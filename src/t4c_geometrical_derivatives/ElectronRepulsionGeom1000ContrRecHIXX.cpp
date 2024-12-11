#include "ElectronRepulsionGeom1000ContrRecHIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_hixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_hixx,
                                            const size_t idx_gixx,
                                            const size_t idx_geom_10_gixx,
                                            const size_t idx_geom_10_gkxx,
                                            const TPoint<double>& r_ab,
                                            const int c_angmom,
                                            const int d_angmom) -> void
{
    const auto nelems = cbuffer.number_of_active_elements();

    const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});

    const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});

    // set up R(AB) distances

    const auto xyz = r_ab.coordinates();

    const auto ab_x = xyz[0];

    const auto ab_y = xyz[1];

    const auto ab_z = xyz[2];

    for (int i = 0; i < ccomps; i++)
    {
        for (int j = 0; j < dcomps; j++)
        {
            /// Set up components of auxilary buffer : GISS

            const auto gi_off = idx_gixx + i * dcomps + j;

            auto g_xxxx_xxxxxx = cbuffer.data(gi_off + 0 * ccomps * dcomps);

            auto g_xxxx_xxxxxy = cbuffer.data(gi_off + 1 * ccomps * dcomps);

            auto g_xxxx_xxxxxz = cbuffer.data(gi_off + 2 * ccomps * dcomps);

            auto g_xxxx_xxxxyy = cbuffer.data(gi_off + 3 * ccomps * dcomps);

            auto g_xxxx_xxxxyz = cbuffer.data(gi_off + 4 * ccomps * dcomps);

            auto g_xxxx_xxxxzz = cbuffer.data(gi_off + 5 * ccomps * dcomps);

            auto g_xxxx_xxxyyy = cbuffer.data(gi_off + 6 * ccomps * dcomps);

            auto g_xxxx_xxxyyz = cbuffer.data(gi_off + 7 * ccomps * dcomps);

            auto g_xxxx_xxxyzz = cbuffer.data(gi_off + 8 * ccomps * dcomps);

            auto g_xxxx_xxxzzz = cbuffer.data(gi_off + 9 * ccomps * dcomps);

            auto g_xxxx_xxyyyy = cbuffer.data(gi_off + 10 * ccomps * dcomps);

            auto g_xxxx_xxyyyz = cbuffer.data(gi_off + 11 * ccomps * dcomps);

            auto g_xxxx_xxyyzz = cbuffer.data(gi_off + 12 * ccomps * dcomps);

            auto g_xxxx_xxyzzz = cbuffer.data(gi_off + 13 * ccomps * dcomps);

            auto g_xxxx_xxzzzz = cbuffer.data(gi_off + 14 * ccomps * dcomps);

            auto g_xxxx_xyyyyy = cbuffer.data(gi_off + 15 * ccomps * dcomps);

            auto g_xxxx_xyyyyz = cbuffer.data(gi_off + 16 * ccomps * dcomps);

            auto g_xxxx_xyyyzz = cbuffer.data(gi_off + 17 * ccomps * dcomps);

            auto g_xxxx_xyyzzz = cbuffer.data(gi_off + 18 * ccomps * dcomps);

            auto g_xxxx_xyzzzz = cbuffer.data(gi_off + 19 * ccomps * dcomps);

            auto g_xxxx_xzzzzz = cbuffer.data(gi_off + 20 * ccomps * dcomps);

            auto g_xxxx_yyyyyy = cbuffer.data(gi_off + 21 * ccomps * dcomps);

            auto g_xxxx_yyyyyz = cbuffer.data(gi_off + 22 * ccomps * dcomps);

            auto g_xxxx_yyyyzz = cbuffer.data(gi_off + 23 * ccomps * dcomps);

            auto g_xxxx_yyyzzz = cbuffer.data(gi_off + 24 * ccomps * dcomps);

            auto g_xxxx_yyzzzz = cbuffer.data(gi_off + 25 * ccomps * dcomps);

            auto g_xxxx_yzzzzz = cbuffer.data(gi_off + 26 * ccomps * dcomps);

            auto g_xxxx_zzzzzz = cbuffer.data(gi_off + 27 * ccomps * dcomps);

            auto g_yyyy_xxxxxx = cbuffer.data(gi_off + 280 * ccomps * dcomps);

            auto g_yyyy_xxxxxy = cbuffer.data(gi_off + 281 * ccomps * dcomps);

            auto g_yyyy_xxxxxz = cbuffer.data(gi_off + 282 * ccomps * dcomps);

            auto g_yyyy_xxxxyy = cbuffer.data(gi_off + 283 * ccomps * dcomps);

            auto g_yyyy_xxxxyz = cbuffer.data(gi_off + 284 * ccomps * dcomps);

            auto g_yyyy_xxxxzz = cbuffer.data(gi_off + 285 * ccomps * dcomps);

            auto g_yyyy_xxxyyy = cbuffer.data(gi_off + 286 * ccomps * dcomps);

            auto g_yyyy_xxxyyz = cbuffer.data(gi_off + 287 * ccomps * dcomps);

            auto g_yyyy_xxxyzz = cbuffer.data(gi_off + 288 * ccomps * dcomps);

            auto g_yyyy_xxxzzz = cbuffer.data(gi_off + 289 * ccomps * dcomps);

            auto g_yyyy_xxyyyy = cbuffer.data(gi_off + 290 * ccomps * dcomps);

            auto g_yyyy_xxyyyz = cbuffer.data(gi_off + 291 * ccomps * dcomps);

            auto g_yyyy_xxyyzz = cbuffer.data(gi_off + 292 * ccomps * dcomps);

            auto g_yyyy_xxyzzz = cbuffer.data(gi_off + 293 * ccomps * dcomps);

            auto g_yyyy_xxzzzz = cbuffer.data(gi_off + 294 * ccomps * dcomps);

            auto g_yyyy_xyyyyy = cbuffer.data(gi_off + 295 * ccomps * dcomps);

            auto g_yyyy_xyyyyz = cbuffer.data(gi_off + 296 * ccomps * dcomps);

            auto g_yyyy_xyyyzz = cbuffer.data(gi_off + 297 * ccomps * dcomps);

            auto g_yyyy_xyyzzz = cbuffer.data(gi_off + 298 * ccomps * dcomps);

            auto g_yyyy_xyzzzz = cbuffer.data(gi_off + 299 * ccomps * dcomps);

            auto g_yyyy_xzzzzz = cbuffer.data(gi_off + 300 * ccomps * dcomps);

            auto g_yyyy_yyyyyy = cbuffer.data(gi_off + 301 * ccomps * dcomps);

            auto g_yyyy_yyyyyz = cbuffer.data(gi_off + 302 * ccomps * dcomps);

            auto g_yyyy_yyyyzz = cbuffer.data(gi_off + 303 * ccomps * dcomps);

            auto g_yyyy_yyyzzz = cbuffer.data(gi_off + 304 * ccomps * dcomps);

            auto g_yyyy_yyzzzz = cbuffer.data(gi_off + 305 * ccomps * dcomps);

            auto g_yyyy_yzzzzz = cbuffer.data(gi_off + 306 * ccomps * dcomps);

            auto g_yyyy_zzzzzz = cbuffer.data(gi_off + 307 * ccomps * dcomps);

            auto g_zzzz_xxxxxx = cbuffer.data(gi_off + 392 * ccomps * dcomps);

            auto g_zzzz_xxxxxy = cbuffer.data(gi_off + 393 * ccomps * dcomps);

            auto g_zzzz_xxxxxz = cbuffer.data(gi_off + 394 * ccomps * dcomps);

            auto g_zzzz_xxxxyy = cbuffer.data(gi_off + 395 * ccomps * dcomps);

            auto g_zzzz_xxxxyz = cbuffer.data(gi_off + 396 * ccomps * dcomps);

            auto g_zzzz_xxxxzz = cbuffer.data(gi_off + 397 * ccomps * dcomps);

            auto g_zzzz_xxxyyy = cbuffer.data(gi_off + 398 * ccomps * dcomps);

            auto g_zzzz_xxxyyz = cbuffer.data(gi_off + 399 * ccomps * dcomps);

            auto g_zzzz_xxxyzz = cbuffer.data(gi_off + 400 * ccomps * dcomps);

            auto g_zzzz_xxxzzz = cbuffer.data(gi_off + 401 * ccomps * dcomps);

            auto g_zzzz_xxyyyy = cbuffer.data(gi_off + 402 * ccomps * dcomps);

            auto g_zzzz_xxyyyz = cbuffer.data(gi_off + 403 * ccomps * dcomps);

            auto g_zzzz_xxyyzz = cbuffer.data(gi_off + 404 * ccomps * dcomps);

            auto g_zzzz_xxyzzz = cbuffer.data(gi_off + 405 * ccomps * dcomps);

            auto g_zzzz_xxzzzz = cbuffer.data(gi_off + 406 * ccomps * dcomps);

            auto g_zzzz_xyyyyy = cbuffer.data(gi_off + 407 * ccomps * dcomps);

            auto g_zzzz_xyyyyz = cbuffer.data(gi_off + 408 * ccomps * dcomps);

            auto g_zzzz_xyyyzz = cbuffer.data(gi_off + 409 * ccomps * dcomps);

            auto g_zzzz_xyyzzz = cbuffer.data(gi_off + 410 * ccomps * dcomps);

            auto g_zzzz_xyzzzz = cbuffer.data(gi_off + 411 * ccomps * dcomps);

            auto g_zzzz_xzzzzz = cbuffer.data(gi_off + 412 * ccomps * dcomps);

            auto g_zzzz_yyyyyy = cbuffer.data(gi_off + 413 * ccomps * dcomps);

            auto g_zzzz_yyyyyz = cbuffer.data(gi_off + 414 * ccomps * dcomps);

            auto g_zzzz_yyyyzz = cbuffer.data(gi_off + 415 * ccomps * dcomps);

            auto g_zzzz_yyyzzz = cbuffer.data(gi_off + 416 * ccomps * dcomps);

            auto g_zzzz_yyzzzz = cbuffer.data(gi_off + 417 * ccomps * dcomps);

            auto g_zzzz_yzzzzz = cbuffer.data(gi_off + 418 * ccomps * dcomps);

            auto g_zzzz_zzzzzz = cbuffer.data(gi_off + 419 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GISS

            const auto gi_geom_10_off = idx_geom_10_gixx + i * dcomps + j;

            auto g_x_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 251 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 279 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 307 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 335 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 356 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 360 * ccomps * dcomps);

            auto g_x_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 362 * ccomps * dcomps);

            auto g_x_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 363 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 364 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 365 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 366 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 369 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 373 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 374 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 378 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 380 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 384 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 387 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 391 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 392 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 395 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 396 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 398 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 401 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 404 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 405 * ccomps * dcomps);

            auto g_x_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 410 * ccomps * dcomps);

            auto g_x_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 416 * ccomps * dcomps);

            auto g_x_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 419 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 434 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 439 * ccomps * dcomps);

            auto g_y_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 449 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 459 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 464 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 469 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 475 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 479 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 489 * ccomps * dcomps);

            auto g_y_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 494 * ccomps * dcomps);

            auto g_y_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 499 * ccomps * dcomps);

            auto g_y_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 503 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 509 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 519 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 524 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 529 * ccomps * dcomps);

            auto g_y_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 531 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 539 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 559 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 560 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 561 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 562 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 563 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 564 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 565 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 566 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 567 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 568 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 569 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 570 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 571 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 572 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 573 * ccomps * dcomps);

            auto g_y_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 574 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 575 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 576 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 577 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 578 * ccomps * dcomps);

            auto g_y_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 579 * ccomps * dcomps);

            auto g_y_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 580 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 581 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 582 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 583 * ccomps * dcomps);

            auto g_y_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 584 * ccomps * dcomps);

            auto g_y_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 585 * ccomps * dcomps);

            auto g_y_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 586 * ccomps * dcomps);

            auto g_y_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 587 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 608 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 609 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 610 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 611 * ccomps * dcomps);

            auto g_y_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 614 * ccomps * dcomps);

            auto g_y_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 615 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 629 * ccomps * dcomps);

            auto g_y_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 630 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 631 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 632 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 633 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 634 * ccomps * dcomps);

            auto g_y_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 635 * ccomps * dcomps);

            auto g_y_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 636 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 637 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 638 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 639 * ccomps * dcomps);

            auto g_y_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 640 * ccomps * dcomps);

            auto g_y_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 641 * ccomps * dcomps);

            auto g_y_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 642 * ccomps * dcomps);

            auto g_y_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 643 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 644 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 645 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 646 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 647 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 648 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 649 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 650 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 651 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 652 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 653 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 654 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 655 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 656 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 657 * ccomps * dcomps);

            auto g_y_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 658 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 659 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 660 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 661 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 662 * ccomps * dcomps);

            auto g_y_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 663 * ccomps * dcomps);

            auto g_y_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 664 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 665 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 666 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 667 * ccomps * dcomps);

            auto g_y_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 668 * ccomps * dcomps);

            auto g_y_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 669 * ccomps * dcomps);

            auto g_y_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 670 * ccomps * dcomps);

            auto g_y_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 671 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 672 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 673 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 674 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 675 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 676 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 677 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 678 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 679 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 680 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 681 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 682 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 683 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 684 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 685 * ccomps * dcomps);

            auto g_y_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 686 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 687 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 688 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 689 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 690 * ccomps * dcomps);

            auto g_y_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 691 * ccomps * dcomps);

            auto g_y_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 692 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 693 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 694 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 695 * ccomps * dcomps);

            auto g_y_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 696 * ccomps * dcomps);

            auto g_y_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 697 * ccomps * dcomps);

            auto g_y_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 698 * ccomps * dcomps);

            auto g_y_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 699 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 700 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 701 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 702 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 703 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 704 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 705 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 706 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 707 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 708 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 709 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 710 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 711 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 712 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 713 * ccomps * dcomps);

            auto g_y_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 714 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 715 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 716 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 717 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 718 * ccomps * dcomps);

            auto g_y_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 719 * ccomps * dcomps);

            auto g_y_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 720 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 721 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 722 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 723 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 724 * ccomps * dcomps);

            auto g_y_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 725 * ccomps * dcomps);

            auto g_y_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 726 * ccomps * dcomps);

            auto g_y_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 727 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 728 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 729 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 730 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 731 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 732 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 733 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 734 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 735 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 736 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 737 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 738 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 739 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 740 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 741 * ccomps * dcomps);

            auto g_y_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 742 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 743 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 744 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 745 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 746 * ccomps * dcomps);

            auto g_y_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 747 * ccomps * dcomps);

            auto g_y_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 748 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 749 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 750 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 751 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 752 * ccomps * dcomps);

            auto g_y_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 753 * ccomps * dcomps);

            auto g_y_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 754 * ccomps * dcomps);

            auto g_y_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 755 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 756 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 757 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 758 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 759 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 760 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 761 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 762 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 763 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 764 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 765 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 766 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 767 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 768 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 769 * ccomps * dcomps);

            auto g_y_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 770 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 771 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 772 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 773 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 774 * ccomps * dcomps);

            auto g_y_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 775 * ccomps * dcomps);

            auto g_y_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 776 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 777 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 778 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 779 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 780 * ccomps * dcomps);

            auto g_y_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 781 * ccomps * dcomps);

            auto g_y_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 782 * ccomps * dcomps);

            auto g_y_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 783 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 784 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 785 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 786 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 787 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 788 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 789 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 790 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 791 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 794 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 797 * ccomps * dcomps);

            auto g_y_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 805 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 809 * ccomps * dcomps);

            auto g_y_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 811 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 818 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 820 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 821 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 822 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 823 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 824 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 825 * ccomps * dcomps);

            auto g_y_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 826 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 827 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 833 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 839 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxx = cbuffer.data(gi_geom_10_off + 840 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxy = cbuffer.data(gi_geom_10_off + 841 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxz = cbuffer.data(gi_geom_10_off + 842 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxyy = cbuffer.data(gi_geom_10_off + 843 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxyz = cbuffer.data(gi_geom_10_off + 844 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxzz = cbuffer.data(gi_geom_10_off + 845 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyyy = cbuffer.data(gi_geom_10_off + 846 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyyz = cbuffer.data(gi_geom_10_off + 847 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyzz = cbuffer.data(gi_geom_10_off + 848 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxzzz = cbuffer.data(gi_geom_10_off + 849 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyyy = cbuffer.data(gi_geom_10_off + 850 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyyz = cbuffer.data(gi_geom_10_off + 851 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyzz = cbuffer.data(gi_geom_10_off + 852 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyzzz = cbuffer.data(gi_geom_10_off + 853 * ccomps * dcomps);

            auto g_z_0_xxxx_xxzzzz = cbuffer.data(gi_geom_10_off + 854 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyyy = cbuffer.data(gi_geom_10_off + 855 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyyz = cbuffer.data(gi_geom_10_off + 856 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyzz = cbuffer.data(gi_geom_10_off + 857 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyzzz = cbuffer.data(gi_geom_10_off + 858 * ccomps * dcomps);

            auto g_z_0_xxxx_xyzzzz = cbuffer.data(gi_geom_10_off + 859 * ccomps * dcomps);

            auto g_z_0_xxxx_xzzzzz = cbuffer.data(gi_geom_10_off + 860 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyyyy = cbuffer.data(gi_geom_10_off + 861 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyyyz = cbuffer.data(gi_geom_10_off + 862 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyyzz = cbuffer.data(gi_geom_10_off + 863 * ccomps * dcomps);

            auto g_z_0_xxxx_yyyzzz = cbuffer.data(gi_geom_10_off + 864 * ccomps * dcomps);

            auto g_z_0_xxxx_yyzzzz = cbuffer.data(gi_geom_10_off + 865 * ccomps * dcomps);

            auto g_z_0_xxxx_yzzzzz = cbuffer.data(gi_geom_10_off + 866 * ccomps * dcomps);

            auto g_z_0_xxxx_zzzzzz = cbuffer.data(gi_geom_10_off + 867 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxx = cbuffer.data(gi_geom_10_off + 868 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxy = cbuffer.data(gi_geom_10_off + 869 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxz = cbuffer.data(gi_geom_10_off + 870 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxyy = cbuffer.data(gi_geom_10_off + 871 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxyz = cbuffer.data(gi_geom_10_off + 872 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxzz = cbuffer.data(gi_geom_10_off + 873 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyyy = cbuffer.data(gi_geom_10_off + 874 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyyz = cbuffer.data(gi_geom_10_off + 875 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyzz = cbuffer.data(gi_geom_10_off + 876 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxzzz = cbuffer.data(gi_geom_10_off + 877 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyyy = cbuffer.data(gi_geom_10_off + 878 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyyz = cbuffer.data(gi_geom_10_off + 879 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyzz = cbuffer.data(gi_geom_10_off + 880 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyzzz = cbuffer.data(gi_geom_10_off + 881 * ccomps * dcomps);

            auto g_z_0_xxxy_xxzzzz = cbuffer.data(gi_geom_10_off + 882 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyyy = cbuffer.data(gi_geom_10_off + 883 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyyz = cbuffer.data(gi_geom_10_off + 884 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyzz = cbuffer.data(gi_geom_10_off + 885 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyzzz = cbuffer.data(gi_geom_10_off + 886 * ccomps * dcomps);

            auto g_z_0_xxxy_xyzzzz = cbuffer.data(gi_geom_10_off + 887 * ccomps * dcomps);

            auto g_z_0_xxxy_xzzzzz = cbuffer.data(gi_geom_10_off + 888 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyyyy = cbuffer.data(gi_geom_10_off + 889 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyyyz = cbuffer.data(gi_geom_10_off + 890 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyyzz = cbuffer.data(gi_geom_10_off + 891 * ccomps * dcomps);

            auto g_z_0_xxxy_yyyzzz = cbuffer.data(gi_geom_10_off + 892 * ccomps * dcomps);

            auto g_z_0_xxxy_yyzzzz = cbuffer.data(gi_geom_10_off + 893 * ccomps * dcomps);

            auto g_z_0_xxxy_yzzzzz = cbuffer.data(gi_geom_10_off + 894 * ccomps * dcomps);

            auto g_z_0_xxxy_zzzzzz = cbuffer.data(gi_geom_10_off + 895 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxx = cbuffer.data(gi_geom_10_off + 896 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxy = cbuffer.data(gi_geom_10_off + 897 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxz = cbuffer.data(gi_geom_10_off + 898 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxyy = cbuffer.data(gi_geom_10_off + 899 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxyz = cbuffer.data(gi_geom_10_off + 900 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxzz = cbuffer.data(gi_geom_10_off + 901 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyyy = cbuffer.data(gi_geom_10_off + 902 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyyz = cbuffer.data(gi_geom_10_off + 903 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyzz = cbuffer.data(gi_geom_10_off + 904 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxzzz = cbuffer.data(gi_geom_10_off + 905 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyyy = cbuffer.data(gi_geom_10_off + 906 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyyz = cbuffer.data(gi_geom_10_off + 907 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyzz = cbuffer.data(gi_geom_10_off + 908 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyzzz = cbuffer.data(gi_geom_10_off + 909 * ccomps * dcomps);

            auto g_z_0_xxxz_xxzzzz = cbuffer.data(gi_geom_10_off + 910 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyyy = cbuffer.data(gi_geom_10_off + 911 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyyz = cbuffer.data(gi_geom_10_off + 912 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyzz = cbuffer.data(gi_geom_10_off + 913 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyzzz = cbuffer.data(gi_geom_10_off + 914 * ccomps * dcomps);

            auto g_z_0_xxxz_xyzzzz = cbuffer.data(gi_geom_10_off + 915 * ccomps * dcomps);

            auto g_z_0_xxxz_xzzzzz = cbuffer.data(gi_geom_10_off + 916 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyyyy = cbuffer.data(gi_geom_10_off + 917 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyyyz = cbuffer.data(gi_geom_10_off + 918 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyyzz = cbuffer.data(gi_geom_10_off + 919 * ccomps * dcomps);

            auto g_z_0_xxxz_yyyzzz = cbuffer.data(gi_geom_10_off + 920 * ccomps * dcomps);

            auto g_z_0_xxxz_yyzzzz = cbuffer.data(gi_geom_10_off + 921 * ccomps * dcomps);

            auto g_z_0_xxxz_yzzzzz = cbuffer.data(gi_geom_10_off + 922 * ccomps * dcomps);

            auto g_z_0_xxxz_zzzzzz = cbuffer.data(gi_geom_10_off + 923 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxx = cbuffer.data(gi_geom_10_off + 924 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxy = cbuffer.data(gi_geom_10_off + 925 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxz = cbuffer.data(gi_geom_10_off + 926 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxyy = cbuffer.data(gi_geom_10_off + 927 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxyz = cbuffer.data(gi_geom_10_off + 928 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxzz = cbuffer.data(gi_geom_10_off + 929 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyyy = cbuffer.data(gi_geom_10_off + 930 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyyz = cbuffer.data(gi_geom_10_off + 931 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyzz = cbuffer.data(gi_geom_10_off + 932 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxzzz = cbuffer.data(gi_geom_10_off + 933 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyyy = cbuffer.data(gi_geom_10_off + 934 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyyz = cbuffer.data(gi_geom_10_off + 935 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyzz = cbuffer.data(gi_geom_10_off + 936 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyzzz = cbuffer.data(gi_geom_10_off + 937 * ccomps * dcomps);

            auto g_z_0_xxyy_xxzzzz = cbuffer.data(gi_geom_10_off + 938 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyyy = cbuffer.data(gi_geom_10_off + 939 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyyz = cbuffer.data(gi_geom_10_off + 940 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyzz = cbuffer.data(gi_geom_10_off + 941 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyzzz = cbuffer.data(gi_geom_10_off + 942 * ccomps * dcomps);

            auto g_z_0_xxyy_xyzzzz = cbuffer.data(gi_geom_10_off + 943 * ccomps * dcomps);

            auto g_z_0_xxyy_xzzzzz = cbuffer.data(gi_geom_10_off + 944 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyyyy = cbuffer.data(gi_geom_10_off + 945 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyyyz = cbuffer.data(gi_geom_10_off + 946 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyyzz = cbuffer.data(gi_geom_10_off + 947 * ccomps * dcomps);

            auto g_z_0_xxyy_yyyzzz = cbuffer.data(gi_geom_10_off + 948 * ccomps * dcomps);

            auto g_z_0_xxyy_yyzzzz = cbuffer.data(gi_geom_10_off + 949 * ccomps * dcomps);

            auto g_z_0_xxyy_yzzzzz = cbuffer.data(gi_geom_10_off + 950 * ccomps * dcomps);

            auto g_z_0_xxyy_zzzzzz = cbuffer.data(gi_geom_10_off + 951 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxx = cbuffer.data(gi_geom_10_off + 952 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxy = cbuffer.data(gi_geom_10_off + 953 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxz = cbuffer.data(gi_geom_10_off + 954 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxyy = cbuffer.data(gi_geom_10_off + 955 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxyz = cbuffer.data(gi_geom_10_off + 956 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxzz = cbuffer.data(gi_geom_10_off + 957 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyyy = cbuffer.data(gi_geom_10_off + 958 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyyz = cbuffer.data(gi_geom_10_off + 959 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyzz = cbuffer.data(gi_geom_10_off + 960 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxzzz = cbuffer.data(gi_geom_10_off + 961 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyyy = cbuffer.data(gi_geom_10_off + 962 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyyz = cbuffer.data(gi_geom_10_off + 963 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyzz = cbuffer.data(gi_geom_10_off + 964 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyzzz = cbuffer.data(gi_geom_10_off + 965 * ccomps * dcomps);

            auto g_z_0_xxyz_xxzzzz = cbuffer.data(gi_geom_10_off + 966 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyyy = cbuffer.data(gi_geom_10_off + 967 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyyz = cbuffer.data(gi_geom_10_off + 968 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyzz = cbuffer.data(gi_geom_10_off + 969 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyzzz = cbuffer.data(gi_geom_10_off + 970 * ccomps * dcomps);

            auto g_z_0_xxyz_xyzzzz = cbuffer.data(gi_geom_10_off + 971 * ccomps * dcomps);

            auto g_z_0_xxyz_xzzzzz = cbuffer.data(gi_geom_10_off + 972 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyyyy = cbuffer.data(gi_geom_10_off + 973 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyyyz = cbuffer.data(gi_geom_10_off + 974 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyyzz = cbuffer.data(gi_geom_10_off + 975 * ccomps * dcomps);

            auto g_z_0_xxyz_yyyzzz = cbuffer.data(gi_geom_10_off + 976 * ccomps * dcomps);

            auto g_z_0_xxyz_yyzzzz = cbuffer.data(gi_geom_10_off + 977 * ccomps * dcomps);

            auto g_z_0_xxyz_yzzzzz = cbuffer.data(gi_geom_10_off + 978 * ccomps * dcomps);

            auto g_z_0_xxyz_zzzzzz = cbuffer.data(gi_geom_10_off + 979 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxx = cbuffer.data(gi_geom_10_off + 980 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxy = cbuffer.data(gi_geom_10_off + 981 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxz = cbuffer.data(gi_geom_10_off + 982 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxyy = cbuffer.data(gi_geom_10_off + 983 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxyz = cbuffer.data(gi_geom_10_off + 984 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxzz = cbuffer.data(gi_geom_10_off + 985 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyyy = cbuffer.data(gi_geom_10_off + 986 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyyz = cbuffer.data(gi_geom_10_off + 987 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyzz = cbuffer.data(gi_geom_10_off + 988 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxzzz = cbuffer.data(gi_geom_10_off + 989 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyyy = cbuffer.data(gi_geom_10_off + 990 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyyz = cbuffer.data(gi_geom_10_off + 991 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyzz = cbuffer.data(gi_geom_10_off + 992 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyzzz = cbuffer.data(gi_geom_10_off + 993 * ccomps * dcomps);

            auto g_z_0_xxzz_xxzzzz = cbuffer.data(gi_geom_10_off + 994 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyyy = cbuffer.data(gi_geom_10_off + 995 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyyz = cbuffer.data(gi_geom_10_off + 996 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyzz = cbuffer.data(gi_geom_10_off + 997 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyzzz = cbuffer.data(gi_geom_10_off + 998 * ccomps * dcomps);

            auto g_z_0_xxzz_xyzzzz = cbuffer.data(gi_geom_10_off + 999 * ccomps * dcomps);

            auto g_z_0_xxzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1000 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1001 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1002 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1003 * ccomps * dcomps);

            auto g_z_0_xxzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1004 * ccomps * dcomps);

            auto g_z_0_xxzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1005 * ccomps * dcomps);

            auto g_z_0_xxzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1006 * ccomps * dcomps);

            auto g_z_0_xxzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1007 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 1008 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 1009 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 1010 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 1011 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 1012 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 1013 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 1014 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 1015 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 1016 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 1017 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 1018 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 1019 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 1020 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 1021 * ccomps * dcomps);

            auto g_z_0_xyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 1022 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 1023 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 1024 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 1025 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 1026 * ccomps * dcomps);

            auto g_z_0_xyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 1027 * ccomps * dcomps);

            auto g_z_0_xyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 1028 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 1029 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 1030 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 1031 * ccomps * dcomps);

            auto g_z_0_xyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 1032 * ccomps * dcomps);

            auto g_z_0_xyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 1033 * ccomps * dcomps);

            auto g_z_0_xyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 1034 * ccomps * dcomps);

            auto g_z_0_xyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 1035 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 1036 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 1037 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 1038 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 1039 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 1040 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 1041 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 1042 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 1043 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 1044 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 1045 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 1046 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 1047 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 1048 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 1049 * ccomps * dcomps);

            auto g_z_0_xyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 1050 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 1051 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 1052 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 1053 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 1054 * ccomps * dcomps);

            auto g_z_0_xyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 1055 * ccomps * dcomps);

            auto g_z_0_xyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 1056 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 1057 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 1058 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 1059 * ccomps * dcomps);

            auto g_z_0_xyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 1060 * ccomps * dcomps);

            auto g_z_0_xyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 1061 * ccomps * dcomps);

            auto g_z_0_xyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 1062 * ccomps * dcomps);

            auto g_z_0_xyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 1063 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 1064 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 1065 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 1066 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 1067 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 1068 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 1069 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 1070 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 1071 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 1072 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 1073 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 1074 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 1075 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 1076 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 1077 * ccomps * dcomps);

            auto g_z_0_xyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 1078 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 1079 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 1080 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 1081 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 1082 * ccomps * dcomps);

            auto g_z_0_xyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 1083 * ccomps * dcomps);

            auto g_z_0_xyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1084 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1085 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1086 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1087 * ccomps * dcomps);

            auto g_z_0_xyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1088 * ccomps * dcomps);

            auto g_z_0_xyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1089 * ccomps * dcomps);

            auto g_z_0_xyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1090 * ccomps * dcomps);

            auto g_z_0_xyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1091 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 1092 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 1093 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 1094 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 1095 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 1096 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 1097 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 1098 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 1099 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 1100 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 1101 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 1102 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 1103 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 1104 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 1105 * ccomps * dcomps);

            auto g_z_0_xzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 1106 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 1107 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 1108 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 1109 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 1110 * ccomps * dcomps);

            auto g_z_0_xzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 1111 * ccomps * dcomps);

            auto g_z_0_xzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1112 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1113 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1114 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1115 * ccomps * dcomps);

            auto g_z_0_xzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1116 * ccomps * dcomps);

            auto g_z_0_xzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1117 * ccomps * dcomps);

            auto g_z_0_xzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1118 * ccomps * dcomps);

            auto g_z_0_xzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1119 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxx = cbuffer.data(gi_geom_10_off + 1120 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxy = cbuffer.data(gi_geom_10_off + 1121 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxz = cbuffer.data(gi_geom_10_off + 1122 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxyy = cbuffer.data(gi_geom_10_off + 1123 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxyz = cbuffer.data(gi_geom_10_off + 1124 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxzz = cbuffer.data(gi_geom_10_off + 1125 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyyy = cbuffer.data(gi_geom_10_off + 1126 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyyz = cbuffer.data(gi_geom_10_off + 1127 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyzz = cbuffer.data(gi_geom_10_off + 1128 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxzzz = cbuffer.data(gi_geom_10_off + 1129 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyyy = cbuffer.data(gi_geom_10_off + 1130 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyyz = cbuffer.data(gi_geom_10_off + 1131 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyzz = cbuffer.data(gi_geom_10_off + 1132 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyzzz = cbuffer.data(gi_geom_10_off + 1133 * ccomps * dcomps);

            auto g_z_0_yyyy_xxzzzz = cbuffer.data(gi_geom_10_off + 1134 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyyy = cbuffer.data(gi_geom_10_off + 1135 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyyz = cbuffer.data(gi_geom_10_off + 1136 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyzz = cbuffer.data(gi_geom_10_off + 1137 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyzzz = cbuffer.data(gi_geom_10_off + 1138 * ccomps * dcomps);

            auto g_z_0_yyyy_xyzzzz = cbuffer.data(gi_geom_10_off + 1139 * ccomps * dcomps);

            auto g_z_0_yyyy_xzzzzz = cbuffer.data(gi_geom_10_off + 1140 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyyy = cbuffer.data(gi_geom_10_off + 1141 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyyz = cbuffer.data(gi_geom_10_off + 1142 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyzz = cbuffer.data(gi_geom_10_off + 1143 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyzzz = cbuffer.data(gi_geom_10_off + 1144 * ccomps * dcomps);

            auto g_z_0_yyyy_yyzzzz = cbuffer.data(gi_geom_10_off + 1145 * ccomps * dcomps);

            auto g_z_0_yyyy_yzzzzz = cbuffer.data(gi_geom_10_off + 1146 * ccomps * dcomps);

            auto g_z_0_yyyy_zzzzzz = cbuffer.data(gi_geom_10_off + 1147 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxx = cbuffer.data(gi_geom_10_off + 1148 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxy = cbuffer.data(gi_geom_10_off + 1149 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxz = cbuffer.data(gi_geom_10_off + 1150 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxyy = cbuffer.data(gi_geom_10_off + 1151 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxyz = cbuffer.data(gi_geom_10_off + 1152 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxzz = cbuffer.data(gi_geom_10_off + 1153 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyyy = cbuffer.data(gi_geom_10_off + 1154 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyyz = cbuffer.data(gi_geom_10_off + 1155 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyzz = cbuffer.data(gi_geom_10_off + 1156 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxzzz = cbuffer.data(gi_geom_10_off + 1157 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyyy = cbuffer.data(gi_geom_10_off + 1158 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyyz = cbuffer.data(gi_geom_10_off + 1159 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyzz = cbuffer.data(gi_geom_10_off + 1160 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyzzz = cbuffer.data(gi_geom_10_off + 1161 * ccomps * dcomps);

            auto g_z_0_yyyz_xxzzzz = cbuffer.data(gi_geom_10_off + 1162 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyyy = cbuffer.data(gi_geom_10_off + 1163 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyyz = cbuffer.data(gi_geom_10_off + 1164 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyzz = cbuffer.data(gi_geom_10_off + 1165 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyzzz = cbuffer.data(gi_geom_10_off + 1166 * ccomps * dcomps);

            auto g_z_0_yyyz_xyzzzz = cbuffer.data(gi_geom_10_off + 1167 * ccomps * dcomps);

            auto g_z_0_yyyz_xzzzzz = cbuffer.data(gi_geom_10_off + 1168 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyyy = cbuffer.data(gi_geom_10_off + 1169 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyyz = cbuffer.data(gi_geom_10_off + 1170 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyzz = cbuffer.data(gi_geom_10_off + 1171 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyzzz = cbuffer.data(gi_geom_10_off + 1172 * ccomps * dcomps);

            auto g_z_0_yyyz_yyzzzz = cbuffer.data(gi_geom_10_off + 1173 * ccomps * dcomps);

            auto g_z_0_yyyz_yzzzzz = cbuffer.data(gi_geom_10_off + 1174 * ccomps * dcomps);

            auto g_z_0_yyyz_zzzzzz = cbuffer.data(gi_geom_10_off + 1175 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxx = cbuffer.data(gi_geom_10_off + 1176 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxy = cbuffer.data(gi_geom_10_off + 1177 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxz = cbuffer.data(gi_geom_10_off + 1178 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxyy = cbuffer.data(gi_geom_10_off + 1179 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxyz = cbuffer.data(gi_geom_10_off + 1180 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxzz = cbuffer.data(gi_geom_10_off + 1181 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyyy = cbuffer.data(gi_geom_10_off + 1182 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyyz = cbuffer.data(gi_geom_10_off + 1183 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyzz = cbuffer.data(gi_geom_10_off + 1184 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxzzz = cbuffer.data(gi_geom_10_off + 1185 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyyy = cbuffer.data(gi_geom_10_off + 1186 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyyz = cbuffer.data(gi_geom_10_off + 1187 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyzz = cbuffer.data(gi_geom_10_off + 1188 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyzzz = cbuffer.data(gi_geom_10_off + 1189 * ccomps * dcomps);

            auto g_z_0_yyzz_xxzzzz = cbuffer.data(gi_geom_10_off + 1190 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyyy = cbuffer.data(gi_geom_10_off + 1191 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyyz = cbuffer.data(gi_geom_10_off + 1192 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyzz = cbuffer.data(gi_geom_10_off + 1193 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyzzz = cbuffer.data(gi_geom_10_off + 1194 * ccomps * dcomps);

            auto g_z_0_yyzz_xyzzzz = cbuffer.data(gi_geom_10_off + 1195 * ccomps * dcomps);

            auto g_z_0_yyzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1196 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1197 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1198 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1199 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1200 * ccomps * dcomps);

            auto g_z_0_yyzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1201 * ccomps * dcomps);

            auto g_z_0_yyzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1202 * ccomps * dcomps);

            auto g_z_0_yyzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1203 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 1204 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 1205 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 1206 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 1207 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 1208 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 1209 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 1210 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 1211 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 1212 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 1213 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 1214 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 1215 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 1216 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 1217 * ccomps * dcomps);

            auto g_z_0_yzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 1218 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 1219 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 1220 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 1221 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 1222 * ccomps * dcomps);

            auto g_z_0_yzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 1223 * ccomps * dcomps);

            auto g_z_0_yzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1224 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1225 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1226 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1227 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1228 * ccomps * dcomps);

            auto g_z_0_yzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1229 * ccomps * dcomps);

            auto g_z_0_yzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1230 * ccomps * dcomps);

            auto g_z_0_yzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1231 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxx = cbuffer.data(gi_geom_10_off + 1232 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxy = cbuffer.data(gi_geom_10_off + 1233 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxz = cbuffer.data(gi_geom_10_off + 1234 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxyy = cbuffer.data(gi_geom_10_off + 1235 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxyz = cbuffer.data(gi_geom_10_off + 1236 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxzz = cbuffer.data(gi_geom_10_off + 1237 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyyy = cbuffer.data(gi_geom_10_off + 1238 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyyz = cbuffer.data(gi_geom_10_off + 1239 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyzz = cbuffer.data(gi_geom_10_off + 1240 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxzzz = cbuffer.data(gi_geom_10_off + 1241 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyyy = cbuffer.data(gi_geom_10_off + 1242 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyyz = cbuffer.data(gi_geom_10_off + 1243 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyzz = cbuffer.data(gi_geom_10_off + 1244 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyzzz = cbuffer.data(gi_geom_10_off + 1245 * ccomps * dcomps);

            auto g_z_0_zzzz_xxzzzz = cbuffer.data(gi_geom_10_off + 1246 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyyy = cbuffer.data(gi_geom_10_off + 1247 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyyz = cbuffer.data(gi_geom_10_off + 1248 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyzz = cbuffer.data(gi_geom_10_off + 1249 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyzzz = cbuffer.data(gi_geom_10_off + 1250 * ccomps * dcomps);

            auto g_z_0_zzzz_xyzzzz = cbuffer.data(gi_geom_10_off + 1251 * ccomps * dcomps);

            auto g_z_0_zzzz_xzzzzz = cbuffer.data(gi_geom_10_off + 1252 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyyy = cbuffer.data(gi_geom_10_off + 1253 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyyz = cbuffer.data(gi_geom_10_off + 1254 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyzz = cbuffer.data(gi_geom_10_off + 1255 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyzzz = cbuffer.data(gi_geom_10_off + 1256 * ccomps * dcomps);

            auto g_z_0_zzzz_yyzzzz = cbuffer.data(gi_geom_10_off + 1257 * ccomps * dcomps);

            auto g_z_0_zzzz_yzzzzz = cbuffer.data(gi_geom_10_off + 1258 * ccomps * dcomps);

            auto g_z_0_zzzz_zzzzzz = cbuffer.data(gi_geom_10_off + 1259 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GKSS

            const auto gk_geom_10_off = idx_geom_10_gkxx + i * dcomps + j;

            auto g_x_0_xxxx_xxxxxxx = cbuffer.data(gk_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxxxy = cbuffer.data(gk_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxxxz = cbuffer.data(gk_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxxyy = cbuffer.data(gk_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxxyz = cbuffer.data(gk_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxxzz = cbuffer.data(gk_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxyyy = cbuffer.data(gk_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxyyz = cbuffer.data(gk_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxyzz = cbuffer.data(gk_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxxzzz = cbuffer.data(gk_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyyyy = cbuffer.data(gk_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyyyz = cbuffer.data(gk_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyyzz = cbuffer.data(gk_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxyzzz = cbuffer.data(gk_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxx_xxxzzzz = cbuffer.data(gk_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyyyy = cbuffer.data(gk_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyyyz = cbuffer.data(gk_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyyzz = cbuffer.data(gk_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyyzzz = cbuffer.data(gk_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxx_xxyzzzz = cbuffer.data(gk_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxx_xxzzzzz = cbuffer.data(gk_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyyyy = cbuffer.data(gk_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyyyz = cbuffer.data(gk_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyyzz = cbuffer.data(gk_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyyzzz = cbuffer.data(gk_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxx_xyyzzzz = cbuffer.data(gk_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxx_xyzzzzz = cbuffer.data(gk_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxx_xzzzzzz = cbuffer.data(gk_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyyyy = cbuffer.data(gk_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyyyz = cbuffer.data(gk_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyyzz = cbuffer.data(gk_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyyzzz = cbuffer.data(gk_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxx_yyyzzzz = cbuffer.data(gk_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxx_yyzzzzz = cbuffer.data(gk_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxx_yzzzzzz = cbuffer.data(gk_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxx_zzzzzzz = cbuffer.data(gk_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxxxy = cbuffer.data(gk_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxxyy = cbuffer.data(gk_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxxyz = cbuffer.data(gk_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxyyy = cbuffer.data(gk_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxyyz = cbuffer.data(gk_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxxyzz = cbuffer.data(gk_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyyyy = cbuffer.data(gk_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyyyz = cbuffer.data(gk_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyyzz = cbuffer.data(gk_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxy_xxxyzzz = cbuffer.data(gk_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyyyy = cbuffer.data(gk_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyyyz = cbuffer.data(gk_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyyzz = cbuffer.data(gk_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyyzzz = cbuffer.data(gk_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxy_xxyzzzz = cbuffer.data(gk_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyyyy = cbuffer.data(gk_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyyyz = cbuffer.data(gk_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyyzz = cbuffer.data(gk_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyyzzz = cbuffer.data(gk_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxy_xyyzzzz = cbuffer.data(gk_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxy_xyzzzzz = cbuffer.data(gk_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyyyy = cbuffer.data(gk_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyyyz = cbuffer.data(gk_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyyzz = cbuffer.data(gk_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyyzzz = cbuffer.data(gk_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxy_yyyzzzz = cbuffer.data(gk_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxy_yyzzzzz = cbuffer.data(gk_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxy_yzzzzzz = cbuffer.data(gk_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxxxy = cbuffer.data(gk_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxxxz = cbuffer.data(gk_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxxyy = cbuffer.data(gk_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxxyz = cbuffer.data(gk_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxxzz = cbuffer.data(gk_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxyyy = cbuffer.data(gk_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxyyz = cbuffer.data(gk_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxyzz = cbuffer.data(gk_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxxzzz = cbuffer.data(gk_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyyyy = cbuffer.data(gk_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyyyz = cbuffer.data(gk_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyyzz = cbuffer.data(gk_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxyzzz = cbuffer.data(gk_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxz_xxxzzzz = cbuffer.data(gk_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyyyy = cbuffer.data(gk_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyyyz = cbuffer.data(gk_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyyzz = cbuffer.data(gk_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyyzzz = cbuffer.data(gk_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxxz_xxyzzzz = cbuffer.data(gk_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxxz_xxzzzzz = cbuffer.data(gk_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyyyy = cbuffer.data(gk_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyyyz = cbuffer.data(gk_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyyzz = cbuffer.data(gk_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyyzzz = cbuffer.data(gk_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxxz_xyyzzzz = cbuffer.data(gk_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxxz_xyzzzzz = cbuffer.data(gk_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxxz_xzzzzzz = cbuffer.data(gk_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyyyy = cbuffer.data(gk_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyyyz = cbuffer.data(gk_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyyzz = cbuffer.data(gk_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyyzzz = cbuffer.data(gk_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxxz_yyyzzzz = cbuffer.data(gk_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxxz_yyzzzzz = cbuffer.data(gk_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxxz_yzzzzzz = cbuffer.data(gk_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxxz_zzzzzzz = cbuffer.data(gk_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xxyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xxyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xxyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xxyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xxyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xxyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xxzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xxzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xxzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xxzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xxzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xxzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xxzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xxzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xxzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xxzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xyyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xyyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xyyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_xyyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_xyyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_xyyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_xyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_xyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 307 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_xyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_xyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_xyzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_xyzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_xyzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 335 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_xzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_xzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_xzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_xzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_xzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_xzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_xzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 356 * ccomps * dcomps);

            auto g_x_0_xzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_xzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_xzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 363 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 364 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 366 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_yyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 373 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 378 * ccomps * dcomps);

            auto g_x_0_yyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 384 * ccomps * dcomps);

            auto g_x_0_yyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_yyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 391 * ccomps * dcomps);

            auto g_x_0_yyyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 392 * ccomps * dcomps);

            auto g_x_0_yyyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_yyyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 404 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_yyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_yyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 419 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 420 * ccomps * dcomps);

            auto g_x_0_yyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 421 * ccomps * dcomps);

            auto g_x_0_yyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 422 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 424 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 425 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 426 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 427 * ccomps * dcomps);

            auto g_x_0_yyyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 428 * ccomps * dcomps);

            auto g_x_0_yyyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 429 * ccomps * dcomps);

            auto g_x_0_yyyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 430 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 433 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 435 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 436 * ccomps * dcomps); 

            auto g_x_0_yyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 438 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 439 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 440 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 442 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 443 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 444 * ccomps * dcomps);

            auto g_x_0_yyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 445 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 447 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 448 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 449 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 450 * ccomps * dcomps);

            auto g_x_0_yyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 451 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 453 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 454 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 455 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 456 * ccomps * dcomps);

            auto g_x_0_yyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 457 * ccomps * dcomps);

            auto g_x_0_yyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 458 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 460 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 461 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 462 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 463 * ccomps * dcomps);

            auto g_x_0_yyzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 464 * ccomps * dcomps);

            auto g_x_0_yyzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 465 * ccomps * dcomps);

            auto g_x_0_yyzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 466 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 469 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 471 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 472 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 474 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 475 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 476 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 478 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 479 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 480 * ccomps * dcomps);

            auto g_x_0_yzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 481 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 483 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 484 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 485 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 486 * ccomps * dcomps);

            auto g_x_0_yzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 487 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 489 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 490 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 491 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 492 * ccomps * dcomps);

            auto g_x_0_yzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 493 * ccomps * dcomps);

            auto g_x_0_yzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 494 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 496 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 497 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 498 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 499 * ccomps * dcomps);

            auto g_x_0_yzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 500 * ccomps * dcomps);

            auto g_x_0_yzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 501 * ccomps * dcomps);

            auto g_x_0_yzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 502 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 505 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 506 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 507 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 508 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 509 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 510 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 511 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 512 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 513 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 514 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 515 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 516 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 517 * ccomps * dcomps);

            auto g_x_0_zzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 518 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 519 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 520 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 521 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 522 * ccomps * dcomps);

            auto g_x_0_zzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 523 * ccomps * dcomps);

            auto g_x_0_zzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 524 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 525 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 526 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 527 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 528 * ccomps * dcomps);

            auto g_x_0_zzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 529 * ccomps * dcomps);

            auto g_x_0_zzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 530 * ccomps * dcomps);

            auto g_x_0_zzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 531 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 532 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 533 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 534 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 535 * ccomps * dcomps);

            auto g_x_0_zzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 536 * ccomps * dcomps);

            auto g_x_0_zzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 537 * ccomps * dcomps);

            auto g_x_0_zzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 538 * ccomps * dcomps);

            auto g_x_0_zzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 539 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxxx = cbuffer.data(gk_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxxy = cbuffer.data(gk_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxxz = cbuffer.data(gk_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxyy = cbuffer.data(gk_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxyz = cbuffer.data(gk_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxxzz = cbuffer.data(gk_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxyyy = cbuffer.data(gk_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxyyz = cbuffer.data(gk_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxyzz = cbuffer.data(gk_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxxzzz = cbuffer.data(gk_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyyyy = cbuffer.data(gk_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyyyz = cbuffer.data(gk_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyyzz = cbuffer.data(gk_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxyzzz = cbuffer.data(gk_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_xxxx_xxxzzzz = cbuffer.data(gk_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyyyy = cbuffer.data(gk_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyyyz = cbuffer.data(gk_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyyzz = cbuffer.data(gk_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyyzzz = cbuffer.data(gk_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_xxxx_xxyzzzz = cbuffer.data(gk_geom_10_off + 559 * ccomps * dcomps);

            auto g_y_0_xxxx_xxzzzzz = cbuffer.data(gk_geom_10_off + 560 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyyyy = cbuffer.data(gk_geom_10_off + 561 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyyyz = cbuffer.data(gk_geom_10_off + 562 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyyzz = cbuffer.data(gk_geom_10_off + 563 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyyzzz = cbuffer.data(gk_geom_10_off + 564 * ccomps * dcomps);

            auto g_y_0_xxxx_xyyzzzz = cbuffer.data(gk_geom_10_off + 565 * ccomps * dcomps);

            auto g_y_0_xxxx_xyzzzzz = cbuffer.data(gk_geom_10_off + 566 * ccomps * dcomps);

            auto g_y_0_xxxx_xzzzzzz = cbuffer.data(gk_geom_10_off + 567 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxxx = cbuffer.data(gk_geom_10_off + 576 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxxy = cbuffer.data(gk_geom_10_off + 577 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxxz = cbuffer.data(gk_geom_10_off + 578 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxyy = cbuffer.data(gk_geom_10_off + 579 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxyz = cbuffer.data(gk_geom_10_off + 580 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxxzz = cbuffer.data(gk_geom_10_off + 581 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxyyy = cbuffer.data(gk_geom_10_off + 582 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxyyz = cbuffer.data(gk_geom_10_off + 583 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxyzz = cbuffer.data(gk_geom_10_off + 584 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxxzzz = cbuffer.data(gk_geom_10_off + 585 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyyyy = cbuffer.data(gk_geom_10_off + 586 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyyyz = cbuffer.data(gk_geom_10_off + 587 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyyzz = cbuffer.data(gk_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxyzzz = cbuffer.data(gk_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_xxxy_xxxzzzz = cbuffer.data(gk_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyyyy = cbuffer.data(gk_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyyyz = cbuffer.data(gk_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyyzz = cbuffer.data(gk_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyyzzz = cbuffer.data(gk_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_xxxy_xxyzzzz = cbuffer.data(gk_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_xxxy_xxzzzzz = cbuffer.data(gk_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyyyy = cbuffer.data(gk_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyyyz = cbuffer.data(gk_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyyzz = cbuffer.data(gk_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyyzzz = cbuffer.data(gk_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_xxxy_xyyzzzz = cbuffer.data(gk_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_xxxy_xyzzzzz = cbuffer.data(gk_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_xxxy_xzzzzzz = cbuffer.data(gk_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxxx = cbuffer.data(gk_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxxy = cbuffer.data(gk_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxxz = cbuffer.data(gk_geom_10_off + 614 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxyy = cbuffer.data(gk_geom_10_off + 615 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxyz = cbuffer.data(gk_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxxzz = cbuffer.data(gk_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxyyy = cbuffer.data(gk_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxyyz = cbuffer.data(gk_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxyzz = cbuffer.data(gk_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxxzzz = cbuffer.data(gk_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyyyy = cbuffer.data(gk_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyyyz = cbuffer.data(gk_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyyzz = cbuffer.data(gk_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxyzzz = cbuffer.data(gk_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_xxxz_xxxzzzz = cbuffer.data(gk_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyyyy = cbuffer.data(gk_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyyyz = cbuffer.data(gk_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyyzz = cbuffer.data(gk_geom_10_off + 629 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyyzzz = cbuffer.data(gk_geom_10_off + 630 * ccomps * dcomps);

            auto g_y_0_xxxz_xxyzzzz = cbuffer.data(gk_geom_10_off + 631 * ccomps * dcomps);

            auto g_y_0_xxxz_xxzzzzz = cbuffer.data(gk_geom_10_off + 632 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyyyy = cbuffer.data(gk_geom_10_off + 633 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyyyz = cbuffer.data(gk_geom_10_off + 634 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyyzz = cbuffer.data(gk_geom_10_off + 635 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyyzzz = cbuffer.data(gk_geom_10_off + 636 * ccomps * dcomps);

            auto g_y_0_xxxz_xyyzzzz = cbuffer.data(gk_geom_10_off + 637 * ccomps * dcomps);

            auto g_y_0_xxxz_xyzzzzz = cbuffer.data(gk_geom_10_off + 638 * ccomps * dcomps);

            auto g_y_0_xxxz_xzzzzzz = cbuffer.data(gk_geom_10_off + 639 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 648 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 649 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 650 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 651 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 652 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 653 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 654 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 655 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 656 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 657 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 658 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 659 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 660 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 661 * ccomps * dcomps);

            auto g_y_0_xxyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 662 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 663 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 664 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 665 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 666 * ccomps * dcomps);

            auto g_y_0_xxyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 667 * ccomps * dcomps);

            auto g_y_0_xxyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 668 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 669 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 670 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 671 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 672 * ccomps * dcomps);

            auto g_y_0_xxyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 673 * ccomps * dcomps);

            auto g_y_0_xxyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 674 * ccomps * dcomps);

            auto g_y_0_xxyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 675 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 684 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 685 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 686 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 687 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 688 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 689 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 690 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 691 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 692 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 693 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 694 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 695 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 696 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 697 * ccomps * dcomps);

            auto g_y_0_xxyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 698 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 699 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 700 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 701 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 702 * ccomps * dcomps);

            auto g_y_0_xxyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 703 * ccomps * dcomps);

            auto g_y_0_xxyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 704 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 705 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 706 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 707 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 708 * ccomps * dcomps);

            auto g_y_0_xxyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 709 * ccomps * dcomps);

            auto g_y_0_xxyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 710 * ccomps * dcomps);

            auto g_y_0_xxyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 711 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 720 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 721 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 722 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 723 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 724 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 725 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 726 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 727 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 728 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 729 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 730 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 731 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 732 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 733 * ccomps * dcomps);

            auto g_y_0_xxzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 734 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 735 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 736 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 737 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 738 * ccomps * dcomps);

            auto g_y_0_xxzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 739 * ccomps * dcomps);

            auto g_y_0_xxzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 740 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 741 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 742 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 743 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 744 * ccomps * dcomps);

            auto g_y_0_xxzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 745 * ccomps * dcomps);

            auto g_y_0_xxzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 746 * ccomps * dcomps);

            auto g_y_0_xxzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 747 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 756 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 757 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 758 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 759 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 760 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 761 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 762 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 763 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 764 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 765 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 766 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 767 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 768 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 769 * ccomps * dcomps);

            auto g_y_0_xyyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 770 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 771 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 772 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 773 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 774 * ccomps * dcomps);

            auto g_y_0_xyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 775 * ccomps * dcomps);

            auto g_y_0_xyyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 776 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 777 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 778 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 779 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 780 * ccomps * dcomps);

            auto g_y_0_xyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 781 * ccomps * dcomps);

            auto g_y_0_xyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 782 * ccomps * dcomps);

            auto g_y_0_xyyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 783 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 794 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 797 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 805 * ccomps * dcomps);

            auto g_y_0_xyyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 809 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_xyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 811 * ccomps * dcomps);

            auto g_y_0_xyyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_xyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_xyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 818 * ccomps * dcomps);

            auto g_y_0_xyyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 833 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 839 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 840 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 841 * ccomps * dcomps);

            auto g_y_0_xyzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 842 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 843 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 844 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 845 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 846 * ccomps * dcomps);

            auto g_y_0_xyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 847 * ccomps * dcomps);

            auto g_y_0_xyzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 848 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 849 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 850 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 851 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 852 * ccomps * dcomps);

            auto g_y_0_xyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 853 * ccomps * dcomps);

            auto g_y_0_xyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 854 * ccomps * dcomps);

            auto g_y_0_xyzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 855 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 864 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 865 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 866 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 867 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 868 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 869 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 870 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 871 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 872 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 873 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 874 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 875 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 876 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 877 * ccomps * dcomps);

            auto g_y_0_xzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 878 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 879 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 880 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 881 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 882 * ccomps * dcomps);

            auto g_y_0_xzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 883 * ccomps * dcomps);

            auto g_y_0_xzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 884 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 885 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 886 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 887 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 888 * ccomps * dcomps);

            auto g_y_0_xzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 889 * ccomps * dcomps);

            auto g_y_0_xzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 890 * ccomps * dcomps);

            auto g_y_0_xzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 891 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 900 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 901 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 902 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 903 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 904 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 905 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 906 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 907 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 908 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 909 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 910 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 911 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 912 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 913 * ccomps * dcomps);

            auto g_y_0_yyyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 914 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 915 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 916 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 917 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 918 * ccomps * dcomps);

            auto g_y_0_yyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 919 * ccomps * dcomps);

            auto g_y_0_yyyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 920 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 921 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 922 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 923 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 924 * ccomps * dcomps);

            auto g_y_0_yyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 925 * ccomps * dcomps);

            auto g_y_0_yyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 926 * ccomps * dcomps);

            auto g_y_0_yyyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 927 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 928 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 929 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 930 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 931 * ccomps * dcomps);

            auto g_y_0_yyyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 932 * ccomps * dcomps);

            auto g_y_0_yyyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 933 * ccomps * dcomps);

            auto g_y_0_yyyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 934 * ccomps * dcomps);

            auto g_y_0_yyyy_zzzzzzz = cbuffer.data(gk_geom_10_off + 935 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 936 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 937 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 938 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 939 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 940 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 941 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 942 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 943 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 944 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 945 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 946 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 947 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 948 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 949 * ccomps * dcomps);

            auto g_y_0_yyyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 950 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 951 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 952 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 953 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 954 * ccomps * dcomps);

            auto g_y_0_yyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 955 * ccomps * dcomps);

            auto g_y_0_yyyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 956 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 957 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 958 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 959 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 960 * ccomps * dcomps);

            auto g_y_0_yyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 961 * ccomps * dcomps);

            auto g_y_0_yyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 962 * ccomps * dcomps);

            auto g_y_0_yyyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 963 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 965 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 966 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 967 * ccomps * dcomps);

            auto g_y_0_yyyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 968 * ccomps * dcomps);

            auto g_y_0_yyyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 969 * ccomps * dcomps);

            auto g_y_0_yyyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 970 * ccomps * dcomps);

            auto g_y_0_yyyz_zzzzzzz = cbuffer.data(gk_geom_10_off + 971 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 972 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 973 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 974 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 975 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 976 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 977 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 978 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 979 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 980 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 981 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 982 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 983 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 984 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 985 * ccomps * dcomps);

            auto g_y_0_yyzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 986 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 987 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 988 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 989 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 990 * ccomps * dcomps);

            auto g_y_0_yyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 991 * ccomps * dcomps);

            auto g_y_0_yyzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 992 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 993 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 994 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 995 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 996 * ccomps * dcomps);

            auto g_y_0_yyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 997 * ccomps * dcomps);

            auto g_y_0_yyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 998 * ccomps * dcomps);

            auto g_y_0_yyzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 999 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1001 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1002 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1003 * ccomps * dcomps);

            auto g_y_0_yyzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1004 * ccomps * dcomps);

            auto g_y_0_yyzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1005 * ccomps * dcomps);

            auto g_y_0_yyzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1006 * ccomps * dcomps);

            auto g_y_0_yyzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1007 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1008 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1009 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1010 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1011 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1012 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1013 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1014 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1015 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1016 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1017 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1018 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1019 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1020 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1021 * ccomps * dcomps);

            auto g_y_0_yzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1022 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1023 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1024 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1025 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1026 * ccomps * dcomps);

            auto g_y_0_yzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1027 * ccomps * dcomps);

            auto g_y_0_yzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1028 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1029 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1030 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1031 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1032 * ccomps * dcomps);

            auto g_y_0_yzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1033 * ccomps * dcomps);

            auto g_y_0_yzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1034 * ccomps * dcomps);

            auto g_y_0_yzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1035 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1037 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1038 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1039 * ccomps * dcomps);

            auto g_y_0_yzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1040 * ccomps * dcomps);

            auto g_y_0_yzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1041 * ccomps * dcomps);

            auto g_y_0_yzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1042 * ccomps * dcomps);

            auto g_y_0_yzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1043 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1044 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1045 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1046 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1047 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1048 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1049 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1050 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1051 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1052 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1053 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1054 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1055 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1056 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1057 * ccomps * dcomps);

            auto g_y_0_zzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1058 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1059 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1060 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1061 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1062 * ccomps * dcomps);

            auto g_y_0_zzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1063 * ccomps * dcomps);

            auto g_y_0_zzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1064 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1065 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1066 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1067 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1068 * ccomps * dcomps);

            auto g_y_0_zzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1069 * ccomps * dcomps);

            auto g_y_0_zzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1070 * ccomps * dcomps);

            auto g_y_0_zzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1071 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1073 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1074 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1075 * ccomps * dcomps);

            auto g_y_0_zzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1076 * ccomps * dcomps);

            auto g_y_0_zzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1077 * ccomps * dcomps);

            auto g_y_0_zzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1078 * ccomps * dcomps);

            auto g_y_0_zzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1079 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxxx = cbuffer.data(gk_geom_10_off + 1080 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxxy = cbuffer.data(gk_geom_10_off + 1081 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxxz = cbuffer.data(gk_geom_10_off + 1082 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxyy = cbuffer.data(gk_geom_10_off + 1083 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxyz = cbuffer.data(gk_geom_10_off + 1084 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxxzz = cbuffer.data(gk_geom_10_off + 1085 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxyyy = cbuffer.data(gk_geom_10_off + 1086 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxyyz = cbuffer.data(gk_geom_10_off + 1087 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxyzz = cbuffer.data(gk_geom_10_off + 1088 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxxzzz = cbuffer.data(gk_geom_10_off + 1089 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyyyy = cbuffer.data(gk_geom_10_off + 1090 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyyyz = cbuffer.data(gk_geom_10_off + 1091 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyyzz = cbuffer.data(gk_geom_10_off + 1092 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxyzzz = cbuffer.data(gk_geom_10_off + 1093 * ccomps * dcomps);

            auto g_z_0_xxxx_xxxzzzz = cbuffer.data(gk_geom_10_off + 1094 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyyyy = cbuffer.data(gk_geom_10_off + 1095 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyyyz = cbuffer.data(gk_geom_10_off + 1096 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyyzz = cbuffer.data(gk_geom_10_off + 1097 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyyzzz = cbuffer.data(gk_geom_10_off + 1098 * ccomps * dcomps);

            auto g_z_0_xxxx_xxyzzzz = cbuffer.data(gk_geom_10_off + 1099 * ccomps * dcomps);

            auto g_z_0_xxxx_xxzzzzz = cbuffer.data(gk_geom_10_off + 1100 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyyyy = cbuffer.data(gk_geom_10_off + 1101 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyyyz = cbuffer.data(gk_geom_10_off + 1102 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyyzz = cbuffer.data(gk_geom_10_off + 1103 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyyzzz = cbuffer.data(gk_geom_10_off + 1104 * ccomps * dcomps);

            auto g_z_0_xxxx_xyyzzzz = cbuffer.data(gk_geom_10_off + 1105 * ccomps * dcomps);

            auto g_z_0_xxxx_xyzzzzz = cbuffer.data(gk_geom_10_off + 1106 * ccomps * dcomps);

            auto g_z_0_xxxx_xzzzzzz = cbuffer.data(gk_geom_10_off + 1107 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxxx = cbuffer.data(gk_geom_10_off + 1116 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxxy = cbuffer.data(gk_geom_10_off + 1117 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxxz = cbuffer.data(gk_geom_10_off + 1118 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxyy = cbuffer.data(gk_geom_10_off + 1119 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxyz = cbuffer.data(gk_geom_10_off + 1120 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxxzz = cbuffer.data(gk_geom_10_off + 1121 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxyyy = cbuffer.data(gk_geom_10_off + 1122 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxyyz = cbuffer.data(gk_geom_10_off + 1123 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxyzz = cbuffer.data(gk_geom_10_off + 1124 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxxzzz = cbuffer.data(gk_geom_10_off + 1125 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyyyy = cbuffer.data(gk_geom_10_off + 1126 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyyyz = cbuffer.data(gk_geom_10_off + 1127 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyyzz = cbuffer.data(gk_geom_10_off + 1128 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxyzzz = cbuffer.data(gk_geom_10_off + 1129 * ccomps * dcomps);

            auto g_z_0_xxxy_xxxzzzz = cbuffer.data(gk_geom_10_off + 1130 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyyyy = cbuffer.data(gk_geom_10_off + 1131 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyyyz = cbuffer.data(gk_geom_10_off + 1132 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyyzz = cbuffer.data(gk_geom_10_off + 1133 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyyzzz = cbuffer.data(gk_geom_10_off + 1134 * ccomps * dcomps);

            auto g_z_0_xxxy_xxyzzzz = cbuffer.data(gk_geom_10_off + 1135 * ccomps * dcomps);

            auto g_z_0_xxxy_xxzzzzz = cbuffer.data(gk_geom_10_off + 1136 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyyyy = cbuffer.data(gk_geom_10_off + 1137 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyyyz = cbuffer.data(gk_geom_10_off + 1138 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyyzz = cbuffer.data(gk_geom_10_off + 1139 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyyzzz = cbuffer.data(gk_geom_10_off + 1140 * ccomps * dcomps);

            auto g_z_0_xxxy_xyyzzzz = cbuffer.data(gk_geom_10_off + 1141 * ccomps * dcomps);

            auto g_z_0_xxxy_xyzzzzz = cbuffer.data(gk_geom_10_off + 1142 * ccomps * dcomps);

            auto g_z_0_xxxy_xzzzzzz = cbuffer.data(gk_geom_10_off + 1143 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1152 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1153 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1154 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1155 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1156 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1157 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1158 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1159 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1160 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1161 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1162 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1163 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1164 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1165 * ccomps * dcomps);

            auto g_z_0_xxxz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1166 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1167 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1168 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1169 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1170 * ccomps * dcomps);

            auto g_z_0_xxxz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1171 * ccomps * dcomps);

            auto g_z_0_xxxz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1172 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1173 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1174 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1175 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1176 * ccomps * dcomps);

            auto g_z_0_xxxz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1177 * ccomps * dcomps);

            auto g_z_0_xxxz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1178 * ccomps * dcomps);

            auto g_z_0_xxxz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1179 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 1188 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 1189 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 1190 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 1191 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 1192 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 1193 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 1194 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 1195 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 1196 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 1197 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 1198 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 1199 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 1200 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 1201 * ccomps * dcomps);

            auto g_z_0_xxyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 1202 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 1203 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 1204 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 1205 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 1206 * ccomps * dcomps);

            auto g_z_0_xxyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 1207 * ccomps * dcomps);

            auto g_z_0_xxyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 1208 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 1209 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 1210 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 1211 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 1212 * ccomps * dcomps);

            auto g_z_0_xxyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 1213 * ccomps * dcomps);

            auto g_z_0_xxyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 1214 * ccomps * dcomps);

            auto g_z_0_xxyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 1215 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1224 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1225 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1226 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1227 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1228 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1229 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1230 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1231 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1232 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1233 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1234 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1235 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1236 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1237 * ccomps * dcomps);

            auto g_z_0_xxyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1238 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1239 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1240 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1241 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1242 * ccomps * dcomps);

            auto g_z_0_xxyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1243 * ccomps * dcomps);

            auto g_z_0_xxyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1244 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1245 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1246 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1247 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1248 * ccomps * dcomps);

            auto g_z_0_xxyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1249 * ccomps * dcomps);

            auto g_z_0_xxyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1250 * ccomps * dcomps);

            auto g_z_0_xxyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1251 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1260 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1261 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1262 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1263 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1264 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1265 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1266 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1267 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1268 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1269 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1270 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1271 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1272 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1273 * ccomps * dcomps);

            auto g_z_0_xxzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1274 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1275 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1276 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1277 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1278 * ccomps * dcomps);

            auto g_z_0_xxzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1279 * ccomps * dcomps);

            auto g_z_0_xxzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1280 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1281 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1282 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1283 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1284 * ccomps * dcomps);

            auto g_z_0_xxzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1285 * ccomps * dcomps);

            auto g_z_0_xxzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1286 * ccomps * dcomps);

            auto g_z_0_xxzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1287 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 1296 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 1297 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 1298 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 1299 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 1300 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 1301 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 1302 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 1303 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 1304 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 1305 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 1306 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 1307 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 1308 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 1309 * ccomps * dcomps);

            auto g_z_0_xyyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 1310 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 1311 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 1312 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 1313 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 1314 * ccomps * dcomps);

            auto g_z_0_xyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 1315 * ccomps * dcomps);

            auto g_z_0_xyyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 1316 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 1317 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 1318 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 1319 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 1320 * ccomps * dcomps);

            auto g_z_0_xyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 1321 * ccomps * dcomps);

            auto g_z_0_xyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 1322 * ccomps * dcomps);

            auto g_z_0_xyyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 1323 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1332 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1333 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1334 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1335 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1336 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1337 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1338 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1339 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1340 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1341 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1342 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1343 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1344 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1345 * ccomps * dcomps);

            auto g_z_0_xyyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1346 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1347 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1348 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1349 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1350 * ccomps * dcomps);

            auto g_z_0_xyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1351 * ccomps * dcomps);

            auto g_z_0_xyyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1352 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1353 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1354 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1355 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1356 * ccomps * dcomps);

            auto g_z_0_xyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1357 * ccomps * dcomps);

            auto g_z_0_xyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1358 * ccomps * dcomps);

            auto g_z_0_xyyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1359 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1368 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1369 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1370 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1371 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1372 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1373 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1374 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1375 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1376 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1377 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1378 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1379 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1380 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1381 * ccomps * dcomps);

            auto g_z_0_xyzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1382 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1383 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1384 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1385 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1386 * ccomps * dcomps);

            auto g_z_0_xyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1387 * ccomps * dcomps);

            auto g_z_0_xyzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1388 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1389 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1390 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1391 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1392 * ccomps * dcomps);

            auto g_z_0_xyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1393 * ccomps * dcomps);

            auto g_z_0_xyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1394 * ccomps * dcomps);

            auto g_z_0_xyzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1395 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1404 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1405 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1406 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1407 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1408 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1409 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1410 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1411 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1412 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1413 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1414 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1415 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1416 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1417 * ccomps * dcomps);

            auto g_z_0_xzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1418 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1419 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1420 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1421 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1422 * ccomps * dcomps);

            auto g_z_0_xzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1423 * ccomps * dcomps);

            auto g_z_0_xzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1424 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1425 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1426 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1427 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1428 * ccomps * dcomps);

            auto g_z_0_xzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1429 * ccomps * dcomps);

            auto g_z_0_xzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1430 * ccomps * dcomps);

            auto g_z_0_xzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1431 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxxx = cbuffer.data(gk_geom_10_off + 1440 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxxy = cbuffer.data(gk_geom_10_off + 1441 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxxz = cbuffer.data(gk_geom_10_off + 1442 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxyy = cbuffer.data(gk_geom_10_off + 1443 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxyz = cbuffer.data(gk_geom_10_off + 1444 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxxzz = cbuffer.data(gk_geom_10_off + 1445 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxyyy = cbuffer.data(gk_geom_10_off + 1446 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxyyz = cbuffer.data(gk_geom_10_off + 1447 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxyzz = cbuffer.data(gk_geom_10_off + 1448 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxxzzz = cbuffer.data(gk_geom_10_off + 1449 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyyyy = cbuffer.data(gk_geom_10_off + 1450 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyyyz = cbuffer.data(gk_geom_10_off + 1451 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyyzz = cbuffer.data(gk_geom_10_off + 1452 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxyzzz = cbuffer.data(gk_geom_10_off + 1453 * ccomps * dcomps);

            auto g_z_0_yyyy_xxxzzzz = cbuffer.data(gk_geom_10_off + 1454 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyyyy = cbuffer.data(gk_geom_10_off + 1455 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyyyz = cbuffer.data(gk_geom_10_off + 1456 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyyzz = cbuffer.data(gk_geom_10_off + 1457 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyyzzz = cbuffer.data(gk_geom_10_off + 1458 * ccomps * dcomps);

            auto g_z_0_yyyy_xxyzzzz = cbuffer.data(gk_geom_10_off + 1459 * ccomps * dcomps);

            auto g_z_0_yyyy_xxzzzzz = cbuffer.data(gk_geom_10_off + 1460 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyyyy = cbuffer.data(gk_geom_10_off + 1461 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyyyz = cbuffer.data(gk_geom_10_off + 1462 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyyzz = cbuffer.data(gk_geom_10_off + 1463 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyyzzz = cbuffer.data(gk_geom_10_off + 1464 * ccomps * dcomps);

            auto g_z_0_yyyy_xyyzzzz = cbuffer.data(gk_geom_10_off + 1465 * ccomps * dcomps);

            auto g_z_0_yyyy_xyzzzzz = cbuffer.data(gk_geom_10_off + 1466 * ccomps * dcomps);

            auto g_z_0_yyyy_xzzzzzz = cbuffer.data(gk_geom_10_off + 1467 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyyyy = cbuffer.data(gk_geom_10_off + 1468 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyyyz = cbuffer.data(gk_geom_10_off + 1469 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyyzz = cbuffer.data(gk_geom_10_off + 1470 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyyzzz = cbuffer.data(gk_geom_10_off + 1471 * ccomps * dcomps);

            auto g_z_0_yyyy_yyyzzzz = cbuffer.data(gk_geom_10_off + 1472 * ccomps * dcomps);

            auto g_z_0_yyyy_yyzzzzz = cbuffer.data(gk_geom_10_off + 1473 * ccomps * dcomps);

            auto g_z_0_yyyy_yzzzzzz = cbuffer.data(gk_geom_10_off + 1474 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1476 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1477 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1478 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1479 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1480 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1481 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1482 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1483 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1484 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1485 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1486 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1487 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1488 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1489 * ccomps * dcomps);

            auto g_z_0_yyyz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1490 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1491 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1492 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1493 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1494 * ccomps * dcomps);

            auto g_z_0_yyyz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1495 * ccomps * dcomps);

            auto g_z_0_yyyz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1496 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1497 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1498 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1499 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1500 * ccomps * dcomps);

            auto g_z_0_yyyz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1501 * ccomps * dcomps);

            auto g_z_0_yyyz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1502 * ccomps * dcomps);

            auto g_z_0_yyyz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1503 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1504 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1505 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1506 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1507 * ccomps * dcomps);

            auto g_z_0_yyyz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1508 * ccomps * dcomps);

            auto g_z_0_yyyz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1509 * ccomps * dcomps);

            auto g_z_0_yyyz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1510 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1512 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1513 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1514 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1515 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1516 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1517 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1518 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1519 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1520 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1521 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1522 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1523 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1524 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1525 * ccomps * dcomps);

            auto g_z_0_yyzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1526 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1527 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1528 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1529 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1530 * ccomps * dcomps);

            auto g_z_0_yyzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1531 * ccomps * dcomps);

            auto g_z_0_yyzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1532 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1533 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1534 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1535 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1536 * ccomps * dcomps);

            auto g_z_0_yyzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1537 * ccomps * dcomps);

            auto g_z_0_yyzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1538 * ccomps * dcomps);

            auto g_z_0_yyzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1539 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1540 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1541 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1542 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1543 * ccomps * dcomps);

            auto g_z_0_yyzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1544 * ccomps * dcomps);

            auto g_z_0_yyzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1545 * ccomps * dcomps);

            auto g_z_0_yyzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1546 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1548 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1549 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1550 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1551 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1552 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1553 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1554 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1555 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1556 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1557 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1558 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1559 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1560 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1561 * ccomps * dcomps);

            auto g_z_0_yzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1562 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1563 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1564 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1565 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1566 * ccomps * dcomps);

            auto g_z_0_yzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1567 * ccomps * dcomps);

            auto g_z_0_yzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1568 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1569 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1570 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1571 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1572 * ccomps * dcomps);

            auto g_z_0_yzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1573 * ccomps * dcomps);

            auto g_z_0_yzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1574 * ccomps * dcomps);

            auto g_z_0_yzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1575 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1576 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1577 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1578 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1579 * ccomps * dcomps);

            auto g_z_0_yzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1580 * ccomps * dcomps);

            auto g_z_0_yzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1581 * ccomps * dcomps);

            auto g_z_0_yzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1582 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxxx = cbuffer.data(gk_geom_10_off + 1584 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxxy = cbuffer.data(gk_geom_10_off + 1585 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxxz = cbuffer.data(gk_geom_10_off + 1586 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxyy = cbuffer.data(gk_geom_10_off + 1587 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxyz = cbuffer.data(gk_geom_10_off + 1588 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxxzz = cbuffer.data(gk_geom_10_off + 1589 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxyyy = cbuffer.data(gk_geom_10_off + 1590 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxyyz = cbuffer.data(gk_geom_10_off + 1591 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxyzz = cbuffer.data(gk_geom_10_off + 1592 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxxzzz = cbuffer.data(gk_geom_10_off + 1593 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyyyy = cbuffer.data(gk_geom_10_off + 1594 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyyyz = cbuffer.data(gk_geom_10_off + 1595 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyyzz = cbuffer.data(gk_geom_10_off + 1596 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxyzzz = cbuffer.data(gk_geom_10_off + 1597 * ccomps * dcomps);

            auto g_z_0_zzzz_xxxzzzz = cbuffer.data(gk_geom_10_off + 1598 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyyyy = cbuffer.data(gk_geom_10_off + 1599 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyyyz = cbuffer.data(gk_geom_10_off + 1600 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyyzz = cbuffer.data(gk_geom_10_off + 1601 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyyzzz = cbuffer.data(gk_geom_10_off + 1602 * ccomps * dcomps);

            auto g_z_0_zzzz_xxyzzzz = cbuffer.data(gk_geom_10_off + 1603 * ccomps * dcomps);

            auto g_z_0_zzzz_xxzzzzz = cbuffer.data(gk_geom_10_off + 1604 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyyyy = cbuffer.data(gk_geom_10_off + 1605 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyyyz = cbuffer.data(gk_geom_10_off + 1606 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyyzz = cbuffer.data(gk_geom_10_off + 1607 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyyzzz = cbuffer.data(gk_geom_10_off + 1608 * ccomps * dcomps);

            auto g_z_0_zzzz_xyyzzzz = cbuffer.data(gk_geom_10_off + 1609 * ccomps * dcomps);

            auto g_z_0_zzzz_xyzzzzz = cbuffer.data(gk_geom_10_off + 1610 * ccomps * dcomps);

            auto g_z_0_zzzz_xzzzzzz = cbuffer.data(gk_geom_10_off + 1611 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyyyy = cbuffer.data(gk_geom_10_off + 1612 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyyyz = cbuffer.data(gk_geom_10_off + 1613 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyyzz = cbuffer.data(gk_geom_10_off + 1614 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyyzzz = cbuffer.data(gk_geom_10_off + 1615 * ccomps * dcomps);

            auto g_z_0_zzzz_yyyzzzz = cbuffer.data(gk_geom_10_off + 1616 * ccomps * dcomps);

            auto g_z_0_zzzz_yyzzzzz = cbuffer.data(gk_geom_10_off + 1617 * ccomps * dcomps);

            auto g_z_0_zzzz_yzzzzzz = cbuffer.data(gk_geom_10_off + 1618 * ccomps * dcomps);

            auto g_z_0_zzzz_zzzzzzz = cbuffer.data(gk_geom_10_off + 1619 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hixx

            const auto hi_geom_10_off = idx_geom_10_hixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyyy = cbuffer.data(hi_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyyz = cbuffer.data(hi_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyzz = cbuffer.data(hi_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyzzz = cbuffer.data(hi_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyzzzz = cbuffer.data(hi_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxx_yzzzzz = cbuffer.data(hi_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxx_zzzzzz = cbuffer.data(hi_geom_10_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxx_xxxxxx, g_x_0_xxxx_xxxxxxx, g_x_0_xxxx_xxxxxxy, g_x_0_xxxx_xxxxxxz, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxxyy, g_x_0_xxxx_xxxxxyz, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxxzz, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyyy, g_x_0_xxxx_xxxxyyz, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxyzz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxxzzz, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyyy, g_x_0_xxxx_xxxyyyz, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyyzz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxyzzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxxzzzz, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyyy, g_x_0_xxxx_xxyyyyz, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyyzz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyyzzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxyzzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xxzzzzz, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyyy, g_x_0_xxxx_xyyyyyz, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyyzz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyyzzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyyzzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xyzzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_xzzzzzz, g_x_0_xxxx_yyyyyy, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_zzzzzz, g_x_0_xxxxx_xxxxxx, g_x_0_xxxxx_xxxxxy, g_x_0_xxxxx_xxxxxz, g_x_0_xxxxx_xxxxyy, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxzz, g_x_0_xxxxx_xxxyyy, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxzzz, g_x_0_xxxxx_xxyyyy, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxzzzz, g_x_0_xxxxx_xyyyyy, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xzzzzz, g_x_0_xxxxx_yyyyyy, g_x_0_xxxxx_yyyyyz, g_x_0_xxxxx_yyyyzz, g_x_0_xxxxx_yyyzzz, g_x_0_xxxxx_yyzzzz, g_x_0_xxxxx_yzzzzz, g_x_0_xxxxx_zzzzzz, g_xxxx_xxxxxx, g_xxxx_xxxxxy, g_xxxx_xxxxxz, g_xxxx_xxxxyy, g_xxxx_xxxxyz, g_xxxx_xxxxzz, g_xxxx_xxxyyy, g_xxxx_xxxyyz, g_xxxx_xxxyzz, g_xxxx_xxxzzz, g_xxxx_xxyyyy, g_xxxx_xxyyyz, g_xxxx_xxyyzz, g_xxxx_xxyzzz, g_xxxx_xxzzzz, g_xxxx_xyyyyy, g_xxxx_xyyyyz, g_xxxx_xyyyzz, g_xxxx_xyyzzz, g_xxxx_xyzzzz, g_xxxx_xzzzzz, g_xxxx_yyyyyy, g_xxxx_yyyyyz, g_xxxx_yyyyzz, g_xxxx_yyyzzz, g_xxxx_yyzzzz, g_xxxx_yzzzzz, g_xxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxx_xxxxxx[k] = -g_xxxx_xxxxxx[k] - g_x_0_xxxx_xxxxxx[k] * ab_x + g_x_0_xxxx_xxxxxxx[k];

                g_x_0_xxxxx_xxxxxy[k] = -g_xxxx_xxxxxy[k] - g_x_0_xxxx_xxxxxy[k] * ab_x + g_x_0_xxxx_xxxxxxy[k];

                g_x_0_xxxxx_xxxxxz[k] = -g_xxxx_xxxxxz[k] - g_x_0_xxxx_xxxxxz[k] * ab_x + g_x_0_xxxx_xxxxxxz[k];

                g_x_0_xxxxx_xxxxyy[k] = -g_xxxx_xxxxyy[k] - g_x_0_xxxx_xxxxyy[k] * ab_x + g_x_0_xxxx_xxxxxyy[k];

                g_x_0_xxxxx_xxxxyz[k] = -g_xxxx_xxxxyz[k] - g_x_0_xxxx_xxxxyz[k] * ab_x + g_x_0_xxxx_xxxxxyz[k];

                g_x_0_xxxxx_xxxxzz[k] = -g_xxxx_xxxxzz[k] - g_x_0_xxxx_xxxxzz[k] * ab_x + g_x_0_xxxx_xxxxxzz[k];

                g_x_0_xxxxx_xxxyyy[k] = -g_xxxx_xxxyyy[k] - g_x_0_xxxx_xxxyyy[k] * ab_x + g_x_0_xxxx_xxxxyyy[k];

                g_x_0_xxxxx_xxxyyz[k] = -g_xxxx_xxxyyz[k] - g_x_0_xxxx_xxxyyz[k] * ab_x + g_x_0_xxxx_xxxxyyz[k];

                g_x_0_xxxxx_xxxyzz[k] = -g_xxxx_xxxyzz[k] - g_x_0_xxxx_xxxyzz[k] * ab_x + g_x_0_xxxx_xxxxyzz[k];

                g_x_0_xxxxx_xxxzzz[k] = -g_xxxx_xxxzzz[k] - g_x_0_xxxx_xxxzzz[k] * ab_x + g_x_0_xxxx_xxxxzzz[k];

                g_x_0_xxxxx_xxyyyy[k] = -g_xxxx_xxyyyy[k] - g_x_0_xxxx_xxyyyy[k] * ab_x + g_x_0_xxxx_xxxyyyy[k];

                g_x_0_xxxxx_xxyyyz[k] = -g_xxxx_xxyyyz[k] - g_x_0_xxxx_xxyyyz[k] * ab_x + g_x_0_xxxx_xxxyyyz[k];

                g_x_0_xxxxx_xxyyzz[k] = -g_xxxx_xxyyzz[k] - g_x_0_xxxx_xxyyzz[k] * ab_x + g_x_0_xxxx_xxxyyzz[k];

                g_x_0_xxxxx_xxyzzz[k] = -g_xxxx_xxyzzz[k] - g_x_0_xxxx_xxyzzz[k] * ab_x + g_x_0_xxxx_xxxyzzz[k];

                g_x_0_xxxxx_xxzzzz[k] = -g_xxxx_xxzzzz[k] - g_x_0_xxxx_xxzzzz[k] * ab_x + g_x_0_xxxx_xxxzzzz[k];

                g_x_0_xxxxx_xyyyyy[k] = -g_xxxx_xyyyyy[k] - g_x_0_xxxx_xyyyyy[k] * ab_x + g_x_0_xxxx_xxyyyyy[k];

                g_x_0_xxxxx_xyyyyz[k] = -g_xxxx_xyyyyz[k] - g_x_0_xxxx_xyyyyz[k] * ab_x + g_x_0_xxxx_xxyyyyz[k];

                g_x_0_xxxxx_xyyyzz[k] = -g_xxxx_xyyyzz[k] - g_x_0_xxxx_xyyyzz[k] * ab_x + g_x_0_xxxx_xxyyyzz[k];

                g_x_0_xxxxx_xyyzzz[k] = -g_xxxx_xyyzzz[k] - g_x_0_xxxx_xyyzzz[k] * ab_x + g_x_0_xxxx_xxyyzzz[k];

                g_x_0_xxxxx_xyzzzz[k] = -g_xxxx_xyzzzz[k] - g_x_0_xxxx_xyzzzz[k] * ab_x + g_x_0_xxxx_xxyzzzz[k];

                g_x_0_xxxxx_xzzzzz[k] = -g_xxxx_xzzzzz[k] - g_x_0_xxxx_xzzzzz[k] * ab_x + g_x_0_xxxx_xxzzzzz[k];

                g_x_0_xxxxx_yyyyyy[k] = -g_xxxx_yyyyyy[k] - g_x_0_xxxx_yyyyyy[k] * ab_x + g_x_0_xxxx_xyyyyyy[k];

                g_x_0_xxxxx_yyyyyz[k] = -g_xxxx_yyyyyz[k] - g_x_0_xxxx_yyyyyz[k] * ab_x + g_x_0_xxxx_xyyyyyz[k];

                g_x_0_xxxxx_yyyyzz[k] = -g_xxxx_yyyyzz[k] - g_x_0_xxxx_yyyyzz[k] * ab_x + g_x_0_xxxx_xyyyyzz[k];

                g_x_0_xxxxx_yyyzzz[k] = -g_xxxx_yyyzzz[k] - g_x_0_xxxx_yyyzzz[k] * ab_x + g_x_0_xxxx_xyyyzzz[k];

                g_x_0_xxxxx_yyzzzz[k] = -g_xxxx_yyzzzz[k] - g_x_0_xxxx_yyzzzz[k] * ab_x + g_x_0_xxxx_xyyzzzz[k];

                g_x_0_xxxxx_yzzzzz[k] = -g_xxxx_yzzzzz[k] - g_x_0_xxxx_yzzzzz[k] * ab_x + g_x_0_xxxx_xyzzzzz[k];

                g_x_0_xxxxx_zzzzzz[k] = -g_xxxx_zzzzzz[k] - g_x_0_xxxx_zzzzzz[k] * ab_x + g_x_0_xxxx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyyy = cbuffer.data(hi_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyyz = cbuffer.data(hi_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyzz = cbuffer.data(hi_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyzzz = cbuffer.data(hi_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyzzzz = cbuffer.data(hi_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxxy_yzzzzz = cbuffer.data(hi_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxxy_zzzzzz = cbuffer.data(hi_geom_10_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxx_xxxxxx, g_x_0_xxxx_xxxxxxy, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxxyy, g_x_0_xxxx_xxxxxyz, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyyy, g_x_0_xxxx_xxxxyyz, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxyzz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyyy, g_x_0_xxxx_xxxyyyz, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyyzz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxyzzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyyy, g_x_0_xxxx_xxyyyyz, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyyzz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyyzzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxyzzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyyy, g_x_0_xxxx_xyyyyyz, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyyzz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyyzzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyyzzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xyzzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_yyyyyy, g_x_0_xxxx_yyyyyyy, g_x_0_xxxx_yyyyyyz, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyyzz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyyzzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyyzzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yyzzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_yzzzzzz, g_x_0_xxxx_zzzzzz, g_x_0_xxxxy_xxxxxx, g_x_0_xxxxy_xxxxxy, g_x_0_xxxxy_xxxxxz, g_x_0_xxxxy_xxxxyy, g_x_0_xxxxy_xxxxyz, g_x_0_xxxxy_xxxxzz, g_x_0_xxxxy_xxxyyy, g_x_0_xxxxy_xxxyyz, g_x_0_xxxxy_xxxyzz, g_x_0_xxxxy_xxxzzz, g_x_0_xxxxy_xxyyyy, g_x_0_xxxxy_xxyyyz, g_x_0_xxxxy_xxyyzz, g_x_0_xxxxy_xxyzzz, g_x_0_xxxxy_xxzzzz, g_x_0_xxxxy_xyyyyy, g_x_0_xxxxy_xyyyyz, g_x_0_xxxxy_xyyyzz, g_x_0_xxxxy_xyyzzz, g_x_0_xxxxy_xyzzzz, g_x_0_xxxxy_xzzzzz, g_x_0_xxxxy_yyyyyy, g_x_0_xxxxy_yyyyyz, g_x_0_xxxxy_yyyyzz, g_x_0_xxxxy_yyyzzz, g_x_0_xxxxy_yyzzzz, g_x_0_xxxxy_yzzzzz, g_x_0_xxxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxy_xxxxxx[k] = -g_x_0_xxxx_xxxxxx[k] * ab_y + g_x_0_xxxx_xxxxxxy[k];

                g_x_0_xxxxy_xxxxxy[k] = -g_x_0_xxxx_xxxxxy[k] * ab_y + g_x_0_xxxx_xxxxxyy[k];

                g_x_0_xxxxy_xxxxxz[k] = -g_x_0_xxxx_xxxxxz[k] * ab_y + g_x_0_xxxx_xxxxxyz[k];

                g_x_0_xxxxy_xxxxyy[k] = -g_x_0_xxxx_xxxxyy[k] * ab_y + g_x_0_xxxx_xxxxyyy[k];

                g_x_0_xxxxy_xxxxyz[k] = -g_x_0_xxxx_xxxxyz[k] * ab_y + g_x_0_xxxx_xxxxyyz[k];

                g_x_0_xxxxy_xxxxzz[k] = -g_x_0_xxxx_xxxxzz[k] * ab_y + g_x_0_xxxx_xxxxyzz[k];

                g_x_0_xxxxy_xxxyyy[k] = -g_x_0_xxxx_xxxyyy[k] * ab_y + g_x_0_xxxx_xxxyyyy[k];

                g_x_0_xxxxy_xxxyyz[k] = -g_x_0_xxxx_xxxyyz[k] * ab_y + g_x_0_xxxx_xxxyyyz[k];

                g_x_0_xxxxy_xxxyzz[k] = -g_x_0_xxxx_xxxyzz[k] * ab_y + g_x_0_xxxx_xxxyyzz[k];

                g_x_0_xxxxy_xxxzzz[k] = -g_x_0_xxxx_xxxzzz[k] * ab_y + g_x_0_xxxx_xxxyzzz[k];

                g_x_0_xxxxy_xxyyyy[k] = -g_x_0_xxxx_xxyyyy[k] * ab_y + g_x_0_xxxx_xxyyyyy[k];

                g_x_0_xxxxy_xxyyyz[k] = -g_x_0_xxxx_xxyyyz[k] * ab_y + g_x_0_xxxx_xxyyyyz[k];

                g_x_0_xxxxy_xxyyzz[k] = -g_x_0_xxxx_xxyyzz[k] * ab_y + g_x_0_xxxx_xxyyyzz[k];

                g_x_0_xxxxy_xxyzzz[k] = -g_x_0_xxxx_xxyzzz[k] * ab_y + g_x_0_xxxx_xxyyzzz[k];

                g_x_0_xxxxy_xxzzzz[k] = -g_x_0_xxxx_xxzzzz[k] * ab_y + g_x_0_xxxx_xxyzzzz[k];

                g_x_0_xxxxy_xyyyyy[k] = -g_x_0_xxxx_xyyyyy[k] * ab_y + g_x_0_xxxx_xyyyyyy[k];

                g_x_0_xxxxy_xyyyyz[k] = -g_x_0_xxxx_xyyyyz[k] * ab_y + g_x_0_xxxx_xyyyyyz[k];

                g_x_0_xxxxy_xyyyzz[k] = -g_x_0_xxxx_xyyyzz[k] * ab_y + g_x_0_xxxx_xyyyyzz[k];

                g_x_0_xxxxy_xyyzzz[k] = -g_x_0_xxxx_xyyzzz[k] * ab_y + g_x_0_xxxx_xyyyzzz[k];

                g_x_0_xxxxy_xyzzzz[k] = -g_x_0_xxxx_xyzzzz[k] * ab_y + g_x_0_xxxx_xyyzzzz[k];

                g_x_0_xxxxy_xzzzzz[k] = -g_x_0_xxxx_xzzzzz[k] * ab_y + g_x_0_xxxx_xyzzzzz[k];

                g_x_0_xxxxy_yyyyyy[k] = -g_x_0_xxxx_yyyyyy[k] * ab_y + g_x_0_xxxx_yyyyyyy[k];

                g_x_0_xxxxy_yyyyyz[k] = -g_x_0_xxxx_yyyyyz[k] * ab_y + g_x_0_xxxx_yyyyyyz[k];

                g_x_0_xxxxy_yyyyzz[k] = -g_x_0_xxxx_yyyyzz[k] * ab_y + g_x_0_xxxx_yyyyyzz[k];

                g_x_0_xxxxy_yyyzzz[k] = -g_x_0_xxxx_yyyzzz[k] * ab_y + g_x_0_xxxx_yyyyzzz[k];

                g_x_0_xxxxy_yyzzzz[k] = -g_x_0_xxxx_yyzzzz[k] * ab_y + g_x_0_xxxx_yyyzzzz[k];

                g_x_0_xxxxy_yzzzzz[k] = -g_x_0_xxxx_yzzzzz[k] * ab_y + g_x_0_xxxx_yyzzzzz[k];

                g_x_0_xxxxy_zzzzzz[k] = -g_x_0_xxxx_zzzzzz[k] * ab_y + g_x_0_xxxx_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyyy = cbuffer.data(hi_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyyz = cbuffer.data(hi_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyzz = cbuffer.data(hi_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyzzz = cbuffer.data(hi_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyzzzz = cbuffer.data(hi_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxxz_yzzzzz = cbuffer.data(hi_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxxz_zzzzzz = cbuffer.data(hi_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxx_xxxxxx, g_x_0_xxxx_xxxxxxz, g_x_0_xxxx_xxxxxy, g_x_0_xxxx_xxxxxyz, g_x_0_xxxx_xxxxxz, g_x_0_xxxx_xxxxxzz, g_x_0_xxxx_xxxxyy, g_x_0_xxxx_xxxxyyz, g_x_0_xxxx_xxxxyz, g_x_0_xxxx_xxxxyzz, g_x_0_xxxx_xxxxzz, g_x_0_xxxx_xxxxzzz, g_x_0_xxxx_xxxyyy, g_x_0_xxxx_xxxyyyz, g_x_0_xxxx_xxxyyz, g_x_0_xxxx_xxxyyzz, g_x_0_xxxx_xxxyzz, g_x_0_xxxx_xxxyzzz, g_x_0_xxxx_xxxzzz, g_x_0_xxxx_xxxzzzz, g_x_0_xxxx_xxyyyy, g_x_0_xxxx_xxyyyyz, g_x_0_xxxx_xxyyyz, g_x_0_xxxx_xxyyyzz, g_x_0_xxxx_xxyyzz, g_x_0_xxxx_xxyyzzz, g_x_0_xxxx_xxyzzz, g_x_0_xxxx_xxyzzzz, g_x_0_xxxx_xxzzzz, g_x_0_xxxx_xxzzzzz, g_x_0_xxxx_xyyyyy, g_x_0_xxxx_xyyyyyz, g_x_0_xxxx_xyyyyz, g_x_0_xxxx_xyyyyzz, g_x_0_xxxx_xyyyzz, g_x_0_xxxx_xyyyzzz, g_x_0_xxxx_xyyzzz, g_x_0_xxxx_xyyzzzz, g_x_0_xxxx_xyzzzz, g_x_0_xxxx_xyzzzzz, g_x_0_xxxx_xzzzzz, g_x_0_xxxx_xzzzzzz, g_x_0_xxxx_yyyyyy, g_x_0_xxxx_yyyyyyz, g_x_0_xxxx_yyyyyz, g_x_0_xxxx_yyyyyzz, g_x_0_xxxx_yyyyzz, g_x_0_xxxx_yyyyzzz, g_x_0_xxxx_yyyzzz, g_x_0_xxxx_yyyzzzz, g_x_0_xxxx_yyzzzz, g_x_0_xxxx_yyzzzzz, g_x_0_xxxx_yzzzzz, g_x_0_xxxx_yzzzzzz, g_x_0_xxxx_zzzzzz, g_x_0_xxxx_zzzzzzz, g_x_0_xxxxz_xxxxxx, g_x_0_xxxxz_xxxxxy, g_x_0_xxxxz_xxxxxz, g_x_0_xxxxz_xxxxyy, g_x_0_xxxxz_xxxxyz, g_x_0_xxxxz_xxxxzz, g_x_0_xxxxz_xxxyyy, g_x_0_xxxxz_xxxyyz, g_x_0_xxxxz_xxxyzz, g_x_0_xxxxz_xxxzzz, g_x_0_xxxxz_xxyyyy, g_x_0_xxxxz_xxyyyz, g_x_0_xxxxz_xxyyzz, g_x_0_xxxxz_xxyzzz, g_x_0_xxxxz_xxzzzz, g_x_0_xxxxz_xyyyyy, g_x_0_xxxxz_xyyyyz, g_x_0_xxxxz_xyyyzz, g_x_0_xxxxz_xyyzzz, g_x_0_xxxxz_xyzzzz, g_x_0_xxxxz_xzzzzz, g_x_0_xxxxz_yyyyyy, g_x_0_xxxxz_yyyyyz, g_x_0_xxxxz_yyyyzz, g_x_0_xxxxz_yyyzzz, g_x_0_xxxxz_yyzzzz, g_x_0_xxxxz_yzzzzz, g_x_0_xxxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxz_xxxxxx[k] = -g_x_0_xxxx_xxxxxx[k] * ab_z + g_x_0_xxxx_xxxxxxz[k];

                g_x_0_xxxxz_xxxxxy[k] = -g_x_0_xxxx_xxxxxy[k] * ab_z + g_x_0_xxxx_xxxxxyz[k];

                g_x_0_xxxxz_xxxxxz[k] = -g_x_0_xxxx_xxxxxz[k] * ab_z + g_x_0_xxxx_xxxxxzz[k];

                g_x_0_xxxxz_xxxxyy[k] = -g_x_0_xxxx_xxxxyy[k] * ab_z + g_x_0_xxxx_xxxxyyz[k];

                g_x_0_xxxxz_xxxxyz[k] = -g_x_0_xxxx_xxxxyz[k] * ab_z + g_x_0_xxxx_xxxxyzz[k];

                g_x_0_xxxxz_xxxxzz[k] = -g_x_0_xxxx_xxxxzz[k] * ab_z + g_x_0_xxxx_xxxxzzz[k];

                g_x_0_xxxxz_xxxyyy[k] = -g_x_0_xxxx_xxxyyy[k] * ab_z + g_x_0_xxxx_xxxyyyz[k];

                g_x_0_xxxxz_xxxyyz[k] = -g_x_0_xxxx_xxxyyz[k] * ab_z + g_x_0_xxxx_xxxyyzz[k];

                g_x_0_xxxxz_xxxyzz[k] = -g_x_0_xxxx_xxxyzz[k] * ab_z + g_x_0_xxxx_xxxyzzz[k];

                g_x_0_xxxxz_xxxzzz[k] = -g_x_0_xxxx_xxxzzz[k] * ab_z + g_x_0_xxxx_xxxzzzz[k];

                g_x_0_xxxxz_xxyyyy[k] = -g_x_0_xxxx_xxyyyy[k] * ab_z + g_x_0_xxxx_xxyyyyz[k];

                g_x_0_xxxxz_xxyyyz[k] = -g_x_0_xxxx_xxyyyz[k] * ab_z + g_x_0_xxxx_xxyyyzz[k];

                g_x_0_xxxxz_xxyyzz[k] = -g_x_0_xxxx_xxyyzz[k] * ab_z + g_x_0_xxxx_xxyyzzz[k];

                g_x_0_xxxxz_xxyzzz[k] = -g_x_0_xxxx_xxyzzz[k] * ab_z + g_x_0_xxxx_xxyzzzz[k];

                g_x_0_xxxxz_xxzzzz[k] = -g_x_0_xxxx_xxzzzz[k] * ab_z + g_x_0_xxxx_xxzzzzz[k];

                g_x_0_xxxxz_xyyyyy[k] = -g_x_0_xxxx_xyyyyy[k] * ab_z + g_x_0_xxxx_xyyyyyz[k];

                g_x_0_xxxxz_xyyyyz[k] = -g_x_0_xxxx_xyyyyz[k] * ab_z + g_x_0_xxxx_xyyyyzz[k];

                g_x_0_xxxxz_xyyyzz[k] = -g_x_0_xxxx_xyyyzz[k] * ab_z + g_x_0_xxxx_xyyyzzz[k];

                g_x_0_xxxxz_xyyzzz[k] = -g_x_0_xxxx_xyyzzz[k] * ab_z + g_x_0_xxxx_xyyzzzz[k];

                g_x_0_xxxxz_xyzzzz[k] = -g_x_0_xxxx_xyzzzz[k] * ab_z + g_x_0_xxxx_xyzzzzz[k];

                g_x_0_xxxxz_xzzzzz[k] = -g_x_0_xxxx_xzzzzz[k] * ab_z + g_x_0_xxxx_xzzzzzz[k];

                g_x_0_xxxxz_yyyyyy[k] = -g_x_0_xxxx_yyyyyy[k] * ab_z + g_x_0_xxxx_yyyyyyz[k];

                g_x_0_xxxxz_yyyyyz[k] = -g_x_0_xxxx_yyyyyz[k] * ab_z + g_x_0_xxxx_yyyyyzz[k];

                g_x_0_xxxxz_yyyyzz[k] = -g_x_0_xxxx_yyyyzz[k] * ab_z + g_x_0_xxxx_yyyyzzz[k];

                g_x_0_xxxxz_yyyzzz[k] = -g_x_0_xxxx_yyyzzz[k] * ab_z + g_x_0_xxxx_yyyzzzz[k];

                g_x_0_xxxxz_yyzzzz[k] = -g_x_0_xxxx_yyzzzz[k] * ab_z + g_x_0_xxxx_yyzzzzz[k];

                g_x_0_xxxxz_yzzzzz[k] = -g_x_0_xxxx_yzzzzz[k] * ab_z + g_x_0_xxxx_yzzzzzz[k];

                g_x_0_xxxxz_zzzzzz[k] = -g_x_0_xxxx_zzzzzz[k] * ab_z + g_x_0_xxxx_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyyy = cbuffer.data(hi_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyyz = cbuffer.data(hi_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyzz = cbuffer.data(hi_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyzzz = cbuffer.data(hi_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyzzzz = cbuffer.data(hi_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxxyy_yzzzzz = cbuffer.data(hi_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxxyy_zzzzzz = cbuffer.data(hi_geom_10_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxy_xxxxxx, g_x_0_xxxy_xxxxxxy, g_x_0_xxxy_xxxxxy, g_x_0_xxxy_xxxxxyy, g_x_0_xxxy_xxxxxyz, g_x_0_xxxy_xxxxxz, g_x_0_xxxy_xxxxyy, g_x_0_xxxy_xxxxyyy, g_x_0_xxxy_xxxxyyz, g_x_0_xxxy_xxxxyz, g_x_0_xxxy_xxxxyzz, g_x_0_xxxy_xxxxzz, g_x_0_xxxy_xxxyyy, g_x_0_xxxy_xxxyyyy, g_x_0_xxxy_xxxyyyz, g_x_0_xxxy_xxxyyz, g_x_0_xxxy_xxxyyzz, g_x_0_xxxy_xxxyzz, g_x_0_xxxy_xxxyzzz, g_x_0_xxxy_xxxzzz, g_x_0_xxxy_xxyyyy, g_x_0_xxxy_xxyyyyy, g_x_0_xxxy_xxyyyyz, g_x_0_xxxy_xxyyyz, g_x_0_xxxy_xxyyyzz, g_x_0_xxxy_xxyyzz, g_x_0_xxxy_xxyyzzz, g_x_0_xxxy_xxyzzz, g_x_0_xxxy_xxyzzzz, g_x_0_xxxy_xxzzzz, g_x_0_xxxy_xyyyyy, g_x_0_xxxy_xyyyyyy, g_x_0_xxxy_xyyyyyz, g_x_0_xxxy_xyyyyz, g_x_0_xxxy_xyyyyzz, g_x_0_xxxy_xyyyzz, g_x_0_xxxy_xyyyzzz, g_x_0_xxxy_xyyzzz, g_x_0_xxxy_xyyzzzz, g_x_0_xxxy_xyzzzz, g_x_0_xxxy_xyzzzzz, g_x_0_xxxy_xzzzzz, g_x_0_xxxy_yyyyyy, g_x_0_xxxy_yyyyyyy, g_x_0_xxxy_yyyyyyz, g_x_0_xxxy_yyyyyz, g_x_0_xxxy_yyyyyzz, g_x_0_xxxy_yyyyzz, g_x_0_xxxy_yyyyzzz, g_x_0_xxxy_yyyzzz, g_x_0_xxxy_yyyzzzz, g_x_0_xxxy_yyzzzz, g_x_0_xxxy_yyzzzzz, g_x_0_xxxy_yzzzzz, g_x_0_xxxy_yzzzzzz, g_x_0_xxxy_zzzzzz, g_x_0_xxxyy_xxxxxx, g_x_0_xxxyy_xxxxxy, g_x_0_xxxyy_xxxxxz, g_x_0_xxxyy_xxxxyy, g_x_0_xxxyy_xxxxyz, g_x_0_xxxyy_xxxxzz, g_x_0_xxxyy_xxxyyy, g_x_0_xxxyy_xxxyyz, g_x_0_xxxyy_xxxyzz, g_x_0_xxxyy_xxxzzz, g_x_0_xxxyy_xxyyyy, g_x_0_xxxyy_xxyyyz, g_x_0_xxxyy_xxyyzz, g_x_0_xxxyy_xxyzzz, g_x_0_xxxyy_xxzzzz, g_x_0_xxxyy_xyyyyy, g_x_0_xxxyy_xyyyyz, g_x_0_xxxyy_xyyyzz, g_x_0_xxxyy_xyyzzz, g_x_0_xxxyy_xyzzzz, g_x_0_xxxyy_xzzzzz, g_x_0_xxxyy_yyyyyy, g_x_0_xxxyy_yyyyyz, g_x_0_xxxyy_yyyyzz, g_x_0_xxxyy_yyyzzz, g_x_0_xxxyy_yyzzzz, g_x_0_xxxyy_yzzzzz, g_x_0_xxxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyy_xxxxxx[k] = -g_x_0_xxxy_xxxxxx[k] * ab_y + g_x_0_xxxy_xxxxxxy[k];

                g_x_0_xxxyy_xxxxxy[k] = -g_x_0_xxxy_xxxxxy[k] * ab_y + g_x_0_xxxy_xxxxxyy[k];

                g_x_0_xxxyy_xxxxxz[k] = -g_x_0_xxxy_xxxxxz[k] * ab_y + g_x_0_xxxy_xxxxxyz[k];

                g_x_0_xxxyy_xxxxyy[k] = -g_x_0_xxxy_xxxxyy[k] * ab_y + g_x_0_xxxy_xxxxyyy[k];

                g_x_0_xxxyy_xxxxyz[k] = -g_x_0_xxxy_xxxxyz[k] * ab_y + g_x_0_xxxy_xxxxyyz[k];

                g_x_0_xxxyy_xxxxzz[k] = -g_x_0_xxxy_xxxxzz[k] * ab_y + g_x_0_xxxy_xxxxyzz[k];

                g_x_0_xxxyy_xxxyyy[k] = -g_x_0_xxxy_xxxyyy[k] * ab_y + g_x_0_xxxy_xxxyyyy[k];

                g_x_0_xxxyy_xxxyyz[k] = -g_x_0_xxxy_xxxyyz[k] * ab_y + g_x_0_xxxy_xxxyyyz[k];

                g_x_0_xxxyy_xxxyzz[k] = -g_x_0_xxxy_xxxyzz[k] * ab_y + g_x_0_xxxy_xxxyyzz[k];

                g_x_0_xxxyy_xxxzzz[k] = -g_x_0_xxxy_xxxzzz[k] * ab_y + g_x_0_xxxy_xxxyzzz[k];

                g_x_0_xxxyy_xxyyyy[k] = -g_x_0_xxxy_xxyyyy[k] * ab_y + g_x_0_xxxy_xxyyyyy[k];

                g_x_0_xxxyy_xxyyyz[k] = -g_x_0_xxxy_xxyyyz[k] * ab_y + g_x_0_xxxy_xxyyyyz[k];

                g_x_0_xxxyy_xxyyzz[k] = -g_x_0_xxxy_xxyyzz[k] * ab_y + g_x_0_xxxy_xxyyyzz[k];

                g_x_0_xxxyy_xxyzzz[k] = -g_x_0_xxxy_xxyzzz[k] * ab_y + g_x_0_xxxy_xxyyzzz[k];

                g_x_0_xxxyy_xxzzzz[k] = -g_x_0_xxxy_xxzzzz[k] * ab_y + g_x_0_xxxy_xxyzzzz[k];

                g_x_0_xxxyy_xyyyyy[k] = -g_x_0_xxxy_xyyyyy[k] * ab_y + g_x_0_xxxy_xyyyyyy[k];

                g_x_0_xxxyy_xyyyyz[k] = -g_x_0_xxxy_xyyyyz[k] * ab_y + g_x_0_xxxy_xyyyyyz[k];

                g_x_0_xxxyy_xyyyzz[k] = -g_x_0_xxxy_xyyyzz[k] * ab_y + g_x_0_xxxy_xyyyyzz[k];

                g_x_0_xxxyy_xyyzzz[k] = -g_x_0_xxxy_xyyzzz[k] * ab_y + g_x_0_xxxy_xyyyzzz[k];

                g_x_0_xxxyy_xyzzzz[k] = -g_x_0_xxxy_xyzzzz[k] * ab_y + g_x_0_xxxy_xyyzzzz[k];

                g_x_0_xxxyy_xzzzzz[k] = -g_x_0_xxxy_xzzzzz[k] * ab_y + g_x_0_xxxy_xyzzzzz[k];

                g_x_0_xxxyy_yyyyyy[k] = -g_x_0_xxxy_yyyyyy[k] * ab_y + g_x_0_xxxy_yyyyyyy[k];

                g_x_0_xxxyy_yyyyyz[k] = -g_x_0_xxxy_yyyyyz[k] * ab_y + g_x_0_xxxy_yyyyyyz[k];

                g_x_0_xxxyy_yyyyzz[k] = -g_x_0_xxxy_yyyyzz[k] * ab_y + g_x_0_xxxy_yyyyyzz[k];

                g_x_0_xxxyy_yyyzzz[k] = -g_x_0_xxxy_yyyzzz[k] * ab_y + g_x_0_xxxy_yyyyzzz[k];

                g_x_0_xxxyy_yyzzzz[k] = -g_x_0_xxxy_yyzzzz[k] * ab_y + g_x_0_xxxy_yyyzzzz[k];

                g_x_0_xxxyy_yzzzzz[k] = -g_x_0_xxxy_yzzzzz[k] * ab_y + g_x_0_xxxy_yyzzzzz[k];

                g_x_0_xxxyy_zzzzzz[k] = -g_x_0_xxxy_zzzzzz[k] * ab_y + g_x_0_xxxy_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyyy = cbuffer.data(hi_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyyz = cbuffer.data(hi_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyzz = cbuffer.data(hi_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyzzz = cbuffer.data(hi_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyzzzz = cbuffer.data(hi_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxxyz_yzzzzz = cbuffer.data(hi_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxxyz_zzzzzz = cbuffer.data(hi_geom_10_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyz_xxxxxx, g_x_0_xxxyz_xxxxxy, g_x_0_xxxyz_xxxxxz, g_x_0_xxxyz_xxxxyy, g_x_0_xxxyz_xxxxyz, g_x_0_xxxyz_xxxxzz, g_x_0_xxxyz_xxxyyy, g_x_0_xxxyz_xxxyyz, g_x_0_xxxyz_xxxyzz, g_x_0_xxxyz_xxxzzz, g_x_0_xxxyz_xxyyyy, g_x_0_xxxyz_xxyyyz, g_x_0_xxxyz_xxyyzz, g_x_0_xxxyz_xxyzzz, g_x_0_xxxyz_xxzzzz, g_x_0_xxxyz_xyyyyy, g_x_0_xxxyz_xyyyyz, g_x_0_xxxyz_xyyyzz, g_x_0_xxxyz_xyyzzz, g_x_0_xxxyz_xyzzzz, g_x_0_xxxyz_xzzzzz, g_x_0_xxxyz_yyyyyy, g_x_0_xxxyz_yyyyyz, g_x_0_xxxyz_yyyyzz, g_x_0_xxxyz_yyyzzz, g_x_0_xxxyz_yyzzzz, g_x_0_xxxyz_yzzzzz, g_x_0_xxxyz_zzzzzz, g_x_0_xxxz_xxxxxx, g_x_0_xxxz_xxxxxxy, g_x_0_xxxz_xxxxxy, g_x_0_xxxz_xxxxxyy, g_x_0_xxxz_xxxxxyz, g_x_0_xxxz_xxxxxz, g_x_0_xxxz_xxxxyy, g_x_0_xxxz_xxxxyyy, g_x_0_xxxz_xxxxyyz, g_x_0_xxxz_xxxxyz, g_x_0_xxxz_xxxxyzz, g_x_0_xxxz_xxxxzz, g_x_0_xxxz_xxxyyy, g_x_0_xxxz_xxxyyyy, g_x_0_xxxz_xxxyyyz, g_x_0_xxxz_xxxyyz, g_x_0_xxxz_xxxyyzz, g_x_0_xxxz_xxxyzz, g_x_0_xxxz_xxxyzzz, g_x_0_xxxz_xxxzzz, g_x_0_xxxz_xxyyyy, g_x_0_xxxz_xxyyyyy, g_x_0_xxxz_xxyyyyz, g_x_0_xxxz_xxyyyz, g_x_0_xxxz_xxyyyzz, g_x_0_xxxz_xxyyzz, g_x_0_xxxz_xxyyzzz, g_x_0_xxxz_xxyzzz, g_x_0_xxxz_xxyzzzz, g_x_0_xxxz_xxzzzz, g_x_0_xxxz_xyyyyy, g_x_0_xxxz_xyyyyyy, g_x_0_xxxz_xyyyyyz, g_x_0_xxxz_xyyyyz, g_x_0_xxxz_xyyyyzz, g_x_0_xxxz_xyyyzz, g_x_0_xxxz_xyyyzzz, g_x_0_xxxz_xyyzzz, g_x_0_xxxz_xyyzzzz, g_x_0_xxxz_xyzzzz, g_x_0_xxxz_xyzzzzz, g_x_0_xxxz_xzzzzz, g_x_0_xxxz_yyyyyy, g_x_0_xxxz_yyyyyyy, g_x_0_xxxz_yyyyyyz, g_x_0_xxxz_yyyyyz, g_x_0_xxxz_yyyyyzz, g_x_0_xxxz_yyyyzz, g_x_0_xxxz_yyyyzzz, g_x_0_xxxz_yyyzzz, g_x_0_xxxz_yyyzzzz, g_x_0_xxxz_yyzzzz, g_x_0_xxxz_yyzzzzz, g_x_0_xxxz_yzzzzz, g_x_0_xxxz_yzzzzzz, g_x_0_xxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyz_xxxxxx[k] = -g_x_0_xxxz_xxxxxx[k] * ab_y + g_x_0_xxxz_xxxxxxy[k];

                g_x_0_xxxyz_xxxxxy[k] = -g_x_0_xxxz_xxxxxy[k] * ab_y + g_x_0_xxxz_xxxxxyy[k];

                g_x_0_xxxyz_xxxxxz[k] = -g_x_0_xxxz_xxxxxz[k] * ab_y + g_x_0_xxxz_xxxxxyz[k];

                g_x_0_xxxyz_xxxxyy[k] = -g_x_0_xxxz_xxxxyy[k] * ab_y + g_x_0_xxxz_xxxxyyy[k];

                g_x_0_xxxyz_xxxxyz[k] = -g_x_0_xxxz_xxxxyz[k] * ab_y + g_x_0_xxxz_xxxxyyz[k];

                g_x_0_xxxyz_xxxxzz[k] = -g_x_0_xxxz_xxxxzz[k] * ab_y + g_x_0_xxxz_xxxxyzz[k];

                g_x_0_xxxyz_xxxyyy[k] = -g_x_0_xxxz_xxxyyy[k] * ab_y + g_x_0_xxxz_xxxyyyy[k];

                g_x_0_xxxyz_xxxyyz[k] = -g_x_0_xxxz_xxxyyz[k] * ab_y + g_x_0_xxxz_xxxyyyz[k];

                g_x_0_xxxyz_xxxyzz[k] = -g_x_0_xxxz_xxxyzz[k] * ab_y + g_x_0_xxxz_xxxyyzz[k];

                g_x_0_xxxyz_xxxzzz[k] = -g_x_0_xxxz_xxxzzz[k] * ab_y + g_x_0_xxxz_xxxyzzz[k];

                g_x_0_xxxyz_xxyyyy[k] = -g_x_0_xxxz_xxyyyy[k] * ab_y + g_x_0_xxxz_xxyyyyy[k];

                g_x_0_xxxyz_xxyyyz[k] = -g_x_0_xxxz_xxyyyz[k] * ab_y + g_x_0_xxxz_xxyyyyz[k];

                g_x_0_xxxyz_xxyyzz[k] = -g_x_0_xxxz_xxyyzz[k] * ab_y + g_x_0_xxxz_xxyyyzz[k];

                g_x_0_xxxyz_xxyzzz[k] = -g_x_0_xxxz_xxyzzz[k] * ab_y + g_x_0_xxxz_xxyyzzz[k];

                g_x_0_xxxyz_xxzzzz[k] = -g_x_0_xxxz_xxzzzz[k] * ab_y + g_x_0_xxxz_xxyzzzz[k];

                g_x_0_xxxyz_xyyyyy[k] = -g_x_0_xxxz_xyyyyy[k] * ab_y + g_x_0_xxxz_xyyyyyy[k];

                g_x_0_xxxyz_xyyyyz[k] = -g_x_0_xxxz_xyyyyz[k] * ab_y + g_x_0_xxxz_xyyyyyz[k];

                g_x_0_xxxyz_xyyyzz[k] = -g_x_0_xxxz_xyyyzz[k] * ab_y + g_x_0_xxxz_xyyyyzz[k];

                g_x_0_xxxyz_xyyzzz[k] = -g_x_0_xxxz_xyyzzz[k] * ab_y + g_x_0_xxxz_xyyyzzz[k];

                g_x_0_xxxyz_xyzzzz[k] = -g_x_0_xxxz_xyzzzz[k] * ab_y + g_x_0_xxxz_xyyzzzz[k];

                g_x_0_xxxyz_xzzzzz[k] = -g_x_0_xxxz_xzzzzz[k] * ab_y + g_x_0_xxxz_xyzzzzz[k];

                g_x_0_xxxyz_yyyyyy[k] = -g_x_0_xxxz_yyyyyy[k] * ab_y + g_x_0_xxxz_yyyyyyy[k];

                g_x_0_xxxyz_yyyyyz[k] = -g_x_0_xxxz_yyyyyz[k] * ab_y + g_x_0_xxxz_yyyyyyz[k];

                g_x_0_xxxyz_yyyyzz[k] = -g_x_0_xxxz_yyyyzz[k] * ab_y + g_x_0_xxxz_yyyyyzz[k];

                g_x_0_xxxyz_yyyzzz[k] = -g_x_0_xxxz_yyyzzz[k] * ab_y + g_x_0_xxxz_yyyyzzz[k];

                g_x_0_xxxyz_yyzzzz[k] = -g_x_0_xxxz_yyzzzz[k] * ab_y + g_x_0_xxxz_yyyzzzz[k];

                g_x_0_xxxyz_yzzzzz[k] = -g_x_0_xxxz_yzzzzz[k] * ab_y + g_x_0_xxxz_yyzzzzz[k];

                g_x_0_xxxyz_zzzzzz[k] = -g_x_0_xxxz_zzzzzz[k] * ab_y + g_x_0_xxxz_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyyy = cbuffer.data(hi_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyyz = cbuffer.data(hi_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyzz = cbuffer.data(hi_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyzzz = cbuffer.data(hi_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyzzzz = cbuffer.data(hi_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxxzz_yzzzzz = cbuffer.data(hi_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxxzz_zzzzzz = cbuffer.data(hi_geom_10_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxz_xxxxxx, g_x_0_xxxz_xxxxxxz, g_x_0_xxxz_xxxxxy, g_x_0_xxxz_xxxxxyz, g_x_0_xxxz_xxxxxz, g_x_0_xxxz_xxxxxzz, g_x_0_xxxz_xxxxyy, g_x_0_xxxz_xxxxyyz, g_x_0_xxxz_xxxxyz, g_x_0_xxxz_xxxxyzz, g_x_0_xxxz_xxxxzz, g_x_0_xxxz_xxxxzzz, g_x_0_xxxz_xxxyyy, g_x_0_xxxz_xxxyyyz, g_x_0_xxxz_xxxyyz, g_x_0_xxxz_xxxyyzz, g_x_0_xxxz_xxxyzz, g_x_0_xxxz_xxxyzzz, g_x_0_xxxz_xxxzzz, g_x_0_xxxz_xxxzzzz, g_x_0_xxxz_xxyyyy, g_x_0_xxxz_xxyyyyz, g_x_0_xxxz_xxyyyz, g_x_0_xxxz_xxyyyzz, g_x_0_xxxz_xxyyzz, g_x_0_xxxz_xxyyzzz, g_x_0_xxxz_xxyzzz, g_x_0_xxxz_xxyzzzz, g_x_0_xxxz_xxzzzz, g_x_0_xxxz_xxzzzzz, g_x_0_xxxz_xyyyyy, g_x_0_xxxz_xyyyyyz, g_x_0_xxxz_xyyyyz, g_x_0_xxxz_xyyyyzz, g_x_0_xxxz_xyyyzz, g_x_0_xxxz_xyyyzzz, g_x_0_xxxz_xyyzzz, g_x_0_xxxz_xyyzzzz, g_x_0_xxxz_xyzzzz, g_x_0_xxxz_xyzzzzz, g_x_0_xxxz_xzzzzz, g_x_0_xxxz_xzzzzzz, g_x_0_xxxz_yyyyyy, g_x_0_xxxz_yyyyyyz, g_x_0_xxxz_yyyyyz, g_x_0_xxxz_yyyyyzz, g_x_0_xxxz_yyyyzz, g_x_0_xxxz_yyyyzzz, g_x_0_xxxz_yyyzzz, g_x_0_xxxz_yyyzzzz, g_x_0_xxxz_yyzzzz, g_x_0_xxxz_yyzzzzz, g_x_0_xxxz_yzzzzz, g_x_0_xxxz_yzzzzzz, g_x_0_xxxz_zzzzzz, g_x_0_xxxz_zzzzzzz, g_x_0_xxxzz_xxxxxx, g_x_0_xxxzz_xxxxxy, g_x_0_xxxzz_xxxxxz, g_x_0_xxxzz_xxxxyy, g_x_0_xxxzz_xxxxyz, g_x_0_xxxzz_xxxxzz, g_x_0_xxxzz_xxxyyy, g_x_0_xxxzz_xxxyyz, g_x_0_xxxzz_xxxyzz, g_x_0_xxxzz_xxxzzz, g_x_0_xxxzz_xxyyyy, g_x_0_xxxzz_xxyyyz, g_x_0_xxxzz_xxyyzz, g_x_0_xxxzz_xxyzzz, g_x_0_xxxzz_xxzzzz, g_x_0_xxxzz_xyyyyy, g_x_0_xxxzz_xyyyyz, g_x_0_xxxzz_xyyyzz, g_x_0_xxxzz_xyyzzz, g_x_0_xxxzz_xyzzzz, g_x_0_xxxzz_xzzzzz, g_x_0_xxxzz_yyyyyy, g_x_0_xxxzz_yyyyyz, g_x_0_xxxzz_yyyyzz, g_x_0_xxxzz_yyyzzz, g_x_0_xxxzz_yyzzzz, g_x_0_xxxzz_yzzzzz, g_x_0_xxxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzz_xxxxxx[k] = -g_x_0_xxxz_xxxxxx[k] * ab_z + g_x_0_xxxz_xxxxxxz[k];

                g_x_0_xxxzz_xxxxxy[k] = -g_x_0_xxxz_xxxxxy[k] * ab_z + g_x_0_xxxz_xxxxxyz[k];

                g_x_0_xxxzz_xxxxxz[k] = -g_x_0_xxxz_xxxxxz[k] * ab_z + g_x_0_xxxz_xxxxxzz[k];

                g_x_0_xxxzz_xxxxyy[k] = -g_x_0_xxxz_xxxxyy[k] * ab_z + g_x_0_xxxz_xxxxyyz[k];

                g_x_0_xxxzz_xxxxyz[k] = -g_x_0_xxxz_xxxxyz[k] * ab_z + g_x_0_xxxz_xxxxyzz[k];

                g_x_0_xxxzz_xxxxzz[k] = -g_x_0_xxxz_xxxxzz[k] * ab_z + g_x_0_xxxz_xxxxzzz[k];

                g_x_0_xxxzz_xxxyyy[k] = -g_x_0_xxxz_xxxyyy[k] * ab_z + g_x_0_xxxz_xxxyyyz[k];

                g_x_0_xxxzz_xxxyyz[k] = -g_x_0_xxxz_xxxyyz[k] * ab_z + g_x_0_xxxz_xxxyyzz[k];

                g_x_0_xxxzz_xxxyzz[k] = -g_x_0_xxxz_xxxyzz[k] * ab_z + g_x_0_xxxz_xxxyzzz[k];

                g_x_0_xxxzz_xxxzzz[k] = -g_x_0_xxxz_xxxzzz[k] * ab_z + g_x_0_xxxz_xxxzzzz[k];

                g_x_0_xxxzz_xxyyyy[k] = -g_x_0_xxxz_xxyyyy[k] * ab_z + g_x_0_xxxz_xxyyyyz[k];

                g_x_0_xxxzz_xxyyyz[k] = -g_x_0_xxxz_xxyyyz[k] * ab_z + g_x_0_xxxz_xxyyyzz[k];

                g_x_0_xxxzz_xxyyzz[k] = -g_x_0_xxxz_xxyyzz[k] * ab_z + g_x_0_xxxz_xxyyzzz[k];

                g_x_0_xxxzz_xxyzzz[k] = -g_x_0_xxxz_xxyzzz[k] * ab_z + g_x_0_xxxz_xxyzzzz[k];

                g_x_0_xxxzz_xxzzzz[k] = -g_x_0_xxxz_xxzzzz[k] * ab_z + g_x_0_xxxz_xxzzzzz[k];

                g_x_0_xxxzz_xyyyyy[k] = -g_x_0_xxxz_xyyyyy[k] * ab_z + g_x_0_xxxz_xyyyyyz[k];

                g_x_0_xxxzz_xyyyyz[k] = -g_x_0_xxxz_xyyyyz[k] * ab_z + g_x_0_xxxz_xyyyyzz[k];

                g_x_0_xxxzz_xyyyzz[k] = -g_x_0_xxxz_xyyyzz[k] * ab_z + g_x_0_xxxz_xyyyzzz[k];

                g_x_0_xxxzz_xyyzzz[k] = -g_x_0_xxxz_xyyzzz[k] * ab_z + g_x_0_xxxz_xyyzzzz[k];

                g_x_0_xxxzz_xyzzzz[k] = -g_x_0_xxxz_xyzzzz[k] * ab_z + g_x_0_xxxz_xyzzzzz[k];

                g_x_0_xxxzz_xzzzzz[k] = -g_x_0_xxxz_xzzzzz[k] * ab_z + g_x_0_xxxz_xzzzzzz[k];

                g_x_0_xxxzz_yyyyyy[k] = -g_x_0_xxxz_yyyyyy[k] * ab_z + g_x_0_xxxz_yyyyyyz[k];

                g_x_0_xxxzz_yyyyyz[k] = -g_x_0_xxxz_yyyyyz[k] * ab_z + g_x_0_xxxz_yyyyyzz[k];

                g_x_0_xxxzz_yyyyzz[k] = -g_x_0_xxxz_yyyyzz[k] * ab_z + g_x_0_xxxz_yyyyzzz[k];

                g_x_0_xxxzz_yyyzzz[k] = -g_x_0_xxxz_yyyzzz[k] * ab_z + g_x_0_xxxz_yyyzzzz[k];

                g_x_0_xxxzz_yyzzzz[k] = -g_x_0_xxxz_yyzzzz[k] * ab_z + g_x_0_xxxz_yyzzzzz[k];

                g_x_0_xxxzz_yzzzzz[k] = -g_x_0_xxxz_yzzzzz[k] * ab_z + g_x_0_xxxz_yzzzzzz[k];

                g_x_0_xxxzz_zzzzzz[k] = -g_x_0_xxxz_zzzzzz[k] * ab_z + g_x_0_xxxz_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xxyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xxyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyy_xxxxxx, g_x_0_xxyy_xxxxxxy, g_x_0_xxyy_xxxxxy, g_x_0_xxyy_xxxxxyy, g_x_0_xxyy_xxxxxyz, g_x_0_xxyy_xxxxxz, g_x_0_xxyy_xxxxyy, g_x_0_xxyy_xxxxyyy, g_x_0_xxyy_xxxxyyz, g_x_0_xxyy_xxxxyz, g_x_0_xxyy_xxxxyzz, g_x_0_xxyy_xxxxzz, g_x_0_xxyy_xxxyyy, g_x_0_xxyy_xxxyyyy, g_x_0_xxyy_xxxyyyz, g_x_0_xxyy_xxxyyz, g_x_0_xxyy_xxxyyzz, g_x_0_xxyy_xxxyzz, g_x_0_xxyy_xxxyzzz, g_x_0_xxyy_xxxzzz, g_x_0_xxyy_xxyyyy, g_x_0_xxyy_xxyyyyy, g_x_0_xxyy_xxyyyyz, g_x_0_xxyy_xxyyyz, g_x_0_xxyy_xxyyyzz, g_x_0_xxyy_xxyyzz, g_x_0_xxyy_xxyyzzz, g_x_0_xxyy_xxyzzz, g_x_0_xxyy_xxyzzzz, g_x_0_xxyy_xxzzzz, g_x_0_xxyy_xyyyyy, g_x_0_xxyy_xyyyyyy, g_x_0_xxyy_xyyyyyz, g_x_0_xxyy_xyyyyz, g_x_0_xxyy_xyyyyzz, g_x_0_xxyy_xyyyzz, g_x_0_xxyy_xyyyzzz, g_x_0_xxyy_xyyzzz, g_x_0_xxyy_xyyzzzz, g_x_0_xxyy_xyzzzz, g_x_0_xxyy_xyzzzzz, g_x_0_xxyy_xzzzzz, g_x_0_xxyy_yyyyyy, g_x_0_xxyy_yyyyyyy, g_x_0_xxyy_yyyyyyz, g_x_0_xxyy_yyyyyz, g_x_0_xxyy_yyyyyzz, g_x_0_xxyy_yyyyzz, g_x_0_xxyy_yyyyzzz, g_x_0_xxyy_yyyzzz, g_x_0_xxyy_yyyzzzz, g_x_0_xxyy_yyzzzz, g_x_0_xxyy_yyzzzzz, g_x_0_xxyy_yzzzzz, g_x_0_xxyy_yzzzzzz, g_x_0_xxyy_zzzzzz, g_x_0_xxyyy_xxxxxx, g_x_0_xxyyy_xxxxxy, g_x_0_xxyyy_xxxxxz, g_x_0_xxyyy_xxxxyy, g_x_0_xxyyy_xxxxyz, g_x_0_xxyyy_xxxxzz, g_x_0_xxyyy_xxxyyy, g_x_0_xxyyy_xxxyyz, g_x_0_xxyyy_xxxyzz, g_x_0_xxyyy_xxxzzz, g_x_0_xxyyy_xxyyyy, g_x_0_xxyyy_xxyyyz, g_x_0_xxyyy_xxyyzz, g_x_0_xxyyy_xxyzzz, g_x_0_xxyyy_xxzzzz, g_x_0_xxyyy_xyyyyy, g_x_0_xxyyy_xyyyyz, g_x_0_xxyyy_xyyyzz, g_x_0_xxyyy_xyyzzz, g_x_0_xxyyy_xyzzzz, g_x_0_xxyyy_xzzzzz, g_x_0_xxyyy_yyyyyy, g_x_0_xxyyy_yyyyyz, g_x_0_xxyyy_yyyyzz, g_x_0_xxyyy_yyyzzz, g_x_0_xxyyy_yyzzzz, g_x_0_xxyyy_yzzzzz, g_x_0_xxyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyy_xxxxxx[k] = -g_x_0_xxyy_xxxxxx[k] * ab_y + g_x_0_xxyy_xxxxxxy[k];

                g_x_0_xxyyy_xxxxxy[k] = -g_x_0_xxyy_xxxxxy[k] * ab_y + g_x_0_xxyy_xxxxxyy[k];

                g_x_0_xxyyy_xxxxxz[k] = -g_x_0_xxyy_xxxxxz[k] * ab_y + g_x_0_xxyy_xxxxxyz[k];

                g_x_0_xxyyy_xxxxyy[k] = -g_x_0_xxyy_xxxxyy[k] * ab_y + g_x_0_xxyy_xxxxyyy[k];

                g_x_0_xxyyy_xxxxyz[k] = -g_x_0_xxyy_xxxxyz[k] * ab_y + g_x_0_xxyy_xxxxyyz[k];

                g_x_0_xxyyy_xxxxzz[k] = -g_x_0_xxyy_xxxxzz[k] * ab_y + g_x_0_xxyy_xxxxyzz[k];

                g_x_0_xxyyy_xxxyyy[k] = -g_x_0_xxyy_xxxyyy[k] * ab_y + g_x_0_xxyy_xxxyyyy[k];

                g_x_0_xxyyy_xxxyyz[k] = -g_x_0_xxyy_xxxyyz[k] * ab_y + g_x_0_xxyy_xxxyyyz[k];

                g_x_0_xxyyy_xxxyzz[k] = -g_x_0_xxyy_xxxyzz[k] * ab_y + g_x_0_xxyy_xxxyyzz[k];

                g_x_0_xxyyy_xxxzzz[k] = -g_x_0_xxyy_xxxzzz[k] * ab_y + g_x_0_xxyy_xxxyzzz[k];

                g_x_0_xxyyy_xxyyyy[k] = -g_x_0_xxyy_xxyyyy[k] * ab_y + g_x_0_xxyy_xxyyyyy[k];

                g_x_0_xxyyy_xxyyyz[k] = -g_x_0_xxyy_xxyyyz[k] * ab_y + g_x_0_xxyy_xxyyyyz[k];

                g_x_0_xxyyy_xxyyzz[k] = -g_x_0_xxyy_xxyyzz[k] * ab_y + g_x_0_xxyy_xxyyyzz[k];

                g_x_0_xxyyy_xxyzzz[k] = -g_x_0_xxyy_xxyzzz[k] * ab_y + g_x_0_xxyy_xxyyzzz[k];

                g_x_0_xxyyy_xxzzzz[k] = -g_x_0_xxyy_xxzzzz[k] * ab_y + g_x_0_xxyy_xxyzzzz[k];

                g_x_0_xxyyy_xyyyyy[k] = -g_x_0_xxyy_xyyyyy[k] * ab_y + g_x_0_xxyy_xyyyyyy[k];

                g_x_0_xxyyy_xyyyyz[k] = -g_x_0_xxyy_xyyyyz[k] * ab_y + g_x_0_xxyy_xyyyyyz[k];

                g_x_0_xxyyy_xyyyzz[k] = -g_x_0_xxyy_xyyyzz[k] * ab_y + g_x_0_xxyy_xyyyyzz[k];

                g_x_0_xxyyy_xyyzzz[k] = -g_x_0_xxyy_xyyzzz[k] * ab_y + g_x_0_xxyy_xyyyzzz[k];

                g_x_0_xxyyy_xyzzzz[k] = -g_x_0_xxyy_xyzzzz[k] * ab_y + g_x_0_xxyy_xyyzzzz[k];

                g_x_0_xxyyy_xzzzzz[k] = -g_x_0_xxyy_xzzzzz[k] * ab_y + g_x_0_xxyy_xyzzzzz[k];

                g_x_0_xxyyy_yyyyyy[k] = -g_x_0_xxyy_yyyyyy[k] * ab_y + g_x_0_xxyy_yyyyyyy[k];

                g_x_0_xxyyy_yyyyyz[k] = -g_x_0_xxyy_yyyyyz[k] * ab_y + g_x_0_xxyy_yyyyyyz[k];

                g_x_0_xxyyy_yyyyzz[k] = -g_x_0_xxyy_yyyyzz[k] * ab_y + g_x_0_xxyy_yyyyyzz[k];

                g_x_0_xxyyy_yyyzzz[k] = -g_x_0_xxyy_yyyzzz[k] * ab_y + g_x_0_xxyy_yyyyzzz[k];

                g_x_0_xxyyy_yyzzzz[k] = -g_x_0_xxyy_yyzzzz[k] * ab_y + g_x_0_xxyy_yyyzzzz[k];

                g_x_0_xxyyy_yzzzzz[k] = -g_x_0_xxyy_yzzzzz[k] * ab_y + g_x_0_xxyy_yyzzzzz[k];

                g_x_0_xxyyy_zzzzzz[k] = -g_x_0_xxyy_zzzzzz[k] * ab_y + g_x_0_xxyy_yzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xxyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xxyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyz_xxxxxx, g_x_0_xxyyz_xxxxxy, g_x_0_xxyyz_xxxxxz, g_x_0_xxyyz_xxxxyy, g_x_0_xxyyz_xxxxyz, g_x_0_xxyyz_xxxxzz, g_x_0_xxyyz_xxxyyy, g_x_0_xxyyz_xxxyyz, g_x_0_xxyyz_xxxyzz, g_x_0_xxyyz_xxxzzz, g_x_0_xxyyz_xxyyyy, g_x_0_xxyyz_xxyyyz, g_x_0_xxyyz_xxyyzz, g_x_0_xxyyz_xxyzzz, g_x_0_xxyyz_xxzzzz, g_x_0_xxyyz_xyyyyy, g_x_0_xxyyz_xyyyyz, g_x_0_xxyyz_xyyyzz, g_x_0_xxyyz_xyyzzz, g_x_0_xxyyz_xyzzzz, g_x_0_xxyyz_xzzzzz, g_x_0_xxyyz_yyyyyy, g_x_0_xxyyz_yyyyyz, g_x_0_xxyyz_yyyyzz, g_x_0_xxyyz_yyyzzz, g_x_0_xxyyz_yyzzzz, g_x_0_xxyyz_yzzzzz, g_x_0_xxyyz_zzzzzz, g_x_0_xxyz_xxxxxx, g_x_0_xxyz_xxxxxxy, g_x_0_xxyz_xxxxxy, g_x_0_xxyz_xxxxxyy, g_x_0_xxyz_xxxxxyz, g_x_0_xxyz_xxxxxz, g_x_0_xxyz_xxxxyy, g_x_0_xxyz_xxxxyyy, g_x_0_xxyz_xxxxyyz, g_x_0_xxyz_xxxxyz, g_x_0_xxyz_xxxxyzz, g_x_0_xxyz_xxxxzz, g_x_0_xxyz_xxxyyy, g_x_0_xxyz_xxxyyyy, g_x_0_xxyz_xxxyyyz, g_x_0_xxyz_xxxyyz, g_x_0_xxyz_xxxyyzz, g_x_0_xxyz_xxxyzz, g_x_0_xxyz_xxxyzzz, g_x_0_xxyz_xxxzzz, g_x_0_xxyz_xxyyyy, g_x_0_xxyz_xxyyyyy, g_x_0_xxyz_xxyyyyz, g_x_0_xxyz_xxyyyz, g_x_0_xxyz_xxyyyzz, g_x_0_xxyz_xxyyzz, g_x_0_xxyz_xxyyzzz, g_x_0_xxyz_xxyzzz, g_x_0_xxyz_xxyzzzz, g_x_0_xxyz_xxzzzz, g_x_0_xxyz_xyyyyy, g_x_0_xxyz_xyyyyyy, g_x_0_xxyz_xyyyyyz, g_x_0_xxyz_xyyyyz, g_x_0_xxyz_xyyyyzz, g_x_0_xxyz_xyyyzz, g_x_0_xxyz_xyyyzzz, g_x_0_xxyz_xyyzzz, g_x_0_xxyz_xyyzzzz, g_x_0_xxyz_xyzzzz, g_x_0_xxyz_xyzzzzz, g_x_0_xxyz_xzzzzz, g_x_0_xxyz_yyyyyy, g_x_0_xxyz_yyyyyyy, g_x_0_xxyz_yyyyyyz, g_x_0_xxyz_yyyyyz, g_x_0_xxyz_yyyyyzz, g_x_0_xxyz_yyyyzz, g_x_0_xxyz_yyyyzzz, g_x_0_xxyz_yyyzzz, g_x_0_xxyz_yyyzzzz, g_x_0_xxyz_yyzzzz, g_x_0_xxyz_yyzzzzz, g_x_0_xxyz_yzzzzz, g_x_0_xxyz_yzzzzzz, g_x_0_xxyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyz_xxxxxx[k] = -g_x_0_xxyz_xxxxxx[k] * ab_y + g_x_0_xxyz_xxxxxxy[k];

                g_x_0_xxyyz_xxxxxy[k] = -g_x_0_xxyz_xxxxxy[k] * ab_y + g_x_0_xxyz_xxxxxyy[k];

                g_x_0_xxyyz_xxxxxz[k] = -g_x_0_xxyz_xxxxxz[k] * ab_y + g_x_0_xxyz_xxxxxyz[k];

                g_x_0_xxyyz_xxxxyy[k] = -g_x_0_xxyz_xxxxyy[k] * ab_y + g_x_0_xxyz_xxxxyyy[k];

                g_x_0_xxyyz_xxxxyz[k] = -g_x_0_xxyz_xxxxyz[k] * ab_y + g_x_0_xxyz_xxxxyyz[k];

                g_x_0_xxyyz_xxxxzz[k] = -g_x_0_xxyz_xxxxzz[k] * ab_y + g_x_0_xxyz_xxxxyzz[k];

                g_x_0_xxyyz_xxxyyy[k] = -g_x_0_xxyz_xxxyyy[k] * ab_y + g_x_0_xxyz_xxxyyyy[k];

                g_x_0_xxyyz_xxxyyz[k] = -g_x_0_xxyz_xxxyyz[k] * ab_y + g_x_0_xxyz_xxxyyyz[k];

                g_x_0_xxyyz_xxxyzz[k] = -g_x_0_xxyz_xxxyzz[k] * ab_y + g_x_0_xxyz_xxxyyzz[k];

                g_x_0_xxyyz_xxxzzz[k] = -g_x_0_xxyz_xxxzzz[k] * ab_y + g_x_0_xxyz_xxxyzzz[k];

                g_x_0_xxyyz_xxyyyy[k] = -g_x_0_xxyz_xxyyyy[k] * ab_y + g_x_0_xxyz_xxyyyyy[k];

                g_x_0_xxyyz_xxyyyz[k] = -g_x_0_xxyz_xxyyyz[k] * ab_y + g_x_0_xxyz_xxyyyyz[k];

                g_x_0_xxyyz_xxyyzz[k] = -g_x_0_xxyz_xxyyzz[k] * ab_y + g_x_0_xxyz_xxyyyzz[k];

                g_x_0_xxyyz_xxyzzz[k] = -g_x_0_xxyz_xxyzzz[k] * ab_y + g_x_0_xxyz_xxyyzzz[k];

                g_x_0_xxyyz_xxzzzz[k] = -g_x_0_xxyz_xxzzzz[k] * ab_y + g_x_0_xxyz_xxyzzzz[k];

                g_x_0_xxyyz_xyyyyy[k] = -g_x_0_xxyz_xyyyyy[k] * ab_y + g_x_0_xxyz_xyyyyyy[k];

                g_x_0_xxyyz_xyyyyz[k] = -g_x_0_xxyz_xyyyyz[k] * ab_y + g_x_0_xxyz_xyyyyyz[k];

                g_x_0_xxyyz_xyyyzz[k] = -g_x_0_xxyz_xyyyzz[k] * ab_y + g_x_0_xxyz_xyyyyzz[k];

                g_x_0_xxyyz_xyyzzz[k] = -g_x_0_xxyz_xyyzzz[k] * ab_y + g_x_0_xxyz_xyyyzzz[k];

                g_x_0_xxyyz_xyzzzz[k] = -g_x_0_xxyz_xyzzzz[k] * ab_y + g_x_0_xxyz_xyyzzzz[k];

                g_x_0_xxyyz_xzzzzz[k] = -g_x_0_xxyz_xzzzzz[k] * ab_y + g_x_0_xxyz_xyzzzzz[k];

                g_x_0_xxyyz_yyyyyy[k] = -g_x_0_xxyz_yyyyyy[k] * ab_y + g_x_0_xxyz_yyyyyyy[k];

                g_x_0_xxyyz_yyyyyz[k] = -g_x_0_xxyz_yyyyyz[k] * ab_y + g_x_0_xxyz_yyyyyyz[k];

                g_x_0_xxyyz_yyyyzz[k] = -g_x_0_xxyz_yyyyzz[k] * ab_y + g_x_0_xxyz_yyyyyzz[k];

                g_x_0_xxyyz_yyyzzz[k] = -g_x_0_xxyz_yyyzzz[k] * ab_y + g_x_0_xxyz_yyyyzzz[k];

                g_x_0_xxyyz_yyzzzz[k] = -g_x_0_xxyz_yyzzzz[k] * ab_y + g_x_0_xxyz_yyyzzzz[k];

                g_x_0_xxyyz_yzzzzz[k] = -g_x_0_xxyz_yzzzzz[k] * ab_y + g_x_0_xxyz_yyzzzzz[k];

                g_x_0_xxyyz_zzzzzz[k] = -g_x_0_xxyz_zzzzzz[k] * ab_y + g_x_0_xxyz_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xxyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xxyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyzz_xxxxxx, g_x_0_xxyzz_xxxxxy, g_x_0_xxyzz_xxxxxz, g_x_0_xxyzz_xxxxyy, g_x_0_xxyzz_xxxxyz, g_x_0_xxyzz_xxxxzz, g_x_0_xxyzz_xxxyyy, g_x_0_xxyzz_xxxyyz, g_x_0_xxyzz_xxxyzz, g_x_0_xxyzz_xxxzzz, g_x_0_xxyzz_xxyyyy, g_x_0_xxyzz_xxyyyz, g_x_0_xxyzz_xxyyzz, g_x_0_xxyzz_xxyzzz, g_x_0_xxyzz_xxzzzz, g_x_0_xxyzz_xyyyyy, g_x_0_xxyzz_xyyyyz, g_x_0_xxyzz_xyyyzz, g_x_0_xxyzz_xyyzzz, g_x_0_xxyzz_xyzzzz, g_x_0_xxyzz_xzzzzz, g_x_0_xxyzz_yyyyyy, g_x_0_xxyzz_yyyyyz, g_x_0_xxyzz_yyyyzz, g_x_0_xxyzz_yyyzzz, g_x_0_xxyzz_yyzzzz, g_x_0_xxyzz_yzzzzz, g_x_0_xxyzz_zzzzzz, g_x_0_xxzz_xxxxxx, g_x_0_xxzz_xxxxxxy, g_x_0_xxzz_xxxxxy, g_x_0_xxzz_xxxxxyy, g_x_0_xxzz_xxxxxyz, g_x_0_xxzz_xxxxxz, g_x_0_xxzz_xxxxyy, g_x_0_xxzz_xxxxyyy, g_x_0_xxzz_xxxxyyz, g_x_0_xxzz_xxxxyz, g_x_0_xxzz_xxxxyzz, g_x_0_xxzz_xxxxzz, g_x_0_xxzz_xxxyyy, g_x_0_xxzz_xxxyyyy, g_x_0_xxzz_xxxyyyz, g_x_0_xxzz_xxxyyz, g_x_0_xxzz_xxxyyzz, g_x_0_xxzz_xxxyzz, g_x_0_xxzz_xxxyzzz, g_x_0_xxzz_xxxzzz, g_x_0_xxzz_xxyyyy, g_x_0_xxzz_xxyyyyy, g_x_0_xxzz_xxyyyyz, g_x_0_xxzz_xxyyyz, g_x_0_xxzz_xxyyyzz, g_x_0_xxzz_xxyyzz, g_x_0_xxzz_xxyyzzz, g_x_0_xxzz_xxyzzz, g_x_0_xxzz_xxyzzzz, g_x_0_xxzz_xxzzzz, g_x_0_xxzz_xyyyyy, g_x_0_xxzz_xyyyyyy, g_x_0_xxzz_xyyyyyz, g_x_0_xxzz_xyyyyz, g_x_0_xxzz_xyyyyzz, g_x_0_xxzz_xyyyzz, g_x_0_xxzz_xyyyzzz, g_x_0_xxzz_xyyzzz, g_x_0_xxzz_xyyzzzz, g_x_0_xxzz_xyzzzz, g_x_0_xxzz_xyzzzzz, g_x_0_xxzz_xzzzzz, g_x_0_xxzz_yyyyyy, g_x_0_xxzz_yyyyyyy, g_x_0_xxzz_yyyyyyz, g_x_0_xxzz_yyyyyz, g_x_0_xxzz_yyyyyzz, g_x_0_xxzz_yyyyzz, g_x_0_xxzz_yyyyzzz, g_x_0_xxzz_yyyzzz, g_x_0_xxzz_yyyzzzz, g_x_0_xxzz_yyzzzz, g_x_0_xxzz_yyzzzzz, g_x_0_xxzz_yzzzzz, g_x_0_xxzz_yzzzzzz, g_x_0_xxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzz_xxxxxx[k] = -g_x_0_xxzz_xxxxxx[k] * ab_y + g_x_0_xxzz_xxxxxxy[k];

                g_x_0_xxyzz_xxxxxy[k] = -g_x_0_xxzz_xxxxxy[k] * ab_y + g_x_0_xxzz_xxxxxyy[k];

                g_x_0_xxyzz_xxxxxz[k] = -g_x_0_xxzz_xxxxxz[k] * ab_y + g_x_0_xxzz_xxxxxyz[k];

                g_x_0_xxyzz_xxxxyy[k] = -g_x_0_xxzz_xxxxyy[k] * ab_y + g_x_0_xxzz_xxxxyyy[k];

                g_x_0_xxyzz_xxxxyz[k] = -g_x_0_xxzz_xxxxyz[k] * ab_y + g_x_0_xxzz_xxxxyyz[k];

                g_x_0_xxyzz_xxxxzz[k] = -g_x_0_xxzz_xxxxzz[k] * ab_y + g_x_0_xxzz_xxxxyzz[k];

                g_x_0_xxyzz_xxxyyy[k] = -g_x_0_xxzz_xxxyyy[k] * ab_y + g_x_0_xxzz_xxxyyyy[k];

                g_x_0_xxyzz_xxxyyz[k] = -g_x_0_xxzz_xxxyyz[k] * ab_y + g_x_0_xxzz_xxxyyyz[k];

                g_x_0_xxyzz_xxxyzz[k] = -g_x_0_xxzz_xxxyzz[k] * ab_y + g_x_0_xxzz_xxxyyzz[k];

                g_x_0_xxyzz_xxxzzz[k] = -g_x_0_xxzz_xxxzzz[k] * ab_y + g_x_0_xxzz_xxxyzzz[k];

                g_x_0_xxyzz_xxyyyy[k] = -g_x_0_xxzz_xxyyyy[k] * ab_y + g_x_0_xxzz_xxyyyyy[k];

                g_x_0_xxyzz_xxyyyz[k] = -g_x_0_xxzz_xxyyyz[k] * ab_y + g_x_0_xxzz_xxyyyyz[k];

                g_x_0_xxyzz_xxyyzz[k] = -g_x_0_xxzz_xxyyzz[k] * ab_y + g_x_0_xxzz_xxyyyzz[k];

                g_x_0_xxyzz_xxyzzz[k] = -g_x_0_xxzz_xxyzzz[k] * ab_y + g_x_0_xxzz_xxyyzzz[k];

                g_x_0_xxyzz_xxzzzz[k] = -g_x_0_xxzz_xxzzzz[k] * ab_y + g_x_0_xxzz_xxyzzzz[k];

                g_x_0_xxyzz_xyyyyy[k] = -g_x_0_xxzz_xyyyyy[k] * ab_y + g_x_0_xxzz_xyyyyyy[k];

                g_x_0_xxyzz_xyyyyz[k] = -g_x_0_xxzz_xyyyyz[k] * ab_y + g_x_0_xxzz_xyyyyyz[k];

                g_x_0_xxyzz_xyyyzz[k] = -g_x_0_xxzz_xyyyzz[k] * ab_y + g_x_0_xxzz_xyyyyzz[k];

                g_x_0_xxyzz_xyyzzz[k] = -g_x_0_xxzz_xyyzzz[k] * ab_y + g_x_0_xxzz_xyyyzzz[k];

                g_x_0_xxyzz_xyzzzz[k] = -g_x_0_xxzz_xyzzzz[k] * ab_y + g_x_0_xxzz_xyyzzzz[k];

                g_x_0_xxyzz_xzzzzz[k] = -g_x_0_xxzz_xzzzzz[k] * ab_y + g_x_0_xxzz_xyzzzzz[k];

                g_x_0_xxyzz_yyyyyy[k] = -g_x_0_xxzz_yyyyyy[k] * ab_y + g_x_0_xxzz_yyyyyyy[k];

                g_x_0_xxyzz_yyyyyz[k] = -g_x_0_xxzz_yyyyyz[k] * ab_y + g_x_0_xxzz_yyyyyyz[k];

                g_x_0_xxyzz_yyyyzz[k] = -g_x_0_xxzz_yyyyzz[k] * ab_y + g_x_0_xxzz_yyyyyzz[k];

                g_x_0_xxyzz_yyyzzz[k] = -g_x_0_xxzz_yyyzzz[k] * ab_y + g_x_0_xxzz_yyyyzzz[k];

                g_x_0_xxyzz_yyzzzz[k] = -g_x_0_xxzz_yyzzzz[k] * ab_y + g_x_0_xxzz_yyyzzzz[k];

                g_x_0_xxyzz_yzzzzz[k] = -g_x_0_xxzz_yzzzzz[k] * ab_y + g_x_0_xxzz_yyzzzzz[k];

                g_x_0_xxyzz_zzzzzz[k] = -g_x_0_xxzz_zzzzzz[k] * ab_y + g_x_0_xxzz_yzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xxzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xxzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxzz_xxxxxx, g_x_0_xxzz_xxxxxxz, g_x_0_xxzz_xxxxxy, g_x_0_xxzz_xxxxxyz, g_x_0_xxzz_xxxxxz, g_x_0_xxzz_xxxxxzz, g_x_0_xxzz_xxxxyy, g_x_0_xxzz_xxxxyyz, g_x_0_xxzz_xxxxyz, g_x_0_xxzz_xxxxyzz, g_x_0_xxzz_xxxxzz, g_x_0_xxzz_xxxxzzz, g_x_0_xxzz_xxxyyy, g_x_0_xxzz_xxxyyyz, g_x_0_xxzz_xxxyyz, g_x_0_xxzz_xxxyyzz, g_x_0_xxzz_xxxyzz, g_x_0_xxzz_xxxyzzz, g_x_0_xxzz_xxxzzz, g_x_0_xxzz_xxxzzzz, g_x_0_xxzz_xxyyyy, g_x_0_xxzz_xxyyyyz, g_x_0_xxzz_xxyyyz, g_x_0_xxzz_xxyyyzz, g_x_0_xxzz_xxyyzz, g_x_0_xxzz_xxyyzzz, g_x_0_xxzz_xxyzzz, g_x_0_xxzz_xxyzzzz, g_x_0_xxzz_xxzzzz, g_x_0_xxzz_xxzzzzz, g_x_0_xxzz_xyyyyy, g_x_0_xxzz_xyyyyyz, g_x_0_xxzz_xyyyyz, g_x_0_xxzz_xyyyyzz, g_x_0_xxzz_xyyyzz, g_x_0_xxzz_xyyyzzz, g_x_0_xxzz_xyyzzz, g_x_0_xxzz_xyyzzzz, g_x_0_xxzz_xyzzzz, g_x_0_xxzz_xyzzzzz, g_x_0_xxzz_xzzzzz, g_x_0_xxzz_xzzzzzz, g_x_0_xxzz_yyyyyy, g_x_0_xxzz_yyyyyyz, g_x_0_xxzz_yyyyyz, g_x_0_xxzz_yyyyyzz, g_x_0_xxzz_yyyyzz, g_x_0_xxzz_yyyyzzz, g_x_0_xxzz_yyyzzz, g_x_0_xxzz_yyyzzzz, g_x_0_xxzz_yyzzzz, g_x_0_xxzz_yyzzzzz, g_x_0_xxzz_yzzzzz, g_x_0_xxzz_yzzzzzz, g_x_0_xxzz_zzzzzz, g_x_0_xxzz_zzzzzzz, g_x_0_xxzzz_xxxxxx, g_x_0_xxzzz_xxxxxy, g_x_0_xxzzz_xxxxxz, g_x_0_xxzzz_xxxxyy, g_x_0_xxzzz_xxxxyz, g_x_0_xxzzz_xxxxzz, g_x_0_xxzzz_xxxyyy, g_x_0_xxzzz_xxxyyz, g_x_0_xxzzz_xxxyzz, g_x_0_xxzzz_xxxzzz, g_x_0_xxzzz_xxyyyy, g_x_0_xxzzz_xxyyyz, g_x_0_xxzzz_xxyyzz, g_x_0_xxzzz_xxyzzz, g_x_0_xxzzz_xxzzzz, g_x_0_xxzzz_xyyyyy, g_x_0_xxzzz_xyyyyz, g_x_0_xxzzz_xyyyzz, g_x_0_xxzzz_xyyzzz, g_x_0_xxzzz_xyzzzz, g_x_0_xxzzz_xzzzzz, g_x_0_xxzzz_yyyyyy, g_x_0_xxzzz_yyyyyz, g_x_0_xxzzz_yyyyzz, g_x_0_xxzzz_yyyzzz, g_x_0_xxzzz_yyzzzz, g_x_0_xxzzz_yzzzzz, g_x_0_xxzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzz_xxxxxx[k] = -g_x_0_xxzz_xxxxxx[k] * ab_z + g_x_0_xxzz_xxxxxxz[k];

                g_x_0_xxzzz_xxxxxy[k] = -g_x_0_xxzz_xxxxxy[k] * ab_z + g_x_0_xxzz_xxxxxyz[k];

                g_x_0_xxzzz_xxxxxz[k] = -g_x_0_xxzz_xxxxxz[k] * ab_z + g_x_0_xxzz_xxxxxzz[k];

                g_x_0_xxzzz_xxxxyy[k] = -g_x_0_xxzz_xxxxyy[k] * ab_z + g_x_0_xxzz_xxxxyyz[k];

                g_x_0_xxzzz_xxxxyz[k] = -g_x_0_xxzz_xxxxyz[k] * ab_z + g_x_0_xxzz_xxxxyzz[k];

                g_x_0_xxzzz_xxxxzz[k] = -g_x_0_xxzz_xxxxzz[k] * ab_z + g_x_0_xxzz_xxxxzzz[k];

                g_x_0_xxzzz_xxxyyy[k] = -g_x_0_xxzz_xxxyyy[k] * ab_z + g_x_0_xxzz_xxxyyyz[k];

                g_x_0_xxzzz_xxxyyz[k] = -g_x_0_xxzz_xxxyyz[k] * ab_z + g_x_0_xxzz_xxxyyzz[k];

                g_x_0_xxzzz_xxxyzz[k] = -g_x_0_xxzz_xxxyzz[k] * ab_z + g_x_0_xxzz_xxxyzzz[k];

                g_x_0_xxzzz_xxxzzz[k] = -g_x_0_xxzz_xxxzzz[k] * ab_z + g_x_0_xxzz_xxxzzzz[k];

                g_x_0_xxzzz_xxyyyy[k] = -g_x_0_xxzz_xxyyyy[k] * ab_z + g_x_0_xxzz_xxyyyyz[k];

                g_x_0_xxzzz_xxyyyz[k] = -g_x_0_xxzz_xxyyyz[k] * ab_z + g_x_0_xxzz_xxyyyzz[k];

                g_x_0_xxzzz_xxyyzz[k] = -g_x_0_xxzz_xxyyzz[k] * ab_z + g_x_0_xxzz_xxyyzzz[k];

                g_x_0_xxzzz_xxyzzz[k] = -g_x_0_xxzz_xxyzzz[k] * ab_z + g_x_0_xxzz_xxyzzzz[k];

                g_x_0_xxzzz_xxzzzz[k] = -g_x_0_xxzz_xxzzzz[k] * ab_z + g_x_0_xxzz_xxzzzzz[k];

                g_x_0_xxzzz_xyyyyy[k] = -g_x_0_xxzz_xyyyyy[k] * ab_z + g_x_0_xxzz_xyyyyyz[k];

                g_x_0_xxzzz_xyyyyz[k] = -g_x_0_xxzz_xyyyyz[k] * ab_z + g_x_0_xxzz_xyyyyzz[k];

                g_x_0_xxzzz_xyyyzz[k] = -g_x_0_xxzz_xyyyzz[k] * ab_z + g_x_0_xxzz_xyyyzzz[k];

                g_x_0_xxzzz_xyyzzz[k] = -g_x_0_xxzz_xyyzzz[k] * ab_z + g_x_0_xxzz_xyyzzzz[k];

                g_x_0_xxzzz_xyzzzz[k] = -g_x_0_xxzz_xyzzzz[k] * ab_z + g_x_0_xxzz_xyzzzzz[k];

                g_x_0_xxzzz_xzzzzz[k] = -g_x_0_xxzz_xzzzzz[k] * ab_z + g_x_0_xxzz_xzzzzzz[k];

                g_x_0_xxzzz_yyyyyy[k] = -g_x_0_xxzz_yyyyyy[k] * ab_z + g_x_0_xxzz_yyyyyyz[k];

                g_x_0_xxzzz_yyyyyz[k] = -g_x_0_xxzz_yyyyyz[k] * ab_z + g_x_0_xxzz_yyyyyzz[k];

                g_x_0_xxzzz_yyyyzz[k] = -g_x_0_xxzz_yyyyzz[k] * ab_z + g_x_0_xxzz_yyyyzzz[k];

                g_x_0_xxzzz_yyyzzz[k] = -g_x_0_xxzz_yyyzzz[k] * ab_z + g_x_0_xxzz_yyyzzzz[k];

                g_x_0_xxzzz_yyzzzz[k] = -g_x_0_xxzz_yyzzzz[k] * ab_z + g_x_0_xxzz_yyzzzzz[k];

                g_x_0_xxzzz_yzzzzz[k] = -g_x_0_xxzz_yzzzzz[k] * ab_z + g_x_0_xxzz_yzzzzzz[k];

                g_x_0_xxzzz_zzzzzz[k] = -g_x_0_xxzz_zzzzzz[k] * ab_z + g_x_0_xxzz_zzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_xyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_xyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 307 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyy_xxxxxx, g_x_0_xyyy_xxxxxxy, g_x_0_xyyy_xxxxxy, g_x_0_xyyy_xxxxxyy, g_x_0_xyyy_xxxxxyz, g_x_0_xyyy_xxxxxz, g_x_0_xyyy_xxxxyy, g_x_0_xyyy_xxxxyyy, g_x_0_xyyy_xxxxyyz, g_x_0_xyyy_xxxxyz, g_x_0_xyyy_xxxxyzz, g_x_0_xyyy_xxxxzz, g_x_0_xyyy_xxxyyy, g_x_0_xyyy_xxxyyyy, g_x_0_xyyy_xxxyyyz, g_x_0_xyyy_xxxyyz, g_x_0_xyyy_xxxyyzz, g_x_0_xyyy_xxxyzz, g_x_0_xyyy_xxxyzzz, g_x_0_xyyy_xxxzzz, g_x_0_xyyy_xxyyyy, g_x_0_xyyy_xxyyyyy, g_x_0_xyyy_xxyyyyz, g_x_0_xyyy_xxyyyz, g_x_0_xyyy_xxyyyzz, g_x_0_xyyy_xxyyzz, g_x_0_xyyy_xxyyzzz, g_x_0_xyyy_xxyzzz, g_x_0_xyyy_xxyzzzz, g_x_0_xyyy_xxzzzz, g_x_0_xyyy_xyyyyy, g_x_0_xyyy_xyyyyyy, g_x_0_xyyy_xyyyyyz, g_x_0_xyyy_xyyyyz, g_x_0_xyyy_xyyyyzz, g_x_0_xyyy_xyyyzz, g_x_0_xyyy_xyyyzzz, g_x_0_xyyy_xyyzzz, g_x_0_xyyy_xyyzzzz, g_x_0_xyyy_xyzzzz, g_x_0_xyyy_xyzzzzz, g_x_0_xyyy_xzzzzz, g_x_0_xyyy_yyyyyy, g_x_0_xyyy_yyyyyyy, g_x_0_xyyy_yyyyyyz, g_x_0_xyyy_yyyyyz, g_x_0_xyyy_yyyyyzz, g_x_0_xyyy_yyyyzz, g_x_0_xyyy_yyyyzzz, g_x_0_xyyy_yyyzzz, g_x_0_xyyy_yyyzzzz, g_x_0_xyyy_yyzzzz, g_x_0_xyyy_yyzzzzz, g_x_0_xyyy_yzzzzz, g_x_0_xyyy_yzzzzzz, g_x_0_xyyy_zzzzzz, g_x_0_xyyyy_xxxxxx, g_x_0_xyyyy_xxxxxy, g_x_0_xyyyy_xxxxxz, g_x_0_xyyyy_xxxxyy, g_x_0_xyyyy_xxxxyz, g_x_0_xyyyy_xxxxzz, g_x_0_xyyyy_xxxyyy, g_x_0_xyyyy_xxxyyz, g_x_0_xyyyy_xxxyzz, g_x_0_xyyyy_xxxzzz, g_x_0_xyyyy_xxyyyy, g_x_0_xyyyy_xxyyyz, g_x_0_xyyyy_xxyyzz, g_x_0_xyyyy_xxyzzz, g_x_0_xyyyy_xxzzzz, g_x_0_xyyyy_xyyyyy, g_x_0_xyyyy_xyyyyz, g_x_0_xyyyy_xyyyzz, g_x_0_xyyyy_xyyzzz, g_x_0_xyyyy_xyzzzz, g_x_0_xyyyy_xzzzzz, g_x_0_xyyyy_yyyyyy, g_x_0_xyyyy_yyyyyz, g_x_0_xyyyy_yyyyzz, g_x_0_xyyyy_yyyzzz, g_x_0_xyyyy_yyzzzz, g_x_0_xyyyy_yzzzzz, g_x_0_xyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyy_xxxxxx[k] = -g_x_0_xyyy_xxxxxx[k] * ab_y + g_x_0_xyyy_xxxxxxy[k];

                g_x_0_xyyyy_xxxxxy[k] = -g_x_0_xyyy_xxxxxy[k] * ab_y + g_x_0_xyyy_xxxxxyy[k];

                g_x_0_xyyyy_xxxxxz[k] = -g_x_0_xyyy_xxxxxz[k] * ab_y + g_x_0_xyyy_xxxxxyz[k];

                g_x_0_xyyyy_xxxxyy[k] = -g_x_0_xyyy_xxxxyy[k] * ab_y + g_x_0_xyyy_xxxxyyy[k];

                g_x_0_xyyyy_xxxxyz[k] = -g_x_0_xyyy_xxxxyz[k] * ab_y + g_x_0_xyyy_xxxxyyz[k];

                g_x_0_xyyyy_xxxxzz[k] = -g_x_0_xyyy_xxxxzz[k] * ab_y + g_x_0_xyyy_xxxxyzz[k];

                g_x_0_xyyyy_xxxyyy[k] = -g_x_0_xyyy_xxxyyy[k] * ab_y + g_x_0_xyyy_xxxyyyy[k];

                g_x_0_xyyyy_xxxyyz[k] = -g_x_0_xyyy_xxxyyz[k] * ab_y + g_x_0_xyyy_xxxyyyz[k];

                g_x_0_xyyyy_xxxyzz[k] = -g_x_0_xyyy_xxxyzz[k] * ab_y + g_x_0_xyyy_xxxyyzz[k];

                g_x_0_xyyyy_xxxzzz[k] = -g_x_0_xyyy_xxxzzz[k] * ab_y + g_x_0_xyyy_xxxyzzz[k];

                g_x_0_xyyyy_xxyyyy[k] = -g_x_0_xyyy_xxyyyy[k] * ab_y + g_x_0_xyyy_xxyyyyy[k];

                g_x_0_xyyyy_xxyyyz[k] = -g_x_0_xyyy_xxyyyz[k] * ab_y + g_x_0_xyyy_xxyyyyz[k];

                g_x_0_xyyyy_xxyyzz[k] = -g_x_0_xyyy_xxyyzz[k] * ab_y + g_x_0_xyyy_xxyyyzz[k];

                g_x_0_xyyyy_xxyzzz[k] = -g_x_0_xyyy_xxyzzz[k] * ab_y + g_x_0_xyyy_xxyyzzz[k];

                g_x_0_xyyyy_xxzzzz[k] = -g_x_0_xyyy_xxzzzz[k] * ab_y + g_x_0_xyyy_xxyzzzz[k];

                g_x_0_xyyyy_xyyyyy[k] = -g_x_0_xyyy_xyyyyy[k] * ab_y + g_x_0_xyyy_xyyyyyy[k];

                g_x_0_xyyyy_xyyyyz[k] = -g_x_0_xyyy_xyyyyz[k] * ab_y + g_x_0_xyyy_xyyyyyz[k];

                g_x_0_xyyyy_xyyyzz[k] = -g_x_0_xyyy_xyyyzz[k] * ab_y + g_x_0_xyyy_xyyyyzz[k];

                g_x_0_xyyyy_xyyzzz[k] = -g_x_0_xyyy_xyyzzz[k] * ab_y + g_x_0_xyyy_xyyyzzz[k];

                g_x_0_xyyyy_xyzzzz[k] = -g_x_0_xyyy_xyzzzz[k] * ab_y + g_x_0_xyyy_xyyzzzz[k];

                g_x_0_xyyyy_xzzzzz[k] = -g_x_0_xyyy_xzzzzz[k] * ab_y + g_x_0_xyyy_xyzzzzz[k];

                g_x_0_xyyyy_yyyyyy[k] = -g_x_0_xyyy_yyyyyy[k] * ab_y + g_x_0_xyyy_yyyyyyy[k];

                g_x_0_xyyyy_yyyyyz[k] = -g_x_0_xyyy_yyyyyz[k] * ab_y + g_x_0_xyyy_yyyyyyz[k];

                g_x_0_xyyyy_yyyyzz[k] = -g_x_0_xyyy_yyyyzz[k] * ab_y + g_x_0_xyyy_yyyyyzz[k];

                g_x_0_xyyyy_yyyzzz[k] = -g_x_0_xyyy_yyyzzz[k] * ab_y + g_x_0_xyyy_yyyyzzz[k];

                g_x_0_xyyyy_yyzzzz[k] = -g_x_0_xyyy_yyzzzz[k] * ab_y + g_x_0_xyyy_yyyzzzz[k];

                g_x_0_xyyyy_yzzzzz[k] = -g_x_0_xyyy_yzzzzz[k] * ab_y + g_x_0_xyyy_yyzzzzz[k];

                g_x_0_xyyyy_zzzzzz[k] = -g_x_0_xyyy_zzzzzz[k] * ab_y + g_x_0_xyyy_yzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_xyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_xyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyz_xxxxxx, g_x_0_xyyyz_xxxxxy, g_x_0_xyyyz_xxxxxz, g_x_0_xyyyz_xxxxyy, g_x_0_xyyyz_xxxxyz, g_x_0_xyyyz_xxxxzz, g_x_0_xyyyz_xxxyyy, g_x_0_xyyyz_xxxyyz, g_x_0_xyyyz_xxxyzz, g_x_0_xyyyz_xxxzzz, g_x_0_xyyyz_xxyyyy, g_x_0_xyyyz_xxyyyz, g_x_0_xyyyz_xxyyzz, g_x_0_xyyyz_xxyzzz, g_x_0_xyyyz_xxzzzz, g_x_0_xyyyz_xyyyyy, g_x_0_xyyyz_xyyyyz, g_x_0_xyyyz_xyyyzz, g_x_0_xyyyz_xyyzzz, g_x_0_xyyyz_xyzzzz, g_x_0_xyyyz_xzzzzz, g_x_0_xyyyz_yyyyyy, g_x_0_xyyyz_yyyyyz, g_x_0_xyyyz_yyyyzz, g_x_0_xyyyz_yyyzzz, g_x_0_xyyyz_yyzzzz, g_x_0_xyyyz_yzzzzz, g_x_0_xyyyz_zzzzzz, g_x_0_xyyz_xxxxxx, g_x_0_xyyz_xxxxxxy, g_x_0_xyyz_xxxxxy, g_x_0_xyyz_xxxxxyy, g_x_0_xyyz_xxxxxyz, g_x_0_xyyz_xxxxxz, g_x_0_xyyz_xxxxyy, g_x_0_xyyz_xxxxyyy, g_x_0_xyyz_xxxxyyz, g_x_0_xyyz_xxxxyz, g_x_0_xyyz_xxxxyzz, g_x_0_xyyz_xxxxzz, g_x_0_xyyz_xxxyyy, g_x_0_xyyz_xxxyyyy, g_x_0_xyyz_xxxyyyz, g_x_0_xyyz_xxxyyz, g_x_0_xyyz_xxxyyzz, g_x_0_xyyz_xxxyzz, g_x_0_xyyz_xxxyzzz, g_x_0_xyyz_xxxzzz, g_x_0_xyyz_xxyyyy, g_x_0_xyyz_xxyyyyy, g_x_0_xyyz_xxyyyyz, g_x_0_xyyz_xxyyyz, g_x_0_xyyz_xxyyyzz, g_x_0_xyyz_xxyyzz, g_x_0_xyyz_xxyyzzz, g_x_0_xyyz_xxyzzz, g_x_0_xyyz_xxyzzzz, g_x_0_xyyz_xxzzzz, g_x_0_xyyz_xyyyyy, g_x_0_xyyz_xyyyyyy, g_x_0_xyyz_xyyyyyz, g_x_0_xyyz_xyyyyz, g_x_0_xyyz_xyyyyzz, g_x_0_xyyz_xyyyzz, g_x_0_xyyz_xyyyzzz, g_x_0_xyyz_xyyzzz, g_x_0_xyyz_xyyzzzz, g_x_0_xyyz_xyzzzz, g_x_0_xyyz_xyzzzzz, g_x_0_xyyz_xzzzzz, g_x_0_xyyz_yyyyyy, g_x_0_xyyz_yyyyyyy, g_x_0_xyyz_yyyyyyz, g_x_0_xyyz_yyyyyz, g_x_0_xyyz_yyyyyzz, g_x_0_xyyz_yyyyzz, g_x_0_xyyz_yyyyzzz, g_x_0_xyyz_yyyzzz, g_x_0_xyyz_yyyzzzz, g_x_0_xyyz_yyzzzz, g_x_0_xyyz_yyzzzzz, g_x_0_xyyz_yzzzzz, g_x_0_xyyz_yzzzzzz, g_x_0_xyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyz_xxxxxx[k] = -g_x_0_xyyz_xxxxxx[k] * ab_y + g_x_0_xyyz_xxxxxxy[k];

                g_x_0_xyyyz_xxxxxy[k] = -g_x_0_xyyz_xxxxxy[k] * ab_y + g_x_0_xyyz_xxxxxyy[k];

                g_x_0_xyyyz_xxxxxz[k] = -g_x_0_xyyz_xxxxxz[k] * ab_y + g_x_0_xyyz_xxxxxyz[k];

                g_x_0_xyyyz_xxxxyy[k] = -g_x_0_xyyz_xxxxyy[k] * ab_y + g_x_0_xyyz_xxxxyyy[k];

                g_x_0_xyyyz_xxxxyz[k] = -g_x_0_xyyz_xxxxyz[k] * ab_y + g_x_0_xyyz_xxxxyyz[k];

                g_x_0_xyyyz_xxxxzz[k] = -g_x_0_xyyz_xxxxzz[k] * ab_y + g_x_0_xyyz_xxxxyzz[k];

                g_x_0_xyyyz_xxxyyy[k] = -g_x_0_xyyz_xxxyyy[k] * ab_y + g_x_0_xyyz_xxxyyyy[k];

                g_x_0_xyyyz_xxxyyz[k] = -g_x_0_xyyz_xxxyyz[k] * ab_y + g_x_0_xyyz_xxxyyyz[k];

                g_x_0_xyyyz_xxxyzz[k] = -g_x_0_xyyz_xxxyzz[k] * ab_y + g_x_0_xyyz_xxxyyzz[k];

                g_x_0_xyyyz_xxxzzz[k] = -g_x_0_xyyz_xxxzzz[k] * ab_y + g_x_0_xyyz_xxxyzzz[k];

                g_x_0_xyyyz_xxyyyy[k] = -g_x_0_xyyz_xxyyyy[k] * ab_y + g_x_0_xyyz_xxyyyyy[k];

                g_x_0_xyyyz_xxyyyz[k] = -g_x_0_xyyz_xxyyyz[k] * ab_y + g_x_0_xyyz_xxyyyyz[k];

                g_x_0_xyyyz_xxyyzz[k] = -g_x_0_xyyz_xxyyzz[k] * ab_y + g_x_0_xyyz_xxyyyzz[k];

                g_x_0_xyyyz_xxyzzz[k] = -g_x_0_xyyz_xxyzzz[k] * ab_y + g_x_0_xyyz_xxyyzzz[k];

                g_x_0_xyyyz_xxzzzz[k] = -g_x_0_xyyz_xxzzzz[k] * ab_y + g_x_0_xyyz_xxyzzzz[k];

                g_x_0_xyyyz_xyyyyy[k] = -g_x_0_xyyz_xyyyyy[k] * ab_y + g_x_0_xyyz_xyyyyyy[k];

                g_x_0_xyyyz_xyyyyz[k] = -g_x_0_xyyz_xyyyyz[k] * ab_y + g_x_0_xyyz_xyyyyyz[k];

                g_x_0_xyyyz_xyyyzz[k] = -g_x_0_xyyz_xyyyzz[k] * ab_y + g_x_0_xyyz_xyyyyzz[k];

                g_x_0_xyyyz_xyyzzz[k] = -g_x_0_xyyz_xyyzzz[k] * ab_y + g_x_0_xyyz_xyyyzzz[k];

                g_x_0_xyyyz_xyzzzz[k] = -g_x_0_xyyz_xyzzzz[k] * ab_y + g_x_0_xyyz_xyyzzzz[k];

                g_x_0_xyyyz_xzzzzz[k] = -g_x_0_xyyz_xzzzzz[k] * ab_y + g_x_0_xyyz_xyzzzzz[k];

                g_x_0_xyyyz_yyyyyy[k] = -g_x_0_xyyz_yyyyyy[k] * ab_y + g_x_0_xyyz_yyyyyyy[k];

                g_x_0_xyyyz_yyyyyz[k] = -g_x_0_xyyz_yyyyyz[k] * ab_y + g_x_0_xyyz_yyyyyyz[k];

                g_x_0_xyyyz_yyyyzz[k] = -g_x_0_xyyz_yyyyzz[k] * ab_y + g_x_0_xyyz_yyyyyzz[k];

                g_x_0_xyyyz_yyyzzz[k] = -g_x_0_xyyz_yyyzzz[k] * ab_y + g_x_0_xyyz_yyyyzzz[k];

                g_x_0_xyyyz_yyzzzz[k] = -g_x_0_xyyz_yyzzzz[k] * ab_y + g_x_0_xyyz_yyyzzzz[k];

                g_x_0_xyyyz_yzzzzz[k] = -g_x_0_xyyz_yzzzzz[k] * ab_y + g_x_0_xyyz_yyzzzzz[k];

                g_x_0_xyyyz_zzzzzz[k] = -g_x_0_xyyz_zzzzzz[k] * ab_y + g_x_0_xyyz_yzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 356 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 360 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_xyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 362 * ccomps * dcomps);

            auto g_x_0_xyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 363 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyzz_xxxxxx, g_x_0_xyyzz_xxxxxy, g_x_0_xyyzz_xxxxxz, g_x_0_xyyzz_xxxxyy, g_x_0_xyyzz_xxxxyz, g_x_0_xyyzz_xxxxzz, g_x_0_xyyzz_xxxyyy, g_x_0_xyyzz_xxxyyz, g_x_0_xyyzz_xxxyzz, g_x_0_xyyzz_xxxzzz, g_x_0_xyyzz_xxyyyy, g_x_0_xyyzz_xxyyyz, g_x_0_xyyzz_xxyyzz, g_x_0_xyyzz_xxyzzz, g_x_0_xyyzz_xxzzzz, g_x_0_xyyzz_xyyyyy, g_x_0_xyyzz_xyyyyz, g_x_0_xyyzz_xyyyzz, g_x_0_xyyzz_xyyzzz, g_x_0_xyyzz_xyzzzz, g_x_0_xyyzz_xzzzzz, g_x_0_xyyzz_yyyyyy, g_x_0_xyyzz_yyyyyz, g_x_0_xyyzz_yyyyzz, g_x_0_xyyzz_yyyzzz, g_x_0_xyyzz_yyzzzz, g_x_0_xyyzz_yzzzzz, g_x_0_xyyzz_zzzzzz, g_x_0_xyzz_xxxxxx, g_x_0_xyzz_xxxxxxy, g_x_0_xyzz_xxxxxy, g_x_0_xyzz_xxxxxyy, g_x_0_xyzz_xxxxxyz, g_x_0_xyzz_xxxxxz, g_x_0_xyzz_xxxxyy, g_x_0_xyzz_xxxxyyy, g_x_0_xyzz_xxxxyyz, g_x_0_xyzz_xxxxyz, g_x_0_xyzz_xxxxyzz, g_x_0_xyzz_xxxxzz, g_x_0_xyzz_xxxyyy, g_x_0_xyzz_xxxyyyy, g_x_0_xyzz_xxxyyyz, g_x_0_xyzz_xxxyyz, g_x_0_xyzz_xxxyyzz, g_x_0_xyzz_xxxyzz, g_x_0_xyzz_xxxyzzz, g_x_0_xyzz_xxxzzz, g_x_0_xyzz_xxyyyy, g_x_0_xyzz_xxyyyyy, g_x_0_xyzz_xxyyyyz, g_x_0_xyzz_xxyyyz, g_x_0_xyzz_xxyyyzz, g_x_0_xyzz_xxyyzz, g_x_0_xyzz_xxyyzzz, g_x_0_xyzz_xxyzzz, g_x_0_xyzz_xxyzzzz, g_x_0_xyzz_xxzzzz, g_x_0_xyzz_xyyyyy, g_x_0_xyzz_xyyyyyy, g_x_0_xyzz_xyyyyyz, g_x_0_xyzz_xyyyyz, g_x_0_xyzz_xyyyyzz, g_x_0_xyzz_xyyyzz, g_x_0_xyzz_xyyyzzz, g_x_0_xyzz_xyyzzz, g_x_0_xyzz_xyyzzzz, g_x_0_xyzz_xyzzzz, g_x_0_xyzz_xyzzzzz, g_x_0_xyzz_xzzzzz, g_x_0_xyzz_yyyyyy, g_x_0_xyzz_yyyyyyy, g_x_0_xyzz_yyyyyyz, g_x_0_xyzz_yyyyyz, g_x_0_xyzz_yyyyyzz, g_x_0_xyzz_yyyyzz, g_x_0_xyzz_yyyyzzz, g_x_0_xyzz_yyyzzz, g_x_0_xyzz_yyyzzzz, g_x_0_xyzz_yyzzzz, g_x_0_xyzz_yyzzzzz, g_x_0_xyzz_yzzzzz, g_x_0_xyzz_yzzzzzz, g_x_0_xyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzz_xxxxxx[k] = -g_x_0_xyzz_xxxxxx[k] * ab_y + g_x_0_xyzz_xxxxxxy[k];

                g_x_0_xyyzz_xxxxxy[k] = -g_x_0_xyzz_xxxxxy[k] * ab_y + g_x_0_xyzz_xxxxxyy[k];

                g_x_0_xyyzz_xxxxxz[k] = -g_x_0_xyzz_xxxxxz[k] * ab_y + g_x_0_xyzz_xxxxxyz[k];

                g_x_0_xyyzz_xxxxyy[k] = -g_x_0_xyzz_xxxxyy[k] * ab_y + g_x_0_xyzz_xxxxyyy[k];

                g_x_0_xyyzz_xxxxyz[k] = -g_x_0_xyzz_xxxxyz[k] * ab_y + g_x_0_xyzz_xxxxyyz[k];

                g_x_0_xyyzz_xxxxzz[k] = -g_x_0_xyzz_xxxxzz[k] * ab_y + g_x_0_xyzz_xxxxyzz[k];

                g_x_0_xyyzz_xxxyyy[k] = -g_x_0_xyzz_xxxyyy[k] * ab_y + g_x_0_xyzz_xxxyyyy[k];

                g_x_0_xyyzz_xxxyyz[k] = -g_x_0_xyzz_xxxyyz[k] * ab_y + g_x_0_xyzz_xxxyyyz[k];

                g_x_0_xyyzz_xxxyzz[k] = -g_x_0_xyzz_xxxyzz[k] * ab_y + g_x_0_xyzz_xxxyyzz[k];

                g_x_0_xyyzz_xxxzzz[k] = -g_x_0_xyzz_xxxzzz[k] * ab_y + g_x_0_xyzz_xxxyzzz[k];

                g_x_0_xyyzz_xxyyyy[k] = -g_x_0_xyzz_xxyyyy[k] * ab_y + g_x_0_xyzz_xxyyyyy[k];

                g_x_0_xyyzz_xxyyyz[k] = -g_x_0_xyzz_xxyyyz[k] * ab_y + g_x_0_xyzz_xxyyyyz[k];

                g_x_0_xyyzz_xxyyzz[k] = -g_x_0_xyzz_xxyyzz[k] * ab_y + g_x_0_xyzz_xxyyyzz[k];

                g_x_0_xyyzz_xxyzzz[k] = -g_x_0_xyzz_xxyzzz[k] * ab_y + g_x_0_xyzz_xxyyzzz[k];

                g_x_0_xyyzz_xxzzzz[k] = -g_x_0_xyzz_xxzzzz[k] * ab_y + g_x_0_xyzz_xxyzzzz[k];

                g_x_0_xyyzz_xyyyyy[k] = -g_x_0_xyzz_xyyyyy[k] * ab_y + g_x_0_xyzz_xyyyyyy[k];

                g_x_0_xyyzz_xyyyyz[k] = -g_x_0_xyzz_xyyyyz[k] * ab_y + g_x_0_xyzz_xyyyyyz[k];

                g_x_0_xyyzz_xyyyzz[k] = -g_x_0_xyzz_xyyyzz[k] * ab_y + g_x_0_xyzz_xyyyyzz[k];

                g_x_0_xyyzz_xyyzzz[k] = -g_x_0_xyzz_xyyzzz[k] * ab_y + g_x_0_xyzz_xyyyzzz[k];

                g_x_0_xyyzz_xyzzzz[k] = -g_x_0_xyzz_xyzzzz[k] * ab_y + g_x_0_xyzz_xyyzzzz[k];

                g_x_0_xyyzz_xzzzzz[k] = -g_x_0_xyzz_xzzzzz[k] * ab_y + g_x_0_xyzz_xyzzzzz[k];

                g_x_0_xyyzz_yyyyyy[k] = -g_x_0_xyzz_yyyyyy[k] * ab_y + g_x_0_xyzz_yyyyyyy[k];

                g_x_0_xyyzz_yyyyyz[k] = -g_x_0_xyzz_yyyyyz[k] * ab_y + g_x_0_xyzz_yyyyyyz[k];

                g_x_0_xyyzz_yyyyzz[k] = -g_x_0_xyzz_yyyyzz[k] * ab_y + g_x_0_xyzz_yyyyyzz[k];

                g_x_0_xyyzz_yyyzzz[k] = -g_x_0_xyzz_yyyzzz[k] * ab_y + g_x_0_xyzz_yyyyzzz[k];

                g_x_0_xyyzz_yyzzzz[k] = -g_x_0_xyzz_yyzzzz[k] * ab_y + g_x_0_xyzz_yyyzzzz[k];

                g_x_0_xyyzz_yzzzzz[k] = -g_x_0_xyzz_yzzzzz[k] * ab_y + g_x_0_xyzz_yyzzzzz[k];

                g_x_0_xyyzz_zzzzzz[k] = -g_x_0_xyzz_zzzzzz[k] * ab_y + g_x_0_xyzz_yzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 364 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 365 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 366 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 369 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 373 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 374 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 378 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 380 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 384 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 387 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_xyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_xyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 391 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzzz_xxxxxx, g_x_0_xyzzz_xxxxxy, g_x_0_xyzzz_xxxxxz, g_x_0_xyzzz_xxxxyy, g_x_0_xyzzz_xxxxyz, g_x_0_xyzzz_xxxxzz, g_x_0_xyzzz_xxxyyy, g_x_0_xyzzz_xxxyyz, g_x_0_xyzzz_xxxyzz, g_x_0_xyzzz_xxxzzz, g_x_0_xyzzz_xxyyyy, g_x_0_xyzzz_xxyyyz, g_x_0_xyzzz_xxyyzz, g_x_0_xyzzz_xxyzzz, g_x_0_xyzzz_xxzzzz, g_x_0_xyzzz_xyyyyy, g_x_0_xyzzz_xyyyyz, g_x_0_xyzzz_xyyyzz, g_x_0_xyzzz_xyyzzz, g_x_0_xyzzz_xyzzzz, g_x_0_xyzzz_xzzzzz, g_x_0_xyzzz_yyyyyy, g_x_0_xyzzz_yyyyyz, g_x_0_xyzzz_yyyyzz, g_x_0_xyzzz_yyyzzz, g_x_0_xyzzz_yyzzzz, g_x_0_xyzzz_yzzzzz, g_x_0_xyzzz_zzzzzz, g_x_0_xzzz_xxxxxx, g_x_0_xzzz_xxxxxxy, g_x_0_xzzz_xxxxxy, g_x_0_xzzz_xxxxxyy, g_x_0_xzzz_xxxxxyz, g_x_0_xzzz_xxxxxz, g_x_0_xzzz_xxxxyy, g_x_0_xzzz_xxxxyyy, g_x_0_xzzz_xxxxyyz, g_x_0_xzzz_xxxxyz, g_x_0_xzzz_xxxxyzz, g_x_0_xzzz_xxxxzz, g_x_0_xzzz_xxxyyy, g_x_0_xzzz_xxxyyyy, g_x_0_xzzz_xxxyyyz, g_x_0_xzzz_xxxyyz, g_x_0_xzzz_xxxyyzz, g_x_0_xzzz_xxxyzz, g_x_0_xzzz_xxxyzzz, g_x_0_xzzz_xxxzzz, g_x_0_xzzz_xxyyyy, g_x_0_xzzz_xxyyyyy, g_x_0_xzzz_xxyyyyz, g_x_0_xzzz_xxyyyz, g_x_0_xzzz_xxyyyzz, g_x_0_xzzz_xxyyzz, g_x_0_xzzz_xxyyzzz, g_x_0_xzzz_xxyzzz, g_x_0_xzzz_xxyzzzz, g_x_0_xzzz_xxzzzz, g_x_0_xzzz_xyyyyy, g_x_0_xzzz_xyyyyyy, g_x_0_xzzz_xyyyyyz, g_x_0_xzzz_xyyyyz, g_x_0_xzzz_xyyyyzz, g_x_0_xzzz_xyyyzz, g_x_0_xzzz_xyyyzzz, g_x_0_xzzz_xyyzzz, g_x_0_xzzz_xyyzzzz, g_x_0_xzzz_xyzzzz, g_x_0_xzzz_xyzzzzz, g_x_0_xzzz_xzzzzz, g_x_0_xzzz_yyyyyy, g_x_0_xzzz_yyyyyyy, g_x_0_xzzz_yyyyyyz, g_x_0_xzzz_yyyyyz, g_x_0_xzzz_yyyyyzz, g_x_0_xzzz_yyyyzz, g_x_0_xzzz_yyyyzzz, g_x_0_xzzz_yyyzzz, g_x_0_xzzz_yyyzzzz, g_x_0_xzzz_yyzzzz, g_x_0_xzzz_yyzzzzz, g_x_0_xzzz_yzzzzz, g_x_0_xzzz_yzzzzzz, g_x_0_xzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzz_xxxxxx[k] = -g_x_0_xzzz_xxxxxx[k] * ab_y + g_x_0_xzzz_xxxxxxy[k];

                g_x_0_xyzzz_xxxxxy[k] = -g_x_0_xzzz_xxxxxy[k] * ab_y + g_x_0_xzzz_xxxxxyy[k];

                g_x_0_xyzzz_xxxxxz[k] = -g_x_0_xzzz_xxxxxz[k] * ab_y + g_x_0_xzzz_xxxxxyz[k];

                g_x_0_xyzzz_xxxxyy[k] = -g_x_0_xzzz_xxxxyy[k] * ab_y + g_x_0_xzzz_xxxxyyy[k];

                g_x_0_xyzzz_xxxxyz[k] = -g_x_0_xzzz_xxxxyz[k] * ab_y + g_x_0_xzzz_xxxxyyz[k];

                g_x_0_xyzzz_xxxxzz[k] = -g_x_0_xzzz_xxxxzz[k] * ab_y + g_x_0_xzzz_xxxxyzz[k];

                g_x_0_xyzzz_xxxyyy[k] = -g_x_0_xzzz_xxxyyy[k] * ab_y + g_x_0_xzzz_xxxyyyy[k];

                g_x_0_xyzzz_xxxyyz[k] = -g_x_0_xzzz_xxxyyz[k] * ab_y + g_x_0_xzzz_xxxyyyz[k];

                g_x_0_xyzzz_xxxyzz[k] = -g_x_0_xzzz_xxxyzz[k] * ab_y + g_x_0_xzzz_xxxyyzz[k];

                g_x_0_xyzzz_xxxzzz[k] = -g_x_0_xzzz_xxxzzz[k] * ab_y + g_x_0_xzzz_xxxyzzz[k];

                g_x_0_xyzzz_xxyyyy[k] = -g_x_0_xzzz_xxyyyy[k] * ab_y + g_x_0_xzzz_xxyyyyy[k];

                g_x_0_xyzzz_xxyyyz[k] = -g_x_0_xzzz_xxyyyz[k] * ab_y + g_x_0_xzzz_xxyyyyz[k];

                g_x_0_xyzzz_xxyyzz[k] = -g_x_0_xzzz_xxyyzz[k] * ab_y + g_x_0_xzzz_xxyyyzz[k];

                g_x_0_xyzzz_xxyzzz[k] = -g_x_0_xzzz_xxyzzz[k] * ab_y + g_x_0_xzzz_xxyyzzz[k];

                g_x_0_xyzzz_xxzzzz[k] = -g_x_0_xzzz_xxzzzz[k] * ab_y + g_x_0_xzzz_xxyzzzz[k];

                g_x_0_xyzzz_xyyyyy[k] = -g_x_0_xzzz_xyyyyy[k] * ab_y + g_x_0_xzzz_xyyyyyy[k];

                g_x_0_xyzzz_xyyyyz[k] = -g_x_0_xzzz_xyyyyz[k] * ab_y + g_x_0_xzzz_xyyyyyz[k];

                g_x_0_xyzzz_xyyyzz[k] = -g_x_0_xzzz_xyyyzz[k] * ab_y + g_x_0_xzzz_xyyyyzz[k];

                g_x_0_xyzzz_xyyzzz[k] = -g_x_0_xzzz_xyyzzz[k] * ab_y + g_x_0_xzzz_xyyyzzz[k];

                g_x_0_xyzzz_xyzzzz[k] = -g_x_0_xzzz_xyzzzz[k] * ab_y + g_x_0_xzzz_xyyzzzz[k];

                g_x_0_xyzzz_xzzzzz[k] = -g_x_0_xzzz_xzzzzz[k] * ab_y + g_x_0_xzzz_xyzzzzz[k];

                g_x_0_xyzzz_yyyyyy[k] = -g_x_0_xzzz_yyyyyy[k] * ab_y + g_x_0_xzzz_yyyyyyy[k];

                g_x_0_xyzzz_yyyyyz[k] = -g_x_0_xzzz_yyyyyz[k] * ab_y + g_x_0_xzzz_yyyyyyz[k];

                g_x_0_xyzzz_yyyyzz[k] = -g_x_0_xzzz_yyyyzz[k] * ab_y + g_x_0_xzzz_yyyyyzz[k];

                g_x_0_xyzzz_yyyzzz[k] = -g_x_0_xzzz_yyyzzz[k] * ab_y + g_x_0_xzzz_yyyyzzz[k];

                g_x_0_xyzzz_yyzzzz[k] = -g_x_0_xzzz_yyzzzz[k] * ab_y + g_x_0_xzzz_yyyzzzz[k];

                g_x_0_xyzzz_yzzzzz[k] = -g_x_0_xzzz_yzzzzz[k] * ab_y + g_x_0_xzzz_yyzzzzz[k];

                g_x_0_xyzzz_zzzzzz[k] = -g_x_0_xzzz_zzzzzz[k] * ab_y + g_x_0_xzzz_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 392 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 395 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 396 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 398 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 401 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 404 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 405 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 410 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 416 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_xzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_xzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzzz_xxxxxx, g_x_0_xzzz_xxxxxxz, g_x_0_xzzz_xxxxxy, g_x_0_xzzz_xxxxxyz, g_x_0_xzzz_xxxxxz, g_x_0_xzzz_xxxxxzz, g_x_0_xzzz_xxxxyy, g_x_0_xzzz_xxxxyyz, g_x_0_xzzz_xxxxyz, g_x_0_xzzz_xxxxyzz, g_x_0_xzzz_xxxxzz, g_x_0_xzzz_xxxxzzz, g_x_0_xzzz_xxxyyy, g_x_0_xzzz_xxxyyyz, g_x_0_xzzz_xxxyyz, g_x_0_xzzz_xxxyyzz, g_x_0_xzzz_xxxyzz, g_x_0_xzzz_xxxyzzz, g_x_0_xzzz_xxxzzz, g_x_0_xzzz_xxxzzzz, g_x_0_xzzz_xxyyyy, g_x_0_xzzz_xxyyyyz, g_x_0_xzzz_xxyyyz, g_x_0_xzzz_xxyyyzz, g_x_0_xzzz_xxyyzz, g_x_0_xzzz_xxyyzzz, g_x_0_xzzz_xxyzzz, g_x_0_xzzz_xxyzzzz, g_x_0_xzzz_xxzzzz, g_x_0_xzzz_xxzzzzz, g_x_0_xzzz_xyyyyy, g_x_0_xzzz_xyyyyyz, g_x_0_xzzz_xyyyyz, g_x_0_xzzz_xyyyyzz, g_x_0_xzzz_xyyyzz, g_x_0_xzzz_xyyyzzz, g_x_0_xzzz_xyyzzz, g_x_0_xzzz_xyyzzzz, g_x_0_xzzz_xyzzzz, g_x_0_xzzz_xyzzzzz, g_x_0_xzzz_xzzzzz, g_x_0_xzzz_xzzzzzz, g_x_0_xzzz_yyyyyy, g_x_0_xzzz_yyyyyyz, g_x_0_xzzz_yyyyyz, g_x_0_xzzz_yyyyyzz, g_x_0_xzzz_yyyyzz, g_x_0_xzzz_yyyyzzz, g_x_0_xzzz_yyyzzz, g_x_0_xzzz_yyyzzzz, g_x_0_xzzz_yyzzzz, g_x_0_xzzz_yyzzzzz, g_x_0_xzzz_yzzzzz, g_x_0_xzzz_yzzzzzz, g_x_0_xzzz_zzzzzz, g_x_0_xzzz_zzzzzzz, g_x_0_xzzzz_xxxxxx, g_x_0_xzzzz_xxxxxy, g_x_0_xzzzz_xxxxxz, g_x_0_xzzzz_xxxxyy, g_x_0_xzzzz_xxxxyz, g_x_0_xzzzz_xxxxzz, g_x_0_xzzzz_xxxyyy, g_x_0_xzzzz_xxxyyz, g_x_0_xzzzz_xxxyzz, g_x_0_xzzzz_xxxzzz, g_x_0_xzzzz_xxyyyy, g_x_0_xzzzz_xxyyyz, g_x_0_xzzzz_xxyyzz, g_x_0_xzzzz_xxyzzz, g_x_0_xzzzz_xxzzzz, g_x_0_xzzzz_xyyyyy, g_x_0_xzzzz_xyyyyz, g_x_0_xzzzz_xyyyzz, g_x_0_xzzzz_xyyzzz, g_x_0_xzzzz_xyzzzz, g_x_0_xzzzz_xzzzzz, g_x_0_xzzzz_yyyyyy, g_x_0_xzzzz_yyyyyz, g_x_0_xzzzz_yyyyzz, g_x_0_xzzzz_yyyzzz, g_x_0_xzzzz_yyzzzz, g_x_0_xzzzz_yzzzzz, g_x_0_xzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzz_xxxxxx[k] = -g_x_0_xzzz_xxxxxx[k] * ab_z + g_x_0_xzzz_xxxxxxz[k];

                g_x_0_xzzzz_xxxxxy[k] = -g_x_0_xzzz_xxxxxy[k] * ab_z + g_x_0_xzzz_xxxxxyz[k];

                g_x_0_xzzzz_xxxxxz[k] = -g_x_0_xzzz_xxxxxz[k] * ab_z + g_x_0_xzzz_xxxxxzz[k];

                g_x_0_xzzzz_xxxxyy[k] = -g_x_0_xzzz_xxxxyy[k] * ab_z + g_x_0_xzzz_xxxxyyz[k];

                g_x_0_xzzzz_xxxxyz[k] = -g_x_0_xzzz_xxxxyz[k] * ab_z + g_x_0_xzzz_xxxxyzz[k];

                g_x_0_xzzzz_xxxxzz[k] = -g_x_0_xzzz_xxxxzz[k] * ab_z + g_x_0_xzzz_xxxxzzz[k];

                g_x_0_xzzzz_xxxyyy[k] = -g_x_0_xzzz_xxxyyy[k] * ab_z + g_x_0_xzzz_xxxyyyz[k];

                g_x_0_xzzzz_xxxyyz[k] = -g_x_0_xzzz_xxxyyz[k] * ab_z + g_x_0_xzzz_xxxyyzz[k];

                g_x_0_xzzzz_xxxyzz[k] = -g_x_0_xzzz_xxxyzz[k] * ab_z + g_x_0_xzzz_xxxyzzz[k];

                g_x_0_xzzzz_xxxzzz[k] = -g_x_0_xzzz_xxxzzz[k] * ab_z + g_x_0_xzzz_xxxzzzz[k];

                g_x_0_xzzzz_xxyyyy[k] = -g_x_0_xzzz_xxyyyy[k] * ab_z + g_x_0_xzzz_xxyyyyz[k];

                g_x_0_xzzzz_xxyyyz[k] = -g_x_0_xzzz_xxyyyz[k] * ab_z + g_x_0_xzzz_xxyyyzz[k];

                g_x_0_xzzzz_xxyyzz[k] = -g_x_0_xzzz_xxyyzz[k] * ab_z + g_x_0_xzzz_xxyyzzz[k];

                g_x_0_xzzzz_xxyzzz[k] = -g_x_0_xzzz_xxyzzz[k] * ab_z + g_x_0_xzzz_xxyzzzz[k];

                g_x_0_xzzzz_xxzzzz[k] = -g_x_0_xzzz_xxzzzz[k] * ab_z + g_x_0_xzzz_xxzzzzz[k];

                g_x_0_xzzzz_xyyyyy[k] = -g_x_0_xzzz_xyyyyy[k] * ab_z + g_x_0_xzzz_xyyyyyz[k];

                g_x_0_xzzzz_xyyyyz[k] = -g_x_0_xzzz_xyyyyz[k] * ab_z + g_x_0_xzzz_xyyyyzz[k];

                g_x_0_xzzzz_xyyyzz[k] = -g_x_0_xzzz_xyyyzz[k] * ab_z + g_x_0_xzzz_xyyyzzz[k];

                g_x_0_xzzzz_xyyzzz[k] = -g_x_0_xzzz_xyyzzz[k] * ab_z + g_x_0_xzzz_xyyzzzz[k];

                g_x_0_xzzzz_xyzzzz[k] = -g_x_0_xzzz_xyzzzz[k] * ab_z + g_x_0_xzzz_xyzzzzz[k];

                g_x_0_xzzzz_xzzzzz[k] = -g_x_0_xzzz_xzzzzz[k] * ab_z + g_x_0_xzzz_xzzzzzz[k];

                g_x_0_xzzzz_yyyyyy[k] = -g_x_0_xzzz_yyyyyy[k] * ab_z + g_x_0_xzzz_yyyyyyz[k];

                g_x_0_xzzzz_yyyyyz[k] = -g_x_0_xzzz_yyyyyz[k] * ab_z + g_x_0_xzzz_yyyyyzz[k];

                g_x_0_xzzzz_yyyyzz[k] = -g_x_0_xzzz_yyyyzz[k] * ab_z + g_x_0_xzzz_yyyyzzz[k];

                g_x_0_xzzzz_yyyzzz[k] = -g_x_0_xzzz_yyyzzz[k] * ab_z + g_x_0_xzzz_yyyzzzz[k];

                g_x_0_xzzzz_yyzzzz[k] = -g_x_0_xzzz_yyzzzz[k] * ab_z + g_x_0_xzzz_yyzzzzz[k];

                g_x_0_xzzzz_yzzzzz[k] = -g_x_0_xzzz_yzzzzz[k] * ab_z + g_x_0_xzzz_yzzzzzz[k];

                g_x_0_xzzzz_zzzzzz[k] = -g_x_0_xzzz_zzzzzz[k] * ab_z + g_x_0_xzzz_zzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 420 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 421 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 422 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 423 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 424 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 425 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 426 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 427 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 428 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 429 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 430 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 431 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 432 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 433 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 434 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 435 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 436 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 437 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 438 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 439 * ccomps * dcomps);

            auto g_x_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 440 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 441 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 442 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 443 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 444 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 445 * ccomps * dcomps);

            auto g_x_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 446 * ccomps * dcomps);

            auto g_x_0_yyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 447 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyy_xxxxxx, g_x_0_yyyy_xxxxxxy, g_x_0_yyyy_xxxxxy, g_x_0_yyyy_xxxxxyy, g_x_0_yyyy_xxxxxyz, g_x_0_yyyy_xxxxxz, g_x_0_yyyy_xxxxyy, g_x_0_yyyy_xxxxyyy, g_x_0_yyyy_xxxxyyz, g_x_0_yyyy_xxxxyz, g_x_0_yyyy_xxxxyzz, g_x_0_yyyy_xxxxzz, g_x_0_yyyy_xxxyyy, g_x_0_yyyy_xxxyyyy, g_x_0_yyyy_xxxyyyz, g_x_0_yyyy_xxxyyz, g_x_0_yyyy_xxxyyzz, g_x_0_yyyy_xxxyzz, g_x_0_yyyy_xxxyzzz, g_x_0_yyyy_xxxzzz, g_x_0_yyyy_xxyyyy, g_x_0_yyyy_xxyyyyy, g_x_0_yyyy_xxyyyyz, g_x_0_yyyy_xxyyyz, g_x_0_yyyy_xxyyyzz, g_x_0_yyyy_xxyyzz, g_x_0_yyyy_xxyyzzz, g_x_0_yyyy_xxyzzz, g_x_0_yyyy_xxyzzzz, g_x_0_yyyy_xxzzzz, g_x_0_yyyy_xyyyyy, g_x_0_yyyy_xyyyyyy, g_x_0_yyyy_xyyyyyz, g_x_0_yyyy_xyyyyz, g_x_0_yyyy_xyyyyzz, g_x_0_yyyy_xyyyzz, g_x_0_yyyy_xyyyzzz, g_x_0_yyyy_xyyzzz, g_x_0_yyyy_xyyzzzz, g_x_0_yyyy_xyzzzz, g_x_0_yyyy_xyzzzzz, g_x_0_yyyy_xzzzzz, g_x_0_yyyy_yyyyyy, g_x_0_yyyy_yyyyyyy, g_x_0_yyyy_yyyyyyz, g_x_0_yyyy_yyyyyz, g_x_0_yyyy_yyyyyzz, g_x_0_yyyy_yyyyzz, g_x_0_yyyy_yyyyzzz, g_x_0_yyyy_yyyzzz, g_x_0_yyyy_yyyzzzz, g_x_0_yyyy_yyzzzz, g_x_0_yyyy_yyzzzzz, g_x_0_yyyy_yzzzzz, g_x_0_yyyy_yzzzzzz, g_x_0_yyyy_zzzzzz, g_x_0_yyyyy_xxxxxx, g_x_0_yyyyy_xxxxxy, g_x_0_yyyyy_xxxxxz, g_x_0_yyyyy_xxxxyy, g_x_0_yyyyy_xxxxyz, g_x_0_yyyyy_xxxxzz, g_x_0_yyyyy_xxxyyy, g_x_0_yyyyy_xxxyyz, g_x_0_yyyyy_xxxyzz, g_x_0_yyyyy_xxxzzz, g_x_0_yyyyy_xxyyyy, g_x_0_yyyyy_xxyyyz, g_x_0_yyyyy_xxyyzz, g_x_0_yyyyy_xxyzzz, g_x_0_yyyyy_xxzzzz, g_x_0_yyyyy_xyyyyy, g_x_0_yyyyy_xyyyyz, g_x_0_yyyyy_xyyyzz, g_x_0_yyyyy_xyyzzz, g_x_0_yyyyy_xyzzzz, g_x_0_yyyyy_xzzzzz, g_x_0_yyyyy_yyyyyy, g_x_0_yyyyy_yyyyyz, g_x_0_yyyyy_yyyyzz, g_x_0_yyyyy_yyyzzz, g_x_0_yyyyy_yyzzzz, g_x_0_yyyyy_yzzzzz, g_x_0_yyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyy_xxxxxx[k] = -g_x_0_yyyy_xxxxxx[k] * ab_y + g_x_0_yyyy_xxxxxxy[k];

                g_x_0_yyyyy_xxxxxy[k] = -g_x_0_yyyy_xxxxxy[k] * ab_y + g_x_0_yyyy_xxxxxyy[k];

                g_x_0_yyyyy_xxxxxz[k] = -g_x_0_yyyy_xxxxxz[k] * ab_y + g_x_0_yyyy_xxxxxyz[k];

                g_x_0_yyyyy_xxxxyy[k] = -g_x_0_yyyy_xxxxyy[k] * ab_y + g_x_0_yyyy_xxxxyyy[k];

                g_x_0_yyyyy_xxxxyz[k] = -g_x_0_yyyy_xxxxyz[k] * ab_y + g_x_0_yyyy_xxxxyyz[k];

                g_x_0_yyyyy_xxxxzz[k] = -g_x_0_yyyy_xxxxzz[k] * ab_y + g_x_0_yyyy_xxxxyzz[k];

                g_x_0_yyyyy_xxxyyy[k] = -g_x_0_yyyy_xxxyyy[k] * ab_y + g_x_0_yyyy_xxxyyyy[k];

                g_x_0_yyyyy_xxxyyz[k] = -g_x_0_yyyy_xxxyyz[k] * ab_y + g_x_0_yyyy_xxxyyyz[k];

                g_x_0_yyyyy_xxxyzz[k] = -g_x_0_yyyy_xxxyzz[k] * ab_y + g_x_0_yyyy_xxxyyzz[k];

                g_x_0_yyyyy_xxxzzz[k] = -g_x_0_yyyy_xxxzzz[k] * ab_y + g_x_0_yyyy_xxxyzzz[k];

                g_x_0_yyyyy_xxyyyy[k] = -g_x_0_yyyy_xxyyyy[k] * ab_y + g_x_0_yyyy_xxyyyyy[k];

                g_x_0_yyyyy_xxyyyz[k] = -g_x_0_yyyy_xxyyyz[k] * ab_y + g_x_0_yyyy_xxyyyyz[k];

                g_x_0_yyyyy_xxyyzz[k] = -g_x_0_yyyy_xxyyzz[k] * ab_y + g_x_0_yyyy_xxyyyzz[k];

                g_x_0_yyyyy_xxyzzz[k] = -g_x_0_yyyy_xxyzzz[k] * ab_y + g_x_0_yyyy_xxyyzzz[k];

                g_x_0_yyyyy_xxzzzz[k] = -g_x_0_yyyy_xxzzzz[k] * ab_y + g_x_0_yyyy_xxyzzzz[k];

                g_x_0_yyyyy_xyyyyy[k] = -g_x_0_yyyy_xyyyyy[k] * ab_y + g_x_0_yyyy_xyyyyyy[k];

                g_x_0_yyyyy_xyyyyz[k] = -g_x_0_yyyy_xyyyyz[k] * ab_y + g_x_0_yyyy_xyyyyyz[k];

                g_x_0_yyyyy_xyyyzz[k] = -g_x_0_yyyy_xyyyzz[k] * ab_y + g_x_0_yyyy_xyyyyzz[k];

                g_x_0_yyyyy_xyyzzz[k] = -g_x_0_yyyy_xyyzzz[k] * ab_y + g_x_0_yyyy_xyyyzzz[k];

                g_x_0_yyyyy_xyzzzz[k] = -g_x_0_yyyy_xyzzzz[k] * ab_y + g_x_0_yyyy_xyyzzzz[k];

                g_x_0_yyyyy_xzzzzz[k] = -g_x_0_yyyy_xzzzzz[k] * ab_y + g_x_0_yyyy_xyzzzzz[k];

                g_x_0_yyyyy_yyyyyy[k] = -g_x_0_yyyy_yyyyyy[k] * ab_y + g_x_0_yyyy_yyyyyyy[k];

                g_x_0_yyyyy_yyyyyz[k] = -g_x_0_yyyy_yyyyyz[k] * ab_y + g_x_0_yyyy_yyyyyyz[k];

                g_x_0_yyyyy_yyyyzz[k] = -g_x_0_yyyy_yyyyzz[k] * ab_y + g_x_0_yyyy_yyyyyzz[k];

                g_x_0_yyyyy_yyyzzz[k] = -g_x_0_yyyy_yyyzzz[k] * ab_y + g_x_0_yyyy_yyyyzzz[k];

                g_x_0_yyyyy_yyzzzz[k] = -g_x_0_yyyy_yyzzzz[k] * ab_y + g_x_0_yyyy_yyyzzzz[k];

                g_x_0_yyyyy_yzzzzz[k] = -g_x_0_yyyy_yzzzzz[k] * ab_y + g_x_0_yyyy_yyzzzzz[k];

                g_x_0_yyyyy_zzzzzz[k] = -g_x_0_yyyy_zzzzzz[k] * ab_y + g_x_0_yyyy_yzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 448 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 449 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 450 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 451 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 452 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 453 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 454 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 455 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 456 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 457 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 458 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 459 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 460 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 461 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 462 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 463 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 464 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 465 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 466 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 467 * ccomps * dcomps);

            auto g_x_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 468 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 469 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 470 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 471 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 472 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 473 * ccomps * dcomps);

            auto g_x_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 474 * ccomps * dcomps);

            auto g_x_0_yyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 475 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyz_xxxxxx, g_x_0_yyyyz_xxxxxy, g_x_0_yyyyz_xxxxxz, g_x_0_yyyyz_xxxxyy, g_x_0_yyyyz_xxxxyz, g_x_0_yyyyz_xxxxzz, g_x_0_yyyyz_xxxyyy, g_x_0_yyyyz_xxxyyz, g_x_0_yyyyz_xxxyzz, g_x_0_yyyyz_xxxzzz, g_x_0_yyyyz_xxyyyy, g_x_0_yyyyz_xxyyyz, g_x_0_yyyyz_xxyyzz, g_x_0_yyyyz_xxyzzz, g_x_0_yyyyz_xxzzzz, g_x_0_yyyyz_xyyyyy, g_x_0_yyyyz_xyyyyz, g_x_0_yyyyz_xyyyzz, g_x_0_yyyyz_xyyzzz, g_x_0_yyyyz_xyzzzz, g_x_0_yyyyz_xzzzzz, g_x_0_yyyyz_yyyyyy, g_x_0_yyyyz_yyyyyz, g_x_0_yyyyz_yyyyzz, g_x_0_yyyyz_yyyzzz, g_x_0_yyyyz_yyzzzz, g_x_0_yyyyz_yzzzzz, g_x_0_yyyyz_zzzzzz, g_x_0_yyyz_xxxxxx, g_x_0_yyyz_xxxxxxy, g_x_0_yyyz_xxxxxy, g_x_0_yyyz_xxxxxyy, g_x_0_yyyz_xxxxxyz, g_x_0_yyyz_xxxxxz, g_x_0_yyyz_xxxxyy, g_x_0_yyyz_xxxxyyy, g_x_0_yyyz_xxxxyyz, g_x_0_yyyz_xxxxyz, g_x_0_yyyz_xxxxyzz, g_x_0_yyyz_xxxxzz, g_x_0_yyyz_xxxyyy, g_x_0_yyyz_xxxyyyy, g_x_0_yyyz_xxxyyyz, g_x_0_yyyz_xxxyyz, g_x_0_yyyz_xxxyyzz, g_x_0_yyyz_xxxyzz, g_x_0_yyyz_xxxyzzz, g_x_0_yyyz_xxxzzz, g_x_0_yyyz_xxyyyy, g_x_0_yyyz_xxyyyyy, g_x_0_yyyz_xxyyyyz, g_x_0_yyyz_xxyyyz, g_x_0_yyyz_xxyyyzz, g_x_0_yyyz_xxyyzz, g_x_0_yyyz_xxyyzzz, g_x_0_yyyz_xxyzzz, g_x_0_yyyz_xxyzzzz, g_x_0_yyyz_xxzzzz, g_x_0_yyyz_xyyyyy, g_x_0_yyyz_xyyyyyy, g_x_0_yyyz_xyyyyyz, g_x_0_yyyz_xyyyyz, g_x_0_yyyz_xyyyyzz, g_x_0_yyyz_xyyyzz, g_x_0_yyyz_xyyyzzz, g_x_0_yyyz_xyyzzz, g_x_0_yyyz_xyyzzzz, g_x_0_yyyz_xyzzzz, g_x_0_yyyz_xyzzzzz, g_x_0_yyyz_xzzzzz, g_x_0_yyyz_yyyyyy, g_x_0_yyyz_yyyyyyy, g_x_0_yyyz_yyyyyyz, g_x_0_yyyz_yyyyyz, g_x_0_yyyz_yyyyyzz, g_x_0_yyyz_yyyyzz, g_x_0_yyyz_yyyyzzz, g_x_0_yyyz_yyyzzz, g_x_0_yyyz_yyyzzzz, g_x_0_yyyz_yyzzzz, g_x_0_yyyz_yyzzzzz, g_x_0_yyyz_yzzzzz, g_x_0_yyyz_yzzzzzz, g_x_0_yyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyz_xxxxxx[k] = -g_x_0_yyyz_xxxxxx[k] * ab_y + g_x_0_yyyz_xxxxxxy[k];

                g_x_0_yyyyz_xxxxxy[k] = -g_x_0_yyyz_xxxxxy[k] * ab_y + g_x_0_yyyz_xxxxxyy[k];

                g_x_0_yyyyz_xxxxxz[k] = -g_x_0_yyyz_xxxxxz[k] * ab_y + g_x_0_yyyz_xxxxxyz[k];

                g_x_0_yyyyz_xxxxyy[k] = -g_x_0_yyyz_xxxxyy[k] * ab_y + g_x_0_yyyz_xxxxyyy[k];

                g_x_0_yyyyz_xxxxyz[k] = -g_x_0_yyyz_xxxxyz[k] * ab_y + g_x_0_yyyz_xxxxyyz[k];

                g_x_0_yyyyz_xxxxzz[k] = -g_x_0_yyyz_xxxxzz[k] * ab_y + g_x_0_yyyz_xxxxyzz[k];

                g_x_0_yyyyz_xxxyyy[k] = -g_x_0_yyyz_xxxyyy[k] * ab_y + g_x_0_yyyz_xxxyyyy[k];

                g_x_0_yyyyz_xxxyyz[k] = -g_x_0_yyyz_xxxyyz[k] * ab_y + g_x_0_yyyz_xxxyyyz[k];

                g_x_0_yyyyz_xxxyzz[k] = -g_x_0_yyyz_xxxyzz[k] * ab_y + g_x_0_yyyz_xxxyyzz[k];

                g_x_0_yyyyz_xxxzzz[k] = -g_x_0_yyyz_xxxzzz[k] * ab_y + g_x_0_yyyz_xxxyzzz[k];

                g_x_0_yyyyz_xxyyyy[k] = -g_x_0_yyyz_xxyyyy[k] * ab_y + g_x_0_yyyz_xxyyyyy[k];

                g_x_0_yyyyz_xxyyyz[k] = -g_x_0_yyyz_xxyyyz[k] * ab_y + g_x_0_yyyz_xxyyyyz[k];

                g_x_0_yyyyz_xxyyzz[k] = -g_x_0_yyyz_xxyyzz[k] * ab_y + g_x_0_yyyz_xxyyyzz[k];

                g_x_0_yyyyz_xxyzzz[k] = -g_x_0_yyyz_xxyzzz[k] * ab_y + g_x_0_yyyz_xxyyzzz[k];

                g_x_0_yyyyz_xxzzzz[k] = -g_x_0_yyyz_xxzzzz[k] * ab_y + g_x_0_yyyz_xxyzzzz[k];

                g_x_0_yyyyz_xyyyyy[k] = -g_x_0_yyyz_xyyyyy[k] * ab_y + g_x_0_yyyz_xyyyyyy[k];

                g_x_0_yyyyz_xyyyyz[k] = -g_x_0_yyyz_xyyyyz[k] * ab_y + g_x_0_yyyz_xyyyyyz[k];

                g_x_0_yyyyz_xyyyzz[k] = -g_x_0_yyyz_xyyyzz[k] * ab_y + g_x_0_yyyz_xyyyyzz[k];

                g_x_0_yyyyz_xyyzzz[k] = -g_x_0_yyyz_xyyzzz[k] * ab_y + g_x_0_yyyz_xyyyzzz[k];

                g_x_0_yyyyz_xyzzzz[k] = -g_x_0_yyyz_xyzzzz[k] * ab_y + g_x_0_yyyz_xyyzzzz[k];

                g_x_0_yyyyz_xzzzzz[k] = -g_x_0_yyyz_xzzzzz[k] * ab_y + g_x_0_yyyz_xyzzzzz[k];

                g_x_0_yyyyz_yyyyyy[k] = -g_x_0_yyyz_yyyyyy[k] * ab_y + g_x_0_yyyz_yyyyyyy[k];

                g_x_0_yyyyz_yyyyyz[k] = -g_x_0_yyyz_yyyyyz[k] * ab_y + g_x_0_yyyz_yyyyyyz[k];

                g_x_0_yyyyz_yyyyzz[k] = -g_x_0_yyyz_yyyyzz[k] * ab_y + g_x_0_yyyz_yyyyyzz[k];

                g_x_0_yyyyz_yyyzzz[k] = -g_x_0_yyyz_yyyzzz[k] * ab_y + g_x_0_yyyz_yyyyzzz[k];

                g_x_0_yyyyz_yyzzzz[k] = -g_x_0_yyyz_yyzzzz[k] * ab_y + g_x_0_yyyz_yyyzzzz[k];

                g_x_0_yyyyz_yzzzzz[k] = -g_x_0_yyyz_yzzzzz[k] * ab_y + g_x_0_yyyz_yyzzzzz[k];

                g_x_0_yyyyz_zzzzzz[k] = -g_x_0_yyyz_zzzzzz[k] * ab_y + g_x_0_yyyz_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 476 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 477 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 478 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 479 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 480 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 481 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 482 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 483 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 484 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 485 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 486 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 487 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 488 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 489 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 490 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 491 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 492 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 493 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 494 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 495 * ccomps * dcomps);

            auto g_x_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 496 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 497 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 498 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 499 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 500 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 501 * ccomps * dcomps);

            auto g_x_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 502 * ccomps * dcomps);

            auto g_x_0_yyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyzz_xxxxxx, g_x_0_yyyzz_xxxxxy, g_x_0_yyyzz_xxxxxz, g_x_0_yyyzz_xxxxyy, g_x_0_yyyzz_xxxxyz, g_x_0_yyyzz_xxxxzz, g_x_0_yyyzz_xxxyyy, g_x_0_yyyzz_xxxyyz, g_x_0_yyyzz_xxxyzz, g_x_0_yyyzz_xxxzzz, g_x_0_yyyzz_xxyyyy, g_x_0_yyyzz_xxyyyz, g_x_0_yyyzz_xxyyzz, g_x_0_yyyzz_xxyzzz, g_x_0_yyyzz_xxzzzz, g_x_0_yyyzz_xyyyyy, g_x_0_yyyzz_xyyyyz, g_x_0_yyyzz_xyyyzz, g_x_0_yyyzz_xyyzzz, g_x_0_yyyzz_xyzzzz, g_x_0_yyyzz_xzzzzz, g_x_0_yyyzz_yyyyyy, g_x_0_yyyzz_yyyyyz, g_x_0_yyyzz_yyyyzz, g_x_0_yyyzz_yyyzzz, g_x_0_yyyzz_yyzzzz, g_x_0_yyyzz_yzzzzz, g_x_0_yyyzz_zzzzzz, g_x_0_yyzz_xxxxxx, g_x_0_yyzz_xxxxxxy, g_x_0_yyzz_xxxxxy, g_x_0_yyzz_xxxxxyy, g_x_0_yyzz_xxxxxyz, g_x_0_yyzz_xxxxxz, g_x_0_yyzz_xxxxyy, g_x_0_yyzz_xxxxyyy, g_x_0_yyzz_xxxxyyz, g_x_0_yyzz_xxxxyz, g_x_0_yyzz_xxxxyzz, g_x_0_yyzz_xxxxzz, g_x_0_yyzz_xxxyyy, g_x_0_yyzz_xxxyyyy, g_x_0_yyzz_xxxyyyz, g_x_0_yyzz_xxxyyz, g_x_0_yyzz_xxxyyzz, g_x_0_yyzz_xxxyzz, g_x_0_yyzz_xxxyzzz, g_x_0_yyzz_xxxzzz, g_x_0_yyzz_xxyyyy, g_x_0_yyzz_xxyyyyy, g_x_0_yyzz_xxyyyyz, g_x_0_yyzz_xxyyyz, g_x_0_yyzz_xxyyyzz, g_x_0_yyzz_xxyyzz, g_x_0_yyzz_xxyyzzz, g_x_0_yyzz_xxyzzz, g_x_0_yyzz_xxyzzzz, g_x_0_yyzz_xxzzzz, g_x_0_yyzz_xyyyyy, g_x_0_yyzz_xyyyyyy, g_x_0_yyzz_xyyyyyz, g_x_0_yyzz_xyyyyz, g_x_0_yyzz_xyyyyzz, g_x_0_yyzz_xyyyzz, g_x_0_yyzz_xyyyzzz, g_x_0_yyzz_xyyzzz, g_x_0_yyzz_xyyzzzz, g_x_0_yyzz_xyzzzz, g_x_0_yyzz_xyzzzzz, g_x_0_yyzz_xzzzzz, g_x_0_yyzz_yyyyyy, g_x_0_yyzz_yyyyyyy, g_x_0_yyzz_yyyyyyz, g_x_0_yyzz_yyyyyz, g_x_0_yyzz_yyyyyzz, g_x_0_yyzz_yyyyzz, g_x_0_yyzz_yyyyzzz, g_x_0_yyzz_yyyzzz, g_x_0_yyzz_yyyzzzz, g_x_0_yyzz_yyzzzz, g_x_0_yyzz_yyzzzzz, g_x_0_yyzz_yzzzzz, g_x_0_yyzz_yzzzzzz, g_x_0_yyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzz_xxxxxx[k] = -g_x_0_yyzz_xxxxxx[k] * ab_y + g_x_0_yyzz_xxxxxxy[k];

                g_x_0_yyyzz_xxxxxy[k] = -g_x_0_yyzz_xxxxxy[k] * ab_y + g_x_0_yyzz_xxxxxyy[k];

                g_x_0_yyyzz_xxxxxz[k] = -g_x_0_yyzz_xxxxxz[k] * ab_y + g_x_0_yyzz_xxxxxyz[k];

                g_x_0_yyyzz_xxxxyy[k] = -g_x_0_yyzz_xxxxyy[k] * ab_y + g_x_0_yyzz_xxxxyyy[k];

                g_x_0_yyyzz_xxxxyz[k] = -g_x_0_yyzz_xxxxyz[k] * ab_y + g_x_0_yyzz_xxxxyyz[k];

                g_x_0_yyyzz_xxxxzz[k] = -g_x_0_yyzz_xxxxzz[k] * ab_y + g_x_0_yyzz_xxxxyzz[k];

                g_x_0_yyyzz_xxxyyy[k] = -g_x_0_yyzz_xxxyyy[k] * ab_y + g_x_0_yyzz_xxxyyyy[k];

                g_x_0_yyyzz_xxxyyz[k] = -g_x_0_yyzz_xxxyyz[k] * ab_y + g_x_0_yyzz_xxxyyyz[k];

                g_x_0_yyyzz_xxxyzz[k] = -g_x_0_yyzz_xxxyzz[k] * ab_y + g_x_0_yyzz_xxxyyzz[k];

                g_x_0_yyyzz_xxxzzz[k] = -g_x_0_yyzz_xxxzzz[k] * ab_y + g_x_0_yyzz_xxxyzzz[k];

                g_x_0_yyyzz_xxyyyy[k] = -g_x_0_yyzz_xxyyyy[k] * ab_y + g_x_0_yyzz_xxyyyyy[k];

                g_x_0_yyyzz_xxyyyz[k] = -g_x_0_yyzz_xxyyyz[k] * ab_y + g_x_0_yyzz_xxyyyyz[k];

                g_x_0_yyyzz_xxyyzz[k] = -g_x_0_yyzz_xxyyzz[k] * ab_y + g_x_0_yyzz_xxyyyzz[k];

                g_x_0_yyyzz_xxyzzz[k] = -g_x_0_yyzz_xxyzzz[k] * ab_y + g_x_0_yyzz_xxyyzzz[k];

                g_x_0_yyyzz_xxzzzz[k] = -g_x_0_yyzz_xxzzzz[k] * ab_y + g_x_0_yyzz_xxyzzzz[k];

                g_x_0_yyyzz_xyyyyy[k] = -g_x_0_yyzz_xyyyyy[k] * ab_y + g_x_0_yyzz_xyyyyyy[k];

                g_x_0_yyyzz_xyyyyz[k] = -g_x_0_yyzz_xyyyyz[k] * ab_y + g_x_0_yyzz_xyyyyyz[k];

                g_x_0_yyyzz_xyyyzz[k] = -g_x_0_yyzz_xyyyzz[k] * ab_y + g_x_0_yyzz_xyyyyzz[k];

                g_x_0_yyyzz_xyyzzz[k] = -g_x_0_yyzz_xyyzzz[k] * ab_y + g_x_0_yyzz_xyyyzzz[k];

                g_x_0_yyyzz_xyzzzz[k] = -g_x_0_yyzz_xyzzzz[k] * ab_y + g_x_0_yyzz_xyyzzzz[k];

                g_x_0_yyyzz_xzzzzz[k] = -g_x_0_yyzz_xzzzzz[k] * ab_y + g_x_0_yyzz_xyzzzzz[k];

                g_x_0_yyyzz_yyyyyy[k] = -g_x_0_yyzz_yyyyyy[k] * ab_y + g_x_0_yyzz_yyyyyyy[k];

                g_x_0_yyyzz_yyyyyz[k] = -g_x_0_yyzz_yyyyyz[k] * ab_y + g_x_0_yyzz_yyyyyyz[k];

                g_x_0_yyyzz_yyyyzz[k] = -g_x_0_yyzz_yyyyzz[k] * ab_y + g_x_0_yyzz_yyyyyzz[k];

                g_x_0_yyyzz_yyyzzz[k] = -g_x_0_yyzz_yyyzzz[k] * ab_y + g_x_0_yyzz_yyyyzzz[k];

                g_x_0_yyyzz_yyzzzz[k] = -g_x_0_yyzz_yyzzzz[k] * ab_y + g_x_0_yyzz_yyyzzzz[k];

                g_x_0_yyyzz_yzzzzz[k] = -g_x_0_yyzz_yzzzzz[k] * ab_y + g_x_0_yyzz_yyzzzzz[k];

                g_x_0_yyyzz_zzzzzz[k] = -g_x_0_yyzz_zzzzzz[k] * ab_y + g_x_0_yyzz_yzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 504 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 505 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 506 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 507 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 508 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 509 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 510 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 511 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 512 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 513 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 514 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 515 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 516 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 517 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 518 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 519 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 520 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 521 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 522 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 523 * ccomps * dcomps);

            auto g_x_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 524 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 525 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 526 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 527 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 528 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 529 * ccomps * dcomps);

            auto g_x_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 530 * ccomps * dcomps);

            auto g_x_0_yyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 531 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzzz_xxxxxx, g_x_0_yyzzz_xxxxxy, g_x_0_yyzzz_xxxxxz, g_x_0_yyzzz_xxxxyy, g_x_0_yyzzz_xxxxyz, g_x_0_yyzzz_xxxxzz, g_x_0_yyzzz_xxxyyy, g_x_0_yyzzz_xxxyyz, g_x_0_yyzzz_xxxyzz, g_x_0_yyzzz_xxxzzz, g_x_0_yyzzz_xxyyyy, g_x_0_yyzzz_xxyyyz, g_x_0_yyzzz_xxyyzz, g_x_0_yyzzz_xxyzzz, g_x_0_yyzzz_xxzzzz, g_x_0_yyzzz_xyyyyy, g_x_0_yyzzz_xyyyyz, g_x_0_yyzzz_xyyyzz, g_x_0_yyzzz_xyyzzz, g_x_0_yyzzz_xyzzzz, g_x_0_yyzzz_xzzzzz, g_x_0_yyzzz_yyyyyy, g_x_0_yyzzz_yyyyyz, g_x_0_yyzzz_yyyyzz, g_x_0_yyzzz_yyyzzz, g_x_0_yyzzz_yyzzzz, g_x_0_yyzzz_yzzzzz, g_x_0_yyzzz_zzzzzz, g_x_0_yzzz_xxxxxx, g_x_0_yzzz_xxxxxxy, g_x_0_yzzz_xxxxxy, g_x_0_yzzz_xxxxxyy, g_x_0_yzzz_xxxxxyz, g_x_0_yzzz_xxxxxz, g_x_0_yzzz_xxxxyy, g_x_0_yzzz_xxxxyyy, g_x_0_yzzz_xxxxyyz, g_x_0_yzzz_xxxxyz, g_x_0_yzzz_xxxxyzz, g_x_0_yzzz_xxxxzz, g_x_0_yzzz_xxxyyy, g_x_0_yzzz_xxxyyyy, g_x_0_yzzz_xxxyyyz, g_x_0_yzzz_xxxyyz, g_x_0_yzzz_xxxyyzz, g_x_0_yzzz_xxxyzz, g_x_0_yzzz_xxxyzzz, g_x_0_yzzz_xxxzzz, g_x_0_yzzz_xxyyyy, g_x_0_yzzz_xxyyyyy, g_x_0_yzzz_xxyyyyz, g_x_0_yzzz_xxyyyz, g_x_0_yzzz_xxyyyzz, g_x_0_yzzz_xxyyzz, g_x_0_yzzz_xxyyzzz, g_x_0_yzzz_xxyzzz, g_x_0_yzzz_xxyzzzz, g_x_0_yzzz_xxzzzz, g_x_0_yzzz_xyyyyy, g_x_0_yzzz_xyyyyyy, g_x_0_yzzz_xyyyyyz, g_x_0_yzzz_xyyyyz, g_x_0_yzzz_xyyyyzz, g_x_0_yzzz_xyyyzz, g_x_0_yzzz_xyyyzzz, g_x_0_yzzz_xyyzzz, g_x_0_yzzz_xyyzzzz, g_x_0_yzzz_xyzzzz, g_x_0_yzzz_xyzzzzz, g_x_0_yzzz_xzzzzz, g_x_0_yzzz_yyyyyy, g_x_0_yzzz_yyyyyyy, g_x_0_yzzz_yyyyyyz, g_x_0_yzzz_yyyyyz, g_x_0_yzzz_yyyyyzz, g_x_0_yzzz_yyyyzz, g_x_0_yzzz_yyyyzzz, g_x_0_yzzz_yyyzzz, g_x_0_yzzz_yyyzzzz, g_x_0_yzzz_yyzzzz, g_x_0_yzzz_yyzzzzz, g_x_0_yzzz_yzzzzz, g_x_0_yzzz_yzzzzzz, g_x_0_yzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzz_xxxxxx[k] = -g_x_0_yzzz_xxxxxx[k] * ab_y + g_x_0_yzzz_xxxxxxy[k];

                g_x_0_yyzzz_xxxxxy[k] = -g_x_0_yzzz_xxxxxy[k] * ab_y + g_x_0_yzzz_xxxxxyy[k];

                g_x_0_yyzzz_xxxxxz[k] = -g_x_0_yzzz_xxxxxz[k] * ab_y + g_x_0_yzzz_xxxxxyz[k];

                g_x_0_yyzzz_xxxxyy[k] = -g_x_0_yzzz_xxxxyy[k] * ab_y + g_x_0_yzzz_xxxxyyy[k];

                g_x_0_yyzzz_xxxxyz[k] = -g_x_0_yzzz_xxxxyz[k] * ab_y + g_x_0_yzzz_xxxxyyz[k];

                g_x_0_yyzzz_xxxxzz[k] = -g_x_0_yzzz_xxxxzz[k] * ab_y + g_x_0_yzzz_xxxxyzz[k];

                g_x_0_yyzzz_xxxyyy[k] = -g_x_0_yzzz_xxxyyy[k] * ab_y + g_x_0_yzzz_xxxyyyy[k];

                g_x_0_yyzzz_xxxyyz[k] = -g_x_0_yzzz_xxxyyz[k] * ab_y + g_x_0_yzzz_xxxyyyz[k];

                g_x_0_yyzzz_xxxyzz[k] = -g_x_0_yzzz_xxxyzz[k] * ab_y + g_x_0_yzzz_xxxyyzz[k];

                g_x_0_yyzzz_xxxzzz[k] = -g_x_0_yzzz_xxxzzz[k] * ab_y + g_x_0_yzzz_xxxyzzz[k];

                g_x_0_yyzzz_xxyyyy[k] = -g_x_0_yzzz_xxyyyy[k] * ab_y + g_x_0_yzzz_xxyyyyy[k];

                g_x_0_yyzzz_xxyyyz[k] = -g_x_0_yzzz_xxyyyz[k] * ab_y + g_x_0_yzzz_xxyyyyz[k];

                g_x_0_yyzzz_xxyyzz[k] = -g_x_0_yzzz_xxyyzz[k] * ab_y + g_x_0_yzzz_xxyyyzz[k];

                g_x_0_yyzzz_xxyzzz[k] = -g_x_0_yzzz_xxyzzz[k] * ab_y + g_x_0_yzzz_xxyyzzz[k];

                g_x_0_yyzzz_xxzzzz[k] = -g_x_0_yzzz_xxzzzz[k] * ab_y + g_x_0_yzzz_xxyzzzz[k];

                g_x_0_yyzzz_xyyyyy[k] = -g_x_0_yzzz_xyyyyy[k] * ab_y + g_x_0_yzzz_xyyyyyy[k];

                g_x_0_yyzzz_xyyyyz[k] = -g_x_0_yzzz_xyyyyz[k] * ab_y + g_x_0_yzzz_xyyyyyz[k];

                g_x_0_yyzzz_xyyyzz[k] = -g_x_0_yzzz_xyyyzz[k] * ab_y + g_x_0_yzzz_xyyyyzz[k];

                g_x_0_yyzzz_xyyzzz[k] = -g_x_0_yzzz_xyyzzz[k] * ab_y + g_x_0_yzzz_xyyyzzz[k];

                g_x_0_yyzzz_xyzzzz[k] = -g_x_0_yzzz_xyzzzz[k] * ab_y + g_x_0_yzzz_xyyzzzz[k];

                g_x_0_yyzzz_xzzzzz[k] = -g_x_0_yzzz_xzzzzz[k] * ab_y + g_x_0_yzzz_xyzzzzz[k];

                g_x_0_yyzzz_yyyyyy[k] = -g_x_0_yzzz_yyyyyy[k] * ab_y + g_x_0_yzzz_yyyyyyy[k];

                g_x_0_yyzzz_yyyyyz[k] = -g_x_0_yzzz_yyyyyz[k] * ab_y + g_x_0_yzzz_yyyyyyz[k];

                g_x_0_yyzzz_yyyyzz[k] = -g_x_0_yzzz_yyyyzz[k] * ab_y + g_x_0_yzzz_yyyyyzz[k];

                g_x_0_yyzzz_yyyzzz[k] = -g_x_0_yzzz_yyyzzz[k] * ab_y + g_x_0_yzzz_yyyyzzz[k];

                g_x_0_yyzzz_yyzzzz[k] = -g_x_0_yzzz_yyzzzz[k] * ab_y + g_x_0_yzzz_yyyzzzz[k];

                g_x_0_yyzzz_yzzzzz[k] = -g_x_0_yzzz_yzzzzz[k] * ab_y + g_x_0_yzzz_yyzzzzz[k];

                g_x_0_yyzzz_zzzzzz[k] = -g_x_0_yzzz_zzzzzz[k] * ab_y + g_x_0_yzzz_yzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 532 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 533 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 534 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 535 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 536 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 537 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 538 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 539 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 540 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 541 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 542 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 543 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 544 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 545 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 546 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 547 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 548 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 549 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 550 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 551 * ccomps * dcomps);

            auto g_x_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 552 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 553 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 554 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 555 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 556 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 557 * ccomps * dcomps);

            auto g_x_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 558 * ccomps * dcomps);

            auto g_x_0_yzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzzz_xxxxxx, g_x_0_yzzzz_xxxxxy, g_x_0_yzzzz_xxxxxz, g_x_0_yzzzz_xxxxyy, g_x_0_yzzzz_xxxxyz, g_x_0_yzzzz_xxxxzz, g_x_0_yzzzz_xxxyyy, g_x_0_yzzzz_xxxyyz, g_x_0_yzzzz_xxxyzz, g_x_0_yzzzz_xxxzzz, g_x_0_yzzzz_xxyyyy, g_x_0_yzzzz_xxyyyz, g_x_0_yzzzz_xxyyzz, g_x_0_yzzzz_xxyzzz, g_x_0_yzzzz_xxzzzz, g_x_0_yzzzz_xyyyyy, g_x_0_yzzzz_xyyyyz, g_x_0_yzzzz_xyyyzz, g_x_0_yzzzz_xyyzzz, g_x_0_yzzzz_xyzzzz, g_x_0_yzzzz_xzzzzz, g_x_0_yzzzz_yyyyyy, g_x_0_yzzzz_yyyyyz, g_x_0_yzzzz_yyyyzz, g_x_0_yzzzz_yyyzzz, g_x_0_yzzzz_yyzzzz, g_x_0_yzzzz_yzzzzz, g_x_0_yzzzz_zzzzzz, g_x_0_zzzz_xxxxxx, g_x_0_zzzz_xxxxxxy, g_x_0_zzzz_xxxxxy, g_x_0_zzzz_xxxxxyy, g_x_0_zzzz_xxxxxyz, g_x_0_zzzz_xxxxxz, g_x_0_zzzz_xxxxyy, g_x_0_zzzz_xxxxyyy, g_x_0_zzzz_xxxxyyz, g_x_0_zzzz_xxxxyz, g_x_0_zzzz_xxxxyzz, g_x_0_zzzz_xxxxzz, g_x_0_zzzz_xxxyyy, g_x_0_zzzz_xxxyyyy, g_x_0_zzzz_xxxyyyz, g_x_0_zzzz_xxxyyz, g_x_0_zzzz_xxxyyzz, g_x_0_zzzz_xxxyzz, g_x_0_zzzz_xxxyzzz, g_x_0_zzzz_xxxzzz, g_x_0_zzzz_xxyyyy, g_x_0_zzzz_xxyyyyy, g_x_0_zzzz_xxyyyyz, g_x_0_zzzz_xxyyyz, g_x_0_zzzz_xxyyyzz, g_x_0_zzzz_xxyyzz, g_x_0_zzzz_xxyyzzz, g_x_0_zzzz_xxyzzz, g_x_0_zzzz_xxyzzzz, g_x_0_zzzz_xxzzzz, g_x_0_zzzz_xyyyyy, g_x_0_zzzz_xyyyyyy, g_x_0_zzzz_xyyyyyz, g_x_0_zzzz_xyyyyz, g_x_0_zzzz_xyyyyzz, g_x_0_zzzz_xyyyzz, g_x_0_zzzz_xyyyzzz, g_x_0_zzzz_xyyzzz, g_x_0_zzzz_xyyzzzz, g_x_0_zzzz_xyzzzz, g_x_0_zzzz_xyzzzzz, g_x_0_zzzz_xzzzzz, g_x_0_zzzz_yyyyyy, g_x_0_zzzz_yyyyyyy, g_x_0_zzzz_yyyyyyz, g_x_0_zzzz_yyyyyz, g_x_0_zzzz_yyyyyzz, g_x_0_zzzz_yyyyzz, g_x_0_zzzz_yyyyzzz, g_x_0_zzzz_yyyzzz, g_x_0_zzzz_yyyzzzz, g_x_0_zzzz_yyzzzz, g_x_0_zzzz_yyzzzzz, g_x_0_zzzz_yzzzzz, g_x_0_zzzz_yzzzzzz, g_x_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzz_xxxxxx[k] = -g_x_0_zzzz_xxxxxx[k] * ab_y + g_x_0_zzzz_xxxxxxy[k];

                g_x_0_yzzzz_xxxxxy[k] = -g_x_0_zzzz_xxxxxy[k] * ab_y + g_x_0_zzzz_xxxxxyy[k];

                g_x_0_yzzzz_xxxxxz[k] = -g_x_0_zzzz_xxxxxz[k] * ab_y + g_x_0_zzzz_xxxxxyz[k];

                g_x_0_yzzzz_xxxxyy[k] = -g_x_0_zzzz_xxxxyy[k] * ab_y + g_x_0_zzzz_xxxxyyy[k];

                g_x_0_yzzzz_xxxxyz[k] = -g_x_0_zzzz_xxxxyz[k] * ab_y + g_x_0_zzzz_xxxxyyz[k];

                g_x_0_yzzzz_xxxxzz[k] = -g_x_0_zzzz_xxxxzz[k] * ab_y + g_x_0_zzzz_xxxxyzz[k];

                g_x_0_yzzzz_xxxyyy[k] = -g_x_0_zzzz_xxxyyy[k] * ab_y + g_x_0_zzzz_xxxyyyy[k];

                g_x_0_yzzzz_xxxyyz[k] = -g_x_0_zzzz_xxxyyz[k] * ab_y + g_x_0_zzzz_xxxyyyz[k];

                g_x_0_yzzzz_xxxyzz[k] = -g_x_0_zzzz_xxxyzz[k] * ab_y + g_x_0_zzzz_xxxyyzz[k];

                g_x_0_yzzzz_xxxzzz[k] = -g_x_0_zzzz_xxxzzz[k] * ab_y + g_x_0_zzzz_xxxyzzz[k];

                g_x_0_yzzzz_xxyyyy[k] = -g_x_0_zzzz_xxyyyy[k] * ab_y + g_x_0_zzzz_xxyyyyy[k];

                g_x_0_yzzzz_xxyyyz[k] = -g_x_0_zzzz_xxyyyz[k] * ab_y + g_x_0_zzzz_xxyyyyz[k];

                g_x_0_yzzzz_xxyyzz[k] = -g_x_0_zzzz_xxyyzz[k] * ab_y + g_x_0_zzzz_xxyyyzz[k];

                g_x_0_yzzzz_xxyzzz[k] = -g_x_0_zzzz_xxyzzz[k] * ab_y + g_x_0_zzzz_xxyyzzz[k];

                g_x_0_yzzzz_xxzzzz[k] = -g_x_0_zzzz_xxzzzz[k] * ab_y + g_x_0_zzzz_xxyzzzz[k];

                g_x_0_yzzzz_xyyyyy[k] = -g_x_0_zzzz_xyyyyy[k] * ab_y + g_x_0_zzzz_xyyyyyy[k];

                g_x_0_yzzzz_xyyyyz[k] = -g_x_0_zzzz_xyyyyz[k] * ab_y + g_x_0_zzzz_xyyyyyz[k];

                g_x_0_yzzzz_xyyyzz[k] = -g_x_0_zzzz_xyyyzz[k] * ab_y + g_x_0_zzzz_xyyyyzz[k];

                g_x_0_yzzzz_xyyzzz[k] = -g_x_0_zzzz_xyyzzz[k] * ab_y + g_x_0_zzzz_xyyyzzz[k];

                g_x_0_yzzzz_xyzzzz[k] = -g_x_0_zzzz_xyzzzz[k] * ab_y + g_x_0_zzzz_xyyzzzz[k];

                g_x_0_yzzzz_xzzzzz[k] = -g_x_0_zzzz_xzzzzz[k] * ab_y + g_x_0_zzzz_xyzzzzz[k];

                g_x_0_yzzzz_yyyyyy[k] = -g_x_0_zzzz_yyyyyy[k] * ab_y + g_x_0_zzzz_yyyyyyy[k];

                g_x_0_yzzzz_yyyyyz[k] = -g_x_0_zzzz_yyyyyz[k] * ab_y + g_x_0_zzzz_yyyyyyz[k];

                g_x_0_yzzzz_yyyyzz[k] = -g_x_0_zzzz_yyyyzz[k] * ab_y + g_x_0_zzzz_yyyyyzz[k];

                g_x_0_yzzzz_yyyzzz[k] = -g_x_0_zzzz_yyyzzz[k] * ab_y + g_x_0_zzzz_yyyyzzz[k];

                g_x_0_yzzzz_yyzzzz[k] = -g_x_0_zzzz_yyzzzz[k] * ab_y + g_x_0_zzzz_yyyzzzz[k];

                g_x_0_yzzzz_yzzzzz[k] = -g_x_0_zzzz_yzzzzz[k] * ab_y + g_x_0_zzzz_yyzzzzz[k];

                g_x_0_yzzzz_zzzzzz[k] = -g_x_0_zzzz_zzzzzz[k] * ab_y + g_x_0_zzzz_yzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 560 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 561 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 562 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 563 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 564 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 565 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 566 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 567 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 568 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 569 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 570 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 571 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 572 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 573 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 574 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 575 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 576 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 577 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 578 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 579 * ccomps * dcomps);

            auto g_x_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 580 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 581 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 582 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 583 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 584 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 585 * ccomps * dcomps);

            auto g_x_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 586 * ccomps * dcomps);

            auto g_x_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzz_xxxxxx, g_x_0_zzzz_xxxxxxz, g_x_0_zzzz_xxxxxy, g_x_0_zzzz_xxxxxyz, g_x_0_zzzz_xxxxxz, g_x_0_zzzz_xxxxxzz, g_x_0_zzzz_xxxxyy, g_x_0_zzzz_xxxxyyz, g_x_0_zzzz_xxxxyz, g_x_0_zzzz_xxxxyzz, g_x_0_zzzz_xxxxzz, g_x_0_zzzz_xxxxzzz, g_x_0_zzzz_xxxyyy, g_x_0_zzzz_xxxyyyz, g_x_0_zzzz_xxxyyz, g_x_0_zzzz_xxxyyzz, g_x_0_zzzz_xxxyzz, g_x_0_zzzz_xxxyzzz, g_x_0_zzzz_xxxzzz, g_x_0_zzzz_xxxzzzz, g_x_0_zzzz_xxyyyy, g_x_0_zzzz_xxyyyyz, g_x_0_zzzz_xxyyyz, g_x_0_zzzz_xxyyyzz, g_x_0_zzzz_xxyyzz, g_x_0_zzzz_xxyyzzz, g_x_0_zzzz_xxyzzz, g_x_0_zzzz_xxyzzzz, g_x_0_zzzz_xxzzzz, g_x_0_zzzz_xxzzzzz, g_x_0_zzzz_xyyyyy, g_x_0_zzzz_xyyyyyz, g_x_0_zzzz_xyyyyz, g_x_0_zzzz_xyyyyzz, g_x_0_zzzz_xyyyzz, g_x_0_zzzz_xyyyzzz, g_x_0_zzzz_xyyzzz, g_x_0_zzzz_xyyzzzz, g_x_0_zzzz_xyzzzz, g_x_0_zzzz_xyzzzzz, g_x_0_zzzz_xzzzzz, g_x_0_zzzz_xzzzzzz, g_x_0_zzzz_yyyyyy, g_x_0_zzzz_yyyyyyz, g_x_0_zzzz_yyyyyz, g_x_0_zzzz_yyyyyzz, g_x_0_zzzz_yyyyzz, g_x_0_zzzz_yyyyzzz, g_x_0_zzzz_yyyzzz, g_x_0_zzzz_yyyzzzz, g_x_0_zzzz_yyzzzz, g_x_0_zzzz_yyzzzzz, g_x_0_zzzz_yzzzzz, g_x_0_zzzz_yzzzzzz, g_x_0_zzzz_zzzzzz, g_x_0_zzzz_zzzzzzz, g_x_0_zzzzz_xxxxxx, g_x_0_zzzzz_xxxxxy, g_x_0_zzzzz_xxxxxz, g_x_0_zzzzz_xxxxyy, g_x_0_zzzzz_xxxxyz, g_x_0_zzzzz_xxxxzz, g_x_0_zzzzz_xxxyyy, g_x_0_zzzzz_xxxyyz, g_x_0_zzzzz_xxxyzz, g_x_0_zzzzz_xxxzzz, g_x_0_zzzzz_xxyyyy, g_x_0_zzzzz_xxyyyz, g_x_0_zzzzz_xxyyzz, g_x_0_zzzzz_xxyzzz, g_x_0_zzzzz_xxzzzz, g_x_0_zzzzz_xyyyyy, g_x_0_zzzzz_xyyyyz, g_x_0_zzzzz_xyyyzz, g_x_0_zzzzz_xyyzzz, g_x_0_zzzzz_xyzzzz, g_x_0_zzzzz_xzzzzz, g_x_0_zzzzz_yyyyyy, g_x_0_zzzzz_yyyyyz, g_x_0_zzzzz_yyyyzz, g_x_0_zzzzz_yyyzzz, g_x_0_zzzzz_yyzzzz, g_x_0_zzzzz_yzzzzz, g_x_0_zzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzz_xxxxxx[k] = -g_x_0_zzzz_xxxxxx[k] * ab_z + g_x_0_zzzz_xxxxxxz[k];

                g_x_0_zzzzz_xxxxxy[k] = -g_x_0_zzzz_xxxxxy[k] * ab_z + g_x_0_zzzz_xxxxxyz[k];

                g_x_0_zzzzz_xxxxxz[k] = -g_x_0_zzzz_xxxxxz[k] * ab_z + g_x_0_zzzz_xxxxxzz[k];

                g_x_0_zzzzz_xxxxyy[k] = -g_x_0_zzzz_xxxxyy[k] * ab_z + g_x_0_zzzz_xxxxyyz[k];

                g_x_0_zzzzz_xxxxyz[k] = -g_x_0_zzzz_xxxxyz[k] * ab_z + g_x_0_zzzz_xxxxyzz[k];

                g_x_0_zzzzz_xxxxzz[k] = -g_x_0_zzzz_xxxxzz[k] * ab_z + g_x_0_zzzz_xxxxzzz[k];

                g_x_0_zzzzz_xxxyyy[k] = -g_x_0_zzzz_xxxyyy[k] * ab_z + g_x_0_zzzz_xxxyyyz[k];

                g_x_0_zzzzz_xxxyyz[k] = -g_x_0_zzzz_xxxyyz[k] * ab_z + g_x_0_zzzz_xxxyyzz[k];

                g_x_0_zzzzz_xxxyzz[k] = -g_x_0_zzzz_xxxyzz[k] * ab_z + g_x_0_zzzz_xxxyzzz[k];

                g_x_0_zzzzz_xxxzzz[k] = -g_x_0_zzzz_xxxzzz[k] * ab_z + g_x_0_zzzz_xxxzzzz[k];

                g_x_0_zzzzz_xxyyyy[k] = -g_x_0_zzzz_xxyyyy[k] * ab_z + g_x_0_zzzz_xxyyyyz[k];

                g_x_0_zzzzz_xxyyyz[k] = -g_x_0_zzzz_xxyyyz[k] * ab_z + g_x_0_zzzz_xxyyyzz[k];

                g_x_0_zzzzz_xxyyzz[k] = -g_x_0_zzzz_xxyyzz[k] * ab_z + g_x_0_zzzz_xxyyzzz[k];

                g_x_0_zzzzz_xxyzzz[k] = -g_x_0_zzzz_xxyzzz[k] * ab_z + g_x_0_zzzz_xxyzzzz[k];

                g_x_0_zzzzz_xxzzzz[k] = -g_x_0_zzzz_xxzzzz[k] * ab_z + g_x_0_zzzz_xxzzzzz[k];

                g_x_0_zzzzz_xyyyyy[k] = -g_x_0_zzzz_xyyyyy[k] * ab_z + g_x_0_zzzz_xyyyyyz[k];

                g_x_0_zzzzz_xyyyyz[k] = -g_x_0_zzzz_xyyyyz[k] * ab_z + g_x_0_zzzz_xyyyyzz[k];

                g_x_0_zzzzz_xyyyzz[k] = -g_x_0_zzzz_xyyyzz[k] * ab_z + g_x_0_zzzz_xyyyzzz[k];

                g_x_0_zzzzz_xyyzzz[k] = -g_x_0_zzzz_xyyzzz[k] * ab_z + g_x_0_zzzz_xyyzzzz[k];

                g_x_0_zzzzz_xyzzzz[k] = -g_x_0_zzzz_xyzzzz[k] * ab_z + g_x_0_zzzz_xyzzzzz[k];

                g_x_0_zzzzz_xzzzzz[k] = -g_x_0_zzzz_xzzzzz[k] * ab_z + g_x_0_zzzz_xzzzzzz[k];

                g_x_0_zzzzz_yyyyyy[k] = -g_x_0_zzzz_yyyyyy[k] * ab_z + g_x_0_zzzz_yyyyyyz[k];

                g_x_0_zzzzz_yyyyyz[k] = -g_x_0_zzzz_yyyyyz[k] * ab_z + g_x_0_zzzz_yyyyyzz[k];

                g_x_0_zzzzz_yyyyzz[k] = -g_x_0_zzzz_yyyyzz[k] * ab_z + g_x_0_zzzz_yyyyzzz[k];

                g_x_0_zzzzz_yyyzzz[k] = -g_x_0_zzzz_yyyzzz[k] * ab_z + g_x_0_zzzz_yyyzzzz[k];

                g_x_0_zzzzz_yyzzzz[k] = -g_x_0_zzzz_yyzzzz[k] * ab_z + g_x_0_zzzz_yyzzzzz[k];

                g_x_0_zzzzz_yzzzzz[k] = -g_x_0_zzzz_yzzzzz[k] * ab_z + g_x_0_zzzz_yzzzzzz[k];

                g_x_0_zzzzz_zzzzzz[k] = -g_x_0_zzzz_zzzzzz[k] * ab_z + g_x_0_zzzz_zzzzzzz[k];
            }

            /// Set up 588-616 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 608 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyyyy = cbuffer.data(hi_geom_10_off + 609 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyyyz = cbuffer.data(hi_geom_10_off + 610 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyyzz = cbuffer.data(hi_geom_10_off + 611 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyzzz = cbuffer.data(hi_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyzzzz = cbuffer.data(hi_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_xxxxx_yzzzzz = cbuffer.data(hi_geom_10_off + 614 * ccomps * dcomps);

            auto g_y_0_xxxxx_zzzzzz = cbuffer.data(hi_geom_10_off + 615 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxx_xxxxxx, g_y_0_xxxx_xxxxxxx, g_y_0_xxxx_xxxxxxy, g_y_0_xxxx_xxxxxxz, g_y_0_xxxx_xxxxxy, g_y_0_xxxx_xxxxxyy, g_y_0_xxxx_xxxxxyz, g_y_0_xxxx_xxxxxz, g_y_0_xxxx_xxxxxzz, g_y_0_xxxx_xxxxyy, g_y_0_xxxx_xxxxyyy, g_y_0_xxxx_xxxxyyz, g_y_0_xxxx_xxxxyz, g_y_0_xxxx_xxxxyzz, g_y_0_xxxx_xxxxzz, g_y_0_xxxx_xxxxzzz, g_y_0_xxxx_xxxyyy, g_y_0_xxxx_xxxyyyy, g_y_0_xxxx_xxxyyyz, g_y_0_xxxx_xxxyyz, g_y_0_xxxx_xxxyyzz, g_y_0_xxxx_xxxyzz, g_y_0_xxxx_xxxyzzz, g_y_0_xxxx_xxxzzz, g_y_0_xxxx_xxxzzzz, g_y_0_xxxx_xxyyyy, g_y_0_xxxx_xxyyyyy, g_y_0_xxxx_xxyyyyz, g_y_0_xxxx_xxyyyz, g_y_0_xxxx_xxyyyzz, g_y_0_xxxx_xxyyzz, g_y_0_xxxx_xxyyzzz, g_y_0_xxxx_xxyzzz, g_y_0_xxxx_xxyzzzz, g_y_0_xxxx_xxzzzz, g_y_0_xxxx_xxzzzzz, g_y_0_xxxx_xyyyyy, g_y_0_xxxx_xyyyyyy, g_y_0_xxxx_xyyyyyz, g_y_0_xxxx_xyyyyz, g_y_0_xxxx_xyyyyzz, g_y_0_xxxx_xyyyzz, g_y_0_xxxx_xyyyzzz, g_y_0_xxxx_xyyzzz, g_y_0_xxxx_xyyzzzz, g_y_0_xxxx_xyzzzz, g_y_0_xxxx_xyzzzzz, g_y_0_xxxx_xzzzzz, g_y_0_xxxx_xzzzzzz, g_y_0_xxxx_yyyyyy, g_y_0_xxxx_yyyyyz, g_y_0_xxxx_yyyyzz, g_y_0_xxxx_yyyzzz, g_y_0_xxxx_yyzzzz, g_y_0_xxxx_yzzzzz, g_y_0_xxxx_zzzzzz, g_y_0_xxxxx_xxxxxx, g_y_0_xxxxx_xxxxxy, g_y_0_xxxxx_xxxxxz, g_y_0_xxxxx_xxxxyy, g_y_0_xxxxx_xxxxyz, g_y_0_xxxxx_xxxxzz, g_y_0_xxxxx_xxxyyy, g_y_0_xxxxx_xxxyyz, g_y_0_xxxxx_xxxyzz, g_y_0_xxxxx_xxxzzz, g_y_0_xxxxx_xxyyyy, g_y_0_xxxxx_xxyyyz, g_y_0_xxxxx_xxyyzz, g_y_0_xxxxx_xxyzzz, g_y_0_xxxxx_xxzzzz, g_y_0_xxxxx_xyyyyy, g_y_0_xxxxx_xyyyyz, g_y_0_xxxxx_xyyyzz, g_y_0_xxxxx_xyyzzz, g_y_0_xxxxx_xyzzzz, g_y_0_xxxxx_xzzzzz, g_y_0_xxxxx_yyyyyy, g_y_0_xxxxx_yyyyyz, g_y_0_xxxxx_yyyyzz, g_y_0_xxxxx_yyyzzz, g_y_0_xxxxx_yyzzzz, g_y_0_xxxxx_yzzzzz, g_y_0_xxxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxx_xxxxxx[k] = -g_y_0_xxxx_xxxxxx[k] * ab_x + g_y_0_xxxx_xxxxxxx[k];

                g_y_0_xxxxx_xxxxxy[k] = -g_y_0_xxxx_xxxxxy[k] * ab_x + g_y_0_xxxx_xxxxxxy[k];

                g_y_0_xxxxx_xxxxxz[k] = -g_y_0_xxxx_xxxxxz[k] * ab_x + g_y_0_xxxx_xxxxxxz[k];

                g_y_0_xxxxx_xxxxyy[k] = -g_y_0_xxxx_xxxxyy[k] * ab_x + g_y_0_xxxx_xxxxxyy[k];

                g_y_0_xxxxx_xxxxyz[k] = -g_y_0_xxxx_xxxxyz[k] * ab_x + g_y_0_xxxx_xxxxxyz[k];

                g_y_0_xxxxx_xxxxzz[k] = -g_y_0_xxxx_xxxxzz[k] * ab_x + g_y_0_xxxx_xxxxxzz[k];

                g_y_0_xxxxx_xxxyyy[k] = -g_y_0_xxxx_xxxyyy[k] * ab_x + g_y_0_xxxx_xxxxyyy[k];

                g_y_0_xxxxx_xxxyyz[k] = -g_y_0_xxxx_xxxyyz[k] * ab_x + g_y_0_xxxx_xxxxyyz[k];

                g_y_0_xxxxx_xxxyzz[k] = -g_y_0_xxxx_xxxyzz[k] * ab_x + g_y_0_xxxx_xxxxyzz[k];

                g_y_0_xxxxx_xxxzzz[k] = -g_y_0_xxxx_xxxzzz[k] * ab_x + g_y_0_xxxx_xxxxzzz[k];

                g_y_0_xxxxx_xxyyyy[k] = -g_y_0_xxxx_xxyyyy[k] * ab_x + g_y_0_xxxx_xxxyyyy[k];

                g_y_0_xxxxx_xxyyyz[k] = -g_y_0_xxxx_xxyyyz[k] * ab_x + g_y_0_xxxx_xxxyyyz[k];

                g_y_0_xxxxx_xxyyzz[k] = -g_y_0_xxxx_xxyyzz[k] * ab_x + g_y_0_xxxx_xxxyyzz[k];

                g_y_0_xxxxx_xxyzzz[k] = -g_y_0_xxxx_xxyzzz[k] * ab_x + g_y_0_xxxx_xxxyzzz[k];

                g_y_0_xxxxx_xxzzzz[k] = -g_y_0_xxxx_xxzzzz[k] * ab_x + g_y_0_xxxx_xxxzzzz[k];

                g_y_0_xxxxx_xyyyyy[k] = -g_y_0_xxxx_xyyyyy[k] * ab_x + g_y_0_xxxx_xxyyyyy[k];

                g_y_0_xxxxx_xyyyyz[k] = -g_y_0_xxxx_xyyyyz[k] * ab_x + g_y_0_xxxx_xxyyyyz[k];

                g_y_0_xxxxx_xyyyzz[k] = -g_y_0_xxxx_xyyyzz[k] * ab_x + g_y_0_xxxx_xxyyyzz[k];

                g_y_0_xxxxx_xyyzzz[k] = -g_y_0_xxxx_xyyzzz[k] * ab_x + g_y_0_xxxx_xxyyzzz[k];

                g_y_0_xxxxx_xyzzzz[k] = -g_y_0_xxxx_xyzzzz[k] * ab_x + g_y_0_xxxx_xxyzzzz[k];

                g_y_0_xxxxx_xzzzzz[k] = -g_y_0_xxxx_xzzzzz[k] * ab_x + g_y_0_xxxx_xxzzzzz[k];

                g_y_0_xxxxx_yyyyyy[k] = -g_y_0_xxxx_yyyyyy[k] * ab_x + g_y_0_xxxx_xyyyyyy[k];

                g_y_0_xxxxx_yyyyyz[k] = -g_y_0_xxxx_yyyyyz[k] * ab_x + g_y_0_xxxx_xyyyyyz[k];

                g_y_0_xxxxx_yyyyzz[k] = -g_y_0_xxxx_yyyyzz[k] * ab_x + g_y_0_xxxx_xyyyyzz[k];

                g_y_0_xxxxx_yyyzzz[k] = -g_y_0_xxxx_yyyzzz[k] * ab_x + g_y_0_xxxx_xyyyzzz[k];

                g_y_0_xxxxx_yyzzzz[k] = -g_y_0_xxxx_yyzzzz[k] * ab_x + g_y_0_xxxx_xyyzzzz[k];

                g_y_0_xxxxx_yzzzzz[k] = -g_y_0_xxxx_yzzzzz[k] * ab_x + g_y_0_xxxx_xyzzzzz[k];

                g_y_0_xxxxx_zzzzzz[k] = -g_y_0_xxxx_zzzzzz[k] * ab_x + g_y_0_xxxx_xzzzzzz[k];
            }

            /// Set up 616-644 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 629 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 630 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 631 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 632 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 633 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 634 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 635 * ccomps * dcomps);

            auto g_y_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 636 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyyyy = cbuffer.data(hi_geom_10_off + 637 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyyyz = cbuffer.data(hi_geom_10_off + 638 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyyzz = cbuffer.data(hi_geom_10_off + 639 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyzzz = cbuffer.data(hi_geom_10_off + 640 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyzzzz = cbuffer.data(hi_geom_10_off + 641 * ccomps * dcomps);

            auto g_y_0_xxxxy_yzzzzz = cbuffer.data(hi_geom_10_off + 642 * ccomps * dcomps);

            auto g_y_0_xxxxy_zzzzzz = cbuffer.data(hi_geom_10_off + 643 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxy_xxxxxx, g_y_0_xxxxy_xxxxxy, g_y_0_xxxxy_xxxxxz, g_y_0_xxxxy_xxxxyy, g_y_0_xxxxy_xxxxyz, g_y_0_xxxxy_xxxxzz, g_y_0_xxxxy_xxxyyy, g_y_0_xxxxy_xxxyyz, g_y_0_xxxxy_xxxyzz, g_y_0_xxxxy_xxxzzz, g_y_0_xxxxy_xxyyyy, g_y_0_xxxxy_xxyyyz, g_y_0_xxxxy_xxyyzz, g_y_0_xxxxy_xxyzzz, g_y_0_xxxxy_xxzzzz, g_y_0_xxxxy_xyyyyy, g_y_0_xxxxy_xyyyyz, g_y_0_xxxxy_xyyyzz, g_y_0_xxxxy_xyyzzz, g_y_0_xxxxy_xyzzzz, g_y_0_xxxxy_xzzzzz, g_y_0_xxxxy_yyyyyy, g_y_0_xxxxy_yyyyyz, g_y_0_xxxxy_yyyyzz, g_y_0_xxxxy_yyyzzz, g_y_0_xxxxy_yyzzzz, g_y_0_xxxxy_yzzzzz, g_y_0_xxxxy_zzzzzz, g_y_0_xxxy_xxxxxx, g_y_0_xxxy_xxxxxxx, g_y_0_xxxy_xxxxxxy, g_y_0_xxxy_xxxxxxz, g_y_0_xxxy_xxxxxy, g_y_0_xxxy_xxxxxyy, g_y_0_xxxy_xxxxxyz, g_y_0_xxxy_xxxxxz, g_y_0_xxxy_xxxxxzz, g_y_0_xxxy_xxxxyy, g_y_0_xxxy_xxxxyyy, g_y_0_xxxy_xxxxyyz, g_y_0_xxxy_xxxxyz, g_y_0_xxxy_xxxxyzz, g_y_0_xxxy_xxxxzz, g_y_0_xxxy_xxxxzzz, g_y_0_xxxy_xxxyyy, g_y_0_xxxy_xxxyyyy, g_y_0_xxxy_xxxyyyz, g_y_0_xxxy_xxxyyz, g_y_0_xxxy_xxxyyzz, g_y_0_xxxy_xxxyzz, g_y_0_xxxy_xxxyzzz, g_y_0_xxxy_xxxzzz, g_y_0_xxxy_xxxzzzz, g_y_0_xxxy_xxyyyy, g_y_0_xxxy_xxyyyyy, g_y_0_xxxy_xxyyyyz, g_y_0_xxxy_xxyyyz, g_y_0_xxxy_xxyyyzz, g_y_0_xxxy_xxyyzz, g_y_0_xxxy_xxyyzzz, g_y_0_xxxy_xxyzzz, g_y_0_xxxy_xxyzzzz, g_y_0_xxxy_xxzzzz, g_y_0_xxxy_xxzzzzz, g_y_0_xxxy_xyyyyy, g_y_0_xxxy_xyyyyyy, g_y_0_xxxy_xyyyyyz, g_y_0_xxxy_xyyyyz, g_y_0_xxxy_xyyyyzz, g_y_0_xxxy_xyyyzz, g_y_0_xxxy_xyyyzzz, g_y_0_xxxy_xyyzzz, g_y_0_xxxy_xyyzzzz, g_y_0_xxxy_xyzzzz, g_y_0_xxxy_xyzzzzz, g_y_0_xxxy_xzzzzz, g_y_0_xxxy_xzzzzzz, g_y_0_xxxy_yyyyyy, g_y_0_xxxy_yyyyyz, g_y_0_xxxy_yyyyzz, g_y_0_xxxy_yyyzzz, g_y_0_xxxy_yyzzzz, g_y_0_xxxy_yzzzzz, g_y_0_xxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxy_xxxxxx[k] = -g_y_0_xxxy_xxxxxx[k] * ab_x + g_y_0_xxxy_xxxxxxx[k];

                g_y_0_xxxxy_xxxxxy[k] = -g_y_0_xxxy_xxxxxy[k] * ab_x + g_y_0_xxxy_xxxxxxy[k];

                g_y_0_xxxxy_xxxxxz[k] = -g_y_0_xxxy_xxxxxz[k] * ab_x + g_y_0_xxxy_xxxxxxz[k];

                g_y_0_xxxxy_xxxxyy[k] = -g_y_0_xxxy_xxxxyy[k] * ab_x + g_y_0_xxxy_xxxxxyy[k];

                g_y_0_xxxxy_xxxxyz[k] = -g_y_0_xxxy_xxxxyz[k] * ab_x + g_y_0_xxxy_xxxxxyz[k];

                g_y_0_xxxxy_xxxxzz[k] = -g_y_0_xxxy_xxxxzz[k] * ab_x + g_y_0_xxxy_xxxxxzz[k];

                g_y_0_xxxxy_xxxyyy[k] = -g_y_0_xxxy_xxxyyy[k] * ab_x + g_y_0_xxxy_xxxxyyy[k];

                g_y_0_xxxxy_xxxyyz[k] = -g_y_0_xxxy_xxxyyz[k] * ab_x + g_y_0_xxxy_xxxxyyz[k];

                g_y_0_xxxxy_xxxyzz[k] = -g_y_0_xxxy_xxxyzz[k] * ab_x + g_y_0_xxxy_xxxxyzz[k];

                g_y_0_xxxxy_xxxzzz[k] = -g_y_0_xxxy_xxxzzz[k] * ab_x + g_y_0_xxxy_xxxxzzz[k];

                g_y_0_xxxxy_xxyyyy[k] = -g_y_0_xxxy_xxyyyy[k] * ab_x + g_y_0_xxxy_xxxyyyy[k];

                g_y_0_xxxxy_xxyyyz[k] = -g_y_0_xxxy_xxyyyz[k] * ab_x + g_y_0_xxxy_xxxyyyz[k];

                g_y_0_xxxxy_xxyyzz[k] = -g_y_0_xxxy_xxyyzz[k] * ab_x + g_y_0_xxxy_xxxyyzz[k];

                g_y_0_xxxxy_xxyzzz[k] = -g_y_0_xxxy_xxyzzz[k] * ab_x + g_y_0_xxxy_xxxyzzz[k];

                g_y_0_xxxxy_xxzzzz[k] = -g_y_0_xxxy_xxzzzz[k] * ab_x + g_y_0_xxxy_xxxzzzz[k];

                g_y_0_xxxxy_xyyyyy[k] = -g_y_0_xxxy_xyyyyy[k] * ab_x + g_y_0_xxxy_xxyyyyy[k];

                g_y_0_xxxxy_xyyyyz[k] = -g_y_0_xxxy_xyyyyz[k] * ab_x + g_y_0_xxxy_xxyyyyz[k];

                g_y_0_xxxxy_xyyyzz[k] = -g_y_0_xxxy_xyyyzz[k] * ab_x + g_y_0_xxxy_xxyyyzz[k];

                g_y_0_xxxxy_xyyzzz[k] = -g_y_0_xxxy_xyyzzz[k] * ab_x + g_y_0_xxxy_xxyyzzz[k];

                g_y_0_xxxxy_xyzzzz[k] = -g_y_0_xxxy_xyzzzz[k] * ab_x + g_y_0_xxxy_xxyzzzz[k];

                g_y_0_xxxxy_xzzzzz[k] = -g_y_0_xxxy_xzzzzz[k] * ab_x + g_y_0_xxxy_xxzzzzz[k];

                g_y_0_xxxxy_yyyyyy[k] = -g_y_0_xxxy_yyyyyy[k] * ab_x + g_y_0_xxxy_xyyyyyy[k];

                g_y_0_xxxxy_yyyyyz[k] = -g_y_0_xxxy_yyyyyz[k] * ab_x + g_y_0_xxxy_xyyyyyz[k];

                g_y_0_xxxxy_yyyyzz[k] = -g_y_0_xxxy_yyyyzz[k] * ab_x + g_y_0_xxxy_xyyyyzz[k];

                g_y_0_xxxxy_yyyzzz[k] = -g_y_0_xxxy_yyyzzz[k] * ab_x + g_y_0_xxxy_xyyyzzz[k];

                g_y_0_xxxxy_yyzzzz[k] = -g_y_0_xxxy_yyzzzz[k] * ab_x + g_y_0_xxxy_xyyzzzz[k];

                g_y_0_xxxxy_yzzzzz[k] = -g_y_0_xxxy_yzzzzz[k] * ab_x + g_y_0_xxxy_xyzzzzz[k];

                g_y_0_xxxxy_zzzzzz[k] = -g_y_0_xxxy_zzzzzz[k] * ab_x + g_y_0_xxxy_xzzzzzz[k];
            }

            /// Set up 644-672 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 644 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 645 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 646 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 647 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 648 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 649 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 650 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 651 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 652 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 653 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 654 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 655 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 656 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 657 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 658 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 659 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 660 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 661 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 662 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 663 * ccomps * dcomps);

            auto g_y_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 664 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyyyy = cbuffer.data(hi_geom_10_off + 665 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyyyz = cbuffer.data(hi_geom_10_off + 666 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyyzz = cbuffer.data(hi_geom_10_off + 667 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyzzz = cbuffer.data(hi_geom_10_off + 668 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyzzzz = cbuffer.data(hi_geom_10_off + 669 * ccomps * dcomps);

            auto g_y_0_xxxxz_yzzzzz = cbuffer.data(hi_geom_10_off + 670 * ccomps * dcomps);

            auto g_y_0_xxxxz_zzzzzz = cbuffer.data(hi_geom_10_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxz_xxxxxx, g_y_0_xxxxz_xxxxxy, g_y_0_xxxxz_xxxxxz, g_y_0_xxxxz_xxxxyy, g_y_0_xxxxz_xxxxyz, g_y_0_xxxxz_xxxxzz, g_y_0_xxxxz_xxxyyy, g_y_0_xxxxz_xxxyyz, g_y_0_xxxxz_xxxyzz, g_y_0_xxxxz_xxxzzz, g_y_0_xxxxz_xxyyyy, g_y_0_xxxxz_xxyyyz, g_y_0_xxxxz_xxyyzz, g_y_0_xxxxz_xxyzzz, g_y_0_xxxxz_xxzzzz, g_y_0_xxxxz_xyyyyy, g_y_0_xxxxz_xyyyyz, g_y_0_xxxxz_xyyyzz, g_y_0_xxxxz_xyyzzz, g_y_0_xxxxz_xyzzzz, g_y_0_xxxxz_xzzzzz, g_y_0_xxxxz_yyyyyy, g_y_0_xxxxz_yyyyyz, g_y_0_xxxxz_yyyyzz, g_y_0_xxxxz_yyyzzz, g_y_0_xxxxz_yyzzzz, g_y_0_xxxxz_yzzzzz, g_y_0_xxxxz_zzzzzz, g_y_0_xxxz_xxxxxx, g_y_0_xxxz_xxxxxxx, g_y_0_xxxz_xxxxxxy, g_y_0_xxxz_xxxxxxz, g_y_0_xxxz_xxxxxy, g_y_0_xxxz_xxxxxyy, g_y_0_xxxz_xxxxxyz, g_y_0_xxxz_xxxxxz, g_y_0_xxxz_xxxxxzz, g_y_0_xxxz_xxxxyy, g_y_0_xxxz_xxxxyyy, g_y_0_xxxz_xxxxyyz, g_y_0_xxxz_xxxxyz, g_y_0_xxxz_xxxxyzz, g_y_0_xxxz_xxxxzz, g_y_0_xxxz_xxxxzzz, g_y_0_xxxz_xxxyyy, g_y_0_xxxz_xxxyyyy, g_y_0_xxxz_xxxyyyz, g_y_0_xxxz_xxxyyz, g_y_0_xxxz_xxxyyzz, g_y_0_xxxz_xxxyzz, g_y_0_xxxz_xxxyzzz, g_y_0_xxxz_xxxzzz, g_y_0_xxxz_xxxzzzz, g_y_0_xxxz_xxyyyy, g_y_0_xxxz_xxyyyyy, g_y_0_xxxz_xxyyyyz, g_y_0_xxxz_xxyyyz, g_y_0_xxxz_xxyyyzz, g_y_0_xxxz_xxyyzz, g_y_0_xxxz_xxyyzzz, g_y_0_xxxz_xxyzzz, g_y_0_xxxz_xxyzzzz, g_y_0_xxxz_xxzzzz, g_y_0_xxxz_xxzzzzz, g_y_0_xxxz_xyyyyy, g_y_0_xxxz_xyyyyyy, g_y_0_xxxz_xyyyyyz, g_y_0_xxxz_xyyyyz, g_y_0_xxxz_xyyyyzz, g_y_0_xxxz_xyyyzz, g_y_0_xxxz_xyyyzzz, g_y_0_xxxz_xyyzzz, g_y_0_xxxz_xyyzzzz, g_y_0_xxxz_xyzzzz, g_y_0_xxxz_xyzzzzz, g_y_0_xxxz_xzzzzz, g_y_0_xxxz_xzzzzzz, g_y_0_xxxz_yyyyyy, g_y_0_xxxz_yyyyyz, g_y_0_xxxz_yyyyzz, g_y_0_xxxz_yyyzzz, g_y_0_xxxz_yyzzzz, g_y_0_xxxz_yzzzzz, g_y_0_xxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxz_xxxxxx[k] = -g_y_0_xxxz_xxxxxx[k] * ab_x + g_y_0_xxxz_xxxxxxx[k];

                g_y_0_xxxxz_xxxxxy[k] = -g_y_0_xxxz_xxxxxy[k] * ab_x + g_y_0_xxxz_xxxxxxy[k];

                g_y_0_xxxxz_xxxxxz[k] = -g_y_0_xxxz_xxxxxz[k] * ab_x + g_y_0_xxxz_xxxxxxz[k];

                g_y_0_xxxxz_xxxxyy[k] = -g_y_0_xxxz_xxxxyy[k] * ab_x + g_y_0_xxxz_xxxxxyy[k];

                g_y_0_xxxxz_xxxxyz[k] = -g_y_0_xxxz_xxxxyz[k] * ab_x + g_y_0_xxxz_xxxxxyz[k];

                g_y_0_xxxxz_xxxxzz[k] = -g_y_0_xxxz_xxxxzz[k] * ab_x + g_y_0_xxxz_xxxxxzz[k];

                g_y_0_xxxxz_xxxyyy[k] = -g_y_0_xxxz_xxxyyy[k] * ab_x + g_y_0_xxxz_xxxxyyy[k];

                g_y_0_xxxxz_xxxyyz[k] = -g_y_0_xxxz_xxxyyz[k] * ab_x + g_y_0_xxxz_xxxxyyz[k];

                g_y_0_xxxxz_xxxyzz[k] = -g_y_0_xxxz_xxxyzz[k] * ab_x + g_y_0_xxxz_xxxxyzz[k];

                g_y_0_xxxxz_xxxzzz[k] = -g_y_0_xxxz_xxxzzz[k] * ab_x + g_y_0_xxxz_xxxxzzz[k];

                g_y_0_xxxxz_xxyyyy[k] = -g_y_0_xxxz_xxyyyy[k] * ab_x + g_y_0_xxxz_xxxyyyy[k];

                g_y_0_xxxxz_xxyyyz[k] = -g_y_0_xxxz_xxyyyz[k] * ab_x + g_y_0_xxxz_xxxyyyz[k];

                g_y_0_xxxxz_xxyyzz[k] = -g_y_0_xxxz_xxyyzz[k] * ab_x + g_y_0_xxxz_xxxyyzz[k];

                g_y_0_xxxxz_xxyzzz[k] = -g_y_0_xxxz_xxyzzz[k] * ab_x + g_y_0_xxxz_xxxyzzz[k];

                g_y_0_xxxxz_xxzzzz[k] = -g_y_0_xxxz_xxzzzz[k] * ab_x + g_y_0_xxxz_xxxzzzz[k];

                g_y_0_xxxxz_xyyyyy[k] = -g_y_0_xxxz_xyyyyy[k] * ab_x + g_y_0_xxxz_xxyyyyy[k];

                g_y_0_xxxxz_xyyyyz[k] = -g_y_0_xxxz_xyyyyz[k] * ab_x + g_y_0_xxxz_xxyyyyz[k];

                g_y_0_xxxxz_xyyyzz[k] = -g_y_0_xxxz_xyyyzz[k] * ab_x + g_y_0_xxxz_xxyyyzz[k];

                g_y_0_xxxxz_xyyzzz[k] = -g_y_0_xxxz_xyyzzz[k] * ab_x + g_y_0_xxxz_xxyyzzz[k];

                g_y_0_xxxxz_xyzzzz[k] = -g_y_0_xxxz_xyzzzz[k] * ab_x + g_y_0_xxxz_xxyzzzz[k];

                g_y_0_xxxxz_xzzzzz[k] = -g_y_0_xxxz_xzzzzz[k] * ab_x + g_y_0_xxxz_xxzzzzz[k];

                g_y_0_xxxxz_yyyyyy[k] = -g_y_0_xxxz_yyyyyy[k] * ab_x + g_y_0_xxxz_xyyyyyy[k];

                g_y_0_xxxxz_yyyyyz[k] = -g_y_0_xxxz_yyyyyz[k] * ab_x + g_y_0_xxxz_xyyyyyz[k];

                g_y_0_xxxxz_yyyyzz[k] = -g_y_0_xxxz_yyyyzz[k] * ab_x + g_y_0_xxxz_xyyyyzz[k];

                g_y_0_xxxxz_yyyzzz[k] = -g_y_0_xxxz_yyyzzz[k] * ab_x + g_y_0_xxxz_xyyyzzz[k];

                g_y_0_xxxxz_yyzzzz[k] = -g_y_0_xxxz_yyzzzz[k] * ab_x + g_y_0_xxxz_xyyzzzz[k];

                g_y_0_xxxxz_yzzzzz[k] = -g_y_0_xxxz_yzzzzz[k] * ab_x + g_y_0_xxxz_xyzzzzz[k];

                g_y_0_xxxxz_zzzzzz[k] = -g_y_0_xxxz_zzzzzz[k] * ab_x + g_y_0_xxxz_xzzzzzz[k];
            }

            /// Set up 672-700 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 672 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 673 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 674 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 675 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 676 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 677 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 678 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 679 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 680 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 681 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 682 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 683 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 684 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 685 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 686 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 687 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 688 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 689 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 690 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 691 * ccomps * dcomps);

            auto g_y_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 692 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyyyy = cbuffer.data(hi_geom_10_off + 693 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyyyz = cbuffer.data(hi_geom_10_off + 694 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyyzz = cbuffer.data(hi_geom_10_off + 695 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyzzz = cbuffer.data(hi_geom_10_off + 696 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyzzzz = cbuffer.data(hi_geom_10_off + 697 * ccomps * dcomps);

            auto g_y_0_xxxyy_yzzzzz = cbuffer.data(hi_geom_10_off + 698 * ccomps * dcomps);

            auto g_y_0_xxxyy_zzzzzz = cbuffer.data(hi_geom_10_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyy_xxxxxx, g_y_0_xxxyy_xxxxxy, g_y_0_xxxyy_xxxxxz, g_y_0_xxxyy_xxxxyy, g_y_0_xxxyy_xxxxyz, g_y_0_xxxyy_xxxxzz, g_y_0_xxxyy_xxxyyy, g_y_0_xxxyy_xxxyyz, g_y_0_xxxyy_xxxyzz, g_y_0_xxxyy_xxxzzz, g_y_0_xxxyy_xxyyyy, g_y_0_xxxyy_xxyyyz, g_y_0_xxxyy_xxyyzz, g_y_0_xxxyy_xxyzzz, g_y_0_xxxyy_xxzzzz, g_y_0_xxxyy_xyyyyy, g_y_0_xxxyy_xyyyyz, g_y_0_xxxyy_xyyyzz, g_y_0_xxxyy_xyyzzz, g_y_0_xxxyy_xyzzzz, g_y_0_xxxyy_xzzzzz, g_y_0_xxxyy_yyyyyy, g_y_0_xxxyy_yyyyyz, g_y_0_xxxyy_yyyyzz, g_y_0_xxxyy_yyyzzz, g_y_0_xxxyy_yyzzzz, g_y_0_xxxyy_yzzzzz, g_y_0_xxxyy_zzzzzz, g_y_0_xxyy_xxxxxx, g_y_0_xxyy_xxxxxxx, g_y_0_xxyy_xxxxxxy, g_y_0_xxyy_xxxxxxz, g_y_0_xxyy_xxxxxy, g_y_0_xxyy_xxxxxyy, g_y_0_xxyy_xxxxxyz, g_y_0_xxyy_xxxxxz, g_y_0_xxyy_xxxxxzz, g_y_0_xxyy_xxxxyy, g_y_0_xxyy_xxxxyyy, g_y_0_xxyy_xxxxyyz, g_y_0_xxyy_xxxxyz, g_y_0_xxyy_xxxxyzz, g_y_0_xxyy_xxxxzz, g_y_0_xxyy_xxxxzzz, g_y_0_xxyy_xxxyyy, g_y_0_xxyy_xxxyyyy, g_y_0_xxyy_xxxyyyz, g_y_0_xxyy_xxxyyz, g_y_0_xxyy_xxxyyzz, g_y_0_xxyy_xxxyzz, g_y_0_xxyy_xxxyzzz, g_y_0_xxyy_xxxzzz, g_y_0_xxyy_xxxzzzz, g_y_0_xxyy_xxyyyy, g_y_0_xxyy_xxyyyyy, g_y_0_xxyy_xxyyyyz, g_y_0_xxyy_xxyyyz, g_y_0_xxyy_xxyyyzz, g_y_0_xxyy_xxyyzz, g_y_0_xxyy_xxyyzzz, g_y_0_xxyy_xxyzzz, g_y_0_xxyy_xxyzzzz, g_y_0_xxyy_xxzzzz, g_y_0_xxyy_xxzzzzz, g_y_0_xxyy_xyyyyy, g_y_0_xxyy_xyyyyyy, g_y_0_xxyy_xyyyyyz, g_y_0_xxyy_xyyyyz, g_y_0_xxyy_xyyyyzz, g_y_0_xxyy_xyyyzz, g_y_0_xxyy_xyyyzzz, g_y_0_xxyy_xyyzzz, g_y_0_xxyy_xyyzzzz, g_y_0_xxyy_xyzzzz, g_y_0_xxyy_xyzzzzz, g_y_0_xxyy_xzzzzz, g_y_0_xxyy_xzzzzzz, g_y_0_xxyy_yyyyyy, g_y_0_xxyy_yyyyyz, g_y_0_xxyy_yyyyzz, g_y_0_xxyy_yyyzzz, g_y_0_xxyy_yyzzzz, g_y_0_xxyy_yzzzzz, g_y_0_xxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyy_xxxxxx[k] = -g_y_0_xxyy_xxxxxx[k] * ab_x + g_y_0_xxyy_xxxxxxx[k];

                g_y_0_xxxyy_xxxxxy[k] = -g_y_0_xxyy_xxxxxy[k] * ab_x + g_y_0_xxyy_xxxxxxy[k];

                g_y_0_xxxyy_xxxxxz[k] = -g_y_0_xxyy_xxxxxz[k] * ab_x + g_y_0_xxyy_xxxxxxz[k];

                g_y_0_xxxyy_xxxxyy[k] = -g_y_0_xxyy_xxxxyy[k] * ab_x + g_y_0_xxyy_xxxxxyy[k];

                g_y_0_xxxyy_xxxxyz[k] = -g_y_0_xxyy_xxxxyz[k] * ab_x + g_y_0_xxyy_xxxxxyz[k];

                g_y_0_xxxyy_xxxxzz[k] = -g_y_0_xxyy_xxxxzz[k] * ab_x + g_y_0_xxyy_xxxxxzz[k];

                g_y_0_xxxyy_xxxyyy[k] = -g_y_0_xxyy_xxxyyy[k] * ab_x + g_y_0_xxyy_xxxxyyy[k];

                g_y_0_xxxyy_xxxyyz[k] = -g_y_0_xxyy_xxxyyz[k] * ab_x + g_y_0_xxyy_xxxxyyz[k];

                g_y_0_xxxyy_xxxyzz[k] = -g_y_0_xxyy_xxxyzz[k] * ab_x + g_y_0_xxyy_xxxxyzz[k];

                g_y_0_xxxyy_xxxzzz[k] = -g_y_0_xxyy_xxxzzz[k] * ab_x + g_y_0_xxyy_xxxxzzz[k];

                g_y_0_xxxyy_xxyyyy[k] = -g_y_0_xxyy_xxyyyy[k] * ab_x + g_y_0_xxyy_xxxyyyy[k];

                g_y_0_xxxyy_xxyyyz[k] = -g_y_0_xxyy_xxyyyz[k] * ab_x + g_y_0_xxyy_xxxyyyz[k];

                g_y_0_xxxyy_xxyyzz[k] = -g_y_0_xxyy_xxyyzz[k] * ab_x + g_y_0_xxyy_xxxyyzz[k];

                g_y_0_xxxyy_xxyzzz[k] = -g_y_0_xxyy_xxyzzz[k] * ab_x + g_y_0_xxyy_xxxyzzz[k];

                g_y_0_xxxyy_xxzzzz[k] = -g_y_0_xxyy_xxzzzz[k] * ab_x + g_y_0_xxyy_xxxzzzz[k];

                g_y_0_xxxyy_xyyyyy[k] = -g_y_0_xxyy_xyyyyy[k] * ab_x + g_y_0_xxyy_xxyyyyy[k];

                g_y_0_xxxyy_xyyyyz[k] = -g_y_0_xxyy_xyyyyz[k] * ab_x + g_y_0_xxyy_xxyyyyz[k];

                g_y_0_xxxyy_xyyyzz[k] = -g_y_0_xxyy_xyyyzz[k] * ab_x + g_y_0_xxyy_xxyyyzz[k];

                g_y_0_xxxyy_xyyzzz[k] = -g_y_0_xxyy_xyyzzz[k] * ab_x + g_y_0_xxyy_xxyyzzz[k];

                g_y_0_xxxyy_xyzzzz[k] = -g_y_0_xxyy_xyzzzz[k] * ab_x + g_y_0_xxyy_xxyzzzz[k];

                g_y_0_xxxyy_xzzzzz[k] = -g_y_0_xxyy_xzzzzz[k] * ab_x + g_y_0_xxyy_xxzzzzz[k];

                g_y_0_xxxyy_yyyyyy[k] = -g_y_0_xxyy_yyyyyy[k] * ab_x + g_y_0_xxyy_xyyyyyy[k];

                g_y_0_xxxyy_yyyyyz[k] = -g_y_0_xxyy_yyyyyz[k] * ab_x + g_y_0_xxyy_xyyyyyz[k];

                g_y_0_xxxyy_yyyyzz[k] = -g_y_0_xxyy_yyyyzz[k] * ab_x + g_y_0_xxyy_xyyyyzz[k];

                g_y_0_xxxyy_yyyzzz[k] = -g_y_0_xxyy_yyyzzz[k] * ab_x + g_y_0_xxyy_xyyyzzz[k];

                g_y_0_xxxyy_yyzzzz[k] = -g_y_0_xxyy_yyzzzz[k] * ab_x + g_y_0_xxyy_xyyzzzz[k];

                g_y_0_xxxyy_yzzzzz[k] = -g_y_0_xxyy_yzzzzz[k] * ab_x + g_y_0_xxyy_xyzzzzz[k];

                g_y_0_xxxyy_zzzzzz[k] = -g_y_0_xxyy_zzzzzz[k] * ab_x + g_y_0_xxyy_xzzzzzz[k];
            }

            /// Set up 700-728 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 700 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 701 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 702 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 703 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 704 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 705 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 706 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 707 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 708 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 709 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 710 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 711 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 712 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 713 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 714 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 715 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 716 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 717 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 718 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 719 * ccomps * dcomps);

            auto g_y_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 720 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyyyy = cbuffer.data(hi_geom_10_off + 721 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyyyz = cbuffer.data(hi_geom_10_off + 722 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyyzz = cbuffer.data(hi_geom_10_off + 723 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyzzz = cbuffer.data(hi_geom_10_off + 724 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyzzzz = cbuffer.data(hi_geom_10_off + 725 * ccomps * dcomps);

            auto g_y_0_xxxyz_yzzzzz = cbuffer.data(hi_geom_10_off + 726 * ccomps * dcomps);

            auto g_y_0_xxxyz_zzzzzz = cbuffer.data(hi_geom_10_off + 727 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyz_xxxxxx, g_y_0_xxxyz_xxxxxy, g_y_0_xxxyz_xxxxxz, g_y_0_xxxyz_xxxxyy, g_y_0_xxxyz_xxxxyz, g_y_0_xxxyz_xxxxzz, g_y_0_xxxyz_xxxyyy, g_y_0_xxxyz_xxxyyz, g_y_0_xxxyz_xxxyzz, g_y_0_xxxyz_xxxzzz, g_y_0_xxxyz_xxyyyy, g_y_0_xxxyz_xxyyyz, g_y_0_xxxyz_xxyyzz, g_y_0_xxxyz_xxyzzz, g_y_0_xxxyz_xxzzzz, g_y_0_xxxyz_xyyyyy, g_y_0_xxxyz_xyyyyz, g_y_0_xxxyz_xyyyzz, g_y_0_xxxyz_xyyzzz, g_y_0_xxxyz_xyzzzz, g_y_0_xxxyz_xzzzzz, g_y_0_xxxyz_yyyyyy, g_y_0_xxxyz_yyyyyz, g_y_0_xxxyz_yyyyzz, g_y_0_xxxyz_yyyzzz, g_y_0_xxxyz_yyzzzz, g_y_0_xxxyz_yzzzzz, g_y_0_xxxyz_zzzzzz, g_y_0_xxyz_xxxxxx, g_y_0_xxyz_xxxxxxx, g_y_0_xxyz_xxxxxxy, g_y_0_xxyz_xxxxxxz, g_y_0_xxyz_xxxxxy, g_y_0_xxyz_xxxxxyy, g_y_0_xxyz_xxxxxyz, g_y_0_xxyz_xxxxxz, g_y_0_xxyz_xxxxxzz, g_y_0_xxyz_xxxxyy, g_y_0_xxyz_xxxxyyy, g_y_0_xxyz_xxxxyyz, g_y_0_xxyz_xxxxyz, g_y_0_xxyz_xxxxyzz, g_y_0_xxyz_xxxxzz, g_y_0_xxyz_xxxxzzz, g_y_0_xxyz_xxxyyy, g_y_0_xxyz_xxxyyyy, g_y_0_xxyz_xxxyyyz, g_y_0_xxyz_xxxyyz, g_y_0_xxyz_xxxyyzz, g_y_0_xxyz_xxxyzz, g_y_0_xxyz_xxxyzzz, g_y_0_xxyz_xxxzzz, g_y_0_xxyz_xxxzzzz, g_y_0_xxyz_xxyyyy, g_y_0_xxyz_xxyyyyy, g_y_0_xxyz_xxyyyyz, g_y_0_xxyz_xxyyyz, g_y_0_xxyz_xxyyyzz, g_y_0_xxyz_xxyyzz, g_y_0_xxyz_xxyyzzz, g_y_0_xxyz_xxyzzz, g_y_0_xxyz_xxyzzzz, g_y_0_xxyz_xxzzzz, g_y_0_xxyz_xxzzzzz, g_y_0_xxyz_xyyyyy, g_y_0_xxyz_xyyyyyy, g_y_0_xxyz_xyyyyyz, g_y_0_xxyz_xyyyyz, g_y_0_xxyz_xyyyyzz, g_y_0_xxyz_xyyyzz, g_y_0_xxyz_xyyyzzz, g_y_0_xxyz_xyyzzz, g_y_0_xxyz_xyyzzzz, g_y_0_xxyz_xyzzzz, g_y_0_xxyz_xyzzzzz, g_y_0_xxyz_xzzzzz, g_y_0_xxyz_xzzzzzz, g_y_0_xxyz_yyyyyy, g_y_0_xxyz_yyyyyz, g_y_0_xxyz_yyyyzz, g_y_0_xxyz_yyyzzz, g_y_0_xxyz_yyzzzz, g_y_0_xxyz_yzzzzz, g_y_0_xxyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyz_xxxxxx[k] = -g_y_0_xxyz_xxxxxx[k] * ab_x + g_y_0_xxyz_xxxxxxx[k];

                g_y_0_xxxyz_xxxxxy[k] = -g_y_0_xxyz_xxxxxy[k] * ab_x + g_y_0_xxyz_xxxxxxy[k];

                g_y_0_xxxyz_xxxxxz[k] = -g_y_0_xxyz_xxxxxz[k] * ab_x + g_y_0_xxyz_xxxxxxz[k];

                g_y_0_xxxyz_xxxxyy[k] = -g_y_0_xxyz_xxxxyy[k] * ab_x + g_y_0_xxyz_xxxxxyy[k];

                g_y_0_xxxyz_xxxxyz[k] = -g_y_0_xxyz_xxxxyz[k] * ab_x + g_y_0_xxyz_xxxxxyz[k];

                g_y_0_xxxyz_xxxxzz[k] = -g_y_0_xxyz_xxxxzz[k] * ab_x + g_y_0_xxyz_xxxxxzz[k];

                g_y_0_xxxyz_xxxyyy[k] = -g_y_0_xxyz_xxxyyy[k] * ab_x + g_y_0_xxyz_xxxxyyy[k];

                g_y_0_xxxyz_xxxyyz[k] = -g_y_0_xxyz_xxxyyz[k] * ab_x + g_y_0_xxyz_xxxxyyz[k];

                g_y_0_xxxyz_xxxyzz[k] = -g_y_0_xxyz_xxxyzz[k] * ab_x + g_y_0_xxyz_xxxxyzz[k];

                g_y_0_xxxyz_xxxzzz[k] = -g_y_0_xxyz_xxxzzz[k] * ab_x + g_y_0_xxyz_xxxxzzz[k];

                g_y_0_xxxyz_xxyyyy[k] = -g_y_0_xxyz_xxyyyy[k] * ab_x + g_y_0_xxyz_xxxyyyy[k];

                g_y_0_xxxyz_xxyyyz[k] = -g_y_0_xxyz_xxyyyz[k] * ab_x + g_y_0_xxyz_xxxyyyz[k];

                g_y_0_xxxyz_xxyyzz[k] = -g_y_0_xxyz_xxyyzz[k] * ab_x + g_y_0_xxyz_xxxyyzz[k];

                g_y_0_xxxyz_xxyzzz[k] = -g_y_0_xxyz_xxyzzz[k] * ab_x + g_y_0_xxyz_xxxyzzz[k];

                g_y_0_xxxyz_xxzzzz[k] = -g_y_0_xxyz_xxzzzz[k] * ab_x + g_y_0_xxyz_xxxzzzz[k];

                g_y_0_xxxyz_xyyyyy[k] = -g_y_0_xxyz_xyyyyy[k] * ab_x + g_y_0_xxyz_xxyyyyy[k];

                g_y_0_xxxyz_xyyyyz[k] = -g_y_0_xxyz_xyyyyz[k] * ab_x + g_y_0_xxyz_xxyyyyz[k];

                g_y_0_xxxyz_xyyyzz[k] = -g_y_0_xxyz_xyyyzz[k] * ab_x + g_y_0_xxyz_xxyyyzz[k];

                g_y_0_xxxyz_xyyzzz[k] = -g_y_0_xxyz_xyyzzz[k] * ab_x + g_y_0_xxyz_xxyyzzz[k];

                g_y_0_xxxyz_xyzzzz[k] = -g_y_0_xxyz_xyzzzz[k] * ab_x + g_y_0_xxyz_xxyzzzz[k];

                g_y_0_xxxyz_xzzzzz[k] = -g_y_0_xxyz_xzzzzz[k] * ab_x + g_y_0_xxyz_xxzzzzz[k];

                g_y_0_xxxyz_yyyyyy[k] = -g_y_0_xxyz_yyyyyy[k] * ab_x + g_y_0_xxyz_xyyyyyy[k];

                g_y_0_xxxyz_yyyyyz[k] = -g_y_0_xxyz_yyyyyz[k] * ab_x + g_y_0_xxyz_xyyyyyz[k];

                g_y_0_xxxyz_yyyyzz[k] = -g_y_0_xxyz_yyyyzz[k] * ab_x + g_y_0_xxyz_xyyyyzz[k];

                g_y_0_xxxyz_yyyzzz[k] = -g_y_0_xxyz_yyyzzz[k] * ab_x + g_y_0_xxyz_xyyyzzz[k];

                g_y_0_xxxyz_yyzzzz[k] = -g_y_0_xxyz_yyzzzz[k] * ab_x + g_y_0_xxyz_xyyzzzz[k];

                g_y_0_xxxyz_yzzzzz[k] = -g_y_0_xxyz_yzzzzz[k] * ab_x + g_y_0_xxyz_xyzzzzz[k];

                g_y_0_xxxyz_zzzzzz[k] = -g_y_0_xxyz_zzzzzz[k] * ab_x + g_y_0_xxyz_xzzzzzz[k];
            }

            /// Set up 728-756 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 728 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 729 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 730 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 731 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 732 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 733 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 734 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 735 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 736 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 737 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 738 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 739 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 740 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 741 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 742 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 743 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 744 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 745 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 746 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 747 * ccomps * dcomps);

            auto g_y_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 748 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyyyy = cbuffer.data(hi_geom_10_off + 749 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyyyz = cbuffer.data(hi_geom_10_off + 750 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyyzz = cbuffer.data(hi_geom_10_off + 751 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyzzz = cbuffer.data(hi_geom_10_off + 752 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyzzzz = cbuffer.data(hi_geom_10_off + 753 * ccomps * dcomps);

            auto g_y_0_xxxzz_yzzzzz = cbuffer.data(hi_geom_10_off + 754 * ccomps * dcomps);

            auto g_y_0_xxxzz_zzzzzz = cbuffer.data(hi_geom_10_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxzz_xxxxxx, g_y_0_xxxzz_xxxxxy, g_y_0_xxxzz_xxxxxz, g_y_0_xxxzz_xxxxyy, g_y_0_xxxzz_xxxxyz, g_y_0_xxxzz_xxxxzz, g_y_0_xxxzz_xxxyyy, g_y_0_xxxzz_xxxyyz, g_y_0_xxxzz_xxxyzz, g_y_0_xxxzz_xxxzzz, g_y_0_xxxzz_xxyyyy, g_y_0_xxxzz_xxyyyz, g_y_0_xxxzz_xxyyzz, g_y_0_xxxzz_xxyzzz, g_y_0_xxxzz_xxzzzz, g_y_0_xxxzz_xyyyyy, g_y_0_xxxzz_xyyyyz, g_y_0_xxxzz_xyyyzz, g_y_0_xxxzz_xyyzzz, g_y_0_xxxzz_xyzzzz, g_y_0_xxxzz_xzzzzz, g_y_0_xxxzz_yyyyyy, g_y_0_xxxzz_yyyyyz, g_y_0_xxxzz_yyyyzz, g_y_0_xxxzz_yyyzzz, g_y_0_xxxzz_yyzzzz, g_y_0_xxxzz_yzzzzz, g_y_0_xxxzz_zzzzzz, g_y_0_xxzz_xxxxxx, g_y_0_xxzz_xxxxxxx, g_y_0_xxzz_xxxxxxy, g_y_0_xxzz_xxxxxxz, g_y_0_xxzz_xxxxxy, g_y_0_xxzz_xxxxxyy, g_y_0_xxzz_xxxxxyz, g_y_0_xxzz_xxxxxz, g_y_0_xxzz_xxxxxzz, g_y_0_xxzz_xxxxyy, g_y_0_xxzz_xxxxyyy, g_y_0_xxzz_xxxxyyz, g_y_0_xxzz_xxxxyz, g_y_0_xxzz_xxxxyzz, g_y_0_xxzz_xxxxzz, g_y_0_xxzz_xxxxzzz, g_y_0_xxzz_xxxyyy, g_y_0_xxzz_xxxyyyy, g_y_0_xxzz_xxxyyyz, g_y_0_xxzz_xxxyyz, g_y_0_xxzz_xxxyyzz, g_y_0_xxzz_xxxyzz, g_y_0_xxzz_xxxyzzz, g_y_0_xxzz_xxxzzz, g_y_0_xxzz_xxxzzzz, g_y_0_xxzz_xxyyyy, g_y_0_xxzz_xxyyyyy, g_y_0_xxzz_xxyyyyz, g_y_0_xxzz_xxyyyz, g_y_0_xxzz_xxyyyzz, g_y_0_xxzz_xxyyzz, g_y_0_xxzz_xxyyzzz, g_y_0_xxzz_xxyzzz, g_y_0_xxzz_xxyzzzz, g_y_0_xxzz_xxzzzz, g_y_0_xxzz_xxzzzzz, g_y_0_xxzz_xyyyyy, g_y_0_xxzz_xyyyyyy, g_y_0_xxzz_xyyyyyz, g_y_0_xxzz_xyyyyz, g_y_0_xxzz_xyyyyzz, g_y_0_xxzz_xyyyzz, g_y_0_xxzz_xyyyzzz, g_y_0_xxzz_xyyzzz, g_y_0_xxzz_xyyzzzz, g_y_0_xxzz_xyzzzz, g_y_0_xxzz_xyzzzzz, g_y_0_xxzz_xzzzzz, g_y_0_xxzz_xzzzzzz, g_y_0_xxzz_yyyyyy, g_y_0_xxzz_yyyyyz, g_y_0_xxzz_yyyyzz, g_y_0_xxzz_yyyzzz, g_y_0_xxzz_yyzzzz, g_y_0_xxzz_yzzzzz, g_y_0_xxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzz_xxxxxx[k] = -g_y_0_xxzz_xxxxxx[k] * ab_x + g_y_0_xxzz_xxxxxxx[k];

                g_y_0_xxxzz_xxxxxy[k] = -g_y_0_xxzz_xxxxxy[k] * ab_x + g_y_0_xxzz_xxxxxxy[k];

                g_y_0_xxxzz_xxxxxz[k] = -g_y_0_xxzz_xxxxxz[k] * ab_x + g_y_0_xxzz_xxxxxxz[k];

                g_y_0_xxxzz_xxxxyy[k] = -g_y_0_xxzz_xxxxyy[k] * ab_x + g_y_0_xxzz_xxxxxyy[k];

                g_y_0_xxxzz_xxxxyz[k] = -g_y_0_xxzz_xxxxyz[k] * ab_x + g_y_0_xxzz_xxxxxyz[k];

                g_y_0_xxxzz_xxxxzz[k] = -g_y_0_xxzz_xxxxzz[k] * ab_x + g_y_0_xxzz_xxxxxzz[k];

                g_y_0_xxxzz_xxxyyy[k] = -g_y_0_xxzz_xxxyyy[k] * ab_x + g_y_0_xxzz_xxxxyyy[k];

                g_y_0_xxxzz_xxxyyz[k] = -g_y_0_xxzz_xxxyyz[k] * ab_x + g_y_0_xxzz_xxxxyyz[k];

                g_y_0_xxxzz_xxxyzz[k] = -g_y_0_xxzz_xxxyzz[k] * ab_x + g_y_0_xxzz_xxxxyzz[k];

                g_y_0_xxxzz_xxxzzz[k] = -g_y_0_xxzz_xxxzzz[k] * ab_x + g_y_0_xxzz_xxxxzzz[k];

                g_y_0_xxxzz_xxyyyy[k] = -g_y_0_xxzz_xxyyyy[k] * ab_x + g_y_0_xxzz_xxxyyyy[k];

                g_y_0_xxxzz_xxyyyz[k] = -g_y_0_xxzz_xxyyyz[k] * ab_x + g_y_0_xxzz_xxxyyyz[k];

                g_y_0_xxxzz_xxyyzz[k] = -g_y_0_xxzz_xxyyzz[k] * ab_x + g_y_0_xxzz_xxxyyzz[k];

                g_y_0_xxxzz_xxyzzz[k] = -g_y_0_xxzz_xxyzzz[k] * ab_x + g_y_0_xxzz_xxxyzzz[k];

                g_y_0_xxxzz_xxzzzz[k] = -g_y_0_xxzz_xxzzzz[k] * ab_x + g_y_0_xxzz_xxxzzzz[k];

                g_y_0_xxxzz_xyyyyy[k] = -g_y_0_xxzz_xyyyyy[k] * ab_x + g_y_0_xxzz_xxyyyyy[k];

                g_y_0_xxxzz_xyyyyz[k] = -g_y_0_xxzz_xyyyyz[k] * ab_x + g_y_0_xxzz_xxyyyyz[k];

                g_y_0_xxxzz_xyyyzz[k] = -g_y_0_xxzz_xyyyzz[k] * ab_x + g_y_0_xxzz_xxyyyzz[k];

                g_y_0_xxxzz_xyyzzz[k] = -g_y_0_xxzz_xyyzzz[k] * ab_x + g_y_0_xxzz_xxyyzzz[k];

                g_y_0_xxxzz_xyzzzz[k] = -g_y_0_xxzz_xyzzzz[k] * ab_x + g_y_0_xxzz_xxyzzzz[k];

                g_y_0_xxxzz_xzzzzz[k] = -g_y_0_xxzz_xzzzzz[k] * ab_x + g_y_0_xxzz_xxzzzzz[k];

                g_y_0_xxxzz_yyyyyy[k] = -g_y_0_xxzz_yyyyyy[k] * ab_x + g_y_0_xxzz_xyyyyyy[k];

                g_y_0_xxxzz_yyyyyz[k] = -g_y_0_xxzz_yyyyyz[k] * ab_x + g_y_0_xxzz_xyyyyyz[k];

                g_y_0_xxxzz_yyyyzz[k] = -g_y_0_xxzz_yyyyzz[k] * ab_x + g_y_0_xxzz_xyyyyzz[k];

                g_y_0_xxxzz_yyyzzz[k] = -g_y_0_xxzz_yyyzzz[k] * ab_x + g_y_0_xxzz_xyyyzzz[k];

                g_y_0_xxxzz_yyzzzz[k] = -g_y_0_xxzz_yyzzzz[k] * ab_x + g_y_0_xxzz_xyyzzzz[k];

                g_y_0_xxxzz_yzzzzz[k] = -g_y_0_xxzz_yzzzzz[k] * ab_x + g_y_0_xxzz_xyzzzzz[k];

                g_y_0_xxxzz_zzzzzz[k] = -g_y_0_xxzz_zzzzzz[k] * ab_x + g_y_0_xxzz_xzzzzzz[k];
            }

            /// Set up 756-784 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 756 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 757 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 758 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 759 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 760 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 761 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 762 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 763 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 764 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 765 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 766 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 767 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 768 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 769 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 770 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 771 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 772 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 773 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 774 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 775 * ccomps * dcomps);

            auto g_y_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 776 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 777 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 778 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 779 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 780 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 781 * ccomps * dcomps);

            auto g_y_0_xxyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 782 * ccomps * dcomps);

            auto g_y_0_xxyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 783 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyy_xxxxxx, g_y_0_xxyyy_xxxxxy, g_y_0_xxyyy_xxxxxz, g_y_0_xxyyy_xxxxyy, g_y_0_xxyyy_xxxxyz, g_y_0_xxyyy_xxxxzz, g_y_0_xxyyy_xxxyyy, g_y_0_xxyyy_xxxyyz, g_y_0_xxyyy_xxxyzz, g_y_0_xxyyy_xxxzzz, g_y_0_xxyyy_xxyyyy, g_y_0_xxyyy_xxyyyz, g_y_0_xxyyy_xxyyzz, g_y_0_xxyyy_xxyzzz, g_y_0_xxyyy_xxzzzz, g_y_0_xxyyy_xyyyyy, g_y_0_xxyyy_xyyyyz, g_y_0_xxyyy_xyyyzz, g_y_0_xxyyy_xyyzzz, g_y_0_xxyyy_xyzzzz, g_y_0_xxyyy_xzzzzz, g_y_0_xxyyy_yyyyyy, g_y_0_xxyyy_yyyyyz, g_y_0_xxyyy_yyyyzz, g_y_0_xxyyy_yyyzzz, g_y_0_xxyyy_yyzzzz, g_y_0_xxyyy_yzzzzz, g_y_0_xxyyy_zzzzzz, g_y_0_xyyy_xxxxxx, g_y_0_xyyy_xxxxxxx, g_y_0_xyyy_xxxxxxy, g_y_0_xyyy_xxxxxxz, g_y_0_xyyy_xxxxxy, g_y_0_xyyy_xxxxxyy, g_y_0_xyyy_xxxxxyz, g_y_0_xyyy_xxxxxz, g_y_0_xyyy_xxxxxzz, g_y_0_xyyy_xxxxyy, g_y_0_xyyy_xxxxyyy, g_y_0_xyyy_xxxxyyz, g_y_0_xyyy_xxxxyz, g_y_0_xyyy_xxxxyzz, g_y_0_xyyy_xxxxzz, g_y_0_xyyy_xxxxzzz, g_y_0_xyyy_xxxyyy, g_y_0_xyyy_xxxyyyy, g_y_0_xyyy_xxxyyyz, g_y_0_xyyy_xxxyyz, g_y_0_xyyy_xxxyyzz, g_y_0_xyyy_xxxyzz, g_y_0_xyyy_xxxyzzz, g_y_0_xyyy_xxxzzz, g_y_0_xyyy_xxxzzzz, g_y_0_xyyy_xxyyyy, g_y_0_xyyy_xxyyyyy, g_y_0_xyyy_xxyyyyz, g_y_0_xyyy_xxyyyz, g_y_0_xyyy_xxyyyzz, g_y_0_xyyy_xxyyzz, g_y_0_xyyy_xxyyzzz, g_y_0_xyyy_xxyzzz, g_y_0_xyyy_xxyzzzz, g_y_0_xyyy_xxzzzz, g_y_0_xyyy_xxzzzzz, g_y_0_xyyy_xyyyyy, g_y_0_xyyy_xyyyyyy, g_y_0_xyyy_xyyyyyz, g_y_0_xyyy_xyyyyz, g_y_0_xyyy_xyyyyzz, g_y_0_xyyy_xyyyzz, g_y_0_xyyy_xyyyzzz, g_y_0_xyyy_xyyzzz, g_y_0_xyyy_xyyzzzz, g_y_0_xyyy_xyzzzz, g_y_0_xyyy_xyzzzzz, g_y_0_xyyy_xzzzzz, g_y_0_xyyy_xzzzzzz, g_y_0_xyyy_yyyyyy, g_y_0_xyyy_yyyyyz, g_y_0_xyyy_yyyyzz, g_y_0_xyyy_yyyzzz, g_y_0_xyyy_yyzzzz, g_y_0_xyyy_yzzzzz, g_y_0_xyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyy_xxxxxx[k] = -g_y_0_xyyy_xxxxxx[k] * ab_x + g_y_0_xyyy_xxxxxxx[k];

                g_y_0_xxyyy_xxxxxy[k] = -g_y_0_xyyy_xxxxxy[k] * ab_x + g_y_0_xyyy_xxxxxxy[k];

                g_y_0_xxyyy_xxxxxz[k] = -g_y_0_xyyy_xxxxxz[k] * ab_x + g_y_0_xyyy_xxxxxxz[k];

                g_y_0_xxyyy_xxxxyy[k] = -g_y_0_xyyy_xxxxyy[k] * ab_x + g_y_0_xyyy_xxxxxyy[k];

                g_y_0_xxyyy_xxxxyz[k] = -g_y_0_xyyy_xxxxyz[k] * ab_x + g_y_0_xyyy_xxxxxyz[k];

                g_y_0_xxyyy_xxxxzz[k] = -g_y_0_xyyy_xxxxzz[k] * ab_x + g_y_0_xyyy_xxxxxzz[k];

                g_y_0_xxyyy_xxxyyy[k] = -g_y_0_xyyy_xxxyyy[k] * ab_x + g_y_0_xyyy_xxxxyyy[k];

                g_y_0_xxyyy_xxxyyz[k] = -g_y_0_xyyy_xxxyyz[k] * ab_x + g_y_0_xyyy_xxxxyyz[k];

                g_y_0_xxyyy_xxxyzz[k] = -g_y_0_xyyy_xxxyzz[k] * ab_x + g_y_0_xyyy_xxxxyzz[k];

                g_y_0_xxyyy_xxxzzz[k] = -g_y_0_xyyy_xxxzzz[k] * ab_x + g_y_0_xyyy_xxxxzzz[k];

                g_y_0_xxyyy_xxyyyy[k] = -g_y_0_xyyy_xxyyyy[k] * ab_x + g_y_0_xyyy_xxxyyyy[k];

                g_y_0_xxyyy_xxyyyz[k] = -g_y_0_xyyy_xxyyyz[k] * ab_x + g_y_0_xyyy_xxxyyyz[k];

                g_y_0_xxyyy_xxyyzz[k] = -g_y_0_xyyy_xxyyzz[k] * ab_x + g_y_0_xyyy_xxxyyzz[k];

                g_y_0_xxyyy_xxyzzz[k] = -g_y_0_xyyy_xxyzzz[k] * ab_x + g_y_0_xyyy_xxxyzzz[k];

                g_y_0_xxyyy_xxzzzz[k] = -g_y_0_xyyy_xxzzzz[k] * ab_x + g_y_0_xyyy_xxxzzzz[k];

                g_y_0_xxyyy_xyyyyy[k] = -g_y_0_xyyy_xyyyyy[k] * ab_x + g_y_0_xyyy_xxyyyyy[k];

                g_y_0_xxyyy_xyyyyz[k] = -g_y_0_xyyy_xyyyyz[k] * ab_x + g_y_0_xyyy_xxyyyyz[k];

                g_y_0_xxyyy_xyyyzz[k] = -g_y_0_xyyy_xyyyzz[k] * ab_x + g_y_0_xyyy_xxyyyzz[k];

                g_y_0_xxyyy_xyyzzz[k] = -g_y_0_xyyy_xyyzzz[k] * ab_x + g_y_0_xyyy_xxyyzzz[k];

                g_y_0_xxyyy_xyzzzz[k] = -g_y_0_xyyy_xyzzzz[k] * ab_x + g_y_0_xyyy_xxyzzzz[k];

                g_y_0_xxyyy_xzzzzz[k] = -g_y_0_xyyy_xzzzzz[k] * ab_x + g_y_0_xyyy_xxzzzzz[k];

                g_y_0_xxyyy_yyyyyy[k] = -g_y_0_xyyy_yyyyyy[k] * ab_x + g_y_0_xyyy_xyyyyyy[k];

                g_y_0_xxyyy_yyyyyz[k] = -g_y_0_xyyy_yyyyyz[k] * ab_x + g_y_0_xyyy_xyyyyyz[k];

                g_y_0_xxyyy_yyyyzz[k] = -g_y_0_xyyy_yyyyzz[k] * ab_x + g_y_0_xyyy_xyyyyzz[k];

                g_y_0_xxyyy_yyyzzz[k] = -g_y_0_xyyy_yyyzzz[k] * ab_x + g_y_0_xyyy_xyyyzzz[k];

                g_y_0_xxyyy_yyzzzz[k] = -g_y_0_xyyy_yyzzzz[k] * ab_x + g_y_0_xyyy_xyyzzzz[k];

                g_y_0_xxyyy_yzzzzz[k] = -g_y_0_xyyy_yzzzzz[k] * ab_x + g_y_0_xyyy_xyzzzzz[k];

                g_y_0_xxyyy_zzzzzz[k] = -g_y_0_xyyy_zzzzzz[k] * ab_x + g_y_0_xyyy_xzzzzzz[k];
            }

            /// Set up 784-812 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 784 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 785 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 786 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 787 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 788 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 789 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 790 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 791 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 794 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 797 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 805 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 809 * ccomps * dcomps);

            auto g_y_0_xxyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_xxyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 811 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyz_xxxxxx, g_y_0_xxyyz_xxxxxy, g_y_0_xxyyz_xxxxxz, g_y_0_xxyyz_xxxxyy, g_y_0_xxyyz_xxxxyz, g_y_0_xxyyz_xxxxzz, g_y_0_xxyyz_xxxyyy, g_y_0_xxyyz_xxxyyz, g_y_0_xxyyz_xxxyzz, g_y_0_xxyyz_xxxzzz, g_y_0_xxyyz_xxyyyy, g_y_0_xxyyz_xxyyyz, g_y_0_xxyyz_xxyyzz, g_y_0_xxyyz_xxyzzz, g_y_0_xxyyz_xxzzzz, g_y_0_xxyyz_xyyyyy, g_y_0_xxyyz_xyyyyz, g_y_0_xxyyz_xyyyzz, g_y_0_xxyyz_xyyzzz, g_y_0_xxyyz_xyzzzz, g_y_0_xxyyz_xzzzzz, g_y_0_xxyyz_yyyyyy, g_y_0_xxyyz_yyyyyz, g_y_0_xxyyz_yyyyzz, g_y_0_xxyyz_yyyzzz, g_y_0_xxyyz_yyzzzz, g_y_0_xxyyz_yzzzzz, g_y_0_xxyyz_zzzzzz, g_y_0_xyyz_xxxxxx, g_y_0_xyyz_xxxxxxx, g_y_0_xyyz_xxxxxxy, g_y_0_xyyz_xxxxxxz, g_y_0_xyyz_xxxxxy, g_y_0_xyyz_xxxxxyy, g_y_0_xyyz_xxxxxyz, g_y_0_xyyz_xxxxxz, g_y_0_xyyz_xxxxxzz, g_y_0_xyyz_xxxxyy, g_y_0_xyyz_xxxxyyy, g_y_0_xyyz_xxxxyyz, g_y_0_xyyz_xxxxyz, g_y_0_xyyz_xxxxyzz, g_y_0_xyyz_xxxxzz, g_y_0_xyyz_xxxxzzz, g_y_0_xyyz_xxxyyy, g_y_0_xyyz_xxxyyyy, g_y_0_xyyz_xxxyyyz, g_y_0_xyyz_xxxyyz, g_y_0_xyyz_xxxyyzz, g_y_0_xyyz_xxxyzz, g_y_0_xyyz_xxxyzzz, g_y_0_xyyz_xxxzzz, g_y_0_xyyz_xxxzzzz, g_y_0_xyyz_xxyyyy, g_y_0_xyyz_xxyyyyy, g_y_0_xyyz_xxyyyyz, g_y_0_xyyz_xxyyyz, g_y_0_xyyz_xxyyyzz, g_y_0_xyyz_xxyyzz, g_y_0_xyyz_xxyyzzz, g_y_0_xyyz_xxyzzz, g_y_0_xyyz_xxyzzzz, g_y_0_xyyz_xxzzzz, g_y_0_xyyz_xxzzzzz, g_y_0_xyyz_xyyyyy, g_y_0_xyyz_xyyyyyy, g_y_0_xyyz_xyyyyyz, g_y_0_xyyz_xyyyyz, g_y_0_xyyz_xyyyyzz, g_y_0_xyyz_xyyyzz, g_y_0_xyyz_xyyyzzz, g_y_0_xyyz_xyyzzz, g_y_0_xyyz_xyyzzzz, g_y_0_xyyz_xyzzzz, g_y_0_xyyz_xyzzzzz, g_y_0_xyyz_xzzzzz, g_y_0_xyyz_xzzzzzz, g_y_0_xyyz_yyyyyy, g_y_0_xyyz_yyyyyz, g_y_0_xyyz_yyyyzz, g_y_0_xyyz_yyyzzz, g_y_0_xyyz_yyzzzz, g_y_0_xyyz_yzzzzz, g_y_0_xyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyz_xxxxxx[k] = -g_y_0_xyyz_xxxxxx[k] * ab_x + g_y_0_xyyz_xxxxxxx[k];

                g_y_0_xxyyz_xxxxxy[k] = -g_y_0_xyyz_xxxxxy[k] * ab_x + g_y_0_xyyz_xxxxxxy[k];

                g_y_0_xxyyz_xxxxxz[k] = -g_y_0_xyyz_xxxxxz[k] * ab_x + g_y_0_xyyz_xxxxxxz[k];

                g_y_0_xxyyz_xxxxyy[k] = -g_y_0_xyyz_xxxxyy[k] * ab_x + g_y_0_xyyz_xxxxxyy[k];

                g_y_0_xxyyz_xxxxyz[k] = -g_y_0_xyyz_xxxxyz[k] * ab_x + g_y_0_xyyz_xxxxxyz[k];

                g_y_0_xxyyz_xxxxzz[k] = -g_y_0_xyyz_xxxxzz[k] * ab_x + g_y_0_xyyz_xxxxxzz[k];

                g_y_0_xxyyz_xxxyyy[k] = -g_y_0_xyyz_xxxyyy[k] * ab_x + g_y_0_xyyz_xxxxyyy[k];

                g_y_0_xxyyz_xxxyyz[k] = -g_y_0_xyyz_xxxyyz[k] * ab_x + g_y_0_xyyz_xxxxyyz[k];

                g_y_0_xxyyz_xxxyzz[k] = -g_y_0_xyyz_xxxyzz[k] * ab_x + g_y_0_xyyz_xxxxyzz[k];

                g_y_0_xxyyz_xxxzzz[k] = -g_y_0_xyyz_xxxzzz[k] * ab_x + g_y_0_xyyz_xxxxzzz[k];

                g_y_0_xxyyz_xxyyyy[k] = -g_y_0_xyyz_xxyyyy[k] * ab_x + g_y_0_xyyz_xxxyyyy[k];

                g_y_0_xxyyz_xxyyyz[k] = -g_y_0_xyyz_xxyyyz[k] * ab_x + g_y_0_xyyz_xxxyyyz[k];

                g_y_0_xxyyz_xxyyzz[k] = -g_y_0_xyyz_xxyyzz[k] * ab_x + g_y_0_xyyz_xxxyyzz[k];

                g_y_0_xxyyz_xxyzzz[k] = -g_y_0_xyyz_xxyzzz[k] * ab_x + g_y_0_xyyz_xxxyzzz[k];

                g_y_0_xxyyz_xxzzzz[k] = -g_y_0_xyyz_xxzzzz[k] * ab_x + g_y_0_xyyz_xxxzzzz[k];

                g_y_0_xxyyz_xyyyyy[k] = -g_y_0_xyyz_xyyyyy[k] * ab_x + g_y_0_xyyz_xxyyyyy[k];

                g_y_0_xxyyz_xyyyyz[k] = -g_y_0_xyyz_xyyyyz[k] * ab_x + g_y_0_xyyz_xxyyyyz[k];

                g_y_0_xxyyz_xyyyzz[k] = -g_y_0_xyyz_xyyyzz[k] * ab_x + g_y_0_xyyz_xxyyyzz[k];

                g_y_0_xxyyz_xyyzzz[k] = -g_y_0_xyyz_xyyzzz[k] * ab_x + g_y_0_xyyz_xxyyzzz[k];

                g_y_0_xxyyz_xyzzzz[k] = -g_y_0_xyyz_xyzzzz[k] * ab_x + g_y_0_xyyz_xxyzzzz[k];

                g_y_0_xxyyz_xzzzzz[k] = -g_y_0_xyyz_xzzzzz[k] * ab_x + g_y_0_xyyz_xxzzzzz[k];

                g_y_0_xxyyz_yyyyyy[k] = -g_y_0_xyyz_yyyyyy[k] * ab_x + g_y_0_xyyz_xyyyyyy[k];

                g_y_0_xxyyz_yyyyyz[k] = -g_y_0_xyyz_yyyyyz[k] * ab_x + g_y_0_xyyz_xyyyyyz[k];

                g_y_0_xxyyz_yyyyzz[k] = -g_y_0_xyyz_yyyyzz[k] * ab_x + g_y_0_xyyz_xyyyyzz[k];

                g_y_0_xxyyz_yyyzzz[k] = -g_y_0_xyyz_yyyzzz[k] * ab_x + g_y_0_xyyz_xyyyzzz[k];

                g_y_0_xxyyz_yyzzzz[k] = -g_y_0_xyyz_yyzzzz[k] * ab_x + g_y_0_xyyz_xyyzzzz[k];

                g_y_0_xxyyz_yzzzzz[k] = -g_y_0_xyyz_yzzzzz[k] * ab_x + g_y_0_xyyz_xyzzzzz[k];

                g_y_0_xxyyz_zzzzzz[k] = -g_y_0_xyyz_zzzzzz[k] * ab_x + g_y_0_xyyz_xzzzzzz[k];
            }

            /// Set up 812-840 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 818 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 820 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 821 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 822 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 823 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 824 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 825 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 826 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 827 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 833 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_xxyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_xxyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyzz_xxxxxx, g_y_0_xxyzz_xxxxxy, g_y_0_xxyzz_xxxxxz, g_y_0_xxyzz_xxxxyy, g_y_0_xxyzz_xxxxyz, g_y_0_xxyzz_xxxxzz, g_y_0_xxyzz_xxxyyy, g_y_0_xxyzz_xxxyyz, g_y_0_xxyzz_xxxyzz, g_y_0_xxyzz_xxxzzz, g_y_0_xxyzz_xxyyyy, g_y_0_xxyzz_xxyyyz, g_y_0_xxyzz_xxyyzz, g_y_0_xxyzz_xxyzzz, g_y_0_xxyzz_xxzzzz, g_y_0_xxyzz_xyyyyy, g_y_0_xxyzz_xyyyyz, g_y_0_xxyzz_xyyyzz, g_y_0_xxyzz_xyyzzz, g_y_0_xxyzz_xyzzzz, g_y_0_xxyzz_xzzzzz, g_y_0_xxyzz_yyyyyy, g_y_0_xxyzz_yyyyyz, g_y_0_xxyzz_yyyyzz, g_y_0_xxyzz_yyyzzz, g_y_0_xxyzz_yyzzzz, g_y_0_xxyzz_yzzzzz, g_y_0_xxyzz_zzzzzz, g_y_0_xyzz_xxxxxx, g_y_0_xyzz_xxxxxxx, g_y_0_xyzz_xxxxxxy, g_y_0_xyzz_xxxxxxz, g_y_0_xyzz_xxxxxy, g_y_0_xyzz_xxxxxyy, g_y_0_xyzz_xxxxxyz, g_y_0_xyzz_xxxxxz, g_y_0_xyzz_xxxxxzz, g_y_0_xyzz_xxxxyy, g_y_0_xyzz_xxxxyyy, g_y_0_xyzz_xxxxyyz, g_y_0_xyzz_xxxxyz, g_y_0_xyzz_xxxxyzz, g_y_0_xyzz_xxxxzz, g_y_0_xyzz_xxxxzzz, g_y_0_xyzz_xxxyyy, g_y_0_xyzz_xxxyyyy, g_y_0_xyzz_xxxyyyz, g_y_0_xyzz_xxxyyz, g_y_0_xyzz_xxxyyzz, g_y_0_xyzz_xxxyzz, g_y_0_xyzz_xxxyzzz, g_y_0_xyzz_xxxzzz, g_y_0_xyzz_xxxzzzz, g_y_0_xyzz_xxyyyy, g_y_0_xyzz_xxyyyyy, g_y_0_xyzz_xxyyyyz, g_y_0_xyzz_xxyyyz, g_y_0_xyzz_xxyyyzz, g_y_0_xyzz_xxyyzz, g_y_0_xyzz_xxyyzzz, g_y_0_xyzz_xxyzzz, g_y_0_xyzz_xxyzzzz, g_y_0_xyzz_xxzzzz, g_y_0_xyzz_xxzzzzz, g_y_0_xyzz_xyyyyy, g_y_0_xyzz_xyyyyyy, g_y_0_xyzz_xyyyyyz, g_y_0_xyzz_xyyyyz, g_y_0_xyzz_xyyyyzz, g_y_0_xyzz_xyyyzz, g_y_0_xyzz_xyyyzzz, g_y_0_xyzz_xyyzzz, g_y_0_xyzz_xyyzzzz, g_y_0_xyzz_xyzzzz, g_y_0_xyzz_xyzzzzz, g_y_0_xyzz_xzzzzz, g_y_0_xyzz_xzzzzzz, g_y_0_xyzz_yyyyyy, g_y_0_xyzz_yyyyyz, g_y_0_xyzz_yyyyzz, g_y_0_xyzz_yyyzzz, g_y_0_xyzz_yyzzzz, g_y_0_xyzz_yzzzzz, g_y_0_xyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzz_xxxxxx[k] = -g_y_0_xyzz_xxxxxx[k] * ab_x + g_y_0_xyzz_xxxxxxx[k];

                g_y_0_xxyzz_xxxxxy[k] = -g_y_0_xyzz_xxxxxy[k] * ab_x + g_y_0_xyzz_xxxxxxy[k];

                g_y_0_xxyzz_xxxxxz[k] = -g_y_0_xyzz_xxxxxz[k] * ab_x + g_y_0_xyzz_xxxxxxz[k];

                g_y_0_xxyzz_xxxxyy[k] = -g_y_0_xyzz_xxxxyy[k] * ab_x + g_y_0_xyzz_xxxxxyy[k];

                g_y_0_xxyzz_xxxxyz[k] = -g_y_0_xyzz_xxxxyz[k] * ab_x + g_y_0_xyzz_xxxxxyz[k];

                g_y_0_xxyzz_xxxxzz[k] = -g_y_0_xyzz_xxxxzz[k] * ab_x + g_y_0_xyzz_xxxxxzz[k];

                g_y_0_xxyzz_xxxyyy[k] = -g_y_0_xyzz_xxxyyy[k] * ab_x + g_y_0_xyzz_xxxxyyy[k];

                g_y_0_xxyzz_xxxyyz[k] = -g_y_0_xyzz_xxxyyz[k] * ab_x + g_y_0_xyzz_xxxxyyz[k];

                g_y_0_xxyzz_xxxyzz[k] = -g_y_0_xyzz_xxxyzz[k] * ab_x + g_y_0_xyzz_xxxxyzz[k];

                g_y_0_xxyzz_xxxzzz[k] = -g_y_0_xyzz_xxxzzz[k] * ab_x + g_y_0_xyzz_xxxxzzz[k];

                g_y_0_xxyzz_xxyyyy[k] = -g_y_0_xyzz_xxyyyy[k] * ab_x + g_y_0_xyzz_xxxyyyy[k];

                g_y_0_xxyzz_xxyyyz[k] = -g_y_0_xyzz_xxyyyz[k] * ab_x + g_y_0_xyzz_xxxyyyz[k];

                g_y_0_xxyzz_xxyyzz[k] = -g_y_0_xyzz_xxyyzz[k] * ab_x + g_y_0_xyzz_xxxyyzz[k];

                g_y_0_xxyzz_xxyzzz[k] = -g_y_0_xyzz_xxyzzz[k] * ab_x + g_y_0_xyzz_xxxyzzz[k];

                g_y_0_xxyzz_xxzzzz[k] = -g_y_0_xyzz_xxzzzz[k] * ab_x + g_y_0_xyzz_xxxzzzz[k];

                g_y_0_xxyzz_xyyyyy[k] = -g_y_0_xyzz_xyyyyy[k] * ab_x + g_y_0_xyzz_xxyyyyy[k];

                g_y_0_xxyzz_xyyyyz[k] = -g_y_0_xyzz_xyyyyz[k] * ab_x + g_y_0_xyzz_xxyyyyz[k];

                g_y_0_xxyzz_xyyyzz[k] = -g_y_0_xyzz_xyyyzz[k] * ab_x + g_y_0_xyzz_xxyyyzz[k];

                g_y_0_xxyzz_xyyzzz[k] = -g_y_0_xyzz_xyyzzz[k] * ab_x + g_y_0_xyzz_xxyyzzz[k];

                g_y_0_xxyzz_xyzzzz[k] = -g_y_0_xyzz_xyzzzz[k] * ab_x + g_y_0_xyzz_xxyzzzz[k];

                g_y_0_xxyzz_xzzzzz[k] = -g_y_0_xyzz_xzzzzz[k] * ab_x + g_y_0_xyzz_xxzzzzz[k];

                g_y_0_xxyzz_yyyyyy[k] = -g_y_0_xyzz_yyyyyy[k] * ab_x + g_y_0_xyzz_xyyyyyy[k];

                g_y_0_xxyzz_yyyyyz[k] = -g_y_0_xyzz_yyyyyz[k] * ab_x + g_y_0_xyzz_xyyyyyz[k];

                g_y_0_xxyzz_yyyyzz[k] = -g_y_0_xyzz_yyyyzz[k] * ab_x + g_y_0_xyzz_xyyyyzz[k];

                g_y_0_xxyzz_yyyzzz[k] = -g_y_0_xyzz_yyyzzz[k] * ab_x + g_y_0_xyzz_xyyyzzz[k];

                g_y_0_xxyzz_yyzzzz[k] = -g_y_0_xyzz_yyzzzz[k] * ab_x + g_y_0_xyzz_xyyzzzz[k];

                g_y_0_xxyzz_yzzzzz[k] = -g_y_0_xyzz_yzzzzz[k] * ab_x + g_y_0_xyzz_xyzzzzz[k];

                g_y_0_xxyzz_zzzzzz[k] = -g_y_0_xyzz_zzzzzz[k] * ab_x + g_y_0_xyzz_xzzzzzz[k];
            }

            /// Set up 840-868 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 840 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 841 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 842 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 843 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 844 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 845 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 846 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 847 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 848 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 849 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 850 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 851 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 852 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 853 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 854 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 855 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 856 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 857 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 858 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 859 * ccomps * dcomps);

            auto g_y_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 860 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 861 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 862 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 863 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 864 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 865 * ccomps * dcomps);

            auto g_y_0_xxzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 866 * ccomps * dcomps);

            auto g_y_0_xxzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 867 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzzz_xxxxxx, g_y_0_xxzzz_xxxxxy, g_y_0_xxzzz_xxxxxz, g_y_0_xxzzz_xxxxyy, g_y_0_xxzzz_xxxxyz, g_y_0_xxzzz_xxxxzz, g_y_0_xxzzz_xxxyyy, g_y_0_xxzzz_xxxyyz, g_y_0_xxzzz_xxxyzz, g_y_0_xxzzz_xxxzzz, g_y_0_xxzzz_xxyyyy, g_y_0_xxzzz_xxyyyz, g_y_0_xxzzz_xxyyzz, g_y_0_xxzzz_xxyzzz, g_y_0_xxzzz_xxzzzz, g_y_0_xxzzz_xyyyyy, g_y_0_xxzzz_xyyyyz, g_y_0_xxzzz_xyyyzz, g_y_0_xxzzz_xyyzzz, g_y_0_xxzzz_xyzzzz, g_y_0_xxzzz_xzzzzz, g_y_0_xxzzz_yyyyyy, g_y_0_xxzzz_yyyyyz, g_y_0_xxzzz_yyyyzz, g_y_0_xxzzz_yyyzzz, g_y_0_xxzzz_yyzzzz, g_y_0_xxzzz_yzzzzz, g_y_0_xxzzz_zzzzzz, g_y_0_xzzz_xxxxxx, g_y_0_xzzz_xxxxxxx, g_y_0_xzzz_xxxxxxy, g_y_0_xzzz_xxxxxxz, g_y_0_xzzz_xxxxxy, g_y_0_xzzz_xxxxxyy, g_y_0_xzzz_xxxxxyz, g_y_0_xzzz_xxxxxz, g_y_0_xzzz_xxxxxzz, g_y_0_xzzz_xxxxyy, g_y_0_xzzz_xxxxyyy, g_y_0_xzzz_xxxxyyz, g_y_0_xzzz_xxxxyz, g_y_0_xzzz_xxxxyzz, g_y_0_xzzz_xxxxzz, g_y_0_xzzz_xxxxzzz, g_y_0_xzzz_xxxyyy, g_y_0_xzzz_xxxyyyy, g_y_0_xzzz_xxxyyyz, g_y_0_xzzz_xxxyyz, g_y_0_xzzz_xxxyyzz, g_y_0_xzzz_xxxyzz, g_y_0_xzzz_xxxyzzz, g_y_0_xzzz_xxxzzz, g_y_0_xzzz_xxxzzzz, g_y_0_xzzz_xxyyyy, g_y_0_xzzz_xxyyyyy, g_y_0_xzzz_xxyyyyz, g_y_0_xzzz_xxyyyz, g_y_0_xzzz_xxyyyzz, g_y_0_xzzz_xxyyzz, g_y_0_xzzz_xxyyzzz, g_y_0_xzzz_xxyzzz, g_y_0_xzzz_xxyzzzz, g_y_0_xzzz_xxzzzz, g_y_0_xzzz_xxzzzzz, g_y_0_xzzz_xyyyyy, g_y_0_xzzz_xyyyyyy, g_y_0_xzzz_xyyyyyz, g_y_0_xzzz_xyyyyz, g_y_0_xzzz_xyyyyzz, g_y_0_xzzz_xyyyzz, g_y_0_xzzz_xyyyzzz, g_y_0_xzzz_xyyzzz, g_y_0_xzzz_xyyzzzz, g_y_0_xzzz_xyzzzz, g_y_0_xzzz_xyzzzzz, g_y_0_xzzz_xzzzzz, g_y_0_xzzz_xzzzzzz, g_y_0_xzzz_yyyyyy, g_y_0_xzzz_yyyyyz, g_y_0_xzzz_yyyyzz, g_y_0_xzzz_yyyzzz, g_y_0_xzzz_yyzzzz, g_y_0_xzzz_yzzzzz, g_y_0_xzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzz_xxxxxx[k] = -g_y_0_xzzz_xxxxxx[k] * ab_x + g_y_0_xzzz_xxxxxxx[k];

                g_y_0_xxzzz_xxxxxy[k] = -g_y_0_xzzz_xxxxxy[k] * ab_x + g_y_0_xzzz_xxxxxxy[k];

                g_y_0_xxzzz_xxxxxz[k] = -g_y_0_xzzz_xxxxxz[k] * ab_x + g_y_0_xzzz_xxxxxxz[k];

                g_y_0_xxzzz_xxxxyy[k] = -g_y_0_xzzz_xxxxyy[k] * ab_x + g_y_0_xzzz_xxxxxyy[k];

                g_y_0_xxzzz_xxxxyz[k] = -g_y_0_xzzz_xxxxyz[k] * ab_x + g_y_0_xzzz_xxxxxyz[k];

                g_y_0_xxzzz_xxxxzz[k] = -g_y_0_xzzz_xxxxzz[k] * ab_x + g_y_0_xzzz_xxxxxzz[k];

                g_y_0_xxzzz_xxxyyy[k] = -g_y_0_xzzz_xxxyyy[k] * ab_x + g_y_0_xzzz_xxxxyyy[k];

                g_y_0_xxzzz_xxxyyz[k] = -g_y_0_xzzz_xxxyyz[k] * ab_x + g_y_0_xzzz_xxxxyyz[k];

                g_y_0_xxzzz_xxxyzz[k] = -g_y_0_xzzz_xxxyzz[k] * ab_x + g_y_0_xzzz_xxxxyzz[k];

                g_y_0_xxzzz_xxxzzz[k] = -g_y_0_xzzz_xxxzzz[k] * ab_x + g_y_0_xzzz_xxxxzzz[k];

                g_y_0_xxzzz_xxyyyy[k] = -g_y_0_xzzz_xxyyyy[k] * ab_x + g_y_0_xzzz_xxxyyyy[k];

                g_y_0_xxzzz_xxyyyz[k] = -g_y_0_xzzz_xxyyyz[k] * ab_x + g_y_0_xzzz_xxxyyyz[k];

                g_y_0_xxzzz_xxyyzz[k] = -g_y_0_xzzz_xxyyzz[k] * ab_x + g_y_0_xzzz_xxxyyzz[k];

                g_y_0_xxzzz_xxyzzz[k] = -g_y_0_xzzz_xxyzzz[k] * ab_x + g_y_0_xzzz_xxxyzzz[k];

                g_y_0_xxzzz_xxzzzz[k] = -g_y_0_xzzz_xxzzzz[k] * ab_x + g_y_0_xzzz_xxxzzzz[k];

                g_y_0_xxzzz_xyyyyy[k] = -g_y_0_xzzz_xyyyyy[k] * ab_x + g_y_0_xzzz_xxyyyyy[k];

                g_y_0_xxzzz_xyyyyz[k] = -g_y_0_xzzz_xyyyyz[k] * ab_x + g_y_0_xzzz_xxyyyyz[k];

                g_y_0_xxzzz_xyyyzz[k] = -g_y_0_xzzz_xyyyzz[k] * ab_x + g_y_0_xzzz_xxyyyzz[k];

                g_y_0_xxzzz_xyyzzz[k] = -g_y_0_xzzz_xyyzzz[k] * ab_x + g_y_0_xzzz_xxyyzzz[k];

                g_y_0_xxzzz_xyzzzz[k] = -g_y_0_xzzz_xyzzzz[k] * ab_x + g_y_0_xzzz_xxyzzzz[k];

                g_y_0_xxzzz_xzzzzz[k] = -g_y_0_xzzz_xzzzzz[k] * ab_x + g_y_0_xzzz_xxzzzzz[k];

                g_y_0_xxzzz_yyyyyy[k] = -g_y_0_xzzz_yyyyyy[k] * ab_x + g_y_0_xzzz_xyyyyyy[k];

                g_y_0_xxzzz_yyyyyz[k] = -g_y_0_xzzz_yyyyyz[k] * ab_x + g_y_0_xzzz_xyyyyyz[k];

                g_y_0_xxzzz_yyyyzz[k] = -g_y_0_xzzz_yyyyzz[k] * ab_x + g_y_0_xzzz_xyyyyzz[k];

                g_y_0_xxzzz_yyyzzz[k] = -g_y_0_xzzz_yyyzzz[k] * ab_x + g_y_0_xzzz_xyyyzzz[k];

                g_y_0_xxzzz_yyzzzz[k] = -g_y_0_xzzz_yyzzzz[k] * ab_x + g_y_0_xzzz_xyyzzzz[k];

                g_y_0_xxzzz_yzzzzz[k] = -g_y_0_xzzz_yzzzzz[k] * ab_x + g_y_0_xzzz_xyzzzzz[k];

                g_y_0_xxzzz_zzzzzz[k] = -g_y_0_xzzz_zzzzzz[k] * ab_x + g_y_0_xzzz_xzzzzzz[k];
            }

            /// Set up 868-896 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 868 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 869 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 870 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 871 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 872 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 873 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 874 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 875 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 876 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 877 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 878 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 879 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 880 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 881 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 882 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 883 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 884 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 885 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 886 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 887 * ccomps * dcomps);

            auto g_y_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 888 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 889 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 890 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 891 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 892 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 893 * ccomps * dcomps);

            auto g_y_0_xyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 894 * ccomps * dcomps);

            auto g_y_0_xyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 895 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyy_xxxxxx, g_y_0_xyyyy_xxxxxy, g_y_0_xyyyy_xxxxxz, g_y_0_xyyyy_xxxxyy, g_y_0_xyyyy_xxxxyz, g_y_0_xyyyy_xxxxzz, g_y_0_xyyyy_xxxyyy, g_y_0_xyyyy_xxxyyz, g_y_0_xyyyy_xxxyzz, g_y_0_xyyyy_xxxzzz, g_y_0_xyyyy_xxyyyy, g_y_0_xyyyy_xxyyyz, g_y_0_xyyyy_xxyyzz, g_y_0_xyyyy_xxyzzz, g_y_0_xyyyy_xxzzzz, g_y_0_xyyyy_xyyyyy, g_y_0_xyyyy_xyyyyz, g_y_0_xyyyy_xyyyzz, g_y_0_xyyyy_xyyzzz, g_y_0_xyyyy_xyzzzz, g_y_0_xyyyy_xzzzzz, g_y_0_xyyyy_yyyyyy, g_y_0_xyyyy_yyyyyz, g_y_0_xyyyy_yyyyzz, g_y_0_xyyyy_yyyzzz, g_y_0_xyyyy_yyzzzz, g_y_0_xyyyy_yzzzzz, g_y_0_xyyyy_zzzzzz, g_y_0_yyyy_xxxxxx, g_y_0_yyyy_xxxxxxx, g_y_0_yyyy_xxxxxxy, g_y_0_yyyy_xxxxxxz, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxxyy, g_y_0_yyyy_xxxxxyz, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxxzz, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyyy, g_y_0_yyyy_xxxxyyz, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxyzz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxxzzz, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyyy, g_y_0_yyyy_xxxyyyz, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyyzz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxyzzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxxzzzz, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyyy, g_y_0_yyyy_xxyyyyz, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyyzz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyyzzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxyzzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xxzzzzz, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyyy, g_y_0_yyyy_xyyyyyz, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyyzz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyyzzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyyzzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xyzzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_xzzzzzz, g_y_0_yyyy_yyyyyy, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyy_xxxxxx[k] = -g_y_0_yyyy_xxxxxx[k] * ab_x + g_y_0_yyyy_xxxxxxx[k];

                g_y_0_xyyyy_xxxxxy[k] = -g_y_0_yyyy_xxxxxy[k] * ab_x + g_y_0_yyyy_xxxxxxy[k];

                g_y_0_xyyyy_xxxxxz[k] = -g_y_0_yyyy_xxxxxz[k] * ab_x + g_y_0_yyyy_xxxxxxz[k];

                g_y_0_xyyyy_xxxxyy[k] = -g_y_0_yyyy_xxxxyy[k] * ab_x + g_y_0_yyyy_xxxxxyy[k];

                g_y_0_xyyyy_xxxxyz[k] = -g_y_0_yyyy_xxxxyz[k] * ab_x + g_y_0_yyyy_xxxxxyz[k];

                g_y_0_xyyyy_xxxxzz[k] = -g_y_0_yyyy_xxxxzz[k] * ab_x + g_y_0_yyyy_xxxxxzz[k];

                g_y_0_xyyyy_xxxyyy[k] = -g_y_0_yyyy_xxxyyy[k] * ab_x + g_y_0_yyyy_xxxxyyy[k];

                g_y_0_xyyyy_xxxyyz[k] = -g_y_0_yyyy_xxxyyz[k] * ab_x + g_y_0_yyyy_xxxxyyz[k];

                g_y_0_xyyyy_xxxyzz[k] = -g_y_0_yyyy_xxxyzz[k] * ab_x + g_y_0_yyyy_xxxxyzz[k];

                g_y_0_xyyyy_xxxzzz[k] = -g_y_0_yyyy_xxxzzz[k] * ab_x + g_y_0_yyyy_xxxxzzz[k];

                g_y_0_xyyyy_xxyyyy[k] = -g_y_0_yyyy_xxyyyy[k] * ab_x + g_y_0_yyyy_xxxyyyy[k];

                g_y_0_xyyyy_xxyyyz[k] = -g_y_0_yyyy_xxyyyz[k] * ab_x + g_y_0_yyyy_xxxyyyz[k];

                g_y_0_xyyyy_xxyyzz[k] = -g_y_0_yyyy_xxyyzz[k] * ab_x + g_y_0_yyyy_xxxyyzz[k];

                g_y_0_xyyyy_xxyzzz[k] = -g_y_0_yyyy_xxyzzz[k] * ab_x + g_y_0_yyyy_xxxyzzz[k];

                g_y_0_xyyyy_xxzzzz[k] = -g_y_0_yyyy_xxzzzz[k] * ab_x + g_y_0_yyyy_xxxzzzz[k];

                g_y_0_xyyyy_xyyyyy[k] = -g_y_0_yyyy_xyyyyy[k] * ab_x + g_y_0_yyyy_xxyyyyy[k];

                g_y_0_xyyyy_xyyyyz[k] = -g_y_0_yyyy_xyyyyz[k] * ab_x + g_y_0_yyyy_xxyyyyz[k];

                g_y_0_xyyyy_xyyyzz[k] = -g_y_0_yyyy_xyyyzz[k] * ab_x + g_y_0_yyyy_xxyyyzz[k];

                g_y_0_xyyyy_xyyzzz[k] = -g_y_0_yyyy_xyyzzz[k] * ab_x + g_y_0_yyyy_xxyyzzz[k];

                g_y_0_xyyyy_xyzzzz[k] = -g_y_0_yyyy_xyzzzz[k] * ab_x + g_y_0_yyyy_xxyzzzz[k];

                g_y_0_xyyyy_xzzzzz[k] = -g_y_0_yyyy_xzzzzz[k] * ab_x + g_y_0_yyyy_xxzzzzz[k];

                g_y_0_xyyyy_yyyyyy[k] = -g_y_0_yyyy_yyyyyy[k] * ab_x + g_y_0_yyyy_xyyyyyy[k];

                g_y_0_xyyyy_yyyyyz[k] = -g_y_0_yyyy_yyyyyz[k] * ab_x + g_y_0_yyyy_xyyyyyz[k];

                g_y_0_xyyyy_yyyyzz[k] = -g_y_0_yyyy_yyyyzz[k] * ab_x + g_y_0_yyyy_xyyyyzz[k];

                g_y_0_xyyyy_yyyzzz[k] = -g_y_0_yyyy_yyyzzz[k] * ab_x + g_y_0_yyyy_xyyyzzz[k];

                g_y_0_xyyyy_yyzzzz[k] = -g_y_0_yyyy_yyzzzz[k] * ab_x + g_y_0_yyyy_xyyzzzz[k];

                g_y_0_xyyyy_yzzzzz[k] = -g_y_0_yyyy_yzzzzz[k] * ab_x + g_y_0_yyyy_xyzzzzz[k];

                g_y_0_xyyyy_zzzzzz[k] = -g_y_0_yyyy_zzzzzz[k] * ab_x + g_y_0_yyyy_xzzzzzz[k];
            }

            /// Set up 896-924 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 896 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 897 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 898 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 899 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 900 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 901 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 902 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 903 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 904 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 905 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 906 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 907 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 908 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 909 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 910 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 911 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 912 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 913 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 914 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 915 * ccomps * dcomps);

            auto g_y_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 916 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 917 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 918 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 919 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 920 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 921 * ccomps * dcomps);

            auto g_y_0_xyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 922 * ccomps * dcomps);

            auto g_y_0_xyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyz_xxxxxx, g_y_0_xyyyz_xxxxxy, g_y_0_xyyyz_xxxxxz, g_y_0_xyyyz_xxxxyy, g_y_0_xyyyz_xxxxyz, g_y_0_xyyyz_xxxxzz, g_y_0_xyyyz_xxxyyy, g_y_0_xyyyz_xxxyyz, g_y_0_xyyyz_xxxyzz, g_y_0_xyyyz_xxxzzz, g_y_0_xyyyz_xxyyyy, g_y_0_xyyyz_xxyyyz, g_y_0_xyyyz_xxyyzz, g_y_0_xyyyz_xxyzzz, g_y_0_xyyyz_xxzzzz, g_y_0_xyyyz_xyyyyy, g_y_0_xyyyz_xyyyyz, g_y_0_xyyyz_xyyyzz, g_y_0_xyyyz_xyyzzz, g_y_0_xyyyz_xyzzzz, g_y_0_xyyyz_xzzzzz, g_y_0_xyyyz_yyyyyy, g_y_0_xyyyz_yyyyyz, g_y_0_xyyyz_yyyyzz, g_y_0_xyyyz_yyyzzz, g_y_0_xyyyz_yyzzzz, g_y_0_xyyyz_yzzzzz, g_y_0_xyyyz_zzzzzz, g_y_0_yyyz_xxxxxx, g_y_0_yyyz_xxxxxxx, g_y_0_yyyz_xxxxxxy, g_y_0_yyyz_xxxxxxz, g_y_0_yyyz_xxxxxy, g_y_0_yyyz_xxxxxyy, g_y_0_yyyz_xxxxxyz, g_y_0_yyyz_xxxxxz, g_y_0_yyyz_xxxxxzz, g_y_0_yyyz_xxxxyy, g_y_0_yyyz_xxxxyyy, g_y_0_yyyz_xxxxyyz, g_y_0_yyyz_xxxxyz, g_y_0_yyyz_xxxxyzz, g_y_0_yyyz_xxxxzz, g_y_0_yyyz_xxxxzzz, g_y_0_yyyz_xxxyyy, g_y_0_yyyz_xxxyyyy, g_y_0_yyyz_xxxyyyz, g_y_0_yyyz_xxxyyz, g_y_0_yyyz_xxxyyzz, g_y_0_yyyz_xxxyzz, g_y_0_yyyz_xxxyzzz, g_y_0_yyyz_xxxzzz, g_y_0_yyyz_xxxzzzz, g_y_0_yyyz_xxyyyy, g_y_0_yyyz_xxyyyyy, g_y_0_yyyz_xxyyyyz, g_y_0_yyyz_xxyyyz, g_y_0_yyyz_xxyyyzz, g_y_0_yyyz_xxyyzz, g_y_0_yyyz_xxyyzzz, g_y_0_yyyz_xxyzzz, g_y_0_yyyz_xxyzzzz, g_y_0_yyyz_xxzzzz, g_y_0_yyyz_xxzzzzz, g_y_0_yyyz_xyyyyy, g_y_0_yyyz_xyyyyyy, g_y_0_yyyz_xyyyyyz, g_y_0_yyyz_xyyyyz, g_y_0_yyyz_xyyyyzz, g_y_0_yyyz_xyyyzz, g_y_0_yyyz_xyyyzzz, g_y_0_yyyz_xyyzzz, g_y_0_yyyz_xyyzzzz, g_y_0_yyyz_xyzzzz, g_y_0_yyyz_xyzzzzz, g_y_0_yyyz_xzzzzz, g_y_0_yyyz_xzzzzzz, g_y_0_yyyz_yyyyyy, g_y_0_yyyz_yyyyyz, g_y_0_yyyz_yyyyzz, g_y_0_yyyz_yyyzzz, g_y_0_yyyz_yyzzzz, g_y_0_yyyz_yzzzzz, g_y_0_yyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyz_xxxxxx[k] = -g_y_0_yyyz_xxxxxx[k] * ab_x + g_y_0_yyyz_xxxxxxx[k];

                g_y_0_xyyyz_xxxxxy[k] = -g_y_0_yyyz_xxxxxy[k] * ab_x + g_y_0_yyyz_xxxxxxy[k];

                g_y_0_xyyyz_xxxxxz[k] = -g_y_0_yyyz_xxxxxz[k] * ab_x + g_y_0_yyyz_xxxxxxz[k];

                g_y_0_xyyyz_xxxxyy[k] = -g_y_0_yyyz_xxxxyy[k] * ab_x + g_y_0_yyyz_xxxxxyy[k];

                g_y_0_xyyyz_xxxxyz[k] = -g_y_0_yyyz_xxxxyz[k] * ab_x + g_y_0_yyyz_xxxxxyz[k];

                g_y_0_xyyyz_xxxxzz[k] = -g_y_0_yyyz_xxxxzz[k] * ab_x + g_y_0_yyyz_xxxxxzz[k];

                g_y_0_xyyyz_xxxyyy[k] = -g_y_0_yyyz_xxxyyy[k] * ab_x + g_y_0_yyyz_xxxxyyy[k];

                g_y_0_xyyyz_xxxyyz[k] = -g_y_0_yyyz_xxxyyz[k] * ab_x + g_y_0_yyyz_xxxxyyz[k];

                g_y_0_xyyyz_xxxyzz[k] = -g_y_0_yyyz_xxxyzz[k] * ab_x + g_y_0_yyyz_xxxxyzz[k];

                g_y_0_xyyyz_xxxzzz[k] = -g_y_0_yyyz_xxxzzz[k] * ab_x + g_y_0_yyyz_xxxxzzz[k];

                g_y_0_xyyyz_xxyyyy[k] = -g_y_0_yyyz_xxyyyy[k] * ab_x + g_y_0_yyyz_xxxyyyy[k];

                g_y_0_xyyyz_xxyyyz[k] = -g_y_0_yyyz_xxyyyz[k] * ab_x + g_y_0_yyyz_xxxyyyz[k];

                g_y_0_xyyyz_xxyyzz[k] = -g_y_0_yyyz_xxyyzz[k] * ab_x + g_y_0_yyyz_xxxyyzz[k];

                g_y_0_xyyyz_xxyzzz[k] = -g_y_0_yyyz_xxyzzz[k] * ab_x + g_y_0_yyyz_xxxyzzz[k];

                g_y_0_xyyyz_xxzzzz[k] = -g_y_0_yyyz_xxzzzz[k] * ab_x + g_y_0_yyyz_xxxzzzz[k];

                g_y_0_xyyyz_xyyyyy[k] = -g_y_0_yyyz_xyyyyy[k] * ab_x + g_y_0_yyyz_xxyyyyy[k];

                g_y_0_xyyyz_xyyyyz[k] = -g_y_0_yyyz_xyyyyz[k] * ab_x + g_y_0_yyyz_xxyyyyz[k];

                g_y_0_xyyyz_xyyyzz[k] = -g_y_0_yyyz_xyyyzz[k] * ab_x + g_y_0_yyyz_xxyyyzz[k];

                g_y_0_xyyyz_xyyzzz[k] = -g_y_0_yyyz_xyyzzz[k] * ab_x + g_y_0_yyyz_xxyyzzz[k];

                g_y_0_xyyyz_xyzzzz[k] = -g_y_0_yyyz_xyzzzz[k] * ab_x + g_y_0_yyyz_xxyzzzz[k];

                g_y_0_xyyyz_xzzzzz[k] = -g_y_0_yyyz_xzzzzz[k] * ab_x + g_y_0_yyyz_xxzzzzz[k];

                g_y_0_xyyyz_yyyyyy[k] = -g_y_0_yyyz_yyyyyy[k] * ab_x + g_y_0_yyyz_xyyyyyy[k];

                g_y_0_xyyyz_yyyyyz[k] = -g_y_0_yyyz_yyyyyz[k] * ab_x + g_y_0_yyyz_xyyyyyz[k];

                g_y_0_xyyyz_yyyyzz[k] = -g_y_0_yyyz_yyyyzz[k] * ab_x + g_y_0_yyyz_xyyyyzz[k];

                g_y_0_xyyyz_yyyzzz[k] = -g_y_0_yyyz_yyyzzz[k] * ab_x + g_y_0_yyyz_xyyyzzz[k];

                g_y_0_xyyyz_yyzzzz[k] = -g_y_0_yyyz_yyzzzz[k] * ab_x + g_y_0_yyyz_xyyzzzz[k];

                g_y_0_xyyyz_yzzzzz[k] = -g_y_0_yyyz_yzzzzz[k] * ab_x + g_y_0_yyyz_xyzzzzz[k];

                g_y_0_xyyyz_zzzzzz[k] = -g_y_0_yyyz_zzzzzz[k] * ab_x + g_y_0_yyyz_xzzzzzz[k];
            }

            /// Set up 924-952 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 924 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 925 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 926 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 927 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 928 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 929 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 930 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 931 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 932 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 933 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 934 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 935 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 936 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 937 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 938 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 939 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 940 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 941 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 942 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 943 * ccomps * dcomps);

            auto g_y_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 944 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 945 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 946 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 947 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 948 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 949 * ccomps * dcomps);

            auto g_y_0_xyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 950 * ccomps * dcomps);

            auto g_y_0_xyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 951 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyzz_xxxxxx, g_y_0_xyyzz_xxxxxy, g_y_0_xyyzz_xxxxxz, g_y_0_xyyzz_xxxxyy, g_y_0_xyyzz_xxxxyz, g_y_0_xyyzz_xxxxzz, g_y_0_xyyzz_xxxyyy, g_y_0_xyyzz_xxxyyz, g_y_0_xyyzz_xxxyzz, g_y_0_xyyzz_xxxzzz, g_y_0_xyyzz_xxyyyy, g_y_0_xyyzz_xxyyyz, g_y_0_xyyzz_xxyyzz, g_y_0_xyyzz_xxyzzz, g_y_0_xyyzz_xxzzzz, g_y_0_xyyzz_xyyyyy, g_y_0_xyyzz_xyyyyz, g_y_0_xyyzz_xyyyzz, g_y_0_xyyzz_xyyzzz, g_y_0_xyyzz_xyzzzz, g_y_0_xyyzz_xzzzzz, g_y_0_xyyzz_yyyyyy, g_y_0_xyyzz_yyyyyz, g_y_0_xyyzz_yyyyzz, g_y_0_xyyzz_yyyzzz, g_y_0_xyyzz_yyzzzz, g_y_0_xyyzz_yzzzzz, g_y_0_xyyzz_zzzzzz, g_y_0_yyzz_xxxxxx, g_y_0_yyzz_xxxxxxx, g_y_0_yyzz_xxxxxxy, g_y_0_yyzz_xxxxxxz, g_y_0_yyzz_xxxxxy, g_y_0_yyzz_xxxxxyy, g_y_0_yyzz_xxxxxyz, g_y_0_yyzz_xxxxxz, g_y_0_yyzz_xxxxxzz, g_y_0_yyzz_xxxxyy, g_y_0_yyzz_xxxxyyy, g_y_0_yyzz_xxxxyyz, g_y_0_yyzz_xxxxyz, g_y_0_yyzz_xxxxyzz, g_y_0_yyzz_xxxxzz, g_y_0_yyzz_xxxxzzz, g_y_0_yyzz_xxxyyy, g_y_0_yyzz_xxxyyyy, g_y_0_yyzz_xxxyyyz, g_y_0_yyzz_xxxyyz, g_y_0_yyzz_xxxyyzz, g_y_0_yyzz_xxxyzz, g_y_0_yyzz_xxxyzzz, g_y_0_yyzz_xxxzzz, g_y_0_yyzz_xxxzzzz, g_y_0_yyzz_xxyyyy, g_y_0_yyzz_xxyyyyy, g_y_0_yyzz_xxyyyyz, g_y_0_yyzz_xxyyyz, g_y_0_yyzz_xxyyyzz, g_y_0_yyzz_xxyyzz, g_y_0_yyzz_xxyyzzz, g_y_0_yyzz_xxyzzz, g_y_0_yyzz_xxyzzzz, g_y_0_yyzz_xxzzzz, g_y_0_yyzz_xxzzzzz, g_y_0_yyzz_xyyyyy, g_y_0_yyzz_xyyyyyy, g_y_0_yyzz_xyyyyyz, g_y_0_yyzz_xyyyyz, g_y_0_yyzz_xyyyyzz, g_y_0_yyzz_xyyyzz, g_y_0_yyzz_xyyyzzz, g_y_0_yyzz_xyyzzz, g_y_0_yyzz_xyyzzzz, g_y_0_yyzz_xyzzzz, g_y_0_yyzz_xyzzzzz, g_y_0_yyzz_xzzzzz, g_y_0_yyzz_xzzzzzz, g_y_0_yyzz_yyyyyy, g_y_0_yyzz_yyyyyz, g_y_0_yyzz_yyyyzz, g_y_0_yyzz_yyyzzz, g_y_0_yyzz_yyzzzz, g_y_0_yyzz_yzzzzz, g_y_0_yyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzz_xxxxxx[k] = -g_y_0_yyzz_xxxxxx[k] * ab_x + g_y_0_yyzz_xxxxxxx[k];

                g_y_0_xyyzz_xxxxxy[k] = -g_y_0_yyzz_xxxxxy[k] * ab_x + g_y_0_yyzz_xxxxxxy[k];

                g_y_0_xyyzz_xxxxxz[k] = -g_y_0_yyzz_xxxxxz[k] * ab_x + g_y_0_yyzz_xxxxxxz[k];

                g_y_0_xyyzz_xxxxyy[k] = -g_y_0_yyzz_xxxxyy[k] * ab_x + g_y_0_yyzz_xxxxxyy[k];

                g_y_0_xyyzz_xxxxyz[k] = -g_y_0_yyzz_xxxxyz[k] * ab_x + g_y_0_yyzz_xxxxxyz[k];

                g_y_0_xyyzz_xxxxzz[k] = -g_y_0_yyzz_xxxxzz[k] * ab_x + g_y_0_yyzz_xxxxxzz[k];

                g_y_0_xyyzz_xxxyyy[k] = -g_y_0_yyzz_xxxyyy[k] * ab_x + g_y_0_yyzz_xxxxyyy[k];

                g_y_0_xyyzz_xxxyyz[k] = -g_y_0_yyzz_xxxyyz[k] * ab_x + g_y_0_yyzz_xxxxyyz[k];

                g_y_0_xyyzz_xxxyzz[k] = -g_y_0_yyzz_xxxyzz[k] * ab_x + g_y_0_yyzz_xxxxyzz[k];

                g_y_0_xyyzz_xxxzzz[k] = -g_y_0_yyzz_xxxzzz[k] * ab_x + g_y_0_yyzz_xxxxzzz[k];

                g_y_0_xyyzz_xxyyyy[k] = -g_y_0_yyzz_xxyyyy[k] * ab_x + g_y_0_yyzz_xxxyyyy[k];

                g_y_0_xyyzz_xxyyyz[k] = -g_y_0_yyzz_xxyyyz[k] * ab_x + g_y_0_yyzz_xxxyyyz[k];

                g_y_0_xyyzz_xxyyzz[k] = -g_y_0_yyzz_xxyyzz[k] * ab_x + g_y_0_yyzz_xxxyyzz[k];

                g_y_0_xyyzz_xxyzzz[k] = -g_y_0_yyzz_xxyzzz[k] * ab_x + g_y_0_yyzz_xxxyzzz[k];

                g_y_0_xyyzz_xxzzzz[k] = -g_y_0_yyzz_xxzzzz[k] * ab_x + g_y_0_yyzz_xxxzzzz[k];

                g_y_0_xyyzz_xyyyyy[k] = -g_y_0_yyzz_xyyyyy[k] * ab_x + g_y_0_yyzz_xxyyyyy[k];

                g_y_0_xyyzz_xyyyyz[k] = -g_y_0_yyzz_xyyyyz[k] * ab_x + g_y_0_yyzz_xxyyyyz[k];

                g_y_0_xyyzz_xyyyzz[k] = -g_y_0_yyzz_xyyyzz[k] * ab_x + g_y_0_yyzz_xxyyyzz[k];

                g_y_0_xyyzz_xyyzzz[k] = -g_y_0_yyzz_xyyzzz[k] * ab_x + g_y_0_yyzz_xxyyzzz[k];

                g_y_0_xyyzz_xyzzzz[k] = -g_y_0_yyzz_xyzzzz[k] * ab_x + g_y_0_yyzz_xxyzzzz[k];

                g_y_0_xyyzz_xzzzzz[k] = -g_y_0_yyzz_xzzzzz[k] * ab_x + g_y_0_yyzz_xxzzzzz[k];

                g_y_0_xyyzz_yyyyyy[k] = -g_y_0_yyzz_yyyyyy[k] * ab_x + g_y_0_yyzz_xyyyyyy[k];

                g_y_0_xyyzz_yyyyyz[k] = -g_y_0_yyzz_yyyyyz[k] * ab_x + g_y_0_yyzz_xyyyyyz[k];

                g_y_0_xyyzz_yyyyzz[k] = -g_y_0_yyzz_yyyyzz[k] * ab_x + g_y_0_yyzz_xyyyyzz[k];

                g_y_0_xyyzz_yyyzzz[k] = -g_y_0_yyzz_yyyzzz[k] * ab_x + g_y_0_yyzz_xyyyzzz[k];

                g_y_0_xyyzz_yyzzzz[k] = -g_y_0_yyzz_yyzzzz[k] * ab_x + g_y_0_yyzz_xyyzzzz[k];

                g_y_0_xyyzz_yzzzzz[k] = -g_y_0_yyzz_yzzzzz[k] * ab_x + g_y_0_yyzz_xyzzzzz[k];

                g_y_0_xyyzz_zzzzzz[k] = -g_y_0_yyzz_zzzzzz[k] * ab_x + g_y_0_yyzz_xzzzzzz[k];
            }

            /// Set up 952-980 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 952 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 953 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 954 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 955 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 956 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 957 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 958 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 959 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 960 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 961 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 962 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 963 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 964 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 965 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 966 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 967 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 968 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 969 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 970 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 971 * ccomps * dcomps);

            auto g_y_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 972 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 973 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 974 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 975 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 976 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 977 * ccomps * dcomps);

            auto g_y_0_xyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 978 * ccomps * dcomps);

            auto g_y_0_xyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 979 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzzz_xxxxxx, g_y_0_xyzzz_xxxxxy, g_y_0_xyzzz_xxxxxz, g_y_0_xyzzz_xxxxyy, g_y_0_xyzzz_xxxxyz, g_y_0_xyzzz_xxxxzz, g_y_0_xyzzz_xxxyyy, g_y_0_xyzzz_xxxyyz, g_y_0_xyzzz_xxxyzz, g_y_0_xyzzz_xxxzzz, g_y_0_xyzzz_xxyyyy, g_y_0_xyzzz_xxyyyz, g_y_0_xyzzz_xxyyzz, g_y_0_xyzzz_xxyzzz, g_y_0_xyzzz_xxzzzz, g_y_0_xyzzz_xyyyyy, g_y_0_xyzzz_xyyyyz, g_y_0_xyzzz_xyyyzz, g_y_0_xyzzz_xyyzzz, g_y_0_xyzzz_xyzzzz, g_y_0_xyzzz_xzzzzz, g_y_0_xyzzz_yyyyyy, g_y_0_xyzzz_yyyyyz, g_y_0_xyzzz_yyyyzz, g_y_0_xyzzz_yyyzzz, g_y_0_xyzzz_yyzzzz, g_y_0_xyzzz_yzzzzz, g_y_0_xyzzz_zzzzzz, g_y_0_yzzz_xxxxxx, g_y_0_yzzz_xxxxxxx, g_y_0_yzzz_xxxxxxy, g_y_0_yzzz_xxxxxxz, g_y_0_yzzz_xxxxxy, g_y_0_yzzz_xxxxxyy, g_y_0_yzzz_xxxxxyz, g_y_0_yzzz_xxxxxz, g_y_0_yzzz_xxxxxzz, g_y_0_yzzz_xxxxyy, g_y_0_yzzz_xxxxyyy, g_y_0_yzzz_xxxxyyz, g_y_0_yzzz_xxxxyz, g_y_0_yzzz_xxxxyzz, g_y_0_yzzz_xxxxzz, g_y_0_yzzz_xxxxzzz, g_y_0_yzzz_xxxyyy, g_y_0_yzzz_xxxyyyy, g_y_0_yzzz_xxxyyyz, g_y_0_yzzz_xxxyyz, g_y_0_yzzz_xxxyyzz, g_y_0_yzzz_xxxyzz, g_y_0_yzzz_xxxyzzz, g_y_0_yzzz_xxxzzz, g_y_0_yzzz_xxxzzzz, g_y_0_yzzz_xxyyyy, g_y_0_yzzz_xxyyyyy, g_y_0_yzzz_xxyyyyz, g_y_0_yzzz_xxyyyz, g_y_0_yzzz_xxyyyzz, g_y_0_yzzz_xxyyzz, g_y_0_yzzz_xxyyzzz, g_y_0_yzzz_xxyzzz, g_y_0_yzzz_xxyzzzz, g_y_0_yzzz_xxzzzz, g_y_0_yzzz_xxzzzzz, g_y_0_yzzz_xyyyyy, g_y_0_yzzz_xyyyyyy, g_y_0_yzzz_xyyyyyz, g_y_0_yzzz_xyyyyz, g_y_0_yzzz_xyyyyzz, g_y_0_yzzz_xyyyzz, g_y_0_yzzz_xyyyzzz, g_y_0_yzzz_xyyzzz, g_y_0_yzzz_xyyzzzz, g_y_0_yzzz_xyzzzz, g_y_0_yzzz_xyzzzzz, g_y_0_yzzz_xzzzzz, g_y_0_yzzz_xzzzzzz, g_y_0_yzzz_yyyyyy, g_y_0_yzzz_yyyyyz, g_y_0_yzzz_yyyyzz, g_y_0_yzzz_yyyzzz, g_y_0_yzzz_yyzzzz, g_y_0_yzzz_yzzzzz, g_y_0_yzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzz_xxxxxx[k] = -g_y_0_yzzz_xxxxxx[k] * ab_x + g_y_0_yzzz_xxxxxxx[k];

                g_y_0_xyzzz_xxxxxy[k] = -g_y_0_yzzz_xxxxxy[k] * ab_x + g_y_0_yzzz_xxxxxxy[k];

                g_y_0_xyzzz_xxxxxz[k] = -g_y_0_yzzz_xxxxxz[k] * ab_x + g_y_0_yzzz_xxxxxxz[k];

                g_y_0_xyzzz_xxxxyy[k] = -g_y_0_yzzz_xxxxyy[k] * ab_x + g_y_0_yzzz_xxxxxyy[k];

                g_y_0_xyzzz_xxxxyz[k] = -g_y_0_yzzz_xxxxyz[k] * ab_x + g_y_0_yzzz_xxxxxyz[k];

                g_y_0_xyzzz_xxxxzz[k] = -g_y_0_yzzz_xxxxzz[k] * ab_x + g_y_0_yzzz_xxxxxzz[k];

                g_y_0_xyzzz_xxxyyy[k] = -g_y_0_yzzz_xxxyyy[k] * ab_x + g_y_0_yzzz_xxxxyyy[k];

                g_y_0_xyzzz_xxxyyz[k] = -g_y_0_yzzz_xxxyyz[k] * ab_x + g_y_0_yzzz_xxxxyyz[k];

                g_y_0_xyzzz_xxxyzz[k] = -g_y_0_yzzz_xxxyzz[k] * ab_x + g_y_0_yzzz_xxxxyzz[k];

                g_y_0_xyzzz_xxxzzz[k] = -g_y_0_yzzz_xxxzzz[k] * ab_x + g_y_0_yzzz_xxxxzzz[k];

                g_y_0_xyzzz_xxyyyy[k] = -g_y_0_yzzz_xxyyyy[k] * ab_x + g_y_0_yzzz_xxxyyyy[k];

                g_y_0_xyzzz_xxyyyz[k] = -g_y_0_yzzz_xxyyyz[k] * ab_x + g_y_0_yzzz_xxxyyyz[k];

                g_y_0_xyzzz_xxyyzz[k] = -g_y_0_yzzz_xxyyzz[k] * ab_x + g_y_0_yzzz_xxxyyzz[k];

                g_y_0_xyzzz_xxyzzz[k] = -g_y_0_yzzz_xxyzzz[k] * ab_x + g_y_0_yzzz_xxxyzzz[k];

                g_y_0_xyzzz_xxzzzz[k] = -g_y_0_yzzz_xxzzzz[k] * ab_x + g_y_0_yzzz_xxxzzzz[k];

                g_y_0_xyzzz_xyyyyy[k] = -g_y_0_yzzz_xyyyyy[k] * ab_x + g_y_0_yzzz_xxyyyyy[k];

                g_y_0_xyzzz_xyyyyz[k] = -g_y_0_yzzz_xyyyyz[k] * ab_x + g_y_0_yzzz_xxyyyyz[k];

                g_y_0_xyzzz_xyyyzz[k] = -g_y_0_yzzz_xyyyzz[k] * ab_x + g_y_0_yzzz_xxyyyzz[k];

                g_y_0_xyzzz_xyyzzz[k] = -g_y_0_yzzz_xyyzzz[k] * ab_x + g_y_0_yzzz_xxyyzzz[k];

                g_y_0_xyzzz_xyzzzz[k] = -g_y_0_yzzz_xyzzzz[k] * ab_x + g_y_0_yzzz_xxyzzzz[k];

                g_y_0_xyzzz_xzzzzz[k] = -g_y_0_yzzz_xzzzzz[k] * ab_x + g_y_0_yzzz_xxzzzzz[k];

                g_y_0_xyzzz_yyyyyy[k] = -g_y_0_yzzz_yyyyyy[k] * ab_x + g_y_0_yzzz_xyyyyyy[k];

                g_y_0_xyzzz_yyyyyz[k] = -g_y_0_yzzz_yyyyyz[k] * ab_x + g_y_0_yzzz_xyyyyyz[k];

                g_y_0_xyzzz_yyyyzz[k] = -g_y_0_yzzz_yyyyzz[k] * ab_x + g_y_0_yzzz_xyyyyzz[k];

                g_y_0_xyzzz_yyyzzz[k] = -g_y_0_yzzz_yyyzzz[k] * ab_x + g_y_0_yzzz_xyyyzzz[k];

                g_y_0_xyzzz_yyzzzz[k] = -g_y_0_yzzz_yyzzzz[k] * ab_x + g_y_0_yzzz_xyyzzzz[k];

                g_y_0_xyzzz_yzzzzz[k] = -g_y_0_yzzz_yzzzzz[k] * ab_x + g_y_0_yzzz_xyzzzzz[k];

                g_y_0_xyzzz_zzzzzz[k] = -g_y_0_yzzz_zzzzzz[k] * ab_x + g_y_0_yzzz_xzzzzzz[k];
            }

            /// Set up 980-1008 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 980 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 981 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 982 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 983 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 984 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 985 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 986 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 987 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 988 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 989 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 990 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 991 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 992 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 993 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 994 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 995 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 996 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 997 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 998 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 999 * ccomps * dcomps);

            auto g_y_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1000 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1001 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1002 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1003 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1004 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1005 * ccomps * dcomps);

            auto g_y_0_xzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1006 * ccomps * dcomps);

            auto g_y_0_xzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzzz_xxxxxx, g_y_0_xzzzz_xxxxxy, g_y_0_xzzzz_xxxxxz, g_y_0_xzzzz_xxxxyy, g_y_0_xzzzz_xxxxyz, g_y_0_xzzzz_xxxxzz, g_y_0_xzzzz_xxxyyy, g_y_0_xzzzz_xxxyyz, g_y_0_xzzzz_xxxyzz, g_y_0_xzzzz_xxxzzz, g_y_0_xzzzz_xxyyyy, g_y_0_xzzzz_xxyyyz, g_y_0_xzzzz_xxyyzz, g_y_0_xzzzz_xxyzzz, g_y_0_xzzzz_xxzzzz, g_y_0_xzzzz_xyyyyy, g_y_0_xzzzz_xyyyyz, g_y_0_xzzzz_xyyyzz, g_y_0_xzzzz_xyyzzz, g_y_0_xzzzz_xyzzzz, g_y_0_xzzzz_xzzzzz, g_y_0_xzzzz_yyyyyy, g_y_0_xzzzz_yyyyyz, g_y_0_xzzzz_yyyyzz, g_y_0_xzzzz_yyyzzz, g_y_0_xzzzz_yyzzzz, g_y_0_xzzzz_yzzzzz, g_y_0_xzzzz_zzzzzz, g_y_0_zzzz_xxxxxx, g_y_0_zzzz_xxxxxxx, g_y_0_zzzz_xxxxxxy, g_y_0_zzzz_xxxxxxz, g_y_0_zzzz_xxxxxy, g_y_0_zzzz_xxxxxyy, g_y_0_zzzz_xxxxxyz, g_y_0_zzzz_xxxxxz, g_y_0_zzzz_xxxxxzz, g_y_0_zzzz_xxxxyy, g_y_0_zzzz_xxxxyyy, g_y_0_zzzz_xxxxyyz, g_y_0_zzzz_xxxxyz, g_y_0_zzzz_xxxxyzz, g_y_0_zzzz_xxxxzz, g_y_0_zzzz_xxxxzzz, g_y_0_zzzz_xxxyyy, g_y_0_zzzz_xxxyyyy, g_y_0_zzzz_xxxyyyz, g_y_0_zzzz_xxxyyz, g_y_0_zzzz_xxxyyzz, g_y_0_zzzz_xxxyzz, g_y_0_zzzz_xxxyzzz, g_y_0_zzzz_xxxzzz, g_y_0_zzzz_xxxzzzz, g_y_0_zzzz_xxyyyy, g_y_0_zzzz_xxyyyyy, g_y_0_zzzz_xxyyyyz, g_y_0_zzzz_xxyyyz, g_y_0_zzzz_xxyyyzz, g_y_0_zzzz_xxyyzz, g_y_0_zzzz_xxyyzzz, g_y_0_zzzz_xxyzzz, g_y_0_zzzz_xxyzzzz, g_y_0_zzzz_xxzzzz, g_y_0_zzzz_xxzzzzz, g_y_0_zzzz_xyyyyy, g_y_0_zzzz_xyyyyyy, g_y_0_zzzz_xyyyyyz, g_y_0_zzzz_xyyyyz, g_y_0_zzzz_xyyyyzz, g_y_0_zzzz_xyyyzz, g_y_0_zzzz_xyyyzzz, g_y_0_zzzz_xyyzzz, g_y_0_zzzz_xyyzzzz, g_y_0_zzzz_xyzzzz, g_y_0_zzzz_xyzzzzz, g_y_0_zzzz_xzzzzz, g_y_0_zzzz_xzzzzzz, g_y_0_zzzz_yyyyyy, g_y_0_zzzz_yyyyyz, g_y_0_zzzz_yyyyzz, g_y_0_zzzz_yyyzzz, g_y_0_zzzz_yyzzzz, g_y_0_zzzz_yzzzzz, g_y_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzz_xxxxxx[k] = -g_y_0_zzzz_xxxxxx[k] * ab_x + g_y_0_zzzz_xxxxxxx[k];

                g_y_0_xzzzz_xxxxxy[k] = -g_y_0_zzzz_xxxxxy[k] * ab_x + g_y_0_zzzz_xxxxxxy[k];

                g_y_0_xzzzz_xxxxxz[k] = -g_y_0_zzzz_xxxxxz[k] * ab_x + g_y_0_zzzz_xxxxxxz[k];

                g_y_0_xzzzz_xxxxyy[k] = -g_y_0_zzzz_xxxxyy[k] * ab_x + g_y_0_zzzz_xxxxxyy[k];

                g_y_0_xzzzz_xxxxyz[k] = -g_y_0_zzzz_xxxxyz[k] * ab_x + g_y_0_zzzz_xxxxxyz[k];

                g_y_0_xzzzz_xxxxzz[k] = -g_y_0_zzzz_xxxxzz[k] * ab_x + g_y_0_zzzz_xxxxxzz[k];

                g_y_0_xzzzz_xxxyyy[k] = -g_y_0_zzzz_xxxyyy[k] * ab_x + g_y_0_zzzz_xxxxyyy[k];

                g_y_0_xzzzz_xxxyyz[k] = -g_y_0_zzzz_xxxyyz[k] * ab_x + g_y_0_zzzz_xxxxyyz[k];

                g_y_0_xzzzz_xxxyzz[k] = -g_y_0_zzzz_xxxyzz[k] * ab_x + g_y_0_zzzz_xxxxyzz[k];

                g_y_0_xzzzz_xxxzzz[k] = -g_y_0_zzzz_xxxzzz[k] * ab_x + g_y_0_zzzz_xxxxzzz[k];

                g_y_0_xzzzz_xxyyyy[k] = -g_y_0_zzzz_xxyyyy[k] * ab_x + g_y_0_zzzz_xxxyyyy[k];

                g_y_0_xzzzz_xxyyyz[k] = -g_y_0_zzzz_xxyyyz[k] * ab_x + g_y_0_zzzz_xxxyyyz[k];

                g_y_0_xzzzz_xxyyzz[k] = -g_y_0_zzzz_xxyyzz[k] * ab_x + g_y_0_zzzz_xxxyyzz[k];

                g_y_0_xzzzz_xxyzzz[k] = -g_y_0_zzzz_xxyzzz[k] * ab_x + g_y_0_zzzz_xxxyzzz[k];

                g_y_0_xzzzz_xxzzzz[k] = -g_y_0_zzzz_xxzzzz[k] * ab_x + g_y_0_zzzz_xxxzzzz[k];

                g_y_0_xzzzz_xyyyyy[k] = -g_y_0_zzzz_xyyyyy[k] * ab_x + g_y_0_zzzz_xxyyyyy[k];

                g_y_0_xzzzz_xyyyyz[k] = -g_y_0_zzzz_xyyyyz[k] * ab_x + g_y_0_zzzz_xxyyyyz[k];

                g_y_0_xzzzz_xyyyzz[k] = -g_y_0_zzzz_xyyyzz[k] * ab_x + g_y_0_zzzz_xxyyyzz[k];

                g_y_0_xzzzz_xyyzzz[k] = -g_y_0_zzzz_xyyzzz[k] * ab_x + g_y_0_zzzz_xxyyzzz[k];

                g_y_0_xzzzz_xyzzzz[k] = -g_y_0_zzzz_xyzzzz[k] * ab_x + g_y_0_zzzz_xxyzzzz[k];

                g_y_0_xzzzz_xzzzzz[k] = -g_y_0_zzzz_xzzzzz[k] * ab_x + g_y_0_zzzz_xxzzzzz[k];

                g_y_0_xzzzz_yyyyyy[k] = -g_y_0_zzzz_yyyyyy[k] * ab_x + g_y_0_zzzz_xyyyyyy[k];

                g_y_0_xzzzz_yyyyyz[k] = -g_y_0_zzzz_yyyyyz[k] * ab_x + g_y_0_zzzz_xyyyyyz[k];

                g_y_0_xzzzz_yyyyzz[k] = -g_y_0_zzzz_yyyyzz[k] * ab_x + g_y_0_zzzz_xyyyyzz[k];

                g_y_0_xzzzz_yyyzzz[k] = -g_y_0_zzzz_yyyzzz[k] * ab_x + g_y_0_zzzz_xyyyzzz[k];

                g_y_0_xzzzz_yyzzzz[k] = -g_y_0_zzzz_yyzzzz[k] * ab_x + g_y_0_zzzz_xyyzzzz[k];

                g_y_0_xzzzz_yzzzzz[k] = -g_y_0_zzzz_yzzzzz[k] * ab_x + g_y_0_zzzz_xyzzzzz[k];

                g_y_0_xzzzz_zzzzzz[k] = -g_y_0_zzzz_zzzzzz[k] * ab_x + g_y_0_zzzz_xzzzzzz[k];
            }

            /// Set up 1008-1036 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1008 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1009 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1010 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1011 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1012 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1013 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1014 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1015 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1016 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1017 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1018 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1019 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1020 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1021 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1022 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1023 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1024 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1025 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1026 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1027 * ccomps * dcomps);

            auto g_y_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1028 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1029 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1030 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1031 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1032 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1033 * ccomps * dcomps);

            auto g_y_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1034 * ccomps * dcomps);

            auto g_y_0_yyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1035 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyy_xxxxxx, g_y_0_yyyy_xxxxxxy, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxxyy, g_y_0_yyyy_xxxxxyz, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyyy, g_y_0_yyyy_xxxxyyz, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxyzz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyyy, g_y_0_yyyy_xxxyyyz, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyyzz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxyzzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyyy, g_y_0_yyyy_xxyyyyz, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyyzz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyyzzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxyzzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyyy, g_y_0_yyyy_xyyyyyz, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyyzz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyyzzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyyzzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xyzzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_yyyyyy, g_y_0_yyyy_yyyyyyy, g_y_0_yyyy_yyyyyyz, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyyzz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyyzzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyyzzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yyzzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_yzzzzzz, g_y_0_yyyy_zzzzzz, g_y_0_yyyyy_xxxxxx, g_y_0_yyyyy_xxxxxy, g_y_0_yyyyy_xxxxxz, g_y_0_yyyyy_xxxxyy, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxzz, g_y_0_yyyyy_xxxyyy, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxzzz, g_y_0_yyyyy_xxyyyy, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxzzzz, g_y_0_yyyyy_xyyyyy, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xzzzzz, g_y_0_yyyyy_yyyyyy, g_y_0_yyyyy_yyyyyz, g_y_0_yyyyy_yyyyzz, g_y_0_yyyyy_yyyzzz, g_y_0_yyyyy_yyzzzz, g_y_0_yyyyy_yzzzzz, g_y_0_yyyyy_zzzzzz, g_yyyy_xxxxxx, g_yyyy_xxxxxy, g_yyyy_xxxxxz, g_yyyy_xxxxyy, g_yyyy_xxxxyz, g_yyyy_xxxxzz, g_yyyy_xxxyyy, g_yyyy_xxxyyz, g_yyyy_xxxyzz, g_yyyy_xxxzzz, g_yyyy_xxyyyy, g_yyyy_xxyyyz, g_yyyy_xxyyzz, g_yyyy_xxyzzz, g_yyyy_xxzzzz, g_yyyy_xyyyyy, g_yyyy_xyyyyz, g_yyyy_xyyyzz, g_yyyy_xyyzzz, g_yyyy_xyzzzz, g_yyyy_xzzzzz, g_yyyy_yyyyyy, g_yyyy_yyyyyz, g_yyyy_yyyyzz, g_yyyy_yyyzzz, g_yyyy_yyzzzz, g_yyyy_yzzzzz, g_yyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyy_xxxxxx[k] = -g_yyyy_xxxxxx[k] - g_y_0_yyyy_xxxxxx[k] * ab_y + g_y_0_yyyy_xxxxxxy[k];

                g_y_0_yyyyy_xxxxxy[k] = -g_yyyy_xxxxxy[k] - g_y_0_yyyy_xxxxxy[k] * ab_y + g_y_0_yyyy_xxxxxyy[k];

                g_y_0_yyyyy_xxxxxz[k] = -g_yyyy_xxxxxz[k] - g_y_0_yyyy_xxxxxz[k] * ab_y + g_y_0_yyyy_xxxxxyz[k];

                g_y_0_yyyyy_xxxxyy[k] = -g_yyyy_xxxxyy[k] - g_y_0_yyyy_xxxxyy[k] * ab_y + g_y_0_yyyy_xxxxyyy[k];

                g_y_0_yyyyy_xxxxyz[k] = -g_yyyy_xxxxyz[k] - g_y_0_yyyy_xxxxyz[k] * ab_y + g_y_0_yyyy_xxxxyyz[k];

                g_y_0_yyyyy_xxxxzz[k] = -g_yyyy_xxxxzz[k] - g_y_0_yyyy_xxxxzz[k] * ab_y + g_y_0_yyyy_xxxxyzz[k];

                g_y_0_yyyyy_xxxyyy[k] = -g_yyyy_xxxyyy[k] - g_y_0_yyyy_xxxyyy[k] * ab_y + g_y_0_yyyy_xxxyyyy[k];

                g_y_0_yyyyy_xxxyyz[k] = -g_yyyy_xxxyyz[k] - g_y_0_yyyy_xxxyyz[k] * ab_y + g_y_0_yyyy_xxxyyyz[k];

                g_y_0_yyyyy_xxxyzz[k] = -g_yyyy_xxxyzz[k] - g_y_0_yyyy_xxxyzz[k] * ab_y + g_y_0_yyyy_xxxyyzz[k];

                g_y_0_yyyyy_xxxzzz[k] = -g_yyyy_xxxzzz[k] - g_y_0_yyyy_xxxzzz[k] * ab_y + g_y_0_yyyy_xxxyzzz[k];

                g_y_0_yyyyy_xxyyyy[k] = -g_yyyy_xxyyyy[k] - g_y_0_yyyy_xxyyyy[k] * ab_y + g_y_0_yyyy_xxyyyyy[k];

                g_y_0_yyyyy_xxyyyz[k] = -g_yyyy_xxyyyz[k] - g_y_0_yyyy_xxyyyz[k] * ab_y + g_y_0_yyyy_xxyyyyz[k];

                g_y_0_yyyyy_xxyyzz[k] = -g_yyyy_xxyyzz[k] - g_y_0_yyyy_xxyyzz[k] * ab_y + g_y_0_yyyy_xxyyyzz[k];

                g_y_0_yyyyy_xxyzzz[k] = -g_yyyy_xxyzzz[k] - g_y_0_yyyy_xxyzzz[k] * ab_y + g_y_0_yyyy_xxyyzzz[k];

                g_y_0_yyyyy_xxzzzz[k] = -g_yyyy_xxzzzz[k] - g_y_0_yyyy_xxzzzz[k] * ab_y + g_y_0_yyyy_xxyzzzz[k];

                g_y_0_yyyyy_xyyyyy[k] = -g_yyyy_xyyyyy[k] - g_y_0_yyyy_xyyyyy[k] * ab_y + g_y_0_yyyy_xyyyyyy[k];

                g_y_0_yyyyy_xyyyyz[k] = -g_yyyy_xyyyyz[k] - g_y_0_yyyy_xyyyyz[k] * ab_y + g_y_0_yyyy_xyyyyyz[k];

                g_y_0_yyyyy_xyyyzz[k] = -g_yyyy_xyyyzz[k] - g_y_0_yyyy_xyyyzz[k] * ab_y + g_y_0_yyyy_xyyyyzz[k];

                g_y_0_yyyyy_xyyzzz[k] = -g_yyyy_xyyzzz[k] - g_y_0_yyyy_xyyzzz[k] * ab_y + g_y_0_yyyy_xyyyzzz[k];

                g_y_0_yyyyy_xyzzzz[k] = -g_yyyy_xyzzzz[k] - g_y_0_yyyy_xyzzzz[k] * ab_y + g_y_0_yyyy_xyyzzzz[k];

                g_y_0_yyyyy_xzzzzz[k] = -g_yyyy_xzzzzz[k] - g_y_0_yyyy_xzzzzz[k] * ab_y + g_y_0_yyyy_xyzzzzz[k];

                g_y_0_yyyyy_yyyyyy[k] = -g_yyyy_yyyyyy[k] - g_y_0_yyyy_yyyyyy[k] * ab_y + g_y_0_yyyy_yyyyyyy[k];

                g_y_0_yyyyy_yyyyyz[k] = -g_yyyy_yyyyyz[k] - g_y_0_yyyy_yyyyyz[k] * ab_y + g_y_0_yyyy_yyyyyyz[k];

                g_y_0_yyyyy_yyyyzz[k] = -g_yyyy_yyyyzz[k] - g_y_0_yyyy_yyyyzz[k] * ab_y + g_y_0_yyyy_yyyyyzz[k];

                g_y_0_yyyyy_yyyzzz[k] = -g_yyyy_yyyzzz[k] - g_y_0_yyyy_yyyzzz[k] * ab_y + g_y_0_yyyy_yyyyzzz[k];

                g_y_0_yyyyy_yyzzzz[k] = -g_yyyy_yyzzzz[k] - g_y_0_yyyy_yyzzzz[k] * ab_y + g_y_0_yyyy_yyyzzzz[k];

                g_y_0_yyyyy_yzzzzz[k] = -g_yyyy_yzzzzz[k] - g_y_0_yyyy_yzzzzz[k] * ab_y + g_y_0_yyyy_yyzzzzz[k];

                g_y_0_yyyyy_zzzzzz[k] = -g_yyyy_zzzzzz[k] - g_y_0_yyyy_zzzzzz[k] * ab_y + g_y_0_yyyy_yzzzzzz[k];
            }

            /// Set up 1036-1064 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1036 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1037 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1038 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1039 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1040 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1041 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1042 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1043 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1044 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1045 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1046 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1047 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1048 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1049 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1050 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1051 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1052 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1053 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1054 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1055 * ccomps * dcomps);

            auto g_y_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1056 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1057 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1058 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1059 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1060 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1061 * ccomps * dcomps);

            auto g_y_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1062 * ccomps * dcomps);

            auto g_y_0_yyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1063 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyy_xxxxxx, g_y_0_yyyy_xxxxxxz, g_y_0_yyyy_xxxxxy, g_y_0_yyyy_xxxxxyz, g_y_0_yyyy_xxxxxz, g_y_0_yyyy_xxxxxzz, g_y_0_yyyy_xxxxyy, g_y_0_yyyy_xxxxyyz, g_y_0_yyyy_xxxxyz, g_y_0_yyyy_xxxxyzz, g_y_0_yyyy_xxxxzz, g_y_0_yyyy_xxxxzzz, g_y_0_yyyy_xxxyyy, g_y_0_yyyy_xxxyyyz, g_y_0_yyyy_xxxyyz, g_y_0_yyyy_xxxyyzz, g_y_0_yyyy_xxxyzz, g_y_0_yyyy_xxxyzzz, g_y_0_yyyy_xxxzzz, g_y_0_yyyy_xxxzzzz, g_y_0_yyyy_xxyyyy, g_y_0_yyyy_xxyyyyz, g_y_0_yyyy_xxyyyz, g_y_0_yyyy_xxyyyzz, g_y_0_yyyy_xxyyzz, g_y_0_yyyy_xxyyzzz, g_y_0_yyyy_xxyzzz, g_y_0_yyyy_xxyzzzz, g_y_0_yyyy_xxzzzz, g_y_0_yyyy_xxzzzzz, g_y_0_yyyy_xyyyyy, g_y_0_yyyy_xyyyyyz, g_y_0_yyyy_xyyyyz, g_y_0_yyyy_xyyyyzz, g_y_0_yyyy_xyyyzz, g_y_0_yyyy_xyyyzzz, g_y_0_yyyy_xyyzzz, g_y_0_yyyy_xyyzzzz, g_y_0_yyyy_xyzzzz, g_y_0_yyyy_xyzzzzz, g_y_0_yyyy_xzzzzz, g_y_0_yyyy_xzzzzzz, g_y_0_yyyy_yyyyyy, g_y_0_yyyy_yyyyyyz, g_y_0_yyyy_yyyyyz, g_y_0_yyyy_yyyyyzz, g_y_0_yyyy_yyyyzz, g_y_0_yyyy_yyyyzzz, g_y_0_yyyy_yyyzzz, g_y_0_yyyy_yyyzzzz, g_y_0_yyyy_yyzzzz, g_y_0_yyyy_yyzzzzz, g_y_0_yyyy_yzzzzz, g_y_0_yyyy_yzzzzzz, g_y_0_yyyy_zzzzzz, g_y_0_yyyy_zzzzzzz, g_y_0_yyyyz_xxxxxx, g_y_0_yyyyz_xxxxxy, g_y_0_yyyyz_xxxxxz, g_y_0_yyyyz_xxxxyy, g_y_0_yyyyz_xxxxyz, g_y_0_yyyyz_xxxxzz, g_y_0_yyyyz_xxxyyy, g_y_0_yyyyz_xxxyyz, g_y_0_yyyyz_xxxyzz, g_y_0_yyyyz_xxxzzz, g_y_0_yyyyz_xxyyyy, g_y_0_yyyyz_xxyyyz, g_y_0_yyyyz_xxyyzz, g_y_0_yyyyz_xxyzzz, g_y_0_yyyyz_xxzzzz, g_y_0_yyyyz_xyyyyy, g_y_0_yyyyz_xyyyyz, g_y_0_yyyyz_xyyyzz, g_y_0_yyyyz_xyyzzz, g_y_0_yyyyz_xyzzzz, g_y_0_yyyyz_xzzzzz, g_y_0_yyyyz_yyyyyy, g_y_0_yyyyz_yyyyyz, g_y_0_yyyyz_yyyyzz, g_y_0_yyyyz_yyyzzz, g_y_0_yyyyz_yyzzzz, g_y_0_yyyyz_yzzzzz, g_y_0_yyyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyz_xxxxxx[k] = -g_y_0_yyyy_xxxxxx[k] * ab_z + g_y_0_yyyy_xxxxxxz[k];

                g_y_0_yyyyz_xxxxxy[k] = -g_y_0_yyyy_xxxxxy[k] * ab_z + g_y_0_yyyy_xxxxxyz[k];

                g_y_0_yyyyz_xxxxxz[k] = -g_y_0_yyyy_xxxxxz[k] * ab_z + g_y_0_yyyy_xxxxxzz[k];

                g_y_0_yyyyz_xxxxyy[k] = -g_y_0_yyyy_xxxxyy[k] * ab_z + g_y_0_yyyy_xxxxyyz[k];

                g_y_0_yyyyz_xxxxyz[k] = -g_y_0_yyyy_xxxxyz[k] * ab_z + g_y_0_yyyy_xxxxyzz[k];

                g_y_0_yyyyz_xxxxzz[k] = -g_y_0_yyyy_xxxxzz[k] * ab_z + g_y_0_yyyy_xxxxzzz[k];

                g_y_0_yyyyz_xxxyyy[k] = -g_y_0_yyyy_xxxyyy[k] * ab_z + g_y_0_yyyy_xxxyyyz[k];

                g_y_0_yyyyz_xxxyyz[k] = -g_y_0_yyyy_xxxyyz[k] * ab_z + g_y_0_yyyy_xxxyyzz[k];

                g_y_0_yyyyz_xxxyzz[k] = -g_y_0_yyyy_xxxyzz[k] * ab_z + g_y_0_yyyy_xxxyzzz[k];

                g_y_0_yyyyz_xxxzzz[k] = -g_y_0_yyyy_xxxzzz[k] * ab_z + g_y_0_yyyy_xxxzzzz[k];

                g_y_0_yyyyz_xxyyyy[k] = -g_y_0_yyyy_xxyyyy[k] * ab_z + g_y_0_yyyy_xxyyyyz[k];

                g_y_0_yyyyz_xxyyyz[k] = -g_y_0_yyyy_xxyyyz[k] * ab_z + g_y_0_yyyy_xxyyyzz[k];

                g_y_0_yyyyz_xxyyzz[k] = -g_y_0_yyyy_xxyyzz[k] * ab_z + g_y_0_yyyy_xxyyzzz[k];

                g_y_0_yyyyz_xxyzzz[k] = -g_y_0_yyyy_xxyzzz[k] * ab_z + g_y_0_yyyy_xxyzzzz[k];

                g_y_0_yyyyz_xxzzzz[k] = -g_y_0_yyyy_xxzzzz[k] * ab_z + g_y_0_yyyy_xxzzzzz[k];

                g_y_0_yyyyz_xyyyyy[k] = -g_y_0_yyyy_xyyyyy[k] * ab_z + g_y_0_yyyy_xyyyyyz[k];

                g_y_0_yyyyz_xyyyyz[k] = -g_y_0_yyyy_xyyyyz[k] * ab_z + g_y_0_yyyy_xyyyyzz[k];

                g_y_0_yyyyz_xyyyzz[k] = -g_y_0_yyyy_xyyyzz[k] * ab_z + g_y_0_yyyy_xyyyzzz[k];

                g_y_0_yyyyz_xyyzzz[k] = -g_y_0_yyyy_xyyzzz[k] * ab_z + g_y_0_yyyy_xyyzzzz[k];

                g_y_0_yyyyz_xyzzzz[k] = -g_y_0_yyyy_xyzzzz[k] * ab_z + g_y_0_yyyy_xyzzzzz[k];

                g_y_0_yyyyz_xzzzzz[k] = -g_y_0_yyyy_xzzzzz[k] * ab_z + g_y_0_yyyy_xzzzzzz[k];

                g_y_0_yyyyz_yyyyyy[k] = -g_y_0_yyyy_yyyyyy[k] * ab_z + g_y_0_yyyy_yyyyyyz[k];

                g_y_0_yyyyz_yyyyyz[k] = -g_y_0_yyyy_yyyyyz[k] * ab_z + g_y_0_yyyy_yyyyyzz[k];

                g_y_0_yyyyz_yyyyzz[k] = -g_y_0_yyyy_yyyyzz[k] * ab_z + g_y_0_yyyy_yyyyzzz[k];

                g_y_0_yyyyz_yyyzzz[k] = -g_y_0_yyyy_yyyzzz[k] * ab_z + g_y_0_yyyy_yyyzzzz[k];

                g_y_0_yyyyz_yyzzzz[k] = -g_y_0_yyyy_yyzzzz[k] * ab_z + g_y_0_yyyy_yyzzzzz[k];

                g_y_0_yyyyz_yzzzzz[k] = -g_y_0_yyyy_yzzzzz[k] * ab_z + g_y_0_yyyy_yzzzzzz[k];

                g_y_0_yyyyz_zzzzzz[k] = -g_y_0_yyyy_zzzzzz[k] * ab_z + g_y_0_yyyy_zzzzzzz[k];
            }

            /// Set up 1064-1092 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1064 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1065 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1066 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1067 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1068 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1069 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1070 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1071 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1072 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1073 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1074 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1075 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1076 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1077 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1078 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1079 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1080 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1081 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1082 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1083 * ccomps * dcomps);

            auto g_y_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1084 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1085 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1086 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1087 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1088 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1089 * ccomps * dcomps);

            auto g_y_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1090 * ccomps * dcomps);

            auto g_y_0_yyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1091 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyz_xxxxxx, g_y_0_yyyz_xxxxxxz, g_y_0_yyyz_xxxxxy, g_y_0_yyyz_xxxxxyz, g_y_0_yyyz_xxxxxz, g_y_0_yyyz_xxxxxzz, g_y_0_yyyz_xxxxyy, g_y_0_yyyz_xxxxyyz, g_y_0_yyyz_xxxxyz, g_y_0_yyyz_xxxxyzz, g_y_0_yyyz_xxxxzz, g_y_0_yyyz_xxxxzzz, g_y_0_yyyz_xxxyyy, g_y_0_yyyz_xxxyyyz, g_y_0_yyyz_xxxyyz, g_y_0_yyyz_xxxyyzz, g_y_0_yyyz_xxxyzz, g_y_0_yyyz_xxxyzzz, g_y_0_yyyz_xxxzzz, g_y_0_yyyz_xxxzzzz, g_y_0_yyyz_xxyyyy, g_y_0_yyyz_xxyyyyz, g_y_0_yyyz_xxyyyz, g_y_0_yyyz_xxyyyzz, g_y_0_yyyz_xxyyzz, g_y_0_yyyz_xxyyzzz, g_y_0_yyyz_xxyzzz, g_y_0_yyyz_xxyzzzz, g_y_0_yyyz_xxzzzz, g_y_0_yyyz_xxzzzzz, g_y_0_yyyz_xyyyyy, g_y_0_yyyz_xyyyyyz, g_y_0_yyyz_xyyyyz, g_y_0_yyyz_xyyyyzz, g_y_0_yyyz_xyyyzz, g_y_0_yyyz_xyyyzzz, g_y_0_yyyz_xyyzzz, g_y_0_yyyz_xyyzzzz, g_y_0_yyyz_xyzzzz, g_y_0_yyyz_xyzzzzz, g_y_0_yyyz_xzzzzz, g_y_0_yyyz_xzzzzzz, g_y_0_yyyz_yyyyyy, g_y_0_yyyz_yyyyyyz, g_y_0_yyyz_yyyyyz, g_y_0_yyyz_yyyyyzz, g_y_0_yyyz_yyyyzz, g_y_0_yyyz_yyyyzzz, g_y_0_yyyz_yyyzzz, g_y_0_yyyz_yyyzzzz, g_y_0_yyyz_yyzzzz, g_y_0_yyyz_yyzzzzz, g_y_0_yyyz_yzzzzz, g_y_0_yyyz_yzzzzzz, g_y_0_yyyz_zzzzzz, g_y_0_yyyz_zzzzzzz, g_y_0_yyyzz_xxxxxx, g_y_0_yyyzz_xxxxxy, g_y_0_yyyzz_xxxxxz, g_y_0_yyyzz_xxxxyy, g_y_0_yyyzz_xxxxyz, g_y_0_yyyzz_xxxxzz, g_y_0_yyyzz_xxxyyy, g_y_0_yyyzz_xxxyyz, g_y_0_yyyzz_xxxyzz, g_y_0_yyyzz_xxxzzz, g_y_0_yyyzz_xxyyyy, g_y_0_yyyzz_xxyyyz, g_y_0_yyyzz_xxyyzz, g_y_0_yyyzz_xxyzzz, g_y_0_yyyzz_xxzzzz, g_y_0_yyyzz_xyyyyy, g_y_0_yyyzz_xyyyyz, g_y_0_yyyzz_xyyyzz, g_y_0_yyyzz_xyyzzz, g_y_0_yyyzz_xyzzzz, g_y_0_yyyzz_xzzzzz, g_y_0_yyyzz_yyyyyy, g_y_0_yyyzz_yyyyyz, g_y_0_yyyzz_yyyyzz, g_y_0_yyyzz_yyyzzz, g_y_0_yyyzz_yyzzzz, g_y_0_yyyzz_yzzzzz, g_y_0_yyyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzz_xxxxxx[k] = -g_y_0_yyyz_xxxxxx[k] * ab_z + g_y_0_yyyz_xxxxxxz[k];

                g_y_0_yyyzz_xxxxxy[k] = -g_y_0_yyyz_xxxxxy[k] * ab_z + g_y_0_yyyz_xxxxxyz[k];

                g_y_0_yyyzz_xxxxxz[k] = -g_y_0_yyyz_xxxxxz[k] * ab_z + g_y_0_yyyz_xxxxxzz[k];

                g_y_0_yyyzz_xxxxyy[k] = -g_y_0_yyyz_xxxxyy[k] * ab_z + g_y_0_yyyz_xxxxyyz[k];

                g_y_0_yyyzz_xxxxyz[k] = -g_y_0_yyyz_xxxxyz[k] * ab_z + g_y_0_yyyz_xxxxyzz[k];

                g_y_0_yyyzz_xxxxzz[k] = -g_y_0_yyyz_xxxxzz[k] * ab_z + g_y_0_yyyz_xxxxzzz[k];

                g_y_0_yyyzz_xxxyyy[k] = -g_y_0_yyyz_xxxyyy[k] * ab_z + g_y_0_yyyz_xxxyyyz[k];

                g_y_0_yyyzz_xxxyyz[k] = -g_y_0_yyyz_xxxyyz[k] * ab_z + g_y_0_yyyz_xxxyyzz[k];

                g_y_0_yyyzz_xxxyzz[k] = -g_y_0_yyyz_xxxyzz[k] * ab_z + g_y_0_yyyz_xxxyzzz[k];

                g_y_0_yyyzz_xxxzzz[k] = -g_y_0_yyyz_xxxzzz[k] * ab_z + g_y_0_yyyz_xxxzzzz[k];

                g_y_0_yyyzz_xxyyyy[k] = -g_y_0_yyyz_xxyyyy[k] * ab_z + g_y_0_yyyz_xxyyyyz[k];

                g_y_0_yyyzz_xxyyyz[k] = -g_y_0_yyyz_xxyyyz[k] * ab_z + g_y_0_yyyz_xxyyyzz[k];

                g_y_0_yyyzz_xxyyzz[k] = -g_y_0_yyyz_xxyyzz[k] * ab_z + g_y_0_yyyz_xxyyzzz[k];

                g_y_0_yyyzz_xxyzzz[k] = -g_y_0_yyyz_xxyzzz[k] * ab_z + g_y_0_yyyz_xxyzzzz[k];

                g_y_0_yyyzz_xxzzzz[k] = -g_y_0_yyyz_xxzzzz[k] * ab_z + g_y_0_yyyz_xxzzzzz[k];

                g_y_0_yyyzz_xyyyyy[k] = -g_y_0_yyyz_xyyyyy[k] * ab_z + g_y_0_yyyz_xyyyyyz[k];

                g_y_0_yyyzz_xyyyyz[k] = -g_y_0_yyyz_xyyyyz[k] * ab_z + g_y_0_yyyz_xyyyyzz[k];

                g_y_0_yyyzz_xyyyzz[k] = -g_y_0_yyyz_xyyyzz[k] * ab_z + g_y_0_yyyz_xyyyzzz[k];

                g_y_0_yyyzz_xyyzzz[k] = -g_y_0_yyyz_xyyzzz[k] * ab_z + g_y_0_yyyz_xyyzzzz[k];

                g_y_0_yyyzz_xyzzzz[k] = -g_y_0_yyyz_xyzzzz[k] * ab_z + g_y_0_yyyz_xyzzzzz[k];

                g_y_0_yyyzz_xzzzzz[k] = -g_y_0_yyyz_xzzzzz[k] * ab_z + g_y_0_yyyz_xzzzzzz[k];

                g_y_0_yyyzz_yyyyyy[k] = -g_y_0_yyyz_yyyyyy[k] * ab_z + g_y_0_yyyz_yyyyyyz[k];

                g_y_0_yyyzz_yyyyyz[k] = -g_y_0_yyyz_yyyyyz[k] * ab_z + g_y_0_yyyz_yyyyyzz[k];

                g_y_0_yyyzz_yyyyzz[k] = -g_y_0_yyyz_yyyyzz[k] * ab_z + g_y_0_yyyz_yyyyzzz[k];

                g_y_0_yyyzz_yyyzzz[k] = -g_y_0_yyyz_yyyzzz[k] * ab_z + g_y_0_yyyz_yyyzzzz[k];

                g_y_0_yyyzz_yyzzzz[k] = -g_y_0_yyyz_yyzzzz[k] * ab_z + g_y_0_yyyz_yyzzzzz[k];

                g_y_0_yyyzz_yzzzzz[k] = -g_y_0_yyyz_yzzzzz[k] * ab_z + g_y_0_yyyz_yzzzzzz[k];

                g_y_0_yyyzz_zzzzzz[k] = -g_y_0_yyyz_zzzzzz[k] * ab_z + g_y_0_yyyz_zzzzzzz[k];
            }

            /// Set up 1092-1120 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1092 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1093 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1094 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1095 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1096 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1097 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1098 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1099 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1100 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1101 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1102 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1103 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1104 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1105 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1106 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1107 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1108 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1109 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1110 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1111 * ccomps * dcomps);

            auto g_y_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1112 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1113 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1114 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1115 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1116 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1117 * ccomps * dcomps);

            auto g_y_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1118 * ccomps * dcomps);

            auto g_y_0_yyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1119 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyzz_xxxxxx, g_y_0_yyzz_xxxxxxz, g_y_0_yyzz_xxxxxy, g_y_0_yyzz_xxxxxyz, g_y_0_yyzz_xxxxxz, g_y_0_yyzz_xxxxxzz, g_y_0_yyzz_xxxxyy, g_y_0_yyzz_xxxxyyz, g_y_0_yyzz_xxxxyz, g_y_0_yyzz_xxxxyzz, g_y_0_yyzz_xxxxzz, g_y_0_yyzz_xxxxzzz, g_y_0_yyzz_xxxyyy, g_y_0_yyzz_xxxyyyz, g_y_0_yyzz_xxxyyz, g_y_0_yyzz_xxxyyzz, g_y_0_yyzz_xxxyzz, g_y_0_yyzz_xxxyzzz, g_y_0_yyzz_xxxzzz, g_y_0_yyzz_xxxzzzz, g_y_0_yyzz_xxyyyy, g_y_0_yyzz_xxyyyyz, g_y_0_yyzz_xxyyyz, g_y_0_yyzz_xxyyyzz, g_y_0_yyzz_xxyyzz, g_y_0_yyzz_xxyyzzz, g_y_0_yyzz_xxyzzz, g_y_0_yyzz_xxyzzzz, g_y_0_yyzz_xxzzzz, g_y_0_yyzz_xxzzzzz, g_y_0_yyzz_xyyyyy, g_y_0_yyzz_xyyyyyz, g_y_0_yyzz_xyyyyz, g_y_0_yyzz_xyyyyzz, g_y_0_yyzz_xyyyzz, g_y_0_yyzz_xyyyzzz, g_y_0_yyzz_xyyzzz, g_y_0_yyzz_xyyzzzz, g_y_0_yyzz_xyzzzz, g_y_0_yyzz_xyzzzzz, g_y_0_yyzz_xzzzzz, g_y_0_yyzz_xzzzzzz, g_y_0_yyzz_yyyyyy, g_y_0_yyzz_yyyyyyz, g_y_0_yyzz_yyyyyz, g_y_0_yyzz_yyyyyzz, g_y_0_yyzz_yyyyzz, g_y_0_yyzz_yyyyzzz, g_y_0_yyzz_yyyzzz, g_y_0_yyzz_yyyzzzz, g_y_0_yyzz_yyzzzz, g_y_0_yyzz_yyzzzzz, g_y_0_yyzz_yzzzzz, g_y_0_yyzz_yzzzzzz, g_y_0_yyzz_zzzzzz, g_y_0_yyzz_zzzzzzz, g_y_0_yyzzz_xxxxxx, g_y_0_yyzzz_xxxxxy, g_y_0_yyzzz_xxxxxz, g_y_0_yyzzz_xxxxyy, g_y_0_yyzzz_xxxxyz, g_y_0_yyzzz_xxxxzz, g_y_0_yyzzz_xxxyyy, g_y_0_yyzzz_xxxyyz, g_y_0_yyzzz_xxxyzz, g_y_0_yyzzz_xxxzzz, g_y_0_yyzzz_xxyyyy, g_y_0_yyzzz_xxyyyz, g_y_0_yyzzz_xxyyzz, g_y_0_yyzzz_xxyzzz, g_y_0_yyzzz_xxzzzz, g_y_0_yyzzz_xyyyyy, g_y_0_yyzzz_xyyyyz, g_y_0_yyzzz_xyyyzz, g_y_0_yyzzz_xyyzzz, g_y_0_yyzzz_xyzzzz, g_y_0_yyzzz_xzzzzz, g_y_0_yyzzz_yyyyyy, g_y_0_yyzzz_yyyyyz, g_y_0_yyzzz_yyyyzz, g_y_0_yyzzz_yyyzzz, g_y_0_yyzzz_yyzzzz, g_y_0_yyzzz_yzzzzz, g_y_0_yyzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzz_xxxxxx[k] = -g_y_0_yyzz_xxxxxx[k] * ab_z + g_y_0_yyzz_xxxxxxz[k];

                g_y_0_yyzzz_xxxxxy[k] = -g_y_0_yyzz_xxxxxy[k] * ab_z + g_y_0_yyzz_xxxxxyz[k];

                g_y_0_yyzzz_xxxxxz[k] = -g_y_0_yyzz_xxxxxz[k] * ab_z + g_y_0_yyzz_xxxxxzz[k];

                g_y_0_yyzzz_xxxxyy[k] = -g_y_0_yyzz_xxxxyy[k] * ab_z + g_y_0_yyzz_xxxxyyz[k];

                g_y_0_yyzzz_xxxxyz[k] = -g_y_0_yyzz_xxxxyz[k] * ab_z + g_y_0_yyzz_xxxxyzz[k];

                g_y_0_yyzzz_xxxxzz[k] = -g_y_0_yyzz_xxxxzz[k] * ab_z + g_y_0_yyzz_xxxxzzz[k];

                g_y_0_yyzzz_xxxyyy[k] = -g_y_0_yyzz_xxxyyy[k] * ab_z + g_y_0_yyzz_xxxyyyz[k];

                g_y_0_yyzzz_xxxyyz[k] = -g_y_0_yyzz_xxxyyz[k] * ab_z + g_y_0_yyzz_xxxyyzz[k];

                g_y_0_yyzzz_xxxyzz[k] = -g_y_0_yyzz_xxxyzz[k] * ab_z + g_y_0_yyzz_xxxyzzz[k];

                g_y_0_yyzzz_xxxzzz[k] = -g_y_0_yyzz_xxxzzz[k] * ab_z + g_y_0_yyzz_xxxzzzz[k];

                g_y_0_yyzzz_xxyyyy[k] = -g_y_0_yyzz_xxyyyy[k] * ab_z + g_y_0_yyzz_xxyyyyz[k];

                g_y_0_yyzzz_xxyyyz[k] = -g_y_0_yyzz_xxyyyz[k] * ab_z + g_y_0_yyzz_xxyyyzz[k];

                g_y_0_yyzzz_xxyyzz[k] = -g_y_0_yyzz_xxyyzz[k] * ab_z + g_y_0_yyzz_xxyyzzz[k];

                g_y_0_yyzzz_xxyzzz[k] = -g_y_0_yyzz_xxyzzz[k] * ab_z + g_y_0_yyzz_xxyzzzz[k];

                g_y_0_yyzzz_xxzzzz[k] = -g_y_0_yyzz_xxzzzz[k] * ab_z + g_y_0_yyzz_xxzzzzz[k];

                g_y_0_yyzzz_xyyyyy[k] = -g_y_0_yyzz_xyyyyy[k] * ab_z + g_y_0_yyzz_xyyyyyz[k];

                g_y_0_yyzzz_xyyyyz[k] = -g_y_0_yyzz_xyyyyz[k] * ab_z + g_y_0_yyzz_xyyyyzz[k];

                g_y_0_yyzzz_xyyyzz[k] = -g_y_0_yyzz_xyyyzz[k] * ab_z + g_y_0_yyzz_xyyyzzz[k];

                g_y_0_yyzzz_xyyzzz[k] = -g_y_0_yyzz_xyyzzz[k] * ab_z + g_y_0_yyzz_xyyzzzz[k];

                g_y_0_yyzzz_xyzzzz[k] = -g_y_0_yyzz_xyzzzz[k] * ab_z + g_y_0_yyzz_xyzzzzz[k];

                g_y_0_yyzzz_xzzzzz[k] = -g_y_0_yyzz_xzzzzz[k] * ab_z + g_y_0_yyzz_xzzzzzz[k];

                g_y_0_yyzzz_yyyyyy[k] = -g_y_0_yyzz_yyyyyy[k] * ab_z + g_y_0_yyzz_yyyyyyz[k];

                g_y_0_yyzzz_yyyyyz[k] = -g_y_0_yyzz_yyyyyz[k] * ab_z + g_y_0_yyzz_yyyyyzz[k];

                g_y_0_yyzzz_yyyyzz[k] = -g_y_0_yyzz_yyyyzz[k] * ab_z + g_y_0_yyzz_yyyyzzz[k];

                g_y_0_yyzzz_yyyzzz[k] = -g_y_0_yyzz_yyyzzz[k] * ab_z + g_y_0_yyzz_yyyzzzz[k];

                g_y_0_yyzzz_yyzzzz[k] = -g_y_0_yyzz_yyzzzz[k] * ab_z + g_y_0_yyzz_yyzzzzz[k];

                g_y_0_yyzzz_yzzzzz[k] = -g_y_0_yyzz_yzzzzz[k] * ab_z + g_y_0_yyzz_yzzzzzz[k];

                g_y_0_yyzzz_zzzzzz[k] = -g_y_0_yyzz_zzzzzz[k] * ab_z + g_y_0_yyzz_zzzzzzz[k];
            }

            /// Set up 1120-1148 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1120 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1121 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1122 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1123 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1124 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1125 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1126 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1127 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1128 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1129 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1130 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1131 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1132 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1133 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1134 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1135 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1136 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1137 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1138 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1139 * ccomps * dcomps);

            auto g_y_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1140 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1141 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1142 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1143 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1144 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1145 * ccomps * dcomps);

            auto g_y_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1146 * ccomps * dcomps);

            auto g_y_0_yzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1147 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzzz_xxxxxx, g_y_0_yzzz_xxxxxxz, g_y_0_yzzz_xxxxxy, g_y_0_yzzz_xxxxxyz, g_y_0_yzzz_xxxxxz, g_y_0_yzzz_xxxxxzz, g_y_0_yzzz_xxxxyy, g_y_0_yzzz_xxxxyyz, g_y_0_yzzz_xxxxyz, g_y_0_yzzz_xxxxyzz, g_y_0_yzzz_xxxxzz, g_y_0_yzzz_xxxxzzz, g_y_0_yzzz_xxxyyy, g_y_0_yzzz_xxxyyyz, g_y_0_yzzz_xxxyyz, g_y_0_yzzz_xxxyyzz, g_y_0_yzzz_xxxyzz, g_y_0_yzzz_xxxyzzz, g_y_0_yzzz_xxxzzz, g_y_0_yzzz_xxxzzzz, g_y_0_yzzz_xxyyyy, g_y_0_yzzz_xxyyyyz, g_y_0_yzzz_xxyyyz, g_y_0_yzzz_xxyyyzz, g_y_0_yzzz_xxyyzz, g_y_0_yzzz_xxyyzzz, g_y_0_yzzz_xxyzzz, g_y_0_yzzz_xxyzzzz, g_y_0_yzzz_xxzzzz, g_y_0_yzzz_xxzzzzz, g_y_0_yzzz_xyyyyy, g_y_0_yzzz_xyyyyyz, g_y_0_yzzz_xyyyyz, g_y_0_yzzz_xyyyyzz, g_y_0_yzzz_xyyyzz, g_y_0_yzzz_xyyyzzz, g_y_0_yzzz_xyyzzz, g_y_0_yzzz_xyyzzzz, g_y_0_yzzz_xyzzzz, g_y_0_yzzz_xyzzzzz, g_y_0_yzzz_xzzzzz, g_y_0_yzzz_xzzzzzz, g_y_0_yzzz_yyyyyy, g_y_0_yzzz_yyyyyyz, g_y_0_yzzz_yyyyyz, g_y_0_yzzz_yyyyyzz, g_y_0_yzzz_yyyyzz, g_y_0_yzzz_yyyyzzz, g_y_0_yzzz_yyyzzz, g_y_0_yzzz_yyyzzzz, g_y_0_yzzz_yyzzzz, g_y_0_yzzz_yyzzzzz, g_y_0_yzzz_yzzzzz, g_y_0_yzzz_yzzzzzz, g_y_0_yzzz_zzzzzz, g_y_0_yzzz_zzzzzzz, g_y_0_yzzzz_xxxxxx, g_y_0_yzzzz_xxxxxy, g_y_0_yzzzz_xxxxxz, g_y_0_yzzzz_xxxxyy, g_y_0_yzzzz_xxxxyz, g_y_0_yzzzz_xxxxzz, g_y_0_yzzzz_xxxyyy, g_y_0_yzzzz_xxxyyz, g_y_0_yzzzz_xxxyzz, g_y_0_yzzzz_xxxzzz, g_y_0_yzzzz_xxyyyy, g_y_0_yzzzz_xxyyyz, g_y_0_yzzzz_xxyyzz, g_y_0_yzzzz_xxyzzz, g_y_0_yzzzz_xxzzzz, g_y_0_yzzzz_xyyyyy, g_y_0_yzzzz_xyyyyz, g_y_0_yzzzz_xyyyzz, g_y_0_yzzzz_xyyzzz, g_y_0_yzzzz_xyzzzz, g_y_0_yzzzz_xzzzzz, g_y_0_yzzzz_yyyyyy, g_y_0_yzzzz_yyyyyz, g_y_0_yzzzz_yyyyzz, g_y_0_yzzzz_yyyzzz, g_y_0_yzzzz_yyzzzz, g_y_0_yzzzz_yzzzzz, g_y_0_yzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzz_xxxxxx[k] = -g_y_0_yzzz_xxxxxx[k] * ab_z + g_y_0_yzzz_xxxxxxz[k];

                g_y_0_yzzzz_xxxxxy[k] = -g_y_0_yzzz_xxxxxy[k] * ab_z + g_y_0_yzzz_xxxxxyz[k];

                g_y_0_yzzzz_xxxxxz[k] = -g_y_0_yzzz_xxxxxz[k] * ab_z + g_y_0_yzzz_xxxxxzz[k];

                g_y_0_yzzzz_xxxxyy[k] = -g_y_0_yzzz_xxxxyy[k] * ab_z + g_y_0_yzzz_xxxxyyz[k];

                g_y_0_yzzzz_xxxxyz[k] = -g_y_0_yzzz_xxxxyz[k] * ab_z + g_y_0_yzzz_xxxxyzz[k];

                g_y_0_yzzzz_xxxxzz[k] = -g_y_0_yzzz_xxxxzz[k] * ab_z + g_y_0_yzzz_xxxxzzz[k];

                g_y_0_yzzzz_xxxyyy[k] = -g_y_0_yzzz_xxxyyy[k] * ab_z + g_y_0_yzzz_xxxyyyz[k];

                g_y_0_yzzzz_xxxyyz[k] = -g_y_0_yzzz_xxxyyz[k] * ab_z + g_y_0_yzzz_xxxyyzz[k];

                g_y_0_yzzzz_xxxyzz[k] = -g_y_0_yzzz_xxxyzz[k] * ab_z + g_y_0_yzzz_xxxyzzz[k];

                g_y_0_yzzzz_xxxzzz[k] = -g_y_0_yzzz_xxxzzz[k] * ab_z + g_y_0_yzzz_xxxzzzz[k];

                g_y_0_yzzzz_xxyyyy[k] = -g_y_0_yzzz_xxyyyy[k] * ab_z + g_y_0_yzzz_xxyyyyz[k];

                g_y_0_yzzzz_xxyyyz[k] = -g_y_0_yzzz_xxyyyz[k] * ab_z + g_y_0_yzzz_xxyyyzz[k];

                g_y_0_yzzzz_xxyyzz[k] = -g_y_0_yzzz_xxyyzz[k] * ab_z + g_y_0_yzzz_xxyyzzz[k];

                g_y_0_yzzzz_xxyzzz[k] = -g_y_0_yzzz_xxyzzz[k] * ab_z + g_y_0_yzzz_xxyzzzz[k];

                g_y_0_yzzzz_xxzzzz[k] = -g_y_0_yzzz_xxzzzz[k] * ab_z + g_y_0_yzzz_xxzzzzz[k];

                g_y_0_yzzzz_xyyyyy[k] = -g_y_0_yzzz_xyyyyy[k] * ab_z + g_y_0_yzzz_xyyyyyz[k];

                g_y_0_yzzzz_xyyyyz[k] = -g_y_0_yzzz_xyyyyz[k] * ab_z + g_y_0_yzzz_xyyyyzz[k];

                g_y_0_yzzzz_xyyyzz[k] = -g_y_0_yzzz_xyyyzz[k] * ab_z + g_y_0_yzzz_xyyyzzz[k];

                g_y_0_yzzzz_xyyzzz[k] = -g_y_0_yzzz_xyyzzz[k] * ab_z + g_y_0_yzzz_xyyzzzz[k];

                g_y_0_yzzzz_xyzzzz[k] = -g_y_0_yzzz_xyzzzz[k] * ab_z + g_y_0_yzzz_xyzzzzz[k];

                g_y_0_yzzzz_xzzzzz[k] = -g_y_0_yzzz_xzzzzz[k] * ab_z + g_y_0_yzzz_xzzzzzz[k];

                g_y_0_yzzzz_yyyyyy[k] = -g_y_0_yzzz_yyyyyy[k] * ab_z + g_y_0_yzzz_yyyyyyz[k];

                g_y_0_yzzzz_yyyyyz[k] = -g_y_0_yzzz_yyyyyz[k] * ab_z + g_y_0_yzzz_yyyyyzz[k];

                g_y_0_yzzzz_yyyyzz[k] = -g_y_0_yzzz_yyyyzz[k] * ab_z + g_y_0_yzzz_yyyyzzz[k];

                g_y_0_yzzzz_yyyzzz[k] = -g_y_0_yzzz_yyyzzz[k] * ab_z + g_y_0_yzzz_yyyzzzz[k];

                g_y_0_yzzzz_yyzzzz[k] = -g_y_0_yzzz_yyzzzz[k] * ab_z + g_y_0_yzzz_yyzzzzz[k];

                g_y_0_yzzzz_yzzzzz[k] = -g_y_0_yzzz_yzzzzz[k] * ab_z + g_y_0_yzzz_yzzzzzz[k];

                g_y_0_yzzzz_zzzzzz[k] = -g_y_0_yzzz_zzzzzz[k] * ab_z + g_y_0_yzzz_zzzzzzz[k];
            }

            /// Set up 1148-1176 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1148 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1149 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1150 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1151 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1152 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1153 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1154 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1155 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1156 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1157 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1158 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1159 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1160 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1161 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1162 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1163 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1164 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1165 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1166 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1167 * ccomps * dcomps);

            auto g_y_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1168 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1169 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1170 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1171 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1172 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1173 * ccomps * dcomps);

            auto g_y_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1174 * ccomps * dcomps);

            auto g_y_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1175 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzz_xxxxxx, g_y_0_zzzz_xxxxxxz, g_y_0_zzzz_xxxxxy, g_y_0_zzzz_xxxxxyz, g_y_0_zzzz_xxxxxz, g_y_0_zzzz_xxxxxzz, g_y_0_zzzz_xxxxyy, g_y_0_zzzz_xxxxyyz, g_y_0_zzzz_xxxxyz, g_y_0_zzzz_xxxxyzz, g_y_0_zzzz_xxxxzz, g_y_0_zzzz_xxxxzzz, g_y_0_zzzz_xxxyyy, g_y_0_zzzz_xxxyyyz, g_y_0_zzzz_xxxyyz, g_y_0_zzzz_xxxyyzz, g_y_0_zzzz_xxxyzz, g_y_0_zzzz_xxxyzzz, g_y_0_zzzz_xxxzzz, g_y_0_zzzz_xxxzzzz, g_y_0_zzzz_xxyyyy, g_y_0_zzzz_xxyyyyz, g_y_0_zzzz_xxyyyz, g_y_0_zzzz_xxyyyzz, g_y_0_zzzz_xxyyzz, g_y_0_zzzz_xxyyzzz, g_y_0_zzzz_xxyzzz, g_y_0_zzzz_xxyzzzz, g_y_0_zzzz_xxzzzz, g_y_0_zzzz_xxzzzzz, g_y_0_zzzz_xyyyyy, g_y_0_zzzz_xyyyyyz, g_y_0_zzzz_xyyyyz, g_y_0_zzzz_xyyyyzz, g_y_0_zzzz_xyyyzz, g_y_0_zzzz_xyyyzzz, g_y_0_zzzz_xyyzzz, g_y_0_zzzz_xyyzzzz, g_y_0_zzzz_xyzzzz, g_y_0_zzzz_xyzzzzz, g_y_0_zzzz_xzzzzz, g_y_0_zzzz_xzzzzzz, g_y_0_zzzz_yyyyyy, g_y_0_zzzz_yyyyyyz, g_y_0_zzzz_yyyyyz, g_y_0_zzzz_yyyyyzz, g_y_0_zzzz_yyyyzz, g_y_0_zzzz_yyyyzzz, g_y_0_zzzz_yyyzzz, g_y_0_zzzz_yyyzzzz, g_y_0_zzzz_yyzzzz, g_y_0_zzzz_yyzzzzz, g_y_0_zzzz_yzzzzz, g_y_0_zzzz_yzzzzzz, g_y_0_zzzz_zzzzzz, g_y_0_zzzz_zzzzzzz, g_y_0_zzzzz_xxxxxx, g_y_0_zzzzz_xxxxxy, g_y_0_zzzzz_xxxxxz, g_y_0_zzzzz_xxxxyy, g_y_0_zzzzz_xxxxyz, g_y_0_zzzzz_xxxxzz, g_y_0_zzzzz_xxxyyy, g_y_0_zzzzz_xxxyyz, g_y_0_zzzzz_xxxyzz, g_y_0_zzzzz_xxxzzz, g_y_0_zzzzz_xxyyyy, g_y_0_zzzzz_xxyyyz, g_y_0_zzzzz_xxyyzz, g_y_0_zzzzz_xxyzzz, g_y_0_zzzzz_xxzzzz, g_y_0_zzzzz_xyyyyy, g_y_0_zzzzz_xyyyyz, g_y_0_zzzzz_xyyyzz, g_y_0_zzzzz_xyyzzz, g_y_0_zzzzz_xyzzzz, g_y_0_zzzzz_xzzzzz, g_y_0_zzzzz_yyyyyy, g_y_0_zzzzz_yyyyyz, g_y_0_zzzzz_yyyyzz, g_y_0_zzzzz_yyyzzz, g_y_0_zzzzz_yyzzzz, g_y_0_zzzzz_yzzzzz, g_y_0_zzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzz_xxxxxx[k] = -g_y_0_zzzz_xxxxxx[k] * ab_z + g_y_0_zzzz_xxxxxxz[k];

                g_y_0_zzzzz_xxxxxy[k] = -g_y_0_zzzz_xxxxxy[k] * ab_z + g_y_0_zzzz_xxxxxyz[k];

                g_y_0_zzzzz_xxxxxz[k] = -g_y_0_zzzz_xxxxxz[k] * ab_z + g_y_0_zzzz_xxxxxzz[k];

                g_y_0_zzzzz_xxxxyy[k] = -g_y_0_zzzz_xxxxyy[k] * ab_z + g_y_0_zzzz_xxxxyyz[k];

                g_y_0_zzzzz_xxxxyz[k] = -g_y_0_zzzz_xxxxyz[k] * ab_z + g_y_0_zzzz_xxxxyzz[k];

                g_y_0_zzzzz_xxxxzz[k] = -g_y_0_zzzz_xxxxzz[k] * ab_z + g_y_0_zzzz_xxxxzzz[k];

                g_y_0_zzzzz_xxxyyy[k] = -g_y_0_zzzz_xxxyyy[k] * ab_z + g_y_0_zzzz_xxxyyyz[k];

                g_y_0_zzzzz_xxxyyz[k] = -g_y_0_zzzz_xxxyyz[k] * ab_z + g_y_0_zzzz_xxxyyzz[k];

                g_y_0_zzzzz_xxxyzz[k] = -g_y_0_zzzz_xxxyzz[k] * ab_z + g_y_0_zzzz_xxxyzzz[k];

                g_y_0_zzzzz_xxxzzz[k] = -g_y_0_zzzz_xxxzzz[k] * ab_z + g_y_0_zzzz_xxxzzzz[k];

                g_y_0_zzzzz_xxyyyy[k] = -g_y_0_zzzz_xxyyyy[k] * ab_z + g_y_0_zzzz_xxyyyyz[k];

                g_y_0_zzzzz_xxyyyz[k] = -g_y_0_zzzz_xxyyyz[k] * ab_z + g_y_0_zzzz_xxyyyzz[k];

                g_y_0_zzzzz_xxyyzz[k] = -g_y_0_zzzz_xxyyzz[k] * ab_z + g_y_0_zzzz_xxyyzzz[k];

                g_y_0_zzzzz_xxyzzz[k] = -g_y_0_zzzz_xxyzzz[k] * ab_z + g_y_0_zzzz_xxyzzzz[k];

                g_y_0_zzzzz_xxzzzz[k] = -g_y_0_zzzz_xxzzzz[k] * ab_z + g_y_0_zzzz_xxzzzzz[k];

                g_y_0_zzzzz_xyyyyy[k] = -g_y_0_zzzz_xyyyyy[k] * ab_z + g_y_0_zzzz_xyyyyyz[k];

                g_y_0_zzzzz_xyyyyz[k] = -g_y_0_zzzz_xyyyyz[k] * ab_z + g_y_0_zzzz_xyyyyzz[k];

                g_y_0_zzzzz_xyyyzz[k] = -g_y_0_zzzz_xyyyzz[k] * ab_z + g_y_0_zzzz_xyyyzzz[k];

                g_y_0_zzzzz_xyyzzz[k] = -g_y_0_zzzz_xyyzzz[k] * ab_z + g_y_0_zzzz_xyyzzzz[k];

                g_y_0_zzzzz_xyzzzz[k] = -g_y_0_zzzz_xyzzzz[k] * ab_z + g_y_0_zzzz_xyzzzzz[k];

                g_y_0_zzzzz_xzzzzz[k] = -g_y_0_zzzz_xzzzzz[k] * ab_z + g_y_0_zzzz_xzzzzzz[k];

                g_y_0_zzzzz_yyyyyy[k] = -g_y_0_zzzz_yyyyyy[k] * ab_z + g_y_0_zzzz_yyyyyyz[k];

                g_y_0_zzzzz_yyyyyz[k] = -g_y_0_zzzz_yyyyyz[k] * ab_z + g_y_0_zzzz_yyyyyzz[k];

                g_y_0_zzzzz_yyyyzz[k] = -g_y_0_zzzz_yyyyzz[k] * ab_z + g_y_0_zzzz_yyyyzzz[k];

                g_y_0_zzzzz_yyyzzz[k] = -g_y_0_zzzz_yyyzzz[k] * ab_z + g_y_0_zzzz_yyyzzzz[k];

                g_y_0_zzzzz_yyzzzz[k] = -g_y_0_zzzz_yyzzzz[k] * ab_z + g_y_0_zzzz_yyzzzzz[k];

                g_y_0_zzzzz_yzzzzz[k] = -g_y_0_zzzz_yzzzzz[k] * ab_z + g_y_0_zzzz_yzzzzzz[k];

                g_y_0_zzzzz_zzzzzz[k] = -g_y_0_zzzz_zzzzzz[k] * ab_z + g_y_0_zzzz_zzzzzzz[k];
            }

            /// Set up 1176-1204 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxx_xxxxxx = cbuffer.data(hi_geom_10_off + 1176 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxy = cbuffer.data(hi_geom_10_off + 1177 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxz = cbuffer.data(hi_geom_10_off + 1178 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxyy = cbuffer.data(hi_geom_10_off + 1179 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxyz = cbuffer.data(hi_geom_10_off + 1180 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxzz = cbuffer.data(hi_geom_10_off + 1181 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyyy = cbuffer.data(hi_geom_10_off + 1182 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyyz = cbuffer.data(hi_geom_10_off + 1183 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyzz = cbuffer.data(hi_geom_10_off + 1184 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxzzz = cbuffer.data(hi_geom_10_off + 1185 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyyy = cbuffer.data(hi_geom_10_off + 1186 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyyz = cbuffer.data(hi_geom_10_off + 1187 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyzz = cbuffer.data(hi_geom_10_off + 1188 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyzzz = cbuffer.data(hi_geom_10_off + 1189 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxzzzz = cbuffer.data(hi_geom_10_off + 1190 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyyy = cbuffer.data(hi_geom_10_off + 1191 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyyz = cbuffer.data(hi_geom_10_off + 1192 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyzz = cbuffer.data(hi_geom_10_off + 1193 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyzzz = cbuffer.data(hi_geom_10_off + 1194 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyzzzz = cbuffer.data(hi_geom_10_off + 1195 * ccomps * dcomps);

            auto g_z_0_xxxxx_xzzzzz = cbuffer.data(hi_geom_10_off + 1196 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyyyy = cbuffer.data(hi_geom_10_off + 1197 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyyyz = cbuffer.data(hi_geom_10_off + 1198 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyyzz = cbuffer.data(hi_geom_10_off + 1199 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyzzz = cbuffer.data(hi_geom_10_off + 1200 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyzzzz = cbuffer.data(hi_geom_10_off + 1201 * ccomps * dcomps);

            auto g_z_0_xxxxx_yzzzzz = cbuffer.data(hi_geom_10_off + 1202 * ccomps * dcomps);

            auto g_z_0_xxxxx_zzzzzz = cbuffer.data(hi_geom_10_off + 1203 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxx_xxxxxx, g_z_0_xxxx_xxxxxxx, g_z_0_xxxx_xxxxxxy, g_z_0_xxxx_xxxxxxz, g_z_0_xxxx_xxxxxy, g_z_0_xxxx_xxxxxyy, g_z_0_xxxx_xxxxxyz, g_z_0_xxxx_xxxxxz, g_z_0_xxxx_xxxxxzz, g_z_0_xxxx_xxxxyy, g_z_0_xxxx_xxxxyyy, g_z_0_xxxx_xxxxyyz, g_z_0_xxxx_xxxxyz, g_z_0_xxxx_xxxxyzz, g_z_0_xxxx_xxxxzz, g_z_0_xxxx_xxxxzzz, g_z_0_xxxx_xxxyyy, g_z_0_xxxx_xxxyyyy, g_z_0_xxxx_xxxyyyz, g_z_0_xxxx_xxxyyz, g_z_0_xxxx_xxxyyzz, g_z_0_xxxx_xxxyzz, g_z_0_xxxx_xxxyzzz, g_z_0_xxxx_xxxzzz, g_z_0_xxxx_xxxzzzz, g_z_0_xxxx_xxyyyy, g_z_0_xxxx_xxyyyyy, g_z_0_xxxx_xxyyyyz, g_z_0_xxxx_xxyyyz, g_z_0_xxxx_xxyyyzz, g_z_0_xxxx_xxyyzz, g_z_0_xxxx_xxyyzzz, g_z_0_xxxx_xxyzzz, g_z_0_xxxx_xxyzzzz, g_z_0_xxxx_xxzzzz, g_z_0_xxxx_xxzzzzz, g_z_0_xxxx_xyyyyy, g_z_0_xxxx_xyyyyyy, g_z_0_xxxx_xyyyyyz, g_z_0_xxxx_xyyyyz, g_z_0_xxxx_xyyyyzz, g_z_0_xxxx_xyyyzz, g_z_0_xxxx_xyyyzzz, g_z_0_xxxx_xyyzzz, g_z_0_xxxx_xyyzzzz, g_z_0_xxxx_xyzzzz, g_z_0_xxxx_xyzzzzz, g_z_0_xxxx_xzzzzz, g_z_0_xxxx_xzzzzzz, g_z_0_xxxx_yyyyyy, g_z_0_xxxx_yyyyyz, g_z_0_xxxx_yyyyzz, g_z_0_xxxx_yyyzzz, g_z_0_xxxx_yyzzzz, g_z_0_xxxx_yzzzzz, g_z_0_xxxx_zzzzzz, g_z_0_xxxxx_xxxxxx, g_z_0_xxxxx_xxxxxy, g_z_0_xxxxx_xxxxxz, g_z_0_xxxxx_xxxxyy, g_z_0_xxxxx_xxxxyz, g_z_0_xxxxx_xxxxzz, g_z_0_xxxxx_xxxyyy, g_z_0_xxxxx_xxxyyz, g_z_0_xxxxx_xxxyzz, g_z_0_xxxxx_xxxzzz, g_z_0_xxxxx_xxyyyy, g_z_0_xxxxx_xxyyyz, g_z_0_xxxxx_xxyyzz, g_z_0_xxxxx_xxyzzz, g_z_0_xxxxx_xxzzzz, g_z_0_xxxxx_xyyyyy, g_z_0_xxxxx_xyyyyz, g_z_0_xxxxx_xyyyzz, g_z_0_xxxxx_xyyzzz, g_z_0_xxxxx_xyzzzz, g_z_0_xxxxx_xzzzzz, g_z_0_xxxxx_yyyyyy, g_z_0_xxxxx_yyyyyz, g_z_0_xxxxx_yyyyzz, g_z_0_xxxxx_yyyzzz, g_z_0_xxxxx_yyzzzz, g_z_0_xxxxx_yzzzzz, g_z_0_xxxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxx_xxxxxx[k] = -g_z_0_xxxx_xxxxxx[k] * ab_x + g_z_0_xxxx_xxxxxxx[k];

                g_z_0_xxxxx_xxxxxy[k] = -g_z_0_xxxx_xxxxxy[k] * ab_x + g_z_0_xxxx_xxxxxxy[k];

                g_z_0_xxxxx_xxxxxz[k] = -g_z_0_xxxx_xxxxxz[k] * ab_x + g_z_0_xxxx_xxxxxxz[k];

                g_z_0_xxxxx_xxxxyy[k] = -g_z_0_xxxx_xxxxyy[k] * ab_x + g_z_0_xxxx_xxxxxyy[k];

                g_z_0_xxxxx_xxxxyz[k] = -g_z_0_xxxx_xxxxyz[k] * ab_x + g_z_0_xxxx_xxxxxyz[k];

                g_z_0_xxxxx_xxxxzz[k] = -g_z_0_xxxx_xxxxzz[k] * ab_x + g_z_0_xxxx_xxxxxzz[k];

                g_z_0_xxxxx_xxxyyy[k] = -g_z_0_xxxx_xxxyyy[k] * ab_x + g_z_0_xxxx_xxxxyyy[k];

                g_z_0_xxxxx_xxxyyz[k] = -g_z_0_xxxx_xxxyyz[k] * ab_x + g_z_0_xxxx_xxxxyyz[k];

                g_z_0_xxxxx_xxxyzz[k] = -g_z_0_xxxx_xxxyzz[k] * ab_x + g_z_0_xxxx_xxxxyzz[k];

                g_z_0_xxxxx_xxxzzz[k] = -g_z_0_xxxx_xxxzzz[k] * ab_x + g_z_0_xxxx_xxxxzzz[k];

                g_z_0_xxxxx_xxyyyy[k] = -g_z_0_xxxx_xxyyyy[k] * ab_x + g_z_0_xxxx_xxxyyyy[k];

                g_z_0_xxxxx_xxyyyz[k] = -g_z_0_xxxx_xxyyyz[k] * ab_x + g_z_0_xxxx_xxxyyyz[k];

                g_z_0_xxxxx_xxyyzz[k] = -g_z_0_xxxx_xxyyzz[k] * ab_x + g_z_0_xxxx_xxxyyzz[k];

                g_z_0_xxxxx_xxyzzz[k] = -g_z_0_xxxx_xxyzzz[k] * ab_x + g_z_0_xxxx_xxxyzzz[k];

                g_z_0_xxxxx_xxzzzz[k] = -g_z_0_xxxx_xxzzzz[k] * ab_x + g_z_0_xxxx_xxxzzzz[k];

                g_z_0_xxxxx_xyyyyy[k] = -g_z_0_xxxx_xyyyyy[k] * ab_x + g_z_0_xxxx_xxyyyyy[k];

                g_z_0_xxxxx_xyyyyz[k] = -g_z_0_xxxx_xyyyyz[k] * ab_x + g_z_0_xxxx_xxyyyyz[k];

                g_z_0_xxxxx_xyyyzz[k] = -g_z_0_xxxx_xyyyzz[k] * ab_x + g_z_0_xxxx_xxyyyzz[k];

                g_z_0_xxxxx_xyyzzz[k] = -g_z_0_xxxx_xyyzzz[k] * ab_x + g_z_0_xxxx_xxyyzzz[k];

                g_z_0_xxxxx_xyzzzz[k] = -g_z_0_xxxx_xyzzzz[k] * ab_x + g_z_0_xxxx_xxyzzzz[k];

                g_z_0_xxxxx_xzzzzz[k] = -g_z_0_xxxx_xzzzzz[k] * ab_x + g_z_0_xxxx_xxzzzzz[k];

                g_z_0_xxxxx_yyyyyy[k] = -g_z_0_xxxx_yyyyyy[k] * ab_x + g_z_0_xxxx_xyyyyyy[k];

                g_z_0_xxxxx_yyyyyz[k] = -g_z_0_xxxx_yyyyyz[k] * ab_x + g_z_0_xxxx_xyyyyyz[k];

                g_z_0_xxxxx_yyyyzz[k] = -g_z_0_xxxx_yyyyzz[k] * ab_x + g_z_0_xxxx_xyyyyzz[k];

                g_z_0_xxxxx_yyyzzz[k] = -g_z_0_xxxx_yyyzzz[k] * ab_x + g_z_0_xxxx_xyyyzzz[k];

                g_z_0_xxxxx_yyzzzz[k] = -g_z_0_xxxx_yyzzzz[k] * ab_x + g_z_0_xxxx_xyyzzzz[k];

                g_z_0_xxxxx_yzzzzz[k] = -g_z_0_xxxx_yzzzzz[k] * ab_x + g_z_0_xxxx_xyzzzzz[k];

                g_z_0_xxxxx_zzzzzz[k] = -g_z_0_xxxx_zzzzzz[k] * ab_x + g_z_0_xxxx_xzzzzzz[k];
            }

            /// Set up 1204-1232 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxy_xxxxxx = cbuffer.data(hi_geom_10_off + 1204 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxy = cbuffer.data(hi_geom_10_off + 1205 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxz = cbuffer.data(hi_geom_10_off + 1206 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxyy = cbuffer.data(hi_geom_10_off + 1207 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxyz = cbuffer.data(hi_geom_10_off + 1208 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxzz = cbuffer.data(hi_geom_10_off + 1209 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyyy = cbuffer.data(hi_geom_10_off + 1210 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyyz = cbuffer.data(hi_geom_10_off + 1211 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyzz = cbuffer.data(hi_geom_10_off + 1212 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxzzz = cbuffer.data(hi_geom_10_off + 1213 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyyy = cbuffer.data(hi_geom_10_off + 1214 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyyz = cbuffer.data(hi_geom_10_off + 1215 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyzz = cbuffer.data(hi_geom_10_off + 1216 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyzzz = cbuffer.data(hi_geom_10_off + 1217 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxzzzz = cbuffer.data(hi_geom_10_off + 1218 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyyy = cbuffer.data(hi_geom_10_off + 1219 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyyz = cbuffer.data(hi_geom_10_off + 1220 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyzz = cbuffer.data(hi_geom_10_off + 1221 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyzzz = cbuffer.data(hi_geom_10_off + 1222 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyzzzz = cbuffer.data(hi_geom_10_off + 1223 * ccomps * dcomps);

            auto g_z_0_xxxxy_xzzzzz = cbuffer.data(hi_geom_10_off + 1224 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyyyy = cbuffer.data(hi_geom_10_off + 1225 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyyyz = cbuffer.data(hi_geom_10_off + 1226 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyyzz = cbuffer.data(hi_geom_10_off + 1227 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyzzz = cbuffer.data(hi_geom_10_off + 1228 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyzzzz = cbuffer.data(hi_geom_10_off + 1229 * ccomps * dcomps);

            auto g_z_0_xxxxy_yzzzzz = cbuffer.data(hi_geom_10_off + 1230 * ccomps * dcomps);

            auto g_z_0_xxxxy_zzzzzz = cbuffer.data(hi_geom_10_off + 1231 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxy_xxxxxx, g_z_0_xxxxy_xxxxxy, g_z_0_xxxxy_xxxxxz, g_z_0_xxxxy_xxxxyy, g_z_0_xxxxy_xxxxyz, g_z_0_xxxxy_xxxxzz, g_z_0_xxxxy_xxxyyy, g_z_0_xxxxy_xxxyyz, g_z_0_xxxxy_xxxyzz, g_z_0_xxxxy_xxxzzz, g_z_0_xxxxy_xxyyyy, g_z_0_xxxxy_xxyyyz, g_z_0_xxxxy_xxyyzz, g_z_0_xxxxy_xxyzzz, g_z_0_xxxxy_xxzzzz, g_z_0_xxxxy_xyyyyy, g_z_0_xxxxy_xyyyyz, g_z_0_xxxxy_xyyyzz, g_z_0_xxxxy_xyyzzz, g_z_0_xxxxy_xyzzzz, g_z_0_xxxxy_xzzzzz, g_z_0_xxxxy_yyyyyy, g_z_0_xxxxy_yyyyyz, g_z_0_xxxxy_yyyyzz, g_z_0_xxxxy_yyyzzz, g_z_0_xxxxy_yyzzzz, g_z_0_xxxxy_yzzzzz, g_z_0_xxxxy_zzzzzz, g_z_0_xxxy_xxxxxx, g_z_0_xxxy_xxxxxxx, g_z_0_xxxy_xxxxxxy, g_z_0_xxxy_xxxxxxz, g_z_0_xxxy_xxxxxy, g_z_0_xxxy_xxxxxyy, g_z_0_xxxy_xxxxxyz, g_z_0_xxxy_xxxxxz, g_z_0_xxxy_xxxxxzz, g_z_0_xxxy_xxxxyy, g_z_0_xxxy_xxxxyyy, g_z_0_xxxy_xxxxyyz, g_z_0_xxxy_xxxxyz, g_z_0_xxxy_xxxxyzz, g_z_0_xxxy_xxxxzz, g_z_0_xxxy_xxxxzzz, g_z_0_xxxy_xxxyyy, g_z_0_xxxy_xxxyyyy, g_z_0_xxxy_xxxyyyz, g_z_0_xxxy_xxxyyz, g_z_0_xxxy_xxxyyzz, g_z_0_xxxy_xxxyzz, g_z_0_xxxy_xxxyzzz, g_z_0_xxxy_xxxzzz, g_z_0_xxxy_xxxzzzz, g_z_0_xxxy_xxyyyy, g_z_0_xxxy_xxyyyyy, g_z_0_xxxy_xxyyyyz, g_z_0_xxxy_xxyyyz, g_z_0_xxxy_xxyyyzz, g_z_0_xxxy_xxyyzz, g_z_0_xxxy_xxyyzzz, g_z_0_xxxy_xxyzzz, g_z_0_xxxy_xxyzzzz, g_z_0_xxxy_xxzzzz, g_z_0_xxxy_xxzzzzz, g_z_0_xxxy_xyyyyy, g_z_0_xxxy_xyyyyyy, g_z_0_xxxy_xyyyyyz, g_z_0_xxxy_xyyyyz, g_z_0_xxxy_xyyyyzz, g_z_0_xxxy_xyyyzz, g_z_0_xxxy_xyyyzzz, g_z_0_xxxy_xyyzzz, g_z_0_xxxy_xyyzzzz, g_z_0_xxxy_xyzzzz, g_z_0_xxxy_xyzzzzz, g_z_0_xxxy_xzzzzz, g_z_0_xxxy_xzzzzzz, g_z_0_xxxy_yyyyyy, g_z_0_xxxy_yyyyyz, g_z_0_xxxy_yyyyzz, g_z_0_xxxy_yyyzzz, g_z_0_xxxy_yyzzzz, g_z_0_xxxy_yzzzzz, g_z_0_xxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxy_xxxxxx[k] = -g_z_0_xxxy_xxxxxx[k] * ab_x + g_z_0_xxxy_xxxxxxx[k];

                g_z_0_xxxxy_xxxxxy[k] = -g_z_0_xxxy_xxxxxy[k] * ab_x + g_z_0_xxxy_xxxxxxy[k];

                g_z_0_xxxxy_xxxxxz[k] = -g_z_0_xxxy_xxxxxz[k] * ab_x + g_z_0_xxxy_xxxxxxz[k];

                g_z_0_xxxxy_xxxxyy[k] = -g_z_0_xxxy_xxxxyy[k] * ab_x + g_z_0_xxxy_xxxxxyy[k];

                g_z_0_xxxxy_xxxxyz[k] = -g_z_0_xxxy_xxxxyz[k] * ab_x + g_z_0_xxxy_xxxxxyz[k];

                g_z_0_xxxxy_xxxxzz[k] = -g_z_0_xxxy_xxxxzz[k] * ab_x + g_z_0_xxxy_xxxxxzz[k];

                g_z_0_xxxxy_xxxyyy[k] = -g_z_0_xxxy_xxxyyy[k] * ab_x + g_z_0_xxxy_xxxxyyy[k];

                g_z_0_xxxxy_xxxyyz[k] = -g_z_0_xxxy_xxxyyz[k] * ab_x + g_z_0_xxxy_xxxxyyz[k];

                g_z_0_xxxxy_xxxyzz[k] = -g_z_0_xxxy_xxxyzz[k] * ab_x + g_z_0_xxxy_xxxxyzz[k];

                g_z_0_xxxxy_xxxzzz[k] = -g_z_0_xxxy_xxxzzz[k] * ab_x + g_z_0_xxxy_xxxxzzz[k];

                g_z_0_xxxxy_xxyyyy[k] = -g_z_0_xxxy_xxyyyy[k] * ab_x + g_z_0_xxxy_xxxyyyy[k];

                g_z_0_xxxxy_xxyyyz[k] = -g_z_0_xxxy_xxyyyz[k] * ab_x + g_z_0_xxxy_xxxyyyz[k];

                g_z_0_xxxxy_xxyyzz[k] = -g_z_0_xxxy_xxyyzz[k] * ab_x + g_z_0_xxxy_xxxyyzz[k];

                g_z_0_xxxxy_xxyzzz[k] = -g_z_0_xxxy_xxyzzz[k] * ab_x + g_z_0_xxxy_xxxyzzz[k];

                g_z_0_xxxxy_xxzzzz[k] = -g_z_0_xxxy_xxzzzz[k] * ab_x + g_z_0_xxxy_xxxzzzz[k];

                g_z_0_xxxxy_xyyyyy[k] = -g_z_0_xxxy_xyyyyy[k] * ab_x + g_z_0_xxxy_xxyyyyy[k];

                g_z_0_xxxxy_xyyyyz[k] = -g_z_0_xxxy_xyyyyz[k] * ab_x + g_z_0_xxxy_xxyyyyz[k];

                g_z_0_xxxxy_xyyyzz[k] = -g_z_0_xxxy_xyyyzz[k] * ab_x + g_z_0_xxxy_xxyyyzz[k];

                g_z_0_xxxxy_xyyzzz[k] = -g_z_0_xxxy_xyyzzz[k] * ab_x + g_z_0_xxxy_xxyyzzz[k];

                g_z_0_xxxxy_xyzzzz[k] = -g_z_0_xxxy_xyzzzz[k] * ab_x + g_z_0_xxxy_xxyzzzz[k];

                g_z_0_xxxxy_xzzzzz[k] = -g_z_0_xxxy_xzzzzz[k] * ab_x + g_z_0_xxxy_xxzzzzz[k];

                g_z_0_xxxxy_yyyyyy[k] = -g_z_0_xxxy_yyyyyy[k] * ab_x + g_z_0_xxxy_xyyyyyy[k];

                g_z_0_xxxxy_yyyyyz[k] = -g_z_0_xxxy_yyyyyz[k] * ab_x + g_z_0_xxxy_xyyyyyz[k];

                g_z_0_xxxxy_yyyyzz[k] = -g_z_0_xxxy_yyyyzz[k] * ab_x + g_z_0_xxxy_xyyyyzz[k];

                g_z_0_xxxxy_yyyzzz[k] = -g_z_0_xxxy_yyyzzz[k] * ab_x + g_z_0_xxxy_xyyyzzz[k];

                g_z_0_xxxxy_yyzzzz[k] = -g_z_0_xxxy_yyzzzz[k] * ab_x + g_z_0_xxxy_xyyzzzz[k];

                g_z_0_xxxxy_yzzzzz[k] = -g_z_0_xxxy_yzzzzz[k] * ab_x + g_z_0_xxxy_xyzzzzz[k];

                g_z_0_xxxxy_zzzzzz[k] = -g_z_0_xxxy_zzzzzz[k] * ab_x + g_z_0_xxxy_xzzzzzz[k];
            }

            /// Set up 1232-1260 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxz_xxxxxx = cbuffer.data(hi_geom_10_off + 1232 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxy = cbuffer.data(hi_geom_10_off + 1233 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxz = cbuffer.data(hi_geom_10_off + 1234 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxyy = cbuffer.data(hi_geom_10_off + 1235 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxyz = cbuffer.data(hi_geom_10_off + 1236 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxzz = cbuffer.data(hi_geom_10_off + 1237 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyyy = cbuffer.data(hi_geom_10_off + 1238 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyyz = cbuffer.data(hi_geom_10_off + 1239 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyzz = cbuffer.data(hi_geom_10_off + 1240 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxzzz = cbuffer.data(hi_geom_10_off + 1241 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyyy = cbuffer.data(hi_geom_10_off + 1242 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyyz = cbuffer.data(hi_geom_10_off + 1243 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyzz = cbuffer.data(hi_geom_10_off + 1244 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyzzz = cbuffer.data(hi_geom_10_off + 1245 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxzzzz = cbuffer.data(hi_geom_10_off + 1246 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyyy = cbuffer.data(hi_geom_10_off + 1247 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyyz = cbuffer.data(hi_geom_10_off + 1248 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyzz = cbuffer.data(hi_geom_10_off + 1249 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyzzz = cbuffer.data(hi_geom_10_off + 1250 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyzzzz = cbuffer.data(hi_geom_10_off + 1251 * ccomps * dcomps);

            auto g_z_0_xxxxz_xzzzzz = cbuffer.data(hi_geom_10_off + 1252 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyyyy = cbuffer.data(hi_geom_10_off + 1253 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyyyz = cbuffer.data(hi_geom_10_off + 1254 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyyzz = cbuffer.data(hi_geom_10_off + 1255 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyzzz = cbuffer.data(hi_geom_10_off + 1256 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyzzzz = cbuffer.data(hi_geom_10_off + 1257 * ccomps * dcomps);

            auto g_z_0_xxxxz_yzzzzz = cbuffer.data(hi_geom_10_off + 1258 * ccomps * dcomps);

            auto g_z_0_xxxxz_zzzzzz = cbuffer.data(hi_geom_10_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxz_xxxxxx, g_z_0_xxxxz_xxxxxy, g_z_0_xxxxz_xxxxxz, g_z_0_xxxxz_xxxxyy, g_z_0_xxxxz_xxxxyz, g_z_0_xxxxz_xxxxzz, g_z_0_xxxxz_xxxyyy, g_z_0_xxxxz_xxxyyz, g_z_0_xxxxz_xxxyzz, g_z_0_xxxxz_xxxzzz, g_z_0_xxxxz_xxyyyy, g_z_0_xxxxz_xxyyyz, g_z_0_xxxxz_xxyyzz, g_z_0_xxxxz_xxyzzz, g_z_0_xxxxz_xxzzzz, g_z_0_xxxxz_xyyyyy, g_z_0_xxxxz_xyyyyz, g_z_0_xxxxz_xyyyzz, g_z_0_xxxxz_xyyzzz, g_z_0_xxxxz_xyzzzz, g_z_0_xxxxz_xzzzzz, g_z_0_xxxxz_yyyyyy, g_z_0_xxxxz_yyyyyz, g_z_0_xxxxz_yyyyzz, g_z_0_xxxxz_yyyzzz, g_z_0_xxxxz_yyzzzz, g_z_0_xxxxz_yzzzzz, g_z_0_xxxxz_zzzzzz, g_z_0_xxxz_xxxxxx, g_z_0_xxxz_xxxxxxx, g_z_0_xxxz_xxxxxxy, g_z_0_xxxz_xxxxxxz, g_z_0_xxxz_xxxxxy, g_z_0_xxxz_xxxxxyy, g_z_0_xxxz_xxxxxyz, g_z_0_xxxz_xxxxxz, g_z_0_xxxz_xxxxxzz, g_z_0_xxxz_xxxxyy, g_z_0_xxxz_xxxxyyy, g_z_0_xxxz_xxxxyyz, g_z_0_xxxz_xxxxyz, g_z_0_xxxz_xxxxyzz, g_z_0_xxxz_xxxxzz, g_z_0_xxxz_xxxxzzz, g_z_0_xxxz_xxxyyy, g_z_0_xxxz_xxxyyyy, g_z_0_xxxz_xxxyyyz, g_z_0_xxxz_xxxyyz, g_z_0_xxxz_xxxyyzz, g_z_0_xxxz_xxxyzz, g_z_0_xxxz_xxxyzzz, g_z_0_xxxz_xxxzzz, g_z_0_xxxz_xxxzzzz, g_z_0_xxxz_xxyyyy, g_z_0_xxxz_xxyyyyy, g_z_0_xxxz_xxyyyyz, g_z_0_xxxz_xxyyyz, g_z_0_xxxz_xxyyyzz, g_z_0_xxxz_xxyyzz, g_z_0_xxxz_xxyyzzz, g_z_0_xxxz_xxyzzz, g_z_0_xxxz_xxyzzzz, g_z_0_xxxz_xxzzzz, g_z_0_xxxz_xxzzzzz, g_z_0_xxxz_xyyyyy, g_z_0_xxxz_xyyyyyy, g_z_0_xxxz_xyyyyyz, g_z_0_xxxz_xyyyyz, g_z_0_xxxz_xyyyyzz, g_z_0_xxxz_xyyyzz, g_z_0_xxxz_xyyyzzz, g_z_0_xxxz_xyyzzz, g_z_0_xxxz_xyyzzzz, g_z_0_xxxz_xyzzzz, g_z_0_xxxz_xyzzzzz, g_z_0_xxxz_xzzzzz, g_z_0_xxxz_xzzzzzz, g_z_0_xxxz_yyyyyy, g_z_0_xxxz_yyyyyz, g_z_0_xxxz_yyyyzz, g_z_0_xxxz_yyyzzz, g_z_0_xxxz_yyzzzz, g_z_0_xxxz_yzzzzz, g_z_0_xxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxz_xxxxxx[k] = -g_z_0_xxxz_xxxxxx[k] * ab_x + g_z_0_xxxz_xxxxxxx[k];

                g_z_0_xxxxz_xxxxxy[k] = -g_z_0_xxxz_xxxxxy[k] * ab_x + g_z_0_xxxz_xxxxxxy[k];

                g_z_0_xxxxz_xxxxxz[k] = -g_z_0_xxxz_xxxxxz[k] * ab_x + g_z_0_xxxz_xxxxxxz[k];

                g_z_0_xxxxz_xxxxyy[k] = -g_z_0_xxxz_xxxxyy[k] * ab_x + g_z_0_xxxz_xxxxxyy[k];

                g_z_0_xxxxz_xxxxyz[k] = -g_z_0_xxxz_xxxxyz[k] * ab_x + g_z_0_xxxz_xxxxxyz[k];

                g_z_0_xxxxz_xxxxzz[k] = -g_z_0_xxxz_xxxxzz[k] * ab_x + g_z_0_xxxz_xxxxxzz[k];

                g_z_0_xxxxz_xxxyyy[k] = -g_z_0_xxxz_xxxyyy[k] * ab_x + g_z_0_xxxz_xxxxyyy[k];

                g_z_0_xxxxz_xxxyyz[k] = -g_z_0_xxxz_xxxyyz[k] * ab_x + g_z_0_xxxz_xxxxyyz[k];

                g_z_0_xxxxz_xxxyzz[k] = -g_z_0_xxxz_xxxyzz[k] * ab_x + g_z_0_xxxz_xxxxyzz[k];

                g_z_0_xxxxz_xxxzzz[k] = -g_z_0_xxxz_xxxzzz[k] * ab_x + g_z_0_xxxz_xxxxzzz[k];

                g_z_0_xxxxz_xxyyyy[k] = -g_z_0_xxxz_xxyyyy[k] * ab_x + g_z_0_xxxz_xxxyyyy[k];

                g_z_0_xxxxz_xxyyyz[k] = -g_z_0_xxxz_xxyyyz[k] * ab_x + g_z_0_xxxz_xxxyyyz[k];

                g_z_0_xxxxz_xxyyzz[k] = -g_z_0_xxxz_xxyyzz[k] * ab_x + g_z_0_xxxz_xxxyyzz[k];

                g_z_0_xxxxz_xxyzzz[k] = -g_z_0_xxxz_xxyzzz[k] * ab_x + g_z_0_xxxz_xxxyzzz[k];

                g_z_0_xxxxz_xxzzzz[k] = -g_z_0_xxxz_xxzzzz[k] * ab_x + g_z_0_xxxz_xxxzzzz[k];

                g_z_0_xxxxz_xyyyyy[k] = -g_z_0_xxxz_xyyyyy[k] * ab_x + g_z_0_xxxz_xxyyyyy[k];

                g_z_0_xxxxz_xyyyyz[k] = -g_z_0_xxxz_xyyyyz[k] * ab_x + g_z_0_xxxz_xxyyyyz[k];

                g_z_0_xxxxz_xyyyzz[k] = -g_z_0_xxxz_xyyyzz[k] * ab_x + g_z_0_xxxz_xxyyyzz[k];

                g_z_0_xxxxz_xyyzzz[k] = -g_z_0_xxxz_xyyzzz[k] * ab_x + g_z_0_xxxz_xxyyzzz[k];

                g_z_0_xxxxz_xyzzzz[k] = -g_z_0_xxxz_xyzzzz[k] * ab_x + g_z_0_xxxz_xxyzzzz[k];

                g_z_0_xxxxz_xzzzzz[k] = -g_z_0_xxxz_xzzzzz[k] * ab_x + g_z_0_xxxz_xxzzzzz[k];

                g_z_0_xxxxz_yyyyyy[k] = -g_z_0_xxxz_yyyyyy[k] * ab_x + g_z_0_xxxz_xyyyyyy[k];

                g_z_0_xxxxz_yyyyyz[k] = -g_z_0_xxxz_yyyyyz[k] * ab_x + g_z_0_xxxz_xyyyyyz[k];

                g_z_0_xxxxz_yyyyzz[k] = -g_z_0_xxxz_yyyyzz[k] * ab_x + g_z_0_xxxz_xyyyyzz[k];

                g_z_0_xxxxz_yyyzzz[k] = -g_z_0_xxxz_yyyzzz[k] * ab_x + g_z_0_xxxz_xyyyzzz[k];

                g_z_0_xxxxz_yyzzzz[k] = -g_z_0_xxxz_yyzzzz[k] * ab_x + g_z_0_xxxz_xyyzzzz[k];

                g_z_0_xxxxz_yzzzzz[k] = -g_z_0_xxxz_yzzzzz[k] * ab_x + g_z_0_xxxz_xyzzzzz[k];

                g_z_0_xxxxz_zzzzzz[k] = -g_z_0_xxxz_zzzzzz[k] * ab_x + g_z_0_xxxz_xzzzzzz[k];
            }

            /// Set up 1260-1288 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1260 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1261 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1262 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1263 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1264 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1265 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1266 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1267 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1268 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1269 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1270 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1271 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1272 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1273 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1274 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1275 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1276 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1277 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1278 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1279 * ccomps * dcomps);

            auto g_z_0_xxxyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1280 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1281 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1282 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1283 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1284 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1285 * ccomps * dcomps);

            auto g_z_0_xxxyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1286 * ccomps * dcomps);

            auto g_z_0_xxxyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1287 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyy_xxxxxx, g_z_0_xxxyy_xxxxxy, g_z_0_xxxyy_xxxxxz, g_z_0_xxxyy_xxxxyy, g_z_0_xxxyy_xxxxyz, g_z_0_xxxyy_xxxxzz, g_z_0_xxxyy_xxxyyy, g_z_0_xxxyy_xxxyyz, g_z_0_xxxyy_xxxyzz, g_z_0_xxxyy_xxxzzz, g_z_0_xxxyy_xxyyyy, g_z_0_xxxyy_xxyyyz, g_z_0_xxxyy_xxyyzz, g_z_0_xxxyy_xxyzzz, g_z_0_xxxyy_xxzzzz, g_z_0_xxxyy_xyyyyy, g_z_0_xxxyy_xyyyyz, g_z_0_xxxyy_xyyyzz, g_z_0_xxxyy_xyyzzz, g_z_0_xxxyy_xyzzzz, g_z_0_xxxyy_xzzzzz, g_z_0_xxxyy_yyyyyy, g_z_0_xxxyy_yyyyyz, g_z_0_xxxyy_yyyyzz, g_z_0_xxxyy_yyyzzz, g_z_0_xxxyy_yyzzzz, g_z_0_xxxyy_yzzzzz, g_z_0_xxxyy_zzzzzz, g_z_0_xxyy_xxxxxx, g_z_0_xxyy_xxxxxxx, g_z_0_xxyy_xxxxxxy, g_z_0_xxyy_xxxxxxz, g_z_0_xxyy_xxxxxy, g_z_0_xxyy_xxxxxyy, g_z_0_xxyy_xxxxxyz, g_z_0_xxyy_xxxxxz, g_z_0_xxyy_xxxxxzz, g_z_0_xxyy_xxxxyy, g_z_0_xxyy_xxxxyyy, g_z_0_xxyy_xxxxyyz, g_z_0_xxyy_xxxxyz, g_z_0_xxyy_xxxxyzz, g_z_0_xxyy_xxxxzz, g_z_0_xxyy_xxxxzzz, g_z_0_xxyy_xxxyyy, g_z_0_xxyy_xxxyyyy, g_z_0_xxyy_xxxyyyz, g_z_0_xxyy_xxxyyz, g_z_0_xxyy_xxxyyzz, g_z_0_xxyy_xxxyzz, g_z_0_xxyy_xxxyzzz, g_z_0_xxyy_xxxzzz, g_z_0_xxyy_xxxzzzz, g_z_0_xxyy_xxyyyy, g_z_0_xxyy_xxyyyyy, g_z_0_xxyy_xxyyyyz, g_z_0_xxyy_xxyyyz, g_z_0_xxyy_xxyyyzz, g_z_0_xxyy_xxyyzz, g_z_0_xxyy_xxyyzzz, g_z_0_xxyy_xxyzzz, g_z_0_xxyy_xxyzzzz, g_z_0_xxyy_xxzzzz, g_z_0_xxyy_xxzzzzz, g_z_0_xxyy_xyyyyy, g_z_0_xxyy_xyyyyyy, g_z_0_xxyy_xyyyyyz, g_z_0_xxyy_xyyyyz, g_z_0_xxyy_xyyyyzz, g_z_0_xxyy_xyyyzz, g_z_0_xxyy_xyyyzzz, g_z_0_xxyy_xyyzzz, g_z_0_xxyy_xyyzzzz, g_z_0_xxyy_xyzzzz, g_z_0_xxyy_xyzzzzz, g_z_0_xxyy_xzzzzz, g_z_0_xxyy_xzzzzzz, g_z_0_xxyy_yyyyyy, g_z_0_xxyy_yyyyyz, g_z_0_xxyy_yyyyzz, g_z_0_xxyy_yyyzzz, g_z_0_xxyy_yyzzzz, g_z_0_xxyy_yzzzzz, g_z_0_xxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyy_xxxxxx[k] = -g_z_0_xxyy_xxxxxx[k] * ab_x + g_z_0_xxyy_xxxxxxx[k];

                g_z_0_xxxyy_xxxxxy[k] = -g_z_0_xxyy_xxxxxy[k] * ab_x + g_z_0_xxyy_xxxxxxy[k];

                g_z_0_xxxyy_xxxxxz[k] = -g_z_0_xxyy_xxxxxz[k] * ab_x + g_z_0_xxyy_xxxxxxz[k];

                g_z_0_xxxyy_xxxxyy[k] = -g_z_0_xxyy_xxxxyy[k] * ab_x + g_z_0_xxyy_xxxxxyy[k];

                g_z_0_xxxyy_xxxxyz[k] = -g_z_0_xxyy_xxxxyz[k] * ab_x + g_z_0_xxyy_xxxxxyz[k];

                g_z_0_xxxyy_xxxxzz[k] = -g_z_0_xxyy_xxxxzz[k] * ab_x + g_z_0_xxyy_xxxxxzz[k];

                g_z_0_xxxyy_xxxyyy[k] = -g_z_0_xxyy_xxxyyy[k] * ab_x + g_z_0_xxyy_xxxxyyy[k];

                g_z_0_xxxyy_xxxyyz[k] = -g_z_0_xxyy_xxxyyz[k] * ab_x + g_z_0_xxyy_xxxxyyz[k];

                g_z_0_xxxyy_xxxyzz[k] = -g_z_0_xxyy_xxxyzz[k] * ab_x + g_z_0_xxyy_xxxxyzz[k];

                g_z_0_xxxyy_xxxzzz[k] = -g_z_0_xxyy_xxxzzz[k] * ab_x + g_z_0_xxyy_xxxxzzz[k];

                g_z_0_xxxyy_xxyyyy[k] = -g_z_0_xxyy_xxyyyy[k] * ab_x + g_z_0_xxyy_xxxyyyy[k];

                g_z_0_xxxyy_xxyyyz[k] = -g_z_0_xxyy_xxyyyz[k] * ab_x + g_z_0_xxyy_xxxyyyz[k];

                g_z_0_xxxyy_xxyyzz[k] = -g_z_0_xxyy_xxyyzz[k] * ab_x + g_z_0_xxyy_xxxyyzz[k];

                g_z_0_xxxyy_xxyzzz[k] = -g_z_0_xxyy_xxyzzz[k] * ab_x + g_z_0_xxyy_xxxyzzz[k];

                g_z_0_xxxyy_xxzzzz[k] = -g_z_0_xxyy_xxzzzz[k] * ab_x + g_z_0_xxyy_xxxzzzz[k];

                g_z_0_xxxyy_xyyyyy[k] = -g_z_0_xxyy_xyyyyy[k] * ab_x + g_z_0_xxyy_xxyyyyy[k];

                g_z_0_xxxyy_xyyyyz[k] = -g_z_0_xxyy_xyyyyz[k] * ab_x + g_z_0_xxyy_xxyyyyz[k];

                g_z_0_xxxyy_xyyyzz[k] = -g_z_0_xxyy_xyyyzz[k] * ab_x + g_z_0_xxyy_xxyyyzz[k];

                g_z_0_xxxyy_xyyzzz[k] = -g_z_0_xxyy_xyyzzz[k] * ab_x + g_z_0_xxyy_xxyyzzz[k];

                g_z_0_xxxyy_xyzzzz[k] = -g_z_0_xxyy_xyzzzz[k] * ab_x + g_z_0_xxyy_xxyzzzz[k];

                g_z_0_xxxyy_xzzzzz[k] = -g_z_0_xxyy_xzzzzz[k] * ab_x + g_z_0_xxyy_xxzzzzz[k];

                g_z_0_xxxyy_yyyyyy[k] = -g_z_0_xxyy_yyyyyy[k] * ab_x + g_z_0_xxyy_xyyyyyy[k];

                g_z_0_xxxyy_yyyyyz[k] = -g_z_0_xxyy_yyyyyz[k] * ab_x + g_z_0_xxyy_xyyyyyz[k];

                g_z_0_xxxyy_yyyyzz[k] = -g_z_0_xxyy_yyyyzz[k] * ab_x + g_z_0_xxyy_xyyyyzz[k];

                g_z_0_xxxyy_yyyzzz[k] = -g_z_0_xxyy_yyyzzz[k] * ab_x + g_z_0_xxyy_xyyyzzz[k];

                g_z_0_xxxyy_yyzzzz[k] = -g_z_0_xxyy_yyzzzz[k] * ab_x + g_z_0_xxyy_xyyzzzz[k];

                g_z_0_xxxyy_yzzzzz[k] = -g_z_0_xxyy_yzzzzz[k] * ab_x + g_z_0_xxyy_xyzzzzz[k];

                g_z_0_xxxyy_zzzzzz[k] = -g_z_0_xxyy_zzzzzz[k] * ab_x + g_z_0_xxyy_xzzzzzz[k];
            }

            /// Set up 1288-1316 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1288 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1289 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1290 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1291 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1292 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1293 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1294 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1295 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1296 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1297 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1298 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1299 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1300 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1301 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1302 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1303 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1304 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1305 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1306 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1307 * ccomps * dcomps);

            auto g_z_0_xxxyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1308 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1309 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1310 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1311 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1312 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1313 * ccomps * dcomps);

            auto g_z_0_xxxyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1314 * ccomps * dcomps);

            auto g_z_0_xxxyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1315 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyz_xxxxxx, g_z_0_xxxyz_xxxxxy, g_z_0_xxxyz_xxxxxz, g_z_0_xxxyz_xxxxyy, g_z_0_xxxyz_xxxxyz, g_z_0_xxxyz_xxxxzz, g_z_0_xxxyz_xxxyyy, g_z_0_xxxyz_xxxyyz, g_z_0_xxxyz_xxxyzz, g_z_0_xxxyz_xxxzzz, g_z_0_xxxyz_xxyyyy, g_z_0_xxxyz_xxyyyz, g_z_0_xxxyz_xxyyzz, g_z_0_xxxyz_xxyzzz, g_z_0_xxxyz_xxzzzz, g_z_0_xxxyz_xyyyyy, g_z_0_xxxyz_xyyyyz, g_z_0_xxxyz_xyyyzz, g_z_0_xxxyz_xyyzzz, g_z_0_xxxyz_xyzzzz, g_z_0_xxxyz_xzzzzz, g_z_0_xxxyz_yyyyyy, g_z_0_xxxyz_yyyyyz, g_z_0_xxxyz_yyyyzz, g_z_0_xxxyz_yyyzzz, g_z_0_xxxyz_yyzzzz, g_z_0_xxxyz_yzzzzz, g_z_0_xxxyz_zzzzzz, g_z_0_xxyz_xxxxxx, g_z_0_xxyz_xxxxxxx, g_z_0_xxyz_xxxxxxy, g_z_0_xxyz_xxxxxxz, g_z_0_xxyz_xxxxxy, g_z_0_xxyz_xxxxxyy, g_z_0_xxyz_xxxxxyz, g_z_0_xxyz_xxxxxz, g_z_0_xxyz_xxxxxzz, g_z_0_xxyz_xxxxyy, g_z_0_xxyz_xxxxyyy, g_z_0_xxyz_xxxxyyz, g_z_0_xxyz_xxxxyz, g_z_0_xxyz_xxxxyzz, g_z_0_xxyz_xxxxzz, g_z_0_xxyz_xxxxzzz, g_z_0_xxyz_xxxyyy, g_z_0_xxyz_xxxyyyy, g_z_0_xxyz_xxxyyyz, g_z_0_xxyz_xxxyyz, g_z_0_xxyz_xxxyyzz, g_z_0_xxyz_xxxyzz, g_z_0_xxyz_xxxyzzz, g_z_0_xxyz_xxxzzz, g_z_0_xxyz_xxxzzzz, g_z_0_xxyz_xxyyyy, g_z_0_xxyz_xxyyyyy, g_z_0_xxyz_xxyyyyz, g_z_0_xxyz_xxyyyz, g_z_0_xxyz_xxyyyzz, g_z_0_xxyz_xxyyzz, g_z_0_xxyz_xxyyzzz, g_z_0_xxyz_xxyzzz, g_z_0_xxyz_xxyzzzz, g_z_0_xxyz_xxzzzz, g_z_0_xxyz_xxzzzzz, g_z_0_xxyz_xyyyyy, g_z_0_xxyz_xyyyyyy, g_z_0_xxyz_xyyyyyz, g_z_0_xxyz_xyyyyz, g_z_0_xxyz_xyyyyzz, g_z_0_xxyz_xyyyzz, g_z_0_xxyz_xyyyzzz, g_z_0_xxyz_xyyzzz, g_z_0_xxyz_xyyzzzz, g_z_0_xxyz_xyzzzz, g_z_0_xxyz_xyzzzzz, g_z_0_xxyz_xzzzzz, g_z_0_xxyz_xzzzzzz, g_z_0_xxyz_yyyyyy, g_z_0_xxyz_yyyyyz, g_z_0_xxyz_yyyyzz, g_z_0_xxyz_yyyzzz, g_z_0_xxyz_yyzzzz, g_z_0_xxyz_yzzzzz, g_z_0_xxyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyz_xxxxxx[k] = -g_z_0_xxyz_xxxxxx[k] * ab_x + g_z_0_xxyz_xxxxxxx[k];

                g_z_0_xxxyz_xxxxxy[k] = -g_z_0_xxyz_xxxxxy[k] * ab_x + g_z_0_xxyz_xxxxxxy[k];

                g_z_0_xxxyz_xxxxxz[k] = -g_z_0_xxyz_xxxxxz[k] * ab_x + g_z_0_xxyz_xxxxxxz[k];

                g_z_0_xxxyz_xxxxyy[k] = -g_z_0_xxyz_xxxxyy[k] * ab_x + g_z_0_xxyz_xxxxxyy[k];

                g_z_0_xxxyz_xxxxyz[k] = -g_z_0_xxyz_xxxxyz[k] * ab_x + g_z_0_xxyz_xxxxxyz[k];

                g_z_0_xxxyz_xxxxzz[k] = -g_z_0_xxyz_xxxxzz[k] * ab_x + g_z_0_xxyz_xxxxxzz[k];

                g_z_0_xxxyz_xxxyyy[k] = -g_z_0_xxyz_xxxyyy[k] * ab_x + g_z_0_xxyz_xxxxyyy[k];

                g_z_0_xxxyz_xxxyyz[k] = -g_z_0_xxyz_xxxyyz[k] * ab_x + g_z_0_xxyz_xxxxyyz[k];

                g_z_0_xxxyz_xxxyzz[k] = -g_z_0_xxyz_xxxyzz[k] * ab_x + g_z_0_xxyz_xxxxyzz[k];

                g_z_0_xxxyz_xxxzzz[k] = -g_z_0_xxyz_xxxzzz[k] * ab_x + g_z_0_xxyz_xxxxzzz[k];

                g_z_0_xxxyz_xxyyyy[k] = -g_z_0_xxyz_xxyyyy[k] * ab_x + g_z_0_xxyz_xxxyyyy[k];

                g_z_0_xxxyz_xxyyyz[k] = -g_z_0_xxyz_xxyyyz[k] * ab_x + g_z_0_xxyz_xxxyyyz[k];

                g_z_0_xxxyz_xxyyzz[k] = -g_z_0_xxyz_xxyyzz[k] * ab_x + g_z_0_xxyz_xxxyyzz[k];

                g_z_0_xxxyz_xxyzzz[k] = -g_z_0_xxyz_xxyzzz[k] * ab_x + g_z_0_xxyz_xxxyzzz[k];

                g_z_0_xxxyz_xxzzzz[k] = -g_z_0_xxyz_xxzzzz[k] * ab_x + g_z_0_xxyz_xxxzzzz[k];

                g_z_0_xxxyz_xyyyyy[k] = -g_z_0_xxyz_xyyyyy[k] * ab_x + g_z_0_xxyz_xxyyyyy[k];

                g_z_0_xxxyz_xyyyyz[k] = -g_z_0_xxyz_xyyyyz[k] * ab_x + g_z_0_xxyz_xxyyyyz[k];

                g_z_0_xxxyz_xyyyzz[k] = -g_z_0_xxyz_xyyyzz[k] * ab_x + g_z_0_xxyz_xxyyyzz[k];

                g_z_0_xxxyz_xyyzzz[k] = -g_z_0_xxyz_xyyzzz[k] * ab_x + g_z_0_xxyz_xxyyzzz[k];

                g_z_0_xxxyz_xyzzzz[k] = -g_z_0_xxyz_xyzzzz[k] * ab_x + g_z_0_xxyz_xxyzzzz[k];

                g_z_0_xxxyz_xzzzzz[k] = -g_z_0_xxyz_xzzzzz[k] * ab_x + g_z_0_xxyz_xxzzzzz[k];

                g_z_0_xxxyz_yyyyyy[k] = -g_z_0_xxyz_yyyyyy[k] * ab_x + g_z_0_xxyz_xyyyyyy[k];

                g_z_0_xxxyz_yyyyyz[k] = -g_z_0_xxyz_yyyyyz[k] * ab_x + g_z_0_xxyz_xyyyyyz[k];

                g_z_0_xxxyz_yyyyzz[k] = -g_z_0_xxyz_yyyyzz[k] * ab_x + g_z_0_xxyz_xyyyyzz[k];

                g_z_0_xxxyz_yyyzzz[k] = -g_z_0_xxyz_yyyzzz[k] * ab_x + g_z_0_xxyz_xyyyzzz[k];

                g_z_0_xxxyz_yyzzzz[k] = -g_z_0_xxyz_yyzzzz[k] * ab_x + g_z_0_xxyz_xyyzzzz[k];

                g_z_0_xxxyz_yzzzzz[k] = -g_z_0_xxyz_yzzzzz[k] * ab_x + g_z_0_xxyz_xyzzzzz[k];

                g_z_0_xxxyz_zzzzzz[k] = -g_z_0_xxyz_zzzzzz[k] * ab_x + g_z_0_xxyz_xzzzzzz[k];
            }

            /// Set up 1316-1344 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1316 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1317 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1318 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1319 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1320 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1321 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1322 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1323 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1324 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1325 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1326 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1327 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1328 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1329 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1330 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1331 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1332 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1333 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1334 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1335 * ccomps * dcomps);

            auto g_z_0_xxxzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1336 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1337 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1338 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1339 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1340 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1341 * ccomps * dcomps);

            auto g_z_0_xxxzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1342 * ccomps * dcomps);

            auto g_z_0_xxxzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1343 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxzz_xxxxxx, g_z_0_xxxzz_xxxxxy, g_z_0_xxxzz_xxxxxz, g_z_0_xxxzz_xxxxyy, g_z_0_xxxzz_xxxxyz, g_z_0_xxxzz_xxxxzz, g_z_0_xxxzz_xxxyyy, g_z_0_xxxzz_xxxyyz, g_z_0_xxxzz_xxxyzz, g_z_0_xxxzz_xxxzzz, g_z_0_xxxzz_xxyyyy, g_z_0_xxxzz_xxyyyz, g_z_0_xxxzz_xxyyzz, g_z_0_xxxzz_xxyzzz, g_z_0_xxxzz_xxzzzz, g_z_0_xxxzz_xyyyyy, g_z_0_xxxzz_xyyyyz, g_z_0_xxxzz_xyyyzz, g_z_0_xxxzz_xyyzzz, g_z_0_xxxzz_xyzzzz, g_z_0_xxxzz_xzzzzz, g_z_0_xxxzz_yyyyyy, g_z_0_xxxzz_yyyyyz, g_z_0_xxxzz_yyyyzz, g_z_0_xxxzz_yyyzzz, g_z_0_xxxzz_yyzzzz, g_z_0_xxxzz_yzzzzz, g_z_0_xxxzz_zzzzzz, g_z_0_xxzz_xxxxxx, g_z_0_xxzz_xxxxxxx, g_z_0_xxzz_xxxxxxy, g_z_0_xxzz_xxxxxxz, g_z_0_xxzz_xxxxxy, g_z_0_xxzz_xxxxxyy, g_z_0_xxzz_xxxxxyz, g_z_0_xxzz_xxxxxz, g_z_0_xxzz_xxxxxzz, g_z_0_xxzz_xxxxyy, g_z_0_xxzz_xxxxyyy, g_z_0_xxzz_xxxxyyz, g_z_0_xxzz_xxxxyz, g_z_0_xxzz_xxxxyzz, g_z_0_xxzz_xxxxzz, g_z_0_xxzz_xxxxzzz, g_z_0_xxzz_xxxyyy, g_z_0_xxzz_xxxyyyy, g_z_0_xxzz_xxxyyyz, g_z_0_xxzz_xxxyyz, g_z_0_xxzz_xxxyyzz, g_z_0_xxzz_xxxyzz, g_z_0_xxzz_xxxyzzz, g_z_0_xxzz_xxxzzz, g_z_0_xxzz_xxxzzzz, g_z_0_xxzz_xxyyyy, g_z_0_xxzz_xxyyyyy, g_z_0_xxzz_xxyyyyz, g_z_0_xxzz_xxyyyz, g_z_0_xxzz_xxyyyzz, g_z_0_xxzz_xxyyzz, g_z_0_xxzz_xxyyzzz, g_z_0_xxzz_xxyzzz, g_z_0_xxzz_xxyzzzz, g_z_0_xxzz_xxzzzz, g_z_0_xxzz_xxzzzzz, g_z_0_xxzz_xyyyyy, g_z_0_xxzz_xyyyyyy, g_z_0_xxzz_xyyyyyz, g_z_0_xxzz_xyyyyz, g_z_0_xxzz_xyyyyzz, g_z_0_xxzz_xyyyzz, g_z_0_xxzz_xyyyzzz, g_z_0_xxzz_xyyzzz, g_z_0_xxzz_xyyzzzz, g_z_0_xxzz_xyzzzz, g_z_0_xxzz_xyzzzzz, g_z_0_xxzz_xzzzzz, g_z_0_xxzz_xzzzzzz, g_z_0_xxzz_yyyyyy, g_z_0_xxzz_yyyyyz, g_z_0_xxzz_yyyyzz, g_z_0_xxzz_yyyzzz, g_z_0_xxzz_yyzzzz, g_z_0_xxzz_yzzzzz, g_z_0_xxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzz_xxxxxx[k] = -g_z_0_xxzz_xxxxxx[k] * ab_x + g_z_0_xxzz_xxxxxxx[k];

                g_z_0_xxxzz_xxxxxy[k] = -g_z_0_xxzz_xxxxxy[k] * ab_x + g_z_0_xxzz_xxxxxxy[k];

                g_z_0_xxxzz_xxxxxz[k] = -g_z_0_xxzz_xxxxxz[k] * ab_x + g_z_0_xxzz_xxxxxxz[k];

                g_z_0_xxxzz_xxxxyy[k] = -g_z_0_xxzz_xxxxyy[k] * ab_x + g_z_0_xxzz_xxxxxyy[k];

                g_z_0_xxxzz_xxxxyz[k] = -g_z_0_xxzz_xxxxyz[k] * ab_x + g_z_0_xxzz_xxxxxyz[k];

                g_z_0_xxxzz_xxxxzz[k] = -g_z_0_xxzz_xxxxzz[k] * ab_x + g_z_0_xxzz_xxxxxzz[k];

                g_z_0_xxxzz_xxxyyy[k] = -g_z_0_xxzz_xxxyyy[k] * ab_x + g_z_0_xxzz_xxxxyyy[k];

                g_z_0_xxxzz_xxxyyz[k] = -g_z_0_xxzz_xxxyyz[k] * ab_x + g_z_0_xxzz_xxxxyyz[k];

                g_z_0_xxxzz_xxxyzz[k] = -g_z_0_xxzz_xxxyzz[k] * ab_x + g_z_0_xxzz_xxxxyzz[k];

                g_z_0_xxxzz_xxxzzz[k] = -g_z_0_xxzz_xxxzzz[k] * ab_x + g_z_0_xxzz_xxxxzzz[k];

                g_z_0_xxxzz_xxyyyy[k] = -g_z_0_xxzz_xxyyyy[k] * ab_x + g_z_0_xxzz_xxxyyyy[k];

                g_z_0_xxxzz_xxyyyz[k] = -g_z_0_xxzz_xxyyyz[k] * ab_x + g_z_0_xxzz_xxxyyyz[k];

                g_z_0_xxxzz_xxyyzz[k] = -g_z_0_xxzz_xxyyzz[k] * ab_x + g_z_0_xxzz_xxxyyzz[k];

                g_z_0_xxxzz_xxyzzz[k] = -g_z_0_xxzz_xxyzzz[k] * ab_x + g_z_0_xxzz_xxxyzzz[k];

                g_z_0_xxxzz_xxzzzz[k] = -g_z_0_xxzz_xxzzzz[k] * ab_x + g_z_0_xxzz_xxxzzzz[k];

                g_z_0_xxxzz_xyyyyy[k] = -g_z_0_xxzz_xyyyyy[k] * ab_x + g_z_0_xxzz_xxyyyyy[k];

                g_z_0_xxxzz_xyyyyz[k] = -g_z_0_xxzz_xyyyyz[k] * ab_x + g_z_0_xxzz_xxyyyyz[k];

                g_z_0_xxxzz_xyyyzz[k] = -g_z_0_xxzz_xyyyzz[k] * ab_x + g_z_0_xxzz_xxyyyzz[k];

                g_z_0_xxxzz_xyyzzz[k] = -g_z_0_xxzz_xyyzzz[k] * ab_x + g_z_0_xxzz_xxyyzzz[k];

                g_z_0_xxxzz_xyzzzz[k] = -g_z_0_xxzz_xyzzzz[k] * ab_x + g_z_0_xxzz_xxyzzzz[k];

                g_z_0_xxxzz_xzzzzz[k] = -g_z_0_xxzz_xzzzzz[k] * ab_x + g_z_0_xxzz_xxzzzzz[k];

                g_z_0_xxxzz_yyyyyy[k] = -g_z_0_xxzz_yyyyyy[k] * ab_x + g_z_0_xxzz_xyyyyyy[k];

                g_z_0_xxxzz_yyyyyz[k] = -g_z_0_xxzz_yyyyyz[k] * ab_x + g_z_0_xxzz_xyyyyyz[k];

                g_z_0_xxxzz_yyyyzz[k] = -g_z_0_xxzz_yyyyzz[k] * ab_x + g_z_0_xxzz_xyyyyzz[k];

                g_z_0_xxxzz_yyyzzz[k] = -g_z_0_xxzz_yyyzzz[k] * ab_x + g_z_0_xxzz_xyyyzzz[k];

                g_z_0_xxxzz_yyzzzz[k] = -g_z_0_xxzz_yyzzzz[k] * ab_x + g_z_0_xxzz_xyyzzzz[k];

                g_z_0_xxxzz_yzzzzz[k] = -g_z_0_xxzz_yzzzzz[k] * ab_x + g_z_0_xxzz_xyzzzzz[k];

                g_z_0_xxxzz_zzzzzz[k] = -g_z_0_xxzz_zzzzzz[k] * ab_x + g_z_0_xxzz_xzzzzzz[k];
            }

            /// Set up 1344-1372 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1344 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1345 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1346 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1347 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1348 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1349 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1350 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1351 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1352 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1353 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1354 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1355 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1356 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1357 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1358 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1359 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1360 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1361 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1362 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1363 * ccomps * dcomps);

            auto g_z_0_xxyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1364 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1365 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1366 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1367 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1368 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1369 * ccomps * dcomps);

            auto g_z_0_xxyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1370 * ccomps * dcomps);

            auto g_z_0_xxyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1371 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyy_xxxxxx, g_z_0_xxyyy_xxxxxy, g_z_0_xxyyy_xxxxxz, g_z_0_xxyyy_xxxxyy, g_z_0_xxyyy_xxxxyz, g_z_0_xxyyy_xxxxzz, g_z_0_xxyyy_xxxyyy, g_z_0_xxyyy_xxxyyz, g_z_0_xxyyy_xxxyzz, g_z_0_xxyyy_xxxzzz, g_z_0_xxyyy_xxyyyy, g_z_0_xxyyy_xxyyyz, g_z_0_xxyyy_xxyyzz, g_z_0_xxyyy_xxyzzz, g_z_0_xxyyy_xxzzzz, g_z_0_xxyyy_xyyyyy, g_z_0_xxyyy_xyyyyz, g_z_0_xxyyy_xyyyzz, g_z_0_xxyyy_xyyzzz, g_z_0_xxyyy_xyzzzz, g_z_0_xxyyy_xzzzzz, g_z_0_xxyyy_yyyyyy, g_z_0_xxyyy_yyyyyz, g_z_0_xxyyy_yyyyzz, g_z_0_xxyyy_yyyzzz, g_z_0_xxyyy_yyzzzz, g_z_0_xxyyy_yzzzzz, g_z_0_xxyyy_zzzzzz, g_z_0_xyyy_xxxxxx, g_z_0_xyyy_xxxxxxx, g_z_0_xyyy_xxxxxxy, g_z_0_xyyy_xxxxxxz, g_z_0_xyyy_xxxxxy, g_z_0_xyyy_xxxxxyy, g_z_0_xyyy_xxxxxyz, g_z_0_xyyy_xxxxxz, g_z_0_xyyy_xxxxxzz, g_z_0_xyyy_xxxxyy, g_z_0_xyyy_xxxxyyy, g_z_0_xyyy_xxxxyyz, g_z_0_xyyy_xxxxyz, g_z_0_xyyy_xxxxyzz, g_z_0_xyyy_xxxxzz, g_z_0_xyyy_xxxxzzz, g_z_0_xyyy_xxxyyy, g_z_0_xyyy_xxxyyyy, g_z_0_xyyy_xxxyyyz, g_z_0_xyyy_xxxyyz, g_z_0_xyyy_xxxyyzz, g_z_0_xyyy_xxxyzz, g_z_0_xyyy_xxxyzzz, g_z_0_xyyy_xxxzzz, g_z_0_xyyy_xxxzzzz, g_z_0_xyyy_xxyyyy, g_z_0_xyyy_xxyyyyy, g_z_0_xyyy_xxyyyyz, g_z_0_xyyy_xxyyyz, g_z_0_xyyy_xxyyyzz, g_z_0_xyyy_xxyyzz, g_z_0_xyyy_xxyyzzz, g_z_0_xyyy_xxyzzz, g_z_0_xyyy_xxyzzzz, g_z_0_xyyy_xxzzzz, g_z_0_xyyy_xxzzzzz, g_z_0_xyyy_xyyyyy, g_z_0_xyyy_xyyyyyy, g_z_0_xyyy_xyyyyyz, g_z_0_xyyy_xyyyyz, g_z_0_xyyy_xyyyyzz, g_z_0_xyyy_xyyyzz, g_z_0_xyyy_xyyyzzz, g_z_0_xyyy_xyyzzz, g_z_0_xyyy_xyyzzzz, g_z_0_xyyy_xyzzzz, g_z_0_xyyy_xyzzzzz, g_z_0_xyyy_xzzzzz, g_z_0_xyyy_xzzzzzz, g_z_0_xyyy_yyyyyy, g_z_0_xyyy_yyyyyz, g_z_0_xyyy_yyyyzz, g_z_0_xyyy_yyyzzz, g_z_0_xyyy_yyzzzz, g_z_0_xyyy_yzzzzz, g_z_0_xyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyy_xxxxxx[k] = -g_z_0_xyyy_xxxxxx[k] * ab_x + g_z_0_xyyy_xxxxxxx[k];

                g_z_0_xxyyy_xxxxxy[k] = -g_z_0_xyyy_xxxxxy[k] * ab_x + g_z_0_xyyy_xxxxxxy[k];

                g_z_0_xxyyy_xxxxxz[k] = -g_z_0_xyyy_xxxxxz[k] * ab_x + g_z_0_xyyy_xxxxxxz[k];

                g_z_0_xxyyy_xxxxyy[k] = -g_z_0_xyyy_xxxxyy[k] * ab_x + g_z_0_xyyy_xxxxxyy[k];

                g_z_0_xxyyy_xxxxyz[k] = -g_z_0_xyyy_xxxxyz[k] * ab_x + g_z_0_xyyy_xxxxxyz[k];

                g_z_0_xxyyy_xxxxzz[k] = -g_z_0_xyyy_xxxxzz[k] * ab_x + g_z_0_xyyy_xxxxxzz[k];

                g_z_0_xxyyy_xxxyyy[k] = -g_z_0_xyyy_xxxyyy[k] * ab_x + g_z_0_xyyy_xxxxyyy[k];

                g_z_0_xxyyy_xxxyyz[k] = -g_z_0_xyyy_xxxyyz[k] * ab_x + g_z_0_xyyy_xxxxyyz[k];

                g_z_0_xxyyy_xxxyzz[k] = -g_z_0_xyyy_xxxyzz[k] * ab_x + g_z_0_xyyy_xxxxyzz[k];

                g_z_0_xxyyy_xxxzzz[k] = -g_z_0_xyyy_xxxzzz[k] * ab_x + g_z_0_xyyy_xxxxzzz[k];

                g_z_0_xxyyy_xxyyyy[k] = -g_z_0_xyyy_xxyyyy[k] * ab_x + g_z_0_xyyy_xxxyyyy[k];

                g_z_0_xxyyy_xxyyyz[k] = -g_z_0_xyyy_xxyyyz[k] * ab_x + g_z_0_xyyy_xxxyyyz[k];

                g_z_0_xxyyy_xxyyzz[k] = -g_z_0_xyyy_xxyyzz[k] * ab_x + g_z_0_xyyy_xxxyyzz[k];

                g_z_0_xxyyy_xxyzzz[k] = -g_z_0_xyyy_xxyzzz[k] * ab_x + g_z_0_xyyy_xxxyzzz[k];

                g_z_0_xxyyy_xxzzzz[k] = -g_z_0_xyyy_xxzzzz[k] * ab_x + g_z_0_xyyy_xxxzzzz[k];

                g_z_0_xxyyy_xyyyyy[k] = -g_z_0_xyyy_xyyyyy[k] * ab_x + g_z_0_xyyy_xxyyyyy[k];

                g_z_0_xxyyy_xyyyyz[k] = -g_z_0_xyyy_xyyyyz[k] * ab_x + g_z_0_xyyy_xxyyyyz[k];

                g_z_0_xxyyy_xyyyzz[k] = -g_z_0_xyyy_xyyyzz[k] * ab_x + g_z_0_xyyy_xxyyyzz[k];

                g_z_0_xxyyy_xyyzzz[k] = -g_z_0_xyyy_xyyzzz[k] * ab_x + g_z_0_xyyy_xxyyzzz[k];

                g_z_0_xxyyy_xyzzzz[k] = -g_z_0_xyyy_xyzzzz[k] * ab_x + g_z_0_xyyy_xxyzzzz[k];

                g_z_0_xxyyy_xzzzzz[k] = -g_z_0_xyyy_xzzzzz[k] * ab_x + g_z_0_xyyy_xxzzzzz[k];

                g_z_0_xxyyy_yyyyyy[k] = -g_z_0_xyyy_yyyyyy[k] * ab_x + g_z_0_xyyy_xyyyyyy[k];

                g_z_0_xxyyy_yyyyyz[k] = -g_z_0_xyyy_yyyyyz[k] * ab_x + g_z_0_xyyy_xyyyyyz[k];

                g_z_0_xxyyy_yyyyzz[k] = -g_z_0_xyyy_yyyyzz[k] * ab_x + g_z_0_xyyy_xyyyyzz[k];

                g_z_0_xxyyy_yyyzzz[k] = -g_z_0_xyyy_yyyzzz[k] * ab_x + g_z_0_xyyy_xyyyzzz[k];

                g_z_0_xxyyy_yyzzzz[k] = -g_z_0_xyyy_yyzzzz[k] * ab_x + g_z_0_xyyy_xyyzzzz[k];

                g_z_0_xxyyy_yzzzzz[k] = -g_z_0_xyyy_yzzzzz[k] * ab_x + g_z_0_xyyy_xyzzzzz[k];

                g_z_0_xxyyy_zzzzzz[k] = -g_z_0_xyyy_zzzzzz[k] * ab_x + g_z_0_xyyy_xzzzzzz[k];
            }

            /// Set up 1372-1400 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1372 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1373 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1374 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1375 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1376 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1377 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1378 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1379 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1380 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1381 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1382 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1383 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1384 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1385 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1386 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1387 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1388 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1389 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1390 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1391 * ccomps * dcomps);

            auto g_z_0_xxyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1392 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1393 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1394 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1395 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1396 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1397 * ccomps * dcomps);

            auto g_z_0_xxyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1398 * ccomps * dcomps);

            auto g_z_0_xxyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1399 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyz_xxxxxx, g_z_0_xxyyz_xxxxxy, g_z_0_xxyyz_xxxxxz, g_z_0_xxyyz_xxxxyy, g_z_0_xxyyz_xxxxyz, g_z_0_xxyyz_xxxxzz, g_z_0_xxyyz_xxxyyy, g_z_0_xxyyz_xxxyyz, g_z_0_xxyyz_xxxyzz, g_z_0_xxyyz_xxxzzz, g_z_0_xxyyz_xxyyyy, g_z_0_xxyyz_xxyyyz, g_z_0_xxyyz_xxyyzz, g_z_0_xxyyz_xxyzzz, g_z_0_xxyyz_xxzzzz, g_z_0_xxyyz_xyyyyy, g_z_0_xxyyz_xyyyyz, g_z_0_xxyyz_xyyyzz, g_z_0_xxyyz_xyyzzz, g_z_0_xxyyz_xyzzzz, g_z_0_xxyyz_xzzzzz, g_z_0_xxyyz_yyyyyy, g_z_0_xxyyz_yyyyyz, g_z_0_xxyyz_yyyyzz, g_z_0_xxyyz_yyyzzz, g_z_0_xxyyz_yyzzzz, g_z_0_xxyyz_yzzzzz, g_z_0_xxyyz_zzzzzz, g_z_0_xyyz_xxxxxx, g_z_0_xyyz_xxxxxxx, g_z_0_xyyz_xxxxxxy, g_z_0_xyyz_xxxxxxz, g_z_0_xyyz_xxxxxy, g_z_0_xyyz_xxxxxyy, g_z_0_xyyz_xxxxxyz, g_z_0_xyyz_xxxxxz, g_z_0_xyyz_xxxxxzz, g_z_0_xyyz_xxxxyy, g_z_0_xyyz_xxxxyyy, g_z_0_xyyz_xxxxyyz, g_z_0_xyyz_xxxxyz, g_z_0_xyyz_xxxxyzz, g_z_0_xyyz_xxxxzz, g_z_0_xyyz_xxxxzzz, g_z_0_xyyz_xxxyyy, g_z_0_xyyz_xxxyyyy, g_z_0_xyyz_xxxyyyz, g_z_0_xyyz_xxxyyz, g_z_0_xyyz_xxxyyzz, g_z_0_xyyz_xxxyzz, g_z_0_xyyz_xxxyzzz, g_z_0_xyyz_xxxzzz, g_z_0_xyyz_xxxzzzz, g_z_0_xyyz_xxyyyy, g_z_0_xyyz_xxyyyyy, g_z_0_xyyz_xxyyyyz, g_z_0_xyyz_xxyyyz, g_z_0_xyyz_xxyyyzz, g_z_0_xyyz_xxyyzz, g_z_0_xyyz_xxyyzzz, g_z_0_xyyz_xxyzzz, g_z_0_xyyz_xxyzzzz, g_z_0_xyyz_xxzzzz, g_z_0_xyyz_xxzzzzz, g_z_0_xyyz_xyyyyy, g_z_0_xyyz_xyyyyyy, g_z_0_xyyz_xyyyyyz, g_z_0_xyyz_xyyyyz, g_z_0_xyyz_xyyyyzz, g_z_0_xyyz_xyyyzz, g_z_0_xyyz_xyyyzzz, g_z_0_xyyz_xyyzzz, g_z_0_xyyz_xyyzzzz, g_z_0_xyyz_xyzzzz, g_z_0_xyyz_xyzzzzz, g_z_0_xyyz_xzzzzz, g_z_0_xyyz_xzzzzzz, g_z_0_xyyz_yyyyyy, g_z_0_xyyz_yyyyyz, g_z_0_xyyz_yyyyzz, g_z_0_xyyz_yyyzzz, g_z_0_xyyz_yyzzzz, g_z_0_xyyz_yzzzzz, g_z_0_xyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyz_xxxxxx[k] = -g_z_0_xyyz_xxxxxx[k] * ab_x + g_z_0_xyyz_xxxxxxx[k];

                g_z_0_xxyyz_xxxxxy[k] = -g_z_0_xyyz_xxxxxy[k] * ab_x + g_z_0_xyyz_xxxxxxy[k];

                g_z_0_xxyyz_xxxxxz[k] = -g_z_0_xyyz_xxxxxz[k] * ab_x + g_z_0_xyyz_xxxxxxz[k];

                g_z_0_xxyyz_xxxxyy[k] = -g_z_0_xyyz_xxxxyy[k] * ab_x + g_z_0_xyyz_xxxxxyy[k];

                g_z_0_xxyyz_xxxxyz[k] = -g_z_0_xyyz_xxxxyz[k] * ab_x + g_z_0_xyyz_xxxxxyz[k];

                g_z_0_xxyyz_xxxxzz[k] = -g_z_0_xyyz_xxxxzz[k] * ab_x + g_z_0_xyyz_xxxxxzz[k];

                g_z_0_xxyyz_xxxyyy[k] = -g_z_0_xyyz_xxxyyy[k] * ab_x + g_z_0_xyyz_xxxxyyy[k];

                g_z_0_xxyyz_xxxyyz[k] = -g_z_0_xyyz_xxxyyz[k] * ab_x + g_z_0_xyyz_xxxxyyz[k];

                g_z_0_xxyyz_xxxyzz[k] = -g_z_0_xyyz_xxxyzz[k] * ab_x + g_z_0_xyyz_xxxxyzz[k];

                g_z_0_xxyyz_xxxzzz[k] = -g_z_0_xyyz_xxxzzz[k] * ab_x + g_z_0_xyyz_xxxxzzz[k];

                g_z_0_xxyyz_xxyyyy[k] = -g_z_0_xyyz_xxyyyy[k] * ab_x + g_z_0_xyyz_xxxyyyy[k];

                g_z_0_xxyyz_xxyyyz[k] = -g_z_0_xyyz_xxyyyz[k] * ab_x + g_z_0_xyyz_xxxyyyz[k];

                g_z_0_xxyyz_xxyyzz[k] = -g_z_0_xyyz_xxyyzz[k] * ab_x + g_z_0_xyyz_xxxyyzz[k];

                g_z_0_xxyyz_xxyzzz[k] = -g_z_0_xyyz_xxyzzz[k] * ab_x + g_z_0_xyyz_xxxyzzz[k];

                g_z_0_xxyyz_xxzzzz[k] = -g_z_0_xyyz_xxzzzz[k] * ab_x + g_z_0_xyyz_xxxzzzz[k];

                g_z_0_xxyyz_xyyyyy[k] = -g_z_0_xyyz_xyyyyy[k] * ab_x + g_z_0_xyyz_xxyyyyy[k];

                g_z_0_xxyyz_xyyyyz[k] = -g_z_0_xyyz_xyyyyz[k] * ab_x + g_z_0_xyyz_xxyyyyz[k];

                g_z_0_xxyyz_xyyyzz[k] = -g_z_0_xyyz_xyyyzz[k] * ab_x + g_z_0_xyyz_xxyyyzz[k];

                g_z_0_xxyyz_xyyzzz[k] = -g_z_0_xyyz_xyyzzz[k] * ab_x + g_z_0_xyyz_xxyyzzz[k];

                g_z_0_xxyyz_xyzzzz[k] = -g_z_0_xyyz_xyzzzz[k] * ab_x + g_z_0_xyyz_xxyzzzz[k];

                g_z_0_xxyyz_xzzzzz[k] = -g_z_0_xyyz_xzzzzz[k] * ab_x + g_z_0_xyyz_xxzzzzz[k];

                g_z_0_xxyyz_yyyyyy[k] = -g_z_0_xyyz_yyyyyy[k] * ab_x + g_z_0_xyyz_xyyyyyy[k];

                g_z_0_xxyyz_yyyyyz[k] = -g_z_0_xyyz_yyyyyz[k] * ab_x + g_z_0_xyyz_xyyyyyz[k];

                g_z_0_xxyyz_yyyyzz[k] = -g_z_0_xyyz_yyyyzz[k] * ab_x + g_z_0_xyyz_xyyyyzz[k];

                g_z_0_xxyyz_yyyzzz[k] = -g_z_0_xyyz_yyyzzz[k] * ab_x + g_z_0_xyyz_xyyyzzz[k];

                g_z_0_xxyyz_yyzzzz[k] = -g_z_0_xyyz_yyzzzz[k] * ab_x + g_z_0_xyyz_xyyzzzz[k];

                g_z_0_xxyyz_yzzzzz[k] = -g_z_0_xyyz_yzzzzz[k] * ab_x + g_z_0_xyyz_xyzzzzz[k];

                g_z_0_xxyyz_zzzzzz[k] = -g_z_0_xyyz_zzzzzz[k] * ab_x + g_z_0_xyyz_xzzzzzz[k];
            }

            /// Set up 1400-1428 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1400 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1401 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1402 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1403 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1404 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1405 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1406 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1407 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1408 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1409 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1410 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1411 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1412 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1413 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1414 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1415 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1416 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1417 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1418 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1419 * ccomps * dcomps);

            auto g_z_0_xxyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1420 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1421 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1422 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1423 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1424 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1425 * ccomps * dcomps);

            auto g_z_0_xxyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1426 * ccomps * dcomps);

            auto g_z_0_xxyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1427 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyzz_xxxxxx, g_z_0_xxyzz_xxxxxy, g_z_0_xxyzz_xxxxxz, g_z_0_xxyzz_xxxxyy, g_z_0_xxyzz_xxxxyz, g_z_0_xxyzz_xxxxzz, g_z_0_xxyzz_xxxyyy, g_z_0_xxyzz_xxxyyz, g_z_0_xxyzz_xxxyzz, g_z_0_xxyzz_xxxzzz, g_z_0_xxyzz_xxyyyy, g_z_0_xxyzz_xxyyyz, g_z_0_xxyzz_xxyyzz, g_z_0_xxyzz_xxyzzz, g_z_0_xxyzz_xxzzzz, g_z_0_xxyzz_xyyyyy, g_z_0_xxyzz_xyyyyz, g_z_0_xxyzz_xyyyzz, g_z_0_xxyzz_xyyzzz, g_z_0_xxyzz_xyzzzz, g_z_0_xxyzz_xzzzzz, g_z_0_xxyzz_yyyyyy, g_z_0_xxyzz_yyyyyz, g_z_0_xxyzz_yyyyzz, g_z_0_xxyzz_yyyzzz, g_z_0_xxyzz_yyzzzz, g_z_0_xxyzz_yzzzzz, g_z_0_xxyzz_zzzzzz, g_z_0_xyzz_xxxxxx, g_z_0_xyzz_xxxxxxx, g_z_0_xyzz_xxxxxxy, g_z_0_xyzz_xxxxxxz, g_z_0_xyzz_xxxxxy, g_z_0_xyzz_xxxxxyy, g_z_0_xyzz_xxxxxyz, g_z_0_xyzz_xxxxxz, g_z_0_xyzz_xxxxxzz, g_z_0_xyzz_xxxxyy, g_z_0_xyzz_xxxxyyy, g_z_0_xyzz_xxxxyyz, g_z_0_xyzz_xxxxyz, g_z_0_xyzz_xxxxyzz, g_z_0_xyzz_xxxxzz, g_z_0_xyzz_xxxxzzz, g_z_0_xyzz_xxxyyy, g_z_0_xyzz_xxxyyyy, g_z_0_xyzz_xxxyyyz, g_z_0_xyzz_xxxyyz, g_z_0_xyzz_xxxyyzz, g_z_0_xyzz_xxxyzz, g_z_0_xyzz_xxxyzzz, g_z_0_xyzz_xxxzzz, g_z_0_xyzz_xxxzzzz, g_z_0_xyzz_xxyyyy, g_z_0_xyzz_xxyyyyy, g_z_0_xyzz_xxyyyyz, g_z_0_xyzz_xxyyyz, g_z_0_xyzz_xxyyyzz, g_z_0_xyzz_xxyyzz, g_z_0_xyzz_xxyyzzz, g_z_0_xyzz_xxyzzz, g_z_0_xyzz_xxyzzzz, g_z_0_xyzz_xxzzzz, g_z_0_xyzz_xxzzzzz, g_z_0_xyzz_xyyyyy, g_z_0_xyzz_xyyyyyy, g_z_0_xyzz_xyyyyyz, g_z_0_xyzz_xyyyyz, g_z_0_xyzz_xyyyyzz, g_z_0_xyzz_xyyyzz, g_z_0_xyzz_xyyyzzz, g_z_0_xyzz_xyyzzz, g_z_0_xyzz_xyyzzzz, g_z_0_xyzz_xyzzzz, g_z_0_xyzz_xyzzzzz, g_z_0_xyzz_xzzzzz, g_z_0_xyzz_xzzzzzz, g_z_0_xyzz_yyyyyy, g_z_0_xyzz_yyyyyz, g_z_0_xyzz_yyyyzz, g_z_0_xyzz_yyyzzz, g_z_0_xyzz_yyzzzz, g_z_0_xyzz_yzzzzz, g_z_0_xyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzz_xxxxxx[k] = -g_z_0_xyzz_xxxxxx[k] * ab_x + g_z_0_xyzz_xxxxxxx[k];

                g_z_0_xxyzz_xxxxxy[k] = -g_z_0_xyzz_xxxxxy[k] * ab_x + g_z_0_xyzz_xxxxxxy[k];

                g_z_0_xxyzz_xxxxxz[k] = -g_z_0_xyzz_xxxxxz[k] * ab_x + g_z_0_xyzz_xxxxxxz[k];

                g_z_0_xxyzz_xxxxyy[k] = -g_z_0_xyzz_xxxxyy[k] * ab_x + g_z_0_xyzz_xxxxxyy[k];

                g_z_0_xxyzz_xxxxyz[k] = -g_z_0_xyzz_xxxxyz[k] * ab_x + g_z_0_xyzz_xxxxxyz[k];

                g_z_0_xxyzz_xxxxzz[k] = -g_z_0_xyzz_xxxxzz[k] * ab_x + g_z_0_xyzz_xxxxxzz[k];

                g_z_0_xxyzz_xxxyyy[k] = -g_z_0_xyzz_xxxyyy[k] * ab_x + g_z_0_xyzz_xxxxyyy[k];

                g_z_0_xxyzz_xxxyyz[k] = -g_z_0_xyzz_xxxyyz[k] * ab_x + g_z_0_xyzz_xxxxyyz[k];

                g_z_0_xxyzz_xxxyzz[k] = -g_z_0_xyzz_xxxyzz[k] * ab_x + g_z_0_xyzz_xxxxyzz[k];

                g_z_0_xxyzz_xxxzzz[k] = -g_z_0_xyzz_xxxzzz[k] * ab_x + g_z_0_xyzz_xxxxzzz[k];

                g_z_0_xxyzz_xxyyyy[k] = -g_z_0_xyzz_xxyyyy[k] * ab_x + g_z_0_xyzz_xxxyyyy[k];

                g_z_0_xxyzz_xxyyyz[k] = -g_z_0_xyzz_xxyyyz[k] * ab_x + g_z_0_xyzz_xxxyyyz[k];

                g_z_0_xxyzz_xxyyzz[k] = -g_z_0_xyzz_xxyyzz[k] * ab_x + g_z_0_xyzz_xxxyyzz[k];

                g_z_0_xxyzz_xxyzzz[k] = -g_z_0_xyzz_xxyzzz[k] * ab_x + g_z_0_xyzz_xxxyzzz[k];

                g_z_0_xxyzz_xxzzzz[k] = -g_z_0_xyzz_xxzzzz[k] * ab_x + g_z_0_xyzz_xxxzzzz[k];

                g_z_0_xxyzz_xyyyyy[k] = -g_z_0_xyzz_xyyyyy[k] * ab_x + g_z_0_xyzz_xxyyyyy[k];

                g_z_0_xxyzz_xyyyyz[k] = -g_z_0_xyzz_xyyyyz[k] * ab_x + g_z_0_xyzz_xxyyyyz[k];

                g_z_0_xxyzz_xyyyzz[k] = -g_z_0_xyzz_xyyyzz[k] * ab_x + g_z_0_xyzz_xxyyyzz[k];

                g_z_0_xxyzz_xyyzzz[k] = -g_z_0_xyzz_xyyzzz[k] * ab_x + g_z_0_xyzz_xxyyzzz[k];

                g_z_0_xxyzz_xyzzzz[k] = -g_z_0_xyzz_xyzzzz[k] * ab_x + g_z_0_xyzz_xxyzzzz[k];

                g_z_0_xxyzz_xzzzzz[k] = -g_z_0_xyzz_xzzzzz[k] * ab_x + g_z_0_xyzz_xxzzzzz[k];

                g_z_0_xxyzz_yyyyyy[k] = -g_z_0_xyzz_yyyyyy[k] * ab_x + g_z_0_xyzz_xyyyyyy[k];

                g_z_0_xxyzz_yyyyyz[k] = -g_z_0_xyzz_yyyyyz[k] * ab_x + g_z_0_xyzz_xyyyyyz[k];

                g_z_0_xxyzz_yyyyzz[k] = -g_z_0_xyzz_yyyyzz[k] * ab_x + g_z_0_xyzz_xyyyyzz[k];

                g_z_0_xxyzz_yyyzzz[k] = -g_z_0_xyzz_yyyzzz[k] * ab_x + g_z_0_xyzz_xyyyzzz[k];

                g_z_0_xxyzz_yyzzzz[k] = -g_z_0_xyzz_yyzzzz[k] * ab_x + g_z_0_xyzz_xyyzzzz[k];

                g_z_0_xxyzz_yzzzzz[k] = -g_z_0_xyzz_yzzzzz[k] * ab_x + g_z_0_xyzz_xyzzzzz[k];

                g_z_0_xxyzz_zzzzzz[k] = -g_z_0_xyzz_zzzzzz[k] * ab_x + g_z_0_xyzz_xzzzzzz[k];
            }

            /// Set up 1428-1456 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1428 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1429 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1430 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1431 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1432 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1433 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1434 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1435 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1436 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1437 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1438 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1439 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1440 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1441 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1442 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1443 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1444 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1445 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1446 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1447 * ccomps * dcomps);

            auto g_z_0_xxzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1448 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1449 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1450 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1451 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1452 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1453 * ccomps * dcomps);

            auto g_z_0_xxzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1454 * ccomps * dcomps);

            auto g_z_0_xxzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1455 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzzz_xxxxxx, g_z_0_xxzzz_xxxxxy, g_z_0_xxzzz_xxxxxz, g_z_0_xxzzz_xxxxyy, g_z_0_xxzzz_xxxxyz, g_z_0_xxzzz_xxxxzz, g_z_0_xxzzz_xxxyyy, g_z_0_xxzzz_xxxyyz, g_z_0_xxzzz_xxxyzz, g_z_0_xxzzz_xxxzzz, g_z_0_xxzzz_xxyyyy, g_z_0_xxzzz_xxyyyz, g_z_0_xxzzz_xxyyzz, g_z_0_xxzzz_xxyzzz, g_z_0_xxzzz_xxzzzz, g_z_0_xxzzz_xyyyyy, g_z_0_xxzzz_xyyyyz, g_z_0_xxzzz_xyyyzz, g_z_0_xxzzz_xyyzzz, g_z_0_xxzzz_xyzzzz, g_z_0_xxzzz_xzzzzz, g_z_0_xxzzz_yyyyyy, g_z_0_xxzzz_yyyyyz, g_z_0_xxzzz_yyyyzz, g_z_0_xxzzz_yyyzzz, g_z_0_xxzzz_yyzzzz, g_z_0_xxzzz_yzzzzz, g_z_0_xxzzz_zzzzzz, g_z_0_xzzz_xxxxxx, g_z_0_xzzz_xxxxxxx, g_z_0_xzzz_xxxxxxy, g_z_0_xzzz_xxxxxxz, g_z_0_xzzz_xxxxxy, g_z_0_xzzz_xxxxxyy, g_z_0_xzzz_xxxxxyz, g_z_0_xzzz_xxxxxz, g_z_0_xzzz_xxxxxzz, g_z_0_xzzz_xxxxyy, g_z_0_xzzz_xxxxyyy, g_z_0_xzzz_xxxxyyz, g_z_0_xzzz_xxxxyz, g_z_0_xzzz_xxxxyzz, g_z_0_xzzz_xxxxzz, g_z_0_xzzz_xxxxzzz, g_z_0_xzzz_xxxyyy, g_z_0_xzzz_xxxyyyy, g_z_0_xzzz_xxxyyyz, g_z_0_xzzz_xxxyyz, g_z_0_xzzz_xxxyyzz, g_z_0_xzzz_xxxyzz, g_z_0_xzzz_xxxyzzz, g_z_0_xzzz_xxxzzz, g_z_0_xzzz_xxxzzzz, g_z_0_xzzz_xxyyyy, g_z_0_xzzz_xxyyyyy, g_z_0_xzzz_xxyyyyz, g_z_0_xzzz_xxyyyz, g_z_0_xzzz_xxyyyzz, g_z_0_xzzz_xxyyzz, g_z_0_xzzz_xxyyzzz, g_z_0_xzzz_xxyzzz, g_z_0_xzzz_xxyzzzz, g_z_0_xzzz_xxzzzz, g_z_0_xzzz_xxzzzzz, g_z_0_xzzz_xyyyyy, g_z_0_xzzz_xyyyyyy, g_z_0_xzzz_xyyyyyz, g_z_0_xzzz_xyyyyz, g_z_0_xzzz_xyyyyzz, g_z_0_xzzz_xyyyzz, g_z_0_xzzz_xyyyzzz, g_z_0_xzzz_xyyzzz, g_z_0_xzzz_xyyzzzz, g_z_0_xzzz_xyzzzz, g_z_0_xzzz_xyzzzzz, g_z_0_xzzz_xzzzzz, g_z_0_xzzz_xzzzzzz, g_z_0_xzzz_yyyyyy, g_z_0_xzzz_yyyyyz, g_z_0_xzzz_yyyyzz, g_z_0_xzzz_yyyzzz, g_z_0_xzzz_yyzzzz, g_z_0_xzzz_yzzzzz, g_z_0_xzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzz_xxxxxx[k] = -g_z_0_xzzz_xxxxxx[k] * ab_x + g_z_0_xzzz_xxxxxxx[k];

                g_z_0_xxzzz_xxxxxy[k] = -g_z_0_xzzz_xxxxxy[k] * ab_x + g_z_0_xzzz_xxxxxxy[k];

                g_z_0_xxzzz_xxxxxz[k] = -g_z_0_xzzz_xxxxxz[k] * ab_x + g_z_0_xzzz_xxxxxxz[k];

                g_z_0_xxzzz_xxxxyy[k] = -g_z_0_xzzz_xxxxyy[k] * ab_x + g_z_0_xzzz_xxxxxyy[k];

                g_z_0_xxzzz_xxxxyz[k] = -g_z_0_xzzz_xxxxyz[k] * ab_x + g_z_0_xzzz_xxxxxyz[k];

                g_z_0_xxzzz_xxxxzz[k] = -g_z_0_xzzz_xxxxzz[k] * ab_x + g_z_0_xzzz_xxxxxzz[k];

                g_z_0_xxzzz_xxxyyy[k] = -g_z_0_xzzz_xxxyyy[k] * ab_x + g_z_0_xzzz_xxxxyyy[k];

                g_z_0_xxzzz_xxxyyz[k] = -g_z_0_xzzz_xxxyyz[k] * ab_x + g_z_0_xzzz_xxxxyyz[k];

                g_z_0_xxzzz_xxxyzz[k] = -g_z_0_xzzz_xxxyzz[k] * ab_x + g_z_0_xzzz_xxxxyzz[k];

                g_z_0_xxzzz_xxxzzz[k] = -g_z_0_xzzz_xxxzzz[k] * ab_x + g_z_0_xzzz_xxxxzzz[k];

                g_z_0_xxzzz_xxyyyy[k] = -g_z_0_xzzz_xxyyyy[k] * ab_x + g_z_0_xzzz_xxxyyyy[k];

                g_z_0_xxzzz_xxyyyz[k] = -g_z_0_xzzz_xxyyyz[k] * ab_x + g_z_0_xzzz_xxxyyyz[k];

                g_z_0_xxzzz_xxyyzz[k] = -g_z_0_xzzz_xxyyzz[k] * ab_x + g_z_0_xzzz_xxxyyzz[k];

                g_z_0_xxzzz_xxyzzz[k] = -g_z_0_xzzz_xxyzzz[k] * ab_x + g_z_0_xzzz_xxxyzzz[k];

                g_z_0_xxzzz_xxzzzz[k] = -g_z_0_xzzz_xxzzzz[k] * ab_x + g_z_0_xzzz_xxxzzzz[k];

                g_z_0_xxzzz_xyyyyy[k] = -g_z_0_xzzz_xyyyyy[k] * ab_x + g_z_0_xzzz_xxyyyyy[k];

                g_z_0_xxzzz_xyyyyz[k] = -g_z_0_xzzz_xyyyyz[k] * ab_x + g_z_0_xzzz_xxyyyyz[k];

                g_z_0_xxzzz_xyyyzz[k] = -g_z_0_xzzz_xyyyzz[k] * ab_x + g_z_0_xzzz_xxyyyzz[k];

                g_z_0_xxzzz_xyyzzz[k] = -g_z_0_xzzz_xyyzzz[k] * ab_x + g_z_0_xzzz_xxyyzzz[k];

                g_z_0_xxzzz_xyzzzz[k] = -g_z_0_xzzz_xyzzzz[k] * ab_x + g_z_0_xzzz_xxyzzzz[k];

                g_z_0_xxzzz_xzzzzz[k] = -g_z_0_xzzz_xzzzzz[k] * ab_x + g_z_0_xzzz_xxzzzzz[k];

                g_z_0_xxzzz_yyyyyy[k] = -g_z_0_xzzz_yyyyyy[k] * ab_x + g_z_0_xzzz_xyyyyyy[k];

                g_z_0_xxzzz_yyyyyz[k] = -g_z_0_xzzz_yyyyyz[k] * ab_x + g_z_0_xzzz_xyyyyyz[k];

                g_z_0_xxzzz_yyyyzz[k] = -g_z_0_xzzz_yyyyzz[k] * ab_x + g_z_0_xzzz_xyyyyzz[k];

                g_z_0_xxzzz_yyyzzz[k] = -g_z_0_xzzz_yyyzzz[k] * ab_x + g_z_0_xzzz_xyyyzzz[k];

                g_z_0_xxzzz_yyzzzz[k] = -g_z_0_xzzz_yyzzzz[k] * ab_x + g_z_0_xzzz_xyyzzzz[k];

                g_z_0_xxzzz_yzzzzz[k] = -g_z_0_xzzz_yzzzzz[k] * ab_x + g_z_0_xzzz_xyzzzzz[k];

                g_z_0_xxzzz_zzzzzz[k] = -g_z_0_xzzz_zzzzzz[k] * ab_x + g_z_0_xzzz_xzzzzzz[k];
            }

            /// Set up 1456-1484 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1456 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1457 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1458 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1459 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1460 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1461 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1462 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1463 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1464 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1465 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1466 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1467 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1468 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1469 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1470 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1471 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1472 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1473 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1474 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1475 * ccomps * dcomps);

            auto g_z_0_xyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1476 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1477 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1478 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1479 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1480 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1481 * ccomps * dcomps);

            auto g_z_0_xyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1482 * ccomps * dcomps);

            auto g_z_0_xyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1483 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyy_xxxxxx, g_z_0_xyyyy_xxxxxy, g_z_0_xyyyy_xxxxxz, g_z_0_xyyyy_xxxxyy, g_z_0_xyyyy_xxxxyz, g_z_0_xyyyy_xxxxzz, g_z_0_xyyyy_xxxyyy, g_z_0_xyyyy_xxxyyz, g_z_0_xyyyy_xxxyzz, g_z_0_xyyyy_xxxzzz, g_z_0_xyyyy_xxyyyy, g_z_0_xyyyy_xxyyyz, g_z_0_xyyyy_xxyyzz, g_z_0_xyyyy_xxyzzz, g_z_0_xyyyy_xxzzzz, g_z_0_xyyyy_xyyyyy, g_z_0_xyyyy_xyyyyz, g_z_0_xyyyy_xyyyzz, g_z_0_xyyyy_xyyzzz, g_z_0_xyyyy_xyzzzz, g_z_0_xyyyy_xzzzzz, g_z_0_xyyyy_yyyyyy, g_z_0_xyyyy_yyyyyz, g_z_0_xyyyy_yyyyzz, g_z_0_xyyyy_yyyzzz, g_z_0_xyyyy_yyzzzz, g_z_0_xyyyy_yzzzzz, g_z_0_xyyyy_zzzzzz, g_z_0_yyyy_xxxxxx, g_z_0_yyyy_xxxxxxx, g_z_0_yyyy_xxxxxxy, g_z_0_yyyy_xxxxxxz, g_z_0_yyyy_xxxxxy, g_z_0_yyyy_xxxxxyy, g_z_0_yyyy_xxxxxyz, g_z_0_yyyy_xxxxxz, g_z_0_yyyy_xxxxxzz, g_z_0_yyyy_xxxxyy, g_z_0_yyyy_xxxxyyy, g_z_0_yyyy_xxxxyyz, g_z_0_yyyy_xxxxyz, g_z_0_yyyy_xxxxyzz, g_z_0_yyyy_xxxxzz, g_z_0_yyyy_xxxxzzz, g_z_0_yyyy_xxxyyy, g_z_0_yyyy_xxxyyyy, g_z_0_yyyy_xxxyyyz, g_z_0_yyyy_xxxyyz, g_z_0_yyyy_xxxyyzz, g_z_0_yyyy_xxxyzz, g_z_0_yyyy_xxxyzzz, g_z_0_yyyy_xxxzzz, g_z_0_yyyy_xxxzzzz, g_z_0_yyyy_xxyyyy, g_z_0_yyyy_xxyyyyy, g_z_0_yyyy_xxyyyyz, g_z_0_yyyy_xxyyyz, g_z_0_yyyy_xxyyyzz, g_z_0_yyyy_xxyyzz, g_z_0_yyyy_xxyyzzz, g_z_0_yyyy_xxyzzz, g_z_0_yyyy_xxyzzzz, g_z_0_yyyy_xxzzzz, g_z_0_yyyy_xxzzzzz, g_z_0_yyyy_xyyyyy, g_z_0_yyyy_xyyyyyy, g_z_0_yyyy_xyyyyyz, g_z_0_yyyy_xyyyyz, g_z_0_yyyy_xyyyyzz, g_z_0_yyyy_xyyyzz, g_z_0_yyyy_xyyyzzz, g_z_0_yyyy_xyyzzz, g_z_0_yyyy_xyyzzzz, g_z_0_yyyy_xyzzzz, g_z_0_yyyy_xyzzzzz, g_z_0_yyyy_xzzzzz, g_z_0_yyyy_xzzzzzz, g_z_0_yyyy_yyyyyy, g_z_0_yyyy_yyyyyz, g_z_0_yyyy_yyyyzz, g_z_0_yyyy_yyyzzz, g_z_0_yyyy_yyzzzz, g_z_0_yyyy_yzzzzz, g_z_0_yyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyy_xxxxxx[k] = -g_z_0_yyyy_xxxxxx[k] * ab_x + g_z_0_yyyy_xxxxxxx[k];

                g_z_0_xyyyy_xxxxxy[k] = -g_z_0_yyyy_xxxxxy[k] * ab_x + g_z_0_yyyy_xxxxxxy[k];

                g_z_0_xyyyy_xxxxxz[k] = -g_z_0_yyyy_xxxxxz[k] * ab_x + g_z_0_yyyy_xxxxxxz[k];

                g_z_0_xyyyy_xxxxyy[k] = -g_z_0_yyyy_xxxxyy[k] * ab_x + g_z_0_yyyy_xxxxxyy[k];

                g_z_0_xyyyy_xxxxyz[k] = -g_z_0_yyyy_xxxxyz[k] * ab_x + g_z_0_yyyy_xxxxxyz[k];

                g_z_0_xyyyy_xxxxzz[k] = -g_z_0_yyyy_xxxxzz[k] * ab_x + g_z_0_yyyy_xxxxxzz[k];

                g_z_0_xyyyy_xxxyyy[k] = -g_z_0_yyyy_xxxyyy[k] * ab_x + g_z_0_yyyy_xxxxyyy[k];

                g_z_0_xyyyy_xxxyyz[k] = -g_z_0_yyyy_xxxyyz[k] * ab_x + g_z_0_yyyy_xxxxyyz[k];

                g_z_0_xyyyy_xxxyzz[k] = -g_z_0_yyyy_xxxyzz[k] * ab_x + g_z_0_yyyy_xxxxyzz[k];

                g_z_0_xyyyy_xxxzzz[k] = -g_z_0_yyyy_xxxzzz[k] * ab_x + g_z_0_yyyy_xxxxzzz[k];

                g_z_0_xyyyy_xxyyyy[k] = -g_z_0_yyyy_xxyyyy[k] * ab_x + g_z_0_yyyy_xxxyyyy[k];

                g_z_0_xyyyy_xxyyyz[k] = -g_z_0_yyyy_xxyyyz[k] * ab_x + g_z_0_yyyy_xxxyyyz[k];

                g_z_0_xyyyy_xxyyzz[k] = -g_z_0_yyyy_xxyyzz[k] * ab_x + g_z_0_yyyy_xxxyyzz[k];

                g_z_0_xyyyy_xxyzzz[k] = -g_z_0_yyyy_xxyzzz[k] * ab_x + g_z_0_yyyy_xxxyzzz[k];

                g_z_0_xyyyy_xxzzzz[k] = -g_z_0_yyyy_xxzzzz[k] * ab_x + g_z_0_yyyy_xxxzzzz[k];

                g_z_0_xyyyy_xyyyyy[k] = -g_z_0_yyyy_xyyyyy[k] * ab_x + g_z_0_yyyy_xxyyyyy[k];

                g_z_0_xyyyy_xyyyyz[k] = -g_z_0_yyyy_xyyyyz[k] * ab_x + g_z_0_yyyy_xxyyyyz[k];

                g_z_0_xyyyy_xyyyzz[k] = -g_z_0_yyyy_xyyyzz[k] * ab_x + g_z_0_yyyy_xxyyyzz[k];

                g_z_0_xyyyy_xyyzzz[k] = -g_z_0_yyyy_xyyzzz[k] * ab_x + g_z_0_yyyy_xxyyzzz[k];

                g_z_0_xyyyy_xyzzzz[k] = -g_z_0_yyyy_xyzzzz[k] * ab_x + g_z_0_yyyy_xxyzzzz[k];

                g_z_0_xyyyy_xzzzzz[k] = -g_z_0_yyyy_xzzzzz[k] * ab_x + g_z_0_yyyy_xxzzzzz[k];

                g_z_0_xyyyy_yyyyyy[k] = -g_z_0_yyyy_yyyyyy[k] * ab_x + g_z_0_yyyy_xyyyyyy[k];

                g_z_0_xyyyy_yyyyyz[k] = -g_z_0_yyyy_yyyyyz[k] * ab_x + g_z_0_yyyy_xyyyyyz[k];

                g_z_0_xyyyy_yyyyzz[k] = -g_z_0_yyyy_yyyyzz[k] * ab_x + g_z_0_yyyy_xyyyyzz[k];

                g_z_0_xyyyy_yyyzzz[k] = -g_z_0_yyyy_yyyzzz[k] * ab_x + g_z_0_yyyy_xyyyzzz[k];

                g_z_0_xyyyy_yyzzzz[k] = -g_z_0_yyyy_yyzzzz[k] * ab_x + g_z_0_yyyy_xyyzzzz[k];

                g_z_0_xyyyy_yzzzzz[k] = -g_z_0_yyyy_yzzzzz[k] * ab_x + g_z_0_yyyy_xyzzzzz[k];

                g_z_0_xyyyy_zzzzzz[k] = -g_z_0_yyyy_zzzzzz[k] * ab_x + g_z_0_yyyy_xzzzzzz[k];
            }

            /// Set up 1484-1512 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1484 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1485 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1486 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1487 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1488 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1489 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1490 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1491 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1492 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1493 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1494 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1495 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1496 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1497 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1498 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1499 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1500 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1501 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1502 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1503 * ccomps * dcomps);

            auto g_z_0_xyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1504 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1505 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1506 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1507 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1508 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1509 * ccomps * dcomps);

            auto g_z_0_xyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1510 * ccomps * dcomps);

            auto g_z_0_xyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1511 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyz_xxxxxx, g_z_0_xyyyz_xxxxxy, g_z_0_xyyyz_xxxxxz, g_z_0_xyyyz_xxxxyy, g_z_0_xyyyz_xxxxyz, g_z_0_xyyyz_xxxxzz, g_z_0_xyyyz_xxxyyy, g_z_0_xyyyz_xxxyyz, g_z_0_xyyyz_xxxyzz, g_z_0_xyyyz_xxxzzz, g_z_0_xyyyz_xxyyyy, g_z_0_xyyyz_xxyyyz, g_z_0_xyyyz_xxyyzz, g_z_0_xyyyz_xxyzzz, g_z_0_xyyyz_xxzzzz, g_z_0_xyyyz_xyyyyy, g_z_0_xyyyz_xyyyyz, g_z_0_xyyyz_xyyyzz, g_z_0_xyyyz_xyyzzz, g_z_0_xyyyz_xyzzzz, g_z_0_xyyyz_xzzzzz, g_z_0_xyyyz_yyyyyy, g_z_0_xyyyz_yyyyyz, g_z_0_xyyyz_yyyyzz, g_z_0_xyyyz_yyyzzz, g_z_0_xyyyz_yyzzzz, g_z_0_xyyyz_yzzzzz, g_z_0_xyyyz_zzzzzz, g_z_0_yyyz_xxxxxx, g_z_0_yyyz_xxxxxxx, g_z_0_yyyz_xxxxxxy, g_z_0_yyyz_xxxxxxz, g_z_0_yyyz_xxxxxy, g_z_0_yyyz_xxxxxyy, g_z_0_yyyz_xxxxxyz, g_z_0_yyyz_xxxxxz, g_z_0_yyyz_xxxxxzz, g_z_0_yyyz_xxxxyy, g_z_0_yyyz_xxxxyyy, g_z_0_yyyz_xxxxyyz, g_z_0_yyyz_xxxxyz, g_z_0_yyyz_xxxxyzz, g_z_0_yyyz_xxxxzz, g_z_0_yyyz_xxxxzzz, g_z_0_yyyz_xxxyyy, g_z_0_yyyz_xxxyyyy, g_z_0_yyyz_xxxyyyz, g_z_0_yyyz_xxxyyz, g_z_0_yyyz_xxxyyzz, g_z_0_yyyz_xxxyzz, g_z_0_yyyz_xxxyzzz, g_z_0_yyyz_xxxzzz, g_z_0_yyyz_xxxzzzz, g_z_0_yyyz_xxyyyy, g_z_0_yyyz_xxyyyyy, g_z_0_yyyz_xxyyyyz, g_z_0_yyyz_xxyyyz, g_z_0_yyyz_xxyyyzz, g_z_0_yyyz_xxyyzz, g_z_0_yyyz_xxyyzzz, g_z_0_yyyz_xxyzzz, g_z_0_yyyz_xxyzzzz, g_z_0_yyyz_xxzzzz, g_z_0_yyyz_xxzzzzz, g_z_0_yyyz_xyyyyy, g_z_0_yyyz_xyyyyyy, g_z_0_yyyz_xyyyyyz, g_z_0_yyyz_xyyyyz, g_z_0_yyyz_xyyyyzz, g_z_0_yyyz_xyyyzz, g_z_0_yyyz_xyyyzzz, g_z_0_yyyz_xyyzzz, g_z_0_yyyz_xyyzzzz, g_z_0_yyyz_xyzzzz, g_z_0_yyyz_xyzzzzz, g_z_0_yyyz_xzzzzz, g_z_0_yyyz_xzzzzzz, g_z_0_yyyz_yyyyyy, g_z_0_yyyz_yyyyyz, g_z_0_yyyz_yyyyzz, g_z_0_yyyz_yyyzzz, g_z_0_yyyz_yyzzzz, g_z_0_yyyz_yzzzzz, g_z_0_yyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyz_xxxxxx[k] = -g_z_0_yyyz_xxxxxx[k] * ab_x + g_z_0_yyyz_xxxxxxx[k];

                g_z_0_xyyyz_xxxxxy[k] = -g_z_0_yyyz_xxxxxy[k] * ab_x + g_z_0_yyyz_xxxxxxy[k];

                g_z_0_xyyyz_xxxxxz[k] = -g_z_0_yyyz_xxxxxz[k] * ab_x + g_z_0_yyyz_xxxxxxz[k];

                g_z_0_xyyyz_xxxxyy[k] = -g_z_0_yyyz_xxxxyy[k] * ab_x + g_z_0_yyyz_xxxxxyy[k];

                g_z_0_xyyyz_xxxxyz[k] = -g_z_0_yyyz_xxxxyz[k] * ab_x + g_z_0_yyyz_xxxxxyz[k];

                g_z_0_xyyyz_xxxxzz[k] = -g_z_0_yyyz_xxxxzz[k] * ab_x + g_z_0_yyyz_xxxxxzz[k];

                g_z_0_xyyyz_xxxyyy[k] = -g_z_0_yyyz_xxxyyy[k] * ab_x + g_z_0_yyyz_xxxxyyy[k];

                g_z_0_xyyyz_xxxyyz[k] = -g_z_0_yyyz_xxxyyz[k] * ab_x + g_z_0_yyyz_xxxxyyz[k];

                g_z_0_xyyyz_xxxyzz[k] = -g_z_0_yyyz_xxxyzz[k] * ab_x + g_z_0_yyyz_xxxxyzz[k];

                g_z_0_xyyyz_xxxzzz[k] = -g_z_0_yyyz_xxxzzz[k] * ab_x + g_z_0_yyyz_xxxxzzz[k];

                g_z_0_xyyyz_xxyyyy[k] = -g_z_0_yyyz_xxyyyy[k] * ab_x + g_z_0_yyyz_xxxyyyy[k];

                g_z_0_xyyyz_xxyyyz[k] = -g_z_0_yyyz_xxyyyz[k] * ab_x + g_z_0_yyyz_xxxyyyz[k];

                g_z_0_xyyyz_xxyyzz[k] = -g_z_0_yyyz_xxyyzz[k] * ab_x + g_z_0_yyyz_xxxyyzz[k];

                g_z_0_xyyyz_xxyzzz[k] = -g_z_0_yyyz_xxyzzz[k] * ab_x + g_z_0_yyyz_xxxyzzz[k];

                g_z_0_xyyyz_xxzzzz[k] = -g_z_0_yyyz_xxzzzz[k] * ab_x + g_z_0_yyyz_xxxzzzz[k];

                g_z_0_xyyyz_xyyyyy[k] = -g_z_0_yyyz_xyyyyy[k] * ab_x + g_z_0_yyyz_xxyyyyy[k];

                g_z_0_xyyyz_xyyyyz[k] = -g_z_0_yyyz_xyyyyz[k] * ab_x + g_z_0_yyyz_xxyyyyz[k];

                g_z_0_xyyyz_xyyyzz[k] = -g_z_0_yyyz_xyyyzz[k] * ab_x + g_z_0_yyyz_xxyyyzz[k];

                g_z_0_xyyyz_xyyzzz[k] = -g_z_0_yyyz_xyyzzz[k] * ab_x + g_z_0_yyyz_xxyyzzz[k];

                g_z_0_xyyyz_xyzzzz[k] = -g_z_0_yyyz_xyzzzz[k] * ab_x + g_z_0_yyyz_xxyzzzz[k];

                g_z_0_xyyyz_xzzzzz[k] = -g_z_0_yyyz_xzzzzz[k] * ab_x + g_z_0_yyyz_xxzzzzz[k];

                g_z_0_xyyyz_yyyyyy[k] = -g_z_0_yyyz_yyyyyy[k] * ab_x + g_z_0_yyyz_xyyyyyy[k];

                g_z_0_xyyyz_yyyyyz[k] = -g_z_0_yyyz_yyyyyz[k] * ab_x + g_z_0_yyyz_xyyyyyz[k];

                g_z_0_xyyyz_yyyyzz[k] = -g_z_0_yyyz_yyyyzz[k] * ab_x + g_z_0_yyyz_xyyyyzz[k];

                g_z_0_xyyyz_yyyzzz[k] = -g_z_0_yyyz_yyyzzz[k] * ab_x + g_z_0_yyyz_xyyyzzz[k];

                g_z_0_xyyyz_yyzzzz[k] = -g_z_0_yyyz_yyzzzz[k] * ab_x + g_z_0_yyyz_xyyzzzz[k];

                g_z_0_xyyyz_yzzzzz[k] = -g_z_0_yyyz_yzzzzz[k] * ab_x + g_z_0_yyyz_xyzzzzz[k];

                g_z_0_xyyyz_zzzzzz[k] = -g_z_0_yyyz_zzzzzz[k] * ab_x + g_z_0_yyyz_xzzzzzz[k];
            }

            /// Set up 1512-1540 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1512 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1513 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1514 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1515 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1516 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1517 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1518 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1519 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1520 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1521 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1522 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1523 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1524 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1525 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1526 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1527 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1528 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1529 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1530 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1531 * ccomps * dcomps);

            auto g_z_0_xyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1532 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1533 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1534 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1535 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1536 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1537 * ccomps * dcomps);

            auto g_z_0_xyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1538 * ccomps * dcomps);

            auto g_z_0_xyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1539 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyzz_xxxxxx, g_z_0_xyyzz_xxxxxy, g_z_0_xyyzz_xxxxxz, g_z_0_xyyzz_xxxxyy, g_z_0_xyyzz_xxxxyz, g_z_0_xyyzz_xxxxzz, g_z_0_xyyzz_xxxyyy, g_z_0_xyyzz_xxxyyz, g_z_0_xyyzz_xxxyzz, g_z_0_xyyzz_xxxzzz, g_z_0_xyyzz_xxyyyy, g_z_0_xyyzz_xxyyyz, g_z_0_xyyzz_xxyyzz, g_z_0_xyyzz_xxyzzz, g_z_0_xyyzz_xxzzzz, g_z_0_xyyzz_xyyyyy, g_z_0_xyyzz_xyyyyz, g_z_0_xyyzz_xyyyzz, g_z_0_xyyzz_xyyzzz, g_z_0_xyyzz_xyzzzz, g_z_0_xyyzz_xzzzzz, g_z_0_xyyzz_yyyyyy, g_z_0_xyyzz_yyyyyz, g_z_0_xyyzz_yyyyzz, g_z_0_xyyzz_yyyzzz, g_z_0_xyyzz_yyzzzz, g_z_0_xyyzz_yzzzzz, g_z_0_xyyzz_zzzzzz, g_z_0_yyzz_xxxxxx, g_z_0_yyzz_xxxxxxx, g_z_0_yyzz_xxxxxxy, g_z_0_yyzz_xxxxxxz, g_z_0_yyzz_xxxxxy, g_z_0_yyzz_xxxxxyy, g_z_0_yyzz_xxxxxyz, g_z_0_yyzz_xxxxxz, g_z_0_yyzz_xxxxxzz, g_z_0_yyzz_xxxxyy, g_z_0_yyzz_xxxxyyy, g_z_0_yyzz_xxxxyyz, g_z_0_yyzz_xxxxyz, g_z_0_yyzz_xxxxyzz, g_z_0_yyzz_xxxxzz, g_z_0_yyzz_xxxxzzz, g_z_0_yyzz_xxxyyy, g_z_0_yyzz_xxxyyyy, g_z_0_yyzz_xxxyyyz, g_z_0_yyzz_xxxyyz, g_z_0_yyzz_xxxyyzz, g_z_0_yyzz_xxxyzz, g_z_0_yyzz_xxxyzzz, g_z_0_yyzz_xxxzzz, g_z_0_yyzz_xxxzzzz, g_z_0_yyzz_xxyyyy, g_z_0_yyzz_xxyyyyy, g_z_0_yyzz_xxyyyyz, g_z_0_yyzz_xxyyyz, g_z_0_yyzz_xxyyyzz, g_z_0_yyzz_xxyyzz, g_z_0_yyzz_xxyyzzz, g_z_0_yyzz_xxyzzz, g_z_0_yyzz_xxyzzzz, g_z_0_yyzz_xxzzzz, g_z_0_yyzz_xxzzzzz, g_z_0_yyzz_xyyyyy, g_z_0_yyzz_xyyyyyy, g_z_0_yyzz_xyyyyyz, g_z_0_yyzz_xyyyyz, g_z_0_yyzz_xyyyyzz, g_z_0_yyzz_xyyyzz, g_z_0_yyzz_xyyyzzz, g_z_0_yyzz_xyyzzz, g_z_0_yyzz_xyyzzzz, g_z_0_yyzz_xyzzzz, g_z_0_yyzz_xyzzzzz, g_z_0_yyzz_xzzzzz, g_z_0_yyzz_xzzzzzz, g_z_0_yyzz_yyyyyy, g_z_0_yyzz_yyyyyz, g_z_0_yyzz_yyyyzz, g_z_0_yyzz_yyyzzz, g_z_0_yyzz_yyzzzz, g_z_0_yyzz_yzzzzz, g_z_0_yyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzz_xxxxxx[k] = -g_z_0_yyzz_xxxxxx[k] * ab_x + g_z_0_yyzz_xxxxxxx[k];

                g_z_0_xyyzz_xxxxxy[k] = -g_z_0_yyzz_xxxxxy[k] * ab_x + g_z_0_yyzz_xxxxxxy[k];

                g_z_0_xyyzz_xxxxxz[k] = -g_z_0_yyzz_xxxxxz[k] * ab_x + g_z_0_yyzz_xxxxxxz[k];

                g_z_0_xyyzz_xxxxyy[k] = -g_z_0_yyzz_xxxxyy[k] * ab_x + g_z_0_yyzz_xxxxxyy[k];

                g_z_0_xyyzz_xxxxyz[k] = -g_z_0_yyzz_xxxxyz[k] * ab_x + g_z_0_yyzz_xxxxxyz[k];

                g_z_0_xyyzz_xxxxzz[k] = -g_z_0_yyzz_xxxxzz[k] * ab_x + g_z_0_yyzz_xxxxxzz[k];

                g_z_0_xyyzz_xxxyyy[k] = -g_z_0_yyzz_xxxyyy[k] * ab_x + g_z_0_yyzz_xxxxyyy[k];

                g_z_0_xyyzz_xxxyyz[k] = -g_z_0_yyzz_xxxyyz[k] * ab_x + g_z_0_yyzz_xxxxyyz[k];

                g_z_0_xyyzz_xxxyzz[k] = -g_z_0_yyzz_xxxyzz[k] * ab_x + g_z_0_yyzz_xxxxyzz[k];

                g_z_0_xyyzz_xxxzzz[k] = -g_z_0_yyzz_xxxzzz[k] * ab_x + g_z_0_yyzz_xxxxzzz[k];

                g_z_0_xyyzz_xxyyyy[k] = -g_z_0_yyzz_xxyyyy[k] * ab_x + g_z_0_yyzz_xxxyyyy[k];

                g_z_0_xyyzz_xxyyyz[k] = -g_z_0_yyzz_xxyyyz[k] * ab_x + g_z_0_yyzz_xxxyyyz[k];

                g_z_0_xyyzz_xxyyzz[k] = -g_z_0_yyzz_xxyyzz[k] * ab_x + g_z_0_yyzz_xxxyyzz[k];

                g_z_0_xyyzz_xxyzzz[k] = -g_z_0_yyzz_xxyzzz[k] * ab_x + g_z_0_yyzz_xxxyzzz[k];

                g_z_0_xyyzz_xxzzzz[k] = -g_z_0_yyzz_xxzzzz[k] * ab_x + g_z_0_yyzz_xxxzzzz[k];

                g_z_0_xyyzz_xyyyyy[k] = -g_z_0_yyzz_xyyyyy[k] * ab_x + g_z_0_yyzz_xxyyyyy[k];

                g_z_0_xyyzz_xyyyyz[k] = -g_z_0_yyzz_xyyyyz[k] * ab_x + g_z_0_yyzz_xxyyyyz[k];

                g_z_0_xyyzz_xyyyzz[k] = -g_z_0_yyzz_xyyyzz[k] * ab_x + g_z_0_yyzz_xxyyyzz[k];

                g_z_0_xyyzz_xyyzzz[k] = -g_z_0_yyzz_xyyzzz[k] * ab_x + g_z_0_yyzz_xxyyzzz[k];

                g_z_0_xyyzz_xyzzzz[k] = -g_z_0_yyzz_xyzzzz[k] * ab_x + g_z_0_yyzz_xxyzzzz[k];

                g_z_0_xyyzz_xzzzzz[k] = -g_z_0_yyzz_xzzzzz[k] * ab_x + g_z_0_yyzz_xxzzzzz[k];

                g_z_0_xyyzz_yyyyyy[k] = -g_z_0_yyzz_yyyyyy[k] * ab_x + g_z_0_yyzz_xyyyyyy[k];

                g_z_0_xyyzz_yyyyyz[k] = -g_z_0_yyzz_yyyyyz[k] * ab_x + g_z_0_yyzz_xyyyyyz[k];

                g_z_0_xyyzz_yyyyzz[k] = -g_z_0_yyzz_yyyyzz[k] * ab_x + g_z_0_yyzz_xyyyyzz[k];

                g_z_0_xyyzz_yyyzzz[k] = -g_z_0_yyzz_yyyzzz[k] * ab_x + g_z_0_yyzz_xyyyzzz[k];

                g_z_0_xyyzz_yyzzzz[k] = -g_z_0_yyzz_yyzzzz[k] * ab_x + g_z_0_yyzz_xyyzzzz[k];

                g_z_0_xyyzz_yzzzzz[k] = -g_z_0_yyzz_yzzzzz[k] * ab_x + g_z_0_yyzz_xyzzzzz[k];

                g_z_0_xyyzz_zzzzzz[k] = -g_z_0_yyzz_zzzzzz[k] * ab_x + g_z_0_yyzz_xzzzzzz[k];
            }

            /// Set up 1540-1568 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1540 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1541 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1542 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1543 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1544 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1545 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1546 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1547 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1548 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1549 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1550 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1551 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1552 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1553 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1554 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1555 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1556 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1557 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1558 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1559 * ccomps * dcomps);

            auto g_z_0_xyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1560 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1561 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1562 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1563 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1564 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1565 * ccomps * dcomps);

            auto g_z_0_xyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1566 * ccomps * dcomps);

            auto g_z_0_xyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1567 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzzz_xxxxxx, g_z_0_xyzzz_xxxxxy, g_z_0_xyzzz_xxxxxz, g_z_0_xyzzz_xxxxyy, g_z_0_xyzzz_xxxxyz, g_z_0_xyzzz_xxxxzz, g_z_0_xyzzz_xxxyyy, g_z_0_xyzzz_xxxyyz, g_z_0_xyzzz_xxxyzz, g_z_0_xyzzz_xxxzzz, g_z_0_xyzzz_xxyyyy, g_z_0_xyzzz_xxyyyz, g_z_0_xyzzz_xxyyzz, g_z_0_xyzzz_xxyzzz, g_z_0_xyzzz_xxzzzz, g_z_0_xyzzz_xyyyyy, g_z_0_xyzzz_xyyyyz, g_z_0_xyzzz_xyyyzz, g_z_0_xyzzz_xyyzzz, g_z_0_xyzzz_xyzzzz, g_z_0_xyzzz_xzzzzz, g_z_0_xyzzz_yyyyyy, g_z_0_xyzzz_yyyyyz, g_z_0_xyzzz_yyyyzz, g_z_0_xyzzz_yyyzzz, g_z_0_xyzzz_yyzzzz, g_z_0_xyzzz_yzzzzz, g_z_0_xyzzz_zzzzzz, g_z_0_yzzz_xxxxxx, g_z_0_yzzz_xxxxxxx, g_z_0_yzzz_xxxxxxy, g_z_0_yzzz_xxxxxxz, g_z_0_yzzz_xxxxxy, g_z_0_yzzz_xxxxxyy, g_z_0_yzzz_xxxxxyz, g_z_0_yzzz_xxxxxz, g_z_0_yzzz_xxxxxzz, g_z_0_yzzz_xxxxyy, g_z_0_yzzz_xxxxyyy, g_z_0_yzzz_xxxxyyz, g_z_0_yzzz_xxxxyz, g_z_0_yzzz_xxxxyzz, g_z_0_yzzz_xxxxzz, g_z_0_yzzz_xxxxzzz, g_z_0_yzzz_xxxyyy, g_z_0_yzzz_xxxyyyy, g_z_0_yzzz_xxxyyyz, g_z_0_yzzz_xxxyyz, g_z_0_yzzz_xxxyyzz, g_z_0_yzzz_xxxyzz, g_z_0_yzzz_xxxyzzz, g_z_0_yzzz_xxxzzz, g_z_0_yzzz_xxxzzzz, g_z_0_yzzz_xxyyyy, g_z_0_yzzz_xxyyyyy, g_z_0_yzzz_xxyyyyz, g_z_0_yzzz_xxyyyz, g_z_0_yzzz_xxyyyzz, g_z_0_yzzz_xxyyzz, g_z_0_yzzz_xxyyzzz, g_z_0_yzzz_xxyzzz, g_z_0_yzzz_xxyzzzz, g_z_0_yzzz_xxzzzz, g_z_0_yzzz_xxzzzzz, g_z_0_yzzz_xyyyyy, g_z_0_yzzz_xyyyyyy, g_z_0_yzzz_xyyyyyz, g_z_0_yzzz_xyyyyz, g_z_0_yzzz_xyyyyzz, g_z_0_yzzz_xyyyzz, g_z_0_yzzz_xyyyzzz, g_z_0_yzzz_xyyzzz, g_z_0_yzzz_xyyzzzz, g_z_0_yzzz_xyzzzz, g_z_0_yzzz_xyzzzzz, g_z_0_yzzz_xzzzzz, g_z_0_yzzz_xzzzzzz, g_z_0_yzzz_yyyyyy, g_z_0_yzzz_yyyyyz, g_z_0_yzzz_yyyyzz, g_z_0_yzzz_yyyzzz, g_z_0_yzzz_yyzzzz, g_z_0_yzzz_yzzzzz, g_z_0_yzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzz_xxxxxx[k] = -g_z_0_yzzz_xxxxxx[k] * ab_x + g_z_0_yzzz_xxxxxxx[k];

                g_z_0_xyzzz_xxxxxy[k] = -g_z_0_yzzz_xxxxxy[k] * ab_x + g_z_0_yzzz_xxxxxxy[k];

                g_z_0_xyzzz_xxxxxz[k] = -g_z_0_yzzz_xxxxxz[k] * ab_x + g_z_0_yzzz_xxxxxxz[k];

                g_z_0_xyzzz_xxxxyy[k] = -g_z_0_yzzz_xxxxyy[k] * ab_x + g_z_0_yzzz_xxxxxyy[k];

                g_z_0_xyzzz_xxxxyz[k] = -g_z_0_yzzz_xxxxyz[k] * ab_x + g_z_0_yzzz_xxxxxyz[k];

                g_z_0_xyzzz_xxxxzz[k] = -g_z_0_yzzz_xxxxzz[k] * ab_x + g_z_0_yzzz_xxxxxzz[k];

                g_z_0_xyzzz_xxxyyy[k] = -g_z_0_yzzz_xxxyyy[k] * ab_x + g_z_0_yzzz_xxxxyyy[k];

                g_z_0_xyzzz_xxxyyz[k] = -g_z_0_yzzz_xxxyyz[k] * ab_x + g_z_0_yzzz_xxxxyyz[k];

                g_z_0_xyzzz_xxxyzz[k] = -g_z_0_yzzz_xxxyzz[k] * ab_x + g_z_0_yzzz_xxxxyzz[k];

                g_z_0_xyzzz_xxxzzz[k] = -g_z_0_yzzz_xxxzzz[k] * ab_x + g_z_0_yzzz_xxxxzzz[k];

                g_z_0_xyzzz_xxyyyy[k] = -g_z_0_yzzz_xxyyyy[k] * ab_x + g_z_0_yzzz_xxxyyyy[k];

                g_z_0_xyzzz_xxyyyz[k] = -g_z_0_yzzz_xxyyyz[k] * ab_x + g_z_0_yzzz_xxxyyyz[k];

                g_z_0_xyzzz_xxyyzz[k] = -g_z_0_yzzz_xxyyzz[k] * ab_x + g_z_0_yzzz_xxxyyzz[k];

                g_z_0_xyzzz_xxyzzz[k] = -g_z_0_yzzz_xxyzzz[k] * ab_x + g_z_0_yzzz_xxxyzzz[k];

                g_z_0_xyzzz_xxzzzz[k] = -g_z_0_yzzz_xxzzzz[k] * ab_x + g_z_0_yzzz_xxxzzzz[k];

                g_z_0_xyzzz_xyyyyy[k] = -g_z_0_yzzz_xyyyyy[k] * ab_x + g_z_0_yzzz_xxyyyyy[k];

                g_z_0_xyzzz_xyyyyz[k] = -g_z_0_yzzz_xyyyyz[k] * ab_x + g_z_0_yzzz_xxyyyyz[k];

                g_z_0_xyzzz_xyyyzz[k] = -g_z_0_yzzz_xyyyzz[k] * ab_x + g_z_0_yzzz_xxyyyzz[k];

                g_z_0_xyzzz_xyyzzz[k] = -g_z_0_yzzz_xyyzzz[k] * ab_x + g_z_0_yzzz_xxyyzzz[k];

                g_z_0_xyzzz_xyzzzz[k] = -g_z_0_yzzz_xyzzzz[k] * ab_x + g_z_0_yzzz_xxyzzzz[k];

                g_z_0_xyzzz_xzzzzz[k] = -g_z_0_yzzz_xzzzzz[k] * ab_x + g_z_0_yzzz_xxzzzzz[k];

                g_z_0_xyzzz_yyyyyy[k] = -g_z_0_yzzz_yyyyyy[k] * ab_x + g_z_0_yzzz_xyyyyyy[k];

                g_z_0_xyzzz_yyyyyz[k] = -g_z_0_yzzz_yyyyyz[k] * ab_x + g_z_0_yzzz_xyyyyyz[k];

                g_z_0_xyzzz_yyyyzz[k] = -g_z_0_yzzz_yyyyzz[k] * ab_x + g_z_0_yzzz_xyyyyzz[k];

                g_z_0_xyzzz_yyyzzz[k] = -g_z_0_yzzz_yyyzzz[k] * ab_x + g_z_0_yzzz_xyyyzzz[k];

                g_z_0_xyzzz_yyzzzz[k] = -g_z_0_yzzz_yyzzzz[k] * ab_x + g_z_0_yzzz_xyyzzzz[k];

                g_z_0_xyzzz_yzzzzz[k] = -g_z_0_yzzz_yzzzzz[k] * ab_x + g_z_0_yzzz_xyzzzzz[k];

                g_z_0_xyzzz_zzzzzz[k] = -g_z_0_yzzz_zzzzzz[k] * ab_x + g_z_0_yzzz_xzzzzzz[k];
            }

            /// Set up 1568-1596 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1568 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1569 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1570 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1571 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1572 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1573 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1574 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1575 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1576 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1577 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1578 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1579 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1580 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1581 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1582 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1583 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1584 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1585 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1586 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1587 * ccomps * dcomps);

            auto g_z_0_xzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1588 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1589 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1590 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1591 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1592 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1593 * ccomps * dcomps);

            auto g_z_0_xzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1594 * ccomps * dcomps);

            auto g_z_0_xzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1595 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzzz_xxxxxx, g_z_0_xzzzz_xxxxxy, g_z_0_xzzzz_xxxxxz, g_z_0_xzzzz_xxxxyy, g_z_0_xzzzz_xxxxyz, g_z_0_xzzzz_xxxxzz, g_z_0_xzzzz_xxxyyy, g_z_0_xzzzz_xxxyyz, g_z_0_xzzzz_xxxyzz, g_z_0_xzzzz_xxxzzz, g_z_0_xzzzz_xxyyyy, g_z_0_xzzzz_xxyyyz, g_z_0_xzzzz_xxyyzz, g_z_0_xzzzz_xxyzzz, g_z_0_xzzzz_xxzzzz, g_z_0_xzzzz_xyyyyy, g_z_0_xzzzz_xyyyyz, g_z_0_xzzzz_xyyyzz, g_z_0_xzzzz_xyyzzz, g_z_0_xzzzz_xyzzzz, g_z_0_xzzzz_xzzzzz, g_z_0_xzzzz_yyyyyy, g_z_0_xzzzz_yyyyyz, g_z_0_xzzzz_yyyyzz, g_z_0_xzzzz_yyyzzz, g_z_0_xzzzz_yyzzzz, g_z_0_xzzzz_yzzzzz, g_z_0_xzzzz_zzzzzz, g_z_0_zzzz_xxxxxx, g_z_0_zzzz_xxxxxxx, g_z_0_zzzz_xxxxxxy, g_z_0_zzzz_xxxxxxz, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxxyy, g_z_0_zzzz_xxxxxyz, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxxzz, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyyy, g_z_0_zzzz_xxxxyyz, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxyzz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxxzzz, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyyy, g_z_0_zzzz_xxxyyyz, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyyzz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxyzzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxxzzzz, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyyy, g_z_0_zzzz_xxyyyyz, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyyzz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyyzzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxyzzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xxzzzzz, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyyy, g_z_0_zzzz_xyyyyyz, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyyzz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyyzzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyyzzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xyzzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_xzzzzzz, g_z_0_zzzz_yyyyyy, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzz_xxxxxx[k] = -g_z_0_zzzz_xxxxxx[k] * ab_x + g_z_0_zzzz_xxxxxxx[k];

                g_z_0_xzzzz_xxxxxy[k] = -g_z_0_zzzz_xxxxxy[k] * ab_x + g_z_0_zzzz_xxxxxxy[k];

                g_z_0_xzzzz_xxxxxz[k] = -g_z_0_zzzz_xxxxxz[k] * ab_x + g_z_0_zzzz_xxxxxxz[k];

                g_z_0_xzzzz_xxxxyy[k] = -g_z_0_zzzz_xxxxyy[k] * ab_x + g_z_0_zzzz_xxxxxyy[k];

                g_z_0_xzzzz_xxxxyz[k] = -g_z_0_zzzz_xxxxyz[k] * ab_x + g_z_0_zzzz_xxxxxyz[k];

                g_z_0_xzzzz_xxxxzz[k] = -g_z_0_zzzz_xxxxzz[k] * ab_x + g_z_0_zzzz_xxxxxzz[k];

                g_z_0_xzzzz_xxxyyy[k] = -g_z_0_zzzz_xxxyyy[k] * ab_x + g_z_0_zzzz_xxxxyyy[k];

                g_z_0_xzzzz_xxxyyz[k] = -g_z_0_zzzz_xxxyyz[k] * ab_x + g_z_0_zzzz_xxxxyyz[k];

                g_z_0_xzzzz_xxxyzz[k] = -g_z_0_zzzz_xxxyzz[k] * ab_x + g_z_0_zzzz_xxxxyzz[k];

                g_z_0_xzzzz_xxxzzz[k] = -g_z_0_zzzz_xxxzzz[k] * ab_x + g_z_0_zzzz_xxxxzzz[k];

                g_z_0_xzzzz_xxyyyy[k] = -g_z_0_zzzz_xxyyyy[k] * ab_x + g_z_0_zzzz_xxxyyyy[k];

                g_z_0_xzzzz_xxyyyz[k] = -g_z_0_zzzz_xxyyyz[k] * ab_x + g_z_0_zzzz_xxxyyyz[k];

                g_z_0_xzzzz_xxyyzz[k] = -g_z_0_zzzz_xxyyzz[k] * ab_x + g_z_0_zzzz_xxxyyzz[k];

                g_z_0_xzzzz_xxyzzz[k] = -g_z_0_zzzz_xxyzzz[k] * ab_x + g_z_0_zzzz_xxxyzzz[k];

                g_z_0_xzzzz_xxzzzz[k] = -g_z_0_zzzz_xxzzzz[k] * ab_x + g_z_0_zzzz_xxxzzzz[k];

                g_z_0_xzzzz_xyyyyy[k] = -g_z_0_zzzz_xyyyyy[k] * ab_x + g_z_0_zzzz_xxyyyyy[k];

                g_z_0_xzzzz_xyyyyz[k] = -g_z_0_zzzz_xyyyyz[k] * ab_x + g_z_0_zzzz_xxyyyyz[k];

                g_z_0_xzzzz_xyyyzz[k] = -g_z_0_zzzz_xyyyzz[k] * ab_x + g_z_0_zzzz_xxyyyzz[k];

                g_z_0_xzzzz_xyyzzz[k] = -g_z_0_zzzz_xyyzzz[k] * ab_x + g_z_0_zzzz_xxyyzzz[k];

                g_z_0_xzzzz_xyzzzz[k] = -g_z_0_zzzz_xyzzzz[k] * ab_x + g_z_0_zzzz_xxyzzzz[k];

                g_z_0_xzzzz_xzzzzz[k] = -g_z_0_zzzz_xzzzzz[k] * ab_x + g_z_0_zzzz_xxzzzzz[k];

                g_z_0_xzzzz_yyyyyy[k] = -g_z_0_zzzz_yyyyyy[k] * ab_x + g_z_0_zzzz_xyyyyyy[k];

                g_z_0_xzzzz_yyyyyz[k] = -g_z_0_zzzz_yyyyyz[k] * ab_x + g_z_0_zzzz_xyyyyyz[k];

                g_z_0_xzzzz_yyyyzz[k] = -g_z_0_zzzz_yyyyzz[k] * ab_x + g_z_0_zzzz_xyyyyzz[k];

                g_z_0_xzzzz_yyyzzz[k] = -g_z_0_zzzz_yyyzzz[k] * ab_x + g_z_0_zzzz_xyyyzzz[k];

                g_z_0_xzzzz_yyzzzz[k] = -g_z_0_zzzz_yyzzzz[k] * ab_x + g_z_0_zzzz_xyyzzzz[k];

                g_z_0_xzzzz_yzzzzz[k] = -g_z_0_zzzz_yzzzzz[k] * ab_x + g_z_0_zzzz_xyzzzzz[k];

                g_z_0_xzzzz_zzzzzz[k] = -g_z_0_zzzz_zzzzzz[k] * ab_x + g_z_0_zzzz_xzzzzzz[k];
            }

            /// Set up 1596-1624 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyy_xxxxxx = cbuffer.data(hi_geom_10_off + 1596 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxy = cbuffer.data(hi_geom_10_off + 1597 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxz = cbuffer.data(hi_geom_10_off + 1598 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxyy = cbuffer.data(hi_geom_10_off + 1599 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxyz = cbuffer.data(hi_geom_10_off + 1600 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxzz = cbuffer.data(hi_geom_10_off + 1601 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyyy = cbuffer.data(hi_geom_10_off + 1602 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyyz = cbuffer.data(hi_geom_10_off + 1603 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyzz = cbuffer.data(hi_geom_10_off + 1604 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxzzz = cbuffer.data(hi_geom_10_off + 1605 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyyy = cbuffer.data(hi_geom_10_off + 1606 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyyz = cbuffer.data(hi_geom_10_off + 1607 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyzz = cbuffer.data(hi_geom_10_off + 1608 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyzzz = cbuffer.data(hi_geom_10_off + 1609 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxzzzz = cbuffer.data(hi_geom_10_off + 1610 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyyy = cbuffer.data(hi_geom_10_off + 1611 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyyz = cbuffer.data(hi_geom_10_off + 1612 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyzz = cbuffer.data(hi_geom_10_off + 1613 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyzzz = cbuffer.data(hi_geom_10_off + 1614 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyzzzz = cbuffer.data(hi_geom_10_off + 1615 * ccomps * dcomps);

            auto g_z_0_yyyyy_xzzzzz = cbuffer.data(hi_geom_10_off + 1616 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyyy = cbuffer.data(hi_geom_10_off + 1617 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyyz = cbuffer.data(hi_geom_10_off + 1618 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyzz = cbuffer.data(hi_geom_10_off + 1619 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyzzz = cbuffer.data(hi_geom_10_off + 1620 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyzzzz = cbuffer.data(hi_geom_10_off + 1621 * ccomps * dcomps);

            auto g_z_0_yyyyy_yzzzzz = cbuffer.data(hi_geom_10_off + 1622 * ccomps * dcomps);

            auto g_z_0_yyyyy_zzzzzz = cbuffer.data(hi_geom_10_off + 1623 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyy_xxxxxx, g_z_0_yyyy_xxxxxxy, g_z_0_yyyy_xxxxxy, g_z_0_yyyy_xxxxxyy, g_z_0_yyyy_xxxxxyz, g_z_0_yyyy_xxxxxz, g_z_0_yyyy_xxxxyy, g_z_0_yyyy_xxxxyyy, g_z_0_yyyy_xxxxyyz, g_z_0_yyyy_xxxxyz, g_z_0_yyyy_xxxxyzz, g_z_0_yyyy_xxxxzz, g_z_0_yyyy_xxxyyy, g_z_0_yyyy_xxxyyyy, g_z_0_yyyy_xxxyyyz, g_z_0_yyyy_xxxyyz, g_z_0_yyyy_xxxyyzz, g_z_0_yyyy_xxxyzz, g_z_0_yyyy_xxxyzzz, g_z_0_yyyy_xxxzzz, g_z_0_yyyy_xxyyyy, g_z_0_yyyy_xxyyyyy, g_z_0_yyyy_xxyyyyz, g_z_0_yyyy_xxyyyz, g_z_0_yyyy_xxyyyzz, g_z_0_yyyy_xxyyzz, g_z_0_yyyy_xxyyzzz, g_z_0_yyyy_xxyzzz, g_z_0_yyyy_xxyzzzz, g_z_0_yyyy_xxzzzz, g_z_0_yyyy_xyyyyy, g_z_0_yyyy_xyyyyyy, g_z_0_yyyy_xyyyyyz, g_z_0_yyyy_xyyyyz, g_z_0_yyyy_xyyyyzz, g_z_0_yyyy_xyyyzz, g_z_0_yyyy_xyyyzzz, g_z_0_yyyy_xyyzzz, g_z_0_yyyy_xyyzzzz, g_z_0_yyyy_xyzzzz, g_z_0_yyyy_xyzzzzz, g_z_0_yyyy_xzzzzz, g_z_0_yyyy_yyyyyy, g_z_0_yyyy_yyyyyyy, g_z_0_yyyy_yyyyyyz, g_z_0_yyyy_yyyyyz, g_z_0_yyyy_yyyyyzz, g_z_0_yyyy_yyyyzz, g_z_0_yyyy_yyyyzzz, g_z_0_yyyy_yyyzzz, g_z_0_yyyy_yyyzzzz, g_z_0_yyyy_yyzzzz, g_z_0_yyyy_yyzzzzz, g_z_0_yyyy_yzzzzz, g_z_0_yyyy_yzzzzzz, g_z_0_yyyy_zzzzzz, g_z_0_yyyyy_xxxxxx, g_z_0_yyyyy_xxxxxy, g_z_0_yyyyy_xxxxxz, g_z_0_yyyyy_xxxxyy, g_z_0_yyyyy_xxxxyz, g_z_0_yyyyy_xxxxzz, g_z_0_yyyyy_xxxyyy, g_z_0_yyyyy_xxxyyz, g_z_0_yyyyy_xxxyzz, g_z_0_yyyyy_xxxzzz, g_z_0_yyyyy_xxyyyy, g_z_0_yyyyy_xxyyyz, g_z_0_yyyyy_xxyyzz, g_z_0_yyyyy_xxyzzz, g_z_0_yyyyy_xxzzzz, g_z_0_yyyyy_xyyyyy, g_z_0_yyyyy_xyyyyz, g_z_0_yyyyy_xyyyzz, g_z_0_yyyyy_xyyzzz, g_z_0_yyyyy_xyzzzz, g_z_0_yyyyy_xzzzzz, g_z_0_yyyyy_yyyyyy, g_z_0_yyyyy_yyyyyz, g_z_0_yyyyy_yyyyzz, g_z_0_yyyyy_yyyzzz, g_z_0_yyyyy_yyzzzz, g_z_0_yyyyy_yzzzzz, g_z_0_yyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyy_xxxxxx[k] = -g_z_0_yyyy_xxxxxx[k] * ab_y + g_z_0_yyyy_xxxxxxy[k];

                g_z_0_yyyyy_xxxxxy[k] = -g_z_0_yyyy_xxxxxy[k] * ab_y + g_z_0_yyyy_xxxxxyy[k];

                g_z_0_yyyyy_xxxxxz[k] = -g_z_0_yyyy_xxxxxz[k] * ab_y + g_z_0_yyyy_xxxxxyz[k];

                g_z_0_yyyyy_xxxxyy[k] = -g_z_0_yyyy_xxxxyy[k] * ab_y + g_z_0_yyyy_xxxxyyy[k];

                g_z_0_yyyyy_xxxxyz[k] = -g_z_0_yyyy_xxxxyz[k] * ab_y + g_z_0_yyyy_xxxxyyz[k];

                g_z_0_yyyyy_xxxxzz[k] = -g_z_0_yyyy_xxxxzz[k] * ab_y + g_z_0_yyyy_xxxxyzz[k];

                g_z_0_yyyyy_xxxyyy[k] = -g_z_0_yyyy_xxxyyy[k] * ab_y + g_z_0_yyyy_xxxyyyy[k];

                g_z_0_yyyyy_xxxyyz[k] = -g_z_0_yyyy_xxxyyz[k] * ab_y + g_z_0_yyyy_xxxyyyz[k];

                g_z_0_yyyyy_xxxyzz[k] = -g_z_0_yyyy_xxxyzz[k] * ab_y + g_z_0_yyyy_xxxyyzz[k];

                g_z_0_yyyyy_xxxzzz[k] = -g_z_0_yyyy_xxxzzz[k] * ab_y + g_z_0_yyyy_xxxyzzz[k];

                g_z_0_yyyyy_xxyyyy[k] = -g_z_0_yyyy_xxyyyy[k] * ab_y + g_z_0_yyyy_xxyyyyy[k];

                g_z_0_yyyyy_xxyyyz[k] = -g_z_0_yyyy_xxyyyz[k] * ab_y + g_z_0_yyyy_xxyyyyz[k];

                g_z_0_yyyyy_xxyyzz[k] = -g_z_0_yyyy_xxyyzz[k] * ab_y + g_z_0_yyyy_xxyyyzz[k];

                g_z_0_yyyyy_xxyzzz[k] = -g_z_0_yyyy_xxyzzz[k] * ab_y + g_z_0_yyyy_xxyyzzz[k];

                g_z_0_yyyyy_xxzzzz[k] = -g_z_0_yyyy_xxzzzz[k] * ab_y + g_z_0_yyyy_xxyzzzz[k];

                g_z_0_yyyyy_xyyyyy[k] = -g_z_0_yyyy_xyyyyy[k] * ab_y + g_z_0_yyyy_xyyyyyy[k];

                g_z_0_yyyyy_xyyyyz[k] = -g_z_0_yyyy_xyyyyz[k] * ab_y + g_z_0_yyyy_xyyyyyz[k];

                g_z_0_yyyyy_xyyyzz[k] = -g_z_0_yyyy_xyyyzz[k] * ab_y + g_z_0_yyyy_xyyyyzz[k];

                g_z_0_yyyyy_xyyzzz[k] = -g_z_0_yyyy_xyyzzz[k] * ab_y + g_z_0_yyyy_xyyyzzz[k];

                g_z_0_yyyyy_xyzzzz[k] = -g_z_0_yyyy_xyzzzz[k] * ab_y + g_z_0_yyyy_xyyzzzz[k];

                g_z_0_yyyyy_xzzzzz[k] = -g_z_0_yyyy_xzzzzz[k] * ab_y + g_z_0_yyyy_xyzzzzz[k];

                g_z_0_yyyyy_yyyyyy[k] = -g_z_0_yyyy_yyyyyy[k] * ab_y + g_z_0_yyyy_yyyyyyy[k];

                g_z_0_yyyyy_yyyyyz[k] = -g_z_0_yyyy_yyyyyz[k] * ab_y + g_z_0_yyyy_yyyyyyz[k];

                g_z_0_yyyyy_yyyyzz[k] = -g_z_0_yyyy_yyyyzz[k] * ab_y + g_z_0_yyyy_yyyyyzz[k];

                g_z_0_yyyyy_yyyzzz[k] = -g_z_0_yyyy_yyyzzz[k] * ab_y + g_z_0_yyyy_yyyyzzz[k];

                g_z_0_yyyyy_yyzzzz[k] = -g_z_0_yyyy_yyzzzz[k] * ab_y + g_z_0_yyyy_yyyzzzz[k];

                g_z_0_yyyyy_yzzzzz[k] = -g_z_0_yyyy_yzzzzz[k] * ab_y + g_z_0_yyyy_yyzzzzz[k];

                g_z_0_yyyyy_zzzzzz[k] = -g_z_0_yyyy_zzzzzz[k] * ab_y + g_z_0_yyyy_yzzzzzz[k];
            }

            /// Set up 1624-1652 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyz_xxxxxx = cbuffer.data(hi_geom_10_off + 1624 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxy = cbuffer.data(hi_geom_10_off + 1625 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxz = cbuffer.data(hi_geom_10_off + 1626 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxyy = cbuffer.data(hi_geom_10_off + 1627 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxyz = cbuffer.data(hi_geom_10_off + 1628 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxzz = cbuffer.data(hi_geom_10_off + 1629 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyyy = cbuffer.data(hi_geom_10_off + 1630 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyyz = cbuffer.data(hi_geom_10_off + 1631 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyzz = cbuffer.data(hi_geom_10_off + 1632 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxzzz = cbuffer.data(hi_geom_10_off + 1633 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyyy = cbuffer.data(hi_geom_10_off + 1634 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyyz = cbuffer.data(hi_geom_10_off + 1635 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyzz = cbuffer.data(hi_geom_10_off + 1636 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyzzz = cbuffer.data(hi_geom_10_off + 1637 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxzzzz = cbuffer.data(hi_geom_10_off + 1638 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyyy = cbuffer.data(hi_geom_10_off + 1639 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyyz = cbuffer.data(hi_geom_10_off + 1640 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyzz = cbuffer.data(hi_geom_10_off + 1641 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyzzz = cbuffer.data(hi_geom_10_off + 1642 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyzzzz = cbuffer.data(hi_geom_10_off + 1643 * ccomps * dcomps);

            auto g_z_0_yyyyz_xzzzzz = cbuffer.data(hi_geom_10_off + 1644 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyyy = cbuffer.data(hi_geom_10_off + 1645 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyyz = cbuffer.data(hi_geom_10_off + 1646 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyzz = cbuffer.data(hi_geom_10_off + 1647 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyzzz = cbuffer.data(hi_geom_10_off + 1648 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyzzzz = cbuffer.data(hi_geom_10_off + 1649 * ccomps * dcomps);

            auto g_z_0_yyyyz_yzzzzz = cbuffer.data(hi_geom_10_off + 1650 * ccomps * dcomps);

            auto g_z_0_yyyyz_zzzzzz = cbuffer.data(hi_geom_10_off + 1651 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyz_xxxxxx, g_z_0_yyyyz_xxxxxy, g_z_0_yyyyz_xxxxxz, g_z_0_yyyyz_xxxxyy, g_z_0_yyyyz_xxxxyz, g_z_0_yyyyz_xxxxzz, g_z_0_yyyyz_xxxyyy, g_z_0_yyyyz_xxxyyz, g_z_0_yyyyz_xxxyzz, g_z_0_yyyyz_xxxzzz, g_z_0_yyyyz_xxyyyy, g_z_0_yyyyz_xxyyyz, g_z_0_yyyyz_xxyyzz, g_z_0_yyyyz_xxyzzz, g_z_0_yyyyz_xxzzzz, g_z_0_yyyyz_xyyyyy, g_z_0_yyyyz_xyyyyz, g_z_0_yyyyz_xyyyzz, g_z_0_yyyyz_xyyzzz, g_z_0_yyyyz_xyzzzz, g_z_0_yyyyz_xzzzzz, g_z_0_yyyyz_yyyyyy, g_z_0_yyyyz_yyyyyz, g_z_0_yyyyz_yyyyzz, g_z_0_yyyyz_yyyzzz, g_z_0_yyyyz_yyzzzz, g_z_0_yyyyz_yzzzzz, g_z_0_yyyyz_zzzzzz, g_z_0_yyyz_xxxxxx, g_z_0_yyyz_xxxxxxy, g_z_0_yyyz_xxxxxy, g_z_0_yyyz_xxxxxyy, g_z_0_yyyz_xxxxxyz, g_z_0_yyyz_xxxxxz, g_z_0_yyyz_xxxxyy, g_z_0_yyyz_xxxxyyy, g_z_0_yyyz_xxxxyyz, g_z_0_yyyz_xxxxyz, g_z_0_yyyz_xxxxyzz, g_z_0_yyyz_xxxxzz, g_z_0_yyyz_xxxyyy, g_z_0_yyyz_xxxyyyy, g_z_0_yyyz_xxxyyyz, g_z_0_yyyz_xxxyyz, g_z_0_yyyz_xxxyyzz, g_z_0_yyyz_xxxyzz, g_z_0_yyyz_xxxyzzz, g_z_0_yyyz_xxxzzz, g_z_0_yyyz_xxyyyy, g_z_0_yyyz_xxyyyyy, g_z_0_yyyz_xxyyyyz, g_z_0_yyyz_xxyyyz, g_z_0_yyyz_xxyyyzz, g_z_0_yyyz_xxyyzz, g_z_0_yyyz_xxyyzzz, g_z_0_yyyz_xxyzzz, g_z_0_yyyz_xxyzzzz, g_z_0_yyyz_xxzzzz, g_z_0_yyyz_xyyyyy, g_z_0_yyyz_xyyyyyy, g_z_0_yyyz_xyyyyyz, g_z_0_yyyz_xyyyyz, g_z_0_yyyz_xyyyyzz, g_z_0_yyyz_xyyyzz, g_z_0_yyyz_xyyyzzz, g_z_0_yyyz_xyyzzz, g_z_0_yyyz_xyyzzzz, g_z_0_yyyz_xyzzzz, g_z_0_yyyz_xyzzzzz, g_z_0_yyyz_xzzzzz, g_z_0_yyyz_yyyyyy, g_z_0_yyyz_yyyyyyy, g_z_0_yyyz_yyyyyyz, g_z_0_yyyz_yyyyyz, g_z_0_yyyz_yyyyyzz, g_z_0_yyyz_yyyyzz, g_z_0_yyyz_yyyyzzz, g_z_0_yyyz_yyyzzz, g_z_0_yyyz_yyyzzzz, g_z_0_yyyz_yyzzzz, g_z_0_yyyz_yyzzzzz, g_z_0_yyyz_yzzzzz, g_z_0_yyyz_yzzzzzz, g_z_0_yyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyz_xxxxxx[k] = -g_z_0_yyyz_xxxxxx[k] * ab_y + g_z_0_yyyz_xxxxxxy[k];

                g_z_0_yyyyz_xxxxxy[k] = -g_z_0_yyyz_xxxxxy[k] * ab_y + g_z_0_yyyz_xxxxxyy[k];

                g_z_0_yyyyz_xxxxxz[k] = -g_z_0_yyyz_xxxxxz[k] * ab_y + g_z_0_yyyz_xxxxxyz[k];

                g_z_0_yyyyz_xxxxyy[k] = -g_z_0_yyyz_xxxxyy[k] * ab_y + g_z_0_yyyz_xxxxyyy[k];

                g_z_0_yyyyz_xxxxyz[k] = -g_z_0_yyyz_xxxxyz[k] * ab_y + g_z_0_yyyz_xxxxyyz[k];

                g_z_0_yyyyz_xxxxzz[k] = -g_z_0_yyyz_xxxxzz[k] * ab_y + g_z_0_yyyz_xxxxyzz[k];

                g_z_0_yyyyz_xxxyyy[k] = -g_z_0_yyyz_xxxyyy[k] * ab_y + g_z_0_yyyz_xxxyyyy[k];

                g_z_0_yyyyz_xxxyyz[k] = -g_z_0_yyyz_xxxyyz[k] * ab_y + g_z_0_yyyz_xxxyyyz[k];

                g_z_0_yyyyz_xxxyzz[k] = -g_z_0_yyyz_xxxyzz[k] * ab_y + g_z_0_yyyz_xxxyyzz[k];

                g_z_0_yyyyz_xxxzzz[k] = -g_z_0_yyyz_xxxzzz[k] * ab_y + g_z_0_yyyz_xxxyzzz[k];

                g_z_0_yyyyz_xxyyyy[k] = -g_z_0_yyyz_xxyyyy[k] * ab_y + g_z_0_yyyz_xxyyyyy[k];

                g_z_0_yyyyz_xxyyyz[k] = -g_z_0_yyyz_xxyyyz[k] * ab_y + g_z_0_yyyz_xxyyyyz[k];

                g_z_0_yyyyz_xxyyzz[k] = -g_z_0_yyyz_xxyyzz[k] * ab_y + g_z_0_yyyz_xxyyyzz[k];

                g_z_0_yyyyz_xxyzzz[k] = -g_z_0_yyyz_xxyzzz[k] * ab_y + g_z_0_yyyz_xxyyzzz[k];

                g_z_0_yyyyz_xxzzzz[k] = -g_z_0_yyyz_xxzzzz[k] * ab_y + g_z_0_yyyz_xxyzzzz[k];

                g_z_0_yyyyz_xyyyyy[k] = -g_z_0_yyyz_xyyyyy[k] * ab_y + g_z_0_yyyz_xyyyyyy[k];

                g_z_0_yyyyz_xyyyyz[k] = -g_z_0_yyyz_xyyyyz[k] * ab_y + g_z_0_yyyz_xyyyyyz[k];

                g_z_0_yyyyz_xyyyzz[k] = -g_z_0_yyyz_xyyyzz[k] * ab_y + g_z_0_yyyz_xyyyyzz[k];

                g_z_0_yyyyz_xyyzzz[k] = -g_z_0_yyyz_xyyzzz[k] * ab_y + g_z_0_yyyz_xyyyzzz[k];

                g_z_0_yyyyz_xyzzzz[k] = -g_z_0_yyyz_xyzzzz[k] * ab_y + g_z_0_yyyz_xyyzzzz[k];

                g_z_0_yyyyz_xzzzzz[k] = -g_z_0_yyyz_xzzzzz[k] * ab_y + g_z_0_yyyz_xyzzzzz[k];

                g_z_0_yyyyz_yyyyyy[k] = -g_z_0_yyyz_yyyyyy[k] * ab_y + g_z_0_yyyz_yyyyyyy[k];

                g_z_0_yyyyz_yyyyyz[k] = -g_z_0_yyyz_yyyyyz[k] * ab_y + g_z_0_yyyz_yyyyyyz[k];

                g_z_0_yyyyz_yyyyzz[k] = -g_z_0_yyyz_yyyyzz[k] * ab_y + g_z_0_yyyz_yyyyyzz[k];

                g_z_0_yyyyz_yyyzzz[k] = -g_z_0_yyyz_yyyzzz[k] * ab_y + g_z_0_yyyz_yyyyzzz[k];

                g_z_0_yyyyz_yyzzzz[k] = -g_z_0_yyyz_yyzzzz[k] * ab_y + g_z_0_yyyz_yyyzzzz[k];

                g_z_0_yyyyz_yzzzzz[k] = -g_z_0_yyyz_yzzzzz[k] * ab_y + g_z_0_yyyz_yyzzzzz[k];

                g_z_0_yyyyz_zzzzzz[k] = -g_z_0_yyyz_zzzzzz[k] * ab_y + g_z_0_yyyz_yzzzzzz[k];
            }

            /// Set up 1652-1680 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1652 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1653 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1654 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1655 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1656 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1657 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1658 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1659 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1660 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1661 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1662 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1663 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1664 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1665 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1666 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1667 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1668 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1669 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1670 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1671 * ccomps * dcomps);

            auto g_z_0_yyyzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1672 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1673 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1674 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1675 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1676 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1677 * ccomps * dcomps);

            auto g_z_0_yyyzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1678 * ccomps * dcomps);

            auto g_z_0_yyyzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1679 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyzz_xxxxxx, g_z_0_yyyzz_xxxxxy, g_z_0_yyyzz_xxxxxz, g_z_0_yyyzz_xxxxyy, g_z_0_yyyzz_xxxxyz, g_z_0_yyyzz_xxxxzz, g_z_0_yyyzz_xxxyyy, g_z_0_yyyzz_xxxyyz, g_z_0_yyyzz_xxxyzz, g_z_0_yyyzz_xxxzzz, g_z_0_yyyzz_xxyyyy, g_z_0_yyyzz_xxyyyz, g_z_0_yyyzz_xxyyzz, g_z_0_yyyzz_xxyzzz, g_z_0_yyyzz_xxzzzz, g_z_0_yyyzz_xyyyyy, g_z_0_yyyzz_xyyyyz, g_z_0_yyyzz_xyyyzz, g_z_0_yyyzz_xyyzzz, g_z_0_yyyzz_xyzzzz, g_z_0_yyyzz_xzzzzz, g_z_0_yyyzz_yyyyyy, g_z_0_yyyzz_yyyyyz, g_z_0_yyyzz_yyyyzz, g_z_0_yyyzz_yyyzzz, g_z_0_yyyzz_yyzzzz, g_z_0_yyyzz_yzzzzz, g_z_0_yyyzz_zzzzzz, g_z_0_yyzz_xxxxxx, g_z_0_yyzz_xxxxxxy, g_z_0_yyzz_xxxxxy, g_z_0_yyzz_xxxxxyy, g_z_0_yyzz_xxxxxyz, g_z_0_yyzz_xxxxxz, g_z_0_yyzz_xxxxyy, g_z_0_yyzz_xxxxyyy, g_z_0_yyzz_xxxxyyz, g_z_0_yyzz_xxxxyz, g_z_0_yyzz_xxxxyzz, g_z_0_yyzz_xxxxzz, g_z_0_yyzz_xxxyyy, g_z_0_yyzz_xxxyyyy, g_z_0_yyzz_xxxyyyz, g_z_0_yyzz_xxxyyz, g_z_0_yyzz_xxxyyzz, g_z_0_yyzz_xxxyzz, g_z_0_yyzz_xxxyzzz, g_z_0_yyzz_xxxzzz, g_z_0_yyzz_xxyyyy, g_z_0_yyzz_xxyyyyy, g_z_0_yyzz_xxyyyyz, g_z_0_yyzz_xxyyyz, g_z_0_yyzz_xxyyyzz, g_z_0_yyzz_xxyyzz, g_z_0_yyzz_xxyyzzz, g_z_0_yyzz_xxyzzz, g_z_0_yyzz_xxyzzzz, g_z_0_yyzz_xxzzzz, g_z_0_yyzz_xyyyyy, g_z_0_yyzz_xyyyyyy, g_z_0_yyzz_xyyyyyz, g_z_0_yyzz_xyyyyz, g_z_0_yyzz_xyyyyzz, g_z_0_yyzz_xyyyzz, g_z_0_yyzz_xyyyzzz, g_z_0_yyzz_xyyzzz, g_z_0_yyzz_xyyzzzz, g_z_0_yyzz_xyzzzz, g_z_0_yyzz_xyzzzzz, g_z_0_yyzz_xzzzzz, g_z_0_yyzz_yyyyyy, g_z_0_yyzz_yyyyyyy, g_z_0_yyzz_yyyyyyz, g_z_0_yyzz_yyyyyz, g_z_0_yyzz_yyyyyzz, g_z_0_yyzz_yyyyzz, g_z_0_yyzz_yyyyzzz, g_z_0_yyzz_yyyzzz, g_z_0_yyzz_yyyzzzz, g_z_0_yyzz_yyzzzz, g_z_0_yyzz_yyzzzzz, g_z_0_yyzz_yzzzzz, g_z_0_yyzz_yzzzzzz, g_z_0_yyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzz_xxxxxx[k] = -g_z_0_yyzz_xxxxxx[k] * ab_y + g_z_0_yyzz_xxxxxxy[k];

                g_z_0_yyyzz_xxxxxy[k] = -g_z_0_yyzz_xxxxxy[k] * ab_y + g_z_0_yyzz_xxxxxyy[k];

                g_z_0_yyyzz_xxxxxz[k] = -g_z_0_yyzz_xxxxxz[k] * ab_y + g_z_0_yyzz_xxxxxyz[k];

                g_z_0_yyyzz_xxxxyy[k] = -g_z_0_yyzz_xxxxyy[k] * ab_y + g_z_0_yyzz_xxxxyyy[k];

                g_z_0_yyyzz_xxxxyz[k] = -g_z_0_yyzz_xxxxyz[k] * ab_y + g_z_0_yyzz_xxxxyyz[k];

                g_z_0_yyyzz_xxxxzz[k] = -g_z_0_yyzz_xxxxzz[k] * ab_y + g_z_0_yyzz_xxxxyzz[k];

                g_z_0_yyyzz_xxxyyy[k] = -g_z_0_yyzz_xxxyyy[k] * ab_y + g_z_0_yyzz_xxxyyyy[k];

                g_z_0_yyyzz_xxxyyz[k] = -g_z_0_yyzz_xxxyyz[k] * ab_y + g_z_0_yyzz_xxxyyyz[k];

                g_z_0_yyyzz_xxxyzz[k] = -g_z_0_yyzz_xxxyzz[k] * ab_y + g_z_0_yyzz_xxxyyzz[k];

                g_z_0_yyyzz_xxxzzz[k] = -g_z_0_yyzz_xxxzzz[k] * ab_y + g_z_0_yyzz_xxxyzzz[k];

                g_z_0_yyyzz_xxyyyy[k] = -g_z_0_yyzz_xxyyyy[k] * ab_y + g_z_0_yyzz_xxyyyyy[k];

                g_z_0_yyyzz_xxyyyz[k] = -g_z_0_yyzz_xxyyyz[k] * ab_y + g_z_0_yyzz_xxyyyyz[k];

                g_z_0_yyyzz_xxyyzz[k] = -g_z_0_yyzz_xxyyzz[k] * ab_y + g_z_0_yyzz_xxyyyzz[k];

                g_z_0_yyyzz_xxyzzz[k] = -g_z_0_yyzz_xxyzzz[k] * ab_y + g_z_0_yyzz_xxyyzzz[k];

                g_z_0_yyyzz_xxzzzz[k] = -g_z_0_yyzz_xxzzzz[k] * ab_y + g_z_0_yyzz_xxyzzzz[k];

                g_z_0_yyyzz_xyyyyy[k] = -g_z_0_yyzz_xyyyyy[k] * ab_y + g_z_0_yyzz_xyyyyyy[k];

                g_z_0_yyyzz_xyyyyz[k] = -g_z_0_yyzz_xyyyyz[k] * ab_y + g_z_0_yyzz_xyyyyyz[k];

                g_z_0_yyyzz_xyyyzz[k] = -g_z_0_yyzz_xyyyzz[k] * ab_y + g_z_0_yyzz_xyyyyzz[k];

                g_z_0_yyyzz_xyyzzz[k] = -g_z_0_yyzz_xyyzzz[k] * ab_y + g_z_0_yyzz_xyyyzzz[k];

                g_z_0_yyyzz_xyzzzz[k] = -g_z_0_yyzz_xyzzzz[k] * ab_y + g_z_0_yyzz_xyyzzzz[k];

                g_z_0_yyyzz_xzzzzz[k] = -g_z_0_yyzz_xzzzzz[k] * ab_y + g_z_0_yyzz_xyzzzzz[k];

                g_z_0_yyyzz_yyyyyy[k] = -g_z_0_yyzz_yyyyyy[k] * ab_y + g_z_0_yyzz_yyyyyyy[k];

                g_z_0_yyyzz_yyyyyz[k] = -g_z_0_yyzz_yyyyyz[k] * ab_y + g_z_0_yyzz_yyyyyyz[k];

                g_z_0_yyyzz_yyyyzz[k] = -g_z_0_yyzz_yyyyzz[k] * ab_y + g_z_0_yyzz_yyyyyzz[k];

                g_z_0_yyyzz_yyyzzz[k] = -g_z_0_yyzz_yyyzzz[k] * ab_y + g_z_0_yyzz_yyyyzzz[k];

                g_z_0_yyyzz_yyzzzz[k] = -g_z_0_yyzz_yyzzzz[k] * ab_y + g_z_0_yyzz_yyyzzzz[k];

                g_z_0_yyyzz_yzzzzz[k] = -g_z_0_yyzz_yzzzzz[k] * ab_y + g_z_0_yyzz_yyzzzzz[k];

                g_z_0_yyyzz_zzzzzz[k] = -g_z_0_yyzz_zzzzzz[k] * ab_y + g_z_0_yyzz_yzzzzzz[k];
            }

            /// Set up 1680-1708 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1680 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1681 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1682 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1683 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1684 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1685 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1686 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1687 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1688 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1689 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1690 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1691 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1692 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1693 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1694 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1695 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1696 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1697 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1698 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1699 * ccomps * dcomps);

            auto g_z_0_yyzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1700 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1701 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1702 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1703 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1704 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1705 * ccomps * dcomps);

            auto g_z_0_yyzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1706 * ccomps * dcomps);

            auto g_z_0_yyzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1707 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzzz_xxxxxx, g_z_0_yyzzz_xxxxxy, g_z_0_yyzzz_xxxxxz, g_z_0_yyzzz_xxxxyy, g_z_0_yyzzz_xxxxyz, g_z_0_yyzzz_xxxxzz, g_z_0_yyzzz_xxxyyy, g_z_0_yyzzz_xxxyyz, g_z_0_yyzzz_xxxyzz, g_z_0_yyzzz_xxxzzz, g_z_0_yyzzz_xxyyyy, g_z_0_yyzzz_xxyyyz, g_z_0_yyzzz_xxyyzz, g_z_0_yyzzz_xxyzzz, g_z_0_yyzzz_xxzzzz, g_z_0_yyzzz_xyyyyy, g_z_0_yyzzz_xyyyyz, g_z_0_yyzzz_xyyyzz, g_z_0_yyzzz_xyyzzz, g_z_0_yyzzz_xyzzzz, g_z_0_yyzzz_xzzzzz, g_z_0_yyzzz_yyyyyy, g_z_0_yyzzz_yyyyyz, g_z_0_yyzzz_yyyyzz, g_z_0_yyzzz_yyyzzz, g_z_0_yyzzz_yyzzzz, g_z_0_yyzzz_yzzzzz, g_z_0_yyzzz_zzzzzz, g_z_0_yzzz_xxxxxx, g_z_0_yzzz_xxxxxxy, g_z_0_yzzz_xxxxxy, g_z_0_yzzz_xxxxxyy, g_z_0_yzzz_xxxxxyz, g_z_0_yzzz_xxxxxz, g_z_0_yzzz_xxxxyy, g_z_0_yzzz_xxxxyyy, g_z_0_yzzz_xxxxyyz, g_z_0_yzzz_xxxxyz, g_z_0_yzzz_xxxxyzz, g_z_0_yzzz_xxxxzz, g_z_0_yzzz_xxxyyy, g_z_0_yzzz_xxxyyyy, g_z_0_yzzz_xxxyyyz, g_z_0_yzzz_xxxyyz, g_z_0_yzzz_xxxyyzz, g_z_0_yzzz_xxxyzz, g_z_0_yzzz_xxxyzzz, g_z_0_yzzz_xxxzzz, g_z_0_yzzz_xxyyyy, g_z_0_yzzz_xxyyyyy, g_z_0_yzzz_xxyyyyz, g_z_0_yzzz_xxyyyz, g_z_0_yzzz_xxyyyzz, g_z_0_yzzz_xxyyzz, g_z_0_yzzz_xxyyzzz, g_z_0_yzzz_xxyzzz, g_z_0_yzzz_xxyzzzz, g_z_0_yzzz_xxzzzz, g_z_0_yzzz_xyyyyy, g_z_0_yzzz_xyyyyyy, g_z_0_yzzz_xyyyyyz, g_z_0_yzzz_xyyyyz, g_z_0_yzzz_xyyyyzz, g_z_0_yzzz_xyyyzz, g_z_0_yzzz_xyyyzzz, g_z_0_yzzz_xyyzzz, g_z_0_yzzz_xyyzzzz, g_z_0_yzzz_xyzzzz, g_z_0_yzzz_xyzzzzz, g_z_0_yzzz_xzzzzz, g_z_0_yzzz_yyyyyy, g_z_0_yzzz_yyyyyyy, g_z_0_yzzz_yyyyyyz, g_z_0_yzzz_yyyyyz, g_z_0_yzzz_yyyyyzz, g_z_0_yzzz_yyyyzz, g_z_0_yzzz_yyyyzzz, g_z_0_yzzz_yyyzzz, g_z_0_yzzz_yyyzzzz, g_z_0_yzzz_yyzzzz, g_z_0_yzzz_yyzzzzz, g_z_0_yzzz_yzzzzz, g_z_0_yzzz_yzzzzzz, g_z_0_yzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzz_xxxxxx[k] = -g_z_0_yzzz_xxxxxx[k] * ab_y + g_z_0_yzzz_xxxxxxy[k];

                g_z_0_yyzzz_xxxxxy[k] = -g_z_0_yzzz_xxxxxy[k] * ab_y + g_z_0_yzzz_xxxxxyy[k];

                g_z_0_yyzzz_xxxxxz[k] = -g_z_0_yzzz_xxxxxz[k] * ab_y + g_z_0_yzzz_xxxxxyz[k];

                g_z_0_yyzzz_xxxxyy[k] = -g_z_0_yzzz_xxxxyy[k] * ab_y + g_z_0_yzzz_xxxxyyy[k];

                g_z_0_yyzzz_xxxxyz[k] = -g_z_0_yzzz_xxxxyz[k] * ab_y + g_z_0_yzzz_xxxxyyz[k];

                g_z_0_yyzzz_xxxxzz[k] = -g_z_0_yzzz_xxxxzz[k] * ab_y + g_z_0_yzzz_xxxxyzz[k];

                g_z_0_yyzzz_xxxyyy[k] = -g_z_0_yzzz_xxxyyy[k] * ab_y + g_z_0_yzzz_xxxyyyy[k];

                g_z_0_yyzzz_xxxyyz[k] = -g_z_0_yzzz_xxxyyz[k] * ab_y + g_z_0_yzzz_xxxyyyz[k];

                g_z_0_yyzzz_xxxyzz[k] = -g_z_0_yzzz_xxxyzz[k] * ab_y + g_z_0_yzzz_xxxyyzz[k];

                g_z_0_yyzzz_xxxzzz[k] = -g_z_0_yzzz_xxxzzz[k] * ab_y + g_z_0_yzzz_xxxyzzz[k];

                g_z_0_yyzzz_xxyyyy[k] = -g_z_0_yzzz_xxyyyy[k] * ab_y + g_z_0_yzzz_xxyyyyy[k];

                g_z_0_yyzzz_xxyyyz[k] = -g_z_0_yzzz_xxyyyz[k] * ab_y + g_z_0_yzzz_xxyyyyz[k];

                g_z_0_yyzzz_xxyyzz[k] = -g_z_0_yzzz_xxyyzz[k] * ab_y + g_z_0_yzzz_xxyyyzz[k];

                g_z_0_yyzzz_xxyzzz[k] = -g_z_0_yzzz_xxyzzz[k] * ab_y + g_z_0_yzzz_xxyyzzz[k];

                g_z_0_yyzzz_xxzzzz[k] = -g_z_0_yzzz_xxzzzz[k] * ab_y + g_z_0_yzzz_xxyzzzz[k];

                g_z_0_yyzzz_xyyyyy[k] = -g_z_0_yzzz_xyyyyy[k] * ab_y + g_z_0_yzzz_xyyyyyy[k];

                g_z_0_yyzzz_xyyyyz[k] = -g_z_0_yzzz_xyyyyz[k] * ab_y + g_z_0_yzzz_xyyyyyz[k];

                g_z_0_yyzzz_xyyyzz[k] = -g_z_0_yzzz_xyyyzz[k] * ab_y + g_z_0_yzzz_xyyyyzz[k];

                g_z_0_yyzzz_xyyzzz[k] = -g_z_0_yzzz_xyyzzz[k] * ab_y + g_z_0_yzzz_xyyyzzz[k];

                g_z_0_yyzzz_xyzzzz[k] = -g_z_0_yzzz_xyzzzz[k] * ab_y + g_z_0_yzzz_xyyzzzz[k];

                g_z_0_yyzzz_xzzzzz[k] = -g_z_0_yzzz_xzzzzz[k] * ab_y + g_z_0_yzzz_xyzzzzz[k];

                g_z_0_yyzzz_yyyyyy[k] = -g_z_0_yzzz_yyyyyy[k] * ab_y + g_z_0_yzzz_yyyyyyy[k];

                g_z_0_yyzzz_yyyyyz[k] = -g_z_0_yzzz_yyyyyz[k] * ab_y + g_z_0_yzzz_yyyyyyz[k];

                g_z_0_yyzzz_yyyyzz[k] = -g_z_0_yzzz_yyyyzz[k] * ab_y + g_z_0_yzzz_yyyyyzz[k];

                g_z_0_yyzzz_yyyzzz[k] = -g_z_0_yzzz_yyyzzz[k] * ab_y + g_z_0_yzzz_yyyyzzz[k];

                g_z_0_yyzzz_yyzzzz[k] = -g_z_0_yzzz_yyzzzz[k] * ab_y + g_z_0_yzzz_yyyzzzz[k];

                g_z_0_yyzzz_yzzzzz[k] = -g_z_0_yzzz_yzzzzz[k] * ab_y + g_z_0_yzzz_yyzzzzz[k];

                g_z_0_yyzzz_zzzzzz[k] = -g_z_0_yzzz_zzzzzz[k] * ab_y + g_z_0_yzzz_yzzzzzz[k];
            }

            /// Set up 1708-1736 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1708 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1709 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1710 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1711 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1712 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1713 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1714 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1715 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1716 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1717 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1718 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1719 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1720 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1721 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1722 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1723 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1724 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1725 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1726 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1727 * ccomps * dcomps);

            auto g_z_0_yzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1728 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1729 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1730 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1731 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1732 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1733 * ccomps * dcomps);

            auto g_z_0_yzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1734 * ccomps * dcomps);

            auto g_z_0_yzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1735 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzzz_xxxxxx, g_z_0_yzzzz_xxxxxy, g_z_0_yzzzz_xxxxxz, g_z_0_yzzzz_xxxxyy, g_z_0_yzzzz_xxxxyz, g_z_0_yzzzz_xxxxzz, g_z_0_yzzzz_xxxyyy, g_z_0_yzzzz_xxxyyz, g_z_0_yzzzz_xxxyzz, g_z_0_yzzzz_xxxzzz, g_z_0_yzzzz_xxyyyy, g_z_0_yzzzz_xxyyyz, g_z_0_yzzzz_xxyyzz, g_z_0_yzzzz_xxyzzz, g_z_0_yzzzz_xxzzzz, g_z_0_yzzzz_xyyyyy, g_z_0_yzzzz_xyyyyz, g_z_0_yzzzz_xyyyzz, g_z_0_yzzzz_xyyzzz, g_z_0_yzzzz_xyzzzz, g_z_0_yzzzz_xzzzzz, g_z_0_yzzzz_yyyyyy, g_z_0_yzzzz_yyyyyz, g_z_0_yzzzz_yyyyzz, g_z_0_yzzzz_yyyzzz, g_z_0_yzzzz_yyzzzz, g_z_0_yzzzz_yzzzzz, g_z_0_yzzzz_zzzzzz, g_z_0_zzzz_xxxxxx, g_z_0_zzzz_xxxxxxy, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxxyy, g_z_0_zzzz_xxxxxyz, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyyy, g_z_0_zzzz_xxxxyyz, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxyzz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyyy, g_z_0_zzzz_xxxyyyz, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyyzz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxyzzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyyy, g_z_0_zzzz_xxyyyyz, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyyzz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyyzzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxyzzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyyy, g_z_0_zzzz_xyyyyyz, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyyzz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyyzzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyyzzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xyzzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_yyyyyy, g_z_0_zzzz_yyyyyyy, g_z_0_zzzz_yyyyyyz, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyyzz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyyzzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyyzzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yyzzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_yzzzzzz, g_z_0_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzz_xxxxxx[k] = -g_z_0_zzzz_xxxxxx[k] * ab_y + g_z_0_zzzz_xxxxxxy[k];

                g_z_0_yzzzz_xxxxxy[k] = -g_z_0_zzzz_xxxxxy[k] * ab_y + g_z_0_zzzz_xxxxxyy[k];

                g_z_0_yzzzz_xxxxxz[k] = -g_z_0_zzzz_xxxxxz[k] * ab_y + g_z_0_zzzz_xxxxxyz[k];

                g_z_0_yzzzz_xxxxyy[k] = -g_z_0_zzzz_xxxxyy[k] * ab_y + g_z_0_zzzz_xxxxyyy[k];

                g_z_0_yzzzz_xxxxyz[k] = -g_z_0_zzzz_xxxxyz[k] * ab_y + g_z_0_zzzz_xxxxyyz[k];

                g_z_0_yzzzz_xxxxzz[k] = -g_z_0_zzzz_xxxxzz[k] * ab_y + g_z_0_zzzz_xxxxyzz[k];

                g_z_0_yzzzz_xxxyyy[k] = -g_z_0_zzzz_xxxyyy[k] * ab_y + g_z_0_zzzz_xxxyyyy[k];

                g_z_0_yzzzz_xxxyyz[k] = -g_z_0_zzzz_xxxyyz[k] * ab_y + g_z_0_zzzz_xxxyyyz[k];

                g_z_0_yzzzz_xxxyzz[k] = -g_z_0_zzzz_xxxyzz[k] * ab_y + g_z_0_zzzz_xxxyyzz[k];

                g_z_0_yzzzz_xxxzzz[k] = -g_z_0_zzzz_xxxzzz[k] * ab_y + g_z_0_zzzz_xxxyzzz[k];

                g_z_0_yzzzz_xxyyyy[k] = -g_z_0_zzzz_xxyyyy[k] * ab_y + g_z_0_zzzz_xxyyyyy[k];

                g_z_0_yzzzz_xxyyyz[k] = -g_z_0_zzzz_xxyyyz[k] * ab_y + g_z_0_zzzz_xxyyyyz[k];

                g_z_0_yzzzz_xxyyzz[k] = -g_z_0_zzzz_xxyyzz[k] * ab_y + g_z_0_zzzz_xxyyyzz[k];

                g_z_0_yzzzz_xxyzzz[k] = -g_z_0_zzzz_xxyzzz[k] * ab_y + g_z_0_zzzz_xxyyzzz[k];

                g_z_0_yzzzz_xxzzzz[k] = -g_z_0_zzzz_xxzzzz[k] * ab_y + g_z_0_zzzz_xxyzzzz[k];

                g_z_0_yzzzz_xyyyyy[k] = -g_z_0_zzzz_xyyyyy[k] * ab_y + g_z_0_zzzz_xyyyyyy[k];

                g_z_0_yzzzz_xyyyyz[k] = -g_z_0_zzzz_xyyyyz[k] * ab_y + g_z_0_zzzz_xyyyyyz[k];

                g_z_0_yzzzz_xyyyzz[k] = -g_z_0_zzzz_xyyyzz[k] * ab_y + g_z_0_zzzz_xyyyyzz[k];

                g_z_0_yzzzz_xyyzzz[k] = -g_z_0_zzzz_xyyzzz[k] * ab_y + g_z_0_zzzz_xyyyzzz[k];

                g_z_0_yzzzz_xyzzzz[k] = -g_z_0_zzzz_xyzzzz[k] * ab_y + g_z_0_zzzz_xyyzzzz[k];

                g_z_0_yzzzz_xzzzzz[k] = -g_z_0_zzzz_xzzzzz[k] * ab_y + g_z_0_zzzz_xyzzzzz[k];

                g_z_0_yzzzz_yyyyyy[k] = -g_z_0_zzzz_yyyyyy[k] * ab_y + g_z_0_zzzz_yyyyyyy[k];

                g_z_0_yzzzz_yyyyyz[k] = -g_z_0_zzzz_yyyyyz[k] * ab_y + g_z_0_zzzz_yyyyyyz[k];

                g_z_0_yzzzz_yyyyzz[k] = -g_z_0_zzzz_yyyyzz[k] * ab_y + g_z_0_zzzz_yyyyyzz[k];

                g_z_0_yzzzz_yyyzzz[k] = -g_z_0_zzzz_yyyzzz[k] * ab_y + g_z_0_zzzz_yyyyzzz[k];

                g_z_0_yzzzz_yyzzzz[k] = -g_z_0_zzzz_yyzzzz[k] * ab_y + g_z_0_zzzz_yyyzzzz[k];

                g_z_0_yzzzz_yzzzzz[k] = -g_z_0_zzzz_yzzzzz[k] * ab_y + g_z_0_zzzz_yyzzzzz[k];

                g_z_0_yzzzz_zzzzzz[k] = -g_z_0_zzzz_zzzzzz[k] * ab_y + g_z_0_zzzz_yzzzzzz[k];
            }

            /// Set up 1736-1764 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzz_xxxxxx = cbuffer.data(hi_geom_10_off + 1736 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxy = cbuffer.data(hi_geom_10_off + 1737 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxz = cbuffer.data(hi_geom_10_off + 1738 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxyy = cbuffer.data(hi_geom_10_off + 1739 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxyz = cbuffer.data(hi_geom_10_off + 1740 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxzz = cbuffer.data(hi_geom_10_off + 1741 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyyy = cbuffer.data(hi_geom_10_off + 1742 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyyz = cbuffer.data(hi_geom_10_off + 1743 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyzz = cbuffer.data(hi_geom_10_off + 1744 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxzzz = cbuffer.data(hi_geom_10_off + 1745 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyyy = cbuffer.data(hi_geom_10_off + 1746 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyyz = cbuffer.data(hi_geom_10_off + 1747 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyzz = cbuffer.data(hi_geom_10_off + 1748 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyzzz = cbuffer.data(hi_geom_10_off + 1749 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxzzzz = cbuffer.data(hi_geom_10_off + 1750 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyyy = cbuffer.data(hi_geom_10_off + 1751 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyyz = cbuffer.data(hi_geom_10_off + 1752 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyzz = cbuffer.data(hi_geom_10_off + 1753 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyzzz = cbuffer.data(hi_geom_10_off + 1754 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyzzzz = cbuffer.data(hi_geom_10_off + 1755 * ccomps * dcomps);

            auto g_z_0_zzzzz_xzzzzz = cbuffer.data(hi_geom_10_off + 1756 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyyy = cbuffer.data(hi_geom_10_off + 1757 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyyz = cbuffer.data(hi_geom_10_off + 1758 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyzz = cbuffer.data(hi_geom_10_off + 1759 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyzzz = cbuffer.data(hi_geom_10_off + 1760 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyzzzz = cbuffer.data(hi_geom_10_off + 1761 * ccomps * dcomps);

            auto g_z_0_zzzzz_yzzzzz = cbuffer.data(hi_geom_10_off + 1762 * ccomps * dcomps);

            auto g_z_0_zzzzz_zzzzzz = cbuffer.data(hi_geom_10_off + 1763 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzz_xxxxxx, g_z_0_zzzz_xxxxxxz, g_z_0_zzzz_xxxxxy, g_z_0_zzzz_xxxxxyz, g_z_0_zzzz_xxxxxz, g_z_0_zzzz_xxxxxzz, g_z_0_zzzz_xxxxyy, g_z_0_zzzz_xxxxyyz, g_z_0_zzzz_xxxxyz, g_z_0_zzzz_xxxxyzz, g_z_0_zzzz_xxxxzz, g_z_0_zzzz_xxxxzzz, g_z_0_zzzz_xxxyyy, g_z_0_zzzz_xxxyyyz, g_z_0_zzzz_xxxyyz, g_z_0_zzzz_xxxyyzz, g_z_0_zzzz_xxxyzz, g_z_0_zzzz_xxxyzzz, g_z_0_zzzz_xxxzzz, g_z_0_zzzz_xxxzzzz, g_z_0_zzzz_xxyyyy, g_z_0_zzzz_xxyyyyz, g_z_0_zzzz_xxyyyz, g_z_0_zzzz_xxyyyzz, g_z_0_zzzz_xxyyzz, g_z_0_zzzz_xxyyzzz, g_z_0_zzzz_xxyzzz, g_z_0_zzzz_xxyzzzz, g_z_0_zzzz_xxzzzz, g_z_0_zzzz_xxzzzzz, g_z_0_zzzz_xyyyyy, g_z_0_zzzz_xyyyyyz, g_z_0_zzzz_xyyyyz, g_z_0_zzzz_xyyyyzz, g_z_0_zzzz_xyyyzz, g_z_0_zzzz_xyyyzzz, g_z_0_zzzz_xyyzzz, g_z_0_zzzz_xyyzzzz, g_z_0_zzzz_xyzzzz, g_z_0_zzzz_xyzzzzz, g_z_0_zzzz_xzzzzz, g_z_0_zzzz_xzzzzzz, g_z_0_zzzz_yyyyyy, g_z_0_zzzz_yyyyyyz, g_z_0_zzzz_yyyyyz, g_z_0_zzzz_yyyyyzz, g_z_0_zzzz_yyyyzz, g_z_0_zzzz_yyyyzzz, g_z_0_zzzz_yyyzzz, g_z_0_zzzz_yyyzzzz, g_z_0_zzzz_yyzzzz, g_z_0_zzzz_yyzzzzz, g_z_0_zzzz_yzzzzz, g_z_0_zzzz_yzzzzzz, g_z_0_zzzz_zzzzzz, g_z_0_zzzz_zzzzzzz, g_z_0_zzzzz_xxxxxx, g_z_0_zzzzz_xxxxxy, g_z_0_zzzzz_xxxxxz, g_z_0_zzzzz_xxxxyy, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxzz, g_z_0_zzzzz_xxxyyy, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxzzz, g_z_0_zzzzz_xxyyyy, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxzzzz, g_z_0_zzzzz_xyyyyy, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xzzzzz, g_z_0_zzzzz_yyyyyy, g_z_0_zzzzz_yyyyyz, g_z_0_zzzzz_yyyyzz, g_z_0_zzzzz_yyyzzz, g_z_0_zzzzz_yyzzzz, g_z_0_zzzzz_yzzzzz, g_z_0_zzzzz_zzzzzz, g_zzzz_xxxxxx, g_zzzz_xxxxxy, g_zzzz_xxxxxz, g_zzzz_xxxxyy, g_zzzz_xxxxyz, g_zzzz_xxxxzz, g_zzzz_xxxyyy, g_zzzz_xxxyyz, g_zzzz_xxxyzz, g_zzzz_xxxzzz, g_zzzz_xxyyyy, g_zzzz_xxyyyz, g_zzzz_xxyyzz, g_zzzz_xxyzzz, g_zzzz_xxzzzz, g_zzzz_xyyyyy, g_zzzz_xyyyyz, g_zzzz_xyyyzz, g_zzzz_xyyzzz, g_zzzz_xyzzzz, g_zzzz_xzzzzz, g_zzzz_yyyyyy, g_zzzz_yyyyyz, g_zzzz_yyyyzz, g_zzzz_yyyzzz, g_zzzz_yyzzzz, g_zzzz_yzzzzz, g_zzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzz_xxxxxx[k] = -g_zzzz_xxxxxx[k] - g_z_0_zzzz_xxxxxx[k] * ab_z + g_z_0_zzzz_xxxxxxz[k];

                g_z_0_zzzzz_xxxxxy[k] = -g_zzzz_xxxxxy[k] - g_z_0_zzzz_xxxxxy[k] * ab_z + g_z_0_zzzz_xxxxxyz[k];

                g_z_0_zzzzz_xxxxxz[k] = -g_zzzz_xxxxxz[k] - g_z_0_zzzz_xxxxxz[k] * ab_z + g_z_0_zzzz_xxxxxzz[k];

                g_z_0_zzzzz_xxxxyy[k] = -g_zzzz_xxxxyy[k] - g_z_0_zzzz_xxxxyy[k] * ab_z + g_z_0_zzzz_xxxxyyz[k];

                g_z_0_zzzzz_xxxxyz[k] = -g_zzzz_xxxxyz[k] - g_z_0_zzzz_xxxxyz[k] * ab_z + g_z_0_zzzz_xxxxyzz[k];

                g_z_0_zzzzz_xxxxzz[k] = -g_zzzz_xxxxzz[k] - g_z_0_zzzz_xxxxzz[k] * ab_z + g_z_0_zzzz_xxxxzzz[k];

                g_z_0_zzzzz_xxxyyy[k] = -g_zzzz_xxxyyy[k] - g_z_0_zzzz_xxxyyy[k] * ab_z + g_z_0_zzzz_xxxyyyz[k];

                g_z_0_zzzzz_xxxyyz[k] = -g_zzzz_xxxyyz[k] - g_z_0_zzzz_xxxyyz[k] * ab_z + g_z_0_zzzz_xxxyyzz[k];

                g_z_0_zzzzz_xxxyzz[k] = -g_zzzz_xxxyzz[k] - g_z_0_zzzz_xxxyzz[k] * ab_z + g_z_0_zzzz_xxxyzzz[k];

                g_z_0_zzzzz_xxxzzz[k] = -g_zzzz_xxxzzz[k] - g_z_0_zzzz_xxxzzz[k] * ab_z + g_z_0_zzzz_xxxzzzz[k];

                g_z_0_zzzzz_xxyyyy[k] = -g_zzzz_xxyyyy[k] - g_z_0_zzzz_xxyyyy[k] * ab_z + g_z_0_zzzz_xxyyyyz[k];

                g_z_0_zzzzz_xxyyyz[k] = -g_zzzz_xxyyyz[k] - g_z_0_zzzz_xxyyyz[k] * ab_z + g_z_0_zzzz_xxyyyzz[k];

                g_z_0_zzzzz_xxyyzz[k] = -g_zzzz_xxyyzz[k] - g_z_0_zzzz_xxyyzz[k] * ab_z + g_z_0_zzzz_xxyyzzz[k];

                g_z_0_zzzzz_xxyzzz[k] = -g_zzzz_xxyzzz[k] - g_z_0_zzzz_xxyzzz[k] * ab_z + g_z_0_zzzz_xxyzzzz[k];

                g_z_0_zzzzz_xxzzzz[k] = -g_zzzz_xxzzzz[k] - g_z_0_zzzz_xxzzzz[k] * ab_z + g_z_0_zzzz_xxzzzzz[k];

                g_z_0_zzzzz_xyyyyy[k] = -g_zzzz_xyyyyy[k] - g_z_0_zzzz_xyyyyy[k] * ab_z + g_z_0_zzzz_xyyyyyz[k];

                g_z_0_zzzzz_xyyyyz[k] = -g_zzzz_xyyyyz[k] - g_z_0_zzzz_xyyyyz[k] * ab_z + g_z_0_zzzz_xyyyyzz[k];

                g_z_0_zzzzz_xyyyzz[k] = -g_zzzz_xyyyzz[k] - g_z_0_zzzz_xyyyzz[k] * ab_z + g_z_0_zzzz_xyyyzzz[k];

                g_z_0_zzzzz_xyyzzz[k] = -g_zzzz_xyyzzz[k] - g_z_0_zzzz_xyyzzz[k] * ab_z + g_z_0_zzzz_xyyzzzz[k];

                g_z_0_zzzzz_xyzzzz[k] = -g_zzzz_xyzzzz[k] - g_z_0_zzzz_xyzzzz[k] * ab_z + g_z_0_zzzz_xyzzzzz[k];

                g_z_0_zzzzz_xzzzzz[k] = -g_zzzz_xzzzzz[k] - g_z_0_zzzz_xzzzzz[k] * ab_z + g_z_0_zzzz_xzzzzzz[k];

                g_z_0_zzzzz_yyyyyy[k] = -g_zzzz_yyyyyy[k] - g_z_0_zzzz_yyyyyy[k] * ab_z + g_z_0_zzzz_yyyyyyz[k];

                g_z_0_zzzzz_yyyyyz[k] = -g_zzzz_yyyyyz[k] - g_z_0_zzzz_yyyyyz[k] * ab_z + g_z_0_zzzz_yyyyyzz[k];

                g_z_0_zzzzz_yyyyzz[k] = -g_zzzz_yyyyzz[k] - g_z_0_zzzz_yyyyzz[k] * ab_z + g_z_0_zzzz_yyyyzzz[k];

                g_z_0_zzzzz_yyyzzz[k] = -g_zzzz_yyyzzz[k] - g_z_0_zzzz_yyyzzz[k] * ab_z + g_z_0_zzzz_yyyzzzz[k];

                g_z_0_zzzzz_yyzzzz[k] = -g_zzzz_yyzzzz[k] - g_z_0_zzzz_yyzzzz[k] * ab_z + g_z_0_zzzz_yyzzzzz[k];

                g_z_0_zzzzz_yzzzzz[k] = -g_zzzz_yzzzzz[k] - g_z_0_zzzz_yzzzzz[k] * ab_z + g_z_0_zzzz_yzzzzzz[k];

                g_z_0_zzzzz_zzzzzz[k] = -g_zzzz_zzzzzz[k] - g_z_0_zzzz_zzzzzz[k] * ab_z + g_z_0_zzzz_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

