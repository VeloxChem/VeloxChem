#include "ElectronRepulsionGeom1000ContrRecIIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_iixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_iixx,
                                            const size_t idx_hixx,
                                            const size_t idx_geom_10_hixx,
                                            const size_t idx_geom_10_hkxx,
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
            /// Set up components of auxilary buffer : HISS

            const auto hi_off = idx_hixx + i * dcomps + j;

            auto g_xxxxx_xxxxxx = cbuffer.data(hi_off + 0 * ccomps * dcomps);

            auto g_xxxxx_xxxxxy = cbuffer.data(hi_off + 1 * ccomps * dcomps);

            auto g_xxxxx_xxxxxz = cbuffer.data(hi_off + 2 * ccomps * dcomps);

            auto g_xxxxx_xxxxyy = cbuffer.data(hi_off + 3 * ccomps * dcomps);

            auto g_xxxxx_xxxxyz = cbuffer.data(hi_off + 4 * ccomps * dcomps);

            auto g_xxxxx_xxxxzz = cbuffer.data(hi_off + 5 * ccomps * dcomps);

            auto g_xxxxx_xxxyyy = cbuffer.data(hi_off + 6 * ccomps * dcomps);

            auto g_xxxxx_xxxyyz = cbuffer.data(hi_off + 7 * ccomps * dcomps);

            auto g_xxxxx_xxxyzz = cbuffer.data(hi_off + 8 * ccomps * dcomps);

            auto g_xxxxx_xxxzzz = cbuffer.data(hi_off + 9 * ccomps * dcomps);

            auto g_xxxxx_xxyyyy = cbuffer.data(hi_off + 10 * ccomps * dcomps);

            auto g_xxxxx_xxyyyz = cbuffer.data(hi_off + 11 * ccomps * dcomps);

            auto g_xxxxx_xxyyzz = cbuffer.data(hi_off + 12 * ccomps * dcomps);

            auto g_xxxxx_xxyzzz = cbuffer.data(hi_off + 13 * ccomps * dcomps);

            auto g_xxxxx_xxzzzz = cbuffer.data(hi_off + 14 * ccomps * dcomps);

            auto g_xxxxx_xyyyyy = cbuffer.data(hi_off + 15 * ccomps * dcomps);

            auto g_xxxxx_xyyyyz = cbuffer.data(hi_off + 16 * ccomps * dcomps);

            auto g_xxxxx_xyyyzz = cbuffer.data(hi_off + 17 * ccomps * dcomps);

            auto g_xxxxx_xyyzzz = cbuffer.data(hi_off + 18 * ccomps * dcomps);

            auto g_xxxxx_xyzzzz = cbuffer.data(hi_off + 19 * ccomps * dcomps);

            auto g_xxxxx_xzzzzz = cbuffer.data(hi_off + 20 * ccomps * dcomps);

            auto g_xxxxx_yyyyyy = cbuffer.data(hi_off + 21 * ccomps * dcomps);

            auto g_xxxxx_yyyyyz = cbuffer.data(hi_off + 22 * ccomps * dcomps);

            auto g_xxxxx_yyyyzz = cbuffer.data(hi_off + 23 * ccomps * dcomps);

            auto g_xxxxx_yyyzzz = cbuffer.data(hi_off + 24 * ccomps * dcomps);

            auto g_xxxxx_yyzzzz = cbuffer.data(hi_off + 25 * ccomps * dcomps);

            auto g_xxxxx_yzzzzz = cbuffer.data(hi_off + 26 * ccomps * dcomps);

            auto g_xxxxx_zzzzzz = cbuffer.data(hi_off + 27 * ccomps * dcomps);

            auto g_yyyyy_xxxxxx = cbuffer.data(hi_off + 420 * ccomps * dcomps);

            auto g_yyyyy_xxxxxy = cbuffer.data(hi_off + 421 * ccomps * dcomps);

            auto g_yyyyy_xxxxxz = cbuffer.data(hi_off + 422 * ccomps * dcomps);

            auto g_yyyyy_xxxxyy = cbuffer.data(hi_off + 423 * ccomps * dcomps);

            auto g_yyyyy_xxxxyz = cbuffer.data(hi_off + 424 * ccomps * dcomps);

            auto g_yyyyy_xxxxzz = cbuffer.data(hi_off + 425 * ccomps * dcomps);

            auto g_yyyyy_xxxyyy = cbuffer.data(hi_off + 426 * ccomps * dcomps);

            auto g_yyyyy_xxxyyz = cbuffer.data(hi_off + 427 * ccomps * dcomps);

            auto g_yyyyy_xxxyzz = cbuffer.data(hi_off + 428 * ccomps * dcomps);

            auto g_yyyyy_xxxzzz = cbuffer.data(hi_off + 429 * ccomps * dcomps);

            auto g_yyyyy_xxyyyy = cbuffer.data(hi_off + 430 * ccomps * dcomps);

            auto g_yyyyy_xxyyyz = cbuffer.data(hi_off + 431 * ccomps * dcomps);

            auto g_yyyyy_xxyyzz = cbuffer.data(hi_off + 432 * ccomps * dcomps);

            auto g_yyyyy_xxyzzz = cbuffer.data(hi_off + 433 * ccomps * dcomps);

            auto g_yyyyy_xxzzzz = cbuffer.data(hi_off + 434 * ccomps * dcomps);

            auto g_yyyyy_xyyyyy = cbuffer.data(hi_off + 435 * ccomps * dcomps);

            auto g_yyyyy_xyyyyz = cbuffer.data(hi_off + 436 * ccomps * dcomps);

            auto g_yyyyy_xyyyzz = cbuffer.data(hi_off + 437 * ccomps * dcomps);

            auto g_yyyyy_xyyzzz = cbuffer.data(hi_off + 438 * ccomps * dcomps);

            auto g_yyyyy_xyzzzz = cbuffer.data(hi_off + 439 * ccomps * dcomps);

            auto g_yyyyy_xzzzzz = cbuffer.data(hi_off + 440 * ccomps * dcomps);

            auto g_yyyyy_yyyyyy = cbuffer.data(hi_off + 441 * ccomps * dcomps);

            auto g_yyyyy_yyyyyz = cbuffer.data(hi_off + 442 * ccomps * dcomps);

            auto g_yyyyy_yyyyzz = cbuffer.data(hi_off + 443 * ccomps * dcomps);

            auto g_yyyyy_yyyzzz = cbuffer.data(hi_off + 444 * ccomps * dcomps);

            auto g_yyyyy_yyzzzz = cbuffer.data(hi_off + 445 * ccomps * dcomps);

            auto g_yyyyy_yzzzzz = cbuffer.data(hi_off + 446 * ccomps * dcomps);

            auto g_yyyyy_zzzzzz = cbuffer.data(hi_off + 447 * ccomps * dcomps);

            auto g_zzzzz_xxxxxx = cbuffer.data(hi_off + 560 * ccomps * dcomps);

            auto g_zzzzz_xxxxxy = cbuffer.data(hi_off + 561 * ccomps * dcomps);

            auto g_zzzzz_xxxxxz = cbuffer.data(hi_off + 562 * ccomps * dcomps);

            auto g_zzzzz_xxxxyy = cbuffer.data(hi_off + 563 * ccomps * dcomps);

            auto g_zzzzz_xxxxyz = cbuffer.data(hi_off + 564 * ccomps * dcomps);

            auto g_zzzzz_xxxxzz = cbuffer.data(hi_off + 565 * ccomps * dcomps);

            auto g_zzzzz_xxxyyy = cbuffer.data(hi_off + 566 * ccomps * dcomps);

            auto g_zzzzz_xxxyyz = cbuffer.data(hi_off + 567 * ccomps * dcomps);

            auto g_zzzzz_xxxyzz = cbuffer.data(hi_off + 568 * ccomps * dcomps);

            auto g_zzzzz_xxxzzz = cbuffer.data(hi_off + 569 * ccomps * dcomps);

            auto g_zzzzz_xxyyyy = cbuffer.data(hi_off + 570 * ccomps * dcomps);

            auto g_zzzzz_xxyyyz = cbuffer.data(hi_off + 571 * ccomps * dcomps);

            auto g_zzzzz_xxyyzz = cbuffer.data(hi_off + 572 * ccomps * dcomps);

            auto g_zzzzz_xxyzzz = cbuffer.data(hi_off + 573 * ccomps * dcomps);

            auto g_zzzzz_xxzzzz = cbuffer.data(hi_off + 574 * ccomps * dcomps);

            auto g_zzzzz_xyyyyy = cbuffer.data(hi_off + 575 * ccomps * dcomps);

            auto g_zzzzz_xyyyyz = cbuffer.data(hi_off + 576 * ccomps * dcomps);

            auto g_zzzzz_xyyyzz = cbuffer.data(hi_off + 577 * ccomps * dcomps);

            auto g_zzzzz_xyyzzz = cbuffer.data(hi_off + 578 * ccomps * dcomps);

            auto g_zzzzz_xyzzzz = cbuffer.data(hi_off + 579 * ccomps * dcomps);

            auto g_zzzzz_xzzzzz = cbuffer.data(hi_off + 580 * ccomps * dcomps);

            auto g_zzzzz_yyyyyy = cbuffer.data(hi_off + 581 * ccomps * dcomps);

            auto g_zzzzz_yyyyyz = cbuffer.data(hi_off + 582 * ccomps * dcomps);

            auto g_zzzzz_yyyyzz = cbuffer.data(hi_off + 583 * ccomps * dcomps);

            auto g_zzzzz_yyyzzz = cbuffer.data(hi_off + 584 * ccomps * dcomps);

            auto g_zzzzz_yyzzzz = cbuffer.data(hi_off + 585 * ccomps * dcomps);

            auto g_zzzzz_yzzzzz = cbuffer.data(hi_off + 586 * ccomps * dcomps);

            auto g_zzzzz_zzzzzz = cbuffer.data(hi_off + 587 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HISS

            const auto hi_geom_10_off = idx_geom_10_hixx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : HKSS

            const auto hk_geom_10_off = idx_geom_10_hkxx + i * dcomps + j;

            auto g_x_0_xxxxx_xxxxxxx = cbuffer.data(hk_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxxxy = cbuffer.data(hk_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxxxz = cbuffer.data(hk_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxxyy = cbuffer.data(hk_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxxyz = cbuffer.data(hk_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxxzz = cbuffer.data(hk_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxyyy = cbuffer.data(hk_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxyyz = cbuffer.data(hk_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxyzz = cbuffer.data(hk_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxzzz = cbuffer.data(hk_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyyyy = cbuffer.data(hk_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyyyz = cbuffer.data(hk_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyyzz = cbuffer.data(hk_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyzzz = cbuffer.data(hk_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxzzzz = cbuffer.data(hk_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyyyy = cbuffer.data(hk_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyyyz = cbuffer.data(hk_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyyzz = cbuffer.data(hk_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyzzz = cbuffer.data(hk_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyzzzz = cbuffer.data(hk_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxzzzzz = cbuffer.data(hk_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyyyy = cbuffer.data(hk_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyyyz = cbuffer.data(hk_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyyzz = cbuffer.data(hk_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyzzz = cbuffer.data(hk_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyzzzz = cbuffer.data(hk_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyzzzzz = cbuffer.data(hk_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxx_xzzzzzz = cbuffer.data(hk_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyyyy = cbuffer.data(hk_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyyyz = cbuffer.data(hk_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyyzz = cbuffer.data(hk_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyzzz = cbuffer.data(hk_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyzzzz = cbuffer.data(hk_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyzzzzz = cbuffer.data(hk_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxxx_yzzzzzz = cbuffer.data(hk_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxxx_zzzzzzz = cbuffer.data(hk_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxxxy = cbuffer.data(hk_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxxyy = cbuffer.data(hk_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxxyz = cbuffer.data(hk_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxyyy = cbuffer.data(hk_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxyyz = cbuffer.data(hk_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxyzz = cbuffer.data(hk_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyyyy = cbuffer.data(hk_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyyyz = cbuffer.data(hk_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyyzz = cbuffer.data(hk_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyzzz = cbuffer.data(hk_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyyyy = cbuffer.data(hk_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyyyz = cbuffer.data(hk_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyyzz = cbuffer.data(hk_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyzzz = cbuffer.data(hk_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyzzzz = cbuffer.data(hk_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyyyy = cbuffer.data(hk_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyyyz = cbuffer.data(hk_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyyzz = cbuffer.data(hk_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyzzz = cbuffer.data(hk_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyzzzz = cbuffer.data(hk_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyzzzzz = cbuffer.data(hk_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyyyy = cbuffer.data(hk_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyyyz = cbuffer.data(hk_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyyzz = cbuffer.data(hk_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyzzz = cbuffer.data(hk_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyzzzz = cbuffer.data(hk_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyzzzzz = cbuffer.data(hk_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxxy_yzzzzzz = cbuffer.data(hk_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxxxy = cbuffer.data(hk_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxxxz = cbuffer.data(hk_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxxyy = cbuffer.data(hk_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxxyz = cbuffer.data(hk_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxxzz = cbuffer.data(hk_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxyyy = cbuffer.data(hk_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxyyz = cbuffer.data(hk_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxyzz = cbuffer.data(hk_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxzzz = cbuffer.data(hk_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyyyy = cbuffer.data(hk_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyyyz = cbuffer.data(hk_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyyzz = cbuffer.data(hk_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyzzz = cbuffer.data(hk_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxzzzz = cbuffer.data(hk_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyyyy = cbuffer.data(hk_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyyyz = cbuffer.data(hk_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyyzz = cbuffer.data(hk_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyzzz = cbuffer.data(hk_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyzzzz = cbuffer.data(hk_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxzzzzz = cbuffer.data(hk_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyyyy = cbuffer.data(hk_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyyyz = cbuffer.data(hk_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyyzz = cbuffer.data(hk_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyzzz = cbuffer.data(hk_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyzzzz = cbuffer.data(hk_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyzzzzz = cbuffer.data(hk_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxxxz_xzzzzzz = cbuffer.data(hk_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyyyy = cbuffer.data(hk_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyyyz = cbuffer.data(hk_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyyzz = cbuffer.data(hk_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyzzz = cbuffer.data(hk_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyzzzz = cbuffer.data(hk_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyzzzzz = cbuffer.data(hk_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxxxz_yzzzzzz = cbuffer.data(hk_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxxxz_zzzzzzz = cbuffer.data(hk_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyyyy = cbuffer.data(hk_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyyyz = cbuffer.data(hk_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyyzz = cbuffer.data(hk_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyzzz = cbuffer.data(hk_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyzzzz = cbuffer.data(hk_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyzzzzz = cbuffer.data(hk_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxxyy_yzzzzzz = cbuffer.data(hk_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyyyy = cbuffer.data(hk_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyyyz = cbuffer.data(hk_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyyzz = cbuffer.data(hk_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyzzz = cbuffer.data(hk_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyzzzz = cbuffer.data(hk_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyzzzzz = cbuffer.data(hk_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xxxyz_yzzzzzz = cbuffer.data(hk_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xxxzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xxxzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xxxzz_zzzzzzz = cbuffer.data(hk_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyyyy = cbuffer.data(hk_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyyyz = cbuffer.data(hk_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyyzz = cbuffer.data(hk_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyzzz = cbuffer.data(hk_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyzzzz = cbuffer.data(hk_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyzzzzz = cbuffer.data(hk_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xxyyy_yzzzzzz = cbuffer.data(hk_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyyyy = cbuffer.data(hk_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyyyz = cbuffer.data(hk_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyyzz = cbuffer.data(hk_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyzzz = cbuffer.data(hk_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyzzzz = cbuffer.data(hk_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyzzzzz = cbuffer.data(hk_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_xxyyz_yzzzzzz = cbuffer.data(hk_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 307 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_xxyzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 335 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_xxzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 356 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_xxzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_xxzzz_zzzzzzz = cbuffer.data(hk_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 363 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 364 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 366 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 373 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 378 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 384 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyyyy = cbuffer.data(hk_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyyyz = cbuffer.data(hk_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyyzz = cbuffer.data(hk_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyzzz = cbuffer.data(hk_geom_10_off + 391 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyzzzz = cbuffer.data(hk_geom_10_off + 392 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyzzzzz = cbuffer.data(hk_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_xyyyy_yzzzzzz = cbuffer.data(hk_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 404 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 419 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 420 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 421 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 422 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyyyy = cbuffer.data(hk_geom_10_off + 424 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyyyz = cbuffer.data(hk_geom_10_off + 425 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyyzz = cbuffer.data(hk_geom_10_off + 426 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyzzz = cbuffer.data(hk_geom_10_off + 427 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyzzzz = cbuffer.data(hk_geom_10_off + 428 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyzzzzz = cbuffer.data(hk_geom_10_off + 429 * ccomps * dcomps);

            auto g_x_0_xyyyz_yzzzzzz = cbuffer.data(hk_geom_10_off + 430 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 433 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 435 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 436 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 438 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 439 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 440 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 442 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 443 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 444 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 445 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 447 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 448 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 449 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 450 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 451 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 453 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 454 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 455 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 456 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 457 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 458 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 460 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 461 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 462 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 463 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 464 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 465 * ccomps * dcomps);

            auto g_x_0_xyyzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 466 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 469 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 471 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 472 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 474 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 475 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 476 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 478 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 479 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 480 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 481 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 483 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 484 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 485 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 486 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 487 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 489 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 490 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 491 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 492 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 493 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 494 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 496 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 497 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 498 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 499 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 500 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 501 * ccomps * dcomps);

            auto g_x_0_xyzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 502 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 505 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 506 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 507 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 508 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 509 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 510 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 511 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 512 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 513 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 514 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 515 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 516 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 517 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 518 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 519 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 520 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 521 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 522 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 523 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 524 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 525 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 526 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 527 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 528 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 529 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 530 * ccomps * dcomps);

            auto g_x_0_xzzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 531 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 532 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 533 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 534 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 535 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 536 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 537 * ccomps * dcomps);

            auto g_x_0_xzzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 538 * ccomps * dcomps);

            auto g_x_0_xzzzz_zzzzzzz = cbuffer.data(hk_geom_10_off + 539 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 541 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 543 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 544 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 546 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 547 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 548 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 550 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 551 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 552 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 553 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 555 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 556 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 557 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 558 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 559 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 561 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 562 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 563 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 564 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 565 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 566 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyyyy = cbuffer.data(hk_geom_10_off + 568 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyyyz = cbuffer.data(hk_geom_10_off + 569 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyyzz = cbuffer.data(hk_geom_10_off + 570 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyzzz = cbuffer.data(hk_geom_10_off + 571 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyzzzz = cbuffer.data(hk_geom_10_off + 572 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyzzzzz = cbuffer.data(hk_geom_10_off + 573 * ccomps * dcomps);

            auto g_x_0_yyyyy_yzzzzzz = cbuffer.data(hk_geom_10_off + 574 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 577 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 579 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 580 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 582 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 583 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 584 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 586 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 587 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 588 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 589 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 591 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 592 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 593 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 594 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 595 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 597 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 598 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 599 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 600 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 601 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 602 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyyyy = cbuffer.data(hk_geom_10_off + 604 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyyyz = cbuffer.data(hk_geom_10_off + 605 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyyzz = cbuffer.data(hk_geom_10_off + 606 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyzzz = cbuffer.data(hk_geom_10_off + 607 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyzzzz = cbuffer.data(hk_geom_10_off + 608 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyzzzzz = cbuffer.data(hk_geom_10_off + 609 * ccomps * dcomps);

            auto g_x_0_yyyyz_yzzzzzz = cbuffer.data(hk_geom_10_off + 610 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 613 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 615 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 616 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 618 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 619 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 620 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 622 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 623 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 624 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 625 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 627 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 628 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 629 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 630 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 631 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 633 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 634 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 635 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 636 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 637 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 638 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 640 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 641 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 642 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 643 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 644 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 645 * ccomps * dcomps);

            auto g_x_0_yyyzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 646 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 649 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 651 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 652 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 654 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 655 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 656 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 658 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 659 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 660 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 661 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 663 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 664 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 665 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 666 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 667 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 669 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 670 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 671 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 672 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 673 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 674 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 676 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 677 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 678 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 679 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 680 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 681 * ccomps * dcomps);

            auto g_x_0_yyzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 682 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 685 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 687 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 688 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 690 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 691 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 692 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 694 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 695 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 696 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 697 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 699 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 700 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 701 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 702 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 703 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 705 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 706 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 707 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 708 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 709 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 710 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 712 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 713 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 714 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 715 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 716 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 717 * ccomps * dcomps);

            auto g_x_0_yzzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 718 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 721 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 722 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 723 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 724 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 725 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 726 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 727 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 728 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 729 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 730 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 731 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 732 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 733 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 734 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 735 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 736 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 737 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 738 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 739 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 740 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 741 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 742 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 743 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 744 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 745 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 746 * ccomps * dcomps);

            auto g_x_0_zzzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 747 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 748 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 749 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 750 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 751 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 752 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 753 * ccomps * dcomps);

            auto g_x_0_zzzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 754 * ccomps * dcomps);

            auto g_x_0_zzzzz_zzzzzzz = cbuffer.data(hk_geom_10_off + 755 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxxx = cbuffer.data(hk_geom_10_off + 756 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxxy = cbuffer.data(hk_geom_10_off + 757 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxxz = cbuffer.data(hk_geom_10_off + 758 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxyy = cbuffer.data(hk_geom_10_off + 759 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxyz = cbuffer.data(hk_geom_10_off + 760 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxxzz = cbuffer.data(hk_geom_10_off + 761 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxyyy = cbuffer.data(hk_geom_10_off + 762 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxyyz = cbuffer.data(hk_geom_10_off + 763 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxyzz = cbuffer.data(hk_geom_10_off + 764 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxzzz = cbuffer.data(hk_geom_10_off + 765 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyyyy = cbuffer.data(hk_geom_10_off + 766 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyyyz = cbuffer.data(hk_geom_10_off + 767 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyyzz = cbuffer.data(hk_geom_10_off + 768 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyzzz = cbuffer.data(hk_geom_10_off + 769 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxzzzz = cbuffer.data(hk_geom_10_off + 770 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyyyy = cbuffer.data(hk_geom_10_off + 771 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyyyz = cbuffer.data(hk_geom_10_off + 772 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyyzz = cbuffer.data(hk_geom_10_off + 773 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyzzz = cbuffer.data(hk_geom_10_off + 774 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyzzzz = cbuffer.data(hk_geom_10_off + 775 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxzzzzz = cbuffer.data(hk_geom_10_off + 776 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyyyy = cbuffer.data(hk_geom_10_off + 777 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyyyz = cbuffer.data(hk_geom_10_off + 778 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyyzz = cbuffer.data(hk_geom_10_off + 779 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyzzz = cbuffer.data(hk_geom_10_off + 780 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyzzzz = cbuffer.data(hk_geom_10_off + 781 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyzzzzz = cbuffer.data(hk_geom_10_off + 782 * ccomps * dcomps);

            auto g_y_0_xxxxx_xzzzzzz = cbuffer.data(hk_geom_10_off + 783 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxxx = cbuffer.data(hk_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxxy = cbuffer.data(hk_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxxz = cbuffer.data(hk_geom_10_off + 794 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxyy = cbuffer.data(hk_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxyz = cbuffer.data(hk_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxxzz = cbuffer.data(hk_geom_10_off + 797 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxyyy = cbuffer.data(hk_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxyyz = cbuffer.data(hk_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxyzz = cbuffer.data(hk_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxzzz = cbuffer.data(hk_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyyyy = cbuffer.data(hk_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyyyz = cbuffer.data(hk_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyyzz = cbuffer.data(hk_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyzzz = cbuffer.data(hk_geom_10_off + 805 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxzzzz = cbuffer.data(hk_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyyyy = cbuffer.data(hk_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyyyz = cbuffer.data(hk_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyyzz = cbuffer.data(hk_geom_10_off + 809 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyzzz = cbuffer.data(hk_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyzzzz = cbuffer.data(hk_geom_10_off + 811 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxzzzzz = cbuffer.data(hk_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyyyy = cbuffer.data(hk_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyyyz = cbuffer.data(hk_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyyzz = cbuffer.data(hk_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyzzz = cbuffer.data(hk_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyzzzz = cbuffer.data(hk_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyzzzzz = cbuffer.data(hk_geom_10_off + 818 * ccomps * dcomps);

            auto g_y_0_xxxxy_xzzzzzz = cbuffer.data(hk_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxxx = cbuffer.data(hk_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxxy = cbuffer.data(hk_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxxz = cbuffer.data(hk_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxyy = cbuffer.data(hk_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxyz = cbuffer.data(hk_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxxzz = cbuffer.data(hk_geom_10_off + 833 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxyyy = cbuffer.data(hk_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxyyz = cbuffer.data(hk_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxyzz = cbuffer.data(hk_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxzzz = cbuffer.data(hk_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyyyy = cbuffer.data(hk_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyyyz = cbuffer.data(hk_geom_10_off + 839 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyyzz = cbuffer.data(hk_geom_10_off + 840 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyzzz = cbuffer.data(hk_geom_10_off + 841 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxzzzz = cbuffer.data(hk_geom_10_off + 842 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyyyy = cbuffer.data(hk_geom_10_off + 843 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyyyz = cbuffer.data(hk_geom_10_off + 844 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyyzz = cbuffer.data(hk_geom_10_off + 845 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyzzz = cbuffer.data(hk_geom_10_off + 846 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyzzzz = cbuffer.data(hk_geom_10_off + 847 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxzzzzz = cbuffer.data(hk_geom_10_off + 848 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyyyy = cbuffer.data(hk_geom_10_off + 849 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyyyz = cbuffer.data(hk_geom_10_off + 850 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyyzz = cbuffer.data(hk_geom_10_off + 851 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyzzz = cbuffer.data(hk_geom_10_off + 852 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyzzzz = cbuffer.data(hk_geom_10_off + 853 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyzzzzz = cbuffer.data(hk_geom_10_off + 854 * ccomps * dcomps);

            auto g_y_0_xxxxz_xzzzzzz = cbuffer.data(hk_geom_10_off + 855 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxxx = cbuffer.data(hk_geom_10_off + 864 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 865 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxxz = cbuffer.data(hk_geom_10_off + 866 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 867 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 868 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxxzz = cbuffer.data(hk_geom_10_off + 869 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 870 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 871 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 872 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxzzz = cbuffer.data(hk_geom_10_off + 873 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 874 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 875 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 876 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 877 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxzzzz = cbuffer.data(hk_geom_10_off + 878 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 879 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 880 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 881 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 882 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 883 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxzzzzz = cbuffer.data(hk_geom_10_off + 884 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 885 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 886 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 887 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 888 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 889 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 890 * ccomps * dcomps);

            auto g_y_0_xxxyy_xzzzzzz = cbuffer.data(hk_geom_10_off + 891 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxxx = cbuffer.data(hk_geom_10_off + 900 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 901 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxxz = cbuffer.data(hk_geom_10_off + 902 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 903 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 904 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxxzz = cbuffer.data(hk_geom_10_off + 905 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 906 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 907 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 908 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxzzz = cbuffer.data(hk_geom_10_off + 909 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 910 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 911 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 912 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 913 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxzzzz = cbuffer.data(hk_geom_10_off + 914 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 915 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 916 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 917 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 918 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 919 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxzzzzz = cbuffer.data(hk_geom_10_off + 920 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 921 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 922 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 923 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 924 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 925 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 926 * ccomps * dcomps);

            auto g_y_0_xxxyz_xzzzzzz = cbuffer.data(hk_geom_10_off + 927 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 936 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 937 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 938 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 939 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 940 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 941 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 942 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 943 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 944 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 945 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 946 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 947 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 948 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 949 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 950 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 951 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 952 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 953 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 954 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 955 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 956 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 957 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 958 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 959 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 960 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 961 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 962 * ccomps * dcomps);

            auto g_y_0_xxxzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 963 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxxx = cbuffer.data(hk_geom_10_off + 972 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 973 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxxz = cbuffer.data(hk_geom_10_off + 974 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 975 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 976 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxxzz = cbuffer.data(hk_geom_10_off + 977 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 978 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 979 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 980 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxzzz = cbuffer.data(hk_geom_10_off + 981 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 982 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 983 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 984 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 985 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxzzzz = cbuffer.data(hk_geom_10_off + 986 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 987 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 988 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 989 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 990 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 991 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxzzzzz = cbuffer.data(hk_geom_10_off + 992 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 993 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 994 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 995 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 996 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 997 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 998 * ccomps * dcomps);

            auto g_y_0_xxyyy_xzzzzzz = cbuffer.data(hk_geom_10_off + 999 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1008 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1009 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1010 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1011 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1012 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1013 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1014 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1015 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1016 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1017 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1018 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1019 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1020 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1021 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1022 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1023 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1024 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1025 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1026 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1027 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1028 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1029 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1030 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1031 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1032 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1033 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1034 * ccomps * dcomps);

            auto g_y_0_xxyyz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1035 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1044 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1045 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1046 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1047 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1048 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1049 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1050 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1051 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1052 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1053 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1054 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1055 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1056 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1057 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1058 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1059 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1060 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1061 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1062 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1063 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1064 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1065 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1066 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1067 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1068 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1069 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1070 * ccomps * dcomps);

            auto g_y_0_xxyzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1071 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1080 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1081 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1082 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1083 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1084 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1085 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1086 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1087 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1088 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1089 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1090 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1091 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1092 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1093 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1094 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1095 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1096 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1097 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1098 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1099 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1100 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1101 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1102 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1103 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1104 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1105 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1106 * ccomps * dcomps);

            auto g_y_0_xxzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1107 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxxx = cbuffer.data(hk_geom_10_off + 1116 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 1117 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxxz = cbuffer.data(hk_geom_10_off + 1118 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 1119 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 1120 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxxzz = cbuffer.data(hk_geom_10_off + 1121 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 1122 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 1123 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 1124 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxzzz = cbuffer.data(hk_geom_10_off + 1125 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 1126 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 1127 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 1128 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 1129 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxzzzz = cbuffer.data(hk_geom_10_off + 1130 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 1131 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 1132 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 1133 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 1134 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 1135 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxzzzzz = cbuffer.data(hk_geom_10_off + 1136 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 1137 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 1138 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 1139 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 1140 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 1141 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 1142 * ccomps * dcomps);

            auto g_y_0_xyyyy_xzzzzzz = cbuffer.data(hk_geom_10_off + 1143 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1152 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1153 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1154 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1155 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1156 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1157 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1158 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1159 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1160 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1161 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1162 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1163 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1164 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1165 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1166 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1167 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1168 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1169 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1170 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1171 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1172 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1173 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1174 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1175 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1176 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1177 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1178 * ccomps * dcomps);

            auto g_y_0_xyyyz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1179 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1188 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1189 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1190 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1191 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1192 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1193 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1194 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1195 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1196 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1197 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1198 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1199 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1200 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1201 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1202 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1203 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1204 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1205 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1206 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1207 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1208 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1209 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1210 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1211 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1212 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1213 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1214 * ccomps * dcomps);

            auto g_y_0_xyyzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1215 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1224 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1225 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1226 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1227 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1228 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1229 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1230 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1231 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1232 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1233 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1234 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1235 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1236 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1237 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1238 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1239 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1240 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1241 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1242 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1243 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1244 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1245 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1246 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1247 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1248 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1249 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1250 * ccomps * dcomps);

            auto g_y_0_xyzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1251 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1260 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1261 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1262 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1263 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1264 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1265 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1266 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1267 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1268 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1269 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1270 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1271 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1272 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1273 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1274 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1275 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1276 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1277 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1278 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1279 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1280 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1281 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1282 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1283 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1284 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1285 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1286 * ccomps * dcomps);

            auto g_y_0_xzzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1287 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxxx = cbuffer.data(hk_geom_10_off + 1296 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 1297 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxxz = cbuffer.data(hk_geom_10_off + 1298 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 1299 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 1300 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxxzz = cbuffer.data(hk_geom_10_off + 1301 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 1302 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 1303 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 1304 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxzzz = cbuffer.data(hk_geom_10_off + 1305 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 1306 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 1307 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 1308 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 1309 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxzzzz = cbuffer.data(hk_geom_10_off + 1310 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 1311 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 1312 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 1313 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 1314 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 1315 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxzzzzz = cbuffer.data(hk_geom_10_off + 1316 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 1317 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 1318 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 1319 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 1320 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 1321 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 1322 * ccomps * dcomps);

            auto g_y_0_yyyyy_xzzzzzz = cbuffer.data(hk_geom_10_off + 1323 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyyyy = cbuffer.data(hk_geom_10_off + 1324 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyyyz = cbuffer.data(hk_geom_10_off + 1325 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyyzz = cbuffer.data(hk_geom_10_off + 1326 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyzzz = cbuffer.data(hk_geom_10_off + 1327 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyzzzz = cbuffer.data(hk_geom_10_off + 1328 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyzzzzz = cbuffer.data(hk_geom_10_off + 1329 * ccomps * dcomps);

            auto g_y_0_yyyyy_yzzzzzz = cbuffer.data(hk_geom_10_off + 1330 * ccomps * dcomps);

            auto g_y_0_yyyyy_zzzzzzz = cbuffer.data(hk_geom_10_off + 1331 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1332 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1333 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1334 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1335 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1336 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1337 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1338 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1339 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1340 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1341 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1342 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1343 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1344 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1345 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1346 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1347 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1348 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1349 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1350 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1351 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1352 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1353 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1354 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1355 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1356 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1357 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1358 * ccomps * dcomps);

            auto g_y_0_yyyyz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1359 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyyyz = cbuffer.data(hk_geom_10_off + 1361 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyyzz = cbuffer.data(hk_geom_10_off + 1362 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyzzz = cbuffer.data(hk_geom_10_off + 1363 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyzzzz = cbuffer.data(hk_geom_10_off + 1364 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyzzzzz = cbuffer.data(hk_geom_10_off + 1365 * ccomps * dcomps);

            auto g_y_0_yyyyz_yzzzzzz = cbuffer.data(hk_geom_10_off + 1366 * ccomps * dcomps);

            auto g_y_0_yyyyz_zzzzzzz = cbuffer.data(hk_geom_10_off + 1367 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1368 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1369 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1370 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1371 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1372 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1373 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1374 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1375 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1376 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1377 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1378 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1379 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1380 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1381 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1382 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1383 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1384 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1385 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1386 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1387 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1388 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1389 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1390 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1391 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1392 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1393 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1394 * ccomps * dcomps);

            auto g_y_0_yyyzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1395 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 1397 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 1398 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 1399 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 1400 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 1401 * ccomps * dcomps);

            auto g_y_0_yyyzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 1402 * ccomps * dcomps);

            auto g_y_0_yyyzz_zzzzzzz = cbuffer.data(hk_geom_10_off + 1403 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1404 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1405 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1406 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1407 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1408 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1409 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1410 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1411 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1412 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1413 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1414 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1415 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1416 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1417 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1418 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1419 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1420 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1421 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1422 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1423 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1424 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1425 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1426 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1427 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1428 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1429 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1430 * ccomps * dcomps);

            auto g_y_0_yyzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1431 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 1433 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 1434 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 1435 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 1436 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 1437 * ccomps * dcomps);

            auto g_y_0_yyzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 1438 * ccomps * dcomps);

            auto g_y_0_yyzzz_zzzzzzz = cbuffer.data(hk_geom_10_off + 1439 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1440 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1441 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1442 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1443 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1444 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1445 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1446 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1447 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1448 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1449 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1450 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1451 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1452 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1453 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1454 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1455 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1456 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1457 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1458 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1459 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1460 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1461 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1462 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1463 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1464 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1465 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1466 * ccomps * dcomps);

            auto g_y_0_yzzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1467 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 1469 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 1470 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 1471 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 1472 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 1473 * ccomps * dcomps);

            auto g_y_0_yzzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 1474 * ccomps * dcomps);

            auto g_y_0_yzzzz_zzzzzzz = cbuffer.data(hk_geom_10_off + 1475 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1476 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1477 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1478 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1479 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1480 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1481 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1482 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1483 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1484 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1485 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1486 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1487 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1488 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1489 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1490 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1491 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1492 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1493 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1494 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1495 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1496 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1497 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1498 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1499 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1500 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1501 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1502 * ccomps * dcomps);

            auto g_y_0_zzzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1503 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 1505 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 1506 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 1507 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 1508 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 1509 * ccomps * dcomps);

            auto g_y_0_zzzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 1510 * ccomps * dcomps);

            auto g_y_0_zzzzz_zzzzzzz = cbuffer.data(hk_geom_10_off + 1511 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxxx = cbuffer.data(hk_geom_10_off + 1512 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxxy = cbuffer.data(hk_geom_10_off + 1513 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxxz = cbuffer.data(hk_geom_10_off + 1514 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxyy = cbuffer.data(hk_geom_10_off + 1515 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxyz = cbuffer.data(hk_geom_10_off + 1516 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxxzz = cbuffer.data(hk_geom_10_off + 1517 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxyyy = cbuffer.data(hk_geom_10_off + 1518 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxyyz = cbuffer.data(hk_geom_10_off + 1519 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxyzz = cbuffer.data(hk_geom_10_off + 1520 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxzzz = cbuffer.data(hk_geom_10_off + 1521 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyyyy = cbuffer.data(hk_geom_10_off + 1522 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyyyz = cbuffer.data(hk_geom_10_off + 1523 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyyzz = cbuffer.data(hk_geom_10_off + 1524 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyzzz = cbuffer.data(hk_geom_10_off + 1525 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxzzzz = cbuffer.data(hk_geom_10_off + 1526 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyyyy = cbuffer.data(hk_geom_10_off + 1527 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyyyz = cbuffer.data(hk_geom_10_off + 1528 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyyzz = cbuffer.data(hk_geom_10_off + 1529 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyzzz = cbuffer.data(hk_geom_10_off + 1530 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyzzzz = cbuffer.data(hk_geom_10_off + 1531 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxzzzzz = cbuffer.data(hk_geom_10_off + 1532 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyyyy = cbuffer.data(hk_geom_10_off + 1533 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyyyz = cbuffer.data(hk_geom_10_off + 1534 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyyzz = cbuffer.data(hk_geom_10_off + 1535 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyzzz = cbuffer.data(hk_geom_10_off + 1536 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyzzzz = cbuffer.data(hk_geom_10_off + 1537 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyzzzzz = cbuffer.data(hk_geom_10_off + 1538 * ccomps * dcomps);

            auto g_z_0_xxxxx_xzzzzzz = cbuffer.data(hk_geom_10_off + 1539 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxxx = cbuffer.data(hk_geom_10_off + 1548 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxxy = cbuffer.data(hk_geom_10_off + 1549 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxxz = cbuffer.data(hk_geom_10_off + 1550 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxyy = cbuffer.data(hk_geom_10_off + 1551 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxyz = cbuffer.data(hk_geom_10_off + 1552 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxxzz = cbuffer.data(hk_geom_10_off + 1553 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxyyy = cbuffer.data(hk_geom_10_off + 1554 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxyyz = cbuffer.data(hk_geom_10_off + 1555 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxyzz = cbuffer.data(hk_geom_10_off + 1556 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxzzz = cbuffer.data(hk_geom_10_off + 1557 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyyyy = cbuffer.data(hk_geom_10_off + 1558 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyyyz = cbuffer.data(hk_geom_10_off + 1559 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyyzz = cbuffer.data(hk_geom_10_off + 1560 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyzzz = cbuffer.data(hk_geom_10_off + 1561 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxzzzz = cbuffer.data(hk_geom_10_off + 1562 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyyyy = cbuffer.data(hk_geom_10_off + 1563 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyyyz = cbuffer.data(hk_geom_10_off + 1564 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyyzz = cbuffer.data(hk_geom_10_off + 1565 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyzzz = cbuffer.data(hk_geom_10_off + 1566 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyzzzz = cbuffer.data(hk_geom_10_off + 1567 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxzzzzz = cbuffer.data(hk_geom_10_off + 1568 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyyyy = cbuffer.data(hk_geom_10_off + 1569 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyyyz = cbuffer.data(hk_geom_10_off + 1570 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyyzz = cbuffer.data(hk_geom_10_off + 1571 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyzzz = cbuffer.data(hk_geom_10_off + 1572 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyzzzz = cbuffer.data(hk_geom_10_off + 1573 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyzzzzz = cbuffer.data(hk_geom_10_off + 1574 * ccomps * dcomps);

            auto g_z_0_xxxxy_xzzzzzz = cbuffer.data(hk_geom_10_off + 1575 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1584 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1585 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1586 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1587 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1588 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1589 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1590 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1591 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1592 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1593 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1594 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1595 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1596 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1597 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1598 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1599 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1600 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1601 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1602 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1603 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1604 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1605 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1606 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1607 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1608 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1609 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1610 * ccomps * dcomps);

            auto g_z_0_xxxxz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1611 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxxx = cbuffer.data(hk_geom_10_off + 1620 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 1621 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxxz = cbuffer.data(hk_geom_10_off + 1622 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 1623 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 1624 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxxzz = cbuffer.data(hk_geom_10_off + 1625 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 1626 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 1627 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 1628 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxzzz = cbuffer.data(hk_geom_10_off + 1629 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 1630 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 1631 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 1632 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 1633 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxzzzz = cbuffer.data(hk_geom_10_off + 1634 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 1635 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 1636 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 1637 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 1638 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 1639 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxzzzzz = cbuffer.data(hk_geom_10_off + 1640 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 1641 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 1642 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 1643 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 1644 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 1645 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 1646 * ccomps * dcomps);

            auto g_z_0_xxxyy_xzzzzzz = cbuffer.data(hk_geom_10_off + 1647 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1656 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1657 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1658 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1659 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1660 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1661 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1662 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1663 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1664 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1665 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1666 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1667 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1668 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1669 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1670 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1671 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1672 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1673 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1674 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1675 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1676 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1677 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1678 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1679 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1680 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1681 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1682 * ccomps * dcomps);

            auto g_z_0_xxxyz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1683 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1692 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1693 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1694 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1695 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1696 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1697 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1698 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1699 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1700 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1701 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1702 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1703 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1704 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1705 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1706 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1707 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1708 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1709 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1710 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1711 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1712 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1713 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1714 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1715 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1716 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1717 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1718 * ccomps * dcomps);

            auto g_z_0_xxxzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1719 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxxx = cbuffer.data(hk_geom_10_off + 1728 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 1729 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxxz = cbuffer.data(hk_geom_10_off + 1730 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 1731 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 1732 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxxzz = cbuffer.data(hk_geom_10_off + 1733 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 1734 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 1735 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 1736 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxzzz = cbuffer.data(hk_geom_10_off + 1737 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 1738 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 1739 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 1740 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 1741 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxzzzz = cbuffer.data(hk_geom_10_off + 1742 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 1743 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 1744 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 1745 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 1746 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 1747 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxzzzzz = cbuffer.data(hk_geom_10_off + 1748 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 1749 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 1750 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 1751 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 1752 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 1753 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 1754 * ccomps * dcomps);

            auto g_z_0_xxyyy_xzzzzzz = cbuffer.data(hk_geom_10_off + 1755 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1764 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1765 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1766 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1767 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1768 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1769 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1770 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1771 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1772 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1773 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1774 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1775 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1776 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1777 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1778 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1779 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1780 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1781 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1782 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1783 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1784 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1785 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1786 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1787 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1788 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1789 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1790 * ccomps * dcomps);

            auto g_z_0_xxyyz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1791 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1800 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1801 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1802 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1803 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1804 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1805 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1806 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1807 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1808 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1809 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1810 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1811 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1812 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1813 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1814 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1815 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1816 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1817 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1818 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1819 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1820 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1821 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1822 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1823 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1824 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1825 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1826 * ccomps * dcomps);

            auto g_z_0_xxyzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1827 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1836 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1837 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1838 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1839 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1840 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1841 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1842 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1843 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1844 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1845 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1846 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1847 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1848 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1849 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1850 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1851 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1852 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1853 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1854 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1855 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1856 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1857 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1858 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1859 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1860 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1861 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1862 * ccomps * dcomps);

            auto g_z_0_xxzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1863 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxxx = cbuffer.data(hk_geom_10_off + 1872 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 1873 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxxz = cbuffer.data(hk_geom_10_off + 1874 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 1875 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 1876 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxxzz = cbuffer.data(hk_geom_10_off + 1877 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 1878 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 1879 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 1880 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxzzz = cbuffer.data(hk_geom_10_off + 1881 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 1882 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 1883 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 1884 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 1885 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxzzzz = cbuffer.data(hk_geom_10_off + 1886 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 1887 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 1888 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 1889 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 1890 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 1891 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxzzzzz = cbuffer.data(hk_geom_10_off + 1892 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 1893 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 1894 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 1895 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 1896 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 1897 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 1898 * ccomps * dcomps);

            auto g_z_0_xyyyy_xzzzzzz = cbuffer.data(hk_geom_10_off + 1899 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1908 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1909 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1910 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1911 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1912 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1913 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1914 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1915 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1916 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1917 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1918 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1919 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1920 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1921 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1922 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1923 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1924 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1925 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1926 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1927 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1928 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1929 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1930 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1931 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1932 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1933 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1934 * ccomps * dcomps);

            auto g_z_0_xyyyz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1935 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1944 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1945 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1946 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1947 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1948 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1949 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1950 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1951 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1952 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1953 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1954 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1955 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1956 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1957 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1958 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1959 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1960 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1961 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1962 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1963 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 1964 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 1965 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 1966 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 1967 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 1968 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 1969 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 1970 * ccomps * dcomps);

            auto g_z_0_xyyzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 1971 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 1980 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 1981 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 1982 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 1983 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 1984 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 1985 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 1986 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 1987 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 1988 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 1989 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 1990 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 1991 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 1992 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 1993 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 1994 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 1995 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 1996 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 1997 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 1998 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 1999 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 2000 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 2001 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 2002 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 2003 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 2004 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 2005 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 2006 * ccomps * dcomps);

            auto g_z_0_xyzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 2007 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 2016 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 2017 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 2018 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 2019 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 2020 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 2021 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 2022 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 2023 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 2024 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 2025 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 2026 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 2027 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 2028 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 2029 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 2030 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 2031 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 2032 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 2033 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 2034 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 2035 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 2036 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 2037 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 2038 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 2039 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 2040 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 2041 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 2042 * ccomps * dcomps);

            auto g_z_0_xzzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 2043 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxxx = cbuffer.data(hk_geom_10_off + 2052 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxxy = cbuffer.data(hk_geom_10_off + 2053 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxxz = cbuffer.data(hk_geom_10_off + 2054 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxyy = cbuffer.data(hk_geom_10_off + 2055 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxyz = cbuffer.data(hk_geom_10_off + 2056 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxxzz = cbuffer.data(hk_geom_10_off + 2057 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxyyy = cbuffer.data(hk_geom_10_off + 2058 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxyyz = cbuffer.data(hk_geom_10_off + 2059 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxyzz = cbuffer.data(hk_geom_10_off + 2060 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxzzz = cbuffer.data(hk_geom_10_off + 2061 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyyyy = cbuffer.data(hk_geom_10_off + 2062 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyyyz = cbuffer.data(hk_geom_10_off + 2063 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyyzz = cbuffer.data(hk_geom_10_off + 2064 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyzzz = cbuffer.data(hk_geom_10_off + 2065 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxzzzz = cbuffer.data(hk_geom_10_off + 2066 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyyyy = cbuffer.data(hk_geom_10_off + 2067 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyyyz = cbuffer.data(hk_geom_10_off + 2068 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyyzz = cbuffer.data(hk_geom_10_off + 2069 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyzzz = cbuffer.data(hk_geom_10_off + 2070 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyzzzz = cbuffer.data(hk_geom_10_off + 2071 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxzzzzz = cbuffer.data(hk_geom_10_off + 2072 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyyyy = cbuffer.data(hk_geom_10_off + 2073 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyyyz = cbuffer.data(hk_geom_10_off + 2074 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyyzz = cbuffer.data(hk_geom_10_off + 2075 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyzzz = cbuffer.data(hk_geom_10_off + 2076 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyzzzz = cbuffer.data(hk_geom_10_off + 2077 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyzzzzz = cbuffer.data(hk_geom_10_off + 2078 * ccomps * dcomps);

            auto g_z_0_yyyyy_xzzzzzz = cbuffer.data(hk_geom_10_off + 2079 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyyyy = cbuffer.data(hk_geom_10_off + 2080 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyyyz = cbuffer.data(hk_geom_10_off + 2081 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyyzz = cbuffer.data(hk_geom_10_off + 2082 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyzzz = cbuffer.data(hk_geom_10_off + 2083 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyzzzz = cbuffer.data(hk_geom_10_off + 2084 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyzzzzz = cbuffer.data(hk_geom_10_off + 2085 * ccomps * dcomps);

            auto g_z_0_yyyyy_yzzzzzz = cbuffer.data(hk_geom_10_off + 2086 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxxx = cbuffer.data(hk_geom_10_off + 2088 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxxy = cbuffer.data(hk_geom_10_off + 2089 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxxz = cbuffer.data(hk_geom_10_off + 2090 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxyy = cbuffer.data(hk_geom_10_off + 2091 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxyz = cbuffer.data(hk_geom_10_off + 2092 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxxzz = cbuffer.data(hk_geom_10_off + 2093 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxyyy = cbuffer.data(hk_geom_10_off + 2094 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxyyz = cbuffer.data(hk_geom_10_off + 2095 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxyzz = cbuffer.data(hk_geom_10_off + 2096 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxzzz = cbuffer.data(hk_geom_10_off + 2097 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyyyy = cbuffer.data(hk_geom_10_off + 2098 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyyyz = cbuffer.data(hk_geom_10_off + 2099 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyyzz = cbuffer.data(hk_geom_10_off + 2100 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyzzz = cbuffer.data(hk_geom_10_off + 2101 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxzzzz = cbuffer.data(hk_geom_10_off + 2102 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyyyy = cbuffer.data(hk_geom_10_off + 2103 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyyyz = cbuffer.data(hk_geom_10_off + 2104 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyyzz = cbuffer.data(hk_geom_10_off + 2105 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyzzz = cbuffer.data(hk_geom_10_off + 2106 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyzzzz = cbuffer.data(hk_geom_10_off + 2107 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxzzzzz = cbuffer.data(hk_geom_10_off + 2108 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyyyy = cbuffer.data(hk_geom_10_off + 2109 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyyyz = cbuffer.data(hk_geom_10_off + 2110 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyyzz = cbuffer.data(hk_geom_10_off + 2111 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyzzz = cbuffer.data(hk_geom_10_off + 2112 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyzzzz = cbuffer.data(hk_geom_10_off + 2113 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyzzzzz = cbuffer.data(hk_geom_10_off + 2114 * ccomps * dcomps);

            auto g_z_0_yyyyz_xzzzzzz = cbuffer.data(hk_geom_10_off + 2115 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyyyy = cbuffer.data(hk_geom_10_off + 2116 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyyyz = cbuffer.data(hk_geom_10_off + 2117 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyyzz = cbuffer.data(hk_geom_10_off + 2118 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyzzz = cbuffer.data(hk_geom_10_off + 2119 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyzzzz = cbuffer.data(hk_geom_10_off + 2120 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyzzzzz = cbuffer.data(hk_geom_10_off + 2121 * ccomps * dcomps);

            auto g_z_0_yyyyz_yzzzzzz = cbuffer.data(hk_geom_10_off + 2122 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 2124 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 2125 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 2126 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 2127 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 2128 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 2129 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 2130 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 2131 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 2132 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 2133 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 2134 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 2135 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 2136 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 2137 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 2138 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 2139 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 2140 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 2141 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 2142 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 2143 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 2144 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 2145 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 2146 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 2147 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 2148 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 2149 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 2150 * ccomps * dcomps);

            auto g_z_0_yyyzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 2151 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 2152 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 2153 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 2154 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 2155 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 2156 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 2157 * ccomps * dcomps);

            auto g_z_0_yyyzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 2158 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 2160 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 2161 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 2162 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 2163 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 2164 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 2165 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 2166 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 2167 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 2168 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 2169 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 2170 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 2171 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 2172 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 2173 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 2174 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 2175 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 2176 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 2177 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 2178 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 2179 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 2180 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 2181 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 2182 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 2183 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 2184 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 2185 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 2186 * ccomps * dcomps);

            auto g_z_0_yyzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 2187 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 2188 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 2189 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 2190 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 2191 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 2192 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 2193 * ccomps * dcomps);

            auto g_z_0_yyzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 2194 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 2196 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 2197 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 2198 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 2199 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 2200 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 2201 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 2202 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 2203 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 2204 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 2205 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 2206 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 2207 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 2208 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 2209 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 2210 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 2211 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 2212 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 2213 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 2214 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 2215 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 2216 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 2217 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 2218 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 2219 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 2220 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 2221 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 2222 * ccomps * dcomps);

            auto g_z_0_yzzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 2223 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 2224 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 2225 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 2226 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 2227 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 2228 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 2229 * ccomps * dcomps);

            auto g_z_0_yzzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 2230 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxxx = cbuffer.data(hk_geom_10_off + 2232 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxxy = cbuffer.data(hk_geom_10_off + 2233 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxxz = cbuffer.data(hk_geom_10_off + 2234 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxyy = cbuffer.data(hk_geom_10_off + 2235 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxyz = cbuffer.data(hk_geom_10_off + 2236 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxxzz = cbuffer.data(hk_geom_10_off + 2237 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxyyy = cbuffer.data(hk_geom_10_off + 2238 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxyyz = cbuffer.data(hk_geom_10_off + 2239 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxyzz = cbuffer.data(hk_geom_10_off + 2240 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxzzz = cbuffer.data(hk_geom_10_off + 2241 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyyyy = cbuffer.data(hk_geom_10_off + 2242 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyyyz = cbuffer.data(hk_geom_10_off + 2243 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyyzz = cbuffer.data(hk_geom_10_off + 2244 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyzzz = cbuffer.data(hk_geom_10_off + 2245 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxzzzz = cbuffer.data(hk_geom_10_off + 2246 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyyyy = cbuffer.data(hk_geom_10_off + 2247 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyyyz = cbuffer.data(hk_geom_10_off + 2248 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyyzz = cbuffer.data(hk_geom_10_off + 2249 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyzzz = cbuffer.data(hk_geom_10_off + 2250 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyzzzz = cbuffer.data(hk_geom_10_off + 2251 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxzzzzz = cbuffer.data(hk_geom_10_off + 2252 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyyyy = cbuffer.data(hk_geom_10_off + 2253 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyyyz = cbuffer.data(hk_geom_10_off + 2254 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyyzz = cbuffer.data(hk_geom_10_off + 2255 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyzzz = cbuffer.data(hk_geom_10_off + 2256 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyzzzz = cbuffer.data(hk_geom_10_off + 2257 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyzzzzz = cbuffer.data(hk_geom_10_off + 2258 * ccomps * dcomps);

            auto g_z_0_zzzzz_xzzzzzz = cbuffer.data(hk_geom_10_off + 2259 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyyyy = cbuffer.data(hk_geom_10_off + 2260 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyyyz = cbuffer.data(hk_geom_10_off + 2261 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyyzz = cbuffer.data(hk_geom_10_off + 2262 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyzzz = cbuffer.data(hk_geom_10_off + 2263 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyzzzz = cbuffer.data(hk_geom_10_off + 2264 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyzzzzz = cbuffer.data(hk_geom_10_off + 2265 * ccomps * dcomps);

            auto g_z_0_zzzzz_yzzzzzz = cbuffer.data(hk_geom_10_off + 2266 * ccomps * dcomps);

            auto g_z_0_zzzzz_zzzzzzz = cbuffer.data(hk_geom_10_off + 2267 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_iixx

            const auto ii_geom_10_off = idx_geom_10_iixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_xxxxxx = cbuffer.data(ii_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxxxy = cbuffer.data(ii_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxxxz = cbuffer.data(ii_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxxyy = cbuffer.data(ii_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxxyz = cbuffer.data(ii_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxxzz = cbuffer.data(ii_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxyyy = cbuffer.data(ii_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxyyz = cbuffer.data(ii_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxyzz = cbuffer.data(ii_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxzzz = cbuffer.data(ii_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxyyyy = cbuffer.data(ii_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxyyyz = cbuffer.data(ii_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxyyzz = cbuffer.data(ii_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxyzzz = cbuffer.data(ii_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxzzzz = cbuffer.data(ii_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyyyyy = cbuffer.data(ii_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyyyyz = cbuffer.data(ii_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyyyzz = cbuffer.data(ii_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyyzzz = cbuffer.data(ii_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyzzzz = cbuffer.data(ii_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xzzzzz = cbuffer.data(ii_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyyyyy = cbuffer.data(ii_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyyyyz = cbuffer.data(ii_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyyyzz = cbuffer.data(ii_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyyzzz = cbuffer.data(ii_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyzzzz = cbuffer.data(ii_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yzzzzz = cbuffer.data(ii_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxxx_zzzzzz = cbuffer.data(ii_geom_10_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxxxxx, g_x_0_xxxxx_xxxxxxx, g_x_0_xxxxx_xxxxxxy, g_x_0_xxxxx_xxxxxxz, g_x_0_xxxxx_xxxxxy, g_x_0_xxxxx_xxxxxyy, g_x_0_xxxxx_xxxxxyz, g_x_0_xxxxx_xxxxxz, g_x_0_xxxxx_xxxxxzz, g_x_0_xxxxx_xxxxyy, g_x_0_xxxxx_xxxxyyy, g_x_0_xxxxx_xxxxyyz, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxyzz, g_x_0_xxxxx_xxxxzz, g_x_0_xxxxx_xxxxzzz, g_x_0_xxxxx_xxxyyy, g_x_0_xxxxx_xxxyyyy, g_x_0_xxxxx_xxxyyyz, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyyzz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxyzzz, g_x_0_xxxxx_xxxzzz, g_x_0_xxxxx_xxxzzzz, g_x_0_xxxxx_xxyyyy, g_x_0_xxxxx_xxyyyyy, g_x_0_xxxxx_xxyyyyz, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyyzz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyyzzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxyzzzz, g_x_0_xxxxx_xxzzzz, g_x_0_xxxxx_xxzzzzz, g_x_0_xxxxx_xyyyyy, g_x_0_xxxxx_xyyyyyy, g_x_0_xxxxx_xyyyyyz, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyyzz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyyzzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyyzzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xyzzzzz, g_x_0_xxxxx_xzzzzz, g_x_0_xxxxx_xzzzzzz, g_x_0_xxxxx_yyyyyy, g_x_0_xxxxx_yyyyyz, g_x_0_xxxxx_yyyyzz, g_x_0_xxxxx_yyyzzz, g_x_0_xxxxx_yyzzzz, g_x_0_xxxxx_yzzzzz, g_x_0_xxxxx_zzzzzz, g_x_0_xxxxxx_xxxxxx, g_x_0_xxxxxx_xxxxxy, g_x_0_xxxxxx_xxxxxz, g_x_0_xxxxxx_xxxxyy, g_x_0_xxxxxx_xxxxyz, g_x_0_xxxxxx_xxxxzz, g_x_0_xxxxxx_xxxyyy, g_x_0_xxxxxx_xxxyyz, g_x_0_xxxxxx_xxxyzz, g_x_0_xxxxxx_xxxzzz, g_x_0_xxxxxx_xxyyyy, g_x_0_xxxxxx_xxyyyz, g_x_0_xxxxxx_xxyyzz, g_x_0_xxxxxx_xxyzzz, g_x_0_xxxxxx_xxzzzz, g_x_0_xxxxxx_xyyyyy, g_x_0_xxxxxx_xyyyyz, g_x_0_xxxxxx_xyyyzz, g_x_0_xxxxxx_xyyzzz, g_x_0_xxxxxx_xyzzzz, g_x_0_xxxxxx_xzzzzz, g_x_0_xxxxxx_yyyyyy, g_x_0_xxxxxx_yyyyyz, g_x_0_xxxxxx_yyyyzz, g_x_0_xxxxxx_yyyzzz, g_x_0_xxxxxx_yyzzzz, g_x_0_xxxxxx_yzzzzz, g_x_0_xxxxxx_zzzzzz, g_xxxxx_xxxxxx, g_xxxxx_xxxxxy, g_xxxxx_xxxxxz, g_xxxxx_xxxxyy, g_xxxxx_xxxxyz, g_xxxxx_xxxxzz, g_xxxxx_xxxyyy, g_xxxxx_xxxyyz, g_xxxxx_xxxyzz, g_xxxxx_xxxzzz, g_xxxxx_xxyyyy, g_xxxxx_xxyyyz, g_xxxxx_xxyyzz, g_xxxxx_xxyzzz, g_xxxxx_xxzzzz, g_xxxxx_xyyyyy, g_xxxxx_xyyyyz, g_xxxxx_xyyyzz, g_xxxxx_xyyzzz, g_xxxxx_xyzzzz, g_xxxxx_xzzzzz, g_xxxxx_yyyyyy, g_xxxxx_yyyyyz, g_xxxxx_yyyyzz, g_xxxxx_yyyzzz, g_xxxxx_yyzzzz, g_xxxxx_yzzzzz, g_xxxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_xxxxxx[k] = -g_xxxxx_xxxxxx[k] - g_x_0_xxxxx_xxxxxx[k] * ab_x + g_x_0_xxxxx_xxxxxxx[k];

                g_x_0_xxxxxx_xxxxxy[k] = -g_xxxxx_xxxxxy[k] - g_x_0_xxxxx_xxxxxy[k] * ab_x + g_x_0_xxxxx_xxxxxxy[k];

                g_x_0_xxxxxx_xxxxxz[k] = -g_xxxxx_xxxxxz[k] - g_x_0_xxxxx_xxxxxz[k] * ab_x + g_x_0_xxxxx_xxxxxxz[k];

                g_x_0_xxxxxx_xxxxyy[k] = -g_xxxxx_xxxxyy[k] - g_x_0_xxxxx_xxxxyy[k] * ab_x + g_x_0_xxxxx_xxxxxyy[k];

                g_x_0_xxxxxx_xxxxyz[k] = -g_xxxxx_xxxxyz[k] - g_x_0_xxxxx_xxxxyz[k] * ab_x + g_x_0_xxxxx_xxxxxyz[k];

                g_x_0_xxxxxx_xxxxzz[k] = -g_xxxxx_xxxxzz[k] - g_x_0_xxxxx_xxxxzz[k] * ab_x + g_x_0_xxxxx_xxxxxzz[k];

                g_x_0_xxxxxx_xxxyyy[k] = -g_xxxxx_xxxyyy[k] - g_x_0_xxxxx_xxxyyy[k] * ab_x + g_x_0_xxxxx_xxxxyyy[k];

                g_x_0_xxxxxx_xxxyyz[k] = -g_xxxxx_xxxyyz[k] - g_x_0_xxxxx_xxxyyz[k] * ab_x + g_x_0_xxxxx_xxxxyyz[k];

                g_x_0_xxxxxx_xxxyzz[k] = -g_xxxxx_xxxyzz[k] - g_x_0_xxxxx_xxxyzz[k] * ab_x + g_x_0_xxxxx_xxxxyzz[k];

                g_x_0_xxxxxx_xxxzzz[k] = -g_xxxxx_xxxzzz[k] - g_x_0_xxxxx_xxxzzz[k] * ab_x + g_x_0_xxxxx_xxxxzzz[k];

                g_x_0_xxxxxx_xxyyyy[k] = -g_xxxxx_xxyyyy[k] - g_x_0_xxxxx_xxyyyy[k] * ab_x + g_x_0_xxxxx_xxxyyyy[k];

                g_x_0_xxxxxx_xxyyyz[k] = -g_xxxxx_xxyyyz[k] - g_x_0_xxxxx_xxyyyz[k] * ab_x + g_x_0_xxxxx_xxxyyyz[k];

                g_x_0_xxxxxx_xxyyzz[k] = -g_xxxxx_xxyyzz[k] - g_x_0_xxxxx_xxyyzz[k] * ab_x + g_x_0_xxxxx_xxxyyzz[k];

                g_x_0_xxxxxx_xxyzzz[k] = -g_xxxxx_xxyzzz[k] - g_x_0_xxxxx_xxyzzz[k] * ab_x + g_x_0_xxxxx_xxxyzzz[k];

                g_x_0_xxxxxx_xxzzzz[k] = -g_xxxxx_xxzzzz[k] - g_x_0_xxxxx_xxzzzz[k] * ab_x + g_x_0_xxxxx_xxxzzzz[k];

                g_x_0_xxxxxx_xyyyyy[k] = -g_xxxxx_xyyyyy[k] - g_x_0_xxxxx_xyyyyy[k] * ab_x + g_x_0_xxxxx_xxyyyyy[k];

                g_x_0_xxxxxx_xyyyyz[k] = -g_xxxxx_xyyyyz[k] - g_x_0_xxxxx_xyyyyz[k] * ab_x + g_x_0_xxxxx_xxyyyyz[k];

                g_x_0_xxxxxx_xyyyzz[k] = -g_xxxxx_xyyyzz[k] - g_x_0_xxxxx_xyyyzz[k] * ab_x + g_x_0_xxxxx_xxyyyzz[k];

                g_x_0_xxxxxx_xyyzzz[k] = -g_xxxxx_xyyzzz[k] - g_x_0_xxxxx_xyyzzz[k] * ab_x + g_x_0_xxxxx_xxyyzzz[k];

                g_x_0_xxxxxx_xyzzzz[k] = -g_xxxxx_xyzzzz[k] - g_x_0_xxxxx_xyzzzz[k] * ab_x + g_x_0_xxxxx_xxyzzzz[k];

                g_x_0_xxxxxx_xzzzzz[k] = -g_xxxxx_xzzzzz[k] - g_x_0_xxxxx_xzzzzz[k] * ab_x + g_x_0_xxxxx_xxzzzzz[k];

                g_x_0_xxxxxx_yyyyyy[k] = -g_xxxxx_yyyyyy[k] - g_x_0_xxxxx_yyyyyy[k] * ab_x + g_x_0_xxxxx_xyyyyyy[k];

                g_x_0_xxxxxx_yyyyyz[k] = -g_xxxxx_yyyyyz[k] - g_x_0_xxxxx_yyyyyz[k] * ab_x + g_x_0_xxxxx_xyyyyyz[k];

                g_x_0_xxxxxx_yyyyzz[k] = -g_xxxxx_yyyyzz[k] - g_x_0_xxxxx_yyyyzz[k] * ab_x + g_x_0_xxxxx_xyyyyzz[k];

                g_x_0_xxxxxx_yyyzzz[k] = -g_xxxxx_yyyzzz[k] - g_x_0_xxxxx_yyyzzz[k] * ab_x + g_x_0_xxxxx_xyyyzzz[k];

                g_x_0_xxxxxx_yyzzzz[k] = -g_xxxxx_yyzzzz[k] - g_x_0_xxxxx_yyzzzz[k] * ab_x + g_x_0_xxxxx_xyyzzzz[k];

                g_x_0_xxxxxx_yzzzzz[k] = -g_xxxxx_yzzzzz[k] - g_x_0_xxxxx_yzzzzz[k] * ab_x + g_x_0_xxxxx_xyzzzzz[k];

                g_x_0_xxxxxx_zzzzzz[k] = -g_xxxxx_zzzzzz[k] - g_x_0_xxxxx_zzzzzz[k] * ab_x + g_x_0_xxxxx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_xxxxxx = cbuffer.data(ii_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxxxy = cbuffer.data(ii_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxxxz = cbuffer.data(ii_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxxyy = cbuffer.data(ii_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxxyz = cbuffer.data(ii_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxxzz = cbuffer.data(ii_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxyyy = cbuffer.data(ii_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxyyz = cbuffer.data(ii_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxyzz = cbuffer.data(ii_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxzzz = cbuffer.data(ii_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxyyyy = cbuffer.data(ii_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxyyyz = cbuffer.data(ii_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxyyzz = cbuffer.data(ii_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxyzzz = cbuffer.data(ii_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxzzzz = cbuffer.data(ii_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyyyyy = cbuffer.data(ii_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyyyyz = cbuffer.data(ii_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyyyzz = cbuffer.data(ii_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyyzzz = cbuffer.data(ii_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyzzzz = cbuffer.data(ii_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xzzzzz = cbuffer.data(ii_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyyyyy = cbuffer.data(ii_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyyyyz = cbuffer.data(ii_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyyyzz = cbuffer.data(ii_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyyzzz = cbuffer.data(ii_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyzzzz = cbuffer.data(ii_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yzzzzz = cbuffer.data(ii_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxxxy_zzzzzz = cbuffer.data(ii_geom_10_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxxxxx, g_x_0_xxxxx_xxxxxxy, g_x_0_xxxxx_xxxxxy, g_x_0_xxxxx_xxxxxyy, g_x_0_xxxxx_xxxxxyz, g_x_0_xxxxx_xxxxxz, g_x_0_xxxxx_xxxxyy, g_x_0_xxxxx_xxxxyyy, g_x_0_xxxxx_xxxxyyz, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxyzz, g_x_0_xxxxx_xxxxzz, g_x_0_xxxxx_xxxyyy, g_x_0_xxxxx_xxxyyyy, g_x_0_xxxxx_xxxyyyz, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyyzz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxyzzz, g_x_0_xxxxx_xxxzzz, g_x_0_xxxxx_xxyyyy, g_x_0_xxxxx_xxyyyyy, g_x_0_xxxxx_xxyyyyz, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyyzz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyyzzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxyzzzz, g_x_0_xxxxx_xxzzzz, g_x_0_xxxxx_xyyyyy, g_x_0_xxxxx_xyyyyyy, g_x_0_xxxxx_xyyyyyz, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyyzz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyyzzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyyzzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xyzzzzz, g_x_0_xxxxx_xzzzzz, g_x_0_xxxxx_yyyyyy, g_x_0_xxxxx_yyyyyyy, g_x_0_xxxxx_yyyyyyz, g_x_0_xxxxx_yyyyyz, g_x_0_xxxxx_yyyyyzz, g_x_0_xxxxx_yyyyzz, g_x_0_xxxxx_yyyyzzz, g_x_0_xxxxx_yyyzzz, g_x_0_xxxxx_yyyzzzz, g_x_0_xxxxx_yyzzzz, g_x_0_xxxxx_yyzzzzz, g_x_0_xxxxx_yzzzzz, g_x_0_xxxxx_yzzzzzz, g_x_0_xxxxx_zzzzzz, g_x_0_xxxxxy_xxxxxx, g_x_0_xxxxxy_xxxxxy, g_x_0_xxxxxy_xxxxxz, g_x_0_xxxxxy_xxxxyy, g_x_0_xxxxxy_xxxxyz, g_x_0_xxxxxy_xxxxzz, g_x_0_xxxxxy_xxxyyy, g_x_0_xxxxxy_xxxyyz, g_x_0_xxxxxy_xxxyzz, g_x_0_xxxxxy_xxxzzz, g_x_0_xxxxxy_xxyyyy, g_x_0_xxxxxy_xxyyyz, g_x_0_xxxxxy_xxyyzz, g_x_0_xxxxxy_xxyzzz, g_x_0_xxxxxy_xxzzzz, g_x_0_xxxxxy_xyyyyy, g_x_0_xxxxxy_xyyyyz, g_x_0_xxxxxy_xyyyzz, g_x_0_xxxxxy_xyyzzz, g_x_0_xxxxxy_xyzzzz, g_x_0_xxxxxy_xzzzzz, g_x_0_xxxxxy_yyyyyy, g_x_0_xxxxxy_yyyyyz, g_x_0_xxxxxy_yyyyzz, g_x_0_xxxxxy_yyyzzz, g_x_0_xxxxxy_yyzzzz, g_x_0_xxxxxy_yzzzzz, g_x_0_xxxxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_xxxxxx[k] = -g_x_0_xxxxx_xxxxxx[k] * ab_y + g_x_0_xxxxx_xxxxxxy[k];

                g_x_0_xxxxxy_xxxxxy[k] = -g_x_0_xxxxx_xxxxxy[k] * ab_y + g_x_0_xxxxx_xxxxxyy[k];

                g_x_0_xxxxxy_xxxxxz[k] = -g_x_0_xxxxx_xxxxxz[k] * ab_y + g_x_0_xxxxx_xxxxxyz[k];

                g_x_0_xxxxxy_xxxxyy[k] = -g_x_0_xxxxx_xxxxyy[k] * ab_y + g_x_0_xxxxx_xxxxyyy[k];

                g_x_0_xxxxxy_xxxxyz[k] = -g_x_0_xxxxx_xxxxyz[k] * ab_y + g_x_0_xxxxx_xxxxyyz[k];

                g_x_0_xxxxxy_xxxxzz[k] = -g_x_0_xxxxx_xxxxzz[k] * ab_y + g_x_0_xxxxx_xxxxyzz[k];

                g_x_0_xxxxxy_xxxyyy[k] = -g_x_0_xxxxx_xxxyyy[k] * ab_y + g_x_0_xxxxx_xxxyyyy[k];

                g_x_0_xxxxxy_xxxyyz[k] = -g_x_0_xxxxx_xxxyyz[k] * ab_y + g_x_0_xxxxx_xxxyyyz[k];

                g_x_0_xxxxxy_xxxyzz[k] = -g_x_0_xxxxx_xxxyzz[k] * ab_y + g_x_0_xxxxx_xxxyyzz[k];

                g_x_0_xxxxxy_xxxzzz[k] = -g_x_0_xxxxx_xxxzzz[k] * ab_y + g_x_0_xxxxx_xxxyzzz[k];

                g_x_0_xxxxxy_xxyyyy[k] = -g_x_0_xxxxx_xxyyyy[k] * ab_y + g_x_0_xxxxx_xxyyyyy[k];

                g_x_0_xxxxxy_xxyyyz[k] = -g_x_0_xxxxx_xxyyyz[k] * ab_y + g_x_0_xxxxx_xxyyyyz[k];

                g_x_0_xxxxxy_xxyyzz[k] = -g_x_0_xxxxx_xxyyzz[k] * ab_y + g_x_0_xxxxx_xxyyyzz[k];

                g_x_0_xxxxxy_xxyzzz[k] = -g_x_0_xxxxx_xxyzzz[k] * ab_y + g_x_0_xxxxx_xxyyzzz[k];

                g_x_0_xxxxxy_xxzzzz[k] = -g_x_0_xxxxx_xxzzzz[k] * ab_y + g_x_0_xxxxx_xxyzzzz[k];

                g_x_0_xxxxxy_xyyyyy[k] = -g_x_0_xxxxx_xyyyyy[k] * ab_y + g_x_0_xxxxx_xyyyyyy[k];

                g_x_0_xxxxxy_xyyyyz[k] = -g_x_0_xxxxx_xyyyyz[k] * ab_y + g_x_0_xxxxx_xyyyyyz[k];

                g_x_0_xxxxxy_xyyyzz[k] = -g_x_0_xxxxx_xyyyzz[k] * ab_y + g_x_0_xxxxx_xyyyyzz[k];

                g_x_0_xxxxxy_xyyzzz[k] = -g_x_0_xxxxx_xyyzzz[k] * ab_y + g_x_0_xxxxx_xyyyzzz[k];

                g_x_0_xxxxxy_xyzzzz[k] = -g_x_0_xxxxx_xyzzzz[k] * ab_y + g_x_0_xxxxx_xyyzzzz[k];

                g_x_0_xxxxxy_xzzzzz[k] = -g_x_0_xxxxx_xzzzzz[k] * ab_y + g_x_0_xxxxx_xyzzzzz[k];

                g_x_0_xxxxxy_yyyyyy[k] = -g_x_0_xxxxx_yyyyyy[k] * ab_y + g_x_0_xxxxx_yyyyyyy[k];

                g_x_0_xxxxxy_yyyyyz[k] = -g_x_0_xxxxx_yyyyyz[k] * ab_y + g_x_0_xxxxx_yyyyyyz[k];

                g_x_0_xxxxxy_yyyyzz[k] = -g_x_0_xxxxx_yyyyzz[k] * ab_y + g_x_0_xxxxx_yyyyyzz[k];

                g_x_0_xxxxxy_yyyzzz[k] = -g_x_0_xxxxx_yyyzzz[k] * ab_y + g_x_0_xxxxx_yyyyzzz[k];

                g_x_0_xxxxxy_yyzzzz[k] = -g_x_0_xxxxx_yyzzzz[k] * ab_y + g_x_0_xxxxx_yyyzzzz[k];

                g_x_0_xxxxxy_yzzzzz[k] = -g_x_0_xxxxx_yzzzzz[k] * ab_y + g_x_0_xxxxx_yyzzzzz[k];

                g_x_0_xxxxxy_zzzzzz[k] = -g_x_0_xxxxx_zzzzzz[k] * ab_y + g_x_0_xxxxx_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_xxxxxx = cbuffer.data(ii_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxxxy = cbuffer.data(ii_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxxxz = cbuffer.data(ii_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxxyy = cbuffer.data(ii_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxxyz = cbuffer.data(ii_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxxzz = cbuffer.data(ii_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxyyy = cbuffer.data(ii_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxyyz = cbuffer.data(ii_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxyzz = cbuffer.data(ii_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxzzz = cbuffer.data(ii_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxyyyy = cbuffer.data(ii_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxyyyz = cbuffer.data(ii_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxyyzz = cbuffer.data(ii_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxyzzz = cbuffer.data(ii_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxzzzz = cbuffer.data(ii_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyyyyy = cbuffer.data(ii_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyyyyz = cbuffer.data(ii_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyyyzz = cbuffer.data(ii_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyyzzz = cbuffer.data(ii_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyzzzz = cbuffer.data(ii_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xzzzzz = cbuffer.data(ii_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyyyyy = cbuffer.data(ii_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyyyyz = cbuffer.data(ii_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyyyzz = cbuffer.data(ii_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyyzzz = cbuffer.data(ii_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyzzzz = cbuffer.data(ii_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yzzzzz = cbuffer.data(ii_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxxxz_zzzzzz = cbuffer.data(ii_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxxxxx, g_x_0_xxxxx_xxxxxxz, g_x_0_xxxxx_xxxxxy, g_x_0_xxxxx_xxxxxyz, g_x_0_xxxxx_xxxxxz, g_x_0_xxxxx_xxxxxzz, g_x_0_xxxxx_xxxxyy, g_x_0_xxxxx_xxxxyyz, g_x_0_xxxxx_xxxxyz, g_x_0_xxxxx_xxxxyzz, g_x_0_xxxxx_xxxxzz, g_x_0_xxxxx_xxxxzzz, g_x_0_xxxxx_xxxyyy, g_x_0_xxxxx_xxxyyyz, g_x_0_xxxxx_xxxyyz, g_x_0_xxxxx_xxxyyzz, g_x_0_xxxxx_xxxyzz, g_x_0_xxxxx_xxxyzzz, g_x_0_xxxxx_xxxzzz, g_x_0_xxxxx_xxxzzzz, g_x_0_xxxxx_xxyyyy, g_x_0_xxxxx_xxyyyyz, g_x_0_xxxxx_xxyyyz, g_x_0_xxxxx_xxyyyzz, g_x_0_xxxxx_xxyyzz, g_x_0_xxxxx_xxyyzzz, g_x_0_xxxxx_xxyzzz, g_x_0_xxxxx_xxyzzzz, g_x_0_xxxxx_xxzzzz, g_x_0_xxxxx_xxzzzzz, g_x_0_xxxxx_xyyyyy, g_x_0_xxxxx_xyyyyyz, g_x_0_xxxxx_xyyyyz, g_x_0_xxxxx_xyyyyzz, g_x_0_xxxxx_xyyyzz, g_x_0_xxxxx_xyyyzzz, g_x_0_xxxxx_xyyzzz, g_x_0_xxxxx_xyyzzzz, g_x_0_xxxxx_xyzzzz, g_x_0_xxxxx_xyzzzzz, g_x_0_xxxxx_xzzzzz, g_x_0_xxxxx_xzzzzzz, g_x_0_xxxxx_yyyyyy, g_x_0_xxxxx_yyyyyyz, g_x_0_xxxxx_yyyyyz, g_x_0_xxxxx_yyyyyzz, g_x_0_xxxxx_yyyyzz, g_x_0_xxxxx_yyyyzzz, g_x_0_xxxxx_yyyzzz, g_x_0_xxxxx_yyyzzzz, g_x_0_xxxxx_yyzzzz, g_x_0_xxxxx_yyzzzzz, g_x_0_xxxxx_yzzzzz, g_x_0_xxxxx_yzzzzzz, g_x_0_xxxxx_zzzzzz, g_x_0_xxxxx_zzzzzzz, g_x_0_xxxxxz_xxxxxx, g_x_0_xxxxxz_xxxxxy, g_x_0_xxxxxz_xxxxxz, g_x_0_xxxxxz_xxxxyy, g_x_0_xxxxxz_xxxxyz, g_x_0_xxxxxz_xxxxzz, g_x_0_xxxxxz_xxxyyy, g_x_0_xxxxxz_xxxyyz, g_x_0_xxxxxz_xxxyzz, g_x_0_xxxxxz_xxxzzz, g_x_0_xxxxxz_xxyyyy, g_x_0_xxxxxz_xxyyyz, g_x_0_xxxxxz_xxyyzz, g_x_0_xxxxxz_xxyzzz, g_x_0_xxxxxz_xxzzzz, g_x_0_xxxxxz_xyyyyy, g_x_0_xxxxxz_xyyyyz, g_x_0_xxxxxz_xyyyzz, g_x_0_xxxxxz_xyyzzz, g_x_0_xxxxxz_xyzzzz, g_x_0_xxxxxz_xzzzzz, g_x_0_xxxxxz_yyyyyy, g_x_0_xxxxxz_yyyyyz, g_x_0_xxxxxz_yyyyzz, g_x_0_xxxxxz_yyyzzz, g_x_0_xxxxxz_yyzzzz, g_x_0_xxxxxz_yzzzzz, g_x_0_xxxxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_xxxxxx[k] = -g_x_0_xxxxx_xxxxxx[k] * ab_z + g_x_0_xxxxx_xxxxxxz[k];

                g_x_0_xxxxxz_xxxxxy[k] = -g_x_0_xxxxx_xxxxxy[k] * ab_z + g_x_0_xxxxx_xxxxxyz[k];

                g_x_0_xxxxxz_xxxxxz[k] = -g_x_0_xxxxx_xxxxxz[k] * ab_z + g_x_0_xxxxx_xxxxxzz[k];

                g_x_0_xxxxxz_xxxxyy[k] = -g_x_0_xxxxx_xxxxyy[k] * ab_z + g_x_0_xxxxx_xxxxyyz[k];

                g_x_0_xxxxxz_xxxxyz[k] = -g_x_0_xxxxx_xxxxyz[k] * ab_z + g_x_0_xxxxx_xxxxyzz[k];

                g_x_0_xxxxxz_xxxxzz[k] = -g_x_0_xxxxx_xxxxzz[k] * ab_z + g_x_0_xxxxx_xxxxzzz[k];

                g_x_0_xxxxxz_xxxyyy[k] = -g_x_0_xxxxx_xxxyyy[k] * ab_z + g_x_0_xxxxx_xxxyyyz[k];

                g_x_0_xxxxxz_xxxyyz[k] = -g_x_0_xxxxx_xxxyyz[k] * ab_z + g_x_0_xxxxx_xxxyyzz[k];

                g_x_0_xxxxxz_xxxyzz[k] = -g_x_0_xxxxx_xxxyzz[k] * ab_z + g_x_0_xxxxx_xxxyzzz[k];

                g_x_0_xxxxxz_xxxzzz[k] = -g_x_0_xxxxx_xxxzzz[k] * ab_z + g_x_0_xxxxx_xxxzzzz[k];

                g_x_0_xxxxxz_xxyyyy[k] = -g_x_0_xxxxx_xxyyyy[k] * ab_z + g_x_0_xxxxx_xxyyyyz[k];

                g_x_0_xxxxxz_xxyyyz[k] = -g_x_0_xxxxx_xxyyyz[k] * ab_z + g_x_0_xxxxx_xxyyyzz[k];

                g_x_0_xxxxxz_xxyyzz[k] = -g_x_0_xxxxx_xxyyzz[k] * ab_z + g_x_0_xxxxx_xxyyzzz[k];

                g_x_0_xxxxxz_xxyzzz[k] = -g_x_0_xxxxx_xxyzzz[k] * ab_z + g_x_0_xxxxx_xxyzzzz[k];

                g_x_0_xxxxxz_xxzzzz[k] = -g_x_0_xxxxx_xxzzzz[k] * ab_z + g_x_0_xxxxx_xxzzzzz[k];

                g_x_0_xxxxxz_xyyyyy[k] = -g_x_0_xxxxx_xyyyyy[k] * ab_z + g_x_0_xxxxx_xyyyyyz[k];

                g_x_0_xxxxxz_xyyyyz[k] = -g_x_0_xxxxx_xyyyyz[k] * ab_z + g_x_0_xxxxx_xyyyyzz[k];

                g_x_0_xxxxxz_xyyyzz[k] = -g_x_0_xxxxx_xyyyzz[k] * ab_z + g_x_0_xxxxx_xyyyzzz[k];

                g_x_0_xxxxxz_xyyzzz[k] = -g_x_0_xxxxx_xyyzzz[k] * ab_z + g_x_0_xxxxx_xyyzzzz[k];

                g_x_0_xxxxxz_xyzzzz[k] = -g_x_0_xxxxx_xyzzzz[k] * ab_z + g_x_0_xxxxx_xyzzzzz[k];

                g_x_0_xxxxxz_xzzzzz[k] = -g_x_0_xxxxx_xzzzzz[k] * ab_z + g_x_0_xxxxx_xzzzzzz[k];

                g_x_0_xxxxxz_yyyyyy[k] = -g_x_0_xxxxx_yyyyyy[k] * ab_z + g_x_0_xxxxx_yyyyyyz[k];

                g_x_0_xxxxxz_yyyyyz[k] = -g_x_0_xxxxx_yyyyyz[k] * ab_z + g_x_0_xxxxx_yyyyyzz[k];

                g_x_0_xxxxxz_yyyyzz[k] = -g_x_0_xxxxx_yyyyzz[k] * ab_z + g_x_0_xxxxx_yyyyzzz[k];

                g_x_0_xxxxxz_yyyzzz[k] = -g_x_0_xxxxx_yyyzzz[k] * ab_z + g_x_0_xxxxx_yyyzzzz[k];

                g_x_0_xxxxxz_yyzzzz[k] = -g_x_0_xxxxx_yyzzzz[k] * ab_z + g_x_0_xxxxx_yyzzzzz[k];

                g_x_0_xxxxxz_yzzzzz[k] = -g_x_0_xxxxx_yzzzzz[k] * ab_z + g_x_0_xxxxx_yzzzzzz[k];

                g_x_0_xxxxxz_zzzzzz[k] = -g_x_0_xxxxx_zzzzzz[k] * ab_z + g_x_0_xxxxx_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_xxxxxx = cbuffer.data(ii_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxxxy = cbuffer.data(ii_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxxxz = cbuffer.data(ii_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxxyy = cbuffer.data(ii_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxxyz = cbuffer.data(ii_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxxzz = cbuffer.data(ii_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxyyy = cbuffer.data(ii_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxyyz = cbuffer.data(ii_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxyzz = cbuffer.data(ii_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxzzz = cbuffer.data(ii_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxyyyy = cbuffer.data(ii_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxyyyz = cbuffer.data(ii_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxyyzz = cbuffer.data(ii_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxyzzz = cbuffer.data(ii_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxzzzz = cbuffer.data(ii_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyyyyy = cbuffer.data(ii_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyyyyz = cbuffer.data(ii_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyyyzz = cbuffer.data(ii_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyyzzz = cbuffer.data(ii_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyzzzz = cbuffer.data(ii_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xzzzzz = cbuffer.data(ii_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyyyyy = cbuffer.data(ii_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyyyyz = cbuffer.data(ii_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyyyzz = cbuffer.data(ii_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyyzzz = cbuffer.data(ii_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyzzzz = cbuffer.data(ii_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yzzzzz = cbuffer.data(ii_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxxxyy_zzzzzz = cbuffer.data(ii_geom_10_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxy_xxxxxx, g_x_0_xxxxy_xxxxxxy, g_x_0_xxxxy_xxxxxy, g_x_0_xxxxy_xxxxxyy, g_x_0_xxxxy_xxxxxyz, g_x_0_xxxxy_xxxxxz, g_x_0_xxxxy_xxxxyy, g_x_0_xxxxy_xxxxyyy, g_x_0_xxxxy_xxxxyyz, g_x_0_xxxxy_xxxxyz, g_x_0_xxxxy_xxxxyzz, g_x_0_xxxxy_xxxxzz, g_x_0_xxxxy_xxxyyy, g_x_0_xxxxy_xxxyyyy, g_x_0_xxxxy_xxxyyyz, g_x_0_xxxxy_xxxyyz, g_x_0_xxxxy_xxxyyzz, g_x_0_xxxxy_xxxyzz, g_x_0_xxxxy_xxxyzzz, g_x_0_xxxxy_xxxzzz, g_x_0_xxxxy_xxyyyy, g_x_0_xxxxy_xxyyyyy, g_x_0_xxxxy_xxyyyyz, g_x_0_xxxxy_xxyyyz, g_x_0_xxxxy_xxyyyzz, g_x_0_xxxxy_xxyyzz, g_x_0_xxxxy_xxyyzzz, g_x_0_xxxxy_xxyzzz, g_x_0_xxxxy_xxyzzzz, g_x_0_xxxxy_xxzzzz, g_x_0_xxxxy_xyyyyy, g_x_0_xxxxy_xyyyyyy, g_x_0_xxxxy_xyyyyyz, g_x_0_xxxxy_xyyyyz, g_x_0_xxxxy_xyyyyzz, g_x_0_xxxxy_xyyyzz, g_x_0_xxxxy_xyyyzzz, g_x_0_xxxxy_xyyzzz, g_x_0_xxxxy_xyyzzzz, g_x_0_xxxxy_xyzzzz, g_x_0_xxxxy_xyzzzzz, g_x_0_xxxxy_xzzzzz, g_x_0_xxxxy_yyyyyy, g_x_0_xxxxy_yyyyyyy, g_x_0_xxxxy_yyyyyyz, g_x_0_xxxxy_yyyyyz, g_x_0_xxxxy_yyyyyzz, g_x_0_xxxxy_yyyyzz, g_x_0_xxxxy_yyyyzzz, g_x_0_xxxxy_yyyzzz, g_x_0_xxxxy_yyyzzzz, g_x_0_xxxxy_yyzzzz, g_x_0_xxxxy_yyzzzzz, g_x_0_xxxxy_yzzzzz, g_x_0_xxxxy_yzzzzzz, g_x_0_xxxxy_zzzzzz, g_x_0_xxxxyy_xxxxxx, g_x_0_xxxxyy_xxxxxy, g_x_0_xxxxyy_xxxxxz, g_x_0_xxxxyy_xxxxyy, g_x_0_xxxxyy_xxxxyz, g_x_0_xxxxyy_xxxxzz, g_x_0_xxxxyy_xxxyyy, g_x_0_xxxxyy_xxxyyz, g_x_0_xxxxyy_xxxyzz, g_x_0_xxxxyy_xxxzzz, g_x_0_xxxxyy_xxyyyy, g_x_0_xxxxyy_xxyyyz, g_x_0_xxxxyy_xxyyzz, g_x_0_xxxxyy_xxyzzz, g_x_0_xxxxyy_xxzzzz, g_x_0_xxxxyy_xyyyyy, g_x_0_xxxxyy_xyyyyz, g_x_0_xxxxyy_xyyyzz, g_x_0_xxxxyy_xyyzzz, g_x_0_xxxxyy_xyzzzz, g_x_0_xxxxyy_xzzzzz, g_x_0_xxxxyy_yyyyyy, g_x_0_xxxxyy_yyyyyz, g_x_0_xxxxyy_yyyyzz, g_x_0_xxxxyy_yyyzzz, g_x_0_xxxxyy_yyzzzz, g_x_0_xxxxyy_yzzzzz, g_x_0_xxxxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_xxxxxx[k] = -g_x_0_xxxxy_xxxxxx[k] * ab_y + g_x_0_xxxxy_xxxxxxy[k];

                g_x_0_xxxxyy_xxxxxy[k] = -g_x_0_xxxxy_xxxxxy[k] * ab_y + g_x_0_xxxxy_xxxxxyy[k];

                g_x_0_xxxxyy_xxxxxz[k] = -g_x_0_xxxxy_xxxxxz[k] * ab_y + g_x_0_xxxxy_xxxxxyz[k];

                g_x_0_xxxxyy_xxxxyy[k] = -g_x_0_xxxxy_xxxxyy[k] * ab_y + g_x_0_xxxxy_xxxxyyy[k];

                g_x_0_xxxxyy_xxxxyz[k] = -g_x_0_xxxxy_xxxxyz[k] * ab_y + g_x_0_xxxxy_xxxxyyz[k];

                g_x_0_xxxxyy_xxxxzz[k] = -g_x_0_xxxxy_xxxxzz[k] * ab_y + g_x_0_xxxxy_xxxxyzz[k];

                g_x_0_xxxxyy_xxxyyy[k] = -g_x_0_xxxxy_xxxyyy[k] * ab_y + g_x_0_xxxxy_xxxyyyy[k];

                g_x_0_xxxxyy_xxxyyz[k] = -g_x_0_xxxxy_xxxyyz[k] * ab_y + g_x_0_xxxxy_xxxyyyz[k];

                g_x_0_xxxxyy_xxxyzz[k] = -g_x_0_xxxxy_xxxyzz[k] * ab_y + g_x_0_xxxxy_xxxyyzz[k];

                g_x_0_xxxxyy_xxxzzz[k] = -g_x_0_xxxxy_xxxzzz[k] * ab_y + g_x_0_xxxxy_xxxyzzz[k];

                g_x_0_xxxxyy_xxyyyy[k] = -g_x_0_xxxxy_xxyyyy[k] * ab_y + g_x_0_xxxxy_xxyyyyy[k];

                g_x_0_xxxxyy_xxyyyz[k] = -g_x_0_xxxxy_xxyyyz[k] * ab_y + g_x_0_xxxxy_xxyyyyz[k];

                g_x_0_xxxxyy_xxyyzz[k] = -g_x_0_xxxxy_xxyyzz[k] * ab_y + g_x_0_xxxxy_xxyyyzz[k];

                g_x_0_xxxxyy_xxyzzz[k] = -g_x_0_xxxxy_xxyzzz[k] * ab_y + g_x_0_xxxxy_xxyyzzz[k];

                g_x_0_xxxxyy_xxzzzz[k] = -g_x_0_xxxxy_xxzzzz[k] * ab_y + g_x_0_xxxxy_xxyzzzz[k];

                g_x_0_xxxxyy_xyyyyy[k] = -g_x_0_xxxxy_xyyyyy[k] * ab_y + g_x_0_xxxxy_xyyyyyy[k];

                g_x_0_xxxxyy_xyyyyz[k] = -g_x_0_xxxxy_xyyyyz[k] * ab_y + g_x_0_xxxxy_xyyyyyz[k];

                g_x_0_xxxxyy_xyyyzz[k] = -g_x_0_xxxxy_xyyyzz[k] * ab_y + g_x_0_xxxxy_xyyyyzz[k];

                g_x_0_xxxxyy_xyyzzz[k] = -g_x_0_xxxxy_xyyzzz[k] * ab_y + g_x_0_xxxxy_xyyyzzz[k];

                g_x_0_xxxxyy_xyzzzz[k] = -g_x_0_xxxxy_xyzzzz[k] * ab_y + g_x_0_xxxxy_xyyzzzz[k];

                g_x_0_xxxxyy_xzzzzz[k] = -g_x_0_xxxxy_xzzzzz[k] * ab_y + g_x_0_xxxxy_xyzzzzz[k];

                g_x_0_xxxxyy_yyyyyy[k] = -g_x_0_xxxxy_yyyyyy[k] * ab_y + g_x_0_xxxxy_yyyyyyy[k];

                g_x_0_xxxxyy_yyyyyz[k] = -g_x_0_xxxxy_yyyyyz[k] * ab_y + g_x_0_xxxxy_yyyyyyz[k];

                g_x_0_xxxxyy_yyyyzz[k] = -g_x_0_xxxxy_yyyyzz[k] * ab_y + g_x_0_xxxxy_yyyyyzz[k];

                g_x_0_xxxxyy_yyyzzz[k] = -g_x_0_xxxxy_yyyzzz[k] * ab_y + g_x_0_xxxxy_yyyyzzz[k];

                g_x_0_xxxxyy_yyzzzz[k] = -g_x_0_xxxxy_yyzzzz[k] * ab_y + g_x_0_xxxxy_yyyzzzz[k];

                g_x_0_xxxxyy_yzzzzz[k] = -g_x_0_xxxxy_yzzzzz[k] * ab_y + g_x_0_xxxxy_yyzzzzz[k];

                g_x_0_xxxxyy_zzzzzz[k] = -g_x_0_xxxxy_zzzzzz[k] * ab_y + g_x_0_xxxxy_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_xxxxxx = cbuffer.data(ii_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxxxy = cbuffer.data(ii_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxxxz = cbuffer.data(ii_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxxyy = cbuffer.data(ii_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxxyz = cbuffer.data(ii_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxxzz = cbuffer.data(ii_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxyyy = cbuffer.data(ii_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxyyz = cbuffer.data(ii_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxyzz = cbuffer.data(ii_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxzzz = cbuffer.data(ii_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxyyyy = cbuffer.data(ii_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxyyyz = cbuffer.data(ii_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxyyzz = cbuffer.data(ii_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxyzzz = cbuffer.data(ii_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxzzzz = cbuffer.data(ii_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyyyyy = cbuffer.data(ii_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyyyyz = cbuffer.data(ii_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyyyzz = cbuffer.data(ii_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyyzzz = cbuffer.data(ii_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyzzzz = cbuffer.data(ii_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xzzzzz = cbuffer.data(ii_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyyyyy = cbuffer.data(ii_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyyyyz = cbuffer.data(ii_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyyyzz = cbuffer.data(ii_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyyzzz = cbuffer.data(ii_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyzzzz = cbuffer.data(ii_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yzzzzz = cbuffer.data(ii_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxxxyz_zzzzzz = cbuffer.data(ii_geom_10_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxyz_xxxxxx, g_x_0_xxxxyz_xxxxxy, g_x_0_xxxxyz_xxxxxz, g_x_0_xxxxyz_xxxxyy, g_x_0_xxxxyz_xxxxyz, g_x_0_xxxxyz_xxxxzz, g_x_0_xxxxyz_xxxyyy, g_x_0_xxxxyz_xxxyyz, g_x_0_xxxxyz_xxxyzz, g_x_0_xxxxyz_xxxzzz, g_x_0_xxxxyz_xxyyyy, g_x_0_xxxxyz_xxyyyz, g_x_0_xxxxyz_xxyyzz, g_x_0_xxxxyz_xxyzzz, g_x_0_xxxxyz_xxzzzz, g_x_0_xxxxyz_xyyyyy, g_x_0_xxxxyz_xyyyyz, g_x_0_xxxxyz_xyyyzz, g_x_0_xxxxyz_xyyzzz, g_x_0_xxxxyz_xyzzzz, g_x_0_xxxxyz_xzzzzz, g_x_0_xxxxyz_yyyyyy, g_x_0_xxxxyz_yyyyyz, g_x_0_xxxxyz_yyyyzz, g_x_0_xxxxyz_yyyzzz, g_x_0_xxxxyz_yyzzzz, g_x_0_xxxxyz_yzzzzz, g_x_0_xxxxyz_zzzzzz, g_x_0_xxxxz_xxxxxx, g_x_0_xxxxz_xxxxxxy, g_x_0_xxxxz_xxxxxy, g_x_0_xxxxz_xxxxxyy, g_x_0_xxxxz_xxxxxyz, g_x_0_xxxxz_xxxxxz, g_x_0_xxxxz_xxxxyy, g_x_0_xxxxz_xxxxyyy, g_x_0_xxxxz_xxxxyyz, g_x_0_xxxxz_xxxxyz, g_x_0_xxxxz_xxxxyzz, g_x_0_xxxxz_xxxxzz, g_x_0_xxxxz_xxxyyy, g_x_0_xxxxz_xxxyyyy, g_x_0_xxxxz_xxxyyyz, g_x_0_xxxxz_xxxyyz, g_x_0_xxxxz_xxxyyzz, g_x_0_xxxxz_xxxyzz, g_x_0_xxxxz_xxxyzzz, g_x_0_xxxxz_xxxzzz, g_x_0_xxxxz_xxyyyy, g_x_0_xxxxz_xxyyyyy, g_x_0_xxxxz_xxyyyyz, g_x_0_xxxxz_xxyyyz, g_x_0_xxxxz_xxyyyzz, g_x_0_xxxxz_xxyyzz, g_x_0_xxxxz_xxyyzzz, g_x_0_xxxxz_xxyzzz, g_x_0_xxxxz_xxyzzzz, g_x_0_xxxxz_xxzzzz, g_x_0_xxxxz_xyyyyy, g_x_0_xxxxz_xyyyyyy, g_x_0_xxxxz_xyyyyyz, g_x_0_xxxxz_xyyyyz, g_x_0_xxxxz_xyyyyzz, g_x_0_xxxxz_xyyyzz, g_x_0_xxxxz_xyyyzzz, g_x_0_xxxxz_xyyzzz, g_x_0_xxxxz_xyyzzzz, g_x_0_xxxxz_xyzzzz, g_x_0_xxxxz_xyzzzzz, g_x_0_xxxxz_xzzzzz, g_x_0_xxxxz_yyyyyy, g_x_0_xxxxz_yyyyyyy, g_x_0_xxxxz_yyyyyyz, g_x_0_xxxxz_yyyyyz, g_x_0_xxxxz_yyyyyzz, g_x_0_xxxxz_yyyyzz, g_x_0_xxxxz_yyyyzzz, g_x_0_xxxxz_yyyzzz, g_x_0_xxxxz_yyyzzzz, g_x_0_xxxxz_yyzzzz, g_x_0_xxxxz_yyzzzzz, g_x_0_xxxxz_yzzzzz, g_x_0_xxxxz_yzzzzzz, g_x_0_xxxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_xxxxxx[k] = -g_x_0_xxxxz_xxxxxx[k] * ab_y + g_x_0_xxxxz_xxxxxxy[k];

                g_x_0_xxxxyz_xxxxxy[k] = -g_x_0_xxxxz_xxxxxy[k] * ab_y + g_x_0_xxxxz_xxxxxyy[k];

                g_x_0_xxxxyz_xxxxxz[k] = -g_x_0_xxxxz_xxxxxz[k] * ab_y + g_x_0_xxxxz_xxxxxyz[k];

                g_x_0_xxxxyz_xxxxyy[k] = -g_x_0_xxxxz_xxxxyy[k] * ab_y + g_x_0_xxxxz_xxxxyyy[k];

                g_x_0_xxxxyz_xxxxyz[k] = -g_x_0_xxxxz_xxxxyz[k] * ab_y + g_x_0_xxxxz_xxxxyyz[k];

                g_x_0_xxxxyz_xxxxzz[k] = -g_x_0_xxxxz_xxxxzz[k] * ab_y + g_x_0_xxxxz_xxxxyzz[k];

                g_x_0_xxxxyz_xxxyyy[k] = -g_x_0_xxxxz_xxxyyy[k] * ab_y + g_x_0_xxxxz_xxxyyyy[k];

                g_x_0_xxxxyz_xxxyyz[k] = -g_x_0_xxxxz_xxxyyz[k] * ab_y + g_x_0_xxxxz_xxxyyyz[k];

                g_x_0_xxxxyz_xxxyzz[k] = -g_x_0_xxxxz_xxxyzz[k] * ab_y + g_x_0_xxxxz_xxxyyzz[k];

                g_x_0_xxxxyz_xxxzzz[k] = -g_x_0_xxxxz_xxxzzz[k] * ab_y + g_x_0_xxxxz_xxxyzzz[k];

                g_x_0_xxxxyz_xxyyyy[k] = -g_x_0_xxxxz_xxyyyy[k] * ab_y + g_x_0_xxxxz_xxyyyyy[k];

                g_x_0_xxxxyz_xxyyyz[k] = -g_x_0_xxxxz_xxyyyz[k] * ab_y + g_x_0_xxxxz_xxyyyyz[k];

                g_x_0_xxxxyz_xxyyzz[k] = -g_x_0_xxxxz_xxyyzz[k] * ab_y + g_x_0_xxxxz_xxyyyzz[k];

                g_x_0_xxxxyz_xxyzzz[k] = -g_x_0_xxxxz_xxyzzz[k] * ab_y + g_x_0_xxxxz_xxyyzzz[k];

                g_x_0_xxxxyz_xxzzzz[k] = -g_x_0_xxxxz_xxzzzz[k] * ab_y + g_x_0_xxxxz_xxyzzzz[k];

                g_x_0_xxxxyz_xyyyyy[k] = -g_x_0_xxxxz_xyyyyy[k] * ab_y + g_x_0_xxxxz_xyyyyyy[k];

                g_x_0_xxxxyz_xyyyyz[k] = -g_x_0_xxxxz_xyyyyz[k] * ab_y + g_x_0_xxxxz_xyyyyyz[k];

                g_x_0_xxxxyz_xyyyzz[k] = -g_x_0_xxxxz_xyyyzz[k] * ab_y + g_x_0_xxxxz_xyyyyzz[k];

                g_x_0_xxxxyz_xyyzzz[k] = -g_x_0_xxxxz_xyyzzz[k] * ab_y + g_x_0_xxxxz_xyyyzzz[k];

                g_x_0_xxxxyz_xyzzzz[k] = -g_x_0_xxxxz_xyzzzz[k] * ab_y + g_x_0_xxxxz_xyyzzzz[k];

                g_x_0_xxxxyz_xzzzzz[k] = -g_x_0_xxxxz_xzzzzz[k] * ab_y + g_x_0_xxxxz_xyzzzzz[k];

                g_x_0_xxxxyz_yyyyyy[k] = -g_x_0_xxxxz_yyyyyy[k] * ab_y + g_x_0_xxxxz_yyyyyyy[k];

                g_x_0_xxxxyz_yyyyyz[k] = -g_x_0_xxxxz_yyyyyz[k] * ab_y + g_x_0_xxxxz_yyyyyyz[k];

                g_x_0_xxxxyz_yyyyzz[k] = -g_x_0_xxxxz_yyyyzz[k] * ab_y + g_x_0_xxxxz_yyyyyzz[k];

                g_x_0_xxxxyz_yyyzzz[k] = -g_x_0_xxxxz_yyyzzz[k] * ab_y + g_x_0_xxxxz_yyyyzzz[k];

                g_x_0_xxxxyz_yyzzzz[k] = -g_x_0_xxxxz_yyzzzz[k] * ab_y + g_x_0_xxxxz_yyyzzzz[k];

                g_x_0_xxxxyz_yzzzzz[k] = -g_x_0_xxxxz_yzzzzz[k] * ab_y + g_x_0_xxxxz_yyzzzzz[k];

                g_x_0_xxxxyz_zzzzzz[k] = -g_x_0_xxxxz_zzzzzz[k] * ab_y + g_x_0_xxxxz_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_xxxxxx = cbuffer.data(ii_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxxxy = cbuffer.data(ii_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxxxz = cbuffer.data(ii_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxxyy = cbuffer.data(ii_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxxyz = cbuffer.data(ii_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxxzz = cbuffer.data(ii_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxyyy = cbuffer.data(ii_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxyyz = cbuffer.data(ii_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxyzz = cbuffer.data(ii_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxzzz = cbuffer.data(ii_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxyyyy = cbuffer.data(ii_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxyyyz = cbuffer.data(ii_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxyyzz = cbuffer.data(ii_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxyzzz = cbuffer.data(ii_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxzzzz = cbuffer.data(ii_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyyyyy = cbuffer.data(ii_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyyyyz = cbuffer.data(ii_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyyyzz = cbuffer.data(ii_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyyzzz = cbuffer.data(ii_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyzzzz = cbuffer.data(ii_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xzzzzz = cbuffer.data(ii_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyyyyy = cbuffer.data(ii_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyyyyz = cbuffer.data(ii_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyyyzz = cbuffer.data(ii_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyyzzz = cbuffer.data(ii_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyzzzz = cbuffer.data(ii_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yzzzzz = cbuffer.data(ii_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxxxzz_zzzzzz = cbuffer.data(ii_geom_10_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxz_xxxxxx, g_x_0_xxxxz_xxxxxxz, g_x_0_xxxxz_xxxxxy, g_x_0_xxxxz_xxxxxyz, g_x_0_xxxxz_xxxxxz, g_x_0_xxxxz_xxxxxzz, g_x_0_xxxxz_xxxxyy, g_x_0_xxxxz_xxxxyyz, g_x_0_xxxxz_xxxxyz, g_x_0_xxxxz_xxxxyzz, g_x_0_xxxxz_xxxxzz, g_x_0_xxxxz_xxxxzzz, g_x_0_xxxxz_xxxyyy, g_x_0_xxxxz_xxxyyyz, g_x_0_xxxxz_xxxyyz, g_x_0_xxxxz_xxxyyzz, g_x_0_xxxxz_xxxyzz, g_x_0_xxxxz_xxxyzzz, g_x_0_xxxxz_xxxzzz, g_x_0_xxxxz_xxxzzzz, g_x_0_xxxxz_xxyyyy, g_x_0_xxxxz_xxyyyyz, g_x_0_xxxxz_xxyyyz, g_x_0_xxxxz_xxyyyzz, g_x_0_xxxxz_xxyyzz, g_x_0_xxxxz_xxyyzzz, g_x_0_xxxxz_xxyzzz, g_x_0_xxxxz_xxyzzzz, g_x_0_xxxxz_xxzzzz, g_x_0_xxxxz_xxzzzzz, g_x_0_xxxxz_xyyyyy, g_x_0_xxxxz_xyyyyyz, g_x_0_xxxxz_xyyyyz, g_x_0_xxxxz_xyyyyzz, g_x_0_xxxxz_xyyyzz, g_x_0_xxxxz_xyyyzzz, g_x_0_xxxxz_xyyzzz, g_x_0_xxxxz_xyyzzzz, g_x_0_xxxxz_xyzzzz, g_x_0_xxxxz_xyzzzzz, g_x_0_xxxxz_xzzzzz, g_x_0_xxxxz_xzzzzzz, g_x_0_xxxxz_yyyyyy, g_x_0_xxxxz_yyyyyyz, g_x_0_xxxxz_yyyyyz, g_x_0_xxxxz_yyyyyzz, g_x_0_xxxxz_yyyyzz, g_x_0_xxxxz_yyyyzzz, g_x_0_xxxxz_yyyzzz, g_x_0_xxxxz_yyyzzzz, g_x_0_xxxxz_yyzzzz, g_x_0_xxxxz_yyzzzzz, g_x_0_xxxxz_yzzzzz, g_x_0_xxxxz_yzzzzzz, g_x_0_xxxxz_zzzzzz, g_x_0_xxxxz_zzzzzzz, g_x_0_xxxxzz_xxxxxx, g_x_0_xxxxzz_xxxxxy, g_x_0_xxxxzz_xxxxxz, g_x_0_xxxxzz_xxxxyy, g_x_0_xxxxzz_xxxxyz, g_x_0_xxxxzz_xxxxzz, g_x_0_xxxxzz_xxxyyy, g_x_0_xxxxzz_xxxyyz, g_x_0_xxxxzz_xxxyzz, g_x_0_xxxxzz_xxxzzz, g_x_0_xxxxzz_xxyyyy, g_x_0_xxxxzz_xxyyyz, g_x_0_xxxxzz_xxyyzz, g_x_0_xxxxzz_xxyzzz, g_x_0_xxxxzz_xxzzzz, g_x_0_xxxxzz_xyyyyy, g_x_0_xxxxzz_xyyyyz, g_x_0_xxxxzz_xyyyzz, g_x_0_xxxxzz_xyyzzz, g_x_0_xxxxzz_xyzzzz, g_x_0_xxxxzz_xzzzzz, g_x_0_xxxxzz_yyyyyy, g_x_0_xxxxzz_yyyyyz, g_x_0_xxxxzz_yyyyzz, g_x_0_xxxxzz_yyyzzz, g_x_0_xxxxzz_yyzzzz, g_x_0_xxxxzz_yzzzzz, g_x_0_xxxxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_xxxxxx[k] = -g_x_0_xxxxz_xxxxxx[k] * ab_z + g_x_0_xxxxz_xxxxxxz[k];

                g_x_0_xxxxzz_xxxxxy[k] = -g_x_0_xxxxz_xxxxxy[k] * ab_z + g_x_0_xxxxz_xxxxxyz[k];

                g_x_0_xxxxzz_xxxxxz[k] = -g_x_0_xxxxz_xxxxxz[k] * ab_z + g_x_0_xxxxz_xxxxxzz[k];

                g_x_0_xxxxzz_xxxxyy[k] = -g_x_0_xxxxz_xxxxyy[k] * ab_z + g_x_0_xxxxz_xxxxyyz[k];

                g_x_0_xxxxzz_xxxxyz[k] = -g_x_0_xxxxz_xxxxyz[k] * ab_z + g_x_0_xxxxz_xxxxyzz[k];

                g_x_0_xxxxzz_xxxxzz[k] = -g_x_0_xxxxz_xxxxzz[k] * ab_z + g_x_0_xxxxz_xxxxzzz[k];

                g_x_0_xxxxzz_xxxyyy[k] = -g_x_0_xxxxz_xxxyyy[k] * ab_z + g_x_0_xxxxz_xxxyyyz[k];

                g_x_0_xxxxzz_xxxyyz[k] = -g_x_0_xxxxz_xxxyyz[k] * ab_z + g_x_0_xxxxz_xxxyyzz[k];

                g_x_0_xxxxzz_xxxyzz[k] = -g_x_0_xxxxz_xxxyzz[k] * ab_z + g_x_0_xxxxz_xxxyzzz[k];

                g_x_0_xxxxzz_xxxzzz[k] = -g_x_0_xxxxz_xxxzzz[k] * ab_z + g_x_0_xxxxz_xxxzzzz[k];

                g_x_0_xxxxzz_xxyyyy[k] = -g_x_0_xxxxz_xxyyyy[k] * ab_z + g_x_0_xxxxz_xxyyyyz[k];

                g_x_0_xxxxzz_xxyyyz[k] = -g_x_0_xxxxz_xxyyyz[k] * ab_z + g_x_0_xxxxz_xxyyyzz[k];

                g_x_0_xxxxzz_xxyyzz[k] = -g_x_0_xxxxz_xxyyzz[k] * ab_z + g_x_0_xxxxz_xxyyzzz[k];

                g_x_0_xxxxzz_xxyzzz[k] = -g_x_0_xxxxz_xxyzzz[k] * ab_z + g_x_0_xxxxz_xxyzzzz[k];

                g_x_0_xxxxzz_xxzzzz[k] = -g_x_0_xxxxz_xxzzzz[k] * ab_z + g_x_0_xxxxz_xxzzzzz[k];

                g_x_0_xxxxzz_xyyyyy[k] = -g_x_0_xxxxz_xyyyyy[k] * ab_z + g_x_0_xxxxz_xyyyyyz[k];

                g_x_0_xxxxzz_xyyyyz[k] = -g_x_0_xxxxz_xyyyyz[k] * ab_z + g_x_0_xxxxz_xyyyyzz[k];

                g_x_0_xxxxzz_xyyyzz[k] = -g_x_0_xxxxz_xyyyzz[k] * ab_z + g_x_0_xxxxz_xyyyzzz[k];

                g_x_0_xxxxzz_xyyzzz[k] = -g_x_0_xxxxz_xyyzzz[k] * ab_z + g_x_0_xxxxz_xyyzzzz[k];

                g_x_0_xxxxzz_xyzzzz[k] = -g_x_0_xxxxz_xyzzzz[k] * ab_z + g_x_0_xxxxz_xyzzzzz[k];

                g_x_0_xxxxzz_xzzzzz[k] = -g_x_0_xxxxz_xzzzzz[k] * ab_z + g_x_0_xxxxz_xzzzzzz[k];

                g_x_0_xxxxzz_yyyyyy[k] = -g_x_0_xxxxz_yyyyyy[k] * ab_z + g_x_0_xxxxz_yyyyyyz[k];

                g_x_0_xxxxzz_yyyyyz[k] = -g_x_0_xxxxz_yyyyyz[k] * ab_z + g_x_0_xxxxz_yyyyyzz[k];

                g_x_0_xxxxzz_yyyyzz[k] = -g_x_0_xxxxz_yyyyzz[k] * ab_z + g_x_0_xxxxz_yyyyzzz[k];

                g_x_0_xxxxzz_yyyzzz[k] = -g_x_0_xxxxz_yyyzzz[k] * ab_z + g_x_0_xxxxz_yyyzzzz[k];

                g_x_0_xxxxzz_yyzzzz[k] = -g_x_0_xxxxz_yyzzzz[k] * ab_z + g_x_0_xxxxz_yyzzzzz[k];

                g_x_0_xxxxzz_yzzzzz[k] = -g_x_0_xxxxz_yzzzzz[k] * ab_z + g_x_0_xxxxz_yzzzzzz[k];

                g_x_0_xxxxzz_zzzzzz[k] = -g_x_0_xxxxz_zzzzzz[k] * ab_z + g_x_0_xxxxz_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xxxyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyy_xxxxxx, g_x_0_xxxyy_xxxxxxy, g_x_0_xxxyy_xxxxxy, g_x_0_xxxyy_xxxxxyy, g_x_0_xxxyy_xxxxxyz, g_x_0_xxxyy_xxxxxz, g_x_0_xxxyy_xxxxyy, g_x_0_xxxyy_xxxxyyy, g_x_0_xxxyy_xxxxyyz, g_x_0_xxxyy_xxxxyz, g_x_0_xxxyy_xxxxyzz, g_x_0_xxxyy_xxxxzz, g_x_0_xxxyy_xxxyyy, g_x_0_xxxyy_xxxyyyy, g_x_0_xxxyy_xxxyyyz, g_x_0_xxxyy_xxxyyz, g_x_0_xxxyy_xxxyyzz, g_x_0_xxxyy_xxxyzz, g_x_0_xxxyy_xxxyzzz, g_x_0_xxxyy_xxxzzz, g_x_0_xxxyy_xxyyyy, g_x_0_xxxyy_xxyyyyy, g_x_0_xxxyy_xxyyyyz, g_x_0_xxxyy_xxyyyz, g_x_0_xxxyy_xxyyyzz, g_x_0_xxxyy_xxyyzz, g_x_0_xxxyy_xxyyzzz, g_x_0_xxxyy_xxyzzz, g_x_0_xxxyy_xxyzzzz, g_x_0_xxxyy_xxzzzz, g_x_0_xxxyy_xyyyyy, g_x_0_xxxyy_xyyyyyy, g_x_0_xxxyy_xyyyyyz, g_x_0_xxxyy_xyyyyz, g_x_0_xxxyy_xyyyyzz, g_x_0_xxxyy_xyyyzz, g_x_0_xxxyy_xyyyzzz, g_x_0_xxxyy_xyyzzz, g_x_0_xxxyy_xyyzzzz, g_x_0_xxxyy_xyzzzz, g_x_0_xxxyy_xyzzzzz, g_x_0_xxxyy_xzzzzz, g_x_0_xxxyy_yyyyyy, g_x_0_xxxyy_yyyyyyy, g_x_0_xxxyy_yyyyyyz, g_x_0_xxxyy_yyyyyz, g_x_0_xxxyy_yyyyyzz, g_x_0_xxxyy_yyyyzz, g_x_0_xxxyy_yyyyzzz, g_x_0_xxxyy_yyyzzz, g_x_0_xxxyy_yyyzzzz, g_x_0_xxxyy_yyzzzz, g_x_0_xxxyy_yyzzzzz, g_x_0_xxxyy_yzzzzz, g_x_0_xxxyy_yzzzzzz, g_x_0_xxxyy_zzzzzz, g_x_0_xxxyyy_xxxxxx, g_x_0_xxxyyy_xxxxxy, g_x_0_xxxyyy_xxxxxz, g_x_0_xxxyyy_xxxxyy, g_x_0_xxxyyy_xxxxyz, g_x_0_xxxyyy_xxxxzz, g_x_0_xxxyyy_xxxyyy, g_x_0_xxxyyy_xxxyyz, g_x_0_xxxyyy_xxxyzz, g_x_0_xxxyyy_xxxzzz, g_x_0_xxxyyy_xxyyyy, g_x_0_xxxyyy_xxyyyz, g_x_0_xxxyyy_xxyyzz, g_x_0_xxxyyy_xxyzzz, g_x_0_xxxyyy_xxzzzz, g_x_0_xxxyyy_xyyyyy, g_x_0_xxxyyy_xyyyyz, g_x_0_xxxyyy_xyyyzz, g_x_0_xxxyyy_xyyzzz, g_x_0_xxxyyy_xyzzzz, g_x_0_xxxyyy_xzzzzz, g_x_0_xxxyyy_yyyyyy, g_x_0_xxxyyy_yyyyyz, g_x_0_xxxyyy_yyyyzz, g_x_0_xxxyyy_yyyzzz, g_x_0_xxxyyy_yyzzzz, g_x_0_xxxyyy_yzzzzz, g_x_0_xxxyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_xxxxxx[k] = -g_x_0_xxxyy_xxxxxx[k] * ab_y + g_x_0_xxxyy_xxxxxxy[k];

                g_x_0_xxxyyy_xxxxxy[k] = -g_x_0_xxxyy_xxxxxy[k] * ab_y + g_x_0_xxxyy_xxxxxyy[k];

                g_x_0_xxxyyy_xxxxxz[k] = -g_x_0_xxxyy_xxxxxz[k] * ab_y + g_x_0_xxxyy_xxxxxyz[k];

                g_x_0_xxxyyy_xxxxyy[k] = -g_x_0_xxxyy_xxxxyy[k] * ab_y + g_x_0_xxxyy_xxxxyyy[k];

                g_x_0_xxxyyy_xxxxyz[k] = -g_x_0_xxxyy_xxxxyz[k] * ab_y + g_x_0_xxxyy_xxxxyyz[k];

                g_x_0_xxxyyy_xxxxzz[k] = -g_x_0_xxxyy_xxxxzz[k] * ab_y + g_x_0_xxxyy_xxxxyzz[k];

                g_x_0_xxxyyy_xxxyyy[k] = -g_x_0_xxxyy_xxxyyy[k] * ab_y + g_x_0_xxxyy_xxxyyyy[k];

                g_x_0_xxxyyy_xxxyyz[k] = -g_x_0_xxxyy_xxxyyz[k] * ab_y + g_x_0_xxxyy_xxxyyyz[k];

                g_x_0_xxxyyy_xxxyzz[k] = -g_x_0_xxxyy_xxxyzz[k] * ab_y + g_x_0_xxxyy_xxxyyzz[k];

                g_x_0_xxxyyy_xxxzzz[k] = -g_x_0_xxxyy_xxxzzz[k] * ab_y + g_x_0_xxxyy_xxxyzzz[k];

                g_x_0_xxxyyy_xxyyyy[k] = -g_x_0_xxxyy_xxyyyy[k] * ab_y + g_x_0_xxxyy_xxyyyyy[k];

                g_x_0_xxxyyy_xxyyyz[k] = -g_x_0_xxxyy_xxyyyz[k] * ab_y + g_x_0_xxxyy_xxyyyyz[k];

                g_x_0_xxxyyy_xxyyzz[k] = -g_x_0_xxxyy_xxyyzz[k] * ab_y + g_x_0_xxxyy_xxyyyzz[k];

                g_x_0_xxxyyy_xxyzzz[k] = -g_x_0_xxxyy_xxyzzz[k] * ab_y + g_x_0_xxxyy_xxyyzzz[k];

                g_x_0_xxxyyy_xxzzzz[k] = -g_x_0_xxxyy_xxzzzz[k] * ab_y + g_x_0_xxxyy_xxyzzzz[k];

                g_x_0_xxxyyy_xyyyyy[k] = -g_x_0_xxxyy_xyyyyy[k] * ab_y + g_x_0_xxxyy_xyyyyyy[k];

                g_x_0_xxxyyy_xyyyyz[k] = -g_x_0_xxxyy_xyyyyz[k] * ab_y + g_x_0_xxxyy_xyyyyyz[k];

                g_x_0_xxxyyy_xyyyzz[k] = -g_x_0_xxxyy_xyyyzz[k] * ab_y + g_x_0_xxxyy_xyyyyzz[k];

                g_x_0_xxxyyy_xyyzzz[k] = -g_x_0_xxxyy_xyyzzz[k] * ab_y + g_x_0_xxxyy_xyyyzzz[k];

                g_x_0_xxxyyy_xyzzzz[k] = -g_x_0_xxxyy_xyzzzz[k] * ab_y + g_x_0_xxxyy_xyyzzzz[k];

                g_x_0_xxxyyy_xzzzzz[k] = -g_x_0_xxxyy_xzzzzz[k] * ab_y + g_x_0_xxxyy_xyzzzzz[k];

                g_x_0_xxxyyy_yyyyyy[k] = -g_x_0_xxxyy_yyyyyy[k] * ab_y + g_x_0_xxxyy_yyyyyyy[k];

                g_x_0_xxxyyy_yyyyyz[k] = -g_x_0_xxxyy_yyyyyz[k] * ab_y + g_x_0_xxxyy_yyyyyyz[k];

                g_x_0_xxxyyy_yyyyzz[k] = -g_x_0_xxxyy_yyyyzz[k] * ab_y + g_x_0_xxxyy_yyyyyzz[k];

                g_x_0_xxxyyy_yyyzzz[k] = -g_x_0_xxxyy_yyyzzz[k] * ab_y + g_x_0_xxxyy_yyyyzzz[k];

                g_x_0_xxxyyy_yyzzzz[k] = -g_x_0_xxxyy_yyzzzz[k] * ab_y + g_x_0_xxxyy_yyyzzzz[k];

                g_x_0_xxxyyy_yzzzzz[k] = -g_x_0_xxxyy_yzzzzz[k] * ab_y + g_x_0_xxxyy_yyzzzzz[k];

                g_x_0_xxxyyy_zzzzzz[k] = -g_x_0_xxxyy_zzzzzz[k] * ab_y + g_x_0_xxxyy_yzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xxxyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyyz_xxxxxx, g_x_0_xxxyyz_xxxxxy, g_x_0_xxxyyz_xxxxxz, g_x_0_xxxyyz_xxxxyy, g_x_0_xxxyyz_xxxxyz, g_x_0_xxxyyz_xxxxzz, g_x_0_xxxyyz_xxxyyy, g_x_0_xxxyyz_xxxyyz, g_x_0_xxxyyz_xxxyzz, g_x_0_xxxyyz_xxxzzz, g_x_0_xxxyyz_xxyyyy, g_x_0_xxxyyz_xxyyyz, g_x_0_xxxyyz_xxyyzz, g_x_0_xxxyyz_xxyzzz, g_x_0_xxxyyz_xxzzzz, g_x_0_xxxyyz_xyyyyy, g_x_0_xxxyyz_xyyyyz, g_x_0_xxxyyz_xyyyzz, g_x_0_xxxyyz_xyyzzz, g_x_0_xxxyyz_xyzzzz, g_x_0_xxxyyz_xzzzzz, g_x_0_xxxyyz_yyyyyy, g_x_0_xxxyyz_yyyyyz, g_x_0_xxxyyz_yyyyzz, g_x_0_xxxyyz_yyyzzz, g_x_0_xxxyyz_yyzzzz, g_x_0_xxxyyz_yzzzzz, g_x_0_xxxyyz_zzzzzz, g_x_0_xxxyz_xxxxxx, g_x_0_xxxyz_xxxxxxy, g_x_0_xxxyz_xxxxxy, g_x_0_xxxyz_xxxxxyy, g_x_0_xxxyz_xxxxxyz, g_x_0_xxxyz_xxxxxz, g_x_0_xxxyz_xxxxyy, g_x_0_xxxyz_xxxxyyy, g_x_0_xxxyz_xxxxyyz, g_x_0_xxxyz_xxxxyz, g_x_0_xxxyz_xxxxyzz, g_x_0_xxxyz_xxxxzz, g_x_0_xxxyz_xxxyyy, g_x_0_xxxyz_xxxyyyy, g_x_0_xxxyz_xxxyyyz, g_x_0_xxxyz_xxxyyz, g_x_0_xxxyz_xxxyyzz, g_x_0_xxxyz_xxxyzz, g_x_0_xxxyz_xxxyzzz, g_x_0_xxxyz_xxxzzz, g_x_0_xxxyz_xxyyyy, g_x_0_xxxyz_xxyyyyy, g_x_0_xxxyz_xxyyyyz, g_x_0_xxxyz_xxyyyz, g_x_0_xxxyz_xxyyyzz, g_x_0_xxxyz_xxyyzz, g_x_0_xxxyz_xxyyzzz, g_x_0_xxxyz_xxyzzz, g_x_0_xxxyz_xxyzzzz, g_x_0_xxxyz_xxzzzz, g_x_0_xxxyz_xyyyyy, g_x_0_xxxyz_xyyyyyy, g_x_0_xxxyz_xyyyyyz, g_x_0_xxxyz_xyyyyz, g_x_0_xxxyz_xyyyyzz, g_x_0_xxxyz_xyyyzz, g_x_0_xxxyz_xyyyzzz, g_x_0_xxxyz_xyyzzz, g_x_0_xxxyz_xyyzzzz, g_x_0_xxxyz_xyzzzz, g_x_0_xxxyz_xyzzzzz, g_x_0_xxxyz_xzzzzz, g_x_0_xxxyz_yyyyyy, g_x_0_xxxyz_yyyyyyy, g_x_0_xxxyz_yyyyyyz, g_x_0_xxxyz_yyyyyz, g_x_0_xxxyz_yyyyyzz, g_x_0_xxxyz_yyyyzz, g_x_0_xxxyz_yyyyzzz, g_x_0_xxxyz_yyyzzz, g_x_0_xxxyz_yyyzzzz, g_x_0_xxxyz_yyzzzz, g_x_0_xxxyz_yyzzzzz, g_x_0_xxxyz_yzzzzz, g_x_0_xxxyz_yzzzzzz, g_x_0_xxxyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_xxxxxx[k] = -g_x_0_xxxyz_xxxxxx[k] * ab_y + g_x_0_xxxyz_xxxxxxy[k];

                g_x_0_xxxyyz_xxxxxy[k] = -g_x_0_xxxyz_xxxxxy[k] * ab_y + g_x_0_xxxyz_xxxxxyy[k];

                g_x_0_xxxyyz_xxxxxz[k] = -g_x_0_xxxyz_xxxxxz[k] * ab_y + g_x_0_xxxyz_xxxxxyz[k];

                g_x_0_xxxyyz_xxxxyy[k] = -g_x_0_xxxyz_xxxxyy[k] * ab_y + g_x_0_xxxyz_xxxxyyy[k];

                g_x_0_xxxyyz_xxxxyz[k] = -g_x_0_xxxyz_xxxxyz[k] * ab_y + g_x_0_xxxyz_xxxxyyz[k];

                g_x_0_xxxyyz_xxxxzz[k] = -g_x_0_xxxyz_xxxxzz[k] * ab_y + g_x_0_xxxyz_xxxxyzz[k];

                g_x_0_xxxyyz_xxxyyy[k] = -g_x_0_xxxyz_xxxyyy[k] * ab_y + g_x_0_xxxyz_xxxyyyy[k];

                g_x_0_xxxyyz_xxxyyz[k] = -g_x_0_xxxyz_xxxyyz[k] * ab_y + g_x_0_xxxyz_xxxyyyz[k];

                g_x_0_xxxyyz_xxxyzz[k] = -g_x_0_xxxyz_xxxyzz[k] * ab_y + g_x_0_xxxyz_xxxyyzz[k];

                g_x_0_xxxyyz_xxxzzz[k] = -g_x_0_xxxyz_xxxzzz[k] * ab_y + g_x_0_xxxyz_xxxyzzz[k];

                g_x_0_xxxyyz_xxyyyy[k] = -g_x_0_xxxyz_xxyyyy[k] * ab_y + g_x_0_xxxyz_xxyyyyy[k];

                g_x_0_xxxyyz_xxyyyz[k] = -g_x_0_xxxyz_xxyyyz[k] * ab_y + g_x_0_xxxyz_xxyyyyz[k];

                g_x_0_xxxyyz_xxyyzz[k] = -g_x_0_xxxyz_xxyyzz[k] * ab_y + g_x_0_xxxyz_xxyyyzz[k];

                g_x_0_xxxyyz_xxyzzz[k] = -g_x_0_xxxyz_xxyzzz[k] * ab_y + g_x_0_xxxyz_xxyyzzz[k];

                g_x_0_xxxyyz_xxzzzz[k] = -g_x_0_xxxyz_xxzzzz[k] * ab_y + g_x_0_xxxyz_xxyzzzz[k];

                g_x_0_xxxyyz_xyyyyy[k] = -g_x_0_xxxyz_xyyyyy[k] * ab_y + g_x_0_xxxyz_xyyyyyy[k];

                g_x_0_xxxyyz_xyyyyz[k] = -g_x_0_xxxyz_xyyyyz[k] * ab_y + g_x_0_xxxyz_xyyyyyz[k];

                g_x_0_xxxyyz_xyyyzz[k] = -g_x_0_xxxyz_xyyyzz[k] * ab_y + g_x_0_xxxyz_xyyyyzz[k];

                g_x_0_xxxyyz_xyyzzz[k] = -g_x_0_xxxyz_xyyzzz[k] * ab_y + g_x_0_xxxyz_xyyyzzz[k];

                g_x_0_xxxyyz_xyzzzz[k] = -g_x_0_xxxyz_xyzzzz[k] * ab_y + g_x_0_xxxyz_xyyzzzz[k];

                g_x_0_xxxyyz_xzzzzz[k] = -g_x_0_xxxyz_xzzzzz[k] * ab_y + g_x_0_xxxyz_xyzzzzz[k];

                g_x_0_xxxyyz_yyyyyy[k] = -g_x_0_xxxyz_yyyyyy[k] * ab_y + g_x_0_xxxyz_yyyyyyy[k];

                g_x_0_xxxyyz_yyyyyz[k] = -g_x_0_xxxyz_yyyyyz[k] * ab_y + g_x_0_xxxyz_yyyyyyz[k];

                g_x_0_xxxyyz_yyyyzz[k] = -g_x_0_xxxyz_yyyyzz[k] * ab_y + g_x_0_xxxyz_yyyyyzz[k];

                g_x_0_xxxyyz_yyyzzz[k] = -g_x_0_xxxyz_yyyzzz[k] * ab_y + g_x_0_xxxyz_yyyyzzz[k];

                g_x_0_xxxyyz_yyzzzz[k] = -g_x_0_xxxyz_yyzzzz[k] * ab_y + g_x_0_xxxyz_yyyzzzz[k];

                g_x_0_xxxyyz_yzzzzz[k] = -g_x_0_xxxyz_yzzzzz[k] * ab_y + g_x_0_xxxyz_yyzzzzz[k];

                g_x_0_xxxyyz_zzzzzz[k] = -g_x_0_xxxyz_zzzzzz[k] * ab_y + g_x_0_xxxyz_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xxxyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyzz_xxxxxx, g_x_0_xxxyzz_xxxxxy, g_x_0_xxxyzz_xxxxxz, g_x_0_xxxyzz_xxxxyy, g_x_0_xxxyzz_xxxxyz, g_x_0_xxxyzz_xxxxzz, g_x_0_xxxyzz_xxxyyy, g_x_0_xxxyzz_xxxyyz, g_x_0_xxxyzz_xxxyzz, g_x_0_xxxyzz_xxxzzz, g_x_0_xxxyzz_xxyyyy, g_x_0_xxxyzz_xxyyyz, g_x_0_xxxyzz_xxyyzz, g_x_0_xxxyzz_xxyzzz, g_x_0_xxxyzz_xxzzzz, g_x_0_xxxyzz_xyyyyy, g_x_0_xxxyzz_xyyyyz, g_x_0_xxxyzz_xyyyzz, g_x_0_xxxyzz_xyyzzz, g_x_0_xxxyzz_xyzzzz, g_x_0_xxxyzz_xzzzzz, g_x_0_xxxyzz_yyyyyy, g_x_0_xxxyzz_yyyyyz, g_x_0_xxxyzz_yyyyzz, g_x_0_xxxyzz_yyyzzz, g_x_0_xxxyzz_yyzzzz, g_x_0_xxxyzz_yzzzzz, g_x_0_xxxyzz_zzzzzz, g_x_0_xxxzz_xxxxxx, g_x_0_xxxzz_xxxxxxy, g_x_0_xxxzz_xxxxxy, g_x_0_xxxzz_xxxxxyy, g_x_0_xxxzz_xxxxxyz, g_x_0_xxxzz_xxxxxz, g_x_0_xxxzz_xxxxyy, g_x_0_xxxzz_xxxxyyy, g_x_0_xxxzz_xxxxyyz, g_x_0_xxxzz_xxxxyz, g_x_0_xxxzz_xxxxyzz, g_x_0_xxxzz_xxxxzz, g_x_0_xxxzz_xxxyyy, g_x_0_xxxzz_xxxyyyy, g_x_0_xxxzz_xxxyyyz, g_x_0_xxxzz_xxxyyz, g_x_0_xxxzz_xxxyyzz, g_x_0_xxxzz_xxxyzz, g_x_0_xxxzz_xxxyzzz, g_x_0_xxxzz_xxxzzz, g_x_0_xxxzz_xxyyyy, g_x_0_xxxzz_xxyyyyy, g_x_0_xxxzz_xxyyyyz, g_x_0_xxxzz_xxyyyz, g_x_0_xxxzz_xxyyyzz, g_x_0_xxxzz_xxyyzz, g_x_0_xxxzz_xxyyzzz, g_x_0_xxxzz_xxyzzz, g_x_0_xxxzz_xxyzzzz, g_x_0_xxxzz_xxzzzz, g_x_0_xxxzz_xyyyyy, g_x_0_xxxzz_xyyyyyy, g_x_0_xxxzz_xyyyyyz, g_x_0_xxxzz_xyyyyz, g_x_0_xxxzz_xyyyyzz, g_x_0_xxxzz_xyyyzz, g_x_0_xxxzz_xyyyzzz, g_x_0_xxxzz_xyyzzz, g_x_0_xxxzz_xyyzzzz, g_x_0_xxxzz_xyzzzz, g_x_0_xxxzz_xyzzzzz, g_x_0_xxxzz_xzzzzz, g_x_0_xxxzz_yyyyyy, g_x_0_xxxzz_yyyyyyy, g_x_0_xxxzz_yyyyyyz, g_x_0_xxxzz_yyyyyz, g_x_0_xxxzz_yyyyyzz, g_x_0_xxxzz_yyyyzz, g_x_0_xxxzz_yyyyzzz, g_x_0_xxxzz_yyyzzz, g_x_0_xxxzz_yyyzzzz, g_x_0_xxxzz_yyzzzz, g_x_0_xxxzz_yyzzzzz, g_x_0_xxxzz_yzzzzz, g_x_0_xxxzz_yzzzzzz, g_x_0_xxxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_xxxxxx[k] = -g_x_0_xxxzz_xxxxxx[k] * ab_y + g_x_0_xxxzz_xxxxxxy[k];

                g_x_0_xxxyzz_xxxxxy[k] = -g_x_0_xxxzz_xxxxxy[k] * ab_y + g_x_0_xxxzz_xxxxxyy[k];

                g_x_0_xxxyzz_xxxxxz[k] = -g_x_0_xxxzz_xxxxxz[k] * ab_y + g_x_0_xxxzz_xxxxxyz[k];

                g_x_0_xxxyzz_xxxxyy[k] = -g_x_0_xxxzz_xxxxyy[k] * ab_y + g_x_0_xxxzz_xxxxyyy[k];

                g_x_0_xxxyzz_xxxxyz[k] = -g_x_0_xxxzz_xxxxyz[k] * ab_y + g_x_0_xxxzz_xxxxyyz[k];

                g_x_0_xxxyzz_xxxxzz[k] = -g_x_0_xxxzz_xxxxzz[k] * ab_y + g_x_0_xxxzz_xxxxyzz[k];

                g_x_0_xxxyzz_xxxyyy[k] = -g_x_0_xxxzz_xxxyyy[k] * ab_y + g_x_0_xxxzz_xxxyyyy[k];

                g_x_0_xxxyzz_xxxyyz[k] = -g_x_0_xxxzz_xxxyyz[k] * ab_y + g_x_0_xxxzz_xxxyyyz[k];

                g_x_0_xxxyzz_xxxyzz[k] = -g_x_0_xxxzz_xxxyzz[k] * ab_y + g_x_0_xxxzz_xxxyyzz[k];

                g_x_0_xxxyzz_xxxzzz[k] = -g_x_0_xxxzz_xxxzzz[k] * ab_y + g_x_0_xxxzz_xxxyzzz[k];

                g_x_0_xxxyzz_xxyyyy[k] = -g_x_0_xxxzz_xxyyyy[k] * ab_y + g_x_0_xxxzz_xxyyyyy[k];

                g_x_0_xxxyzz_xxyyyz[k] = -g_x_0_xxxzz_xxyyyz[k] * ab_y + g_x_0_xxxzz_xxyyyyz[k];

                g_x_0_xxxyzz_xxyyzz[k] = -g_x_0_xxxzz_xxyyzz[k] * ab_y + g_x_0_xxxzz_xxyyyzz[k];

                g_x_0_xxxyzz_xxyzzz[k] = -g_x_0_xxxzz_xxyzzz[k] * ab_y + g_x_0_xxxzz_xxyyzzz[k];

                g_x_0_xxxyzz_xxzzzz[k] = -g_x_0_xxxzz_xxzzzz[k] * ab_y + g_x_0_xxxzz_xxyzzzz[k];

                g_x_0_xxxyzz_xyyyyy[k] = -g_x_0_xxxzz_xyyyyy[k] * ab_y + g_x_0_xxxzz_xyyyyyy[k];

                g_x_0_xxxyzz_xyyyyz[k] = -g_x_0_xxxzz_xyyyyz[k] * ab_y + g_x_0_xxxzz_xyyyyyz[k];

                g_x_0_xxxyzz_xyyyzz[k] = -g_x_0_xxxzz_xyyyzz[k] * ab_y + g_x_0_xxxzz_xyyyyzz[k];

                g_x_0_xxxyzz_xyyzzz[k] = -g_x_0_xxxzz_xyyzzz[k] * ab_y + g_x_0_xxxzz_xyyyzzz[k];

                g_x_0_xxxyzz_xyzzzz[k] = -g_x_0_xxxzz_xyzzzz[k] * ab_y + g_x_0_xxxzz_xyyzzzz[k];

                g_x_0_xxxyzz_xzzzzz[k] = -g_x_0_xxxzz_xzzzzz[k] * ab_y + g_x_0_xxxzz_xyzzzzz[k];

                g_x_0_xxxyzz_yyyyyy[k] = -g_x_0_xxxzz_yyyyyy[k] * ab_y + g_x_0_xxxzz_yyyyyyy[k];

                g_x_0_xxxyzz_yyyyyz[k] = -g_x_0_xxxzz_yyyyyz[k] * ab_y + g_x_0_xxxzz_yyyyyyz[k];

                g_x_0_xxxyzz_yyyyzz[k] = -g_x_0_xxxzz_yyyyzz[k] * ab_y + g_x_0_xxxzz_yyyyyzz[k];

                g_x_0_xxxyzz_yyyzzz[k] = -g_x_0_xxxzz_yyyzzz[k] * ab_y + g_x_0_xxxzz_yyyyzzz[k];

                g_x_0_xxxyzz_yyzzzz[k] = -g_x_0_xxxzz_yyzzzz[k] * ab_y + g_x_0_xxxzz_yyyzzzz[k];

                g_x_0_xxxyzz_yzzzzz[k] = -g_x_0_xxxzz_yzzzzz[k] * ab_y + g_x_0_xxxzz_yyzzzzz[k];

                g_x_0_xxxyzz_zzzzzz[k] = -g_x_0_xxxzz_zzzzzz[k] * ab_y + g_x_0_xxxzz_yzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xxxzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxzz_xxxxxx, g_x_0_xxxzz_xxxxxxz, g_x_0_xxxzz_xxxxxy, g_x_0_xxxzz_xxxxxyz, g_x_0_xxxzz_xxxxxz, g_x_0_xxxzz_xxxxxzz, g_x_0_xxxzz_xxxxyy, g_x_0_xxxzz_xxxxyyz, g_x_0_xxxzz_xxxxyz, g_x_0_xxxzz_xxxxyzz, g_x_0_xxxzz_xxxxzz, g_x_0_xxxzz_xxxxzzz, g_x_0_xxxzz_xxxyyy, g_x_0_xxxzz_xxxyyyz, g_x_0_xxxzz_xxxyyz, g_x_0_xxxzz_xxxyyzz, g_x_0_xxxzz_xxxyzz, g_x_0_xxxzz_xxxyzzz, g_x_0_xxxzz_xxxzzz, g_x_0_xxxzz_xxxzzzz, g_x_0_xxxzz_xxyyyy, g_x_0_xxxzz_xxyyyyz, g_x_0_xxxzz_xxyyyz, g_x_0_xxxzz_xxyyyzz, g_x_0_xxxzz_xxyyzz, g_x_0_xxxzz_xxyyzzz, g_x_0_xxxzz_xxyzzz, g_x_0_xxxzz_xxyzzzz, g_x_0_xxxzz_xxzzzz, g_x_0_xxxzz_xxzzzzz, g_x_0_xxxzz_xyyyyy, g_x_0_xxxzz_xyyyyyz, g_x_0_xxxzz_xyyyyz, g_x_0_xxxzz_xyyyyzz, g_x_0_xxxzz_xyyyzz, g_x_0_xxxzz_xyyyzzz, g_x_0_xxxzz_xyyzzz, g_x_0_xxxzz_xyyzzzz, g_x_0_xxxzz_xyzzzz, g_x_0_xxxzz_xyzzzzz, g_x_0_xxxzz_xzzzzz, g_x_0_xxxzz_xzzzzzz, g_x_0_xxxzz_yyyyyy, g_x_0_xxxzz_yyyyyyz, g_x_0_xxxzz_yyyyyz, g_x_0_xxxzz_yyyyyzz, g_x_0_xxxzz_yyyyzz, g_x_0_xxxzz_yyyyzzz, g_x_0_xxxzz_yyyzzz, g_x_0_xxxzz_yyyzzzz, g_x_0_xxxzz_yyzzzz, g_x_0_xxxzz_yyzzzzz, g_x_0_xxxzz_yzzzzz, g_x_0_xxxzz_yzzzzzz, g_x_0_xxxzz_zzzzzz, g_x_0_xxxzz_zzzzzzz, g_x_0_xxxzzz_xxxxxx, g_x_0_xxxzzz_xxxxxy, g_x_0_xxxzzz_xxxxxz, g_x_0_xxxzzz_xxxxyy, g_x_0_xxxzzz_xxxxyz, g_x_0_xxxzzz_xxxxzz, g_x_0_xxxzzz_xxxyyy, g_x_0_xxxzzz_xxxyyz, g_x_0_xxxzzz_xxxyzz, g_x_0_xxxzzz_xxxzzz, g_x_0_xxxzzz_xxyyyy, g_x_0_xxxzzz_xxyyyz, g_x_0_xxxzzz_xxyyzz, g_x_0_xxxzzz_xxyzzz, g_x_0_xxxzzz_xxzzzz, g_x_0_xxxzzz_xyyyyy, g_x_0_xxxzzz_xyyyyz, g_x_0_xxxzzz_xyyyzz, g_x_0_xxxzzz_xyyzzz, g_x_0_xxxzzz_xyzzzz, g_x_0_xxxzzz_xzzzzz, g_x_0_xxxzzz_yyyyyy, g_x_0_xxxzzz_yyyyyz, g_x_0_xxxzzz_yyyyzz, g_x_0_xxxzzz_yyyzzz, g_x_0_xxxzzz_yyzzzz, g_x_0_xxxzzz_yzzzzz, g_x_0_xxxzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_xxxxxx[k] = -g_x_0_xxxzz_xxxxxx[k] * ab_z + g_x_0_xxxzz_xxxxxxz[k];

                g_x_0_xxxzzz_xxxxxy[k] = -g_x_0_xxxzz_xxxxxy[k] * ab_z + g_x_0_xxxzz_xxxxxyz[k];

                g_x_0_xxxzzz_xxxxxz[k] = -g_x_0_xxxzz_xxxxxz[k] * ab_z + g_x_0_xxxzz_xxxxxzz[k];

                g_x_0_xxxzzz_xxxxyy[k] = -g_x_0_xxxzz_xxxxyy[k] * ab_z + g_x_0_xxxzz_xxxxyyz[k];

                g_x_0_xxxzzz_xxxxyz[k] = -g_x_0_xxxzz_xxxxyz[k] * ab_z + g_x_0_xxxzz_xxxxyzz[k];

                g_x_0_xxxzzz_xxxxzz[k] = -g_x_0_xxxzz_xxxxzz[k] * ab_z + g_x_0_xxxzz_xxxxzzz[k];

                g_x_0_xxxzzz_xxxyyy[k] = -g_x_0_xxxzz_xxxyyy[k] * ab_z + g_x_0_xxxzz_xxxyyyz[k];

                g_x_0_xxxzzz_xxxyyz[k] = -g_x_0_xxxzz_xxxyyz[k] * ab_z + g_x_0_xxxzz_xxxyyzz[k];

                g_x_0_xxxzzz_xxxyzz[k] = -g_x_0_xxxzz_xxxyzz[k] * ab_z + g_x_0_xxxzz_xxxyzzz[k];

                g_x_0_xxxzzz_xxxzzz[k] = -g_x_0_xxxzz_xxxzzz[k] * ab_z + g_x_0_xxxzz_xxxzzzz[k];

                g_x_0_xxxzzz_xxyyyy[k] = -g_x_0_xxxzz_xxyyyy[k] * ab_z + g_x_0_xxxzz_xxyyyyz[k];

                g_x_0_xxxzzz_xxyyyz[k] = -g_x_0_xxxzz_xxyyyz[k] * ab_z + g_x_0_xxxzz_xxyyyzz[k];

                g_x_0_xxxzzz_xxyyzz[k] = -g_x_0_xxxzz_xxyyzz[k] * ab_z + g_x_0_xxxzz_xxyyzzz[k];

                g_x_0_xxxzzz_xxyzzz[k] = -g_x_0_xxxzz_xxyzzz[k] * ab_z + g_x_0_xxxzz_xxyzzzz[k];

                g_x_0_xxxzzz_xxzzzz[k] = -g_x_0_xxxzz_xxzzzz[k] * ab_z + g_x_0_xxxzz_xxzzzzz[k];

                g_x_0_xxxzzz_xyyyyy[k] = -g_x_0_xxxzz_xyyyyy[k] * ab_z + g_x_0_xxxzz_xyyyyyz[k];

                g_x_0_xxxzzz_xyyyyz[k] = -g_x_0_xxxzz_xyyyyz[k] * ab_z + g_x_0_xxxzz_xyyyyzz[k];

                g_x_0_xxxzzz_xyyyzz[k] = -g_x_0_xxxzz_xyyyzz[k] * ab_z + g_x_0_xxxzz_xyyyzzz[k];

                g_x_0_xxxzzz_xyyzzz[k] = -g_x_0_xxxzz_xyyzzz[k] * ab_z + g_x_0_xxxzz_xyyzzzz[k];

                g_x_0_xxxzzz_xyzzzz[k] = -g_x_0_xxxzz_xyzzzz[k] * ab_z + g_x_0_xxxzz_xyzzzzz[k];

                g_x_0_xxxzzz_xzzzzz[k] = -g_x_0_xxxzz_xzzzzz[k] * ab_z + g_x_0_xxxzz_xzzzzzz[k];

                g_x_0_xxxzzz_yyyyyy[k] = -g_x_0_xxxzz_yyyyyy[k] * ab_z + g_x_0_xxxzz_yyyyyyz[k];

                g_x_0_xxxzzz_yyyyyz[k] = -g_x_0_xxxzz_yyyyyz[k] * ab_z + g_x_0_xxxzz_yyyyyzz[k];

                g_x_0_xxxzzz_yyyyzz[k] = -g_x_0_xxxzz_yyyyzz[k] * ab_z + g_x_0_xxxzz_yyyyzzz[k];

                g_x_0_xxxzzz_yyyzzz[k] = -g_x_0_xxxzz_yyyzzz[k] * ab_z + g_x_0_xxxzz_yyyzzzz[k];

                g_x_0_xxxzzz_yyzzzz[k] = -g_x_0_xxxzz_yyzzzz[k] * ab_z + g_x_0_xxxzz_yyzzzzz[k];

                g_x_0_xxxzzz_yzzzzz[k] = -g_x_0_xxxzz_yzzzzz[k] * ab_z + g_x_0_xxxzz_yzzzzzz[k];

                g_x_0_xxxzzz_zzzzzz[k] = -g_x_0_xxxzz_zzzzzz[k] * ab_z + g_x_0_xxxzz_zzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_xxyyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 307 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyy_xxxxxx, g_x_0_xxyyy_xxxxxxy, g_x_0_xxyyy_xxxxxy, g_x_0_xxyyy_xxxxxyy, g_x_0_xxyyy_xxxxxyz, g_x_0_xxyyy_xxxxxz, g_x_0_xxyyy_xxxxyy, g_x_0_xxyyy_xxxxyyy, g_x_0_xxyyy_xxxxyyz, g_x_0_xxyyy_xxxxyz, g_x_0_xxyyy_xxxxyzz, g_x_0_xxyyy_xxxxzz, g_x_0_xxyyy_xxxyyy, g_x_0_xxyyy_xxxyyyy, g_x_0_xxyyy_xxxyyyz, g_x_0_xxyyy_xxxyyz, g_x_0_xxyyy_xxxyyzz, g_x_0_xxyyy_xxxyzz, g_x_0_xxyyy_xxxyzzz, g_x_0_xxyyy_xxxzzz, g_x_0_xxyyy_xxyyyy, g_x_0_xxyyy_xxyyyyy, g_x_0_xxyyy_xxyyyyz, g_x_0_xxyyy_xxyyyz, g_x_0_xxyyy_xxyyyzz, g_x_0_xxyyy_xxyyzz, g_x_0_xxyyy_xxyyzzz, g_x_0_xxyyy_xxyzzz, g_x_0_xxyyy_xxyzzzz, g_x_0_xxyyy_xxzzzz, g_x_0_xxyyy_xyyyyy, g_x_0_xxyyy_xyyyyyy, g_x_0_xxyyy_xyyyyyz, g_x_0_xxyyy_xyyyyz, g_x_0_xxyyy_xyyyyzz, g_x_0_xxyyy_xyyyzz, g_x_0_xxyyy_xyyyzzz, g_x_0_xxyyy_xyyzzz, g_x_0_xxyyy_xyyzzzz, g_x_0_xxyyy_xyzzzz, g_x_0_xxyyy_xyzzzzz, g_x_0_xxyyy_xzzzzz, g_x_0_xxyyy_yyyyyy, g_x_0_xxyyy_yyyyyyy, g_x_0_xxyyy_yyyyyyz, g_x_0_xxyyy_yyyyyz, g_x_0_xxyyy_yyyyyzz, g_x_0_xxyyy_yyyyzz, g_x_0_xxyyy_yyyyzzz, g_x_0_xxyyy_yyyzzz, g_x_0_xxyyy_yyyzzzz, g_x_0_xxyyy_yyzzzz, g_x_0_xxyyy_yyzzzzz, g_x_0_xxyyy_yzzzzz, g_x_0_xxyyy_yzzzzzz, g_x_0_xxyyy_zzzzzz, g_x_0_xxyyyy_xxxxxx, g_x_0_xxyyyy_xxxxxy, g_x_0_xxyyyy_xxxxxz, g_x_0_xxyyyy_xxxxyy, g_x_0_xxyyyy_xxxxyz, g_x_0_xxyyyy_xxxxzz, g_x_0_xxyyyy_xxxyyy, g_x_0_xxyyyy_xxxyyz, g_x_0_xxyyyy_xxxyzz, g_x_0_xxyyyy_xxxzzz, g_x_0_xxyyyy_xxyyyy, g_x_0_xxyyyy_xxyyyz, g_x_0_xxyyyy_xxyyzz, g_x_0_xxyyyy_xxyzzz, g_x_0_xxyyyy_xxzzzz, g_x_0_xxyyyy_xyyyyy, g_x_0_xxyyyy_xyyyyz, g_x_0_xxyyyy_xyyyzz, g_x_0_xxyyyy_xyyzzz, g_x_0_xxyyyy_xyzzzz, g_x_0_xxyyyy_xzzzzz, g_x_0_xxyyyy_yyyyyy, g_x_0_xxyyyy_yyyyyz, g_x_0_xxyyyy_yyyyzz, g_x_0_xxyyyy_yyyzzz, g_x_0_xxyyyy_yyzzzz, g_x_0_xxyyyy_yzzzzz, g_x_0_xxyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_xxxxxx[k] = -g_x_0_xxyyy_xxxxxx[k] * ab_y + g_x_0_xxyyy_xxxxxxy[k];

                g_x_0_xxyyyy_xxxxxy[k] = -g_x_0_xxyyy_xxxxxy[k] * ab_y + g_x_0_xxyyy_xxxxxyy[k];

                g_x_0_xxyyyy_xxxxxz[k] = -g_x_0_xxyyy_xxxxxz[k] * ab_y + g_x_0_xxyyy_xxxxxyz[k];

                g_x_0_xxyyyy_xxxxyy[k] = -g_x_0_xxyyy_xxxxyy[k] * ab_y + g_x_0_xxyyy_xxxxyyy[k];

                g_x_0_xxyyyy_xxxxyz[k] = -g_x_0_xxyyy_xxxxyz[k] * ab_y + g_x_0_xxyyy_xxxxyyz[k];

                g_x_0_xxyyyy_xxxxzz[k] = -g_x_0_xxyyy_xxxxzz[k] * ab_y + g_x_0_xxyyy_xxxxyzz[k];

                g_x_0_xxyyyy_xxxyyy[k] = -g_x_0_xxyyy_xxxyyy[k] * ab_y + g_x_0_xxyyy_xxxyyyy[k];

                g_x_0_xxyyyy_xxxyyz[k] = -g_x_0_xxyyy_xxxyyz[k] * ab_y + g_x_0_xxyyy_xxxyyyz[k];

                g_x_0_xxyyyy_xxxyzz[k] = -g_x_0_xxyyy_xxxyzz[k] * ab_y + g_x_0_xxyyy_xxxyyzz[k];

                g_x_0_xxyyyy_xxxzzz[k] = -g_x_0_xxyyy_xxxzzz[k] * ab_y + g_x_0_xxyyy_xxxyzzz[k];

                g_x_0_xxyyyy_xxyyyy[k] = -g_x_0_xxyyy_xxyyyy[k] * ab_y + g_x_0_xxyyy_xxyyyyy[k];

                g_x_0_xxyyyy_xxyyyz[k] = -g_x_0_xxyyy_xxyyyz[k] * ab_y + g_x_0_xxyyy_xxyyyyz[k];

                g_x_0_xxyyyy_xxyyzz[k] = -g_x_0_xxyyy_xxyyzz[k] * ab_y + g_x_0_xxyyy_xxyyyzz[k];

                g_x_0_xxyyyy_xxyzzz[k] = -g_x_0_xxyyy_xxyzzz[k] * ab_y + g_x_0_xxyyy_xxyyzzz[k];

                g_x_0_xxyyyy_xxzzzz[k] = -g_x_0_xxyyy_xxzzzz[k] * ab_y + g_x_0_xxyyy_xxyzzzz[k];

                g_x_0_xxyyyy_xyyyyy[k] = -g_x_0_xxyyy_xyyyyy[k] * ab_y + g_x_0_xxyyy_xyyyyyy[k];

                g_x_0_xxyyyy_xyyyyz[k] = -g_x_0_xxyyy_xyyyyz[k] * ab_y + g_x_0_xxyyy_xyyyyyz[k];

                g_x_0_xxyyyy_xyyyzz[k] = -g_x_0_xxyyy_xyyyzz[k] * ab_y + g_x_0_xxyyy_xyyyyzz[k];

                g_x_0_xxyyyy_xyyzzz[k] = -g_x_0_xxyyy_xyyzzz[k] * ab_y + g_x_0_xxyyy_xyyyzzz[k];

                g_x_0_xxyyyy_xyzzzz[k] = -g_x_0_xxyyy_xyzzzz[k] * ab_y + g_x_0_xxyyy_xyyzzzz[k];

                g_x_0_xxyyyy_xzzzzz[k] = -g_x_0_xxyyy_xzzzzz[k] * ab_y + g_x_0_xxyyy_xyzzzzz[k];

                g_x_0_xxyyyy_yyyyyy[k] = -g_x_0_xxyyy_yyyyyy[k] * ab_y + g_x_0_xxyyy_yyyyyyy[k];

                g_x_0_xxyyyy_yyyyyz[k] = -g_x_0_xxyyy_yyyyyz[k] * ab_y + g_x_0_xxyyy_yyyyyyz[k];

                g_x_0_xxyyyy_yyyyzz[k] = -g_x_0_xxyyy_yyyyzz[k] * ab_y + g_x_0_xxyyy_yyyyyzz[k];

                g_x_0_xxyyyy_yyyzzz[k] = -g_x_0_xxyyy_yyyzzz[k] * ab_y + g_x_0_xxyyy_yyyyzzz[k];

                g_x_0_xxyyyy_yyzzzz[k] = -g_x_0_xxyyy_yyzzzz[k] * ab_y + g_x_0_xxyyy_yyyzzzz[k];

                g_x_0_xxyyyy_yzzzzz[k] = -g_x_0_xxyyy_yzzzzz[k] * ab_y + g_x_0_xxyyy_yyzzzzz[k];

                g_x_0_xxyyyy_zzzzzz[k] = -g_x_0_xxyyy_zzzzzz[k] * ab_y + g_x_0_xxyyy_yzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_xxyyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyyz_xxxxxx, g_x_0_xxyyyz_xxxxxy, g_x_0_xxyyyz_xxxxxz, g_x_0_xxyyyz_xxxxyy, g_x_0_xxyyyz_xxxxyz, g_x_0_xxyyyz_xxxxzz, g_x_0_xxyyyz_xxxyyy, g_x_0_xxyyyz_xxxyyz, g_x_0_xxyyyz_xxxyzz, g_x_0_xxyyyz_xxxzzz, g_x_0_xxyyyz_xxyyyy, g_x_0_xxyyyz_xxyyyz, g_x_0_xxyyyz_xxyyzz, g_x_0_xxyyyz_xxyzzz, g_x_0_xxyyyz_xxzzzz, g_x_0_xxyyyz_xyyyyy, g_x_0_xxyyyz_xyyyyz, g_x_0_xxyyyz_xyyyzz, g_x_0_xxyyyz_xyyzzz, g_x_0_xxyyyz_xyzzzz, g_x_0_xxyyyz_xzzzzz, g_x_0_xxyyyz_yyyyyy, g_x_0_xxyyyz_yyyyyz, g_x_0_xxyyyz_yyyyzz, g_x_0_xxyyyz_yyyzzz, g_x_0_xxyyyz_yyzzzz, g_x_0_xxyyyz_yzzzzz, g_x_0_xxyyyz_zzzzzz, g_x_0_xxyyz_xxxxxx, g_x_0_xxyyz_xxxxxxy, g_x_0_xxyyz_xxxxxy, g_x_0_xxyyz_xxxxxyy, g_x_0_xxyyz_xxxxxyz, g_x_0_xxyyz_xxxxxz, g_x_0_xxyyz_xxxxyy, g_x_0_xxyyz_xxxxyyy, g_x_0_xxyyz_xxxxyyz, g_x_0_xxyyz_xxxxyz, g_x_0_xxyyz_xxxxyzz, g_x_0_xxyyz_xxxxzz, g_x_0_xxyyz_xxxyyy, g_x_0_xxyyz_xxxyyyy, g_x_0_xxyyz_xxxyyyz, g_x_0_xxyyz_xxxyyz, g_x_0_xxyyz_xxxyyzz, g_x_0_xxyyz_xxxyzz, g_x_0_xxyyz_xxxyzzz, g_x_0_xxyyz_xxxzzz, g_x_0_xxyyz_xxyyyy, g_x_0_xxyyz_xxyyyyy, g_x_0_xxyyz_xxyyyyz, g_x_0_xxyyz_xxyyyz, g_x_0_xxyyz_xxyyyzz, g_x_0_xxyyz_xxyyzz, g_x_0_xxyyz_xxyyzzz, g_x_0_xxyyz_xxyzzz, g_x_0_xxyyz_xxyzzzz, g_x_0_xxyyz_xxzzzz, g_x_0_xxyyz_xyyyyy, g_x_0_xxyyz_xyyyyyy, g_x_0_xxyyz_xyyyyyz, g_x_0_xxyyz_xyyyyz, g_x_0_xxyyz_xyyyyzz, g_x_0_xxyyz_xyyyzz, g_x_0_xxyyz_xyyyzzz, g_x_0_xxyyz_xyyzzz, g_x_0_xxyyz_xyyzzzz, g_x_0_xxyyz_xyzzzz, g_x_0_xxyyz_xyzzzzz, g_x_0_xxyyz_xzzzzz, g_x_0_xxyyz_yyyyyy, g_x_0_xxyyz_yyyyyyy, g_x_0_xxyyz_yyyyyyz, g_x_0_xxyyz_yyyyyz, g_x_0_xxyyz_yyyyyzz, g_x_0_xxyyz_yyyyzz, g_x_0_xxyyz_yyyyzzz, g_x_0_xxyyz_yyyzzz, g_x_0_xxyyz_yyyzzzz, g_x_0_xxyyz_yyzzzz, g_x_0_xxyyz_yyzzzzz, g_x_0_xxyyz_yzzzzz, g_x_0_xxyyz_yzzzzzz, g_x_0_xxyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_xxxxxx[k] = -g_x_0_xxyyz_xxxxxx[k] * ab_y + g_x_0_xxyyz_xxxxxxy[k];

                g_x_0_xxyyyz_xxxxxy[k] = -g_x_0_xxyyz_xxxxxy[k] * ab_y + g_x_0_xxyyz_xxxxxyy[k];

                g_x_0_xxyyyz_xxxxxz[k] = -g_x_0_xxyyz_xxxxxz[k] * ab_y + g_x_0_xxyyz_xxxxxyz[k];

                g_x_0_xxyyyz_xxxxyy[k] = -g_x_0_xxyyz_xxxxyy[k] * ab_y + g_x_0_xxyyz_xxxxyyy[k];

                g_x_0_xxyyyz_xxxxyz[k] = -g_x_0_xxyyz_xxxxyz[k] * ab_y + g_x_0_xxyyz_xxxxyyz[k];

                g_x_0_xxyyyz_xxxxzz[k] = -g_x_0_xxyyz_xxxxzz[k] * ab_y + g_x_0_xxyyz_xxxxyzz[k];

                g_x_0_xxyyyz_xxxyyy[k] = -g_x_0_xxyyz_xxxyyy[k] * ab_y + g_x_0_xxyyz_xxxyyyy[k];

                g_x_0_xxyyyz_xxxyyz[k] = -g_x_0_xxyyz_xxxyyz[k] * ab_y + g_x_0_xxyyz_xxxyyyz[k];

                g_x_0_xxyyyz_xxxyzz[k] = -g_x_0_xxyyz_xxxyzz[k] * ab_y + g_x_0_xxyyz_xxxyyzz[k];

                g_x_0_xxyyyz_xxxzzz[k] = -g_x_0_xxyyz_xxxzzz[k] * ab_y + g_x_0_xxyyz_xxxyzzz[k];

                g_x_0_xxyyyz_xxyyyy[k] = -g_x_0_xxyyz_xxyyyy[k] * ab_y + g_x_0_xxyyz_xxyyyyy[k];

                g_x_0_xxyyyz_xxyyyz[k] = -g_x_0_xxyyz_xxyyyz[k] * ab_y + g_x_0_xxyyz_xxyyyyz[k];

                g_x_0_xxyyyz_xxyyzz[k] = -g_x_0_xxyyz_xxyyzz[k] * ab_y + g_x_0_xxyyz_xxyyyzz[k];

                g_x_0_xxyyyz_xxyzzz[k] = -g_x_0_xxyyz_xxyzzz[k] * ab_y + g_x_0_xxyyz_xxyyzzz[k];

                g_x_0_xxyyyz_xxzzzz[k] = -g_x_0_xxyyz_xxzzzz[k] * ab_y + g_x_0_xxyyz_xxyzzzz[k];

                g_x_0_xxyyyz_xyyyyy[k] = -g_x_0_xxyyz_xyyyyy[k] * ab_y + g_x_0_xxyyz_xyyyyyy[k];

                g_x_0_xxyyyz_xyyyyz[k] = -g_x_0_xxyyz_xyyyyz[k] * ab_y + g_x_0_xxyyz_xyyyyyz[k];

                g_x_0_xxyyyz_xyyyzz[k] = -g_x_0_xxyyz_xyyyzz[k] * ab_y + g_x_0_xxyyz_xyyyyzz[k];

                g_x_0_xxyyyz_xyyzzz[k] = -g_x_0_xxyyz_xyyzzz[k] * ab_y + g_x_0_xxyyz_xyyyzzz[k];

                g_x_0_xxyyyz_xyzzzz[k] = -g_x_0_xxyyz_xyzzzz[k] * ab_y + g_x_0_xxyyz_xyyzzzz[k];

                g_x_0_xxyyyz_xzzzzz[k] = -g_x_0_xxyyz_xzzzzz[k] * ab_y + g_x_0_xxyyz_xyzzzzz[k];

                g_x_0_xxyyyz_yyyyyy[k] = -g_x_0_xxyyz_yyyyyy[k] * ab_y + g_x_0_xxyyz_yyyyyyy[k];

                g_x_0_xxyyyz_yyyyyz[k] = -g_x_0_xxyyz_yyyyyz[k] * ab_y + g_x_0_xxyyz_yyyyyyz[k];

                g_x_0_xxyyyz_yyyyzz[k] = -g_x_0_xxyyz_yyyyzz[k] * ab_y + g_x_0_xxyyz_yyyyyzz[k];

                g_x_0_xxyyyz_yyyzzz[k] = -g_x_0_xxyyz_yyyzzz[k] * ab_y + g_x_0_xxyyz_yyyyzzz[k];

                g_x_0_xxyyyz_yyzzzz[k] = -g_x_0_xxyyz_yyzzzz[k] * ab_y + g_x_0_xxyyz_yyyzzzz[k];

                g_x_0_xxyyyz_yzzzzz[k] = -g_x_0_xxyyz_yzzzzz[k] * ab_y + g_x_0_xxyyz_yyzzzzz[k];

                g_x_0_xxyyyz_zzzzzz[k] = -g_x_0_xxyyz_zzzzzz[k] * ab_y + g_x_0_xxyyz_yzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 356 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 360 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 362 * ccomps * dcomps);

            auto g_x_0_xxyyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 363 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyzz_xxxxxx, g_x_0_xxyyzz_xxxxxy, g_x_0_xxyyzz_xxxxxz, g_x_0_xxyyzz_xxxxyy, g_x_0_xxyyzz_xxxxyz, g_x_0_xxyyzz_xxxxzz, g_x_0_xxyyzz_xxxyyy, g_x_0_xxyyzz_xxxyyz, g_x_0_xxyyzz_xxxyzz, g_x_0_xxyyzz_xxxzzz, g_x_0_xxyyzz_xxyyyy, g_x_0_xxyyzz_xxyyyz, g_x_0_xxyyzz_xxyyzz, g_x_0_xxyyzz_xxyzzz, g_x_0_xxyyzz_xxzzzz, g_x_0_xxyyzz_xyyyyy, g_x_0_xxyyzz_xyyyyz, g_x_0_xxyyzz_xyyyzz, g_x_0_xxyyzz_xyyzzz, g_x_0_xxyyzz_xyzzzz, g_x_0_xxyyzz_xzzzzz, g_x_0_xxyyzz_yyyyyy, g_x_0_xxyyzz_yyyyyz, g_x_0_xxyyzz_yyyyzz, g_x_0_xxyyzz_yyyzzz, g_x_0_xxyyzz_yyzzzz, g_x_0_xxyyzz_yzzzzz, g_x_0_xxyyzz_zzzzzz, g_x_0_xxyzz_xxxxxx, g_x_0_xxyzz_xxxxxxy, g_x_0_xxyzz_xxxxxy, g_x_0_xxyzz_xxxxxyy, g_x_0_xxyzz_xxxxxyz, g_x_0_xxyzz_xxxxxz, g_x_0_xxyzz_xxxxyy, g_x_0_xxyzz_xxxxyyy, g_x_0_xxyzz_xxxxyyz, g_x_0_xxyzz_xxxxyz, g_x_0_xxyzz_xxxxyzz, g_x_0_xxyzz_xxxxzz, g_x_0_xxyzz_xxxyyy, g_x_0_xxyzz_xxxyyyy, g_x_0_xxyzz_xxxyyyz, g_x_0_xxyzz_xxxyyz, g_x_0_xxyzz_xxxyyzz, g_x_0_xxyzz_xxxyzz, g_x_0_xxyzz_xxxyzzz, g_x_0_xxyzz_xxxzzz, g_x_0_xxyzz_xxyyyy, g_x_0_xxyzz_xxyyyyy, g_x_0_xxyzz_xxyyyyz, g_x_0_xxyzz_xxyyyz, g_x_0_xxyzz_xxyyyzz, g_x_0_xxyzz_xxyyzz, g_x_0_xxyzz_xxyyzzz, g_x_0_xxyzz_xxyzzz, g_x_0_xxyzz_xxyzzzz, g_x_0_xxyzz_xxzzzz, g_x_0_xxyzz_xyyyyy, g_x_0_xxyzz_xyyyyyy, g_x_0_xxyzz_xyyyyyz, g_x_0_xxyzz_xyyyyz, g_x_0_xxyzz_xyyyyzz, g_x_0_xxyzz_xyyyzz, g_x_0_xxyzz_xyyyzzz, g_x_0_xxyzz_xyyzzz, g_x_0_xxyzz_xyyzzzz, g_x_0_xxyzz_xyzzzz, g_x_0_xxyzz_xyzzzzz, g_x_0_xxyzz_xzzzzz, g_x_0_xxyzz_yyyyyy, g_x_0_xxyzz_yyyyyyy, g_x_0_xxyzz_yyyyyyz, g_x_0_xxyzz_yyyyyz, g_x_0_xxyzz_yyyyyzz, g_x_0_xxyzz_yyyyzz, g_x_0_xxyzz_yyyyzzz, g_x_0_xxyzz_yyyzzz, g_x_0_xxyzz_yyyzzzz, g_x_0_xxyzz_yyzzzz, g_x_0_xxyzz_yyzzzzz, g_x_0_xxyzz_yzzzzz, g_x_0_xxyzz_yzzzzzz, g_x_0_xxyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_xxxxxx[k] = -g_x_0_xxyzz_xxxxxx[k] * ab_y + g_x_0_xxyzz_xxxxxxy[k];

                g_x_0_xxyyzz_xxxxxy[k] = -g_x_0_xxyzz_xxxxxy[k] * ab_y + g_x_0_xxyzz_xxxxxyy[k];

                g_x_0_xxyyzz_xxxxxz[k] = -g_x_0_xxyzz_xxxxxz[k] * ab_y + g_x_0_xxyzz_xxxxxyz[k];

                g_x_0_xxyyzz_xxxxyy[k] = -g_x_0_xxyzz_xxxxyy[k] * ab_y + g_x_0_xxyzz_xxxxyyy[k];

                g_x_0_xxyyzz_xxxxyz[k] = -g_x_0_xxyzz_xxxxyz[k] * ab_y + g_x_0_xxyzz_xxxxyyz[k];

                g_x_0_xxyyzz_xxxxzz[k] = -g_x_0_xxyzz_xxxxzz[k] * ab_y + g_x_0_xxyzz_xxxxyzz[k];

                g_x_0_xxyyzz_xxxyyy[k] = -g_x_0_xxyzz_xxxyyy[k] * ab_y + g_x_0_xxyzz_xxxyyyy[k];

                g_x_0_xxyyzz_xxxyyz[k] = -g_x_0_xxyzz_xxxyyz[k] * ab_y + g_x_0_xxyzz_xxxyyyz[k];

                g_x_0_xxyyzz_xxxyzz[k] = -g_x_0_xxyzz_xxxyzz[k] * ab_y + g_x_0_xxyzz_xxxyyzz[k];

                g_x_0_xxyyzz_xxxzzz[k] = -g_x_0_xxyzz_xxxzzz[k] * ab_y + g_x_0_xxyzz_xxxyzzz[k];

                g_x_0_xxyyzz_xxyyyy[k] = -g_x_0_xxyzz_xxyyyy[k] * ab_y + g_x_0_xxyzz_xxyyyyy[k];

                g_x_0_xxyyzz_xxyyyz[k] = -g_x_0_xxyzz_xxyyyz[k] * ab_y + g_x_0_xxyzz_xxyyyyz[k];

                g_x_0_xxyyzz_xxyyzz[k] = -g_x_0_xxyzz_xxyyzz[k] * ab_y + g_x_0_xxyzz_xxyyyzz[k];

                g_x_0_xxyyzz_xxyzzz[k] = -g_x_0_xxyzz_xxyzzz[k] * ab_y + g_x_0_xxyzz_xxyyzzz[k];

                g_x_0_xxyyzz_xxzzzz[k] = -g_x_0_xxyzz_xxzzzz[k] * ab_y + g_x_0_xxyzz_xxyzzzz[k];

                g_x_0_xxyyzz_xyyyyy[k] = -g_x_0_xxyzz_xyyyyy[k] * ab_y + g_x_0_xxyzz_xyyyyyy[k];

                g_x_0_xxyyzz_xyyyyz[k] = -g_x_0_xxyzz_xyyyyz[k] * ab_y + g_x_0_xxyzz_xyyyyyz[k];

                g_x_0_xxyyzz_xyyyzz[k] = -g_x_0_xxyzz_xyyyzz[k] * ab_y + g_x_0_xxyzz_xyyyyzz[k];

                g_x_0_xxyyzz_xyyzzz[k] = -g_x_0_xxyzz_xyyzzz[k] * ab_y + g_x_0_xxyzz_xyyyzzz[k];

                g_x_0_xxyyzz_xyzzzz[k] = -g_x_0_xxyzz_xyzzzz[k] * ab_y + g_x_0_xxyzz_xyyzzzz[k];

                g_x_0_xxyyzz_xzzzzz[k] = -g_x_0_xxyzz_xzzzzz[k] * ab_y + g_x_0_xxyzz_xyzzzzz[k];

                g_x_0_xxyyzz_yyyyyy[k] = -g_x_0_xxyzz_yyyyyy[k] * ab_y + g_x_0_xxyzz_yyyyyyy[k];

                g_x_0_xxyyzz_yyyyyz[k] = -g_x_0_xxyzz_yyyyyz[k] * ab_y + g_x_0_xxyzz_yyyyyyz[k];

                g_x_0_xxyyzz_yyyyzz[k] = -g_x_0_xxyzz_yyyyzz[k] * ab_y + g_x_0_xxyzz_yyyyyzz[k];

                g_x_0_xxyyzz_yyyzzz[k] = -g_x_0_xxyzz_yyyzzz[k] * ab_y + g_x_0_xxyzz_yyyyzzz[k];

                g_x_0_xxyyzz_yyzzzz[k] = -g_x_0_xxyzz_yyzzzz[k] * ab_y + g_x_0_xxyzz_yyyzzzz[k];

                g_x_0_xxyyzz_yzzzzz[k] = -g_x_0_xxyzz_yzzzzz[k] * ab_y + g_x_0_xxyzz_yyzzzzz[k];

                g_x_0_xxyyzz_zzzzzz[k] = -g_x_0_xxyzz_zzzzzz[k] * ab_y + g_x_0_xxyzz_yzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 364 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 365 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 366 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 369 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 373 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 374 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 378 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 380 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 384 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 387 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_xxyzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 391 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyzzz_xxxxxx, g_x_0_xxyzzz_xxxxxy, g_x_0_xxyzzz_xxxxxz, g_x_0_xxyzzz_xxxxyy, g_x_0_xxyzzz_xxxxyz, g_x_0_xxyzzz_xxxxzz, g_x_0_xxyzzz_xxxyyy, g_x_0_xxyzzz_xxxyyz, g_x_0_xxyzzz_xxxyzz, g_x_0_xxyzzz_xxxzzz, g_x_0_xxyzzz_xxyyyy, g_x_0_xxyzzz_xxyyyz, g_x_0_xxyzzz_xxyyzz, g_x_0_xxyzzz_xxyzzz, g_x_0_xxyzzz_xxzzzz, g_x_0_xxyzzz_xyyyyy, g_x_0_xxyzzz_xyyyyz, g_x_0_xxyzzz_xyyyzz, g_x_0_xxyzzz_xyyzzz, g_x_0_xxyzzz_xyzzzz, g_x_0_xxyzzz_xzzzzz, g_x_0_xxyzzz_yyyyyy, g_x_0_xxyzzz_yyyyyz, g_x_0_xxyzzz_yyyyzz, g_x_0_xxyzzz_yyyzzz, g_x_0_xxyzzz_yyzzzz, g_x_0_xxyzzz_yzzzzz, g_x_0_xxyzzz_zzzzzz, g_x_0_xxzzz_xxxxxx, g_x_0_xxzzz_xxxxxxy, g_x_0_xxzzz_xxxxxy, g_x_0_xxzzz_xxxxxyy, g_x_0_xxzzz_xxxxxyz, g_x_0_xxzzz_xxxxxz, g_x_0_xxzzz_xxxxyy, g_x_0_xxzzz_xxxxyyy, g_x_0_xxzzz_xxxxyyz, g_x_0_xxzzz_xxxxyz, g_x_0_xxzzz_xxxxyzz, g_x_0_xxzzz_xxxxzz, g_x_0_xxzzz_xxxyyy, g_x_0_xxzzz_xxxyyyy, g_x_0_xxzzz_xxxyyyz, g_x_0_xxzzz_xxxyyz, g_x_0_xxzzz_xxxyyzz, g_x_0_xxzzz_xxxyzz, g_x_0_xxzzz_xxxyzzz, g_x_0_xxzzz_xxxzzz, g_x_0_xxzzz_xxyyyy, g_x_0_xxzzz_xxyyyyy, g_x_0_xxzzz_xxyyyyz, g_x_0_xxzzz_xxyyyz, g_x_0_xxzzz_xxyyyzz, g_x_0_xxzzz_xxyyzz, g_x_0_xxzzz_xxyyzzz, g_x_0_xxzzz_xxyzzz, g_x_0_xxzzz_xxyzzzz, g_x_0_xxzzz_xxzzzz, g_x_0_xxzzz_xyyyyy, g_x_0_xxzzz_xyyyyyy, g_x_0_xxzzz_xyyyyyz, g_x_0_xxzzz_xyyyyz, g_x_0_xxzzz_xyyyyzz, g_x_0_xxzzz_xyyyzz, g_x_0_xxzzz_xyyyzzz, g_x_0_xxzzz_xyyzzz, g_x_0_xxzzz_xyyzzzz, g_x_0_xxzzz_xyzzzz, g_x_0_xxzzz_xyzzzzz, g_x_0_xxzzz_xzzzzz, g_x_0_xxzzz_yyyyyy, g_x_0_xxzzz_yyyyyyy, g_x_0_xxzzz_yyyyyyz, g_x_0_xxzzz_yyyyyz, g_x_0_xxzzz_yyyyyzz, g_x_0_xxzzz_yyyyzz, g_x_0_xxzzz_yyyyzzz, g_x_0_xxzzz_yyyzzz, g_x_0_xxzzz_yyyzzzz, g_x_0_xxzzz_yyzzzz, g_x_0_xxzzz_yyzzzzz, g_x_0_xxzzz_yzzzzz, g_x_0_xxzzz_yzzzzzz, g_x_0_xxzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_xxxxxx[k] = -g_x_0_xxzzz_xxxxxx[k] * ab_y + g_x_0_xxzzz_xxxxxxy[k];

                g_x_0_xxyzzz_xxxxxy[k] = -g_x_0_xxzzz_xxxxxy[k] * ab_y + g_x_0_xxzzz_xxxxxyy[k];

                g_x_0_xxyzzz_xxxxxz[k] = -g_x_0_xxzzz_xxxxxz[k] * ab_y + g_x_0_xxzzz_xxxxxyz[k];

                g_x_0_xxyzzz_xxxxyy[k] = -g_x_0_xxzzz_xxxxyy[k] * ab_y + g_x_0_xxzzz_xxxxyyy[k];

                g_x_0_xxyzzz_xxxxyz[k] = -g_x_0_xxzzz_xxxxyz[k] * ab_y + g_x_0_xxzzz_xxxxyyz[k];

                g_x_0_xxyzzz_xxxxzz[k] = -g_x_0_xxzzz_xxxxzz[k] * ab_y + g_x_0_xxzzz_xxxxyzz[k];

                g_x_0_xxyzzz_xxxyyy[k] = -g_x_0_xxzzz_xxxyyy[k] * ab_y + g_x_0_xxzzz_xxxyyyy[k];

                g_x_0_xxyzzz_xxxyyz[k] = -g_x_0_xxzzz_xxxyyz[k] * ab_y + g_x_0_xxzzz_xxxyyyz[k];

                g_x_0_xxyzzz_xxxyzz[k] = -g_x_0_xxzzz_xxxyzz[k] * ab_y + g_x_0_xxzzz_xxxyyzz[k];

                g_x_0_xxyzzz_xxxzzz[k] = -g_x_0_xxzzz_xxxzzz[k] * ab_y + g_x_0_xxzzz_xxxyzzz[k];

                g_x_0_xxyzzz_xxyyyy[k] = -g_x_0_xxzzz_xxyyyy[k] * ab_y + g_x_0_xxzzz_xxyyyyy[k];

                g_x_0_xxyzzz_xxyyyz[k] = -g_x_0_xxzzz_xxyyyz[k] * ab_y + g_x_0_xxzzz_xxyyyyz[k];

                g_x_0_xxyzzz_xxyyzz[k] = -g_x_0_xxzzz_xxyyzz[k] * ab_y + g_x_0_xxzzz_xxyyyzz[k];

                g_x_0_xxyzzz_xxyzzz[k] = -g_x_0_xxzzz_xxyzzz[k] * ab_y + g_x_0_xxzzz_xxyyzzz[k];

                g_x_0_xxyzzz_xxzzzz[k] = -g_x_0_xxzzz_xxzzzz[k] * ab_y + g_x_0_xxzzz_xxyzzzz[k];

                g_x_0_xxyzzz_xyyyyy[k] = -g_x_0_xxzzz_xyyyyy[k] * ab_y + g_x_0_xxzzz_xyyyyyy[k];

                g_x_0_xxyzzz_xyyyyz[k] = -g_x_0_xxzzz_xyyyyz[k] * ab_y + g_x_0_xxzzz_xyyyyyz[k];

                g_x_0_xxyzzz_xyyyzz[k] = -g_x_0_xxzzz_xyyyzz[k] * ab_y + g_x_0_xxzzz_xyyyyzz[k];

                g_x_0_xxyzzz_xyyzzz[k] = -g_x_0_xxzzz_xyyzzz[k] * ab_y + g_x_0_xxzzz_xyyyzzz[k];

                g_x_0_xxyzzz_xyzzzz[k] = -g_x_0_xxzzz_xyzzzz[k] * ab_y + g_x_0_xxzzz_xyyzzzz[k];

                g_x_0_xxyzzz_xzzzzz[k] = -g_x_0_xxzzz_xzzzzz[k] * ab_y + g_x_0_xxzzz_xyzzzzz[k];

                g_x_0_xxyzzz_yyyyyy[k] = -g_x_0_xxzzz_yyyyyy[k] * ab_y + g_x_0_xxzzz_yyyyyyy[k];

                g_x_0_xxyzzz_yyyyyz[k] = -g_x_0_xxzzz_yyyyyz[k] * ab_y + g_x_0_xxzzz_yyyyyyz[k];

                g_x_0_xxyzzz_yyyyzz[k] = -g_x_0_xxzzz_yyyyzz[k] * ab_y + g_x_0_xxzzz_yyyyyzz[k];

                g_x_0_xxyzzz_yyyzzz[k] = -g_x_0_xxzzz_yyyzzz[k] * ab_y + g_x_0_xxzzz_yyyyzzz[k];

                g_x_0_xxyzzz_yyzzzz[k] = -g_x_0_xxzzz_yyzzzz[k] * ab_y + g_x_0_xxzzz_yyyzzzz[k];

                g_x_0_xxyzzz_yzzzzz[k] = -g_x_0_xxzzz_yzzzzz[k] * ab_y + g_x_0_xxzzz_yyzzzzz[k];

                g_x_0_xxyzzz_zzzzzz[k] = -g_x_0_xxzzz_zzzzzz[k] * ab_y + g_x_0_xxzzz_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 392 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 395 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 396 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 398 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 401 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 404 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 405 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 410 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 416 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_xxzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxzzz_xxxxxx, g_x_0_xxzzz_xxxxxxz, g_x_0_xxzzz_xxxxxy, g_x_0_xxzzz_xxxxxyz, g_x_0_xxzzz_xxxxxz, g_x_0_xxzzz_xxxxxzz, g_x_0_xxzzz_xxxxyy, g_x_0_xxzzz_xxxxyyz, g_x_0_xxzzz_xxxxyz, g_x_0_xxzzz_xxxxyzz, g_x_0_xxzzz_xxxxzz, g_x_0_xxzzz_xxxxzzz, g_x_0_xxzzz_xxxyyy, g_x_0_xxzzz_xxxyyyz, g_x_0_xxzzz_xxxyyz, g_x_0_xxzzz_xxxyyzz, g_x_0_xxzzz_xxxyzz, g_x_0_xxzzz_xxxyzzz, g_x_0_xxzzz_xxxzzz, g_x_0_xxzzz_xxxzzzz, g_x_0_xxzzz_xxyyyy, g_x_0_xxzzz_xxyyyyz, g_x_0_xxzzz_xxyyyz, g_x_0_xxzzz_xxyyyzz, g_x_0_xxzzz_xxyyzz, g_x_0_xxzzz_xxyyzzz, g_x_0_xxzzz_xxyzzz, g_x_0_xxzzz_xxyzzzz, g_x_0_xxzzz_xxzzzz, g_x_0_xxzzz_xxzzzzz, g_x_0_xxzzz_xyyyyy, g_x_0_xxzzz_xyyyyyz, g_x_0_xxzzz_xyyyyz, g_x_0_xxzzz_xyyyyzz, g_x_0_xxzzz_xyyyzz, g_x_0_xxzzz_xyyyzzz, g_x_0_xxzzz_xyyzzz, g_x_0_xxzzz_xyyzzzz, g_x_0_xxzzz_xyzzzz, g_x_0_xxzzz_xyzzzzz, g_x_0_xxzzz_xzzzzz, g_x_0_xxzzz_xzzzzzz, g_x_0_xxzzz_yyyyyy, g_x_0_xxzzz_yyyyyyz, g_x_0_xxzzz_yyyyyz, g_x_0_xxzzz_yyyyyzz, g_x_0_xxzzz_yyyyzz, g_x_0_xxzzz_yyyyzzz, g_x_0_xxzzz_yyyzzz, g_x_0_xxzzz_yyyzzzz, g_x_0_xxzzz_yyzzzz, g_x_0_xxzzz_yyzzzzz, g_x_0_xxzzz_yzzzzz, g_x_0_xxzzz_yzzzzzz, g_x_0_xxzzz_zzzzzz, g_x_0_xxzzz_zzzzzzz, g_x_0_xxzzzz_xxxxxx, g_x_0_xxzzzz_xxxxxy, g_x_0_xxzzzz_xxxxxz, g_x_0_xxzzzz_xxxxyy, g_x_0_xxzzzz_xxxxyz, g_x_0_xxzzzz_xxxxzz, g_x_0_xxzzzz_xxxyyy, g_x_0_xxzzzz_xxxyyz, g_x_0_xxzzzz_xxxyzz, g_x_0_xxzzzz_xxxzzz, g_x_0_xxzzzz_xxyyyy, g_x_0_xxzzzz_xxyyyz, g_x_0_xxzzzz_xxyyzz, g_x_0_xxzzzz_xxyzzz, g_x_0_xxzzzz_xxzzzz, g_x_0_xxzzzz_xyyyyy, g_x_0_xxzzzz_xyyyyz, g_x_0_xxzzzz_xyyyzz, g_x_0_xxzzzz_xyyzzz, g_x_0_xxzzzz_xyzzzz, g_x_0_xxzzzz_xzzzzz, g_x_0_xxzzzz_yyyyyy, g_x_0_xxzzzz_yyyyyz, g_x_0_xxzzzz_yyyyzz, g_x_0_xxzzzz_yyyzzz, g_x_0_xxzzzz_yyzzzz, g_x_0_xxzzzz_yzzzzz, g_x_0_xxzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_xxxxxx[k] = -g_x_0_xxzzz_xxxxxx[k] * ab_z + g_x_0_xxzzz_xxxxxxz[k];

                g_x_0_xxzzzz_xxxxxy[k] = -g_x_0_xxzzz_xxxxxy[k] * ab_z + g_x_0_xxzzz_xxxxxyz[k];

                g_x_0_xxzzzz_xxxxxz[k] = -g_x_0_xxzzz_xxxxxz[k] * ab_z + g_x_0_xxzzz_xxxxxzz[k];

                g_x_0_xxzzzz_xxxxyy[k] = -g_x_0_xxzzz_xxxxyy[k] * ab_z + g_x_0_xxzzz_xxxxyyz[k];

                g_x_0_xxzzzz_xxxxyz[k] = -g_x_0_xxzzz_xxxxyz[k] * ab_z + g_x_0_xxzzz_xxxxyzz[k];

                g_x_0_xxzzzz_xxxxzz[k] = -g_x_0_xxzzz_xxxxzz[k] * ab_z + g_x_0_xxzzz_xxxxzzz[k];

                g_x_0_xxzzzz_xxxyyy[k] = -g_x_0_xxzzz_xxxyyy[k] * ab_z + g_x_0_xxzzz_xxxyyyz[k];

                g_x_0_xxzzzz_xxxyyz[k] = -g_x_0_xxzzz_xxxyyz[k] * ab_z + g_x_0_xxzzz_xxxyyzz[k];

                g_x_0_xxzzzz_xxxyzz[k] = -g_x_0_xxzzz_xxxyzz[k] * ab_z + g_x_0_xxzzz_xxxyzzz[k];

                g_x_0_xxzzzz_xxxzzz[k] = -g_x_0_xxzzz_xxxzzz[k] * ab_z + g_x_0_xxzzz_xxxzzzz[k];

                g_x_0_xxzzzz_xxyyyy[k] = -g_x_0_xxzzz_xxyyyy[k] * ab_z + g_x_0_xxzzz_xxyyyyz[k];

                g_x_0_xxzzzz_xxyyyz[k] = -g_x_0_xxzzz_xxyyyz[k] * ab_z + g_x_0_xxzzz_xxyyyzz[k];

                g_x_0_xxzzzz_xxyyzz[k] = -g_x_0_xxzzz_xxyyzz[k] * ab_z + g_x_0_xxzzz_xxyyzzz[k];

                g_x_0_xxzzzz_xxyzzz[k] = -g_x_0_xxzzz_xxyzzz[k] * ab_z + g_x_0_xxzzz_xxyzzzz[k];

                g_x_0_xxzzzz_xxzzzz[k] = -g_x_0_xxzzz_xxzzzz[k] * ab_z + g_x_0_xxzzz_xxzzzzz[k];

                g_x_0_xxzzzz_xyyyyy[k] = -g_x_0_xxzzz_xyyyyy[k] * ab_z + g_x_0_xxzzz_xyyyyyz[k];

                g_x_0_xxzzzz_xyyyyz[k] = -g_x_0_xxzzz_xyyyyz[k] * ab_z + g_x_0_xxzzz_xyyyyzz[k];

                g_x_0_xxzzzz_xyyyzz[k] = -g_x_0_xxzzz_xyyyzz[k] * ab_z + g_x_0_xxzzz_xyyyzzz[k];

                g_x_0_xxzzzz_xyyzzz[k] = -g_x_0_xxzzz_xyyzzz[k] * ab_z + g_x_0_xxzzz_xyyzzzz[k];

                g_x_0_xxzzzz_xyzzzz[k] = -g_x_0_xxzzz_xyzzzz[k] * ab_z + g_x_0_xxzzz_xyzzzzz[k];

                g_x_0_xxzzzz_xzzzzz[k] = -g_x_0_xxzzz_xzzzzz[k] * ab_z + g_x_0_xxzzz_xzzzzzz[k];

                g_x_0_xxzzzz_yyyyyy[k] = -g_x_0_xxzzz_yyyyyy[k] * ab_z + g_x_0_xxzzz_yyyyyyz[k];

                g_x_0_xxzzzz_yyyyyz[k] = -g_x_0_xxzzz_yyyyyz[k] * ab_z + g_x_0_xxzzz_yyyyyzz[k];

                g_x_0_xxzzzz_yyyyzz[k] = -g_x_0_xxzzz_yyyyzz[k] * ab_z + g_x_0_xxzzz_yyyyzzz[k];

                g_x_0_xxzzzz_yyyzzz[k] = -g_x_0_xxzzz_yyyzzz[k] * ab_z + g_x_0_xxzzz_yyyzzzz[k];

                g_x_0_xxzzzz_yyzzzz[k] = -g_x_0_xxzzz_yyzzzz[k] * ab_z + g_x_0_xxzzz_yyzzzzz[k];

                g_x_0_xxzzzz_yzzzzz[k] = -g_x_0_xxzzz_yzzzzz[k] * ab_z + g_x_0_xxzzz_yzzzzzz[k];

                g_x_0_xxzzzz_zzzzzz[k] = -g_x_0_xxzzz_zzzzzz[k] * ab_z + g_x_0_xxzzz_zzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 420 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 421 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 422 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 423 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 424 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 425 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 426 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 427 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 428 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 429 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 430 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 431 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 432 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 433 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 434 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 435 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 436 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 437 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 438 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 439 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 440 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 441 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 442 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 443 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 444 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 445 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 446 * ccomps * dcomps);

            auto g_x_0_xyyyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 447 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyy_xxxxxx, g_x_0_xyyyy_xxxxxxy, g_x_0_xyyyy_xxxxxy, g_x_0_xyyyy_xxxxxyy, g_x_0_xyyyy_xxxxxyz, g_x_0_xyyyy_xxxxxz, g_x_0_xyyyy_xxxxyy, g_x_0_xyyyy_xxxxyyy, g_x_0_xyyyy_xxxxyyz, g_x_0_xyyyy_xxxxyz, g_x_0_xyyyy_xxxxyzz, g_x_0_xyyyy_xxxxzz, g_x_0_xyyyy_xxxyyy, g_x_0_xyyyy_xxxyyyy, g_x_0_xyyyy_xxxyyyz, g_x_0_xyyyy_xxxyyz, g_x_0_xyyyy_xxxyyzz, g_x_0_xyyyy_xxxyzz, g_x_0_xyyyy_xxxyzzz, g_x_0_xyyyy_xxxzzz, g_x_0_xyyyy_xxyyyy, g_x_0_xyyyy_xxyyyyy, g_x_0_xyyyy_xxyyyyz, g_x_0_xyyyy_xxyyyz, g_x_0_xyyyy_xxyyyzz, g_x_0_xyyyy_xxyyzz, g_x_0_xyyyy_xxyyzzz, g_x_0_xyyyy_xxyzzz, g_x_0_xyyyy_xxyzzzz, g_x_0_xyyyy_xxzzzz, g_x_0_xyyyy_xyyyyy, g_x_0_xyyyy_xyyyyyy, g_x_0_xyyyy_xyyyyyz, g_x_0_xyyyy_xyyyyz, g_x_0_xyyyy_xyyyyzz, g_x_0_xyyyy_xyyyzz, g_x_0_xyyyy_xyyyzzz, g_x_0_xyyyy_xyyzzz, g_x_0_xyyyy_xyyzzzz, g_x_0_xyyyy_xyzzzz, g_x_0_xyyyy_xyzzzzz, g_x_0_xyyyy_xzzzzz, g_x_0_xyyyy_yyyyyy, g_x_0_xyyyy_yyyyyyy, g_x_0_xyyyy_yyyyyyz, g_x_0_xyyyy_yyyyyz, g_x_0_xyyyy_yyyyyzz, g_x_0_xyyyy_yyyyzz, g_x_0_xyyyy_yyyyzzz, g_x_0_xyyyy_yyyzzz, g_x_0_xyyyy_yyyzzzz, g_x_0_xyyyy_yyzzzz, g_x_0_xyyyy_yyzzzzz, g_x_0_xyyyy_yzzzzz, g_x_0_xyyyy_yzzzzzz, g_x_0_xyyyy_zzzzzz, g_x_0_xyyyyy_xxxxxx, g_x_0_xyyyyy_xxxxxy, g_x_0_xyyyyy_xxxxxz, g_x_0_xyyyyy_xxxxyy, g_x_0_xyyyyy_xxxxyz, g_x_0_xyyyyy_xxxxzz, g_x_0_xyyyyy_xxxyyy, g_x_0_xyyyyy_xxxyyz, g_x_0_xyyyyy_xxxyzz, g_x_0_xyyyyy_xxxzzz, g_x_0_xyyyyy_xxyyyy, g_x_0_xyyyyy_xxyyyz, g_x_0_xyyyyy_xxyyzz, g_x_0_xyyyyy_xxyzzz, g_x_0_xyyyyy_xxzzzz, g_x_0_xyyyyy_xyyyyy, g_x_0_xyyyyy_xyyyyz, g_x_0_xyyyyy_xyyyzz, g_x_0_xyyyyy_xyyzzz, g_x_0_xyyyyy_xyzzzz, g_x_0_xyyyyy_xzzzzz, g_x_0_xyyyyy_yyyyyy, g_x_0_xyyyyy_yyyyyz, g_x_0_xyyyyy_yyyyzz, g_x_0_xyyyyy_yyyzzz, g_x_0_xyyyyy_yyzzzz, g_x_0_xyyyyy_yzzzzz, g_x_0_xyyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_xxxxxx[k] = -g_x_0_xyyyy_xxxxxx[k] * ab_y + g_x_0_xyyyy_xxxxxxy[k];

                g_x_0_xyyyyy_xxxxxy[k] = -g_x_0_xyyyy_xxxxxy[k] * ab_y + g_x_0_xyyyy_xxxxxyy[k];

                g_x_0_xyyyyy_xxxxxz[k] = -g_x_0_xyyyy_xxxxxz[k] * ab_y + g_x_0_xyyyy_xxxxxyz[k];

                g_x_0_xyyyyy_xxxxyy[k] = -g_x_0_xyyyy_xxxxyy[k] * ab_y + g_x_0_xyyyy_xxxxyyy[k];

                g_x_0_xyyyyy_xxxxyz[k] = -g_x_0_xyyyy_xxxxyz[k] * ab_y + g_x_0_xyyyy_xxxxyyz[k];

                g_x_0_xyyyyy_xxxxzz[k] = -g_x_0_xyyyy_xxxxzz[k] * ab_y + g_x_0_xyyyy_xxxxyzz[k];

                g_x_0_xyyyyy_xxxyyy[k] = -g_x_0_xyyyy_xxxyyy[k] * ab_y + g_x_0_xyyyy_xxxyyyy[k];

                g_x_0_xyyyyy_xxxyyz[k] = -g_x_0_xyyyy_xxxyyz[k] * ab_y + g_x_0_xyyyy_xxxyyyz[k];

                g_x_0_xyyyyy_xxxyzz[k] = -g_x_0_xyyyy_xxxyzz[k] * ab_y + g_x_0_xyyyy_xxxyyzz[k];

                g_x_0_xyyyyy_xxxzzz[k] = -g_x_0_xyyyy_xxxzzz[k] * ab_y + g_x_0_xyyyy_xxxyzzz[k];

                g_x_0_xyyyyy_xxyyyy[k] = -g_x_0_xyyyy_xxyyyy[k] * ab_y + g_x_0_xyyyy_xxyyyyy[k];

                g_x_0_xyyyyy_xxyyyz[k] = -g_x_0_xyyyy_xxyyyz[k] * ab_y + g_x_0_xyyyy_xxyyyyz[k];

                g_x_0_xyyyyy_xxyyzz[k] = -g_x_0_xyyyy_xxyyzz[k] * ab_y + g_x_0_xyyyy_xxyyyzz[k];

                g_x_0_xyyyyy_xxyzzz[k] = -g_x_0_xyyyy_xxyzzz[k] * ab_y + g_x_0_xyyyy_xxyyzzz[k];

                g_x_0_xyyyyy_xxzzzz[k] = -g_x_0_xyyyy_xxzzzz[k] * ab_y + g_x_0_xyyyy_xxyzzzz[k];

                g_x_0_xyyyyy_xyyyyy[k] = -g_x_0_xyyyy_xyyyyy[k] * ab_y + g_x_0_xyyyy_xyyyyyy[k];

                g_x_0_xyyyyy_xyyyyz[k] = -g_x_0_xyyyy_xyyyyz[k] * ab_y + g_x_0_xyyyy_xyyyyyz[k];

                g_x_0_xyyyyy_xyyyzz[k] = -g_x_0_xyyyy_xyyyzz[k] * ab_y + g_x_0_xyyyy_xyyyyzz[k];

                g_x_0_xyyyyy_xyyzzz[k] = -g_x_0_xyyyy_xyyzzz[k] * ab_y + g_x_0_xyyyy_xyyyzzz[k];

                g_x_0_xyyyyy_xyzzzz[k] = -g_x_0_xyyyy_xyzzzz[k] * ab_y + g_x_0_xyyyy_xyyzzzz[k];

                g_x_0_xyyyyy_xzzzzz[k] = -g_x_0_xyyyy_xzzzzz[k] * ab_y + g_x_0_xyyyy_xyzzzzz[k];

                g_x_0_xyyyyy_yyyyyy[k] = -g_x_0_xyyyy_yyyyyy[k] * ab_y + g_x_0_xyyyy_yyyyyyy[k];

                g_x_0_xyyyyy_yyyyyz[k] = -g_x_0_xyyyy_yyyyyz[k] * ab_y + g_x_0_xyyyy_yyyyyyz[k];

                g_x_0_xyyyyy_yyyyzz[k] = -g_x_0_xyyyy_yyyyzz[k] * ab_y + g_x_0_xyyyy_yyyyyzz[k];

                g_x_0_xyyyyy_yyyzzz[k] = -g_x_0_xyyyy_yyyzzz[k] * ab_y + g_x_0_xyyyy_yyyyzzz[k];

                g_x_0_xyyyyy_yyzzzz[k] = -g_x_0_xyyyy_yyzzzz[k] * ab_y + g_x_0_xyyyy_yyyzzzz[k];

                g_x_0_xyyyyy_yzzzzz[k] = -g_x_0_xyyyy_yzzzzz[k] * ab_y + g_x_0_xyyyy_yyzzzzz[k];

                g_x_0_xyyyyy_zzzzzz[k] = -g_x_0_xyyyy_zzzzzz[k] * ab_y + g_x_0_xyyyy_yzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 448 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 449 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 450 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 451 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 452 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 453 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 454 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 455 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 456 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 457 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 458 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 459 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 460 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 461 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 462 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 463 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 464 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 465 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 466 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 467 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 468 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 469 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 470 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 471 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 472 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 473 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 474 * ccomps * dcomps);

            auto g_x_0_xyyyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 475 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyyz_xxxxxx, g_x_0_xyyyyz_xxxxxy, g_x_0_xyyyyz_xxxxxz, g_x_0_xyyyyz_xxxxyy, g_x_0_xyyyyz_xxxxyz, g_x_0_xyyyyz_xxxxzz, g_x_0_xyyyyz_xxxyyy, g_x_0_xyyyyz_xxxyyz, g_x_0_xyyyyz_xxxyzz, g_x_0_xyyyyz_xxxzzz, g_x_0_xyyyyz_xxyyyy, g_x_0_xyyyyz_xxyyyz, g_x_0_xyyyyz_xxyyzz, g_x_0_xyyyyz_xxyzzz, g_x_0_xyyyyz_xxzzzz, g_x_0_xyyyyz_xyyyyy, g_x_0_xyyyyz_xyyyyz, g_x_0_xyyyyz_xyyyzz, g_x_0_xyyyyz_xyyzzz, g_x_0_xyyyyz_xyzzzz, g_x_0_xyyyyz_xzzzzz, g_x_0_xyyyyz_yyyyyy, g_x_0_xyyyyz_yyyyyz, g_x_0_xyyyyz_yyyyzz, g_x_0_xyyyyz_yyyzzz, g_x_0_xyyyyz_yyzzzz, g_x_0_xyyyyz_yzzzzz, g_x_0_xyyyyz_zzzzzz, g_x_0_xyyyz_xxxxxx, g_x_0_xyyyz_xxxxxxy, g_x_0_xyyyz_xxxxxy, g_x_0_xyyyz_xxxxxyy, g_x_0_xyyyz_xxxxxyz, g_x_0_xyyyz_xxxxxz, g_x_0_xyyyz_xxxxyy, g_x_0_xyyyz_xxxxyyy, g_x_0_xyyyz_xxxxyyz, g_x_0_xyyyz_xxxxyz, g_x_0_xyyyz_xxxxyzz, g_x_0_xyyyz_xxxxzz, g_x_0_xyyyz_xxxyyy, g_x_0_xyyyz_xxxyyyy, g_x_0_xyyyz_xxxyyyz, g_x_0_xyyyz_xxxyyz, g_x_0_xyyyz_xxxyyzz, g_x_0_xyyyz_xxxyzz, g_x_0_xyyyz_xxxyzzz, g_x_0_xyyyz_xxxzzz, g_x_0_xyyyz_xxyyyy, g_x_0_xyyyz_xxyyyyy, g_x_0_xyyyz_xxyyyyz, g_x_0_xyyyz_xxyyyz, g_x_0_xyyyz_xxyyyzz, g_x_0_xyyyz_xxyyzz, g_x_0_xyyyz_xxyyzzz, g_x_0_xyyyz_xxyzzz, g_x_0_xyyyz_xxyzzzz, g_x_0_xyyyz_xxzzzz, g_x_0_xyyyz_xyyyyy, g_x_0_xyyyz_xyyyyyy, g_x_0_xyyyz_xyyyyyz, g_x_0_xyyyz_xyyyyz, g_x_0_xyyyz_xyyyyzz, g_x_0_xyyyz_xyyyzz, g_x_0_xyyyz_xyyyzzz, g_x_0_xyyyz_xyyzzz, g_x_0_xyyyz_xyyzzzz, g_x_0_xyyyz_xyzzzz, g_x_0_xyyyz_xyzzzzz, g_x_0_xyyyz_xzzzzz, g_x_0_xyyyz_yyyyyy, g_x_0_xyyyz_yyyyyyy, g_x_0_xyyyz_yyyyyyz, g_x_0_xyyyz_yyyyyz, g_x_0_xyyyz_yyyyyzz, g_x_0_xyyyz_yyyyzz, g_x_0_xyyyz_yyyyzzz, g_x_0_xyyyz_yyyzzz, g_x_0_xyyyz_yyyzzzz, g_x_0_xyyyz_yyzzzz, g_x_0_xyyyz_yyzzzzz, g_x_0_xyyyz_yzzzzz, g_x_0_xyyyz_yzzzzzz, g_x_0_xyyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_xxxxxx[k] = -g_x_0_xyyyz_xxxxxx[k] * ab_y + g_x_0_xyyyz_xxxxxxy[k];

                g_x_0_xyyyyz_xxxxxy[k] = -g_x_0_xyyyz_xxxxxy[k] * ab_y + g_x_0_xyyyz_xxxxxyy[k];

                g_x_0_xyyyyz_xxxxxz[k] = -g_x_0_xyyyz_xxxxxz[k] * ab_y + g_x_0_xyyyz_xxxxxyz[k];

                g_x_0_xyyyyz_xxxxyy[k] = -g_x_0_xyyyz_xxxxyy[k] * ab_y + g_x_0_xyyyz_xxxxyyy[k];

                g_x_0_xyyyyz_xxxxyz[k] = -g_x_0_xyyyz_xxxxyz[k] * ab_y + g_x_0_xyyyz_xxxxyyz[k];

                g_x_0_xyyyyz_xxxxzz[k] = -g_x_0_xyyyz_xxxxzz[k] * ab_y + g_x_0_xyyyz_xxxxyzz[k];

                g_x_0_xyyyyz_xxxyyy[k] = -g_x_0_xyyyz_xxxyyy[k] * ab_y + g_x_0_xyyyz_xxxyyyy[k];

                g_x_0_xyyyyz_xxxyyz[k] = -g_x_0_xyyyz_xxxyyz[k] * ab_y + g_x_0_xyyyz_xxxyyyz[k];

                g_x_0_xyyyyz_xxxyzz[k] = -g_x_0_xyyyz_xxxyzz[k] * ab_y + g_x_0_xyyyz_xxxyyzz[k];

                g_x_0_xyyyyz_xxxzzz[k] = -g_x_0_xyyyz_xxxzzz[k] * ab_y + g_x_0_xyyyz_xxxyzzz[k];

                g_x_0_xyyyyz_xxyyyy[k] = -g_x_0_xyyyz_xxyyyy[k] * ab_y + g_x_0_xyyyz_xxyyyyy[k];

                g_x_0_xyyyyz_xxyyyz[k] = -g_x_0_xyyyz_xxyyyz[k] * ab_y + g_x_0_xyyyz_xxyyyyz[k];

                g_x_0_xyyyyz_xxyyzz[k] = -g_x_0_xyyyz_xxyyzz[k] * ab_y + g_x_0_xyyyz_xxyyyzz[k];

                g_x_0_xyyyyz_xxyzzz[k] = -g_x_0_xyyyz_xxyzzz[k] * ab_y + g_x_0_xyyyz_xxyyzzz[k];

                g_x_0_xyyyyz_xxzzzz[k] = -g_x_0_xyyyz_xxzzzz[k] * ab_y + g_x_0_xyyyz_xxyzzzz[k];

                g_x_0_xyyyyz_xyyyyy[k] = -g_x_0_xyyyz_xyyyyy[k] * ab_y + g_x_0_xyyyz_xyyyyyy[k];

                g_x_0_xyyyyz_xyyyyz[k] = -g_x_0_xyyyz_xyyyyz[k] * ab_y + g_x_0_xyyyz_xyyyyyz[k];

                g_x_0_xyyyyz_xyyyzz[k] = -g_x_0_xyyyz_xyyyzz[k] * ab_y + g_x_0_xyyyz_xyyyyzz[k];

                g_x_0_xyyyyz_xyyzzz[k] = -g_x_0_xyyyz_xyyzzz[k] * ab_y + g_x_0_xyyyz_xyyyzzz[k];

                g_x_0_xyyyyz_xyzzzz[k] = -g_x_0_xyyyz_xyzzzz[k] * ab_y + g_x_0_xyyyz_xyyzzzz[k];

                g_x_0_xyyyyz_xzzzzz[k] = -g_x_0_xyyyz_xzzzzz[k] * ab_y + g_x_0_xyyyz_xyzzzzz[k];

                g_x_0_xyyyyz_yyyyyy[k] = -g_x_0_xyyyz_yyyyyy[k] * ab_y + g_x_0_xyyyz_yyyyyyy[k];

                g_x_0_xyyyyz_yyyyyz[k] = -g_x_0_xyyyz_yyyyyz[k] * ab_y + g_x_0_xyyyz_yyyyyyz[k];

                g_x_0_xyyyyz_yyyyzz[k] = -g_x_0_xyyyz_yyyyzz[k] * ab_y + g_x_0_xyyyz_yyyyyzz[k];

                g_x_0_xyyyyz_yyyzzz[k] = -g_x_0_xyyyz_yyyzzz[k] * ab_y + g_x_0_xyyyz_yyyyzzz[k];

                g_x_0_xyyyyz_yyzzzz[k] = -g_x_0_xyyyz_yyzzzz[k] * ab_y + g_x_0_xyyyz_yyyzzzz[k];

                g_x_0_xyyyyz_yzzzzz[k] = -g_x_0_xyyyz_yzzzzz[k] * ab_y + g_x_0_xyyyz_yyzzzzz[k];

                g_x_0_xyyyyz_zzzzzz[k] = -g_x_0_xyyyz_zzzzzz[k] * ab_y + g_x_0_xyyyz_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 476 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 477 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 478 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 479 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 480 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 481 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 482 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 483 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 484 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 485 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 486 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 487 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 488 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 489 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 490 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 491 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 492 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 493 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 494 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 495 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 496 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 497 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 498 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 499 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 500 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 501 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 502 * ccomps * dcomps);

            auto g_x_0_xyyyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyzz_xxxxxx, g_x_0_xyyyzz_xxxxxy, g_x_0_xyyyzz_xxxxxz, g_x_0_xyyyzz_xxxxyy, g_x_0_xyyyzz_xxxxyz, g_x_0_xyyyzz_xxxxzz, g_x_0_xyyyzz_xxxyyy, g_x_0_xyyyzz_xxxyyz, g_x_0_xyyyzz_xxxyzz, g_x_0_xyyyzz_xxxzzz, g_x_0_xyyyzz_xxyyyy, g_x_0_xyyyzz_xxyyyz, g_x_0_xyyyzz_xxyyzz, g_x_0_xyyyzz_xxyzzz, g_x_0_xyyyzz_xxzzzz, g_x_0_xyyyzz_xyyyyy, g_x_0_xyyyzz_xyyyyz, g_x_0_xyyyzz_xyyyzz, g_x_0_xyyyzz_xyyzzz, g_x_0_xyyyzz_xyzzzz, g_x_0_xyyyzz_xzzzzz, g_x_0_xyyyzz_yyyyyy, g_x_0_xyyyzz_yyyyyz, g_x_0_xyyyzz_yyyyzz, g_x_0_xyyyzz_yyyzzz, g_x_0_xyyyzz_yyzzzz, g_x_0_xyyyzz_yzzzzz, g_x_0_xyyyzz_zzzzzz, g_x_0_xyyzz_xxxxxx, g_x_0_xyyzz_xxxxxxy, g_x_0_xyyzz_xxxxxy, g_x_0_xyyzz_xxxxxyy, g_x_0_xyyzz_xxxxxyz, g_x_0_xyyzz_xxxxxz, g_x_0_xyyzz_xxxxyy, g_x_0_xyyzz_xxxxyyy, g_x_0_xyyzz_xxxxyyz, g_x_0_xyyzz_xxxxyz, g_x_0_xyyzz_xxxxyzz, g_x_0_xyyzz_xxxxzz, g_x_0_xyyzz_xxxyyy, g_x_0_xyyzz_xxxyyyy, g_x_0_xyyzz_xxxyyyz, g_x_0_xyyzz_xxxyyz, g_x_0_xyyzz_xxxyyzz, g_x_0_xyyzz_xxxyzz, g_x_0_xyyzz_xxxyzzz, g_x_0_xyyzz_xxxzzz, g_x_0_xyyzz_xxyyyy, g_x_0_xyyzz_xxyyyyy, g_x_0_xyyzz_xxyyyyz, g_x_0_xyyzz_xxyyyz, g_x_0_xyyzz_xxyyyzz, g_x_0_xyyzz_xxyyzz, g_x_0_xyyzz_xxyyzzz, g_x_0_xyyzz_xxyzzz, g_x_0_xyyzz_xxyzzzz, g_x_0_xyyzz_xxzzzz, g_x_0_xyyzz_xyyyyy, g_x_0_xyyzz_xyyyyyy, g_x_0_xyyzz_xyyyyyz, g_x_0_xyyzz_xyyyyz, g_x_0_xyyzz_xyyyyzz, g_x_0_xyyzz_xyyyzz, g_x_0_xyyzz_xyyyzzz, g_x_0_xyyzz_xyyzzz, g_x_0_xyyzz_xyyzzzz, g_x_0_xyyzz_xyzzzz, g_x_0_xyyzz_xyzzzzz, g_x_0_xyyzz_xzzzzz, g_x_0_xyyzz_yyyyyy, g_x_0_xyyzz_yyyyyyy, g_x_0_xyyzz_yyyyyyz, g_x_0_xyyzz_yyyyyz, g_x_0_xyyzz_yyyyyzz, g_x_0_xyyzz_yyyyzz, g_x_0_xyyzz_yyyyzzz, g_x_0_xyyzz_yyyzzz, g_x_0_xyyzz_yyyzzzz, g_x_0_xyyzz_yyzzzz, g_x_0_xyyzz_yyzzzzz, g_x_0_xyyzz_yzzzzz, g_x_0_xyyzz_yzzzzzz, g_x_0_xyyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_xxxxxx[k] = -g_x_0_xyyzz_xxxxxx[k] * ab_y + g_x_0_xyyzz_xxxxxxy[k];

                g_x_0_xyyyzz_xxxxxy[k] = -g_x_0_xyyzz_xxxxxy[k] * ab_y + g_x_0_xyyzz_xxxxxyy[k];

                g_x_0_xyyyzz_xxxxxz[k] = -g_x_0_xyyzz_xxxxxz[k] * ab_y + g_x_0_xyyzz_xxxxxyz[k];

                g_x_0_xyyyzz_xxxxyy[k] = -g_x_0_xyyzz_xxxxyy[k] * ab_y + g_x_0_xyyzz_xxxxyyy[k];

                g_x_0_xyyyzz_xxxxyz[k] = -g_x_0_xyyzz_xxxxyz[k] * ab_y + g_x_0_xyyzz_xxxxyyz[k];

                g_x_0_xyyyzz_xxxxzz[k] = -g_x_0_xyyzz_xxxxzz[k] * ab_y + g_x_0_xyyzz_xxxxyzz[k];

                g_x_0_xyyyzz_xxxyyy[k] = -g_x_0_xyyzz_xxxyyy[k] * ab_y + g_x_0_xyyzz_xxxyyyy[k];

                g_x_0_xyyyzz_xxxyyz[k] = -g_x_0_xyyzz_xxxyyz[k] * ab_y + g_x_0_xyyzz_xxxyyyz[k];

                g_x_0_xyyyzz_xxxyzz[k] = -g_x_0_xyyzz_xxxyzz[k] * ab_y + g_x_0_xyyzz_xxxyyzz[k];

                g_x_0_xyyyzz_xxxzzz[k] = -g_x_0_xyyzz_xxxzzz[k] * ab_y + g_x_0_xyyzz_xxxyzzz[k];

                g_x_0_xyyyzz_xxyyyy[k] = -g_x_0_xyyzz_xxyyyy[k] * ab_y + g_x_0_xyyzz_xxyyyyy[k];

                g_x_0_xyyyzz_xxyyyz[k] = -g_x_0_xyyzz_xxyyyz[k] * ab_y + g_x_0_xyyzz_xxyyyyz[k];

                g_x_0_xyyyzz_xxyyzz[k] = -g_x_0_xyyzz_xxyyzz[k] * ab_y + g_x_0_xyyzz_xxyyyzz[k];

                g_x_0_xyyyzz_xxyzzz[k] = -g_x_0_xyyzz_xxyzzz[k] * ab_y + g_x_0_xyyzz_xxyyzzz[k];

                g_x_0_xyyyzz_xxzzzz[k] = -g_x_0_xyyzz_xxzzzz[k] * ab_y + g_x_0_xyyzz_xxyzzzz[k];

                g_x_0_xyyyzz_xyyyyy[k] = -g_x_0_xyyzz_xyyyyy[k] * ab_y + g_x_0_xyyzz_xyyyyyy[k];

                g_x_0_xyyyzz_xyyyyz[k] = -g_x_0_xyyzz_xyyyyz[k] * ab_y + g_x_0_xyyzz_xyyyyyz[k];

                g_x_0_xyyyzz_xyyyzz[k] = -g_x_0_xyyzz_xyyyzz[k] * ab_y + g_x_0_xyyzz_xyyyyzz[k];

                g_x_0_xyyyzz_xyyzzz[k] = -g_x_0_xyyzz_xyyzzz[k] * ab_y + g_x_0_xyyzz_xyyyzzz[k];

                g_x_0_xyyyzz_xyzzzz[k] = -g_x_0_xyyzz_xyzzzz[k] * ab_y + g_x_0_xyyzz_xyyzzzz[k];

                g_x_0_xyyyzz_xzzzzz[k] = -g_x_0_xyyzz_xzzzzz[k] * ab_y + g_x_0_xyyzz_xyzzzzz[k];

                g_x_0_xyyyzz_yyyyyy[k] = -g_x_0_xyyzz_yyyyyy[k] * ab_y + g_x_0_xyyzz_yyyyyyy[k];

                g_x_0_xyyyzz_yyyyyz[k] = -g_x_0_xyyzz_yyyyyz[k] * ab_y + g_x_0_xyyzz_yyyyyyz[k];

                g_x_0_xyyyzz_yyyyzz[k] = -g_x_0_xyyzz_yyyyzz[k] * ab_y + g_x_0_xyyzz_yyyyyzz[k];

                g_x_0_xyyyzz_yyyzzz[k] = -g_x_0_xyyzz_yyyzzz[k] * ab_y + g_x_0_xyyzz_yyyyzzz[k];

                g_x_0_xyyyzz_yyzzzz[k] = -g_x_0_xyyzz_yyzzzz[k] * ab_y + g_x_0_xyyzz_yyyzzzz[k];

                g_x_0_xyyyzz_yzzzzz[k] = -g_x_0_xyyzz_yzzzzz[k] * ab_y + g_x_0_xyyzz_yyzzzzz[k];

                g_x_0_xyyyzz_zzzzzz[k] = -g_x_0_xyyzz_zzzzzz[k] * ab_y + g_x_0_xyyzz_yzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 504 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 505 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 506 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 507 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 508 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 509 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 510 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 511 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 512 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 513 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 514 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 515 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 516 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 517 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 518 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 519 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 520 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 521 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 522 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 523 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 524 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 525 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 526 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 527 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 528 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 529 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 530 * ccomps * dcomps);

            auto g_x_0_xyyzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 531 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyzzz_xxxxxx, g_x_0_xyyzzz_xxxxxy, g_x_0_xyyzzz_xxxxxz, g_x_0_xyyzzz_xxxxyy, g_x_0_xyyzzz_xxxxyz, g_x_0_xyyzzz_xxxxzz, g_x_0_xyyzzz_xxxyyy, g_x_0_xyyzzz_xxxyyz, g_x_0_xyyzzz_xxxyzz, g_x_0_xyyzzz_xxxzzz, g_x_0_xyyzzz_xxyyyy, g_x_0_xyyzzz_xxyyyz, g_x_0_xyyzzz_xxyyzz, g_x_0_xyyzzz_xxyzzz, g_x_0_xyyzzz_xxzzzz, g_x_0_xyyzzz_xyyyyy, g_x_0_xyyzzz_xyyyyz, g_x_0_xyyzzz_xyyyzz, g_x_0_xyyzzz_xyyzzz, g_x_0_xyyzzz_xyzzzz, g_x_0_xyyzzz_xzzzzz, g_x_0_xyyzzz_yyyyyy, g_x_0_xyyzzz_yyyyyz, g_x_0_xyyzzz_yyyyzz, g_x_0_xyyzzz_yyyzzz, g_x_0_xyyzzz_yyzzzz, g_x_0_xyyzzz_yzzzzz, g_x_0_xyyzzz_zzzzzz, g_x_0_xyzzz_xxxxxx, g_x_0_xyzzz_xxxxxxy, g_x_0_xyzzz_xxxxxy, g_x_0_xyzzz_xxxxxyy, g_x_0_xyzzz_xxxxxyz, g_x_0_xyzzz_xxxxxz, g_x_0_xyzzz_xxxxyy, g_x_0_xyzzz_xxxxyyy, g_x_0_xyzzz_xxxxyyz, g_x_0_xyzzz_xxxxyz, g_x_0_xyzzz_xxxxyzz, g_x_0_xyzzz_xxxxzz, g_x_0_xyzzz_xxxyyy, g_x_0_xyzzz_xxxyyyy, g_x_0_xyzzz_xxxyyyz, g_x_0_xyzzz_xxxyyz, g_x_0_xyzzz_xxxyyzz, g_x_0_xyzzz_xxxyzz, g_x_0_xyzzz_xxxyzzz, g_x_0_xyzzz_xxxzzz, g_x_0_xyzzz_xxyyyy, g_x_0_xyzzz_xxyyyyy, g_x_0_xyzzz_xxyyyyz, g_x_0_xyzzz_xxyyyz, g_x_0_xyzzz_xxyyyzz, g_x_0_xyzzz_xxyyzz, g_x_0_xyzzz_xxyyzzz, g_x_0_xyzzz_xxyzzz, g_x_0_xyzzz_xxyzzzz, g_x_0_xyzzz_xxzzzz, g_x_0_xyzzz_xyyyyy, g_x_0_xyzzz_xyyyyyy, g_x_0_xyzzz_xyyyyyz, g_x_0_xyzzz_xyyyyz, g_x_0_xyzzz_xyyyyzz, g_x_0_xyzzz_xyyyzz, g_x_0_xyzzz_xyyyzzz, g_x_0_xyzzz_xyyzzz, g_x_0_xyzzz_xyyzzzz, g_x_0_xyzzz_xyzzzz, g_x_0_xyzzz_xyzzzzz, g_x_0_xyzzz_xzzzzz, g_x_0_xyzzz_yyyyyy, g_x_0_xyzzz_yyyyyyy, g_x_0_xyzzz_yyyyyyz, g_x_0_xyzzz_yyyyyz, g_x_0_xyzzz_yyyyyzz, g_x_0_xyzzz_yyyyzz, g_x_0_xyzzz_yyyyzzz, g_x_0_xyzzz_yyyzzz, g_x_0_xyzzz_yyyzzzz, g_x_0_xyzzz_yyzzzz, g_x_0_xyzzz_yyzzzzz, g_x_0_xyzzz_yzzzzz, g_x_0_xyzzz_yzzzzzz, g_x_0_xyzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_xxxxxx[k] = -g_x_0_xyzzz_xxxxxx[k] * ab_y + g_x_0_xyzzz_xxxxxxy[k];

                g_x_0_xyyzzz_xxxxxy[k] = -g_x_0_xyzzz_xxxxxy[k] * ab_y + g_x_0_xyzzz_xxxxxyy[k];

                g_x_0_xyyzzz_xxxxxz[k] = -g_x_0_xyzzz_xxxxxz[k] * ab_y + g_x_0_xyzzz_xxxxxyz[k];

                g_x_0_xyyzzz_xxxxyy[k] = -g_x_0_xyzzz_xxxxyy[k] * ab_y + g_x_0_xyzzz_xxxxyyy[k];

                g_x_0_xyyzzz_xxxxyz[k] = -g_x_0_xyzzz_xxxxyz[k] * ab_y + g_x_0_xyzzz_xxxxyyz[k];

                g_x_0_xyyzzz_xxxxzz[k] = -g_x_0_xyzzz_xxxxzz[k] * ab_y + g_x_0_xyzzz_xxxxyzz[k];

                g_x_0_xyyzzz_xxxyyy[k] = -g_x_0_xyzzz_xxxyyy[k] * ab_y + g_x_0_xyzzz_xxxyyyy[k];

                g_x_0_xyyzzz_xxxyyz[k] = -g_x_0_xyzzz_xxxyyz[k] * ab_y + g_x_0_xyzzz_xxxyyyz[k];

                g_x_0_xyyzzz_xxxyzz[k] = -g_x_0_xyzzz_xxxyzz[k] * ab_y + g_x_0_xyzzz_xxxyyzz[k];

                g_x_0_xyyzzz_xxxzzz[k] = -g_x_0_xyzzz_xxxzzz[k] * ab_y + g_x_0_xyzzz_xxxyzzz[k];

                g_x_0_xyyzzz_xxyyyy[k] = -g_x_0_xyzzz_xxyyyy[k] * ab_y + g_x_0_xyzzz_xxyyyyy[k];

                g_x_0_xyyzzz_xxyyyz[k] = -g_x_0_xyzzz_xxyyyz[k] * ab_y + g_x_0_xyzzz_xxyyyyz[k];

                g_x_0_xyyzzz_xxyyzz[k] = -g_x_0_xyzzz_xxyyzz[k] * ab_y + g_x_0_xyzzz_xxyyyzz[k];

                g_x_0_xyyzzz_xxyzzz[k] = -g_x_0_xyzzz_xxyzzz[k] * ab_y + g_x_0_xyzzz_xxyyzzz[k];

                g_x_0_xyyzzz_xxzzzz[k] = -g_x_0_xyzzz_xxzzzz[k] * ab_y + g_x_0_xyzzz_xxyzzzz[k];

                g_x_0_xyyzzz_xyyyyy[k] = -g_x_0_xyzzz_xyyyyy[k] * ab_y + g_x_0_xyzzz_xyyyyyy[k];

                g_x_0_xyyzzz_xyyyyz[k] = -g_x_0_xyzzz_xyyyyz[k] * ab_y + g_x_0_xyzzz_xyyyyyz[k];

                g_x_0_xyyzzz_xyyyzz[k] = -g_x_0_xyzzz_xyyyzz[k] * ab_y + g_x_0_xyzzz_xyyyyzz[k];

                g_x_0_xyyzzz_xyyzzz[k] = -g_x_0_xyzzz_xyyzzz[k] * ab_y + g_x_0_xyzzz_xyyyzzz[k];

                g_x_0_xyyzzz_xyzzzz[k] = -g_x_0_xyzzz_xyzzzz[k] * ab_y + g_x_0_xyzzz_xyyzzzz[k];

                g_x_0_xyyzzz_xzzzzz[k] = -g_x_0_xyzzz_xzzzzz[k] * ab_y + g_x_0_xyzzz_xyzzzzz[k];

                g_x_0_xyyzzz_yyyyyy[k] = -g_x_0_xyzzz_yyyyyy[k] * ab_y + g_x_0_xyzzz_yyyyyyy[k];

                g_x_0_xyyzzz_yyyyyz[k] = -g_x_0_xyzzz_yyyyyz[k] * ab_y + g_x_0_xyzzz_yyyyyyz[k];

                g_x_0_xyyzzz_yyyyzz[k] = -g_x_0_xyzzz_yyyyzz[k] * ab_y + g_x_0_xyzzz_yyyyyzz[k];

                g_x_0_xyyzzz_yyyzzz[k] = -g_x_0_xyzzz_yyyzzz[k] * ab_y + g_x_0_xyzzz_yyyyzzz[k];

                g_x_0_xyyzzz_yyzzzz[k] = -g_x_0_xyzzz_yyzzzz[k] * ab_y + g_x_0_xyzzz_yyyzzzz[k];

                g_x_0_xyyzzz_yzzzzz[k] = -g_x_0_xyzzz_yzzzzz[k] * ab_y + g_x_0_xyzzz_yyzzzzz[k];

                g_x_0_xyyzzz_zzzzzz[k] = -g_x_0_xyzzz_zzzzzz[k] * ab_y + g_x_0_xyzzz_yzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 532 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 533 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 534 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 535 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 536 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 537 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 538 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 539 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 540 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 541 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 542 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 543 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 544 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 545 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 546 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 547 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 548 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 549 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 550 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 551 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 552 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 553 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 554 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 555 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 556 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 557 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 558 * ccomps * dcomps);

            auto g_x_0_xyzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzzzz_xxxxxx, g_x_0_xyzzzz_xxxxxy, g_x_0_xyzzzz_xxxxxz, g_x_0_xyzzzz_xxxxyy, g_x_0_xyzzzz_xxxxyz, g_x_0_xyzzzz_xxxxzz, g_x_0_xyzzzz_xxxyyy, g_x_0_xyzzzz_xxxyyz, g_x_0_xyzzzz_xxxyzz, g_x_0_xyzzzz_xxxzzz, g_x_0_xyzzzz_xxyyyy, g_x_0_xyzzzz_xxyyyz, g_x_0_xyzzzz_xxyyzz, g_x_0_xyzzzz_xxyzzz, g_x_0_xyzzzz_xxzzzz, g_x_0_xyzzzz_xyyyyy, g_x_0_xyzzzz_xyyyyz, g_x_0_xyzzzz_xyyyzz, g_x_0_xyzzzz_xyyzzz, g_x_0_xyzzzz_xyzzzz, g_x_0_xyzzzz_xzzzzz, g_x_0_xyzzzz_yyyyyy, g_x_0_xyzzzz_yyyyyz, g_x_0_xyzzzz_yyyyzz, g_x_0_xyzzzz_yyyzzz, g_x_0_xyzzzz_yyzzzz, g_x_0_xyzzzz_yzzzzz, g_x_0_xyzzzz_zzzzzz, g_x_0_xzzzz_xxxxxx, g_x_0_xzzzz_xxxxxxy, g_x_0_xzzzz_xxxxxy, g_x_0_xzzzz_xxxxxyy, g_x_0_xzzzz_xxxxxyz, g_x_0_xzzzz_xxxxxz, g_x_0_xzzzz_xxxxyy, g_x_0_xzzzz_xxxxyyy, g_x_0_xzzzz_xxxxyyz, g_x_0_xzzzz_xxxxyz, g_x_0_xzzzz_xxxxyzz, g_x_0_xzzzz_xxxxzz, g_x_0_xzzzz_xxxyyy, g_x_0_xzzzz_xxxyyyy, g_x_0_xzzzz_xxxyyyz, g_x_0_xzzzz_xxxyyz, g_x_0_xzzzz_xxxyyzz, g_x_0_xzzzz_xxxyzz, g_x_0_xzzzz_xxxyzzz, g_x_0_xzzzz_xxxzzz, g_x_0_xzzzz_xxyyyy, g_x_0_xzzzz_xxyyyyy, g_x_0_xzzzz_xxyyyyz, g_x_0_xzzzz_xxyyyz, g_x_0_xzzzz_xxyyyzz, g_x_0_xzzzz_xxyyzz, g_x_0_xzzzz_xxyyzzz, g_x_0_xzzzz_xxyzzz, g_x_0_xzzzz_xxyzzzz, g_x_0_xzzzz_xxzzzz, g_x_0_xzzzz_xyyyyy, g_x_0_xzzzz_xyyyyyy, g_x_0_xzzzz_xyyyyyz, g_x_0_xzzzz_xyyyyz, g_x_0_xzzzz_xyyyyzz, g_x_0_xzzzz_xyyyzz, g_x_0_xzzzz_xyyyzzz, g_x_0_xzzzz_xyyzzz, g_x_0_xzzzz_xyyzzzz, g_x_0_xzzzz_xyzzzz, g_x_0_xzzzz_xyzzzzz, g_x_0_xzzzz_xzzzzz, g_x_0_xzzzz_yyyyyy, g_x_0_xzzzz_yyyyyyy, g_x_0_xzzzz_yyyyyyz, g_x_0_xzzzz_yyyyyz, g_x_0_xzzzz_yyyyyzz, g_x_0_xzzzz_yyyyzz, g_x_0_xzzzz_yyyyzzz, g_x_0_xzzzz_yyyzzz, g_x_0_xzzzz_yyyzzzz, g_x_0_xzzzz_yyzzzz, g_x_0_xzzzz_yyzzzzz, g_x_0_xzzzz_yzzzzz, g_x_0_xzzzz_yzzzzzz, g_x_0_xzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_xxxxxx[k] = -g_x_0_xzzzz_xxxxxx[k] * ab_y + g_x_0_xzzzz_xxxxxxy[k];

                g_x_0_xyzzzz_xxxxxy[k] = -g_x_0_xzzzz_xxxxxy[k] * ab_y + g_x_0_xzzzz_xxxxxyy[k];

                g_x_0_xyzzzz_xxxxxz[k] = -g_x_0_xzzzz_xxxxxz[k] * ab_y + g_x_0_xzzzz_xxxxxyz[k];

                g_x_0_xyzzzz_xxxxyy[k] = -g_x_0_xzzzz_xxxxyy[k] * ab_y + g_x_0_xzzzz_xxxxyyy[k];

                g_x_0_xyzzzz_xxxxyz[k] = -g_x_0_xzzzz_xxxxyz[k] * ab_y + g_x_0_xzzzz_xxxxyyz[k];

                g_x_0_xyzzzz_xxxxzz[k] = -g_x_0_xzzzz_xxxxzz[k] * ab_y + g_x_0_xzzzz_xxxxyzz[k];

                g_x_0_xyzzzz_xxxyyy[k] = -g_x_0_xzzzz_xxxyyy[k] * ab_y + g_x_0_xzzzz_xxxyyyy[k];

                g_x_0_xyzzzz_xxxyyz[k] = -g_x_0_xzzzz_xxxyyz[k] * ab_y + g_x_0_xzzzz_xxxyyyz[k];

                g_x_0_xyzzzz_xxxyzz[k] = -g_x_0_xzzzz_xxxyzz[k] * ab_y + g_x_0_xzzzz_xxxyyzz[k];

                g_x_0_xyzzzz_xxxzzz[k] = -g_x_0_xzzzz_xxxzzz[k] * ab_y + g_x_0_xzzzz_xxxyzzz[k];

                g_x_0_xyzzzz_xxyyyy[k] = -g_x_0_xzzzz_xxyyyy[k] * ab_y + g_x_0_xzzzz_xxyyyyy[k];

                g_x_0_xyzzzz_xxyyyz[k] = -g_x_0_xzzzz_xxyyyz[k] * ab_y + g_x_0_xzzzz_xxyyyyz[k];

                g_x_0_xyzzzz_xxyyzz[k] = -g_x_0_xzzzz_xxyyzz[k] * ab_y + g_x_0_xzzzz_xxyyyzz[k];

                g_x_0_xyzzzz_xxyzzz[k] = -g_x_0_xzzzz_xxyzzz[k] * ab_y + g_x_0_xzzzz_xxyyzzz[k];

                g_x_0_xyzzzz_xxzzzz[k] = -g_x_0_xzzzz_xxzzzz[k] * ab_y + g_x_0_xzzzz_xxyzzzz[k];

                g_x_0_xyzzzz_xyyyyy[k] = -g_x_0_xzzzz_xyyyyy[k] * ab_y + g_x_0_xzzzz_xyyyyyy[k];

                g_x_0_xyzzzz_xyyyyz[k] = -g_x_0_xzzzz_xyyyyz[k] * ab_y + g_x_0_xzzzz_xyyyyyz[k];

                g_x_0_xyzzzz_xyyyzz[k] = -g_x_0_xzzzz_xyyyzz[k] * ab_y + g_x_0_xzzzz_xyyyyzz[k];

                g_x_0_xyzzzz_xyyzzz[k] = -g_x_0_xzzzz_xyyzzz[k] * ab_y + g_x_0_xzzzz_xyyyzzz[k];

                g_x_0_xyzzzz_xyzzzz[k] = -g_x_0_xzzzz_xyzzzz[k] * ab_y + g_x_0_xzzzz_xyyzzzz[k];

                g_x_0_xyzzzz_xzzzzz[k] = -g_x_0_xzzzz_xzzzzz[k] * ab_y + g_x_0_xzzzz_xyzzzzz[k];

                g_x_0_xyzzzz_yyyyyy[k] = -g_x_0_xzzzz_yyyyyy[k] * ab_y + g_x_0_xzzzz_yyyyyyy[k];

                g_x_0_xyzzzz_yyyyyz[k] = -g_x_0_xzzzz_yyyyyz[k] * ab_y + g_x_0_xzzzz_yyyyyyz[k];

                g_x_0_xyzzzz_yyyyzz[k] = -g_x_0_xzzzz_yyyyzz[k] * ab_y + g_x_0_xzzzz_yyyyyzz[k];

                g_x_0_xyzzzz_yyyzzz[k] = -g_x_0_xzzzz_yyyzzz[k] * ab_y + g_x_0_xzzzz_yyyyzzz[k];

                g_x_0_xyzzzz_yyzzzz[k] = -g_x_0_xzzzz_yyzzzz[k] * ab_y + g_x_0_xzzzz_yyyzzzz[k];

                g_x_0_xyzzzz_yzzzzz[k] = -g_x_0_xzzzz_yzzzzz[k] * ab_y + g_x_0_xzzzz_yyzzzzz[k];

                g_x_0_xyzzzz_zzzzzz[k] = -g_x_0_xzzzz_zzzzzz[k] * ab_y + g_x_0_xzzzz_yzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 560 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 561 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 562 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 563 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 564 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 565 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 566 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 567 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 568 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 569 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 570 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 571 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 572 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 573 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 574 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 575 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 576 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 577 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 578 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 579 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 580 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 581 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 582 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 583 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 584 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 585 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 586 * ccomps * dcomps);

            auto g_x_0_xzzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzzzz_xxxxxx, g_x_0_xzzzz_xxxxxxz, g_x_0_xzzzz_xxxxxy, g_x_0_xzzzz_xxxxxyz, g_x_0_xzzzz_xxxxxz, g_x_0_xzzzz_xxxxxzz, g_x_0_xzzzz_xxxxyy, g_x_0_xzzzz_xxxxyyz, g_x_0_xzzzz_xxxxyz, g_x_0_xzzzz_xxxxyzz, g_x_0_xzzzz_xxxxzz, g_x_0_xzzzz_xxxxzzz, g_x_0_xzzzz_xxxyyy, g_x_0_xzzzz_xxxyyyz, g_x_0_xzzzz_xxxyyz, g_x_0_xzzzz_xxxyyzz, g_x_0_xzzzz_xxxyzz, g_x_0_xzzzz_xxxyzzz, g_x_0_xzzzz_xxxzzz, g_x_0_xzzzz_xxxzzzz, g_x_0_xzzzz_xxyyyy, g_x_0_xzzzz_xxyyyyz, g_x_0_xzzzz_xxyyyz, g_x_0_xzzzz_xxyyyzz, g_x_0_xzzzz_xxyyzz, g_x_0_xzzzz_xxyyzzz, g_x_0_xzzzz_xxyzzz, g_x_0_xzzzz_xxyzzzz, g_x_0_xzzzz_xxzzzz, g_x_0_xzzzz_xxzzzzz, g_x_0_xzzzz_xyyyyy, g_x_0_xzzzz_xyyyyyz, g_x_0_xzzzz_xyyyyz, g_x_0_xzzzz_xyyyyzz, g_x_0_xzzzz_xyyyzz, g_x_0_xzzzz_xyyyzzz, g_x_0_xzzzz_xyyzzz, g_x_0_xzzzz_xyyzzzz, g_x_0_xzzzz_xyzzzz, g_x_0_xzzzz_xyzzzzz, g_x_0_xzzzz_xzzzzz, g_x_0_xzzzz_xzzzzzz, g_x_0_xzzzz_yyyyyy, g_x_0_xzzzz_yyyyyyz, g_x_0_xzzzz_yyyyyz, g_x_0_xzzzz_yyyyyzz, g_x_0_xzzzz_yyyyzz, g_x_0_xzzzz_yyyyzzz, g_x_0_xzzzz_yyyzzz, g_x_0_xzzzz_yyyzzzz, g_x_0_xzzzz_yyzzzz, g_x_0_xzzzz_yyzzzzz, g_x_0_xzzzz_yzzzzz, g_x_0_xzzzz_yzzzzzz, g_x_0_xzzzz_zzzzzz, g_x_0_xzzzz_zzzzzzz, g_x_0_xzzzzz_xxxxxx, g_x_0_xzzzzz_xxxxxy, g_x_0_xzzzzz_xxxxxz, g_x_0_xzzzzz_xxxxyy, g_x_0_xzzzzz_xxxxyz, g_x_0_xzzzzz_xxxxzz, g_x_0_xzzzzz_xxxyyy, g_x_0_xzzzzz_xxxyyz, g_x_0_xzzzzz_xxxyzz, g_x_0_xzzzzz_xxxzzz, g_x_0_xzzzzz_xxyyyy, g_x_0_xzzzzz_xxyyyz, g_x_0_xzzzzz_xxyyzz, g_x_0_xzzzzz_xxyzzz, g_x_0_xzzzzz_xxzzzz, g_x_0_xzzzzz_xyyyyy, g_x_0_xzzzzz_xyyyyz, g_x_0_xzzzzz_xyyyzz, g_x_0_xzzzzz_xyyzzz, g_x_0_xzzzzz_xyzzzz, g_x_0_xzzzzz_xzzzzz, g_x_0_xzzzzz_yyyyyy, g_x_0_xzzzzz_yyyyyz, g_x_0_xzzzzz_yyyyzz, g_x_0_xzzzzz_yyyzzz, g_x_0_xzzzzz_yyzzzz, g_x_0_xzzzzz_yzzzzz, g_x_0_xzzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_xxxxxx[k] = -g_x_0_xzzzz_xxxxxx[k] * ab_z + g_x_0_xzzzz_xxxxxxz[k];

                g_x_0_xzzzzz_xxxxxy[k] = -g_x_0_xzzzz_xxxxxy[k] * ab_z + g_x_0_xzzzz_xxxxxyz[k];

                g_x_0_xzzzzz_xxxxxz[k] = -g_x_0_xzzzz_xxxxxz[k] * ab_z + g_x_0_xzzzz_xxxxxzz[k];

                g_x_0_xzzzzz_xxxxyy[k] = -g_x_0_xzzzz_xxxxyy[k] * ab_z + g_x_0_xzzzz_xxxxyyz[k];

                g_x_0_xzzzzz_xxxxyz[k] = -g_x_0_xzzzz_xxxxyz[k] * ab_z + g_x_0_xzzzz_xxxxyzz[k];

                g_x_0_xzzzzz_xxxxzz[k] = -g_x_0_xzzzz_xxxxzz[k] * ab_z + g_x_0_xzzzz_xxxxzzz[k];

                g_x_0_xzzzzz_xxxyyy[k] = -g_x_0_xzzzz_xxxyyy[k] * ab_z + g_x_0_xzzzz_xxxyyyz[k];

                g_x_0_xzzzzz_xxxyyz[k] = -g_x_0_xzzzz_xxxyyz[k] * ab_z + g_x_0_xzzzz_xxxyyzz[k];

                g_x_0_xzzzzz_xxxyzz[k] = -g_x_0_xzzzz_xxxyzz[k] * ab_z + g_x_0_xzzzz_xxxyzzz[k];

                g_x_0_xzzzzz_xxxzzz[k] = -g_x_0_xzzzz_xxxzzz[k] * ab_z + g_x_0_xzzzz_xxxzzzz[k];

                g_x_0_xzzzzz_xxyyyy[k] = -g_x_0_xzzzz_xxyyyy[k] * ab_z + g_x_0_xzzzz_xxyyyyz[k];

                g_x_0_xzzzzz_xxyyyz[k] = -g_x_0_xzzzz_xxyyyz[k] * ab_z + g_x_0_xzzzz_xxyyyzz[k];

                g_x_0_xzzzzz_xxyyzz[k] = -g_x_0_xzzzz_xxyyzz[k] * ab_z + g_x_0_xzzzz_xxyyzzz[k];

                g_x_0_xzzzzz_xxyzzz[k] = -g_x_0_xzzzz_xxyzzz[k] * ab_z + g_x_0_xzzzz_xxyzzzz[k];

                g_x_0_xzzzzz_xxzzzz[k] = -g_x_0_xzzzz_xxzzzz[k] * ab_z + g_x_0_xzzzz_xxzzzzz[k];

                g_x_0_xzzzzz_xyyyyy[k] = -g_x_0_xzzzz_xyyyyy[k] * ab_z + g_x_0_xzzzz_xyyyyyz[k];

                g_x_0_xzzzzz_xyyyyz[k] = -g_x_0_xzzzz_xyyyyz[k] * ab_z + g_x_0_xzzzz_xyyyyzz[k];

                g_x_0_xzzzzz_xyyyzz[k] = -g_x_0_xzzzz_xyyyzz[k] * ab_z + g_x_0_xzzzz_xyyyzzz[k];

                g_x_0_xzzzzz_xyyzzz[k] = -g_x_0_xzzzz_xyyzzz[k] * ab_z + g_x_0_xzzzz_xyyzzzz[k];

                g_x_0_xzzzzz_xyzzzz[k] = -g_x_0_xzzzz_xyzzzz[k] * ab_z + g_x_0_xzzzz_xyzzzzz[k];

                g_x_0_xzzzzz_xzzzzz[k] = -g_x_0_xzzzz_xzzzzz[k] * ab_z + g_x_0_xzzzz_xzzzzzz[k];

                g_x_0_xzzzzz_yyyyyy[k] = -g_x_0_xzzzz_yyyyyy[k] * ab_z + g_x_0_xzzzz_yyyyyyz[k];

                g_x_0_xzzzzz_yyyyyz[k] = -g_x_0_xzzzz_yyyyyz[k] * ab_z + g_x_0_xzzzz_yyyyyzz[k];

                g_x_0_xzzzzz_yyyyzz[k] = -g_x_0_xzzzz_yyyyzz[k] * ab_z + g_x_0_xzzzz_yyyyzzz[k];

                g_x_0_xzzzzz_yyyzzz[k] = -g_x_0_xzzzz_yyyzzz[k] * ab_z + g_x_0_xzzzz_yyyzzzz[k];

                g_x_0_xzzzzz_yyzzzz[k] = -g_x_0_xzzzz_yyzzzz[k] * ab_z + g_x_0_xzzzz_yyzzzzz[k];

                g_x_0_xzzzzz_yzzzzz[k] = -g_x_0_xzzzz_yzzzzz[k] * ab_z + g_x_0_xzzzz_yzzzzzz[k];

                g_x_0_xzzzzz_zzzzzz[k] = -g_x_0_xzzzz_zzzzzz[k] * ab_z + g_x_0_xzzzz_zzzzzzz[k];
            }

            /// Set up 588-616 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 588 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 589 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 590 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 591 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 592 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 593 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 594 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 595 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 596 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 597 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 598 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 599 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 600 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 601 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 602 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 603 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 604 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 605 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 606 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 607 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 608 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 609 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 610 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 611 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 612 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 613 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 614 * ccomps * dcomps);

            auto g_x_0_yyyyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 615 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyy_xxxxxx, g_x_0_yyyyy_xxxxxxy, g_x_0_yyyyy_xxxxxy, g_x_0_yyyyy_xxxxxyy, g_x_0_yyyyy_xxxxxyz, g_x_0_yyyyy_xxxxxz, g_x_0_yyyyy_xxxxyy, g_x_0_yyyyy_xxxxyyy, g_x_0_yyyyy_xxxxyyz, g_x_0_yyyyy_xxxxyz, g_x_0_yyyyy_xxxxyzz, g_x_0_yyyyy_xxxxzz, g_x_0_yyyyy_xxxyyy, g_x_0_yyyyy_xxxyyyy, g_x_0_yyyyy_xxxyyyz, g_x_0_yyyyy_xxxyyz, g_x_0_yyyyy_xxxyyzz, g_x_0_yyyyy_xxxyzz, g_x_0_yyyyy_xxxyzzz, g_x_0_yyyyy_xxxzzz, g_x_0_yyyyy_xxyyyy, g_x_0_yyyyy_xxyyyyy, g_x_0_yyyyy_xxyyyyz, g_x_0_yyyyy_xxyyyz, g_x_0_yyyyy_xxyyyzz, g_x_0_yyyyy_xxyyzz, g_x_0_yyyyy_xxyyzzz, g_x_0_yyyyy_xxyzzz, g_x_0_yyyyy_xxyzzzz, g_x_0_yyyyy_xxzzzz, g_x_0_yyyyy_xyyyyy, g_x_0_yyyyy_xyyyyyy, g_x_0_yyyyy_xyyyyyz, g_x_0_yyyyy_xyyyyz, g_x_0_yyyyy_xyyyyzz, g_x_0_yyyyy_xyyyzz, g_x_0_yyyyy_xyyyzzz, g_x_0_yyyyy_xyyzzz, g_x_0_yyyyy_xyyzzzz, g_x_0_yyyyy_xyzzzz, g_x_0_yyyyy_xyzzzzz, g_x_0_yyyyy_xzzzzz, g_x_0_yyyyy_yyyyyy, g_x_0_yyyyy_yyyyyyy, g_x_0_yyyyy_yyyyyyz, g_x_0_yyyyy_yyyyyz, g_x_0_yyyyy_yyyyyzz, g_x_0_yyyyy_yyyyzz, g_x_0_yyyyy_yyyyzzz, g_x_0_yyyyy_yyyzzz, g_x_0_yyyyy_yyyzzzz, g_x_0_yyyyy_yyzzzz, g_x_0_yyyyy_yyzzzzz, g_x_0_yyyyy_yzzzzz, g_x_0_yyyyy_yzzzzzz, g_x_0_yyyyy_zzzzzz, g_x_0_yyyyyy_xxxxxx, g_x_0_yyyyyy_xxxxxy, g_x_0_yyyyyy_xxxxxz, g_x_0_yyyyyy_xxxxyy, g_x_0_yyyyyy_xxxxyz, g_x_0_yyyyyy_xxxxzz, g_x_0_yyyyyy_xxxyyy, g_x_0_yyyyyy_xxxyyz, g_x_0_yyyyyy_xxxyzz, g_x_0_yyyyyy_xxxzzz, g_x_0_yyyyyy_xxyyyy, g_x_0_yyyyyy_xxyyyz, g_x_0_yyyyyy_xxyyzz, g_x_0_yyyyyy_xxyzzz, g_x_0_yyyyyy_xxzzzz, g_x_0_yyyyyy_xyyyyy, g_x_0_yyyyyy_xyyyyz, g_x_0_yyyyyy_xyyyzz, g_x_0_yyyyyy_xyyzzz, g_x_0_yyyyyy_xyzzzz, g_x_0_yyyyyy_xzzzzz, g_x_0_yyyyyy_yyyyyy, g_x_0_yyyyyy_yyyyyz, g_x_0_yyyyyy_yyyyzz, g_x_0_yyyyyy_yyyzzz, g_x_0_yyyyyy_yyzzzz, g_x_0_yyyyyy_yzzzzz, g_x_0_yyyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_xxxxxx[k] = -g_x_0_yyyyy_xxxxxx[k] * ab_y + g_x_0_yyyyy_xxxxxxy[k];

                g_x_0_yyyyyy_xxxxxy[k] = -g_x_0_yyyyy_xxxxxy[k] * ab_y + g_x_0_yyyyy_xxxxxyy[k];

                g_x_0_yyyyyy_xxxxxz[k] = -g_x_0_yyyyy_xxxxxz[k] * ab_y + g_x_0_yyyyy_xxxxxyz[k];

                g_x_0_yyyyyy_xxxxyy[k] = -g_x_0_yyyyy_xxxxyy[k] * ab_y + g_x_0_yyyyy_xxxxyyy[k];

                g_x_0_yyyyyy_xxxxyz[k] = -g_x_0_yyyyy_xxxxyz[k] * ab_y + g_x_0_yyyyy_xxxxyyz[k];

                g_x_0_yyyyyy_xxxxzz[k] = -g_x_0_yyyyy_xxxxzz[k] * ab_y + g_x_0_yyyyy_xxxxyzz[k];

                g_x_0_yyyyyy_xxxyyy[k] = -g_x_0_yyyyy_xxxyyy[k] * ab_y + g_x_0_yyyyy_xxxyyyy[k];

                g_x_0_yyyyyy_xxxyyz[k] = -g_x_0_yyyyy_xxxyyz[k] * ab_y + g_x_0_yyyyy_xxxyyyz[k];

                g_x_0_yyyyyy_xxxyzz[k] = -g_x_0_yyyyy_xxxyzz[k] * ab_y + g_x_0_yyyyy_xxxyyzz[k];

                g_x_0_yyyyyy_xxxzzz[k] = -g_x_0_yyyyy_xxxzzz[k] * ab_y + g_x_0_yyyyy_xxxyzzz[k];

                g_x_0_yyyyyy_xxyyyy[k] = -g_x_0_yyyyy_xxyyyy[k] * ab_y + g_x_0_yyyyy_xxyyyyy[k];

                g_x_0_yyyyyy_xxyyyz[k] = -g_x_0_yyyyy_xxyyyz[k] * ab_y + g_x_0_yyyyy_xxyyyyz[k];

                g_x_0_yyyyyy_xxyyzz[k] = -g_x_0_yyyyy_xxyyzz[k] * ab_y + g_x_0_yyyyy_xxyyyzz[k];

                g_x_0_yyyyyy_xxyzzz[k] = -g_x_0_yyyyy_xxyzzz[k] * ab_y + g_x_0_yyyyy_xxyyzzz[k];

                g_x_0_yyyyyy_xxzzzz[k] = -g_x_0_yyyyy_xxzzzz[k] * ab_y + g_x_0_yyyyy_xxyzzzz[k];

                g_x_0_yyyyyy_xyyyyy[k] = -g_x_0_yyyyy_xyyyyy[k] * ab_y + g_x_0_yyyyy_xyyyyyy[k];

                g_x_0_yyyyyy_xyyyyz[k] = -g_x_0_yyyyy_xyyyyz[k] * ab_y + g_x_0_yyyyy_xyyyyyz[k];

                g_x_0_yyyyyy_xyyyzz[k] = -g_x_0_yyyyy_xyyyzz[k] * ab_y + g_x_0_yyyyy_xyyyyzz[k];

                g_x_0_yyyyyy_xyyzzz[k] = -g_x_0_yyyyy_xyyzzz[k] * ab_y + g_x_0_yyyyy_xyyyzzz[k];

                g_x_0_yyyyyy_xyzzzz[k] = -g_x_0_yyyyy_xyzzzz[k] * ab_y + g_x_0_yyyyy_xyyzzzz[k];

                g_x_0_yyyyyy_xzzzzz[k] = -g_x_0_yyyyy_xzzzzz[k] * ab_y + g_x_0_yyyyy_xyzzzzz[k];

                g_x_0_yyyyyy_yyyyyy[k] = -g_x_0_yyyyy_yyyyyy[k] * ab_y + g_x_0_yyyyy_yyyyyyy[k];

                g_x_0_yyyyyy_yyyyyz[k] = -g_x_0_yyyyy_yyyyyz[k] * ab_y + g_x_0_yyyyy_yyyyyyz[k];

                g_x_0_yyyyyy_yyyyzz[k] = -g_x_0_yyyyy_yyyyzz[k] * ab_y + g_x_0_yyyyy_yyyyyzz[k];

                g_x_0_yyyyyy_yyyzzz[k] = -g_x_0_yyyyy_yyyzzz[k] * ab_y + g_x_0_yyyyy_yyyyzzz[k];

                g_x_0_yyyyyy_yyzzzz[k] = -g_x_0_yyyyy_yyzzzz[k] * ab_y + g_x_0_yyyyy_yyyzzzz[k];

                g_x_0_yyyyyy_yzzzzz[k] = -g_x_0_yyyyy_yzzzzz[k] * ab_y + g_x_0_yyyyy_yyzzzzz[k];

                g_x_0_yyyyyy_zzzzzz[k] = -g_x_0_yyyyy_zzzzzz[k] * ab_y + g_x_0_yyyyy_yzzzzzz[k];
            }

            /// Set up 616-644 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 616 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 617 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 618 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 619 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 620 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 621 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 622 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 623 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 624 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 625 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 626 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 627 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 628 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 629 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 630 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 631 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 632 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 633 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 634 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 635 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 636 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 637 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 638 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 639 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 640 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 641 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 642 * ccomps * dcomps);

            auto g_x_0_yyyyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 643 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyyz_xxxxxx, g_x_0_yyyyyz_xxxxxy, g_x_0_yyyyyz_xxxxxz, g_x_0_yyyyyz_xxxxyy, g_x_0_yyyyyz_xxxxyz, g_x_0_yyyyyz_xxxxzz, g_x_0_yyyyyz_xxxyyy, g_x_0_yyyyyz_xxxyyz, g_x_0_yyyyyz_xxxyzz, g_x_0_yyyyyz_xxxzzz, g_x_0_yyyyyz_xxyyyy, g_x_0_yyyyyz_xxyyyz, g_x_0_yyyyyz_xxyyzz, g_x_0_yyyyyz_xxyzzz, g_x_0_yyyyyz_xxzzzz, g_x_0_yyyyyz_xyyyyy, g_x_0_yyyyyz_xyyyyz, g_x_0_yyyyyz_xyyyzz, g_x_0_yyyyyz_xyyzzz, g_x_0_yyyyyz_xyzzzz, g_x_0_yyyyyz_xzzzzz, g_x_0_yyyyyz_yyyyyy, g_x_0_yyyyyz_yyyyyz, g_x_0_yyyyyz_yyyyzz, g_x_0_yyyyyz_yyyzzz, g_x_0_yyyyyz_yyzzzz, g_x_0_yyyyyz_yzzzzz, g_x_0_yyyyyz_zzzzzz, g_x_0_yyyyz_xxxxxx, g_x_0_yyyyz_xxxxxxy, g_x_0_yyyyz_xxxxxy, g_x_0_yyyyz_xxxxxyy, g_x_0_yyyyz_xxxxxyz, g_x_0_yyyyz_xxxxxz, g_x_0_yyyyz_xxxxyy, g_x_0_yyyyz_xxxxyyy, g_x_0_yyyyz_xxxxyyz, g_x_0_yyyyz_xxxxyz, g_x_0_yyyyz_xxxxyzz, g_x_0_yyyyz_xxxxzz, g_x_0_yyyyz_xxxyyy, g_x_0_yyyyz_xxxyyyy, g_x_0_yyyyz_xxxyyyz, g_x_0_yyyyz_xxxyyz, g_x_0_yyyyz_xxxyyzz, g_x_0_yyyyz_xxxyzz, g_x_0_yyyyz_xxxyzzz, g_x_0_yyyyz_xxxzzz, g_x_0_yyyyz_xxyyyy, g_x_0_yyyyz_xxyyyyy, g_x_0_yyyyz_xxyyyyz, g_x_0_yyyyz_xxyyyz, g_x_0_yyyyz_xxyyyzz, g_x_0_yyyyz_xxyyzz, g_x_0_yyyyz_xxyyzzz, g_x_0_yyyyz_xxyzzz, g_x_0_yyyyz_xxyzzzz, g_x_0_yyyyz_xxzzzz, g_x_0_yyyyz_xyyyyy, g_x_0_yyyyz_xyyyyyy, g_x_0_yyyyz_xyyyyyz, g_x_0_yyyyz_xyyyyz, g_x_0_yyyyz_xyyyyzz, g_x_0_yyyyz_xyyyzz, g_x_0_yyyyz_xyyyzzz, g_x_0_yyyyz_xyyzzz, g_x_0_yyyyz_xyyzzzz, g_x_0_yyyyz_xyzzzz, g_x_0_yyyyz_xyzzzzz, g_x_0_yyyyz_xzzzzz, g_x_0_yyyyz_yyyyyy, g_x_0_yyyyz_yyyyyyy, g_x_0_yyyyz_yyyyyyz, g_x_0_yyyyz_yyyyyz, g_x_0_yyyyz_yyyyyzz, g_x_0_yyyyz_yyyyzz, g_x_0_yyyyz_yyyyzzz, g_x_0_yyyyz_yyyzzz, g_x_0_yyyyz_yyyzzzz, g_x_0_yyyyz_yyzzzz, g_x_0_yyyyz_yyzzzzz, g_x_0_yyyyz_yzzzzz, g_x_0_yyyyz_yzzzzzz, g_x_0_yyyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_xxxxxx[k] = -g_x_0_yyyyz_xxxxxx[k] * ab_y + g_x_0_yyyyz_xxxxxxy[k];

                g_x_0_yyyyyz_xxxxxy[k] = -g_x_0_yyyyz_xxxxxy[k] * ab_y + g_x_0_yyyyz_xxxxxyy[k];

                g_x_0_yyyyyz_xxxxxz[k] = -g_x_0_yyyyz_xxxxxz[k] * ab_y + g_x_0_yyyyz_xxxxxyz[k];

                g_x_0_yyyyyz_xxxxyy[k] = -g_x_0_yyyyz_xxxxyy[k] * ab_y + g_x_0_yyyyz_xxxxyyy[k];

                g_x_0_yyyyyz_xxxxyz[k] = -g_x_0_yyyyz_xxxxyz[k] * ab_y + g_x_0_yyyyz_xxxxyyz[k];

                g_x_0_yyyyyz_xxxxzz[k] = -g_x_0_yyyyz_xxxxzz[k] * ab_y + g_x_0_yyyyz_xxxxyzz[k];

                g_x_0_yyyyyz_xxxyyy[k] = -g_x_0_yyyyz_xxxyyy[k] * ab_y + g_x_0_yyyyz_xxxyyyy[k];

                g_x_0_yyyyyz_xxxyyz[k] = -g_x_0_yyyyz_xxxyyz[k] * ab_y + g_x_0_yyyyz_xxxyyyz[k];

                g_x_0_yyyyyz_xxxyzz[k] = -g_x_0_yyyyz_xxxyzz[k] * ab_y + g_x_0_yyyyz_xxxyyzz[k];

                g_x_0_yyyyyz_xxxzzz[k] = -g_x_0_yyyyz_xxxzzz[k] * ab_y + g_x_0_yyyyz_xxxyzzz[k];

                g_x_0_yyyyyz_xxyyyy[k] = -g_x_0_yyyyz_xxyyyy[k] * ab_y + g_x_0_yyyyz_xxyyyyy[k];

                g_x_0_yyyyyz_xxyyyz[k] = -g_x_0_yyyyz_xxyyyz[k] * ab_y + g_x_0_yyyyz_xxyyyyz[k];

                g_x_0_yyyyyz_xxyyzz[k] = -g_x_0_yyyyz_xxyyzz[k] * ab_y + g_x_0_yyyyz_xxyyyzz[k];

                g_x_0_yyyyyz_xxyzzz[k] = -g_x_0_yyyyz_xxyzzz[k] * ab_y + g_x_0_yyyyz_xxyyzzz[k];

                g_x_0_yyyyyz_xxzzzz[k] = -g_x_0_yyyyz_xxzzzz[k] * ab_y + g_x_0_yyyyz_xxyzzzz[k];

                g_x_0_yyyyyz_xyyyyy[k] = -g_x_0_yyyyz_xyyyyy[k] * ab_y + g_x_0_yyyyz_xyyyyyy[k];

                g_x_0_yyyyyz_xyyyyz[k] = -g_x_0_yyyyz_xyyyyz[k] * ab_y + g_x_0_yyyyz_xyyyyyz[k];

                g_x_0_yyyyyz_xyyyzz[k] = -g_x_0_yyyyz_xyyyzz[k] * ab_y + g_x_0_yyyyz_xyyyyzz[k];

                g_x_0_yyyyyz_xyyzzz[k] = -g_x_0_yyyyz_xyyzzz[k] * ab_y + g_x_0_yyyyz_xyyyzzz[k];

                g_x_0_yyyyyz_xyzzzz[k] = -g_x_0_yyyyz_xyzzzz[k] * ab_y + g_x_0_yyyyz_xyyzzzz[k];

                g_x_0_yyyyyz_xzzzzz[k] = -g_x_0_yyyyz_xzzzzz[k] * ab_y + g_x_0_yyyyz_xyzzzzz[k];

                g_x_0_yyyyyz_yyyyyy[k] = -g_x_0_yyyyz_yyyyyy[k] * ab_y + g_x_0_yyyyz_yyyyyyy[k];

                g_x_0_yyyyyz_yyyyyz[k] = -g_x_0_yyyyz_yyyyyz[k] * ab_y + g_x_0_yyyyz_yyyyyyz[k];

                g_x_0_yyyyyz_yyyyzz[k] = -g_x_0_yyyyz_yyyyzz[k] * ab_y + g_x_0_yyyyz_yyyyyzz[k];

                g_x_0_yyyyyz_yyyzzz[k] = -g_x_0_yyyyz_yyyzzz[k] * ab_y + g_x_0_yyyyz_yyyyzzz[k];

                g_x_0_yyyyyz_yyzzzz[k] = -g_x_0_yyyyz_yyzzzz[k] * ab_y + g_x_0_yyyyz_yyyzzzz[k];

                g_x_0_yyyyyz_yzzzzz[k] = -g_x_0_yyyyz_yzzzzz[k] * ab_y + g_x_0_yyyyz_yyzzzzz[k];

                g_x_0_yyyyyz_zzzzzz[k] = -g_x_0_yyyyz_zzzzzz[k] * ab_y + g_x_0_yyyyz_yzzzzzz[k];
            }

            /// Set up 644-672 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 644 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 645 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 646 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 647 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 648 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 649 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 650 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 651 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 652 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 653 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 654 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 655 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 656 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 657 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 658 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 659 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 660 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 661 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 662 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 663 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 664 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 665 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 666 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 667 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 668 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 669 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 670 * ccomps * dcomps);

            auto g_x_0_yyyyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyzz_xxxxxx, g_x_0_yyyyzz_xxxxxy, g_x_0_yyyyzz_xxxxxz, g_x_0_yyyyzz_xxxxyy, g_x_0_yyyyzz_xxxxyz, g_x_0_yyyyzz_xxxxzz, g_x_0_yyyyzz_xxxyyy, g_x_0_yyyyzz_xxxyyz, g_x_0_yyyyzz_xxxyzz, g_x_0_yyyyzz_xxxzzz, g_x_0_yyyyzz_xxyyyy, g_x_0_yyyyzz_xxyyyz, g_x_0_yyyyzz_xxyyzz, g_x_0_yyyyzz_xxyzzz, g_x_0_yyyyzz_xxzzzz, g_x_0_yyyyzz_xyyyyy, g_x_0_yyyyzz_xyyyyz, g_x_0_yyyyzz_xyyyzz, g_x_0_yyyyzz_xyyzzz, g_x_0_yyyyzz_xyzzzz, g_x_0_yyyyzz_xzzzzz, g_x_0_yyyyzz_yyyyyy, g_x_0_yyyyzz_yyyyyz, g_x_0_yyyyzz_yyyyzz, g_x_0_yyyyzz_yyyzzz, g_x_0_yyyyzz_yyzzzz, g_x_0_yyyyzz_yzzzzz, g_x_0_yyyyzz_zzzzzz, g_x_0_yyyzz_xxxxxx, g_x_0_yyyzz_xxxxxxy, g_x_0_yyyzz_xxxxxy, g_x_0_yyyzz_xxxxxyy, g_x_0_yyyzz_xxxxxyz, g_x_0_yyyzz_xxxxxz, g_x_0_yyyzz_xxxxyy, g_x_0_yyyzz_xxxxyyy, g_x_0_yyyzz_xxxxyyz, g_x_0_yyyzz_xxxxyz, g_x_0_yyyzz_xxxxyzz, g_x_0_yyyzz_xxxxzz, g_x_0_yyyzz_xxxyyy, g_x_0_yyyzz_xxxyyyy, g_x_0_yyyzz_xxxyyyz, g_x_0_yyyzz_xxxyyz, g_x_0_yyyzz_xxxyyzz, g_x_0_yyyzz_xxxyzz, g_x_0_yyyzz_xxxyzzz, g_x_0_yyyzz_xxxzzz, g_x_0_yyyzz_xxyyyy, g_x_0_yyyzz_xxyyyyy, g_x_0_yyyzz_xxyyyyz, g_x_0_yyyzz_xxyyyz, g_x_0_yyyzz_xxyyyzz, g_x_0_yyyzz_xxyyzz, g_x_0_yyyzz_xxyyzzz, g_x_0_yyyzz_xxyzzz, g_x_0_yyyzz_xxyzzzz, g_x_0_yyyzz_xxzzzz, g_x_0_yyyzz_xyyyyy, g_x_0_yyyzz_xyyyyyy, g_x_0_yyyzz_xyyyyyz, g_x_0_yyyzz_xyyyyz, g_x_0_yyyzz_xyyyyzz, g_x_0_yyyzz_xyyyzz, g_x_0_yyyzz_xyyyzzz, g_x_0_yyyzz_xyyzzz, g_x_0_yyyzz_xyyzzzz, g_x_0_yyyzz_xyzzzz, g_x_0_yyyzz_xyzzzzz, g_x_0_yyyzz_xzzzzz, g_x_0_yyyzz_yyyyyy, g_x_0_yyyzz_yyyyyyy, g_x_0_yyyzz_yyyyyyz, g_x_0_yyyzz_yyyyyz, g_x_0_yyyzz_yyyyyzz, g_x_0_yyyzz_yyyyzz, g_x_0_yyyzz_yyyyzzz, g_x_0_yyyzz_yyyzzz, g_x_0_yyyzz_yyyzzzz, g_x_0_yyyzz_yyzzzz, g_x_0_yyyzz_yyzzzzz, g_x_0_yyyzz_yzzzzz, g_x_0_yyyzz_yzzzzzz, g_x_0_yyyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_xxxxxx[k] = -g_x_0_yyyzz_xxxxxx[k] * ab_y + g_x_0_yyyzz_xxxxxxy[k];

                g_x_0_yyyyzz_xxxxxy[k] = -g_x_0_yyyzz_xxxxxy[k] * ab_y + g_x_0_yyyzz_xxxxxyy[k];

                g_x_0_yyyyzz_xxxxxz[k] = -g_x_0_yyyzz_xxxxxz[k] * ab_y + g_x_0_yyyzz_xxxxxyz[k];

                g_x_0_yyyyzz_xxxxyy[k] = -g_x_0_yyyzz_xxxxyy[k] * ab_y + g_x_0_yyyzz_xxxxyyy[k];

                g_x_0_yyyyzz_xxxxyz[k] = -g_x_0_yyyzz_xxxxyz[k] * ab_y + g_x_0_yyyzz_xxxxyyz[k];

                g_x_0_yyyyzz_xxxxzz[k] = -g_x_0_yyyzz_xxxxzz[k] * ab_y + g_x_0_yyyzz_xxxxyzz[k];

                g_x_0_yyyyzz_xxxyyy[k] = -g_x_0_yyyzz_xxxyyy[k] * ab_y + g_x_0_yyyzz_xxxyyyy[k];

                g_x_0_yyyyzz_xxxyyz[k] = -g_x_0_yyyzz_xxxyyz[k] * ab_y + g_x_0_yyyzz_xxxyyyz[k];

                g_x_0_yyyyzz_xxxyzz[k] = -g_x_0_yyyzz_xxxyzz[k] * ab_y + g_x_0_yyyzz_xxxyyzz[k];

                g_x_0_yyyyzz_xxxzzz[k] = -g_x_0_yyyzz_xxxzzz[k] * ab_y + g_x_0_yyyzz_xxxyzzz[k];

                g_x_0_yyyyzz_xxyyyy[k] = -g_x_0_yyyzz_xxyyyy[k] * ab_y + g_x_0_yyyzz_xxyyyyy[k];

                g_x_0_yyyyzz_xxyyyz[k] = -g_x_0_yyyzz_xxyyyz[k] * ab_y + g_x_0_yyyzz_xxyyyyz[k];

                g_x_0_yyyyzz_xxyyzz[k] = -g_x_0_yyyzz_xxyyzz[k] * ab_y + g_x_0_yyyzz_xxyyyzz[k];

                g_x_0_yyyyzz_xxyzzz[k] = -g_x_0_yyyzz_xxyzzz[k] * ab_y + g_x_0_yyyzz_xxyyzzz[k];

                g_x_0_yyyyzz_xxzzzz[k] = -g_x_0_yyyzz_xxzzzz[k] * ab_y + g_x_0_yyyzz_xxyzzzz[k];

                g_x_0_yyyyzz_xyyyyy[k] = -g_x_0_yyyzz_xyyyyy[k] * ab_y + g_x_0_yyyzz_xyyyyyy[k];

                g_x_0_yyyyzz_xyyyyz[k] = -g_x_0_yyyzz_xyyyyz[k] * ab_y + g_x_0_yyyzz_xyyyyyz[k];

                g_x_0_yyyyzz_xyyyzz[k] = -g_x_0_yyyzz_xyyyzz[k] * ab_y + g_x_0_yyyzz_xyyyyzz[k];

                g_x_0_yyyyzz_xyyzzz[k] = -g_x_0_yyyzz_xyyzzz[k] * ab_y + g_x_0_yyyzz_xyyyzzz[k];

                g_x_0_yyyyzz_xyzzzz[k] = -g_x_0_yyyzz_xyzzzz[k] * ab_y + g_x_0_yyyzz_xyyzzzz[k];

                g_x_0_yyyyzz_xzzzzz[k] = -g_x_0_yyyzz_xzzzzz[k] * ab_y + g_x_0_yyyzz_xyzzzzz[k];

                g_x_0_yyyyzz_yyyyyy[k] = -g_x_0_yyyzz_yyyyyy[k] * ab_y + g_x_0_yyyzz_yyyyyyy[k];

                g_x_0_yyyyzz_yyyyyz[k] = -g_x_0_yyyzz_yyyyyz[k] * ab_y + g_x_0_yyyzz_yyyyyyz[k];

                g_x_0_yyyyzz_yyyyzz[k] = -g_x_0_yyyzz_yyyyzz[k] * ab_y + g_x_0_yyyzz_yyyyyzz[k];

                g_x_0_yyyyzz_yyyzzz[k] = -g_x_0_yyyzz_yyyzzz[k] * ab_y + g_x_0_yyyzz_yyyyzzz[k];

                g_x_0_yyyyzz_yyzzzz[k] = -g_x_0_yyyzz_yyzzzz[k] * ab_y + g_x_0_yyyzz_yyyzzzz[k];

                g_x_0_yyyyzz_yzzzzz[k] = -g_x_0_yyyzz_yzzzzz[k] * ab_y + g_x_0_yyyzz_yyzzzzz[k];

                g_x_0_yyyyzz_zzzzzz[k] = -g_x_0_yyyzz_zzzzzz[k] * ab_y + g_x_0_yyyzz_yzzzzzz[k];
            }

            /// Set up 672-700 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 672 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 673 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 674 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 675 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 676 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 677 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 678 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 679 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 680 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 681 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 682 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 683 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 684 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 685 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 686 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 687 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 688 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 689 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 690 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 691 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 692 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 693 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 694 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 695 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 696 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 697 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 698 * ccomps * dcomps);

            auto g_x_0_yyyzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyzzz_xxxxxx, g_x_0_yyyzzz_xxxxxy, g_x_0_yyyzzz_xxxxxz, g_x_0_yyyzzz_xxxxyy, g_x_0_yyyzzz_xxxxyz, g_x_0_yyyzzz_xxxxzz, g_x_0_yyyzzz_xxxyyy, g_x_0_yyyzzz_xxxyyz, g_x_0_yyyzzz_xxxyzz, g_x_0_yyyzzz_xxxzzz, g_x_0_yyyzzz_xxyyyy, g_x_0_yyyzzz_xxyyyz, g_x_0_yyyzzz_xxyyzz, g_x_0_yyyzzz_xxyzzz, g_x_0_yyyzzz_xxzzzz, g_x_0_yyyzzz_xyyyyy, g_x_0_yyyzzz_xyyyyz, g_x_0_yyyzzz_xyyyzz, g_x_0_yyyzzz_xyyzzz, g_x_0_yyyzzz_xyzzzz, g_x_0_yyyzzz_xzzzzz, g_x_0_yyyzzz_yyyyyy, g_x_0_yyyzzz_yyyyyz, g_x_0_yyyzzz_yyyyzz, g_x_0_yyyzzz_yyyzzz, g_x_0_yyyzzz_yyzzzz, g_x_0_yyyzzz_yzzzzz, g_x_0_yyyzzz_zzzzzz, g_x_0_yyzzz_xxxxxx, g_x_0_yyzzz_xxxxxxy, g_x_0_yyzzz_xxxxxy, g_x_0_yyzzz_xxxxxyy, g_x_0_yyzzz_xxxxxyz, g_x_0_yyzzz_xxxxxz, g_x_0_yyzzz_xxxxyy, g_x_0_yyzzz_xxxxyyy, g_x_0_yyzzz_xxxxyyz, g_x_0_yyzzz_xxxxyz, g_x_0_yyzzz_xxxxyzz, g_x_0_yyzzz_xxxxzz, g_x_0_yyzzz_xxxyyy, g_x_0_yyzzz_xxxyyyy, g_x_0_yyzzz_xxxyyyz, g_x_0_yyzzz_xxxyyz, g_x_0_yyzzz_xxxyyzz, g_x_0_yyzzz_xxxyzz, g_x_0_yyzzz_xxxyzzz, g_x_0_yyzzz_xxxzzz, g_x_0_yyzzz_xxyyyy, g_x_0_yyzzz_xxyyyyy, g_x_0_yyzzz_xxyyyyz, g_x_0_yyzzz_xxyyyz, g_x_0_yyzzz_xxyyyzz, g_x_0_yyzzz_xxyyzz, g_x_0_yyzzz_xxyyzzz, g_x_0_yyzzz_xxyzzz, g_x_0_yyzzz_xxyzzzz, g_x_0_yyzzz_xxzzzz, g_x_0_yyzzz_xyyyyy, g_x_0_yyzzz_xyyyyyy, g_x_0_yyzzz_xyyyyyz, g_x_0_yyzzz_xyyyyz, g_x_0_yyzzz_xyyyyzz, g_x_0_yyzzz_xyyyzz, g_x_0_yyzzz_xyyyzzz, g_x_0_yyzzz_xyyzzz, g_x_0_yyzzz_xyyzzzz, g_x_0_yyzzz_xyzzzz, g_x_0_yyzzz_xyzzzzz, g_x_0_yyzzz_xzzzzz, g_x_0_yyzzz_yyyyyy, g_x_0_yyzzz_yyyyyyy, g_x_0_yyzzz_yyyyyyz, g_x_0_yyzzz_yyyyyz, g_x_0_yyzzz_yyyyyzz, g_x_0_yyzzz_yyyyzz, g_x_0_yyzzz_yyyyzzz, g_x_0_yyzzz_yyyzzz, g_x_0_yyzzz_yyyzzzz, g_x_0_yyzzz_yyzzzz, g_x_0_yyzzz_yyzzzzz, g_x_0_yyzzz_yzzzzz, g_x_0_yyzzz_yzzzzzz, g_x_0_yyzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_xxxxxx[k] = -g_x_0_yyzzz_xxxxxx[k] * ab_y + g_x_0_yyzzz_xxxxxxy[k];

                g_x_0_yyyzzz_xxxxxy[k] = -g_x_0_yyzzz_xxxxxy[k] * ab_y + g_x_0_yyzzz_xxxxxyy[k];

                g_x_0_yyyzzz_xxxxxz[k] = -g_x_0_yyzzz_xxxxxz[k] * ab_y + g_x_0_yyzzz_xxxxxyz[k];

                g_x_0_yyyzzz_xxxxyy[k] = -g_x_0_yyzzz_xxxxyy[k] * ab_y + g_x_0_yyzzz_xxxxyyy[k];

                g_x_0_yyyzzz_xxxxyz[k] = -g_x_0_yyzzz_xxxxyz[k] * ab_y + g_x_0_yyzzz_xxxxyyz[k];

                g_x_0_yyyzzz_xxxxzz[k] = -g_x_0_yyzzz_xxxxzz[k] * ab_y + g_x_0_yyzzz_xxxxyzz[k];

                g_x_0_yyyzzz_xxxyyy[k] = -g_x_0_yyzzz_xxxyyy[k] * ab_y + g_x_0_yyzzz_xxxyyyy[k];

                g_x_0_yyyzzz_xxxyyz[k] = -g_x_0_yyzzz_xxxyyz[k] * ab_y + g_x_0_yyzzz_xxxyyyz[k];

                g_x_0_yyyzzz_xxxyzz[k] = -g_x_0_yyzzz_xxxyzz[k] * ab_y + g_x_0_yyzzz_xxxyyzz[k];

                g_x_0_yyyzzz_xxxzzz[k] = -g_x_0_yyzzz_xxxzzz[k] * ab_y + g_x_0_yyzzz_xxxyzzz[k];

                g_x_0_yyyzzz_xxyyyy[k] = -g_x_0_yyzzz_xxyyyy[k] * ab_y + g_x_0_yyzzz_xxyyyyy[k];

                g_x_0_yyyzzz_xxyyyz[k] = -g_x_0_yyzzz_xxyyyz[k] * ab_y + g_x_0_yyzzz_xxyyyyz[k];

                g_x_0_yyyzzz_xxyyzz[k] = -g_x_0_yyzzz_xxyyzz[k] * ab_y + g_x_0_yyzzz_xxyyyzz[k];

                g_x_0_yyyzzz_xxyzzz[k] = -g_x_0_yyzzz_xxyzzz[k] * ab_y + g_x_0_yyzzz_xxyyzzz[k];

                g_x_0_yyyzzz_xxzzzz[k] = -g_x_0_yyzzz_xxzzzz[k] * ab_y + g_x_0_yyzzz_xxyzzzz[k];

                g_x_0_yyyzzz_xyyyyy[k] = -g_x_0_yyzzz_xyyyyy[k] * ab_y + g_x_0_yyzzz_xyyyyyy[k];

                g_x_0_yyyzzz_xyyyyz[k] = -g_x_0_yyzzz_xyyyyz[k] * ab_y + g_x_0_yyzzz_xyyyyyz[k];

                g_x_0_yyyzzz_xyyyzz[k] = -g_x_0_yyzzz_xyyyzz[k] * ab_y + g_x_0_yyzzz_xyyyyzz[k];

                g_x_0_yyyzzz_xyyzzz[k] = -g_x_0_yyzzz_xyyzzz[k] * ab_y + g_x_0_yyzzz_xyyyzzz[k];

                g_x_0_yyyzzz_xyzzzz[k] = -g_x_0_yyzzz_xyzzzz[k] * ab_y + g_x_0_yyzzz_xyyzzzz[k];

                g_x_0_yyyzzz_xzzzzz[k] = -g_x_0_yyzzz_xzzzzz[k] * ab_y + g_x_0_yyzzz_xyzzzzz[k];

                g_x_0_yyyzzz_yyyyyy[k] = -g_x_0_yyzzz_yyyyyy[k] * ab_y + g_x_0_yyzzz_yyyyyyy[k];

                g_x_0_yyyzzz_yyyyyz[k] = -g_x_0_yyzzz_yyyyyz[k] * ab_y + g_x_0_yyzzz_yyyyyyz[k];

                g_x_0_yyyzzz_yyyyzz[k] = -g_x_0_yyzzz_yyyyzz[k] * ab_y + g_x_0_yyzzz_yyyyyzz[k];

                g_x_0_yyyzzz_yyyzzz[k] = -g_x_0_yyzzz_yyyzzz[k] * ab_y + g_x_0_yyzzz_yyyyzzz[k];

                g_x_0_yyyzzz_yyzzzz[k] = -g_x_0_yyzzz_yyzzzz[k] * ab_y + g_x_0_yyzzz_yyyzzzz[k];

                g_x_0_yyyzzz_yzzzzz[k] = -g_x_0_yyzzz_yzzzzz[k] * ab_y + g_x_0_yyzzz_yyzzzzz[k];

                g_x_0_yyyzzz_zzzzzz[k] = -g_x_0_yyzzz_zzzzzz[k] * ab_y + g_x_0_yyzzz_yzzzzzz[k];
            }

            /// Set up 700-728 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 700 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 701 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 702 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 703 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 704 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 705 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 706 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 707 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 708 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 709 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 710 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 711 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 712 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 713 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 714 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 715 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 716 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 717 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 718 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 719 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 720 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 721 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 722 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 723 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 724 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 725 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 726 * ccomps * dcomps);

            auto g_x_0_yyzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 727 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzzzz_xxxxxx, g_x_0_yyzzzz_xxxxxy, g_x_0_yyzzzz_xxxxxz, g_x_0_yyzzzz_xxxxyy, g_x_0_yyzzzz_xxxxyz, g_x_0_yyzzzz_xxxxzz, g_x_0_yyzzzz_xxxyyy, g_x_0_yyzzzz_xxxyyz, g_x_0_yyzzzz_xxxyzz, g_x_0_yyzzzz_xxxzzz, g_x_0_yyzzzz_xxyyyy, g_x_0_yyzzzz_xxyyyz, g_x_0_yyzzzz_xxyyzz, g_x_0_yyzzzz_xxyzzz, g_x_0_yyzzzz_xxzzzz, g_x_0_yyzzzz_xyyyyy, g_x_0_yyzzzz_xyyyyz, g_x_0_yyzzzz_xyyyzz, g_x_0_yyzzzz_xyyzzz, g_x_0_yyzzzz_xyzzzz, g_x_0_yyzzzz_xzzzzz, g_x_0_yyzzzz_yyyyyy, g_x_0_yyzzzz_yyyyyz, g_x_0_yyzzzz_yyyyzz, g_x_0_yyzzzz_yyyzzz, g_x_0_yyzzzz_yyzzzz, g_x_0_yyzzzz_yzzzzz, g_x_0_yyzzzz_zzzzzz, g_x_0_yzzzz_xxxxxx, g_x_0_yzzzz_xxxxxxy, g_x_0_yzzzz_xxxxxy, g_x_0_yzzzz_xxxxxyy, g_x_0_yzzzz_xxxxxyz, g_x_0_yzzzz_xxxxxz, g_x_0_yzzzz_xxxxyy, g_x_0_yzzzz_xxxxyyy, g_x_0_yzzzz_xxxxyyz, g_x_0_yzzzz_xxxxyz, g_x_0_yzzzz_xxxxyzz, g_x_0_yzzzz_xxxxzz, g_x_0_yzzzz_xxxyyy, g_x_0_yzzzz_xxxyyyy, g_x_0_yzzzz_xxxyyyz, g_x_0_yzzzz_xxxyyz, g_x_0_yzzzz_xxxyyzz, g_x_0_yzzzz_xxxyzz, g_x_0_yzzzz_xxxyzzz, g_x_0_yzzzz_xxxzzz, g_x_0_yzzzz_xxyyyy, g_x_0_yzzzz_xxyyyyy, g_x_0_yzzzz_xxyyyyz, g_x_0_yzzzz_xxyyyz, g_x_0_yzzzz_xxyyyzz, g_x_0_yzzzz_xxyyzz, g_x_0_yzzzz_xxyyzzz, g_x_0_yzzzz_xxyzzz, g_x_0_yzzzz_xxyzzzz, g_x_0_yzzzz_xxzzzz, g_x_0_yzzzz_xyyyyy, g_x_0_yzzzz_xyyyyyy, g_x_0_yzzzz_xyyyyyz, g_x_0_yzzzz_xyyyyz, g_x_0_yzzzz_xyyyyzz, g_x_0_yzzzz_xyyyzz, g_x_0_yzzzz_xyyyzzz, g_x_0_yzzzz_xyyzzz, g_x_0_yzzzz_xyyzzzz, g_x_0_yzzzz_xyzzzz, g_x_0_yzzzz_xyzzzzz, g_x_0_yzzzz_xzzzzz, g_x_0_yzzzz_yyyyyy, g_x_0_yzzzz_yyyyyyy, g_x_0_yzzzz_yyyyyyz, g_x_0_yzzzz_yyyyyz, g_x_0_yzzzz_yyyyyzz, g_x_0_yzzzz_yyyyzz, g_x_0_yzzzz_yyyyzzz, g_x_0_yzzzz_yyyzzz, g_x_0_yzzzz_yyyzzzz, g_x_0_yzzzz_yyzzzz, g_x_0_yzzzz_yyzzzzz, g_x_0_yzzzz_yzzzzz, g_x_0_yzzzz_yzzzzzz, g_x_0_yzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_xxxxxx[k] = -g_x_0_yzzzz_xxxxxx[k] * ab_y + g_x_0_yzzzz_xxxxxxy[k];

                g_x_0_yyzzzz_xxxxxy[k] = -g_x_0_yzzzz_xxxxxy[k] * ab_y + g_x_0_yzzzz_xxxxxyy[k];

                g_x_0_yyzzzz_xxxxxz[k] = -g_x_0_yzzzz_xxxxxz[k] * ab_y + g_x_0_yzzzz_xxxxxyz[k];

                g_x_0_yyzzzz_xxxxyy[k] = -g_x_0_yzzzz_xxxxyy[k] * ab_y + g_x_0_yzzzz_xxxxyyy[k];

                g_x_0_yyzzzz_xxxxyz[k] = -g_x_0_yzzzz_xxxxyz[k] * ab_y + g_x_0_yzzzz_xxxxyyz[k];

                g_x_0_yyzzzz_xxxxzz[k] = -g_x_0_yzzzz_xxxxzz[k] * ab_y + g_x_0_yzzzz_xxxxyzz[k];

                g_x_0_yyzzzz_xxxyyy[k] = -g_x_0_yzzzz_xxxyyy[k] * ab_y + g_x_0_yzzzz_xxxyyyy[k];

                g_x_0_yyzzzz_xxxyyz[k] = -g_x_0_yzzzz_xxxyyz[k] * ab_y + g_x_0_yzzzz_xxxyyyz[k];

                g_x_0_yyzzzz_xxxyzz[k] = -g_x_0_yzzzz_xxxyzz[k] * ab_y + g_x_0_yzzzz_xxxyyzz[k];

                g_x_0_yyzzzz_xxxzzz[k] = -g_x_0_yzzzz_xxxzzz[k] * ab_y + g_x_0_yzzzz_xxxyzzz[k];

                g_x_0_yyzzzz_xxyyyy[k] = -g_x_0_yzzzz_xxyyyy[k] * ab_y + g_x_0_yzzzz_xxyyyyy[k];

                g_x_0_yyzzzz_xxyyyz[k] = -g_x_0_yzzzz_xxyyyz[k] * ab_y + g_x_0_yzzzz_xxyyyyz[k];

                g_x_0_yyzzzz_xxyyzz[k] = -g_x_0_yzzzz_xxyyzz[k] * ab_y + g_x_0_yzzzz_xxyyyzz[k];

                g_x_0_yyzzzz_xxyzzz[k] = -g_x_0_yzzzz_xxyzzz[k] * ab_y + g_x_0_yzzzz_xxyyzzz[k];

                g_x_0_yyzzzz_xxzzzz[k] = -g_x_0_yzzzz_xxzzzz[k] * ab_y + g_x_0_yzzzz_xxyzzzz[k];

                g_x_0_yyzzzz_xyyyyy[k] = -g_x_0_yzzzz_xyyyyy[k] * ab_y + g_x_0_yzzzz_xyyyyyy[k];

                g_x_0_yyzzzz_xyyyyz[k] = -g_x_0_yzzzz_xyyyyz[k] * ab_y + g_x_0_yzzzz_xyyyyyz[k];

                g_x_0_yyzzzz_xyyyzz[k] = -g_x_0_yzzzz_xyyyzz[k] * ab_y + g_x_0_yzzzz_xyyyyzz[k];

                g_x_0_yyzzzz_xyyzzz[k] = -g_x_0_yzzzz_xyyzzz[k] * ab_y + g_x_0_yzzzz_xyyyzzz[k];

                g_x_0_yyzzzz_xyzzzz[k] = -g_x_0_yzzzz_xyzzzz[k] * ab_y + g_x_0_yzzzz_xyyzzzz[k];

                g_x_0_yyzzzz_xzzzzz[k] = -g_x_0_yzzzz_xzzzzz[k] * ab_y + g_x_0_yzzzz_xyzzzzz[k];

                g_x_0_yyzzzz_yyyyyy[k] = -g_x_0_yzzzz_yyyyyy[k] * ab_y + g_x_0_yzzzz_yyyyyyy[k];

                g_x_0_yyzzzz_yyyyyz[k] = -g_x_0_yzzzz_yyyyyz[k] * ab_y + g_x_0_yzzzz_yyyyyyz[k];

                g_x_0_yyzzzz_yyyyzz[k] = -g_x_0_yzzzz_yyyyzz[k] * ab_y + g_x_0_yzzzz_yyyyyzz[k];

                g_x_0_yyzzzz_yyyzzz[k] = -g_x_0_yzzzz_yyyzzz[k] * ab_y + g_x_0_yzzzz_yyyyzzz[k];

                g_x_0_yyzzzz_yyzzzz[k] = -g_x_0_yzzzz_yyzzzz[k] * ab_y + g_x_0_yzzzz_yyyzzzz[k];

                g_x_0_yyzzzz_yzzzzz[k] = -g_x_0_yzzzz_yzzzzz[k] * ab_y + g_x_0_yzzzz_yyzzzzz[k];

                g_x_0_yyzzzz_zzzzzz[k] = -g_x_0_yzzzz_zzzzzz[k] * ab_y + g_x_0_yzzzz_yzzzzzz[k];
            }

            /// Set up 728-756 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 728 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 729 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 730 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 731 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 732 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 733 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 734 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 735 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 736 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 737 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 738 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 739 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 740 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 741 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 742 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 743 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 744 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 745 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 746 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 747 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 748 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 749 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 750 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 751 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 752 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 753 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 754 * ccomps * dcomps);

            auto g_x_0_yzzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzzzz_xxxxxx, g_x_0_yzzzzz_xxxxxy, g_x_0_yzzzzz_xxxxxz, g_x_0_yzzzzz_xxxxyy, g_x_0_yzzzzz_xxxxyz, g_x_0_yzzzzz_xxxxzz, g_x_0_yzzzzz_xxxyyy, g_x_0_yzzzzz_xxxyyz, g_x_0_yzzzzz_xxxyzz, g_x_0_yzzzzz_xxxzzz, g_x_0_yzzzzz_xxyyyy, g_x_0_yzzzzz_xxyyyz, g_x_0_yzzzzz_xxyyzz, g_x_0_yzzzzz_xxyzzz, g_x_0_yzzzzz_xxzzzz, g_x_0_yzzzzz_xyyyyy, g_x_0_yzzzzz_xyyyyz, g_x_0_yzzzzz_xyyyzz, g_x_0_yzzzzz_xyyzzz, g_x_0_yzzzzz_xyzzzz, g_x_0_yzzzzz_xzzzzz, g_x_0_yzzzzz_yyyyyy, g_x_0_yzzzzz_yyyyyz, g_x_0_yzzzzz_yyyyzz, g_x_0_yzzzzz_yyyzzz, g_x_0_yzzzzz_yyzzzz, g_x_0_yzzzzz_yzzzzz, g_x_0_yzzzzz_zzzzzz, g_x_0_zzzzz_xxxxxx, g_x_0_zzzzz_xxxxxxy, g_x_0_zzzzz_xxxxxy, g_x_0_zzzzz_xxxxxyy, g_x_0_zzzzz_xxxxxyz, g_x_0_zzzzz_xxxxxz, g_x_0_zzzzz_xxxxyy, g_x_0_zzzzz_xxxxyyy, g_x_0_zzzzz_xxxxyyz, g_x_0_zzzzz_xxxxyz, g_x_0_zzzzz_xxxxyzz, g_x_0_zzzzz_xxxxzz, g_x_0_zzzzz_xxxyyy, g_x_0_zzzzz_xxxyyyy, g_x_0_zzzzz_xxxyyyz, g_x_0_zzzzz_xxxyyz, g_x_0_zzzzz_xxxyyzz, g_x_0_zzzzz_xxxyzz, g_x_0_zzzzz_xxxyzzz, g_x_0_zzzzz_xxxzzz, g_x_0_zzzzz_xxyyyy, g_x_0_zzzzz_xxyyyyy, g_x_0_zzzzz_xxyyyyz, g_x_0_zzzzz_xxyyyz, g_x_0_zzzzz_xxyyyzz, g_x_0_zzzzz_xxyyzz, g_x_0_zzzzz_xxyyzzz, g_x_0_zzzzz_xxyzzz, g_x_0_zzzzz_xxyzzzz, g_x_0_zzzzz_xxzzzz, g_x_0_zzzzz_xyyyyy, g_x_0_zzzzz_xyyyyyy, g_x_0_zzzzz_xyyyyyz, g_x_0_zzzzz_xyyyyz, g_x_0_zzzzz_xyyyyzz, g_x_0_zzzzz_xyyyzz, g_x_0_zzzzz_xyyyzzz, g_x_0_zzzzz_xyyzzz, g_x_0_zzzzz_xyyzzzz, g_x_0_zzzzz_xyzzzz, g_x_0_zzzzz_xyzzzzz, g_x_0_zzzzz_xzzzzz, g_x_0_zzzzz_yyyyyy, g_x_0_zzzzz_yyyyyyy, g_x_0_zzzzz_yyyyyyz, g_x_0_zzzzz_yyyyyz, g_x_0_zzzzz_yyyyyzz, g_x_0_zzzzz_yyyyzz, g_x_0_zzzzz_yyyyzzz, g_x_0_zzzzz_yyyzzz, g_x_0_zzzzz_yyyzzzz, g_x_0_zzzzz_yyzzzz, g_x_0_zzzzz_yyzzzzz, g_x_0_zzzzz_yzzzzz, g_x_0_zzzzz_yzzzzzz, g_x_0_zzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_xxxxxx[k] = -g_x_0_zzzzz_xxxxxx[k] * ab_y + g_x_0_zzzzz_xxxxxxy[k];

                g_x_0_yzzzzz_xxxxxy[k] = -g_x_0_zzzzz_xxxxxy[k] * ab_y + g_x_0_zzzzz_xxxxxyy[k];

                g_x_0_yzzzzz_xxxxxz[k] = -g_x_0_zzzzz_xxxxxz[k] * ab_y + g_x_0_zzzzz_xxxxxyz[k];

                g_x_0_yzzzzz_xxxxyy[k] = -g_x_0_zzzzz_xxxxyy[k] * ab_y + g_x_0_zzzzz_xxxxyyy[k];

                g_x_0_yzzzzz_xxxxyz[k] = -g_x_0_zzzzz_xxxxyz[k] * ab_y + g_x_0_zzzzz_xxxxyyz[k];

                g_x_0_yzzzzz_xxxxzz[k] = -g_x_0_zzzzz_xxxxzz[k] * ab_y + g_x_0_zzzzz_xxxxyzz[k];

                g_x_0_yzzzzz_xxxyyy[k] = -g_x_0_zzzzz_xxxyyy[k] * ab_y + g_x_0_zzzzz_xxxyyyy[k];

                g_x_0_yzzzzz_xxxyyz[k] = -g_x_0_zzzzz_xxxyyz[k] * ab_y + g_x_0_zzzzz_xxxyyyz[k];

                g_x_0_yzzzzz_xxxyzz[k] = -g_x_0_zzzzz_xxxyzz[k] * ab_y + g_x_0_zzzzz_xxxyyzz[k];

                g_x_0_yzzzzz_xxxzzz[k] = -g_x_0_zzzzz_xxxzzz[k] * ab_y + g_x_0_zzzzz_xxxyzzz[k];

                g_x_0_yzzzzz_xxyyyy[k] = -g_x_0_zzzzz_xxyyyy[k] * ab_y + g_x_0_zzzzz_xxyyyyy[k];

                g_x_0_yzzzzz_xxyyyz[k] = -g_x_0_zzzzz_xxyyyz[k] * ab_y + g_x_0_zzzzz_xxyyyyz[k];

                g_x_0_yzzzzz_xxyyzz[k] = -g_x_0_zzzzz_xxyyzz[k] * ab_y + g_x_0_zzzzz_xxyyyzz[k];

                g_x_0_yzzzzz_xxyzzz[k] = -g_x_0_zzzzz_xxyzzz[k] * ab_y + g_x_0_zzzzz_xxyyzzz[k];

                g_x_0_yzzzzz_xxzzzz[k] = -g_x_0_zzzzz_xxzzzz[k] * ab_y + g_x_0_zzzzz_xxyzzzz[k];

                g_x_0_yzzzzz_xyyyyy[k] = -g_x_0_zzzzz_xyyyyy[k] * ab_y + g_x_0_zzzzz_xyyyyyy[k];

                g_x_0_yzzzzz_xyyyyz[k] = -g_x_0_zzzzz_xyyyyz[k] * ab_y + g_x_0_zzzzz_xyyyyyz[k];

                g_x_0_yzzzzz_xyyyzz[k] = -g_x_0_zzzzz_xyyyzz[k] * ab_y + g_x_0_zzzzz_xyyyyzz[k];

                g_x_0_yzzzzz_xyyzzz[k] = -g_x_0_zzzzz_xyyzzz[k] * ab_y + g_x_0_zzzzz_xyyyzzz[k];

                g_x_0_yzzzzz_xyzzzz[k] = -g_x_0_zzzzz_xyzzzz[k] * ab_y + g_x_0_zzzzz_xyyzzzz[k];

                g_x_0_yzzzzz_xzzzzz[k] = -g_x_0_zzzzz_xzzzzz[k] * ab_y + g_x_0_zzzzz_xyzzzzz[k];

                g_x_0_yzzzzz_yyyyyy[k] = -g_x_0_zzzzz_yyyyyy[k] * ab_y + g_x_0_zzzzz_yyyyyyy[k];

                g_x_0_yzzzzz_yyyyyz[k] = -g_x_0_zzzzz_yyyyyz[k] * ab_y + g_x_0_zzzzz_yyyyyyz[k];

                g_x_0_yzzzzz_yyyyzz[k] = -g_x_0_zzzzz_yyyyzz[k] * ab_y + g_x_0_zzzzz_yyyyyzz[k];

                g_x_0_yzzzzz_yyyzzz[k] = -g_x_0_zzzzz_yyyzzz[k] * ab_y + g_x_0_zzzzz_yyyyzzz[k];

                g_x_0_yzzzzz_yyzzzz[k] = -g_x_0_zzzzz_yyzzzz[k] * ab_y + g_x_0_zzzzz_yyyzzzz[k];

                g_x_0_yzzzzz_yzzzzz[k] = -g_x_0_zzzzz_yzzzzz[k] * ab_y + g_x_0_zzzzz_yyzzzzz[k];

                g_x_0_yzzzzz_zzzzzz[k] = -g_x_0_zzzzz_zzzzzz[k] * ab_y + g_x_0_zzzzz_yzzzzzz[k];
            }

            /// Set up 756-784 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 756 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 757 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 758 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 759 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 760 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 761 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 762 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 763 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 764 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 765 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 766 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 767 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 768 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 769 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 770 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 771 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 772 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 773 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 774 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 775 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 776 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 777 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 778 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 779 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 780 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 781 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 782 * ccomps * dcomps);

            auto g_x_0_zzzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 783 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzzz_xxxxxx, g_x_0_zzzzz_xxxxxxz, g_x_0_zzzzz_xxxxxy, g_x_0_zzzzz_xxxxxyz, g_x_0_zzzzz_xxxxxz, g_x_0_zzzzz_xxxxxzz, g_x_0_zzzzz_xxxxyy, g_x_0_zzzzz_xxxxyyz, g_x_0_zzzzz_xxxxyz, g_x_0_zzzzz_xxxxyzz, g_x_0_zzzzz_xxxxzz, g_x_0_zzzzz_xxxxzzz, g_x_0_zzzzz_xxxyyy, g_x_0_zzzzz_xxxyyyz, g_x_0_zzzzz_xxxyyz, g_x_0_zzzzz_xxxyyzz, g_x_0_zzzzz_xxxyzz, g_x_0_zzzzz_xxxyzzz, g_x_0_zzzzz_xxxzzz, g_x_0_zzzzz_xxxzzzz, g_x_0_zzzzz_xxyyyy, g_x_0_zzzzz_xxyyyyz, g_x_0_zzzzz_xxyyyz, g_x_0_zzzzz_xxyyyzz, g_x_0_zzzzz_xxyyzz, g_x_0_zzzzz_xxyyzzz, g_x_0_zzzzz_xxyzzz, g_x_0_zzzzz_xxyzzzz, g_x_0_zzzzz_xxzzzz, g_x_0_zzzzz_xxzzzzz, g_x_0_zzzzz_xyyyyy, g_x_0_zzzzz_xyyyyyz, g_x_0_zzzzz_xyyyyz, g_x_0_zzzzz_xyyyyzz, g_x_0_zzzzz_xyyyzz, g_x_0_zzzzz_xyyyzzz, g_x_0_zzzzz_xyyzzz, g_x_0_zzzzz_xyyzzzz, g_x_0_zzzzz_xyzzzz, g_x_0_zzzzz_xyzzzzz, g_x_0_zzzzz_xzzzzz, g_x_0_zzzzz_xzzzzzz, g_x_0_zzzzz_yyyyyy, g_x_0_zzzzz_yyyyyyz, g_x_0_zzzzz_yyyyyz, g_x_0_zzzzz_yyyyyzz, g_x_0_zzzzz_yyyyzz, g_x_0_zzzzz_yyyyzzz, g_x_0_zzzzz_yyyzzz, g_x_0_zzzzz_yyyzzzz, g_x_0_zzzzz_yyzzzz, g_x_0_zzzzz_yyzzzzz, g_x_0_zzzzz_yzzzzz, g_x_0_zzzzz_yzzzzzz, g_x_0_zzzzz_zzzzzz, g_x_0_zzzzz_zzzzzzz, g_x_0_zzzzzz_xxxxxx, g_x_0_zzzzzz_xxxxxy, g_x_0_zzzzzz_xxxxxz, g_x_0_zzzzzz_xxxxyy, g_x_0_zzzzzz_xxxxyz, g_x_0_zzzzzz_xxxxzz, g_x_0_zzzzzz_xxxyyy, g_x_0_zzzzzz_xxxyyz, g_x_0_zzzzzz_xxxyzz, g_x_0_zzzzzz_xxxzzz, g_x_0_zzzzzz_xxyyyy, g_x_0_zzzzzz_xxyyyz, g_x_0_zzzzzz_xxyyzz, g_x_0_zzzzzz_xxyzzz, g_x_0_zzzzzz_xxzzzz, g_x_0_zzzzzz_xyyyyy, g_x_0_zzzzzz_xyyyyz, g_x_0_zzzzzz_xyyyzz, g_x_0_zzzzzz_xyyzzz, g_x_0_zzzzzz_xyzzzz, g_x_0_zzzzzz_xzzzzz, g_x_0_zzzzzz_yyyyyy, g_x_0_zzzzzz_yyyyyz, g_x_0_zzzzzz_yyyyzz, g_x_0_zzzzzz_yyyzzz, g_x_0_zzzzzz_yyzzzz, g_x_0_zzzzzz_yzzzzz, g_x_0_zzzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_xxxxxx[k] = -g_x_0_zzzzz_xxxxxx[k] * ab_z + g_x_0_zzzzz_xxxxxxz[k];

                g_x_0_zzzzzz_xxxxxy[k] = -g_x_0_zzzzz_xxxxxy[k] * ab_z + g_x_0_zzzzz_xxxxxyz[k];

                g_x_0_zzzzzz_xxxxxz[k] = -g_x_0_zzzzz_xxxxxz[k] * ab_z + g_x_0_zzzzz_xxxxxzz[k];

                g_x_0_zzzzzz_xxxxyy[k] = -g_x_0_zzzzz_xxxxyy[k] * ab_z + g_x_0_zzzzz_xxxxyyz[k];

                g_x_0_zzzzzz_xxxxyz[k] = -g_x_0_zzzzz_xxxxyz[k] * ab_z + g_x_0_zzzzz_xxxxyzz[k];

                g_x_0_zzzzzz_xxxxzz[k] = -g_x_0_zzzzz_xxxxzz[k] * ab_z + g_x_0_zzzzz_xxxxzzz[k];

                g_x_0_zzzzzz_xxxyyy[k] = -g_x_0_zzzzz_xxxyyy[k] * ab_z + g_x_0_zzzzz_xxxyyyz[k];

                g_x_0_zzzzzz_xxxyyz[k] = -g_x_0_zzzzz_xxxyyz[k] * ab_z + g_x_0_zzzzz_xxxyyzz[k];

                g_x_0_zzzzzz_xxxyzz[k] = -g_x_0_zzzzz_xxxyzz[k] * ab_z + g_x_0_zzzzz_xxxyzzz[k];

                g_x_0_zzzzzz_xxxzzz[k] = -g_x_0_zzzzz_xxxzzz[k] * ab_z + g_x_0_zzzzz_xxxzzzz[k];

                g_x_0_zzzzzz_xxyyyy[k] = -g_x_0_zzzzz_xxyyyy[k] * ab_z + g_x_0_zzzzz_xxyyyyz[k];

                g_x_0_zzzzzz_xxyyyz[k] = -g_x_0_zzzzz_xxyyyz[k] * ab_z + g_x_0_zzzzz_xxyyyzz[k];

                g_x_0_zzzzzz_xxyyzz[k] = -g_x_0_zzzzz_xxyyzz[k] * ab_z + g_x_0_zzzzz_xxyyzzz[k];

                g_x_0_zzzzzz_xxyzzz[k] = -g_x_0_zzzzz_xxyzzz[k] * ab_z + g_x_0_zzzzz_xxyzzzz[k];

                g_x_0_zzzzzz_xxzzzz[k] = -g_x_0_zzzzz_xxzzzz[k] * ab_z + g_x_0_zzzzz_xxzzzzz[k];

                g_x_0_zzzzzz_xyyyyy[k] = -g_x_0_zzzzz_xyyyyy[k] * ab_z + g_x_0_zzzzz_xyyyyyz[k];

                g_x_0_zzzzzz_xyyyyz[k] = -g_x_0_zzzzz_xyyyyz[k] * ab_z + g_x_0_zzzzz_xyyyyzz[k];

                g_x_0_zzzzzz_xyyyzz[k] = -g_x_0_zzzzz_xyyyzz[k] * ab_z + g_x_0_zzzzz_xyyyzzz[k];

                g_x_0_zzzzzz_xyyzzz[k] = -g_x_0_zzzzz_xyyzzz[k] * ab_z + g_x_0_zzzzz_xyyzzzz[k];

                g_x_0_zzzzzz_xyzzzz[k] = -g_x_0_zzzzz_xyzzzz[k] * ab_z + g_x_0_zzzzz_xyzzzzz[k];

                g_x_0_zzzzzz_xzzzzz[k] = -g_x_0_zzzzz_xzzzzz[k] * ab_z + g_x_0_zzzzz_xzzzzzz[k];

                g_x_0_zzzzzz_yyyyyy[k] = -g_x_0_zzzzz_yyyyyy[k] * ab_z + g_x_0_zzzzz_yyyyyyz[k];

                g_x_0_zzzzzz_yyyyyz[k] = -g_x_0_zzzzz_yyyyyz[k] * ab_z + g_x_0_zzzzz_yyyyyzz[k];

                g_x_0_zzzzzz_yyyyzz[k] = -g_x_0_zzzzz_yyyyzz[k] * ab_z + g_x_0_zzzzz_yyyyzzz[k];

                g_x_0_zzzzzz_yyyzzz[k] = -g_x_0_zzzzz_yyyzzz[k] * ab_z + g_x_0_zzzzz_yyyzzzz[k];

                g_x_0_zzzzzz_yyzzzz[k] = -g_x_0_zzzzz_yyzzzz[k] * ab_z + g_x_0_zzzzz_yyzzzzz[k];

                g_x_0_zzzzzz_yzzzzz[k] = -g_x_0_zzzzz_yzzzzz[k] * ab_z + g_x_0_zzzzz_yzzzzzz[k];

                g_x_0_zzzzzz_zzzzzz[k] = -g_x_0_zzzzz_zzzzzz[k] * ab_z + g_x_0_zzzzz_zzzzzzz[k];
            }

            /// Set up 784-812 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_xxxxxx = cbuffer.data(ii_geom_10_off + 784 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxxxy = cbuffer.data(ii_geom_10_off + 785 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxxxz = cbuffer.data(ii_geom_10_off + 786 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxxyy = cbuffer.data(ii_geom_10_off + 787 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxxyz = cbuffer.data(ii_geom_10_off + 788 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxxzz = cbuffer.data(ii_geom_10_off + 789 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxyyy = cbuffer.data(ii_geom_10_off + 790 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxyyz = cbuffer.data(ii_geom_10_off + 791 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxyzz = cbuffer.data(ii_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxzzz = cbuffer.data(ii_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxyyyy = cbuffer.data(ii_geom_10_off + 794 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxyyyz = cbuffer.data(ii_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxyyzz = cbuffer.data(ii_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxyzzz = cbuffer.data(ii_geom_10_off + 797 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxzzzz = cbuffer.data(ii_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyyyyy = cbuffer.data(ii_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyyyyz = cbuffer.data(ii_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyyyzz = cbuffer.data(ii_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyyzzz = cbuffer.data(ii_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyzzzz = cbuffer.data(ii_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xzzzzz = cbuffer.data(ii_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyyyyy = cbuffer.data(ii_geom_10_off + 805 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyyyyz = cbuffer.data(ii_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyyyzz = cbuffer.data(ii_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyyzzz = cbuffer.data(ii_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyzzzz = cbuffer.data(ii_geom_10_off + 809 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yzzzzz = cbuffer.data(ii_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_xxxxxx_zzzzzz = cbuffer.data(ii_geom_10_off + 811 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxx_xxxxxx, g_y_0_xxxxx_xxxxxxx, g_y_0_xxxxx_xxxxxxy, g_y_0_xxxxx_xxxxxxz, g_y_0_xxxxx_xxxxxy, g_y_0_xxxxx_xxxxxyy, g_y_0_xxxxx_xxxxxyz, g_y_0_xxxxx_xxxxxz, g_y_0_xxxxx_xxxxxzz, g_y_0_xxxxx_xxxxyy, g_y_0_xxxxx_xxxxyyy, g_y_0_xxxxx_xxxxyyz, g_y_0_xxxxx_xxxxyz, g_y_0_xxxxx_xxxxyzz, g_y_0_xxxxx_xxxxzz, g_y_0_xxxxx_xxxxzzz, g_y_0_xxxxx_xxxyyy, g_y_0_xxxxx_xxxyyyy, g_y_0_xxxxx_xxxyyyz, g_y_0_xxxxx_xxxyyz, g_y_0_xxxxx_xxxyyzz, g_y_0_xxxxx_xxxyzz, g_y_0_xxxxx_xxxyzzz, g_y_0_xxxxx_xxxzzz, g_y_0_xxxxx_xxxzzzz, g_y_0_xxxxx_xxyyyy, g_y_0_xxxxx_xxyyyyy, g_y_0_xxxxx_xxyyyyz, g_y_0_xxxxx_xxyyyz, g_y_0_xxxxx_xxyyyzz, g_y_0_xxxxx_xxyyzz, g_y_0_xxxxx_xxyyzzz, g_y_0_xxxxx_xxyzzz, g_y_0_xxxxx_xxyzzzz, g_y_0_xxxxx_xxzzzz, g_y_0_xxxxx_xxzzzzz, g_y_0_xxxxx_xyyyyy, g_y_0_xxxxx_xyyyyyy, g_y_0_xxxxx_xyyyyyz, g_y_0_xxxxx_xyyyyz, g_y_0_xxxxx_xyyyyzz, g_y_0_xxxxx_xyyyzz, g_y_0_xxxxx_xyyyzzz, g_y_0_xxxxx_xyyzzz, g_y_0_xxxxx_xyyzzzz, g_y_0_xxxxx_xyzzzz, g_y_0_xxxxx_xyzzzzz, g_y_0_xxxxx_xzzzzz, g_y_0_xxxxx_xzzzzzz, g_y_0_xxxxx_yyyyyy, g_y_0_xxxxx_yyyyyz, g_y_0_xxxxx_yyyyzz, g_y_0_xxxxx_yyyzzz, g_y_0_xxxxx_yyzzzz, g_y_0_xxxxx_yzzzzz, g_y_0_xxxxx_zzzzzz, g_y_0_xxxxxx_xxxxxx, g_y_0_xxxxxx_xxxxxy, g_y_0_xxxxxx_xxxxxz, g_y_0_xxxxxx_xxxxyy, g_y_0_xxxxxx_xxxxyz, g_y_0_xxxxxx_xxxxzz, g_y_0_xxxxxx_xxxyyy, g_y_0_xxxxxx_xxxyyz, g_y_0_xxxxxx_xxxyzz, g_y_0_xxxxxx_xxxzzz, g_y_0_xxxxxx_xxyyyy, g_y_0_xxxxxx_xxyyyz, g_y_0_xxxxxx_xxyyzz, g_y_0_xxxxxx_xxyzzz, g_y_0_xxxxxx_xxzzzz, g_y_0_xxxxxx_xyyyyy, g_y_0_xxxxxx_xyyyyz, g_y_0_xxxxxx_xyyyzz, g_y_0_xxxxxx_xyyzzz, g_y_0_xxxxxx_xyzzzz, g_y_0_xxxxxx_xzzzzz, g_y_0_xxxxxx_yyyyyy, g_y_0_xxxxxx_yyyyyz, g_y_0_xxxxxx_yyyyzz, g_y_0_xxxxxx_yyyzzz, g_y_0_xxxxxx_yyzzzz, g_y_0_xxxxxx_yzzzzz, g_y_0_xxxxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_xxxxxx[k] = -g_y_0_xxxxx_xxxxxx[k] * ab_x + g_y_0_xxxxx_xxxxxxx[k];

                g_y_0_xxxxxx_xxxxxy[k] = -g_y_0_xxxxx_xxxxxy[k] * ab_x + g_y_0_xxxxx_xxxxxxy[k];

                g_y_0_xxxxxx_xxxxxz[k] = -g_y_0_xxxxx_xxxxxz[k] * ab_x + g_y_0_xxxxx_xxxxxxz[k];

                g_y_0_xxxxxx_xxxxyy[k] = -g_y_0_xxxxx_xxxxyy[k] * ab_x + g_y_0_xxxxx_xxxxxyy[k];

                g_y_0_xxxxxx_xxxxyz[k] = -g_y_0_xxxxx_xxxxyz[k] * ab_x + g_y_0_xxxxx_xxxxxyz[k];

                g_y_0_xxxxxx_xxxxzz[k] = -g_y_0_xxxxx_xxxxzz[k] * ab_x + g_y_0_xxxxx_xxxxxzz[k];

                g_y_0_xxxxxx_xxxyyy[k] = -g_y_0_xxxxx_xxxyyy[k] * ab_x + g_y_0_xxxxx_xxxxyyy[k];

                g_y_0_xxxxxx_xxxyyz[k] = -g_y_0_xxxxx_xxxyyz[k] * ab_x + g_y_0_xxxxx_xxxxyyz[k];

                g_y_0_xxxxxx_xxxyzz[k] = -g_y_0_xxxxx_xxxyzz[k] * ab_x + g_y_0_xxxxx_xxxxyzz[k];

                g_y_0_xxxxxx_xxxzzz[k] = -g_y_0_xxxxx_xxxzzz[k] * ab_x + g_y_0_xxxxx_xxxxzzz[k];

                g_y_0_xxxxxx_xxyyyy[k] = -g_y_0_xxxxx_xxyyyy[k] * ab_x + g_y_0_xxxxx_xxxyyyy[k];

                g_y_0_xxxxxx_xxyyyz[k] = -g_y_0_xxxxx_xxyyyz[k] * ab_x + g_y_0_xxxxx_xxxyyyz[k];

                g_y_0_xxxxxx_xxyyzz[k] = -g_y_0_xxxxx_xxyyzz[k] * ab_x + g_y_0_xxxxx_xxxyyzz[k];

                g_y_0_xxxxxx_xxyzzz[k] = -g_y_0_xxxxx_xxyzzz[k] * ab_x + g_y_0_xxxxx_xxxyzzz[k];

                g_y_0_xxxxxx_xxzzzz[k] = -g_y_0_xxxxx_xxzzzz[k] * ab_x + g_y_0_xxxxx_xxxzzzz[k];

                g_y_0_xxxxxx_xyyyyy[k] = -g_y_0_xxxxx_xyyyyy[k] * ab_x + g_y_0_xxxxx_xxyyyyy[k];

                g_y_0_xxxxxx_xyyyyz[k] = -g_y_0_xxxxx_xyyyyz[k] * ab_x + g_y_0_xxxxx_xxyyyyz[k];

                g_y_0_xxxxxx_xyyyzz[k] = -g_y_0_xxxxx_xyyyzz[k] * ab_x + g_y_0_xxxxx_xxyyyzz[k];

                g_y_0_xxxxxx_xyyzzz[k] = -g_y_0_xxxxx_xyyzzz[k] * ab_x + g_y_0_xxxxx_xxyyzzz[k];

                g_y_0_xxxxxx_xyzzzz[k] = -g_y_0_xxxxx_xyzzzz[k] * ab_x + g_y_0_xxxxx_xxyzzzz[k];

                g_y_0_xxxxxx_xzzzzz[k] = -g_y_0_xxxxx_xzzzzz[k] * ab_x + g_y_0_xxxxx_xxzzzzz[k];

                g_y_0_xxxxxx_yyyyyy[k] = -g_y_0_xxxxx_yyyyyy[k] * ab_x + g_y_0_xxxxx_xyyyyyy[k];

                g_y_0_xxxxxx_yyyyyz[k] = -g_y_0_xxxxx_yyyyyz[k] * ab_x + g_y_0_xxxxx_xyyyyyz[k];

                g_y_0_xxxxxx_yyyyzz[k] = -g_y_0_xxxxx_yyyyzz[k] * ab_x + g_y_0_xxxxx_xyyyyzz[k];

                g_y_0_xxxxxx_yyyzzz[k] = -g_y_0_xxxxx_yyyzzz[k] * ab_x + g_y_0_xxxxx_xyyyzzz[k];

                g_y_0_xxxxxx_yyzzzz[k] = -g_y_0_xxxxx_yyzzzz[k] * ab_x + g_y_0_xxxxx_xyyzzzz[k];

                g_y_0_xxxxxx_yzzzzz[k] = -g_y_0_xxxxx_yzzzzz[k] * ab_x + g_y_0_xxxxx_xyzzzzz[k];

                g_y_0_xxxxxx_zzzzzz[k] = -g_y_0_xxxxx_zzzzzz[k] * ab_x + g_y_0_xxxxx_xzzzzzz[k];
            }

            /// Set up 812-840 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_xxxxxx = cbuffer.data(ii_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxxxy = cbuffer.data(ii_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxxxz = cbuffer.data(ii_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxxyy = cbuffer.data(ii_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxxyz = cbuffer.data(ii_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxxzz = cbuffer.data(ii_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxyyy = cbuffer.data(ii_geom_10_off + 818 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxyyz = cbuffer.data(ii_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxyzz = cbuffer.data(ii_geom_10_off + 820 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxzzz = cbuffer.data(ii_geom_10_off + 821 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxyyyy = cbuffer.data(ii_geom_10_off + 822 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxyyyz = cbuffer.data(ii_geom_10_off + 823 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxyyzz = cbuffer.data(ii_geom_10_off + 824 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxyzzz = cbuffer.data(ii_geom_10_off + 825 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxzzzz = cbuffer.data(ii_geom_10_off + 826 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyyyyy = cbuffer.data(ii_geom_10_off + 827 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyyyyz = cbuffer.data(ii_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyyyzz = cbuffer.data(ii_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyyzzz = cbuffer.data(ii_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyzzzz = cbuffer.data(ii_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xzzzzz = cbuffer.data(ii_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyyyyy = cbuffer.data(ii_geom_10_off + 833 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyyyyz = cbuffer.data(ii_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyyyzz = cbuffer.data(ii_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyyzzz = cbuffer.data(ii_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyzzzz = cbuffer.data(ii_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yzzzzz = cbuffer.data(ii_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_xxxxxy_zzzzzz = cbuffer.data(ii_geom_10_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxxy_xxxxxx, g_y_0_xxxxxy_xxxxxy, g_y_0_xxxxxy_xxxxxz, g_y_0_xxxxxy_xxxxyy, g_y_0_xxxxxy_xxxxyz, g_y_0_xxxxxy_xxxxzz, g_y_0_xxxxxy_xxxyyy, g_y_0_xxxxxy_xxxyyz, g_y_0_xxxxxy_xxxyzz, g_y_0_xxxxxy_xxxzzz, g_y_0_xxxxxy_xxyyyy, g_y_0_xxxxxy_xxyyyz, g_y_0_xxxxxy_xxyyzz, g_y_0_xxxxxy_xxyzzz, g_y_0_xxxxxy_xxzzzz, g_y_0_xxxxxy_xyyyyy, g_y_0_xxxxxy_xyyyyz, g_y_0_xxxxxy_xyyyzz, g_y_0_xxxxxy_xyyzzz, g_y_0_xxxxxy_xyzzzz, g_y_0_xxxxxy_xzzzzz, g_y_0_xxxxxy_yyyyyy, g_y_0_xxxxxy_yyyyyz, g_y_0_xxxxxy_yyyyzz, g_y_0_xxxxxy_yyyzzz, g_y_0_xxxxxy_yyzzzz, g_y_0_xxxxxy_yzzzzz, g_y_0_xxxxxy_zzzzzz, g_y_0_xxxxy_xxxxxx, g_y_0_xxxxy_xxxxxxx, g_y_0_xxxxy_xxxxxxy, g_y_0_xxxxy_xxxxxxz, g_y_0_xxxxy_xxxxxy, g_y_0_xxxxy_xxxxxyy, g_y_0_xxxxy_xxxxxyz, g_y_0_xxxxy_xxxxxz, g_y_0_xxxxy_xxxxxzz, g_y_0_xxxxy_xxxxyy, g_y_0_xxxxy_xxxxyyy, g_y_0_xxxxy_xxxxyyz, g_y_0_xxxxy_xxxxyz, g_y_0_xxxxy_xxxxyzz, g_y_0_xxxxy_xxxxzz, g_y_0_xxxxy_xxxxzzz, g_y_0_xxxxy_xxxyyy, g_y_0_xxxxy_xxxyyyy, g_y_0_xxxxy_xxxyyyz, g_y_0_xxxxy_xxxyyz, g_y_0_xxxxy_xxxyyzz, g_y_0_xxxxy_xxxyzz, g_y_0_xxxxy_xxxyzzz, g_y_0_xxxxy_xxxzzz, g_y_0_xxxxy_xxxzzzz, g_y_0_xxxxy_xxyyyy, g_y_0_xxxxy_xxyyyyy, g_y_0_xxxxy_xxyyyyz, g_y_0_xxxxy_xxyyyz, g_y_0_xxxxy_xxyyyzz, g_y_0_xxxxy_xxyyzz, g_y_0_xxxxy_xxyyzzz, g_y_0_xxxxy_xxyzzz, g_y_0_xxxxy_xxyzzzz, g_y_0_xxxxy_xxzzzz, g_y_0_xxxxy_xxzzzzz, g_y_0_xxxxy_xyyyyy, g_y_0_xxxxy_xyyyyyy, g_y_0_xxxxy_xyyyyyz, g_y_0_xxxxy_xyyyyz, g_y_0_xxxxy_xyyyyzz, g_y_0_xxxxy_xyyyzz, g_y_0_xxxxy_xyyyzzz, g_y_0_xxxxy_xyyzzz, g_y_0_xxxxy_xyyzzzz, g_y_0_xxxxy_xyzzzz, g_y_0_xxxxy_xyzzzzz, g_y_0_xxxxy_xzzzzz, g_y_0_xxxxy_xzzzzzz, g_y_0_xxxxy_yyyyyy, g_y_0_xxxxy_yyyyyz, g_y_0_xxxxy_yyyyzz, g_y_0_xxxxy_yyyzzz, g_y_0_xxxxy_yyzzzz, g_y_0_xxxxy_yzzzzz, g_y_0_xxxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_xxxxxx[k] = -g_y_0_xxxxy_xxxxxx[k] * ab_x + g_y_0_xxxxy_xxxxxxx[k];

                g_y_0_xxxxxy_xxxxxy[k] = -g_y_0_xxxxy_xxxxxy[k] * ab_x + g_y_0_xxxxy_xxxxxxy[k];

                g_y_0_xxxxxy_xxxxxz[k] = -g_y_0_xxxxy_xxxxxz[k] * ab_x + g_y_0_xxxxy_xxxxxxz[k];

                g_y_0_xxxxxy_xxxxyy[k] = -g_y_0_xxxxy_xxxxyy[k] * ab_x + g_y_0_xxxxy_xxxxxyy[k];

                g_y_0_xxxxxy_xxxxyz[k] = -g_y_0_xxxxy_xxxxyz[k] * ab_x + g_y_0_xxxxy_xxxxxyz[k];

                g_y_0_xxxxxy_xxxxzz[k] = -g_y_0_xxxxy_xxxxzz[k] * ab_x + g_y_0_xxxxy_xxxxxzz[k];

                g_y_0_xxxxxy_xxxyyy[k] = -g_y_0_xxxxy_xxxyyy[k] * ab_x + g_y_0_xxxxy_xxxxyyy[k];

                g_y_0_xxxxxy_xxxyyz[k] = -g_y_0_xxxxy_xxxyyz[k] * ab_x + g_y_0_xxxxy_xxxxyyz[k];

                g_y_0_xxxxxy_xxxyzz[k] = -g_y_0_xxxxy_xxxyzz[k] * ab_x + g_y_0_xxxxy_xxxxyzz[k];

                g_y_0_xxxxxy_xxxzzz[k] = -g_y_0_xxxxy_xxxzzz[k] * ab_x + g_y_0_xxxxy_xxxxzzz[k];

                g_y_0_xxxxxy_xxyyyy[k] = -g_y_0_xxxxy_xxyyyy[k] * ab_x + g_y_0_xxxxy_xxxyyyy[k];

                g_y_0_xxxxxy_xxyyyz[k] = -g_y_0_xxxxy_xxyyyz[k] * ab_x + g_y_0_xxxxy_xxxyyyz[k];

                g_y_0_xxxxxy_xxyyzz[k] = -g_y_0_xxxxy_xxyyzz[k] * ab_x + g_y_0_xxxxy_xxxyyzz[k];

                g_y_0_xxxxxy_xxyzzz[k] = -g_y_0_xxxxy_xxyzzz[k] * ab_x + g_y_0_xxxxy_xxxyzzz[k];

                g_y_0_xxxxxy_xxzzzz[k] = -g_y_0_xxxxy_xxzzzz[k] * ab_x + g_y_0_xxxxy_xxxzzzz[k];

                g_y_0_xxxxxy_xyyyyy[k] = -g_y_0_xxxxy_xyyyyy[k] * ab_x + g_y_0_xxxxy_xxyyyyy[k];

                g_y_0_xxxxxy_xyyyyz[k] = -g_y_0_xxxxy_xyyyyz[k] * ab_x + g_y_0_xxxxy_xxyyyyz[k];

                g_y_0_xxxxxy_xyyyzz[k] = -g_y_0_xxxxy_xyyyzz[k] * ab_x + g_y_0_xxxxy_xxyyyzz[k];

                g_y_0_xxxxxy_xyyzzz[k] = -g_y_0_xxxxy_xyyzzz[k] * ab_x + g_y_0_xxxxy_xxyyzzz[k];

                g_y_0_xxxxxy_xyzzzz[k] = -g_y_0_xxxxy_xyzzzz[k] * ab_x + g_y_0_xxxxy_xxyzzzz[k];

                g_y_0_xxxxxy_xzzzzz[k] = -g_y_0_xxxxy_xzzzzz[k] * ab_x + g_y_0_xxxxy_xxzzzzz[k];

                g_y_0_xxxxxy_yyyyyy[k] = -g_y_0_xxxxy_yyyyyy[k] * ab_x + g_y_0_xxxxy_xyyyyyy[k];

                g_y_0_xxxxxy_yyyyyz[k] = -g_y_0_xxxxy_yyyyyz[k] * ab_x + g_y_0_xxxxy_xyyyyyz[k];

                g_y_0_xxxxxy_yyyyzz[k] = -g_y_0_xxxxy_yyyyzz[k] * ab_x + g_y_0_xxxxy_xyyyyzz[k];

                g_y_0_xxxxxy_yyyzzz[k] = -g_y_0_xxxxy_yyyzzz[k] * ab_x + g_y_0_xxxxy_xyyyzzz[k];

                g_y_0_xxxxxy_yyzzzz[k] = -g_y_0_xxxxy_yyzzzz[k] * ab_x + g_y_0_xxxxy_xyyzzzz[k];

                g_y_0_xxxxxy_yzzzzz[k] = -g_y_0_xxxxy_yzzzzz[k] * ab_x + g_y_0_xxxxy_xyzzzzz[k];

                g_y_0_xxxxxy_zzzzzz[k] = -g_y_0_xxxxy_zzzzzz[k] * ab_x + g_y_0_xxxxy_xzzzzzz[k];
            }

            /// Set up 840-868 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_xxxxxx = cbuffer.data(ii_geom_10_off + 840 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxxxy = cbuffer.data(ii_geom_10_off + 841 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxxxz = cbuffer.data(ii_geom_10_off + 842 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxxyy = cbuffer.data(ii_geom_10_off + 843 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxxyz = cbuffer.data(ii_geom_10_off + 844 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxxzz = cbuffer.data(ii_geom_10_off + 845 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxyyy = cbuffer.data(ii_geom_10_off + 846 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxyyz = cbuffer.data(ii_geom_10_off + 847 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxyzz = cbuffer.data(ii_geom_10_off + 848 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxzzz = cbuffer.data(ii_geom_10_off + 849 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxyyyy = cbuffer.data(ii_geom_10_off + 850 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxyyyz = cbuffer.data(ii_geom_10_off + 851 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxyyzz = cbuffer.data(ii_geom_10_off + 852 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxyzzz = cbuffer.data(ii_geom_10_off + 853 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxzzzz = cbuffer.data(ii_geom_10_off + 854 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyyyyy = cbuffer.data(ii_geom_10_off + 855 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyyyyz = cbuffer.data(ii_geom_10_off + 856 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyyyzz = cbuffer.data(ii_geom_10_off + 857 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyyzzz = cbuffer.data(ii_geom_10_off + 858 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyzzzz = cbuffer.data(ii_geom_10_off + 859 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xzzzzz = cbuffer.data(ii_geom_10_off + 860 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyyyyy = cbuffer.data(ii_geom_10_off + 861 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyyyyz = cbuffer.data(ii_geom_10_off + 862 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyyyzz = cbuffer.data(ii_geom_10_off + 863 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyyzzz = cbuffer.data(ii_geom_10_off + 864 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyzzzz = cbuffer.data(ii_geom_10_off + 865 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yzzzzz = cbuffer.data(ii_geom_10_off + 866 * ccomps * dcomps);

            auto g_y_0_xxxxxz_zzzzzz = cbuffer.data(ii_geom_10_off + 867 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxxz_xxxxxx, g_y_0_xxxxxz_xxxxxy, g_y_0_xxxxxz_xxxxxz, g_y_0_xxxxxz_xxxxyy, g_y_0_xxxxxz_xxxxyz, g_y_0_xxxxxz_xxxxzz, g_y_0_xxxxxz_xxxyyy, g_y_0_xxxxxz_xxxyyz, g_y_0_xxxxxz_xxxyzz, g_y_0_xxxxxz_xxxzzz, g_y_0_xxxxxz_xxyyyy, g_y_0_xxxxxz_xxyyyz, g_y_0_xxxxxz_xxyyzz, g_y_0_xxxxxz_xxyzzz, g_y_0_xxxxxz_xxzzzz, g_y_0_xxxxxz_xyyyyy, g_y_0_xxxxxz_xyyyyz, g_y_0_xxxxxz_xyyyzz, g_y_0_xxxxxz_xyyzzz, g_y_0_xxxxxz_xyzzzz, g_y_0_xxxxxz_xzzzzz, g_y_0_xxxxxz_yyyyyy, g_y_0_xxxxxz_yyyyyz, g_y_0_xxxxxz_yyyyzz, g_y_0_xxxxxz_yyyzzz, g_y_0_xxxxxz_yyzzzz, g_y_0_xxxxxz_yzzzzz, g_y_0_xxxxxz_zzzzzz, g_y_0_xxxxz_xxxxxx, g_y_0_xxxxz_xxxxxxx, g_y_0_xxxxz_xxxxxxy, g_y_0_xxxxz_xxxxxxz, g_y_0_xxxxz_xxxxxy, g_y_0_xxxxz_xxxxxyy, g_y_0_xxxxz_xxxxxyz, g_y_0_xxxxz_xxxxxz, g_y_0_xxxxz_xxxxxzz, g_y_0_xxxxz_xxxxyy, g_y_0_xxxxz_xxxxyyy, g_y_0_xxxxz_xxxxyyz, g_y_0_xxxxz_xxxxyz, g_y_0_xxxxz_xxxxyzz, g_y_0_xxxxz_xxxxzz, g_y_0_xxxxz_xxxxzzz, g_y_0_xxxxz_xxxyyy, g_y_0_xxxxz_xxxyyyy, g_y_0_xxxxz_xxxyyyz, g_y_0_xxxxz_xxxyyz, g_y_0_xxxxz_xxxyyzz, g_y_0_xxxxz_xxxyzz, g_y_0_xxxxz_xxxyzzz, g_y_0_xxxxz_xxxzzz, g_y_0_xxxxz_xxxzzzz, g_y_0_xxxxz_xxyyyy, g_y_0_xxxxz_xxyyyyy, g_y_0_xxxxz_xxyyyyz, g_y_0_xxxxz_xxyyyz, g_y_0_xxxxz_xxyyyzz, g_y_0_xxxxz_xxyyzz, g_y_0_xxxxz_xxyyzzz, g_y_0_xxxxz_xxyzzz, g_y_0_xxxxz_xxyzzzz, g_y_0_xxxxz_xxzzzz, g_y_0_xxxxz_xxzzzzz, g_y_0_xxxxz_xyyyyy, g_y_0_xxxxz_xyyyyyy, g_y_0_xxxxz_xyyyyyz, g_y_0_xxxxz_xyyyyz, g_y_0_xxxxz_xyyyyzz, g_y_0_xxxxz_xyyyzz, g_y_0_xxxxz_xyyyzzz, g_y_0_xxxxz_xyyzzz, g_y_0_xxxxz_xyyzzzz, g_y_0_xxxxz_xyzzzz, g_y_0_xxxxz_xyzzzzz, g_y_0_xxxxz_xzzzzz, g_y_0_xxxxz_xzzzzzz, g_y_0_xxxxz_yyyyyy, g_y_0_xxxxz_yyyyyz, g_y_0_xxxxz_yyyyzz, g_y_0_xxxxz_yyyzzz, g_y_0_xxxxz_yyzzzz, g_y_0_xxxxz_yzzzzz, g_y_0_xxxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_xxxxxx[k] = -g_y_0_xxxxz_xxxxxx[k] * ab_x + g_y_0_xxxxz_xxxxxxx[k];

                g_y_0_xxxxxz_xxxxxy[k] = -g_y_0_xxxxz_xxxxxy[k] * ab_x + g_y_0_xxxxz_xxxxxxy[k];

                g_y_0_xxxxxz_xxxxxz[k] = -g_y_0_xxxxz_xxxxxz[k] * ab_x + g_y_0_xxxxz_xxxxxxz[k];

                g_y_0_xxxxxz_xxxxyy[k] = -g_y_0_xxxxz_xxxxyy[k] * ab_x + g_y_0_xxxxz_xxxxxyy[k];

                g_y_0_xxxxxz_xxxxyz[k] = -g_y_0_xxxxz_xxxxyz[k] * ab_x + g_y_0_xxxxz_xxxxxyz[k];

                g_y_0_xxxxxz_xxxxzz[k] = -g_y_0_xxxxz_xxxxzz[k] * ab_x + g_y_0_xxxxz_xxxxxzz[k];

                g_y_0_xxxxxz_xxxyyy[k] = -g_y_0_xxxxz_xxxyyy[k] * ab_x + g_y_0_xxxxz_xxxxyyy[k];

                g_y_0_xxxxxz_xxxyyz[k] = -g_y_0_xxxxz_xxxyyz[k] * ab_x + g_y_0_xxxxz_xxxxyyz[k];

                g_y_0_xxxxxz_xxxyzz[k] = -g_y_0_xxxxz_xxxyzz[k] * ab_x + g_y_0_xxxxz_xxxxyzz[k];

                g_y_0_xxxxxz_xxxzzz[k] = -g_y_0_xxxxz_xxxzzz[k] * ab_x + g_y_0_xxxxz_xxxxzzz[k];

                g_y_0_xxxxxz_xxyyyy[k] = -g_y_0_xxxxz_xxyyyy[k] * ab_x + g_y_0_xxxxz_xxxyyyy[k];

                g_y_0_xxxxxz_xxyyyz[k] = -g_y_0_xxxxz_xxyyyz[k] * ab_x + g_y_0_xxxxz_xxxyyyz[k];

                g_y_0_xxxxxz_xxyyzz[k] = -g_y_0_xxxxz_xxyyzz[k] * ab_x + g_y_0_xxxxz_xxxyyzz[k];

                g_y_0_xxxxxz_xxyzzz[k] = -g_y_0_xxxxz_xxyzzz[k] * ab_x + g_y_0_xxxxz_xxxyzzz[k];

                g_y_0_xxxxxz_xxzzzz[k] = -g_y_0_xxxxz_xxzzzz[k] * ab_x + g_y_0_xxxxz_xxxzzzz[k];

                g_y_0_xxxxxz_xyyyyy[k] = -g_y_0_xxxxz_xyyyyy[k] * ab_x + g_y_0_xxxxz_xxyyyyy[k];

                g_y_0_xxxxxz_xyyyyz[k] = -g_y_0_xxxxz_xyyyyz[k] * ab_x + g_y_0_xxxxz_xxyyyyz[k];

                g_y_0_xxxxxz_xyyyzz[k] = -g_y_0_xxxxz_xyyyzz[k] * ab_x + g_y_0_xxxxz_xxyyyzz[k];

                g_y_0_xxxxxz_xyyzzz[k] = -g_y_0_xxxxz_xyyzzz[k] * ab_x + g_y_0_xxxxz_xxyyzzz[k];

                g_y_0_xxxxxz_xyzzzz[k] = -g_y_0_xxxxz_xyzzzz[k] * ab_x + g_y_0_xxxxz_xxyzzzz[k];

                g_y_0_xxxxxz_xzzzzz[k] = -g_y_0_xxxxz_xzzzzz[k] * ab_x + g_y_0_xxxxz_xxzzzzz[k];

                g_y_0_xxxxxz_yyyyyy[k] = -g_y_0_xxxxz_yyyyyy[k] * ab_x + g_y_0_xxxxz_xyyyyyy[k];

                g_y_0_xxxxxz_yyyyyz[k] = -g_y_0_xxxxz_yyyyyz[k] * ab_x + g_y_0_xxxxz_xyyyyyz[k];

                g_y_0_xxxxxz_yyyyzz[k] = -g_y_0_xxxxz_yyyyzz[k] * ab_x + g_y_0_xxxxz_xyyyyzz[k];

                g_y_0_xxxxxz_yyyzzz[k] = -g_y_0_xxxxz_yyyzzz[k] * ab_x + g_y_0_xxxxz_xyyyzzz[k];

                g_y_0_xxxxxz_yyzzzz[k] = -g_y_0_xxxxz_yyzzzz[k] * ab_x + g_y_0_xxxxz_xyyzzzz[k];

                g_y_0_xxxxxz_yzzzzz[k] = -g_y_0_xxxxz_yzzzzz[k] * ab_x + g_y_0_xxxxz_xyzzzzz[k];

                g_y_0_xxxxxz_zzzzzz[k] = -g_y_0_xxxxz_zzzzzz[k] * ab_x + g_y_0_xxxxz_xzzzzzz[k];
            }

            /// Set up 868-896 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_xxxxxx = cbuffer.data(ii_geom_10_off + 868 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxxxy = cbuffer.data(ii_geom_10_off + 869 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxxxz = cbuffer.data(ii_geom_10_off + 870 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxxyy = cbuffer.data(ii_geom_10_off + 871 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxxyz = cbuffer.data(ii_geom_10_off + 872 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxxzz = cbuffer.data(ii_geom_10_off + 873 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxyyy = cbuffer.data(ii_geom_10_off + 874 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxyyz = cbuffer.data(ii_geom_10_off + 875 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxyzz = cbuffer.data(ii_geom_10_off + 876 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxzzz = cbuffer.data(ii_geom_10_off + 877 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxyyyy = cbuffer.data(ii_geom_10_off + 878 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxyyyz = cbuffer.data(ii_geom_10_off + 879 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxyyzz = cbuffer.data(ii_geom_10_off + 880 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxyzzz = cbuffer.data(ii_geom_10_off + 881 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxzzzz = cbuffer.data(ii_geom_10_off + 882 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyyyyy = cbuffer.data(ii_geom_10_off + 883 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyyyyz = cbuffer.data(ii_geom_10_off + 884 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyyyzz = cbuffer.data(ii_geom_10_off + 885 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyyzzz = cbuffer.data(ii_geom_10_off + 886 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyzzzz = cbuffer.data(ii_geom_10_off + 887 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xzzzzz = cbuffer.data(ii_geom_10_off + 888 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyyyyy = cbuffer.data(ii_geom_10_off + 889 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyyyyz = cbuffer.data(ii_geom_10_off + 890 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyyyzz = cbuffer.data(ii_geom_10_off + 891 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyyzzz = cbuffer.data(ii_geom_10_off + 892 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyzzzz = cbuffer.data(ii_geom_10_off + 893 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yzzzzz = cbuffer.data(ii_geom_10_off + 894 * ccomps * dcomps);

            auto g_y_0_xxxxyy_zzzzzz = cbuffer.data(ii_geom_10_off + 895 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxyy_xxxxxx, g_y_0_xxxxyy_xxxxxy, g_y_0_xxxxyy_xxxxxz, g_y_0_xxxxyy_xxxxyy, g_y_0_xxxxyy_xxxxyz, g_y_0_xxxxyy_xxxxzz, g_y_0_xxxxyy_xxxyyy, g_y_0_xxxxyy_xxxyyz, g_y_0_xxxxyy_xxxyzz, g_y_0_xxxxyy_xxxzzz, g_y_0_xxxxyy_xxyyyy, g_y_0_xxxxyy_xxyyyz, g_y_0_xxxxyy_xxyyzz, g_y_0_xxxxyy_xxyzzz, g_y_0_xxxxyy_xxzzzz, g_y_0_xxxxyy_xyyyyy, g_y_0_xxxxyy_xyyyyz, g_y_0_xxxxyy_xyyyzz, g_y_0_xxxxyy_xyyzzz, g_y_0_xxxxyy_xyzzzz, g_y_0_xxxxyy_xzzzzz, g_y_0_xxxxyy_yyyyyy, g_y_0_xxxxyy_yyyyyz, g_y_0_xxxxyy_yyyyzz, g_y_0_xxxxyy_yyyzzz, g_y_0_xxxxyy_yyzzzz, g_y_0_xxxxyy_yzzzzz, g_y_0_xxxxyy_zzzzzz, g_y_0_xxxyy_xxxxxx, g_y_0_xxxyy_xxxxxxx, g_y_0_xxxyy_xxxxxxy, g_y_0_xxxyy_xxxxxxz, g_y_0_xxxyy_xxxxxy, g_y_0_xxxyy_xxxxxyy, g_y_0_xxxyy_xxxxxyz, g_y_0_xxxyy_xxxxxz, g_y_0_xxxyy_xxxxxzz, g_y_0_xxxyy_xxxxyy, g_y_0_xxxyy_xxxxyyy, g_y_0_xxxyy_xxxxyyz, g_y_0_xxxyy_xxxxyz, g_y_0_xxxyy_xxxxyzz, g_y_0_xxxyy_xxxxzz, g_y_0_xxxyy_xxxxzzz, g_y_0_xxxyy_xxxyyy, g_y_0_xxxyy_xxxyyyy, g_y_0_xxxyy_xxxyyyz, g_y_0_xxxyy_xxxyyz, g_y_0_xxxyy_xxxyyzz, g_y_0_xxxyy_xxxyzz, g_y_0_xxxyy_xxxyzzz, g_y_0_xxxyy_xxxzzz, g_y_0_xxxyy_xxxzzzz, g_y_0_xxxyy_xxyyyy, g_y_0_xxxyy_xxyyyyy, g_y_0_xxxyy_xxyyyyz, g_y_0_xxxyy_xxyyyz, g_y_0_xxxyy_xxyyyzz, g_y_0_xxxyy_xxyyzz, g_y_0_xxxyy_xxyyzzz, g_y_0_xxxyy_xxyzzz, g_y_0_xxxyy_xxyzzzz, g_y_0_xxxyy_xxzzzz, g_y_0_xxxyy_xxzzzzz, g_y_0_xxxyy_xyyyyy, g_y_0_xxxyy_xyyyyyy, g_y_0_xxxyy_xyyyyyz, g_y_0_xxxyy_xyyyyz, g_y_0_xxxyy_xyyyyzz, g_y_0_xxxyy_xyyyzz, g_y_0_xxxyy_xyyyzzz, g_y_0_xxxyy_xyyzzz, g_y_0_xxxyy_xyyzzzz, g_y_0_xxxyy_xyzzzz, g_y_0_xxxyy_xyzzzzz, g_y_0_xxxyy_xzzzzz, g_y_0_xxxyy_xzzzzzz, g_y_0_xxxyy_yyyyyy, g_y_0_xxxyy_yyyyyz, g_y_0_xxxyy_yyyyzz, g_y_0_xxxyy_yyyzzz, g_y_0_xxxyy_yyzzzz, g_y_0_xxxyy_yzzzzz, g_y_0_xxxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_xxxxxx[k] = -g_y_0_xxxyy_xxxxxx[k] * ab_x + g_y_0_xxxyy_xxxxxxx[k];

                g_y_0_xxxxyy_xxxxxy[k] = -g_y_0_xxxyy_xxxxxy[k] * ab_x + g_y_0_xxxyy_xxxxxxy[k];

                g_y_0_xxxxyy_xxxxxz[k] = -g_y_0_xxxyy_xxxxxz[k] * ab_x + g_y_0_xxxyy_xxxxxxz[k];

                g_y_0_xxxxyy_xxxxyy[k] = -g_y_0_xxxyy_xxxxyy[k] * ab_x + g_y_0_xxxyy_xxxxxyy[k];

                g_y_0_xxxxyy_xxxxyz[k] = -g_y_0_xxxyy_xxxxyz[k] * ab_x + g_y_0_xxxyy_xxxxxyz[k];

                g_y_0_xxxxyy_xxxxzz[k] = -g_y_0_xxxyy_xxxxzz[k] * ab_x + g_y_0_xxxyy_xxxxxzz[k];

                g_y_0_xxxxyy_xxxyyy[k] = -g_y_0_xxxyy_xxxyyy[k] * ab_x + g_y_0_xxxyy_xxxxyyy[k];

                g_y_0_xxxxyy_xxxyyz[k] = -g_y_0_xxxyy_xxxyyz[k] * ab_x + g_y_0_xxxyy_xxxxyyz[k];

                g_y_0_xxxxyy_xxxyzz[k] = -g_y_0_xxxyy_xxxyzz[k] * ab_x + g_y_0_xxxyy_xxxxyzz[k];

                g_y_0_xxxxyy_xxxzzz[k] = -g_y_0_xxxyy_xxxzzz[k] * ab_x + g_y_0_xxxyy_xxxxzzz[k];

                g_y_0_xxxxyy_xxyyyy[k] = -g_y_0_xxxyy_xxyyyy[k] * ab_x + g_y_0_xxxyy_xxxyyyy[k];

                g_y_0_xxxxyy_xxyyyz[k] = -g_y_0_xxxyy_xxyyyz[k] * ab_x + g_y_0_xxxyy_xxxyyyz[k];

                g_y_0_xxxxyy_xxyyzz[k] = -g_y_0_xxxyy_xxyyzz[k] * ab_x + g_y_0_xxxyy_xxxyyzz[k];

                g_y_0_xxxxyy_xxyzzz[k] = -g_y_0_xxxyy_xxyzzz[k] * ab_x + g_y_0_xxxyy_xxxyzzz[k];

                g_y_0_xxxxyy_xxzzzz[k] = -g_y_0_xxxyy_xxzzzz[k] * ab_x + g_y_0_xxxyy_xxxzzzz[k];

                g_y_0_xxxxyy_xyyyyy[k] = -g_y_0_xxxyy_xyyyyy[k] * ab_x + g_y_0_xxxyy_xxyyyyy[k];

                g_y_0_xxxxyy_xyyyyz[k] = -g_y_0_xxxyy_xyyyyz[k] * ab_x + g_y_0_xxxyy_xxyyyyz[k];

                g_y_0_xxxxyy_xyyyzz[k] = -g_y_0_xxxyy_xyyyzz[k] * ab_x + g_y_0_xxxyy_xxyyyzz[k];

                g_y_0_xxxxyy_xyyzzz[k] = -g_y_0_xxxyy_xyyzzz[k] * ab_x + g_y_0_xxxyy_xxyyzzz[k];

                g_y_0_xxxxyy_xyzzzz[k] = -g_y_0_xxxyy_xyzzzz[k] * ab_x + g_y_0_xxxyy_xxyzzzz[k];

                g_y_0_xxxxyy_xzzzzz[k] = -g_y_0_xxxyy_xzzzzz[k] * ab_x + g_y_0_xxxyy_xxzzzzz[k];

                g_y_0_xxxxyy_yyyyyy[k] = -g_y_0_xxxyy_yyyyyy[k] * ab_x + g_y_0_xxxyy_xyyyyyy[k];

                g_y_0_xxxxyy_yyyyyz[k] = -g_y_0_xxxyy_yyyyyz[k] * ab_x + g_y_0_xxxyy_xyyyyyz[k];

                g_y_0_xxxxyy_yyyyzz[k] = -g_y_0_xxxyy_yyyyzz[k] * ab_x + g_y_0_xxxyy_xyyyyzz[k];

                g_y_0_xxxxyy_yyyzzz[k] = -g_y_0_xxxyy_yyyzzz[k] * ab_x + g_y_0_xxxyy_xyyyzzz[k];

                g_y_0_xxxxyy_yyzzzz[k] = -g_y_0_xxxyy_yyzzzz[k] * ab_x + g_y_0_xxxyy_xyyzzzz[k];

                g_y_0_xxxxyy_yzzzzz[k] = -g_y_0_xxxyy_yzzzzz[k] * ab_x + g_y_0_xxxyy_xyzzzzz[k];

                g_y_0_xxxxyy_zzzzzz[k] = -g_y_0_xxxyy_zzzzzz[k] * ab_x + g_y_0_xxxyy_xzzzzzz[k];
            }

            /// Set up 896-924 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_xxxxxx = cbuffer.data(ii_geom_10_off + 896 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxxxy = cbuffer.data(ii_geom_10_off + 897 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxxxz = cbuffer.data(ii_geom_10_off + 898 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxxyy = cbuffer.data(ii_geom_10_off + 899 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxxyz = cbuffer.data(ii_geom_10_off + 900 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxxzz = cbuffer.data(ii_geom_10_off + 901 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxyyy = cbuffer.data(ii_geom_10_off + 902 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxyyz = cbuffer.data(ii_geom_10_off + 903 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxyzz = cbuffer.data(ii_geom_10_off + 904 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxzzz = cbuffer.data(ii_geom_10_off + 905 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxyyyy = cbuffer.data(ii_geom_10_off + 906 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxyyyz = cbuffer.data(ii_geom_10_off + 907 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxyyzz = cbuffer.data(ii_geom_10_off + 908 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxyzzz = cbuffer.data(ii_geom_10_off + 909 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxzzzz = cbuffer.data(ii_geom_10_off + 910 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyyyyy = cbuffer.data(ii_geom_10_off + 911 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyyyyz = cbuffer.data(ii_geom_10_off + 912 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyyyzz = cbuffer.data(ii_geom_10_off + 913 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyyzzz = cbuffer.data(ii_geom_10_off + 914 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyzzzz = cbuffer.data(ii_geom_10_off + 915 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xzzzzz = cbuffer.data(ii_geom_10_off + 916 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyyyyy = cbuffer.data(ii_geom_10_off + 917 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyyyyz = cbuffer.data(ii_geom_10_off + 918 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyyyzz = cbuffer.data(ii_geom_10_off + 919 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyyzzz = cbuffer.data(ii_geom_10_off + 920 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyzzzz = cbuffer.data(ii_geom_10_off + 921 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yzzzzz = cbuffer.data(ii_geom_10_off + 922 * ccomps * dcomps);

            auto g_y_0_xxxxyz_zzzzzz = cbuffer.data(ii_geom_10_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxyz_xxxxxx, g_y_0_xxxxyz_xxxxxy, g_y_0_xxxxyz_xxxxxz, g_y_0_xxxxyz_xxxxyy, g_y_0_xxxxyz_xxxxyz, g_y_0_xxxxyz_xxxxzz, g_y_0_xxxxyz_xxxyyy, g_y_0_xxxxyz_xxxyyz, g_y_0_xxxxyz_xxxyzz, g_y_0_xxxxyz_xxxzzz, g_y_0_xxxxyz_xxyyyy, g_y_0_xxxxyz_xxyyyz, g_y_0_xxxxyz_xxyyzz, g_y_0_xxxxyz_xxyzzz, g_y_0_xxxxyz_xxzzzz, g_y_0_xxxxyz_xyyyyy, g_y_0_xxxxyz_xyyyyz, g_y_0_xxxxyz_xyyyzz, g_y_0_xxxxyz_xyyzzz, g_y_0_xxxxyz_xyzzzz, g_y_0_xxxxyz_xzzzzz, g_y_0_xxxxyz_yyyyyy, g_y_0_xxxxyz_yyyyyz, g_y_0_xxxxyz_yyyyzz, g_y_0_xxxxyz_yyyzzz, g_y_0_xxxxyz_yyzzzz, g_y_0_xxxxyz_yzzzzz, g_y_0_xxxxyz_zzzzzz, g_y_0_xxxyz_xxxxxx, g_y_0_xxxyz_xxxxxxx, g_y_0_xxxyz_xxxxxxy, g_y_0_xxxyz_xxxxxxz, g_y_0_xxxyz_xxxxxy, g_y_0_xxxyz_xxxxxyy, g_y_0_xxxyz_xxxxxyz, g_y_0_xxxyz_xxxxxz, g_y_0_xxxyz_xxxxxzz, g_y_0_xxxyz_xxxxyy, g_y_0_xxxyz_xxxxyyy, g_y_0_xxxyz_xxxxyyz, g_y_0_xxxyz_xxxxyz, g_y_0_xxxyz_xxxxyzz, g_y_0_xxxyz_xxxxzz, g_y_0_xxxyz_xxxxzzz, g_y_0_xxxyz_xxxyyy, g_y_0_xxxyz_xxxyyyy, g_y_0_xxxyz_xxxyyyz, g_y_0_xxxyz_xxxyyz, g_y_0_xxxyz_xxxyyzz, g_y_0_xxxyz_xxxyzz, g_y_0_xxxyz_xxxyzzz, g_y_0_xxxyz_xxxzzz, g_y_0_xxxyz_xxxzzzz, g_y_0_xxxyz_xxyyyy, g_y_0_xxxyz_xxyyyyy, g_y_0_xxxyz_xxyyyyz, g_y_0_xxxyz_xxyyyz, g_y_0_xxxyz_xxyyyzz, g_y_0_xxxyz_xxyyzz, g_y_0_xxxyz_xxyyzzz, g_y_0_xxxyz_xxyzzz, g_y_0_xxxyz_xxyzzzz, g_y_0_xxxyz_xxzzzz, g_y_0_xxxyz_xxzzzzz, g_y_0_xxxyz_xyyyyy, g_y_0_xxxyz_xyyyyyy, g_y_0_xxxyz_xyyyyyz, g_y_0_xxxyz_xyyyyz, g_y_0_xxxyz_xyyyyzz, g_y_0_xxxyz_xyyyzz, g_y_0_xxxyz_xyyyzzz, g_y_0_xxxyz_xyyzzz, g_y_0_xxxyz_xyyzzzz, g_y_0_xxxyz_xyzzzz, g_y_0_xxxyz_xyzzzzz, g_y_0_xxxyz_xzzzzz, g_y_0_xxxyz_xzzzzzz, g_y_0_xxxyz_yyyyyy, g_y_0_xxxyz_yyyyyz, g_y_0_xxxyz_yyyyzz, g_y_0_xxxyz_yyyzzz, g_y_0_xxxyz_yyzzzz, g_y_0_xxxyz_yzzzzz, g_y_0_xxxyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_xxxxxx[k] = -g_y_0_xxxyz_xxxxxx[k] * ab_x + g_y_0_xxxyz_xxxxxxx[k];

                g_y_0_xxxxyz_xxxxxy[k] = -g_y_0_xxxyz_xxxxxy[k] * ab_x + g_y_0_xxxyz_xxxxxxy[k];

                g_y_0_xxxxyz_xxxxxz[k] = -g_y_0_xxxyz_xxxxxz[k] * ab_x + g_y_0_xxxyz_xxxxxxz[k];

                g_y_0_xxxxyz_xxxxyy[k] = -g_y_0_xxxyz_xxxxyy[k] * ab_x + g_y_0_xxxyz_xxxxxyy[k];

                g_y_0_xxxxyz_xxxxyz[k] = -g_y_0_xxxyz_xxxxyz[k] * ab_x + g_y_0_xxxyz_xxxxxyz[k];

                g_y_0_xxxxyz_xxxxzz[k] = -g_y_0_xxxyz_xxxxzz[k] * ab_x + g_y_0_xxxyz_xxxxxzz[k];

                g_y_0_xxxxyz_xxxyyy[k] = -g_y_0_xxxyz_xxxyyy[k] * ab_x + g_y_0_xxxyz_xxxxyyy[k];

                g_y_0_xxxxyz_xxxyyz[k] = -g_y_0_xxxyz_xxxyyz[k] * ab_x + g_y_0_xxxyz_xxxxyyz[k];

                g_y_0_xxxxyz_xxxyzz[k] = -g_y_0_xxxyz_xxxyzz[k] * ab_x + g_y_0_xxxyz_xxxxyzz[k];

                g_y_0_xxxxyz_xxxzzz[k] = -g_y_0_xxxyz_xxxzzz[k] * ab_x + g_y_0_xxxyz_xxxxzzz[k];

                g_y_0_xxxxyz_xxyyyy[k] = -g_y_0_xxxyz_xxyyyy[k] * ab_x + g_y_0_xxxyz_xxxyyyy[k];

                g_y_0_xxxxyz_xxyyyz[k] = -g_y_0_xxxyz_xxyyyz[k] * ab_x + g_y_0_xxxyz_xxxyyyz[k];

                g_y_0_xxxxyz_xxyyzz[k] = -g_y_0_xxxyz_xxyyzz[k] * ab_x + g_y_0_xxxyz_xxxyyzz[k];

                g_y_0_xxxxyz_xxyzzz[k] = -g_y_0_xxxyz_xxyzzz[k] * ab_x + g_y_0_xxxyz_xxxyzzz[k];

                g_y_0_xxxxyz_xxzzzz[k] = -g_y_0_xxxyz_xxzzzz[k] * ab_x + g_y_0_xxxyz_xxxzzzz[k];

                g_y_0_xxxxyz_xyyyyy[k] = -g_y_0_xxxyz_xyyyyy[k] * ab_x + g_y_0_xxxyz_xxyyyyy[k];

                g_y_0_xxxxyz_xyyyyz[k] = -g_y_0_xxxyz_xyyyyz[k] * ab_x + g_y_0_xxxyz_xxyyyyz[k];

                g_y_0_xxxxyz_xyyyzz[k] = -g_y_0_xxxyz_xyyyzz[k] * ab_x + g_y_0_xxxyz_xxyyyzz[k];

                g_y_0_xxxxyz_xyyzzz[k] = -g_y_0_xxxyz_xyyzzz[k] * ab_x + g_y_0_xxxyz_xxyyzzz[k];

                g_y_0_xxxxyz_xyzzzz[k] = -g_y_0_xxxyz_xyzzzz[k] * ab_x + g_y_0_xxxyz_xxyzzzz[k];

                g_y_0_xxxxyz_xzzzzz[k] = -g_y_0_xxxyz_xzzzzz[k] * ab_x + g_y_0_xxxyz_xxzzzzz[k];

                g_y_0_xxxxyz_yyyyyy[k] = -g_y_0_xxxyz_yyyyyy[k] * ab_x + g_y_0_xxxyz_xyyyyyy[k];

                g_y_0_xxxxyz_yyyyyz[k] = -g_y_0_xxxyz_yyyyyz[k] * ab_x + g_y_0_xxxyz_xyyyyyz[k];

                g_y_0_xxxxyz_yyyyzz[k] = -g_y_0_xxxyz_yyyyzz[k] * ab_x + g_y_0_xxxyz_xyyyyzz[k];

                g_y_0_xxxxyz_yyyzzz[k] = -g_y_0_xxxyz_yyyzzz[k] * ab_x + g_y_0_xxxyz_xyyyzzz[k];

                g_y_0_xxxxyz_yyzzzz[k] = -g_y_0_xxxyz_yyzzzz[k] * ab_x + g_y_0_xxxyz_xyyzzzz[k];

                g_y_0_xxxxyz_yzzzzz[k] = -g_y_0_xxxyz_yzzzzz[k] * ab_x + g_y_0_xxxyz_xyzzzzz[k];

                g_y_0_xxxxyz_zzzzzz[k] = -g_y_0_xxxyz_zzzzzz[k] * ab_x + g_y_0_xxxyz_xzzzzzz[k];
            }

            /// Set up 924-952 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_xxxxxx = cbuffer.data(ii_geom_10_off + 924 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxxxy = cbuffer.data(ii_geom_10_off + 925 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxxxz = cbuffer.data(ii_geom_10_off + 926 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxxyy = cbuffer.data(ii_geom_10_off + 927 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxxyz = cbuffer.data(ii_geom_10_off + 928 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxxzz = cbuffer.data(ii_geom_10_off + 929 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxyyy = cbuffer.data(ii_geom_10_off + 930 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxyyz = cbuffer.data(ii_geom_10_off + 931 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxyzz = cbuffer.data(ii_geom_10_off + 932 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxzzz = cbuffer.data(ii_geom_10_off + 933 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxyyyy = cbuffer.data(ii_geom_10_off + 934 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxyyyz = cbuffer.data(ii_geom_10_off + 935 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxyyzz = cbuffer.data(ii_geom_10_off + 936 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxyzzz = cbuffer.data(ii_geom_10_off + 937 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxzzzz = cbuffer.data(ii_geom_10_off + 938 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyyyyy = cbuffer.data(ii_geom_10_off + 939 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyyyyz = cbuffer.data(ii_geom_10_off + 940 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyyyzz = cbuffer.data(ii_geom_10_off + 941 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyyzzz = cbuffer.data(ii_geom_10_off + 942 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyzzzz = cbuffer.data(ii_geom_10_off + 943 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xzzzzz = cbuffer.data(ii_geom_10_off + 944 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyyyyy = cbuffer.data(ii_geom_10_off + 945 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyyyyz = cbuffer.data(ii_geom_10_off + 946 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyyyzz = cbuffer.data(ii_geom_10_off + 947 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyyzzz = cbuffer.data(ii_geom_10_off + 948 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyzzzz = cbuffer.data(ii_geom_10_off + 949 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yzzzzz = cbuffer.data(ii_geom_10_off + 950 * ccomps * dcomps);

            auto g_y_0_xxxxzz_zzzzzz = cbuffer.data(ii_geom_10_off + 951 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxzz_xxxxxx, g_y_0_xxxxzz_xxxxxy, g_y_0_xxxxzz_xxxxxz, g_y_0_xxxxzz_xxxxyy, g_y_0_xxxxzz_xxxxyz, g_y_0_xxxxzz_xxxxzz, g_y_0_xxxxzz_xxxyyy, g_y_0_xxxxzz_xxxyyz, g_y_0_xxxxzz_xxxyzz, g_y_0_xxxxzz_xxxzzz, g_y_0_xxxxzz_xxyyyy, g_y_0_xxxxzz_xxyyyz, g_y_0_xxxxzz_xxyyzz, g_y_0_xxxxzz_xxyzzz, g_y_0_xxxxzz_xxzzzz, g_y_0_xxxxzz_xyyyyy, g_y_0_xxxxzz_xyyyyz, g_y_0_xxxxzz_xyyyzz, g_y_0_xxxxzz_xyyzzz, g_y_0_xxxxzz_xyzzzz, g_y_0_xxxxzz_xzzzzz, g_y_0_xxxxzz_yyyyyy, g_y_0_xxxxzz_yyyyyz, g_y_0_xxxxzz_yyyyzz, g_y_0_xxxxzz_yyyzzz, g_y_0_xxxxzz_yyzzzz, g_y_0_xxxxzz_yzzzzz, g_y_0_xxxxzz_zzzzzz, g_y_0_xxxzz_xxxxxx, g_y_0_xxxzz_xxxxxxx, g_y_0_xxxzz_xxxxxxy, g_y_0_xxxzz_xxxxxxz, g_y_0_xxxzz_xxxxxy, g_y_0_xxxzz_xxxxxyy, g_y_0_xxxzz_xxxxxyz, g_y_0_xxxzz_xxxxxz, g_y_0_xxxzz_xxxxxzz, g_y_0_xxxzz_xxxxyy, g_y_0_xxxzz_xxxxyyy, g_y_0_xxxzz_xxxxyyz, g_y_0_xxxzz_xxxxyz, g_y_0_xxxzz_xxxxyzz, g_y_0_xxxzz_xxxxzz, g_y_0_xxxzz_xxxxzzz, g_y_0_xxxzz_xxxyyy, g_y_0_xxxzz_xxxyyyy, g_y_0_xxxzz_xxxyyyz, g_y_0_xxxzz_xxxyyz, g_y_0_xxxzz_xxxyyzz, g_y_0_xxxzz_xxxyzz, g_y_0_xxxzz_xxxyzzz, g_y_0_xxxzz_xxxzzz, g_y_0_xxxzz_xxxzzzz, g_y_0_xxxzz_xxyyyy, g_y_0_xxxzz_xxyyyyy, g_y_0_xxxzz_xxyyyyz, g_y_0_xxxzz_xxyyyz, g_y_0_xxxzz_xxyyyzz, g_y_0_xxxzz_xxyyzz, g_y_0_xxxzz_xxyyzzz, g_y_0_xxxzz_xxyzzz, g_y_0_xxxzz_xxyzzzz, g_y_0_xxxzz_xxzzzz, g_y_0_xxxzz_xxzzzzz, g_y_0_xxxzz_xyyyyy, g_y_0_xxxzz_xyyyyyy, g_y_0_xxxzz_xyyyyyz, g_y_0_xxxzz_xyyyyz, g_y_0_xxxzz_xyyyyzz, g_y_0_xxxzz_xyyyzz, g_y_0_xxxzz_xyyyzzz, g_y_0_xxxzz_xyyzzz, g_y_0_xxxzz_xyyzzzz, g_y_0_xxxzz_xyzzzz, g_y_0_xxxzz_xyzzzzz, g_y_0_xxxzz_xzzzzz, g_y_0_xxxzz_xzzzzzz, g_y_0_xxxzz_yyyyyy, g_y_0_xxxzz_yyyyyz, g_y_0_xxxzz_yyyyzz, g_y_0_xxxzz_yyyzzz, g_y_0_xxxzz_yyzzzz, g_y_0_xxxzz_yzzzzz, g_y_0_xxxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_xxxxxx[k] = -g_y_0_xxxzz_xxxxxx[k] * ab_x + g_y_0_xxxzz_xxxxxxx[k];

                g_y_0_xxxxzz_xxxxxy[k] = -g_y_0_xxxzz_xxxxxy[k] * ab_x + g_y_0_xxxzz_xxxxxxy[k];

                g_y_0_xxxxzz_xxxxxz[k] = -g_y_0_xxxzz_xxxxxz[k] * ab_x + g_y_0_xxxzz_xxxxxxz[k];

                g_y_0_xxxxzz_xxxxyy[k] = -g_y_0_xxxzz_xxxxyy[k] * ab_x + g_y_0_xxxzz_xxxxxyy[k];

                g_y_0_xxxxzz_xxxxyz[k] = -g_y_0_xxxzz_xxxxyz[k] * ab_x + g_y_0_xxxzz_xxxxxyz[k];

                g_y_0_xxxxzz_xxxxzz[k] = -g_y_0_xxxzz_xxxxzz[k] * ab_x + g_y_0_xxxzz_xxxxxzz[k];

                g_y_0_xxxxzz_xxxyyy[k] = -g_y_0_xxxzz_xxxyyy[k] * ab_x + g_y_0_xxxzz_xxxxyyy[k];

                g_y_0_xxxxzz_xxxyyz[k] = -g_y_0_xxxzz_xxxyyz[k] * ab_x + g_y_0_xxxzz_xxxxyyz[k];

                g_y_0_xxxxzz_xxxyzz[k] = -g_y_0_xxxzz_xxxyzz[k] * ab_x + g_y_0_xxxzz_xxxxyzz[k];

                g_y_0_xxxxzz_xxxzzz[k] = -g_y_0_xxxzz_xxxzzz[k] * ab_x + g_y_0_xxxzz_xxxxzzz[k];

                g_y_0_xxxxzz_xxyyyy[k] = -g_y_0_xxxzz_xxyyyy[k] * ab_x + g_y_0_xxxzz_xxxyyyy[k];

                g_y_0_xxxxzz_xxyyyz[k] = -g_y_0_xxxzz_xxyyyz[k] * ab_x + g_y_0_xxxzz_xxxyyyz[k];

                g_y_0_xxxxzz_xxyyzz[k] = -g_y_0_xxxzz_xxyyzz[k] * ab_x + g_y_0_xxxzz_xxxyyzz[k];

                g_y_0_xxxxzz_xxyzzz[k] = -g_y_0_xxxzz_xxyzzz[k] * ab_x + g_y_0_xxxzz_xxxyzzz[k];

                g_y_0_xxxxzz_xxzzzz[k] = -g_y_0_xxxzz_xxzzzz[k] * ab_x + g_y_0_xxxzz_xxxzzzz[k];

                g_y_0_xxxxzz_xyyyyy[k] = -g_y_0_xxxzz_xyyyyy[k] * ab_x + g_y_0_xxxzz_xxyyyyy[k];

                g_y_0_xxxxzz_xyyyyz[k] = -g_y_0_xxxzz_xyyyyz[k] * ab_x + g_y_0_xxxzz_xxyyyyz[k];

                g_y_0_xxxxzz_xyyyzz[k] = -g_y_0_xxxzz_xyyyzz[k] * ab_x + g_y_0_xxxzz_xxyyyzz[k];

                g_y_0_xxxxzz_xyyzzz[k] = -g_y_0_xxxzz_xyyzzz[k] * ab_x + g_y_0_xxxzz_xxyyzzz[k];

                g_y_0_xxxxzz_xyzzzz[k] = -g_y_0_xxxzz_xyzzzz[k] * ab_x + g_y_0_xxxzz_xxyzzzz[k];

                g_y_0_xxxxzz_xzzzzz[k] = -g_y_0_xxxzz_xzzzzz[k] * ab_x + g_y_0_xxxzz_xxzzzzz[k];

                g_y_0_xxxxzz_yyyyyy[k] = -g_y_0_xxxzz_yyyyyy[k] * ab_x + g_y_0_xxxzz_xyyyyyy[k];

                g_y_0_xxxxzz_yyyyyz[k] = -g_y_0_xxxzz_yyyyyz[k] * ab_x + g_y_0_xxxzz_xyyyyyz[k];

                g_y_0_xxxxzz_yyyyzz[k] = -g_y_0_xxxzz_yyyyzz[k] * ab_x + g_y_0_xxxzz_xyyyyzz[k];

                g_y_0_xxxxzz_yyyzzz[k] = -g_y_0_xxxzz_yyyzzz[k] * ab_x + g_y_0_xxxzz_xyyyzzz[k];

                g_y_0_xxxxzz_yyzzzz[k] = -g_y_0_xxxzz_yyzzzz[k] * ab_x + g_y_0_xxxzz_xyyzzzz[k];

                g_y_0_xxxxzz_yzzzzz[k] = -g_y_0_xxxzz_yzzzzz[k] * ab_x + g_y_0_xxxzz_xyzzzzz[k];

                g_y_0_xxxxzz_zzzzzz[k] = -g_y_0_xxxzz_zzzzzz[k] * ab_x + g_y_0_xxxzz_xzzzzzz[k];
            }

            /// Set up 952-980 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 952 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 953 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 954 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 955 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 956 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 957 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 958 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 959 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 960 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 961 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 962 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 963 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 964 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 965 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 966 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 967 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 968 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 969 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 970 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 971 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 972 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 973 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 974 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 975 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 976 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 977 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 978 * ccomps * dcomps);

            auto g_y_0_xxxyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 979 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyyy_xxxxxx, g_y_0_xxxyyy_xxxxxy, g_y_0_xxxyyy_xxxxxz, g_y_0_xxxyyy_xxxxyy, g_y_0_xxxyyy_xxxxyz, g_y_0_xxxyyy_xxxxzz, g_y_0_xxxyyy_xxxyyy, g_y_0_xxxyyy_xxxyyz, g_y_0_xxxyyy_xxxyzz, g_y_0_xxxyyy_xxxzzz, g_y_0_xxxyyy_xxyyyy, g_y_0_xxxyyy_xxyyyz, g_y_0_xxxyyy_xxyyzz, g_y_0_xxxyyy_xxyzzz, g_y_0_xxxyyy_xxzzzz, g_y_0_xxxyyy_xyyyyy, g_y_0_xxxyyy_xyyyyz, g_y_0_xxxyyy_xyyyzz, g_y_0_xxxyyy_xyyzzz, g_y_0_xxxyyy_xyzzzz, g_y_0_xxxyyy_xzzzzz, g_y_0_xxxyyy_yyyyyy, g_y_0_xxxyyy_yyyyyz, g_y_0_xxxyyy_yyyyzz, g_y_0_xxxyyy_yyyzzz, g_y_0_xxxyyy_yyzzzz, g_y_0_xxxyyy_yzzzzz, g_y_0_xxxyyy_zzzzzz, g_y_0_xxyyy_xxxxxx, g_y_0_xxyyy_xxxxxxx, g_y_0_xxyyy_xxxxxxy, g_y_0_xxyyy_xxxxxxz, g_y_0_xxyyy_xxxxxy, g_y_0_xxyyy_xxxxxyy, g_y_0_xxyyy_xxxxxyz, g_y_0_xxyyy_xxxxxz, g_y_0_xxyyy_xxxxxzz, g_y_0_xxyyy_xxxxyy, g_y_0_xxyyy_xxxxyyy, g_y_0_xxyyy_xxxxyyz, g_y_0_xxyyy_xxxxyz, g_y_0_xxyyy_xxxxyzz, g_y_0_xxyyy_xxxxzz, g_y_0_xxyyy_xxxxzzz, g_y_0_xxyyy_xxxyyy, g_y_0_xxyyy_xxxyyyy, g_y_0_xxyyy_xxxyyyz, g_y_0_xxyyy_xxxyyz, g_y_0_xxyyy_xxxyyzz, g_y_0_xxyyy_xxxyzz, g_y_0_xxyyy_xxxyzzz, g_y_0_xxyyy_xxxzzz, g_y_0_xxyyy_xxxzzzz, g_y_0_xxyyy_xxyyyy, g_y_0_xxyyy_xxyyyyy, g_y_0_xxyyy_xxyyyyz, g_y_0_xxyyy_xxyyyz, g_y_0_xxyyy_xxyyyzz, g_y_0_xxyyy_xxyyzz, g_y_0_xxyyy_xxyyzzz, g_y_0_xxyyy_xxyzzz, g_y_0_xxyyy_xxyzzzz, g_y_0_xxyyy_xxzzzz, g_y_0_xxyyy_xxzzzzz, g_y_0_xxyyy_xyyyyy, g_y_0_xxyyy_xyyyyyy, g_y_0_xxyyy_xyyyyyz, g_y_0_xxyyy_xyyyyz, g_y_0_xxyyy_xyyyyzz, g_y_0_xxyyy_xyyyzz, g_y_0_xxyyy_xyyyzzz, g_y_0_xxyyy_xyyzzz, g_y_0_xxyyy_xyyzzzz, g_y_0_xxyyy_xyzzzz, g_y_0_xxyyy_xyzzzzz, g_y_0_xxyyy_xzzzzz, g_y_0_xxyyy_xzzzzzz, g_y_0_xxyyy_yyyyyy, g_y_0_xxyyy_yyyyyz, g_y_0_xxyyy_yyyyzz, g_y_0_xxyyy_yyyzzz, g_y_0_xxyyy_yyzzzz, g_y_0_xxyyy_yzzzzz, g_y_0_xxyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_xxxxxx[k] = -g_y_0_xxyyy_xxxxxx[k] * ab_x + g_y_0_xxyyy_xxxxxxx[k];

                g_y_0_xxxyyy_xxxxxy[k] = -g_y_0_xxyyy_xxxxxy[k] * ab_x + g_y_0_xxyyy_xxxxxxy[k];

                g_y_0_xxxyyy_xxxxxz[k] = -g_y_0_xxyyy_xxxxxz[k] * ab_x + g_y_0_xxyyy_xxxxxxz[k];

                g_y_0_xxxyyy_xxxxyy[k] = -g_y_0_xxyyy_xxxxyy[k] * ab_x + g_y_0_xxyyy_xxxxxyy[k];

                g_y_0_xxxyyy_xxxxyz[k] = -g_y_0_xxyyy_xxxxyz[k] * ab_x + g_y_0_xxyyy_xxxxxyz[k];

                g_y_0_xxxyyy_xxxxzz[k] = -g_y_0_xxyyy_xxxxzz[k] * ab_x + g_y_0_xxyyy_xxxxxzz[k];

                g_y_0_xxxyyy_xxxyyy[k] = -g_y_0_xxyyy_xxxyyy[k] * ab_x + g_y_0_xxyyy_xxxxyyy[k];

                g_y_0_xxxyyy_xxxyyz[k] = -g_y_0_xxyyy_xxxyyz[k] * ab_x + g_y_0_xxyyy_xxxxyyz[k];

                g_y_0_xxxyyy_xxxyzz[k] = -g_y_0_xxyyy_xxxyzz[k] * ab_x + g_y_0_xxyyy_xxxxyzz[k];

                g_y_0_xxxyyy_xxxzzz[k] = -g_y_0_xxyyy_xxxzzz[k] * ab_x + g_y_0_xxyyy_xxxxzzz[k];

                g_y_0_xxxyyy_xxyyyy[k] = -g_y_0_xxyyy_xxyyyy[k] * ab_x + g_y_0_xxyyy_xxxyyyy[k];

                g_y_0_xxxyyy_xxyyyz[k] = -g_y_0_xxyyy_xxyyyz[k] * ab_x + g_y_0_xxyyy_xxxyyyz[k];

                g_y_0_xxxyyy_xxyyzz[k] = -g_y_0_xxyyy_xxyyzz[k] * ab_x + g_y_0_xxyyy_xxxyyzz[k];

                g_y_0_xxxyyy_xxyzzz[k] = -g_y_0_xxyyy_xxyzzz[k] * ab_x + g_y_0_xxyyy_xxxyzzz[k];

                g_y_0_xxxyyy_xxzzzz[k] = -g_y_0_xxyyy_xxzzzz[k] * ab_x + g_y_0_xxyyy_xxxzzzz[k];

                g_y_0_xxxyyy_xyyyyy[k] = -g_y_0_xxyyy_xyyyyy[k] * ab_x + g_y_0_xxyyy_xxyyyyy[k];

                g_y_0_xxxyyy_xyyyyz[k] = -g_y_0_xxyyy_xyyyyz[k] * ab_x + g_y_0_xxyyy_xxyyyyz[k];

                g_y_0_xxxyyy_xyyyzz[k] = -g_y_0_xxyyy_xyyyzz[k] * ab_x + g_y_0_xxyyy_xxyyyzz[k];

                g_y_0_xxxyyy_xyyzzz[k] = -g_y_0_xxyyy_xyyzzz[k] * ab_x + g_y_0_xxyyy_xxyyzzz[k];

                g_y_0_xxxyyy_xyzzzz[k] = -g_y_0_xxyyy_xyzzzz[k] * ab_x + g_y_0_xxyyy_xxyzzzz[k];

                g_y_0_xxxyyy_xzzzzz[k] = -g_y_0_xxyyy_xzzzzz[k] * ab_x + g_y_0_xxyyy_xxzzzzz[k];

                g_y_0_xxxyyy_yyyyyy[k] = -g_y_0_xxyyy_yyyyyy[k] * ab_x + g_y_0_xxyyy_xyyyyyy[k];

                g_y_0_xxxyyy_yyyyyz[k] = -g_y_0_xxyyy_yyyyyz[k] * ab_x + g_y_0_xxyyy_xyyyyyz[k];

                g_y_0_xxxyyy_yyyyzz[k] = -g_y_0_xxyyy_yyyyzz[k] * ab_x + g_y_0_xxyyy_xyyyyzz[k];

                g_y_0_xxxyyy_yyyzzz[k] = -g_y_0_xxyyy_yyyzzz[k] * ab_x + g_y_0_xxyyy_xyyyzzz[k];

                g_y_0_xxxyyy_yyzzzz[k] = -g_y_0_xxyyy_yyzzzz[k] * ab_x + g_y_0_xxyyy_xyyzzzz[k];

                g_y_0_xxxyyy_yzzzzz[k] = -g_y_0_xxyyy_yzzzzz[k] * ab_x + g_y_0_xxyyy_xyzzzzz[k];

                g_y_0_xxxyyy_zzzzzz[k] = -g_y_0_xxyyy_zzzzzz[k] * ab_x + g_y_0_xxyyy_xzzzzzz[k];
            }

            /// Set up 980-1008 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 980 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 981 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 982 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 983 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 984 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 985 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 986 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 987 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 988 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 989 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 990 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 991 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 992 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 993 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 994 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 995 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 996 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 997 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 998 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 999 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 1000 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 1001 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 1002 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 1003 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 1004 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 1005 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 1006 * ccomps * dcomps);

            auto g_y_0_xxxyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyyz_xxxxxx, g_y_0_xxxyyz_xxxxxy, g_y_0_xxxyyz_xxxxxz, g_y_0_xxxyyz_xxxxyy, g_y_0_xxxyyz_xxxxyz, g_y_0_xxxyyz_xxxxzz, g_y_0_xxxyyz_xxxyyy, g_y_0_xxxyyz_xxxyyz, g_y_0_xxxyyz_xxxyzz, g_y_0_xxxyyz_xxxzzz, g_y_0_xxxyyz_xxyyyy, g_y_0_xxxyyz_xxyyyz, g_y_0_xxxyyz_xxyyzz, g_y_0_xxxyyz_xxyzzz, g_y_0_xxxyyz_xxzzzz, g_y_0_xxxyyz_xyyyyy, g_y_0_xxxyyz_xyyyyz, g_y_0_xxxyyz_xyyyzz, g_y_0_xxxyyz_xyyzzz, g_y_0_xxxyyz_xyzzzz, g_y_0_xxxyyz_xzzzzz, g_y_0_xxxyyz_yyyyyy, g_y_0_xxxyyz_yyyyyz, g_y_0_xxxyyz_yyyyzz, g_y_0_xxxyyz_yyyzzz, g_y_0_xxxyyz_yyzzzz, g_y_0_xxxyyz_yzzzzz, g_y_0_xxxyyz_zzzzzz, g_y_0_xxyyz_xxxxxx, g_y_0_xxyyz_xxxxxxx, g_y_0_xxyyz_xxxxxxy, g_y_0_xxyyz_xxxxxxz, g_y_0_xxyyz_xxxxxy, g_y_0_xxyyz_xxxxxyy, g_y_0_xxyyz_xxxxxyz, g_y_0_xxyyz_xxxxxz, g_y_0_xxyyz_xxxxxzz, g_y_0_xxyyz_xxxxyy, g_y_0_xxyyz_xxxxyyy, g_y_0_xxyyz_xxxxyyz, g_y_0_xxyyz_xxxxyz, g_y_0_xxyyz_xxxxyzz, g_y_0_xxyyz_xxxxzz, g_y_0_xxyyz_xxxxzzz, g_y_0_xxyyz_xxxyyy, g_y_0_xxyyz_xxxyyyy, g_y_0_xxyyz_xxxyyyz, g_y_0_xxyyz_xxxyyz, g_y_0_xxyyz_xxxyyzz, g_y_0_xxyyz_xxxyzz, g_y_0_xxyyz_xxxyzzz, g_y_0_xxyyz_xxxzzz, g_y_0_xxyyz_xxxzzzz, g_y_0_xxyyz_xxyyyy, g_y_0_xxyyz_xxyyyyy, g_y_0_xxyyz_xxyyyyz, g_y_0_xxyyz_xxyyyz, g_y_0_xxyyz_xxyyyzz, g_y_0_xxyyz_xxyyzz, g_y_0_xxyyz_xxyyzzz, g_y_0_xxyyz_xxyzzz, g_y_0_xxyyz_xxyzzzz, g_y_0_xxyyz_xxzzzz, g_y_0_xxyyz_xxzzzzz, g_y_0_xxyyz_xyyyyy, g_y_0_xxyyz_xyyyyyy, g_y_0_xxyyz_xyyyyyz, g_y_0_xxyyz_xyyyyz, g_y_0_xxyyz_xyyyyzz, g_y_0_xxyyz_xyyyzz, g_y_0_xxyyz_xyyyzzz, g_y_0_xxyyz_xyyzzz, g_y_0_xxyyz_xyyzzzz, g_y_0_xxyyz_xyzzzz, g_y_0_xxyyz_xyzzzzz, g_y_0_xxyyz_xzzzzz, g_y_0_xxyyz_xzzzzzz, g_y_0_xxyyz_yyyyyy, g_y_0_xxyyz_yyyyyz, g_y_0_xxyyz_yyyyzz, g_y_0_xxyyz_yyyzzz, g_y_0_xxyyz_yyzzzz, g_y_0_xxyyz_yzzzzz, g_y_0_xxyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_xxxxxx[k] = -g_y_0_xxyyz_xxxxxx[k] * ab_x + g_y_0_xxyyz_xxxxxxx[k];

                g_y_0_xxxyyz_xxxxxy[k] = -g_y_0_xxyyz_xxxxxy[k] * ab_x + g_y_0_xxyyz_xxxxxxy[k];

                g_y_0_xxxyyz_xxxxxz[k] = -g_y_0_xxyyz_xxxxxz[k] * ab_x + g_y_0_xxyyz_xxxxxxz[k];

                g_y_0_xxxyyz_xxxxyy[k] = -g_y_0_xxyyz_xxxxyy[k] * ab_x + g_y_0_xxyyz_xxxxxyy[k];

                g_y_0_xxxyyz_xxxxyz[k] = -g_y_0_xxyyz_xxxxyz[k] * ab_x + g_y_0_xxyyz_xxxxxyz[k];

                g_y_0_xxxyyz_xxxxzz[k] = -g_y_0_xxyyz_xxxxzz[k] * ab_x + g_y_0_xxyyz_xxxxxzz[k];

                g_y_0_xxxyyz_xxxyyy[k] = -g_y_0_xxyyz_xxxyyy[k] * ab_x + g_y_0_xxyyz_xxxxyyy[k];

                g_y_0_xxxyyz_xxxyyz[k] = -g_y_0_xxyyz_xxxyyz[k] * ab_x + g_y_0_xxyyz_xxxxyyz[k];

                g_y_0_xxxyyz_xxxyzz[k] = -g_y_0_xxyyz_xxxyzz[k] * ab_x + g_y_0_xxyyz_xxxxyzz[k];

                g_y_0_xxxyyz_xxxzzz[k] = -g_y_0_xxyyz_xxxzzz[k] * ab_x + g_y_0_xxyyz_xxxxzzz[k];

                g_y_0_xxxyyz_xxyyyy[k] = -g_y_0_xxyyz_xxyyyy[k] * ab_x + g_y_0_xxyyz_xxxyyyy[k];

                g_y_0_xxxyyz_xxyyyz[k] = -g_y_0_xxyyz_xxyyyz[k] * ab_x + g_y_0_xxyyz_xxxyyyz[k];

                g_y_0_xxxyyz_xxyyzz[k] = -g_y_0_xxyyz_xxyyzz[k] * ab_x + g_y_0_xxyyz_xxxyyzz[k];

                g_y_0_xxxyyz_xxyzzz[k] = -g_y_0_xxyyz_xxyzzz[k] * ab_x + g_y_0_xxyyz_xxxyzzz[k];

                g_y_0_xxxyyz_xxzzzz[k] = -g_y_0_xxyyz_xxzzzz[k] * ab_x + g_y_0_xxyyz_xxxzzzz[k];

                g_y_0_xxxyyz_xyyyyy[k] = -g_y_0_xxyyz_xyyyyy[k] * ab_x + g_y_0_xxyyz_xxyyyyy[k];

                g_y_0_xxxyyz_xyyyyz[k] = -g_y_0_xxyyz_xyyyyz[k] * ab_x + g_y_0_xxyyz_xxyyyyz[k];

                g_y_0_xxxyyz_xyyyzz[k] = -g_y_0_xxyyz_xyyyzz[k] * ab_x + g_y_0_xxyyz_xxyyyzz[k];

                g_y_0_xxxyyz_xyyzzz[k] = -g_y_0_xxyyz_xyyzzz[k] * ab_x + g_y_0_xxyyz_xxyyzzz[k];

                g_y_0_xxxyyz_xyzzzz[k] = -g_y_0_xxyyz_xyzzzz[k] * ab_x + g_y_0_xxyyz_xxyzzzz[k];

                g_y_0_xxxyyz_xzzzzz[k] = -g_y_0_xxyyz_xzzzzz[k] * ab_x + g_y_0_xxyyz_xxzzzzz[k];

                g_y_0_xxxyyz_yyyyyy[k] = -g_y_0_xxyyz_yyyyyy[k] * ab_x + g_y_0_xxyyz_xyyyyyy[k];

                g_y_0_xxxyyz_yyyyyz[k] = -g_y_0_xxyyz_yyyyyz[k] * ab_x + g_y_0_xxyyz_xyyyyyz[k];

                g_y_0_xxxyyz_yyyyzz[k] = -g_y_0_xxyyz_yyyyzz[k] * ab_x + g_y_0_xxyyz_xyyyyzz[k];

                g_y_0_xxxyyz_yyyzzz[k] = -g_y_0_xxyyz_yyyzzz[k] * ab_x + g_y_0_xxyyz_xyyyzzz[k];

                g_y_0_xxxyyz_yyzzzz[k] = -g_y_0_xxyyz_yyzzzz[k] * ab_x + g_y_0_xxyyz_xyyzzzz[k];

                g_y_0_xxxyyz_yzzzzz[k] = -g_y_0_xxyyz_yzzzzz[k] * ab_x + g_y_0_xxyyz_xyzzzzz[k];

                g_y_0_xxxyyz_zzzzzz[k] = -g_y_0_xxyyz_zzzzzz[k] * ab_x + g_y_0_xxyyz_xzzzzzz[k];
            }

            /// Set up 1008-1036 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1008 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1009 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1010 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1011 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1012 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1013 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1014 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1015 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1016 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1017 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1018 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1019 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1020 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1021 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1022 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1023 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1024 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1025 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1026 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1027 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1028 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1029 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1030 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1031 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1032 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1033 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1034 * ccomps * dcomps);

            auto g_y_0_xxxyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1035 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyzz_xxxxxx, g_y_0_xxxyzz_xxxxxy, g_y_0_xxxyzz_xxxxxz, g_y_0_xxxyzz_xxxxyy, g_y_0_xxxyzz_xxxxyz, g_y_0_xxxyzz_xxxxzz, g_y_0_xxxyzz_xxxyyy, g_y_0_xxxyzz_xxxyyz, g_y_0_xxxyzz_xxxyzz, g_y_0_xxxyzz_xxxzzz, g_y_0_xxxyzz_xxyyyy, g_y_0_xxxyzz_xxyyyz, g_y_0_xxxyzz_xxyyzz, g_y_0_xxxyzz_xxyzzz, g_y_0_xxxyzz_xxzzzz, g_y_0_xxxyzz_xyyyyy, g_y_0_xxxyzz_xyyyyz, g_y_0_xxxyzz_xyyyzz, g_y_0_xxxyzz_xyyzzz, g_y_0_xxxyzz_xyzzzz, g_y_0_xxxyzz_xzzzzz, g_y_0_xxxyzz_yyyyyy, g_y_0_xxxyzz_yyyyyz, g_y_0_xxxyzz_yyyyzz, g_y_0_xxxyzz_yyyzzz, g_y_0_xxxyzz_yyzzzz, g_y_0_xxxyzz_yzzzzz, g_y_0_xxxyzz_zzzzzz, g_y_0_xxyzz_xxxxxx, g_y_0_xxyzz_xxxxxxx, g_y_0_xxyzz_xxxxxxy, g_y_0_xxyzz_xxxxxxz, g_y_0_xxyzz_xxxxxy, g_y_0_xxyzz_xxxxxyy, g_y_0_xxyzz_xxxxxyz, g_y_0_xxyzz_xxxxxz, g_y_0_xxyzz_xxxxxzz, g_y_0_xxyzz_xxxxyy, g_y_0_xxyzz_xxxxyyy, g_y_0_xxyzz_xxxxyyz, g_y_0_xxyzz_xxxxyz, g_y_0_xxyzz_xxxxyzz, g_y_0_xxyzz_xxxxzz, g_y_0_xxyzz_xxxxzzz, g_y_0_xxyzz_xxxyyy, g_y_0_xxyzz_xxxyyyy, g_y_0_xxyzz_xxxyyyz, g_y_0_xxyzz_xxxyyz, g_y_0_xxyzz_xxxyyzz, g_y_0_xxyzz_xxxyzz, g_y_0_xxyzz_xxxyzzz, g_y_0_xxyzz_xxxzzz, g_y_0_xxyzz_xxxzzzz, g_y_0_xxyzz_xxyyyy, g_y_0_xxyzz_xxyyyyy, g_y_0_xxyzz_xxyyyyz, g_y_0_xxyzz_xxyyyz, g_y_0_xxyzz_xxyyyzz, g_y_0_xxyzz_xxyyzz, g_y_0_xxyzz_xxyyzzz, g_y_0_xxyzz_xxyzzz, g_y_0_xxyzz_xxyzzzz, g_y_0_xxyzz_xxzzzz, g_y_0_xxyzz_xxzzzzz, g_y_0_xxyzz_xyyyyy, g_y_0_xxyzz_xyyyyyy, g_y_0_xxyzz_xyyyyyz, g_y_0_xxyzz_xyyyyz, g_y_0_xxyzz_xyyyyzz, g_y_0_xxyzz_xyyyzz, g_y_0_xxyzz_xyyyzzz, g_y_0_xxyzz_xyyzzz, g_y_0_xxyzz_xyyzzzz, g_y_0_xxyzz_xyzzzz, g_y_0_xxyzz_xyzzzzz, g_y_0_xxyzz_xzzzzz, g_y_0_xxyzz_xzzzzzz, g_y_0_xxyzz_yyyyyy, g_y_0_xxyzz_yyyyyz, g_y_0_xxyzz_yyyyzz, g_y_0_xxyzz_yyyzzz, g_y_0_xxyzz_yyzzzz, g_y_0_xxyzz_yzzzzz, g_y_0_xxyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_xxxxxx[k] = -g_y_0_xxyzz_xxxxxx[k] * ab_x + g_y_0_xxyzz_xxxxxxx[k];

                g_y_0_xxxyzz_xxxxxy[k] = -g_y_0_xxyzz_xxxxxy[k] * ab_x + g_y_0_xxyzz_xxxxxxy[k];

                g_y_0_xxxyzz_xxxxxz[k] = -g_y_0_xxyzz_xxxxxz[k] * ab_x + g_y_0_xxyzz_xxxxxxz[k];

                g_y_0_xxxyzz_xxxxyy[k] = -g_y_0_xxyzz_xxxxyy[k] * ab_x + g_y_0_xxyzz_xxxxxyy[k];

                g_y_0_xxxyzz_xxxxyz[k] = -g_y_0_xxyzz_xxxxyz[k] * ab_x + g_y_0_xxyzz_xxxxxyz[k];

                g_y_0_xxxyzz_xxxxzz[k] = -g_y_0_xxyzz_xxxxzz[k] * ab_x + g_y_0_xxyzz_xxxxxzz[k];

                g_y_0_xxxyzz_xxxyyy[k] = -g_y_0_xxyzz_xxxyyy[k] * ab_x + g_y_0_xxyzz_xxxxyyy[k];

                g_y_0_xxxyzz_xxxyyz[k] = -g_y_0_xxyzz_xxxyyz[k] * ab_x + g_y_0_xxyzz_xxxxyyz[k];

                g_y_0_xxxyzz_xxxyzz[k] = -g_y_0_xxyzz_xxxyzz[k] * ab_x + g_y_0_xxyzz_xxxxyzz[k];

                g_y_0_xxxyzz_xxxzzz[k] = -g_y_0_xxyzz_xxxzzz[k] * ab_x + g_y_0_xxyzz_xxxxzzz[k];

                g_y_0_xxxyzz_xxyyyy[k] = -g_y_0_xxyzz_xxyyyy[k] * ab_x + g_y_0_xxyzz_xxxyyyy[k];

                g_y_0_xxxyzz_xxyyyz[k] = -g_y_0_xxyzz_xxyyyz[k] * ab_x + g_y_0_xxyzz_xxxyyyz[k];

                g_y_0_xxxyzz_xxyyzz[k] = -g_y_0_xxyzz_xxyyzz[k] * ab_x + g_y_0_xxyzz_xxxyyzz[k];

                g_y_0_xxxyzz_xxyzzz[k] = -g_y_0_xxyzz_xxyzzz[k] * ab_x + g_y_0_xxyzz_xxxyzzz[k];

                g_y_0_xxxyzz_xxzzzz[k] = -g_y_0_xxyzz_xxzzzz[k] * ab_x + g_y_0_xxyzz_xxxzzzz[k];

                g_y_0_xxxyzz_xyyyyy[k] = -g_y_0_xxyzz_xyyyyy[k] * ab_x + g_y_0_xxyzz_xxyyyyy[k];

                g_y_0_xxxyzz_xyyyyz[k] = -g_y_0_xxyzz_xyyyyz[k] * ab_x + g_y_0_xxyzz_xxyyyyz[k];

                g_y_0_xxxyzz_xyyyzz[k] = -g_y_0_xxyzz_xyyyzz[k] * ab_x + g_y_0_xxyzz_xxyyyzz[k];

                g_y_0_xxxyzz_xyyzzz[k] = -g_y_0_xxyzz_xyyzzz[k] * ab_x + g_y_0_xxyzz_xxyyzzz[k];

                g_y_0_xxxyzz_xyzzzz[k] = -g_y_0_xxyzz_xyzzzz[k] * ab_x + g_y_0_xxyzz_xxyzzzz[k];

                g_y_0_xxxyzz_xzzzzz[k] = -g_y_0_xxyzz_xzzzzz[k] * ab_x + g_y_0_xxyzz_xxzzzzz[k];

                g_y_0_xxxyzz_yyyyyy[k] = -g_y_0_xxyzz_yyyyyy[k] * ab_x + g_y_0_xxyzz_xyyyyyy[k];

                g_y_0_xxxyzz_yyyyyz[k] = -g_y_0_xxyzz_yyyyyz[k] * ab_x + g_y_0_xxyzz_xyyyyyz[k];

                g_y_0_xxxyzz_yyyyzz[k] = -g_y_0_xxyzz_yyyyzz[k] * ab_x + g_y_0_xxyzz_xyyyyzz[k];

                g_y_0_xxxyzz_yyyzzz[k] = -g_y_0_xxyzz_yyyzzz[k] * ab_x + g_y_0_xxyzz_xyyyzzz[k];

                g_y_0_xxxyzz_yyzzzz[k] = -g_y_0_xxyzz_yyzzzz[k] * ab_x + g_y_0_xxyzz_xyyzzzz[k];

                g_y_0_xxxyzz_yzzzzz[k] = -g_y_0_xxyzz_yzzzzz[k] * ab_x + g_y_0_xxyzz_xyzzzzz[k];

                g_y_0_xxxyzz_zzzzzz[k] = -g_y_0_xxyzz_zzzzzz[k] * ab_x + g_y_0_xxyzz_xzzzzzz[k];
            }

            /// Set up 1036-1064 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1036 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1037 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1038 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1039 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1040 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1041 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1042 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1043 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1044 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1045 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1046 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1047 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1048 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1049 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1050 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1051 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1052 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1053 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1054 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1055 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1056 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1057 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1058 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1059 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1060 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1061 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1062 * ccomps * dcomps);

            auto g_y_0_xxxzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1063 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxzzz_xxxxxx, g_y_0_xxxzzz_xxxxxy, g_y_0_xxxzzz_xxxxxz, g_y_0_xxxzzz_xxxxyy, g_y_0_xxxzzz_xxxxyz, g_y_0_xxxzzz_xxxxzz, g_y_0_xxxzzz_xxxyyy, g_y_0_xxxzzz_xxxyyz, g_y_0_xxxzzz_xxxyzz, g_y_0_xxxzzz_xxxzzz, g_y_0_xxxzzz_xxyyyy, g_y_0_xxxzzz_xxyyyz, g_y_0_xxxzzz_xxyyzz, g_y_0_xxxzzz_xxyzzz, g_y_0_xxxzzz_xxzzzz, g_y_0_xxxzzz_xyyyyy, g_y_0_xxxzzz_xyyyyz, g_y_0_xxxzzz_xyyyzz, g_y_0_xxxzzz_xyyzzz, g_y_0_xxxzzz_xyzzzz, g_y_0_xxxzzz_xzzzzz, g_y_0_xxxzzz_yyyyyy, g_y_0_xxxzzz_yyyyyz, g_y_0_xxxzzz_yyyyzz, g_y_0_xxxzzz_yyyzzz, g_y_0_xxxzzz_yyzzzz, g_y_0_xxxzzz_yzzzzz, g_y_0_xxxzzz_zzzzzz, g_y_0_xxzzz_xxxxxx, g_y_0_xxzzz_xxxxxxx, g_y_0_xxzzz_xxxxxxy, g_y_0_xxzzz_xxxxxxz, g_y_0_xxzzz_xxxxxy, g_y_0_xxzzz_xxxxxyy, g_y_0_xxzzz_xxxxxyz, g_y_0_xxzzz_xxxxxz, g_y_0_xxzzz_xxxxxzz, g_y_0_xxzzz_xxxxyy, g_y_0_xxzzz_xxxxyyy, g_y_0_xxzzz_xxxxyyz, g_y_0_xxzzz_xxxxyz, g_y_0_xxzzz_xxxxyzz, g_y_0_xxzzz_xxxxzz, g_y_0_xxzzz_xxxxzzz, g_y_0_xxzzz_xxxyyy, g_y_0_xxzzz_xxxyyyy, g_y_0_xxzzz_xxxyyyz, g_y_0_xxzzz_xxxyyz, g_y_0_xxzzz_xxxyyzz, g_y_0_xxzzz_xxxyzz, g_y_0_xxzzz_xxxyzzz, g_y_0_xxzzz_xxxzzz, g_y_0_xxzzz_xxxzzzz, g_y_0_xxzzz_xxyyyy, g_y_0_xxzzz_xxyyyyy, g_y_0_xxzzz_xxyyyyz, g_y_0_xxzzz_xxyyyz, g_y_0_xxzzz_xxyyyzz, g_y_0_xxzzz_xxyyzz, g_y_0_xxzzz_xxyyzzz, g_y_0_xxzzz_xxyzzz, g_y_0_xxzzz_xxyzzzz, g_y_0_xxzzz_xxzzzz, g_y_0_xxzzz_xxzzzzz, g_y_0_xxzzz_xyyyyy, g_y_0_xxzzz_xyyyyyy, g_y_0_xxzzz_xyyyyyz, g_y_0_xxzzz_xyyyyz, g_y_0_xxzzz_xyyyyzz, g_y_0_xxzzz_xyyyzz, g_y_0_xxzzz_xyyyzzz, g_y_0_xxzzz_xyyzzz, g_y_0_xxzzz_xyyzzzz, g_y_0_xxzzz_xyzzzz, g_y_0_xxzzz_xyzzzzz, g_y_0_xxzzz_xzzzzz, g_y_0_xxzzz_xzzzzzz, g_y_0_xxzzz_yyyyyy, g_y_0_xxzzz_yyyyyz, g_y_0_xxzzz_yyyyzz, g_y_0_xxzzz_yyyzzz, g_y_0_xxzzz_yyzzzz, g_y_0_xxzzz_yzzzzz, g_y_0_xxzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_xxxxxx[k] = -g_y_0_xxzzz_xxxxxx[k] * ab_x + g_y_0_xxzzz_xxxxxxx[k];

                g_y_0_xxxzzz_xxxxxy[k] = -g_y_0_xxzzz_xxxxxy[k] * ab_x + g_y_0_xxzzz_xxxxxxy[k];

                g_y_0_xxxzzz_xxxxxz[k] = -g_y_0_xxzzz_xxxxxz[k] * ab_x + g_y_0_xxzzz_xxxxxxz[k];

                g_y_0_xxxzzz_xxxxyy[k] = -g_y_0_xxzzz_xxxxyy[k] * ab_x + g_y_0_xxzzz_xxxxxyy[k];

                g_y_0_xxxzzz_xxxxyz[k] = -g_y_0_xxzzz_xxxxyz[k] * ab_x + g_y_0_xxzzz_xxxxxyz[k];

                g_y_0_xxxzzz_xxxxzz[k] = -g_y_0_xxzzz_xxxxzz[k] * ab_x + g_y_0_xxzzz_xxxxxzz[k];

                g_y_0_xxxzzz_xxxyyy[k] = -g_y_0_xxzzz_xxxyyy[k] * ab_x + g_y_0_xxzzz_xxxxyyy[k];

                g_y_0_xxxzzz_xxxyyz[k] = -g_y_0_xxzzz_xxxyyz[k] * ab_x + g_y_0_xxzzz_xxxxyyz[k];

                g_y_0_xxxzzz_xxxyzz[k] = -g_y_0_xxzzz_xxxyzz[k] * ab_x + g_y_0_xxzzz_xxxxyzz[k];

                g_y_0_xxxzzz_xxxzzz[k] = -g_y_0_xxzzz_xxxzzz[k] * ab_x + g_y_0_xxzzz_xxxxzzz[k];

                g_y_0_xxxzzz_xxyyyy[k] = -g_y_0_xxzzz_xxyyyy[k] * ab_x + g_y_0_xxzzz_xxxyyyy[k];

                g_y_0_xxxzzz_xxyyyz[k] = -g_y_0_xxzzz_xxyyyz[k] * ab_x + g_y_0_xxzzz_xxxyyyz[k];

                g_y_0_xxxzzz_xxyyzz[k] = -g_y_0_xxzzz_xxyyzz[k] * ab_x + g_y_0_xxzzz_xxxyyzz[k];

                g_y_0_xxxzzz_xxyzzz[k] = -g_y_0_xxzzz_xxyzzz[k] * ab_x + g_y_0_xxzzz_xxxyzzz[k];

                g_y_0_xxxzzz_xxzzzz[k] = -g_y_0_xxzzz_xxzzzz[k] * ab_x + g_y_0_xxzzz_xxxzzzz[k];

                g_y_0_xxxzzz_xyyyyy[k] = -g_y_0_xxzzz_xyyyyy[k] * ab_x + g_y_0_xxzzz_xxyyyyy[k];

                g_y_0_xxxzzz_xyyyyz[k] = -g_y_0_xxzzz_xyyyyz[k] * ab_x + g_y_0_xxzzz_xxyyyyz[k];

                g_y_0_xxxzzz_xyyyzz[k] = -g_y_0_xxzzz_xyyyzz[k] * ab_x + g_y_0_xxzzz_xxyyyzz[k];

                g_y_0_xxxzzz_xyyzzz[k] = -g_y_0_xxzzz_xyyzzz[k] * ab_x + g_y_0_xxzzz_xxyyzzz[k];

                g_y_0_xxxzzz_xyzzzz[k] = -g_y_0_xxzzz_xyzzzz[k] * ab_x + g_y_0_xxzzz_xxyzzzz[k];

                g_y_0_xxxzzz_xzzzzz[k] = -g_y_0_xxzzz_xzzzzz[k] * ab_x + g_y_0_xxzzz_xxzzzzz[k];

                g_y_0_xxxzzz_yyyyyy[k] = -g_y_0_xxzzz_yyyyyy[k] * ab_x + g_y_0_xxzzz_xyyyyyy[k];

                g_y_0_xxxzzz_yyyyyz[k] = -g_y_0_xxzzz_yyyyyz[k] * ab_x + g_y_0_xxzzz_xyyyyyz[k];

                g_y_0_xxxzzz_yyyyzz[k] = -g_y_0_xxzzz_yyyyzz[k] * ab_x + g_y_0_xxzzz_xyyyyzz[k];

                g_y_0_xxxzzz_yyyzzz[k] = -g_y_0_xxzzz_yyyzzz[k] * ab_x + g_y_0_xxzzz_xyyyzzz[k];

                g_y_0_xxxzzz_yyzzzz[k] = -g_y_0_xxzzz_yyzzzz[k] * ab_x + g_y_0_xxzzz_xyyzzzz[k];

                g_y_0_xxxzzz_yzzzzz[k] = -g_y_0_xxzzz_yzzzzz[k] * ab_x + g_y_0_xxzzz_xyzzzzz[k];

                g_y_0_xxxzzz_zzzzzz[k] = -g_y_0_xxzzz_zzzzzz[k] * ab_x + g_y_0_xxzzz_xzzzzzz[k];
            }

            /// Set up 1064-1092 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 1064 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 1065 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 1066 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 1067 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 1068 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 1069 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 1070 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 1071 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 1072 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 1073 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 1074 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 1075 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 1076 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 1077 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 1078 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 1079 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 1080 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 1081 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 1082 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 1083 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 1084 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 1085 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 1086 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 1087 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 1088 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 1089 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 1090 * ccomps * dcomps);

            auto g_y_0_xxyyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 1091 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyyy_xxxxxx, g_y_0_xxyyyy_xxxxxy, g_y_0_xxyyyy_xxxxxz, g_y_0_xxyyyy_xxxxyy, g_y_0_xxyyyy_xxxxyz, g_y_0_xxyyyy_xxxxzz, g_y_0_xxyyyy_xxxyyy, g_y_0_xxyyyy_xxxyyz, g_y_0_xxyyyy_xxxyzz, g_y_0_xxyyyy_xxxzzz, g_y_0_xxyyyy_xxyyyy, g_y_0_xxyyyy_xxyyyz, g_y_0_xxyyyy_xxyyzz, g_y_0_xxyyyy_xxyzzz, g_y_0_xxyyyy_xxzzzz, g_y_0_xxyyyy_xyyyyy, g_y_0_xxyyyy_xyyyyz, g_y_0_xxyyyy_xyyyzz, g_y_0_xxyyyy_xyyzzz, g_y_0_xxyyyy_xyzzzz, g_y_0_xxyyyy_xzzzzz, g_y_0_xxyyyy_yyyyyy, g_y_0_xxyyyy_yyyyyz, g_y_0_xxyyyy_yyyyzz, g_y_0_xxyyyy_yyyzzz, g_y_0_xxyyyy_yyzzzz, g_y_0_xxyyyy_yzzzzz, g_y_0_xxyyyy_zzzzzz, g_y_0_xyyyy_xxxxxx, g_y_0_xyyyy_xxxxxxx, g_y_0_xyyyy_xxxxxxy, g_y_0_xyyyy_xxxxxxz, g_y_0_xyyyy_xxxxxy, g_y_0_xyyyy_xxxxxyy, g_y_0_xyyyy_xxxxxyz, g_y_0_xyyyy_xxxxxz, g_y_0_xyyyy_xxxxxzz, g_y_0_xyyyy_xxxxyy, g_y_0_xyyyy_xxxxyyy, g_y_0_xyyyy_xxxxyyz, g_y_0_xyyyy_xxxxyz, g_y_0_xyyyy_xxxxyzz, g_y_0_xyyyy_xxxxzz, g_y_0_xyyyy_xxxxzzz, g_y_0_xyyyy_xxxyyy, g_y_0_xyyyy_xxxyyyy, g_y_0_xyyyy_xxxyyyz, g_y_0_xyyyy_xxxyyz, g_y_0_xyyyy_xxxyyzz, g_y_0_xyyyy_xxxyzz, g_y_0_xyyyy_xxxyzzz, g_y_0_xyyyy_xxxzzz, g_y_0_xyyyy_xxxzzzz, g_y_0_xyyyy_xxyyyy, g_y_0_xyyyy_xxyyyyy, g_y_0_xyyyy_xxyyyyz, g_y_0_xyyyy_xxyyyz, g_y_0_xyyyy_xxyyyzz, g_y_0_xyyyy_xxyyzz, g_y_0_xyyyy_xxyyzzz, g_y_0_xyyyy_xxyzzz, g_y_0_xyyyy_xxyzzzz, g_y_0_xyyyy_xxzzzz, g_y_0_xyyyy_xxzzzzz, g_y_0_xyyyy_xyyyyy, g_y_0_xyyyy_xyyyyyy, g_y_0_xyyyy_xyyyyyz, g_y_0_xyyyy_xyyyyz, g_y_0_xyyyy_xyyyyzz, g_y_0_xyyyy_xyyyzz, g_y_0_xyyyy_xyyyzzz, g_y_0_xyyyy_xyyzzz, g_y_0_xyyyy_xyyzzzz, g_y_0_xyyyy_xyzzzz, g_y_0_xyyyy_xyzzzzz, g_y_0_xyyyy_xzzzzz, g_y_0_xyyyy_xzzzzzz, g_y_0_xyyyy_yyyyyy, g_y_0_xyyyy_yyyyyz, g_y_0_xyyyy_yyyyzz, g_y_0_xyyyy_yyyzzz, g_y_0_xyyyy_yyzzzz, g_y_0_xyyyy_yzzzzz, g_y_0_xyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_xxxxxx[k] = -g_y_0_xyyyy_xxxxxx[k] * ab_x + g_y_0_xyyyy_xxxxxxx[k];

                g_y_0_xxyyyy_xxxxxy[k] = -g_y_0_xyyyy_xxxxxy[k] * ab_x + g_y_0_xyyyy_xxxxxxy[k];

                g_y_0_xxyyyy_xxxxxz[k] = -g_y_0_xyyyy_xxxxxz[k] * ab_x + g_y_0_xyyyy_xxxxxxz[k];

                g_y_0_xxyyyy_xxxxyy[k] = -g_y_0_xyyyy_xxxxyy[k] * ab_x + g_y_0_xyyyy_xxxxxyy[k];

                g_y_0_xxyyyy_xxxxyz[k] = -g_y_0_xyyyy_xxxxyz[k] * ab_x + g_y_0_xyyyy_xxxxxyz[k];

                g_y_0_xxyyyy_xxxxzz[k] = -g_y_0_xyyyy_xxxxzz[k] * ab_x + g_y_0_xyyyy_xxxxxzz[k];

                g_y_0_xxyyyy_xxxyyy[k] = -g_y_0_xyyyy_xxxyyy[k] * ab_x + g_y_0_xyyyy_xxxxyyy[k];

                g_y_0_xxyyyy_xxxyyz[k] = -g_y_0_xyyyy_xxxyyz[k] * ab_x + g_y_0_xyyyy_xxxxyyz[k];

                g_y_0_xxyyyy_xxxyzz[k] = -g_y_0_xyyyy_xxxyzz[k] * ab_x + g_y_0_xyyyy_xxxxyzz[k];

                g_y_0_xxyyyy_xxxzzz[k] = -g_y_0_xyyyy_xxxzzz[k] * ab_x + g_y_0_xyyyy_xxxxzzz[k];

                g_y_0_xxyyyy_xxyyyy[k] = -g_y_0_xyyyy_xxyyyy[k] * ab_x + g_y_0_xyyyy_xxxyyyy[k];

                g_y_0_xxyyyy_xxyyyz[k] = -g_y_0_xyyyy_xxyyyz[k] * ab_x + g_y_0_xyyyy_xxxyyyz[k];

                g_y_0_xxyyyy_xxyyzz[k] = -g_y_0_xyyyy_xxyyzz[k] * ab_x + g_y_0_xyyyy_xxxyyzz[k];

                g_y_0_xxyyyy_xxyzzz[k] = -g_y_0_xyyyy_xxyzzz[k] * ab_x + g_y_0_xyyyy_xxxyzzz[k];

                g_y_0_xxyyyy_xxzzzz[k] = -g_y_0_xyyyy_xxzzzz[k] * ab_x + g_y_0_xyyyy_xxxzzzz[k];

                g_y_0_xxyyyy_xyyyyy[k] = -g_y_0_xyyyy_xyyyyy[k] * ab_x + g_y_0_xyyyy_xxyyyyy[k];

                g_y_0_xxyyyy_xyyyyz[k] = -g_y_0_xyyyy_xyyyyz[k] * ab_x + g_y_0_xyyyy_xxyyyyz[k];

                g_y_0_xxyyyy_xyyyzz[k] = -g_y_0_xyyyy_xyyyzz[k] * ab_x + g_y_0_xyyyy_xxyyyzz[k];

                g_y_0_xxyyyy_xyyzzz[k] = -g_y_0_xyyyy_xyyzzz[k] * ab_x + g_y_0_xyyyy_xxyyzzz[k];

                g_y_0_xxyyyy_xyzzzz[k] = -g_y_0_xyyyy_xyzzzz[k] * ab_x + g_y_0_xyyyy_xxyzzzz[k];

                g_y_0_xxyyyy_xzzzzz[k] = -g_y_0_xyyyy_xzzzzz[k] * ab_x + g_y_0_xyyyy_xxzzzzz[k];

                g_y_0_xxyyyy_yyyyyy[k] = -g_y_0_xyyyy_yyyyyy[k] * ab_x + g_y_0_xyyyy_xyyyyyy[k];

                g_y_0_xxyyyy_yyyyyz[k] = -g_y_0_xyyyy_yyyyyz[k] * ab_x + g_y_0_xyyyy_xyyyyyz[k];

                g_y_0_xxyyyy_yyyyzz[k] = -g_y_0_xyyyy_yyyyzz[k] * ab_x + g_y_0_xyyyy_xyyyyzz[k];

                g_y_0_xxyyyy_yyyzzz[k] = -g_y_0_xyyyy_yyyzzz[k] * ab_x + g_y_0_xyyyy_xyyyzzz[k];

                g_y_0_xxyyyy_yyzzzz[k] = -g_y_0_xyyyy_yyzzzz[k] * ab_x + g_y_0_xyyyy_xyyzzzz[k];

                g_y_0_xxyyyy_yzzzzz[k] = -g_y_0_xyyyy_yzzzzz[k] * ab_x + g_y_0_xyyyy_xyzzzzz[k];

                g_y_0_xxyyyy_zzzzzz[k] = -g_y_0_xyyyy_zzzzzz[k] * ab_x + g_y_0_xyyyy_xzzzzzz[k];
            }

            /// Set up 1092-1120 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 1092 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 1093 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 1094 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 1095 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 1096 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 1097 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 1098 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 1099 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 1100 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 1101 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 1102 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 1103 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 1104 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 1105 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 1106 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 1107 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 1108 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 1109 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 1110 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 1111 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 1112 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 1113 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 1114 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 1115 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 1116 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 1117 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 1118 * ccomps * dcomps);

            auto g_y_0_xxyyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 1119 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyyz_xxxxxx, g_y_0_xxyyyz_xxxxxy, g_y_0_xxyyyz_xxxxxz, g_y_0_xxyyyz_xxxxyy, g_y_0_xxyyyz_xxxxyz, g_y_0_xxyyyz_xxxxzz, g_y_0_xxyyyz_xxxyyy, g_y_0_xxyyyz_xxxyyz, g_y_0_xxyyyz_xxxyzz, g_y_0_xxyyyz_xxxzzz, g_y_0_xxyyyz_xxyyyy, g_y_0_xxyyyz_xxyyyz, g_y_0_xxyyyz_xxyyzz, g_y_0_xxyyyz_xxyzzz, g_y_0_xxyyyz_xxzzzz, g_y_0_xxyyyz_xyyyyy, g_y_0_xxyyyz_xyyyyz, g_y_0_xxyyyz_xyyyzz, g_y_0_xxyyyz_xyyzzz, g_y_0_xxyyyz_xyzzzz, g_y_0_xxyyyz_xzzzzz, g_y_0_xxyyyz_yyyyyy, g_y_0_xxyyyz_yyyyyz, g_y_0_xxyyyz_yyyyzz, g_y_0_xxyyyz_yyyzzz, g_y_0_xxyyyz_yyzzzz, g_y_0_xxyyyz_yzzzzz, g_y_0_xxyyyz_zzzzzz, g_y_0_xyyyz_xxxxxx, g_y_0_xyyyz_xxxxxxx, g_y_0_xyyyz_xxxxxxy, g_y_0_xyyyz_xxxxxxz, g_y_0_xyyyz_xxxxxy, g_y_0_xyyyz_xxxxxyy, g_y_0_xyyyz_xxxxxyz, g_y_0_xyyyz_xxxxxz, g_y_0_xyyyz_xxxxxzz, g_y_0_xyyyz_xxxxyy, g_y_0_xyyyz_xxxxyyy, g_y_0_xyyyz_xxxxyyz, g_y_0_xyyyz_xxxxyz, g_y_0_xyyyz_xxxxyzz, g_y_0_xyyyz_xxxxzz, g_y_0_xyyyz_xxxxzzz, g_y_0_xyyyz_xxxyyy, g_y_0_xyyyz_xxxyyyy, g_y_0_xyyyz_xxxyyyz, g_y_0_xyyyz_xxxyyz, g_y_0_xyyyz_xxxyyzz, g_y_0_xyyyz_xxxyzz, g_y_0_xyyyz_xxxyzzz, g_y_0_xyyyz_xxxzzz, g_y_0_xyyyz_xxxzzzz, g_y_0_xyyyz_xxyyyy, g_y_0_xyyyz_xxyyyyy, g_y_0_xyyyz_xxyyyyz, g_y_0_xyyyz_xxyyyz, g_y_0_xyyyz_xxyyyzz, g_y_0_xyyyz_xxyyzz, g_y_0_xyyyz_xxyyzzz, g_y_0_xyyyz_xxyzzz, g_y_0_xyyyz_xxyzzzz, g_y_0_xyyyz_xxzzzz, g_y_0_xyyyz_xxzzzzz, g_y_0_xyyyz_xyyyyy, g_y_0_xyyyz_xyyyyyy, g_y_0_xyyyz_xyyyyyz, g_y_0_xyyyz_xyyyyz, g_y_0_xyyyz_xyyyyzz, g_y_0_xyyyz_xyyyzz, g_y_0_xyyyz_xyyyzzz, g_y_0_xyyyz_xyyzzz, g_y_0_xyyyz_xyyzzzz, g_y_0_xyyyz_xyzzzz, g_y_0_xyyyz_xyzzzzz, g_y_0_xyyyz_xzzzzz, g_y_0_xyyyz_xzzzzzz, g_y_0_xyyyz_yyyyyy, g_y_0_xyyyz_yyyyyz, g_y_0_xyyyz_yyyyzz, g_y_0_xyyyz_yyyzzz, g_y_0_xyyyz_yyzzzz, g_y_0_xyyyz_yzzzzz, g_y_0_xyyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_xxxxxx[k] = -g_y_0_xyyyz_xxxxxx[k] * ab_x + g_y_0_xyyyz_xxxxxxx[k];

                g_y_0_xxyyyz_xxxxxy[k] = -g_y_0_xyyyz_xxxxxy[k] * ab_x + g_y_0_xyyyz_xxxxxxy[k];

                g_y_0_xxyyyz_xxxxxz[k] = -g_y_0_xyyyz_xxxxxz[k] * ab_x + g_y_0_xyyyz_xxxxxxz[k];

                g_y_0_xxyyyz_xxxxyy[k] = -g_y_0_xyyyz_xxxxyy[k] * ab_x + g_y_0_xyyyz_xxxxxyy[k];

                g_y_0_xxyyyz_xxxxyz[k] = -g_y_0_xyyyz_xxxxyz[k] * ab_x + g_y_0_xyyyz_xxxxxyz[k];

                g_y_0_xxyyyz_xxxxzz[k] = -g_y_0_xyyyz_xxxxzz[k] * ab_x + g_y_0_xyyyz_xxxxxzz[k];

                g_y_0_xxyyyz_xxxyyy[k] = -g_y_0_xyyyz_xxxyyy[k] * ab_x + g_y_0_xyyyz_xxxxyyy[k];

                g_y_0_xxyyyz_xxxyyz[k] = -g_y_0_xyyyz_xxxyyz[k] * ab_x + g_y_0_xyyyz_xxxxyyz[k];

                g_y_0_xxyyyz_xxxyzz[k] = -g_y_0_xyyyz_xxxyzz[k] * ab_x + g_y_0_xyyyz_xxxxyzz[k];

                g_y_0_xxyyyz_xxxzzz[k] = -g_y_0_xyyyz_xxxzzz[k] * ab_x + g_y_0_xyyyz_xxxxzzz[k];

                g_y_0_xxyyyz_xxyyyy[k] = -g_y_0_xyyyz_xxyyyy[k] * ab_x + g_y_0_xyyyz_xxxyyyy[k];

                g_y_0_xxyyyz_xxyyyz[k] = -g_y_0_xyyyz_xxyyyz[k] * ab_x + g_y_0_xyyyz_xxxyyyz[k];

                g_y_0_xxyyyz_xxyyzz[k] = -g_y_0_xyyyz_xxyyzz[k] * ab_x + g_y_0_xyyyz_xxxyyzz[k];

                g_y_0_xxyyyz_xxyzzz[k] = -g_y_0_xyyyz_xxyzzz[k] * ab_x + g_y_0_xyyyz_xxxyzzz[k];

                g_y_0_xxyyyz_xxzzzz[k] = -g_y_0_xyyyz_xxzzzz[k] * ab_x + g_y_0_xyyyz_xxxzzzz[k];

                g_y_0_xxyyyz_xyyyyy[k] = -g_y_0_xyyyz_xyyyyy[k] * ab_x + g_y_0_xyyyz_xxyyyyy[k];

                g_y_0_xxyyyz_xyyyyz[k] = -g_y_0_xyyyz_xyyyyz[k] * ab_x + g_y_0_xyyyz_xxyyyyz[k];

                g_y_0_xxyyyz_xyyyzz[k] = -g_y_0_xyyyz_xyyyzz[k] * ab_x + g_y_0_xyyyz_xxyyyzz[k];

                g_y_0_xxyyyz_xyyzzz[k] = -g_y_0_xyyyz_xyyzzz[k] * ab_x + g_y_0_xyyyz_xxyyzzz[k];

                g_y_0_xxyyyz_xyzzzz[k] = -g_y_0_xyyyz_xyzzzz[k] * ab_x + g_y_0_xyyyz_xxyzzzz[k];

                g_y_0_xxyyyz_xzzzzz[k] = -g_y_0_xyyyz_xzzzzz[k] * ab_x + g_y_0_xyyyz_xxzzzzz[k];

                g_y_0_xxyyyz_yyyyyy[k] = -g_y_0_xyyyz_yyyyyy[k] * ab_x + g_y_0_xyyyz_xyyyyyy[k];

                g_y_0_xxyyyz_yyyyyz[k] = -g_y_0_xyyyz_yyyyyz[k] * ab_x + g_y_0_xyyyz_xyyyyyz[k];

                g_y_0_xxyyyz_yyyyzz[k] = -g_y_0_xyyyz_yyyyzz[k] * ab_x + g_y_0_xyyyz_xyyyyzz[k];

                g_y_0_xxyyyz_yyyzzz[k] = -g_y_0_xyyyz_yyyzzz[k] * ab_x + g_y_0_xyyyz_xyyyzzz[k];

                g_y_0_xxyyyz_yyzzzz[k] = -g_y_0_xyyyz_yyzzzz[k] * ab_x + g_y_0_xyyyz_xyyzzzz[k];

                g_y_0_xxyyyz_yzzzzz[k] = -g_y_0_xyyyz_yzzzzz[k] * ab_x + g_y_0_xyyyz_xyzzzzz[k];

                g_y_0_xxyyyz_zzzzzz[k] = -g_y_0_xyyyz_zzzzzz[k] * ab_x + g_y_0_xyyyz_xzzzzzz[k];
            }

            /// Set up 1120-1148 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1120 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1121 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1122 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1123 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1124 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1125 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1126 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1127 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1128 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1129 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1130 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1131 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1132 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1133 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1134 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1135 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1136 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1137 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1138 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1139 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1140 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1141 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1142 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1143 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1144 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1145 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1146 * ccomps * dcomps);

            auto g_y_0_xxyyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1147 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyzz_xxxxxx, g_y_0_xxyyzz_xxxxxy, g_y_0_xxyyzz_xxxxxz, g_y_0_xxyyzz_xxxxyy, g_y_0_xxyyzz_xxxxyz, g_y_0_xxyyzz_xxxxzz, g_y_0_xxyyzz_xxxyyy, g_y_0_xxyyzz_xxxyyz, g_y_0_xxyyzz_xxxyzz, g_y_0_xxyyzz_xxxzzz, g_y_0_xxyyzz_xxyyyy, g_y_0_xxyyzz_xxyyyz, g_y_0_xxyyzz_xxyyzz, g_y_0_xxyyzz_xxyzzz, g_y_0_xxyyzz_xxzzzz, g_y_0_xxyyzz_xyyyyy, g_y_0_xxyyzz_xyyyyz, g_y_0_xxyyzz_xyyyzz, g_y_0_xxyyzz_xyyzzz, g_y_0_xxyyzz_xyzzzz, g_y_0_xxyyzz_xzzzzz, g_y_0_xxyyzz_yyyyyy, g_y_0_xxyyzz_yyyyyz, g_y_0_xxyyzz_yyyyzz, g_y_0_xxyyzz_yyyzzz, g_y_0_xxyyzz_yyzzzz, g_y_0_xxyyzz_yzzzzz, g_y_0_xxyyzz_zzzzzz, g_y_0_xyyzz_xxxxxx, g_y_0_xyyzz_xxxxxxx, g_y_0_xyyzz_xxxxxxy, g_y_0_xyyzz_xxxxxxz, g_y_0_xyyzz_xxxxxy, g_y_0_xyyzz_xxxxxyy, g_y_0_xyyzz_xxxxxyz, g_y_0_xyyzz_xxxxxz, g_y_0_xyyzz_xxxxxzz, g_y_0_xyyzz_xxxxyy, g_y_0_xyyzz_xxxxyyy, g_y_0_xyyzz_xxxxyyz, g_y_0_xyyzz_xxxxyz, g_y_0_xyyzz_xxxxyzz, g_y_0_xyyzz_xxxxzz, g_y_0_xyyzz_xxxxzzz, g_y_0_xyyzz_xxxyyy, g_y_0_xyyzz_xxxyyyy, g_y_0_xyyzz_xxxyyyz, g_y_0_xyyzz_xxxyyz, g_y_0_xyyzz_xxxyyzz, g_y_0_xyyzz_xxxyzz, g_y_0_xyyzz_xxxyzzz, g_y_0_xyyzz_xxxzzz, g_y_0_xyyzz_xxxzzzz, g_y_0_xyyzz_xxyyyy, g_y_0_xyyzz_xxyyyyy, g_y_0_xyyzz_xxyyyyz, g_y_0_xyyzz_xxyyyz, g_y_0_xyyzz_xxyyyzz, g_y_0_xyyzz_xxyyzz, g_y_0_xyyzz_xxyyzzz, g_y_0_xyyzz_xxyzzz, g_y_0_xyyzz_xxyzzzz, g_y_0_xyyzz_xxzzzz, g_y_0_xyyzz_xxzzzzz, g_y_0_xyyzz_xyyyyy, g_y_0_xyyzz_xyyyyyy, g_y_0_xyyzz_xyyyyyz, g_y_0_xyyzz_xyyyyz, g_y_0_xyyzz_xyyyyzz, g_y_0_xyyzz_xyyyzz, g_y_0_xyyzz_xyyyzzz, g_y_0_xyyzz_xyyzzz, g_y_0_xyyzz_xyyzzzz, g_y_0_xyyzz_xyzzzz, g_y_0_xyyzz_xyzzzzz, g_y_0_xyyzz_xzzzzz, g_y_0_xyyzz_xzzzzzz, g_y_0_xyyzz_yyyyyy, g_y_0_xyyzz_yyyyyz, g_y_0_xyyzz_yyyyzz, g_y_0_xyyzz_yyyzzz, g_y_0_xyyzz_yyzzzz, g_y_0_xyyzz_yzzzzz, g_y_0_xyyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_xxxxxx[k] = -g_y_0_xyyzz_xxxxxx[k] * ab_x + g_y_0_xyyzz_xxxxxxx[k];

                g_y_0_xxyyzz_xxxxxy[k] = -g_y_0_xyyzz_xxxxxy[k] * ab_x + g_y_0_xyyzz_xxxxxxy[k];

                g_y_0_xxyyzz_xxxxxz[k] = -g_y_0_xyyzz_xxxxxz[k] * ab_x + g_y_0_xyyzz_xxxxxxz[k];

                g_y_0_xxyyzz_xxxxyy[k] = -g_y_0_xyyzz_xxxxyy[k] * ab_x + g_y_0_xyyzz_xxxxxyy[k];

                g_y_0_xxyyzz_xxxxyz[k] = -g_y_0_xyyzz_xxxxyz[k] * ab_x + g_y_0_xyyzz_xxxxxyz[k];

                g_y_0_xxyyzz_xxxxzz[k] = -g_y_0_xyyzz_xxxxzz[k] * ab_x + g_y_0_xyyzz_xxxxxzz[k];

                g_y_0_xxyyzz_xxxyyy[k] = -g_y_0_xyyzz_xxxyyy[k] * ab_x + g_y_0_xyyzz_xxxxyyy[k];

                g_y_0_xxyyzz_xxxyyz[k] = -g_y_0_xyyzz_xxxyyz[k] * ab_x + g_y_0_xyyzz_xxxxyyz[k];

                g_y_0_xxyyzz_xxxyzz[k] = -g_y_0_xyyzz_xxxyzz[k] * ab_x + g_y_0_xyyzz_xxxxyzz[k];

                g_y_0_xxyyzz_xxxzzz[k] = -g_y_0_xyyzz_xxxzzz[k] * ab_x + g_y_0_xyyzz_xxxxzzz[k];

                g_y_0_xxyyzz_xxyyyy[k] = -g_y_0_xyyzz_xxyyyy[k] * ab_x + g_y_0_xyyzz_xxxyyyy[k];

                g_y_0_xxyyzz_xxyyyz[k] = -g_y_0_xyyzz_xxyyyz[k] * ab_x + g_y_0_xyyzz_xxxyyyz[k];

                g_y_0_xxyyzz_xxyyzz[k] = -g_y_0_xyyzz_xxyyzz[k] * ab_x + g_y_0_xyyzz_xxxyyzz[k];

                g_y_0_xxyyzz_xxyzzz[k] = -g_y_0_xyyzz_xxyzzz[k] * ab_x + g_y_0_xyyzz_xxxyzzz[k];

                g_y_0_xxyyzz_xxzzzz[k] = -g_y_0_xyyzz_xxzzzz[k] * ab_x + g_y_0_xyyzz_xxxzzzz[k];

                g_y_0_xxyyzz_xyyyyy[k] = -g_y_0_xyyzz_xyyyyy[k] * ab_x + g_y_0_xyyzz_xxyyyyy[k];

                g_y_0_xxyyzz_xyyyyz[k] = -g_y_0_xyyzz_xyyyyz[k] * ab_x + g_y_0_xyyzz_xxyyyyz[k];

                g_y_0_xxyyzz_xyyyzz[k] = -g_y_0_xyyzz_xyyyzz[k] * ab_x + g_y_0_xyyzz_xxyyyzz[k];

                g_y_0_xxyyzz_xyyzzz[k] = -g_y_0_xyyzz_xyyzzz[k] * ab_x + g_y_0_xyyzz_xxyyzzz[k];

                g_y_0_xxyyzz_xyzzzz[k] = -g_y_0_xyyzz_xyzzzz[k] * ab_x + g_y_0_xyyzz_xxyzzzz[k];

                g_y_0_xxyyzz_xzzzzz[k] = -g_y_0_xyyzz_xzzzzz[k] * ab_x + g_y_0_xyyzz_xxzzzzz[k];

                g_y_0_xxyyzz_yyyyyy[k] = -g_y_0_xyyzz_yyyyyy[k] * ab_x + g_y_0_xyyzz_xyyyyyy[k];

                g_y_0_xxyyzz_yyyyyz[k] = -g_y_0_xyyzz_yyyyyz[k] * ab_x + g_y_0_xyyzz_xyyyyyz[k];

                g_y_0_xxyyzz_yyyyzz[k] = -g_y_0_xyyzz_yyyyzz[k] * ab_x + g_y_0_xyyzz_xyyyyzz[k];

                g_y_0_xxyyzz_yyyzzz[k] = -g_y_0_xyyzz_yyyzzz[k] * ab_x + g_y_0_xyyzz_xyyyzzz[k];

                g_y_0_xxyyzz_yyzzzz[k] = -g_y_0_xyyzz_yyzzzz[k] * ab_x + g_y_0_xyyzz_xyyzzzz[k];

                g_y_0_xxyyzz_yzzzzz[k] = -g_y_0_xyyzz_yzzzzz[k] * ab_x + g_y_0_xyyzz_xyzzzzz[k];

                g_y_0_xxyyzz_zzzzzz[k] = -g_y_0_xyyzz_zzzzzz[k] * ab_x + g_y_0_xyyzz_xzzzzzz[k];
            }

            /// Set up 1148-1176 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1148 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1149 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1150 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1151 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1152 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1153 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1154 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1155 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1156 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1157 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1158 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1159 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1160 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1161 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1162 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1163 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1164 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1165 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1166 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1167 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1168 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1169 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1170 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1171 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1172 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1173 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1174 * ccomps * dcomps);

            auto g_y_0_xxyzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1175 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyzzz_xxxxxx, g_y_0_xxyzzz_xxxxxy, g_y_0_xxyzzz_xxxxxz, g_y_0_xxyzzz_xxxxyy, g_y_0_xxyzzz_xxxxyz, g_y_0_xxyzzz_xxxxzz, g_y_0_xxyzzz_xxxyyy, g_y_0_xxyzzz_xxxyyz, g_y_0_xxyzzz_xxxyzz, g_y_0_xxyzzz_xxxzzz, g_y_0_xxyzzz_xxyyyy, g_y_0_xxyzzz_xxyyyz, g_y_0_xxyzzz_xxyyzz, g_y_0_xxyzzz_xxyzzz, g_y_0_xxyzzz_xxzzzz, g_y_0_xxyzzz_xyyyyy, g_y_0_xxyzzz_xyyyyz, g_y_0_xxyzzz_xyyyzz, g_y_0_xxyzzz_xyyzzz, g_y_0_xxyzzz_xyzzzz, g_y_0_xxyzzz_xzzzzz, g_y_0_xxyzzz_yyyyyy, g_y_0_xxyzzz_yyyyyz, g_y_0_xxyzzz_yyyyzz, g_y_0_xxyzzz_yyyzzz, g_y_0_xxyzzz_yyzzzz, g_y_0_xxyzzz_yzzzzz, g_y_0_xxyzzz_zzzzzz, g_y_0_xyzzz_xxxxxx, g_y_0_xyzzz_xxxxxxx, g_y_0_xyzzz_xxxxxxy, g_y_0_xyzzz_xxxxxxz, g_y_0_xyzzz_xxxxxy, g_y_0_xyzzz_xxxxxyy, g_y_0_xyzzz_xxxxxyz, g_y_0_xyzzz_xxxxxz, g_y_0_xyzzz_xxxxxzz, g_y_0_xyzzz_xxxxyy, g_y_0_xyzzz_xxxxyyy, g_y_0_xyzzz_xxxxyyz, g_y_0_xyzzz_xxxxyz, g_y_0_xyzzz_xxxxyzz, g_y_0_xyzzz_xxxxzz, g_y_0_xyzzz_xxxxzzz, g_y_0_xyzzz_xxxyyy, g_y_0_xyzzz_xxxyyyy, g_y_0_xyzzz_xxxyyyz, g_y_0_xyzzz_xxxyyz, g_y_0_xyzzz_xxxyyzz, g_y_0_xyzzz_xxxyzz, g_y_0_xyzzz_xxxyzzz, g_y_0_xyzzz_xxxzzz, g_y_0_xyzzz_xxxzzzz, g_y_0_xyzzz_xxyyyy, g_y_0_xyzzz_xxyyyyy, g_y_0_xyzzz_xxyyyyz, g_y_0_xyzzz_xxyyyz, g_y_0_xyzzz_xxyyyzz, g_y_0_xyzzz_xxyyzz, g_y_0_xyzzz_xxyyzzz, g_y_0_xyzzz_xxyzzz, g_y_0_xyzzz_xxyzzzz, g_y_0_xyzzz_xxzzzz, g_y_0_xyzzz_xxzzzzz, g_y_0_xyzzz_xyyyyy, g_y_0_xyzzz_xyyyyyy, g_y_0_xyzzz_xyyyyyz, g_y_0_xyzzz_xyyyyz, g_y_0_xyzzz_xyyyyzz, g_y_0_xyzzz_xyyyzz, g_y_0_xyzzz_xyyyzzz, g_y_0_xyzzz_xyyzzz, g_y_0_xyzzz_xyyzzzz, g_y_0_xyzzz_xyzzzz, g_y_0_xyzzz_xyzzzzz, g_y_0_xyzzz_xzzzzz, g_y_0_xyzzz_xzzzzzz, g_y_0_xyzzz_yyyyyy, g_y_0_xyzzz_yyyyyz, g_y_0_xyzzz_yyyyzz, g_y_0_xyzzz_yyyzzz, g_y_0_xyzzz_yyzzzz, g_y_0_xyzzz_yzzzzz, g_y_0_xyzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_xxxxxx[k] = -g_y_0_xyzzz_xxxxxx[k] * ab_x + g_y_0_xyzzz_xxxxxxx[k];

                g_y_0_xxyzzz_xxxxxy[k] = -g_y_0_xyzzz_xxxxxy[k] * ab_x + g_y_0_xyzzz_xxxxxxy[k];

                g_y_0_xxyzzz_xxxxxz[k] = -g_y_0_xyzzz_xxxxxz[k] * ab_x + g_y_0_xyzzz_xxxxxxz[k];

                g_y_0_xxyzzz_xxxxyy[k] = -g_y_0_xyzzz_xxxxyy[k] * ab_x + g_y_0_xyzzz_xxxxxyy[k];

                g_y_0_xxyzzz_xxxxyz[k] = -g_y_0_xyzzz_xxxxyz[k] * ab_x + g_y_0_xyzzz_xxxxxyz[k];

                g_y_0_xxyzzz_xxxxzz[k] = -g_y_0_xyzzz_xxxxzz[k] * ab_x + g_y_0_xyzzz_xxxxxzz[k];

                g_y_0_xxyzzz_xxxyyy[k] = -g_y_0_xyzzz_xxxyyy[k] * ab_x + g_y_0_xyzzz_xxxxyyy[k];

                g_y_0_xxyzzz_xxxyyz[k] = -g_y_0_xyzzz_xxxyyz[k] * ab_x + g_y_0_xyzzz_xxxxyyz[k];

                g_y_0_xxyzzz_xxxyzz[k] = -g_y_0_xyzzz_xxxyzz[k] * ab_x + g_y_0_xyzzz_xxxxyzz[k];

                g_y_0_xxyzzz_xxxzzz[k] = -g_y_0_xyzzz_xxxzzz[k] * ab_x + g_y_0_xyzzz_xxxxzzz[k];

                g_y_0_xxyzzz_xxyyyy[k] = -g_y_0_xyzzz_xxyyyy[k] * ab_x + g_y_0_xyzzz_xxxyyyy[k];

                g_y_0_xxyzzz_xxyyyz[k] = -g_y_0_xyzzz_xxyyyz[k] * ab_x + g_y_0_xyzzz_xxxyyyz[k];

                g_y_0_xxyzzz_xxyyzz[k] = -g_y_0_xyzzz_xxyyzz[k] * ab_x + g_y_0_xyzzz_xxxyyzz[k];

                g_y_0_xxyzzz_xxyzzz[k] = -g_y_0_xyzzz_xxyzzz[k] * ab_x + g_y_0_xyzzz_xxxyzzz[k];

                g_y_0_xxyzzz_xxzzzz[k] = -g_y_0_xyzzz_xxzzzz[k] * ab_x + g_y_0_xyzzz_xxxzzzz[k];

                g_y_0_xxyzzz_xyyyyy[k] = -g_y_0_xyzzz_xyyyyy[k] * ab_x + g_y_0_xyzzz_xxyyyyy[k];

                g_y_0_xxyzzz_xyyyyz[k] = -g_y_0_xyzzz_xyyyyz[k] * ab_x + g_y_0_xyzzz_xxyyyyz[k];

                g_y_0_xxyzzz_xyyyzz[k] = -g_y_0_xyzzz_xyyyzz[k] * ab_x + g_y_0_xyzzz_xxyyyzz[k];

                g_y_0_xxyzzz_xyyzzz[k] = -g_y_0_xyzzz_xyyzzz[k] * ab_x + g_y_0_xyzzz_xxyyzzz[k];

                g_y_0_xxyzzz_xyzzzz[k] = -g_y_0_xyzzz_xyzzzz[k] * ab_x + g_y_0_xyzzz_xxyzzzz[k];

                g_y_0_xxyzzz_xzzzzz[k] = -g_y_0_xyzzz_xzzzzz[k] * ab_x + g_y_0_xyzzz_xxzzzzz[k];

                g_y_0_xxyzzz_yyyyyy[k] = -g_y_0_xyzzz_yyyyyy[k] * ab_x + g_y_0_xyzzz_xyyyyyy[k];

                g_y_0_xxyzzz_yyyyyz[k] = -g_y_0_xyzzz_yyyyyz[k] * ab_x + g_y_0_xyzzz_xyyyyyz[k];

                g_y_0_xxyzzz_yyyyzz[k] = -g_y_0_xyzzz_yyyyzz[k] * ab_x + g_y_0_xyzzz_xyyyyzz[k];

                g_y_0_xxyzzz_yyyzzz[k] = -g_y_0_xyzzz_yyyzzz[k] * ab_x + g_y_0_xyzzz_xyyyzzz[k];

                g_y_0_xxyzzz_yyzzzz[k] = -g_y_0_xyzzz_yyzzzz[k] * ab_x + g_y_0_xyzzz_xyyzzzz[k];

                g_y_0_xxyzzz_yzzzzz[k] = -g_y_0_xyzzz_yzzzzz[k] * ab_x + g_y_0_xyzzz_xyzzzzz[k];

                g_y_0_xxyzzz_zzzzzz[k] = -g_y_0_xyzzz_zzzzzz[k] * ab_x + g_y_0_xyzzz_xzzzzzz[k];
            }

            /// Set up 1176-1204 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1176 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1177 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1178 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1179 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1180 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1181 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1182 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1183 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1184 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1185 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1186 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1187 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1188 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1189 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1190 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1191 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1192 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1193 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1194 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1195 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1196 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1197 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1198 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1199 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1200 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1201 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1202 * ccomps * dcomps);

            auto g_y_0_xxzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1203 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzzzz_xxxxxx, g_y_0_xxzzzz_xxxxxy, g_y_0_xxzzzz_xxxxxz, g_y_0_xxzzzz_xxxxyy, g_y_0_xxzzzz_xxxxyz, g_y_0_xxzzzz_xxxxzz, g_y_0_xxzzzz_xxxyyy, g_y_0_xxzzzz_xxxyyz, g_y_0_xxzzzz_xxxyzz, g_y_0_xxzzzz_xxxzzz, g_y_0_xxzzzz_xxyyyy, g_y_0_xxzzzz_xxyyyz, g_y_0_xxzzzz_xxyyzz, g_y_0_xxzzzz_xxyzzz, g_y_0_xxzzzz_xxzzzz, g_y_0_xxzzzz_xyyyyy, g_y_0_xxzzzz_xyyyyz, g_y_0_xxzzzz_xyyyzz, g_y_0_xxzzzz_xyyzzz, g_y_0_xxzzzz_xyzzzz, g_y_0_xxzzzz_xzzzzz, g_y_0_xxzzzz_yyyyyy, g_y_0_xxzzzz_yyyyyz, g_y_0_xxzzzz_yyyyzz, g_y_0_xxzzzz_yyyzzz, g_y_0_xxzzzz_yyzzzz, g_y_0_xxzzzz_yzzzzz, g_y_0_xxzzzz_zzzzzz, g_y_0_xzzzz_xxxxxx, g_y_0_xzzzz_xxxxxxx, g_y_0_xzzzz_xxxxxxy, g_y_0_xzzzz_xxxxxxz, g_y_0_xzzzz_xxxxxy, g_y_0_xzzzz_xxxxxyy, g_y_0_xzzzz_xxxxxyz, g_y_0_xzzzz_xxxxxz, g_y_0_xzzzz_xxxxxzz, g_y_0_xzzzz_xxxxyy, g_y_0_xzzzz_xxxxyyy, g_y_0_xzzzz_xxxxyyz, g_y_0_xzzzz_xxxxyz, g_y_0_xzzzz_xxxxyzz, g_y_0_xzzzz_xxxxzz, g_y_0_xzzzz_xxxxzzz, g_y_0_xzzzz_xxxyyy, g_y_0_xzzzz_xxxyyyy, g_y_0_xzzzz_xxxyyyz, g_y_0_xzzzz_xxxyyz, g_y_0_xzzzz_xxxyyzz, g_y_0_xzzzz_xxxyzz, g_y_0_xzzzz_xxxyzzz, g_y_0_xzzzz_xxxzzz, g_y_0_xzzzz_xxxzzzz, g_y_0_xzzzz_xxyyyy, g_y_0_xzzzz_xxyyyyy, g_y_0_xzzzz_xxyyyyz, g_y_0_xzzzz_xxyyyz, g_y_0_xzzzz_xxyyyzz, g_y_0_xzzzz_xxyyzz, g_y_0_xzzzz_xxyyzzz, g_y_0_xzzzz_xxyzzz, g_y_0_xzzzz_xxyzzzz, g_y_0_xzzzz_xxzzzz, g_y_0_xzzzz_xxzzzzz, g_y_0_xzzzz_xyyyyy, g_y_0_xzzzz_xyyyyyy, g_y_0_xzzzz_xyyyyyz, g_y_0_xzzzz_xyyyyz, g_y_0_xzzzz_xyyyyzz, g_y_0_xzzzz_xyyyzz, g_y_0_xzzzz_xyyyzzz, g_y_0_xzzzz_xyyzzz, g_y_0_xzzzz_xyyzzzz, g_y_0_xzzzz_xyzzzz, g_y_0_xzzzz_xyzzzzz, g_y_0_xzzzz_xzzzzz, g_y_0_xzzzz_xzzzzzz, g_y_0_xzzzz_yyyyyy, g_y_0_xzzzz_yyyyyz, g_y_0_xzzzz_yyyyzz, g_y_0_xzzzz_yyyzzz, g_y_0_xzzzz_yyzzzz, g_y_0_xzzzz_yzzzzz, g_y_0_xzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_xxxxxx[k] = -g_y_0_xzzzz_xxxxxx[k] * ab_x + g_y_0_xzzzz_xxxxxxx[k];

                g_y_0_xxzzzz_xxxxxy[k] = -g_y_0_xzzzz_xxxxxy[k] * ab_x + g_y_0_xzzzz_xxxxxxy[k];

                g_y_0_xxzzzz_xxxxxz[k] = -g_y_0_xzzzz_xxxxxz[k] * ab_x + g_y_0_xzzzz_xxxxxxz[k];

                g_y_0_xxzzzz_xxxxyy[k] = -g_y_0_xzzzz_xxxxyy[k] * ab_x + g_y_0_xzzzz_xxxxxyy[k];

                g_y_0_xxzzzz_xxxxyz[k] = -g_y_0_xzzzz_xxxxyz[k] * ab_x + g_y_0_xzzzz_xxxxxyz[k];

                g_y_0_xxzzzz_xxxxzz[k] = -g_y_0_xzzzz_xxxxzz[k] * ab_x + g_y_0_xzzzz_xxxxxzz[k];

                g_y_0_xxzzzz_xxxyyy[k] = -g_y_0_xzzzz_xxxyyy[k] * ab_x + g_y_0_xzzzz_xxxxyyy[k];

                g_y_0_xxzzzz_xxxyyz[k] = -g_y_0_xzzzz_xxxyyz[k] * ab_x + g_y_0_xzzzz_xxxxyyz[k];

                g_y_0_xxzzzz_xxxyzz[k] = -g_y_0_xzzzz_xxxyzz[k] * ab_x + g_y_0_xzzzz_xxxxyzz[k];

                g_y_0_xxzzzz_xxxzzz[k] = -g_y_0_xzzzz_xxxzzz[k] * ab_x + g_y_0_xzzzz_xxxxzzz[k];

                g_y_0_xxzzzz_xxyyyy[k] = -g_y_0_xzzzz_xxyyyy[k] * ab_x + g_y_0_xzzzz_xxxyyyy[k];

                g_y_0_xxzzzz_xxyyyz[k] = -g_y_0_xzzzz_xxyyyz[k] * ab_x + g_y_0_xzzzz_xxxyyyz[k];

                g_y_0_xxzzzz_xxyyzz[k] = -g_y_0_xzzzz_xxyyzz[k] * ab_x + g_y_0_xzzzz_xxxyyzz[k];

                g_y_0_xxzzzz_xxyzzz[k] = -g_y_0_xzzzz_xxyzzz[k] * ab_x + g_y_0_xzzzz_xxxyzzz[k];

                g_y_0_xxzzzz_xxzzzz[k] = -g_y_0_xzzzz_xxzzzz[k] * ab_x + g_y_0_xzzzz_xxxzzzz[k];

                g_y_0_xxzzzz_xyyyyy[k] = -g_y_0_xzzzz_xyyyyy[k] * ab_x + g_y_0_xzzzz_xxyyyyy[k];

                g_y_0_xxzzzz_xyyyyz[k] = -g_y_0_xzzzz_xyyyyz[k] * ab_x + g_y_0_xzzzz_xxyyyyz[k];

                g_y_0_xxzzzz_xyyyzz[k] = -g_y_0_xzzzz_xyyyzz[k] * ab_x + g_y_0_xzzzz_xxyyyzz[k];

                g_y_0_xxzzzz_xyyzzz[k] = -g_y_0_xzzzz_xyyzzz[k] * ab_x + g_y_0_xzzzz_xxyyzzz[k];

                g_y_0_xxzzzz_xyzzzz[k] = -g_y_0_xzzzz_xyzzzz[k] * ab_x + g_y_0_xzzzz_xxyzzzz[k];

                g_y_0_xxzzzz_xzzzzz[k] = -g_y_0_xzzzz_xzzzzz[k] * ab_x + g_y_0_xzzzz_xxzzzzz[k];

                g_y_0_xxzzzz_yyyyyy[k] = -g_y_0_xzzzz_yyyyyy[k] * ab_x + g_y_0_xzzzz_xyyyyyy[k];

                g_y_0_xxzzzz_yyyyyz[k] = -g_y_0_xzzzz_yyyyyz[k] * ab_x + g_y_0_xzzzz_xyyyyyz[k];

                g_y_0_xxzzzz_yyyyzz[k] = -g_y_0_xzzzz_yyyyzz[k] * ab_x + g_y_0_xzzzz_xyyyyzz[k];

                g_y_0_xxzzzz_yyyzzz[k] = -g_y_0_xzzzz_yyyzzz[k] * ab_x + g_y_0_xzzzz_xyyyzzz[k];

                g_y_0_xxzzzz_yyzzzz[k] = -g_y_0_xzzzz_yyzzzz[k] * ab_x + g_y_0_xzzzz_xyyzzzz[k];

                g_y_0_xxzzzz_yzzzzz[k] = -g_y_0_xzzzz_yzzzzz[k] * ab_x + g_y_0_xzzzz_xyzzzzz[k];

                g_y_0_xxzzzz_zzzzzz[k] = -g_y_0_xzzzz_zzzzzz[k] * ab_x + g_y_0_xzzzz_xzzzzzz[k];
            }

            /// Set up 1204-1232 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 1204 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 1205 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 1206 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 1207 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 1208 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 1209 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 1210 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 1211 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 1212 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 1213 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 1214 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 1215 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 1216 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 1217 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 1218 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 1219 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 1220 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 1221 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 1222 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 1223 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 1224 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 1225 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 1226 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 1227 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 1228 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 1229 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 1230 * ccomps * dcomps);

            auto g_y_0_xyyyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 1231 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyyy_xxxxxx, g_y_0_xyyyyy_xxxxxy, g_y_0_xyyyyy_xxxxxz, g_y_0_xyyyyy_xxxxyy, g_y_0_xyyyyy_xxxxyz, g_y_0_xyyyyy_xxxxzz, g_y_0_xyyyyy_xxxyyy, g_y_0_xyyyyy_xxxyyz, g_y_0_xyyyyy_xxxyzz, g_y_0_xyyyyy_xxxzzz, g_y_0_xyyyyy_xxyyyy, g_y_0_xyyyyy_xxyyyz, g_y_0_xyyyyy_xxyyzz, g_y_0_xyyyyy_xxyzzz, g_y_0_xyyyyy_xxzzzz, g_y_0_xyyyyy_xyyyyy, g_y_0_xyyyyy_xyyyyz, g_y_0_xyyyyy_xyyyzz, g_y_0_xyyyyy_xyyzzz, g_y_0_xyyyyy_xyzzzz, g_y_0_xyyyyy_xzzzzz, g_y_0_xyyyyy_yyyyyy, g_y_0_xyyyyy_yyyyyz, g_y_0_xyyyyy_yyyyzz, g_y_0_xyyyyy_yyyzzz, g_y_0_xyyyyy_yyzzzz, g_y_0_xyyyyy_yzzzzz, g_y_0_xyyyyy_zzzzzz, g_y_0_yyyyy_xxxxxx, g_y_0_yyyyy_xxxxxxx, g_y_0_yyyyy_xxxxxxy, g_y_0_yyyyy_xxxxxxz, g_y_0_yyyyy_xxxxxy, g_y_0_yyyyy_xxxxxyy, g_y_0_yyyyy_xxxxxyz, g_y_0_yyyyy_xxxxxz, g_y_0_yyyyy_xxxxxzz, g_y_0_yyyyy_xxxxyy, g_y_0_yyyyy_xxxxyyy, g_y_0_yyyyy_xxxxyyz, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxyzz, g_y_0_yyyyy_xxxxzz, g_y_0_yyyyy_xxxxzzz, g_y_0_yyyyy_xxxyyy, g_y_0_yyyyy_xxxyyyy, g_y_0_yyyyy_xxxyyyz, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyyzz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxyzzz, g_y_0_yyyyy_xxxzzz, g_y_0_yyyyy_xxxzzzz, g_y_0_yyyyy_xxyyyy, g_y_0_yyyyy_xxyyyyy, g_y_0_yyyyy_xxyyyyz, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyyzz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyyzzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxyzzzz, g_y_0_yyyyy_xxzzzz, g_y_0_yyyyy_xxzzzzz, g_y_0_yyyyy_xyyyyy, g_y_0_yyyyy_xyyyyyy, g_y_0_yyyyy_xyyyyyz, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyyzz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyyzzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyyzzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xyzzzzz, g_y_0_yyyyy_xzzzzz, g_y_0_yyyyy_xzzzzzz, g_y_0_yyyyy_yyyyyy, g_y_0_yyyyy_yyyyyz, g_y_0_yyyyy_yyyyzz, g_y_0_yyyyy_yyyzzz, g_y_0_yyyyy_yyzzzz, g_y_0_yyyyy_yzzzzz, g_y_0_yyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_xxxxxx[k] = -g_y_0_yyyyy_xxxxxx[k] * ab_x + g_y_0_yyyyy_xxxxxxx[k];

                g_y_0_xyyyyy_xxxxxy[k] = -g_y_0_yyyyy_xxxxxy[k] * ab_x + g_y_0_yyyyy_xxxxxxy[k];

                g_y_0_xyyyyy_xxxxxz[k] = -g_y_0_yyyyy_xxxxxz[k] * ab_x + g_y_0_yyyyy_xxxxxxz[k];

                g_y_0_xyyyyy_xxxxyy[k] = -g_y_0_yyyyy_xxxxyy[k] * ab_x + g_y_0_yyyyy_xxxxxyy[k];

                g_y_0_xyyyyy_xxxxyz[k] = -g_y_0_yyyyy_xxxxyz[k] * ab_x + g_y_0_yyyyy_xxxxxyz[k];

                g_y_0_xyyyyy_xxxxzz[k] = -g_y_0_yyyyy_xxxxzz[k] * ab_x + g_y_0_yyyyy_xxxxxzz[k];

                g_y_0_xyyyyy_xxxyyy[k] = -g_y_0_yyyyy_xxxyyy[k] * ab_x + g_y_0_yyyyy_xxxxyyy[k];

                g_y_0_xyyyyy_xxxyyz[k] = -g_y_0_yyyyy_xxxyyz[k] * ab_x + g_y_0_yyyyy_xxxxyyz[k];

                g_y_0_xyyyyy_xxxyzz[k] = -g_y_0_yyyyy_xxxyzz[k] * ab_x + g_y_0_yyyyy_xxxxyzz[k];

                g_y_0_xyyyyy_xxxzzz[k] = -g_y_0_yyyyy_xxxzzz[k] * ab_x + g_y_0_yyyyy_xxxxzzz[k];

                g_y_0_xyyyyy_xxyyyy[k] = -g_y_0_yyyyy_xxyyyy[k] * ab_x + g_y_0_yyyyy_xxxyyyy[k];

                g_y_0_xyyyyy_xxyyyz[k] = -g_y_0_yyyyy_xxyyyz[k] * ab_x + g_y_0_yyyyy_xxxyyyz[k];

                g_y_0_xyyyyy_xxyyzz[k] = -g_y_0_yyyyy_xxyyzz[k] * ab_x + g_y_0_yyyyy_xxxyyzz[k];

                g_y_0_xyyyyy_xxyzzz[k] = -g_y_0_yyyyy_xxyzzz[k] * ab_x + g_y_0_yyyyy_xxxyzzz[k];

                g_y_0_xyyyyy_xxzzzz[k] = -g_y_0_yyyyy_xxzzzz[k] * ab_x + g_y_0_yyyyy_xxxzzzz[k];

                g_y_0_xyyyyy_xyyyyy[k] = -g_y_0_yyyyy_xyyyyy[k] * ab_x + g_y_0_yyyyy_xxyyyyy[k];

                g_y_0_xyyyyy_xyyyyz[k] = -g_y_0_yyyyy_xyyyyz[k] * ab_x + g_y_0_yyyyy_xxyyyyz[k];

                g_y_0_xyyyyy_xyyyzz[k] = -g_y_0_yyyyy_xyyyzz[k] * ab_x + g_y_0_yyyyy_xxyyyzz[k];

                g_y_0_xyyyyy_xyyzzz[k] = -g_y_0_yyyyy_xyyzzz[k] * ab_x + g_y_0_yyyyy_xxyyzzz[k];

                g_y_0_xyyyyy_xyzzzz[k] = -g_y_0_yyyyy_xyzzzz[k] * ab_x + g_y_0_yyyyy_xxyzzzz[k];

                g_y_0_xyyyyy_xzzzzz[k] = -g_y_0_yyyyy_xzzzzz[k] * ab_x + g_y_0_yyyyy_xxzzzzz[k];

                g_y_0_xyyyyy_yyyyyy[k] = -g_y_0_yyyyy_yyyyyy[k] * ab_x + g_y_0_yyyyy_xyyyyyy[k];

                g_y_0_xyyyyy_yyyyyz[k] = -g_y_0_yyyyy_yyyyyz[k] * ab_x + g_y_0_yyyyy_xyyyyyz[k];

                g_y_0_xyyyyy_yyyyzz[k] = -g_y_0_yyyyy_yyyyzz[k] * ab_x + g_y_0_yyyyy_xyyyyzz[k];

                g_y_0_xyyyyy_yyyzzz[k] = -g_y_0_yyyyy_yyyzzz[k] * ab_x + g_y_0_yyyyy_xyyyzzz[k];

                g_y_0_xyyyyy_yyzzzz[k] = -g_y_0_yyyyy_yyzzzz[k] * ab_x + g_y_0_yyyyy_xyyzzzz[k];

                g_y_0_xyyyyy_yzzzzz[k] = -g_y_0_yyyyy_yzzzzz[k] * ab_x + g_y_0_yyyyy_xyzzzzz[k];

                g_y_0_xyyyyy_zzzzzz[k] = -g_y_0_yyyyy_zzzzzz[k] * ab_x + g_y_0_yyyyy_xzzzzzz[k];
            }

            /// Set up 1232-1260 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 1232 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 1233 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 1234 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 1235 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 1236 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 1237 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 1238 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 1239 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 1240 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 1241 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 1242 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 1243 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 1244 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 1245 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 1246 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 1247 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 1248 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 1249 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 1250 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 1251 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 1252 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 1253 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 1254 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 1255 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 1256 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 1257 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 1258 * ccomps * dcomps);

            auto g_y_0_xyyyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyyz_xxxxxx, g_y_0_xyyyyz_xxxxxy, g_y_0_xyyyyz_xxxxxz, g_y_0_xyyyyz_xxxxyy, g_y_0_xyyyyz_xxxxyz, g_y_0_xyyyyz_xxxxzz, g_y_0_xyyyyz_xxxyyy, g_y_0_xyyyyz_xxxyyz, g_y_0_xyyyyz_xxxyzz, g_y_0_xyyyyz_xxxzzz, g_y_0_xyyyyz_xxyyyy, g_y_0_xyyyyz_xxyyyz, g_y_0_xyyyyz_xxyyzz, g_y_0_xyyyyz_xxyzzz, g_y_0_xyyyyz_xxzzzz, g_y_0_xyyyyz_xyyyyy, g_y_0_xyyyyz_xyyyyz, g_y_0_xyyyyz_xyyyzz, g_y_0_xyyyyz_xyyzzz, g_y_0_xyyyyz_xyzzzz, g_y_0_xyyyyz_xzzzzz, g_y_0_xyyyyz_yyyyyy, g_y_0_xyyyyz_yyyyyz, g_y_0_xyyyyz_yyyyzz, g_y_0_xyyyyz_yyyzzz, g_y_0_xyyyyz_yyzzzz, g_y_0_xyyyyz_yzzzzz, g_y_0_xyyyyz_zzzzzz, g_y_0_yyyyz_xxxxxx, g_y_0_yyyyz_xxxxxxx, g_y_0_yyyyz_xxxxxxy, g_y_0_yyyyz_xxxxxxz, g_y_0_yyyyz_xxxxxy, g_y_0_yyyyz_xxxxxyy, g_y_0_yyyyz_xxxxxyz, g_y_0_yyyyz_xxxxxz, g_y_0_yyyyz_xxxxxzz, g_y_0_yyyyz_xxxxyy, g_y_0_yyyyz_xxxxyyy, g_y_0_yyyyz_xxxxyyz, g_y_0_yyyyz_xxxxyz, g_y_0_yyyyz_xxxxyzz, g_y_0_yyyyz_xxxxzz, g_y_0_yyyyz_xxxxzzz, g_y_0_yyyyz_xxxyyy, g_y_0_yyyyz_xxxyyyy, g_y_0_yyyyz_xxxyyyz, g_y_0_yyyyz_xxxyyz, g_y_0_yyyyz_xxxyyzz, g_y_0_yyyyz_xxxyzz, g_y_0_yyyyz_xxxyzzz, g_y_0_yyyyz_xxxzzz, g_y_0_yyyyz_xxxzzzz, g_y_0_yyyyz_xxyyyy, g_y_0_yyyyz_xxyyyyy, g_y_0_yyyyz_xxyyyyz, g_y_0_yyyyz_xxyyyz, g_y_0_yyyyz_xxyyyzz, g_y_0_yyyyz_xxyyzz, g_y_0_yyyyz_xxyyzzz, g_y_0_yyyyz_xxyzzz, g_y_0_yyyyz_xxyzzzz, g_y_0_yyyyz_xxzzzz, g_y_0_yyyyz_xxzzzzz, g_y_0_yyyyz_xyyyyy, g_y_0_yyyyz_xyyyyyy, g_y_0_yyyyz_xyyyyyz, g_y_0_yyyyz_xyyyyz, g_y_0_yyyyz_xyyyyzz, g_y_0_yyyyz_xyyyzz, g_y_0_yyyyz_xyyyzzz, g_y_0_yyyyz_xyyzzz, g_y_0_yyyyz_xyyzzzz, g_y_0_yyyyz_xyzzzz, g_y_0_yyyyz_xyzzzzz, g_y_0_yyyyz_xzzzzz, g_y_0_yyyyz_xzzzzzz, g_y_0_yyyyz_yyyyyy, g_y_0_yyyyz_yyyyyz, g_y_0_yyyyz_yyyyzz, g_y_0_yyyyz_yyyzzz, g_y_0_yyyyz_yyzzzz, g_y_0_yyyyz_yzzzzz, g_y_0_yyyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_xxxxxx[k] = -g_y_0_yyyyz_xxxxxx[k] * ab_x + g_y_0_yyyyz_xxxxxxx[k];

                g_y_0_xyyyyz_xxxxxy[k] = -g_y_0_yyyyz_xxxxxy[k] * ab_x + g_y_0_yyyyz_xxxxxxy[k];

                g_y_0_xyyyyz_xxxxxz[k] = -g_y_0_yyyyz_xxxxxz[k] * ab_x + g_y_0_yyyyz_xxxxxxz[k];

                g_y_0_xyyyyz_xxxxyy[k] = -g_y_0_yyyyz_xxxxyy[k] * ab_x + g_y_0_yyyyz_xxxxxyy[k];

                g_y_0_xyyyyz_xxxxyz[k] = -g_y_0_yyyyz_xxxxyz[k] * ab_x + g_y_0_yyyyz_xxxxxyz[k];

                g_y_0_xyyyyz_xxxxzz[k] = -g_y_0_yyyyz_xxxxzz[k] * ab_x + g_y_0_yyyyz_xxxxxzz[k];

                g_y_0_xyyyyz_xxxyyy[k] = -g_y_0_yyyyz_xxxyyy[k] * ab_x + g_y_0_yyyyz_xxxxyyy[k];

                g_y_0_xyyyyz_xxxyyz[k] = -g_y_0_yyyyz_xxxyyz[k] * ab_x + g_y_0_yyyyz_xxxxyyz[k];

                g_y_0_xyyyyz_xxxyzz[k] = -g_y_0_yyyyz_xxxyzz[k] * ab_x + g_y_0_yyyyz_xxxxyzz[k];

                g_y_0_xyyyyz_xxxzzz[k] = -g_y_0_yyyyz_xxxzzz[k] * ab_x + g_y_0_yyyyz_xxxxzzz[k];

                g_y_0_xyyyyz_xxyyyy[k] = -g_y_0_yyyyz_xxyyyy[k] * ab_x + g_y_0_yyyyz_xxxyyyy[k];

                g_y_0_xyyyyz_xxyyyz[k] = -g_y_0_yyyyz_xxyyyz[k] * ab_x + g_y_0_yyyyz_xxxyyyz[k];

                g_y_0_xyyyyz_xxyyzz[k] = -g_y_0_yyyyz_xxyyzz[k] * ab_x + g_y_0_yyyyz_xxxyyzz[k];

                g_y_0_xyyyyz_xxyzzz[k] = -g_y_0_yyyyz_xxyzzz[k] * ab_x + g_y_0_yyyyz_xxxyzzz[k];

                g_y_0_xyyyyz_xxzzzz[k] = -g_y_0_yyyyz_xxzzzz[k] * ab_x + g_y_0_yyyyz_xxxzzzz[k];

                g_y_0_xyyyyz_xyyyyy[k] = -g_y_0_yyyyz_xyyyyy[k] * ab_x + g_y_0_yyyyz_xxyyyyy[k];

                g_y_0_xyyyyz_xyyyyz[k] = -g_y_0_yyyyz_xyyyyz[k] * ab_x + g_y_0_yyyyz_xxyyyyz[k];

                g_y_0_xyyyyz_xyyyzz[k] = -g_y_0_yyyyz_xyyyzz[k] * ab_x + g_y_0_yyyyz_xxyyyzz[k];

                g_y_0_xyyyyz_xyyzzz[k] = -g_y_0_yyyyz_xyyzzz[k] * ab_x + g_y_0_yyyyz_xxyyzzz[k];

                g_y_0_xyyyyz_xyzzzz[k] = -g_y_0_yyyyz_xyzzzz[k] * ab_x + g_y_0_yyyyz_xxyzzzz[k];

                g_y_0_xyyyyz_xzzzzz[k] = -g_y_0_yyyyz_xzzzzz[k] * ab_x + g_y_0_yyyyz_xxzzzzz[k];

                g_y_0_xyyyyz_yyyyyy[k] = -g_y_0_yyyyz_yyyyyy[k] * ab_x + g_y_0_yyyyz_xyyyyyy[k];

                g_y_0_xyyyyz_yyyyyz[k] = -g_y_0_yyyyz_yyyyyz[k] * ab_x + g_y_0_yyyyz_xyyyyyz[k];

                g_y_0_xyyyyz_yyyyzz[k] = -g_y_0_yyyyz_yyyyzz[k] * ab_x + g_y_0_yyyyz_xyyyyzz[k];

                g_y_0_xyyyyz_yyyzzz[k] = -g_y_0_yyyyz_yyyzzz[k] * ab_x + g_y_0_yyyyz_xyyyzzz[k];

                g_y_0_xyyyyz_yyzzzz[k] = -g_y_0_yyyyz_yyzzzz[k] * ab_x + g_y_0_yyyyz_xyyzzzz[k];

                g_y_0_xyyyyz_yzzzzz[k] = -g_y_0_yyyyz_yzzzzz[k] * ab_x + g_y_0_yyyyz_xyzzzzz[k];

                g_y_0_xyyyyz_zzzzzz[k] = -g_y_0_yyyyz_zzzzzz[k] * ab_x + g_y_0_yyyyz_xzzzzzz[k];
            }

            /// Set up 1260-1288 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1260 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1261 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1262 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1263 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1264 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1265 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1266 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1267 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1268 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1269 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1270 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1271 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1272 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1273 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1274 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1275 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1276 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1277 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1278 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1279 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1280 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1281 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1282 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1283 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1284 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1285 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1286 * ccomps * dcomps);

            auto g_y_0_xyyyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1287 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyzz_xxxxxx, g_y_0_xyyyzz_xxxxxy, g_y_0_xyyyzz_xxxxxz, g_y_0_xyyyzz_xxxxyy, g_y_0_xyyyzz_xxxxyz, g_y_0_xyyyzz_xxxxzz, g_y_0_xyyyzz_xxxyyy, g_y_0_xyyyzz_xxxyyz, g_y_0_xyyyzz_xxxyzz, g_y_0_xyyyzz_xxxzzz, g_y_0_xyyyzz_xxyyyy, g_y_0_xyyyzz_xxyyyz, g_y_0_xyyyzz_xxyyzz, g_y_0_xyyyzz_xxyzzz, g_y_0_xyyyzz_xxzzzz, g_y_0_xyyyzz_xyyyyy, g_y_0_xyyyzz_xyyyyz, g_y_0_xyyyzz_xyyyzz, g_y_0_xyyyzz_xyyzzz, g_y_0_xyyyzz_xyzzzz, g_y_0_xyyyzz_xzzzzz, g_y_0_xyyyzz_yyyyyy, g_y_0_xyyyzz_yyyyyz, g_y_0_xyyyzz_yyyyzz, g_y_0_xyyyzz_yyyzzz, g_y_0_xyyyzz_yyzzzz, g_y_0_xyyyzz_yzzzzz, g_y_0_xyyyzz_zzzzzz, g_y_0_yyyzz_xxxxxx, g_y_0_yyyzz_xxxxxxx, g_y_0_yyyzz_xxxxxxy, g_y_0_yyyzz_xxxxxxz, g_y_0_yyyzz_xxxxxy, g_y_0_yyyzz_xxxxxyy, g_y_0_yyyzz_xxxxxyz, g_y_0_yyyzz_xxxxxz, g_y_0_yyyzz_xxxxxzz, g_y_0_yyyzz_xxxxyy, g_y_0_yyyzz_xxxxyyy, g_y_0_yyyzz_xxxxyyz, g_y_0_yyyzz_xxxxyz, g_y_0_yyyzz_xxxxyzz, g_y_0_yyyzz_xxxxzz, g_y_0_yyyzz_xxxxzzz, g_y_0_yyyzz_xxxyyy, g_y_0_yyyzz_xxxyyyy, g_y_0_yyyzz_xxxyyyz, g_y_0_yyyzz_xxxyyz, g_y_0_yyyzz_xxxyyzz, g_y_0_yyyzz_xxxyzz, g_y_0_yyyzz_xxxyzzz, g_y_0_yyyzz_xxxzzz, g_y_0_yyyzz_xxxzzzz, g_y_0_yyyzz_xxyyyy, g_y_0_yyyzz_xxyyyyy, g_y_0_yyyzz_xxyyyyz, g_y_0_yyyzz_xxyyyz, g_y_0_yyyzz_xxyyyzz, g_y_0_yyyzz_xxyyzz, g_y_0_yyyzz_xxyyzzz, g_y_0_yyyzz_xxyzzz, g_y_0_yyyzz_xxyzzzz, g_y_0_yyyzz_xxzzzz, g_y_0_yyyzz_xxzzzzz, g_y_0_yyyzz_xyyyyy, g_y_0_yyyzz_xyyyyyy, g_y_0_yyyzz_xyyyyyz, g_y_0_yyyzz_xyyyyz, g_y_0_yyyzz_xyyyyzz, g_y_0_yyyzz_xyyyzz, g_y_0_yyyzz_xyyyzzz, g_y_0_yyyzz_xyyzzz, g_y_0_yyyzz_xyyzzzz, g_y_0_yyyzz_xyzzzz, g_y_0_yyyzz_xyzzzzz, g_y_0_yyyzz_xzzzzz, g_y_0_yyyzz_xzzzzzz, g_y_0_yyyzz_yyyyyy, g_y_0_yyyzz_yyyyyz, g_y_0_yyyzz_yyyyzz, g_y_0_yyyzz_yyyzzz, g_y_0_yyyzz_yyzzzz, g_y_0_yyyzz_yzzzzz, g_y_0_yyyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_xxxxxx[k] = -g_y_0_yyyzz_xxxxxx[k] * ab_x + g_y_0_yyyzz_xxxxxxx[k];

                g_y_0_xyyyzz_xxxxxy[k] = -g_y_0_yyyzz_xxxxxy[k] * ab_x + g_y_0_yyyzz_xxxxxxy[k];

                g_y_0_xyyyzz_xxxxxz[k] = -g_y_0_yyyzz_xxxxxz[k] * ab_x + g_y_0_yyyzz_xxxxxxz[k];

                g_y_0_xyyyzz_xxxxyy[k] = -g_y_0_yyyzz_xxxxyy[k] * ab_x + g_y_0_yyyzz_xxxxxyy[k];

                g_y_0_xyyyzz_xxxxyz[k] = -g_y_0_yyyzz_xxxxyz[k] * ab_x + g_y_0_yyyzz_xxxxxyz[k];

                g_y_0_xyyyzz_xxxxzz[k] = -g_y_0_yyyzz_xxxxzz[k] * ab_x + g_y_0_yyyzz_xxxxxzz[k];

                g_y_0_xyyyzz_xxxyyy[k] = -g_y_0_yyyzz_xxxyyy[k] * ab_x + g_y_0_yyyzz_xxxxyyy[k];

                g_y_0_xyyyzz_xxxyyz[k] = -g_y_0_yyyzz_xxxyyz[k] * ab_x + g_y_0_yyyzz_xxxxyyz[k];

                g_y_0_xyyyzz_xxxyzz[k] = -g_y_0_yyyzz_xxxyzz[k] * ab_x + g_y_0_yyyzz_xxxxyzz[k];

                g_y_0_xyyyzz_xxxzzz[k] = -g_y_0_yyyzz_xxxzzz[k] * ab_x + g_y_0_yyyzz_xxxxzzz[k];

                g_y_0_xyyyzz_xxyyyy[k] = -g_y_0_yyyzz_xxyyyy[k] * ab_x + g_y_0_yyyzz_xxxyyyy[k];

                g_y_0_xyyyzz_xxyyyz[k] = -g_y_0_yyyzz_xxyyyz[k] * ab_x + g_y_0_yyyzz_xxxyyyz[k];

                g_y_0_xyyyzz_xxyyzz[k] = -g_y_0_yyyzz_xxyyzz[k] * ab_x + g_y_0_yyyzz_xxxyyzz[k];

                g_y_0_xyyyzz_xxyzzz[k] = -g_y_0_yyyzz_xxyzzz[k] * ab_x + g_y_0_yyyzz_xxxyzzz[k];

                g_y_0_xyyyzz_xxzzzz[k] = -g_y_0_yyyzz_xxzzzz[k] * ab_x + g_y_0_yyyzz_xxxzzzz[k];

                g_y_0_xyyyzz_xyyyyy[k] = -g_y_0_yyyzz_xyyyyy[k] * ab_x + g_y_0_yyyzz_xxyyyyy[k];

                g_y_0_xyyyzz_xyyyyz[k] = -g_y_0_yyyzz_xyyyyz[k] * ab_x + g_y_0_yyyzz_xxyyyyz[k];

                g_y_0_xyyyzz_xyyyzz[k] = -g_y_0_yyyzz_xyyyzz[k] * ab_x + g_y_0_yyyzz_xxyyyzz[k];

                g_y_0_xyyyzz_xyyzzz[k] = -g_y_0_yyyzz_xyyzzz[k] * ab_x + g_y_0_yyyzz_xxyyzzz[k];

                g_y_0_xyyyzz_xyzzzz[k] = -g_y_0_yyyzz_xyzzzz[k] * ab_x + g_y_0_yyyzz_xxyzzzz[k];

                g_y_0_xyyyzz_xzzzzz[k] = -g_y_0_yyyzz_xzzzzz[k] * ab_x + g_y_0_yyyzz_xxzzzzz[k];

                g_y_0_xyyyzz_yyyyyy[k] = -g_y_0_yyyzz_yyyyyy[k] * ab_x + g_y_0_yyyzz_xyyyyyy[k];

                g_y_0_xyyyzz_yyyyyz[k] = -g_y_0_yyyzz_yyyyyz[k] * ab_x + g_y_0_yyyzz_xyyyyyz[k];

                g_y_0_xyyyzz_yyyyzz[k] = -g_y_0_yyyzz_yyyyzz[k] * ab_x + g_y_0_yyyzz_xyyyyzz[k];

                g_y_0_xyyyzz_yyyzzz[k] = -g_y_0_yyyzz_yyyzzz[k] * ab_x + g_y_0_yyyzz_xyyyzzz[k];

                g_y_0_xyyyzz_yyzzzz[k] = -g_y_0_yyyzz_yyzzzz[k] * ab_x + g_y_0_yyyzz_xyyzzzz[k];

                g_y_0_xyyyzz_yzzzzz[k] = -g_y_0_yyyzz_yzzzzz[k] * ab_x + g_y_0_yyyzz_xyzzzzz[k];

                g_y_0_xyyyzz_zzzzzz[k] = -g_y_0_yyyzz_zzzzzz[k] * ab_x + g_y_0_yyyzz_xzzzzzz[k];
            }

            /// Set up 1288-1316 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1288 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1289 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1290 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1291 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1292 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1293 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1294 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1295 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1296 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1297 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1298 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1299 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1300 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1301 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1302 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1303 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1304 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1305 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1306 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1307 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1308 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1309 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1310 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1311 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1312 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1313 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1314 * ccomps * dcomps);

            auto g_y_0_xyyzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1315 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyzzz_xxxxxx, g_y_0_xyyzzz_xxxxxy, g_y_0_xyyzzz_xxxxxz, g_y_0_xyyzzz_xxxxyy, g_y_0_xyyzzz_xxxxyz, g_y_0_xyyzzz_xxxxzz, g_y_0_xyyzzz_xxxyyy, g_y_0_xyyzzz_xxxyyz, g_y_0_xyyzzz_xxxyzz, g_y_0_xyyzzz_xxxzzz, g_y_0_xyyzzz_xxyyyy, g_y_0_xyyzzz_xxyyyz, g_y_0_xyyzzz_xxyyzz, g_y_0_xyyzzz_xxyzzz, g_y_0_xyyzzz_xxzzzz, g_y_0_xyyzzz_xyyyyy, g_y_0_xyyzzz_xyyyyz, g_y_0_xyyzzz_xyyyzz, g_y_0_xyyzzz_xyyzzz, g_y_0_xyyzzz_xyzzzz, g_y_0_xyyzzz_xzzzzz, g_y_0_xyyzzz_yyyyyy, g_y_0_xyyzzz_yyyyyz, g_y_0_xyyzzz_yyyyzz, g_y_0_xyyzzz_yyyzzz, g_y_0_xyyzzz_yyzzzz, g_y_0_xyyzzz_yzzzzz, g_y_0_xyyzzz_zzzzzz, g_y_0_yyzzz_xxxxxx, g_y_0_yyzzz_xxxxxxx, g_y_0_yyzzz_xxxxxxy, g_y_0_yyzzz_xxxxxxz, g_y_0_yyzzz_xxxxxy, g_y_0_yyzzz_xxxxxyy, g_y_0_yyzzz_xxxxxyz, g_y_0_yyzzz_xxxxxz, g_y_0_yyzzz_xxxxxzz, g_y_0_yyzzz_xxxxyy, g_y_0_yyzzz_xxxxyyy, g_y_0_yyzzz_xxxxyyz, g_y_0_yyzzz_xxxxyz, g_y_0_yyzzz_xxxxyzz, g_y_0_yyzzz_xxxxzz, g_y_0_yyzzz_xxxxzzz, g_y_0_yyzzz_xxxyyy, g_y_0_yyzzz_xxxyyyy, g_y_0_yyzzz_xxxyyyz, g_y_0_yyzzz_xxxyyz, g_y_0_yyzzz_xxxyyzz, g_y_0_yyzzz_xxxyzz, g_y_0_yyzzz_xxxyzzz, g_y_0_yyzzz_xxxzzz, g_y_0_yyzzz_xxxzzzz, g_y_0_yyzzz_xxyyyy, g_y_0_yyzzz_xxyyyyy, g_y_0_yyzzz_xxyyyyz, g_y_0_yyzzz_xxyyyz, g_y_0_yyzzz_xxyyyzz, g_y_0_yyzzz_xxyyzz, g_y_0_yyzzz_xxyyzzz, g_y_0_yyzzz_xxyzzz, g_y_0_yyzzz_xxyzzzz, g_y_0_yyzzz_xxzzzz, g_y_0_yyzzz_xxzzzzz, g_y_0_yyzzz_xyyyyy, g_y_0_yyzzz_xyyyyyy, g_y_0_yyzzz_xyyyyyz, g_y_0_yyzzz_xyyyyz, g_y_0_yyzzz_xyyyyzz, g_y_0_yyzzz_xyyyzz, g_y_0_yyzzz_xyyyzzz, g_y_0_yyzzz_xyyzzz, g_y_0_yyzzz_xyyzzzz, g_y_0_yyzzz_xyzzzz, g_y_0_yyzzz_xyzzzzz, g_y_0_yyzzz_xzzzzz, g_y_0_yyzzz_xzzzzzz, g_y_0_yyzzz_yyyyyy, g_y_0_yyzzz_yyyyyz, g_y_0_yyzzz_yyyyzz, g_y_0_yyzzz_yyyzzz, g_y_0_yyzzz_yyzzzz, g_y_0_yyzzz_yzzzzz, g_y_0_yyzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_xxxxxx[k] = -g_y_0_yyzzz_xxxxxx[k] * ab_x + g_y_0_yyzzz_xxxxxxx[k];

                g_y_0_xyyzzz_xxxxxy[k] = -g_y_0_yyzzz_xxxxxy[k] * ab_x + g_y_0_yyzzz_xxxxxxy[k];

                g_y_0_xyyzzz_xxxxxz[k] = -g_y_0_yyzzz_xxxxxz[k] * ab_x + g_y_0_yyzzz_xxxxxxz[k];

                g_y_0_xyyzzz_xxxxyy[k] = -g_y_0_yyzzz_xxxxyy[k] * ab_x + g_y_0_yyzzz_xxxxxyy[k];

                g_y_0_xyyzzz_xxxxyz[k] = -g_y_0_yyzzz_xxxxyz[k] * ab_x + g_y_0_yyzzz_xxxxxyz[k];

                g_y_0_xyyzzz_xxxxzz[k] = -g_y_0_yyzzz_xxxxzz[k] * ab_x + g_y_0_yyzzz_xxxxxzz[k];

                g_y_0_xyyzzz_xxxyyy[k] = -g_y_0_yyzzz_xxxyyy[k] * ab_x + g_y_0_yyzzz_xxxxyyy[k];

                g_y_0_xyyzzz_xxxyyz[k] = -g_y_0_yyzzz_xxxyyz[k] * ab_x + g_y_0_yyzzz_xxxxyyz[k];

                g_y_0_xyyzzz_xxxyzz[k] = -g_y_0_yyzzz_xxxyzz[k] * ab_x + g_y_0_yyzzz_xxxxyzz[k];

                g_y_0_xyyzzz_xxxzzz[k] = -g_y_0_yyzzz_xxxzzz[k] * ab_x + g_y_0_yyzzz_xxxxzzz[k];

                g_y_0_xyyzzz_xxyyyy[k] = -g_y_0_yyzzz_xxyyyy[k] * ab_x + g_y_0_yyzzz_xxxyyyy[k];

                g_y_0_xyyzzz_xxyyyz[k] = -g_y_0_yyzzz_xxyyyz[k] * ab_x + g_y_0_yyzzz_xxxyyyz[k];

                g_y_0_xyyzzz_xxyyzz[k] = -g_y_0_yyzzz_xxyyzz[k] * ab_x + g_y_0_yyzzz_xxxyyzz[k];

                g_y_0_xyyzzz_xxyzzz[k] = -g_y_0_yyzzz_xxyzzz[k] * ab_x + g_y_0_yyzzz_xxxyzzz[k];

                g_y_0_xyyzzz_xxzzzz[k] = -g_y_0_yyzzz_xxzzzz[k] * ab_x + g_y_0_yyzzz_xxxzzzz[k];

                g_y_0_xyyzzz_xyyyyy[k] = -g_y_0_yyzzz_xyyyyy[k] * ab_x + g_y_0_yyzzz_xxyyyyy[k];

                g_y_0_xyyzzz_xyyyyz[k] = -g_y_0_yyzzz_xyyyyz[k] * ab_x + g_y_0_yyzzz_xxyyyyz[k];

                g_y_0_xyyzzz_xyyyzz[k] = -g_y_0_yyzzz_xyyyzz[k] * ab_x + g_y_0_yyzzz_xxyyyzz[k];

                g_y_0_xyyzzz_xyyzzz[k] = -g_y_0_yyzzz_xyyzzz[k] * ab_x + g_y_0_yyzzz_xxyyzzz[k];

                g_y_0_xyyzzz_xyzzzz[k] = -g_y_0_yyzzz_xyzzzz[k] * ab_x + g_y_0_yyzzz_xxyzzzz[k];

                g_y_0_xyyzzz_xzzzzz[k] = -g_y_0_yyzzz_xzzzzz[k] * ab_x + g_y_0_yyzzz_xxzzzzz[k];

                g_y_0_xyyzzz_yyyyyy[k] = -g_y_0_yyzzz_yyyyyy[k] * ab_x + g_y_0_yyzzz_xyyyyyy[k];

                g_y_0_xyyzzz_yyyyyz[k] = -g_y_0_yyzzz_yyyyyz[k] * ab_x + g_y_0_yyzzz_xyyyyyz[k];

                g_y_0_xyyzzz_yyyyzz[k] = -g_y_0_yyzzz_yyyyzz[k] * ab_x + g_y_0_yyzzz_xyyyyzz[k];

                g_y_0_xyyzzz_yyyzzz[k] = -g_y_0_yyzzz_yyyzzz[k] * ab_x + g_y_0_yyzzz_xyyyzzz[k];

                g_y_0_xyyzzz_yyzzzz[k] = -g_y_0_yyzzz_yyzzzz[k] * ab_x + g_y_0_yyzzz_xyyzzzz[k];

                g_y_0_xyyzzz_yzzzzz[k] = -g_y_0_yyzzz_yzzzzz[k] * ab_x + g_y_0_yyzzz_xyzzzzz[k];

                g_y_0_xyyzzz_zzzzzz[k] = -g_y_0_yyzzz_zzzzzz[k] * ab_x + g_y_0_yyzzz_xzzzzzz[k];
            }

            /// Set up 1316-1344 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1316 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1317 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1318 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1319 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1320 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1321 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1322 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1323 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1324 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1325 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1326 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1327 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1328 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1329 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1330 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1331 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1332 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1333 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1334 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1335 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1336 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1337 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1338 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1339 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1340 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1341 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1342 * ccomps * dcomps);

            auto g_y_0_xyzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1343 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzzzz_xxxxxx, g_y_0_xyzzzz_xxxxxy, g_y_0_xyzzzz_xxxxxz, g_y_0_xyzzzz_xxxxyy, g_y_0_xyzzzz_xxxxyz, g_y_0_xyzzzz_xxxxzz, g_y_0_xyzzzz_xxxyyy, g_y_0_xyzzzz_xxxyyz, g_y_0_xyzzzz_xxxyzz, g_y_0_xyzzzz_xxxzzz, g_y_0_xyzzzz_xxyyyy, g_y_0_xyzzzz_xxyyyz, g_y_0_xyzzzz_xxyyzz, g_y_0_xyzzzz_xxyzzz, g_y_0_xyzzzz_xxzzzz, g_y_0_xyzzzz_xyyyyy, g_y_0_xyzzzz_xyyyyz, g_y_0_xyzzzz_xyyyzz, g_y_0_xyzzzz_xyyzzz, g_y_0_xyzzzz_xyzzzz, g_y_0_xyzzzz_xzzzzz, g_y_0_xyzzzz_yyyyyy, g_y_0_xyzzzz_yyyyyz, g_y_0_xyzzzz_yyyyzz, g_y_0_xyzzzz_yyyzzz, g_y_0_xyzzzz_yyzzzz, g_y_0_xyzzzz_yzzzzz, g_y_0_xyzzzz_zzzzzz, g_y_0_yzzzz_xxxxxx, g_y_0_yzzzz_xxxxxxx, g_y_0_yzzzz_xxxxxxy, g_y_0_yzzzz_xxxxxxz, g_y_0_yzzzz_xxxxxy, g_y_0_yzzzz_xxxxxyy, g_y_0_yzzzz_xxxxxyz, g_y_0_yzzzz_xxxxxz, g_y_0_yzzzz_xxxxxzz, g_y_0_yzzzz_xxxxyy, g_y_0_yzzzz_xxxxyyy, g_y_0_yzzzz_xxxxyyz, g_y_0_yzzzz_xxxxyz, g_y_0_yzzzz_xxxxyzz, g_y_0_yzzzz_xxxxzz, g_y_0_yzzzz_xxxxzzz, g_y_0_yzzzz_xxxyyy, g_y_0_yzzzz_xxxyyyy, g_y_0_yzzzz_xxxyyyz, g_y_0_yzzzz_xxxyyz, g_y_0_yzzzz_xxxyyzz, g_y_0_yzzzz_xxxyzz, g_y_0_yzzzz_xxxyzzz, g_y_0_yzzzz_xxxzzz, g_y_0_yzzzz_xxxzzzz, g_y_0_yzzzz_xxyyyy, g_y_0_yzzzz_xxyyyyy, g_y_0_yzzzz_xxyyyyz, g_y_0_yzzzz_xxyyyz, g_y_0_yzzzz_xxyyyzz, g_y_0_yzzzz_xxyyzz, g_y_0_yzzzz_xxyyzzz, g_y_0_yzzzz_xxyzzz, g_y_0_yzzzz_xxyzzzz, g_y_0_yzzzz_xxzzzz, g_y_0_yzzzz_xxzzzzz, g_y_0_yzzzz_xyyyyy, g_y_0_yzzzz_xyyyyyy, g_y_0_yzzzz_xyyyyyz, g_y_0_yzzzz_xyyyyz, g_y_0_yzzzz_xyyyyzz, g_y_0_yzzzz_xyyyzz, g_y_0_yzzzz_xyyyzzz, g_y_0_yzzzz_xyyzzz, g_y_0_yzzzz_xyyzzzz, g_y_0_yzzzz_xyzzzz, g_y_0_yzzzz_xyzzzzz, g_y_0_yzzzz_xzzzzz, g_y_0_yzzzz_xzzzzzz, g_y_0_yzzzz_yyyyyy, g_y_0_yzzzz_yyyyyz, g_y_0_yzzzz_yyyyzz, g_y_0_yzzzz_yyyzzz, g_y_0_yzzzz_yyzzzz, g_y_0_yzzzz_yzzzzz, g_y_0_yzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_xxxxxx[k] = -g_y_0_yzzzz_xxxxxx[k] * ab_x + g_y_0_yzzzz_xxxxxxx[k];

                g_y_0_xyzzzz_xxxxxy[k] = -g_y_0_yzzzz_xxxxxy[k] * ab_x + g_y_0_yzzzz_xxxxxxy[k];

                g_y_0_xyzzzz_xxxxxz[k] = -g_y_0_yzzzz_xxxxxz[k] * ab_x + g_y_0_yzzzz_xxxxxxz[k];

                g_y_0_xyzzzz_xxxxyy[k] = -g_y_0_yzzzz_xxxxyy[k] * ab_x + g_y_0_yzzzz_xxxxxyy[k];

                g_y_0_xyzzzz_xxxxyz[k] = -g_y_0_yzzzz_xxxxyz[k] * ab_x + g_y_0_yzzzz_xxxxxyz[k];

                g_y_0_xyzzzz_xxxxzz[k] = -g_y_0_yzzzz_xxxxzz[k] * ab_x + g_y_0_yzzzz_xxxxxzz[k];

                g_y_0_xyzzzz_xxxyyy[k] = -g_y_0_yzzzz_xxxyyy[k] * ab_x + g_y_0_yzzzz_xxxxyyy[k];

                g_y_0_xyzzzz_xxxyyz[k] = -g_y_0_yzzzz_xxxyyz[k] * ab_x + g_y_0_yzzzz_xxxxyyz[k];

                g_y_0_xyzzzz_xxxyzz[k] = -g_y_0_yzzzz_xxxyzz[k] * ab_x + g_y_0_yzzzz_xxxxyzz[k];

                g_y_0_xyzzzz_xxxzzz[k] = -g_y_0_yzzzz_xxxzzz[k] * ab_x + g_y_0_yzzzz_xxxxzzz[k];

                g_y_0_xyzzzz_xxyyyy[k] = -g_y_0_yzzzz_xxyyyy[k] * ab_x + g_y_0_yzzzz_xxxyyyy[k];

                g_y_0_xyzzzz_xxyyyz[k] = -g_y_0_yzzzz_xxyyyz[k] * ab_x + g_y_0_yzzzz_xxxyyyz[k];

                g_y_0_xyzzzz_xxyyzz[k] = -g_y_0_yzzzz_xxyyzz[k] * ab_x + g_y_0_yzzzz_xxxyyzz[k];

                g_y_0_xyzzzz_xxyzzz[k] = -g_y_0_yzzzz_xxyzzz[k] * ab_x + g_y_0_yzzzz_xxxyzzz[k];

                g_y_0_xyzzzz_xxzzzz[k] = -g_y_0_yzzzz_xxzzzz[k] * ab_x + g_y_0_yzzzz_xxxzzzz[k];

                g_y_0_xyzzzz_xyyyyy[k] = -g_y_0_yzzzz_xyyyyy[k] * ab_x + g_y_0_yzzzz_xxyyyyy[k];

                g_y_0_xyzzzz_xyyyyz[k] = -g_y_0_yzzzz_xyyyyz[k] * ab_x + g_y_0_yzzzz_xxyyyyz[k];

                g_y_0_xyzzzz_xyyyzz[k] = -g_y_0_yzzzz_xyyyzz[k] * ab_x + g_y_0_yzzzz_xxyyyzz[k];

                g_y_0_xyzzzz_xyyzzz[k] = -g_y_0_yzzzz_xyyzzz[k] * ab_x + g_y_0_yzzzz_xxyyzzz[k];

                g_y_0_xyzzzz_xyzzzz[k] = -g_y_0_yzzzz_xyzzzz[k] * ab_x + g_y_0_yzzzz_xxyzzzz[k];

                g_y_0_xyzzzz_xzzzzz[k] = -g_y_0_yzzzz_xzzzzz[k] * ab_x + g_y_0_yzzzz_xxzzzzz[k];

                g_y_0_xyzzzz_yyyyyy[k] = -g_y_0_yzzzz_yyyyyy[k] * ab_x + g_y_0_yzzzz_xyyyyyy[k];

                g_y_0_xyzzzz_yyyyyz[k] = -g_y_0_yzzzz_yyyyyz[k] * ab_x + g_y_0_yzzzz_xyyyyyz[k];

                g_y_0_xyzzzz_yyyyzz[k] = -g_y_0_yzzzz_yyyyzz[k] * ab_x + g_y_0_yzzzz_xyyyyzz[k];

                g_y_0_xyzzzz_yyyzzz[k] = -g_y_0_yzzzz_yyyzzz[k] * ab_x + g_y_0_yzzzz_xyyyzzz[k];

                g_y_0_xyzzzz_yyzzzz[k] = -g_y_0_yzzzz_yyzzzz[k] * ab_x + g_y_0_yzzzz_xyyzzzz[k];

                g_y_0_xyzzzz_yzzzzz[k] = -g_y_0_yzzzz_yzzzzz[k] * ab_x + g_y_0_yzzzz_xyzzzzz[k];

                g_y_0_xyzzzz_zzzzzz[k] = -g_y_0_yzzzz_zzzzzz[k] * ab_x + g_y_0_yzzzz_xzzzzzz[k];
            }

            /// Set up 1344-1372 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1344 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1345 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1346 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1347 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1348 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1349 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1350 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1351 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1352 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1353 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1354 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1355 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1356 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1357 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1358 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1359 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1360 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1361 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1362 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1363 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1364 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1365 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1366 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1367 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1368 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1369 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1370 * ccomps * dcomps);

            auto g_y_0_xzzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1371 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzzzz_xxxxxx, g_y_0_xzzzzz_xxxxxy, g_y_0_xzzzzz_xxxxxz, g_y_0_xzzzzz_xxxxyy, g_y_0_xzzzzz_xxxxyz, g_y_0_xzzzzz_xxxxzz, g_y_0_xzzzzz_xxxyyy, g_y_0_xzzzzz_xxxyyz, g_y_0_xzzzzz_xxxyzz, g_y_0_xzzzzz_xxxzzz, g_y_0_xzzzzz_xxyyyy, g_y_0_xzzzzz_xxyyyz, g_y_0_xzzzzz_xxyyzz, g_y_0_xzzzzz_xxyzzz, g_y_0_xzzzzz_xxzzzz, g_y_0_xzzzzz_xyyyyy, g_y_0_xzzzzz_xyyyyz, g_y_0_xzzzzz_xyyyzz, g_y_0_xzzzzz_xyyzzz, g_y_0_xzzzzz_xyzzzz, g_y_0_xzzzzz_xzzzzz, g_y_0_xzzzzz_yyyyyy, g_y_0_xzzzzz_yyyyyz, g_y_0_xzzzzz_yyyyzz, g_y_0_xzzzzz_yyyzzz, g_y_0_xzzzzz_yyzzzz, g_y_0_xzzzzz_yzzzzz, g_y_0_xzzzzz_zzzzzz, g_y_0_zzzzz_xxxxxx, g_y_0_zzzzz_xxxxxxx, g_y_0_zzzzz_xxxxxxy, g_y_0_zzzzz_xxxxxxz, g_y_0_zzzzz_xxxxxy, g_y_0_zzzzz_xxxxxyy, g_y_0_zzzzz_xxxxxyz, g_y_0_zzzzz_xxxxxz, g_y_0_zzzzz_xxxxxzz, g_y_0_zzzzz_xxxxyy, g_y_0_zzzzz_xxxxyyy, g_y_0_zzzzz_xxxxyyz, g_y_0_zzzzz_xxxxyz, g_y_0_zzzzz_xxxxyzz, g_y_0_zzzzz_xxxxzz, g_y_0_zzzzz_xxxxzzz, g_y_0_zzzzz_xxxyyy, g_y_0_zzzzz_xxxyyyy, g_y_0_zzzzz_xxxyyyz, g_y_0_zzzzz_xxxyyz, g_y_0_zzzzz_xxxyyzz, g_y_0_zzzzz_xxxyzz, g_y_0_zzzzz_xxxyzzz, g_y_0_zzzzz_xxxzzz, g_y_0_zzzzz_xxxzzzz, g_y_0_zzzzz_xxyyyy, g_y_0_zzzzz_xxyyyyy, g_y_0_zzzzz_xxyyyyz, g_y_0_zzzzz_xxyyyz, g_y_0_zzzzz_xxyyyzz, g_y_0_zzzzz_xxyyzz, g_y_0_zzzzz_xxyyzzz, g_y_0_zzzzz_xxyzzz, g_y_0_zzzzz_xxyzzzz, g_y_0_zzzzz_xxzzzz, g_y_0_zzzzz_xxzzzzz, g_y_0_zzzzz_xyyyyy, g_y_0_zzzzz_xyyyyyy, g_y_0_zzzzz_xyyyyyz, g_y_0_zzzzz_xyyyyz, g_y_0_zzzzz_xyyyyzz, g_y_0_zzzzz_xyyyzz, g_y_0_zzzzz_xyyyzzz, g_y_0_zzzzz_xyyzzz, g_y_0_zzzzz_xyyzzzz, g_y_0_zzzzz_xyzzzz, g_y_0_zzzzz_xyzzzzz, g_y_0_zzzzz_xzzzzz, g_y_0_zzzzz_xzzzzzz, g_y_0_zzzzz_yyyyyy, g_y_0_zzzzz_yyyyyz, g_y_0_zzzzz_yyyyzz, g_y_0_zzzzz_yyyzzz, g_y_0_zzzzz_yyzzzz, g_y_0_zzzzz_yzzzzz, g_y_0_zzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_xxxxxx[k] = -g_y_0_zzzzz_xxxxxx[k] * ab_x + g_y_0_zzzzz_xxxxxxx[k];

                g_y_0_xzzzzz_xxxxxy[k] = -g_y_0_zzzzz_xxxxxy[k] * ab_x + g_y_0_zzzzz_xxxxxxy[k];

                g_y_0_xzzzzz_xxxxxz[k] = -g_y_0_zzzzz_xxxxxz[k] * ab_x + g_y_0_zzzzz_xxxxxxz[k];

                g_y_0_xzzzzz_xxxxyy[k] = -g_y_0_zzzzz_xxxxyy[k] * ab_x + g_y_0_zzzzz_xxxxxyy[k];

                g_y_0_xzzzzz_xxxxyz[k] = -g_y_0_zzzzz_xxxxyz[k] * ab_x + g_y_0_zzzzz_xxxxxyz[k];

                g_y_0_xzzzzz_xxxxzz[k] = -g_y_0_zzzzz_xxxxzz[k] * ab_x + g_y_0_zzzzz_xxxxxzz[k];

                g_y_0_xzzzzz_xxxyyy[k] = -g_y_0_zzzzz_xxxyyy[k] * ab_x + g_y_0_zzzzz_xxxxyyy[k];

                g_y_0_xzzzzz_xxxyyz[k] = -g_y_0_zzzzz_xxxyyz[k] * ab_x + g_y_0_zzzzz_xxxxyyz[k];

                g_y_0_xzzzzz_xxxyzz[k] = -g_y_0_zzzzz_xxxyzz[k] * ab_x + g_y_0_zzzzz_xxxxyzz[k];

                g_y_0_xzzzzz_xxxzzz[k] = -g_y_0_zzzzz_xxxzzz[k] * ab_x + g_y_0_zzzzz_xxxxzzz[k];

                g_y_0_xzzzzz_xxyyyy[k] = -g_y_0_zzzzz_xxyyyy[k] * ab_x + g_y_0_zzzzz_xxxyyyy[k];

                g_y_0_xzzzzz_xxyyyz[k] = -g_y_0_zzzzz_xxyyyz[k] * ab_x + g_y_0_zzzzz_xxxyyyz[k];

                g_y_0_xzzzzz_xxyyzz[k] = -g_y_0_zzzzz_xxyyzz[k] * ab_x + g_y_0_zzzzz_xxxyyzz[k];

                g_y_0_xzzzzz_xxyzzz[k] = -g_y_0_zzzzz_xxyzzz[k] * ab_x + g_y_0_zzzzz_xxxyzzz[k];

                g_y_0_xzzzzz_xxzzzz[k] = -g_y_0_zzzzz_xxzzzz[k] * ab_x + g_y_0_zzzzz_xxxzzzz[k];

                g_y_0_xzzzzz_xyyyyy[k] = -g_y_0_zzzzz_xyyyyy[k] * ab_x + g_y_0_zzzzz_xxyyyyy[k];

                g_y_0_xzzzzz_xyyyyz[k] = -g_y_0_zzzzz_xyyyyz[k] * ab_x + g_y_0_zzzzz_xxyyyyz[k];

                g_y_0_xzzzzz_xyyyzz[k] = -g_y_0_zzzzz_xyyyzz[k] * ab_x + g_y_0_zzzzz_xxyyyzz[k];

                g_y_0_xzzzzz_xyyzzz[k] = -g_y_0_zzzzz_xyyzzz[k] * ab_x + g_y_0_zzzzz_xxyyzzz[k];

                g_y_0_xzzzzz_xyzzzz[k] = -g_y_0_zzzzz_xyzzzz[k] * ab_x + g_y_0_zzzzz_xxyzzzz[k];

                g_y_0_xzzzzz_xzzzzz[k] = -g_y_0_zzzzz_xzzzzz[k] * ab_x + g_y_0_zzzzz_xxzzzzz[k];

                g_y_0_xzzzzz_yyyyyy[k] = -g_y_0_zzzzz_yyyyyy[k] * ab_x + g_y_0_zzzzz_xyyyyyy[k];

                g_y_0_xzzzzz_yyyyyz[k] = -g_y_0_zzzzz_yyyyyz[k] * ab_x + g_y_0_zzzzz_xyyyyyz[k];

                g_y_0_xzzzzz_yyyyzz[k] = -g_y_0_zzzzz_yyyyzz[k] * ab_x + g_y_0_zzzzz_xyyyyzz[k];

                g_y_0_xzzzzz_yyyzzz[k] = -g_y_0_zzzzz_yyyzzz[k] * ab_x + g_y_0_zzzzz_xyyyzzz[k];

                g_y_0_xzzzzz_yyzzzz[k] = -g_y_0_zzzzz_yyzzzz[k] * ab_x + g_y_0_zzzzz_xyyzzzz[k];

                g_y_0_xzzzzz_yzzzzz[k] = -g_y_0_zzzzz_yzzzzz[k] * ab_x + g_y_0_zzzzz_xyzzzzz[k];

                g_y_0_xzzzzz_zzzzzz[k] = -g_y_0_zzzzz_zzzzzz[k] * ab_x + g_y_0_zzzzz_xzzzzzz[k];
            }

            /// Set up 1372-1400 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 1372 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 1373 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 1374 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 1375 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 1376 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 1377 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 1378 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 1379 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 1380 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 1381 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 1382 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 1383 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 1384 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 1385 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 1386 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 1387 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 1388 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 1389 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 1390 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 1391 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 1392 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 1393 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 1394 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 1395 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 1396 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 1397 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 1398 * ccomps * dcomps);

            auto g_y_0_yyyyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 1399 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_xxxxxx, g_y_0_yyyyy_xxxxxxy, g_y_0_yyyyy_xxxxxy, g_y_0_yyyyy_xxxxxyy, g_y_0_yyyyy_xxxxxyz, g_y_0_yyyyy_xxxxxz, g_y_0_yyyyy_xxxxyy, g_y_0_yyyyy_xxxxyyy, g_y_0_yyyyy_xxxxyyz, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxyzz, g_y_0_yyyyy_xxxxzz, g_y_0_yyyyy_xxxyyy, g_y_0_yyyyy_xxxyyyy, g_y_0_yyyyy_xxxyyyz, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyyzz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxyzzz, g_y_0_yyyyy_xxxzzz, g_y_0_yyyyy_xxyyyy, g_y_0_yyyyy_xxyyyyy, g_y_0_yyyyy_xxyyyyz, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyyzz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyyzzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxyzzzz, g_y_0_yyyyy_xxzzzz, g_y_0_yyyyy_xyyyyy, g_y_0_yyyyy_xyyyyyy, g_y_0_yyyyy_xyyyyyz, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyyzz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyyzzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyyzzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xyzzzzz, g_y_0_yyyyy_xzzzzz, g_y_0_yyyyy_yyyyyy, g_y_0_yyyyy_yyyyyyy, g_y_0_yyyyy_yyyyyyz, g_y_0_yyyyy_yyyyyz, g_y_0_yyyyy_yyyyyzz, g_y_0_yyyyy_yyyyzz, g_y_0_yyyyy_yyyyzzz, g_y_0_yyyyy_yyyzzz, g_y_0_yyyyy_yyyzzzz, g_y_0_yyyyy_yyzzzz, g_y_0_yyyyy_yyzzzzz, g_y_0_yyyyy_yzzzzz, g_y_0_yyyyy_yzzzzzz, g_y_0_yyyyy_zzzzzz, g_y_0_yyyyyy_xxxxxx, g_y_0_yyyyyy_xxxxxy, g_y_0_yyyyyy_xxxxxz, g_y_0_yyyyyy_xxxxyy, g_y_0_yyyyyy_xxxxyz, g_y_0_yyyyyy_xxxxzz, g_y_0_yyyyyy_xxxyyy, g_y_0_yyyyyy_xxxyyz, g_y_0_yyyyyy_xxxyzz, g_y_0_yyyyyy_xxxzzz, g_y_0_yyyyyy_xxyyyy, g_y_0_yyyyyy_xxyyyz, g_y_0_yyyyyy_xxyyzz, g_y_0_yyyyyy_xxyzzz, g_y_0_yyyyyy_xxzzzz, g_y_0_yyyyyy_xyyyyy, g_y_0_yyyyyy_xyyyyz, g_y_0_yyyyyy_xyyyzz, g_y_0_yyyyyy_xyyzzz, g_y_0_yyyyyy_xyzzzz, g_y_0_yyyyyy_xzzzzz, g_y_0_yyyyyy_yyyyyy, g_y_0_yyyyyy_yyyyyz, g_y_0_yyyyyy_yyyyzz, g_y_0_yyyyyy_yyyzzz, g_y_0_yyyyyy_yyzzzz, g_y_0_yyyyyy_yzzzzz, g_y_0_yyyyyy_zzzzzz, g_yyyyy_xxxxxx, g_yyyyy_xxxxxy, g_yyyyy_xxxxxz, g_yyyyy_xxxxyy, g_yyyyy_xxxxyz, g_yyyyy_xxxxzz, g_yyyyy_xxxyyy, g_yyyyy_xxxyyz, g_yyyyy_xxxyzz, g_yyyyy_xxxzzz, g_yyyyy_xxyyyy, g_yyyyy_xxyyyz, g_yyyyy_xxyyzz, g_yyyyy_xxyzzz, g_yyyyy_xxzzzz, g_yyyyy_xyyyyy, g_yyyyy_xyyyyz, g_yyyyy_xyyyzz, g_yyyyy_xyyzzz, g_yyyyy_xyzzzz, g_yyyyy_xzzzzz, g_yyyyy_yyyyyy, g_yyyyy_yyyyyz, g_yyyyy_yyyyzz, g_yyyyy_yyyzzz, g_yyyyy_yyzzzz, g_yyyyy_yzzzzz, g_yyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_xxxxxx[k] = -g_yyyyy_xxxxxx[k] - g_y_0_yyyyy_xxxxxx[k] * ab_y + g_y_0_yyyyy_xxxxxxy[k];

                g_y_0_yyyyyy_xxxxxy[k] = -g_yyyyy_xxxxxy[k] - g_y_0_yyyyy_xxxxxy[k] * ab_y + g_y_0_yyyyy_xxxxxyy[k];

                g_y_0_yyyyyy_xxxxxz[k] = -g_yyyyy_xxxxxz[k] - g_y_0_yyyyy_xxxxxz[k] * ab_y + g_y_0_yyyyy_xxxxxyz[k];

                g_y_0_yyyyyy_xxxxyy[k] = -g_yyyyy_xxxxyy[k] - g_y_0_yyyyy_xxxxyy[k] * ab_y + g_y_0_yyyyy_xxxxyyy[k];

                g_y_0_yyyyyy_xxxxyz[k] = -g_yyyyy_xxxxyz[k] - g_y_0_yyyyy_xxxxyz[k] * ab_y + g_y_0_yyyyy_xxxxyyz[k];

                g_y_0_yyyyyy_xxxxzz[k] = -g_yyyyy_xxxxzz[k] - g_y_0_yyyyy_xxxxzz[k] * ab_y + g_y_0_yyyyy_xxxxyzz[k];

                g_y_0_yyyyyy_xxxyyy[k] = -g_yyyyy_xxxyyy[k] - g_y_0_yyyyy_xxxyyy[k] * ab_y + g_y_0_yyyyy_xxxyyyy[k];

                g_y_0_yyyyyy_xxxyyz[k] = -g_yyyyy_xxxyyz[k] - g_y_0_yyyyy_xxxyyz[k] * ab_y + g_y_0_yyyyy_xxxyyyz[k];

                g_y_0_yyyyyy_xxxyzz[k] = -g_yyyyy_xxxyzz[k] - g_y_0_yyyyy_xxxyzz[k] * ab_y + g_y_0_yyyyy_xxxyyzz[k];

                g_y_0_yyyyyy_xxxzzz[k] = -g_yyyyy_xxxzzz[k] - g_y_0_yyyyy_xxxzzz[k] * ab_y + g_y_0_yyyyy_xxxyzzz[k];

                g_y_0_yyyyyy_xxyyyy[k] = -g_yyyyy_xxyyyy[k] - g_y_0_yyyyy_xxyyyy[k] * ab_y + g_y_0_yyyyy_xxyyyyy[k];

                g_y_0_yyyyyy_xxyyyz[k] = -g_yyyyy_xxyyyz[k] - g_y_0_yyyyy_xxyyyz[k] * ab_y + g_y_0_yyyyy_xxyyyyz[k];

                g_y_0_yyyyyy_xxyyzz[k] = -g_yyyyy_xxyyzz[k] - g_y_0_yyyyy_xxyyzz[k] * ab_y + g_y_0_yyyyy_xxyyyzz[k];

                g_y_0_yyyyyy_xxyzzz[k] = -g_yyyyy_xxyzzz[k] - g_y_0_yyyyy_xxyzzz[k] * ab_y + g_y_0_yyyyy_xxyyzzz[k];

                g_y_0_yyyyyy_xxzzzz[k] = -g_yyyyy_xxzzzz[k] - g_y_0_yyyyy_xxzzzz[k] * ab_y + g_y_0_yyyyy_xxyzzzz[k];

                g_y_0_yyyyyy_xyyyyy[k] = -g_yyyyy_xyyyyy[k] - g_y_0_yyyyy_xyyyyy[k] * ab_y + g_y_0_yyyyy_xyyyyyy[k];

                g_y_0_yyyyyy_xyyyyz[k] = -g_yyyyy_xyyyyz[k] - g_y_0_yyyyy_xyyyyz[k] * ab_y + g_y_0_yyyyy_xyyyyyz[k];

                g_y_0_yyyyyy_xyyyzz[k] = -g_yyyyy_xyyyzz[k] - g_y_0_yyyyy_xyyyzz[k] * ab_y + g_y_0_yyyyy_xyyyyzz[k];

                g_y_0_yyyyyy_xyyzzz[k] = -g_yyyyy_xyyzzz[k] - g_y_0_yyyyy_xyyzzz[k] * ab_y + g_y_0_yyyyy_xyyyzzz[k];

                g_y_0_yyyyyy_xyzzzz[k] = -g_yyyyy_xyzzzz[k] - g_y_0_yyyyy_xyzzzz[k] * ab_y + g_y_0_yyyyy_xyyzzzz[k];

                g_y_0_yyyyyy_xzzzzz[k] = -g_yyyyy_xzzzzz[k] - g_y_0_yyyyy_xzzzzz[k] * ab_y + g_y_0_yyyyy_xyzzzzz[k];

                g_y_0_yyyyyy_yyyyyy[k] = -g_yyyyy_yyyyyy[k] - g_y_0_yyyyy_yyyyyy[k] * ab_y + g_y_0_yyyyy_yyyyyyy[k];

                g_y_0_yyyyyy_yyyyyz[k] = -g_yyyyy_yyyyyz[k] - g_y_0_yyyyy_yyyyyz[k] * ab_y + g_y_0_yyyyy_yyyyyyz[k];

                g_y_0_yyyyyy_yyyyzz[k] = -g_yyyyy_yyyyzz[k] - g_y_0_yyyyy_yyyyzz[k] * ab_y + g_y_0_yyyyy_yyyyyzz[k];

                g_y_0_yyyyyy_yyyzzz[k] = -g_yyyyy_yyyzzz[k] - g_y_0_yyyyy_yyyzzz[k] * ab_y + g_y_0_yyyyy_yyyyzzz[k];

                g_y_0_yyyyyy_yyzzzz[k] = -g_yyyyy_yyzzzz[k] - g_y_0_yyyyy_yyzzzz[k] * ab_y + g_y_0_yyyyy_yyyzzzz[k];

                g_y_0_yyyyyy_yzzzzz[k] = -g_yyyyy_yzzzzz[k] - g_y_0_yyyyy_yzzzzz[k] * ab_y + g_y_0_yyyyy_yyzzzzz[k];

                g_y_0_yyyyyy_zzzzzz[k] = -g_yyyyy_zzzzzz[k] - g_y_0_yyyyy_zzzzzz[k] * ab_y + g_y_0_yyyyy_yzzzzzz[k];
            }

            /// Set up 1400-1428 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 1400 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 1401 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 1402 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 1403 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 1404 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 1405 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 1406 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 1407 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 1408 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 1409 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 1410 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 1411 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 1412 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 1413 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 1414 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 1415 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 1416 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 1417 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 1418 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 1419 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 1420 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 1421 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 1422 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 1423 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 1424 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 1425 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 1426 * ccomps * dcomps);

            auto g_y_0_yyyyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 1427 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_xxxxxx, g_y_0_yyyyy_xxxxxxz, g_y_0_yyyyy_xxxxxy, g_y_0_yyyyy_xxxxxyz, g_y_0_yyyyy_xxxxxz, g_y_0_yyyyy_xxxxxzz, g_y_0_yyyyy_xxxxyy, g_y_0_yyyyy_xxxxyyz, g_y_0_yyyyy_xxxxyz, g_y_0_yyyyy_xxxxyzz, g_y_0_yyyyy_xxxxzz, g_y_0_yyyyy_xxxxzzz, g_y_0_yyyyy_xxxyyy, g_y_0_yyyyy_xxxyyyz, g_y_0_yyyyy_xxxyyz, g_y_0_yyyyy_xxxyyzz, g_y_0_yyyyy_xxxyzz, g_y_0_yyyyy_xxxyzzz, g_y_0_yyyyy_xxxzzz, g_y_0_yyyyy_xxxzzzz, g_y_0_yyyyy_xxyyyy, g_y_0_yyyyy_xxyyyyz, g_y_0_yyyyy_xxyyyz, g_y_0_yyyyy_xxyyyzz, g_y_0_yyyyy_xxyyzz, g_y_0_yyyyy_xxyyzzz, g_y_0_yyyyy_xxyzzz, g_y_0_yyyyy_xxyzzzz, g_y_0_yyyyy_xxzzzz, g_y_0_yyyyy_xxzzzzz, g_y_0_yyyyy_xyyyyy, g_y_0_yyyyy_xyyyyyz, g_y_0_yyyyy_xyyyyz, g_y_0_yyyyy_xyyyyzz, g_y_0_yyyyy_xyyyzz, g_y_0_yyyyy_xyyyzzz, g_y_0_yyyyy_xyyzzz, g_y_0_yyyyy_xyyzzzz, g_y_0_yyyyy_xyzzzz, g_y_0_yyyyy_xyzzzzz, g_y_0_yyyyy_xzzzzz, g_y_0_yyyyy_xzzzzzz, g_y_0_yyyyy_yyyyyy, g_y_0_yyyyy_yyyyyyz, g_y_0_yyyyy_yyyyyz, g_y_0_yyyyy_yyyyyzz, g_y_0_yyyyy_yyyyzz, g_y_0_yyyyy_yyyyzzz, g_y_0_yyyyy_yyyzzz, g_y_0_yyyyy_yyyzzzz, g_y_0_yyyyy_yyzzzz, g_y_0_yyyyy_yyzzzzz, g_y_0_yyyyy_yzzzzz, g_y_0_yyyyy_yzzzzzz, g_y_0_yyyyy_zzzzzz, g_y_0_yyyyy_zzzzzzz, g_y_0_yyyyyz_xxxxxx, g_y_0_yyyyyz_xxxxxy, g_y_0_yyyyyz_xxxxxz, g_y_0_yyyyyz_xxxxyy, g_y_0_yyyyyz_xxxxyz, g_y_0_yyyyyz_xxxxzz, g_y_0_yyyyyz_xxxyyy, g_y_0_yyyyyz_xxxyyz, g_y_0_yyyyyz_xxxyzz, g_y_0_yyyyyz_xxxzzz, g_y_0_yyyyyz_xxyyyy, g_y_0_yyyyyz_xxyyyz, g_y_0_yyyyyz_xxyyzz, g_y_0_yyyyyz_xxyzzz, g_y_0_yyyyyz_xxzzzz, g_y_0_yyyyyz_xyyyyy, g_y_0_yyyyyz_xyyyyz, g_y_0_yyyyyz_xyyyzz, g_y_0_yyyyyz_xyyzzz, g_y_0_yyyyyz_xyzzzz, g_y_0_yyyyyz_xzzzzz, g_y_0_yyyyyz_yyyyyy, g_y_0_yyyyyz_yyyyyz, g_y_0_yyyyyz_yyyyzz, g_y_0_yyyyyz_yyyzzz, g_y_0_yyyyyz_yyzzzz, g_y_0_yyyyyz_yzzzzz, g_y_0_yyyyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_xxxxxx[k] = -g_y_0_yyyyy_xxxxxx[k] * ab_z + g_y_0_yyyyy_xxxxxxz[k];

                g_y_0_yyyyyz_xxxxxy[k] = -g_y_0_yyyyy_xxxxxy[k] * ab_z + g_y_0_yyyyy_xxxxxyz[k];

                g_y_0_yyyyyz_xxxxxz[k] = -g_y_0_yyyyy_xxxxxz[k] * ab_z + g_y_0_yyyyy_xxxxxzz[k];

                g_y_0_yyyyyz_xxxxyy[k] = -g_y_0_yyyyy_xxxxyy[k] * ab_z + g_y_0_yyyyy_xxxxyyz[k];

                g_y_0_yyyyyz_xxxxyz[k] = -g_y_0_yyyyy_xxxxyz[k] * ab_z + g_y_0_yyyyy_xxxxyzz[k];

                g_y_0_yyyyyz_xxxxzz[k] = -g_y_0_yyyyy_xxxxzz[k] * ab_z + g_y_0_yyyyy_xxxxzzz[k];

                g_y_0_yyyyyz_xxxyyy[k] = -g_y_0_yyyyy_xxxyyy[k] * ab_z + g_y_0_yyyyy_xxxyyyz[k];

                g_y_0_yyyyyz_xxxyyz[k] = -g_y_0_yyyyy_xxxyyz[k] * ab_z + g_y_0_yyyyy_xxxyyzz[k];

                g_y_0_yyyyyz_xxxyzz[k] = -g_y_0_yyyyy_xxxyzz[k] * ab_z + g_y_0_yyyyy_xxxyzzz[k];

                g_y_0_yyyyyz_xxxzzz[k] = -g_y_0_yyyyy_xxxzzz[k] * ab_z + g_y_0_yyyyy_xxxzzzz[k];

                g_y_0_yyyyyz_xxyyyy[k] = -g_y_0_yyyyy_xxyyyy[k] * ab_z + g_y_0_yyyyy_xxyyyyz[k];

                g_y_0_yyyyyz_xxyyyz[k] = -g_y_0_yyyyy_xxyyyz[k] * ab_z + g_y_0_yyyyy_xxyyyzz[k];

                g_y_0_yyyyyz_xxyyzz[k] = -g_y_0_yyyyy_xxyyzz[k] * ab_z + g_y_0_yyyyy_xxyyzzz[k];

                g_y_0_yyyyyz_xxyzzz[k] = -g_y_0_yyyyy_xxyzzz[k] * ab_z + g_y_0_yyyyy_xxyzzzz[k];

                g_y_0_yyyyyz_xxzzzz[k] = -g_y_0_yyyyy_xxzzzz[k] * ab_z + g_y_0_yyyyy_xxzzzzz[k];

                g_y_0_yyyyyz_xyyyyy[k] = -g_y_0_yyyyy_xyyyyy[k] * ab_z + g_y_0_yyyyy_xyyyyyz[k];

                g_y_0_yyyyyz_xyyyyz[k] = -g_y_0_yyyyy_xyyyyz[k] * ab_z + g_y_0_yyyyy_xyyyyzz[k];

                g_y_0_yyyyyz_xyyyzz[k] = -g_y_0_yyyyy_xyyyzz[k] * ab_z + g_y_0_yyyyy_xyyyzzz[k];

                g_y_0_yyyyyz_xyyzzz[k] = -g_y_0_yyyyy_xyyzzz[k] * ab_z + g_y_0_yyyyy_xyyzzzz[k];

                g_y_0_yyyyyz_xyzzzz[k] = -g_y_0_yyyyy_xyzzzz[k] * ab_z + g_y_0_yyyyy_xyzzzzz[k];

                g_y_0_yyyyyz_xzzzzz[k] = -g_y_0_yyyyy_xzzzzz[k] * ab_z + g_y_0_yyyyy_xzzzzzz[k];

                g_y_0_yyyyyz_yyyyyy[k] = -g_y_0_yyyyy_yyyyyy[k] * ab_z + g_y_0_yyyyy_yyyyyyz[k];

                g_y_0_yyyyyz_yyyyyz[k] = -g_y_0_yyyyy_yyyyyz[k] * ab_z + g_y_0_yyyyy_yyyyyzz[k];

                g_y_0_yyyyyz_yyyyzz[k] = -g_y_0_yyyyy_yyyyzz[k] * ab_z + g_y_0_yyyyy_yyyyzzz[k];

                g_y_0_yyyyyz_yyyzzz[k] = -g_y_0_yyyyy_yyyzzz[k] * ab_z + g_y_0_yyyyy_yyyzzzz[k];

                g_y_0_yyyyyz_yyzzzz[k] = -g_y_0_yyyyy_yyzzzz[k] * ab_z + g_y_0_yyyyy_yyzzzzz[k];

                g_y_0_yyyyyz_yzzzzz[k] = -g_y_0_yyyyy_yzzzzz[k] * ab_z + g_y_0_yyyyy_yzzzzzz[k];

                g_y_0_yyyyyz_zzzzzz[k] = -g_y_0_yyyyy_zzzzzz[k] * ab_z + g_y_0_yyyyy_zzzzzzz[k];
            }

            /// Set up 1428-1456 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1428 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1429 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1430 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1431 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1432 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1433 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1434 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1435 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1436 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1437 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1438 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1439 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1440 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1441 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1442 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1443 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1444 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1445 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1446 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1447 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1448 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1449 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1450 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1451 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1452 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1453 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1454 * ccomps * dcomps);

            auto g_y_0_yyyyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1455 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyz_xxxxxx, g_y_0_yyyyz_xxxxxxz, g_y_0_yyyyz_xxxxxy, g_y_0_yyyyz_xxxxxyz, g_y_0_yyyyz_xxxxxz, g_y_0_yyyyz_xxxxxzz, g_y_0_yyyyz_xxxxyy, g_y_0_yyyyz_xxxxyyz, g_y_0_yyyyz_xxxxyz, g_y_0_yyyyz_xxxxyzz, g_y_0_yyyyz_xxxxzz, g_y_0_yyyyz_xxxxzzz, g_y_0_yyyyz_xxxyyy, g_y_0_yyyyz_xxxyyyz, g_y_0_yyyyz_xxxyyz, g_y_0_yyyyz_xxxyyzz, g_y_0_yyyyz_xxxyzz, g_y_0_yyyyz_xxxyzzz, g_y_0_yyyyz_xxxzzz, g_y_0_yyyyz_xxxzzzz, g_y_0_yyyyz_xxyyyy, g_y_0_yyyyz_xxyyyyz, g_y_0_yyyyz_xxyyyz, g_y_0_yyyyz_xxyyyzz, g_y_0_yyyyz_xxyyzz, g_y_0_yyyyz_xxyyzzz, g_y_0_yyyyz_xxyzzz, g_y_0_yyyyz_xxyzzzz, g_y_0_yyyyz_xxzzzz, g_y_0_yyyyz_xxzzzzz, g_y_0_yyyyz_xyyyyy, g_y_0_yyyyz_xyyyyyz, g_y_0_yyyyz_xyyyyz, g_y_0_yyyyz_xyyyyzz, g_y_0_yyyyz_xyyyzz, g_y_0_yyyyz_xyyyzzz, g_y_0_yyyyz_xyyzzz, g_y_0_yyyyz_xyyzzzz, g_y_0_yyyyz_xyzzzz, g_y_0_yyyyz_xyzzzzz, g_y_0_yyyyz_xzzzzz, g_y_0_yyyyz_xzzzzzz, g_y_0_yyyyz_yyyyyy, g_y_0_yyyyz_yyyyyyz, g_y_0_yyyyz_yyyyyz, g_y_0_yyyyz_yyyyyzz, g_y_0_yyyyz_yyyyzz, g_y_0_yyyyz_yyyyzzz, g_y_0_yyyyz_yyyzzz, g_y_0_yyyyz_yyyzzzz, g_y_0_yyyyz_yyzzzz, g_y_0_yyyyz_yyzzzzz, g_y_0_yyyyz_yzzzzz, g_y_0_yyyyz_yzzzzzz, g_y_0_yyyyz_zzzzzz, g_y_0_yyyyz_zzzzzzz, g_y_0_yyyyzz_xxxxxx, g_y_0_yyyyzz_xxxxxy, g_y_0_yyyyzz_xxxxxz, g_y_0_yyyyzz_xxxxyy, g_y_0_yyyyzz_xxxxyz, g_y_0_yyyyzz_xxxxzz, g_y_0_yyyyzz_xxxyyy, g_y_0_yyyyzz_xxxyyz, g_y_0_yyyyzz_xxxyzz, g_y_0_yyyyzz_xxxzzz, g_y_0_yyyyzz_xxyyyy, g_y_0_yyyyzz_xxyyyz, g_y_0_yyyyzz_xxyyzz, g_y_0_yyyyzz_xxyzzz, g_y_0_yyyyzz_xxzzzz, g_y_0_yyyyzz_xyyyyy, g_y_0_yyyyzz_xyyyyz, g_y_0_yyyyzz_xyyyzz, g_y_0_yyyyzz_xyyzzz, g_y_0_yyyyzz_xyzzzz, g_y_0_yyyyzz_xzzzzz, g_y_0_yyyyzz_yyyyyy, g_y_0_yyyyzz_yyyyyz, g_y_0_yyyyzz_yyyyzz, g_y_0_yyyyzz_yyyzzz, g_y_0_yyyyzz_yyzzzz, g_y_0_yyyyzz_yzzzzz, g_y_0_yyyyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_xxxxxx[k] = -g_y_0_yyyyz_xxxxxx[k] * ab_z + g_y_0_yyyyz_xxxxxxz[k];

                g_y_0_yyyyzz_xxxxxy[k] = -g_y_0_yyyyz_xxxxxy[k] * ab_z + g_y_0_yyyyz_xxxxxyz[k];

                g_y_0_yyyyzz_xxxxxz[k] = -g_y_0_yyyyz_xxxxxz[k] * ab_z + g_y_0_yyyyz_xxxxxzz[k];

                g_y_0_yyyyzz_xxxxyy[k] = -g_y_0_yyyyz_xxxxyy[k] * ab_z + g_y_0_yyyyz_xxxxyyz[k];

                g_y_0_yyyyzz_xxxxyz[k] = -g_y_0_yyyyz_xxxxyz[k] * ab_z + g_y_0_yyyyz_xxxxyzz[k];

                g_y_0_yyyyzz_xxxxzz[k] = -g_y_0_yyyyz_xxxxzz[k] * ab_z + g_y_0_yyyyz_xxxxzzz[k];

                g_y_0_yyyyzz_xxxyyy[k] = -g_y_0_yyyyz_xxxyyy[k] * ab_z + g_y_0_yyyyz_xxxyyyz[k];

                g_y_0_yyyyzz_xxxyyz[k] = -g_y_0_yyyyz_xxxyyz[k] * ab_z + g_y_0_yyyyz_xxxyyzz[k];

                g_y_0_yyyyzz_xxxyzz[k] = -g_y_0_yyyyz_xxxyzz[k] * ab_z + g_y_0_yyyyz_xxxyzzz[k];

                g_y_0_yyyyzz_xxxzzz[k] = -g_y_0_yyyyz_xxxzzz[k] * ab_z + g_y_0_yyyyz_xxxzzzz[k];

                g_y_0_yyyyzz_xxyyyy[k] = -g_y_0_yyyyz_xxyyyy[k] * ab_z + g_y_0_yyyyz_xxyyyyz[k];

                g_y_0_yyyyzz_xxyyyz[k] = -g_y_0_yyyyz_xxyyyz[k] * ab_z + g_y_0_yyyyz_xxyyyzz[k];

                g_y_0_yyyyzz_xxyyzz[k] = -g_y_0_yyyyz_xxyyzz[k] * ab_z + g_y_0_yyyyz_xxyyzzz[k];

                g_y_0_yyyyzz_xxyzzz[k] = -g_y_0_yyyyz_xxyzzz[k] * ab_z + g_y_0_yyyyz_xxyzzzz[k];

                g_y_0_yyyyzz_xxzzzz[k] = -g_y_0_yyyyz_xxzzzz[k] * ab_z + g_y_0_yyyyz_xxzzzzz[k];

                g_y_0_yyyyzz_xyyyyy[k] = -g_y_0_yyyyz_xyyyyy[k] * ab_z + g_y_0_yyyyz_xyyyyyz[k];

                g_y_0_yyyyzz_xyyyyz[k] = -g_y_0_yyyyz_xyyyyz[k] * ab_z + g_y_0_yyyyz_xyyyyzz[k];

                g_y_0_yyyyzz_xyyyzz[k] = -g_y_0_yyyyz_xyyyzz[k] * ab_z + g_y_0_yyyyz_xyyyzzz[k];

                g_y_0_yyyyzz_xyyzzz[k] = -g_y_0_yyyyz_xyyzzz[k] * ab_z + g_y_0_yyyyz_xyyzzzz[k];

                g_y_0_yyyyzz_xyzzzz[k] = -g_y_0_yyyyz_xyzzzz[k] * ab_z + g_y_0_yyyyz_xyzzzzz[k];

                g_y_0_yyyyzz_xzzzzz[k] = -g_y_0_yyyyz_xzzzzz[k] * ab_z + g_y_0_yyyyz_xzzzzzz[k];

                g_y_0_yyyyzz_yyyyyy[k] = -g_y_0_yyyyz_yyyyyy[k] * ab_z + g_y_0_yyyyz_yyyyyyz[k];

                g_y_0_yyyyzz_yyyyyz[k] = -g_y_0_yyyyz_yyyyyz[k] * ab_z + g_y_0_yyyyz_yyyyyzz[k];

                g_y_0_yyyyzz_yyyyzz[k] = -g_y_0_yyyyz_yyyyzz[k] * ab_z + g_y_0_yyyyz_yyyyzzz[k];

                g_y_0_yyyyzz_yyyzzz[k] = -g_y_0_yyyyz_yyyzzz[k] * ab_z + g_y_0_yyyyz_yyyzzzz[k];

                g_y_0_yyyyzz_yyzzzz[k] = -g_y_0_yyyyz_yyzzzz[k] * ab_z + g_y_0_yyyyz_yyzzzzz[k];

                g_y_0_yyyyzz_yzzzzz[k] = -g_y_0_yyyyz_yzzzzz[k] * ab_z + g_y_0_yyyyz_yzzzzzz[k];

                g_y_0_yyyyzz_zzzzzz[k] = -g_y_0_yyyyz_zzzzzz[k] * ab_z + g_y_0_yyyyz_zzzzzzz[k];
            }

            /// Set up 1456-1484 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1456 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1457 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1458 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1459 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1460 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1461 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1462 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1463 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1464 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1465 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1466 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1467 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1468 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1469 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1470 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1471 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1472 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1473 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1474 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1475 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1476 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1477 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1478 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1479 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1480 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1481 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1482 * ccomps * dcomps);

            auto g_y_0_yyyzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1483 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyzz_xxxxxx, g_y_0_yyyzz_xxxxxxz, g_y_0_yyyzz_xxxxxy, g_y_0_yyyzz_xxxxxyz, g_y_0_yyyzz_xxxxxz, g_y_0_yyyzz_xxxxxzz, g_y_0_yyyzz_xxxxyy, g_y_0_yyyzz_xxxxyyz, g_y_0_yyyzz_xxxxyz, g_y_0_yyyzz_xxxxyzz, g_y_0_yyyzz_xxxxzz, g_y_0_yyyzz_xxxxzzz, g_y_0_yyyzz_xxxyyy, g_y_0_yyyzz_xxxyyyz, g_y_0_yyyzz_xxxyyz, g_y_0_yyyzz_xxxyyzz, g_y_0_yyyzz_xxxyzz, g_y_0_yyyzz_xxxyzzz, g_y_0_yyyzz_xxxzzz, g_y_0_yyyzz_xxxzzzz, g_y_0_yyyzz_xxyyyy, g_y_0_yyyzz_xxyyyyz, g_y_0_yyyzz_xxyyyz, g_y_0_yyyzz_xxyyyzz, g_y_0_yyyzz_xxyyzz, g_y_0_yyyzz_xxyyzzz, g_y_0_yyyzz_xxyzzz, g_y_0_yyyzz_xxyzzzz, g_y_0_yyyzz_xxzzzz, g_y_0_yyyzz_xxzzzzz, g_y_0_yyyzz_xyyyyy, g_y_0_yyyzz_xyyyyyz, g_y_0_yyyzz_xyyyyz, g_y_0_yyyzz_xyyyyzz, g_y_0_yyyzz_xyyyzz, g_y_0_yyyzz_xyyyzzz, g_y_0_yyyzz_xyyzzz, g_y_0_yyyzz_xyyzzzz, g_y_0_yyyzz_xyzzzz, g_y_0_yyyzz_xyzzzzz, g_y_0_yyyzz_xzzzzz, g_y_0_yyyzz_xzzzzzz, g_y_0_yyyzz_yyyyyy, g_y_0_yyyzz_yyyyyyz, g_y_0_yyyzz_yyyyyz, g_y_0_yyyzz_yyyyyzz, g_y_0_yyyzz_yyyyzz, g_y_0_yyyzz_yyyyzzz, g_y_0_yyyzz_yyyzzz, g_y_0_yyyzz_yyyzzzz, g_y_0_yyyzz_yyzzzz, g_y_0_yyyzz_yyzzzzz, g_y_0_yyyzz_yzzzzz, g_y_0_yyyzz_yzzzzzz, g_y_0_yyyzz_zzzzzz, g_y_0_yyyzz_zzzzzzz, g_y_0_yyyzzz_xxxxxx, g_y_0_yyyzzz_xxxxxy, g_y_0_yyyzzz_xxxxxz, g_y_0_yyyzzz_xxxxyy, g_y_0_yyyzzz_xxxxyz, g_y_0_yyyzzz_xxxxzz, g_y_0_yyyzzz_xxxyyy, g_y_0_yyyzzz_xxxyyz, g_y_0_yyyzzz_xxxyzz, g_y_0_yyyzzz_xxxzzz, g_y_0_yyyzzz_xxyyyy, g_y_0_yyyzzz_xxyyyz, g_y_0_yyyzzz_xxyyzz, g_y_0_yyyzzz_xxyzzz, g_y_0_yyyzzz_xxzzzz, g_y_0_yyyzzz_xyyyyy, g_y_0_yyyzzz_xyyyyz, g_y_0_yyyzzz_xyyyzz, g_y_0_yyyzzz_xyyzzz, g_y_0_yyyzzz_xyzzzz, g_y_0_yyyzzz_xzzzzz, g_y_0_yyyzzz_yyyyyy, g_y_0_yyyzzz_yyyyyz, g_y_0_yyyzzz_yyyyzz, g_y_0_yyyzzz_yyyzzz, g_y_0_yyyzzz_yyzzzz, g_y_0_yyyzzz_yzzzzz, g_y_0_yyyzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_xxxxxx[k] = -g_y_0_yyyzz_xxxxxx[k] * ab_z + g_y_0_yyyzz_xxxxxxz[k];

                g_y_0_yyyzzz_xxxxxy[k] = -g_y_0_yyyzz_xxxxxy[k] * ab_z + g_y_0_yyyzz_xxxxxyz[k];

                g_y_0_yyyzzz_xxxxxz[k] = -g_y_0_yyyzz_xxxxxz[k] * ab_z + g_y_0_yyyzz_xxxxxzz[k];

                g_y_0_yyyzzz_xxxxyy[k] = -g_y_0_yyyzz_xxxxyy[k] * ab_z + g_y_0_yyyzz_xxxxyyz[k];

                g_y_0_yyyzzz_xxxxyz[k] = -g_y_0_yyyzz_xxxxyz[k] * ab_z + g_y_0_yyyzz_xxxxyzz[k];

                g_y_0_yyyzzz_xxxxzz[k] = -g_y_0_yyyzz_xxxxzz[k] * ab_z + g_y_0_yyyzz_xxxxzzz[k];

                g_y_0_yyyzzz_xxxyyy[k] = -g_y_0_yyyzz_xxxyyy[k] * ab_z + g_y_0_yyyzz_xxxyyyz[k];

                g_y_0_yyyzzz_xxxyyz[k] = -g_y_0_yyyzz_xxxyyz[k] * ab_z + g_y_0_yyyzz_xxxyyzz[k];

                g_y_0_yyyzzz_xxxyzz[k] = -g_y_0_yyyzz_xxxyzz[k] * ab_z + g_y_0_yyyzz_xxxyzzz[k];

                g_y_0_yyyzzz_xxxzzz[k] = -g_y_0_yyyzz_xxxzzz[k] * ab_z + g_y_0_yyyzz_xxxzzzz[k];

                g_y_0_yyyzzz_xxyyyy[k] = -g_y_0_yyyzz_xxyyyy[k] * ab_z + g_y_0_yyyzz_xxyyyyz[k];

                g_y_0_yyyzzz_xxyyyz[k] = -g_y_0_yyyzz_xxyyyz[k] * ab_z + g_y_0_yyyzz_xxyyyzz[k];

                g_y_0_yyyzzz_xxyyzz[k] = -g_y_0_yyyzz_xxyyzz[k] * ab_z + g_y_0_yyyzz_xxyyzzz[k];

                g_y_0_yyyzzz_xxyzzz[k] = -g_y_0_yyyzz_xxyzzz[k] * ab_z + g_y_0_yyyzz_xxyzzzz[k];

                g_y_0_yyyzzz_xxzzzz[k] = -g_y_0_yyyzz_xxzzzz[k] * ab_z + g_y_0_yyyzz_xxzzzzz[k];

                g_y_0_yyyzzz_xyyyyy[k] = -g_y_0_yyyzz_xyyyyy[k] * ab_z + g_y_0_yyyzz_xyyyyyz[k];

                g_y_0_yyyzzz_xyyyyz[k] = -g_y_0_yyyzz_xyyyyz[k] * ab_z + g_y_0_yyyzz_xyyyyzz[k];

                g_y_0_yyyzzz_xyyyzz[k] = -g_y_0_yyyzz_xyyyzz[k] * ab_z + g_y_0_yyyzz_xyyyzzz[k];

                g_y_0_yyyzzz_xyyzzz[k] = -g_y_0_yyyzz_xyyzzz[k] * ab_z + g_y_0_yyyzz_xyyzzzz[k];

                g_y_0_yyyzzz_xyzzzz[k] = -g_y_0_yyyzz_xyzzzz[k] * ab_z + g_y_0_yyyzz_xyzzzzz[k];

                g_y_0_yyyzzz_xzzzzz[k] = -g_y_0_yyyzz_xzzzzz[k] * ab_z + g_y_0_yyyzz_xzzzzzz[k];

                g_y_0_yyyzzz_yyyyyy[k] = -g_y_0_yyyzz_yyyyyy[k] * ab_z + g_y_0_yyyzz_yyyyyyz[k];

                g_y_0_yyyzzz_yyyyyz[k] = -g_y_0_yyyzz_yyyyyz[k] * ab_z + g_y_0_yyyzz_yyyyyzz[k];

                g_y_0_yyyzzz_yyyyzz[k] = -g_y_0_yyyzz_yyyyzz[k] * ab_z + g_y_0_yyyzz_yyyyzzz[k];

                g_y_0_yyyzzz_yyyzzz[k] = -g_y_0_yyyzz_yyyzzz[k] * ab_z + g_y_0_yyyzz_yyyzzzz[k];

                g_y_0_yyyzzz_yyzzzz[k] = -g_y_0_yyyzz_yyzzzz[k] * ab_z + g_y_0_yyyzz_yyzzzzz[k];

                g_y_0_yyyzzz_yzzzzz[k] = -g_y_0_yyyzz_yzzzzz[k] * ab_z + g_y_0_yyyzz_yzzzzzz[k];

                g_y_0_yyyzzz_zzzzzz[k] = -g_y_0_yyyzz_zzzzzz[k] * ab_z + g_y_0_yyyzz_zzzzzzz[k];
            }

            /// Set up 1484-1512 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1484 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1485 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1486 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1487 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1488 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1489 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1490 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1491 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1492 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1493 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1494 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1495 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1496 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1497 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1498 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1499 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1500 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1501 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1502 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1503 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1504 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1505 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1506 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1507 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1508 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1509 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1510 * ccomps * dcomps);

            auto g_y_0_yyzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1511 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyzzz_xxxxxx, g_y_0_yyzzz_xxxxxxz, g_y_0_yyzzz_xxxxxy, g_y_0_yyzzz_xxxxxyz, g_y_0_yyzzz_xxxxxz, g_y_0_yyzzz_xxxxxzz, g_y_0_yyzzz_xxxxyy, g_y_0_yyzzz_xxxxyyz, g_y_0_yyzzz_xxxxyz, g_y_0_yyzzz_xxxxyzz, g_y_0_yyzzz_xxxxzz, g_y_0_yyzzz_xxxxzzz, g_y_0_yyzzz_xxxyyy, g_y_0_yyzzz_xxxyyyz, g_y_0_yyzzz_xxxyyz, g_y_0_yyzzz_xxxyyzz, g_y_0_yyzzz_xxxyzz, g_y_0_yyzzz_xxxyzzz, g_y_0_yyzzz_xxxzzz, g_y_0_yyzzz_xxxzzzz, g_y_0_yyzzz_xxyyyy, g_y_0_yyzzz_xxyyyyz, g_y_0_yyzzz_xxyyyz, g_y_0_yyzzz_xxyyyzz, g_y_0_yyzzz_xxyyzz, g_y_0_yyzzz_xxyyzzz, g_y_0_yyzzz_xxyzzz, g_y_0_yyzzz_xxyzzzz, g_y_0_yyzzz_xxzzzz, g_y_0_yyzzz_xxzzzzz, g_y_0_yyzzz_xyyyyy, g_y_0_yyzzz_xyyyyyz, g_y_0_yyzzz_xyyyyz, g_y_0_yyzzz_xyyyyzz, g_y_0_yyzzz_xyyyzz, g_y_0_yyzzz_xyyyzzz, g_y_0_yyzzz_xyyzzz, g_y_0_yyzzz_xyyzzzz, g_y_0_yyzzz_xyzzzz, g_y_0_yyzzz_xyzzzzz, g_y_0_yyzzz_xzzzzz, g_y_0_yyzzz_xzzzzzz, g_y_0_yyzzz_yyyyyy, g_y_0_yyzzz_yyyyyyz, g_y_0_yyzzz_yyyyyz, g_y_0_yyzzz_yyyyyzz, g_y_0_yyzzz_yyyyzz, g_y_0_yyzzz_yyyyzzz, g_y_0_yyzzz_yyyzzz, g_y_0_yyzzz_yyyzzzz, g_y_0_yyzzz_yyzzzz, g_y_0_yyzzz_yyzzzzz, g_y_0_yyzzz_yzzzzz, g_y_0_yyzzz_yzzzzzz, g_y_0_yyzzz_zzzzzz, g_y_0_yyzzz_zzzzzzz, g_y_0_yyzzzz_xxxxxx, g_y_0_yyzzzz_xxxxxy, g_y_0_yyzzzz_xxxxxz, g_y_0_yyzzzz_xxxxyy, g_y_0_yyzzzz_xxxxyz, g_y_0_yyzzzz_xxxxzz, g_y_0_yyzzzz_xxxyyy, g_y_0_yyzzzz_xxxyyz, g_y_0_yyzzzz_xxxyzz, g_y_0_yyzzzz_xxxzzz, g_y_0_yyzzzz_xxyyyy, g_y_0_yyzzzz_xxyyyz, g_y_0_yyzzzz_xxyyzz, g_y_0_yyzzzz_xxyzzz, g_y_0_yyzzzz_xxzzzz, g_y_0_yyzzzz_xyyyyy, g_y_0_yyzzzz_xyyyyz, g_y_0_yyzzzz_xyyyzz, g_y_0_yyzzzz_xyyzzz, g_y_0_yyzzzz_xyzzzz, g_y_0_yyzzzz_xzzzzz, g_y_0_yyzzzz_yyyyyy, g_y_0_yyzzzz_yyyyyz, g_y_0_yyzzzz_yyyyzz, g_y_0_yyzzzz_yyyzzz, g_y_0_yyzzzz_yyzzzz, g_y_0_yyzzzz_yzzzzz, g_y_0_yyzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_xxxxxx[k] = -g_y_0_yyzzz_xxxxxx[k] * ab_z + g_y_0_yyzzz_xxxxxxz[k];

                g_y_0_yyzzzz_xxxxxy[k] = -g_y_0_yyzzz_xxxxxy[k] * ab_z + g_y_0_yyzzz_xxxxxyz[k];

                g_y_0_yyzzzz_xxxxxz[k] = -g_y_0_yyzzz_xxxxxz[k] * ab_z + g_y_0_yyzzz_xxxxxzz[k];

                g_y_0_yyzzzz_xxxxyy[k] = -g_y_0_yyzzz_xxxxyy[k] * ab_z + g_y_0_yyzzz_xxxxyyz[k];

                g_y_0_yyzzzz_xxxxyz[k] = -g_y_0_yyzzz_xxxxyz[k] * ab_z + g_y_0_yyzzz_xxxxyzz[k];

                g_y_0_yyzzzz_xxxxzz[k] = -g_y_0_yyzzz_xxxxzz[k] * ab_z + g_y_0_yyzzz_xxxxzzz[k];

                g_y_0_yyzzzz_xxxyyy[k] = -g_y_0_yyzzz_xxxyyy[k] * ab_z + g_y_0_yyzzz_xxxyyyz[k];

                g_y_0_yyzzzz_xxxyyz[k] = -g_y_0_yyzzz_xxxyyz[k] * ab_z + g_y_0_yyzzz_xxxyyzz[k];

                g_y_0_yyzzzz_xxxyzz[k] = -g_y_0_yyzzz_xxxyzz[k] * ab_z + g_y_0_yyzzz_xxxyzzz[k];

                g_y_0_yyzzzz_xxxzzz[k] = -g_y_0_yyzzz_xxxzzz[k] * ab_z + g_y_0_yyzzz_xxxzzzz[k];

                g_y_0_yyzzzz_xxyyyy[k] = -g_y_0_yyzzz_xxyyyy[k] * ab_z + g_y_0_yyzzz_xxyyyyz[k];

                g_y_0_yyzzzz_xxyyyz[k] = -g_y_0_yyzzz_xxyyyz[k] * ab_z + g_y_0_yyzzz_xxyyyzz[k];

                g_y_0_yyzzzz_xxyyzz[k] = -g_y_0_yyzzz_xxyyzz[k] * ab_z + g_y_0_yyzzz_xxyyzzz[k];

                g_y_0_yyzzzz_xxyzzz[k] = -g_y_0_yyzzz_xxyzzz[k] * ab_z + g_y_0_yyzzz_xxyzzzz[k];

                g_y_0_yyzzzz_xxzzzz[k] = -g_y_0_yyzzz_xxzzzz[k] * ab_z + g_y_0_yyzzz_xxzzzzz[k];

                g_y_0_yyzzzz_xyyyyy[k] = -g_y_0_yyzzz_xyyyyy[k] * ab_z + g_y_0_yyzzz_xyyyyyz[k];

                g_y_0_yyzzzz_xyyyyz[k] = -g_y_0_yyzzz_xyyyyz[k] * ab_z + g_y_0_yyzzz_xyyyyzz[k];

                g_y_0_yyzzzz_xyyyzz[k] = -g_y_0_yyzzz_xyyyzz[k] * ab_z + g_y_0_yyzzz_xyyyzzz[k];

                g_y_0_yyzzzz_xyyzzz[k] = -g_y_0_yyzzz_xyyzzz[k] * ab_z + g_y_0_yyzzz_xyyzzzz[k];

                g_y_0_yyzzzz_xyzzzz[k] = -g_y_0_yyzzz_xyzzzz[k] * ab_z + g_y_0_yyzzz_xyzzzzz[k];

                g_y_0_yyzzzz_xzzzzz[k] = -g_y_0_yyzzz_xzzzzz[k] * ab_z + g_y_0_yyzzz_xzzzzzz[k];

                g_y_0_yyzzzz_yyyyyy[k] = -g_y_0_yyzzz_yyyyyy[k] * ab_z + g_y_0_yyzzz_yyyyyyz[k];

                g_y_0_yyzzzz_yyyyyz[k] = -g_y_0_yyzzz_yyyyyz[k] * ab_z + g_y_0_yyzzz_yyyyyzz[k];

                g_y_0_yyzzzz_yyyyzz[k] = -g_y_0_yyzzz_yyyyzz[k] * ab_z + g_y_0_yyzzz_yyyyzzz[k];

                g_y_0_yyzzzz_yyyzzz[k] = -g_y_0_yyzzz_yyyzzz[k] * ab_z + g_y_0_yyzzz_yyyzzzz[k];

                g_y_0_yyzzzz_yyzzzz[k] = -g_y_0_yyzzz_yyzzzz[k] * ab_z + g_y_0_yyzzz_yyzzzzz[k];

                g_y_0_yyzzzz_yzzzzz[k] = -g_y_0_yyzzz_yzzzzz[k] * ab_z + g_y_0_yyzzz_yzzzzzz[k];

                g_y_0_yyzzzz_zzzzzz[k] = -g_y_0_yyzzz_zzzzzz[k] * ab_z + g_y_0_yyzzz_zzzzzzz[k];
            }

            /// Set up 1512-1540 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1512 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1513 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1514 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1515 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1516 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1517 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1518 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1519 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1520 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1521 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1522 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1523 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1524 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1525 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1526 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1527 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1528 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1529 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1530 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1531 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1532 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1533 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1534 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1535 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1536 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1537 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1538 * ccomps * dcomps);

            auto g_y_0_yzzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1539 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzzzz_xxxxxx, g_y_0_yzzzz_xxxxxxz, g_y_0_yzzzz_xxxxxy, g_y_0_yzzzz_xxxxxyz, g_y_0_yzzzz_xxxxxz, g_y_0_yzzzz_xxxxxzz, g_y_0_yzzzz_xxxxyy, g_y_0_yzzzz_xxxxyyz, g_y_0_yzzzz_xxxxyz, g_y_0_yzzzz_xxxxyzz, g_y_0_yzzzz_xxxxzz, g_y_0_yzzzz_xxxxzzz, g_y_0_yzzzz_xxxyyy, g_y_0_yzzzz_xxxyyyz, g_y_0_yzzzz_xxxyyz, g_y_0_yzzzz_xxxyyzz, g_y_0_yzzzz_xxxyzz, g_y_0_yzzzz_xxxyzzz, g_y_0_yzzzz_xxxzzz, g_y_0_yzzzz_xxxzzzz, g_y_0_yzzzz_xxyyyy, g_y_0_yzzzz_xxyyyyz, g_y_0_yzzzz_xxyyyz, g_y_0_yzzzz_xxyyyzz, g_y_0_yzzzz_xxyyzz, g_y_0_yzzzz_xxyyzzz, g_y_0_yzzzz_xxyzzz, g_y_0_yzzzz_xxyzzzz, g_y_0_yzzzz_xxzzzz, g_y_0_yzzzz_xxzzzzz, g_y_0_yzzzz_xyyyyy, g_y_0_yzzzz_xyyyyyz, g_y_0_yzzzz_xyyyyz, g_y_0_yzzzz_xyyyyzz, g_y_0_yzzzz_xyyyzz, g_y_0_yzzzz_xyyyzzz, g_y_0_yzzzz_xyyzzz, g_y_0_yzzzz_xyyzzzz, g_y_0_yzzzz_xyzzzz, g_y_0_yzzzz_xyzzzzz, g_y_0_yzzzz_xzzzzz, g_y_0_yzzzz_xzzzzzz, g_y_0_yzzzz_yyyyyy, g_y_0_yzzzz_yyyyyyz, g_y_0_yzzzz_yyyyyz, g_y_0_yzzzz_yyyyyzz, g_y_0_yzzzz_yyyyzz, g_y_0_yzzzz_yyyyzzz, g_y_0_yzzzz_yyyzzz, g_y_0_yzzzz_yyyzzzz, g_y_0_yzzzz_yyzzzz, g_y_0_yzzzz_yyzzzzz, g_y_0_yzzzz_yzzzzz, g_y_0_yzzzz_yzzzzzz, g_y_0_yzzzz_zzzzzz, g_y_0_yzzzz_zzzzzzz, g_y_0_yzzzzz_xxxxxx, g_y_0_yzzzzz_xxxxxy, g_y_0_yzzzzz_xxxxxz, g_y_0_yzzzzz_xxxxyy, g_y_0_yzzzzz_xxxxyz, g_y_0_yzzzzz_xxxxzz, g_y_0_yzzzzz_xxxyyy, g_y_0_yzzzzz_xxxyyz, g_y_0_yzzzzz_xxxyzz, g_y_0_yzzzzz_xxxzzz, g_y_0_yzzzzz_xxyyyy, g_y_0_yzzzzz_xxyyyz, g_y_0_yzzzzz_xxyyzz, g_y_0_yzzzzz_xxyzzz, g_y_0_yzzzzz_xxzzzz, g_y_0_yzzzzz_xyyyyy, g_y_0_yzzzzz_xyyyyz, g_y_0_yzzzzz_xyyyzz, g_y_0_yzzzzz_xyyzzz, g_y_0_yzzzzz_xyzzzz, g_y_0_yzzzzz_xzzzzz, g_y_0_yzzzzz_yyyyyy, g_y_0_yzzzzz_yyyyyz, g_y_0_yzzzzz_yyyyzz, g_y_0_yzzzzz_yyyzzz, g_y_0_yzzzzz_yyzzzz, g_y_0_yzzzzz_yzzzzz, g_y_0_yzzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_xxxxxx[k] = -g_y_0_yzzzz_xxxxxx[k] * ab_z + g_y_0_yzzzz_xxxxxxz[k];

                g_y_0_yzzzzz_xxxxxy[k] = -g_y_0_yzzzz_xxxxxy[k] * ab_z + g_y_0_yzzzz_xxxxxyz[k];

                g_y_0_yzzzzz_xxxxxz[k] = -g_y_0_yzzzz_xxxxxz[k] * ab_z + g_y_0_yzzzz_xxxxxzz[k];

                g_y_0_yzzzzz_xxxxyy[k] = -g_y_0_yzzzz_xxxxyy[k] * ab_z + g_y_0_yzzzz_xxxxyyz[k];

                g_y_0_yzzzzz_xxxxyz[k] = -g_y_0_yzzzz_xxxxyz[k] * ab_z + g_y_0_yzzzz_xxxxyzz[k];

                g_y_0_yzzzzz_xxxxzz[k] = -g_y_0_yzzzz_xxxxzz[k] * ab_z + g_y_0_yzzzz_xxxxzzz[k];

                g_y_0_yzzzzz_xxxyyy[k] = -g_y_0_yzzzz_xxxyyy[k] * ab_z + g_y_0_yzzzz_xxxyyyz[k];

                g_y_0_yzzzzz_xxxyyz[k] = -g_y_0_yzzzz_xxxyyz[k] * ab_z + g_y_0_yzzzz_xxxyyzz[k];

                g_y_0_yzzzzz_xxxyzz[k] = -g_y_0_yzzzz_xxxyzz[k] * ab_z + g_y_0_yzzzz_xxxyzzz[k];

                g_y_0_yzzzzz_xxxzzz[k] = -g_y_0_yzzzz_xxxzzz[k] * ab_z + g_y_0_yzzzz_xxxzzzz[k];

                g_y_0_yzzzzz_xxyyyy[k] = -g_y_0_yzzzz_xxyyyy[k] * ab_z + g_y_0_yzzzz_xxyyyyz[k];

                g_y_0_yzzzzz_xxyyyz[k] = -g_y_0_yzzzz_xxyyyz[k] * ab_z + g_y_0_yzzzz_xxyyyzz[k];

                g_y_0_yzzzzz_xxyyzz[k] = -g_y_0_yzzzz_xxyyzz[k] * ab_z + g_y_0_yzzzz_xxyyzzz[k];

                g_y_0_yzzzzz_xxyzzz[k] = -g_y_0_yzzzz_xxyzzz[k] * ab_z + g_y_0_yzzzz_xxyzzzz[k];

                g_y_0_yzzzzz_xxzzzz[k] = -g_y_0_yzzzz_xxzzzz[k] * ab_z + g_y_0_yzzzz_xxzzzzz[k];

                g_y_0_yzzzzz_xyyyyy[k] = -g_y_0_yzzzz_xyyyyy[k] * ab_z + g_y_0_yzzzz_xyyyyyz[k];

                g_y_0_yzzzzz_xyyyyz[k] = -g_y_0_yzzzz_xyyyyz[k] * ab_z + g_y_0_yzzzz_xyyyyzz[k];

                g_y_0_yzzzzz_xyyyzz[k] = -g_y_0_yzzzz_xyyyzz[k] * ab_z + g_y_0_yzzzz_xyyyzzz[k];

                g_y_0_yzzzzz_xyyzzz[k] = -g_y_0_yzzzz_xyyzzz[k] * ab_z + g_y_0_yzzzz_xyyzzzz[k];

                g_y_0_yzzzzz_xyzzzz[k] = -g_y_0_yzzzz_xyzzzz[k] * ab_z + g_y_0_yzzzz_xyzzzzz[k];

                g_y_0_yzzzzz_xzzzzz[k] = -g_y_0_yzzzz_xzzzzz[k] * ab_z + g_y_0_yzzzz_xzzzzzz[k];

                g_y_0_yzzzzz_yyyyyy[k] = -g_y_0_yzzzz_yyyyyy[k] * ab_z + g_y_0_yzzzz_yyyyyyz[k];

                g_y_0_yzzzzz_yyyyyz[k] = -g_y_0_yzzzz_yyyyyz[k] * ab_z + g_y_0_yzzzz_yyyyyzz[k];

                g_y_0_yzzzzz_yyyyzz[k] = -g_y_0_yzzzz_yyyyzz[k] * ab_z + g_y_0_yzzzz_yyyyzzz[k];

                g_y_0_yzzzzz_yyyzzz[k] = -g_y_0_yzzzz_yyyzzz[k] * ab_z + g_y_0_yzzzz_yyyzzzz[k];

                g_y_0_yzzzzz_yyzzzz[k] = -g_y_0_yzzzz_yyzzzz[k] * ab_z + g_y_0_yzzzz_yyzzzzz[k];

                g_y_0_yzzzzz_yzzzzz[k] = -g_y_0_yzzzz_yzzzzz[k] * ab_z + g_y_0_yzzzz_yzzzzzz[k];

                g_y_0_yzzzzz_zzzzzz[k] = -g_y_0_yzzzz_zzzzzz[k] * ab_z + g_y_0_yzzzz_zzzzzzz[k];
            }

            /// Set up 1540-1568 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1540 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1541 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1542 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1543 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1544 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1545 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1546 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1547 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1548 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1549 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1550 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1551 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1552 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1553 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1554 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1555 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1556 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1557 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1558 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1559 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1560 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1561 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1562 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1563 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1564 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1565 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1566 * ccomps * dcomps);

            auto g_y_0_zzzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1567 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzzz_xxxxxx, g_y_0_zzzzz_xxxxxxz, g_y_0_zzzzz_xxxxxy, g_y_0_zzzzz_xxxxxyz, g_y_0_zzzzz_xxxxxz, g_y_0_zzzzz_xxxxxzz, g_y_0_zzzzz_xxxxyy, g_y_0_zzzzz_xxxxyyz, g_y_0_zzzzz_xxxxyz, g_y_0_zzzzz_xxxxyzz, g_y_0_zzzzz_xxxxzz, g_y_0_zzzzz_xxxxzzz, g_y_0_zzzzz_xxxyyy, g_y_0_zzzzz_xxxyyyz, g_y_0_zzzzz_xxxyyz, g_y_0_zzzzz_xxxyyzz, g_y_0_zzzzz_xxxyzz, g_y_0_zzzzz_xxxyzzz, g_y_0_zzzzz_xxxzzz, g_y_0_zzzzz_xxxzzzz, g_y_0_zzzzz_xxyyyy, g_y_0_zzzzz_xxyyyyz, g_y_0_zzzzz_xxyyyz, g_y_0_zzzzz_xxyyyzz, g_y_0_zzzzz_xxyyzz, g_y_0_zzzzz_xxyyzzz, g_y_0_zzzzz_xxyzzz, g_y_0_zzzzz_xxyzzzz, g_y_0_zzzzz_xxzzzz, g_y_0_zzzzz_xxzzzzz, g_y_0_zzzzz_xyyyyy, g_y_0_zzzzz_xyyyyyz, g_y_0_zzzzz_xyyyyz, g_y_0_zzzzz_xyyyyzz, g_y_0_zzzzz_xyyyzz, g_y_0_zzzzz_xyyyzzz, g_y_0_zzzzz_xyyzzz, g_y_0_zzzzz_xyyzzzz, g_y_0_zzzzz_xyzzzz, g_y_0_zzzzz_xyzzzzz, g_y_0_zzzzz_xzzzzz, g_y_0_zzzzz_xzzzzzz, g_y_0_zzzzz_yyyyyy, g_y_0_zzzzz_yyyyyyz, g_y_0_zzzzz_yyyyyz, g_y_0_zzzzz_yyyyyzz, g_y_0_zzzzz_yyyyzz, g_y_0_zzzzz_yyyyzzz, g_y_0_zzzzz_yyyzzz, g_y_0_zzzzz_yyyzzzz, g_y_0_zzzzz_yyzzzz, g_y_0_zzzzz_yyzzzzz, g_y_0_zzzzz_yzzzzz, g_y_0_zzzzz_yzzzzzz, g_y_0_zzzzz_zzzzzz, g_y_0_zzzzz_zzzzzzz, g_y_0_zzzzzz_xxxxxx, g_y_0_zzzzzz_xxxxxy, g_y_0_zzzzzz_xxxxxz, g_y_0_zzzzzz_xxxxyy, g_y_0_zzzzzz_xxxxyz, g_y_0_zzzzzz_xxxxzz, g_y_0_zzzzzz_xxxyyy, g_y_0_zzzzzz_xxxyyz, g_y_0_zzzzzz_xxxyzz, g_y_0_zzzzzz_xxxzzz, g_y_0_zzzzzz_xxyyyy, g_y_0_zzzzzz_xxyyyz, g_y_0_zzzzzz_xxyyzz, g_y_0_zzzzzz_xxyzzz, g_y_0_zzzzzz_xxzzzz, g_y_0_zzzzzz_xyyyyy, g_y_0_zzzzzz_xyyyyz, g_y_0_zzzzzz_xyyyzz, g_y_0_zzzzzz_xyyzzz, g_y_0_zzzzzz_xyzzzz, g_y_0_zzzzzz_xzzzzz, g_y_0_zzzzzz_yyyyyy, g_y_0_zzzzzz_yyyyyz, g_y_0_zzzzzz_yyyyzz, g_y_0_zzzzzz_yyyzzz, g_y_0_zzzzzz_yyzzzz, g_y_0_zzzzzz_yzzzzz, g_y_0_zzzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_xxxxxx[k] = -g_y_0_zzzzz_xxxxxx[k] * ab_z + g_y_0_zzzzz_xxxxxxz[k];

                g_y_0_zzzzzz_xxxxxy[k] = -g_y_0_zzzzz_xxxxxy[k] * ab_z + g_y_0_zzzzz_xxxxxyz[k];

                g_y_0_zzzzzz_xxxxxz[k] = -g_y_0_zzzzz_xxxxxz[k] * ab_z + g_y_0_zzzzz_xxxxxzz[k];

                g_y_0_zzzzzz_xxxxyy[k] = -g_y_0_zzzzz_xxxxyy[k] * ab_z + g_y_0_zzzzz_xxxxyyz[k];

                g_y_0_zzzzzz_xxxxyz[k] = -g_y_0_zzzzz_xxxxyz[k] * ab_z + g_y_0_zzzzz_xxxxyzz[k];

                g_y_0_zzzzzz_xxxxzz[k] = -g_y_0_zzzzz_xxxxzz[k] * ab_z + g_y_0_zzzzz_xxxxzzz[k];

                g_y_0_zzzzzz_xxxyyy[k] = -g_y_0_zzzzz_xxxyyy[k] * ab_z + g_y_0_zzzzz_xxxyyyz[k];

                g_y_0_zzzzzz_xxxyyz[k] = -g_y_0_zzzzz_xxxyyz[k] * ab_z + g_y_0_zzzzz_xxxyyzz[k];

                g_y_0_zzzzzz_xxxyzz[k] = -g_y_0_zzzzz_xxxyzz[k] * ab_z + g_y_0_zzzzz_xxxyzzz[k];

                g_y_0_zzzzzz_xxxzzz[k] = -g_y_0_zzzzz_xxxzzz[k] * ab_z + g_y_0_zzzzz_xxxzzzz[k];

                g_y_0_zzzzzz_xxyyyy[k] = -g_y_0_zzzzz_xxyyyy[k] * ab_z + g_y_0_zzzzz_xxyyyyz[k];

                g_y_0_zzzzzz_xxyyyz[k] = -g_y_0_zzzzz_xxyyyz[k] * ab_z + g_y_0_zzzzz_xxyyyzz[k];

                g_y_0_zzzzzz_xxyyzz[k] = -g_y_0_zzzzz_xxyyzz[k] * ab_z + g_y_0_zzzzz_xxyyzzz[k];

                g_y_0_zzzzzz_xxyzzz[k] = -g_y_0_zzzzz_xxyzzz[k] * ab_z + g_y_0_zzzzz_xxyzzzz[k];

                g_y_0_zzzzzz_xxzzzz[k] = -g_y_0_zzzzz_xxzzzz[k] * ab_z + g_y_0_zzzzz_xxzzzzz[k];

                g_y_0_zzzzzz_xyyyyy[k] = -g_y_0_zzzzz_xyyyyy[k] * ab_z + g_y_0_zzzzz_xyyyyyz[k];

                g_y_0_zzzzzz_xyyyyz[k] = -g_y_0_zzzzz_xyyyyz[k] * ab_z + g_y_0_zzzzz_xyyyyzz[k];

                g_y_0_zzzzzz_xyyyzz[k] = -g_y_0_zzzzz_xyyyzz[k] * ab_z + g_y_0_zzzzz_xyyyzzz[k];

                g_y_0_zzzzzz_xyyzzz[k] = -g_y_0_zzzzz_xyyzzz[k] * ab_z + g_y_0_zzzzz_xyyzzzz[k];

                g_y_0_zzzzzz_xyzzzz[k] = -g_y_0_zzzzz_xyzzzz[k] * ab_z + g_y_0_zzzzz_xyzzzzz[k];

                g_y_0_zzzzzz_xzzzzz[k] = -g_y_0_zzzzz_xzzzzz[k] * ab_z + g_y_0_zzzzz_xzzzzzz[k];

                g_y_0_zzzzzz_yyyyyy[k] = -g_y_0_zzzzz_yyyyyy[k] * ab_z + g_y_0_zzzzz_yyyyyyz[k];

                g_y_0_zzzzzz_yyyyyz[k] = -g_y_0_zzzzz_yyyyyz[k] * ab_z + g_y_0_zzzzz_yyyyyzz[k];

                g_y_0_zzzzzz_yyyyzz[k] = -g_y_0_zzzzz_yyyyzz[k] * ab_z + g_y_0_zzzzz_yyyyzzz[k];

                g_y_0_zzzzzz_yyyzzz[k] = -g_y_0_zzzzz_yyyzzz[k] * ab_z + g_y_0_zzzzz_yyyzzzz[k];

                g_y_0_zzzzzz_yyzzzz[k] = -g_y_0_zzzzz_yyzzzz[k] * ab_z + g_y_0_zzzzz_yyzzzzz[k];

                g_y_0_zzzzzz_yzzzzz[k] = -g_y_0_zzzzz_yzzzzz[k] * ab_z + g_y_0_zzzzz_yzzzzzz[k];

                g_y_0_zzzzzz_zzzzzz[k] = -g_y_0_zzzzz_zzzzzz[k] * ab_z + g_y_0_zzzzz_zzzzzzz[k];
            }

            /// Set up 1568-1596 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_xxxxxx = cbuffer.data(ii_geom_10_off + 1568 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxxxy = cbuffer.data(ii_geom_10_off + 1569 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxxxz = cbuffer.data(ii_geom_10_off + 1570 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxxyy = cbuffer.data(ii_geom_10_off + 1571 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxxyz = cbuffer.data(ii_geom_10_off + 1572 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxxzz = cbuffer.data(ii_geom_10_off + 1573 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxyyy = cbuffer.data(ii_geom_10_off + 1574 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxyyz = cbuffer.data(ii_geom_10_off + 1575 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxyzz = cbuffer.data(ii_geom_10_off + 1576 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxzzz = cbuffer.data(ii_geom_10_off + 1577 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxyyyy = cbuffer.data(ii_geom_10_off + 1578 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxyyyz = cbuffer.data(ii_geom_10_off + 1579 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxyyzz = cbuffer.data(ii_geom_10_off + 1580 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxyzzz = cbuffer.data(ii_geom_10_off + 1581 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxzzzz = cbuffer.data(ii_geom_10_off + 1582 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyyyyy = cbuffer.data(ii_geom_10_off + 1583 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyyyyz = cbuffer.data(ii_geom_10_off + 1584 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyyyzz = cbuffer.data(ii_geom_10_off + 1585 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyyzzz = cbuffer.data(ii_geom_10_off + 1586 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyzzzz = cbuffer.data(ii_geom_10_off + 1587 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xzzzzz = cbuffer.data(ii_geom_10_off + 1588 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyyyyy = cbuffer.data(ii_geom_10_off + 1589 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyyyyz = cbuffer.data(ii_geom_10_off + 1590 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyyyzz = cbuffer.data(ii_geom_10_off + 1591 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyyzzz = cbuffer.data(ii_geom_10_off + 1592 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyzzzz = cbuffer.data(ii_geom_10_off + 1593 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yzzzzz = cbuffer.data(ii_geom_10_off + 1594 * ccomps * dcomps);

            auto g_z_0_xxxxxx_zzzzzz = cbuffer.data(ii_geom_10_off + 1595 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxx_xxxxxx, g_z_0_xxxxx_xxxxxxx, g_z_0_xxxxx_xxxxxxy, g_z_0_xxxxx_xxxxxxz, g_z_0_xxxxx_xxxxxy, g_z_0_xxxxx_xxxxxyy, g_z_0_xxxxx_xxxxxyz, g_z_0_xxxxx_xxxxxz, g_z_0_xxxxx_xxxxxzz, g_z_0_xxxxx_xxxxyy, g_z_0_xxxxx_xxxxyyy, g_z_0_xxxxx_xxxxyyz, g_z_0_xxxxx_xxxxyz, g_z_0_xxxxx_xxxxyzz, g_z_0_xxxxx_xxxxzz, g_z_0_xxxxx_xxxxzzz, g_z_0_xxxxx_xxxyyy, g_z_0_xxxxx_xxxyyyy, g_z_0_xxxxx_xxxyyyz, g_z_0_xxxxx_xxxyyz, g_z_0_xxxxx_xxxyyzz, g_z_0_xxxxx_xxxyzz, g_z_0_xxxxx_xxxyzzz, g_z_0_xxxxx_xxxzzz, g_z_0_xxxxx_xxxzzzz, g_z_0_xxxxx_xxyyyy, g_z_0_xxxxx_xxyyyyy, g_z_0_xxxxx_xxyyyyz, g_z_0_xxxxx_xxyyyz, g_z_0_xxxxx_xxyyyzz, g_z_0_xxxxx_xxyyzz, g_z_0_xxxxx_xxyyzzz, g_z_0_xxxxx_xxyzzz, g_z_0_xxxxx_xxyzzzz, g_z_0_xxxxx_xxzzzz, g_z_0_xxxxx_xxzzzzz, g_z_0_xxxxx_xyyyyy, g_z_0_xxxxx_xyyyyyy, g_z_0_xxxxx_xyyyyyz, g_z_0_xxxxx_xyyyyz, g_z_0_xxxxx_xyyyyzz, g_z_0_xxxxx_xyyyzz, g_z_0_xxxxx_xyyyzzz, g_z_0_xxxxx_xyyzzz, g_z_0_xxxxx_xyyzzzz, g_z_0_xxxxx_xyzzzz, g_z_0_xxxxx_xyzzzzz, g_z_0_xxxxx_xzzzzz, g_z_0_xxxxx_xzzzzzz, g_z_0_xxxxx_yyyyyy, g_z_0_xxxxx_yyyyyz, g_z_0_xxxxx_yyyyzz, g_z_0_xxxxx_yyyzzz, g_z_0_xxxxx_yyzzzz, g_z_0_xxxxx_yzzzzz, g_z_0_xxxxx_zzzzzz, g_z_0_xxxxxx_xxxxxx, g_z_0_xxxxxx_xxxxxy, g_z_0_xxxxxx_xxxxxz, g_z_0_xxxxxx_xxxxyy, g_z_0_xxxxxx_xxxxyz, g_z_0_xxxxxx_xxxxzz, g_z_0_xxxxxx_xxxyyy, g_z_0_xxxxxx_xxxyyz, g_z_0_xxxxxx_xxxyzz, g_z_0_xxxxxx_xxxzzz, g_z_0_xxxxxx_xxyyyy, g_z_0_xxxxxx_xxyyyz, g_z_0_xxxxxx_xxyyzz, g_z_0_xxxxxx_xxyzzz, g_z_0_xxxxxx_xxzzzz, g_z_0_xxxxxx_xyyyyy, g_z_0_xxxxxx_xyyyyz, g_z_0_xxxxxx_xyyyzz, g_z_0_xxxxxx_xyyzzz, g_z_0_xxxxxx_xyzzzz, g_z_0_xxxxxx_xzzzzz, g_z_0_xxxxxx_yyyyyy, g_z_0_xxxxxx_yyyyyz, g_z_0_xxxxxx_yyyyzz, g_z_0_xxxxxx_yyyzzz, g_z_0_xxxxxx_yyzzzz, g_z_0_xxxxxx_yzzzzz, g_z_0_xxxxxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_xxxxxx[k] = -g_z_0_xxxxx_xxxxxx[k] * ab_x + g_z_0_xxxxx_xxxxxxx[k];

                g_z_0_xxxxxx_xxxxxy[k] = -g_z_0_xxxxx_xxxxxy[k] * ab_x + g_z_0_xxxxx_xxxxxxy[k];

                g_z_0_xxxxxx_xxxxxz[k] = -g_z_0_xxxxx_xxxxxz[k] * ab_x + g_z_0_xxxxx_xxxxxxz[k];

                g_z_0_xxxxxx_xxxxyy[k] = -g_z_0_xxxxx_xxxxyy[k] * ab_x + g_z_0_xxxxx_xxxxxyy[k];

                g_z_0_xxxxxx_xxxxyz[k] = -g_z_0_xxxxx_xxxxyz[k] * ab_x + g_z_0_xxxxx_xxxxxyz[k];

                g_z_0_xxxxxx_xxxxzz[k] = -g_z_0_xxxxx_xxxxzz[k] * ab_x + g_z_0_xxxxx_xxxxxzz[k];

                g_z_0_xxxxxx_xxxyyy[k] = -g_z_0_xxxxx_xxxyyy[k] * ab_x + g_z_0_xxxxx_xxxxyyy[k];

                g_z_0_xxxxxx_xxxyyz[k] = -g_z_0_xxxxx_xxxyyz[k] * ab_x + g_z_0_xxxxx_xxxxyyz[k];

                g_z_0_xxxxxx_xxxyzz[k] = -g_z_0_xxxxx_xxxyzz[k] * ab_x + g_z_0_xxxxx_xxxxyzz[k];

                g_z_0_xxxxxx_xxxzzz[k] = -g_z_0_xxxxx_xxxzzz[k] * ab_x + g_z_0_xxxxx_xxxxzzz[k];

                g_z_0_xxxxxx_xxyyyy[k] = -g_z_0_xxxxx_xxyyyy[k] * ab_x + g_z_0_xxxxx_xxxyyyy[k];

                g_z_0_xxxxxx_xxyyyz[k] = -g_z_0_xxxxx_xxyyyz[k] * ab_x + g_z_0_xxxxx_xxxyyyz[k];

                g_z_0_xxxxxx_xxyyzz[k] = -g_z_0_xxxxx_xxyyzz[k] * ab_x + g_z_0_xxxxx_xxxyyzz[k];

                g_z_0_xxxxxx_xxyzzz[k] = -g_z_0_xxxxx_xxyzzz[k] * ab_x + g_z_0_xxxxx_xxxyzzz[k];

                g_z_0_xxxxxx_xxzzzz[k] = -g_z_0_xxxxx_xxzzzz[k] * ab_x + g_z_0_xxxxx_xxxzzzz[k];

                g_z_0_xxxxxx_xyyyyy[k] = -g_z_0_xxxxx_xyyyyy[k] * ab_x + g_z_0_xxxxx_xxyyyyy[k];

                g_z_0_xxxxxx_xyyyyz[k] = -g_z_0_xxxxx_xyyyyz[k] * ab_x + g_z_0_xxxxx_xxyyyyz[k];

                g_z_0_xxxxxx_xyyyzz[k] = -g_z_0_xxxxx_xyyyzz[k] * ab_x + g_z_0_xxxxx_xxyyyzz[k];

                g_z_0_xxxxxx_xyyzzz[k] = -g_z_0_xxxxx_xyyzzz[k] * ab_x + g_z_0_xxxxx_xxyyzzz[k];

                g_z_0_xxxxxx_xyzzzz[k] = -g_z_0_xxxxx_xyzzzz[k] * ab_x + g_z_0_xxxxx_xxyzzzz[k];

                g_z_0_xxxxxx_xzzzzz[k] = -g_z_0_xxxxx_xzzzzz[k] * ab_x + g_z_0_xxxxx_xxzzzzz[k];

                g_z_0_xxxxxx_yyyyyy[k] = -g_z_0_xxxxx_yyyyyy[k] * ab_x + g_z_0_xxxxx_xyyyyyy[k];

                g_z_0_xxxxxx_yyyyyz[k] = -g_z_0_xxxxx_yyyyyz[k] * ab_x + g_z_0_xxxxx_xyyyyyz[k];

                g_z_0_xxxxxx_yyyyzz[k] = -g_z_0_xxxxx_yyyyzz[k] * ab_x + g_z_0_xxxxx_xyyyyzz[k];

                g_z_0_xxxxxx_yyyzzz[k] = -g_z_0_xxxxx_yyyzzz[k] * ab_x + g_z_0_xxxxx_xyyyzzz[k];

                g_z_0_xxxxxx_yyzzzz[k] = -g_z_0_xxxxx_yyzzzz[k] * ab_x + g_z_0_xxxxx_xyyzzzz[k];

                g_z_0_xxxxxx_yzzzzz[k] = -g_z_0_xxxxx_yzzzzz[k] * ab_x + g_z_0_xxxxx_xyzzzzz[k];

                g_z_0_xxxxxx_zzzzzz[k] = -g_z_0_xxxxx_zzzzzz[k] * ab_x + g_z_0_xxxxx_xzzzzzz[k];
            }

            /// Set up 1596-1624 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_xxxxxx = cbuffer.data(ii_geom_10_off + 1596 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxxxy = cbuffer.data(ii_geom_10_off + 1597 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxxxz = cbuffer.data(ii_geom_10_off + 1598 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxxyy = cbuffer.data(ii_geom_10_off + 1599 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxxyz = cbuffer.data(ii_geom_10_off + 1600 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxxzz = cbuffer.data(ii_geom_10_off + 1601 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxyyy = cbuffer.data(ii_geom_10_off + 1602 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxyyz = cbuffer.data(ii_geom_10_off + 1603 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxyzz = cbuffer.data(ii_geom_10_off + 1604 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxzzz = cbuffer.data(ii_geom_10_off + 1605 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxyyyy = cbuffer.data(ii_geom_10_off + 1606 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxyyyz = cbuffer.data(ii_geom_10_off + 1607 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxyyzz = cbuffer.data(ii_geom_10_off + 1608 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxyzzz = cbuffer.data(ii_geom_10_off + 1609 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxzzzz = cbuffer.data(ii_geom_10_off + 1610 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyyyyy = cbuffer.data(ii_geom_10_off + 1611 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyyyyz = cbuffer.data(ii_geom_10_off + 1612 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyyyzz = cbuffer.data(ii_geom_10_off + 1613 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyyzzz = cbuffer.data(ii_geom_10_off + 1614 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyzzzz = cbuffer.data(ii_geom_10_off + 1615 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xzzzzz = cbuffer.data(ii_geom_10_off + 1616 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyyyyy = cbuffer.data(ii_geom_10_off + 1617 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyyyyz = cbuffer.data(ii_geom_10_off + 1618 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyyyzz = cbuffer.data(ii_geom_10_off + 1619 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyyzzz = cbuffer.data(ii_geom_10_off + 1620 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyzzzz = cbuffer.data(ii_geom_10_off + 1621 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yzzzzz = cbuffer.data(ii_geom_10_off + 1622 * ccomps * dcomps);

            auto g_z_0_xxxxxy_zzzzzz = cbuffer.data(ii_geom_10_off + 1623 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxxy_xxxxxx, g_z_0_xxxxxy_xxxxxy, g_z_0_xxxxxy_xxxxxz, g_z_0_xxxxxy_xxxxyy, g_z_0_xxxxxy_xxxxyz, g_z_0_xxxxxy_xxxxzz, g_z_0_xxxxxy_xxxyyy, g_z_0_xxxxxy_xxxyyz, g_z_0_xxxxxy_xxxyzz, g_z_0_xxxxxy_xxxzzz, g_z_0_xxxxxy_xxyyyy, g_z_0_xxxxxy_xxyyyz, g_z_0_xxxxxy_xxyyzz, g_z_0_xxxxxy_xxyzzz, g_z_0_xxxxxy_xxzzzz, g_z_0_xxxxxy_xyyyyy, g_z_0_xxxxxy_xyyyyz, g_z_0_xxxxxy_xyyyzz, g_z_0_xxxxxy_xyyzzz, g_z_0_xxxxxy_xyzzzz, g_z_0_xxxxxy_xzzzzz, g_z_0_xxxxxy_yyyyyy, g_z_0_xxxxxy_yyyyyz, g_z_0_xxxxxy_yyyyzz, g_z_0_xxxxxy_yyyzzz, g_z_0_xxxxxy_yyzzzz, g_z_0_xxxxxy_yzzzzz, g_z_0_xxxxxy_zzzzzz, g_z_0_xxxxy_xxxxxx, g_z_0_xxxxy_xxxxxxx, g_z_0_xxxxy_xxxxxxy, g_z_0_xxxxy_xxxxxxz, g_z_0_xxxxy_xxxxxy, g_z_0_xxxxy_xxxxxyy, g_z_0_xxxxy_xxxxxyz, g_z_0_xxxxy_xxxxxz, g_z_0_xxxxy_xxxxxzz, g_z_0_xxxxy_xxxxyy, g_z_0_xxxxy_xxxxyyy, g_z_0_xxxxy_xxxxyyz, g_z_0_xxxxy_xxxxyz, g_z_0_xxxxy_xxxxyzz, g_z_0_xxxxy_xxxxzz, g_z_0_xxxxy_xxxxzzz, g_z_0_xxxxy_xxxyyy, g_z_0_xxxxy_xxxyyyy, g_z_0_xxxxy_xxxyyyz, g_z_0_xxxxy_xxxyyz, g_z_0_xxxxy_xxxyyzz, g_z_0_xxxxy_xxxyzz, g_z_0_xxxxy_xxxyzzz, g_z_0_xxxxy_xxxzzz, g_z_0_xxxxy_xxxzzzz, g_z_0_xxxxy_xxyyyy, g_z_0_xxxxy_xxyyyyy, g_z_0_xxxxy_xxyyyyz, g_z_0_xxxxy_xxyyyz, g_z_0_xxxxy_xxyyyzz, g_z_0_xxxxy_xxyyzz, g_z_0_xxxxy_xxyyzzz, g_z_0_xxxxy_xxyzzz, g_z_0_xxxxy_xxyzzzz, g_z_0_xxxxy_xxzzzz, g_z_0_xxxxy_xxzzzzz, g_z_0_xxxxy_xyyyyy, g_z_0_xxxxy_xyyyyyy, g_z_0_xxxxy_xyyyyyz, g_z_0_xxxxy_xyyyyz, g_z_0_xxxxy_xyyyyzz, g_z_0_xxxxy_xyyyzz, g_z_0_xxxxy_xyyyzzz, g_z_0_xxxxy_xyyzzz, g_z_0_xxxxy_xyyzzzz, g_z_0_xxxxy_xyzzzz, g_z_0_xxxxy_xyzzzzz, g_z_0_xxxxy_xzzzzz, g_z_0_xxxxy_xzzzzzz, g_z_0_xxxxy_yyyyyy, g_z_0_xxxxy_yyyyyz, g_z_0_xxxxy_yyyyzz, g_z_0_xxxxy_yyyzzz, g_z_0_xxxxy_yyzzzz, g_z_0_xxxxy_yzzzzz, g_z_0_xxxxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_xxxxxx[k] = -g_z_0_xxxxy_xxxxxx[k] * ab_x + g_z_0_xxxxy_xxxxxxx[k];

                g_z_0_xxxxxy_xxxxxy[k] = -g_z_0_xxxxy_xxxxxy[k] * ab_x + g_z_0_xxxxy_xxxxxxy[k];

                g_z_0_xxxxxy_xxxxxz[k] = -g_z_0_xxxxy_xxxxxz[k] * ab_x + g_z_0_xxxxy_xxxxxxz[k];

                g_z_0_xxxxxy_xxxxyy[k] = -g_z_0_xxxxy_xxxxyy[k] * ab_x + g_z_0_xxxxy_xxxxxyy[k];

                g_z_0_xxxxxy_xxxxyz[k] = -g_z_0_xxxxy_xxxxyz[k] * ab_x + g_z_0_xxxxy_xxxxxyz[k];

                g_z_0_xxxxxy_xxxxzz[k] = -g_z_0_xxxxy_xxxxzz[k] * ab_x + g_z_0_xxxxy_xxxxxzz[k];

                g_z_0_xxxxxy_xxxyyy[k] = -g_z_0_xxxxy_xxxyyy[k] * ab_x + g_z_0_xxxxy_xxxxyyy[k];

                g_z_0_xxxxxy_xxxyyz[k] = -g_z_0_xxxxy_xxxyyz[k] * ab_x + g_z_0_xxxxy_xxxxyyz[k];

                g_z_0_xxxxxy_xxxyzz[k] = -g_z_0_xxxxy_xxxyzz[k] * ab_x + g_z_0_xxxxy_xxxxyzz[k];

                g_z_0_xxxxxy_xxxzzz[k] = -g_z_0_xxxxy_xxxzzz[k] * ab_x + g_z_0_xxxxy_xxxxzzz[k];

                g_z_0_xxxxxy_xxyyyy[k] = -g_z_0_xxxxy_xxyyyy[k] * ab_x + g_z_0_xxxxy_xxxyyyy[k];

                g_z_0_xxxxxy_xxyyyz[k] = -g_z_0_xxxxy_xxyyyz[k] * ab_x + g_z_0_xxxxy_xxxyyyz[k];

                g_z_0_xxxxxy_xxyyzz[k] = -g_z_0_xxxxy_xxyyzz[k] * ab_x + g_z_0_xxxxy_xxxyyzz[k];

                g_z_0_xxxxxy_xxyzzz[k] = -g_z_0_xxxxy_xxyzzz[k] * ab_x + g_z_0_xxxxy_xxxyzzz[k];

                g_z_0_xxxxxy_xxzzzz[k] = -g_z_0_xxxxy_xxzzzz[k] * ab_x + g_z_0_xxxxy_xxxzzzz[k];

                g_z_0_xxxxxy_xyyyyy[k] = -g_z_0_xxxxy_xyyyyy[k] * ab_x + g_z_0_xxxxy_xxyyyyy[k];

                g_z_0_xxxxxy_xyyyyz[k] = -g_z_0_xxxxy_xyyyyz[k] * ab_x + g_z_0_xxxxy_xxyyyyz[k];

                g_z_0_xxxxxy_xyyyzz[k] = -g_z_0_xxxxy_xyyyzz[k] * ab_x + g_z_0_xxxxy_xxyyyzz[k];

                g_z_0_xxxxxy_xyyzzz[k] = -g_z_0_xxxxy_xyyzzz[k] * ab_x + g_z_0_xxxxy_xxyyzzz[k];

                g_z_0_xxxxxy_xyzzzz[k] = -g_z_0_xxxxy_xyzzzz[k] * ab_x + g_z_0_xxxxy_xxyzzzz[k];

                g_z_0_xxxxxy_xzzzzz[k] = -g_z_0_xxxxy_xzzzzz[k] * ab_x + g_z_0_xxxxy_xxzzzzz[k];

                g_z_0_xxxxxy_yyyyyy[k] = -g_z_0_xxxxy_yyyyyy[k] * ab_x + g_z_0_xxxxy_xyyyyyy[k];

                g_z_0_xxxxxy_yyyyyz[k] = -g_z_0_xxxxy_yyyyyz[k] * ab_x + g_z_0_xxxxy_xyyyyyz[k];

                g_z_0_xxxxxy_yyyyzz[k] = -g_z_0_xxxxy_yyyyzz[k] * ab_x + g_z_0_xxxxy_xyyyyzz[k];

                g_z_0_xxxxxy_yyyzzz[k] = -g_z_0_xxxxy_yyyzzz[k] * ab_x + g_z_0_xxxxy_xyyyzzz[k];

                g_z_0_xxxxxy_yyzzzz[k] = -g_z_0_xxxxy_yyzzzz[k] * ab_x + g_z_0_xxxxy_xyyzzzz[k];

                g_z_0_xxxxxy_yzzzzz[k] = -g_z_0_xxxxy_yzzzzz[k] * ab_x + g_z_0_xxxxy_xyzzzzz[k];

                g_z_0_xxxxxy_zzzzzz[k] = -g_z_0_xxxxy_zzzzzz[k] * ab_x + g_z_0_xxxxy_xzzzzzz[k];
            }

            /// Set up 1624-1652 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_xxxxxx = cbuffer.data(ii_geom_10_off + 1624 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxxxy = cbuffer.data(ii_geom_10_off + 1625 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxxxz = cbuffer.data(ii_geom_10_off + 1626 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxxyy = cbuffer.data(ii_geom_10_off + 1627 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxxyz = cbuffer.data(ii_geom_10_off + 1628 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxxzz = cbuffer.data(ii_geom_10_off + 1629 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxyyy = cbuffer.data(ii_geom_10_off + 1630 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxyyz = cbuffer.data(ii_geom_10_off + 1631 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxyzz = cbuffer.data(ii_geom_10_off + 1632 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxzzz = cbuffer.data(ii_geom_10_off + 1633 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxyyyy = cbuffer.data(ii_geom_10_off + 1634 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxyyyz = cbuffer.data(ii_geom_10_off + 1635 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxyyzz = cbuffer.data(ii_geom_10_off + 1636 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxyzzz = cbuffer.data(ii_geom_10_off + 1637 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxzzzz = cbuffer.data(ii_geom_10_off + 1638 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyyyyy = cbuffer.data(ii_geom_10_off + 1639 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyyyyz = cbuffer.data(ii_geom_10_off + 1640 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyyyzz = cbuffer.data(ii_geom_10_off + 1641 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyyzzz = cbuffer.data(ii_geom_10_off + 1642 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyzzzz = cbuffer.data(ii_geom_10_off + 1643 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xzzzzz = cbuffer.data(ii_geom_10_off + 1644 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyyyyy = cbuffer.data(ii_geom_10_off + 1645 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyyyyz = cbuffer.data(ii_geom_10_off + 1646 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyyyzz = cbuffer.data(ii_geom_10_off + 1647 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyyzzz = cbuffer.data(ii_geom_10_off + 1648 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyzzzz = cbuffer.data(ii_geom_10_off + 1649 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yzzzzz = cbuffer.data(ii_geom_10_off + 1650 * ccomps * dcomps);

            auto g_z_0_xxxxxz_zzzzzz = cbuffer.data(ii_geom_10_off + 1651 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxxz_xxxxxx, g_z_0_xxxxxz_xxxxxy, g_z_0_xxxxxz_xxxxxz, g_z_0_xxxxxz_xxxxyy, g_z_0_xxxxxz_xxxxyz, g_z_0_xxxxxz_xxxxzz, g_z_0_xxxxxz_xxxyyy, g_z_0_xxxxxz_xxxyyz, g_z_0_xxxxxz_xxxyzz, g_z_0_xxxxxz_xxxzzz, g_z_0_xxxxxz_xxyyyy, g_z_0_xxxxxz_xxyyyz, g_z_0_xxxxxz_xxyyzz, g_z_0_xxxxxz_xxyzzz, g_z_0_xxxxxz_xxzzzz, g_z_0_xxxxxz_xyyyyy, g_z_0_xxxxxz_xyyyyz, g_z_0_xxxxxz_xyyyzz, g_z_0_xxxxxz_xyyzzz, g_z_0_xxxxxz_xyzzzz, g_z_0_xxxxxz_xzzzzz, g_z_0_xxxxxz_yyyyyy, g_z_0_xxxxxz_yyyyyz, g_z_0_xxxxxz_yyyyzz, g_z_0_xxxxxz_yyyzzz, g_z_0_xxxxxz_yyzzzz, g_z_0_xxxxxz_yzzzzz, g_z_0_xxxxxz_zzzzzz, g_z_0_xxxxz_xxxxxx, g_z_0_xxxxz_xxxxxxx, g_z_0_xxxxz_xxxxxxy, g_z_0_xxxxz_xxxxxxz, g_z_0_xxxxz_xxxxxy, g_z_0_xxxxz_xxxxxyy, g_z_0_xxxxz_xxxxxyz, g_z_0_xxxxz_xxxxxz, g_z_0_xxxxz_xxxxxzz, g_z_0_xxxxz_xxxxyy, g_z_0_xxxxz_xxxxyyy, g_z_0_xxxxz_xxxxyyz, g_z_0_xxxxz_xxxxyz, g_z_0_xxxxz_xxxxyzz, g_z_0_xxxxz_xxxxzz, g_z_0_xxxxz_xxxxzzz, g_z_0_xxxxz_xxxyyy, g_z_0_xxxxz_xxxyyyy, g_z_0_xxxxz_xxxyyyz, g_z_0_xxxxz_xxxyyz, g_z_0_xxxxz_xxxyyzz, g_z_0_xxxxz_xxxyzz, g_z_0_xxxxz_xxxyzzz, g_z_0_xxxxz_xxxzzz, g_z_0_xxxxz_xxxzzzz, g_z_0_xxxxz_xxyyyy, g_z_0_xxxxz_xxyyyyy, g_z_0_xxxxz_xxyyyyz, g_z_0_xxxxz_xxyyyz, g_z_0_xxxxz_xxyyyzz, g_z_0_xxxxz_xxyyzz, g_z_0_xxxxz_xxyyzzz, g_z_0_xxxxz_xxyzzz, g_z_0_xxxxz_xxyzzzz, g_z_0_xxxxz_xxzzzz, g_z_0_xxxxz_xxzzzzz, g_z_0_xxxxz_xyyyyy, g_z_0_xxxxz_xyyyyyy, g_z_0_xxxxz_xyyyyyz, g_z_0_xxxxz_xyyyyz, g_z_0_xxxxz_xyyyyzz, g_z_0_xxxxz_xyyyzz, g_z_0_xxxxz_xyyyzzz, g_z_0_xxxxz_xyyzzz, g_z_0_xxxxz_xyyzzzz, g_z_0_xxxxz_xyzzzz, g_z_0_xxxxz_xyzzzzz, g_z_0_xxxxz_xzzzzz, g_z_0_xxxxz_xzzzzzz, g_z_0_xxxxz_yyyyyy, g_z_0_xxxxz_yyyyyz, g_z_0_xxxxz_yyyyzz, g_z_0_xxxxz_yyyzzz, g_z_0_xxxxz_yyzzzz, g_z_0_xxxxz_yzzzzz, g_z_0_xxxxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_xxxxxx[k] = -g_z_0_xxxxz_xxxxxx[k] * ab_x + g_z_0_xxxxz_xxxxxxx[k];

                g_z_0_xxxxxz_xxxxxy[k] = -g_z_0_xxxxz_xxxxxy[k] * ab_x + g_z_0_xxxxz_xxxxxxy[k];

                g_z_0_xxxxxz_xxxxxz[k] = -g_z_0_xxxxz_xxxxxz[k] * ab_x + g_z_0_xxxxz_xxxxxxz[k];

                g_z_0_xxxxxz_xxxxyy[k] = -g_z_0_xxxxz_xxxxyy[k] * ab_x + g_z_0_xxxxz_xxxxxyy[k];

                g_z_0_xxxxxz_xxxxyz[k] = -g_z_0_xxxxz_xxxxyz[k] * ab_x + g_z_0_xxxxz_xxxxxyz[k];

                g_z_0_xxxxxz_xxxxzz[k] = -g_z_0_xxxxz_xxxxzz[k] * ab_x + g_z_0_xxxxz_xxxxxzz[k];

                g_z_0_xxxxxz_xxxyyy[k] = -g_z_0_xxxxz_xxxyyy[k] * ab_x + g_z_0_xxxxz_xxxxyyy[k];

                g_z_0_xxxxxz_xxxyyz[k] = -g_z_0_xxxxz_xxxyyz[k] * ab_x + g_z_0_xxxxz_xxxxyyz[k];

                g_z_0_xxxxxz_xxxyzz[k] = -g_z_0_xxxxz_xxxyzz[k] * ab_x + g_z_0_xxxxz_xxxxyzz[k];

                g_z_0_xxxxxz_xxxzzz[k] = -g_z_0_xxxxz_xxxzzz[k] * ab_x + g_z_0_xxxxz_xxxxzzz[k];

                g_z_0_xxxxxz_xxyyyy[k] = -g_z_0_xxxxz_xxyyyy[k] * ab_x + g_z_0_xxxxz_xxxyyyy[k];

                g_z_0_xxxxxz_xxyyyz[k] = -g_z_0_xxxxz_xxyyyz[k] * ab_x + g_z_0_xxxxz_xxxyyyz[k];

                g_z_0_xxxxxz_xxyyzz[k] = -g_z_0_xxxxz_xxyyzz[k] * ab_x + g_z_0_xxxxz_xxxyyzz[k];

                g_z_0_xxxxxz_xxyzzz[k] = -g_z_0_xxxxz_xxyzzz[k] * ab_x + g_z_0_xxxxz_xxxyzzz[k];

                g_z_0_xxxxxz_xxzzzz[k] = -g_z_0_xxxxz_xxzzzz[k] * ab_x + g_z_0_xxxxz_xxxzzzz[k];

                g_z_0_xxxxxz_xyyyyy[k] = -g_z_0_xxxxz_xyyyyy[k] * ab_x + g_z_0_xxxxz_xxyyyyy[k];

                g_z_0_xxxxxz_xyyyyz[k] = -g_z_0_xxxxz_xyyyyz[k] * ab_x + g_z_0_xxxxz_xxyyyyz[k];

                g_z_0_xxxxxz_xyyyzz[k] = -g_z_0_xxxxz_xyyyzz[k] * ab_x + g_z_0_xxxxz_xxyyyzz[k];

                g_z_0_xxxxxz_xyyzzz[k] = -g_z_0_xxxxz_xyyzzz[k] * ab_x + g_z_0_xxxxz_xxyyzzz[k];

                g_z_0_xxxxxz_xyzzzz[k] = -g_z_0_xxxxz_xyzzzz[k] * ab_x + g_z_0_xxxxz_xxyzzzz[k];

                g_z_0_xxxxxz_xzzzzz[k] = -g_z_0_xxxxz_xzzzzz[k] * ab_x + g_z_0_xxxxz_xxzzzzz[k];

                g_z_0_xxxxxz_yyyyyy[k] = -g_z_0_xxxxz_yyyyyy[k] * ab_x + g_z_0_xxxxz_xyyyyyy[k];

                g_z_0_xxxxxz_yyyyyz[k] = -g_z_0_xxxxz_yyyyyz[k] * ab_x + g_z_0_xxxxz_xyyyyyz[k];

                g_z_0_xxxxxz_yyyyzz[k] = -g_z_0_xxxxz_yyyyzz[k] * ab_x + g_z_0_xxxxz_xyyyyzz[k];

                g_z_0_xxxxxz_yyyzzz[k] = -g_z_0_xxxxz_yyyzzz[k] * ab_x + g_z_0_xxxxz_xyyyzzz[k];

                g_z_0_xxxxxz_yyzzzz[k] = -g_z_0_xxxxz_yyzzzz[k] * ab_x + g_z_0_xxxxz_xyyzzzz[k];

                g_z_0_xxxxxz_yzzzzz[k] = -g_z_0_xxxxz_yzzzzz[k] * ab_x + g_z_0_xxxxz_xyzzzzz[k];

                g_z_0_xxxxxz_zzzzzz[k] = -g_z_0_xxxxz_zzzzzz[k] * ab_x + g_z_0_xxxxz_xzzzzzz[k];
            }

            /// Set up 1652-1680 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_xxxxxx = cbuffer.data(ii_geom_10_off + 1652 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxxxy = cbuffer.data(ii_geom_10_off + 1653 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxxxz = cbuffer.data(ii_geom_10_off + 1654 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxxyy = cbuffer.data(ii_geom_10_off + 1655 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxxyz = cbuffer.data(ii_geom_10_off + 1656 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxxzz = cbuffer.data(ii_geom_10_off + 1657 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxyyy = cbuffer.data(ii_geom_10_off + 1658 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxyyz = cbuffer.data(ii_geom_10_off + 1659 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxyzz = cbuffer.data(ii_geom_10_off + 1660 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxzzz = cbuffer.data(ii_geom_10_off + 1661 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxyyyy = cbuffer.data(ii_geom_10_off + 1662 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxyyyz = cbuffer.data(ii_geom_10_off + 1663 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxyyzz = cbuffer.data(ii_geom_10_off + 1664 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxyzzz = cbuffer.data(ii_geom_10_off + 1665 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxzzzz = cbuffer.data(ii_geom_10_off + 1666 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyyyyy = cbuffer.data(ii_geom_10_off + 1667 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyyyyz = cbuffer.data(ii_geom_10_off + 1668 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyyyzz = cbuffer.data(ii_geom_10_off + 1669 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyyzzz = cbuffer.data(ii_geom_10_off + 1670 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyzzzz = cbuffer.data(ii_geom_10_off + 1671 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xzzzzz = cbuffer.data(ii_geom_10_off + 1672 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyyyyy = cbuffer.data(ii_geom_10_off + 1673 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyyyyz = cbuffer.data(ii_geom_10_off + 1674 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyyyzz = cbuffer.data(ii_geom_10_off + 1675 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyyzzz = cbuffer.data(ii_geom_10_off + 1676 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyzzzz = cbuffer.data(ii_geom_10_off + 1677 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yzzzzz = cbuffer.data(ii_geom_10_off + 1678 * ccomps * dcomps);

            auto g_z_0_xxxxyy_zzzzzz = cbuffer.data(ii_geom_10_off + 1679 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxyy_xxxxxx, g_z_0_xxxxyy_xxxxxy, g_z_0_xxxxyy_xxxxxz, g_z_0_xxxxyy_xxxxyy, g_z_0_xxxxyy_xxxxyz, g_z_0_xxxxyy_xxxxzz, g_z_0_xxxxyy_xxxyyy, g_z_0_xxxxyy_xxxyyz, g_z_0_xxxxyy_xxxyzz, g_z_0_xxxxyy_xxxzzz, g_z_0_xxxxyy_xxyyyy, g_z_0_xxxxyy_xxyyyz, g_z_0_xxxxyy_xxyyzz, g_z_0_xxxxyy_xxyzzz, g_z_0_xxxxyy_xxzzzz, g_z_0_xxxxyy_xyyyyy, g_z_0_xxxxyy_xyyyyz, g_z_0_xxxxyy_xyyyzz, g_z_0_xxxxyy_xyyzzz, g_z_0_xxxxyy_xyzzzz, g_z_0_xxxxyy_xzzzzz, g_z_0_xxxxyy_yyyyyy, g_z_0_xxxxyy_yyyyyz, g_z_0_xxxxyy_yyyyzz, g_z_0_xxxxyy_yyyzzz, g_z_0_xxxxyy_yyzzzz, g_z_0_xxxxyy_yzzzzz, g_z_0_xxxxyy_zzzzzz, g_z_0_xxxyy_xxxxxx, g_z_0_xxxyy_xxxxxxx, g_z_0_xxxyy_xxxxxxy, g_z_0_xxxyy_xxxxxxz, g_z_0_xxxyy_xxxxxy, g_z_0_xxxyy_xxxxxyy, g_z_0_xxxyy_xxxxxyz, g_z_0_xxxyy_xxxxxz, g_z_0_xxxyy_xxxxxzz, g_z_0_xxxyy_xxxxyy, g_z_0_xxxyy_xxxxyyy, g_z_0_xxxyy_xxxxyyz, g_z_0_xxxyy_xxxxyz, g_z_0_xxxyy_xxxxyzz, g_z_0_xxxyy_xxxxzz, g_z_0_xxxyy_xxxxzzz, g_z_0_xxxyy_xxxyyy, g_z_0_xxxyy_xxxyyyy, g_z_0_xxxyy_xxxyyyz, g_z_0_xxxyy_xxxyyz, g_z_0_xxxyy_xxxyyzz, g_z_0_xxxyy_xxxyzz, g_z_0_xxxyy_xxxyzzz, g_z_0_xxxyy_xxxzzz, g_z_0_xxxyy_xxxzzzz, g_z_0_xxxyy_xxyyyy, g_z_0_xxxyy_xxyyyyy, g_z_0_xxxyy_xxyyyyz, g_z_0_xxxyy_xxyyyz, g_z_0_xxxyy_xxyyyzz, g_z_0_xxxyy_xxyyzz, g_z_0_xxxyy_xxyyzzz, g_z_0_xxxyy_xxyzzz, g_z_0_xxxyy_xxyzzzz, g_z_0_xxxyy_xxzzzz, g_z_0_xxxyy_xxzzzzz, g_z_0_xxxyy_xyyyyy, g_z_0_xxxyy_xyyyyyy, g_z_0_xxxyy_xyyyyyz, g_z_0_xxxyy_xyyyyz, g_z_0_xxxyy_xyyyyzz, g_z_0_xxxyy_xyyyzz, g_z_0_xxxyy_xyyyzzz, g_z_0_xxxyy_xyyzzz, g_z_0_xxxyy_xyyzzzz, g_z_0_xxxyy_xyzzzz, g_z_0_xxxyy_xyzzzzz, g_z_0_xxxyy_xzzzzz, g_z_0_xxxyy_xzzzzzz, g_z_0_xxxyy_yyyyyy, g_z_0_xxxyy_yyyyyz, g_z_0_xxxyy_yyyyzz, g_z_0_xxxyy_yyyzzz, g_z_0_xxxyy_yyzzzz, g_z_0_xxxyy_yzzzzz, g_z_0_xxxyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_xxxxxx[k] = -g_z_0_xxxyy_xxxxxx[k] * ab_x + g_z_0_xxxyy_xxxxxxx[k];

                g_z_0_xxxxyy_xxxxxy[k] = -g_z_0_xxxyy_xxxxxy[k] * ab_x + g_z_0_xxxyy_xxxxxxy[k];

                g_z_0_xxxxyy_xxxxxz[k] = -g_z_0_xxxyy_xxxxxz[k] * ab_x + g_z_0_xxxyy_xxxxxxz[k];

                g_z_0_xxxxyy_xxxxyy[k] = -g_z_0_xxxyy_xxxxyy[k] * ab_x + g_z_0_xxxyy_xxxxxyy[k];

                g_z_0_xxxxyy_xxxxyz[k] = -g_z_0_xxxyy_xxxxyz[k] * ab_x + g_z_0_xxxyy_xxxxxyz[k];

                g_z_0_xxxxyy_xxxxzz[k] = -g_z_0_xxxyy_xxxxzz[k] * ab_x + g_z_0_xxxyy_xxxxxzz[k];

                g_z_0_xxxxyy_xxxyyy[k] = -g_z_0_xxxyy_xxxyyy[k] * ab_x + g_z_0_xxxyy_xxxxyyy[k];

                g_z_0_xxxxyy_xxxyyz[k] = -g_z_0_xxxyy_xxxyyz[k] * ab_x + g_z_0_xxxyy_xxxxyyz[k];

                g_z_0_xxxxyy_xxxyzz[k] = -g_z_0_xxxyy_xxxyzz[k] * ab_x + g_z_0_xxxyy_xxxxyzz[k];

                g_z_0_xxxxyy_xxxzzz[k] = -g_z_0_xxxyy_xxxzzz[k] * ab_x + g_z_0_xxxyy_xxxxzzz[k];

                g_z_0_xxxxyy_xxyyyy[k] = -g_z_0_xxxyy_xxyyyy[k] * ab_x + g_z_0_xxxyy_xxxyyyy[k];

                g_z_0_xxxxyy_xxyyyz[k] = -g_z_0_xxxyy_xxyyyz[k] * ab_x + g_z_0_xxxyy_xxxyyyz[k];

                g_z_0_xxxxyy_xxyyzz[k] = -g_z_0_xxxyy_xxyyzz[k] * ab_x + g_z_0_xxxyy_xxxyyzz[k];

                g_z_0_xxxxyy_xxyzzz[k] = -g_z_0_xxxyy_xxyzzz[k] * ab_x + g_z_0_xxxyy_xxxyzzz[k];

                g_z_0_xxxxyy_xxzzzz[k] = -g_z_0_xxxyy_xxzzzz[k] * ab_x + g_z_0_xxxyy_xxxzzzz[k];

                g_z_0_xxxxyy_xyyyyy[k] = -g_z_0_xxxyy_xyyyyy[k] * ab_x + g_z_0_xxxyy_xxyyyyy[k];

                g_z_0_xxxxyy_xyyyyz[k] = -g_z_0_xxxyy_xyyyyz[k] * ab_x + g_z_0_xxxyy_xxyyyyz[k];

                g_z_0_xxxxyy_xyyyzz[k] = -g_z_0_xxxyy_xyyyzz[k] * ab_x + g_z_0_xxxyy_xxyyyzz[k];

                g_z_0_xxxxyy_xyyzzz[k] = -g_z_0_xxxyy_xyyzzz[k] * ab_x + g_z_0_xxxyy_xxyyzzz[k];

                g_z_0_xxxxyy_xyzzzz[k] = -g_z_0_xxxyy_xyzzzz[k] * ab_x + g_z_0_xxxyy_xxyzzzz[k];

                g_z_0_xxxxyy_xzzzzz[k] = -g_z_0_xxxyy_xzzzzz[k] * ab_x + g_z_0_xxxyy_xxzzzzz[k];

                g_z_0_xxxxyy_yyyyyy[k] = -g_z_0_xxxyy_yyyyyy[k] * ab_x + g_z_0_xxxyy_xyyyyyy[k];

                g_z_0_xxxxyy_yyyyyz[k] = -g_z_0_xxxyy_yyyyyz[k] * ab_x + g_z_0_xxxyy_xyyyyyz[k];

                g_z_0_xxxxyy_yyyyzz[k] = -g_z_0_xxxyy_yyyyzz[k] * ab_x + g_z_0_xxxyy_xyyyyzz[k];

                g_z_0_xxxxyy_yyyzzz[k] = -g_z_0_xxxyy_yyyzzz[k] * ab_x + g_z_0_xxxyy_xyyyzzz[k];

                g_z_0_xxxxyy_yyzzzz[k] = -g_z_0_xxxyy_yyzzzz[k] * ab_x + g_z_0_xxxyy_xyyzzzz[k];

                g_z_0_xxxxyy_yzzzzz[k] = -g_z_0_xxxyy_yzzzzz[k] * ab_x + g_z_0_xxxyy_xyzzzzz[k];

                g_z_0_xxxxyy_zzzzzz[k] = -g_z_0_xxxyy_zzzzzz[k] * ab_x + g_z_0_xxxyy_xzzzzzz[k];
            }

            /// Set up 1680-1708 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_xxxxxx = cbuffer.data(ii_geom_10_off + 1680 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxxxy = cbuffer.data(ii_geom_10_off + 1681 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxxxz = cbuffer.data(ii_geom_10_off + 1682 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxxyy = cbuffer.data(ii_geom_10_off + 1683 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxxyz = cbuffer.data(ii_geom_10_off + 1684 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxxzz = cbuffer.data(ii_geom_10_off + 1685 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxyyy = cbuffer.data(ii_geom_10_off + 1686 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxyyz = cbuffer.data(ii_geom_10_off + 1687 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxyzz = cbuffer.data(ii_geom_10_off + 1688 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxzzz = cbuffer.data(ii_geom_10_off + 1689 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxyyyy = cbuffer.data(ii_geom_10_off + 1690 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxyyyz = cbuffer.data(ii_geom_10_off + 1691 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxyyzz = cbuffer.data(ii_geom_10_off + 1692 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxyzzz = cbuffer.data(ii_geom_10_off + 1693 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxzzzz = cbuffer.data(ii_geom_10_off + 1694 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyyyyy = cbuffer.data(ii_geom_10_off + 1695 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyyyyz = cbuffer.data(ii_geom_10_off + 1696 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyyyzz = cbuffer.data(ii_geom_10_off + 1697 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyyzzz = cbuffer.data(ii_geom_10_off + 1698 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyzzzz = cbuffer.data(ii_geom_10_off + 1699 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xzzzzz = cbuffer.data(ii_geom_10_off + 1700 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyyyyy = cbuffer.data(ii_geom_10_off + 1701 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyyyyz = cbuffer.data(ii_geom_10_off + 1702 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyyyzz = cbuffer.data(ii_geom_10_off + 1703 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyyzzz = cbuffer.data(ii_geom_10_off + 1704 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyzzzz = cbuffer.data(ii_geom_10_off + 1705 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yzzzzz = cbuffer.data(ii_geom_10_off + 1706 * ccomps * dcomps);

            auto g_z_0_xxxxyz_zzzzzz = cbuffer.data(ii_geom_10_off + 1707 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxyz_xxxxxx, g_z_0_xxxxyz_xxxxxy, g_z_0_xxxxyz_xxxxxz, g_z_0_xxxxyz_xxxxyy, g_z_0_xxxxyz_xxxxyz, g_z_0_xxxxyz_xxxxzz, g_z_0_xxxxyz_xxxyyy, g_z_0_xxxxyz_xxxyyz, g_z_0_xxxxyz_xxxyzz, g_z_0_xxxxyz_xxxzzz, g_z_0_xxxxyz_xxyyyy, g_z_0_xxxxyz_xxyyyz, g_z_0_xxxxyz_xxyyzz, g_z_0_xxxxyz_xxyzzz, g_z_0_xxxxyz_xxzzzz, g_z_0_xxxxyz_xyyyyy, g_z_0_xxxxyz_xyyyyz, g_z_0_xxxxyz_xyyyzz, g_z_0_xxxxyz_xyyzzz, g_z_0_xxxxyz_xyzzzz, g_z_0_xxxxyz_xzzzzz, g_z_0_xxxxyz_yyyyyy, g_z_0_xxxxyz_yyyyyz, g_z_0_xxxxyz_yyyyzz, g_z_0_xxxxyz_yyyzzz, g_z_0_xxxxyz_yyzzzz, g_z_0_xxxxyz_yzzzzz, g_z_0_xxxxyz_zzzzzz, g_z_0_xxxyz_xxxxxx, g_z_0_xxxyz_xxxxxxx, g_z_0_xxxyz_xxxxxxy, g_z_0_xxxyz_xxxxxxz, g_z_0_xxxyz_xxxxxy, g_z_0_xxxyz_xxxxxyy, g_z_0_xxxyz_xxxxxyz, g_z_0_xxxyz_xxxxxz, g_z_0_xxxyz_xxxxxzz, g_z_0_xxxyz_xxxxyy, g_z_0_xxxyz_xxxxyyy, g_z_0_xxxyz_xxxxyyz, g_z_0_xxxyz_xxxxyz, g_z_0_xxxyz_xxxxyzz, g_z_0_xxxyz_xxxxzz, g_z_0_xxxyz_xxxxzzz, g_z_0_xxxyz_xxxyyy, g_z_0_xxxyz_xxxyyyy, g_z_0_xxxyz_xxxyyyz, g_z_0_xxxyz_xxxyyz, g_z_0_xxxyz_xxxyyzz, g_z_0_xxxyz_xxxyzz, g_z_0_xxxyz_xxxyzzz, g_z_0_xxxyz_xxxzzz, g_z_0_xxxyz_xxxzzzz, g_z_0_xxxyz_xxyyyy, g_z_0_xxxyz_xxyyyyy, g_z_0_xxxyz_xxyyyyz, g_z_0_xxxyz_xxyyyz, g_z_0_xxxyz_xxyyyzz, g_z_0_xxxyz_xxyyzz, g_z_0_xxxyz_xxyyzzz, g_z_0_xxxyz_xxyzzz, g_z_0_xxxyz_xxyzzzz, g_z_0_xxxyz_xxzzzz, g_z_0_xxxyz_xxzzzzz, g_z_0_xxxyz_xyyyyy, g_z_0_xxxyz_xyyyyyy, g_z_0_xxxyz_xyyyyyz, g_z_0_xxxyz_xyyyyz, g_z_0_xxxyz_xyyyyzz, g_z_0_xxxyz_xyyyzz, g_z_0_xxxyz_xyyyzzz, g_z_0_xxxyz_xyyzzz, g_z_0_xxxyz_xyyzzzz, g_z_0_xxxyz_xyzzzz, g_z_0_xxxyz_xyzzzzz, g_z_0_xxxyz_xzzzzz, g_z_0_xxxyz_xzzzzzz, g_z_0_xxxyz_yyyyyy, g_z_0_xxxyz_yyyyyz, g_z_0_xxxyz_yyyyzz, g_z_0_xxxyz_yyyzzz, g_z_0_xxxyz_yyzzzz, g_z_0_xxxyz_yzzzzz, g_z_0_xxxyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_xxxxxx[k] = -g_z_0_xxxyz_xxxxxx[k] * ab_x + g_z_0_xxxyz_xxxxxxx[k];

                g_z_0_xxxxyz_xxxxxy[k] = -g_z_0_xxxyz_xxxxxy[k] * ab_x + g_z_0_xxxyz_xxxxxxy[k];

                g_z_0_xxxxyz_xxxxxz[k] = -g_z_0_xxxyz_xxxxxz[k] * ab_x + g_z_0_xxxyz_xxxxxxz[k];

                g_z_0_xxxxyz_xxxxyy[k] = -g_z_0_xxxyz_xxxxyy[k] * ab_x + g_z_0_xxxyz_xxxxxyy[k];

                g_z_0_xxxxyz_xxxxyz[k] = -g_z_0_xxxyz_xxxxyz[k] * ab_x + g_z_0_xxxyz_xxxxxyz[k];

                g_z_0_xxxxyz_xxxxzz[k] = -g_z_0_xxxyz_xxxxzz[k] * ab_x + g_z_0_xxxyz_xxxxxzz[k];

                g_z_0_xxxxyz_xxxyyy[k] = -g_z_0_xxxyz_xxxyyy[k] * ab_x + g_z_0_xxxyz_xxxxyyy[k];

                g_z_0_xxxxyz_xxxyyz[k] = -g_z_0_xxxyz_xxxyyz[k] * ab_x + g_z_0_xxxyz_xxxxyyz[k];

                g_z_0_xxxxyz_xxxyzz[k] = -g_z_0_xxxyz_xxxyzz[k] * ab_x + g_z_0_xxxyz_xxxxyzz[k];

                g_z_0_xxxxyz_xxxzzz[k] = -g_z_0_xxxyz_xxxzzz[k] * ab_x + g_z_0_xxxyz_xxxxzzz[k];

                g_z_0_xxxxyz_xxyyyy[k] = -g_z_0_xxxyz_xxyyyy[k] * ab_x + g_z_0_xxxyz_xxxyyyy[k];

                g_z_0_xxxxyz_xxyyyz[k] = -g_z_0_xxxyz_xxyyyz[k] * ab_x + g_z_0_xxxyz_xxxyyyz[k];

                g_z_0_xxxxyz_xxyyzz[k] = -g_z_0_xxxyz_xxyyzz[k] * ab_x + g_z_0_xxxyz_xxxyyzz[k];

                g_z_0_xxxxyz_xxyzzz[k] = -g_z_0_xxxyz_xxyzzz[k] * ab_x + g_z_0_xxxyz_xxxyzzz[k];

                g_z_0_xxxxyz_xxzzzz[k] = -g_z_0_xxxyz_xxzzzz[k] * ab_x + g_z_0_xxxyz_xxxzzzz[k];

                g_z_0_xxxxyz_xyyyyy[k] = -g_z_0_xxxyz_xyyyyy[k] * ab_x + g_z_0_xxxyz_xxyyyyy[k];

                g_z_0_xxxxyz_xyyyyz[k] = -g_z_0_xxxyz_xyyyyz[k] * ab_x + g_z_0_xxxyz_xxyyyyz[k];

                g_z_0_xxxxyz_xyyyzz[k] = -g_z_0_xxxyz_xyyyzz[k] * ab_x + g_z_0_xxxyz_xxyyyzz[k];

                g_z_0_xxxxyz_xyyzzz[k] = -g_z_0_xxxyz_xyyzzz[k] * ab_x + g_z_0_xxxyz_xxyyzzz[k];

                g_z_0_xxxxyz_xyzzzz[k] = -g_z_0_xxxyz_xyzzzz[k] * ab_x + g_z_0_xxxyz_xxyzzzz[k];

                g_z_0_xxxxyz_xzzzzz[k] = -g_z_0_xxxyz_xzzzzz[k] * ab_x + g_z_0_xxxyz_xxzzzzz[k];

                g_z_0_xxxxyz_yyyyyy[k] = -g_z_0_xxxyz_yyyyyy[k] * ab_x + g_z_0_xxxyz_xyyyyyy[k];

                g_z_0_xxxxyz_yyyyyz[k] = -g_z_0_xxxyz_yyyyyz[k] * ab_x + g_z_0_xxxyz_xyyyyyz[k];

                g_z_0_xxxxyz_yyyyzz[k] = -g_z_0_xxxyz_yyyyzz[k] * ab_x + g_z_0_xxxyz_xyyyyzz[k];

                g_z_0_xxxxyz_yyyzzz[k] = -g_z_0_xxxyz_yyyzzz[k] * ab_x + g_z_0_xxxyz_xyyyzzz[k];

                g_z_0_xxxxyz_yyzzzz[k] = -g_z_0_xxxyz_yyzzzz[k] * ab_x + g_z_0_xxxyz_xyyzzzz[k];

                g_z_0_xxxxyz_yzzzzz[k] = -g_z_0_xxxyz_yzzzzz[k] * ab_x + g_z_0_xxxyz_xyzzzzz[k];

                g_z_0_xxxxyz_zzzzzz[k] = -g_z_0_xxxyz_zzzzzz[k] * ab_x + g_z_0_xxxyz_xzzzzzz[k];
            }

            /// Set up 1708-1736 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1708 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1709 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1710 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1711 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1712 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1713 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1714 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1715 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1716 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1717 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1718 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1719 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1720 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1721 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1722 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1723 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1724 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1725 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1726 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1727 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1728 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1729 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1730 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1731 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1732 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1733 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1734 * ccomps * dcomps);

            auto g_z_0_xxxxzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1735 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxzz_xxxxxx, g_z_0_xxxxzz_xxxxxy, g_z_0_xxxxzz_xxxxxz, g_z_0_xxxxzz_xxxxyy, g_z_0_xxxxzz_xxxxyz, g_z_0_xxxxzz_xxxxzz, g_z_0_xxxxzz_xxxyyy, g_z_0_xxxxzz_xxxyyz, g_z_0_xxxxzz_xxxyzz, g_z_0_xxxxzz_xxxzzz, g_z_0_xxxxzz_xxyyyy, g_z_0_xxxxzz_xxyyyz, g_z_0_xxxxzz_xxyyzz, g_z_0_xxxxzz_xxyzzz, g_z_0_xxxxzz_xxzzzz, g_z_0_xxxxzz_xyyyyy, g_z_0_xxxxzz_xyyyyz, g_z_0_xxxxzz_xyyyzz, g_z_0_xxxxzz_xyyzzz, g_z_0_xxxxzz_xyzzzz, g_z_0_xxxxzz_xzzzzz, g_z_0_xxxxzz_yyyyyy, g_z_0_xxxxzz_yyyyyz, g_z_0_xxxxzz_yyyyzz, g_z_0_xxxxzz_yyyzzz, g_z_0_xxxxzz_yyzzzz, g_z_0_xxxxzz_yzzzzz, g_z_0_xxxxzz_zzzzzz, g_z_0_xxxzz_xxxxxx, g_z_0_xxxzz_xxxxxxx, g_z_0_xxxzz_xxxxxxy, g_z_0_xxxzz_xxxxxxz, g_z_0_xxxzz_xxxxxy, g_z_0_xxxzz_xxxxxyy, g_z_0_xxxzz_xxxxxyz, g_z_0_xxxzz_xxxxxz, g_z_0_xxxzz_xxxxxzz, g_z_0_xxxzz_xxxxyy, g_z_0_xxxzz_xxxxyyy, g_z_0_xxxzz_xxxxyyz, g_z_0_xxxzz_xxxxyz, g_z_0_xxxzz_xxxxyzz, g_z_0_xxxzz_xxxxzz, g_z_0_xxxzz_xxxxzzz, g_z_0_xxxzz_xxxyyy, g_z_0_xxxzz_xxxyyyy, g_z_0_xxxzz_xxxyyyz, g_z_0_xxxzz_xxxyyz, g_z_0_xxxzz_xxxyyzz, g_z_0_xxxzz_xxxyzz, g_z_0_xxxzz_xxxyzzz, g_z_0_xxxzz_xxxzzz, g_z_0_xxxzz_xxxzzzz, g_z_0_xxxzz_xxyyyy, g_z_0_xxxzz_xxyyyyy, g_z_0_xxxzz_xxyyyyz, g_z_0_xxxzz_xxyyyz, g_z_0_xxxzz_xxyyyzz, g_z_0_xxxzz_xxyyzz, g_z_0_xxxzz_xxyyzzz, g_z_0_xxxzz_xxyzzz, g_z_0_xxxzz_xxyzzzz, g_z_0_xxxzz_xxzzzz, g_z_0_xxxzz_xxzzzzz, g_z_0_xxxzz_xyyyyy, g_z_0_xxxzz_xyyyyyy, g_z_0_xxxzz_xyyyyyz, g_z_0_xxxzz_xyyyyz, g_z_0_xxxzz_xyyyyzz, g_z_0_xxxzz_xyyyzz, g_z_0_xxxzz_xyyyzzz, g_z_0_xxxzz_xyyzzz, g_z_0_xxxzz_xyyzzzz, g_z_0_xxxzz_xyzzzz, g_z_0_xxxzz_xyzzzzz, g_z_0_xxxzz_xzzzzz, g_z_0_xxxzz_xzzzzzz, g_z_0_xxxzz_yyyyyy, g_z_0_xxxzz_yyyyyz, g_z_0_xxxzz_yyyyzz, g_z_0_xxxzz_yyyzzz, g_z_0_xxxzz_yyzzzz, g_z_0_xxxzz_yzzzzz, g_z_0_xxxzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_xxxxxx[k] = -g_z_0_xxxzz_xxxxxx[k] * ab_x + g_z_0_xxxzz_xxxxxxx[k];

                g_z_0_xxxxzz_xxxxxy[k] = -g_z_0_xxxzz_xxxxxy[k] * ab_x + g_z_0_xxxzz_xxxxxxy[k];

                g_z_0_xxxxzz_xxxxxz[k] = -g_z_0_xxxzz_xxxxxz[k] * ab_x + g_z_0_xxxzz_xxxxxxz[k];

                g_z_0_xxxxzz_xxxxyy[k] = -g_z_0_xxxzz_xxxxyy[k] * ab_x + g_z_0_xxxzz_xxxxxyy[k];

                g_z_0_xxxxzz_xxxxyz[k] = -g_z_0_xxxzz_xxxxyz[k] * ab_x + g_z_0_xxxzz_xxxxxyz[k];

                g_z_0_xxxxzz_xxxxzz[k] = -g_z_0_xxxzz_xxxxzz[k] * ab_x + g_z_0_xxxzz_xxxxxzz[k];

                g_z_0_xxxxzz_xxxyyy[k] = -g_z_0_xxxzz_xxxyyy[k] * ab_x + g_z_0_xxxzz_xxxxyyy[k];

                g_z_0_xxxxzz_xxxyyz[k] = -g_z_0_xxxzz_xxxyyz[k] * ab_x + g_z_0_xxxzz_xxxxyyz[k];

                g_z_0_xxxxzz_xxxyzz[k] = -g_z_0_xxxzz_xxxyzz[k] * ab_x + g_z_0_xxxzz_xxxxyzz[k];

                g_z_0_xxxxzz_xxxzzz[k] = -g_z_0_xxxzz_xxxzzz[k] * ab_x + g_z_0_xxxzz_xxxxzzz[k];

                g_z_0_xxxxzz_xxyyyy[k] = -g_z_0_xxxzz_xxyyyy[k] * ab_x + g_z_0_xxxzz_xxxyyyy[k];

                g_z_0_xxxxzz_xxyyyz[k] = -g_z_0_xxxzz_xxyyyz[k] * ab_x + g_z_0_xxxzz_xxxyyyz[k];

                g_z_0_xxxxzz_xxyyzz[k] = -g_z_0_xxxzz_xxyyzz[k] * ab_x + g_z_0_xxxzz_xxxyyzz[k];

                g_z_0_xxxxzz_xxyzzz[k] = -g_z_0_xxxzz_xxyzzz[k] * ab_x + g_z_0_xxxzz_xxxyzzz[k];

                g_z_0_xxxxzz_xxzzzz[k] = -g_z_0_xxxzz_xxzzzz[k] * ab_x + g_z_0_xxxzz_xxxzzzz[k];

                g_z_0_xxxxzz_xyyyyy[k] = -g_z_0_xxxzz_xyyyyy[k] * ab_x + g_z_0_xxxzz_xxyyyyy[k];

                g_z_0_xxxxzz_xyyyyz[k] = -g_z_0_xxxzz_xyyyyz[k] * ab_x + g_z_0_xxxzz_xxyyyyz[k];

                g_z_0_xxxxzz_xyyyzz[k] = -g_z_0_xxxzz_xyyyzz[k] * ab_x + g_z_0_xxxzz_xxyyyzz[k];

                g_z_0_xxxxzz_xyyzzz[k] = -g_z_0_xxxzz_xyyzzz[k] * ab_x + g_z_0_xxxzz_xxyyzzz[k];

                g_z_0_xxxxzz_xyzzzz[k] = -g_z_0_xxxzz_xyzzzz[k] * ab_x + g_z_0_xxxzz_xxyzzzz[k];

                g_z_0_xxxxzz_xzzzzz[k] = -g_z_0_xxxzz_xzzzzz[k] * ab_x + g_z_0_xxxzz_xxzzzzz[k];

                g_z_0_xxxxzz_yyyyyy[k] = -g_z_0_xxxzz_yyyyyy[k] * ab_x + g_z_0_xxxzz_xyyyyyy[k];

                g_z_0_xxxxzz_yyyyyz[k] = -g_z_0_xxxzz_yyyyyz[k] * ab_x + g_z_0_xxxzz_xyyyyyz[k];

                g_z_0_xxxxzz_yyyyzz[k] = -g_z_0_xxxzz_yyyyzz[k] * ab_x + g_z_0_xxxzz_xyyyyzz[k];

                g_z_0_xxxxzz_yyyzzz[k] = -g_z_0_xxxzz_yyyzzz[k] * ab_x + g_z_0_xxxzz_xyyyzzz[k];

                g_z_0_xxxxzz_yyzzzz[k] = -g_z_0_xxxzz_yyzzzz[k] * ab_x + g_z_0_xxxzz_xyyzzzz[k];

                g_z_0_xxxxzz_yzzzzz[k] = -g_z_0_xxxzz_yzzzzz[k] * ab_x + g_z_0_xxxzz_xyzzzzz[k];

                g_z_0_xxxxzz_zzzzzz[k] = -g_z_0_xxxzz_zzzzzz[k] * ab_x + g_z_0_xxxzz_xzzzzzz[k];
            }

            /// Set up 1736-1764 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 1736 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 1737 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 1738 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 1739 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 1740 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 1741 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 1742 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 1743 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 1744 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 1745 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 1746 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 1747 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 1748 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 1749 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 1750 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 1751 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 1752 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 1753 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 1754 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 1755 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 1756 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 1757 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 1758 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 1759 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 1760 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 1761 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 1762 * ccomps * dcomps);

            auto g_z_0_xxxyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 1763 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyyy_xxxxxx, g_z_0_xxxyyy_xxxxxy, g_z_0_xxxyyy_xxxxxz, g_z_0_xxxyyy_xxxxyy, g_z_0_xxxyyy_xxxxyz, g_z_0_xxxyyy_xxxxzz, g_z_0_xxxyyy_xxxyyy, g_z_0_xxxyyy_xxxyyz, g_z_0_xxxyyy_xxxyzz, g_z_0_xxxyyy_xxxzzz, g_z_0_xxxyyy_xxyyyy, g_z_0_xxxyyy_xxyyyz, g_z_0_xxxyyy_xxyyzz, g_z_0_xxxyyy_xxyzzz, g_z_0_xxxyyy_xxzzzz, g_z_0_xxxyyy_xyyyyy, g_z_0_xxxyyy_xyyyyz, g_z_0_xxxyyy_xyyyzz, g_z_0_xxxyyy_xyyzzz, g_z_0_xxxyyy_xyzzzz, g_z_0_xxxyyy_xzzzzz, g_z_0_xxxyyy_yyyyyy, g_z_0_xxxyyy_yyyyyz, g_z_0_xxxyyy_yyyyzz, g_z_0_xxxyyy_yyyzzz, g_z_0_xxxyyy_yyzzzz, g_z_0_xxxyyy_yzzzzz, g_z_0_xxxyyy_zzzzzz, g_z_0_xxyyy_xxxxxx, g_z_0_xxyyy_xxxxxxx, g_z_0_xxyyy_xxxxxxy, g_z_0_xxyyy_xxxxxxz, g_z_0_xxyyy_xxxxxy, g_z_0_xxyyy_xxxxxyy, g_z_0_xxyyy_xxxxxyz, g_z_0_xxyyy_xxxxxz, g_z_0_xxyyy_xxxxxzz, g_z_0_xxyyy_xxxxyy, g_z_0_xxyyy_xxxxyyy, g_z_0_xxyyy_xxxxyyz, g_z_0_xxyyy_xxxxyz, g_z_0_xxyyy_xxxxyzz, g_z_0_xxyyy_xxxxzz, g_z_0_xxyyy_xxxxzzz, g_z_0_xxyyy_xxxyyy, g_z_0_xxyyy_xxxyyyy, g_z_0_xxyyy_xxxyyyz, g_z_0_xxyyy_xxxyyz, g_z_0_xxyyy_xxxyyzz, g_z_0_xxyyy_xxxyzz, g_z_0_xxyyy_xxxyzzz, g_z_0_xxyyy_xxxzzz, g_z_0_xxyyy_xxxzzzz, g_z_0_xxyyy_xxyyyy, g_z_0_xxyyy_xxyyyyy, g_z_0_xxyyy_xxyyyyz, g_z_0_xxyyy_xxyyyz, g_z_0_xxyyy_xxyyyzz, g_z_0_xxyyy_xxyyzz, g_z_0_xxyyy_xxyyzzz, g_z_0_xxyyy_xxyzzz, g_z_0_xxyyy_xxyzzzz, g_z_0_xxyyy_xxzzzz, g_z_0_xxyyy_xxzzzzz, g_z_0_xxyyy_xyyyyy, g_z_0_xxyyy_xyyyyyy, g_z_0_xxyyy_xyyyyyz, g_z_0_xxyyy_xyyyyz, g_z_0_xxyyy_xyyyyzz, g_z_0_xxyyy_xyyyzz, g_z_0_xxyyy_xyyyzzz, g_z_0_xxyyy_xyyzzz, g_z_0_xxyyy_xyyzzzz, g_z_0_xxyyy_xyzzzz, g_z_0_xxyyy_xyzzzzz, g_z_0_xxyyy_xzzzzz, g_z_0_xxyyy_xzzzzzz, g_z_0_xxyyy_yyyyyy, g_z_0_xxyyy_yyyyyz, g_z_0_xxyyy_yyyyzz, g_z_0_xxyyy_yyyzzz, g_z_0_xxyyy_yyzzzz, g_z_0_xxyyy_yzzzzz, g_z_0_xxyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_xxxxxx[k] = -g_z_0_xxyyy_xxxxxx[k] * ab_x + g_z_0_xxyyy_xxxxxxx[k];

                g_z_0_xxxyyy_xxxxxy[k] = -g_z_0_xxyyy_xxxxxy[k] * ab_x + g_z_0_xxyyy_xxxxxxy[k];

                g_z_0_xxxyyy_xxxxxz[k] = -g_z_0_xxyyy_xxxxxz[k] * ab_x + g_z_0_xxyyy_xxxxxxz[k];

                g_z_0_xxxyyy_xxxxyy[k] = -g_z_0_xxyyy_xxxxyy[k] * ab_x + g_z_0_xxyyy_xxxxxyy[k];

                g_z_0_xxxyyy_xxxxyz[k] = -g_z_0_xxyyy_xxxxyz[k] * ab_x + g_z_0_xxyyy_xxxxxyz[k];

                g_z_0_xxxyyy_xxxxzz[k] = -g_z_0_xxyyy_xxxxzz[k] * ab_x + g_z_0_xxyyy_xxxxxzz[k];

                g_z_0_xxxyyy_xxxyyy[k] = -g_z_0_xxyyy_xxxyyy[k] * ab_x + g_z_0_xxyyy_xxxxyyy[k];

                g_z_0_xxxyyy_xxxyyz[k] = -g_z_0_xxyyy_xxxyyz[k] * ab_x + g_z_0_xxyyy_xxxxyyz[k];

                g_z_0_xxxyyy_xxxyzz[k] = -g_z_0_xxyyy_xxxyzz[k] * ab_x + g_z_0_xxyyy_xxxxyzz[k];

                g_z_0_xxxyyy_xxxzzz[k] = -g_z_0_xxyyy_xxxzzz[k] * ab_x + g_z_0_xxyyy_xxxxzzz[k];

                g_z_0_xxxyyy_xxyyyy[k] = -g_z_0_xxyyy_xxyyyy[k] * ab_x + g_z_0_xxyyy_xxxyyyy[k];

                g_z_0_xxxyyy_xxyyyz[k] = -g_z_0_xxyyy_xxyyyz[k] * ab_x + g_z_0_xxyyy_xxxyyyz[k];

                g_z_0_xxxyyy_xxyyzz[k] = -g_z_0_xxyyy_xxyyzz[k] * ab_x + g_z_0_xxyyy_xxxyyzz[k];

                g_z_0_xxxyyy_xxyzzz[k] = -g_z_0_xxyyy_xxyzzz[k] * ab_x + g_z_0_xxyyy_xxxyzzz[k];

                g_z_0_xxxyyy_xxzzzz[k] = -g_z_0_xxyyy_xxzzzz[k] * ab_x + g_z_0_xxyyy_xxxzzzz[k];

                g_z_0_xxxyyy_xyyyyy[k] = -g_z_0_xxyyy_xyyyyy[k] * ab_x + g_z_0_xxyyy_xxyyyyy[k];

                g_z_0_xxxyyy_xyyyyz[k] = -g_z_0_xxyyy_xyyyyz[k] * ab_x + g_z_0_xxyyy_xxyyyyz[k];

                g_z_0_xxxyyy_xyyyzz[k] = -g_z_0_xxyyy_xyyyzz[k] * ab_x + g_z_0_xxyyy_xxyyyzz[k];

                g_z_0_xxxyyy_xyyzzz[k] = -g_z_0_xxyyy_xyyzzz[k] * ab_x + g_z_0_xxyyy_xxyyzzz[k];

                g_z_0_xxxyyy_xyzzzz[k] = -g_z_0_xxyyy_xyzzzz[k] * ab_x + g_z_0_xxyyy_xxyzzzz[k];

                g_z_0_xxxyyy_xzzzzz[k] = -g_z_0_xxyyy_xzzzzz[k] * ab_x + g_z_0_xxyyy_xxzzzzz[k];

                g_z_0_xxxyyy_yyyyyy[k] = -g_z_0_xxyyy_yyyyyy[k] * ab_x + g_z_0_xxyyy_xyyyyyy[k];

                g_z_0_xxxyyy_yyyyyz[k] = -g_z_0_xxyyy_yyyyyz[k] * ab_x + g_z_0_xxyyy_xyyyyyz[k];

                g_z_0_xxxyyy_yyyyzz[k] = -g_z_0_xxyyy_yyyyzz[k] * ab_x + g_z_0_xxyyy_xyyyyzz[k];

                g_z_0_xxxyyy_yyyzzz[k] = -g_z_0_xxyyy_yyyzzz[k] * ab_x + g_z_0_xxyyy_xyyyzzz[k];

                g_z_0_xxxyyy_yyzzzz[k] = -g_z_0_xxyyy_yyzzzz[k] * ab_x + g_z_0_xxyyy_xyyzzzz[k];

                g_z_0_xxxyyy_yzzzzz[k] = -g_z_0_xxyyy_yzzzzz[k] * ab_x + g_z_0_xxyyy_xyzzzzz[k];

                g_z_0_xxxyyy_zzzzzz[k] = -g_z_0_xxyyy_zzzzzz[k] * ab_x + g_z_0_xxyyy_xzzzzzz[k];
            }

            /// Set up 1764-1792 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 1764 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 1765 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 1766 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 1767 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 1768 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 1769 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 1770 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 1771 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 1772 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 1773 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 1774 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 1775 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 1776 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 1777 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 1778 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 1779 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 1780 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 1781 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 1782 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 1783 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 1784 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 1785 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 1786 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 1787 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 1788 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 1789 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 1790 * ccomps * dcomps);

            auto g_z_0_xxxyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 1791 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyyz_xxxxxx, g_z_0_xxxyyz_xxxxxy, g_z_0_xxxyyz_xxxxxz, g_z_0_xxxyyz_xxxxyy, g_z_0_xxxyyz_xxxxyz, g_z_0_xxxyyz_xxxxzz, g_z_0_xxxyyz_xxxyyy, g_z_0_xxxyyz_xxxyyz, g_z_0_xxxyyz_xxxyzz, g_z_0_xxxyyz_xxxzzz, g_z_0_xxxyyz_xxyyyy, g_z_0_xxxyyz_xxyyyz, g_z_0_xxxyyz_xxyyzz, g_z_0_xxxyyz_xxyzzz, g_z_0_xxxyyz_xxzzzz, g_z_0_xxxyyz_xyyyyy, g_z_0_xxxyyz_xyyyyz, g_z_0_xxxyyz_xyyyzz, g_z_0_xxxyyz_xyyzzz, g_z_0_xxxyyz_xyzzzz, g_z_0_xxxyyz_xzzzzz, g_z_0_xxxyyz_yyyyyy, g_z_0_xxxyyz_yyyyyz, g_z_0_xxxyyz_yyyyzz, g_z_0_xxxyyz_yyyzzz, g_z_0_xxxyyz_yyzzzz, g_z_0_xxxyyz_yzzzzz, g_z_0_xxxyyz_zzzzzz, g_z_0_xxyyz_xxxxxx, g_z_0_xxyyz_xxxxxxx, g_z_0_xxyyz_xxxxxxy, g_z_0_xxyyz_xxxxxxz, g_z_0_xxyyz_xxxxxy, g_z_0_xxyyz_xxxxxyy, g_z_0_xxyyz_xxxxxyz, g_z_0_xxyyz_xxxxxz, g_z_0_xxyyz_xxxxxzz, g_z_0_xxyyz_xxxxyy, g_z_0_xxyyz_xxxxyyy, g_z_0_xxyyz_xxxxyyz, g_z_0_xxyyz_xxxxyz, g_z_0_xxyyz_xxxxyzz, g_z_0_xxyyz_xxxxzz, g_z_0_xxyyz_xxxxzzz, g_z_0_xxyyz_xxxyyy, g_z_0_xxyyz_xxxyyyy, g_z_0_xxyyz_xxxyyyz, g_z_0_xxyyz_xxxyyz, g_z_0_xxyyz_xxxyyzz, g_z_0_xxyyz_xxxyzz, g_z_0_xxyyz_xxxyzzz, g_z_0_xxyyz_xxxzzz, g_z_0_xxyyz_xxxzzzz, g_z_0_xxyyz_xxyyyy, g_z_0_xxyyz_xxyyyyy, g_z_0_xxyyz_xxyyyyz, g_z_0_xxyyz_xxyyyz, g_z_0_xxyyz_xxyyyzz, g_z_0_xxyyz_xxyyzz, g_z_0_xxyyz_xxyyzzz, g_z_0_xxyyz_xxyzzz, g_z_0_xxyyz_xxyzzzz, g_z_0_xxyyz_xxzzzz, g_z_0_xxyyz_xxzzzzz, g_z_0_xxyyz_xyyyyy, g_z_0_xxyyz_xyyyyyy, g_z_0_xxyyz_xyyyyyz, g_z_0_xxyyz_xyyyyz, g_z_0_xxyyz_xyyyyzz, g_z_0_xxyyz_xyyyzz, g_z_0_xxyyz_xyyyzzz, g_z_0_xxyyz_xyyzzz, g_z_0_xxyyz_xyyzzzz, g_z_0_xxyyz_xyzzzz, g_z_0_xxyyz_xyzzzzz, g_z_0_xxyyz_xzzzzz, g_z_0_xxyyz_xzzzzzz, g_z_0_xxyyz_yyyyyy, g_z_0_xxyyz_yyyyyz, g_z_0_xxyyz_yyyyzz, g_z_0_xxyyz_yyyzzz, g_z_0_xxyyz_yyzzzz, g_z_0_xxyyz_yzzzzz, g_z_0_xxyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_xxxxxx[k] = -g_z_0_xxyyz_xxxxxx[k] * ab_x + g_z_0_xxyyz_xxxxxxx[k];

                g_z_0_xxxyyz_xxxxxy[k] = -g_z_0_xxyyz_xxxxxy[k] * ab_x + g_z_0_xxyyz_xxxxxxy[k];

                g_z_0_xxxyyz_xxxxxz[k] = -g_z_0_xxyyz_xxxxxz[k] * ab_x + g_z_0_xxyyz_xxxxxxz[k];

                g_z_0_xxxyyz_xxxxyy[k] = -g_z_0_xxyyz_xxxxyy[k] * ab_x + g_z_0_xxyyz_xxxxxyy[k];

                g_z_0_xxxyyz_xxxxyz[k] = -g_z_0_xxyyz_xxxxyz[k] * ab_x + g_z_0_xxyyz_xxxxxyz[k];

                g_z_0_xxxyyz_xxxxzz[k] = -g_z_0_xxyyz_xxxxzz[k] * ab_x + g_z_0_xxyyz_xxxxxzz[k];

                g_z_0_xxxyyz_xxxyyy[k] = -g_z_0_xxyyz_xxxyyy[k] * ab_x + g_z_0_xxyyz_xxxxyyy[k];

                g_z_0_xxxyyz_xxxyyz[k] = -g_z_0_xxyyz_xxxyyz[k] * ab_x + g_z_0_xxyyz_xxxxyyz[k];

                g_z_0_xxxyyz_xxxyzz[k] = -g_z_0_xxyyz_xxxyzz[k] * ab_x + g_z_0_xxyyz_xxxxyzz[k];

                g_z_0_xxxyyz_xxxzzz[k] = -g_z_0_xxyyz_xxxzzz[k] * ab_x + g_z_0_xxyyz_xxxxzzz[k];

                g_z_0_xxxyyz_xxyyyy[k] = -g_z_0_xxyyz_xxyyyy[k] * ab_x + g_z_0_xxyyz_xxxyyyy[k];

                g_z_0_xxxyyz_xxyyyz[k] = -g_z_0_xxyyz_xxyyyz[k] * ab_x + g_z_0_xxyyz_xxxyyyz[k];

                g_z_0_xxxyyz_xxyyzz[k] = -g_z_0_xxyyz_xxyyzz[k] * ab_x + g_z_0_xxyyz_xxxyyzz[k];

                g_z_0_xxxyyz_xxyzzz[k] = -g_z_0_xxyyz_xxyzzz[k] * ab_x + g_z_0_xxyyz_xxxyzzz[k];

                g_z_0_xxxyyz_xxzzzz[k] = -g_z_0_xxyyz_xxzzzz[k] * ab_x + g_z_0_xxyyz_xxxzzzz[k];

                g_z_0_xxxyyz_xyyyyy[k] = -g_z_0_xxyyz_xyyyyy[k] * ab_x + g_z_0_xxyyz_xxyyyyy[k];

                g_z_0_xxxyyz_xyyyyz[k] = -g_z_0_xxyyz_xyyyyz[k] * ab_x + g_z_0_xxyyz_xxyyyyz[k];

                g_z_0_xxxyyz_xyyyzz[k] = -g_z_0_xxyyz_xyyyzz[k] * ab_x + g_z_0_xxyyz_xxyyyzz[k];

                g_z_0_xxxyyz_xyyzzz[k] = -g_z_0_xxyyz_xyyzzz[k] * ab_x + g_z_0_xxyyz_xxyyzzz[k];

                g_z_0_xxxyyz_xyzzzz[k] = -g_z_0_xxyyz_xyzzzz[k] * ab_x + g_z_0_xxyyz_xxyzzzz[k];

                g_z_0_xxxyyz_xzzzzz[k] = -g_z_0_xxyyz_xzzzzz[k] * ab_x + g_z_0_xxyyz_xxzzzzz[k];

                g_z_0_xxxyyz_yyyyyy[k] = -g_z_0_xxyyz_yyyyyy[k] * ab_x + g_z_0_xxyyz_xyyyyyy[k];

                g_z_0_xxxyyz_yyyyyz[k] = -g_z_0_xxyyz_yyyyyz[k] * ab_x + g_z_0_xxyyz_xyyyyyz[k];

                g_z_0_xxxyyz_yyyyzz[k] = -g_z_0_xxyyz_yyyyzz[k] * ab_x + g_z_0_xxyyz_xyyyyzz[k];

                g_z_0_xxxyyz_yyyzzz[k] = -g_z_0_xxyyz_yyyzzz[k] * ab_x + g_z_0_xxyyz_xyyyzzz[k];

                g_z_0_xxxyyz_yyzzzz[k] = -g_z_0_xxyyz_yyzzzz[k] * ab_x + g_z_0_xxyyz_xyyzzzz[k];

                g_z_0_xxxyyz_yzzzzz[k] = -g_z_0_xxyyz_yzzzzz[k] * ab_x + g_z_0_xxyyz_xyzzzzz[k];

                g_z_0_xxxyyz_zzzzzz[k] = -g_z_0_xxyyz_zzzzzz[k] * ab_x + g_z_0_xxyyz_xzzzzzz[k];
            }

            /// Set up 1792-1820 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1792 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1793 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1794 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1795 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1796 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1797 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1798 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1799 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1800 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1801 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1802 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1803 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1804 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1805 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1806 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1807 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1808 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1809 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1810 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1811 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1812 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1813 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1814 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1815 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1816 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1817 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1818 * ccomps * dcomps);

            auto g_z_0_xxxyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1819 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyzz_xxxxxx, g_z_0_xxxyzz_xxxxxy, g_z_0_xxxyzz_xxxxxz, g_z_0_xxxyzz_xxxxyy, g_z_0_xxxyzz_xxxxyz, g_z_0_xxxyzz_xxxxzz, g_z_0_xxxyzz_xxxyyy, g_z_0_xxxyzz_xxxyyz, g_z_0_xxxyzz_xxxyzz, g_z_0_xxxyzz_xxxzzz, g_z_0_xxxyzz_xxyyyy, g_z_0_xxxyzz_xxyyyz, g_z_0_xxxyzz_xxyyzz, g_z_0_xxxyzz_xxyzzz, g_z_0_xxxyzz_xxzzzz, g_z_0_xxxyzz_xyyyyy, g_z_0_xxxyzz_xyyyyz, g_z_0_xxxyzz_xyyyzz, g_z_0_xxxyzz_xyyzzz, g_z_0_xxxyzz_xyzzzz, g_z_0_xxxyzz_xzzzzz, g_z_0_xxxyzz_yyyyyy, g_z_0_xxxyzz_yyyyyz, g_z_0_xxxyzz_yyyyzz, g_z_0_xxxyzz_yyyzzz, g_z_0_xxxyzz_yyzzzz, g_z_0_xxxyzz_yzzzzz, g_z_0_xxxyzz_zzzzzz, g_z_0_xxyzz_xxxxxx, g_z_0_xxyzz_xxxxxxx, g_z_0_xxyzz_xxxxxxy, g_z_0_xxyzz_xxxxxxz, g_z_0_xxyzz_xxxxxy, g_z_0_xxyzz_xxxxxyy, g_z_0_xxyzz_xxxxxyz, g_z_0_xxyzz_xxxxxz, g_z_0_xxyzz_xxxxxzz, g_z_0_xxyzz_xxxxyy, g_z_0_xxyzz_xxxxyyy, g_z_0_xxyzz_xxxxyyz, g_z_0_xxyzz_xxxxyz, g_z_0_xxyzz_xxxxyzz, g_z_0_xxyzz_xxxxzz, g_z_0_xxyzz_xxxxzzz, g_z_0_xxyzz_xxxyyy, g_z_0_xxyzz_xxxyyyy, g_z_0_xxyzz_xxxyyyz, g_z_0_xxyzz_xxxyyz, g_z_0_xxyzz_xxxyyzz, g_z_0_xxyzz_xxxyzz, g_z_0_xxyzz_xxxyzzz, g_z_0_xxyzz_xxxzzz, g_z_0_xxyzz_xxxzzzz, g_z_0_xxyzz_xxyyyy, g_z_0_xxyzz_xxyyyyy, g_z_0_xxyzz_xxyyyyz, g_z_0_xxyzz_xxyyyz, g_z_0_xxyzz_xxyyyzz, g_z_0_xxyzz_xxyyzz, g_z_0_xxyzz_xxyyzzz, g_z_0_xxyzz_xxyzzz, g_z_0_xxyzz_xxyzzzz, g_z_0_xxyzz_xxzzzz, g_z_0_xxyzz_xxzzzzz, g_z_0_xxyzz_xyyyyy, g_z_0_xxyzz_xyyyyyy, g_z_0_xxyzz_xyyyyyz, g_z_0_xxyzz_xyyyyz, g_z_0_xxyzz_xyyyyzz, g_z_0_xxyzz_xyyyzz, g_z_0_xxyzz_xyyyzzz, g_z_0_xxyzz_xyyzzz, g_z_0_xxyzz_xyyzzzz, g_z_0_xxyzz_xyzzzz, g_z_0_xxyzz_xyzzzzz, g_z_0_xxyzz_xzzzzz, g_z_0_xxyzz_xzzzzzz, g_z_0_xxyzz_yyyyyy, g_z_0_xxyzz_yyyyyz, g_z_0_xxyzz_yyyyzz, g_z_0_xxyzz_yyyzzz, g_z_0_xxyzz_yyzzzz, g_z_0_xxyzz_yzzzzz, g_z_0_xxyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_xxxxxx[k] = -g_z_0_xxyzz_xxxxxx[k] * ab_x + g_z_0_xxyzz_xxxxxxx[k];

                g_z_0_xxxyzz_xxxxxy[k] = -g_z_0_xxyzz_xxxxxy[k] * ab_x + g_z_0_xxyzz_xxxxxxy[k];

                g_z_0_xxxyzz_xxxxxz[k] = -g_z_0_xxyzz_xxxxxz[k] * ab_x + g_z_0_xxyzz_xxxxxxz[k];

                g_z_0_xxxyzz_xxxxyy[k] = -g_z_0_xxyzz_xxxxyy[k] * ab_x + g_z_0_xxyzz_xxxxxyy[k];

                g_z_0_xxxyzz_xxxxyz[k] = -g_z_0_xxyzz_xxxxyz[k] * ab_x + g_z_0_xxyzz_xxxxxyz[k];

                g_z_0_xxxyzz_xxxxzz[k] = -g_z_0_xxyzz_xxxxzz[k] * ab_x + g_z_0_xxyzz_xxxxxzz[k];

                g_z_0_xxxyzz_xxxyyy[k] = -g_z_0_xxyzz_xxxyyy[k] * ab_x + g_z_0_xxyzz_xxxxyyy[k];

                g_z_0_xxxyzz_xxxyyz[k] = -g_z_0_xxyzz_xxxyyz[k] * ab_x + g_z_0_xxyzz_xxxxyyz[k];

                g_z_0_xxxyzz_xxxyzz[k] = -g_z_0_xxyzz_xxxyzz[k] * ab_x + g_z_0_xxyzz_xxxxyzz[k];

                g_z_0_xxxyzz_xxxzzz[k] = -g_z_0_xxyzz_xxxzzz[k] * ab_x + g_z_0_xxyzz_xxxxzzz[k];

                g_z_0_xxxyzz_xxyyyy[k] = -g_z_0_xxyzz_xxyyyy[k] * ab_x + g_z_0_xxyzz_xxxyyyy[k];

                g_z_0_xxxyzz_xxyyyz[k] = -g_z_0_xxyzz_xxyyyz[k] * ab_x + g_z_0_xxyzz_xxxyyyz[k];

                g_z_0_xxxyzz_xxyyzz[k] = -g_z_0_xxyzz_xxyyzz[k] * ab_x + g_z_0_xxyzz_xxxyyzz[k];

                g_z_0_xxxyzz_xxyzzz[k] = -g_z_0_xxyzz_xxyzzz[k] * ab_x + g_z_0_xxyzz_xxxyzzz[k];

                g_z_0_xxxyzz_xxzzzz[k] = -g_z_0_xxyzz_xxzzzz[k] * ab_x + g_z_0_xxyzz_xxxzzzz[k];

                g_z_0_xxxyzz_xyyyyy[k] = -g_z_0_xxyzz_xyyyyy[k] * ab_x + g_z_0_xxyzz_xxyyyyy[k];

                g_z_0_xxxyzz_xyyyyz[k] = -g_z_0_xxyzz_xyyyyz[k] * ab_x + g_z_0_xxyzz_xxyyyyz[k];

                g_z_0_xxxyzz_xyyyzz[k] = -g_z_0_xxyzz_xyyyzz[k] * ab_x + g_z_0_xxyzz_xxyyyzz[k];

                g_z_0_xxxyzz_xyyzzz[k] = -g_z_0_xxyzz_xyyzzz[k] * ab_x + g_z_0_xxyzz_xxyyzzz[k];

                g_z_0_xxxyzz_xyzzzz[k] = -g_z_0_xxyzz_xyzzzz[k] * ab_x + g_z_0_xxyzz_xxyzzzz[k];

                g_z_0_xxxyzz_xzzzzz[k] = -g_z_0_xxyzz_xzzzzz[k] * ab_x + g_z_0_xxyzz_xxzzzzz[k];

                g_z_0_xxxyzz_yyyyyy[k] = -g_z_0_xxyzz_yyyyyy[k] * ab_x + g_z_0_xxyzz_xyyyyyy[k];

                g_z_0_xxxyzz_yyyyyz[k] = -g_z_0_xxyzz_yyyyyz[k] * ab_x + g_z_0_xxyzz_xyyyyyz[k];

                g_z_0_xxxyzz_yyyyzz[k] = -g_z_0_xxyzz_yyyyzz[k] * ab_x + g_z_0_xxyzz_xyyyyzz[k];

                g_z_0_xxxyzz_yyyzzz[k] = -g_z_0_xxyzz_yyyzzz[k] * ab_x + g_z_0_xxyzz_xyyyzzz[k];

                g_z_0_xxxyzz_yyzzzz[k] = -g_z_0_xxyzz_yyzzzz[k] * ab_x + g_z_0_xxyzz_xyyzzzz[k];

                g_z_0_xxxyzz_yzzzzz[k] = -g_z_0_xxyzz_yzzzzz[k] * ab_x + g_z_0_xxyzz_xyzzzzz[k];

                g_z_0_xxxyzz_zzzzzz[k] = -g_z_0_xxyzz_zzzzzz[k] * ab_x + g_z_0_xxyzz_xzzzzzz[k];
            }

            /// Set up 1820-1848 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1820 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1821 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1822 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1823 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1824 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1825 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1826 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1827 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1828 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1829 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1830 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1831 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1832 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1833 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1834 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1835 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1836 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1837 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1838 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1839 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1840 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1841 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1842 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1843 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1844 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1845 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1846 * ccomps * dcomps);

            auto g_z_0_xxxzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1847 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxzzz_xxxxxx, g_z_0_xxxzzz_xxxxxy, g_z_0_xxxzzz_xxxxxz, g_z_0_xxxzzz_xxxxyy, g_z_0_xxxzzz_xxxxyz, g_z_0_xxxzzz_xxxxzz, g_z_0_xxxzzz_xxxyyy, g_z_0_xxxzzz_xxxyyz, g_z_0_xxxzzz_xxxyzz, g_z_0_xxxzzz_xxxzzz, g_z_0_xxxzzz_xxyyyy, g_z_0_xxxzzz_xxyyyz, g_z_0_xxxzzz_xxyyzz, g_z_0_xxxzzz_xxyzzz, g_z_0_xxxzzz_xxzzzz, g_z_0_xxxzzz_xyyyyy, g_z_0_xxxzzz_xyyyyz, g_z_0_xxxzzz_xyyyzz, g_z_0_xxxzzz_xyyzzz, g_z_0_xxxzzz_xyzzzz, g_z_0_xxxzzz_xzzzzz, g_z_0_xxxzzz_yyyyyy, g_z_0_xxxzzz_yyyyyz, g_z_0_xxxzzz_yyyyzz, g_z_0_xxxzzz_yyyzzz, g_z_0_xxxzzz_yyzzzz, g_z_0_xxxzzz_yzzzzz, g_z_0_xxxzzz_zzzzzz, g_z_0_xxzzz_xxxxxx, g_z_0_xxzzz_xxxxxxx, g_z_0_xxzzz_xxxxxxy, g_z_0_xxzzz_xxxxxxz, g_z_0_xxzzz_xxxxxy, g_z_0_xxzzz_xxxxxyy, g_z_0_xxzzz_xxxxxyz, g_z_0_xxzzz_xxxxxz, g_z_0_xxzzz_xxxxxzz, g_z_0_xxzzz_xxxxyy, g_z_0_xxzzz_xxxxyyy, g_z_0_xxzzz_xxxxyyz, g_z_0_xxzzz_xxxxyz, g_z_0_xxzzz_xxxxyzz, g_z_0_xxzzz_xxxxzz, g_z_0_xxzzz_xxxxzzz, g_z_0_xxzzz_xxxyyy, g_z_0_xxzzz_xxxyyyy, g_z_0_xxzzz_xxxyyyz, g_z_0_xxzzz_xxxyyz, g_z_0_xxzzz_xxxyyzz, g_z_0_xxzzz_xxxyzz, g_z_0_xxzzz_xxxyzzz, g_z_0_xxzzz_xxxzzz, g_z_0_xxzzz_xxxzzzz, g_z_0_xxzzz_xxyyyy, g_z_0_xxzzz_xxyyyyy, g_z_0_xxzzz_xxyyyyz, g_z_0_xxzzz_xxyyyz, g_z_0_xxzzz_xxyyyzz, g_z_0_xxzzz_xxyyzz, g_z_0_xxzzz_xxyyzzz, g_z_0_xxzzz_xxyzzz, g_z_0_xxzzz_xxyzzzz, g_z_0_xxzzz_xxzzzz, g_z_0_xxzzz_xxzzzzz, g_z_0_xxzzz_xyyyyy, g_z_0_xxzzz_xyyyyyy, g_z_0_xxzzz_xyyyyyz, g_z_0_xxzzz_xyyyyz, g_z_0_xxzzz_xyyyyzz, g_z_0_xxzzz_xyyyzz, g_z_0_xxzzz_xyyyzzz, g_z_0_xxzzz_xyyzzz, g_z_0_xxzzz_xyyzzzz, g_z_0_xxzzz_xyzzzz, g_z_0_xxzzz_xyzzzzz, g_z_0_xxzzz_xzzzzz, g_z_0_xxzzz_xzzzzzz, g_z_0_xxzzz_yyyyyy, g_z_0_xxzzz_yyyyyz, g_z_0_xxzzz_yyyyzz, g_z_0_xxzzz_yyyzzz, g_z_0_xxzzz_yyzzzz, g_z_0_xxzzz_yzzzzz, g_z_0_xxzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_xxxxxx[k] = -g_z_0_xxzzz_xxxxxx[k] * ab_x + g_z_0_xxzzz_xxxxxxx[k];

                g_z_0_xxxzzz_xxxxxy[k] = -g_z_0_xxzzz_xxxxxy[k] * ab_x + g_z_0_xxzzz_xxxxxxy[k];

                g_z_0_xxxzzz_xxxxxz[k] = -g_z_0_xxzzz_xxxxxz[k] * ab_x + g_z_0_xxzzz_xxxxxxz[k];

                g_z_0_xxxzzz_xxxxyy[k] = -g_z_0_xxzzz_xxxxyy[k] * ab_x + g_z_0_xxzzz_xxxxxyy[k];

                g_z_0_xxxzzz_xxxxyz[k] = -g_z_0_xxzzz_xxxxyz[k] * ab_x + g_z_0_xxzzz_xxxxxyz[k];

                g_z_0_xxxzzz_xxxxzz[k] = -g_z_0_xxzzz_xxxxzz[k] * ab_x + g_z_0_xxzzz_xxxxxzz[k];

                g_z_0_xxxzzz_xxxyyy[k] = -g_z_0_xxzzz_xxxyyy[k] * ab_x + g_z_0_xxzzz_xxxxyyy[k];

                g_z_0_xxxzzz_xxxyyz[k] = -g_z_0_xxzzz_xxxyyz[k] * ab_x + g_z_0_xxzzz_xxxxyyz[k];

                g_z_0_xxxzzz_xxxyzz[k] = -g_z_0_xxzzz_xxxyzz[k] * ab_x + g_z_0_xxzzz_xxxxyzz[k];

                g_z_0_xxxzzz_xxxzzz[k] = -g_z_0_xxzzz_xxxzzz[k] * ab_x + g_z_0_xxzzz_xxxxzzz[k];

                g_z_0_xxxzzz_xxyyyy[k] = -g_z_0_xxzzz_xxyyyy[k] * ab_x + g_z_0_xxzzz_xxxyyyy[k];

                g_z_0_xxxzzz_xxyyyz[k] = -g_z_0_xxzzz_xxyyyz[k] * ab_x + g_z_0_xxzzz_xxxyyyz[k];

                g_z_0_xxxzzz_xxyyzz[k] = -g_z_0_xxzzz_xxyyzz[k] * ab_x + g_z_0_xxzzz_xxxyyzz[k];

                g_z_0_xxxzzz_xxyzzz[k] = -g_z_0_xxzzz_xxyzzz[k] * ab_x + g_z_0_xxzzz_xxxyzzz[k];

                g_z_0_xxxzzz_xxzzzz[k] = -g_z_0_xxzzz_xxzzzz[k] * ab_x + g_z_0_xxzzz_xxxzzzz[k];

                g_z_0_xxxzzz_xyyyyy[k] = -g_z_0_xxzzz_xyyyyy[k] * ab_x + g_z_0_xxzzz_xxyyyyy[k];

                g_z_0_xxxzzz_xyyyyz[k] = -g_z_0_xxzzz_xyyyyz[k] * ab_x + g_z_0_xxzzz_xxyyyyz[k];

                g_z_0_xxxzzz_xyyyzz[k] = -g_z_0_xxzzz_xyyyzz[k] * ab_x + g_z_0_xxzzz_xxyyyzz[k];

                g_z_0_xxxzzz_xyyzzz[k] = -g_z_0_xxzzz_xyyzzz[k] * ab_x + g_z_0_xxzzz_xxyyzzz[k];

                g_z_0_xxxzzz_xyzzzz[k] = -g_z_0_xxzzz_xyzzzz[k] * ab_x + g_z_0_xxzzz_xxyzzzz[k];

                g_z_0_xxxzzz_xzzzzz[k] = -g_z_0_xxzzz_xzzzzz[k] * ab_x + g_z_0_xxzzz_xxzzzzz[k];

                g_z_0_xxxzzz_yyyyyy[k] = -g_z_0_xxzzz_yyyyyy[k] * ab_x + g_z_0_xxzzz_xyyyyyy[k];

                g_z_0_xxxzzz_yyyyyz[k] = -g_z_0_xxzzz_yyyyyz[k] * ab_x + g_z_0_xxzzz_xyyyyyz[k];

                g_z_0_xxxzzz_yyyyzz[k] = -g_z_0_xxzzz_yyyyzz[k] * ab_x + g_z_0_xxzzz_xyyyyzz[k];

                g_z_0_xxxzzz_yyyzzz[k] = -g_z_0_xxzzz_yyyzzz[k] * ab_x + g_z_0_xxzzz_xyyyzzz[k];

                g_z_0_xxxzzz_yyzzzz[k] = -g_z_0_xxzzz_yyzzzz[k] * ab_x + g_z_0_xxzzz_xyyzzzz[k];

                g_z_0_xxxzzz_yzzzzz[k] = -g_z_0_xxzzz_yzzzzz[k] * ab_x + g_z_0_xxzzz_xyzzzzz[k];

                g_z_0_xxxzzz_zzzzzz[k] = -g_z_0_xxzzz_zzzzzz[k] * ab_x + g_z_0_xxzzz_xzzzzzz[k];
            }

            /// Set up 1848-1876 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 1848 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 1849 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 1850 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 1851 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 1852 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 1853 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 1854 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 1855 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 1856 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 1857 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 1858 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 1859 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 1860 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 1861 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 1862 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 1863 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 1864 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 1865 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 1866 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 1867 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 1868 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 1869 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 1870 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 1871 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 1872 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 1873 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 1874 * ccomps * dcomps);

            auto g_z_0_xxyyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 1875 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyyy_xxxxxx, g_z_0_xxyyyy_xxxxxy, g_z_0_xxyyyy_xxxxxz, g_z_0_xxyyyy_xxxxyy, g_z_0_xxyyyy_xxxxyz, g_z_0_xxyyyy_xxxxzz, g_z_0_xxyyyy_xxxyyy, g_z_0_xxyyyy_xxxyyz, g_z_0_xxyyyy_xxxyzz, g_z_0_xxyyyy_xxxzzz, g_z_0_xxyyyy_xxyyyy, g_z_0_xxyyyy_xxyyyz, g_z_0_xxyyyy_xxyyzz, g_z_0_xxyyyy_xxyzzz, g_z_0_xxyyyy_xxzzzz, g_z_0_xxyyyy_xyyyyy, g_z_0_xxyyyy_xyyyyz, g_z_0_xxyyyy_xyyyzz, g_z_0_xxyyyy_xyyzzz, g_z_0_xxyyyy_xyzzzz, g_z_0_xxyyyy_xzzzzz, g_z_0_xxyyyy_yyyyyy, g_z_0_xxyyyy_yyyyyz, g_z_0_xxyyyy_yyyyzz, g_z_0_xxyyyy_yyyzzz, g_z_0_xxyyyy_yyzzzz, g_z_0_xxyyyy_yzzzzz, g_z_0_xxyyyy_zzzzzz, g_z_0_xyyyy_xxxxxx, g_z_0_xyyyy_xxxxxxx, g_z_0_xyyyy_xxxxxxy, g_z_0_xyyyy_xxxxxxz, g_z_0_xyyyy_xxxxxy, g_z_0_xyyyy_xxxxxyy, g_z_0_xyyyy_xxxxxyz, g_z_0_xyyyy_xxxxxz, g_z_0_xyyyy_xxxxxzz, g_z_0_xyyyy_xxxxyy, g_z_0_xyyyy_xxxxyyy, g_z_0_xyyyy_xxxxyyz, g_z_0_xyyyy_xxxxyz, g_z_0_xyyyy_xxxxyzz, g_z_0_xyyyy_xxxxzz, g_z_0_xyyyy_xxxxzzz, g_z_0_xyyyy_xxxyyy, g_z_0_xyyyy_xxxyyyy, g_z_0_xyyyy_xxxyyyz, g_z_0_xyyyy_xxxyyz, g_z_0_xyyyy_xxxyyzz, g_z_0_xyyyy_xxxyzz, g_z_0_xyyyy_xxxyzzz, g_z_0_xyyyy_xxxzzz, g_z_0_xyyyy_xxxzzzz, g_z_0_xyyyy_xxyyyy, g_z_0_xyyyy_xxyyyyy, g_z_0_xyyyy_xxyyyyz, g_z_0_xyyyy_xxyyyz, g_z_0_xyyyy_xxyyyzz, g_z_0_xyyyy_xxyyzz, g_z_0_xyyyy_xxyyzzz, g_z_0_xyyyy_xxyzzz, g_z_0_xyyyy_xxyzzzz, g_z_0_xyyyy_xxzzzz, g_z_0_xyyyy_xxzzzzz, g_z_0_xyyyy_xyyyyy, g_z_0_xyyyy_xyyyyyy, g_z_0_xyyyy_xyyyyyz, g_z_0_xyyyy_xyyyyz, g_z_0_xyyyy_xyyyyzz, g_z_0_xyyyy_xyyyzz, g_z_0_xyyyy_xyyyzzz, g_z_0_xyyyy_xyyzzz, g_z_0_xyyyy_xyyzzzz, g_z_0_xyyyy_xyzzzz, g_z_0_xyyyy_xyzzzzz, g_z_0_xyyyy_xzzzzz, g_z_0_xyyyy_xzzzzzz, g_z_0_xyyyy_yyyyyy, g_z_0_xyyyy_yyyyyz, g_z_0_xyyyy_yyyyzz, g_z_0_xyyyy_yyyzzz, g_z_0_xyyyy_yyzzzz, g_z_0_xyyyy_yzzzzz, g_z_0_xyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_xxxxxx[k] = -g_z_0_xyyyy_xxxxxx[k] * ab_x + g_z_0_xyyyy_xxxxxxx[k];

                g_z_0_xxyyyy_xxxxxy[k] = -g_z_0_xyyyy_xxxxxy[k] * ab_x + g_z_0_xyyyy_xxxxxxy[k];

                g_z_0_xxyyyy_xxxxxz[k] = -g_z_0_xyyyy_xxxxxz[k] * ab_x + g_z_0_xyyyy_xxxxxxz[k];

                g_z_0_xxyyyy_xxxxyy[k] = -g_z_0_xyyyy_xxxxyy[k] * ab_x + g_z_0_xyyyy_xxxxxyy[k];

                g_z_0_xxyyyy_xxxxyz[k] = -g_z_0_xyyyy_xxxxyz[k] * ab_x + g_z_0_xyyyy_xxxxxyz[k];

                g_z_0_xxyyyy_xxxxzz[k] = -g_z_0_xyyyy_xxxxzz[k] * ab_x + g_z_0_xyyyy_xxxxxzz[k];

                g_z_0_xxyyyy_xxxyyy[k] = -g_z_0_xyyyy_xxxyyy[k] * ab_x + g_z_0_xyyyy_xxxxyyy[k];

                g_z_0_xxyyyy_xxxyyz[k] = -g_z_0_xyyyy_xxxyyz[k] * ab_x + g_z_0_xyyyy_xxxxyyz[k];

                g_z_0_xxyyyy_xxxyzz[k] = -g_z_0_xyyyy_xxxyzz[k] * ab_x + g_z_0_xyyyy_xxxxyzz[k];

                g_z_0_xxyyyy_xxxzzz[k] = -g_z_0_xyyyy_xxxzzz[k] * ab_x + g_z_0_xyyyy_xxxxzzz[k];

                g_z_0_xxyyyy_xxyyyy[k] = -g_z_0_xyyyy_xxyyyy[k] * ab_x + g_z_0_xyyyy_xxxyyyy[k];

                g_z_0_xxyyyy_xxyyyz[k] = -g_z_0_xyyyy_xxyyyz[k] * ab_x + g_z_0_xyyyy_xxxyyyz[k];

                g_z_0_xxyyyy_xxyyzz[k] = -g_z_0_xyyyy_xxyyzz[k] * ab_x + g_z_0_xyyyy_xxxyyzz[k];

                g_z_0_xxyyyy_xxyzzz[k] = -g_z_0_xyyyy_xxyzzz[k] * ab_x + g_z_0_xyyyy_xxxyzzz[k];

                g_z_0_xxyyyy_xxzzzz[k] = -g_z_0_xyyyy_xxzzzz[k] * ab_x + g_z_0_xyyyy_xxxzzzz[k];

                g_z_0_xxyyyy_xyyyyy[k] = -g_z_0_xyyyy_xyyyyy[k] * ab_x + g_z_0_xyyyy_xxyyyyy[k];

                g_z_0_xxyyyy_xyyyyz[k] = -g_z_0_xyyyy_xyyyyz[k] * ab_x + g_z_0_xyyyy_xxyyyyz[k];

                g_z_0_xxyyyy_xyyyzz[k] = -g_z_0_xyyyy_xyyyzz[k] * ab_x + g_z_0_xyyyy_xxyyyzz[k];

                g_z_0_xxyyyy_xyyzzz[k] = -g_z_0_xyyyy_xyyzzz[k] * ab_x + g_z_0_xyyyy_xxyyzzz[k];

                g_z_0_xxyyyy_xyzzzz[k] = -g_z_0_xyyyy_xyzzzz[k] * ab_x + g_z_0_xyyyy_xxyzzzz[k];

                g_z_0_xxyyyy_xzzzzz[k] = -g_z_0_xyyyy_xzzzzz[k] * ab_x + g_z_0_xyyyy_xxzzzzz[k];

                g_z_0_xxyyyy_yyyyyy[k] = -g_z_0_xyyyy_yyyyyy[k] * ab_x + g_z_0_xyyyy_xyyyyyy[k];

                g_z_0_xxyyyy_yyyyyz[k] = -g_z_0_xyyyy_yyyyyz[k] * ab_x + g_z_0_xyyyy_xyyyyyz[k];

                g_z_0_xxyyyy_yyyyzz[k] = -g_z_0_xyyyy_yyyyzz[k] * ab_x + g_z_0_xyyyy_xyyyyzz[k];

                g_z_0_xxyyyy_yyyzzz[k] = -g_z_0_xyyyy_yyyzzz[k] * ab_x + g_z_0_xyyyy_xyyyzzz[k];

                g_z_0_xxyyyy_yyzzzz[k] = -g_z_0_xyyyy_yyzzzz[k] * ab_x + g_z_0_xyyyy_xyyzzzz[k];

                g_z_0_xxyyyy_yzzzzz[k] = -g_z_0_xyyyy_yzzzzz[k] * ab_x + g_z_0_xyyyy_xyzzzzz[k];

                g_z_0_xxyyyy_zzzzzz[k] = -g_z_0_xyyyy_zzzzzz[k] * ab_x + g_z_0_xyyyy_xzzzzzz[k];
            }

            /// Set up 1876-1904 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 1876 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 1877 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 1878 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 1879 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 1880 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 1881 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 1882 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 1883 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 1884 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 1885 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 1886 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 1887 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 1888 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 1889 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 1890 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 1891 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 1892 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 1893 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 1894 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 1895 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 1896 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 1897 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 1898 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 1899 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 1900 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 1901 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 1902 * ccomps * dcomps);

            auto g_z_0_xxyyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 1903 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyyz_xxxxxx, g_z_0_xxyyyz_xxxxxy, g_z_0_xxyyyz_xxxxxz, g_z_0_xxyyyz_xxxxyy, g_z_0_xxyyyz_xxxxyz, g_z_0_xxyyyz_xxxxzz, g_z_0_xxyyyz_xxxyyy, g_z_0_xxyyyz_xxxyyz, g_z_0_xxyyyz_xxxyzz, g_z_0_xxyyyz_xxxzzz, g_z_0_xxyyyz_xxyyyy, g_z_0_xxyyyz_xxyyyz, g_z_0_xxyyyz_xxyyzz, g_z_0_xxyyyz_xxyzzz, g_z_0_xxyyyz_xxzzzz, g_z_0_xxyyyz_xyyyyy, g_z_0_xxyyyz_xyyyyz, g_z_0_xxyyyz_xyyyzz, g_z_0_xxyyyz_xyyzzz, g_z_0_xxyyyz_xyzzzz, g_z_0_xxyyyz_xzzzzz, g_z_0_xxyyyz_yyyyyy, g_z_0_xxyyyz_yyyyyz, g_z_0_xxyyyz_yyyyzz, g_z_0_xxyyyz_yyyzzz, g_z_0_xxyyyz_yyzzzz, g_z_0_xxyyyz_yzzzzz, g_z_0_xxyyyz_zzzzzz, g_z_0_xyyyz_xxxxxx, g_z_0_xyyyz_xxxxxxx, g_z_0_xyyyz_xxxxxxy, g_z_0_xyyyz_xxxxxxz, g_z_0_xyyyz_xxxxxy, g_z_0_xyyyz_xxxxxyy, g_z_0_xyyyz_xxxxxyz, g_z_0_xyyyz_xxxxxz, g_z_0_xyyyz_xxxxxzz, g_z_0_xyyyz_xxxxyy, g_z_0_xyyyz_xxxxyyy, g_z_0_xyyyz_xxxxyyz, g_z_0_xyyyz_xxxxyz, g_z_0_xyyyz_xxxxyzz, g_z_0_xyyyz_xxxxzz, g_z_0_xyyyz_xxxxzzz, g_z_0_xyyyz_xxxyyy, g_z_0_xyyyz_xxxyyyy, g_z_0_xyyyz_xxxyyyz, g_z_0_xyyyz_xxxyyz, g_z_0_xyyyz_xxxyyzz, g_z_0_xyyyz_xxxyzz, g_z_0_xyyyz_xxxyzzz, g_z_0_xyyyz_xxxzzz, g_z_0_xyyyz_xxxzzzz, g_z_0_xyyyz_xxyyyy, g_z_0_xyyyz_xxyyyyy, g_z_0_xyyyz_xxyyyyz, g_z_0_xyyyz_xxyyyz, g_z_0_xyyyz_xxyyyzz, g_z_0_xyyyz_xxyyzz, g_z_0_xyyyz_xxyyzzz, g_z_0_xyyyz_xxyzzz, g_z_0_xyyyz_xxyzzzz, g_z_0_xyyyz_xxzzzz, g_z_0_xyyyz_xxzzzzz, g_z_0_xyyyz_xyyyyy, g_z_0_xyyyz_xyyyyyy, g_z_0_xyyyz_xyyyyyz, g_z_0_xyyyz_xyyyyz, g_z_0_xyyyz_xyyyyzz, g_z_0_xyyyz_xyyyzz, g_z_0_xyyyz_xyyyzzz, g_z_0_xyyyz_xyyzzz, g_z_0_xyyyz_xyyzzzz, g_z_0_xyyyz_xyzzzz, g_z_0_xyyyz_xyzzzzz, g_z_0_xyyyz_xzzzzz, g_z_0_xyyyz_xzzzzzz, g_z_0_xyyyz_yyyyyy, g_z_0_xyyyz_yyyyyz, g_z_0_xyyyz_yyyyzz, g_z_0_xyyyz_yyyzzz, g_z_0_xyyyz_yyzzzz, g_z_0_xyyyz_yzzzzz, g_z_0_xyyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_xxxxxx[k] = -g_z_0_xyyyz_xxxxxx[k] * ab_x + g_z_0_xyyyz_xxxxxxx[k];

                g_z_0_xxyyyz_xxxxxy[k] = -g_z_0_xyyyz_xxxxxy[k] * ab_x + g_z_0_xyyyz_xxxxxxy[k];

                g_z_0_xxyyyz_xxxxxz[k] = -g_z_0_xyyyz_xxxxxz[k] * ab_x + g_z_0_xyyyz_xxxxxxz[k];

                g_z_0_xxyyyz_xxxxyy[k] = -g_z_0_xyyyz_xxxxyy[k] * ab_x + g_z_0_xyyyz_xxxxxyy[k];

                g_z_0_xxyyyz_xxxxyz[k] = -g_z_0_xyyyz_xxxxyz[k] * ab_x + g_z_0_xyyyz_xxxxxyz[k];

                g_z_0_xxyyyz_xxxxzz[k] = -g_z_0_xyyyz_xxxxzz[k] * ab_x + g_z_0_xyyyz_xxxxxzz[k];

                g_z_0_xxyyyz_xxxyyy[k] = -g_z_0_xyyyz_xxxyyy[k] * ab_x + g_z_0_xyyyz_xxxxyyy[k];

                g_z_0_xxyyyz_xxxyyz[k] = -g_z_0_xyyyz_xxxyyz[k] * ab_x + g_z_0_xyyyz_xxxxyyz[k];

                g_z_0_xxyyyz_xxxyzz[k] = -g_z_0_xyyyz_xxxyzz[k] * ab_x + g_z_0_xyyyz_xxxxyzz[k];

                g_z_0_xxyyyz_xxxzzz[k] = -g_z_0_xyyyz_xxxzzz[k] * ab_x + g_z_0_xyyyz_xxxxzzz[k];

                g_z_0_xxyyyz_xxyyyy[k] = -g_z_0_xyyyz_xxyyyy[k] * ab_x + g_z_0_xyyyz_xxxyyyy[k];

                g_z_0_xxyyyz_xxyyyz[k] = -g_z_0_xyyyz_xxyyyz[k] * ab_x + g_z_0_xyyyz_xxxyyyz[k];

                g_z_0_xxyyyz_xxyyzz[k] = -g_z_0_xyyyz_xxyyzz[k] * ab_x + g_z_0_xyyyz_xxxyyzz[k];

                g_z_0_xxyyyz_xxyzzz[k] = -g_z_0_xyyyz_xxyzzz[k] * ab_x + g_z_0_xyyyz_xxxyzzz[k];

                g_z_0_xxyyyz_xxzzzz[k] = -g_z_0_xyyyz_xxzzzz[k] * ab_x + g_z_0_xyyyz_xxxzzzz[k];

                g_z_0_xxyyyz_xyyyyy[k] = -g_z_0_xyyyz_xyyyyy[k] * ab_x + g_z_0_xyyyz_xxyyyyy[k];

                g_z_0_xxyyyz_xyyyyz[k] = -g_z_0_xyyyz_xyyyyz[k] * ab_x + g_z_0_xyyyz_xxyyyyz[k];

                g_z_0_xxyyyz_xyyyzz[k] = -g_z_0_xyyyz_xyyyzz[k] * ab_x + g_z_0_xyyyz_xxyyyzz[k];

                g_z_0_xxyyyz_xyyzzz[k] = -g_z_0_xyyyz_xyyzzz[k] * ab_x + g_z_0_xyyyz_xxyyzzz[k];

                g_z_0_xxyyyz_xyzzzz[k] = -g_z_0_xyyyz_xyzzzz[k] * ab_x + g_z_0_xyyyz_xxyzzzz[k];

                g_z_0_xxyyyz_xzzzzz[k] = -g_z_0_xyyyz_xzzzzz[k] * ab_x + g_z_0_xyyyz_xxzzzzz[k];

                g_z_0_xxyyyz_yyyyyy[k] = -g_z_0_xyyyz_yyyyyy[k] * ab_x + g_z_0_xyyyz_xyyyyyy[k];

                g_z_0_xxyyyz_yyyyyz[k] = -g_z_0_xyyyz_yyyyyz[k] * ab_x + g_z_0_xyyyz_xyyyyyz[k];

                g_z_0_xxyyyz_yyyyzz[k] = -g_z_0_xyyyz_yyyyzz[k] * ab_x + g_z_0_xyyyz_xyyyyzz[k];

                g_z_0_xxyyyz_yyyzzz[k] = -g_z_0_xyyyz_yyyzzz[k] * ab_x + g_z_0_xyyyz_xyyyzzz[k];

                g_z_0_xxyyyz_yyzzzz[k] = -g_z_0_xyyyz_yyzzzz[k] * ab_x + g_z_0_xyyyz_xyyzzzz[k];

                g_z_0_xxyyyz_yzzzzz[k] = -g_z_0_xyyyz_yzzzzz[k] * ab_x + g_z_0_xyyyz_xyzzzzz[k];

                g_z_0_xxyyyz_zzzzzz[k] = -g_z_0_xyyyz_zzzzzz[k] * ab_x + g_z_0_xyyyz_xzzzzzz[k];
            }

            /// Set up 1904-1932 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1904 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1905 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1906 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1907 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1908 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1909 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1910 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1911 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1912 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1913 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1914 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1915 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1916 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1917 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1918 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1919 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1920 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1921 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1922 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1923 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1924 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1925 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1926 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1927 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1928 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1929 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1930 * ccomps * dcomps);

            auto g_z_0_xxyyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1931 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyzz_xxxxxx, g_z_0_xxyyzz_xxxxxy, g_z_0_xxyyzz_xxxxxz, g_z_0_xxyyzz_xxxxyy, g_z_0_xxyyzz_xxxxyz, g_z_0_xxyyzz_xxxxzz, g_z_0_xxyyzz_xxxyyy, g_z_0_xxyyzz_xxxyyz, g_z_0_xxyyzz_xxxyzz, g_z_0_xxyyzz_xxxzzz, g_z_0_xxyyzz_xxyyyy, g_z_0_xxyyzz_xxyyyz, g_z_0_xxyyzz_xxyyzz, g_z_0_xxyyzz_xxyzzz, g_z_0_xxyyzz_xxzzzz, g_z_0_xxyyzz_xyyyyy, g_z_0_xxyyzz_xyyyyz, g_z_0_xxyyzz_xyyyzz, g_z_0_xxyyzz_xyyzzz, g_z_0_xxyyzz_xyzzzz, g_z_0_xxyyzz_xzzzzz, g_z_0_xxyyzz_yyyyyy, g_z_0_xxyyzz_yyyyyz, g_z_0_xxyyzz_yyyyzz, g_z_0_xxyyzz_yyyzzz, g_z_0_xxyyzz_yyzzzz, g_z_0_xxyyzz_yzzzzz, g_z_0_xxyyzz_zzzzzz, g_z_0_xyyzz_xxxxxx, g_z_0_xyyzz_xxxxxxx, g_z_0_xyyzz_xxxxxxy, g_z_0_xyyzz_xxxxxxz, g_z_0_xyyzz_xxxxxy, g_z_0_xyyzz_xxxxxyy, g_z_0_xyyzz_xxxxxyz, g_z_0_xyyzz_xxxxxz, g_z_0_xyyzz_xxxxxzz, g_z_0_xyyzz_xxxxyy, g_z_0_xyyzz_xxxxyyy, g_z_0_xyyzz_xxxxyyz, g_z_0_xyyzz_xxxxyz, g_z_0_xyyzz_xxxxyzz, g_z_0_xyyzz_xxxxzz, g_z_0_xyyzz_xxxxzzz, g_z_0_xyyzz_xxxyyy, g_z_0_xyyzz_xxxyyyy, g_z_0_xyyzz_xxxyyyz, g_z_0_xyyzz_xxxyyz, g_z_0_xyyzz_xxxyyzz, g_z_0_xyyzz_xxxyzz, g_z_0_xyyzz_xxxyzzz, g_z_0_xyyzz_xxxzzz, g_z_0_xyyzz_xxxzzzz, g_z_0_xyyzz_xxyyyy, g_z_0_xyyzz_xxyyyyy, g_z_0_xyyzz_xxyyyyz, g_z_0_xyyzz_xxyyyz, g_z_0_xyyzz_xxyyyzz, g_z_0_xyyzz_xxyyzz, g_z_0_xyyzz_xxyyzzz, g_z_0_xyyzz_xxyzzz, g_z_0_xyyzz_xxyzzzz, g_z_0_xyyzz_xxzzzz, g_z_0_xyyzz_xxzzzzz, g_z_0_xyyzz_xyyyyy, g_z_0_xyyzz_xyyyyyy, g_z_0_xyyzz_xyyyyyz, g_z_0_xyyzz_xyyyyz, g_z_0_xyyzz_xyyyyzz, g_z_0_xyyzz_xyyyzz, g_z_0_xyyzz_xyyyzzz, g_z_0_xyyzz_xyyzzz, g_z_0_xyyzz_xyyzzzz, g_z_0_xyyzz_xyzzzz, g_z_0_xyyzz_xyzzzzz, g_z_0_xyyzz_xzzzzz, g_z_0_xyyzz_xzzzzzz, g_z_0_xyyzz_yyyyyy, g_z_0_xyyzz_yyyyyz, g_z_0_xyyzz_yyyyzz, g_z_0_xyyzz_yyyzzz, g_z_0_xyyzz_yyzzzz, g_z_0_xyyzz_yzzzzz, g_z_0_xyyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_xxxxxx[k] = -g_z_0_xyyzz_xxxxxx[k] * ab_x + g_z_0_xyyzz_xxxxxxx[k];

                g_z_0_xxyyzz_xxxxxy[k] = -g_z_0_xyyzz_xxxxxy[k] * ab_x + g_z_0_xyyzz_xxxxxxy[k];

                g_z_0_xxyyzz_xxxxxz[k] = -g_z_0_xyyzz_xxxxxz[k] * ab_x + g_z_0_xyyzz_xxxxxxz[k];

                g_z_0_xxyyzz_xxxxyy[k] = -g_z_0_xyyzz_xxxxyy[k] * ab_x + g_z_0_xyyzz_xxxxxyy[k];

                g_z_0_xxyyzz_xxxxyz[k] = -g_z_0_xyyzz_xxxxyz[k] * ab_x + g_z_0_xyyzz_xxxxxyz[k];

                g_z_0_xxyyzz_xxxxzz[k] = -g_z_0_xyyzz_xxxxzz[k] * ab_x + g_z_0_xyyzz_xxxxxzz[k];

                g_z_0_xxyyzz_xxxyyy[k] = -g_z_0_xyyzz_xxxyyy[k] * ab_x + g_z_0_xyyzz_xxxxyyy[k];

                g_z_0_xxyyzz_xxxyyz[k] = -g_z_0_xyyzz_xxxyyz[k] * ab_x + g_z_0_xyyzz_xxxxyyz[k];

                g_z_0_xxyyzz_xxxyzz[k] = -g_z_0_xyyzz_xxxyzz[k] * ab_x + g_z_0_xyyzz_xxxxyzz[k];

                g_z_0_xxyyzz_xxxzzz[k] = -g_z_0_xyyzz_xxxzzz[k] * ab_x + g_z_0_xyyzz_xxxxzzz[k];

                g_z_0_xxyyzz_xxyyyy[k] = -g_z_0_xyyzz_xxyyyy[k] * ab_x + g_z_0_xyyzz_xxxyyyy[k];

                g_z_0_xxyyzz_xxyyyz[k] = -g_z_0_xyyzz_xxyyyz[k] * ab_x + g_z_0_xyyzz_xxxyyyz[k];

                g_z_0_xxyyzz_xxyyzz[k] = -g_z_0_xyyzz_xxyyzz[k] * ab_x + g_z_0_xyyzz_xxxyyzz[k];

                g_z_0_xxyyzz_xxyzzz[k] = -g_z_0_xyyzz_xxyzzz[k] * ab_x + g_z_0_xyyzz_xxxyzzz[k];

                g_z_0_xxyyzz_xxzzzz[k] = -g_z_0_xyyzz_xxzzzz[k] * ab_x + g_z_0_xyyzz_xxxzzzz[k];

                g_z_0_xxyyzz_xyyyyy[k] = -g_z_0_xyyzz_xyyyyy[k] * ab_x + g_z_0_xyyzz_xxyyyyy[k];

                g_z_0_xxyyzz_xyyyyz[k] = -g_z_0_xyyzz_xyyyyz[k] * ab_x + g_z_0_xyyzz_xxyyyyz[k];

                g_z_0_xxyyzz_xyyyzz[k] = -g_z_0_xyyzz_xyyyzz[k] * ab_x + g_z_0_xyyzz_xxyyyzz[k];

                g_z_0_xxyyzz_xyyzzz[k] = -g_z_0_xyyzz_xyyzzz[k] * ab_x + g_z_0_xyyzz_xxyyzzz[k];

                g_z_0_xxyyzz_xyzzzz[k] = -g_z_0_xyyzz_xyzzzz[k] * ab_x + g_z_0_xyyzz_xxyzzzz[k];

                g_z_0_xxyyzz_xzzzzz[k] = -g_z_0_xyyzz_xzzzzz[k] * ab_x + g_z_0_xyyzz_xxzzzzz[k];

                g_z_0_xxyyzz_yyyyyy[k] = -g_z_0_xyyzz_yyyyyy[k] * ab_x + g_z_0_xyyzz_xyyyyyy[k];

                g_z_0_xxyyzz_yyyyyz[k] = -g_z_0_xyyzz_yyyyyz[k] * ab_x + g_z_0_xyyzz_xyyyyyz[k];

                g_z_0_xxyyzz_yyyyzz[k] = -g_z_0_xyyzz_yyyyzz[k] * ab_x + g_z_0_xyyzz_xyyyyzz[k];

                g_z_0_xxyyzz_yyyzzz[k] = -g_z_0_xyyzz_yyyzzz[k] * ab_x + g_z_0_xyyzz_xyyyzzz[k];

                g_z_0_xxyyzz_yyzzzz[k] = -g_z_0_xyyzz_yyzzzz[k] * ab_x + g_z_0_xyyzz_xyyzzzz[k];

                g_z_0_xxyyzz_yzzzzz[k] = -g_z_0_xyyzz_yzzzzz[k] * ab_x + g_z_0_xyyzz_xyzzzzz[k];

                g_z_0_xxyyzz_zzzzzz[k] = -g_z_0_xyyzz_zzzzzz[k] * ab_x + g_z_0_xyyzz_xzzzzzz[k];
            }

            /// Set up 1932-1960 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1932 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1933 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1934 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1935 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1936 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1937 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1938 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1939 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1940 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1941 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1942 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1943 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1944 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1945 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1946 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1947 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1948 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1949 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1950 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1951 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1952 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1953 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1954 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1955 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1956 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1957 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1958 * ccomps * dcomps);

            auto g_z_0_xxyzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1959 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyzzz_xxxxxx, g_z_0_xxyzzz_xxxxxy, g_z_0_xxyzzz_xxxxxz, g_z_0_xxyzzz_xxxxyy, g_z_0_xxyzzz_xxxxyz, g_z_0_xxyzzz_xxxxzz, g_z_0_xxyzzz_xxxyyy, g_z_0_xxyzzz_xxxyyz, g_z_0_xxyzzz_xxxyzz, g_z_0_xxyzzz_xxxzzz, g_z_0_xxyzzz_xxyyyy, g_z_0_xxyzzz_xxyyyz, g_z_0_xxyzzz_xxyyzz, g_z_0_xxyzzz_xxyzzz, g_z_0_xxyzzz_xxzzzz, g_z_0_xxyzzz_xyyyyy, g_z_0_xxyzzz_xyyyyz, g_z_0_xxyzzz_xyyyzz, g_z_0_xxyzzz_xyyzzz, g_z_0_xxyzzz_xyzzzz, g_z_0_xxyzzz_xzzzzz, g_z_0_xxyzzz_yyyyyy, g_z_0_xxyzzz_yyyyyz, g_z_0_xxyzzz_yyyyzz, g_z_0_xxyzzz_yyyzzz, g_z_0_xxyzzz_yyzzzz, g_z_0_xxyzzz_yzzzzz, g_z_0_xxyzzz_zzzzzz, g_z_0_xyzzz_xxxxxx, g_z_0_xyzzz_xxxxxxx, g_z_0_xyzzz_xxxxxxy, g_z_0_xyzzz_xxxxxxz, g_z_0_xyzzz_xxxxxy, g_z_0_xyzzz_xxxxxyy, g_z_0_xyzzz_xxxxxyz, g_z_0_xyzzz_xxxxxz, g_z_0_xyzzz_xxxxxzz, g_z_0_xyzzz_xxxxyy, g_z_0_xyzzz_xxxxyyy, g_z_0_xyzzz_xxxxyyz, g_z_0_xyzzz_xxxxyz, g_z_0_xyzzz_xxxxyzz, g_z_0_xyzzz_xxxxzz, g_z_0_xyzzz_xxxxzzz, g_z_0_xyzzz_xxxyyy, g_z_0_xyzzz_xxxyyyy, g_z_0_xyzzz_xxxyyyz, g_z_0_xyzzz_xxxyyz, g_z_0_xyzzz_xxxyyzz, g_z_0_xyzzz_xxxyzz, g_z_0_xyzzz_xxxyzzz, g_z_0_xyzzz_xxxzzz, g_z_0_xyzzz_xxxzzzz, g_z_0_xyzzz_xxyyyy, g_z_0_xyzzz_xxyyyyy, g_z_0_xyzzz_xxyyyyz, g_z_0_xyzzz_xxyyyz, g_z_0_xyzzz_xxyyyzz, g_z_0_xyzzz_xxyyzz, g_z_0_xyzzz_xxyyzzz, g_z_0_xyzzz_xxyzzz, g_z_0_xyzzz_xxyzzzz, g_z_0_xyzzz_xxzzzz, g_z_0_xyzzz_xxzzzzz, g_z_0_xyzzz_xyyyyy, g_z_0_xyzzz_xyyyyyy, g_z_0_xyzzz_xyyyyyz, g_z_0_xyzzz_xyyyyz, g_z_0_xyzzz_xyyyyzz, g_z_0_xyzzz_xyyyzz, g_z_0_xyzzz_xyyyzzz, g_z_0_xyzzz_xyyzzz, g_z_0_xyzzz_xyyzzzz, g_z_0_xyzzz_xyzzzz, g_z_0_xyzzz_xyzzzzz, g_z_0_xyzzz_xzzzzz, g_z_0_xyzzz_xzzzzzz, g_z_0_xyzzz_yyyyyy, g_z_0_xyzzz_yyyyyz, g_z_0_xyzzz_yyyyzz, g_z_0_xyzzz_yyyzzz, g_z_0_xyzzz_yyzzzz, g_z_0_xyzzz_yzzzzz, g_z_0_xyzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_xxxxxx[k] = -g_z_0_xyzzz_xxxxxx[k] * ab_x + g_z_0_xyzzz_xxxxxxx[k];

                g_z_0_xxyzzz_xxxxxy[k] = -g_z_0_xyzzz_xxxxxy[k] * ab_x + g_z_0_xyzzz_xxxxxxy[k];

                g_z_0_xxyzzz_xxxxxz[k] = -g_z_0_xyzzz_xxxxxz[k] * ab_x + g_z_0_xyzzz_xxxxxxz[k];

                g_z_0_xxyzzz_xxxxyy[k] = -g_z_0_xyzzz_xxxxyy[k] * ab_x + g_z_0_xyzzz_xxxxxyy[k];

                g_z_0_xxyzzz_xxxxyz[k] = -g_z_0_xyzzz_xxxxyz[k] * ab_x + g_z_0_xyzzz_xxxxxyz[k];

                g_z_0_xxyzzz_xxxxzz[k] = -g_z_0_xyzzz_xxxxzz[k] * ab_x + g_z_0_xyzzz_xxxxxzz[k];

                g_z_0_xxyzzz_xxxyyy[k] = -g_z_0_xyzzz_xxxyyy[k] * ab_x + g_z_0_xyzzz_xxxxyyy[k];

                g_z_0_xxyzzz_xxxyyz[k] = -g_z_0_xyzzz_xxxyyz[k] * ab_x + g_z_0_xyzzz_xxxxyyz[k];

                g_z_0_xxyzzz_xxxyzz[k] = -g_z_0_xyzzz_xxxyzz[k] * ab_x + g_z_0_xyzzz_xxxxyzz[k];

                g_z_0_xxyzzz_xxxzzz[k] = -g_z_0_xyzzz_xxxzzz[k] * ab_x + g_z_0_xyzzz_xxxxzzz[k];

                g_z_0_xxyzzz_xxyyyy[k] = -g_z_0_xyzzz_xxyyyy[k] * ab_x + g_z_0_xyzzz_xxxyyyy[k];

                g_z_0_xxyzzz_xxyyyz[k] = -g_z_0_xyzzz_xxyyyz[k] * ab_x + g_z_0_xyzzz_xxxyyyz[k];

                g_z_0_xxyzzz_xxyyzz[k] = -g_z_0_xyzzz_xxyyzz[k] * ab_x + g_z_0_xyzzz_xxxyyzz[k];

                g_z_0_xxyzzz_xxyzzz[k] = -g_z_0_xyzzz_xxyzzz[k] * ab_x + g_z_0_xyzzz_xxxyzzz[k];

                g_z_0_xxyzzz_xxzzzz[k] = -g_z_0_xyzzz_xxzzzz[k] * ab_x + g_z_0_xyzzz_xxxzzzz[k];

                g_z_0_xxyzzz_xyyyyy[k] = -g_z_0_xyzzz_xyyyyy[k] * ab_x + g_z_0_xyzzz_xxyyyyy[k];

                g_z_0_xxyzzz_xyyyyz[k] = -g_z_0_xyzzz_xyyyyz[k] * ab_x + g_z_0_xyzzz_xxyyyyz[k];

                g_z_0_xxyzzz_xyyyzz[k] = -g_z_0_xyzzz_xyyyzz[k] * ab_x + g_z_0_xyzzz_xxyyyzz[k];

                g_z_0_xxyzzz_xyyzzz[k] = -g_z_0_xyzzz_xyyzzz[k] * ab_x + g_z_0_xyzzz_xxyyzzz[k];

                g_z_0_xxyzzz_xyzzzz[k] = -g_z_0_xyzzz_xyzzzz[k] * ab_x + g_z_0_xyzzz_xxyzzzz[k];

                g_z_0_xxyzzz_xzzzzz[k] = -g_z_0_xyzzz_xzzzzz[k] * ab_x + g_z_0_xyzzz_xxzzzzz[k];

                g_z_0_xxyzzz_yyyyyy[k] = -g_z_0_xyzzz_yyyyyy[k] * ab_x + g_z_0_xyzzz_xyyyyyy[k];

                g_z_0_xxyzzz_yyyyyz[k] = -g_z_0_xyzzz_yyyyyz[k] * ab_x + g_z_0_xyzzz_xyyyyyz[k];

                g_z_0_xxyzzz_yyyyzz[k] = -g_z_0_xyzzz_yyyyzz[k] * ab_x + g_z_0_xyzzz_xyyyyzz[k];

                g_z_0_xxyzzz_yyyzzz[k] = -g_z_0_xyzzz_yyyzzz[k] * ab_x + g_z_0_xyzzz_xyyyzzz[k];

                g_z_0_xxyzzz_yyzzzz[k] = -g_z_0_xyzzz_yyzzzz[k] * ab_x + g_z_0_xyzzz_xyyzzzz[k];

                g_z_0_xxyzzz_yzzzzz[k] = -g_z_0_xyzzz_yzzzzz[k] * ab_x + g_z_0_xyzzz_xyzzzzz[k];

                g_z_0_xxyzzz_zzzzzz[k] = -g_z_0_xyzzz_zzzzzz[k] * ab_x + g_z_0_xyzzz_xzzzzzz[k];
            }

            /// Set up 1960-1988 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 1960 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 1961 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 1962 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 1963 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 1964 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 1965 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 1966 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 1967 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 1968 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 1969 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 1970 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 1971 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 1972 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 1973 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 1974 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 1975 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 1976 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 1977 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 1978 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 1979 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 1980 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 1981 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 1982 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 1983 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 1984 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 1985 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 1986 * ccomps * dcomps);

            auto g_z_0_xxzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 1987 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzzzz_xxxxxx, g_z_0_xxzzzz_xxxxxy, g_z_0_xxzzzz_xxxxxz, g_z_0_xxzzzz_xxxxyy, g_z_0_xxzzzz_xxxxyz, g_z_0_xxzzzz_xxxxzz, g_z_0_xxzzzz_xxxyyy, g_z_0_xxzzzz_xxxyyz, g_z_0_xxzzzz_xxxyzz, g_z_0_xxzzzz_xxxzzz, g_z_0_xxzzzz_xxyyyy, g_z_0_xxzzzz_xxyyyz, g_z_0_xxzzzz_xxyyzz, g_z_0_xxzzzz_xxyzzz, g_z_0_xxzzzz_xxzzzz, g_z_0_xxzzzz_xyyyyy, g_z_0_xxzzzz_xyyyyz, g_z_0_xxzzzz_xyyyzz, g_z_0_xxzzzz_xyyzzz, g_z_0_xxzzzz_xyzzzz, g_z_0_xxzzzz_xzzzzz, g_z_0_xxzzzz_yyyyyy, g_z_0_xxzzzz_yyyyyz, g_z_0_xxzzzz_yyyyzz, g_z_0_xxzzzz_yyyzzz, g_z_0_xxzzzz_yyzzzz, g_z_0_xxzzzz_yzzzzz, g_z_0_xxzzzz_zzzzzz, g_z_0_xzzzz_xxxxxx, g_z_0_xzzzz_xxxxxxx, g_z_0_xzzzz_xxxxxxy, g_z_0_xzzzz_xxxxxxz, g_z_0_xzzzz_xxxxxy, g_z_0_xzzzz_xxxxxyy, g_z_0_xzzzz_xxxxxyz, g_z_0_xzzzz_xxxxxz, g_z_0_xzzzz_xxxxxzz, g_z_0_xzzzz_xxxxyy, g_z_0_xzzzz_xxxxyyy, g_z_0_xzzzz_xxxxyyz, g_z_0_xzzzz_xxxxyz, g_z_0_xzzzz_xxxxyzz, g_z_0_xzzzz_xxxxzz, g_z_0_xzzzz_xxxxzzz, g_z_0_xzzzz_xxxyyy, g_z_0_xzzzz_xxxyyyy, g_z_0_xzzzz_xxxyyyz, g_z_0_xzzzz_xxxyyz, g_z_0_xzzzz_xxxyyzz, g_z_0_xzzzz_xxxyzz, g_z_0_xzzzz_xxxyzzz, g_z_0_xzzzz_xxxzzz, g_z_0_xzzzz_xxxzzzz, g_z_0_xzzzz_xxyyyy, g_z_0_xzzzz_xxyyyyy, g_z_0_xzzzz_xxyyyyz, g_z_0_xzzzz_xxyyyz, g_z_0_xzzzz_xxyyyzz, g_z_0_xzzzz_xxyyzz, g_z_0_xzzzz_xxyyzzz, g_z_0_xzzzz_xxyzzz, g_z_0_xzzzz_xxyzzzz, g_z_0_xzzzz_xxzzzz, g_z_0_xzzzz_xxzzzzz, g_z_0_xzzzz_xyyyyy, g_z_0_xzzzz_xyyyyyy, g_z_0_xzzzz_xyyyyyz, g_z_0_xzzzz_xyyyyz, g_z_0_xzzzz_xyyyyzz, g_z_0_xzzzz_xyyyzz, g_z_0_xzzzz_xyyyzzz, g_z_0_xzzzz_xyyzzz, g_z_0_xzzzz_xyyzzzz, g_z_0_xzzzz_xyzzzz, g_z_0_xzzzz_xyzzzzz, g_z_0_xzzzz_xzzzzz, g_z_0_xzzzz_xzzzzzz, g_z_0_xzzzz_yyyyyy, g_z_0_xzzzz_yyyyyz, g_z_0_xzzzz_yyyyzz, g_z_0_xzzzz_yyyzzz, g_z_0_xzzzz_yyzzzz, g_z_0_xzzzz_yzzzzz, g_z_0_xzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_xxxxxx[k] = -g_z_0_xzzzz_xxxxxx[k] * ab_x + g_z_0_xzzzz_xxxxxxx[k];

                g_z_0_xxzzzz_xxxxxy[k] = -g_z_0_xzzzz_xxxxxy[k] * ab_x + g_z_0_xzzzz_xxxxxxy[k];

                g_z_0_xxzzzz_xxxxxz[k] = -g_z_0_xzzzz_xxxxxz[k] * ab_x + g_z_0_xzzzz_xxxxxxz[k];

                g_z_0_xxzzzz_xxxxyy[k] = -g_z_0_xzzzz_xxxxyy[k] * ab_x + g_z_0_xzzzz_xxxxxyy[k];

                g_z_0_xxzzzz_xxxxyz[k] = -g_z_0_xzzzz_xxxxyz[k] * ab_x + g_z_0_xzzzz_xxxxxyz[k];

                g_z_0_xxzzzz_xxxxzz[k] = -g_z_0_xzzzz_xxxxzz[k] * ab_x + g_z_0_xzzzz_xxxxxzz[k];

                g_z_0_xxzzzz_xxxyyy[k] = -g_z_0_xzzzz_xxxyyy[k] * ab_x + g_z_0_xzzzz_xxxxyyy[k];

                g_z_0_xxzzzz_xxxyyz[k] = -g_z_0_xzzzz_xxxyyz[k] * ab_x + g_z_0_xzzzz_xxxxyyz[k];

                g_z_0_xxzzzz_xxxyzz[k] = -g_z_0_xzzzz_xxxyzz[k] * ab_x + g_z_0_xzzzz_xxxxyzz[k];

                g_z_0_xxzzzz_xxxzzz[k] = -g_z_0_xzzzz_xxxzzz[k] * ab_x + g_z_0_xzzzz_xxxxzzz[k];

                g_z_0_xxzzzz_xxyyyy[k] = -g_z_0_xzzzz_xxyyyy[k] * ab_x + g_z_0_xzzzz_xxxyyyy[k];

                g_z_0_xxzzzz_xxyyyz[k] = -g_z_0_xzzzz_xxyyyz[k] * ab_x + g_z_0_xzzzz_xxxyyyz[k];

                g_z_0_xxzzzz_xxyyzz[k] = -g_z_0_xzzzz_xxyyzz[k] * ab_x + g_z_0_xzzzz_xxxyyzz[k];

                g_z_0_xxzzzz_xxyzzz[k] = -g_z_0_xzzzz_xxyzzz[k] * ab_x + g_z_0_xzzzz_xxxyzzz[k];

                g_z_0_xxzzzz_xxzzzz[k] = -g_z_0_xzzzz_xxzzzz[k] * ab_x + g_z_0_xzzzz_xxxzzzz[k];

                g_z_0_xxzzzz_xyyyyy[k] = -g_z_0_xzzzz_xyyyyy[k] * ab_x + g_z_0_xzzzz_xxyyyyy[k];

                g_z_0_xxzzzz_xyyyyz[k] = -g_z_0_xzzzz_xyyyyz[k] * ab_x + g_z_0_xzzzz_xxyyyyz[k];

                g_z_0_xxzzzz_xyyyzz[k] = -g_z_0_xzzzz_xyyyzz[k] * ab_x + g_z_0_xzzzz_xxyyyzz[k];

                g_z_0_xxzzzz_xyyzzz[k] = -g_z_0_xzzzz_xyyzzz[k] * ab_x + g_z_0_xzzzz_xxyyzzz[k];

                g_z_0_xxzzzz_xyzzzz[k] = -g_z_0_xzzzz_xyzzzz[k] * ab_x + g_z_0_xzzzz_xxyzzzz[k];

                g_z_0_xxzzzz_xzzzzz[k] = -g_z_0_xzzzz_xzzzzz[k] * ab_x + g_z_0_xzzzz_xxzzzzz[k];

                g_z_0_xxzzzz_yyyyyy[k] = -g_z_0_xzzzz_yyyyyy[k] * ab_x + g_z_0_xzzzz_xyyyyyy[k];

                g_z_0_xxzzzz_yyyyyz[k] = -g_z_0_xzzzz_yyyyyz[k] * ab_x + g_z_0_xzzzz_xyyyyyz[k];

                g_z_0_xxzzzz_yyyyzz[k] = -g_z_0_xzzzz_yyyyzz[k] * ab_x + g_z_0_xzzzz_xyyyyzz[k];

                g_z_0_xxzzzz_yyyzzz[k] = -g_z_0_xzzzz_yyyzzz[k] * ab_x + g_z_0_xzzzz_xyyyzzz[k];

                g_z_0_xxzzzz_yyzzzz[k] = -g_z_0_xzzzz_yyzzzz[k] * ab_x + g_z_0_xzzzz_xyyzzzz[k];

                g_z_0_xxzzzz_yzzzzz[k] = -g_z_0_xzzzz_yzzzzz[k] * ab_x + g_z_0_xzzzz_xyzzzzz[k];

                g_z_0_xxzzzz_zzzzzz[k] = -g_z_0_xzzzz_zzzzzz[k] * ab_x + g_z_0_xzzzz_xzzzzzz[k];
            }

            /// Set up 1988-2016 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 1988 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 1989 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 1990 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 1991 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 1992 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 1993 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 1994 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 1995 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 1996 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 1997 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 1998 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 1999 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 2000 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 2001 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 2002 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 2003 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 2004 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 2005 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 2006 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 2007 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 2008 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 2009 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 2010 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 2011 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 2012 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 2013 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 2014 * ccomps * dcomps);

            auto g_z_0_xyyyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 2015 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyyy_xxxxxx, g_z_0_xyyyyy_xxxxxy, g_z_0_xyyyyy_xxxxxz, g_z_0_xyyyyy_xxxxyy, g_z_0_xyyyyy_xxxxyz, g_z_0_xyyyyy_xxxxzz, g_z_0_xyyyyy_xxxyyy, g_z_0_xyyyyy_xxxyyz, g_z_0_xyyyyy_xxxyzz, g_z_0_xyyyyy_xxxzzz, g_z_0_xyyyyy_xxyyyy, g_z_0_xyyyyy_xxyyyz, g_z_0_xyyyyy_xxyyzz, g_z_0_xyyyyy_xxyzzz, g_z_0_xyyyyy_xxzzzz, g_z_0_xyyyyy_xyyyyy, g_z_0_xyyyyy_xyyyyz, g_z_0_xyyyyy_xyyyzz, g_z_0_xyyyyy_xyyzzz, g_z_0_xyyyyy_xyzzzz, g_z_0_xyyyyy_xzzzzz, g_z_0_xyyyyy_yyyyyy, g_z_0_xyyyyy_yyyyyz, g_z_0_xyyyyy_yyyyzz, g_z_0_xyyyyy_yyyzzz, g_z_0_xyyyyy_yyzzzz, g_z_0_xyyyyy_yzzzzz, g_z_0_xyyyyy_zzzzzz, g_z_0_yyyyy_xxxxxx, g_z_0_yyyyy_xxxxxxx, g_z_0_yyyyy_xxxxxxy, g_z_0_yyyyy_xxxxxxz, g_z_0_yyyyy_xxxxxy, g_z_0_yyyyy_xxxxxyy, g_z_0_yyyyy_xxxxxyz, g_z_0_yyyyy_xxxxxz, g_z_0_yyyyy_xxxxxzz, g_z_0_yyyyy_xxxxyy, g_z_0_yyyyy_xxxxyyy, g_z_0_yyyyy_xxxxyyz, g_z_0_yyyyy_xxxxyz, g_z_0_yyyyy_xxxxyzz, g_z_0_yyyyy_xxxxzz, g_z_0_yyyyy_xxxxzzz, g_z_0_yyyyy_xxxyyy, g_z_0_yyyyy_xxxyyyy, g_z_0_yyyyy_xxxyyyz, g_z_0_yyyyy_xxxyyz, g_z_0_yyyyy_xxxyyzz, g_z_0_yyyyy_xxxyzz, g_z_0_yyyyy_xxxyzzz, g_z_0_yyyyy_xxxzzz, g_z_0_yyyyy_xxxzzzz, g_z_0_yyyyy_xxyyyy, g_z_0_yyyyy_xxyyyyy, g_z_0_yyyyy_xxyyyyz, g_z_0_yyyyy_xxyyyz, g_z_0_yyyyy_xxyyyzz, g_z_0_yyyyy_xxyyzz, g_z_0_yyyyy_xxyyzzz, g_z_0_yyyyy_xxyzzz, g_z_0_yyyyy_xxyzzzz, g_z_0_yyyyy_xxzzzz, g_z_0_yyyyy_xxzzzzz, g_z_0_yyyyy_xyyyyy, g_z_0_yyyyy_xyyyyyy, g_z_0_yyyyy_xyyyyyz, g_z_0_yyyyy_xyyyyz, g_z_0_yyyyy_xyyyyzz, g_z_0_yyyyy_xyyyzz, g_z_0_yyyyy_xyyyzzz, g_z_0_yyyyy_xyyzzz, g_z_0_yyyyy_xyyzzzz, g_z_0_yyyyy_xyzzzz, g_z_0_yyyyy_xyzzzzz, g_z_0_yyyyy_xzzzzz, g_z_0_yyyyy_xzzzzzz, g_z_0_yyyyy_yyyyyy, g_z_0_yyyyy_yyyyyz, g_z_0_yyyyy_yyyyzz, g_z_0_yyyyy_yyyzzz, g_z_0_yyyyy_yyzzzz, g_z_0_yyyyy_yzzzzz, g_z_0_yyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_xxxxxx[k] = -g_z_0_yyyyy_xxxxxx[k] * ab_x + g_z_0_yyyyy_xxxxxxx[k];

                g_z_0_xyyyyy_xxxxxy[k] = -g_z_0_yyyyy_xxxxxy[k] * ab_x + g_z_0_yyyyy_xxxxxxy[k];

                g_z_0_xyyyyy_xxxxxz[k] = -g_z_0_yyyyy_xxxxxz[k] * ab_x + g_z_0_yyyyy_xxxxxxz[k];

                g_z_0_xyyyyy_xxxxyy[k] = -g_z_0_yyyyy_xxxxyy[k] * ab_x + g_z_0_yyyyy_xxxxxyy[k];

                g_z_0_xyyyyy_xxxxyz[k] = -g_z_0_yyyyy_xxxxyz[k] * ab_x + g_z_0_yyyyy_xxxxxyz[k];

                g_z_0_xyyyyy_xxxxzz[k] = -g_z_0_yyyyy_xxxxzz[k] * ab_x + g_z_0_yyyyy_xxxxxzz[k];

                g_z_0_xyyyyy_xxxyyy[k] = -g_z_0_yyyyy_xxxyyy[k] * ab_x + g_z_0_yyyyy_xxxxyyy[k];

                g_z_0_xyyyyy_xxxyyz[k] = -g_z_0_yyyyy_xxxyyz[k] * ab_x + g_z_0_yyyyy_xxxxyyz[k];

                g_z_0_xyyyyy_xxxyzz[k] = -g_z_0_yyyyy_xxxyzz[k] * ab_x + g_z_0_yyyyy_xxxxyzz[k];

                g_z_0_xyyyyy_xxxzzz[k] = -g_z_0_yyyyy_xxxzzz[k] * ab_x + g_z_0_yyyyy_xxxxzzz[k];

                g_z_0_xyyyyy_xxyyyy[k] = -g_z_0_yyyyy_xxyyyy[k] * ab_x + g_z_0_yyyyy_xxxyyyy[k];

                g_z_0_xyyyyy_xxyyyz[k] = -g_z_0_yyyyy_xxyyyz[k] * ab_x + g_z_0_yyyyy_xxxyyyz[k];

                g_z_0_xyyyyy_xxyyzz[k] = -g_z_0_yyyyy_xxyyzz[k] * ab_x + g_z_0_yyyyy_xxxyyzz[k];

                g_z_0_xyyyyy_xxyzzz[k] = -g_z_0_yyyyy_xxyzzz[k] * ab_x + g_z_0_yyyyy_xxxyzzz[k];

                g_z_0_xyyyyy_xxzzzz[k] = -g_z_0_yyyyy_xxzzzz[k] * ab_x + g_z_0_yyyyy_xxxzzzz[k];

                g_z_0_xyyyyy_xyyyyy[k] = -g_z_0_yyyyy_xyyyyy[k] * ab_x + g_z_0_yyyyy_xxyyyyy[k];

                g_z_0_xyyyyy_xyyyyz[k] = -g_z_0_yyyyy_xyyyyz[k] * ab_x + g_z_0_yyyyy_xxyyyyz[k];

                g_z_0_xyyyyy_xyyyzz[k] = -g_z_0_yyyyy_xyyyzz[k] * ab_x + g_z_0_yyyyy_xxyyyzz[k];

                g_z_0_xyyyyy_xyyzzz[k] = -g_z_0_yyyyy_xyyzzz[k] * ab_x + g_z_0_yyyyy_xxyyzzz[k];

                g_z_0_xyyyyy_xyzzzz[k] = -g_z_0_yyyyy_xyzzzz[k] * ab_x + g_z_0_yyyyy_xxyzzzz[k];

                g_z_0_xyyyyy_xzzzzz[k] = -g_z_0_yyyyy_xzzzzz[k] * ab_x + g_z_0_yyyyy_xxzzzzz[k];

                g_z_0_xyyyyy_yyyyyy[k] = -g_z_0_yyyyy_yyyyyy[k] * ab_x + g_z_0_yyyyy_xyyyyyy[k];

                g_z_0_xyyyyy_yyyyyz[k] = -g_z_0_yyyyy_yyyyyz[k] * ab_x + g_z_0_yyyyy_xyyyyyz[k];

                g_z_0_xyyyyy_yyyyzz[k] = -g_z_0_yyyyy_yyyyzz[k] * ab_x + g_z_0_yyyyy_xyyyyzz[k];

                g_z_0_xyyyyy_yyyzzz[k] = -g_z_0_yyyyy_yyyzzz[k] * ab_x + g_z_0_yyyyy_xyyyzzz[k];

                g_z_0_xyyyyy_yyzzzz[k] = -g_z_0_yyyyy_yyzzzz[k] * ab_x + g_z_0_yyyyy_xyyzzzz[k];

                g_z_0_xyyyyy_yzzzzz[k] = -g_z_0_yyyyy_yzzzzz[k] * ab_x + g_z_0_yyyyy_xyzzzzz[k];

                g_z_0_xyyyyy_zzzzzz[k] = -g_z_0_yyyyy_zzzzzz[k] * ab_x + g_z_0_yyyyy_xzzzzzz[k];
            }

            /// Set up 2016-2044 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 2016 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 2017 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 2018 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 2019 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 2020 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 2021 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 2022 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 2023 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 2024 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 2025 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 2026 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 2027 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 2028 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 2029 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 2030 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 2031 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 2032 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 2033 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 2034 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 2035 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 2036 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 2037 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 2038 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 2039 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 2040 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 2041 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 2042 * ccomps * dcomps);

            auto g_z_0_xyyyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 2043 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyyz_xxxxxx, g_z_0_xyyyyz_xxxxxy, g_z_0_xyyyyz_xxxxxz, g_z_0_xyyyyz_xxxxyy, g_z_0_xyyyyz_xxxxyz, g_z_0_xyyyyz_xxxxzz, g_z_0_xyyyyz_xxxyyy, g_z_0_xyyyyz_xxxyyz, g_z_0_xyyyyz_xxxyzz, g_z_0_xyyyyz_xxxzzz, g_z_0_xyyyyz_xxyyyy, g_z_0_xyyyyz_xxyyyz, g_z_0_xyyyyz_xxyyzz, g_z_0_xyyyyz_xxyzzz, g_z_0_xyyyyz_xxzzzz, g_z_0_xyyyyz_xyyyyy, g_z_0_xyyyyz_xyyyyz, g_z_0_xyyyyz_xyyyzz, g_z_0_xyyyyz_xyyzzz, g_z_0_xyyyyz_xyzzzz, g_z_0_xyyyyz_xzzzzz, g_z_0_xyyyyz_yyyyyy, g_z_0_xyyyyz_yyyyyz, g_z_0_xyyyyz_yyyyzz, g_z_0_xyyyyz_yyyzzz, g_z_0_xyyyyz_yyzzzz, g_z_0_xyyyyz_yzzzzz, g_z_0_xyyyyz_zzzzzz, g_z_0_yyyyz_xxxxxx, g_z_0_yyyyz_xxxxxxx, g_z_0_yyyyz_xxxxxxy, g_z_0_yyyyz_xxxxxxz, g_z_0_yyyyz_xxxxxy, g_z_0_yyyyz_xxxxxyy, g_z_0_yyyyz_xxxxxyz, g_z_0_yyyyz_xxxxxz, g_z_0_yyyyz_xxxxxzz, g_z_0_yyyyz_xxxxyy, g_z_0_yyyyz_xxxxyyy, g_z_0_yyyyz_xxxxyyz, g_z_0_yyyyz_xxxxyz, g_z_0_yyyyz_xxxxyzz, g_z_0_yyyyz_xxxxzz, g_z_0_yyyyz_xxxxzzz, g_z_0_yyyyz_xxxyyy, g_z_0_yyyyz_xxxyyyy, g_z_0_yyyyz_xxxyyyz, g_z_0_yyyyz_xxxyyz, g_z_0_yyyyz_xxxyyzz, g_z_0_yyyyz_xxxyzz, g_z_0_yyyyz_xxxyzzz, g_z_0_yyyyz_xxxzzz, g_z_0_yyyyz_xxxzzzz, g_z_0_yyyyz_xxyyyy, g_z_0_yyyyz_xxyyyyy, g_z_0_yyyyz_xxyyyyz, g_z_0_yyyyz_xxyyyz, g_z_0_yyyyz_xxyyyzz, g_z_0_yyyyz_xxyyzz, g_z_0_yyyyz_xxyyzzz, g_z_0_yyyyz_xxyzzz, g_z_0_yyyyz_xxyzzzz, g_z_0_yyyyz_xxzzzz, g_z_0_yyyyz_xxzzzzz, g_z_0_yyyyz_xyyyyy, g_z_0_yyyyz_xyyyyyy, g_z_0_yyyyz_xyyyyyz, g_z_0_yyyyz_xyyyyz, g_z_0_yyyyz_xyyyyzz, g_z_0_yyyyz_xyyyzz, g_z_0_yyyyz_xyyyzzz, g_z_0_yyyyz_xyyzzz, g_z_0_yyyyz_xyyzzzz, g_z_0_yyyyz_xyzzzz, g_z_0_yyyyz_xyzzzzz, g_z_0_yyyyz_xzzzzz, g_z_0_yyyyz_xzzzzzz, g_z_0_yyyyz_yyyyyy, g_z_0_yyyyz_yyyyyz, g_z_0_yyyyz_yyyyzz, g_z_0_yyyyz_yyyzzz, g_z_0_yyyyz_yyzzzz, g_z_0_yyyyz_yzzzzz, g_z_0_yyyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_xxxxxx[k] = -g_z_0_yyyyz_xxxxxx[k] * ab_x + g_z_0_yyyyz_xxxxxxx[k];

                g_z_0_xyyyyz_xxxxxy[k] = -g_z_0_yyyyz_xxxxxy[k] * ab_x + g_z_0_yyyyz_xxxxxxy[k];

                g_z_0_xyyyyz_xxxxxz[k] = -g_z_0_yyyyz_xxxxxz[k] * ab_x + g_z_0_yyyyz_xxxxxxz[k];

                g_z_0_xyyyyz_xxxxyy[k] = -g_z_0_yyyyz_xxxxyy[k] * ab_x + g_z_0_yyyyz_xxxxxyy[k];

                g_z_0_xyyyyz_xxxxyz[k] = -g_z_0_yyyyz_xxxxyz[k] * ab_x + g_z_0_yyyyz_xxxxxyz[k];

                g_z_0_xyyyyz_xxxxzz[k] = -g_z_0_yyyyz_xxxxzz[k] * ab_x + g_z_0_yyyyz_xxxxxzz[k];

                g_z_0_xyyyyz_xxxyyy[k] = -g_z_0_yyyyz_xxxyyy[k] * ab_x + g_z_0_yyyyz_xxxxyyy[k];

                g_z_0_xyyyyz_xxxyyz[k] = -g_z_0_yyyyz_xxxyyz[k] * ab_x + g_z_0_yyyyz_xxxxyyz[k];

                g_z_0_xyyyyz_xxxyzz[k] = -g_z_0_yyyyz_xxxyzz[k] * ab_x + g_z_0_yyyyz_xxxxyzz[k];

                g_z_0_xyyyyz_xxxzzz[k] = -g_z_0_yyyyz_xxxzzz[k] * ab_x + g_z_0_yyyyz_xxxxzzz[k];

                g_z_0_xyyyyz_xxyyyy[k] = -g_z_0_yyyyz_xxyyyy[k] * ab_x + g_z_0_yyyyz_xxxyyyy[k];

                g_z_0_xyyyyz_xxyyyz[k] = -g_z_0_yyyyz_xxyyyz[k] * ab_x + g_z_0_yyyyz_xxxyyyz[k];

                g_z_0_xyyyyz_xxyyzz[k] = -g_z_0_yyyyz_xxyyzz[k] * ab_x + g_z_0_yyyyz_xxxyyzz[k];

                g_z_0_xyyyyz_xxyzzz[k] = -g_z_0_yyyyz_xxyzzz[k] * ab_x + g_z_0_yyyyz_xxxyzzz[k];

                g_z_0_xyyyyz_xxzzzz[k] = -g_z_0_yyyyz_xxzzzz[k] * ab_x + g_z_0_yyyyz_xxxzzzz[k];

                g_z_0_xyyyyz_xyyyyy[k] = -g_z_0_yyyyz_xyyyyy[k] * ab_x + g_z_0_yyyyz_xxyyyyy[k];

                g_z_0_xyyyyz_xyyyyz[k] = -g_z_0_yyyyz_xyyyyz[k] * ab_x + g_z_0_yyyyz_xxyyyyz[k];

                g_z_0_xyyyyz_xyyyzz[k] = -g_z_0_yyyyz_xyyyzz[k] * ab_x + g_z_0_yyyyz_xxyyyzz[k];

                g_z_0_xyyyyz_xyyzzz[k] = -g_z_0_yyyyz_xyyzzz[k] * ab_x + g_z_0_yyyyz_xxyyzzz[k];

                g_z_0_xyyyyz_xyzzzz[k] = -g_z_0_yyyyz_xyzzzz[k] * ab_x + g_z_0_yyyyz_xxyzzzz[k];

                g_z_0_xyyyyz_xzzzzz[k] = -g_z_0_yyyyz_xzzzzz[k] * ab_x + g_z_0_yyyyz_xxzzzzz[k];

                g_z_0_xyyyyz_yyyyyy[k] = -g_z_0_yyyyz_yyyyyy[k] * ab_x + g_z_0_yyyyz_xyyyyyy[k];

                g_z_0_xyyyyz_yyyyyz[k] = -g_z_0_yyyyz_yyyyyz[k] * ab_x + g_z_0_yyyyz_xyyyyyz[k];

                g_z_0_xyyyyz_yyyyzz[k] = -g_z_0_yyyyz_yyyyzz[k] * ab_x + g_z_0_yyyyz_xyyyyzz[k];

                g_z_0_xyyyyz_yyyzzz[k] = -g_z_0_yyyyz_yyyzzz[k] * ab_x + g_z_0_yyyyz_xyyyzzz[k];

                g_z_0_xyyyyz_yyzzzz[k] = -g_z_0_yyyyz_yyzzzz[k] * ab_x + g_z_0_yyyyz_xyyzzzz[k];

                g_z_0_xyyyyz_yzzzzz[k] = -g_z_0_yyyyz_yzzzzz[k] * ab_x + g_z_0_yyyyz_xyzzzzz[k];

                g_z_0_xyyyyz_zzzzzz[k] = -g_z_0_yyyyz_zzzzzz[k] * ab_x + g_z_0_yyyyz_xzzzzzz[k];
            }

            /// Set up 2044-2072 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 2044 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 2045 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 2046 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 2047 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 2048 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 2049 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 2050 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 2051 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 2052 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 2053 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 2054 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 2055 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 2056 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 2057 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 2058 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 2059 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 2060 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 2061 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 2062 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 2063 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 2064 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 2065 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 2066 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 2067 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 2068 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 2069 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 2070 * ccomps * dcomps);

            auto g_z_0_xyyyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 2071 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyzz_xxxxxx, g_z_0_xyyyzz_xxxxxy, g_z_0_xyyyzz_xxxxxz, g_z_0_xyyyzz_xxxxyy, g_z_0_xyyyzz_xxxxyz, g_z_0_xyyyzz_xxxxzz, g_z_0_xyyyzz_xxxyyy, g_z_0_xyyyzz_xxxyyz, g_z_0_xyyyzz_xxxyzz, g_z_0_xyyyzz_xxxzzz, g_z_0_xyyyzz_xxyyyy, g_z_0_xyyyzz_xxyyyz, g_z_0_xyyyzz_xxyyzz, g_z_0_xyyyzz_xxyzzz, g_z_0_xyyyzz_xxzzzz, g_z_0_xyyyzz_xyyyyy, g_z_0_xyyyzz_xyyyyz, g_z_0_xyyyzz_xyyyzz, g_z_0_xyyyzz_xyyzzz, g_z_0_xyyyzz_xyzzzz, g_z_0_xyyyzz_xzzzzz, g_z_0_xyyyzz_yyyyyy, g_z_0_xyyyzz_yyyyyz, g_z_0_xyyyzz_yyyyzz, g_z_0_xyyyzz_yyyzzz, g_z_0_xyyyzz_yyzzzz, g_z_0_xyyyzz_yzzzzz, g_z_0_xyyyzz_zzzzzz, g_z_0_yyyzz_xxxxxx, g_z_0_yyyzz_xxxxxxx, g_z_0_yyyzz_xxxxxxy, g_z_0_yyyzz_xxxxxxz, g_z_0_yyyzz_xxxxxy, g_z_0_yyyzz_xxxxxyy, g_z_0_yyyzz_xxxxxyz, g_z_0_yyyzz_xxxxxz, g_z_0_yyyzz_xxxxxzz, g_z_0_yyyzz_xxxxyy, g_z_0_yyyzz_xxxxyyy, g_z_0_yyyzz_xxxxyyz, g_z_0_yyyzz_xxxxyz, g_z_0_yyyzz_xxxxyzz, g_z_0_yyyzz_xxxxzz, g_z_0_yyyzz_xxxxzzz, g_z_0_yyyzz_xxxyyy, g_z_0_yyyzz_xxxyyyy, g_z_0_yyyzz_xxxyyyz, g_z_0_yyyzz_xxxyyz, g_z_0_yyyzz_xxxyyzz, g_z_0_yyyzz_xxxyzz, g_z_0_yyyzz_xxxyzzz, g_z_0_yyyzz_xxxzzz, g_z_0_yyyzz_xxxzzzz, g_z_0_yyyzz_xxyyyy, g_z_0_yyyzz_xxyyyyy, g_z_0_yyyzz_xxyyyyz, g_z_0_yyyzz_xxyyyz, g_z_0_yyyzz_xxyyyzz, g_z_0_yyyzz_xxyyzz, g_z_0_yyyzz_xxyyzzz, g_z_0_yyyzz_xxyzzz, g_z_0_yyyzz_xxyzzzz, g_z_0_yyyzz_xxzzzz, g_z_0_yyyzz_xxzzzzz, g_z_0_yyyzz_xyyyyy, g_z_0_yyyzz_xyyyyyy, g_z_0_yyyzz_xyyyyyz, g_z_0_yyyzz_xyyyyz, g_z_0_yyyzz_xyyyyzz, g_z_0_yyyzz_xyyyzz, g_z_0_yyyzz_xyyyzzz, g_z_0_yyyzz_xyyzzz, g_z_0_yyyzz_xyyzzzz, g_z_0_yyyzz_xyzzzz, g_z_0_yyyzz_xyzzzzz, g_z_0_yyyzz_xzzzzz, g_z_0_yyyzz_xzzzzzz, g_z_0_yyyzz_yyyyyy, g_z_0_yyyzz_yyyyyz, g_z_0_yyyzz_yyyyzz, g_z_0_yyyzz_yyyzzz, g_z_0_yyyzz_yyzzzz, g_z_0_yyyzz_yzzzzz, g_z_0_yyyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_xxxxxx[k] = -g_z_0_yyyzz_xxxxxx[k] * ab_x + g_z_0_yyyzz_xxxxxxx[k];

                g_z_0_xyyyzz_xxxxxy[k] = -g_z_0_yyyzz_xxxxxy[k] * ab_x + g_z_0_yyyzz_xxxxxxy[k];

                g_z_0_xyyyzz_xxxxxz[k] = -g_z_0_yyyzz_xxxxxz[k] * ab_x + g_z_0_yyyzz_xxxxxxz[k];

                g_z_0_xyyyzz_xxxxyy[k] = -g_z_0_yyyzz_xxxxyy[k] * ab_x + g_z_0_yyyzz_xxxxxyy[k];

                g_z_0_xyyyzz_xxxxyz[k] = -g_z_0_yyyzz_xxxxyz[k] * ab_x + g_z_0_yyyzz_xxxxxyz[k];

                g_z_0_xyyyzz_xxxxzz[k] = -g_z_0_yyyzz_xxxxzz[k] * ab_x + g_z_0_yyyzz_xxxxxzz[k];

                g_z_0_xyyyzz_xxxyyy[k] = -g_z_0_yyyzz_xxxyyy[k] * ab_x + g_z_0_yyyzz_xxxxyyy[k];

                g_z_0_xyyyzz_xxxyyz[k] = -g_z_0_yyyzz_xxxyyz[k] * ab_x + g_z_0_yyyzz_xxxxyyz[k];

                g_z_0_xyyyzz_xxxyzz[k] = -g_z_0_yyyzz_xxxyzz[k] * ab_x + g_z_0_yyyzz_xxxxyzz[k];

                g_z_0_xyyyzz_xxxzzz[k] = -g_z_0_yyyzz_xxxzzz[k] * ab_x + g_z_0_yyyzz_xxxxzzz[k];

                g_z_0_xyyyzz_xxyyyy[k] = -g_z_0_yyyzz_xxyyyy[k] * ab_x + g_z_0_yyyzz_xxxyyyy[k];

                g_z_0_xyyyzz_xxyyyz[k] = -g_z_0_yyyzz_xxyyyz[k] * ab_x + g_z_0_yyyzz_xxxyyyz[k];

                g_z_0_xyyyzz_xxyyzz[k] = -g_z_0_yyyzz_xxyyzz[k] * ab_x + g_z_0_yyyzz_xxxyyzz[k];

                g_z_0_xyyyzz_xxyzzz[k] = -g_z_0_yyyzz_xxyzzz[k] * ab_x + g_z_0_yyyzz_xxxyzzz[k];

                g_z_0_xyyyzz_xxzzzz[k] = -g_z_0_yyyzz_xxzzzz[k] * ab_x + g_z_0_yyyzz_xxxzzzz[k];

                g_z_0_xyyyzz_xyyyyy[k] = -g_z_0_yyyzz_xyyyyy[k] * ab_x + g_z_0_yyyzz_xxyyyyy[k];

                g_z_0_xyyyzz_xyyyyz[k] = -g_z_0_yyyzz_xyyyyz[k] * ab_x + g_z_0_yyyzz_xxyyyyz[k];

                g_z_0_xyyyzz_xyyyzz[k] = -g_z_0_yyyzz_xyyyzz[k] * ab_x + g_z_0_yyyzz_xxyyyzz[k];

                g_z_0_xyyyzz_xyyzzz[k] = -g_z_0_yyyzz_xyyzzz[k] * ab_x + g_z_0_yyyzz_xxyyzzz[k];

                g_z_0_xyyyzz_xyzzzz[k] = -g_z_0_yyyzz_xyzzzz[k] * ab_x + g_z_0_yyyzz_xxyzzzz[k];

                g_z_0_xyyyzz_xzzzzz[k] = -g_z_0_yyyzz_xzzzzz[k] * ab_x + g_z_0_yyyzz_xxzzzzz[k];

                g_z_0_xyyyzz_yyyyyy[k] = -g_z_0_yyyzz_yyyyyy[k] * ab_x + g_z_0_yyyzz_xyyyyyy[k];

                g_z_0_xyyyzz_yyyyyz[k] = -g_z_0_yyyzz_yyyyyz[k] * ab_x + g_z_0_yyyzz_xyyyyyz[k];

                g_z_0_xyyyzz_yyyyzz[k] = -g_z_0_yyyzz_yyyyzz[k] * ab_x + g_z_0_yyyzz_xyyyyzz[k];

                g_z_0_xyyyzz_yyyzzz[k] = -g_z_0_yyyzz_yyyzzz[k] * ab_x + g_z_0_yyyzz_xyyyzzz[k];

                g_z_0_xyyyzz_yyzzzz[k] = -g_z_0_yyyzz_yyzzzz[k] * ab_x + g_z_0_yyyzz_xyyzzzz[k];

                g_z_0_xyyyzz_yzzzzz[k] = -g_z_0_yyyzz_yzzzzz[k] * ab_x + g_z_0_yyyzz_xyzzzzz[k];

                g_z_0_xyyyzz_zzzzzz[k] = -g_z_0_yyyzz_zzzzzz[k] * ab_x + g_z_0_yyyzz_xzzzzzz[k];
            }

            /// Set up 2072-2100 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 2072 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 2073 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 2074 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 2075 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 2076 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 2077 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 2078 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 2079 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 2080 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 2081 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 2082 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 2083 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 2084 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 2085 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 2086 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 2087 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 2088 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 2089 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 2090 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 2091 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 2092 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 2093 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 2094 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 2095 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 2096 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 2097 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 2098 * ccomps * dcomps);

            auto g_z_0_xyyzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 2099 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyzzz_xxxxxx, g_z_0_xyyzzz_xxxxxy, g_z_0_xyyzzz_xxxxxz, g_z_0_xyyzzz_xxxxyy, g_z_0_xyyzzz_xxxxyz, g_z_0_xyyzzz_xxxxzz, g_z_0_xyyzzz_xxxyyy, g_z_0_xyyzzz_xxxyyz, g_z_0_xyyzzz_xxxyzz, g_z_0_xyyzzz_xxxzzz, g_z_0_xyyzzz_xxyyyy, g_z_0_xyyzzz_xxyyyz, g_z_0_xyyzzz_xxyyzz, g_z_0_xyyzzz_xxyzzz, g_z_0_xyyzzz_xxzzzz, g_z_0_xyyzzz_xyyyyy, g_z_0_xyyzzz_xyyyyz, g_z_0_xyyzzz_xyyyzz, g_z_0_xyyzzz_xyyzzz, g_z_0_xyyzzz_xyzzzz, g_z_0_xyyzzz_xzzzzz, g_z_0_xyyzzz_yyyyyy, g_z_0_xyyzzz_yyyyyz, g_z_0_xyyzzz_yyyyzz, g_z_0_xyyzzz_yyyzzz, g_z_0_xyyzzz_yyzzzz, g_z_0_xyyzzz_yzzzzz, g_z_0_xyyzzz_zzzzzz, g_z_0_yyzzz_xxxxxx, g_z_0_yyzzz_xxxxxxx, g_z_0_yyzzz_xxxxxxy, g_z_0_yyzzz_xxxxxxz, g_z_0_yyzzz_xxxxxy, g_z_0_yyzzz_xxxxxyy, g_z_0_yyzzz_xxxxxyz, g_z_0_yyzzz_xxxxxz, g_z_0_yyzzz_xxxxxzz, g_z_0_yyzzz_xxxxyy, g_z_0_yyzzz_xxxxyyy, g_z_0_yyzzz_xxxxyyz, g_z_0_yyzzz_xxxxyz, g_z_0_yyzzz_xxxxyzz, g_z_0_yyzzz_xxxxzz, g_z_0_yyzzz_xxxxzzz, g_z_0_yyzzz_xxxyyy, g_z_0_yyzzz_xxxyyyy, g_z_0_yyzzz_xxxyyyz, g_z_0_yyzzz_xxxyyz, g_z_0_yyzzz_xxxyyzz, g_z_0_yyzzz_xxxyzz, g_z_0_yyzzz_xxxyzzz, g_z_0_yyzzz_xxxzzz, g_z_0_yyzzz_xxxzzzz, g_z_0_yyzzz_xxyyyy, g_z_0_yyzzz_xxyyyyy, g_z_0_yyzzz_xxyyyyz, g_z_0_yyzzz_xxyyyz, g_z_0_yyzzz_xxyyyzz, g_z_0_yyzzz_xxyyzz, g_z_0_yyzzz_xxyyzzz, g_z_0_yyzzz_xxyzzz, g_z_0_yyzzz_xxyzzzz, g_z_0_yyzzz_xxzzzz, g_z_0_yyzzz_xxzzzzz, g_z_0_yyzzz_xyyyyy, g_z_0_yyzzz_xyyyyyy, g_z_0_yyzzz_xyyyyyz, g_z_0_yyzzz_xyyyyz, g_z_0_yyzzz_xyyyyzz, g_z_0_yyzzz_xyyyzz, g_z_0_yyzzz_xyyyzzz, g_z_0_yyzzz_xyyzzz, g_z_0_yyzzz_xyyzzzz, g_z_0_yyzzz_xyzzzz, g_z_0_yyzzz_xyzzzzz, g_z_0_yyzzz_xzzzzz, g_z_0_yyzzz_xzzzzzz, g_z_0_yyzzz_yyyyyy, g_z_0_yyzzz_yyyyyz, g_z_0_yyzzz_yyyyzz, g_z_0_yyzzz_yyyzzz, g_z_0_yyzzz_yyzzzz, g_z_0_yyzzz_yzzzzz, g_z_0_yyzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_xxxxxx[k] = -g_z_0_yyzzz_xxxxxx[k] * ab_x + g_z_0_yyzzz_xxxxxxx[k];

                g_z_0_xyyzzz_xxxxxy[k] = -g_z_0_yyzzz_xxxxxy[k] * ab_x + g_z_0_yyzzz_xxxxxxy[k];

                g_z_0_xyyzzz_xxxxxz[k] = -g_z_0_yyzzz_xxxxxz[k] * ab_x + g_z_0_yyzzz_xxxxxxz[k];

                g_z_0_xyyzzz_xxxxyy[k] = -g_z_0_yyzzz_xxxxyy[k] * ab_x + g_z_0_yyzzz_xxxxxyy[k];

                g_z_0_xyyzzz_xxxxyz[k] = -g_z_0_yyzzz_xxxxyz[k] * ab_x + g_z_0_yyzzz_xxxxxyz[k];

                g_z_0_xyyzzz_xxxxzz[k] = -g_z_0_yyzzz_xxxxzz[k] * ab_x + g_z_0_yyzzz_xxxxxzz[k];

                g_z_0_xyyzzz_xxxyyy[k] = -g_z_0_yyzzz_xxxyyy[k] * ab_x + g_z_0_yyzzz_xxxxyyy[k];

                g_z_0_xyyzzz_xxxyyz[k] = -g_z_0_yyzzz_xxxyyz[k] * ab_x + g_z_0_yyzzz_xxxxyyz[k];

                g_z_0_xyyzzz_xxxyzz[k] = -g_z_0_yyzzz_xxxyzz[k] * ab_x + g_z_0_yyzzz_xxxxyzz[k];

                g_z_0_xyyzzz_xxxzzz[k] = -g_z_0_yyzzz_xxxzzz[k] * ab_x + g_z_0_yyzzz_xxxxzzz[k];

                g_z_0_xyyzzz_xxyyyy[k] = -g_z_0_yyzzz_xxyyyy[k] * ab_x + g_z_0_yyzzz_xxxyyyy[k];

                g_z_0_xyyzzz_xxyyyz[k] = -g_z_0_yyzzz_xxyyyz[k] * ab_x + g_z_0_yyzzz_xxxyyyz[k];

                g_z_0_xyyzzz_xxyyzz[k] = -g_z_0_yyzzz_xxyyzz[k] * ab_x + g_z_0_yyzzz_xxxyyzz[k];

                g_z_0_xyyzzz_xxyzzz[k] = -g_z_0_yyzzz_xxyzzz[k] * ab_x + g_z_0_yyzzz_xxxyzzz[k];

                g_z_0_xyyzzz_xxzzzz[k] = -g_z_0_yyzzz_xxzzzz[k] * ab_x + g_z_0_yyzzz_xxxzzzz[k];

                g_z_0_xyyzzz_xyyyyy[k] = -g_z_0_yyzzz_xyyyyy[k] * ab_x + g_z_0_yyzzz_xxyyyyy[k];

                g_z_0_xyyzzz_xyyyyz[k] = -g_z_0_yyzzz_xyyyyz[k] * ab_x + g_z_0_yyzzz_xxyyyyz[k];

                g_z_0_xyyzzz_xyyyzz[k] = -g_z_0_yyzzz_xyyyzz[k] * ab_x + g_z_0_yyzzz_xxyyyzz[k];

                g_z_0_xyyzzz_xyyzzz[k] = -g_z_0_yyzzz_xyyzzz[k] * ab_x + g_z_0_yyzzz_xxyyzzz[k];

                g_z_0_xyyzzz_xyzzzz[k] = -g_z_0_yyzzz_xyzzzz[k] * ab_x + g_z_0_yyzzz_xxyzzzz[k];

                g_z_0_xyyzzz_xzzzzz[k] = -g_z_0_yyzzz_xzzzzz[k] * ab_x + g_z_0_yyzzz_xxzzzzz[k];

                g_z_0_xyyzzz_yyyyyy[k] = -g_z_0_yyzzz_yyyyyy[k] * ab_x + g_z_0_yyzzz_xyyyyyy[k];

                g_z_0_xyyzzz_yyyyyz[k] = -g_z_0_yyzzz_yyyyyz[k] * ab_x + g_z_0_yyzzz_xyyyyyz[k];

                g_z_0_xyyzzz_yyyyzz[k] = -g_z_0_yyzzz_yyyyzz[k] * ab_x + g_z_0_yyzzz_xyyyyzz[k];

                g_z_0_xyyzzz_yyyzzz[k] = -g_z_0_yyzzz_yyyzzz[k] * ab_x + g_z_0_yyzzz_xyyyzzz[k];

                g_z_0_xyyzzz_yyzzzz[k] = -g_z_0_yyzzz_yyzzzz[k] * ab_x + g_z_0_yyzzz_xyyzzzz[k];

                g_z_0_xyyzzz_yzzzzz[k] = -g_z_0_yyzzz_yzzzzz[k] * ab_x + g_z_0_yyzzz_xyzzzzz[k];

                g_z_0_xyyzzz_zzzzzz[k] = -g_z_0_yyzzz_zzzzzz[k] * ab_x + g_z_0_yyzzz_xzzzzzz[k];
            }

            /// Set up 2100-2128 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 2100 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 2101 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 2102 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 2103 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 2104 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 2105 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 2106 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 2107 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 2108 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 2109 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 2110 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 2111 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 2112 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 2113 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 2114 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 2115 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 2116 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 2117 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 2118 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 2119 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 2120 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 2121 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 2122 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 2123 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 2124 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 2125 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 2126 * ccomps * dcomps);

            auto g_z_0_xyzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 2127 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzzzz_xxxxxx, g_z_0_xyzzzz_xxxxxy, g_z_0_xyzzzz_xxxxxz, g_z_0_xyzzzz_xxxxyy, g_z_0_xyzzzz_xxxxyz, g_z_0_xyzzzz_xxxxzz, g_z_0_xyzzzz_xxxyyy, g_z_0_xyzzzz_xxxyyz, g_z_0_xyzzzz_xxxyzz, g_z_0_xyzzzz_xxxzzz, g_z_0_xyzzzz_xxyyyy, g_z_0_xyzzzz_xxyyyz, g_z_0_xyzzzz_xxyyzz, g_z_0_xyzzzz_xxyzzz, g_z_0_xyzzzz_xxzzzz, g_z_0_xyzzzz_xyyyyy, g_z_0_xyzzzz_xyyyyz, g_z_0_xyzzzz_xyyyzz, g_z_0_xyzzzz_xyyzzz, g_z_0_xyzzzz_xyzzzz, g_z_0_xyzzzz_xzzzzz, g_z_0_xyzzzz_yyyyyy, g_z_0_xyzzzz_yyyyyz, g_z_0_xyzzzz_yyyyzz, g_z_0_xyzzzz_yyyzzz, g_z_0_xyzzzz_yyzzzz, g_z_0_xyzzzz_yzzzzz, g_z_0_xyzzzz_zzzzzz, g_z_0_yzzzz_xxxxxx, g_z_0_yzzzz_xxxxxxx, g_z_0_yzzzz_xxxxxxy, g_z_0_yzzzz_xxxxxxz, g_z_0_yzzzz_xxxxxy, g_z_0_yzzzz_xxxxxyy, g_z_0_yzzzz_xxxxxyz, g_z_0_yzzzz_xxxxxz, g_z_0_yzzzz_xxxxxzz, g_z_0_yzzzz_xxxxyy, g_z_0_yzzzz_xxxxyyy, g_z_0_yzzzz_xxxxyyz, g_z_0_yzzzz_xxxxyz, g_z_0_yzzzz_xxxxyzz, g_z_0_yzzzz_xxxxzz, g_z_0_yzzzz_xxxxzzz, g_z_0_yzzzz_xxxyyy, g_z_0_yzzzz_xxxyyyy, g_z_0_yzzzz_xxxyyyz, g_z_0_yzzzz_xxxyyz, g_z_0_yzzzz_xxxyyzz, g_z_0_yzzzz_xxxyzz, g_z_0_yzzzz_xxxyzzz, g_z_0_yzzzz_xxxzzz, g_z_0_yzzzz_xxxzzzz, g_z_0_yzzzz_xxyyyy, g_z_0_yzzzz_xxyyyyy, g_z_0_yzzzz_xxyyyyz, g_z_0_yzzzz_xxyyyz, g_z_0_yzzzz_xxyyyzz, g_z_0_yzzzz_xxyyzz, g_z_0_yzzzz_xxyyzzz, g_z_0_yzzzz_xxyzzz, g_z_0_yzzzz_xxyzzzz, g_z_0_yzzzz_xxzzzz, g_z_0_yzzzz_xxzzzzz, g_z_0_yzzzz_xyyyyy, g_z_0_yzzzz_xyyyyyy, g_z_0_yzzzz_xyyyyyz, g_z_0_yzzzz_xyyyyz, g_z_0_yzzzz_xyyyyzz, g_z_0_yzzzz_xyyyzz, g_z_0_yzzzz_xyyyzzz, g_z_0_yzzzz_xyyzzz, g_z_0_yzzzz_xyyzzzz, g_z_0_yzzzz_xyzzzz, g_z_0_yzzzz_xyzzzzz, g_z_0_yzzzz_xzzzzz, g_z_0_yzzzz_xzzzzzz, g_z_0_yzzzz_yyyyyy, g_z_0_yzzzz_yyyyyz, g_z_0_yzzzz_yyyyzz, g_z_0_yzzzz_yyyzzz, g_z_0_yzzzz_yyzzzz, g_z_0_yzzzz_yzzzzz, g_z_0_yzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_xxxxxx[k] = -g_z_0_yzzzz_xxxxxx[k] * ab_x + g_z_0_yzzzz_xxxxxxx[k];

                g_z_0_xyzzzz_xxxxxy[k] = -g_z_0_yzzzz_xxxxxy[k] * ab_x + g_z_0_yzzzz_xxxxxxy[k];

                g_z_0_xyzzzz_xxxxxz[k] = -g_z_0_yzzzz_xxxxxz[k] * ab_x + g_z_0_yzzzz_xxxxxxz[k];

                g_z_0_xyzzzz_xxxxyy[k] = -g_z_0_yzzzz_xxxxyy[k] * ab_x + g_z_0_yzzzz_xxxxxyy[k];

                g_z_0_xyzzzz_xxxxyz[k] = -g_z_0_yzzzz_xxxxyz[k] * ab_x + g_z_0_yzzzz_xxxxxyz[k];

                g_z_0_xyzzzz_xxxxzz[k] = -g_z_0_yzzzz_xxxxzz[k] * ab_x + g_z_0_yzzzz_xxxxxzz[k];

                g_z_0_xyzzzz_xxxyyy[k] = -g_z_0_yzzzz_xxxyyy[k] * ab_x + g_z_0_yzzzz_xxxxyyy[k];

                g_z_0_xyzzzz_xxxyyz[k] = -g_z_0_yzzzz_xxxyyz[k] * ab_x + g_z_0_yzzzz_xxxxyyz[k];

                g_z_0_xyzzzz_xxxyzz[k] = -g_z_0_yzzzz_xxxyzz[k] * ab_x + g_z_0_yzzzz_xxxxyzz[k];

                g_z_0_xyzzzz_xxxzzz[k] = -g_z_0_yzzzz_xxxzzz[k] * ab_x + g_z_0_yzzzz_xxxxzzz[k];

                g_z_0_xyzzzz_xxyyyy[k] = -g_z_0_yzzzz_xxyyyy[k] * ab_x + g_z_0_yzzzz_xxxyyyy[k];

                g_z_0_xyzzzz_xxyyyz[k] = -g_z_0_yzzzz_xxyyyz[k] * ab_x + g_z_0_yzzzz_xxxyyyz[k];

                g_z_0_xyzzzz_xxyyzz[k] = -g_z_0_yzzzz_xxyyzz[k] * ab_x + g_z_0_yzzzz_xxxyyzz[k];

                g_z_0_xyzzzz_xxyzzz[k] = -g_z_0_yzzzz_xxyzzz[k] * ab_x + g_z_0_yzzzz_xxxyzzz[k];

                g_z_0_xyzzzz_xxzzzz[k] = -g_z_0_yzzzz_xxzzzz[k] * ab_x + g_z_0_yzzzz_xxxzzzz[k];

                g_z_0_xyzzzz_xyyyyy[k] = -g_z_0_yzzzz_xyyyyy[k] * ab_x + g_z_0_yzzzz_xxyyyyy[k];

                g_z_0_xyzzzz_xyyyyz[k] = -g_z_0_yzzzz_xyyyyz[k] * ab_x + g_z_0_yzzzz_xxyyyyz[k];

                g_z_0_xyzzzz_xyyyzz[k] = -g_z_0_yzzzz_xyyyzz[k] * ab_x + g_z_0_yzzzz_xxyyyzz[k];

                g_z_0_xyzzzz_xyyzzz[k] = -g_z_0_yzzzz_xyyzzz[k] * ab_x + g_z_0_yzzzz_xxyyzzz[k];

                g_z_0_xyzzzz_xyzzzz[k] = -g_z_0_yzzzz_xyzzzz[k] * ab_x + g_z_0_yzzzz_xxyzzzz[k];

                g_z_0_xyzzzz_xzzzzz[k] = -g_z_0_yzzzz_xzzzzz[k] * ab_x + g_z_0_yzzzz_xxzzzzz[k];

                g_z_0_xyzzzz_yyyyyy[k] = -g_z_0_yzzzz_yyyyyy[k] * ab_x + g_z_0_yzzzz_xyyyyyy[k];

                g_z_0_xyzzzz_yyyyyz[k] = -g_z_0_yzzzz_yyyyyz[k] * ab_x + g_z_0_yzzzz_xyyyyyz[k];

                g_z_0_xyzzzz_yyyyzz[k] = -g_z_0_yzzzz_yyyyzz[k] * ab_x + g_z_0_yzzzz_xyyyyzz[k];

                g_z_0_xyzzzz_yyyzzz[k] = -g_z_0_yzzzz_yyyzzz[k] * ab_x + g_z_0_yzzzz_xyyyzzz[k];

                g_z_0_xyzzzz_yyzzzz[k] = -g_z_0_yzzzz_yyzzzz[k] * ab_x + g_z_0_yzzzz_xyyzzzz[k];

                g_z_0_xyzzzz_yzzzzz[k] = -g_z_0_yzzzz_yzzzzz[k] * ab_x + g_z_0_yzzzz_xyzzzzz[k];

                g_z_0_xyzzzz_zzzzzz[k] = -g_z_0_yzzzz_zzzzzz[k] * ab_x + g_z_0_yzzzz_xzzzzzz[k];
            }

            /// Set up 2128-2156 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 2128 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 2129 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 2130 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 2131 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 2132 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 2133 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 2134 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 2135 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 2136 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 2137 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 2138 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 2139 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 2140 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 2141 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 2142 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 2143 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 2144 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 2145 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 2146 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 2147 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 2148 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 2149 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 2150 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 2151 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 2152 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 2153 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 2154 * ccomps * dcomps);

            auto g_z_0_xzzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 2155 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzzzz_xxxxxx, g_z_0_xzzzzz_xxxxxy, g_z_0_xzzzzz_xxxxxz, g_z_0_xzzzzz_xxxxyy, g_z_0_xzzzzz_xxxxyz, g_z_0_xzzzzz_xxxxzz, g_z_0_xzzzzz_xxxyyy, g_z_0_xzzzzz_xxxyyz, g_z_0_xzzzzz_xxxyzz, g_z_0_xzzzzz_xxxzzz, g_z_0_xzzzzz_xxyyyy, g_z_0_xzzzzz_xxyyyz, g_z_0_xzzzzz_xxyyzz, g_z_0_xzzzzz_xxyzzz, g_z_0_xzzzzz_xxzzzz, g_z_0_xzzzzz_xyyyyy, g_z_0_xzzzzz_xyyyyz, g_z_0_xzzzzz_xyyyzz, g_z_0_xzzzzz_xyyzzz, g_z_0_xzzzzz_xyzzzz, g_z_0_xzzzzz_xzzzzz, g_z_0_xzzzzz_yyyyyy, g_z_0_xzzzzz_yyyyyz, g_z_0_xzzzzz_yyyyzz, g_z_0_xzzzzz_yyyzzz, g_z_0_xzzzzz_yyzzzz, g_z_0_xzzzzz_yzzzzz, g_z_0_xzzzzz_zzzzzz, g_z_0_zzzzz_xxxxxx, g_z_0_zzzzz_xxxxxxx, g_z_0_zzzzz_xxxxxxy, g_z_0_zzzzz_xxxxxxz, g_z_0_zzzzz_xxxxxy, g_z_0_zzzzz_xxxxxyy, g_z_0_zzzzz_xxxxxyz, g_z_0_zzzzz_xxxxxz, g_z_0_zzzzz_xxxxxzz, g_z_0_zzzzz_xxxxyy, g_z_0_zzzzz_xxxxyyy, g_z_0_zzzzz_xxxxyyz, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxyzz, g_z_0_zzzzz_xxxxzz, g_z_0_zzzzz_xxxxzzz, g_z_0_zzzzz_xxxyyy, g_z_0_zzzzz_xxxyyyy, g_z_0_zzzzz_xxxyyyz, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyyzz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxyzzz, g_z_0_zzzzz_xxxzzz, g_z_0_zzzzz_xxxzzzz, g_z_0_zzzzz_xxyyyy, g_z_0_zzzzz_xxyyyyy, g_z_0_zzzzz_xxyyyyz, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyyzz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyyzzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxyzzzz, g_z_0_zzzzz_xxzzzz, g_z_0_zzzzz_xxzzzzz, g_z_0_zzzzz_xyyyyy, g_z_0_zzzzz_xyyyyyy, g_z_0_zzzzz_xyyyyyz, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyyzz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyyzzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyyzzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xyzzzzz, g_z_0_zzzzz_xzzzzz, g_z_0_zzzzz_xzzzzzz, g_z_0_zzzzz_yyyyyy, g_z_0_zzzzz_yyyyyz, g_z_0_zzzzz_yyyyzz, g_z_0_zzzzz_yyyzzz, g_z_0_zzzzz_yyzzzz, g_z_0_zzzzz_yzzzzz, g_z_0_zzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_xxxxxx[k] = -g_z_0_zzzzz_xxxxxx[k] * ab_x + g_z_0_zzzzz_xxxxxxx[k];

                g_z_0_xzzzzz_xxxxxy[k] = -g_z_0_zzzzz_xxxxxy[k] * ab_x + g_z_0_zzzzz_xxxxxxy[k];

                g_z_0_xzzzzz_xxxxxz[k] = -g_z_0_zzzzz_xxxxxz[k] * ab_x + g_z_0_zzzzz_xxxxxxz[k];

                g_z_0_xzzzzz_xxxxyy[k] = -g_z_0_zzzzz_xxxxyy[k] * ab_x + g_z_0_zzzzz_xxxxxyy[k];

                g_z_0_xzzzzz_xxxxyz[k] = -g_z_0_zzzzz_xxxxyz[k] * ab_x + g_z_0_zzzzz_xxxxxyz[k];

                g_z_0_xzzzzz_xxxxzz[k] = -g_z_0_zzzzz_xxxxzz[k] * ab_x + g_z_0_zzzzz_xxxxxzz[k];

                g_z_0_xzzzzz_xxxyyy[k] = -g_z_0_zzzzz_xxxyyy[k] * ab_x + g_z_0_zzzzz_xxxxyyy[k];

                g_z_0_xzzzzz_xxxyyz[k] = -g_z_0_zzzzz_xxxyyz[k] * ab_x + g_z_0_zzzzz_xxxxyyz[k];

                g_z_0_xzzzzz_xxxyzz[k] = -g_z_0_zzzzz_xxxyzz[k] * ab_x + g_z_0_zzzzz_xxxxyzz[k];

                g_z_0_xzzzzz_xxxzzz[k] = -g_z_0_zzzzz_xxxzzz[k] * ab_x + g_z_0_zzzzz_xxxxzzz[k];

                g_z_0_xzzzzz_xxyyyy[k] = -g_z_0_zzzzz_xxyyyy[k] * ab_x + g_z_0_zzzzz_xxxyyyy[k];

                g_z_0_xzzzzz_xxyyyz[k] = -g_z_0_zzzzz_xxyyyz[k] * ab_x + g_z_0_zzzzz_xxxyyyz[k];

                g_z_0_xzzzzz_xxyyzz[k] = -g_z_0_zzzzz_xxyyzz[k] * ab_x + g_z_0_zzzzz_xxxyyzz[k];

                g_z_0_xzzzzz_xxyzzz[k] = -g_z_0_zzzzz_xxyzzz[k] * ab_x + g_z_0_zzzzz_xxxyzzz[k];

                g_z_0_xzzzzz_xxzzzz[k] = -g_z_0_zzzzz_xxzzzz[k] * ab_x + g_z_0_zzzzz_xxxzzzz[k];

                g_z_0_xzzzzz_xyyyyy[k] = -g_z_0_zzzzz_xyyyyy[k] * ab_x + g_z_0_zzzzz_xxyyyyy[k];

                g_z_0_xzzzzz_xyyyyz[k] = -g_z_0_zzzzz_xyyyyz[k] * ab_x + g_z_0_zzzzz_xxyyyyz[k];

                g_z_0_xzzzzz_xyyyzz[k] = -g_z_0_zzzzz_xyyyzz[k] * ab_x + g_z_0_zzzzz_xxyyyzz[k];

                g_z_0_xzzzzz_xyyzzz[k] = -g_z_0_zzzzz_xyyzzz[k] * ab_x + g_z_0_zzzzz_xxyyzzz[k];

                g_z_0_xzzzzz_xyzzzz[k] = -g_z_0_zzzzz_xyzzzz[k] * ab_x + g_z_0_zzzzz_xxyzzzz[k];

                g_z_0_xzzzzz_xzzzzz[k] = -g_z_0_zzzzz_xzzzzz[k] * ab_x + g_z_0_zzzzz_xxzzzzz[k];

                g_z_0_xzzzzz_yyyyyy[k] = -g_z_0_zzzzz_yyyyyy[k] * ab_x + g_z_0_zzzzz_xyyyyyy[k];

                g_z_0_xzzzzz_yyyyyz[k] = -g_z_0_zzzzz_yyyyyz[k] * ab_x + g_z_0_zzzzz_xyyyyyz[k];

                g_z_0_xzzzzz_yyyyzz[k] = -g_z_0_zzzzz_yyyyzz[k] * ab_x + g_z_0_zzzzz_xyyyyzz[k];

                g_z_0_xzzzzz_yyyzzz[k] = -g_z_0_zzzzz_yyyzzz[k] * ab_x + g_z_0_zzzzz_xyyyzzz[k];

                g_z_0_xzzzzz_yyzzzz[k] = -g_z_0_zzzzz_yyzzzz[k] * ab_x + g_z_0_zzzzz_xyyzzzz[k];

                g_z_0_xzzzzz_yzzzzz[k] = -g_z_0_zzzzz_yzzzzz[k] * ab_x + g_z_0_zzzzz_xyzzzzz[k];

                g_z_0_xzzzzz_zzzzzz[k] = -g_z_0_zzzzz_zzzzzz[k] * ab_x + g_z_0_zzzzz_xzzzzzz[k];
            }

            /// Set up 2156-2184 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_xxxxxx = cbuffer.data(ii_geom_10_off + 2156 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxxxy = cbuffer.data(ii_geom_10_off + 2157 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxxxz = cbuffer.data(ii_geom_10_off + 2158 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxxyy = cbuffer.data(ii_geom_10_off + 2159 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxxyz = cbuffer.data(ii_geom_10_off + 2160 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxxzz = cbuffer.data(ii_geom_10_off + 2161 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxyyy = cbuffer.data(ii_geom_10_off + 2162 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxyyz = cbuffer.data(ii_geom_10_off + 2163 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxyzz = cbuffer.data(ii_geom_10_off + 2164 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxzzz = cbuffer.data(ii_geom_10_off + 2165 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxyyyy = cbuffer.data(ii_geom_10_off + 2166 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxyyyz = cbuffer.data(ii_geom_10_off + 2167 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxyyzz = cbuffer.data(ii_geom_10_off + 2168 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxyzzz = cbuffer.data(ii_geom_10_off + 2169 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxzzzz = cbuffer.data(ii_geom_10_off + 2170 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyyyyy = cbuffer.data(ii_geom_10_off + 2171 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyyyyz = cbuffer.data(ii_geom_10_off + 2172 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyyyzz = cbuffer.data(ii_geom_10_off + 2173 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyyzzz = cbuffer.data(ii_geom_10_off + 2174 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyzzzz = cbuffer.data(ii_geom_10_off + 2175 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xzzzzz = cbuffer.data(ii_geom_10_off + 2176 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyyyyy = cbuffer.data(ii_geom_10_off + 2177 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyyyyz = cbuffer.data(ii_geom_10_off + 2178 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyyyzz = cbuffer.data(ii_geom_10_off + 2179 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyyzzz = cbuffer.data(ii_geom_10_off + 2180 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyzzzz = cbuffer.data(ii_geom_10_off + 2181 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yzzzzz = cbuffer.data(ii_geom_10_off + 2182 * ccomps * dcomps);

            auto g_z_0_yyyyyy_zzzzzz = cbuffer.data(ii_geom_10_off + 2183 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyy_xxxxxx, g_z_0_yyyyy_xxxxxxy, g_z_0_yyyyy_xxxxxy, g_z_0_yyyyy_xxxxxyy, g_z_0_yyyyy_xxxxxyz, g_z_0_yyyyy_xxxxxz, g_z_0_yyyyy_xxxxyy, g_z_0_yyyyy_xxxxyyy, g_z_0_yyyyy_xxxxyyz, g_z_0_yyyyy_xxxxyz, g_z_0_yyyyy_xxxxyzz, g_z_0_yyyyy_xxxxzz, g_z_0_yyyyy_xxxyyy, g_z_0_yyyyy_xxxyyyy, g_z_0_yyyyy_xxxyyyz, g_z_0_yyyyy_xxxyyz, g_z_0_yyyyy_xxxyyzz, g_z_0_yyyyy_xxxyzz, g_z_0_yyyyy_xxxyzzz, g_z_0_yyyyy_xxxzzz, g_z_0_yyyyy_xxyyyy, g_z_0_yyyyy_xxyyyyy, g_z_0_yyyyy_xxyyyyz, g_z_0_yyyyy_xxyyyz, g_z_0_yyyyy_xxyyyzz, g_z_0_yyyyy_xxyyzz, g_z_0_yyyyy_xxyyzzz, g_z_0_yyyyy_xxyzzz, g_z_0_yyyyy_xxyzzzz, g_z_0_yyyyy_xxzzzz, g_z_0_yyyyy_xyyyyy, g_z_0_yyyyy_xyyyyyy, g_z_0_yyyyy_xyyyyyz, g_z_0_yyyyy_xyyyyz, g_z_0_yyyyy_xyyyyzz, g_z_0_yyyyy_xyyyzz, g_z_0_yyyyy_xyyyzzz, g_z_0_yyyyy_xyyzzz, g_z_0_yyyyy_xyyzzzz, g_z_0_yyyyy_xyzzzz, g_z_0_yyyyy_xyzzzzz, g_z_0_yyyyy_xzzzzz, g_z_0_yyyyy_yyyyyy, g_z_0_yyyyy_yyyyyyy, g_z_0_yyyyy_yyyyyyz, g_z_0_yyyyy_yyyyyz, g_z_0_yyyyy_yyyyyzz, g_z_0_yyyyy_yyyyzz, g_z_0_yyyyy_yyyyzzz, g_z_0_yyyyy_yyyzzz, g_z_0_yyyyy_yyyzzzz, g_z_0_yyyyy_yyzzzz, g_z_0_yyyyy_yyzzzzz, g_z_0_yyyyy_yzzzzz, g_z_0_yyyyy_yzzzzzz, g_z_0_yyyyy_zzzzzz, g_z_0_yyyyyy_xxxxxx, g_z_0_yyyyyy_xxxxxy, g_z_0_yyyyyy_xxxxxz, g_z_0_yyyyyy_xxxxyy, g_z_0_yyyyyy_xxxxyz, g_z_0_yyyyyy_xxxxzz, g_z_0_yyyyyy_xxxyyy, g_z_0_yyyyyy_xxxyyz, g_z_0_yyyyyy_xxxyzz, g_z_0_yyyyyy_xxxzzz, g_z_0_yyyyyy_xxyyyy, g_z_0_yyyyyy_xxyyyz, g_z_0_yyyyyy_xxyyzz, g_z_0_yyyyyy_xxyzzz, g_z_0_yyyyyy_xxzzzz, g_z_0_yyyyyy_xyyyyy, g_z_0_yyyyyy_xyyyyz, g_z_0_yyyyyy_xyyyzz, g_z_0_yyyyyy_xyyzzz, g_z_0_yyyyyy_xyzzzz, g_z_0_yyyyyy_xzzzzz, g_z_0_yyyyyy_yyyyyy, g_z_0_yyyyyy_yyyyyz, g_z_0_yyyyyy_yyyyzz, g_z_0_yyyyyy_yyyzzz, g_z_0_yyyyyy_yyzzzz, g_z_0_yyyyyy_yzzzzz, g_z_0_yyyyyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_xxxxxx[k] = -g_z_0_yyyyy_xxxxxx[k] * ab_y + g_z_0_yyyyy_xxxxxxy[k];

                g_z_0_yyyyyy_xxxxxy[k] = -g_z_0_yyyyy_xxxxxy[k] * ab_y + g_z_0_yyyyy_xxxxxyy[k];

                g_z_0_yyyyyy_xxxxxz[k] = -g_z_0_yyyyy_xxxxxz[k] * ab_y + g_z_0_yyyyy_xxxxxyz[k];

                g_z_0_yyyyyy_xxxxyy[k] = -g_z_0_yyyyy_xxxxyy[k] * ab_y + g_z_0_yyyyy_xxxxyyy[k];

                g_z_0_yyyyyy_xxxxyz[k] = -g_z_0_yyyyy_xxxxyz[k] * ab_y + g_z_0_yyyyy_xxxxyyz[k];

                g_z_0_yyyyyy_xxxxzz[k] = -g_z_0_yyyyy_xxxxzz[k] * ab_y + g_z_0_yyyyy_xxxxyzz[k];

                g_z_0_yyyyyy_xxxyyy[k] = -g_z_0_yyyyy_xxxyyy[k] * ab_y + g_z_0_yyyyy_xxxyyyy[k];

                g_z_0_yyyyyy_xxxyyz[k] = -g_z_0_yyyyy_xxxyyz[k] * ab_y + g_z_0_yyyyy_xxxyyyz[k];

                g_z_0_yyyyyy_xxxyzz[k] = -g_z_0_yyyyy_xxxyzz[k] * ab_y + g_z_0_yyyyy_xxxyyzz[k];

                g_z_0_yyyyyy_xxxzzz[k] = -g_z_0_yyyyy_xxxzzz[k] * ab_y + g_z_0_yyyyy_xxxyzzz[k];

                g_z_0_yyyyyy_xxyyyy[k] = -g_z_0_yyyyy_xxyyyy[k] * ab_y + g_z_0_yyyyy_xxyyyyy[k];

                g_z_0_yyyyyy_xxyyyz[k] = -g_z_0_yyyyy_xxyyyz[k] * ab_y + g_z_0_yyyyy_xxyyyyz[k];

                g_z_0_yyyyyy_xxyyzz[k] = -g_z_0_yyyyy_xxyyzz[k] * ab_y + g_z_0_yyyyy_xxyyyzz[k];

                g_z_0_yyyyyy_xxyzzz[k] = -g_z_0_yyyyy_xxyzzz[k] * ab_y + g_z_0_yyyyy_xxyyzzz[k];

                g_z_0_yyyyyy_xxzzzz[k] = -g_z_0_yyyyy_xxzzzz[k] * ab_y + g_z_0_yyyyy_xxyzzzz[k];

                g_z_0_yyyyyy_xyyyyy[k] = -g_z_0_yyyyy_xyyyyy[k] * ab_y + g_z_0_yyyyy_xyyyyyy[k];

                g_z_0_yyyyyy_xyyyyz[k] = -g_z_0_yyyyy_xyyyyz[k] * ab_y + g_z_0_yyyyy_xyyyyyz[k];

                g_z_0_yyyyyy_xyyyzz[k] = -g_z_0_yyyyy_xyyyzz[k] * ab_y + g_z_0_yyyyy_xyyyyzz[k];

                g_z_0_yyyyyy_xyyzzz[k] = -g_z_0_yyyyy_xyyzzz[k] * ab_y + g_z_0_yyyyy_xyyyzzz[k];

                g_z_0_yyyyyy_xyzzzz[k] = -g_z_0_yyyyy_xyzzzz[k] * ab_y + g_z_0_yyyyy_xyyzzzz[k];

                g_z_0_yyyyyy_xzzzzz[k] = -g_z_0_yyyyy_xzzzzz[k] * ab_y + g_z_0_yyyyy_xyzzzzz[k];

                g_z_0_yyyyyy_yyyyyy[k] = -g_z_0_yyyyy_yyyyyy[k] * ab_y + g_z_0_yyyyy_yyyyyyy[k];

                g_z_0_yyyyyy_yyyyyz[k] = -g_z_0_yyyyy_yyyyyz[k] * ab_y + g_z_0_yyyyy_yyyyyyz[k];

                g_z_0_yyyyyy_yyyyzz[k] = -g_z_0_yyyyy_yyyyzz[k] * ab_y + g_z_0_yyyyy_yyyyyzz[k];

                g_z_0_yyyyyy_yyyzzz[k] = -g_z_0_yyyyy_yyyzzz[k] * ab_y + g_z_0_yyyyy_yyyyzzz[k];

                g_z_0_yyyyyy_yyzzzz[k] = -g_z_0_yyyyy_yyzzzz[k] * ab_y + g_z_0_yyyyy_yyyzzzz[k];

                g_z_0_yyyyyy_yzzzzz[k] = -g_z_0_yyyyy_yzzzzz[k] * ab_y + g_z_0_yyyyy_yyzzzzz[k];

                g_z_0_yyyyyy_zzzzzz[k] = -g_z_0_yyyyy_zzzzzz[k] * ab_y + g_z_0_yyyyy_yzzzzzz[k];
            }

            /// Set up 2184-2212 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_xxxxxx = cbuffer.data(ii_geom_10_off + 2184 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxxxy = cbuffer.data(ii_geom_10_off + 2185 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxxxz = cbuffer.data(ii_geom_10_off + 2186 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxxyy = cbuffer.data(ii_geom_10_off + 2187 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxxyz = cbuffer.data(ii_geom_10_off + 2188 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxxzz = cbuffer.data(ii_geom_10_off + 2189 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxyyy = cbuffer.data(ii_geom_10_off + 2190 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxyyz = cbuffer.data(ii_geom_10_off + 2191 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxyzz = cbuffer.data(ii_geom_10_off + 2192 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxzzz = cbuffer.data(ii_geom_10_off + 2193 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxyyyy = cbuffer.data(ii_geom_10_off + 2194 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxyyyz = cbuffer.data(ii_geom_10_off + 2195 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxyyzz = cbuffer.data(ii_geom_10_off + 2196 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxyzzz = cbuffer.data(ii_geom_10_off + 2197 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxzzzz = cbuffer.data(ii_geom_10_off + 2198 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyyyyy = cbuffer.data(ii_geom_10_off + 2199 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyyyyz = cbuffer.data(ii_geom_10_off + 2200 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyyyzz = cbuffer.data(ii_geom_10_off + 2201 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyyzzz = cbuffer.data(ii_geom_10_off + 2202 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyzzzz = cbuffer.data(ii_geom_10_off + 2203 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xzzzzz = cbuffer.data(ii_geom_10_off + 2204 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyyyyy = cbuffer.data(ii_geom_10_off + 2205 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyyyyz = cbuffer.data(ii_geom_10_off + 2206 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyyyzz = cbuffer.data(ii_geom_10_off + 2207 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyyzzz = cbuffer.data(ii_geom_10_off + 2208 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyzzzz = cbuffer.data(ii_geom_10_off + 2209 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yzzzzz = cbuffer.data(ii_geom_10_off + 2210 * ccomps * dcomps);

            auto g_z_0_yyyyyz_zzzzzz = cbuffer.data(ii_geom_10_off + 2211 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyyz_xxxxxx, g_z_0_yyyyyz_xxxxxy, g_z_0_yyyyyz_xxxxxz, g_z_0_yyyyyz_xxxxyy, g_z_0_yyyyyz_xxxxyz, g_z_0_yyyyyz_xxxxzz, g_z_0_yyyyyz_xxxyyy, g_z_0_yyyyyz_xxxyyz, g_z_0_yyyyyz_xxxyzz, g_z_0_yyyyyz_xxxzzz, g_z_0_yyyyyz_xxyyyy, g_z_0_yyyyyz_xxyyyz, g_z_0_yyyyyz_xxyyzz, g_z_0_yyyyyz_xxyzzz, g_z_0_yyyyyz_xxzzzz, g_z_0_yyyyyz_xyyyyy, g_z_0_yyyyyz_xyyyyz, g_z_0_yyyyyz_xyyyzz, g_z_0_yyyyyz_xyyzzz, g_z_0_yyyyyz_xyzzzz, g_z_0_yyyyyz_xzzzzz, g_z_0_yyyyyz_yyyyyy, g_z_0_yyyyyz_yyyyyz, g_z_0_yyyyyz_yyyyzz, g_z_0_yyyyyz_yyyzzz, g_z_0_yyyyyz_yyzzzz, g_z_0_yyyyyz_yzzzzz, g_z_0_yyyyyz_zzzzzz, g_z_0_yyyyz_xxxxxx, g_z_0_yyyyz_xxxxxxy, g_z_0_yyyyz_xxxxxy, g_z_0_yyyyz_xxxxxyy, g_z_0_yyyyz_xxxxxyz, g_z_0_yyyyz_xxxxxz, g_z_0_yyyyz_xxxxyy, g_z_0_yyyyz_xxxxyyy, g_z_0_yyyyz_xxxxyyz, g_z_0_yyyyz_xxxxyz, g_z_0_yyyyz_xxxxyzz, g_z_0_yyyyz_xxxxzz, g_z_0_yyyyz_xxxyyy, g_z_0_yyyyz_xxxyyyy, g_z_0_yyyyz_xxxyyyz, g_z_0_yyyyz_xxxyyz, g_z_0_yyyyz_xxxyyzz, g_z_0_yyyyz_xxxyzz, g_z_0_yyyyz_xxxyzzz, g_z_0_yyyyz_xxxzzz, g_z_0_yyyyz_xxyyyy, g_z_0_yyyyz_xxyyyyy, g_z_0_yyyyz_xxyyyyz, g_z_0_yyyyz_xxyyyz, g_z_0_yyyyz_xxyyyzz, g_z_0_yyyyz_xxyyzz, g_z_0_yyyyz_xxyyzzz, g_z_0_yyyyz_xxyzzz, g_z_0_yyyyz_xxyzzzz, g_z_0_yyyyz_xxzzzz, g_z_0_yyyyz_xyyyyy, g_z_0_yyyyz_xyyyyyy, g_z_0_yyyyz_xyyyyyz, g_z_0_yyyyz_xyyyyz, g_z_0_yyyyz_xyyyyzz, g_z_0_yyyyz_xyyyzz, g_z_0_yyyyz_xyyyzzz, g_z_0_yyyyz_xyyzzz, g_z_0_yyyyz_xyyzzzz, g_z_0_yyyyz_xyzzzz, g_z_0_yyyyz_xyzzzzz, g_z_0_yyyyz_xzzzzz, g_z_0_yyyyz_yyyyyy, g_z_0_yyyyz_yyyyyyy, g_z_0_yyyyz_yyyyyyz, g_z_0_yyyyz_yyyyyz, g_z_0_yyyyz_yyyyyzz, g_z_0_yyyyz_yyyyzz, g_z_0_yyyyz_yyyyzzz, g_z_0_yyyyz_yyyzzz, g_z_0_yyyyz_yyyzzzz, g_z_0_yyyyz_yyzzzz, g_z_0_yyyyz_yyzzzzz, g_z_0_yyyyz_yzzzzz, g_z_0_yyyyz_yzzzzzz, g_z_0_yyyyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_xxxxxx[k] = -g_z_0_yyyyz_xxxxxx[k] * ab_y + g_z_0_yyyyz_xxxxxxy[k];

                g_z_0_yyyyyz_xxxxxy[k] = -g_z_0_yyyyz_xxxxxy[k] * ab_y + g_z_0_yyyyz_xxxxxyy[k];

                g_z_0_yyyyyz_xxxxxz[k] = -g_z_0_yyyyz_xxxxxz[k] * ab_y + g_z_0_yyyyz_xxxxxyz[k];

                g_z_0_yyyyyz_xxxxyy[k] = -g_z_0_yyyyz_xxxxyy[k] * ab_y + g_z_0_yyyyz_xxxxyyy[k];

                g_z_0_yyyyyz_xxxxyz[k] = -g_z_0_yyyyz_xxxxyz[k] * ab_y + g_z_0_yyyyz_xxxxyyz[k];

                g_z_0_yyyyyz_xxxxzz[k] = -g_z_0_yyyyz_xxxxzz[k] * ab_y + g_z_0_yyyyz_xxxxyzz[k];

                g_z_0_yyyyyz_xxxyyy[k] = -g_z_0_yyyyz_xxxyyy[k] * ab_y + g_z_0_yyyyz_xxxyyyy[k];

                g_z_0_yyyyyz_xxxyyz[k] = -g_z_0_yyyyz_xxxyyz[k] * ab_y + g_z_0_yyyyz_xxxyyyz[k];

                g_z_0_yyyyyz_xxxyzz[k] = -g_z_0_yyyyz_xxxyzz[k] * ab_y + g_z_0_yyyyz_xxxyyzz[k];

                g_z_0_yyyyyz_xxxzzz[k] = -g_z_0_yyyyz_xxxzzz[k] * ab_y + g_z_0_yyyyz_xxxyzzz[k];

                g_z_0_yyyyyz_xxyyyy[k] = -g_z_0_yyyyz_xxyyyy[k] * ab_y + g_z_0_yyyyz_xxyyyyy[k];

                g_z_0_yyyyyz_xxyyyz[k] = -g_z_0_yyyyz_xxyyyz[k] * ab_y + g_z_0_yyyyz_xxyyyyz[k];

                g_z_0_yyyyyz_xxyyzz[k] = -g_z_0_yyyyz_xxyyzz[k] * ab_y + g_z_0_yyyyz_xxyyyzz[k];

                g_z_0_yyyyyz_xxyzzz[k] = -g_z_0_yyyyz_xxyzzz[k] * ab_y + g_z_0_yyyyz_xxyyzzz[k];

                g_z_0_yyyyyz_xxzzzz[k] = -g_z_0_yyyyz_xxzzzz[k] * ab_y + g_z_0_yyyyz_xxyzzzz[k];

                g_z_0_yyyyyz_xyyyyy[k] = -g_z_0_yyyyz_xyyyyy[k] * ab_y + g_z_0_yyyyz_xyyyyyy[k];

                g_z_0_yyyyyz_xyyyyz[k] = -g_z_0_yyyyz_xyyyyz[k] * ab_y + g_z_0_yyyyz_xyyyyyz[k];

                g_z_0_yyyyyz_xyyyzz[k] = -g_z_0_yyyyz_xyyyzz[k] * ab_y + g_z_0_yyyyz_xyyyyzz[k];

                g_z_0_yyyyyz_xyyzzz[k] = -g_z_0_yyyyz_xyyzzz[k] * ab_y + g_z_0_yyyyz_xyyyzzz[k];

                g_z_0_yyyyyz_xyzzzz[k] = -g_z_0_yyyyz_xyzzzz[k] * ab_y + g_z_0_yyyyz_xyyzzzz[k];

                g_z_0_yyyyyz_xzzzzz[k] = -g_z_0_yyyyz_xzzzzz[k] * ab_y + g_z_0_yyyyz_xyzzzzz[k];

                g_z_0_yyyyyz_yyyyyy[k] = -g_z_0_yyyyz_yyyyyy[k] * ab_y + g_z_0_yyyyz_yyyyyyy[k];

                g_z_0_yyyyyz_yyyyyz[k] = -g_z_0_yyyyz_yyyyyz[k] * ab_y + g_z_0_yyyyz_yyyyyyz[k];

                g_z_0_yyyyyz_yyyyzz[k] = -g_z_0_yyyyz_yyyyzz[k] * ab_y + g_z_0_yyyyz_yyyyyzz[k];

                g_z_0_yyyyyz_yyyzzz[k] = -g_z_0_yyyyz_yyyzzz[k] * ab_y + g_z_0_yyyyz_yyyyzzz[k];

                g_z_0_yyyyyz_yyzzzz[k] = -g_z_0_yyyyz_yyzzzz[k] * ab_y + g_z_0_yyyyz_yyyzzzz[k];

                g_z_0_yyyyyz_yzzzzz[k] = -g_z_0_yyyyz_yzzzzz[k] * ab_y + g_z_0_yyyyz_yyzzzzz[k];

                g_z_0_yyyyyz_zzzzzz[k] = -g_z_0_yyyyz_zzzzzz[k] * ab_y + g_z_0_yyyyz_yzzzzzz[k];
            }

            /// Set up 2212-2240 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_xxxxxx = cbuffer.data(ii_geom_10_off + 2212 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxxxy = cbuffer.data(ii_geom_10_off + 2213 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxxxz = cbuffer.data(ii_geom_10_off + 2214 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxxyy = cbuffer.data(ii_geom_10_off + 2215 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxxyz = cbuffer.data(ii_geom_10_off + 2216 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxxzz = cbuffer.data(ii_geom_10_off + 2217 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxyyy = cbuffer.data(ii_geom_10_off + 2218 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxyyz = cbuffer.data(ii_geom_10_off + 2219 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxyzz = cbuffer.data(ii_geom_10_off + 2220 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxzzz = cbuffer.data(ii_geom_10_off + 2221 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxyyyy = cbuffer.data(ii_geom_10_off + 2222 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxyyyz = cbuffer.data(ii_geom_10_off + 2223 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxyyzz = cbuffer.data(ii_geom_10_off + 2224 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxyzzz = cbuffer.data(ii_geom_10_off + 2225 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxzzzz = cbuffer.data(ii_geom_10_off + 2226 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyyyyy = cbuffer.data(ii_geom_10_off + 2227 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyyyyz = cbuffer.data(ii_geom_10_off + 2228 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyyyzz = cbuffer.data(ii_geom_10_off + 2229 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyyzzz = cbuffer.data(ii_geom_10_off + 2230 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyzzzz = cbuffer.data(ii_geom_10_off + 2231 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xzzzzz = cbuffer.data(ii_geom_10_off + 2232 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyyyyy = cbuffer.data(ii_geom_10_off + 2233 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyyyyz = cbuffer.data(ii_geom_10_off + 2234 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyyyzz = cbuffer.data(ii_geom_10_off + 2235 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyyzzz = cbuffer.data(ii_geom_10_off + 2236 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyzzzz = cbuffer.data(ii_geom_10_off + 2237 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yzzzzz = cbuffer.data(ii_geom_10_off + 2238 * ccomps * dcomps);

            auto g_z_0_yyyyzz_zzzzzz = cbuffer.data(ii_geom_10_off + 2239 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyzz_xxxxxx, g_z_0_yyyyzz_xxxxxy, g_z_0_yyyyzz_xxxxxz, g_z_0_yyyyzz_xxxxyy, g_z_0_yyyyzz_xxxxyz, g_z_0_yyyyzz_xxxxzz, g_z_0_yyyyzz_xxxyyy, g_z_0_yyyyzz_xxxyyz, g_z_0_yyyyzz_xxxyzz, g_z_0_yyyyzz_xxxzzz, g_z_0_yyyyzz_xxyyyy, g_z_0_yyyyzz_xxyyyz, g_z_0_yyyyzz_xxyyzz, g_z_0_yyyyzz_xxyzzz, g_z_0_yyyyzz_xxzzzz, g_z_0_yyyyzz_xyyyyy, g_z_0_yyyyzz_xyyyyz, g_z_0_yyyyzz_xyyyzz, g_z_0_yyyyzz_xyyzzz, g_z_0_yyyyzz_xyzzzz, g_z_0_yyyyzz_xzzzzz, g_z_0_yyyyzz_yyyyyy, g_z_0_yyyyzz_yyyyyz, g_z_0_yyyyzz_yyyyzz, g_z_0_yyyyzz_yyyzzz, g_z_0_yyyyzz_yyzzzz, g_z_0_yyyyzz_yzzzzz, g_z_0_yyyyzz_zzzzzz, g_z_0_yyyzz_xxxxxx, g_z_0_yyyzz_xxxxxxy, g_z_0_yyyzz_xxxxxy, g_z_0_yyyzz_xxxxxyy, g_z_0_yyyzz_xxxxxyz, g_z_0_yyyzz_xxxxxz, g_z_0_yyyzz_xxxxyy, g_z_0_yyyzz_xxxxyyy, g_z_0_yyyzz_xxxxyyz, g_z_0_yyyzz_xxxxyz, g_z_0_yyyzz_xxxxyzz, g_z_0_yyyzz_xxxxzz, g_z_0_yyyzz_xxxyyy, g_z_0_yyyzz_xxxyyyy, g_z_0_yyyzz_xxxyyyz, g_z_0_yyyzz_xxxyyz, g_z_0_yyyzz_xxxyyzz, g_z_0_yyyzz_xxxyzz, g_z_0_yyyzz_xxxyzzz, g_z_0_yyyzz_xxxzzz, g_z_0_yyyzz_xxyyyy, g_z_0_yyyzz_xxyyyyy, g_z_0_yyyzz_xxyyyyz, g_z_0_yyyzz_xxyyyz, g_z_0_yyyzz_xxyyyzz, g_z_0_yyyzz_xxyyzz, g_z_0_yyyzz_xxyyzzz, g_z_0_yyyzz_xxyzzz, g_z_0_yyyzz_xxyzzzz, g_z_0_yyyzz_xxzzzz, g_z_0_yyyzz_xyyyyy, g_z_0_yyyzz_xyyyyyy, g_z_0_yyyzz_xyyyyyz, g_z_0_yyyzz_xyyyyz, g_z_0_yyyzz_xyyyyzz, g_z_0_yyyzz_xyyyzz, g_z_0_yyyzz_xyyyzzz, g_z_0_yyyzz_xyyzzz, g_z_0_yyyzz_xyyzzzz, g_z_0_yyyzz_xyzzzz, g_z_0_yyyzz_xyzzzzz, g_z_0_yyyzz_xzzzzz, g_z_0_yyyzz_yyyyyy, g_z_0_yyyzz_yyyyyyy, g_z_0_yyyzz_yyyyyyz, g_z_0_yyyzz_yyyyyz, g_z_0_yyyzz_yyyyyzz, g_z_0_yyyzz_yyyyzz, g_z_0_yyyzz_yyyyzzz, g_z_0_yyyzz_yyyzzz, g_z_0_yyyzz_yyyzzzz, g_z_0_yyyzz_yyzzzz, g_z_0_yyyzz_yyzzzzz, g_z_0_yyyzz_yzzzzz, g_z_0_yyyzz_yzzzzzz, g_z_0_yyyzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_xxxxxx[k] = -g_z_0_yyyzz_xxxxxx[k] * ab_y + g_z_0_yyyzz_xxxxxxy[k];

                g_z_0_yyyyzz_xxxxxy[k] = -g_z_0_yyyzz_xxxxxy[k] * ab_y + g_z_0_yyyzz_xxxxxyy[k];

                g_z_0_yyyyzz_xxxxxz[k] = -g_z_0_yyyzz_xxxxxz[k] * ab_y + g_z_0_yyyzz_xxxxxyz[k];

                g_z_0_yyyyzz_xxxxyy[k] = -g_z_0_yyyzz_xxxxyy[k] * ab_y + g_z_0_yyyzz_xxxxyyy[k];

                g_z_0_yyyyzz_xxxxyz[k] = -g_z_0_yyyzz_xxxxyz[k] * ab_y + g_z_0_yyyzz_xxxxyyz[k];

                g_z_0_yyyyzz_xxxxzz[k] = -g_z_0_yyyzz_xxxxzz[k] * ab_y + g_z_0_yyyzz_xxxxyzz[k];

                g_z_0_yyyyzz_xxxyyy[k] = -g_z_0_yyyzz_xxxyyy[k] * ab_y + g_z_0_yyyzz_xxxyyyy[k];

                g_z_0_yyyyzz_xxxyyz[k] = -g_z_0_yyyzz_xxxyyz[k] * ab_y + g_z_0_yyyzz_xxxyyyz[k];

                g_z_0_yyyyzz_xxxyzz[k] = -g_z_0_yyyzz_xxxyzz[k] * ab_y + g_z_0_yyyzz_xxxyyzz[k];

                g_z_0_yyyyzz_xxxzzz[k] = -g_z_0_yyyzz_xxxzzz[k] * ab_y + g_z_0_yyyzz_xxxyzzz[k];

                g_z_0_yyyyzz_xxyyyy[k] = -g_z_0_yyyzz_xxyyyy[k] * ab_y + g_z_0_yyyzz_xxyyyyy[k];

                g_z_0_yyyyzz_xxyyyz[k] = -g_z_0_yyyzz_xxyyyz[k] * ab_y + g_z_0_yyyzz_xxyyyyz[k];

                g_z_0_yyyyzz_xxyyzz[k] = -g_z_0_yyyzz_xxyyzz[k] * ab_y + g_z_0_yyyzz_xxyyyzz[k];

                g_z_0_yyyyzz_xxyzzz[k] = -g_z_0_yyyzz_xxyzzz[k] * ab_y + g_z_0_yyyzz_xxyyzzz[k];

                g_z_0_yyyyzz_xxzzzz[k] = -g_z_0_yyyzz_xxzzzz[k] * ab_y + g_z_0_yyyzz_xxyzzzz[k];

                g_z_0_yyyyzz_xyyyyy[k] = -g_z_0_yyyzz_xyyyyy[k] * ab_y + g_z_0_yyyzz_xyyyyyy[k];

                g_z_0_yyyyzz_xyyyyz[k] = -g_z_0_yyyzz_xyyyyz[k] * ab_y + g_z_0_yyyzz_xyyyyyz[k];

                g_z_0_yyyyzz_xyyyzz[k] = -g_z_0_yyyzz_xyyyzz[k] * ab_y + g_z_0_yyyzz_xyyyyzz[k];

                g_z_0_yyyyzz_xyyzzz[k] = -g_z_0_yyyzz_xyyzzz[k] * ab_y + g_z_0_yyyzz_xyyyzzz[k];

                g_z_0_yyyyzz_xyzzzz[k] = -g_z_0_yyyzz_xyzzzz[k] * ab_y + g_z_0_yyyzz_xyyzzzz[k];

                g_z_0_yyyyzz_xzzzzz[k] = -g_z_0_yyyzz_xzzzzz[k] * ab_y + g_z_0_yyyzz_xyzzzzz[k];

                g_z_0_yyyyzz_yyyyyy[k] = -g_z_0_yyyzz_yyyyyy[k] * ab_y + g_z_0_yyyzz_yyyyyyy[k];

                g_z_0_yyyyzz_yyyyyz[k] = -g_z_0_yyyzz_yyyyyz[k] * ab_y + g_z_0_yyyzz_yyyyyyz[k];

                g_z_0_yyyyzz_yyyyzz[k] = -g_z_0_yyyzz_yyyyzz[k] * ab_y + g_z_0_yyyzz_yyyyyzz[k];

                g_z_0_yyyyzz_yyyzzz[k] = -g_z_0_yyyzz_yyyzzz[k] * ab_y + g_z_0_yyyzz_yyyyzzz[k];

                g_z_0_yyyyzz_yyzzzz[k] = -g_z_0_yyyzz_yyzzzz[k] * ab_y + g_z_0_yyyzz_yyyzzzz[k];

                g_z_0_yyyyzz_yzzzzz[k] = -g_z_0_yyyzz_yzzzzz[k] * ab_y + g_z_0_yyyzz_yyzzzzz[k];

                g_z_0_yyyyzz_zzzzzz[k] = -g_z_0_yyyzz_zzzzzz[k] * ab_y + g_z_0_yyyzz_yzzzzzz[k];
            }

            /// Set up 2240-2268 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 2240 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 2241 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 2242 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 2243 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 2244 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 2245 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 2246 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 2247 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 2248 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 2249 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 2250 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 2251 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 2252 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 2253 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 2254 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 2255 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 2256 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 2257 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 2258 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 2259 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 2260 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 2261 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 2262 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 2263 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 2264 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 2265 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 2266 * ccomps * dcomps);

            auto g_z_0_yyyzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 2267 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyzzz_xxxxxx, g_z_0_yyyzzz_xxxxxy, g_z_0_yyyzzz_xxxxxz, g_z_0_yyyzzz_xxxxyy, g_z_0_yyyzzz_xxxxyz, g_z_0_yyyzzz_xxxxzz, g_z_0_yyyzzz_xxxyyy, g_z_0_yyyzzz_xxxyyz, g_z_0_yyyzzz_xxxyzz, g_z_0_yyyzzz_xxxzzz, g_z_0_yyyzzz_xxyyyy, g_z_0_yyyzzz_xxyyyz, g_z_0_yyyzzz_xxyyzz, g_z_0_yyyzzz_xxyzzz, g_z_0_yyyzzz_xxzzzz, g_z_0_yyyzzz_xyyyyy, g_z_0_yyyzzz_xyyyyz, g_z_0_yyyzzz_xyyyzz, g_z_0_yyyzzz_xyyzzz, g_z_0_yyyzzz_xyzzzz, g_z_0_yyyzzz_xzzzzz, g_z_0_yyyzzz_yyyyyy, g_z_0_yyyzzz_yyyyyz, g_z_0_yyyzzz_yyyyzz, g_z_0_yyyzzz_yyyzzz, g_z_0_yyyzzz_yyzzzz, g_z_0_yyyzzz_yzzzzz, g_z_0_yyyzzz_zzzzzz, g_z_0_yyzzz_xxxxxx, g_z_0_yyzzz_xxxxxxy, g_z_0_yyzzz_xxxxxy, g_z_0_yyzzz_xxxxxyy, g_z_0_yyzzz_xxxxxyz, g_z_0_yyzzz_xxxxxz, g_z_0_yyzzz_xxxxyy, g_z_0_yyzzz_xxxxyyy, g_z_0_yyzzz_xxxxyyz, g_z_0_yyzzz_xxxxyz, g_z_0_yyzzz_xxxxyzz, g_z_0_yyzzz_xxxxzz, g_z_0_yyzzz_xxxyyy, g_z_0_yyzzz_xxxyyyy, g_z_0_yyzzz_xxxyyyz, g_z_0_yyzzz_xxxyyz, g_z_0_yyzzz_xxxyyzz, g_z_0_yyzzz_xxxyzz, g_z_0_yyzzz_xxxyzzz, g_z_0_yyzzz_xxxzzz, g_z_0_yyzzz_xxyyyy, g_z_0_yyzzz_xxyyyyy, g_z_0_yyzzz_xxyyyyz, g_z_0_yyzzz_xxyyyz, g_z_0_yyzzz_xxyyyzz, g_z_0_yyzzz_xxyyzz, g_z_0_yyzzz_xxyyzzz, g_z_0_yyzzz_xxyzzz, g_z_0_yyzzz_xxyzzzz, g_z_0_yyzzz_xxzzzz, g_z_0_yyzzz_xyyyyy, g_z_0_yyzzz_xyyyyyy, g_z_0_yyzzz_xyyyyyz, g_z_0_yyzzz_xyyyyz, g_z_0_yyzzz_xyyyyzz, g_z_0_yyzzz_xyyyzz, g_z_0_yyzzz_xyyyzzz, g_z_0_yyzzz_xyyzzz, g_z_0_yyzzz_xyyzzzz, g_z_0_yyzzz_xyzzzz, g_z_0_yyzzz_xyzzzzz, g_z_0_yyzzz_xzzzzz, g_z_0_yyzzz_yyyyyy, g_z_0_yyzzz_yyyyyyy, g_z_0_yyzzz_yyyyyyz, g_z_0_yyzzz_yyyyyz, g_z_0_yyzzz_yyyyyzz, g_z_0_yyzzz_yyyyzz, g_z_0_yyzzz_yyyyzzz, g_z_0_yyzzz_yyyzzz, g_z_0_yyzzz_yyyzzzz, g_z_0_yyzzz_yyzzzz, g_z_0_yyzzz_yyzzzzz, g_z_0_yyzzz_yzzzzz, g_z_0_yyzzz_yzzzzzz, g_z_0_yyzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_xxxxxx[k] = -g_z_0_yyzzz_xxxxxx[k] * ab_y + g_z_0_yyzzz_xxxxxxy[k];

                g_z_0_yyyzzz_xxxxxy[k] = -g_z_0_yyzzz_xxxxxy[k] * ab_y + g_z_0_yyzzz_xxxxxyy[k];

                g_z_0_yyyzzz_xxxxxz[k] = -g_z_0_yyzzz_xxxxxz[k] * ab_y + g_z_0_yyzzz_xxxxxyz[k];

                g_z_0_yyyzzz_xxxxyy[k] = -g_z_0_yyzzz_xxxxyy[k] * ab_y + g_z_0_yyzzz_xxxxyyy[k];

                g_z_0_yyyzzz_xxxxyz[k] = -g_z_0_yyzzz_xxxxyz[k] * ab_y + g_z_0_yyzzz_xxxxyyz[k];

                g_z_0_yyyzzz_xxxxzz[k] = -g_z_0_yyzzz_xxxxzz[k] * ab_y + g_z_0_yyzzz_xxxxyzz[k];

                g_z_0_yyyzzz_xxxyyy[k] = -g_z_0_yyzzz_xxxyyy[k] * ab_y + g_z_0_yyzzz_xxxyyyy[k];

                g_z_0_yyyzzz_xxxyyz[k] = -g_z_0_yyzzz_xxxyyz[k] * ab_y + g_z_0_yyzzz_xxxyyyz[k];

                g_z_0_yyyzzz_xxxyzz[k] = -g_z_0_yyzzz_xxxyzz[k] * ab_y + g_z_0_yyzzz_xxxyyzz[k];

                g_z_0_yyyzzz_xxxzzz[k] = -g_z_0_yyzzz_xxxzzz[k] * ab_y + g_z_0_yyzzz_xxxyzzz[k];

                g_z_0_yyyzzz_xxyyyy[k] = -g_z_0_yyzzz_xxyyyy[k] * ab_y + g_z_0_yyzzz_xxyyyyy[k];

                g_z_0_yyyzzz_xxyyyz[k] = -g_z_0_yyzzz_xxyyyz[k] * ab_y + g_z_0_yyzzz_xxyyyyz[k];

                g_z_0_yyyzzz_xxyyzz[k] = -g_z_0_yyzzz_xxyyzz[k] * ab_y + g_z_0_yyzzz_xxyyyzz[k];

                g_z_0_yyyzzz_xxyzzz[k] = -g_z_0_yyzzz_xxyzzz[k] * ab_y + g_z_0_yyzzz_xxyyzzz[k];

                g_z_0_yyyzzz_xxzzzz[k] = -g_z_0_yyzzz_xxzzzz[k] * ab_y + g_z_0_yyzzz_xxyzzzz[k];

                g_z_0_yyyzzz_xyyyyy[k] = -g_z_0_yyzzz_xyyyyy[k] * ab_y + g_z_0_yyzzz_xyyyyyy[k];

                g_z_0_yyyzzz_xyyyyz[k] = -g_z_0_yyzzz_xyyyyz[k] * ab_y + g_z_0_yyzzz_xyyyyyz[k];

                g_z_0_yyyzzz_xyyyzz[k] = -g_z_0_yyzzz_xyyyzz[k] * ab_y + g_z_0_yyzzz_xyyyyzz[k];

                g_z_0_yyyzzz_xyyzzz[k] = -g_z_0_yyzzz_xyyzzz[k] * ab_y + g_z_0_yyzzz_xyyyzzz[k];

                g_z_0_yyyzzz_xyzzzz[k] = -g_z_0_yyzzz_xyzzzz[k] * ab_y + g_z_0_yyzzz_xyyzzzz[k];

                g_z_0_yyyzzz_xzzzzz[k] = -g_z_0_yyzzz_xzzzzz[k] * ab_y + g_z_0_yyzzz_xyzzzzz[k];

                g_z_0_yyyzzz_yyyyyy[k] = -g_z_0_yyzzz_yyyyyy[k] * ab_y + g_z_0_yyzzz_yyyyyyy[k];

                g_z_0_yyyzzz_yyyyyz[k] = -g_z_0_yyzzz_yyyyyz[k] * ab_y + g_z_0_yyzzz_yyyyyyz[k];

                g_z_0_yyyzzz_yyyyzz[k] = -g_z_0_yyzzz_yyyyzz[k] * ab_y + g_z_0_yyzzz_yyyyyzz[k];

                g_z_0_yyyzzz_yyyzzz[k] = -g_z_0_yyzzz_yyyzzz[k] * ab_y + g_z_0_yyzzz_yyyyzzz[k];

                g_z_0_yyyzzz_yyzzzz[k] = -g_z_0_yyzzz_yyzzzz[k] * ab_y + g_z_0_yyzzz_yyyzzzz[k];

                g_z_0_yyyzzz_yzzzzz[k] = -g_z_0_yyzzz_yzzzzz[k] * ab_y + g_z_0_yyzzz_yyzzzzz[k];

                g_z_0_yyyzzz_zzzzzz[k] = -g_z_0_yyzzz_zzzzzz[k] * ab_y + g_z_0_yyzzz_yzzzzzz[k];
            }

            /// Set up 2268-2296 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 2268 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 2269 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 2270 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 2271 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 2272 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 2273 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 2274 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 2275 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 2276 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 2277 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 2278 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 2279 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 2280 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 2281 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 2282 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 2283 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 2284 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 2285 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 2286 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 2287 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 2288 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 2289 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 2290 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 2291 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 2292 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 2293 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 2294 * ccomps * dcomps);

            auto g_z_0_yyzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 2295 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzzzz_xxxxxx, g_z_0_yyzzzz_xxxxxy, g_z_0_yyzzzz_xxxxxz, g_z_0_yyzzzz_xxxxyy, g_z_0_yyzzzz_xxxxyz, g_z_0_yyzzzz_xxxxzz, g_z_0_yyzzzz_xxxyyy, g_z_0_yyzzzz_xxxyyz, g_z_0_yyzzzz_xxxyzz, g_z_0_yyzzzz_xxxzzz, g_z_0_yyzzzz_xxyyyy, g_z_0_yyzzzz_xxyyyz, g_z_0_yyzzzz_xxyyzz, g_z_0_yyzzzz_xxyzzz, g_z_0_yyzzzz_xxzzzz, g_z_0_yyzzzz_xyyyyy, g_z_0_yyzzzz_xyyyyz, g_z_0_yyzzzz_xyyyzz, g_z_0_yyzzzz_xyyzzz, g_z_0_yyzzzz_xyzzzz, g_z_0_yyzzzz_xzzzzz, g_z_0_yyzzzz_yyyyyy, g_z_0_yyzzzz_yyyyyz, g_z_0_yyzzzz_yyyyzz, g_z_0_yyzzzz_yyyzzz, g_z_0_yyzzzz_yyzzzz, g_z_0_yyzzzz_yzzzzz, g_z_0_yyzzzz_zzzzzz, g_z_0_yzzzz_xxxxxx, g_z_0_yzzzz_xxxxxxy, g_z_0_yzzzz_xxxxxy, g_z_0_yzzzz_xxxxxyy, g_z_0_yzzzz_xxxxxyz, g_z_0_yzzzz_xxxxxz, g_z_0_yzzzz_xxxxyy, g_z_0_yzzzz_xxxxyyy, g_z_0_yzzzz_xxxxyyz, g_z_0_yzzzz_xxxxyz, g_z_0_yzzzz_xxxxyzz, g_z_0_yzzzz_xxxxzz, g_z_0_yzzzz_xxxyyy, g_z_0_yzzzz_xxxyyyy, g_z_0_yzzzz_xxxyyyz, g_z_0_yzzzz_xxxyyz, g_z_0_yzzzz_xxxyyzz, g_z_0_yzzzz_xxxyzz, g_z_0_yzzzz_xxxyzzz, g_z_0_yzzzz_xxxzzz, g_z_0_yzzzz_xxyyyy, g_z_0_yzzzz_xxyyyyy, g_z_0_yzzzz_xxyyyyz, g_z_0_yzzzz_xxyyyz, g_z_0_yzzzz_xxyyyzz, g_z_0_yzzzz_xxyyzz, g_z_0_yzzzz_xxyyzzz, g_z_0_yzzzz_xxyzzz, g_z_0_yzzzz_xxyzzzz, g_z_0_yzzzz_xxzzzz, g_z_0_yzzzz_xyyyyy, g_z_0_yzzzz_xyyyyyy, g_z_0_yzzzz_xyyyyyz, g_z_0_yzzzz_xyyyyz, g_z_0_yzzzz_xyyyyzz, g_z_0_yzzzz_xyyyzz, g_z_0_yzzzz_xyyyzzz, g_z_0_yzzzz_xyyzzz, g_z_0_yzzzz_xyyzzzz, g_z_0_yzzzz_xyzzzz, g_z_0_yzzzz_xyzzzzz, g_z_0_yzzzz_xzzzzz, g_z_0_yzzzz_yyyyyy, g_z_0_yzzzz_yyyyyyy, g_z_0_yzzzz_yyyyyyz, g_z_0_yzzzz_yyyyyz, g_z_0_yzzzz_yyyyyzz, g_z_0_yzzzz_yyyyzz, g_z_0_yzzzz_yyyyzzz, g_z_0_yzzzz_yyyzzz, g_z_0_yzzzz_yyyzzzz, g_z_0_yzzzz_yyzzzz, g_z_0_yzzzz_yyzzzzz, g_z_0_yzzzz_yzzzzz, g_z_0_yzzzz_yzzzzzz, g_z_0_yzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_xxxxxx[k] = -g_z_0_yzzzz_xxxxxx[k] * ab_y + g_z_0_yzzzz_xxxxxxy[k];

                g_z_0_yyzzzz_xxxxxy[k] = -g_z_0_yzzzz_xxxxxy[k] * ab_y + g_z_0_yzzzz_xxxxxyy[k];

                g_z_0_yyzzzz_xxxxxz[k] = -g_z_0_yzzzz_xxxxxz[k] * ab_y + g_z_0_yzzzz_xxxxxyz[k];

                g_z_0_yyzzzz_xxxxyy[k] = -g_z_0_yzzzz_xxxxyy[k] * ab_y + g_z_0_yzzzz_xxxxyyy[k];

                g_z_0_yyzzzz_xxxxyz[k] = -g_z_0_yzzzz_xxxxyz[k] * ab_y + g_z_0_yzzzz_xxxxyyz[k];

                g_z_0_yyzzzz_xxxxzz[k] = -g_z_0_yzzzz_xxxxzz[k] * ab_y + g_z_0_yzzzz_xxxxyzz[k];

                g_z_0_yyzzzz_xxxyyy[k] = -g_z_0_yzzzz_xxxyyy[k] * ab_y + g_z_0_yzzzz_xxxyyyy[k];

                g_z_0_yyzzzz_xxxyyz[k] = -g_z_0_yzzzz_xxxyyz[k] * ab_y + g_z_0_yzzzz_xxxyyyz[k];

                g_z_0_yyzzzz_xxxyzz[k] = -g_z_0_yzzzz_xxxyzz[k] * ab_y + g_z_0_yzzzz_xxxyyzz[k];

                g_z_0_yyzzzz_xxxzzz[k] = -g_z_0_yzzzz_xxxzzz[k] * ab_y + g_z_0_yzzzz_xxxyzzz[k];

                g_z_0_yyzzzz_xxyyyy[k] = -g_z_0_yzzzz_xxyyyy[k] * ab_y + g_z_0_yzzzz_xxyyyyy[k];

                g_z_0_yyzzzz_xxyyyz[k] = -g_z_0_yzzzz_xxyyyz[k] * ab_y + g_z_0_yzzzz_xxyyyyz[k];

                g_z_0_yyzzzz_xxyyzz[k] = -g_z_0_yzzzz_xxyyzz[k] * ab_y + g_z_0_yzzzz_xxyyyzz[k];

                g_z_0_yyzzzz_xxyzzz[k] = -g_z_0_yzzzz_xxyzzz[k] * ab_y + g_z_0_yzzzz_xxyyzzz[k];

                g_z_0_yyzzzz_xxzzzz[k] = -g_z_0_yzzzz_xxzzzz[k] * ab_y + g_z_0_yzzzz_xxyzzzz[k];

                g_z_0_yyzzzz_xyyyyy[k] = -g_z_0_yzzzz_xyyyyy[k] * ab_y + g_z_0_yzzzz_xyyyyyy[k];

                g_z_0_yyzzzz_xyyyyz[k] = -g_z_0_yzzzz_xyyyyz[k] * ab_y + g_z_0_yzzzz_xyyyyyz[k];

                g_z_0_yyzzzz_xyyyzz[k] = -g_z_0_yzzzz_xyyyzz[k] * ab_y + g_z_0_yzzzz_xyyyyzz[k];

                g_z_0_yyzzzz_xyyzzz[k] = -g_z_0_yzzzz_xyyzzz[k] * ab_y + g_z_0_yzzzz_xyyyzzz[k];

                g_z_0_yyzzzz_xyzzzz[k] = -g_z_0_yzzzz_xyzzzz[k] * ab_y + g_z_0_yzzzz_xyyzzzz[k];

                g_z_0_yyzzzz_xzzzzz[k] = -g_z_0_yzzzz_xzzzzz[k] * ab_y + g_z_0_yzzzz_xyzzzzz[k];

                g_z_0_yyzzzz_yyyyyy[k] = -g_z_0_yzzzz_yyyyyy[k] * ab_y + g_z_0_yzzzz_yyyyyyy[k];

                g_z_0_yyzzzz_yyyyyz[k] = -g_z_0_yzzzz_yyyyyz[k] * ab_y + g_z_0_yzzzz_yyyyyyz[k];

                g_z_0_yyzzzz_yyyyzz[k] = -g_z_0_yzzzz_yyyyzz[k] * ab_y + g_z_0_yzzzz_yyyyyzz[k];

                g_z_0_yyzzzz_yyyzzz[k] = -g_z_0_yzzzz_yyyzzz[k] * ab_y + g_z_0_yzzzz_yyyyzzz[k];

                g_z_0_yyzzzz_yyzzzz[k] = -g_z_0_yzzzz_yyzzzz[k] * ab_y + g_z_0_yzzzz_yyyzzzz[k];

                g_z_0_yyzzzz_yzzzzz[k] = -g_z_0_yzzzz_yzzzzz[k] * ab_y + g_z_0_yzzzz_yyzzzzz[k];

                g_z_0_yyzzzz_zzzzzz[k] = -g_z_0_yzzzz_zzzzzz[k] * ab_y + g_z_0_yzzzz_yzzzzzz[k];
            }

            /// Set up 2296-2324 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 2296 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 2297 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 2298 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 2299 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 2300 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 2301 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 2302 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 2303 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 2304 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 2305 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 2306 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 2307 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 2308 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 2309 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 2310 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 2311 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 2312 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 2313 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 2314 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 2315 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 2316 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 2317 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 2318 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 2319 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 2320 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 2321 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 2322 * ccomps * dcomps);

            auto g_z_0_yzzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 2323 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzzzz_xxxxxx, g_z_0_yzzzzz_xxxxxy, g_z_0_yzzzzz_xxxxxz, g_z_0_yzzzzz_xxxxyy, g_z_0_yzzzzz_xxxxyz, g_z_0_yzzzzz_xxxxzz, g_z_0_yzzzzz_xxxyyy, g_z_0_yzzzzz_xxxyyz, g_z_0_yzzzzz_xxxyzz, g_z_0_yzzzzz_xxxzzz, g_z_0_yzzzzz_xxyyyy, g_z_0_yzzzzz_xxyyyz, g_z_0_yzzzzz_xxyyzz, g_z_0_yzzzzz_xxyzzz, g_z_0_yzzzzz_xxzzzz, g_z_0_yzzzzz_xyyyyy, g_z_0_yzzzzz_xyyyyz, g_z_0_yzzzzz_xyyyzz, g_z_0_yzzzzz_xyyzzz, g_z_0_yzzzzz_xyzzzz, g_z_0_yzzzzz_xzzzzz, g_z_0_yzzzzz_yyyyyy, g_z_0_yzzzzz_yyyyyz, g_z_0_yzzzzz_yyyyzz, g_z_0_yzzzzz_yyyzzz, g_z_0_yzzzzz_yyzzzz, g_z_0_yzzzzz_yzzzzz, g_z_0_yzzzzz_zzzzzz, g_z_0_zzzzz_xxxxxx, g_z_0_zzzzz_xxxxxxy, g_z_0_zzzzz_xxxxxy, g_z_0_zzzzz_xxxxxyy, g_z_0_zzzzz_xxxxxyz, g_z_0_zzzzz_xxxxxz, g_z_0_zzzzz_xxxxyy, g_z_0_zzzzz_xxxxyyy, g_z_0_zzzzz_xxxxyyz, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxyzz, g_z_0_zzzzz_xxxxzz, g_z_0_zzzzz_xxxyyy, g_z_0_zzzzz_xxxyyyy, g_z_0_zzzzz_xxxyyyz, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyyzz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxyzzz, g_z_0_zzzzz_xxxzzz, g_z_0_zzzzz_xxyyyy, g_z_0_zzzzz_xxyyyyy, g_z_0_zzzzz_xxyyyyz, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyyzz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyyzzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxyzzzz, g_z_0_zzzzz_xxzzzz, g_z_0_zzzzz_xyyyyy, g_z_0_zzzzz_xyyyyyy, g_z_0_zzzzz_xyyyyyz, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyyzz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyyzzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyyzzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xyzzzzz, g_z_0_zzzzz_xzzzzz, g_z_0_zzzzz_yyyyyy, g_z_0_zzzzz_yyyyyyy, g_z_0_zzzzz_yyyyyyz, g_z_0_zzzzz_yyyyyz, g_z_0_zzzzz_yyyyyzz, g_z_0_zzzzz_yyyyzz, g_z_0_zzzzz_yyyyzzz, g_z_0_zzzzz_yyyzzz, g_z_0_zzzzz_yyyzzzz, g_z_0_zzzzz_yyzzzz, g_z_0_zzzzz_yyzzzzz, g_z_0_zzzzz_yzzzzz, g_z_0_zzzzz_yzzzzzz, g_z_0_zzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_xxxxxx[k] = -g_z_0_zzzzz_xxxxxx[k] * ab_y + g_z_0_zzzzz_xxxxxxy[k];

                g_z_0_yzzzzz_xxxxxy[k] = -g_z_0_zzzzz_xxxxxy[k] * ab_y + g_z_0_zzzzz_xxxxxyy[k];

                g_z_0_yzzzzz_xxxxxz[k] = -g_z_0_zzzzz_xxxxxz[k] * ab_y + g_z_0_zzzzz_xxxxxyz[k];

                g_z_0_yzzzzz_xxxxyy[k] = -g_z_0_zzzzz_xxxxyy[k] * ab_y + g_z_0_zzzzz_xxxxyyy[k];

                g_z_0_yzzzzz_xxxxyz[k] = -g_z_0_zzzzz_xxxxyz[k] * ab_y + g_z_0_zzzzz_xxxxyyz[k];

                g_z_0_yzzzzz_xxxxzz[k] = -g_z_0_zzzzz_xxxxzz[k] * ab_y + g_z_0_zzzzz_xxxxyzz[k];

                g_z_0_yzzzzz_xxxyyy[k] = -g_z_0_zzzzz_xxxyyy[k] * ab_y + g_z_0_zzzzz_xxxyyyy[k];

                g_z_0_yzzzzz_xxxyyz[k] = -g_z_0_zzzzz_xxxyyz[k] * ab_y + g_z_0_zzzzz_xxxyyyz[k];

                g_z_0_yzzzzz_xxxyzz[k] = -g_z_0_zzzzz_xxxyzz[k] * ab_y + g_z_0_zzzzz_xxxyyzz[k];

                g_z_0_yzzzzz_xxxzzz[k] = -g_z_0_zzzzz_xxxzzz[k] * ab_y + g_z_0_zzzzz_xxxyzzz[k];

                g_z_0_yzzzzz_xxyyyy[k] = -g_z_0_zzzzz_xxyyyy[k] * ab_y + g_z_0_zzzzz_xxyyyyy[k];

                g_z_0_yzzzzz_xxyyyz[k] = -g_z_0_zzzzz_xxyyyz[k] * ab_y + g_z_0_zzzzz_xxyyyyz[k];

                g_z_0_yzzzzz_xxyyzz[k] = -g_z_0_zzzzz_xxyyzz[k] * ab_y + g_z_0_zzzzz_xxyyyzz[k];

                g_z_0_yzzzzz_xxyzzz[k] = -g_z_0_zzzzz_xxyzzz[k] * ab_y + g_z_0_zzzzz_xxyyzzz[k];

                g_z_0_yzzzzz_xxzzzz[k] = -g_z_0_zzzzz_xxzzzz[k] * ab_y + g_z_0_zzzzz_xxyzzzz[k];

                g_z_0_yzzzzz_xyyyyy[k] = -g_z_0_zzzzz_xyyyyy[k] * ab_y + g_z_0_zzzzz_xyyyyyy[k];

                g_z_0_yzzzzz_xyyyyz[k] = -g_z_0_zzzzz_xyyyyz[k] * ab_y + g_z_0_zzzzz_xyyyyyz[k];

                g_z_0_yzzzzz_xyyyzz[k] = -g_z_0_zzzzz_xyyyzz[k] * ab_y + g_z_0_zzzzz_xyyyyzz[k];

                g_z_0_yzzzzz_xyyzzz[k] = -g_z_0_zzzzz_xyyzzz[k] * ab_y + g_z_0_zzzzz_xyyyzzz[k];

                g_z_0_yzzzzz_xyzzzz[k] = -g_z_0_zzzzz_xyzzzz[k] * ab_y + g_z_0_zzzzz_xyyzzzz[k];

                g_z_0_yzzzzz_xzzzzz[k] = -g_z_0_zzzzz_xzzzzz[k] * ab_y + g_z_0_zzzzz_xyzzzzz[k];

                g_z_0_yzzzzz_yyyyyy[k] = -g_z_0_zzzzz_yyyyyy[k] * ab_y + g_z_0_zzzzz_yyyyyyy[k];

                g_z_0_yzzzzz_yyyyyz[k] = -g_z_0_zzzzz_yyyyyz[k] * ab_y + g_z_0_zzzzz_yyyyyyz[k];

                g_z_0_yzzzzz_yyyyzz[k] = -g_z_0_zzzzz_yyyyzz[k] * ab_y + g_z_0_zzzzz_yyyyyzz[k];

                g_z_0_yzzzzz_yyyzzz[k] = -g_z_0_zzzzz_yyyzzz[k] * ab_y + g_z_0_zzzzz_yyyyzzz[k];

                g_z_0_yzzzzz_yyzzzz[k] = -g_z_0_zzzzz_yyzzzz[k] * ab_y + g_z_0_zzzzz_yyyzzzz[k];

                g_z_0_yzzzzz_yzzzzz[k] = -g_z_0_zzzzz_yzzzzz[k] * ab_y + g_z_0_zzzzz_yyzzzzz[k];

                g_z_0_yzzzzz_zzzzzz[k] = -g_z_0_zzzzz_zzzzzz[k] * ab_y + g_z_0_zzzzz_yzzzzzz[k];
            }

            /// Set up 2324-2352 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_xxxxxx = cbuffer.data(ii_geom_10_off + 2324 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxxxy = cbuffer.data(ii_geom_10_off + 2325 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxxxz = cbuffer.data(ii_geom_10_off + 2326 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxxyy = cbuffer.data(ii_geom_10_off + 2327 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxxyz = cbuffer.data(ii_geom_10_off + 2328 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxxzz = cbuffer.data(ii_geom_10_off + 2329 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxyyy = cbuffer.data(ii_geom_10_off + 2330 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxyyz = cbuffer.data(ii_geom_10_off + 2331 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxyzz = cbuffer.data(ii_geom_10_off + 2332 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxzzz = cbuffer.data(ii_geom_10_off + 2333 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxyyyy = cbuffer.data(ii_geom_10_off + 2334 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxyyyz = cbuffer.data(ii_geom_10_off + 2335 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxyyzz = cbuffer.data(ii_geom_10_off + 2336 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxyzzz = cbuffer.data(ii_geom_10_off + 2337 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxzzzz = cbuffer.data(ii_geom_10_off + 2338 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyyyyy = cbuffer.data(ii_geom_10_off + 2339 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyyyyz = cbuffer.data(ii_geom_10_off + 2340 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyyyzz = cbuffer.data(ii_geom_10_off + 2341 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyyzzz = cbuffer.data(ii_geom_10_off + 2342 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyzzzz = cbuffer.data(ii_geom_10_off + 2343 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xzzzzz = cbuffer.data(ii_geom_10_off + 2344 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyyyyy = cbuffer.data(ii_geom_10_off + 2345 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyyyyz = cbuffer.data(ii_geom_10_off + 2346 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyyyzz = cbuffer.data(ii_geom_10_off + 2347 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyyzzz = cbuffer.data(ii_geom_10_off + 2348 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyzzzz = cbuffer.data(ii_geom_10_off + 2349 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yzzzzz = cbuffer.data(ii_geom_10_off + 2350 * ccomps * dcomps);

            auto g_z_0_zzzzzz_zzzzzz = cbuffer.data(ii_geom_10_off + 2351 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzzz_xxxxxx, g_z_0_zzzzz_xxxxxxz, g_z_0_zzzzz_xxxxxy, g_z_0_zzzzz_xxxxxyz, g_z_0_zzzzz_xxxxxz, g_z_0_zzzzz_xxxxxzz, g_z_0_zzzzz_xxxxyy, g_z_0_zzzzz_xxxxyyz, g_z_0_zzzzz_xxxxyz, g_z_0_zzzzz_xxxxyzz, g_z_0_zzzzz_xxxxzz, g_z_0_zzzzz_xxxxzzz, g_z_0_zzzzz_xxxyyy, g_z_0_zzzzz_xxxyyyz, g_z_0_zzzzz_xxxyyz, g_z_0_zzzzz_xxxyyzz, g_z_0_zzzzz_xxxyzz, g_z_0_zzzzz_xxxyzzz, g_z_0_zzzzz_xxxzzz, g_z_0_zzzzz_xxxzzzz, g_z_0_zzzzz_xxyyyy, g_z_0_zzzzz_xxyyyyz, g_z_0_zzzzz_xxyyyz, g_z_0_zzzzz_xxyyyzz, g_z_0_zzzzz_xxyyzz, g_z_0_zzzzz_xxyyzzz, g_z_0_zzzzz_xxyzzz, g_z_0_zzzzz_xxyzzzz, g_z_0_zzzzz_xxzzzz, g_z_0_zzzzz_xxzzzzz, g_z_0_zzzzz_xyyyyy, g_z_0_zzzzz_xyyyyyz, g_z_0_zzzzz_xyyyyz, g_z_0_zzzzz_xyyyyzz, g_z_0_zzzzz_xyyyzz, g_z_0_zzzzz_xyyyzzz, g_z_0_zzzzz_xyyzzz, g_z_0_zzzzz_xyyzzzz, g_z_0_zzzzz_xyzzzz, g_z_0_zzzzz_xyzzzzz, g_z_0_zzzzz_xzzzzz, g_z_0_zzzzz_xzzzzzz, g_z_0_zzzzz_yyyyyy, g_z_0_zzzzz_yyyyyyz, g_z_0_zzzzz_yyyyyz, g_z_0_zzzzz_yyyyyzz, g_z_0_zzzzz_yyyyzz, g_z_0_zzzzz_yyyyzzz, g_z_0_zzzzz_yyyzzz, g_z_0_zzzzz_yyyzzzz, g_z_0_zzzzz_yyzzzz, g_z_0_zzzzz_yyzzzzz, g_z_0_zzzzz_yzzzzz, g_z_0_zzzzz_yzzzzzz, g_z_0_zzzzz_zzzzzz, g_z_0_zzzzz_zzzzzzz, g_z_0_zzzzzz_xxxxxx, g_z_0_zzzzzz_xxxxxy, g_z_0_zzzzzz_xxxxxz, g_z_0_zzzzzz_xxxxyy, g_z_0_zzzzzz_xxxxyz, g_z_0_zzzzzz_xxxxzz, g_z_0_zzzzzz_xxxyyy, g_z_0_zzzzzz_xxxyyz, g_z_0_zzzzzz_xxxyzz, g_z_0_zzzzzz_xxxzzz, g_z_0_zzzzzz_xxyyyy, g_z_0_zzzzzz_xxyyyz, g_z_0_zzzzzz_xxyyzz, g_z_0_zzzzzz_xxyzzz, g_z_0_zzzzzz_xxzzzz, g_z_0_zzzzzz_xyyyyy, g_z_0_zzzzzz_xyyyyz, g_z_0_zzzzzz_xyyyzz, g_z_0_zzzzzz_xyyzzz, g_z_0_zzzzzz_xyzzzz, g_z_0_zzzzzz_xzzzzz, g_z_0_zzzzzz_yyyyyy, g_z_0_zzzzzz_yyyyyz, g_z_0_zzzzzz_yyyyzz, g_z_0_zzzzzz_yyyzzz, g_z_0_zzzzzz_yyzzzz, g_z_0_zzzzzz_yzzzzz, g_z_0_zzzzzz_zzzzzz, g_zzzzz_xxxxxx, g_zzzzz_xxxxxy, g_zzzzz_xxxxxz, g_zzzzz_xxxxyy, g_zzzzz_xxxxyz, g_zzzzz_xxxxzz, g_zzzzz_xxxyyy, g_zzzzz_xxxyyz, g_zzzzz_xxxyzz, g_zzzzz_xxxzzz, g_zzzzz_xxyyyy, g_zzzzz_xxyyyz, g_zzzzz_xxyyzz, g_zzzzz_xxyzzz, g_zzzzz_xxzzzz, g_zzzzz_xyyyyy, g_zzzzz_xyyyyz, g_zzzzz_xyyyzz, g_zzzzz_xyyzzz, g_zzzzz_xyzzzz, g_zzzzz_xzzzzz, g_zzzzz_yyyyyy, g_zzzzz_yyyyyz, g_zzzzz_yyyyzz, g_zzzzz_yyyzzz, g_zzzzz_yyzzzz, g_zzzzz_yzzzzz, g_zzzzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_xxxxxx[k] = -g_zzzzz_xxxxxx[k] - g_z_0_zzzzz_xxxxxx[k] * ab_z + g_z_0_zzzzz_xxxxxxz[k];

                g_z_0_zzzzzz_xxxxxy[k] = -g_zzzzz_xxxxxy[k] - g_z_0_zzzzz_xxxxxy[k] * ab_z + g_z_0_zzzzz_xxxxxyz[k];

                g_z_0_zzzzzz_xxxxxz[k] = -g_zzzzz_xxxxxz[k] - g_z_0_zzzzz_xxxxxz[k] * ab_z + g_z_0_zzzzz_xxxxxzz[k];

                g_z_0_zzzzzz_xxxxyy[k] = -g_zzzzz_xxxxyy[k] - g_z_0_zzzzz_xxxxyy[k] * ab_z + g_z_0_zzzzz_xxxxyyz[k];

                g_z_0_zzzzzz_xxxxyz[k] = -g_zzzzz_xxxxyz[k] - g_z_0_zzzzz_xxxxyz[k] * ab_z + g_z_0_zzzzz_xxxxyzz[k];

                g_z_0_zzzzzz_xxxxzz[k] = -g_zzzzz_xxxxzz[k] - g_z_0_zzzzz_xxxxzz[k] * ab_z + g_z_0_zzzzz_xxxxzzz[k];

                g_z_0_zzzzzz_xxxyyy[k] = -g_zzzzz_xxxyyy[k] - g_z_0_zzzzz_xxxyyy[k] * ab_z + g_z_0_zzzzz_xxxyyyz[k];

                g_z_0_zzzzzz_xxxyyz[k] = -g_zzzzz_xxxyyz[k] - g_z_0_zzzzz_xxxyyz[k] * ab_z + g_z_0_zzzzz_xxxyyzz[k];

                g_z_0_zzzzzz_xxxyzz[k] = -g_zzzzz_xxxyzz[k] - g_z_0_zzzzz_xxxyzz[k] * ab_z + g_z_0_zzzzz_xxxyzzz[k];

                g_z_0_zzzzzz_xxxzzz[k] = -g_zzzzz_xxxzzz[k] - g_z_0_zzzzz_xxxzzz[k] * ab_z + g_z_0_zzzzz_xxxzzzz[k];

                g_z_0_zzzzzz_xxyyyy[k] = -g_zzzzz_xxyyyy[k] - g_z_0_zzzzz_xxyyyy[k] * ab_z + g_z_0_zzzzz_xxyyyyz[k];

                g_z_0_zzzzzz_xxyyyz[k] = -g_zzzzz_xxyyyz[k] - g_z_0_zzzzz_xxyyyz[k] * ab_z + g_z_0_zzzzz_xxyyyzz[k];

                g_z_0_zzzzzz_xxyyzz[k] = -g_zzzzz_xxyyzz[k] - g_z_0_zzzzz_xxyyzz[k] * ab_z + g_z_0_zzzzz_xxyyzzz[k];

                g_z_0_zzzzzz_xxyzzz[k] = -g_zzzzz_xxyzzz[k] - g_z_0_zzzzz_xxyzzz[k] * ab_z + g_z_0_zzzzz_xxyzzzz[k];

                g_z_0_zzzzzz_xxzzzz[k] = -g_zzzzz_xxzzzz[k] - g_z_0_zzzzz_xxzzzz[k] * ab_z + g_z_0_zzzzz_xxzzzzz[k];

                g_z_0_zzzzzz_xyyyyy[k] = -g_zzzzz_xyyyyy[k] - g_z_0_zzzzz_xyyyyy[k] * ab_z + g_z_0_zzzzz_xyyyyyz[k];

                g_z_0_zzzzzz_xyyyyz[k] = -g_zzzzz_xyyyyz[k] - g_z_0_zzzzz_xyyyyz[k] * ab_z + g_z_0_zzzzz_xyyyyzz[k];

                g_z_0_zzzzzz_xyyyzz[k] = -g_zzzzz_xyyyzz[k] - g_z_0_zzzzz_xyyyzz[k] * ab_z + g_z_0_zzzzz_xyyyzzz[k];

                g_z_0_zzzzzz_xyyzzz[k] = -g_zzzzz_xyyzzz[k] - g_z_0_zzzzz_xyyzzz[k] * ab_z + g_z_0_zzzzz_xyyzzzz[k];

                g_z_0_zzzzzz_xyzzzz[k] = -g_zzzzz_xyzzzz[k] - g_z_0_zzzzz_xyzzzz[k] * ab_z + g_z_0_zzzzz_xyzzzzz[k];

                g_z_0_zzzzzz_xzzzzz[k] = -g_zzzzz_xzzzzz[k] - g_z_0_zzzzz_xzzzzz[k] * ab_z + g_z_0_zzzzz_xzzzzzz[k];

                g_z_0_zzzzzz_yyyyyy[k] = -g_zzzzz_yyyyyy[k] - g_z_0_zzzzz_yyyyyy[k] * ab_z + g_z_0_zzzzz_yyyyyyz[k];

                g_z_0_zzzzzz_yyyyyz[k] = -g_zzzzz_yyyyyz[k] - g_z_0_zzzzz_yyyyyz[k] * ab_z + g_z_0_zzzzz_yyyyyzz[k];

                g_z_0_zzzzzz_yyyyzz[k] = -g_zzzzz_yyyyzz[k] - g_z_0_zzzzz_yyyyzz[k] * ab_z + g_z_0_zzzzz_yyyyzzz[k];

                g_z_0_zzzzzz_yyyzzz[k] = -g_zzzzz_yyyzzz[k] - g_z_0_zzzzz_yyyzzz[k] * ab_z + g_z_0_zzzzz_yyyzzzz[k];

                g_z_0_zzzzzz_yyzzzz[k] = -g_zzzzz_yyzzzz[k] - g_z_0_zzzzz_yyzzzz[k] * ab_z + g_z_0_zzzzz_yyzzzzz[k];

                g_z_0_zzzzzz_yzzzzz[k] = -g_zzzzz_yzzzzz[k] - g_z_0_zzzzz_yzzzzz[k] * ab_z + g_z_0_zzzzz_yzzzzzz[k];

                g_z_0_zzzzzz_zzzzzz[k] = -g_zzzzz_zzzzzz[k] - g_z_0_zzzzz_zzzzzz[k] * ab_z + g_z_0_zzzzz_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

