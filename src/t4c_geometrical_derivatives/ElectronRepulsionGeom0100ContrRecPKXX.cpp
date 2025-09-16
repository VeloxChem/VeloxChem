#include "ElectronRepulsionGeom0100ContrRecPKXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_pkxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_pkxx,
                                            const size_t idx_skxx,
                                            const size_t idx_geom_01_skxx,
                                            const size_t idx_geom_01_slxx,
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

            /// Set up components of auxilary buffer : SKSS

            const auto sk_geom_01_off = idx_geom_01_skxx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : SLSS

            const auto sl_geom_01_off = idx_geom_01_slxx + i * dcomps + j;

            auto g_0_x_0_xxxxxxxx = cbuffer.data(sl_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_0_xxxxxxxy = cbuffer.data(sl_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_0_xxxxxxxz = cbuffer.data(sl_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_0_xxxxxxyy = cbuffer.data(sl_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_0_xxxxxxyz = cbuffer.data(sl_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_0_xxxxxxzz = cbuffer.data(sl_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_0_xxxxxyyy = cbuffer.data(sl_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_0_xxxxxyyz = cbuffer.data(sl_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_0_xxxxxyzz = cbuffer.data(sl_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_0_xxxxxzzz = cbuffer.data(sl_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_0_xxxxyyyy = cbuffer.data(sl_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_0_xxxxyyyz = cbuffer.data(sl_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_0_xxxxyyzz = cbuffer.data(sl_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_0_xxxxyzzz = cbuffer.data(sl_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_0_xxxxzzzz = cbuffer.data(sl_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_0_xxxyyyyy = cbuffer.data(sl_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_0_xxxyyyyz = cbuffer.data(sl_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_0_xxxyyyzz = cbuffer.data(sl_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_0_xxxyyzzz = cbuffer.data(sl_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_0_xxxyzzzz = cbuffer.data(sl_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_0_xxxzzzzz = cbuffer.data(sl_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_0_xxyyyyyy = cbuffer.data(sl_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_0_xxyyyyyz = cbuffer.data(sl_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_0_xxyyyyzz = cbuffer.data(sl_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_0_xxyyyzzz = cbuffer.data(sl_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_0_xxyyzzzz = cbuffer.data(sl_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_0_xxyzzzzz = cbuffer.data(sl_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_0_xxzzzzzz = cbuffer.data(sl_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_0_xyyyyyyy = cbuffer.data(sl_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_0_xyyyyyyz = cbuffer.data(sl_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_0_xyyyyyzz = cbuffer.data(sl_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_0_xyyyyzzz = cbuffer.data(sl_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_0_xyyyzzzz = cbuffer.data(sl_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_0_xyyzzzzz = cbuffer.data(sl_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_0_xyzzzzzz = cbuffer.data(sl_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_0_xzzzzzzz = cbuffer.data(sl_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_0_yyyyyyyy = cbuffer.data(sl_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_0_yyyyyyyz = cbuffer.data(sl_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_0_yyyyyyzz = cbuffer.data(sl_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_0_yyyyyzzz = cbuffer.data(sl_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_0_yyyyzzzz = cbuffer.data(sl_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_0_yyyzzzzz = cbuffer.data(sl_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_0_yyzzzzzz = cbuffer.data(sl_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_0_yzzzzzzz = cbuffer.data(sl_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_0_zzzzzzzz = cbuffer.data(sl_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_y_0_xxxxxxxx = cbuffer.data(sl_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_0_xxxxxxxy = cbuffer.data(sl_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_0_xxxxxxxz = cbuffer.data(sl_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_y_0_xxxxxxyy = cbuffer.data(sl_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_0_xxxxxxyz = cbuffer.data(sl_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_0_xxxxxxzz = cbuffer.data(sl_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_y_0_xxxxxyyy = cbuffer.data(sl_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_0_xxxxxyyz = cbuffer.data(sl_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_0_xxxxxyzz = cbuffer.data(sl_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_y_0_xxxxxzzz = cbuffer.data(sl_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_0_xxxxyyyy = cbuffer.data(sl_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_y_0_xxxxyyyz = cbuffer.data(sl_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_y_0_xxxxyyzz = cbuffer.data(sl_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_y_0_xxxxyzzz = cbuffer.data(sl_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_y_0_xxxxzzzz = cbuffer.data(sl_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_y_0_xxxyyyyy = cbuffer.data(sl_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_y_0_xxxyyyyz = cbuffer.data(sl_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_y_0_xxxyyyzz = cbuffer.data(sl_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_y_0_xxxyyzzz = cbuffer.data(sl_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_0_xxxyzzzz = cbuffer.data(sl_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_0_xxxzzzzz = cbuffer.data(sl_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_y_0_xxyyyyyy = cbuffer.data(sl_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_0_xxyyyyyz = cbuffer.data(sl_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_0_xxyyyyzz = cbuffer.data(sl_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_y_0_xxyyyzzz = cbuffer.data(sl_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_0_xxyyzzzz = cbuffer.data(sl_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_0_xxyzzzzz = cbuffer.data(sl_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_y_0_xxzzzzzz = cbuffer.data(sl_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_y_0_xyyyyyyy = cbuffer.data(sl_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_y_0_xyyyyyyz = cbuffer.data(sl_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_y_0_xyyyyyzz = cbuffer.data(sl_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_y_0_xyyyyzzz = cbuffer.data(sl_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_y_0_xyyyzzzz = cbuffer.data(sl_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_y_0_xyyzzzzz = cbuffer.data(sl_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_y_0_xyzzzzzz = cbuffer.data(sl_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_y_0_xzzzzzzz = cbuffer.data(sl_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_y_0_yyyyyyyy = cbuffer.data(sl_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_y_0_yyyyyyyz = cbuffer.data(sl_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_y_0_yyyyyyzz = cbuffer.data(sl_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_y_0_yyyyyzzz = cbuffer.data(sl_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_0_yyyyzzzz = cbuffer.data(sl_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_0_yyyzzzzz = cbuffer.data(sl_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_0_yyzzzzzz = cbuffer.data(sl_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_0_yzzzzzzz = cbuffer.data(sl_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_0_zzzzzzzz = cbuffer.data(sl_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_z_0_xxxxxxxx = cbuffer.data(sl_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_z_0_xxxxxxxy = cbuffer.data(sl_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_z_0_xxxxxxxz = cbuffer.data(sl_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_z_0_xxxxxxyy = cbuffer.data(sl_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_z_0_xxxxxxyz = cbuffer.data(sl_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_z_0_xxxxxxzz = cbuffer.data(sl_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_z_0_xxxxxyyy = cbuffer.data(sl_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_z_0_xxxxxyyz = cbuffer.data(sl_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_z_0_xxxxxyzz = cbuffer.data(sl_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_z_0_xxxxxzzz = cbuffer.data(sl_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_z_0_xxxxyyyy = cbuffer.data(sl_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_z_0_xxxxyyyz = cbuffer.data(sl_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_z_0_xxxxyyzz = cbuffer.data(sl_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_z_0_xxxxyzzz = cbuffer.data(sl_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_z_0_xxxxzzzz = cbuffer.data(sl_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_z_0_xxxyyyyy = cbuffer.data(sl_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_z_0_xxxyyyyz = cbuffer.data(sl_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_z_0_xxxyyyzz = cbuffer.data(sl_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_z_0_xxxyyzzz = cbuffer.data(sl_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_z_0_xxxyzzzz = cbuffer.data(sl_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_z_0_xxxzzzzz = cbuffer.data(sl_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_z_0_xxyyyyyy = cbuffer.data(sl_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_z_0_xxyyyyyz = cbuffer.data(sl_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_z_0_xxyyyyzz = cbuffer.data(sl_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_z_0_xxyyyzzz = cbuffer.data(sl_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_z_0_xxyyzzzz = cbuffer.data(sl_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_z_0_xxyzzzzz = cbuffer.data(sl_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_z_0_xxzzzzzz = cbuffer.data(sl_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_z_0_xyyyyyyy = cbuffer.data(sl_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_z_0_xyyyyyyz = cbuffer.data(sl_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_z_0_xyyyyyzz = cbuffer.data(sl_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_z_0_xyyyyzzz = cbuffer.data(sl_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_z_0_xyyyzzzz = cbuffer.data(sl_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_z_0_xyyzzzzz = cbuffer.data(sl_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_z_0_xyzzzzzz = cbuffer.data(sl_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_z_0_xzzzzzzz = cbuffer.data(sl_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_z_0_yyyyyyyy = cbuffer.data(sl_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_z_0_yyyyyyyz = cbuffer.data(sl_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_z_0_yyyyyyzz = cbuffer.data(sl_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_z_0_yyyyyzzz = cbuffer.data(sl_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_z_0_yyyyzzzz = cbuffer.data(sl_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_z_0_yyyzzzzz = cbuffer.data(sl_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_z_0_yyzzzzzz = cbuffer.data(sl_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_z_0_yzzzzzzz = cbuffer.data(sl_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_z_0_zzzzzzzz = cbuffer.data(sl_geom_01_off + 134 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pkxx

            const auto pk_geom_01_off = idx_geom_01_pkxx + i * dcomps + j;

            /// Set up 0-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_x_xxxxxxx = cbuffer.data(pk_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xxxxxxy = cbuffer.data(pk_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xxxxxxz = cbuffer.data(pk_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_xxxxxyy = cbuffer.data(pk_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_xxxxxyz = cbuffer.data(pk_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_xxxxxzz = cbuffer.data(pk_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_x_xxxxyyy = cbuffer.data(pk_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_x_xxxxyyz = cbuffer.data(pk_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_x_xxxxyzz = cbuffer.data(pk_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_x_xxxxzzz = cbuffer.data(pk_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_x_xxxyyyy = cbuffer.data(pk_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_x_xxxyyyz = cbuffer.data(pk_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_x_xxxyyzz = cbuffer.data(pk_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_x_xxxyzzz = cbuffer.data(pk_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_x_xxxzzzz = cbuffer.data(pk_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_x_xxyyyyy = cbuffer.data(pk_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_x_xxyyyyz = cbuffer.data(pk_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_x_xxyyyzz = cbuffer.data(pk_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_x_xxyyzzz = cbuffer.data(pk_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_x_xxyzzzz = cbuffer.data(pk_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_x_xxzzzzz = cbuffer.data(pk_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_x_xyyyyyy = cbuffer.data(pk_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_x_xyyyyyz = cbuffer.data(pk_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_x_xyyyyzz = cbuffer.data(pk_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_x_xyyyzzz = cbuffer.data(pk_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_x_xyyzzzz = cbuffer.data(pk_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_x_xyzzzzz = cbuffer.data(pk_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_x_xzzzzzz = cbuffer.data(pk_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_x_yyyyyyy = cbuffer.data(pk_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_x_yyyyyyz = cbuffer.data(pk_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_x_yyyyyzz = cbuffer.data(pk_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_x_yyyyzzz = cbuffer.data(pk_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_x_yyyzzzz = cbuffer.data(pk_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_x_yyzzzzz = cbuffer.data(pk_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_x_yzzzzzz = cbuffer.data(pk_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_x_zzzzzzz = cbuffer.data(pk_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxxx, g_0_x_0_xxxxxxxx, g_0_x_0_xxxxxxxy, g_0_x_0_xxxxxxxz, g_0_x_0_xxxxxxy, g_0_x_0_xxxxxxyy, g_0_x_0_xxxxxxyz, g_0_x_0_xxxxxxz, g_0_x_0_xxxxxxzz, g_0_x_0_xxxxxyy, g_0_x_0_xxxxxyyy, g_0_x_0_xxxxxyyz, g_0_x_0_xxxxxyz, g_0_x_0_xxxxxyzz, g_0_x_0_xxxxxzz, g_0_x_0_xxxxxzzz, g_0_x_0_xxxxyyy, g_0_x_0_xxxxyyyy, g_0_x_0_xxxxyyyz, g_0_x_0_xxxxyyz, g_0_x_0_xxxxyyzz, g_0_x_0_xxxxyzz, g_0_x_0_xxxxyzzz, g_0_x_0_xxxxzzz, g_0_x_0_xxxxzzzz, g_0_x_0_xxxyyyy, g_0_x_0_xxxyyyyy, g_0_x_0_xxxyyyyz, g_0_x_0_xxxyyyz, g_0_x_0_xxxyyyzz, g_0_x_0_xxxyyzz, g_0_x_0_xxxyyzzz, g_0_x_0_xxxyzzz, g_0_x_0_xxxyzzzz, g_0_x_0_xxxzzzz, g_0_x_0_xxxzzzzz, g_0_x_0_xxyyyyy, g_0_x_0_xxyyyyyy, g_0_x_0_xxyyyyyz, g_0_x_0_xxyyyyz, g_0_x_0_xxyyyyzz, g_0_x_0_xxyyyzz, g_0_x_0_xxyyyzzz, g_0_x_0_xxyyzzz, g_0_x_0_xxyyzzzz, g_0_x_0_xxyzzzz, g_0_x_0_xxyzzzzz, g_0_x_0_xxzzzzz, g_0_x_0_xxzzzzzz, g_0_x_0_xyyyyyy, g_0_x_0_xyyyyyyy, g_0_x_0_xyyyyyyz, g_0_x_0_xyyyyyz, g_0_x_0_xyyyyyzz, g_0_x_0_xyyyyzz, g_0_x_0_xyyyyzzz, g_0_x_0_xyyyzzz, g_0_x_0_xyyyzzzz, g_0_x_0_xyyzzzz, g_0_x_0_xyyzzzzz, g_0_x_0_xyzzzzz, g_0_x_0_xyzzzzzz, g_0_x_0_xzzzzzz, g_0_x_0_xzzzzzzz, g_0_x_0_yyyyyyy, g_0_x_0_yyyyyyz, g_0_x_0_yyyyyzz, g_0_x_0_yyyyzzz, g_0_x_0_yyyzzzz, g_0_x_0_yyzzzzz, g_0_x_0_yzzzzzz, g_0_x_0_zzzzzzz, g_0_x_x_xxxxxxx, g_0_x_x_xxxxxxy, g_0_x_x_xxxxxxz, g_0_x_x_xxxxxyy, g_0_x_x_xxxxxyz, g_0_x_x_xxxxxzz, g_0_x_x_xxxxyyy, g_0_x_x_xxxxyyz, g_0_x_x_xxxxyzz, g_0_x_x_xxxxzzz, g_0_x_x_xxxyyyy, g_0_x_x_xxxyyyz, g_0_x_x_xxxyyzz, g_0_x_x_xxxyzzz, g_0_x_x_xxxzzzz, g_0_x_x_xxyyyyy, g_0_x_x_xxyyyyz, g_0_x_x_xxyyyzz, g_0_x_x_xxyyzzz, g_0_x_x_xxyzzzz, g_0_x_x_xxzzzzz, g_0_x_x_xyyyyyy, g_0_x_x_xyyyyyz, g_0_x_x_xyyyyzz, g_0_x_x_xyyyzzz, g_0_x_x_xyyzzzz, g_0_x_x_xyzzzzz, g_0_x_x_xzzzzzz, g_0_x_x_yyyyyyy, g_0_x_x_yyyyyyz, g_0_x_x_yyyyyzz, g_0_x_x_yyyyzzz, g_0_x_x_yyyzzzz, g_0_x_x_yyzzzzz, g_0_x_x_yzzzzzz, g_0_x_x_zzzzzzz, g_0_xxxxxxx, g_0_xxxxxxy, g_0_xxxxxxz, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxzz, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyzz, g_0_xxxxzzz, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyzz, g_0_xxxyzzz, g_0_xxxzzzz, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyzz, g_0_xxyyzzz, g_0_xxyzzzz, g_0_xxzzzzz, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyzz, g_0_xyyyzzz, g_0_xyyzzzz, g_0_xyzzzzz, g_0_xzzzzzz, g_0_yyyyyyy, g_0_yyyyyyz, g_0_yyyyyzz, g_0_yyyyzzz, g_0_yyyzzzz, g_0_yyzzzzz, g_0_yzzzzzz, g_0_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_x_xxxxxxx[k] = g_0_xxxxxxx[k] - g_0_x_0_xxxxxxx[k] * ab_x + g_0_x_0_xxxxxxxx[k];

                g_0_x_x_xxxxxxy[k] = g_0_xxxxxxy[k] - g_0_x_0_xxxxxxy[k] * ab_x + g_0_x_0_xxxxxxxy[k];

                g_0_x_x_xxxxxxz[k] = g_0_xxxxxxz[k] - g_0_x_0_xxxxxxz[k] * ab_x + g_0_x_0_xxxxxxxz[k];

                g_0_x_x_xxxxxyy[k] = g_0_xxxxxyy[k] - g_0_x_0_xxxxxyy[k] * ab_x + g_0_x_0_xxxxxxyy[k];

                g_0_x_x_xxxxxyz[k] = g_0_xxxxxyz[k] - g_0_x_0_xxxxxyz[k] * ab_x + g_0_x_0_xxxxxxyz[k];

                g_0_x_x_xxxxxzz[k] = g_0_xxxxxzz[k] - g_0_x_0_xxxxxzz[k] * ab_x + g_0_x_0_xxxxxxzz[k];

                g_0_x_x_xxxxyyy[k] = g_0_xxxxyyy[k] - g_0_x_0_xxxxyyy[k] * ab_x + g_0_x_0_xxxxxyyy[k];

                g_0_x_x_xxxxyyz[k] = g_0_xxxxyyz[k] - g_0_x_0_xxxxyyz[k] * ab_x + g_0_x_0_xxxxxyyz[k];

                g_0_x_x_xxxxyzz[k] = g_0_xxxxyzz[k] - g_0_x_0_xxxxyzz[k] * ab_x + g_0_x_0_xxxxxyzz[k];

                g_0_x_x_xxxxzzz[k] = g_0_xxxxzzz[k] - g_0_x_0_xxxxzzz[k] * ab_x + g_0_x_0_xxxxxzzz[k];

                g_0_x_x_xxxyyyy[k] = g_0_xxxyyyy[k] - g_0_x_0_xxxyyyy[k] * ab_x + g_0_x_0_xxxxyyyy[k];

                g_0_x_x_xxxyyyz[k] = g_0_xxxyyyz[k] - g_0_x_0_xxxyyyz[k] * ab_x + g_0_x_0_xxxxyyyz[k];

                g_0_x_x_xxxyyzz[k] = g_0_xxxyyzz[k] - g_0_x_0_xxxyyzz[k] * ab_x + g_0_x_0_xxxxyyzz[k];

                g_0_x_x_xxxyzzz[k] = g_0_xxxyzzz[k] - g_0_x_0_xxxyzzz[k] * ab_x + g_0_x_0_xxxxyzzz[k];

                g_0_x_x_xxxzzzz[k] = g_0_xxxzzzz[k] - g_0_x_0_xxxzzzz[k] * ab_x + g_0_x_0_xxxxzzzz[k];

                g_0_x_x_xxyyyyy[k] = g_0_xxyyyyy[k] - g_0_x_0_xxyyyyy[k] * ab_x + g_0_x_0_xxxyyyyy[k];

                g_0_x_x_xxyyyyz[k] = g_0_xxyyyyz[k] - g_0_x_0_xxyyyyz[k] * ab_x + g_0_x_0_xxxyyyyz[k];

                g_0_x_x_xxyyyzz[k] = g_0_xxyyyzz[k] - g_0_x_0_xxyyyzz[k] * ab_x + g_0_x_0_xxxyyyzz[k];

                g_0_x_x_xxyyzzz[k] = g_0_xxyyzzz[k] - g_0_x_0_xxyyzzz[k] * ab_x + g_0_x_0_xxxyyzzz[k];

                g_0_x_x_xxyzzzz[k] = g_0_xxyzzzz[k] - g_0_x_0_xxyzzzz[k] * ab_x + g_0_x_0_xxxyzzzz[k];

                g_0_x_x_xxzzzzz[k] = g_0_xxzzzzz[k] - g_0_x_0_xxzzzzz[k] * ab_x + g_0_x_0_xxxzzzzz[k];

                g_0_x_x_xyyyyyy[k] = g_0_xyyyyyy[k] - g_0_x_0_xyyyyyy[k] * ab_x + g_0_x_0_xxyyyyyy[k];

                g_0_x_x_xyyyyyz[k] = g_0_xyyyyyz[k] - g_0_x_0_xyyyyyz[k] * ab_x + g_0_x_0_xxyyyyyz[k];

                g_0_x_x_xyyyyzz[k] = g_0_xyyyyzz[k] - g_0_x_0_xyyyyzz[k] * ab_x + g_0_x_0_xxyyyyzz[k];

                g_0_x_x_xyyyzzz[k] = g_0_xyyyzzz[k] - g_0_x_0_xyyyzzz[k] * ab_x + g_0_x_0_xxyyyzzz[k];

                g_0_x_x_xyyzzzz[k] = g_0_xyyzzzz[k] - g_0_x_0_xyyzzzz[k] * ab_x + g_0_x_0_xxyyzzzz[k];

                g_0_x_x_xyzzzzz[k] = g_0_xyzzzzz[k] - g_0_x_0_xyzzzzz[k] * ab_x + g_0_x_0_xxyzzzzz[k];

                g_0_x_x_xzzzzzz[k] = g_0_xzzzzzz[k] - g_0_x_0_xzzzzzz[k] * ab_x + g_0_x_0_xxzzzzzz[k];

                g_0_x_x_yyyyyyy[k] = g_0_yyyyyyy[k] - g_0_x_0_yyyyyyy[k] * ab_x + g_0_x_0_xyyyyyyy[k];

                g_0_x_x_yyyyyyz[k] = g_0_yyyyyyz[k] - g_0_x_0_yyyyyyz[k] * ab_x + g_0_x_0_xyyyyyyz[k];

                g_0_x_x_yyyyyzz[k] = g_0_yyyyyzz[k] - g_0_x_0_yyyyyzz[k] * ab_x + g_0_x_0_xyyyyyzz[k];

                g_0_x_x_yyyyzzz[k] = g_0_yyyyzzz[k] - g_0_x_0_yyyyzzz[k] * ab_x + g_0_x_0_xyyyyzzz[k];

                g_0_x_x_yyyzzzz[k] = g_0_yyyzzzz[k] - g_0_x_0_yyyzzzz[k] * ab_x + g_0_x_0_xyyyzzzz[k];

                g_0_x_x_yyzzzzz[k] = g_0_yyzzzzz[k] - g_0_x_0_yyzzzzz[k] * ab_x + g_0_x_0_xyyzzzzz[k];

                g_0_x_x_yzzzzzz[k] = g_0_yzzzzzz[k] - g_0_x_0_yzzzzzz[k] * ab_x + g_0_x_0_xyzzzzzz[k];

                g_0_x_x_zzzzzzz[k] = g_0_zzzzzzz[k] - g_0_x_0_zzzzzzz[k] * ab_x + g_0_x_0_xzzzzzzz[k];
            }

            /// Set up 36-72 components of targeted buffer : cbuffer.data(

            auto g_0_x_y_xxxxxxx = cbuffer.data(pk_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxy = cbuffer.data(pk_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_y_xxxxxxz = cbuffer.data(pk_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_y_xxxxxyy = cbuffer.data(pk_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_y_xxxxxyz = cbuffer.data(pk_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_y_xxxxxzz = cbuffer.data(pk_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_y_xxxxyyy = cbuffer.data(pk_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_y_xxxxyyz = cbuffer.data(pk_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_y_xxxxyzz = cbuffer.data(pk_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_y_xxxxzzz = cbuffer.data(pk_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_y_xxxyyyy = cbuffer.data(pk_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_y_xxxyyyz = cbuffer.data(pk_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_y_xxxyyzz = cbuffer.data(pk_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_y_xxxyzzz = cbuffer.data(pk_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_y_xxxzzzz = cbuffer.data(pk_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_y_xxyyyyy = cbuffer.data(pk_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_y_xxyyyyz = cbuffer.data(pk_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_y_xxyyyzz = cbuffer.data(pk_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_y_xxyyzzz = cbuffer.data(pk_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_y_xxyzzzz = cbuffer.data(pk_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_y_xxzzzzz = cbuffer.data(pk_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_y_xyyyyyy = cbuffer.data(pk_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_y_xyyyyyz = cbuffer.data(pk_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_y_xyyyyzz = cbuffer.data(pk_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_y_xyyyzzz = cbuffer.data(pk_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_y_xyyzzzz = cbuffer.data(pk_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_y_xyzzzzz = cbuffer.data(pk_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_y_xzzzzzz = cbuffer.data(pk_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_y_yyyyyyy = cbuffer.data(pk_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_y_yyyyyyz = cbuffer.data(pk_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_y_yyyyyzz = cbuffer.data(pk_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_y_yyyyzzz = cbuffer.data(pk_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_y_yyyzzzz = cbuffer.data(pk_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_y_yyzzzzz = cbuffer.data(pk_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_y_yzzzzzz = cbuffer.data(pk_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_y_zzzzzzz = cbuffer.data(pk_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxxx, g_0_x_0_xxxxxxxy, g_0_x_0_xxxxxxy, g_0_x_0_xxxxxxyy, g_0_x_0_xxxxxxyz, g_0_x_0_xxxxxxz, g_0_x_0_xxxxxyy, g_0_x_0_xxxxxyyy, g_0_x_0_xxxxxyyz, g_0_x_0_xxxxxyz, g_0_x_0_xxxxxyzz, g_0_x_0_xxxxxzz, g_0_x_0_xxxxyyy, g_0_x_0_xxxxyyyy, g_0_x_0_xxxxyyyz, g_0_x_0_xxxxyyz, g_0_x_0_xxxxyyzz, g_0_x_0_xxxxyzz, g_0_x_0_xxxxyzzz, g_0_x_0_xxxxzzz, g_0_x_0_xxxyyyy, g_0_x_0_xxxyyyyy, g_0_x_0_xxxyyyyz, g_0_x_0_xxxyyyz, g_0_x_0_xxxyyyzz, g_0_x_0_xxxyyzz, g_0_x_0_xxxyyzzz, g_0_x_0_xxxyzzz, g_0_x_0_xxxyzzzz, g_0_x_0_xxxzzzz, g_0_x_0_xxyyyyy, g_0_x_0_xxyyyyyy, g_0_x_0_xxyyyyyz, g_0_x_0_xxyyyyz, g_0_x_0_xxyyyyzz, g_0_x_0_xxyyyzz, g_0_x_0_xxyyyzzz, g_0_x_0_xxyyzzz, g_0_x_0_xxyyzzzz, g_0_x_0_xxyzzzz, g_0_x_0_xxyzzzzz, g_0_x_0_xxzzzzz, g_0_x_0_xyyyyyy, g_0_x_0_xyyyyyyy, g_0_x_0_xyyyyyyz, g_0_x_0_xyyyyyz, g_0_x_0_xyyyyyzz, g_0_x_0_xyyyyzz, g_0_x_0_xyyyyzzz, g_0_x_0_xyyyzzz, g_0_x_0_xyyyzzzz, g_0_x_0_xyyzzzz, g_0_x_0_xyyzzzzz, g_0_x_0_xyzzzzz, g_0_x_0_xyzzzzzz, g_0_x_0_xzzzzzz, g_0_x_0_yyyyyyy, g_0_x_0_yyyyyyyy, g_0_x_0_yyyyyyyz, g_0_x_0_yyyyyyz, g_0_x_0_yyyyyyzz, g_0_x_0_yyyyyzz, g_0_x_0_yyyyyzzz, g_0_x_0_yyyyzzz, g_0_x_0_yyyyzzzz, g_0_x_0_yyyzzzz, g_0_x_0_yyyzzzzz, g_0_x_0_yyzzzzz, g_0_x_0_yyzzzzzz, g_0_x_0_yzzzzzz, g_0_x_0_yzzzzzzz, g_0_x_0_zzzzzzz, g_0_x_y_xxxxxxx, g_0_x_y_xxxxxxy, g_0_x_y_xxxxxxz, g_0_x_y_xxxxxyy, g_0_x_y_xxxxxyz, g_0_x_y_xxxxxzz, g_0_x_y_xxxxyyy, g_0_x_y_xxxxyyz, g_0_x_y_xxxxyzz, g_0_x_y_xxxxzzz, g_0_x_y_xxxyyyy, g_0_x_y_xxxyyyz, g_0_x_y_xxxyyzz, g_0_x_y_xxxyzzz, g_0_x_y_xxxzzzz, g_0_x_y_xxyyyyy, g_0_x_y_xxyyyyz, g_0_x_y_xxyyyzz, g_0_x_y_xxyyzzz, g_0_x_y_xxyzzzz, g_0_x_y_xxzzzzz, g_0_x_y_xyyyyyy, g_0_x_y_xyyyyyz, g_0_x_y_xyyyyzz, g_0_x_y_xyyyzzz, g_0_x_y_xyyzzzz, g_0_x_y_xyzzzzz, g_0_x_y_xzzzzzz, g_0_x_y_yyyyyyy, g_0_x_y_yyyyyyz, g_0_x_y_yyyyyzz, g_0_x_y_yyyyzzz, g_0_x_y_yyyzzzz, g_0_x_y_yyzzzzz, g_0_x_y_yzzzzzz, g_0_x_y_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_y_xxxxxxx[k] = -g_0_x_0_xxxxxxx[k] * ab_y + g_0_x_0_xxxxxxxy[k];

                g_0_x_y_xxxxxxy[k] = -g_0_x_0_xxxxxxy[k] * ab_y + g_0_x_0_xxxxxxyy[k];

                g_0_x_y_xxxxxxz[k] = -g_0_x_0_xxxxxxz[k] * ab_y + g_0_x_0_xxxxxxyz[k];

                g_0_x_y_xxxxxyy[k] = -g_0_x_0_xxxxxyy[k] * ab_y + g_0_x_0_xxxxxyyy[k];

                g_0_x_y_xxxxxyz[k] = -g_0_x_0_xxxxxyz[k] * ab_y + g_0_x_0_xxxxxyyz[k];

                g_0_x_y_xxxxxzz[k] = -g_0_x_0_xxxxxzz[k] * ab_y + g_0_x_0_xxxxxyzz[k];

                g_0_x_y_xxxxyyy[k] = -g_0_x_0_xxxxyyy[k] * ab_y + g_0_x_0_xxxxyyyy[k];

                g_0_x_y_xxxxyyz[k] = -g_0_x_0_xxxxyyz[k] * ab_y + g_0_x_0_xxxxyyyz[k];

                g_0_x_y_xxxxyzz[k] = -g_0_x_0_xxxxyzz[k] * ab_y + g_0_x_0_xxxxyyzz[k];

                g_0_x_y_xxxxzzz[k] = -g_0_x_0_xxxxzzz[k] * ab_y + g_0_x_0_xxxxyzzz[k];

                g_0_x_y_xxxyyyy[k] = -g_0_x_0_xxxyyyy[k] * ab_y + g_0_x_0_xxxyyyyy[k];

                g_0_x_y_xxxyyyz[k] = -g_0_x_0_xxxyyyz[k] * ab_y + g_0_x_0_xxxyyyyz[k];

                g_0_x_y_xxxyyzz[k] = -g_0_x_0_xxxyyzz[k] * ab_y + g_0_x_0_xxxyyyzz[k];

                g_0_x_y_xxxyzzz[k] = -g_0_x_0_xxxyzzz[k] * ab_y + g_0_x_0_xxxyyzzz[k];

                g_0_x_y_xxxzzzz[k] = -g_0_x_0_xxxzzzz[k] * ab_y + g_0_x_0_xxxyzzzz[k];

                g_0_x_y_xxyyyyy[k] = -g_0_x_0_xxyyyyy[k] * ab_y + g_0_x_0_xxyyyyyy[k];

                g_0_x_y_xxyyyyz[k] = -g_0_x_0_xxyyyyz[k] * ab_y + g_0_x_0_xxyyyyyz[k];

                g_0_x_y_xxyyyzz[k] = -g_0_x_0_xxyyyzz[k] * ab_y + g_0_x_0_xxyyyyzz[k];

                g_0_x_y_xxyyzzz[k] = -g_0_x_0_xxyyzzz[k] * ab_y + g_0_x_0_xxyyyzzz[k];

                g_0_x_y_xxyzzzz[k] = -g_0_x_0_xxyzzzz[k] * ab_y + g_0_x_0_xxyyzzzz[k];

                g_0_x_y_xxzzzzz[k] = -g_0_x_0_xxzzzzz[k] * ab_y + g_0_x_0_xxyzzzzz[k];

                g_0_x_y_xyyyyyy[k] = -g_0_x_0_xyyyyyy[k] * ab_y + g_0_x_0_xyyyyyyy[k];

                g_0_x_y_xyyyyyz[k] = -g_0_x_0_xyyyyyz[k] * ab_y + g_0_x_0_xyyyyyyz[k];

                g_0_x_y_xyyyyzz[k] = -g_0_x_0_xyyyyzz[k] * ab_y + g_0_x_0_xyyyyyzz[k];

                g_0_x_y_xyyyzzz[k] = -g_0_x_0_xyyyzzz[k] * ab_y + g_0_x_0_xyyyyzzz[k];

                g_0_x_y_xyyzzzz[k] = -g_0_x_0_xyyzzzz[k] * ab_y + g_0_x_0_xyyyzzzz[k];

                g_0_x_y_xyzzzzz[k] = -g_0_x_0_xyzzzzz[k] * ab_y + g_0_x_0_xyyzzzzz[k];

                g_0_x_y_xzzzzzz[k] = -g_0_x_0_xzzzzzz[k] * ab_y + g_0_x_0_xyzzzzzz[k];

                g_0_x_y_yyyyyyy[k] = -g_0_x_0_yyyyyyy[k] * ab_y + g_0_x_0_yyyyyyyy[k];

                g_0_x_y_yyyyyyz[k] = -g_0_x_0_yyyyyyz[k] * ab_y + g_0_x_0_yyyyyyyz[k];

                g_0_x_y_yyyyyzz[k] = -g_0_x_0_yyyyyzz[k] * ab_y + g_0_x_0_yyyyyyzz[k];

                g_0_x_y_yyyyzzz[k] = -g_0_x_0_yyyyzzz[k] * ab_y + g_0_x_0_yyyyyzzz[k];

                g_0_x_y_yyyzzzz[k] = -g_0_x_0_yyyzzzz[k] * ab_y + g_0_x_0_yyyyzzzz[k];

                g_0_x_y_yyzzzzz[k] = -g_0_x_0_yyzzzzz[k] * ab_y + g_0_x_0_yyyzzzzz[k];

                g_0_x_y_yzzzzzz[k] = -g_0_x_0_yzzzzzz[k] * ab_y + g_0_x_0_yyzzzzzz[k];

                g_0_x_y_zzzzzzz[k] = -g_0_x_0_zzzzzzz[k] * ab_y + g_0_x_0_yzzzzzzz[k];
            }

            /// Set up 72-108 components of targeted buffer : cbuffer.data(

            auto g_0_x_z_xxxxxxx = cbuffer.data(pk_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxy = cbuffer.data(pk_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_z_xxxxxxz = cbuffer.data(pk_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_z_xxxxxyy = cbuffer.data(pk_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_z_xxxxxyz = cbuffer.data(pk_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_z_xxxxxzz = cbuffer.data(pk_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_z_xxxxyyy = cbuffer.data(pk_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_z_xxxxyyz = cbuffer.data(pk_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_z_xxxxyzz = cbuffer.data(pk_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_z_xxxxzzz = cbuffer.data(pk_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_z_xxxyyyy = cbuffer.data(pk_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_z_xxxyyyz = cbuffer.data(pk_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_z_xxxyyzz = cbuffer.data(pk_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_z_xxxyzzz = cbuffer.data(pk_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_z_xxxzzzz = cbuffer.data(pk_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_z_xxyyyyy = cbuffer.data(pk_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_z_xxyyyyz = cbuffer.data(pk_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_z_xxyyyzz = cbuffer.data(pk_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_z_xxyyzzz = cbuffer.data(pk_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_z_xxyzzzz = cbuffer.data(pk_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_z_xxzzzzz = cbuffer.data(pk_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_z_xyyyyyy = cbuffer.data(pk_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_z_xyyyyyz = cbuffer.data(pk_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_z_xyyyyzz = cbuffer.data(pk_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_z_xyyyzzz = cbuffer.data(pk_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_z_xyyzzzz = cbuffer.data(pk_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_z_xyzzzzz = cbuffer.data(pk_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_z_xzzzzzz = cbuffer.data(pk_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_z_yyyyyyy = cbuffer.data(pk_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_z_yyyyyyz = cbuffer.data(pk_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_z_yyyyyzz = cbuffer.data(pk_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_z_yyyyzzz = cbuffer.data(pk_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_z_yyyzzzz = cbuffer.data(pk_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_z_yyzzzzz = cbuffer.data(pk_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_z_yzzzzzz = cbuffer.data(pk_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_z_zzzzzzz = cbuffer.data(pk_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxxx, g_0_x_0_xxxxxxxz, g_0_x_0_xxxxxxy, g_0_x_0_xxxxxxyz, g_0_x_0_xxxxxxz, g_0_x_0_xxxxxxzz, g_0_x_0_xxxxxyy, g_0_x_0_xxxxxyyz, g_0_x_0_xxxxxyz, g_0_x_0_xxxxxyzz, g_0_x_0_xxxxxzz, g_0_x_0_xxxxxzzz, g_0_x_0_xxxxyyy, g_0_x_0_xxxxyyyz, g_0_x_0_xxxxyyz, g_0_x_0_xxxxyyzz, g_0_x_0_xxxxyzz, g_0_x_0_xxxxyzzz, g_0_x_0_xxxxzzz, g_0_x_0_xxxxzzzz, g_0_x_0_xxxyyyy, g_0_x_0_xxxyyyyz, g_0_x_0_xxxyyyz, g_0_x_0_xxxyyyzz, g_0_x_0_xxxyyzz, g_0_x_0_xxxyyzzz, g_0_x_0_xxxyzzz, g_0_x_0_xxxyzzzz, g_0_x_0_xxxzzzz, g_0_x_0_xxxzzzzz, g_0_x_0_xxyyyyy, g_0_x_0_xxyyyyyz, g_0_x_0_xxyyyyz, g_0_x_0_xxyyyyzz, g_0_x_0_xxyyyzz, g_0_x_0_xxyyyzzz, g_0_x_0_xxyyzzz, g_0_x_0_xxyyzzzz, g_0_x_0_xxyzzzz, g_0_x_0_xxyzzzzz, g_0_x_0_xxzzzzz, g_0_x_0_xxzzzzzz, g_0_x_0_xyyyyyy, g_0_x_0_xyyyyyyz, g_0_x_0_xyyyyyz, g_0_x_0_xyyyyyzz, g_0_x_0_xyyyyzz, g_0_x_0_xyyyyzzz, g_0_x_0_xyyyzzz, g_0_x_0_xyyyzzzz, g_0_x_0_xyyzzzz, g_0_x_0_xyyzzzzz, g_0_x_0_xyzzzzz, g_0_x_0_xyzzzzzz, g_0_x_0_xzzzzzz, g_0_x_0_xzzzzzzz, g_0_x_0_yyyyyyy, g_0_x_0_yyyyyyyz, g_0_x_0_yyyyyyz, g_0_x_0_yyyyyyzz, g_0_x_0_yyyyyzz, g_0_x_0_yyyyyzzz, g_0_x_0_yyyyzzz, g_0_x_0_yyyyzzzz, g_0_x_0_yyyzzzz, g_0_x_0_yyyzzzzz, g_0_x_0_yyzzzzz, g_0_x_0_yyzzzzzz, g_0_x_0_yzzzzzz, g_0_x_0_yzzzzzzz, g_0_x_0_zzzzzzz, g_0_x_0_zzzzzzzz, g_0_x_z_xxxxxxx, g_0_x_z_xxxxxxy, g_0_x_z_xxxxxxz, g_0_x_z_xxxxxyy, g_0_x_z_xxxxxyz, g_0_x_z_xxxxxzz, g_0_x_z_xxxxyyy, g_0_x_z_xxxxyyz, g_0_x_z_xxxxyzz, g_0_x_z_xxxxzzz, g_0_x_z_xxxyyyy, g_0_x_z_xxxyyyz, g_0_x_z_xxxyyzz, g_0_x_z_xxxyzzz, g_0_x_z_xxxzzzz, g_0_x_z_xxyyyyy, g_0_x_z_xxyyyyz, g_0_x_z_xxyyyzz, g_0_x_z_xxyyzzz, g_0_x_z_xxyzzzz, g_0_x_z_xxzzzzz, g_0_x_z_xyyyyyy, g_0_x_z_xyyyyyz, g_0_x_z_xyyyyzz, g_0_x_z_xyyyzzz, g_0_x_z_xyyzzzz, g_0_x_z_xyzzzzz, g_0_x_z_xzzzzzz, g_0_x_z_yyyyyyy, g_0_x_z_yyyyyyz, g_0_x_z_yyyyyzz, g_0_x_z_yyyyzzz, g_0_x_z_yyyzzzz, g_0_x_z_yyzzzzz, g_0_x_z_yzzzzzz, g_0_x_z_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_z_xxxxxxx[k] = -g_0_x_0_xxxxxxx[k] * ab_z + g_0_x_0_xxxxxxxz[k];

                g_0_x_z_xxxxxxy[k] = -g_0_x_0_xxxxxxy[k] * ab_z + g_0_x_0_xxxxxxyz[k];

                g_0_x_z_xxxxxxz[k] = -g_0_x_0_xxxxxxz[k] * ab_z + g_0_x_0_xxxxxxzz[k];

                g_0_x_z_xxxxxyy[k] = -g_0_x_0_xxxxxyy[k] * ab_z + g_0_x_0_xxxxxyyz[k];

                g_0_x_z_xxxxxyz[k] = -g_0_x_0_xxxxxyz[k] * ab_z + g_0_x_0_xxxxxyzz[k];

                g_0_x_z_xxxxxzz[k] = -g_0_x_0_xxxxxzz[k] * ab_z + g_0_x_0_xxxxxzzz[k];

                g_0_x_z_xxxxyyy[k] = -g_0_x_0_xxxxyyy[k] * ab_z + g_0_x_0_xxxxyyyz[k];

                g_0_x_z_xxxxyyz[k] = -g_0_x_0_xxxxyyz[k] * ab_z + g_0_x_0_xxxxyyzz[k];

                g_0_x_z_xxxxyzz[k] = -g_0_x_0_xxxxyzz[k] * ab_z + g_0_x_0_xxxxyzzz[k];

                g_0_x_z_xxxxzzz[k] = -g_0_x_0_xxxxzzz[k] * ab_z + g_0_x_0_xxxxzzzz[k];

                g_0_x_z_xxxyyyy[k] = -g_0_x_0_xxxyyyy[k] * ab_z + g_0_x_0_xxxyyyyz[k];

                g_0_x_z_xxxyyyz[k] = -g_0_x_0_xxxyyyz[k] * ab_z + g_0_x_0_xxxyyyzz[k];

                g_0_x_z_xxxyyzz[k] = -g_0_x_0_xxxyyzz[k] * ab_z + g_0_x_0_xxxyyzzz[k];

                g_0_x_z_xxxyzzz[k] = -g_0_x_0_xxxyzzz[k] * ab_z + g_0_x_0_xxxyzzzz[k];

                g_0_x_z_xxxzzzz[k] = -g_0_x_0_xxxzzzz[k] * ab_z + g_0_x_0_xxxzzzzz[k];

                g_0_x_z_xxyyyyy[k] = -g_0_x_0_xxyyyyy[k] * ab_z + g_0_x_0_xxyyyyyz[k];

                g_0_x_z_xxyyyyz[k] = -g_0_x_0_xxyyyyz[k] * ab_z + g_0_x_0_xxyyyyzz[k];

                g_0_x_z_xxyyyzz[k] = -g_0_x_0_xxyyyzz[k] * ab_z + g_0_x_0_xxyyyzzz[k];

                g_0_x_z_xxyyzzz[k] = -g_0_x_0_xxyyzzz[k] * ab_z + g_0_x_0_xxyyzzzz[k];

                g_0_x_z_xxyzzzz[k] = -g_0_x_0_xxyzzzz[k] * ab_z + g_0_x_0_xxyzzzzz[k];

                g_0_x_z_xxzzzzz[k] = -g_0_x_0_xxzzzzz[k] * ab_z + g_0_x_0_xxzzzzzz[k];

                g_0_x_z_xyyyyyy[k] = -g_0_x_0_xyyyyyy[k] * ab_z + g_0_x_0_xyyyyyyz[k];

                g_0_x_z_xyyyyyz[k] = -g_0_x_0_xyyyyyz[k] * ab_z + g_0_x_0_xyyyyyzz[k];

                g_0_x_z_xyyyyzz[k] = -g_0_x_0_xyyyyzz[k] * ab_z + g_0_x_0_xyyyyzzz[k];

                g_0_x_z_xyyyzzz[k] = -g_0_x_0_xyyyzzz[k] * ab_z + g_0_x_0_xyyyzzzz[k];

                g_0_x_z_xyyzzzz[k] = -g_0_x_0_xyyzzzz[k] * ab_z + g_0_x_0_xyyzzzzz[k];

                g_0_x_z_xyzzzzz[k] = -g_0_x_0_xyzzzzz[k] * ab_z + g_0_x_0_xyzzzzzz[k];

                g_0_x_z_xzzzzzz[k] = -g_0_x_0_xzzzzzz[k] * ab_z + g_0_x_0_xzzzzzzz[k];

                g_0_x_z_yyyyyyy[k] = -g_0_x_0_yyyyyyy[k] * ab_z + g_0_x_0_yyyyyyyz[k];

                g_0_x_z_yyyyyyz[k] = -g_0_x_0_yyyyyyz[k] * ab_z + g_0_x_0_yyyyyyzz[k];

                g_0_x_z_yyyyyzz[k] = -g_0_x_0_yyyyyzz[k] * ab_z + g_0_x_0_yyyyyzzz[k];

                g_0_x_z_yyyyzzz[k] = -g_0_x_0_yyyyzzz[k] * ab_z + g_0_x_0_yyyyzzzz[k];

                g_0_x_z_yyyzzzz[k] = -g_0_x_0_yyyzzzz[k] * ab_z + g_0_x_0_yyyzzzzz[k];

                g_0_x_z_yyzzzzz[k] = -g_0_x_0_yyzzzzz[k] * ab_z + g_0_x_0_yyzzzzzz[k];

                g_0_x_z_yzzzzzz[k] = -g_0_x_0_yzzzzzz[k] * ab_z + g_0_x_0_yzzzzzzz[k];

                g_0_x_z_zzzzzzz[k] = -g_0_x_0_zzzzzzz[k] * ab_z + g_0_x_0_zzzzzzzz[k];
            }

            /// Set up 108-144 components of targeted buffer : cbuffer.data(

            auto g_0_y_x_xxxxxxx = cbuffer.data(pk_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxy = cbuffer.data(pk_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_x_xxxxxxz = cbuffer.data(pk_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_x_xxxxxyy = cbuffer.data(pk_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_x_xxxxxyz = cbuffer.data(pk_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_x_xxxxxzz = cbuffer.data(pk_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_x_xxxxyyy = cbuffer.data(pk_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_x_xxxxyyz = cbuffer.data(pk_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_x_xxxxyzz = cbuffer.data(pk_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_x_xxxxzzz = cbuffer.data(pk_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_x_xxxyyyy = cbuffer.data(pk_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_x_xxxyyyz = cbuffer.data(pk_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_x_xxxyyzz = cbuffer.data(pk_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_x_xxxyzzz = cbuffer.data(pk_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_x_xxxzzzz = cbuffer.data(pk_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_x_xxyyyyy = cbuffer.data(pk_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_x_xxyyyyz = cbuffer.data(pk_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_x_xxyyyzz = cbuffer.data(pk_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_x_xxyyzzz = cbuffer.data(pk_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_x_xxyzzzz = cbuffer.data(pk_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_x_xxzzzzz = cbuffer.data(pk_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_x_xyyyyyy = cbuffer.data(pk_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_x_xyyyyyz = cbuffer.data(pk_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_x_xyyyyzz = cbuffer.data(pk_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_x_xyyyzzz = cbuffer.data(pk_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_x_xyyzzzz = cbuffer.data(pk_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_x_xyzzzzz = cbuffer.data(pk_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_x_xzzzzzz = cbuffer.data(pk_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_x_yyyyyyy = cbuffer.data(pk_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_x_yyyyyyz = cbuffer.data(pk_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_x_yyyyyzz = cbuffer.data(pk_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_x_yyyyzzz = cbuffer.data(pk_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_x_yyyzzzz = cbuffer.data(pk_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_x_yyzzzzz = cbuffer.data(pk_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_x_yzzzzzz = cbuffer.data(pk_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_x_zzzzzzz = cbuffer.data(pk_geom_01_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxxxx, g_0_y_0_xxxxxxxx, g_0_y_0_xxxxxxxy, g_0_y_0_xxxxxxxz, g_0_y_0_xxxxxxy, g_0_y_0_xxxxxxyy, g_0_y_0_xxxxxxyz, g_0_y_0_xxxxxxz, g_0_y_0_xxxxxxzz, g_0_y_0_xxxxxyy, g_0_y_0_xxxxxyyy, g_0_y_0_xxxxxyyz, g_0_y_0_xxxxxyz, g_0_y_0_xxxxxyzz, g_0_y_0_xxxxxzz, g_0_y_0_xxxxxzzz, g_0_y_0_xxxxyyy, g_0_y_0_xxxxyyyy, g_0_y_0_xxxxyyyz, g_0_y_0_xxxxyyz, g_0_y_0_xxxxyyzz, g_0_y_0_xxxxyzz, g_0_y_0_xxxxyzzz, g_0_y_0_xxxxzzz, g_0_y_0_xxxxzzzz, g_0_y_0_xxxyyyy, g_0_y_0_xxxyyyyy, g_0_y_0_xxxyyyyz, g_0_y_0_xxxyyyz, g_0_y_0_xxxyyyzz, g_0_y_0_xxxyyzz, g_0_y_0_xxxyyzzz, g_0_y_0_xxxyzzz, g_0_y_0_xxxyzzzz, g_0_y_0_xxxzzzz, g_0_y_0_xxxzzzzz, g_0_y_0_xxyyyyy, g_0_y_0_xxyyyyyy, g_0_y_0_xxyyyyyz, g_0_y_0_xxyyyyz, g_0_y_0_xxyyyyzz, g_0_y_0_xxyyyzz, g_0_y_0_xxyyyzzz, g_0_y_0_xxyyzzz, g_0_y_0_xxyyzzzz, g_0_y_0_xxyzzzz, g_0_y_0_xxyzzzzz, g_0_y_0_xxzzzzz, g_0_y_0_xxzzzzzz, g_0_y_0_xyyyyyy, g_0_y_0_xyyyyyyy, g_0_y_0_xyyyyyyz, g_0_y_0_xyyyyyz, g_0_y_0_xyyyyyzz, g_0_y_0_xyyyyzz, g_0_y_0_xyyyyzzz, g_0_y_0_xyyyzzz, g_0_y_0_xyyyzzzz, g_0_y_0_xyyzzzz, g_0_y_0_xyyzzzzz, g_0_y_0_xyzzzzz, g_0_y_0_xyzzzzzz, g_0_y_0_xzzzzzz, g_0_y_0_xzzzzzzz, g_0_y_0_yyyyyyy, g_0_y_0_yyyyyyz, g_0_y_0_yyyyyzz, g_0_y_0_yyyyzzz, g_0_y_0_yyyzzzz, g_0_y_0_yyzzzzz, g_0_y_0_yzzzzzz, g_0_y_0_zzzzzzz, g_0_y_x_xxxxxxx, g_0_y_x_xxxxxxy, g_0_y_x_xxxxxxz, g_0_y_x_xxxxxyy, g_0_y_x_xxxxxyz, g_0_y_x_xxxxxzz, g_0_y_x_xxxxyyy, g_0_y_x_xxxxyyz, g_0_y_x_xxxxyzz, g_0_y_x_xxxxzzz, g_0_y_x_xxxyyyy, g_0_y_x_xxxyyyz, g_0_y_x_xxxyyzz, g_0_y_x_xxxyzzz, g_0_y_x_xxxzzzz, g_0_y_x_xxyyyyy, g_0_y_x_xxyyyyz, g_0_y_x_xxyyyzz, g_0_y_x_xxyyzzz, g_0_y_x_xxyzzzz, g_0_y_x_xxzzzzz, g_0_y_x_xyyyyyy, g_0_y_x_xyyyyyz, g_0_y_x_xyyyyzz, g_0_y_x_xyyyzzz, g_0_y_x_xyyzzzz, g_0_y_x_xyzzzzz, g_0_y_x_xzzzzzz, g_0_y_x_yyyyyyy, g_0_y_x_yyyyyyz, g_0_y_x_yyyyyzz, g_0_y_x_yyyyzzz, g_0_y_x_yyyzzzz, g_0_y_x_yyzzzzz, g_0_y_x_yzzzzzz, g_0_y_x_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_x_xxxxxxx[k] = -g_0_y_0_xxxxxxx[k] * ab_x + g_0_y_0_xxxxxxxx[k];

                g_0_y_x_xxxxxxy[k] = -g_0_y_0_xxxxxxy[k] * ab_x + g_0_y_0_xxxxxxxy[k];

                g_0_y_x_xxxxxxz[k] = -g_0_y_0_xxxxxxz[k] * ab_x + g_0_y_0_xxxxxxxz[k];

                g_0_y_x_xxxxxyy[k] = -g_0_y_0_xxxxxyy[k] * ab_x + g_0_y_0_xxxxxxyy[k];

                g_0_y_x_xxxxxyz[k] = -g_0_y_0_xxxxxyz[k] * ab_x + g_0_y_0_xxxxxxyz[k];

                g_0_y_x_xxxxxzz[k] = -g_0_y_0_xxxxxzz[k] * ab_x + g_0_y_0_xxxxxxzz[k];

                g_0_y_x_xxxxyyy[k] = -g_0_y_0_xxxxyyy[k] * ab_x + g_0_y_0_xxxxxyyy[k];

                g_0_y_x_xxxxyyz[k] = -g_0_y_0_xxxxyyz[k] * ab_x + g_0_y_0_xxxxxyyz[k];

                g_0_y_x_xxxxyzz[k] = -g_0_y_0_xxxxyzz[k] * ab_x + g_0_y_0_xxxxxyzz[k];

                g_0_y_x_xxxxzzz[k] = -g_0_y_0_xxxxzzz[k] * ab_x + g_0_y_0_xxxxxzzz[k];

                g_0_y_x_xxxyyyy[k] = -g_0_y_0_xxxyyyy[k] * ab_x + g_0_y_0_xxxxyyyy[k];

                g_0_y_x_xxxyyyz[k] = -g_0_y_0_xxxyyyz[k] * ab_x + g_0_y_0_xxxxyyyz[k];

                g_0_y_x_xxxyyzz[k] = -g_0_y_0_xxxyyzz[k] * ab_x + g_0_y_0_xxxxyyzz[k];

                g_0_y_x_xxxyzzz[k] = -g_0_y_0_xxxyzzz[k] * ab_x + g_0_y_0_xxxxyzzz[k];

                g_0_y_x_xxxzzzz[k] = -g_0_y_0_xxxzzzz[k] * ab_x + g_0_y_0_xxxxzzzz[k];

                g_0_y_x_xxyyyyy[k] = -g_0_y_0_xxyyyyy[k] * ab_x + g_0_y_0_xxxyyyyy[k];

                g_0_y_x_xxyyyyz[k] = -g_0_y_0_xxyyyyz[k] * ab_x + g_0_y_0_xxxyyyyz[k];

                g_0_y_x_xxyyyzz[k] = -g_0_y_0_xxyyyzz[k] * ab_x + g_0_y_0_xxxyyyzz[k];

                g_0_y_x_xxyyzzz[k] = -g_0_y_0_xxyyzzz[k] * ab_x + g_0_y_0_xxxyyzzz[k];

                g_0_y_x_xxyzzzz[k] = -g_0_y_0_xxyzzzz[k] * ab_x + g_0_y_0_xxxyzzzz[k];

                g_0_y_x_xxzzzzz[k] = -g_0_y_0_xxzzzzz[k] * ab_x + g_0_y_0_xxxzzzzz[k];

                g_0_y_x_xyyyyyy[k] = -g_0_y_0_xyyyyyy[k] * ab_x + g_0_y_0_xxyyyyyy[k];

                g_0_y_x_xyyyyyz[k] = -g_0_y_0_xyyyyyz[k] * ab_x + g_0_y_0_xxyyyyyz[k];

                g_0_y_x_xyyyyzz[k] = -g_0_y_0_xyyyyzz[k] * ab_x + g_0_y_0_xxyyyyzz[k];

                g_0_y_x_xyyyzzz[k] = -g_0_y_0_xyyyzzz[k] * ab_x + g_0_y_0_xxyyyzzz[k];

                g_0_y_x_xyyzzzz[k] = -g_0_y_0_xyyzzzz[k] * ab_x + g_0_y_0_xxyyzzzz[k];

                g_0_y_x_xyzzzzz[k] = -g_0_y_0_xyzzzzz[k] * ab_x + g_0_y_0_xxyzzzzz[k];

                g_0_y_x_xzzzzzz[k] = -g_0_y_0_xzzzzzz[k] * ab_x + g_0_y_0_xxzzzzzz[k];

                g_0_y_x_yyyyyyy[k] = -g_0_y_0_yyyyyyy[k] * ab_x + g_0_y_0_xyyyyyyy[k];

                g_0_y_x_yyyyyyz[k] = -g_0_y_0_yyyyyyz[k] * ab_x + g_0_y_0_xyyyyyyz[k];

                g_0_y_x_yyyyyzz[k] = -g_0_y_0_yyyyyzz[k] * ab_x + g_0_y_0_xyyyyyzz[k];

                g_0_y_x_yyyyzzz[k] = -g_0_y_0_yyyyzzz[k] * ab_x + g_0_y_0_xyyyyzzz[k];

                g_0_y_x_yyyzzzz[k] = -g_0_y_0_yyyzzzz[k] * ab_x + g_0_y_0_xyyyzzzz[k];

                g_0_y_x_yyzzzzz[k] = -g_0_y_0_yyzzzzz[k] * ab_x + g_0_y_0_xyyzzzzz[k];

                g_0_y_x_yzzzzzz[k] = -g_0_y_0_yzzzzzz[k] * ab_x + g_0_y_0_xyzzzzzz[k];

                g_0_y_x_zzzzzzz[k] = -g_0_y_0_zzzzzzz[k] * ab_x + g_0_y_0_xzzzzzzz[k];
            }

            /// Set up 144-180 components of targeted buffer : cbuffer.data(

            auto g_0_y_y_xxxxxxx = cbuffer.data(pk_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxy = cbuffer.data(pk_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_y_xxxxxxz = cbuffer.data(pk_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_y_xxxxxyy = cbuffer.data(pk_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_y_xxxxxyz = cbuffer.data(pk_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_y_xxxxxzz = cbuffer.data(pk_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_y_xxxxyyy = cbuffer.data(pk_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_y_xxxxyyz = cbuffer.data(pk_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_y_xxxxyzz = cbuffer.data(pk_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_y_xxxxzzz = cbuffer.data(pk_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_y_xxxyyyy = cbuffer.data(pk_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_y_xxxyyyz = cbuffer.data(pk_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_y_xxxyyzz = cbuffer.data(pk_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_y_xxxyzzz = cbuffer.data(pk_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_y_xxxzzzz = cbuffer.data(pk_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_y_xxyyyyy = cbuffer.data(pk_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_y_xxyyyyz = cbuffer.data(pk_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_y_xxyyyzz = cbuffer.data(pk_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_y_xxyyzzz = cbuffer.data(pk_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_y_xxyzzzz = cbuffer.data(pk_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_y_xxzzzzz = cbuffer.data(pk_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_y_xyyyyyy = cbuffer.data(pk_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_y_xyyyyyz = cbuffer.data(pk_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_y_xyyyyzz = cbuffer.data(pk_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_y_xyyyzzz = cbuffer.data(pk_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_y_xyyzzzz = cbuffer.data(pk_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_y_xyzzzzz = cbuffer.data(pk_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_y_xzzzzzz = cbuffer.data(pk_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_y_yyyyyyy = cbuffer.data(pk_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_y_yyyyyyz = cbuffer.data(pk_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_y_yyyyyzz = cbuffer.data(pk_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_y_yyyyzzz = cbuffer.data(pk_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_y_yyyzzzz = cbuffer.data(pk_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_y_yyzzzzz = cbuffer.data(pk_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_y_yzzzzzz = cbuffer.data(pk_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_y_zzzzzzz = cbuffer.data(pk_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxxxx, g_0_xxxxxxy, g_0_xxxxxxz, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxzz, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyzz, g_0_xxxxzzz, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyzz, g_0_xxxyzzz, g_0_xxxzzzz, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyzz, g_0_xxyyzzz, g_0_xxyzzzz, g_0_xxzzzzz, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyzz, g_0_xyyyzzz, g_0_xyyzzzz, g_0_xyzzzzz, g_0_xzzzzzz, g_0_y_0_xxxxxxx, g_0_y_0_xxxxxxxy, g_0_y_0_xxxxxxy, g_0_y_0_xxxxxxyy, g_0_y_0_xxxxxxyz, g_0_y_0_xxxxxxz, g_0_y_0_xxxxxyy, g_0_y_0_xxxxxyyy, g_0_y_0_xxxxxyyz, g_0_y_0_xxxxxyz, g_0_y_0_xxxxxyzz, g_0_y_0_xxxxxzz, g_0_y_0_xxxxyyy, g_0_y_0_xxxxyyyy, g_0_y_0_xxxxyyyz, g_0_y_0_xxxxyyz, g_0_y_0_xxxxyyzz, g_0_y_0_xxxxyzz, g_0_y_0_xxxxyzzz, g_0_y_0_xxxxzzz, g_0_y_0_xxxyyyy, g_0_y_0_xxxyyyyy, g_0_y_0_xxxyyyyz, g_0_y_0_xxxyyyz, g_0_y_0_xxxyyyzz, g_0_y_0_xxxyyzz, g_0_y_0_xxxyyzzz, g_0_y_0_xxxyzzz, g_0_y_0_xxxyzzzz, g_0_y_0_xxxzzzz, g_0_y_0_xxyyyyy, g_0_y_0_xxyyyyyy, g_0_y_0_xxyyyyyz, g_0_y_0_xxyyyyz, g_0_y_0_xxyyyyzz, g_0_y_0_xxyyyzz, g_0_y_0_xxyyyzzz, g_0_y_0_xxyyzzz, g_0_y_0_xxyyzzzz, g_0_y_0_xxyzzzz, g_0_y_0_xxyzzzzz, g_0_y_0_xxzzzzz, g_0_y_0_xyyyyyy, g_0_y_0_xyyyyyyy, g_0_y_0_xyyyyyyz, g_0_y_0_xyyyyyz, g_0_y_0_xyyyyyzz, g_0_y_0_xyyyyzz, g_0_y_0_xyyyyzzz, g_0_y_0_xyyyzzz, g_0_y_0_xyyyzzzz, g_0_y_0_xyyzzzz, g_0_y_0_xyyzzzzz, g_0_y_0_xyzzzzz, g_0_y_0_xyzzzzzz, g_0_y_0_xzzzzzz, g_0_y_0_yyyyyyy, g_0_y_0_yyyyyyyy, g_0_y_0_yyyyyyyz, g_0_y_0_yyyyyyz, g_0_y_0_yyyyyyzz, g_0_y_0_yyyyyzz, g_0_y_0_yyyyyzzz, g_0_y_0_yyyyzzz, g_0_y_0_yyyyzzzz, g_0_y_0_yyyzzzz, g_0_y_0_yyyzzzzz, g_0_y_0_yyzzzzz, g_0_y_0_yyzzzzzz, g_0_y_0_yzzzzzz, g_0_y_0_yzzzzzzz, g_0_y_0_zzzzzzz, g_0_y_y_xxxxxxx, g_0_y_y_xxxxxxy, g_0_y_y_xxxxxxz, g_0_y_y_xxxxxyy, g_0_y_y_xxxxxyz, g_0_y_y_xxxxxzz, g_0_y_y_xxxxyyy, g_0_y_y_xxxxyyz, g_0_y_y_xxxxyzz, g_0_y_y_xxxxzzz, g_0_y_y_xxxyyyy, g_0_y_y_xxxyyyz, g_0_y_y_xxxyyzz, g_0_y_y_xxxyzzz, g_0_y_y_xxxzzzz, g_0_y_y_xxyyyyy, g_0_y_y_xxyyyyz, g_0_y_y_xxyyyzz, g_0_y_y_xxyyzzz, g_0_y_y_xxyzzzz, g_0_y_y_xxzzzzz, g_0_y_y_xyyyyyy, g_0_y_y_xyyyyyz, g_0_y_y_xyyyyzz, g_0_y_y_xyyyzzz, g_0_y_y_xyyzzzz, g_0_y_y_xyzzzzz, g_0_y_y_xzzzzzz, g_0_y_y_yyyyyyy, g_0_y_y_yyyyyyz, g_0_y_y_yyyyyzz, g_0_y_y_yyyyzzz, g_0_y_y_yyyzzzz, g_0_y_y_yyzzzzz, g_0_y_y_yzzzzzz, g_0_y_y_zzzzzzz, g_0_yyyyyyy, g_0_yyyyyyz, g_0_yyyyyzz, g_0_yyyyzzz, g_0_yyyzzzz, g_0_yyzzzzz, g_0_yzzzzzz, g_0_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_y_xxxxxxx[k] = g_0_xxxxxxx[k] - g_0_y_0_xxxxxxx[k] * ab_y + g_0_y_0_xxxxxxxy[k];

                g_0_y_y_xxxxxxy[k] = g_0_xxxxxxy[k] - g_0_y_0_xxxxxxy[k] * ab_y + g_0_y_0_xxxxxxyy[k];

                g_0_y_y_xxxxxxz[k] = g_0_xxxxxxz[k] - g_0_y_0_xxxxxxz[k] * ab_y + g_0_y_0_xxxxxxyz[k];

                g_0_y_y_xxxxxyy[k] = g_0_xxxxxyy[k] - g_0_y_0_xxxxxyy[k] * ab_y + g_0_y_0_xxxxxyyy[k];

                g_0_y_y_xxxxxyz[k] = g_0_xxxxxyz[k] - g_0_y_0_xxxxxyz[k] * ab_y + g_0_y_0_xxxxxyyz[k];

                g_0_y_y_xxxxxzz[k] = g_0_xxxxxzz[k] - g_0_y_0_xxxxxzz[k] * ab_y + g_0_y_0_xxxxxyzz[k];

                g_0_y_y_xxxxyyy[k] = g_0_xxxxyyy[k] - g_0_y_0_xxxxyyy[k] * ab_y + g_0_y_0_xxxxyyyy[k];

                g_0_y_y_xxxxyyz[k] = g_0_xxxxyyz[k] - g_0_y_0_xxxxyyz[k] * ab_y + g_0_y_0_xxxxyyyz[k];

                g_0_y_y_xxxxyzz[k] = g_0_xxxxyzz[k] - g_0_y_0_xxxxyzz[k] * ab_y + g_0_y_0_xxxxyyzz[k];

                g_0_y_y_xxxxzzz[k] = g_0_xxxxzzz[k] - g_0_y_0_xxxxzzz[k] * ab_y + g_0_y_0_xxxxyzzz[k];

                g_0_y_y_xxxyyyy[k] = g_0_xxxyyyy[k] - g_0_y_0_xxxyyyy[k] * ab_y + g_0_y_0_xxxyyyyy[k];

                g_0_y_y_xxxyyyz[k] = g_0_xxxyyyz[k] - g_0_y_0_xxxyyyz[k] * ab_y + g_0_y_0_xxxyyyyz[k];

                g_0_y_y_xxxyyzz[k] = g_0_xxxyyzz[k] - g_0_y_0_xxxyyzz[k] * ab_y + g_0_y_0_xxxyyyzz[k];

                g_0_y_y_xxxyzzz[k] = g_0_xxxyzzz[k] - g_0_y_0_xxxyzzz[k] * ab_y + g_0_y_0_xxxyyzzz[k];

                g_0_y_y_xxxzzzz[k] = g_0_xxxzzzz[k] - g_0_y_0_xxxzzzz[k] * ab_y + g_0_y_0_xxxyzzzz[k];

                g_0_y_y_xxyyyyy[k] = g_0_xxyyyyy[k] - g_0_y_0_xxyyyyy[k] * ab_y + g_0_y_0_xxyyyyyy[k];

                g_0_y_y_xxyyyyz[k] = g_0_xxyyyyz[k] - g_0_y_0_xxyyyyz[k] * ab_y + g_0_y_0_xxyyyyyz[k];

                g_0_y_y_xxyyyzz[k] = g_0_xxyyyzz[k] - g_0_y_0_xxyyyzz[k] * ab_y + g_0_y_0_xxyyyyzz[k];

                g_0_y_y_xxyyzzz[k] = g_0_xxyyzzz[k] - g_0_y_0_xxyyzzz[k] * ab_y + g_0_y_0_xxyyyzzz[k];

                g_0_y_y_xxyzzzz[k] = g_0_xxyzzzz[k] - g_0_y_0_xxyzzzz[k] * ab_y + g_0_y_0_xxyyzzzz[k];

                g_0_y_y_xxzzzzz[k] = g_0_xxzzzzz[k] - g_0_y_0_xxzzzzz[k] * ab_y + g_0_y_0_xxyzzzzz[k];

                g_0_y_y_xyyyyyy[k] = g_0_xyyyyyy[k] - g_0_y_0_xyyyyyy[k] * ab_y + g_0_y_0_xyyyyyyy[k];

                g_0_y_y_xyyyyyz[k] = g_0_xyyyyyz[k] - g_0_y_0_xyyyyyz[k] * ab_y + g_0_y_0_xyyyyyyz[k];

                g_0_y_y_xyyyyzz[k] = g_0_xyyyyzz[k] - g_0_y_0_xyyyyzz[k] * ab_y + g_0_y_0_xyyyyyzz[k];

                g_0_y_y_xyyyzzz[k] = g_0_xyyyzzz[k] - g_0_y_0_xyyyzzz[k] * ab_y + g_0_y_0_xyyyyzzz[k];

                g_0_y_y_xyyzzzz[k] = g_0_xyyzzzz[k] - g_0_y_0_xyyzzzz[k] * ab_y + g_0_y_0_xyyyzzzz[k];

                g_0_y_y_xyzzzzz[k] = g_0_xyzzzzz[k] - g_0_y_0_xyzzzzz[k] * ab_y + g_0_y_0_xyyzzzzz[k];

                g_0_y_y_xzzzzzz[k] = g_0_xzzzzzz[k] - g_0_y_0_xzzzzzz[k] * ab_y + g_0_y_0_xyzzzzzz[k];

                g_0_y_y_yyyyyyy[k] = g_0_yyyyyyy[k] - g_0_y_0_yyyyyyy[k] * ab_y + g_0_y_0_yyyyyyyy[k];

                g_0_y_y_yyyyyyz[k] = g_0_yyyyyyz[k] - g_0_y_0_yyyyyyz[k] * ab_y + g_0_y_0_yyyyyyyz[k];

                g_0_y_y_yyyyyzz[k] = g_0_yyyyyzz[k] - g_0_y_0_yyyyyzz[k] * ab_y + g_0_y_0_yyyyyyzz[k];

                g_0_y_y_yyyyzzz[k] = g_0_yyyyzzz[k] - g_0_y_0_yyyyzzz[k] * ab_y + g_0_y_0_yyyyyzzz[k];

                g_0_y_y_yyyzzzz[k] = g_0_yyyzzzz[k] - g_0_y_0_yyyzzzz[k] * ab_y + g_0_y_0_yyyyzzzz[k];

                g_0_y_y_yyzzzzz[k] = g_0_yyzzzzz[k] - g_0_y_0_yyzzzzz[k] * ab_y + g_0_y_0_yyyzzzzz[k];

                g_0_y_y_yzzzzzz[k] = g_0_yzzzzzz[k] - g_0_y_0_yzzzzzz[k] * ab_y + g_0_y_0_yyzzzzzz[k];

                g_0_y_y_zzzzzzz[k] = g_0_zzzzzzz[k] - g_0_y_0_zzzzzzz[k] * ab_y + g_0_y_0_yzzzzzzz[k];
            }

            /// Set up 180-216 components of targeted buffer : cbuffer.data(

            auto g_0_y_z_xxxxxxx = cbuffer.data(pk_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxy = cbuffer.data(pk_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_z_xxxxxxz = cbuffer.data(pk_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_z_xxxxxyy = cbuffer.data(pk_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_z_xxxxxyz = cbuffer.data(pk_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_z_xxxxxzz = cbuffer.data(pk_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_z_xxxxyyy = cbuffer.data(pk_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_z_xxxxyyz = cbuffer.data(pk_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_z_xxxxyzz = cbuffer.data(pk_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_z_xxxxzzz = cbuffer.data(pk_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_z_xxxyyyy = cbuffer.data(pk_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_z_xxxyyyz = cbuffer.data(pk_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_z_xxxyyzz = cbuffer.data(pk_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_z_xxxyzzz = cbuffer.data(pk_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_z_xxxzzzz = cbuffer.data(pk_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_z_xxyyyyy = cbuffer.data(pk_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_z_xxyyyyz = cbuffer.data(pk_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_z_xxyyyzz = cbuffer.data(pk_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_z_xxyyzzz = cbuffer.data(pk_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_z_xxyzzzz = cbuffer.data(pk_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_z_xxzzzzz = cbuffer.data(pk_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_z_xyyyyyy = cbuffer.data(pk_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_z_xyyyyyz = cbuffer.data(pk_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_z_xyyyyzz = cbuffer.data(pk_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_z_xyyyzzz = cbuffer.data(pk_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_z_xyyzzzz = cbuffer.data(pk_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_z_xyzzzzz = cbuffer.data(pk_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_z_xzzzzzz = cbuffer.data(pk_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_z_yyyyyyy = cbuffer.data(pk_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_z_yyyyyyz = cbuffer.data(pk_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_z_yyyyyzz = cbuffer.data(pk_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_z_yyyyzzz = cbuffer.data(pk_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_z_yyyzzzz = cbuffer.data(pk_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_z_yyzzzzz = cbuffer.data(pk_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_z_yzzzzzz = cbuffer.data(pk_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_z_zzzzzzz = cbuffer.data(pk_geom_01_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxxxx, g_0_y_0_xxxxxxxz, g_0_y_0_xxxxxxy, g_0_y_0_xxxxxxyz, g_0_y_0_xxxxxxz, g_0_y_0_xxxxxxzz, g_0_y_0_xxxxxyy, g_0_y_0_xxxxxyyz, g_0_y_0_xxxxxyz, g_0_y_0_xxxxxyzz, g_0_y_0_xxxxxzz, g_0_y_0_xxxxxzzz, g_0_y_0_xxxxyyy, g_0_y_0_xxxxyyyz, g_0_y_0_xxxxyyz, g_0_y_0_xxxxyyzz, g_0_y_0_xxxxyzz, g_0_y_0_xxxxyzzz, g_0_y_0_xxxxzzz, g_0_y_0_xxxxzzzz, g_0_y_0_xxxyyyy, g_0_y_0_xxxyyyyz, g_0_y_0_xxxyyyz, g_0_y_0_xxxyyyzz, g_0_y_0_xxxyyzz, g_0_y_0_xxxyyzzz, g_0_y_0_xxxyzzz, g_0_y_0_xxxyzzzz, g_0_y_0_xxxzzzz, g_0_y_0_xxxzzzzz, g_0_y_0_xxyyyyy, g_0_y_0_xxyyyyyz, g_0_y_0_xxyyyyz, g_0_y_0_xxyyyyzz, g_0_y_0_xxyyyzz, g_0_y_0_xxyyyzzz, g_0_y_0_xxyyzzz, g_0_y_0_xxyyzzzz, g_0_y_0_xxyzzzz, g_0_y_0_xxyzzzzz, g_0_y_0_xxzzzzz, g_0_y_0_xxzzzzzz, g_0_y_0_xyyyyyy, g_0_y_0_xyyyyyyz, g_0_y_0_xyyyyyz, g_0_y_0_xyyyyyzz, g_0_y_0_xyyyyzz, g_0_y_0_xyyyyzzz, g_0_y_0_xyyyzzz, g_0_y_0_xyyyzzzz, g_0_y_0_xyyzzzz, g_0_y_0_xyyzzzzz, g_0_y_0_xyzzzzz, g_0_y_0_xyzzzzzz, g_0_y_0_xzzzzzz, g_0_y_0_xzzzzzzz, g_0_y_0_yyyyyyy, g_0_y_0_yyyyyyyz, g_0_y_0_yyyyyyz, g_0_y_0_yyyyyyzz, g_0_y_0_yyyyyzz, g_0_y_0_yyyyyzzz, g_0_y_0_yyyyzzz, g_0_y_0_yyyyzzzz, g_0_y_0_yyyzzzz, g_0_y_0_yyyzzzzz, g_0_y_0_yyzzzzz, g_0_y_0_yyzzzzzz, g_0_y_0_yzzzzzz, g_0_y_0_yzzzzzzz, g_0_y_0_zzzzzzz, g_0_y_0_zzzzzzzz, g_0_y_z_xxxxxxx, g_0_y_z_xxxxxxy, g_0_y_z_xxxxxxz, g_0_y_z_xxxxxyy, g_0_y_z_xxxxxyz, g_0_y_z_xxxxxzz, g_0_y_z_xxxxyyy, g_0_y_z_xxxxyyz, g_0_y_z_xxxxyzz, g_0_y_z_xxxxzzz, g_0_y_z_xxxyyyy, g_0_y_z_xxxyyyz, g_0_y_z_xxxyyzz, g_0_y_z_xxxyzzz, g_0_y_z_xxxzzzz, g_0_y_z_xxyyyyy, g_0_y_z_xxyyyyz, g_0_y_z_xxyyyzz, g_0_y_z_xxyyzzz, g_0_y_z_xxyzzzz, g_0_y_z_xxzzzzz, g_0_y_z_xyyyyyy, g_0_y_z_xyyyyyz, g_0_y_z_xyyyyzz, g_0_y_z_xyyyzzz, g_0_y_z_xyyzzzz, g_0_y_z_xyzzzzz, g_0_y_z_xzzzzzz, g_0_y_z_yyyyyyy, g_0_y_z_yyyyyyz, g_0_y_z_yyyyyzz, g_0_y_z_yyyyzzz, g_0_y_z_yyyzzzz, g_0_y_z_yyzzzzz, g_0_y_z_yzzzzzz, g_0_y_z_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_z_xxxxxxx[k] = -g_0_y_0_xxxxxxx[k] * ab_z + g_0_y_0_xxxxxxxz[k];

                g_0_y_z_xxxxxxy[k] = -g_0_y_0_xxxxxxy[k] * ab_z + g_0_y_0_xxxxxxyz[k];

                g_0_y_z_xxxxxxz[k] = -g_0_y_0_xxxxxxz[k] * ab_z + g_0_y_0_xxxxxxzz[k];

                g_0_y_z_xxxxxyy[k] = -g_0_y_0_xxxxxyy[k] * ab_z + g_0_y_0_xxxxxyyz[k];

                g_0_y_z_xxxxxyz[k] = -g_0_y_0_xxxxxyz[k] * ab_z + g_0_y_0_xxxxxyzz[k];

                g_0_y_z_xxxxxzz[k] = -g_0_y_0_xxxxxzz[k] * ab_z + g_0_y_0_xxxxxzzz[k];

                g_0_y_z_xxxxyyy[k] = -g_0_y_0_xxxxyyy[k] * ab_z + g_0_y_0_xxxxyyyz[k];

                g_0_y_z_xxxxyyz[k] = -g_0_y_0_xxxxyyz[k] * ab_z + g_0_y_0_xxxxyyzz[k];

                g_0_y_z_xxxxyzz[k] = -g_0_y_0_xxxxyzz[k] * ab_z + g_0_y_0_xxxxyzzz[k];

                g_0_y_z_xxxxzzz[k] = -g_0_y_0_xxxxzzz[k] * ab_z + g_0_y_0_xxxxzzzz[k];

                g_0_y_z_xxxyyyy[k] = -g_0_y_0_xxxyyyy[k] * ab_z + g_0_y_0_xxxyyyyz[k];

                g_0_y_z_xxxyyyz[k] = -g_0_y_0_xxxyyyz[k] * ab_z + g_0_y_0_xxxyyyzz[k];

                g_0_y_z_xxxyyzz[k] = -g_0_y_0_xxxyyzz[k] * ab_z + g_0_y_0_xxxyyzzz[k];

                g_0_y_z_xxxyzzz[k] = -g_0_y_0_xxxyzzz[k] * ab_z + g_0_y_0_xxxyzzzz[k];

                g_0_y_z_xxxzzzz[k] = -g_0_y_0_xxxzzzz[k] * ab_z + g_0_y_0_xxxzzzzz[k];

                g_0_y_z_xxyyyyy[k] = -g_0_y_0_xxyyyyy[k] * ab_z + g_0_y_0_xxyyyyyz[k];

                g_0_y_z_xxyyyyz[k] = -g_0_y_0_xxyyyyz[k] * ab_z + g_0_y_0_xxyyyyzz[k];

                g_0_y_z_xxyyyzz[k] = -g_0_y_0_xxyyyzz[k] * ab_z + g_0_y_0_xxyyyzzz[k];

                g_0_y_z_xxyyzzz[k] = -g_0_y_0_xxyyzzz[k] * ab_z + g_0_y_0_xxyyzzzz[k];

                g_0_y_z_xxyzzzz[k] = -g_0_y_0_xxyzzzz[k] * ab_z + g_0_y_0_xxyzzzzz[k];

                g_0_y_z_xxzzzzz[k] = -g_0_y_0_xxzzzzz[k] * ab_z + g_0_y_0_xxzzzzzz[k];

                g_0_y_z_xyyyyyy[k] = -g_0_y_0_xyyyyyy[k] * ab_z + g_0_y_0_xyyyyyyz[k];

                g_0_y_z_xyyyyyz[k] = -g_0_y_0_xyyyyyz[k] * ab_z + g_0_y_0_xyyyyyzz[k];

                g_0_y_z_xyyyyzz[k] = -g_0_y_0_xyyyyzz[k] * ab_z + g_0_y_0_xyyyyzzz[k];

                g_0_y_z_xyyyzzz[k] = -g_0_y_0_xyyyzzz[k] * ab_z + g_0_y_0_xyyyzzzz[k];

                g_0_y_z_xyyzzzz[k] = -g_0_y_0_xyyzzzz[k] * ab_z + g_0_y_0_xyyzzzzz[k];

                g_0_y_z_xyzzzzz[k] = -g_0_y_0_xyzzzzz[k] * ab_z + g_0_y_0_xyzzzzzz[k];

                g_0_y_z_xzzzzzz[k] = -g_0_y_0_xzzzzzz[k] * ab_z + g_0_y_0_xzzzzzzz[k];

                g_0_y_z_yyyyyyy[k] = -g_0_y_0_yyyyyyy[k] * ab_z + g_0_y_0_yyyyyyyz[k];

                g_0_y_z_yyyyyyz[k] = -g_0_y_0_yyyyyyz[k] * ab_z + g_0_y_0_yyyyyyzz[k];

                g_0_y_z_yyyyyzz[k] = -g_0_y_0_yyyyyzz[k] * ab_z + g_0_y_0_yyyyyzzz[k];

                g_0_y_z_yyyyzzz[k] = -g_0_y_0_yyyyzzz[k] * ab_z + g_0_y_0_yyyyzzzz[k];

                g_0_y_z_yyyzzzz[k] = -g_0_y_0_yyyzzzz[k] * ab_z + g_0_y_0_yyyzzzzz[k];

                g_0_y_z_yyzzzzz[k] = -g_0_y_0_yyzzzzz[k] * ab_z + g_0_y_0_yyzzzzzz[k];

                g_0_y_z_yzzzzzz[k] = -g_0_y_0_yzzzzzz[k] * ab_z + g_0_y_0_yzzzzzzz[k];

                g_0_y_z_zzzzzzz[k] = -g_0_y_0_zzzzzzz[k] * ab_z + g_0_y_0_zzzzzzzz[k];
            }

            /// Set up 216-252 components of targeted buffer : cbuffer.data(

            auto g_0_z_x_xxxxxxx = cbuffer.data(pk_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxy = cbuffer.data(pk_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_x_xxxxxxz = cbuffer.data(pk_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_x_xxxxxyy = cbuffer.data(pk_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_x_xxxxxyz = cbuffer.data(pk_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_x_xxxxxzz = cbuffer.data(pk_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_x_xxxxyyy = cbuffer.data(pk_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_x_xxxxyyz = cbuffer.data(pk_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_x_xxxxyzz = cbuffer.data(pk_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_x_xxxxzzz = cbuffer.data(pk_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_x_xxxyyyy = cbuffer.data(pk_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_x_xxxyyyz = cbuffer.data(pk_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_x_xxxyyzz = cbuffer.data(pk_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_x_xxxyzzz = cbuffer.data(pk_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_x_xxxzzzz = cbuffer.data(pk_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_x_xxyyyyy = cbuffer.data(pk_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_x_xxyyyyz = cbuffer.data(pk_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_x_xxyyyzz = cbuffer.data(pk_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_x_xxyyzzz = cbuffer.data(pk_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_x_xxyzzzz = cbuffer.data(pk_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_x_xxzzzzz = cbuffer.data(pk_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_x_xyyyyyy = cbuffer.data(pk_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_x_xyyyyyz = cbuffer.data(pk_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_x_xyyyyzz = cbuffer.data(pk_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_x_xyyyzzz = cbuffer.data(pk_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_x_xyyzzzz = cbuffer.data(pk_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_x_xyzzzzz = cbuffer.data(pk_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_x_xzzzzzz = cbuffer.data(pk_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_x_yyyyyyy = cbuffer.data(pk_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_x_yyyyyyz = cbuffer.data(pk_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_x_yyyyyzz = cbuffer.data(pk_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_x_yyyyzzz = cbuffer.data(pk_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_x_yyyzzzz = cbuffer.data(pk_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_x_yyzzzzz = cbuffer.data(pk_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_x_yzzzzzz = cbuffer.data(pk_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_x_zzzzzzz = cbuffer.data(pk_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxxxx, g_0_z_0_xxxxxxxx, g_0_z_0_xxxxxxxy, g_0_z_0_xxxxxxxz, g_0_z_0_xxxxxxy, g_0_z_0_xxxxxxyy, g_0_z_0_xxxxxxyz, g_0_z_0_xxxxxxz, g_0_z_0_xxxxxxzz, g_0_z_0_xxxxxyy, g_0_z_0_xxxxxyyy, g_0_z_0_xxxxxyyz, g_0_z_0_xxxxxyz, g_0_z_0_xxxxxyzz, g_0_z_0_xxxxxzz, g_0_z_0_xxxxxzzz, g_0_z_0_xxxxyyy, g_0_z_0_xxxxyyyy, g_0_z_0_xxxxyyyz, g_0_z_0_xxxxyyz, g_0_z_0_xxxxyyzz, g_0_z_0_xxxxyzz, g_0_z_0_xxxxyzzz, g_0_z_0_xxxxzzz, g_0_z_0_xxxxzzzz, g_0_z_0_xxxyyyy, g_0_z_0_xxxyyyyy, g_0_z_0_xxxyyyyz, g_0_z_0_xxxyyyz, g_0_z_0_xxxyyyzz, g_0_z_0_xxxyyzz, g_0_z_0_xxxyyzzz, g_0_z_0_xxxyzzz, g_0_z_0_xxxyzzzz, g_0_z_0_xxxzzzz, g_0_z_0_xxxzzzzz, g_0_z_0_xxyyyyy, g_0_z_0_xxyyyyyy, g_0_z_0_xxyyyyyz, g_0_z_0_xxyyyyz, g_0_z_0_xxyyyyzz, g_0_z_0_xxyyyzz, g_0_z_0_xxyyyzzz, g_0_z_0_xxyyzzz, g_0_z_0_xxyyzzzz, g_0_z_0_xxyzzzz, g_0_z_0_xxyzzzzz, g_0_z_0_xxzzzzz, g_0_z_0_xxzzzzzz, g_0_z_0_xyyyyyy, g_0_z_0_xyyyyyyy, g_0_z_0_xyyyyyyz, g_0_z_0_xyyyyyz, g_0_z_0_xyyyyyzz, g_0_z_0_xyyyyzz, g_0_z_0_xyyyyzzz, g_0_z_0_xyyyzzz, g_0_z_0_xyyyzzzz, g_0_z_0_xyyzzzz, g_0_z_0_xyyzzzzz, g_0_z_0_xyzzzzz, g_0_z_0_xyzzzzzz, g_0_z_0_xzzzzzz, g_0_z_0_xzzzzzzz, g_0_z_0_yyyyyyy, g_0_z_0_yyyyyyz, g_0_z_0_yyyyyzz, g_0_z_0_yyyyzzz, g_0_z_0_yyyzzzz, g_0_z_0_yyzzzzz, g_0_z_0_yzzzzzz, g_0_z_0_zzzzzzz, g_0_z_x_xxxxxxx, g_0_z_x_xxxxxxy, g_0_z_x_xxxxxxz, g_0_z_x_xxxxxyy, g_0_z_x_xxxxxyz, g_0_z_x_xxxxxzz, g_0_z_x_xxxxyyy, g_0_z_x_xxxxyyz, g_0_z_x_xxxxyzz, g_0_z_x_xxxxzzz, g_0_z_x_xxxyyyy, g_0_z_x_xxxyyyz, g_0_z_x_xxxyyzz, g_0_z_x_xxxyzzz, g_0_z_x_xxxzzzz, g_0_z_x_xxyyyyy, g_0_z_x_xxyyyyz, g_0_z_x_xxyyyzz, g_0_z_x_xxyyzzz, g_0_z_x_xxyzzzz, g_0_z_x_xxzzzzz, g_0_z_x_xyyyyyy, g_0_z_x_xyyyyyz, g_0_z_x_xyyyyzz, g_0_z_x_xyyyzzz, g_0_z_x_xyyzzzz, g_0_z_x_xyzzzzz, g_0_z_x_xzzzzzz, g_0_z_x_yyyyyyy, g_0_z_x_yyyyyyz, g_0_z_x_yyyyyzz, g_0_z_x_yyyyzzz, g_0_z_x_yyyzzzz, g_0_z_x_yyzzzzz, g_0_z_x_yzzzzzz, g_0_z_x_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_x_xxxxxxx[k] = -g_0_z_0_xxxxxxx[k] * ab_x + g_0_z_0_xxxxxxxx[k];

                g_0_z_x_xxxxxxy[k] = -g_0_z_0_xxxxxxy[k] * ab_x + g_0_z_0_xxxxxxxy[k];

                g_0_z_x_xxxxxxz[k] = -g_0_z_0_xxxxxxz[k] * ab_x + g_0_z_0_xxxxxxxz[k];

                g_0_z_x_xxxxxyy[k] = -g_0_z_0_xxxxxyy[k] * ab_x + g_0_z_0_xxxxxxyy[k];

                g_0_z_x_xxxxxyz[k] = -g_0_z_0_xxxxxyz[k] * ab_x + g_0_z_0_xxxxxxyz[k];

                g_0_z_x_xxxxxzz[k] = -g_0_z_0_xxxxxzz[k] * ab_x + g_0_z_0_xxxxxxzz[k];

                g_0_z_x_xxxxyyy[k] = -g_0_z_0_xxxxyyy[k] * ab_x + g_0_z_0_xxxxxyyy[k];

                g_0_z_x_xxxxyyz[k] = -g_0_z_0_xxxxyyz[k] * ab_x + g_0_z_0_xxxxxyyz[k];

                g_0_z_x_xxxxyzz[k] = -g_0_z_0_xxxxyzz[k] * ab_x + g_0_z_0_xxxxxyzz[k];

                g_0_z_x_xxxxzzz[k] = -g_0_z_0_xxxxzzz[k] * ab_x + g_0_z_0_xxxxxzzz[k];

                g_0_z_x_xxxyyyy[k] = -g_0_z_0_xxxyyyy[k] * ab_x + g_0_z_0_xxxxyyyy[k];

                g_0_z_x_xxxyyyz[k] = -g_0_z_0_xxxyyyz[k] * ab_x + g_0_z_0_xxxxyyyz[k];

                g_0_z_x_xxxyyzz[k] = -g_0_z_0_xxxyyzz[k] * ab_x + g_0_z_0_xxxxyyzz[k];

                g_0_z_x_xxxyzzz[k] = -g_0_z_0_xxxyzzz[k] * ab_x + g_0_z_0_xxxxyzzz[k];

                g_0_z_x_xxxzzzz[k] = -g_0_z_0_xxxzzzz[k] * ab_x + g_0_z_0_xxxxzzzz[k];

                g_0_z_x_xxyyyyy[k] = -g_0_z_0_xxyyyyy[k] * ab_x + g_0_z_0_xxxyyyyy[k];

                g_0_z_x_xxyyyyz[k] = -g_0_z_0_xxyyyyz[k] * ab_x + g_0_z_0_xxxyyyyz[k];

                g_0_z_x_xxyyyzz[k] = -g_0_z_0_xxyyyzz[k] * ab_x + g_0_z_0_xxxyyyzz[k];

                g_0_z_x_xxyyzzz[k] = -g_0_z_0_xxyyzzz[k] * ab_x + g_0_z_0_xxxyyzzz[k];

                g_0_z_x_xxyzzzz[k] = -g_0_z_0_xxyzzzz[k] * ab_x + g_0_z_0_xxxyzzzz[k];

                g_0_z_x_xxzzzzz[k] = -g_0_z_0_xxzzzzz[k] * ab_x + g_0_z_0_xxxzzzzz[k];

                g_0_z_x_xyyyyyy[k] = -g_0_z_0_xyyyyyy[k] * ab_x + g_0_z_0_xxyyyyyy[k];

                g_0_z_x_xyyyyyz[k] = -g_0_z_0_xyyyyyz[k] * ab_x + g_0_z_0_xxyyyyyz[k];

                g_0_z_x_xyyyyzz[k] = -g_0_z_0_xyyyyzz[k] * ab_x + g_0_z_0_xxyyyyzz[k];

                g_0_z_x_xyyyzzz[k] = -g_0_z_0_xyyyzzz[k] * ab_x + g_0_z_0_xxyyyzzz[k];

                g_0_z_x_xyyzzzz[k] = -g_0_z_0_xyyzzzz[k] * ab_x + g_0_z_0_xxyyzzzz[k];

                g_0_z_x_xyzzzzz[k] = -g_0_z_0_xyzzzzz[k] * ab_x + g_0_z_0_xxyzzzzz[k];

                g_0_z_x_xzzzzzz[k] = -g_0_z_0_xzzzzzz[k] * ab_x + g_0_z_0_xxzzzzzz[k];

                g_0_z_x_yyyyyyy[k] = -g_0_z_0_yyyyyyy[k] * ab_x + g_0_z_0_xyyyyyyy[k];

                g_0_z_x_yyyyyyz[k] = -g_0_z_0_yyyyyyz[k] * ab_x + g_0_z_0_xyyyyyyz[k];

                g_0_z_x_yyyyyzz[k] = -g_0_z_0_yyyyyzz[k] * ab_x + g_0_z_0_xyyyyyzz[k];

                g_0_z_x_yyyyzzz[k] = -g_0_z_0_yyyyzzz[k] * ab_x + g_0_z_0_xyyyyzzz[k];

                g_0_z_x_yyyzzzz[k] = -g_0_z_0_yyyzzzz[k] * ab_x + g_0_z_0_xyyyzzzz[k];

                g_0_z_x_yyzzzzz[k] = -g_0_z_0_yyzzzzz[k] * ab_x + g_0_z_0_xyyzzzzz[k];

                g_0_z_x_yzzzzzz[k] = -g_0_z_0_yzzzzzz[k] * ab_x + g_0_z_0_xyzzzzzz[k];

                g_0_z_x_zzzzzzz[k] = -g_0_z_0_zzzzzzz[k] * ab_x + g_0_z_0_xzzzzzzz[k];
            }

            /// Set up 252-288 components of targeted buffer : cbuffer.data(

            auto g_0_z_y_xxxxxxx = cbuffer.data(pk_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxy = cbuffer.data(pk_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_y_xxxxxxz = cbuffer.data(pk_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_z_y_xxxxxyy = cbuffer.data(pk_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_y_xxxxxyz = cbuffer.data(pk_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_y_xxxxxzz = cbuffer.data(pk_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_z_y_xxxxyyy = cbuffer.data(pk_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_y_xxxxyyz = cbuffer.data(pk_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_y_xxxxyzz = cbuffer.data(pk_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_z_y_xxxxzzz = cbuffer.data(pk_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_y_xxxyyyy = cbuffer.data(pk_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_y_xxxyyyz = cbuffer.data(pk_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_z_y_xxxyyzz = cbuffer.data(pk_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_y_xxxyzzz = cbuffer.data(pk_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_y_xxxzzzz = cbuffer.data(pk_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_z_y_xxyyyyy = cbuffer.data(pk_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_y_xxyyyyz = cbuffer.data(pk_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_y_xxyyyzz = cbuffer.data(pk_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_z_y_xxyyzzz = cbuffer.data(pk_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_y_xxyzzzz = cbuffer.data(pk_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_y_xxzzzzz = cbuffer.data(pk_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_z_y_xyyyyyy = cbuffer.data(pk_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_y_xyyyyyz = cbuffer.data(pk_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_y_xyyyyzz = cbuffer.data(pk_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_z_y_xyyyzzz = cbuffer.data(pk_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_y_xyyzzzz = cbuffer.data(pk_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_y_xyzzzzz = cbuffer.data(pk_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_z_y_xzzzzzz = cbuffer.data(pk_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_y_yyyyyyy = cbuffer.data(pk_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_y_yyyyyyz = cbuffer.data(pk_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_z_y_yyyyyzz = cbuffer.data(pk_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_y_yyyyzzz = cbuffer.data(pk_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_y_yyyzzzz = cbuffer.data(pk_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_z_y_yyzzzzz = cbuffer.data(pk_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_y_yzzzzzz = cbuffer.data(pk_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_y_zzzzzzz = cbuffer.data(pk_geom_01_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxxxx, g_0_z_0_xxxxxxxy, g_0_z_0_xxxxxxy, g_0_z_0_xxxxxxyy, g_0_z_0_xxxxxxyz, g_0_z_0_xxxxxxz, g_0_z_0_xxxxxyy, g_0_z_0_xxxxxyyy, g_0_z_0_xxxxxyyz, g_0_z_0_xxxxxyz, g_0_z_0_xxxxxyzz, g_0_z_0_xxxxxzz, g_0_z_0_xxxxyyy, g_0_z_0_xxxxyyyy, g_0_z_0_xxxxyyyz, g_0_z_0_xxxxyyz, g_0_z_0_xxxxyyzz, g_0_z_0_xxxxyzz, g_0_z_0_xxxxyzzz, g_0_z_0_xxxxzzz, g_0_z_0_xxxyyyy, g_0_z_0_xxxyyyyy, g_0_z_0_xxxyyyyz, g_0_z_0_xxxyyyz, g_0_z_0_xxxyyyzz, g_0_z_0_xxxyyzz, g_0_z_0_xxxyyzzz, g_0_z_0_xxxyzzz, g_0_z_0_xxxyzzzz, g_0_z_0_xxxzzzz, g_0_z_0_xxyyyyy, g_0_z_0_xxyyyyyy, g_0_z_0_xxyyyyyz, g_0_z_0_xxyyyyz, g_0_z_0_xxyyyyzz, g_0_z_0_xxyyyzz, g_0_z_0_xxyyyzzz, g_0_z_0_xxyyzzz, g_0_z_0_xxyyzzzz, g_0_z_0_xxyzzzz, g_0_z_0_xxyzzzzz, g_0_z_0_xxzzzzz, g_0_z_0_xyyyyyy, g_0_z_0_xyyyyyyy, g_0_z_0_xyyyyyyz, g_0_z_0_xyyyyyz, g_0_z_0_xyyyyyzz, g_0_z_0_xyyyyzz, g_0_z_0_xyyyyzzz, g_0_z_0_xyyyzzz, g_0_z_0_xyyyzzzz, g_0_z_0_xyyzzzz, g_0_z_0_xyyzzzzz, g_0_z_0_xyzzzzz, g_0_z_0_xyzzzzzz, g_0_z_0_xzzzzzz, g_0_z_0_yyyyyyy, g_0_z_0_yyyyyyyy, g_0_z_0_yyyyyyyz, g_0_z_0_yyyyyyz, g_0_z_0_yyyyyyzz, g_0_z_0_yyyyyzz, g_0_z_0_yyyyyzzz, g_0_z_0_yyyyzzz, g_0_z_0_yyyyzzzz, g_0_z_0_yyyzzzz, g_0_z_0_yyyzzzzz, g_0_z_0_yyzzzzz, g_0_z_0_yyzzzzzz, g_0_z_0_yzzzzzz, g_0_z_0_yzzzzzzz, g_0_z_0_zzzzzzz, g_0_z_y_xxxxxxx, g_0_z_y_xxxxxxy, g_0_z_y_xxxxxxz, g_0_z_y_xxxxxyy, g_0_z_y_xxxxxyz, g_0_z_y_xxxxxzz, g_0_z_y_xxxxyyy, g_0_z_y_xxxxyyz, g_0_z_y_xxxxyzz, g_0_z_y_xxxxzzz, g_0_z_y_xxxyyyy, g_0_z_y_xxxyyyz, g_0_z_y_xxxyyzz, g_0_z_y_xxxyzzz, g_0_z_y_xxxzzzz, g_0_z_y_xxyyyyy, g_0_z_y_xxyyyyz, g_0_z_y_xxyyyzz, g_0_z_y_xxyyzzz, g_0_z_y_xxyzzzz, g_0_z_y_xxzzzzz, g_0_z_y_xyyyyyy, g_0_z_y_xyyyyyz, g_0_z_y_xyyyyzz, g_0_z_y_xyyyzzz, g_0_z_y_xyyzzzz, g_0_z_y_xyzzzzz, g_0_z_y_xzzzzzz, g_0_z_y_yyyyyyy, g_0_z_y_yyyyyyz, g_0_z_y_yyyyyzz, g_0_z_y_yyyyzzz, g_0_z_y_yyyzzzz, g_0_z_y_yyzzzzz, g_0_z_y_yzzzzzz, g_0_z_y_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_y_xxxxxxx[k] = -g_0_z_0_xxxxxxx[k] * ab_y + g_0_z_0_xxxxxxxy[k];

                g_0_z_y_xxxxxxy[k] = -g_0_z_0_xxxxxxy[k] * ab_y + g_0_z_0_xxxxxxyy[k];

                g_0_z_y_xxxxxxz[k] = -g_0_z_0_xxxxxxz[k] * ab_y + g_0_z_0_xxxxxxyz[k];

                g_0_z_y_xxxxxyy[k] = -g_0_z_0_xxxxxyy[k] * ab_y + g_0_z_0_xxxxxyyy[k];

                g_0_z_y_xxxxxyz[k] = -g_0_z_0_xxxxxyz[k] * ab_y + g_0_z_0_xxxxxyyz[k];

                g_0_z_y_xxxxxzz[k] = -g_0_z_0_xxxxxzz[k] * ab_y + g_0_z_0_xxxxxyzz[k];

                g_0_z_y_xxxxyyy[k] = -g_0_z_0_xxxxyyy[k] * ab_y + g_0_z_0_xxxxyyyy[k];

                g_0_z_y_xxxxyyz[k] = -g_0_z_0_xxxxyyz[k] * ab_y + g_0_z_0_xxxxyyyz[k];

                g_0_z_y_xxxxyzz[k] = -g_0_z_0_xxxxyzz[k] * ab_y + g_0_z_0_xxxxyyzz[k];

                g_0_z_y_xxxxzzz[k] = -g_0_z_0_xxxxzzz[k] * ab_y + g_0_z_0_xxxxyzzz[k];

                g_0_z_y_xxxyyyy[k] = -g_0_z_0_xxxyyyy[k] * ab_y + g_0_z_0_xxxyyyyy[k];

                g_0_z_y_xxxyyyz[k] = -g_0_z_0_xxxyyyz[k] * ab_y + g_0_z_0_xxxyyyyz[k];

                g_0_z_y_xxxyyzz[k] = -g_0_z_0_xxxyyzz[k] * ab_y + g_0_z_0_xxxyyyzz[k];

                g_0_z_y_xxxyzzz[k] = -g_0_z_0_xxxyzzz[k] * ab_y + g_0_z_0_xxxyyzzz[k];

                g_0_z_y_xxxzzzz[k] = -g_0_z_0_xxxzzzz[k] * ab_y + g_0_z_0_xxxyzzzz[k];

                g_0_z_y_xxyyyyy[k] = -g_0_z_0_xxyyyyy[k] * ab_y + g_0_z_0_xxyyyyyy[k];

                g_0_z_y_xxyyyyz[k] = -g_0_z_0_xxyyyyz[k] * ab_y + g_0_z_0_xxyyyyyz[k];

                g_0_z_y_xxyyyzz[k] = -g_0_z_0_xxyyyzz[k] * ab_y + g_0_z_0_xxyyyyzz[k];

                g_0_z_y_xxyyzzz[k] = -g_0_z_0_xxyyzzz[k] * ab_y + g_0_z_0_xxyyyzzz[k];

                g_0_z_y_xxyzzzz[k] = -g_0_z_0_xxyzzzz[k] * ab_y + g_0_z_0_xxyyzzzz[k];

                g_0_z_y_xxzzzzz[k] = -g_0_z_0_xxzzzzz[k] * ab_y + g_0_z_0_xxyzzzzz[k];

                g_0_z_y_xyyyyyy[k] = -g_0_z_0_xyyyyyy[k] * ab_y + g_0_z_0_xyyyyyyy[k];

                g_0_z_y_xyyyyyz[k] = -g_0_z_0_xyyyyyz[k] * ab_y + g_0_z_0_xyyyyyyz[k];

                g_0_z_y_xyyyyzz[k] = -g_0_z_0_xyyyyzz[k] * ab_y + g_0_z_0_xyyyyyzz[k];

                g_0_z_y_xyyyzzz[k] = -g_0_z_0_xyyyzzz[k] * ab_y + g_0_z_0_xyyyyzzz[k];

                g_0_z_y_xyyzzzz[k] = -g_0_z_0_xyyzzzz[k] * ab_y + g_0_z_0_xyyyzzzz[k];

                g_0_z_y_xyzzzzz[k] = -g_0_z_0_xyzzzzz[k] * ab_y + g_0_z_0_xyyzzzzz[k];

                g_0_z_y_xzzzzzz[k] = -g_0_z_0_xzzzzzz[k] * ab_y + g_0_z_0_xyzzzzzz[k];

                g_0_z_y_yyyyyyy[k] = -g_0_z_0_yyyyyyy[k] * ab_y + g_0_z_0_yyyyyyyy[k];

                g_0_z_y_yyyyyyz[k] = -g_0_z_0_yyyyyyz[k] * ab_y + g_0_z_0_yyyyyyyz[k];

                g_0_z_y_yyyyyzz[k] = -g_0_z_0_yyyyyzz[k] * ab_y + g_0_z_0_yyyyyyzz[k];

                g_0_z_y_yyyyzzz[k] = -g_0_z_0_yyyyzzz[k] * ab_y + g_0_z_0_yyyyyzzz[k];

                g_0_z_y_yyyzzzz[k] = -g_0_z_0_yyyzzzz[k] * ab_y + g_0_z_0_yyyyzzzz[k];

                g_0_z_y_yyzzzzz[k] = -g_0_z_0_yyzzzzz[k] * ab_y + g_0_z_0_yyyzzzzz[k];

                g_0_z_y_yzzzzzz[k] = -g_0_z_0_yzzzzzz[k] * ab_y + g_0_z_0_yyzzzzzz[k];

                g_0_z_y_zzzzzzz[k] = -g_0_z_0_zzzzzzz[k] * ab_y + g_0_z_0_yzzzzzzz[k];
            }

            /// Set up 288-324 components of targeted buffer : cbuffer.data(

            auto g_0_z_z_xxxxxxx = cbuffer.data(pk_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxy = cbuffer.data(pk_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_z_xxxxxxz = cbuffer.data(pk_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_z_z_xxxxxyy = cbuffer.data(pk_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_z_xxxxxyz = cbuffer.data(pk_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_z_xxxxxzz = cbuffer.data(pk_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_z_z_xxxxyyy = cbuffer.data(pk_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_z_xxxxyyz = cbuffer.data(pk_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_z_xxxxyzz = cbuffer.data(pk_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_z_z_xxxxzzz = cbuffer.data(pk_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_z_xxxyyyy = cbuffer.data(pk_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_z_xxxyyyz = cbuffer.data(pk_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_z_z_xxxyyzz = cbuffer.data(pk_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_z_xxxyzzz = cbuffer.data(pk_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_z_xxxzzzz = cbuffer.data(pk_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_z_z_xxyyyyy = cbuffer.data(pk_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_z_xxyyyyz = cbuffer.data(pk_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_z_xxyyyzz = cbuffer.data(pk_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_z_z_xxyyzzz = cbuffer.data(pk_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_z_xxyzzzz = cbuffer.data(pk_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_z_xxzzzzz = cbuffer.data(pk_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_z_z_xyyyyyy = cbuffer.data(pk_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_z_xyyyyyz = cbuffer.data(pk_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_z_xyyyyzz = cbuffer.data(pk_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_z_z_xyyyzzz = cbuffer.data(pk_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_z_xyyzzzz = cbuffer.data(pk_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_z_xyzzzzz = cbuffer.data(pk_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_z_z_xzzzzzz = cbuffer.data(pk_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_z_yyyyyyy = cbuffer.data(pk_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_z_yyyyyyz = cbuffer.data(pk_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_z_z_yyyyyzz = cbuffer.data(pk_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_z_yyyyzzz = cbuffer.data(pk_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_z_yyyzzzz = cbuffer.data(pk_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_z_z_yyzzzzz = cbuffer.data(pk_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_z_yzzzzzz = cbuffer.data(pk_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_z_zzzzzzz = cbuffer.data(pk_geom_01_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxxxx, g_0_xxxxxxy, g_0_xxxxxxz, g_0_xxxxxyy, g_0_xxxxxyz, g_0_xxxxxzz, g_0_xxxxyyy, g_0_xxxxyyz, g_0_xxxxyzz, g_0_xxxxzzz, g_0_xxxyyyy, g_0_xxxyyyz, g_0_xxxyyzz, g_0_xxxyzzz, g_0_xxxzzzz, g_0_xxyyyyy, g_0_xxyyyyz, g_0_xxyyyzz, g_0_xxyyzzz, g_0_xxyzzzz, g_0_xxzzzzz, g_0_xyyyyyy, g_0_xyyyyyz, g_0_xyyyyzz, g_0_xyyyzzz, g_0_xyyzzzz, g_0_xyzzzzz, g_0_xzzzzzz, g_0_yyyyyyy, g_0_yyyyyyz, g_0_yyyyyzz, g_0_yyyyzzz, g_0_yyyzzzz, g_0_yyzzzzz, g_0_yzzzzzz, g_0_z_0_xxxxxxx, g_0_z_0_xxxxxxxz, g_0_z_0_xxxxxxy, g_0_z_0_xxxxxxyz, g_0_z_0_xxxxxxz, g_0_z_0_xxxxxxzz, g_0_z_0_xxxxxyy, g_0_z_0_xxxxxyyz, g_0_z_0_xxxxxyz, g_0_z_0_xxxxxyzz, g_0_z_0_xxxxxzz, g_0_z_0_xxxxxzzz, g_0_z_0_xxxxyyy, g_0_z_0_xxxxyyyz, g_0_z_0_xxxxyyz, g_0_z_0_xxxxyyzz, g_0_z_0_xxxxyzz, g_0_z_0_xxxxyzzz, g_0_z_0_xxxxzzz, g_0_z_0_xxxxzzzz, g_0_z_0_xxxyyyy, g_0_z_0_xxxyyyyz, g_0_z_0_xxxyyyz, g_0_z_0_xxxyyyzz, g_0_z_0_xxxyyzz, g_0_z_0_xxxyyzzz, g_0_z_0_xxxyzzz, g_0_z_0_xxxyzzzz, g_0_z_0_xxxzzzz, g_0_z_0_xxxzzzzz, g_0_z_0_xxyyyyy, g_0_z_0_xxyyyyyz, g_0_z_0_xxyyyyz, g_0_z_0_xxyyyyzz, g_0_z_0_xxyyyzz, g_0_z_0_xxyyyzzz, g_0_z_0_xxyyzzz, g_0_z_0_xxyyzzzz, g_0_z_0_xxyzzzz, g_0_z_0_xxyzzzzz, g_0_z_0_xxzzzzz, g_0_z_0_xxzzzzzz, g_0_z_0_xyyyyyy, g_0_z_0_xyyyyyyz, g_0_z_0_xyyyyyz, g_0_z_0_xyyyyyzz, g_0_z_0_xyyyyzz, g_0_z_0_xyyyyzzz, g_0_z_0_xyyyzzz, g_0_z_0_xyyyzzzz, g_0_z_0_xyyzzzz, g_0_z_0_xyyzzzzz, g_0_z_0_xyzzzzz, g_0_z_0_xyzzzzzz, g_0_z_0_xzzzzzz, g_0_z_0_xzzzzzzz, g_0_z_0_yyyyyyy, g_0_z_0_yyyyyyyz, g_0_z_0_yyyyyyz, g_0_z_0_yyyyyyzz, g_0_z_0_yyyyyzz, g_0_z_0_yyyyyzzz, g_0_z_0_yyyyzzz, g_0_z_0_yyyyzzzz, g_0_z_0_yyyzzzz, g_0_z_0_yyyzzzzz, g_0_z_0_yyzzzzz, g_0_z_0_yyzzzzzz, g_0_z_0_yzzzzzz, g_0_z_0_yzzzzzzz, g_0_z_0_zzzzzzz, g_0_z_0_zzzzzzzz, g_0_z_z_xxxxxxx, g_0_z_z_xxxxxxy, g_0_z_z_xxxxxxz, g_0_z_z_xxxxxyy, g_0_z_z_xxxxxyz, g_0_z_z_xxxxxzz, g_0_z_z_xxxxyyy, g_0_z_z_xxxxyyz, g_0_z_z_xxxxyzz, g_0_z_z_xxxxzzz, g_0_z_z_xxxyyyy, g_0_z_z_xxxyyyz, g_0_z_z_xxxyyzz, g_0_z_z_xxxyzzz, g_0_z_z_xxxzzzz, g_0_z_z_xxyyyyy, g_0_z_z_xxyyyyz, g_0_z_z_xxyyyzz, g_0_z_z_xxyyzzz, g_0_z_z_xxyzzzz, g_0_z_z_xxzzzzz, g_0_z_z_xyyyyyy, g_0_z_z_xyyyyyz, g_0_z_z_xyyyyzz, g_0_z_z_xyyyzzz, g_0_z_z_xyyzzzz, g_0_z_z_xyzzzzz, g_0_z_z_xzzzzzz, g_0_z_z_yyyyyyy, g_0_z_z_yyyyyyz, g_0_z_z_yyyyyzz, g_0_z_z_yyyyzzz, g_0_z_z_yyyzzzz, g_0_z_z_yyzzzzz, g_0_z_z_yzzzzzz, g_0_z_z_zzzzzzz, g_0_zzzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_z_xxxxxxx[k] = g_0_xxxxxxx[k] - g_0_z_0_xxxxxxx[k] * ab_z + g_0_z_0_xxxxxxxz[k];

                g_0_z_z_xxxxxxy[k] = g_0_xxxxxxy[k] - g_0_z_0_xxxxxxy[k] * ab_z + g_0_z_0_xxxxxxyz[k];

                g_0_z_z_xxxxxxz[k] = g_0_xxxxxxz[k] - g_0_z_0_xxxxxxz[k] * ab_z + g_0_z_0_xxxxxxzz[k];

                g_0_z_z_xxxxxyy[k] = g_0_xxxxxyy[k] - g_0_z_0_xxxxxyy[k] * ab_z + g_0_z_0_xxxxxyyz[k];

                g_0_z_z_xxxxxyz[k] = g_0_xxxxxyz[k] - g_0_z_0_xxxxxyz[k] * ab_z + g_0_z_0_xxxxxyzz[k];

                g_0_z_z_xxxxxzz[k] = g_0_xxxxxzz[k] - g_0_z_0_xxxxxzz[k] * ab_z + g_0_z_0_xxxxxzzz[k];

                g_0_z_z_xxxxyyy[k] = g_0_xxxxyyy[k] - g_0_z_0_xxxxyyy[k] * ab_z + g_0_z_0_xxxxyyyz[k];

                g_0_z_z_xxxxyyz[k] = g_0_xxxxyyz[k] - g_0_z_0_xxxxyyz[k] * ab_z + g_0_z_0_xxxxyyzz[k];

                g_0_z_z_xxxxyzz[k] = g_0_xxxxyzz[k] - g_0_z_0_xxxxyzz[k] * ab_z + g_0_z_0_xxxxyzzz[k];

                g_0_z_z_xxxxzzz[k] = g_0_xxxxzzz[k] - g_0_z_0_xxxxzzz[k] * ab_z + g_0_z_0_xxxxzzzz[k];

                g_0_z_z_xxxyyyy[k] = g_0_xxxyyyy[k] - g_0_z_0_xxxyyyy[k] * ab_z + g_0_z_0_xxxyyyyz[k];

                g_0_z_z_xxxyyyz[k] = g_0_xxxyyyz[k] - g_0_z_0_xxxyyyz[k] * ab_z + g_0_z_0_xxxyyyzz[k];

                g_0_z_z_xxxyyzz[k] = g_0_xxxyyzz[k] - g_0_z_0_xxxyyzz[k] * ab_z + g_0_z_0_xxxyyzzz[k];

                g_0_z_z_xxxyzzz[k] = g_0_xxxyzzz[k] - g_0_z_0_xxxyzzz[k] * ab_z + g_0_z_0_xxxyzzzz[k];

                g_0_z_z_xxxzzzz[k] = g_0_xxxzzzz[k] - g_0_z_0_xxxzzzz[k] * ab_z + g_0_z_0_xxxzzzzz[k];

                g_0_z_z_xxyyyyy[k] = g_0_xxyyyyy[k] - g_0_z_0_xxyyyyy[k] * ab_z + g_0_z_0_xxyyyyyz[k];

                g_0_z_z_xxyyyyz[k] = g_0_xxyyyyz[k] - g_0_z_0_xxyyyyz[k] * ab_z + g_0_z_0_xxyyyyzz[k];

                g_0_z_z_xxyyyzz[k] = g_0_xxyyyzz[k] - g_0_z_0_xxyyyzz[k] * ab_z + g_0_z_0_xxyyyzzz[k];

                g_0_z_z_xxyyzzz[k] = g_0_xxyyzzz[k] - g_0_z_0_xxyyzzz[k] * ab_z + g_0_z_0_xxyyzzzz[k];

                g_0_z_z_xxyzzzz[k] = g_0_xxyzzzz[k] - g_0_z_0_xxyzzzz[k] * ab_z + g_0_z_0_xxyzzzzz[k];

                g_0_z_z_xxzzzzz[k] = g_0_xxzzzzz[k] - g_0_z_0_xxzzzzz[k] * ab_z + g_0_z_0_xxzzzzzz[k];

                g_0_z_z_xyyyyyy[k] = g_0_xyyyyyy[k] - g_0_z_0_xyyyyyy[k] * ab_z + g_0_z_0_xyyyyyyz[k];

                g_0_z_z_xyyyyyz[k] = g_0_xyyyyyz[k] - g_0_z_0_xyyyyyz[k] * ab_z + g_0_z_0_xyyyyyzz[k];

                g_0_z_z_xyyyyzz[k] = g_0_xyyyyzz[k] - g_0_z_0_xyyyyzz[k] * ab_z + g_0_z_0_xyyyyzzz[k];

                g_0_z_z_xyyyzzz[k] = g_0_xyyyzzz[k] - g_0_z_0_xyyyzzz[k] * ab_z + g_0_z_0_xyyyzzzz[k];

                g_0_z_z_xyyzzzz[k] = g_0_xyyzzzz[k] - g_0_z_0_xyyzzzz[k] * ab_z + g_0_z_0_xyyzzzzz[k];

                g_0_z_z_xyzzzzz[k] = g_0_xyzzzzz[k] - g_0_z_0_xyzzzzz[k] * ab_z + g_0_z_0_xyzzzzzz[k];

                g_0_z_z_xzzzzzz[k] = g_0_xzzzzzz[k] - g_0_z_0_xzzzzzz[k] * ab_z + g_0_z_0_xzzzzzzz[k];

                g_0_z_z_yyyyyyy[k] = g_0_yyyyyyy[k] - g_0_z_0_yyyyyyy[k] * ab_z + g_0_z_0_yyyyyyyz[k];

                g_0_z_z_yyyyyyz[k] = g_0_yyyyyyz[k] - g_0_z_0_yyyyyyz[k] * ab_z + g_0_z_0_yyyyyyzz[k];

                g_0_z_z_yyyyyzz[k] = g_0_yyyyyzz[k] - g_0_z_0_yyyyyzz[k] * ab_z + g_0_z_0_yyyyyzzz[k];

                g_0_z_z_yyyyzzz[k] = g_0_yyyyzzz[k] - g_0_z_0_yyyyzzz[k] * ab_z + g_0_z_0_yyyyzzzz[k];

                g_0_z_z_yyyzzzz[k] = g_0_yyyzzzz[k] - g_0_z_0_yyyzzzz[k] * ab_z + g_0_z_0_yyyzzzzz[k];

                g_0_z_z_yyzzzzz[k] = g_0_yyzzzzz[k] - g_0_z_0_yyzzzzz[k] * ab_z + g_0_z_0_yyzzzzzz[k];

                g_0_z_z_yzzzzzz[k] = g_0_yzzzzzz[k] - g_0_z_0_yzzzzzz[k] * ab_z + g_0_z_0_yzzzzzzz[k];

                g_0_z_z_zzzzzzz[k] = g_0_zzzzzzz[k] - g_0_z_0_zzzzzzz[k] * ab_z + g_0_z_0_zzzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

