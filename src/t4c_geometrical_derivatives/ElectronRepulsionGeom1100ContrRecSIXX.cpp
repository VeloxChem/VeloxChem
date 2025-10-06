#include "ElectronRepulsionGeom1100ContrRecSIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_sixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_sixx,
                                            const size_t idx_sixx,
                                            const size_t idx_geom_01_sixx,
                                            const size_t idx_geom_01_skxx,
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

            /// Set up components of auxilary buffer : SISS

            const auto si_geom_01_off = idx_geom_01_sixx + i * dcomps + j;

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

            /// set up bra offset for contr_buffer_sixx

            const auto si_geom_11_off = idx_geom_11_sixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_x_0_xxxxxx = cbuffer.data(si_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_0_xxxxxy = cbuffer.data(si_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_0_xxxxxz = cbuffer.data(si_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_0_xxxxyy = cbuffer.data(si_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_0_xxxxyz = cbuffer.data(si_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_0_xxxxzz = cbuffer.data(si_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_0_xxxyyy = cbuffer.data(si_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_0_xxxyyz = cbuffer.data(si_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_0_xxxyzz = cbuffer.data(si_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_0_xxxzzz = cbuffer.data(si_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_0_xxyyyy = cbuffer.data(si_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_0_xxyyyz = cbuffer.data(si_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_0_xxyyzz = cbuffer.data(si_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_0_xxyzzz = cbuffer.data(si_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_0_xxzzzz = cbuffer.data(si_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_0_xyyyyy = cbuffer.data(si_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_0_xyyyyz = cbuffer.data(si_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_0_xyyyzz = cbuffer.data(si_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_0_xyyzzz = cbuffer.data(si_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_0_xyzzzz = cbuffer.data(si_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_0_xzzzzz = cbuffer.data(si_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_0_yyyyyy = cbuffer.data(si_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_0_yyyyyz = cbuffer.data(si_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_0_yyyyzz = cbuffer.data(si_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_0_yyyzzz = cbuffer.data(si_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_0_yyzzzz = cbuffer.data(si_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_0_yzzzzz = cbuffer.data(si_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_0_zzzzzz = cbuffer.data(si_geom_11_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxx, g_0_x_0_xxxxxxx, g_0_x_0_xxxxxxy, g_0_x_0_xxxxxxz, g_0_x_0_xxxxxy, g_0_x_0_xxxxxyy, g_0_x_0_xxxxxyz, g_0_x_0_xxxxxz, g_0_x_0_xxxxxzz, g_0_x_0_xxxxyy, g_0_x_0_xxxxyyy, g_0_x_0_xxxxyyz, g_0_x_0_xxxxyz, g_0_x_0_xxxxyzz, g_0_x_0_xxxxzz, g_0_x_0_xxxxzzz, g_0_x_0_xxxyyy, g_0_x_0_xxxyyyy, g_0_x_0_xxxyyyz, g_0_x_0_xxxyyz, g_0_x_0_xxxyyzz, g_0_x_0_xxxyzz, g_0_x_0_xxxyzzz, g_0_x_0_xxxzzz, g_0_x_0_xxxzzzz, g_0_x_0_xxyyyy, g_0_x_0_xxyyyyy, g_0_x_0_xxyyyyz, g_0_x_0_xxyyyz, g_0_x_0_xxyyyzz, g_0_x_0_xxyyzz, g_0_x_0_xxyyzzz, g_0_x_0_xxyzzz, g_0_x_0_xxyzzzz, g_0_x_0_xxzzzz, g_0_x_0_xxzzzzz, g_0_x_0_xyyyyy, g_0_x_0_xyyyyyy, g_0_x_0_xyyyyyz, g_0_x_0_xyyyyz, g_0_x_0_xyyyyzz, g_0_x_0_xyyyzz, g_0_x_0_xyyyzzz, g_0_x_0_xyyzzz, g_0_x_0_xyyzzzz, g_0_x_0_xyzzzz, g_0_x_0_xyzzzzz, g_0_x_0_xzzzzz, g_0_x_0_xzzzzzz, g_0_x_0_yyyyyy, g_0_x_0_yyyyyz, g_0_x_0_yyyyzz, g_0_x_0_yyyzzz, g_0_x_0_yyzzzz, g_0_x_0_yzzzzz, g_0_x_0_zzzzzz, g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyzzz, g_0_xyzzzz, g_0_xzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_zzzzzz, g_x_x_0_xxxxxx, g_x_x_0_xxxxxy, g_x_x_0_xxxxxz, g_x_x_0_xxxxyy, g_x_x_0_xxxxyz, g_x_x_0_xxxxzz, g_x_x_0_xxxyyy, g_x_x_0_xxxyyz, g_x_x_0_xxxyzz, g_x_x_0_xxxzzz, g_x_x_0_xxyyyy, g_x_x_0_xxyyyz, g_x_x_0_xxyyzz, g_x_x_0_xxyzzz, g_x_x_0_xxzzzz, g_x_x_0_xyyyyy, g_x_x_0_xyyyyz, g_x_x_0_xyyyzz, g_x_x_0_xyyzzz, g_x_x_0_xyzzzz, g_x_x_0_xzzzzz, g_x_x_0_yyyyyy, g_x_x_0_yyyyyz, g_x_x_0_yyyyzz, g_x_x_0_yyyzzz, g_x_x_0_yyzzzz, g_x_x_0_yzzzzz, g_x_x_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_0_xxxxxx[k] = g_0_xxxxxx[k] - g_0_x_0_xxxxxx[k] * ab_x + g_0_x_0_xxxxxxx[k];

                g_x_x_0_xxxxxy[k] = g_0_xxxxxy[k] - g_0_x_0_xxxxxy[k] * ab_x + g_0_x_0_xxxxxxy[k];

                g_x_x_0_xxxxxz[k] = g_0_xxxxxz[k] - g_0_x_0_xxxxxz[k] * ab_x + g_0_x_0_xxxxxxz[k];

                g_x_x_0_xxxxyy[k] = g_0_xxxxyy[k] - g_0_x_0_xxxxyy[k] * ab_x + g_0_x_0_xxxxxyy[k];

                g_x_x_0_xxxxyz[k] = g_0_xxxxyz[k] - g_0_x_0_xxxxyz[k] * ab_x + g_0_x_0_xxxxxyz[k];

                g_x_x_0_xxxxzz[k] = g_0_xxxxzz[k] - g_0_x_0_xxxxzz[k] * ab_x + g_0_x_0_xxxxxzz[k];

                g_x_x_0_xxxyyy[k] = g_0_xxxyyy[k] - g_0_x_0_xxxyyy[k] * ab_x + g_0_x_0_xxxxyyy[k];

                g_x_x_0_xxxyyz[k] = g_0_xxxyyz[k] - g_0_x_0_xxxyyz[k] * ab_x + g_0_x_0_xxxxyyz[k];

                g_x_x_0_xxxyzz[k] = g_0_xxxyzz[k] - g_0_x_0_xxxyzz[k] * ab_x + g_0_x_0_xxxxyzz[k];

                g_x_x_0_xxxzzz[k] = g_0_xxxzzz[k] - g_0_x_0_xxxzzz[k] * ab_x + g_0_x_0_xxxxzzz[k];

                g_x_x_0_xxyyyy[k] = g_0_xxyyyy[k] - g_0_x_0_xxyyyy[k] * ab_x + g_0_x_0_xxxyyyy[k];

                g_x_x_0_xxyyyz[k] = g_0_xxyyyz[k] - g_0_x_0_xxyyyz[k] * ab_x + g_0_x_0_xxxyyyz[k];

                g_x_x_0_xxyyzz[k] = g_0_xxyyzz[k] - g_0_x_0_xxyyzz[k] * ab_x + g_0_x_0_xxxyyzz[k];

                g_x_x_0_xxyzzz[k] = g_0_xxyzzz[k] - g_0_x_0_xxyzzz[k] * ab_x + g_0_x_0_xxxyzzz[k];

                g_x_x_0_xxzzzz[k] = g_0_xxzzzz[k] - g_0_x_0_xxzzzz[k] * ab_x + g_0_x_0_xxxzzzz[k];

                g_x_x_0_xyyyyy[k] = g_0_xyyyyy[k] - g_0_x_0_xyyyyy[k] * ab_x + g_0_x_0_xxyyyyy[k];

                g_x_x_0_xyyyyz[k] = g_0_xyyyyz[k] - g_0_x_0_xyyyyz[k] * ab_x + g_0_x_0_xxyyyyz[k];

                g_x_x_0_xyyyzz[k] = g_0_xyyyzz[k] - g_0_x_0_xyyyzz[k] * ab_x + g_0_x_0_xxyyyzz[k];

                g_x_x_0_xyyzzz[k] = g_0_xyyzzz[k] - g_0_x_0_xyyzzz[k] * ab_x + g_0_x_0_xxyyzzz[k];

                g_x_x_0_xyzzzz[k] = g_0_xyzzzz[k] - g_0_x_0_xyzzzz[k] * ab_x + g_0_x_0_xxyzzzz[k];

                g_x_x_0_xzzzzz[k] = g_0_xzzzzz[k] - g_0_x_0_xzzzzz[k] * ab_x + g_0_x_0_xxzzzzz[k];

                g_x_x_0_yyyyyy[k] = g_0_yyyyyy[k] - g_0_x_0_yyyyyy[k] * ab_x + g_0_x_0_xyyyyyy[k];

                g_x_x_0_yyyyyz[k] = g_0_yyyyyz[k] - g_0_x_0_yyyyyz[k] * ab_x + g_0_x_0_xyyyyyz[k];

                g_x_x_0_yyyyzz[k] = g_0_yyyyzz[k] - g_0_x_0_yyyyzz[k] * ab_x + g_0_x_0_xyyyyzz[k];

                g_x_x_0_yyyzzz[k] = g_0_yyyzzz[k] - g_0_x_0_yyyzzz[k] * ab_x + g_0_x_0_xyyyzzz[k];

                g_x_x_0_yyzzzz[k] = g_0_yyzzzz[k] - g_0_x_0_yyzzzz[k] * ab_x + g_0_x_0_xyyzzzz[k];

                g_x_x_0_yzzzzz[k] = g_0_yzzzzz[k] - g_0_x_0_yzzzzz[k] * ab_x + g_0_x_0_xyzzzzz[k];

                g_x_x_0_zzzzzz[k] = g_0_zzzzzz[k] - g_0_x_0_zzzzzz[k] * ab_x + g_0_x_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_y_0_xxxxxx = cbuffer.data(si_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_y_0_xxxxxy = cbuffer.data(si_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_y_0_xxxxxz = cbuffer.data(si_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_y_0_xxxxyy = cbuffer.data(si_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_y_0_xxxxyz = cbuffer.data(si_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_y_0_xxxxzz = cbuffer.data(si_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_y_0_xxxyyy = cbuffer.data(si_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_y_0_xxxyyz = cbuffer.data(si_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_y_0_xxxyzz = cbuffer.data(si_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_y_0_xxxzzz = cbuffer.data(si_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_y_0_xxyyyy = cbuffer.data(si_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_y_0_xxyyyz = cbuffer.data(si_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_y_0_xxyyzz = cbuffer.data(si_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_y_0_xxyzzz = cbuffer.data(si_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_y_0_xxzzzz = cbuffer.data(si_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_y_0_xyyyyy = cbuffer.data(si_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_y_0_xyyyyz = cbuffer.data(si_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_0_xyyyzz = cbuffer.data(si_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_0_xyyzzz = cbuffer.data(si_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_0_xyzzzz = cbuffer.data(si_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_y_0_xzzzzz = cbuffer.data(si_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_0_yyyyyy = cbuffer.data(si_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_0_yyyyyz = cbuffer.data(si_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_0_yyyyzz = cbuffer.data(si_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_0_yyyzzz = cbuffer.data(si_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_0_yyzzzz = cbuffer.data(si_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_y_0_yzzzzz = cbuffer.data(si_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_0_zzzzzz = cbuffer.data(si_geom_11_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxxx, g_0_y_0_xxxxxxx, g_0_y_0_xxxxxxy, g_0_y_0_xxxxxxz, g_0_y_0_xxxxxy, g_0_y_0_xxxxxyy, g_0_y_0_xxxxxyz, g_0_y_0_xxxxxz, g_0_y_0_xxxxxzz, g_0_y_0_xxxxyy, g_0_y_0_xxxxyyy, g_0_y_0_xxxxyyz, g_0_y_0_xxxxyz, g_0_y_0_xxxxyzz, g_0_y_0_xxxxzz, g_0_y_0_xxxxzzz, g_0_y_0_xxxyyy, g_0_y_0_xxxyyyy, g_0_y_0_xxxyyyz, g_0_y_0_xxxyyz, g_0_y_0_xxxyyzz, g_0_y_0_xxxyzz, g_0_y_0_xxxyzzz, g_0_y_0_xxxzzz, g_0_y_0_xxxzzzz, g_0_y_0_xxyyyy, g_0_y_0_xxyyyyy, g_0_y_0_xxyyyyz, g_0_y_0_xxyyyz, g_0_y_0_xxyyyzz, g_0_y_0_xxyyzz, g_0_y_0_xxyyzzz, g_0_y_0_xxyzzz, g_0_y_0_xxyzzzz, g_0_y_0_xxzzzz, g_0_y_0_xxzzzzz, g_0_y_0_xyyyyy, g_0_y_0_xyyyyyy, g_0_y_0_xyyyyyz, g_0_y_0_xyyyyz, g_0_y_0_xyyyyzz, g_0_y_0_xyyyzz, g_0_y_0_xyyyzzz, g_0_y_0_xyyzzz, g_0_y_0_xyyzzzz, g_0_y_0_xyzzzz, g_0_y_0_xyzzzzz, g_0_y_0_xzzzzz, g_0_y_0_xzzzzzz, g_0_y_0_yyyyyy, g_0_y_0_yyyyyz, g_0_y_0_yyyyzz, g_0_y_0_yyyzzz, g_0_y_0_yyzzzz, g_0_y_0_yzzzzz, g_0_y_0_zzzzzz, g_x_y_0_xxxxxx, g_x_y_0_xxxxxy, g_x_y_0_xxxxxz, g_x_y_0_xxxxyy, g_x_y_0_xxxxyz, g_x_y_0_xxxxzz, g_x_y_0_xxxyyy, g_x_y_0_xxxyyz, g_x_y_0_xxxyzz, g_x_y_0_xxxzzz, g_x_y_0_xxyyyy, g_x_y_0_xxyyyz, g_x_y_0_xxyyzz, g_x_y_0_xxyzzz, g_x_y_0_xxzzzz, g_x_y_0_xyyyyy, g_x_y_0_xyyyyz, g_x_y_0_xyyyzz, g_x_y_0_xyyzzz, g_x_y_0_xyzzzz, g_x_y_0_xzzzzz, g_x_y_0_yyyyyy, g_x_y_0_yyyyyz, g_x_y_0_yyyyzz, g_x_y_0_yyyzzz, g_x_y_0_yyzzzz, g_x_y_0_yzzzzz, g_x_y_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_0_xxxxxx[k] = -g_0_y_0_xxxxxx[k] * ab_x + g_0_y_0_xxxxxxx[k];

                g_x_y_0_xxxxxy[k] = -g_0_y_0_xxxxxy[k] * ab_x + g_0_y_0_xxxxxxy[k];

                g_x_y_0_xxxxxz[k] = -g_0_y_0_xxxxxz[k] * ab_x + g_0_y_0_xxxxxxz[k];

                g_x_y_0_xxxxyy[k] = -g_0_y_0_xxxxyy[k] * ab_x + g_0_y_0_xxxxxyy[k];

                g_x_y_0_xxxxyz[k] = -g_0_y_0_xxxxyz[k] * ab_x + g_0_y_0_xxxxxyz[k];

                g_x_y_0_xxxxzz[k] = -g_0_y_0_xxxxzz[k] * ab_x + g_0_y_0_xxxxxzz[k];

                g_x_y_0_xxxyyy[k] = -g_0_y_0_xxxyyy[k] * ab_x + g_0_y_0_xxxxyyy[k];

                g_x_y_0_xxxyyz[k] = -g_0_y_0_xxxyyz[k] * ab_x + g_0_y_0_xxxxyyz[k];

                g_x_y_0_xxxyzz[k] = -g_0_y_0_xxxyzz[k] * ab_x + g_0_y_0_xxxxyzz[k];

                g_x_y_0_xxxzzz[k] = -g_0_y_0_xxxzzz[k] * ab_x + g_0_y_0_xxxxzzz[k];

                g_x_y_0_xxyyyy[k] = -g_0_y_0_xxyyyy[k] * ab_x + g_0_y_0_xxxyyyy[k];

                g_x_y_0_xxyyyz[k] = -g_0_y_0_xxyyyz[k] * ab_x + g_0_y_0_xxxyyyz[k];

                g_x_y_0_xxyyzz[k] = -g_0_y_0_xxyyzz[k] * ab_x + g_0_y_0_xxxyyzz[k];

                g_x_y_0_xxyzzz[k] = -g_0_y_0_xxyzzz[k] * ab_x + g_0_y_0_xxxyzzz[k];

                g_x_y_0_xxzzzz[k] = -g_0_y_0_xxzzzz[k] * ab_x + g_0_y_0_xxxzzzz[k];

                g_x_y_0_xyyyyy[k] = -g_0_y_0_xyyyyy[k] * ab_x + g_0_y_0_xxyyyyy[k];

                g_x_y_0_xyyyyz[k] = -g_0_y_0_xyyyyz[k] * ab_x + g_0_y_0_xxyyyyz[k];

                g_x_y_0_xyyyzz[k] = -g_0_y_0_xyyyzz[k] * ab_x + g_0_y_0_xxyyyzz[k];

                g_x_y_0_xyyzzz[k] = -g_0_y_0_xyyzzz[k] * ab_x + g_0_y_0_xxyyzzz[k];

                g_x_y_0_xyzzzz[k] = -g_0_y_0_xyzzzz[k] * ab_x + g_0_y_0_xxyzzzz[k];

                g_x_y_0_xzzzzz[k] = -g_0_y_0_xzzzzz[k] * ab_x + g_0_y_0_xxzzzzz[k];

                g_x_y_0_yyyyyy[k] = -g_0_y_0_yyyyyy[k] * ab_x + g_0_y_0_xyyyyyy[k];

                g_x_y_0_yyyyyz[k] = -g_0_y_0_yyyyyz[k] * ab_x + g_0_y_0_xyyyyyz[k];

                g_x_y_0_yyyyzz[k] = -g_0_y_0_yyyyzz[k] * ab_x + g_0_y_0_xyyyyzz[k];

                g_x_y_0_yyyzzz[k] = -g_0_y_0_yyyzzz[k] * ab_x + g_0_y_0_xyyyzzz[k];

                g_x_y_0_yyzzzz[k] = -g_0_y_0_yyzzzz[k] * ab_x + g_0_y_0_xyyzzzz[k];

                g_x_y_0_yzzzzz[k] = -g_0_y_0_yzzzzz[k] * ab_x + g_0_y_0_xyzzzzz[k];

                g_x_y_0_zzzzzz[k] = -g_0_y_0_zzzzzz[k] * ab_x + g_0_y_0_xzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_z_0_xxxxxx = cbuffer.data(si_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_z_0_xxxxxy = cbuffer.data(si_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_z_0_xxxxxz = cbuffer.data(si_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_z_0_xxxxyy = cbuffer.data(si_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_z_0_xxxxyz = cbuffer.data(si_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_z_0_xxxxzz = cbuffer.data(si_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_z_0_xxxyyy = cbuffer.data(si_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_z_0_xxxyyz = cbuffer.data(si_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_z_0_xxxyzz = cbuffer.data(si_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_z_0_xxxzzz = cbuffer.data(si_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_z_0_xxyyyy = cbuffer.data(si_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_z_0_xxyyyz = cbuffer.data(si_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_z_0_xxyyzz = cbuffer.data(si_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_z_0_xxyzzz = cbuffer.data(si_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_z_0_xxzzzz = cbuffer.data(si_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_z_0_xyyyyy = cbuffer.data(si_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_z_0_xyyyyz = cbuffer.data(si_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_z_0_xyyyzz = cbuffer.data(si_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_z_0_xyyzzz = cbuffer.data(si_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_z_0_xyzzzz = cbuffer.data(si_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_z_0_xzzzzz = cbuffer.data(si_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_z_0_yyyyyy = cbuffer.data(si_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_z_0_yyyyyz = cbuffer.data(si_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_z_0_yyyyzz = cbuffer.data(si_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_z_0_yyyzzz = cbuffer.data(si_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_z_0_yyzzzz = cbuffer.data(si_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_z_0_yzzzzz = cbuffer.data(si_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_z_0_zzzzzz = cbuffer.data(si_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxxx, g_0_z_0_xxxxxxx, g_0_z_0_xxxxxxy, g_0_z_0_xxxxxxz, g_0_z_0_xxxxxy, g_0_z_0_xxxxxyy, g_0_z_0_xxxxxyz, g_0_z_0_xxxxxz, g_0_z_0_xxxxxzz, g_0_z_0_xxxxyy, g_0_z_0_xxxxyyy, g_0_z_0_xxxxyyz, g_0_z_0_xxxxyz, g_0_z_0_xxxxyzz, g_0_z_0_xxxxzz, g_0_z_0_xxxxzzz, g_0_z_0_xxxyyy, g_0_z_0_xxxyyyy, g_0_z_0_xxxyyyz, g_0_z_0_xxxyyz, g_0_z_0_xxxyyzz, g_0_z_0_xxxyzz, g_0_z_0_xxxyzzz, g_0_z_0_xxxzzz, g_0_z_0_xxxzzzz, g_0_z_0_xxyyyy, g_0_z_0_xxyyyyy, g_0_z_0_xxyyyyz, g_0_z_0_xxyyyz, g_0_z_0_xxyyyzz, g_0_z_0_xxyyzz, g_0_z_0_xxyyzzz, g_0_z_0_xxyzzz, g_0_z_0_xxyzzzz, g_0_z_0_xxzzzz, g_0_z_0_xxzzzzz, g_0_z_0_xyyyyy, g_0_z_0_xyyyyyy, g_0_z_0_xyyyyyz, g_0_z_0_xyyyyz, g_0_z_0_xyyyyzz, g_0_z_0_xyyyzz, g_0_z_0_xyyyzzz, g_0_z_0_xyyzzz, g_0_z_0_xyyzzzz, g_0_z_0_xyzzzz, g_0_z_0_xyzzzzz, g_0_z_0_xzzzzz, g_0_z_0_xzzzzzz, g_0_z_0_yyyyyy, g_0_z_0_yyyyyz, g_0_z_0_yyyyzz, g_0_z_0_yyyzzz, g_0_z_0_yyzzzz, g_0_z_0_yzzzzz, g_0_z_0_zzzzzz, g_x_z_0_xxxxxx, g_x_z_0_xxxxxy, g_x_z_0_xxxxxz, g_x_z_0_xxxxyy, g_x_z_0_xxxxyz, g_x_z_0_xxxxzz, g_x_z_0_xxxyyy, g_x_z_0_xxxyyz, g_x_z_0_xxxyzz, g_x_z_0_xxxzzz, g_x_z_0_xxyyyy, g_x_z_0_xxyyyz, g_x_z_0_xxyyzz, g_x_z_0_xxyzzz, g_x_z_0_xxzzzz, g_x_z_0_xyyyyy, g_x_z_0_xyyyyz, g_x_z_0_xyyyzz, g_x_z_0_xyyzzz, g_x_z_0_xyzzzz, g_x_z_0_xzzzzz, g_x_z_0_yyyyyy, g_x_z_0_yyyyyz, g_x_z_0_yyyyzz, g_x_z_0_yyyzzz, g_x_z_0_yyzzzz, g_x_z_0_yzzzzz, g_x_z_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_0_xxxxxx[k] = -g_0_z_0_xxxxxx[k] * ab_x + g_0_z_0_xxxxxxx[k];

                g_x_z_0_xxxxxy[k] = -g_0_z_0_xxxxxy[k] * ab_x + g_0_z_0_xxxxxxy[k];

                g_x_z_0_xxxxxz[k] = -g_0_z_0_xxxxxz[k] * ab_x + g_0_z_0_xxxxxxz[k];

                g_x_z_0_xxxxyy[k] = -g_0_z_0_xxxxyy[k] * ab_x + g_0_z_0_xxxxxyy[k];

                g_x_z_0_xxxxyz[k] = -g_0_z_0_xxxxyz[k] * ab_x + g_0_z_0_xxxxxyz[k];

                g_x_z_0_xxxxzz[k] = -g_0_z_0_xxxxzz[k] * ab_x + g_0_z_0_xxxxxzz[k];

                g_x_z_0_xxxyyy[k] = -g_0_z_0_xxxyyy[k] * ab_x + g_0_z_0_xxxxyyy[k];

                g_x_z_0_xxxyyz[k] = -g_0_z_0_xxxyyz[k] * ab_x + g_0_z_0_xxxxyyz[k];

                g_x_z_0_xxxyzz[k] = -g_0_z_0_xxxyzz[k] * ab_x + g_0_z_0_xxxxyzz[k];

                g_x_z_0_xxxzzz[k] = -g_0_z_0_xxxzzz[k] * ab_x + g_0_z_0_xxxxzzz[k];

                g_x_z_0_xxyyyy[k] = -g_0_z_0_xxyyyy[k] * ab_x + g_0_z_0_xxxyyyy[k];

                g_x_z_0_xxyyyz[k] = -g_0_z_0_xxyyyz[k] * ab_x + g_0_z_0_xxxyyyz[k];

                g_x_z_0_xxyyzz[k] = -g_0_z_0_xxyyzz[k] * ab_x + g_0_z_0_xxxyyzz[k];

                g_x_z_0_xxyzzz[k] = -g_0_z_0_xxyzzz[k] * ab_x + g_0_z_0_xxxyzzz[k];

                g_x_z_0_xxzzzz[k] = -g_0_z_0_xxzzzz[k] * ab_x + g_0_z_0_xxxzzzz[k];

                g_x_z_0_xyyyyy[k] = -g_0_z_0_xyyyyy[k] * ab_x + g_0_z_0_xxyyyyy[k];

                g_x_z_0_xyyyyz[k] = -g_0_z_0_xyyyyz[k] * ab_x + g_0_z_0_xxyyyyz[k];

                g_x_z_0_xyyyzz[k] = -g_0_z_0_xyyyzz[k] * ab_x + g_0_z_0_xxyyyzz[k];

                g_x_z_0_xyyzzz[k] = -g_0_z_0_xyyzzz[k] * ab_x + g_0_z_0_xxyyzzz[k];

                g_x_z_0_xyzzzz[k] = -g_0_z_0_xyzzzz[k] * ab_x + g_0_z_0_xxyzzzz[k];

                g_x_z_0_xzzzzz[k] = -g_0_z_0_xzzzzz[k] * ab_x + g_0_z_0_xxzzzzz[k];

                g_x_z_0_yyyyyy[k] = -g_0_z_0_yyyyyy[k] * ab_x + g_0_z_0_xyyyyyy[k];

                g_x_z_0_yyyyyz[k] = -g_0_z_0_yyyyyz[k] * ab_x + g_0_z_0_xyyyyyz[k];

                g_x_z_0_yyyyzz[k] = -g_0_z_0_yyyyzz[k] * ab_x + g_0_z_0_xyyyyzz[k];

                g_x_z_0_yyyzzz[k] = -g_0_z_0_yyyzzz[k] * ab_x + g_0_z_0_xyyyzzz[k];

                g_x_z_0_yyzzzz[k] = -g_0_z_0_yyzzzz[k] * ab_x + g_0_z_0_xyyzzzz[k];

                g_x_z_0_yzzzzz[k] = -g_0_z_0_yzzzzz[k] * ab_x + g_0_z_0_xyzzzzz[k];

                g_x_z_0_zzzzzz[k] = -g_0_z_0_zzzzzz[k] * ab_x + g_0_z_0_xzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_y_x_0_xxxxxx = cbuffer.data(si_geom_11_off + 84 * ccomps * dcomps);

            auto g_y_x_0_xxxxxy = cbuffer.data(si_geom_11_off + 85 * ccomps * dcomps);

            auto g_y_x_0_xxxxxz = cbuffer.data(si_geom_11_off + 86 * ccomps * dcomps);

            auto g_y_x_0_xxxxyy = cbuffer.data(si_geom_11_off + 87 * ccomps * dcomps);

            auto g_y_x_0_xxxxyz = cbuffer.data(si_geom_11_off + 88 * ccomps * dcomps);

            auto g_y_x_0_xxxxzz = cbuffer.data(si_geom_11_off + 89 * ccomps * dcomps);

            auto g_y_x_0_xxxyyy = cbuffer.data(si_geom_11_off + 90 * ccomps * dcomps);

            auto g_y_x_0_xxxyyz = cbuffer.data(si_geom_11_off + 91 * ccomps * dcomps);

            auto g_y_x_0_xxxyzz = cbuffer.data(si_geom_11_off + 92 * ccomps * dcomps);

            auto g_y_x_0_xxxzzz = cbuffer.data(si_geom_11_off + 93 * ccomps * dcomps);

            auto g_y_x_0_xxyyyy = cbuffer.data(si_geom_11_off + 94 * ccomps * dcomps);

            auto g_y_x_0_xxyyyz = cbuffer.data(si_geom_11_off + 95 * ccomps * dcomps);

            auto g_y_x_0_xxyyzz = cbuffer.data(si_geom_11_off + 96 * ccomps * dcomps);

            auto g_y_x_0_xxyzzz = cbuffer.data(si_geom_11_off + 97 * ccomps * dcomps);

            auto g_y_x_0_xxzzzz = cbuffer.data(si_geom_11_off + 98 * ccomps * dcomps);

            auto g_y_x_0_xyyyyy = cbuffer.data(si_geom_11_off + 99 * ccomps * dcomps);

            auto g_y_x_0_xyyyyz = cbuffer.data(si_geom_11_off + 100 * ccomps * dcomps);

            auto g_y_x_0_xyyyzz = cbuffer.data(si_geom_11_off + 101 * ccomps * dcomps);

            auto g_y_x_0_xyyzzz = cbuffer.data(si_geom_11_off + 102 * ccomps * dcomps);

            auto g_y_x_0_xyzzzz = cbuffer.data(si_geom_11_off + 103 * ccomps * dcomps);

            auto g_y_x_0_xzzzzz = cbuffer.data(si_geom_11_off + 104 * ccomps * dcomps);

            auto g_y_x_0_yyyyyy = cbuffer.data(si_geom_11_off + 105 * ccomps * dcomps);

            auto g_y_x_0_yyyyyz = cbuffer.data(si_geom_11_off + 106 * ccomps * dcomps);

            auto g_y_x_0_yyyyzz = cbuffer.data(si_geom_11_off + 107 * ccomps * dcomps);

            auto g_y_x_0_yyyzzz = cbuffer.data(si_geom_11_off + 108 * ccomps * dcomps);

            auto g_y_x_0_yyzzzz = cbuffer.data(si_geom_11_off + 109 * ccomps * dcomps);

            auto g_y_x_0_yzzzzz = cbuffer.data(si_geom_11_off + 110 * ccomps * dcomps);

            auto g_y_x_0_zzzzzz = cbuffer.data(si_geom_11_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxx, g_0_x_0_xxxxxxy, g_0_x_0_xxxxxy, g_0_x_0_xxxxxyy, g_0_x_0_xxxxxyz, g_0_x_0_xxxxxz, g_0_x_0_xxxxyy, g_0_x_0_xxxxyyy, g_0_x_0_xxxxyyz, g_0_x_0_xxxxyz, g_0_x_0_xxxxyzz, g_0_x_0_xxxxzz, g_0_x_0_xxxyyy, g_0_x_0_xxxyyyy, g_0_x_0_xxxyyyz, g_0_x_0_xxxyyz, g_0_x_0_xxxyyzz, g_0_x_0_xxxyzz, g_0_x_0_xxxyzzz, g_0_x_0_xxxzzz, g_0_x_0_xxyyyy, g_0_x_0_xxyyyyy, g_0_x_0_xxyyyyz, g_0_x_0_xxyyyz, g_0_x_0_xxyyyzz, g_0_x_0_xxyyzz, g_0_x_0_xxyyzzz, g_0_x_0_xxyzzz, g_0_x_0_xxyzzzz, g_0_x_0_xxzzzz, g_0_x_0_xyyyyy, g_0_x_0_xyyyyyy, g_0_x_0_xyyyyyz, g_0_x_0_xyyyyz, g_0_x_0_xyyyyzz, g_0_x_0_xyyyzz, g_0_x_0_xyyyzzz, g_0_x_0_xyyzzz, g_0_x_0_xyyzzzz, g_0_x_0_xyzzzz, g_0_x_0_xyzzzzz, g_0_x_0_xzzzzz, g_0_x_0_yyyyyy, g_0_x_0_yyyyyyy, g_0_x_0_yyyyyyz, g_0_x_0_yyyyyz, g_0_x_0_yyyyyzz, g_0_x_0_yyyyzz, g_0_x_0_yyyyzzz, g_0_x_0_yyyzzz, g_0_x_0_yyyzzzz, g_0_x_0_yyzzzz, g_0_x_0_yyzzzzz, g_0_x_0_yzzzzz, g_0_x_0_yzzzzzz, g_0_x_0_zzzzzz, g_y_x_0_xxxxxx, g_y_x_0_xxxxxy, g_y_x_0_xxxxxz, g_y_x_0_xxxxyy, g_y_x_0_xxxxyz, g_y_x_0_xxxxzz, g_y_x_0_xxxyyy, g_y_x_0_xxxyyz, g_y_x_0_xxxyzz, g_y_x_0_xxxzzz, g_y_x_0_xxyyyy, g_y_x_0_xxyyyz, g_y_x_0_xxyyzz, g_y_x_0_xxyzzz, g_y_x_0_xxzzzz, g_y_x_0_xyyyyy, g_y_x_0_xyyyyz, g_y_x_0_xyyyzz, g_y_x_0_xyyzzz, g_y_x_0_xyzzzz, g_y_x_0_xzzzzz, g_y_x_0_yyyyyy, g_y_x_0_yyyyyz, g_y_x_0_yyyyzz, g_y_x_0_yyyzzz, g_y_x_0_yyzzzz, g_y_x_0_yzzzzz, g_y_x_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_0_xxxxxx[k] = -g_0_x_0_xxxxxx[k] * ab_y + g_0_x_0_xxxxxxy[k];

                g_y_x_0_xxxxxy[k] = -g_0_x_0_xxxxxy[k] * ab_y + g_0_x_0_xxxxxyy[k];

                g_y_x_0_xxxxxz[k] = -g_0_x_0_xxxxxz[k] * ab_y + g_0_x_0_xxxxxyz[k];

                g_y_x_0_xxxxyy[k] = -g_0_x_0_xxxxyy[k] * ab_y + g_0_x_0_xxxxyyy[k];

                g_y_x_0_xxxxyz[k] = -g_0_x_0_xxxxyz[k] * ab_y + g_0_x_0_xxxxyyz[k];

                g_y_x_0_xxxxzz[k] = -g_0_x_0_xxxxzz[k] * ab_y + g_0_x_0_xxxxyzz[k];

                g_y_x_0_xxxyyy[k] = -g_0_x_0_xxxyyy[k] * ab_y + g_0_x_0_xxxyyyy[k];

                g_y_x_0_xxxyyz[k] = -g_0_x_0_xxxyyz[k] * ab_y + g_0_x_0_xxxyyyz[k];

                g_y_x_0_xxxyzz[k] = -g_0_x_0_xxxyzz[k] * ab_y + g_0_x_0_xxxyyzz[k];

                g_y_x_0_xxxzzz[k] = -g_0_x_0_xxxzzz[k] * ab_y + g_0_x_0_xxxyzzz[k];

                g_y_x_0_xxyyyy[k] = -g_0_x_0_xxyyyy[k] * ab_y + g_0_x_0_xxyyyyy[k];

                g_y_x_0_xxyyyz[k] = -g_0_x_0_xxyyyz[k] * ab_y + g_0_x_0_xxyyyyz[k];

                g_y_x_0_xxyyzz[k] = -g_0_x_0_xxyyzz[k] * ab_y + g_0_x_0_xxyyyzz[k];

                g_y_x_0_xxyzzz[k] = -g_0_x_0_xxyzzz[k] * ab_y + g_0_x_0_xxyyzzz[k];

                g_y_x_0_xxzzzz[k] = -g_0_x_0_xxzzzz[k] * ab_y + g_0_x_0_xxyzzzz[k];

                g_y_x_0_xyyyyy[k] = -g_0_x_0_xyyyyy[k] * ab_y + g_0_x_0_xyyyyyy[k];

                g_y_x_0_xyyyyz[k] = -g_0_x_0_xyyyyz[k] * ab_y + g_0_x_0_xyyyyyz[k];

                g_y_x_0_xyyyzz[k] = -g_0_x_0_xyyyzz[k] * ab_y + g_0_x_0_xyyyyzz[k];

                g_y_x_0_xyyzzz[k] = -g_0_x_0_xyyzzz[k] * ab_y + g_0_x_0_xyyyzzz[k];

                g_y_x_0_xyzzzz[k] = -g_0_x_0_xyzzzz[k] * ab_y + g_0_x_0_xyyzzzz[k];

                g_y_x_0_xzzzzz[k] = -g_0_x_0_xzzzzz[k] * ab_y + g_0_x_0_xyzzzzz[k];

                g_y_x_0_yyyyyy[k] = -g_0_x_0_yyyyyy[k] * ab_y + g_0_x_0_yyyyyyy[k];

                g_y_x_0_yyyyyz[k] = -g_0_x_0_yyyyyz[k] * ab_y + g_0_x_0_yyyyyyz[k];

                g_y_x_0_yyyyzz[k] = -g_0_x_0_yyyyzz[k] * ab_y + g_0_x_0_yyyyyzz[k];

                g_y_x_0_yyyzzz[k] = -g_0_x_0_yyyzzz[k] * ab_y + g_0_x_0_yyyyzzz[k];

                g_y_x_0_yyzzzz[k] = -g_0_x_0_yyzzzz[k] * ab_y + g_0_x_0_yyyzzzz[k];

                g_y_x_0_yzzzzz[k] = -g_0_x_0_yzzzzz[k] * ab_y + g_0_x_0_yyzzzzz[k];

                g_y_x_0_zzzzzz[k] = -g_0_x_0_zzzzzz[k] * ab_y + g_0_x_0_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_y_y_0_xxxxxx = cbuffer.data(si_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_y_0_xxxxxy = cbuffer.data(si_geom_11_off + 113 * ccomps * dcomps);

            auto g_y_y_0_xxxxxz = cbuffer.data(si_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_y_0_xxxxyy = cbuffer.data(si_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_y_0_xxxxyz = cbuffer.data(si_geom_11_off + 116 * ccomps * dcomps);

            auto g_y_y_0_xxxxzz = cbuffer.data(si_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_y_0_xxxyyy = cbuffer.data(si_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_y_0_xxxyyz = cbuffer.data(si_geom_11_off + 119 * ccomps * dcomps);

            auto g_y_y_0_xxxyzz = cbuffer.data(si_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_y_0_xxxzzz = cbuffer.data(si_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_y_0_xxyyyy = cbuffer.data(si_geom_11_off + 122 * ccomps * dcomps);

            auto g_y_y_0_xxyyyz = cbuffer.data(si_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_y_0_xxyyzz = cbuffer.data(si_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_y_0_xxyzzz = cbuffer.data(si_geom_11_off + 125 * ccomps * dcomps);

            auto g_y_y_0_xxzzzz = cbuffer.data(si_geom_11_off + 126 * ccomps * dcomps);

            auto g_y_y_0_xyyyyy = cbuffer.data(si_geom_11_off + 127 * ccomps * dcomps);

            auto g_y_y_0_xyyyyz = cbuffer.data(si_geom_11_off + 128 * ccomps * dcomps);

            auto g_y_y_0_xyyyzz = cbuffer.data(si_geom_11_off + 129 * ccomps * dcomps);

            auto g_y_y_0_xyyzzz = cbuffer.data(si_geom_11_off + 130 * ccomps * dcomps);

            auto g_y_y_0_xyzzzz = cbuffer.data(si_geom_11_off + 131 * ccomps * dcomps);

            auto g_y_y_0_xzzzzz = cbuffer.data(si_geom_11_off + 132 * ccomps * dcomps);

            auto g_y_y_0_yyyyyy = cbuffer.data(si_geom_11_off + 133 * ccomps * dcomps);

            auto g_y_y_0_yyyyyz = cbuffer.data(si_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_y_0_yyyyzz = cbuffer.data(si_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_y_0_yyyzzz = cbuffer.data(si_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_y_0_yyzzzz = cbuffer.data(si_geom_11_off + 137 * ccomps * dcomps);

            auto g_y_y_0_yzzzzz = cbuffer.data(si_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_y_0_zzzzzz = cbuffer.data(si_geom_11_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyzzz, g_0_xyzzzz, g_0_xzzzzz, g_0_y_0_xxxxxx, g_0_y_0_xxxxxxy, g_0_y_0_xxxxxy, g_0_y_0_xxxxxyy, g_0_y_0_xxxxxyz, g_0_y_0_xxxxxz, g_0_y_0_xxxxyy, g_0_y_0_xxxxyyy, g_0_y_0_xxxxyyz, g_0_y_0_xxxxyz, g_0_y_0_xxxxyzz, g_0_y_0_xxxxzz, g_0_y_0_xxxyyy, g_0_y_0_xxxyyyy, g_0_y_0_xxxyyyz, g_0_y_0_xxxyyz, g_0_y_0_xxxyyzz, g_0_y_0_xxxyzz, g_0_y_0_xxxyzzz, g_0_y_0_xxxzzz, g_0_y_0_xxyyyy, g_0_y_0_xxyyyyy, g_0_y_0_xxyyyyz, g_0_y_0_xxyyyz, g_0_y_0_xxyyyzz, g_0_y_0_xxyyzz, g_0_y_0_xxyyzzz, g_0_y_0_xxyzzz, g_0_y_0_xxyzzzz, g_0_y_0_xxzzzz, g_0_y_0_xyyyyy, g_0_y_0_xyyyyyy, g_0_y_0_xyyyyyz, g_0_y_0_xyyyyz, g_0_y_0_xyyyyzz, g_0_y_0_xyyyzz, g_0_y_0_xyyyzzz, g_0_y_0_xyyzzz, g_0_y_0_xyyzzzz, g_0_y_0_xyzzzz, g_0_y_0_xyzzzzz, g_0_y_0_xzzzzz, g_0_y_0_yyyyyy, g_0_y_0_yyyyyyy, g_0_y_0_yyyyyyz, g_0_y_0_yyyyyz, g_0_y_0_yyyyyzz, g_0_y_0_yyyyzz, g_0_y_0_yyyyzzz, g_0_y_0_yyyzzz, g_0_y_0_yyyzzzz, g_0_y_0_yyzzzz, g_0_y_0_yyzzzzz, g_0_y_0_yzzzzz, g_0_y_0_yzzzzzz, g_0_y_0_zzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_zzzzzz, g_y_y_0_xxxxxx, g_y_y_0_xxxxxy, g_y_y_0_xxxxxz, g_y_y_0_xxxxyy, g_y_y_0_xxxxyz, g_y_y_0_xxxxzz, g_y_y_0_xxxyyy, g_y_y_0_xxxyyz, g_y_y_0_xxxyzz, g_y_y_0_xxxzzz, g_y_y_0_xxyyyy, g_y_y_0_xxyyyz, g_y_y_0_xxyyzz, g_y_y_0_xxyzzz, g_y_y_0_xxzzzz, g_y_y_0_xyyyyy, g_y_y_0_xyyyyz, g_y_y_0_xyyyzz, g_y_y_0_xyyzzz, g_y_y_0_xyzzzz, g_y_y_0_xzzzzz, g_y_y_0_yyyyyy, g_y_y_0_yyyyyz, g_y_y_0_yyyyzz, g_y_y_0_yyyzzz, g_y_y_0_yyzzzz, g_y_y_0_yzzzzz, g_y_y_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_0_xxxxxx[k] = g_0_xxxxxx[k] - g_0_y_0_xxxxxx[k] * ab_y + g_0_y_0_xxxxxxy[k];

                g_y_y_0_xxxxxy[k] = g_0_xxxxxy[k] - g_0_y_0_xxxxxy[k] * ab_y + g_0_y_0_xxxxxyy[k];

                g_y_y_0_xxxxxz[k] = g_0_xxxxxz[k] - g_0_y_0_xxxxxz[k] * ab_y + g_0_y_0_xxxxxyz[k];

                g_y_y_0_xxxxyy[k] = g_0_xxxxyy[k] - g_0_y_0_xxxxyy[k] * ab_y + g_0_y_0_xxxxyyy[k];

                g_y_y_0_xxxxyz[k] = g_0_xxxxyz[k] - g_0_y_0_xxxxyz[k] * ab_y + g_0_y_0_xxxxyyz[k];

                g_y_y_0_xxxxzz[k] = g_0_xxxxzz[k] - g_0_y_0_xxxxzz[k] * ab_y + g_0_y_0_xxxxyzz[k];

                g_y_y_0_xxxyyy[k] = g_0_xxxyyy[k] - g_0_y_0_xxxyyy[k] * ab_y + g_0_y_0_xxxyyyy[k];

                g_y_y_0_xxxyyz[k] = g_0_xxxyyz[k] - g_0_y_0_xxxyyz[k] * ab_y + g_0_y_0_xxxyyyz[k];

                g_y_y_0_xxxyzz[k] = g_0_xxxyzz[k] - g_0_y_0_xxxyzz[k] * ab_y + g_0_y_0_xxxyyzz[k];

                g_y_y_0_xxxzzz[k] = g_0_xxxzzz[k] - g_0_y_0_xxxzzz[k] * ab_y + g_0_y_0_xxxyzzz[k];

                g_y_y_0_xxyyyy[k] = g_0_xxyyyy[k] - g_0_y_0_xxyyyy[k] * ab_y + g_0_y_0_xxyyyyy[k];

                g_y_y_0_xxyyyz[k] = g_0_xxyyyz[k] - g_0_y_0_xxyyyz[k] * ab_y + g_0_y_0_xxyyyyz[k];

                g_y_y_0_xxyyzz[k] = g_0_xxyyzz[k] - g_0_y_0_xxyyzz[k] * ab_y + g_0_y_0_xxyyyzz[k];

                g_y_y_0_xxyzzz[k] = g_0_xxyzzz[k] - g_0_y_0_xxyzzz[k] * ab_y + g_0_y_0_xxyyzzz[k];

                g_y_y_0_xxzzzz[k] = g_0_xxzzzz[k] - g_0_y_0_xxzzzz[k] * ab_y + g_0_y_0_xxyzzzz[k];

                g_y_y_0_xyyyyy[k] = g_0_xyyyyy[k] - g_0_y_0_xyyyyy[k] * ab_y + g_0_y_0_xyyyyyy[k];

                g_y_y_0_xyyyyz[k] = g_0_xyyyyz[k] - g_0_y_0_xyyyyz[k] * ab_y + g_0_y_0_xyyyyyz[k];

                g_y_y_0_xyyyzz[k] = g_0_xyyyzz[k] - g_0_y_0_xyyyzz[k] * ab_y + g_0_y_0_xyyyyzz[k];

                g_y_y_0_xyyzzz[k] = g_0_xyyzzz[k] - g_0_y_0_xyyzzz[k] * ab_y + g_0_y_0_xyyyzzz[k];

                g_y_y_0_xyzzzz[k] = g_0_xyzzzz[k] - g_0_y_0_xyzzzz[k] * ab_y + g_0_y_0_xyyzzzz[k];

                g_y_y_0_xzzzzz[k] = g_0_xzzzzz[k] - g_0_y_0_xzzzzz[k] * ab_y + g_0_y_0_xyzzzzz[k];

                g_y_y_0_yyyyyy[k] = g_0_yyyyyy[k] - g_0_y_0_yyyyyy[k] * ab_y + g_0_y_0_yyyyyyy[k];

                g_y_y_0_yyyyyz[k] = g_0_yyyyyz[k] - g_0_y_0_yyyyyz[k] * ab_y + g_0_y_0_yyyyyyz[k];

                g_y_y_0_yyyyzz[k] = g_0_yyyyzz[k] - g_0_y_0_yyyyzz[k] * ab_y + g_0_y_0_yyyyyzz[k];

                g_y_y_0_yyyzzz[k] = g_0_yyyzzz[k] - g_0_y_0_yyyzzz[k] * ab_y + g_0_y_0_yyyyzzz[k];

                g_y_y_0_yyzzzz[k] = g_0_yyzzzz[k] - g_0_y_0_yyzzzz[k] * ab_y + g_0_y_0_yyyzzzz[k];

                g_y_y_0_yzzzzz[k] = g_0_yzzzzz[k] - g_0_y_0_yzzzzz[k] * ab_y + g_0_y_0_yyzzzzz[k];

                g_y_y_0_zzzzzz[k] = g_0_zzzzzz[k] - g_0_y_0_zzzzzz[k] * ab_y + g_0_y_0_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_y_z_0_xxxxxx = cbuffer.data(si_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_z_0_xxxxxy = cbuffer.data(si_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_z_0_xxxxxz = cbuffer.data(si_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_z_0_xxxxyy = cbuffer.data(si_geom_11_off + 143 * ccomps * dcomps);

            auto g_y_z_0_xxxxyz = cbuffer.data(si_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_z_0_xxxxzz = cbuffer.data(si_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_z_0_xxxyyy = cbuffer.data(si_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_z_0_xxxyyz = cbuffer.data(si_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_z_0_xxxyzz = cbuffer.data(si_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_z_0_xxxzzz = cbuffer.data(si_geom_11_off + 149 * ccomps * dcomps);

            auto g_y_z_0_xxyyyy = cbuffer.data(si_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_z_0_xxyyyz = cbuffer.data(si_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_z_0_xxyyzz = cbuffer.data(si_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_z_0_xxyzzz = cbuffer.data(si_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_z_0_xxzzzz = cbuffer.data(si_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_z_0_xyyyyy = cbuffer.data(si_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_z_0_xyyyyz = cbuffer.data(si_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_z_0_xyyyzz = cbuffer.data(si_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_z_0_xyyzzz = cbuffer.data(si_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_z_0_xyzzzz = cbuffer.data(si_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_z_0_xzzzzz = cbuffer.data(si_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_z_0_yyyyyy = cbuffer.data(si_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_z_0_yyyyyz = cbuffer.data(si_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_z_0_yyyyzz = cbuffer.data(si_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_z_0_yyyzzz = cbuffer.data(si_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_z_0_yyzzzz = cbuffer.data(si_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_z_0_yzzzzz = cbuffer.data(si_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_z_0_zzzzzz = cbuffer.data(si_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxxx, g_0_z_0_xxxxxxy, g_0_z_0_xxxxxy, g_0_z_0_xxxxxyy, g_0_z_0_xxxxxyz, g_0_z_0_xxxxxz, g_0_z_0_xxxxyy, g_0_z_0_xxxxyyy, g_0_z_0_xxxxyyz, g_0_z_0_xxxxyz, g_0_z_0_xxxxyzz, g_0_z_0_xxxxzz, g_0_z_0_xxxyyy, g_0_z_0_xxxyyyy, g_0_z_0_xxxyyyz, g_0_z_0_xxxyyz, g_0_z_0_xxxyyzz, g_0_z_0_xxxyzz, g_0_z_0_xxxyzzz, g_0_z_0_xxxzzz, g_0_z_0_xxyyyy, g_0_z_0_xxyyyyy, g_0_z_0_xxyyyyz, g_0_z_0_xxyyyz, g_0_z_0_xxyyyzz, g_0_z_0_xxyyzz, g_0_z_0_xxyyzzz, g_0_z_0_xxyzzz, g_0_z_0_xxyzzzz, g_0_z_0_xxzzzz, g_0_z_0_xyyyyy, g_0_z_0_xyyyyyy, g_0_z_0_xyyyyyz, g_0_z_0_xyyyyz, g_0_z_0_xyyyyzz, g_0_z_0_xyyyzz, g_0_z_0_xyyyzzz, g_0_z_0_xyyzzz, g_0_z_0_xyyzzzz, g_0_z_0_xyzzzz, g_0_z_0_xyzzzzz, g_0_z_0_xzzzzz, g_0_z_0_yyyyyy, g_0_z_0_yyyyyyy, g_0_z_0_yyyyyyz, g_0_z_0_yyyyyz, g_0_z_0_yyyyyzz, g_0_z_0_yyyyzz, g_0_z_0_yyyyzzz, g_0_z_0_yyyzzz, g_0_z_0_yyyzzzz, g_0_z_0_yyzzzz, g_0_z_0_yyzzzzz, g_0_z_0_yzzzzz, g_0_z_0_yzzzzzz, g_0_z_0_zzzzzz, g_y_z_0_xxxxxx, g_y_z_0_xxxxxy, g_y_z_0_xxxxxz, g_y_z_0_xxxxyy, g_y_z_0_xxxxyz, g_y_z_0_xxxxzz, g_y_z_0_xxxyyy, g_y_z_0_xxxyyz, g_y_z_0_xxxyzz, g_y_z_0_xxxzzz, g_y_z_0_xxyyyy, g_y_z_0_xxyyyz, g_y_z_0_xxyyzz, g_y_z_0_xxyzzz, g_y_z_0_xxzzzz, g_y_z_0_xyyyyy, g_y_z_0_xyyyyz, g_y_z_0_xyyyzz, g_y_z_0_xyyzzz, g_y_z_0_xyzzzz, g_y_z_0_xzzzzz, g_y_z_0_yyyyyy, g_y_z_0_yyyyyz, g_y_z_0_yyyyzz, g_y_z_0_yyyzzz, g_y_z_0_yyzzzz, g_y_z_0_yzzzzz, g_y_z_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_0_xxxxxx[k] = -g_0_z_0_xxxxxx[k] * ab_y + g_0_z_0_xxxxxxy[k];

                g_y_z_0_xxxxxy[k] = -g_0_z_0_xxxxxy[k] * ab_y + g_0_z_0_xxxxxyy[k];

                g_y_z_0_xxxxxz[k] = -g_0_z_0_xxxxxz[k] * ab_y + g_0_z_0_xxxxxyz[k];

                g_y_z_0_xxxxyy[k] = -g_0_z_0_xxxxyy[k] * ab_y + g_0_z_0_xxxxyyy[k];

                g_y_z_0_xxxxyz[k] = -g_0_z_0_xxxxyz[k] * ab_y + g_0_z_0_xxxxyyz[k];

                g_y_z_0_xxxxzz[k] = -g_0_z_0_xxxxzz[k] * ab_y + g_0_z_0_xxxxyzz[k];

                g_y_z_0_xxxyyy[k] = -g_0_z_0_xxxyyy[k] * ab_y + g_0_z_0_xxxyyyy[k];

                g_y_z_0_xxxyyz[k] = -g_0_z_0_xxxyyz[k] * ab_y + g_0_z_0_xxxyyyz[k];

                g_y_z_0_xxxyzz[k] = -g_0_z_0_xxxyzz[k] * ab_y + g_0_z_0_xxxyyzz[k];

                g_y_z_0_xxxzzz[k] = -g_0_z_0_xxxzzz[k] * ab_y + g_0_z_0_xxxyzzz[k];

                g_y_z_0_xxyyyy[k] = -g_0_z_0_xxyyyy[k] * ab_y + g_0_z_0_xxyyyyy[k];

                g_y_z_0_xxyyyz[k] = -g_0_z_0_xxyyyz[k] * ab_y + g_0_z_0_xxyyyyz[k];

                g_y_z_0_xxyyzz[k] = -g_0_z_0_xxyyzz[k] * ab_y + g_0_z_0_xxyyyzz[k];

                g_y_z_0_xxyzzz[k] = -g_0_z_0_xxyzzz[k] * ab_y + g_0_z_0_xxyyzzz[k];

                g_y_z_0_xxzzzz[k] = -g_0_z_0_xxzzzz[k] * ab_y + g_0_z_0_xxyzzzz[k];

                g_y_z_0_xyyyyy[k] = -g_0_z_0_xyyyyy[k] * ab_y + g_0_z_0_xyyyyyy[k];

                g_y_z_0_xyyyyz[k] = -g_0_z_0_xyyyyz[k] * ab_y + g_0_z_0_xyyyyyz[k];

                g_y_z_0_xyyyzz[k] = -g_0_z_0_xyyyzz[k] * ab_y + g_0_z_0_xyyyyzz[k];

                g_y_z_0_xyyzzz[k] = -g_0_z_0_xyyzzz[k] * ab_y + g_0_z_0_xyyyzzz[k];

                g_y_z_0_xyzzzz[k] = -g_0_z_0_xyzzzz[k] * ab_y + g_0_z_0_xyyzzzz[k];

                g_y_z_0_xzzzzz[k] = -g_0_z_0_xzzzzz[k] * ab_y + g_0_z_0_xyzzzzz[k];

                g_y_z_0_yyyyyy[k] = -g_0_z_0_yyyyyy[k] * ab_y + g_0_z_0_yyyyyyy[k];

                g_y_z_0_yyyyyz[k] = -g_0_z_0_yyyyyz[k] * ab_y + g_0_z_0_yyyyyyz[k];

                g_y_z_0_yyyyzz[k] = -g_0_z_0_yyyyzz[k] * ab_y + g_0_z_0_yyyyyzz[k];

                g_y_z_0_yyyzzz[k] = -g_0_z_0_yyyzzz[k] * ab_y + g_0_z_0_yyyyzzz[k];

                g_y_z_0_yyzzzz[k] = -g_0_z_0_yyzzzz[k] * ab_y + g_0_z_0_yyyzzzz[k];

                g_y_z_0_yzzzzz[k] = -g_0_z_0_yzzzzz[k] * ab_y + g_0_z_0_yyzzzzz[k];

                g_y_z_0_zzzzzz[k] = -g_0_z_0_zzzzzz[k] * ab_y + g_0_z_0_yzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_z_x_0_xxxxxx = cbuffer.data(si_geom_11_off + 168 * ccomps * dcomps);

            auto g_z_x_0_xxxxxy = cbuffer.data(si_geom_11_off + 169 * ccomps * dcomps);

            auto g_z_x_0_xxxxxz = cbuffer.data(si_geom_11_off + 170 * ccomps * dcomps);

            auto g_z_x_0_xxxxyy = cbuffer.data(si_geom_11_off + 171 * ccomps * dcomps);

            auto g_z_x_0_xxxxyz = cbuffer.data(si_geom_11_off + 172 * ccomps * dcomps);

            auto g_z_x_0_xxxxzz = cbuffer.data(si_geom_11_off + 173 * ccomps * dcomps);

            auto g_z_x_0_xxxyyy = cbuffer.data(si_geom_11_off + 174 * ccomps * dcomps);

            auto g_z_x_0_xxxyyz = cbuffer.data(si_geom_11_off + 175 * ccomps * dcomps);

            auto g_z_x_0_xxxyzz = cbuffer.data(si_geom_11_off + 176 * ccomps * dcomps);

            auto g_z_x_0_xxxzzz = cbuffer.data(si_geom_11_off + 177 * ccomps * dcomps);

            auto g_z_x_0_xxyyyy = cbuffer.data(si_geom_11_off + 178 * ccomps * dcomps);

            auto g_z_x_0_xxyyyz = cbuffer.data(si_geom_11_off + 179 * ccomps * dcomps);

            auto g_z_x_0_xxyyzz = cbuffer.data(si_geom_11_off + 180 * ccomps * dcomps);

            auto g_z_x_0_xxyzzz = cbuffer.data(si_geom_11_off + 181 * ccomps * dcomps);

            auto g_z_x_0_xxzzzz = cbuffer.data(si_geom_11_off + 182 * ccomps * dcomps);

            auto g_z_x_0_xyyyyy = cbuffer.data(si_geom_11_off + 183 * ccomps * dcomps);

            auto g_z_x_0_xyyyyz = cbuffer.data(si_geom_11_off + 184 * ccomps * dcomps);

            auto g_z_x_0_xyyyzz = cbuffer.data(si_geom_11_off + 185 * ccomps * dcomps);

            auto g_z_x_0_xyyzzz = cbuffer.data(si_geom_11_off + 186 * ccomps * dcomps);

            auto g_z_x_0_xyzzzz = cbuffer.data(si_geom_11_off + 187 * ccomps * dcomps);

            auto g_z_x_0_xzzzzz = cbuffer.data(si_geom_11_off + 188 * ccomps * dcomps);

            auto g_z_x_0_yyyyyy = cbuffer.data(si_geom_11_off + 189 * ccomps * dcomps);

            auto g_z_x_0_yyyyyz = cbuffer.data(si_geom_11_off + 190 * ccomps * dcomps);

            auto g_z_x_0_yyyyzz = cbuffer.data(si_geom_11_off + 191 * ccomps * dcomps);

            auto g_z_x_0_yyyzzz = cbuffer.data(si_geom_11_off + 192 * ccomps * dcomps);

            auto g_z_x_0_yyzzzz = cbuffer.data(si_geom_11_off + 193 * ccomps * dcomps);

            auto g_z_x_0_yzzzzz = cbuffer.data(si_geom_11_off + 194 * ccomps * dcomps);

            auto g_z_x_0_zzzzzz = cbuffer.data(si_geom_11_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxx, g_0_x_0_xxxxxxz, g_0_x_0_xxxxxy, g_0_x_0_xxxxxyz, g_0_x_0_xxxxxz, g_0_x_0_xxxxxzz, g_0_x_0_xxxxyy, g_0_x_0_xxxxyyz, g_0_x_0_xxxxyz, g_0_x_0_xxxxyzz, g_0_x_0_xxxxzz, g_0_x_0_xxxxzzz, g_0_x_0_xxxyyy, g_0_x_0_xxxyyyz, g_0_x_0_xxxyyz, g_0_x_0_xxxyyzz, g_0_x_0_xxxyzz, g_0_x_0_xxxyzzz, g_0_x_0_xxxzzz, g_0_x_0_xxxzzzz, g_0_x_0_xxyyyy, g_0_x_0_xxyyyyz, g_0_x_0_xxyyyz, g_0_x_0_xxyyyzz, g_0_x_0_xxyyzz, g_0_x_0_xxyyzzz, g_0_x_0_xxyzzz, g_0_x_0_xxyzzzz, g_0_x_0_xxzzzz, g_0_x_0_xxzzzzz, g_0_x_0_xyyyyy, g_0_x_0_xyyyyyz, g_0_x_0_xyyyyz, g_0_x_0_xyyyyzz, g_0_x_0_xyyyzz, g_0_x_0_xyyyzzz, g_0_x_0_xyyzzz, g_0_x_0_xyyzzzz, g_0_x_0_xyzzzz, g_0_x_0_xyzzzzz, g_0_x_0_xzzzzz, g_0_x_0_xzzzzzz, g_0_x_0_yyyyyy, g_0_x_0_yyyyyyz, g_0_x_0_yyyyyz, g_0_x_0_yyyyyzz, g_0_x_0_yyyyzz, g_0_x_0_yyyyzzz, g_0_x_0_yyyzzz, g_0_x_0_yyyzzzz, g_0_x_0_yyzzzz, g_0_x_0_yyzzzzz, g_0_x_0_yzzzzz, g_0_x_0_yzzzzzz, g_0_x_0_zzzzzz, g_0_x_0_zzzzzzz, g_z_x_0_xxxxxx, g_z_x_0_xxxxxy, g_z_x_0_xxxxxz, g_z_x_0_xxxxyy, g_z_x_0_xxxxyz, g_z_x_0_xxxxzz, g_z_x_0_xxxyyy, g_z_x_0_xxxyyz, g_z_x_0_xxxyzz, g_z_x_0_xxxzzz, g_z_x_0_xxyyyy, g_z_x_0_xxyyyz, g_z_x_0_xxyyzz, g_z_x_0_xxyzzz, g_z_x_0_xxzzzz, g_z_x_0_xyyyyy, g_z_x_0_xyyyyz, g_z_x_0_xyyyzz, g_z_x_0_xyyzzz, g_z_x_0_xyzzzz, g_z_x_0_xzzzzz, g_z_x_0_yyyyyy, g_z_x_0_yyyyyz, g_z_x_0_yyyyzz, g_z_x_0_yyyzzz, g_z_x_0_yyzzzz, g_z_x_0_yzzzzz, g_z_x_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_0_xxxxxx[k] = -g_0_x_0_xxxxxx[k] * ab_z + g_0_x_0_xxxxxxz[k];

                g_z_x_0_xxxxxy[k] = -g_0_x_0_xxxxxy[k] * ab_z + g_0_x_0_xxxxxyz[k];

                g_z_x_0_xxxxxz[k] = -g_0_x_0_xxxxxz[k] * ab_z + g_0_x_0_xxxxxzz[k];

                g_z_x_0_xxxxyy[k] = -g_0_x_0_xxxxyy[k] * ab_z + g_0_x_0_xxxxyyz[k];

                g_z_x_0_xxxxyz[k] = -g_0_x_0_xxxxyz[k] * ab_z + g_0_x_0_xxxxyzz[k];

                g_z_x_0_xxxxzz[k] = -g_0_x_0_xxxxzz[k] * ab_z + g_0_x_0_xxxxzzz[k];

                g_z_x_0_xxxyyy[k] = -g_0_x_0_xxxyyy[k] * ab_z + g_0_x_0_xxxyyyz[k];

                g_z_x_0_xxxyyz[k] = -g_0_x_0_xxxyyz[k] * ab_z + g_0_x_0_xxxyyzz[k];

                g_z_x_0_xxxyzz[k] = -g_0_x_0_xxxyzz[k] * ab_z + g_0_x_0_xxxyzzz[k];

                g_z_x_0_xxxzzz[k] = -g_0_x_0_xxxzzz[k] * ab_z + g_0_x_0_xxxzzzz[k];

                g_z_x_0_xxyyyy[k] = -g_0_x_0_xxyyyy[k] * ab_z + g_0_x_0_xxyyyyz[k];

                g_z_x_0_xxyyyz[k] = -g_0_x_0_xxyyyz[k] * ab_z + g_0_x_0_xxyyyzz[k];

                g_z_x_0_xxyyzz[k] = -g_0_x_0_xxyyzz[k] * ab_z + g_0_x_0_xxyyzzz[k];

                g_z_x_0_xxyzzz[k] = -g_0_x_0_xxyzzz[k] * ab_z + g_0_x_0_xxyzzzz[k];

                g_z_x_0_xxzzzz[k] = -g_0_x_0_xxzzzz[k] * ab_z + g_0_x_0_xxzzzzz[k];

                g_z_x_0_xyyyyy[k] = -g_0_x_0_xyyyyy[k] * ab_z + g_0_x_0_xyyyyyz[k];

                g_z_x_0_xyyyyz[k] = -g_0_x_0_xyyyyz[k] * ab_z + g_0_x_0_xyyyyzz[k];

                g_z_x_0_xyyyzz[k] = -g_0_x_0_xyyyzz[k] * ab_z + g_0_x_0_xyyyzzz[k];

                g_z_x_0_xyyzzz[k] = -g_0_x_0_xyyzzz[k] * ab_z + g_0_x_0_xyyzzzz[k];

                g_z_x_0_xyzzzz[k] = -g_0_x_0_xyzzzz[k] * ab_z + g_0_x_0_xyzzzzz[k];

                g_z_x_0_xzzzzz[k] = -g_0_x_0_xzzzzz[k] * ab_z + g_0_x_0_xzzzzzz[k];

                g_z_x_0_yyyyyy[k] = -g_0_x_0_yyyyyy[k] * ab_z + g_0_x_0_yyyyyyz[k];

                g_z_x_0_yyyyyz[k] = -g_0_x_0_yyyyyz[k] * ab_z + g_0_x_0_yyyyyzz[k];

                g_z_x_0_yyyyzz[k] = -g_0_x_0_yyyyzz[k] * ab_z + g_0_x_0_yyyyzzz[k];

                g_z_x_0_yyyzzz[k] = -g_0_x_0_yyyzzz[k] * ab_z + g_0_x_0_yyyzzzz[k];

                g_z_x_0_yyzzzz[k] = -g_0_x_0_yyzzzz[k] * ab_z + g_0_x_0_yyzzzzz[k];

                g_z_x_0_yzzzzz[k] = -g_0_x_0_yzzzzz[k] * ab_z + g_0_x_0_yzzzzzz[k];

                g_z_x_0_zzzzzz[k] = -g_0_x_0_zzzzzz[k] * ab_z + g_0_x_0_zzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_z_y_0_xxxxxx = cbuffer.data(si_geom_11_off + 196 * ccomps * dcomps);

            auto g_z_y_0_xxxxxy = cbuffer.data(si_geom_11_off + 197 * ccomps * dcomps);

            auto g_z_y_0_xxxxxz = cbuffer.data(si_geom_11_off + 198 * ccomps * dcomps);

            auto g_z_y_0_xxxxyy = cbuffer.data(si_geom_11_off + 199 * ccomps * dcomps);

            auto g_z_y_0_xxxxyz = cbuffer.data(si_geom_11_off + 200 * ccomps * dcomps);

            auto g_z_y_0_xxxxzz = cbuffer.data(si_geom_11_off + 201 * ccomps * dcomps);

            auto g_z_y_0_xxxyyy = cbuffer.data(si_geom_11_off + 202 * ccomps * dcomps);

            auto g_z_y_0_xxxyyz = cbuffer.data(si_geom_11_off + 203 * ccomps * dcomps);

            auto g_z_y_0_xxxyzz = cbuffer.data(si_geom_11_off + 204 * ccomps * dcomps);

            auto g_z_y_0_xxxzzz = cbuffer.data(si_geom_11_off + 205 * ccomps * dcomps);

            auto g_z_y_0_xxyyyy = cbuffer.data(si_geom_11_off + 206 * ccomps * dcomps);

            auto g_z_y_0_xxyyyz = cbuffer.data(si_geom_11_off + 207 * ccomps * dcomps);

            auto g_z_y_0_xxyyzz = cbuffer.data(si_geom_11_off + 208 * ccomps * dcomps);

            auto g_z_y_0_xxyzzz = cbuffer.data(si_geom_11_off + 209 * ccomps * dcomps);

            auto g_z_y_0_xxzzzz = cbuffer.data(si_geom_11_off + 210 * ccomps * dcomps);

            auto g_z_y_0_xyyyyy = cbuffer.data(si_geom_11_off + 211 * ccomps * dcomps);

            auto g_z_y_0_xyyyyz = cbuffer.data(si_geom_11_off + 212 * ccomps * dcomps);

            auto g_z_y_0_xyyyzz = cbuffer.data(si_geom_11_off + 213 * ccomps * dcomps);

            auto g_z_y_0_xyyzzz = cbuffer.data(si_geom_11_off + 214 * ccomps * dcomps);

            auto g_z_y_0_xyzzzz = cbuffer.data(si_geom_11_off + 215 * ccomps * dcomps);

            auto g_z_y_0_xzzzzz = cbuffer.data(si_geom_11_off + 216 * ccomps * dcomps);

            auto g_z_y_0_yyyyyy = cbuffer.data(si_geom_11_off + 217 * ccomps * dcomps);

            auto g_z_y_0_yyyyyz = cbuffer.data(si_geom_11_off + 218 * ccomps * dcomps);

            auto g_z_y_0_yyyyzz = cbuffer.data(si_geom_11_off + 219 * ccomps * dcomps);

            auto g_z_y_0_yyyzzz = cbuffer.data(si_geom_11_off + 220 * ccomps * dcomps);

            auto g_z_y_0_yyzzzz = cbuffer.data(si_geom_11_off + 221 * ccomps * dcomps);

            auto g_z_y_0_yzzzzz = cbuffer.data(si_geom_11_off + 222 * ccomps * dcomps);

            auto g_z_y_0_zzzzzz = cbuffer.data(si_geom_11_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxxx, g_0_y_0_xxxxxxz, g_0_y_0_xxxxxy, g_0_y_0_xxxxxyz, g_0_y_0_xxxxxz, g_0_y_0_xxxxxzz, g_0_y_0_xxxxyy, g_0_y_0_xxxxyyz, g_0_y_0_xxxxyz, g_0_y_0_xxxxyzz, g_0_y_0_xxxxzz, g_0_y_0_xxxxzzz, g_0_y_0_xxxyyy, g_0_y_0_xxxyyyz, g_0_y_0_xxxyyz, g_0_y_0_xxxyyzz, g_0_y_0_xxxyzz, g_0_y_0_xxxyzzz, g_0_y_0_xxxzzz, g_0_y_0_xxxzzzz, g_0_y_0_xxyyyy, g_0_y_0_xxyyyyz, g_0_y_0_xxyyyz, g_0_y_0_xxyyyzz, g_0_y_0_xxyyzz, g_0_y_0_xxyyzzz, g_0_y_0_xxyzzz, g_0_y_0_xxyzzzz, g_0_y_0_xxzzzz, g_0_y_0_xxzzzzz, g_0_y_0_xyyyyy, g_0_y_0_xyyyyyz, g_0_y_0_xyyyyz, g_0_y_0_xyyyyzz, g_0_y_0_xyyyzz, g_0_y_0_xyyyzzz, g_0_y_0_xyyzzz, g_0_y_0_xyyzzzz, g_0_y_0_xyzzzz, g_0_y_0_xyzzzzz, g_0_y_0_xzzzzz, g_0_y_0_xzzzzzz, g_0_y_0_yyyyyy, g_0_y_0_yyyyyyz, g_0_y_0_yyyyyz, g_0_y_0_yyyyyzz, g_0_y_0_yyyyzz, g_0_y_0_yyyyzzz, g_0_y_0_yyyzzz, g_0_y_0_yyyzzzz, g_0_y_0_yyzzzz, g_0_y_0_yyzzzzz, g_0_y_0_yzzzzz, g_0_y_0_yzzzzzz, g_0_y_0_zzzzzz, g_0_y_0_zzzzzzz, g_z_y_0_xxxxxx, g_z_y_0_xxxxxy, g_z_y_0_xxxxxz, g_z_y_0_xxxxyy, g_z_y_0_xxxxyz, g_z_y_0_xxxxzz, g_z_y_0_xxxyyy, g_z_y_0_xxxyyz, g_z_y_0_xxxyzz, g_z_y_0_xxxzzz, g_z_y_0_xxyyyy, g_z_y_0_xxyyyz, g_z_y_0_xxyyzz, g_z_y_0_xxyzzz, g_z_y_0_xxzzzz, g_z_y_0_xyyyyy, g_z_y_0_xyyyyz, g_z_y_0_xyyyzz, g_z_y_0_xyyzzz, g_z_y_0_xyzzzz, g_z_y_0_xzzzzz, g_z_y_0_yyyyyy, g_z_y_0_yyyyyz, g_z_y_0_yyyyzz, g_z_y_0_yyyzzz, g_z_y_0_yyzzzz, g_z_y_0_yzzzzz, g_z_y_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_0_xxxxxx[k] = -g_0_y_0_xxxxxx[k] * ab_z + g_0_y_0_xxxxxxz[k];

                g_z_y_0_xxxxxy[k] = -g_0_y_0_xxxxxy[k] * ab_z + g_0_y_0_xxxxxyz[k];

                g_z_y_0_xxxxxz[k] = -g_0_y_0_xxxxxz[k] * ab_z + g_0_y_0_xxxxxzz[k];

                g_z_y_0_xxxxyy[k] = -g_0_y_0_xxxxyy[k] * ab_z + g_0_y_0_xxxxyyz[k];

                g_z_y_0_xxxxyz[k] = -g_0_y_0_xxxxyz[k] * ab_z + g_0_y_0_xxxxyzz[k];

                g_z_y_0_xxxxzz[k] = -g_0_y_0_xxxxzz[k] * ab_z + g_0_y_0_xxxxzzz[k];

                g_z_y_0_xxxyyy[k] = -g_0_y_0_xxxyyy[k] * ab_z + g_0_y_0_xxxyyyz[k];

                g_z_y_0_xxxyyz[k] = -g_0_y_0_xxxyyz[k] * ab_z + g_0_y_0_xxxyyzz[k];

                g_z_y_0_xxxyzz[k] = -g_0_y_0_xxxyzz[k] * ab_z + g_0_y_0_xxxyzzz[k];

                g_z_y_0_xxxzzz[k] = -g_0_y_0_xxxzzz[k] * ab_z + g_0_y_0_xxxzzzz[k];

                g_z_y_0_xxyyyy[k] = -g_0_y_0_xxyyyy[k] * ab_z + g_0_y_0_xxyyyyz[k];

                g_z_y_0_xxyyyz[k] = -g_0_y_0_xxyyyz[k] * ab_z + g_0_y_0_xxyyyzz[k];

                g_z_y_0_xxyyzz[k] = -g_0_y_0_xxyyzz[k] * ab_z + g_0_y_0_xxyyzzz[k];

                g_z_y_0_xxyzzz[k] = -g_0_y_0_xxyzzz[k] * ab_z + g_0_y_0_xxyzzzz[k];

                g_z_y_0_xxzzzz[k] = -g_0_y_0_xxzzzz[k] * ab_z + g_0_y_0_xxzzzzz[k];

                g_z_y_0_xyyyyy[k] = -g_0_y_0_xyyyyy[k] * ab_z + g_0_y_0_xyyyyyz[k];

                g_z_y_0_xyyyyz[k] = -g_0_y_0_xyyyyz[k] * ab_z + g_0_y_0_xyyyyzz[k];

                g_z_y_0_xyyyzz[k] = -g_0_y_0_xyyyzz[k] * ab_z + g_0_y_0_xyyyzzz[k];

                g_z_y_0_xyyzzz[k] = -g_0_y_0_xyyzzz[k] * ab_z + g_0_y_0_xyyzzzz[k];

                g_z_y_0_xyzzzz[k] = -g_0_y_0_xyzzzz[k] * ab_z + g_0_y_0_xyzzzzz[k];

                g_z_y_0_xzzzzz[k] = -g_0_y_0_xzzzzz[k] * ab_z + g_0_y_0_xzzzzzz[k];

                g_z_y_0_yyyyyy[k] = -g_0_y_0_yyyyyy[k] * ab_z + g_0_y_0_yyyyyyz[k];

                g_z_y_0_yyyyyz[k] = -g_0_y_0_yyyyyz[k] * ab_z + g_0_y_0_yyyyyzz[k];

                g_z_y_0_yyyyzz[k] = -g_0_y_0_yyyyzz[k] * ab_z + g_0_y_0_yyyyzzz[k];

                g_z_y_0_yyyzzz[k] = -g_0_y_0_yyyzzz[k] * ab_z + g_0_y_0_yyyzzzz[k];

                g_z_y_0_yyzzzz[k] = -g_0_y_0_yyzzzz[k] * ab_z + g_0_y_0_yyzzzzz[k];

                g_z_y_0_yzzzzz[k] = -g_0_y_0_yzzzzz[k] * ab_z + g_0_y_0_yzzzzzz[k];

                g_z_y_0_zzzzzz[k] = -g_0_y_0_zzzzzz[k] * ab_z + g_0_y_0_zzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_z_z_0_xxxxxx = cbuffer.data(si_geom_11_off + 224 * ccomps * dcomps);

            auto g_z_z_0_xxxxxy = cbuffer.data(si_geom_11_off + 225 * ccomps * dcomps);

            auto g_z_z_0_xxxxxz = cbuffer.data(si_geom_11_off + 226 * ccomps * dcomps);

            auto g_z_z_0_xxxxyy = cbuffer.data(si_geom_11_off + 227 * ccomps * dcomps);

            auto g_z_z_0_xxxxyz = cbuffer.data(si_geom_11_off + 228 * ccomps * dcomps);

            auto g_z_z_0_xxxxzz = cbuffer.data(si_geom_11_off + 229 * ccomps * dcomps);

            auto g_z_z_0_xxxyyy = cbuffer.data(si_geom_11_off + 230 * ccomps * dcomps);

            auto g_z_z_0_xxxyyz = cbuffer.data(si_geom_11_off + 231 * ccomps * dcomps);

            auto g_z_z_0_xxxyzz = cbuffer.data(si_geom_11_off + 232 * ccomps * dcomps);

            auto g_z_z_0_xxxzzz = cbuffer.data(si_geom_11_off + 233 * ccomps * dcomps);

            auto g_z_z_0_xxyyyy = cbuffer.data(si_geom_11_off + 234 * ccomps * dcomps);

            auto g_z_z_0_xxyyyz = cbuffer.data(si_geom_11_off + 235 * ccomps * dcomps);

            auto g_z_z_0_xxyyzz = cbuffer.data(si_geom_11_off + 236 * ccomps * dcomps);

            auto g_z_z_0_xxyzzz = cbuffer.data(si_geom_11_off + 237 * ccomps * dcomps);

            auto g_z_z_0_xxzzzz = cbuffer.data(si_geom_11_off + 238 * ccomps * dcomps);

            auto g_z_z_0_xyyyyy = cbuffer.data(si_geom_11_off + 239 * ccomps * dcomps);

            auto g_z_z_0_xyyyyz = cbuffer.data(si_geom_11_off + 240 * ccomps * dcomps);

            auto g_z_z_0_xyyyzz = cbuffer.data(si_geom_11_off + 241 * ccomps * dcomps);

            auto g_z_z_0_xyyzzz = cbuffer.data(si_geom_11_off + 242 * ccomps * dcomps);

            auto g_z_z_0_xyzzzz = cbuffer.data(si_geom_11_off + 243 * ccomps * dcomps);

            auto g_z_z_0_xzzzzz = cbuffer.data(si_geom_11_off + 244 * ccomps * dcomps);

            auto g_z_z_0_yyyyyy = cbuffer.data(si_geom_11_off + 245 * ccomps * dcomps);

            auto g_z_z_0_yyyyyz = cbuffer.data(si_geom_11_off + 246 * ccomps * dcomps);

            auto g_z_z_0_yyyyzz = cbuffer.data(si_geom_11_off + 247 * ccomps * dcomps);

            auto g_z_z_0_yyyzzz = cbuffer.data(si_geom_11_off + 248 * ccomps * dcomps);

            auto g_z_z_0_yyzzzz = cbuffer.data(si_geom_11_off + 249 * ccomps * dcomps);

            auto g_z_z_0_yzzzzz = cbuffer.data(si_geom_11_off + 250 * ccomps * dcomps);

            auto g_z_z_0_zzzzzz = cbuffer.data(si_geom_11_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyzzz, g_0_xyzzzz, g_0_xzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_z_0_xxxxxx, g_0_z_0_xxxxxxz, g_0_z_0_xxxxxy, g_0_z_0_xxxxxyz, g_0_z_0_xxxxxz, g_0_z_0_xxxxxzz, g_0_z_0_xxxxyy, g_0_z_0_xxxxyyz, g_0_z_0_xxxxyz, g_0_z_0_xxxxyzz, g_0_z_0_xxxxzz, g_0_z_0_xxxxzzz, g_0_z_0_xxxyyy, g_0_z_0_xxxyyyz, g_0_z_0_xxxyyz, g_0_z_0_xxxyyzz, g_0_z_0_xxxyzz, g_0_z_0_xxxyzzz, g_0_z_0_xxxzzz, g_0_z_0_xxxzzzz, g_0_z_0_xxyyyy, g_0_z_0_xxyyyyz, g_0_z_0_xxyyyz, g_0_z_0_xxyyyzz, g_0_z_0_xxyyzz, g_0_z_0_xxyyzzz, g_0_z_0_xxyzzz, g_0_z_0_xxyzzzz, g_0_z_0_xxzzzz, g_0_z_0_xxzzzzz, g_0_z_0_xyyyyy, g_0_z_0_xyyyyyz, g_0_z_0_xyyyyz, g_0_z_0_xyyyyzz, g_0_z_0_xyyyzz, g_0_z_0_xyyyzzz, g_0_z_0_xyyzzz, g_0_z_0_xyyzzzz, g_0_z_0_xyzzzz, g_0_z_0_xyzzzzz, g_0_z_0_xzzzzz, g_0_z_0_xzzzzzz, g_0_z_0_yyyyyy, g_0_z_0_yyyyyyz, g_0_z_0_yyyyyz, g_0_z_0_yyyyyzz, g_0_z_0_yyyyzz, g_0_z_0_yyyyzzz, g_0_z_0_yyyzzz, g_0_z_0_yyyzzzz, g_0_z_0_yyzzzz, g_0_z_0_yyzzzzz, g_0_z_0_yzzzzz, g_0_z_0_yzzzzzz, g_0_z_0_zzzzzz, g_0_z_0_zzzzzzz, g_0_zzzzzz, g_z_z_0_xxxxxx, g_z_z_0_xxxxxy, g_z_z_0_xxxxxz, g_z_z_0_xxxxyy, g_z_z_0_xxxxyz, g_z_z_0_xxxxzz, g_z_z_0_xxxyyy, g_z_z_0_xxxyyz, g_z_z_0_xxxyzz, g_z_z_0_xxxzzz, g_z_z_0_xxyyyy, g_z_z_0_xxyyyz, g_z_z_0_xxyyzz, g_z_z_0_xxyzzz, g_z_z_0_xxzzzz, g_z_z_0_xyyyyy, g_z_z_0_xyyyyz, g_z_z_0_xyyyzz, g_z_z_0_xyyzzz, g_z_z_0_xyzzzz, g_z_z_0_xzzzzz, g_z_z_0_yyyyyy, g_z_z_0_yyyyyz, g_z_z_0_yyyyzz, g_z_z_0_yyyzzz, g_z_z_0_yyzzzz, g_z_z_0_yzzzzz, g_z_z_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_0_xxxxxx[k] = g_0_xxxxxx[k] - g_0_z_0_xxxxxx[k] * ab_z + g_0_z_0_xxxxxxz[k];

                g_z_z_0_xxxxxy[k] = g_0_xxxxxy[k] - g_0_z_0_xxxxxy[k] * ab_z + g_0_z_0_xxxxxyz[k];

                g_z_z_0_xxxxxz[k] = g_0_xxxxxz[k] - g_0_z_0_xxxxxz[k] * ab_z + g_0_z_0_xxxxxzz[k];

                g_z_z_0_xxxxyy[k] = g_0_xxxxyy[k] - g_0_z_0_xxxxyy[k] * ab_z + g_0_z_0_xxxxyyz[k];

                g_z_z_0_xxxxyz[k] = g_0_xxxxyz[k] - g_0_z_0_xxxxyz[k] * ab_z + g_0_z_0_xxxxyzz[k];

                g_z_z_0_xxxxzz[k] = g_0_xxxxzz[k] - g_0_z_0_xxxxzz[k] * ab_z + g_0_z_0_xxxxzzz[k];

                g_z_z_0_xxxyyy[k] = g_0_xxxyyy[k] - g_0_z_0_xxxyyy[k] * ab_z + g_0_z_0_xxxyyyz[k];

                g_z_z_0_xxxyyz[k] = g_0_xxxyyz[k] - g_0_z_0_xxxyyz[k] * ab_z + g_0_z_0_xxxyyzz[k];

                g_z_z_0_xxxyzz[k] = g_0_xxxyzz[k] - g_0_z_0_xxxyzz[k] * ab_z + g_0_z_0_xxxyzzz[k];

                g_z_z_0_xxxzzz[k] = g_0_xxxzzz[k] - g_0_z_0_xxxzzz[k] * ab_z + g_0_z_0_xxxzzzz[k];

                g_z_z_0_xxyyyy[k] = g_0_xxyyyy[k] - g_0_z_0_xxyyyy[k] * ab_z + g_0_z_0_xxyyyyz[k];

                g_z_z_0_xxyyyz[k] = g_0_xxyyyz[k] - g_0_z_0_xxyyyz[k] * ab_z + g_0_z_0_xxyyyzz[k];

                g_z_z_0_xxyyzz[k] = g_0_xxyyzz[k] - g_0_z_0_xxyyzz[k] * ab_z + g_0_z_0_xxyyzzz[k];

                g_z_z_0_xxyzzz[k] = g_0_xxyzzz[k] - g_0_z_0_xxyzzz[k] * ab_z + g_0_z_0_xxyzzzz[k];

                g_z_z_0_xxzzzz[k] = g_0_xxzzzz[k] - g_0_z_0_xxzzzz[k] * ab_z + g_0_z_0_xxzzzzz[k];

                g_z_z_0_xyyyyy[k] = g_0_xyyyyy[k] - g_0_z_0_xyyyyy[k] * ab_z + g_0_z_0_xyyyyyz[k];

                g_z_z_0_xyyyyz[k] = g_0_xyyyyz[k] - g_0_z_0_xyyyyz[k] * ab_z + g_0_z_0_xyyyyzz[k];

                g_z_z_0_xyyyzz[k] = g_0_xyyyzz[k] - g_0_z_0_xyyyzz[k] * ab_z + g_0_z_0_xyyyzzz[k];

                g_z_z_0_xyyzzz[k] = g_0_xyyzzz[k] - g_0_z_0_xyyzzz[k] * ab_z + g_0_z_0_xyyzzzz[k];

                g_z_z_0_xyzzzz[k] = g_0_xyzzzz[k] - g_0_z_0_xyzzzz[k] * ab_z + g_0_z_0_xyzzzzz[k];

                g_z_z_0_xzzzzz[k] = g_0_xzzzzz[k] - g_0_z_0_xzzzzz[k] * ab_z + g_0_z_0_xzzzzzz[k];

                g_z_z_0_yyyyyy[k] = g_0_yyyyyy[k] - g_0_z_0_yyyyyy[k] * ab_z + g_0_z_0_yyyyyyz[k];

                g_z_z_0_yyyyyz[k] = g_0_yyyyyz[k] - g_0_z_0_yyyyyz[k] * ab_z + g_0_z_0_yyyyyzz[k];

                g_z_z_0_yyyyzz[k] = g_0_yyyyzz[k] - g_0_z_0_yyyyzz[k] * ab_z + g_0_z_0_yyyyzzz[k];

                g_z_z_0_yyyzzz[k] = g_0_yyyzzz[k] - g_0_z_0_yyyzzz[k] * ab_z + g_0_z_0_yyyzzzz[k];

                g_z_z_0_yyzzzz[k] = g_0_yyzzzz[k] - g_0_z_0_yyzzzz[k] * ab_z + g_0_z_0_yyzzzzz[k];

                g_z_z_0_yzzzzz[k] = g_0_yzzzzz[k] - g_0_z_0_yzzzzz[k] * ab_z + g_0_z_0_yzzzzzz[k];

                g_z_z_0_zzzzzz[k] = g_0_zzzzzz[k] - g_0_z_0_zzzzzz[k] * ab_z + g_0_z_0_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

