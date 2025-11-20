#include "ElectronRepulsionGeom0100ContrRecPIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_pixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_pixx,
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

            /// set up bra offset for contr_buffer_pixx

            const auto pi_geom_01_off = idx_geom_01_pixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_0_x_x_xxxxxx = cbuffer.data(pi_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xxxxxy = cbuffer.data(pi_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xxxxxz = cbuffer.data(pi_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_xxxxyy = cbuffer.data(pi_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_xxxxyz = cbuffer.data(pi_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_xxxxzz = cbuffer.data(pi_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_x_xxxyyy = cbuffer.data(pi_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_x_xxxyyz = cbuffer.data(pi_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_x_xxxyzz = cbuffer.data(pi_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_x_xxxzzz = cbuffer.data(pi_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_x_xxyyyy = cbuffer.data(pi_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_x_xxyyyz = cbuffer.data(pi_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_x_xxyyzz = cbuffer.data(pi_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_x_xxyzzz = cbuffer.data(pi_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_x_xxzzzz = cbuffer.data(pi_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_x_xyyyyy = cbuffer.data(pi_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_x_xyyyyz = cbuffer.data(pi_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_x_xyyyzz = cbuffer.data(pi_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_x_xyyzzz = cbuffer.data(pi_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_x_xyzzzz = cbuffer.data(pi_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_x_xzzzzz = cbuffer.data(pi_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_x_yyyyyy = cbuffer.data(pi_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_x_yyyyyz = cbuffer.data(pi_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_x_yyyyzz = cbuffer.data(pi_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_x_yyyzzz = cbuffer.data(pi_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_x_yyzzzz = cbuffer.data(pi_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_x_yzzzzz = cbuffer.data(pi_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_x_zzzzzz = cbuffer.data(pi_geom_01_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxx, g_0_x_0_xxxxxxx, g_0_x_0_xxxxxxy, g_0_x_0_xxxxxxz, g_0_x_0_xxxxxy, g_0_x_0_xxxxxyy, g_0_x_0_xxxxxyz, g_0_x_0_xxxxxz, g_0_x_0_xxxxxzz, g_0_x_0_xxxxyy, g_0_x_0_xxxxyyy, g_0_x_0_xxxxyyz, g_0_x_0_xxxxyz, g_0_x_0_xxxxyzz, g_0_x_0_xxxxzz, g_0_x_0_xxxxzzz, g_0_x_0_xxxyyy, g_0_x_0_xxxyyyy, g_0_x_0_xxxyyyz, g_0_x_0_xxxyyz, g_0_x_0_xxxyyzz, g_0_x_0_xxxyzz, g_0_x_0_xxxyzzz, g_0_x_0_xxxzzz, g_0_x_0_xxxzzzz, g_0_x_0_xxyyyy, g_0_x_0_xxyyyyy, g_0_x_0_xxyyyyz, g_0_x_0_xxyyyz, g_0_x_0_xxyyyzz, g_0_x_0_xxyyzz, g_0_x_0_xxyyzzz, g_0_x_0_xxyzzz, g_0_x_0_xxyzzzz, g_0_x_0_xxzzzz, g_0_x_0_xxzzzzz, g_0_x_0_xyyyyy, g_0_x_0_xyyyyyy, g_0_x_0_xyyyyyz, g_0_x_0_xyyyyz, g_0_x_0_xyyyyzz, g_0_x_0_xyyyzz, g_0_x_0_xyyyzzz, g_0_x_0_xyyzzz, g_0_x_0_xyyzzzz, g_0_x_0_xyzzzz, g_0_x_0_xyzzzzz, g_0_x_0_xzzzzz, g_0_x_0_xzzzzzz, g_0_x_0_yyyyyy, g_0_x_0_yyyyyz, g_0_x_0_yyyyzz, g_0_x_0_yyyzzz, g_0_x_0_yyzzzz, g_0_x_0_yzzzzz, g_0_x_0_zzzzzz, g_0_x_x_xxxxxx, g_0_x_x_xxxxxy, g_0_x_x_xxxxxz, g_0_x_x_xxxxyy, g_0_x_x_xxxxyz, g_0_x_x_xxxxzz, g_0_x_x_xxxyyy, g_0_x_x_xxxyyz, g_0_x_x_xxxyzz, g_0_x_x_xxxzzz, g_0_x_x_xxyyyy, g_0_x_x_xxyyyz, g_0_x_x_xxyyzz, g_0_x_x_xxyzzz, g_0_x_x_xxzzzz, g_0_x_x_xyyyyy, g_0_x_x_xyyyyz, g_0_x_x_xyyyzz, g_0_x_x_xyyzzz, g_0_x_x_xyzzzz, g_0_x_x_xzzzzz, g_0_x_x_yyyyyy, g_0_x_x_yyyyyz, g_0_x_x_yyyyzz, g_0_x_x_yyyzzz, g_0_x_x_yyzzzz, g_0_x_x_yzzzzz, g_0_x_x_zzzzzz, g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyzzz, g_0_xyzzzz, g_0_xzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_x_xxxxxx[k] = g_0_xxxxxx[k] - g_0_x_0_xxxxxx[k] * ab_x + g_0_x_0_xxxxxxx[k];

                g_0_x_x_xxxxxy[k] = g_0_xxxxxy[k] - g_0_x_0_xxxxxy[k] * ab_x + g_0_x_0_xxxxxxy[k];

                g_0_x_x_xxxxxz[k] = g_0_xxxxxz[k] - g_0_x_0_xxxxxz[k] * ab_x + g_0_x_0_xxxxxxz[k];

                g_0_x_x_xxxxyy[k] = g_0_xxxxyy[k] - g_0_x_0_xxxxyy[k] * ab_x + g_0_x_0_xxxxxyy[k];

                g_0_x_x_xxxxyz[k] = g_0_xxxxyz[k] - g_0_x_0_xxxxyz[k] * ab_x + g_0_x_0_xxxxxyz[k];

                g_0_x_x_xxxxzz[k] = g_0_xxxxzz[k] - g_0_x_0_xxxxzz[k] * ab_x + g_0_x_0_xxxxxzz[k];

                g_0_x_x_xxxyyy[k] = g_0_xxxyyy[k] - g_0_x_0_xxxyyy[k] * ab_x + g_0_x_0_xxxxyyy[k];

                g_0_x_x_xxxyyz[k] = g_0_xxxyyz[k] - g_0_x_0_xxxyyz[k] * ab_x + g_0_x_0_xxxxyyz[k];

                g_0_x_x_xxxyzz[k] = g_0_xxxyzz[k] - g_0_x_0_xxxyzz[k] * ab_x + g_0_x_0_xxxxyzz[k];

                g_0_x_x_xxxzzz[k] = g_0_xxxzzz[k] - g_0_x_0_xxxzzz[k] * ab_x + g_0_x_0_xxxxzzz[k];

                g_0_x_x_xxyyyy[k] = g_0_xxyyyy[k] - g_0_x_0_xxyyyy[k] * ab_x + g_0_x_0_xxxyyyy[k];

                g_0_x_x_xxyyyz[k] = g_0_xxyyyz[k] - g_0_x_0_xxyyyz[k] * ab_x + g_0_x_0_xxxyyyz[k];

                g_0_x_x_xxyyzz[k] = g_0_xxyyzz[k] - g_0_x_0_xxyyzz[k] * ab_x + g_0_x_0_xxxyyzz[k];

                g_0_x_x_xxyzzz[k] = g_0_xxyzzz[k] - g_0_x_0_xxyzzz[k] * ab_x + g_0_x_0_xxxyzzz[k];

                g_0_x_x_xxzzzz[k] = g_0_xxzzzz[k] - g_0_x_0_xxzzzz[k] * ab_x + g_0_x_0_xxxzzzz[k];

                g_0_x_x_xyyyyy[k] = g_0_xyyyyy[k] - g_0_x_0_xyyyyy[k] * ab_x + g_0_x_0_xxyyyyy[k];

                g_0_x_x_xyyyyz[k] = g_0_xyyyyz[k] - g_0_x_0_xyyyyz[k] * ab_x + g_0_x_0_xxyyyyz[k];

                g_0_x_x_xyyyzz[k] = g_0_xyyyzz[k] - g_0_x_0_xyyyzz[k] * ab_x + g_0_x_0_xxyyyzz[k];

                g_0_x_x_xyyzzz[k] = g_0_xyyzzz[k] - g_0_x_0_xyyzzz[k] * ab_x + g_0_x_0_xxyyzzz[k];

                g_0_x_x_xyzzzz[k] = g_0_xyzzzz[k] - g_0_x_0_xyzzzz[k] * ab_x + g_0_x_0_xxyzzzz[k];

                g_0_x_x_xzzzzz[k] = g_0_xzzzzz[k] - g_0_x_0_xzzzzz[k] * ab_x + g_0_x_0_xxzzzzz[k];

                g_0_x_x_yyyyyy[k] = g_0_yyyyyy[k] - g_0_x_0_yyyyyy[k] * ab_x + g_0_x_0_xyyyyyy[k];

                g_0_x_x_yyyyyz[k] = g_0_yyyyyz[k] - g_0_x_0_yyyyyz[k] * ab_x + g_0_x_0_xyyyyyz[k];

                g_0_x_x_yyyyzz[k] = g_0_yyyyzz[k] - g_0_x_0_yyyyzz[k] * ab_x + g_0_x_0_xyyyyzz[k];

                g_0_x_x_yyyzzz[k] = g_0_yyyzzz[k] - g_0_x_0_yyyzzz[k] * ab_x + g_0_x_0_xyyyzzz[k];

                g_0_x_x_yyzzzz[k] = g_0_yyzzzz[k] - g_0_x_0_yyzzzz[k] * ab_x + g_0_x_0_xyyzzzz[k];

                g_0_x_x_yzzzzz[k] = g_0_yzzzzz[k] - g_0_x_0_yzzzzz[k] * ab_x + g_0_x_0_xyzzzzz[k];

                g_0_x_x_zzzzzz[k] = g_0_zzzzzz[k] - g_0_x_0_zzzzzz[k] * ab_x + g_0_x_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_0_x_y_xxxxxx = cbuffer.data(pi_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_y_xxxxxy = cbuffer.data(pi_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_y_xxxxxz = cbuffer.data(pi_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_y_xxxxyy = cbuffer.data(pi_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_y_xxxxyz = cbuffer.data(pi_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_y_xxxxzz = cbuffer.data(pi_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_y_xxxyyy = cbuffer.data(pi_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_y_xxxyyz = cbuffer.data(pi_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_y_xxxyzz = cbuffer.data(pi_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_y_xxxzzz = cbuffer.data(pi_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_y_xxyyyy = cbuffer.data(pi_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_y_xxyyyz = cbuffer.data(pi_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_y_xxyyzz = cbuffer.data(pi_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_y_xxyzzz = cbuffer.data(pi_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_y_xxzzzz = cbuffer.data(pi_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_y_xyyyyy = cbuffer.data(pi_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_y_xyyyyz = cbuffer.data(pi_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_y_xyyyzz = cbuffer.data(pi_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_y_xyyzzz = cbuffer.data(pi_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_y_xyzzzz = cbuffer.data(pi_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_y_xzzzzz = cbuffer.data(pi_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_y_yyyyyy = cbuffer.data(pi_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_y_yyyyyz = cbuffer.data(pi_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_y_yyyyzz = cbuffer.data(pi_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_y_yyyzzz = cbuffer.data(pi_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_y_yyzzzz = cbuffer.data(pi_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_y_yzzzzz = cbuffer.data(pi_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_y_zzzzzz = cbuffer.data(pi_geom_01_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxx, g_0_x_0_xxxxxxy, g_0_x_0_xxxxxy, g_0_x_0_xxxxxyy, g_0_x_0_xxxxxyz, g_0_x_0_xxxxxz, g_0_x_0_xxxxyy, g_0_x_0_xxxxyyy, g_0_x_0_xxxxyyz, g_0_x_0_xxxxyz, g_0_x_0_xxxxyzz, g_0_x_0_xxxxzz, g_0_x_0_xxxyyy, g_0_x_0_xxxyyyy, g_0_x_0_xxxyyyz, g_0_x_0_xxxyyz, g_0_x_0_xxxyyzz, g_0_x_0_xxxyzz, g_0_x_0_xxxyzzz, g_0_x_0_xxxzzz, g_0_x_0_xxyyyy, g_0_x_0_xxyyyyy, g_0_x_0_xxyyyyz, g_0_x_0_xxyyyz, g_0_x_0_xxyyyzz, g_0_x_0_xxyyzz, g_0_x_0_xxyyzzz, g_0_x_0_xxyzzz, g_0_x_0_xxyzzzz, g_0_x_0_xxzzzz, g_0_x_0_xyyyyy, g_0_x_0_xyyyyyy, g_0_x_0_xyyyyyz, g_0_x_0_xyyyyz, g_0_x_0_xyyyyzz, g_0_x_0_xyyyzz, g_0_x_0_xyyyzzz, g_0_x_0_xyyzzz, g_0_x_0_xyyzzzz, g_0_x_0_xyzzzz, g_0_x_0_xyzzzzz, g_0_x_0_xzzzzz, g_0_x_0_yyyyyy, g_0_x_0_yyyyyyy, g_0_x_0_yyyyyyz, g_0_x_0_yyyyyz, g_0_x_0_yyyyyzz, g_0_x_0_yyyyzz, g_0_x_0_yyyyzzz, g_0_x_0_yyyzzz, g_0_x_0_yyyzzzz, g_0_x_0_yyzzzz, g_0_x_0_yyzzzzz, g_0_x_0_yzzzzz, g_0_x_0_yzzzzzz, g_0_x_0_zzzzzz, g_0_x_y_xxxxxx, g_0_x_y_xxxxxy, g_0_x_y_xxxxxz, g_0_x_y_xxxxyy, g_0_x_y_xxxxyz, g_0_x_y_xxxxzz, g_0_x_y_xxxyyy, g_0_x_y_xxxyyz, g_0_x_y_xxxyzz, g_0_x_y_xxxzzz, g_0_x_y_xxyyyy, g_0_x_y_xxyyyz, g_0_x_y_xxyyzz, g_0_x_y_xxyzzz, g_0_x_y_xxzzzz, g_0_x_y_xyyyyy, g_0_x_y_xyyyyz, g_0_x_y_xyyyzz, g_0_x_y_xyyzzz, g_0_x_y_xyzzzz, g_0_x_y_xzzzzz, g_0_x_y_yyyyyy, g_0_x_y_yyyyyz, g_0_x_y_yyyyzz, g_0_x_y_yyyzzz, g_0_x_y_yyzzzz, g_0_x_y_yzzzzz, g_0_x_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_y_xxxxxx[k] = -g_0_x_0_xxxxxx[k] * ab_y + g_0_x_0_xxxxxxy[k];

                g_0_x_y_xxxxxy[k] = -g_0_x_0_xxxxxy[k] * ab_y + g_0_x_0_xxxxxyy[k];

                g_0_x_y_xxxxxz[k] = -g_0_x_0_xxxxxz[k] * ab_y + g_0_x_0_xxxxxyz[k];

                g_0_x_y_xxxxyy[k] = -g_0_x_0_xxxxyy[k] * ab_y + g_0_x_0_xxxxyyy[k];

                g_0_x_y_xxxxyz[k] = -g_0_x_0_xxxxyz[k] * ab_y + g_0_x_0_xxxxyyz[k];

                g_0_x_y_xxxxzz[k] = -g_0_x_0_xxxxzz[k] * ab_y + g_0_x_0_xxxxyzz[k];

                g_0_x_y_xxxyyy[k] = -g_0_x_0_xxxyyy[k] * ab_y + g_0_x_0_xxxyyyy[k];

                g_0_x_y_xxxyyz[k] = -g_0_x_0_xxxyyz[k] * ab_y + g_0_x_0_xxxyyyz[k];

                g_0_x_y_xxxyzz[k] = -g_0_x_0_xxxyzz[k] * ab_y + g_0_x_0_xxxyyzz[k];

                g_0_x_y_xxxzzz[k] = -g_0_x_0_xxxzzz[k] * ab_y + g_0_x_0_xxxyzzz[k];

                g_0_x_y_xxyyyy[k] = -g_0_x_0_xxyyyy[k] * ab_y + g_0_x_0_xxyyyyy[k];

                g_0_x_y_xxyyyz[k] = -g_0_x_0_xxyyyz[k] * ab_y + g_0_x_0_xxyyyyz[k];

                g_0_x_y_xxyyzz[k] = -g_0_x_0_xxyyzz[k] * ab_y + g_0_x_0_xxyyyzz[k];

                g_0_x_y_xxyzzz[k] = -g_0_x_0_xxyzzz[k] * ab_y + g_0_x_0_xxyyzzz[k];

                g_0_x_y_xxzzzz[k] = -g_0_x_0_xxzzzz[k] * ab_y + g_0_x_0_xxyzzzz[k];

                g_0_x_y_xyyyyy[k] = -g_0_x_0_xyyyyy[k] * ab_y + g_0_x_0_xyyyyyy[k];

                g_0_x_y_xyyyyz[k] = -g_0_x_0_xyyyyz[k] * ab_y + g_0_x_0_xyyyyyz[k];

                g_0_x_y_xyyyzz[k] = -g_0_x_0_xyyyzz[k] * ab_y + g_0_x_0_xyyyyzz[k];

                g_0_x_y_xyyzzz[k] = -g_0_x_0_xyyzzz[k] * ab_y + g_0_x_0_xyyyzzz[k];

                g_0_x_y_xyzzzz[k] = -g_0_x_0_xyzzzz[k] * ab_y + g_0_x_0_xyyzzzz[k];

                g_0_x_y_xzzzzz[k] = -g_0_x_0_xzzzzz[k] * ab_y + g_0_x_0_xyzzzzz[k];

                g_0_x_y_yyyyyy[k] = -g_0_x_0_yyyyyy[k] * ab_y + g_0_x_0_yyyyyyy[k];

                g_0_x_y_yyyyyz[k] = -g_0_x_0_yyyyyz[k] * ab_y + g_0_x_0_yyyyyyz[k];

                g_0_x_y_yyyyzz[k] = -g_0_x_0_yyyyzz[k] * ab_y + g_0_x_0_yyyyyzz[k];

                g_0_x_y_yyyzzz[k] = -g_0_x_0_yyyzzz[k] * ab_y + g_0_x_0_yyyyzzz[k];

                g_0_x_y_yyzzzz[k] = -g_0_x_0_yyzzzz[k] * ab_y + g_0_x_0_yyyzzzz[k];

                g_0_x_y_yzzzzz[k] = -g_0_x_0_yzzzzz[k] * ab_y + g_0_x_0_yyzzzzz[k];

                g_0_x_y_zzzzzz[k] = -g_0_x_0_zzzzzz[k] * ab_y + g_0_x_0_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_0_x_z_xxxxxx = cbuffer.data(pi_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_z_xxxxxy = cbuffer.data(pi_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_z_xxxxxz = cbuffer.data(pi_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_z_xxxxyy = cbuffer.data(pi_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_z_xxxxyz = cbuffer.data(pi_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_z_xxxxzz = cbuffer.data(pi_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_z_xxxyyy = cbuffer.data(pi_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_z_xxxyyz = cbuffer.data(pi_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_z_xxxyzz = cbuffer.data(pi_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_z_xxxzzz = cbuffer.data(pi_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_z_xxyyyy = cbuffer.data(pi_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_z_xxyyyz = cbuffer.data(pi_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_z_xxyyzz = cbuffer.data(pi_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_z_xxyzzz = cbuffer.data(pi_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_z_xxzzzz = cbuffer.data(pi_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_z_xyyyyy = cbuffer.data(pi_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_z_xyyyyz = cbuffer.data(pi_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_z_xyyyzz = cbuffer.data(pi_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_z_xyyzzz = cbuffer.data(pi_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_z_xyzzzz = cbuffer.data(pi_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_z_xzzzzz = cbuffer.data(pi_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_z_yyyyyy = cbuffer.data(pi_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_z_yyyyyz = cbuffer.data(pi_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_z_yyyyzz = cbuffer.data(pi_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_z_yyyzzz = cbuffer.data(pi_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_z_yyzzzz = cbuffer.data(pi_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_z_yzzzzz = cbuffer.data(pi_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_z_zzzzzz = cbuffer.data(pi_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxx, g_0_x_0_xxxxxxz, g_0_x_0_xxxxxy, g_0_x_0_xxxxxyz, g_0_x_0_xxxxxz, g_0_x_0_xxxxxzz, g_0_x_0_xxxxyy, g_0_x_0_xxxxyyz, g_0_x_0_xxxxyz, g_0_x_0_xxxxyzz, g_0_x_0_xxxxzz, g_0_x_0_xxxxzzz, g_0_x_0_xxxyyy, g_0_x_0_xxxyyyz, g_0_x_0_xxxyyz, g_0_x_0_xxxyyzz, g_0_x_0_xxxyzz, g_0_x_0_xxxyzzz, g_0_x_0_xxxzzz, g_0_x_0_xxxzzzz, g_0_x_0_xxyyyy, g_0_x_0_xxyyyyz, g_0_x_0_xxyyyz, g_0_x_0_xxyyyzz, g_0_x_0_xxyyzz, g_0_x_0_xxyyzzz, g_0_x_0_xxyzzz, g_0_x_0_xxyzzzz, g_0_x_0_xxzzzz, g_0_x_0_xxzzzzz, g_0_x_0_xyyyyy, g_0_x_0_xyyyyyz, g_0_x_0_xyyyyz, g_0_x_0_xyyyyzz, g_0_x_0_xyyyzz, g_0_x_0_xyyyzzz, g_0_x_0_xyyzzz, g_0_x_0_xyyzzzz, g_0_x_0_xyzzzz, g_0_x_0_xyzzzzz, g_0_x_0_xzzzzz, g_0_x_0_xzzzzzz, g_0_x_0_yyyyyy, g_0_x_0_yyyyyyz, g_0_x_0_yyyyyz, g_0_x_0_yyyyyzz, g_0_x_0_yyyyzz, g_0_x_0_yyyyzzz, g_0_x_0_yyyzzz, g_0_x_0_yyyzzzz, g_0_x_0_yyzzzz, g_0_x_0_yyzzzzz, g_0_x_0_yzzzzz, g_0_x_0_yzzzzzz, g_0_x_0_zzzzzz, g_0_x_0_zzzzzzz, g_0_x_z_xxxxxx, g_0_x_z_xxxxxy, g_0_x_z_xxxxxz, g_0_x_z_xxxxyy, g_0_x_z_xxxxyz, g_0_x_z_xxxxzz, g_0_x_z_xxxyyy, g_0_x_z_xxxyyz, g_0_x_z_xxxyzz, g_0_x_z_xxxzzz, g_0_x_z_xxyyyy, g_0_x_z_xxyyyz, g_0_x_z_xxyyzz, g_0_x_z_xxyzzz, g_0_x_z_xxzzzz, g_0_x_z_xyyyyy, g_0_x_z_xyyyyz, g_0_x_z_xyyyzz, g_0_x_z_xyyzzz, g_0_x_z_xyzzzz, g_0_x_z_xzzzzz, g_0_x_z_yyyyyy, g_0_x_z_yyyyyz, g_0_x_z_yyyyzz, g_0_x_z_yyyzzz, g_0_x_z_yyzzzz, g_0_x_z_yzzzzz, g_0_x_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_z_xxxxxx[k] = -g_0_x_0_xxxxxx[k] * ab_z + g_0_x_0_xxxxxxz[k];

                g_0_x_z_xxxxxy[k] = -g_0_x_0_xxxxxy[k] * ab_z + g_0_x_0_xxxxxyz[k];

                g_0_x_z_xxxxxz[k] = -g_0_x_0_xxxxxz[k] * ab_z + g_0_x_0_xxxxxzz[k];

                g_0_x_z_xxxxyy[k] = -g_0_x_0_xxxxyy[k] * ab_z + g_0_x_0_xxxxyyz[k];

                g_0_x_z_xxxxyz[k] = -g_0_x_0_xxxxyz[k] * ab_z + g_0_x_0_xxxxyzz[k];

                g_0_x_z_xxxxzz[k] = -g_0_x_0_xxxxzz[k] * ab_z + g_0_x_0_xxxxzzz[k];

                g_0_x_z_xxxyyy[k] = -g_0_x_0_xxxyyy[k] * ab_z + g_0_x_0_xxxyyyz[k];

                g_0_x_z_xxxyyz[k] = -g_0_x_0_xxxyyz[k] * ab_z + g_0_x_0_xxxyyzz[k];

                g_0_x_z_xxxyzz[k] = -g_0_x_0_xxxyzz[k] * ab_z + g_0_x_0_xxxyzzz[k];

                g_0_x_z_xxxzzz[k] = -g_0_x_0_xxxzzz[k] * ab_z + g_0_x_0_xxxzzzz[k];

                g_0_x_z_xxyyyy[k] = -g_0_x_0_xxyyyy[k] * ab_z + g_0_x_0_xxyyyyz[k];

                g_0_x_z_xxyyyz[k] = -g_0_x_0_xxyyyz[k] * ab_z + g_0_x_0_xxyyyzz[k];

                g_0_x_z_xxyyzz[k] = -g_0_x_0_xxyyzz[k] * ab_z + g_0_x_0_xxyyzzz[k];

                g_0_x_z_xxyzzz[k] = -g_0_x_0_xxyzzz[k] * ab_z + g_0_x_0_xxyzzzz[k];

                g_0_x_z_xxzzzz[k] = -g_0_x_0_xxzzzz[k] * ab_z + g_0_x_0_xxzzzzz[k];

                g_0_x_z_xyyyyy[k] = -g_0_x_0_xyyyyy[k] * ab_z + g_0_x_0_xyyyyyz[k];

                g_0_x_z_xyyyyz[k] = -g_0_x_0_xyyyyz[k] * ab_z + g_0_x_0_xyyyyzz[k];

                g_0_x_z_xyyyzz[k] = -g_0_x_0_xyyyzz[k] * ab_z + g_0_x_0_xyyyzzz[k];

                g_0_x_z_xyyzzz[k] = -g_0_x_0_xyyzzz[k] * ab_z + g_0_x_0_xyyzzzz[k];

                g_0_x_z_xyzzzz[k] = -g_0_x_0_xyzzzz[k] * ab_z + g_0_x_0_xyzzzzz[k];

                g_0_x_z_xzzzzz[k] = -g_0_x_0_xzzzzz[k] * ab_z + g_0_x_0_xzzzzzz[k];

                g_0_x_z_yyyyyy[k] = -g_0_x_0_yyyyyy[k] * ab_z + g_0_x_0_yyyyyyz[k];

                g_0_x_z_yyyyyz[k] = -g_0_x_0_yyyyyz[k] * ab_z + g_0_x_0_yyyyyzz[k];

                g_0_x_z_yyyyzz[k] = -g_0_x_0_yyyyzz[k] * ab_z + g_0_x_0_yyyyzzz[k];

                g_0_x_z_yyyzzz[k] = -g_0_x_0_yyyzzz[k] * ab_z + g_0_x_0_yyyzzzz[k];

                g_0_x_z_yyzzzz[k] = -g_0_x_0_yyzzzz[k] * ab_z + g_0_x_0_yyzzzzz[k];

                g_0_x_z_yzzzzz[k] = -g_0_x_0_yzzzzz[k] * ab_z + g_0_x_0_yzzzzzz[k];

                g_0_x_z_zzzzzz[k] = -g_0_x_0_zzzzzz[k] * ab_z + g_0_x_0_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_0_y_x_xxxxxx = cbuffer.data(pi_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_x_xxxxxy = cbuffer.data(pi_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_x_xxxxxz = cbuffer.data(pi_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_x_xxxxyy = cbuffer.data(pi_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_x_xxxxyz = cbuffer.data(pi_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_x_xxxxzz = cbuffer.data(pi_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_x_xxxyyy = cbuffer.data(pi_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_x_xxxyyz = cbuffer.data(pi_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_x_xxxyzz = cbuffer.data(pi_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_x_xxxzzz = cbuffer.data(pi_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_x_xxyyyy = cbuffer.data(pi_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_x_xxyyyz = cbuffer.data(pi_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_x_xxyyzz = cbuffer.data(pi_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_x_xxyzzz = cbuffer.data(pi_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_x_xxzzzz = cbuffer.data(pi_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_x_xyyyyy = cbuffer.data(pi_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_x_xyyyyz = cbuffer.data(pi_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_x_xyyyzz = cbuffer.data(pi_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_x_xyyzzz = cbuffer.data(pi_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_x_xyzzzz = cbuffer.data(pi_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_x_xzzzzz = cbuffer.data(pi_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_x_yyyyyy = cbuffer.data(pi_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_x_yyyyyz = cbuffer.data(pi_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_x_yyyyzz = cbuffer.data(pi_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_x_yyyzzz = cbuffer.data(pi_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_x_yyzzzz = cbuffer.data(pi_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_x_yzzzzz = cbuffer.data(pi_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_x_zzzzzz = cbuffer.data(pi_geom_01_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxxx, g_0_y_0_xxxxxxx, g_0_y_0_xxxxxxy, g_0_y_0_xxxxxxz, g_0_y_0_xxxxxy, g_0_y_0_xxxxxyy, g_0_y_0_xxxxxyz, g_0_y_0_xxxxxz, g_0_y_0_xxxxxzz, g_0_y_0_xxxxyy, g_0_y_0_xxxxyyy, g_0_y_0_xxxxyyz, g_0_y_0_xxxxyz, g_0_y_0_xxxxyzz, g_0_y_0_xxxxzz, g_0_y_0_xxxxzzz, g_0_y_0_xxxyyy, g_0_y_0_xxxyyyy, g_0_y_0_xxxyyyz, g_0_y_0_xxxyyz, g_0_y_0_xxxyyzz, g_0_y_0_xxxyzz, g_0_y_0_xxxyzzz, g_0_y_0_xxxzzz, g_0_y_0_xxxzzzz, g_0_y_0_xxyyyy, g_0_y_0_xxyyyyy, g_0_y_0_xxyyyyz, g_0_y_0_xxyyyz, g_0_y_0_xxyyyzz, g_0_y_0_xxyyzz, g_0_y_0_xxyyzzz, g_0_y_0_xxyzzz, g_0_y_0_xxyzzzz, g_0_y_0_xxzzzz, g_0_y_0_xxzzzzz, g_0_y_0_xyyyyy, g_0_y_0_xyyyyyy, g_0_y_0_xyyyyyz, g_0_y_0_xyyyyz, g_0_y_0_xyyyyzz, g_0_y_0_xyyyzz, g_0_y_0_xyyyzzz, g_0_y_0_xyyzzz, g_0_y_0_xyyzzzz, g_0_y_0_xyzzzz, g_0_y_0_xyzzzzz, g_0_y_0_xzzzzz, g_0_y_0_xzzzzzz, g_0_y_0_yyyyyy, g_0_y_0_yyyyyz, g_0_y_0_yyyyzz, g_0_y_0_yyyzzz, g_0_y_0_yyzzzz, g_0_y_0_yzzzzz, g_0_y_0_zzzzzz, g_0_y_x_xxxxxx, g_0_y_x_xxxxxy, g_0_y_x_xxxxxz, g_0_y_x_xxxxyy, g_0_y_x_xxxxyz, g_0_y_x_xxxxzz, g_0_y_x_xxxyyy, g_0_y_x_xxxyyz, g_0_y_x_xxxyzz, g_0_y_x_xxxzzz, g_0_y_x_xxyyyy, g_0_y_x_xxyyyz, g_0_y_x_xxyyzz, g_0_y_x_xxyzzz, g_0_y_x_xxzzzz, g_0_y_x_xyyyyy, g_0_y_x_xyyyyz, g_0_y_x_xyyyzz, g_0_y_x_xyyzzz, g_0_y_x_xyzzzz, g_0_y_x_xzzzzz, g_0_y_x_yyyyyy, g_0_y_x_yyyyyz, g_0_y_x_yyyyzz, g_0_y_x_yyyzzz, g_0_y_x_yyzzzz, g_0_y_x_yzzzzz, g_0_y_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_x_xxxxxx[k] = -g_0_y_0_xxxxxx[k] * ab_x + g_0_y_0_xxxxxxx[k];

                g_0_y_x_xxxxxy[k] = -g_0_y_0_xxxxxy[k] * ab_x + g_0_y_0_xxxxxxy[k];

                g_0_y_x_xxxxxz[k] = -g_0_y_0_xxxxxz[k] * ab_x + g_0_y_0_xxxxxxz[k];

                g_0_y_x_xxxxyy[k] = -g_0_y_0_xxxxyy[k] * ab_x + g_0_y_0_xxxxxyy[k];

                g_0_y_x_xxxxyz[k] = -g_0_y_0_xxxxyz[k] * ab_x + g_0_y_0_xxxxxyz[k];

                g_0_y_x_xxxxzz[k] = -g_0_y_0_xxxxzz[k] * ab_x + g_0_y_0_xxxxxzz[k];

                g_0_y_x_xxxyyy[k] = -g_0_y_0_xxxyyy[k] * ab_x + g_0_y_0_xxxxyyy[k];

                g_0_y_x_xxxyyz[k] = -g_0_y_0_xxxyyz[k] * ab_x + g_0_y_0_xxxxyyz[k];

                g_0_y_x_xxxyzz[k] = -g_0_y_0_xxxyzz[k] * ab_x + g_0_y_0_xxxxyzz[k];

                g_0_y_x_xxxzzz[k] = -g_0_y_0_xxxzzz[k] * ab_x + g_0_y_0_xxxxzzz[k];

                g_0_y_x_xxyyyy[k] = -g_0_y_0_xxyyyy[k] * ab_x + g_0_y_0_xxxyyyy[k];

                g_0_y_x_xxyyyz[k] = -g_0_y_0_xxyyyz[k] * ab_x + g_0_y_0_xxxyyyz[k];

                g_0_y_x_xxyyzz[k] = -g_0_y_0_xxyyzz[k] * ab_x + g_0_y_0_xxxyyzz[k];

                g_0_y_x_xxyzzz[k] = -g_0_y_0_xxyzzz[k] * ab_x + g_0_y_0_xxxyzzz[k];

                g_0_y_x_xxzzzz[k] = -g_0_y_0_xxzzzz[k] * ab_x + g_0_y_0_xxxzzzz[k];

                g_0_y_x_xyyyyy[k] = -g_0_y_0_xyyyyy[k] * ab_x + g_0_y_0_xxyyyyy[k];

                g_0_y_x_xyyyyz[k] = -g_0_y_0_xyyyyz[k] * ab_x + g_0_y_0_xxyyyyz[k];

                g_0_y_x_xyyyzz[k] = -g_0_y_0_xyyyzz[k] * ab_x + g_0_y_0_xxyyyzz[k];

                g_0_y_x_xyyzzz[k] = -g_0_y_0_xyyzzz[k] * ab_x + g_0_y_0_xxyyzzz[k];

                g_0_y_x_xyzzzz[k] = -g_0_y_0_xyzzzz[k] * ab_x + g_0_y_0_xxyzzzz[k];

                g_0_y_x_xzzzzz[k] = -g_0_y_0_xzzzzz[k] * ab_x + g_0_y_0_xxzzzzz[k];

                g_0_y_x_yyyyyy[k] = -g_0_y_0_yyyyyy[k] * ab_x + g_0_y_0_xyyyyyy[k];

                g_0_y_x_yyyyyz[k] = -g_0_y_0_yyyyyz[k] * ab_x + g_0_y_0_xyyyyyz[k];

                g_0_y_x_yyyyzz[k] = -g_0_y_0_yyyyzz[k] * ab_x + g_0_y_0_xyyyyzz[k];

                g_0_y_x_yyyzzz[k] = -g_0_y_0_yyyzzz[k] * ab_x + g_0_y_0_xyyyzzz[k];

                g_0_y_x_yyzzzz[k] = -g_0_y_0_yyzzzz[k] * ab_x + g_0_y_0_xyyzzzz[k];

                g_0_y_x_yzzzzz[k] = -g_0_y_0_yzzzzz[k] * ab_x + g_0_y_0_xyzzzzz[k];

                g_0_y_x_zzzzzz[k] = -g_0_y_0_zzzzzz[k] * ab_x + g_0_y_0_xzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_0_y_y_xxxxxx = cbuffer.data(pi_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_y_xxxxxy = cbuffer.data(pi_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_y_xxxxxz = cbuffer.data(pi_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_y_xxxxyy = cbuffer.data(pi_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_y_xxxxyz = cbuffer.data(pi_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_y_xxxxzz = cbuffer.data(pi_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_y_xxxyyy = cbuffer.data(pi_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_y_xxxyyz = cbuffer.data(pi_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_y_xxxyzz = cbuffer.data(pi_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_y_xxxzzz = cbuffer.data(pi_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_y_xxyyyy = cbuffer.data(pi_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_y_xxyyyz = cbuffer.data(pi_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_y_xxyyzz = cbuffer.data(pi_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_y_xxyzzz = cbuffer.data(pi_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_y_xxzzzz = cbuffer.data(pi_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_y_xyyyyy = cbuffer.data(pi_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_y_xyyyyz = cbuffer.data(pi_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_y_xyyyzz = cbuffer.data(pi_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_y_xyyzzz = cbuffer.data(pi_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_y_xyzzzz = cbuffer.data(pi_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_y_xzzzzz = cbuffer.data(pi_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_y_yyyyyy = cbuffer.data(pi_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_y_yyyyyz = cbuffer.data(pi_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_y_yyyyzz = cbuffer.data(pi_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_y_yyyzzz = cbuffer.data(pi_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_y_yyzzzz = cbuffer.data(pi_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_y_yzzzzz = cbuffer.data(pi_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_y_zzzzzz = cbuffer.data(pi_geom_01_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyzzz, g_0_xyzzzz, g_0_xzzzzz, g_0_y_0_xxxxxx, g_0_y_0_xxxxxxy, g_0_y_0_xxxxxy, g_0_y_0_xxxxxyy, g_0_y_0_xxxxxyz, g_0_y_0_xxxxxz, g_0_y_0_xxxxyy, g_0_y_0_xxxxyyy, g_0_y_0_xxxxyyz, g_0_y_0_xxxxyz, g_0_y_0_xxxxyzz, g_0_y_0_xxxxzz, g_0_y_0_xxxyyy, g_0_y_0_xxxyyyy, g_0_y_0_xxxyyyz, g_0_y_0_xxxyyz, g_0_y_0_xxxyyzz, g_0_y_0_xxxyzz, g_0_y_0_xxxyzzz, g_0_y_0_xxxzzz, g_0_y_0_xxyyyy, g_0_y_0_xxyyyyy, g_0_y_0_xxyyyyz, g_0_y_0_xxyyyz, g_0_y_0_xxyyyzz, g_0_y_0_xxyyzz, g_0_y_0_xxyyzzz, g_0_y_0_xxyzzz, g_0_y_0_xxyzzzz, g_0_y_0_xxzzzz, g_0_y_0_xyyyyy, g_0_y_0_xyyyyyy, g_0_y_0_xyyyyyz, g_0_y_0_xyyyyz, g_0_y_0_xyyyyzz, g_0_y_0_xyyyzz, g_0_y_0_xyyyzzz, g_0_y_0_xyyzzz, g_0_y_0_xyyzzzz, g_0_y_0_xyzzzz, g_0_y_0_xyzzzzz, g_0_y_0_xzzzzz, g_0_y_0_yyyyyy, g_0_y_0_yyyyyyy, g_0_y_0_yyyyyyz, g_0_y_0_yyyyyz, g_0_y_0_yyyyyzz, g_0_y_0_yyyyzz, g_0_y_0_yyyyzzz, g_0_y_0_yyyzzz, g_0_y_0_yyyzzzz, g_0_y_0_yyzzzz, g_0_y_0_yyzzzzz, g_0_y_0_yzzzzz, g_0_y_0_yzzzzzz, g_0_y_0_zzzzzz, g_0_y_y_xxxxxx, g_0_y_y_xxxxxy, g_0_y_y_xxxxxz, g_0_y_y_xxxxyy, g_0_y_y_xxxxyz, g_0_y_y_xxxxzz, g_0_y_y_xxxyyy, g_0_y_y_xxxyyz, g_0_y_y_xxxyzz, g_0_y_y_xxxzzz, g_0_y_y_xxyyyy, g_0_y_y_xxyyyz, g_0_y_y_xxyyzz, g_0_y_y_xxyzzz, g_0_y_y_xxzzzz, g_0_y_y_xyyyyy, g_0_y_y_xyyyyz, g_0_y_y_xyyyzz, g_0_y_y_xyyzzz, g_0_y_y_xyzzzz, g_0_y_y_xzzzzz, g_0_y_y_yyyyyy, g_0_y_y_yyyyyz, g_0_y_y_yyyyzz, g_0_y_y_yyyzzz, g_0_y_y_yyzzzz, g_0_y_y_yzzzzz, g_0_y_y_zzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_y_xxxxxx[k] = g_0_xxxxxx[k] - g_0_y_0_xxxxxx[k] * ab_y + g_0_y_0_xxxxxxy[k];

                g_0_y_y_xxxxxy[k] = g_0_xxxxxy[k] - g_0_y_0_xxxxxy[k] * ab_y + g_0_y_0_xxxxxyy[k];

                g_0_y_y_xxxxxz[k] = g_0_xxxxxz[k] - g_0_y_0_xxxxxz[k] * ab_y + g_0_y_0_xxxxxyz[k];

                g_0_y_y_xxxxyy[k] = g_0_xxxxyy[k] - g_0_y_0_xxxxyy[k] * ab_y + g_0_y_0_xxxxyyy[k];

                g_0_y_y_xxxxyz[k] = g_0_xxxxyz[k] - g_0_y_0_xxxxyz[k] * ab_y + g_0_y_0_xxxxyyz[k];

                g_0_y_y_xxxxzz[k] = g_0_xxxxzz[k] - g_0_y_0_xxxxzz[k] * ab_y + g_0_y_0_xxxxyzz[k];

                g_0_y_y_xxxyyy[k] = g_0_xxxyyy[k] - g_0_y_0_xxxyyy[k] * ab_y + g_0_y_0_xxxyyyy[k];

                g_0_y_y_xxxyyz[k] = g_0_xxxyyz[k] - g_0_y_0_xxxyyz[k] * ab_y + g_0_y_0_xxxyyyz[k];

                g_0_y_y_xxxyzz[k] = g_0_xxxyzz[k] - g_0_y_0_xxxyzz[k] * ab_y + g_0_y_0_xxxyyzz[k];

                g_0_y_y_xxxzzz[k] = g_0_xxxzzz[k] - g_0_y_0_xxxzzz[k] * ab_y + g_0_y_0_xxxyzzz[k];

                g_0_y_y_xxyyyy[k] = g_0_xxyyyy[k] - g_0_y_0_xxyyyy[k] * ab_y + g_0_y_0_xxyyyyy[k];

                g_0_y_y_xxyyyz[k] = g_0_xxyyyz[k] - g_0_y_0_xxyyyz[k] * ab_y + g_0_y_0_xxyyyyz[k];

                g_0_y_y_xxyyzz[k] = g_0_xxyyzz[k] - g_0_y_0_xxyyzz[k] * ab_y + g_0_y_0_xxyyyzz[k];

                g_0_y_y_xxyzzz[k] = g_0_xxyzzz[k] - g_0_y_0_xxyzzz[k] * ab_y + g_0_y_0_xxyyzzz[k];

                g_0_y_y_xxzzzz[k] = g_0_xxzzzz[k] - g_0_y_0_xxzzzz[k] * ab_y + g_0_y_0_xxyzzzz[k];

                g_0_y_y_xyyyyy[k] = g_0_xyyyyy[k] - g_0_y_0_xyyyyy[k] * ab_y + g_0_y_0_xyyyyyy[k];

                g_0_y_y_xyyyyz[k] = g_0_xyyyyz[k] - g_0_y_0_xyyyyz[k] * ab_y + g_0_y_0_xyyyyyz[k];

                g_0_y_y_xyyyzz[k] = g_0_xyyyzz[k] - g_0_y_0_xyyyzz[k] * ab_y + g_0_y_0_xyyyyzz[k];

                g_0_y_y_xyyzzz[k] = g_0_xyyzzz[k] - g_0_y_0_xyyzzz[k] * ab_y + g_0_y_0_xyyyzzz[k];

                g_0_y_y_xyzzzz[k] = g_0_xyzzzz[k] - g_0_y_0_xyzzzz[k] * ab_y + g_0_y_0_xyyzzzz[k];

                g_0_y_y_xzzzzz[k] = g_0_xzzzzz[k] - g_0_y_0_xzzzzz[k] * ab_y + g_0_y_0_xyzzzzz[k];

                g_0_y_y_yyyyyy[k] = g_0_yyyyyy[k] - g_0_y_0_yyyyyy[k] * ab_y + g_0_y_0_yyyyyyy[k];

                g_0_y_y_yyyyyz[k] = g_0_yyyyyz[k] - g_0_y_0_yyyyyz[k] * ab_y + g_0_y_0_yyyyyyz[k];

                g_0_y_y_yyyyzz[k] = g_0_yyyyzz[k] - g_0_y_0_yyyyzz[k] * ab_y + g_0_y_0_yyyyyzz[k];

                g_0_y_y_yyyzzz[k] = g_0_yyyzzz[k] - g_0_y_0_yyyzzz[k] * ab_y + g_0_y_0_yyyyzzz[k];

                g_0_y_y_yyzzzz[k] = g_0_yyzzzz[k] - g_0_y_0_yyzzzz[k] * ab_y + g_0_y_0_yyyzzzz[k];

                g_0_y_y_yzzzzz[k] = g_0_yzzzzz[k] - g_0_y_0_yzzzzz[k] * ab_y + g_0_y_0_yyzzzzz[k];

                g_0_y_y_zzzzzz[k] = g_0_zzzzzz[k] - g_0_y_0_zzzzzz[k] * ab_y + g_0_y_0_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_0_y_z_xxxxxx = cbuffer.data(pi_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_z_xxxxxy = cbuffer.data(pi_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_z_xxxxxz = cbuffer.data(pi_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_z_xxxxyy = cbuffer.data(pi_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_z_xxxxyz = cbuffer.data(pi_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_z_xxxxzz = cbuffer.data(pi_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_z_xxxyyy = cbuffer.data(pi_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_z_xxxyyz = cbuffer.data(pi_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_z_xxxyzz = cbuffer.data(pi_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_z_xxxzzz = cbuffer.data(pi_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_z_xxyyyy = cbuffer.data(pi_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_z_xxyyyz = cbuffer.data(pi_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_z_xxyyzz = cbuffer.data(pi_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_z_xxyzzz = cbuffer.data(pi_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_z_xxzzzz = cbuffer.data(pi_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_z_xyyyyy = cbuffer.data(pi_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_z_xyyyyz = cbuffer.data(pi_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_z_xyyyzz = cbuffer.data(pi_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_z_xyyzzz = cbuffer.data(pi_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_z_xyzzzz = cbuffer.data(pi_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_z_xzzzzz = cbuffer.data(pi_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_z_yyyyyy = cbuffer.data(pi_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_z_yyyyyz = cbuffer.data(pi_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_z_yyyyzz = cbuffer.data(pi_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_z_yyyzzz = cbuffer.data(pi_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_z_yyzzzz = cbuffer.data(pi_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_z_yzzzzz = cbuffer.data(pi_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_z_zzzzzz = cbuffer.data(pi_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxxx, g_0_y_0_xxxxxxz, g_0_y_0_xxxxxy, g_0_y_0_xxxxxyz, g_0_y_0_xxxxxz, g_0_y_0_xxxxxzz, g_0_y_0_xxxxyy, g_0_y_0_xxxxyyz, g_0_y_0_xxxxyz, g_0_y_0_xxxxyzz, g_0_y_0_xxxxzz, g_0_y_0_xxxxzzz, g_0_y_0_xxxyyy, g_0_y_0_xxxyyyz, g_0_y_0_xxxyyz, g_0_y_0_xxxyyzz, g_0_y_0_xxxyzz, g_0_y_0_xxxyzzz, g_0_y_0_xxxzzz, g_0_y_0_xxxzzzz, g_0_y_0_xxyyyy, g_0_y_0_xxyyyyz, g_0_y_0_xxyyyz, g_0_y_0_xxyyyzz, g_0_y_0_xxyyzz, g_0_y_0_xxyyzzz, g_0_y_0_xxyzzz, g_0_y_0_xxyzzzz, g_0_y_0_xxzzzz, g_0_y_0_xxzzzzz, g_0_y_0_xyyyyy, g_0_y_0_xyyyyyz, g_0_y_0_xyyyyz, g_0_y_0_xyyyyzz, g_0_y_0_xyyyzz, g_0_y_0_xyyyzzz, g_0_y_0_xyyzzz, g_0_y_0_xyyzzzz, g_0_y_0_xyzzzz, g_0_y_0_xyzzzzz, g_0_y_0_xzzzzz, g_0_y_0_xzzzzzz, g_0_y_0_yyyyyy, g_0_y_0_yyyyyyz, g_0_y_0_yyyyyz, g_0_y_0_yyyyyzz, g_0_y_0_yyyyzz, g_0_y_0_yyyyzzz, g_0_y_0_yyyzzz, g_0_y_0_yyyzzzz, g_0_y_0_yyzzzz, g_0_y_0_yyzzzzz, g_0_y_0_yzzzzz, g_0_y_0_yzzzzzz, g_0_y_0_zzzzzz, g_0_y_0_zzzzzzz, g_0_y_z_xxxxxx, g_0_y_z_xxxxxy, g_0_y_z_xxxxxz, g_0_y_z_xxxxyy, g_0_y_z_xxxxyz, g_0_y_z_xxxxzz, g_0_y_z_xxxyyy, g_0_y_z_xxxyyz, g_0_y_z_xxxyzz, g_0_y_z_xxxzzz, g_0_y_z_xxyyyy, g_0_y_z_xxyyyz, g_0_y_z_xxyyzz, g_0_y_z_xxyzzz, g_0_y_z_xxzzzz, g_0_y_z_xyyyyy, g_0_y_z_xyyyyz, g_0_y_z_xyyyzz, g_0_y_z_xyyzzz, g_0_y_z_xyzzzz, g_0_y_z_xzzzzz, g_0_y_z_yyyyyy, g_0_y_z_yyyyyz, g_0_y_z_yyyyzz, g_0_y_z_yyyzzz, g_0_y_z_yyzzzz, g_0_y_z_yzzzzz, g_0_y_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_z_xxxxxx[k] = -g_0_y_0_xxxxxx[k] * ab_z + g_0_y_0_xxxxxxz[k];

                g_0_y_z_xxxxxy[k] = -g_0_y_0_xxxxxy[k] * ab_z + g_0_y_0_xxxxxyz[k];

                g_0_y_z_xxxxxz[k] = -g_0_y_0_xxxxxz[k] * ab_z + g_0_y_0_xxxxxzz[k];

                g_0_y_z_xxxxyy[k] = -g_0_y_0_xxxxyy[k] * ab_z + g_0_y_0_xxxxyyz[k];

                g_0_y_z_xxxxyz[k] = -g_0_y_0_xxxxyz[k] * ab_z + g_0_y_0_xxxxyzz[k];

                g_0_y_z_xxxxzz[k] = -g_0_y_0_xxxxzz[k] * ab_z + g_0_y_0_xxxxzzz[k];

                g_0_y_z_xxxyyy[k] = -g_0_y_0_xxxyyy[k] * ab_z + g_0_y_0_xxxyyyz[k];

                g_0_y_z_xxxyyz[k] = -g_0_y_0_xxxyyz[k] * ab_z + g_0_y_0_xxxyyzz[k];

                g_0_y_z_xxxyzz[k] = -g_0_y_0_xxxyzz[k] * ab_z + g_0_y_0_xxxyzzz[k];

                g_0_y_z_xxxzzz[k] = -g_0_y_0_xxxzzz[k] * ab_z + g_0_y_0_xxxzzzz[k];

                g_0_y_z_xxyyyy[k] = -g_0_y_0_xxyyyy[k] * ab_z + g_0_y_0_xxyyyyz[k];

                g_0_y_z_xxyyyz[k] = -g_0_y_0_xxyyyz[k] * ab_z + g_0_y_0_xxyyyzz[k];

                g_0_y_z_xxyyzz[k] = -g_0_y_0_xxyyzz[k] * ab_z + g_0_y_0_xxyyzzz[k];

                g_0_y_z_xxyzzz[k] = -g_0_y_0_xxyzzz[k] * ab_z + g_0_y_0_xxyzzzz[k];

                g_0_y_z_xxzzzz[k] = -g_0_y_0_xxzzzz[k] * ab_z + g_0_y_0_xxzzzzz[k];

                g_0_y_z_xyyyyy[k] = -g_0_y_0_xyyyyy[k] * ab_z + g_0_y_0_xyyyyyz[k];

                g_0_y_z_xyyyyz[k] = -g_0_y_0_xyyyyz[k] * ab_z + g_0_y_0_xyyyyzz[k];

                g_0_y_z_xyyyzz[k] = -g_0_y_0_xyyyzz[k] * ab_z + g_0_y_0_xyyyzzz[k];

                g_0_y_z_xyyzzz[k] = -g_0_y_0_xyyzzz[k] * ab_z + g_0_y_0_xyyzzzz[k];

                g_0_y_z_xyzzzz[k] = -g_0_y_0_xyzzzz[k] * ab_z + g_0_y_0_xyzzzzz[k];

                g_0_y_z_xzzzzz[k] = -g_0_y_0_xzzzzz[k] * ab_z + g_0_y_0_xzzzzzz[k];

                g_0_y_z_yyyyyy[k] = -g_0_y_0_yyyyyy[k] * ab_z + g_0_y_0_yyyyyyz[k];

                g_0_y_z_yyyyyz[k] = -g_0_y_0_yyyyyz[k] * ab_z + g_0_y_0_yyyyyzz[k];

                g_0_y_z_yyyyzz[k] = -g_0_y_0_yyyyzz[k] * ab_z + g_0_y_0_yyyyzzz[k];

                g_0_y_z_yyyzzz[k] = -g_0_y_0_yyyzzz[k] * ab_z + g_0_y_0_yyyzzzz[k];

                g_0_y_z_yyzzzz[k] = -g_0_y_0_yyzzzz[k] * ab_z + g_0_y_0_yyzzzzz[k];

                g_0_y_z_yzzzzz[k] = -g_0_y_0_yzzzzz[k] * ab_z + g_0_y_0_yzzzzzz[k];

                g_0_y_z_zzzzzz[k] = -g_0_y_0_zzzzzz[k] * ab_z + g_0_y_0_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_0_z_x_xxxxxx = cbuffer.data(pi_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_z_x_xxxxxy = cbuffer.data(pi_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_z_x_xxxxxz = cbuffer.data(pi_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_z_x_xxxxyy = cbuffer.data(pi_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_z_x_xxxxyz = cbuffer.data(pi_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_z_x_xxxxzz = cbuffer.data(pi_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_z_x_xxxyyy = cbuffer.data(pi_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_z_x_xxxyyz = cbuffer.data(pi_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_z_x_xxxyzz = cbuffer.data(pi_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_z_x_xxxzzz = cbuffer.data(pi_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_z_x_xxyyyy = cbuffer.data(pi_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_z_x_xxyyyz = cbuffer.data(pi_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_x_xxyyzz = cbuffer.data(pi_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_x_xxyzzz = cbuffer.data(pi_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_x_xxzzzz = cbuffer.data(pi_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_x_xyyyyy = cbuffer.data(pi_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_x_xyyyyz = cbuffer.data(pi_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_x_xyyyzz = cbuffer.data(pi_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_x_xyyzzz = cbuffer.data(pi_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_x_xyzzzz = cbuffer.data(pi_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_x_xzzzzz = cbuffer.data(pi_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_z_x_yyyyyy = cbuffer.data(pi_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_z_x_yyyyyz = cbuffer.data(pi_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_z_x_yyyyzz = cbuffer.data(pi_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_z_x_yyyzzz = cbuffer.data(pi_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_z_x_yyzzzz = cbuffer.data(pi_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_z_x_yzzzzz = cbuffer.data(pi_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_z_x_zzzzzz = cbuffer.data(pi_geom_01_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxxx, g_0_z_0_xxxxxxx, g_0_z_0_xxxxxxy, g_0_z_0_xxxxxxz, g_0_z_0_xxxxxy, g_0_z_0_xxxxxyy, g_0_z_0_xxxxxyz, g_0_z_0_xxxxxz, g_0_z_0_xxxxxzz, g_0_z_0_xxxxyy, g_0_z_0_xxxxyyy, g_0_z_0_xxxxyyz, g_0_z_0_xxxxyz, g_0_z_0_xxxxyzz, g_0_z_0_xxxxzz, g_0_z_0_xxxxzzz, g_0_z_0_xxxyyy, g_0_z_0_xxxyyyy, g_0_z_0_xxxyyyz, g_0_z_0_xxxyyz, g_0_z_0_xxxyyzz, g_0_z_0_xxxyzz, g_0_z_0_xxxyzzz, g_0_z_0_xxxzzz, g_0_z_0_xxxzzzz, g_0_z_0_xxyyyy, g_0_z_0_xxyyyyy, g_0_z_0_xxyyyyz, g_0_z_0_xxyyyz, g_0_z_0_xxyyyzz, g_0_z_0_xxyyzz, g_0_z_0_xxyyzzz, g_0_z_0_xxyzzz, g_0_z_0_xxyzzzz, g_0_z_0_xxzzzz, g_0_z_0_xxzzzzz, g_0_z_0_xyyyyy, g_0_z_0_xyyyyyy, g_0_z_0_xyyyyyz, g_0_z_0_xyyyyz, g_0_z_0_xyyyyzz, g_0_z_0_xyyyzz, g_0_z_0_xyyyzzz, g_0_z_0_xyyzzz, g_0_z_0_xyyzzzz, g_0_z_0_xyzzzz, g_0_z_0_xyzzzzz, g_0_z_0_xzzzzz, g_0_z_0_xzzzzzz, g_0_z_0_yyyyyy, g_0_z_0_yyyyyz, g_0_z_0_yyyyzz, g_0_z_0_yyyzzz, g_0_z_0_yyzzzz, g_0_z_0_yzzzzz, g_0_z_0_zzzzzz, g_0_z_x_xxxxxx, g_0_z_x_xxxxxy, g_0_z_x_xxxxxz, g_0_z_x_xxxxyy, g_0_z_x_xxxxyz, g_0_z_x_xxxxzz, g_0_z_x_xxxyyy, g_0_z_x_xxxyyz, g_0_z_x_xxxyzz, g_0_z_x_xxxzzz, g_0_z_x_xxyyyy, g_0_z_x_xxyyyz, g_0_z_x_xxyyzz, g_0_z_x_xxyzzz, g_0_z_x_xxzzzz, g_0_z_x_xyyyyy, g_0_z_x_xyyyyz, g_0_z_x_xyyyzz, g_0_z_x_xyyzzz, g_0_z_x_xyzzzz, g_0_z_x_xzzzzz, g_0_z_x_yyyyyy, g_0_z_x_yyyyyz, g_0_z_x_yyyyzz, g_0_z_x_yyyzzz, g_0_z_x_yyzzzz, g_0_z_x_yzzzzz, g_0_z_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_x_xxxxxx[k] = -g_0_z_0_xxxxxx[k] * ab_x + g_0_z_0_xxxxxxx[k];

                g_0_z_x_xxxxxy[k] = -g_0_z_0_xxxxxy[k] * ab_x + g_0_z_0_xxxxxxy[k];

                g_0_z_x_xxxxxz[k] = -g_0_z_0_xxxxxz[k] * ab_x + g_0_z_0_xxxxxxz[k];

                g_0_z_x_xxxxyy[k] = -g_0_z_0_xxxxyy[k] * ab_x + g_0_z_0_xxxxxyy[k];

                g_0_z_x_xxxxyz[k] = -g_0_z_0_xxxxyz[k] * ab_x + g_0_z_0_xxxxxyz[k];

                g_0_z_x_xxxxzz[k] = -g_0_z_0_xxxxzz[k] * ab_x + g_0_z_0_xxxxxzz[k];

                g_0_z_x_xxxyyy[k] = -g_0_z_0_xxxyyy[k] * ab_x + g_0_z_0_xxxxyyy[k];

                g_0_z_x_xxxyyz[k] = -g_0_z_0_xxxyyz[k] * ab_x + g_0_z_0_xxxxyyz[k];

                g_0_z_x_xxxyzz[k] = -g_0_z_0_xxxyzz[k] * ab_x + g_0_z_0_xxxxyzz[k];

                g_0_z_x_xxxzzz[k] = -g_0_z_0_xxxzzz[k] * ab_x + g_0_z_0_xxxxzzz[k];

                g_0_z_x_xxyyyy[k] = -g_0_z_0_xxyyyy[k] * ab_x + g_0_z_0_xxxyyyy[k];

                g_0_z_x_xxyyyz[k] = -g_0_z_0_xxyyyz[k] * ab_x + g_0_z_0_xxxyyyz[k];

                g_0_z_x_xxyyzz[k] = -g_0_z_0_xxyyzz[k] * ab_x + g_0_z_0_xxxyyzz[k];

                g_0_z_x_xxyzzz[k] = -g_0_z_0_xxyzzz[k] * ab_x + g_0_z_0_xxxyzzz[k];

                g_0_z_x_xxzzzz[k] = -g_0_z_0_xxzzzz[k] * ab_x + g_0_z_0_xxxzzzz[k];

                g_0_z_x_xyyyyy[k] = -g_0_z_0_xyyyyy[k] * ab_x + g_0_z_0_xxyyyyy[k];

                g_0_z_x_xyyyyz[k] = -g_0_z_0_xyyyyz[k] * ab_x + g_0_z_0_xxyyyyz[k];

                g_0_z_x_xyyyzz[k] = -g_0_z_0_xyyyzz[k] * ab_x + g_0_z_0_xxyyyzz[k];

                g_0_z_x_xyyzzz[k] = -g_0_z_0_xyyzzz[k] * ab_x + g_0_z_0_xxyyzzz[k];

                g_0_z_x_xyzzzz[k] = -g_0_z_0_xyzzzz[k] * ab_x + g_0_z_0_xxyzzzz[k];

                g_0_z_x_xzzzzz[k] = -g_0_z_0_xzzzzz[k] * ab_x + g_0_z_0_xxzzzzz[k];

                g_0_z_x_yyyyyy[k] = -g_0_z_0_yyyyyy[k] * ab_x + g_0_z_0_xyyyyyy[k];

                g_0_z_x_yyyyyz[k] = -g_0_z_0_yyyyyz[k] * ab_x + g_0_z_0_xyyyyyz[k];

                g_0_z_x_yyyyzz[k] = -g_0_z_0_yyyyzz[k] * ab_x + g_0_z_0_xyyyyzz[k];

                g_0_z_x_yyyzzz[k] = -g_0_z_0_yyyzzz[k] * ab_x + g_0_z_0_xyyyzzz[k];

                g_0_z_x_yyzzzz[k] = -g_0_z_0_yyzzzz[k] * ab_x + g_0_z_0_xyyzzzz[k];

                g_0_z_x_yzzzzz[k] = -g_0_z_0_yzzzzz[k] * ab_x + g_0_z_0_xyzzzzz[k];

                g_0_z_x_zzzzzz[k] = -g_0_z_0_zzzzzz[k] * ab_x + g_0_z_0_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_0_z_y_xxxxxx = cbuffer.data(pi_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_z_y_xxxxxy = cbuffer.data(pi_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_z_y_xxxxxz = cbuffer.data(pi_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_z_y_xxxxyy = cbuffer.data(pi_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_y_xxxxyz = cbuffer.data(pi_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_y_xxxxzz = cbuffer.data(pi_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_y_xxxyyy = cbuffer.data(pi_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_y_xxxyyz = cbuffer.data(pi_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_z_y_xxxyzz = cbuffer.data(pi_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_y_xxxzzz = cbuffer.data(pi_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_y_xxyyyy = cbuffer.data(pi_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_y_xxyyyz = cbuffer.data(pi_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_y_xxyyzz = cbuffer.data(pi_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_y_xxyzzz = cbuffer.data(pi_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_z_y_xxzzzz = cbuffer.data(pi_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_y_xyyyyy = cbuffer.data(pi_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_y_xyyyyz = cbuffer.data(pi_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_y_xyyyzz = cbuffer.data(pi_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_y_xyyzzz = cbuffer.data(pi_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_y_xyzzzz = cbuffer.data(pi_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_y_xzzzzz = cbuffer.data(pi_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_y_yyyyyy = cbuffer.data(pi_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_y_yyyyyz = cbuffer.data(pi_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_y_yyyyzz = cbuffer.data(pi_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_y_yyyzzz = cbuffer.data(pi_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_y_yyzzzz = cbuffer.data(pi_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_y_yzzzzz = cbuffer.data(pi_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_y_zzzzzz = cbuffer.data(pi_geom_01_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxxx, g_0_z_0_xxxxxxy, g_0_z_0_xxxxxy, g_0_z_0_xxxxxyy, g_0_z_0_xxxxxyz, g_0_z_0_xxxxxz, g_0_z_0_xxxxyy, g_0_z_0_xxxxyyy, g_0_z_0_xxxxyyz, g_0_z_0_xxxxyz, g_0_z_0_xxxxyzz, g_0_z_0_xxxxzz, g_0_z_0_xxxyyy, g_0_z_0_xxxyyyy, g_0_z_0_xxxyyyz, g_0_z_0_xxxyyz, g_0_z_0_xxxyyzz, g_0_z_0_xxxyzz, g_0_z_0_xxxyzzz, g_0_z_0_xxxzzz, g_0_z_0_xxyyyy, g_0_z_0_xxyyyyy, g_0_z_0_xxyyyyz, g_0_z_0_xxyyyz, g_0_z_0_xxyyyzz, g_0_z_0_xxyyzz, g_0_z_0_xxyyzzz, g_0_z_0_xxyzzz, g_0_z_0_xxyzzzz, g_0_z_0_xxzzzz, g_0_z_0_xyyyyy, g_0_z_0_xyyyyyy, g_0_z_0_xyyyyyz, g_0_z_0_xyyyyz, g_0_z_0_xyyyyzz, g_0_z_0_xyyyzz, g_0_z_0_xyyyzzz, g_0_z_0_xyyzzz, g_0_z_0_xyyzzzz, g_0_z_0_xyzzzz, g_0_z_0_xyzzzzz, g_0_z_0_xzzzzz, g_0_z_0_yyyyyy, g_0_z_0_yyyyyyy, g_0_z_0_yyyyyyz, g_0_z_0_yyyyyz, g_0_z_0_yyyyyzz, g_0_z_0_yyyyzz, g_0_z_0_yyyyzzz, g_0_z_0_yyyzzz, g_0_z_0_yyyzzzz, g_0_z_0_yyzzzz, g_0_z_0_yyzzzzz, g_0_z_0_yzzzzz, g_0_z_0_yzzzzzz, g_0_z_0_zzzzzz, g_0_z_y_xxxxxx, g_0_z_y_xxxxxy, g_0_z_y_xxxxxz, g_0_z_y_xxxxyy, g_0_z_y_xxxxyz, g_0_z_y_xxxxzz, g_0_z_y_xxxyyy, g_0_z_y_xxxyyz, g_0_z_y_xxxyzz, g_0_z_y_xxxzzz, g_0_z_y_xxyyyy, g_0_z_y_xxyyyz, g_0_z_y_xxyyzz, g_0_z_y_xxyzzz, g_0_z_y_xxzzzz, g_0_z_y_xyyyyy, g_0_z_y_xyyyyz, g_0_z_y_xyyyzz, g_0_z_y_xyyzzz, g_0_z_y_xyzzzz, g_0_z_y_xzzzzz, g_0_z_y_yyyyyy, g_0_z_y_yyyyyz, g_0_z_y_yyyyzz, g_0_z_y_yyyzzz, g_0_z_y_yyzzzz, g_0_z_y_yzzzzz, g_0_z_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_y_xxxxxx[k] = -g_0_z_0_xxxxxx[k] * ab_y + g_0_z_0_xxxxxxy[k];

                g_0_z_y_xxxxxy[k] = -g_0_z_0_xxxxxy[k] * ab_y + g_0_z_0_xxxxxyy[k];

                g_0_z_y_xxxxxz[k] = -g_0_z_0_xxxxxz[k] * ab_y + g_0_z_0_xxxxxyz[k];

                g_0_z_y_xxxxyy[k] = -g_0_z_0_xxxxyy[k] * ab_y + g_0_z_0_xxxxyyy[k];

                g_0_z_y_xxxxyz[k] = -g_0_z_0_xxxxyz[k] * ab_y + g_0_z_0_xxxxyyz[k];

                g_0_z_y_xxxxzz[k] = -g_0_z_0_xxxxzz[k] * ab_y + g_0_z_0_xxxxyzz[k];

                g_0_z_y_xxxyyy[k] = -g_0_z_0_xxxyyy[k] * ab_y + g_0_z_0_xxxyyyy[k];

                g_0_z_y_xxxyyz[k] = -g_0_z_0_xxxyyz[k] * ab_y + g_0_z_0_xxxyyyz[k];

                g_0_z_y_xxxyzz[k] = -g_0_z_0_xxxyzz[k] * ab_y + g_0_z_0_xxxyyzz[k];

                g_0_z_y_xxxzzz[k] = -g_0_z_0_xxxzzz[k] * ab_y + g_0_z_0_xxxyzzz[k];

                g_0_z_y_xxyyyy[k] = -g_0_z_0_xxyyyy[k] * ab_y + g_0_z_0_xxyyyyy[k];

                g_0_z_y_xxyyyz[k] = -g_0_z_0_xxyyyz[k] * ab_y + g_0_z_0_xxyyyyz[k];

                g_0_z_y_xxyyzz[k] = -g_0_z_0_xxyyzz[k] * ab_y + g_0_z_0_xxyyyzz[k];

                g_0_z_y_xxyzzz[k] = -g_0_z_0_xxyzzz[k] * ab_y + g_0_z_0_xxyyzzz[k];

                g_0_z_y_xxzzzz[k] = -g_0_z_0_xxzzzz[k] * ab_y + g_0_z_0_xxyzzzz[k];

                g_0_z_y_xyyyyy[k] = -g_0_z_0_xyyyyy[k] * ab_y + g_0_z_0_xyyyyyy[k];

                g_0_z_y_xyyyyz[k] = -g_0_z_0_xyyyyz[k] * ab_y + g_0_z_0_xyyyyyz[k];

                g_0_z_y_xyyyzz[k] = -g_0_z_0_xyyyzz[k] * ab_y + g_0_z_0_xyyyyzz[k];

                g_0_z_y_xyyzzz[k] = -g_0_z_0_xyyzzz[k] * ab_y + g_0_z_0_xyyyzzz[k];

                g_0_z_y_xyzzzz[k] = -g_0_z_0_xyzzzz[k] * ab_y + g_0_z_0_xyyzzzz[k];

                g_0_z_y_xzzzzz[k] = -g_0_z_0_xzzzzz[k] * ab_y + g_0_z_0_xyzzzzz[k];

                g_0_z_y_yyyyyy[k] = -g_0_z_0_yyyyyy[k] * ab_y + g_0_z_0_yyyyyyy[k];

                g_0_z_y_yyyyyz[k] = -g_0_z_0_yyyyyz[k] * ab_y + g_0_z_0_yyyyyyz[k];

                g_0_z_y_yyyyzz[k] = -g_0_z_0_yyyyzz[k] * ab_y + g_0_z_0_yyyyyzz[k];

                g_0_z_y_yyyzzz[k] = -g_0_z_0_yyyzzz[k] * ab_y + g_0_z_0_yyyyzzz[k];

                g_0_z_y_yyzzzz[k] = -g_0_z_0_yyzzzz[k] * ab_y + g_0_z_0_yyyzzzz[k];

                g_0_z_y_yzzzzz[k] = -g_0_z_0_yzzzzz[k] * ab_y + g_0_z_0_yyzzzzz[k];

                g_0_z_y_zzzzzz[k] = -g_0_z_0_zzzzzz[k] * ab_y + g_0_z_0_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_0_z_z_xxxxxx = cbuffer.data(pi_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_z_xxxxxy = cbuffer.data(pi_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_z_xxxxxz = cbuffer.data(pi_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_z_xxxxyy = cbuffer.data(pi_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_z_xxxxyz = cbuffer.data(pi_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_z_xxxxzz = cbuffer.data(pi_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_z_xxxyyy = cbuffer.data(pi_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_z_xxxyyz = cbuffer.data(pi_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_z_xxxyzz = cbuffer.data(pi_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_z_xxxzzz = cbuffer.data(pi_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_z_xxyyyy = cbuffer.data(pi_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_z_xxyyyz = cbuffer.data(pi_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_z_xxyyzz = cbuffer.data(pi_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_z_xxyzzz = cbuffer.data(pi_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_z_xxzzzz = cbuffer.data(pi_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_z_xyyyyy = cbuffer.data(pi_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_z_xyyyyz = cbuffer.data(pi_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_z_xyyyzz = cbuffer.data(pi_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_z_xyyzzz = cbuffer.data(pi_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_z_xyzzzz = cbuffer.data(pi_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_z_xzzzzz = cbuffer.data(pi_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_z_yyyyyy = cbuffer.data(pi_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_z_yyyyyz = cbuffer.data(pi_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_z_yyyyzz = cbuffer.data(pi_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_z_yyyzzz = cbuffer.data(pi_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_z_yyzzzz = cbuffer.data(pi_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_z_yzzzzz = cbuffer.data(pi_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_z_zzzzzz = cbuffer.data(pi_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_xxxxxx, g_0_xxxxxy, g_0_xxxxxz, g_0_xxxxyy, g_0_xxxxyz, g_0_xxxxzz, g_0_xxxyyy, g_0_xxxyyz, g_0_xxxyzz, g_0_xxxzzz, g_0_xxyyyy, g_0_xxyyyz, g_0_xxyyzz, g_0_xxyzzz, g_0_xxzzzz, g_0_xyyyyy, g_0_xyyyyz, g_0_xyyyzz, g_0_xyyzzz, g_0_xyzzzz, g_0_xzzzzz, g_0_yyyyyy, g_0_yyyyyz, g_0_yyyyzz, g_0_yyyzzz, g_0_yyzzzz, g_0_yzzzzz, g_0_z_0_xxxxxx, g_0_z_0_xxxxxxz, g_0_z_0_xxxxxy, g_0_z_0_xxxxxyz, g_0_z_0_xxxxxz, g_0_z_0_xxxxxzz, g_0_z_0_xxxxyy, g_0_z_0_xxxxyyz, g_0_z_0_xxxxyz, g_0_z_0_xxxxyzz, g_0_z_0_xxxxzz, g_0_z_0_xxxxzzz, g_0_z_0_xxxyyy, g_0_z_0_xxxyyyz, g_0_z_0_xxxyyz, g_0_z_0_xxxyyzz, g_0_z_0_xxxyzz, g_0_z_0_xxxyzzz, g_0_z_0_xxxzzz, g_0_z_0_xxxzzzz, g_0_z_0_xxyyyy, g_0_z_0_xxyyyyz, g_0_z_0_xxyyyz, g_0_z_0_xxyyyzz, g_0_z_0_xxyyzz, g_0_z_0_xxyyzzz, g_0_z_0_xxyzzz, g_0_z_0_xxyzzzz, g_0_z_0_xxzzzz, g_0_z_0_xxzzzzz, g_0_z_0_xyyyyy, g_0_z_0_xyyyyyz, g_0_z_0_xyyyyz, g_0_z_0_xyyyyzz, g_0_z_0_xyyyzz, g_0_z_0_xyyyzzz, g_0_z_0_xyyzzz, g_0_z_0_xyyzzzz, g_0_z_0_xyzzzz, g_0_z_0_xyzzzzz, g_0_z_0_xzzzzz, g_0_z_0_xzzzzzz, g_0_z_0_yyyyyy, g_0_z_0_yyyyyyz, g_0_z_0_yyyyyz, g_0_z_0_yyyyyzz, g_0_z_0_yyyyzz, g_0_z_0_yyyyzzz, g_0_z_0_yyyzzz, g_0_z_0_yyyzzzz, g_0_z_0_yyzzzz, g_0_z_0_yyzzzzz, g_0_z_0_yzzzzz, g_0_z_0_yzzzzzz, g_0_z_0_zzzzzz, g_0_z_0_zzzzzzz, g_0_z_z_xxxxxx, g_0_z_z_xxxxxy, g_0_z_z_xxxxxz, g_0_z_z_xxxxyy, g_0_z_z_xxxxyz, g_0_z_z_xxxxzz, g_0_z_z_xxxyyy, g_0_z_z_xxxyyz, g_0_z_z_xxxyzz, g_0_z_z_xxxzzz, g_0_z_z_xxyyyy, g_0_z_z_xxyyyz, g_0_z_z_xxyyzz, g_0_z_z_xxyzzz, g_0_z_z_xxzzzz, g_0_z_z_xyyyyy, g_0_z_z_xyyyyz, g_0_z_z_xyyyzz, g_0_z_z_xyyzzz, g_0_z_z_xyzzzz, g_0_z_z_xzzzzz, g_0_z_z_yyyyyy, g_0_z_z_yyyyyz, g_0_z_z_yyyyzz, g_0_z_z_yyyzzz, g_0_z_z_yyzzzz, g_0_z_z_yzzzzz, g_0_z_z_zzzzzz, g_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_z_xxxxxx[k] = g_0_xxxxxx[k] - g_0_z_0_xxxxxx[k] * ab_z + g_0_z_0_xxxxxxz[k];

                g_0_z_z_xxxxxy[k] = g_0_xxxxxy[k] - g_0_z_0_xxxxxy[k] * ab_z + g_0_z_0_xxxxxyz[k];

                g_0_z_z_xxxxxz[k] = g_0_xxxxxz[k] - g_0_z_0_xxxxxz[k] * ab_z + g_0_z_0_xxxxxzz[k];

                g_0_z_z_xxxxyy[k] = g_0_xxxxyy[k] - g_0_z_0_xxxxyy[k] * ab_z + g_0_z_0_xxxxyyz[k];

                g_0_z_z_xxxxyz[k] = g_0_xxxxyz[k] - g_0_z_0_xxxxyz[k] * ab_z + g_0_z_0_xxxxyzz[k];

                g_0_z_z_xxxxzz[k] = g_0_xxxxzz[k] - g_0_z_0_xxxxzz[k] * ab_z + g_0_z_0_xxxxzzz[k];

                g_0_z_z_xxxyyy[k] = g_0_xxxyyy[k] - g_0_z_0_xxxyyy[k] * ab_z + g_0_z_0_xxxyyyz[k];

                g_0_z_z_xxxyyz[k] = g_0_xxxyyz[k] - g_0_z_0_xxxyyz[k] * ab_z + g_0_z_0_xxxyyzz[k];

                g_0_z_z_xxxyzz[k] = g_0_xxxyzz[k] - g_0_z_0_xxxyzz[k] * ab_z + g_0_z_0_xxxyzzz[k];

                g_0_z_z_xxxzzz[k] = g_0_xxxzzz[k] - g_0_z_0_xxxzzz[k] * ab_z + g_0_z_0_xxxzzzz[k];

                g_0_z_z_xxyyyy[k] = g_0_xxyyyy[k] - g_0_z_0_xxyyyy[k] * ab_z + g_0_z_0_xxyyyyz[k];

                g_0_z_z_xxyyyz[k] = g_0_xxyyyz[k] - g_0_z_0_xxyyyz[k] * ab_z + g_0_z_0_xxyyyzz[k];

                g_0_z_z_xxyyzz[k] = g_0_xxyyzz[k] - g_0_z_0_xxyyzz[k] * ab_z + g_0_z_0_xxyyzzz[k];

                g_0_z_z_xxyzzz[k] = g_0_xxyzzz[k] - g_0_z_0_xxyzzz[k] * ab_z + g_0_z_0_xxyzzzz[k];

                g_0_z_z_xxzzzz[k] = g_0_xxzzzz[k] - g_0_z_0_xxzzzz[k] * ab_z + g_0_z_0_xxzzzzz[k];

                g_0_z_z_xyyyyy[k] = g_0_xyyyyy[k] - g_0_z_0_xyyyyy[k] * ab_z + g_0_z_0_xyyyyyz[k];

                g_0_z_z_xyyyyz[k] = g_0_xyyyyz[k] - g_0_z_0_xyyyyz[k] * ab_z + g_0_z_0_xyyyyzz[k];

                g_0_z_z_xyyyzz[k] = g_0_xyyyzz[k] - g_0_z_0_xyyyzz[k] * ab_z + g_0_z_0_xyyyzzz[k];

                g_0_z_z_xyyzzz[k] = g_0_xyyzzz[k] - g_0_z_0_xyyzzz[k] * ab_z + g_0_z_0_xyyzzzz[k];

                g_0_z_z_xyzzzz[k] = g_0_xyzzzz[k] - g_0_z_0_xyzzzz[k] * ab_z + g_0_z_0_xyzzzzz[k];

                g_0_z_z_xzzzzz[k] = g_0_xzzzzz[k] - g_0_z_0_xzzzzz[k] * ab_z + g_0_z_0_xzzzzzz[k];

                g_0_z_z_yyyyyy[k] = g_0_yyyyyy[k] - g_0_z_0_yyyyyy[k] * ab_z + g_0_z_0_yyyyyyz[k];

                g_0_z_z_yyyyyz[k] = g_0_yyyyyz[k] - g_0_z_0_yyyyyz[k] * ab_z + g_0_z_0_yyyyyzz[k];

                g_0_z_z_yyyyzz[k] = g_0_yyyyzz[k] - g_0_z_0_yyyyzz[k] * ab_z + g_0_z_0_yyyyzzz[k];

                g_0_z_z_yyyzzz[k] = g_0_yyyzzz[k] - g_0_z_0_yyyzzz[k] * ab_z + g_0_z_0_yyyzzzz[k];

                g_0_z_z_yyzzzz[k] = g_0_yyzzzz[k] - g_0_z_0_yyzzzz[k] * ab_z + g_0_z_0_yyzzzzz[k];

                g_0_z_z_yzzzzz[k] = g_0_yzzzzz[k] - g_0_z_0_yzzzzz[k] * ab_z + g_0_z_0_yzzzzzz[k];

                g_0_z_z_zzzzzz[k] = g_0_zzzzzz[k] - g_0_z_0_zzzzzz[k] * ab_z + g_0_z_0_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

