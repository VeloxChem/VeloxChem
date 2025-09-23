#include "ElectronRepulsionGeom2000ContrRecPIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_pixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_pixx,
                                            const size_t idx_geom_10_sixx,
                                            const size_t idx_geom_20_sixx,
                                            const size_t idx_geom_20_skxx,
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

            const auto si_geom_10_off = idx_geom_10_sixx + i * dcomps + j;

            auto g_x_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 27 * ccomps * dcomps);

            auto g_y_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 28 * ccomps * dcomps);

            auto g_y_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 29 * ccomps * dcomps);

            auto g_y_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 30 * ccomps * dcomps);

            auto g_y_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 31 * ccomps * dcomps);

            auto g_y_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 32 * ccomps * dcomps);

            auto g_y_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 33 * ccomps * dcomps);

            auto g_y_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 34 * ccomps * dcomps);

            auto g_y_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 35 * ccomps * dcomps);

            auto g_y_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 36 * ccomps * dcomps);

            auto g_y_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 37 * ccomps * dcomps);

            auto g_y_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 38 * ccomps * dcomps);

            auto g_y_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 39 * ccomps * dcomps);

            auto g_y_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 40 * ccomps * dcomps);

            auto g_y_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 41 * ccomps * dcomps);

            auto g_y_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 42 * ccomps * dcomps);

            auto g_y_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 43 * ccomps * dcomps);

            auto g_y_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 44 * ccomps * dcomps);

            auto g_y_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 45 * ccomps * dcomps);

            auto g_y_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 46 * ccomps * dcomps);

            auto g_y_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 47 * ccomps * dcomps);

            auto g_y_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 48 * ccomps * dcomps);

            auto g_y_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 49 * ccomps * dcomps);

            auto g_y_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 50 * ccomps * dcomps);

            auto g_y_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 51 * ccomps * dcomps);

            auto g_y_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 52 * ccomps * dcomps);

            auto g_y_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 53 * ccomps * dcomps);

            auto g_y_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 54 * ccomps * dcomps);

            auto g_y_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 55 * ccomps * dcomps);

            auto g_z_0_0_xxxxxx = cbuffer.data(si_geom_10_off + 56 * ccomps * dcomps);

            auto g_z_0_0_xxxxxy = cbuffer.data(si_geom_10_off + 57 * ccomps * dcomps);

            auto g_z_0_0_xxxxxz = cbuffer.data(si_geom_10_off + 58 * ccomps * dcomps);

            auto g_z_0_0_xxxxyy = cbuffer.data(si_geom_10_off + 59 * ccomps * dcomps);

            auto g_z_0_0_xxxxyz = cbuffer.data(si_geom_10_off + 60 * ccomps * dcomps);

            auto g_z_0_0_xxxxzz = cbuffer.data(si_geom_10_off + 61 * ccomps * dcomps);

            auto g_z_0_0_xxxyyy = cbuffer.data(si_geom_10_off + 62 * ccomps * dcomps);

            auto g_z_0_0_xxxyyz = cbuffer.data(si_geom_10_off + 63 * ccomps * dcomps);

            auto g_z_0_0_xxxyzz = cbuffer.data(si_geom_10_off + 64 * ccomps * dcomps);

            auto g_z_0_0_xxxzzz = cbuffer.data(si_geom_10_off + 65 * ccomps * dcomps);

            auto g_z_0_0_xxyyyy = cbuffer.data(si_geom_10_off + 66 * ccomps * dcomps);

            auto g_z_0_0_xxyyyz = cbuffer.data(si_geom_10_off + 67 * ccomps * dcomps);

            auto g_z_0_0_xxyyzz = cbuffer.data(si_geom_10_off + 68 * ccomps * dcomps);

            auto g_z_0_0_xxyzzz = cbuffer.data(si_geom_10_off + 69 * ccomps * dcomps);

            auto g_z_0_0_xxzzzz = cbuffer.data(si_geom_10_off + 70 * ccomps * dcomps);

            auto g_z_0_0_xyyyyy = cbuffer.data(si_geom_10_off + 71 * ccomps * dcomps);

            auto g_z_0_0_xyyyyz = cbuffer.data(si_geom_10_off + 72 * ccomps * dcomps);

            auto g_z_0_0_xyyyzz = cbuffer.data(si_geom_10_off + 73 * ccomps * dcomps);

            auto g_z_0_0_xyyzzz = cbuffer.data(si_geom_10_off + 74 * ccomps * dcomps);

            auto g_z_0_0_xyzzzz = cbuffer.data(si_geom_10_off + 75 * ccomps * dcomps);

            auto g_z_0_0_xzzzzz = cbuffer.data(si_geom_10_off + 76 * ccomps * dcomps);

            auto g_z_0_0_yyyyyy = cbuffer.data(si_geom_10_off + 77 * ccomps * dcomps);

            auto g_z_0_0_yyyyyz = cbuffer.data(si_geom_10_off + 78 * ccomps * dcomps);

            auto g_z_0_0_yyyyzz = cbuffer.data(si_geom_10_off + 79 * ccomps * dcomps);

            auto g_z_0_0_yyyzzz = cbuffer.data(si_geom_10_off + 80 * ccomps * dcomps);

            auto g_z_0_0_yyzzzz = cbuffer.data(si_geom_10_off + 81 * ccomps * dcomps);

            auto g_z_0_0_yzzzzz = cbuffer.data(si_geom_10_off + 82 * ccomps * dcomps);

            auto g_z_0_0_zzzzzz = cbuffer.data(si_geom_10_off + 83 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SISS

            const auto si_geom_20_off = idx_geom_20_sixx + i * dcomps + j;

            auto g_xx_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 27 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 28 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 29 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 30 * ccomps * dcomps);

            auto g_xy_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 31 * ccomps * dcomps);

            auto g_xy_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 32 * ccomps * dcomps);

            auto g_xy_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 33 * ccomps * dcomps);

            auto g_xy_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 34 * ccomps * dcomps);

            auto g_xy_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 35 * ccomps * dcomps);

            auto g_xy_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 36 * ccomps * dcomps);

            auto g_xy_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 37 * ccomps * dcomps);

            auto g_xy_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 38 * ccomps * dcomps);

            auto g_xy_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 39 * ccomps * dcomps);

            auto g_xy_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 40 * ccomps * dcomps);

            auto g_xy_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 41 * ccomps * dcomps);

            auto g_xy_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 42 * ccomps * dcomps);

            auto g_xy_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 43 * ccomps * dcomps);

            auto g_xy_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 44 * ccomps * dcomps);

            auto g_xy_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 47 * ccomps * dcomps);

            auto g_xy_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 50 * ccomps * dcomps);

            auto g_xy_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 53 * ccomps * dcomps);

            auto g_xy_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 55 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 56 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 57 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 58 * ccomps * dcomps);

            auto g_xz_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 59 * ccomps * dcomps);

            auto g_xz_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 60 * ccomps * dcomps);

            auto g_xz_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 61 * ccomps * dcomps);

            auto g_xz_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 62 * ccomps * dcomps);

            auto g_xz_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 63 * ccomps * dcomps);

            auto g_xz_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 64 * ccomps * dcomps);

            auto g_xz_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 65 * ccomps * dcomps);

            auto g_xz_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 66 * ccomps * dcomps);

            auto g_xz_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 67 * ccomps * dcomps);

            auto g_xz_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 68 * ccomps * dcomps);

            auto g_xz_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 69 * ccomps * dcomps);

            auto g_xz_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 70 * ccomps * dcomps);

            auto g_xz_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 71 * ccomps * dcomps);

            auto g_xz_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 72 * ccomps * dcomps);

            auto g_xz_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 73 * ccomps * dcomps);

            auto g_xz_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 74 * ccomps * dcomps);

            auto g_xz_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 75 * ccomps * dcomps);

            auto g_xz_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 76 * ccomps * dcomps);

            auto g_xz_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 77 * ccomps * dcomps);

            auto g_xz_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 78 * ccomps * dcomps);

            auto g_xz_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 79 * ccomps * dcomps);

            auto g_xz_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 80 * ccomps * dcomps);

            auto g_xz_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 81 * ccomps * dcomps);

            auto g_xz_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 82 * ccomps * dcomps);

            auto g_xz_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 83 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 84 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 85 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 86 * ccomps * dcomps);

            auto g_yy_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 87 * ccomps * dcomps);

            auto g_yy_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 88 * ccomps * dcomps);

            auto g_yy_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 89 * ccomps * dcomps);

            auto g_yy_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 90 * ccomps * dcomps);

            auto g_yy_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 91 * ccomps * dcomps);

            auto g_yy_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 92 * ccomps * dcomps);

            auto g_yy_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 93 * ccomps * dcomps);

            auto g_yy_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 94 * ccomps * dcomps);

            auto g_yy_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 95 * ccomps * dcomps);

            auto g_yy_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 96 * ccomps * dcomps);

            auto g_yy_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 97 * ccomps * dcomps);

            auto g_yy_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 98 * ccomps * dcomps);

            auto g_yy_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 99 * ccomps * dcomps);

            auto g_yy_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 100 * ccomps * dcomps);

            auto g_yy_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 101 * ccomps * dcomps);

            auto g_yy_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 102 * ccomps * dcomps);

            auto g_yy_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 103 * ccomps * dcomps);

            auto g_yy_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 104 * ccomps * dcomps);

            auto g_yy_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 105 * ccomps * dcomps);

            auto g_yy_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 106 * ccomps * dcomps);

            auto g_yy_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 107 * ccomps * dcomps);

            auto g_yy_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 108 * ccomps * dcomps);

            auto g_yy_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 109 * ccomps * dcomps);

            auto g_yy_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 110 * ccomps * dcomps);

            auto g_yy_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 111 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 112 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 113 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 114 * ccomps * dcomps);

            auto g_yz_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 115 * ccomps * dcomps);

            auto g_yz_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 116 * ccomps * dcomps);

            auto g_yz_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 117 * ccomps * dcomps);

            auto g_yz_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 118 * ccomps * dcomps);

            auto g_yz_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 119 * ccomps * dcomps);

            auto g_yz_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 120 * ccomps * dcomps);

            auto g_yz_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 121 * ccomps * dcomps);

            auto g_yz_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 122 * ccomps * dcomps);

            auto g_yz_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 123 * ccomps * dcomps);

            auto g_yz_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 124 * ccomps * dcomps);

            auto g_yz_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 125 * ccomps * dcomps);

            auto g_yz_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 126 * ccomps * dcomps);

            auto g_yz_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 127 * ccomps * dcomps);

            auto g_yz_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 128 * ccomps * dcomps);

            auto g_yz_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 129 * ccomps * dcomps);

            auto g_yz_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 130 * ccomps * dcomps);

            auto g_yz_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 131 * ccomps * dcomps);

            auto g_yz_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 132 * ccomps * dcomps);

            auto g_yz_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 133 * ccomps * dcomps);

            auto g_yz_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 134 * ccomps * dcomps);

            auto g_yz_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 135 * ccomps * dcomps);

            auto g_yz_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 136 * ccomps * dcomps);

            auto g_yz_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 137 * ccomps * dcomps);

            auto g_yz_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 138 * ccomps * dcomps);

            auto g_yz_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 139 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxx = cbuffer.data(si_geom_20_off + 140 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxy = cbuffer.data(si_geom_20_off + 141 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxz = cbuffer.data(si_geom_20_off + 142 * ccomps * dcomps);

            auto g_zz_0_0_xxxxyy = cbuffer.data(si_geom_20_off + 143 * ccomps * dcomps);

            auto g_zz_0_0_xxxxyz = cbuffer.data(si_geom_20_off + 144 * ccomps * dcomps);

            auto g_zz_0_0_xxxxzz = cbuffer.data(si_geom_20_off + 145 * ccomps * dcomps);

            auto g_zz_0_0_xxxyyy = cbuffer.data(si_geom_20_off + 146 * ccomps * dcomps);

            auto g_zz_0_0_xxxyyz = cbuffer.data(si_geom_20_off + 147 * ccomps * dcomps);

            auto g_zz_0_0_xxxyzz = cbuffer.data(si_geom_20_off + 148 * ccomps * dcomps);

            auto g_zz_0_0_xxxzzz = cbuffer.data(si_geom_20_off + 149 * ccomps * dcomps);

            auto g_zz_0_0_xxyyyy = cbuffer.data(si_geom_20_off + 150 * ccomps * dcomps);

            auto g_zz_0_0_xxyyyz = cbuffer.data(si_geom_20_off + 151 * ccomps * dcomps);

            auto g_zz_0_0_xxyyzz = cbuffer.data(si_geom_20_off + 152 * ccomps * dcomps);

            auto g_zz_0_0_xxyzzz = cbuffer.data(si_geom_20_off + 153 * ccomps * dcomps);

            auto g_zz_0_0_xxzzzz = cbuffer.data(si_geom_20_off + 154 * ccomps * dcomps);

            auto g_zz_0_0_xyyyyy = cbuffer.data(si_geom_20_off + 155 * ccomps * dcomps);

            auto g_zz_0_0_xyyyyz = cbuffer.data(si_geom_20_off + 156 * ccomps * dcomps);

            auto g_zz_0_0_xyyyzz = cbuffer.data(si_geom_20_off + 157 * ccomps * dcomps);

            auto g_zz_0_0_xyyzzz = cbuffer.data(si_geom_20_off + 158 * ccomps * dcomps);

            auto g_zz_0_0_xyzzzz = cbuffer.data(si_geom_20_off + 159 * ccomps * dcomps);

            auto g_zz_0_0_xzzzzz = cbuffer.data(si_geom_20_off + 160 * ccomps * dcomps);

            auto g_zz_0_0_yyyyyy = cbuffer.data(si_geom_20_off + 161 * ccomps * dcomps);

            auto g_zz_0_0_yyyyyz = cbuffer.data(si_geom_20_off + 162 * ccomps * dcomps);

            auto g_zz_0_0_yyyyzz = cbuffer.data(si_geom_20_off + 163 * ccomps * dcomps);

            auto g_zz_0_0_yyyzzz = cbuffer.data(si_geom_20_off + 164 * ccomps * dcomps);

            auto g_zz_0_0_yyzzzz = cbuffer.data(si_geom_20_off + 165 * ccomps * dcomps);

            auto g_zz_0_0_yzzzzz = cbuffer.data(si_geom_20_off + 166 * ccomps * dcomps);

            auto g_zz_0_0_zzzzzz = cbuffer.data(si_geom_20_off + 167 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SKSS

            const auto sk_geom_20_off = idx_geom_20_skxx + i * dcomps + j;

            auto g_xx_0_0_xxxxxxx = cbuffer.data(sk_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_0_xxxxxxy = cbuffer.data(sk_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_0_xxxxxxz = cbuffer.data(sk_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_0_xxxxxyy = cbuffer.data(sk_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_0_xxxxxyz = cbuffer.data(sk_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_0_xxxxxzz = cbuffer.data(sk_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_0_xxxxyyy = cbuffer.data(sk_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_0_xxxxyyz = cbuffer.data(sk_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_0_xxxxyzz = cbuffer.data(sk_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_0_xxxxzzz = cbuffer.data(sk_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_0_xxxyyyy = cbuffer.data(sk_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_0_xxxyyyz = cbuffer.data(sk_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_0_xxxyyzz = cbuffer.data(sk_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_0_xxxyzzz = cbuffer.data(sk_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_0_xxxzzzz = cbuffer.data(sk_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_0_xxyyyyy = cbuffer.data(sk_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_0_xxyyyyz = cbuffer.data(sk_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_0_xxyyyzz = cbuffer.data(sk_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_0_xxyyzzz = cbuffer.data(sk_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_0_xxyzzzz = cbuffer.data(sk_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_0_xxzzzzz = cbuffer.data(sk_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_0_xyyyyyy = cbuffer.data(sk_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_0_xyyyyyz = cbuffer.data(sk_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_0_xyyyyzz = cbuffer.data(sk_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_0_xyyyzzz = cbuffer.data(sk_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_0_xyyzzzz = cbuffer.data(sk_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_0_xyzzzzz = cbuffer.data(sk_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_0_xzzzzzz = cbuffer.data(sk_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_0_yyyyyyy = cbuffer.data(sk_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_0_yyyyyyz = cbuffer.data(sk_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_0_yyyyyzz = cbuffer.data(sk_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_0_yyyyzzz = cbuffer.data(sk_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_0_yyyzzzz = cbuffer.data(sk_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_0_yyzzzzz = cbuffer.data(sk_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_0_yzzzzzz = cbuffer.data(sk_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_0_zzzzzzz = cbuffer.data(sk_geom_20_off + 35 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxxx = cbuffer.data(sk_geom_20_off + 36 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxxy = cbuffer.data(sk_geom_20_off + 37 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxxz = cbuffer.data(sk_geom_20_off + 38 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxyy = cbuffer.data(sk_geom_20_off + 39 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxyz = cbuffer.data(sk_geom_20_off + 40 * ccomps * dcomps);

            auto g_xy_0_0_xxxxxzz = cbuffer.data(sk_geom_20_off + 41 * ccomps * dcomps);

            auto g_xy_0_0_xxxxyyy = cbuffer.data(sk_geom_20_off + 42 * ccomps * dcomps);

            auto g_xy_0_0_xxxxyyz = cbuffer.data(sk_geom_20_off + 43 * ccomps * dcomps);

            auto g_xy_0_0_xxxxyzz = cbuffer.data(sk_geom_20_off + 44 * ccomps * dcomps);

            auto g_xy_0_0_xxxxzzz = cbuffer.data(sk_geom_20_off + 45 * ccomps * dcomps);

            auto g_xy_0_0_xxxyyyy = cbuffer.data(sk_geom_20_off + 46 * ccomps * dcomps);

            auto g_xy_0_0_xxxyyyz = cbuffer.data(sk_geom_20_off + 47 * ccomps * dcomps);

            auto g_xy_0_0_xxxyyzz = cbuffer.data(sk_geom_20_off + 48 * ccomps * dcomps);

            auto g_xy_0_0_xxxyzzz = cbuffer.data(sk_geom_20_off + 49 * ccomps * dcomps);

            auto g_xy_0_0_xxxzzzz = cbuffer.data(sk_geom_20_off + 50 * ccomps * dcomps);

            auto g_xy_0_0_xxyyyyy = cbuffer.data(sk_geom_20_off + 51 * ccomps * dcomps);

            auto g_xy_0_0_xxyyyyz = cbuffer.data(sk_geom_20_off + 52 * ccomps * dcomps);

            auto g_xy_0_0_xxyyyzz = cbuffer.data(sk_geom_20_off + 53 * ccomps * dcomps);

            auto g_xy_0_0_xxyyzzz = cbuffer.data(sk_geom_20_off + 54 * ccomps * dcomps);

            auto g_xy_0_0_xxyzzzz = cbuffer.data(sk_geom_20_off + 55 * ccomps * dcomps);

            auto g_xy_0_0_xxzzzzz = cbuffer.data(sk_geom_20_off + 56 * ccomps * dcomps);

            auto g_xy_0_0_xyyyyyy = cbuffer.data(sk_geom_20_off + 57 * ccomps * dcomps);

            auto g_xy_0_0_xyyyyyz = cbuffer.data(sk_geom_20_off + 58 * ccomps * dcomps);

            auto g_xy_0_0_xyyyyzz = cbuffer.data(sk_geom_20_off + 59 * ccomps * dcomps);

            auto g_xy_0_0_xyyyzzz = cbuffer.data(sk_geom_20_off + 60 * ccomps * dcomps);

            auto g_xy_0_0_xyyzzzz = cbuffer.data(sk_geom_20_off + 61 * ccomps * dcomps);

            auto g_xy_0_0_xyzzzzz = cbuffer.data(sk_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_0_xzzzzzz = cbuffer.data(sk_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_0_yyyyyyy = cbuffer.data(sk_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_0_yyyyyyz = cbuffer.data(sk_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_0_yyyyyzz = cbuffer.data(sk_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_0_yyyyzzz = cbuffer.data(sk_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_0_yyyzzzz = cbuffer.data(sk_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_0_yyzzzzz = cbuffer.data(sk_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_0_yzzzzzz = cbuffer.data(sk_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_0_zzzzzzz = cbuffer.data(sk_geom_20_off + 71 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxxx = cbuffer.data(sk_geom_20_off + 72 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxxy = cbuffer.data(sk_geom_20_off + 73 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxxz = cbuffer.data(sk_geom_20_off + 74 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxyy = cbuffer.data(sk_geom_20_off + 75 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxyz = cbuffer.data(sk_geom_20_off + 76 * ccomps * dcomps);

            auto g_xz_0_0_xxxxxzz = cbuffer.data(sk_geom_20_off + 77 * ccomps * dcomps);

            auto g_xz_0_0_xxxxyyy = cbuffer.data(sk_geom_20_off + 78 * ccomps * dcomps);

            auto g_xz_0_0_xxxxyyz = cbuffer.data(sk_geom_20_off + 79 * ccomps * dcomps);

            auto g_xz_0_0_xxxxyzz = cbuffer.data(sk_geom_20_off + 80 * ccomps * dcomps);

            auto g_xz_0_0_xxxxzzz = cbuffer.data(sk_geom_20_off + 81 * ccomps * dcomps);

            auto g_xz_0_0_xxxyyyy = cbuffer.data(sk_geom_20_off + 82 * ccomps * dcomps);

            auto g_xz_0_0_xxxyyyz = cbuffer.data(sk_geom_20_off + 83 * ccomps * dcomps);

            auto g_xz_0_0_xxxyyzz = cbuffer.data(sk_geom_20_off + 84 * ccomps * dcomps);

            auto g_xz_0_0_xxxyzzz = cbuffer.data(sk_geom_20_off + 85 * ccomps * dcomps);

            auto g_xz_0_0_xxxzzzz = cbuffer.data(sk_geom_20_off + 86 * ccomps * dcomps);

            auto g_xz_0_0_xxyyyyy = cbuffer.data(sk_geom_20_off + 87 * ccomps * dcomps);

            auto g_xz_0_0_xxyyyyz = cbuffer.data(sk_geom_20_off + 88 * ccomps * dcomps);

            auto g_xz_0_0_xxyyyzz = cbuffer.data(sk_geom_20_off + 89 * ccomps * dcomps);

            auto g_xz_0_0_xxyyzzz = cbuffer.data(sk_geom_20_off + 90 * ccomps * dcomps);

            auto g_xz_0_0_xxyzzzz = cbuffer.data(sk_geom_20_off + 91 * ccomps * dcomps);

            auto g_xz_0_0_xxzzzzz = cbuffer.data(sk_geom_20_off + 92 * ccomps * dcomps);

            auto g_xz_0_0_xyyyyyy = cbuffer.data(sk_geom_20_off + 93 * ccomps * dcomps);

            auto g_xz_0_0_xyyyyyz = cbuffer.data(sk_geom_20_off + 94 * ccomps * dcomps);

            auto g_xz_0_0_xyyyyzz = cbuffer.data(sk_geom_20_off + 95 * ccomps * dcomps);

            auto g_xz_0_0_xyyyzzz = cbuffer.data(sk_geom_20_off + 96 * ccomps * dcomps);

            auto g_xz_0_0_xyyzzzz = cbuffer.data(sk_geom_20_off + 97 * ccomps * dcomps);

            auto g_xz_0_0_xyzzzzz = cbuffer.data(sk_geom_20_off + 98 * ccomps * dcomps);

            auto g_xz_0_0_xzzzzzz = cbuffer.data(sk_geom_20_off + 99 * ccomps * dcomps);

            auto g_xz_0_0_yyyyyyy = cbuffer.data(sk_geom_20_off + 100 * ccomps * dcomps);

            auto g_xz_0_0_yyyyyyz = cbuffer.data(sk_geom_20_off + 101 * ccomps * dcomps);

            auto g_xz_0_0_yyyyyzz = cbuffer.data(sk_geom_20_off + 102 * ccomps * dcomps);

            auto g_xz_0_0_yyyyzzz = cbuffer.data(sk_geom_20_off + 103 * ccomps * dcomps);

            auto g_xz_0_0_yyyzzzz = cbuffer.data(sk_geom_20_off + 104 * ccomps * dcomps);

            auto g_xz_0_0_yyzzzzz = cbuffer.data(sk_geom_20_off + 105 * ccomps * dcomps);

            auto g_xz_0_0_yzzzzzz = cbuffer.data(sk_geom_20_off + 106 * ccomps * dcomps);

            auto g_xz_0_0_zzzzzzz = cbuffer.data(sk_geom_20_off + 107 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxxx = cbuffer.data(sk_geom_20_off + 108 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxxy = cbuffer.data(sk_geom_20_off + 109 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxxz = cbuffer.data(sk_geom_20_off + 110 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxyy = cbuffer.data(sk_geom_20_off + 111 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxyz = cbuffer.data(sk_geom_20_off + 112 * ccomps * dcomps);

            auto g_yy_0_0_xxxxxzz = cbuffer.data(sk_geom_20_off + 113 * ccomps * dcomps);

            auto g_yy_0_0_xxxxyyy = cbuffer.data(sk_geom_20_off + 114 * ccomps * dcomps);

            auto g_yy_0_0_xxxxyyz = cbuffer.data(sk_geom_20_off + 115 * ccomps * dcomps);

            auto g_yy_0_0_xxxxyzz = cbuffer.data(sk_geom_20_off + 116 * ccomps * dcomps);

            auto g_yy_0_0_xxxxzzz = cbuffer.data(sk_geom_20_off + 117 * ccomps * dcomps);

            auto g_yy_0_0_xxxyyyy = cbuffer.data(sk_geom_20_off + 118 * ccomps * dcomps);

            auto g_yy_0_0_xxxyyyz = cbuffer.data(sk_geom_20_off + 119 * ccomps * dcomps);

            auto g_yy_0_0_xxxyyzz = cbuffer.data(sk_geom_20_off + 120 * ccomps * dcomps);

            auto g_yy_0_0_xxxyzzz = cbuffer.data(sk_geom_20_off + 121 * ccomps * dcomps);

            auto g_yy_0_0_xxxzzzz = cbuffer.data(sk_geom_20_off + 122 * ccomps * dcomps);

            auto g_yy_0_0_xxyyyyy = cbuffer.data(sk_geom_20_off + 123 * ccomps * dcomps);

            auto g_yy_0_0_xxyyyyz = cbuffer.data(sk_geom_20_off + 124 * ccomps * dcomps);

            auto g_yy_0_0_xxyyyzz = cbuffer.data(sk_geom_20_off + 125 * ccomps * dcomps);

            auto g_yy_0_0_xxyyzzz = cbuffer.data(sk_geom_20_off + 126 * ccomps * dcomps);

            auto g_yy_0_0_xxyzzzz = cbuffer.data(sk_geom_20_off + 127 * ccomps * dcomps);

            auto g_yy_0_0_xxzzzzz = cbuffer.data(sk_geom_20_off + 128 * ccomps * dcomps);

            auto g_yy_0_0_xyyyyyy = cbuffer.data(sk_geom_20_off + 129 * ccomps * dcomps);

            auto g_yy_0_0_xyyyyyz = cbuffer.data(sk_geom_20_off + 130 * ccomps * dcomps);

            auto g_yy_0_0_xyyyyzz = cbuffer.data(sk_geom_20_off + 131 * ccomps * dcomps);

            auto g_yy_0_0_xyyyzzz = cbuffer.data(sk_geom_20_off + 132 * ccomps * dcomps);

            auto g_yy_0_0_xyyzzzz = cbuffer.data(sk_geom_20_off + 133 * ccomps * dcomps);

            auto g_yy_0_0_xyzzzzz = cbuffer.data(sk_geom_20_off + 134 * ccomps * dcomps);

            auto g_yy_0_0_xzzzzzz = cbuffer.data(sk_geom_20_off + 135 * ccomps * dcomps);

            auto g_yy_0_0_yyyyyyy = cbuffer.data(sk_geom_20_off + 136 * ccomps * dcomps);

            auto g_yy_0_0_yyyyyyz = cbuffer.data(sk_geom_20_off + 137 * ccomps * dcomps);

            auto g_yy_0_0_yyyyyzz = cbuffer.data(sk_geom_20_off + 138 * ccomps * dcomps);

            auto g_yy_0_0_yyyyzzz = cbuffer.data(sk_geom_20_off + 139 * ccomps * dcomps);

            auto g_yy_0_0_yyyzzzz = cbuffer.data(sk_geom_20_off + 140 * ccomps * dcomps);

            auto g_yy_0_0_yyzzzzz = cbuffer.data(sk_geom_20_off + 141 * ccomps * dcomps);

            auto g_yy_0_0_yzzzzzz = cbuffer.data(sk_geom_20_off + 142 * ccomps * dcomps);

            auto g_yy_0_0_zzzzzzz = cbuffer.data(sk_geom_20_off + 143 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxxx = cbuffer.data(sk_geom_20_off + 144 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxxy = cbuffer.data(sk_geom_20_off + 145 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxxz = cbuffer.data(sk_geom_20_off + 146 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxyy = cbuffer.data(sk_geom_20_off + 147 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxyz = cbuffer.data(sk_geom_20_off + 148 * ccomps * dcomps);

            auto g_yz_0_0_xxxxxzz = cbuffer.data(sk_geom_20_off + 149 * ccomps * dcomps);

            auto g_yz_0_0_xxxxyyy = cbuffer.data(sk_geom_20_off + 150 * ccomps * dcomps);

            auto g_yz_0_0_xxxxyyz = cbuffer.data(sk_geom_20_off + 151 * ccomps * dcomps);

            auto g_yz_0_0_xxxxyzz = cbuffer.data(sk_geom_20_off + 152 * ccomps * dcomps);

            auto g_yz_0_0_xxxxzzz = cbuffer.data(sk_geom_20_off + 153 * ccomps * dcomps);

            auto g_yz_0_0_xxxyyyy = cbuffer.data(sk_geom_20_off + 154 * ccomps * dcomps);

            auto g_yz_0_0_xxxyyyz = cbuffer.data(sk_geom_20_off + 155 * ccomps * dcomps);

            auto g_yz_0_0_xxxyyzz = cbuffer.data(sk_geom_20_off + 156 * ccomps * dcomps);

            auto g_yz_0_0_xxxyzzz = cbuffer.data(sk_geom_20_off + 157 * ccomps * dcomps);

            auto g_yz_0_0_xxxzzzz = cbuffer.data(sk_geom_20_off + 158 * ccomps * dcomps);

            auto g_yz_0_0_xxyyyyy = cbuffer.data(sk_geom_20_off + 159 * ccomps * dcomps);

            auto g_yz_0_0_xxyyyyz = cbuffer.data(sk_geom_20_off + 160 * ccomps * dcomps);

            auto g_yz_0_0_xxyyyzz = cbuffer.data(sk_geom_20_off + 161 * ccomps * dcomps);

            auto g_yz_0_0_xxyyzzz = cbuffer.data(sk_geom_20_off + 162 * ccomps * dcomps);

            auto g_yz_0_0_xxyzzzz = cbuffer.data(sk_geom_20_off + 163 * ccomps * dcomps);

            auto g_yz_0_0_xxzzzzz = cbuffer.data(sk_geom_20_off + 164 * ccomps * dcomps);

            auto g_yz_0_0_xyyyyyy = cbuffer.data(sk_geom_20_off + 165 * ccomps * dcomps);

            auto g_yz_0_0_xyyyyyz = cbuffer.data(sk_geom_20_off + 166 * ccomps * dcomps);

            auto g_yz_0_0_xyyyyzz = cbuffer.data(sk_geom_20_off + 167 * ccomps * dcomps);

            auto g_yz_0_0_xyyyzzz = cbuffer.data(sk_geom_20_off + 168 * ccomps * dcomps);

            auto g_yz_0_0_xyyzzzz = cbuffer.data(sk_geom_20_off + 169 * ccomps * dcomps);

            auto g_yz_0_0_xyzzzzz = cbuffer.data(sk_geom_20_off + 170 * ccomps * dcomps);

            auto g_yz_0_0_xzzzzzz = cbuffer.data(sk_geom_20_off + 171 * ccomps * dcomps);

            auto g_yz_0_0_yyyyyyy = cbuffer.data(sk_geom_20_off + 172 * ccomps * dcomps);

            auto g_yz_0_0_yyyyyyz = cbuffer.data(sk_geom_20_off + 173 * ccomps * dcomps);

            auto g_yz_0_0_yyyyyzz = cbuffer.data(sk_geom_20_off + 174 * ccomps * dcomps);

            auto g_yz_0_0_yyyyzzz = cbuffer.data(sk_geom_20_off + 175 * ccomps * dcomps);

            auto g_yz_0_0_yyyzzzz = cbuffer.data(sk_geom_20_off + 176 * ccomps * dcomps);

            auto g_yz_0_0_yyzzzzz = cbuffer.data(sk_geom_20_off + 177 * ccomps * dcomps);

            auto g_yz_0_0_yzzzzzz = cbuffer.data(sk_geom_20_off + 178 * ccomps * dcomps);

            auto g_yz_0_0_zzzzzzz = cbuffer.data(sk_geom_20_off + 179 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxxx = cbuffer.data(sk_geom_20_off + 180 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxxy = cbuffer.data(sk_geom_20_off + 181 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxxz = cbuffer.data(sk_geom_20_off + 182 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxyy = cbuffer.data(sk_geom_20_off + 183 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxyz = cbuffer.data(sk_geom_20_off + 184 * ccomps * dcomps);

            auto g_zz_0_0_xxxxxzz = cbuffer.data(sk_geom_20_off + 185 * ccomps * dcomps);

            auto g_zz_0_0_xxxxyyy = cbuffer.data(sk_geom_20_off + 186 * ccomps * dcomps);

            auto g_zz_0_0_xxxxyyz = cbuffer.data(sk_geom_20_off + 187 * ccomps * dcomps);

            auto g_zz_0_0_xxxxyzz = cbuffer.data(sk_geom_20_off + 188 * ccomps * dcomps);

            auto g_zz_0_0_xxxxzzz = cbuffer.data(sk_geom_20_off + 189 * ccomps * dcomps);

            auto g_zz_0_0_xxxyyyy = cbuffer.data(sk_geom_20_off + 190 * ccomps * dcomps);

            auto g_zz_0_0_xxxyyyz = cbuffer.data(sk_geom_20_off + 191 * ccomps * dcomps);

            auto g_zz_0_0_xxxyyzz = cbuffer.data(sk_geom_20_off + 192 * ccomps * dcomps);

            auto g_zz_0_0_xxxyzzz = cbuffer.data(sk_geom_20_off + 193 * ccomps * dcomps);

            auto g_zz_0_0_xxxzzzz = cbuffer.data(sk_geom_20_off + 194 * ccomps * dcomps);

            auto g_zz_0_0_xxyyyyy = cbuffer.data(sk_geom_20_off + 195 * ccomps * dcomps);

            auto g_zz_0_0_xxyyyyz = cbuffer.data(sk_geom_20_off + 196 * ccomps * dcomps);

            auto g_zz_0_0_xxyyyzz = cbuffer.data(sk_geom_20_off + 197 * ccomps * dcomps);

            auto g_zz_0_0_xxyyzzz = cbuffer.data(sk_geom_20_off + 198 * ccomps * dcomps);

            auto g_zz_0_0_xxyzzzz = cbuffer.data(sk_geom_20_off + 199 * ccomps * dcomps);

            auto g_zz_0_0_xxzzzzz = cbuffer.data(sk_geom_20_off + 200 * ccomps * dcomps);

            auto g_zz_0_0_xyyyyyy = cbuffer.data(sk_geom_20_off + 201 * ccomps * dcomps);

            auto g_zz_0_0_xyyyyyz = cbuffer.data(sk_geom_20_off + 202 * ccomps * dcomps);

            auto g_zz_0_0_xyyyyzz = cbuffer.data(sk_geom_20_off + 203 * ccomps * dcomps);

            auto g_zz_0_0_xyyyzzz = cbuffer.data(sk_geom_20_off + 204 * ccomps * dcomps);

            auto g_zz_0_0_xyyzzzz = cbuffer.data(sk_geom_20_off + 205 * ccomps * dcomps);

            auto g_zz_0_0_xyzzzzz = cbuffer.data(sk_geom_20_off + 206 * ccomps * dcomps);

            auto g_zz_0_0_xzzzzzz = cbuffer.data(sk_geom_20_off + 207 * ccomps * dcomps);

            auto g_zz_0_0_yyyyyyy = cbuffer.data(sk_geom_20_off + 208 * ccomps * dcomps);

            auto g_zz_0_0_yyyyyyz = cbuffer.data(sk_geom_20_off + 209 * ccomps * dcomps);

            auto g_zz_0_0_yyyyyzz = cbuffer.data(sk_geom_20_off + 210 * ccomps * dcomps);

            auto g_zz_0_0_yyyyzzz = cbuffer.data(sk_geom_20_off + 211 * ccomps * dcomps);

            auto g_zz_0_0_yyyzzzz = cbuffer.data(sk_geom_20_off + 212 * ccomps * dcomps);

            auto g_zz_0_0_yyzzzzz = cbuffer.data(sk_geom_20_off + 213 * ccomps * dcomps);

            auto g_zz_0_0_yzzzzzz = cbuffer.data(sk_geom_20_off + 214 * ccomps * dcomps);

            auto g_zz_0_0_zzzzzzz = cbuffer.data(sk_geom_20_off + 215 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pixx

            const auto pi_geom_20_off = idx_geom_20_pixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_xx_0_x_xxxxxx = cbuffer.data(pi_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_x_xxxxxy = cbuffer.data(pi_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_x_xxxxxz = cbuffer.data(pi_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_x_xxxxyy = cbuffer.data(pi_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_x_xxxxyz = cbuffer.data(pi_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_x_xxxxzz = cbuffer.data(pi_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_x_xxxyyy = cbuffer.data(pi_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_x_xxxyyz = cbuffer.data(pi_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_x_xxxyzz = cbuffer.data(pi_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_x_xxxzzz = cbuffer.data(pi_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_x_xxyyyy = cbuffer.data(pi_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_x_xxyyyz = cbuffer.data(pi_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_x_xxyyzz = cbuffer.data(pi_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_x_xxyzzz = cbuffer.data(pi_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_x_xxzzzz = cbuffer.data(pi_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_x_xyyyyy = cbuffer.data(pi_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_x_xyyyyz = cbuffer.data(pi_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_x_xyyyzz = cbuffer.data(pi_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_x_xyyzzz = cbuffer.data(pi_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_x_xyzzzz = cbuffer.data(pi_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_x_xzzzzz = cbuffer.data(pi_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_x_yyyyyy = cbuffer.data(pi_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_x_yyyyyz = cbuffer.data(pi_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_x_yyyyzz = cbuffer.data(pi_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_x_yyyzzz = cbuffer.data(pi_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_x_yyzzzz = cbuffer.data(pi_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_x_yzzzzz = cbuffer.data(pi_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_x_zzzzzz = cbuffer.data(pi_geom_20_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxxz, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxzz, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzzz, g_xx_0_0_xxxxxx, g_xx_0_0_xxxxxxx, g_xx_0_0_xxxxxxy, g_xx_0_0_xxxxxxz, g_xx_0_0_xxxxxy, g_xx_0_0_xxxxxyy, g_xx_0_0_xxxxxyz, g_xx_0_0_xxxxxz, g_xx_0_0_xxxxxzz, g_xx_0_0_xxxxyy, g_xx_0_0_xxxxyyy, g_xx_0_0_xxxxyyz, g_xx_0_0_xxxxyz, g_xx_0_0_xxxxyzz, g_xx_0_0_xxxxzz, g_xx_0_0_xxxxzzz, g_xx_0_0_xxxyyy, g_xx_0_0_xxxyyyy, g_xx_0_0_xxxyyyz, g_xx_0_0_xxxyyz, g_xx_0_0_xxxyyzz, g_xx_0_0_xxxyzz, g_xx_0_0_xxxyzzz, g_xx_0_0_xxxzzz, g_xx_0_0_xxxzzzz, g_xx_0_0_xxyyyy, g_xx_0_0_xxyyyyy, g_xx_0_0_xxyyyyz, g_xx_0_0_xxyyyz, g_xx_0_0_xxyyyzz, g_xx_0_0_xxyyzz, g_xx_0_0_xxyyzzz, g_xx_0_0_xxyzzz, g_xx_0_0_xxyzzzz, g_xx_0_0_xxzzzz, g_xx_0_0_xxzzzzz, g_xx_0_0_xyyyyy, g_xx_0_0_xyyyyyy, g_xx_0_0_xyyyyyz, g_xx_0_0_xyyyyz, g_xx_0_0_xyyyyzz, g_xx_0_0_xyyyzz, g_xx_0_0_xyyyzzz, g_xx_0_0_xyyzzz, g_xx_0_0_xyyzzzz, g_xx_0_0_xyzzzz, g_xx_0_0_xyzzzzz, g_xx_0_0_xzzzzz, g_xx_0_0_xzzzzzz, g_xx_0_0_yyyyyy, g_xx_0_0_yyyyyz, g_xx_0_0_yyyyzz, g_xx_0_0_yyyzzz, g_xx_0_0_yyzzzz, g_xx_0_0_yzzzzz, g_xx_0_0_zzzzzz, g_xx_0_x_xxxxxx, g_xx_0_x_xxxxxy, g_xx_0_x_xxxxxz, g_xx_0_x_xxxxyy, g_xx_0_x_xxxxyz, g_xx_0_x_xxxxzz, g_xx_0_x_xxxyyy, g_xx_0_x_xxxyyz, g_xx_0_x_xxxyzz, g_xx_0_x_xxxzzz, g_xx_0_x_xxyyyy, g_xx_0_x_xxyyyz, g_xx_0_x_xxyyzz, g_xx_0_x_xxyzzz, g_xx_0_x_xxzzzz, g_xx_0_x_xyyyyy, g_xx_0_x_xyyyyz, g_xx_0_x_xyyyzz, g_xx_0_x_xyyzzz, g_xx_0_x_xyzzzz, g_xx_0_x_xzzzzz, g_xx_0_x_yyyyyy, g_xx_0_x_yyyyyz, g_xx_0_x_yyyyzz, g_xx_0_x_yyyzzz, g_xx_0_x_yyzzzz, g_xx_0_x_yzzzzz, g_xx_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_x_xxxxxx[k] = -2.0 * g_x_0_0_xxxxxx[k] - g_xx_0_0_xxxxxx[k] * ab_x + g_xx_0_0_xxxxxxx[k];

                g_xx_0_x_xxxxxy[k] = -2.0 * g_x_0_0_xxxxxy[k] - g_xx_0_0_xxxxxy[k] * ab_x + g_xx_0_0_xxxxxxy[k];

                g_xx_0_x_xxxxxz[k] = -2.0 * g_x_0_0_xxxxxz[k] - g_xx_0_0_xxxxxz[k] * ab_x + g_xx_0_0_xxxxxxz[k];

                g_xx_0_x_xxxxyy[k] = -2.0 * g_x_0_0_xxxxyy[k] - g_xx_0_0_xxxxyy[k] * ab_x + g_xx_0_0_xxxxxyy[k];

                g_xx_0_x_xxxxyz[k] = -2.0 * g_x_0_0_xxxxyz[k] - g_xx_0_0_xxxxyz[k] * ab_x + g_xx_0_0_xxxxxyz[k];

                g_xx_0_x_xxxxzz[k] = -2.0 * g_x_0_0_xxxxzz[k] - g_xx_0_0_xxxxzz[k] * ab_x + g_xx_0_0_xxxxxzz[k];

                g_xx_0_x_xxxyyy[k] = -2.0 * g_x_0_0_xxxyyy[k] - g_xx_0_0_xxxyyy[k] * ab_x + g_xx_0_0_xxxxyyy[k];

                g_xx_0_x_xxxyyz[k] = -2.0 * g_x_0_0_xxxyyz[k] - g_xx_0_0_xxxyyz[k] * ab_x + g_xx_0_0_xxxxyyz[k];

                g_xx_0_x_xxxyzz[k] = -2.0 * g_x_0_0_xxxyzz[k] - g_xx_0_0_xxxyzz[k] * ab_x + g_xx_0_0_xxxxyzz[k];

                g_xx_0_x_xxxzzz[k] = -2.0 * g_x_0_0_xxxzzz[k] - g_xx_0_0_xxxzzz[k] * ab_x + g_xx_0_0_xxxxzzz[k];

                g_xx_0_x_xxyyyy[k] = -2.0 * g_x_0_0_xxyyyy[k] - g_xx_0_0_xxyyyy[k] * ab_x + g_xx_0_0_xxxyyyy[k];

                g_xx_0_x_xxyyyz[k] = -2.0 * g_x_0_0_xxyyyz[k] - g_xx_0_0_xxyyyz[k] * ab_x + g_xx_0_0_xxxyyyz[k];

                g_xx_0_x_xxyyzz[k] = -2.0 * g_x_0_0_xxyyzz[k] - g_xx_0_0_xxyyzz[k] * ab_x + g_xx_0_0_xxxyyzz[k];

                g_xx_0_x_xxyzzz[k] = -2.0 * g_x_0_0_xxyzzz[k] - g_xx_0_0_xxyzzz[k] * ab_x + g_xx_0_0_xxxyzzz[k];

                g_xx_0_x_xxzzzz[k] = -2.0 * g_x_0_0_xxzzzz[k] - g_xx_0_0_xxzzzz[k] * ab_x + g_xx_0_0_xxxzzzz[k];

                g_xx_0_x_xyyyyy[k] = -2.0 * g_x_0_0_xyyyyy[k] - g_xx_0_0_xyyyyy[k] * ab_x + g_xx_0_0_xxyyyyy[k];

                g_xx_0_x_xyyyyz[k] = -2.0 * g_x_0_0_xyyyyz[k] - g_xx_0_0_xyyyyz[k] * ab_x + g_xx_0_0_xxyyyyz[k];

                g_xx_0_x_xyyyzz[k] = -2.0 * g_x_0_0_xyyyzz[k] - g_xx_0_0_xyyyzz[k] * ab_x + g_xx_0_0_xxyyyzz[k];

                g_xx_0_x_xyyzzz[k] = -2.0 * g_x_0_0_xyyzzz[k] - g_xx_0_0_xyyzzz[k] * ab_x + g_xx_0_0_xxyyzzz[k];

                g_xx_0_x_xyzzzz[k] = -2.0 * g_x_0_0_xyzzzz[k] - g_xx_0_0_xyzzzz[k] * ab_x + g_xx_0_0_xxyzzzz[k];

                g_xx_0_x_xzzzzz[k] = -2.0 * g_x_0_0_xzzzzz[k] - g_xx_0_0_xzzzzz[k] * ab_x + g_xx_0_0_xxzzzzz[k];

                g_xx_0_x_yyyyyy[k] = -2.0 * g_x_0_0_yyyyyy[k] - g_xx_0_0_yyyyyy[k] * ab_x + g_xx_0_0_xyyyyyy[k];

                g_xx_0_x_yyyyyz[k] = -2.0 * g_x_0_0_yyyyyz[k] - g_xx_0_0_yyyyyz[k] * ab_x + g_xx_0_0_xyyyyyz[k];

                g_xx_0_x_yyyyzz[k] = -2.0 * g_x_0_0_yyyyzz[k] - g_xx_0_0_yyyyzz[k] * ab_x + g_xx_0_0_xyyyyzz[k];

                g_xx_0_x_yyyzzz[k] = -2.0 * g_x_0_0_yyyzzz[k] - g_xx_0_0_yyyzzz[k] * ab_x + g_xx_0_0_xyyyzzz[k];

                g_xx_0_x_yyzzzz[k] = -2.0 * g_x_0_0_yyzzzz[k] - g_xx_0_0_yyzzzz[k] * ab_x + g_xx_0_0_xyyzzzz[k];

                g_xx_0_x_yzzzzz[k] = -2.0 * g_x_0_0_yzzzzz[k] - g_xx_0_0_yzzzzz[k] * ab_x + g_xx_0_0_xyzzzzz[k];

                g_xx_0_x_zzzzzz[k] = -2.0 * g_x_0_0_zzzzzz[k] - g_xx_0_0_zzzzzz[k] * ab_x + g_xx_0_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_xx_0_y_xxxxxx = cbuffer.data(pi_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_y_xxxxxy = cbuffer.data(pi_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_y_xxxxxz = cbuffer.data(pi_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_y_xxxxyy = cbuffer.data(pi_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_y_xxxxyz = cbuffer.data(pi_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_y_xxxxzz = cbuffer.data(pi_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_y_xxxyyy = cbuffer.data(pi_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_y_xxxyyz = cbuffer.data(pi_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_y_xxxyzz = cbuffer.data(pi_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_y_xxxzzz = cbuffer.data(pi_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_y_xxyyyy = cbuffer.data(pi_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_y_xxyyyz = cbuffer.data(pi_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_y_xxyyzz = cbuffer.data(pi_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_y_xxyzzz = cbuffer.data(pi_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_y_xxzzzz = cbuffer.data(pi_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_y_xyyyyy = cbuffer.data(pi_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_y_xyyyyz = cbuffer.data(pi_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_y_xyyyzz = cbuffer.data(pi_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_y_xyyzzz = cbuffer.data(pi_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_y_xyzzzz = cbuffer.data(pi_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_y_xzzzzz = cbuffer.data(pi_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_y_yyyyyy = cbuffer.data(pi_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_y_yyyyyz = cbuffer.data(pi_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_y_yyyyzz = cbuffer.data(pi_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_y_yyyzzz = cbuffer.data(pi_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_y_yyzzzz = cbuffer.data(pi_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_y_yzzzzz = cbuffer.data(pi_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_y_zzzzzz = cbuffer.data(pi_geom_20_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_0_xxxxxx, g_xx_0_0_xxxxxxy, g_xx_0_0_xxxxxy, g_xx_0_0_xxxxxyy, g_xx_0_0_xxxxxyz, g_xx_0_0_xxxxxz, g_xx_0_0_xxxxyy, g_xx_0_0_xxxxyyy, g_xx_0_0_xxxxyyz, g_xx_0_0_xxxxyz, g_xx_0_0_xxxxyzz, g_xx_0_0_xxxxzz, g_xx_0_0_xxxyyy, g_xx_0_0_xxxyyyy, g_xx_0_0_xxxyyyz, g_xx_0_0_xxxyyz, g_xx_0_0_xxxyyzz, g_xx_0_0_xxxyzz, g_xx_0_0_xxxyzzz, g_xx_0_0_xxxzzz, g_xx_0_0_xxyyyy, g_xx_0_0_xxyyyyy, g_xx_0_0_xxyyyyz, g_xx_0_0_xxyyyz, g_xx_0_0_xxyyyzz, g_xx_0_0_xxyyzz, g_xx_0_0_xxyyzzz, g_xx_0_0_xxyzzz, g_xx_0_0_xxyzzzz, g_xx_0_0_xxzzzz, g_xx_0_0_xyyyyy, g_xx_0_0_xyyyyyy, g_xx_0_0_xyyyyyz, g_xx_0_0_xyyyyz, g_xx_0_0_xyyyyzz, g_xx_0_0_xyyyzz, g_xx_0_0_xyyyzzz, g_xx_0_0_xyyzzz, g_xx_0_0_xyyzzzz, g_xx_0_0_xyzzzz, g_xx_0_0_xyzzzzz, g_xx_0_0_xzzzzz, g_xx_0_0_yyyyyy, g_xx_0_0_yyyyyyy, g_xx_0_0_yyyyyyz, g_xx_0_0_yyyyyz, g_xx_0_0_yyyyyzz, g_xx_0_0_yyyyzz, g_xx_0_0_yyyyzzz, g_xx_0_0_yyyzzz, g_xx_0_0_yyyzzzz, g_xx_0_0_yyzzzz, g_xx_0_0_yyzzzzz, g_xx_0_0_yzzzzz, g_xx_0_0_yzzzzzz, g_xx_0_0_zzzzzz, g_xx_0_y_xxxxxx, g_xx_0_y_xxxxxy, g_xx_0_y_xxxxxz, g_xx_0_y_xxxxyy, g_xx_0_y_xxxxyz, g_xx_0_y_xxxxzz, g_xx_0_y_xxxyyy, g_xx_0_y_xxxyyz, g_xx_0_y_xxxyzz, g_xx_0_y_xxxzzz, g_xx_0_y_xxyyyy, g_xx_0_y_xxyyyz, g_xx_0_y_xxyyzz, g_xx_0_y_xxyzzz, g_xx_0_y_xxzzzz, g_xx_0_y_xyyyyy, g_xx_0_y_xyyyyz, g_xx_0_y_xyyyzz, g_xx_0_y_xyyzzz, g_xx_0_y_xyzzzz, g_xx_0_y_xzzzzz, g_xx_0_y_yyyyyy, g_xx_0_y_yyyyyz, g_xx_0_y_yyyyzz, g_xx_0_y_yyyzzz, g_xx_0_y_yyzzzz, g_xx_0_y_yzzzzz, g_xx_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_y_xxxxxx[k] = -g_xx_0_0_xxxxxx[k] * ab_y + g_xx_0_0_xxxxxxy[k];

                g_xx_0_y_xxxxxy[k] = -g_xx_0_0_xxxxxy[k] * ab_y + g_xx_0_0_xxxxxyy[k];

                g_xx_0_y_xxxxxz[k] = -g_xx_0_0_xxxxxz[k] * ab_y + g_xx_0_0_xxxxxyz[k];

                g_xx_0_y_xxxxyy[k] = -g_xx_0_0_xxxxyy[k] * ab_y + g_xx_0_0_xxxxyyy[k];

                g_xx_0_y_xxxxyz[k] = -g_xx_0_0_xxxxyz[k] * ab_y + g_xx_0_0_xxxxyyz[k];

                g_xx_0_y_xxxxzz[k] = -g_xx_0_0_xxxxzz[k] * ab_y + g_xx_0_0_xxxxyzz[k];

                g_xx_0_y_xxxyyy[k] = -g_xx_0_0_xxxyyy[k] * ab_y + g_xx_0_0_xxxyyyy[k];

                g_xx_0_y_xxxyyz[k] = -g_xx_0_0_xxxyyz[k] * ab_y + g_xx_0_0_xxxyyyz[k];

                g_xx_0_y_xxxyzz[k] = -g_xx_0_0_xxxyzz[k] * ab_y + g_xx_0_0_xxxyyzz[k];

                g_xx_0_y_xxxzzz[k] = -g_xx_0_0_xxxzzz[k] * ab_y + g_xx_0_0_xxxyzzz[k];

                g_xx_0_y_xxyyyy[k] = -g_xx_0_0_xxyyyy[k] * ab_y + g_xx_0_0_xxyyyyy[k];

                g_xx_0_y_xxyyyz[k] = -g_xx_0_0_xxyyyz[k] * ab_y + g_xx_0_0_xxyyyyz[k];

                g_xx_0_y_xxyyzz[k] = -g_xx_0_0_xxyyzz[k] * ab_y + g_xx_0_0_xxyyyzz[k];

                g_xx_0_y_xxyzzz[k] = -g_xx_0_0_xxyzzz[k] * ab_y + g_xx_0_0_xxyyzzz[k];

                g_xx_0_y_xxzzzz[k] = -g_xx_0_0_xxzzzz[k] * ab_y + g_xx_0_0_xxyzzzz[k];

                g_xx_0_y_xyyyyy[k] = -g_xx_0_0_xyyyyy[k] * ab_y + g_xx_0_0_xyyyyyy[k];

                g_xx_0_y_xyyyyz[k] = -g_xx_0_0_xyyyyz[k] * ab_y + g_xx_0_0_xyyyyyz[k];

                g_xx_0_y_xyyyzz[k] = -g_xx_0_0_xyyyzz[k] * ab_y + g_xx_0_0_xyyyyzz[k];

                g_xx_0_y_xyyzzz[k] = -g_xx_0_0_xyyzzz[k] * ab_y + g_xx_0_0_xyyyzzz[k];

                g_xx_0_y_xyzzzz[k] = -g_xx_0_0_xyzzzz[k] * ab_y + g_xx_0_0_xyyzzzz[k];

                g_xx_0_y_xzzzzz[k] = -g_xx_0_0_xzzzzz[k] * ab_y + g_xx_0_0_xyzzzzz[k];

                g_xx_0_y_yyyyyy[k] = -g_xx_0_0_yyyyyy[k] * ab_y + g_xx_0_0_yyyyyyy[k];

                g_xx_0_y_yyyyyz[k] = -g_xx_0_0_yyyyyz[k] * ab_y + g_xx_0_0_yyyyyyz[k];

                g_xx_0_y_yyyyzz[k] = -g_xx_0_0_yyyyzz[k] * ab_y + g_xx_0_0_yyyyyzz[k];

                g_xx_0_y_yyyzzz[k] = -g_xx_0_0_yyyzzz[k] * ab_y + g_xx_0_0_yyyyzzz[k];

                g_xx_0_y_yyzzzz[k] = -g_xx_0_0_yyzzzz[k] * ab_y + g_xx_0_0_yyyzzzz[k];

                g_xx_0_y_yzzzzz[k] = -g_xx_0_0_yzzzzz[k] * ab_y + g_xx_0_0_yyzzzzz[k];

                g_xx_0_y_zzzzzz[k] = -g_xx_0_0_zzzzzz[k] * ab_y + g_xx_0_0_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_xx_0_z_xxxxxx = cbuffer.data(pi_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_z_xxxxxy = cbuffer.data(pi_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_z_xxxxxz = cbuffer.data(pi_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_z_xxxxyy = cbuffer.data(pi_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_z_xxxxyz = cbuffer.data(pi_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_z_xxxxzz = cbuffer.data(pi_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_z_xxxyyy = cbuffer.data(pi_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_z_xxxyyz = cbuffer.data(pi_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_z_xxxyzz = cbuffer.data(pi_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_z_xxxzzz = cbuffer.data(pi_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_z_xxyyyy = cbuffer.data(pi_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_z_xxyyyz = cbuffer.data(pi_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_z_xxyyzz = cbuffer.data(pi_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_z_xxyzzz = cbuffer.data(pi_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_z_xxzzzz = cbuffer.data(pi_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_z_xyyyyy = cbuffer.data(pi_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_z_xyyyyz = cbuffer.data(pi_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_z_xyyyzz = cbuffer.data(pi_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_z_xyyzzz = cbuffer.data(pi_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_z_xyzzzz = cbuffer.data(pi_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_z_xzzzzz = cbuffer.data(pi_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_z_yyyyyy = cbuffer.data(pi_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_z_yyyyyz = cbuffer.data(pi_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_z_yyyyzz = cbuffer.data(pi_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_z_yyyzzz = cbuffer.data(pi_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_z_yyzzzz = cbuffer.data(pi_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_z_yzzzzz = cbuffer.data(pi_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_z_zzzzzz = cbuffer.data(pi_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_0_xxxxxx, g_xx_0_0_xxxxxxz, g_xx_0_0_xxxxxy, g_xx_0_0_xxxxxyz, g_xx_0_0_xxxxxz, g_xx_0_0_xxxxxzz, g_xx_0_0_xxxxyy, g_xx_0_0_xxxxyyz, g_xx_0_0_xxxxyz, g_xx_0_0_xxxxyzz, g_xx_0_0_xxxxzz, g_xx_0_0_xxxxzzz, g_xx_0_0_xxxyyy, g_xx_0_0_xxxyyyz, g_xx_0_0_xxxyyz, g_xx_0_0_xxxyyzz, g_xx_0_0_xxxyzz, g_xx_0_0_xxxyzzz, g_xx_0_0_xxxzzz, g_xx_0_0_xxxzzzz, g_xx_0_0_xxyyyy, g_xx_0_0_xxyyyyz, g_xx_0_0_xxyyyz, g_xx_0_0_xxyyyzz, g_xx_0_0_xxyyzz, g_xx_0_0_xxyyzzz, g_xx_0_0_xxyzzz, g_xx_0_0_xxyzzzz, g_xx_0_0_xxzzzz, g_xx_0_0_xxzzzzz, g_xx_0_0_xyyyyy, g_xx_0_0_xyyyyyz, g_xx_0_0_xyyyyz, g_xx_0_0_xyyyyzz, g_xx_0_0_xyyyzz, g_xx_0_0_xyyyzzz, g_xx_0_0_xyyzzz, g_xx_0_0_xyyzzzz, g_xx_0_0_xyzzzz, g_xx_0_0_xyzzzzz, g_xx_0_0_xzzzzz, g_xx_0_0_xzzzzzz, g_xx_0_0_yyyyyy, g_xx_0_0_yyyyyyz, g_xx_0_0_yyyyyz, g_xx_0_0_yyyyyzz, g_xx_0_0_yyyyzz, g_xx_0_0_yyyyzzz, g_xx_0_0_yyyzzz, g_xx_0_0_yyyzzzz, g_xx_0_0_yyzzzz, g_xx_0_0_yyzzzzz, g_xx_0_0_yzzzzz, g_xx_0_0_yzzzzzz, g_xx_0_0_zzzzzz, g_xx_0_0_zzzzzzz, g_xx_0_z_xxxxxx, g_xx_0_z_xxxxxy, g_xx_0_z_xxxxxz, g_xx_0_z_xxxxyy, g_xx_0_z_xxxxyz, g_xx_0_z_xxxxzz, g_xx_0_z_xxxyyy, g_xx_0_z_xxxyyz, g_xx_0_z_xxxyzz, g_xx_0_z_xxxzzz, g_xx_0_z_xxyyyy, g_xx_0_z_xxyyyz, g_xx_0_z_xxyyzz, g_xx_0_z_xxyzzz, g_xx_0_z_xxzzzz, g_xx_0_z_xyyyyy, g_xx_0_z_xyyyyz, g_xx_0_z_xyyyzz, g_xx_0_z_xyyzzz, g_xx_0_z_xyzzzz, g_xx_0_z_xzzzzz, g_xx_0_z_yyyyyy, g_xx_0_z_yyyyyz, g_xx_0_z_yyyyzz, g_xx_0_z_yyyzzz, g_xx_0_z_yyzzzz, g_xx_0_z_yzzzzz, g_xx_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_z_xxxxxx[k] = -g_xx_0_0_xxxxxx[k] * ab_z + g_xx_0_0_xxxxxxz[k];

                g_xx_0_z_xxxxxy[k] = -g_xx_0_0_xxxxxy[k] * ab_z + g_xx_0_0_xxxxxyz[k];

                g_xx_0_z_xxxxxz[k] = -g_xx_0_0_xxxxxz[k] * ab_z + g_xx_0_0_xxxxxzz[k];

                g_xx_0_z_xxxxyy[k] = -g_xx_0_0_xxxxyy[k] * ab_z + g_xx_0_0_xxxxyyz[k];

                g_xx_0_z_xxxxyz[k] = -g_xx_0_0_xxxxyz[k] * ab_z + g_xx_0_0_xxxxyzz[k];

                g_xx_0_z_xxxxzz[k] = -g_xx_0_0_xxxxzz[k] * ab_z + g_xx_0_0_xxxxzzz[k];

                g_xx_0_z_xxxyyy[k] = -g_xx_0_0_xxxyyy[k] * ab_z + g_xx_0_0_xxxyyyz[k];

                g_xx_0_z_xxxyyz[k] = -g_xx_0_0_xxxyyz[k] * ab_z + g_xx_0_0_xxxyyzz[k];

                g_xx_0_z_xxxyzz[k] = -g_xx_0_0_xxxyzz[k] * ab_z + g_xx_0_0_xxxyzzz[k];

                g_xx_0_z_xxxzzz[k] = -g_xx_0_0_xxxzzz[k] * ab_z + g_xx_0_0_xxxzzzz[k];

                g_xx_0_z_xxyyyy[k] = -g_xx_0_0_xxyyyy[k] * ab_z + g_xx_0_0_xxyyyyz[k];

                g_xx_0_z_xxyyyz[k] = -g_xx_0_0_xxyyyz[k] * ab_z + g_xx_0_0_xxyyyzz[k];

                g_xx_0_z_xxyyzz[k] = -g_xx_0_0_xxyyzz[k] * ab_z + g_xx_0_0_xxyyzzz[k];

                g_xx_0_z_xxyzzz[k] = -g_xx_0_0_xxyzzz[k] * ab_z + g_xx_0_0_xxyzzzz[k];

                g_xx_0_z_xxzzzz[k] = -g_xx_0_0_xxzzzz[k] * ab_z + g_xx_0_0_xxzzzzz[k];

                g_xx_0_z_xyyyyy[k] = -g_xx_0_0_xyyyyy[k] * ab_z + g_xx_0_0_xyyyyyz[k];

                g_xx_0_z_xyyyyz[k] = -g_xx_0_0_xyyyyz[k] * ab_z + g_xx_0_0_xyyyyzz[k];

                g_xx_0_z_xyyyzz[k] = -g_xx_0_0_xyyyzz[k] * ab_z + g_xx_0_0_xyyyzzz[k];

                g_xx_0_z_xyyzzz[k] = -g_xx_0_0_xyyzzz[k] * ab_z + g_xx_0_0_xyyzzzz[k];

                g_xx_0_z_xyzzzz[k] = -g_xx_0_0_xyzzzz[k] * ab_z + g_xx_0_0_xyzzzzz[k];

                g_xx_0_z_xzzzzz[k] = -g_xx_0_0_xzzzzz[k] * ab_z + g_xx_0_0_xzzzzzz[k];

                g_xx_0_z_yyyyyy[k] = -g_xx_0_0_yyyyyy[k] * ab_z + g_xx_0_0_yyyyyyz[k];

                g_xx_0_z_yyyyyz[k] = -g_xx_0_0_yyyyyz[k] * ab_z + g_xx_0_0_yyyyyzz[k];

                g_xx_0_z_yyyyzz[k] = -g_xx_0_0_yyyyzz[k] * ab_z + g_xx_0_0_yyyyzzz[k];

                g_xx_0_z_yyyzzz[k] = -g_xx_0_0_yyyzzz[k] * ab_z + g_xx_0_0_yyyzzzz[k];

                g_xx_0_z_yyzzzz[k] = -g_xx_0_0_yyzzzz[k] * ab_z + g_xx_0_0_yyzzzzz[k];

                g_xx_0_z_yzzzzz[k] = -g_xx_0_0_yzzzzz[k] * ab_z + g_xx_0_0_yzzzzzz[k];

                g_xx_0_z_zzzzzz[k] = -g_xx_0_0_zzzzzz[k] * ab_z + g_xx_0_0_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_xy_0_x_xxxxxx = cbuffer.data(pi_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_x_xxxxxy = cbuffer.data(pi_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_x_xxxxxz = cbuffer.data(pi_geom_20_off + 86 * ccomps * dcomps);

            auto g_xy_0_x_xxxxyy = cbuffer.data(pi_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_x_xxxxyz = cbuffer.data(pi_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_x_xxxxzz = cbuffer.data(pi_geom_20_off + 89 * ccomps * dcomps);

            auto g_xy_0_x_xxxyyy = cbuffer.data(pi_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_x_xxxyyz = cbuffer.data(pi_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_x_xxxyzz = cbuffer.data(pi_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_x_xxxzzz = cbuffer.data(pi_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_x_xxyyyy = cbuffer.data(pi_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_x_xxyyyz = cbuffer.data(pi_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_x_xxyyzz = cbuffer.data(pi_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_x_xxyzzz = cbuffer.data(pi_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_x_xxzzzz = cbuffer.data(pi_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_x_xyyyyy = cbuffer.data(pi_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_x_xyyyyz = cbuffer.data(pi_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_x_xyyyzz = cbuffer.data(pi_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_x_xyyzzz = cbuffer.data(pi_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_x_xyzzzz = cbuffer.data(pi_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_x_xzzzzz = cbuffer.data(pi_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_x_yyyyyy = cbuffer.data(pi_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_x_yyyyyz = cbuffer.data(pi_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_x_yyyyzz = cbuffer.data(pi_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_x_yyyzzz = cbuffer.data(pi_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_x_yyzzzz = cbuffer.data(pi_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_x_yzzzzz = cbuffer.data(pi_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_x_zzzzzz = cbuffer.data(pi_geom_20_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_0_xxxxxx, g_xy_0_0_xxxxxxx, g_xy_0_0_xxxxxxy, g_xy_0_0_xxxxxxz, g_xy_0_0_xxxxxy, g_xy_0_0_xxxxxyy, g_xy_0_0_xxxxxyz, g_xy_0_0_xxxxxz, g_xy_0_0_xxxxxzz, g_xy_0_0_xxxxyy, g_xy_0_0_xxxxyyy, g_xy_0_0_xxxxyyz, g_xy_0_0_xxxxyz, g_xy_0_0_xxxxyzz, g_xy_0_0_xxxxzz, g_xy_0_0_xxxxzzz, g_xy_0_0_xxxyyy, g_xy_0_0_xxxyyyy, g_xy_0_0_xxxyyyz, g_xy_0_0_xxxyyz, g_xy_0_0_xxxyyzz, g_xy_0_0_xxxyzz, g_xy_0_0_xxxyzzz, g_xy_0_0_xxxzzz, g_xy_0_0_xxxzzzz, g_xy_0_0_xxyyyy, g_xy_0_0_xxyyyyy, g_xy_0_0_xxyyyyz, g_xy_0_0_xxyyyz, g_xy_0_0_xxyyyzz, g_xy_0_0_xxyyzz, g_xy_0_0_xxyyzzz, g_xy_0_0_xxyzzz, g_xy_0_0_xxyzzzz, g_xy_0_0_xxzzzz, g_xy_0_0_xxzzzzz, g_xy_0_0_xyyyyy, g_xy_0_0_xyyyyyy, g_xy_0_0_xyyyyyz, g_xy_0_0_xyyyyz, g_xy_0_0_xyyyyzz, g_xy_0_0_xyyyzz, g_xy_0_0_xyyyzzz, g_xy_0_0_xyyzzz, g_xy_0_0_xyyzzzz, g_xy_0_0_xyzzzz, g_xy_0_0_xyzzzzz, g_xy_0_0_xzzzzz, g_xy_0_0_xzzzzzz, g_xy_0_0_yyyyyy, g_xy_0_0_yyyyyz, g_xy_0_0_yyyyzz, g_xy_0_0_yyyzzz, g_xy_0_0_yyzzzz, g_xy_0_0_yzzzzz, g_xy_0_0_zzzzzz, g_xy_0_x_xxxxxx, g_xy_0_x_xxxxxy, g_xy_0_x_xxxxxz, g_xy_0_x_xxxxyy, g_xy_0_x_xxxxyz, g_xy_0_x_xxxxzz, g_xy_0_x_xxxyyy, g_xy_0_x_xxxyyz, g_xy_0_x_xxxyzz, g_xy_0_x_xxxzzz, g_xy_0_x_xxyyyy, g_xy_0_x_xxyyyz, g_xy_0_x_xxyyzz, g_xy_0_x_xxyzzz, g_xy_0_x_xxzzzz, g_xy_0_x_xyyyyy, g_xy_0_x_xyyyyz, g_xy_0_x_xyyyzz, g_xy_0_x_xyyzzz, g_xy_0_x_xyzzzz, g_xy_0_x_xzzzzz, g_xy_0_x_yyyyyy, g_xy_0_x_yyyyyz, g_xy_0_x_yyyyzz, g_xy_0_x_yyyzzz, g_xy_0_x_yyzzzz, g_xy_0_x_yzzzzz, g_xy_0_x_zzzzzz, g_y_0_0_xxxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxxz, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxzz, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_x_xxxxxx[k] = -g_y_0_0_xxxxxx[k] - g_xy_0_0_xxxxxx[k] * ab_x + g_xy_0_0_xxxxxxx[k];

                g_xy_0_x_xxxxxy[k] = -g_y_0_0_xxxxxy[k] - g_xy_0_0_xxxxxy[k] * ab_x + g_xy_0_0_xxxxxxy[k];

                g_xy_0_x_xxxxxz[k] = -g_y_0_0_xxxxxz[k] - g_xy_0_0_xxxxxz[k] * ab_x + g_xy_0_0_xxxxxxz[k];

                g_xy_0_x_xxxxyy[k] = -g_y_0_0_xxxxyy[k] - g_xy_0_0_xxxxyy[k] * ab_x + g_xy_0_0_xxxxxyy[k];

                g_xy_0_x_xxxxyz[k] = -g_y_0_0_xxxxyz[k] - g_xy_0_0_xxxxyz[k] * ab_x + g_xy_0_0_xxxxxyz[k];

                g_xy_0_x_xxxxzz[k] = -g_y_0_0_xxxxzz[k] - g_xy_0_0_xxxxzz[k] * ab_x + g_xy_0_0_xxxxxzz[k];

                g_xy_0_x_xxxyyy[k] = -g_y_0_0_xxxyyy[k] - g_xy_0_0_xxxyyy[k] * ab_x + g_xy_0_0_xxxxyyy[k];

                g_xy_0_x_xxxyyz[k] = -g_y_0_0_xxxyyz[k] - g_xy_0_0_xxxyyz[k] * ab_x + g_xy_0_0_xxxxyyz[k];

                g_xy_0_x_xxxyzz[k] = -g_y_0_0_xxxyzz[k] - g_xy_0_0_xxxyzz[k] * ab_x + g_xy_0_0_xxxxyzz[k];

                g_xy_0_x_xxxzzz[k] = -g_y_0_0_xxxzzz[k] - g_xy_0_0_xxxzzz[k] * ab_x + g_xy_0_0_xxxxzzz[k];

                g_xy_0_x_xxyyyy[k] = -g_y_0_0_xxyyyy[k] - g_xy_0_0_xxyyyy[k] * ab_x + g_xy_0_0_xxxyyyy[k];

                g_xy_0_x_xxyyyz[k] = -g_y_0_0_xxyyyz[k] - g_xy_0_0_xxyyyz[k] * ab_x + g_xy_0_0_xxxyyyz[k];

                g_xy_0_x_xxyyzz[k] = -g_y_0_0_xxyyzz[k] - g_xy_0_0_xxyyzz[k] * ab_x + g_xy_0_0_xxxyyzz[k];

                g_xy_0_x_xxyzzz[k] = -g_y_0_0_xxyzzz[k] - g_xy_0_0_xxyzzz[k] * ab_x + g_xy_0_0_xxxyzzz[k];

                g_xy_0_x_xxzzzz[k] = -g_y_0_0_xxzzzz[k] - g_xy_0_0_xxzzzz[k] * ab_x + g_xy_0_0_xxxzzzz[k];

                g_xy_0_x_xyyyyy[k] = -g_y_0_0_xyyyyy[k] - g_xy_0_0_xyyyyy[k] * ab_x + g_xy_0_0_xxyyyyy[k];

                g_xy_0_x_xyyyyz[k] = -g_y_0_0_xyyyyz[k] - g_xy_0_0_xyyyyz[k] * ab_x + g_xy_0_0_xxyyyyz[k];

                g_xy_0_x_xyyyzz[k] = -g_y_0_0_xyyyzz[k] - g_xy_0_0_xyyyzz[k] * ab_x + g_xy_0_0_xxyyyzz[k];

                g_xy_0_x_xyyzzz[k] = -g_y_0_0_xyyzzz[k] - g_xy_0_0_xyyzzz[k] * ab_x + g_xy_0_0_xxyyzzz[k];

                g_xy_0_x_xyzzzz[k] = -g_y_0_0_xyzzzz[k] - g_xy_0_0_xyzzzz[k] * ab_x + g_xy_0_0_xxyzzzz[k];

                g_xy_0_x_xzzzzz[k] = -g_y_0_0_xzzzzz[k] - g_xy_0_0_xzzzzz[k] * ab_x + g_xy_0_0_xxzzzzz[k];

                g_xy_0_x_yyyyyy[k] = -g_y_0_0_yyyyyy[k] - g_xy_0_0_yyyyyy[k] * ab_x + g_xy_0_0_xyyyyyy[k];

                g_xy_0_x_yyyyyz[k] = -g_y_0_0_yyyyyz[k] - g_xy_0_0_yyyyyz[k] * ab_x + g_xy_0_0_xyyyyyz[k];

                g_xy_0_x_yyyyzz[k] = -g_y_0_0_yyyyzz[k] - g_xy_0_0_yyyyzz[k] * ab_x + g_xy_0_0_xyyyyzz[k];

                g_xy_0_x_yyyzzz[k] = -g_y_0_0_yyyzzz[k] - g_xy_0_0_yyyzzz[k] * ab_x + g_xy_0_0_xyyyzzz[k];

                g_xy_0_x_yyzzzz[k] = -g_y_0_0_yyzzzz[k] - g_xy_0_0_yyzzzz[k] * ab_x + g_xy_0_0_xyyzzzz[k];

                g_xy_0_x_yzzzzz[k] = -g_y_0_0_yzzzzz[k] - g_xy_0_0_yzzzzz[k] * ab_x + g_xy_0_0_xyzzzzz[k];

                g_xy_0_x_zzzzzz[k] = -g_y_0_0_zzzzzz[k] - g_xy_0_0_zzzzzz[k] * ab_x + g_xy_0_0_xzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_xy_0_y_xxxxxx = cbuffer.data(pi_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_y_xxxxxy = cbuffer.data(pi_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_y_xxxxxz = cbuffer.data(pi_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_y_xxxxyy = cbuffer.data(pi_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_y_xxxxyz = cbuffer.data(pi_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_y_xxxxzz = cbuffer.data(pi_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_y_xxxyyy = cbuffer.data(pi_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_y_xxxyyz = cbuffer.data(pi_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_y_xxxyzz = cbuffer.data(pi_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_y_xxxzzz = cbuffer.data(pi_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_y_xxyyyy = cbuffer.data(pi_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_y_xxyyyz = cbuffer.data(pi_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_y_xxyyzz = cbuffer.data(pi_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_y_xxyzzz = cbuffer.data(pi_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_y_xxzzzz = cbuffer.data(pi_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_y_xyyyyy = cbuffer.data(pi_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_y_xyyyyz = cbuffer.data(pi_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_y_xyyyzz = cbuffer.data(pi_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_y_xyyzzz = cbuffer.data(pi_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_y_xyzzzz = cbuffer.data(pi_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_y_xzzzzz = cbuffer.data(pi_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_y_yyyyyy = cbuffer.data(pi_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_y_yyyyyz = cbuffer.data(pi_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_y_yyyyzz = cbuffer.data(pi_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_y_yyyzzz = cbuffer.data(pi_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_y_yyzzzz = cbuffer.data(pi_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_y_yzzzzz = cbuffer.data(pi_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_y_zzzzzz = cbuffer.data(pi_geom_20_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxxz, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxzz, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzzz, g_xy_0_0_xxxxxx, g_xy_0_0_xxxxxxy, g_xy_0_0_xxxxxy, g_xy_0_0_xxxxxyy, g_xy_0_0_xxxxxyz, g_xy_0_0_xxxxxz, g_xy_0_0_xxxxyy, g_xy_0_0_xxxxyyy, g_xy_0_0_xxxxyyz, g_xy_0_0_xxxxyz, g_xy_0_0_xxxxyzz, g_xy_0_0_xxxxzz, g_xy_0_0_xxxyyy, g_xy_0_0_xxxyyyy, g_xy_0_0_xxxyyyz, g_xy_0_0_xxxyyz, g_xy_0_0_xxxyyzz, g_xy_0_0_xxxyzz, g_xy_0_0_xxxyzzz, g_xy_0_0_xxxzzz, g_xy_0_0_xxyyyy, g_xy_0_0_xxyyyyy, g_xy_0_0_xxyyyyz, g_xy_0_0_xxyyyz, g_xy_0_0_xxyyyzz, g_xy_0_0_xxyyzz, g_xy_0_0_xxyyzzz, g_xy_0_0_xxyzzz, g_xy_0_0_xxyzzzz, g_xy_0_0_xxzzzz, g_xy_0_0_xyyyyy, g_xy_0_0_xyyyyyy, g_xy_0_0_xyyyyyz, g_xy_0_0_xyyyyz, g_xy_0_0_xyyyyzz, g_xy_0_0_xyyyzz, g_xy_0_0_xyyyzzz, g_xy_0_0_xyyzzz, g_xy_0_0_xyyzzzz, g_xy_0_0_xyzzzz, g_xy_0_0_xyzzzzz, g_xy_0_0_xzzzzz, g_xy_0_0_yyyyyy, g_xy_0_0_yyyyyyy, g_xy_0_0_yyyyyyz, g_xy_0_0_yyyyyz, g_xy_0_0_yyyyyzz, g_xy_0_0_yyyyzz, g_xy_0_0_yyyyzzz, g_xy_0_0_yyyzzz, g_xy_0_0_yyyzzzz, g_xy_0_0_yyzzzz, g_xy_0_0_yyzzzzz, g_xy_0_0_yzzzzz, g_xy_0_0_yzzzzzz, g_xy_0_0_zzzzzz, g_xy_0_y_xxxxxx, g_xy_0_y_xxxxxy, g_xy_0_y_xxxxxz, g_xy_0_y_xxxxyy, g_xy_0_y_xxxxyz, g_xy_0_y_xxxxzz, g_xy_0_y_xxxyyy, g_xy_0_y_xxxyyz, g_xy_0_y_xxxyzz, g_xy_0_y_xxxzzz, g_xy_0_y_xxyyyy, g_xy_0_y_xxyyyz, g_xy_0_y_xxyyzz, g_xy_0_y_xxyzzz, g_xy_0_y_xxzzzz, g_xy_0_y_xyyyyy, g_xy_0_y_xyyyyz, g_xy_0_y_xyyyzz, g_xy_0_y_xyyzzz, g_xy_0_y_xyzzzz, g_xy_0_y_xzzzzz, g_xy_0_y_yyyyyy, g_xy_0_y_yyyyyz, g_xy_0_y_yyyyzz, g_xy_0_y_yyyzzz, g_xy_0_y_yyzzzz, g_xy_0_y_yzzzzz, g_xy_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_y_xxxxxx[k] = -g_x_0_0_xxxxxx[k] - g_xy_0_0_xxxxxx[k] * ab_y + g_xy_0_0_xxxxxxy[k];

                g_xy_0_y_xxxxxy[k] = -g_x_0_0_xxxxxy[k] - g_xy_0_0_xxxxxy[k] * ab_y + g_xy_0_0_xxxxxyy[k];

                g_xy_0_y_xxxxxz[k] = -g_x_0_0_xxxxxz[k] - g_xy_0_0_xxxxxz[k] * ab_y + g_xy_0_0_xxxxxyz[k];

                g_xy_0_y_xxxxyy[k] = -g_x_0_0_xxxxyy[k] - g_xy_0_0_xxxxyy[k] * ab_y + g_xy_0_0_xxxxyyy[k];

                g_xy_0_y_xxxxyz[k] = -g_x_0_0_xxxxyz[k] - g_xy_0_0_xxxxyz[k] * ab_y + g_xy_0_0_xxxxyyz[k];

                g_xy_0_y_xxxxzz[k] = -g_x_0_0_xxxxzz[k] - g_xy_0_0_xxxxzz[k] * ab_y + g_xy_0_0_xxxxyzz[k];

                g_xy_0_y_xxxyyy[k] = -g_x_0_0_xxxyyy[k] - g_xy_0_0_xxxyyy[k] * ab_y + g_xy_0_0_xxxyyyy[k];

                g_xy_0_y_xxxyyz[k] = -g_x_0_0_xxxyyz[k] - g_xy_0_0_xxxyyz[k] * ab_y + g_xy_0_0_xxxyyyz[k];

                g_xy_0_y_xxxyzz[k] = -g_x_0_0_xxxyzz[k] - g_xy_0_0_xxxyzz[k] * ab_y + g_xy_0_0_xxxyyzz[k];

                g_xy_0_y_xxxzzz[k] = -g_x_0_0_xxxzzz[k] - g_xy_0_0_xxxzzz[k] * ab_y + g_xy_0_0_xxxyzzz[k];

                g_xy_0_y_xxyyyy[k] = -g_x_0_0_xxyyyy[k] - g_xy_0_0_xxyyyy[k] * ab_y + g_xy_0_0_xxyyyyy[k];

                g_xy_0_y_xxyyyz[k] = -g_x_0_0_xxyyyz[k] - g_xy_0_0_xxyyyz[k] * ab_y + g_xy_0_0_xxyyyyz[k];

                g_xy_0_y_xxyyzz[k] = -g_x_0_0_xxyyzz[k] - g_xy_0_0_xxyyzz[k] * ab_y + g_xy_0_0_xxyyyzz[k];

                g_xy_0_y_xxyzzz[k] = -g_x_0_0_xxyzzz[k] - g_xy_0_0_xxyzzz[k] * ab_y + g_xy_0_0_xxyyzzz[k];

                g_xy_0_y_xxzzzz[k] = -g_x_0_0_xxzzzz[k] - g_xy_0_0_xxzzzz[k] * ab_y + g_xy_0_0_xxyzzzz[k];

                g_xy_0_y_xyyyyy[k] = -g_x_0_0_xyyyyy[k] - g_xy_0_0_xyyyyy[k] * ab_y + g_xy_0_0_xyyyyyy[k];

                g_xy_0_y_xyyyyz[k] = -g_x_0_0_xyyyyz[k] - g_xy_0_0_xyyyyz[k] * ab_y + g_xy_0_0_xyyyyyz[k];

                g_xy_0_y_xyyyzz[k] = -g_x_0_0_xyyyzz[k] - g_xy_0_0_xyyyzz[k] * ab_y + g_xy_0_0_xyyyyzz[k];

                g_xy_0_y_xyyzzz[k] = -g_x_0_0_xyyzzz[k] - g_xy_0_0_xyyzzz[k] * ab_y + g_xy_0_0_xyyyzzz[k];

                g_xy_0_y_xyzzzz[k] = -g_x_0_0_xyzzzz[k] - g_xy_0_0_xyzzzz[k] * ab_y + g_xy_0_0_xyyzzzz[k];

                g_xy_0_y_xzzzzz[k] = -g_x_0_0_xzzzzz[k] - g_xy_0_0_xzzzzz[k] * ab_y + g_xy_0_0_xyzzzzz[k];

                g_xy_0_y_yyyyyy[k] = -g_x_0_0_yyyyyy[k] - g_xy_0_0_yyyyyy[k] * ab_y + g_xy_0_0_yyyyyyy[k];

                g_xy_0_y_yyyyyz[k] = -g_x_0_0_yyyyyz[k] - g_xy_0_0_yyyyyz[k] * ab_y + g_xy_0_0_yyyyyyz[k];

                g_xy_0_y_yyyyzz[k] = -g_x_0_0_yyyyzz[k] - g_xy_0_0_yyyyzz[k] * ab_y + g_xy_0_0_yyyyyzz[k];

                g_xy_0_y_yyyzzz[k] = -g_x_0_0_yyyzzz[k] - g_xy_0_0_yyyzzz[k] * ab_y + g_xy_0_0_yyyyzzz[k];

                g_xy_0_y_yyzzzz[k] = -g_x_0_0_yyzzzz[k] - g_xy_0_0_yyzzzz[k] * ab_y + g_xy_0_0_yyyzzzz[k];

                g_xy_0_y_yzzzzz[k] = -g_x_0_0_yzzzzz[k] - g_xy_0_0_yzzzzz[k] * ab_y + g_xy_0_0_yyzzzzz[k];

                g_xy_0_y_zzzzzz[k] = -g_x_0_0_zzzzzz[k] - g_xy_0_0_zzzzzz[k] * ab_y + g_xy_0_0_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_xy_0_z_xxxxxx = cbuffer.data(pi_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_z_xxxxxy = cbuffer.data(pi_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_z_xxxxxz = cbuffer.data(pi_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_z_xxxxyy = cbuffer.data(pi_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_z_xxxxyz = cbuffer.data(pi_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_z_xxxxzz = cbuffer.data(pi_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_z_xxxyyy = cbuffer.data(pi_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_z_xxxyyz = cbuffer.data(pi_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_z_xxxyzz = cbuffer.data(pi_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_z_xxxzzz = cbuffer.data(pi_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_z_xxyyyy = cbuffer.data(pi_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_z_xxyyyz = cbuffer.data(pi_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_z_xxyyzz = cbuffer.data(pi_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_z_xxyzzz = cbuffer.data(pi_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_z_xxzzzz = cbuffer.data(pi_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_z_xyyyyy = cbuffer.data(pi_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_z_xyyyyz = cbuffer.data(pi_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_z_xyyyzz = cbuffer.data(pi_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_z_xyyzzz = cbuffer.data(pi_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_z_xyzzzz = cbuffer.data(pi_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_z_xzzzzz = cbuffer.data(pi_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_z_yyyyyy = cbuffer.data(pi_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_z_yyyyyz = cbuffer.data(pi_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_z_yyyyzz = cbuffer.data(pi_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_z_yyyzzz = cbuffer.data(pi_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_z_yyzzzz = cbuffer.data(pi_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_z_yzzzzz = cbuffer.data(pi_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_z_zzzzzz = cbuffer.data(pi_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_0_xxxxxx, g_xy_0_0_xxxxxxz, g_xy_0_0_xxxxxy, g_xy_0_0_xxxxxyz, g_xy_0_0_xxxxxz, g_xy_0_0_xxxxxzz, g_xy_0_0_xxxxyy, g_xy_0_0_xxxxyyz, g_xy_0_0_xxxxyz, g_xy_0_0_xxxxyzz, g_xy_0_0_xxxxzz, g_xy_0_0_xxxxzzz, g_xy_0_0_xxxyyy, g_xy_0_0_xxxyyyz, g_xy_0_0_xxxyyz, g_xy_0_0_xxxyyzz, g_xy_0_0_xxxyzz, g_xy_0_0_xxxyzzz, g_xy_0_0_xxxzzz, g_xy_0_0_xxxzzzz, g_xy_0_0_xxyyyy, g_xy_0_0_xxyyyyz, g_xy_0_0_xxyyyz, g_xy_0_0_xxyyyzz, g_xy_0_0_xxyyzz, g_xy_0_0_xxyyzzz, g_xy_0_0_xxyzzz, g_xy_0_0_xxyzzzz, g_xy_0_0_xxzzzz, g_xy_0_0_xxzzzzz, g_xy_0_0_xyyyyy, g_xy_0_0_xyyyyyz, g_xy_0_0_xyyyyz, g_xy_0_0_xyyyyzz, g_xy_0_0_xyyyzz, g_xy_0_0_xyyyzzz, g_xy_0_0_xyyzzz, g_xy_0_0_xyyzzzz, g_xy_0_0_xyzzzz, g_xy_0_0_xyzzzzz, g_xy_0_0_xzzzzz, g_xy_0_0_xzzzzzz, g_xy_0_0_yyyyyy, g_xy_0_0_yyyyyyz, g_xy_0_0_yyyyyz, g_xy_0_0_yyyyyzz, g_xy_0_0_yyyyzz, g_xy_0_0_yyyyzzz, g_xy_0_0_yyyzzz, g_xy_0_0_yyyzzzz, g_xy_0_0_yyzzzz, g_xy_0_0_yyzzzzz, g_xy_0_0_yzzzzz, g_xy_0_0_yzzzzzz, g_xy_0_0_zzzzzz, g_xy_0_0_zzzzzzz, g_xy_0_z_xxxxxx, g_xy_0_z_xxxxxy, g_xy_0_z_xxxxxz, g_xy_0_z_xxxxyy, g_xy_0_z_xxxxyz, g_xy_0_z_xxxxzz, g_xy_0_z_xxxyyy, g_xy_0_z_xxxyyz, g_xy_0_z_xxxyzz, g_xy_0_z_xxxzzz, g_xy_0_z_xxyyyy, g_xy_0_z_xxyyyz, g_xy_0_z_xxyyzz, g_xy_0_z_xxyzzz, g_xy_0_z_xxzzzz, g_xy_0_z_xyyyyy, g_xy_0_z_xyyyyz, g_xy_0_z_xyyyzz, g_xy_0_z_xyyzzz, g_xy_0_z_xyzzzz, g_xy_0_z_xzzzzz, g_xy_0_z_yyyyyy, g_xy_0_z_yyyyyz, g_xy_0_z_yyyyzz, g_xy_0_z_yyyzzz, g_xy_0_z_yyzzzz, g_xy_0_z_yzzzzz, g_xy_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_z_xxxxxx[k] = -g_xy_0_0_xxxxxx[k] * ab_z + g_xy_0_0_xxxxxxz[k];

                g_xy_0_z_xxxxxy[k] = -g_xy_0_0_xxxxxy[k] * ab_z + g_xy_0_0_xxxxxyz[k];

                g_xy_0_z_xxxxxz[k] = -g_xy_0_0_xxxxxz[k] * ab_z + g_xy_0_0_xxxxxzz[k];

                g_xy_0_z_xxxxyy[k] = -g_xy_0_0_xxxxyy[k] * ab_z + g_xy_0_0_xxxxyyz[k];

                g_xy_0_z_xxxxyz[k] = -g_xy_0_0_xxxxyz[k] * ab_z + g_xy_0_0_xxxxyzz[k];

                g_xy_0_z_xxxxzz[k] = -g_xy_0_0_xxxxzz[k] * ab_z + g_xy_0_0_xxxxzzz[k];

                g_xy_0_z_xxxyyy[k] = -g_xy_0_0_xxxyyy[k] * ab_z + g_xy_0_0_xxxyyyz[k];

                g_xy_0_z_xxxyyz[k] = -g_xy_0_0_xxxyyz[k] * ab_z + g_xy_0_0_xxxyyzz[k];

                g_xy_0_z_xxxyzz[k] = -g_xy_0_0_xxxyzz[k] * ab_z + g_xy_0_0_xxxyzzz[k];

                g_xy_0_z_xxxzzz[k] = -g_xy_0_0_xxxzzz[k] * ab_z + g_xy_0_0_xxxzzzz[k];

                g_xy_0_z_xxyyyy[k] = -g_xy_0_0_xxyyyy[k] * ab_z + g_xy_0_0_xxyyyyz[k];

                g_xy_0_z_xxyyyz[k] = -g_xy_0_0_xxyyyz[k] * ab_z + g_xy_0_0_xxyyyzz[k];

                g_xy_0_z_xxyyzz[k] = -g_xy_0_0_xxyyzz[k] * ab_z + g_xy_0_0_xxyyzzz[k];

                g_xy_0_z_xxyzzz[k] = -g_xy_0_0_xxyzzz[k] * ab_z + g_xy_0_0_xxyzzzz[k];

                g_xy_0_z_xxzzzz[k] = -g_xy_0_0_xxzzzz[k] * ab_z + g_xy_0_0_xxzzzzz[k];

                g_xy_0_z_xyyyyy[k] = -g_xy_0_0_xyyyyy[k] * ab_z + g_xy_0_0_xyyyyyz[k];

                g_xy_0_z_xyyyyz[k] = -g_xy_0_0_xyyyyz[k] * ab_z + g_xy_0_0_xyyyyzz[k];

                g_xy_0_z_xyyyzz[k] = -g_xy_0_0_xyyyzz[k] * ab_z + g_xy_0_0_xyyyzzz[k];

                g_xy_0_z_xyyzzz[k] = -g_xy_0_0_xyyzzz[k] * ab_z + g_xy_0_0_xyyzzzz[k];

                g_xy_0_z_xyzzzz[k] = -g_xy_0_0_xyzzzz[k] * ab_z + g_xy_0_0_xyzzzzz[k];

                g_xy_0_z_xzzzzz[k] = -g_xy_0_0_xzzzzz[k] * ab_z + g_xy_0_0_xzzzzzz[k];

                g_xy_0_z_yyyyyy[k] = -g_xy_0_0_yyyyyy[k] * ab_z + g_xy_0_0_yyyyyyz[k];

                g_xy_0_z_yyyyyz[k] = -g_xy_0_0_yyyyyz[k] * ab_z + g_xy_0_0_yyyyyzz[k];

                g_xy_0_z_yyyyzz[k] = -g_xy_0_0_yyyyzz[k] * ab_z + g_xy_0_0_yyyyzzz[k];

                g_xy_0_z_yyyzzz[k] = -g_xy_0_0_yyyzzz[k] * ab_z + g_xy_0_0_yyyzzzz[k];

                g_xy_0_z_yyzzzz[k] = -g_xy_0_0_yyzzzz[k] * ab_z + g_xy_0_0_yyzzzzz[k];

                g_xy_0_z_yzzzzz[k] = -g_xy_0_0_yzzzzz[k] * ab_z + g_xy_0_0_yzzzzzz[k];

                g_xy_0_z_zzzzzz[k] = -g_xy_0_0_zzzzzz[k] * ab_z + g_xy_0_0_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_xz_0_x_xxxxxx = cbuffer.data(pi_geom_20_off + 168 * ccomps * dcomps);

            auto g_xz_0_x_xxxxxy = cbuffer.data(pi_geom_20_off + 169 * ccomps * dcomps);

            auto g_xz_0_x_xxxxxz = cbuffer.data(pi_geom_20_off + 170 * ccomps * dcomps);

            auto g_xz_0_x_xxxxyy = cbuffer.data(pi_geom_20_off + 171 * ccomps * dcomps);

            auto g_xz_0_x_xxxxyz = cbuffer.data(pi_geom_20_off + 172 * ccomps * dcomps);

            auto g_xz_0_x_xxxxzz = cbuffer.data(pi_geom_20_off + 173 * ccomps * dcomps);

            auto g_xz_0_x_xxxyyy = cbuffer.data(pi_geom_20_off + 174 * ccomps * dcomps);

            auto g_xz_0_x_xxxyyz = cbuffer.data(pi_geom_20_off + 175 * ccomps * dcomps);

            auto g_xz_0_x_xxxyzz = cbuffer.data(pi_geom_20_off + 176 * ccomps * dcomps);

            auto g_xz_0_x_xxxzzz = cbuffer.data(pi_geom_20_off + 177 * ccomps * dcomps);

            auto g_xz_0_x_xxyyyy = cbuffer.data(pi_geom_20_off + 178 * ccomps * dcomps);

            auto g_xz_0_x_xxyyyz = cbuffer.data(pi_geom_20_off + 179 * ccomps * dcomps);

            auto g_xz_0_x_xxyyzz = cbuffer.data(pi_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_x_xxyzzz = cbuffer.data(pi_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_x_xxzzzz = cbuffer.data(pi_geom_20_off + 182 * ccomps * dcomps);

            auto g_xz_0_x_xyyyyy = cbuffer.data(pi_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_x_xyyyyz = cbuffer.data(pi_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_x_xyyyzz = cbuffer.data(pi_geom_20_off + 185 * ccomps * dcomps);

            auto g_xz_0_x_xyyzzz = cbuffer.data(pi_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_x_xyzzzz = cbuffer.data(pi_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_x_xzzzzz = cbuffer.data(pi_geom_20_off + 188 * ccomps * dcomps);

            auto g_xz_0_x_yyyyyy = cbuffer.data(pi_geom_20_off + 189 * ccomps * dcomps);

            auto g_xz_0_x_yyyyyz = cbuffer.data(pi_geom_20_off + 190 * ccomps * dcomps);

            auto g_xz_0_x_yyyyzz = cbuffer.data(pi_geom_20_off + 191 * ccomps * dcomps);

            auto g_xz_0_x_yyyzzz = cbuffer.data(pi_geom_20_off + 192 * ccomps * dcomps);

            auto g_xz_0_x_yyzzzz = cbuffer.data(pi_geom_20_off + 193 * ccomps * dcomps);

            auto g_xz_0_x_yzzzzz = cbuffer.data(pi_geom_20_off + 194 * ccomps * dcomps);

            auto g_xz_0_x_zzzzzz = cbuffer.data(pi_geom_20_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_0_xxxxxx, g_xz_0_0_xxxxxxx, g_xz_0_0_xxxxxxy, g_xz_0_0_xxxxxxz, g_xz_0_0_xxxxxy, g_xz_0_0_xxxxxyy, g_xz_0_0_xxxxxyz, g_xz_0_0_xxxxxz, g_xz_0_0_xxxxxzz, g_xz_0_0_xxxxyy, g_xz_0_0_xxxxyyy, g_xz_0_0_xxxxyyz, g_xz_0_0_xxxxyz, g_xz_0_0_xxxxyzz, g_xz_0_0_xxxxzz, g_xz_0_0_xxxxzzz, g_xz_0_0_xxxyyy, g_xz_0_0_xxxyyyy, g_xz_0_0_xxxyyyz, g_xz_0_0_xxxyyz, g_xz_0_0_xxxyyzz, g_xz_0_0_xxxyzz, g_xz_0_0_xxxyzzz, g_xz_0_0_xxxzzz, g_xz_0_0_xxxzzzz, g_xz_0_0_xxyyyy, g_xz_0_0_xxyyyyy, g_xz_0_0_xxyyyyz, g_xz_0_0_xxyyyz, g_xz_0_0_xxyyyzz, g_xz_0_0_xxyyzz, g_xz_0_0_xxyyzzz, g_xz_0_0_xxyzzz, g_xz_0_0_xxyzzzz, g_xz_0_0_xxzzzz, g_xz_0_0_xxzzzzz, g_xz_0_0_xyyyyy, g_xz_0_0_xyyyyyy, g_xz_0_0_xyyyyyz, g_xz_0_0_xyyyyz, g_xz_0_0_xyyyyzz, g_xz_0_0_xyyyzz, g_xz_0_0_xyyyzzz, g_xz_0_0_xyyzzz, g_xz_0_0_xyyzzzz, g_xz_0_0_xyzzzz, g_xz_0_0_xyzzzzz, g_xz_0_0_xzzzzz, g_xz_0_0_xzzzzzz, g_xz_0_0_yyyyyy, g_xz_0_0_yyyyyz, g_xz_0_0_yyyyzz, g_xz_0_0_yyyzzz, g_xz_0_0_yyzzzz, g_xz_0_0_yzzzzz, g_xz_0_0_zzzzzz, g_xz_0_x_xxxxxx, g_xz_0_x_xxxxxy, g_xz_0_x_xxxxxz, g_xz_0_x_xxxxyy, g_xz_0_x_xxxxyz, g_xz_0_x_xxxxzz, g_xz_0_x_xxxyyy, g_xz_0_x_xxxyyz, g_xz_0_x_xxxyzz, g_xz_0_x_xxxzzz, g_xz_0_x_xxyyyy, g_xz_0_x_xxyyyz, g_xz_0_x_xxyyzz, g_xz_0_x_xxyzzz, g_xz_0_x_xxzzzz, g_xz_0_x_xyyyyy, g_xz_0_x_xyyyyz, g_xz_0_x_xyyyzz, g_xz_0_x_xyyzzz, g_xz_0_x_xyzzzz, g_xz_0_x_xzzzzz, g_xz_0_x_yyyyyy, g_xz_0_x_yyyyyz, g_xz_0_x_yyyyzz, g_xz_0_x_yyyzzz, g_xz_0_x_yyzzzz, g_xz_0_x_yzzzzz, g_xz_0_x_zzzzzz, g_z_0_0_xxxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxxz, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxzz, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_x_xxxxxx[k] = -g_z_0_0_xxxxxx[k] - g_xz_0_0_xxxxxx[k] * ab_x + g_xz_0_0_xxxxxxx[k];

                g_xz_0_x_xxxxxy[k] = -g_z_0_0_xxxxxy[k] - g_xz_0_0_xxxxxy[k] * ab_x + g_xz_0_0_xxxxxxy[k];

                g_xz_0_x_xxxxxz[k] = -g_z_0_0_xxxxxz[k] - g_xz_0_0_xxxxxz[k] * ab_x + g_xz_0_0_xxxxxxz[k];

                g_xz_0_x_xxxxyy[k] = -g_z_0_0_xxxxyy[k] - g_xz_0_0_xxxxyy[k] * ab_x + g_xz_0_0_xxxxxyy[k];

                g_xz_0_x_xxxxyz[k] = -g_z_0_0_xxxxyz[k] - g_xz_0_0_xxxxyz[k] * ab_x + g_xz_0_0_xxxxxyz[k];

                g_xz_0_x_xxxxzz[k] = -g_z_0_0_xxxxzz[k] - g_xz_0_0_xxxxzz[k] * ab_x + g_xz_0_0_xxxxxzz[k];

                g_xz_0_x_xxxyyy[k] = -g_z_0_0_xxxyyy[k] - g_xz_0_0_xxxyyy[k] * ab_x + g_xz_0_0_xxxxyyy[k];

                g_xz_0_x_xxxyyz[k] = -g_z_0_0_xxxyyz[k] - g_xz_0_0_xxxyyz[k] * ab_x + g_xz_0_0_xxxxyyz[k];

                g_xz_0_x_xxxyzz[k] = -g_z_0_0_xxxyzz[k] - g_xz_0_0_xxxyzz[k] * ab_x + g_xz_0_0_xxxxyzz[k];

                g_xz_0_x_xxxzzz[k] = -g_z_0_0_xxxzzz[k] - g_xz_0_0_xxxzzz[k] * ab_x + g_xz_0_0_xxxxzzz[k];

                g_xz_0_x_xxyyyy[k] = -g_z_0_0_xxyyyy[k] - g_xz_0_0_xxyyyy[k] * ab_x + g_xz_0_0_xxxyyyy[k];

                g_xz_0_x_xxyyyz[k] = -g_z_0_0_xxyyyz[k] - g_xz_0_0_xxyyyz[k] * ab_x + g_xz_0_0_xxxyyyz[k];

                g_xz_0_x_xxyyzz[k] = -g_z_0_0_xxyyzz[k] - g_xz_0_0_xxyyzz[k] * ab_x + g_xz_0_0_xxxyyzz[k];

                g_xz_0_x_xxyzzz[k] = -g_z_0_0_xxyzzz[k] - g_xz_0_0_xxyzzz[k] * ab_x + g_xz_0_0_xxxyzzz[k];

                g_xz_0_x_xxzzzz[k] = -g_z_0_0_xxzzzz[k] - g_xz_0_0_xxzzzz[k] * ab_x + g_xz_0_0_xxxzzzz[k];

                g_xz_0_x_xyyyyy[k] = -g_z_0_0_xyyyyy[k] - g_xz_0_0_xyyyyy[k] * ab_x + g_xz_0_0_xxyyyyy[k];

                g_xz_0_x_xyyyyz[k] = -g_z_0_0_xyyyyz[k] - g_xz_0_0_xyyyyz[k] * ab_x + g_xz_0_0_xxyyyyz[k];

                g_xz_0_x_xyyyzz[k] = -g_z_0_0_xyyyzz[k] - g_xz_0_0_xyyyzz[k] * ab_x + g_xz_0_0_xxyyyzz[k];

                g_xz_0_x_xyyzzz[k] = -g_z_0_0_xyyzzz[k] - g_xz_0_0_xyyzzz[k] * ab_x + g_xz_0_0_xxyyzzz[k];

                g_xz_0_x_xyzzzz[k] = -g_z_0_0_xyzzzz[k] - g_xz_0_0_xyzzzz[k] * ab_x + g_xz_0_0_xxyzzzz[k];

                g_xz_0_x_xzzzzz[k] = -g_z_0_0_xzzzzz[k] - g_xz_0_0_xzzzzz[k] * ab_x + g_xz_0_0_xxzzzzz[k];

                g_xz_0_x_yyyyyy[k] = -g_z_0_0_yyyyyy[k] - g_xz_0_0_yyyyyy[k] * ab_x + g_xz_0_0_xyyyyyy[k];

                g_xz_0_x_yyyyyz[k] = -g_z_0_0_yyyyyz[k] - g_xz_0_0_yyyyyz[k] * ab_x + g_xz_0_0_xyyyyyz[k];

                g_xz_0_x_yyyyzz[k] = -g_z_0_0_yyyyzz[k] - g_xz_0_0_yyyyzz[k] * ab_x + g_xz_0_0_xyyyyzz[k];

                g_xz_0_x_yyyzzz[k] = -g_z_0_0_yyyzzz[k] - g_xz_0_0_yyyzzz[k] * ab_x + g_xz_0_0_xyyyzzz[k];

                g_xz_0_x_yyzzzz[k] = -g_z_0_0_yyzzzz[k] - g_xz_0_0_yyzzzz[k] * ab_x + g_xz_0_0_xyyzzzz[k];

                g_xz_0_x_yzzzzz[k] = -g_z_0_0_yzzzzz[k] - g_xz_0_0_yzzzzz[k] * ab_x + g_xz_0_0_xyzzzzz[k];

                g_xz_0_x_zzzzzz[k] = -g_z_0_0_zzzzzz[k] - g_xz_0_0_zzzzzz[k] * ab_x + g_xz_0_0_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_xz_0_y_xxxxxx = cbuffer.data(pi_geom_20_off + 196 * ccomps * dcomps);

            auto g_xz_0_y_xxxxxy = cbuffer.data(pi_geom_20_off + 197 * ccomps * dcomps);

            auto g_xz_0_y_xxxxxz = cbuffer.data(pi_geom_20_off + 198 * ccomps * dcomps);

            auto g_xz_0_y_xxxxyy = cbuffer.data(pi_geom_20_off + 199 * ccomps * dcomps);

            auto g_xz_0_y_xxxxyz = cbuffer.data(pi_geom_20_off + 200 * ccomps * dcomps);

            auto g_xz_0_y_xxxxzz = cbuffer.data(pi_geom_20_off + 201 * ccomps * dcomps);

            auto g_xz_0_y_xxxyyy = cbuffer.data(pi_geom_20_off + 202 * ccomps * dcomps);

            auto g_xz_0_y_xxxyyz = cbuffer.data(pi_geom_20_off + 203 * ccomps * dcomps);

            auto g_xz_0_y_xxxyzz = cbuffer.data(pi_geom_20_off + 204 * ccomps * dcomps);

            auto g_xz_0_y_xxxzzz = cbuffer.data(pi_geom_20_off + 205 * ccomps * dcomps);

            auto g_xz_0_y_xxyyyy = cbuffer.data(pi_geom_20_off + 206 * ccomps * dcomps);

            auto g_xz_0_y_xxyyyz = cbuffer.data(pi_geom_20_off + 207 * ccomps * dcomps);

            auto g_xz_0_y_xxyyzz = cbuffer.data(pi_geom_20_off + 208 * ccomps * dcomps);

            auto g_xz_0_y_xxyzzz = cbuffer.data(pi_geom_20_off + 209 * ccomps * dcomps);

            auto g_xz_0_y_xxzzzz = cbuffer.data(pi_geom_20_off + 210 * ccomps * dcomps);

            auto g_xz_0_y_xyyyyy = cbuffer.data(pi_geom_20_off + 211 * ccomps * dcomps);

            auto g_xz_0_y_xyyyyz = cbuffer.data(pi_geom_20_off + 212 * ccomps * dcomps);

            auto g_xz_0_y_xyyyzz = cbuffer.data(pi_geom_20_off + 213 * ccomps * dcomps);

            auto g_xz_0_y_xyyzzz = cbuffer.data(pi_geom_20_off + 214 * ccomps * dcomps);

            auto g_xz_0_y_xyzzzz = cbuffer.data(pi_geom_20_off + 215 * ccomps * dcomps);

            auto g_xz_0_y_xzzzzz = cbuffer.data(pi_geom_20_off + 216 * ccomps * dcomps);

            auto g_xz_0_y_yyyyyy = cbuffer.data(pi_geom_20_off + 217 * ccomps * dcomps);

            auto g_xz_0_y_yyyyyz = cbuffer.data(pi_geom_20_off + 218 * ccomps * dcomps);

            auto g_xz_0_y_yyyyzz = cbuffer.data(pi_geom_20_off + 219 * ccomps * dcomps);

            auto g_xz_0_y_yyyzzz = cbuffer.data(pi_geom_20_off + 220 * ccomps * dcomps);

            auto g_xz_0_y_yyzzzz = cbuffer.data(pi_geom_20_off + 221 * ccomps * dcomps);

            auto g_xz_0_y_yzzzzz = cbuffer.data(pi_geom_20_off + 222 * ccomps * dcomps);

            auto g_xz_0_y_zzzzzz = cbuffer.data(pi_geom_20_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_0_xxxxxx, g_xz_0_0_xxxxxxy, g_xz_0_0_xxxxxy, g_xz_0_0_xxxxxyy, g_xz_0_0_xxxxxyz, g_xz_0_0_xxxxxz, g_xz_0_0_xxxxyy, g_xz_0_0_xxxxyyy, g_xz_0_0_xxxxyyz, g_xz_0_0_xxxxyz, g_xz_0_0_xxxxyzz, g_xz_0_0_xxxxzz, g_xz_0_0_xxxyyy, g_xz_0_0_xxxyyyy, g_xz_0_0_xxxyyyz, g_xz_0_0_xxxyyz, g_xz_0_0_xxxyyzz, g_xz_0_0_xxxyzz, g_xz_0_0_xxxyzzz, g_xz_0_0_xxxzzz, g_xz_0_0_xxyyyy, g_xz_0_0_xxyyyyy, g_xz_0_0_xxyyyyz, g_xz_0_0_xxyyyz, g_xz_0_0_xxyyyzz, g_xz_0_0_xxyyzz, g_xz_0_0_xxyyzzz, g_xz_0_0_xxyzzz, g_xz_0_0_xxyzzzz, g_xz_0_0_xxzzzz, g_xz_0_0_xyyyyy, g_xz_0_0_xyyyyyy, g_xz_0_0_xyyyyyz, g_xz_0_0_xyyyyz, g_xz_0_0_xyyyyzz, g_xz_0_0_xyyyzz, g_xz_0_0_xyyyzzz, g_xz_0_0_xyyzzz, g_xz_0_0_xyyzzzz, g_xz_0_0_xyzzzz, g_xz_0_0_xyzzzzz, g_xz_0_0_xzzzzz, g_xz_0_0_yyyyyy, g_xz_0_0_yyyyyyy, g_xz_0_0_yyyyyyz, g_xz_0_0_yyyyyz, g_xz_0_0_yyyyyzz, g_xz_0_0_yyyyzz, g_xz_0_0_yyyyzzz, g_xz_0_0_yyyzzz, g_xz_0_0_yyyzzzz, g_xz_0_0_yyzzzz, g_xz_0_0_yyzzzzz, g_xz_0_0_yzzzzz, g_xz_0_0_yzzzzzz, g_xz_0_0_zzzzzz, g_xz_0_y_xxxxxx, g_xz_0_y_xxxxxy, g_xz_0_y_xxxxxz, g_xz_0_y_xxxxyy, g_xz_0_y_xxxxyz, g_xz_0_y_xxxxzz, g_xz_0_y_xxxyyy, g_xz_0_y_xxxyyz, g_xz_0_y_xxxyzz, g_xz_0_y_xxxzzz, g_xz_0_y_xxyyyy, g_xz_0_y_xxyyyz, g_xz_0_y_xxyyzz, g_xz_0_y_xxyzzz, g_xz_0_y_xxzzzz, g_xz_0_y_xyyyyy, g_xz_0_y_xyyyyz, g_xz_0_y_xyyyzz, g_xz_0_y_xyyzzz, g_xz_0_y_xyzzzz, g_xz_0_y_xzzzzz, g_xz_0_y_yyyyyy, g_xz_0_y_yyyyyz, g_xz_0_y_yyyyzz, g_xz_0_y_yyyzzz, g_xz_0_y_yyzzzz, g_xz_0_y_yzzzzz, g_xz_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_y_xxxxxx[k] = -g_xz_0_0_xxxxxx[k] * ab_y + g_xz_0_0_xxxxxxy[k];

                g_xz_0_y_xxxxxy[k] = -g_xz_0_0_xxxxxy[k] * ab_y + g_xz_0_0_xxxxxyy[k];

                g_xz_0_y_xxxxxz[k] = -g_xz_0_0_xxxxxz[k] * ab_y + g_xz_0_0_xxxxxyz[k];

                g_xz_0_y_xxxxyy[k] = -g_xz_0_0_xxxxyy[k] * ab_y + g_xz_0_0_xxxxyyy[k];

                g_xz_0_y_xxxxyz[k] = -g_xz_0_0_xxxxyz[k] * ab_y + g_xz_0_0_xxxxyyz[k];

                g_xz_0_y_xxxxzz[k] = -g_xz_0_0_xxxxzz[k] * ab_y + g_xz_0_0_xxxxyzz[k];

                g_xz_0_y_xxxyyy[k] = -g_xz_0_0_xxxyyy[k] * ab_y + g_xz_0_0_xxxyyyy[k];

                g_xz_0_y_xxxyyz[k] = -g_xz_0_0_xxxyyz[k] * ab_y + g_xz_0_0_xxxyyyz[k];

                g_xz_0_y_xxxyzz[k] = -g_xz_0_0_xxxyzz[k] * ab_y + g_xz_0_0_xxxyyzz[k];

                g_xz_0_y_xxxzzz[k] = -g_xz_0_0_xxxzzz[k] * ab_y + g_xz_0_0_xxxyzzz[k];

                g_xz_0_y_xxyyyy[k] = -g_xz_0_0_xxyyyy[k] * ab_y + g_xz_0_0_xxyyyyy[k];

                g_xz_0_y_xxyyyz[k] = -g_xz_0_0_xxyyyz[k] * ab_y + g_xz_0_0_xxyyyyz[k];

                g_xz_0_y_xxyyzz[k] = -g_xz_0_0_xxyyzz[k] * ab_y + g_xz_0_0_xxyyyzz[k];

                g_xz_0_y_xxyzzz[k] = -g_xz_0_0_xxyzzz[k] * ab_y + g_xz_0_0_xxyyzzz[k];

                g_xz_0_y_xxzzzz[k] = -g_xz_0_0_xxzzzz[k] * ab_y + g_xz_0_0_xxyzzzz[k];

                g_xz_0_y_xyyyyy[k] = -g_xz_0_0_xyyyyy[k] * ab_y + g_xz_0_0_xyyyyyy[k];

                g_xz_0_y_xyyyyz[k] = -g_xz_0_0_xyyyyz[k] * ab_y + g_xz_0_0_xyyyyyz[k];

                g_xz_0_y_xyyyzz[k] = -g_xz_0_0_xyyyzz[k] * ab_y + g_xz_0_0_xyyyyzz[k];

                g_xz_0_y_xyyzzz[k] = -g_xz_0_0_xyyzzz[k] * ab_y + g_xz_0_0_xyyyzzz[k];

                g_xz_0_y_xyzzzz[k] = -g_xz_0_0_xyzzzz[k] * ab_y + g_xz_0_0_xyyzzzz[k];

                g_xz_0_y_xzzzzz[k] = -g_xz_0_0_xzzzzz[k] * ab_y + g_xz_0_0_xyzzzzz[k];

                g_xz_0_y_yyyyyy[k] = -g_xz_0_0_yyyyyy[k] * ab_y + g_xz_0_0_yyyyyyy[k];

                g_xz_0_y_yyyyyz[k] = -g_xz_0_0_yyyyyz[k] * ab_y + g_xz_0_0_yyyyyyz[k];

                g_xz_0_y_yyyyzz[k] = -g_xz_0_0_yyyyzz[k] * ab_y + g_xz_0_0_yyyyyzz[k];

                g_xz_0_y_yyyzzz[k] = -g_xz_0_0_yyyzzz[k] * ab_y + g_xz_0_0_yyyyzzz[k];

                g_xz_0_y_yyzzzz[k] = -g_xz_0_0_yyzzzz[k] * ab_y + g_xz_0_0_yyyzzzz[k];

                g_xz_0_y_yzzzzz[k] = -g_xz_0_0_yzzzzz[k] * ab_y + g_xz_0_0_yyzzzzz[k];

                g_xz_0_y_zzzzzz[k] = -g_xz_0_0_zzzzzz[k] * ab_y + g_xz_0_0_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_xz_0_z_xxxxxx = cbuffer.data(pi_geom_20_off + 224 * ccomps * dcomps);

            auto g_xz_0_z_xxxxxy = cbuffer.data(pi_geom_20_off + 225 * ccomps * dcomps);

            auto g_xz_0_z_xxxxxz = cbuffer.data(pi_geom_20_off + 226 * ccomps * dcomps);

            auto g_xz_0_z_xxxxyy = cbuffer.data(pi_geom_20_off + 227 * ccomps * dcomps);

            auto g_xz_0_z_xxxxyz = cbuffer.data(pi_geom_20_off + 228 * ccomps * dcomps);

            auto g_xz_0_z_xxxxzz = cbuffer.data(pi_geom_20_off + 229 * ccomps * dcomps);

            auto g_xz_0_z_xxxyyy = cbuffer.data(pi_geom_20_off + 230 * ccomps * dcomps);

            auto g_xz_0_z_xxxyyz = cbuffer.data(pi_geom_20_off + 231 * ccomps * dcomps);

            auto g_xz_0_z_xxxyzz = cbuffer.data(pi_geom_20_off + 232 * ccomps * dcomps);

            auto g_xz_0_z_xxxzzz = cbuffer.data(pi_geom_20_off + 233 * ccomps * dcomps);

            auto g_xz_0_z_xxyyyy = cbuffer.data(pi_geom_20_off + 234 * ccomps * dcomps);

            auto g_xz_0_z_xxyyyz = cbuffer.data(pi_geom_20_off + 235 * ccomps * dcomps);

            auto g_xz_0_z_xxyyzz = cbuffer.data(pi_geom_20_off + 236 * ccomps * dcomps);

            auto g_xz_0_z_xxyzzz = cbuffer.data(pi_geom_20_off + 237 * ccomps * dcomps);

            auto g_xz_0_z_xxzzzz = cbuffer.data(pi_geom_20_off + 238 * ccomps * dcomps);

            auto g_xz_0_z_xyyyyy = cbuffer.data(pi_geom_20_off + 239 * ccomps * dcomps);

            auto g_xz_0_z_xyyyyz = cbuffer.data(pi_geom_20_off + 240 * ccomps * dcomps);

            auto g_xz_0_z_xyyyzz = cbuffer.data(pi_geom_20_off + 241 * ccomps * dcomps);

            auto g_xz_0_z_xyyzzz = cbuffer.data(pi_geom_20_off + 242 * ccomps * dcomps);

            auto g_xz_0_z_xyzzzz = cbuffer.data(pi_geom_20_off + 243 * ccomps * dcomps);

            auto g_xz_0_z_xzzzzz = cbuffer.data(pi_geom_20_off + 244 * ccomps * dcomps);

            auto g_xz_0_z_yyyyyy = cbuffer.data(pi_geom_20_off + 245 * ccomps * dcomps);

            auto g_xz_0_z_yyyyyz = cbuffer.data(pi_geom_20_off + 246 * ccomps * dcomps);

            auto g_xz_0_z_yyyyzz = cbuffer.data(pi_geom_20_off + 247 * ccomps * dcomps);

            auto g_xz_0_z_yyyzzz = cbuffer.data(pi_geom_20_off + 248 * ccomps * dcomps);

            auto g_xz_0_z_yyzzzz = cbuffer.data(pi_geom_20_off + 249 * ccomps * dcomps);

            auto g_xz_0_z_yzzzzz = cbuffer.data(pi_geom_20_off + 250 * ccomps * dcomps);

            auto g_xz_0_z_zzzzzz = cbuffer.data(pi_geom_20_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxxz, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxzz, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzzz, g_xz_0_0_xxxxxx, g_xz_0_0_xxxxxxz, g_xz_0_0_xxxxxy, g_xz_0_0_xxxxxyz, g_xz_0_0_xxxxxz, g_xz_0_0_xxxxxzz, g_xz_0_0_xxxxyy, g_xz_0_0_xxxxyyz, g_xz_0_0_xxxxyz, g_xz_0_0_xxxxyzz, g_xz_0_0_xxxxzz, g_xz_0_0_xxxxzzz, g_xz_0_0_xxxyyy, g_xz_0_0_xxxyyyz, g_xz_0_0_xxxyyz, g_xz_0_0_xxxyyzz, g_xz_0_0_xxxyzz, g_xz_0_0_xxxyzzz, g_xz_0_0_xxxzzz, g_xz_0_0_xxxzzzz, g_xz_0_0_xxyyyy, g_xz_0_0_xxyyyyz, g_xz_0_0_xxyyyz, g_xz_0_0_xxyyyzz, g_xz_0_0_xxyyzz, g_xz_0_0_xxyyzzz, g_xz_0_0_xxyzzz, g_xz_0_0_xxyzzzz, g_xz_0_0_xxzzzz, g_xz_0_0_xxzzzzz, g_xz_0_0_xyyyyy, g_xz_0_0_xyyyyyz, g_xz_0_0_xyyyyz, g_xz_0_0_xyyyyzz, g_xz_0_0_xyyyzz, g_xz_0_0_xyyyzzz, g_xz_0_0_xyyzzz, g_xz_0_0_xyyzzzz, g_xz_0_0_xyzzzz, g_xz_0_0_xyzzzzz, g_xz_0_0_xzzzzz, g_xz_0_0_xzzzzzz, g_xz_0_0_yyyyyy, g_xz_0_0_yyyyyyz, g_xz_0_0_yyyyyz, g_xz_0_0_yyyyyzz, g_xz_0_0_yyyyzz, g_xz_0_0_yyyyzzz, g_xz_0_0_yyyzzz, g_xz_0_0_yyyzzzz, g_xz_0_0_yyzzzz, g_xz_0_0_yyzzzzz, g_xz_0_0_yzzzzz, g_xz_0_0_yzzzzzz, g_xz_0_0_zzzzzz, g_xz_0_0_zzzzzzz, g_xz_0_z_xxxxxx, g_xz_0_z_xxxxxy, g_xz_0_z_xxxxxz, g_xz_0_z_xxxxyy, g_xz_0_z_xxxxyz, g_xz_0_z_xxxxzz, g_xz_0_z_xxxyyy, g_xz_0_z_xxxyyz, g_xz_0_z_xxxyzz, g_xz_0_z_xxxzzz, g_xz_0_z_xxyyyy, g_xz_0_z_xxyyyz, g_xz_0_z_xxyyzz, g_xz_0_z_xxyzzz, g_xz_0_z_xxzzzz, g_xz_0_z_xyyyyy, g_xz_0_z_xyyyyz, g_xz_0_z_xyyyzz, g_xz_0_z_xyyzzz, g_xz_0_z_xyzzzz, g_xz_0_z_xzzzzz, g_xz_0_z_yyyyyy, g_xz_0_z_yyyyyz, g_xz_0_z_yyyyzz, g_xz_0_z_yyyzzz, g_xz_0_z_yyzzzz, g_xz_0_z_yzzzzz, g_xz_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_z_xxxxxx[k] = -g_x_0_0_xxxxxx[k] - g_xz_0_0_xxxxxx[k] * ab_z + g_xz_0_0_xxxxxxz[k];

                g_xz_0_z_xxxxxy[k] = -g_x_0_0_xxxxxy[k] - g_xz_0_0_xxxxxy[k] * ab_z + g_xz_0_0_xxxxxyz[k];

                g_xz_0_z_xxxxxz[k] = -g_x_0_0_xxxxxz[k] - g_xz_0_0_xxxxxz[k] * ab_z + g_xz_0_0_xxxxxzz[k];

                g_xz_0_z_xxxxyy[k] = -g_x_0_0_xxxxyy[k] - g_xz_0_0_xxxxyy[k] * ab_z + g_xz_0_0_xxxxyyz[k];

                g_xz_0_z_xxxxyz[k] = -g_x_0_0_xxxxyz[k] - g_xz_0_0_xxxxyz[k] * ab_z + g_xz_0_0_xxxxyzz[k];

                g_xz_0_z_xxxxzz[k] = -g_x_0_0_xxxxzz[k] - g_xz_0_0_xxxxzz[k] * ab_z + g_xz_0_0_xxxxzzz[k];

                g_xz_0_z_xxxyyy[k] = -g_x_0_0_xxxyyy[k] - g_xz_0_0_xxxyyy[k] * ab_z + g_xz_0_0_xxxyyyz[k];

                g_xz_0_z_xxxyyz[k] = -g_x_0_0_xxxyyz[k] - g_xz_0_0_xxxyyz[k] * ab_z + g_xz_0_0_xxxyyzz[k];

                g_xz_0_z_xxxyzz[k] = -g_x_0_0_xxxyzz[k] - g_xz_0_0_xxxyzz[k] * ab_z + g_xz_0_0_xxxyzzz[k];

                g_xz_0_z_xxxzzz[k] = -g_x_0_0_xxxzzz[k] - g_xz_0_0_xxxzzz[k] * ab_z + g_xz_0_0_xxxzzzz[k];

                g_xz_0_z_xxyyyy[k] = -g_x_0_0_xxyyyy[k] - g_xz_0_0_xxyyyy[k] * ab_z + g_xz_0_0_xxyyyyz[k];

                g_xz_0_z_xxyyyz[k] = -g_x_0_0_xxyyyz[k] - g_xz_0_0_xxyyyz[k] * ab_z + g_xz_0_0_xxyyyzz[k];

                g_xz_0_z_xxyyzz[k] = -g_x_0_0_xxyyzz[k] - g_xz_0_0_xxyyzz[k] * ab_z + g_xz_0_0_xxyyzzz[k];

                g_xz_0_z_xxyzzz[k] = -g_x_0_0_xxyzzz[k] - g_xz_0_0_xxyzzz[k] * ab_z + g_xz_0_0_xxyzzzz[k];

                g_xz_0_z_xxzzzz[k] = -g_x_0_0_xxzzzz[k] - g_xz_0_0_xxzzzz[k] * ab_z + g_xz_0_0_xxzzzzz[k];

                g_xz_0_z_xyyyyy[k] = -g_x_0_0_xyyyyy[k] - g_xz_0_0_xyyyyy[k] * ab_z + g_xz_0_0_xyyyyyz[k];

                g_xz_0_z_xyyyyz[k] = -g_x_0_0_xyyyyz[k] - g_xz_0_0_xyyyyz[k] * ab_z + g_xz_0_0_xyyyyzz[k];

                g_xz_0_z_xyyyzz[k] = -g_x_0_0_xyyyzz[k] - g_xz_0_0_xyyyzz[k] * ab_z + g_xz_0_0_xyyyzzz[k];

                g_xz_0_z_xyyzzz[k] = -g_x_0_0_xyyzzz[k] - g_xz_0_0_xyyzzz[k] * ab_z + g_xz_0_0_xyyzzzz[k];

                g_xz_0_z_xyzzzz[k] = -g_x_0_0_xyzzzz[k] - g_xz_0_0_xyzzzz[k] * ab_z + g_xz_0_0_xyzzzzz[k];

                g_xz_0_z_xzzzzz[k] = -g_x_0_0_xzzzzz[k] - g_xz_0_0_xzzzzz[k] * ab_z + g_xz_0_0_xzzzzzz[k];

                g_xz_0_z_yyyyyy[k] = -g_x_0_0_yyyyyy[k] - g_xz_0_0_yyyyyy[k] * ab_z + g_xz_0_0_yyyyyyz[k];

                g_xz_0_z_yyyyyz[k] = -g_x_0_0_yyyyyz[k] - g_xz_0_0_yyyyyz[k] * ab_z + g_xz_0_0_yyyyyzz[k];

                g_xz_0_z_yyyyzz[k] = -g_x_0_0_yyyyzz[k] - g_xz_0_0_yyyyzz[k] * ab_z + g_xz_0_0_yyyyzzz[k];

                g_xz_0_z_yyyzzz[k] = -g_x_0_0_yyyzzz[k] - g_xz_0_0_yyyzzz[k] * ab_z + g_xz_0_0_yyyzzzz[k];

                g_xz_0_z_yyzzzz[k] = -g_x_0_0_yyzzzz[k] - g_xz_0_0_yyzzzz[k] * ab_z + g_xz_0_0_yyzzzzz[k];

                g_xz_0_z_yzzzzz[k] = -g_x_0_0_yzzzzz[k] - g_xz_0_0_yzzzzz[k] * ab_z + g_xz_0_0_yzzzzzz[k];

                g_xz_0_z_zzzzzz[k] = -g_x_0_0_zzzzzz[k] - g_xz_0_0_zzzzzz[k] * ab_z + g_xz_0_0_zzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_yy_0_x_xxxxxx = cbuffer.data(pi_geom_20_off + 252 * ccomps * dcomps);

            auto g_yy_0_x_xxxxxy = cbuffer.data(pi_geom_20_off + 253 * ccomps * dcomps);

            auto g_yy_0_x_xxxxxz = cbuffer.data(pi_geom_20_off + 254 * ccomps * dcomps);

            auto g_yy_0_x_xxxxyy = cbuffer.data(pi_geom_20_off + 255 * ccomps * dcomps);

            auto g_yy_0_x_xxxxyz = cbuffer.data(pi_geom_20_off + 256 * ccomps * dcomps);

            auto g_yy_0_x_xxxxzz = cbuffer.data(pi_geom_20_off + 257 * ccomps * dcomps);

            auto g_yy_0_x_xxxyyy = cbuffer.data(pi_geom_20_off + 258 * ccomps * dcomps);

            auto g_yy_0_x_xxxyyz = cbuffer.data(pi_geom_20_off + 259 * ccomps * dcomps);

            auto g_yy_0_x_xxxyzz = cbuffer.data(pi_geom_20_off + 260 * ccomps * dcomps);

            auto g_yy_0_x_xxxzzz = cbuffer.data(pi_geom_20_off + 261 * ccomps * dcomps);

            auto g_yy_0_x_xxyyyy = cbuffer.data(pi_geom_20_off + 262 * ccomps * dcomps);

            auto g_yy_0_x_xxyyyz = cbuffer.data(pi_geom_20_off + 263 * ccomps * dcomps);

            auto g_yy_0_x_xxyyzz = cbuffer.data(pi_geom_20_off + 264 * ccomps * dcomps);

            auto g_yy_0_x_xxyzzz = cbuffer.data(pi_geom_20_off + 265 * ccomps * dcomps);

            auto g_yy_0_x_xxzzzz = cbuffer.data(pi_geom_20_off + 266 * ccomps * dcomps);

            auto g_yy_0_x_xyyyyy = cbuffer.data(pi_geom_20_off + 267 * ccomps * dcomps);

            auto g_yy_0_x_xyyyyz = cbuffer.data(pi_geom_20_off + 268 * ccomps * dcomps);

            auto g_yy_0_x_xyyyzz = cbuffer.data(pi_geom_20_off + 269 * ccomps * dcomps);

            auto g_yy_0_x_xyyzzz = cbuffer.data(pi_geom_20_off + 270 * ccomps * dcomps);

            auto g_yy_0_x_xyzzzz = cbuffer.data(pi_geom_20_off + 271 * ccomps * dcomps);

            auto g_yy_0_x_xzzzzz = cbuffer.data(pi_geom_20_off + 272 * ccomps * dcomps);

            auto g_yy_0_x_yyyyyy = cbuffer.data(pi_geom_20_off + 273 * ccomps * dcomps);

            auto g_yy_0_x_yyyyyz = cbuffer.data(pi_geom_20_off + 274 * ccomps * dcomps);

            auto g_yy_0_x_yyyyzz = cbuffer.data(pi_geom_20_off + 275 * ccomps * dcomps);

            auto g_yy_0_x_yyyzzz = cbuffer.data(pi_geom_20_off + 276 * ccomps * dcomps);

            auto g_yy_0_x_yyzzzz = cbuffer.data(pi_geom_20_off + 277 * ccomps * dcomps);

            auto g_yy_0_x_yzzzzz = cbuffer.data(pi_geom_20_off + 278 * ccomps * dcomps);

            auto g_yy_0_x_zzzzzz = cbuffer.data(pi_geom_20_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_0_xxxxxx, g_yy_0_0_xxxxxxx, g_yy_0_0_xxxxxxy, g_yy_0_0_xxxxxxz, g_yy_0_0_xxxxxy, g_yy_0_0_xxxxxyy, g_yy_0_0_xxxxxyz, g_yy_0_0_xxxxxz, g_yy_0_0_xxxxxzz, g_yy_0_0_xxxxyy, g_yy_0_0_xxxxyyy, g_yy_0_0_xxxxyyz, g_yy_0_0_xxxxyz, g_yy_0_0_xxxxyzz, g_yy_0_0_xxxxzz, g_yy_0_0_xxxxzzz, g_yy_0_0_xxxyyy, g_yy_0_0_xxxyyyy, g_yy_0_0_xxxyyyz, g_yy_0_0_xxxyyz, g_yy_0_0_xxxyyzz, g_yy_0_0_xxxyzz, g_yy_0_0_xxxyzzz, g_yy_0_0_xxxzzz, g_yy_0_0_xxxzzzz, g_yy_0_0_xxyyyy, g_yy_0_0_xxyyyyy, g_yy_0_0_xxyyyyz, g_yy_0_0_xxyyyz, g_yy_0_0_xxyyyzz, g_yy_0_0_xxyyzz, g_yy_0_0_xxyyzzz, g_yy_0_0_xxyzzz, g_yy_0_0_xxyzzzz, g_yy_0_0_xxzzzz, g_yy_0_0_xxzzzzz, g_yy_0_0_xyyyyy, g_yy_0_0_xyyyyyy, g_yy_0_0_xyyyyyz, g_yy_0_0_xyyyyz, g_yy_0_0_xyyyyzz, g_yy_0_0_xyyyzz, g_yy_0_0_xyyyzzz, g_yy_0_0_xyyzzz, g_yy_0_0_xyyzzzz, g_yy_0_0_xyzzzz, g_yy_0_0_xyzzzzz, g_yy_0_0_xzzzzz, g_yy_0_0_xzzzzzz, g_yy_0_0_yyyyyy, g_yy_0_0_yyyyyz, g_yy_0_0_yyyyzz, g_yy_0_0_yyyzzz, g_yy_0_0_yyzzzz, g_yy_0_0_yzzzzz, g_yy_0_0_zzzzzz, g_yy_0_x_xxxxxx, g_yy_0_x_xxxxxy, g_yy_0_x_xxxxxz, g_yy_0_x_xxxxyy, g_yy_0_x_xxxxyz, g_yy_0_x_xxxxzz, g_yy_0_x_xxxyyy, g_yy_0_x_xxxyyz, g_yy_0_x_xxxyzz, g_yy_0_x_xxxzzz, g_yy_0_x_xxyyyy, g_yy_0_x_xxyyyz, g_yy_0_x_xxyyzz, g_yy_0_x_xxyzzz, g_yy_0_x_xxzzzz, g_yy_0_x_xyyyyy, g_yy_0_x_xyyyyz, g_yy_0_x_xyyyzz, g_yy_0_x_xyyzzz, g_yy_0_x_xyzzzz, g_yy_0_x_xzzzzz, g_yy_0_x_yyyyyy, g_yy_0_x_yyyyyz, g_yy_0_x_yyyyzz, g_yy_0_x_yyyzzz, g_yy_0_x_yyzzzz, g_yy_0_x_yzzzzz, g_yy_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_x_xxxxxx[k] = -g_yy_0_0_xxxxxx[k] * ab_x + g_yy_0_0_xxxxxxx[k];

                g_yy_0_x_xxxxxy[k] = -g_yy_0_0_xxxxxy[k] * ab_x + g_yy_0_0_xxxxxxy[k];

                g_yy_0_x_xxxxxz[k] = -g_yy_0_0_xxxxxz[k] * ab_x + g_yy_0_0_xxxxxxz[k];

                g_yy_0_x_xxxxyy[k] = -g_yy_0_0_xxxxyy[k] * ab_x + g_yy_0_0_xxxxxyy[k];

                g_yy_0_x_xxxxyz[k] = -g_yy_0_0_xxxxyz[k] * ab_x + g_yy_0_0_xxxxxyz[k];

                g_yy_0_x_xxxxzz[k] = -g_yy_0_0_xxxxzz[k] * ab_x + g_yy_0_0_xxxxxzz[k];

                g_yy_0_x_xxxyyy[k] = -g_yy_0_0_xxxyyy[k] * ab_x + g_yy_0_0_xxxxyyy[k];

                g_yy_0_x_xxxyyz[k] = -g_yy_0_0_xxxyyz[k] * ab_x + g_yy_0_0_xxxxyyz[k];

                g_yy_0_x_xxxyzz[k] = -g_yy_0_0_xxxyzz[k] * ab_x + g_yy_0_0_xxxxyzz[k];

                g_yy_0_x_xxxzzz[k] = -g_yy_0_0_xxxzzz[k] * ab_x + g_yy_0_0_xxxxzzz[k];

                g_yy_0_x_xxyyyy[k] = -g_yy_0_0_xxyyyy[k] * ab_x + g_yy_0_0_xxxyyyy[k];

                g_yy_0_x_xxyyyz[k] = -g_yy_0_0_xxyyyz[k] * ab_x + g_yy_0_0_xxxyyyz[k];

                g_yy_0_x_xxyyzz[k] = -g_yy_0_0_xxyyzz[k] * ab_x + g_yy_0_0_xxxyyzz[k];

                g_yy_0_x_xxyzzz[k] = -g_yy_0_0_xxyzzz[k] * ab_x + g_yy_0_0_xxxyzzz[k];

                g_yy_0_x_xxzzzz[k] = -g_yy_0_0_xxzzzz[k] * ab_x + g_yy_0_0_xxxzzzz[k];

                g_yy_0_x_xyyyyy[k] = -g_yy_0_0_xyyyyy[k] * ab_x + g_yy_0_0_xxyyyyy[k];

                g_yy_0_x_xyyyyz[k] = -g_yy_0_0_xyyyyz[k] * ab_x + g_yy_0_0_xxyyyyz[k];

                g_yy_0_x_xyyyzz[k] = -g_yy_0_0_xyyyzz[k] * ab_x + g_yy_0_0_xxyyyzz[k];

                g_yy_0_x_xyyzzz[k] = -g_yy_0_0_xyyzzz[k] * ab_x + g_yy_0_0_xxyyzzz[k];

                g_yy_0_x_xyzzzz[k] = -g_yy_0_0_xyzzzz[k] * ab_x + g_yy_0_0_xxyzzzz[k];

                g_yy_0_x_xzzzzz[k] = -g_yy_0_0_xzzzzz[k] * ab_x + g_yy_0_0_xxzzzzz[k];

                g_yy_0_x_yyyyyy[k] = -g_yy_0_0_yyyyyy[k] * ab_x + g_yy_0_0_xyyyyyy[k];

                g_yy_0_x_yyyyyz[k] = -g_yy_0_0_yyyyyz[k] * ab_x + g_yy_0_0_xyyyyyz[k];

                g_yy_0_x_yyyyzz[k] = -g_yy_0_0_yyyyzz[k] * ab_x + g_yy_0_0_xyyyyzz[k];

                g_yy_0_x_yyyzzz[k] = -g_yy_0_0_yyyzzz[k] * ab_x + g_yy_0_0_xyyyzzz[k];

                g_yy_0_x_yyzzzz[k] = -g_yy_0_0_yyzzzz[k] * ab_x + g_yy_0_0_xyyzzzz[k];

                g_yy_0_x_yzzzzz[k] = -g_yy_0_0_yzzzzz[k] * ab_x + g_yy_0_0_xyzzzzz[k];

                g_yy_0_x_zzzzzz[k] = -g_yy_0_0_zzzzzz[k] * ab_x + g_yy_0_0_xzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_yy_0_y_xxxxxx = cbuffer.data(pi_geom_20_off + 280 * ccomps * dcomps);

            auto g_yy_0_y_xxxxxy = cbuffer.data(pi_geom_20_off + 281 * ccomps * dcomps);

            auto g_yy_0_y_xxxxxz = cbuffer.data(pi_geom_20_off + 282 * ccomps * dcomps);

            auto g_yy_0_y_xxxxyy = cbuffer.data(pi_geom_20_off + 283 * ccomps * dcomps);

            auto g_yy_0_y_xxxxyz = cbuffer.data(pi_geom_20_off + 284 * ccomps * dcomps);

            auto g_yy_0_y_xxxxzz = cbuffer.data(pi_geom_20_off + 285 * ccomps * dcomps);

            auto g_yy_0_y_xxxyyy = cbuffer.data(pi_geom_20_off + 286 * ccomps * dcomps);

            auto g_yy_0_y_xxxyyz = cbuffer.data(pi_geom_20_off + 287 * ccomps * dcomps);

            auto g_yy_0_y_xxxyzz = cbuffer.data(pi_geom_20_off + 288 * ccomps * dcomps);

            auto g_yy_0_y_xxxzzz = cbuffer.data(pi_geom_20_off + 289 * ccomps * dcomps);

            auto g_yy_0_y_xxyyyy = cbuffer.data(pi_geom_20_off + 290 * ccomps * dcomps);

            auto g_yy_0_y_xxyyyz = cbuffer.data(pi_geom_20_off + 291 * ccomps * dcomps);

            auto g_yy_0_y_xxyyzz = cbuffer.data(pi_geom_20_off + 292 * ccomps * dcomps);

            auto g_yy_0_y_xxyzzz = cbuffer.data(pi_geom_20_off + 293 * ccomps * dcomps);

            auto g_yy_0_y_xxzzzz = cbuffer.data(pi_geom_20_off + 294 * ccomps * dcomps);

            auto g_yy_0_y_xyyyyy = cbuffer.data(pi_geom_20_off + 295 * ccomps * dcomps);

            auto g_yy_0_y_xyyyyz = cbuffer.data(pi_geom_20_off + 296 * ccomps * dcomps);

            auto g_yy_0_y_xyyyzz = cbuffer.data(pi_geom_20_off + 297 * ccomps * dcomps);

            auto g_yy_0_y_xyyzzz = cbuffer.data(pi_geom_20_off + 298 * ccomps * dcomps);

            auto g_yy_0_y_xyzzzz = cbuffer.data(pi_geom_20_off + 299 * ccomps * dcomps);

            auto g_yy_0_y_xzzzzz = cbuffer.data(pi_geom_20_off + 300 * ccomps * dcomps);

            auto g_yy_0_y_yyyyyy = cbuffer.data(pi_geom_20_off + 301 * ccomps * dcomps);

            auto g_yy_0_y_yyyyyz = cbuffer.data(pi_geom_20_off + 302 * ccomps * dcomps);

            auto g_yy_0_y_yyyyzz = cbuffer.data(pi_geom_20_off + 303 * ccomps * dcomps);

            auto g_yy_0_y_yyyzzz = cbuffer.data(pi_geom_20_off + 304 * ccomps * dcomps);

            auto g_yy_0_y_yyzzzz = cbuffer.data(pi_geom_20_off + 305 * ccomps * dcomps);

            auto g_yy_0_y_yzzzzz = cbuffer.data(pi_geom_20_off + 306 * ccomps * dcomps);

            auto g_yy_0_y_zzzzzz = cbuffer.data(pi_geom_20_off + 307 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxxz, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxzz, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzzz, g_yy_0_0_xxxxxx, g_yy_0_0_xxxxxxy, g_yy_0_0_xxxxxy, g_yy_0_0_xxxxxyy, g_yy_0_0_xxxxxyz, g_yy_0_0_xxxxxz, g_yy_0_0_xxxxyy, g_yy_0_0_xxxxyyy, g_yy_0_0_xxxxyyz, g_yy_0_0_xxxxyz, g_yy_0_0_xxxxyzz, g_yy_0_0_xxxxzz, g_yy_0_0_xxxyyy, g_yy_0_0_xxxyyyy, g_yy_0_0_xxxyyyz, g_yy_0_0_xxxyyz, g_yy_0_0_xxxyyzz, g_yy_0_0_xxxyzz, g_yy_0_0_xxxyzzz, g_yy_0_0_xxxzzz, g_yy_0_0_xxyyyy, g_yy_0_0_xxyyyyy, g_yy_0_0_xxyyyyz, g_yy_0_0_xxyyyz, g_yy_0_0_xxyyyzz, g_yy_0_0_xxyyzz, g_yy_0_0_xxyyzzz, g_yy_0_0_xxyzzz, g_yy_0_0_xxyzzzz, g_yy_0_0_xxzzzz, g_yy_0_0_xyyyyy, g_yy_0_0_xyyyyyy, g_yy_0_0_xyyyyyz, g_yy_0_0_xyyyyz, g_yy_0_0_xyyyyzz, g_yy_0_0_xyyyzz, g_yy_0_0_xyyyzzz, g_yy_0_0_xyyzzz, g_yy_0_0_xyyzzzz, g_yy_0_0_xyzzzz, g_yy_0_0_xyzzzzz, g_yy_0_0_xzzzzz, g_yy_0_0_yyyyyy, g_yy_0_0_yyyyyyy, g_yy_0_0_yyyyyyz, g_yy_0_0_yyyyyz, g_yy_0_0_yyyyyzz, g_yy_0_0_yyyyzz, g_yy_0_0_yyyyzzz, g_yy_0_0_yyyzzz, g_yy_0_0_yyyzzzz, g_yy_0_0_yyzzzz, g_yy_0_0_yyzzzzz, g_yy_0_0_yzzzzz, g_yy_0_0_yzzzzzz, g_yy_0_0_zzzzzz, g_yy_0_y_xxxxxx, g_yy_0_y_xxxxxy, g_yy_0_y_xxxxxz, g_yy_0_y_xxxxyy, g_yy_0_y_xxxxyz, g_yy_0_y_xxxxzz, g_yy_0_y_xxxyyy, g_yy_0_y_xxxyyz, g_yy_0_y_xxxyzz, g_yy_0_y_xxxzzz, g_yy_0_y_xxyyyy, g_yy_0_y_xxyyyz, g_yy_0_y_xxyyzz, g_yy_0_y_xxyzzz, g_yy_0_y_xxzzzz, g_yy_0_y_xyyyyy, g_yy_0_y_xyyyyz, g_yy_0_y_xyyyzz, g_yy_0_y_xyyzzz, g_yy_0_y_xyzzzz, g_yy_0_y_xzzzzz, g_yy_0_y_yyyyyy, g_yy_0_y_yyyyyz, g_yy_0_y_yyyyzz, g_yy_0_y_yyyzzz, g_yy_0_y_yyzzzz, g_yy_0_y_yzzzzz, g_yy_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_y_xxxxxx[k] = -2.0 * g_y_0_0_xxxxxx[k] - g_yy_0_0_xxxxxx[k] * ab_y + g_yy_0_0_xxxxxxy[k];

                g_yy_0_y_xxxxxy[k] = -2.0 * g_y_0_0_xxxxxy[k] - g_yy_0_0_xxxxxy[k] * ab_y + g_yy_0_0_xxxxxyy[k];

                g_yy_0_y_xxxxxz[k] = -2.0 * g_y_0_0_xxxxxz[k] - g_yy_0_0_xxxxxz[k] * ab_y + g_yy_0_0_xxxxxyz[k];

                g_yy_0_y_xxxxyy[k] = -2.0 * g_y_0_0_xxxxyy[k] - g_yy_0_0_xxxxyy[k] * ab_y + g_yy_0_0_xxxxyyy[k];

                g_yy_0_y_xxxxyz[k] = -2.0 * g_y_0_0_xxxxyz[k] - g_yy_0_0_xxxxyz[k] * ab_y + g_yy_0_0_xxxxyyz[k];

                g_yy_0_y_xxxxzz[k] = -2.0 * g_y_0_0_xxxxzz[k] - g_yy_0_0_xxxxzz[k] * ab_y + g_yy_0_0_xxxxyzz[k];

                g_yy_0_y_xxxyyy[k] = -2.0 * g_y_0_0_xxxyyy[k] - g_yy_0_0_xxxyyy[k] * ab_y + g_yy_0_0_xxxyyyy[k];

                g_yy_0_y_xxxyyz[k] = -2.0 * g_y_0_0_xxxyyz[k] - g_yy_0_0_xxxyyz[k] * ab_y + g_yy_0_0_xxxyyyz[k];

                g_yy_0_y_xxxyzz[k] = -2.0 * g_y_0_0_xxxyzz[k] - g_yy_0_0_xxxyzz[k] * ab_y + g_yy_0_0_xxxyyzz[k];

                g_yy_0_y_xxxzzz[k] = -2.0 * g_y_0_0_xxxzzz[k] - g_yy_0_0_xxxzzz[k] * ab_y + g_yy_0_0_xxxyzzz[k];

                g_yy_0_y_xxyyyy[k] = -2.0 * g_y_0_0_xxyyyy[k] - g_yy_0_0_xxyyyy[k] * ab_y + g_yy_0_0_xxyyyyy[k];

                g_yy_0_y_xxyyyz[k] = -2.0 * g_y_0_0_xxyyyz[k] - g_yy_0_0_xxyyyz[k] * ab_y + g_yy_0_0_xxyyyyz[k];

                g_yy_0_y_xxyyzz[k] = -2.0 * g_y_0_0_xxyyzz[k] - g_yy_0_0_xxyyzz[k] * ab_y + g_yy_0_0_xxyyyzz[k];

                g_yy_0_y_xxyzzz[k] = -2.0 * g_y_0_0_xxyzzz[k] - g_yy_0_0_xxyzzz[k] * ab_y + g_yy_0_0_xxyyzzz[k];

                g_yy_0_y_xxzzzz[k] = -2.0 * g_y_0_0_xxzzzz[k] - g_yy_0_0_xxzzzz[k] * ab_y + g_yy_0_0_xxyzzzz[k];

                g_yy_0_y_xyyyyy[k] = -2.0 * g_y_0_0_xyyyyy[k] - g_yy_0_0_xyyyyy[k] * ab_y + g_yy_0_0_xyyyyyy[k];

                g_yy_0_y_xyyyyz[k] = -2.0 * g_y_0_0_xyyyyz[k] - g_yy_0_0_xyyyyz[k] * ab_y + g_yy_0_0_xyyyyyz[k];

                g_yy_0_y_xyyyzz[k] = -2.0 * g_y_0_0_xyyyzz[k] - g_yy_0_0_xyyyzz[k] * ab_y + g_yy_0_0_xyyyyzz[k];

                g_yy_0_y_xyyzzz[k] = -2.0 * g_y_0_0_xyyzzz[k] - g_yy_0_0_xyyzzz[k] * ab_y + g_yy_0_0_xyyyzzz[k];

                g_yy_0_y_xyzzzz[k] = -2.0 * g_y_0_0_xyzzzz[k] - g_yy_0_0_xyzzzz[k] * ab_y + g_yy_0_0_xyyzzzz[k];

                g_yy_0_y_xzzzzz[k] = -2.0 * g_y_0_0_xzzzzz[k] - g_yy_0_0_xzzzzz[k] * ab_y + g_yy_0_0_xyzzzzz[k];

                g_yy_0_y_yyyyyy[k] = -2.0 * g_y_0_0_yyyyyy[k] - g_yy_0_0_yyyyyy[k] * ab_y + g_yy_0_0_yyyyyyy[k];

                g_yy_0_y_yyyyyz[k] = -2.0 * g_y_0_0_yyyyyz[k] - g_yy_0_0_yyyyyz[k] * ab_y + g_yy_0_0_yyyyyyz[k];

                g_yy_0_y_yyyyzz[k] = -2.0 * g_y_0_0_yyyyzz[k] - g_yy_0_0_yyyyzz[k] * ab_y + g_yy_0_0_yyyyyzz[k];

                g_yy_0_y_yyyzzz[k] = -2.0 * g_y_0_0_yyyzzz[k] - g_yy_0_0_yyyzzz[k] * ab_y + g_yy_0_0_yyyyzzz[k];

                g_yy_0_y_yyzzzz[k] = -2.0 * g_y_0_0_yyzzzz[k] - g_yy_0_0_yyzzzz[k] * ab_y + g_yy_0_0_yyyzzzz[k];

                g_yy_0_y_yzzzzz[k] = -2.0 * g_y_0_0_yzzzzz[k] - g_yy_0_0_yzzzzz[k] * ab_y + g_yy_0_0_yyzzzzz[k];

                g_yy_0_y_zzzzzz[k] = -2.0 * g_y_0_0_zzzzzz[k] - g_yy_0_0_zzzzzz[k] * ab_y + g_yy_0_0_yzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_yy_0_z_xxxxxx = cbuffer.data(pi_geom_20_off + 308 * ccomps * dcomps);

            auto g_yy_0_z_xxxxxy = cbuffer.data(pi_geom_20_off + 309 * ccomps * dcomps);

            auto g_yy_0_z_xxxxxz = cbuffer.data(pi_geom_20_off + 310 * ccomps * dcomps);

            auto g_yy_0_z_xxxxyy = cbuffer.data(pi_geom_20_off + 311 * ccomps * dcomps);

            auto g_yy_0_z_xxxxyz = cbuffer.data(pi_geom_20_off + 312 * ccomps * dcomps);

            auto g_yy_0_z_xxxxzz = cbuffer.data(pi_geom_20_off + 313 * ccomps * dcomps);

            auto g_yy_0_z_xxxyyy = cbuffer.data(pi_geom_20_off + 314 * ccomps * dcomps);

            auto g_yy_0_z_xxxyyz = cbuffer.data(pi_geom_20_off + 315 * ccomps * dcomps);

            auto g_yy_0_z_xxxyzz = cbuffer.data(pi_geom_20_off + 316 * ccomps * dcomps);

            auto g_yy_0_z_xxxzzz = cbuffer.data(pi_geom_20_off + 317 * ccomps * dcomps);

            auto g_yy_0_z_xxyyyy = cbuffer.data(pi_geom_20_off + 318 * ccomps * dcomps);

            auto g_yy_0_z_xxyyyz = cbuffer.data(pi_geom_20_off + 319 * ccomps * dcomps);

            auto g_yy_0_z_xxyyzz = cbuffer.data(pi_geom_20_off + 320 * ccomps * dcomps);

            auto g_yy_0_z_xxyzzz = cbuffer.data(pi_geom_20_off + 321 * ccomps * dcomps);

            auto g_yy_0_z_xxzzzz = cbuffer.data(pi_geom_20_off + 322 * ccomps * dcomps);

            auto g_yy_0_z_xyyyyy = cbuffer.data(pi_geom_20_off + 323 * ccomps * dcomps);

            auto g_yy_0_z_xyyyyz = cbuffer.data(pi_geom_20_off + 324 * ccomps * dcomps);

            auto g_yy_0_z_xyyyzz = cbuffer.data(pi_geom_20_off + 325 * ccomps * dcomps);

            auto g_yy_0_z_xyyzzz = cbuffer.data(pi_geom_20_off + 326 * ccomps * dcomps);

            auto g_yy_0_z_xyzzzz = cbuffer.data(pi_geom_20_off + 327 * ccomps * dcomps);

            auto g_yy_0_z_xzzzzz = cbuffer.data(pi_geom_20_off + 328 * ccomps * dcomps);

            auto g_yy_0_z_yyyyyy = cbuffer.data(pi_geom_20_off + 329 * ccomps * dcomps);

            auto g_yy_0_z_yyyyyz = cbuffer.data(pi_geom_20_off + 330 * ccomps * dcomps);

            auto g_yy_0_z_yyyyzz = cbuffer.data(pi_geom_20_off + 331 * ccomps * dcomps);

            auto g_yy_0_z_yyyzzz = cbuffer.data(pi_geom_20_off + 332 * ccomps * dcomps);

            auto g_yy_0_z_yyzzzz = cbuffer.data(pi_geom_20_off + 333 * ccomps * dcomps);

            auto g_yy_0_z_yzzzzz = cbuffer.data(pi_geom_20_off + 334 * ccomps * dcomps);

            auto g_yy_0_z_zzzzzz = cbuffer.data(pi_geom_20_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_0_xxxxxx, g_yy_0_0_xxxxxxz, g_yy_0_0_xxxxxy, g_yy_0_0_xxxxxyz, g_yy_0_0_xxxxxz, g_yy_0_0_xxxxxzz, g_yy_0_0_xxxxyy, g_yy_0_0_xxxxyyz, g_yy_0_0_xxxxyz, g_yy_0_0_xxxxyzz, g_yy_0_0_xxxxzz, g_yy_0_0_xxxxzzz, g_yy_0_0_xxxyyy, g_yy_0_0_xxxyyyz, g_yy_0_0_xxxyyz, g_yy_0_0_xxxyyzz, g_yy_0_0_xxxyzz, g_yy_0_0_xxxyzzz, g_yy_0_0_xxxzzz, g_yy_0_0_xxxzzzz, g_yy_0_0_xxyyyy, g_yy_0_0_xxyyyyz, g_yy_0_0_xxyyyz, g_yy_0_0_xxyyyzz, g_yy_0_0_xxyyzz, g_yy_0_0_xxyyzzz, g_yy_0_0_xxyzzz, g_yy_0_0_xxyzzzz, g_yy_0_0_xxzzzz, g_yy_0_0_xxzzzzz, g_yy_0_0_xyyyyy, g_yy_0_0_xyyyyyz, g_yy_0_0_xyyyyz, g_yy_0_0_xyyyyzz, g_yy_0_0_xyyyzz, g_yy_0_0_xyyyzzz, g_yy_0_0_xyyzzz, g_yy_0_0_xyyzzzz, g_yy_0_0_xyzzzz, g_yy_0_0_xyzzzzz, g_yy_0_0_xzzzzz, g_yy_0_0_xzzzzzz, g_yy_0_0_yyyyyy, g_yy_0_0_yyyyyyz, g_yy_0_0_yyyyyz, g_yy_0_0_yyyyyzz, g_yy_0_0_yyyyzz, g_yy_0_0_yyyyzzz, g_yy_0_0_yyyzzz, g_yy_0_0_yyyzzzz, g_yy_0_0_yyzzzz, g_yy_0_0_yyzzzzz, g_yy_0_0_yzzzzz, g_yy_0_0_yzzzzzz, g_yy_0_0_zzzzzz, g_yy_0_0_zzzzzzz, g_yy_0_z_xxxxxx, g_yy_0_z_xxxxxy, g_yy_0_z_xxxxxz, g_yy_0_z_xxxxyy, g_yy_0_z_xxxxyz, g_yy_0_z_xxxxzz, g_yy_0_z_xxxyyy, g_yy_0_z_xxxyyz, g_yy_0_z_xxxyzz, g_yy_0_z_xxxzzz, g_yy_0_z_xxyyyy, g_yy_0_z_xxyyyz, g_yy_0_z_xxyyzz, g_yy_0_z_xxyzzz, g_yy_0_z_xxzzzz, g_yy_0_z_xyyyyy, g_yy_0_z_xyyyyz, g_yy_0_z_xyyyzz, g_yy_0_z_xyyzzz, g_yy_0_z_xyzzzz, g_yy_0_z_xzzzzz, g_yy_0_z_yyyyyy, g_yy_0_z_yyyyyz, g_yy_0_z_yyyyzz, g_yy_0_z_yyyzzz, g_yy_0_z_yyzzzz, g_yy_0_z_yzzzzz, g_yy_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_z_xxxxxx[k] = -g_yy_0_0_xxxxxx[k] * ab_z + g_yy_0_0_xxxxxxz[k];

                g_yy_0_z_xxxxxy[k] = -g_yy_0_0_xxxxxy[k] * ab_z + g_yy_0_0_xxxxxyz[k];

                g_yy_0_z_xxxxxz[k] = -g_yy_0_0_xxxxxz[k] * ab_z + g_yy_0_0_xxxxxzz[k];

                g_yy_0_z_xxxxyy[k] = -g_yy_0_0_xxxxyy[k] * ab_z + g_yy_0_0_xxxxyyz[k];

                g_yy_0_z_xxxxyz[k] = -g_yy_0_0_xxxxyz[k] * ab_z + g_yy_0_0_xxxxyzz[k];

                g_yy_0_z_xxxxzz[k] = -g_yy_0_0_xxxxzz[k] * ab_z + g_yy_0_0_xxxxzzz[k];

                g_yy_0_z_xxxyyy[k] = -g_yy_0_0_xxxyyy[k] * ab_z + g_yy_0_0_xxxyyyz[k];

                g_yy_0_z_xxxyyz[k] = -g_yy_0_0_xxxyyz[k] * ab_z + g_yy_0_0_xxxyyzz[k];

                g_yy_0_z_xxxyzz[k] = -g_yy_0_0_xxxyzz[k] * ab_z + g_yy_0_0_xxxyzzz[k];

                g_yy_0_z_xxxzzz[k] = -g_yy_0_0_xxxzzz[k] * ab_z + g_yy_0_0_xxxzzzz[k];

                g_yy_0_z_xxyyyy[k] = -g_yy_0_0_xxyyyy[k] * ab_z + g_yy_0_0_xxyyyyz[k];

                g_yy_0_z_xxyyyz[k] = -g_yy_0_0_xxyyyz[k] * ab_z + g_yy_0_0_xxyyyzz[k];

                g_yy_0_z_xxyyzz[k] = -g_yy_0_0_xxyyzz[k] * ab_z + g_yy_0_0_xxyyzzz[k];

                g_yy_0_z_xxyzzz[k] = -g_yy_0_0_xxyzzz[k] * ab_z + g_yy_0_0_xxyzzzz[k];

                g_yy_0_z_xxzzzz[k] = -g_yy_0_0_xxzzzz[k] * ab_z + g_yy_0_0_xxzzzzz[k];

                g_yy_0_z_xyyyyy[k] = -g_yy_0_0_xyyyyy[k] * ab_z + g_yy_0_0_xyyyyyz[k];

                g_yy_0_z_xyyyyz[k] = -g_yy_0_0_xyyyyz[k] * ab_z + g_yy_0_0_xyyyyzz[k];

                g_yy_0_z_xyyyzz[k] = -g_yy_0_0_xyyyzz[k] * ab_z + g_yy_0_0_xyyyzzz[k];

                g_yy_0_z_xyyzzz[k] = -g_yy_0_0_xyyzzz[k] * ab_z + g_yy_0_0_xyyzzzz[k];

                g_yy_0_z_xyzzzz[k] = -g_yy_0_0_xyzzzz[k] * ab_z + g_yy_0_0_xyzzzzz[k];

                g_yy_0_z_xzzzzz[k] = -g_yy_0_0_xzzzzz[k] * ab_z + g_yy_0_0_xzzzzzz[k];

                g_yy_0_z_yyyyyy[k] = -g_yy_0_0_yyyyyy[k] * ab_z + g_yy_0_0_yyyyyyz[k];

                g_yy_0_z_yyyyyz[k] = -g_yy_0_0_yyyyyz[k] * ab_z + g_yy_0_0_yyyyyzz[k];

                g_yy_0_z_yyyyzz[k] = -g_yy_0_0_yyyyzz[k] * ab_z + g_yy_0_0_yyyyzzz[k];

                g_yy_0_z_yyyzzz[k] = -g_yy_0_0_yyyzzz[k] * ab_z + g_yy_0_0_yyyzzzz[k];

                g_yy_0_z_yyzzzz[k] = -g_yy_0_0_yyzzzz[k] * ab_z + g_yy_0_0_yyzzzzz[k];

                g_yy_0_z_yzzzzz[k] = -g_yy_0_0_yzzzzz[k] * ab_z + g_yy_0_0_yzzzzzz[k];

                g_yy_0_z_zzzzzz[k] = -g_yy_0_0_zzzzzz[k] * ab_z + g_yy_0_0_zzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_yz_0_x_xxxxxx = cbuffer.data(pi_geom_20_off + 336 * ccomps * dcomps);

            auto g_yz_0_x_xxxxxy = cbuffer.data(pi_geom_20_off + 337 * ccomps * dcomps);

            auto g_yz_0_x_xxxxxz = cbuffer.data(pi_geom_20_off + 338 * ccomps * dcomps);

            auto g_yz_0_x_xxxxyy = cbuffer.data(pi_geom_20_off + 339 * ccomps * dcomps);

            auto g_yz_0_x_xxxxyz = cbuffer.data(pi_geom_20_off + 340 * ccomps * dcomps);

            auto g_yz_0_x_xxxxzz = cbuffer.data(pi_geom_20_off + 341 * ccomps * dcomps);

            auto g_yz_0_x_xxxyyy = cbuffer.data(pi_geom_20_off + 342 * ccomps * dcomps);

            auto g_yz_0_x_xxxyyz = cbuffer.data(pi_geom_20_off + 343 * ccomps * dcomps);

            auto g_yz_0_x_xxxyzz = cbuffer.data(pi_geom_20_off + 344 * ccomps * dcomps);

            auto g_yz_0_x_xxxzzz = cbuffer.data(pi_geom_20_off + 345 * ccomps * dcomps);

            auto g_yz_0_x_xxyyyy = cbuffer.data(pi_geom_20_off + 346 * ccomps * dcomps);

            auto g_yz_0_x_xxyyyz = cbuffer.data(pi_geom_20_off + 347 * ccomps * dcomps);

            auto g_yz_0_x_xxyyzz = cbuffer.data(pi_geom_20_off + 348 * ccomps * dcomps);

            auto g_yz_0_x_xxyzzz = cbuffer.data(pi_geom_20_off + 349 * ccomps * dcomps);

            auto g_yz_0_x_xxzzzz = cbuffer.data(pi_geom_20_off + 350 * ccomps * dcomps);

            auto g_yz_0_x_xyyyyy = cbuffer.data(pi_geom_20_off + 351 * ccomps * dcomps);

            auto g_yz_0_x_xyyyyz = cbuffer.data(pi_geom_20_off + 352 * ccomps * dcomps);

            auto g_yz_0_x_xyyyzz = cbuffer.data(pi_geom_20_off + 353 * ccomps * dcomps);

            auto g_yz_0_x_xyyzzz = cbuffer.data(pi_geom_20_off + 354 * ccomps * dcomps);

            auto g_yz_0_x_xyzzzz = cbuffer.data(pi_geom_20_off + 355 * ccomps * dcomps);

            auto g_yz_0_x_xzzzzz = cbuffer.data(pi_geom_20_off + 356 * ccomps * dcomps);

            auto g_yz_0_x_yyyyyy = cbuffer.data(pi_geom_20_off + 357 * ccomps * dcomps);

            auto g_yz_0_x_yyyyyz = cbuffer.data(pi_geom_20_off + 358 * ccomps * dcomps);

            auto g_yz_0_x_yyyyzz = cbuffer.data(pi_geom_20_off + 359 * ccomps * dcomps);

            auto g_yz_0_x_yyyzzz = cbuffer.data(pi_geom_20_off + 360 * ccomps * dcomps);

            auto g_yz_0_x_yyzzzz = cbuffer.data(pi_geom_20_off + 361 * ccomps * dcomps);

            auto g_yz_0_x_yzzzzz = cbuffer.data(pi_geom_20_off + 362 * ccomps * dcomps);

            auto g_yz_0_x_zzzzzz = cbuffer.data(pi_geom_20_off + 363 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_0_xxxxxx, g_yz_0_0_xxxxxxx, g_yz_0_0_xxxxxxy, g_yz_0_0_xxxxxxz, g_yz_0_0_xxxxxy, g_yz_0_0_xxxxxyy, g_yz_0_0_xxxxxyz, g_yz_0_0_xxxxxz, g_yz_0_0_xxxxxzz, g_yz_0_0_xxxxyy, g_yz_0_0_xxxxyyy, g_yz_0_0_xxxxyyz, g_yz_0_0_xxxxyz, g_yz_0_0_xxxxyzz, g_yz_0_0_xxxxzz, g_yz_0_0_xxxxzzz, g_yz_0_0_xxxyyy, g_yz_0_0_xxxyyyy, g_yz_0_0_xxxyyyz, g_yz_0_0_xxxyyz, g_yz_0_0_xxxyyzz, g_yz_0_0_xxxyzz, g_yz_0_0_xxxyzzz, g_yz_0_0_xxxzzz, g_yz_0_0_xxxzzzz, g_yz_0_0_xxyyyy, g_yz_0_0_xxyyyyy, g_yz_0_0_xxyyyyz, g_yz_0_0_xxyyyz, g_yz_0_0_xxyyyzz, g_yz_0_0_xxyyzz, g_yz_0_0_xxyyzzz, g_yz_0_0_xxyzzz, g_yz_0_0_xxyzzzz, g_yz_0_0_xxzzzz, g_yz_0_0_xxzzzzz, g_yz_0_0_xyyyyy, g_yz_0_0_xyyyyyy, g_yz_0_0_xyyyyyz, g_yz_0_0_xyyyyz, g_yz_0_0_xyyyyzz, g_yz_0_0_xyyyzz, g_yz_0_0_xyyyzzz, g_yz_0_0_xyyzzz, g_yz_0_0_xyyzzzz, g_yz_0_0_xyzzzz, g_yz_0_0_xyzzzzz, g_yz_0_0_xzzzzz, g_yz_0_0_xzzzzzz, g_yz_0_0_yyyyyy, g_yz_0_0_yyyyyz, g_yz_0_0_yyyyzz, g_yz_0_0_yyyzzz, g_yz_0_0_yyzzzz, g_yz_0_0_yzzzzz, g_yz_0_0_zzzzzz, g_yz_0_x_xxxxxx, g_yz_0_x_xxxxxy, g_yz_0_x_xxxxxz, g_yz_0_x_xxxxyy, g_yz_0_x_xxxxyz, g_yz_0_x_xxxxzz, g_yz_0_x_xxxyyy, g_yz_0_x_xxxyyz, g_yz_0_x_xxxyzz, g_yz_0_x_xxxzzz, g_yz_0_x_xxyyyy, g_yz_0_x_xxyyyz, g_yz_0_x_xxyyzz, g_yz_0_x_xxyzzz, g_yz_0_x_xxzzzz, g_yz_0_x_xyyyyy, g_yz_0_x_xyyyyz, g_yz_0_x_xyyyzz, g_yz_0_x_xyyzzz, g_yz_0_x_xyzzzz, g_yz_0_x_xzzzzz, g_yz_0_x_yyyyyy, g_yz_0_x_yyyyyz, g_yz_0_x_yyyyzz, g_yz_0_x_yyyzzz, g_yz_0_x_yyzzzz, g_yz_0_x_yzzzzz, g_yz_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_x_xxxxxx[k] = -g_yz_0_0_xxxxxx[k] * ab_x + g_yz_0_0_xxxxxxx[k];

                g_yz_0_x_xxxxxy[k] = -g_yz_0_0_xxxxxy[k] * ab_x + g_yz_0_0_xxxxxxy[k];

                g_yz_0_x_xxxxxz[k] = -g_yz_0_0_xxxxxz[k] * ab_x + g_yz_0_0_xxxxxxz[k];

                g_yz_0_x_xxxxyy[k] = -g_yz_0_0_xxxxyy[k] * ab_x + g_yz_0_0_xxxxxyy[k];

                g_yz_0_x_xxxxyz[k] = -g_yz_0_0_xxxxyz[k] * ab_x + g_yz_0_0_xxxxxyz[k];

                g_yz_0_x_xxxxzz[k] = -g_yz_0_0_xxxxzz[k] * ab_x + g_yz_0_0_xxxxxzz[k];

                g_yz_0_x_xxxyyy[k] = -g_yz_0_0_xxxyyy[k] * ab_x + g_yz_0_0_xxxxyyy[k];

                g_yz_0_x_xxxyyz[k] = -g_yz_0_0_xxxyyz[k] * ab_x + g_yz_0_0_xxxxyyz[k];

                g_yz_0_x_xxxyzz[k] = -g_yz_0_0_xxxyzz[k] * ab_x + g_yz_0_0_xxxxyzz[k];

                g_yz_0_x_xxxzzz[k] = -g_yz_0_0_xxxzzz[k] * ab_x + g_yz_0_0_xxxxzzz[k];

                g_yz_0_x_xxyyyy[k] = -g_yz_0_0_xxyyyy[k] * ab_x + g_yz_0_0_xxxyyyy[k];

                g_yz_0_x_xxyyyz[k] = -g_yz_0_0_xxyyyz[k] * ab_x + g_yz_0_0_xxxyyyz[k];

                g_yz_0_x_xxyyzz[k] = -g_yz_0_0_xxyyzz[k] * ab_x + g_yz_0_0_xxxyyzz[k];

                g_yz_0_x_xxyzzz[k] = -g_yz_0_0_xxyzzz[k] * ab_x + g_yz_0_0_xxxyzzz[k];

                g_yz_0_x_xxzzzz[k] = -g_yz_0_0_xxzzzz[k] * ab_x + g_yz_0_0_xxxzzzz[k];

                g_yz_0_x_xyyyyy[k] = -g_yz_0_0_xyyyyy[k] * ab_x + g_yz_0_0_xxyyyyy[k];

                g_yz_0_x_xyyyyz[k] = -g_yz_0_0_xyyyyz[k] * ab_x + g_yz_0_0_xxyyyyz[k];

                g_yz_0_x_xyyyzz[k] = -g_yz_0_0_xyyyzz[k] * ab_x + g_yz_0_0_xxyyyzz[k];

                g_yz_0_x_xyyzzz[k] = -g_yz_0_0_xyyzzz[k] * ab_x + g_yz_0_0_xxyyzzz[k];

                g_yz_0_x_xyzzzz[k] = -g_yz_0_0_xyzzzz[k] * ab_x + g_yz_0_0_xxyzzzz[k];

                g_yz_0_x_xzzzzz[k] = -g_yz_0_0_xzzzzz[k] * ab_x + g_yz_0_0_xxzzzzz[k];

                g_yz_0_x_yyyyyy[k] = -g_yz_0_0_yyyyyy[k] * ab_x + g_yz_0_0_xyyyyyy[k];

                g_yz_0_x_yyyyyz[k] = -g_yz_0_0_yyyyyz[k] * ab_x + g_yz_0_0_xyyyyyz[k];

                g_yz_0_x_yyyyzz[k] = -g_yz_0_0_yyyyzz[k] * ab_x + g_yz_0_0_xyyyyzz[k];

                g_yz_0_x_yyyzzz[k] = -g_yz_0_0_yyyzzz[k] * ab_x + g_yz_0_0_xyyyzzz[k];

                g_yz_0_x_yyzzzz[k] = -g_yz_0_0_yyzzzz[k] * ab_x + g_yz_0_0_xyyzzzz[k];

                g_yz_0_x_yzzzzz[k] = -g_yz_0_0_yzzzzz[k] * ab_x + g_yz_0_0_xyzzzzz[k];

                g_yz_0_x_zzzzzz[k] = -g_yz_0_0_zzzzzz[k] * ab_x + g_yz_0_0_xzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_yz_0_y_xxxxxx = cbuffer.data(pi_geom_20_off + 364 * ccomps * dcomps);

            auto g_yz_0_y_xxxxxy = cbuffer.data(pi_geom_20_off + 365 * ccomps * dcomps);

            auto g_yz_0_y_xxxxxz = cbuffer.data(pi_geom_20_off + 366 * ccomps * dcomps);

            auto g_yz_0_y_xxxxyy = cbuffer.data(pi_geom_20_off + 367 * ccomps * dcomps);

            auto g_yz_0_y_xxxxyz = cbuffer.data(pi_geom_20_off + 368 * ccomps * dcomps);

            auto g_yz_0_y_xxxxzz = cbuffer.data(pi_geom_20_off + 369 * ccomps * dcomps);

            auto g_yz_0_y_xxxyyy = cbuffer.data(pi_geom_20_off + 370 * ccomps * dcomps);

            auto g_yz_0_y_xxxyyz = cbuffer.data(pi_geom_20_off + 371 * ccomps * dcomps);

            auto g_yz_0_y_xxxyzz = cbuffer.data(pi_geom_20_off + 372 * ccomps * dcomps);

            auto g_yz_0_y_xxxzzz = cbuffer.data(pi_geom_20_off + 373 * ccomps * dcomps);

            auto g_yz_0_y_xxyyyy = cbuffer.data(pi_geom_20_off + 374 * ccomps * dcomps);

            auto g_yz_0_y_xxyyyz = cbuffer.data(pi_geom_20_off + 375 * ccomps * dcomps);

            auto g_yz_0_y_xxyyzz = cbuffer.data(pi_geom_20_off + 376 * ccomps * dcomps);

            auto g_yz_0_y_xxyzzz = cbuffer.data(pi_geom_20_off + 377 * ccomps * dcomps);

            auto g_yz_0_y_xxzzzz = cbuffer.data(pi_geom_20_off + 378 * ccomps * dcomps);

            auto g_yz_0_y_xyyyyy = cbuffer.data(pi_geom_20_off + 379 * ccomps * dcomps);

            auto g_yz_0_y_xyyyyz = cbuffer.data(pi_geom_20_off + 380 * ccomps * dcomps);

            auto g_yz_0_y_xyyyzz = cbuffer.data(pi_geom_20_off + 381 * ccomps * dcomps);

            auto g_yz_0_y_xyyzzz = cbuffer.data(pi_geom_20_off + 382 * ccomps * dcomps);

            auto g_yz_0_y_xyzzzz = cbuffer.data(pi_geom_20_off + 383 * ccomps * dcomps);

            auto g_yz_0_y_xzzzzz = cbuffer.data(pi_geom_20_off + 384 * ccomps * dcomps);

            auto g_yz_0_y_yyyyyy = cbuffer.data(pi_geom_20_off + 385 * ccomps * dcomps);

            auto g_yz_0_y_yyyyyz = cbuffer.data(pi_geom_20_off + 386 * ccomps * dcomps);

            auto g_yz_0_y_yyyyzz = cbuffer.data(pi_geom_20_off + 387 * ccomps * dcomps);

            auto g_yz_0_y_yyyzzz = cbuffer.data(pi_geom_20_off + 388 * ccomps * dcomps);

            auto g_yz_0_y_yyzzzz = cbuffer.data(pi_geom_20_off + 389 * ccomps * dcomps);

            auto g_yz_0_y_yzzzzz = cbuffer.data(pi_geom_20_off + 390 * ccomps * dcomps);

            auto g_yz_0_y_zzzzzz = cbuffer.data(pi_geom_20_off + 391 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_0_xxxxxx, g_yz_0_0_xxxxxxy, g_yz_0_0_xxxxxy, g_yz_0_0_xxxxxyy, g_yz_0_0_xxxxxyz, g_yz_0_0_xxxxxz, g_yz_0_0_xxxxyy, g_yz_0_0_xxxxyyy, g_yz_0_0_xxxxyyz, g_yz_0_0_xxxxyz, g_yz_0_0_xxxxyzz, g_yz_0_0_xxxxzz, g_yz_0_0_xxxyyy, g_yz_0_0_xxxyyyy, g_yz_0_0_xxxyyyz, g_yz_0_0_xxxyyz, g_yz_0_0_xxxyyzz, g_yz_0_0_xxxyzz, g_yz_0_0_xxxyzzz, g_yz_0_0_xxxzzz, g_yz_0_0_xxyyyy, g_yz_0_0_xxyyyyy, g_yz_0_0_xxyyyyz, g_yz_0_0_xxyyyz, g_yz_0_0_xxyyyzz, g_yz_0_0_xxyyzz, g_yz_0_0_xxyyzzz, g_yz_0_0_xxyzzz, g_yz_0_0_xxyzzzz, g_yz_0_0_xxzzzz, g_yz_0_0_xyyyyy, g_yz_0_0_xyyyyyy, g_yz_0_0_xyyyyyz, g_yz_0_0_xyyyyz, g_yz_0_0_xyyyyzz, g_yz_0_0_xyyyzz, g_yz_0_0_xyyyzzz, g_yz_0_0_xyyzzz, g_yz_0_0_xyyzzzz, g_yz_0_0_xyzzzz, g_yz_0_0_xyzzzzz, g_yz_0_0_xzzzzz, g_yz_0_0_yyyyyy, g_yz_0_0_yyyyyyy, g_yz_0_0_yyyyyyz, g_yz_0_0_yyyyyz, g_yz_0_0_yyyyyzz, g_yz_0_0_yyyyzz, g_yz_0_0_yyyyzzz, g_yz_0_0_yyyzzz, g_yz_0_0_yyyzzzz, g_yz_0_0_yyzzzz, g_yz_0_0_yyzzzzz, g_yz_0_0_yzzzzz, g_yz_0_0_yzzzzzz, g_yz_0_0_zzzzzz, g_yz_0_y_xxxxxx, g_yz_0_y_xxxxxy, g_yz_0_y_xxxxxz, g_yz_0_y_xxxxyy, g_yz_0_y_xxxxyz, g_yz_0_y_xxxxzz, g_yz_0_y_xxxyyy, g_yz_0_y_xxxyyz, g_yz_0_y_xxxyzz, g_yz_0_y_xxxzzz, g_yz_0_y_xxyyyy, g_yz_0_y_xxyyyz, g_yz_0_y_xxyyzz, g_yz_0_y_xxyzzz, g_yz_0_y_xxzzzz, g_yz_0_y_xyyyyy, g_yz_0_y_xyyyyz, g_yz_0_y_xyyyzz, g_yz_0_y_xyyzzz, g_yz_0_y_xyzzzz, g_yz_0_y_xzzzzz, g_yz_0_y_yyyyyy, g_yz_0_y_yyyyyz, g_yz_0_y_yyyyzz, g_yz_0_y_yyyzzz, g_yz_0_y_yyzzzz, g_yz_0_y_yzzzzz, g_yz_0_y_zzzzzz, g_z_0_0_xxxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxxz, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxzz, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_y_xxxxxx[k] = -g_z_0_0_xxxxxx[k] - g_yz_0_0_xxxxxx[k] * ab_y + g_yz_0_0_xxxxxxy[k];

                g_yz_0_y_xxxxxy[k] = -g_z_0_0_xxxxxy[k] - g_yz_0_0_xxxxxy[k] * ab_y + g_yz_0_0_xxxxxyy[k];

                g_yz_0_y_xxxxxz[k] = -g_z_0_0_xxxxxz[k] - g_yz_0_0_xxxxxz[k] * ab_y + g_yz_0_0_xxxxxyz[k];

                g_yz_0_y_xxxxyy[k] = -g_z_0_0_xxxxyy[k] - g_yz_0_0_xxxxyy[k] * ab_y + g_yz_0_0_xxxxyyy[k];

                g_yz_0_y_xxxxyz[k] = -g_z_0_0_xxxxyz[k] - g_yz_0_0_xxxxyz[k] * ab_y + g_yz_0_0_xxxxyyz[k];

                g_yz_0_y_xxxxzz[k] = -g_z_0_0_xxxxzz[k] - g_yz_0_0_xxxxzz[k] * ab_y + g_yz_0_0_xxxxyzz[k];

                g_yz_0_y_xxxyyy[k] = -g_z_0_0_xxxyyy[k] - g_yz_0_0_xxxyyy[k] * ab_y + g_yz_0_0_xxxyyyy[k];

                g_yz_0_y_xxxyyz[k] = -g_z_0_0_xxxyyz[k] - g_yz_0_0_xxxyyz[k] * ab_y + g_yz_0_0_xxxyyyz[k];

                g_yz_0_y_xxxyzz[k] = -g_z_0_0_xxxyzz[k] - g_yz_0_0_xxxyzz[k] * ab_y + g_yz_0_0_xxxyyzz[k];

                g_yz_0_y_xxxzzz[k] = -g_z_0_0_xxxzzz[k] - g_yz_0_0_xxxzzz[k] * ab_y + g_yz_0_0_xxxyzzz[k];

                g_yz_0_y_xxyyyy[k] = -g_z_0_0_xxyyyy[k] - g_yz_0_0_xxyyyy[k] * ab_y + g_yz_0_0_xxyyyyy[k];

                g_yz_0_y_xxyyyz[k] = -g_z_0_0_xxyyyz[k] - g_yz_0_0_xxyyyz[k] * ab_y + g_yz_0_0_xxyyyyz[k];

                g_yz_0_y_xxyyzz[k] = -g_z_0_0_xxyyzz[k] - g_yz_0_0_xxyyzz[k] * ab_y + g_yz_0_0_xxyyyzz[k];

                g_yz_0_y_xxyzzz[k] = -g_z_0_0_xxyzzz[k] - g_yz_0_0_xxyzzz[k] * ab_y + g_yz_0_0_xxyyzzz[k];

                g_yz_0_y_xxzzzz[k] = -g_z_0_0_xxzzzz[k] - g_yz_0_0_xxzzzz[k] * ab_y + g_yz_0_0_xxyzzzz[k];

                g_yz_0_y_xyyyyy[k] = -g_z_0_0_xyyyyy[k] - g_yz_0_0_xyyyyy[k] * ab_y + g_yz_0_0_xyyyyyy[k];

                g_yz_0_y_xyyyyz[k] = -g_z_0_0_xyyyyz[k] - g_yz_0_0_xyyyyz[k] * ab_y + g_yz_0_0_xyyyyyz[k];

                g_yz_0_y_xyyyzz[k] = -g_z_0_0_xyyyzz[k] - g_yz_0_0_xyyyzz[k] * ab_y + g_yz_0_0_xyyyyzz[k];

                g_yz_0_y_xyyzzz[k] = -g_z_0_0_xyyzzz[k] - g_yz_0_0_xyyzzz[k] * ab_y + g_yz_0_0_xyyyzzz[k];

                g_yz_0_y_xyzzzz[k] = -g_z_0_0_xyzzzz[k] - g_yz_0_0_xyzzzz[k] * ab_y + g_yz_0_0_xyyzzzz[k];

                g_yz_0_y_xzzzzz[k] = -g_z_0_0_xzzzzz[k] - g_yz_0_0_xzzzzz[k] * ab_y + g_yz_0_0_xyzzzzz[k];

                g_yz_0_y_yyyyyy[k] = -g_z_0_0_yyyyyy[k] - g_yz_0_0_yyyyyy[k] * ab_y + g_yz_0_0_yyyyyyy[k];

                g_yz_0_y_yyyyyz[k] = -g_z_0_0_yyyyyz[k] - g_yz_0_0_yyyyyz[k] * ab_y + g_yz_0_0_yyyyyyz[k];

                g_yz_0_y_yyyyzz[k] = -g_z_0_0_yyyyzz[k] - g_yz_0_0_yyyyzz[k] * ab_y + g_yz_0_0_yyyyyzz[k];

                g_yz_0_y_yyyzzz[k] = -g_z_0_0_yyyzzz[k] - g_yz_0_0_yyyzzz[k] * ab_y + g_yz_0_0_yyyyzzz[k];

                g_yz_0_y_yyzzzz[k] = -g_z_0_0_yyzzzz[k] - g_yz_0_0_yyzzzz[k] * ab_y + g_yz_0_0_yyyzzzz[k];

                g_yz_0_y_yzzzzz[k] = -g_z_0_0_yzzzzz[k] - g_yz_0_0_yzzzzz[k] * ab_y + g_yz_0_0_yyzzzzz[k];

                g_yz_0_y_zzzzzz[k] = -g_z_0_0_zzzzzz[k] - g_yz_0_0_zzzzzz[k] * ab_y + g_yz_0_0_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_yz_0_z_xxxxxx = cbuffer.data(pi_geom_20_off + 392 * ccomps * dcomps);

            auto g_yz_0_z_xxxxxy = cbuffer.data(pi_geom_20_off + 393 * ccomps * dcomps);

            auto g_yz_0_z_xxxxxz = cbuffer.data(pi_geom_20_off + 394 * ccomps * dcomps);

            auto g_yz_0_z_xxxxyy = cbuffer.data(pi_geom_20_off + 395 * ccomps * dcomps);

            auto g_yz_0_z_xxxxyz = cbuffer.data(pi_geom_20_off + 396 * ccomps * dcomps);

            auto g_yz_0_z_xxxxzz = cbuffer.data(pi_geom_20_off + 397 * ccomps * dcomps);

            auto g_yz_0_z_xxxyyy = cbuffer.data(pi_geom_20_off + 398 * ccomps * dcomps);

            auto g_yz_0_z_xxxyyz = cbuffer.data(pi_geom_20_off + 399 * ccomps * dcomps);

            auto g_yz_0_z_xxxyzz = cbuffer.data(pi_geom_20_off + 400 * ccomps * dcomps);

            auto g_yz_0_z_xxxzzz = cbuffer.data(pi_geom_20_off + 401 * ccomps * dcomps);

            auto g_yz_0_z_xxyyyy = cbuffer.data(pi_geom_20_off + 402 * ccomps * dcomps);

            auto g_yz_0_z_xxyyyz = cbuffer.data(pi_geom_20_off + 403 * ccomps * dcomps);

            auto g_yz_0_z_xxyyzz = cbuffer.data(pi_geom_20_off + 404 * ccomps * dcomps);

            auto g_yz_0_z_xxyzzz = cbuffer.data(pi_geom_20_off + 405 * ccomps * dcomps);

            auto g_yz_0_z_xxzzzz = cbuffer.data(pi_geom_20_off + 406 * ccomps * dcomps);

            auto g_yz_0_z_xyyyyy = cbuffer.data(pi_geom_20_off + 407 * ccomps * dcomps);

            auto g_yz_0_z_xyyyyz = cbuffer.data(pi_geom_20_off + 408 * ccomps * dcomps);

            auto g_yz_0_z_xyyyzz = cbuffer.data(pi_geom_20_off + 409 * ccomps * dcomps);

            auto g_yz_0_z_xyyzzz = cbuffer.data(pi_geom_20_off + 410 * ccomps * dcomps);

            auto g_yz_0_z_xyzzzz = cbuffer.data(pi_geom_20_off + 411 * ccomps * dcomps);

            auto g_yz_0_z_xzzzzz = cbuffer.data(pi_geom_20_off + 412 * ccomps * dcomps);

            auto g_yz_0_z_yyyyyy = cbuffer.data(pi_geom_20_off + 413 * ccomps * dcomps);

            auto g_yz_0_z_yyyyyz = cbuffer.data(pi_geom_20_off + 414 * ccomps * dcomps);

            auto g_yz_0_z_yyyyzz = cbuffer.data(pi_geom_20_off + 415 * ccomps * dcomps);

            auto g_yz_0_z_yyyzzz = cbuffer.data(pi_geom_20_off + 416 * ccomps * dcomps);

            auto g_yz_0_z_yyzzzz = cbuffer.data(pi_geom_20_off + 417 * ccomps * dcomps);

            auto g_yz_0_z_yzzzzz = cbuffer.data(pi_geom_20_off + 418 * ccomps * dcomps);

            auto g_yz_0_z_zzzzzz = cbuffer.data(pi_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxxz, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxzz, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzzz, g_yz_0_0_xxxxxx, g_yz_0_0_xxxxxxz, g_yz_0_0_xxxxxy, g_yz_0_0_xxxxxyz, g_yz_0_0_xxxxxz, g_yz_0_0_xxxxxzz, g_yz_0_0_xxxxyy, g_yz_0_0_xxxxyyz, g_yz_0_0_xxxxyz, g_yz_0_0_xxxxyzz, g_yz_0_0_xxxxzz, g_yz_0_0_xxxxzzz, g_yz_0_0_xxxyyy, g_yz_0_0_xxxyyyz, g_yz_0_0_xxxyyz, g_yz_0_0_xxxyyzz, g_yz_0_0_xxxyzz, g_yz_0_0_xxxyzzz, g_yz_0_0_xxxzzz, g_yz_0_0_xxxzzzz, g_yz_0_0_xxyyyy, g_yz_0_0_xxyyyyz, g_yz_0_0_xxyyyz, g_yz_0_0_xxyyyzz, g_yz_0_0_xxyyzz, g_yz_0_0_xxyyzzz, g_yz_0_0_xxyzzz, g_yz_0_0_xxyzzzz, g_yz_0_0_xxzzzz, g_yz_0_0_xxzzzzz, g_yz_0_0_xyyyyy, g_yz_0_0_xyyyyyz, g_yz_0_0_xyyyyz, g_yz_0_0_xyyyyzz, g_yz_0_0_xyyyzz, g_yz_0_0_xyyyzzz, g_yz_0_0_xyyzzz, g_yz_0_0_xyyzzzz, g_yz_0_0_xyzzzz, g_yz_0_0_xyzzzzz, g_yz_0_0_xzzzzz, g_yz_0_0_xzzzzzz, g_yz_0_0_yyyyyy, g_yz_0_0_yyyyyyz, g_yz_0_0_yyyyyz, g_yz_0_0_yyyyyzz, g_yz_0_0_yyyyzz, g_yz_0_0_yyyyzzz, g_yz_0_0_yyyzzz, g_yz_0_0_yyyzzzz, g_yz_0_0_yyzzzz, g_yz_0_0_yyzzzzz, g_yz_0_0_yzzzzz, g_yz_0_0_yzzzzzz, g_yz_0_0_zzzzzz, g_yz_0_0_zzzzzzz, g_yz_0_z_xxxxxx, g_yz_0_z_xxxxxy, g_yz_0_z_xxxxxz, g_yz_0_z_xxxxyy, g_yz_0_z_xxxxyz, g_yz_0_z_xxxxzz, g_yz_0_z_xxxyyy, g_yz_0_z_xxxyyz, g_yz_0_z_xxxyzz, g_yz_0_z_xxxzzz, g_yz_0_z_xxyyyy, g_yz_0_z_xxyyyz, g_yz_0_z_xxyyzz, g_yz_0_z_xxyzzz, g_yz_0_z_xxzzzz, g_yz_0_z_xyyyyy, g_yz_0_z_xyyyyz, g_yz_0_z_xyyyzz, g_yz_0_z_xyyzzz, g_yz_0_z_xyzzzz, g_yz_0_z_xzzzzz, g_yz_0_z_yyyyyy, g_yz_0_z_yyyyyz, g_yz_0_z_yyyyzz, g_yz_0_z_yyyzzz, g_yz_0_z_yyzzzz, g_yz_0_z_yzzzzz, g_yz_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_z_xxxxxx[k] = -g_y_0_0_xxxxxx[k] - g_yz_0_0_xxxxxx[k] * ab_z + g_yz_0_0_xxxxxxz[k];

                g_yz_0_z_xxxxxy[k] = -g_y_0_0_xxxxxy[k] - g_yz_0_0_xxxxxy[k] * ab_z + g_yz_0_0_xxxxxyz[k];

                g_yz_0_z_xxxxxz[k] = -g_y_0_0_xxxxxz[k] - g_yz_0_0_xxxxxz[k] * ab_z + g_yz_0_0_xxxxxzz[k];

                g_yz_0_z_xxxxyy[k] = -g_y_0_0_xxxxyy[k] - g_yz_0_0_xxxxyy[k] * ab_z + g_yz_0_0_xxxxyyz[k];

                g_yz_0_z_xxxxyz[k] = -g_y_0_0_xxxxyz[k] - g_yz_0_0_xxxxyz[k] * ab_z + g_yz_0_0_xxxxyzz[k];

                g_yz_0_z_xxxxzz[k] = -g_y_0_0_xxxxzz[k] - g_yz_0_0_xxxxzz[k] * ab_z + g_yz_0_0_xxxxzzz[k];

                g_yz_0_z_xxxyyy[k] = -g_y_0_0_xxxyyy[k] - g_yz_0_0_xxxyyy[k] * ab_z + g_yz_0_0_xxxyyyz[k];

                g_yz_0_z_xxxyyz[k] = -g_y_0_0_xxxyyz[k] - g_yz_0_0_xxxyyz[k] * ab_z + g_yz_0_0_xxxyyzz[k];

                g_yz_0_z_xxxyzz[k] = -g_y_0_0_xxxyzz[k] - g_yz_0_0_xxxyzz[k] * ab_z + g_yz_0_0_xxxyzzz[k];

                g_yz_0_z_xxxzzz[k] = -g_y_0_0_xxxzzz[k] - g_yz_0_0_xxxzzz[k] * ab_z + g_yz_0_0_xxxzzzz[k];

                g_yz_0_z_xxyyyy[k] = -g_y_0_0_xxyyyy[k] - g_yz_0_0_xxyyyy[k] * ab_z + g_yz_0_0_xxyyyyz[k];

                g_yz_0_z_xxyyyz[k] = -g_y_0_0_xxyyyz[k] - g_yz_0_0_xxyyyz[k] * ab_z + g_yz_0_0_xxyyyzz[k];

                g_yz_0_z_xxyyzz[k] = -g_y_0_0_xxyyzz[k] - g_yz_0_0_xxyyzz[k] * ab_z + g_yz_0_0_xxyyzzz[k];

                g_yz_0_z_xxyzzz[k] = -g_y_0_0_xxyzzz[k] - g_yz_0_0_xxyzzz[k] * ab_z + g_yz_0_0_xxyzzzz[k];

                g_yz_0_z_xxzzzz[k] = -g_y_0_0_xxzzzz[k] - g_yz_0_0_xxzzzz[k] * ab_z + g_yz_0_0_xxzzzzz[k];

                g_yz_0_z_xyyyyy[k] = -g_y_0_0_xyyyyy[k] - g_yz_0_0_xyyyyy[k] * ab_z + g_yz_0_0_xyyyyyz[k];

                g_yz_0_z_xyyyyz[k] = -g_y_0_0_xyyyyz[k] - g_yz_0_0_xyyyyz[k] * ab_z + g_yz_0_0_xyyyyzz[k];

                g_yz_0_z_xyyyzz[k] = -g_y_0_0_xyyyzz[k] - g_yz_0_0_xyyyzz[k] * ab_z + g_yz_0_0_xyyyzzz[k];

                g_yz_0_z_xyyzzz[k] = -g_y_0_0_xyyzzz[k] - g_yz_0_0_xyyzzz[k] * ab_z + g_yz_0_0_xyyzzzz[k];

                g_yz_0_z_xyzzzz[k] = -g_y_0_0_xyzzzz[k] - g_yz_0_0_xyzzzz[k] * ab_z + g_yz_0_0_xyzzzzz[k];

                g_yz_0_z_xzzzzz[k] = -g_y_0_0_xzzzzz[k] - g_yz_0_0_xzzzzz[k] * ab_z + g_yz_0_0_xzzzzzz[k];

                g_yz_0_z_yyyyyy[k] = -g_y_0_0_yyyyyy[k] - g_yz_0_0_yyyyyy[k] * ab_z + g_yz_0_0_yyyyyyz[k];

                g_yz_0_z_yyyyyz[k] = -g_y_0_0_yyyyyz[k] - g_yz_0_0_yyyyyz[k] * ab_z + g_yz_0_0_yyyyyzz[k];

                g_yz_0_z_yyyyzz[k] = -g_y_0_0_yyyyzz[k] - g_yz_0_0_yyyyzz[k] * ab_z + g_yz_0_0_yyyyzzz[k];

                g_yz_0_z_yyyzzz[k] = -g_y_0_0_yyyzzz[k] - g_yz_0_0_yyyzzz[k] * ab_z + g_yz_0_0_yyyzzzz[k];

                g_yz_0_z_yyzzzz[k] = -g_y_0_0_yyzzzz[k] - g_yz_0_0_yyzzzz[k] * ab_z + g_yz_0_0_yyzzzzz[k];

                g_yz_0_z_yzzzzz[k] = -g_y_0_0_yzzzzz[k] - g_yz_0_0_yzzzzz[k] * ab_z + g_yz_0_0_yzzzzzz[k];

                g_yz_0_z_zzzzzz[k] = -g_y_0_0_zzzzzz[k] - g_yz_0_0_zzzzzz[k] * ab_z + g_yz_0_0_zzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_zz_0_x_xxxxxx = cbuffer.data(pi_geom_20_off + 420 * ccomps * dcomps);

            auto g_zz_0_x_xxxxxy = cbuffer.data(pi_geom_20_off + 421 * ccomps * dcomps);

            auto g_zz_0_x_xxxxxz = cbuffer.data(pi_geom_20_off + 422 * ccomps * dcomps);

            auto g_zz_0_x_xxxxyy = cbuffer.data(pi_geom_20_off + 423 * ccomps * dcomps);

            auto g_zz_0_x_xxxxyz = cbuffer.data(pi_geom_20_off + 424 * ccomps * dcomps);

            auto g_zz_0_x_xxxxzz = cbuffer.data(pi_geom_20_off + 425 * ccomps * dcomps);

            auto g_zz_0_x_xxxyyy = cbuffer.data(pi_geom_20_off + 426 * ccomps * dcomps);

            auto g_zz_0_x_xxxyyz = cbuffer.data(pi_geom_20_off + 427 * ccomps * dcomps);

            auto g_zz_0_x_xxxyzz = cbuffer.data(pi_geom_20_off + 428 * ccomps * dcomps);

            auto g_zz_0_x_xxxzzz = cbuffer.data(pi_geom_20_off + 429 * ccomps * dcomps);

            auto g_zz_0_x_xxyyyy = cbuffer.data(pi_geom_20_off + 430 * ccomps * dcomps);

            auto g_zz_0_x_xxyyyz = cbuffer.data(pi_geom_20_off + 431 * ccomps * dcomps);

            auto g_zz_0_x_xxyyzz = cbuffer.data(pi_geom_20_off + 432 * ccomps * dcomps);

            auto g_zz_0_x_xxyzzz = cbuffer.data(pi_geom_20_off + 433 * ccomps * dcomps);

            auto g_zz_0_x_xxzzzz = cbuffer.data(pi_geom_20_off + 434 * ccomps * dcomps);

            auto g_zz_0_x_xyyyyy = cbuffer.data(pi_geom_20_off + 435 * ccomps * dcomps);

            auto g_zz_0_x_xyyyyz = cbuffer.data(pi_geom_20_off + 436 * ccomps * dcomps);

            auto g_zz_0_x_xyyyzz = cbuffer.data(pi_geom_20_off + 437 * ccomps * dcomps);

            auto g_zz_0_x_xyyzzz = cbuffer.data(pi_geom_20_off + 438 * ccomps * dcomps);

            auto g_zz_0_x_xyzzzz = cbuffer.data(pi_geom_20_off + 439 * ccomps * dcomps);

            auto g_zz_0_x_xzzzzz = cbuffer.data(pi_geom_20_off + 440 * ccomps * dcomps);

            auto g_zz_0_x_yyyyyy = cbuffer.data(pi_geom_20_off + 441 * ccomps * dcomps);

            auto g_zz_0_x_yyyyyz = cbuffer.data(pi_geom_20_off + 442 * ccomps * dcomps);

            auto g_zz_0_x_yyyyzz = cbuffer.data(pi_geom_20_off + 443 * ccomps * dcomps);

            auto g_zz_0_x_yyyzzz = cbuffer.data(pi_geom_20_off + 444 * ccomps * dcomps);

            auto g_zz_0_x_yyzzzz = cbuffer.data(pi_geom_20_off + 445 * ccomps * dcomps);

            auto g_zz_0_x_yzzzzz = cbuffer.data(pi_geom_20_off + 446 * ccomps * dcomps);

            auto g_zz_0_x_zzzzzz = cbuffer.data(pi_geom_20_off + 447 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_0_xxxxxx, g_zz_0_0_xxxxxxx, g_zz_0_0_xxxxxxy, g_zz_0_0_xxxxxxz, g_zz_0_0_xxxxxy, g_zz_0_0_xxxxxyy, g_zz_0_0_xxxxxyz, g_zz_0_0_xxxxxz, g_zz_0_0_xxxxxzz, g_zz_0_0_xxxxyy, g_zz_0_0_xxxxyyy, g_zz_0_0_xxxxyyz, g_zz_0_0_xxxxyz, g_zz_0_0_xxxxyzz, g_zz_0_0_xxxxzz, g_zz_0_0_xxxxzzz, g_zz_0_0_xxxyyy, g_zz_0_0_xxxyyyy, g_zz_0_0_xxxyyyz, g_zz_0_0_xxxyyz, g_zz_0_0_xxxyyzz, g_zz_0_0_xxxyzz, g_zz_0_0_xxxyzzz, g_zz_0_0_xxxzzz, g_zz_0_0_xxxzzzz, g_zz_0_0_xxyyyy, g_zz_0_0_xxyyyyy, g_zz_0_0_xxyyyyz, g_zz_0_0_xxyyyz, g_zz_0_0_xxyyyzz, g_zz_0_0_xxyyzz, g_zz_0_0_xxyyzzz, g_zz_0_0_xxyzzz, g_zz_0_0_xxyzzzz, g_zz_0_0_xxzzzz, g_zz_0_0_xxzzzzz, g_zz_0_0_xyyyyy, g_zz_0_0_xyyyyyy, g_zz_0_0_xyyyyyz, g_zz_0_0_xyyyyz, g_zz_0_0_xyyyyzz, g_zz_0_0_xyyyzz, g_zz_0_0_xyyyzzz, g_zz_0_0_xyyzzz, g_zz_0_0_xyyzzzz, g_zz_0_0_xyzzzz, g_zz_0_0_xyzzzzz, g_zz_0_0_xzzzzz, g_zz_0_0_xzzzzzz, g_zz_0_0_yyyyyy, g_zz_0_0_yyyyyz, g_zz_0_0_yyyyzz, g_zz_0_0_yyyzzz, g_zz_0_0_yyzzzz, g_zz_0_0_yzzzzz, g_zz_0_0_zzzzzz, g_zz_0_x_xxxxxx, g_zz_0_x_xxxxxy, g_zz_0_x_xxxxxz, g_zz_0_x_xxxxyy, g_zz_0_x_xxxxyz, g_zz_0_x_xxxxzz, g_zz_0_x_xxxyyy, g_zz_0_x_xxxyyz, g_zz_0_x_xxxyzz, g_zz_0_x_xxxzzz, g_zz_0_x_xxyyyy, g_zz_0_x_xxyyyz, g_zz_0_x_xxyyzz, g_zz_0_x_xxyzzz, g_zz_0_x_xxzzzz, g_zz_0_x_xyyyyy, g_zz_0_x_xyyyyz, g_zz_0_x_xyyyzz, g_zz_0_x_xyyzzz, g_zz_0_x_xyzzzz, g_zz_0_x_xzzzzz, g_zz_0_x_yyyyyy, g_zz_0_x_yyyyyz, g_zz_0_x_yyyyzz, g_zz_0_x_yyyzzz, g_zz_0_x_yyzzzz, g_zz_0_x_yzzzzz, g_zz_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_x_xxxxxx[k] = -g_zz_0_0_xxxxxx[k] * ab_x + g_zz_0_0_xxxxxxx[k];

                g_zz_0_x_xxxxxy[k] = -g_zz_0_0_xxxxxy[k] * ab_x + g_zz_0_0_xxxxxxy[k];

                g_zz_0_x_xxxxxz[k] = -g_zz_0_0_xxxxxz[k] * ab_x + g_zz_0_0_xxxxxxz[k];

                g_zz_0_x_xxxxyy[k] = -g_zz_0_0_xxxxyy[k] * ab_x + g_zz_0_0_xxxxxyy[k];

                g_zz_0_x_xxxxyz[k] = -g_zz_0_0_xxxxyz[k] * ab_x + g_zz_0_0_xxxxxyz[k];

                g_zz_0_x_xxxxzz[k] = -g_zz_0_0_xxxxzz[k] * ab_x + g_zz_0_0_xxxxxzz[k];

                g_zz_0_x_xxxyyy[k] = -g_zz_0_0_xxxyyy[k] * ab_x + g_zz_0_0_xxxxyyy[k];

                g_zz_0_x_xxxyyz[k] = -g_zz_0_0_xxxyyz[k] * ab_x + g_zz_0_0_xxxxyyz[k];

                g_zz_0_x_xxxyzz[k] = -g_zz_0_0_xxxyzz[k] * ab_x + g_zz_0_0_xxxxyzz[k];

                g_zz_0_x_xxxzzz[k] = -g_zz_0_0_xxxzzz[k] * ab_x + g_zz_0_0_xxxxzzz[k];

                g_zz_0_x_xxyyyy[k] = -g_zz_0_0_xxyyyy[k] * ab_x + g_zz_0_0_xxxyyyy[k];

                g_zz_0_x_xxyyyz[k] = -g_zz_0_0_xxyyyz[k] * ab_x + g_zz_0_0_xxxyyyz[k];

                g_zz_0_x_xxyyzz[k] = -g_zz_0_0_xxyyzz[k] * ab_x + g_zz_0_0_xxxyyzz[k];

                g_zz_0_x_xxyzzz[k] = -g_zz_0_0_xxyzzz[k] * ab_x + g_zz_0_0_xxxyzzz[k];

                g_zz_0_x_xxzzzz[k] = -g_zz_0_0_xxzzzz[k] * ab_x + g_zz_0_0_xxxzzzz[k];

                g_zz_0_x_xyyyyy[k] = -g_zz_0_0_xyyyyy[k] * ab_x + g_zz_0_0_xxyyyyy[k];

                g_zz_0_x_xyyyyz[k] = -g_zz_0_0_xyyyyz[k] * ab_x + g_zz_0_0_xxyyyyz[k];

                g_zz_0_x_xyyyzz[k] = -g_zz_0_0_xyyyzz[k] * ab_x + g_zz_0_0_xxyyyzz[k];

                g_zz_0_x_xyyzzz[k] = -g_zz_0_0_xyyzzz[k] * ab_x + g_zz_0_0_xxyyzzz[k];

                g_zz_0_x_xyzzzz[k] = -g_zz_0_0_xyzzzz[k] * ab_x + g_zz_0_0_xxyzzzz[k];

                g_zz_0_x_xzzzzz[k] = -g_zz_0_0_xzzzzz[k] * ab_x + g_zz_0_0_xxzzzzz[k];

                g_zz_0_x_yyyyyy[k] = -g_zz_0_0_yyyyyy[k] * ab_x + g_zz_0_0_xyyyyyy[k];

                g_zz_0_x_yyyyyz[k] = -g_zz_0_0_yyyyyz[k] * ab_x + g_zz_0_0_xyyyyyz[k];

                g_zz_0_x_yyyyzz[k] = -g_zz_0_0_yyyyzz[k] * ab_x + g_zz_0_0_xyyyyzz[k];

                g_zz_0_x_yyyzzz[k] = -g_zz_0_0_yyyzzz[k] * ab_x + g_zz_0_0_xyyyzzz[k];

                g_zz_0_x_yyzzzz[k] = -g_zz_0_0_yyzzzz[k] * ab_x + g_zz_0_0_xyyzzzz[k];

                g_zz_0_x_yzzzzz[k] = -g_zz_0_0_yzzzzz[k] * ab_x + g_zz_0_0_xyzzzzz[k];

                g_zz_0_x_zzzzzz[k] = -g_zz_0_0_zzzzzz[k] * ab_x + g_zz_0_0_xzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_zz_0_y_xxxxxx = cbuffer.data(pi_geom_20_off + 448 * ccomps * dcomps);

            auto g_zz_0_y_xxxxxy = cbuffer.data(pi_geom_20_off + 449 * ccomps * dcomps);

            auto g_zz_0_y_xxxxxz = cbuffer.data(pi_geom_20_off + 450 * ccomps * dcomps);

            auto g_zz_0_y_xxxxyy = cbuffer.data(pi_geom_20_off + 451 * ccomps * dcomps);

            auto g_zz_0_y_xxxxyz = cbuffer.data(pi_geom_20_off + 452 * ccomps * dcomps);

            auto g_zz_0_y_xxxxzz = cbuffer.data(pi_geom_20_off + 453 * ccomps * dcomps);

            auto g_zz_0_y_xxxyyy = cbuffer.data(pi_geom_20_off + 454 * ccomps * dcomps);

            auto g_zz_0_y_xxxyyz = cbuffer.data(pi_geom_20_off + 455 * ccomps * dcomps);

            auto g_zz_0_y_xxxyzz = cbuffer.data(pi_geom_20_off + 456 * ccomps * dcomps);

            auto g_zz_0_y_xxxzzz = cbuffer.data(pi_geom_20_off + 457 * ccomps * dcomps);

            auto g_zz_0_y_xxyyyy = cbuffer.data(pi_geom_20_off + 458 * ccomps * dcomps);

            auto g_zz_0_y_xxyyyz = cbuffer.data(pi_geom_20_off + 459 * ccomps * dcomps);

            auto g_zz_0_y_xxyyzz = cbuffer.data(pi_geom_20_off + 460 * ccomps * dcomps);

            auto g_zz_0_y_xxyzzz = cbuffer.data(pi_geom_20_off + 461 * ccomps * dcomps);

            auto g_zz_0_y_xxzzzz = cbuffer.data(pi_geom_20_off + 462 * ccomps * dcomps);

            auto g_zz_0_y_xyyyyy = cbuffer.data(pi_geom_20_off + 463 * ccomps * dcomps);

            auto g_zz_0_y_xyyyyz = cbuffer.data(pi_geom_20_off + 464 * ccomps * dcomps);

            auto g_zz_0_y_xyyyzz = cbuffer.data(pi_geom_20_off + 465 * ccomps * dcomps);

            auto g_zz_0_y_xyyzzz = cbuffer.data(pi_geom_20_off + 466 * ccomps * dcomps);

            auto g_zz_0_y_xyzzzz = cbuffer.data(pi_geom_20_off + 467 * ccomps * dcomps);

            auto g_zz_0_y_xzzzzz = cbuffer.data(pi_geom_20_off + 468 * ccomps * dcomps);

            auto g_zz_0_y_yyyyyy = cbuffer.data(pi_geom_20_off + 469 * ccomps * dcomps);

            auto g_zz_0_y_yyyyyz = cbuffer.data(pi_geom_20_off + 470 * ccomps * dcomps);

            auto g_zz_0_y_yyyyzz = cbuffer.data(pi_geom_20_off + 471 * ccomps * dcomps);

            auto g_zz_0_y_yyyzzz = cbuffer.data(pi_geom_20_off + 472 * ccomps * dcomps);

            auto g_zz_0_y_yyzzzz = cbuffer.data(pi_geom_20_off + 473 * ccomps * dcomps);

            auto g_zz_0_y_yzzzzz = cbuffer.data(pi_geom_20_off + 474 * ccomps * dcomps);

            auto g_zz_0_y_zzzzzz = cbuffer.data(pi_geom_20_off + 475 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_0_xxxxxx, g_zz_0_0_xxxxxxy, g_zz_0_0_xxxxxy, g_zz_0_0_xxxxxyy, g_zz_0_0_xxxxxyz, g_zz_0_0_xxxxxz, g_zz_0_0_xxxxyy, g_zz_0_0_xxxxyyy, g_zz_0_0_xxxxyyz, g_zz_0_0_xxxxyz, g_zz_0_0_xxxxyzz, g_zz_0_0_xxxxzz, g_zz_0_0_xxxyyy, g_zz_0_0_xxxyyyy, g_zz_0_0_xxxyyyz, g_zz_0_0_xxxyyz, g_zz_0_0_xxxyyzz, g_zz_0_0_xxxyzz, g_zz_0_0_xxxyzzz, g_zz_0_0_xxxzzz, g_zz_0_0_xxyyyy, g_zz_0_0_xxyyyyy, g_zz_0_0_xxyyyyz, g_zz_0_0_xxyyyz, g_zz_0_0_xxyyyzz, g_zz_0_0_xxyyzz, g_zz_0_0_xxyyzzz, g_zz_0_0_xxyzzz, g_zz_0_0_xxyzzzz, g_zz_0_0_xxzzzz, g_zz_0_0_xyyyyy, g_zz_0_0_xyyyyyy, g_zz_0_0_xyyyyyz, g_zz_0_0_xyyyyz, g_zz_0_0_xyyyyzz, g_zz_0_0_xyyyzz, g_zz_0_0_xyyyzzz, g_zz_0_0_xyyzzz, g_zz_0_0_xyyzzzz, g_zz_0_0_xyzzzz, g_zz_0_0_xyzzzzz, g_zz_0_0_xzzzzz, g_zz_0_0_yyyyyy, g_zz_0_0_yyyyyyy, g_zz_0_0_yyyyyyz, g_zz_0_0_yyyyyz, g_zz_0_0_yyyyyzz, g_zz_0_0_yyyyzz, g_zz_0_0_yyyyzzz, g_zz_0_0_yyyzzz, g_zz_0_0_yyyzzzz, g_zz_0_0_yyzzzz, g_zz_0_0_yyzzzzz, g_zz_0_0_yzzzzz, g_zz_0_0_yzzzzzz, g_zz_0_0_zzzzzz, g_zz_0_y_xxxxxx, g_zz_0_y_xxxxxy, g_zz_0_y_xxxxxz, g_zz_0_y_xxxxyy, g_zz_0_y_xxxxyz, g_zz_0_y_xxxxzz, g_zz_0_y_xxxyyy, g_zz_0_y_xxxyyz, g_zz_0_y_xxxyzz, g_zz_0_y_xxxzzz, g_zz_0_y_xxyyyy, g_zz_0_y_xxyyyz, g_zz_0_y_xxyyzz, g_zz_0_y_xxyzzz, g_zz_0_y_xxzzzz, g_zz_0_y_xyyyyy, g_zz_0_y_xyyyyz, g_zz_0_y_xyyyzz, g_zz_0_y_xyyzzz, g_zz_0_y_xyzzzz, g_zz_0_y_xzzzzz, g_zz_0_y_yyyyyy, g_zz_0_y_yyyyyz, g_zz_0_y_yyyyzz, g_zz_0_y_yyyzzz, g_zz_0_y_yyzzzz, g_zz_0_y_yzzzzz, g_zz_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_y_xxxxxx[k] = -g_zz_0_0_xxxxxx[k] * ab_y + g_zz_0_0_xxxxxxy[k];

                g_zz_0_y_xxxxxy[k] = -g_zz_0_0_xxxxxy[k] * ab_y + g_zz_0_0_xxxxxyy[k];

                g_zz_0_y_xxxxxz[k] = -g_zz_0_0_xxxxxz[k] * ab_y + g_zz_0_0_xxxxxyz[k];

                g_zz_0_y_xxxxyy[k] = -g_zz_0_0_xxxxyy[k] * ab_y + g_zz_0_0_xxxxyyy[k];

                g_zz_0_y_xxxxyz[k] = -g_zz_0_0_xxxxyz[k] * ab_y + g_zz_0_0_xxxxyyz[k];

                g_zz_0_y_xxxxzz[k] = -g_zz_0_0_xxxxzz[k] * ab_y + g_zz_0_0_xxxxyzz[k];

                g_zz_0_y_xxxyyy[k] = -g_zz_0_0_xxxyyy[k] * ab_y + g_zz_0_0_xxxyyyy[k];

                g_zz_0_y_xxxyyz[k] = -g_zz_0_0_xxxyyz[k] * ab_y + g_zz_0_0_xxxyyyz[k];

                g_zz_0_y_xxxyzz[k] = -g_zz_0_0_xxxyzz[k] * ab_y + g_zz_0_0_xxxyyzz[k];

                g_zz_0_y_xxxzzz[k] = -g_zz_0_0_xxxzzz[k] * ab_y + g_zz_0_0_xxxyzzz[k];

                g_zz_0_y_xxyyyy[k] = -g_zz_0_0_xxyyyy[k] * ab_y + g_zz_0_0_xxyyyyy[k];

                g_zz_0_y_xxyyyz[k] = -g_zz_0_0_xxyyyz[k] * ab_y + g_zz_0_0_xxyyyyz[k];

                g_zz_0_y_xxyyzz[k] = -g_zz_0_0_xxyyzz[k] * ab_y + g_zz_0_0_xxyyyzz[k];

                g_zz_0_y_xxyzzz[k] = -g_zz_0_0_xxyzzz[k] * ab_y + g_zz_0_0_xxyyzzz[k];

                g_zz_0_y_xxzzzz[k] = -g_zz_0_0_xxzzzz[k] * ab_y + g_zz_0_0_xxyzzzz[k];

                g_zz_0_y_xyyyyy[k] = -g_zz_0_0_xyyyyy[k] * ab_y + g_zz_0_0_xyyyyyy[k];

                g_zz_0_y_xyyyyz[k] = -g_zz_0_0_xyyyyz[k] * ab_y + g_zz_0_0_xyyyyyz[k];

                g_zz_0_y_xyyyzz[k] = -g_zz_0_0_xyyyzz[k] * ab_y + g_zz_0_0_xyyyyzz[k];

                g_zz_0_y_xyyzzz[k] = -g_zz_0_0_xyyzzz[k] * ab_y + g_zz_0_0_xyyyzzz[k];

                g_zz_0_y_xyzzzz[k] = -g_zz_0_0_xyzzzz[k] * ab_y + g_zz_0_0_xyyzzzz[k];

                g_zz_0_y_xzzzzz[k] = -g_zz_0_0_xzzzzz[k] * ab_y + g_zz_0_0_xyzzzzz[k];

                g_zz_0_y_yyyyyy[k] = -g_zz_0_0_yyyyyy[k] * ab_y + g_zz_0_0_yyyyyyy[k];

                g_zz_0_y_yyyyyz[k] = -g_zz_0_0_yyyyyz[k] * ab_y + g_zz_0_0_yyyyyyz[k];

                g_zz_0_y_yyyyzz[k] = -g_zz_0_0_yyyyzz[k] * ab_y + g_zz_0_0_yyyyyzz[k];

                g_zz_0_y_yyyzzz[k] = -g_zz_0_0_yyyzzz[k] * ab_y + g_zz_0_0_yyyyzzz[k];

                g_zz_0_y_yyzzzz[k] = -g_zz_0_0_yyzzzz[k] * ab_y + g_zz_0_0_yyyzzzz[k];

                g_zz_0_y_yzzzzz[k] = -g_zz_0_0_yzzzzz[k] * ab_y + g_zz_0_0_yyzzzzz[k];

                g_zz_0_y_zzzzzz[k] = -g_zz_0_0_zzzzzz[k] * ab_y + g_zz_0_0_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_zz_0_z_xxxxxx = cbuffer.data(pi_geom_20_off + 476 * ccomps * dcomps);

            auto g_zz_0_z_xxxxxy = cbuffer.data(pi_geom_20_off + 477 * ccomps * dcomps);

            auto g_zz_0_z_xxxxxz = cbuffer.data(pi_geom_20_off + 478 * ccomps * dcomps);

            auto g_zz_0_z_xxxxyy = cbuffer.data(pi_geom_20_off + 479 * ccomps * dcomps);

            auto g_zz_0_z_xxxxyz = cbuffer.data(pi_geom_20_off + 480 * ccomps * dcomps);

            auto g_zz_0_z_xxxxzz = cbuffer.data(pi_geom_20_off + 481 * ccomps * dcomps);

            auto g_zz_0_z_xxxyyy = cbuffer.data(pi_geom_20_off + 482 * ccomps * dcomps);

            auto g_zz_0_z_xxxyyz = cbuffer.data(pi_geom_20_off + 483 * ccomps * dcomps);

            auto g_zz_0_z_xxxyzz = cbuffer.data(pi_geom_20_off + 484 * ccomps * dcomps);

            auto g_zz_0_z_xxxzzz = cbuffer.data(pi_geom_20_off + 485 * ccomps * dcomps);

            auto g_zz_0_z_xxyyyy = cbuffer.data(pi_geom_20_off + 486 * ccomps * dcomps);

            auto g_zz_0_z_xxyyyz = cbuffer.data(pi_geom_20_off + 487 * ccomps * dcomps);

            auto g_zz_0_z_xxyyzz = cbuffer.data(pi_geom_20_off + 488 * ccomps * dcomps);

            auto g_zz_0_z_xxyzzz = cbuffer.data(pi_geom_20_off + 489 * ccomps * dcomps);

            auto g_zz_0_z_xxzzzz = cbuffer.data(pi_geom_20_off + 490 * ccomps * dcomps);

            auto g_zz_0_z_xyyyyy = cbuffer.data(pi_geom_20_off + 491 * ccomps * dcomps);

            auto g_zz_0_z_xyyyyz = cbuffer.data(pi_geom_20_off + 492 * ccomps * dcomps);

            auto g_zz_0_z_xyyyzz = cbuffer.data(pi_geom_20_off + 493 * ccomps * dcomps);

            auto g_zz_0_z_xyyzzz = cbuffer.data(pi_geom_20_off + 494 * ccomps * dcomps);

            auto g_zz_0_z_xyzzzz = cbuffer.data(pi_geom_20_off + 495 * ccomps * dcomps);

            auto g_zz_0_z_xzzzzz = cbuffer.data(pi_geom_20_off + 496 * ccomps * dcomps);

            auto g_zz_0_z_yyyyyy = cbuffer.data(pi_geom_20_off + 497 * ccomps * dcomps);

            auto g_zz_0_z_yyyyyz = cbuffer.data(pi_geom_20_off + 498 * ccomps * dcomps);

            auto g_zz_0_z_yyyyzz = cbuffer.data(pi_geom_20_off + 499 * ccomps * dcomps);

            auto g_zz_0_z_yyyzzz = cbuffer.data(pi_geom_20_off + 500 * ccomps * dcomps);

            auto g_zz_0_z_yyzzzz = cbuffer.data(pi_geom_20_off + 501 * ccomps * dcomps);

            auto g_zz_0_z_yzzzzz = cbuffer.data(pi_geom_20_off + 502 * ccomps * dcomps);

            auto g_zz_0_z_zzzzzz = cbuffer.data(pi_geom_20_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxxz, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxzz, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzzz, g_zz_0_0_xxxxxx, g_zz_0_0_xxxxxxz, g_zz_0_0_xxxxxy, g_zz_0_0_xxxxxyz, g_zz_0_0_xxxxxz, g_zz_0_0_xxxxxzz, g_zz_0_0_xxxxyy, g_zz_0_0_xxxxyyz, g_zz_0_0_xxxxyz, g_zz_0_0_xxxxyzz, g_zz_0_0_xxxxzz, g_zz_0_0_xxxxzzz, g_zz_0_0_xxxyyy, g_zz_0_0_xxxyyyz, g_zz_0_0_xxxyyz, g_zz_0_0_xxxyyzz, g_zz_0_0_xxxyzz, g_zz_0_0_xxxyzzz, g_zz_0_0_xxxzzz, g_zz_0_0_xxxzzzz, g_zz_0_0_xxyyyy, g_zz_0_0_xxyyyyz, g_zz_0_0_xxyyyz, g_zz_0_0_xxyyyzz, g_zz_0_0_xxyyzz, g_zz_0_0_xxyyzzz, g_zz_0_0_xxyzzz, g_zz_0_0_xxyzzzz, g_zz_0_0_xxzzzz, g_zz_0_0_xxzzzzz, g_zz_0_0_xyyyyy, g_zz_0_0_xyyyyyz, g_zz_0_0_xyyyyz, g_zz_0_0_xyyyyzz, g_zz_0_0_xyyyzz, g_zz_0_0_xyyyzzz, g_zz_0_0_xyyzzz, g_zz_0_0_xyyzzzz, g_zz_0_0_xyzzzz, g_zz_0_0_xyzzzzz, g_zz_0_0_xzzzzz, g_zz_0_0_xzzzzzz, g_zz_0_0_yyyyyy, g_zz_0_0_yyyyyyz, g_zz_0_0_yyyyyz, g_zz_0_0_yyyyyzz, g_zz_0_0_yyyyzz, g_zz_0_0_yyyyzzz, g_zz_0_0_yyyzzz, g_zz_0_0_yyyzzzz, g_zz_0_0_yyzzzz, g_zz_0_0_yyzzzzz, g_zz_0_0_yzzzzz, g_zz_0_0_yzzzzzz, g_zz_0_0_zzzzzz, g_zz_0_0_zzzzzzz, g_zz_0_z_xxxxxx, g_zz_0_z_xxxxxy, g_zz_0_z_xxxxxz, g_zz_0_z_xxxxyy, g_zz_0_z_xxxxyz, g_zz_0_z_xxxxzz, g_zz_0_z_xxxyyy, g_zz_0_z_xxxyyz, g_zz_0_z_xxxyzz, g_zz_0_z_xxxzzz, g_zz_0_z_xxyyyy, g_zz_0_z_xxyyyz, g_zz_0_z_xxyyzz, g_zz_0_z_xxyzzz, g_zz_0_z_xxzzzz, g_zz_0_z_xyyyyy, g_zz_0_z_xyyyyz, g_zz_0_z_xyyyzz, g_zz_0_z_xyyzzz, g_zz_0_z_xyzzzz, g_zz_0_z_xzzzzz, g_zz_0_z_yyyyyy, g_zz_0_z_yyyyyz, g_zz_0_z_yyyyzz, g_zz_0_z_yyyzzz, g_zz_0_z_yyzzzz, g_zz_0_z_yzzzzz, g_zz_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_z_xxxxxx[k] = -2.0 * g_z_0_0_xxxxxx[k] - g_zz_0_0_xxxxxx[k] * ab_z + g_zz_0_0_xxxxxxz[k];

                g_zz_0_z_xxxxxy[k] = -2.0 * g_z_0_0_xxxxxy[k] - g_zz_0_0_xxxxxy[k] * ab_z + g_zz_0_0_xxxxxyz[k];

                g_zz_0_z_xxxxxz[k] = -2.0 * g_z_0_0_xxxxxz[k] - g_zz_0_0_xxxxxz[k] * ab_z + g_zz_0_0_xxxxxzz[k];

                g_zz_0_z_xxxxyy[k] = -2.0 * g_z_0_0_xxxxyy[k] - g_zz_0_0_xxxxyy[k] * ab_z + g_zz_0_0_xxxxyyz[k];

                g_zz_0_z_xxxxyz[k] = -2.0 * g_z_0_0_xxxxyz[k] - g_zz_0_0_xxxxyz[k] * ab_z + g_zz_0_0_xxxxyzz[k];

                g_zz_0_z_xxxxzz[k] = -2.0 * g_z_0_0_xxxxzz[k] - g_zz_0_0_xxxxzz[k] * ab_z + g_zz_0_0_xxxxzzz[k];

                g_zz_0_z_xxxyyy[k] = -2.0 * g_z_0_0_xxxyyy[k] - g_zz_0_0_xxxyyy[k] * ab_z + g_zz_0_0_xxxyyyz[k];

                g_zz_0_z_xxxyyz[k] = -2.0 * g_z_0_0_xxxyyz[k] - g_zz_0_0_xxxyyz[k] * ab_z + g_zz_0_0_xxxyyzz[k];

                g_zz_0_z_xxxyzz[k] = -2.0 * g_z_0_0_xxxyzz[k] - g_zz_0_0_xxxyzz[k] * ab_z + g_zz_0_0_xxxyzzz[k];

                g_zz_0_z_xxxzzz[k] = -2.0 * g_z_0_0_xxxzzz[k] - g_zz_0_0_xxxzzz[k] * ab_z + g_zz_0_0_xxxzzzz[k];

                g_zz_0_z_xxyyyy[k] = -2.0 * g_z_0_0_xxyyyy[k] - g_zz_0_0_xxyyyy[k] * ab_z + g_zz_0_0_xxyyyyz[k];

                g_zz_0_z_xxyyyz[k] = -2.0 * g_z_0_0_xxyyyz[k] - g_zz_0_0_xxyyyz[k] * ab_z + g_zz_0_0_xxyyyzz[k];

                g_zz_0_z_xxyyzz[k] = -2.0 * g_z_0_0_xxyyzz[k] - g_zz_0_0_xxyyzz[k] * ab_z + g_zz_0_0_xxyyzzz[k];

                g_zz_0_z_xxyzzz[k] = -2.0 * g_z_0_0_xxyzzz[k] - g_zz_0_0_xxyzzz[k] * ab_z + g_zz_0_0_xxyzzzz[k];

                g_zz_0_z_xxzzzz[k] = -2.0 * g_z_0_0_xxzzzz[k] - g_zz_0_0_xxzzzz[k] * ab_z + g_zz_0_0_xxzzzzz[k];

                g_zz_0_z_xyyyyy[k] = -2.0 * g_z_0_0_xyyyyy[k] - g_zz_0_0_xyyyyy[k] * ab_z + g_zz_0_0_xyyyyyz[k];

                g_zz_0_z_xyyyyz[k] = -2.0 * g_z_0_0_xyyyyz[k] - g_zz_0_0_xyyyyz[k] * ab_z + g_zz_0_0_xyyyyzz[k];

                g_zz_0_z_xyyyzz[k] = -2.0 * g_z_0_0_xyyyzz[k] - g_zz_0_0_xyyyzz[k] * ab_z + g_zz_0_0_xyyyzzz[k];

                g_zz_0_z_xyyzzz[k] = -2.0 * g_z_0_0_xyyzzz[k] - g_zz_0_0_xyyzzz[k] * ab_z + g_zz_0_0_xyyzzzz[k];

                g_zz_0_z_xyzzzz[k] = -2.0 * g_z_0_0_xyzzzz[k] - g_zz_0_0_xyzzzz[k] * ab_z + g_zz_0_0_xyzzzzz[k];

                g_zz_0_z_xzzzzz[k] = -2.0 * g_z_0_0_xzzzzz[k] - g_zz_0_0_xzzzzz[k] * ab_z + g_zz_0_0_xzzzzzz[k];

                g_zz_0_z_yyyyyy[k] = -2.0 * g_z_0_0_yyyyyy[k] - g_zz_0_0_yyyyyy[k] * ab_z + g_zz_0_0_yyyyyyz[k];

                g_zz_0_z_yyyyyz[k] = -2.0 * g_z_0_0_yyyyyz[k] - g_zz_0_0_yyyyyz[k] * ab_z + g_zz_0_0_yyyyyzz[k];

                g_zz_0_z_yyyyzz[k] = -2.0 * g_z_0_0_yyyyzz[k] - g_zz_0_0_yyyyzz[k] * ab_z + g_zz_0_0_yyyyzzz[k];

                g_zz_0_z_yyyzzz[k] = -2.0 * g_z_0_0_yyyzzz[k] - g_zz_0_0_yyyzzz[k] * ab_z + g_zz_0_0_yyyzzzz[k];

                g_zz_0_z_yyzzzz[k] = -2.0 * g_z_0_0_yyzzzz[k] - g_zz_0_0_yyzzzz[k] * ab_z + g_zz_0_0_yyzzzzz[k];

                g_zz_0_z_yzzzzz[k] = -2.0 * g_z_0_0_yzzzzz[k] - g_zz_0_0_yzzzzz[k] * ab_z + g_zz_0_0_yzzzzzz[k];

                g_zz_0_z_zzzzzz[k] = -2.0 * g_z_0_0_zzzzzz[k] - g_zz_0_0_zzzzzz[k] * ab_z + g_zz_0_0_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

