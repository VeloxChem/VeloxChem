#include "ElectronRepulsionGeom1000ContrRecDIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_dixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_dixx,
                                            const size_t idx_pixx,
                                            const size_t idx_geom_10_pixx,
                                            const size_t idx_geom_10_pkxx,
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
            /// Set up components of auxilary buffer : PISS

            const auto pi_off = idx_pixx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : PISS

            const auto pi_geom_10_off = idx_geom_10_pixx + i * dcomps + j;

            auto g_x_0_x_xxxxxx = cbuffer.data(pi_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxxxxy = cbuffer.data(pi_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxxxxz = cbuffer.data(pi_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xxxxyy = cbuffer.data(pi_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xxxxyz = cbuffer.data(pi_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xxxxzz = cbuffer.data(pi_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_xxxyyy = cbuffer.data(pi_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_xxxyyz = cbuffer.data(pi_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_xxxyzz = cbuffer.data(pi_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_xxxzzz = cbuffer.data(pi_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_x_xxyyyy = cbuffer.data(pi_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_x_xxyyyz = cbuffer.data(pi_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_x_xxyyzz = cbuffer.data(pi_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_x_xxyzzz = cbuffer.data(pi_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_x_xxzzzz = cbuffer.data(pi_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_x_xyyyyy = cbuffer.data(pi_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_x_xyyyyz = cbuffer.data(pi_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_x_xyyyzz = cbuffer.data(pi_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_x_xyyzzz = cbuffer.data(pi_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_x_xyzzzz = cbuffer.data(pi_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_x_xzzzzz = cbuffer.data(pi_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_x_yyyyyy = cbuffer.data(pi_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_x_yyyyyz = cbuffer.data(pi_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_x_yyyyzz = cbuffer.data(pi_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_x_yyyzzz = cbuffer.data(pi_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_x_yyzzzz = cbuffer.data(pi_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_x_yzzzzz = cbuffer.data(pi_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_x_zzzzzz = cbuffer.data(pi_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_y_xxxxxx = cbuffer.data(pi_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_y_xxxxxy = cbuffer.data(pi_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_y_xxxxxz = cbuffer.data(pi_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_y_xxxxyy = cbuffer.data(pi_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_y_xxxxyz = cbuffer.data(pi_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_y_xxxxzz = cbuffer.data(pi_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_y_xxxyyy = cbuffer.data(pi_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_y_xxxyyz = cbuffer.data(pi_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_y_xxxyzz = cbuffer.data(pi_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_y_xxxzzz = cbuffer.data(pi_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_y_xxyyyy = cbuffer.data(pi_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_y_xxyyyz = cbuffer.data(pi_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_y_xxyyzz = cbuffer.data(pi_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_y_xxyzzz = cbuffer.data(pi_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_y_xxzzzz = cbuffer.data(pi_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_y_xyyyyy = cbuffer.data(pi_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_y_xyyyyz = cbuffer.data(pi_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_y_xyyyzz = cbuffer.data(pi_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_y_xyyzzz = cbuffer.data(pi_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_y_xyzzzz = cbuffer.data(pi_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_y_xzzzzz = cbuffer.data(pi_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_y_yyyyyy = cbuffer.data(pi_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_y_yyyyyz = cbuffer.data(pi_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_y_yyyyzz = cbuffer.data(pi_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_y_yyyzzz = cbuffer.data(pi_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_y_yyzzzz = cbuffer.data(pi_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_y_yzzzzz = cbuffer.data(pi_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_y_zzzzzz = cbuffer.data(pi_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_z_xxxxxx = cbuffer.data(pi_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_z_xxxxxy = cbuffer.data(pi_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_z_xxxxxz = cbuffer.data(pi_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_z_xxxxyy = cbuffer.data(pi_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_z_xxxxyz = cbuffer.data(pi_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_z_xxxxzz = cbuffer.data(pi_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_z_xxxyyy = cbuffer.data(pi_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_z_xxxyyz = cbuffer.data(pi_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_z_xxxyzz = cbuffer.data(pi_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_z_xxxzzz = cbuffer.data(pi_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_z_xxyyyy = cbuffer.data(pi_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_z_xxyyyz = cbuffer.data(pi_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_z_xxyyzz = cbuffer.data(pi_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_z_xxyzzz = cbuffer.data(pi_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_z_xxzzzz = cbuffer.data(pi_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_z_xyyyyy = cbuffer.data(pi_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_z_xyyyyz = cbuffer.data(pi_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_z_xyyyzz = cbuffer.data(pi_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_z_xyyzzz = cbuffer.data(pi_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_z_xyzzzz = cbuffer.data(pi_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_z_xzzzzz = cbuffer.data(pi_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_z_yyyyyy = cbuffer.data(pi_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_z_yyyyyz = cbuffer.data(pi_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_z_yyyyzz = cbuffer.data(pi_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_z_yyyzzz = cbuffer.data(pi_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_z_yyzzzz = cbuffer.data(pi_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_z_yzzzzz = cbuffer.data(pi_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_z_zzzzzz = cbuffer.data(pi_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_x_xxxxxx = cbuffer.data(pi_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_x_xxxxxy = cbuffer.data(pi_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_x_xxxxxz = cbuffer.data(pi_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_x_xxxxyy = cbuffer.data(pi_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_x_xxxxyz = cbuffer.data(pi_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_x_xxxxzz = cbuffer.data(pi_geom_10_off + 89 * ccomps * dcomps);

            auto g_y_0_x_xxxyyy = cbuffer.data(pi_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_x_xxxyyz = cbuffer.data(pi_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_x_xxxyzz = cbuffer.data(pi_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_x_xxxzzz = cbuffer.data(pi_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_x_xxyyyy = cbuffer.data(pi_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_x_xxyyyz = cbuffer.data(pi_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_x_xxyyzz = cbuffer.data(pi_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_x_xxyzzz = cbuffer.data(pi_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_x_xxzzzz = cbuffer.data(pi_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_x_xyyyyy = cbuffer.data(pi_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_x_xyyyyz = cbuffer.data(pi_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_x_xyyyzz = cbuffer.data(pi_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_x_xyyzzz = cbuffer.data(pi_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_x_xyzzzz = cbuffer.data(pi_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_x_xzzzzz = cbuffer.data(pi_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_x_yyyyyy = cbuffer.data(pi_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_x_yyyyyz = cbuffer.data(pi_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_x_yyyyzz = cbuffer.data(pi_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_x_yyyzzz = cbuffer.data(pi_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_x_yyzzzz = cbuffer.data(pi_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_x_yzzzzz = cbuffer.data(pi_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_x_zzzzzz = cbuffer.data(pi_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_y_xxxxxx = cbuffer.data(pi_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_y_xxxxxy = cbuffer.data(pi_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_y_xxxxxz = cbuffer.data(pi_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_y_xxxxyy = cbuffer.data(pi_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_y_xxxxyz = cbuffer.data(pi_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_y_xxxxzz = cbuffer.data(pi_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_y_xxxyyy = cbuffer.data(pi_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_y_xxxyyz = cbuffer.data(pi_geom_10_off + 119 * ccomps * dcomps);

            auto g_y_0_y_xxxyzz = cbuffer.data(pi_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_y_xxxzzz = cbuffer.data(pi_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_y_xxyyyy = cbuffer.data(pi_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_y_xxyyyz = cbuffer.data(pi_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_y_xxyyzz = cbuffer.data(pi_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_y_xxyzzz = cbuffer.data(pi_geom_10_off + 125 * ccomps * dcomps);

            auto g_y_0_y_xxzzzz = cbuffer.data(pi_geom_10_off + 126 * ccomps * dcomps);

            auto g_y_0_y_xyyyyy = cbuffer.data(pi_geom_10_off + 127 * ccomps * dcomps);

            auto g_y_0_y_xyyyyz = cbuffer.data(pi_geom_10_off + 128 * ccomps * dcomps);

            auto g_y_0_y_xyyyzz = cbuffer.data(pi_geom_10_off + 129 * ccomps * dcomps);

            auto g_y_0_y_xyyzzz = cbuffer.data(pi_geom_10_off + 130 * ccomps * dcomps);

            auto g_y_0_y_xyzzzz = cbuffer.data(pi_geom_10_off + 131 * ccomps * dcomps);

            auto g_y_0_y_xzzzzz = cbuffer.data(pi_geom_10_off + 132 * ccomps * dcomps);

            auto g_y_0_y_yyyyyy = cbuffer.data(pi_geom_10_off + 133 * ccomps * dcomps);

            auto g_y_0_y_yyyyyz = cbuffer.data(pi_geom_10_off + 134 * ccomps * dcomps);

            auto g_y_0_y_yyyyzz = cbuffer.data(pi_geom_10_off + 135 * ccomps * dcomps);

            auto g_y_0_y_yyyzzz = cbuffer.data(pi_geom_10_off + 136 * ccomps * dcomps);

            auto g_y_0_y_yyzzzz = cbuffer.data(pi_geom_10_off + 137 * ccomps * dcomps);

            auto g_y_0_y_yzzzzz = cbuffer.data(pi_geom_10_off + 138 * ccomps * dcomps);

            auto g_y_0_y_zzzzzz = cbuffer.data(pi_geom_10_off + 139 * ccomps * dcomps);

            auto g_y_0_z_xxxxxx = cbuffer.data(pi_geom_10_off + 140 * ccomps * dcomps);

            auto g_y_0_z_xxxxxy = cbuffer.data(pi_geom_10_off + 141 * ccomps * dcomps);

            auto g_y_0_z_xxxxxz = cbuffer.data(pi_geom_10_off + 142 * ccomps * dcomps);

            auto g_y_0_z_xxxxyy = cbuffer.data(pi_geom_10_off + 143 * ccomps * dcomps);

            auto g_y_0_z_xxxxyz = cbuffer.data(pi_geom_10_off + 144 * ccomps * dcomps);

            auto g_y_0_z_xxxxzz = cbuffer.data(pi_geom_10_off + 145 * ccomps * dcomps);

            auto g_y_0_z_xxxyyy = cbuffer.data(pi_geom_10_off + 146 * ccomps * dcomps);

            auto g_y_0_z_xxxyyz = cbuffer.data(pi_geom_10_off + 147 * ccomps * dcomps);

            auto g_y_0_z_xxxyzz = cbuffer.data(pi_geom_10_off + 148 * ccomps * dcomps);

            auto g_y_0_z_xxxzzz = cbuffer.data(pi_geom_10_off + 149 * ccomps * dcomps);

            auto g_y_0_z_xxyyyy = cbuffer.data(pi_geom_10_off + 150 * ccomps * dcomps);

            auto g_y_0_z_xxyyyz = cbuffer.data(pi_geom_10_off + 151 * ccomps * dcomps);

            auto g_y_0_z_xxyyzz = cbuffer.data(pi_geom_10_off + 152 * ccomps * dcomps);

            auto g_y_0_z_xxyzzz = cbuffer.data(pi_geom_10_off + 153 * ccomps * dcomps);

            auto g_y_0_z_xxzzzz = cbuffer.data(pi_geom_10_off + 154 * ccomps * dcomps);

            auto g_y_0_z_xyyyyy = cbuffer.data(pi_geom_10_off + 155 * ccomps * dcomps);

            auto g_y_0_z_xyyyyz = cbuffer.data(pi_geom_10_off + 156 * ccomps * dcomps);

            auto g_y_0_z_xyyyzz = cbuffer.data(pi_geom_10_off + 157 * ccomps * dcomps);

            auto g_y_0_z_xyyzzz = cbuffer.data(pi_geom_10_off + 158 * ccomps * dcomps);

            auto g_y_0_z_xyzzzz = cbuffer.data(pi_geom_10_off + 159 * ccomps * dcomps);

            auto g_y_0_z_xzzzzz = cbuffer.data(pi_geom_10_off + 160 * ccomps * dcomps);

            auto g_y_0_z_yyyyyy = cbuffer.data(pi_geom_10_off + 161 * ccomps * dcomps);

            auto g_y_0_z_yyyyyz = cbuffer.data(pi_geom_10_off + 162 * ccomps * dcomps);

            auto g_y_0_z_yyyyzz = cbuffer.data(pi_geom_10_off + 163 * ccomps * dcomps);

            auto g_y_0_z_yyyzzz = cbuffer.data(pi_geom_10_off + 164 * ccomps * dcomps);

            auto g_y_0_z_yyzzzz = cbuffer.data(pi_geom_10_off + 165 * ccomps * dcomps);

            auto g_y_0_z_yzzzzz = cbuffer.data(pi_geom_10_off + 166 * ccomps * dcomps);

            auto g_y_0_z_zzzzzz = cbuffer.data(pi_geom_10_off + 167 * ccomps * dcomps);

            auto g_z_0_x_xxxxxx = cbuffer.data(pi_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_x_xxxxxy = cbuffer.data(pi_geom_10_off + 169 * ccomps * dcomps);

            auto g_z_0_x_xxxxxz = cbuffer.data(pi_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_x_xxxxyy = cbuffer.data(pi_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_x_xxxxyz = cbuffer.data(pi_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_x_xxxxzz = cbuffer.data(pi_geom_10_off + 173 * ccomps * dcomps);

            auto g_z_0_x_xxxyyy = cbuffer.data(pi_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_x_xxxyyz = cbuffer.data(pi_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_x_xxxyzz = cbuffer.data(pi_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_x_xxxzzz = cbuffer.data(pi_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_x_xxyyyy = cbuffer.data(pi_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_x_xxyyyz = cbuffer.data(pi_geom_10_off + 179 * ccomps * dcomps);

            auto g_z_0_x_xxyyzz = cbuffer.data(pi_geom_10_off + 180 * ccomps * dcomps);

            auto g_z_0_x_xxyzzz = cbuffer.data(pi_geom_10_off + 181 * ccomps * dcomps);

            auto g_z_0_x_xxzzzz = cbuffer.data(pi_geom_10_off + 182 * ccomps * dcomps);

            auto g_z_0_x_xyyyyy = cbuffer.data(pi_geom_10_off + 183 * ccomps * dcomps);

            auto g_z_0_x_xyyyyz = cbuffer.data(pi_geom_10_off + 184 * ccomps * dcomps);

            auto g_z_0_x_xyyyzz = cbuffer.data(pi_geom_10_off + 185 * ccomps * dcomps);

            auto g_z_0_x_xyyzzz = cbuffer.data(pi_geom_10_off + 186 * ccomps * dcomps);

            auto g_z_0_x_xyzzzz = cbuffer.data(pi_geom_10_off + 187 * ccomps * dcomps);

            auto g_z_0_x_xzzzzz = cbuffer.data(pi_geom_10_off + 188 * ccomps * dcomps);

            auto g_z_0_x_yyyyyy = cbuffer.data(pi_geom_10_off + 189 * ccomps * dcomps);

            auto g_z_0_x_yyyyyz = cbuffer.data(pi_geom_10_off + 190 * ccomps * dcomps);

            auto g_z_0_x_yyyyzz = cbuffer.data(pi_geom_10_off + 191 * ccomps * dcomps);

            auto g_z_0_x_yyyzzz = cbuffer.data(pi_geom_10_off + 192 * ccomps * dcomps);

            auto g_z_0_x_yyzzzz = cbuffer.data(pi_geom_10_off + 193 * ccomps * dcomps);

            auto g_z_0_x_yzzzzz = cbuffer.data(pi_geom_10_off + 194 * ccomps * dcomps);

            auto g_z_0_x_zzzzzz = cbuffer.data(pi_geom_10_off + 195 * ccomps * dcomps);

            auto g_z_0_y_xxxxxx = cbuffer.data(pi_geom_10_off + 196 * ccomps * dcomps);

            auto g_z_0_y_xxxxxy = cbuffer.data(pi_geom_10_off + 197 * ccomps * dcomps);

            auto g_z_0_y_xxxxxz = cbuffer.data(pi_geom_10_off + 198 * ccomps * dcomps);

            auto g_z_0_y_xxxxyy = cbuffer.data(pi_geom_10_off + 199 * ccomps * dcomps);

            auto g_z_0_y_xxxxyz = cbuffer.data(pi_geom_10_off + 200 * ccomps * dcomps);

            auto g_z_0_y_xxxxzz = cbuffer.data(pi_geom_10_off + 201 * ccomps * dcomps);

            auto g_z_0_y_xxxyyy = cbuffer.data(pi_geom_10_off + 202 * ccomps * dcomps);

            auto g_z_0_y_xxxyyz = cbuffer.data(pi_geom_10_off + 203 * ccomps * dcomps);

            auto g_z_0_y_xxxyzz = cbuffer.data(pi_geom_10_off + 204 * ccomps * dcomps);

            auto g_z_0_y_xxxzzz = cbuffer.data(pi_geom_10_off + 205 * ccomps * dcomps);

            auto g_z_0_y_xxyyyy = cbuffer.data(pi_geom_10_off + 206 * ccomps * dcomps);

            auto g_z_0_y_xxyyyz = cbuffer.data(pi_geom_10_off + 207 * ccomps * dcomps);

            auto g_z_0_y_xxyyzz = cbuffer.data(pi_geom_10_off + 208 * ccomps * dcomps);

            auto g_z_0_y_xxyzzz = cbuffer.data(pi_geom_10_off + 209 * ccomps * dcomps);

            auto g_z_0_y_xxzzzz = cbuffer.data(pi_geom_10_off + 210 * ccomps * dcomps);

            auto g_z_0_y_xyyyyy = cbuffer.data(pi_geom_10_off + 211 * ccomps * dcomps);

            auto g_z_0_y_xyyyyz = cbuffer.data(pi_geom_10_off + 212 * ccomps * dcomps);

            auto g_z_0_y_xyyyzz = cbuffer.data(pi_geom_10_off + 213 * ccomps * dcomps);

            auto g_z_0_y_xyyzzz = cbuffer.data(pi_geom_10_off + 214 * ccomps * dcomps);

            auto g_z_0_y_xyzzzz = cbuffer.data(pi_geom_10_off + 215 * ccomps * dcomps);

            auto g_z_0_y_xzzzzz = cbuffer.data(pi_geom_10_off + 216 * ccomps * dcomps);

            auto g_z_0_y_yyyyyy = cbuffer.data(pi_geom_10_off + 217 * ccomps * dcomps);

            auto g_z_0_y_yyyyyz = cbuffer.data(pi_geom_10_off + 218 * ccomps * dcomps);

            auto g_z_0_y_yyyyzz = cbuffer.data(pi_geom_10_off + 219 * ccomps * dcomps);

            auto g_z_0_y_yyyzzz = cbuffer.data(pi_geom_10_off + 220 * ccomps * dcomps);

            auto g_z_0_y_yyzzzz = cbuffer.data(pi_geom_10_off + 221 * ccomps * dcomps);

            auto g_z_0_y_yzzzzz = cbuffer.data(pi_geom_10_off + 222 * ccomps * dcomps);

            auto g_z_0_y_zzzzzz = cbuffer.data(pi_geom_10_off + 223 * ccomps * dcomps);

            auto g_z_0_z_xxxxxx = cbuffer.data(pi_geom_10_off + 224 * ccomps * dcomps);

            auto g_z_0_z_xxxxxy = cbuffer.data(pi_geom_10_off + 225 * ccomps * dcomps);

            auto g_z_0_z_xxxxxz = cbuffer.data(pi_geom_10_off + 226 * ccomps * dcomps);

            auto g_z_0_z_xxxxyy = cbuffer.data(pi_geom_10_off + 227 * ccomps * dcomps);

            auto g_z_0_z_xxxxyz = cbuffer.data(pi_geom_10_off + 228 * ccomps * dcomps);

            auto g_z_0_z_xxxxzz = cbuffer.data(pi_geom_10_off + 229 * ccomps * dcomps);

            auto g_z_0_z_xxxyyy = cbuffer.data(pi_geom_10_off + 230 * ccomps * dcomps);

            auto g_z_0_z_xxxyyz = cbuffer.data(pi_geom_10_off + 231 * ccomps * dcomps);

            auto g_z_0_z_xxxyzz = cbuffer.data(pi_geom_10_off + 232 * ccomps * dcomps);

            auto g_z_0_z_xxxzzz = cbuffer.data(pi_geom_10_off + 233 * ccomps * dcomps);

            auto g_z_0_z_xxyyyy = cbuffer.data(pi_geom_10_off + 234 * ccomps * dcomps);

            auto g_z_0_z_xxyyyz = cbuffer.data(pi_geom_10_off + 235 * ccomps * dcomps);

            auto g_z_0_z_xxyyzz = cbuffer.data(pi_geom_10_off + 236 * ccomps * dcomps);

            auto g_z_0_z_xxyzzz = cbuffer.data(pi_geom_10_off + 237 * ccomps * dcomps);

            auto g_z_0_z_xxzzzz = cbuffer.data(pi_geom_10_off + 238 * ccomps * dcomps);

            auto g_z_0_z_xyyyyy = cbuffer.data(pi_geom_10_off + 239 * ccomps * dcomps);

            auto g_z_0_z_xyyyyz = cbuffer.data(pi_geom_10_off + 240 * ccomps * dcomps);

            auto g_z_0_z_xyyyzz = cbuffer.data(pi_geom_10_off + 241 * ccomps * dcomps);

            auto g_z_0_z_xyyzzz = cbuffer.data(pi_geom_10_off + 242 * ccomps * dcomps);

            auto g_z_0_z_xyzzzz = cbuffer.data(pi_geom_10_off + 243 * ccomps * dcomps);

            auto g_z_0_z_xzzzzz = cbuffer.data(pi_geom_10_off + 244 * ccomps * dcomps);

            auto g_z_0_z_yyyyyy = cbuffer.data(pi_geom_10_off + 245 * ccomps * dcomps);

            auto g_z_0_z_yyyyyz = cbuffer.data(pi_geom_10_off + 246 * ccomps * dcomps);

            auto g_z_0_z_yyyyzz = cbuffer.data(pi_geom_10_off + 247 * ccomps * dcomps);

            auto g_z_0_z_yyyzzz = cbuffer.data(pi_geom_10_off + 248 * ccomps * dcomps);

            auto g_z_0_z_yyzzzz = cbuffer.data(pi_geom_10_off + 249 * ccomps * dcomps);

            auto g_z_0_z_yzzzzz = cbuffer.data(pi_geom_10_off + 250 * ccomps * dcomps);

            auto g_z_0_z_zzzzzz = cbuffer.data(pi_geom_10_off + 251 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PKSS

            const auto pk_geom_10_off = idx_geom_10_pkxx + i * dcomps + j;

            auto g_x_0_x_xxxxxxx = cbuffer.data(pk_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxxxxxy = cbuffer.data(pk_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxxxxxz = cbuffer.data(pk_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xxxxxyy = cbuffer.data(pk_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xxxxxyz = cbuffer.data(pk_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xxxxxzz = cbuffer.data(pk_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_xxxxyyy = cbuffer.data(pk_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_xxxxyyz = cbuffer.data(pk_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_xxxxyzz = cbuffer.data(pk_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_xxxxzzz = cbuffer.data(pk_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_x_xxxyyyy = cbuffer.data(pk_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_x_xxxyyyz = cbuffer.data(pk_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_x_xxxyyzz = cbuffer.data(pk_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_x_xxxyzzz = cbuffer.data(pk_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_x_xxxzzzz = cbuffer.data(pk_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_x_xxyyyyy = cbuffer.data(pk_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_x_xxyyyyz = cbuffer.data(pk_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_x_xxyyyzz = cbuffer.data(pk_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_x_xxyyzzz = cbuffer.data(pk_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_x_xxyzzzz = cbuffer.data(pk_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_x_xxzzzzz = cbuffer.data(pk_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_x_xyyyyyy = cbuffer.data(pk_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_x_xyyyyyz = cbuffer.data(pk_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_x_xyyyyzz = cbuffer.data(pk_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_x_xyyyzzz = cbuffer.data(pk_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_x_xyyzzzz = cbuffer.data(pk_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_x_xyzzzzz = cbuffer.data(pk_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_x_xzzzzzz = cbuffer.data(pk_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_x_yyyyyyy = cbuffer.data(pk_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_x_yyyyyyz = cbuffer.data(pk_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_x_yyyyyzz = cbuffer.data(pk_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_x_yyyyzzz = cbuffer.data(pk_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_x_yyyzzzz = cbuffer.data(pk_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_x_yyzzzzz = cbuffer.data(pk_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_x_yzzzzzz = cbuffer.data(pk_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_x_zzzzzzz = cbuffer.data(pk_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_y_xxxxxxy = cbuffer.data(pk_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_y_xxxxxyy = cbuffer.data(pk_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_y_xxxxxyz = cbuffer.data(pk_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_y_xxxxyyy = cbuffer.data(pk_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_y_xxxxyyz = cbuffer.data(pk_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_y_xxxxyzz = cbuffer.data(pk_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_y_xxxyyyy = cbuffer.data(pk_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_y_xxxyyyz = cbuffer.data(pk_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_y_xxxyyzz = cbuffer.data(pk_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_y_xxxyzzz = cbuffer.data(pk_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_y_xxyyyyy = cbuffer.data(pk_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_y_xxyyyyz = cbuffer.data(pk_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_y_xxyyyzz = cbuffer.data(pk_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_y_xxyyzzz = cbuffer.data(pk_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_y_xxyzzzz = cbuffer.data(pk_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_y_xyyyyyy = cbuffer.data(pk_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_y_xyyyyyz = cbuffer.data(pk_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_y_xyyyyzz = cbuffer.data(pk_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_y_xyyyzzz = cbuffer.data(pk_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_y_xyyzzzz = cbuffer.data(pk_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_y_xyzzzzz = cbuffer.data(pk_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_y_yyyyyyy = cbuffer.data(pk_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_y_yyyyyyz = cbuffer.data(pk_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_y_yyyyyzz = cbuffer.data(pk_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_y_yyyyzzz = cbuffer.data(pk_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_y_yyyzzzz = cbuffer.data(pk_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_y_yyzzzzz = cbuffer.data(pk_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_y_yzzzzzz = cbuffer.data(pk_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_z_xxxxxxy = cbuffer.data(pk_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_z_xxxxxxz = cbuffer.data(pk_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_z_xxxxxyy = cbuffer.data(pk_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_z_xxxxxyz = cbuffer.data(pk_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_z_xxxxxzz = cbuffer.data(pk_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_z_xxxxyyy = cbuffer.data(pk_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_z_xxxxyyz = cbuffer.data(pk_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_z_xxxxyzz = cbuffer.data(pk_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_z_xxxxzzz = cbuffer.data(pk_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_z_xxxyyyy = cbuffer.data(pk_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_z_xxxyyyz = cbuffer.data(pk_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_z_xxxyyzz = cbuffer.data(pk_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_z_xxxyzzz = cbuffer.data(pk_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_z_xxxzzzz = cbuffer.data(pk_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_z_xxyyyyy = cbuffer.data(pk_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_z_xxyyyyz = cbuffer.data(pk_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_z_xxyyyzz = cbuffer.data(pk_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_z_xxyyzzz = cbuffer.data(pk_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_z_xxyzzzz = cbuffer.data(pk_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_z_xxzzzzz = cbuffer.data(pk_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_z_xyyyyyy = cbuffer.data(pk_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_z_xyyyyyz = cbuffer.data(pk_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_z_xyyyyzz = cbuffer.data(pk_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_z_xyyyzzz = cbuffer.data(pk_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_z_xyyzzzz = cbuffer.data(pk_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_z_xyzzzzz = cbuffer.data(pk_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_z_xzzzzzz = cbuffer.data(pk_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_z_yyyyyyy = cbuffer.data(pk_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_z_yyyyyyz = cbuffer.data(pk_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_z_yyyyyzz = cbuffer.data(pk_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_z_yyyyzzz = cbuffer.data(pk_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_z_yyyzzzz = cbuffer.data(pk_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_z_yyzzzzz = cbuffer.data(pk_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_z_yzzzzzz = cbuffer.data(pk_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_z_zzzzzzz = cbuffer.data(pk_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_x_xxxxxxx = cbuffer.data(pk_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_x_xxxxxxy = cbuffer.data(pk_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_x_xxxxxxz = cbuffer.data(pk_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_x_xxxxxyy = cbuffer.data(pk_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_x_xxxxxyz = cbuffer.data(pk_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_x_xxxxxzz = cbuffer.data(pk_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_x_xxxxyyy = cbuffer.data(pk_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_x_xxxxyyz = cbuffer.data(pk_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_x_xxxxyzz = cbuffer.data(pk_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_x_xxxxzzz = cbuffer.data(pk_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_x_xxxyyyy = cbuffer.data(pk_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_x_xxxyyyz = cbuffer.data(pk_geom_10_off + 119 * ccomps * dcomps);

            auto g_y_0_x_xxxyyzz = cbuffer.data(pk_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_x_xxxyzzz = cbuffer.data(pk_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_x_xxxzzzz = cbuffer.data(pk_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_x_xxyyyyy = cbuffer.data(pk_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_x_xxyyyyz = cbuffer.data(pk_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_x_xxyyyzz = cbuffer.data(pk_geom_10_off + 125 * ccomps * dcomps);

            auto g_y_0_x_xxyyzzz = cbuffer.data(pk_geom_10_off + 126 * ccomps * dcomps);

            auto g_y_0_x_xxyzzzz = cbuffer.data(pk_geom_10_off + 127 * ccomps * dcomps);

            auto g_y_0_x_xxzzzzz = cbuffer.data(pk_geom_10_off + 128 * ccomps * dcomps);

            auto g_y_0_x_xyyyyyy = cbuffer.data(pk_geom_10_off + 129 * ccomps * dcomps);

            auto g_y_0_x_xyyyyyz = cbuffer.data(pk_geom_10_off + 130 * ccomps * dcomps);

            auto g_y_0_x_xyyyyzz = cbuffer.data(pk_geom_10_off + 131 * ccomps * dcomps);

            auto g_y_0_x_xyyyzzz = cbuffer.data(pk_geom_10_off + 132 * ccomps * dcomps);

            auto g_y_0_x_xyyzzzz = cbuffer.data(pk_geom_10_off + 133 * ccomps * dcomps);

            auto g_y_0_x_xyzzzzz = cbuffer.data(pk_geom_10_off + 134 * ccomps * dcomps);

            auto g_y_0_x_xzzzzzz = cbuffer.data(pk_geom_10_off + 135 * ccomps * dcomps);

            auto g_y_0_y_xxxxxxx = cbuffer.data(pk_geom_10_off + 144 * ccomps * dcomps);

            auto g_y_0_y_xxxxxxy = cbuffer.data(pk_geom_10_off + 145 * ccomps * dcomps);

            auto g_y_0_y_xxxxxxz = cbuffer.data(pk_geom_10_off + 146 * ccomps * dcomps);

            auto g_y_0_y_xxxxxyy = cbuffer.data(pk_geom_10_off + 147 * ccomps * dcomps);

            auto g_y_0_y_xxxxxyz = cbuffer.data(pk_geom_10_off + 148 * ccomps * dcomps);

            auto g_y_0_y_xxxxxzz = cbuffer.data(pk_geom_10_off + 149 * ccomps * dcomps);

            auto g_y_0_y_xxxxyyy = cbuffer.data(pk_geom_10_off + 150 * ccomps * dcomps);

            auto g_y_0_y_xxxxyyz = cbuffer.data(pk_geom_10_off + 151 * ccomps * dcomps);

            auto g_y_0_y_xxxxyzz = cbuffer.data(pk_geom_10_off + 152 * ccomps * dcomps);

            auto g_y_0_y_xxxxzzz = cbuffer.data(pk_geom_10_off + 153 * ccomps * dcomps);

            auto g_y_0_y_xxxyyyy = cbuffer.data(pk_geom_10_off + 154 * ccomps * dcomps);

            auto g_y_0_y_xxxyyyz = cbuffer.data(pk_geom_10_off + 155 * ccomps * dcomps);

            auto g_y_0_y_xxxyyzz = cbuffer.data(pk_geom_10_off + 156 * ccomps * dcomps);

            auto g_y_0_y_xxxyzzz = cbuffer.data(pk_geom_10_off + 157 * ccomps * dcomps);

            auto g_y_0_y_xxxzzzz = cbuffer.data(pk_geom_10_off + 158 * ccomps * dcomps);

            auto g_y_0_y_xxyyyyy = cbuffer.data(pk_geom_10_off + 159 * ccomps * dcomps);

            auto g_y_0_y_xxyyyyz = cbuffer.data(pk_geom_10_off + 160 * ccomps * dcomps);

            auto g_y_0_y_xxyyyzz = cbuffer.data(pk_geom_10_off + 161 * ccomps * dcomps);

            auto g_y_0_y_xxyyzzz = cbuffer.data(pk_geom_10_off + 162 * ccomps * dcomps);

            auto g_y_0_y_xxyzzzz = cbuffer.data(pk_geom_10_off + 163 * ccomps * dcomps);

            auto g_y_0_y_xxzzzzz = cbuffer.data(pk_geom_10_off + 164 * ccomps * dcomps);

            auto g_y_0_y_xyyyyyy = cbuffer.data(pk_geom_10_off + 165 * ccomps * dcomps);

            auto g_y_0_y_xyyyyyz = cbuffer.data(pk_geom_10_off + 166 * ccomps * dcomps);

            auto g_y_0_y_xyyyyzz = cbuffer.data(pk_geom_10_off + 167 * ccomps * dcomps);

            auto g_y_0_y_xyyyzzz = cbuffer.data(pk_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_y_xyyzzzz = cbuffer.data(pk_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_y_xyzzzzz = cbuffer.data(pk_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_y_xzzzzzz = cbuffer.data(pk_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_y_yyyyyyy = cbuffer.data(pk_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_y_yyyyyyz = cbuffer.data(pk_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_y_yyyyyzz = cbuffer.data(pk_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_y_yyyyzzz = cbuffer.data(pk_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_y_yyyzzzz = cbuffer.data(pk_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_y_yyzzzzz = cbuffer.data(pk_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_y_yzzzzzz = cbuffer.data(pk_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_y_zzzzzzz = cbuffer.data(pk_geom_10_off + 179 * ccomps * dcomps);

            auto g_y_0_z_xxxxxxx = cbuffer.data(pk_geom_10_off + 180 * ccomps * dcomps);

            auto g_y_0_z_xxxxxxy = cbuffer.data(pk_geom_10_off + 181 * ccomps * dcomps);

            auto g_y_0_z_xxxxxxz = cbuffer.data(pk_geom_10_off + 182 * ccomps * dcomps);

            auto g_y_0_z_xxxxxyy = cbuffer.data(pk_geom_10_off + 183 * ccomps * dcomps);

            auto g_y_0_z_xxxxxyz = cbuffer.data(pk_geom_10_off + 184 * ccomps * dcomps);

            auto g_y_0_z_xxxxxzz = cbuffer.data(pk_geom_10_off + 185 * ccomps * dcomps);

            auto g_y_0_z_xxxxyyy = cbuffer.data(pk_geom_10_off + 186 * ccomps * dcomps);

            auto g_y_0_z_xxxxyyz = cbuffer.data(pk_geom_10_off + 187 * ccomps * dcomps);

            auto g_y_0_z_xxxxyzz = cbuffer.data(pk_geom_10_off + 188 * ccomps * dcomps);

            auto g_y_0_z_xxxxzzz = cbuffer.data(pk_geom_10_off + 189 * ccomps * dcomps);

            auto g_y_0_z_xxxyyyy = cbuffer.data(pk_geom_10_off + 190 * ccomps * dcomps);

            auto g_y_0_z_xxxyyyz = cbuffer.data(pk_geom_10_off + 191 * ccomps * dcomps);

            auto g_y_0_z_xxxyyzz = cbuffer.data(pk_geom_10_off + 192 * ccomps * dcomps);

            auto g_y_0_z_xxxyzzz = cbuffer.data(pk_geom_10_off + 193 * ccomps * dcomps);

            auto g_y_0_z_xxxzzzz = cbuffer.data(pk_geom_10_off + 194 * ccomps * dcomps);

            auto g_y_0_z_xxyyyyy = cbuffer.data(pk_geom_10_off + 195 * ccomps * dcomps);

            auto g_y_0_z_xxyyyyz = cbuffer.data(pk_geom_10_off + 196 * ccomps * dcomps);

            auto g_y_0_z_xxyyyzz = cbuffer.data(pk_geom_10_off + 197 * ccomps * dcomps);

            auto g_y_0_z_xxyyzzz = cbuffer.data(pk_geom_10_off + 198 * ccomps * dcomps);

            auto g_y_0_z_xxyzzzz = cbuffer.data(pk_geom_10_off + 199 * ccomps * dcomps);

            auto g_y_0_z_xxzzzzz = cbuffer.data(pk_geom_10_off + 200 * ccomps * dcomps);

            auto g_y_0_z_xyyyyyy = cbuffer.data(pk_geom_10_off + 201 * ccomps * dcomps);

            auto g_y_0_z_xyyyyyz = cbuffer.data(pk_geom_10_off + 202 * ccomps * dcomps);

            auto g_y_0_z_xyyyyzz = cbuffer.data(pk_geom_10_off + 203 * ccomps * dcomps);

            auto g_y_0_z_xyyyzzz = cbuffer.data(pk_geom_10_off + 204 * ccomps * dcomps);

            auto g_y_0_z_xyyzzzz = cbuffer.data(pk_geom_10_off + 205 * ccomps * dcomps);

            auto g_y_0_z_xyzzzzz = cbuffer.data(pk_geom_10_off + 206 * ccomps * dcomps);

            auto g_y_0_z_xzzzzzz = cbuffer.data(pk_geom_10_off + 207 * ccomps * dcomps);

            auto g_y_0_z_yyyyyyz = cbuffer.data(pk_geom_10_off + 209 * ccomps * dcomps);

            auto g_y_0_z_yyyyyzz = cbuffer.data(pk_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_z_yyyyzzz = cbuffer.data(pk_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_z_yyyzzzz = cbuffer.data(pk_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_z_yyzzzzz = cbuffer.data(pk_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_z_yzzzzzz = cbuffer.data(pk_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_z_zzzzzzz = cbuffer.data(pk_geom_10_off + 215 * ccomps * dcomps);

            auto g_z_0_x_xxxxxxx = cbuffer.data(pk_geom_10_off + 216 * ccomps * dcomps);

            auto g_z_0_x_xxxxxxy = cbuffer.data(pk_geom_10_off + 217 * ccomps * dcomps);

            auto g_z_0_x_xxxxxxz = cbuffer.data(pk_geom_10_off + 218 * ccomps * dcomps);

            auto g_z_0_x_xxxxxyy = cbuffer.data(pk_geom_10_off + 219 * ccomps * dcomps);

            auto g_z_0_x_xxxxxyz = cbuffer.data(pk_geom_10_off + 220 * ccomps * dcomps);

            auto g_z_0_x_xxxxxzz = cbuffer.data(pk_geom_10_off + 221 * ccomps * dcomps);

            auto g_z_0_x_xxxxyyy = cbuffer.data(pk_geom_10_off + 222 * ccomps * dcomps);

            auto g_z_0_x_xxxxyyz = cbuffer.data(pk_geom_10_off + 223 * ccomps * dcomps);

            auto g_z_0_x_xxxxyzz = cbuffer.data(pk_geom_10_off + 224 * ccomps * dcomps);

            auto g_z_0_x_xxxxzzz = cbuffer.data(pk_geom_10_off + 225 * ccomps * dcomps);

            auto g_z_0_x_xxxyyyy = cbuffer.data(pk_geom_10_off + 226 * ccomps * dcomps);

            auto g_z_0_x_xxxyyyz = cbuffer.data(pk_geom_10_off + 227 * ccomps * dcomps);

            auto g_z_0_x_xxxyyzz = cbuffer.data(pk_geom_10_off + 228 * ccomps * dcomps);

            auto g_z_0_x_xxxyzzz = cbuffer.data(pk_geom_10_off + 229 * ccomps * dcomps);

            auto g_z_0_x_xxxzzzz = cbuffer.data(pk_geom_10_off + 230 * ccomps * dcomps);

            auto g_z_0_x_xxyyyyy = cbuffer.data(pk_geom_10_off + 231 * ccomps * dcomps);

            auto g_z_0_x_xxyyyyz = cbuffer.data(pk_geom_10_off + 232 * ccomps * dcomps);

            auto g_z_0_x_xxyyyzz = cbuffer.data(pk_geom_10_off + 233 * ccomps * dcomps);

            auto g_z_0_x_xxyyzzz = cbuffer.data(pk_geom_10_off + 234 * ccomps * dcomps);

            auto g_z_0_x_xxyzzzz = cbuffer.data(pk_geom_10_off + 235 * ccomps * dcomps);

            auto g_z_0_x_xxzzzzz = cbuffer.data(pk_geom_10_off + 236 * ccomps * dcomps);

            auto g_z_0_x_xyyyyyy = cbuffer.data(pk_geom_10_off + 237 * ccomps * dcomps);

            auto g_z_0_x_xyyyyyz = cbuffer.data(pk_geom_10_off + 238 * ccomps * dcomps);

            auto g_z_0_x_xyyyyzz = cbuffer.data(pk_geom_10_off + 239 * ccomps * dcomps);

            auto g_z_0_x_xyyyzzz = cbuffer.data(pk_geom_10_off + 240 * ccomps * dcomps);

            auto g_z_0_x_xyyzzzz = cbuffer.data(pk_geom_10_off + 241 * ccomps * dcomps);

            auto g_z_0_x_xyzzzzz = cbuffer.data(pk_geom_10_off + 242 * ccomps * dcomps);

            auto g_z_0_x_xzzzzzz = cbuffer.data(pk_geom_10_off + 243 * ccomps * dcomps);

            auto g_z_0_y_xxxxxxx = cbuffer.data(pk_geom_10_off + 252 * ccomps * dcomps);

            auto g_z_0_y_xxxxxxy = cbuffer.data(pk_geom_10_off + 253 * ccomps * dcomps);

            auto g_z_0_y_xxxxxxz = cbuffer.data(pk_geom_10_off + 254 * ccomps * dcomps);

            auto g_z_0_y_xxxxxyy = cbuffer.data(pk_geom_10_off + 255 * ccomps * dcomps);

            auto g_z_0_y_xxxxxyz = cbuffer.data(pk_geom_10_off + 256 * ccomps * dcomps);

            auto g_z_0_y_xxxxxzz = cbuffer.data(pk_geom_10_off + 257 * ccomps * dcomps);

            auto g_z_0_y_xxxxyyy = cbuffer.data(pk_geom_10_off + 258 * ccomps * dcomps);

            auto g_z_0_y_xxxxyyz = cbuffer.data(pk_geom_10_off + 259 * ccomps * dcomps);

            auto g_z_0_y_xxxxyzz = cbuffer.data(pk_geom_10_off + 260 * ccomps * dcomps);

            auto g_z_0_y_xxxxzzz = cbuffer.data(pk_geom_10_off + 261 * ccomps * dcomps);

            auto g_z_0_y_xxxyyyy = cbuffer.data(pk_geom_10_off + 262 * ccomps * dcomps);

            auto g_z_0_y_xxxyyyz = cbuffer.data(pk_geom_10_off + 263 * ccomps * dcomps);

            auto g_z_0_y_xxxyyzz = cbuffer.data(pk_geom_10_off + 264 * ccomps * dcomps);

            auto g_z_0_y_xxxyzzz = cbuffer.data(pk_geom_10_off + 265 * ccomps * dcomps);

            auto g_z_0_y_xxxzzzz = cbuffer.data(pk_geom_10_off + 266 * ccomps * dcomps);

            auto g_z_0_y_xxyyyyy = cbuffer.data(pk_geom_10_off + 267 * ccomps * dcomps);

            auto g_z_0_y_xxyyyyz = cbuffer.data(pk_geom_10_off + 268 * ccomps * dcomps);

            auto g_z_0_y_xxyyyzz = cbuffer.data(pk_geom_10_off + 269 * ccomps * dcomps);

            auto g_z_0_y_xxyyzzz = cbuffer.data(pk_geom_10_off + 270 * ccomps * dcomps);

            auto g_z_0_y_xxyzzzz = cbuffer.data(pk_geom_10_off + 271 * ccomps * dcomps);

            auto g_z_0_y_xxzzzzz = cbuffer.data(pk_geom_10_off + 272 * ccomps * dcomps);

            auto g_z_0_y_xyyyyyy = cbuffer.data(pk_geom_10_off + 273 * ccomps * dcomps);

            auto g_z_0_y_xyyyyyz = cbuffer.data(pk_geom_10_off + 274 * ccomps * dcomps);

            auto g_z_0_y_xyyyyzz = cbuffer.data(pk_geom_10_off + 275 * ccomps * dcomps);

            auto g_z_0_y_xyyyzzz = cbuffer.data(pk_geom_10_off + 276 * ccomps * dcomps);

            auto g_z_0_y_xyyzzzz = cbuffer.data(pk_geom_10_off + 277 * ccomps * dcomps);

            auto g_z_0_y_xyzzzzz = cbuffer.data(pk_geom_10_off + 278 * ccomps * dcomps);

            auto g_z_0_y_xzzzzzz = cbuffer.data(pk_geom_10_off + 279 * ccomps * dcomps);

            auto g_z_0_y_yyyyyyy = cbuffer.data(pk_geom_10_off + 280 * ccomps * dcomps);

            auto g_z_0_y_yyyyyyz = cbuffer.data(pk_geom_10_off + 281 * ccomps * dcomps);

            auto g_z_0_y_yyyyyzz = cbuffer.data(pk_geom_10_off + 282 * ccomps * dcomps);

            auto g_z_0_y_yyyyzzz = cbuffer.data(pk_geom_10_off + 283 * ccomps * dcomps);

            auto g_z_0_y_yyyzzzz = cbuffer.data(pk_geom_10_off + 284 * ccomps * dcomps);

            auto g_z_0_y_yyzzzzz = cbuffer.data(pk_geom_10_off + 285 * ccomps * dcomps);

            auto g_z_0_y_yzzzzzz = cbuffer.data(pk_geom_10_off + 286 * ccomps * dcomps);

            auto g_z_0_z_xxxxxxx = cbuffer.data(pk_geom_10_off + 288 * ccomps * dcomps);

            auto g_z_0_z_xxxxxxy = cbuffer.data(pk_geom_10_off + 289 * ccomps * dcomps);

            auto g_z_0_z_xxxxxxz = cbuffer.data(pk_geom_10_off + 290 * ccomps * dcomps);

            auto g_z_0_z_xxxxxyy = cbuffer.data(pk_geom_10_off + 291 * ccomps * dcomps);

            auto g_z_0_z_xxxxxyz = cbuffer.data(pk_geom_10_off + 292 * ccomps * dcomps);

            auto g_z_0_z_xxxxxzz = cbuffer.data(pk_geom_10_off + 293 * ccomps * dcomps);

            auto g_z_0_z_xxxxyyy = cbuffer.data(pk_geom_10_off + 294 * ccomps * dcomps);

            auto g_z_0_z_xxxxyyz = cbuffer.data(pk_geom_10_off + 295 * ccomps * dcomps);

            auto g_z_0_z_xxxxyzz = cbuffer.data(pk_geom_10_off + 296 * ccomps * dcomps);

            auto g_z_0_z_xxxxzzz = cbuffer.data(pk_geom_10_off + 297 * ccomps * dcomps);

            auto g_z_0_z_xxxyyyy = cbuffer.data(pk_geom_10_off + 298 * ccomps * dcomps);

            auto g_z_0_z_xxxyyyz = cbuffer.data(pk_geom_10_off + 299 * ccomps * dcomps);

            auto g_z_0_z_xxxyyzz = cbuffer.data(pk_geom_10_off + 300 * ccomps * dcomps);

            auto g_z_0_z_xxxyzzz = cbuffer.data(pk_geom_10_off + 301 * ccomps * dcomps);

            auto g_z_0_z_xxxzzzz = cbuffer.data(pk_geom_10_off + 302 * ccomps * dcomps);

            auto g_z_0_z_xxyyyyy = cbuffer.data(pk_geom_10_off + 303 * ccomps * dcomps);

            auto g_z_0_z_xxyyyyz = cbuffer.data(pk_geom_10_off + 304 * ccomps * dcomps);

            auto g_z_0_z_xxyyyzz = cbuffer.data(pk_geom_10_off + 305 * ccomps * dcomps);

            auto g_z_0_z_xxyyzzz = cbuffer.data(pk_geom_10_off + 306 * ccomps * dcomps);

            auto g_z_0_z_xxyzzzz = cbuffer.data(pk_geom_10_off + 307 * ccomps * dcomps);

            auto g_z_0_z_xxzzzzz = cbuffer.data(pk_geom_10_off + 308 * ccomps * dcomps);

            auto g_z_0_z_xyyyyyy = cbuffer.data(pk_geom_10_off + 309 * ccomps * dcomps);

            auto g_z_0_z_xyyyyyz = cbuffer.data(pk_geom_10_off + 310 * ccomps * dcomps);

            auto g_z_0_z_xyyyyzz = cbuffer.data(pk_geom_10_off + 311 * ccomps * dcomps);

            auto g_z_0_z_xyyyzzz = cbuffer.data(pk_geom_10_off + 312 * ccomps * dcomps);

            auto g_z_0_z_xyyzzzz = cbuffer.data(pk_geom_10_off + 313 * ccomps * dcomps);

            auto g_z_0_z_xyzzzzz = cbuffer.data(pk_geom_10_off + 314 * ccomps * dcomps);

            auto g_z_0_z_xzzzzzz = cbuffer.data(pk_geom_10_off + 315 * ccomps * dcomps);

            auto g_z_0_z_yyyyyyy = cbuffer.data(pk_geom_10_off + 316 * ccomps * dcomps);

            auto g_z_0_z_yyyyyyz = cbuffer.data(pk_geom_10_off + 317 * ccomps * dcomps);

            auto g_z_0_z_yyyyyzz = cbuffer.data(pk_geom_10_off + 318 * ccomps * dcomps);

            auto g_z_0_z_yyyyzzz = cbuffer.data(pk_geom_10_off + 319 * ccomps * dcomps);

            auto g_z_0_z_yyyzzzz = cbuffer.data(pk_geom_10_off + 320 * ccomps * dcomps);

            auto g_z_0_z_yyzzzzz = cbuffer.data(pk_geom_10_off + 321 * ccomps * dcomps);

            auto g_z_0_z_yzzzzzz = cbuffer.data(pk_geom_10_off + 322 * ccomps * dcomps);

            auto g_z_0_z_zzzzzzz = cbuffer.data(pk_geom_10_off + 323 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dixx

            const auto di_geom_10_off = idx_geom_10_dixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xx_yyyyyy = cbuffer.data(di_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xx_yyyyyz = cbuffer.data(di_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xx_yyyyzz = cbuffer.data(di_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xx_yyyzzz = cbuffer.data(di_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xx_yyzzzz = cbuffer.data(di_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xx_yzzzzz = cbuffer.data(di_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xx_zzzzzz = cbuffer.data(di_geom_10_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxxxxx, g_x_0_x_xxxxxxx, g_x_0_x_xxxxxxy, g_x_0_x_xxxxxxz, g_x_0_x_xxxxxy, g_x_0_x_xxxxxyy, g_x_0_x_xxxxxyz, g_x_0_x_xxxxxz, g_x_0_x_xxxxxzz, g_x_0_x_xxxxyy, g_x_0_x_xxxxyyy, g_x_0_x_xxxxyyz, g_x_0_x_xxxxyz, g_x_0_x_xxxxyzz, g_x_0_x_xxxxzz, g_x_0_x_xxxxzzz, g_x_0_x_xxxyyy, g_x_0_x_xxxyyyy, g_x_0_x_xxxyyyz, g_x_0_x_xxxyyz, g_x_0_x_xxxyyzz, g_x_0_x_xxxyzz, g_x_0_x_xxxyzzz, g_x_0_x_xxxzzz, g_x_0_x_xxxzzzz, g_x_0_x_xxyyyy, g_x_0_x_xxyyyyy, g_x_0_x_xxyyyyz, g_x_0_x_xxyyyz, g_x_0_x_xxyyyzz, g_x_0_x_xxyyzz, g_x_0_x_xxyyzzz, g_x_0_x_xxyzzz, g_x_0_x_xxyzzzz, g_x_0_x_xxzzzz, g_x_0_x_xxzzzzz, g_x_0_x_xyyyyy, g_x_0_x_xyyyyyy, g_x_0_x_xyyyyyz, g_x_0_x_xyyyyz, g_x_0_x_xyyyyzz, g_x_0_x_xyyyzz, g_x_0_x_xyyyzzz, g_x_0_x_xyyzzz, g_x_0_x_xyyzzzz, g_x_0_x_xyzzzz, g_x_0_x_xyzzzzz, g_x_0_x_xzzzzz, g_x_0_x_xzzzzzz, g_x_0_x_yyyyyy, g_x_0_x_yyyyyz, g_x_0_x_yyyyzz, g_x_0_x_yyyzzz, g_x_0_x_yyzzzz, g_x_0_x_yzzzzz, g_x_0_x_zzzzzz, g_x_0_xx_xxxxxx, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_yyyyyy, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_zzzzzz, g_x_xxxxxx, g_x_xxxxxy, g_x_xxxxxz, g_x_xxxxyy, g_x_xxxxyz, g_x_xxxxzz, g_x_xxxyyy, g_x_xxxyyz, g_x_xxxyzz, g_x_xxxzzz, g_x_xxyyyy, g_x_xxyyyz, g_x_xxyyzz, g_x_xxyzzz, g_x_xxzzzz, g_x_xyyyyy, g_x_xyyyyz, g_x_xyyyzz, g_x_xyyzzz, g_x_xyzzzz, g_x_xzzzzz, g_x_yyyyyy, g_x_yyyyyz, g_x_yyyyzz, g_x_yyyzzz, g_x_yyzzzz, g_x_yzzzzz, g_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xx_xxxxxx[k] = -g_x_xxxxxx[k] - g_x_0_x_xxxxxx[k] * ab_x + g_x_0_x_xxxxxxx[k];

                g_x_0_xx_xxxxxy[k] = -g_x_xxxxxy[k] - g_x_0_x_xxxxxy[k] * ab_x + g_x_0_x_xxxxxxy[k];

                g_x_0_xx_xxxxxz[k] = -g_x_xxxxxz[k] - g_x_0_x_xxxxxz[k] * ab_x + g_x_0_x_xxxxxxz[k];

                g_x_0_xx_xxxxyy[k] = -g_x_xxxxyy[k] - g_x_0_x_xxxxyy[k] * ab_x + g_x_0_x_xxxxxyy[k];

                g_x_0_xx_xxxxyz[k] = -g_x_xxxxyz[k] - g_x_0_x_xxxxyz[k] * ab_x + g_x_0_x_xxxxxyz[k];

                g_x_0_xx_xxxxzz[k] = -g_x_xxxxzz[k] - g_x_0_x_xxxxzz[k] * ab_x + g_x_0_x_xxxxxzz[k];

                g_x_0_xx_xxxyyy[k] = -g_x_xxxyyy[k] - g_x_0_x_xxxyyy[k] * ab_x + g_x_0_x_xxxxyyy[k];

                g_x_0_xx_xxxyyz[k] = -g_x_xxxyyz[k] - g_x_0_x_xxxyyz[k] * ab_x + g_x_0_x_xxxxyyz[k];

                g_x_0_xx_xxxyzz[k] = -g_x_xxxyzz[k] - g_x_0_x_xxxyzz[k] * ab_x + g_x_0_x_xxxxyzz[k];

                g_x_0_xx_xxxzzz[k] = -g_x_xxxzzz[k] - g_x_0_x_xxxzzz[k] * ab_x + g_x_0_x_xxxxzzz[k];

                g_x_0_xx_xxyyyy[k] = -g_x_xxyyyy[k] - g_x_0_x_xxyyyy[k] * ab_x + g_x_0_x_xxxyyyy[k];

                g_x_0_xx_xxyyyz[k] = -g_x_xxyyyz[k] - g_x_0_x_xxyyyz[k] * ab_x + g_x_0_x_xxxyyyz[k];

                g_x_0_xx_xxyyzz[k] = -g_x_xxyyzz[k] - g_x_0_x_xxyyzz[k] * ab_x + g_x_0_x_xxxyyzz[k];

                g_x_0_xx_xxyzzz[k] = -g_x_xxyzzz[k] - g_x_0_x_xxyzzz[k] * ab_x + g_x_0_x_xxxyzzz[k];

                g_x_0_xx_xxzzzz[k] = -g_x_xxzzzz[k] - g_x_0_x_xxzzzz[k] * ab_x + g_x_0_x_xxxzzzz[k];

                g_x_0_xx_xyyyyy[k] = -g_x_xyyyyy[k] - g_x_0_x_xyyyyy[k] * ab_x + g_x_0_x_xxyyyyy[k];

                g_x_0_xx_xyyyyz[k] = -g_x_xyyyyz[k] - g_x_0_x_xyyyyz[k] * ab_x + g_x_0_x_xxyyyyz[k];

                g_x_0_xx_xyyyzz[k] = -g_x_xyyyzz[k] - g_x_0_x_xyyyzz[k] * ab_x + g_x_0_x_xxyyyzz[k];

                g_x_0_xx_xyyzzz[k] = -g_x_xyyzzz[k] - g_x_0_x_xyyzzz[k] * ab_x + g_x_0_x_xxyyzzz[k];

                g_x_0_xx_xyzzzz[k] = -g_x_xyzzzz[k] - g_x_0_x_xyzzzz[k] * ab_x + g_x_0_x_xxyzzzz[k];

                g_x_0_xx_xzzzzz[k] = -g_x_xzzzzz[k] - g_x_0_x_xzzzzz[k] * ab_x + g_x_0_x_xxzzzzz[k];

                g_x_0_xx_yyyyyy[k] = -g_x_yyyyyy[k] - g_x_0_x_yyyyyy[k] * ab_x + g_x_0_x_xyyyyyy[k];

                g_x_0_xx_yyyyyz[k] = -g_x_yyyyyz[k] - g_x_0_x_yyyyyz[k] * ab_x + g_x_0_x_xyyyyyz[k];

                g_x_0_xx_yyyyzz[k] = -g_x_yyyyzz[k] - g_x_0_x_yyyyzz[k] * ab_x + g_x_0_x_xyyyyzz[k];

                g_x_0_xx_yyyzzz[k] = -g_x_yyyzzz[k] - g_x_0_x_yyyzzz[k] * ab_x + g_x_0_x_xyyyzzz[k];

                g_x_0_xx_yyzzzz[k] = -g_x_yyzzzz[k] - g_x_0_x_yyzzzz[k] * ab_x + g_x_0_x_xyyzzzz[k];

                g_x_0_xx_yzzzzz[k] = -g_x_yzzzzz[k] - g_x_0_x_yzzzzz[k] * ab_x + g_x_0_x_xyzzzzz[k];

                g_x_0_xx_zzzzzz[k] = -g_x_zzzzzz[k] - g_x_0_x_zzzzzz[k] * ab_x + g_x_0_x_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_0_xy_xxxxxx = cbuffer.data(di_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xy_xxxxxz = cbuffer.data(di_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xy_xxxxzz = cbuffer.data(di_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xy_xxxzzz = cbuffer.data(di_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xy_xxzzzz = cbuffer.data(di_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xy_xzzzzz = cbuffer.data(di_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xy_yyyyyy = cbuffer.data(di_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xy_yyyyyz = cbuffer.data(di_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xy_yyyyzz = cbuffer.data(di_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xy_yyyzzz = cbuffer.data(di_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xy_yyzzzz = cbuffer.data(di_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xy_yzzzzz = cbuffer.data(di_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xy_zzzzzz = cbuffer.data(di_geom_10_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxxxxx, g_x_0_x_xxxxxxy, g_x_0_x_xxxxxy, g_x_0_x_xxxxxyy, g_x_0_x_xxxxxyz, g_x_0_x_xxxxxz, g_x_0_x_xxxxyy, g_x_0_x_xxxxyyy, g_x_0_x_xxxxyyz, g_x_0_x_xxxxyz, g_x_0_x_xxxxyzz, g_x_0_x_xxxxzz, g_x_0_x_xxxyyy, g_x_0_x_xxxyyyy, g_x_0_x_xxxyyyz, g_x_0_x_xxxyyz, g_x_0_x_xxxyyzz, g_x_0_x_xxxyzz, g_x_0_x_xxxyzzz, g_x_0_x_xxxzzz, g_x_0_x_xxyyyy, g_x_0_x_xxyyyyy, g_x_0_x_xxyyyyz, g_x_0_x_xxyyyz, g_x_0_x_xxyyyzz, g_x_0_x_xxyyzz, g_x_0_x_xxyyzzz, g_x_0_x_xxyzzz, g_x_0_x_xxyzzzz, g_x_0_x_xxzzzz, g_x_0_x_xyyyyy, g_x_0_x_xyyyyyy, g_x_0_x_xyyyyyz, g_x_0_x_xyyyyz, g_x_0_x_xyyyyzz, g_x_0_x_xyyyzz, g_x_0_x_xyyyzzz, g_x_0_x_xyyzzz, g_x_0_x_xyyzzzz, g_x_0_x_xyzzzz, g_x_0_x_xyzzzzz, g_x_0_x_xzzzzz, g_x_0_x_yyyyyy, g_x_0_x_yyyyyyy, g_x_0_x_yyyyyyz, g_x_0_x_yyyyyz, g_x_0_x_yyyyyzz, g_x_0_x_yyyyzz, g_x_0_x_yyyyzzz, g_x_0_x_yyyzzz, g_x_0_x_yyyzzzz, g_x_0_x_yyzzzz, g_x_0_x_yyzzzzz, g_x_0_x_yzzzzz, g_x_0_x_yzzzzzz, g_x_0_x_zzzzzz, g_x_0_xy_xxxxxx, g_x_0_xy_xxxxxy, g_x_0_xy_xxxxxz, g_x_0_xy_xxxxyy, g_x_0_xy_xxxxyz, g_x_0_xy_xxxxzz, g_x_0_xy_xxxyyy, g_x_0_xy_xxxyyz, g_x_0_xy_xxxyzz, g_x_0_xy_xxxzzz, g_x_0_xy_xxyyyy, g_x_0_xy_xxyyyz, g_x_0_xy_xxyyzz, g_x_0_xy_xxyzzz, g_x_0_xy_xxzzzz, g_x_0_xy_xyyyyy, g_x_0_xy_xyyyyz, g_x_0_xy_xyyyzz, g_x_0_xy_xyyzzz, g_x_0_xy_xyzzzz, g_x_0_xy_xzzzzz, g_x_0_xy_yyyyyy, g_x_0_xy_yyyyyz, g_x_0_xy_yyyyzz, g_x_0_xy_yyyzzz, g_x_0_xy_yyzzzz, g_x_0_xy_yzzzzz, g_x_0_xy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xy_xxxxxx[k] = -g_x_0_x_xxxxxx[k] * ab_y + g_x_0_x_xxxxxxy[k];

                g_x_0_xy_xxxxxy[k] = -g_x_0_x_xxxxxy[k] * ab_y + g_x_0_x_xxxxxyy[k];

                g_x_0_xy_xxxxxz[k] = -g_x_0_x_xxxxxz[k] * ab_y + g_x_0_x_xxxxxyz[k];

                g_x_0_xy_xxxxyy[k] = -g_x_0_x_xxxxyy[k] * ab_y + g_x_0_x_xxxxyyy[k];

                g_x_0_xy_xxxxyz[k] = -g_x_0_x_xxxxyz[k] * ab_y + g_x_0_x_xxxxyyz[k];

                g_x_0_xy_xxxxzz[k] = -g_x_0_x_xxxxzz[k] * ab_y + g_x_0_x_xxxxyzz[k];

                g_x_0_xy_xxxyyy[k] = -g_x_0_x_xxxyyy[k] * ab_y + g_x_0_x_xxxyyyy[k];

                g_x_0_xy_xxxyyz[k] = -g_x_0_x_xxxyyz[k] * ab_y + g_x_0_x_xxxyyyz[k];

                g_x_0_xy_xxxyzz[k] = -g_x_0_x_xxxyzz[k] * ab_y + g_x_0_x_xxxyyzz[k];

                g_x_0_xy_xxxzzz[k] = -g_x_0_x_xxxzzz[k] * ab_y + g_x_0_x_xxxyzzz[k];

                g_x_0_xy_xxyyyy[k] = -g_x_0_x_xxyyyy[k] * ab_y + g_x_0_x_xxyyyyy[k];

                g_x_0_xy_xxyyyz[k] = -g_x_0_x_xxyyyz[k] * ab_y + g_x_0_x_xxyyyyz[k];

                g_x_0_xy_xxyyzz[k] = -g_x_0_x_xxyyzz[k] * ab_y + g_x_0_x_xxyyyzz[k];

                g_x_0_xy_xxyzzz[k] = -g_x_0_x_xxyzzz[k] * ab_y + g_x_0_x_xxyyzzz[k];

                g_x_0_xy_xxzzzz[k] = -g_x_0_x_xxzzzz[k] * ab_y + g_x_0_x_xxyzzzz[k];

                g_x_0_xy_xyyyyy[k] = -g_x_0_x_xyyyyy[k] * ab_y + g_x_0_x_xyyyyyy[k];

                g_x_0_xy_xyyyyz[k] = -g_x_0_x_xyyyyz[k] * ab_y + g_x_0_x_xyyyyyz[k];

                g_x_0_xy_xyyyzz[k] = -g_x_0_x_xyyyzz[k] * ab_y + g_x_0_x_xyyyyzz[k];

                g_x_0_xy_xyyzzz[k] = -g_x_0_x_xyyzzz[k] * ab_y + g_x_0_x_xyyyzzz[k];

                g_x_0_xy_xyzzzz[k] = -g_x_0_x_xyzzzz[k] * ab_y + g_x_0_x_xyyzzzz[k];

                g_x_0_xy_xzzzzz[k] = -g_x_0_x_xzzzzz[k] * ab_y + g_x_0_x_xyzzzzz[k];

                g_x_0_xy_yyyyyy[k] = -g_x_0_x_yyyyyy[k] * ab_y + g_x_0_x_yyyyyyy[k];

                g_x_0_xy_yyyyyz[k] = -g_x_0_x_yyyyyz[k] * ab_y + g_x_0_x_yyyyyyz[k];

                g_x_0_xy_yyyyzz[k] = -g_x_0_x_yyyyzz[k] * ab_y + g_x_0_x_yyyyyzz[k];

                g_x_0_xy_yyyzzz[k] = -g_x_0_x_yyyzzz[k] * ab_y + g_x_0_x_yyyyzzz[k];

                g_x_0_xy_yyzzzz[k] = -g_x_0_x_yyzzzz[k] * ab_y + g_x_0_x_yyyzzzz[k];

                g_x_0_xy_yzzzzz[k] = -g_x_0_x_yzzzzz[k] * ab_y + g_x_0_x_yyzzzzz[k];

                g_x_0_xy_zzzzzz[k] = -g_x_0_x_zzzzzz[k] * ab_y + g_x_0_x_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xz_xxxxxx = cbuffer.data(di_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xz_yyyyyy = cbuffer.data(di_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xz_yyyyyz = cbuffer.data(di_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xz_yyyyzz = cbuffer.data(di_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xz_yyyzzz = cbuffer.data(di_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xz_yyzzzz = cbuffer.data(di_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xz_yzzzzz = cbuffer.data(di_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xz_zzzzzz = cbuffer.data(di_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxxxxx, g_x_0_x_xxxxxxz, g_x_0_x_xxxxxy, g_x_0_x_xxxxxyz, g_x_0_x_xxxxxz, g_x_0_x_xxxxxzz, g_x_0_x_xxxxyy, g_x_0_x_xxxxyyz, g_x_0_x_xxxxyz, g_x_0_x_xxxxyzz, g_x_0_x_xxxxzz, g_x_0_x_xxxxzzz, g_x_0_x_xxxyyy, g_x_0_x_xxxyyyz, g_x_0_x_xxxyyz, g_x_0_x_xxxyyzz, g_x_0_x_xxxyzz, g_x_0_x_xxxyzzz, g_x_0_x_xxxzzz, g_x_0_x_xxxzzzz, g_x_0_x_xxyyyy, g_x_0_x_xxyyyyz, g_x_0_x_xxyyyz, g_x_0_x_xxyyyzz, g_x_0_x_xxyyzz, g_x_0_x_xxyyzzz, g_x_0_x_xxyzzz, g_x_0_x_xxyzzzz, g_x_0_x_xxzzzz, g_x_0_x_xxzzzzz, g_x_0_x_xyyyyy, g_x_0_x_xyyyyyz, g_x_0_x_xyyyyz, g_x_0_x_xyyyyzz, g_x_0_x_xyyyzz, g_x_0_x_xyyyzzz, g_x_0_x_xyyzzz, g_x_0_x_xyyzzzz, g_x_0_x_xyzzzz, g_x_0_x_xyzzzzz, g_x_0_x_xzzzzz, g_x_0_x_xzzzzzz, g_x_0_x_yyyyyy, g_x_0_x_yyyyyyz, g_x_0_x_yyyyyz, g_x_0_x_yyyyyzz, g_x_0_x_yyyyzz, g_x_0_x_yyyyzzz, g_x_0_x_yyyzzz, g_x_0_x_yyyzzzz, g_x_0_x_yyzzzz, g_x_0_x_yyzzzzz, g_x_0_x_yzzzzz, g_x_0_x_yzzzzzz, g_x_0_x_zzzzzz, g_x_0_x_zzzzzzz, g_x_0_xz_xxxxxx, g_x_0_xz_xxxxxy, g_x_0_xz_xxxxxz, g_x_0_xz_xxxxyy, g_x_0_xz_xxxxyz, g_x_0_xz_xxxxzz, g_x_0_xz_xxxyyy, g_x_0_xz_xxxyyz, g_x_0_xz_xxxyzz, g_x_0_xz_xxxzzz, g_x_0_xz_xxyyyy, g_x_0_xz_xxyyyz, g_x_0_xz_xxyyzz, g_x_0_xz_xxyzzz, g_x_0_xz_xxzzzz, g_x_0_xz_xyyyyy, g_x_0_xz_xyyyyz, g_x_0_xz_xyyyzz, g_x_0_xz_xyyzzz, g_x_0_xz_xyzzzz, g_x_0_xz_xzzzzz, g_x_0_xz_yyyyyy, g_x_0_xz_yyyyyz, g_x_0_xz_yyyyzz, g_x_0_xz_yyyzzz, g_x_0_xz_yyzzzz, g_x_0_xz_yzzzzz, g_x_0_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xz_xxxxxx[k] = -g_x_0_x_xxxxxx[k] * ab_z + g_x_0_x_xxxxxxz[k];

                g_x_0_xz_xxxxxy[k] = -g_x_0_x_xxxxxy[k] * ab_z + g_x_0_x_xxxxxyz[k];

                g_x_0_xz_xxxxxz[k] = -g_x_0_x_xxxxxz[k] * ab_z + g_x_0_x_xxxxxzz[k];

                g_x_0_xz_xxxxyy[k] = -g_x_0_x_xxxxyy[k] * ab_z + g_x_0_x_xxxxyyz[k];

                g_x_0_xz_xxxxyz[k] = -g_x_0_x_xxxxyz[k] * ab_z + g_x_0_x_xxxxyzz[k];

                g_x_0_xz_xxxxzz[k] = -g_x_0_x_xxxxzz[k] * ab_z + g_x_0_x_xxxxzzz[k];

                g_x_0_xz_xxxyyy[k] = -g_x_0_x_xxxyyy[k] * ab_z + g_x_0_x_xxxyyyz[k];

                g_x_0_xz_xxxyyz[k] = -g_x_0_x_xxxyyz[k] * ab_z + g_x_0_x_xxxyyzz[k];

                g_x_0_xz_xxxyzz[k] = -g_x_0_x_xxxyzz[k] * ab_z + g_x_0_x_xxxyzzz[k];

                g_x_0_xz_xxxzzz[k] = -g_x_0_x_xxxzzz[k] * ab_z + g_x_0_x_xxxzzzz[k];

                g_x_0_xz_xxyyyy[k] = -g_x_0_x_xxyyyy[k] * ab_z + g_x_0_x_xxyyyyz[k];

                g_x_0_xz_xxyyyz[k] = -g_x_0_x_xxyyyz[k] * ab_z + g_x_0_x_xxyyyzz[k];

                g_x_0_xz_xxyyzz[k] = -g_x_0_x_xxyyzz[k] * ab_z + g_x_0_x_xxyyzzz[k];

                g_x_0_xz_xxyzzz[k] = -g_x_0_x_xxyzzz[k] * ab_z + g_x_0_x_xxyzzzz[k];

                g_x_0_xz_xxzzzz[k] = -g_x_0_x_xxzzzz[k] * ab_z + g_x_0_x_xxzzzzz[k];

                g_x_0_xz_xyyyyy[k] = -g_x_0_x_xyyyyy[k] * ab_z + g_x_0_x_xyyyyyz[k];

                g_x_0_xz_xyyyyz[k] = -g_x_0_x_xyyyyz[k] * ab_z + g_x_0_x_xyyyyzz[k];

                g_x_0_xz_xyyyzz[k] = -g_x_0_x_xyyyzz[k] * ab_z + g_x_0_x_xyyyzzz[k];

                g_x_0_xz_xyyzzz[k] = -g_x_0_x_xyyzzz[k] * ab_z + g_x_0_x_xyyzzzz[k];

                g_x_0_xz_xyzzzz[k] = -g_x_0_x_xyzzzz[k] * ab_z + g_x_0_x_xyzzzzz[k];

                g_x_0_xz_xzzzzz[k] = -g_x_0_x_xzzzzz[k] * ab_z + g_x_0_x_xzzzzzz[k];

                g_x_0_xz_yyyyyy[k] = -g_x_0_x_yyyyyy[k] * ab_z + g_x_0_x_yyyyyyz[k];

                g_x_0_xz_yyyyyz[k] = -g_x_0_x_yyyyyz[k] * ab_z + g_x_0_x_yyyyyzz[k];

                g_x_0_xz_yyyyzz[k] = -g_x_0_x_yyyyzz[k] * ab_z + g_x_0_x_yyyyzzz[k];

                g_x_0_xz_yyyzzz[k] = -g_x_0_x_yyyzzz[k] * ab_z + g_x_0_x_yyyzzzz[k];

                g_x_0_xz_yyzzzz[k] = -g_x_0_x_yyzzzz[k] * ab_z + g_x_0_x_yyzzzzz[k];

                g_x_0_xz_yzzzzz[k] = -g_x_0_x_yzzzzz[k] * ab_z + g_x_0_x_yzzzzzz[k];

                g_x_0_xz_zzzzzz[k] = -g_x_0_x_zzzzzz[k] * ab_z + g_x_0_x_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_x_0_yy_xxxxxx = cbuffer.data(di_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_yy_xxxxxz = cbuffer.data(di_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_yy_xxxxzz = cbuffer.data(di_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_yy_xxxzzz = cbuffer.data(di_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_yy_xxzzzz = cbuffer.data(di_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_yy_xzzzzz = cbuffer.data(di_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_yy_zzzzzz = cbuffer.data(di_geom_10_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxxxxx, g_x_0_y_xxxxxxy, g_x_0_y_xxxxxy, g_x_0_y_xxxxxyy, g_x_0_y_xxxxxyz, g_x_0_y_xxxxxz, g_x_0_y_xxxxyy, g_x_0_y_xxxxyyy, g_x_0_y_xxxxyyz, g_x_0_y_xxxxyz, g_x_0_y_xxxxyzz, g_x_0_y_xxxxzz, g_x_0_y_xxxyyy, g_x_0_y_xxxyyyy, g_x_0_y_xxxyyyz, g_x_0_y_xxxyyz, g_x_0_y_xxxyyzz, g_x_0_y_xxxyzz, g_x_0_y_xxxyzzz, g_x_0_y_xxxzzz, g_x_0_y_xxyyyy, g_x_0_y_xxyyyyy, g_x_0_y_xxyyyyz, g_x_0_y_xxyyyz, g_x_0_y_xxyyyzz, g_x_0_y_xxyyzz, g_x_0_y_xxyyzzz, g_x_0_y_xxyzzz, g_x_0_y_xxyzzzz, g_x_0_y_xxzzzz, g_x_0_y_xyyyyy, g_x_0_y_xyyyyyy, g_x_0_y_xyyyyyz, g_x_0_y_xyyyyz, g_x_0_y_xyyyyzz, g_x_0_y_xyyyzz, g_x_0_y_xyyyzzz, g_x_0_y_xyyzzz, g_x_0_y_xyyzzzz, g_x_0_y_xyzzzz, g_x_0_y_xyzzzzz, g_x_0_y_xzzzzz, g_x_0_y_yyyyyy, g_x_0_y_yyyyyyy, g_x_0_y_yyyyyyz, g_x_0_y_yyyyyz, g_x_0_y_yyyyyzz, g_x_0_y_yyyyzz, g_x_0_y_yyyyzzz, g_x_0_y_yyyzzz, g_x_0_y_yyyzzzz, g_x_0_y_yyzzzz, g_x_0_y_yyzzzzz, g_x_0_y_yzzzzz, g_x_0_y_yzzzzzz, g_x_0_y_zzzzzz, g_x_0_yy_xxxxxx, g_x_0_yy_xxxxxy, g_x_0_yy_xxxxxz, g_x_0_yy_xxxxyy, g_x_0_yy_xxxxyz, g_x_0_yy_xxxxzz, g_x_0_yy_xxxyyy, g_x_0_yy_xxxyyz, g_x_0_yy_xxxyzz, g_x_0_yy_xxxzzz, g_x_0_yy_xxyyyy, g_x_0_yy_xxyyyz, g_x_0_yy_xxyyzz, g_x_0_yy_xxyzzz, g_x_0_yy_xxzzzz, g_x_0_yy_xyyyyy, g_x_0_yy_xyyyyz, g_x_0_yy_xyyyzz, g_x_0_yy_xyyzzz, g_x_0_yy_xyzzzz, g_x_0_yy_xzzzzz, g_x_0_yy_yyyyyy, g_x_0_yy_yyyyyz, g_x_0_yy_yyyyzz, g_x_0_yy_yyyzzz, g_x_0_yy_yyzzzz, g_x_0_yy_yzzzzz, g_x_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yy_xxxxxx[k] = -g_x_0_y_xxxxxx[k] * ab_y + g_x_0_y_xxxxxxy[k];

                g_x_0_yy_xxxxxy[k] = -g_x_0_y_xxxxxy[k] * ab_y + g_x_0_y_xxxxxyy[k];

                g_x_0_yy_xxxxxz[k] = -g_x_0_y_xxxxxz[k] * ab_y + g_x_0_y_xxxxxyz[k];

                g_x_0_yy_xxxxyy[k] = -g_x_0_y_xxxxyy[k] * ab_y + g_x_0_y_xxxxyyy[k];

                g_x_0_yy_xxxxyz[k] = -g_x_0_y_xxxxyz[k] * ab_y + g_x_0_y_xxxxyyz[k];

                g_x_0_yy_xxxxzz[k] = -g_x_0_y_xxxxzz[k] * ab_y + g_x_0_y_xxxxyzz[k];

                g_x_0_yy_xxxyyy[k] = -g_x_0_y_xxxyyy[k] * ab_y + g_x_0_y_xxxyyyy[k];

                g_x_0_yy_xxxyyz[k] = -g_x_0_y_xxxyyz[k] * ab_y + g_x_0_y_xxxyyyz[k];

                g_x_0_yy_xxxyzz[k] = -g_x_0_y_xxxyzz[k] * ab_y + g_x_0_y_xxxyyzz[k];

                g_x_0_yy_xxxzzz[k] = -g_x_0_y_xxxzzz[k] * ab_y + g_x_0_y_xxxyzzz[k];

                g_x_0_yy_xxyyyy[k] = -g_x_0_y_xxyyyy[k] * ab_y + g_x_0_y_xxyyyyy[k];

                g_x_0_yy_xxyyyz[k] = -g_x_0_y_xxyyyz[k] * ab_y + g_x_0_y_xxyyyyz[k];

                g_x_0_yy_xxyyzz[k] = -g_x_0_y_xxyyzz[k] * ab_y + g_x_0_y_xxyyyzz[k];

                g_x_0_yy_xxyzzz[k] = -g_x_0_y_xxyzzz[k] * ab_y + g_x_0_y_xxyyzzz[k];

                g_x_0_yy_xxzzzz[k] = -g_x_0_y_xxzzzz[k] * ab_y + g_x_0_y_xxyzzzz[k];

                g_x_0_yy_xyyyyy[k] = -g_x_0_y_xyyyyy[k] * ab_y + g_x_0_y_xyyyyyy[k];

                g_x_0_yy_xyyyyz[k] = -g_x_0_y_xyyyyz[k] * ab_y + g_x_0_y_xyyyyyz[k];

                g_x_0_yy_xyyyzz[k] = -g_x_0_y_xyyyzz[k] * ab_y + g_x_0_y_xyyyyzz[k];

                g_x_0_yy_xyyzzz[k] = -g_x_0_y_xyyzzz[k] * ab_y + g_x_0_y_xyyyzzz[k];

                g_x_0_yy_xyzzzz[k] = -g_x_0_y_xyzzzz[k] * ab_y + g_x_0_y_xyyzzzz[k];

                g_x_0_yy_xzzzzz[k] = -g_x_0_y_xzzzzz[k] * ab_y + g_x_0_y_xyzzzzz[k];

                g_x_0_yy_yyyyyy[k] = -g_x_0_y_yyyyyy[k] * ab_y + g_x_0_y_yyyyyyy[k];

                g_x_0_yy_yyyyyz[k] = -g_x_0_y_yyyyyz[k] * ab_y + g_x_0_y_yyyyyyz[k];

                g_x_0_yy_yyyyzz[k] = -g_x_0_y_yyyyzz[k] * ab_y + g_x_0_y_yyyyyzz[k];

                g_x_0_yy_yyyzzz[k] = -g_x_0_y_yyyzzz[k] * ab_y + g_x_0_y_yyyyzzz[k];

                g_x_0_yy_yyzzzz[k] = -g_x_0_y_yyzzzz[k] * ab_y + g_x_0_y_yyyzzzz[k];

                g_x_0_yy_yzzzzz[k] = -g_x_0_y_yzzzzz[k] * ab_y + g_x_0_y_yyzzzzz[k];

                g_x_0_yy_zzzzzz[k] = -g_x_0_y_zzzzzz[k] * ab_y + g_x_0_y_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_yz_xxxxxx = cbuffer.data(di_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_yz_xxxxxz = cbuffer.data(di_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_yz_xxxxzz = cbuffer.data(di_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_yz_xxxzzz = cbuffer.data(di_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_yz_xxzzzz = cbuffer.data(di_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_yz_xzzzzz = cbuffer.data(di_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_yz_yyyyyy = cbuffer.data(di_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_yz_zzzzzz = cbuffer.data(di_geom_10_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yz_xxxxxx, g_x_0_yz_xxxxxy, g_x_0_yz_xxxxxz, g_x_0_yz_xxxxyy, g_x_0_yz_xxxxyz, g_x_0_yz_xxxxzz, g_x_0_yz_xxxyyy, g_x_0_yz_xxxyyz, g_x_0_yz_xxxyzz, g_x_0_yz_xxxzzz, g_x_0_yz_xxyyyy, g_x_0_yz_xxyyyz, g_x_0_yz_xxyyzz, g_x_0_yz_xxyzzz, g_x_0_yz_xxzzzz, g_x_0_yz_xyyyyy, g_x_0_yz_xyyyyz, g_x_0_yz_xyyyzz, g_x_0_yz_xyyzzz, g_x_0_yz_xyzzzz, g_x_0_yz_xzzzzz, g_x_0_yz_yyyyyy, g_x_0_yz_yyyyyz, g_x_0_yz_yyyyzz, g_x_0_yz_yyyzzz, g_x_0_yz_yyzzzz, g_x_0_yz_yzzzzz, g_x_0_yz_zzzzzz, g_x_0_z_xxxxxx, g_x_0_z_xxxxxxy, g_x_0_z_xxxxxy, g_x_0_z_xxxxxyy, g_x_0_z_xxxxxyz, g_x_0_z_xxxxxz, g_x_0_z_xxxxyy, g_x_0_z_xxxxyyy, g_x_0_z_xxxxyyz, g_x_0_z_xxxxyz, g_x_0_z_xxxxyzz, g_x_0_z_xxxxzz, g_x_0_z_xxxyyy, g_x_0_z_xxxyyyy, g_x_0_z_xxxyyyz, g_x_0_z_xxxyyz, g_x_0_z_xxxyyzz, g_x_0_z_xxxyzz, g_x_0_z_xxxyzzz, g_x_0_z_xxxzzz, g_x_0_z_xxyyyy, g_x_0_z_xxyyyyy, g_x_0_z_xxyyyyz, g_x_0_z_xxyyyz, g_x_0_z_xxyyyzz, g_x_0_z_xxyyzz, g_x_0_z_xxyyzzz, g_x_0_z_xxyzzz, g_x_0_z_xxyzzzz, g_x_0_z_xxzzzz, g_x_0_z_xyyyyy, g_x_0_z_xyyyyyy, g_x_0_z_xyyyyyz, g_x_0_z_xyyyyz, g_x_0_z_xyyyyzz, g_x_0_z_xyyyzz, g_x_0_z_xyyyzzz, g_x_0_z_xyyzzz, g_x_0_z_xyyzzzz, g_x_0_z_xyzzzz, g_x_0_z_xyzzzzz, g_x_0_z_xzzzzz, g_x_0_z_yyyyyy, g_x_0_z_yyyyyyy, g_x_0_z_yyyyyyz, g_x_0_z_yyyyyz, g_x_0_z_yyyyyzz, g_x_0_z_yyyyzz, g_x_0_z_yyyyzzz, g_x_0_z_yyyzzz, g_x_0_z_yyyzzzz, g_x_0_z_yyzzzz, g_x_0_z_yyzzzzz, g_x_0_z_yzzzzz, g_x_0_z_yzzzzzz, g_x_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yz_xxxxxx[k] = -g_x_0_z_xxxxxx[k] * ab_y + g_x_0_z_xxxxxxy[k];

                g_x_0_yz_xxxxxy[k] = -g_x_0_z_xxxxxy[k] * ab_y + g_x_0_z_xxxxxyy[k];

                g_x_0_yz_xxxxxz[k] = -g_x_0_z_xxxxxz[k] * ab_y + g_x_0_z_xxxxxyz[k];

                g_x_0_yz_xxxxyy[k] = -g_x_0_z_xxxxyy[k] * ab_y + g_x_0_z_xxxxyyy[k];

                g_x_0_yz_xxxxyz[k] = -g_x_0_z_xxxxyz[k] * ab_y + g_x_0_z_xxxxyyz[k];

                g_x_0_yz_xxxxzz[k] = -g_x_0_z_xxxxzz[k] * ab_y + g_x_0_z_xxxxyzz[k];

                g_x_0_yz_xxxyyy[k] = -g_x_0_z_xxxyyy[k] * ab_y + g_x_0_z_xxxyyyy[k];

                g_x_0_yz_xxxyyz[k] = -g_x_0_z_xxxyyz[k] * ab_y + g_x_0_z_xxxyyyz[k];

                g_x_0_yz_xxxyzz[k] = -g_x_0_z_xxxyzz[k] * ab_y + g_x_0_z_xxxyyzz[k];

                g_x_0_yz_xxxzzz[k] = -g_x_0_z_xxxzzz[k] * ab_y + g_x_0_z_xxxyzzz[k];

                g_x_0_yz_xxyyyy[k] = -g_x_0_z_xxyyyy[k] * ab_y + g_x_0_z_xxyyyyy[k];

                g_x_0_yz_xxyyyz[k] = -g_x_0_z_xxyyyz[k] * ab_y + g_x_0_z_xxyyyyz[k];

                g_x_0_yz_xxyyzz[k] = -g_x_0_z_xxyyzz[k] * ab_y + g_x_0_z_xxyyyzz[k];

                g_x_0_yz_xxyzzz[k] = -g_x_0_z_xxyzzz[k] * ab_y + g_x_0_z_xxyyzzz[k];

                g_x_0_yz_xxzzzz[k] = -g_x_0_z_xxzzzz[k] * ab_y + g_x_0_z_xxyzzzz[k];

                g_x_0_yz_xyyyyy[k] = -g_x_0_z_xyyyyy[k] * ab_y + g_x_0_z_xyyyyyy[k];

                g_x_0_yz_xyyyyz[k] = -g_x_0_z_xyyyyz[k] * ab_y + g_x_0_z_xyyyyyz[k];

                g_x_0_yz_xyyyzz[k] = -g_x_0_z_xyyyzz[k] * ab_y + g_x_0_z_xyyyyzz[k];

                g_x_0_yz_xyyzzz[k] = -g_x_0_z_xyyzzz[k] * ab_y + g_x_0_z_xyyyzzz[k];

                g_x_0_yz_xyzzzz[k] = -g_x_0_z_xyzzzz[k] * ab_y + g_x_0_z_xyyzzzz[k];

                g_x_0_yz_xzzzzz[k] = -g_x_0_z_xzzzzz[k] * ab_y + g_x_0_z_xyzzzzz[k];

                g_x_0_yz_yyyyyy[k] = -g_x_0_z_yyyyyy[k] * ab_y + g_x_0_z_yyyyyyy[k];

                g_x_0_yz_yyyyyz[k] = -g_x_0_z_yyyyyz[k] * ab_y + g_x_0_z_yyyyyyz[k];

                g_x_0_yz_yyyyzz[k] = -g_x_0_z_yyyyzz[k] * ab_y + g_x_0_z_yyyyyzz[k];

                g_x_0_yz_yyyzzz[k] = -g_x_0_z_yyyzzz[k] * ab_y + g_x_0_z_yyyyzzz[k];

                g_x_0_yz_yyzzzz[k] = -g_x_0_z_yyzzzz[k] * ab_y + g_x_0_z_yyyzzzz[k];

                g_x_0_yz_yzzzzz[k] = -g_x_0_z_yzzzzz[k] * ab_y + g_x_0_z_yyzzzzz[k];

                g_x_0_yz_zzzzzz[k] = -g_x_0_z_zzzzzz[k] * ab_y + g_x_0_z_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_zz_xxxxxx = cbuffer.data(di_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_zz_yyyyyy = cbuffer.data(di_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxxxxx, g_x_0_z_xxxxxxz, g_x_0_z_xxxxxy, g_x_0_z_xxxxxyz, g_x_0_z_xxxxxz, g_x_0_z_xxxxxzz, g_x_0_z_xxxxyy, g_x_0_z_xxxxyyz, g_x_0_z_xxxxyz, g_x_0_z_xxxxyzz, g_x_0_z_xxxxzz, g_x_0_z_xxxxzzz, g_x_0_z_xxxyyy, g_x_0_z_xxxyyyz, g_x_0_z_xxxyyz, g_x_0_z_xxxyyzz, g_x_0_z_xxxyzz, g_x_0_z_xxxyzzz, g_x_0_z_xxxzzz, g_x_0_z_xxxzzzz, g_x_0_z_xxyyyy, g_x_0_z_xxyyyyz, g_x_0_z_xxyyyz, g_x_0_z_xxyyyzz, g_x_0_z_xxyyzz, g_x_0_z_xxyyzzz, g_x_0_z_xxyzzz, g_x_0_z_xxyzzzz, g_x_0_z_xxzzzz, g_x_0_z_xxzzzzz, g_x_0_z_xyyyyy, g_x_0_z_xyyyyyz, g_x_0_z_xyyyyz, g_x_0_z_xyyyyzz, g_x_0_z_xyyyzz, g_x_0_z_xyyyzzz, g_x_0_z_xyyzzz, g_x_0_z_xyyzzzz, g_x_0_z_xyzzzz, g_x_0_z_xyzzzzz, g_x_0_z_xzzzzz, g_x_0_z_xzzzzzz, g_x_0_z_yyyyyy, g_x_0_z_yyyyyyz, g_x_0_z_yyyyyz, g_x_0_z_yyyyyzz, g_x_0_z_yyyyzz, g_x_0_z_yyyyzzz, g_x_0_z_yyyzzz, g_x_0_z_yyyzzzz, g_x_0_z_yyzzzz, g_x_0_z_yyzzzzz, g_x_0_z_yzzzzz, g_x_0_z_yzzzzzz, g_x_0_z_zzzzzz, g_x_0_z_zzzzzzz, g_x_0_zz_xxxxxx, g_x_0_zz_xxxxxy, g_x_0_zz_xxxxxz, g_x_0_zz_xxxxyy, g_x_0_zz_xxxxyz, g_x_0_zz_xxxxzz, g_x_0_zz_xxxyyy, g_x_0_zz_xxxyyz, g_x_0_zz_xxxyzz, g_x_0_zz_xxxzzz, g_x_0_zz_xxyyyy, g_x_0_zz_xxyyyz, g_x_0_zz_xxyyzz, g_x_0_zz_xxyzzz, g_x_0_zz_xxzzzz, g_x_0_zz_xyyyyy, g_x_0_zz_xyyyyz, g_x_0_zz_xyyyzz, g_x_0_zz_xyyzzz, g_x_0_zz_xyzzzz, g_x_0_zz_xzzzzz, g_x_0_zz_yyyyyy, g_x_0_zz_yyyyyz, g_x_0_zz_yyyyzz, g_x_0_zz_yyyzzz, g_x_0_zz_yyzzzz, g_x_0_zz_yzzzzz, g_x_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zz_xxxxxx[k] = -g_x_0_z_xxxxxx[k] * ab_z + g_x_0_z_xxxxxxz[k];

                g_x_0_zz_xxxxxy[k] = -g_x_0_z_xxxxxy[k] * ab_z + g_x_0_z_xxxxxyz[k];

                g_x_0_zz_xxxxxz[k] = -g_x_0_z_xxxxxz[k] * ab_z + g_x_0_z_xxxxxzz[k];

                g_x_0_zz_xxxxyy[k] = -g_x_0_z_xxxxyy[k] * ab_z + g_x_0_z_xxxxyyz[k];

                g_x_0_zz_xxxxyz[k] = -g_x_0_z_xxxxyz[k] * ab_z + g_x_0_z_xxxxyzz[k];

                g_x_0_zz_xxxxzz[k] = -g_x_0_z_xxxxzz[k] * ab_z + g_x_0_z_xxxxzzz[k];

                g_x_0_zz_xxxyyy[k] = -g_x_0_z_xxxyyy[k] * ab_z + g_x_0_z_xxxyyyz[k];

                g_x_0_zz_xxxyyz[k] = -g_x_0_z_xxxyyz[k] * ab_z + g_x_0_z_xxxyyzz[k];

                g_x_0_zz_xxxyzz[k] = -g_x_0_z_xxxyzz[k] * ab_z + g_x_0_z_xxxyzzz[k];

                g_x_0_zz_xxxzzz[k] = -g_x_0_z_xxxzzz[k] * ab_z + g_x_0_z_xxxzzzz[k];

                g_x_0_zz_xxyyyy[k] = -g_x_0_z_xxyyyy[k] * ab_z + g_x_0_z_xxyyyyz[k];

                g_x_0_zz_xxyyyz[k] = -g_x_0_z_xxyyyz[k] * ab_z + g_x_0_z_xxyyyzz[k];

                g_x_0_zz_xxyyzz[k] = -g_x_0_z_xxyyzz[k] * ab_z + g_x_0_z_xxyyzzz[k];

                g_x_0_zz_xxyzzz[k] = -g_x_0_z_xxyzzz[k] * ab_z + g_x_0_z_xxyzzzz[k];

                g_x_0_zz_xxzzzz[k] = -g_x_0_z_xxzzzz[k] * ab_z + g_x_0_z_xxzzzzz[k];

                g_x_0_zz_xyyyyy[k] = -g_x_0_z_xyyyyy[k] * ab_z + g_x_0_z_xyyyyyz[k];

                g_x_0_zz_xyyyyz[k] = -g_x_0_z_xyyyyz[k] * ab_z + g_x_0_z_xyyyyzz[k];

                g_x_0_zz_xyyyzz[k] = -g_x_0_z_xyyyzz[k] * ab_z + g_x_0_z_xyyyzzz[k];

                g_x_0_zz_xyyzzz[k] = -g_x_0_z_xyyzzz[k] * ab_z + g_x_0_z_xyyzzzz[k];

                g_x_0_zz_xyzzzz[k] = -g_x_0_z_xyzzzz[k] * ab_z + g_x_0_z_xyzzzzz[k];

                g_x_0_zz_xzzzzz[k] = -g_x_0_z_xzzzzz[k] * ab_z + g_x_0_z_xzzzzzz[k];

                g_x_0_zz_yyyyyy[k] = -g_x_0_z_yyyyyy[k] * ab_z + g_x_0_z_yyyyyyz[k];

                g_x_0_zz_yyyyyz[k] = -g_x_0_z_yyyyyz[k] * ab_z + g_x_0_z_yyyyyzz[k];

                g_x_0_zz_yyyyzz[k] = -g_x_0_z_yyyyzz[k] * ab_z + g_x_0_z_yyyyzzz[k];

                g_x_0_zz_yyyzzz[k] = -g_x_0_z_yyyzzz[k] * ab_z + g_x_0_z_yyyzzzz[k];

                g_x_0_zz_yyzzzz[k] = -g_x_0_z_yyzzzz[k] * ab_z + g_x_0_z_yyzzzzz[k];

                g_x_0_zz_yzzzzz[k] = -g_x_0_z_yzzzzz[k] * ab_z + g_x_0_z_yzzzzzz[k];

                g_x_0_zz_zzzzzz[k] = -g_x_0_z_zzzzzz[k] * ab_z + g_x_0_z_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_y_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 168 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 169 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 170 * ccomps * dcomps);

            auto g_y_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 171 * ccomps * dcomps);

            auto g_y_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 172 * ccomps * dcomps);

            auto g_y_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 173 * ccomps * dcomps);

            auto g_y_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 174 * ccomps * dcomps);

            auto g_y_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 175 * ccomps * dcomps);

            auto g_y_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 176 * ccomps * dcomps);

            auto g_y_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 177 * ccomps * dcomps);

            auto g_y_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 178 * ccomps * dcomps);

            auto g_y_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 179 * ccomps * dcomps);

            auto g_y_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 180 * ccomps * dcomps);

            auto g_y_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 181 * ccomps * dcomps);

            auto g_y_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 182 * ccomps * dcomps);

            auto g_y_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 183 * ccomps * dcomps);

            auto g_y_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 184 * ccomps * dcomps);

            auto g_y_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 185 * ccomps * dcomps);

            auto g_y_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 186 * ccomps * dcomps);

            auto g_y_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 187 * ccomps * dcomps);

            auto g_y_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 188 * ccomps * dcomps);

            auto g_y_0_xx_yyyyyy = cbuffer.data(di_geom_10_off + 189 * ccomps * dcomps);

            auto g_y_0_xx_yyyyyz = cbuffer.data(di_geom_10_off + 190 * ccomps * dcomps);

            auto g_y_0_xx_yyyyzz = cbuffer.data(di_geom_10_off + 191 * ccomps * dcomps);

            auto g_y_0_xx_yyyzzz = cbuffer.data(di_geom_10_off + 192 * ccomps * dcomps);

            auto g_y_0_xx_yyzzzz = cbuffer.data(di_geom_10_off + 193 * ccomps * dcomps);

            auto g_y_0_xx_yzzzzz = cbuffer.data(di_geom_10_off + 194 * ccomps * dcomps);

            auto g_y_0_xx_zzzzzz = cbuffer.data(di_geom_10_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_xxxxxx, g_y_0_x_xxxxxxx, g_y_0_x_xxxxxxy, g_y_0_x_xxxxxxz, g_y_0_x_xxxxxy, g_y_0_x_xxxxxyy, g_y_0_x_xxxxxyz, g_y_0_x_xxxxxz, g_y_0_x_xxxxxzz, g_y_0_x_xxxxyy, g_y_0_x_xxxxyyy, g_y_0_x_xxxxyyz, g_y_0_x_xxxxyz, g_y_0_x_xxxxyzz, g_y_0_x_xxxxzz, g_y_0_x_xxxxzzz, g_y_0_x_xxxyyy, g_y_0_x_xxxyyyy, g_y_0_x_xxxyyyz, g_y_0_x_xxxyyz, g_y_0_x_xxxyyzz, g_y_0_x_xxxyzz, g_y_0_x_xxxyzzz, g_y_0_x_xxxzzz, g_y_0_x_xxxzzzz, g_y_0_x_xxyyyy, g_y_0_x_xxyyyyy, g_y_0_x_xxyyyyz, g_y_0_x_xxyyyz, g_y_0_x_xxyyyzz, g_y_0_x_xxyyzz, g_y_0_x_xxyyzzz, g_y_0_x_xxyzzz, g_y_0_x_xxyzzzz, g_y_0_x_xxzzzz, g_y_0_x_xxzzzzz, g_y_0_x_xyyyyy, g_y_0_x_xyyyyyy, g_y_0_x_xyyyyyz, g_y_0_x_xyyyyz, g_y_0_x_xyyyyzz, g_y_0_x_xyyyzz, g_y_0_x_xyyyzzz, g_y_0_x_xyyzzz, g_y_0_x_xyyzzzz, g_y_0_x_xyzzzz, g_y_0_x_xyzzzzz, g_y_0_x_xzzzzz, g_y_0_x_xzzzzzz, g_y_0_x_yyyyyy, g_y_0_x_yyyyyz, g_y_0_x_yyyyzz, g_y_0_x_yyyzzz, g_y_0_x_yyzzzz, g_y_0_x_yzzzzz, g_y_0_x_zzzzzz, g_y_0_xx_xxxxxx, g_y_0_xx_xxxxxy, g_y_0_xx_xxxxxz, g_y_0_xx_xxxxyy, g_y_0_xx_xxxxyz, g_y_0_xx_xxxxzz, g_y_0_xx_xxxyyy, g_y_0_xx_xxxyyz, g_y_0_xx_xxxyzz, g_y_0_xx_xxxzzz, g_y_0_xx_xxyyyy, g_y_0_xx_xxyyyz, g_y_0_xx_xxyyzz, g_y_0_xx_xxyzzz, g_y_0_xx_xxzzzz, g_y_0_xx_xyyyyy, g_y_0_xx_xyyyyz, g_y_0_xx_xyyyzz, g_y_0_xx_xyyzzz, g_y_0_xx_xyzzzz, g_y_0_xx_xzzzzz, g_y_0_xx_yyyyyy, g_y_0_xx_yyyyyz, g_y_0_xx_yyyyzz, g_y_0_xx_yyyzzz, g_y_0_xx_yyzzzz, g_y_0_xx_yzzzzz, g_y_0_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xx_xxxxxx[k] = -g_y_0_x_xxxxxx[k] * ab_x + g_y_0_x_xxxxxxx[k];

                g_y_0_xx_xxxxxy[k] = -g_y_0_x_xxxxxy[k] * ab_x + g_y_0_x_xxxxxxy[k];

                g_y_0_xx_xxxxxz[k] = -g_y_0_x_xxxxxz[k] * ab_x + g_y_0_x_xxxxxxz[k];

                g_y_0_xx_xxxxyy[k] = -g_y_0_x_xxxxyy[k] * ab_x + g_y_0_x_xxxxxyy[k];

                g_y_0_xx_xxxxyz[k] = -g_y_0_x_xxxxyz[k] * ab_x + g_y_0_x_xxxxxyz[k];

                g_y_0_xx_xxxxzz[k] = -g_y_0_x_xxxxzz[k] * ab_x + g_y_0_x_xxxxxzz[k];

                g_y_0_xx_xxxyyy[k] = -g_y_0_x_xxxyyy[k] * ab_x + g_y_0_x_xxxxyyy[k];

                g_y_0_xx_xxxyyz[k] = -g_y_0_x_xxxyyz[k] * ab_x + g_y_0_x_xxxxyyz[k];

                g_y_0_xx_xxxyzz[k] = -g_y_0_x_xxxyzz[k] * ab_x + g_y_0_x_xxxxyzz[k];

                g_y_0_xx_xxxzzz[k] = -g_y_0_x_xxxzzz[k] * ab_x + g_y_0_x_xxxxzzz[k];

                g_y_0_xx_xxyyyy[k] = -g_y_0_x_xxyyyy[k] * ab_x + g_y_0_x_xxxyyyy[k];

                g_y_0_xx_xxyyyz[k] = -g_y_0_x_xxyyyz[k] * ab_x + g_y_0_x_xxxyyyz[k];

                g_y_0_xx_xxyyzz[k] = -g_y_0_x_xxyyzz[k] * ab_x + g_y_0_x_xxxyyzz[k];

                g_y_0_xx_xxyzzz[k] = -g_y_0_x_xxyzzz[k] * ab_x + g_y_0_x_xxxyzzz[k];

                g_y_0_xx_xxzzzz[k] = -g_y_0_x_xxzzzz[k] * ab_x + g_y_0_x_xxxzzzz[k];

                g_y_0_xx_xyyyyy[k] = -g_y_0_x_xyyyyy[k] * ab_x + g_y_0_x_xxyyyyy[k];

                g_y_0_xx_xyyyyz[k] = -g_y_0_x_xyyyyz[k] * ab_x + g_y_0_x_xxyyyyz[k];

                g_y_0_xx_xyyyzz[k] = -g_y_0_x_xyyyzz[k] * ab_x + g_y_0_x_xxyyyzz[k];

                g_y_0_xx_xyyzzz[k] = -g_y_0_x_xyyzzz[k] * ab_x + g_y_0_x_xxyyzzz[k];

                g_y_0_xx_xyzzzz[k] = -g_y_0_x_xyzzzz[k] * ab_x + g_y_0_x_xxyzzzz[k];

                g_y_0_xx_xzzzzz[k] = -g_y_0_x_xzzzzz[k] * ab_x + g_y_0_x_xxzzzzz[k];

                g_y_0_xx_yyyyyy[k] = -g_y_0_x_yyyyyy[k] * ab_x + g_y_0_x_xyyyyyy[k];

                g_y_0_xx_yyyyyz[k] = -g_y_0_x_yyyyyz[k] * ab_x + g_y_0_x_xyyyyyz[k];

                g_y_0_xx_yyyyzz[k] = -g_y_0_x_yyyyzz[k] * ab_x + g_y_0_x_xyyyyzz[k];

                g_y_0_xx_yyyzzz[k] = -g_y_0_x_yyyzzz[k] * ab_x + g_y_0_x_xyyyzzz[k];

                g_y_0_xx_yyzzzz[k] = -g_y_0_x_yyzzzz[k] * ab_x + g_y_0_x_xyyzzzz[k];

                g_y_0_xx_yzzzzz[k] = -g_y_0_x_yzzzzz[k] * ab_x + g_y_0_x_xyzzzzz[k];

                g_y_0_xx_zzzzzz[k] = -g_y_0_x_zzzzzz[k] * ab_x + g_y_0_x_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_y_0_xy_xxxxxx = cbuffer.data(di_geom_10_off + 196 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 197 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxz = cbuffer.data(di_geom_10_off + 198 * ccomps * dcomps);

            auto g_y_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 199 * ccomps * dcomps);

            auto g_y_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 200 * ccomps * dcomps);

            auto g_y_0_xy_xxxxzz = cbuffer.data(di_geom_10_off + 201 * ccomps * dcomps);

            auto g_y_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 202 * ccomps * dcomps);

            auto g_y_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 203 * ccomps * dcomps);

            auto g_y_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 204 * ccomps * dcomps);

            auto g_y_0_xy_xxxzzz = cbuffer.data(di_geom_10_off + 205 * ccomps * dcomps);

            auto g_y_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 206 * ccomps * dcomps);

            auto g_y_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 207 * ccomps * dcomps);

            auto g_y_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 208 * ccomps * dcomps);

            auto g_y_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 209 * ccomps * dcomps);

            auto g_y_0_xy_xxzzzz = cbuffer.data(di_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_xy_xzzzzz = cbuffer.data(di_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_xy_yyyyyy = cbuffer.data(di_geom_10_off + 217 * ccomps * dcomps);

            auto g_y_0_xy_yyyyyz = cbuffer.data(di_geom_10_off + 218 * ccomps * dcomps);

            auto g_y_0_xy_yyyyzz = cbuffer.data(di_geom_10_off + 219 * ccomps * dcomps);

            auto g_y_0_xy_yyyzzz = cbuffer.data(di_geom_10_off + 220 * ccomps * dcomps);

            auto g_y_0_xy_yyzzzz = cbuffer.data(di_geom_10_off + 221 * ccomps * dcomps);

            auto g_y_0_xy_yzzzzz = cbuffer.data(di_geom_10_off + 222 * ccomps * dcomps);

            auto g_y_0_xy_zzzzzz = cbuffer.data(di_geom_10_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xy_xxxxxx, g_y_0_xy_xxxxxy, g_y_0_xy_xxxxxz, g_y_0_xy_xxxxyy, g_y_0_xy_xxxxyz, g_y_0_xy_xxxxzz, g_y_0_xy_xxxyyy, g_y_0_xy_xxxyyz, g_y_0_xy_xxxyzz, g_y_0_xy_xxxzzz, g_y_0_xy_xxyyyy, g_y_0_xy_xxyyyz, g_y_0_xy_xxyyzz, g_y_0_xy_xxyzzz, g_y_0_xy_xxzzzz, g_y_0_xy_xyyyyy, g_y_0_xy_xyyyyz, g_y_0_xy_xyyyzz, g_y_0_xy_xyyzzz, g_y_0_xy_xyzzzz, g_y_0_xy_xzzzzz, g_y_0_xy_yyyyyy, g_y_0_xy_yyyyyz, g_y_0_xy_yyyyzz, g_y_0_xy_yyyzzz, g_y_0_xy_yyzzzz, g_y_0_xy_yzzzzz, g_y_0_xy_zzzzzz, g_y_0_y_xxxxxx, g_y_0_y_xxxxxxx, g_y_0_y_xxxxxxy, g_y_0_y_xxxxxxz, g_y_0_y_xxxxxy, g_y_0_y_xxxxxyy, g_y_0_y_xxxxxyz, g_y_0_y_xxxxxz, g_y_0_y_xxxxxzz, g_y_0_y_xxxxyy, g_y_0_y_xxxxyyy, g_y_0_y_xxxxyyz, g_y_0_y_xxxxyz, g_y_0_y_xxxxyzz, g_y_0_y_xxxxzz, g_y_0_y_xxxxzzz, g_y_0_y_xxxyyy, g_y_0_y_xxxyyyy, g_y_0_y_xxxyyyz, g_y_0_y_xxxyyz, g_y_0_y_xxxyyzz, g_y_0_y_xxxyzz, g_y_0_y_xxxyzzz, g_y_0_y_xxxzzz, g_y_0_y_xxxzzzz, g_y_0_y_xxyyyy, g_y_0_y_xxyyyyy, g_y_0_y_xxyyyyz, g_y_0_y_xxyyyz, g_y_0_y_xxyyyzz, g_y_0_y_xxyyzz, g_y_0_y_xxyyzzz, g_y_0_y_xxyzzz, g_y_0_y_xxyzzzz, g_y_0_y_xxzzzz, g_y_0_y_xxzzzzz, g_y_0_y_xyyyyy, g_y_0_y_xyyyyyy, g_y_0_y_xyyyyyz, g_y_0_y_xyyyyz, g_y_0_y_xyyyyzz, g_y_0_y_xyyyzz, g_y_0_y_xyyyzzz, g_y_0_y_xyyzzz, g_y_0_y_xyyzzzz, g_y_0_y_xyzzzz, g_y_0_y_xyzzzzz, g_y_0_y_xzzzzz, g_y_0_y_xzzzzzz, g_y_0_y_yyyyyy, g_y_0_y_yyyyyz, g_y_0_y_yyyyzz, g_y_0_y_yyyzzz, g_y_0_y_yyzzzz, g_y_0_y_yzzzzz, g_y_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xy_xxxxxx[k] = -g_y_0_y_xxxxxx[k] * ab_x + g_y_0_y_xxxxxxx[k];

                g_y_0_xy_xxxxxy[k] = -g_y_0_y_xxxxxy[k] * ab_x + g_y_0_y_xxxxxxy[k];

                g_y_0_xy_xxxxxz[k] = -g_y_0_y_xxxxxz[k] * ab_x + g_y_0_y_xxxxxxz[k];

                g_y_0_xy_xxxxyy[k] = -g_y_0_y_xxxxyy[k] * ab_x + g_y_0_y_xxxxxyy[k];

                g_y_0_xy_xxxxyz[k] = -g_y_0_y_xxxxyz[k] * ab_x + g_y_0_y_xxxxxyz[k];

                g_y_0_xy_xxxxzz[k] = -g_y_0_y_xxxxzz[k] * ab_x + g_y_0_y_xxxxxzz[k];

                g_y_0_xy_xxxyyy[k] = -g_y_0_y_xxxyyy[k] * ab_x + g_y_0_y_xxxxyyy[k];

                g_y_0_xy_xxxyyz[k] = -g_y_0_y_xxxyyz[k] * ab_x + g_y_0_y_xxxxyyz[k];

                g_y_0_xy_xxxyzz[k] = -g_y_0_y_xxxyzz[k] * ab_x + g_y_0_y_xxxxyzz[k];

                g_y_0_xy_xxxzzz[k] = -g_y_0_y_xxxzzz[k] * ab_x + g_y_0_y_xxxxzzz[k];

                g_y_0_xy_xxyyyy[k] = -g_y_0_y_xxyyyy[k] * ab_x + g_y_0_y_xxxyyyy[k];

                g_y_0_xy_xxyyyz[k] = -g_y_0_y_xxyyyz[k] * ab_x + g_y_0_y_xxxyyyz[k];

                g_y_0_xy_xxyyzz[k] = -g_y_0_y_xxyyzz[k] * ab_x + g_y_0_y_xxxyyzz[k];

                g_y_0_xy_xxyzzz[k] = -g_y_0_y_xxyzzz[k] * ab_x + g_y_0_y_xxxyzzz[k];

                g_y_0_xy_xxzzzz[k] = -g_y_0_y_xxzzzz[k] * ab_x + g_y_0_y_xxxzzzz[k];

                g_y_0_xy_xyyyyy[k] = -g_y_0_y_xyyyyy[k] * ab_x + g_y_0_y_xxyyyyy[k];

                g_y_0_xy_xyyyyz[k] = -g_y_0_y_xyyyyz[k] * ab_x + g_y_0_y_xxyyyyz[k];

                g_y_0_xy_xyyyzz[k] = -g_y_0_y_xyyyzz[k] * ab_x + g_y_0_y_xxyyyzz[k];

                g_y_0_xy_xyyzzz[k] = -g_y_0_y_xyyzzz[k] * ab_x + g_y_0_y_xxyyzzz[k];

                g_y_0_xy_xyzzzz[k] = -g_y_0_y_xyzzzz[k] * ab_x + g_y_0_y_xxyzzzz[k];

                g_y_0_xy_xzzzzz[k] = -g_y_0_y_xzzzzz[k] * ab_x + g_y_0_y_xxzzzzz[k];

                g_y_0_xy_yyyyyy[k] = -g_y_0_y_yyyyyy[k] * ab_x + g_y_0_y_xyyyyyy[k];

                g_y_0_xy_yyyyyz[k] = -g_y_0_y_yyyyyz[k] * ab_x + g_y_0_y_xyyyyyz[k];

                g_y_0_xy_yyyyzz[k] = -g_y_0_y_yyyyzz[k] * ab_x + g_y_0_y_xyyyyzz[k];

                g_y_0_xy_yyyzzz[k] = -g_y_0_y_yyyzzz[k] * ab_x + g_y_0_y_xyyyzzz[k];

                g_y_0_xy_yyzzzz[k] = -g_y_0_y_yyzzzz[k] * ab_x + g_y_0_y_xyyzzzz[k];

                g_y_0_xy_yzzzzz[k] = -g_y_0_y_yzzzzz[k] * ab_x + g_y_0_y_xyzzzzz[k];

                g_y_0_xy_zzzzzz[k] = -g_y_0_y_zzzzzz[k] * ab_x + g_y_0_y_xzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_y_0_xz_xxxxxx = cbuffer.data(di_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_xz_yyyyyy = cbuffer.data(di_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_xz_yyyyyz = cbuffer.data(di_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_xz_yyyyzz = cbuffer.data(di_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_xz_yyyzzz = cbuffer.data(di_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_xz_yyzzzz = cbuffer.data(di_geom_10_off + 249 * ccomps * dcomps);

            auto g_y_0_xz_yzzzzz = cbuffer.data(di_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_xz_zzzzzz = cbuffer.data(di_geom_10_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xz_xxxxxx, g_y_0_xz_xxxxxy, g_y_0_xz_xxxxxz, g_y_0_xz_xxxxyy, g_y_0_xz_xxxxyz, g_y_0_xz_xxxxzz, g_y_0_xz_xxxyyy, g_y_0_xz_xxxyyz, g_y_0_xz_xxxyzz, g_y_0_xz_xxxzzz, g_y_0_xz_xxyyyy, g_y_0_xz_xxyyyz, g_y_0_xz_xxyyzz, g_y_0_xz_xxyzzz, g_y_0_xz_xxzzzz, g_y_0_xz_xyyyyy, g_y_0_xz_xyyyyz, g_y_0_xz_xyyyzz, g_y_0_xz_xyyzzz, g_y_0_xz_xyzzzz, g_y_0_xz_xzzzzz, g_y_0_xz_yyyyyy, g_y_0_xz_yyyyyz, g_y_0_xz_yyyyzz, g_y_0_xz_yyyzzz, g_y_0_xz_yyzzzz, g_y_0_xz_yzzzzz, g_y_0_xz_zzzzzz, g_y_0_z_xxxxxx, g_y_0_z_xxxxxxx, g_y_0_z_xxxxxxy, g_y_0_z_xxxxxxz, g_y_0_z_xxxxxy, g_y_0_z_xxxxxyy, g_y_0_z_xxxxxyz, g_y_0_z_xxxxxz, g_y_0_z_xxxxxzz, g_y_0_z_xxxxyy, g_y_0_z_xxxxyyy, g_y_0_z_xxxxyyz, g_y_0_z_xxxxyz, g_y_0_z_xxxxyzz, g_y_0_z_xxxxzz, g_y_0_z_xxxxzzz, g_y_0_z_xxxyyy, g_y_0_z_xxxyyyy, g_y_0_z_xxxyyyz, g_y_0_z_xxxyyz, g_y_0_z_xxxyyzz, g_y_0_z_xxxyzz, g_y_0_z_xxxyzzz, g_y_0_z_xxxzzz, g_y_0_z_xxxzzzz, g_y_0_z_xxyyyy, g_y_0_z_xxyyyyy, g_y_0_z_xxyyyyz, g_y_0_z_xxyyyz, g_y_0_z_xxyyyzz, g_y_0_z_xxyyzz, g_y_0_z_xxyyzzz, g_y_0_z_xxyzzz, g_y_0_z_xxyzzzz, g_y_0_z_xxzzzz, g_y_0_z_xxzzzzz, g_y_0_z_xyyyyy, g_y_0_z_xyyyyyy, g_y_0_z_xyyyyyz, g_y_0_z_xyyyyz, g_y_0_z_xyyyyzz, g_y_0_z_xyyyzz, g_y_0_z_xyyyzzz, g_y_0_z_xyyzzz, g_y_0_z_xyyzzzz, g_y_0_z_xyzzzz, g_y_0_z_xyzzzzz, g_y_0_z_xzzzzz, g_y_0_z_xzzzzzz, g_y_0_z_yyyyyy, g_y_0_z_yyyyyz, g_y_0_z_yyyyzz, g_y_0_z_yyyzzz, g_y_0_z_yyzzzz, g_y_0_z_yzzzzz, g_y_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xz_xxxxxx[k] = -g_y_0_z_xxxxxx[k] * ab_x + g_y_0_z_xxxxxxx[k];

                g_y_0_xz_xxxxxy[k] = -g_y_0_z_xxxxxy[k] * ab_x + g_y_0_z_xxxxxxy[k];

                g_y_0_xz_xxxxxz[k] = -g_y_0_z_xxxxxz[k] * ab_x + g_y_0_z_xxxxxxz[k];

                g_y_0_xz_xxxxyy[k] = -g_y_0_z_xxxxyy[k] * ab_x + g_y_0_z_xxxxxyy[k];

                g_y_0_xz_xxxxyz[k] = -g_y_0_z_xxxxyz[k] * ab_x + g_y_0_z_xxxxxyz[k];

                g_y_0_xz_xxxxzz[k] = -g_y_0_z_xxxxzz[k] * ab_x + g_y_0_z_xxxxxzz[k];

                g_y_0_xz_xxxyyy[k] = -g_y_0_z_xxxyyy[k] * ab_x + g_y_0_z_xxxxyyy[k];

                g_y_0_xz_xxxyyz[k] = -g_y_0_z_xxxyyz[k] * ab_x + g_y_0_z_xxxxyyz[k];

                g_y_0_xz_xxxyzz[k] = -g_y_0_z_xxxyzz[k] * ab_x + g_y_0_z_xxxxyzz[k];

                g_y_0_xz_xxxzzz[k] = -g_y_0_z_xxxzzz[k] * ab_x + g_y_0_z_xxxxzzz[k];

                g_y_0_xz_xxyyyy[k] = -g_y_0_z_xxyyyy[k] * ab_x + g_y_0_z_xxxyyyy[k];

                g_y_0_xz_xxyyyz[k] = -g_y_0_z_xxyyyz[k] * ab_x + g_y_0_z_xxxyyyz[k];

                g_y_0_xz_xxyyzz[k] = -g_y_0_z_xxyyzz[k] * ab_x + g_y_0_z_xxxyyzz[k];

                g_y_0_xz_xxyzzz[k] = -g_y_0_z_xxyzzz[k] * ab_x + g_y_0_z_xxxyzzz[k];

                g_y_0_xz_xxzzzz[k] = -g_y_0_z_xxzzzz[k] * ab_x + g_y_0_z_xxxzzzz[k];

                g_y_0_xz_xyyyyy[k] = -g_y_0_z_xyyyyy[k] * ab_x + g_y_0_z_xxyyyyy[k];

                g_y_0_xz_xyyyyz[k] = -g_y_0_z_xyyyyz[k] * ab_x + g_y_0_z_xxyyyyz[k];

                g_y_0_xz_xyyyzz[k] = -g_y_0_z_xyyyzz[k] * ab_x + g_y_0_z_xxyyyzz[k];

                g_y_0_xz_xyyzzz[k] = -g_y_0_z_xyyzzz[k] * ab_x + g_y_0_z_xxyyzzz[k];

                g_y_0_xz_xyzzzz[k] = -g_y_0_z_xyzzzz[k] * ab_x + g_y_0_z_xxyzzzz[k];

                g_y_0_xz_xzzzzz[k] = -g_y_0_z_xzzzzz[k] * ab_x + g_y_0_z_xxzzzzz[k];

                g_y_0_xz_yyyyyy[k] = -g_y_0_z_yyyyyy[k] * ab_x + g_y_0_z_xyyyyyy[k];

                g_y_0_xz_yyyyyz[k] = -g_y_0_z_yyyyyz[k] * ab_x + g_y_0_z_xyyyyyz[k];

                g_y_0_xz_yyyyzz[k] = -g_y_0_z_yyyyzz[k] * ab_x + g_y_0_z_xyyyyzz[k];

                g_y_0_xz_yyyzzz[k] = -g_y_0_z_yyyzzz[k] * ab_x + g_y_0_z_xyyyzzz[k];

                g_y_0_xz_yyzzzz[k] = -g_y_0_z_yyzzzz[k] * ab_x + g_y_0_z_xyyzzzz[k];

                g_y_0_xz_yzzzzz[k] = -g_y_0_z_yzzzzz[k] * ab_x + g_y_0_z_xyzzzzz[k];

                g_y_0_xz_zzzzzz[k] = -g_y_0_z_zzzzzz[k] * ab_x + g_y_0_z_xzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_y_0_yy_xxxxxx = cbuffer.data(di_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxz = cbuffer.data(di_geom_10_off + 254 * ccomps * dcomps);

            auto g_y_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_yy_xxxxzz = cbuffer.data(di_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 259 * ccomps * dcomps);

            auto g_y_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_yy_xxxzzz = cbuffer.data(di_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_yy_xxzzzz = cbuffer.data(di_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 269 * ccomps * dcomps);

            auto g_y_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_yy_xzzzzz = cbuffer.data(di_geom_10_off + 272 * ccomps * dcomps);

            auto g_y_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_yy_zzzzzz = cbuffer.data(di_geom_10_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxxxxx, g_y_0_y_xxxxxxy, g_y_0_y_xxxxxy, g_y_0_y_xxxxxyy, g_y_0_y_xxxxxyz, g_y_0_y_xxxxxz, g_y_0_y_xxxxyy, g_y_0_y_xxxxyyy, g_y_0_y_xxxxyyz, g_y_0_y_xxxxyz, g_y_0_y_xxxxyzz, g_y_0_y_xxxxzz, g_y_0_y_xxxyyy, g_y_0_y_xxxyyyy, g_y_0_y_xxxyyyz, g_y_0_y_xxxyyz, g_y_0_y_xxxyyzz, g_y_0_y_xxxyzz, g_y_0_y_xxxyzzz, g_y_0_y_xxxzzz, g_y_0_y_xxyyyy, g_y_0_y_xxyyyyy, g_y_0_y_xxyyyyz, g_y_0_y_xxyyyz, g_y_0_y_xxyyyzz, g_y_0_y_xxyyzz, g_y_0_y_xxyyzzz, g_y_0_y_xxyzzz, g_y_0_y_xxyzzzz, g_y_0_y_xxzzzz, g_y_0_y_xyyyyy, g_y_0_y_xyyyyyy, g_y_0_y_xyyyyyz, g_y_0_y_xyyyyz, g_y_0_y_xyyyyzz, g_y_0_y_xyyyzz, g_y_0_y_xyyyzzz, g_y_0_y_xyyzzz, g_y_0_y_xyyzzzz, g_y_0_y_xyzzzz, g_y_0_y_xyzzzzz, g_y_0_y_xzzzzz, g_y_0_y_yyyyyy, g_y_0_y_yyyyyyy, g_y_0_y_yyyyyyz, g_y_0_y_yyyyyz, g_y_0_y_yyyyyzz, g_y_0_y_yyyyzz, g_y_0_y_yyyyzzz, g_y_0_y_yyyzzz, g_y_0_y_yyyzzzz, g_y_0_y_yyzzzz, g_y_0_y_yyzzzzz, g_y_0_y_yzzzzz, g_y_0_y_yzzzzzz, g_y_0_y_zzzzzz, g_y_0_yy_xxxxxx, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_yyyyyy, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_zzzzzz, g_y_xxxxxx, g_y_xxxxxy, g_y_xxxxxz, g_y_xxxxyy, g_y_xxxxyz, g_y_xxxxzz, g_y_xxxyyy, g_y_xxxyyz, g_y_xxxyzz, g_y_xxxzzz, g_y_xxyyyy, g_y_xxyyyz, g_y_xxyyzz, g_y_xxyzzz, g_y_xxzzzz, g_y_xyyyyy, g_y_xyyyyz, g_y_xyyyzz, g_y_xyyzzz, g_y_xyzzzz, g_y_xzzzzz, g_y_yyyyyy, g_y_yyyyyz, g_y_yyyyzz, g_y_yyyzzz, g_y_yyzzzz, g_y_yzzzzz, g_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yy_xxxxxx[k] = -g_y_xxxxxx[k] - g_y_0_y_xxxxxx[k] * ab_y + g_y_0_y_xxxxxxy[k];

                g_y_0_yy_xxxxxy[k] = -g_y_xxxxxy[k] - g_y_0_y_xxxxxy[k] * ab_y + g_y_0_y_xxxxxyy[k];

                g_y_0_yy_xxxxxz[k] = -g_y_xxxxxz[k] - g_y_0_y_xxxxxz[k] * ab_y + g_y_0_y_xxxxxyz[k];

                g_y_0_yy_xxxxyy[k] = -g_y_xxxxyy[k] - g_y_0_y_xxxxyy[k] * ab_y + g_y_0_y_xxxxyyy[k];

                g_y_0_yy_xxxxyz[k] = -g_y_xxxxyz[k] - g_y_0_y_xxxxyz[k] * ab_y + g_y_0_y_xxxxyyz[k];

                g_y_0_yy_xxxxzz[k] = -g_y_xxxxzz[k] - g_y_0_y_xxxxzz[k] * ab_y + g_y_0_y_xxxxyzz[k];

                g_y_0_yy_xxxyyy[k] = -g_y_xxxyyy[k] - g_y_0_y_xxxyyy[k] * ab_y + g_y_0_y_xxxyyyy[k];

                g_y_0_yy_xxxyyz[k] = -g_y_xxxyyz[k] - g_y_0_y_xxxyyz[k] * ab_y + g_y_0_y_xxxyyyz[k];

                g_y_0_yy_xxxyzz[k] = -g_y_xxxyzz[k] - g_y_0_y_xxxyzz[k] * ab_y + g_y_0_y_xxxyyzz[k];

                g_y_0_yy_xxxzzz[k] = -g_y_xxxzzz[k] - g_y_0_y_xxxzzz[k] * ab_y + g_y_0_y_xxxyzzz[k];

                g_y_0_yy_xxyyyy[k] = -g_y_xxyyyy[k] - g_y_0_y_xxyyyy[k] * ab_y + g_y_0_y_xxyyyyy[k];

                g_y_0_yy_xxyyyz[k] = -g_y_xxyyyz[k] - g_y_0_y_xxyyyz[k] * ab_y + g_y_0_y_xxyyyyz[k];

                g_y_0_yy_xxyyzz[k] = -g_y_xxyyzz[k] - g_y_0_y_xxyyzz[k] * ab_y + g_y_0_y_xxyyyzz[k];

                g_y_0_yy_xxyzzz[k] = -g_y_xxyzzz[k] - g_y_0_y_xxyzzz[k] * ab_y + g_y_0_y_xxyyzzz[k];

                g_y_0_yy_xxzzzz[k] = -g_y_xxzzzz[k] - g_y_0_y_xxzzzz[k] * ab_y + g_y_0_y_xxyzzzz[k];

                g_y_0_yy_xyyyyy[k] = -g_y_xyyyyy[k] - g_y_0_y_xyyyyy[k] * ab_y + g_y_0_y_xyyyyyy[k];

                g_y_0_yy_xyyyyz[k] = -g_y_xyyyyz[k] - g_y_0_y_xyyyyz[k] * ab_y + g_y_0_y_xyyyyyz[k];

                g_y_0_yy_xyyyzz[k] = -g_y_xyyyzz[k] - g_y_0_y_xyyyzz[k] * ab_y + g_y_0_y_xyyyyzz[k];

                g_y_0_yy_xyyzzz[k] = -g_y_xyyzzz[k] - g_y_0_y_xyyzzz[k] * ab_y + g_y_0_y_xyyyzzz[k];

                g_y_0_yy_xyzzzz[k] = -g_y_xyzzzz[k] - g_y_0_y_xyzzzz[k] * ab_y + g_y_0_y_xyyzzzz[k];

                g_y_0_yy_xzzzzz[k] = -g_y_xzzzzz[k] - g_y_0_y_xzzzzz[k] * ab_y + g_y_0_y_xyzzzzz[k];

                g_y_0_yy_yyyyyy[k] = -g_y_yyyyyy[k] - g_y_0_y_yyyyyy[k] * ab_y + g_y_0_y_yyyyyyy[k];

                g_y_0_yy_yyyyyz[k] = -g_y_yyyyyz[k] - g_y_0_y_yyyyyz[k] * ab_y + g_y_0_y_yyyyyyz[k];

                g_y_0_yy_yyyyzz[k] = -g_y_yyyyzz[k] - g_y_0_y_yyyyzz[k] * ab_y + g_y_0_y_yyyyyzz[k];

                g_y_0_yy_yyyzzz[k] = -g_y_yyyzzz[k] - g_y_0_y_yyyzzz[k] * ab_y + g_y_0_y_yyyyzzz[k];

                g_y_0_yy_yyzzzz[k] = -g_y_yyzzzz[k] - g_y_0_y_yyzzzz[k] * ab_y + g_y_0_y_yyyzzzz[k];

                g_y_0_yy_yzzzzz[k] = -g_y_yzzzzz[k] - g_y_0_y_yzzzzz[k] * ab_y + g_y_0_y_yyzzzzz[k];

                g_y_0_yy_zzzzzz[k] = -g_y_zzzzzz[k] - g_y_0_y_zzzzzz[k] * ab_y + g_y_0_y_yzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_y_0_yz_xxxxxx = cbuffer.data(di_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxz = cbuffer.data(di_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_yz_xxxxzz = cbuffer.data(di_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_yz_xxxzzz = cbuffer.data(di_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_yz_xxzzzz = cbuffer.data(di_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 299 * ccomps * dcomps);

            auto g_y_0_yz_xzzzzz = cbuffer.data(di_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_yz_yyyyyy = cbuffer.data(di_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_yz_zzzzzz = cbuffer.data(di_geom_10_off + 307 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxxxxx, g_y_0_y_xxxxxxz, g_y_0_y_xxxxxy, g_y_0_y_xxxxxyz, g_y_0_y_xxxxxz, g_y_0_y_xxxxxzz, g_y_0_y_xxxxyy, g_y_0_y_xxxxyyz, g_y_0_y_xxxxyz, g_y_0_y_xxxxyzz, g_y_0_y_xxxxzz, g_y_0_y_xxxxzzz, g_y_0_y_xxxyyy, g_y_0_y_xxxyyyz, g_y_0_y_xxxyyz, g_y_0_y_xxxyyzz, g_y_0_y_xxxyzz, g_y_0_y_xxxyzzz, g_y_0_y_xxxzzz, g_y_0_y_xxxzzzz, g_y_0_y_xxyyyy, g_y_0_y_xxyyyyz, g_y_0_y_xxyyyz, g_y_0_y_xxyyyzz, g_y_0_y_xxyyzz, g_y_0_y_xxyyzzz, g_y_0_y_xxyzzz, g_y_0_y_xxyzzzz, g_y_0_y_xxzzzz, g_y_0_y_xxzzzzz, g_y_0_y_xyyyyy, g_y_0_y_xyyyyyz, g_y_0_y_xyyyyz, g_y_0_y_xyyyyzz, g_y_0_y_xyyyzz, g_y_0_y_xyyyzzz, g_y_0_y_xyyzzz, g_y_0_y_xyyzzzz, g_y_0_y_xyzzzz, g_y_0_y_xyzzzzz, g_y_0_y_xzzzzz, g_y_0_y_xzzzzzz, g_y_0_y_yyyyyy, g_y_0_y_yyyyyyz, g_y_0_y_yyyyyz, g_y_0_y_yyyyyzz, g_y_0_y_yyyyzz, g_y_0_y_yyyyzzz, g_y_0_y_yyyzzz, g_y_0_y_yyyzzzz, g_y_0_y_yyzzzz, g_y_0_y_yyzzzzz, g_y_0_y_yzzzzz, g_y_0_y_yzzzzzz, g_y_0_y_zzzzzz, g_y_0_y_zzzzzzz, g_y_0_yz_xxxxxx, g_y_0_yz_xxxxxy, g_y_0_yz_xxxxxz, g_y_0_yz_xxxxyy, g_y_0_yz_xxxxyz, g_y_0_yz_xxxxzz, g_y_0_yz_xxxyyy, g_y_0_yz_xxxyyz, g_y_0_yz_xxxyzz, g_y_0_yz_xxxzzz, g_y_0_yz_xxyyyy, g_y_0_yz_xxyyyz, g_y_0_yz_xxyyzz, g_y_0_yz_xxyzzz, g_y_0_yz_xxzzzz, g_y_0_yz_xyyyyy, g_y_0_yz_xyyyyz, g_y_0_yz_xyyyzz, g_y_0_yz_xyyzzz, g_y_0_yz_xyzzzz, g_y_0_yz_xzzzzz, g_y_0_yz_yyyyyy, g_y_0_yz_yyyyyz, g_y_0_yz_yyyyzz, g_y_0_yz_yyyzzz, g_y_0_yz_yyzzzz, g_y_0_yz_yzzzzz, g_y_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yz_xxxxxx[k] = -g_y_0_y_xxxxxx[k] * ab_z + g_y_0_y_xxxxxxz[k];

                g_y_0_yz_xxxxxy[k] = -g_y_0_y_xxxxxy[k] * ab_z + g_y_0_y_xxxxxyz[k];

                g_y_0_yz_xxxxxz[k] = -g_y_0_y_xxxxxz[k] * ab_z + g_y_0_y_xxxxxzz[k];

                g_y_0_yz_xxxxyy[k] = -g_y_0_y_xxxxyy[k] * ab_z + g_y_0_y_xxxxyyz[k];

                g_y_0_yz_xxxxyz[k] = -g_y_0_y_xxxxyz[k] * ab_z + g_y_0_y_xxxxyzz[k];

                g_y_0_yz_xxxxzz[k] = -g_y_0_y_xxxxzz[k] * ab_z + g_y_0_y_xxxxzzz[k];

                g_y_0_yz_xxxyyy[k] = -g_y_0_y_xxxyyy[k] * ab_z + g_y_0_y_xxxyyyz[k];

                g_y_0_yz_xxxyyz[k] = -g_y_0_y_xxxyyz[k] * ab_z + g_y_0_y_xxxyyzz[k];

                g_y_0_yz_xxxyzz[k] = -g_y_0_y_xxxyzz[k] * ab_z + g_y_0_y_xxxyzzz[k];

                g_y_0_yz_xxxzzz[k] = -g_y_0_y_xxxzzz[k] * ab_z + g_y_0_y_xxxzzzz[k];

                g_y_0_yz_xxyyyy[k] = -g_y_0_y_xxyyyy[k] * ab_z + g_y_0_y_xxyyyyz[k];

                g_y_0_yz_xxyyyz[k] = -g_y_0_y_xxyyyz[k] * ab_z + g_y_0_y_xxyyyzz[k];

                g_y_0_yz_xxyyzz[k] = -g_y_0_y_xxyyzz[k] * ab_z + g_y_0_y_xxyyzzz[k];

                g_y_0_yz_xxyzzz[k] = -g_y_0_y_xxyzzz[k] * ab_z + g_y_0_y_xxyzzzz[k];

                g_y_0_yz_xxzzzz[k] = -g_y_0_y_xxzzzz[k] * ab_z + g_y_0_y_xxzzzzz[k];

                g_y_0_yz_xyyyyy[k] = -g_y_0_y_xyyyyy[k] * ab_z + g_y_0_y_xyyyyyz[k];

                g_y_0_yz_xyyyyz[k] = -g_y_0_y_xyyyyz[k] * ab_z + g_y_0_y_xyyyyzz[k];

                g_y_0_yz_xyyyzz[k] = -g_y_0_y_xyyyzz[k] * ab_z + g_y_0_y_xyyyzzz[k];

                g_y_0_yz_xyyzzz[k] = -g_y_0_y_xyyzzz[k] * ab_z + g_y_0_y_xyyzzzz[k];

                g_y_0_yz_xyzzzz[k] = -g_y_0_y_xyzzzz[k] * ab_z + g_y_0_y_xyzzzzz[k];

                g_y_0_yz_xzzzzz[k] = -g_y_0_y_xzzzzz[k] * ab_z + g_y_0_y_xzzzzzz[k];

                g_y_0_yz_yyyyyy[k] = -g_y_0_y_yyyyyy[k] * ab_z + g_y_0_y_yyyyyyz[k];

                g_y_0_yz_yyyyyz[k] = -g_y_0_y_yyyyyz[k] * ab_z + g_y_0_y_yyyyyzz[k];

                g_y_0_yz_yyyyzz[k] = -g_y_0_y_yyyyzz[k] * ab_z + g_y_0_y_yyyyzzz[k];

                g_y_0_yz_yyyzzz[k] = -g_y_0_y_yyyzzz[k] * ab_z + g_y_0_y_yyyzzzz[k];

                g_y_0_yz_yyzzzz[k] = -g_y_0_y_yyzzzz[k] * ab_z + g_y_0_y_yyzzzzz[k];

                g_y_0_yz_yzzzzz[k] = -g_y_0_y_yzzzzz[k] * ab_z + g_y_0_y_yzzzzzz[k];

                g_y_0_yz_zzzzzz[k] = -g_y_0_y_zzzzzz[k] * ab_z + g_y_0_y_zzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_y_0_zz_xxxxxx = cbuffer.data(di_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 309 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_zz_yyyyyy = cbuffer.data(di_geom_10_off + 329 * ccomps * dcomps);

            auto g_y_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxxxxx, g_y_0_z_xxxxxxz, g_y_0_z_xxxxxy, g_y_0_z_xxxxxyz, g_y_0_z_xxxxxz, g_y_0_z_xxxxxzz, g_y_0_z_xxxxyy, g_y_0_z_xxxxyyz, g_y_0_z_xxxxyz, g_y_0_z_xxxxyzz, g_y_0_z_xxxxzz, g_y_0_z_xxxxzzz, g_y_0_z_xxxyyy, g_y_0_z_xxxyyyz, g_y_0_z_xxxyyz, g_y_0_z_xxxyyzz, g_y_0_z_xxxyzz, g_y_0_z_xxxyzzz, g_y_0_z_xxxzzz, g_y_0_z_xxxzzzz, g_y_0_z_xxyyyy, g_y_0_z_xxyyyyz, g_y_0_z_xxyyyz, g_y_0_z_xxyyyzz, g_y_0_z_xxyyzz, g_y_0_z_xxyyzzz, g_y_0_z_xxyzzz, g_y_0_z_xxyzzzz, g_y_0_z_xxzzzz, g_y_0_z_xxzzzzz, g_y_0_z_xyyyyy, g_y_0_z_xyyyyyz, g_y_0_z_xyyyyz, g_y_0_z_xyyyyzz, g_y_0_z_xyyyzz, g_y_0_z_xyyyzzz, g_y_0_z_xyyzzz, g_y_0_z_xyyzzzz, g_y_0_z_xyzzzz, g_y_0_z_xyzzzzz, g_y_0_z_xzzzzz, g_y_0_z_xzzzzzz, g_y_0_z_yyyyyy, g_y_0_z_yyyyyyz, g_y_0_z_yyyyyz, g_y_0_z_yyyyyzz, g_y_0_z_yyyyzz, g_y_0_z_yyyyzzz, g_y_0_z_yyyzzz, g_y_0_z_yyyzzzz, g_y_0_z_yyzzzz, g_y_0_z_yyzzzzz, g_y_0_z_yzzzzz, g_y_0_z_yzzzzzz, g_y_0_z_zzzzzz, g_y_0_z_zzzzzzz, g_y_0_zz_xxxxxx, g_y_0_zz_xxxxxy, g_y_0_zz_xxxxxz, g_y_0_zz_xxxxyy, g_y_0_zz_xxxxyz, g_y_0_zz_xxxxzz, g_y_0_zz_xxxyyy, g_y_0_zz_xxxyyz, g_y_0_zz_xxxyzz, g_y_0_zz_xxxzzz, g_y_0_zz_xxyyyy, g_y_0_zz_xxyyyz, g_y_0_zz_xxyyzz, g_y_0_zz_xxyzzz, g_y_0_zz_xxzzzz, g_y_0_zz_xyyyyy, g_y_0_zz_xyyyyz, g_y_0_zz_xyyyzz, g_y_0_zz_xyyzzz, g_y_0_zz_xyzzzz, g_y_0_zz_xzzzzz, g_y_0_zz_yyyyyy, g_y_0_zz_yyyyyz, g_y_0_zz_yyyyzz, g_y_0_zz_yyyzzz, g_y_0_zz_yyzzzz, g_y_0_zz_yzzzzz, g_y_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zz_xxxxxx[k] = -g_y_0_z_xxxxxx[k] * ab_z + g_y_0_z_xxxxxxz[k];

                g_y_0_zz_xxxxxy[k] = -g_y_0_z_xxxxxy[k] * ab_z + g_y_0_z_xxxxxyz[k];

                g_y_0_zz_xxxxxz[k] = -g_y_0_z_xxxxxz[k] * ab_z + g_y_0_z_xxxxxzz[k];

                g_y_0_zz_xxxxyy[k] = -g_y_0_z_xxxxyy[k] * ab_z + g_y_0_z_xxxxyyz[k];

                g_y_0_zz_xxxxyz[k] = -g_y_0_z_xxxxyz[k] * ab_z + g_y_0_z_xxxxyzz[k];

                g_y_0_zz_xxxxzz[k] = -g_y_0_z_xxxxzz[k] * ab_z + g_y_0_z_xxxxzzz[k];

                g_y_0_zz_xxxyyy[k] = -g_y_0_z_xxxyyy[k] * ab_z + g_y_0_z_xxxyyyz[k];

                g_y_0_zz_xxxyyz[k] = -g_y_0_z_xxxyyz[k] * ab_z + g_y_0_z_xxxyyzz[k];

                g_y_0_zz_xxxyzz[k] = -g_y_0_z_xxxyzz[k] * ab_z + g_y_0_z_xxxyzzz[k];

                g_y_0_zz_xxxzzz[k] = -g_y_0_z_xxxzzz[k] * ab_z + g_y_0_z_xxxzzzz[k];

                g_y_0_zz_xxyyyy[k] = -g_y_0_z_xxyyyy[k] * ab_z + g_y_0_z_xxyyyyz[k];

                g_y_0_zz_xxyyyz[k] = -g_y_0_z_xxyyyz[k] * ab_z + g_y_0_z_xxyyyzz[k];

                g_y_0_zz_xxyyzz[k] = -g_y_0_z_xxyyzz[k] * ab_z + g_y_0_z_xxyyzzz[k];

                g_y_0_zz_xxyzzz[k] = -g_y_0_z_xxyzzz[k] * ab_z + g_y_0_z_xxyzzzz[k];

                g_y_0_zz_xxzzzz[k] = -g_y_0_z_xxzzzz[k] * ab_z + g_y_0_z_xxzzzzz[k];

                g_y_0_zz_xyyyyy[k] = -g_y_0_z_xyyyyy[k] * ab_z + g_y_0_z_xyyyyyz[k];

                g_y_0_zz_xyyyyz[k] = -g_y_0_z_xyyyyz[k] * ab_z + g_y_0_z_xyyyyzz[k];

                g_y_0_zz_xyyyzz[k] = -g_y_0_z_xyyyzz[k] * ab_z + g_y_0_z_xyyyzzz[k];

                g_y_0_zz_xyyzzz[k] = -g_y_0_z_xyyzzz[k] * ab_z + g_y_0_z_xyyzzzz[k];

                g_y_0_zz_xyzzzz[k] = -g_y_0_z_xyzzzz[k] * ab_z + g_y_0_z_xyzzzzz[k];

                g_y_0_zz_xzzzzz[k] = -g_y_0_z_xzzzzz[k] * ab_z + g_y_0_z_xzzzzzz[k];

                g_y_0_zz_yyyyyy[k] = -g_y_0_z_yyyyyy[k] * ab_z + g_y_0_z_yyyyyyz[k];

                g_y_0_zz_yyyyyz[k] = -g_y_0_z_yyyyyz[k] * ab_z + g_y_0_z_yyyyyzz[k];

                g_y_0_zz_yyyyzz[k] = -g_y_0_z_yyyyzz[k] * ab_z + g_y_0_z_yyyyzzz[k];

                g_y_0_zz_yyyzzz[k] = -g_y_0_z_yyyzzz[k] * ab_z + g_y_0_z_yyyzzzz[k];

                g_y_0_zz_yyzzzz[k] = -g_y_0_z_yyzzzz[k] * ab_z + g_y_0_z_yyzzzzz[k];

                g_y_0_zz_yzzzzz[k] = -g_y_0_z_yzzzzz[k] * ab_z + g_y_0_z_yzzzzzz[k];

                g_y_0_zz_zzzzzz[k] = -g_y_0_z_zzzzzz[k] * ab_z + g_y_0_z_zzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_z_0_xx_xxxxxx = cbuffer.data(di_geom_10_off + 336 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxy = cbuffer.data(di_geom_10_off + 337 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxz = cbuffer.data(di_geom_10_off + 338 * ccomps * dcomps);

            auto g_z_0_xx_xxxxyy = cbuffer.data(di_geom_10_off + 339 * ccomps * dcomps);

            auto g_z_0_xx_xxxxyz = cbuffer.data(di_geom_10_off + 340 * ccomps * dcomps);

            auto g_z_0_xx_xxxxzz = cbuffer.data(di_geom_10_off + 341 * ccomps * dcomps);

            auto g_z_0_xx_xxxyyy = cbuffer.data(di_geom_10_off + 342 * ccomps * dcomps);

            auto g_z_0_xx_xxxyyz = cbuffer.data(di_geom_10_off + 343 * ccomps * dcomps);

            auto g_z_0_xx_xxxyzz = cbuffer.data(di_geom_10_off + 344 * ccomps * dcomps);

            auto g_z_0_xx_xxxzzz = cbuffer.data(di_geom_10_off + 345 * ccomps * dcomps);

            auto g_z_0_xx_xxyyyy = cbuffer.data(di_geom_10_off + 346 * ccomps * dcomps);

            auto g_z_0_xx_xxyyyz = cbuffer.data(di_geom_10_off + 347 * ccomps * dcomps);

            auto g_z_0_xx_xxyyzz = cbuffer.data(di_geom_10_off + 348 * ccomps * dcomps);

            auto g_z_0_xx_xxyzzz = cbuffer.data(di_geom_10_off + 349 * ccomps * dcomps);

            auto g_z_0_xx_xxzzzz = cbuffer.data(di_geom_10_off + 350 * ccomps * dcomps);

            auto g_z_0_xx_xyyyyy = cbuffer.data(di_geom_10_off + 351 * ccomps * dcomps);

            auto g_z_0_xx_xyyyyz = cbuffer.data(di_geom_10_off + 352 * ccomps * dcomps);

            auto g_z_0_xx_xyyyzz = cbuffer.data(di_geom_10_off + 353 * ccomps * dcomps);

            auto g_z_0_xx_xyyzzz = cbuffer.data(di_geom_10_off + 354 * ccomps * dcomps);

            auto g_z_0_xx_xyzzzz = cbuffer.data(di_geom_10_off + 355 * ccomps * dcomps);

            auto g_z_0_xx_xzzzzz = cbuffer.data(di_geom_10_off + 356 * ccomps * dcomps);

            auto g_z_0_xx_yyyyyy = cbuffer.data(di_geom_10_off + 357 * ccomps * dcomps);

            auto g_z_0_xx_yyyyyz = cbuffer.data(di_geom_10_off + 358 * ccomps * dcomps);

            auto g_z_0_xx_yyyyzz = cbuffer.data(di_geom_10_off + 359 * ccomps * dcomps);

            auto g_z_0_xx_yyyzzz = cbuffer.data(di_geom_10_off + 360 * ccomps * dcomps);

            auto g_z_0_xx_yyzzzz = cbuffer.data(di_geom_10_off + 361 * ccomps * dcomps);

            auto g_z_0_xx_yzzzzz = cbuffer.data(di_geom_10_off + 362 * ccomps * dcomps);

            auto g_z_0_xx_zzzzzz = cbuffer.data(di_geom_10_off + 363 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_xxxxxx, g_z_0_x_xxxxxxx, g_z_0_x_xxxxxxy, g_z_0_x_xxxxxxz, g_z_0_x_xxxxxy, g_z_0_x_xxxxxyy, g_z_0_x_xxxxxyz, g_z_0_x_xxxxxz, g_z_0_x_xxxxxzz, g_z_0_x_xxxxyy, g_z_0_x_xxxxyyy, g_z_0_x_xxxxyyz, g_z_0_x_xxxxyz, g_z_0_x_xxxxyzz, g_z_0_x_xxxxzz, g_z_0_x_xxxxzzz, g_z_0_x_xxxyyy, g_z_0_x_xxxyyyy, g_z_0_x_xxxyyyz, g_z_0_x_xxxyyz, g_z_0_x_xxxyyzz, g_z_0_x_xxxyzz, g_z_0_x_xxxyzzz, g_z_0_x_xxxzzz, g_z_0_x_xxxzzzz, g_z_0_x_xxyyyy, g_z_0_x_xxyyyyy, g_z_0_x_xxyyyyz, g_z_0_x_xxyyyz, g_z_0_x_xxyyyzz, g_z_0_x_xxyyzz, g_z_0_x_xxyyzzz, g_z_0_x_xxyzzz, g_z_0_x_xxyzzzz, g_z_0_x_xxzzzz, g_z_0_x_xxzzzzz, g_z_0_x_xyyyyy, g_z_0_x_xyyyyyy, g_z_0_x_xyyyyyz, g_z_0_x_xyyyyz, g_z_0_x_xyyyyzz, g_z_0_x_xyyyzz, g_z_0_x_xyyyzzz, g_z_0_x_xyyzzz, g_z_0_x_xyyzzzz, g_z_0_x_xyzzzz, g_z_0_x_xyzzzzz, g_z_0_x_xzzzzz, g_z_0_x_xzzzzzz, g_z_0_x_yyyyyy, g_z_0_x_yyyyyz, g_z_0_x_yyyyzz, g_z_0_x_yyyzzz, g_z_0_x_yyzzzz, g_z_0_x_yzzzzz, g_z_0_x_zzzzzz, g_z_0_xx_xxxxxx, g_z_0_xx_xxxxxy, g_z_0_xx_xxxxxz, g_z_0_xx_xxxxyy, g_z_0_xx_xxxxyz, g_z_0_xx_xxxxzz, g_z_0_xx_xxxyyy, g_z_0_xx_xxxyyz, g_z_0_xx_xxxyzz, g_z_0_xx_xxxzzz, g_z_0_xx_xxyyyy, g_z_0_xx_xxyyyz, g_z_0_xx_xxyyzz, g_z_0_xx_xxyzzz, g_z_0_xx_xxzzzz, g_z_0_xx_xyyyyy, g_z_0_xx_xyyyyz, g_z_0_xx_xyyyzz, g_z_0_xx_xyyzzz, g_z_0_xx_xyzzzz, g_z_0_xx_xzzzzz, g_z_0_xx_yyyyyy, g_z_0_xx_yyyyyz, g_z_0_xx_yyyyzz, g_z_0_xx_yyyzzz, g_z_0_xx_yyzzzz, g_z_0_xx_yzzzzz, g_z_0_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xx_xxxxxx[k] = -g_z_0_x_xxxxxx[k] * ab_x + g_z_0_x_xxxxxxx[k];

                g_z_0_xx_xxxxxy[k] = -g_z_0_x_xxxxxy[k] * ab_x + g_z_0_x_xxxxxxy[k];

                g_z_0_xx_xxxxxz[k] = -g_z_0_x_xxxxxz[k] * ab_x + g_z_0_x_xxxxxxz[k];

                g_z_0_xx_xxxxyy[k] = -g_z_0_x_xxxxyy[k] * ab_x + g_z_0_x_xxxxxyy[k];

                g_z_0_xx_xxxxyz[k] = -g_z_0_x_xxxxyz[k] * ab_x + g_z_0_x_xxxxxyz[k];

                g_z_0_xx_xxxxzz[k] = -g_z_0_x_xxxxzz[k] * ab_x + g_z_0_x_xxxxxzz[k];

                g_z_0_xx_xxxyyy[k] = -g_z_0_x_xxxyyy[k] * ab_x + g_z_0_x_xxxxyyy[k];

                g_z_0_xx_xxxyyz[k] = -g_z_0_x_xxxyyz[k] * ab_x + g_z_0_x_xxxxyyz[k];

                g_z_0_xx_xxxyzz[k] = -g_z_0_x_xxxyzz[k] * ab_x + g_z_0_x_xxxxyzz[k];

                g_z_0_xx_xxxzzz[k] = -g_z_0_x_xxxzzz[k] * ab_x + g_z_0_x_xxxxzzz[k];

                g_z_0_xx_xxyyyy[k] = -g_z_0_x_xxyyyy[k] * ab_x + g_z_0_x_xxxyyyy[k];

                g_z_0_xx_xxyyyz[k] = -g_z_0_x_xxyyyz[k] * ab_x + g_z_0_x_xxxyyyz[k];

                g_z_0_xx_xxyyzz[k] = -g_z_0_x_xxyyzz[k] * ab_x + g_z_0_x_xxxyyzz[k];

                g_z_0_xx_xxyzzz[k] = -g_z_0_x_xxyzzz[k] * ab_x + g_z_0_x_xxxyzzz[k];

                g_z_0_xx_xxzzzz[k] = -g_z_0_x_xxzzzz[k] * ab_x + g_z_0_x_xxxzzzz[k];

                g_z_0_xx_xyyyyy[k] = -g_z_0_x_xyyyyy[k] * ab_x + g_z_0_x_xxyyyyy[k];

                g_z_0_xx_xyyyyz[k] = -g_z_0_x_xyyyyz[k] * ab_x + g_z_0_x_xxyyyyz[k];

                g_z_0_xx_xyyyzz[k] = -g_z_0_x_xyyyzz[k] * ab_x + g_z_0_x_xxyyyzz[k];

                g_z_0_xx_xyyzzz[k] = -g_z_0_x_xyyzzz[k] * ab_x + g_z_0_x_xxyyzzz[k];

                g_z_0_xx_xyzzzz[k] = -g_z_0_x_xyzzzz[k] * ab_x + g_z_0_x_xxyzzzz[k];

                g_z_0_xx_xzzzzz[k] = -g_z_0_x_xzzzzz[k] * ab_x + g_z_0_x_xxzzzzz[k];

                g_z_0_xx_yyyyyy[k] = -g_z_0_x_yyyyyy[k] * ab_x + g_z_0_x_xyyyyyy[k];

                g_z_0_xx_yyyyyz[k] = -g_z_0_x_yyyyyz[k] * ab_x + g_z_0_x_xyyyyyz[k];

                g_z_0_xx_yyyyzz[k] = -g_z_0_x_yyyyzz[k] * ab_x + g_z_0_x_xyyyyzz[k];

                g_z_0_xx_yyyzzz[k] = -g_z_0_x_yyyzzz[k] * ab_x + g_z_0_x_xyyyzzz[k];

                g_z_0_xx_yyzzzz[k] = -g_z_0_x_yyzzzz[k] * ab_x + g_z_0_x_xyyzzzz[k];

                g_z_0_xx_yzzzzz[k] = -g_z_0_x_yzzzzz[k] * ab_x + g_z_0_x_xyzzzzz[k];

                g_z_0_xx_zzzzzz[k] = -g_z_0_x_zzzzzz[k] * ab_x + g_z_0_x_xzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_z_0_xy_xxxxxx = cbuffer.data(di_geom_10_off + 364 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxy = cbuffer.data(di_geom_10_off + 365 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxz = cbuffer.data(di_geom_10_off + 366 * ccomps * dcomps);

            auto g_z_0_xy_xxxxyy = cbuffer.data(di_geom_10_off + 367 * ccomps * dcomps);

            auto g_z_0_xy_xxxxyz = cbuffer.data(di_geom_10_off + 368 * ccomps * dcomps);

            auto g_z_0_xy_xxxxzz = cbuffer.data(di_geom_10_off + 369 * ccomps * dcomps);

            auto g_z_0_xy_xxxyyy = cbuffer.data(di_geom_10_off + 370 * ccomps * dcomps);

            auto g_z_0_xy_xxxyyz = cbuffer.data(di_geom_10_off + 371 * ccomps * dcomps);

            auto g_z_0_xy_xxxyzz = cbuffer.data(di_geom_10_off + 372 * ccomps * dcomps);

            auto g_z_0_xy_xxxzzz = cbuffer.data(di_geom_10_off + 373 * ccomps * dcomps);

            auto g_z_0_xy_xxyyyy = cbuffer.data(di_geom_10_off + 374 * ccomps * dcomps);

            auto g_z_0_xy_xxyyyz = cbuffer.data(di_geom_10_off + 375 * ccomps * dcomps);

            auto g_z_0_xy_xxyyzz = cbuffer.data(di_geom_10_off + 376 * ccomps * dcomps);

            auto g_z_0_xy_xxyzzz = cbuffer.data(di_geom_10_off + 377 * ccomps * dcomps);

            auto g_z_0_xy_xxzzzz = cbuffer.data(di_geom_10_off + 378 * ccomps * dcomps);

            auto g_z_0_xy_xyyyyy = cbuffer.data(di_geom_10_off + 379 * ccomps * dcomps);

            auto g_z_0_xy_xyyyyz = cbuffer.data(di_geom_10_off + 380 * ccomps * dcomps);

            auto g_z_0_xy_xyyyzz = cbuffer.data(di_geom_10_off + 381 * ccomps * dcomps);

            auto g_z_0_xy_xyyzzz = cbuffer.data(di_geom_10_off + 382 * ccomps * dcomps);

            auto g_z_0_xy_xyzzzz = cbuffer.data(di_geom_10_off + 383 * ccomps * dcomps);

            auto g_z_0_xy_xzzzzz = cbuffer.data(di_geom_10_off + 384 * ccomps * dcomps);

            auto g_z_0_xy_yyyyyy = cbuffer.data(di_geom_10_off + 385 * ccomps * dcomps);

            auto g_z_0_xy_yyyyyz = cbuffer.data(di_geom_10_off + 386 * ccomps * dcomps);

            auto g_z_0_xy_yyyyzz = cbuffer.data(di_geom_10_off + 387 * ccomps * dcomps);

            auto g_z_0_xy_yyyzzz = cbuffer.data(di_geom_10_off + 388 * ccomps * dcomps);

            auto g_z_0_xy_yyzzzz = cbuffer.data(di_geom_10_off + 389 * ccomps * dcomps);

            auto g_z_0_xy_yzzzzz = cbuffer.data(di_geom_10_off + 390 * ccomps * dcomps);

            auto g_z_0_xy_zzzzzz = cbuffer.data(di_geom_10_off + 391 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xy_xxxxxx, g_z_0_xy_xxxxxy, g_z_0_xy_xxxxxz, g_z_0_xy_xxxxyy, g_z_0_xy_xxxxyz, g_z_0_xy_xxxxzz, g_z_0_xy_xxxyyy, g_z_0_xy_xxxyyz, g_z_0_xy_xxxyzz, g_z_0_xy_xxxzzz, g_z_0_xy_xxyyyy, g_z_0_xy_xxyyyz, g_z_0_xy_xxyyzz, g_z_0_xy_xxyzzz, g_z_0_xy_xxzzzz, g_z_0_xy_xyyyyy, g_z_0_xy_xyyyyz, g_z_0_xy_xyyyzz, g_z_0_xy_xyyzzz, g_z_0_xy_xyzzzz, g_z_0_xy_xzzzzz, g_z_0_xy_yyyyyy, g_z_0_xy_yyyyyz, g_z_0_xy_yyyyzz, g_z_0_xy_yyyzzz, g_z_0_xy_yyzzzz, g_z_0_xy_yzzzzz, g_z_0_xy_zzzzzz, g_z_0_y_xxxxxx, g_z_0_y_xxxxxxx, g_z_0_y_xxxxxxy, g_z_0_y_xxxxxxz, g_z_0_y_xxxxxy, g_z_0_y_xxxxxyy, g_z_0_y_xxxxxyz, g_z_0_y_xxxxxz, g_z_0_y_xxxxxzz, g_z_0_y_xxxxyy, g_z_0_y_xxxxyyy, g_z_0_y_xxxxyyz, g_z_0_y_xxxxyz, g_z_0_y_xxxxyzz, g_z_0_y_xxxxzz, g_z_0_y_xxxxzzz, g_z_0_y_xxxyyy, g_z_0_y_xxxyyyy, g_z_0_y_xxxyyyz, g_z_0_y_xxxyyz, g_z_0_y_xxxyyzz, g_z_0_y_xxxyzz, g_z_0_y_xxxyzzz, g_z_0_y_xxxzzz, g_z_0_y_xxxzzzz, g_z_0_y_xxyyyy, g_z_0_y_xxyyyyy, g_z_0_y_xxyyyyz, g_z_0_y_xxyyyz, g_z_0_y_xxyyyzz, g_z_0_y_xxyyzz, g_z_0_y_xxyyzzz, g_z_0_y_xxyzzz, g_z_0_y_xxyzzzz, g_z_0_y_xxzzzz, g_z_0_y_xxzzzzz, g_z_0_y_xyyyyy, g_z_0_y_xyyyyyy, g_z_0_y_xyyyyyz, g_z_0_y_xyyyyz, g_z_0_y_xyyyyzz, g_z_0_y_xyyyzz, g_z_0_y_xyyyzzz, g_z_0_y_xyyzzz, g_z_0_y_xyyzzzz, g_z_0_y_xyzzzz, g_z_0_y_xyzzzzz, g_z_0_y_xzzzzz, g_z_0_y_xzzzzzz, g_z_0_y_yyyyyy, g_z_0_y_yyyyyz, g_z_0_y_yyyyzz, g_z_0_y_yyyzzz, g_z_0_y_yyzzzz, g_z_0_y_yzzzzz, g_z_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xy_xxxxxx[k] = -g_z_0_y_xxxxxx[k] * ab_x + g_z_0_y_xxxxxxx[k];

                g_z_0_xy_xxxxxy[k] = -g_z_0_y_xxxxxy[k] * ab_x + g_z_0_y_xxxxxxy[k];

                g_z_0_xy_xxxxxz[k] = -g_z_0_y_xxxxxz[k] * ab_x + g_z_0_y_xxxxxxz[k];

                g_z_0_xy_xxxxyy[k] = -g_z_0_y_xxxxyy[k] * ab_x + g_z_0_y_xxxxxyy[k];

                g_z_0_xy_xxxxyz[k] = -g_z_0_y_xxxxyz[k] * ab_x + g_z_0_y_xxxxxyz[k];

                g_z_0_xy_xxxxzz[k] = -g_z_0_y_xxxxzz[k] * ab_x + g_z_0_y_xxxxxzz[k];

                g_z_0_xy_xxxyyy[k] = -g_z_0_y_xxxyyy[k] * ab_x + g_z_0_y_xxxxyyy[k];

                g_z_0_xy_xxxyyz[k] = -g_z_0_y_xxxyyz[k] * ab_x + g_z_0_y_xxxxyyz[k];

                g_z_0_xy_xxxyzz[k] = -g_z_0_y_xxxyzz[k] * ab_x + g_z_0_y_xxxxyzz[k];

                g_z_0_xy_xxxzzz[k] = -g_z_0_y_xxxzzz[k] * ab_x + g_z_0_y_xxxxzzz[k];

                g_z_0_xy_xxyyyy[k] = -g_z_0_y_xxyyyy[k] * ab_x + g_z_0_y_xxxyyyy[k];

                g_z_0_xy_xxyyyz[k] = -g_z_0_y_xxyyyz[k] * ab_x + g_z_0_y_xxxyyyz[k];

                g_z_0_xy_xxyyzz[k] = -g_z_0_y_xxyyzz[k] * ab_x + g_z_0_y_xxxyyzz[k];

                g_z_0_xy_xxyzzz[k] = -g_z_0_y_xxyzzz[k] * ab_x + g_z_0_y_xxxyzzz[k];

                g_z_0_xy_xxzzzz[k] = -g_z_0_y_xxzzzz[k] * ab_x + g_z_0_y_xxxzzzz[k];

                g_z_0_xy_xyyyyy[k] = -g_z_0_y_xyyyyy[k] * ab_x + g_z_0_y_xxyyyyy[k];

                g_z_0_xy_xyyyyz[k] = -g_z_0_y_xyyyyz[k] * ab_x + g_z_0_y_xxyyyyz[k];

                g_z_0_xy_xyyyzz[k] = -g_z_0_y_xyyyzz[k] * ab_x + g_z_0_y_xxyyyzz[k];

                g_z_0_xy_xyyzzz[k] = -g_z_0_y_xyyzzz[k] * ab_x + g_z_0_y_xxyyzzz[k];

                g_z_0_xy_xyzzzz[k] = -g_z_0_y_xyzzzz[k] * ab_x + g_z_0_y_xxyzzzz[k];

                g_z_0_xy_xzzzzz[k] = -g_z_0_y_xzzzzz[k] * ab_x + g_z_0_y_xxzzzzz[k];

                g_z_0_xy_yyyyyy[k] = -g_z_0_y_yyyyyy[k] * ab_x + g_z_0_y_xyyyyyy[k];

                g_z_0_xy_yyyyyz[k] = -g_z_0_y_yyyyyz[k] * ab_x + g_z_0_y_xyyyyyz[k];

                g_z_0_xy_yyyyzz[k] = -g_z_0_y_yyyyzz[k] * ab_x + g_z_0_y_xyyyyzz[k];

                g_z_0_xy_yyyzzz[k] = -g_z_0_y_yyyzzz[k] * ab_x + g_z_0_y_xyyyzzz[k];

                g_z_0_xy_yyzzzz[k] = -g_z_0_y_yyzzzz[k] * ab_x + g_z_0_y_xyyzzzz[k];

                g_z_0_xy_yzzzzz[k] = -g_z_0_y_yzzzzz[k] * ab_x + g_z_0_y_xyzzzzz[k];

                g_z_0_xy_zzzzzz[k] = -g_z_0_y_zzzzzz[k] * ab_x + g_z_0_y_xzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_z_0_xz_xxxxxx = cbuffer.data(di_geom_10_off + 392 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxy = cbuffer.data(di_geom_10_off + 393 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxz = cbuffer.data(di_geom_10_off + 394 * ccomps * dcomps);

            auto g_z_0_xz_xxxxyy = cbuffer.data(di_geom_10_off + 395 * ccomps * dcomps);

            auto g_z_0_xz_xxxxyz = cbuffer.data(di_geom_10_off + 396 * ccomps * dcomps);

            auto g_z_0_xz_xxxxzz = cbuffer.data(di_geom_10_off + 397 * ccomps * dcomps);

            auto g_z_0_xz_xxxyyy = cbuffer.data(di_geom_10_off + 398 * ccomps * dcomps);

            auto g_z_0_xz_xxxyyz = cbuffer.data(di_geom_10_off + 399 * ccomps * dcomps);

            auto g_z_0_xz_xxxyzz = cbuffer.data(di_geom_10_off + 400 * ccomps * dcomps);

            auto g_z_0_xz_xxxzzz = cbuffer.data(di_geom_10_off + 401 * ccomps * dcomps);

            auto g_z_0_xz_xxyyyy = cbuffer.data(di_geom_10_off + 402 * ccomps * dcomps);

            auto g_z_0_xz_xxyyyz = cbuffer.data(di_geom_10_off + 403 * ccomps * dcomps);

            auto g_z_0_xz_xxyyzz = cbuffer.data(di_geom_10_off + 404 * ccomps * dcomps);

            auto g_z_0_xz_xxyzzz = cbuffer.data(di_geom_10_off + 405 * ccomps * dcomps);

            auto g_z_0_xz_xxzzzz = cbuffer.data(di_geom_10_off + 406 * ccomps * dcomps);

            auto g_z_0_xz_xyyyyy = cbuffer.data(di_geom_10_off + 407 * ccomps * dcomps);

            auto g_z_0_xz_xyyyyz = cbuffer.data(di_geom_10_off + 408 * ccomps * dcomps);

            auto g_z_0_xz_xyyyzz = cbuffer.data(di_geom_10_off + 409 * ccomps * dcomps);

            auto g_z_0_xz_xyyzzz = cbuffer.data(di_geom_10_off + 410 * ccomps * dcomps);

            auto g_z_0_xz_xyzzzz = cbuffer.data(di_geom_10_off + 411 * ccomps * dcomps);

            auto g_z_0_xz_xzzzzz = cbuffer.data(di_geom_10_off + 412 * ccomps * dcomps);

            auto g_z_0_xz_yyyyyy = cbuffer.data(di_geom_10_off + 413 * ccomps * dcomps);

            auto g_z_0_xz_yyyyyz = cbuffer.data(di_geom_10_off + 414 * ccomps * dcomps);

            auto g_z_0_xz_yyyyzz = cbuffer.data(di_geom_10_off + 415 * ccomps * dcomps);

            auto g_z_0_xz_yyyzzz = cbuffer.data(di_geom_10_off + 416 * ccomps * dcomps);

            auto g_z_0_xz_yyzzzz = cbuffer.data(di_geom_10_off + 417 * ccomps * dcomps);

            auto g_z_0_xz_yzzzzz = cbuffer.data(di_geom_10_off + 418 * ccomps * dcomps);

            auto g_z_0_xz_zzzzzz = cbuffer.data(di_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xz_xxxxxx, g_z_0_xz_xxxxxy, g_z_0_xz_xxxxxz, g_z_0_xz_xxxxyy, g_z_0_xz_xxxxyz, g_z_0_xz_xxxxzz, g_z_0_xz_xxxyyy, g_z_0_xz_xxxyyz, g_z_0_xz_xxxyzz, g_z_0_xz_xxxzzz, g_z_0_xz_xxyyyy, g_z_0_xz_xxyyyz, g_z_0_xz_xxyyzz, g_z_0_xz_xxyzzz, g_z_0_xz_xxzzzz, g_z_0_xz_xyyyyy, g_z_0_xz_xyyyyz, g_z_0_xz_xyyyzz, g_z_0_xz_xyyzzz, g_z_0_xz_xyzzzz, g_z_0_xz_xzzzzz, g_z_0_xz_yyyyyy, g_z_0_xz_yyyyyz, g_z_0_xz_yyyyzz, g_z_0_xz_yyyzzz, g_z_0_xz_yyzzzz, g_z_0_xz_yzzzzz, g_z_0_xz_zzzzzz, g_z_0_z_xxxxxx, g_z_0_z_xxxxxxx, g_z_0_z_xxxxxxy, g_z_0_z_xxxxxxz, g_z_0_z_xxxxxy, g_z_0_z_xxxxxyy, g_z_0_z_xxxxxyz, g_z_0_z_xxxxxz, g_z_0_z_xxxxxzz, g_z_0_z_xxxxyy, g_z_0_z_xxxxyyy, g_z_0_z_xxxxyyz, g_z_0_z_xxxxyz, g_z_0_z_xxxxyzz, g_z_0_z_xxxxzz, g_z_0_z_xxxxzzz, g_z_0_z_xxxyyy, g_z_0_z_xxxyyyy, g_z_0_z_xxxyyyz, g_z_0_z_xxxyyz, g_z_0_z_xxxyyzz, g_z_0_z_xxxyzz, g_z_0_z_xxxyzzz, g_z_0_z_xxxzzz, g_z_0_z_xxxzzzz, g_z_0_z_xxyyyy, g_z_0_z_xxyyyyy, g_z_0_z_xxyyyyz, g_z_0_z_xxyyyz, g_z_0_z_xxyyyzz, g_z_0_z_xxyyzz, g_z_0_z_xxyyzzz, g_z_0_z_xxyzzz, g_z_0_z_xxyzzzz, g_z_0_z_xxzzzz, g_z_0_z_xxzzzzz, g_z_0_z_xyyyyy, g_z_0_z_xyyyyyy, g_z_0_z_xyyyyyz, g_z_0_z_xyyyyz, g_z_0_z_xyyyyzz, g_z_0_z_xyyyzz, g_z_0_z_xyyyzzz, g_z_0_z_xyyzzz, g_z_0_z_xyyzzzz, g_z_0_z_xyzzzz, g_z_0_z_xyzzzzz, g_z_0_z_xzzzzz, g_z_0_z_xzzzzzz, g_z_0_z_yyyyyy, g_z_0_z_yyyyyz, g_z_0_z_yyyyzz, g_z_0_z_yyyzzz, g_z_0_z_yyzzzz, g_z_0_z_yzzzzz, g_z_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xz_xxxxxx[k] = -g_z_0_z_xxxxxx[k] * ab_x + g_z_0_z_xxxxxxx[k];

                g_z_0_xz_xxxxxy[k] = -g_z_0_z_xxxxxy[k] * ab_x + g_z_0_z_xxxxxxy[k];

                g_z_0_xz_xxxxxz[k] = -g_z_0_z_xxxxxz[k] * ab_x + g_z_0_z_xxxxxxz[k];

                g_z_0_xz_xxxxyy[k] = -g_z_0_z_xxxxyy[k] * ab_x + g_z_0_z_xxxxxyy[k];

                g_z_0_xz_xxxxyz[k] = -g_z_0_z_xxxxyz[k] * ab_x + g_z_0_z_xxxxxyz[k];

                g_z_0_xz_xxxxzz[k] = -g_z_0_z_xxxxzz[k] * ab_x + g_z_0_z_xxxxxzz[k];

                g_z_0_xz_xxxyyy[k] = -g_z_0_z_xxxyyy[k] * ab_x + g_z_0_z_xxxxyyy[k];

                g_z_0_xz_xxxyyz[k] = -g_z_0_z_xxxyyz[k] * ab_x + g_z_0_z_xxxxyyz[k];

                g_z_0_xz_xxxyzz[k] = -g_z_0_z_xxxyzz[k] * ab_x + g_z_0_z_xxxxyzz[k];

                g_z_0_xz_xxxzzz[k] = -g_z_0_z_xxxzzz[k] * ab_x + g_z_0_z_xxxxzzz[k];

                g_z_0_xz_xxyyyy[k] = -g_z_0_z_xxyyyy[k] * ab_x + g_z_0_z_xxxyyyy[k];

                g_z_0_xz_xxyyyz[k] = -g_z_0_z_xxyyyz[k] * ab_x + g_z_0_z_xxxyyyz[k];

                g_z_0_xz_xxyyzz[k] = -g_z_0_z_xxyyzz[k] * ab_x + g_z_0_z_xxxyyzz[k];

                g_z_0_xz_xxyzzz[k] = -g_z_0_z_xxyzzz[k] * ab_x + g_z_0_z_xxxyzzz[k];

                g_z_0_xz_xxzzzz[k] = -g_z_0_z_xxzzzz[k] * ab_x + g_z_0_z_xxxzzzz[k];

                g_z_0_xz_xyyyyy[k] = -g_z_0_z_xyyyyy[k] * ab_x + g_z_0_z_xxyyyyy[k];

                g_z_0_xz_xyyyyz[k] = -g_z_0_z_xyyyyz[k] * ab_x + g_z_0_z_xxyyyyz[k];

                g_z_0_xz_xyyyzz[k] = -g_z_0_z_xyyyzz[k] * ab_x + g_z_0_z_xxyyyzz[k];

                g_z_0_xz_xyyzzz[k] = -g_z_0_z_xyyzzz[k] * ab_x + g_z_0_z_xxyyzzz[k];

                g_z_0_xz_xyzzzz[k] = -g_z_0_z_xyzzzz[k] * ab_x + g_z_0_z_xxyzzzz[k];

                g_z_0_xz_xzzzzz[k] = -g_z_0_z_xzzzzz[k] * ab_x + g_z_0_z_xxzzzzz[k];

                g_z_0_xz_yyyyyy[k] = -g_z_0_z_yyyyyy[k] * ab_x + g_z_0_z_xyyyyyy[k];

                g_z_0_xz_yyyyyz[k] = -g_z_0_z_yyyyyz[k] * ab_x + g_z_0_z_xyyyyyz[k];

                g_z_0_xz_yyyyzz[k] = -g_z_0_z_yyyyzz[k] * ab_x + g_z_0_z_xyyyyzz[k];

                g_z_0_xz_yyyzzz[k] = -g_z_0_z_yyyzzz[k] * ab_x + g_z_0_z_xyyyzzz[k];

                g_z_0_xz_yyzzzz[k] = -g_z_0_z_yyzzzz[k] * ab_x + g_z_0_z_xyyzzzz[k];

                g_z_0_xz_yzzzzz[k] = -g_z_0_z_yzzzzz[k] * ab_x + g_z_0_z_xyzzzzz[k];

                g_z_0_xz_zzzzzz[k] = -g_z_0_z_zzzzzz[k] * ab_x + g_z_0_z_xzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_z_0_yy_xxxxxx = cbuffer.data(di_geom_10_off + 420 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxy = cbuffer.data(di_geom_10_off + 421 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxz = cbuffer.data(di_geom_10_off + 422 * ccomps * dcomps);

            auto g_z_0_yy_xxxxyy = cbuffer.data(di_geom_10_off + 423 * ccomps * dcomps);

            auto g_z_0_yy_xxxxyz = cbuffer.data(di_geom_10_off + 424 * ccomps * dcomps);

            auto g_z_0_yy_xxxxzz = cbuffer.data(di_geom_10_off + 425 * ccomps * dcomps);

            auto g_z_0_yy_xxxyyy = cbuffer.data(di_geom_10_off + 426 * ccomps * dcomps);

            auto g_z_0_yy_xxxyyz = cbuffer.data(di_geom_10_off + 427 * ccomps * dcomps);

            auto g_z_0_yy_xxxyzz = cbuffer.data(di_geom_10_off + 428 * ccomps * dcomps);

            auto g_z_0_yy_xxxzzz = cbuffer.data(di_geom_10_off + 429 * ccomps * dcomps);

            auto g_z_0_yy_xxyyyy = cbuffer.data(di_geom_10_off + 430 * ccomps * dcomps);

            auto g_z_0_yy_xxyyyz = cbuffer.data(di_geom_10_off + 431 * ccomps * dcomps);

            auto g_z_0_yy_xxyyzz = cbuffer.data(di_geom_10_off + 432 * ccomps * dcomps);

            auto g_z_0_yy_xxyzzz = cbuffer.data(di_geom_10_off + 433 * ccomps * dcomps);

            auto g_z_0_yy_xxzzzz = cbuffer.data(di_geom_10_off + 434 * ccomps * dcomps);

            auto g_z_0_yy_xyyyyy = cbuffer.data(di_geom_10_off + 435 * ccomps * dcomps);

            auto g_z_0_yy_xyyyyz = cbuffer.data(di_geom_10_off + 436 * ccomps * dcomps);

            auto g_z_0_yy_xyyyzz = cbuffer.data(di_geom_10_off + 437 * ccomps * dcomps);

            auto g_z_0_yy_xyyzzz = cbuffer.data(di_geom_10_off + 438 * ccomps * dcomps);

            auto g_z_0_yy_xyzzzz = cbuffer.data(di_geom_10_off + 439 * ccomps * dcomps);

            auto g_z_0_yy_xzzzzz = cbuffer.data(di_geom_10_off + 440 * ccomps * dcomps);

            auto g_z_0_yy_yyyyyy = cbuffer.data(di_geom_10_off + 441 * ccomps * dcomps);

            auto g_z_0_yy_yyyyyz = cbuffer.data(di_geom_10_off + 442 * ccomps * dcomps);

            auto g_z_0_yy_yyyyzz = cbuffer.data(di_geom_10_off + 443 * ccomps * dcomps);

            auto g_z_0_yy_yyyzzz = cbuffer.data(di_geom_10_off + 444 * ccomps * dcomps);

            auto g_z_0_yy_yyzzzz = cbuffer.data(di_geom_10_off + 445 * ccomps * dcomps);

            auto g_z_0_yy_yzzzzz = cbuffer.data(di_geom_10_off + 446 * ccomps * dcomps);

            auto g_z_0_yy_zzzzzz = cbuffer.data(di_geom_10_off + 447 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_xxxxxx, g_z_0_y_xxxxxxy, g_z_0_y_xxxxxy, g_z_0_y_xxxxxyy, g_z_0_y_xxxxxyz, g_z_0_y_xxxxxz, g_z_0_y_xxxxyy, g_z_0_y_xxxxyyy, g_z_0_y_xxxxyyz, g_z_0_y_xxxxyz, g_z_0_y_xxxxyzz, g_z_0_y_xxxxzz, g_z_0_y_xxxyyy, g_z_0_y_xxxyyyy, g_z_0_y_xxxyyyz, g_z_0_y_xxxyyz, g_z_0_y_xxxyyzz, g_z_0_y_xxxyzz, g_z_0_y_xxxyzzz, g_z_0_y_xxxzzz, g_z_0_y_xxyyyy, g_z_0_y_xxyyyyy, g_z_0_y_xxyyyyz, g_z_0_y_xxyyyz, g_z_0_y_xxyyyzz, g_z_0_y_xxyyzz, g_z_0_y_xxyyzzz, g_z_0_y_xxyzzz, g_z_0_y_xxyzzzz, g_z_0_y_xxzzzz, g_z_0_y_xyyyyy, g_z_0_y_xyyyyyy, g_z_0_y_xyyyyyz, g_z_0_y_xyyyyz, g_z_0_y_xyyyyzz, g_z_0_y_xyyyzz, g_z_0_y_xyyyzzz, g_z_0_y_xyyzzz, g_z_0_y_xyyzzzz, g_z_0_y_xyzzzz, g_z_0_y_xyzzzzz, g_z_0_y_xzzzzz, g_z_0_y_yyyyyy, g_z_0_y_yyyyyyy, g_z_0_y_yyyyyyz, g_z_0_y_yyyyyz, g_z_0_y_yyyyyzz, g_z_0_y_yyyyzz, g_z_0_y_yyyyzzz, g_z_0_y_yyyzzz, g_z_0_y_yyyzzzz, g_z_0_y_yyzzzz, g_z_0_y_yyzzzzz, g_z_0_y_yzzzzz, g_z_0_y_yzzzzzz, g_z_0_y_zzzzzz, g_z_0_yy_xxxxxx, g_z_0_yy_xxxxxy, g_z_0_yy_xxxxxz, g_z_0_yy_xxxxyy, g_z_0_yy_xxxxyz, g_z_0_yy_xxxxzz, g_z_0_yy_xxxyyy, g_z_0_yy_xxxyyz, g_z_0_yy_xxxyzz, g_z_0_yy_xxxzzz, g_z_0_yy_xxyyyy, g_z_0_yy_xxyyyz, g_z_0_yy_xxyyzz, g_z_0_yy_xxyzzz, g_z_0_yy_xxzzzz, g_z_0_yy_xyyyyy, g_z_0_yy_xyyyyz, g_z_0_yy_xyyyzz, g_z_0_yy_xyyzzz, g_z_0_yy_xyzzzz, g_z_0_yy_xzzzzz, g_z_0_yy_yyyyyy, g_z_0_yy_yyyyyz, g_z_0_yy_yyyyzz, g_z_0_yy_yyyzzz, g_z_0_yy_yyzzzz, g_z_0_yy_yzzzzz, g_z_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yy_xxxxxx[k] = -g_z_0_y_xxxxxx[k] * ab_y + g_z_0_y_xxxxxxy[k];

                g_z_0_yy_xxxxxy[k] = -g_z_0_y_xxxxxy[k] * ab_y + g_z_0_y_xxxxxyy[k];

                g_z_0_yy_xxxxxz[k] = -g_z_0_y_xxxxxz[k] * ab_y + g_z_0_y_xxxxxyz[k];

                g_z_0_yy_xxxxyy[k] = -g_z_0_y_xxxxyy[k] * ab_y + g_z_0_y_xxxxyyy[k];

                g_z_0_yy_xxxxyz[k] = -g_z_0_y_xxxxyz[k] * ab_y + g_z_0_y_xxxxyyz[k];

                g_z_0_yy_xxxxzz[k] = -g_z_0_y_xxxxzz[k] * ab_y + g_z_0_y_xxxxyzz[k];

                g_z_0_yy_xxxyyy[k] = -g_z_0_y_xxxyyy[k] * ab_y + g_z_0_y_xxxyyyy[k];

                g_z_0_yy_xxxyyz[k] = -g_z_0_y_xxxyyz[k] * ab_y + g_z_0_y_xxxyyyz[k];

                g_z_0_yy_xxxyzz[k] = -g_z_0_y_xxxyzz[k] * ab_y + g_z_0_y_xxxyyzz[k];

                g_z_0_yy_xxxzzz[k] = -g_z_0_y_xxxzzz[k] * ab_y + g_z_0_y_xxxyzzz[k];

                g_z_0_yy_xxyyyy[k] = -g_z_0_y_xxyyyy[k] * ab_y + g_z_0_y_xxyyyyy[k];

                g_z_0_yy_xxyyyz[k] = -g_z_0_y_xxyyyz[k] * ab_y + g_z_0_y_xxyyyyz[k];

                g_z_0_yy_xxyyzz[k] = -g_z_0_y_xxyyzz[k] * ab_y + g_z_0_y_xxyyyzz[k];

                g_z_0_yy_xxyzzz[k] = -g_z_0_y_xxyzzz[k] * ab_y + g_z_0_y_xxyyzzz[k];

                g_z_0_yy_xxzzzz[k] = -g_z_0_y_xxzzzz[k] * ab_y + g_z_0_y_xxyzzzz[k];

                g_z_0_yy_xyyyyy[k] = -g_z_0_y_xyyyyy[k] * ab_y + g_z_0_y_xyyyyyy[k];

                g_z_0_yy_xyyyyz[k] = -g_z_0_y_xyyyyz[k] * ab_y + g_z_0_y_xyyyyyz[k];

                g_z_0_yy_xyyyzz[k] = -g_z_0_y_xyyyzz[k] * ab_y + g_z_0_y_xyyyyzz[k];

                g_z_0_yy_xyyzzz[k] = -g_z_0_y_xyyzzz[k] * ab_y + g_z_0_y_xyyyzzz[k];

                g_z_0_yy_xyzzzz[k] = -g_z_0_y_xyzzzz[k] * ab_y + g_z_0_y_xyyzzzz[k];

                g_z_0_yy_xzzzzz[k] = -g_z_0_y_xzzzzz[k] * ab_y + g_z_0_y_xyzzzzz[k];

                g_z_0_yy_yyyyyy[k] = -g_z_0_y_yyyyyy[k] * ab_y + g_z_0_y_yyyyyyy[k];

                g_z_0_yy_yyyyyz[k] = -g_z_0_y_yyyyyz[k] * ab_y + g_z_0_y_yyyyyyz[k];

                g_z_0_yy_yyyyzz[k] = -g_z_0_y_yyyyzz[k] * ab_y + g_z_0_y_yyyyyzz[k];

                g_z_0_yy_yyyzzz[k] = -g_z_0_y_yyyzzz[k] * ab_y + g_z_0_y_yyyyzzz[k];

                g_z_0_yy_yyzzzz[k] = -g_z_0_y_yyzzzz[k] * ab_y + g_z_0_y_yyyzzzz[k];

                g_z_0_yy_yzzzzz[k] = -g_z_0_y_yzzzzz[k] * ab_y + g_z_0_y_yyzzzzz[k];

                g_z_0_yy_zzzzzz[k] = -g_z_0_y_zzzzzz[k] * ab_y + g_z_0_y_yzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_z_0_yz_xxxxxx = cbuffer.data(di_geom_10_off + 448 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxy = cbuffer.data(di_geom_10_off + 449 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxz = cbuffer.data(di_geom_10_off + 450 * ccomps * dcomps);

            auto g_z_0_yz_xxxxyy = cbuffer.data(di_geom_10_off + 451 * ccomps * dcomps);

            auto g_z_0_yz_xxxxyz = cbuffer.data(di_geom_10_off + 452 * ccomps * dcomps);

            auto g_z_0_yz_xxxxzz = cbuffer.data(di_geom_10_off + 453 * ccomps * dcomps);

            auto g_z_0_yz_xxxyyy = cbuffer.data(di_geom_10_off + 454 * ccomps * dcomps);

            auto g_z_0_yz_xxxyyz = cbuffer.data(di_geom_10_off + 455 * ccomps * dcomps);

            auto g_z_0_yz_xxxyzz = cbuffer.data(di_geom_10_off + 456 * ccomps * dcomps);

            auto g_z_0_yz_xxxzzz = cbuffer.data(di_geom_10_off + 457 * ccomps * dcomps);

            auto g_z_0_yz_xxyyyy = cbuffer.data(di_geom_10_off + 458 * ccomps * dcomps);

            auto g_z_0_yz_xxyyyz = cbuffer.data(di_geom_10_off + 459 * ccomps * dcomps);

            auto g_z_0_yz_xxyyzz = cbuffer.data(di_geom_10_off + 460 * ccomps * dcomps);

            auto g_z_0_yz_xxyzzz = cbuffer.data(di_geom_10_off + 461 * ccomps * dcomps);

            auto g_z_0_yz_xxzzzz = cbuffer.data(di_geom_10_off + 462 * ccomps * dcomps);

            auto g_z_0_yz_xyyyyy = cbuffer.data(di_geom_10_off + 463 * ccomps * dcomps);

            auto g_z_0_yz_xyyyyz = cbuffer.data(di_geom_10_off + 464 * ccomps * dcomps);

            auto g_z_0_yz_xyyyzz = cbuffer.data(di_geom_10_off + 465 * ccomps * dcomps);

            auto g_z_0_yz_xyyzzz = cbuffer.data(di_geom_10_off + 466 * ccomps * dcomps);

            auto g_z_0_yz_xyzzzz = cbuffer.data(di_geom_10_off + 467 * ccomps * dcomps);

            auto g_z_0_yz_xzzzzz = cbuffer.data(di_geom_10_off + 468 * ccomps * dcomps);

            auto g_z_0_yz_yyyyyy = cbuffer.data(di_geom_10_off + 469 * ccomps * dcomps);

            auto g_z_0_yz_yyyyyz = cbuffer.data(di_geom_10_off + 470 * ccomps * dcomps);

            auto g_z_0_yz_yyyyzz = cbuffer.data(di_geom_10_off + 471 * ccomps * dcomps);

            auto g_z_0_yz_yyyzzz = cbuffer.data(di_geom_10_off + 472 * ccomps * dcomps);

            auto g_z_0_yz_yyzzzz = cbuffer.data(di_geom_10_off + 473 * ccomps * dcomps);

            auto g_z_0_yz_yzzzzz = cbuffer.data(di_geom_10_off + 474 * ccomps * dcomps);

            auto g_z_0_yz_zzzzzz = cbuffer.data(di_geom_10_off + 475 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yz_xxxxxx, g_z_0_yz_xxxxxy, g_z_0_yz_xxxxxz, g_z_0_yz_xxxxyy, g_z_0_yz_xxxxyz, g_z_0_yz_xxxxzz, g_z_0_yz_xxxyyy, g_z_0_yz_xxxyyz, g_z_0_yz_xxxyzz, g_z_0_yz_xxxzzz, g_z_0_yz_xxyyyy, g_z_0_yz_xxyyyz, g_z_0_yz_xxyyzz, g_z_0_yz_xxyzzz, g_z_0_yz_xxzzzz, g_z_0_yz_xyyyyy, g_z_0_yz_xyyyyz, g_z_0_yz_xyyyzz, g_z_0_yz_xyyzzz, g_z_0_yz_xyzzzz, g_z_0_yz_xzzzzz, g_z_0_yz_yyyyyy, g_z_0_yz_yyyyyz, g_z_0_yz_yyyyzz, g_z_0_yz_yyyzzz, g_z_0_yz_yyzzzz, g_z_0_yz_yzzzzz, g_z_0_yz_zzzzzz, g_z_0_z_xxxxxx, g_z_0_z_xxxxxxy, g_z_0_z_xxxxxy, g_z_0_z_xxxxxyy, g_z_0_z_xxxxxyz, g_z_0_z_xxxxxz, g_z_0_z_xxxxyy, g_z_0_z_xxxxyyy, g_z_0_z_xxxxyyz, g_z_0_z_xxxxyz, g_z_0_z_xxxxyzz, g_z_0_z_xxxxzz, g_z_0_z_xxxyyy, g_z_0_z_xxxyyyy, g_z_0_z_xxxyyyz, g_z_0_z_xxxyyz, g_z_0_z_xxxyyzz, g_z_0_z_xxxyzz, g_z_0_z_xxxyzzz, g_z_0_z_xxxzzz, g_z_0_z_xxyyyy, g_z_0_z_xxyyyyy, g_z_0_z_xxyyyyz, g_z_0_z_xxyyyz, g_z_0_z_xxyyyzz, g_z_0_z_xxyyzz, g_z_0_z_xxyyzzz, g_z_0_z_xxyzzz, g_z_0_z_xxyzzzz, g_z_0_z_xxzzzz, g_z_0_z_xyyyyy, g_z_0_z_xyyyyyy, g_z_0_z_xyyyyyz, g_z_0_z_xyyyyz, g_z_0_z_xyyyyzz, g_z_0_z_xyyyzz, g_z_0_z_xyyyzzz, g_z_0_z_xyyzzz, g_z_0_z_xyyzzzz, g_z_0_z_xyzzzz, g_z_0_z_xyzzzzz, g_z_0_z_xzzzzz, g_z_0_z_yyyyyy, g_z_0_z_yyyyyyy, g_z_0_z_yyyyyyz, g_z_0_z_yyyyyz, g_z_0_z_yyyyyzz, g_z_0_z_yyyyzz, g_z_0_z_yyyyzzz, g_z_0_z_yyyzzz, g_z_0_z_yyyzzzz, g_z_0_z_yyzzzz, g_z_0_z_yyzzzzz, g_z_0_z_yzzzzz, g_z_0_z_yzzzzzz, g_z_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yz_xxxxxx[k] = -g_z_0_z_xxxxxx[k] * ab_y + g_z_0_z_xxxxxxy[k];

                g_z_0_yz_xxxxxy[k] = -g_z_0_z_xxxxxy[k] * ab_y + g_z_0_z_xxxxxyy[k];

                g_z_0_yz_xxxxxz[k] = -g_z_0_z_xxxxxz[k] * ab_y + g_z_0_z_xxxxxyz[k];

                g_z_0_yz_xxxxyy[k] = -g_z_0_z_xxxxyy[k] * ab_y + g_z_0_z_xxxxyyy[k];

                g_z_0_yz_xxxxyz[k] = -g_z_0_z_xxxxyz[k] * ab_y + g_z_0_z_xxxxyyz[k];

                g_z_0_yz_xxxxzz[k] = -g_z_0_z_xxxxzz[k] * ab_y + g_z_0_z_xxxxyzz[k];

                g_z_0_yz_xxxyyy[k] = -g_z_0_z_xxxyyy[k] * ab_y + g_z_0_z_xxxyyyy[k];

                g_z_0_yz_xxxyyz[k] = -g_z_0_z_xxxyyz[k] * ab_y + g_z_0_z_xxxyyyz[k];

                g_z_0_yz_xxxyzz[k] = -g_z_0_z_xxxyzz[k] * ab_y + g_z_0_z_xxxyyzz[k];

                g_z_0_yz_xxxzzz[k] = -g_z_0_z_xxxzzz[k] * ab_y + g_z_0_z_xxxyzzz[k];

                g_z_0_yz_xxyyyy[k] = -g_z_0_z_xxyyyy[k] * ab_y + g_z_0_z_xxyyyyy[k];

                g_z_0_yz_xxyyyz[k] = -g_z_0_z_xxyyyz[k] * ab_y + g_z_0_z_xxyyyyz[k];

                g_z_0_yz_xxyyzz[k] = -g_z_0_z_xxyyzz[k] * ab_y + g_z_0_z_xxyyyzz[k];

                g_z_0_yz_xxyzzz[k] = -g_z_0_z_xxyzzz[k] * ab_y + g_z_0_z_xxyyzzz[k];

                g_z_0_yz_xxzzzz[k] = -g_z_0_z_xxzzzz[k] * ab_y + g_z_0_z_xxyzzzz[k];

                g_z_0_yz_xyyyyy[k] = -g_z_0_z_xyyyyy[k] * ab_y + g_z_0_z_xyyyyyy[k];

                g_z_0_yz_xyyyyz[k] = -g_z_0_z_xyyyyz[k] * ab_y + g_z_0_z_xyyyyyz[k];

                g_z_0_yz_xyyyzz[k] = -g_z_0_z_xyyyzz[k] * ab_y + g_z_0_z_xyyyyzz[k];

                g_z_0_yz_xyyzzz[k] = -g_z_0_z_xyyzzz[k] * ab_y + g_z_0_z_xyyyzzz[k];

                g_z_0_yz_xyzzzz[k] = -g_z_0_z_xyzzzz[k] * ab_y + g_z_0_z_xyyzzzz[k];

                g_z_0_yz_xzzzzz[k] = -g_z_0_z_xzzzzz[k] * ab_y + g_z_0_z_xyzzzzz[k];

                g_z_0_yz_yyyyyy[k] = -g_z_0_z_yyyyyy[k] * ab_y + g_z_0_z_yyyyyyy[k];

                g_z_0_yz_yyyyyz[k] = -g_z_0_z_yyyyyz[k] * ab_y + g_z_0_z_yyyyyyz[k];

                g_z_0_yz_yyyyzz[k] = -g_z_0_z_yyyyzz[k] * ab_y + g_z_0_z_yyyyyzz[k];

                g_z_0_yz_yyyzzz[k] = -g_z_0_z_yyyzzz[k] * ab_y + g_z_0_z_yyyyzzz[k];

                g_z_0_yz_yyzzzz[k] = -g_z_0_z_yyzzzz[k] * ab_y + g_z_0_z_yyyzzzz[k];

                g_z_0_yz_yzzzzz[k] = -g_z_0_z_yzzzzz[k] * ab_y + g_z_0_z_yyzzzzz[k];

                g_z_0_yz_zzzzzz[k] = -g_z_0_z_zzzzzz[k] * ab_y + g_z_0_z_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_z_0_zz_xxxxxx = cbuffer.data(di_geom_10_off + 476 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxy = cbuffer.data(di_geom_10_off + 477 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxz = cbuffer.data(di_geom_10_off + 478 * ccomps * dcomps);

            auto g_z_0_zz_xxxxyy = cbuffer.data(di_geom_10_off + 479 * ccomps * dcomps);

            auto g_z_0_zz_xxxxyz = cbuffer.data(di_geom_10_off + 480 * ccomps * dcomps);

            auto g_z_0_zz_xxxxzz = cbuffer.data(di_geom_10_off + 481 * ccomps * dcomps);

            auto g_z_0_zz_xxxyyy = cbuffer.data(di_geom_10_off + 482 * ccomps * dcomps);

            auto g_z_0_zz_xxxyyz = cbuffer.data(di_geom_10_off + 483 * ccomps * dcomps);

            auto g_z_0_zz_xxxyzz = cbuffer.data(di_geom_10_off + 484 * ccomps * dcomps);

            auto g_z_0_zz_xxxzzz = cbuffer.data(di_geom_10_off + 485 * ccomps * dcomps);

            auto g_z_0_zz_xxyyyy = cbuffer.data(di_geom_10_off + 486 * ccomps * dcomps);

            auto g_z_0_zz_xxyyyz = cbuffer.data(di_geom_10_off + 487 * ccomps * dcomps);

            auto g_z_0_zz_xxyyzz = cbuffer.data(di_geom_10_off + 488 * ccomps * dcomps);

            auto g_z_0_zz_xxyzzz = cbuffer.data(di_geom_10_off + 489 * ccomps * dcomps);

            auto g_z_0_zz_xxzzzz = cbuffer.data(di_geom_10_off + 490 * ccomps * dcomps);

            auto g_z_0_zz_xyyyyy = cbuffer.data(di_geom_10_off + 491 * ccomps * dcomps);

            auto g_z_0_zz_xyyyyz = cbuffer.data(di_geom_10_off + 492 * ccomps * dcomps);

            auto g_z_0_zz_xyyyzz = cbuffer.data(di_geom_10_off + 493 * ccomps * dcomps);

            auto g_z_0_zz_xyyzzz = cbuffer.data(di_geom_10_off + 494 * ccomps * dcomps);

            auto g_z_0_zz_xyzzzz = cbuffer.data(di_geom_10_off + 495 * ccomps * dcomps);

            auto g_z_0_zz_xzzzzz = cbuffer.data(di_geom_10_off + 496 * ccomps * dcomps);

            auto g_z_0_zz_yyyyyy = cbuffer.data(di_geom_10_off + 497 * ccomps * dcomps);

            auto g_z_0_zz_yyyyyz = cbuffer.data(di_geom_10_off + 498 * ccomps * dcomps);

            auto g_z_0_zz_yyyyzz = cbuffer.data(di_geom_10_off + 499 * ccomps * dcomps);

            auto g_z_0_zz_yyyzzz = cbuffer.data(di_geom_10_off + 500 * ccomps * dcomps);

            auto g_z_0_zz_yyzzzz = cbuffer.data(di_geom_10_off + 501 * ccomps * dcomps);

            auto g_z_0_zz_yzzzzz = cbuffer.data(di_geom_10_off + 502 * ccomps * dcomps);

            auto g_z_0_zz_zzzzzz = cbuffer.data(di_geom_10_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxxxx, g_z_0_z_xxxxxxz, g_z_0_z_xxxxxy, g_z_0_z_xxxxxyz, g_z_0_z_xxxxxz, g_z_0_z_xxxxxzz, g_z_0_z_xxxxyy, g_z_0_z_xxxxyyz, g_z_0_z_xxxxyz, g_z_0_z_xxxxyzz, g_z_0_z_xxxxzz, g_z_0_z_xxxxzzz, g_z_0_z_xxxyyy, g_z_0_z_xxxyyyz, g_z_0_z_xxxyyz, g_z_0_z_xxxyyzz, g_z_0_z_xxxyzz, g_z_0_z_xxxyzzz, g_z_0_z_xxxzzz, g_z_0_z_xxxzzzz, g_z_0_z_xxyyyy, g_z_0_z_xxyyyyz, g_z_0_z_xxyyyz, g_z_0_z_xxyyyzz, g_z_0_z_xxyyzz, g_z_0_z_xxyyzzz, g_z_0_z_xxyzzz, g_z_0_z_xxyzzzz, g_z_0_z_xxzzzz, g_z_0_z_xxzzzzz, g_z_0_z_xyyyyy, g_z_0_z_xyyyyyz, g_z_0_z_xyyyyz, g_z_0_z_xyyyyzz, g_z_0_z_xyyyzz, g_z_0_z_xyyyzzz, g_z_0_z_xyyzzz, g_z_0_z_xyyzzzz, g_z_0_z_xyzzzz, g_z_0_z_xyzzzzz, g_z_0_z_xzzzzz, g_z_0_z_xzzzzzz, g_z_0_z_yyyyyy, g_z_0_z_yyyyyyz, g_z_0_z_yyyyyz, g_z_0_z_yyyyyzz, g_z_0_z_yyyyzz, g_z_0_z_yyyyzzz, g_z_0_z_yyyzzz, g_z_0_z_yyyzzzz, g_z_0_z_yyzzzz, g_z_0_z_yyzzzzz, g_z_0_z_yzzzzz, g_z_0_z_yzzzzzz, g_z_0_z_zzzzzz, g_z_0_z_zzzzzzz, g_z_0_zz_xxxxxx, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_yyyyyy, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_zzzzzz, g_z_xxxxxx, g_z_xxxxxy, g_z_xxxxxz, g_z_xxxxyy, g_z_xxxxyz, g_z_xxxxzz, g_z_xxxyyy, g_z_xxxyyz, g_z_xxxyzz, g_z_xxxzzz, g_z_xxyyyy, g_z_xxyyyz, g_z_xxyyzz, g_z_xxyzzz, g_z_xxzzzz, g_z_xyyyyy, g_z_xyyyyz, g_z_xyyyzz, g_z_xyyzzz, g_z_xyzzzz, g_z_xzzzzz, g_z_yyyyyy, g_z_yyyyyz, g_z_yyyyzz, g_z_yyyzzz, g_z_yyzzzz, g_z_yzzzzz, g_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zz_xxxxxx[k] = -g_z_xxxxxx[k] - g_z_0_z_xxxxxx[k] * ab_z + g_z_0_z_xxxxxxz[k];

                g_z_0_zz_xxxxxy[k] = -g_z_xxxxxy[k] - g_z_0_z_xxxxxy[k] * ab_z + g_z_0_z_xxxxxyz[k];

                g_z_0_zz_xxxxxz[k] = -g_z_xxxxxz[k] - g_z_0_z_xxxxxz[k] * ab_z + g_z_0_z_xxxxxzz[k];

                g_z_0_zz_xxxxyy[k] = -g_z_xxxxyy[k] - g_z_0_z_xxxxyy[k] * ab_z + g_z_0_z_xxxxyyz[k];

                g_z_0_zz_xxxxyz[k] = -g_z_xxxxyz[k] - g_z_0_z_xxxxyz[k] * ab_z + g_z_0_z_xxxxyzz[k];

                g_z_0_zz_xxxxzz[k] = -g_z_xxxxzz[k] - g_z_0_z_xxxxzz[k] * ab_z + g_z_0_z_xxxxzzz[k];

                g_z_0_zz_xxxyyy[k] = -g_z_xxxyyy[k] - g_z_0_z_xxxyyy[k] * ab_z + g_z_0_z_xxxyyyz[k];

                g_z_0_zz_xxxyyz[k] = -g_z_xxxyyz[k] - g_z_0_z_xxxyyz[k] * ab_z + g_z_0_z_xxxyyzz[k];

                g_z_0_zz_xxxyzz[k] = -g_z_xxxyzz[k] - g_z_0_z_xxxyzz[k] * ab_z + g_z_0_z_xxxyzzz[k];

                g_z_0_zz_xxxzzz[k] = -g_z_xxxzzz[k] - g_z_0_z_xxxzzz[k] * ab_z + g_z_0_z_xxxzzzz[k];

                g_z_0_zz_xxyyyy[k] = -g_z_xxyyyy[k] - g_z_0_z_xxyyyy[k] * ab_z + g_z_0_z_xxyyyyz[k];

                g_z_0_zz_xxyyyz[k] = -g_z_xxyyyz[k] - g_z_0_z_xxyyyz[k] * ab_z + g_z_0_z_xxyyyzz[k];

                g_z_0_zz_xxyyzz[k] = -g_z_xxyyzz[k] - g_z_0_z_xxyyzz[k] * ab_z + g_z_0_z_xxyyzzz[k];

                g_z_0_zz_xxyzzz[k] = -g_z_xxyzzz[k] - g_z_0_z_xxyzzz[k] * ab_z + g_z_0_z_xxyzzzz[k];

                g_z_0_zz_xxzzzz[k] = -g_z_xxzzzz[k] - g_z_0_z_xxzzzz[k] * ab_z + g_z_0_z_xxzzzzz[k];

                g_z_0_zz_xyyyyy[k] = -g_z_xyyyyy[k] - g_z_0_z_xyyyyy[k] * ab_z + g_z_0_z_xyyyyyz[k];

                g_z_0_zz_xyyyyz[k] = -g_z_xyyyyz[k] - g_z_0_z_xyyyyz[k] * ab_z + g_z_0_z_xyyyyzz[k];

                g_z_0_zz_xyyyzz[k] = -g_z_xyyyzz[k] - g_z_0_z_xyyyzz[k] * ab_z + g_z_0_z_xyyyzzz[k];

                g_z_0_zz_xyyzzz[k] = -g_z_xyyzzz[k] - g_z_0_z_xyyzzz[k] * ab_z + g_z_0_z_xyyzzzz[k];

                g_z_0_zz_xyzzzz[k] = -g_z_xyzzzz[k] - g_z_0_z_xyzzzz[k] * ab_z + g_z_0_z_xyzzzzz[k];

                g_z_0_zz_xzzzzz[k] = -g_z_xzzzzz[k] - g_z_0_z_xzzzzz[k] * ab_z + g_z_0_z_xzzzzzz[k];

                g_z_0_zz_yyyyyy[k] = -g_z_yyyyyy[k] - g_z_0_z_yyyyyy[k] * ab_z + g_z_0_z_yyyyyyz[k];

                g_z_0_zz_yyyyyz[k] = -g_z_yyyyyz[k] - g_z_0_z_yyyyyz[k] * ab_z + g_z_0_z_yyyyyzz[k];

                g_z_0_zz_yyyyzz[k] = -g_z_yyyyzz[k] - g_z_0_z_yyyyzz[k] * ab_z + g_z_0_z_yyyyzzz[k];

                g_z_0_zz_yyyzzz[k] = -g_z_yyyzzz[k] - g_z_0_z_yyyzzz[k] * ab_z + g_z_0_z_yyyzzzz[k];

                g_z_0_zz_yyzzzz[k] = -g_z_yyzzzz[k] - g_z_0_z_yyzzzz[k] * ab_z + g_z_0_z_yyzzzzz[k];

                g_z_0_zz_yzzzzz[k] = -g_z_yzzzzz[k] - g_z_0_z_yzzzzz[k] * ab_z + g_z_0_z_yzzzzzz[k];

                g_z_0_zz_zzzzzz[k] = -g_z_zzzzzz[k] - g_z_0_z_zzzzzz[k] * ab_z + g_z_0_z_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

