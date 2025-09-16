#include "ElectronRepulsionGeom1100ContrRecPIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_pixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_pixx,
                                            const size_t idx_geom_01_sixx,
                                            const size_t idx_geom_10_sixx,
                                            const size_t idx_geom_11_sixx,
                                            const size_t idx_geom_11_skxx,
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

            const auto si_geom_11_off = idx_geom_11_sixx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : SKSS

            const auto sk_geom_11_off = idx_geom_11_skxx + i * dcomps + j;

            auto g_x_x_0_xxxxxxx = cbuffer.data(sk_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_0_xxxxxxy = cbuffer.data(sk_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_0_xxxxxxz = cbuffer.data(sk_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_0_xxxxxyy = cbuffer.data(sk_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_0_xxxxxyz = cbuffer.data(sk_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_0_xxxxxzz = cbuffer.data(sk_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_0_xxxxyyy = cbuffer.data(sk_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_0_xxxxyyz = cbuffer.data(sk_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_0_xxxxyzz = cbuffer.data(sk_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_0_xxxxzzz = cbuffer.data(sk_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_0_xxxyyyy = cbuffer.data(sk_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_0_xxxyyyz = cbuffer.data(sk_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_0_xxxyyzz = cbuffer.data(sk_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_0_xxxyzzz = cbuffer.data(sk_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_0_xxxzzzz = cbuffer.data(sk_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_0_xxyyyyy = cbuffer.data(sk_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_0_xxyyyyz = cbuffer.data(sk_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_0_xxyyyzz = cbuffer.data(sk_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_0_xxyyzzz = cbuffer.data(sk_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_0_xxyzzzz = cbuffer.data(sk_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_0_xxzzzzz = cbuffer.data(sk_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_0_xyyyyyy = cbuffer.data(sk_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_0_xyyyyyz = cbuffer.data(sk_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_0_xyyyyzz = cbuffer.data(sk_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_0_xyyyzzz = cbuffer.data(sk_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_0_xyyzzzz = cbuffer.data(sk_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_0_xyzzzzz = cbuffer.data(sk_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_0_xzzzzzz = cbuffer.data(sk_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_0_yyyyyyy = cbuffer.data(sk_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_0_yyyyyyz = cbuffer.data(sk_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_0_yyyyyzz = cbuffer.data(sk_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_0_yyyyzzz = cbuffer.data(sk_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_0_yyyzzzz = cbuffer.data(sk_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_0_yyzzzzz = cbuffer.data(sk_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_0_yzzzzzz = cbuffer.data(sk_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_0_zzzzzzz = cbuffer.data(sk_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_y_0_xxxxxxx = cbuffer.data(sk_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_y_0_xxxxxxy = cbuffer.data(sk_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_y_0_xxxxxxz = cbuffer.data(sk_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_y_0_xxxxxyy = cbuffer.data(sk_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_y_0_xxxxxyz = cbuffer.data(sk_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_y_0_xxxxxzz = cbuffer.data(sk_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_y_0_xxxxyyy = cbuffer.data(sk_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_y_0_xxxxyyz = cbuffer.data(sk_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_y_0_xxxxyzz = cbuffer.data(sk_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_y_0_xxxxzzz = cbuffer.data(sk_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_y_0_xxxyyyy = cbuffer.data(sk_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_y_0_xxxyyyz = cbuffer.data(sk_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_y_0_xxxyyzz = cbuffer.data(sk_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_y_0_xxxyzzz = cbuffer.data(sk_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_y_0_xxxzzzz = cbuffer.data(sk_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_y_0_xxyyyyy = cbuffer.data(sk_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_y_0_xxyyyyz = cbuffer.data(sk_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_y_0_xxyyyzz = cbuffer.data(sk_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_y_0_xxyyzzz = cbuffer.data(sk_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_y_0_xxyzzzz = cbuffer.data(sk_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_y_0_xxzzzzz = cbuffer.data(sk_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_y_0_xyyyyyy = cbuffer.data(sk_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_y_0_xyyyyyz = cbuffer.data(sk_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_y_0_xyyyyzz = cbuffer.data(sk_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_y_0_xyyyzzz = cbuffer.data(sk_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_y_0_xyyzzzz = cbuffer.data(sk_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_y_0_xyzzzzz = cbuffer.data(sk_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_0_xzzzzzz = cbuffer.data(sk_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_0_yyyyyyy = cbuffer.data(sk_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_0_yyyyyyz = cbuffer.data(sk_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_0_yyyyyzz = cbuffer.data(sk_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_0_yyyyzzz = cbuffer.data(sk_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_0_yyyzzzz = cbuffer.data(sk_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_0_yyzzzzz = cbuffer.data(sk_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_0_yzzzzzz = cbuffer.data(sk_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_0_zzzzzzz = cbuffer.data(sk_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_z_0_xxxxxxx = cbuffer.data(sk_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_z_0_xxxxxxy = cbuffer.data(sk_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_z_0_xxxxxxz = cbuffer.data(sk_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_z_0_xxxxxyy = cbuffer.data(sk_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_z_0_xxxxxyz = cbuffer.data(sk_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_z_0_xxxxxzz = cbuffer.data(sk_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_z_0_xxxxyyy = cbuffer.data(sk_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_z_0_xxxxyyz = cbuffer.data(sk_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_z_0_xxxxyzz = cbuffer.data(sk_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_z_0_xxxxzzz = cbuffer.data(sk_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_z_0_xxxyyyy = cbuffer.data(sk_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_z_0_xxxyyyz = cbuffer.data(sk_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_z_0_xxxyyzz = cbuffer.data(sk_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_z_0_xxxyzzz = cbuffer.data(sk_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_z_0_xxxzzzz = cbuffer.data(sk_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_z_0_xxyyyyy = cbuffer.data(sk_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_z_0_xxyyyyz = cbuffer.data(sk_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_z_0_xxyyyzz = cbuffer.data(sk_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_z_0_xxyyzzz = cbuffer.data(sk_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_z_0_xxyzzzz = cbuffer.data(sk_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_z_0_xxzzzzz = cbuffer.data(sk_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_z_0_xyyyyyy = cbuffer.data(sk_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_z_0_xyyyyyz = cbuffer.data(sk_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_z_0_xyyyyzz = cbuffer.data(sk_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_z_0_xyyyzzz = cbuffer.data(sk_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_z_0_xyyzzzz = cbuffer.data(sk_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_z_0_xyzzzzz = cbuffer.data(sk_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_z_0_xzzzzzz = cbuffer.data(sk_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_z_0_yyyyyyy = cbuffer.data(sk_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_z_0_yyyyyyz = cbuffer.data(sk_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_z_0_yyyyyzz = cbuffer.data(sk_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_z_0_yyyyzzz = cbuffer.data(sk_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_z_0_yyyzzzz = cbuffer.data(sk_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_z_0_yyzzzzz = cbuffer.data(sk_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_z_0_yzzzzzz = cbuffer.data(sk_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_z_0_zzzzzzz = cbuffer.data(sk_geom_11_off + 107 * ccomps * dcomps);

            auto g_y_x_0_xxxxxxx = cbuffer.data(sk_geom_11_off + 108 * ccomps * dcomps);

            auto g_y_x_0_xxxxxxy = cbuffer.data(sk_geom_11_off + 109 * ccomps * dcomps);

            auto g_y_x_0_xxxxxxz = cbuffer.data(sk_geom_11_off + 110 * ccomps * dcomps);

            auto g_y_x_0_xxxxxyy = cbuffer.data(sk_geom_11_off + 111 * ccomps * dcomps);

            auto g_y_x_0_xxxxxyz = cbuffer.data(sk_geom_11_off + 112 * ccomps * dcomps);

            auto g_y_x_0_xxxxxzz = cbuffer.data(sk_geom_11_off + 113 * ccomps * dcomps);

            auto g_y_x_0_xxxxyyy = cbuffer.data(sk_geom_11_off + 114 * ccomps * dcomps);

            auto g_y_x_0_xxxxyyz = cbuffer.data(sk_geom_11_off + 115 * ccomps * dcomps);

            auto g_y_x_0_xxxxyzz = cbuffer.data(sk_geom_11_off + 116 * ccomps * dcomps);

            auto g_y_x_0_xxxxzzz = cbuffer.data(sk_geom_11_off + 117 * ccomps * dcomps);

            auto g_y_x_0_xxxyyyy = cbuffer.data(sk_geom_11_off + 118 * ccomps * dcomps);

            auto g_y_x_0_xxxyyyz = cbuffer.data(sk_geom_11_off + 119 * ccomps * dcomps);

            auto g_y_x_0_xxxyyzz = cbuffer.data(sk_geom_11_off + 120 * ccomps * dcomps);

            auto g_y_x_0_xxxyzzz = cbuffer.data(sk_geom_11_off + 121 * ccomps * dcomps);

            auto g_y_x_0_xxxzzzz = cbuffer.data(sk_geom_11_off + 122 * ccomps * dcomps);

            auto g_y_x_0_xxyyyyy = cbuffer.data(sk_geom_11_off + 123 * ccomps * dcomps);

            auto g_y_x_0_xxyyyyz = cbuffer.data(sk_geom_11_off + 124 * ccomps * dcomps);

            auto g_y_x_0_xxyyyzz = cbuffer.data(sk_geom_11_off + 125 * ccomps * dcomps);

            auto g_y_x_0_xxyyzzz = cbuffer.data(sk_geom_11_off + 126 * ccomps * dcomps);

            auto g_y_x_0_xxyzzzz = cbuffer.data(sk_geom_11_off + 127 * ccomps * dcomps);

            auto g_y_x_0_xxzzzzz = cbuffer.data(sk_geom_11_off + 128 * ccomps * dcomps);

            auto g_y_x_0_xyyyyyy = cbuffer.data(sk_geom_11_off + 129 * ccomps * dcomps);

            auto g_y_x_0_xyyyyyz = cbuffer.data(sk_geom_11_off + 130 * ccomps * dcomps);

            auto g_y_x_0_xyyyyzz = cbuffer.data(sk_geom_11_off + 131 * ccomps * dcomps);

            auto g_y_x_0_xyyyzzz = cbuffer.data(sk_geom_11_off + 132 * ccomps * dcomps);

            auto g_y_x_0_xyyzzzz = cbuffer.data(sk_geom_11_off + 133 * ccomps * dcomps);

            auto g_y_x_0_xyzzzzz = cbuffer.data(sk_geom_11_off + 134 * ccomps * dcomps);

            auto g_y_x_0_xzzzzzz = cbuffer.data(sk_geom_11_off + 135 * ccomps * dcomps);

            auto g_y_x_0_yyyyyyy = cbuffer.data(sk_geom_11_off + 136 * ccomps * dcomps);

            auto g_y_x_0_yyyyyyz = cbuffer.data(sk_geom_11_off + 137 * ccomps * dcomps);

            auto g_y_x_0_yyyyyzz = cbuffer.data(sk_geom_11_off + 138 * ccomps * dcomps);

            auto g_y_x_0_yyyyzzz = cbuffer.data(sk_geom_11_off + 139 * ccomps * dcomps);

            auto g_y_x_0_yyyzzzz = cbuffer.data(sk_geom_11_off + 140 * ccomps * dcomps);

            auto g_y_x_0_yyzzzzz = cbuffer.data(sk_geom_11_off + 141 * ccomps * dcomps);

            auto g_y_x_0_yzzzzzz = cbuffer.data(sk_geom_11_off + 142 * ccomps * dcomps);

            auto g_y_x_0_zzzzzzz = cbuffer.data(sk_geom_11_off + 143 * ccomps * dcomps);

            auto g_y_y_0_xxxxxxx = cbuffer.data(sk_geom_11_off + 144 * ccomps * dcomps);

            auto g_y_y_0_xxxxxxy = cbuffer.data(sk_geom_11_off + 145 * ccomps * dcomps);

            auto g_y_y_0_xxxxxxz = cbuffer.data(sk_geom_11_off + 146 * ccomps * dcomps);

            auto g_y_y_0_xxxxxyy = cbuffer.data(sk_geom_11_off + 147 * ccomps * dcomps);

            auto g_y_y_0_xxxxxyz = cbuffer.data(sk_geom_11_off + 148 * ccomps * dcomps);

            auto g_y_y_0_xxxxxzz = cbuffer.data(sk_geom_11_off + 149 * ccomps * dcomps);

            auto g_y_y_0_xxxxyyy = cbuffer.data(sk_geom_11_off + 150 * ccomps * dcomps);

            auto g_y_y_0_xxxxyyz = cbuffer.data(sk_geom_11_off + 151 * ccomps * dcomps);

            auto g_y_y_0_xxxxyzz = cbuffer.data(sk_geom_11_off + 152 * ccomps * dcomps);

            auto g_y_y_0_xxxxzzz = cbuffer.data(sk_geom_11_off + 153 * ccomps * dcomps);

            auto g_y_y_0_xxxyyyy = cbuffer.data(sk_geom_11_off + 154 * ccomps * dcomps);

            auto g_y_y_0_xxxyyyz = cbuffer.data(sk_geom_11_off + 155 * ccomps * dcomps);

            auto g_y_y_0_xxxyyzz = cbuffer.data(sk_geom_11_off + 156 * ccomps * dcomps);

            auto g_y_y_0_xxxyzzz = cbuffer.data(sk_geom_11_off + 157 * ccomps * dcomps);

            auto g_y_y_0_xxxzzzz = cbuffer.data(sk_geom_11_off + 158 * ccomps * dcomps);

            auto g_y_y_0_xxyyyyy = cbuffer.data(sk_geom_11_off + 159 * ccomps * dcomps);

            auto g_y_y_0_xxyyyyz = cbuffer.data(sk_geom_11_off + 160 * ccomps * dcomps);

            auto g_y_y_0_xxyyyzz = cbuffer.data(sk_geom_11_off + 161 * ccomps * dcomps);

            auto g_y_y_0_xxyyzzz = cbuffer.data(sk_geom_11_off + 162 * ccomps * dcomps);

            auto g_y_y_0_xxyzzzz = cbuffer.data(sk_geom_11_off + 163 * ccomps * dcomps);

            auto g_y_y_0_xxzzzzz = cbuffer.data(sk_geom_11_off + 164 * ccomps * dcomps);

            auto g_y_y_0_xyyyyyy = cbuffer.data(sk_geom_11_off + 165 * ccomps * dcomps);

            auto g_y_y_0_xyyyyyz = cbuffer.data(sk_geom_11_off + 166 * ccomps * dcomps);

            auto g_y_y_0_xyyyyzz = cbuffer.data(sk_geom_11_off + 167 * ccomps * dcomps);

            auto g_y_y_0_xyyyzzz = cbuffer.data(sk_geom_11_off + 168 * ccomps * dcomps);

            auto g_y_y_0_xyyzzzz = cbuffer.data(sk_geom_11_off + 169 * ccomps * dcomps);

            auto g_y_y_0_xyzzzzz = cbuffer.data(sk_geom_11_off + 170 * ccomps * dcomps);

            auto g_y_y_0_xzzzzzz = cbuffer.data(sk_geom_11_off + 171 * ccomps * dcomps);

            auto g_y_y_0_yyyyyyy = cbuffer.data(sk_geom_11_off + 172 * ccomps * dcomps);

            auto g_y_y_0_yyyyyyz = cbuffer.data(sk_geom_11_off + 173 * ccomps * dcomps);

            auto g_y_y_0_yyyyyzz = cbuffer.data(sk_geom_11_off + 174 * ccomps * dcomps);

            auto g_y_y_0_yyyyzzz = cbuffer.data(sk_geom_11_off + 175 * ccomps * dcomps);

            auto g_y_y_0_yyyzzzz = cbuffer.data(sk_geom_11_off + 176 * ccomps * dcomps);

            auto g_y_y_0_yyzzzzz = cbuffer.data(sk_geom_11_off + 177 * ccomps * dcomps);

            auto g_y_y_0_yzzzzzz = cbuffer.data(sk_geom_11_off + 178 * ccomps * dcomps);

            auto g_y_y_0_zzzzzzz = cbuffer.data(sk_geom_11_off + 179 * ccomps * dcomps);

            auto g_y_z_0_xxxxxxx = cbuffer.data(sk_geom_11_off + 180 * ccomps * dcomps);

            auto g_y_z_0_xxxxxxy = cbuffer.data(sk_geom_11_off + 181 * ccomps * dcomps);

            auto g_y_z_0_xxxxxxz = cbuffer.data(sk_geom_11_off + 182 * ccomps * dcomps);

            auto g_y_z_0_xxxxxyy = cbuffer.data(sk_geom_11_off + 183 * ccomps * dcomps);

            auto g_y_z_0_xxxxxyz = cbuffer.data(sk_geom_11_off + 184 * ccomps * dcomps);

            auto g_y_z_0_xxxxxzz = cbuffer.data(sk_geom_11_off + 185 * ccomps * dcomps);

            auto g_y_z_0_xxxxyyy = cbuffer.data(sk_geom_11_off + 186 * ccomps * dcomps);

            auto g_y_z_0_xxxxyyz = cbuffer.data(sk_geom_11_off + 187 * ccomps * dcomps);

            auto g_y_z_0_xxxxyzz = cbuffer.data(sk_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_z_0_xxxxzzz = cbuffer.data(sk_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_z_0_xxxyyyy = cbuffer.data(sk_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_z_0_xxxyyyz = cbuffer.data(sk_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_z_0_xxxyyzz = cbuffer.data(sk_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_z_0_xxxyzzz = cbuffer.data(sk_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_z_0_xxxzzzz = cbuffer.data(sk_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_z_0_xxyyyyy = cbuffer.data(sk_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_z_0_xxyyyyz = cbuffer.data(sk_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_z_0_xxyyyzz = cbuffer.data(sk_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_z_0_xxyyzzz = cbuffer.data(sk_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_z_0_xxyzzzz = cbuffer.data(sk_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_z_0_xxzzzzz = cbuffer.data(sk_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_z_0_xyyyyyy = cbuffer.data(sk_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_z_0_xyyyyyz = cbuffer.data(sk_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_z_0_xyyyyzz = cbuffer.data(sk_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_z_0_xyyyzzz = cbuffer.data(sk_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_z_0_xyyzzzz = cbuffer.data(sk_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_z_0_xyzzzzz = cbuffer.data(sk_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_z_0_xzzzzzz = cbuffer.data(sk_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_z_0_yyyyyyy = cbuffer.data(sk_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_z_0_yyyyyyz = cbuffer.data(sk_geom_11_off + 209 * ccomps * dcomps);

            auto g_y_z_0_yyyyyzz = cbuffer.data(sk_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_z_0_yyyyzzz = cbuffer.data(sk_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_z_0_yyyzzzz = cbuffer.data(sk_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_z_0_yyzzzzz = cbuffer.data(sk_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_z_0_yzzzzzz = cbuffer.data(sk_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_z_0_zzzzzzz = cbuffer.data(sk_geom_11_off + 215 * ccomps * dcomps);

            auto g_z_x_0_xxxxxxx = cbuffer.data(sk_geom_11_off + 216 * ccomps * dcomps);

            auto g_z_x_0_xxxxxxy = cbuffer.data(sk_geom_11_off + 217 * ccomps * dcomps);

            auto g_z_x_0_xxxxxxz = cbuffer.data(sk_geom_11_off + 218 * ccomps * dcomps);

            auto g_z_x_0_xxxxxyy = cbuffer.data(sk_geom_11_off + 219 * ccomps * dcomps);

            auto g_z_x_0_xxxxxyz = cbuffer.data(sk_geom_11_off + 220 * ccomps * dcomps);

            auto g_z_x_0_xxxxxzz = cbuffer.data(sk_geom_11_off + 221 * ccomps * dcomps);

            auto g_z_x_0_xxxxyyy = cbuffer.data(sk_geom_11_off + 222 * ccomps * dcomps);

            auto g_z_x_0_xxxxyyz = cbuffer.data(sk_geom_11_off + 223 * ccomps * dcomps);

            auto g_z_x_0_xxxxyzz = cbuffer.data(sk_geom_11_off + 224 * ccomps * dcomps);

            auto g_z_x_0_xxxxzzz = cbuffer.data(sk_geom_11_off + 225 * ccomps * dcomps);

            auto g_z_x_0_xxxyyyy = cbuffer.data(sk_geom_11_off + 226 * ccomps * dcomps);

            auto g_z_x_0_xxxyyyz = cbuffer.data(sk_geom_11_off + 227 * ccomps * dcomps);

            auto g_z_x_0_xxxyyzz = cbuffer.data(sk_geom_11_off + 228 * ccomps * dcomps);

            auto g_z_x_0_xxxyzzz = cbuffer.data(sk_geom_11_off + 229 * ccomps * dcomps);

            auto g_z_x_0_xxxzzzz = cbuffer.data(sk_geom_11_off + 230 * ccomps * dcomps);

            auto g_z_x_0_xxyyyyy = cbuffer.data(sk_geom_11_off + 231 * ccomps * dcomps);

            auto g_z_x_0_xxyyyyz = cbuffer.data(sk_geom_11_off + 232 * ccomps * dcomps);

            auto g_z_x_0_xxyyyzz = cbuffer.data(sk_geom_11_off + 233 * ccomps * dcomps);

            auto g_z_x_0_xxyyzzz = cbuffer.data(sk_geom_11_off + 234 * ccomps * dcomps);

            auto g_z_x_0_xxyzzzz = cbuffer.data(sk_geom_11_off + 235 * ccomps * dcomps);

            auto g_z_x_0_xxzzzzz = cbuffer.data(sk_geom_11_off + 236 * ccomps * dcomps);

            auto g_z_x_0_xyyyyyy = cbuffer.data(sk_geom_11_off + 237 * ccomps * dcomps);

            auto g_z_x_0_xyyyyyz = cbuffer.data(sk_geom_11_off + 238 * ccomps * dcomps);

            auto g_z_x_0_xyyyyzz = cbuffer.data(sk_geom_11_off + 239 * ccomps * dcomps);

            auto g_z_x_0_xyyyzzz = cbuffer.data(sk_geom_11_off + 240 * ccomps * dcomps);

            auto g_z_x_0_xyyzzzz = cbuffer.data(sk_geom_11_off + 241 * ccomps * dcomps);

            auto g_z_x_0_xyzzzzz = cbuffer.data(sk_geom_11_off + 242 * ccomps * dcomps);

            auto g_z_x_0_xzzzzzz = cbuffer.data(sk_geom_11_off + 243 * ccomps * dcomps);

            auto g_z_x_0_yyyyyyy = cbuffer.data(sk_geom_11_off + 244 * ccomps * dcomps);

            auto g_z_x_0_yyyyyyz = cbuffer.data(sk_geom_11_off + 245 * ccomps * dcomps);

            auto g_z_x_0_yyyyyzz = cbuffer.data(sk_geom_11_off + 246 * ccomps * dcomps);

            auto g_z_x_0_yyyyzzz = cbuffer.data(sk_geom_11_off + 247 * ccomps * dcomps);

            auto g_z_x_0_yyyzzzz = cbuffer.data(sk_geom_11_off + 248 * ccomps * dcomps);

            auto g_z_x_0_yyzzzzz = cbuffer.data(sk_geom_11_off + 249 * ccomps * dcomps);

            auto g_z_x_0_yzzzzzz = cbuffer.data(sk_geom_11_off + 250 * ccomps * dcomps);

            auto g_z_x_0_zzzzzzz = cbuffer.data(sk_geom_11_off + 251 * ccomps * dcomps);

            auto g_z_y_0_xxxxxxx = cbuffer.data(sk_geom_11_off + 252 * ccomps * dcomps);

            auto g_z_y_0_xxxxxxy = cbuffer.data(sk_geom_11_off + 253 * ccomps * dcomps);

            auto g_z_y_0_xxxxxxz = cbuffer.data(sk_geom_11_off + 254 * ccomps * dcomps);

            auto g_z_y_0_xxxxxyy = cbuffer.data(sk_geom_11_off + 255 * ccomps * dcomps);

            auto g_z_y_0_xxxxxyz = cbuffer.data(sk_geom_11_off + 256 * ccomps * dcomps);

            auto g_z_y_0_xxxxxzz = cbuffer.data(sk_geom_11_off + 257 * ccomps * dcomps);

            auto g_z_y_0_xxxxyyy = cbuffer.data(sk_geom_11_off + 258 * ccomps * dcomps);

            auto g_z_y_0_xxxxyyz = cbuffer.data(sk_geom_11_off + 259 * ccomps * dcomps);

            auto g_z_y_0_xxxxyzz = cbuffer.data(sk_geom_11_off + 260 * ccomps * dcomps);

            auto g_z_y_0_xxxxzzz = cbuffer.data(sk_geom_11_off + 261 * ccomps * dcomps);

            auto g_z_y_0_xxxyyyy = cbuffer.data(sk_geom_11_off + 262 * ccomps * dcomps);

            auto g_z_y_0_xxxyyyz = cbuffer.data(sk_geom_11_off + 263 * ccomps * dcomps);

            auto g_z_y_0_xxxyyzz = cbuffer.data(sk_geom_11_off + 264 * ccomps * dcomps);

            auto g_z_y_0_xxxyzzz = cbuffer.data(sk_geom_11_off + 265 * ccomps * dcomps);

            auto g_z_y_0_xxxzzzz = cbuffer.data(sk_geom_11_off + 266 * ccomps * dcomps);

            auto g_z_y_0_xxyyyyy = cbuffer.data(sk_geom_11_off + 267 * ccomps * dcomps);

            auto g_z_y_0_xxyyyyz = cbuffer.data(sk_geom_11_off + 268 * ccomps * dcomps);

            auto g_z_y_0_xxyyyzz = cbuffer.data(sk_geom_11_off + 269 * ccomps * dcomps);

            auto g_z_y_0_xxyyzzz = cbuffer.data(sk_geom_11_off + 270 * ccomps * dcomps);

            auto g_z_y_0_xxyzzzz = cbuffer.data(sk_geom_11_off + 271 * ccomps * dcomps);

            auto g_z_y_0_xxzzzzz = cbuffer.data(sk_geom_11_off + 272 * ccomps * dcomps);

            auto g_z_y_0_xyyyyyy = cbuffer.data(sk_geom_11_off + 273 * ccomps * dcomps);

            auto g_z_y_0_xyyyyyz = cbuffer.data(sk_geom_11_off + 274 * ccomps * dcomps);

            auto g_z_y_0_xyyyyzz = cbuffer.data(sk_geom_11_off + 275 * ccomps * dcomps);

            auto g_z_y_0_xyyyzzz = cbuffer.data(sk_geom_11_off + 276 * ccomps * dcomps);

            auto g_z_y_0_xyyzzzz = cbuffer.data(sk_geom_11_off + 277 * ccomps * dcomps);

            auto g_z_y_0_xyzzzzz = cbuffer.data(sk_geom_11_off + 278 * ccomps * dcomps);

            auto g_z_y_0_xzzzzzz = cbuffer.data(sk_geom_11_off + 279 * ccomps * dcomps);

            auto g_z_y_0_yyyyyyy = cbuffer.data(sk_geom_11_off + 280 * ccomps * dcomps);

            auto g_z_y_0_yyyyyyz = cbuffer.data(sk_geom_11_off + 281 * ccomps * dcomps);

            auto g_z_y_0_yyyyyzz = cbuffer.data(sk_geom_11_off + 282 * ccomps * dcomps);

            auto g_z_y_0_yyyyzzz = cbuffer.data(sk_geom_11_off + 283 * ccomps * dcomps);

            auto g_z_y_0_yyyzzzz = cbuffer.data(sk_geom_11_off + 284 * ccomps * dcomps);

            auto g_z_y_0_yyzzzzz = cbuffer.data(sk_geom_11_off + 285 * ccomps * dcomps);

            auto g_z_y_0_yzzzzzz = cbuffer.data(sk_geom_11_off + 286 * ccomps * dcomps);

            auto g_z_y_0_zzzzzzz = cbuffer.data(sk_geom_11_off + 287 * ccomps * dcomps);

            auto g_z_z_0_xxxxxxx = cbuffer.data(sk_geom_11_off + 288 * ccomps * dcomps);

            auto g_z_z_0_xxxxxxy = cbuffer.data(sk_geom_11_off + 289 * ccomps * dcomps);

            auto g_z_z_0_xxxxxxz = cbuffer.data(sk_geom_11_off + 290 * ccomps * dcomps);

            auto g_z_z_0_xxxxxyy = cbuffer.data(sk_geom_11_off + 291 * ccomps * dcomps);

            auto g_z_z_0_xxxxxyz = cbuffer.data(sk_geom_11_off + 292 * ccomps * dcomps);

            auto g_z_z_0_xxxxxzz = cbuffer.data(sk_geom_11_off + 293 * ccomps * dcomps);

            auto g_z_z_0_xxxxyyy = cbuffer.data(sk_geom_11_off + 294 * ccomps * dcomps);

            auto g_z_z_0_xxxxyyz = cbuffer.data(sk_geom_11_off + 295 * ccomps * dcomps);

            auto g_z_z_0_xxxxyzz = cbuffer.data(sk_geom_11_off + 296 * ccomps * dcomps);

            auto g_z_z_0_xxxxzzz = cbuffer.data(sk_geom_11_off + 297 * ccomps * dcomps);

            auto g_z_z_0_xxxyyyy = cbuffer.data(sk_geom_11_off + 298 * ccomps * dcomps);

            auto g_z_z_0_xxxyyyz = cbuffer.data(sk_geom_11_off + 299 * ccomps * dcomps);

            auto g_z_z_0_xxxyyzz = cbuffer.data(sk_geom_11_off + 300 * ccomps * dcomps);

            auto g_z_z_0_xxxyzzz = cbuffer.data(sk_geom_11_off + 301 * ccomps * dcomps);

            auto g_z_z_0_xxxzzzz = cbuffer.data(sk_geom_11_off + 302 * ccomps * dcomps);

            auto g_z_z_0_xxyyyyy = cbuffer.data(sk_geom_11_off + 303 * ccomps * dcomps);

            auto g_z_z_0_xxyyyyz = cbuffer.data(sk_geom_11_off + 304 * ccomps * dcomps);

            auto g_z_z_0_xxyyyzz = cbuffer.data(sk_geom_11_off + 305 * ccomps * dcomps);

            auto g_z_z_0_xxyyzzz = cbuffer.data(sk_geom_11_off + 306 * ccomps * dcomps);

            auto g_z_z_0_xxyzzzz = cbuffer.data(sk_geom_11_off + 307 * ccomps * dcomps);

            auto g_z_z_0_xxzzzzz = cbuffer.data(sk_geom_11_off + 308 * ccomps * dcomps);

            auto g_z_z_0_xyyyyyy = cbuffer.data(sk_geom_11_off + 309 * ccomps * dcomps);

            auto g_z_z_0_xyyyyyz = cbuffer.data(sk_geom_11_off + 310 * ccomps * dcomps);

            auto g_z_z_0_xyyyyzz = cbuffer.data(sk_geom_11_off + 311 * ccomps * dcomps);

            auto g_z_z_0_xyyyzzz = cbuffer.data(sk_geom_11_off + 312 * ccomps * dcomps);

            auto g_z_z_0_xyyzzzz = cbuffer.data(sk_geom_11_off + 313 * ccomps * dcomps);

            auto g_z_z_0_xyzzzzz = cbuffer.data(sk_geom_11_off + 314 * ccomps * dcomps);

            auto g_z_z_0_xzzzzzz = cbuffer.data(sk_geom_11_off + 315 * ccomps * dcomps);

            auto g_z_z_0_yyyyyyy = cbuffer.data(sk_geom_11_off + 316 * ccomps * dcomps);

            auto g_z_z_0_yyyyyyz = cbuffer.data(sk_geom_11_off + 317 * ccomps * dcomps);

            auto g_z_z_0_yyyyyzz = cbuffer.data(sk_geom_11_off + 318 * ccomps * dcomps);

            auto g_z_z_0_yyyyzzz = cbuffer.data(sk_geom_11_off + 319 * ccomps * dcomps);

            auto g_z_z_0_yyyzzzz = cbuffer.data(sk_geom_11_off + 320 * ccomps * dcomps);

            auto g_z_z_0_yyzzzzz = cbuffer.data(sk_geom_11_off + 321 * ccomps * dcomps);

            auto g_z_z_0_yzzzzzz = cbuffer.data(sk_geom_11_off + 322 * ccomps * dcomps);

            auto g_z_z_0_zzzzzzz = cbuffer.data(sk_geom_11_off + 323 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pixx

            const auto pi_geom_11_off = idx_geom_11_pixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_x_x_xxxxxx = cbuffer.data(pi_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xxxxxy = cbuffer.data(pi_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xxxxxz = cbuffer.data(pi_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_xxxxyy = cbuffer.data(pi_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_xxxxyz = cbuffer.data(pi_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_xxxxzz = cbuffer.data(pi_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_x_xxxyyy = cbuffer.data(pi_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_x_xxxyyz = cbuffer.data(pi_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_x_xxxyzz = cbuffer.data(pi_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_x_xxxzzz = cbuffer.data(pi_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_x_xxyyyy = cbuffer.data(pi_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_x_xxyyyz = cbuffer.data(pi_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_x_xxyyzz = cbuffer.data(pi_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_x_xxyzzz = cbuffer.data(pi_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_x_xxzzzz = cbuffer.data(pi_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_x_xyyyyy = cbuffer.data(pi_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_x_xyyyyz = cbuffer.data(pi_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_x_xyyyzz = cbuffer.data(pi_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_x_xyyzzz = cbuffer.data(pi_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_x_xyzzzz = cbuffer.data(pi_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_x_xzzzzz = cbuffer.data(pi_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_x_yyyyyy = cbuffer.data(pi_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_x_yyyyyz = cbuffer.data(pi_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_x_yyyyzz = cbuffer.data(pi_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_x_yyyzzz = cbuffer.data(pi_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_x_yyzzzz = cbuffer.data(pi_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_x_yzzzzz = cbuffer.data(pi_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_x_zzzzzz = cbuffer.data(pi_geom_11_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxx, g_0_x_0_xxxxxy, g_0_x_0_xxxxxz, g_0_x_0_xxxxyy, g_0_x_0_xxxxyz, g_0_x_0_xxxxzz, g_0_x_0_xxxyyy, g_0_x_0_xxxyyz, g_0_x_0_xxxyzz, g_0_x_0_xxxzzz, g_0_x_0_xxyyyy, g_0_x_0_xxyyyz, g_0_x_0_xxyyzz, g_0_x_0_xxyzzz, g_0_x_0_xxzzzz, g_0_x_0_xyyyyy, g_0_x_0_xyyyyz, g_0_x_0_xyyyzz, g_0_x_0_xyyzzz, g_0_x_0_xyzzzz, g_0_x_0_xzzzzz, g_0_x_0_yyyyyy, g_0_x_0_yyyyyz, g_0_x_0_yyyyzz, g_0_x_0_yyyzzz, g_0_x_0_yyzzzz, g_0_x_0_yzzzzz, g_0_x_0_zzzzzz, g_x_0_0_xxxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxxz, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxzz, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzzz, g_x_x_0_xxxxxx, g_x_x_0_xxxxxxx, g_x_x_0_xxxxxxy, g_x_x_0_xxxxxxz, g_x_x_0_xxxxxy, g_x_x_0_xxxxxyy, g_x_x_0_xxxxxyz, g_x_x_0_xxxxxz, g_x_x_0_xxxxxzz, g_x_x_0_xxxxyy, g_x_x_0_xxxxyyy, g_x_x_0_xxxxyyz, g_x_x_0_xxxxyz, g_x_x_0_xxxxyzz, g_x_x_0_xxxxzz, g_x_x_0_xxxxzzz, g_x_x_0_xxxyyy, g_x_x_0_xxxyyyy, g_x_x_0_xxxyyyz, g_x_x_0_xxxyyz, g_x_x_0_xxxyyzz, g_x_x_0_xxxyzz, g_x_x_0_xxxyzzz, g_x_x_0_xxxzzz, g_x_x_0_xxxzzzz, g_x_x_0_xxyyyy, g_x_x_0_xxyyyyy, g_x_x_0_xxyyyyz, g_x_x_0_xxyyyz, g_x_x_0_xxyyyzz, g_x_x_0_xxyyzz, g_x_x_0_xxyyzzz, g_x_x_0_xxyzzz, g_x_x_0_xxyzzzz, g_x_x_0_xxzzzz, g_x_x_0_xxzzzzz, g_x_x_0_xyyyyy, g_x_x_0_xyyyyyy, g_x_x_0_xyyyyyz, g_x_x_0_xyyyyz, g_x_x_0_xyyyyzz, g_x_x_0_xyyyzz, g_x_x_0_xyyyzzz, g_x_x_0_xyyzzz, g_x_x_0_xyyzzzz, g_x_x_0_xyzzzz, g_x_x_0_xyzzzzz, g_x_x_0_xzzzzz, g_x_x_0_xzzzzzz, g_x_x_0_yyyyyy, g_x_x_0_yyyyyz, g_x_x_0_yyyyzz, g_x_x_0_yyyzzz, g_x_x_0_yyzzzz, g_x_x_0_yzzzzz, g_x_x_0_zzzzzz, g_x_x_x_xxxxxx, g_x_x_x_xxxxxy, g_x_x_x_xxxxxz, g_x_x_x_xxxxyy, g_x_x_x_xxxxyz, g_x_x_x_xxxxzz, g_x_x_x_xxxyyy, g_x_x_x_xxxyyz, g_x_x_x_xxxyzz, g_x_x_x_xxxzzz, g_x_x_x_xxyyyy, g_x_x_x_xxyyyz, g_x_x_x_xxyyzz, g_x_x_x_xxyzzz, g_x_x_x_xxzzzz, g_x_x_x_xyyyyy, g_x_x_x_xyyyyz, g_x_x_x_xyyyzz, g_x_x_x_xyyzzz, g_x_x_x_xyzzzz, g_x_x_x_xzzzzz, g_x_x_x_yyyyyy, g_x_x_x_yyyyyz, g_x_x_x_yyyyzz, g_x_x_x_yyyzzz, g_x_x_x_yyzzzz, g_x_x_x_yzzzzz, g_x_x_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_x_xxxxxx[k] = -g_0_x_0_xxxxxx[k] + g_x_0_0_xxxxxx[k] - g_x_x_0_xxxxxx[k] * ab_x + g_x_x_0_xxxxxxx[k];

                g_x_x_x_xxxxxy[k] = -g_0_x_0_xxxxxy[k] + g_x_0_0_xxxxxy[k] - g_x_x_0_xxxxxy[k] * ab_x + g_x_x_0_xxxxxxy[k];

                g_x_x_x_xxxxxz[k] = -g_0_x_0_xxxxxz[k] + g_x_0_0_xxxxxz[k] - g_x_x_0_xxxxxz[k] * ab_x + g_x_x_0_xxxxxxz[k];

                g_x_x_x_xxxxyy[k] = -g_0_x_0_xxxxyy[k] + g_x_0_0_xxxxyy[k] - g_x_x_0_xxxxyy[k] * ab_x + g_x_x_0_xxxxxyy[k];

                g_x_x_x_xxxxyz[k] = -g_0_x_0_xxxxyz[k] + g_x_0_0_xxxxyz[k] - g_x_x_0_xxxxyz[k] * ab_x + g_x_x_0_xxxxxyz[k];

                g_x_x_x_xxxxzz[k] = -g_0_x_0_xxxxzz[k] + g_x_0_0_xxxxzz[k] - g_x_x_0_xxxxzz[k] * ab_x + g_x_x_0_xxxxxzz[k];

                g_x_x_x_xxxyyy[k] = -g_0_x_0_xxxyyy[k] + g_x_0_0_xxxyyy[k] - g_x_x_0_xxxyyy[k] * ab_x + g_x_x_0_xxxxyyy[k];

                g_x_x_x_xxxyyz[k] = -g_0_x_0_xxxyyz[k] + g_x_0_0_xxxyyz[k] - g_x_x_0_xxxyyz[k] * ab_x + g_x_x_0_xxxxyyz[k];

                g_x_x_x_xxxyzz[k] = -g_0_x_0_xxxyzz[k] + g_x_0_0_xxxyzz[k] - g_x_x_0_xxxyzz[k] * ab_x + g_x_x_0_xxxxyzz[k];

                g_x_x_x_xxxzzz[k] = -g_0_x_0_xxxzzz[k] + g_x_0_0_xxxzzz[k] - g_x_x_0_xxxzzz[k] * ab_x + g_x_x_0_xxxxzzz[k];

                g_x_x_x_xxyyyy[k] = -g_0_x_0_xxyyyy[k] + g_x_0_0_xxyyyy[k] - g_x_x_0_xxyyyy[k] * ab_x + g_x_x_0_xxxyyyy[k];

                g_x_x_x_xxyyyz[k] = -g_0_x_0_xxyyyz[k] + g_x_0_0_xxyyyz[k] - g_x_x_0_xxyyyz[k] * ab_x + g_x_x_0_xxxyyyz[k];

                g_x_x_x_xxyyzz[k] = -g_0_x_0_xxyyzz[k] + g_x_0_0_xxyyzz[k] - g_x_x_0_xxyyzz[k] * ab_x + g_x_x_0_xxxyyzz[k];

                g_x_x_x_xxyzzz[k] = -g_0_x_0_xxyzzz[k] + g_x_0_0_xxyzzz[k] - g_x_x_0_xxyzzz[k] * ab_x + g_x_x_0_xxxyzzz[k];

                g_x_x_x_xxzzzz[k] = -g_0_x_0_xxzzzz[k] + g_x_0_0_xxzzzz[k] - g_x_x_0_xxzzzz[k] * ab_x + g_x_x_0_xxxzzzz[k];

                g_x_x_x_xyyyyy[k] = -g_0_x_0_xyyyyy[k] + g_x_0_0_xyyyyy[k] - g_x_x_0_xyyyyy[k] * ab_x + g_x_x_0_xxyyyyy[k];

                g_x_x_x_xyyyyz[k] = -g_0_x_0_xyyyyz[k] + g_x_0_0_xyyyyz[k] - g_x_x_0_xyyyyz[k] * ab_x + g_x_x_0_xxyyyyz[k];

                g_x_x_x_xyyyzz[k] = -g_0_x_0_xyyyzz[k] + g_x_0_0_xyyyzz[k] - g_x_x_0_xyyyzz[k] * ab_x + g_x_x_0_xxyyyzz[k];

                g_x_x_x_xyyzzz[k] = -g_0_x_0_xyyzzz[k] + g_x_0_0_xyyzzz[k] - g_x_x_0_xyyzzz[k] * ab_x + g_x_x_0_xxyyzzz[k];

                g_x_x_x_xyzzzz[k] = -g_0_x_0_xyzzzz[k] + g_x_0_0_xyzzzz[k] - g_x_x_0_xyzzzz[k] * ab_x + g_x_x_0_xxyzzzz[k];

                g_x_x_x_xzzzzz[k] = -g_0_x_0_xzzzzz[k] + g_x_0_0_xzzzzz[k] - g_x_x_0_xzzzzz[k] * ab_x + g_x_x_0_xxzzzzz[k];

                g_x_x_x_yyyyyy[k] = -g_0_x_0_yyyyyy[k] + g_x_0_0_yyyyyy[k] - g_x_x_0_yyyyyy[k] * ab_x + g_x_x_0_xyyyyyy[k];

                g_x_x_x_yyyyyz[k] = -g_0_x_0_yyyyyz[k] + g_x_0_0_yyyyyz[k] - g_x_x_0_yyyyyz[k] * ab_x + g_x_x_0_xyyyyyz[k];

                g_x_x_x_yyyyzz[k] = -g_0_x_0_yyyyzz[k] + g_x_0_0_yyyyzz[k] - g_x_x_0_yyyyzz[k] * ab_x + g_x_x_0_xyyyyzz[k];

                g_x_x_x_yyyzzz[k] = -g_0_x_0_yyyzzz[k] + g_x_0_0_yyyzzz[k] - g_x_x_0_yyyzzz[k] * ab_x + g_x_x_0_xyyyzzz[k];

                g_x_x_x_yyzzzz[k] = -g_0_x_0_yyzzzz[k] + g_x_0_0_yyzzzz[k] - g_x_x_0_yyzzzz[k] * ab_x + g_x_x_0_xyyzzzz[k];

                g_x_x_x_yzzzzz[k] = -g_0_x_0_yzzzzz[k] + g_x_0_0_yzzzzz[k] - g_x_x_0_yzzzzz[k] * ab_x + g_x_x_0_xyzzzzz[k];

                g_x_x_x_zzzzzz[k] = -g_0_x_0_zzzzzz[k] + g_x_0_0_zzzzzz[k] - g_x_x_0_zzzzzz[k] * ab_x + g_x_x_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_x_y_xxxxxx = cbuffer.data(pi_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_y_xxxxxy = cbuffer.data(pi_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_y_xxxxxz = cbuffer.data(pi_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_y_xxxxyy = cbuffer.data(pi_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_y_xxxxyz = cbuffer.data(pi_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_y_xxxxzz = cbuffer.data(pi_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_y_xxxyyy = cbuffer.data(pi_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_y_xxxyyz = cbuffer.data(pi_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_y_xxxyzz = cbuffer.data(pi_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_y_xxxzzz = cbuffer.data(pi_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_y_xxyyyy = cbuffer.data(pi_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_y_xxyyyz = cbuffer.data(pi_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_y_xxyyzz = cbuffer.data(pi_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_y_xxyzzz = cbuffer.data(pi_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_y_xxzzzz = cbuffer.data(pi_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_y_xyyyyy = cbuffer.data(pi_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_y_xyyyyz = cbuffer.data(pi_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_y_xyyyzz = cbuffer.data(pi_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_y_xyyzzz = cbuffer.data(pi_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_y_xyzzzz = cbuffer.data(pi_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_y_xzzzzz = cbuffer.data(pi_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_y_yyyyyy = cbuffer.data(pi_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_y_yyyyyz = cbuffer.data(pi_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_y_yyyyzz = cbuffer.data(pi_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_y_yyyzzz = cbuffer.data(pi_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_y_yyzzzz = cbuffer.data(pi_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_y_yzzzzz = cbuffer.data(pi_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_y_zzzzzz = cbuffer.data(pi_geom_11_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_0_xxxxxx, g_x_x_0_xxxxxxy, g_x_x_0_xxxxxy, g_x_x_0_xxxxxyy, g_x_x_0_xxxxxyz, g_x_x_0_xxxxxz, g_x_x_0_xxxxyy, g_x_x_0_xxxxyyy, g_x_x_0_xxxxyyz, g_x_x_0_xxxxyz, g_x_x_0_xxxxyzz, g_x_x_0_xxxxzz, g_x_x_0_xxxyyy, g_x_x_0_xxxyyyy, g_x_x_0_xxxyyyz, g_x_x_0_xxxyyz, g_x_x_0_xxxyyzz, g_x_x_0_xxxyzz, g_x_x_0_xxxyzzz, g_x_x_0_xxxzzz, g_x_x_0_xxyyyy, g_x_x_0_xxyyyyy, g_x_x_0_xxyyyyz, g_x_x_0_xxyyyz, g_x_x_0_xxyyyzz, g_x_x_0_xxyyzz, g_x_x_0_xxyyzzz, g_x_x_0_xxyzzz, g_x_x_0_xxyzzzz, g_x_x_0_xxzzzz, g_x_x_0_xyyyyy, g_x_x_0_xyyyyyy, g_x_x_0_xyyyyyz, g_x_x_0_xyyyyz, g_x_x_0_xyyyyzz, g_x_x_0_xyyyzz, g_x_x_0_xyyyzzz, g_x_x_0_xyyzzz, g_x_x_0_xyyzzzz, g_x_x_0_xyzzzz, g_x_x_0_xyzzzzz, g_x_x_0_xzzzzz, g_x_x_0_yyyyyy, g_x_x_0_yyyyyyy, g_x_x_0_yyyyyyz, g_x_x_0_yyyyyz, g_x_x_0_yyyyyzz, g_x_x_0_yyyyzz, g_x_x_0_yyyyzzz, g_x_x_0_yyyzzz, g_x_x_0_yyyzzzz, g_x_x_0_yyzzzz, g_x_x_0_yyzzzzz, g_x_x_0_yzzzzz, g_x_x_0_yzzzzzz, g_x_x_0_zzzzzz, g_x_x_y_xxxxxx, g_x_x_y_xxxxxy, g_x_x_y_xxxxxz, g_x_x_y_xxxxyy, g_x_x_y_xxxxyz, g_x_x_y_xxxxzz, g_x_x_y_xxxyyy, g_x_x_y_xxxyyz, g_x_x_y_xxxyzz, g_x_x_y_xxxzzz, g_x_x_y_xxyyyy, g_x_x_y_xxyyyz, g_x_x_y_xxyyzz, g_x_x_y_xxyzzz, g_x_x_y_xxzzzz, g_x_x_y_xyyyyy, g_x_x_y_xyyyyz, g_x_x_y_xyyyzz, g_x_x_y_xyyzzz, g_x_x_y_xyzzzz, g_x_x_y_xzzzzz, g_x_x_y_yyyyyy, g_x_x_y_yyyyyz, g_x_x_y_yyyyzz, g_x_x_y_yyyzzz, g_x_x_y_yyzzzz, g_x_x_y_yzzzzz, g_x_x_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_y_xxxxxx[k] = -g_x_x_0_xxxxxx[k] * ab_y + g_x_x_0_xxxxxxy[k];

                g_x_x_y_xxxxxy[k] = -g_x_x_0_xxxxxy[k] * ab_y + g_x_x_0_xxxxxyy[k];

                g_x_x_y_xxxxxz[k] = -g_x_x_0_xxxxxz[k] * ab_y + g_x_x_0_xxxxxyz[k];

                g_x_x_y_xxxxyy[k] = -g_x_x_0_xxxxyy[k] * ab_y + g_x_x_0_xxxxyyy[k];

                g_x_x_y_xxxxyz[k] = -g_x_x_0_xxxxyz[k] * ab_y + g_x_x_0_xxxxyyz[k];

                g_x_x_y_xxxxzz[k] = -g_x_x_0_xxxxzz[k] * ab_y + g_x_x_0_xxxxyzz[k];

                g_x_x_y_xxxyyy[k] = -g_x_x_0_xxxyyy[k] * ab_y + g_x_x_0_xxxyyyy[k];

                g_x_x_y_xxxyyz[k] = -g_x_x_0_xxxyyz[k] * ab_y + g_x_x_0_xxxyyyz[k];

                g_x_x_y_xxxyzz[k] = -g_x_x_0_xxxyzz[k] * ab_y + g_x_x_0_xxxyyzz[k];

                g_x_x_y_xxxzzz[k] = -g_x_x_0_xxxzzz[k] * ab_y + g_x_x_0_xxxyzzz[k];

                g_x_x_y_xxyyyy[k] = -g_x_x_0_xxyyyy[k] * ab_y + g_x_x_0_xxyyyyy[k];

                g_x_x_y_xxyyyz[k] = -g_x_x_0_xxyyyz[k] * ab_y + g_x_x_0_xxyyyyz[k];

                g_x_x_y_xxyyzz[k] = -g_x_x_0_xxyyzz[k] * ab_y + g_x_x_0_xxyyyzz[k];

                g_x_x_y_xxyzzz[k] = -g_x_x_0_xxyzzz[k] * ab_y + g_x_x_0_xxyyzzz[k];

                g_x_x_y_xxzzzz[k] = -g_x_x_0_xxzzzz[k] * ab_y + g_x_x_0_xxyzzzz[k];

                g_x_x_y_xyyyyy[k] = -g_x_x_0_xyyyyy[k] * ab_y + g_x_x_0_xyyyyyy[k];

                g_x_x_y_xyyyyz[k] = -g_x_x_0_xyyyyz[k] * ab_y + g_x_x_0_xyyyyyz[k];

                g_x_x_y_xyyyzz[k] = -g_x_x_0_xyyyzz[k] * ab_y + g_x_x_0_xyyyyzz[k];

                g_x_x_y_xyyzzz[k] = -g_x_x_0_xyyzzz[k] * ab_y + g_x_x_0_xyyyzzz[k];

                g_x_x_y_xyzzzz[k] = -g_x_x_0_xyzzzz[k] * ab_y + g_x_x_0_xyyzzzz[k];

                g_x_x_y_xzzzzz[k] = -g_x_x_0_xzzzzz[k] * ab_y + g_x_x_0_xyzzzzz[k];

                g_x_x_y_yyyyyy[k] = -g_x_x_0_yyyyyy[k] * ab_y + g_x_x_0_yyyyyyy[k];

                g_x_x_y_yyyyyz[k] = -g_x_x_0_yyyyyz[k] * ab_y + g_x_x_0_yyyyyyz[k];

                g_x_x_y_yyyyzz[k] = -g_x_x_0_yyyyzz[k] * ab_y + g_x_x_0_yyyyyzz[k];

                g_x_x_y_yyyzzz[k] = -g_x_x_0_yyyzzz[k] * ab_y + g_x_x_0_yyyyzzz[k];

                g_x_x_y_yyzzzz[k] = -g_x_x_0_yyzzzz[k] * ab_y + g_x_x_0_yyyzzzz[k];

                g_x_x_y_yzzzzz[k] = -g_x_x_0_yzzzzz[k] * ab_y + g_x_x_0_yyzzzzz[k];

                g_x_x_y_zzzzzz[k] = -g_x_x_0_zzzzzz[k] * ab_y + g_x_x_0_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_x_z_xxxxxx = cbuffer.data(pi_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_z_xxxxxy = cbuffer.data(pi_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_z_xxxxxz = cbuffer.data(pi_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_z_xxxxyy = cbuffer.data(pi_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_z_xxxxyz = cbuffer.data(pi_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_z_xxxxzz = cbuffer.data(pi_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_z_xxxyyy = cbuffer.data(pi_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_z_xxxyyz = cbuffer.data(pi_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_z_xxxyzz = cbuffer.data(pi_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_z_xxxzzz = cbuffer.data(pi_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_z_xxyyyy = cbuffer.data(pi_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_z_xxyyyz = cbuffer.data(pi_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_z_xxyyzz = cbuffer.data(pi_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_z_xxyzzz = cbuffer.data(pi_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_z_xxzzzz = cbuffer.data(pi_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_z_xyyyyy = cbuffer.data(pi_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_z_xyyyyz = cbuffer.data(pi_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_z_xyyyzz = cbuffer.data(pi_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_z_xyyzzz = cbuffer.data(pi_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_z_xyzzzz = cbuffer.data(pi_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_z_xzzzzz = cbuffer.data(pi_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_z_yyyyyy = cbuffer.data(pi_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_z_yyyyyz = cbuffer.data(pi_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_z_yyyyzz = cbuffer.data(pi_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_z_yyyzzz = cbuffer.data(pi_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_z_yyzzzz = cbuffer.data(pi_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_z_yzzzzz = cbuffer.data(pi_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_z_zzzzzz = cbuffer.data(pi_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_0_xxxxxx, g_x_x_0_xxxxxxz, g_x_x_0_xxxxxy, g_x_x_0_xxxxxyz, g_x_x_0_xxxxxz, g_x_x_0_xxxxxzz, g_x_x_0_xxxxyy, g_x_x_0_xxxxyyz, g_x_x_0_xxxxyz, g_x_x_0_xxxxyzz, g_x_x_0_xxxxzz, g_x_x_0_xxxxzzz, g_x_x_0_xxxyyy, g_x_x_0_xxxyyyz, g_x_x_0_xxxyyz, g_x_x_0_xxxyyzz, g_x_x_0_xxxyzz, g_x_x_0_xxxyzzz, g_x_x_0_xxxzzz, g_x_x_0_xxxzzzz, g_x_x_0_xxyyyy, g_x_x_0_xxyyyyz, g_x_x_0_xxyyyz, g_x_x_0_xxyyyzz, g_x_x_0_xxyyzz, g_x_x_0_xxyyzzz, g_x_x_0_xxyzzz, g_x_x_0_xxyzzzz, g_x_x_0_xxzzzz, g_x_x_0_xxzzzzz, g_x_x_0_xyyyyy, g_x_x_0_xyyyyyz, g_x_x_0_xyyyyz, g_x_x_0_xyyyyzz, g_x_x_0_xyyyzz, g_x_x_0_xyyyzzz, g_x_x_0_xyyzzz, g_x_x_0_xyyzzzz, g_x_x_0_xyzzzz, g_x_x_0_xyzzzzz, g_x_x_0_xzzzzz, g_x_x_0_xzzzzzz, g_x_x_0_yyyyyy, g_x_x_0_yyyyyyz, g_x_x_0_yyyyyz, g_x_x_0_yyyyyzz, g_x_x_0_yyyyzz, g_x_x_0_yyyyzzz, g_x_x_0_yyyzzz, g_x_x_0_yyyzzzz, g_x_x_0_yyzzzz, g_x_x_0_yyzzzzz, g_x_x_0_yzzzzz, g_x_x_0_yzzzzzz, g_x_x_0_zzzzzz, g_x_x_0_zzzzzzz, g_x_x_z_xxxxxx, g_x_x_z_xxxxxy, g_x_x_z_xxxxxz, g_x_x_z_xxxxyy, g_x_x_z_xxxxyz, g_x_x_z_xxxxzz, g_x_x_z_xxxyyy, g_x_x_z_xxxyyz, g_x_x_z_xxxyzz, g_x_x_z_xxxzzz, g_x_x_z_xxyyyy, g_x_x_z_xxyyyz, g_x_x_z_xxyyzz, g_x_x_z_xxyzzz, g_x_x_z_xxzzzz, g_x_x_z_xyyyyy, g_x_x_z_xyyyyz, g_x_x_z_xyyyzz, g_x_x_z_xyyzzz, g_x_x_z_xyzzzz, g_x_x_z_xzzzzz, g_x_x_z_yyyyyy, g_x_x_z_yyyyyz, g_x_x_z_yyyyzz, g_x_x_z_yyyzzz, g_x_x_z_yyzzzz, g_x_x_z_yzzzzz, g_x_x_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_z_xxxxxx[k] = -g_x_x_0_xxxxxx[k] * ab_z + g_x_x_0_xxxxxxz[k];

                g_x_x_z_xxxxxy[k] = -g_x_x_0_xxxxxy[k] * ab_z + g_x_x_0_xxxxxyz[k];

                g_x_x_z_xxxxxz[k] = -g_x_x_0_xxxxxz[k] * ab_z + g_x_x_0_xxxxxzz[k];

                g_x_x_z_xxxxyy[k] = -g_x_x_0_xxxxyy[k] * ab_z + g_x_x_0_xxxxyyz[k];

                g_x_x_z_xxxxyz[k] = -g_x_x_0_xxxxyz[k] * ab_z + g_x_x_0_xxxxyzz[k];

                g_x_x_z_xxxxzz[k] = -g_x_x_0_xxxxzz[k] * ab_z + g_x_x_0_xxxxzzz[k];

                g_x_x_z_xxxyyy[k] = -g_x_x_0_xxxyyy[k] * ab_z + g_x_x_0_xxxyyyz[k];

                g_x_x_z_xxxyyz[k] = -g_x_x_0_xxxyyz[k] * ab_z + g_x_x_0_xxxyyzz[k];

                g_x_x_z_xxxyzz[k] = -g_x_x_0_xxxyzz[k] * ab_z + g_x_x_0_xxxyzzz[k];

                g_x_x_z_xxxzzz[k] = -g_x_x_0_xxxzzz[k] * ab_z + g_x_x_0_xxxzzzz[k];

                g_x_x_z_xxyyyy[k] = -g_x_x_0_xxyyyy[k] * ab_z + g_x_x_0_xxyyyyz[k];

                g_x_x_z_xxyyyz[k] = -g_x_x_0_xxyyyz[k] * ab_z + g_x_x_0_xxyyyzz[k];

                g_x_x_z_xxyyzz[k] = -g_x_x_0_xxyyzz[k] * ab_z + g_x_x_0_xxyyzzz[k];

                g_x_x_z_xxyzzz[k] = -g_x_x_0_xxyzzz[k] * ab_z + g_x_x_0_xxyzzzz[k];

                g_x_x_z_xxzzzz[k] = -g_x_x_0_xxzzzz[k] * ab_z + g_x_x_0_xxzzzzz[k];

                g_x_x_z_xyyyyy[k] = -g_x_x_0_xyyyyy[k] * ab_z + g_x_x_0_xyyyyyz[k];

                g_x_x_z_xyyyyz[k] = -g_x_x_0_xyyyyz[k] * ab_z + g_x_x_0_xyyyyzz[k];

                g_x_x_z_xyyyzz[k] = -g_x_x_0_xyyyzz[k] * ab_z + g_x_x_0_xyyyzzz[k];

                g_x_x_z_xyyzzz[k] = -g_x_x_0_xyyzzz[k] * ab_z + g_x_x_0_xyyzzzz[k];

                g_x_x_z_xyzzzz[k] = -g_x_x_0_xyzzzz[k] * ab_z + g_x_x_0_xyzzzzz[k];

                g_x_x_z_xzzzzz[k] = -g_x_x_0_xzzzzz[k] * ab_z + g_x_x_0_xzzzzzz[k];

                g_x_x_z_yyyyyy[k] = -g_x_x_0_yyyyyy[k] * ab_z + g_x_x_0_yyyyyyz[k];

                g_x_x_z_yyyyyz[k] = -g_x_x_0_yyyyyz[k] * ab_z + g_x_x_0_yyyyyzz[k];

                g_x_x_z_yyyyzz[k] = -g_x_x_0_yyyyzz[k] * ab_z + g_x_x_0_yyyyzzz[k];

                g_x_x_z_yyyzzz[k] = -g_x_x_0_yyyzzz[k] * ab_z + g_x_x_0_yyyzzzz[k];

                g_x_x_z_yyzzzz[k] = -g_x_x_0_yyzzzz[k] * ab_z + g_x_x_0_yyzzzzz[k];

                g_x_x_z_yzzzzz[k] = -g_x_x_0_yzzzzz[k] * ab_z + g_x_x_0_yzzzzzz[k];

                g_x_x_z_zzzzzz[k] = -g_x_x_0_zzzzzz[k] * ab_z + g_x_x_0_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_x_y_x_xxxxxx = cbuffer.data(pi_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_x_xxxxxy = cbuffer.data(pi_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_x_xxxxxz = cbuffer.data(pi_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_y_x_xxxxyy = cbuffer.data(pi_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_x_xxxxyz = cbuffer.data(pi_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_x_xxxxzz = cbuffer.data(pi_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_y_x_xxxyyy = cbuffer.data(pi_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_x_xxxyyz = cbuffer.data(pi_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_x_xxxyzz = cbuffer.data(pi_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_x_xxxzzz = cbuffer.data(pi_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_x_xxyyyy = cbuffer.data(pi_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_x_xxyyyz = cbuffer.data(pi_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_y_x_xxyyzz = cbuffer.data(pi_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_x_xxyzzz = cbuffer.data(pi_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_x_xxzzzz = cbuffer.data(pi_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_x_xyyyyy = cbuffer.data(pi_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_y_x_xyyyyz = cbuffer.data(pi_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_x_xyyyzz = cbuffer.data(pi_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_x_xyyzzz = cbuffer.data(pi_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_x_xyzzzz = cbuffer.data(pi_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_x_xzzzzz = cbuffer.data(pi_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_y_x_yyyyyy = cbuffer.data(pi_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_x_yyyyyz = cbuffer.data(pi_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_x_yyyyzz = cbuffer.data(pi_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_x_yyyzzz = cbuffer.data(pi_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_x_yyzzzz = cbuffer.data(pi_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_x_yzzzzz = cbuffer.data(pi_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_x_zzzzzz = cbuffer.data(pi_geom_11_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxxx, g_0_y_0_xxxxxy, g_0_y_0_xxxxxz, g_0_y_0_xxxxyy, g_0_y_0_xxxxyz, g_0_y_0_xxxxzz, g_0_y_0_xxxyyy, g_0_y_0_xxxyyz, g_0_y_0_xxxyzz, g_0_y_0_xxxzzz, g_0_y_0_xxyyyy, g_0_y_0_xxyyyz, g_0_y_0_xxyyzz, g_0_y_0_xxyzzz, g_0_y_0_xxzzzz, g_0_y_0_xyyyyy, g_0_y_0_xyyyyz, g_0_y_0_xyyyzz, g_0_y_0_xyyzzz, g_0_y_0_xyzzzz, g_0_y_0_xzzzzz, g_0_y_0_yyyyyy, g_0_y_0_yyyyyz, g_0_y_0_yyyyzz, g_0_y_0_yyyzzz, g_0_y_0_yyzzzz, g_0_y_0_yzzzzz, g_0_y_0_zzzzzz, g_x_y_0_xxxxxx, g_x_y_0_xxxxxxx, g_x_y_0_xxxxxxy, g_x_y_0_xxxxxxz, g_x_y_0_xxxxxy, g_x_y_0_xxxxxyy, g_x_y_0_xxxxxyz, g_x_y_0_xxxxxz, g_x_y_0_xxxxxzz, g_x_y_0_xxxxyy, g_x_y_0_xxxxyyy, g_x_y_0_xxxxyyz, g_x_y_0_xxxxyz, g_x_y_0_xxxxyzz, g_x_y_0_xxxxzz, g_x_y_0_xxxxzzz, g_x_y_0_xxxyyy, g_x_y_0_xxxyyyy, g_x_y_0_xxxyyyz, g_x_y_0_xxxyyz, g_x_y_0_xxxyyzz, g_x_y_0_xxxyzz, g_x_y_0_xxxyzzz, g_x_y_0_xxxzzz, g_x_y_0_xxxzzzz, g_x_y_0_xxyyyy, g_x_y_0_xxyyyyy, g_x_y_0_xxyyyyz, g_x_y_0_xxyyyz, g_x_y_0_xxyyyzz, g_x_y_0_xxyyzz, g_x_y_0_xxyyzzz, g_x_y_0_xxyzzz, g_x_y_0_xxyzzzz, g_x_y_0_xxzzzz, g_x_y_0_xxzzzzz, g_x_y_0_xyyyyy, g_x_y_0_xyyyyyy, g_x_y_0_xyyyyyz, g_x_y_0_xyyyyz, g_x_y_0_xyyyyzz, g_x_y_0_xyyyzz, g_x_y_0_xyyyzzz, g_x_y_0_xyyzzz, g_x_y_0_xyyzzzz, g_x_y_0_xyzzzz, g_x_y_0_xyzzzzz, g_x_y_0_xzzzzz, g_x_y_0_xzzzzzz, g_x_y_0_yyyyyy, g_x_y_0_yyyyyz, g_x_y_0_yyyyzz, g_x_y_0_yyyzzz, g_x_y_0_yyzzzz, g_x_y_0_yzzzzz, g_x_y_0_zzzzzz, g_x_y_x_xxxxxx, g_x_y_x_xxxxxy, g_x_y_x_xxxxxz, g_x_y_x_xxxxyy, g_x_y_x_xxxxyz, g_x_y_x_xxxxzz, g_x_y_x_xxxyyy, g_x_y_x_xxxyyz, g_x_y_x_xxxyzz, g_x_y_x_xxxzzz, g_x_y_x_xxyyyy, g_x_y_x_xxyyyz, g_x_y_x_xxyyzz, g_x_y_x_xxyzzz, g_x_y_x_xxzzzz, g_x_y_x_xyyyyy, g_x_y_x_xyyyyz, g_x_y_x_xyyyzz, g_x_y_x_xyyzzz, g_x_y_x_xyzzzz, g_x_y_x_xzzzzz, g_x_y_x_yyyyyy, g_x_y_x_yyyyyz, g_x_y_x_yyyyzz, g_x_y_x_yyyzzz, g_x_y_x_yyzzzz, g_x_y_x_yzzzzz, g_x_y_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_x_xxxxxx[k] = -g_0_y_0_xxxxxx[k] - g_x_y_0_xxxxxx[k] * ab_x + g_x_y_0_xxxxxxx[k];

                g_x_y_x_xxxxxy[k] = -g_0_y_0_xxxxxy[k] - g_x_y_0_xxxxxy[k] * ab_x + g_x_y_0_xxxxxxy[k];

                g_x_y_x_xxxxxz[k] = -g_0_y_0_xxxxxz[k] - g_x_y_0_xxxxxz[k] * ab_x + g_x_y_0_xxxxxxz[k];

                g_x_y_x_xxxxyy[k] = -g_0_y_0_xxxxyy[k] - g_x_y_0_xxxxyy[k] * ab_x + g_x_y_0_xxxxxyy[k];

                g_x_y_x_xxxxyz[k] = -g_0_y_0_xxxxyz[k] - g_x_y_0_xxxxyz[k] * ab_x + g_x_y_0_xxxxxyz[k];

                g_x_y_x_xxxxzz[k] = -g_0_y_0_xxxxzz[k] - g_x_y_0_xxxxzz[k] * ab_x + g_x_y_0_xxxxxzz[k];

                g_x_y_x_xxxyyy[k] = -g_0_y_0_xxxyyy[k] - g_x_y_0_xxxyyy[k] * ab_x + g_x_y_0_xxxxyyy[k];

                g_x_y_x_xxxyyz[k] = -g_0_y_0_xxxyyz[k] - g_x_y_0_xxxyyz[k] * ab_x + g_x_y_0_xxxxyyz[k];

                g_x_y_x_xxxyzz[k] = -g_0_y_0_xxxyzz[k] - g_x_y_0_xxxyzz[k] * ab_x + g_x_y_0_xxxxyzz[k];

                g_x_y_x_xxxzzz[k] = -g_0_y_0_xxxzzz[k] - g_x_y_0_xxxzzz[k] * ab_x + g_x_y_0_xxxxzzz[k];

                g_x_y_x_xxyyyy[k] = -g_0_y_0_xxyyyy[k] - g_x_y_0_xxyyyy[k] * ab_x + g_x_y_0_xxxyyyy[k];

                g_x_y_x_xxyyyz[k] = -g_0_y_0_xxyyyz[k] - g_x_y_0_xxyyyz[k] * ab_x + g_x_y_0_xxxyyyz[k];

                g_x_y_x_xxyyzz[k] = -g_0_y_0_xxyyzz[k] - g_x_y_0_xxyyzz[k] * ab_x + g_x_y_0_xxxyyzz[k];

                g_x_y_x_xxyzzz[k] = -g_0_y_0_xxyzzz[k] - g_x_y_0_xxyzzz[k] * ab_x + g_x_y_0_xxxyzzz[k];

                g_x_y_x_xxzzzz[k] = -g_0_y_0_xxzzzz[k] - g_x_y_0_xxzzzz[k] * ab_x + g_x_y_0_xxxzzzz[k];

                g_x_y_x_xyyyyy[k] = -g_0_y_0_xyyyyy[k] - g_x_y_0_xyyyyy[k] * ab_x + g_x_y_0_xxyyyyy[k];

                g_x_y_x_xyyyyz[k] = -g_0_y_0_xyyyyz[k] - g_x_y_0_xyyyyz[k] * ab_x + g_x_y_0_xxyyyyz[k];

                g_x_y_x_xyyyzz[k] = -g_0_y_0_xyyyzz[k] - g_x_y_0_xyyyzz[k] * ab_x + g_x_y_0_xxyyyzz[k];

                g_x_y_x_xyyzzz[k] = -g_0_y_0_xyyzzz[k] - g_x_y_0_xyyzzz[k] * ab_x + g_x_y_0_xxyyzzz[k];

                g_x_y_x_xyzzzz[k] = -g_0_y_0_xyzzzz[k] - g_x_y_0_xyzzzz[k] * ab_x + g_x_y_0_xxyzzzz[k];

                g_x_y_x_xzzzzz[k] = -g_0_y_0_xzzzzz[k] - g_x_y_0_xzzzzz[k] * ab_x + g_x_y_0_xxzzzzz[k];

                g_x_y_x_yyyyyy[k] = -g_0_y_0_yyyyyy[k] - g_x_y_0_yyyyyy[k] * ab_x + g_x_y_0_xyyyyyy[k];

                g_x_y_x_yyyyyz[k] = -g_0_y_0_yyyyyz[k] - g_x_y_0_yyyyyz[k] * ab_x + g_x_y_0_xyyyyyz[k];

                g_x_y_x_yyyyzz[k] = -g_0_y_0_yyyyzz[k] - g_x_y_0_yyyyzz[k] * ab_x + g_x_y_0_xyyyyzz[k];

                g_x_y_x_yyyzzz[k] = -g_0_y_0_yyyzzz[k] - g_x_y_0_yyyzzz[k] * ab_x + g_x_y_0_xyyyzzz[k];

                g_x_y_x_yyzzzz[k] = -g_0_y_0_yyzzzz[k] - g_x_y_0_yyzzzz[k] * ab_x + g_x_y_0_xyyzzzz[k];

                g_x_y_x_yzzzzz[k] = -g_0_y_0_yzzzzz[k] - g_x_y_0_yzzzzz[k] * ab_x + g_x_y_0_xyzzzzz[k];

                g_x_y_x_zzzzzz[k] = -g_0_y_0_zzzzzz[k] - g_x_y_0_zzzzzz[k] * ab_x + g_x_y_0_xzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_x_y_y_xxxxxx = cbuffer.data(pi_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_y_xxxxxy = cbuffer.data(pi_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_y_xxxxxz = cbuffer.data(pi_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_y_xxxxyy = cbuffer.data(pi_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_y_xxxxyz = cbuffer.data(pi_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_y_xxxxzz = cbuffer.data(pi_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_y_xxxyyy = cbuffer.data(pi_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_y_xxxyyz = cbuffer.data(pi_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_y_y_xxxyzz = cbuffer.data(pi_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_y_y_xxxzzz = cbuffer.data(pi_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_y_y_xxyyyy = cbuffer.data(pi_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_y_y_xxyyyz = cbuffer.data(pi_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_y_y_xxyyzz = cbuffer.data(pi_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_y_y_xxyzzz = cbuffer.data(pi_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_y_y_xxzzzz = cbuffer.data(pi_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_y_y_xyyyyy = cbuffer.data(pi_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_y_y_xyyyyz = cbuffer.data(pi_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_y_y_xyyyzz = cbuffer.data(pi_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_y_y_xyyzzz = cbuffer.data(pi_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_y_y_xyzzzz = cbuffer.data(pi_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_y_y_xzzzzz = cbuffer.data(pi_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_y_y_yyyyyy = cbuffer.data(pi_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_y_y_yyyyyz = cbuffer.data(pi_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_y_y_yyyyzz = cbuffer.data(pi_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_y_y_yyyzzz = cbuffer.data(pi_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_y_y_yyzzzz = cbuffer.data(pi_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_y_y_yzzzzz = cbuffer.data(pi_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_y_y_zzzzzz = cbuffer.data(pi_geom_11_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxxz, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxzz, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzzz, g_x_y_0_xxxxxx, g_x_y_0_xxxxxxy, g_x_y_0_xxxxxy, g_x_y_0_xxxxxyy, g_x_y_0_xxxxxyz, g_x_y_0_xxxxxz, g_x_y_0_xxxxyy, g_x_y_0_xxxxyyy, g_x_y_0_xxxxyyz, g_x_y_0_xxxxyz, g_x_y_0_xxxxyzz, g_x_y_0_xxxxzz, g_x_y_0_xxxyyy, g_x_y_0_xxxyyyy, g_x_y_0_xxxyyyz, g_x_y_0_xxxyyz, g_x_y_0_xxxyyzz, g_x_y_0_xxxyzz, g_x_y_0_xxxyzzz, g_x_y_0_xxxzzz, g_x_y_0_xxyyyy, g_x_y_0_xxyyyyy, g_x_y_0_xxyyyyz, g_x_y_0_xxyyyz, g_x_y_0_xxyyyzz, g_x_y_0_xxyyzz, g_x_y_0_xxyyzzz, g_x_y_0_xxyzzz, g_x_y_0_xxyzzzz, g_x_y_0_xxzzzz, g_x_y_0_xyyyyy, g_x_y_0_xyyyyyy, g_x_y_0_xyyyyyz, g_x_y_0_xyyyyz, g_x_y_0_xyyyyzz, g_x_y_0_xyyyzz, g_x_y_0_xyyyzzz, g_x_y_0_xyyzzz, g_x_y_0_xyyzzzz, g_x_y_0_xyzzzz, g_x_y_0_xyzzzzz, g_x_y_0_xzzzzz, g_x_y_0_yyyyyy, g_x_y_0_yyyyyyy, g_x_y_0_yyyyyyz, g_x_y_0_yyyyyz, g_x_y_0_yyyyyzz, g_x_y_0_yyyyzz, g_x_y_0_yyyyzzz, g_x_y_0_yyyzzz, g_x_y_0_yyyzzzz, g_x_y_0_yyzzzz, g_x_y_0_yyzzzzz, g_x_y_0_yzzzzz, g_x_y_0_yzzzzzz, g_x_y_0_zzzzzz, g_x_y_y_xxxxxx, g_x_y_y_xxxxxy, g_x_y_y_xxxxxz, g_x_y_y_xxxxyy, g_x_y_y_xxxxyz, g_x_y_y_xxxxzz, g_x_y_y_xxxyyy, g_x_y_y_xxxyyz, g_x_y_y_xxxyzz, g_x_y_y_xxxzzz, g_x_y_y_xxyyyy, g_x_y_y_xxyyyz, g_x_y_y_xxyyzz, g_x_y_y_xxyzzz, g_x_y_y_xxzzzz, g_x_y_y_xyyyyy, g_x_y_y_xyyyyz, g_x_y_y_xyyyzz, g_x_y_y_xyyzzz, g_x_y_y_xyzzzz, g_x_y_y_xzzzzz, g_x_y_y_yyyyyy, g_x_y_y_yyyyyz, g_x_y_y_yyyyzz, g_x_y_y_yyyzzz, g_x_y_y_yyzzzz, g_x_y_y_yzzzzz, g_x_y_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_y_xxxxxx[k] = g_x_0_0_xxxxxx[k] - g_x_y_0_xxxxxx[k] * ab_y + g_x_y_0_xxxxxxy[k];

                g_x_y_y_xxxxxy[k] = g_x_0_0_xxxxxy[k] - g_x_y_0_xxxxxy[k] * ab_y + g_x_y_0_xxxxxyy[k];

                g_x_y_y_xxxxxz[k] = g_x_0_0_xxxxxz[k] - g_x_y_0_xxxxxz[k] * ab_y + g_x_y_0_xxxxxyz[k];

                g_x_y_y_xxxxyy[k] = g_x_0_0_xxxxyy[k] - g_x_y_0_xxxxyy[k] * ab_y + g_x_y_0_xxxxyyy[k];

                g_x_y_y_xxxxyz[k] = g_x_0_0_xxxxyz[k] - g_x_y_0_xxxxyz[k] * ab_y + g_x_y_0_xxxxyyz[k];

                g_x_y_y_xxxxzz[k] = g_x_0_0_xxxxzz[k] - g_x_y_0_xxxxzz[k] * ab_y + g_x_y_0_xxxxyzz[k];

                g_x_y_y_xxxyyy[k] = g_x_0_0_xxxyyy[k] - g_x_y_0_xxxyyy[k] * ab_y + g_x_y_0_xxxyyyy[k];

                g_x_y_y_xxxyyz[k] = g_x_0_0_xxxyyz[k] - g_x_y_0_xxxyyz[k] * ab_y + g_x_y_0_xxxyyyz[k];

                g_x_y_y_xxxyzz[k] = g_x_0_0_xxxyzz[k] - g_x_y_0_xxxyzz[k] * ab_y + g_x_y_0_xxxyyzz[k];

                g_x_y_y_xxxzzz[k] = g_x_0_0_xxxzzz[k] - g_x_y_0_xxxzzz[k] * ab_y + g_x_y_0_xxxyzzz[k];

                g_x_y_y_xxyyyy[k] = g_x_0_0_xxyyyy[k] - g_x_y_0_xxyyyy[k] * ab_y + g_x_y_0_xxyyyyy[k];

                g_x_y_y_xxyyyz[k] = g_x_0_0_xxyyyz[k] - g_x_y_0_xxyyyz[k] * ab_y + g_x_y_0_xxyyyyz[k];

                g_x_y_y_xxyyzz[k] = g_x_0_0_xxyyzz[k] - g_x_y_0_xxyyzz[k] * ab_y + g_x_y_0_xxyyyzz[k];

                g_x_y_y_xxyzzz[k] = g_x_0_0_xxyzzz[k] - g_x_y_0_xxyzzz[k] * ab_y + g_x_y_0_xxyyzzz[k];

                g_x_y_y_xxzzzz[k] = g_x_0_0_xxzzzz[k] - g_x_y_0_xxzzzz[k] * ab_y + g_x_y_0_xxyzzzz[k];

                g_x_y_y_xyyyyy[k] = g_x_0_0_xyyyyy[k] - g_x_y_0_xyyyyy[k] * ab_y + g_x_y_0_xyyyyyy[k];

                g_x_y_y_xyyyyz[k] = g_x_0_0_xyyyyz[k] - g_x_y_0_xyyyyz[k] * ab_y + g_x_y_0_xyyyyyz[k];

                g_x_y_y_xyyyzz[k] = g_x_0_0_xyyyzz[k] - g_x_y_0_xyyyzz[k] * ab_y + g_x_y_0_xyyyyzz[k];

                g_x_y_y_xyyzzz[k] = g_x_0_0_xyyzzz[k] - g_x_y_0_xyyzzz[k] * ab_y + g_x_y_0_xyyyzzz[k];

                g_x_y_y_xyzzzz[k] = g_x_0_0_xyzzzz[k] - g_x_y_0_xyzzzz[k] * ab_y + g_x_y_0_xyyzzzz[k];

                g_x_y_y_xzzzzz[k] = g_x_0_0_xzzzzz[k] - g_x_y_0_xzzzzz[k] * ab_y + g_x_y_0_xyzzzzz[k];

                g_x_y_y_yyyyyy[k] = g_x_0_0_yyyyyy[k] - g_x_y_0_yyyyyy[k] * ab_y + g_x_y_0_yyyyyyy[k];

                g_x_y_y_yyyyyz[k] = g_x_0_0_yyyyyz[k] - g_x_y_0_yyyyyz[k] * ab_y + g_x_y_0_yyyyyyz[k];

                g_x_y_y_yyyyzz[k] = g_x_0_0_yyyyzz[k] - g_x_y_0_yyyyzz[k] * ab_y + g_x_y_0_yyyyyzz[k];

                g_x_y_y_yyyzzz[k] = g_x_0_0_yyyzzz[k] - g_x_y_0_yyyzzz[k] * ab_y + g_x_y_0_yyyyzzz[k];

                g_x_y_y_yyzzzz[k] = g_x_0_0_yyzzzz[k] - g_x_y_0_yyzzzz[k] * ab_y + g_x_y_0_yyyzzzz[k];

                g_x_y_y_yzzzzz[k] = g_x_0_0_yzzzzz[k] - g_x_y_0_yzzzzz[k] * ab_y + g_x_y_0_yyzzzzz[k];

                g_x_y_y_zzzzzz[k] = g_x_0_0_zzzzzz[k] - g_x_y_0_zzzzzz[k] * ab_y + g_x_y_0_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_x_y_z_xxxxxx = cbuffer.data(pi_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_y_z_xxxxxy = cbuffer.data(pi_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_y_z_xxxxxz = cbuffer.data(pi_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_y_z_xxxxyy = cbuffer.data(pi_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_y_z_xxxxyz = cbuffer.data(pi_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_y_z_xxxxzz = cbuffer.data(pi_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_y_z_xxxyyy = cbuffer.data(pi_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_y_z_xxxyyz = cbuffer.data(pi_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_y_z_xxxyzz = cbuffer.data(pi_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_y_z_xxxzzz = cbuffer.data(pi_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_y_z_xxyyyy = cbuffer.data(pi_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_z_xxyyyz = cbuffer.data(pi_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_z_xxyyzz = cbuffer.data(pi_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_z_xxyzzz = cbuffer.data(pi_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_z_xxzzzz = cbuffer.data(pi_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_z_xyyyyy = cbuffer.data(pi_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_y_z_xyyyyz = cbuffer.data(pi_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_z_xyyyzz = cbuffer.data(pi_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_z_xyyzzz = cbuffer.data(pi_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_z_xyzzzz = cbuffer.data(pi_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_z_xzzzzz = cbuffer.data(pi_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_z_yyyyyy = cbuffer.data(pi_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_y_z_yyyyyz = cbuffer.data(pi_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_z_yyyyzz = cbuffer.data(pi_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_z_yyyzzz = cbuffer.data(pi_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_y_z_yyzzzz = cbuffer.data(pi_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_z_yzzzzz = cbuffer.data(pi_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_z_zzzzzz = cbuffer.data(pi_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_0_xxxxxx, g_x_y_0_xxxxxxz, g_x_y_0_xxxxxy, g_x_y_0_xxxxxyz, g_x_y_0_xxxxxz, g_x_y_0_xxxxxzz, g_x_y_0_xxxxyy, g_x_y_0_xxxxyyz, g_x_y_0_xxxxyz, g_x_y_0_xxxxyzz, g_x_y_0_xxxxzz, g_x_y_0_xxxxzzz, g_x_y_0_xxxyyy, g_x_y_0_xxxyyyz, g_x_y_0_xxxyyz, g_x_y_0_xxxyyzz, g_x_y_0_xxxyzz, g_x_y_0_xxxyzzz, g_x_y_0_xxxzzz, g_x_y_0_xxxzzzz, g_x_y_0_xxyyyy, g_x_y_0_xxyyyyz, g_x_y_0_xxyyyz, g_x_y_0_xxyyyzz, g_x_y_0_xxyyzz, g_x_y_0_xxyyzzz, g_x_y_0_xxyzzz, g_x_y_0_xxyzzzz, g_x_y_0_xxzzzz, g_x_y_0_xxzzzzz, g_x_y_0_xyyyyy, g_x_y_0_xyyyyyz, g_x_y_0_xyyyyz, g_x_y_0_xyyyyzz, g_x_y_0_xyyyzz, g_x_y_0_xyyyzzz, g_x_y_0_xyyzzz, g_x_y_0_xyyzzzz, g_x_y_0_xyzzzz, g_x_y_0_xyzzzzz, g_x_y_0_xzzzzz, g_x_y_0_xzzzzzz, g_x_y_0_yyyyyy, g_x_y_0_yyyyyyz, g_x_y_0_yyyyyz, g_x_y_0_yyyyyzz, g_x_y_0_yyyyzz, g_x_y_0_yyyyzzz, g_x_y_0_yyyzzz, g_x_y_0_yyyzzzz, g_x_y_0_yyzzzz, g_x_y_0_yyzzzzz, g_x_y_0_yzzzzz, g_x_y_0_yzzzzzz, g_x_y_0_zzzzzz, g_x_y_0_zzzzzzz, g_x_y_z_xxxxxx, g_x_y_z_xxxxxy, g_x_y_z_xxxxxz, g_x_y_z_xxxxyy, g_x_y_z_xxxxyz, g_x_y_z_xxxxzz, g_x_y_z_xxxyyy, g_x_y_z_xxxyyz, g_x_y_z_xxxyzz, g_x_y_z_xxxzzz, g_x_y_z_xxyyyy, g_x_y_z_xxyyyz, g_x_y_z_xxyyzz, g_x_y_z_xxyzzz, g_x_y_z_xxzzzz, g_x_y_z_xyyyyy, g_x_y_z_xyyyyz, g_x_y_z_xyyyzz, g_x_y_z_xyyzzz, g_x_y_z_xyzzzz, g_x_y_z_xzzzzz, g_x_y_z_yyyyyy, g_x_y_z_yyyyyz, g_x_y_z_yyyyzz, g_x_y_z_yyyzzz, g_x_y_z_yyzzzz, g_x_y_z_yzzzzz, g_x_y_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_z_xxxxxx[k] = -g_x_y_0_xxxxxx[k] * ab_z + g_x_y_0_xxxxxxz[k];

                g_x_y_z_xxxxxy[k] = -g_x_y_0_xxxxxy[k] * ab_z + g_x_y_0_xxxxxyz[k];

                g_x_y_z_xxxxxz[k] = -g_x_y_0_xxxxxz[k] * ab_z + g_x_y_0_xxxxxzz[k];

                g_x_y_z_xxxxyy[k] = -g_x_y_0_xxxxyy[k] * ab_z + g_x_y_0_xxxxyyz[k];

                g_x_y_z_xxxxyz[k] = -g_x_y_0_xxxxyz[k] * ab_z + g_x_y_0_xxxxyzz[k];

                g_x_y_z_xxxxzz[k] = -g_x_y_0_xxxxzz[k] * ab_z + g_x_y_0_xxxxzzz[k];

                g_x_y_z_xxxyyy[k] = -g_x_y_0_xxxyyy[k] * ab_z + g_x_y_0_xxxyyyz[k];

                g_x_y_z_xxxyyz[k] = -g_x_y_0_xxxyyz[k] * ab_z + g_x_y_0_xxxyyzz[k];

                g_x_y_z_xxxyzz[k] = -g_x_y_0_xxxyzz[k] * ab_z + g_x_y_0_xxxyzzz[k];

                g_x_y_z_xxxzzz[k] = -g_x_y_0_xxxzzz[k] * ab_z + g_x_y_0_xxxzzzz[k];

                g_x_y_z_xxyyyy[k] = -g_x_y_0_xxyyyy[k] * ab_z + g_x_y_0_xxyyyyz[k];

                g_x_y_z_xxyyyz[k] = -g_x_y_0_xxyyyz[k] * ab_z + g_x_y_0_xxyyyzz[k];

                g_x_y_z_xxyyzz[k] = -g_x_y_0_xxyyzz[k] * ab_z + g_x_y_0_xxyyzzz[k];

                g_x_y_z_xxyzzz[k] = -g_x_y_0_xxyzzz[k] * ab_z + g_x_y_0_xxyzzzz[k];

                g_x_y_z_xxzzzz[k] = -g_x_y_0_xxzzzz[k] * ab_z + g_x_y_0_xxzzzzz[k];

                g_x_y_z_xyyyyy[k] = -g_x_y_0_xyyyyy[k] * ab_z + g_x_y_0_xyyyyyz[k];

                g_x_y_z_xyyyyz[k] = -g_x_y_0_xyyyyz[k] * ab_z + g_x_y_0_xyyyyzz[k];

                g_x_y_z_xyyyzz[k] = -g_x_y_0_xyyyzz[k] * ab_z + g_x_y_0_xyyyzzz[k];

                g_x_y_z_xyyzzz[k] = -g_x_y_0_xyyzzz[k] * ab_z + g_x_y_0_xyyzzzz[k];

                g_x_y_z_xyzzzz[k] = -g_x_y_0_xyzzzz[k] * ab_z + g_x_y_0_xyzzzzz[k];

                g_x_y_z_xzzzzz[k] = -g_x_y_0_xzzzzz[k] * ab_z + g_x_y_0_xzzzzzz[k];

                g_x_y_z_yyyyyy[k] = -g_x_y_0_yyyyyy[k] * ab_z + g_x_y_0_yyyyyyz[k];

                g_x_y_z_yyyyyz[k] = -g_x_y_0_yyyyyz[k] * ab_z + g_x_y_0_yyyyyzz[k];

                g_x_y_z_yyyyzz[k] = -g_x_y_0_yyyyzz[k] * ab_z + g_x_y_0_yyyyzzz[k];

                g_x_y_z_yyyzzz[k] = -g_x_y_0_yyyzzz[k] * ab_z + g_x_y_0_yyyzzzz[k];

                g_x_y_z_yyzzzz[k] = -g_x_y_0_yyzzzz[k] * ab_z + g_x_y_0_yyzzzzz[k];

                g_x_y_z_yzzzzz[k] = -g_x_y_0_yzzzzz[k] * ab_z + g_x_y_0_yzzzzzz[k];

                g_x_y_z_zzzzzz[k] = -g_x_y_0_zzzzzz[k] * ab_z + g_x_y_0_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_x_z_x_xxxxxx = cbuffer.data(pi_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_z_x_xxxxxy = cbuffer.data(pi_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_z_x_xxxxxz = cbuffer.data(pi_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_z_x_xxxxyy = cbuffer.data(pi_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_z_x_xxxxyz = cbuffer.data(pi_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_z_x_xxxxzz = cbuffer.data(pi_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_z_x_xxxyyy = cbuffer.data(pi_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_z_x_xxxyyz = cbuffer.data(pi_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_z_x_xxxyzz = cbuffer.data(pi_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_z_x_xxxzzz = cbuffer.data(pi_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_z_x_xxyyyy = cbuffer.data(pi_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_z_x_xxyyyz = cbuffer.data(pi_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_z_x_xxyyzz = cbuffer.data(pi_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_z_x_xxyzzz = cbuffer.data(pi_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_z_x_xxzzzz = cbuffer.data(pi_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_z_x_xyyyyy = cbuffer.data(pi_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_z_x_xyyyyz = cbuffer.data(pi_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_z_x_xyyyzz = cbuffer.data(pi_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_z_x_xyyzzz = cbuffer.data(pi_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_z_x_xyzzzz = cbuffer.data(pi_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_z_x_xzzzzz = cbuffer.data(pi_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_z_x_yyyyyy = cbuffer.data(pi_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_z_x_yyyyyz = cbuffer.data(pi_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_z_x_yyyyzz = cbuffer.data(pi_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_z_x_yyyzzz = cbuffer.data(pi_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_z_x_yyzzzz = cbuffer.data(pi_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_z_x_yzzzzz = cbuffer.data(pi_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_z_x_zzzzzz = cbuffer.data(pi_geom_11_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxxx, g_0_z_0_xxxxxy, g_0_z_0_xxxxxz, g_0_z_0_xxxxyy, g_0_z_0_xxxxyz, g_0_z_0_xxxxzz, g_0_z_0_xxxyyy, g_0_z_0_xxxyyz, g_0_z_0_xxxyzz, g_0_z_0_xxxzzz, g_0_z_0_xxyyyy, g_0_z_0_xxyyyz, g_0_z_0_xxyyzz, g_0_z_0_xxyzzz, g_0_z_0_xxzzzz, g_0_z_0_xyyyyy, g_0_z_0_xyyyyz, g_0_z_0_xyyyzz, g_0_z_0_xyyzzz, g_0_z_0_xyzzzz, g_0_z_0_xzzzzz, g_0_z_0_yyyyyy, g_0_z_0_yyyyyz, g_0_z_0_yyyyzz, g_0_z_0_yyyzzz, g_0_z_0_yyzzzz, g_0_z_0_yzzzzz, g_0_z_0_zzzzzz, g_x_z_0_xxxxxx, g_x_z_0_xxxxxxx, g_x_z_0_xxxxxxy, g_x_z_0_xxxxxxz, g_x_z_0_xxxxxy, g_x_z_0_xxxxxyy, g_x_z_0_xxxxxyz, g_x_z_0_xxxxxz, g_x_z_0_xxxxxzz, g_x_z_0_xxxxyy, g_x_z_0_xxxxyyy, g_x_z_0_xxxxyyz, g_x_z_0_xxxxyz, g_x_z_0_xxxxyzz, g_x_z_0_xxxxzz, g_x_z_0_xxxxzzz, g_x_z_0_xxxyyy, g_x_z_0_xxxyyyy, g_x_z_0_xxxyyyz, g_x_z_0_xxxyyz, g_x_z_0_xxxyyzz, g_x_z_0_xxxyzz, g_x_z_0_xxxyzzz, g_x_z_0_xxxzzz, g_x_z_0_xxxzzzz, g_x_z_0_xxyyyy, g_x_z_0_xxyyyyy, g_x_z_0_xxyyyyz, g_x_z_0_xxyyyz, g_x_z_0_xxyyyzz, g_x_z_0_xxyyzz, g_x_z_0_xxyyzzz, g_x_z_0_xxyzzz, g_x_z_0_xxyzzzz, g_x_z_0_xxzzzz, g_x_z_0_xxzzzzz, g_x_z_0_xyyyyy, g_x_z_0_xyyyyyy, g_x_z_0_xyyyyyz, g_x_z_0_xyyyyz, g_x_z_0_xyyyyzz, g_x_z_0_xyyyzz, g_x_z_0_xyyyzzz, g_x_z_0_xyyzzz, g_x_z_0_xyyzzzz, g_x_z_0_xyzzzz, g_x_z_0_xyzzzzz, g_x_z_0_xzzzzz, g_x_z_0_xzzzzzz, g_x_z_0_yyyyyy, g_x_z_0_yyyyyz, g_x_z_0_yyyyzz, g_x_z_0_yyyzzz, g_x_z_0_yyzzzz, g_x_z_0_yzzzzz, g_x_z_0_zzzzzz, g_x_z_x_xxxxxx, g_x_z_x_xxxxxy, g_x_z_x_xxxxxz, g_x_z_x_xxxxyy, g_x_z_x_xxxxyz, g_x_z_x_xxxxzz, g_x_z_x_xxxyyy, g_x_z_x_xxxyyz, g_x_z_x_xxxyzz, g_x_z_x_xxxzzz, g_x_z_x_xxyyyy, g_x_z_x_xxyyyz, g_x_z_x_xxyyzz, g_x_z_x_xxyzzz, g_x_z_x_xxzzzz, g_x_z_x_xyyyyy, g_x_z_x_xyyyyz, g_x_z_x_xyyyzz, g_x_z_x_xyyzzz, g_x_z_x_xyzzzz, g_x_z_x_xzzzzz, g_x_z_x_yyyyyy, g_x_z_x_yyyyyz, g_x_z_x_yyyyzz, g_x_z_x_yyyzzz, g_x_z_x_yyzzzz, g_x_z_x_yzzzzz, g_x_z_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_x_xxxxxx[k] = -g_0_z_0_xxxxxx[k] - g_x_z_0_xxxxxx[k] * ab_x + g_x_z_0_xxxxxxx[k];

                g_x_z_x_xxxxxy[k] = -g_0_z_0_xxxxxy[k] - g_x_z_0_xxxxxy[k] * ab_x + g_x_z_0_xxxxxxy[k];

                g_x_z_x_xxxxxz[k] = -g_0_z_0_xxxxxz[k] - g_x_z_0_xxxxxz[k] * ab_x + g_x_z_0_xxxxxxz[k];

                g_x_z_x_xxxxyy[k] = -g_0_z_0_xxxxyy[k] - g_x_z_0_xxxxyy[k] * ab_x + g_x_z_0_xxxxxyy[k];

                g_x_z_x_xxxxyz[k] = -g_0_z_0_xxxxyz[k] - g_x_z_0_xxxxyz[k] * ab_x + g_x_z_0_xxxxxyz[k];

                g_x_z_x_xxxxzz[k] = -g_0_z_0_xxxxzz[k] - g_x_z_0_xxxxzz[k] * ab_x + g_x_z_0_xxxxxzz[k];

                g_x_z_x_xxxyyy[k] = -g_0_z_0_xxxyyy[k] - g_x_z_0_xxxyyy[k] * ab_x + g_x_z_0_xxxxyyy[k];

                g_x_z_x_xxxyyz[k] = -g_0_z_0_xxxyyz[k] - g_x_z_0_xxxyyz[k] * ab_x + g_x_z_0_xxxxyyz[k];

                g_x_z_x_xxxyzz[k] = -g_0_z_0_xxxyzz[k] - g_x_z_0_xxxyzz[k] * ab_x + g_x_z_0_xxxxyzz[k];

                g_x_z_x_xxxzzz[k] = -g_0_z_0_xxxzzz[k] - g_x_z_0_xxxzzz[k] * ab_x + g_x_z_0_xxxxzzz[k];

                g_x_z_x_xxyyyy[k] = -g_0_z_0_xxyyyy[k] - g_x_z_0_xxyyyy[k] * ab_x + g_x_z_0_xxxyyyy[k];

                g_x_z_x_xxyyyz[k] = -g_0_z_0_xxyyyz[k] - g_x_z_0_xxyyyz[k] * ab_x + g_x_z_0_xxxyyyz[k];

                g_x_z_x_xxyyzz[k] = -g_0_z_0_xxyyzz[k] - g_x_z_0_xxyyzz[k] * ab_x + g_x_z_0_xxxyyzz[k];

                g_x_z_x_xxyzzz[k] = -g_0_z_0_xxyzzz[k] - g_x_z_0_xxyzzz[k] * ab_x + g_x_z_0_xxxyzzz[k];

                g_x_z_x_xxzzzz[k] = -g_0_z_0_xxzzzz[k] - g_x_z_0_xxzzzz[k] * ab_x + g_x_z_0_xxxzzzz[k];

                g_x_z_x_xyyyyy[k] = -g_0_z_0_xyyyyy[k] - g_x_z_0_xyyyyy[k] * ab_x + g_x_z_0_xxyyyyy[k];

                g_x_z_x_xyyyyz[k] = -g_0_z_0_xyyyyz[k] - g_x_z_0_xyyyyz[k] * ab_x + g_x_z_0_xxyyyyz[k];

                g_x_z_x_xyyyzz[k] = -g_0_z_0_xyyyzz[k] - g_x_z_0_xyyyzz[k] * ab_x + g_x_z_0_xxyyyzz[k];

                g_x_z_x_xyyzzz[k] = -g_0_z_0_xyyzzz[k] - g_x_z_0_xyyzzz[k] * ab_x + g_x_z_0_xxyyzzz[k];

                g_x_z_x_xyzzzz[k] = -g_0_z_0_xyzzzz[k] - g_x_z_0_xyzzzz[k] * ab_x + g_x_z_0_xxyzzzz[k];

                g_x_z_x_xzzzzz[k] = -g_0_z_0_xzzzzz[k] - g_x_z_0_xzzzzz[k] * ab_x + g_x_z_0_xxzzzzz[k];

                g_x_z_x_yyyyyy[k] = -g_0_z_0_yyyyyy[k] - g_x_z_0_yyyyyy[k] * ab_x + g_x_z_0_xyyyyyy[k];

                g_x_z_x_yyyyyz[k] = -g_0_z_0_yyyyyz[k] - g_x_z_0_yyyyyz[k] * ab_x + g_x_z_0_xyyyyyz[k];

                g_x_z_x_yyyyzz[k] = -g_0_z_0_yyyyzz[k] - g_x_z_0_yyyyzz[k] * ab_x + g_x_z_0_xyyyyzz[k];

                g_x_z_x_yyyzzz[k] = -g_0_z_0_yyyzzz[k] - g_x_z_0_yyyzzz[k] * ab_x + g_x_z_0_xyyyzzz[k];

                g_x_z_x_yyzzzz[k] = -g_0_z_0_yyzzzz[k] - g_x_z_0_yyzzzz[k] * ab_x + g_x_z_0_xyyzzzz[k];

                g_x_z_x_yzzzzz[k] = -g_0_z_0_yzzzzz[k] - g_x_z_0_yzzzzz[k] * ab_x + g_x_z_0_xyzzzzz[k];

                g_x_z_x_zzzzzz[k] = -g_0_z_0_zzzzzz[k] - g_x_z_0_zzzzzz[k] * ab_x + g_x_z_0_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_x_z_y_xxxxxx = cbuffer.data(pi_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_z_y_xxxxxy = cbuffer.data(pi_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_z_y_xxxxxz = cbuffer.data(pi_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_z_y_xxxxyy = cbuffer.data(pi_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_z_y_xxxxyz = cbuffer.data(pi_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_z_y_xxxxzz = cbuffer.data(pi_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_z_y_xxxyyy = cbuffer.data(pi_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_z_y_xxxyyz = cbuffer.data(pi_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_z_y_xxxyzz = cbuffer.data(pi_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_z_y_xxxzzz = cbuffer.data(pi_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_z_y_xxyyyy = cbuffer.data(pi_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_z_y_xxyyyz = cbuffer.data(pi_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_z_y_xxyyzz = cbuffer.data(pi_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_z_y_xxyzzz = cbuffer.data(pi_geom_11_off + 209 * ccomps * dcomps);

            auto g_x_z_y_xxzzzz = cbuffer.data(pi_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_z_y_xyyyyy = cbuffer.data(pi_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_z_y_xyyyyz = cbuffer.data(pi_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_z_y_xyyyzz = cbuffer.data(pi_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_z_y_xyyzzz = cbuffer.data(pi_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_z_y_xyzzzz = cbuffer.data(pi_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_z_y_xzzzzz = cbuffer.data(pi_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_z_y_yyyyyy = cbuffer.data(pi_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_z_y_yyyyyz = cbuffer.data(pi_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_z_y_yyyyzz = cbuffer.data(pi_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_z_y_yyyzzz = cbuffer.data(pi_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_z_y_yyzzzz = cbuffer.data(pi_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_z_y_yzzzzz = cbuffer.data(pi_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_z_y_zzzzzz = cbuffer.data(pi_geom_11_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_0_xxxxxx, g_x_z_0_xxxxxxy, g_x_z_0_xxxxxy, g_x_z_0_xxxxxyy, g_x_z_0_xxxxxyz, g_x_z_0_xxxxxz, g_x_z_0_xxxxyy, g_x_z_0_xxxxyyy, g_x_z_0_xxxxyyz, g_x_z_0_xxxxyz, g_x_z_0_xxxxyzz, g_x_z_0_xxxxzz, g_x_z_0_xxxyyy, g_x_z_0_xxxyyyy, g_x_z_0_xxxyyyz, g_x_z_0_xxxyyz, g_x_z_0_xxxyyzz, g_x_z_0_xxxyzz, g_x_z_0_xxxyzzz, g_x_z_0_xxxzzz, g_x_z_0_xxyyyy, g_x_z_0_xxyyyyy, g_x_z_0_xxyyyyz, g_x_z_0_xxyyyz, g_x_z_0_xxyyyzz, g_x_z_0_xxyyzz, g_x_z_0_xxyyzzz, g_x_z_0_xxyzzz, g_x_z_0_xxyzzzz, g_x_z_0_xxzzzz, g_x_z_0_xyyyyy, g_x_z_0_xyyyyyy, g_x_z_0_xyyyyyz, g_x_z_0_xyyyyz, g_x_z_0_xyyyyzz, g_x_z_0_xyyyzz, g_x_z_0_xyyyzzz, g_x_z_0_xyyzzz, g_x_z_0_xyyzzzz, g_x_z_0_xyzzzz, g_x_z_0_xyzzzzz, g_x_z_0_xzzzzz, g_x_z_0_yyyyyy, g_x_z_0_yyyyyyy, g_x_z_0_yyyyyyz, g_x_z_0_yyyyyz, g_x_z_0_yyyyyzz, g_x_z_0_yyyyzz, g_x_z_0_yyyyzzz, g_x_z_0_yyyzzz, g_x_z_0_yyyzzzz, g_x_z_0_yyzzzz, g_x_z_0_yyzzzzz, g_x_z_0_yzzzzz, g_x_z_0_yzzzzzz, g_x_z_0_zzzzzz, g_x_z_y_xxxxxx, g_x_z_y_xxxxxy, g_x_z_y_xxxxxz, g_x_z_y_xxxxyy, g_x_z_y_xxxxyz, g_x_z_y_xxxxzz, g_x_z_y_xxxyyy, g_x_z_y_xxxyyz, g_x_z_y_xxxyzz, g_x_z_y_xxxzzz, g_x_z_y_xxyyyy, g_x_z_y_xxyyyz, g_x_z_y_xxyyzz, g_x_z_y_xxyzzz, g_x_z_y_xxzzzz, g_x_z_y_xyyyyy, g_x_z_y_xyyyyz, g_x_z_y_xyyyzz, g_x_z_y_xyyzzz, g_x_z_y_xyzzzz, g_x_z_y_xzzzzz, g_x_z_y_yyyyyy, g_x_z_y_yyyyyz, g_x_z_y_yyyyzz, g_x_z_y_yyyzzz, g_x_z_y_yyzzzz, g_x_z_y_yzzzzz, g_x_z_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_y_xxxxxx[k] = -g_x_z_0_xxxxxx[k] * ab_y + g_x_z_0_xxxxxxy[k];

                g_x_z_y_xxxxxy[k] = -g_x_z_0_xxxxxy[k] * ab_y + g_x_z_0_xxxxxyy[k];

                g_x_z_y_xxxxxz[k] = -g_x_z_0_xxxxxz[k] * ab_y + g_x_z_0_xxxxxyz[k];

                g_x_z_y_xxxxyy[k] = -g_x_z_0_xxxxyy[k] * ab_y + g_x_z_0_xxxxyyy[k];

                g_x_z_y_xxxxyz[k] = -g_x_z_0_xxxxyz[k] * ab_y + g_x_z_0_xxxxyyz[k];

                g_x_z_y_xxxxzz[k] = -g_x_z_0_xxxxzz[k] * ab_y + g_x_z_0_xxxxyzz[k];

                g_x_z_y_xxxyyy[k] = -g_x_z_0_xxxyyy[k] * ab_y + g_x_z_0_xxxyyyy[k];

                g_x_z_y_xxxyyz[k] = -g_x_z_0_xxxyyz[k] * ab_y + g_x_z_0_xxxyyyz[k];

                g_x_z_y_xxxyzz[k] = -g_x_z_0_xxxyzz[k] * ab_y + g_x_z_0_xxxyyzz[k];

                g_x_z_y_xxxzzz[k] = -g_x_z_0_xxxzzz[k] * ab_y + g_x_z_0_xxxyzzz[k];

                g_x_z_y_xxyyyy[k] = -g_x_z_0_xxyyyy[k] * ab_y + g_x_z_0_xxyyyyy[k];

                g_x_z_y_xxyyyz[k] = -g_x_z_0_xxyyyz[k] * ab_y + g_x_z_0_xxyyyyz[k];

                g_x_z_y_xxyyzz[k] = -g_x_z_0_xxyyzz[k] * ab_y + g_x_z_0_xxyyyzz[k];

                g_x_z_y_xxyzzz[k] = -g_x_z_0_xxyzzz[k] * ab_y + g_x_z_0_xxyyzzz[k];

                g_x_z_y_xxzzzz[k] = -g_x_z_0_xxzzzz[k] * ab_y + g_x_z_0_xxyzzzz[k];

                g_x_z_y_xyyyyy[k] = -g_x_z_0_xyyyyy[k] * ab_y + g_x_z_0_xyyyyyy[k];

                g_x_z_y_xyyyyz[k] = -g_x_z_0_xyyyyz[k] * ab_y + g_x_z_0_xyyyyyz[k];

                g_x_z_y_xyyyzz[k] = -g_x_z_0_xyyyzz[k] * ab_y + g_x_z_0_xyyyyzz[k];

                g_x_z_y_xyyzzz[k] = -g_x_z_0_xyyzzz[k] * ab_y + g_x_z_0_xyyyzzz[k];

                g_x_z_y_xyzzzz[k] = -g_x_z_0_xyzzzz[k] * ab_y + g_x_z_0_xyyzzzz[k];

                g_x_z_y_xzzzzz[k] = -g_x_z_0_xzzzzz[k] * ab_y + g_x_z_0_xyzzzzz[k];

                g_x_z_y_yyyyyy[k] = -g_x_z_0_yyyyyy[k] * ab_y + g_x_z_0_yyyyyyy[k];

                g_x_z_y_yyyyyz[k] = -g_x_z_0_yyyyyz[k] * ab_y + g_x_z_0_yyyyyyz[k];

                g_x_z_y_yyyyzz[k] = -g_x_z_0_yyyyzz[k] * ab_y + g_x_z_0_yyyyyzz[k];

                g_x_z_y_yyyzzz[k] = -g_x_z_0_yyyzzz[k] * ab_y + g_x_z_0_yyyyzzz[k];

                g_x_z_y_yyzzzz[k] = -g_x_z_0_yyzzzz[k] * ab_y + g_x_z_0_yyyzzzz[k];

                g_x_z_y_yzzzzz[k] = -g_x_z_0_yzzzzz[k] * ab_y + g_x_z_0_yyzzzzz[k];

                g_x_z_y_zzzzzz[k] = -g_x_z_0_zzzzzz[k] * ab_y + g_x_z_0_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_x_z_z_xxxxxx = cbuffer.data(pi_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_z_z_xxxxxy = cbuffer.data(pi_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_z_z_xxxxxz = cbuffer.data(pi_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_z_z_xxxxyy = cbuffer.data(pi_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_z_z_xxxxyz = cbuffer.data(pi_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_z_z_xxxxzz = cbuffer.data(pi_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_z_z_xxxyyy = cbuffer.data(pi_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_z_z_xxxyyz = cbuffer.data(pi_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_z_z_xxxyzz = cbuffer.data(pi_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_z_z_xxxzzz = cbuffer.data(pi_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_z_z_xxyyyy = cbuffer.data(pi_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_z_z_xxyyyz = cbuffer.data(pi_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_z_z_xxyyzz = cbuffer.data(pi_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_z_z_xxyzzz = cbuffer.data(pi_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_z_z_xxzzzz = cbuffer.data(pi_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_z_z_xyyyyy = cbuffer.data(pi_geom_11_off + 239 * ccomps * dcomps);

            auto g_x_z_z_xyyyyz = cbuffer.data(pi_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_z_z_xyyyzz = cbuffer.data(pi_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_z_z_xyyzzz = cbuffer.data(pi_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_z_z_xyzzzz = cbuffer.data(pi_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_z_z_xzzzzz = cbuffer.data(pi_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_z_z_yyyyyy = cbuffer.data(pi_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_z_z_yyyyyz = cbuffer.data(pi_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_z_z_yyyyzz = cbuffer.data(pi_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_z_z_yyyzzz = cbuffer.data(pi_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_z_z_yyzzzz = cbuffer.data(pi_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_z_z_yzzzzz = cbuffer.data(pi_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_z_z_zzzzzz = cbuffer.data(pi_geom_11_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_0_xxxxxx, g_x_0_0_xxxxxy, g_x_0_0_xxxxxz, g_x_0_0_xxxxyy, g_x_0_0_xxxxyz, g_x_0_0_xxxxzz, g_x_0_0_xxxyyy, g_x_0_0_xxxyyz, g_x_0_0_xxxyzz, g_x_0_0_xxxzzz, g_x_0_0_xxyyyy, g_x_0_0_xxyyyz, g_x_0_0_xxyyzz, g_x_0_0_xxyzzz, g_x_0_0_xxzzzz, g_x_0_0_xyyyyy, g_x_0_0_xyyyyz, g_x_0_0_xyyyzz, g_x_0_0_xyyzzz, g_x_0_0_xyzzzz, g_x_0_0_xzzzzz, g_x_0_0_yyyyyy, g_x_0_0_yyyyyz, g_x_0_0_yyyyzz, g_x_0_0_yyyzzz, g_x_0_0_yyzzzz, g_x_0_0_yzzzzz, g_x_0_0_zzzzzz, g_x_z_0_xxxxxx, g_x_z_0_xxxxxxz, g_x_z_0_xxxxxy, g_x_z_0_xxxxxyz, g_x_z_0_xxxxxz, g_x_z_0_xxxxxzz, g_x_z_0_xxxxyy, g_x_z_0_xxxxyyz, g_x_z_0_xxxxyz, g_x_z_0_xxxxyzz, g_x_z_0_xxxxzz, g_x_z_0_xxxxzzz, g_x_z_0_xxxyyy, g_x_z_0_xxxyyyz, g_x_z_0_xxxyyz, g_x_z_0_xxxyyzz, g_x_z_0_xxxyzz, g_x_z_0_xxxyzzz, g_x_z_0_xxxzzz, g_x_z_0_xxxzzzz, g_x_z_0_xxyyyy, g_x_z_0_xxyyyyz, g_x_z_0_xxyyyz, g_x_z_0_xxyyyzz, g_x_z_0_xxyyzz, g_x_z_0_xxyyzzz, g_x_z_0_xxyzzz, g_x_z_0_xxyzzzz, g_x_z_0_xxzzzz, g_x_z_0_xxzzzzz, g_x_z_0_xyyyyy, g_x_z_0_xyyyyyz, g_x_z_0_xyyyyz, g_x_z_0_xyyyyzz, g_x_z_0_xyyyzz, g_x_z_0_xyyyzzz, g_x_z_0_xyyzzz, g_x_z_0_xyyzzzz, g_x_z_0_xyzzzz, g_x_z_0_xyzzzzz, g_x_z_0_xzzzzz, g_x_z_0_xzzzzzz, g_x_z_0_yyyyyy, g_x_z_0_yyyyyyz, g_x_z_0_yyyyyz, g_x_z_0_yyyyyzz, g_x_z_0_yyyyzz, g_x_z_0_yyyyzzz, g_x_z_0_yyyzzz, g_x_z_0_yyyzzzz, g_x_z_0_yyzzzz, g_x_z_0_yyzzzzz, g_x_z_0_yzzzzz, g_x_z_0_yzzzzzz, g_x_z_0_zzzzzz, g_x_z_0_zzzzzzz, g_x_z_z_xxxxxx, g_x_z_z_xxxxxy, g_x_z_z_xxxxxz, g_x_z_z_xxxxyy, g_x_z_z_xxxxyz, g_x_z_z_xxxxzz, g_x_z_z_xxxyyy, g_x_z_z_xxxyyz, g_x_z_z_xxxyzz, g_x_z_z_xxxzzz, g_x_z_z_xxyyyy, g_x_z_z_xxyyyz, g_x_z_z_xxyyzz, g_x_z_z_xxyzzz, g_x_z_z_xxzzzz, g_x_z_z_xyyyyy, g_x_z_z_xyyyyz, g_x_z_z_xyyyzz, g_x_z_z_xyyzzz, g_x_z_z_xyzzzz, g_x_z_z_xzzzzz, g_x_z_z_yyyyyy, g_x_z_z_yyyyyz, g_x_z_z_yyyyzz, g_x_z_z_yyyzzz, g_x_z_z_yyzzzz, g_x_z_z_yzzzzz, g_x_z_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_z_xxxxxx[k] = g_x_0_0_xxxxxx[k] - g_x_z_0_xxxxxx[k] * ab_z + g_x_z_0_xxxxxxz[k];

                g_x_z_z_xxxxxy[k] = g_x_0_0_xxxxxy[k] - g_x_z_0_xxxxxy[k] * ab_z + g_x_z_0_xxxxxyz[k];

                g_x_z_z_xxxxxz[k] = g_x_0_0_xxxxxz[k] - g_x_z_0_xxxxxz[k] * ab_z + g_x_z_0_xxxxxzz[k];

                g_x_z_z_xxxxyy[k] = g_x_0_0_xxxxyy[k] - g_x_z_0_xxxxyy[k] * ab_z + g_x_z_0_xxxxyyz[k];

                g_x_z_z_xxxxyz[k] = g_x_0_0_xxxxyz[k] - g_x_z_0_xxxxyz[k] * ab_z + g_x_z_0_xxxxyzz[k];

                g_x_z_z_xxxxzz[k] = g_x_0_0_xxxxzz[k] - g_x_z_0_xxxxzz[k] * ab_z + g_x_z_0_xxxxzzz[k];

                g_x_z_z_xxxyyy[k] = g_x_0_0_xxxyyy[k] - g_x_z_0_xxxyyy[k] * ab_z + g_x_z_0_xxxyyyz[k];

                g_x_z_z_xxxyyz[k] = g_x_0_0_xxxyyz[k] - g_x_z_0_xxxyyz[k] * ab_z + g_x_z_0_xxxyyzz[k];

                g_x_z_z_xxxyzz[k] = g_x_0_0_xxxyzz[k] - g_x_z_0_xxxyzz[k] * ab_z + g_x_z_0_xxxyzzz[k];

                g_x_z_z_xxxzzz[k] = g_x_0_0_xxxzzz[k] - g_x_z_0_xxxzzz[k] * ab_z + g_x_z_0_xxxzzzz[k];

                g_x_z_z_xxyyyy[k] = g_x_0_0_xxyyyy[k] - g_x_z_0_xxyyyy[k] * ab_z + g_x_z_0_xxyyyyz[k];

                g_x_z_z_xxyyyz[k] = g_x_0_0_xxyyyz[k] - g_x_z_0_xxyyyz[k] * ab_z + g_x_z_0_xxyyyzz[k];

                g_x_z_z_xxyyzz[k] = g_x_0_0_xxyyzz[k] - g_x_z_0_xxyyzz[k] * ab_z + g_x_z_0_xxyyzzz[k];

                g_x_z_z_xxyzzz[k] = g_x_0_0_xxyzzz[k] - g_x_z_0_xxyzzz[k] * ab_z + g_x_z_0_xxyzzzz[k];

                g_x_z_z_xxzzzz[k] = g_x_0_0_xxzzzz[k] - g_x_z_0_xxzzzz[k] * ab_z + g_x_z_0_xxzzzzz[k];

                g_x_z_z_xyyyyy[k] = g_x_0_0_xyyyyy[k] - g_x_z_0_xyyyyy[k] * ab_z + g_x_z_0_xyyyyyz[k];

                g_x_z_z_xyyyyz[k] = g_x_0_0_xyyyyz[k] - g_x_z_0_xyyyyz[k] * ab_z + g_x_z_0_xyyyyzz[k];

                g_x_z_z_xyyyzz[k] = g_x_0_0_xyyyzz[k] - g_x_z_0_xyyyzz[k] * ab_z + g_x_z_0_xyyyzzz[k];

                g_x_z_z_xyyzzz[k] = g_x_0_0_xyyzzz[k] - g_x_z_0_xyyzzz[k] * ab_z + g_x_z_0_xyyzzzz[k];

                g_x_z_z_xyzzzz[k] = g_x_0_0_xyzzzz[k] - g_x_z_0_xyzzzz[k] * ab_z + g_x_z_0_xyzzzzz[k];

                g_x_z_z_xzzzzz[k] = g_x_0_0_xzzzzz[k] - g_x_z_0_xzzzzz[k] * ab_z + g_x_z_0_xzzzzzz[k];

                g_x_z_z_yyyyyy[k] = g_x_0_0_yyyyyy[k] - g_x_z_0_yyyyyy[k] * ab_z + g_x_z_0_yyyyyyz[k];

                g_x_z_z_yyyyyz[k] = g_x_0_0_yyyyyz[k] - g_x_z_0_yyyyyz[k] * ab_z + g_x_z_0_yyyyyzz[k];

                g_x_z_z_yyyyzz[k] = g_x_0_0_yyyyzz[k] - g_x_z_0_yyyyzz[k] * ab_z + g_x_z_0_yyyyzzz[k];

                g_x_z_z_yyyzzz[k] = g_x_0_0_yyyzzz[k] - g_x_z_0_yyyzzz[k] * ab_z + g_x_z_0_yyyzzzz[k];

                g_x_z_z_yyzzzz[k] = g_x_0_0_yyzzzz[k] - g_x_z_0_yyzzzz[k] * ab_z + g_x_z_0_yyzzzzz[k];

                g_x_z_z_yzzzzz[k] = g_x_0_0_yzzzzz[k] - g_x_z_0_yzzzzz[k] * ab_z + g_x_z_0_yzzzzzz[k];

                g_x_z_z_zzzzzz[k] = g_x_0_0_zzzzzz[k] - g_x_z_0_zzzzzz[k] * ab_z + g_x_z_0_zzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_y_x_x_xxxxxx = cbuffer.data(pi_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_x_x_xxxxxy = cbuffer.data(pi_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_x_x_xxxxxz = cbuffer.data(pi_geom_11_off + 254 * ccomps * dcomps);

            auto g_y_x_x_xxxxyy = cbuffer.data(pi_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_x_x_xxxxyz = cbuffer.data(pi_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_x_x_xxxxzz = cbuffer.data(pi_geom_11_off + 257 * ccomps * dcomps);

            auto g_y_x_x_xxxyyy = cbuffer.data(pi_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_x_x_xxxyyz = cbuffer.data(pi_geom_11_off + 259 * ccomps * dcomps);

            auto g_y_x_x_xxxyzz = cbuffer.data(pi_geom_11_off + 260 * ccomps * dcomps);

            auto g_y_x_x_xxxzzz = cbuffer.data(pi_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_x_x_xxyyyy = cbuffer.data(pi_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_x_x_xxyyyz = cbuffer.data(pi_geom_11_off + 263 * ccomps * dcomps);

            auto g_y_x_x_xxyyzz = cbuffer.data(pi_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_x_x_xxyzzz = cbuffer.data(pi_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_x_x_xxzzzz = cbuffer.data(pi_geom_11_off + 266 * ccomps * dcomps);

            auto g_y_x_x_xyyyyy = cbuffer.data(pi_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_x_x_xyyyyz = cbuffer.data(pi_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_x_x_xyyyzz = cbuffer.data(pi_geom_11_off + 269 * ccomps * dcomps);

            auto g_y_x_x_xyyzzz = cbuffer.data(pi_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_x_x_xyzzzz = cbuffer.data(pi_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_x_x_xzzzzz = cbuffer.data(pi_geom_11_off + 272 * ccomps * dcomps);

            auto g_y_x_x_yyyyyy = cbuffer.data(pi_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_x_x_yyyyyz = cbuffer.data(pi_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_x_x_yyyyzz = cbuffer.data(pi_geom_11_off + 275 * ccomps * dcomps);

            auto g_y_x_x_yyyzzz = cbuffer.data(pi_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_x_x_yyzzzz = cbuffer.data(pi_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_x_x_yzzzzz = cbuffer.data(pi_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_x_x_zzzzzz = cbuffer.data(pi_geom_11_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxxz, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxzz, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzzz, g_y_x_0_xxxxxx, g_y_x_0_xxxxxxx, g_y_x_0_xxxxxxy, g_y_x_0_xxxxxxz, g_y_x_0_xxxxxy, g_y_x_0_xxxxxyy, g_y_x_0_xxxxxyz, g_y_x_0_xxxxxz, g_y_x_0_xxxxxzz, g_y_x_0_xxxxyy, g_y_x_0_xxxxyyy, g_y_x_0_xxxxyyz, g_y_x_0_xxxxyz, g_y_x_0_xxxxyzz, g_y_x_0_xxxxzz, g_y_x_0_xxxxzzz, g_y_x_0_xxxyyy, g_y_x_0_xxxyyyy, g_y_x_0_xxxyyyz, g_y_x_0_xxxyyz, g_y_x_0_xxxyyzz, g_y_x_0_xxxyzz, g_y_x_0_xxxyzzz, g_y_x_0_xxxzzz, g_y_x_0_xxxzzzz, g_y_x_0_xxyyyy, g_y_x_0_xxyyyyy, g_y_x_0_xxyyyyz, g_y_x_0_xxyyyz, g_y_x_0_xxyyyzz, g_y_x_0_xxyyzz, g_y_x_0_xxyyzzz, g_y_x_0_xxyzzz, g_y_x_0_xxyzzzz, g_y_x_0_xxzzzz, g_y_x_0_xxzzzzz, g_y_x_0_xyyyyy, g_y_x_0_xyyyyyy, g_y_x_0_xyyyyyz, g_y_x_0_xyyyyz, g_y_x_0_xyyyyzz, g_y_x_0_xyyyzz, g_y_x_0_xyyyzzz, g_y_x_0_xyyzzz, g_y_x_0_xyyzzzz, g_y_x_0_xyzzzz, g_y_x_0_xyzzzzz, g_y_x_0_xzzzzz, g_y_x_0_xzzzzzz, g_y_x_0_yyyyyy, g_y_x_0_yyyyyz, g_y_x_0_yyyyzz, g_y_x_0_yyyzzz, g_y_x_0_yyzzzz, g_y_x_0_yzzzzz, g_y_x_0_zzzzzz, g_y_x_x_xxxxxx, g_y_x_x_xxxxxy, g_y_x_x_xxxxxz, g_y_x_x_xxxxyy, g_y_x_x_xxxxyz, g_y_x_x_xxxxzz, g_y_x_x_xxxyyy, g_y_x_x_xxxyyz, g_y_x_x_xxxyzz, g_y_x_x_xxxzzz, g_y_x_x_xxyyyy, g_y_x_x_xxyyyz, g_y_x_x_xxyyzz, g_y_x_x_xxyzzz, g_y_x_x_xxzzzz, g_y_x_x_xyyyyy, g_y_x_x_xyyyyz, g_y_x_x_xyyyzz, g_y_x_x_xyyzzz, g_y_x_x_xyzzzz, g_y_x_x_xzzzzz, g_y_x_x_yyyyyy, g_y_x_x_yyyyyz, g_y_x_x_yyyyzz, g_y_x_x_yyyzzz, g_y_x_x_yyzzzz, g_y_x_x_yzzzzz, g_y_x_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_x_xxxxxx[k] = g_y_0_0_xxxxxx[k] - g_y_x_0_xxxxxx[k] * ab_x + g_y_x_0_xxxxxxx[k];

                g_y_x_x_xxxxxy[k] = g_y_0_0_xxxxxy[k] - g_y_x_0_xxxxxy[k] * ab_x + g_y_x_0_xxxxxxy[k];

                g_y_x_x_xxxxxz[k] = g_y_0_0_xxxxxz[k] - g_y_x_0_xxxxxz[k] * ab_x + g_y_x_0_xxxxxxz[k];

                g_y_x_x_xxxxyy[k] = g_y_0_0_xxxxyy[k] - g_y_x_0_xxxxyy[k] * ab_x + g_y_x_0_xxxxxyy[k];

                g_y_x_x_xxxxyz[k] = g_y_0_0_xxxxyz[k] - g_y_x_0_xxxxyz[k] * ab_x + g_y_x_0_xxxxxyz[k];

                g_y_x_x_xxxxzz[k] = g_y_0_0_xxxxzz[k] - g_y_x_0_xxxxzz[k] * ab_x + g_y_x_0_xxxxxzz[k];

                g_y_x_x_xxxyyy[k] = g_y_0_0_xxxyyy[k] - g_y_x_0_xxxyyy[k] * ab_x + g_y_x_0_xxxxyyy[k];

                g_y_x_x_xxxyyz[k] = g_y_0_0_xxxyyz[k] - g_y_x_0_xxxyyz[k] * ab_x + g_y_x_0_xxxxyyz[k];

                g_y_x_x_xxxyzz[k] = g_y_0_0_xxxyzz[k] - g_y_x_0_xxxyzz[k] * ab_x + g_y_x_0_xxxxyzz[k];

                g_y_x_x_xxxzzz[k] = g_y_0_0_xxxzzz[k] - g_y_x_0_xxxzzz[k] * ab_x + g_y_x_0_xxxxzzz[k];

                g_y_x_x_xxyyyy[k] = g_y_0_0_xxyyyy[k] - g_y_x_0_xxyyyy[k] * ab_x + g_y_x_0_xxxyyyy[k];

                g_y_x_x_xxyyyz[k] = g_y_0_0_xxyyyz[k] - g_y_x_0_xxyyyz[k] * ab_x + g_y_x_0_xxxyyyz[k];

                g_y_x_x_xxyyzz[k] = g_y_0_0_xxyyzz[k] - g_y_x_0_xxyyzz[k] * ab_x + g_y_x_0_xxxyyzz[k];

                g_y_x_x_xxyzzz[k] = g_y_0_0_xxyzzz[k] - g_y_x_0_xxyzzz[k] * ab_x + g_y_x_0_xxxyzzz[k];

                g_y_x_x_xxzzzz[k] = g_y_0_0_xxzzzz[k] - g_y_x_0_xxzzzz[k] * ab_x + g_y_x_0_xxxzzzz[k];

                g_y_x_x_xyyyyy[k] = g_y_0_0_xyyyyy[k] - g_y_x_0_xyyyyy[k] * ab_x + g_y_x_0_xxyyyyy[k];

                g_y_x_x_xyyyyz[k] = g_y_0_0_xyyyyz[k] - g_y_x_0_xyyyyz[k] * ab_x + g_y_x_0_xxyyyyz[k];

                g_y_x_x_xyyyzz[k] = g_y_0_0_xyyyzz[k] - g_y_x_0_xyyyzz[k] * ab_x + g_y_x_0_xxyyyzz[k];

                g_y_x_x_xyyzzz[k] = g_y_0_0_xyyzzz[k] - g_y_x_0_xyyzzz[k] * ab_x + g_y_x_0_xxyyzzz[k];

                g_y_x_x_xyzzzz[k] = g_y_0_0_xyzzzz[k] - g_y_x_0_xyzzzz[k] * ab_x + g_y_x_0_xxyzzzz[k];

                g_y_x_x_xzzzzz[k] = g_y_0_0_xzzzzz[k] - g_y_x_0_xzzzzz[k] * ab_x + g_y_x_0_xxzzzzz[k];

                g_y_x_x_yyyyyy[k] = g_y_0_0_yyyyyy[k] - g_y_x_0_yyyyyy[k] * ab_x + g_y_x_0_xyyyyyy[k];

                g_y_x_x_yyyyyz[k] = g_y_0_0_yyyyyz[k] - g_y_x_0_yyyyyz[k] * ab_x + g_y_x_0_xyyyyyz[k];

                g_y_x_x_yyyyzz[k] = g_y_0_0_yyyyzz[k] - g_y_x_0_yyyyzz[k] * ab_x + g_y_x_0_xyyyyzz[k];

                g_y_x_x_yyyzzz[k] = g_y_0_0_yyyzzz[k] - g_y_x_0_yyyzzz[k] * ab_x + g_y_x_0_xyyyzzz[k];

                g_y_x_x_yyzzzz[k] = g_y_0_0_yyzzzz[k] - g_y_x_0_yyzzzz[k] * ab_x + g_y_x_0_xyyzzzz[k];

                g_y_x_x_yzzzzz[k] = g_y_0_0_yzzzzz[k] - g_y_x_0_yzzzzz[k] * ab_x + g_y_x_0_xyzzzzz[k];

                g_y_x_x_zzzzzz[k] = g_y_0_0_zzzzzz[k] - g_y_x_0_zzzzzz[k] * ab_x + g_y_x_0_xzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_y_x_y_xxxxxx = cbuffer.data(pi_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_x_y_xxxxxy = cbuffer.data(pi_geom_11_off + 281 * ccomps * dcomps);

            auto g_y_x_y_xxxxxz = cbuffer.data(pi_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_x_y_xxxxyy = cbuffer.data(pi_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_x_y_xxxxyz = cbuffer.data(pi_geom_11_off + 284 * ccomps * dcomps);

            auto g_y_x_y_xxxxzz = cbuffer.data(pi_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_x_y_xxxyyy = cbuffer.data(pi_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_x_y_xxxyyz = cbuffer.data(pi_geom_11_off + 287 * ccomps * dcomps);

            auto g_y_x_y_xxxyzz = cbuffer.data(pi_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_x_y_xxxzzz = cbuffer.data(pi_geom_11_off + 289 * ccomps * dcomps);

            auto g_y_x_y_xxyyyy = cbuffer.data(pi_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_x_y_xxyyyz = cbuffer.data(pi_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_x_y_xxyyzz = cbuffer.data(pi_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_x_y_xxyzzz = cbuffer.data(pi_geom_11_off + 293 * ccomps * dcomps);

            auto g_y_x_y_xxzzzz = cbuffer.data(pi_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_x_y_xyyyyy = cbuffer.data(pi_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_x_y_xyyyyz = cbuffer.data(pi_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_x_y_xyyyzz = cbuffer.data(pi_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_x_y_xyyzzz = cbuffer.data(pi_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_x_y_xyzzzz = cbuffer.data(pi_geom_11_off + 299 * ccomps * dcomps);

            auto g_y_x_y_xzzzzz = cbuffer.data(pi_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_x_y_yyyyyy = cbuffer.data(pi_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_x_y_yyyyyz = cbuffer.data(pi_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_x_y_yyyyzz = cbuffer.data(pi_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_x_y_yyyzzz = cbuffer.data(pi_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_x_y_yyzzzz = cbuffer.data(pi_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_x_y_yzzzzz = cbuffer.data(pi_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_x_y_zzzzzz = cbuffer.data(pi_geom_11_off + 307 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxx, g_0_x_0_xxxxxy, g_0_x_0_xxxxxz, g_0_x_0_xxxxyy, g_0_x_0_xxxxyz, g_0_x_0_xxxxzz, g_0_x_0_xxxyyy, g_0_x_0_xxxyyz, g_0_x_0_xxxyzz, g_0_x_0_xxxzzz, g_0_x_0_xxyyyy, g_0_x_0_xxyyyz, g_0_x_0_xxyyzz, g_0_x_0_xxyzzz, g_0_x_0_xxzzzz, g_0_x_0_xyyyyy, g_0_x_0_xyyyyz, g_0_x_0_xyyyzz, g_0_x_0_xyyzzz, g_0_x_0_xyzzzz, g_0_x_0_xzzzzz, g_0_x_0_yyyyyy, g_0_x_0_yyyyyz, g_0_x_0_yyyyzz, g_0_x_0_yyyzzz, g_0_x_0_yyzzzz, g_0_x_0_yzzzzz, g_0_x_0_zzzzzz, g_y_x_0_xxxxxx, g_y_x_0_xxxxxxy, g_y_x_0_xxxxxy, g_y_x_0_xxxxxyy, g_y_x_0_xxxxxyz, g_y_x_0_xxxxxz, g_y_x_0_xxxxyy, g_y_x_0_xxxxyyy, g_y_x_0_xxxxyyz, g_y_x_0_xxxxyz, g_y_x_0_xxxxyzz, g_y_x_0_xxxxzz, g_y_x_0_xxxyyy, g_y_x_0_xxxyyyy, g_y_x_0_xxxyyyz, g_y_x_0_xxxyyz, g_y_x_0_xxxyyzz, g_y_x_0_xxxyzz, g_y_x_0_xxxyzzz, g_y_x_0_xxxzzz, g_y_x_0_xxyyyy, g_y_x_0_xxyyyyy, g_y_x_0_xxyyyyz, g_y_x_0_xxyyyz, g_y_x_0_xxyyyzz, g_y_x_0_xxyyzz, g_y_x_0_xxyyzzz, g_y_x_0_xxyzzz, g_y_x_0_xxyzzzz, g_y_x_0_xxzzzz, g_y_x_0_xyyyyy, g_y_x_0_xyyyyyy, g_y_x_0_xyyyyyz, g_y_x_0_xyyyyz, g_y_x_0_xyyyyzz, g_y_x_0_xyyyzz, g_y_x_0_xyyyzzz, g_y_x_0_xyyzzz, g_y_x_0_xyyzzzz, g_y_x_0_xyzzzz, g_y_x_0_xyzzzzz, g_y_x_0_xzzzzz, g_y_x_0_yyyyyy, g_y_x_0_yyyyyyy, g_y_x_0_yyyyyyz, g_y_x_0_yyyyyz, g_y_x_0_yyyyyzz, g_y_x_0_yyyyzz, g_y_x_0_yyyyzzz, g_y_x_0_yyyzzz, g_y_x_0_yyyzzzz, g_y_x_0_yyzzzz, g_y_x_0_yyzzzzz, g_y_x_0_yzzzzz, g_y_x_0_yzzzzzz, g_y_x_0_zzzzzz, g_y_x_y_xxxxxx, g_y_x_y_xxxxxy, g_y_x_y_xxxxxz, g_y_x_y_xxxxyy, g_y_x_y_xxxxyz, g_y_x_y_xxxxzz, g_y_x_y_xxxyyy, g_y_x_y_xxxyyz, g_y_x_y_xxxyzz, g_y_x_y_xxxzzz, g_y_x_y_xxyyyy, g_y_x_y_xxyyyz, g_y_x_y_xxyyzz, g_y_x_y_xxyzzz, g_y_x_y_xxzzzz, g_y_x_y_xyyyyy, g_y_x_y_xyyyyz, g_y_x_y_xyyyzz, g_y_x_y_xyyzzz, g_y_x_y_xyzzzz, g_y_x_y_xzzzzz, g_y_x_y_yyyyyy, g_y_x_y_yyyyyz, g_y_x_y_yyyyzz, g_y_x_y_yyyzzz, g_y_x_y_yyzzzz, g_y_x_y_yzzzzz, g_y_x_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_y_xxxxxx[k] = -g_0_x_0_xxxxxx[k] - g_y_x_0_xxxxxx[k] * ab_y + g_y_x_0_xxxxxxy[k];

                g_y_x_y_xxxxxy[k] = -g_0_x_0_xxxxxy[k] - g_y_x_0_xxxxxy[k] * ab_y + g_y_x_0_xxxxxyy[k];

                g_y_x_y_xxxxxz[k] = -g_0_x_0_xxxxxz[k] - g_y_x_0_xxxxxz[k] * ab_y + g_y_x_0_xxxxxyz[k];

                g_y_x_y_xxxxyy[k] = -g_0_x_0_xxxxyy[k] - g_y_x_0_xxxxyy[k] * ab_y + g_y_x_0_xxxxyyy[k];

                g_y_x_y_xxxxyz[k] = -g_0_x_0_xxxxyz[k] - g_y_x_0_xxxxyz[k] * ab_y + g_y_x_0_xxxxyyz[k];

                g_y_x_y_xxxxzz[k] = -g_0_x_0_xxxxzz[k] - g_y_x_0_xxxxzz[k] * ab_y + g_y_x_0_xxxxyzz[k];

                g_y_x_y_xxxyyy[k] = -g_0_x_0_xxxyyy[k] - g_y_x_0_xxxyyy[k] * ab_y + g_y_x_0_xxxyyyy[k];

                g_y_x_y_xxxyyz[k] = -g_0_x_0_xxxyyz[k] - g_y_x_0_xxxyyz[k] * ab_y + g_y_x_0_xxxyyyz[k];

                g_y_x_y_xxxyzz[k] = -g_0_x_0_xxxyzz[k] - g_y_x_0_xxxyzz[k] * ab_y + g_y_x_0_xxxyyzz[k];

                g_y_x_y_xxxzzz[k] = -g_0_x_0_xxxzzz[k] - g_y_x_0_xxxzzz[k] * ab_y + g_y_x_0_xxxyzzz[k];

                g_y_x_y_xxyyyy[k] = -g_0_x_0_xxyyyy[k] - g_y_x_0_xxyyyy[k] * ab_y + g_y_x_0_xxyyyyy[k];

                g_y_x_y_xxyyyz[k] = -g_0_x_0_xxyyyz[k] - g_y_x_0_xxyyyz[k] * ab_y + g_y_x_0_xxyyyyz[k];

                g_y_x_y_xxyyzz[k] = -g_0_x_0_xxyyzz[k] - g_y_x_0_xxyyzz[k] * ab_y + g_y_x_0_xxyyyzz[k];

                g_y_x_y_xxyzzz[k] = -g_0_x_0_xxyzzz[k] - g_y_x_0_xxyzzz[k] * ab_y + g_y_x_0_xxyyzzz[k];

                g_y_x_y_xxzzzz[k] = -g_0_x_0_xxzzzz[k] - g_y_x_0_xxzzzz[k] * ab_y + g_y_x_0_xxyzzzz[k];

                g_y_x_y_xyyyyy[k] = -g_0_x_0_xyyyyy[k] - g_y_x_0_xyyyyy[k] * ab_y + g_y_x_0_xyyyyyy[k];

                g_y_x_y_xyyyyz[k] = -g_0_x_0_xyyyyz[k] - g_y_x_0_xyyyyz[k] * ab_y + g_y_x_0_xyyyyyz[k];

                g_y_x_y_xyyyzz[k] = -g_0_x_0_xyyyzz[k] - g_y_x_0_xyyyzz[k] * ab_y + g_y_x_0_xyyyyzz[k];

                g_y_x_y_xyyzzz[k] = -g_0_x_0_xyyzzz[k] - g_y_x_0_xyyzzz[k] * ab_y + g_y_x_0_xyyyzzz[k];

                g_y_x_y_xyzzzz[k] = -g_0_x_0_xyzzzz[k] - g_y_x_0_xyzzzz[k] * ab_y + g_y_x_0_xyyzzzz[k];

                g_y_x_y_xzzzzz[k] = -g_0_x_0_xzzzzz[k] - g_y_x_0_xzzzzz[k] * ab_y + g_y_x_0_xyzzzzz[k];

                g_y_x_y_yyyyyy[k] = -g_0_x_0_yyyyyy[k] - g_y_x_0_yyyyyy[k] * ab_y + g_y_x_0_yyyyyyy[k];

                g_y_x_y_yyyyyz[k] = -g_0_x_0_yyyyyz[k] - g_y_x_0_yyyyyz[k] * ab_y + g_y_x_0_yyyyyyz[k];

                g_y_x_y_yyyyzz[k] = -g_0_x_0_yyyyzz[k] - g_y_x_0_yyyyzz[k] * ab_y + g_y_x_0_yyyyyzz[k];

                g_y_x_y_yyyzzz[k] = -g_0_x_0_yyyzzz[k] - g_y_x_0_yyyzzz[k] * ab_y + g_y_x_0_yyyyzzz[k];

                g_y_x_y_yyzzzz[k] = -g_0_x_0_yyzzzz[k] - g_y_x_0_yyzzzz[k] * ab_y + g_y_x_0_yyyzzzz[k];

                g_y_x_y_yzzzzz[k] = -g_0_x_0_yzzzzz[k] - g_y_x_0_yzzzzz[k] * ab_y + g_y_x_0_yyzzzzz[k];

                g_y_x_y_zzzzzz[k] = -g_0_x_0_zzzzzz[k] - g_y_x_0_zzzzzz[k] * ab_y + g_y_x_0_yzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_y_x_z_xxxxxx = cbuffer.data(pi_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_x_z_xxxxxy = cbuffer.data(pi_geom_11_off + 309 * ccomps * dcomps);

            auto g_y_x_z_xxxxxz = cbuffer.data(pi_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_x_z_xxxxyy = cbuffer.data(pi_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_x_z_xxxxyz = cbuffer.data(pi_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_x_z_xxxxzz = cbuffer.data(pi_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_x_z_xxxyyy = cbuffer.data(pi_geom_11_off + 314 * ccomps * dcomps);

            auto g_y_x_z_xxxyyz = cbuffer.data(pi_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_x_z_xxxyzz = cbuffer.data(pi_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_x_z_xxxzzz = cbuffer.data(pi_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_x_z_xxyyyy = cbuffer.data(pi_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_x_z_xxyyyz = cbuffer.data(pi_geom_11_off + 319 * ccomps * dcomps);

            auto g_y_x_z_xxyyzz = cbuffer.data(pi_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_x_z_xxyzzz = cbuffer.data(pi_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_x_z_xxzzzz = cbuffer.data(pi_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_x_z_xyyyyy = cbuffer.data(pi_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_x_z_xyyyyz = cbuffer.data(pi_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_x_z_xyyyzz = cbuffer.data(pi_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_x_z_xyyzzz = cbuffer.data(pi_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_x_z_xyzzzz = cbuffer.data(pi_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_x_z_xzzzzz = cbuffer.data(pi_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_x_z_yyyyyy = cbuffer.data(pi_geom_11_off + 329 * ccomps * dcomps);

            auto g_y_x_z_yyyyyz = cbuffer.data(pi_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_x_z_yyyyzz = cbuffer.data(pi_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_x_z_yyyzzz = cbuffer.data(pi_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_x_z_yyzzzz = cbuffer.data(pi_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_x_z_yzzzzz = cbuffer.data(pi_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_x_z_zzzzzz = cbuffer.data(pi_geom_11_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_0_xxxxxx, g_y_x_0_xxxxxxz, g_y_x_0_xxxxxy, g_y_x_0_xxxxxyz, g_y_x_0_xxxxxz, g_y_x_0_xxxxxzz, g_y_x_0_xxxxyy, g_y_x_0_xxxxyyz, g_y_x_0_xxxxyz, g_y_x_0_xxxxyzz, g_y_x_0_xxxxzz, g_y_x_0_xxxxzzz, g_y_x_0_xxxyyy, g_y_x_0_xxxyyyz, g_y_x_0_xxxyyz, g_y_x_0_xxxyyzz, g_y_x_0_xxxyzz, g_y_x_0_xxxyzzz, g_y_x_0_xxxzzz, g_y_x_0_xxxzzzz, g_y_x_0_xxyyyy, g_y_x_0_xxyyyyz, g_y_x_0_xxyyyz, g_y_x_0_xxyyyzz, g_y_x_0_xxyyzz, g_y_x_0_xxyyzzz, g_y_x_0_xxyzzz, g_y_x_0_xxyzzzz, g_y_x_0_xxzzzz, g_y_x_0_xxzzzzz, g_y_x_0_xyyyyy, g_y_x_0_xyyyyyz, g_y_x_0_xyyyyz, g_y_x_0_xyyyyzz, g_y_x_0_xyyyzz, g_y_x_0_xyyyzzz, g_y_x_0_xyyzzz, g_y_x_0_xyyzzzz, g_y_x_0_xyzzzz, g_y_x_0_xyzzzzz, g_y_x_0_xzzzzz, g_y_x_0_xzzzzzz, g_y_x_0_yyyyyy, g_y_x_0_yyyyyyz, g_y_x_0_yyyyyz, g_y_x_0_yyyyyzz, g_y_x_0_yyyyzz, g_y_x_0_yyyyzzz, g_y_x_0_yyyzzz, g_y_x_0_yyyzzzz, g_y_x_0_yyzzzz, g_y_x_0_yyzzzzz, g_y_x_0_yzzzzz, g_y_x_0_yzzzzzz, g_y_x_0_zzzzzz, g_y_x_0_zzzzzzz, g_y_x_z_xxxxxx, g_y_x_z_xxxxxy, g_y_x_z_xxxxxz, g_y_x_z_xxxxyy, g_y_x_z_xxxxyz, g_y_x_z_xxxxzz, g_y_x_z_xxxyyy, g_y_x_z_xxxyyz, g_y_x_z_xxxyzz, g_y_x_z_xxxzzz, g_y_x_z_xxyyyy, g_y_x_z_xxyyyz, g_y_x_z_xxyyzz, g_y_x_z_xxyzzz, g_y_x_z_xxzzzz, g_y_x_z_xyyyyy, g_y_x_z_xyyyyz, g_y_x_z_xyyyzz, g_y_x_z_xyyzzz, g_y_x_z_xyzzzz, g_y_x_z_xzzzzz, g_y_x_z_yyyyyy, g_y_x_z_yyyyyz, g_y_x_z_yyyyzz, g_y_x_z_yyyzzz, g_y_x_z_yyzzzz, g_y_x_z_yzzzzz, g_y_x_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_z_xxxxxx[k] = -g_y_x_0_xxxxxx[k] * ab_z + g_y_x_0_xxxxxxz[k];

                g_y_x_z_xxxxxy[k] = -g_y_x_0_xxxxxy[k] * ab_z + g_y_x_0_xxxxxyz[k];

                g_y_x_z_xxxxxz[k] = -g_y_x_0_xxxxxz[k] * ab_z + g_y_x_0_xxxxxzz[k];

                g_y_x_z_xxxxyy[k] = -g_y_x_0_xxxxyy[k] * ab_z + g_y_x_0_xxxxyyz[k];

                g_y_x_z_xxxxyz[k] = -g_y_x_0_xxxxyz[k] * ab_z + g_y_x_0_xxxxyzz[k];

                g_y_x_z_xxxxzz[k] = -g_y_x_0_xxxxzz[k] * ab_z + g_y_x_0_xxxxzzz[k];

                g_y_x_z_xxxyyy[k] = -g_y_x_0_xxxyyy[k] * ab_z + g_y_x_0_xxxyyyz[k];

                g_y_x_z_xxxyyz[k] = -g_y_x_0_xxxyyz[k] * ab_z + g_y_x_0_xxxyyzz[k];

                g_y_x_z_xxxyzz[k] = -g_y_x_0_xxxyzz[k] * ab_z + g_y_x_0_xxxyzzz[k];

                g_y_x_z_xxxzzz[k] = -g_y_x_0_xxxzzz[k] * ab_z + g_y_x_0_xxxzzzz[k];

                g_y_x_z_xxyyyy[k] = -g_y_x_0_xxyyyy[k] * ab_z + g_y_x_0_xxyyyyz[k];

                g_y_x_z_xxyyyz[k] = -g_y_x_0_xxyyyz[k] * ab_z + g_y_x_0_xxyyyzz[k];

                g_y_x_z_xxyyzz[k] = -g_y_x_0_xxyyzz[k] * ab_z + g_y_x_0_xxyyzzz[k];

                g_y_x_z_xxyzzz[k] = -g_y_x_0_xxyzzz[k] * ab_z + g_y_x_0_xxyzzzz[k];

                g_y_x_z_xxzzzz[k] = -g_y_x_0_xxzzzz[k] * ab_z + g_y_x_0_xxzzzzz[k];

                g_y_x_z_xyyyyy[k] = -g_y_x_0_xyyyyy[k] * ab_z + g_y_x_0_xyyyyyz[k];

                g_y_x_z_xyyyyz[k] = -g_y_x_0_xyyyyz[k] * ab_z + g_y_x_0_xyyyyzz[k];

                g_y_x_z_xyyyzz[k] = -g_y_x_0_xyyyzz[k] * ab_z + g_y_x_0_xyyyzzz[k];

                g_y_x_z_xyyzzz[k] = -g_y_x_0_xyyzzz[k] * ab_z + g_y_x_0_xyyzzzz[k];

                g_y_x_z_xyzzzz[k] = -g_y_x_0_xyzzzz[k] * ab_z + g_y_x_0_xyzzzzz[k];

                g_y_x_z_xzzzzz[k] = -g_y_x_0_xzzzzz[k] * ab_z + g_y_x_0_xzzzzzz[k];

                g_y_x_z_yyyyyy[k] = -g_y_x_0_yyyyyy[k] * ab_z + g_y_x_0_yyyyyyz[k];

                g_y_x_z_yyyyyz[k] = -g_y_x_0_yyyyyz[k] * ab_z + g_y_x_0_yyyyyzz[k];

                g_y_x_z_yyyyzz[k] = -g_y_x_0_yyyyzz[k] * ab_z + g_y_x_0_yyyyzzz[k];

                g_y_x_z_yyyzzz[k] = -g_y_x_0_yyyzzz[k] * ab_z + g_y_x_0_yyyzzzz[k];

                g_y_x_z_yyzzzz[k] = -g_y_x_0_yyzzzz[k] * ab_z + g_y_x_0_yyzzzzz[k];

                g_y_x_z_yzzzzz[k] = -g_y_x_0_yzzzzz[k] * ab_z + g_y_x_0_yzzzzzz[k];

                g_y_x_z_zzzzzz[k] = -g_y_x_0_zzzzzz[k] * ab_z + g_y_x_0_zzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_y_y_x_xxxxxx = cbuffer.data(pi_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_y_x_xxxxxy = cbuffer.data(pi_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_y_x_xxxxxz = cbuffer.data(pi_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_y_x_xxxxyy = cbuffer.data(pi_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_y_x_xxxxyz = cbuffer.data(pi_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_y_x_xxxxzz = cbuffer.data(pi_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_y_x_xxxyyy = cbuffer.data(pi_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_y_x_xxxyyz = cbuffer.data(pi_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_y_x_xxxyzz = cbuffer.data(pi_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_y_x_xxxzzz = cbuffer.data(pi_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_y_x_xxyyyy = cbuffer.data(pi_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_y_x_xxyyyz = cbuffer.data(pi_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_y_x_xxyyzz = cbuffer.data(pi_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_y_x_xxyzzz = cbuffer.data(pi_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_y_x_xxzzzz = cbuffer.data(pi_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_y_x_xyyyyy = cbuffer.data(pi_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_y_x_xyyyyz = cbuffer.data(pi_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_y_x_xyyyzz = cbuffer.data(pi_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_y_x_xyyzzz = cbuffer.data(pi_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_y_x_xyzzzz = cbuffer.data(pi_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_y_x_xzzzzz = cbuffer.data(pi_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_y_x_yyyyyy = cbuffer.data(pi_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_y_x_yyyyyz = cbuffer.data(pi_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_y_x_yyyyzz = cbuffer.data(pi_geom_11_off + 359 * ccomps * dcomps);

            auto g_y_y_x_yyyzzz = cbuffer.data(pi_geom_11_off + 360 * ccomps * dcomps);

            auto g_y_y_x_yyzzzz = cbuffer.data(pi_geom_11_off + 361 * ccomps * dcomps);

            auto g_y_y_x_yzzzzz = cbuffer.data(pi_geom_11_off + 362 * ccomps * dcomps);

            auto g_y_y_x_zzzzzz = cbuffer.data(pi_geom_11_off + 363 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_0_xxxxxx, g_y_y_0_xxxxxxx, g_y_y_0_xxxxxxy, g_y_y_0_xxxxxxz, g_y_y_0_xxxxxy, g_y_y_0_xxxxxyy, g_y_y_0_xxxxxyz, g_y_y_0_xxxxxz, g_y_y_0_xxxxxzz, g_y_y_0_xxxxyy, g_y_y_0_xxxxyyy, g_y_y_0_xxxxyyz, g_y_y_0_xxxxyz, g_y_y_0_xxxxyzz, g_y_y_0_xxxxzz, g_y_y_0_xxxxzzz, g_y_y_0_xxxyyy, g_y_y_0_xxxyyyy, g_y_y_0_xxxyyyz, g_y_y_0_xxxyyz, g_y_y_0_xxxyyzz, g_y_y_0_xxxyzz, g_y_y_0_xxxyzzz, g_y_y_0_xxxzzz, g_y_y_0_xxxzzzz, g_y_y_0_xxyyyy, g_y_y_0_xxyyyyy, g_y_y_0_xxyyyyz, g_y_y_0_xxyyyz, g_y_y_0_xxyyyzz, g_y_y_0_xxyyzz, g_y_y_0_xxyyzzz, g_y_y_0_xxyzzz, g_y_y_0_xxyzzzz, g_y_y_0_xxzzzz, g_y_y_0_xxzzzzz, g_y_y_0_xyyyyy, g_y_y_0_xyyyyyy, g_y_y_0_xyyyyyz, g_y_y_0_xyyyyz, g_y_y_0_xyyyyzz, g_y_y_0_xyyyzz, g_y_y_0_xyyyzzz, g_y_y_0_xyyzzz, g_y_y_0_xyyzzzz, g_y_y_0_xyzzzz, g_y_y_0_xyzzzzz, g_y_y_0_xzzzzz, g_y_y_0_xzzzzzz, g_y_y_0_yyyyyy, g_y_y_0_yyyyyz, g_y_y_0_yyyyzz, g_y_y_0_yyyzzz, g_y_y_0_yyzzzz, g_y_y_0_yzzzzz, g_y_y_0_zzzzzz, g_y_y_x_xxxxxx, g_y_y_x_xxxxxy, g_y_y_x_xxxxxz, g_y_y_x_xxxxyy, g_y_y_x_xxxxyz, g_y_y_x_xxxxzz, g_y_y_x_xxxyyy, g_y_y_x_xxxyyz, g_y_y_x_xxxyzz, g_y_y_x_xxxzzz, g_y_y_x_xxyyyy, g_y_y_x_xxyyyz, g_y_y_x_xxyyzz, g_y_y_x_xxyzzz, g_y_y_x_xxzzzz, g_y_y_x_xyyyyy, g_y_y_x_xyyyyz, g_y_y_x_xyyyzz, g_y_y_x_xyyzzz, g_y_y_x_xyzzzz, g_y_y_x_xzzzzz, g_y_y_x_yyyyyy, g_y_y_x_yyyyyz, g_y_y_x_yyyyzz, g_y_y_x_yyyzzz, g_y_y_x_yyzzzz, g_y_y_x_yzzzzz, g_y_y_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_x_xxxxxx[k] = -g_y_y_0_xxxxxx[k] * ab_x + g_y_y_0_xxxxxxx[k];

                g_y_y_x_xxxxxy[k] = -g_y_y_0_xxxxxy[k] * ab_x + g_y_y_0_xxxxxxy[k];

                g_y_y_x_xxxxxz[k] = -g_y_y_0_xxxxxz[k] * ab_x + g_y_y_0_xxxxxxz[k];

                g_y_y_x_xxxxyy[k] = -g_y_y_0_xxxxyy[k] * ab_x + g_y_y_0_xxxxxyy[k];

                g_y_y_x_xxxxyz[k] = -g_y_y_0_xxxxyz[k] * ab_x + g_y_y_0_xxxxxyz[k];

                g_y_y_x_xxxxzz[k] = -g_y_y_0_xxxxzz[k] * ab_x + g_y_y_0_xxxxxzz[k];

                g_y_y_x_xxxyyy[k] = -g_y_y_0_xxxyyy[k] * ab_x + g_y_y_0_xxxxyyy[k];

                g_y_y_x_xxxyyz[k] = -g_y_y_0_xxxyyz[k] * ab_x + g_y_y_0_xxxxyyz[k];

                g_y_y_x_xxxyzz[k] = -g_y_y_0_xxxyzz[k] * ab_x + g_y_y_0_xxxxyzz[k];

                g_y_y_x_xxxzzz[k] = -g_y_y_0_xxxzzz[k] * ab_x + g_y_y_0_xxxxzzz[k];

                g_y_y_x_xxyyyy[k] = -g_y_y_0_xxyyyy[k] * ab_x + g_y_y_0_xxxyyyy[k];

                g_y_y_x_xxyyyz[k] = -g_y_y_0_xxyyyz[k] * ab_x + g_y_y_0_xxxyyyz[k];

                g_y_y_x_xxyyzz[k] = -g_y_y_0_xxyyzz[k] * ab_x + g_y_y_0_xxxyyzz[k];

                g_y_y_x_xxyzzz[k] = -g_y_y_0_xxyzzz[k] * ab_x + g_y_y_0_xxxyzzz[k];

                g_y_y_x_xxzzzz[k] = -g_y_y_0_xxzzzz[k] * ab_x + g_y_y_0_xxxzzzz[k];

                g_y_y_x_xyyyyy[k] = -g_y_y_0_xyyyyy[k] * ab_x + g_y_y_0_xxyyyyy[k];

                g_y_y_x_xyyyyz[k] = -g_y_y_0_xyyyyz[k] * ab_x + g_y_y_0_xxyyyyz[k];

                g_y_y_x_xyyyzz[k] = -g_y_y_0_xyyyzz[k] * ab_x + g_y_y_0_xxyyyzz[k];

                g_y_y_x_xyyzzz[k] = -g_y_y_0_xyyzzz[k] * ab_x + g_y_y_0_xxyyzzz[k];

                g_y_y_x_xyzzzz[k] = -g_y_y_0_xyzzzz[k] * ab_x + g_y_y_0_xxyzzzz[k];

                g_y_y_x_xzzzzz[k] = -g_y_y_0_xzzzzz[k] * ab_x + g_y_y_0_xxzzzzz[k];

                g_y_y_x_yyyyyy[k] = -g_y_y_0_yyyyyy[k] * ab_x + g_y_y_0_xyyyyyy[k];

                g_y_y_x_yyyyyz[k] = -g_y_y_0_yyyyyz[k] * ab_x + g_y_y_0_xyyyyyz[k];

                g_y_y_x_yyyyzz[k] = -g_y_y_0_yyyyzz[k] * ab_x + g_y_y_0_xyyyyzz[k];

                g_y_y_x_yyyzzz[k] = -g_y_y_0_yyyzzz[k] * ab_x + g_y_y_0_xyyyzzz[k];

                g_y_y_x_yyzzzz[k] = -g_y_y_0_yyzzzz[k] * ab_x + g_y_y_0_xyyzzzz[k];

                g_y_y_x_yzzzzz[k] = -g_y_y_0_yzzzzz[k] * ab_x + g_y_y_0_xyzzzzz[k];

                g_y_y_x_zzzzzz[k] = -g_y_y_0_zzzzzz[k] * ab_x + g_y_y_0_xzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_y_y_y_xxxxxx = cbuffer.data(pi_geom_11_off + 364 * ccomps * dcomps);

            auto g_y_y_y_xxxxxy = cbuffer.data(pi_geom_11_off + 365 * ccomps * dcomps);

            auto g_y_y_y_xxxxxz = cbuffer.data(pi_geom_11_off + 366 * ccomps * dcomps);

            auto g_y_y_y_xxxxyy = cbuffer.data(pi_geom_11_off + 367 * ccomps * dcomps);

            auto g_y_y_y_xxxxyz = cbuffer.data(pi_geom_11_off + 368 * ccomps * dcomps);

            auto g_y_y_y_xxxxzz = cbuffer.data(pi_geom_11_off + 369 * ccomps * dcomps);

            auto g_y_y_y_xxxyyy = cbuffer.data(pi_geom_11_off + 370 * ccomps * dcomps);

            auto g_y_y_y_xxxyyz = cbuffer.data(pi_geom_11_off + 371 * ccomps * dcomps);

            auto g_y_y_y_xxxyzz = cbuffer.data(pi_geom_11_off + 372 * ccomps * dcomps);

            auto g_y_y_y_xxxzzz = cbuffer.data(pi_geom_11_off + 373 * ccomps * dcomps);

            auto g_y_y_y_xxyyyy = cbuffer.data(pi_geom_11_off + 374 * ccomps * dcomps);

            auto g_y_y_y_xxyyyz = cbuffer.data(pi_geom_11_off + 375 * ccomps * dcomps);

            auto g_y_y_y_xxyyzz = cbuffer.data(pi_geom_11_off + 376 * ccomps * dcomps);

            auto g_y_y_y_xxyzzz = cbuffer.data(pi_geom_11_off + 377 * ccomps * dcomps);

            auto g_y_y_y_xxzzzz = cbuffer.data(pi_geom_11_off + 378 * ccomps * dcomps);

            auto g_y_y_y_xyyyyy = cbuffer.data(pi_geom_11_off + 379 * ccomps * dcomps);

            auto g_y_y_y_xyyyyz = cbuffer.data(pi_geom_11_off + 380 * ccomps * dcomps);

            auto g_y_y_y_xyyyzz = cbuffer.data(pi_geom_11_off + 381 * ccomps * dcomps);

            auto g_y_y_y_xyyzzz = cbuffer.data(pi_geom_11_off + 382 * ccomps * dcomps);

            auto g_y_y_y_xyzzzz = cbuffer.data(pi_geom_11_off + 383 * ccomps * dcomps);

            auto g_y_y_y_xzzzzz = cbuffer.data(pi_geom_11_off + 384 * ccomps * dcomps);

            auto g_y_y_y_yyyyyy = cbuffer.data(pi_geom_11_off + 385 * ccomps * dcomps);

            auto g_y_y_y_yyyyyz = cbuffer.data(pi_geom_11_off + 386 * ccomps * dcomps);

            auto g_y_y_y_yyyyzz = cbuffer.data(pi_geom_11_off + 387 * ccomps * dcomps);

            auto g_y_y_y_yyyzzz = cbuffer.data(pi_geom_11_off + 388 * ccomps * dcomps);

            auto g_y_y_y_yyzzzz = cbuffer.data(pi_geom_11_off + 389 * ccomps * dcomps);

            auto g_y_y_y_yzzzzz = cbuffer.data(pi_geom_11_off + 390 * ccomps * dcomps);

            auto g_y_y_y_zzzzzz = cbuffer.data(pi_geom_11_off + 391 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxxx, g_0_y_0_xxxxxy, g_0_y_0_xxxxxz, g_0_y_0_xxxxyy, g_0_y_0_xxxxyz, g_0_y_0_xxxxzz, g_0_y_0_xxxyyy, g_0_y_0_xxxyyz, g_0_y_0_xxxyzz, g_0_y_0_xxxzzz, g_0_y_0_xxyyyy, g_0_y_0_xxyyyz, g_0_y_0_xxyyzz, g_0_y_0_xxyzzz, g_0_y_0_xxzzzz, g_0_y_0_xyyyyy, g_0_y_0_xyyyyz, g_0_y_0_xyyyzz, g_0_y_0_xyyzzz, g_0_y_0_xyzzzz, g_0_y_0_xzzzzz, g_0_y_0_yyyyyy, g_0_y_0_yyyyyz, g_0_y_0_yyyyzz, g_0_y_0_yyyzzz, g_0_y_0_yyzzzz, g_0_y_0_yzzzzz, g_0_y_0_zzzzzz, g_y_0_0_xxxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxxz, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxzz, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzzz, g_y_y_0_xxxxxx, g_y_y_0_xxxxxxy, g_y_y_0_xxxxxy, g_y_y_0_xxxxxyy, g_y_y_0_xxxxxyz, g_y_y_0_xxxxxz, g_y_y_0_xxxxyy, g_y_y_0_xxxxyyy, g_y_y_0_xxxxyyz, g_y_y_0_xxxxyz, g_y_y_0_xxxxyzz, g_y_y_0_xxxxzz, g_y_y_0_xxxyyy, g_y_y_0_xxxyyyy, g_y_y_0_xxxyyyz, g_y_y_0_xxxyyz, g_y_y_0_xxxyyzz, g_y_y_0_xxxyzz, g_y_y_0_xxxyzzz, g_y_y_0_xxxzzz, g_y_y_0_xxyyyy, g_y_y_0_xxyyyyy, g_y_y_0_xxyyyyz, g_y_y_0_xxyyyz, g_y_y_0_xxyyyzz, g_y_y_0_xxyyzz, g_y_y_0_xxyyzzz, g_y_y_0_xxyzzz, g_y_y_0_xxyzzzz, g_y_y_0_xxzzzz, g_y_y_0_xyyyyy, g_y_y_0_xyyyyyy, g_y_y_0_xyyyyyz, g_y_y_0_xyyyyz, g_y_y_0_xyyyyzz, g_y_y_0_xyyyzz, g_y_y_0_xyyyzzz, g_y_y_0_xyyzzz, g_y_y_0_xyyzzzz, g_y_y_0_xyzzzz, g_y_y_0_xyzzzzz, g_y_y_0_xzzzzz, g_y_y_0_yyyyyy, g_y_y_0_yyyyyyy, g_y_y_0_yyyyyyz, g_y_y_0_yyyyyz, g_y_y_0_yyyyyzz, g_y_y_0_yyyyzz, g_y_y_0_yyyyzzz, g_y_y_0_yyyzzz, g_y_y_0_yyyzzzz, g_y_y_0_yyzzzz, g_y_y_0_yyzzzzz, g_y_y_0_yzzzzz, g_y_y_0_yzzzzzz, g_y_y_0_zzzzzz, g_y_y_y_xxxxxx, g_y_y_y_xxxxxy, g_y_y_y_xxxxxz, g_y_y_y_xxxxyy, g_y_y_y_xxxxyz, g_y_y_y_xxxxzz, g_y_y_y_xxxyyy, g_y_y_y_xxxyyz, g_y_y_y_xxxyzz, g_y_y_y_xxxzzz, g_y_y_y_xxyyyy, g_y_y_y_xxyyyz, g_y_y_y_xxyyzz, g_y_y_y_xxyzzz, g_y_y_y_xxzzzz, g_y_y_y_xyyyyy, g_y_y_y_xyyyyz, g_y_y_y_xyyyzz, g_y_y_y_xyyzzz, g_y_y_y_xyzzzz, g_y_y_y_xzzzzz, g_y_y_y_yyyyyy, g_y_y_y_yyyyyz, g_y_y_y_yyyyzz, g_y_y_y_yyyzzz, g_y_y_y_yyzzzz, g_y_y_y_yzzzzz, g_y_y_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_y_xxxxxx[k] = -g_0_y_0_xxxxxx[k] + g_y_0_0_xxxxxx[k] - g_y_y_0_xxxxxx[k] * ab_y + g_y_y_0_xxxxxxy[k];

                g_y_y_y_xxxxxy[k] = -g_0_y_0_xxxxxy[k] + g_y_0_0_xxxxxy[k] - g_y_y_0_xxxxxy[k] * ab_y + g_y_y_0_xxxxxyy[k];

                g_y_y_y_xxxxxz[k] = -g_0_y_0_xxxxxz[k] + g_y_0_0_xxxxxz[k] - g_y_y_0_xxxxxz[k] * ab_y + g_y_y_0_xxxxxyz[k];

                g_y_y_y_xxxxyy[k] = -g_0_y_0_xxxxyy[k] + g_y_0_0_xxxxyy[k] - g_y_y_0_xxxxyy[k] * ab_y + g_y_y_0_xxxxyyy[k];

                g_y_y_y_xxxxyz[k] = -g_0_y_0_xxxxyz[k] + g_y_0_0_xxxxyz[k] - g_y_y_0_xxxxyz[k] * ab_y + g_y_y_0_xxxxyyz[k];

                g_y_y_y_xxxxzz[k] = -g_0_y_0_xxxxzz[k] + g_y_0_0_xxxxzz[k] - g_y_y_0_xxxxzz[k] * ab_y + g_y_y_0_xxxxyzz[k];

                g_y_y_y_xxxyyy[k] = -g_0_y_0_xxxyyy[k] + g_y_0_0_xxxyyy[k] - g_y_y_0_xxxyyy[k] * ab_y + g_y_y_0_xxxyyyy[k];

                g_y_y_y_xxxyyz[k] = -g_0_y_0_xxxyyz[k] + g_y_0_0_xxxyyz[k] - g_y_y_0_xxxyyz[k] * ab_y + g_y_y_0_xxxyyyz[k];

                g_y_y_y_xxxyzz[k] = -g_0_y_0_xxxyzz[k] + g_y_0_0_xxxyzz[k] - g_y_y_0_xxxyzz[k] * ab_y + g_y_y_0_xxxyyzz[k];

                g_y_y_y_xxxzzz[k] = -g_0_y_0_xxxzzz[k] + g_y_0_0_xxxzzz[k] - g_y_y_0_xxxzzz[k] * ab_y + g_y_y_0_xxxyzzz[k];

                g_y_y_y_xxyyyy[k] = -g_0_y_0_xxyyyy[k] + g_y_0_0_xxyyyy[k] - g_y_y_0_xxyyyy[k] * ab_y + g_y_y_0_xxyyyyy[k];

                g_y_y_y_xxyyyz[k] = -g_0_y_0_xxyyyz[k] + g_y_0_0_xxyyyz[k] - g_y_y_0_xxyyyz[k] * ab_y + g_y_y_0_xxyyyyz[k];

                g_y_y_y_xxyyzz[k] = -g_0_y_0_xxyyzz[k] + g_y_0_0_xxyyzz[k] - g_y_y_0_xxyyzz[k] * ab_y + g_y_y_0_xxyyyzz[k];

                g_y_y_y_xxyzzz[k] = -g_0_y_0_xxyzzz[k] + g_y_0_0_xxyzzz[k] - g_y_y_0_xxyzzz[k] * ab_y + g_y_y_0_xxyyzzz[k];

                g_y_y_y_xxzzzz[k] = -g_0_y_0_xxzzzz[k] + g_y_0_0_xxzzzz[k] - g_y_y_0_xxzzzz[k] * ab_y + g_y_y_0_xxyzzzz[k];

                g_y_y_y_xyyyyy[k] = -g_0_y_0_xyyyyy[k] + g_y_0_0_xyyyyy[k] - g_y_y_0_xyyyyy[k] * ab_y + g_y_y_0_xyyyyyy[k];

                g_y_y_y_xyyyyz[k] = -g_0_y_0_xyyyyz[k] + g_y_0_0_xyyyyz[k] - g_y_y_0_xyyyyz[k] * ab_y + g_y_y_0_xyyyyyz[k];

                g_y_y_y_xyyyzz[k] = -g_0_y_0_xyyyzz[k] + g_y_0_0_xyyyzz[k] - g_y_y_0_xyyyzz[k] * ab_y + g_y_y_0_xyyyyzz[k];

                g_y_y_y_xyyzzz[k] = -g_0_y_0_xyyzzz[k] + g_y_0_0_xyyzzz[k] - g_y_y_0_xyyzzz[k] * ab_y + g_y_y_0_xyyyzzz[k];

                g_y_y_y_xyzzzz[k] = -g_0_y_0_xyzzzz[k] + g_y_0_0_xyzzzz[k] - g_y_y_0_xyzzzz[k] * ab_y + g_y_y_0_xyyzzzz[k];

                g_y_y_y_xzzzzz[k] = -g_0_y_0_xzzzzz[k] + g_y_0_0_xzzzzz[k] - g_y_y_0_xzzzzz[k] * ab_y + g_y_y_0_xyzzzzz[k];

                g_y_y_y_yyyyyy[k] = -g_0_y_0_yyyyyy[k] + g_y_0_0_yyyyyy[k] - g_y_y_0_yyyyyy[k] * ab_y + g_y_y_0_yyyyyyy[k];

                g_y_y_y_yyyyyz[k] = -g_0_y_0_yyyyyz[k] + g_y_0_0_yyyyyz[k] - g_y_y_0_yyyyyz[k] * ab_y + g_y_y_0_yyyyyyz[k];

                g_y_y_y_yyyyzz[k] = -g_0_y_0_yyyyzz[k] + g_y_0_0_yyyyzz[k] - g_y_y_0_yyyyzz[k] * ab_y + g_y_y_0_yyyyyzz[k];

                g_y_y_y_yyyzzz[k] = -g_0_y_0_yyyzzz[k] + g_y_0_0_yyyzzz[k] - g_y_y_0_yyyzzz[k] * ab_y + g_y_y_0_yyyyzzz[k];

                g_y_y_y_yyzzzz[k] = -g_0_y_0_yyzzzz[k] + g_y_0_0_yyzzzz[k] - g_y_y_0_yyzzzz[k] * ab_y + g_y_y_0_yyyzzzz[k];

                g_y_y_y_yzzzzz[k] = -g_0_y_0_yzzzzz[k] + g_y_0_0_yzzzzz[k] - g_y_y_0_yzzzzz[k] * ab_y + g_y_y_0_yyzzzzz[k];

                g_y_y_y_zzzzzz[k] = -g_0_y_0_zzzzzz[k] + g_y_0_0_zzzzzz[k] - g_y_y_0_zzzzzz[k] * ab_y + g_y_y_0_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_y_y_z_xxxxxx = cbuffer.data(pi_geom_11_off + 392 * ccomps * dcomps);

            auto g_y_y_z_xxxxxy = cbuffer.data(pi_geom_11_off + 393 * ccomps * dcomps);

            auto g_y_y_z_xxxxxz = cbuffer.data(pi_geom_11_off + 394 * ccomps * dcomps);

            auto g_y_y_z_xxxxyy = cbuffer.data(pi_geom_11_off + 395 * ccomps * dcomps);

            auto g_y_y_z_xxxxyz = cbuffer.data(pi_geom_11_off + 396 * ccomps * dcomps);

            auto g_y_y_z_xxxxzz = cbuffer.data(pi_geom_11_off + 397 * ccomps * dcomps);

            auto g_y_y_z_xxxyyy = cbuffer.data(pi_geom_11_off + 398 * ccomps * dcomps);

            auto g_y_y_z_xxxyyz = cbuffer.data(pi_geom_11_off + 399 * ccomps * dcomps);

            auto g_y_y_z_xxxyzz = cbuffer.data(pi_geom_11_off + 400 * ccomps * dcomps);

            auto g_y_y_z_xxxzzz = cbuffer.data(pi_geom_11_off + 401 * ccomps * dcomps);

            auto g_y_y_z_xxyyyy = cbuffer.data(pi_geom_11_off + 402 * ccomps * dcomps);

            auto g_y_y_z_xxyyyz = cbuffer.data(pi_geom_11_off + 403 * ccomps * dcomps);

            auto g_y_y_z_xxyyzz = cbuffer.data(pi_geom_11_off + 404 * ccomps * dcomps);

            auto g_y_y_z_xxyzzz = cbuffer.data(pi_geom_11_off + 405 * ccomps * dcomps);

            auto g_y_y_z_xxzzzz = cbuffer.data(pi_geom_11_off + 406 * ccomps * dcomps);

            auto g_y_y_z_xyyyyy = cbuffer.data(pi_geom_11_off + 407 * ccomps * dcomps);

            auto g_y_y_z_xyyyyz = cbuffer.data(pi_geom_11_off + 408 * ccomps * dcomps);

            auto g_y_y_z_xyyyzz = cbuffer.data(pi_geom_11_off + 409 * ccomps * dcomps);

            auto g_y_y_z_xyyzzz = cbuffer.data(pi_geom_11_off + 410 * ccomps * dcomps);

            auto g_y_y_z_xyzzzz = cbuffer.data(pi_geom_11_off + 411 * ccomps * dcomps);

            auto g_y_y_z_xzzzzz = cbuffer.data(pi_geom_11_off + 412 * ccomps * dcomps);

            auto g_y_y_z_yyyyyy = cbuffer.data(pi_geom_11_off + 413 * ccomps * dcomps);

            auto g_y_y_z_yyyyyz = cbuffer.data(pi_geom_11_off + 414 * ccomps * dcomps);

            auto g_y_y_z_yyyyzz = cbuffer.data(pi_geom_11_off + 415 * ccomps * dcomps);

            auto g_y_y_z_yyyzzz = cbuffer.data(pi_geom_11_off + 416 * ccomps * dcomps);

            auto g_y_y_z_yyzzzz = cbuffer.data(pi_geom_11_off + 417 * ccomps * dcomps);

            auto g_y_y_z_yzzzzz = cbuffer.data(pi_geom_11_off + 418 * ccomps * dcomps);

            auto g_y_y_z_zzzzzz = cbuffer.data(pi_geom_11_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_0_xxxxxx, g_y_y_0_xxxxxxz, g_y_y_0_xxxxxy, g_y_y_0_xxxxxyz, g_y_y_0_xxxxxz, g_y_y_0_xxxxxzz, g_y_y_0_xxxxyy, g_y_y_0_xxxxyyz, g_y_y_0_xxxxyz, g_y_y_0_xxxxyzz, g_y_y_0_xxxxzz, g_y_y_0_xxxxzzz, g_y_y_0_xxxyyy, g_y_y_0_xxxyyyz, g_y_y_0_xxxyyz, g_y_y_0_xxxyyzz, g_y_y_0_xxxyzz, g_y_y_0_xxxyzzz, g_y_y_0_xxxzzz, g_y_y_0_xxxzzzz, g_y_y_0_xxyyyy, g_y_y_0_xxyyyyz, g_y_y_0_xxyyyz, g_y_y_0_xxyyyzz, g_y_y_0_xxyyzz, g_y_y_0_xxyyzzz, g_y_y_0_xxyzzz, g_y_y_0_xxyzzzz, g_y_y_0_xxzzzz, g_y_y_0_xxzzzzz, g_y_y_0_xyyyyy, g_y_y_0_xyyyyyz, g_y_y_0_xyyyyz, g_y_y_0_xyyyyzz, g_y_y_0_xyyyzz, g_y_y_0_xyyyzzz, g_y_y_0_xyyzzz, g_y_y_0_xyyzzzz, g_y_y_0_xyzzzz, g_y_y_0_xyzzzzz, g_y_y_0_xzzzzz, g_y_y_0_xzzzzzz, g_y_y_0_yyyyyy, g_y_y_0_yyyyyyz, g_y_y_0_yyyyyz, g_y_y_0_yyyyyzz, g_y_y_0_yyyyzz, g_y_y_0_yyyyzzz, g_y_y_0_yyyzzz, g_y_y_0_yyyzzzz, g_y_y_0_yyzzzz, g_y_y_0_yyzzzzz, g_y_y_0_yzzzzz, g_y_y_0_yzzzzzz, g_y_y_0_zzzzzz, g_y_y_0_zzzzzzz, g_y_y_z_xxxxxx, g_y_y_z_xxxxxy, g_y_y_z_xxxxxz, g_y_y_z_xxxxyy, g_y_y_z_xxxxyz, g_y_y_z_xxxxzz, g_y_y_z_xxxyyy, g_y_y_z_xxxyyz, g_y_y_z_xxxyzz, g_y_y_z_xxxzzz, g_y_y_z_xxyyyy, g_y_y_z_xxyyyz, g_y_y_z_xxyyzz, g_y_y_z_xxyzzz, g_y_y_z_xxzzzz, g_y_y_z_xyyyyy, g_y_y_z_xyyyyz, g_y_y_z_xyyyzz, g_y_y_z_xyyzzz, g_y_y_z_xyzzzz, g_y_y_z_xzzzzz, g_y_y_z_yyyyyy, g_y_y_z_yyyyyz, g_y_y_z_yyyyzz, g_y_y_z_yyyzzz, g_y_y_z_yyzzzz, g_y_y_z_yzzzzz, g_y_y_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_z_xxxxxx[k] = -g_y_y_0_xxxxxx[k] * ab_z + g_y_y_0_xxxxxxz[k];

                g_y_y_z_xxxxxy[k] = -g_y_y_0_xxxxxy[k] * ab_z + g_y_y_0_xxxxxyz[k];

                g_y_y_z_xxxxxz[k] = -g_y_y_0_xxxxxz[k] * ab_z + g_y_y_0_xxxxxzz[k];

                g_y_y_z_xxxxyy[k] = -g_y_y_0_xxxxyy[k] * ab_z + g_y_y_0_xxxxyyz[k];

                g_y_y_z_xxxxyz[k] = -g_y_y_0_xxxxyz[k] * ab_z + g_y_y_0_xxxxyzz[k];

                g_y_y_z_xxxxzz[k] = -g_y_y_0_xxxxzz[k] * ab_z + g_y_y_0_xxxxzzz[k];

                g_y_y_z_xxxyyy[k] = -g_y_y_0_xxxyyy[k] * ab_z + g_y_y_0_xxxyyyz[k];

                g_y_y_z_xxxyyz[k] = -g_y_y_0_xxxyyz[k] * ab_z + g_y_y_0_xxxyyzz[k];

                g_y_y_z_xxxyzz[k] = -g_y_y_0_xxxyzz[k] * ab_z + g_y_y_0_xxxyzzz[k];

                g_y_y_z_xxxzzz[k] = -g_y_y_0_xxxzzz[k] * ab_z + g_y_y_0_xxxzzzz[k];

                g_y_y_z_xxyyyy[k] = -g_y_y_0_xxyyyy[k] * ab_z + g_y_y_0_xxyyyyz[k];

                g_y_y_z_xxyyyz[k] = -g_y_y_0_xxyyyz[k] * ab_z + g_y_y_0_xxyyyzz[k];

                g_y_y_z_xxyyzz[k] = -g_y_y_0_xxyyzz[k] * ab_z + g_y_y_0_xxyyzzz[k];

                g_y_y_z_xxyzzz[k] = -g_y_y_0_xxyzzz[k] * ab_z + g_y_y_0_xxyzzzz[k];

                g_y_y_z_xxzzzz[k] = -g_y_y_0_xxzzzz[k] * ab_z + g_y_y_0_xxzzzzz[k];

                g_y_y_z_xyyyyy[k] = -g_y_y_0_xyyyyy[k] * ab_z + g_y_y_0_xyyyyyz[k];

                g_y_y_z_xyyyyz[k] = -g_y_y_0_xyyyyz[k] * ab_z + g_y_y_0_xyyyyzz[k];

                g_y_y_z_xyyyzz[k] = -g_y_y_0_xyyyzz[k] * ab_z + g_y_y_0_xyyyzzz[k];

                g_y_y_z_xyyzzz[k] = -g_y_y_0_xyyzzz[k] * ab_z + g_y_y_0_xyyzzzz[k];

                g_y_y_z_xyzzzz[k] = -g_y_y_0_xyzzzz[k] * ab_z + g_y_y_0_xyzzzzz[k];

                g_y_y_z_xzzzzz[k] = -g_y_y_0_xzzzzz[k] * ab_z + g_y_y_0_xzzzzzz[k];

                g_y_y_z_yyyyyy[k] = -g_y_y_0_yyyyyy[k] * ab_z + g_y_y_0_yyyyyyz[k];

                g_y_y_z_yyyyyz[k] = -g_y_y_0_yyyyyz[k] * ab_z + g_y_y_0_yyyyyzz[k];

                g_y_y_z_yyyyzz[k] = -g_y_y_0_yyyyzz[k] * ab_z + g_y_y_0_yyyyzzz[k];

                g_y_y_z_yyyzzz[k] = -g_y_y_0_yyyzzz[k] * ab_z + g_y_y_0_yyyzzzz[k];

                g_y_y_z_yyzzzz[k] = -g_y_y_0_yyzzzz[k] * ab_z + g_y_y_0_yyzzzzz[k];

                g_y_y_z_yzzzzz[k] = -g_y_y_0_yzzzzz[k] * ab_z + g_y_y_0_yzzzzzz[k];

                g_y_y_z_zzzzzz[k] = -g_y_y_0_zzzzzz[k] * ab_z + g_y_y_0_zzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_y_z_x_xxxxxx = cbuffer.data(pi_geom_11_off + 420 * ccomps * dcomps);

            auto g_y_z_x_xxxxxy = cbuffer.data(pi_geom_11_off + 421 * ccomps * dcomps);

            auto g_y_z_x_xxxxxz = cbuffer.data(pi_geom_11_off + 422 * ccomps * dcomps);

            auto g_y_z_x_xxxxyy = cbuffer.data(pi_geom_11_off + 423 * ccomps * dcomps);

            auto g_y_z_x_xxxxyz = cbuffer.data(pi_geom_11_off + 424 * ccomps * dcomps);

            auto g_y_z_x_xxxxzz = cbuffer.data(pi_geom_11_off + 425 * ccomps * dcomps);

            auto g_y_z_x_xxxyyy = cbuffer.data(pi_geom_11_off + 426 * ccomps * dcomps);

            auto g_y_z_x_xxxyyz = cbuffer.data(pi_geom_11_off + 427 * ccomps * dcomps);

            auto g_y_z_x_xxxyzz = cbuffer.data(pi_geom_11_off + 428 * ccomps * dcomps);

            auto g_y_z_x_xxxzzz = cbuffer.data(pi_geom_11_off + 429 * ccomps * dcomps);

            auto g_y_z_x_xxyyyy = cbuffer.data(pi_geom_11_off + 430 * ccomps * dcomps);

            auto g_y_z_x_xxyyyz = cbuffer.data(pi_geom_11_off + 431 * ccomps * dcomps);

            auto g_y_z_x_xxyyzz = cbuffer.data(pi_geom_11_off + 432 * ccomps * dcomps);

            auto g_y_z_x_xxyzzz = cbuffer.data(pi_geom_11_off + 433 * ccomps * dcomps);

            auto g_y_z_x_xxzzzz = cbuffer.data(pi_geom_11_off + 434 * ccomps * dcomps);

            auto g_y_z_x_xyyyyy = cbuffer.data(pi_geom_11_off + 435 * ccomps * dcomps);

            auto g_y_z_x_xyyyyz = cbuffer.data(pi_geom_11_off + 436 * ccomps * dcomps);

            auto g_y_z_x_xyyyzz = cbuffer.data(pi_geom_11_off + 437 * ccomps * dcomps);

            auto g_y_z_x_xyyzzz = cbuffer.data(pi_geom_11_off + 438 * ccomps * dcomps);

            auto g_y_z_x_xyzzzz = cbuffer.data(pi_geom_11_off + 439 * ccomps * dcomps);

            auto g_y_z_x_xzzzzz = cbuffer.data(pi_geom_11_off + 440 * ccomps * dcomps);

            auto g_y_z_x_yyyyyy = cbuffer.data(pi_geom_11_off + 441 * ccomps * dcomps);

            auto g_y_z_x_yyyyyz = cbuffer.data(pi_geom_11_off + 442 * ccomps * dcomps);

            auto g_y_z_x_yyyyzz = cbuffer.data(pi_geom_11_off + 443 * ccomps * dcomps);

            auto g_y_z_x_yyyzzz = cbuffer.data(pi_geom_11_off + 444 * ccomps * dcomps);

            auto g_y_z_x_yyzzzz = cbuffer.data(pi_geom_11_off + 445 * ccomps * dcomps);

            auto g_y_z_x_yzzzzz = cbuffer.data(pi_geom_11_off + 446 * ccomps * dcomps);

            auto g_y_z_x_zzzzzz = cbuffer.data(pi_geom_11_off + 447 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_0_xxxxxx, g_y_z_0_xxxxxxx, g_y_z_0_xxxxxxy, g_y_z_0_xxxxxxz, g_y_z_0_xxxxxy, g_y_z_0_xxxxxyy, g_y_z_0_xxxxxyz, g_y_z_0_xxxxxz, g_y_z_0_xxxxxzz, g_y_z_0_xxxxyy, g_y_z_0_xxxxyyy, g_y_z_0_xxxxyyz, g_y_z_0_xxxxyz, g_y_z_0_xxxxyzz, g_y_z_0_xxxxzz, g_y_z_0_xxxxzzz, g_y_z_0_xxxyyy, g_y_z_0_xxxyyyy, g_y_z_0_xxxyyyz, g_y_z_0_xxxyyz, g_y_z_0_xxxyyzz, g_y_z_0_xxxyzz, g_y_z_0_xxxyzzz, g_y_z_0_xxxzzz, g_y_z_0_xxxzzzz, g_y_z_0_xxyyyy, g_y_z_0_xxyyyyy, g_y_z_0_xxyyyyz, g_y_z_0_xxyyyz, g_y_z_0_xxyyyzz, g_y_z_0_xxyyzz, g_y_z_0_xxyyzzz, g_y_z_0_xxyzzz, g_y_z_0_xxyzzzz, g_y_z_0_xxzzzz, g_y_z_0_xxzzzzz, g_y_z_0_xyyyyy, g_y_z_0_xyyyyyy, g_y_z_0_xyyyyyz, g_y_z_0_xyyyyz, g_y_z_0_xyyyyzz, g_y_z_0_xyyyzz, g_y_z_0_xyyyzzz, g_y_z_0_xyyzzz, g_y_z_0_xyyzzzz, g_y_z_0_xyzzzz, g_y_z_0_xyzzzzz, g_y_z_0_xzzzzz, g_y_z_0_xzzzzzz, g_y_z_0_yyyyyy, g_y_z_0_yyyyyz, g_y_z_0_yyyyzz, g_y_z_0_yyyzzz, g_y_z_0_yyzzzz, g_y_z_0_yzzzzz, g_y_z_0_zzzzzz, g_y_z_x_xxxxxx, g_y_z_x_xxxxxy, g_y_z_x_xxxxxz, g_y_z_x_xxxxyy, g_y_z_x_xxxxyz, g_y_z_x_xxxxzz, g_y_z_x_xxxyyy, g_y_z_x_xxxyyz, g_y_z_x_xxxyzz, g_y_z_x_xxxzzz, g_y_z_x_xxyyyy, g_y_z_x_xxyyyz, g_y_z_x_xxyyzz, g_y_z_x_xxyzzz, g_y_z_x_xxzzzz, g_y_z_x_xyyyyy, g_y_z_x_xyyyyz, g_y_z_x_xyyyzz, g_y_z_x_xyyzzz, g_y_z_x_xyzzzz, g_y_z_x_xzzzzz, g_y_z_x_yyyyyy, g_y_z_x_yyyyyz, g_y_z_x_yyyyzz, g_y_z_x_yyyzzz, g_y_z_x_yyzzzz, g_y_z_x_yzzzzz, g_y_z_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_x_xxxxxx[k] = -g_y_z_0_xxxxxx[k] * ab_x + g_y_z_0_xxxxxxx[k];

                g_y_z_x_xxxxxy[k] = -g_y_z_0_xxxxxy[k] * ab_x + g_y_z_0_xxxxxxy[k];

                g_y_z_x_xxxxxz[k] = -g_y_z_0_xxxxxz[k] * ab_x + g_y_z_0_xxxxxxz[k];

                g_y_z_x_xxxxyy[k] = -g_y_z_0_xxxxyy[k] * ab_x + g_y_z_0_xxxxxyy[k];

                g_y_z_x_xxxxyz[k] = -g_y_z_0_xxxxyz[k] * ab_x + g_y_z_0_xxxxxyz[k];

                g_y_z_x_xxxxzz[k] = -g_y_z_0_xxxxzz[k] * ab_x + g_y_z_0_xxxxxzz[k];

                g_y_z_x_xxxyyy[k] = -g_y_z_0_xxxyyy[k] * ab_x + g_y_z_0_xxxxyyy[k];

                g_y_z_x_xxxyyz[k] = -g_y_z_0_xxxyyz[k] * ab_x + g_y_z_0_xxxxyyz[k];

                g_y_z_x_xxxyzz[k] = -g_y_z_0_xxxyzz[k] * ab_x + g_y_z_0_xxxxyzz[k];

                g_y_z_x_xxxzzz[k] = -g_y_z_0_xxxzzz[k] * ab_x + g_y_z_0_xxxxzzz[k];

                g_y_z_x_xxyyyy[k] = -g_y_z_0_xxyyyy[k] * ab_x + g_y_z_0_xxxyyyy[k];

                g_y_z_x_xxyyyz[k] = -g_y_z_0_xxyyyz[k] * ab_x + g_y_z_0_xxxyyyz[k];

                g_y_z_x_xxyyzz[k] = -g_y_z_0_xxyyzz[k] * ab_x + g_y_z_0_xxxyyzz[k];

                g_y_z_x_xxyzzz[k] = -g_y_z_0_xxyzzz[k] * ab_x + g_y_z_0_xxxyzzz[k];

                g_y_z_x_xxzzzz[k] = -g_y_z_0_xxzzzz[k] * ab_x + g_y_z_0_xxxzzzz[k];

                g_y_z_x_xyyyyy[k] = -g_y_z_0_xyyyyy[k] * ab_x + g_y_z_0_xxyyyyy[k];

                g_y_z_x_xyyyyz[k] = -g_y_z_0_xyyyyz[k] * ab_x + g_y_z_0_xxyyyyz[k];

                g_y_z_x_xyyyzz[k] = -g_y_z_0_xyyyzz[k] * ab_x + g_y_z_0_xxyyyzz[k];

                g_y_z_x_xyyzzz[k] = -g_y_z_0_xyyzzz[k] * ab_x + g_y_z_0_xxyyzzz[k];

                g_y_z_x_xyzzzz[k] = -g_y_z_0_xyzzzz[k] * ab_x + g_y_z_0_xxyzzzz[k];

                g_y_z_x_xzzzzz[k] = -g_y_z_0_xzzzzz[k] * ab_x + g_y_z_0_xxzzzzz[k];

                g_y_z_x_yyyyyy[k] = -g_y_z_0_yyyyyy[k] * ab_x + g_y_z_0_xyyyyyy[k];

                g_y_z_x_yyyyyz[k] = -g_y_z_0_yyyyyz[k] * ab_x + g_y_z_0_xyyyyyz[k];

                g_y_z_x_yyyyzz[k] = -g_y_z_0_yyyyzz[k] * ab_x + g_y_z_0_xyyyyzz[k];

                g_y_z_x_yyyzzz[k] = -g_y_z_0_yyyzzz[k] * ab_x + g_y_z_0_xyyyzzz[k];

                g_y_z_x_yyzzzz[k] = -g_y_z_0_yyzzzz[k] * ab_x + g_y_z_0_xyyzzzz[k];

                g_y_z_x_yzzzzz[k] = -g_y_z_0_yzzzzz[k] * ab_x + g_y_z_0_xyzzzzz[k];

                g_y_z_x_zzzzzz[k] = -g_y_z_0_zzzzzz[k] * ab_x + g_y_z_0_xzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_y_z_y_xxxxxx = cbuffer.data(pi_geom_11_off + 448 * ccomps * dcomps);

            auto g_y_z_y_xxxxxy = cbuffer.data(pi_geom_11_off + 449 * ccomps * dcomps);

            auto g_y_z_y_xxxxxz = cbuffer.data(pi_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_z_y_xxxxyy = cbuffer.data(pi_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_z_y_xxxxyz = cbuffer.data(pi_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_z_y_xxxxzz = cbuffer.data(pi_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_z_y_xxxyyy = cbuffer.data(pi_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_z_y_xxxyyz = cbuffer.data(pi_geom_11_off + 455 * ccomps * dcomps);

            auto g_y_z_y_xxxyzz = cbuffer.data(pi_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_z_y_xxxzzz = cbuffer.data(pi_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_z_y_xxyyyy = cbuffer.data(pi_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_z_y_xxyyyz = cbuffer.data(pi_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_z_y_xxyyzz = cbuffer.data(pi_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_z_y_xxyzzz = cbuffer.data(pi_geom_11_off + 461 * ccomps * dcomps);

            auto g_y_z_y_xxzzzz = cbuffer.data(pi_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_z_y_xyyyyy = cbuffer.data(pi_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_z_y_xyyyyz = cbuffer.data(pi_geom_11_off + 464 * ccomps * dcomps);

            auto g_y_z_y_xyyyzz = cbuffer.data(pi_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_z_y_xyyzzz = cbuffer.data(pi_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_z_y_xyzzzz = cbuffer.data(pi_geom_11_off + 467 * ccomps * dcomps);

            auto g_y_z_y_xzzzzz = cbuffer.data(pi_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_z_y_yyyyyy = cbuffer.data(pi_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_z_y_yyyyyz = cbuffer.data(pi_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_z_y_yyyyzz = cbuffer.data(pi_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_z_y_yyyzzz = cbuffer.data(pi_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_z_y_yyzzzz = cbuffer.data(pi_geom_11_off + 473 * ccomps * dcomps);

            auto g_y_z_y_yzzzzz = cbuffer.data(pi_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_z_y_zzzzzz = cbuffer.data(pi_geom_11_off + 475 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxxx, g_0_z_0_xxxxxy, g_0_z_0_xxxxxz, g_0_z_0_xxxxyy, g_0_z_0_xxxxyz, g_0_z_0_xxxxzz, g_0_z_0_xxxyyy, g_0_z_0_xxxyyz, g_0_z_0_xxxyzz, g_0_z_0_xxxzzz, g_0_z_0_xxyyyy, g_0_z_0_xxyyyz, g_0_z_0_xxyyzz, g_0_z_0_xxyzzz, g_0_z_0_xxzzzz, g_0_z_0_xyyyyy, g_0_z_0_xyyyyz, g_0_z_0_xyyyzz, g_0_z_0_xyyzzz, g_0_z_0_xyzzzz, g_0_z_0_xzzzzz, g_0_z_0_yyyyyy, g_0_z_0_yyyyyz, g_0_z_0_yyyyzz, g_0_z_0_yyyzzz, g_0_z_0_yyzzzz, g_0_z_0_yzzzzz, g_0_z_0_zzzzzz, g_y_z_0_xxxxxx, g_y_z_0_xxxxxxy, g_y_z_0_xxxxxy, g_y_z_0_xxxxxyy, g_y_z_0_xxxxxyz, g_y_z_0_xxxxxz, g_y_z_0_xxxxyy, g_y_z_0_xxxxyyy, g_y_z_0_xxxxyyz, g_y_z_0_xxxxyz, g_y_z_0_xxxxyzz, g_y_z_0_xxxxzz, g_y_z_0_xxxyyy, g_y_z_0_xxxyyyy, g_y_z_0_xxxyyyz, g_y_z_0_xxxyyz, g_y_z_0_xxxyyzz, g_y_z_0_xxxyzz, g_y_z_0_xxxyzzz, g_y_z_0_xxxzzz, g_y_z_0_xxyyyy, g_y_z_0_xxyyyyy, g_y_z_0_xxyyyyz, g_y_z_0_xxyyyz, g_y_z_0_xxyyyzz, g_y_z_0_xxyyzz, g_y_z_0_xxyyzzz, g_y_z_0_xxyzzz, g_y_z_0_xxyzzzz, g_y_z_0_xxzzzz, g_y_z_0_xyyyyy, g_y_z_0_xyyyyyy, g_y_z_0_xyyyyyz, g_y_z_0_xyyyyz, g_y_z_0_xyyyyzz, g_y_z_0_xyyyzz, g_y_z_0_xyyyzzz, g_y_z_0_xyyzzz, g_y_z_0_xyyzzzz, g_y_z_0_xyzzzz, g_y_z_0_xyzzzzz, g_y_z_0_xzzzzz, g_y_z_0_yyyyyy, g_y_z_0_yyyyyyy, g_y_z_0_yyyyyyz, g_y_z_0_yyyyyz, g_y_z_0_yyyyyzz, g_y_z_0_yyyyzz, g_y_z_0_yyyyzzz, g_y_z_0_yyyzzz, g_y_z_0_yyyzzzz, g_y_z_0_yyzzzz, g_y_z_0_yyzzzzz, g_y_z_0_yzzzzz, g_y_z_0_yzzzzzz, g_y_z_0_zzzzzz, g_y_z_y_xxxxxx, g_y_z_y_xxxxxy, g_y_z_y_xxxxxz, g_y_z_y_xxxxyy, g_y_z_y_xxxxyz, g_y_z_y_xxxxzz, g_y_z_y_xxxyyy, g_y_z_y_xxxyyz, g_y_z_y_xxxyzz, g_y_z_y_xxxzzz, g_y_z_y_xxyyyy, g_y_z_y_xxyyyz, g_y_z_y_xxyyzz, g_y_z_y_xxyzzz, g_y_z_y_xxzzzz, g_y_z_y_xyyyyy, g_y_z_y_xyyyyz, g_y_z_y_xyyyzz, g_y_z_y_xyyzzz, g_y_z_y_xyzzzz, g_y_z_y_xzzzzz, g_y_z_y_yyyyyy, g_y_z_y_yyyyyz, g_y_z_y_yyyyzz, g_y_z_y_yyyzzz, g_y_z_y_yyzzzz, g_y_z_y_yzzzzz, g_y_z_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_y_xxxxxx[k] = -g_0_z_0_xxxxxx[k] - g_y_z_0_xxxxxx[k] * ab_y + g_y_z_0_xxxxxxy[k];

                g_y_z_y_xxxxxy[k] = -g_0_z_0_xxxxxy[k] - g_y_z_0_xxxxxy[k] * ab_y + g_y_z_0_xxxxxyy[k];

                g_y_z_y_xxxxxz[k] = -g_0_z_0_xxxxxz[k] - g_y_z_0_xxxxxz[k] * ab_y + g_y_z_0_xxxxxyz[k];

                g_y_z_y_xxxxyy[k] = -g_0_z_0_xxxxyy[k] - g_y_z_0_xxxxyy[k] * ab_y + g_y_z_0_xxxxyyy[k];

                g_y_z_y_xxxxyz[k] = -g_0_z_0_xxxxyz[k] - g_y_z_0_xxxxyz[k] * ab_y + g_y_z_0_xxxxyyz[k];

                g_y_z_y_xxxxzz[k] = -g_0_z_0_xxxxzz[k] - g_y_z_0_xxxxzz[k] * ab_y + g_y_z_0_xxxxyzz[k];

                g_y_z_y_xxxyyy[k] = -g_0_z_0_xxxyyy[k] - g_y_z_0_xxxyyy[k] * ab_y + g_y_z_0_xxxyyyy[k];

                g_y_z_y_xxxyyz[k] = -g_0_z_0_xxxyyz[k] - g_y_z_0_xxxyyz[k] * ab_y + g_y_z_0_xxxyyyz[k];

                g_y_z_y_xxxyzz[k] = -g_0_z_0_xxxyzz[k] - g_y_z_0_xxxyzz[k] * ab_y + g_y_z_0_xxxyyzz[k];

                g_y_z_y_xxxzzz[k] = -g_0_z_0_xxxzzz[k] - g_y_z_0_xxxzzz[k] * ab_y + g_y_z_0_xxxyzzz[k];

                g_y_z_y_xxyyyy[k] = -g_0_z_0_xxyyyy[k] - g_y_z_0_xxyyyy[k] * ab_y + g_y_z_0_xxyyyyy[k];

                g_y_z_y_xxyyyz[k] = -g_0_z_0_xxyyyz[k] - g_y_z_0_xxyyyz[k] * ab_y + g_y_z_0_xxyyyyz[k];

                g_y_z_y_xxyyzz[k] = -g_0_z_0_xxyyzz[k] - g_y_z_0_xxyyzz[k] * ab_y + g_y_z_0_xxyyyzz[k];

                g_y_z_y_xxyzzz[k] = -g_0_z_0_xxyzzz[k] - g_y_z_0_xxyzzz[k] * ab_y + g_y_z_0_xxyyzzz[k];

                g_y_z_y_xxzzzz[k] = -g_0_z_0_xxzzzz[k] - g_y_z_0_xxzzzz[k] * ab_y + g_y_z_0_xxyzzzz[k];

                g_y_z_y_xyyyyy[k] = -g_0_z_0_xyyyyy[k] - g_y_z_0_xyyyyy[k] * ab_y + g_y_z_0_xyyyyyy[k];

                g_y_z_y_xyyyyz[k] = -g_0_z_0_xyyyyz[k] - g_y_z_0_xyyyyz[k] * ab_y + g_y_z_0_xyyyyyz[k];

                g_y_z_y_xyyyzz[k] = -g_0_z_0_xyyyzz[k] - g_y_z_0_xyyyzz[k] * ab_y + g_y_z_0_xyyyyzz[k];

                g_y_z_y_xyyzzz[k] = -g_0_z_0_xyyzzz[k] - g_y_z_0_xyyzzz[k] * ab_y + g_y_z_0_xyyyzzz[k];

                g_y_z_y_xyzzzz[k] = -g_0_z_0_xyzzzz[k] - g_y_z_0_xyzzzz[k] * ab_y + g_y_z_0_xyyzzzz[k];

                g_y_z_y_xzzzzz[k] = -g_0_z_0_xzzzzz[k] - g_y_z_0_xzzzzz[k] * ab_y + g_y_z_0_xyzzzzz[k];

                g_y_z_y_yyyyyy[k] = -g_0_z_0_yyyyyy[k] - g_y_z_0_yyyyyy[k] * ab_y + g_y_z_0_yyyyyyy[k];

                g_y_z_y_yyyyyz[k] = -g_0_z_0_yyyyyz[k] - g_y_z_0_yyyyyz[k] * ab_y + g_y_z_0_yyyyyyz[k];

                g_y_z_y_yyyyzz[k] = -g_0_z_0_yyyyzz[k] - g_y_z_0_yyyyzz[k] * ab_y + g_y_z_0_yyyyyzz[k];

                g_y_z_y_yyyzzz[k] = -g_0_z_0_yyyzzz[k] - g_y_z_0_yyyzzz[k] * ab_y + g_y_z_0_yyyyzzz[k];

                g_y_z_y_yyzzzz[k] = -g_0_z_0_yyzzzz[k] - g_y_z_0_yyzzzz[k] * ab_y + g_y_z_0_yyyzzzz[k];

                g_y_z_y_yzzzzz[k] = -g_0_z_0_yzzzzz[k] - g_y_z_0_yzzzzz[k] * ab_y + g_y_z_0_yyzzzzz[k];

                g_y_z_y_zzzzzz[k] = -g_0_z_0_zzzzzz[k] - g_y_z_0_zzzzzz[k] * ab_y + g_y_z_0_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_y_z_z_xxxxxx = cbuffer.data(pi_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_z_z_xxxxxy = cbuffer.data(pi_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_z_z_xxxxxz = cbuffer.data(pi_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_z_z_xxxxyy = cbuffer.data(pi_geom_11_off + 479 * ccomps * dcomps);

            auto g_y_z_z_xxxxyz = cbuffer.data(pi_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_z_z_xxxxzz = cbuffer.data(pi_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_z_z_xxxyyy = cbuffer.data(pi_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_z_z_xxxyyz = cbuffer.data(pi_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_z_z_xxxyzz = cbuffer.data(pi_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_z_z_xxxzzz = cbuffer.data(pi_geom_11_off + 485 * ccomps * dcomps);

            auto g_y_z_z_xxyyyy = cbuffer.data(pi_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_z_z_xxyyyz = cbuffer.data(pi_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_z_z_xxyyzz = cbuffer.data(pi_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_z_z_xxyzzz = cbuffer.data(pi_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_z_z_xxzzzz = cbuffer.data(pi_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_z_z_xyyyyy = cbuffer.data(pi_geom_11_off + 491 * ccomps * dcomps);

            auto g_y_z_z_xyyyyz = cbuffer.data(pi_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_z_z_xyyyzz = cbuffer.data(pi_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_z_z_xyyzzz = cbuffer.data(pi_geom_11_off + 494 * ccomps * dcomps);

            auto g_y_z_z_xyzzzz = cbuffer.data(pi_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_z_z_xzzzzz = cbuffer.data(pi_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_z_z_yyyyyy = cbuffer.data(pi_geom_11_off + 497 * ccomps * dcomps);

            auto g_y_z_z_yyyyyz = cbuffer.data(pi_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_z_z_yyyyzz = cbuffer.data(pi_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_z_z_yyyzzz = cbuffer.data(pi_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_z_z_yyzzzz = cbuffer.data(pi_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_z_z_yzzzzz = cbuffer.data(pi_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_z_z_zzzzzz = cbuffer.data(pi_geom_11_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_0_xxxxxx, g_y_0_0_xxxxxy, g_y_0_0_xxxxxz, g_y_0_0_xxxxyy, g_y_0_0_xxxxyz, g_y_0_0_xxxxzz, g_y_0_0_xxxyyy, g_y_0_0_xxxyyz, g_y_0_0_xxxyzz, g_y_0_0_xxxzzz, g_y_0_0_xxyyyy, g_y_0_0_xxyyyz, g_y_0_0_xxyyzz, g_y_0_0_xxyzzz, g_y_0_0_xxzzzz, g_y_0_0_xyyyyy, g_y_0_0_xyyyyz, g_y_0_0_xyyyzz, g_y_0_0_xyyzzz, g_y_0_0_xyzzzz, g_y_0_0_xzzzzz, g_y_0_0_yyyyyy, g_y_0_0_yyyyyz, g_y_0_0_yyyyzz, g_y_0_0_yyyzzz, g_y_0_0_yyzzzz, g_y_0_0_yzzzzz, g_y_0_0_zzzzzz, g_y_z_0_xxxxxx, g_y_z_0_xxxxxxz, g_y_z_0_xxxxxy, g_y_z_0_xxxxxyz, g_y_z_0_xxxxxz, g_y_z_0_xxxxxzz, g_y_z_0_xxxxyy, g_y_z_0_xxxxyyz, g_y_z_0_xxxxyz, g_y_z_0_xxxxyzz, g_y_z_0_xxxxzz, g_y_z_0_xxxxzzz, g_y_z_0_xxxyyy, g_y_z_0_xxxyyyz, g_y_z_0_xxxyyz, g_y_z_0_xxxyyzz, g_y_z_0_xxxyzz, g_y_z_0_xxxyzzz, g_y_z_0_xxxzzz, g_y_z_0_xxxzzzz, g_y_z_0_xxyyyy, g_y_z_0_xxyyyyz, g_y_z_0_xxyyyz, g_y_z_0_xxyyyzz, g_y_z_0_xxyyzz, g_y_z_0_xxyyzzz, g_y_z_0_xxyzzz, g_y_z_0_xxyzzzz, g_y_z_0_xxzzzz, g_y_z_0_xxzzzzz, g_y_z_0_xyyyyy, g_y_z_0_xyyyyyz, g_y_z_0_xyyyyz, g_y_z_0_xyyyyzz, g_y_z_0_xyyyzz, g_y_z_0_xyyyzzz, g_y_z_0_xyyzzz, g_y_z_0_xyyzzzz, g_y_z_0_xyzzzz, g_y_z_0_xyzzzzz, g_y_z_0_xzzzzz, g_y_z_0_xzzzzzz, g_y_z_0_yyyyyy, g_y_z_0_yyyyyyz, g_y_z_0_yyyyyz, g_y_z_0_yyyyyzz, g_y_z_0_yyyyzz, g_y_z_0_yyyyzzz, g_y_z_0_yyyzzz, g_y_z_0_yyyzzzz, g_y_z_0_yyzzzz, g_y_z_0_yyzzzzz, g_y_z_0_yzzzzz, g_y_z_0_yzzzzzz, g_y_z_0_zzzzzz, g_y_z_0_zzzzzzz, g_y_z_z_xxxxxx, g_y_z_z_xxxxxy, g_y_z_z_xxxxxz, g_y_z_z_xxxxyy, g_y_z_z_xxxxyz, g_y_z_z_xxxxzz, g_y_z_z_xxxyyy, g_y_z_z_xxxyyz, g_y_z_z_xxxyzz, g_y_z_z_xxxzzz, g_y_z_z_xxyyyy, g_y_z_z_xxyyyz, g_y_z_z_xxyyzz, g_y_z_z_xxyzzz, g_y_z_z_xxzzzz, g_y_z_z_xyyyyy, g_y_z_z_xyyyyz, g_y_z_z_xyyyzz, g_y_z_z_xyyzzz, g_y_z_z_xyzzzz, g_y_z_z_xzzzzz, g_y_z_z_yyyyyy, g_y_z_z_yyyyyz, g_y_z_z_yyyyzz, g_y_z_z_yyyzzz, g_y_z_z_yyzzzz, g_y_z_z_yzzzzz, g_y_z_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_z_xxxxxx[k] = g_y_0_0_xxxxxx[k] - g_y_z_0_xxxxxx[k] * ab_z + g_y_z_0_xxxxxxz[k];

                g_y_z_z_xxxxxy[k] = g_y_0_0_xxxxxy[k] - g_y_z_0_xxxxxy[k] * ab_z + g_y_z_0_xxxxxyz[k];

                g_y_z_z_xxxxxz[k] = g_y_0_0_xxxxxz[k] - g_y_z_0_xxxxxz[k] * ab_z + g_y_z_0_xxxxxzz[k];

                g_y_z_z_xxxxyy[k] = g_y_0_0_xxxxyy[k] - g_y_z_0_xxxxyy[k] * ab_z + g_y_z_0_xxxxyyz[k];

                g_y_z_z_xxxxyz[k] = g_y_0_0_xxxxyz[k] - g_y_z_0_xxxxyz[k] * ab_z + g_y_z_0_xxxxyzz[k];

                g_y_z_z_xxxxzz[k] = g_y_0_0_xxxxzz[k] - g_y_z_0_xxxxzz[k] * ab_z + g_y_z_0_xxxxzzz[k];

                g_y_z_z_xxxyyy[k] = g_y_0_0_xxxyyy[k] - g_y_z_0_xxxyyy[k] * ab_z + g_y_z_0_xxxyyyz[k];

                g_y_z_z_xxxyyz[k] = g_y_0_0_xxxyyz[k] - g_y_z_0_xxxyyz[k] * ab_z + g_y_z_0_xxxyyzz[k];

                g_y_z_z_xxxyzz[k] = g_y_0_0_xxxyzz[k] - g_y_z_0_xxxyzz[k] * ab_z + g_y_z_0_xxxyzzz[k];

                g_y_z_z_xxxzzz[k] = g_y_0_0_xxxzzz[k] - g_y_z_0_xxxzzz[k] * ab_z + g_y_z_0_xxxzzzz[k];

                g_y_z_z_xxyyyy[k] = g_y_0_0_xxyyyy[k] - g_y_z_0_xxyyyy[k] * ab_z + g_y_z_0_xxyyyyz[k];

                g_y_z_z_xxyyyz[k] = g_y_0_0_xxyyyz[k] - g_y_z_0_xxyyyz[k] * ab_z + g_y_z_0_xxyyyzz[k];

                g_y_z_z_xxyyzz[k] = g_y_0_0_xxyyzz[k] - g_y_z_0_xxyyzz[k] * ab_z + g_y_z_0_xxyyzzz[k];

                g_y_z_z_xxyzzz[k] = g_y_0_0_xxyzzz[k] - g_y_z_0_xxyzzz[k] * ab_z + g_y_z_0_xxyzzzz[k];

                g_y_z_z_xxzzzz[k] = g_y_0_0_xxzzzz[k] - g_y_z_0_xxzzzz[k] * ab_z + g_y_z_0_xxzzzzz[k];

                g_y_z_z_xyyyyy[k] = g_y_0_0_xyyyyy[k] - g_y_z_0_xyyyyy[k] * ab_z + g_y_z_0_xyyyyyz[k];

                g_y_z_z_xyyyyz[k] = g_y_0_0_xyyyyz[k] - g_y_z_0_xyyyyz[k] * ab_z + g_y_z_0_xyyyyzz[k];

                g_y_z_z_xyyyzz[k] = g_y_0_0_xyyyzz[k] - g_y_z_0_xyyyzz[k] * ab_z + g_y_z_0_xyyyzzz[k];

                g_y_z_z_xyyzzz[k] = g_y_0_0_xyyzzz[k] - g_y_z_0_xyyzzz[k] * ab_z + g_y_z_0_xyyzzzz[k];

                g_y_z_z_xyzzzz[k] = g_y_0_0_xyzzzz[k] - g_y_z_0_xyzzzz[k] * ab_z + g_y_z_0_xyzzzzz[k];

                g_y_z_z_xzzzzz[k] = g_y_0_0_xzzzzz[k] - g_y_z_0_xzzzzz[k] * ab_z + g_y_z_0_xzzzzzz[k];

                g_y_z_z_yyyyyy[k] = g_y_0_0_yyyyyy[k] - g_y_z_0_yyyyyy[k] * ab_z + g_y_z_0_yyyyyyz[k];

                g_y_z_z_yyyyyz[k] = g_y_0_0_yyyyyz[k] - g_y_z_0_yyyyyz[k] * ab_z + g_y_z_0_yyyyyzz[k];

                g_y_z_z_yyyyzz[k] = g_y_0_0_yyyyzz[k] - g_y_z_0_yyyyzz[k] * ab_z + g_y_z_0_yyyyzzz[k];

                g_y_z_z_yyyzzz[k] = g_y_0_0_yyyzzz[k] - g_y_z_0_yyyzzz[k] * ab_z + g_y_z_0_yyyzzzz[k];

                g_y_z_z_yyzzzz[k] = g_y_0_0_yyzzzz[k] - g_y_z_0_yyzzzz[k] * ab_z + g_y_z_0_yyzzzzz[k];

                g_y_z_z_yzzzzz[k] = g_y_0_0_yzzzzz[k] - g_y_z_0_yzzzzz[k] * ab_z + g_y_z_0_yzzzzzz[k];

                g_y_z_z_zzzzzz[k] = g_y_0_0_zzzzzz[k] - g_y_z_0_zzzzzz[k] * ab_z + g_y_z_0_zzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_z_x_x_xxxxxx = cbuffer.data(pi_geom_11_off + 504 * ccomps * dcomps);

            auto g_z_x_x_xxxxxy = cbuffer.data(pi_geom_11_off + 505 * ccomps * dcomps);

            auto g_z_x_x_xxxxxz = cbuffer.data(pi_geom_11_off + 506 * ccomps * dcomps);

            auto g_z_x_x_xxxxyy = cbuffer.data(pi_geom_11_off + 507 * ccomps * dcomps);

            auto g_z_x_x_xxxxyz = cbuffer.data(pi_geom_11_off + 508 * ccomps * dcomps);

            auto g_z_x_x_xxxxzz = cbuffer.data(pi_geom_11_off + 509 * ccomps * dcomps);

            auto g_z_x_x_xxxyyy = cbuffer.data(pi_geom_11_off + 510 * ccomps * dcomps);

            auto g_z_x_x_xxxyyz = cbuffer.data(pi_geom_11_off + 511 * ccomps * dcomps);

            auto g_z_x_x_xxxyzz = cbuffer.data(pi_geom_11_off + 512 * ccomps * dcomps);

            auto g_z_x_x_xxxzzz = cbuffer.data(pi_geom_11_off + 513 * ccomps * dcomps);

            auto g_z_x_x_xxyyyy = cbuffer.data(pi_geom_11_off + 514 * ccomps * dcomps);

            auto g_z_x_x_xxyyyz = cbuffer.data(pi_geom_11_off + 515 * ccomps * dcomps);

            auto g_z_x_x_xxyyzz = cbuffer.data(pi_geom_11_off + 516 * ccomps * dcomps);

            auto g_z_x_x_xxyzzz = cbuffer.data(pi_geom_11_off + 517 * ccomps * dcomps);

            auto g_z_x_x_xxzzzz = cbuffer.data(pi_geom_11_off + 518 * ccomps * dcomps);

            auto g_z_x_x_xyyyyy = cbuffer.data(pi_geom_11_off + 519 * ccomps * dcomps);

            auto g_z_x_x_xyyyyz = cbuffer.data(pi_geom_11_off + 520 * ccomps * dcomps);

            auto g_z_x_x_xyyyzz = cbuffer.data(pi_geom_11_off + 521 * ccomps * dcomps);

            auto g_z_x_x_xyyzzz = cbuffer.data(pi_geom_11_off + 522 * ccomps * dcomps);

            auto g_z_x_x_xyzzzz = cbuffer.data(pi_geom_11_off + 523 * ccomps * dcomps);

            auto g_z_x_x_xzzzzz = cbuffer.data(pi_geom_11_off + 524 * ccomps * dcomps);

            auto g_z_x_x_yyyyyy = cbuffer.data(pi_geom_11_off + 525 * ccomps * dcomps);

            auto g_z_x_x_yyyyyz = cbuffer.data(pi_geom_11_off + 526 * ccomps * dcomps);

            auto g_z_x_x_yyyyzz = cbuffer.data(pi_geom_11_off + 527 * ccomps * dcomps);

            auto g_z_x_x_yyyzzz = cbuffer.data(pi_geom_11_off + 528 * ccomps * dcomps);

            auto g_z_x_x_yyzzzz = cbuffer.data(pi_geom_11_off + 529 * ccomps * dcomps);

            auto g_z_x_x_yzzzzz = cbuffer.data(pi_geom_11_off + 530 * ccomps * dcomps);

            auto g_z_x_x_zzzzzz = cbuffer.data(pi_geom_11_off + 531 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxxz, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxzz, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzzz, g_z_x_0_xxxxxx, g_z_x_0_xxxxxxx, g_z_x_0_xxxxxxy, g_z_x_0_xxxxxxz, g_z_x_0_xxxxxy, g_z_x_0_xxxxxyy, g_z_x_0_xxxxxyz, g_z_x_0_xxxxxz, g_z_x_0_xxxxxzz, g_z_x_0_xxxxyy, g_z_x_0_xxxxyyy, g_z_x_0_xxxxyyz, g_z_x_0_xxxxyz, g_z_x_0_xxxxyzz, g_z_x_0_xxxxzz, g_z_x_0_xxxxzzz, g_z_x_0_xxxyyy, g_z_x_0_xxxyyyy, g_z_x_0_xxxyyyz, g_z_x_0_xxxyyz, g_z_x_0_xxxyyzz, g_z_x_0_xxxyzz, g_z_x_0_xxxyzzz, g_z_x_0_xxxzzz, g_z_x_0_xxxzzzz, g_z_x_0_xxyyyy, g_z_x_0_xxyyyyy, g_z_x_0_xxyyyyz, g_z_x_0_xxyyyz, g_z_x_0_xxyyyzz, g_z_x_0_xxyyzz, g_z_x_0_xxyyzzz, g_z_x_0_xxyzzz, g_z_x_0_xxyzzzz, g_z_x_0_xxzzzz, g_z_x_0_xxzzzzz, g_z_x_0_xyyyyy, g_z_x_0_xyyyyyy, g_z_x_0_xyyyyyz, g_z_x_0_xyyyyz, g_z_x_0_xyyyyzz, g_z_x_0_xyyyzz, g_z_x_0_xyyyzzz, g_z_x_0_xyyzzz, g_z_x_0_xyyzzzz, g_z_x_0_xyzzzz, g_z_x_0_xyzzzzz, g_z_x_0_xzzzzz, g_z_x_0_xzzzzzz, g_z_x_0_yyyyyy, g_z_x_0_yyyyyz, g_z_x_0_yyyyzz, g_z_x_0_yyyzzz, g_z_x_0_yyzzzz, g_z_x_0_yzzzzz, g_z_x_0_zzzzzz, g_z_x_x_xxxxxx, g_z_x_x_xxxxxy, g_z_x_x_xxxxxz, g_z_x_x_xxxxyy, g_z_x_x_xxxxyz, g_z_x_x_xxxxzz, g_z_x_x_xxxyyy, g_z_x_x_xxxyyz, g_z_x_x_xxxyzz, g_z_x_x_xxxzzz, g_z_x_x_xxyyyy, g_z_x_x_xxyyyz, g_z_x_x_xxyyzz, g_z_x_x_xxyzzz, g_z_x_x_xxzzzz, g_z_x_x_xyyyyy, g_z_x_x_xyyyyz, g_z_x_x_xyyyzz, g_z_x_x_xyyzzz, g_z_x_x_xyzzzz, g_z_x_x_xzzzzz, g_z_x_x_yyyyyy, g_z_x_x_yyyyyz, g_z_x_x_yyyyzz, g_z_x_x_yyyzzz, g_z_x_x_yyzzzz, g_z_x_x_yzzzzz, g_z_x_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_x_xxxxxx[k] = g_z_0_0_xxxxxx[k] - g_z_x_0_xxxxxx[k] * ab_x + g_z_x_0_xxxxxxx[k];

                g_z_x_x_xxxxxy[k] = g_z_0_0_xxxxxy[k] - g_z_x_0_xxxxxy[k] * ab_x + g_z_x_0_xxxxxxy[k];

                g_z_x_x_xxxxxz[k] = g_z_0_0_xxxxxz[k] - g_z_x_0_xxxxxz[k] * ab_x + g_z_x_0_xxxxxxz[k];

                g_z_x_x_xxxxyy[k] = g_z_0_0_xxxxyy[k] - g_z_x_0_xxxxyy[k] * ab_x + g_z_x_0_xxxxxyy[k];

                g_z_x_x_xxxxyz[k] = g_z_0_0_xxxxyz[k] - g_z_x_0_xxxxyz[k] * ab_x + g_z_x_0_xxxxxyz[k];

                g_z_x_x_xxxxzz[k] = g_z_0_0_xxxxzz[k] - g_z_x_0_xxxxzz[k] * ab_x + g_z_x_0_xxxxxzz[k];

                g_z_x_x_xxxyyy[k] = g_z_0_0_xxxyyy[k] - g_z_x_0_xxxyyy[k] * ab_x + g_z_x_0_xxxxyyy[k];

                g_z_x_x_xxxyyz[k] = g_z_0_0_xxxyyz[k] - g_z_x_0_xxxyyz[k] * ab_x + g_z_x_0_xxxxyyz[k];

                g_z_x_x_xxxyzz[k] = g_z_0_0_xxxyzz[k] - g_z_x_0_xxxyzz[k] * ab_x + g_z_x_0_xxxxyzz[k];

                g_z_x_x_xxxzzz[k] = g_z_0_0_xxxzzz[k] - g_z_x_0_xxxzzz[k] * ab_x + g_z_x_0_xxxxzzz[k];

                g_z_x_x_xxyyyy[k] = g_z_0_0_xxyyyy[k] - g_z_x_0_xxyyyy[k] * ab_x + g_z_x_0_xxxyyyy[k];

                g_z_x_x_xxyyyz[k] = g_z_0_0_xxyyyz[k] - g_z_x_0_xxyyyz[k] * ab_x + g_z_x_0_xxxyyyz[k];

                g_z_x_x_xxyyzz[k] = g_z_0_0_xxyyzz[k] - g_z_x_0_xxyyzz[k] * ab_x + g_z_x_0_xxxyyzz[k];

                g_z_x_x_xxyzzz[k] = g_z_0_0_xxyzzz[k] - g_z_x_0_xxyzzz[k] * ab_x + g_z_x_0_xxxyzzz[k];

                g_z_x_x_xxzzzz[k] = g_z_0_0_xxzzzz[k] - g_z_x_0_xxzzzz[k] * ab_x + g_z_x_0_xxxzzzz[k];

                g_z_x_x_xyyyyy[k] = g_z_0_0_xyyyyy[k] - g_z_x_0_xyyyyy[k] * ab_x + g_z_x_0_xxyyyyy[k];

                g_z_x_x_xyyyyz[k] = g_z_0_0_xyyyyz[k] - g_z_x_0_xyyyyz[k] * ab_x + g_z_x_0_xxyyyyz[k];

                g_z_x_x_xyyyzz[k] = g_z_0_0_xyyyzz[k] - g_z_x_0_xyyyzz[k] * ab_x + g_z_x_0_xxyyyzz[k];

                g_z_x_x_xyyzzz[k] = g_z_0_0_xyyzzz[k] - g_z_x_0_xyyzzz[k] * ab_x + g_z_x_0_xxyyzzz[k];

                g_z_x_x_xyzzzz[k] = g_z_0_0_xyzzzz[k] - g_z_x_0_xyzzzz[k] * ab_x + g_z_x_0_xxyzzzz[k];

                g_z_x_x_xzzzzz[k] = g_z_0_0_xzzzzz[k] - g_z_x_0_xzzzzz[k] * ab_x + g_z_x_0_xxzzzzz[k];

                g_z_x_x_yyyyyy[k] = g_z_0_0_yyyyyy[k] - g_z_x_0_yyyyyy[k] * ab_x + g_z_x_0_xyyyyyy[k];

                g_z_x_x_yyyyyz[k] = g_z_0_0_yyyyyz[k] - g_z_x_0_yyyyyz[k] * ab_x + g_z_x_0_xyyyyyz[k];

                g_z_x_x_yyyyzz[k] = g_z_0_0_yyyyzz[k] - g_z_x_0_yyyyzz[k] * ab_x + g_z_x_0_xyyyyzz[k];

                g_z_x_x_yyyzzz[k] = g_z_0_0_yyyzzz[k] - g_z_x_0_yyyzzz[k] * ab_x + g_z_x_0_xyyyzzz[k];

                g_z_x_x_yyzzzz[k] = g_z_0_0_yyzzzz[k] - g_z_x_0_yyzzzz[k] * ab_x + g_z_x_0_xyyzzzz[k];

                g_z_x_x_yzzzzz[k] = g_z_0_0_yzzzzz[k] - g_z_x_0_yzzzzz[k] * ab_x + g_z_x_0_xyzzzzz[k];

                g_z_x_x_zzzzzz[k] = g_z_0_0_zzzzzz[k] - g_z_x_0_zzzzzz[k] * ab_x + g_z_x_0_xzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_z_x_y_xxxxxx = cbuffer.data(pi_geom_11_off + 532 * ccomps * dcomps);

            auto g_z_x_y_xxxxxy = cbuffer.data(pi_geom_11_off + 533 * ccomps * dcomps);

            auto g_z_x_y_xxxxxz = cbuffer.data(pi_geom_11_off + 534 * ccomps * dcomps);

            auto g_z_x_y_xxxxyy = cbuffer.data(pi_geom_11_off + 535 * ccomps * dcomps);

            auto g_z_x_y_xxxxyz = cbuffer.data(pi_geom_11_off + 536 * ccomps * dcomps);

            auto g_z_x_y_xxxxzz = cbuffer.data(pi_geom_11_off + 537 * ccomps * dcomps);

            auto g_z_x_y_xxxyyy = cbuffer.data(pi_geom_11_off + 538 * ccomps * dcomps);

            auto g_z_x_y_xxxyyz = cbuffer.data(pi_geom_11_off + 539 * ccomps * dcomps);

            auto g_z_x_y_xxxyzz = cbuffer.data(pi_geom_11_off + 540 * ccomps * dcomps);

            auto g_z_x_y_xxxzzz = cbuffer.data(pi_geom_11_off + 541 * ccomps * dcomps);

            auto g_z_x_y_xxyyyy = cbuffer.data(pi_geom_11_off + 542 * ccomps * dcomps);

            auto g_z_x_y_xxyyyz = cbuffer.data(pi_geom_11_off + 543 * ccomps * dcomps);

            auto g_z_x_y_xxyyzz = cbuffer.data(pi_geom_11_off + 544 * ccomps * dcomps);

            auto g_z_x_y_xxyzzz = cbuffer.data(pi_geom_11_off + 545 * ccomps * dcomps);

            auto g_z_x_y_xxzzzz = cbuffer.data(pi_geom_11_off + 546 * ccomps * dcomps);

            auto g_z_x_y_xyyyyy = cbuffer.data(pi_geom_11_off + 547 * ccomps * dcomps);

            auto g_z_x_y_xyyyyz = cbuffer.data(pi_geom_11_off + 548 * ccomps * dcomps);

            auto g_z_x_y_xyyyzz = cbuffer.data(pi_geom_11_off + 549 * ccomps * dcomps);

            auto g_z_x_y_xyyzzz = cbuffer.data(pi_geom_11_off + 550 * ccomps * dcomps);

            auto g_z_x_y_xyzzzz = cbuffer.data(pi_geom_11_off + 551 * ccomps * dcomps);

            auto g_z_x_y_xzzzzz = cbuffer.data(pi_geom_11_off + 552 * ccomps * dcomps);

            auto g_z_x_y_yyyyyy = cbuffer.data(pi_geom_11_off + 553 * ccomps * dcomps);

            auto g_z_x_y_yyyyyz = cbuffer.data(pi_geom_11_off + 554 * ccomps * dcomps);

            auto g_z_x_y_yyyyzz = cbuffer.data(pi_geom_11_off + 555 * ccomps * dcomps);

            auto g_z_x_y_yyyzzz = cbuffer.data(pi_geom_11_off + 556 * ccomps * dcomps);

            auto g_z_x_y_yyzzzz = cbuffer.data(pi_geom_11_off + 557 * ccomps * dcomps);

            auto g_z_x_y_yzzzzz = cbuffer.data(pi_geom_11_off + 558 * ccomps * dcomps);

            auto g_z_x_y_zzzzzz = cbuffer.data(pi_geom_11_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_0_xxxxxx, g_z_x_0_xxxxxxy, g_z_x_0_xxxxxy, g_z_x_0_xxxxxyy, g_z_x_0_xxxxxyz, g_z_x_0_xxxxxz, g_z_x_0_xxxxyy, g_z_x_0_xxxxyyy, g_z_x_0_xxxxyyz, g_z_x_0_xxxxyz, g_z_x_0_xxxxyzz, g_z_x_0_xxxxzz, g_z_x_0_xxxyyy, g_z_x_0_xxxyyyy, g_z_x_0_xxxyyyz, g_z_x_0_xxxyyz, g_z_x_0_xxxyyzz, g_z_x_0_xxxyzz, g_z_x_0_xxxyzzz, g_z_x_0_xxxzzz, g_z_x_0_xxyyyy, g_z_x_0_xxyyyyy, g_z_x_0_xxyyyyz, g_z_x_0_xxyyyz, g_z_x_0_xxyyyzz, g_z_x_0_xxyyzz, g_z_x_0_xxyyzzz, g_z_x_0_xxyzzz, g_z_x_0_xxyzzzz, g_z_x_0_xxzzzz, g_z_x_0_xyyyyy, g_z_x_0_xyyyyyy, g_z_x_0_xyyyyyz, g_z_x_0_xyyyyz, g_z_x_0_xyyyyzz, g_z_x_0_xyyyzz, g_z_x_0_xyyyzzz, g_z_x_0_xyyzzz, g_z_x_0_xyyzzzz, g_z_x_0_xyzzzz, g_z_x_0_xyzzzzz, g_z_x_0_xzzzzz, g_z_x_0_yyyyyy, g_z_x_0_yyyyyyy, g_z_x_0_yyyyyyz, g_z_x_0_yyyyyz, g_z_x_0_yyyyyzz, g_z_x_0_yyyyzz, g_z_x_0_yyyyzzz, g_z_x_0_yyyzzz, g_z_x_0_yyyzzzz, g_z_x_0_yyzzzz, g_z_x_0_yyzzzzz, g_z_x_0_yzzzzz, g_z_x_0_yzzzzzz, g_z_x_0_zzzzzz, g_z_x_y_xxxxxx, g_z_x_y_xxxxxy, g_z_x_y_xxxxxz, g_z_x_y_xxxxyy, g_z_x_y_xxxxyz, g_z_x_y_xxxxzz, g_z_x_y_xxxyyy, g_z_x_y_xxxyyz, g_z_x_y_xxxyzz, g_z_x_y_xxxzzz, g_z_x_y_xxyyyy, g_z_x_y_xxyyyz, g_z_x_y_xxyyzz, g_z_x_y_xxyzzz, g_z_x_y_xxzzzz, g_z_x_y_xyyyyy, g_z_x_y_xyyyyz, g_z_x_y_xyyyzz, g_z_x_y_xyyzzz, g_z_x_y_xyzzzz, g_z_x_y_xzzzzz, g_z_x_y_yyyyyy, g_z_x_y_yyyyyz, g_z_x_y_yyyyzz, g_z_x_y_yyyzzz, g_z_x_y_yyzzzz, g_z_x_y_yzzzzz, g_z_x_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_y_xxxxxx[k] = -g_z_x_0_xxxxxx[k] * ab_y + g_z_x_0_xxxxxxy[k];

                g_z_x_y_xxxxxy[k] = -g_z_x_0_xxxxxy[k] * ab_y + g_z_x_0_xxxxxyy[k];

                g_z_x_y_xxxxxz[k] = -g_z_x_0_xxxxxz[k] * ab_y + g_z_x_0_xxxxxyz[k];

                g_z_x_y_xxxxyy[k] = -g_z_x_0_xxxxyy[k] * ab_y + g_z_x_0_xxxxyyy[k];

                g_z_x_y_xxxxyz[k] = -g_z_x_0_xxxxyz[k] * ab_y + g_z_x_0_xxxxyyz[k];

                g_z_x_y_xxxxzz[k] = -g_z_x_0_xxxxzz[k] * ab_y + g_z_x_0_xxxxyzz[k];

                g_z_x_y_xxxyyy[k] = -g_z_x_0_xxxyyy[k] * ab_y + g_z_x_0_xxxyyyy[k];

                g_z_x_y_xxxyyz[k] = -g_z_x_0_xxxyyz[k] * ab_y + g_z_x_0_xxxyyyz[k];

                g_z_x_y_xxxyzz[k] = -g_z_x_0_xxxyzz[k] * ab_y + g_z_x_0_xxxyyzz[k];

                g_z_x_y_xxxzzz[k] = -g_z_x_0_xxxzzz[k] * ab_y + g_z_x_0_xxxyzzz[k];

                g_z_x_y_xxyyyy[k] = -g_z_x_0_xxyyyy[k] * ab_y + g_z_x_0_xxyyyyy[k];

                g_z_x_y_xxyyyz[k] = -g_z_x_0_xxyyyz[k] * ab_y + g_z_x_0_xxyyyyz[k];

                g_z_x_y_xxyyzz[k] = -g_z_x_0_xxyyzz[k] * ab_y + g_z_x_0_xxyyyzz[k];

                g_z_x_y_xxyzzz[k] = -g_z_x_0_xxyzzz[k] * ab_y + g_z_x_0_xxyyzzz[k];

                g_z_x_y_xxzzzz[k] = -g_z_x_0_xxzzzz[k] * ab_y + g_z_x_0_xxyzzzz[k];

                g_z_x_y_xyyyyy[k] = -g_z_x_0_xyyyyy[k] * ab_y + g_z_x_0_xyyyyyy[k];

                g_z_x_y_xyyyyz[k] = -g_z_x_0_xyyyyz[k] * ab_y + g_z_x_0_xyyyyyz[k];

                g_z_x_y_xyyyzz[k] = -g_z_x_0_xyyyzz[k] * ab_y + g_z_x_0_xyyyyzz[k];

                g_z_x_y_xyyzzz[k] = -g_z_x_0_xyyzzz[k] * ab_y + g_z_x_0_xyyyzzz[k];

                g_z_x_y_xyzzzz[k] = -g_z_x_0_xyzzzz[k] * ab_y + g_z_x_0_xyyzzzz[k];

                g_z_x_y_xzzzzz[k] = -g_z_x_0_xzzzzz[k] * ab_y + g_z_x_0_xyzzzzz[k];

                g_z_x_y_yyyyyy[k] = -g_z_x_0_yyyyyy[k] * ab_y + g_z_x_0_yyyyyyy[k];

                g_z_x_y_yyyyyz[k] = -g_z_x_0_yyyyyz[k] * ab_y + g_z_x_0_yyyyyyz[k];

                g_z_x_y_yyyyzz[k] = -g_z_x_0_yyyyzz[k] * ab_y + g_z_x_0_yyyyyzz[k];

                g_z_x_y_yyyzzz[k] = -g_z_x_0_yyyzzz[k] * ab_y + g_z_x_0_yyyyzzz[k];

                g_z_x_y_yyzzzz[k] = -g_z_x_0_yyzzzz[k] * ab_y + g_z_x_0_yyyzzzz[k];

                g_z_x_y_yzzzzz[k] = -g_z_x_0_yzzzzz[k] * ab_y + g_z_x_0_yyzzzzz[k];

                g_z_x_y_zzzzzz[k] = -g_z_x_0_zzzzzz[k] * ab_y + g_z_x_0_yzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_z_x_z_xxxxxx = cbuffer.data(pi_geom_11_off + 560 * ccomps * dcomps);

            auto g_z_x_z_xxxxxy = cbuffer.data(pi_geom_11_off + 561 * ccomps * dcomps);

            auto g_z_x_z_xxxxxz = cbuffer.data(pi_geom_11_off + 562 * ccomps * dcomps);

            auto g_z_x_z_xxxxyy = cbuffer.data(pi_geom_11_off + 563 * ccomps * dcomps);

            auto g_z_x_z_xxxxyz = cbuffer.data(pi_geom_11_off + 564 * ccomps * dcomps);

            auto g_z_x_z_xxxxzz = cbuffer.data(pi_geom_11_off + 565 * ccomps * dcomps);

            auto g_z_x_z_xxxyyy = cbuffer.data(pi_geom_11_off + 566 * ccomps * dcomps);

            auto g_z_x_z_xxxyyz = cbuffer.data(pi_geom_11_off + 567 * ccomps * dcomps);

            auto g_z_x_z_xxxyzz = cbuffer.data(pi_geom_11_off + 568 * ccomps * dcomps);

            auto g_z_x_z_xxxzzz = cbuffer.data(pi_geom_11_off + 569 * ccomps * dcomps);

            auto g_z_x_z_xxyyyy = cbuffer.data(pi_geom_11_off + 570 * ccomps * dcomps);

            auto g_z_x_z_xxyyyz = cbuffer.data(pi_geom_11_off + 571 * ccomps * dcomps);

            auto g_z_x_z_xxyyzz = cbuffer.data(pi_geom_11_off + 572 * ccomps * dcomps);

            auto g_z_x_z_xxyzzz = cbuffer.data(pi_geom_11_off + 573 * ccomps * dcomps);

            auto g_z_x_z_xxzzzz = cbuffer.data(pi_geom_11_off + 574 * ccomps * dcomps);

            auto g_z_x_z_xyyyyy = cbuffer.data(pi_geom_11_off + 575 * ccomps * dcomps);

            auto g_z_x_z_xyyyyz = cbuffer.data(pi_geom_11_off + 576 * ccomps * dcomps);

            auto g_z_x_z_xyyyzz = cbuffer.data(pi_geom_11_off + 577 * ccomps * dcomps);

            auto g_z_x_z_xyyzzz = cbuffer.data(pi_geom_11_off + 578 * ccomps * dcomps);

            auto g_z_x_z_xyzzzz = cbuffer.data(pi_geom_11_off + 579 * ccomps * dcomps);

            auto g_z_x_z_xzzzzz = cbuffer.data(pi_geom_11_off + 580 * ccomps * dcomps);

            auto g_z_x_z_yyyyyy = cbuffer.data(pi_geom_11_off + 581 * ccomps * dcomps);

            auto g_z_x_z_yyyyyz = cbuffer.data(pi_geom_11_off + 582 * ccomps * dcomps);

            auto g_z_x_z_yyyyzz = cbuffer.data(pi_geom_11_off + 583 * ccomps * dcomps);

            auto g_z_x_z_yyyzzz = cbuffer.data(pi_geom_11_off + 584 * ccomps * dcomps);

            auto g_z_x_z_yyzzzz = cbuffer.data(pi_geom_11_off + 585 * ccomps * dcomps);

            auto g_z_x_z_yzzzzz = cbuffer.data(pi_geom_11_off + 586 * ccomps * dcomps);

            auto g_z_x_z_zzzzzz = cbuffer.data(pi_geom_11_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_0_xxxxxx, g_0_x_0_xxxxxy, g_0_x_0_xxxxxz, g_0_x_0_xxxxyy, g_0_x_0_xxxxyz, g_0_x_0_xxxxzz, g_0_x_0_xxxyyy, g_0_x_0_xxxyyz, g_0_x_0_xxxyzz, g_0_x_0_xxxzzz, g_0_x_0_xxyyyy, g_0_x_0_xxyyyz, g_0_x_0_xxyyzz, g_0_x_0_xxyzzz, g_0_x_0_xxzzzz, g_0_x_0_xyyyyy, g_0_x_0_xyyyyz, g_0_x_0_xyyyzz, g_0_x_0_xyyzzz, g_0_x_0_xyzzzz, g_0_x_0_xzzzzz, g_0_x_0_yyyyyy, g_0_x_0_yyyyyz, g_0_x_0_yyyyzz, g_0_x_0_yyyzzz, g_0_x_0_yyzzzz, g_0_x_0_yzzzzz, g_0_x_0_zzzzzz, g_z_x_0_xxxxxx, g_z_x_0_xxxxxxz, g_z_x_0_xxxxxy, g_z_x_0_xxxxxyz, g_z_x_0_xxxxxz, g_z_x_0_xxxxxzz, g_z_x_0_xxxxyy, g_z_x_0_xxxxyyz, g_z_x_0_xxxxyz, g_z_x_0_xxxxyzz, g_z_x_0_xxxxzz, g_z_x_0_xxxxzzz, g_z_x_0_xxxyyy, g_z_x_0_xxxyyyz, g_z_x_0_xxxyyz, g_z_x_0_xxxyyzz, g_z_x_0_xxxyzz, g_z_x_0_xxxyzzz, g_z_x_0_xxxzzz, g_z_x_0_xxxzzzz, g_z_x_0_xxyyyy, g_z_x_0_xxyyyyz, g_z_x_0_xxyyyz, g_z_x_0_xxyyyzz, g_z_x_0_xxyyzz, g_z_x_0_xxyyzzz, g_z_x_0_xxyzzz, g_z_x_0_xxyzzzz, g_z_x_0_xxzzzz, g_z_x_0_xxzzzzz, g_z_x_0_xyyyyy, g_z_x_0_xyyyyyz, g_z_x_0_xyyyyz, g_z_x_0_xyyyyzz, g_z_x_0_xyyyzz, g_z_x_0_xyyyzzz, g_z_x_0_xyyzzz, g_z_x_0_xyyzzzz, g_z_x_0_xyzzzz, g_z_x_0_xyzzzzz, g_z_x_0_xzzzzz, g_z_x_0_xzzzzzz, g_z_x_0_yyyyyy, g_z_x_0_yyyyyyz, g_z_x_0_yyyyyz, g_z_x_0_yyyyyzz, g_z_x_0_yyyyzz, g_z_x_0_yyyyzzz, g_z_x_0_yyyzzz, g_z_x_0_yyyzzzz, g_z_x_0_yyzzzz, g_z_x_0_yyzzzzz, g_z_x_0_yzzzzz, g_z_x_0_yzzzzzz, g_z_x_0_zzzzzz, g_z_x_0_zzzzzzz, g_z_x_z_xxxxxx, g_z_x_z_xxxxxy, g_z_x_z_xxxxxz, g_z_x_z_xxxxyy, g_z_x_z_xxxxyz, g_z_x_z_xxxxzz, g_z_x_z_xxxyyy, g_z_x_z_xxxyyz, g_z_x_z_xxxyzz, g_z_x_z_xxxzzz, g_z_x_z_xxyyyy, g_z_x_z_xxyyyz, g_z_x_z_xxyyzz, g_z_x_z_xxyzzz, g_z_x_z_xxzzzz, g_z_x_z_xyyyyy, g_z_x_z_xyyyyz, g_z_x_z_xyyyzz, g_z_x_z_xyyzzz, g_z_x_z_xyzzzz, g_z_x_z_xzzzzz, g_z_x_z_yyyyyy, g_z_x_z_yyyyyz, g_z_x_z_yyyyzz, g_z_x_z_yyyzzz, g_z_x_z_yyzzzz, g_z_x_z_yzzzzz, g_z_x_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_z_xxxxxx[k] = -g_0_x_0_xxxxxx[k] - g_z_x_0_xxxxxx[k] * ab_z + g_z_x_0_xxxxxxz[k];

                g_z_x_z_xxxxxy[k] = -g_0_x_0_xxxxxy[k] - g_z_x_0_xxxxxy[k] * ab_z + g_z_x_0_xxxxxyz[k];

                g_z_x_z_xxxxxz[k] = -g_0_x_0_xxxxxz[k] - g_z_x_0_xxxxxz[k] * ab_z + g_z_x_0_xxxxxzz[k];

                g_z_x_z_xxxxyy[k] = -g_0_x_0_xxxxyy[k] - g_z_x_0_xxxxyy[k] * ab_z + g_z_x_0_xxxxyyz[k];

                g_z_x_z_xxxxyz[k] = -g_0_x_0_xxxxyz[k] - g_z_x_0_xxxxyz[k] * ab_z + g_z_x_0_xxxxyzz[k];

                g_z_x_z_xxxxzz[k] = -g_0_x_0_xxxxzz[k] - g_z_x_0_xxxxzz[k] * ab_z + g_z_x_0_xxxxzzz[k];

                g_z_x_z_xxxyyy[k] = -g_0_x_0_xxxyyy[k] - g_z_x_0_xxxyyy[k] * ab_z + g_z_x_0_xxxyyyz[k];

                g_z_x_z_xxxyyz[k] = -g_0_x_0_xxxyyz[k] - g_z_x_0_xxxyyz[k] * ab_z + g_z_x_0_xxxyyzz[k];

                g_z_x_z_xxxyzz[k] = -g_0_x_0_xxxyzz[k] - g_z_x_0_xxxyzz[k] * ab_z + g_z_x_0_xxxyzzz[k];

                g_z_x_z_xxxzzz[k] = -g_0_x_0_xxxzzz[k] - g_z_x_0_xxxzzz[k] * ab_z + g_z_x_0_xxxzzzz[k];

                g_z_x_z_xxyyyy[k] = -g_0_x_0_xxyyyy[k] - g_z_x_0_xxyyyy[k] * ab_z + g_z_x_0_xxyyyyz[k];

                g_z_x_z_xxyyyz[k] = -g_0_x_0_xxyyyz[k] - g_z_x_0_xxyyyz[k] * ab_z + g_z_x_0_xxyyyzz[k];

                g_z_x_z_xxyyzz[k] = -g_0_x_0_xxyyzz[k] - g_z_x_0_xxyyzz[k] * ab_z + g_z_x_0_xxyyzzz[k];

                g_z_x_z_xxyzzz[k] = -g_0_x_0_xxyzzz[k] - g_z_x_0_xxyzzz[k] * ab_z + g_z_x_0_xxyzzzz[k];

                g_z_x_z_xxzzzz[k] = -g_0_x_0_xxzzzz[k] - g_z_x_0_xxzzzz[k] * ab_z + g_z_x_0_xxzzzzz[k];

                g_z_x_z_xyyyyy[k] = -g_0_x_0_xyyyyy[k] - g_z_x_0_xyyyyy[k] * ab_z + g_z_x_0_xyyyyyz[k];

                g_z_x_z_xyyyyz[k] = -g_0_x_0_xyyyyz[k] - g_z_x_0_xyyyyz[k] * ab_z + g_z_x_0_xyyyyzz[k];

                g_z_x_z_xyyyzz[k] = -g_0_x_0_xyyyzz[k] - g_z_x_0_xyyyzz[k] * ab_z + g_z_x_0_xyyyzzz[k];

                g_z_x_z_xyyzzz[k] = -g_0_x_0_xyyzzz[k] - g_z_x_0_xyyzzz[k] * ab_z + g_z_x_0_xyyzzzz[k];

                g_z_x_z_xyzzzz[k] = -g_0_x_0_xyzzzz[k] - g_z_x_0_xyzzzz[k] * ab_z + g_z_x_0_xyzzzzz[k];

                g_z_x_z_xzzzzz[k] = -g_0_x_0_xzzzzz[k] - g_z_x_0_xzzzzz[k] * ab_z + g_z_x_0_xzzzzzz[k];

                g_z_x_z_yyyyyy[k] = -g_0_x_0_yyyyyy[k] - g_z_x_0_yyyyyy[k] * ab_z + g_z_x_0_yyyyyyz[k];

                g_z_x_z_yyyyyz[k] = -g_0_x_0_yyyyyz[k] - g_z_x_0_yyyyyz[k] * ab_z + g_z_x_0_yyyyyzz[k];

                g_z_x_z_yyyyzz[k] = -g_0_x_0_yyyyzz[k] - g_z_x_0_yyyyzz[k] * ab_z + g_z_x_0_yyyyzzz[k];

                g_z_x_z_yyyzzz[k] = -g_0_x_0_yyyzzz[k] - g_z_x_0_yyyzzz[k] * ab_z + g_z_x_0_yyyzzzz[k];

                g_z_x_z_yyzzzz[k] = -g_0_x_0_yyzzzz[k] - g_z_x_0_yyzzzz[k] * ab_z + g_z_x_0_yyzzzzz[k];

                g_z_x_z_yzzzzz[k] = -g_0_x_0_yzzzzz[k] - g_z_x_0_yzzzzz[k] * ab_z + g_z_x_0_yzzzzzz[k];

                g_z_x_z_zzzzzz[k] = -g_0_x_0_zzzzzz[k] - g_z_x_0_zzzzzz[k] * ab_z + g_z_x_0_zzzzzzz[k];
            }

            /// Set up 588-616 components of targeted buffer : cbuffer.data(

            auto g_z_y_x_xxxxxx = cbuffer.data(pi_geom_11_off + 588 * ccomps * dcomps);

            auto g_z_y_x_xxxxxy = cbuffer.data(pi_geom_11_off + 589 * ccomps * dcomps);

            auto g_z_y_x_xxxxxz = cbuffer.data(pi_geom_11_off + 590 * ccomps * dcomps);

            auto g_z_y_x_xxxxyy = cbuffer.data(pi_geom_11_off + 591 * ccomps * dcomps);

            auto g_z_y_x_xxxxyz = cbuffer.data(pi_geom_11_off + 592 * ccomps * dcomps);

            auto g_z_y_x_xxxxzz = cbuffer.data(pi_geom_11_off + 593 * ccomps * dcomps);

            auto g_z_y_x_xxxyyy = cbuffer.data(pi_geom_11_off + 594 * ccomps * dcomps);

            auto g_z_y_x_xxxyyz = cbuffer.data(pi_geom_11_off + 595 * ccomps * dcomps);

            auto g_z_y_x_xxxyzz = cbuffer.data(pi_geom_11_off + 596 * ccomps * dcomps);

            auto g_z_y_x_xxxzzz = cbuffer.data(pi_geom_11_off + 597 * ccomps * dcomps);

            auto g_z_y_x_xxyyyy = cbuffer.data(pi_geom_11_off + 598 * ccomps * dcomps);

            auto g_z_y_x_xxyyyz = cbuffer.data(pi_geom_11_off + 599 * ccomps * dcomps);

            auto g_z_y_x_xxyyzz = cbuffer.data(pi_geom_11_off + 600 * ccomps * dcomps);

            auto g_z_y_x_xxyzzz = cbuffer.data(pi_geom_11_off + 601 * ccomps * dcomps);

            auto g_z_y_x_xxzzzz = cbuffer.data(pi_geom_11_off + 602 * ccomps * dcomps);

            auto g_z_y_x_xyyyyy = cbuffer.data(pi_geom_11_off + 603 * ccomps * dcomps);

            auto g_z_y_x_xyyyyz = cbuffer.data(pi_geom_11_off + 604 * ccomps * dcomps);

            auto g_z_y_x_xyyyzz = cbuffer.data(pi_geom_11_off + 605 * ccomps * dcomps);

            auto g_z_y_x_xyyzzz = cbuffer.data(pi_geom_11_off + 606 * ccomps * dcomps);

            auto g_z_y_x_xyzzzz = cbuffer.data(pi_geom_11_off + 607 * ccomps * dcomps);

            auto g_z_y_x_xzzzzz = cbuffer.data(pi_geom_11_off + 608 * ccomps * dcomps);

            auto g_z_y_x_yyyyyy = cbuffer.data(pi_geom_11_off + 609 * ccomps * dcomps);

            auto g_z_y_x_yyyyyz = cbuffer.data(pi_geom_11_off + 610 * ccomps * dcomps);

            auto g_z_y_x_yyyyzz = cbuffer.data(pi_geom_11_off + 611 * ccomps * dcomps);

            auto g_z_y_x_yyyzzz = cbuffer.data(pi_geom_11_off + 612 * ccomps * dcomps);

            auto g_z_y_x_yyzzzz = cbuffer.data(pi_geom_11_off + 613 * ccomps * dcomps);

            auto g_z_y_x_yzzzzz = cbuffer.data(pi_geom_11_off + 614 * ccomps * dcomps);

            auto g_z_y_x_zzzzzz = cbuffer.data(pi_geom_11_off + 615 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_0_xxxxxx, g_z_y_0_xxxxxxx, g_z_y_0_xxxxxxy, g_z_y_0_xxxxxxz, g_z_y_0_xxxxxy, g_z_y_0_xxxxxyy, g_z_y_0_xxxxxyz, g_z_y_0_xxxxxz, g_z_y_0_xxxxxzz, g_z_y_0_xxxxyy, g_z_y_0_xxxxyyy, g_z_y_0_xxxxyyz, g_z_y_0_xxxxyz, g_z_y_0_xxxxyzz, g_z_y_0_xxxxzz, g_z_y_0_xxxxzzz, g_z_y_0_xxxyyy, g_z_y_0_xxxyyyy, g_z_y_0_xxxyyyz, g_z_y_0_xxxyyz, g_z_y_0_xxxyyzz, g_z_y_0_xxxyzz, g_z_y_0_xxxyzzz, g_z_y_0_xxxzzz, g_z_y_0_xxxzzzz, g_z_y_0_xxyyyy, g_z_y_0_xxyyyyy, g_z_y_0_xxyyyyz, g_z_y_0_xxyyyz, g_z_y_0_xxyyyzz, g_z_y_0_xxyyzz, g_z_y_0_xxyyzzz, g_z_y_0_xxyzzz, g_z_y_0_xxyzzzz, g_z_y_0_xxzzzz, g_z_y_0_xxzzzzz, g_z_y_0_xyyyyy, g_z_y_0_xyyyyyy, g_z_y_0_xyyyyyz, g_z_y_0_xyyyyz, g_z_y_0_xyyyyzz, g_z_y_0_xyyyzz, g_z_y_0_xyyyzzz, g_z_y_0_xyyzzz, g_z_y_0_xyyzzzz, g_z_y_0_xyzzzz, g_z_y_0_xyzzzzz, g_z_y_0_xzzzzz, g_z_y_0_xzzzzzz, g_z_y_0_yyyyyy, g_z_y_0_yyyyyz, g_z_y_0_yyyyzz, g_z_y_0_yyyzzz, g_z_y_0_yyzzzz, g_z_y_0_yzzzzz, g_z_y_0_zzzzzz, g_z_y_x_xxxxxx, g_z_y_x_xxxxxy, g_z_y_x_xxxxxz, g_z_y_x_xxxxyy, g_z_y_x_xxxxyz, g_z_y_x_xxxxzz, g_z_y_x_xxxyyy, g_z_y_x_xxxyyz, g_z_y_x_xxxyzz, g_z_y_x_xxxzzz, g_z_y_x_xxyyyy, g_z_y_x_xxyyyz, g_z_y_x_xxyyzz, g_z_y_x_xxyzzz, g_z_y_x_xxzzzz, g_z_y_x_xyyyyy, g_z_y_x_xyyyyz, g_z_y_x_xyyyzz, g_z_y_x_xyyzzz, g_z_y_x_xyzzzz, g_z_y_x_xzzzzz, g_z_y_x_yyyyyy, g_z_y_x_yyyyyz, g_z_y_x_yyyyzz, g_z_y_x_yyyzzz, g_z_y_x_yyzzzz, g_z_y_x_yzzzzz, g_z_y_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_x_xxxxxx[k] = -g_z_y_0_xxxxxx[k] * ab_x + g_z_y_0_xxxxxxx[k];

                g_z_y_x_xxxxxy[k] = -g_z_y_0_xxxxxy[k] * ab_x + g_z_y_0_xxxxxxy[k];

                g_z_y_x_xxxxxz[k] = -g_z_y_0_xxxxxz[k] * ab_x + g_z_y_0_xxxxxxz[k];

                g_z_y_x_xxxxyy[k] = -g_z_y_0_xxxxyy[k] * ab_x + g_z_y_0_xxxxxyy[k];

                g_z_y_x_xxxxyz[k] = -g_z_y_0_xxxxyz[k] * ab_x + g_z_y_0_xxxxxyz[k];

                g_z_y_x_xxxxzz[k] = -g_z_y_0_xxxxzz[k] * ab_x + g_z_y_0_xxxxxzz[k];

                g_z_y_x_xxxyyy[k] = -g_z_y_0_xxxyyy[k] * ab_x + g_z_y_0_xxxxyyy[k];

                g_z_y_x_xxxyyz[k] = -g_z_y_0_xxxyyz[k] * ab_x + g_z_y_0_xxxxyyz[k];

                g_z_y_x_xxxyzz[k] = -g_z_y_0_xxxyzz[k] * ab_x + g_z_y_0_xxxxyzz[k];

                g_z_y_x_xxxzzz[k] = -g_z_y_0_xxxzzz[k] * ab_x + g_z_y_0_xxxxzzz[k];

                g_z_y_x_xxyyyy[k] = -g_z_y_0_xxyyyy[k] * ab_x + g_z_y_0_xxxyyyy[k];

                g_z_y_x_xxyyyz[k] = -g_z_y_0_xxyyyz[k] * ab_x + g_z_y_0_xxxyyyz[k];

                g_z_y_x_xxyyzz[k] = -g_z_y_0_xxyyzz[k] * ab_x + g_z_y_0_xxxyyzz[k];

                g_z_y_x_xxyzzz[k] = -g_z_y_0_xxyzzz[k] * ab_x + g_z_y_0_xxxyzzz[k];

                g_z_y_x_xxzzzz[k] = -g_z_y_0_xxzzzz[k] * ab_x + g_z_y_0_xxxzzzz[k];

                g_z_y_x_xyyyyy[k] = -g_z_y_0_xyyyyy[k] * ab_x + g_z_y_0_xxyyyyy[k];

                g_z_y_x_xyyyyz[k] = -g_z_y_0_xyyyyz[k] * ab_x + g_z_y_0_xxyyyyz[k];

                g_z_y_x_xyyyzz[k] = -g_z_y_0_xyyyzz[k] * ab_x + g_z_y_0_xxyyyzz[k];

                g_z_y_x_xyyzzz[k] = -g_z_y_0_xyyzzz[k] * ab_x + g_z_y_0_xxyyzzz[k];

                g_z_y_x_xyzzzz[k] = -g_z_y_0_xyzzzz[k] * ab_x + g_z_y_0_xxyzzzz[k];

                g_z_y_x_xzzzzz[k] = -g_z_y_0_xzzzzz[k] * ab_x + g_z_y_0_xxzzzzz[k];

                g_z_y_x_yyyyyy[k] = -g_z_y_0_yyyyyy[k] * ab_x + g_z_y_0_xyyyyyy[k];

                g_z_y_x_yyyyyz[k] = -g_z_y_0_yyyyyz[k] * ab_x + g_z_y_0_xyyyyyz[k];

                g_z_y_x_yyyyzz[k] = -g_z_y_0_yyyyzz[k] * ab_x + g_z_y_0_xyyyyzz[k];

                g_z_y_x_yyyzzz[k] = -g_z_y_0_yyyzzz[k] * ab_x + g_z_y_0_xyyyzzz[k];

                g_z_y_x_yyzzzz[k] = -g_z_y_0_yyzzzz[k] * ab_x + g_z_y_0_xyyzzzz[k];

                g_z_y_x_yzzzzz[k] = -g_z_y_0_yzzzzz[k] * ab_x + g_z_y_0_xyzzzzz[k];

                g_z_y_x_zzzzzz[k] = -g_z_y_0_zzzzzz[k] * ab_x + g_z_y_0_xzzzzzz[k];
            }

            /// Set up 616-644 components of targeted buffer : cbuffer.data(

            auto g_z_y_y_xxxxxx = cbuffer.data(pi_geom_11_off + 616 * ccomps * dcomps);

            auto g_z_y_y_xxxxxy = cbuffer.data(pi_geom_11_off + 617 * ccomps * dcomps);

            auto g_z_y_y_xxxxxz = cbuffer.data(pi_geom_11_off + 618 * ccomps * dcomps);

            auto g_z_y_y_xxxxyy = cbuffer.data(pi_geom_11_off + 619 * ccomps * dcomps);

            auto g_z_y_y_xxxxyz = cbuffer.data(pi_geom_11_off + 620 * ccomps * dcomps);

            auto g_z_y_y_xxxxzz = cbuffer.data(pi_geom_11_off + 621 * ccomps * dcomps);

            auto g_z_y_y_xxxyyy = cbuffer.data(pi_geom_11_off + 622 * ccomps * dcomps);

            auto g_z_y_y_xxxyyz = cbuffer.data(pi_geom_11_off + 623 * ccomps * dcomps);

            auto g_z_y_y_xxxyzz = cbuffer.data(pi_geom_11_off + 624 * ccomps * dcomps);

            auto g_z_y_y_xxxzzz = cbuffer.data(pi_geom_11_off + 625 * ccomps * dcomps);

            auto g_z_y_y_xxyyyy = cbuffer.data(pi_geom_11_off + 626 * ccomps * dcomps);

            auto g_z_y_y_xxyyyz = cbuffer.data(pi_geom_11_off + 627 * ccomps * dcomps);

            auto g_z_y_y_xxyyzz = cbuffer.data(pi_geom_11_off + 628 * ccomps * dcomps);

            auto g_z_y_y_xxyzzz = cbuffer.data(pi_geom_11_off + 629 * ccomps * dcomps);

            auto g_z_y_y_xxzzzz = cbuffer.data(pi_geom_11_off + 630 * ccomps * dcomps);

            auto g_z_y_y_xyyyyy = cbuffer.data(pi_geom_11_off + 631 * ccomps * dcomps);

            auto g_z_y_y_xyyyyz = cbuffer.data(pi_geom_11_off + 632 * ccomps * dcomps);

            auto g_z_y_y_xyyyzz = cbuffer.data(pi_geom_11_off + 633 * ccomps * dcomps);

            auto g_z_y_y_xyyzzz = cbuffer.data(pi_geom_11_off + 634 * ccomps * dcomps);

            auto g_z_y_y_xyzzzz = cbuffer.data(pi_geom_11_off + 635 * ccomps * dcomps);

            auto g_z_y_y_xzzzzz = cbuffer.data(pi_geom_11_off + 636 * ccomps * dcomps);

            auto g_z_y_y_yyyyyy = cbuffer.data(pi_geom_11_off + 637 * ccomps * dcomps);

            auto g_z_y_y_yyyyyz = cbuffer.data(pi_geom_11_off + 638 * ccomps * dcomps);

            auto g_z_y_y_yyyyzz = cbuffer.data(pi_geom_11_off + 639 * ccomps * dcomps);

            auto g_z_y_y_yyyzzz = cbuffer.data(pi_geom_11_off + 640 * ccomps * dcomps);

            auto g_z_y_y_yyzzzz = cbuffer.data(pi_geom_11_off + 641 * ccomps * dcomps);

            auto g_z_y_y_yzzzzz = cbuffer.data(pi_geom_11_off + 642 * ccomps * dcomps);

            auto g_z_y_y_zzzzzz = cbuffer.data(pi_geom_11_off + 643 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_0_xxxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxxz, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxzz, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzzz, g_z_y_0_xxxxxx, g_z_y_0_xxxxxxy, g_z_y_0_xxxxxy, g_z_y_0_xxxxxyy, g_z_y_0_xxxxxyz, g_z_y_0_xxxxxz, g_z_y_0_xxxxyy, g_z_y_0_xxxxyyy, g_z_y_0_xxxxyyz, g_z_y_0_xxxxyz, g_z_y_0_xxxxyzz, g_z_y_0_xxxxzz, g_z_y_0_xxxyyy, g_z_y_0_xxxyyyy, g_z_y_0_xxxyyyz, g_z_y_0_xxxyyz, g_z_y_0_xxxyyzz, g_z_y_0_xxxyzz, g_z_y_0_xxxyzzz, g_z_y_0_xxxzzz, g_z_y_0_xxyyyy, g_z_y_0_xxyyyyy, g_z_y_0_xxyyyyz, g_z_y_0_xxyyyz, g_z_y_0_xxyyyzz, g_z_y_0_xxyyzz, g_z_y_0_xxyyzzz, g_z_y_0_xxyzzz, g_z_y_0_xxyzzzz, g_z_y_0_xxzzzz, g_z_y_0_xyyyyy, g_z_y_0_xyyyyyy, g_z_y_0_xyyyyyz, g_z_y_0_xyyyyz, g_z_y_0_xyyyyzz, g_z_y_0_xyyyzz, g_z_y_0_xyyyzzz, g_z_y_0_xyyzzz, g_z_y_0_xyyzzzz, g_z_y_0_xyzzzz, g_z_y_0_xyzzzzz, g_z_y_0_xzzzzz, g_z_y_0_yyyyyy, g_z_y_0_yyyyyyy, g_z_y_0_yyyyyyz, g_z_y_0_yyyyyz, g_z_y_0_yyyyyzz, g_z_y_0_yyyyzz, g_z_y_0_yyyyzzz, g_z_y_0_yyyzzz, g_z_y_0_yyyzzzz, g_z_y_0_yyzzzz, g_z_y_0_yyzzzzz, g_z_y_0_yzzzzz, g_z_y_0_yzzzzzz, g_z_y_0_zzzzzz, g_z_y_y_xxxxxx, g_z_y_y_xxxxxy, g_z_y_y_xxxxxz, g_z_y_y_xxxxyy, g_z_y_y_xxxxyz, g_z_y_y_xxxxzz, g_z_y_y_xxxyyy, g_z_y_y_xxxyyz, g_z_y_y_xxxyzz, g_z_y_y_xxxzzz, g_z_y_y_xxyyyy, g_z_y_y_xxyyyz, g_z_y_y_xxyyzz, g_z_y_y_xxyzzz, g_z_y_y_xxzzzz, g_z_y_y_xyyyyy, g_z_y_y_xyyyyz, g_z_y_y_xyyyzz, g_z_y_y_xyyzzz, g_z_y_y_xyzzzz, g_z_y_y_xzzzzz, g_z_y_y_yyyyyy, g_z_y_y_yyyyyz, g_z_y_y_yyyyzz, g_z_y_y_yyyzzz, g_z_y_y_yyzzzz, g_z_y_y_yzzzzz, g_z_y_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_y_xxxxxx[k] = g_z_0_0_xxxxxx[k] - g_z_y_0_xxxxxx[k] * ab_y + g_z_y_0_xxxxxxy[k];

                g_z_y_y_xxxxxy[k] = g_z_0_0_xxxxxy[k] - g_z_y_0_xxxxxy[k] * ab_y + g_z_y_0_xxxxxyy[k];

                g_z_y_y_xxxxxz[k] = g_z_0_0_xxxxxz[k] - g_z_y_0_xxxxxz[k] * ab_y + g_z_y_0_xxxxxyz[k];

                g_z_y_y_xxxxyy[k] = g_z_0_0_xxxxyy[k] - g_z_y_0_xxxxyy[k] * ab_y + g_z_y_0_xxxxyyy[k];

                g_z_y_y_xxxxyz[k] = g_z_0_0_xxxxyz[k] - g_z_y_0_xxxxyz[k] * ab_y + g_z_y_0_xxxxyyz[k];

                g_z_y_y_xxxxzz[k] = g_z_0_0_xxxxzz[k] - g_z_y_0_xxxxzz[k] * ab_y + g_z_y_0_xxxxyzz[k];

                g_z_y_y_xxxyyy[k] = g_z_0_0_xxxyyy[k] - g_z_y_0_xxxyyy[k] * ab_y + g_z_y_0_xxxyyyy[k];

                g_z_y_y_xxxyyz[k] = g_z_0_0_xxxyyz[k] - g_z_y_0_xxxyyz[k] * ab_y + g_z_y_0_xxxyyyz[k];

                g_z_y_y_xxxyzz[k] = g_z_0_0_xxxyzz[k] - g_z_y_0_xxxyzz[k] * ab_y + g_z_y_0_xxxyyzz[k];

                g_z_y_y_xxxzzz[k] = g_z_0_0_xxxzzz[k] - g_z_y_0_xxxzzz[k] * ab_y + g_z_y_0_xxxyzzz[k];

                g_z_y_y_xxyyyy[k] = g_z_0_0_xxyyyy[k] - g_z_y_0_xxyyyy[k] * ab_y + g_z_y_0_xxyyyyy[k];

                g_z_y_y_xxyyyz[k] = g_z_0_0_xxyyyz[k] - g_z_y_0_xxyyyz[k] * ab_y + g_z_y_0_xxyyyyz[k];

                g_z_y_y_xxyyzz[k] = g_z_0_0_xxyyzz[k] - g_z_y_0_xxyyzz[k] * ab_y + g_z_y_0_xxyyyzz[k];

                g_z_y_y_xxyzzz[k] = g_z_0_0_xxyzzz[k] - g_z_y_0_xxyzzz[k] * ab_y + g_z_y_0_xxyyzzz[k];

                g_z_y_y_xxzzzz[k] = g_z_0_0_xxzzzz[k] - g_z_y_0_xxzzzz[k] * ab_y + g_z_y_0_xxyzzzz[k];

                g_z_y_y_xyyyyy[k] = g_z_0_0_xyyyyy[k] - g_z_y_0_xyyyyy[k] * ab_y + g_z_y_0_xyyyyyy[k];

                g_z_y_y_xyyyyz[k] = g_z_0_0_xyyyyz[k] - g_z_y_0_xyyyyz[k] * ab_y + g_z_y_0_xyyyyyz[k];

                g_z_y_y_xyyyzz[k] = g_z_0_0_xyyyzz[k] - g_z_y_0_xyyyzz[k] * ab_y + g_z_y_0_xyyyyzz[k];

                g_z_y_y_xyyzzz[k] = g_z_0_0_xyyzzz[k] - g_z_y_0_xyyzzz[k] * ab_y + g_z_y_0_xyyyzzz[k];

                g_z_y_y_xyzzzz[k] = g_z_0_0_xyzzzz[k] - g_z_y_0_xyzzzz[k] * ab_y + g_z_y_0_xyyzzzz[k];

                g_z_y_y_xzzzzz[k] = g_z_0_0_xzzzzz[k] - g_z_y_0_xzzzzz[k] * ab_y + g_z_y_0_xyzzzzz[k];

                g_z_y_y_yyyyyy[k] = g_z_0_0_yyyyyy[k] - g_z_y_0_yyyyyy[k] * ab_y + g_z_y_0_yyyyyyy[k];

                g_z_y_y_yyyyyz[k] = g_z_0_0_yyyyyz[k] - g_z_y_0_yyyyyz[k] * ab_y + g_z_y_0_yyyyyyz[k];

                g_z_y_y_yyyyzz[k] = g_z_0_0_yyyyzz[k] - g_z_y_0_yyyyzz[k] * ab_y + g_z_y_0_yyyyyzz[k];

                g_z_y_y_yyyzzz[k] = g_z_0_0_yyyzzz[k] - g_z_y_0_yyyzzz[k] * ab_y + g_z_y_0_yyyyzzz[k];

                g_z_y_y_yyzzzz[k] = g_z_0_0_yyzzzz[k] - g_z_y_0_yyzzzz[k] * ab_y + g_z_y_0_yyyzzzz[k];

                g_z_y_y_yzzzzz[k] = g_z_0_0_yzzzzz[k] - g_z_y_0_yzzzzz[k] * ab_y + g_z_y_0_yyzzzzz[k];

                g_z_y_y_zzzzzz[k] = g_z_0_0_zzzzzz[k] - g_z_y_0_zzzzzz[k] * ab_y + g_z_y_0_yzzzzzz[k];
            }

            /// Set up 644-672 components of targeted buffer : cbuffer.data(

            auto g_z_y_z_xxxxxx = cbuffer.data(pi_geom_11_off + 644 * ccomps * dcomps);

            auto g_z_y_z_xxxxxy = cbuffer.data(pi_geom_11_off + 645 * ccomps * dcomps);

            auto g_z_y_z_xxxxxz = cbuffer.data(pi_geom_11_off + 646 * ccomps * dcomps);

            auto g_z_y_z_xxxxyy = cbuffer.data(pi_geom_11_off + 647 * ccomps * dcomps);

            auto g_z_y_z_xxxxyz = cbuffer.data(pi_geom_11_off + 648 * ccomps * dcomps);

            auto g_z_y_z_xxxxzz = cbuffer.data(pi_geom_11_off + 649 * ccomps * dcomps);

            auto g_z_y_z_xxxyyy = cbuffer.data(pi_geom_11_off + 650 * ccomps * dcomps);

            auto g_z_y_z_xxxyyz = cbuffer.data(pi_geom_11_off + 651 * ccomps * dcomps);

            auto g_z_y_z_xxxyzz = cbuffer.data(pi_geom_11_off + 652 * ccomps * dcomps);

            auto g_z_y_z_xxxzzz = cbuffer.data(pi_geom_11_off + 653 * ccomps * dcomps);

            auto g_z_y_z_xxyyyy = cbuffer.data(pi_geom_11_off + 654 * ccomps * dcomps);

            auto g_z_y_z_xxyyyz = cbuffer.data(pi_geom_11_off + 655 * ccomps * dcomps);

            auto g_z_y_z_xxyyzz = cbuffer.data(pi_geom_11_off + 656 * ccomps * dcomps);

            auto g_z_y_z_xxyzzz = cbuffer.data(pi_geom_11_off + 657 * ccomps * dcomps);

            auto g_z_y_z_xxzzzz = cbuffer.data(pi_geom_11_off + 658 * ccomps * dcomps);

            auto g_z_y_z_xyyyyy = cbuffer.data(pi_geom_11_off + 659 * ccomps * dcomps);

            auto g_z_y_z_xyyyyz = cbuffer.data(pi_geom_11_off + 660 * ccomps * dcomps);

            auto g_z_y_z_xyyyzz = cbuffer.data(pi_geom_11_off + 661 * ccomps * dcomps);

            auto g_z_y_z_xyyzzz = cbuffer.data(pi_geom_11_off + 662 * ccomps * dcomps);

            auto g_z_y_z_xyzzzz = cbuffer.data(pi_geom_11_off + 663 * ccomps * dcomps);

            auto g_z_y_z_xzzzzz = cbuffer.data(pi_geom_11_off + 664 * ccomps * dcomps);

            auto g_z_y_z_yyyyyy = cbuffer.data(pi_geom_11_off + 665 * ccomps * dcomps);

            auto g_z_y_z_yyyyyz = cbuffer.data(pi_geom_11_off + 666 * ccomps * dcomps);

            auto g_z_y_z_yyyyzz = cbuffer.data(pi_geom_11_off + 667 * ccomps * dcomps);

            auto g_z_y_z_yyyzzz = cbuffer.data(pi_geom_11_off + 668 * ccomps * dcomps);

            auto g_z_y_z_yyzzzz = cbuffer.data(pi_geom_11_off + 669 * ccomps * dcomps);

            auto g_z_y_z_yzzzzz = cbuffer.data(pi_geom_11_off + 670 * ccomps * dcomps);

            auto g_z_y_z_zzzzzz = cbuffer.data(pi_geom_11_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_0_xxxxxx, g_0_y_0_xxxxxy, g_0_y_0_xxxxxz, g_0_y_0_xxxxyy, g_0_y_0_xxxxyz, g_0_y_0_xxxxzz, g_0_y_0_xxxyyy, g_0_y_0_xxxyyz, g_0_y_0_xxxyzz, g_0_y_0_xxxzzz, g_0_y_0_xxyyyy, g_0_y_0_xxyyyz, g_0_y_0_xxyyzz, g_0_y_0_xxyzzz, g_0_y_0_xxzzzz, g_0_y_0_xyyyyy, g_0_y_0_xyyyyz, g_0_y_0_xyyyzz, g_0_y_0_xyyzzz, g_0_y_0_xyzzzz, g_0_y_0_xzzzzz, g_0_y_0_yyyyyy, g_0_y_0_yyyyyz, g_0_y_0_yyyyzz, g_0_y_0_yyyzzz, g_0_y_0_yyzzzz, g_0_y_0_yzzzzz, g_0_y_0_zzzzzz, g_z_y_0_xxxxxx, g_z_y_0_xxxxxxz, g_z_y_0_xxxxxy, g_z_y_0_xxxxxyz, g_z_y_0_xxxxxz, g_z_y_0_xxxxxzz, g_z_y_0_xxxxyy, g_z_y_0_xxxxyyz, g_z_y_0_xxxxyz, g_z_y_0_xxxxyzz, g_z_y_0_xxxxzz, g_z_y_0_xxxxzzz, g_z_y_0_xxxyyy, g_z_y_0_xxxyyyz, g_z_y_0_xxxyyz, g_z_y_0_xxxyyzz, g_z_y_0_xxxyzz, g_z_y_0_xxxyzzz, g_z_y_0_xxxzzz, g_z_y_0_xxxzzzz, g_z_y_0_xxyyyy, g_z_y_0_xxyyyyz, g_z_y_0_xxyyyz, g_z_y_0_xxyyyzz, g_z_y_0_xxyyzz, g_z_y_0_xxyyzzz, g_z_y_0_xxyzzz, g_z_y_0_xxyzzzz, g_z_y_0_xxzzzz, g_z_y_0_xxzzzzz, g_z_y_0_xyyyyy, g_z_y_0_xyyyyyz, g_z_y_0_xyyyyz, g_z_y_0_xyyyyzz, g_z_y_0_xyyyzz, g_z_y_0_xyyyzzz, g_z_y_0_xyyzzz, g_z_y_0_xyyzzzz, g_z_y_0_xyzzzz, g_z_y_0_xyzzzzz, g_z_y_0_xzzzzz, g_z_y_0_xzzzzzz, g_z_y_0_yyyyyy, g_z_y_0_yyyyyyz, g_z_y_0_yyyyyz, g_z_y_0_yyyyyzz, g_z_y_0_yyyyzz, g_z_y_0_yyyyzzz, g_z_y_0_yyyzzz, g_z_y_0_yyyzzzz, g_z_y_0_yyzzzz, g_z_y_0_yyzzzzz, g_z_y_0_yzzzzz, g_z_y_0_yzzzzzz, g_z_y_0_zzzzzz, g_z_y_0_zzzzzzz, g_z_y_z_xxxxxx, g_z_y_z_xxxxxy, g_z_y_z_xxxxxz, g_z_y_z_xxxxyy, g_z_y_z_xxxxyz, g_z_y_z_xxxxzz, g_z_y_z_xxxyyy, g_z_y_z_xxxyyz, g_z_y_z_xxxyzz, g_z_y_z_xxxzzz, g_z_y_z_xxyyyy, g_z_y_z_xxyyyz, g_z_y_z_xxyyzz, g_z_y_z_xxyzzz, g_z_y_z_xxzzzz, g_z_y_z_xyyyyy, g_z_y_z_xyyyyz, g_z_y_z_xyyyzz, g_z_y_z_xyyzzz, g_z_y_z_xyzzzz, g_z_y_z_xzzzzz, g_z_y_z_yyyyyy, g_z_y_z_yyyyyz, g_z_y_z_yyyyzz, g_z_y_z_yyyzzz, g_z_y_z_yyzzzz, g_z_y_z_yzzzzz, g_z_y_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_z_xxxxxx[k] = -g_0_y_0_xxxxxx[k] - g_z_y_0_xxxxxx[k] * ab_z + g_z_y_0_xxxxxxz[k];

                g_z_y_z_xxxxxy[k] = -g_0_y_0_xxxxxy[k] - g_z_y_0_xxxxxy[k] * ab_z + g_z_y_0_xxxxxyz[k];

                g_z_y_z_xxxxxz[k] = -g_0_y_0_xxxxxz[k] - g_z_y_0_xxxxxz[k] * ab_z + g_z_y_0_xxxxxzz[k];

                g_z_y_z_xxxxyy[k] = -g_0_y_0_xxxxyy[k] - g_z_y_0_xxxxyy[k] * ab_z + g_z_y_0_xxxxyyz[k];

                g_z_y_z_xxxxyz[k] = -g_0_y_0_xxxxyz[k] - g_z_y_0_xxxxyz[k] * ab_z + g_z_y_0_xxxxyzz[k];

                g_z_y_z_xxxxzz[k] = -g_0_y_0_xxxxzz[k] - g_z_y_0_xxxxzz[k] * ab_z + g_z_y_0_xxxxzzz[k];

                g_z_y_z_xxxyyy[k] = -g_0_y_0_xxxyyy[k] - g_z_y_0_xxxyyy[k] * ab_z + g_z_y_0_xxxyyyz[k];

                g_z_y_z_xxxyyz[k] = -g_0_y_0_xxxyyz[k] - g_z_y_0_xxxyyz[k] * ab_z + g_z_y_0_xxxyyzz[k];

                g_z_y_z_xxxyzz[k] = -g_0_y_0_xxxyzz[k] - g_z_y_0_xxxyzz[k] * ab_z + g_z_y_0_xxxyzzz[k];

                g_z_y_z_xxxzzz[k] = -g_0_y_0_xxxzzz[k] - g_z_y_0_xxxzzz[k] * ab_z + g_z_y_0_xxxzzzz[k];

                g_z_y_z_xxyyyy[k] = -g_0_y_0_xxyyyy[k] - g_z_y_0_xxyyyy[k] * ab_z + g_z_y_0_xxyyyyz[k];

                g_z_y_z_xxyyyz[k] = -g_0_y_0_xxyyyz[k] - g_z_y_0_xxyyyz[k] * ab_z + g_z_y_0_xxyyyzz[k];

                g_z_y_z_xxyyzz[k] = -g_0_y_0_xxyyzz[k] - g_z_y_0_xxyyzz[k] * ab_z + g_z_y_0_xxyyzzz[k];

                g_z_y_z_xxyzzz[k] = -g_0_y_0_xxyzzz[k] - g_z_y_0_xxyzzz[k] * ab_z + g_z_y_0_xxyzzzz[k];

                g_z_y_z_xxzzzz[k] = -g_0_y_0_xxzzzz[k] - g_z_y_0_xxzzzz[k] * ab_z + g_z_y_0_xxzzzzz[k];

                g_z_y_z_xyyyyy[k] = -g_0_y_0_xyyyyy[k] - g_z_y_0_xyyyyy[k] * ab_z + g_z_y_0_xyyyyyz[k];

                g_z_y_z_xyyyyz[k] = -g_0_y_0_xyyyyz[k] - g_z_y_0_xyyyyz[k] * ab_z + g_z_y_0_xyyyyzz[k];

                g_z_y_z_xyyyzz[k] = -g_0_y_0_xyyyzz[k] - g_z_y_0_xyyyzz[k] * ab_z + g_z_y_0_xyyyzzz[k];

                g_z_y_z_xyyzzz[k] = -g_0_y_0_xyyzzz[k] - g_z_y_0_xyyzzz[k] * ab_z + g_z_y_0_xyyzzzz[k];

                g_z_y_z_xyzzzz[k] = -g_0_y_0_xyzzzz[k] - g_z_y_0_xyzzzz[k] * ab_z + g_z_y_0_xyzzzzz[k];

                g_z_y_z_xzzzzz[k] = -g_0_y_0_xzzzzz[k] - g_z_y_0_xzzzzz[k] * ab_z + g_z_y_0_xzzzzzz[k];

                g_z_y_z_yyyyyy[k] = -g_0_y_0_yyyyyy[k] - g_z_y_0_yyyyyy[k] * ab_z + g_z_y_0_yyyyyyz[k];

                g_z_y_z_yyyyyz[k] = -g_0_y_0_yyyyyz[k] - g_z_y_0_yyyyyz[k] * ab_z + g_z_y_0_yyyyyzz[k];

                g_z_y_z_yyyyzz[k] = -g_0_y_0_yyyyzz[k] - g_z_y_0_yyyyzz[k] * ab_z + g_z_y_0_yyyyzzz[k];

                g_z_y_z_yyyzzz[k] = -g_0_y_0_yyyzzz[k] - g_z_y_0_yyyzzz[k] * ab_z + g_z_y_0_yyyzzzz[k];

                g_z_y_z_yyzzzz[k] = -g_0_y_0_yyzzzz[k] - g_z_y_0_yyzzzz[k] * ab_z + g_z_y_0_yyzzzzz[k];

                g_z_y_z_yzzzzz[k] = -g_0_y_0_yzzzzz[k] - g_z_y_0_yzzzzz[k] * ab_z + g_z_y_0_yzzzzzz[k];

                g_z_y_z_zzzzzz[k] = -g_0_y_0_zzzzzz[k] - g_z_y_0_zzzzzz[k] * ab_z + g_z_y_0_zzzzzzz[k];
            }

            /// Set up 672-700 components of targeted buffer : cbuffer.data(

            auto g_z_z_x_xxxxxx = cbuffer.data(pi_geom_11_off + 672 * ccomps * dcomps);

            auto g_z_z_x_xxxxxy = cbuffer.data(pi_geom_11_off + 673 * ccomps * dcomps);

            auto g_z_z_x_xxxxxz = cbuffer.data(pi_geom_11_off + 674 * ccomps * dcomps);

            auto g_z_z_x_xxxxyy = cbuffer.data(pi_geom_11_off + 675 * ccomps * dcomps);

            auto g_z_z_x_xxxxyz = cbuffer.data(pi_geom_11_off + 676 * ccomps * dcomps);

            auto g_z_z_x_xxxxzz = cbuffer.data(pi_geom_11_off + 677 * ccomps * dcomps);

            auto g_z_z_x_xxxyyy = cbuffer.data(pi_geom_11_off + 678 * ccomps * dcomps);

            auto g_z_z_x_xxxyyz = cbuffer.data(pi_geom_11_off + 679 * ccomps * dcomps);

            auto g_z_z_x_xxxyzz = cbuffer.data(pi_geom_11_off + 680 * ccomps * dcomps);

            auto g_z_z_x_xxxzzz = cbuffer.data(pi_geom_11_off + 681 * ccomps * dcomps);

            auto g_z_z_x_xxyyyy = cbuffer.data(pi_geom_11_off + 682 * ccomps * dcomps);

            auto g_z_z_x_xxyyyz = cbuffer.data(pi_geom_11_off + 683 * ccomps * dcomps);

            auto g_z_z_x_xxyyzz = cbuffer.data(pi_geom_11_off + 684 * ccomps * dcomps);

            auto g_z_z_x_xxyzzz = cbuffer.data(pi_geom_11_off + 685 * ccomps * dcomps);

            auto g_z_z_x_xxzzzz = cbuffer.data(pi_geom_11_off + 686 * ccomps * dcomps);

            auto g_z_z_x_xyyyyy = cbuffer.data(pi_geom_11_off + 687 * ccomps * dcomps);

            auto g_z_z_x_xyyyyz = cbuffer.data(pi_geom_11_off + 688 * ccomps * dcomps);

            auto g_z_z_x_xyyyzz = cbuffer.data(pi_geom_11_off + 689 * ccomps * dcomps);

            auto g_z_z_x_xyyzzz = cbuffer.data(pi_geom_11_off + 690 * ccomps * dcomps);

            auto g_z_z_x_xyzzzz = cbuffer.data(pi_geom_11_off + 691 * ccomps * dcomps);

            auto g_z_z_x_xzzzzz = cbuffer.data(pi_geom_11_off + 692 * ccomps * dcomps);

            auto g_z_z_x_yyyyyy = cbuffer.data(pi_geom_11_off + 693 * ccomps * dcomps);

            auto g_z_z_x_yyyyyz = cbuffer.data(pi_geom_11_off + 694 * ccomps * dcomps);

            auto g_z_z_x_yyyyzz = cbuffer.data(pi_geom_11_off + 695 * ccomps * dcomps);

            auto g_z_z_x_yyyzzz = cbuffer.data(pi_geom_11_off + 696 * ccomps * dcomps);

            auto g_z_z_x_yyzzzz = cbuffer.data(pi_geom_11_off + 697 * ccomps * dcomps);

            auto g_z_z_x_yzzzzz = cbuffer.data(pi_geom_11_off + 698 * ccomps * dcomps);

            auto g_z_z_x_zzzzzz = cbuffer.data(pi_geom_11_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_0_xxxxxx, g_z_z_0_xxxxxxx, g_z_z_0_xxxxxxy, g_z_z_0_xxxxxxz, g_z_z_0_xxxxxy, g_z_z_0_xxxxxyy, g_z_z_0_xxxxxyz, g_z_z_0_xxxxxz, g_z_z_0_xxxxxzz, g_z_z_0_xxxxyy, g_z_z_0_xxxxyyy, g_z_z_0_xxxxyyz, g_z_z_0_xxxxyz, g_z_z_0_xxxxyzz, g_z_z_0_xxxxzz, g_z_z_0_xxxxzzz, g_z_z_0_xxxyyy, g_z_z_0_xxxyyyy, g_z_z_0_xxxyyyz, g_z_z_0_xxxyyz, g_z_z_0_xxxyyzz, g_z_z_0_xxxyzz, g_z_z_0_xxxyzzz, g_z_z_0_xxxzzz, g_z_z_0_xxxzzzz, g_z_z_0_xxyyyy, g_z_z_0_xxyyyyy, g_z_z_0_xxyyyyz, g_z_z_0_xxyyyz, g_z_z_0_xxyyyzz, g_z_z_0_xxyyzz, g_z_z_0_xxyyzzz, g_z_z_0_xxyzzz, g_z_z_0_xxyzzzz, g_z_z_0_xxzzzz, g_z_z_0_xxzzzzz, g_z_z_0_xyyyyy, g_z_z_0_xyyyyyy, g_z_z_0_xyyyyyz, g_z_z_0_xyyyyz, g_z_z_0_xyyyyzz, g_z_z_0_xyyyzz, g_z_z_0_xyyyzzz, g_z_z_0_xyyzzz, g_z_z_0_xyyzzzz, g_z_z_0_xyzzzz, g_z_z_0_xyzzzzz, g_z_z_0_xzzzzz, g_z_z_0_xzzzzzz, g_z_z_0_yyyyyy, g_z_z_0_yyyyyz, g_z_z_0_yyyyzz, g_z_z_0_yyyzzz, g_z_z_0_yyzzzz, g_z_z_0_yzzzzz, g_z_z_0_zzzzzz, g_z_z_x_xxxxxx, g_z_z_x_xxxxxy, g_z_z_x_xxxxxz, g_z_z_x_xxxxyy, g_z_z_x_xxxxyz, g_z_z_x_xxxxzz, g_z_z_x_xxxyyy, g_z_z_x_xxxyyz, g_z_z_x_xxxyzz, g_z_z_x_xxxzzz, g_z_z_x_xxyyyy, g_z_z_x_xxyyyz, g_z_z_x_xxyyzz, g_z_z_x_xxyzzz, g_z_z_x_xxzzzz, g_z_z_x_xyyyyy, g_z_z_x_xyyyyz, g_z_z_x_xyyyzz, g_z_z_x_xyyzzz, g_z_z_x_xyzzzz, g_z_z_x_xzzzzz, g_z_z_x_yyyyyy, g_z_z_x_yyyyyz, g_z_z_x_yyyyzz, g_z_z_x_yyyzzz, g_z_z_x_yyzzzz, g_z_z_x_yzzzzz, g_z_z_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_x_xxxxxx[k] = -g_z_z_0_xxxxxx[k] * ab_x + g_z_z_0_xxxxxxx[k];

                g_z_z_x_xxxxxy[k] = -g_z_z_0_xxxxxy[k] * ab_x + g_z_z_0_xxxxxxy[k];

                g_z_z_x_xxxxxz[k] = -g_z_z_0_xxxxxz[k] * ab_x + g_z_z_0_xxxxxxz[k];

                g_z_z_x_xxxxyy[k] = -g_z_z_0_xxxxyy[k] * ab_x + g_z_z_0_xxxxxyy[k];

                g_z_z_x_xxxxyz[k] = -g_z_z_0_xxxxyz[k] * ab_x + g_z_z_0_xxxxxyz[k];

                g_z_z_x_xxxxzz[k] = -g_z_z_0_xxxxzz[k] * ab_x + g_z_z_0_xxxxxzz[k];

                g_z_z_x_xxxyyy[k] = -g_z_z_0_xxxyyy[k] * ab_x + g_z_z_0_xxxxyyy[k];

                g_z_z_x_xxxyyz[k] = -g_z_z_0_xxxyyz[k] * ab_x + g_z_z_0_xxxxyyz[k];

                g_z_z_x_xxxyzz[k] = -g_z_z_0_xxxyzz[k] * ab_x + g_z_z_0_xxxxyzz[k];

                g_z_z_x_xxxzzz[k] = -g_z_z_0_xxxzzz[k] * ab_x + g_z_z_0_xxxxzzz[k];

                g_z_z_x_xxyyyy[k] = -g_z_z_0_xxyyyy[k] * ab_x + g_z_z_0_xxxyyyy[k];

                g_z_z_x_xxyyyz[k] = -g_z_z_0_xxyyyz[k] * ab_x + g_z_z_0_xxxyyyz[k];

                g_z_z_x_xxyyzz[k] = -g_z_z_0_xxyyzz[k] * ab_x + g_z_z_0_xxxyyzz[k];

                g_z_z_x_xxyzzz[k] = -g_z_z_0_xxyzzz[k] * ab_x + g_z_z_0_xxxyzzz[k];

                g_z_z_x_xxzzzz[k] = -g_z_z_0_xxzzzz[k] * ab_x + g_z_z_0_xxxzzzz[k];

                g_z_z_x_xyyyyy[k] = -g_z_z_0_xyyyyy[k] * ab_x + g_z_z_0_xxyyyyy[k];

                g_z_z_x_xyyyyz[k] = -g_z_z_0_xyyyyz[k] * ab_x + g_z_z_0_xxyyyyz[k];

                g_z_z_x_xyyyzz[k] = -g_z_z_0_xyyyzz[k] * ab_x + g_z_z_0_xxyyyzz[k];

                g_z_z_x_xyyzzz[k] = -g_z_z_0_xyyzzz[k] * ab_x + g_z_z_0_xxyyzzz[k];

                g_z_z_x_xyzzzz[k] = -g_z_z_0_xyzzzz[k] * ab_x + g_z_z_0_xxyzzzz[k];

                g_z_z_x_xzzzzz[k] = -g_z_z_0_xzzzzz[k] * ab_x + g_z_z_0_xxzzzzz[k];

                g_z_z_x_yyyyyy[k] = -g_z_z_0_yyyyyy[k] * ab_x + g_z_z_0_xyyyyyy[k];

                g_z_z_x_yyyyyz[k] = -g_z_z_0_yyyyyz[k] * ab_x + g_z_z_0_xyyyyyz[k];

                g_z_z_x_yyyyzz[k] = -g_z_z_0_yyyyzz[k] * ab_x + g_z_z_0_xyyyyzz[k];

                g_z_z_x_yyyzzz[k] = -g_z_z_0_yyyzzz[k] * ab_x + g_z_z_0_xyyyzzz[k];

                g_z_z_x_yyzzzz[k] = -g_z_z_0_yyzzzz[k] * ab_x + g_z_z_0_xyyzzzz[k];

                g_z_z_x_yzzzzz[k] = -g_z_z_0_yzzzzz[k] * ab_x + g_z_z_0_xyzzzzz[k];

                g_z_z_x_zzzzzz[k] = -g_z_z_0_zzzzzz[k] * ab_x + g_z_z_0_xzzzzzz[k];
            }

            /// Set up 700-728 components of targeted buffer : cbuffer.data(

            auto g_z_z_y_xxxxxx = cbuffer.data(pi_geom_11_off + 700 * ccomps * dcomps);

            auto g_z_z_y_xxxxxy = cbuffer.data(pi_geom_11_off + 701 * ccomps * dcomps);

            auto g_z_z_y_xxxxxz = cbuffer.data(pi_geom_11_off + 702 * ccomps * dcomps);

            auto g_z_z_y_xxxxyy = cbuffer.data(pi_geom_11_off + 703 * ccomps * dcomps);

            auto g_z_z_y_xxxxyz = cbuffer.data(pi_geom_11_off + 704 * ccomps * dcomps);

            auto g_z_z_y_xxxxzz = cbuffer.data(pi_geom_11_off + 705 * ccomps * dcomps);

            auto g_z_z_y_xxxyyy = cbuffer.data(pi_geom_11_off + 706 * ccomps * dcomps);

            auto g_z_z_y_xxxyyz = cbuffer.data(pi_geom_11_off + 707 * ccomps * dcomps);

            auto g_z_z_y_xxxyzz = cbuffer.data(pi_geom_11_off + 708 * ccomps * dcomps);

            auto g_z_z_y_xxxzzz = cbuffer.data(pi_geom_11_off + 709 * ccomps * dcomps);

            auto g_z_z_y_xxyyyy = cbuffer.data(pi_geom_11_off + 710 * ccomps * dcomps);

            auto g_z_z_y_xxyyyz = cbuffer.data(pi_geom_11_off + 711 * ccomps * dcomps);

            auto g_z_z_y_xxyyzz = cbuffer.data(pi_geom_11_off + 712 * ccomps * dcomps);

            auto g_z_z_y_xxyzzz = cbuffer.data(pi_geom_11_off + 713 * ccomps * dcomps);

            auto g_z_z_y_xxzzzz = cbuffer.data(pi_geom_11_off + 714 * ccomps * dcomps);

            auto g_z_z_y_xyyyyy = cbuffer.data(pi_geom_11_off + 715 * ccomps * dcomps);

            auto g_z_z_y_xyyyyz = cbuffer.data(pi_geom_11_off + 716 * ccomps * dcomps);

            auto g_z_z_y_xyyyzz = cbuffer.data(pi_geom_11_off + 717 * ccomps * dcomps);

            auto g_z_z_y_xyyzzz = cbuffer.data(pi_geom_11_off + 718 * ccomps * dcomps);

            auto g_z_z_y_xyzzzz = cbuffer.data(pi_geom_11_off + 719 * ccomps * dcomps);

            auto g_z_z_y_xzzzzz = cbuffer.data(pi_geom_11_off + 720 * ccomps * dcomps);

            auto g_z_z_y_yyyyyy = cbuffer.data(pi_geom_11_off + 721 * ccomps * dcomps);

            auto g_z_z_y_yyyyyz = cbuffer.data(pi_geom_11_off + 722 * ccomps * dcomps);

            auto g_z_z_y_yyyyzz = cbuffer.data(pi_geom_11_off + 723 * ccomps * dcomps);

            auto g_z_z_y_yyyzzz = cbuffer.data(pi_geom_11_off + 724 * ccomps * dcomps);

            auto g_z_z_y_yyzzzz = cbuffer.data(pi_geom_11_off + 725 * ccomps * dcomps);

            auto g_z_z_y_yzzzzz = cbuffer.data(pi_geom_11_off + 726 * ccomps * dcomps);

            auto g_z_z_y_zzzzzz = cbuffer.data(pi_geom_11_off + 727 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_0_xxxxxx, g_z_z_0_xxxxxxy, g_z_z_0_xxxxxy, g_z_z_0_xxxxxyy, g_z_z_0_xxxxxyz, g_z_z_0_xxxxxz, g_z_z_0_xxxxyy, g_z_z_0_xxxxyyy, g_z_z_0_xxxxyyz, g_z_z_0_xxxxyz, g_z_z_0_xxxxyzz, g_z_z_0_xxxxzz, g_z_z_0_xxxyyy, g_z_z_0_xxxyyyy, g_z_z_0_xxxyyyz, g_z_z_0_xxxyyz, g_z_z_0_xxxyyzz, g_z_z_0_xxxyzz, g_z_z_0_xxxyzzz, g_z_z_0_xxxzzz, g_z_z_0_xxyyyy, g_z_z_0_xxyyyyy, g_z_z_0_xxyyyyz, g_z_z_0_xxyyyz, g_z_z_0_xxyyyzz, g_z_z_0_xxyyzz, g_z_z_0_xxyyzzz, g_z_z_0_xxyzzz, g_z_z_0_xxyzzzz, g_z_z_0_xxzzzz, g_z_z_0_xyyyyy, g_z_z_0_xyyyyyy, g_z_z_0_xyyyyyz, g_z_z_0_xyyyyz, g_z_z_0_xyyyyzz, g_z_z_0_xyyyzz, g_z_z_0_xyyyzzz, g_z_z_0_xyyzzz, g_z_z_0_xyyzzzz, g_z_z_0_xyzzzz, g_z_z_0_xyzzzzz, g_z_z_0_xzzzzz, g_z_z_0_yyyyyy, g_z_z_0_yyyyyyy, g_z_z_0_yyyyyyz, g_z_z_0_yyyyyz, g_z_z_0_yyyyyzz, g_z_z_0_yyyyzz, g_z_z_0_yyyyzzz, g_z_z_0_yyyzzz, g_z_z_0_yyyzzzz, g_z_z_0_yyzzzz, g_z_z_0_yyzzzzz, g_z_z_0_yzzzzz, g_z_z_0_yzzzzzz, g_z_z_0_zzzzzz, g_z_z_y_xxxxxx, g_z_z_y_xxxxxy, g_z_z_y_xxxxxz, g_z_z_y_xxxxyy, g_z_z_y_xxxxyz, g_z_z_y_xxxxzz, g_z_z_y_xxxyyy, g_z_z_y_xxxyyz, g_z_z_y_xxxyzz, g_z_z_y_xxxzzz, g_z_z_y_xxyyyy, g_z_z_y_xxyyyz, g_z_z_y_xxyyzz, g_z_z_y_xxyzzz, g_z_z_y_xxzzzz, g_z_z_y_xyyyyy, g_z_z_y_xyyyyz, g_z_z_y_xyyyzz, g_z_z_y_xyyzzz, g_z_z_y_xyzzzz, g_z_z_y_xzzzzz, g_z_z_y_yyyyyy, g_z_z_y_yyyyyz, g_z_z_y_yyyyzz, g_z_z_y_yyyzzz, g_z_z_y_yyzzzz, g_z_z_y_yzzzzz, g_z_z_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_y_xxxxxx[k] = -g_z_z_0_xxxxxx[k] * ab_y + g_z_z_0_xxxxxxy[k];

                g_z_z_y_xxxxxy[k] = -g_z_z_0_xxxxxy[k] * ab_y + g_z_z_0_xxxxxyy[k];

                g_z_z_y_xxxxxz[k] = -g_z_z_0_xxxxxz[k] * ab_y + g_z_z_0_xxxxxyz[k];

                g_z_z_y_xxxxyy[k] = -g_z_z_0_xxxxyy[k] * ab_y + g_z_z_0_xxxxyyy[k];

                g_z_z_y_xxxxyz[k] = -g_z_z_0_xxxxyz[k] * ab_y + g_z_z_0_xxxxyyz[k];

                g_z_z_y_xxxxzz[k] = -g_z_z_0_xxxxzz[k] * ab_y + g_z_z_0_xxxxyzz[k];

                g_z_z_y_xxxyyy[k] = -g_z_z_0_xxxyyy[k] * ab_y + g_z_z_0_xxxyyyy[k];

                g_z_z_y_xxxyyz[k] = -g_z_z_0_xxxyyz[k] * ab_y + g_z_z_0_xxxyyyz[k];

                g_z_z_y_xxxyzz[k] = -g_z_z_0_xxxyzz[k] * ab_y + g_z_z_0_xxxyyzz[k];

                g_z_z_y_xxxzzz[k] = -g_z_z_0_xxxzzz[k] * ab_y + g_z_z_0_xxxyzzz[k];

                g_z_z_y_xxyyyy[k] = -g_z_z_0_xxyyyy[k] * ab_y + g_z_z_0_xxyyyyy[k];

                g_z_z_y_xxyyyz[k] = -g_z_z_0_xxyyyz[k] * ab_y + g_z_z_0_xxyyyyz[k];

                g_z_z_y_xxyyzz[k] = -g_z_z_0_xxyyzz[k] * ab_y + g_z_z_0_xxyyyzz[k];

                g_z_z_y_xxyzzz[k] = -g_z_z_0_xxyzzz[k] * ab_y + g_z_z_0_xxyyzzz[k];

                g_z_z_y_xxzzzz[k] = -g_z_z_0_xxzzzz[k] * ab_y + g_z_z_0_xxyzzzz[k];

                g_z_z_y_xyyyyy[k] = -g_z_z_0_xyyyyy[k] * ab_y + g_z_z_0_xyyyyyy[k];

                g_z_z_y_xyyyyz[k] = -g_z_z_0_xyyyyz[k] * ab_y + g_z_z_0_xyyyyyz[k];

                g_z_z_y_xyyyzz[k] = -g_z_z_0_xyyyzz[k] * ab_y + g_z_z_0_xyyyyzz[k];

                g_z_z_y_xyyzzz[k] = -g_z_z_0_xyyzzz[k] * ab_y + g_z_z_0_xyyyzzz[k];

                g_z_z_y_xyzzzz[k] = -g_z_z_0_xyzzzz[k] * ab_y + g_z_z_0_xyyzzzz[k];

                g_z_z_y_xzzzzz[k] = -g_z_z_0_xzzzzz[k] * ab_y + g_z_z_0_xyzzzzz[k];

                g_z_z_y_yyyyyy[k] = -g_z_z_0_yyyyyy[k] * ab_y + g_z_z_0_yyyyyyy[k];

                g_z_z_y_yyyyyz[k] = -g_z_z_0_yyyyyz[k] * ab_y + g_z_z_0_yyyyyyz[k];

                g_z_z_y_yyyyzz[k] = -g_z_z_0_yyyyzz[k] * ab_y + g_z_z_0_yyyyyzz[k];

                g_z_z_y_yyyzzz[k] = -g_z_z_0_yyyzzz[k] * ab_y + g_z_z_0_yyyyzzz[k];

                g_z_z_y_yyzzzz[k] = -g_z_z_0_yyzzzz[k] * ab_y + g_z_z_0_yyyzzzz[k];

                g_z_z_y_yzzzzz[k] = -g_z_z_0_yzzzzz[k] * ab_y + g_z_z_0_yyzzzzz[k];

                g_z_z_y_zzzzzz[k] = -g_z_z_0_zzzzzz[k] * ab_y + g_z_z_0_yzzzzzz[k];
            }

            /// Set up 728-756 components of targeted buffer : cbuffer.data(

            auto g_z_z_z_xxxxxx = cbuffer.data(pi_geom_11_off + 728 * ccomps * dcomps);

            auto g_z_z_z_xxxxxy = cbuffer.data(pi_geom_11_off + 729 * ccomps * dcomps);

            auto g_z_z_z_xxxxxz = cbuffer.data(pi_geom_11_off + 730 * ccomps * dcomps);

            auto g_z_z_z_xxxxyy = cbuffer.data(pi_geom_11_off + 731 * ccomps * dcomps);

            auto g_z_z_z_xxxxyz = cbuffer.data(pi_geom_11_off + 732 * ccomps * dcomps);

            auto g_z_z_z_xxxxzz = cbuffer.data(pi_geom_11_off + 733 * ccomps * dcomps);

            auto g_z_z_z_xxxyyy = cbuffer.data(pi_geom_11_off + 734 * ccomps * dcomps);

            auto g_z_z_z_xxxyyz = cbuffer.data(pi_geom_11_off + 735 * ccomps * dcomps);

            auto g_z_z_z_xxxyzz = cbuffer.data(pi_geom_11_off + 736 * ccomps * dcomps);

            auto g_z_z_z_xxxzzz = cbuffer.data(pi_geom_11_off + 737 * ccomps * dcomps);

            auto g_z_z_z_xxyyyy = cbuffer.data(pi_geom_11_off + 738 * ccomps * dcomps);

            auto g_z_z_z_xxyyyz = cbuffer.data(pi_geom_11_off + 739 * ccomps * dcomps);

            auto g_z_z_z_xxyyzz = cbuffer.data(pi_geom_11_off + 740 * ccomps * dcomps);

            auto g_z_z_z_xxyzzz = cbuffer.data(pi_geom_11_off + 741 * ccomps * dcomps);

            auto g_z_z_z_xxzzzz = cbuffer.data(pi_geom_11_off + 742 * ccomps * dcomps);

            auto g_z_z_z_xyyyyy = cbuffer.data(pi_geom_11_off + 743 * ccomps * dcomps);

            auto g_z_z_z_xyyyyz = cbuffer.data(pi_geom_11_off + 744 * ccomps * dcomps);

            auto g_z_z_z_xyyyzz = cbuffer.data(pi_geom_11_off + 745 * ccomps * dcomps);

            auto g_z_z_z_xyyzzz = cbuffer.data(pi_geom_11_off + 746 * ccomps * dcomps);

            auto g_z_z_z_xyzzzz = cbuffer.data(pi_geom_11_off + 747 * ccomps * dcomps);

            auto g_z_z_z_xzzzzz = cbuffer.data(pi_geom_11_off + 748 * ccomps * dcomps);

            auto g_z_z_z_yyyyyy = cbuffer.data(pi_geom_11_off + 749 * ccomps * dcomps);

            auto g_z_z_z_yyyyyz = cbuffer.data(pi_geom_11_off + 750 * ccomps * dcomps);

            auto g_z_z_z_yyyyzz = cbuffer.data(pi_geom_11_off + 751 * ccomps * dcomps);

            auto g_z_z_z_yyyzzz = cbuffer.data(pi_geom_11_off + 752 * ccomps * dcomps);

            auto g_z_z_z_yyzzzz = cbuffer.data(pi_geom_11_off + 753 * ccomps * dcomps);

            auto g_z_z_z_yzzzzz = cbuffer.data(pi_geom_11_off + 754 * ccomps * dcomps);

            auto g_z_z_z_zzzzzz = cbuffer.data(pi_geom_11_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_0_xxxxxx, g_0_z_0_xxxxxy, g_0_z_0_xxxxxz, g_0_z_0_xxxxyy, g_0_z_0_xxxxyz, g_0_z_0_xxxxzz, g_0_z_0_xxxyyy, g_0_z_0_xxxyyz, g_0_z_0_xxxyzz, g_0_z_0_xxxzzz, g_0_z_0_xxyyyy, g_0_z_0_xxyyyz, g_0_z_0_xxyyzz, g_0_z_0_xxyzzz, g_0_z_0_xxzzzz, g_0_z_0_xyyyyy, g_0_z_0_xyyyyz, g_0_z_0_xyyyzz, g_0_z_0_xyyzzz, g_0_z_0_xyzzzz, g_0_z_0_xzzzzz, g_0_z_0_yyyyyy, g_0_z_0_yyyyyz, g_0_z_0_yyyyzz, g_0_z_0_yyyzzz, g_0_z_0_yyzzzz, g_0_z_0_yzzzzz, g_0_z_0_zzzzzz, g_z_0_0_xxxxxx, g_z_0_0_xxxxxy, g_z_0_0_xxxxxz, g_z_0_0_xxxxyy, g_z_0_0_xxxxyz, g_z_0_0_xxxxzz, g_z_0_0_xxxyyy, g_z_0_0_xxxyyz, g_z_0_0_xxxyzz, g_z_0_0_xxxzzz, g_z_0_0_xxyyyy, g_z_0_0_xxyyyz, g_z_0_0_xxyyzz, g_z_0_0_xxyzzz, g_z_0_0_xxzzzz, g_z_0_0_xyyyyy, g_z_0_0_xyyyyz, g_z_0_0_xyyyzz, g_z_0_0_xyyzzz, g_z_0_0_xyzzzz, g_z_0_0_xzzzzz, g_z_0_0_yyyyyy, g_z_0_0_yyyyyz, g_z_0_0_yyyyzz, g_z_0_0_yyyzzz, g_z_0_0_yyzzzz, g_z_0_0_yzzzzz, g_z_0_0_zzzzzz, g_z_z_0_xxxxxx, g_z_z_0_xxxxxxz, g_z_z_0_xxxxxy, g_z_z_0_xxxxxyz, g_z_z_0_xxxxxz, g_z_z_0_xxxxxzz, g_z_z_0_xxxxyy, g_z_z_0_xxxxyyz, g_z_z_0_xxxxyz, g_z_z_0_xxxxyzz, g_z_z_0_xxxxzz, g_z_z_0_xxxxzzz, g_z_z_0_xxxyyy, g_z_z_0_xxxyyyz, g_z_z_0_xxxyyz, g_z_z_0_xxxyyzz, g_z_z_0_xxxyzz, g_z_z_0_xxxyzzz, g_z_z_0_xxxzzz, g_z_z_0_xxxzzzz, g_z_z_0_xxyyyy, g_z_z_0_xxyyyyz, g_z_z_0_xxyyyz, g_z_z_0_xxyyyzz, g_z_z_0_xxyyzz, g_z_z_0_xxyyzzz, g_z_z_0_xxyzzz, g_z_z_0_xxyzzzz, g_z_z_0_xxzzzz, g_z_z_0_xxzzzzz, g_z_z_0_xyyyyy, g_z_z_0_xyyyyyz, g_z_z_0_xyyyyz, g_z_z_0_xyyyyzz, g_z_z_0_xyyyzz, g_z_z_0_xyyyzzz, g_z_z_0_xyyzzz, g_z_z_0_xyyzzzz, g_z_z_0_xyzzzz, g_z_z_0_xyzzzzz, g_z_z_0_xzzzzz, g_z_z_0_xzzzzzz, g_z_z_0_yyyyyy, g_z_z_0_yyyyyyz, g_z_z_0_yyyyyz, g_z_z_0_yyyyyzz, g_z_z_0_yyyyzz, g_z_z_0_yyyyzzz, g_z_z_0_yyyzzz, g_z_z_0_yyyzzzz, g_z_z_0_yyzzzz, g_z_z_0_yyzzzzz, g_z_z_0_yzzzzz, g_z_z_0_yzzzzzz, g_z_z_0_zzzzzz, g_z_z_0_zzzzzzz, g_z_z_z_xxxxxx, g_z_z_z_xxxxxy, g_z_z_z_xxxxxz, g_z_z_z_xxxxyy, g_z_z_z_xxxxyz, g_z_z_z_xxxxzz, g_z_z_z_xxxyyy, g_z_z_z_xxxyyz, g_z_z_z_xxxyzz, g_z_z_z_xxxzzz, g_z_z_z_xxyyyy, g_z_z_z_xxyyyz, g_z_z_z_xxyyzz, g_z_z_z_xxyzzz, g_z_z_z_xxzzzz, g_z_z_z_xyyyyy, g_z_z_z_xyyyyz, g_z_z_z_xyyyzz, g_z_z_z_xyyzzz, g_z_z_z_xyzzzz, g_z_z_z_xzzzzz, g_z_z_z_yyyyyy, g_z_z_z_yyyyyz, g_z_z_z_yyyyzz, g_z_z_z_yyyzzz, g_z_z_z_yyzzzz, g_z_z_z_yzzzzz, g_z_z_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_z_xxxxxx[k] = -g_0_z_0_xxxxxx[k] + g_z_0_0_xxxxxx[k] - g_z_z_0_xxxxxx[k] * ab_z + g_z_z_0_xxxxxxz[k];

                g_z_z_z_xxxxxy[k] = -g_0_z_0_xxxxxy[k] + g_z_0_0_xxxxxy[k] - g_z_z_0_xxxxxy[k] * ab_z + g_z_z_0_xxxxxyz[k];

                g_z_z_z_xxxxxz[k] = -g_0_z_0_xxxxxz[k] + g_z_0_0_xxxxxz[k] - g_z_z_0_xxxxxz[k] * ab_z + g_z_z_0_xxxxxzz[k];

                g_z_z_z_xxxxyy[k] = -g_0_z_0_xxxxyy[k] + g_z_0_0_xxxxyy[k] - g_z_z_0_xxxxyy[k] * ab_z + g_z_z_0_xxxxyyz[k];

                g_z_z_z_xxxxyz[k] = -g_0_z_0_xxxxyz[k] + g_z_0_0_xxxxyz[k] - g_z_z_0_xxxxyz[k] * ab_z + g_z_z_0_xxxxyzz[k];

                g_z_z_z_xxxxzz[k] = -g_0_z_0_xxxxzz[k] + g_z_0_0_xxxxzz[k] - g_z_z_0_xxxxzz[k] * ab_z + g_z_z_0_xxxxzzz[k];

                g_z_z_z_xxxyyy[k] = -g_0_z_0_xxxyyy[k] + g_z_0_0_xxxyyy[k] - g_z_z_0_xxxyyy[k] * ab_z + g_z_z_0_xxxyyyz[k];

                g_z_z_z_xxxyyz[k] = -g_0_z_0_xxxyyz[k] + g_z_0_0_xxxyyz[k] - g_z_z_0_xxxyyz[k] * ab_z + g_z_z_0_xxxyyzz[k];

                g_z_z_z_xxxyzz[k] = -g_0_z_0_xxxyzz[k] + g_z_0_0_xxxyzz[k] - g_z_z_0_xxxyzz[k] * ab_z + g_z_z_0_xxxyzzz[k];

                g_z_z_z_xxxzzz[k] = -g_0_z_0_xxxzzz[k] + g_z_0_0_xxxzzz[k] - g_z_z_0_xxxzzz[k] * ab_z + g_z_z_0_xxxzzzz[k];

                g_z_z_z_xxyyyy[k] = -g_0_z_0_xxyyyy[k] + g_z_0_0_xxyyyy[k] - g_z_z_0_xxyyyy[k] * ab_z + g_z_z_0_xxyyyyz[k];

                g_z_z_z_xxyyyz[k] = -g_0_z_0_xxyyyz[k] + g_z_0_0_xxyyyz[k] - g_z_z_0_xxyyyz[k] * ab_z + g_z_z_0_xxyyyzz[k];

                g_z_z_z_xxyyzz[k] = -g_0_z_0_xxyyzz[k] + g_z_0_0_xxyyzz[k] - g_z_z_0_xxyyzz[k] * ab_z + g_z_z_0_xxyyzzz[k];

                g_z_z_z_xxyzzz[k] = -g_0_z_0_xxyzzz[k] + g_z_0_0_xxyzzz[k] - g_z_z_0_xxyzzz[k] * ab_z + g_z_z_0_xxyzzzz[k];

                g_z_z_z_xxzzzz[k] = -g_0_z_0_xxzzzz[k] + g_z_0_0_xxzzzz[k] - g_z_z_0_xxzzzz[k] * ab_z + g_z_z_0_xxzzzzz[k];

                g_z_z_z_xyyyyy[k] = -g_0_z_0_xyyyyy[k] + g_z_0_0_xyyyyy[k] - g_z_z_0_xyyyyy[k] * ab_z + g_z_z_0_xyyyyyz[k];

                g_z_z_z_xyyyyz[k] = -g_0_z_0_xyyyyz[k] + g_z_0_0_xyyyyz[k] - g_z_z_0_xyyyyz[k] * ab_z + g_z_z_0_xyyyyzz[k];

                g_z_z_z_xyyyzz[k] = -g_0_z_0_xyyyzz[k] + g_z_0_0_xyyyzz[k] - g_z_z_0_xyyyzz[k] * ab_z + g_z_z_0_xyyyzzz[k];

                g_z_z_z_xyyzzz[k] = -g_0_z_0_xyyzzz[k] + g_z_0_0_xyyzzz[k] - g_z_z_0_xyyzzz[k] * ab_z + g_z_z_0_xyyzzzz[k];

                g_z_z_z_xyzzzz[k] = -g_0_z_0_xyzzzz[k] + g_z_0_0_xyzzzz[k] - g_z_z_0_xyzzzz[k] * ab_z + g_z_z_0_xyzzzzz[k];

                g_z_z_z_xzzzzz[k] = -g_0_z_0_xzzzzz[k] + g_z_0_0_xzzzzz[k] - g_z_z_0_xzzzzz[k] * ab_z + g_z_z_0_xzzzzzz[k];

                g_z_z_z_yyyyyy[k] = -g_0_z_0_yyyyyy[k] + g_z_0_0_yyyyyy[k] - g_z_z_0_yyyyyy[k] * ab_z + g_z_z_0_yyyyyyz[k];

                g_z_z_z_yyyyyz[k] = -g_0_z_0_yyyyyz[k] + g_z_0_0_yyyyyz[k] - g_z_z_0_yyyyyz[k] * ab_z + g_z_z_0_yyyyyzz[k];

                g_z_z_z_yyyyzz[k] = -g_0_z_0_yyyyzz[k] + g_z_0_0_yyyyzz[k] - g_z_z_0_yyyyzz[k] * ab_z + g_z_z_0_yyyyzzz[k];

                g_z_z_z_yyyzzz[k] = -g_0_z_0_yyyzzz[k] + g_z_0_0_yyyzzz[k] - g_z_z_0_yyyzzz[k] * ab_z + g_z_z_0_yyyzzzz[k];

                g_z_z_z_yyzzzz[k] = -g_0_z_0_yyzzzz[k] + g_z_0_0_yyzzzz[k] - g_z_z_0_yyzzzz[k] * ab_z + g_z_z_0_yyzzzzz[k];

                g_z_z_z_yzzzzz[k] = -g_0_z_0_yzzzzz[k] + g_z_0_0_yzzzzz[k] - g_z_z_0_yzzzzz[k] * ab_z + g_z_z_0_yzzzzzz[k];

                g_z_z_z_zzzzzz[k] = -g_0_z_0_zzzzzz[k] + g_z_0_0_zzzzzz[k] - g_z_z_0_zzzzzz[k] * ab_z + g_z_z_0_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

