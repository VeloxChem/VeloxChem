#include "ElectronRepulsionGeom1010ContrRecPIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom1010_hrr_electron_repulsion_pixx(CSimdArray<double>& cbuffer,
                                              const size_t idx_geom_1010_pixx,
                                              const size_t idx_geom_0010_sixx,
                                              const size_t idx_geom_1010_sixx,
                                              const size_t idx_geom_1010_skxx,
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

            const auto si_geom_0010_off = idx_geom_0010_sixx + i * dcomps + j;

            auto g_0_0_x_0_0_xxxxxx = cbuffer.data(si_geom_0010_off + 0 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxxxy = cbuffer.data(si_geom_0010_off + 1 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxxxz = cbuffer.data(si_geom_0010_off + 2 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxxyy = cbuffer.data(si_geom_0010_off + 3 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxxyz = cbuffer.data(si_geom_0010_off + 4 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxxzz = cbuffer.data(si_geom_0010_off + 5 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxyyy = cbuffer.data(si_geom_0010_off + 6 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxyyz = cbuffer.data(si_geom_0010_off + 7 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxyzz = cbuffer.data(si_geom_0010_off + 8 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxxzzz = cbuffer.data(si_geom_0010_off + 9 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxyyyy = cbuffer.data(si_geom_0010_off + 10 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxyyyz = cbuffer.data(si_geom_0010_off + 11 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxyyzz = cbuffer.data(si_geom_0010_off + 12 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxyzzz = cbuffer.data(si_geom_0010_off + 13 * ccomps * dcomps);

            auto g_0_0_x_0_0_xxzzzz = cbuffer.data(si_geom_0010_off + 14 * ccomps * dcomps);

            auto g_0_0_x_0_0_xyyyyy = cbuffer.data(si_geom_0010_off + 15 * ccomps * dcomps);

            auto g_0_0_x_0_0_xyyyyz = cbuffer.data(si_geom_0010_off + 16 * ccomps * dcomps);

            auto g_0_0_x_0_0_xyyyzz = cbuffer.data(si_geom_0010_off + 17 * ccomps * dcomps);

            auto g_0_0_x_0_0_xyyzzz = cbuffer.data(si_geom_0010_off + 18 * ccomps * dcomps);

            auto g_0_0_x_0_0_xyzzzz = cbuffer.data(si_geom_0010_off + 19 * ccomps * dcomps);

            auto g_0_0_x_0_0_xzzzzz = cbuffer.data(si_geom_0010_off + 20 * ccomps * dcomps);

            auto g_0_0_x_0_0_yyyyyy = cbuffer.data(si_geom_0010_off + 21 * ccomps * dcomps);

            auto g_0_0_x_0_0_yyyyyz = cbuffer.data(si_geom_0010_off + 22 * ccomps * dcomps);

            auto g_0_0_x_0_0_yyyyzz = cbuffer.data(si_geom_0010_off + 23 * ccomps * dcomps);

            auto g_0_0_x_0_0_yyyzzz = cbuffer.data(si_geom_0010_off + 24 * ccomps * dcomps);

            auto g_0_0_x_0_0_yyzzzz = cbuffer.data(si_geom_0010_off + 25 * ccomps * dcomps);

            auto g_0_0_x_0_0_yzzzzz = cbuffer.data(si_geom_0010_off + 26 * ccomps * dcomps);

            auto g_0_0_x_0_0_zzzzzz = cbuffer.data(si_geom_0010_off + 27 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxxxx = cbuffer.data(si_geom_0010_off + 28 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxxxy = cbuffer.data(si_geom_0010_off + 29 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxxxz = cbuffer.data(si_geom_0010_off + 30 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxxyy = cbuffer.data(si_geom_0010_off + 31 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxxyz = cbuffer.data(si_geom_0010_off + 32 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxxzz = cbuffer.data(si_geom_0010_off + 33 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxyyy = cbuffer.data(si_geom_0010_off + 34 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxyyz = cbuffer.data(si_geom_0010_off + 35 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxyzz = cbuffer.data(si_geom_0010_off + 36 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxxzzz = cbuffer.data(si_geom_0010_off + 37 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxyyyy = cbuffer.data(si_geom_0010_off + 38 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxyyyz = cbuffer.data(si_geom_0010_off + 39 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxyyzz = cbuffer.data(si_geom_0010_off + 40 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxyzzz = cbuffer.data(si_geom_0010_off + 41 * ccomps * dcomps);

            auto g_0_0_y_0_0_xxzzzz = cbuffer.data(si_geom_0010_off + 42 * ccomps * dcomps);

            auto g_0_0_y_0_0_xyyyyy = cbuffer.data(si_geom_0010_off + 43 * ccomps * dcomps);

            auto g_0_0_y_0_0_xyyyyz = cbuffer.data(si_geom_0010_off + 44 * ccomps * dcomps);

            auto g_0_0_y_0_0_xyyyzz = cbuffer.data(si_geom_0010_off + 45 * ccomps * dcomps);

            auto g_0_0_y_0_0_xyyzzz = cbuffer.data(si_geom_0010_off + 46 * ccomps * dcomps);

            auto g_0_0_y_0_0_xyzzzz = cbuffer.data(si_geom_0010_off + 47 * ccomps * dcomps);

            auto g_0_0_y_0_0_xzzzzz = cbuffer.data(si_geom_0010_off + 48 * ccomps * dcomps);

            auto g_0_0_y_0_0_yyyyyy = cbuffer.data(si_geom_0010_off + 49 * ccomps * dcomps);

            auto g_0_0_y_0_0_yyyyyz = cbuffer.data(si_geom_0010_off + 50 * ccomps * dcomps);

            auto g_0_0_y_0_0_yyyyzz = cbuffer.data(si_geom_0010_off + 51 * ccomps * dcomps);

            auto g_0_0_y_0_0_yyyzzz = cbuffer.data(si_geom_0010_off + 52 * ccomps * dcomps);

            auto g_0_0_y_0_0_yyzzzz = cbuffer.data(si_geom_0010_off + 53 * ccomps * dcomps);

            auto g_0_0_y_0_0_yzzzzz = cbuffer.data(si_geom_0010_off + 54 * ccomps * dcomps);

            auto g_0_0_y_0_0_zzzzzz = cbuffer.data(si_geom_0010_off + 55 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxxxx = cbuffer.data(si_geom_0010_off + 56 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxxxy = cbuffer.data(si_geom_0010_off + 57 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxxxz = cbuffer.data(si_geom_0010_off + 58 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxxyy = cbuffer.data(si_geom_0010_off + 59 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxxyz = cbuffer.data(si_geom_0010_off + 60 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxxzz = cbuffer.data(si_geom_0010_off + 61 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxyyy = cbuffer.data(si_geom_0010_off + 62 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxyyz = cbuffer.data(si_geom_0010_off + 63 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxyzz = cbuffer.data(si_geom_0010_off + 64 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxxzzz = cbuffer.data(si_geom_0010_off + 65 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxyyyy = cbuffer.data(si_geom_0010_off + 66 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxyyyz = cbuffer.data(si_geom_0010_off + 67 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxyyzz = cbuffer.data(si_geom_0010_off + 68 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxyzzz = cbuffer.data(si_geom_0010_off + 69 * ccomps * dcomps);

            auto g_0_0_z_0_0_xxzzzz = cbuffer.data(si_geom_0010_off + 70 * ccomps * dcomps);

            auto g_0_0_z_0_0_xyyyyy = cbuffer.data(si_geom_0010_off + 71 * ccomps * dcomps);

            auto g_0_0_z_0_0_xyyyyz = cbuffer.data(si_geom_0010_off + 72 * ccomps * dcomps);

            auto g_0_0_z_0_0_xyyyzz = cbuffer.data(si_geom_0010_off + 73 * ccomps * dcomps);

            auto g_0_0_z_0_0_xyyzzz = cbuffer.data(si_geom_0010_off + 74 * ccomps * dcomps);

            auto g_0_0_z_0_0_xyzzzz = cbuffer.data(si_geom_0010_off + 75 * ccomps * dcomps);

            auto g_0_0_z_0_0_xzzzzz = cbuffer.data(si_geom_0010_off + 76 * ccomps * dcomps);

            auto g_0_0_z_0_0_yyyyyy = cbuffer.data(si_geom_0010_off + 77 * ccomps * dcomps);

            auto g_0_0_z_0_0_yyyyyz = cbuffer.data(si_geom_0010_off + 78 * ccomps * dcomps);

            auto g_0_0_z_0_0_yyyyzz = cbuffer.data(si_geom_0010_off + 79 * ccomps * dcomps);

            auto g_0_0_z_0_0_yyyzzz = cbuffer.data(si_geom_0010_off + 80 * ccomps * dcomps);

            auto g_0_0_z_0_0_yyzzzz = cbuffer.data(si_geom_0010_off + 81 * ccomps * dcomps);

            auto g_0_0_z_0_0_yzzzzz = cbuffer.data(si_geom_0010_off + 82 * ccomps * dcomps);

            auto g_0_0_z_0_0_zzzzzz = cbuffer.data(si_geom_0010_off + 83 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SISS

            const auto si_geom_1010_off = idx_geom_1010_sixx + i * dcomps + j;

            auto g_x_0_x_0_0_xxxxxx = cbuffer.data(si_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxxy = cbuffer.data(si_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxxz = cbuffer.data(si_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxyy = cbuffer.data(si_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxyz = cbuffer.data(si_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxzz = cbuffer.data(si_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxyyy = cbuffer.data(si_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxyyz = cbuffer.data(si_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxyzz = cbuffer.data(si_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxzzz = cbuffer.data(si_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyyyy = cbuffer.data(si_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyyyz = cbuffer.data(si_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyyzz = cbuffer.data(si_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyzzz = cbuffer.data(si_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxzzzz = cbuffer.data(si_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyyyy = cbuffer.data(si_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyyyz = cbuffer.data(si_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyyzz = cbuffer.data(si_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyzzz = cbuffer.data(si_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyzzzz = cbuffer.data(si_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_0_xzzzzz = cbuffer.data(si_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyyyy = cbuffer.data(si_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyyyz = cbuffer.data(si_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyyzz = cbuffer.data(si_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyzzz = cbuffer.data(si_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyzzzz = cbuffer.data(si_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_0_yzzzzz = cbuffer.data(si_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_0_zzzzzz = cbuffer.data(si_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxxx = cbuffer.data(si_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxxy = cbuffer.data(si_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxxz = cbuffer.data(si_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxyy = cbuffer.data(si_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxyz = cbuffer.data(si_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxzz = cbuffer.data(si_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxyyy = cbuffer.data(si_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxyyz = cbuffer.data(si_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxyzz = cbuffer.data(si_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxzzz = cbuffer.data(si_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyyyy = cbuffer.data(si_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyyyz = cbuffer.data(si_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyyzz = cbuffer.data(si_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyzzz = cbuffer.data(si_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxzzzz = cbuffer.data(si_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyyyy = cbuffer.data(si_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyyyz = cbuffer.data(si_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyyzz = cbuffer.data(si_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyzzz = cbuffer.data(si_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyzzzz = cbuffer.data(si_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_y_0_0_xzzzzz = cbuffer.data(si_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyyyy = cbuffer.data(si_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyyyz = cbuffer.data(si_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyyzz = cbuffer.data(si_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyzzz = cbuffer.data(si_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyzzzz = cbuffer.data(si_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_y_0_0_yzzzzz = cbuffer.data(si_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_y_0_0_zzzzzz = cbuffer.data(si_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxxx = cbuffer.data(si_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxxy = cbuffer.data(si_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxxz = cbuffer.data(si_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxyy = cbuffer.data(si_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxyz = cbuffer.data(si_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxzz = cbuffer.data(si_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxyyy = cbuffer.data(si_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxyyz = cbuffer.data(si_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxyzz = cbuffer.data(si_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxzzz = cbuffer.data(si_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyyyy = cbuffer.data(si_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyyyz = cbuffer.data(si_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyyzz = cbuffer.data(si_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyzzz = cbuffer.data(si_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxzzzz = cbuffer.data(si_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyyyy = cbuffer.data(si_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyyyz = cbuffer.data(si_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyyzz = cbuffer.data(si_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyzzz = cbuffer.data(si_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyzzzz = cbuffer.data(si_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_z_0_0_xzzzzz = cbuffer.data(si_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyyyy = cbuffer.data(si_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyyyz = cbuffer.data(si_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyyzz = cbuffer.data(si_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyzzz = cbuffer.data(si_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyzzzz = cbuffer.data(si_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_z_0_0_yzzzzz = cbuffer.data(si_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_z_0_0_zzzzzz = cbuffer.data(si_geom_1010_off + 83 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxxx = cbuffer.data(si_geom_1010_off + 84 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxxy = cbuffer.data(si_geom_1010_off + 85 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxxz = cbuffer.data(si_geom_1010_off + 86 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxyy = cbuffer.data(si_geom_1010_off + 87 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxyz = cbuffer.data(si_geom_1010_off + 88 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxzz = cbuffer.data(si_geom_1010_off + 89 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxyyy = cbuffer.data(si_geom_1010_off + 90 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxyyz = cbuffer.data(si_geom_1010_off + 91 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxyzz = cbuffer.data(si_geom_1010_off + 92 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxzzz = cbuffer.data(si_geom_1010_off + 93 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyyyy = cbuffer.data(si_geom_1010_off + 94 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyyyz = cbuffer.data(si_geom_1010_off + 95 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyyzz = cbuffer.data(si_geom_1010_off + 96 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyzzz = cbuffer.data(si_geom_1010_off + 97 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxzzzz = cbuffer.data(si_geom_1010_off + 98 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyyyy = cbuffer.data(si_geom_1010_off + 99 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyyyz = cbuffer.data(si_geom_1010_off + 100 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyyzz = cbuffer.data(si_geom_1010_off + 101 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyzzz = cbuffer.data(si_geom_1010_off + 102 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyzzzz = cbuffer.data(si_geom_1010_off + 103 * ccomps * dcomps);

            auto g_y_0_x_0_0_xzzzzz = cbuffer.data(si_geom_1010_off + 104 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyyyy = cbuffer.data(si_geom_1010_off + 105 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyyyz = cbuffer.data(si_geom_1010_off + 106 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyyzz = cbuffer.data(si_geom_1010_off + 107 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyzzz = cbuffer.data(si_geom_1010_off + 108 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyzzzz = cbuffer.data(si_geom_1010_off + 109 * ccomps * dcomps);

            auto g_y_0_x_0_0_yzzzzz = cbuffer.data(si_geom_1010_off + 110 * ccomps * dcomps);

            auto g_y_0_x_0_0_zzzzzz = cbuffer.data(si_geom_1010_off + 111 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxxx = cbuffer.data(si_geom_1010_off + 112 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxxy = cbuffer.data(si_geom_1010_off + 113 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxxz = cbuffer.data(si_geom_1010_off + 114 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxyy = cbuffer.data(si_geom_1010_off + 115 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxyz = cbuffer.data(si_geom_1010_off + 116 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxzz = cbuffer.data(si_geom_1010_off + 117 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxyyy = cbuffer.data(si_geom_1010_off + 118 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxyyz = cbuffer.data(si_geom_1010_off + 119 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxyzz = cbuffer.data(si_geom_1010_off + 120 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxzzz = cbuffer.data(si_geom_1010_off + 121 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyyyy = cbuffer.data(si_geom_1010_off + 122 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyyyz = cbuffer.data(si_geom_1010_off + 123 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyyzz = cbuffer.data(si_geom_1010_off + 124 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyzzz = cbuffer.data(si_geom_1010_off + 125 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxzzzz = cbuffer.data(si_geom_1010_off + 126 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyyyy = cbuffer.data(si_geom_1010_off + 127 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyyyz = cbuffer.data(si_geom_1010_off + 128 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyyzz = cbuffer.data(si_geom_1010_off + 129 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyzzz = cbuffer.data(si_geom_1010_off + 130 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyzzzz = cbuffer.data(si_geom_1010_off + 131 * ccomps * dcomps);

            auto g_y_0_y_0_0_xzzzzz = cbuffer.data(si_geom_1010_off + 132 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyyyy = cbuffer.data(si_geom_1010_off + 133 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyyyz = cbuffer.data(si_geom_1010_off + 134 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyyzz = cbuffer.data(si_geom_1010_off + 135 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyzzz = cbuffer.data(si_geom_1010_off + 136 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyzzzz = cbuffer.data(si_geom_1010_off + 137 * ccomps * dcomps);

            auto g_y_0_y_0_0_yzzzzz = cbuffer.data(si_geom_1010_off + 138 * ccomps * dcomps);

            auto g_y_0_y_0_0_zzzzzz = cbuffer.data(si_geom_1010_off + 139 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxxx = cbuffer.data(si_geom_1010_off + 140 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxxy = cbuffer.data(si_geom_1010_off + 141 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxxz = cbuffer.data(si_geom_1010_off + 142 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxyy = cbuffer.data(si_geom_1010_off + 143 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxyz = cbuffer.data(si_geom_1010_off + 144 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxzz = cbuffer.data(si_geom_1010_off + 145 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxyyy = cbuffer.data(si_geom_1010_off + 146 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxyyz = cbuffer.data(si_geom_1010_off + 147 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxyzz = cbuffer.data(si_geom_1010_off + 148 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxzzz = cbuffer.data(si_geom_1010_off + 149 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyyyy = cbuffer.data(si_geom_1010_off + 150 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyyyz = cbuffer.data(si_geom_1010_off + 151 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyyzz = cbuffer.data(si_geom_1010_off + 152 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyzzz = cbuffer.data(si_geom_1010_off + 153 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxzzzz = cbuffer.data(si_geom_1010_off + 154 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyyyy = cbuffer.data(si_geom_1010_off + 155 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyyyz = cbuffer.data(si_geom_1010_off + 156 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyyzz = cbuffer.data(si_geom_1010_off + 157 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyzzz = cbuffer.data(si_geom_1010_off + 158 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyzzzz = cbuffer.data(si_geom_1010_off + 159 * ccomps * dcomps);

            auto g_y_0_z_0_0_xzzzzz = cbuffer.data(si_geom_1010_off + 160 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyyyy = cbuffer.data(si_geom_1010_off + 161 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyyyz = cbuffer.data(si_geom_1010_off + 162 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyyzz = cbuffer.data(si_geom_1010_off + 163 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyzzz = cbuffer.data(si_geom_1010_off + 164 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyzzzz = cbuffer.data(si_geom_1010_off + 165 * ccomps * dcomps);

            auto g_y_0_z_0_0_yzzzzz = cbuffer.data(si_geom_1010_off + 166 * ccomps * dcomps);

            auto g_y_0_z_0_0_zzzzzz = cbuffer.data(si_geom_1010_off + 167 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxxx = cbuffer.data(si_geom_1010_off + 168 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxxy = cbuffer.data(si_geom_1010_off + 169 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxxz = cbuffer.data(si_geom_1010_off + 170 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxyy = cbuffer.data(si_geom_1010_off + 171 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxyz = cbuffer.data(si_geom_1010_off + 172 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxzz = cbuffer.data(si_geom_1010_off + 173 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxyyy = cbuffer.data(si_geom_1010_off + 174 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxyyz = cbuffer.data(si_geom_1010_off + 175 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxyzz = cbuffer.data(si_geom_1010_off + 176 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxzzz = cbuffer.data(si_geom_1010_off + 177 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyyyy = cbuffer.data(si_geom_1010_off + 178 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyyyz = cbuffer.data(si_geom_1010_off + 179 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyyzz = cbuffer.data(si_geom_1010_off + 180 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyzzz = cbuffer.data(si_geom_1010_off + 181 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxzzzz = cbuffer.data(si_geom_1010_off + 182 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyyyy = cbuffer.data(si_geom_1010_off + 183 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyyyz = cbuffer.data(si_geom_1010_off + 184 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyyzz = cbuffer.data(si_geom_1010_off + 185 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyzzz = cbuffer.data(si_geom_1010_off + 186 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyzzzz = cbuffer.data(si_geom_1010_off + 187 * ccomps * dcomps);

            auto g_z_0_x_0_0_xzzzzz = cbuffer.data(si_geom_1010_off + 188 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyyyy = cbuffer.data(si_geom_1010_off + 189 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyyyz = cbuffer.data(si_geom_1010_off + 190 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyyzz = cbuffer.data(si_geom_1010_off + 191 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyzzz = cbuffer.data(si_geom_1010_off + 192 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyzzzz = cbuffer.data(si_geom_1010_off + 193 * ccomps * dcomps);

            auto g_z_0_x_0_0_yzzzzz = cbuffer.data(si_geom_1010_off + 194 * ccomps * dcomps);

            auto g_z_0_x_0_0_zzzzzz = cbuffer.data(si_geom_1010_off + 195 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxxx = cbuffer.data(si_geom_1010_off + 196 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxxy = cbuffer.data(si_geom_1010_off + 197 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxxz = cbuffer.data(si_geom_1010_off + 198 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxyy = cbuffer.data(si_geom_1010_off + 199 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxyz = cbuffer.data(si_geom_1010_off + 200 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxzz = cbuffer.data(si_geom_1010_off + 201 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxyyy = cbuffer.data(si_geom_1010_off + 202 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxyyz = cbuffer.data(si_geom_1010_off + 203 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxyzz = cbuffer.data(si_geom_1010_off + 204 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxzzz = cbuffer.data(si_geom_1010_off + 205 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyyyy = cbuffer.data(si_geom_1010_off + 206 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyyyz = cbuffer.data(si_geom_1010_off + 207 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyyzz = cbuffer.data(si_geom_1010_off + 208 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyzzz = cbuffer.data(si_geom_1010_off + 209 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxzzzz = cbuffer.data(si_geom_1010_off + 210 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyyyy = cbuffer.data(si_geom_1010_off + 211 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyyyz = cbuffer.data(si_geom_1010_off + 212 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyyzz = cbuffer.data(si_geom_1010_off + 213 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyzzz = cbuffer.data(si_geom_1010_off + 214 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyzzzz = cbuffer.data(si_geom_1010_off + 215 * ccomps * dcomps);

            auto g_z_0_y_0_0_xzzzzz = cbuffer.data(si_geom_1010_off + 216 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyyyy = cbuffer.data(si_geom_1010_off + 217 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyyyz = cbuffer.data(si_geom_1010_off + 218 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyyzz = cbuffer.data(si_geom_1010_off + 219 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyzzz = cbuffer.data(si_geom_1010_off + 220 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyzzzz = cbuffer.data(si_geom_1010_off + 221 * ccomps * dcomps);

            auto g_z_0_y_0_0_yzzzzz = cbuffer.data(si_geom_1010_off + 222 * ccomps * dcomps);

            auto g_z_0_y_0_0_zzzzzz = cbuffer.data(si_geom_1010_off + 223 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxxx = cbuffer.data(si_geom_1010_off + 224 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxxy = cbuffer.data(si_geom_1010_off + 225 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxxz = cbuffer.data(si_geom_1010_off + 226 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxyy = cbuffer.data(si_geom_1010_off + 227 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxyz = cbuffer.data(si_geom_1010_off + 228 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxzz = cbuffer.data(si_geom_1010_off + 229 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxyyy = cbuffer.data(si_geom_1010_off + 230 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxyyz = cbuffer.data(si_geom_1010_off + 231 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxyzz = cbuffer.data(si_geom_1010_off + 232 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxzzz = cbuffer.data(si_geom_1010_off + 233 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyyyy = cbuffer.data(si_geom_1010_off + 234 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyyyz = cbuffer.data(si_geom_1010_off + 235 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyyzz = cbuffer.data(si_geom_1010_off + 236 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyzzz = cbuffer.data(si_geom_1010_off + 237 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxzzzz = cbuffer.data(si_geom_1010_off + 238 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyyyy = cbuffer.data(si_geom_1010_off + 239 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyyyz = cbuffer.data(si_geom_1010_off + 240 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyyzz = cbuffer.data(si_geom_1010_off + 241 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyzzz = cbuffer.data(si_geom_1010_off + 242 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyzzzz = cbuffer.data(si_geom_1010_off + 243 * ccomps * dcomps);

            auto g_z_0_z_0_0_xzzzzz = cbuffer.data(si_geom_1010_off + 244 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyyyy = cbuffer.data(si_geom_1010_off + 245 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyyyz = cbuffer.data(si_geom_1010_off + 246 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyyzz = cbuffer.data(si_geom_1010_off + 247 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyzzz = cbuffer.data(si_geom_1010_off + 248 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyzzzz = cbuffer.data(si_geom_1010_off + 249 * ccomps * dcomps);

            auto g_z_0_z_0_0_yzzzzz = cbuffer.data(si_geom_1010_off + 250 * ccomps * dcomps);

            auto g_z_0_z_0_0_zzzzzz = cbuffer.data(si_geom_1010_off + 251 * ccomps * dcomps);

            /// Set up components of auxilary buffer : SKSS

            const auto sk_geom_1010_off = idx_geom_1010_skxx + i * dcomps + j;

            auto g_x_0_x_0_0_xxxxxxx = cbuffer.data(sk_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxxxy = cbuffer.data(sk_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxxxz = cbuffer.data(sk_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxxyy = cbuffer.data(sk_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxxyz = cbuffer.data(sk_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxxzz = cbuffer.data(sk_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxyyy = cbuffer.data(sk_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxyyz = cbuffer.data(sk_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxyzz = cbuffer.data(sk_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxxzzz = cbuffer.data(sk_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxyyyy = cbuffer.data(sk_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxyyyz = cbuffer.data(sk_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxyyzz = cbuffer.data(sk_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxyzzz = cbuffer.data(sk_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxxzzzz = cbuffer.data(sk_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyyyyy = cbuffer.data(sk_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyyyyz = cbuffer.data(sk_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyyyzz = cbuffer.data(sk_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyyzzz = cbuffer.data(sk_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxyzzzz = cbuffer.data(sk_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_0_xxzzzzz = cbuffer.data(sk_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyyyyy = cbuffer.data(sk_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyyyyz = cbuffer.data(sk_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyyyzz = cbuffer.data(sk_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyyzzz = cbuffer.data(sk_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyyzzzz = cbuffer.data(sk_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_0_xyzzzzz = cbuffer.data(sk_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_0_xzzzzzz = cbuffer.data(sk_geom_1010_off + 27 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyyyyy = cbuffer.data(sk_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyyyyz = cbuffer.data(sk_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyyyzz = cbuffer.data(sk_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyyzzz = cbuffer.data(sk_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyyzzzz = cbuffer.data(sk_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_0_yyzzzzz = cbuffer.data(sk_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_0_yzzzzzz = cbuffer.data(sk_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_0_zzzzzzz = cbuffer.data(sk_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxxxx = cbuffer.data(sk_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxxxy = cbuffer.data(sk_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxxxz = cbuffer.data(sk_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxxyy = cbuffer.data(sk_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxxyz = cbuffer.data(sk_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxxzz = cbuffer.data(sk_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxyyy = cbuffer.data(sk_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxyyz = cbuffer.data(sk_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxyzz = cbuffer.data(sk_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxxzzz = cbuffer.data(sk_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxyyyy = cbuffer.data(sk_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxyyyz = cbuffer.data(sk_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxyyzz = cbuffer.data(sk_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxyzzz = cbuffer.data(sk_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxxzzzz = cbuffer.data(sk_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyyyyy = cbuffer.data(sk_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyyyyz = cbuffer.data(sk_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyyyzz = cbuffer.data(sk_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyyzzz = cbuffer.data(sk_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxyzzzz = cbuffer.data(sk_geom_1010_off + 55 * ccomps * dcomps);

            auto g_x_0_y_0_0_xxzzzzz = cbuffer.data(sk_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyyyyy = cbuffer.data(sk_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyyyyz = cbuffer.data(sk_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyyyzz = cbuffer.data(sk_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyyzzz = cbuffer.data(sk_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyyzzzz = cbuffer.data(sk_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_y_0_0_xyzzzzz = cbuffer.data(sk_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_y_0_0_xzzzzzz = cbuffer.data(sk_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyyyyy = cbuffer.data(sk_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyyyyz = cbuffer.data(sk_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyyyzz = cbuffer.data(sk_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyyzzz = cbuffer.data(sk_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyyzzzz = cbuffer.data(sk_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_y_0_0_yyzzzzz = cbuffer.data(sk_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_y_0_0_yzzzzzz = cbuffer.data(sk_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_y_0_0_zzzzzzz = cbuffer.data(sk_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxxxx = cbuffer.data(sk_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxxxy = cbuffer.data(sk_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxxxz = cbuffer.data(sk_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxxyy = cbuffer.data(sk_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxxyz = cbuffer.data(sk_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxxzz = cbuffer.data(sk_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxyyy = cbuffer.data(sk_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxyyz = cbuffer.data(sk_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxyzz = cbuffer.data(sk_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxxzzz = cbuffer.data(sk_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxyyyy = cbuffer.data(sk_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxyyyz = cbuffer.data(sk_geom_1010_off + 83 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxyyzz = cbuffer.data(sk_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxyzzz = cbuffer.data(sk_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxxzzzz = cbuffer.data(sk_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyyyyy = cbuffer.data(sk_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyyyyz = cbuffer.data(sk_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyyyzz = cbuffer.data(sk_geom_1010_off + 89 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyyzzz = cbuffer.data(sk_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxyzzzz = cbuffer.data(sk_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_z_0_0_xxzzzzz = cbuffer.data(sk_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyyyyy = cbuffer.data(sk_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyyyyz = cbuffer.data(sk_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyyyzz = cbuffer.data(sk_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyyzzz = cbuffer.data(sk_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyyzzzz = cbuffer.data(sk_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_z_0_0_xyzzzzz = cbuffer.data(sk_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_z_0_0_xzzzzzz = cbuffer.data(sk_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyyyyy = cbuffer.data(sk_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyyyyz = cbuffer.data(sk_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyyyzz = cbuffer.data(sk_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyyzzz = cbuffer.data(sk_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyyzzzz = cbuffer.data(sk_geom_1010_off + 104 * ccomps * dcomps);

            auto g_x_0_z_0_0_yyzzzzz = cbuffer.data(sk_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_z_0_0_yzzzzzz = cbuffer.data(sk_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_z_0_0_zzzzzzz = cbuffer.data(sk_geom_1010_off + 107 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxxxx = cbuffer.data(sk_geom_1010_off + 108 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxxxy = cbuffer.data(sk_geom_1010_off + 109 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxxxz = cbuffer.data(sk_geom_1010_off + 110 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxxyy = cbuffer.data(sk_geom_1010_off + 111 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxxyz = cbuffer.data(sk_geom_1010_off + 112 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxxzz = cbuffer.data(sk_geom_1010_off + 113 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxyyy = cbuffer.data(sk_geom_1010_off + 114 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxyyz = cbuffer.data(sk_geom_1010_off + 115 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxyzz = cbuffer.data(sk_geom_1010_off + 116 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxxzzz = cbuffer.data(sk_geom_1010_off + 117 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxyyyy = cbuffer.data(sk_geom_1010_off + 118 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxyyyz = cbuffer.data(sk_geom_1010_off + 119 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxyyzz = cbuffer.data(sk_geom_1010_off + 120 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxyzzz = cbuffer.data(sk_geom_1010_off + 121 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxxzzzz = cbuffer.data(sk_geom_1010_off + 122 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyyyyy = cbuffer.data(sk_geom_1010_off + 123 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyyyyz = cbuffer.data(sk_geom_1010_off + 124 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyyyzz = cbuffer.data(sk_geom_1010_off + 125 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyyzzz = cbuffer.data(sk_geom_1010_off + 126 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxyzzzz = cbuffer.data(sk_geom_1010_off + 127 * ccomps * dcomps);

            auto g_y_0_x_0_0_xxzzzzz = cbuffer.data(sk_geom_1010_off + 128 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyyyyy = cbuffer.data(sk_geom_1010_off + 129 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyyyyz = cbuffer.data(sk_geom_1010_off + 130 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyyyzz = cbuffer.data(sk_geom_1010_off + 131 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyyzzz = cbuffer.data(sk_geom_1010_off + 132 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyyzzzz = cbuffer.data(sk_geom_1010_off + 133 * ccomps * dcomps);

            auto g_y_0_x_0_0_xyzzzzz = cbuffer.data(sk_geom_1010_off + 134 * ccomps * dcomps);

            auto g_y_0_x_0_0_xzzzzzz = cbuffer.data(sk_geom_1010_off + 135 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyyyyy = cbuffer.data(sk_geom_1010_off + 136 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyyyyz = cbuffer.data(sk_geom_1010_off + 137 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyyyzz = cbuffer.data(sk_geom_1010_off + 138 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyyzzz = cbuffer.data(sk_geom_1010_off + 139 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyyzzzz = cbuffer.data(sk_geom_1010_off + 140 * ccomps * dcomps);

            auto g_y_0_x_0_0_yyzzzzz = cbuffer.data(sk_geom_1010_off + 141 * ccomps * dcomps);

            auto g_y_0_x_0_0_yzzzzzz = cbuffer.data(sk_geom_1010_off + 142 * ccomps * dcomps);

            auto g_y_0_x_0_0_zzzzzzz = cbuffer.data(sk_geom_1010_off + 143 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxxxx = cbuffer.data(sk_geom_1010_off + 144 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxxxy = cbuffer.data(sk_geom_1010_off + 145 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxxxz = cbuffer.data(sk_geom_1010_off + 146 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxxyy = cbuffer.data(sk_geom_1010_off + 147 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxxyz = cbuffer.data(sk_geom_1010_off + 148 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxxzz = cbuffer.data(sk_geom_1010_off + 149 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxyyy = cbuffer.data(sk_geom_1010_off + 150 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxyyz = cbuffer.data(sk_geom_1010_off + 151 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxyzz = cbuffer.data(sk_geom_1010_off + 152 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxxzzz = cbuffer.data(sk_geom_1010_off + 153 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxyyyy = cbuffer.data(sk_geom_1010_off + 154 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxyyyz = cbuffer.data(sk_geom_1010_off + 155 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxyyzz = cbuffer.data(sk_geom_1010_off + 156 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxyzzz = cbuffer.data(sk_geom_1010_off + 157 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxxzzzz = cbuffer.data(sk_geom_1010_off + 158 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyyyyy = cbuffer.data(sk_geom_1010_off + 159 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyyyyz = cbuffer.data(sk_geom_1010_off + 160 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyyyzz = cbuffer.data(sk_geom_1010_off + 161 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyyzzz = cbuffer.data(sk_geom_1010_off + 162 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxyzzzz = cbuffer.data(sk_geom_1010_off + 163 * ccomps * dcomps);

            auto g_y_0_y_0_0_xxzzzzz = cbuffer.data(sk_geom_1010_off + 164 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyyyyy = cbuffer.data(sk_geom_1010_off + 165 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyyyyz = cbuffer.data(sk_geom_1010_off + 166 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyyyzz = cbuffer.data(sk_geom_1010_off + 167 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyyzzz = cbuffer.data(sk_geom_1010_off + 168 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyyzzzz = cbuffer.data(sk_geom_1010_off + 169 * ccomps * dcomps);

            auto g_y_0_y_0_0_xyzzzzz = cbuffer.data(sk_geom_1010_off + 170 * ccomps * dcomps);

            auto g_y_0_y_0_0_xzzzzzz = cbuffer.data(sk_geom_1010_off + 171 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyyyyy = cbuffer.data(sk_geom_1010_off + 172 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyyyyz = cbuffer.data(sk_geom_1010_off + 173 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyyyzz = cbuffer.data(sk_geom_1010_off + 174 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyyzzz = cbuffer.data(sk_geom_1010_off + 175 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyyzzzz = cbuffer.data(sk_geom_1010_off + 176 * ccomps * dcomps);

            auto g_y_0_y_0_0_yyzzzzz = cbuffer.data(sk_geom_1010_off + 177 * ccomps * dcomps);

            auto g_y_0_y_0_0_yzzzzzz = cbuffer.data(sk_geom_1010_off + 178 * ccomps * dcomps);

            auto g_y_0_y_0_0_zzzzzzz = cbuffer.data(sk_geom_1010_off + 179 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxxxx = cbuffer.data(sk_geom_1010_off + 180 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxxxy = cbuffer.data(sk_geom_1010_off + 181 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxxxz = cbuffer.data(sk_geom_1010_off + 182 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxxyy = cbuffer.data(sk_geom_1010_off + 183 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxxyz = cbuffer.data(sk_geom_1010_off + 184 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxxzz = cbuffer.data(sk_geom_1010_off + 185 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxyyy = cbuffer.data(sk_geom_1010_off + 186 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxyyz = cbuffer.data(sk_geom_1010_off + 187 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxyzz = cbuffer.data(sk_geom_1010_off + 188 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxxzzz = cbuffer.data(sk_geom_1010_off + 189 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxyyyy = cbuffer.data(sk_geom_1010_off + 190 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxyyyz = cbuffer.data(sk_geom_1010_off + 191 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxyyzz = cbuffer.data(sk_geom_1010_off + 192 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxyzzz = cbuffer.data(sk_geom_1010_off + 193 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxxzzzz = cbuffer.data(sk_geom_1010_off + 194 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyyyyy = cbuffer.data(sk_geom_1010_off + 195 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyyyyz = cbuffer.data(sk_geom_1010_off + 196 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyyyzz = cbuffer.data(sk_geom_1010_off + 197 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyyzzz = cbuffer.data(sk_geom_1010_off + 198 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxyzzzz = cbuffer.data(sk_geom_1010_off + 199 * ccomps * dcomps);

            auto g_y_0_z_0_0_xxzzzzz = cbuffer.data(sk_geom_1010_off + 200 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyyyyy = cbuffer.data(sk_geom_1010_off + 201 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyyyyz = cbuffer.data(sk_geom_1010_off + 202 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyyyzz = cbuffer.data(sk_geom_1010_off + 203 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyyzzz = cbuffer.data(sk_geom_1010_off + 204 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyyzzzz = cbuffer.data(sk_geom_1010_off + 205 * ccomps * dcomps);

            auto g_y_0_z_0_0_xyzzzzz = cbuffer.data(sk_geom_1010_off + 206 * ccomps * dcomps);

            auto g_y_0_z_0_0_xzzzzzz = cbuffer.data(sk_geom_1010_off + 207 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyyyyy = cbuffer.data(sk_geom_1010_off + 208 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyyyyz = cbuffer.data(sk_geom_1010_off + 209 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyyyzz = cbuffer.data(sk_geom_1010_off + 210 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyyzzz = cbuffer.data(sk_geom_1010_off + 211 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyyzzzz = cbuffer.data(sk_geom_1010_off + 212 * ccomps * dcomps);

            auto g_y_0_z_0_0_yyzzzzz = cbuffer.data(sk_geom_1010_off + 213 * ccomps * dcomps);

            auto g_y_0_z_0_0_yzzzzzz = cbuffer.data(sk_geom_1010_off + 214 * ccomps * dcomps);

            auto g_y_0_z_0_0_zzzzzzz = cbuffer.data(sk_geom_1010_off + 215 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxxxx = cbuffer.data(sk_geom_1010_off + 216 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxxxy = cbuffer.data(sk_geom_1010_off + 217 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxxxz = cbuffer.data(sk_geom_1010_off + 218 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxxyy = cbuffer.data(sk_geom_1010_off + 219 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxxyz = cbuffer.data(sk_geom_1010_off + 220 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxxzz = cbuffer.data(sk_geom_1010_off + 221 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxyyy = cbuffer.data(sk_geom_1010_off + 222 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxyyz = cbuffer.data(sk_geom_1010_off + 223 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxyzz = cbuffer.data(sk_geom_1010_off + 224 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxxzzz = cbuffer.data(sk_geom_1010_off + 225 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxyyyy = cbuffer.data(sk_geom_1010_off + 226 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxyyyz = cbuffer.data(sk_geom_1010_off + 227 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxyyzz = cbuffer.data(sk_geom_1010_off + 228 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxyzzz = cbuffer.data(sk_geom_1010_off + 229 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxxzzzz = cbuffer.data(sk_geom_1010_off + 230 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyyyyy = cbuffer.data(sk_geom_1010_off + 231 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyyyyz = cbuffer.data(sk_geom_1010_off + 232 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyyyzz = cbuffer.data(sk_geom_1010_off + 233 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyyzzz = cbuffer.data(sk_geom_1010_off + 234 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxyzzzz = cbuffer.data(sk_geom_1010_off + 235 * ccomps * dcomps);

            auto g_z_0_x_0_0_xxzzzzz = cbuffer.data(sk_geom_1010_off + 236 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyyyyy = cbuffer.data(sk_geom_1010_off + 237 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyyyyz = cbuffer.data(sk_geom_1010_off + 238 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyyyzz = cbuffer.data(sk_geom_1010_off + 239 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyyzzz = cbuffer.data(sk_geom_1010_off + 240 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyyzzzz = cbuffer.data(sk_geom_1010_off + 241 * ccomps * dcomps);

            auto g_z_0_x_0_0_xyzzzzz = cbuffer.data(sk_geom_1010_off + 242 * ccomps * dcomps);

            auto g_z_0_x_0_0_xzzzzzz = cbuffer.data(sk_geom_1010_off + 243 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyyyyy = cbuffer.data(sk_geom_1010_off + 244 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyyyyz = cbuffer.data(sk_geom_1010_off + 245 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyyyzz = cbuffer.data(sk_geom_1010_off + 246 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyyzzz = cbuffer.data(sk_geom_1010_off + 247 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyyzzzz = cbuffer.data(sk_geom_1010_off + 248 * ccomps * dcomps);

            auto g_z_0_x_0_0_yyzzzzz = cbuffer.data(sk_geom_1010_off + 249 * ccomps * dcomps);

            auto g_z_0_x_0_0_yzzzzzz = cbuffer.data(sk_geom_1010_off + 250 * ccomps * dcomps);

            auto g_z_0_x_0_0_zzzzzzz = cbuffer.data(sk_geom_1010_off + 251 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxxxx = cbuffer.data(sk_geom_1010_off + 252 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxxxy = cbuffer.data(sk_geom_1010_off + 253 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxxxz = cbuffer.data(sk_geom_1010_off + 254 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxxyy = cbuffer.data(sk_geom_1010_off + 255 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxxyz = cbuffer.data(sk_geom_1010_off + 256 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxxzz = cbuffer.data(sk_geom_1010_off + 257 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxyyy = cbuffer.data(sk_geom_1010_off + 258 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxyyz = cbuffer.data(sk_geom_1010_off + 259 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxyzz = cbuffer.data(sk_geom_1010_off + 260 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxxzzz = cbuffer.data(sk_geom_1010_off + 261 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxyyyy = cbuffer.data(sk_geom_1010_off + 262 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxyyyz = cbuffer.data(sk_geom_1010_off + 263 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxyyzz = cbuffer.data(sk_geom_1010_off + 264 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxyzzz = cbuffer.data(sk_geom_1010_off + 265 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxxzzzz = cbuffer.data(sk_geom_1010_off + 266 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyyyyy = cbuffer.data(sk_geom_1010_off + 267 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyyyyz = cbuffer.data(sk_geom_1010_off + 268 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyyyzz = cbuffer.data(sk_geom_1010_off + 269 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyyzzz = cbuffer.data(sk_geom_1010_off + 270 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxyzzzz = cbuffer.data(sk_geom_1010_off + 271 * ccomps * dcomps);

            auto g_z_0_y_0_0_xxzzzzz = cbuffer.data(sk_geom_1010_off + 272 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyyyyy = cbuffer.data(sk_geom_1010_off + 273 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyyyyz = cbuffer.data(sk_geom_1010_off + 274 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyyyzz = cbuffer.data(sk_geom_1010_off + 275 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyyzzz = cbuffer.data(sk_geom_1010_off + 276 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyyzzzz = cbuffer.data(sk_geom_1010_off + 277 * ccomps * dcomps);

            auto g_z_0_y_0_0_xyzzzzz = cbuffer.data(sk_geom_1010_off + 278 * ccomps * dcomps);

            auto g_z_0_y_0_0_xzzzzzz = cbuffer.data(sk_geom_1010_off + 279 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyyyyy = cbuffer.data(sk_geom_1010_off + 280 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyyyyz = cbuffer.data(sk_geom_1010_off + 281 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyyyzz = cbuffer.data(sk_geom_1010_off + 282 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyyzzz = cbuffer.data(sk_geom_1010_off + 283 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyyzzzz = cbuffer.data(sk_geom_1010_off + 284 * ccomps * dcomps);

            auto g_z_0_y_0_0_yyzzzzz = cbuffer.data(sk_geom_1010_off + 285 * ccomps * dcomps);

            auto g_z_0_y_0_0_yzzzzzz = cbuffer.data(sk_geom_1010_off + 286 * ccomps * dcomps);

            auto g_z_0_y_0_0_zzzzzzz = cbuffer.data(sk_geom_1010_off + 287 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxxxx = cbuffer.data(sk_geom_1010_off + 288 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxxxy = cbuffer.data(sk_geom_1010_off + 289 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxxxz = cbuffer.data(sk_geom_1010_off + 290 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxxyy = cbuffer.data(sk_geom_1010_off + 291 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxxyz = cbuffer.data(sk_geom_1010_off + 292 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxxzz = cbuffer.data(sk_geom_1010_off + 293 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxyyy = cbuffer.data(sk_geom_1010_off + 294 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxyyz = cbuffer.data(sk_geom_1010_off + 295 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxyzz = cbuffer.data(sk_geom_1010_off + 296 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxxzzz = cbuffer.data(sk_geom_1010_off + 297 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxyyyy = cbuffer.data(sk_geom_1010_off + 298 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxyyyz = cbuffer.data(sk_geom_1010_off + 299 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxyyzz = cbuffer.data(sk_geom_1010_off + 300 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxyzzz = cbuffer.data(sk_geom_1010_off + 301 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxxzzzz = cbuffer.data(sk_geom_1010_off + 302 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyyyyy = cbuffer.data(sk_geom_1010_off + 303 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyyyyz = cbuffer.data(sk_geom_1010_off + 304 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyyyzz = cbuffer.data(sk_geom_1010_off + 305 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyyzzz = cbuffer.data(sk_geom_1010_off + 306 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxyzzzz = cbuffer.data(sk_geom_1010_off + 307 * ccomps * dcomps);

            auto g_z_0_z_0_0_xxzzzzz = cbuffer.data(sk_geom_1010_off + 308 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyyyyy = cbuffer.data(sk_geom_1010_off + 309 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyyyyz = cbuffer.data(sk_geom_1010_off + 310 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyyyzz = cbuffer.data(sk_geom_1010_off + 311 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyyzzz = cbuffer.data(sk_geom_1010_off + 312 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyyzzzz = cbuffer.data(sk_geom_1010_off + 313 * ccomps * dcomps);

            auto g_z_0_z_0_0_xyzzzzz = cbuffer.data(sk_geom_1010_off + 314 * ccomps * dcomps);

            auto g_z_0_z_0_0_xzzzzzz = cbuffer.data(sk_geom_1010_off + 315 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyyyyy = cbuffer.data(sk_geom_1010_off + 316 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyyyyz = cbuffer.data(sk_geom_1010_off + 317 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyyyzz = cbuffer.data(sk_geom_1010_off + 318 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyyzzz = cbuffer.data(sk_geom_1010_off + 319 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyyzzzz = cbuffer.data(sk_geom_1010_off + 320 * ccomps * dcomps);

            auto g_z_0_z_0_0_yyzzzzz = cbuffer.data(sk_geom_1010_off + 321 * ccomps * dcomps);

            auto g_z_0_z_0_0_yzzzzzz = cbuffer.data(sk_geom_1010_off + 322 * ccomps * dcomps);

            auto g_z_0_z_0_0_zzzzzzz = cbuffer.data(sk_geom_1010_off + 323 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_pixx

            const auto pi_geom_1010_off = idx_geom_1010_pixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_x_xxxxxx = cbuffer.data(pi_geom_1010_off + 0 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxxxy = cbuffer.data(pi_geom_1010_off + 1 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxxxz = cbuffer.data(pi_geom_1010_off + 2 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxxyy = cbuffer.data(pi_geom_1010_off + 3 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxxyz = cbuffer.data(pi_geom_1010_off + 4 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxxzz = cbuffer.data(pi_geom_1010_off + 5 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxyyy = cbuffer.data(pi_geom_1010_off + 6 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxyyz = cbuffer.data(pi_geom_1010_off + 7 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxyzz = cbuffer.data(pi_geom_1010_off + 8 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxxzzz = cbuffer.data(pi_geom_1010_off + 9 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxyyyy = cbuffer.data(pi_geom_1010_off + 10 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxyyyz = cbuffer.data(pi_geom_1010_off + 11 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxyyzz = cbuffer.data(pi_geom_1010_off + 12 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxyzzz = cbuffer.data(pi_geom_1010_off + 13 * ccomps * dcomps);

            auto g_x_0_x_0_x_xxzzzz = cbuffer.data(pi_geom_1010_off + 14 * ccomps * dcomps);

            auto g_x_0_x_0_x_xyyyyy = cbuffer.data(pi_geom_1010_off + 15 * ccomps * dcomps);

            auto g_x_0_x_0_x_xyyyyz = cbuffer.data(pi_geom_1010_off + 16 * ccomps * dcomps);

            auto g_x_0_x_0_x_xyyyzz = cbuffer.data(pi_geom_1010_off + 17 * ccomps * dcomps);

            auto g_x_0_x_0_x_xyyzzz = cbuffer.data(pi_geom_1010_off + 18 * ccomps * dcomps);

            auto g_x_0_x_0_x_xyzzzz = cbuffer.data(pi_geom_1010_off + 19 * ccomps * dcomps);

            auto g_x_0_x_0_x_xzzzzz = cbuffer.data(pi_geom_1010_off + 20 * ccomps * dcomps);

            auto g_x_0_x_0_x_yyyyyy = cbuffer.data(pi_geom_1010_off + 21 * ccomps * dcomps);

            auto g_x_0_x_0_x_yyyyyz = cbuffer.data(pi_geom_1010_off + 22 * ccomps * dcomps);

            auto g_x_0_x_0_x_yyyyzz = cbuffer.data(pi_geom_1010_off + 23 * ccomps * dcomps);

            auto g_x_0_x_0_x_yyyzzz = cbuffer.data(pi_geom_1010_off + 24 * ccomps * dcomps);

            auto g_x_0_x_0_x_yyzzzz = cbuffer.data(pi_geom_1010_off + 25 * ccomps * dcomps);

            auto g_x_0_x_0_x_yzzzzz = cbuffer.data(pi_geom_1010_off + 26 * ccomps * dcomps);

            auto g_x_0_x_0_x_zzzzzz = cbuffer.data(pi_geom_1010_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_xxxxxx, g_0_0_x_0_0_xxxxxy, g_0_0_x_0_0_xxxxxz, g_0_0_x_0_0_xxxxyy, g_0_0_x_0_0_xxxxyz, g_0_0_x_0_0_xxxxzz, g_0_0_x_0_0_xxxyyy, g_0_0_x_0_0_xxxyyz, g_0_0_x_0_0_xxxyzz, g_0_0_x_0_0_xxxzzz, g_0_0_x_0_0_xxyyyy, g_0_0_x_0_0_xxyyyz, g_0_0_x_0_0_xxyyzz, g_0_0_x_0_0_xxyzzz, g_0_0_x_0_0_xxzzzz, g_0_0_x_0_0_xyyyyy, g_0_0_x_0_0_xyyyyz, g_0_0_x_0_0_xyyyzz, g_0_0_x_0_0_xyyzzz, g_0_0_x_0_0_xyzzzz, g_0_0_x_0_0_xzzzzz, g_0_0_x_0_0_yyyyyy, g_0_0_x_0_0_yyyyyz, g_0_0_x_0_0_yyyyzz, g_0_0_x_0_0_yyyzzz, g_0_0_x_0_0_yyzzzz, g_0_0_x_0_0_yzzzzz, g_0_0_x_0_0_zzzzzz, g_x_0_x_0_0_xxxxxx, g_x_0_x_0_0_xxxxxxx, g_x_0_x_0_0_xxxxxxy, g_x_0_x_0_0_xxxxxxz, g_x_0_x_0_0_xxxxxy, g_x_0_x_0_0_xxxxxyy, g_x_0_x_0_0_xxxxxyz, g_x_0_x_0_0_xxxxxz, g_x_0_x_0_0_xxxxxzz, g_x_0_x_0_0_xxxxyy, g_x_0_x_0_0_xxxxyyy, g_x_0_x_0_0_xxxxyyz, g_x_0_x_0_0_xxxxyz, g_x_0_x_0_0_xxxxyzz, g_x_0_x_0_0_xxxxzz, g_x_0_x_0_0_xxxxzzz, g_x_0_x_0_0_xxxyyy, g_x_0_x_0_0_xxxyyyy, g_x_0_x_0_0_xxxyyyz, g_x_0_x_0_0_xxxyyz, g_x_0_x_0_0_xxxyyzz, g_x_0_x_0_0_xxxyzz, g_x_0_x_0_0_xxxyzzz, g_x_0_x_0_0_xxxzzz, g_x_0_x_0_0_xxxzzzz, g_x_0_x_0_0_xxyyyy, g_x_0_x_0_0_xxyyyyy, g_x_0_x_0_0_xxyyyyz, g_x_0_x_0_0_xxyyyz, g_x_0_x_0_0_xxyyyzz, g_x_0_x_0_0_xxyyzz, g_x_0_x_0_0_xxyyzzz, g_x_0_x_0_0_xxyzzz, g_x_0_x_0_0_xxyzzzz, g_x_0_x_0_0_xxzzzz, g_x_0_x_0_0_xxzzzzz, g_x_0_x_0_0_xyyyyy, g_x_0_x_0_0_xyyyyyy, g_x_0_x_0_0_xyyyyyz, g_x_0_x_0_0_xyyyyz, g_x_0_x_0_0_xyyyyzz, g_x_0_x_0_0_xyyyzz, g_x_0_x_0_0_xyyyzzz, g_x_0_x_0_0_xyyzzz, g_x_0_x_0_0_xyyzzzz, g_x_0_x_0_0_xyzzzz, g_x_0_x_0_0_xyzzzzz, g_x_0_x_0_0_xzzzzz, g_x_0_x_0_0_xzzzzzz, g_x_0_x_0_0_yyyyyy, g_x_0_x_0_0_yyyyyz, g_x_0_x_0_0_yyyyzz, g_x_0_x_0_0_yyyzzz, g_x_0_x_0_0_yyzzzz, g_x_0_x_0_0_yzzzzz, g_x_0_x_0_0_zzzzzz, g_x_0_x_0_x_xxxxxx, g_x_0_x_0_x_xxxxxy, g_x_0_x_0_x_xxxxxz, g_x_0_x_0_x_xxxxyy, g_x_0_x_0_x_xxxxyz, g_x_0_x_0_x_xxxxzz, g_x_0_x_0_x_xxxyyy, g_x_0_x_0_x_xxxyyz, g_x_0_x_0_x_xxxyzz, g_x_0_x_0_x_xxxzzz, g_x_0_x_0_x_xxyyyy, g_x_0_x_0_x_xxyyyz, g_x_0_x_0_x_xxyyzz, g_x_0_x_0_x_xxyzzz, g_x_0_x_0_x_xxzzzz, g_x_0_x_0_x_xyyyyy, g_x_0_x_0_x_xyyyyz, g_x_0_x_0_x_xyyyzz, g_x_0_x_0_x_xyyzzz, g_x_0_x_0_x_xyzzzz, g_x_0_x_0_x_xzzzzz, g_x_0_x_0_x_yyyyyy, g_x_0_x_0_x_yyyyyz, g_x_0_x_0_x_yyyyzz, g_x_0_x_0_x_yyyzzz, g_x_0_x_0_x_yyzzzz, g_x_0_x_0_x_yzzzzz, g_x_0_x_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_x_xxxxxx[k] = -g_0_0_x_0_0_xxxxxx[k] - g_x_0_x_0_0_xxxxxx[k] * ab_x + g_x_0_x_0_0_xxxxxxx[k];

                g_x_0_x_0_x_xxxxxy[k] = -g_0_0_x_0_0_xxxxxy[k] - g_x_0_x_0_0_xxxxxy[k] * ab_x + g_x_0_x_0_0_xxxxxxy[k];

                g_x_0_x_0_x_xxxxxz[k] = -g_0_0_x_0_0_xxxxxz[k] - g_x_0_x_0_0_xxxxxz[k] * ab_x + g_x_0_x_0_0_xxxxxxz[k];

                g_x_0_x_0_x_xxxxyy[k] = -g_0_0_x_0_0_xxxxyy[k] - g_x_0_x_0_0_xxxxyy[k] * ab_x + g_x_0_x_0_0_xxxxxyy[k];

                g_x_0_x_0_x_xxxxyz[k] = -g_0_0_x_0_0_xxxxyz[k] - g_x_0_x_0_0_xxxxyz[k] * ab_x + g_x_0_x_0_0_xxxxxyz[k];

                g_x_0_x_0_x_xxxxzz[k] = -g_0_0_x_0_0_xxxxzz[k] - g_x_0_x_0_0_xxxxzz[k] * ab_x + g_x_0_x_0_0_xxxxxzz[k];

                g_x_0_x_0_x_xxxyyy[k] = -g_0_0_x_0_0_xxxyyy[k] - g_x_0_x_0_0_xxxyyy[k] * ab_x + g_x_0_x_0_0_xxxxyyy[k];

                g_x_0_x_0_x_xxxyyz[k] = -g_0_0_x_0_0_xxxyyz[k] - g_x_0_x_0_0_xxxyyz[k] * ab_x + g_x_0_x_0_0_xxxxyyz[k];

                g_x_0_x_0_x_xxxyzz[k] = -g_0_0_x_0_0_xxxyzz[k] - g_x_0_x_0_0_xxxyzz[k] * ab_x + g_x_0_x_0_0_xxxxyzz[k];

                g_x_0_x_0_x_xxxzzz[k] = -g_0_0_x_0_0_xxxzzz[k] - g_x_0_x_0_0_xxxzzz[k] * ab_x + g_x_0_x_0_0_xxxxzzz[k];

                g_x_0_x_0_x_xxyyyy[k] = -g_0_0_x_0_0_xxyyyy[k] - g_x_0_x_0_0_xxyyyy[k] * ab_x + g_x_0_x_0_0_xxxyyyy[k];

                g_x_0_x_0_x_xxyyyz[k] = -g_0_0_x_0_0_xxyyyz[k] - g_x_0_x_0_0_xxyyyz[k] * ab_x + g_x_0_x_0_0_xxxyyyz[k];

                g_x_0_x_0_x_xxyyzz[k] = -g_0_0_x_0_0_xxyyzz[k] - g_x_0_x_0_0_xxyyzz[k] * ab_x + g_x_0_x_0_0_xxxyyzz[k];

                g_x_0_x_0_x_xxyzzz[k] = -g_0_0_x_0_0_xxyzzz[k] - g_x_0_x_0_0_xxyzzz[k] * ab_x + g_x_0_x_0_0_xxxyzzz[k];

                g_x_0_x_0_x_xxzzzz[k] = -g_0_0_x_0_0_xxzzzz[k] - g_x_0_x_0_0_xxzzzz[k] * ab_x + g_x_0_x_0_0_xxxzzzz[k];

                g_x_0_x_0_x_xyyyyy[k] = -g_0_0_x_0_0_xyyyyy[k] - g_x_0_x_0_0_xyyyyy[k] * ab_x + g_x_0_x_0_0_xxyyyyy[k];

                g_x_0_x_0_x_xyyyyz[k] = -g_0_0_x_0_0_xyyyyz[k] - g_x_0_x_0_0_xyyyyz[k] * ab_x + g_x_0_x_0_0_xxyyyyz[k];

                g_x_0_x_0_x_xyyyzz[k] = -g_0_0_x_0_0_xyyyzz[k] - g_x_0_x_0_0_xyyyzz[k] * ab_x + g_x_0_x_0_0_xxyyyzz[k];

                g_x_0_x_0_x_xyyzzz[k] = -g_0_0_x_0_0_xyyzzz[k] - g_x_0_x_0_0_xyyzzz[k] * ab_x + g_x_0_x_0_0_xxyyzzz[k];

                g_x_0_x_0_x_xyzzzz[k] = -g_0_0_x_0_0_xyzzzz[k] - g_x_0_x_0_0_xyzzzz[k] * ab_x + g_x_0_x_0_0_xxyzzzz[k];

                g_x_0_x_0_x_xzzzzz[k] = -g_0_0_x_0_0_xzzzzz[k] - g_x_0_x_0_0_xzzzzz[k] * ab_x + g_x_0_x_0_0_xxzzzzz[k];

                g_x_0_x_0_x_yyyyyy[k] = -g_0_0_x_0_0_yyyyyy[k] - g_x_0_x_0_0_yyyyyy[k] * ab_x + g_x_0_x_0_0_xyyyyyy[k];

                g_x_0_x_0_x_yyyyyz[k] = -g_0_0_x_0_0_yyyyyz[k] - g_x_0_x_0_0_yyyyyz[k] * ab_x + g_x_0_x_0_0_xyyyyyz[k];

                g_x_0_x_0_x_yyyyzz[k] = -g_0_0_x_0_0_yyyyzz[k] - g_x_0_x_0_0_yyyyzz[k] * ab_x + g_x_0_x_0_0_xyyyyzz[k];

                g_x_0_x_0_x_yyyzzz[k] = -g_0_0_x_0_0_yyyzzz[k] - g_x_0_x_0_0_yyyzzz[k] * ab_x + g_x_0_x_0_0_xyyyzzz[k];

                g_x_0_x_0_x_yyzzzz[k] = -g_0_0_x_0_0_yyzzzz[k] - g_x_0_x_0_0_yyzzzz[k] * ab_x + g_x_0_x_0_0_xyyzzzz[k];

                g_x_0_x_0_x_yzzzzz[k] = -g_0_0_x_0_0_yzzzzz[k] - g_x_0_x_0_0_yzzzzz[k] * ab_x + g_x_0_x_0_0_xyzzzzz[k];

                g_x_0_x_0_x_zzzzzz[k] = -g_0_0_x_0_0_zzzzzz[k] - g_x_0_x_0_0_zzzzzz[k] * ab_x + g_x_0_x_0_0_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_y_xxxxxx = cbuffer.data(pi_geom_1010_off + 28 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxxxy = cbuffer.data(pi_geom_1010_off + 29 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxxxz = cbuffer.data(pi_geom_1010_off + 30 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxxyy = cbuffer.data(pi_geom_1010_off + 31 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxxyz = cbuffer.data(pi_geom_1010_off + 32 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxxzz = cbuffer.data(pi_geom_1010_off + 33 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxyyy = cbuffer.data(pi_geom_1010_off + 34 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxyyz = cbuffer.data(pi_geom_1010_off + 35 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxyzz = cbuffer.data(pi_geom_1010_off + 36 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxxzzz = cbuffer.data(pi_geom_1010_off + 37 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxyyyy = cbuffer.data(pi_geom_1010_off + 38 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxyyyz = cbuffer.data(pi_geom_1010_off + 39 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxyyzz = cbuffer.data(pi_geom_1010_off + 40 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxyzzz = cbuffer.data(pi_geom_1010_off + 41 * ccomps * dcomps);

            auto g_x_0_x_0_y_xxzzzz = cbuffer.data(pi_geom_1010_off + 42 * ccomps * dcomps);

            auto g_x_0_x_0_y_xyyyyy = cbuffer.data(pi_geom_1010_off + 43 * ccomps * dcomps);

            auto g_x_0_x_0_y_xyyyyz = cbuffer.data(pi_geom_1010_off + 44 * ccomps * dcomps);

            auto g_x_0_x_0_y_xyyyzz = cbuffer.data(pi_geom_1010_off + 45 * ccomps * dcomps);

            auto g_x_0_x_0_y_xyyzzz = cbuffer.data(pi_geom_1010_off + 46 * ccomps * dcomps);

            auto g_x_0_x_0_y_xyzzzz = cbuffer.data(pi_geom_1010_off + 47 * ccomps * dcomps);

            auto g_x_0_x_0_y_xzzzzz = cbuffer.data(pi_geom_1010_off + 48 * ccomps * dcomps);

            auto g_x_0_x_0_y_yyyyyy = cbuffer.data(pi_geom_1010_off + 49 * ccomps * dcomps);

            auto g_x_0_x_0_y_yyyyyz = cbuffer.data(pi_geom_1010_off + 50 * ccomps * dcomps);

            auto g_x_0_x_0_y_yyyyzz = cbuffer.data(pi_geom_1010_off + 51 * ccomps * dcomps);

            auto g_x_0_x_0_y_yyyzzz = cbuffer.data(pi_geom_1010_off + 52 * ccomps * dcomps);

            auto g_x_0_x_0_y_yyzzzz = cbuffer.data(pi_geom_1010_off + 53 * ccomps * dcomps);

            auto g_x_0_x_0_y_yzzzzz = cbuffer.data(pi_geom_1010_off + 54 * ccomps * dcomps);

            auto g_x_0_x_0_y_zzzzzz = cbuffer.data(pi_geom_1010_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_0_xxxxxx, g_x_0_x_0_0_xxxxxxy, g_x_0_x_0_0_xxxxxy, g_x_0_x_0_0_xxxxxyy, g_x_0_x_0_0_xxxxxyz, g_x_0_x_0_0_xxxxxz, g_x_0_x_0_0_xxxxyy, g_x_0_x_0_0_xxxxyyy, g_x_0_x_0_0_xxxxyyz, g_x_0_x_0_0_xxxxyz, g_x_0_x_0_0_xxxxyzz, g_x_0_x_0_0_xxxxzz, g_x_0_x_0_0_xxxyyy, g_x_0_x_0_0_xxxyyyy, g_x_0_x_0_0_xxxyyyz, g_x_0_x_0_0_xxxyyz, g_x_0_x_0_0_xxxyyzz, g_x_0_x_0_0_xxxyzz, g_x_0_x_0_0_xxxyzzz, g_x_0_x_0_0_xxxzzz, g_x_0_x_0_0_xxyyyy, g_x_0_x_0_0_xxyyyyy, g_x_0_x_0_0_xxyyyyz, g_x_0_x_0_0_xxyyyz, g_x_0_x_0_0_xxyyyzz, g_x_0_x_0_0_xxyyzz, g_x_0_x_0_0_xxyyzzz, g_x_0_x_0_0_xxyzzz, g_x_0_x_0_0_xxyzzzz, g_x_0_x_0_0_xxzzzz, g_x_0_x_0_0_xyyyyy, g_x_0_x_0_0_xyyyyyy, g_x_0_x_0_0_xyyyyyz, g_x_0_x_0_0_xyyyyz, g_x_0_x_0_0_xyyyyzz, g_x_0_x_0_0_xyyyzz, g_x_0_x_0_0_xyyyzzz, g_x_0_x_0_0_xyyzzz, g_x_0_x_0_0_xyyzzzz, g_x_0_x_0_0_xyzzzz, g_x_0_x_0_0_xyzzzzz, g_x_0_x_0_0_xzzzzz, g_x_0_x_0_0_yyyyyy, g_x_0_x_0_0_yyyyyyy, g_x_0_x_0_0_yyyyyyz, g_x_0_x_0_0_yyyyyz, g_x_0_x_0_0_yyyyyzz, g_x_0_x_0_0_yyyyzz, g_x_0_x_0_0_yyyyzzz, g_x_0_x_0_0_yyyzzz, g_x_0_x_0_0_yyyzzzz, g_x_0_x_0_0_yyzzzz, g_x_0_x_0_0_yyzzzzz, g_x_0_x_0_0_yzzzzz, g_x_0_x_0_0_yzzzzzz, g_x_0_x_0_0_zzzzzz, g_x_0_x_0_y_xxxxxx, g_x_0_x_0_y_xxxxxy, g_x_0_x_0_y_xxxxxz, g_x_0_x_0_y_xxxxyy, g_x_0_x_0_y_xxxxyz, g_x_0_x_0_y_xxxxzz, g_x_0_x_0_y_xxxyyy, g_x_0_x_0_y_xxxyyz, g_x_0_x_0_y_xxxyzz, g_x_0_x_0_y_xxxzzz, g_x_0_x_0_y_xxyyyy, g_x_0_x_0_y_xxyyyz, g_x_0_x_0_y_xxyyzz, g_x_0_x_0_y_xxyzzz, g_x_0_x_0_y_xxzzzz, g_x_0_x_0_y_xyyyyy, g_x_0_x_0_y_xyyyyz, g_x_0_x_0_y_xyyyzz, g_x_0_x_0_y_xyyzzz, g_x_0_x_0_y_xyzzzz, g_x_0_x_0_y_xzzzzz, g_x_0_x_0_y_yyyyyy, g_x_0_x_0_y_yyyyyz, g_x_0_x_0_y_yyyyzz, g_x_0_x_0_y_yyyzzz, g_x_0_x_0_y_yyzzzz, g_x_0_x_0_y_yzzzzz, g_x_0_x_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_y_xxxxxx[k] = -g_x_0_x_0_0_xxxxxx[k] * ab_y + g_x_0_x_0_0_xxxxxxy[k];

                g_x_0_x_0_y_xxxxxy[k] = -g_x_0_x_0_0_xxxxxy[k] * ab_y + g_x_0_x_0_0_xxxxxyy[k];

                g_x_0_x_0_y_xxxxxz[k] = -g_x_0_x_0_0_xxxxxz[k] * ab_y + g_x_0_x_0_0_xxxxxyz[k];

                g_x_0_x_0_y_xxxxyy[k] = -g_x_0_x_0_0_xxxxyy[k] * ab_y + g_x_0_x_0_0_xxxxyyy[k];

                g_x_0_x_0_y_xxxxyz[k] = -g_x_0_x_0_0_xxxxyz[k] * ab_y + g_x_0_x_0_0_xxxxyyz[k];

                g_x_0_x_0_y_xxxxzz[k] = -g_x_0_x_0_0_xxxxzz[k] * ab_y + g_x_0_x_0_0_xxxxyzz[k];

                g_x_0_x_0_y_xxxyyy[k] = -g_x_0_x_0_0_xxxyyy[k] * ab_y + g_x_0_x_0_0_xxxyyyy[k];

                g_x_0_x_0_y_xxxyyz[k] = -g_x_0_x_0_0_xxxyyz[k] * ab_y + g_x_0_x_0_0_xxxyyyz[k];

                g_x_0_x_0_y_xxxyzz[k] = -g_x_0_x_0_0_xxxyzz[k] * ab_y + g_x_0_x_0_0_xxxyyzz[k];

                g_x_0_x_0_y_xxxzzz[k] = -g_x_0_x_0_0_xxxzzz[k] * ab_y + g_x_0_x_0_0_xxxyzzz[k];

                g_x_0_x_0_y_xxyyyy[k] = -g_x_0_x_0_0_xxyyyy[k] * ab_y + g_x_0_x_0_0_xxyyyyy[k];

                g_x_0_x_0_y_xxyyyz[k] = -g_x_0_x_0_0_xxyyyz[k] * ab_y + g_x_0_x_0_0_xxyyyyz[k];

                g_x_0_x_0_y_xxyyzz[k] = -g_x_0_x_0_0_xxyyzz[k] * ab_y + g_x_0_x_0_0_xxyyyzz[k];

                g_x_0_x_0_y_xxyzzz[k] = -g_x_0_x_0_0_xxyzzz[k] * ab_y + g_x_0_x_0_0_xxyyzzz[k];

                g_x_0_x_0_y_xxzzzz[k] = -g_x_0_x_0_0_xxzzzz[k] * ab_y + g_x_0_x_0_0_xxyzzzz[k];

                g_x_0_x_0_y_xyyyyy[k] = -g_x_0_x_0_0_xyyyyy[k] * ab_y + g_x_0_x_0_0_xyyyyyy[k];

                g_x_0_x_0_y_xyyyyz[k] = -g_x_0_x_0_0_xyyyyz[k] * ab_y + g_x_0_x_0_0_xyyyyyz[k];

                g_x_0_x_0_y_xyyyzz[k] = -g_x_0_x_0_0_xyyyzz[k] * ab_y + g_x_0_x_0_0_xyyyyzz[k];

                g_x_0_x_0_y_xyyzzz[k] = -g_x_0_x_0_0_xyyzzz[k] * ab_y + g_x_0_x_0_0_xyyyzzz[k];

                g_x_0_x_0_y_xyzzzz[k] = -g_x_0_x_0_0_xyzzzz[k] * ab_y + g_x_0_x_0_0_xyyzzzz[k];

                g_x_0_x_0_y_xzzzzz[k] = -g_x_0_x_0_0_xzzzzz[k] * ab_y + g_x_0_x_0_0_xyzzzzz[k];

                g_x_0_x_0_y_yyyyyy[k] = -g_x_0_x_0_0_yyyyyy[k] * ab_y + g_x_0_x_0_0_yyyyyyy[k];

                g_x_0_x_0_y_yyyyyz[k] = -g_x_0_x_0_0_yyyyyz[k] * ab_y + g_x_0_x_0_0_yyyyyyz[k];

                g_x_0_x_0_y_yyyyzz[k] = -g_x_0_x_0_0_yyyyzz[k] * ab_y + g_x_0_x_0_0_yyyyyzz[k];

                g_x_0_x_0_y_yyyzzz[k] = -g_x_0_x_0_0_yyyzzz[k] * ab_y + g_x_0_x_0_0_yyyyzzz[k];

                g_x_0_x_0_y_yyzzzz[k] = -g_x_0_x_0_0_yyzzzz[k] * ab_y + g_x_0_x_0_0_yyyzzzz[k];

                g_x_0_x_0_y_yzzzzz[k] = -g_x_0_x_0_0_yzzzzz[k] * ab_y + g_x_0_x_0_0_yyzzzzz[k];

                g_x_0_x_0_y_zzzzzz[k] = -g_x_0_x_0_0_zzzzzz[k] * ab_y + g_x_0_x_0_0_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_x_0_z_xxxxxx = cbuffer.data(pi_geom_1010_off + 56 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxxxy = cbuffer.data(pi_geom_1010_off + 57 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxxxz = cbuffer.data(pi_geom_1010_off + 58 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxxyy = cbuffer.data(pi_geom_1010_off + 59 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxxyz = cbuffer.data(pi_geom_1010_off + 60 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxxzz = cbuffer.data(pi_geom_1010_off + 61 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxyyy = cbuffer.data(pi_geom_1010_off + 62 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxyyz = cbuffer.data(pi_geom_1010_off + 63 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxyzz = cbuffer.data(pi_geom_1010_off + 64 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxxzzz = cbuffer.data(pi_geom_1010_off + 65 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxyyyy = cbuffer.data(pi_geom_1010_off + 66 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxyyyz = cbuffer.data(pi_geom_1010_off + 67 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxyyzz = cbuffer.data(pi_geom_1010_off + 68 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxyzzz = cbuffer.data(pi_geom_1010_off + 69 * ccomps * dcomps);

            auto g_x_0_x_0_z_xxzzzz = cbuffer.data(pi_geom_1010_off + 70 * ccomps * dcomps);

            auto g_x_0_x_0_z_xyyyyy = cbuffer.data(pi_geom_1010_off + 71 * ccomps * dcomps);

            auto g_x_0_x_0_z_xyyyyz = cbuffer.data(pi_geom_1010_off + 72 * ccomps * dcomps);

            auto g_x_0_x_0_z_xyyyzz = cbuffer.data(pi_geom_1010_off + 73 * ccomps * dcomps);

            auto g_x_0_x_0_z_xyyzzz = cbuffer.data(pi_geom_1010_off + 74 * ccomps * dcomps);

            auto g_x_0_x_0_z_xyzzzz = cbuffer.data(pi_geom_1010_off + 75 * ccomps * dcomps);

            auto g_x_0_x_0_z_xzzzzz = cbuffer.data(pi_geom_1010_off + 76 * ccomps * dcomps);

            auto g_x_0_x_0_z_yyyyyy = cbuffer.data(pi_geom_1010_off + 77 * ccomps * dcomps);

            auto g_x_0_x_0_z_yyyyyz = cbuffer.data(pi_geom_1010_off + 78 * ccomps * dcomps);

            auto g_x_0_x_0_z_yyyyzz = cbuffer.data(pi_geom_1010_off + 79 * ccomps * dcomps);

            auto g_x_0_x_0_z_yyyzzz = cbuffer.data(pi_geom_1010_off + 80 * ccomps * dcomps);

            auto g_x_0_x_0_z_yyzzzz = cbuffer.data(pi_geom_1010_off + 81 * ccomps * dcomps);

            auto g_x_0_x_0_z_yzzzzz = cbuffer.data(pi_geom_1010_off + 82 * ccomps * dcomps);

            auto g_x_0_x_0_z_zzzzzz = cbuffer.data(pi_geom_1010_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_0_0_xxxxxx, g_x_0_x_0_0_xxxxxxz, g_x_0_x_0_0_xxxxxy, g_x_0_x_0_0_xxxxxyz, g_x_0_x_0_0_xxxxxz, g_x_0_x_0_0_xxxxxzz, g_x_0_x_0_0_xxxxyy, g_x_0_x_0_0_xxxxyyz, g_x_0_x_0_0_xxxxyz, g_x_0_x_0_0_xxxxyzz, g_x_0_x_0_0_xxxxzz, g_x_0_x_0_0_xxxxzzz, g_x_0_x_0_0_xxxyyy, g_x_0_x_0_0_xxxyyyz, g_x_0_x_0_0_xxxyyz, g_x_0_x_0_0_xxxyyzz, g_x_0_x_0_0_xxxyzz, g_x_0_x_0_0_xxxyzzz, g_x_0_x_0_0_xxxzzz, g_x_0_x_0_0_xxxzzzz, g_x_0_x_0_0_xxyyyy, g_x_0_x_0_0_xxyyyyz, g_x_0_x_0_0_xxyyyz, g_x_0_x_0_0_xxyyyzz, g_x_0_x_0_0_xxyyzz, g_x_0_x_0_0_xxyyzzz, g_x_0_x_0_0_xxyzzz, g_x_0_x_0_0_xxyzzzz, g_x_0_x_0_0_xxzzzz, g_x_0_x_0_0_xxzzzzz, g_x_0_x_0_0_xyyyyy, g_x_0_x_0_0_xyyyyyz, g_x_0_x_0_0_xyyyyz, g_x_0_x_0_0_xyyyyzz, g_x_0_x_0_0_xyyyzz, g_x_0_x_0_0_xyyyzzz, g_x_0_x_0_0_xyyzzz, g_x_0_x_0_0_xyyzzzz, g_x_0_x_0_0_xyzzzz, g_x_0_x_0_0_xyzzzzz, g_x_0_x_0_0_xzzzzz, g_x_0_x_0_0_xzzzzzz, g_x_0_x_0_0_yyyyyy, g_x_0_x_0_0_yyyyyyz, g_x_0_x_0_0_yyyyyz, g_x_0_x_0_0_yyyyyzz, g_x_0_x_0_0_yyyyzz, g_x_0_x_0_0_yyyyzzz, g_x_0_x_0_0_yyyzzz, g_x_0_x_0_0_yyyzzzz, g_x_0_x_0_0_yyzzzz, g_x_0_x_0_0_yyzzzzz, g_x_0_x_0_0_yzzzzz, g_x_0_x_0_0_yzzzzzz, g_x_0_x_0_0_zzzzzz, g_x_0_x_0_0_zzzzzzz, g_x_0_x_0_z_xxxxxx, g_x_0_x_0_z_xxxxxy, g_x_0_x_0_z_xxxxxz, g_x_0_x_0_z_xxxxyy, g_x_0_x_0_z_xxxxyz, g_x_0_x_0_z_xxxxzz, g_x_0_x_0_z_xxxyyy, g_x_0_x_0_z_xxxyyz, g_x_0_x_0_z_xxxyzz, g_x_0_x_0_z_xxxzzz, g_x_0_x_0_z_xxyyyy, g_x_0_x_0_z_xxyyyz, g_x_0_x_0_z_xxyyzz, g_x_0_x_0_z_xxyzzz, g_x_0_x_0_z_xxzzzz, g_x_0_x_0_z_xyyyyy, g_x_0_x_0_z_xyyyyz, g_x_0_x_0_z_xyyyzz, g_x_0_x_0_z_xyyzzz, g_x_0_x_0_z_xyzzzz, g_x_0_x_0_z_xzzzzz, g_x_0_x_0_z_yyyyyy, g_x_0_x_0_z_yyyyyz, g_x_0_x_0_z_yyyyzz, g_x_0_x_0_z_yyyzzz, g_x_0_x_0_z_yyzzzz, g_x_0_x_0_z_yzzzzz, g_x_0_x_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_x_0_z_xxxxxx[k] = -g_x_0_x_0_0_xxxxxx[k] * ab_z + g_x_0_x_0_0_xxxxxxz[k];

                g_x_0_x_0_z_xxxxxy[k] = -g_x_0_x_0_0_xxxxxy[k] * ab_z + g_x_0_x_0_0_xxxxxyz[k];

                g_x_0_x_0_z_xxxxxz[k] = -g_x_0_x_0_0_xxxxxz[k] * ab_z + g_x_0_x_0_0_xxxxxzz[k];

                g_x_0_x_0_z_xxxxyy[k] = -g_x_0_x_0_0_xxxxyy[k] * ab_z + g_x_0_x_0_0_xxxxyyz[k];

                g_x_0_x_0_z_xxxxyz[k] = -g_x_0_x_0_0_xxxxyz[k] * ab_z + g_x_0_x_0_0_xxxxyzz[k];

                g_x_0_x_0_z_xxxxzz[k] = -g_x_0_x_0_0_xxxxzz[k] * ab_z + g_x_0_x_0_0_xxxxzzz[k];

                g_x_0_x_0_z_xxxyyy[k] = -g_x_0_x_0_0_xxxyyy[k] * ab_z + g_x_0_x_0_0_xxxyyyz[k];

                g_x_0_x_0_z_xxxyyz[k] = -g_x_0_x_0_0_xxxyyz[k] * ab_z + g_x_0_x_0_0_xxxyyzz[k];

                g_x_0_x_0_z_xxxyzz[k] = -g_x_0_x_0_0_xxxyzz[k] * ab_z + g_x_0_x_0_0_xxxyzzz[k];

                g_x_0_x_0_z_xxxzzz[k] = -g_x_0_x_0_0_xxxzzz[k] * ab_z + g_x_0_x_0_0_xxxzzzz[k];

                g_x_0_x_0_z_xxyyyy[k] = -g_x_0_x_0_0_xxyyyy[k] * ab_z + g_x_0_x_0_0_xxyyyyz[k];

                g_x_0_x_0_z_xxyyyz[k] = -g_x_0_x_0_0_xxyyyz[k] * ab_z + g_x_0_x_0_0_xxyyyzz[k];

                g_x_0_x_0_z_xxyyzz[k] = -g_x_0_x_0_0_xxyyzz[k] * ab_z + g_x_0_x_0_0_xxyyzzz[k];

                g_x_0_x_0_z_xxyzzz[k] = -g_x_0_x_0_0_xxyzzz[k] * ab_z + g_x_0_x_0_0_xxyzzzz[k];

                g_x_0_x_0_z_xxzzzz[k] = -g_x_0_x_0_0_xxzzzz[k] * ab_z + g_x_0_x_0_0_xxzzzzz[k];

                g_x_0_x_0_z_xyyyyy[k] = -g_x_0_x_0_0_xyyyyy[k] * ab_z + g_x_0_x_0_0_xyyyyyz[k];

                g_x_0_x_0_z_xyyyyz[k] = -g_x_0_x_0_0_xyyyyz[k] * ab_z + g_x_0_x_0_0_xyyyyzz[k];

                g_x_0_x_0_z_xyyyzz[k] = -g_x_0_x_0_0_xyyyzz[k] * ab_z + g_x_0_x_0_0_xyyyzzz[k];

                g_x_0_x_0_z_xyyzzz[k] = -g_x_0_x_0_0_xyyzzz[k] * ab_z + g_x_0_x_0_0_xyyzzzz[k];

                g_x_0_x_0_z_xyzzzz[k] = -g_x_0_x_0_0_xyzzzz[k] * ab_z + g_x_0_x_0_0_xyzzzzz[k];

                g_x_0_x_0_z_xzzzzz[k] = -g_x_0_x_0_0_xzzzzz[k] * ab_z + g_x_0_x_0_0_xzzzzzz[k];

                g_x_0_x_0_z_yyyyyy[k] = -g_x_0_x_0_0_yyyyyy[k] * ab_z + g_x_0_x_0_0_yyyyyyz[k];

                g_x_0_x_0_z_yyyyyz[k] = -g_x_0_x_0_0_yyyyyz[k] * ab_z + g_x_0_x_0_0_yyyyyzz[k];

                g_x_0_x_0_z_yyyyzz[k] = -g_x_0_x_0_0_yyyyzz[k] * ab_z + g_x_0_x_0_0_yyyyzzz[k];

                g_x_0_x_0_z_yyyzzz[k] = -g_x_0_x_0_0_yyyzzz[k] * ab_z + g_x_0_x_0_0_yyyzzzz[k];

                g_x_0_x_0_z_yyzzzz[k] = -g_x_0_x_0_0_yyzzzz[k] * ab_z + g_x_0_x_0_0_yyzzzzz[k];

                g_x_0_x_0_z_yzzzzz[k] = -g_x_0_x_0_0_yzzzzz[k] * ab_z + g_x_0_x_0_0_yzzzzzz[k];

                g_x_0_x_0_z_zzzzzz[k] = -g_x_0_x_0_0_zzzzzz[k] * ab_z + g_x_0_x_0_0_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_x_xxxxxx = cbuffer.data(pi_geom_1010_off + 84 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxxxy = cbuffer.data(pi_geom_1010_off + 85 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxxxz = cbuffer.data(pi_geom_1010_off + 86 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxxyy = cbuffer.data(pi_geom_1010_off + 87 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxxyz = cbuffer.data(pi_geom_1010_off + 88 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxxzz = cbuffer.data(pi_geom_1010_off + 89 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxyyy = cbuffer.data(pi_geom_1010_off + 90 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxyyz = cbuffer.data(pi_geom_1010_off + 91 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxyzz = cbuffer.data(pi_geom_1010_off + 92 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxxzzz = cbuffer.data(pi_geom_1010_off + 93 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxyyyy = cbuffer.data(pi_geom_1010_off + 94 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxyyyz = cbuffer.data(pi_geom_1010_off + 95 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxyyzz = cbuffer.data(pi_geom_1010_off + 96 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxyzzz = cbuffer.data(pi_geom_1010_off + 97 * ccomps * dcomps);

            auto g_x_0_y_0_x_xxzzzz = cbuffer.data(pi_geom_1010_off + 98 * ccomps * dcomps);

            auto g_x_0_y_0_x_xyyyyy = cbuffer.data(pi_geom_1010_off + 99 * ccomps * dcomps);

            auto g_x_0_y_0_x_xyyyyz = cbuffer.data(pi_geom_1010_off + 100 * ccomps * dcomps);

            auto g_x_0_y_0_x_xyyyzz = cbuffer.data(pi_geom_1010_off + 101 * ccomps * dcomps);

            auto g_x_0_y_0_x_xyyzzz = cbuffer.data(pi_geom_1010_off + 102 * ccomps * dcomps);

            auto g_x_0_y_0_x_xyzzzz = cbuffer.data(pi_geom_1010_off + 103 * ccomps * dcomps);

            auto g_x_0_y_0_x_xzzzzz = cbuffer.data(pi_geom_1010_off + 104 * ccomps * dcomps);

            auto g_x_0_y_0_x_yyyyyy = cbuffer.data(pi_geom_1010_off + 105 * ccomps * dcomps);

            auto g_x_0_y_0_x_yyyyyz = cbuffer.data(pi_geom_1010_off + 106 * ccomps * dcomps);

            auto g_x_0_y_0_x_yyyyzz = cbuffer.data(pi_geom_1010_off + 107 * ccomps * dcomps);

            auto g_x_0_y_0_x_yyyzzz = cbuffer.data(pi_geom_1010_off + 108 * ccomps * dcomps);

            auto g_x_0_y_0_x_yyzzzz = cbuffer.data(pi_geom_1010_off + 109 * ccomps * dcomps);

            auto g_x_0_y_0_x_yzzzzz = cbuffer.data(pi_geom_1010_off + 110 * ccomps * dcomps);

            auto g_x_0_y_0_x_zzzzzz = cbuffer.data(pi_geom_1010_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_xxxxxx, g_0_0_y_0_0_xxxxxy, g_0_0_y_0_0_xxxxxz, g_0_0_y_0_0_xxxxyy, g_0_0_y_0_0_xxxxyz, g_0_0_y_0_0_xxxxzz, g_0_0_y_0_0_xxxyyy, g_0_0_y_0_0_xxxyyz, g_0_0_y_0_0_xxxyzz, g_0_0_y_0_0_xxxzzz, g_0_0_y_0_0_xxyyyy, g_0_0_y_0_0_xxyyyz, g_0_0_y_0_0_xxyyzz, g_0_0_y_0_0_xxyzzz, g_0_0_y_0_0_xxzzzz, g_0_0_y_0_0_xyyyyy, g_0_0_y_0_0_xyyyyz, g_0_0_y_0_0_xyyyzz, g_0_0_y_0_0_xyyzzz, g_0_0_y_0_0_xyzzzz, g_0_0_y_0_0_xzzzzz, g_0_0_y_0_0_yyyyyy, g_0_0_y_0_0_yyyyyz, g_0_0_y_0_0_yyyyzz, g_0_0_y_0_0_yyyzzz, g_0_0_y_0_0_yyzzzz, g_0_0_y_0_0_yzzzzz, g_0_0_y_0_0_zzzzzz, g_x_0_y_0_0_xxxxxx, g_x_0_y_0_0_xxxxxxx, g_x_0_y_0_0_xxxxxxy, g_x_0_y_0_0_xxxxxxz, g_x_0_y_0_0_xxxxxy, g_x_0_y_0_0_xxxxxyy, g_x_0_y_0_0_xxxxxyz, g_x_0_y_0_0_xxxxxz, g_x_0_y_0_0_xxxxxzz, g_x_0_y_0_0_xxxxyy, g_x_0_y_0_0_xxxxyyy, g_x_0_y_0_0_xxxxyyz, g_x_0_y_0_0_xxxxyz, g_x_0_y_0_0_xxxxyzz, g_x_0_y_0_0_xxxxzz, g_x_0_y_0_0_xxxxzzz, g_x_0_y_0_0_xxxyyy, g_x_0_y_0_0_xxxyyyy, g_x_0_y_0_0_xxxyyyz, g_x_0_y_0_0_xxxyyz, g_x_0_y_0_0_xxxyyzz, g_x_0_y_0_0_xxxyzz, g_x_0_y_0_0_xxxyzzz, g_x_0_y_0_0_xxxzzz, g_x_0_y_0_0_xxxzzzz, g_x_0_y_0_0_xxyyyy, g_x_0_y_0_0_xxyyyyy, g_x_0_y_0_0_xxyyyyz, g_x_0_y_0_0_xxyyyz, g_x_0_y_0_0_xxyyyzz, g_x_0_y_0_0_xxyyzz, g_x_0_y_0_0_xxyyzzz, g_x_0_y_0_0_xxyzzz, g_x_0_y_0_0_xxyzzzz, g_x_0_y_0_0_xxzzzz, g_x_0_y_0_0_xxzzzzz, g_x_0_y_0_0_xyyyyy, g_x_0_y_0_0_xyyyyyy, g_x_0_y_0_0_xyyyyyz, g_x_0_y_0_0_xyyyyz, g_x_0_y_0_0_xyyyyzz, g_x_0_y_0_0_xyyyzz, g_x_0_y_0_0_xyyyzzz, g_x_0_y_0_0_xyyzzz, g_x_0_y_0_0_xyyzzzz, g_x_0_y_0_0_xyzzzz, g_x_0_y_0_0_xyzzzzz, g_x_0_y_0_0_xzzzzz, g_x_0_y_0_0_xzzzzzz, g_x_0_y_0_0_yyyyyy, g_x_0_y_0_0_yyyyyz, g_x_0_y_0_0_yyyyzz, g_x_0_y_0_0_yyyzzz, g_x_0_y_0_0_yyzzzz, g_x_0_y_0_0_yzzzzz, g_x_0_y_0_0_zzzzzz, g_x_0_y_0_x_xxxxxx, g_x_0_y_0_x_xxxxxy, g_x_0_y_0_x_xxxxxz, g_x_0_y_0_x_xxxxyy, g_x_0_y_0_x_xxxxyz, g_x_0_y_0_x_xxxxzz, g_x_0_y_0_x_xxxyyy, g_x_0_y_0_x_xxxyyz, g_x_0_y_0_x_xxxyzz, g_x_0_y_0_x_xxxzzz, g_x_0_y_0_x_xxyyyy, g_x_0_y_0_x_xxyyyz, g_x_0_y_0_x_xxyyzz, g_x_0_y_0_x_xxyzzz, g_x_0_y_0_x_xxzzzz, g_x_0_y_0_x_xyyyyy, g_x_0_y_0_x_xyyyyz, g_x_0_y_0_x_xyyyzz, g_x_0_y_0_x_xyyzzz, g_x_0_y_0_x_xyzzzz, g_x_0_y_0_x_xzzzzz, g_x_0_y_0_x_yyyyyy, g_x_0_y_0_x_yyyyyz, g_x_0_y_0_x_yyyyzz, g_x_0_y_0_x_yyyzzz, g_x_0_y_0_x_yyzzzz, g_x_0_y_0_x_yzzzzz, g_x_0_y_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_x_xxxxxx[k] = -g_0_0_y_0_0_xxxxxx[k] - g_x_0_y_0_0_xxxxxx[k] * ab_x + g_x_0_y_0_0_xxxxxxx[k];

                g_x_0_y_0_x_xxxxxy[k] = -g_0_0_y_0_0_xxxxxy[k] - g_x_0_y_0_0_xxxxxy[k] * ab_x + g_x_0_y_0_0_xxxxxxy[k];

                g_x_0_y_0_x_xxxxxz[k] = -g_0_0_y_0_0_xxxxxz[k] - g_x_0_y_0_0_xxxxxz[k] * ab_x + g_x_0_y_0_0_xxxxxxz[k];

                g_x_0_y_0_x_xxxxyy[k] = -g_0_0_y_0_0_xxxxyy[k] - g_x_0_y_0_0_xxxxyy[k] * ab_x + g_x_0_y_0_0_xxxxxyy[k];

                g_x_0_y_0_x_xxxxyz[k] = -g_0_0_y_0_0_xxxxyz[k] - g_x_0_y_0_0_xxxxyz[k] * ab_x + g_x_0_y_0_0_xxxxxyz[k];

                g_x_0_y_0_x_xxxxzz[k] = -g_0_0_y_0_0_xxxxzz[k] - g_x_0_y_0_0_xxxxzz[k] * ab_x + g_x_0_y_0_0_xxxxxzz[k];

                g_x_0_y_0_x_xxxyyy[k] = -g_0_0_y_0_0_xxxyyy[k] - g_x_0_y_0_0_xxxyyy[k] * ab_x + g_x_0_y_0_0_xxxxyyy[k];

                g_x_0_y_0_x_xxxyyz[k] = -g_0_0_y_0_0_xxxyyz[k] - g_x_0_y_0_0_xxxyyz[k] * ab_x + g_x_0_y_0_0_xxxxyyz[k];

                g_x_0_y_0_x_xxxyzz[k] = -g_0_0_y_0_0_xxxyzz[k] - g_x_0_y_0_0_xxxyzz[k] * ab_x + g_x_0_y_0_0_xxxxyzz[k];

                g_x_0_y_0_x_xxxzzz[k] = -g_0_0_y_0_0_xxxzzz[k] - g_x_0_y_0_0_xxxzzz[k] * ab_x + g_x_0_y_0_0_xxxxzzz[k];

                g_x_0_y_0_x_xxyyyy[k] = -g_0_0_y_0_0_xxyyyy[k] - g_x_0_y_0_0_xxyyyy[k] * ab_x + g_x_0_y_0_0_xxxyyyy[k];

                g_x_0_y_0_x_xxyyyz[k] = -g_0_0_y_0_0_xxyyyz[k] - g_x_0_y_0_0_xxyyyz[k] * ab_x + g_x_0_y_0_0_xxxyyyz[k];

                g_x_0_y_0_x_xxyyzz[k] = -g_0_0_y_0_0_xxyyzz[k] - g_x_0_y_0_0_xxyyzz[k] * ab_x + g_x_0_y_0_0_xxxyyzz[k];

                g_x_0_y_0_x_xxyzzz[k] = -g_0_0_y_0_0_xxyzzz[k] - g_x_0_y_0_0_xxyzzz[k] * ab_x + g_x_0_y_0_0_xxxyzzz[k];

                g_x_0_y_0_x_xxzzzz[k] = -g_0_0_y_0_0_xxzzzz[k] - g_x_0_y_0_0_xxzzzz[k] * ab_x + g_x_0_y_0_0_xxxzzzz[k];

                g_x_0_y_0_x_xyyyyy[k] = -g_0_0_y_0_0_xyyyyy[k] - g_x_0_y_0_0_xyyyyy[k] * ab_x + g_x_0_y_0_0_xxyyyyy[k];

                g_x_0_y_0_x_xyyyyz[k] = -g_0_0_y_0_0_xyyyyz[k] - g_x_0_y_0_0_xyyyyz[k] * ab_x + g_x_0_y_0_0_xxyyyyz[k];

                g_x_0_y_0_x_xyyyzz[k] = -g_0_0_y_0_0_xyyyzz[k] - g_x_0_y_0_0_xyyyzz[k] * ab_x + g_x_0_y_0_0_xxyyyzz[k];

                g_x_0_y_0_x_xyyzzz[k] = -g_0_0_y_0_0_xyyzzz[k] - g_x_0_y_0_0_xyyzzz[k] * ab_x + g_x_0_y_0_0_xxyyzzz[k];

                g_x_0_y_0_x_xyzzzz[k] = -g_0_0_y_0_0_xyzzzz[k] - g_x_0_y_0_0_xyzzzz[k] * ab_x + g_x_0_y_0_0_xxyzzzz[k];

                g_x_0_y_0_x_xzzzzz[k] = -g_0_0_y_0_0_xzzzzz[k] - g_x_0_y_0_0_xzzzzz[k] * ab_x + g_x_0_y_0_0_xxzzzzz[k];

                g_x_0_y_0_x_yyyyyy[k] = -g_0_0_y_0_0_yyyyyy[k] - g_x_0_y_0_0_yyyyyy[k] * ab_x + g_x_0_y_0_0_xyyyyyy[k];

                g_x_0_y_0_x_yyyyyz[k] = -g_0_0_y_0_0_yyyyyz[k] - g_x_0_y_0_0_yyyyyz[k] * ab_x + g_x_0_y_0_0_xyyyyyz[k];

                g_x_0_y_0_x_yyyyzz[k] = -g_0_0_y_0_0_yyyyzz[k] - g_x_0_y_0_0_yyyyzz[k] * ab_x + g_x_0_y_0_0_xyyyyzz[k];

                g_x_0_y_0_x_yyyzzz[k] = -g_0_0_y_0_0_yyyzzz[k] - g_x_0_y_0_0_yyyzzz[k] * ab_x + g_x_0_y_0_0_xyyyzzz[k];

                g_x_0_y_0_x_yyzzzz[k] = -g_0_0_y_0_0_yyzzzz[k] - g_x_0_y_0_0_yyzzzz[k] * ab_x + g_x_0_y_0_0_xyyzzzz[k];

                g_x_0_y_0_x_yzzzzz[k] = -g_0_0_y_0_0_yzzzzz[k] - g_x_0_y_0_0_yzzzzz[k] * ab_x + g_x_0_y_0_0_xyzzzzz[k];

                g_x_0_y_0_x_zzzzzz[k] = -g_0_0_y_0_0_zzzzzz[k] - g_x_0_y_0_0_zzzzzz[k] * ab_x + g_x_0_y_0_0_xzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_y_xxxxxx = cbuffer.data(pi_geom_1010_off + 112 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxxxy = cbuffer.data(pi_geom_1010_off + 113 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxxxz = cbuffer.data(pi_geom_1010_off + 114 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxxyy = cbuffer.data(pi_geom_1010_off + 115 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxxyz = cbuffer.data(pi_geom_1010_off + 116 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxxzz = cbuffer.data(pi_geom_1010_off + 117 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxyyy = cbuffer.data(pi_geom_1010_off + 118 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxyyz = cbuffer.data(pi_geom_1010_off + 119 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxyzz = cbuffer.data(pi_geom_1010_off + 120 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxxzzz = cbuffer.data(pi_geom_1010_off + 121 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxyyyy = cbuffer.data(pi_geom_1010_off + 122 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxyyyz = cbuffer.data(pi_geom_1010_off + 123 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxyyzz = cbuffer.data(pi_geom_1010_off + 124 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxyzzz = cbuffer.data(pi_geom_1010_off + 125 * ccomps * dcomps);

            auto g_x_0_y_0_y_xxzzzz = cbuffer.data(pi_geom_1010_off + 126 * ccomps * dcomps);

            auto g_x_0_y_0_y_xyyyyy = cbuffer.data(pi_geom_1010_off + 127 * ccomps * dcomps);

            auto g_x_0_y_0_y_xyyyyz = cbuffer.data(pi_geom_1010_off + 128 * ccomps * dcomps);

            auto g_x_0_y_0_y_xyyyzz = cbuffer.data(pi_geom_1010_off + 129 * ccomps * dcomps);

            auto g_x_0_y_0_y_xyyzzz = cbuffer.data(pi_geom_1010_off + 130 * ccomps * dcomps);

            auto g_x_0_y_0_y_xyzzzz = cbuffer.data(pi_geom_1010_off + 131 * ccomps * dcomps);

            auto g_x_0_y_0_y_xzzzzz = cbuffer.data(pi_geom_1010_off + 132 * ccomps * dcomps);

            auto g_x_0_y_0_y_yyyyyy = cbuffer.data(pi_geom_1010_off + 133 * ccomps * dcomps);

            auto g_x_0_y_0_y_yyyyyz = cbuffer.data(pi_geom_1010_off + 134 * ccomps * dcomps);

            auto g_x_0_y_0_y_yyyyzz = cbuffer.data(pi_geom_1010_off + 135 * ccomps * dcomps);

            auto g_x_0_y_0_y_yyyzzz = cbuffer.data(pi_geom_1010_off + 136 * ccomps * dcomps);

            auto g_x_0_y_0_y_yyzzzz = cbuffer.data(pi_geom_1010_off + 137 * ccomps * dcomps);

            auto g_x_0_y_0_y_yzzzzz = cbuffer.data(pi_geom_1010_off + 138 * ccomps * dcomps);

            auto g_x_0_y_0_y_zzzzzz = cbuffer.data(pi_geom_1010_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_0_xxxxxx, g_x_0_y_0_0_xxxxxxy, g_x_0_y_0_0_xxxxxy, g_x_0_y_0_0_xxxxxyy, g_x_0_y_0_0_xxxxxyz, g_x_0_y_0_0_xxxxxz, g_x_0_y_0_0_xxxxyy, g_x_0_y_0_0_xxxxyyy, g_x_0_y_0_0_xxxxyyz, g_x_0_y_0_0_xxxxyz, g_x_0_y_0_0_xxxxyzz, g_x_0_y_0_0_xxxxzz, g_x_0_y_0_0_xxxyyy, g_x_0_y_0_0_xxxyyyy, g_x_0_y_0_0_xxxyyyz, g_x_0_y_0_0_xxxyyz, g_x_0_y_0_0_xxxyyzz, g_x_0_y_0_0_xxxyzz, g_x_0_y_0_0_xxxyzzz, g_x_0_y_0_0_xxxzzz, g_x_0_y_0_0_xxyyyy, g_x_0_y_0_0_xxyyyyy, g_x_0_y_0_0_xxyyyyz, g_x_0_y_0_0_xxyyyz, g_x_0_y_0_0_xxyyyzz, g_x_0_y_0_0_xxyyzz, g_x_0_y_0_0_xxyyzzz, g_x_0_y_0_0_xxyzzz, g_x_0_y_0_0_xxyzzzz, g_x_0_y_0_0_xxzzzz, g_x_0_y_0_0_xyyyyy, g_x_0_y_0_0_xyyyyyy, g_x_0_y_0_0_xyyyyyz, g_x_0_y_0_0_xyyyyz, g_x_0_y_0_0_xyyyyzz, g_x_0_y_0_0_xyyyzz, g_x_0_y_0_0_xyyyzzz, g_x_0_y_0_0_xyyzzz, g_x_0_y_0_0_xyyzzzz, g_x_0_y_0_0_xyzzzz, g_x_0_y_0_0_xyzzzzz, g_x_0_y_0_0_xzzzzz, g_x_0_y_0_0_yyyyyy, g_x_0_y_0_0_yyyyyyy, g_x_0_y_0_0_yyyyyyz, g_x_0_y_0_0_yyyyyz, g_x_0_y_0_0_yyyyyzz, g_x_0_y_0_0_yyyyzz, g_x_0_y_0_0_yyyyzzz, g_x_0_y_0_0_yyyzzz, g_x_0_y_0_0_yyyzzzz, g_x_0_y_0_0_yyzzzz, g_x_0_y_0_0_yyzzzzz, g_x_0_y_0_0_yzzzzz, g_x_0_y_0_0_yzzzzzz, g_x_0_y_0_0_zzzzzz, g_x_0_y_0_y_xxxxxx, g_x_0_y_0_y_xxxxxy, g_x_0_y_0_y_xxxxxz, g_x_0_y_0_y_xxxxyy, g_x_0_y_0_y_xxxxyz, g_x_0_y_0_y_xxxxzz, g_x_0_y_0_y_xxxyyy, g_x_0_y_0_y_xxxyyz, g_x_0_y_0_y_xxxyzz, g_x_0_y_0_y_xxxzzz, g_x_0_y_0_y_xxyyyy, g_x_0_y_0_y_xxyyyz, g_x_0_y_0_y_xxyyzz, g_x_0_y_0_y_xxyzzz, g_x_0_y_0_y_xxzzzz, g_x_0_y_0_y_xyyyyy, g_x_0_y_0_y_xyyyyz, g_x_0_y_0_y_xyyyzz, g_x_0_y_0_y_xyyzzz, g_x_0_y_0_y_xyzzzz, g_x_0_y_0_y_xzzzzz, g_x_0_y_0_y_yyyyyy, g_x_0_y_0_y_yyyyyz, g_x_0_y_0_y_yyyyzz, g_x_0_y_0_y_yyyzzz, g_x_0_y_0_y_yyzzzz, g_x_0_y_0_y_yzzzzz, g_x_0_y_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_y_xxxxxx[k] = -g_x_0_y_0_0_xxxxxx[k] * ab_y + g_x_0_y_0_0_xxxxxxy[k];

                g_x_0_y_0_y_xxxxxy[k] = -g_x_0_y_0_0_xxxxxy[k] * ab_y + g_x_0_y_0_0_xxxxxyy[k];

                g_x_0_y_0_y_xxxxxz[k] = -g_x_0_y_0_0_xxxxxz[k] * ab_y + g_x_0_y_0_0_xxxxxyz[k];

                g_x_0_y_0_y_xxxxyy[k] = -g_x_0_y_0_0_xxxxyy[k] * ab_y + g_x_0_y_0_0_xxxxyyy[k];

                g_x_0_y_0_y_xxxxyz[k] = -g_x_0_y_0_0_xxxxyz[k] * ab_y + g_x_0_y_0_0_xxxxyyz[k];

                g_x_0_y_0_y_xxxxzz[k] = -g_x_0_y_0_0_xxxxzz[k] * ab_y + g_x_0_y_0_0_xxxxyzz[k];

                g_x_0_y_0_y_xxxyyy[k] = -g_x_0_y_0_0_xxxyyy[k] * ab_y + g_x_0_y_0_0_xxxyyyy[k];

                g_x_0_y_0_y_xxxyyz[k] = -g_x_0_y_0_0_xxxyyz[k] * ab_y + g_x_0_y_0_0_xxxyyyz[k];

                g_x_0_y_0_y_xxxyzz[k] = -g_x_0_y_0_0_xxxyzz[k] * ab_y + g_x_0_y_0_0_xxxyyzz[k];

                g_x_0_y_0_y_xxxzzz[k] = -g_x_0_y_0_0_xxxzzz[k] * ab_y + g_x_0_y_0_0_xxxyzzz[k];

                g_x_0_y_0_y_xxyyyy[k] = -g_x_0_y_0_0_xxyyyy[k] * ab_y + g_x_0_y_0_0_xxyyyyy[k];

                g_x_0_y_0_y_xxyyyz[k] = -g_x_0_y_0_0_xxyyyz[k] * ab_y + g_x_0_y_0_0_xxyyyyz[k];

                g_x_0_y_0_y_xxyyzz[k] = -g_x_0_y_0_0_xxyyzz[k] * ab_y + g_x_0_y_0_0_xxyyyzz[k];

                g_x_0_y_0_y_xxyzzz[k] = -g_x_0_y_0_0_xxyzzz[k] * ab_y + g_x_0_y_0_0_xxyyzzz[k];

                g_x_0_y_0_y_xxzzzz[k] = -g_x_0_y_0_0_xxzzzz[k] * ab_y + g_x_0_y_0_0_xxyzzzz[k];

                g_x_0_y_0_y_xyyyyy[k] = -g_x_0_y_0_0_xyyyyy[k] * ab_y + g_x_0_y_0_0_xyyyyyy[k];

                g_x_0_y_0_y_xyyyyz[k] = -g_x_0_y_0_0_xyyyyz[k] * ab_y + g_x_0_y_0_0_xyyyyyz[k];

                g_x_0_y_0_y_xyyyzz[k] = -g_x_0_y_0_0_xyyyzz[k] * ab_y + g_x_0_y_0_0_xyyyyzz[k];

                g_x_0_y_0_y_xyyzzz[k] = -g_x_0_y_0_0_xyyzzz[k] * ab_y + g_x_0_y_0_0_xyyyzzz[k];

                g_x_0_y_0_y_xyzzzz[k] = -g_x_0_y_0_0_xyzzzz[k] * ab_y + g_x_0_y_0_0_xyyzzzz[k];

                g_x_0_y_0_y_xzzzzz[k] = -g_x_0_y_0_0_xzzzzz[k] * ab_y + g_x_0_y_0_0_xyzzzzz[k];

                g_x_0_y_0_y_yyyyyy[k] = -g_x_0_y_0_0_yyyyyy[k] * ab_y + g_x_0_y_0_0_yyyyyyy[k];

                g_x_0_y_0_y_yyyyyz[k] = -g_x_0_y_0_0_yyyyyz[k] * ab_y + g_x_0_y_0_0_yyyyyyz[k];

                g_x_0_y_0_y_yyyyzz[k] = -g_x_0_y_0_0_yyyyzz[k] * ab_y + g_x_0_y_0_0_yyyyyzz[k];

                g_x_0_y_0_y_yyyzzz[k] = -g_x_0_y_0_0_yyyzzz[k] * ab_y + g_x_0_y_0_0_yyyyzzz[k];

                g_x_0_y_0_y_yyzzzz[k] = -g_x_0_y_0_0_yyzzzz[k] * ab_y + g_x_0_y_0_0_yyyzzzz[k];

                g_x_0_y_0_y_yzzzzz[k] = -g_x_0_y_0_0_yzzzzz[k] * ab_y + g_x_0_y_0_0_yyzzzzz[k];

                g_x_0_y_0_y_zzzzzz[k] = -g_x_0_y_0_0_zzzzzz[k] * ab_y + g_x_0_y_0_0_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_y_0_z_xxxxxx = cbuffer.data(pi_geom_1010_off + 140 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxxxy = cbuffer.data(pi_geom_1010_off + 141 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxxxz = cbuffer.data(pi_geom_1010_off + 142 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxxyy = cbuffer.data(pi_geom_1010_off + 143 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxxyz = cbuffer.data(pi_geom_1010_off + 144 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxxzz = cbuffer.data(pi_geom_1010_off + 145 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxyyy = cbuffer.data(pi_geom_1010_off + 146 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxyyz = cbuffer.data(pi_geom_1010_off + 147 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxyzz = cbuffer.data(pi_geom_1010_off + 148 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxxzzz = cbuffer.data(pi_geom_1010_off + 149 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxyyyy = cbuffer.data(pi_geom_1010_off + 150 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxyyyz = cbuffer.data(pi_geom_1010_off + 151 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxyyzz = cbuffer.data(pi_geom_1010_off + 152 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxyzzz = cbuffer.data(pi_geom_1010_off + 153 * ccomps * dcomps);

            auto g_x_0_y_0_z_xxzzzz = cbuffer.data(pi_geom_1010_off + 154 * ccomps * dcomps);

            auto g_x_0_y_0_z_xyyyyy = cbuffer.data(pi_geom_1010_off + 155 * ccomps * dcomps);

            auto g_x_0_y_0_z_xyyyyz = cbuffer.data(pi_geom_1010_off + 156 * ccomps * dcomps);

            auto g_x_0_y_0_z_xyyyzz = cbuffer.data(pi_geom_1010_off + 157 * ccomps * dcomps);

            auto g_x_0_y_0_z_xyyzzz = cbuffer.data(pi_geom_1010_off + 158 * ccomps * dcomps);

            auto g_x_0_y_0_z_xyzzzz = cbuffer.data(pi_geom_1010_off + 159 * ccomps * dcomps);

            auto g_x_0_y_0_z_xzzzzz = cbuffer.data(pi_geom_1010_off + 160 * ccomps * dcomps);

            auto g_x_0_y_0_z_yyyyyy = cbuffer.data(pi_geom_1010_off + 161 * ccomps * dcomps);

            auto g_x_0_y_0_z_yyyyyz = cbuffer.data(pi_geom_1010_off + 162 * ccomps * dcomps);

            auto g_x_0_y_0_z_yyyyzz = cbuffer.data(pi_geom_1010_off + 163 * ccomps * dcomps);

            auto g_x_0_y_0_z_yyyzzz = cbuffer.data(pi_geom_1010_off + 164 * ccomps * dcomps);

            auto g_x_0_y_0_z_yyzzzz = cbuffer.data(pi_geom_1010_off + 165 * ccomps * dcomps);

            auto g_x_0_y_0_z_yzzzzz = cbuffer.data(pi_geom_1010_off + 166 * ccomps * dcomps);

            auto g_x_0_y_0_z_zzzzzz = cbuffer.data(pi_geom_1010_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_0_0_xxxxxx, g_x_0_y_0_0_xxxxxxz, g_x_0_y_0_0_xxxxxy, g_x_0_y_0_0_xxxxxyz, g_x_0_y_0_0_xxxxxz, g_x_0_y_0_0_xxxxxzz, g_x_0_y_0_0_xxxxyy, g_x_0_y_0_0_xxxxyyz, g_x_0_y_0_0_xxxxyz, g_x_0_y_0_0_xxxxyzz, g_x_0_y_0_0_xxxxzz, g_x_0_y_0_0_xxxxzzz, g_x_0_y_0_0_xxxyyy, g_x_0_y_0_0_xxxyyyz, g_x_0_y_0_0_xxxyyz, g_x_0_y_0_0_xxxyyzz, g_x_0_y_0_0_xxxyzz, g_x_0_y_0_0_xxxyzzz, g_x_0_y_0_0_xxxzzz, g_x_0_y_0_0_xxxzzzz, g_x_0_y_0_0_xxyyyy, g_x_0_y_0_0_xxyyyyz, g_x_0_y_0_0_xxyyyz, g_x_0_y_0_0_xxyyyzz, g_x_0_y_0_0_xxyyzz, g_x_0_y_0_0_xxyyzzz, g_x_0_y_0_0_xxyzzz, g_x_0_y_0_0_xxyzzzz, g_x_0_y_0_0_xxzzzz, g_x_0_y_0_0_xxzzzzz, g_x_0_y_0_0_xyyyyy, g_x_0_y_0_0_xyyyyyz, g_x_0_y_0_0_xyyyyz, g_x_0_y_0_0_xyyyyzz, g_x_0_y_0_0_xyyyzz, g_x_0_y_0_0_xyyyzzz, g_x_0_y_0_0_xyyzzz, g_x_0_y_0_0_xyyzzzz, g_x_0_y_0_0_xyzzzz, g_x_0_y_0_0_xyzzzzz, g_x_0_y_0_0_xzzzzz, g_x_0_y_0_0_xzzzzzz, g_x_0_y_0_0_yyyyyy, g_x_0_y_0_0_yyyyyyz, g_x_0_y_0_0_yyyyyz, g_x_0_y_0_0_yyyyyzz, g_x_0_y_0_0_yyyyzz, g_x_0_y_0_0_yyyyzzz, g_x_0_y_0_0_yyyzzz, g_x_0_y_0_0_yyyzzzz, g_x_0_y_0_0_yyzzzz, g_x_0_y_0_0_yyzzzzz, g_x_0_y_0_0_yzzzzz, g_x_0_y_0_0_yzzzzzz, g_x_0_y_0_0_zzzzzz, g_x_0_y_0_0_zzzzzzz, g_x_0_y_0_z_xxxxxx, g_x_0_y_0_z_xxxxxy, g_x_0_y_0_z_xxxxxz, g_x_0_y_0_z_xxxxyy, g_x_0_y_0_z_xxxxyz, g_x_0_y_0_z_xxxxzz, g_x_0_y_0_z_xxxyyy, g_x_0_y_0_z_xxxyyz, g_x_0_y_0_z_xxxyzz, g_x_0_y_0_z_xxxzzz, g_x_0_y_0_z_xxyyyy, g_x_0_y_0_z_xxyyyz, g_x_0_y_0_z_xxyyzz, g_x_0_y_0_z_xxyzzz, g_x_0_y_0_z_xxzzzz, g_x_0_y_0_z_xyyyyy, g_x_0_y_0_z_xyyyyz, g_x_0_y_0_z_xyyyzz, g_x_0_y_0_z_xyyzzz, g_x_0_y_0_z_xyzzzz, g_x_0_y_0_z_xzzzzz, g_x_0_y_0_z_yyyyyy, g_x_0_y_0_z_yyyyyz, g_x_0_y_0_z_yyyyzz, g_x_0_y_0_z_yyyzzz, g_x_0_y_0_z_yyzzzz, g_x_0_y_0_z_yzzzzz, g_x_0_y_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_y_0_z_xxxxxx[k] = -g_x_0_y_0_0_xxxxxx[k] * ab_z + g_x_0_y_0_0_xxxxxxz[k];

                g_x_0_y_0_z_xxxxxy[k] = -g_x_0_y_0_0_xxxxxy[k] * ab_z + g_x_0_y_0_0_xxxxxyz[k];

                g_x_0_y_0_z_xxxxxz[k] = -g_x_0_y_0_0_xxxxxz[k] * ab_z + g_x_0_y_0_0_xxxxxzz[k];

                g_x_0_y_0_z_xxxxyy[k] = -g_x_0_y_0_0_xxxxyy[k] * ab_z + g_x_0_y_0_0_xxxxyyz[k];

                g_x_0_y_0_z_xxxxyz[k] = -g_x_0_y_0_0_xxxxyz[k] * ab_z + g_x_0_y_0_0_xxxxyzz[k];

                g_x_0_y_0_z_xxxxzz[k] = -g_x_0_y_0_0_xxxxzz[k] * ab_z + g_x_0_y_0_0_xxxxzzz[k];

                g_x_0_y_0_z_xxxyyy[k] = -g_x_0_y_0_0_xxxyyy[k] * ab_z + g_x_0_y_0_0_xxxyyyz[k];

                g_x_0_y_0_z_xxxyyz[k] = -g_x_0_y_0_0_xxxyyz[k] * ab_z + g_x_0_y_0_0_xxxyyzz[k];

                g_x_0_y_0_z_xxxyzz[k] = -g_x_0_y_0_0_xxxyzz[k] * ab_z + g_x_0_y_0_0_xxxyzzz[k];

                g_x_0_y_0_z_xxxzzz[k] = -g_x_0_y_0_0_xxxzzz[k] * ab_z + g_x_0_y_0_0_xxxzzzz[k];

                g_x_0_y_0_z_xxyyyy[k] = -g_x_0_y_0_0_xxyyyy[k] * ab_z + g_x_0_y_0_0_xxyyyyz[k];

                g_x_0_y_0_z_xxyyyz[k] = -g_x_0_y_0_0_xxyyyz[k] * ab_z + g_x_0_y_0_0_xxyyyzz[k];

                g_x_0_y_0_z_xxyyzz[k] = -g_x_0_y_0_0_xxyyzz[k] * ab_z + g_x_0_y_0_0_xxyyzzz[k];

                g_x_0_y_0_z_xxyzzz[k] = -g_x_0_y_0_0_xxyzzz[k] * ab_z + g_x_0_y_0_0_xxyzzzz[k];

                g_x_0_y_0_z_xxzzzz[k] = -g_x_0_y_0_0_xxzzzz[k] * ab_z + g_x_0_y_0_0_xxzzzzz[k];

                g_x_0_y_0_z_xyyyyy[k] = -g_x_0_y_0_0_xyyyyy[k] * ab_z + g_x_0_y_0_0_xyyyyyz[k];

                g_x_0_y_0_z_xyyyyz[k] = -g_x_0_y_0_0_xyyyyz[k] * ab_z + g_x_0_y_0_0_xyyyyzz[k];

                g_x_0_y_0_z_xyyyzz[k] = -g_x_0_y_0_0_xyyyzz[k] * ab_z + g_x_0_y_0_0_xyyyzzz[k];

                g_x_0_y_0_z_xyyzzz[k] = -g_x_0_y_0_0_xyyzzz[k] * ab_z + g_x_0_y_0_0_xyyzzzz[k];

                g_x_0_y_0_z_xyzzzz[k] = -g_x_0_y_0_0_xyzzzz[k] * ab_z + g_x_0_y_0_0_xyzzzzz[k];

                g_x_0_y_0_z_xzzzzz[k] = -g_x_0_y_0_0_xzzzzz[k] * ab_z + g_x_0_y_0_0_xzzzzzz[k];

                g_x_0_y_0_z_yyyyyy[k] = -g_x_0_y_0_0_yyyyyy[k] * ab_z + g_x_0_y_0_0_yyyyyyz[k];

                g_x_0_y_0_z_yyyyyz[k] = -g_x_0_y_0_0_yyyyyz[k] * ab_z + g_x_0_y_0_0_yyyyyzz[k];

                g_x_0_y_0_z_yyyyzz[k] = -g_x_0_y_0_0_yyyyzz[k] * ab_z + g_x_0_y_0_0_yyyyzzz[k];

                g_x_0_y_0_z_yyyzzz[k] = -g_x_0_y_0_0_yyyzzz[k] * ab_z + g_x_0_y_0_0_yyyzzzz[k];

                g_x_0_y_0_z_yyzzzz[k] = -g_x_0_y_0_0_yyzzzz[k] * ab_z + g_x_0_y_0_0_yyzzzzz[k];

                g_x_0_y_0_z_yzzzzz[k] = -g_x_0_y_0_0_yzzzzz[k] * ab_z + g_x_0_y_0_0_yzzzzzz[k];

                g_x_0_y_0_z_zzzzzz[k] = -g_x_0_y_0_0_zzzzzz[k] * ab_z + g_x_0_y_0_0_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_x_xxxxxx = cbuffer.data(pi_geom_1010_off + 168 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxxxy = cbuffer.data(pi_geom_1010_off + 169 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxxxz = cbuffer.data(pi_geom_1010_off + 170 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxxyy = cbuffer.data(pi_geom_1010_off + 171 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxxyz = cbuffer.data(pi_geom_1010_off + 172 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxxzz = cbuffer.data(pi_geom_1010_off + 173 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxyyy = cbuffer.data(pi_geom_1010_off + 174 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxyyz = cbuffer.data(pi_geom_1010_off + 175 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxyzz = cbuffer.data(pi_geom_1010_off + 176 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxxzzz = cbuffer.data(pi_geom_1010_off + 177 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxyyyy = cbuffer.data(pi_geom_1010_off + 178 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxyyyz = cbuffer.data(pi_geom_1010_off + 179 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxyyzz = cbuffer.data(pi_geom_1010_off + 180 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxyzzz = cbuffer.data(pi_geom_1010_off + 181 * ccomps * dcomps);

            auto g_x_0_z_0_x_xxzzzz = cbuffer.data(pi_geom_1010_off + 182 * ccomps * dcomps);

            auto g_x_0_z_0_x_xyyyyy = cbuffer.data(pi_geom_1010_off + 183 * ccomps * dcomps);

            auto g_x_0_z_0_x_xyyyyz = cbuffer.data(pi_geom_1010_off + 184 * ccomps * dcomps);

            auto g_x_0_z_0_x_xyyyzz = cbuffer.data(pi_geom_1010_off + 185 * ccomps * dcomps);

            auto g_x_0_z_0_x_xyyzzz = cbuffer.data(pi_geom_1010_off + 186 * ccomps * dcomps);

            auto g_x_0_z_0_x_xyzzzz = cbuffer.data(pi_geom_1010_off + 187 * ccomps * dcomps);

            auto g_x_0_z_0_x_xzzzzz = cbuffer.data(pi_geom_1010_off + 188 * ccomps * dcomps);

            auto g_x_0_z_0_x_yyyyyy = cbuffer.data(pi_geom_1010_off + 189 * ccomps * dcomps);

            auto g_x_0_z_0_x_yyyyyz = cbuffer.data(pi_geom_1010_off + 190 * ccomps * dcomps);

            auto g_x_0_z_0_x_yyyyzz = cbuffer.data(pi_geom_1010_off + 191 * ccomps * dcomps);

            auto g_x_0_z_0_x_yyyzzz = cbuffer.data(pi_geom_1010_off + 192 * ccomps * dcomps);

            auto g_x_0_z_0_x_yyzzzz = cbuffer.data(pi_geom_1010_off + 193 * ccomps * dcomps);

            auto g_x_0_z_0_x_yzzzzz = cbuffer.data(pi_geom_1010_off + 194 * ccomps * dcomps);

            auto g_x_0_z_0_x_zzzzzz = cbuffer.data(pi_geom_1010_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_xxxxxx, g_0_0_z_0_0_xxxxxy, g_0_0_z_0_0_xxxxxz, g_0_0_z_0_0_xxxxyy, g_0_0_z_0_0_xxxxyz, g_0_0_z_0_0_xxxxzz, g_0_0_z_0_0_xxxyyy, g_0_0_z_0_0_xxxyyz, g_0_0_z_0_0_xxxyzz, g_0_0_z_0_0_xxxzzz, g_0_0_z_0_0_xxyyyy, g_0_0_z_0_0_xxyyyz, g_0_0_z_0_0_xxyyzz, g_0_0_z_0_0_xxyzzz, g_0_0_z_0_0_xxzzzz, g_0_0_z_0_0_xyyyyy, g_0_0_z_0_0_xyyyyz, g_0_0_z_0_0_xyyyzz, g_0_0_z_0_0_xyyzzz, g_0_0_z_0_0_xyzzzz, g_0_0_z_0_0_xzzzzz, g_0_0_z_0_0_yyyyyy, g_0_0_z_0_0_yyyyyz, g_0_0_z_0_0_yyyyzz, g_0_0_z_0_0_yyyzzz, g_0_0_z_0_0_yyzzzz, g_0_0_z_0_0_yzzzzz, g_0_0_z_0_0_zzzzzz, g_x_0_z_0_0_xxxxxx, g_x_0_z_0_0_xxxxxxx, g_x_0_z_0_0_xxxxxxy, g_x_0_z_0_0_xxxxxxz, g_x_0_z_0_0_xxxxxy, g_x_0_z_0_0_xxxxxyy, g_x_0_z_0_0_xxxxxyz, g_x_0_z_0_0_xxxxxz, g_x_0_z_0_0_xxxxxzz, g_x_0_z_0_0_xxxxyy, g_x_0_z_0_0_xxxxyyy, g_x_0_z_0_0_xxxxyyz, g_x_0_z_0_0_xxxxyz, g_x_0_z_0_0_xxxxyzz, g_x_0_z_0_0_xxxxzz, g_x_0_z_0_0_xxxxzzz, g_x_0_z_0_0_xxxyyy, g_x_0_z_0_0_xxxyyyy, g_x_0_z_0_0_xxxyyyz, g_x_0_z_0_0_xxxyyz, g_x_0_z_0_0_xxxyyzz, g_x_0_z_0_0_xxxyzz, g_x_0_z_0_0_xxxyzzz, g_x_0_z_0_0_xxxzzz, g_x_0_z_0_0_xxxzzzz, g_x_0_z_0_0_xxyyyy, g_x_0_z_0_0_xxyyyyy, g_x_0_z_0_0_xxyyyyz, g_x_0_z_0_0_xxyyyz, g_x_0_z_0_0_xxyyyzz, g_x_0_z_0_0_xxyyzz, g_x_0_z_0_0_xxyyzzz, g_x_0_z_0_0_xxyzzz, g_x_0_z_0_0_xxyzzzz, g_x_0_z_0_0_xxzzzz, g_x_0_z_0_0_xxzzzzz, g_x_0_z_0_0_xyyyyy, g_x_0_z_0_0_xyyyyyy, g_x_0_z_0_0_xyyyyyz, g_x_0_z_0_0_xyyyyz, g_x_0_z_0_0_xyyyyzz, g_x_0_z_0_0_xyyyzz, g_x_0_z_0_0_xyyyzzz, g_x_0_z_0_0_xyyzzz, g_x_0_z_0_0_xyyzzzz, g_x_0_z_0_0_xyzzzz, g_x_0_z_0_0_xyzzzzz, g_x_0_z_0_0_xzzzzz, g_x_0_z_0_0_xzzzzzz, g_x_0_z_0_0_yyyyyy, g_x_0_z_0_0_yyyyyz, g_x_0_z_0_0_yyyyzz, g_x_0_z_0_0_yyyzzz, g_x_0_z_0_0_yyzzzz, g_x_0_z_0_0_yzzzzz, g_x_0_z_0_0_zzzzzz, g_x_0_z_0_x_xxxxxx, g_x_0_z_0_x_xxxxxy, g_x_0_z_0_x_xxxxxz, g_x_0_z_0_x_xxxxyy, g_x_0_z_0_x_xxxxyz, g_x_0_z_0_x_xxxxzz, g_x_0_z_0_x_xxxyyy, g_x_0_z_0_x_xxxyyz, g_x_0_z_0_x_xxxyzz, g_x_0_z_0_x_xxxzzz, g_x_0_z_0_x_xxyyyy, g_x_0_z_0_x_xxyyyz, g_x_0_z_0_x_xxyyzz, g_x_0_z_0_x_xxyzzz, g_x_0_z_0_x_xxzzzz, g_x_0_z_0_x_xyyyyy, g_x_0_z_0_x_xyyyyz, g_x_0_z_0_x_xyyyzz, g_x_0_z_0_x_xyyzzz, g_x_0_z_0_x_xyzzzz, g_x_0_z_0_x_xzzzzz, g_x_0_z_0_x_yyyyyy, g_x_0_z_0_x_yyyyyz, g_x_0_z_0_x_yyyyzz, g_x_0_z_0_x_yyyzzz, g_x_0_z_0_x_yyzzzz, g_x_0_z_0_x_yzzzzz, g_x_0_z_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_x_xxxxxx[k] = -g_0_0_z_0_0_xxxxxx[k] - g_x_0_z_0_0_xxxxxx[k] * ab_x + g_x_0_z_0_0_xxxxxxx[k];

                g_x_0_z_0_x_xxxxxy[k] = -g_0_0_z_0_0_xxxxxy[k] - g_x_0_z_0_0_xxxxxy[k] * ab_x + g_x_0_z_0_0_xxxxxxy[k];

                g_x_0_z_0_x_xxxxxz[k] = -g_0_0_z_0_0_xxxxxz[k] - g_x_0_z_0_0_xxxxxz[k] * ab_x + g_x_0_z_0_0_xxxxxxz[k];

                g_x_0_z_0_x_xxxxyy[k] = -g_0_0_z_0_0_xxxxyy[k] - g_x_0_z_0_0_xxxxyy[k] * ab_x + g_x_0_z_0_0_xxxxxyy[k];

                g_x_0_z_0_x_xxxxyz[k] = -g_0_0_z_0_0_xxxxyz[k] - g_x_0_z_0_0_xxxxyz[k] * ab_x + g_x_0_z_0_0_xxxxxyz[k];

                g_x_0_z_0_x_xxxxzz[k] = -g_0_0_z_0_0_xxxxzz[k] - g_x_0_z_0_0_xxxxzz[k] * ab_x + g_x_0_z_0_0_xxxxxzz[k];

                g_x_0_z_0_x_xxxyyy[k] = -g_0_0_z_0_0_xxxyyy[k] - g_x_0_z_0_0_xxxyyy[k] * ab_x + g_x_0_z_0_0_xxxxyyy[k];

                g_x_0_z_0_x_xxxyyz[k] = -g_0_0_z_0_0_xxxyyz[k] - g_x_0_z_0_0_xxxyyz[k] * ab_x + g_x_0_z_0_0_xxxxyyz[k];

                g_x_0_z_0_x_xxxyzz[k] = -g_0_0_z_0_0_xxxyzz[k] - g_x_0_z_0_0_xxxyzz[k] * ab_x + g_x_0_z_0_0_xxxxyzz[k];

                g_x_0_z_0_x_xxxzzz[k] = -g_0_0_z_0_0_xxxzzz[k] - g_x_0_z_0_0_xxxzzz[k] * ab_x + g_x_0_z_0_0_xxxxzzz[k];

                g_x_0_z_0_x_xxyyyy[k] = -g_0_0_z_0_0_xxyyyy[k] - g_x_0_z_0_0_xxyyyy[k] * ab_x + g_x_0_z_0_0_xxxyyyy[k];

                g_x_0_z_0_x_xxyyyz[k] = -g_0_0_z_0_0_xxyyyz[k] - g_x_0_z_0_0_xxyyyz[k] * ab_x + g_x_0_z_0_0_xxxyyyz[k];

                g_x_0_z_0_x_xxyyzz[k] = -g_0_0_z_0_0_xxyyzz[k] - g_x_0_z_0_0_xxyyzz[k] * ab_x + g_x_0_z_0_0_xxxyyzz[k];

                g_x_0_z_0_x_xxyzzz[k] = -g_0_0_z_0_0_xxyzzz[k] - g_x_0_z_0_0_xxyzzz[k] * ab_x + g_x_0_z_0_0_xxxyzzz[k];

                g_x_0_z_0_x_xxzzzz[k] = -g_0_0_z_0_0_xxzzzz[k] - g_x_0_z_0_0_xxzzzz[k] * ab_x + g_x_0_z_0_0_xxxzzzz[k];

                g_x_0_z_0_x_xyyyyy[k] = -g_0_0_z_0_0_xyyyyy[k] - g_x_0_z_0_0_xyyyyy[k] * ab_x + g_x_0_z_0_0_xxyyyyy[k];

                g_x_0_z_0_x_xyyyyz[k] = -g_0_0_z_0_0_xyyyyz[k] - g_x_0_z_0_0_xyyyyz[k] * ab_x + g_x_0_z_0_0_xxyyyyz[k];

                g_x_0_z_0_x_xyyyzz[k] = -g_0_0_z_0_0_xyyyzz[k] - g_x_0_z_0_0_xyyyzz[k] * ab_x + g_x_0_z_0_0_xxyyyzz[k];

                g_x_0_z_0_x_xyyzzz[k] = -g_0_0_z_0_0_xyyzzz[k] - g_x_0_z_0_0_xyyzzz[k] * ab_x + g_x_0_z_0_0_xxyyzzz[k];

                g_x_0_z_0_x_xyzzzz[k] = -g_0_0_z_0_0_xyzzzz[k] - g_x_0_z_0_0_xyzzzz[k] * ab_x + g_x_0_z_0_0_xxyzzzz[k];

                g_x_0_z_0_x_xzzzzz[k] = -g_0_0_z_0_0_xzzzzz[k] - g_x_0_z_0_0_xzzzzz[k] * ab_x + g_x_0_z_0_0_xxzzzzz[k];

                g_x_0_z_0_x_yyyyyy[k] = -g_0_0_z_0_0_yyyyyy[k] - g_x_0_z_0_0_yyyyyy[k] * ab_x + g_x_0_z_0_0_xyyyyyy[k];

                g_x_0_z_0_x_yyyyyz[k] = -g_0_0_z_0_0_yyyyyz[k] - g_x_0_z_0_0_yyyyyz[k] * ab_x + g_x_0_z_0_0_xyyyyyz[k];

                g_x_0_z_0_x_yyyyzz[k] = -g_0_0_z_0_0_yyyyzz[k] - g_x_0_z_0_0_yyyyzz[k] * ab_x + g_x_0_z_0_0_xyyyyzz[k];

                g_x_0_z_0_x_yyyzzz[k] = -g_0_0_z_0_0_yyyzzz[k] - g_x_0_z_0_0_yyyzzz[k] * ab_x + g_x_0_z_0_0_xyyyzzz[k];

                g_x_0_z_0_x_yyzzzz[k] = -g_0_0_z_0_0_yyzzzz[k] - g_x_0_z_0_0_yyzzzz[k] * ab_x + g_x_0_z_0_0_xyyzzzz[k];

                g_x_0_z_0_x_yzzzzz[k] = -g_0_0_z_0_0_yzzzzz[k] - g_x_0_z_0_0_yzzzzz[k] * ab_x + g_x_0_z_0_0_xyzzzzz[k];

                g_x_0_z_0_x_zzzzzz[k] = -g_0_0_z_0_0_zzzzzz[k] - g_x_0_z_0_0_zzzzzz[k] * ab_x + g_x_0_z_0_0_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_y_xxxxxx = cbuffer.data(pi_geom_1010_off + 196 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxxxy = cbuffer.data(pi_geom_1010_off + 197 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxxxz = cbuffer.data(pi_geom_1010_off + 198 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxxyy = cbuffer.data(pi_geom_1010_off + 199 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxxyz = cbuffer.data(pi_geom_1010_off + 200 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxxzz = cbuffer.data(pi_geom_1010_off + 201 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxyyy = cbuffer.data(pi_geom_1010_off + 202 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxyyz = cbuffer.data(pi_geom_1010_off + 203 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxyzz = cbuffer.data(pi_geom_1010_off + 204 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxxzzz = cbuffer.data(pi_geom_1010_off + 205 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxyyyy = cbuffer.data(pi_geom_1010_off + 206 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxyyyz = cbuffer.data(pi_geom_1010_off + 207 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxyyzz = cbuffer.data(pi_geom_1010_off + 208 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxyzzz = cbuffer.data(pi_geom_1010_off + 209 * ccomps * dcomps);

            auto g_x_0_z_0_y_xxzzzz = cbuffer.data(pi_geom_1010_off + 210 * ccomps * dcomps);

            auto g_x_0_z_0_y_xyyyyy = cbuffer.data(pi_geom_1010_off + 211 * ccomps * dcomps);

            auto g_x_0_z_0_y_xyyyyz = cbuffer.data(pi_geom_1010_off + 212 * ccomps * dcomps);

            auto g_x_0_z_0_y_xyyyzz = cbuffer.data(pi_geom_1010_off + 213 * ccomps * dcomps);

            auto g_x_0_z_0_y_xyyzzz = cbuffer.data(pi_geom_1010_off + 214 * ccomps * dcomps);

            auto g_x_0_z_0_y_xyzzzz = cbuffer.data(pi_geom_1010_off + 215 * ccomps * dcomps);

            auto g_x_0_z_0_y_xzzzzz = cbuffer.data(pi_geom_1010_off + 216 * ccomps * dcomps);

            auto g_x_0_z_0_y_yyyyyy = cbuffer.data(pi_geom_1010_off + 217 * ccomps * dcomps);

            auto g_x_0_z_0_y_yyyyyz = cbuffer.data(pi_geom_1010_off + 218 * ccomps * dcomps);

            auto g_x_0_z_0_y_yyyyzz = cbuffer.data(pi_geom_1010_off + 219 * ccomps * dcomps);

            auto g_x_0_z_0_y_yyyzzz = cbuffer.data(pi_geom_1010_off + 220 * ccomps * dcomps);

            auto g_x_0_z_0_y_yyzzzz = cbuffer.data(pi_geom_1010_off + 221 * ccomps * dcomps);

            auto g_x_0_z_0_y_yzzzzz = cbuffer.data(pi_geom_1010_off + 222 * ccomps * dcomps);

            auto g_x_0_z_0_y_zzzzzz = cbuffer.data(pi_geom_1010_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_0_xxxxxx, g_x_0_z_0_0_xxxxxxy, g_x_0_z_0_0_xxxxxy, g_x_0_z_0_0_xxxxxyy, g_x_0_z_0_0_xxxxxyz, g_x_0_z_0_0_xxxxxz, g_x_0_z_0_0_xxxxyy, g_x_0_z_0_0_xxxxyyy, g_x_0_z_0_0_xxxxyyz, g_x_0_z_0_0_xxxxyz, g_x_0_z_0_0_xxxxyzz, g_x_0_z_0_0_xxxxzz, g_x_0_z_0_0_xxxyyy, g_x_0_z_0_0_xxxyyyy, g_x_0_z_0_0_xxxyyyz, g_x_0_z_0_0_xxxyyz, g_x_0_z_0_0_xxxyyzz, g_x_0_z_0_0_xxxyzz, g_x_0_z_0_0_xxxyzzz, g_x_0_z_0_0_xxxzzz, g_x_0_z_0_0_xxyyyy, g_x_0_z_0_0_xxyyyyy, g_x_0_z_0_0_xxyyyyz, g_x_0_z_0_0_xxyyyz, g_x_0_z_0_0_xxyyyzz, g_x_0_z_0_0_xxyyzz, g_x_0_z_0_0_xxyyzzz, g_x_0_z_0_0_xxyzzz, g_x_0_z_0_0_xxyzzzz, g_x_0_z_0_0_xxzzzz, g_x_0_z_0_0_xyyyyy, g_x_0_z_0_0_xyyyyyy, g_x_0_z_0_0_xyyyyyz, g_x_0_z_0_0_xyyyyz, g_x_0_z_0_0_xyyyyzz, g_x_0_z_0_0_xyyyzz, g_x_0_z_0_0_xyyyzzz, g_x_0_z_0_0_xyyzzz, g_x_0_z_0_0_xyyzzzz, g_x_0_z_0_0_xyzzzz, g_x_0_z_0_0_xyzzzzz, g_x_0_z_0_0_xzzzzz, g_x_0_z_0_0_yyyyyy, g_x_0_z_0_0_yyyyyyy, g_x_0_z_0_0_yyyyyyz, g_x_0_z_0_0_yyyyyz, g_x_0_z_0_0_yyyyyzz, g_x_0_z_0_0_yyyyzz, g_x_0_z_0_0_yyyyzzz, g_x_0_z_0_0_yyyzzz, g_x_0_z_0_0_yyyzzzz, g_x_0_z_0_0_yyzzzz, g_x_0_z_0_0_yyzzzzz, g_x_0_z_0_0_yzzzzz, g_x_0_z_0_0_yzzzzzz, g_x_0_z_0_0_zzzzzz, g_x_0_z_0_y_xxxxxx, g_x_0_z_0_y_xxxxxy, g_x_0_z_0_y_xxxxxz, g_x_0_z_0_y_xxxxyy, g_x_0_z_0_y_xxxxyz, g_x_0_z_0_y_xxxxzz, g_x_0_z_0_y_xxxyyy, g_x_0_z_0_y_xxxyyz, g_x_0_z_0_y_xxxyzz, g_x_0_z_0_y_xxxzzz, g_x_0_z_0_y_xxyyyy, g_x_0_z_0_y_xxyyyz, g_x_0_z_0_y_xxyyzz, g_x_0_z_0_y_xxyzzz, g_x_0_z_0_y_xxzzzz, g_x_0_z_0_y_xyyyyy, g_x_0_z_0_y_xyyyyz, g_x_0_z_0_y_xyyyzz, g_x_0_z_0_y_xyyzzz, g_x_0_z_0_y_xyzzzz, g_x_0_z_0_y_xzzzzz, g_x_0_z_0_y_yyyyyy, g_x_0_z_0_y_yyyyyz, g_x_0_z_0_y_yyyyzz, g_x_0_z_0_y_yyyzzz, g_x_0_z_0_y_yyzzzz, g_x_0_z_0_y_yzzzzz, g_x_0_z_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_y_xxxxxx[k] = -g_x_0_z_0_0_xxxxxx[k] * ab_y + g_x_0_z_0_0_xxxxxxy[k];

                g_x_0_z_0_y_xxxxxy[k] = -g_x_0_z_0_0_xxxxxy[k] * ab_y + g_x_0_z_0_0_xxxxxyy[k];

                g_x_0_z_0_y_xxxxxz[k] = -g_x_0_z_0_0_xxxxxz[k] * ab_y + g_x_0_z_0_0_xxxxxyz[k];

                g_x_0_z_0_y_xxxxyy[k] = -g_x_0_z_0_0_xxxxyy[k] * ab_y + g_x_0_z_0_0_xxxxyyy[k];

                g_x_0_z_0_y_xxxxyz[k] = -g_x_0_z_0_0_xxxxyz[k] * ab_y + g_x_0_z_0_0_xxxxyyz[k];

                g_x_0_z_0_y_xxxxzz[k] = -g_x_0_z_0_0_xxxxzz[k] * ab_y + g_x_0_z_0_0_xxxxyzz[k];

                g_x_0_z_0_y_xxxyyy[k] = -g_x_0_z_0_0_xxxyyy[k] * ab_y + g_x_0_z_0_0_xxxyyyy[k];

                g_x_0_z_0_y_xxxyyz[k] = -g_x_0_z_0_0_xxxyyz[k] * ab_y + g_x_0_z_0_0_xxxyyyz[k];

                g_x_0_z_0_y_xxxyzz[k] = -g_x_0_z_0_0_xxxyzz[k] * ab_y + g_x_0_z_0_0_xxxyyzz[k];

                g_x_0_z_0_y_xxxzzz[k] = -g_x_0_z_0_0_xxxzzz[k] * ab_y + g_x_0_z_0_0_xxxyzzz[k];

                g_x_0_z_0_y_xxyyyy[k] = -g_x_0_z_0_0_xxyyyy[k] * ab_y + g_x_0_z_0_0_xxyyyyy[k];

                g_x_0_z_0_y_xxyyyz[k] = -g_x_0_z_0_0_xxyyyz[k] * ab_y + g_x_0_z_0_0_xxyyyyz[k];

                g_x_0_z_0_y_xxyyzz[k] = -g_x_0_z_0_0_xxyyzz[k] * ab_y + g_x_0_z_0_0_xxyyyzz[k];

                g_x_0_z_0_y_xxyzzz[k] = -g_x_0_z_0_0_xxyzzz[k] * ab_y + g_x_0_z_0_0_xxyyzzz[k];

                g_x_0_z_0_y_xxzzzz[k] = -g_x_0_z_0_0_xxzzzz[k] * ab_y + g_x_0_z_0_0_xxyzzzz[k];

                g_x_0_z_0_y_xyyyyy[k] = -g_x_0_z_0_0_xyyyyy[k] * ab_y + g_x_0_z_0_0_xyyyyyy[k];

                g_x_0_z_0_y_xyyyyz[k] = -g_x_0_z_0_0_xyyyyz[k] * ab_y + g_x_0_z_0_0_xyyyyyz[k];

                g_x_0_z_0_y_xyyyzz[k] = -g_x_0_z_0_0_xyyyzz[k] * ab_y + g_x_0_z_0_0_xyyyyzz[k];

                g_x_0_z_0_y_xyyzzz[k] = -g_x_0_z_0_0_xyyzzz[k] * ab_y + g_x_0_z_0_0_xyyyzzz[k];

                g_x_0_z_0_y_xyzzzz[k] = -g_x_0_z_0_0_xyzzzz[k] * ab_y + g_x_0_z_0_0_xyyzzzz[k];

                g_x_0_z_0_y_xzzzzz[k] = -g_x_0_z_0_0_xzzzzz[k] * ab_y + g_x_0_z_0_0_xyzzzzz[k];

                g_x_0_z_0_y_yyyyyy[k] = -g_x_0_z_0_0_yyyyyy[k] * ab_y + g_x_0_z_0_0_yyyyyyy[k];

                g_x_0_z_0_y_yyyyyz[k] = -g_x_0_z_0_0_yyyyyz[k] * ab_y + g_x_0_z_0_0_yyyyyyz[k];

                g_x_0_z_0_y_yyyyzz[k] = -g_x_0_z_0_0_yyyyzz[k] * ab_y + g_x_0_z_0_0_yyyyyzz[k];

                g_x_0_z_0_y_yyyzzz[k] = -g_x_0_z_0_0_yyyzzz[k] * ab_y + g_x_0_z_0_0_yyyyzzz[k];

                g_x_0_z_0_y_yyzzzz[k] = -g_x_0_z_0_0_yyzzzz[k] * ab_y + g_x_0_z_0_0_yyyzzzz[k];

                g_x_0_z_0_y_yzzzzz[k] = -g_x_0_z_0_0_yzzzzz[k] * ab_y + g_x_0_z_0_0_yyzzzzz[k];

                g_x_0_z_0_y_zzzzzz[k] = -g_x_0_z_0_0_zzzzzz[k] * ab_y + g_x_0_z_0_0_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_z_0_z_xxxxxx = cbuffer.data(pi_geom_1010_off + 224 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxxxy = cbuffer.data(pi_geom_1010_off + 225 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxxxz = cbuffer.data(pi_geom_1010_off + 226 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxxyy = cbuffer.data(pi_geom_1010_off + 227 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxxyz = cbuffer.data(pi_geom_1010_off + 228 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxxzz = cbuffer.data(pi_geom_1010_off + 229 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxyyy = cbuffer.data(pi_geom_1010_off + 230 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxyyz = cbuffer.data(pi_geom_1010_off + 231 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxyzz = cbuffer.data(pi_geom_1010_off + 232 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxxzzz = cbuffer.data(pi_geom_1010_off + 233 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxyyyy = cbuffer.data(pi_geom_1010_off + 234 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxyyyz = cbuffer.data(pi_geom_1010_off + 235 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxyyzz = cbuffer.data(pi_geom_1010_off + 236 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxyzzz = cbuffer.data(pi_geom_1010_off + 237 * ccomps * dcomps);

            auto g_x_0_z_0_z_xxzzzz = cbuffer.data(pi_geom_1010_off + 238 * ccomps * dcomps);

            auto g_x_0_z_0_z_xyyyyy = cbuffer.data(pi_geom_1010_off + 239 * ccomps * dcomps);

            auto g_x_0_z_0_z_xyyyyz = cbuffer.data(pi_geom_1010_off + 240 * ccomps * dcomps);

            auto g_x_0_z_0_z_xyyyzz = cbuffer.data(pi_geom_1010_off + 241 * ccomps * dcomps);

            auto g_x_0_z_0_z_xyyzzz = cbuffer.data(pi_geom_1010_off + 242 * ccomps * dcomps);

            auto g_x_0_z_0_z_xyzzzz = cbuffer.data(pi_geom_1010_off + 243 * ccomps * dcomps);

            auto g_x_0_z_0_z_xzzzzz = cbuffer.data(pi_geom_1010_off + 244 * ccomps * dcomps);

            auto g_x_0_z_0_z_yyyyyy = cbuffer.data(pi_geom_1010_off + 245 * ccomps * dcomps);

            auto g_x_0_z_0_z_yyyyyz = cbuffer.data(pi_geom_1010_off + 246 * ccomps * dcomps);

            auto g_x_0_z_0_z_yyyyzz = cbuffer.data(pi_geom_1010_off + 247 * ccomps * dcomps);

            auto g_x_0_z_0_z_yyyzzz = cbuffer.data(pi_geom_1010_off + 248 * ccomps * dcomps);

            auto g_x_0_z_0_z_yyzzzz = cbuffer.data(pi_geom_1010_off + 249 * ccomps * dcomps);

            auto g_x_0_z_0_z_yzzzzz = cbuffer.data(pi_geom_1010_off + 250 * ccomps * dcomps);

            auto g_x_0_z_0_z_zzzzzz = cbuffer.data(pi_geom_1010_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_0_0_xxxxxx, g_x_0_z_0_0_xxxxxxz, g_x_0_z_0_0_xxxxxy, g_x_0_z_0_0_xxxxxyz, g_x_0_z_0_0_xxxxxz, g_x_0_z_0_0_xxxxxzz, g_x_0_z_0_0_xxxxyy, g_x_0_z_0_0_xxxxyyz, g_x_0_z_0_0_xxxxyz, g_x_0_z_0_0_xxxxyzz, g_x_0_z_0_0_xxxxzz, g_x_0_z_0_0_xxxxzzz, g_x_0_z_0_0_xxxyyy, g_x_0_z_0_0_xxxyyyz, g_x_0_z_0_0_xxxyyz, g_x_0_z_0_0_xxxyyzz, g_x_0_z_0_0_xxxyzz, g_x_0_z_0_0_xxxyzzz, g_x_0_z_0_0_xxxzzz, g_x_0_z_0_0_xxxzzzz, g_x_0_z_0_0_xxyyyy, g_x_0_z_0_0_xxyyyyz, g_x_0_z_0_0_xxyyyz, g_x_0_z_0_0_xxyyyzz, g_x_0_z_0_0_xxyyzz, g_x_0_z_0_0_xxyyzzz, g_x_0_z_0_0_xxyzzz, g_x_0_z_0_0_xxyzzzz, g_x_0_z_0_0_xxzzzz, g_x_0_z_0_0_xxzzzzz, g_x_0_z_0_0_xyyyyy, g_x_0_z_0_0_xyyyyyz, g_x_0_z_0_0_xyyyyz, g_x_0_z_0_0_xyyyyzz, g_x_0_z_0_0_xyyyzz, g_x_0_z_0_0_xyyyzzz, g_x_0_z_0_0_xyyzzz, g_x_0_z_0_0_xyyzzzz, g_x_0_z_0_0_xyzzzz, g_x_0_z_0_0_xyzzzzz, g_x_0_z_0_0_xzzzzz, g_x_0_z_0_0_xzzzzzz, g_x_0_z_0_0_yyyyyy, g_x_0_z_0_0_yyyyyyz, g_x_0_z_0_0_yyyyyz, g_x_0_z_0_0_yyyyyzz, g_x_0_z_0_0_yyyyzz, g_x_0_z_0_0_yyyyzzz, g_x_0_z_0_0_yyyzzz, g_x_0_z_0_0_yyyzzzz, g_x_0_z_0_0_yyzzzz, g_x_0_z_0_0_yyzzzzz, g_x_0_z_0_0_yzzzzz, g_x_0_z_0_0_yzzzzzz, g_x_0_z_0_0_zzzzzz, g_x_0_z_0_0_zzzzzzz, g_x_0_z_0_z_xxxxxx, g_x_0_z_0_z_xxxxxy, g_x_0_z_0_z_xxxxxz, g_x_0_z_0_z_xxxxyy, g_x_0_z_0_z_xxxxyz, g_x_0_z_0_z_xxxxzz, g_x_0_z_0_z_xxxyyy, g_x_0_z_0_z_xxxyyz, g_x_0_z_0_z_xxxyzz, g_x_0_z_0_z_xxxzzz, g_x_0_z_0_z_xxyyyy, g_x_0_z_0_z_xxyyyz, g_x_0_z_0_z_xxyyzz, g_x_0_z_0_z_xxyzzz, g_x_0_z_0_z_xxzzzz, g_x_0_z_0_z_xyyyyy, g_x_0_z_0_z_xyyyyz, g_x_0_z_0_z_xyyyzz, g_x_0_z_0_z_xyyzzz, g_x_0_z_0_z_xyzzzz, g_x_0_z_0_z_xzzzzz, g_x_0_z_0_z_yyyyyy, g_x_0_z_0_z_yyyyyz, g_x_0_z_0_z_yyyyzz, g_x_0_z_0_z_yyyzzz, g_x_0_z_0_z_yyzzzz, g_x_0_z_0_z_yzzzzz, g_x_0_z_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_z_0_z_xxxxxx[k] = -g_x_0_z_0_0_xxxxxx[k] * ab_z + g_x_0_z_0_0_xxxxxxz[k];

                g_x_0_z_0_z_xxxxxy[k] = -g_x_0_z_0_0_xxxxxy[k] * ab_z + g_x_0_z_0_0_xxxxxyz[k];

                g_x_0_z_0_z_xxxxxz[k] = -g_x_0_z_0_0_xxxxxz[k] * ab_z + g_x_0_z_0_0_xxxxxzz[k];

                g_x_0_z_0_z_xxxxyy[k] = -g_x_0_z_0_0_xxxxyy[k] * ab_z + g_x_0_z_0_0_xxxxyyz[k];

                g_x_0_z_0_z_xxxxyz[k] = -g_x_0_z_0_0_xxxxyz[k] * ab_z + g_x_0_z_0_0_xxxxyzz[k];

                g_x_0_z_0_z_xxxxzz[k] = -g_x_0_z_0_0_xxxxzz[k] * ab_z + g_x_0_z_0_0_xxxxzzz[k];

                g_x_0_z_0_z_xxxyyy[k] = -g_x_0_z_0_0_xxxyyy[k] * ab_z + g_x_0_z_0_0_xxxyyyz[k];

                g_x_0_z_0_z_xxxyyz[k] = -g_x_0_z_0_0_xxxyyz[k] * ab_z + g_x_0_z_0_0_xxxyyzz[k];

                g_x_0_z_0_z_xxxyzz[k] = -g_x_0_z_0_0_xxxyzz[k] * ab_z + g_x_0_z_0_0_xxxyzzz[k];

                g_x_0_z_0_z_xxxzzz[k] = -g_x_0_z_0_0_xxxzzz[k] * ab_z + g_x_0_z_0_0_xxxzzzz[k];

                g_x_0_z_0_z_xxyyyy[k] = -g_x_0_z_0_0_xxyyyy[k] * ab_z + g_x_0_z_0_0_xxyyyyz[k];

                g_x_0_z_0_z_xxyyyz[k] = -g_x_0_z_0_0_xxyyyz[k] * ab_z + g_x_0_z_0_0_xxyyyzz[k];

                g_x_0_z_0_z_xxyyzz[k] = -g_x_0_z_0_0_xxyyzz[k] * ab_z + g_x_0_z_0_0_xxyyzzz[k];

                g_x_0_z_0_z_xxyzzz[k] = -g_x_0_z_0_0_xxyzzz[k] * ab_z + g_x_0_z_0_0_xxyzzzz[k];

                g_x_0_z_0_z_xxzzzz[k] = -g_x_0_z_0_0_xxzzzz[k] * ab_z + g_x_0_z_0_0_xxzzzzz[k];

                g_x_0_z_0_z_xyyyyy[k] = -g_x_0_z_0_0_xyyyyy[k] * ab_z + g_x_0_z_0_0_xyyyyyz[k];

                g_x_0_z_0_z_xyyyyz[k] = -g_x_0_z_0_0_xyyyyz[k] * ab_z + g_x_0_z_0_0_xyyyyzz[k];

                g_x_0_z_0_z_xyyyzz[k] = -g_x_0_z_0_0_xyyyzz[k] * ab_z + g_x_0_z_0_0_xyyyzzz[k];

                g_x_0_z_0_z_xyyzzz[k] = -g_x_0_z_0_0_xyyzzz[k] * ab_z + g_x_0_z_0_0_xyyzzzz[k];

                g_x_0_z_0_z_xyzzzz[k] = -g_x_0_z_0_0_xyzzzz[k] * ab_z + g_x_0_z_0_0_xyzzzzz[k];

                g_x_0_z_0_z_xzzzzz[k] = -g_x_0_z_0_0_xzzzzz[k] * ab_z + g_x_0_z_0_0_xzzzzzz[k];

                g_x_0_z_0_z_yyyyyy[k] = -g_x_0_z_0_0_yyyyyy[k] * ab_z + g_x_0_z_0_0_yyyyyyz[k];

                g_x_0_z_0_z_yyyyyz[k] = -g_x_0_z_0_0_yyyyyz[k] * ab_z + g_x_0_z_0_0_yyyyyzz[k];

                g_x_0_z_0_z_yyyyzz[k] = -g_x_0_z_0_0_yyyyzz[k] * ab_z + g_x_0_z_0_0_yyyyzzz[k];

                g_x_0_z_0_z_yyyzzz[k] = -g_x_0_z_0_0_yyyzzz[k] * ab_z + g_x_0_z_0_0_yyyzzzz[k];

                g_x_0_z_0_z_yyzzzz[k] = -g_x_0_z_0_0_yyzzzz[k] * ab_z + g_x_0_z_0_0_yyzzzzz[k];

                g_x_0_z_0_z_yzzzzz[k] = -g_x_0_z_0_0_yzzzzz[k] * ab_z + g_x_0_z_0_0_yzzzzzz[k];

                g_x_0_z_0_z_zzzzzz[k] = -g_x_0_z_0_0_zzzzzz[k] * ab_z + g_x_0_z_0_0_zzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_x_xxxxxx = cbuffer.data(pi_geom_1010_off + 252 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxxxy = cbuffer.data(pi_geom_1010_off + 253 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxxxz = cbuffer.data(pi_geom_1010_off + 254 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxxyy = cbuffer.data(pi_geom_1010_off + 255 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxxyz = cbuffer.data(pi_geom_1010_off + 256 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxxzz = cbuffer.data(pi_geom_1010_off + 257 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxyyy = cbuffer.data(pi_geom_1010_off + 258 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxyyz = cbuffer.data(pi_geom_1010_off + 259 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxyzz = cbuffer.data(pi_geom_1010_off + 260 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxxzzz = cbuffer.data(pi_geom_1010_off + 261 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxyyyy = cbuffer.data(pi_geom_1010_off + 262 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxyyyz = cbuffer.data(pi_geom_1010_off + 263 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxyyzz = cbuffer.data(pi_geom_1010_off + 264 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxyzzz = cbuffer.data(pi_geom_1010_off + 265 * ccomps * dcomps);

            auto g_y_0_x_0_x_xxzzzz = cbuffer.data(pi_geom_1010_off + 266 * ccomps * dcomps);

            auto g_y_0_x_0_x_xyyyyy = cbuffer.data(pi_geom_1010_off + 267 * ccomps * dcomps);

            auto g_y_0_x_0_x_xyyyyz = cbuffer.data(pi_geom_1010_off + 268 * ccomps * dcomps);

            auto g_y_0_x_0_x_xyyyzz = cbuffer.data(pi_geom_1010_off + 269 * ccomps * dcomps);

            auto g_y_0_x_0_x_xyyzzz = cbuffer.data(pi_geom_1010_off + 270 * ccomps * dcomps);

            auto g_y_0_x_0_x_xyzzzz = cbuffer.data(pi_geom_1010_off + 271 * ccomps * dcomps);

            auto g_y_0_x_0_x_xzzzzz = cbuffer.data(pi_geom_1010_off + 272 * ccomps * dcomps);

            auto g_y_0_x_0_x_yyyyyy = cbuffer.data(pi_geom_1010_off + 273 * ccomps * dcomps);

            auto g_y_0_x_0_x_yyyyyz = cbuffer.data(pi_geom_1010_off + 274 * ccomps * dcomps);

            auto g_y_0_x_0_x_yyyyzz = cbuffer.data(pi_geom_1010_off + 275 * ccomps * dcomps);

            auto g_y_0_x_0_x_yyyzzz = cbuffer.data(pi_geom_1010_off + 276 * ccomps * dcomps);

            auto g_y_0_x_0_x_yyzzzz = cbuffer.data(pi_geom_1010_off + 277 * ccomps * dcomps);

            auto g_y_0_x_0_x_yzzzzz = cbuffer.data(pi_geom_1010_off + 278 * ccomps * dcomps);

            auto g_y_0_x_0_x_zzzzzz = cbuffer.data(pi_geom_1010_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_0_xxxxxx, g_y_0_x_0_0_xxxxxxx, g_y_0_x_0_0_xxxxxxy, g_y_0_x_0_0_xxxxxxz, g_y_0_x_0_0_xxxxxy, g_y_0_x_0_0_xxxxxyy, g_y_0_x_0_0_xxxxxyz, g_y_0_x_0_0_xxxxxz, g_y_0_x_0_0_xxxxxzz, g_y_0_x_0_0_xxxxyy, g_y_0_x_0_0_xxxxyyy, g_y_0_x_0_0_xxxxyyz, g_y_0_x_0_0_xxxxyz, g_y_0_x_0_0_xxxxyzz, g_y_0_x_0_0_xxxxzz, g_y_0_x_0_0_xxxxzzz, g_y_0_x_0_0_xxxyyy, g_y_0_x_0_0_xxxyyyy, g_y_0_x_0_0_xxxyyyz, g_y_0_x_0_0_xxxyyz, g_y_0_x_0_0_xxxyyzz, g_y_0_x_0_0_xxxyzz, g_y_0_x_0_0_xxxyzzz, g_y_0_x_0_0_xxxzzz, g_y_0_x_0_0_xxxzzzz, g_y_0_x_0_0_xxyyyy, g_y_0_x_0_0_xxyyyyy, g_y_0_x_0_0_xxyyyyz, g_y_0_x_0_0_xxyyyz, g_y_0_x_0_0_xxyyyzz, g_y_0_x_0_0_xxyyzz, g_y_0_x_0_0_xxyyzzz, g_y_0_x_0_0_xxyzzz, g_y_0_x_0_0_xxyzzzz, g_y_0_x_0_0_xxzzzz, g_y_0_x_0_0_xxzzzzz, g_y_0_x_0_0_xyyyyy, g_y_0_x_0_0_xyyyyyy, g_y_0_x_0_0_xyyyyyz, g_y_0_x_0_0_xyyyyz, g_y_0_x_0_0_xyyyyzz, g_y_0_x_0_0_xyyyzz, g_y_0_x_0_0_xyyyzzz, g_y_0_x_0_0_xyyzzz, g_y_0_x_0_0_xyyzzzz, g_y_0_x_0_0_xyzzzz, g_y_0_x_0_0_xyzzzzz, g_y_0_x_0_0_xzzzzz, g_y_0_x_0_0_xzzzzzz, g_y_0_x_0_0_yyyyyy, g_y_0_x_0_0_yyyyyz, g_y_0_x_0_0_yyyyzz, g_y_0_x_0_0_yyyzzz, g_y_0_x_0_0_yyzzzz, g_y_0_x_0_0_yzzzzz, g_y_0_x_0_0_zzzzzz, g_y_0_x_0_x_xxxxxx, g_y_0_x_0_x_xxxxxy, g_y_0_x_0_x_xxxxxz, g_y_0_x_0_x_xxxxyy, g_y_0_x_0_x_xxxxyz, g_y_0_x_0_x_xxxxzz, g_y_0_x_0_x_xxxyyy, g_y_0_x_0_x_xxxyyz, g_y_0_x_0_x_xxxyzz, g_y_0_x_0_x_xxxzzz, g_y_0_x_0_x_xxyyyy, g_y_0_x_0_x_xxyyyz, g_y_0_x_0_x_xxyyzz, g_y_0_x_0_x_xxyzzz, g_y_0_x_0_x_xxzzzz, g_y_0_x_0_x_xyyyyy, g_y_0_x_0_x_xyyyyz, g_y_0_x_0_x_xyyyzz, g_y_0_x_0_x_xyyzzz, g_y_0_x_0_x_xyzzzz, g_y_0_x_0_x_xzzzzz, g_y_0_x_0_x_yyyyyy, g_y_0_x_0_x_yyyyyz, g_y_0_x_0_x_yyyyzz, g_y_0_x_0_x_yyyzzz, g_y_0_x_0_x_yyzzzz, g_y_0_x_0_x_yzzzzz, g_y_0_x_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_x_xxxxxx[k] = -g_y_0_x_0_0_xxxxxx[k] * ab_x + g_y_0_x_0_0_xxxxxxx[k];

                g_y_0_x_0_x_xxxxxy[k] = -g_y_0_x_0_0_xxxxxy[k] * ab_x + g_y_0_x_0_0_xxxxxxy[k];

                g_y_0_x_0_x_xxxxxz[k] = -g_y_0_x_0_0_xxxxxz[k] * ab_x + g_y_0_x_0_0_xxxxxxz[k];

                g_y_0_x_0_x_xxxxyy[k] = -g_y_0_x_0_0_xxxxyy[k] * ab_x + g_y_0_x_0_0_xxxxxyy[k];

                g_y_0_x_0_x_xxxxyz[k] = -g_y_0_x_0_0_xxxxyz[k] * ab_x + g_y_0_x_0_0_xxxxxyz[k];

                g_y_0_x_0_x_xxxxzz[k] = -g_y_0_x_0_0_xxxxzz[k] * ab_x + g_y_0_x_0_0_xxxxxzz[k];

                g_y_0_x_0_x_xxxyyy[k] = -g_y_0_x_0_0_xxxyyy[k] * ab_x + g_y_0_x_0_0_xxxxyyy[k];

                g_y_0_x_0_x_xxxyyz[k] = -g_y_0_x_0_0_xxxyyz[k] * ab_x + g_y_0_x_0_0_xxxxyyz[k];

                g_y_0_x_0_x_xxxyzz[k] = -g_y_0_x_0_0_xxxyzz[k] * ab_x + g_y_0_x_0_0_xxxxyzz[k];

                g_y_0_x_0_x_xxxzzz[k] = -g_y_0_x_0_0_xxxzzz[k] * ab_x + g_y_0_x_0_0_xxxxzzz[k];

                g_y_0_x_0_x_xxyyyy[k] = -g_y_0_x_0_0_xxyyyy[k] * ab_x + g_y_0_x_0_0_xxxyyyy[k];

                g_y_0_x_0_x_xxyyyz[k] = -g_y_0_x_0_0_xxyyyz[k] * ab_x + g_y_0_x_0_0_xxxyyyz[k];

                g_y_0_x_0_x_xxyyzz[k] = -g_y_0_x_0_0_xxyyzz[k] * ab_x + g_y_0_x_0_0_xxxyyzz[k];

                g_y_0_x_0_x_xxyzzz[k] = -g_y_0_x_0_0_xxyzzz[k] * ab_x + g_y_0_x_0_0_xxxyzzz[k];

                g_y_0_x_0_x_xxzzzz[k] = -g_y_0_x_0_0_xxzzzz[k] * ab_x + g_y_0_x_0_0_xxxzzzz[k];

                g_y_0_x_0_x_xyyyyy[k] = -g_y_0_x_0_0_xyyyyy[k] * ab_x + g_y_0_x_0_0_xxyyyyy[k];

                g_y_0_x_0_x_xyyyyz[k] = -g_y_0_x_0_0_xyyyyz[k] * ab_x + g_y_0_x_0_0_xxyyyyz[k];

                g_y_0_x_0_x_xyyyzz[k] = -g_y_0_x_0_0_xyyyzz[k] * ab_x + g_y_0_x_0_0_xxyyyzz[k];

                g_y_0_x_0_x_xyyzzz[k] = -g_y_0_x_0_0_xyyzzz[k] * ab_x + g_y_0_x_0_0_xxyyzzz[k];

                g_y_0_x_0_x_xyzzzz[k] = -g_y_0_x_0_0_xyzzzz[k] * ab_x + g_y_0_x_0_0_xxyzzzz[k];

                g_y_0_x_0_x_xzzzzz[k] = -g_y_0_x_0_0_xzzzzz[k] * ab_x + g_y_0_x_0_0_xxzzzzz[k];

                g_y_0_x_0_x_yyyyyy[k] = -g_y_0_x_0_0_yyyyyy[k] * ab_x + g_y_0_x_0_0_xyyyyyy[k];

                g_y_0_x_0_x_yyyyyz[k] = -g_y_0_x_0_0_yyyyyz[k] * ab_x + g_y_0_x_0_0_xyyyyyz[k];

                g_y_0_x_0_x_yyyyzz[k] = -g_y_0_x_0_0_yyyyzz[k] * ab_x + g_y_0_x_0_0_xyyyyzz[k];

                g_y_0_x_0_x_yyyzzz[k] = -g_y_0_x_0_0_yyyzzz[k] * ab_x + g_y_0_x_0_0_xyyyzzz[k];

                g_y_0_x_0_x_yyzzzz[k] = -g_y_0_x_0_0_yyzzzz[k] * ab_x + g_y_0_x_0_0_xyyzzzz[k];

                g_y_0_x_0_x_yzzzzz[k] = -g_y_0_x_0_0_yzzzzz[k] * ab_x + g_y_0_x_0_0_xyzzzzz[k];

                g_y_0_x_0_x_zzzzzz[k] = -g_y_0_x_0_0_zzzzzz[k] * ab_x + g_y_0_x_0_0_xzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_y_xxxxxx = cbuffer.data(pi_geom_1010_off + 280 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxxxy = cbuffer.data(pi_geom_1010_off + 281 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxxxz = cbuffer.data(pi_geom_1010_off + 282 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxxyy = cbuffer.data(pi_geom_1010_off + 283 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxxyz = cbuffer.data(pi_geom_1010_off + 284 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxxzz = cbuffer.data(pi_geom_1010_off + 285 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxyyy = cbuffer.data(pi_geom_1010_off + 286 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxyyz = cbuffer.data(pi_geom_1010_off + 287 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxyzz = cbuffer.data(pi_geom_1010_off + 288 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxxzzz = cbuffer.data(pi_geom_1010_off + 289 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxyyyy = cbuffer.data(pi_geom_1010_off + 290 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxyyyz = cbuffer.data(pi_geom_1010_off + 291 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxyyzz = cbuffer.data(pi_geom_1010_off + 292 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxyzzz = cbuffer.data(pi_geom_1010_off + 293 * ccomps * dcomps);

            auto g_y_0_x_0_y_xxzzzz = cbuffer.data(pi_geom_1010_off + 294 * ccomps * dcomps);

            auto g_y_0_x_0_y_xyyyyy = cbuffer.data(pi_geom_1010_off + 295 * ccomps * dcomps);

            auto g_y_0_x_0_y_xyyyyz = cbuffer.data(pi_geom_1010_off + 296 * ccomps * dcomps);

            auto g_y_0_x_0_y_xyyyzz = cbuffer.data(pi_geom_1010_off + 297 * ccomps * dcomps);

            auto g_y_0_x_0_y_xyyzzz = cbuffer.data(pi_geom_1010_off + 298 * ccomps * dcomps);

            auto g_y_0_x_0_y_xyzzzz = cbuffer.data(pi_geom_1010_off + 299 * ccomps * dcomps);

            auto g_y_0_x_0_y_xzzzzz = cbuffer.data(pi_geom_1010_off + 300 * ccomps * dcomps);

            auto g_y_0_x_0_y_yyyyyy = cbuffer.data(pi_geom_1010_off + 301 * ccomps * dcomps);

            auto g_y_0_x_0_y_yyyyyz = cbuffer.data(pi_geom_1010_off + 302 * ccomps * dcomps);

            auto g_y_0_x_0_y_yyyyzz = cbuffer.data(pi_geom_1010_off + 303 * ccomps * dcomps);

            auto g_y_0_x_0_y_yyyzzz = cbuffer.data(pi_geom_1010_off + 304 * ccomps * dcomps);

            auto g_y_0_x_0_y_yyzzzz = cbuffer.data(pi_geom_1010_off + 305 * ccomps * dcomps);

            auto g_y_0_x_0_y_yzzzzz = cbuffer.data(pi_geom_1010_off + 306 * ccomps * dcomps);

            auto g_y_0_x_0_y_zzzzzz = cbuffer.data(pi_geom_1010_off + 307 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_xxxxxx, g_0_0_x_0_0_xxxxxy, g_0_0_x_0_0_xxxxxz, g_0_0_x_0_0_xxxxyy, g_0_0_x_0_0_xxxxyz, g_0_0_x_0_0_xxxxzz, g_0_0_x_0_0_xxxyyy, g_0_0_x_0_0_xxxyyz, g_0_0_x_0_0_xxxyzz, g_0_0_x_0_0_xxxzzz, g_0_0_x_0_0_xxyyyy, g_0_0_x_0_0_xxyyyz, g_0_0_x_0_0_xxyyzz, g_0_0_x_0_0_xxyzzz, g_0_0_x_0_0_xxzzzz, g_0_0_x_0_0_xyyyyy, g_0_0_x_0_0_xyyyyz, g_0_0_x_0_0_xyyyzz, g_0_0_x_0_0_xyyzzz, g_0_0_x_0_0_xyzzzz, g_0_0_x_0_0_xzzzzz, g_0_0_x_0_0_yyyyyy, g_0_0_x_0_0_yyyyyz, g_0_0_x_0_0_yyyyzz, g_0_0_x_0_0_yyyzzz, g_0_0_x_0_0_yyzzzz, g_0_0_x_0_0_yzzzzz, g_0_0_x_0_0_zzzzzz, g_y_0_x_0_0_xxxxxx, g_y_0_x_0_0_xxxxxxy, g_y_0_x_0_0_xxxxxy, g_y_0_x_0_0_xxxxxyy, g_y_0_x_0_0_xxxxxyz, g_y_0_x_0_0_xxxxxz, g_y_0_x_0_0_xxxxyy, g_y_0_x_0_0_xxxxyyy, g_y_0_x_0_0_xxxxyyz, g_y_0_x_0_0_xxxxyz, g_y_0_x_0_0_xxxxyzz, g_y_0_x_0_0_xxxxzz, g_y_0_x_0_0_xxxyyy, g_y_0_x_0_0_xxxyyyy, g_y_0_x_0_0_xxxyyyz, g_y_0_x_0_0_xxxyyz, g_y_0_x_0_0_xxxyyzz, g_y_0_x_0_0_xxxyzz, g_y_0_x_0_0_xxxyzzz, g_y_0_x_0_0_xxxzzz, g_y_0_x_0_0_xxyyyy, g_y_0_x_0_0_xxyyyyy, g_y_0_x_0_0_xxyyyyz, g_y_0_x_0_0_xxyyyz, g_y_0_x_0_0_xxyyyzz, g_y_0_x_0_0_xxyyzz, g_y_0_x_0_0_xxyyzzz, g_y_0_x_0_0_xxyzzz, g_y_0_x_0_0_xxyzzzz, g_y_0_x_0_0_xxzzzz, g_y_0_x_0_0_xyyyyy, g_y_0_x_0_0_xyyyyyy, g_y_0_x_0_0_xyyyyyz, g_y_0_x_0_0_xyyyyz, g_y_0_x_0_0_xyyyyzz, g_y_0_x_0_0_xyyyzz, g_y_0_x_0_0_xyyyzzz, g_y_0_x_0_0_xyyzzz, g_y_0_x_0_0_xyyzzzz, g_y_0_x_0_0_xyzzzz, g_y_0_x_0_0_xyzzzzz, g_y_0_x_0_0_xzzzzz, g_y_0_x_0_0_yyyyyy, g_y_0_x_0_0_yyyyyyy, g_y_0_x_0_0_yyyyyyz, g_y_0_x_0_0_yyyyyz, g_y_0_x_0_0_yyyyyzz, g_y_0_x_0_0_yyyyzz, g_y_0_x_0_0_yyyyzzz, g_y_0_x_0_0_yyyzzz, g_y_0_x_0_0_yyyzzzz, g_y_0_x_0_0_yyzzzz, g_y_0_x_0_0_yyzzzzz, g_y_0_x_0_0_yzzzzz, g_y_0_x_0_0_yzzzzzz, g_y_0_x_0_0_zzzzzz, g_y_0_x_0_y_xxxxxx, g_y_0_x_0_y_xxxxxy, g_y_0_x_0_y_xxxxxz, g_y_0_x_0_y_xxxxyy, g_y_0_x_0_y_xxxxyz, g_y_0_x_0_y_xxxxzz, g_y_0_x_0_y_xxxyyy, g_y_0_x_0_y_xxxyyz, g_y_0_x_0_y_xxxyzz, g_y_0_x_0_y_xxxzzz, g_y_0_x_0_y_xxyyyy, g_y_0_x_0_y_xxyyyz, g_y_0_x_0_y_xxyyzz, g_y_0_x_0_y_xxyzzz, g_y_0_x_0_y_xxzzzz, g_y_0_x_0_y_xyyyyy, g_y_0_x_0_y_xyyyyz, g_y_0_x_0_y_xyyyzz, g_y_0_x_0_y_xyyzzz, g_y_0_x_0_y_xyzzzz, g_y_0_x_0_y_xzzzzz, g_y_0_x_0_y_yyyyyy, g_y_0_x_0_y_yyyyyz, g_y_0_x_0_y_yyyyzz, g_y_0_x_0_y_yyyzzz, g_y_0_x_0_y_yyzzzz, g_y_0_x_0_y_yzzzzz, g_y_0_x_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_y_xxxxxx[k] = -g_0_0_x_0_0_xxxxxx[k] - g_y_0_x_0_0_xxxxxx[k] * ab_y + g_y_0_x_0_0_xxxxxxy[k];

                g_y_0_x_0_y_xxxxxy[k] = -g_0_0_x_0_0_xxxxxy[k] - g_y_0_x_0_0_xxxxxy[k] * ab_y + g_y_0_x_0_0_xxxxxyy[k];

                g_y_0_x_0_y_xxxxxz[k] = -g_0_0_x_0_0_xxxxxz[k] - g_y_0_x_0_0_xxxxxz[k] * ab_y + g_y_0_x_0_0_xxxxxyz[k];

                g_y_0_x_0_y_xxxxyy[k] = -g_0_0_x_0_0_xxxxyy[k] - g_y_0_x_0_0_xxxxyy[k] * ab_y + g_y_0_x_0_0_xxxxyyy[k];

                g_y_0_x_0_y_xxxxyz[k] = -g_0_0_x_0_0_xxxxyz[k] - g_y_0_x_0_0_xxxxyz[k] * ab_y + g_y_0_x_0_0_xxxxyyz[k];

                g_y_0_x_0_y_xxxxzz[k] = -g_0_0_x_0_0_xxxxzz[k] - g_y_0_x_0_0_xxxxzz[k] * ab_y + g_y_0_x_0_0_xxxxyzz[k];

                g_y_0_x_0_y_xxxyyy[k] = -g_0_0_x_0_0_xxxyyy[k] - g_y_0_x_0_0_xxxyyy[k] * ab_y + g_y_0_x_0_0_xxxyyyy[k];

                g_y_0_x_0_y_xxxyyz[k] = -g_0_0_x_0_0_xxxyyz[k] - g_y_0_x_0_0_xxxyyz[k] * ab_y + g_y_0_x_0_0_xxxyyyz[k];

                g_y_0_x_0_y_xxxyzz[k] = -g_0_0_x_0_0_xxxyzz[k] - g_y_0_x_0_0_xxxyzz[k] * ab_y + g_y_0_x_0_0_xxxyyzz[k];

                g_y_0_x_0_y_xxxzzz[k] = -g_0_0_x_0_0_xxxzzz[k] - g_y_0_x_0_0_xxxzzz[k] * ab_y + g_y_0_x_0_0_xxxyzzz[k];

                g_y_0_x_0_y_xxyyyy[k] = -g_0_0_x_0_0_xxyyyy[k] - g_y_0_x_0_0_xxyyyy[k] * ab_y + g_y_0_x_0_0_xxyyyyy[k];

                g_y_0_x_0_y_xxyyyz[k] = -g_0_0_x_0_0_xxyyyz[k] - g_y_0_x_0_0_xxyyyz[k] * ab_y + g_y_0_x_0_0_xxyyyyz[k];

                g_y_0_x_0_y_xxyyzz[k] = -g_0_0_x_0_0_xxyyzz[k] - g_y_0_x_0_0_xxyyzz[k] * ab_y + g_y_0_x_0_0_xxyyyzz[k];

                g_y_0_x_0_y_xxyzzz[k] = -g_0_0_x_0_0_xxyzzz[k] - g_y_0_x_0_0_xxyzzz[k] * ab_y + g_y_0_x_0_0_xxyyzzz[k];

                g_y_0_x_0_y_xxzzzz[k] = -g_0_0_x_0_0_xxzzzz[k] - g_y_0_x_0_0_xxzzzz[k] * ab_y + g_y_0_x_0_0_xxyzzzz[k];

                g_y_0_x_0_y_xyyyyy[k] = -g_0_0_x_0_0_xyyyyy[k] - g_y_0_x_0_0_xyyyyy[k] * ab_y + g_y_0_x_0_0_xyyyyyy[k];

                g_y_0_x_0_y_xyyyyz[k] = -g_0_0_x_0_0_xyyyyz[k] - g_y_0_x_0_0_xyyyyz[k] * ab_y + g_y_0_x_0_0_xyyyyyz[k];

                g_y_0_x_0_y_xyyyzz[k] = -g_0_0_x_0_0_xyyyzz[k] - g_y_0_x_0_0_xyyyzz[k] * ab_y + g_y_0_x_0_0_xyyyyzz[k];

                g_y_0_x_0_y_xyyzzz[k] = -g_0_0_x_0_0_xyyzzz[k] - g_y_0_x_0_0_xyyzzz[k] * ab_y + g_y_0_x_0_0_xyyyzzz[k];

                g_y_0_x_0_y_xyzzzz[k] = -g_0_0_x_0_0_xyzzzz[k] - g_y_0_x_0_0_xyzzzz[k] * ab_y + g_y_0_x_0_0_xyyzzzz[k];

                g_y_0_x_0_y_xzzzzz[k] = -g_0_0_x_0_0_xzzzzz[k] - g_y_0_x_0_0_xzzzzz[k] * ab_y + g_y_0_x_0_0_xyzzzzz[k];

                g_y_0_x_0_y_yyyyyy[k] = -g_0_0_x_0_0_yyyyyy[k] - g_y_0_x_0_0_yyyyyy[k] * ab_y + g_y_0_x_0_0_yyyyyyy[k];

                g_y_0_x_0_y_yyyyyz[k] = -g_0_0_x_0_0_yyyyyz[k] - g_y_0_x_0_0_yyyyyz[k] * ab_y + g_y_0_x_0_0_yyyyyyz[k];

                g_y_0_x_0_y_yyyyzz[k] = -g_0_0_x_0_0_yyyyzz[k] - g_y_0_x_0_0_yyyyzz[k] * ab_y + g_y_0_x_0_0_yyyyyzz[k];

                g_y_0_x_0_y_yyyzzz[k] = -g_0_0_x_0_0_yyyzzz[k] - g_y_0_x_0_0_yyyzzz[k] * ab_y + g_y_0_x_0_0_yyyyzzz[k];

                g_y_0_x_0_y_yyzzzz[k] = -g_0_0_x_0_0_yyzzzz[k] - g_y_0_x_0_0_yyzzzz[k] * ab_y + g_y_0_x_0_0_yyyzzzz[k];

                g_y_0_x_0_y_yzzzzz[k] = -g_0_0_x_0_0_yzzzzz[k] - g_y_0_x_0_0_yzzzzz[k] * ab_y + g_y_0_x_0_0_yyzzzzz[k];

                g_y_0_x_0_y_zzzzzz[k] = -g_0_0_x_0_0_zzzzzz[k] - g_y_0_x_0_0_zzzzzz[k] * ab_y + g_y_0_x_0_0_yzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_y_0_x_0_z_xxxxxx = cbuffer.data(pi_geom_1010_off + 308 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxxxy = cbuffer.data(pi_geom_1010_off + 309 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxxxz = cbuffer.data(pi_geom_1010_off + 310 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxxyy = cbuffer.data(pi_geom_1010_off + 311 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxxyz = cbuffer.data(pi_geom_1010_off + 312 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxxzz = cbuffer.data(pi_geom_1010_off + 313 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxyyy = cbuffer.data(pi_geom_1010_off + 314 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxyyz = cbuffer.data(pi_geom_1010_off + 315 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxyzz = cbuffer.data(pi_geom_1010_off + 316 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxxzzz = cbuffer.data(pi_geom_1010_off + 317 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxyyyy = cbuffer.data(pi_geom_1010_off + 318 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxyyyz = cbuffer.data(pi_geom_1010_off + 319 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxyyzz = cbuffer.data(pi_geom_1010_off + 320 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxyzzz = cbuffer.data(pi_geom_1010_off + 321 * ccomps * dcomps);

            auto g_y_0_x_0_z_xxzzzz = cbuffer.data(pi_geom_1010_off + 322 * ccomps * dcomps);

            auto g_y_0_x_0_z_xyyyyy = cbuffer.data(pi_geom_1010_off + 323 * ccomps * dcomps);

            auto g_y_0_x_0_z_xyyyyz = cbuffer.data(pi_geom_1010_off + 324 * ccomps * dcomps);

            auto g_y_0_x_0_z_xyyyzz = cbuffer.data(pi_geom_1010_off + 325 * ccomps * dcomps);

            auto g_y_0_x_0_z_xyyzzz = cbuffer.data(pi_geom_1010_off + 326 * ccomps * dcomps);

            auto g_y_0_x_0_z_xyzzzz = cbuffer.data(pi_geom_1010_off + 327 * ccomps * dcomps);

            auto g_y_0_x_0_z_xzzzzz = cbuffer.data(pi_geom_1010_off + 328 * ccomps * dcomps);

            auto g_y_0_x_0_z_yyyyyy = cbuffer.data(pi_geom_1010_off + 329 * ccomps * dcomps);

            auto g_y_0_x_0_z_yyyyyz = cbuffer.data(pi_geom_1010_off + 330 * ccomps * dcomps);

            auto g_y_0_x_0_z_yyyyzz = cbuffer.data(pi_geom_1010_off + 331 * ccomps * dcomps);

            auto g_y_0_x_0_z_yyyzzz = cbuffer.data(pi_geom_1010_off + 332 * ccomps * dcomps);

            auto g_y_0_x_0_z_yyzzzz = cbuffer.data(pi_geom_1010_off + 333 * ccomps * dcomps);

            auto g_y_0_x_0_z_yzzzzz = cbuffer.data(pi_geom_1010_off + 334 * ccomps * dcomps);

            auto g_y_0_x_0_z_zzzzzz = cbuffer.data(pi_geom_1010_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_0_0_xxxxxx, g_y_0_x_0_0_xxxxxxz, g_y_0_x_0_0_xxxxxy, g_y_0_x_0_0_xxxxxyz, g_y_0_x_0_0_xxxxxz, g_y_0_x_0_0_xxxxxzz, g_y_0_x_0_0_xxxxyy, g_y_0_x_0_0_xxxxyyz, g_y_0_x_0_0_xxxxyz, g_y_0_x_0_0_xxxxyzz, g_y_0_x_0_0_xxxxzz, g_y_0_x_0_0_xxxxzzz, g_y_0_x_0_0_xxxyyy, g_y_0_x_0_0_xxxyyyz, g_y_0_x_0_0_xxxyyz, g_y_0_x_0_0_xxxyyzz, g_y_0_x_0_0_xxxyzz, g_y_0_x_0_0_xxxyzzz, g_y_0_x_0_0_xxxzzz, g_y_0_x_0_0_xxxzzzz, g_y_0_x_0_0_xxyyyy, g_y_0_x_0_0_xxyyyyz, g_y_0_x_0_0_xxyyyz, g_y_0_x_0_0_xxyyyzz, g_y_0_x_0_0_xxyyzz, g_y_0_x_0_0_xxyyzzz, g_y_0_x_0_0_xxyzzz, g_y_0_x_0_0_xxyzzzz, g_y_0_x_0_0_xxzzzz, g_y_0_x_0_0_xxzzzzz, g_y_0_x_0_0_xyyyyy, g_y_0_x_0_0_xyyyyyz, g_y_0_x_0_0_xyyyyz, g_y_0_x_0_0_xyyyyzz, g_y_0_x_0_0_xyyyzz, g_y_0_x_0_0_xyyyzzz, g_y_0_x_0_0_xyyzzz, g_y_0_x_0_0_xyyzzzz, g_y_0_x_0_0_xyzzzz, g_y_0_x_0_0_xyzzzzz, g_y_0_x_0_0_xzzzzz, g_y_0_x_0_0_xzzzzzz, g_y_0_x_0_0_yyyyyy, g_y_0_x_0_0_yyyyyyz, g_y_0_x_0_0_yyyyyz, g_y_0_x_0_0_yyyyyzz, g_y_0_x_0_0_yyyyzz, g_y_0_x_0_0_yyyyzzz, g_y_0_x_0_0_yyyzzz, g_y_0_x_0_0_yyyzzzz, g_y_0_x_0_0_yyzzzz, g_y_0_x_0_0_yyzzzzz, g_y_0_x_0_0_yzzzzz, g_y_0_x_0_0_yzzzzzz, g_y_0_x_0_0_zzzzzz, g_y_0_x_0_0_zzzzzzz, g_y_0_x_0_z_xxxxxx, g_y_0_x_0_z_xxxxxy, g_y_0_x_0_z_xxxxxz, g_y_0_x_0_z_xxxxyy, g_y_0_x_0_z_xxxxyz, g_y_0_x_0_z_xxxxzz, g_y_0_x_0_z_xxxyyy, g_y_0_x_0_z_xxxyyz, g_y_0_x_0_z_xxxyzz, g_y_0_x_0_z_xxxzzz, g_y_0_x_0_z_xxyyyy, g_y_0_x_0_z_xxyyyz, g_y_0_x_0_z_xxyyzz, g_y_0_x_0_z_xxyzzz, g_y_0_x_0_z_xxzzzz, g_y_0_x_0_z_xyyyyy, g_y_0_x_0_z_xyyyyz, g_y_0_x_0_z_xyyyzz, g_y_0_x_0_z_xyyzzz, g_y_0_x_0_z_xyzzzz, g_y_0_x_0_z_xzzzzz, g_y_0_x_0_z_yyyyyy, g_y_0_x_0_z_yyyyyz, g_y_0_x_0_z_yyyyzz, g_y_0_x_0_z_yyyzzz, g_y_0_x_0_z_yyzzzz, g_y_0_x_0_z_yzzzzz, g_y_0_x_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_x_0_z_xxxxxx[k] = -g_y_0_x_0_0_xxxxxx[k] * ab_z + g_y_0_x_0_0_xxxxxxz[k];

                g_y_0_x_0_z_xxxxxy[k] = -g_y_0_x_0_0_xxxxxy[k] * ab_z + g_y_0_x_0_0_xxxxxyz[k];

                g_y_0_x_0_z_xxxxxz[k] = -g_y_0_x_0_0_xxxxxz[k] * ab_z + g_y_0_x_0_0_xxxxxzz[k];

                g_y_0_x_0_z_xxxxyy[k] = -g_y_0_x_0_0_xxxxyy[k] * ab_z + g_y_0_x_0_0_xxxxyyz[k];

                g_y_0_x_0_z_xxxxyz[k] = -g_y_0_x_0_0_xxxxyz[k] * ab_z + g_y_0_x_0_0_xxxxyzz[k];

                g_y_0_x_0_z_xxxxzz[k] = -g_y_0_x_0_0_xxxxzz[k] * ab_z + g_y_0_x_0_0_xxxxzzz[k];

                g_y_0_x_0_z_xxxyyy[k] = -g_y_0_x_0_0_xxxyyy[k] * ab_z + g_y_0_x_0_0_xxxyyyz[k];

                g_y_0_x_0_z_xxxyyz[k] = -g_y_0_x_0_0_xxxyyz[k] * ab_z + g_y_0_x_0_0_xxxyyzz[k];

                g_y_0_x_0_z_xxxyzz[k] = -g_y_0_x_0_0_xxxyzz[k] * ab_z + g_y_0_x_0_0_xxxyzzz[k];

                g_y_0_x_0_z_xxxzzz[k] = -g_y_0_x_0_0_xxxzzz[k] * ab_z + g_y_0_x_0_0_xxxzzzz[k];

                g_y_0_x_0_z_xxyyyy[k] = -g_y_0_x_0_0_xxyyyy[k] * ab_z + g_y_0_x_0_0_xxyyyyz[k];

                g_y_0_x_0_z_xxyyyz[k] = -g_y_0_x_0_0_xxyyyz[k] * ab_z + g_y_0_x_0_0_xxyyyzz[k];

                g_y_0_x_0_z_xxyyzz[k] = -g_y_0_x_0_0_xxyyzz[k] * ab_z + g_y_0_x_0_0_xxyyzzz[k];

                g_y_0_x_0_z_xxyzzz[k] = -g_y_0_x_0_0_xxyzzz[k] * ab_z + g_y_0_x_0_0_xxyzzzz[k];

                g_y_0_x_0_z_xxzzzz[k] = -g_y_0_x_0_0_xxzzzz[k] * ab_z + g_y_0_x_0_0_xxzzzzz[k];

                g_y_0_x_0_z_xyyyyy[k] = -g_y_0_x_0_0_xyyyyy[k] * ab_z + g_y_0_x_0_0_xyyyyyz[k];

                g_y_0_x_0_z_xyyyyz[k] = -g_y_0_x_0_0_xyyyyz[k] * ab_z + g_y_0_x_0_0_xyyyyzz[k];

                g_y_0_x_0_z_xyyyzz[k] = -g_y_0_x_0_0_xyyyzz[k] * ab_z + g_y_0_x_0_0_xyyyzzz[k];

                g_y_0_x_0_z_xyyzzz[k] = -g_y_0_x_0_0_xyyzzz[k] * ab_z + g_y_0_x_0_0_xyyzzzz[k];

                g_y_0_x_0_z_xyzzzz[k] = -g_y_0_x_0_0_xyzzzz[k] * ab_z + g_y_0_x_0_0_xyzzzzz[k];

                g_y_0_x_0_z_xzzzzz[k] = -g_y_0_x_0_0_xzzzzz[k] * ab_z + g_y_0_x_0_0_xzzzzzz[k];

                g_y_0_x_0_z_yyyyyy[k] = -g_y_0_x_0_0_yyyyyy[k] * ab_z + g_y_0_x_0_0_yyyyyyz[k];

                g_y_0_x_0_z_yyyyyz[k] = -g_y_0_x_0_0_yyyyyz[k] * ab_z + g_y_0_x_0_0_yyyyyzz[k];

                g_y_0_x_0_z_yyyyzz[k] = -g_y_0_x_0_0_yyyyzz[k] * ab_z + g_y_0_x_0_0_yyyyzzz[k];

                g_y_0_x_0_z_yyyzzz[k] = -g_y_0_x_0_0_yyyzzz[k] * ab_z + g_y_0_x_0_0_yyyzzzz[k];

                g_y_0_x_0_z_yyzzzz[k] = -g_y_0_x_0_0_yyzzzz[k] * ab_z + g_y_0_x_0_0_yyzzzzz[k];

                g_y_0_x_0_z_yzzzzz[k] = -g_y_0_x_0_0_yzzzzz[k] * ab_z + g_y_0_x_0_0_yzzzzzz[k];

                g_y_0_x_0_z_zzzzzz[k] = -g_y_0_x_0_0_zzzzzz[k] * ab_z + g_y_0_x_0_0_zzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_x_xxxxxx = cbuffer.data(pi_geom_1010_off + 336 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxxxy = cbuffer.data(pi_geom_1010_off + 337 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxxxz = cbuffer.data(pi_geom_1010_off + 338 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxxyy = cbuffer.data(pi_geom_1010_off + 339 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxxyz = cbuffer.data(pi_geom_1010_off + 340 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxxzz = cbuffer.data(pi_geom_1010_off + 341 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxyyy = cbuffer.data(pi_geom_1010_off + 342 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxyyz = cbuffer.data(pi_geom_1010_off + 343 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxyzz = cbuffer.data(pi_geom_1010_off + 344 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxxzzz = cbuffer.data(pi_geom_1010_off + 345 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxyyyy = cbuffer.data(pi_geom_1010_off + 346 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxyyyz = cbuffer.data(pi_geom_1010_off + 347 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxyyzz = cbuffer.data(pi_geom_1010_off + 348 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxyzzz = cbuffer.data(pi_geom_1010_off + 349 * ccomps * dcomps);

            auto g_y_0_y_0_x_xxzzzz = cbuffer.data(pi_geom_1010_off + 350 * ccomps * dcomps);

            auto g_y_0_y_0_x_xyyyyy = cbuffer.data(pi_geom_1010_off + 351 * ccomps * dcomps);

            auto g_y_0_y_0_x_xyyyyz = cbuffer.data(pi_geom_1010_off + 352 * ccomps * dcomps);

            auto g_y_0_y_0_x_xyyyzz = cbuffer.data(pi_geom_1010_off + 353 * ccomps * dcomps);

            auto g_y_0_y_0_x_xyyzzz = cbuffer.data(pi_geom_1010_off + 354 * ccomps * dcomps);

            auto g_y_0_y_0_x_xyzzzz = cbuffer.data(pi_geom_1010_off + 355 * ccomps * dcomps);

            auto g_y_0_y_0_x_xzzzzz = cbuffer.data(pi_geom_1010_off + 356 * ccomps * dcomps);

            auto g_y_0_y_0_x_yyyyyy = cbuffer.data(pi_geom_1010_off + 357 * ccomps * dcomps);

            auto g_y_0_y_0_x_yyyyyz = cbuffer.data(pi_geom_1010_off + 358 * ccomps * dcomps);

            auto g_y_0_y_0_x_yyyyzz = cbuffer.data(pi_geom_1010_off + 359 * ccomps * dcomps);

            auto g_y_0_y_0_x_yyyzzz = cbuffer.data(pi_geom_1010_off + 360 * ccomps * dcomps);

            auto g_y_0_y_0_x_yyzzzz = cbuffer.data(pi_geom_1010_off + 361 * ccomps * dcomps);

            auto g_y_0_y_0_x_yzzzzz = cbuffer.data(pi_geom_1010_off + 362 * ccomps * dcomps);

            auto g_y_0_y_0_x_zzzzzz = cbuffer.data(pi_geom_1010_off + 363 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_0_xxxxxx, g_y_0_y_0_0_xxxxxxx, g_y_0_y_0_0_xxxxxxy, g_y_0_y_0_0_xxxxxxz, g_y_0_y_0_0_xxxxxy, g_y_0_y_0_0_xxxxxyy, g_y_0_y_0_0_xxxxxyz, g_y_0_y_0_0_xxxxxz, g_y_0_y_0_0_xxxxxzz, g_y_0_y_0_0_xxxxyy, g_y_0_y_0_0_xxxxyyy, g_y_0_y_0_0_xxxxyyz, g_y_0_y_0_0_xxxxyz, g_y_0_y_0_0_xxxxyzz, g_y_0_y_0_0_xxxxzz, g_y_0_y_0_0_xxxxzzz, g_y_0_y_0_0_xxxyyy, g_y_0_y_0_0_xxxyyyy, g_y_0_y_0_0_xxxyyyz, g_y_0_y_0_0_xxxyyz, g_y_0_y_0_0_xxxyyzz, g_y_0_y_0_0_xxxyzz, g_y_0_y_0_0_xxxyzzz, g_y_0_y_0_0_xxxzzz, g_y_0_y_0_0_xxxzzzz, g_y_0_y_0_0_xxyyyy, g_y_0_y_0_0_xxyyyyy, g_y_0_y_0_0_xxyyyyz, g_y_0_y_0_0_xxyyyz, g_y_0_y_0_0_xxyyyzz, g_y_0_y_0_0_xxyyzz, g_y_0_y_0_0_xxyyzzz, g_y_0_y_0_0_xxyzzz, g_y_0_y_0_0_xxyzzzz, g_y_0_y_0_0_xxzzzz, g_y_0_y_0_0_xxzzzzz, g_y_0_y_0_0_xyyyyy, g_y_0_y_0_0_xyyyyyy, g_y_0_y_0_0_xyyyyyz, g_y_0_y_0_0_xyyyyz, g_y_0_y_0_0_xyyyyzz, g_y_0_y_0_0_xyyyzz, g_y_0_y_0_0_xyyyzzz, g_y_0_y_0_0_xyyzzz, g_y_0_y_0_0_xyyzzzz, g_y_0_y_0_0_xyzzzz, g_y_0_y_0_0_xyzzzzz, g_y_0_y_0_0_xzzzzz, g_y_0_y_0_0_xzzzzzz, g_y_0_y_0_0_yyyyyy, g_y_0_y_0_0_yyyyyz, g_y_0_y_0_0_yyyyzz, g_y_0_y_0_0_yyyzzz, g_y_0_y_0_0_yyzzzz, g_y_0_y_0_0_yzzzzz, g_y_0_y_0_0_zzzzzz, g_y_0_y_0_x_xxxxxx, g_y_0_y_0_x_xxxxxy, g_y_0_y_0_x_xxxxxz, g_y_0_y_0_x_xxxxyy, g_y_0_y_0_x_xxxxyz, g_y_0_y_0_x_xxxxzz, g_y_0_y_0_x_xxxyyy, g_y_0_y_0_x_xxxyyz, g_y_0_y_0_x_xxxyzz, g_y_0_y_0_x_xxxzzz, g_y_0_y_0_x_xxyyyy, g_y_0_y_0_x_xxyyyz, g_y_0_y_0_x_xxyyzz, g_y_0_y_0_x_xxyzzz, g_y_0_y_0_x_xxzzzz, g_y_0_y_0_x_xyyyyy, g_y_0_y_0_x_xyyyyz, g_y_0_y_0_x_xyyyzz, g_y_0_y_0_x_xyyzzz, g_y_0_y_0_x_xyzzzz, g_y_0_y_0_x_xzzzzz, g_y_0_y_0_x_yyyyyy, g_y_0_y_0_x_yyyyyz, g_y_0_y_0_x_yyyyzz, g_y_0_y_0_x_yyyzzz, g_y_0_y_0_x_yyzzzz, g_y_0_y_0_x_yzzzzz, g_y_0_y_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_x_xxxxxx[k] = -g_y_0_y_0_0_xxxxxx[k] * ab_x + g_y_0_y_0_0_xxxxxxx[k];

                g_y_0_y_0_x_xxxxxy[k] = -g_y_0_y_0_0_xxxxxy[k] * ab_x + g_y_0_y_0_0_xxxxxxy[k];

                g_y_0_y_0_x_xxxxxz[k] = -g_y_0_y_0_0_xxxxxz[k] * ab_x + g_y_0_y_0_0_xxxxxxz[k];

                g_y_0_y_0_x_xxxxyy[k] = -g_y_0_y_0_0_xxxxyy[k] * ab_x + g_y_0_y_0_0_xxxxxyy[k];

                g_y_0_y_0_x_xxxxyz[k] = -g_y_0_y_0_0_xxxxyz[k] * ab_x + g_y_0_y_0_0_xxxxxyz[k];

                g_y_0_y_0_x_xxxxzz[k] = -g_y_0_y_0_0_xxxxzz[k] * ab_x + g_y_0_y_0_0_xxxxxzz[k];

                g_y_0_y_0_x_xxxyyy[k] = -g_y_0_y_0_0_xxxyyy[k] * ab_x + g_y_0_y_0_0_xxxxyyy[k];

                g_y_0_y_0_x_xxxyyz[k] = -g_y_0_y_0_0_xxxyyz[k] * ab_x + g_y_0_y_0_0_xxxxyyz[k];

                g_y_0_y_0_x_xxxyzz[k] = -g_y_0_y_0_0_xxxyzz[k] * ab_x + g_y_0_y_0_0_xxxxyzz[k];

                g_y_0_y_0_x_xxxzzz[k] = -g_y_0_y_0_0_xxxzzz[k] * ab_x + g_y_0_y_0_0_xxxxzzz[k];

                g_y_0_y_0_x_xxyyyy[k] = -g_y_0_y_0_0_xxyyyy[k] * ab_x + g_y_0_y_0_0_xxxyyyy[k];

                g_y_0_y_0_x_xxyyyz[k] = -g_y_0_y_0_0_xxyyyz[k] * ab_x + g_y_0_y_0_0_xxxyyyz[k];

                g_y_0_y_0_x_xxyyzz[k] = -g_y_0_y_0_0_xxyyzz[k] * ab_x + g_y_0_y_0_0_xxxyyzz[k];

                g_y_0_y_0_x_xxyzzz[k] = -g_y_0_y_0_0_xxyzzz[k] * ab_x + g_y_0_y_0_0_xxxyzzz[k];

                g_y_0_y_0_x_xxzzzz[k] = -g_y_0_y_0_0_xxzzzz[k] * ab_x + g_y_0_y_0_0_xxxzzzz[k];

                g_y_0_y_0_x_xyyyyy[k] = -g_y_0_y_0_0_xyyyyy[k] * ab_x + g_y_0_y_0_0_xxyyyyy[k];

                g_y_0_y_0_x_xyyyyz[k] = -g_y_0_y_0_0_xyyyyz[k] * ab_x + g_y_0_y_0_0_xxyyyyz[k];

                g_y_0_y_0_x_xyyyzz[k] = -g_y_0_y_0_0_xyyyzz[k] * ab_x + g_y_0_y_0_0_xxyyyzz[k];

                g_y_0_y_0_x_xyyzzz[k] = -g_y_0_y_0_0_xyyzzz[k] * ab_x + g_y_0_y_0_0_xxyyzzz[k];

                g_y_0_y_0_x_xyzzzz[k] = -g_y_0_y_0_0_xyzzzz[k] * ab_x + g_y_0_y_0_0_xxyzzzz[k];

                g_y_0_y_0_x_xzzzzz[k] = -g_y_0_y_0_0_xzzzzz[k] * ab_x + g_y_0_y_0_0_xxzzzzz[k];

                g_y_0_y_0_x_yyyyyy[k] = -g_y_0_y_0_0_yyyyyy[k] * ab_x + g_y_0_y_0_0_xyyyyyy[k];

                g_y_0_y_0_x_yyyyyz[k] = -g_y_0_y_0_0_yyyyyz[k] * ab_x + g_y_0_y_0_0_xyyyyyz[k];

                g_y_0_y_0_x_yyyyzz[k] = -g_y_0_y_0_0_yyyyzz[k] * ab_x + g_y_0_y_0_0_xyyyyzz[k];

                g_y_0_y_0_x_yyyzzz[k] = -g_y_0_y_0_0_yyyzzz[k] * ab_x + g_y_0_y_0_0_xyyyzzz[k];

                g_y_0_y_0_x_yyzzzz[k] = -g_y_0_y_0_0_yyzzzz[k] * ab_x + g_y_0_y_0_0_xyyzzzz[k];

                g_y_0_y_0_x_yzzzzz[k] = -g_y_0_y_0_0_yzzzzz[k] * ab_x + g_y_0_y_0_0_xyzzzzz[k];

                g_y_0_y_0_x_zzzzzz[k] = -g_y_0_y_0_0_zzzzzz[k] * ab_x + g_y_0_y_0_0_xzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_y_xxxxxx = cbuffer.data(pi_geom_1010_off + 364 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxxxy = cbuffer.data(pi_geom_1010_off + 365 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxxxz = cbuffer.data(pi_geom_1010_off + 366 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxxyy = cbuffer.data(pi_geom_1010_off + 367 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxxyz = cbuffer.data(pi_geom_1010_off + 368 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxxzz = cbuffer.data(pi_geom_1010_off + 369 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxyyy = cbuffer.data(pi_geom_1010_off + 370 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxyyz = cbuffer.data(pi_geom_1010_off + 371 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxyzz = cbuffer.data(pi_geom_1010_off + 372 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxxzzz = cbuffer.data(pi_geom_1010_off + 373 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxyyyy = cbuffer.data(pi_geom_1010_off + 374 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxyyyz = cbuffer.data(pi_geom_1010_off + 375 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxyyzz = cbuffer.data(pi_geom_1010_off + 376 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxyzzz = cbuffer.data(pi_geom_1010_off + 377 * ccomps * dcomps);

            auto g_y_0_y_0_y_xxzzzz = cbuffer.data(pi_geom_1010_off + 378 * ccomps * dcomps);

            auto g_y_0_y_0_y_xyyyyy = cbuffer.data(pi_geom_1010_off + 379 * ccomps * dcomps);

            auto g_y_0_y_0_y_xyyyyz = cbuffer.data(pi_geom_1010_off + 380 * ccomps * dcomps);

            auto g_y_0_y_0_y_xyyyzz = cbuffer.data(pi_geom_1010_off + 381 * ccomps * dcomps);

            auto g_y_0_y_0_y_xyyzzz = cbuffer.data(pi_geom_1010_off + 382 * ccomps * dcomps);

            auto g_y_0_y_0_y_xyzzzz = cbuffer.data(pi_geom_1010_off + 383 * ccomps * dcomps);

            auto g_y_0_y_0_y_xzzzzz = cbuffer.data(pi_geom_1010_off + 384 * ccomps * dcomps);

            auto g_y_0_y_0_y_yyyyyy = cbuffer.data(pi_geom_1010_off + 385 * ccomps * dcomps);

            auto g_y_0_y_0_y_yyyyyz = cbuffer.data(pi_geom_1010_off + 386 * ccomps * dcomps);

            auto g_y_0_y_0_y_yyyyzz = cbuffer.data(pi_geom_1010_off + 387 * ccomps * dcomps);

            auto g_y_0_y_0_y_yyyzzz = cbuffer.data(pi_geom_1010_off + 388 * ccomps * dcomps);

            auto g_y_0_y_0_y_yyzzzz = cbuffer.data(pi_geom_1010_off + 389 * ccomps * dcomps);

            auto g_y_0_y_0_y_yzzzzz = cbuffer.data(pi_geom_1010_off + 390 * ccomps * dcomps);

            auto g_y_0_y_0_y_zzzzzz = cbuffer.data(pi_geom_1010_off + 391 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_xxxxxx, g_0_0_y_0_0_xxxxxy, g_0_0_y_0_0_xxxxxz, g_0_0_y_0_0_xxxxyy, g_0_0_y_0_0_xxxxyz, g_0_0_y_0_0_xxxxzz, g_0_0_y_0_0_xxxyyy, g_0_0_y_0_0_xxxyyz, g_0_0_y_0_0_xxxyzz, g_0_0_y_0_0_xxxzzz, g_0_0_y_0_0_xxyyyy, g_0_0_y_0_0_xxyyyz, g_0_0_y_0_0_xxyyzz, g_0_0_y_0_0_xxyzzz, g_0_0_y_0_0_xxzzzz, g_0_0_y_0_0_xyyyyy, g_0_0_y_0_0_xyyyyz, g_0_0_y_0_0_xyyyzz, g_0_0_y_0_0_xyyzzz, g_0_0_y_0_0_xyzzzz, g_0_0_y_0_0_xzzzzz, g_0_0_y_0_0_yyyyyy, g_0_0_y_0_0_yyyyyz, g_0_0_y_0_0_yyyyzz, g_0_0_y_0_0_yyyzzz, g_0_0_y_0_0_yyzzzz, g_0_0_y_0_0_yzzzzz, g_0_0_y_0_0_zzzzzz, g_y_0_y_0_0_xxxxxx, g_y_0_y_0_0_xxxxxxy, g_y_0_y_0_0_xxxxxy, g_y_0_y_0_0_xxxxxyy, g_y_0_y_0_0_xxxxxyz, g_y_0_y_0_0_xxxxxz, g_y_0_y_0_0_xxxxyy, g_y_0_y_0_0_xxxxyyy, g_y_0_y_0_0_xxxxyyz, g_y_0_y_0_0_xxxxyz, g_y_0_y_0_0_xxxxyzz, g_y_0_y_0_0_xxxxzz, g_y_0_y_0_0_xxxyyy, g_y_0_y_0_0_xxxyyyy, g_y_0_y_0_0_xxxyyyz, g_y_0_y_0_0_xxxyyz, g_y_0_y_0_0_xxxyyzz, g_y_0_y_0_0_xxxyzz, g_y_0_y_0_0_xxxyzzz, g_y_0_y_0_0_xxxzzz, g_y_0_y_0_0_xxyyyy, g_y_0_y_0_0_xxyyyyy, g_y_0_y_0_0_xxyyyyz, g_y_0_y_0_0_xxyyyz, g_y_0_y_0_0_xxyyyzz, g_y_0_y_0_0_xxyyzz, g_y_0_y_0_0_xxyyzzz, g_y_0_y_0_0_xxyzzz, g_y_0_y_0_0_xxyzzzz, g_y_0_y_0_0_xxzzzz, g_y_0_y_0_0_xyyyyy, g_y_0_y_0_0_xyyyyyy, g_y_0_y_0_0_xyyyyyz, g_y_0_y_0_0_xyyyyz, g_y_0_y_0_0_xyyyyzz, g_y_0_y_0_0_xyyyzz, g_y_0_y_0_0_xyyyzzz, g_y_0_y_0_0_xyyzzz, g_y_0_y_0_0_xyyzzzz, g_y_0_y_0_0_xyzzzz, g_y_0_y_0_0_xyzzzzz, g_y_0_y_0_0_xzzzzz, g_y_0_y_0_0_yyyyyy, g_y_0_y_0_0_yyyyyyy, g_y_0_y_0_0_yyyyyyz, g_y_0_y_0_0_yyyyyz, g_y_0_y_0_0_yyyyyzz, g_y_0_y_0_0_yyyyzz, g_y_0_y_0_0_yyyyzzz, g_y_0_y_0_0_yyyzzz, g_y_0_y_0_0_yyyzzzz, g_y_0_y_0_0_yyzzzz, g_y_0_y_0_0_yyzzzzz, g_y_0_y_0_0_yzzzzz, g_y_0_y_0_0_yzzzzzz, g_y_0_y_0_0_zzzzzz, g_y_0_y_0_y_xxxxxx, g_y_0_y_0_y_xxxxxy, g_y_0_y_0_y_xxxxxz, g_y_0_y_0_y_xxxxyy, g_y_0_y_0_y_xxxxyz, g_y_0_y_0_y_xxxxzz, g_y_0_y_0_y_xxxyyy, g_y_0_y_0_y_xxxyyz, g_y_0_y_0_y_xxxyzz, g_y_0_y_0_y_xxxzzz, g_y_0_y_0_y_xxyyyy, g_y_0_y_0_y_xxyyyz, g_y_0_y_0_y_xxyyzz, g_y_0_y_0_y_xxyzzz, g_y_0_y_0_y_xxzzzz, g_y_0_y_0_y_xyyyyy, g_y_0_y_0_y_xyyyyz, g_y_0_y_0_y_xyyyzz, g_y_0_y_0_y_xyyzzz, g_y_0_y_0_y_xyzzzz, g_y_0_y_0_y_xzzzzz, g_y_0_y_0_y_yyyyyy, g_y_0_y_0_y_yyyyyz, g_y_0_y_0_y_yyyyzz, g_y_0_y_0_y_yyyzzz, g_y_0_y_0_y_yyzzzz, g_y_0_y_0_y_yzzzzz, g_y_0_y_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_y_xxxxxx[k] = -g_0_0_y_0_0_xxxxxx[k] - g_y_0_y_0_0_xxxxxx[k] * ab_y + g_y_0_y_0_0_xxxxxxy[k];

                g_y_0_y_0_y_xxxxxy[k] = -g_0_0_y_0_0_xxxxxy[k] - g_y_0_y_0_0_xxxxxy[k] * ab_y + g_y_0_y_0_0_xxxxxyy[k];

                g_y_0_y_0_y_xxxxxz[k] = -g_0_0_y_0_0_xxxxxz[k] - g_y_0_y_0_0_xxxxxz[k] * ab_y + g_y_0_y_0_0_xxxxxyz[k];

                g_y_0_y_0_y_xxxxyy[k] = -g_0_0_y_0_0_xxxxyy[k] - g_y_0_y_0_0_xxxxyy[k] * ab_y + g_y_0_y_0_0_xxxxyyy[k];

                g_y_0_y_0_y_xxxxyz[k] = -g_0_0_y_0_0_xxxxyz[k] - g_y_0_y_0_0_xxxxyz[k] * ab_y + g_y_0_y_0_0_xxxxyyz[k];

                g_y_0_y_0_y_xxxxzz[k] = -g_0_0_y_0_0_xxxxzz[k] - g_y_0_y_0_0_xxxxzz[k] * ab_y + g_y_0_y_0_0_xxxxyzz[k];

                g_y_0_y_0_y_xxxyyy[k] = -g_0_0_y_0_0_xxxyyy[k] - g_y_0_y_0_0_xxxyyy[k] * ab_y + g_y_0_y_0_0_xxxyyyy[k];

                g_y_0_y_0_y_xxxyyz[k] = -g_0_0_y_0_0_xxxyyz[k] - g_y_0_y_0_0_xxxyyz[k] * ab_y + g_y_0_y_0_0_xxxyyyz[k];

                g_y_0_y_0_y_xxxyzz[k] = -g_0_0_y_0_0_xxxyzz[k] - g_y_0_y_0_0_xxxyzz[k] * ab_y + g_y_0_y_0_0_xxxyyzz[k];

                g_y_0_y_0_y_xxxzzz[k] = -g_0_0_y_0_0_xxxzzz[k] - g_y_0_y_0_0_xxxzzz[k] * ab_y + g_y_0_y_0_0_xxxyzzz[k];

                g_y_0_y_0_y_xxyyyy[k] = -g_0_0_y_0_0_xxyyyy[k] - g_y_0_y_0_0_xxyyyy[k] * ab_y + g_y_0_y_0_0_xxyyyyy[k];

                g_y_0_y_0_y_xxyyyz[k] = -g_0_0_y_0_0_xxyyyz[k] - g_y_0_y_0_0_xxyyyz[k] * ab_y + g_y_0_y_0_0_xxyyyyz[k];

                g_y_0_y_0_y_xxyyzz[k] = -g_0_0_y_0_0_xxyyzz[k] - g_y_0_y_0_0_xxyyzz[k] * ab_y + g_y_0_y_0_0_xxyyyzz[k];

                g_y_0_y_0_y_xxyzzz[k] = -g_0_0_y_0_0_xxyzzz[k] - g_y_0_y_0_0_xxyzzz[k] * ab_y + g_y_0_y_0_0_xxyyzzz[k];

                g_y_0_y_0_y_xxzzzz[k] = -g_0_0_y_0_0_xxzzzz[k] - g_y_0_y_0_0_xxzzzz[k] * ab_y + g_y_0_y_0_0_xxyzzzz[k];

                g_y_0_y_0_y_xyyyyy[k] = -g_0_0_y_0_0_xyyyyy[k] - g_y_0_y_0_0_xyyyyy[k] * ab_y + g_y_0_y_0_0_xyyyyyy[k];

                g_y_0_y_0_y_xyyyyz[k] = -g_0_0_y_0_0_xyyyyz[k] - g_y_0_y_0_0_xyyyyz[k] * ab_y + g_y_0_y_0_0_xyyyyyz[k];

                g_y_0_y_0_y_xyyyzz[k] = -g_0_0_y_0_0_xyyyzz[k] - g_y_0_y_0_0_xyyyzz[k] * ab_y + g_y_0_y_0_0_xyyyyzz[k];

                g_y_0_y_0_y_xyyzzz[k] = -g_0_0_y_0_0_xyyzzz[k] - g_y_0_y_0_0_xyyzzz[k] * ab_y + g_y_0_y_0_0_xyyyzzz[k];

                g_y_0_y_0_y_xyzzzz[k] = -g_0_0_y_0_0_xyzzzz[k] - g_y_0_y_0_0_xyzzzz[k] * ab_y + g_y_0_y_0_0_xyyzzzz[k];

                g_y_0_y_0_y_xzzzzz[k] = -g_0_0_y_0_0_xzzzzz[k] - g_y_0_y_0_0_xzzzzz[k] * ab_y + g_y_0_y_0_0_xyzzzzz[k];

                g_y_0_y_0_y_yyyyyy[k] = -g_0_0_y_0_0_yyyyyy[k] - g_y_0_y_0_0_yyyyyy[k] * ab_y + g_y_0_y_0_0_yyyyyyy[k];

                g_y_0_y_0_y_yyyyyz[k] = -g_0_0_y_0_0_yyyyyz[k] - g_y_0_y_0_0_yyyyyz[k] * ab_y + g_y_0_y_0_0_yyyyyyz[k];

                g_y_0_y_0_y_yyyyzz[k] = -g_0_0_y_0_0_yyyyzz[k] - g_y_0_y_0_0_yyyyzz[k] * ab_y + g_y_0_y_0_0_yyyyyzz[k];

                g_y_0_y_0_y_yyyzzz[k] = -g_0_0_y_0_0_yyyzzz[k] - g_y_0_y_0_0_yyyzzz[k] * ab_y + g_y_0_y_0_0_yyyyzzz[k];

                g_y_0_y_0_y_yyzzzz[k] = -g_0_0_y_0_0_yyzzzz[k] - g_y_0_y_0_0_yyzzzz[k] * ab_y + g_y_0_y_0_0_yyyzzzz[k];

                g_y_0_y_0_y_yzzzzz[k] = -g_0_0_y_0_0_yzzzzz[k] - g_y_0_y_0_0_yzzzzz[k] * ab_y + g_y_0_y_0_0_yyzzzzz[k];

                g_y_0_y_0_y_zzzzzz[k] = -g_0_0_y_0_0_zzzzzz[k] - g_y_0_y_0_0_zzzzzz[k] * ab_y + g_y_0_y_0_0_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_y_0_z_xxxxxx = cbuffer.data(pi_geom_1010_off + 392 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxxxy = cbuffer.data(pi_geom_1010_off + 393 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxxxz = cbuffer.data(pi_geom_1010_off + 394 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxxyy = cbuffer.data(pi_geom_1010_off + 395 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxxyz = cbuffer.data(pi_geom_1010_off + 396 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxxzz = cbuffer.data(pi_geom_1010_off + 397 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxyyy = cbuffer.data(pi_geom_1010_off + 398 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxyyz = cbuffer.data(pi_geom_1010_off + 399 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxyzz = cbuffer.data(pi_geom_1010_off + 400 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxxzzz = cbuffer.data(pi_geom_1010_off + 401 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxyyyy = cbuffer.data(pi_geom_1010_off + 402 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxyyyz = cbuffer.data(pi_geom_1010_off + 403 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxyyzz = cbuffer.data(pi_geom_1010_off + 404 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxyzzz = cbuffer.data(pi_geom_1010_off + 405 * ccomps * dcomps);

            auto g_y_0_y_0_z_xxzzzz = cbuffer.data(pi_geom_1010_off + 406 * ccomps * dcomps);

            auto g_y_0_y_0_z_xyyyyy = cbuffer.data(pi_geom_1010_off + 407 * ccomps * dcomps);

            auto g_y_0_y_0_z_xyyyyz = cbuffer.data(pi_geom_1010_off + 408 * ccomps * dcomps);

            auto g_y_0_y_0_z_xyyyzz = cbuffer.data(pi_geom_1010_off + 409 * ccomps * dcomps);

            auto g_y_0_y_0_z_xyyzzz = cbuffer.data(pi_geom_1010_off + 410 * ccomps * dcomps);

            auto g_y_0_y_0_z_xyzzzz = cbuffer.data(pi_geom_1010_off + 411 * ccomps * dcomps);

            auto g_y_0_y_0_z_xzzzzz = cbuffer.data(pi_geom_1010_off + 412 * ccomps * dcomps);

            auto g_y_0_y_0_z_yyyyyy = cbuffer.data(pi_geom_1010_off + 413 * ccomps * dcomps);

            auto g_y_0_y_0_z_yyyyyz = cbuffer.data(pi_geom_1010_off + 414 * ccomps * dcomps);

            auto g_y_0_y_0_z_yyyyzz = cbuffer.data(pi_geom_1010_off + 415 * ccomps * dcomps);

            auto g_y_0_y_0_z_yyyzzz = cbuffer.data(pi_geom_1010_off + 416 * ccomps * dcomps);

            auto g_y_0_y_0_z_yyzzzz = cbuffer.data(pi_geom_1010_off + 417 * ccomps * dcomps);

            auto g_y_0_y_0_z_yzzzzz = cbuffer.data(pi_geom_1010_off + 418 * ccomps * dcomps);

            auto g_y_0_y_0_z_zzzzzz = cbuffer.data(pi_geom_1010_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_0_0_xxxxxx, g_y_0_y_0_0_xxxxxxz, g_y_0_y_0_0_xxxxxy, g_y_0_y_0_0_xxxxxyz, g_y_0_y_0_0_xxxxxz, g_y_0_y_0_0_xxxxxzz, g_y_0_y_0_0_xxxxyy, g_y_0_y_0_0_xxxxyyz, g_y_0_y_0_0_xxxxyz, g_y_0_y_0_0_xxxxyzz, g_y_0_y_0_0_xxxxzz, g_y_0_y_0_0_xxxxzzz, g_y_0_y_0_0_xxxyyy, g_y_0_y_0_0_xxxyyyz, g_y_0_y_0_0_xxxyyz, g_y_0_y_0_0_xxxyyzz, g_y_0_y_0_0_xxxyzz, g_y_0_y_0_0_xxxyzzz, g_y_0_y_0_0_xxxzzz, g_y_0_y_0_0_xxxzzzz, g_y_0_y_0_0_xxyyyy, g_y_0_y_0_0_xxyyyyz, g_y_0_y_0_0_xxyyyz, g_y_0_y_0_0_xxyyyzz, g_y_0_y_0_0_xxyyzz, g_y_0_y_0_0_xxyyzzz, g_y_0_y_0_0_xxyzzz, g_y_0_y_0_0_xxyzzzz, g_y_0_y_0_0_xxzzzz, g_y_0_y_0_0_xxzzzzz, g_y_0_y_0_0_xyyyyy, g_y_0_y_0_0_xyyyyyz, g_y_0_y_0_0_xyyyyz, g_y_0_y_0_0_xyyyyzz, g_y_0_y_0_0_xyyyzz, g_y_0_y_0_0_xyyyzzz, g_y_0_y_0_0_xyyzzz, g_y_0_y_0_0_xyyzzzz, g_y_0_y_0_0_xyzzzz, g_y_0_y_0_0_xyzzzzz, g_y_0_y_0_0_xzzzzz, g_y_0_y_0_0_xzzzzzz, g_y_0_y_0_0_yyyyyy, g_y_0_y_0_0_yyyyyyz, g_y_0_y_0_0_yyyyyz, g_y_0_y_0_0_yyyyyzz, g_y_0_y_0_0_yyyyzz, g_y_0_y_0_0_yyyyzzz, g_y_0_y_0_0_yyyzzz, g_y_0_y_0_0_yyyzzzz, g_y_0_y_0_0_yyzzzz, g_y_0_y_0_0_yyzzzzz, g_y_0_y_0_0_yzzzzz, g_y_0_y_0_0_yzzzzzz, g_y_0_y_0_0_zzzzzz, g_y_0_y_0_0_zzzzzzz, g_y_0_y_0_z_xxxxxx, g_y_0_y_0_z_xxxxxy, g_y_0_y_0_z_xxxxxz, g_y_0_y_0_z_xxxxyy, g_y_0_y_0_z_xxxxyz, g_y_0_y_0_z_xxxxzz, g_y_0_y_0_z_xxxyyy, g_y_0_y_0_z_xxxyyz, g_y_0_y_0_z_xxxyzz, g_y_0_y_0_z_xxxzzz, g_y_0_y_0_z_xxyyyy, g_y_0_y_0_z_xxyyyz, g_y_0_y_0_z_xxyyzz, g_y_0_y_0_z_xxyzzz, g_y_0_y_0_z_xxzzzz, g_y_0_y_0_z_xyyyyy, g_y_0_y_0_z_xyyyyz, g_y_0_y_0_z_xyyyzz, g_y_0_y_0_z_xyyzzz, g_y_0_y_0_z_xyzzzz, g_y_0_y_0_z_xzzzzz, g_y_0_y_0_z_yyyyyy, g_y_0_y_0_z_yyyyyz, g_y_0_y_0_z_yyyyzz, g_y_0_y_0_z_yyyzzz, g_y_0_y_0_z_yyzzzz, g_y_0_y_0_z_yzzzzz, g_y_0_y_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_y_0_z_xxxxxx[k] = -g_y_0_y_0_0_xxxxxx[k] * ab_z + g_y_0_y_0_0_xxxxxxz[k];

                g_y_0_y_0_z_xxxxxy[k] = -g_y_0_y_0_0_xxxxxy[k] * ab_z + g_y_0_y_0_0_xxxxxyz[k];

                g_y_0_y_0_z_xxxxxz[k] = -g_y_0_y_0_0_xxxxxz[k] * ab_z + g_y_0_y_0_0_xxxxxzz[k];

                g_y_0_y_0_z_xxxxyy[k] = -g_y_0_y_0_0_xxxxyy[k] * ab_z + g_y_0_y_0_0_xxxxyyz[k];

                g_y_0_y_0_z_xxxxyz[k] = -g_y_0_y_0_0_xxxxyz[k] * ab_z + g_y_0_y_0_0_xxxxyzz[k];

                g_y_0_y_0_z_xxxxzz[k] = -g_y_0_y_0_0_xxxxzz[k] * ab_z + g_y_0_y_0_0_xxxxzzz[k];

                g_y_0_y_0_z_xxxyyy[k] = -g_y_0_y_0_0_xxxyyy[k] * ab_z + g_y_0_y_0_0_xxxyyyz[k];

                g_y_0_y_0_z_xxxyyz[k] = -g_y_0_y_0_0_xxxyyz[k] * ab_z + g_y_0_y_0_0_xxxyyzz[k];

                g_y_0_y_0_z_xxxyzz[k] = -g_y_0_y_0_0_xxxyzz[k] * ab_z + g_y_0_y_0_0_xxxyzzz[k];

                g_y_0_y_0_z_xxxzzz[k] = -g_y_0_y_0_0_xxxzzz[k] * ab_z + g_y_0_y_0_0_xxxzzzz[k];

                g_y_0_y_0_z_xxyyyy[k] = -g_y_0_y_0_0_xxyyyy[k] * ab_z + g_y_0_y_0_0_xxyyyyz[k];

                g_y_0_y_0_z_xxyyyz[k] = -g_y_0_y_0_0_xxyyyz[k] * ab_z + g_y_0_y_0_0_xxyyyzz[k];

                g_y_0_y_0_z_xxyyzz[k] = -g_y_0_y_0_0_xxyyzz[k] * ab_z + g_y_0_y_0_0_xxyyzzz[k];

                g_y_0_y_0_z_xxyzzz[k] = -g_y_0_y_0_0_xxyzzz[k] * ab_z + g_y_0_y_0_0_xxyzzzz[k];

                g_y_0_y_0_z_xxzzzz[k] = -g_y_0_y_0_0_xxzzzz[k] * ab_z + g_y_0_y_0_0_xxzzzzz[k];

                g_y_0_y_0_z_xyyyyy[k] = -g_y_0_y_0_0_xyyyyy[k] * ab_z + g_y_0_y_0_0_xyyyyyz[k];

                g_y_0_y_0_z_xyyyyz[k] = -g_y_0_y_0_0_xyyyyz[k] * ab_z + g_y_0_y_0_0_xyyyyzz[k];

                g_y_0_y_0_z_xyyyzz[k] = -g_y_0_y_0_0_xyyyzz[k] * ab_z + g_y_0_y_0_0_xyyyzzz[k];

                g_y_0_y_0_z_xyyzzz[k] = -g_y_0_y_0_0_xyyzzz[k] * ab_z + g_y_0_y_0_0_xyyzzzz[k];

                g_y_0_y_0_z_xyzzzz[k] = -g_y_0_y_0_0_xyzzzz[k] * ab_z + g_y_0_y_0_0_xyzzzzz[k];

                g_y_0_y_0_z_xzzzzz[k] = -g_y_0_y_0_0_xzzzzz[k] * ab_z + g_y_0_y_0_0_xzzzzzz[k];

                g_y_0_y_0_z_yyyyyy[k] = -g_y_0_y_0_0_yyyyyy[k] * ab_z + g_y_0_y_0_0_yyyyyyz[k];

                g_y_0_y_0_z_yyyyyz[k] = -g_y_0_y_0_0_yyyyyz[k] * ab_z + g_y_0_y_0_0_yyyyyzz[k];

                g_y_0_y_0_z_yyyyzz[k] = -g_y_0_y_0_0_yyyyzz[k] * ab_z + g_y_0_y_0_0_yyyyzzz[k];

                g_y_0_y_0_z_yyyzzz[k] = -g_y_0_y_0_0_yyyzzz[k] * ab_z + g_y_0_y_0_0_yyyzzzz[k];

                g_y_0_y_0_z_yyzzzz[k] = -g_y_0_y_0_0_yyzzzz[k] * ab_z + g_y_0_y_0_0_yyzzzzz[k];

                g_y_0_y_0_z_yzzzzz[k] = -g_y_0_y_0_0_yzzzzz[k] * ab_z + g_y_0_y_0_0_yzzzzzz[k];

                g_y_0_y_0_z_zzzzzz[k] = -g_y_0_y_0_0_zzzzzz[k] * ab_z + g_y_0_y_0_0_zzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_x_xxxxxx = cbuffer.data(pi_geom_1010_off + 420 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxxxy = cbuffer.data(pi_geom_1010_off + 421 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxxxz = cbuffer.data(pi_geom_1010_off + 422 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxxyy = cbuffer.data(pi_geom_1010_off + 423 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxxyz = cbuffer.data(pi_geom_1010_off + 424 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxxzz = cbuffer.data(pi_geom_1010_off + 425 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxyyy = cbuffer.data(pi_geom_1010_off + 426 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxyyz = cbuffer.data(pi_geom_1010_off + 427 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxyzz = cbuffer.data(pi_geom_1010_off + 428 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxxzzz = cbuffer.data(pi_geom_1010_off + 429 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxyyyy = cbuffer.data(pi_geom_1010_off + 430 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxyyyz = cbuffer.data(pi_geom_1010_off + 431 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxyyzz = cbuffer.data(pi_geom_1010_off + 432 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxyzzz = cbuffer.data(pi_geom_1010_off + 433 * ccomps * dcomps);

            auto g_y_0_z_0_x_xxzzzz = cbuffer.data(pi_geom_1010_off + 434 * ccomps * dcomps);

            auto g_y_0_z_0_x_xyyyyy = cbuffer.data(pi_geom_1010_off + 435 * ccomps * dcomps);

            auto g_y_0_z_0_x_xyyyyz = cbuffer.data(pi_geom_1010_off + 436 * ccomps * dcomps);

            auto g_y_0_z_0_x_xyyyzz = cbuffer.data(pi_geom_1010_off + 437 * ccomps * dcomps);

            auto g_y_0_z_0_x_xyyzzz = cbuffer.data(pi_geom_1010_off + 438 * ccomps * dcomps);

            auto g_y_0_z_0_x_xyzzzz = cbuffer.data(pi_geom_1010_off + 439 * ccomps * dcomps);

            auto g_y_0_z_0_x_xzzzzz = cbuffer.data(pi_geom_1010_off + 440 * ccomps * dcomps);

            auto g_y_0_z_0_x_yyyyyy = cbuffer.data(pi_geom_1010_off + 441 * ccomps * dcomps);

            auto g_y_0_z_0_x_yyyyyz = cbuffer.data(pi_geom_1010_off + 442 * ccomps * dcomps);

            auto g_y_0_z_0_x_yyyyzz = cbuffer.data(pi_geom_1010_off + 443 * ccomps * dcomps);

            auto g_y_0_z_0_x_yyyzzz = cbuffer.data(pi_geom_1010_off + 444 * ccomps * dcomps);

            auto g_y_0_z_0_x_yyzzzz = cbuffer.data(pi_geom_1010_off + 445 * ccomps * dcomps);

            auto g_y_0_z_0_x_yzzzzz = cbuffer.data(pi_geom_1010_off + 446 * ccomps * dcomps);

            auto g_y_0_z_0_x_zzzzzz = cbuffer.data(pi_geom_1010_off + 447 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_0_xxxxxx, g_y_0_z_0_0_xxxxxxx, g_y_0_z_0_0_xxxxxxy, g_y_0_z_0_0_xxxxxxz, g_y_0_z_0_0_xxxxxy, g_y_0_z_0_0_xxxxxyy, g_y_0_z_0_0_xxxxxyz, g_y_0_z_0_0_xxxxxz, g_y_0_z_0_0_xxxxxzz, g_y_0_z_0_0_xxxxyy, g_y_0_z_0_0_xxxxyyy, g_y_0_z_0_0_xxxxyyz, g_y_0_z_0_0_xxxxyz, g_y_0_z_0_0_xxxxyzz, g_y_0_z_0_0_xxxxzz, g_y_0_z_0_0_xxxxzzz, g_y_0_z_0_0_xxxyyy, g_y_0_z_0_0_xxxyyyy, g_y_0_z_0_0_xxxyyyz, g_y_0_z_0_0_xxxyyz, g_y_0_z_0_0_xxxyyzz, g_y_0_z_0_0_xxxyzz, g_y_0_z_0_0_xxxyzzz, g_y_0_z_0_0_xxxzzz, g_y_0_z_0_0_xxxzzzz, g_y_0_z_0_0_xxyyyy, g_y_0_z_0_0_xxyyyyy, g_y_0_z_0_0_xxyyyyz, g_y_0_z_0_0_xxyyyz, g_y_0_z_0_0_xxyyyzz, g_y_0_z_0_0_xxyyzz, g_y_0_z_0_0_xxyyzzz, g_y_0_z_0_0_xxyzzz, g_y_0_z_0_0_xxyzzzz, g_y_0_z_0_0_xxzzzz, g_y_0_z_0_0_xxzzzzz, g_y_0_z_0_0_xyyyyy, g_y_0_z_0_0_xyyyyyy, g_y_0_z_0_0_xyyyyyz, g_y_0_z_0_0_xyyyyz, g_y_0_z_0_0_xyyyyzz, g_y_0_z_0_0_xyyyzz, g_y_0_z_0_0_xyyyzzz, g_y_0_z_0_0_xyyzzz, g_y_0_z_0_0_xyyzzzz, g_y_0_z_0_0_xyzzzz, g_y_0_z_0_0_xyzzzzz, g_y_0_z_0_0_xzzzzz, g_y_0_z_0_0_xzzzzzz, g_y_0_z_0_0_yyyyyy, g_y_0_z_0_0_yyyyyz, g_y_0_z_0_0_yyyyzz, g_y_0_z_0_0_yyyzzz, g_y_0_z_0_0_yyzzzz, g_y_0_z_0_0_yzzzzz, g_y_0_z_0_0_zzzzzz, g_y_0_z_0_x_xxxxxx, g_y_0_z_0_x_xxxxxy, g_y_0_z_0_x_xxxxxz, g_y_0_z_0_x_xxxxyy, g_y_0_z_0_x_xxxxyz, g_y_0_z_0_x_xxxxzz, g_y_0_z_0_x_xxxyyy, g_y_0_z_0_x_xxxyyz, g_y_0_z_0_x_xxxyzz, g_y_0_z_0_x_xxxzzz, g_y_0_z_0_x_xxyyyy, g_y_0_z_0_x_xxyyyz, g_y_0_z_0_x_xxyyzz, g_y_0_z_0_x_xxyzzz, g_y_0_z_0_x_xxzzzz, g_y_0_z_0_x_xyyyyy, g_y_0_z_0_x_xyyyyz, g_y_0_z_0_x_xyyyzz, g_y_0_z_0_x_xyyzzz, g_y_0_z_0_x_xyzzzz, g_y_0_z_0_x_xzzzzz, g_y_0_z_0_x_yyyyyy, g_y_0_z_0_x_yyyyyz, g_y_0_z_0_x_yyyyzz, g_y_0_z_0_x_yyyzzz, g_y_0_z_0_x_yyzzzz, g_y_0_z_0_x_yzzzzz, g_y_0_z_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_x_xxxxxx[k] = -g_y_0_z_0_0_xxxxxx[k] * ab_x + g_y_0_z_0_0_xxxxxxx[k];

                g_y_0_z_0_x_xxxxxy[k] = -g_y_0_z_0_0_xxxxxy[k] * ab_x + g_y_0_z_0_0_xxxxxxy[k];

                g_y_0_z_0_x_xxxxxz[k] = -g_y_0_z_0_0_xxxxxz[k] * ab_x + g_y_0_z_0_0_xxxxxxz[k];

                g_y_0_z_0_x_xxxxyy[k] = -g_y_0_z_0_0_xxxxyy[k] * ab_x + g_y_0_z_0_0_xxxxxyy[k];

                g_y_0_z_0_x_xxxxyz[k] = -g_y_0_z_0_0_xxxxyz[k] * ab_x + g_y_0_z_0_0_xxxxxyz[k];

                g_y_0_z_0_x_xxxxzz[k] = -g_y_0_z_0_0_xxxxzz[k] * ab_x + g_y_0_z_0_0_xxxxxzz[k];

                g_y_0_z_0_x_xxxyyy[k] = -g_y_0_z_0_0_xxxyyy[k] * ab_x + g_y_0_z_0_0_xxxxyyy[k];

                g_y_0_z_0_x_xxxyyz[k] = -g_y_0_z_0_0_xxxyyz[k] * ab_x + g_y_0_z_0_0_xxxxyyz[k];

                g_y_0_z_0_x_xxxyzz[k] = -g_y_0_z_0_0_xxxyzz[k] * ab_x + g_y_0_z_0_0_xxxxyzz[k];

                g_y_0_z_0_x_xxxzzz[k] = -g_y_0_z_0_0_xxxzzz[k] * ab_x + g_y_0_z_0_0_xxxxzzz[k];

                g_y_0_z_0_x_xxyyyy[k] = -g_y_0_z_0_0_xxyyyy[k] * ab_x + g_y_0_z_0_0_xxxyyyy[k];

                g_y_0_z_0_x_xxyyyz[k] = -g_y_0_z_0_0_xxyyyz[k] * ab_x + g_y_0_z_0_0_xxxyyyz[k];

                g_y_0_z_0_x_xxyyzz[k] = -g_y_0_z_0_0_xxyyzz[k] * ab_x + g_y_0_z_0_0_xxxyyzz[k];

                g_y_0_z_0_x_xxyzzz[k] = -g_y_0_z_0_0_xxyzzz[k] * ab_x + g_y_0_z_0_0_xxxyzzz[k];

                g_y_0_z_0_x_xxzzzz[k] = -g_y_0_z_0_0_xxzzzz[k] * ab_x + g_y_0_z_0_0_xxxzzzz[k];

                g_y_0_z_0_x_xyyyyy[k] = -g_y_0_z_0_0_xyyyyy[k] * ab_x + g_y_0_z_0_0_xxyyyyy[k];

                g_y_0_z_0_x_xyyyyz[k] = -g_y_0_z_0_0_xyyyyz[k] * ab_x + g_y_0_z_0_0_xxyyyyz[k];

                g_y_0_z_0_x_xyyyzz[k] = -g_y_0_z_0_0_xyyyzz[k] * ab_x + g_y_0_z_0_0_xxyyyzz[k];

                g_y_0_z_0_x_xyyzzz[k] = -g_y_0_z_0_0_xyyzzz[k] * ab_x + g_y_0_z_0_0_xxyyzzz[k];

                g_y_0_z_0_x_xyzzzz[k] = -g_y_0_z_0_0_xyzzzz[k] * ab_x + g_y_0_z_0_0_xxyzzzz[k];

                g_y_0_z_0_x_xzzzzz[k] = -g_y_0_z_0_0_xzzzzz[k] * ab_x + g_y_0_z_0_0_xxzzzzz[k];

                g_y_0_z_0_x_yyyyyy[k] = -g_y_0_z_0_0_yyyyyy[k] * ab_x + g_y_0_z_0_0_xyyyyyy[k];

                g_y_0_z_0_x_yyyyyz[k] = -g_y_0_z_0_0_yyyyyz[k] * ab_x + g_y_0_z_0_0_xyyyyyz[k];

                g_y_0_z_0_x_yyyyzz[k] = -g_y_0_z_0_0_yyyyzz[k] * ab_x + g_y_0_z_0_0_xyyyyzz[k];

                g_y_0_z_0_x_yyyzzz[k] = -g_y_0_z_0_0_yyyzzz[k] * ab_x + g_y_0_z_0_0_xyyyzzz[k];

                g_y_0_z_0_x_yyzzzz[k] = -g_y_0_z_0_0_yyzzzz[k] * ab_x + g_y_0_z_0_0_xyyzzzz[k];

                g_y_0_z_0_x_yzzzzz[k] = -g_y_0_z_0_0_yzzzzz[k] * ab_x + g_y_0_z_0_0_xyzzzzz[k];

                g_y_0_z_0_x_zzzzzz[k] = -g_y_0_z_0_0_zzzzzz[k] * ab_x + g_y_0_z_0_0_xzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_y_xxxxxx = cbuffer.data(pi_geom_1010_off + 448 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxxxy = cbuffer.data(pi_geom_1010_off + 449 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxxxz = cbuffer.data(pi_geom_1010_off + 450 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxxyy = cbuffer.data(pi_geom_1010_off + 451 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxxyz = cbuffer.data(pi_geom_1010_off + 452 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxxzz = cbuffer.data(pi_geom_1010_off + 453 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxyyy = cbuffer.data(pi_geom_1010_off + 454 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxyyz = cbuffer.data(pi_geom_1010_off + 455 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxyzz = cbuffer.data(pi_geom_1010_off + 456 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxxzzz = cbuffer.data(pi_geom_1010_off + 457 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxyyyy = cbuffer.data(pi_geom_1010_off + 458 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxyyyz = cbuffer.data(pi_geom_1010_off + 459 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxyyzz = cbuffer.data(pi_geom_1010_off + 460 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxyzzz = cbuffer.data(pi_geom_1010_off + 461 * ccomps * dcomps);

            auto g_y_0_z_0_y_xxzzzz = cbuffer.data(pi_geom_1010_off + 462 * ccomps * dcomps);

            auto g_y_0_z_0_y_xyyyyy = cbuffer.data(pi_geom_1010_off + 463 * ccomps * dcomps);

            auto g_y_0_z_0_y_xyyyyz = cbuffer.data(pi_geom_1010_off + 464 * ccomps * dcomps);

            auto g_y_0_z_0_y_xyyyzz = cbuffer.data(pi_geom_1010_off + 465 * ccomps * dcomps);

            auto g_y_0_z_0_y_xyyzzz = cbuffer.data(pi_geom_1010_off + 466 * ccomps * dcomps);

            auto g_y_0_z_0_y_xyzzzz = cbuffer.data(pi_geom_1010_off + 467 * ccomps * dcomps);

            auto g_y_0_z_0_y_xzzzzz = cbuffer.data(pi_geom_1010_off + 468 * ccomps * dcomps);

            auto g_y_0_z_0_y_yyyyyy = cbuffer.data(pi_geom_1010_off + 469 * ccomps * dcomps);

            auto g_y_0_z_0_y_yyyyyz = cbuffer.data(pi_geom_1010_off + 470 * ccomps * dcomps);

            auto g_y_0_z_0_y_yyyyzz = cbuffer.data(pi_geom_1010_off + 471 * ccomps * dcomps);

            auto g_y_0_z_0_y_yyyzzz = cbuffer.data(pi_geom_1010_off + 472 * ccomps * dcomps);

            auto g_y_0_z_0_y_yyzzzz = cbuffer.data(pi_geom_1010_off + 473 * ccomps * dcomps);

            auto g_y_0_z_0_y_yzzzzz = cbuffer.data(pi_geom_1010_off + 474 * ccomps * dcomps);

            auto g_y_0_z_0_y_zzzzzz = cbuffer.data(pi_geom_1010_off + 475 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_xxxxxx, g_0_0_z_0_0_xxxxxy, g_0_0_z_0_0_xxxxxz, g_0_0_z_0_0_xxxxyy, g_0_0_z_0_0_xxxxyz, g_0_0_z_0_0_xxxxzz, g_0_0_z_0_0_xxxyyy, g_0_0_z_0_0_xxxyyz, g_0_0_z_0_0_xxxyzz, g_0_0_z_0_0_xxxzzz, g_0_0_z_0_0_xxyyyy, g_0_0_z_0_0_xxyyyz, g_0_0_z_0_0_xxyyzz, g_0_0_z_0_0_xxyzzz, g_0_0_z_0_0_xxzzzz, g_0_0_z_0_0_xyyyyy, g_0_0_z_0_0_xyyyyz, g_0_0_z_0_0_xyyyzz, g_0_0_z_0_0_xyyzzz, g_0_0_z_0_0_xyzzzz, g_0_0_z_0_0_xzzzzz, g_0_0_z_0_0_yyyyyy, g_0_0_z_0_0_yyyyyz, g_0_0_z_0_0_yyyyzz, g_0_0_z_0_0_yyyzzz, g_0_0_z_0_0_yyzzzz, g_0_0_z_0_0_yzzzzz, g_0_0_z_0_0_zzzzzz, g_y_0_z_0_0_xxxxxx, g_y_0_z_0_0_xxxxxxy, g_y_0_z_0_0_xxxxxy, g_y_0_z_0_0_xxxxxyy, g_y_0_z_0_0_xxxxxyz, g_y_0_z_0_0_xxxxxz, g_y_0_z_0_0_xxxxyy, g_y_0_z_0_0_xxxxyyy, g_y_0_z_0_0_xxxxyyz, g_y_0_z_0_0_xxxxyz, g_y_0_z_0_0_xxxxyzz, g_y_0_z_0_0_xxxxzz, g_y_0_z_0_0_xxxyyy, g_y_0_z_0_0_xxxyyyy, g_y_0_z_0_0_xxxyyyz, g_y_0_z_0_0_xxxyyz, g_y_0_z_0_0_xxxyyzz, g_y_0_z_0_0_xxxyzz, g_y_0_z_0_0_xxxyzzz, g_y_0_z_0_0_xxxzzz, g_y_0_z_0_0_xxyyyy, g_y_0_z_0_0_xxyyyyy, g_y_0_z_0_0_xxyyyyz, g_y_0_z_0_0_xxyyyz, g_y_0_z_0_0_xxyyyzz, g_y_0_z_0_0_xxyyzz, g_y_0_z_0_0_xxyyzzz, g_y_0_z_0_0_xxyzzz, g_y_0_z_0_0_xxyzzzz, g_y_0_z_0_0_xxzzzz, g_y_0_z_0_0_xyyyyy, g_y_0_z_0_0_xyyyyyy, g_y_0_z_0_0_xyyyyyz, g_y_0_z_0_0_xyyyyz, g_y_0_z_0_0_xyyyyzz, g_y_0_z_0_0_xyyyzz, g_y_0_z_0_0_xyyyzzz, g_y_0_z_0_0_xyyzzz, g_y_0_z_0_0_xyyzzzz, g_y_0_z_0_0_xyzzzz, g_y_0_z_0_0_xyzzzzz, g_y_0_z_0_0_xzzzzz, g_y_0_z_0_0_yyyyyy, g_y_0_z_0_0_yyyyyyy, g_y_0_z_0_0_yyyyyyz, g_y_0_z_0_0_yyyyyz, g_y_0_z_0_0_yyyyyzz, g_y_0_z_0_0_yyyyzz, g_y_0_z_0_0_yyyyzzz, g_y_0_z_0_0_yyyzzz, g_y_0_z_0_0_yyyzzzz, g_y_0_z_0_0_yyzzzz, g_y_0_z_0_0_yyzzzzz, g_y_0_z_0_0_yzzzzz, g_y_0_z_0_0_yzzzzzz, g_y_0_z_0_0_zzzzzz, g_y_0_z_0_y_xxxxxx, g_y_0_z_0_y_xxxxxy, g_y_0_z_0_y_xxxxxz, g_y_0_z_0_y_xxxxyy, g_y_0_z_0_y_xxxxyz, g_y_0_z_0_y_xxxxzz, g_y_0_z_0_y_xxxyyy, g_y_0_z_0_y_xxxyyz, g_y_0_z_0_y_xxxyzz, g_y_0_z_0_y_xxxzzz, g_y_0_z_0_y_xxyyyy, g_y_0_z_0_y_xxyyyz, g_y_0_z_0_y_xxyyzz, g_y_0_z_0_y_xxyzzz, g_y_0_z_0_y_xxzzzz, g_y_0_z_0_y_xyyyyy, g_y_0_z_0_y_xyyyyz, g_y_0_z_0_y_xyyyzz, g_y_0_z_0_y_xyyzzz, g_y_0_z_0_y_xyzzzz, g_y_0_z_0_y_xzzzzz, g_y_0_z_0_y_yyyyyy, g_y_0_z_0_y_yyyyyz, g_y_0_z_0_y_yyyyzz, g_y_0_z_0_y_yyyzzz, g_y_0_z_0_y_yyzzzz, g_y_0_z_0_y_yzzzzz, g_y_0_z_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_y_xxxxxx[k] = -g_0_0_z_0_0_xxxxxx[k] - g_y_0_z_0_0_xxxxxx[k] * ab_y + g_y_0_z_0_0_xxxxxxy[k];

                g_y_0_z_0_y_xxxxxy[k] = -g_0_0_z_0_0_xxxxxy[k] - g_y_0_z_0_0_xxxxxy[k] * ab_y + g_y_0_z_0_0_xxxxxyy[k];

                g_y_0_z_0_y_xxxxxz[k] = -g_0_0_z_0_0_xxxxxz[k] - g_y_0_z_0_0_xxxxxz[k] * ab_y + g_y_0_z_0_0_xxxxxyz[k];

                g_y_0_z_0_y_xxxxyy[k] = -g_0_0_z_0_0_xxxxyy[k] - g_y_0_z_0_0_xxxxyy[k] * ab_y + g_y_0_z_0_0_xxxxyyy[k];

                g_y_0_z_0_y_xxxxyz[k] = -g_0_0_z_0_0_xxxxyz[k] - g_y_0_z_0_0_xxxxyz[k] * ab_y + g_y_0_z_0_0_xxxxyyz[k];

                g_y_0_z_0_y_xxxxzz[k] = -g_0_0_z_0_0_xxxxzz[k] - g_y_0_z_0_0_xxxxzz[k] * ab_y + g_y_0_z_0_0_xxxxyzz[k];

                g_y_0_z_0_y_xxxyyy[k] = -g_0_0_z_0_0_xxxyyy[k] - g_y_0_z_0_0_xxxyyy[k] * ab_y + g_y_0_z_0_0_xxxyyyy[k];

                g_y_0_z_0_y_xxxyyz[k] = -g_0_0_z_0_0_xxxyyz[k] - g_y_0_z_0_0_xxxyyz[k] * ab_y + g_y_0_z_0_0_xxxyyyz[k];

                g_y_0_z_0_y_xxxyzz[k] = -g_0_0_z_0_0_xxxyzz[k] - g_y_0_z_0_0_xxxyzz[k] * ab_y + g_y_0_z_0_0_xxxyyzz[k];

                g_y_0_z_0_y_xxxzzz[k] = -g_0_0_z_0_0_xxxzzz[k] - g_y_0_z_0_0_xxxzzz[k] * ab_y + g_y_0_z_0_0_xxxyzzz[k];

                g_y_0_z_0_y_xxyyyy[k] = -g_0_0_z_0_0_xxyyyy[k] - g_y_0_z_0_0_xxyyyy[k] * ab_y + g_y_0_z_0_0_xxyyyyy[k];

                g_y_0_z_0_y_xxyyyz[k] = -g_0_0_z_0_0_xxyyyz[k] - g_y_0_z_0_0_xxyyyz[k] * ab_y + g_y_0_z_0_0_xxyyyyz[k];

                g_y_0_z_0_y_xxyyzz[k] = -g_0_0_z_0_0_xxyyzz[k] - g_y_0_z_0_0_xxyyzz[k] * ab_y + g_y_0_z_0_0_xxyyyzz[k];

                g_y_0_z_0_y_xxyzzz[k] = -g_0_0_z_0_0_xxyzzz[k] - g_y_0_z_0_0_xxyzzz[k] * ab_y + g_y_0_z_0_0_xxyyzzz[k];

                g_y_0_z_0_y_xxzzzz[k] = -g_0_0_z_0_0_xxzzzz[k] - g_y_0_z_0_0_xxzzzz[k] * ab_y + g_y_0_z_0_0_xxyzzzz[k];

                g_y_0_z_0_y_xyyyyy[k] = -g_0_0_z_0_0_xyyyyy[k] - g_y_0_z_0_0_xyyyyy[k] * ab_y + g_y_0_z_0_0_xyyyyyy[k];

                g_y_0_z_0_y_xyyyyz[k] = -g_0_0_z_0_0_xyyyyz[k] - g_y_0_z_0_0_xyyyyz[k] * ab_y + g_y_0_z_0_0_xyyyyyz[k];

                g_y_0_z_0_y_xyyyzz[k] = -g_0_0_z_0_0_xyyyzz[k] - g_y_0_z_0_0_xyyyzz[k] * ab_y + g_y_0_z_0_0_xyyyyzz[k];

                g_y_0_z_0_y_xyyzzz[k] = -g_0_0_z_0_0_xyyzzz[k] - g_y_0_z_0_0_xyyzzz[k] * ab_y + g_y_0_z_0_0_xyyyzzz[k];

                g_y_0_z_0_y_xyzzzz[k] = -g_0_0_z_0_0_xyzzzz[k] - g_y_0_z_0_0_xyzzzz[k] * ab_y + g_y_0_z_0_0_xyyzzzz[k];

                g_y_0_z_0_y_xzzzzz[k] = -g_0_0_z_0_0_xzzzzz[k] - g_y_0_z_0_0_xzzzzz[k] * ab_y + g_y_0_z_0_0_xyzzzzz[k];

                g_y_0_z_0_y_yyyyyy[k] = -g_0_0_z_0_0_yyyyyy[k] - g_y_0_z_0_0_yyyyyy[k] * ab_y + g_y_0_z_0_0_yyyyyyy[k];

                g_y_0_z_0_y_yyyyyz[k] = -g_0_0_z_0_0_yyyyyz[k] - g_y_0_z_0_0_yyyyyz[k] * ab_y + g_y_0_z_0_0_yyyyyyz[k];

                g_y_0_z_0_y_yyyyzz[k] = -g_0_0_z_0_0_yyyyzz[k] - g_y_0_z_0_0_yyyyzz[k] * ab_y + g_y_0_z_0_0_yyyyyzz[k];

                g_y_0_z_0_y_yyyzzz[k] = -g_0_0_z_0_0_yyyzzz[k] - g_y_0_z_0_0_yyyzzz[k] * ab_y + g_y_0_z_0_0_yyyyzzz[k];

                g_y_0_z_0_y_yyzzzz[k] = -g_0_0_z_0_0_yyzzzz[k] - g_y_0_z_0_0_yyzzzz[k] * ab_y + g_y_0_z_0_0_yyyzzzz[k];

                g_y_0_z_0_y_yzzzzz[k] = -g_0_0_z_0_0_yzzzzz[k] - g_y_0_z_0_0_yzzzzz[k] * ab_y + g_y_0_z_0_0_yyzzzzz[k];

                g_y_0_z_0_y_zzzzzz[k] = -g_0_0_z_0_0_zzzzzz[k] - g_y_0_z_0_0_zzzzzz[k] * ab_y + g_y_0_z_0_0_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_y_0_z_0_z_xxxxxx = cbuffer.data(pi_geom_1010_off + 476 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxxxy = cbuffer.data(pi_geom_1010_off + 477 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxxxz = cbuffer.data(pi_geom_1010_off + 478 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxxyy = cbuffer.data(pi_geom_1010_off + 479 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxxyz = cbuffer.data(pi_geom_1010_off + 480 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxxzz = cbuffer.data(pi_geom_1010_off + 481 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxyyy = cbuffer.data(pi_geom_1010_off + 482 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxyyz = cbuffer.data(pi_geom_1010_off + 483 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxyzz = cbuffer.data(pi_geom_1010_off + 484 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxxzzz = cbuffer.data(pi_geom_1010_off + 485 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxyyyy = cbuffer.data(pi_geom_1010_off + 486 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxyyyz = cbuffer.data(pi_geom_1010_off + 487 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxyyzz = cbuffer.data(pi_geom_1010_off + 488 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxyzzz = cbuffer.data(pi_geom_1010_off + 489 * ccomps * dcomps);

            auto g_y_0_z_0_z_xxzzzz = cbuffer.data(pi_geom_1010_off + 490 * ccomps * dcomps);

            auto g_y_0_z_0_z_xyyyyy = cbuffer.data(pi_geom_1010_off + 491 * ccomps * dcomps);

            auto g_y_0_z_0_z_xyyyyz = cbuffer.data(pi_geom_1010_off + 492 * ccomps * dcomps);

            auto g_y_0_z_0_z_xyyyzz = cbuffer.data(pi_geom_1010_off + 493 * ccomps * dcomps);

            auto g_y_0_z_0_z_xyyzzz = cbuffer.data(pi_geom_1010_off + 494 * ccomps * dcomps);

            auto g_y_0_z_0_z_xyzzzz = cbuffer.data(pi_geom_1010_off + 495 * ccomps * dcomps);

            auto g_y_0_z_0_z_xzzzzz = cbuffer.data(pi_geom_1010_off + 496 * ccomps * dcomps);

            auto g_y_0_z_0_z_yyyyyy = cbuffer.data(pi_geom_1010_off + 497 * ccomps * dcomps);

            auto g_y_0_z_0_z_yyyyyz = cbuffer.data(pi_geom_1010_off + 498 * ccomps * dcomps);

            auto g_y_0_z_0_z_yyyyzz = cbuffer.data(pi_geom_1010_off + 499 * ccomps * dcomps);

            auto g_y_0_z_0_z_yyyzzz = cbuffer.data(pi_geom_1010_off + 500 * ccomps * dcomps);

            auto g_y_0_z_0_z_yyzzzz = cbuffer.data(pi_geom_1010_off + 501 * ccomps * dcomps);

            auto g_y_0_z_0_z_yzzzzz = cbuffer.data(pi_geom_1010_off + 502 * ccomps * dcomps);

            auto g_y_0_z_0_z_zzzzzz = cbuffer.data(pi_geom_1010_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_0_0_xxxxxx, g_y_0_z_0_0_xxxxxxz, g_y_0_z_0_0_xxxxxy, g_y_0_z_0_0_xxxxxyz, g_y_0_z_0_0_xxxxxz, g_y_0_z_0_0_xxxxxzz, g_y_0_z_0_0_xxxxyy, g_y_0_z_0_0_xxxxyyz, g_y_0_z_0_0_xxxxyz, g_y_0_z_0_0_xxxxyzz, g_y_0_z_0_0_xxxxzz, g_y_0_z_0_0_xxxxzzz, g_y_0_z_0_0_xxxyyy, g_y_0_z_0_0_xxxyyyz, g_y_0_z_0_0_xxxyyz, g_y_0_z_0_0_xxxyyzz, g_y_0_z_0_0_xxxyzz, g_y_0_z_0_0_xxxyzzz, g_y_0_z_0_0_xxxzzz, g_y_0_z_0_0_xxxzzzz, g_y_0_z_0_0_xxyyyy, g_y_0_z_0_0_xxyyyyz, g_y_0_z_0_0_xxyyyz, g_y_0_z_0_0_xxyyyzz, g_y_0_z_0_0_xxyyzz, g_y_0_z_0_0_xxyyzzz, g_y_0_z_0_0_xxyzzz, g_y_0_z_0_0_xxyzzzz, g_y_0_z_0_0_xxzzzz, g_y_0_z_0_0_xxzzzzz, g_y_0_z_0_0_xyyyyy, g_y_0_z_0_0_xyyyyyz, g_y_0_z_0_0_xyyyyz, g_y_0_z_0_0_xyyyyzz, g_y_0_z_0_0_xyyyzz, g_y_0_z_0_0_xyyyzzz, g_y_0_z_0_0_xyyzzz, g_y_0_z_0_0_xyyzzzz, g_y_0_z_0_0_xyzzzz, g_y_0_z_0_0_xyzzzzz, g_y_0_z_0_0_xzzzzz, g_y_0_z_0_0_xzzzzzz, g_y_0_z_0_0_yyyyyy, g_y_0_z_0_0_yyyyyyz, g_y_0_z_0_0_yyyyyz, g_y_0_z_0_0_yyyyyzz, g_y_0_z_0_0_yyyyzz, g_y_0_z_0_0_yyyyzzz, g_y_0_z_0_0_yyyzzz, g_y_0_z_0_0_yyyzzzz, g_y_0_z_0_0_yyzzzz, g_y_0_z_0_0_yyzzzzz, g_y_0_z_0_0_yzzzzz, g_y_0_z_0_0_yzzzzzz, g_y_0_z_0_0_zzzzzz, g_y_0_z_0_0_zzzzzzz, g_y_0_z_0_z_xxxxxx, g_y_0_z_0_z_xxxxxy, g_y_0_z_0_z_xxxxxz, g_y_0_z_0_z_xxxxyy, g_y_0_z_0_z_xxxxyz, g_y_0_z_0_z_xxxxzz, g_y_0_z_0_z_xxxyyy, g_y_0_z_0_z_xxxyyz, g_y_0_z_0_z_xxxyzz, g_y_0_z_0_z_xxxzzz, g_y_0_z_0_z_xxyyyy, g_y_0_z_0_z_xxyyyz, g_y_0_z_0_z_xxyyzz, g_y_0_z_0_z_xxyzzz, g_y_0_z_0_z_xxzzzz, g_y_0_z_0_z_xyyyyy, g_y_0_z_0_z_xyyyyz, g_y_0_z_0_z_xyyyzz, g_y_0_z_0_z_xyyzzz, g_y_0_z_0_z_xyzzzz, g_y_0_z_0_z_xzzzzz, g_y_0_z_0_z_yyyyyy, g_y_0_z_0_z_yyyyyz, g_y_0_z_0_z_yyyyzz, g_y_0_z_0_z_yyyzzz, g_y_0_z_0_z_yyzzzz, g_y_0_z_0_z_yzzzzz, g_y_0_z_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_z_0_z_xxxxxx[k] = -g_y_0_z_0_0_xxxxxx[k] * ab_z + g_y_0_z_0_0_xxxxxxz[k];

                g_y_0_z_0_z_xxxxxy[k] = -g_y_0_z_0_0_xxxxxy[k] * ab_z + g_y_0_z_0_0_xxxxxyz[k];

                g_y_0_z_0_z_xxxxxz[k] = -g_y_0_z_0_0_xxxxxz[k] * ab_z + g_y_0_z_0_0_xxxxxzz[k];

                g_y_0_z_0_z_xxxxyy[k] = -g_y_0_z_0_0_xxxxyy[k] * ab_z + g_y_0_z_0_0_xxxxyyz[k];

                g_y_0_z_0_z_xxxxyz[k] = -g_y_0_z_0_0_xxxxyz[k] * ab_z + g_y_0_z_0_0_xxxxyzz[k];

                g_y_0_z_0_z_xxxxzz[k] = -g_y_0_z_0_0_xxxxzz[k] * ab_z + g_y_0_z_0_0_xxxxzzz[k];

                g_y_0_z_0_z_xxxyyy[k] = -g_y_0_z_0_0_xxxyyy[k] * ab_z + g_y_0_z_0_0_xxxyyyz[k];

                g_y_0_z_0_z_xxxyyz[k] = -g_y_0_z_0_0_xxxyyz[k] * ab_z + g_y_0_z_0_0_xxxyyzz[k];

                g_y_0_z_0_z_xxxyzz[k] = -g_y_0_z_0_0_xxxyzz[k] * ab_z + g_y_0_z_0_0_xxxyzzz[k];

                g_y_0_z_0_z_xxxzzz[k] = -g_y_0_z_0_0_xxxzzz[k] * ab_z + g_y_0_z_0_0_xxxzzzz[k];

                g_y_0_z_0_z_xxyyyy[k] = -g_y_0_z_0_0_xxyyyy[k] * ab_z + g_y_0_z_0_0_xxyyyyz[k];

                g_y_0_z_0_z_xxyyyz[k] = -g_y_0_z_0_0_xxyyyz[k] * ab_z + g_y_0_z_0_0_xxyyyzz[k];

                g_y_0_z_0_z_xxyyzz[k] = -g_y_0_z_0_0_xxyyzz[k] * ab_z + g_y_0_z_0_0_xxyyzzz[k];

                g_y_0_z_0_z_xxyzzz[k] = -g_y_0_z_0_0_xxyzzz[k] * ab_z + g_y_0_z_0_0_xxyzzzz[k];

                g_y_0_z_0_z_xxzzzz[k] = -g_y_0_z_0_0_xxzzzz[k] * ab_z + g_y_0_z_0_0_xxzzzzz[k];

                g_y_0_z_0_z_xyyyyy[k] = -g_y_0_z_0_0_xyyyyy[k] * ab_z + g_y_0_z_0_0_xyyyyyz[k];

                g_y_0_z_0_z_xyyyyz[k] = -g_y_0_z_0_0_xyyyyz[k] * ab_z + g_y_0_z_0_0_xyyyyzz[k];

                g_y_0_z_0_z_xyyyzz[k] = -g_y_0_z_0_0_xyyyzz[k] * ab_z + g_y_0_z_0_0_xyyyzzz[k];

                g_y_0_z_0_z_xyyzzz[k] = -g_y_0_z_0_0_xyyzzz[k] * ab_z + g_y_0_z_0_0_xyyzzzz[k];

                g_y_0_z_0_z_xyzzzz[k] = -g_y_0_z_0_0_xyzzzz[k] * ab_z + g_y_0_z_0_0_xyzzzzz[k];

                g_y_0_z_0_z_xzzzzz[k] = -g_y_0_z_0_0_xzzzzz[k] * ab_z + g_y_0_z_0_0_xzzzzzz[k];

                g_y_0_z_0_z_yyyyyy[k] = -g_y_0_z_0_0_yyyyyy[k] * ab_z + g_y_0_z_0_0_yyyyyyz[k];

                g_y_0_z_0_z_yyyyyz[k] = -g_y_0_z_0_0_yyyyyz[k] * ab_z + g_y_0_z_0_0_yyyyyzz[k];

                g_y_0_z_0_z_yyyyzz[k] = -g_y_0_z_0_0_yyyyzz[k] * ab_z + g_y_0_z_0_0_yyyyzzz[k];

                g_y_0_z_0_z_yyyzzz[k] = -g_y_0_z_0_0_yyyzzz[k] * ab_z + g_y_0_z_0_0_yyyzzzz[k];

                g_y_0_z_0_z_yyzzzz[k] = -g_y_0_z_0_0_yyzzzz[k] * ab_z + g_y_0_z_0_0_yyzzzzz[k];

                g_y_0_z_0_z_yzzzzz[k] = -g_y_0_z_0_0_yzzzzz[k] * ab_z + g_y_0_z_0_0_yzzzzzz[k];

                g_y_0_z_0_z_zzzzzz[k] = -g_y_0_z_0_0_zzzzzz[k] * ab_z + g_y_0_z_0_0_zzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_x_xxxxxx = cbuffer.data(pi_geom_1010_off + 504 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxxxy = cbuffer.data(pi_geom_1010_off + 505 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxxxz = cbuffer.data(pi_geom_1010_off + 506 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxxyy = cbuffer.data(pi_geom_1010_off + 507 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxxyz = cbuffer.data(pi_geom_1010_off + 508 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxxzz = cbuffer.data(pi_geom_1010_off + 509 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxyyy = cbuffer.data(pi_geom_1010_off + 510 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxyyz = cbuffer.data(pi_geom_1010_off + 511 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxyzz = cbuffer.data(pi_geom_1010_off + 512 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxxzzz = cbuffer.data(pi_geom_1010_off + 513 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxyyyy = cbuffer.data(pi_geom_1010_off + 514 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxyyyz = cbuffer.data(pi_geom_1010_off + 515 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxyyzz = cbuffer.data(pi_geom_1010_off + 516 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxyzzz = cbuffer.data(pi_geom_1010_off + 517 * ccomps * dcomps);

            auto g_z_0_x_0_x_xxzzzz = cbuffer.data(pi_geom_1010_off + 518 * ccomps * dcomps);

            auto g_z_0_x_0_x_xyyyyy = cbuffer.data(pi_geom_1010_off + 519 * ccomps * dcomps);

            auto g_z_0_x_0_x_xyyyyz = cbuffer.data(pi_geom_1010_off + 520 * ccomps * dcomps);

            auto g_z_0_x_0_x_xyyyzz = cbuffer.data(pi_geom_1010_off + 521 * ccomps * dcomps);

            auto g_z_0_x_0_x_xyyzzz = cbuffer.data(pi_geom_1010_off + 522 * ccomps * dcomps);

            auto g_z_0_x_0_x_xyzzzz = cbuffer.data(pi_geom_1010_off + 523 * ccomps * dcomps);

            auto g_z_0_x_0_x_xzzzzz = cbuffer.data(pi_geom_1010_off + 524 * ccomps * dcomps);

            auto g_z_0_x_0_x_yyyyyy = cbuffer.data(pi_geom_1010_off + 525 * ccomps * dcomps);

            auto g_z_0_x_0_x_yyyyyz = cbuffer.data(pi_geom_1010_off + 526 * ccomps * dcomps);

            auto g_z_0_x_0_x_yyyyzz = cbuffer.data(pi_geom_1010_off + 527 * ccomps * dcomps);

            auto g_z_0_x_0_x_yyyzzz = cbuffer.data(pi_geom_1010_off + 528 * ccomps * dcomps);

            auto g_z_0_x_0_x_yyzzzz = cbuffer.data(pi_geom_1010_off + 529 * ccomps * dcomps);

            auto g_z_0_x_0_x_yzzzzz = cbuffer.data(pi_geom_1010_off + 530 * ccomps * dcomps);

            auto g_z_0_x_0_x_zzzzzz = cbuffer.data(pi_geom_1010_off + 531 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_0_xxxxxx, g_z_0_x_0_0_xxxxxxx, g_z_0_x_0_0_xxxxxxy, g_z_0_x_0_0_xxxxxxz, g_z_0_x_0_0_xxxxxy, g_z_0_x_0_0_xxxxxyy, g_z_0_x_0_0_xxxxxyz, g_z_0_x_0_0_xxxxxz, g_z_0_x_0_0_xxxxxzz, g_z_0_x_0_0_xxxxyy, g_z_0_x_0_0_xxxxyyy, g_z_0_x_0_0_xxxxyyz, g_z_0_x_0_0_xxxxyz, g_z_0_x_0_0_xxxxyzz, g_z_0_x_0_0_xxxxzz, g_z_0_x_0_0_xxxxzzz, g_z_0_x_0_0_xxxyyy, g_z_0_x_0_0_xxxyyyy, g_z_0_x_0_0_xxxyyyz, g_z_0_x_0_0_xxxyyz, g_z_0_x_0_0_xxxyyzz, g_z_0_x_0_0_xxxyzz, g_z_0_x_0_0_xxxyzzz, g_z_0_x_0_0_xxxzzz, g_z_0_x_0_0_xxxzzzz, g_z_0_x_0_0_xxyyyy, g_z_0_x_0_0_xxyyyyy, g_z_0_x_0_0_xxyyyyz, g_z_0_x_0_0_xxyyyz, g_z_0_x_0_0_xxyyyzz, g_z_0_x_0_0_xxyyzz, g_z_0_x_0_0_xxyyzzz, g_z_0_x_0_0_xxyzzz, g_z_0_x_0_0_xxyzzzz, g_z_0_x_0_0_xxzzzz, g_z_0_x_0_0_xxzzzzz, g_z_0_x_0_0_xyyyyy, g_z_0_x_0_0_xyyyyyy, g_z_0_x_0_0_xyyyyyz, g_z_0_x_0_0_xyyyyz, g_z_0_x_0_0_xyyyyzz, g_z_0_x_0_0_xyyyzz, g_z_0_x_0_0_xyyyzzz, g_z_0_x_0_0_xyyzzz, g_z_0_x_0_0_xyyzzzz, g_z_0_x_0_0_xyzzzz, g_z_0_x_0_0_xyzzzzz, g_z_0_x_0_0_xzzzzz, g_z_0_x_0_0_xzzzzzz, g_z_0_x_0_0_yyyyyy, g_z_0_x_0_0_yyyyyz, g_z_0_x_0_0_yyyyzz, g_z_0_x_0_0_yyyzzz, g_z_0_x_0_0_yyzzzz, g_z_0_x_0_0_yzzzzz, g_z_0_x_0_0_zzzzzz, g_z_0_x_0_x_xxxxxx, g_z_0_x_0_x_xxxxxy, g_z_0_x_0_x_xxxxxz, g_z_0_x_0_x_xxxxyy, g_z_0_x_0_x_xxxxyz, g_z_0_x_0_x_xxxxzz, g_z_0_x_0_x_xxxyyy, g_z_0_x_0_x_xxxyyz, g_z_0_x_0_x_xxxyzz, g_z_0_x_0_x_xxxzzz, g_z_0_x_0_x_xxyyyy, g_z_0_x_0_x_xxyyyz, g_z_0_x_0_x_xxyyzz, g_z_0_x_0_x_xxyzzz, g_z_0_x_0_x_xxzzzz, g_z_0_x_0_x_xyyyyy, g_z_0_x_0_x_xyyyyz, g_z_0_x_0_x_xyyyzz, g_z_0_x_0_x_xyyzzz, g_z_0_x_0_x_xyzzzz, g_z_0_x_0_x_xzzzzz, g_z_0_x_0_x_yyyyyy, g_z_0_x_0_x_yyyyyz, g_z_0_x_0_x_yyyyzz, g_z_0_x_0_x_yyyzzz, g_z_0_x_0_x_yyzzzz, g_z_0_x_0_x_yzzzzz, g_z_0_x_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_x_xxxxxx[k] = -g_z_0_x_0_0_xxxxxx[k] * ab_x + g_z_0_x_0_0_xxxxxxx[k];

                g_z_0_x_0_x_xxxxxy[k] = -g_z_0_x_0_0_xxxxxy[k] * ab_x + g_z_0_x_0_0_xxxxxxy[k];

                g_z_0_x_0_x_xxxxxz[k] = -g_z_0_x_0_0_xxxxxz[k] * ab_x + g_z_0_x_0_0_xxxxxxz[k];

                g_z_0_x_0_x_xxxxyy[k] = -g_z_0_x_0_0_xxxxyy[k] * ab_x + g_z_0_x_0_0_xxxxxyy[k];

                g_z_0_x_0_x_xxxxyz[k] = -g_z_0_x_0_0_xxxxyz[k] * ab_x + g_z_0_x_0_0_xxxxxyz[k];

                g_z_0_x_0_x_xxxxzz[k] = -g_z_0_x_0_0_xxxxzz[k] * ab_x + g_z_0_x_0_0_xxxxxzz[k];

                g_z_0_x_0_x_xxxyyy[k] = -g_z_0_x_0_0_xxxyyy[k] * ab_x + g_z_0_x_0_0_xxxxyyy[k];

                g_z_0_x_0_x_xxxyyz[k] = -g_z_0_x_0_0_xxxyyz[k] * ab_x + g_z_0_x_0_0_xxxxyyz[k];

                g_z_0_x_0_x_xxxyzz[k] = -g_z_0_x_0_0_xxxyzz[k] * ab_x + g_z_0_x_0_0_xxxxyzz[k];

                g_z_0_x_0_x_xxxzzz[k] = -g_z_0_x_0_0_xxxzzz[k] * ab_x + g_z_0_x_0_0_xxxxzzz[k];

                g_z_0_x_0_x_xxyyyy[k] = -g_z_0_x_0_0_xxyyyy[k] * ab_x + g_z_0_x_0_0_xxxyyyy[k];

                g_z_0_x_0_x_xxyyyz[k] = -g_z_0_x_0_0_xxyyyz[k] * ab_x + g_z_0_x_0_0_xxxyyyz[k];

                g_z_0_x_0_x_xxyyzz[k] = -g_z_0_x_0_0_xxyyzz[k] * ab_x + g_z_0_x_0_0_xxxyyzz[k];

                g_z_0_x_0_x_xxyzzz[k] = -g_z_0_x_0_0_xxyzzz[k] * ab_x + g_z_0_x_0_0_xxxyzzz[k];

                g_z_0_x_0_x_xxzzzz[k] = -g_z_0_x_0_0_xxzzzz[k] * ab_x + g_z_0_x_0_0_xxxzzzz[k];

                g_z_0_x_0_x_xyyyyy[k] = -g_z_0_x_0_0_xyyyyy[k] * ab_x + g_z_0_x_0_0_xxyyyyy[k];

                g_z_0_x_0_x_xyyyyz[k] = -g_z_0_x_0_0_xyyyyz[k] * ab_x + g_z_0_x_0_0_xxyyyyz[k];

                g_z_0_x_0_x_xyyyzz[k] = -g_z_0_x_0_0_xyyyzz[k] * ab_x + g_z_0_x_0_0_xxyyyzz[k];

                g_z_0_x_0_x_xyyzzz[k] = -g_z_0_x_0_0_xyyzzz[k] * ab_x + g_z_0_x_0_0_xxyyzzz[k];

                g_z_0_x_0_x_xyzzzz[k] = -g_z_0_x_0_0_xyzzzz[k] * ab_x + g_z_0_x_0_0_xxyzzzz[k];

                g_z_0_x_0_x_xzzzzz[k] = -g_z_0_x_0_0_xzzzzz[k] * ab_x + g_z_0_x_0_0_xxzzzzz[k];

                g_z_0_x_0_x_yyyyyy[k] = -g_z_0_x_0_0_yyyyyy[k] * ab_x + g_z_0_x_0_0_xyyyyyy[k];

                g_z_0_x_0_x_yyyyyz[k] = -g_z_0_x_0_0_yyyyyz[k] * ab_x + g_z_0_x_0_0_xyyyyyz[k];

                g_z_0_x_0_x_yyyyzz[k] = -g_z_0_x_0_0_yyyyzz[k] * ab_x + g_z_0_x_0_0_xyyyyzz[k];

                g_z_0_x_0_x_yyyzzz[k] = -g_z_0_x_0_0_yyyzzz[k] * ab_x + g_z_0_x_0_0_xyyyzzz[k];

                g_z_0_x_0_x_yyzzzz[k] = -g_z_0_x_0_0_yyzzzz[k] * ab_x + g_z_0_x_0_0_xyyzzzz[k];

                g_z_0_x_0_x_yzzzzz[k] = -g_z_0_x_0_0_yzzzzz[k] * ab_x + g_z_0_x_0_0_xyzzzzz[k];

                g_z_0_x_0_x_zzzzzz[k] = -g_z_0_x_0_0_zzzzzz[k] * ab_x + g_z_0_x_0_0_xzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_y_xxxxxx = cbuffer.data(pi_geom_1010_off + 532 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxxxy = cbuffer.data(pi_geom_1010_off + 533 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxxxz = cbuffer.data(pi_geom_1010_off + 534 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxxyy = cbuffer.data(pi_geom_1010_off + 535 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxxyz = cbuffer.data(pi_geom_1010_off + 536 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxxzz = cbuffer.data(pi_geom_1010_off + 537 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxyyy = cbuffer.data(pi_geom_1010_off + 538 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxyyz = cbuffer.data(pi_geom_1010_off + 539 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxyzz = cbuffer.data(pi_geom_1010_off + 540 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxxzzz = cbuffer.data(pi_geom_1010_off + 541 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxyyyy = cbuffer.data(pi_geom_1010_off + 542 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxyyyz = cbuffer.data(pi_geom_1010_off + 543 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxyyzz = cbuffer.data(pi_geom_1010_off + 544 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxyzzz = cbuffer.data(pi_geom_1010_off + 545 * ccomps * dcomps);

            auto g_z_0_x_0_y_xxzzzz = cbuffer.data(pi_geom_1010_off + 546 * ccomps * dcomps);

            auto g_z_0_x_0_y_xyyyyy = cbuffer.data(pi_geom_1010_off + 547 * ccomps * dcomps);

            auto g_z_0_x_0_y_xyyyyz = cbuffer.data(pi_geom_1010_off + 548 * ccomps * dcomps);

            auto g_z_0_x_0_y_xyyyzz = cbuffer.data(pi_geom_1010_off + 549 * ccomps * dcomps);

            auto g_z_0_x_0_y_xyyzzz = cbuffer.data(pi_geom_1010_off + 550 * ccomps * dcomps);

            auto g_z_0_x_0_y_xyzzzz = cbuffer.data(pi_geom_1010_off + 551 * ccomps * dcomps);

            auto g_z_0_x_0_y_xzzzzz = cbuffer.data(pi_geom_1010_off + 552 * ccomps * dcomps);

            auto g_z_0_x_0_y_yyyyyy = cbuffer.data(pi_geom_1010_off + 553 * ccomps * dcomps);

            auto g_z_0_x_0_y_yyyyyz = cbuffer.data(pi_geom_1010_off + 554 * ccomps * dcomps);

            auto g_z_0_x_0_y_yyyyzz = cbuffer.data(pi_geom_1010_off + 555 * ccomps * dcomps);

            auto g_z_0_x_0_y_yyyzzz = cbuffer.data(pi_geom_1010_off + 556 * ccomps * dcomps);

            auto g_z_0_x_0_y_yyzzzz = cbuffer.data(pi_geom_1010_off + 557 * ccomps * dcomps);

            auto g_z_0_x_0_y_yzzzzz = cbuffer.data(pi_geom_1010_off + 558 * ccomps * dcomps);

            auto g_z_0_x_0_y_zzzzzz = cbuffer.data(pi_geom_1010_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_0_0_xxxxxx, g_z_0_x_0_0_xxxxxxy, g_z_0_x_0_0_xxxxxy, g_z_0_x_0_0_xxxxxyy, g_z_0_x_0_0_xxxxxyz, g_z_0_x_0_0_xxxxxz, g_z_0_x_0_0_xxxxyy, g_z_0_x_0_0_xxxxyyy, g_z_0_x_0_0_xxxxyyz, g_z_0_x_0_0_xxxxyz, g_z_0_x_0_0_xxxxyzz, g_z_0_x_0_0_xxxxzz, g_z_0_x_0_0_xxxyyy, g_z_0_x_0_0_xxxyyyy, g_z_0_x_0_0_xxxyyyz, g_z_0_x_0_0_xxxyyz, g_z_0_x_0_0_xxxyyzz, g_z_0_x_0_0_xxxyzz, g_z_0_x_0_0_xxxyzzz, g_z_0_x_0_0_xxxzzz, g_z_0_x_0_0_xxyyyy, g_z_0_x_0_0_xxyyyyy, g_z_0_x_0_0_xxyyyyz, g_z_0_x_0_0_xxyyyz, g_z_0_x_0_0_xxyyyzz, g_z_0_x_0_0_xxyyzz, g_z_0_x_0_0_xxyyzzz, g_z_0_x_0_0_xxyzzz, g_z_0_x_0_0_xxyzzzz, g_z_0_x_0_0_xxzzzz, g_z_0_x_0_0_xyyyyy, g_z_0_x_0_0_xyyyyyy, g_z_0_x_0_0_xyyyyyz, g_z_0_x_0_0_xyyyyz, g_z_0_x_0_0_xyyyyzz, g_z_0_x_0_0_xyyyzz, g_z_0_x_0_0_xyyyzzz, g_z_0_x_0_0_xyyzzz, g_z_0_x_0_0_xyyzzzz, g_z_0_x_0_0_xyzzzz, g_z_0_x_0_0_xyzzzzz, g_z_0_x_0_0_xzzzzz, g_z_0_x_0_0_yyyyyy, g_z_0_x_0_0_yyyyyyy, g_z_0_x_0_0_yyyyyyz, g_z_0_x_0_0_yyyyyz, g_z_0_x_0_0_yyyyyzz, g_z_0_x_0_0_yyyyzz, g_z_0_x_0_0_yyyyzzz, g_z_0_x_0_0_yyyzzz, g_z_0_x_0_0_yyyzzzz, g_z_0_x_0_0_yyzzzz, g_z_0_x_0_0_yyzzzzz, g_z_0_x_0_0_yzzzzz, g_z_0_x_0_0_yzzzzzz, g_z_0_x_0_0_zzzzzz, g_z_0_x_0_y_xxxxxx, g_z_0_x_0_y_xxxxxy, g_z_0_x_0_y_xxxxxz, g_z_0_x_0_y_xxxxyy, g_z_0_x_0_y_xxxxyz, g_z_0_x_0_y_xxxxzz, g_z_0_x_0_y_xxxyyy, g_z_0_x_0_y_xxxyyz, g_z_0_x_0_y_xxxyzz, g_z_0_x_0_y_xxxzzz, g_z_0_x_0_y_xxyyyy, g_z_0_x_0_y_xxyyyz, g_z_0_x_0_y_xxyyzz, g_z_0_x_0_y_xxyzzz, g_z_0_x_0_y_xxzzzz, g_z_0_x_0_y_xyyyyy, g_z_0_x_0_y_xyyyyz, g_z_0_x_0_y_xyyyzz, g_z_0_x_0_y_xyyzzz, g_z_0_x_0_y_xyzzzz, g_z_0_x_0_y_xzzzzz, g_z_0_x_0_y_yyyyyy, g_z_0_x_0_y_yyyyyz, g_z_0_x_0_y_yyyyzz, g_z_0_x_0_y_yyyzzz, g_z_0_x_0_y_yyzzzz, g_z_0_x_0_y_yzzzzz, g_z_0_x_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_y_xxxxxx[k] = -g_z_0_x_0_0_xxxxxx[k] * ab_y + g_z_0_x_0_0_xxxxxxy[k];

                g_z_0_x_0_y_xxxxxy[k] = -g_z_0_x_0_0_xxxxxy[k] * ab_y + g_z_0_x_0_0_xxxxxyy[k];

                g_z_0_x_0_y_xxxxxz[k] = -g_z_0_x_0_0_xxxxxz[k] * ab_y + g_z_0_x_0_0_xxxxxyz[k];

                g_z_0_x_0_y_xxxxyy[k] = -g_z_0_x_0_0_xxxxyy[k] * ab_y + g_z_0_x_0_0_xxxxyyy[k];

                g_z_0_x_0_y_xxxxyz[k] = -g_z_0_x_0_0_xxxxyz[k] * ab_y + g_z_0_x_0_0_xxxxyyz[k];

                g_z_0_x_0_y_xxxxzz[k] = -g_z_0_x_0_0_xxxxzz[k] * ab_y + g_z_0_x_0_0_xxxxyzz[k];

                g_z_0_x_0_y_xxxyyy[k] = -g_z_0_x_0_0_xxxyyy[k] * ab_y + g_z_0_x_0_0_xxxyyyy[k];

                g_z_0_x_0_y_xxxyyz[k] = -g_z_0_x_0_0_xxxyyz[k] * ab_y + g_z_0_x_0_0_xxxyyyz[k];

                g_z_0_x_0_y_xxxyzz[k] = -g_z_0_x_0_0_xxxyzz[k] * ab_y + g_z_0_x_0_0_xxxyyzz[k];

                g_z_0_x_0_y_xxxzzz[k] = -g_z_0_x_0_0_xxxzzz[k] * ab_y + g_z_0_x_0_0_xxxyzzz[k];

                g_z_0_x_0_y_xxyyyy[k] = -g_z_0_x_0_0_xxyyyy[k] * ab_y + g_z_0_x_0_0_xxyyyyy[k];

                g_z_0_x_0_y_xxyyyz[k] = -g_z_0_x_0_0_xxyyyz[k] * ab_y + g_z_0_x_0_0_xxyyyyz[k];

                g_z_0_x_0_y_xxyyzz[k] = -g_z_0_x_0_0_xxyyzz[k] * ab_y + g_z_0_x_0_0_xxyyyzz[k];

                g_z_0_x_0_y_xxyzzz[k] = -g_z_0_x_0_0_xxyzzz[k] * ab_y + g_z_0_x_0_0_xxyyzzz[k];

                g_z_0_x_0_y_xxzzzz[k] = -g_z_0_x_0_0_xxzzzz[k] * ab_y + g_z_0_x_0_0_xxyzzzz[k];

                g_z_0_x_0_y_xyyyyy[k] = -g_z_0_x_0_0_xyyyyy[k] * ab_y + g_z_0_x_0_0_xyyyyyy[k];

                g_z_0_x_0_y_xyyyyz[k] = -g_z_0_x_0_0_xyyyyz[k] * ab_y + g_z_0_x_0_0_xyyyyyz[k];

                g_z_0_x_0_y_xyyyzz[k] = -g_z_0_x_0_0_xyyyzz[k] * ab_y + g_z_0_x_0_0_xyyyyzz[k];

                g_z_0_x_0_y_xyyzzz[k] = -g_z_0_x_0_0_xyyzzz[k] * ab_y + g_z_0_x_0_0_xyyyzzz[k];

                g_z_0_x_0_y_xyzzzz[k] = -g_z_0_x_0_0_xyzzzz[k] * ab_y + g_z_0_x_0_0_xyyzzzz[k];

                g_z_0_x_0_y_xzzzzz[k] = -g_z_0_x_0_0_xzzzzz[k] * ab_y + g_z_0_x_0_0_xyzzzzz[k];

                g_z_0_x_0_y_yyyyyy[k] = -g_z_0_x_0_0_yyyyyy[k] * ab_y + g_z_0_x_0_0_yyyyyyy[k];

                g_z_0_x_0_y_yyyyyz[k] = -g_z_0_x_0_0_yyyyyz[k] * ab_y + g_z_0_x_0_0_yyyyyyz[k];

                g_z_0_x_0_y_yyyyzz[k] = -g_z_0_x_0_0_yyyyzz[k] * ab_y + g_z_0_x_0_0_yyyyyzz[k];

                g_z_0_x_0_y_yyyzzz[k] = -g_z_0_x_0_0_yyyzzz[k] * ab_y + g_z_0_x_0_0_yyyyzzz[k];

                g_z_0_x_0_y_yyzzzz[k] = -g_z_0_x_0_0_yyzzzz[k] * ab_y + g_z_0_x_0_0_yyyzzzz[k];

                g_z_0_x_0_y_yzzzzz[k] = -g_z_0_x_0_0_yzzzzz[k] * ab_y + g_z_0_x_0_0_yyzzzzz[k];

                g_z_0_x_0_y_zzzzzz[k] = -g_z_0_x_0_0_zzzzzz[k] * ab_y + g_z_0_x_0_0_yzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_z_0_x_0_z_xxxxxx = cbuffer.data(pi_geom_1010_off + 560 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxxxy = cbuffer.data(pi_geom_1010_off + 561 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxxxz = cbuffer.data(pi_geom_1010_off + 562 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxxyy = cbuffer.data(pi_geom_1010_off + 563 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxxyz = cbuffer.data(pi_geom_1010_off + 564 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxxzz = cbuffer.data(pi_geom_1010_off + 565 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxyyy = cbuffer.data(pi_geom_1010_off + 566 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxyyz = cbuffer.data(pi_geom_1010_off + 567 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxyzz = cbuffer.data(pi_geom_1010_off + 568 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxxzzz = cbuffer.data(pi_geom_1010_off + 569 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxyyyy = cbuffer.data(pi_geom_1010_off + 570 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxyyyz = cbuffer.data(pi_geom_1010_off + 571 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxyyzz = cbuffer.data(pi_geom_1010_off + 572 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxyzzz = cbuffer.data(pi_geom_1010_off + 573 * ccomps * dcomps);

            auto g_z_0_x_0_z_xxzzzz = cbuffer.data(pi_geom_1010_off + 574 * ccomps * dcomps);

            auto g_z_0_x_0_z_xyyyyy = cbuffer.data(pi_geom_1010_off + 575 * ccomps * dcomps);

            auto g_z_0_x_0_z_xyyyyz = cbuffer.data(pi_geom_1010_off + 576 * ccomps * dcomps);

            auto g_z_0_x_0_z_xyyyzz = cbuffer.data(pi_geom_1010_off + 577 * ccomps * dcomps);

            auto g_z_0_x_0_z_xyyzzz = cbuffer.data(pi_geom_1010_off + 578 * ccomps * dcomps);

            auto g_z_0_x_0_z_xyzzzz = cbuffer.data(pi_geom_1010_off + 579 * ccomps * dcomps);

            auto g_z_0_x_0_z_xzzzzz = cbuffer.data(pi_geom_1010_off + 580 * ccomps * dcomps);

            auto g_z_0_x_0_z_yyyyyy = cbuffer.data(pi_geom_1010_off + 581 * ccomps * dcomps);

            auto g_z_0_x_0_z_yyyyyz = cbuffer.data(pi_geom_1010_off + 582 * ccomps * dcomps);

            auto g_z_0_x_0_z_yyyyzz = cbuffer.data(pi_geom_1010_off + 583 * ccomps * dcomps);

            auto g_z_0_x_0_z_yyyzzz = cbuffer.data(pi_geom_1010_off + 584 * ccomps * dcomps);

            auto g_z_0_x_0_z_yyzzzz = cbuffer.data(pi_geom_1010_off + 585 * ccomps * dcomps);

            auto g_z_0_x_0_z_yzzzzz = cbuffer.data(pi_geom_1010_off + 586 * ccomps * dcomps);

            auto g_z_0_x_0_z_zzzzzz = cbuffer.data(pi_geom_1010_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_x_0_0_xxxxxx, g_0_0_x_0_0_xxxxxy, g_0_0_x_0_0_xxxxxz, g_0_0_x_0_0_xxxxyy, g_0_0_x_0_0_xxxxyz, g_0_0_x_0_0_xxxxzz, g_0_0_x_0_0_xxxyyy, g_0_0_x_0_0_xxxyyz, g_0_0_x_0_0_xxxyzz, g_0_0_x_0_0_xxxzzz, g_0_0_x_0_0_xxyyyy, g_0_0_x_0_0_xxyyyz, g_0_0_x_0_0_xxyyzz, g_0_0_x_0_0_xxyzzz, g_0_0_x_0_0_xxzzzz, g_0_0_x_0_0_xyyyyy, g_0_0_x_0_0_xyyyyz, g_0_0_x_0_0_xyyyzz, g_0_0_x_0_0_xyyzzz, g_0_0_x_0_0_xyzzzz, g_0_0_x_0_0_xzzzzz, g_0_0_x_0_0_yyyyyy, g_0_0_x_0_0_yyyyyz, g_0_0_x_0_0_yyyyzz, g_0_0_x_0_0_yyyzzz, g_0_0_x_0_0_yyzzzz, g_0_0_x_0_0_yzzzzz, g_0_0_x_0_0_zzzzzz, g_z_0_x_0_0_xxxxxx, g_z_0_x_0_0_xxxxxxz, g_z_0_x_0_0_xxxxxy, g_z_0_x_0_0_xxxxxyz, g_z_0_x_0_0_xxxxxz, g_z_0_x_0_0_xxxxxzz, g_z_0_x_0_0_xxxxyy, g_z_0_x_0_0_xxxxyyz, g_z_0_x_0_0_xxxxyz, g_z_0_x_0_0_xxxxyzz, g_z_0_x_0_0_xxxxzz, g_z_0_x_0_0_xxxxzzz, g_z_0_x_0_0_xxxyyy, g_z_0_x_0_0_xxxyyyz, g_z_0_x_0_0_xxxyyz, g_z_0_x_0_0_xxxyyzz, g_z_0_x_0_0_xxxyzz, g_z_0_x_0_0_xxxyzzz, g_z_0_x_0_0_xxxzzz, g_z_0_x_0_0_xxxzzzz, g_z_0_x_0_0_xxyyyy, g_z_0_x_0_0_xxyyyyz, g_z_0_x_0_0_xxyyyz, g_z_0_x_0_0_xxyyyzz, g_z_0_x_0_0_xxyyzz, g_z_0_x_0_0_xxyyzzz, g_z_0_x_0_0_xxyzzz, g_z_0_x_0_0_xxyzzzz, g_z_0_x_0_0_xxzzzz, g_z_0_x_0_0_xxzzzzz, g_z_0_x_0_0_xyyyyy, g_z_0_x_0_0_xyyyyyz, g_z_0_x_0_0_xyyyyz, g_z_0_x_0_0_xyyyyzz, g_z_0_x_0_0_xyyyzz, g_z_0_x_0_0_xyyyzzz, g_z_0_x_0_0_xyyzzz, g_z_0_x_0_0_xyyzzzz, g_z_0_x_0_0_xyzzzz, g_z_0_x_0_0_xyzzzzz, g_z_0_x_0_0_xzzzzz, g_z_0_x_0_0_xzzzzzz, g_z_0_x_0_0_yyyyyy, g_z_0_x_0_0_yyyyyyz, g_z_0_x_0_0_yyyyyz, g_z_0_x_0_0_yyyyyzz, g_z_0_x_0_0_yyyyzz, g_z_0_x_0_0_yyyyzzz, g_z_0_x_0_0_yyyzzz, g_z_0_x_0_0_yyyzzzz, g_z_0_x_0_0_yyzzzz, g_z_0_x_0_0_yyzzzzz, g_z_0_x_0_0_yzzzzz, g_z_0_x_0_0_yzzzzzz, g_z_0_x_0_0_zzzzzz, g_z_0_x_0_0_zzzzzzz, g_z_0_x_0_z_xxxxxx, g_z_0_x_0_z_xxxxxy, g_z_0_x_0_z_xxxxxz, g_z_0_x_0_z_xxxxyy, g_z_0_x_0_z_xxxxyz, g_z_0_x_0_z_xxxxzz, g_z_0_x_0_z_xxxyyy, g_z_0_x_0_z_xxxyyz, g_z_0_x_0_z_xxxyzz, g_z_0_x_0_z_xxxzzz, g_z_0_x_0_z_xxyyyy, g_z_0_x_0_z_xxyyyz, g_z_0_x_0_z_xxyyzz, g_z_0_x_0_z_xxyzzz, g_z_0_x_0_z_xxzzzz, g_z_0_x_0_z_xyyyyy, g_z_0_x_0_z_xyyyyz, g_z_0_x_0_z_xyyyzz, g_z_0_x_0_z_xyyzzz, g_z_0_x_0_z_xyzzzz, g_z_0_x_0_z_xzzzzz, g_z_0_x_0_z_yyyyyy, g_z_0_x_0_z_yyyyyz, g_z_0_x_0_z_yyyyzz, g_z_0_x_0_z_yyyzzz, g_z_0_x_0_z_yyzzzz, g_z_0_x_0_z_yzzzzz, g_z_0_x_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_x_0_z_xxxxxx[k] = -g_0_0_x_0_0_xxxxxx[k] - g_z_0_x_0_0_xxxxxx[k] * ab_z + g_z_0_x_0_0_xxxxxxz[k];

                g_z_0_x_0_z_xxxxxy[k] = -g_0_0_x_0_0_xxxxxy[k] - g_z_0_x_0_0_xxxxxy[k] * ab_z + g_z_0_x_0_0_xxxxxyz[k];

                g_z_0_x_0_z_xxxxxz[k] = -g_0_0_x_0_0_xxxxxz[k] - g_z_0_x_0_0_xxxxxz[k] * ab_z + g_z_0_x_0_0_xxxxxzz[k];

                g_z_0_x_0_z_xxxxyy[k] = -g_0_0_x_0_0_xxxxyy[k] - g_z_0_x_0_0_xxxxyy[k] * ab_z + g_z_0_x_0_0_xxxxyyz[k];

                g_z_0_x_0_z_xxxxyz[k] = -g_0_0_x_0_0_xxxxyz[k] - g_z_0_x_0_0_xxxxyz[k] * ab_z + g_z_0_x_0_0_xxxxyzz[k];

                g_z_0_x_0_z_xxxxzz[k] = -g_0_0_x_0_0_xxxxzz[k] - g_z_0_x_0_0_xxxxzz[k] * ab_z + g_z_0_x_0_0_xxxxzzz[k];

                g_z_0_x_0_z_xxxyyy[k] = -g_0_0_x_0_0_xxxyyy[k] - g_z_0_x_0_0_xxxyyy[k] * ab_z + g_z_0_x_0_0_xxxyyyz[k];

                g_z_0_x_0_z_xxxyyz[k] = -g_0_0_x_0_0_xxxyyz[k] - g_z_0_x_0_0_xxxyyz[k] * ab_z + g_z_0_x_0_0_xxxyyzz[k];

                g_z_0_x_0_z_xxxyzz[k] = -g_0_0_x_0_0_xxxyzz[k] - g_z_0_x_0_0_xxxyzz[k] * ab_z + g_z_0_x_0_0_xxxyzzz[k];

                g_z_0_x_0_z_xxxzzz[k] = -g_0_0_x_0_0_xxxzzz[k] - g_z_0_x_0_0_xxxzzz[k] * ab_z + g_z_0_x_0_0_xxxzzzz[k];

                g_z_0_x_0_z_xxyyyy[k] = -g_0_0_x_0_0_xxyyyy[k] - g_z_0_x_0_0_xxyyyy[k] * ab_z + g_z_0_x_0_0_xxyyyyz[k];

                g_z_0_x_0_z_xxyyyz[k] = -g_0_0_x_0_0_xxyyyz[k] - g_z_0_x_0_0_xxyyyz[k] * ab_z + g_z_0_x_0_0_xxyyyzz[k];

                g_z_0_x_0_z_xxyyzz[k] = -g_0_0_x_0_0_xxyyzz[k] - g_z_0_x_0_0_xxyyzz[k] * ab_z + g_z_0_x_0_0_xxyyzzz[k];

                g_z_0_x_0_z_xxyzzz[k] = -g_0_0_x_0_0_xxyzzz[k] - g_z_0_x_0_0_xxyzzz[k] * ab_z + g_z_0_x_0_0_xxyzzzz[k];

                g_z_0_x_0_z_xxzzzz[k] = -g_0_0_x_0_0_xxzzzz[k] - g_z_0_x_0_0_xxzzzz[k] * ab_z + g_z_0_x_0_0_xxzzzzz[k];

                g_z_0_x_0_z_xyyyyy[k] = -g_0_0_x_0_0_xyyyyy[k] - g_z_0_x_0_0_xyyyyy[k] * ab_z + g_z_0_x_0_0_xyyyyyz[k];

                g_z_0_x_0_z_xyyyyz[k] = -g_0_0_x_0_0_xyyyyz[k] - g_z_0_x_0_0_xyyyyz[k] * ab_z + g_z_0_x_0_0_xyyyyzz[k];

                g_z_0_x_0_z_xyyyzz[k] = -g_0_0_x_0_0_xyyyzz[k] - g_z_0_x_0_0_xyyyzz[k] * ab_z + g_z_0_x_0_0_xyyyzzz[k];

                g_z_0_x_0_z_xyyzzz[k] = -g_0_0_x_0_0_xyyzzz[k] - g_z_0_x_0_0_xyyzzz[k] * ab_z + g_z_0_x_0_0_xyyzzzz[k];

                g_z_0_x_0_z_xyzzzz[k] = -g_0_0_x_0_0_xyzzzz[k] - g_z_0_x_0_0_xyzzzz[k] * ab_z + g_z_0_x_0_0_xyzzzzz[k];

                g_z_0_x_0_z_xzzzzz[k] = -g_0_0_x_0_0_xzzzzz[k] - g_z_0_x_0_0_xzzzzz[k] * ab_z + g_z_0_x_0_0_xzzzzzz[k];

                g_z_0_x_0_z_yyyyyy[k] = -g_0_0_x_0_0_yyyyyy[k] - g_z_0_x_0_0_yyyyyy[k] * ab_z + g_z_0_x_0_0_yyyyyyz[k];

                g_z_0_x_0_z_yyyyyz[k] = -g_0_0_x_0_0_yyyyyz[k] - g_z_0_x_0_0_yyyyyz[k] * ab_z + g_z_0_x_0_0_yyyyyzz[k];

                g_z_0_x_0_z_yyyyzz[k] = -g_0_0_x_0_0_yyyyzz[k] - g_z_0_x_0_0_yyyyzz[k] * ab_z + g_z_0_x_0_0_yyyyzzz[k];

                g_z_0_x_0_z_yyyzzz[k] = -g_0_0_x_0_0_yyyzzz[k] - g_z_0_x_0_0_yyyzzz[k] * ab_z + g_z_0_x_0_0_yyyzzzz[k];

                g_z_0_x_0_z_yyzzzz[k] = -g_0_0_x_0_0_yyzzzz[k] - g_z_0_x_0_0_yyzzzz[k] * ab_z + g_z_0_x_0_0_yyzzzzz[k];

                g_z_0_x_0_z_yzzzzz[k] = -g_0_0_x_0_0_yzzzzz[k] - g_z_0_x_0_0_yzzzzz[k] * ab_z + g_z_0_x_0_0_yzzzzzz[k];

                g_z_0_x_0_z_zzzzzz[k] = -g_0_0_x_0_0_zzzzzz[k] - g_z_0_x_0_0_zzzzzz[k] * ab_z + g_z_0_x_0_0_zzzzzzz[k];
            }

            /// Set up 588-616 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_x_xxxxxx = cbuffer.data(pi_geom_1010_off + 588 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxxxy = cbuffer.data(pi_geom_1010_off + 589 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxxxz = cbuffer.data(pi_geom_1010_off + 590 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxxyy = cbuffer.data(pi_geom_1010_off + 591 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxxyz = cbuffer.data(pi_geom_1010_off + 592 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxxzz = cbuffer.data(pi_geom_1010_off + 593 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxyyy = cbuffer.data(pi_geom_1010_off + 594 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxyyz = cbuffer.data(pi_geom_1010_off + 595 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxyzz = cbuffer.data(pi_geom_1010_off + 596 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxxzzz = cbuffer.data(pi_geom_1010_off + 597 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxyyyy = cbuffer.data(pi_geom_1010_off + 598 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxyyyz = cbuffer.data(pi_geom_1010_off + 599 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxyyzz = cbuffer.data(pi_geom_1010_off + 600 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxyzzz = cbuffer.data(pi_geom_1010_off + 601 * ccomps * dcomps);

            auto g_z_0_y_0_x_xxzzzz = cbuffer.data(pi_geom_1010_off + 602 * ccomps * dcomps);

            auto g_z_0_y_0_x_xyyyyy = cbuffer.data(pi_geom_1010_off + 603 * ccomps * dcomps);

            auto g_z_0_y_0_x_xyyyyz = cbuffer.data(pi_geom_1010_off + 604 * ccomps * dcomps);

            auto g_z_0_y_0_x_xyyyzz = cbuffer.data(pi_geom_1010_off + 605 * ccomps * dcomps);

            auto g_z_0_y_0_x_xyyzzz = cbuffer.data(pi_geom_1010_off + 606 * ccomps * dcomps);

            auto g_z_0_y_0_x_xyzzzz = cbuffer.data(pi_geom_1010_off + 607 * ccomps * dcomps);

            auto g_z_0_y_0_x_xzzzzz = cbuffer.data(pi_geom_1010_off + 608 * ccomps * dcomps);

            auto g_z_0_y_0_x_yyyyyy = cbuffer.data(pi_geom_1010_off + 609 * ccomps * dcomps);

            auto g_z_0_y_0_x_yyyyyz = cbuffer.data(pi_geom_1010_off + 610 * ccomps * dcomps);

            auto g_z_0_y_0_x_yyyyzz = cbuffer.data(pi_geom_1010_off + 611 * ccomps * dcomps);

            auto g_z_0_y_0_x_yyyzzz = cbuffer.data(pi_geom_1010_off + 612 * ccomps * dcomps);

            auto g_z_0_y_0_x_yyzzzz = cbuffer.data(pi_geom_1010_off + 613 * ccomps * dcomps);

            auto g_z_0_y_0_x_yzzzzz = cbuffer.data(pi_geom_1010_off + 614 * ccomps * dcomps);

            auto g_z_0_y_0_x_zzzzzz = cbuffer.data(pi_geom_1010_off + 615 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_0_xxxxxx, g_z_0_y_0_0_xxxxxxx, g_z_0_y_0_0_xxxxxxy, g_z_0_y_0_0_xxxxxxz, g_z_0_y_0_0_xxxxxy, g_z_0_y_0_0_xxxxxyy, g_z_0_y_0_0_xxxxxyz, g_z_0_y_0_0_xxxxxz, g_z_0_y_0_0_xxxxxzz, g_z_0_y_0_0_xxxxyy, g_z_0_y_0_0_xxxxyyy, g_z_0_y_0_0_xxxxyyz, g_z_0_y_0_0_xxxxyz, g_z_0_y_0_0_xxxxyzz, g_z_0_y_0_0_xxxxzz, g_z_0_y_0_0_xxxxzzz, g_z_0_y_0_0_xxxyyy, g_z_0_y_0_0_xxxyyyy, g_z_0_y_0_0_xxxyyyz, g_z_0_y_0_0_xxxyyz, g_z_0_y_0_0_xxxyyzz, g_z_0_y_0_0_xxxyzz, g_z_0_y_0_0_xxxyzzz, g_z_0_y_0_0_xxxzzz, g_z_0_y_0_0_xxxzzzz, g_z_0_y_0_0_xxyyyy, g_z_0_y_0_0_xxyyyyy, g_z_0_y_0_0_xxyyyyz, g_z_0_y_0_0_xxyyyz, g_z_0_y_0_0_xxyyyzz, g_z_0_y_0_0_xxyyzz, g_z_0_y_0_0_xxyyzzz, g_z_0_y_0_0_xxyzzz, g_z_0_y_0_0_xxyzzzz, g_z_0_y_0_0_xxzzzz, g_z_0_y_0_0_xxzzzzz, g_z_0_y_0_0_xyyyyy, g_z_0_y_0_0_xyyyyyy, g_z_0_y_0_0_xyyyyyz, g_z_0_y_0_0_xyyyyz, g_z_0_y_0_0_xyyyyzz, g_z_0_y_0_0_xyyyzz, g_z_0_y_0_0_xyyyzzz, g_z_0_y_0_0_xyyzzz, g_z_0_y_0_0_xyyzzzz, g_z_0_y_0_0_xyzzzz, g_z_0_y_0_0_xyzzzzz, g_z_0_y_0_0_xzzzzz, g_z_0_y_0_0_xzzzzzz, g_z_0_y_0_0_yyyyyy, g_z_0_y_0_0_yyyyyz, g_z_0_y_0_0_yyyyzz, g_z_0_y_0_0_yyyzzz, g_z_0_y_0_0_yyzzzz, g_z_0_y_0_0_yzzzzz, g_z_0_y_0_0_zzzzzz, g_z_0_y_0_x_xxxxxx, g_z_0_y_0_x_xxxxxy, g_z_0_y_0_x_xxxxxz, g_z_0_y_0_x_xxxxyy, g_z_0_y_0_x_xxxxyz, g_z_0_y_0_x_xxxxzz, g_z_0_y_0_x_xxxyyy, g_z_0_y_0_x_xxxyyz, g_z_0_y_0_x_xxxyzz, g_z_0_y_0_x_xxxzzz, g_z_0_y_0_x_xxyyyy, g_z_0_y_0_x_xxyyyz, g_z_0_y_0_x_xxyyzz, g_z_0_y_0_x_xxyzzz, g_z_0_y_0_x_xxzzzz, g_z_0_y_0_x_xyyyyy, g_z_0_y_0_x_xyyyyz, g_z_0_y_0_x_xyyyzz, g_z_0_y_0_x_xyyzzz, g_z_0_y_0_x_xyzzzz, g_z_0_y_0_x_xzzzzz, g_z_0_y_0_x_yyyyyy, g_z_0_y_0_x_yyyyyz, g_z_0_y_0_x_yyyyzz, g_z_0_y_0_x_yyyzzz, g_z_0_y_0_x_yyzzzz, g_z_0_y_0_x_yzzzzz, g_z_0_y_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_x_xxxxxx[k] = -g_z_0_y_0_0_xxxxxx[k] * ab_x + g_z_0_y_0_0_xxxxxxx[k];

                g_z_0_y_0_x_xxxxxy[k] = -g_z_0_y_0_0_xxxxxy[k] * ab_x + g_z_0_y_0_0_xxxxxxy[k];

                g_z_0_y_0_x_xxxxxz[k] = -g_z_0_y_0_0_xxxxxz[k] * ab_x + g_z_0_y_0_0_xxxxxxz[k];

                g_z_0_y_0_x_xxxxyy[k] = -g_z_0_y_0_0_xxxxyy[k] * ab_x + g_z_0_y_0_0_xxxxxyy[k];

                g_z_0_y_0_x_xxxxyz[k] = -g_z_0_y_0_0_xxxxyz[k] * ab_x + g_z_0_y_0_0_xxxxxyz[k];

                g_z_0_y_0_x_xxxxzz[k] = -g_z_0_y_0_0_xxxxzz[k] * ab_x + g_z_0_y_0_0_xxxxxzz[k];

                g_z_0_y_0_x_xxxyyy[k] = -g_z_0_y_0_0_xxxyyy[k] * ab_x + g_z_0_y_0_0_xxxxyyy[k];

                g_z_0_y_0_x_xxxyyz[k] = -g_z_0_y_0_0_xxxyyz[k] * ab_x + g_z_0_y_0_0_xxxxyyz[k];

                g_z_0_y_0_x_xxxyzz[k] = -g_z_0_y_0_0_xxxyzz[k] * ab_x + g_z_0_y_0_0_xxxxyzz[k];

                g_z_0_y_0_x_xxxzzz[k] = -g_z_0_y_0_0_xxxzzz[k] * ab_x + g_z_0_y_0_0_xxxxzzz[k];

                g_z_0_y_0_x_xxyyyy[k] = -g_z_0_y_0_0_xxyyyy[k] * ab_x + g_z_0_y_0_0_xxxyyyy[k];

                g_z_0_y_0_x_xxyyyz[k] = -g_z_0_y_0_0_xxyyyz[k] * ab_x + g_z_0_y_0_0_xxxyyyz[k];

                g_z_0_y_0_x_xxyyzz[k] = -g_z_0_y_0_0_xxyyzz[k] * ab_x + g_z_0_y_0_0_xxxyyzz[k];

                g_z_0_y_0_x_xxyzzz[k] = -g_z_0_y_0_0_xxyzzz[k] * ab_x + g_z_0_y_0_0_xxxyzzz[k];

                g_z_0_y_0_x_xxzzzz[k] = -g_z_0_y_0_0_xxzzzz[k] * ab_x + g_z_0_y_0_0_xxxzzzz[k];

                g_z_0_y_0_x_xyyyyy[k] = -g_z_0_y_0_0_xyyyyy[k] * ab_x + g_z_0_y_0_0_xxyyyyy[k];

                g_z_0_y_0_x_xyyyyz[k] = -g_z_0_y_0_0_xyyyyz[k] * ab_x + g_z_0_y_0_0_xxyyyyz[k];

                g_z_0_y_0_x_xyyyzz[k] = -g_z_0_y_0_0_xyyyzz[k] * ab_x + g_z_0_y_0_0_xxyyyzz[k];

                g_z_0_y_0_x_xyyzzz[k] = -g_z_0_y_0_0_xyyzzz[k] * ab_x + g_z_0_y_0_0_xxyyzzz[k];

                g_z_0_y_0_x_xyzzzz[k] = -g_z_0_y_0_0_xyzzzz[k] * ab_x + g_z_0_y_0_0_xxyzzzz[k];

                g_z_0_y_0_x_xzzzzz[k] = -g_z_0_y_0_0_xzzzzz[k] * ab_x + g_z_0_y_0_0_xxzzzzz[k];

                g_z_0_y_0_x_yyyyyy[k] = -g_z_0_y_0_0_yyyyyy[k] * ab_x + g_z_0_y_0_0_xyyyyyy[k];

                g_z_0_y_0_x_yyyyyz[k] = -g_z_0_y_0_0_yyyyyz[k] * ab_x + g_z_0_y_0_0_xyyyyyz[k];

                g_z_0_y_0_x_yyyyzz[k] = -g_z_0_y_0_0_yyyyzz[k] * ab_x + g_z_0_y_0_0_xyyyyzz[k];

                g_z_0_y_0_x_yyyzzz[k] = -g_z_0_y_0_0_yyyzzz[k] * ab_x + g_z_0_y_0_0_xyyyzzz[k];

                g_z_0_y_0_x_yyzzzz[k] = -g_z_0_y_0_0_yyzzzz[k] * ab_x + g_z_0_y_0_0_xyyzzzz[k];

                g_z_0_y_0_x_yzzzzz[k] = -g_z_0_y_0_0_yzzzzz[k] * ab_x + g_z_0_y_0_0_xyzzzzz[k];

                g_z_0_y_0_x_zzzzzz[k] = -g_z_0_y_0_0_zzzzzz[k] * ab_x + g_z_0_y_0_0_xzzzzzz[k];
            }

            /// Set up 616-644 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_y_xxxxxx = cbuffer.data(pi_geom_1010_off + 616 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxxxy = cbuffer.data(pi_geom_1010_off + 617 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxxxz = cbuffer.data(pi_geom_1010_off + 618 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxxyy = cbuffer.data(pi_geom_1010_off + 619 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxxyz = cbuffer.data(pi_geom_1010_off + 620 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxxzz = cbuffer.data(pi_geom_1010_off + 621 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxyyy = cbuffer.data(pi_geom_1010_off + 622 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxyyz = cbuffer.data(pi_geom_1010_off + 623 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxyzz = cbuffer.data(pi_geom_1010_off + 624 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxxzzz = cbuffer.data(pi_geom_1010_off + 625 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxyyyy = cbuffer.data(pi_geom_1010_off + 626 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxyyyz = cbuffer.data(pi_geom_1010_off + 627 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxyyzz = cbuffer.data(pi_geom_1010_off + 628 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxyzzz = cbuffer.data(pi_geom_1010_off + 629 * ccomps * dcomps);

            auto g_z_0_y_0_y_xxzzzz = cbuffer.data(pi_geom_1010_off + 630 * ccomps * dcomps);

            auto g_z_0_y_0_y_xyyyyy = cbuffer.data(pi_geom_1010_off + 631 * ccomps * dcomps);

            auto g_z_0_y_0_y_xyyyyz = cbuffer.data(pi_geom_1010_off + 632 * ccomps * dcomps);

            auto g_z_0_y_0_y_xyyyzz = cbuffer.data(pi_geom_1010_off + 633 * ccomps * dcomps);

            auto g_z_0_y_0_y_xyyzzz = cbuffer.data(pi_geom_1010_off + 634 * ccomps * dcomps);

            auto g_z_0_y_0_y_xyzzzz = cbuffer.data(pi_geom_1010_off + 635 * ccomps * dcomps);

            auto g_z_0_y_0_y_xzzzzz = cbuffer.data(pi_geom_1010_off + 636 * ccomps * dcomps);

            auto g_z_0_y_0_y_yyyyyy = cbuffer.data(pi_geom_1010_off + 637 * ccomps * dcomps);

            auto g_z_0_y_0_y_yyyyyz = cbuffer.data(pi_geom_1010_off + 638 * ccomps * dcomps);

            auto g_z_0_y_0_y_yyyyzz = cbuffer.data(pi_geom_1010_off + 639 * ccomps * dcomps);

            auto g_z_0_y_0_y_yyyzzz = cbuffer.data(pi_geom_1010_off + 640 * ccomps * dcomps);

            auto g_z_0_y_0_y_yyzzzz = cbuffer.data(pi_geom_1010_off + 641 * ccomps * dcomps);

            auto g_z_0_y_0_y_yzzzzz = cbuffer.data(pi_geom_1010_off + 642 * ccomps * dcomps);

            auto g_z_0_y_0_y_zzzzzz = cbuffer.data(pi_geom_1010_off + 643 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_0_0_xxxxxx, g_z_0_y_0_0_xxxxxxy, g_z_0_y_0_0_xxxxxy, g_z_0_y_0_0_xxxxxyy, g_z_0_y_0_0_xxxxxyz, g_z_0_y_0_0_xxxxxz, g_z_0_y_0_0_xxxxyy, g_z_0_y_0_0_xxxxyyy, g_z_0_y_0_0_xxxxyyz, g_z_0_y_0_0_xxxxyz, g_z_0_y_0_0_xxxxyzz, g_z_0_y_0_0_xxxxzz, g_z_0_y_0_0_xxxyyy, g_z_0_y_0_0_xxxyyyy, g_z_0_y_0_0_xxxyyyz, g_z_0_y_0_0_xxxyyz, g_z_0_y_0_0_xxxyyzz, g_z_0_y_0_0_xxxyzz, g_z_0_y_0_0_xxxyzzz, g_z_0_y_0_0_xxxzzz, g_z_0_y_0_0_xxyyyy, g_z_0_y_0_0_xxyyyyy, g_z_0_y_0_0_xxyyyyz, g_z_0_y_0_0_xxyyyz, g_z_0_y_0_0_xxyyyzz, g_z_0_y_0_0_xxyyzz, g_z_0_y_0_0_xxyyzzz, g_z_0_y_0_0_xxyzzz, g_z_0_y_0_0_xxyzzzz, g_z_0_y_0_0_xxzzzz, g_z_0_y_0_0_xyyyyy, g_z_0_y_0_0_xyyyyyy, g_z_0_y_0_0_xyyyyyz, g_z_0_y_0_0_xyyyyz, g_z_0_y_0_0_xyyyyzz, g_z_0_y_0_0_xyyyzz, g_z_0_y_0_0_xyyyzzz, g_z_0_y_0_0_xyyzzz, g_z_0_y_0_0_xyyzzzz, g_z_0_y_0_0_xyzzzz, g_z_0_y_0_0_xyzzzzz, g_z_0_y_0_0_xzzzzz, g_z_0_y_0_0_yyyyyy, g_z_0_y_0_0_yyyyyyy, g_z_0_y_0_0_yyyyyyz, g_z_0_y_0_0_yyyyyz, g_z_0_y_0_0_yyyyyzz, g_z_0_y_0_0_yyyyzz, g_z_0_y_0_0_yyyyzzz, g_z_0_y_0_0_yyyzzz, g_z_0_y_0_0_yyyzzzz, g_z_0_y_0_0_yyzzzz, g_z_0_y_0_0_yyzzzzz, g_z_0_y_0_0_yzzzzz, g_z_0_y_0_0_yzzzzzz, g_z_0_y_0_0_zzzzzz, g_z_0_y_0_y_xxxxxx, g_z_0_y_0_y_xxxxxy, g_z_0_y_0_y_xxxxxz, g_z_0_y_0_y_xxxxyy, g_z_0_y_0_y_xxxxyz, g_z_0_y_0_y_xxxxzz, g_z_0_y_0_y_xxxyyy, g_z_0_y_0_y_xxxyyz, g_z_0_y_0_y_xxxyzz, g_z_0_y_0_y_xxxzzz, g_z_0_y_0_y_xxyyyy, g_z_0_y_0_y_xxyyyz, g_z_0_y_0_y_xxyyzz, g_z_0_y_0_y_xxyzzz, g_z_0_y_0_y_xxzzzz, g_z_0_y_0_y_xyyyyy, g_z_0_y_0_y_xyyyyz, g_z_0_y_0_y_xyyyzz, g_z_0_y_0_y_xyyzzz, g_z_0_y_0_y_xyzzzz, g_z_0_y_0_y_xzzzzz, g_z_0_y_0_y_yyyyyy, g_z_0_y_0_y_yyyyyz, g_z_0_y_0_y_yyyyzz, g_z_0_y_0_y_yyyzzz, g_z_0_y_0_y_yyzzzz, g_z_0_y_0_y_yzzzzz, g_z_0_y_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_y_xxxxxx[k] = -g_z_0_y_0_0_xxxxxx[k] * ab_y + g_z_0_y_0_0_xxxxxxy[k];

                g_z_0_y_0_y_xxxxxy[k] = -g_z_0_y_0_0_xxxxxy[k] * ab_y + g_z_0_y_0_0_xxxxxyy[k];

                g_z_0_y_0_y_xxxxxz[k] = -g_z_0_y_0_0_xxxxxz[k] * ab_y + g_z_0_y_0_0_xxxxxyz[k];

                g_z_0_y_0_y_xxxxyy[k] = -g_z_0_y_0_0_xxxxyy[k] * ab_y + g_z_0_y_0_0_xxxxyyy[k];

                g_z_0_y_0_y_xxxxyz[k] = -g_z_0_y_0_0_xxxxyz[k] * ab_y + g_z_0_y_0_0_xxxxyyz[k];

                g_z_0_y_0_y_xxxxzz[k] = -g_z_0_y_0_0_xxxxzz[k] * ab_y + g_z_0_y_0_0_xxxxyzz[k];

                g_z_0_y_0_y_xxxyyy[k] = -g_z_0_y_0_0_xxxyyy[k] * ab_y + g_z_0_y_0_0_xxxyyyy[k];

                g_z_0_y_0_y_xxxyyz[k] = -g_z_0_y_0_0_xxxyyz[k] * ab_y + g_z_0_y_0_0_xxxyyyz[k];

                g_z_0_y_0_y_xxxyzz[k] = -g_z_0_y_0_0_xxxyzz[k] * ab_y + g_z_0_y_0_0_xxxyyzz[k];

                g_z_0_y_0_y_xxxzzz[k] = -g_z_0_y_0_0_xxxzzz[k] * ab_y + g_z_0_y_0_0_xxxyzzz[k];

                g_z_0_y_0_y_xxyyyy[k] = -g_z_0_y_0_0_xxyyyy[k] * ab_y + g_z_0_y_0_0_xxyyyyy[k];

                g_z_0_y_0_y_xxyyyz[k] = -g_z_0_y_0_0_xxyyyz[k] * ab_y + g_z_0_y_0_0_xxyyyyz[k];

                g_z_0_y_0_y_xxyyzz[k] = -g_z_0_y_0_0_xxyyzz[k] * ab_y + g_z_0_y_0_0_xxyyyzz[k];

                g_z_0_y_0_y_xxyzzz[k] = -g_z_0_y_0_0_xxyzzz[k] * ab_y + g_z_0_y_0_0_xxyyzzz[k];

                g_z_0_y_0_y_xxzzzz[k] = -g_z_0_y_0_0_xxzzzz[k] * ab_y + g_z_0_y_0_0_xxyzzzz[k];

                g_z_0_y_0_y_xyyyyy[k] = -g_z_0_y_0_0_xyyyyy[k] * ab_y + g_z_0_y_0_0_xyyyyyy[k];

                g_z_0_y_0_y_xyyyyz[k] = -g_z_0_y_0_0_xyyyyz[k] * ab_y + g_z_0_y_0_0_xyyyyyz[k];

                g_z_0_y_0_y_xyyyzz[k] = -g_z_0_y_0_0_xyyyzz[k] * ab_y + g_z_0_y_0_0_xyyyyzz[k];

                g_z_0_y_0_y_xyyzzz[k] = -g_z_0_y_0_0_xyyzzz[k] * ab_y + g_z_0_y_0_0_xyyyzzz[k];

                g_z_0_y_0_y_xyzzzz[k] = -g_z_0_y_0_0_xyzzzz[k] * ab_y + g_z_0_y_0_0_xyyzzzz[k];

                g_z_0_y_0_y_xzzzzz[k] = -g_z_0_y_0_0_xzzzzz[k] * ab_y + g_z_0_y_0_0_xyzzzzz[k];

                g_z_0_y_0_y_yyyyyy[k] = -g_z_0_y_0_0_yyyyyy[k] * ab_y + g_z_0_y_0_0_yyyyyyy[k];

                g_z_0_y_0_y_yyyyyz[k] = -g_z_0_y_0_0_yyyyyz[k] * ab_y + g_z_0_y_0_0_yyyyyyz[k];

                g_z_0_y_0_y_yyyyzz[k] = -g_z_0_y_0_0_yyyyzz[k] * ab_y + g_z_0_y_0_0_yyyyyzz[k];

                g_z_0_y_0_y_yyyzzz[k] = -g_z_0_y_0_0_yyyzzz[k] * ab_y + g_z_0_y_0_0_yyyyzzz[k];

                g_z_0_y_0_y_yyzzzz[k] = -g_z_0_y_0_0_yyzzzz[k] * ab_y + g_z_0_y_0_0_yyyzzzz[k];

                g_z_0_y_0_y_yzzzzz[k] = -g_z_0_y_0_0_yzzzzz[k] * ab_y + g_z_0_y_0_0_yyzzzzz[k];

                g_z_0_y_0_y_zzzzzz[k] = -g_z_0_y_0_0_zzzzzz[k] * ab_y + g_z_0_y_0_0_yzzzzzz[k];
            }

            /// Set up 644-672 components of targeted buffer : cbuffer.data(

            auto g_z_0_y_0_z_xxxxxx = cbuffer.data(pi_geom_1010_off + 644 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxxxy = cbuffer.data(pi_geom_1010_off + 645 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxxxz = cbuffer.data(pi_geom_1010_off + 646 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxxyy = cbuffer.data(pi_geom_1010_off + 647 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxxyz = cbuffer.data(pi_geom_1010_off + 648 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxxzz = cbuffer.data(pi_geom_1010_off + 649 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxyyy = cbuffer.data(pi_geom_1010_off + 650 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxyyz = cbuffer.data(pi_geom_1010_off + 651 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxyzz = cbuffer.data(pi_geom_1010_off + 652 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxxzzz = cbuffer.data(pi_geom_1010_off + 653 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxyyyy = cbuffer.data(pi_geom_1010_off + 654 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxyyyz = cbuffer.data(pi_geom_1010_off + 655 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxyyzz = cbuffer.data(pi_geom_1010_off + 656 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxyzzz = cbuffer.data(pi_geom_1010_off + 657 * ccomps * dcomps);

            auto g_z_0_y_0_z_xxzzzz = cbuffer.data(pi_geom_1010_off + 658 * ccomps * dcomps);

            auto g_z_0_y_0_z_xyyyyy = cbuffer.data(pi_geom_1010_off + 659 * ccomps * dcomps);

            auto g_z_0_y_0_z_xyyyyz = cbuffer.data(pi_geom_1010_off + 660 * ccomps * dcomps);

            auto g_z_0_y_0_z_xyyyzz = cbuffer.data(pi_geom_1010_off + 661 * ccomps * dcomps);

            auto g_z_0_y_0_z_xyyzzz = cbuffer.data(pi_geom_1010_off + 662 * ccomps * dcomps);

            auto g_z_0_y_0_z_xyzzzz = cbuffer.data(pi_geom_1010_off + 663 * ccomps * dcomps);

            auto g_z_0_y_0_z_xzzzzz = cbuffer.data(pi_geom_1010_off + 664 * ccomps * dcomps);

            auto g_z_0_y_0_z_yyyyyy = cbuffer.data(pi_geom_1010_off + 665 * ccomps * dcomps);

            auto g_z_0_y_0_z_yyyyyz = cbuffer.data(pi_geom_1010_off + 666 * ccomps * dcomps);

            auto g_z_0_y_0_z_yyyyzz = cbuffer.data(pi_geom_1010_off + 667 * ccomps * dcomps);

            auto g_z_0_y_0_z_yyyzzz = cbuffer.data(pi_geom_1010_off + 668 * ccomps * dcomps);

            auto g_z_0_y_0_z_yyzzzz = cbuffer.data(pi_geom_1010_off + 669 * ccomps * dcomps);

            auto g_z_0_y_0_z_yzzzzz = cbuffer.data(pi_geom_1010_off + 670 * ccomps * dcomps);

            auto g_z_0_y_0_z_zzzzzz = cbuffer.data(pi_geom_1010_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_y_0_0_xxxxxx, g_0_0_y_0_0_xxxxxy, g_0_0_y_0_0_xxxxxz, g_0_0_y_0_0_xxxxyy, g_0_0_y_0_0_xxxxyz, g_0_0_y_0_0_xxxxzz, g_0_0_y_0_0_xxxyyy, g_0_0_y_0_0_xxxyyz, g_0_0_y_0_0_xxxyzz, g_0_0_y_0_0_xxxzzz, g_0_0_y_0_0_xxyyyy, g_0_0_y_0_0_xxyyyz, g_0_0_y_0_0_xxyyzz, g_0_0_y_0_0_xxyzzz, g_0_0_y_0_0_xxzzzz, g_0_0_y_0_0_xyyyyy, g_0_0_y_0_0_xyyyyz, g_0_0_y_0_0_xyyyzz, g_0_0_y_0_0_xyyzzz, g_0_0_y_0_0_xyzzzz, g_0_0_y_0_0_xzzzzz, g_0_0_y_0_0_yyyyyy, g_0_0_y_0_0_yyyyyz, g_0_0_y_0_0_yyyyzz, g_0_0_y_0_0_yyyzzz, g_0_0_y_0_0_yyzzzz, g_0_0_y_0_0_yzzzzz, g_0_0_y_0_0_zzzzzz, g_z_0_y_0_0_xxxxxx, g_z_0_y_0_0_xxxxxxz, g_z_0_y_0_0_xxxxxy, g_z_0_y_0_0_xxxxxyz, g_z_0_y_0_0_xxxxxz, g_z_0_y_0_0_xxxxxzz, g_z_0_y_0_0_xxxxyy, g_z_0_y_0_0_xxxxyyz, g_z_0_y_0_0_xxxxyz, g_z_0_y_0_0_xxxxyzz, g_z_0_y_0_0_xxxxzz, g_z_0_y_0_0_xxxxzzz, g_z_0_y_0_0_xxxyyy, g_z_0_y_0_0_xxxyyyz, g_z_0_y_0_0_xxxyyz, g_z_0_y_0_0_xxxyyzz, g_z_0_y_0_0_xxxyzz, g_z_0_y_0_0_xxxyzzz, g_z_0_y_0_0_xxxzzz, g_z_0_y_0_0_xxxzzzz, g_z_0_y_0_0_xxyyyy, g_z_0_y_0_0_xxyyyyz, g_z_0_y_0_0_xxyyyz, g_z_0_y_0_0_xxyyyzz, g_z_0_y_0_0_xxyyzz, g_z_0_y_0_0_xxyyzzz, g_z_0_y_0_0_xxyzzz, g_z_0_y_0_0_xxyzzzz, g_z_0_y_0_0_xxzzzz, g_z_0_y_0_0_xxzzzzz, g_z_0_y_0_0_xyyyyy, g_z_0_y_0_0_xyyyyyz, g_z_0_y_0_0_xyyyyz, g_z_0_y_0_0_xyyyyzz, g_z_0_y_0_0_xyyyzz, g_z_0_y_0_0_xyyyzzz, g_z_0_y_0_0_xyyzzz, g_z_0_y_0_0_xyyzzzz, g_z_0_y_0_0_xyzzzz, g_z_0_y_0_0_xyzzzzz, g_z_0_y_0_0_xzzzzz, g_z_0_y_0_0_xzzzzzz, g_z_0_y_0_0_yyyyyy, g_z_0_y_0_0_yyyyyyz, g_z_0_y_0_0_yyyyyz, g_z_0_y_0_0_yyyyyzz, g_z_0_y_0_0_yyyyzz, g_z_0_y_0_0_yyyyzzz, g_z_0_y_0_0_yyyzzz, g_z_0_y_0_0_yyyzzzz, g_z_0_y_0_0_yyzzzz, g_z_0_y_0_0_yyzzzzz, g_z_0_y_0_0_yzzzzz, g_z_0_y_0_0_yzzzzzz, g_z_0_y_0_0_zzzzzz, g_z_0_y_0_0_zzzzzzz, g_z_0_y_0_z_xxxxxx, g_z_0_y_0_z_xxxxxy, g_z_0_y_0_z_xxxxxz, g_z_0_y_0_z_xxxxyy, g_z_0_y_0_z_xxxxyz, g_z_0_y_0_z_xxxxzz, g_z_0_y_0_z_xxxyyy, g_z_0_y_0_z_xxxyyz, g_z_0_y_0_z_xxxyzz, g_z_0_y_0_z_xxxzzz, g_z_0_y_0_z_xxyyyy, g_z_0_y_0_z_xxyyyz, g_z_0_y_0_z_xxyyzz, g_z_0_y_0_z_xxyzzz, g_z_0_y_0_z_xxzzzz, g_z_0_y_0_z_xyyyyy, g_z_0_y_0_z_xyyyyz, g_z_0_y_0_z_xyyyzz, g_z_0_y_0_z_xyyzzz, g_z_0_y_0_z_xyzzzz, g_z_0_y_0_z_xzzzzz, g_z_0_y_0_z_yyyyyy, g_z_0_y_0_z_yyyyyz, g_z_0_y_0_z_yyyyzz, g_z_0_y_0_z_yyyzzz, g_z_0_y_0_z_yyzzzz, g_z_0_y_0_z_yzzzzz, g_z_0_y_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_y_0_z_xxxxxx[k] = -g_0_0_y_0_0_xxxxxx[k] - g_z_0_y_0_0_xxxxxx[k] * ab_z + g_z_0_y_0_0_xxxxxxz[k];

                g_z_0_y_0_z_xxxxxy[k] = -g_0_0_y_0_0_xxxxxy[k] - g_z_0_y_0_0_xxxxxy[k] * ab_z + g_z_0_y_0_0_xxxxxyz[k];

                g_z_0_y_0_z_xxxxxz[k] = -g_0_0_y_0_0_xxxxxz[k] - g_z_0_y_0_0_xxxxxz[k] * ab_z + g_z_0_y_0_0_xxxxxzz[k];

                g_z_0_y_0_z_xxxxyy[k] = -g_0_0_y_0_0_xxxxyy[k] - g_z_0_y_0_0_xxxxyy[k] * ab_z + g_z_0_y_0_0_xxxxyyz[k];

                g_z_0_y_0_z_xxxxyz[k] = -g_0_0_y_0_0_xxxxyz[k] - g_z_0_y_0_0_xxxxyz[k] * ab_z + g_z_0_y_0_0_xxxxyzz[k];

                g_z_0_y_0_z_xxxxzz[k] = -g_0_0_y_0_0_xxxxzz[k] - g_z_0_y_0_0_xxxxzz[k] * ab_z + g_z_0_y_0_0_xxxxzzz[k];

                g_z_0_y_0_z_xxxyyy[k] = -g_0_0_y_0_0_xxxyyy[k] - g_z_0_y_0_0_xxxyyy[k] * ab_z + g_z_0_y_0_0_xxxyyyz[k];

                g_z_0_y_0_z_xxxyyz[k] = -g_0_0_y_0_0_xxxyyz[k] - g_z_0_y_0_0_xxxyyz[k] * ab_z + g_z_0_y_0_0_xxxyyzz[k];

                g_z_0_y_0_z_xxxyzz[k] = -g_0_0_y_0_0_xxxyzz[k] - g_z_0_y_0_0_xxxyzz[k] * ab_z + g_z_0_y_0_0_xxxyzzz[k];

                g_z_0_y_0_z_xxxzzz[k] = -g_0_0_y_0_0_xxxzzz[k] - g_z_0_y_0_0_xxxzzz[k] * ab_z + g_z_0_y_0_0_xxxzzzz[k];

                g_z_0_y_0_z_xxyyyy[k] = -g_0_0_y_0_0_xxyyyy[k] - g_z_0_y_0_0_xxyyyy[k] * ab_z + g_z_0_y_0_0_xxyyyyz[k];

                g_z_0_y_0_z_xxyyyz[k] = -g_0_0_y_0_0_xxyyyz[k] - g_z_0_y_0_0_xxyyyz[k] * ab_z + g_z_0_y_0_0_xxyyyzz[k];

                g_z_0_y_0_z_xxyyzz[k] = -g_0_0_y_0_0_xxyyzz[k] - g_z_0_y_0_0_xxyyzz[k] * ab_z + g_z_0_y_0_0_xxyyzzz[k];

                g_z_0_y_0_z_xxyzzz[k] = -g_0_0_y_0_0_xxyzzz[k] - g_z_0_y_0_0_xxyzzz[k] * ab_z + g_z_0_y_0_0_xxyzzzz[k];

                g_z_0_y_0_z_xxzzzz[k] = -g_0_0_y_0_0_xxzzzz[k] - g_z_0_y_0_0_xxzzzz[k] * ab_z + g_z_0_y_0_0_xxzzzzz[k];

                g_z_0_y_0_z_xyyyyy[k] = -g_0_0_y_0_0_xyyyyy[k] - g_z_0_y_0_0_xyyyyy[k] * ab_z + g_z_0_y_0_0_xyyyyyz[k];

                g_z_0_y_0_z_xyyyyz[k] = -g_0_0_y_0_0_xyyyyz[k] - g_z_0_y_0_0_xyyyyz[k] * ab_z + g_z_0_y_0_0_xyyyyzz[k];

                g_z_0_y_0_z_xyyyzz[k] = -g_0_0_y_0_0_xyyyzz[k] - g_z_0_y_0_0_xyyyzz[k] * ab_z + g_z_0_y_0_0_xyyyzzz[k];

                g_z_0_y_0_z_xyyzzz[k] = -g_0_0_y_0_0_xyyzzz[k] - g_z_0_y_0_0_xyyzzz[k] * ab_z + g_z_0_y_0_0_xyyzzzz[k];

                g_z_0_y_0_z_xyzzzz[k] = -g_0_0_y_0_0_xyzzzz[k] - g_z_0_y_0_0_xyzzzz[k] * ab_z + g_z_0_y_0_0_xyzzzzz[k];

                g_z_0_y_0_z_xzzzzz[k] = -g_0_0_y_0_0_xzzzzz[k] - g_z_0_y_0_0_xzzzzz[k] * ab_z + g_z_0_y_0_0_xzzzzzz[k];

                g_z_0_y_0_z_yyyyyy[k] = -g_0_0_y_0_0_yyyyyy[k] - g_z_0_y_0_0_yyyyyy[k] * ab_z + g_z_0_y_0_0_yyyyyyz[k];

                g_z_0_y_0_z_yyyyyz[k] = -g_0_0_y_0_0_yyyyyz[k] - g_z_0_y_0_0_yyyyyz[k] * ab_z + g_z_0_y_0_0_yyyyyzz[k];

                g_z_0_y_0_z_yyyyzz[k] = -g_0_0_y_0_0_yyyyzz[k] - g_z_0_y_0_0_yyyyzz[k] * ab_z + g_z_0_y_0_0_yyyyzzz[k];

                g_z_0_y_0_z_yyyzzz[k] = -g_0_0_y_0_0_yyyzzz[k] - g_z_0_y_0_0_yyyzzz[k] * ab_z + g_z_0_y_0_0_yyyzzzz[k];

                g_z_0_y_0_z_yyzzzz[k] = -g_0_0_y_0_0_yyzzzz[k] - g_z_0_y_0_0_yyzzzz[k] * ab_z + g_z_0_y_0_0_yyzzzzz[k];

                g_z_0_y_0_z_yzzzzz[k] = -g_0_0_y_0_0_yzzzzz[k] - g_z_0_y_0_0_yzzzzz[k] * ab_z + g_z_0_y_0_0_yzzzzzz[k];

                g_z_0_y_0_z_zzzzzz[k] = -g_0_0_y_0_0_zzzzzz[k] - g_z_0_y_0_0_zzzzzz[k] * ab_z + g_z_0_y_0_0_zzzzzzz[k];
            }

            /// Set up 672-700 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_x_xxxxxx = cbuffer.data(pi_geom_1010_off + 672 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxxxy = cbuffer.data(pi_geom_1010_off + 673 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxxxz = cbuffer.data(pi_geom_1010_off + 674 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxxyy = cbuffer.data(pi_geom_1010_off + 675 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxxyz = cbuffer.data(pi_geom_1010_off + 676 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxxzz = cbuffer.data(pi_geom_1010_off + 677 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxyyy = cbuffer.data(pi_geom_1010_off + 678 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxyyz = cbuffer.data(pi_geom_1010_off + 679 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxyzz = cbuffer.data(pi_geom_1010_off + 680 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxxzzz = cbuffer.data(pi_geom_1010_off + 681 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxyyyy = cbuffer.data(pi_geom_1010_off + 682 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxyyyz = cbuffer.data(pi_geom_1010_off + 683 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxyyzz = cbuffer.data(pi_geom_1010_off + 684 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxyzzz = cbuffer.data(pi_geom_1010_off + 685 * ccomps * dcomps);

            auto g_z_0_z_0_x_xxzzzz = cbuffer.data(pi_geom_1010_off + 686 * ccomps * dcomps);

            auto g_z_0_z_0_x_xyyyyy = cbuffer.data(pi_geom_1010_off + 687 * ccomps * dcomps);

            auto g_z_0_z_0_x_xyyyyz = cbuffer.data(pi_geom_1010_off + 688 * ccomps * dcomps);

            auto g_z_0_z_0_x_xyyyzz = cbuffer.data(pi_geom_1010_off + 689 * ccomps * dcomps);

            auto g_z_0_z_0_x_xyyzzz = cbuffer.data(pi_geom_1010_off + 690 * ccomps * dcomps);

            auto g_z_0_z_0_x_xyzzzz = cbuffer.data(pi_geom_1010_off + 691 * ccomps * dcomps);

            auto g_z_0_z_0_x_xzzzzz = cbuffer.data(pi_geom_1010_off + 692 * ccomps * dcomps);

            auto g_z_0_z_0_x_yyyyyy = cbuffer.data(pi_geom_1010_off + 693 * ccomps * dcomps);

            auto g_z_0_z_0_x_yyyyyz = cbuffer.data(pi_geom_1010_off + 694 * ccomps * dcomps);

            auto g_z_0_z_0_x_yyyyzz = cbuffer.data(pi_geom_1010_off + 695 * ccomps * dcomps);

            auto g_z_0_z_0_x_yyyzzz = cbuffer.data(pi_geom_1010_off + 696 * ccomps * dcomps);

            auto g_z_0_z_0_x_yyzzzz = cbuffer.data(pi_geom_1010_off + 697 * ccomps * dcomps);

            auto g_z_0_z_0_x_yzzzzz = cbuffer.data(pi_geom_1010_off + 698 * ccomps * dcomps);

            auto g_z_0_z_0_x_zzzzzz = cbuffer.data(pi_geom_1010_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_0_xxxxxx, g_z_0_z_0_0_xxxxxxx, g_z_0_z_0_0_xxxxxxy, g_z_0_z_0_0_xxxxxxz, g_z_0_z_0_0_xxxxxy, g_z_0_z_0_0_xxxxxyy, g_z_0_z_0_0_xxxxxyz, g_z_0_z_0_0_xxxxxz, g_z_0_z_0_0_xxxxxzz, g_z_0_z_0_0_xxxxyy, g_z_0_z_0_0_xxxxyyy, g_z_0_z_0_0_xxxxyyz, g_z_0_z_0_0_xxxxyz, g_z_0_z_0_0_xxxxyzz, g_z_0_z_0_0_xxxxzz, g_z_0_z_0_0_xxxxzzz, g_z_0_z_0_0_xxxyyy, g_z_0_z_0_0_xxxyyyy, g_z_0_z_0_0_xxxyyyz, g_z_0_z_0_0_xxxyyz, g_z_0_z_0_0_xxxyyzz, g_z_0_z_0_0_xxxyzz, g_z_0_z_0_0_xxxyzzz, g_z_0_z_0_0_xxxzzz, g_z_0_z_0_0_xxxzzzz, g_z_0_z_0_0_xxyyyy, g_z_0_z_0_0_xxyyyyy, g_z_0_z_0_0_xxyyyyz, g_z_0_z_0_0_xxyyyz, g_z_0_z_0_0_xxyyyzz, g_z_0_z_0_0_xxyyzz, g_z_0_z_0_0_xxyyzzz, g_z_0_z_0_0_xxyzzz, g_z_0_z_0_0_xxyzzzz, g_z_0_z_0_0_xxzzzz, g_z_0_z_0_0_xxzzzzz, g_z_0_z_0_0_xyyyyy, g_z_0_z_0_0_xyyyyyy, g_z_0_z_0_0_xyyyyyz, g_z_0_z_0_0_xyyyyz, g_z_0_z_0_0_xyyyyzz, g_z_0_z_0_0_xyyyzz, g_z_0_z_0_0_xyyyzzz, g_z_0_z_0_0_xyyzzz, g_z_0_z_0_0_xyyzzzz, g_z_0_z_0_0_xyzzzz, g_z_0_z_0_0_xyzzzzz, g_z_0_z_0_0_xzzzzz, g_z_0_z_0_0_xzzzzzz, g_z_0_z_0_0_yyyyyy, g_z_0_z_0_0_yyyyyz, g_z_0_z_0_0_yyyyzz, g_z_0_z_0_0_yyyzzz, g_z_0_z_0_0_yyzzzz, g_z_0_z_0_0_yzzzzz, g_z_0_z_0_0_zzzzzz, g_z_0_z_0_x_xxxxxx, g_z_0_z_0_x_xxxxxy, g_z_0_z_0_x_xxxxxz, g_z_0_z_0_x_xxxxyy, g_z_0_z_0_x_xxxxyz, g_z_0_z_0_x_xxxxzz, g_z_0_z_0_x_xxxyyy, g_z_0_z_0_x_xxxyyz, g_z_0_z_0_x_xxxyzz, g_z_0_z_0_x_xxxzzz, g_z_0_z_0_x_xxyyyy, g_z_0_z_0_x_xxyyyz, g_z_0_z_0_x_xxyyzz, g_z_0_z_0_x_xxyzzz, g_z_0_z_0_x_xxzzzz, g_z_0_z_0_x_xyyyyy, g_z_0_z_0_x_xyyyyz, g_z_0_z_0_x_xyyyzz, g_z_0_z_0_x_xyyzzz, g_z_0_z_0_x_xyzzzz, g_z_0_z_0_x_xzzzzz, g_z_0_z_0_x_yyyyyy, g_z_0_z_0_x_yyyyyz, g_z_0_z_0_x_yyyyzz, g_z_0_z_0_x_yyyzzz, g_z_0_z_0_x_yyzzzz, g_z_0_z_0_x_yzzzzz, g_z_0_z_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_x_xxxxxx[k] = -g_z_0_z_0_0_xxxxxx[k] * ab_x + g_z_0_z_0_0_xxxxxxx[k];

                g_z_0_z_0_x_xxxxxy[k] = -g_z_0_z_0_0_xxxxxy[k] * ab_x + g_z_0_z_0_0_xxxxxxy[k];

                g_z_0_z_0_x_xxxxxz[k] = -g_z_0_z_0_0_xxxxxz[k] * ab_x + g_z_0_z_0_0_xxxxxxz[k];

                g_z_0_z_0_x_xxxxyy[k] = -g_z_0_z_0_0_xxxxyy[k] * ab_x + g_z_0_z_0_0_xxxxxyy[k];

                g_z_0_z_0_x_xxxxyz[k] = -g_z_0_z_0_0_xxxxyz[k] * ab_x + g_z_0_z_0_0_xxxxxyz[k];

                g_z_0_z_0_x_xxxxzz[k] = -g_z_0_z_0_0_xxxxzz[k] * ab_x + g_z_0_z_0_0_xxxxxzz[k];

                g_z_0_z_0_x_xxxyyy[k] = -g_z_0_z_0_0_xxxyyy[k] * ab_x + g_z_0_z_0_0_xxxxyyy[k];

                g_z_0_z_0_x_xxxyyz[k] = -g_z_0_z_0_0_xxxyyz[k] * ab_x + g_z_0_z_0_0_xxxxyyz[k];

                g_z_0_z_0_x_xxxyzz[k] = -g_z_0_z_0_0_xxxyzz[k] * ab_x + g_z_0_z_0_0_xxxxyzz[k];

                g_z_0_z_0_x_xxxzzz[k] = -g_z_0_z_0_0_xxxzzz[k] * ab_x + g_z_0_z_0_0_xxxxzzz[k];

                g_z_0_z_0_x_xxyyyy[k] = -g_z_0_z_0_0_xxyyyy[k] * ab_x + g_z_0_z_0_0_xxxyyyy[k];

                g_z_0_z_0_x_xxyyyz[k] = -g_z_0_z_0_0_xxyyyz[k] * ab_x + g_z_0_z_0_0_xxxyyyz[k];

                g_z_0_z_0_x_xxyyzz[k] = -g_z_0_z_0_0_xxyyzz[k] * ab_x + g_z_0_z_0_0_xxxyyzz[k];

                g_z_0_z_0_x_xxyzzz[k] = -g_z_0_z_0_0_xxyzzz[k] * ab_x + g_z_0_z_0_0_xxxyzzz[k];

                g_z_0_z_0_x_xxzzzz[k] = -g_z_0_z_0_0_xxzzzz[k] * ab_x + g_z_0_z_0_0_xxxzzzz[k];

                g_z_0_z_0_x_xyyyyy[k] = -g_z_0_z_0_0_xyyyyy[k] * ab_x + g_z_0_z_0_0_xxyyyyy[k];

                g_z_0_z_0_x_xyyyyz[k] = -g_z_0_z_0_0_xyyyyz[k] * ab_x + g_z_0_z_0_0_xxyyyyz[k];

                g_z_0_z_0_x_xyyyzz[k] = -g_z_0_z_0_0_xyyyzz[k] * ab_x + g_z_0_z_0_0_xxyyyzz[k];

                g_z_0_z_0_x_xyyzzz[k] = -g_z_0_z_0_0_xyyzzz[k] * ab_x + g_z_0_z_0_0_xxyyzzz[k];

                g_z_0_z_0_x_xyzzzz[k] = -g_z_0_z_0_0_xyzzzz[k] * ab_x + g_z_0_z_0_0_xxyzzzz[k];

                g_z_0_z_0_x_xzzzzz[k] = -g_z_0_z_0_0_xzzzzz[k] * ab_x + g_z_0_z_0_0_xxzzzzz[k];

                g_z_0_z_0_x_yyyyyy[k] = -g_z_0_z_0_0_yyyyyy[k] * ab_x + g_z_0_z_0_0_xyyyyyy[k];

                g_z_0_z_0_x_yyyyyz[k] = -g_z_0_z_0_0_yyyyyz[k] * ab_x + g_z_0_z_0_0_xyyyyyz[k];

                g_z_0_z_0_x_yyyyzz[k] = -g_z_0_z_0_0_yyyyzz[k] * ab_x + g_z_0_z_0_0_xyyyyzz[k];

                g_z_0_z_0_x_yyyzzz[k] = -g_z_0_z_0_0_yyyzzz[k] * ab_x + g_z_0_z_0_0_xyyyzzz[k];

                g_z_0_z_0_x_yyzzzz[k] = -g_z_0_z_0_0_yyzzzz[k] * ab_x + g_z_0_z_0_0_xyyzzzz[k];

                g_z_0_z_0_x_yzzzzz[k] = -g_z_0_z_0_0_yzzzzz[k] * ab_x + g_z_0_z_0_0_xyzzzzz[k];

                g_z_0_z_0_x_zzzzzz[k] = -g_z_0_z_0_0_zzzzzz[k] * ab_x + g_z_0_z_0_0_xzzzzzz[k];
            }

            /// Set up 700-728 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_y_xxxxxx = cbuffer.data(pi_geom_1010_off + 700 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxxxy = cbuffer.data(pi_geom_1010_off + 701 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxxxz = cbuffer.data(pi_geom_1010_off + 702 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxxyy = cbuffer.data(pi_geom_1010_off + 703 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxxyz = cbuffer.data(pi_geom_1010_off + 704 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxxzz = cbuffer.data(pi_geom_1010_off + 705 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxyyy = cbuffer.data(pi_geom_1010_off + 706 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxyyz = cbuffer.data(pi_geom_1010_off + 707 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxyzz = cbuffer.data(pi_geom_1010_off + 708 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxxzzz = cbuffer.data(pi_geom_1010_off + 709 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxyyyy = cbuffer.data(pi_geom_1010_off + 710 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxyyyz = cbuffer.data(pi_geom_1010_off + 711 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxyyzz = cbuffer.data(pi_geom_1010_off + 712 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxyzzz = cbuffer.data(pi_geom_1010_off + 713 * ccomps * dcomps);

            auto g_z_0_z_0_y_xxzzzz = cbuffer.data(pi_geom_1010_off + 714 * ccomps * dcomps);

            auto g_z_0_z_0_y_xyyyyy = cbuffer.data(pi_geom_1010_off + 715 * ccomps * dcomps);

            auto g_z_0_z_0_y_xyyyyz = cbuffer.data(pi_geom_1010_off + 716 * ccomps * dcomps);

            auto g_z_0_z_0_y_xyyyzz = cbuffer.data(pi_geom_1010_off + 717 * ccomps * dcomps);

            auto g_z_0_z_0_y_xyyzzz = cbuffer.data(pi_geom_1010_off + 718 * ccomps * dcomps);

            auto g_z_0_z_0_y_xyzzzz = cbuffer.data(pi_geom_1010_off + 719 * ccomps * dcomps);

            auto g_z_0_z_0_y_xzzzzz = cbuffer.data(pi_geom_1010_off + 720 * ccomps * dcomps);

            auto g_z_0_z_0_y_yyyyyy = cbuffer.data(pi_geom_1010_off + 721 * ccomps * dcomps);

            auto g_z_0_z_0_y_yyyyyz = cbuffer.data(pi_geom_1010_off + 722 * ccomps * dcomps);

            auto g_z_0_z_0_y_yyyyzz = cbuffer.data(pi_geom_1010_off + 723 * ccomps * dcomps);

            auto g_z_0_z_0_y_yyyzzz = cbuffer.data(pi_geom_1010_off + 724 * ccomps * dcomps);

            auto g_z_0_z_0_y_yyzzzz = cbuffer.data(pi_geom_1010_off + 725 * ccomps * dcomps);

            auto g_z_0_z_0_y_yzzzzz = cbuffer.data(pi_geom_1010_off + 726 * ccomps * dcomps);

            auto g_z_0_z_0_y_zzzzzz = cbuffer.data(pi_geom_1010_off + 727 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_0_0_xxxxxx, g_z_0_z_0_0_xxxxxxy, g_z_0_z_0_0_xxxxxy, g_z_0_z_0_0_xxxxxyy, g_z_0_z_0_0_xxxxxyz, g_z_0_z_0_0_xxxxxz, g_z_0_z_0_0_xxxxyy, g_z_0_z_0_0_xxxxyyy, g_z_0_z_0_0_xxxxyyz, g_z_0_z_0_0_xxxxyz, g_z_0_z_0_0_xxxxyzz, g_z_0_z_0_0_xxxxzz, g_z_0_z_0_0_xxxyyy, g_z_0_z_0_0_xxxyyyy, g_z_0_z_0_0_xxxyyyz, g_z_0_z_0_0_xxxyyz, g_z_0_z_0_0_xxxyyzz, g_z_0_z_0_0_xxxyzz, g_z_0_z_0_0_xxxyzzz, g_z_0_z_0_0_xxxzzz, g_z_0_z_0_0_xxyyyy, g_z_0_z_0_0_xxyyyyy, g_z_0_z_0_0_xxyyyyz, g_z_0_z_0_0_xxyyyz, g_z_0_z_0_0_xxyyyzz, g_z_0_z_0_0_xxyyzz, g_z_0_z_0_0_xxyyzzz, g_z_0_z_0_0_xxyzzz, g_z_0_z_0_0_xxyzzzz, g_z_0_z_0_0_xxzzzz, g_z_0_z_0_0_xyyyyy, g_z_0_z_0_0_xyyyyyy, g_z_0_z_0_0_xyyyyyz, g_z_0_z_0_0_xyyyyz, g_z_0_z_0_0_xyyyyzz, g_z_0_z_0_0_xyyyzz, g_z_0_z_0_0_xyyyzzz, g_z_0_z_0_0_xyyzzz, g_z_0_z_0_0_xyyzzzz, g_z_0_z_0_0_xyzzzz, g_z_0_z_0_0_xyzzzzz, g_z_0_z_0_0_xzzzzz, g_z_0_z_0_0_yyyyyy, g_z_0_z_0_0_yyyyyyy, g_z_0_z_0_0_yyyyyyz, g_z_0_z_0_0_yyyyyz, g_z_0_z_0_0_yyyyyzz, g_z_0_z_0_0_yyyyzz, g_z_0_z_0_0_yyyyzzz, g_z_0_z_0_0_yyyzzz, g_z_0_z_0_0_yyyzzzz, g_z_0_z_0_0_yyzzzz, g_z_0_z_0_0_yyzzzzz, g_z_0_z_0_0_yzzzzz, g_z_0_z_0_0_yzzzzzz, g_z_0_z_0_0_zzzzzz, g_z_0_z_0_y_xxxxxx, g_z_0_z_0_y_xxxxxy, g_z_0_z_0_y_xxxxxz, g_z_0_z_0_y_xxxxyy, g_z_0_z_0_y_xxxxyz, g_z_0_z_0_y_xxxxzz, g_z_0_z_0_y_xxxyyy, g_z_0_z_0_y_xxxyyz, g_z_0_z_0_y_xxxyzz, g_z_0_z_0_y_xxxzzz, g_z_0_z_0_y_xxyyyy, g_z_0_z_0_y_xxyyyz, g_z_0_z_0_y_xxyyzz, g_z_0_z_0_y_xxyzzz, g_z_0_z_0_y_xxzzzz, g_z_0_z_0_y_xyyyyy, g_z_0_z_0_y_xyyyyz, g_z_0_z_0_y_xyyyzz, g_z_0_z_0_y_xyyzzz, g_z_0_z_0_y_xyzzzz, g_z_0_z_0_y_xzzzzz, g_z_0_z_0_y_yyyyyy, g_z_0_z_0_y_yyyyyz, g_z_0_z_0_y_yyyyzz, g_z_0_z_0_y_yyyzzz, g_z_0_z_0_y_yyzzzz, g_z_0_z_0_y_yzzzzz, g_z_0_z_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_y_xxxxxx[k] = -g_z_0_z_0_0_xxxxxx[k] * ab_y + g_z_0_z_0_0_xxxxxxy[k];

                g_z_0_z_0_y_xxxxxy[k] = -g_z_0_z_0_0_xxxxxy[k] * ab_y + g_z_0_z_0_0_xxxxxyy[k];

                g_z_0_z_0_y_xxxxxz[k] = -g_z_0_z_0_0_xxxxxz[k] * ab_y + g_z_0_z_0_0_xxxxxyz[k];

                g_z_0_z_0_y_xxxxyy[k] = -g_z_0_z_0_0_xxxxyy[k] * ab_y + g_z_0_z_0_0_xxxxyyy[k];

                g_z_0_z_0_y_xxxxyz[k] = -g_z_0_z_0_0_xxxxyz[k] * ab_y + g_z_0_z_0_0_xxxxyyz[k];

                g_z_0_z_0_y_xxxxzz[k] = -g_z_0_z_0_0_xxxxzz[k] * ab_y + g_z_0_z_0_0_xxxxyzz[k];

                g_z_0_z_0_y_xxxyyy[k] = -g_z_0_z_0_0_xxxyyy[k] * ab_y + g_z_0_z_0_0_xxxyyyy[k];

                g_z_0_z_0_y_xxxyyz[k] = -g_z_0_z_0_0_xxxyyz[k] * ab_y + g_z_0_z_0_0_xxxyyyz[k];

                g_z_0_z_0_y_xxxyzz[k] = -g_z_0_z_0_0_xxxyzz[k] * ab_y + g_z_0_z_0_0_xxxyyzz[k];

                g_z_0_z_0_y_xxxzzz[k] = -g_z_0_z_0_0_xxxzzz[k] * ab_y + g_z_0_z_0_0_xxxyzzz[k];

                g_z_0_z_0_y_xxyyyy[k] = -g_z_0_z_0_0_xxyyyy[k] * ab_y + g_z_0_z_0_0_xxyyyyy[k];

                g_z_0_z_0_y_xxyyyz[k] = -g_z_0_z_0_0_xxyyyz[k] * ab_y + g_z_0_z_0_0_xxyyyyz[k];

                g_z_0_z_0_y_xxyyzz[k] = -g_z_0_z_0_0_xxyyzz[k] * ab_y + g_z_0_z_0_0_xxyyyzz[k];

                g_z_0_z_0_y_xxyzzz[k] = -g_z_0_z_0_0_xxyzzz[k] * ab_y + g_z_0_z_0_0_xxyyzzz[k];

                g_z_0_z_0_y_xxzzzz[k] = -g_z_0_z_0_0_xxzzzz[k] * ab_y + g_z_0_z_0_0_xxyzzzz[k];

                g_z_0_z_0_y_xyyyyy[k] = -g_z_0_z_0_0_xyyyyy[k] * ab_y + g_z_0_z_0_0_xyyyyyy[k];

                g_z_0_z_0_y_xyyyyz[k] = -g_z_0_z_0_0_xyyyyz[k] * ab_y + g_z_0_z_0_0_xyyyyyz[k];

                g_z_0_z_0_y_xyyyzz[k] = -g_z_0_z_0_0_xyyyzz[k] * ab_y + g_z_0_z_0_0_xyyyyzz[k];

                g_z_0_z_0_y_xyyzzz[k] = -g_z_0_z_0_0_xyyzzz[k] * ab_y + g_z_0_z_0_0_xyyyzzz[k];

                g_z_0_z_0_y_xyzzzz[k] = -g_z_0_z_0_0_xyzzzz[k] * ab_y + g_z_0_z_0_0_xyyzzzz[k];

                g_z_0_z_0_y_xzzzzz[k] = -g_z_0_z_0_0_xzzzzz[k] * ab_y + g_z_0_z_0_0_xyzzzzz[k];

                g_z_0_z_0_y_yyyyyy[k] = -g_z_0_z_0_0_yyyyyy[k] * ab_y + g_z_0_z_0_0_yyyyyyy[k];

                g_z_0_z_0_y_yyyyyz[k] = -g_z_0_z_0_0_yyyyyz[k] * ab_y + g_z_0_z_0_0_yyyyyyz[k];

                g_z_0_z_0_y_yyyyzz[k] = -g_z_0_z_0_0_yyyyzz[k] * ab_y + g_z_0_z_0_0_yyyyyzz[k];

                g_z_0_z_0_y_yyyzzz[k] = -g_z_0_z_0_0_yyyzzz[k] * ab_y + g_z_0_z_0_0_yyyyzzz[k];

                g_z_0_z_0_y_yyzzzz[k] = -g_z_0_z_0_0_yyzzzz[k] * ab_y + g_z_0_z_0_0_yyyzzzz[k];

                g_z_0_z_0_y_yzzzzz[k] = -g_z_0_z_0_0_yzzzzz[k] * ab_y + g_z_0_z_0_0_yyzzzzz[k];

                g_z_0_z_0_y_zzzzzz[k] = -g_z_0_z_0_0_zzzzzz[k] * ab_y + g_z_0_z_0_0_yzzzzzz[k];
            }

            /// Set up 728-756 components of targeted buffer : cbuffer.data(

            auto g_z_0_z_0_z_xxxxxx = cbuffer.data(pi_geom_1010_off + 728 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxxxy = cbuffer.data(pi_geom_1010_off + 729 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxxxz = cbuffer.data(pi_geom_1010_off + 730 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxxyy = cbuffer.data(pi_geom_1010_off + 731 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxxyz = cbuffer.data(pi_geom_1010_off + 732 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxxzz = cbuffer.data(pi_geom_1010_off + 733 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxyyy = cbuffer.data(pi_geom_1010_off + 734 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxyyz = cbuffer.data(pi_geom_1010_off + 735 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxyzz = cbuffer.data(pi_geom_1010_off + 736 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxxzzz = cbuffer.data(pi_geom_1010_off + 737 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxyyyy = cbuffer.data(pi_geom_1010_off + 738 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxyyyz = cbuffer.data(pi_geom_1010_off + 739 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxyyzz = cbuffer.data(pi_geom_1010_off + 740 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxyzzz = cbuffer.data(pi_geom_1010_off + 741 * ccomps * dcomps);

            auto g_z_0_z_0_z_xxzzzz = cbuffer.data(pi_geom_1010_off + 742 * ccomps * dcomps);

            auto g_z_0_z_0_z_xyyyyy = cbuffer.data(pi_geom_1010_off + 743 * ccomps * dcomps);

            auto g_z_0_z_0_z_xyyyyz = cbuffer.data(pi_geom_1010_off + 744 * ccomps * dcomps);

            auto g_z_0_z_0_z_xyyyzz = cbuffer.data(pi_geom_1010_off + 745 * ccomps * dcomps);

            auto g_z_0_z_0_z_xyyzzz = cbuffer.data(pi_geom_1010_off + 746 * ccomps * dcomps);

            auto g_z_0_z_0_z_xyzzzz = cbuffer.data(pi_geom_1010_off + 747 * ccomps * dcomps);

            auto g_z_0_z_0_z_xzzzzz = cbuffer.data(pi_geom_1010_off + 748 * ccomps * dcomps);

            auto g_z_0_z_0_z_yyyyyy = cbuffer.data(pi_geom_1010_off + 749 * ccomps * dcomps);

            auto g_z_0_z_0_z_yyyyyz = cbuffer.data(pi_geom_1010_off + 750 * ccomps * dcomps);

            auto g_z_0_z_0_z_yyyyzz = cbuffer.data(pi_geom_1010_off + 751 * ccomps * dcomps);

            auto g_z_0_z_0_z_yyyzzz = cbuffer.data(pi_geom_1010_off + 752 * ccomps * dcomps);

            auto g_z_0_z_0_z_yyzzzz = cbuffer.data(pi_geom_1010_off + 753 * ccomps * dcomps);

            auto g_z_0_z_0_z_yzzzzz = cbuffer.data(pi_geom_1010_off + 754 * ccomps * dcomps);

            auto g_z_0_z_0_z_zzzzzz = cbuffer.data(pi_geom_1010_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_0_z_0_0_xxxxxx, g_0_0_z_0_0_xxxxxy, g_0_0_z_0_0_xxxxxz, g_0_0_z_0_0_xxxxyy, g_0_0_z_0_0_xxxxyz, g_0_0_z_0_0_xxxxzz, g_0_0_z_0_0_xxxyyy, g_0_0_z_0_0_xxxyyz, g_0_0_z_0_0_xxxyzz, g_0_0_z_0_0_xxxzzz, g_0_0_z_0_0_xxyyyy, g_0_0_z_0_0_xxyyyz, g_0_0_z_0_0_xxyyzz, g_0_0_z_0_0_xxyzzz, g_0_0_z_0_0_xxzzzz, g_0_0_z_0_0_xyyyyy, g_0_0_z_0_0_xyyyyz, g_0_0_z_0_0_xyyyzz, g_0_0_z_0_0_xyyzzz, g_0_0_z_0_0_xyzzzz, g_0_0_z_0_0_xzzzzz, g_0_0_z_0_0_yyyyyy, g_0_0_z_0_0_yyyyyz, g_0_0_z_0_0_yyyyzz, g_0_0_z_0_0_yyyzzz, g_0_0_z_0_0_yyzzzz, g_0_0_z_0_0_yzzzzz, g_0_0_z_0_0_zzzzzz, g_z_0_z_0_0_xxxxxx, g_z_0_z_0_0_xxxxxxz, g_z_0_z_0_0_xxxxxy, g_z_0_z_0_0_xxxxxyz, g_z_0_z_0_0_xxxxxz, g_z_0_z_0_0_xxxxxzz, g_z_0_z_0_0_xxxxyy, g_z_0_z_0_0_xxxxyyz, g_z_0_z_0_0_xxxxyz, g_z_0_z_0_0_xxxxyzz, g_z_0_z_0_0_xxxxzz, g_z_0_z_0_0_xxxxzzz, g_z_0_z_0_0_xxxyyy, g_z_0_z_0_0_xxxyyyz, g_z_0_z_0_0_xxxyyz, g_z_0_z_0_0_xxxyyzz, g_z_0_z_0_0_xxxyzz, g_z_0_z_0_0_xxxyzzz, g_z_0_z_0_0_xxxzzz, g_z_0_z_0_0_xxxzzzz, g_z_0_z_0_0_xxyyyy, g_z_0_z_0_0_xxyyyyz, g_z_0_z_0_0_xxyyyz, g_z_0_z_0_0_xxyyyzz, g_z_0_z_0_0_xxyyzz, g_z_0_z_0_0_xxyyzzz, g_z_0_z_0_0_xxyzzz, g_z_0_z_0_0_xxyzzzz, g_z_0_z_0_0_xxzzzz, g_z_0_z_0_0_xxzzzzz, g_z_0_z_0_0_xyyyyy, g_z_0_z_0_0_xyyyyyz, g_z_0_z_0_0_xyyyyz, g_z_0_z_0_0_xyyyyzz, g_z_0_z_0_0_xyyyzz, g_z_0_z_0_0_xyyyzzz, g_z_0_z_0_0_xyyzzz, g_z_0_z_0_0_xyyzzzz, g_z_0_z_0_0_xyzzzz, g_z_0_z_0_0_xyzzzzz, g_z_0_z_0_0_xzzzzz, g_z_0_z_0_0_xzzzzzz, g_z_0_z_0_0_yyyyyy, g_z_0_z_0_0_yyyyyyz, g_z_0_z_0_0_yyyyyz, g_z_0_z_0_0_yyyyyzz, g_z_0_z_0_0_yyyyzz, g_z_0_z_0_0_yyyyzzz, g_z_0_z_0_0_yyyzzz, g_z_0_z_0_0_yyyzzzz, g_z_0_z_0_0_yyzzzz, g_z_0_z_0_0_yyzzzzz, g_z_0_z_0_0_yzzzzz, g_z_0_z_0_0_yzzzzzz, g_z_0_z_0_0_zzzzzz, g_z_0_z_0_0_zzzzzzz, g_z_0_z_0_z_xxxxxx, g_z_0_z_0_z_xxxxxy, g_z_0_z_0_z_xxxxxz, g_z_0_z_0_z_xxxxyy, g_z_0_z_0_z_xxxxyz, g_z_0_z_0_z_xxxxzz, g_z_0_z_0_z_xxxyyy, g_z_0_z_0_z_xxxyyz, g_z_0_z_0_z_xxxyzz, g_z_0_z_0_z_xxxzzz, g_z_0_z_0_z_xxyyyy, g_z_0_z_0_z_xxyyyz, g_z_0_z_0_z_xxyyzz, g_z_0_z_0_z_xxyzzz, g_z_0_z_0_z_xxzzzz, g_z_0_z_0_z_xyyyyy, g_z_0_z_0_z_xyyyyz, g_z_0_z_0_z_xyyyzz, g_z_0_z_0_z_xyyzzz, g_z_0_z_0_z_xyzzzz, g_z_0_z_0_z_xzzzzz, g_z_0_z_0_z_yyyyyy, g_z_0_z_0_z_yyyyyz, g_z_0_z_0_z_yyyyzz, g_z_0_z_0_z_yyyzzz, g_z_0_z_0_z_yyzzzz, g_z_0_z_0_z_yzzzzz, g_z_0_z_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_z_0_z_xxxxxx[k] = -g_0_0_z_0_0_xxxxxx[k] - g_z_0_z_0_0_xxxxxx[k] * ab_z + g_z_0_z_0_0_xxxxxxz[k];

                g_z_0_z_0_z_xxxxxy[k] = -g_0_0_z_0_0_xxxxxy[k] - g_z_0_z_0_0_xxxxxy[k] * ab_z + g_z_0_z_0_0_xxxxxyz[k];

                g_z_0_z_0_z_xxxxxz[k] = -g_0_0_z_0_0_xxxxxz[k] - g_z_0_z_0_0_xxxxxz[k] * ab_z + g_z_0_z_0_0_xxxxxzz[k];

                g_z_0_z_0_z_xxxxyy[k] = -g_0_0_z_0_0_xxxxyy[k] - g_z_0_z_0_0_xxxxyy[k] * ab_z + g_z_0_z_0_0_xxxxyyz[k];

                g_z_0_z_0_z_xxxxyz[k] = -g_0_0_z_0_0_xxxxyz[k] - g_z_0_z_0_0_xxxxyz[k] * ab_z + g_z_0_z_0_0_xxxxyzz[k];

                g_z_0_z_0_z_xxxxzz[k] = -g_0_0_z_0_0_xxxxzz[k] - g_z_0_z_0_0_xxxxzz[k] * ab_z + g_z_0_z_0_0_xxxxzzz[k];

                g_z_0_z_0_z_xxxyyy[k] = -g_0_0_z_0_0_xxxyyy[k] - g_z_0_z_0_0_xxxyyy[k] * ab_z + g_z_0_z_0_0_xxxyyyz[k];

                g_z_0_z_0_z_xxxyyz[k] = -g_0_0_z_0_0_xxxyyz[k] - g_z_0_z_0_0_xxxyyz[k] * ab_z + g_z_0_z_0_0_xxxyyzz[k];

                g_z_0_z_0_z_xxxyzz[k] = -g_0_0_z_0_0_xxxyzz[k] - g_z_0_z_0_0_xxxyzz[k] * ab_z + g_z_0_z_0_0_xxxyzzz[k];

                g_z_0_z_0_z_xxxzzz[k] = -g_0_0_z_0_0_xxxzzz[k] - g_z_0_z_0_0_xxxzzz[k] * ab_z + g_z_0_z_0_0_xxxzzzz[k];

                g_z_0_z_0_z_xxyyyy[k] = -g_0_0_z_0_0_xxyyyy[k] - g_z_0_z_0_0_xxyyyy[k] * ab_z + g_z_0_z_0_0_xxyyyyz[k];

                g_z_0_z_0_z_xxyyyz[k] = -g_0_0_z_0_0_xxyyyz[k] - g_z_0_z_0_0_xxyyyz[k] * ab_z + g_z_0_z_0_0_xxyyyzz[k];

                g_z_0_z_0_z_xxyyzz[k] = -g_0_0_z_0_0_xxyyzz[k] - g_z_0_z_0_0_xxyyzz[k] * ab_z + g_z_0_z_0_0_xxyyzzz[k];

                g_z_0_z_0_z_xxyzzz[k] = -g_0_0_z_0_0_xxyzzz[k] - g_z_0_z_0_0_xxyzzz[k] * ab_z + g_z_0_z_0_0_xxyzzzz[k];

                g_z_0_z_0_z_xxzzzz[k] = -g_0_0_z_0_0_xxzzzz[k] - g_z_0_z_0_0_xxzzzz[k] * ab_z + g_z_0_z_0_0_xxzzzzz[k];

                g_z_0_z_0_z_xyyyyy[k] = -g_0_0_z_0_0_xyyyyy[k] - g_z_0_z_0_0_xyyyyy[k] * ab_z + g_z_0_z_0_0_xyyyyyz[k];

                g_z_0_z_0_z_xyyyyz[k] = -g_0_0_z_0_0_xyyyyz[k] - g_z_0_z_0_0_xyyyyz[k] * ab_z + g_z_0_z_0_0_xyyyyzz[k];

                g_z_0_z_0_z_xyyyzz[k] = -g_0_0_z_0_0_xyyyzz[k] - g_z_0_z_0_0_xyyyzz[k] * ab_z + g_z_0_z_0_0_xyyyzzz[k];

                g_z_0_z_0_z_xyyzzz[k] = -g_0_0_z_0_0_xyyzzz[k] - g_z_0_z_0_0_xyyzzz[k] * ab_z + g_z_0_z_0_0_xyyzzzz[k];

                g_z_0_z_0_z_xyzzzz[k] = -g_0_0_z_0_0_xyzzzz[k] - g_z_0_z_0_0_xyzzzz[k] * ab_z + g_z_0_z_0_0_xyzzzzz[k];

                g_z_0_z_0_z_xzzzzz[k] = -g_0_0_z_0_0_xzzzzz[k] - g_z_0_z_0_0_xzzzzz[k] * ab_z + g_z_0_z_0_0_xzzzzzz[k];

                g_z_0_z_0_z_yyyyyy[k] = -g_0_0_z_0_0_yyyyyy[k] - g_z_0_z_0_0_yyyyyy[k] * ab_z + g_z_0_z_0_0_yyyyyyz[k];

                g_z_0_z_0_z_yyyyyz[k] = -g_0_0_z_0_0_yyyyyz[k] - g_z_0_z_0_0_yyyyyz[k] * ab_z + g_z_0_z_0_0_yyyyyzz[k];

                g_z_0_z_0_z_yyyyzz[k] = -g_0_0_z_0_0_yyyyzz[k] - g_z_0_z_0_0_yyyyzz[k] * ab_z + g_z_0_z_0_0_yyyyzzz[k];

                g_z_0_z_0_z_yyyzzz[k] = -g_0_0_z_0_0_yyyzzz[k] - g_z_0_z_0_0_yyyzzz[k] * ab_z + g_z_0_z_0_0_yyyzzzz[k];

                g_z_0_z_0_z_yyzzzz[k] = -g_0_0_z_0_0_yyzzzz[k] - g_z_0_z_0_0_yyzzzz[k] * ab_z + g_z_0_z_0_0_yyzzzzz[k];

                g_z_0_z_0_z_yzzzzz[k] = -g_0_0_z_0_0_yzzzzz[k] - g_z_0_z_0_0_yzzzzz[k] * ab_z + g_z_0_z_0_0_yzzzzzz[k];

                g_z_0_z_0_z_zzzzzz[k] = -g_0_0_z_0_0_zzzzzz[k] - g_z_0_z_0_0_zzzzzz[k] * ab_z + g_z_0_z_0_0_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

