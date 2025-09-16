#include "ElectronRepulsionGeom0100ContrRecKPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_kpxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_kpxx,
                                            const size_t idx_ipxx,
                                            const size_t idx_geom_01_ipxx,
                                            const size_t idx_geom_01_idxx,
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
            /// Set up components of auxilary buffer : IPSS

            const auto ip_off = idx_ipxx + i * dcomps + j;

            auto g_xxxxxx_x = cbuffer.data(ip_off + 0 * ccomps * dcomps);

            auto g_xxxxxx_y = cbuffer.data(ip_off + 1 * ccomps * dcomps);

            auto g_xxxxxx_z = cbuffer.data(ip_off + 2 * ccomps * dcomps);

            auto g_xxxxxy_x = cbuffer.data(ip_off + 3 * ccomps * dcomps);

            auto g_xxxxxy_y = cbuffer.data(ip_off + 4 * ccomps * dcomps);

            auto g_xxxxxy_z = cbuffer.data(ip_off + 5 * ccomps * dcomps);

            auto g_xxxxxz_x = cbuffer.data(ip_off + 6 * ccomps * dcomps);

            auto g_xxxxxz_y = cbuffer.data(ip_off + 7 * ccomps * dcomps);

            auto g_xxxxxz_z = cbuffer.data(ip_off + 8 * ccomps * dcomps);

            auto g_xxxxyy_x = cbuffer.data(ip_off + 9 * ccomps * dcomps);

            auto g_xxxxyy_y = cbuffer.data(ip_off + 10 * ccomps * dcomps);

            auto g_xxxxyy_z = cbuffer.data(ip_off + 11 * ccomps * dcomps);

            auto g_xxxxyz_x = cbuffer.data(ip_off + 12 * ccomps * dcomps);

            auto g_xxxxyz_y = cbuffer.data(ip_off + 13 * ccomps * dcomps);

            auto g_xxxxyz_z = cbuffer.data(ip_off + 14 * ccomps * dcomps);

            auto g_xxxxzz_x = cbuffer.data(ip_off + 15 * ccomps * dcomps);

            auto g_xxxxzz_y = cbuffer.data(ip_off + 16 * ccomps * dcomps);

            auto g_xxxxzz_z = cbuffer.data(ip_off + 17 * ccomps * dcomps);

            auto g_xxxyyy_x = cbuffer.data(ip_off + 18 * ccomps * dcomps);

            auto g_xxxyyy_y = cbuffer.data(ip_off + 19 * ccomps * dcomps);

            auto g_xxxyyy_z = cbuffer.data(ip_off + 20 * ccomps * dcomps);

            auto g_xxxyyz_x = cbuffer.data(ip_off + 21 * ccomps * dcomps);

            auto g_xxxyyz_y = cbuffer.data(ip_off + 22 * ccomps * dcomps);

            auto g_xxxyyz_z = cbuffer.data(ip_off + 23 * ccomps * dcomps);

            auto g_xxxyzz_x = cbuffer.data(ip_off + 24 * ccomps * dcomps);

            auto g_xxxyzz_y = cbuffer.data(ip_off + 25 * ccomps * dcomps);

            auto g_xxxyzz_z = cbuffer.data(ip_off + 26 * ccomps * dcomps);

            auto g_xxxzzz_x = cbuffer.data(ip_off + 27 * ccomps * dcomps);

            auto g_xxxzzz_y = cbuffer.data(ip_off + 28 * ccomps * dcomps);

            auto g_xxxzzz_z = cbuffer.data(ip_off + 29 * ccomps * dcomps);

            auto g_xxyyyy_x = cbuffer.data(ip_off + 30 * ccomps * dcomps);

            auto g_xxyyyy_y = cbuffer.data(ip_off + 31 * ccomps * dcomps);

            auto g_xxyyyy_z = cbuffer.data(ip_off + 32 * ccomps * dcomps);

            auto g_xxyyyz_x = cbuffer.data(ip_off + 33 * ccomps * dcomps);

            auto g_xxyyyz_y = cbuffer.data(ip_off + 34 * ccomps * dcomps);

            auto g_xxyyyz_z = cbuffer.data(ip_off + 35 * ccomps * dcomps);

            auto g_xxyyzz_x = cbuffer.data(ip_off + 36 * ccomps * dcomps);

            auto g_xxyyzz_y = cbuffer.data(ip_off + 37 * ccomps * dcomps);

            auto g_xxyyzz_z = cbuffer.data(ip_off + 38 * ccomps * dcomps);

            auto g_xxyzzz_x = cbuffer.data(ip_off + 39 * ccomps * dcomps);

            auto g_xxyzzz_y = cbuffer.data(ip_off + 40 * ccomps * dcomps);

            auto g_xxyzzz_z = cbuffer.data(ip_off + 41 * ccomps * dcomps);

            auto g_xxzzzz_x = cbuffer.data(ip_off + 42 * ccomps * dcomps);

            auto g_xxzzzz_y = cbuffer.data(ip_off + 43 * ccomps * dcomps);

            auto g_xxzzzz_z = cbuffer.data(ip_off + 44 * ccomps * dcomps);

            auto g_xyyyyy_x = cbuffer.data(ip_off + 45 * ccomps * dcomps);

            auto g_xyyyyy_y = cbuffer.data(ip_off + 46 * ccomps * dcomps);

            auto g_xyyyyy_z = cbuffer.data(ip_off + 47 * ccomps * dcomps);

            auto g_xyyyyz_x = cbuffer.data(ip_off + 48 * ccomps * dcomps);

            auto g_xyyyyz_y = cbuffer.data(ip_off + 49 * ccomps * dcomps);

            auto g_xyyyyz_z = cbuffer.data(ip_off + 50 * ccomps * dcomps);

            auto g_xyyyzz_x = cbuffer.data(ip_off + 51 * ccomps * dcomps);

            auto g_xyyyzz_y = cbuffer.data(ip_off + 52 * ccomps * dcomps);

            auto g_xyyyzz_z = cbuffer.data(ip_off + 53 * ccomps * dcomps);

            auto g_xyyzzz_x = cbuffer.data(ip_off + 54 * ccomps * dcomps);

            auto g_xyyzzz_y = cbuffer.data(ip_off + 55 * ccomps * dcomps);

            auto g_xyyzzz_z = cbuffer.data(ip_off + 56 * ccomps * dcomps);

            auto g_xyzzzz_x = cbuffer.data(ip_off + 57 * ccomps * dcomps);

            auto g_xyzzzz_y = cbuffer.data(ip_off + 58 * ccomps * dcomps);

            auto g_xyzzzz_z = cbuffer.data(ip_off + 59 * ccomps * dcomps);

            auto g_xzzzzz_x = cbuffer.data(ip_off + 60 * ccomps * dcomps);

            auto g_xzzzzz_y = cbuffer.data(ip_off + 61 * ccomps * dcomps);

            auto g_xzzzzz_z = cbuffer.data(ip_off + 62 * ccomps * dcomps);

            auto g_yyyyyy_x = cbuffer.data(ip_off + 63 * ccomps * dcomps);

            auto g_yyyyyy_y = cbuffer.data(ip_off + 64 * ccomps * dcomps);

            auto g_yyyyyy_z = cbuffer.data(ip_off + 65 * ccomps * dcomps);

            auto g_yyyyyz_x = cbuffer.data(ip_off + 66 * ccomps * dcomps);

            auto g_yyyyyz_y = cbuffer.data(ip_off + 67 * ccomps * dcomps);

            auto g_yyyyyz_z = cbuffer.data(ip_off + 68 * ccomps * dcomps);

            auto g_yyyyzz_x = cbuffer.data(ip_off + 69 * ccomps * dcomps);

            auto g_yyyyzz_y = cbuffer.data(ip_off + 70 * ccomps * dcomps);

            auto g_yyyyzz_z = cbuffer.data(ip_off + 71 * ccomps * dcomps);

            auto g_yyyzzz_x = cbuffer.data(ip_off + 72 * ccomps * dcomps);

            auto g_yyyzzz_y = cbuffer.data(ip_off + 73 * ccomps * dcomps);

            auto g_yyyzzz_z = cbuffer.data(ip_off + 74 * ccomps * dcomps);

            auto g_yyzzzz_x = cbuffer.data(ip_off + 75 * ccomps * dcomps);

            auto g_yyzzzz_y = cbuffer.data(ip_off + 76 * ccomps * dcomps);

            auto g_yyzzzz_z = cbuffer.data(ip_off + 77 * ccomps * dcomps);

            auto g_yzzzzz_x = cbuffer.data(ip_off + 78 * ccomps * dcomps);

            auto g_yzzzzz_y = cbuffer.data(ip_off + 79 * ccomps * dcomps);

            auto g_yzzzzz_z = cbuffer.data(ip_off + 80 * ccomps * dcomps);

            auto g_zzzzzz_x = cbuffer.data(ip_off + 81 * ccomps * dcomps);

            auto g_zzzzzz_y = cbuffer.data(ip_off + 82 * ccomps * dcomps);

            auto g_zzzzzz_z = cbuffer.data(ip_off + 83 * ccomps * dcomps);

            /// Set up components of auxilary buffer : IPSS

            const auto ip_geom_01_off = idx_geom_01_ipxx + i * dcomps + j;

            auto g_0_x_xxxxxx_x = cbuffer.data(ip_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxx_y = cbuffer.data(ip_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxx_z = cbuffer.data(ip_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxy_x = cbuffer.data(ip_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxy_y = cbuffer.data(ip_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxy_z = cbuffer.data(ip_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxxz_x = cbuffer.data(ip_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxz_y = cbuffer.data(ip_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxz_z = cbuffer.data(ip_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxyy_x = cbuffer.data(ip_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxyy_y = cbuffer.data(ip_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxyy_z = cbuffer.data(ip_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxyz_x = cbuffer.data(ip_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxyz_y = cbuffer.data(ip_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxyz_z = cbuffer.data(ip_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxzz_x = cbuffer.data(ip_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxzz_y = cbuffer.data(ip_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxzz_z = cbuffer.data(ip_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxyyy_x = cbuffer.data(ip_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxyyy_y = cbuffer.data(ip_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxyyy_z = cbuffer.data(ip_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxyyz_x = cbuffer.data(ip_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxyyz_y = cbuffer.data(ip_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxyyz_z = cbuffer.data(ip_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxyzz_x = cbuffer.data(ip_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxyzz_y = cbuffer.data(ip_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxyzz_z = cbuffer.data(ip_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxzzz_x = cbuffer.data(ip_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxzzz_y = cbuffer.data(ip_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxzzz_z = cbuffer.data(ip_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxyyyy_x = cbuffer.data(ip_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxyyyy_y = cbuffer.data(ip_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxyyyy_z = cbuffer.data(ip_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxyyyz_x = cbuffer.data(ip_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxyyyz_y = cbuffer.data(ip_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxyyyz_z = cbuffer.data(ip_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxyyzz_x = cbuffer.data(ip_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxyyzz_y = cbuffer.data(ip_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxyyzz_z = cbuffer.data(ip_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxyzzz_x = cbuffer.data(ip_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxyzzz_y = cbuffer.data(ip_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxyzzz_z = cbuffer.data(ip_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxzzzz_x = cbuffer.data(ip_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxzzzz_y = cbuffer.data(ip_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxzzzz_z = cbuffer.data(ip_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xyyyyy_x = cbuffer.data(ip_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xyyyyy_y = cbuffer.data(ip_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xyyyyy_z = cbuffer.data(ip_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xyyyyz_x = cbuffer.data(ip_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xyyyyz_y = cbuffer.data(ip_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xyyyyz_z = cbuffer.data(ip_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xyyyzz_x = cbuffer.data(ip_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xyyyzz_y = cbuffer.data(ip_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xyyyzz_z = cbuffer.data(ip_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xyyzzz_x = cbuffer.data(ip_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xyyzzz_y = cbuffer.data(ip_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xyyzzz_z = cbuffer.data(ip_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xyzzzz_x = cbuffer.data(ip_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xyzzzz_y = cbuffer.data(ip_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xyzzzz_z = cbuffer.data(ip_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xzzzzz_x = cbuffer.data(ip_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xzzzzz_y = cbuffer.data(ip_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xzzzzz_z = cbuffer.data(ip_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_yyyyyy_x = cbuffer.data(ip_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_yyyyyy_y = cbuffer.data(ip_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_yyyyyy_z = cbuffer.data(ip_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_yyyyyz_x = cbuffer.data(ip_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_yyyyyz_y = cbuffer.data(ip_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_yyyyyz_z = cbuffer.data(ip_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_yyyyzz_x = cbuffer.data(ip_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_yyyyzz_y = cbuffer.data(ip_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_yyyyzz_z = cbuffer.data(ip_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_yyyzzz_x = cbuffer.data(ip_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_yyyzzz_y = cbuffer.data(ip_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_yyyzzz_z = cbuffer.data(ip_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_yyzzzz_x = cbuffer.data(ip_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_yyzzzz_y = cbuffer.data(ip_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_yyzzzz_z = cbuffer.data(ip_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_yzzzzz_x = cbuffer.data(ip_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_yzzzzz_y = cbuffer.data(ip_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_yzzzzz_z = cbuffer.data(ip_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_zzzzzz_x = cbuffer.data(ip_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_zzzzzz_y = cbuffer.data(ip_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_zzzzzz_z = cbuffer.data(ip_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_y_xxxxxx_x = cbuffer.data(ip_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_xxxxxx_y = cbuffer.data(ip_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_xxxxxx_z = cbuffer.data(ip_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_xxxxxy_x = cbuffer.data(ip_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_xxxxxy_y = cbuffer.data(ip_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_xxxxxy_z = cbuffer.data(ip_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_xxxxxz_x = cbuffer.data(ip_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_xxxxxz_y = cbuffer.data(ip_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_xxxxxz_z = cbuffer.data(ip_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_xxxxyy_x = cbuffer.data(ip_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_xxxxyy_y = cbuffer.data(ip_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_xxxxyy_z = cbuffer.data(ip_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_xxxxyz_x = cbuffer.data(ip_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_xxxxyz_y = cbuffer.data(ip_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_xxxxyz_z = cbuffer.data(ip_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_xxxxzz_x = cbuffer.data(ip_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_xxxxzz_y = cbuffer.data(ip_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_xxxxzz_z = cbuffer.data(ip_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_xxxyyy_x = cbuffer.data(ip_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_xxxyyy_y = cbuffer.data(ip_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_xxxyyy_z = cbuffer.data(ip_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_xxxyyz_x = cbuffer.data(ip_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_xxxyyz_y = cbuffer.data(ip_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_xxxyyz_z = cbuffer.data(ip_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_xxxyzz_x = cbuffer.data(ip_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xxxyzz_y = cbuffer.data(ip_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xxxyzz_z = cbuffer.data(ip_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_xxxzzz_x = cbuffer.data(ip_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xxxzzz_y = cbuffer.data(ip_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xxxzzz_z = cbuffer.data(ip_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_xxyyyy_x = cbuffer.data(ip_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xxyyyy_y = cbuffer.data(ip_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xxyyyy_z = cbuffer.data(ip_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_xxyyyz_x = cbuffer.data(ip_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xxyyyz_y = cbuffer.data(ip_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xxyyyz_z = cbuffer.data(ip_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_xxyyzz_x = cbuffer.data(ip_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xxyyzz_y = cbuffer.data(ip_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xxyyzz_z = cbuffer.data(ip_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_xxyzzz_x = cbuffer.data(ip_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xxyzzz_y = cbuffer.data(ip_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xxyzzz_z = cbuffer.data(ip_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_y_xxzzzz_x = cbuffer.data(ip_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xxzzzz_y = cbuffer.data(ip_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xxzzzz_z = cbuffer.data(ip_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_y_xyyyyy_x = cbuffer.data(ip_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xyyyyy_y = cbuffer.data(ip_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xyyyyy_z = cbuffer.data(ip_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_y_xyyyyz_x = cbuffer.data(ip_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xyyyyz_y = cbuffer.data(ip_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xyyyyz_z = cbuffer.data(ip_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_y_xyyyzz_x = cbuffer.data(ip_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xyyyzz_y = cbuffer.data(ip_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xyyyzz_z = cbuffer.data(ip_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_y_xyyzzz_x = cbuffer.data(ip_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xyyzzz_y = cbuffer.data(ip_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xyyzzz_z = cbuffer.data(ip_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_y_xyzzzz_x = cbuffer.data(ip_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xyzzzz_y = cbuffer.data(ip_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xyzzzz_z = cbuffer.data(ip_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_y_xzzzzz_x = cbuffer.data(ip_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xzzzzz_y = cbuffer.data(ip_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xzzzzz_z = cbuffer.data(ip_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_y_yyyyyy_x = cbuffer.data(ip_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_yyyyyy_y = cbuffer.data(ip_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_yyyyyy_z = cbuffer.data(ip_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_y_yyyyyz_x = cbuffer.data(ip_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_yyyyyz_y = cbuffer.data(ip_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_yyyyyz_z = cbuffer.data(ip_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_y_yyyyzz_x = cbuffer.data(ip_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_yyyyzz_y = cbuffer.data(ip_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_yyyyzz_z = cbuffer.data(ip_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_y_yyyzzz_x = cbuffer.data(ip_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_yyyzzz_y = cbuffer.data(ip_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_yyyzzz_z = cbuffer.data(ip_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_y_yyzzzz_x = cbuffer.data(ip_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_yyzzzz_y = cbuffer.data(ip_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_yyzzzz_z = cbuffer.data(ip_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_y_yzzzzz_x = cbuffer.data(ip_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_yzzzzz_y = cbuffer.data(ip_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_yzzzzz_z = cbuffer.data(ip_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_y_zzzzzz_x = cbuffer.data(ip_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_zzzzzz_y = cbuffer.data(ip_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_zzzzzz_z = cbuffer.data(ip_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_z_xxxxxx_x = cbuffer.data(ip_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_z_xxxxxx_y = cbuffer.data(ip_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_z_xxxxxx_z = cbuffer.data(ip_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_z_xxxxxy_x = cbuffer.data(ip_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_z_xxxxxy_y = cbuffer.data(ip_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_z_xxxxxy_z = cbuffer.data(ip_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_z_xxxxxz_x = cbuffer.data(ip_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_z_xxxxxz_y = cbuffer.data(ip_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_z_xxxxxz_z = cbuffer.data(ip_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_z_xxxxyy_x = cbuffer.data(ip_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_z_xxxxyy_y = cbuffer.data(ip_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_z_xxxxyy_z = cbuffer.data(ip_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_xxxxyz_x = cbuffer.data(ip_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_xxxxyz_y = cbuffer.data(ip_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_xxxxyz_z = cbuffer.data(ip_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_xxxxzz_x = cbuffer.data(ip_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_xxxxzz_y = cbuffer.data(ip_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_xxxxzz_z = cbuffer.data(ip_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_xxxyyy_x = cbuffer.data(ip_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_xxxyyy_y = cbuffer.data(ip_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_xxxyyy_z = cbuffer.data(ip_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_z_xxxyyz_x = cbuffer.data(ip_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_z_xxxyyz_y = cbuffer.data(ip_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_z_xxxyyz_z = cbuffer.data(ip_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_z_xxxyzz_x = cbuffer.data(ip_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_z_xxxyzz_y = cbuffer.data(ip_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_z_xxxyzz_z = cbuffer.data(ip_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_z_xxxzzz_x = cbuffer.data(ip_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_z_xxxzzz_y = cbuffer.data(ip_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_z_xxxzzz_z = cbuffer.data(ip_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_z_xxyyyy_x = cbuffer.data(ip_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_z_xxyyyy_y = cbuffer.data(ip_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_z_xxyyyy_z = cbuffer.data(ip_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_z_xxyyyz_x = cbuffer.data(ip_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_z_xxyyyz_y = cbuffer.data(ip_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_z_xxyyyz_z = cbuffer.data(ip_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_z_xxyyzz_x = cbuffer.data(ip_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_z_xxyyzz_y = cbuffer.data(ip_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_z_xxyyzz_z = cbuffer.data(ip_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_z_xxyzzz_x = cbuffer.data(ip_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_z_xxyzzz_y = cbuffer.data(ip_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_z_xxyzzz_z = cbuffer.data(ip_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_z_xxzzzz_x = cbuffer.data(ip_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_z_xxzzzz_y = cbuffer.data(ip_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_z_xxzzzz_z = cbuffer.data(ip_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_z_xyyyyy_x = cbuffer.data(ip_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_z_xyyyyy_y = cbuffer.data(ip_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_z_xyyyyy_z = cbuffer.data(ip_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_z_xyyyyz_x = cbuffer.data(ip_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xyyyyz_y = cbuffer.data(ip_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xyyyyz_z = cbuffer.data(ip_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_z_xyyyzz_x = cbuffer.data(ip_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xyyyzz_y = cbuffer.data(ip_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xyyyzz_z = cbuffer.data(ip_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_z_xyyzzz_x = cbuffer.data(ip_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xyyzzz_y = cbuffer.data(ip_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xyyzzz_z = cbuffer.data(ip_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_z_xyzzzz_x = cbuffer.data(ip_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_xyzzzz_y = cbuffer.data(ip_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_xyzzzz_z = cbuffer.data(ip_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_z_xzzzzz_x = cbuffer.data(ip_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_xzzzzz_y = cbuffer.data(ip_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_xzzzzz_z = cbuffer.data(ip_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_z_yyyyyy_x = cbuffer.data(ip_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_yyyyyy_y = cbuffer.data(ip_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_yyyyyy_z = cbuffer.data(ip_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_z_yyyyyz_x = cbuffer.data(ip_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_yyyyyz_y = cbuffer.data(ip_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_yyyyyz_z = cbuffer.data(ip_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_z_yyyyzz_x = cbuffer.data(ip_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_yyyyzz_y = cbuffer.data(ip_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_yyyyzz_z = cbuffer.data(ip_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_z_yyyzzz_x = cbuffer.data(ip_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_yyyzzz_y = cbuffer.data(ip_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_yyyzzz_z = cbuffer.data(ip_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_z_yyzzzz_x = cbuffer.data(ip_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_yyzzzz_y = cbuffer.data(ip_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_yyzzzz_z = cbuffer.data(ip_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_z_yzzzzz_x = cbuffer.data(ip_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_yzzzzz_y = cbuffer.data(ip_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_yzzzzz_z = cbuffer.data(ip_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_z_zzzzzz_x = cbuffer.data(ip_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_zzzzzz_y = cbuffer.data(ip_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_zzzzzz_z = cbuffer.data(ip_geom_01_off + 251 * ccomps * dcomps);

            /// Set up components of auxilary buffer : IDSS

            const auto id_geom_01_off = idx_geom_01_idxx + i * dcomps + j;

            auto g_0_x_xxxxxx_xx = cbuffer.data(id_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xy = cbuffer.data(id_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xz = cbuffer.data(id_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yy = cbuffer.data(id_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yz = cbuffer.data(id_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxx_zz = cbuffer.data(id_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xx = cbuffer.data(id_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xy = cbuffer.data(id_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xz = cbuffer.data(id_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yy = cbuffer.data(id_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yz = cbuffer.data(id_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxy_zz = cbuffer.data(id_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xx = cbuffer.data(id_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xy = cbuffer.data(id_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xz = cbuffer.data(id_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yy = cbuffer.data(id_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yz = cbuffer.data(id_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxz_zz = cbuffer.data(id_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xx = cbuffer.data(id_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xy = cbuffer.data(id_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xz = cbuffer.data(id_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yy = cbuffer.data(id_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yz = cbuffer.data(id_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxyy_zz = cbuffer.data(id_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xx = cbuffer.data(id_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xy = cbuffer.data(id_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xz = cbuffer.data(id_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yy = cbuffer.data(id_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yz = cbuffer.data(id_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxyz_zz = cbuffer.data(id_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xx = cbuffer.data(id_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xy = cbuffer.data(id_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xz = cbuffer.data(id_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yy = cbuffer.data(id_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yz = cbuffer.data(id_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxzz_zz = cbuffer.data(id_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xx = cbuffer.data(id_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xy = cbuffer.data(id_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xz = cbuffer.data(id_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yy = cbuffer.data(id_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yz = cbuffer.data(id_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxyyy_zz = cbuffer.data(id_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xx = cbuffer.data(id_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xy = cbuffer.data(id_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xz = cbuffer.data(id_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yy = cbuffer.data(id_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yz = cbuffer.data(id_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxyyz_zz = cbuffer.data(id_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xx = cbuffer.data(id_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xy = cbuffer.data(id_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xz = cbuffer.data(id_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yy = cbuffer.data(id_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yz = cbuffer.data(id_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxyzz_zz = cbuffer.data(id_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xx = cbuffer.data(id_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xy = cbuffer.data(id_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xz = cbuffer.data(id_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yy = cbuffer.data(id_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yz = cbuffer.data(id_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxzzz_zz = cbuffer.data(id_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xx = cbuffer.data(id_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xy = cbuffer.data(id_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xz = cbuffer.data(id_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yy = cbuffer.data(id_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yz = cbuffer.data(id_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxyyyy_zz = cbuffer.data(id_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xx = cbuffer.data(id_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xy = cbuffer.data(id_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xz = cbuffer.data(id_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yy = cbuffer.data(id_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yz = cbuffer.data(id_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxyyyz_zz = cbuffer.data(id_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xx = cbuffer.data(id_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xy = cbuffer.data(id_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xz = cbuffer.data(id_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yy = cbuffer.data(id_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yz = cbuffer.data(id_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxyyzz_zz = cbuffer.data(id_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xx = cbuffer.data(id_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xy = cbuffer.data(id_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xz = cbuffer.data(id_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yy = cbuffer.data(id_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yz = cbuffer.data(id_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxyzzz_zz = cbuffer.data(id_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xx = cbuffer.data(id_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xy = cbuffer.data(id_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xz = cbuffer.data(id_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yy = cbuffer.data(id_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yz = cbuffer.data(id_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxzzzz_zz = cbuffer.data(id_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xx = cbuffer.data(id_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xy = cbuffer.data(id_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xz = cbuffer.data(id_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yy = cbuffer.data(id_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yz = cbuffer.data(id_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xyyyyy_zz = cbuffer.data(id_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xx = cbuffer.data(id_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xy = cbuffer.data(id_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xz = cbuffer.data(id_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yy = cbuffer.data(id_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yz = cbuffer.data(id_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xyyyyz_zz = cbuffer.data(id_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xx = cbuffer.data(id_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xy = cbuffer.data(id_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xz = cbuffer.data(id_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yy = cbuffer.data(id_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yz = cbuffer.data(id_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xyyyzz_zz = cbuffer.data(id_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xx = cbuffer.data(id_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xy = cbuffer.data(id_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xz = cbuffer.data(id_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yy = cbuffer.data(id_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yz = cbuffer.data(id_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xyyzzz_zz = cbuffer.data(id_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xx = cbuffer.data(id_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xy = cbuffer.data(id_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xz = cbuffer.data(id_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yy = cbuffer.data(id_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yz = cbuffer.data(id_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xyzzzz_zz = cbuffer.data(id_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xx = cbuffer.data(id_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xy = cbuffer.data(id_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xz = cbuffer.data(id_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yy = cbuffer.data(id_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yz = cbuffer.data(id_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xzzzzz_zz = cbuffer.data(id_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xx = cbuffer.data(id_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xy = cbuffer.data(id_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xz = cbuffer.data(id_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yy = cbuffer.data(id_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yz = cbuffer.data(id_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_yyyyyy_zz = cbuffer.data(id_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xx = cbuffer.data(id_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xy = cbuffer.data(id_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xz = cbuffer.data(id_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yy = cbuffer.data(id_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yz = cbuffer.data(id_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_yyyyyz_zz = cbuffer.data(id_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xx = cbuffer.data(id_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xy = cbuffer.data(id_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xz = cbuffer.data(id_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yy = cbuffer.data(id_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yz = cbuffer.data(id_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_yyyyzz_zz = cbuffer.data(id_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xx = cbuffer.data(id_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xy = cbuffer.data(id_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xz = cbuffer.data(id_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yy = cbuffer.data(id_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yz = cbuffer.data(id_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_yyyzzz_zz = cbuffer.data(id_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xx = cbuffer.data(id_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xy = cbuffer.data(id_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xz = cbuffer.data(id_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yy = cbuffer.data(id_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yz = cbuffer.data(id_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_yyzzzz_zz = cbuffer.data(id_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xx = cbuffer.data(id_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xy = cbuffer.data(id_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xz = cbuffer.data(id_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yy = cbuffer.data(id_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yz = cbuffer.data(id_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_yzzzzz_zz = cbuffer.data(id_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xx = cbuffer.data(id_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xy = cbuffer.data(id_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xz = cbuffer.data(id_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yy = cbuffer.data(id_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yz = cbuffer.data(id_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_zzzzzz_zz = cbuffer.data(id_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xx = cbuffer.data(id_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xy = cbuffer.data(id_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xz = cbuffer.data(id_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yy = cbuffer.data(id_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yz = cbuffer.data(id_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xxxxxx_zz = cbuffer.data(id_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xx = cbuffer.data(id_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xy = cbuffer.data(id_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xz = cbuffer.data(id_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yy = cbuffer.data(id_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yz = cbuffer.data(id_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xxxxxy_zz = cbuffer.data(id_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xx = cbuffer.data(id_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xy = cbuffer.data(id_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xz = cbuffer.data(id_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yy = cbuffer.data(id_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yz = cbuffer.data(id_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xxxxxz_zz = cbuffer.data(id_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xx = cbuffer.data(id_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xy = cbuffer.data(id_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xz = cbuffer.data(id_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yy = cbuffer.data(id_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yz = cbuffer.data(id_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xxxxyy_zz = cbuffer.data(id_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xx = cbuffer.data(id_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xy = cbuffer.data(id_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xz = cbuffer.data(id_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yy = cbuffer.data(id_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yz = cbuffer.data(id_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_xxxxyz_zz = cbuffer.data(id_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xx = cbuffer.data(id_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xy = cbuffer.data(id_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xz = cbuffer.data(id_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yy = cbuffer.data(id_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yz = cbuffer.data(id_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_xxxxzz_zz = cbuffer.data(id_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xx = cbuffer.data(id_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xy = cbuffer.data(id_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xz = cbuffer.data(id_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yy = cbuffer.data(id_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yz = cbuffer.data(id_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_xxxyyy_zz = cbuffer.data(id_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xx = cbuffer.data(id_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xy = cbuffer.data(id_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xz = cbuffer.data(id_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yy = cbuffer.data(id_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yz = cbuffer.data(id_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_xxxyyz_zz = cbuffer.data(id_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xx = cbuffer.data(id_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xy = cbuffer.data(id_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xz = cbuffer.data(id_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yy = cbuffer.data(id_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yz = cbuffer.data(id_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_y_xxxyzz_zz = cbuffer.data(id_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xx = cbuffer.data(id_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xy = cbuffer.data(id_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xz = cbuffer.data(id_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yy = cbuffer.data(id_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yz = cbuffer.data(id_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_y_xxxzzz_zz = cbuffer.data(id_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xx = cbuffer.data(id_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xy = cbuffer.data(id_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xz = cbuffer.data(id_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yy = cbuffer.data(id_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yz = cbuffer.data(id_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_y_xxyyyy_zz = cbuffer.data(id_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xx = cbuffer.data(id_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xy = cbuffer.data(id_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xz = cbuffer.data(id_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yy = cbuffer.data(id_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yz = cbuffer.data(id_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_y_xxyyyz_zz = cbuffer.data(id_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xx = cbuffer.data(id_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xy = cbuffer.data(id_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xz = cbuffer.data(id_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yy = cbuffer.data(id_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yz = cbuffer.data(id_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_y_xxyyzz_zz = cbuffer.data(id_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xx = cbuffer.data(id_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xy = cbuffer.data(id_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xz = cbuffer.data(id_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yy = cbuffer.data(id_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yz = cbuffer.data(id_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_y_xxyzzz_zz = cbuffer.data(id_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xx = cbuffer.data(id_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xy = cbuffer.data(id_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xz = cbuffer.data(id_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yy = cbuffer.data(id_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yz = cbuffer.data(id_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_y_xxzzzz_zz = cbuffer.data(id_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xx = cbuffer.data(id_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xy = cbuffer.data(id_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xz = cbuffer.data(id_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yy = cbuffer.data(id_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yz = cbuffer.data(id_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_y_xyyyyy_zz = cbuffer.data(id_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xx = cbuffer.data(id_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xy = cbuffer.data(id_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xz = cbuffer.data(id_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yy = cbuffer.data(id_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yz = cbuffer.data(id_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_y_xyyyyz_zz = cbuffer.data(id_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xx = cbuffer.data(id_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xy = cbuffer.data(id_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xz = cbuffer.data(id_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yy = cbuffer.data(id_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yz = cbuffer.data(id_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_y_xyyyzz_zz = cbuffer.data(id_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xx = cbuffer.data(id_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xy = cbuffer.data(id_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xz = cbuffer.data(id_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yy = cbuffer.data(id_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yz = cbuffer.data(id_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_y_xyyzzz_zz = cbuffer.data(id_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xx = cbuffer.data(id_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xy = cbuffer.data(id_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xz = cbuffer.data(id_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yy = cbuffer.data(id_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yz = cbuffer.data(id_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_y_xyzzzz_zz = cbuffer.data(id_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xx = cbuffer.data(id_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xy = cbuffer.data(id_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xz = cbuffer.data(id_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yy = cbuffer.data(id_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yz = cbuffer.data(id_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_y_xzzzzz_zz = cbuffer.data(id_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xx = cbuffer.data(id_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xy = cbuffer.data(id_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xz = cbuffer.data(id_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yy = cbuffer.data(id_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yz = cbuffer.data(id_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_y_yyyyyy_zz = cbuffer.data(id_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xx = cbuffer.data(id_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xy = cbuffer.data(id_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xz = cbuffer.data(id_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yy = cbuffer.data(id_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yz = cbuffer.data(id_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_y_yyyyyz_zz = cbuffer.data(id_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xx = cbuffer.data(id_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xy = cbuffer.data(id_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xz = cbuffer.data(id_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yy = cbuffer.data(id_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yz = cbuffer.data(id_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_y_yyyyzz_zz = cbuffer.data(id_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xx = cbuffer.data(id_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xy = cbuffer.data(id_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xz = cbuffer.data(id_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yy = cbuffer.data(id_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yz = cbuffer.data(id_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_yyyzzz_zz = cbuffer.data(id_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xx = cbuffer.data(id_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xy = cbuffer.data(id_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xz = cbuffer.data(id_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yy = cbuffer.data(id_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yz = cbuffer.data(id_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_yyzzzz_zz = cbuffer.data(id_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xx = cbuffer.data(id_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xy = cbuffer.data(id_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xz = cbuffer.data(id_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yy = cbuffer.data(id_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yz = cbuffer.data(id_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_yzzzzz_zz = cbuffer.data(id_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xx = cbuffer.data(id_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xy = cbuffer.data(id_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xz = cbuffer.data(id_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yy = cbuffer.data(id_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yz = cbuffer.data(id_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_zzzzzz_zz = cbuffer.data(id_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xx = cbuffer.data(id_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xy = cbuffer.data(id_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xz = cbuffer.data(id_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yy = cbuffer.data(id_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yz = cbuffer.data(id_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_z_xxxxxx_zz = cbuffer.data(id_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xx = cbuffer.data(id_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xy = cbuffer.data(id_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xz = cbuffer.data(id_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yy = cbuffer.data(id_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yz = cbuffer.data(id_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_z_xxxxxy_zz = cbuffer.data(id_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xx = cbuffer.data(id_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xy = cbuffer.data(id_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xz = cbuffer.data(id_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yy = cbuffer.data(id_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yz = cbuffer.data(id_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_z_xxxxxz_zz = cbuffer.data(id_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xx = cbuffer.data(id_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xy = cbuffer.data(id_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xz = cbuffer.data(id_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yy = cbuffer.data(id_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yz = cbuffer.data(id_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_z_xxxxyy_zz = cbuffer.data(id_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xx = cbuffer.data(id_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xy = cbuffer.data(id_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xz = cbuffer.data(id_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yy = cbuffer.data(id_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yz = cbuffer.data(id_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_z_xxxxyz_zz = cbuffer.data(id_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xx = cbuffer.data(id_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xy = cbuffer.data(id_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xz = cbuffer.data(id_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yy = cbuffer.data(id_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yz = cbuffer.data(id_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_z_xxxxzz_zz = cbuffer.data(id_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xx = cbuffer.data(id_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xy = cbuffer.data(id_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xz = cbuffer.data(id_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yy = cbuffer.data(id_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yz = cbuffer.data(id_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_z_xxxyyy_zz = cbuffer.data(id_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xx = cbuffer.data(id_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xy = cbuffer.data(id_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xz = cbuffer.data(id_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yy = cbuffer.data(id_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yz = cbuffer.data(id_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_z_xxxyyz_zz = cbuffer.data(id_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xx = cbuffer.data(id_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xy = cbuffer.data(id_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xz = cbuffer.data(id_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yy = cbuffer.data(id_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yz = cbuffer.data(id_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_z_xxxyzz_zz = cbuffer.data(id_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xx = cbuffer.data(id_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xy = cbuffer.data(id_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xz = cbuffer.data(id_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yy = cbuffer.data(id_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yz = cbuffer.data(id_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_z_xxxzzz_zz = cbuffer.data(id_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xx = cbuffer.data(id_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xy = cbuffer.data(id_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xz = cbuffer.data(id_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yy = cbuffer.data(id_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yz = cbuffer.data(id_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_z_xxyyyy_zz = cbuffer.data(id_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xx = cbuffer.data(id_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xy = cbuffer.data(id_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xz = cbuffer.data(id_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yy = cbuffer.data(id_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yz = cbuffer.data(id_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_z_xxyyyz_zz = cbuffer.data(id_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xx = cbuffer.data(id_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xy = cbuffer.data(id_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xz = cbuffer.data(id_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yy = cbuffer.data(id_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yz = cbuffer.data(id_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_z_xxyyzz_zz = cbuffer.data(id_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xx = cbuffer.data(id_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xy = cbuffer.data(id_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xz = cbuffer.data(id_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yy = cbuffer.data(id_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yz = cbuffer.data(id_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_z_xxyzzz_zz = cbuffer.data(id_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xx = cbuffer.data(id_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xy = cbuffer.data(id_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xz = cbuffer.data(id_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yy = cbuffer.data(id_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yz = cbuffer.data(id_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_z_xxzzzz_zz = cbuffer.data(id_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xx = cbuffer.data(id_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xy = cbuffer.data(id_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xz = cbuffer.data(id_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yy = cbuffer.data(id_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yz = cbuffer.data(id_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_z_xyyyyy_zz = cbuffer.data(id_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xx = cbuffer.data(id_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xy = cbuffer.data(id_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xz = cbuffer.data(id_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yy = cbuffer.data(id_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yz = cbuffer.data(id_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_z_xyyyyz_zz = cbuffer.data(id_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xx = cbuffer.data(id_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xy = cbuffer.data(id_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xz = cbuffer.data(id_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yy = cbuffer.data(id_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yz = cbuffer.data(id_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_z_xyyyzz_zz = cbuffer.data(id_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xx = cbuffer.data(id_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xy = cbuffer.data(id_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xz = cbuffer.data(id_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yy = cbuffer.data(id_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yz = cbuffer.data(id_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_z_xyyzzz_zz = cbuffer.data(id_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xx = cbuffer.data(id_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xy = cbuffer.data(id_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xz = cbuffer.data(id_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yy = cbuffer.data(id_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yz = cbuffer.data(id_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_z_xyzzzz_zz = cbuffer.data(id_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xx = cbuffer.data(id_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xy = cbuffer.data(id_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xz = cbuffer.data(id_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yy = cbuffer.data(id_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yz = cbuffer.data(id_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_z_xzzzzz_zz = cbuffer.data(id_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xx = cbuffer.data(id_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xy = cbuffer.data(id_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xz = cbuffer.data(id_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yy = cbuffer.data(id_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yz = cbuffer.data(id_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_z_yyyyyy_zz = cbuffer.data(id_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xx = cbuffer.data(id_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xy = cbuffer.data(id_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xz = cbuffer.data(id_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yy = cbuffer.data(id_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yz = cbuffer.data(id_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_z_yyyyyz_zz = cbuffer.data(id_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xx = cbuffer.data(id_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xy = cbuffer.data(id_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xz = cbuffer.data(id_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yy = cbuffer.data(id_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yz = cbuffer.data(id_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_z_yyyyzz_zz = cbuffer.data(id_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xx = cbuffer.data(id_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xy = cbuffer.data(id_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xz = cbuffer.data(id_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yy = cbuffer.data(id_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yz = cbuffer.data(id_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_z_yyyzzz_zz = cbuffer.data(id_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xx = cbuffer.data(id_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xy = cbuffer.data(id_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xz = cbuffer.data(id_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yy = cbuffer.data(id_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yz = cbuffer.data(id_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_z_yyzzzz_zz = cbuffer.data(id_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xx = cbuffer.data(id_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xy = cbuffer.data(id_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xz = cbuffer.data(id_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yy = cbuffer.data(id_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yz = cbuffer.data(id_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_z_yzzzzz_zz = cbuffer.data(id_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xx = cbuffer.data(id_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xy = cbuffer.data(id_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xz = cbuffer.data(id_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yy = cbuffer.data(id_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yz = cbuffer.data(id_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_z_zzzzzz_zz = cbuffer.data(id_geom_01_off + 503 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_kpxx

            const auto kp_geom_01_off = idx_geom_01_kpxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxx_x = cbuffer.data(kp_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_y = cbuffer.data(kp_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_z = cbuffer.data(kp_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_x, g_0_x_xxxxxx_xx, g_0_x_xxxxxx_xy, g_0_x_xxxxxx_xz, g_0_x_xxxxxx_y, g_0_x_xxxxxx_z, g_0_x_xxxxxxx_x, g_0_x_xxxxxxx_y, g_0_x_xxxxxxx_z, g_xxxxxx_x, g_xxxxxx_y, g_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxx_x[k] = g_xxxxxx_x[k] - g_0_x_xxxxxx_x[k] * ab_x + g_0_x_xxxxxx_xx[k];

                g_0_x_xxxxxxx_y[k] = g_xxxxxx_y[k] - g_0_x_xxxxxx_y[k] * ab_x + g_0_x_xxxxxx_xy[k];

                g_0_x_xxxxxxx_z[k] = g_xxxxxx_z[k] - g_0_x_xxxxxx_z[k] * ab_x + g_0_x_xxxxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxy_x = cbuffer.data(kp_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_y = cbuffer.data(kp_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_z = cbuffer.data(kp_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_x, g_0_x_xxxxxx_xy, g_0_x_xxxxxx_y, g_0_x_xxxxxx_yy, g_0_x_xxxxxx_yz, g_0_x_xxxxxx_z, g_0_x_xxxxxxy_x, g_0_x_xxxxxxy_y, g_0_x_xxxxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxy_x[k] = -g_0_x_xxxxxx_x[k] * ab_y + g_0_x_xxxxxx_xy[k];

                g_0_x_xxxxxxy_y[k] = -g_0_x_xxxxxx_y[k] * ab_y + g_0_x_xxxxxx_yy[k];

                g_0_x_xxxxxxy_z[k] = -g_0_x_xxxxxx_z[k] * ab_y + g_0_x_xxxxxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxz_x = cbuffer.data(kp_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_y = cbuffer.data(kp_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_z = cbuffer.data(kp_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_x, g_0_x_xxxxxx_xz, g_0_x_xxxxxx_y, g_0_x_xxxxxx_yz, g_0_x_xxxxxx_z, g_0_x_xxxxxx_zz, g_0_x_xxxxxxz_x, g_0_x_xxxxxxz_y, g_0_x_xxxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxz_x[k] = -g_0_x_xxxxxx_x[k] * ab_z + g_0_x_xxxxxx_xz[k];

                g_0_x_xxxxxxz_y[k] = -g_0_x_xxxxxx_y[k] * ab_z + g_0_x_xxxxxx_yz[k];

                g_0_x_xxxxxxz_z[k] = -g_0_x_xxxxxx_z[k] * ab_z + g_0_x_xxxxxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxyy_x = cbuffer.data(kp_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_y = cbuffer.data(kp_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_z = cbuffer.data(kp_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxy_x, g_0_x_xxxxxy_xy, g_0_x_xxxxxy_y, g_0_x_xxxxxy_yy, g_0_x_xxxxxy_yz, g_0_x_xxxxxy_z, g_0_x_xxxxxyy_x, g_0_x_xxxxxyy_y, g_0_x_xxxxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxyy_x[k] = -g_0_x_xxxxxy_x[k] * ab_y + g_0_x_xxxxxy_xy[k];

                g_0_x_xxxxxyy_y[k] = -g_0_x_xxxxxy_y[k] * ab_y + g_0_x_xxxxxy_yy[k];

                g_0_x_xxxxxyy_z[k] = -g_0_x_xxxxxy_z[k] * ab_y + g_0_x_xxxxxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxyz_x = cbuffer.data(kp_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_y = cbuffer.data(kp_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_z = cbuffer.data(kp_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxyz_x, g_0_x_xxxxxyz_y, g_0_x_xxxxxyz_z, g_0_x_xxxxxz_x, g_0_x_xxxxxz_xy, g_0_x_xxxxxz_y, g_0_x_xxxxxz_yy, g_0_x_xxxxxz_yz, g_0_x_xxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxyz_x[k] = -g_0_x_xxxxxz_x[k] * ab_y + g_0_x_xxxxxz_xy[k];

                g_0_x_xxxxxyz_y[k] = -g_0_x_xxxxxz_y[k] * ab_y + g_0_x_xxxxxz_yy[k];

                g_0_x_xxxxxyz_z[k] = -g_0_x_xxxxxz_z[k] * ab_y + g_0_x_xxxxxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxzz_x = cbuffer.data(kp_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_y = cbuffer.data(kp_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_z = cbuffer.data(kp_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxz_x, g_0_x_xxxxxz_xz, g_0_x_xxxxxz_y, g_0_x_xxxxxz_yz, g_0_x_xxxxxz_z, g_0_x_xxxxxz_zz, g_0_x_xxxxxzz_x, g_0_x_xxxxxzz_y, g_0_x_xxxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxzz_x[k] = -g_0_x_xxxxxz_x[k] * ab_z + g_0_x_xxxxxz_xz[k];

                g_0_x_xxxxxzz_y[k] = -g_0_x_xxxxxz_y[k] * ab_z + g_0_x_xxxxxz_yz[k];

                g_0_x_xxxxxzz_z[k] = -g_0_x_xxxxxz_z[k] * ab_z + g_0_x_xxxxxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyyy_x = cbuffer.data(kp_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_y = cbuffer.data(kp_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_z = cbuffer.data(kp_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyy_x, g_0_x_xxxxyy_xy, g_0_x_xxxxyy_y, g_0_x_xxxxyy_yy, g_0_x_xxxxyy_yz, g_0_x_xxxxyy_z, g_0_x_xxxxyyy_x, g_0_x_xxxxyyy_y, g_0_x_xxxxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyyy_x[k] = -g_0_x_xxxxyy_x[k] * ab_y + g_0_x_xxxxyy_xy[k];

                g_0_x_xxxxyyy_y[k] = -g_0_x_xxxxyy_y[k] * ab_y + g_0_x_xxxxyy_yy[k];

                g_0_x_xxxxyyy_z[k] = -g_0_x_xxxxyy_z[k] * ab_y + g_0_x_xxxxyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyyz_x = cbuffer.data(kp_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_y = cbuffer.data(kp_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_z = cbuffer.data(kp_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyyz_x, g_0_x_xxxxyyz_y, g_0_x_xxxxyyz_z, g_0_x_xxxxyz_x, g_0_x_xxxxyz_xy, g_0_x_xxxxyz_y, g_0_x_xxxxyz_yy, g_0_x_xxxxyz_yz, g_0_x_xxxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyyz_x[k] = -g_0_x_xxxxyz_x[k] * ab_y + g_0_x_xxxxyz_xy[k];

                g_0_x_xxxxyyz_y[k] = -g_0_x_xxxxyz_y[k] * ab_y + g_0_x_xxxxyz_yy[k];

                g_0_x_xxxxyyz_z[k] = -g_0_x_xxxxyz_z[k] * ab_y + g_0_x_xxxxyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyzz_x = cbuffer.data(kp_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_y = cbuffer.data(kp_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_z = cbuffer.data(kp_geom_01_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyzz_x, g_0_x_xxxxyzz_y, g_0_x_xxxxyzz_z, g_0_x_xxxxzz_x, g_0_x_xxxxzz_xy, g_0_x_xxxxzz_y, g_0_x_xxxxzz_yy, g_0_x_xxxxzz_yz, g_0_x_xxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyzz_x[k] = -g_0_x_xxxxzz_x[k] * ab_y + g_0_x_xxxxzz_xy[k];

                g_0_x_xxxxyzz_y[k] = -g_0_x_xxxxzz_y[k] * ab_y + g_0_x_xxxxzz_yy[k];

                g_0_x_xxxxyzz_z[k] = -g_0_x_xxxxzz_z[k] * ab_y + g_0_x_xxxxzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxzzz_x = cbuffer.data(kp_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_y = cbuffer.data(kp_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_z = cbuffer.data(kp_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxzz_x, g_0_x_xxxxzz_xz, g_0_x_xxxxzz_y, g_0_x_xxxxzz_yz, g_0_x_xxxxzz_z, g_0_x_xxxxzz_zz, g_0_x_xxxxzzz_x, g_0_x_xxxxzzz_y, g_0_x_xxxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxzzz_x[k] = -g_0_x_xxxxzz_x[k] * ab_z + g_0_x_xxxxzz_xz[k];

                g_0_x_xxxxzzz_y[k] = -g_0_x_xxxxzz_y[k] * ab_z + g_0_x_xxxxzz_yz[k];

                g_0_x_xxxxzzz_z[k] = -g_0_x_xxxxzz_z[k] * ab_z + g_0_x_xxxxzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyyy_x = cbuffer.data(kp_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_y = cbuffer.data(kp_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_z = cbuffer.data(kp_geom_01_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyy_x, g_0_x_xxxyyy_xy, g_0_x_xxxyyy_y, g_0_x_xxxyyy_yy, g_0_x_xxxyyy_yz, g_0_x_xxxyyy_z, g_0_x_xxxyyyy_x, g_0_x_xxxyyyy_y, g_0_x_xxxyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyyy_x[k] = -g_0_x_xxxyyy_x[k] * ab_y + g_0_x_xxxyyy_xy[k];

                g_0_x_xxxyyyy_y[k] = -g_0_x_xxxyyy_y[k] * ab_y + g_0_x_xxxyyy_yy[k];

                g_0_x_xxxyyyy_z[k] = -g_0_x_xxxyyy_z[k] * ab_y + g_0_x_xxxyyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyyz_x = cbuffer.data(kp_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_y = cbuffer.data(kp_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_z = cbuffer.data(kp_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyyz_x, g_0_x_xxxyyyz_y, g_0_x_xxxyyyz_z, g_0_x_xxxyyz_x, g_0_x_xxxyyz_xy, g_0_x_xxxyyz_y, g_0_x_xxxyyz_yy, g_0_x_xxxyyz_yz, g_0_x_xxxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyyz_x[k] = -g_0_x_xxxyyz_x[k] * ab_y + g_0_x_xxxyyz_xy[k];

                g_0_x_xxxyyyz_y[k] = -g_0_x_xxxyyz_y[k] * ab_y + g_0_x_xxxyyz_yy[k];

                g_0_x_xxxyyyz_z[k] = -g_0_x_xxxyyz_z[k] * ab_y + g_0_x_xxxyyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyzz_x = cbuffer.data(kp_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_y = cbuffer.data(kp_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_z = cbuffer.data(kp_geom_01_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyzz_x, g_0_x_xxxyyzz_y, g_0_x_xxxyyzz_z, g_0_x_xxxyzz_x, g_0_x_xxxyzz_xy, g_0_x_xxxyzz_y, g_0_x_xxxyzz_yy, g_0_x_xxxyzz_yz, g_0_x_xxxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyzz_x[k] = -g_0_x_xxxyzz_x[k] * ab_y + g_0_x_xxxyzz_xy[k];

                g_0_x_xxxyyzz_y[k] = -g_0_x_xxxyzz_y[k] * ab_y + g_0_x_xxxyzz_yy[k];

                g_0_x_xxxyyzz_z[k] = -g_0_x_xxxyzz_z[k] * ab_y + g_0_x_xxxyzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyzzz_x = cbuffer.data(kp_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_y = cbuffer.data(kp_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_z = cbuffer.data(kp_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyzzz_x, g_0_x_xxxyzzz_y, g_0_x_xxxyzzz_z, g_0_x_xxxzzz_x, g_0_x_xxxzzz_xy, g_0_x_xxxzzz_y, g_0_x_xxxzzz_yy, g_0_x_xxxzzz_yz, g_0_x_xxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyzzz_x[k] = -g_0_x_xxxzzz_x[k] * ab_y + g_0_x_xxxzzz_xy[k];

                g_0_x_xxxyzzz_y[k] = -g_0_x_xxxzzz_y[k] * ab_y + g_0_x_xxxzzz_yy[k];

                g_0_x_xxxyzzz_z[k] = -g_0_x_xxxzzz_z[k] * ab_y + g_0_x_xxxzzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzzzz_x = cbuffer.data(kp_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_y = cbuffer.data(kp_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_z = cbuffer.data(kp_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxzzz_x, g_0_x_xxxzzz_xz, g_0_x_xxxzzz_y, g_0_x_xxxzzz_yz, g_0_x_xxxzzz_z, g_0_x_xxxzzz_zz, g_0_x_xxxzzzz_x, g_0_x_xxxzzzz_y, g_0_x_xxxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzzzz_x[k] = -g_0_x_xxxzzz_x[k] * ab_z + g_0_x_xxxzzz_xz[k];

                g_0_x_xxxzzzz_y[k] = -g_0_x_xxxzzz_y[k] * ab_z + g_0_x_xxxzzz_yz[k];

                g_0_x_xxxzzzz_z[k] = -g_0_x_xxxzzz_z[k] * ab_z + g_0_x_xxxzzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyyy_x = cbuffer.data(kp_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_y = cbuffer.data(kp_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_z = cbuffer.data(kp_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyy_x, g_0_x_xxyyyy_xy, g_0_x_xxyyyy_y, g_0_x_xxyyyy_yy, g_0_x_xxyyyy_yz, g_0_x_xxyyyy_z, g_0_x_xxyyyyy_x, g_0_x_xxyyyyy_y, g_0_x_xxyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyyy_x[k] = -g_0_x_xxyyyy_x[k] * ab_y + g_0_x_xxyyyy_xy[k];

                g_0_x_xxyyyyy_y[k] = -g_0_x_xxyyyy_y[k] * ab_y + g_0_x_xxyyyy_yy[k];

                g_0_x_xxyyyyy_z[k] = -g_0_x_xxyyyy_z[k] * ab_y + g_0_x_xxyyyy_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyyz_x = cbuffer.data(kp_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_y = cbuffer.data(kp_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_z = cbuffer.data(kp_geom_01_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyyz_x, g_0_x_xxyyyyz_y, g_0_x_xxyyyyz_z, g_0_x_xxyyyz_x, g_0_x_xxyyyz_xy, g_0_x_xxyyyz_y, g_0_x_xxyyyz_yy, g_0_x_xxyyyz_yz, g_0_x_xxyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyyz_x[k] = -g_0_x_xxyyyz_x[k] * ab_y + g_0_x_xxyyyz_xy[k];

                g_0_x_xxyyyyz_y[k] = -g_0_x_xxyyyz_y[k] * ab_y + g_0_x_xxyyyz_yy[k];

                g_0_x_xxyyyyz_z[k] = -g_0_x_xxyyyz_z[k] * ab_y + g_0_x_xxyyyz_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyzz_x = cbuffer.data(kp_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_y = cbuffer.data(kp_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_z = cbuffer.data(kp_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyzz_x, g_0_x_xxyyyzz_y, g_0_x_xxyyyzz_z, g_0_x_xxyyzz_x, g_0_x_xxyyzz_xy, g_0_x_xxyyzz_y, g_0_x_xxyyzz_yy, g_0_x_xxyyzz_yz, g_0_x_xxyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyzz_x[k] = -g_0_x_xxyyzz_x[k] * ab_y + g_0_x_xxyyzz_xy[k];

                g_0_x_xxyyyzz_y[k] = -g_0_x_xxyyzz_y[k] * ab_y + g_0_x_xxyyzz_yy[k];

                g_0_x_xxyyyzz_z[k] = -g_0_x_xxyyzz_z[k] * ab_y + g_0_x_xxyyzz_yz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyzzz_x = cbuffer.data(kp_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_y = cbuffer.data(kp_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_z = cbuffer.data(kp_geom_01_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyzzz_x, g_0_x_xxyyzzz_y, g_0_x_xxyyzzz_z, g_0_x_xxyzzz_x, g_0_x_xxyzzz_xy, g_0_x_xxyzzz_y, g_0_x_xxyzzz_yy, g_0_x_xxyzzz_yz, g_0_x_xxyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyzzz_x[k] = -g_0_x_xxyzzz_x[k] * ab_y + g_0_x_xxyzzz_xy[k];

                g_0_x_xxyyzzz_y[k] = -g_0_x_xxyzzz_y[k] * ab_y + g_0_x_xxyzzz_yy[k];

                g_0_x_xxyyzzz_z[k] = -g_0_x_xxyzzz_z[k] * ab_y + g_0_x_xxyzzz_yz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzzzz_x = cbuffer.data(kp_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_y = cbuffer.data(kp_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_z = cbuffer.data(kp_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzzzz_x, g_0_x_xxyzzzz_y, g_0_x_xxyzzzz_z, g_0_x_xxzzzz_x, g_0_x_xxzzzz_xy, g_0_x_xxzzzz_y, g_0_x_xxzzzz_yy, g_0_x_xxzzzz_yz, g_0_x_xxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzzzz_x[k] = -g_0_x_xxzzzz_x[k] * ab_y + g_0_x_xxzzzz_xy[k];

                g_0_x_xxyzzzz_y[k] = -g_0_x_xxzzzz_y[k] * ab_y + g_0_x_xxzzzz_yy[k];

                g_0_x_xxyzzzz_z[k] = -g_0_x_xxzzzz_z[k] * ab_y + g_0_x_xxzzzz_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzzzz_x = cbuffer.data(kp_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_y = cbuffer.data(kp_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_z = cbuffer.data(kp_geom_01_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzzzz_x, g_0_x_xxzzzz_xz, g_0_x_xxzzzz_y, g_0_x_xxzzzz_yz, g_0_x_xxzzzz_z, g_0_x_xxzzzz_zz, g_0_x_xxzzzzz_x, g_0_x_xxzzzzz_y, g_0_x_xxzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzzzz_x[k] = -g_0_x_xxzzzz_x[k] * ab_z + g_0_x_xxzzzz_xz[k];

                g_0_x_xxzzzzz_y[k] = -g_0_x_xxzzzz_y[k] * ab_z + g_0_x_xxzzzz_yz[k];

                g_0_x_xxzzzzz_z[k] = -g_0_x_xxzzzz_z[k] * ab_z + g_0_x_xxzzzz_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyyy_x = cbuffer.data(kp_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_y = cbuffer.data(kp_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_z = cbuffer.data(kp_geom_01_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyy_x, g_0_x_xyyyyy_xy, g_0_x_xyyyyy_y, g_0_x_xyyyyy_yy, g_0_x_xyyyyy_yz, g_0_x_xyyyyy_z, g_0_x_xyyyyyy_x, g_0_x_xyyyyyy_y, g_0_x_xyyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyyy_x[k] = -g_0_x_xyyyyy_x[k] * ab_y + g_0_x_xyyyyy_xy[k];

                g_0_x_xyyyyyy_y[k] = -g_0_x_xyyyyy_y[k] * ab_y + g_0_x_xyyyyy_yy[k];

                g_0_x_xyyyyyy_z[k] = -g_0_x_xyyyyy_z[k] * ab_y + g_0_x_xyyyyy_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyyz_x = cbuffer.data(kp_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_y = cbuffer.data(kp_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_z = cbuffer.data(kp_geom_01_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyyz_x, g_0_x_xyyyyyz_y, g_0_x_xyyyyyz_z, g_0_x_xyyyyz_x, g_0_x_xyyyyz_xy, g_0_x_xyyyyz_y, g_0_x_xyyyyz_yy, g_0_x_xyyyyz_yz, g_0_x_xyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyyz_x[k] = -g_0_x_xyyyyz_x[k] * ab_y + g_0_x_xyyyyz_xy[k];

                g_0_x_xyyyyyz_y[k] = -g_0_x_xyyyyz_y[k] * ab_y + g_0_x_xyyyyz_yy[k];

                g_0_x_xyyyyyz_z[k] = -g_0_x_xyyyyz_z[k] * ab_y + g_0_x_xyyyyz_yz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyzz_x = cbuffer.data(kp_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_y = cbuffer.data(kp_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_z = cbuffer.data(kp_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyzz_x, g_0_x_xyyyyzz_y, g_0_x_xyyyyzz_z, g_0_x_xyyyzz_x, g_0_x_xyyyzz_xy, g_0_x_xyyyzz_y, g_0_x_xyyyzz_yy, g_0_x_xyyyzz_yz, g_0_x_xyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyzz_x[k] = -g_0_x_xyyyzz_x[k] * ab_y + g_0_x_xyyyzz_xy[k];

                g_0_x_xyyyyzz_y[k] = -g_0_x_xyyyzz_y[k] * ab_y + g_0_x_xyyyzz_yy[k];

                g_0_x_xyyyyzz_z[k] = -g_0_x_xyyyzz_z[k] * ab_y + g_0_x_xyyyzz_yz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyzzz_x = cbuffer.data(kp_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_y = cbuffer.data(kp_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_z = cbuffer.data(kp_geom_01_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyzzz_x, g_0_x_xyyyzzz_y, g_0_x_xyyyzzz_z, g_0_x_xyyzzz_x, g_0_x_xyyzzz_xy, g_0_x_xyyzzz_y, g_0_x_xyyzzz_yy, g_0_x_xyyzzz_yz, g_0_x_xyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyzzz_x[k] = -g_0_x_xyyzzz_x[k] * ab_y + g_0_x_xyyzzz_xy[k];

                g_0_x_xyyyzzz_y[k] = -g_0_x_xyyzzz_y[k] * ab_y + g_0_x_xyyzzz_yy[k];

                g_0_x_xyyyzzz_z[k] = -g_0_x_xyyzzz_z[k] * ab_y + g_0_x_xyyzzz_yz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzzzz_x = cbuffer.data(kp_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_y = cbuffer.data(kp_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_z = cbuffer.data(kp_geom_01_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzzzz_x, g_0_x_xyyzzzz_y, g_0_x_xyyzzzz_z, g_0_x_xyzzzz_x, g_0_x_xyzzzz_xy, g_0_x_xyzzzz_y, g_0_x_xyzzzz_yy, g_0_x_xyzzzz_yz, g_0_x_xyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzzzz_x[k] = -g_0_x_xyzzzz_x[k] * ab_y + g_0_x_xyzzzz_xy[k];

                g_0_x_xyyzzzz_y[k] = -g_0_x_xyzzzz_y[k] * ab_y + g_0_x_xyzzzz_yy[k];

                g_0_x_xyyzzzz_z[k] = -g_0_x_xyzzzz_z[k] * ab_y + g_0_x_xyzzzz_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzzzz_x = cbuffer.data(kp_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_y = cbuffer.data(kp_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_z = cbuffer.data(kp_geom_01_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzzzz_x, g_0_x_xyzzzzz_y, g_0_x_xyzzzzz_z, g_0_x_xzzzzz_x, g_0_x_xzzzzz_xy, g_0_x_xzzzzz_y, g_0_x_xzzzzz_yy, g_0_x_xzzzzz_yz, g_0_x_xzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzzzz_x[k] = -g_0_x_xzzzzz_x[k] * ab_y + g_0_x_xzzzzz_xy[k];

                g_0_x_xyzzzzz_y[k] = -g_0_x_xzzzzz_y[k] * ab_y + g_0_x_xzzzzz_yy[k];

                g_0_x_xyzzzzz_z[k] = -g_0_x_xzzzzz_z[k] * ab_y + g_0_x_xzzzzz_yz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzzzz_x = cbuffer.data(kp_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_y = cbuffer.data(kp_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_z = cbuffer.data(kp_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzzzz_x, g_0_x_xzzzzz_xz, g_0_x_xzzzzz_y, g_0_x_xzzzzz_yz, g_0_x_xzzzzz_z, g_0_x_xzzzzz_zz, g_0_x_xzzzzzz_x, g_0_x_xzzzzzz_y, g_0_x_xzzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzzzz_x[k] = -g_0_x_xzzzzz_x[k] * ab_z + g_0_x_xzzzzz_xz[k];

                g_0_x_xzzzzzz_y[k] = -g_0_x_xzzzzz_y[k] * ab_z + g_0_x_xzzzzz_yz[k];

                g_0_x_xzzzzzz_z[k] = -g_0_x_xzzzzz_z[k] * ab_z + g_0_x_xzzzzz_zz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyyy_x = cbuffer.data(kp_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_y = cbuffer.data(kp_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_z = cbuffer.data(kp_geom_01_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyy_x, g_0_x_yyyyyy_xy, g_0_x_yyyyyy_y, g_0_x_yyyyyy_yy, g_0_x_yyyyyy_yz, g_0_x_yyyyyy_z, g_0_x_yyyyyyy_x, g_0_x_yyyyyyy_y, g_0_x_yyyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyyy_x[k] = -g_0_x_yyyyyy_x[k] * ab_y + g_0_x_yyyyyy_xy[k];

                g_0_x_yyyyyyy_y[k] = -g_0_x_yyyyyy_y[k] * ab_y + g_0_x_yyyyyy_yy[k];

                g_0_x_yyyyyyy_z[k] = -g_0_x_yyyyyy_z[k] * ab_y + g_0_x_yyyyyy_yz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyyz_x = cbuffer.data(kp_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_y = cbuffer.data(kp_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_z = cbuffer.data(kp_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyyz_x, g_0_x_yyyyyyz_y, g_0_x_yyyyyyz_z, g_0_x_yyyyyz_x, g_0_x_yyyyyz_xy, g_0_x_yyyyyz_y, g_0_x_yyyyyz_yy, g_0_x_yyyyyz_yz, g_0_x_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyyz_x[k] = -g_0_x_yyyyyz_x[k] * ab_y + g_0_x_yyyyyz_xy[k];

                g_0_x_yyyyyyz_y[k] = -g_0_x_yyyyyz_y[k] * ab_y + g_0_x_yyyyyz_yy[k];

                g_0_x_yyyyyyz_z[k] = -g_0_x_yyyyyz_z[k] * ab_y + g_0_x_yyyyyz_yz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyzz_x = cbuffer.data(kp_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_y = cbuffer.data(kp_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_z = cbuffer.data(kp_geom_01_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyzz_x, g_0_x_yyyyyzz_y, g_0_x_yyyyyzz_z, g_0_x_yyyyzz_x, g_0_x_yyyyzz_xy, g_0_x_yyyyzz_y, g_0_x_yyyyzz_yy, g_0_x_yyyyzz_yz, g_0_x_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyzz_x[k] = -g_0_x_yyyyzz_x[k] * ab_y + g_0_x_yyyyzz_xy[k];

                g_0_x_yyyyyzz_y[k] = -g_0_x_yyyyzz_y[k] * ab_y + g_0_x_yyyyzz_yy[k];

                g_0_x_yyyyyzz_z[k] = -g_0_x_yyyyzz_z[k] * ab_y + g_0_x_yyyyzz_yz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyzzz_x = cbuffer.data(kp_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_y = cbuffer.data(kp_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_z = cbuffer.data(kp_geom_01_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyzzz_x, g_0_x_yyyyzzz_y, g_0_x_yyyyzzz_z, g_0_x_yyyzzz_x, g_0_x_yyyzzz_xy, g_0_x_yyyzzz_y, g_0_x_yyyzzz_yy, g_0_x_yyyzzz_yz, g_0_x_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyzzz_x[k] = -g_0_x_yyyzzz_x[k] * ab_y + g_0_x_yyyzzz_xy[k];

                g_0_x_yyyyzzz_y[k] = -g_0_x_yyyzzz_y[k] * ab_y + g_0_x_yyyzzz_yy[k];

                g_0_x_yyyyzzz_z[k] = -g_0_x_yyyzzz_z[k] * ab_y + g_0_x_yyyzzz_yz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzzzz_x = cbuffer.data(kp_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_y = cbuffer.data(kp_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_z = cbuffer.data(kp_geom_01_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzzzz_x, g_0_x_yyyzzzz_y, g_0_x_yyyzzzz_z, g_0_x_yyzzzz_x, g_0_x_yyzzzz_xy, g_0_x_yyzzzz_y, g_0_x_yyzzzz_yy, g_0_x_yyzzzz_yz, g_0_x_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzzzz_x[k] = -g_0_x_yyzzzz_x[k] * ab_y + g_0_x_yyzzzz_xy[k];

                g_0_x_yyyzzzz_y[k] = -g_0_x_yyzzzz_y[k] * ab_y + g_0_x_yyzzzz_yy[k];

                g_0_x_yyyzzzz_z[k] = -g_0_x_yyzzzz_z[k] * ab_y + g_0_x_yyzzzz_yz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzzzz_x = cbuffer.data(kp_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_y = cbuffer.data(kp_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_z = cbuffer.data(kp_geom_01_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzzzz_x, g_0_x_yyzzzzz_y, g_0_x_yyzzzzz_z, g_0_x_yzzzzz_x, g_0_x_yzzzzz_xy, g_0_x_yzzzzz_y, g_0_x_yzzzzz_yy, g_0_x_yzzzzz_yz, g_0_x_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzzzz_x[k] = -g_0_x_yzzzzz_x[k] * ab_y + g_0_x_yzzzzz_xy[k];

                g_0_x_yyzzzzz_y[k] = -g_0_x_yzzzzz_y[k] * ab_y + g_0_x_yzzzzz_yy[k];

                g_0_x_yyzzzzz_z[k] = -g_0_x_yzzzzz_z[k] * ab_y + g_0_x_yzzzzz_yz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzzzz_x = cbuffer.data(kp_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_y = cbuffer.data(kp_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_z = cbuffer.data(kp_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzzzz_x, g_0_x_yzzzzzz_y, g_0_x_yzzzzzz_z, g_0_x_zzzzzz_x, g_0_x_zzzzzz_xy, g_0_x_zzzzzz_y, g_0_x_zzzzzz_yy, g_0_x_zzzzzz_yz, g_0_x_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzzzz_x[k] = -g_0_x_zzzzzz_x[k] * ab_y + g_0_x_zzzzzz_xy[k];

                g_0_x_yzzzzzz_y[k] = -g_0_x_zzzzzz_y[k] * ab_y + g_0_x_zzzzzz_yy[k];

                g_0_x_yzzzzzz_z[k] = -g_0_x_zzzzzz_z[k] * ab_y + g_0_x_zzzzzz_yz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzzzz_x = cbuffer.data(kp_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_y = cbuffer.data(kp_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_z = cbuffer.data(kp_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzzzz_x, g_0_x_zzzzzz_xz, g_0_x_zzzzzz_y, g_0_x_zzzzzz_yz, g_0_x_zzzzzz_z, g_0_x_zzzzzz_zz, g_0_x_zzzzzzz_x, g_0_x_zzzzzzz_y, g_0_x_zzzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzzzz_x[k] = -g_0_x_zzzzzz_x[k] * ab_z + g_0_x_zzzzzz_xz[k];

                g_0_x_zzzzzzz_y[k] = -g_0_x_zzzzzz_y[k] * ab_z + g_0_x_zzzzzz_yz[k];

                g_0_x_zzzzzzz_z[k] = -g_0_x_zzzzzz_z[k] * ab_z + g_0_x_zzzzzz_zz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxx_x = cbuffer.data(kp_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_y = cbuffer.data(kp_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_z = cbuffer.data(kp_geom_01_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxx_x, g_0_y_xxxxxx_xx, g_0_y_xxxxxx_xy, g_0_y_xxxxxx_xz, g_0_y_xxxxxx_y, g_0_y_xxxxxx_z, g_0_y_xxxxxxx_x, g_0_y_xxxxxxx_y, g_0_y_xxxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxx_x[k] = -g_0_y_xxxxxx_x[k] * ab_x + g_0_y_xxxxxx_xx[k];

                g_0_y_xxxxxxx_y[k] = -g_0_y_xxxxxx_y[k] * ab_x + g_0_y_xxxxxx_xy[k];

                g_0_y_xxxxxxx_z[k] = -g_0_y_xxxxxx_z[k] * ab_x + g_0_y_xxxxxx_xz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxy_x = cbuffer.data(kp_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_y = cbuffer.data(kp_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_z = cbuffer.data(kp_geom_01_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxxy_x, g_0_y_xxxxxxy_y, g_0_y_xxxxxxy_z, g_0_y_xxxxxy_x, g_0_y_xxxxxy_xx, g_0_y_xxxxxy_xy, g_0_y_xxxxxy_xz, g_0_y_xxxxxy_y, g_0_y_xxxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxy_x[k] = -g_0_y_xxxxxy_x[k] * ab_x + g_0_y_xxxxxy_xx[k];

                g_0_y_xxxxxxy_y[k] = -g_0_y_xxxxxy_y[k] * ab_x + g_0_y_xxxxxy_xy[k];

                g_0_y_xxxxxxy_z[k] = -g_0_y_xxxxxy_z[k] * ab_x + g_0_y_xxxxxy_xz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxz_x = cbuffer.data(kp_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_y = cbuffer.data(kp_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_z = cbuffer.data(kp_geom_01_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxxz_x, g_0_y_xxxxxxz_y, g_0_y_xxxxxxz_z, g_0_y_xxxxxz_x, g_0_y_xxxxxz_xx, g_0_y_xxxxxz_xy, g_0_y_xxxxxz_xz, g_0_y_xxxxxz_y, g_0_y_xxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxz_x[k] = -g_0_y_xxxxxz_x[k] * ab_x + g_0_y_xxxxxz_xx[k];

                g_0_y_xxxxxxz_y[k] = -g_0_y_xxxxxz_y[k] * ab_x + g_0_y_xxxxxz_xy[k];

                g_0_y_xxxxxxz_z[k] = -g_0_y_xxxxxz_z[k] * ab_x + g_0_y_xxxxxz_xz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxyy_x = cbuffer.data(kp_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_y = cbuffer.data(kp_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_z = cbuffer.data(kp_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxyy_x, g_0_y_xxxxxyy_y, g_0_y_xxxxxyy_z, g_0_y_xxxxyy_x, g_0_y_xxxxyy_xx, g_0_y_xxxxyy_xy, g_0_y_xxxxyy_xz, g_0_y_xxxxyy_y, g_0_y_xxxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxyy_x[k] = -g_0_y_xxxxyy_x[k] * ab_x + g_0_y_xxxxyy_xx[k];

                g_0_y_xxxxxyy_y[k] = -g_0_y_xxxxyy_y[k] * ab_x + g_0_y_xxxxyy_xy[k];

                g_0_y_xxxxxyy_z[k] = -g_0_y_xxxxyy_z[k] * ab_x + g_0_y_xxxxyy_xz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxyz_x = cbuffer.data(kp_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_y = cbuffer.data(kp_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_z = cbuffer.data(kp_geom_01_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxyz_x, g_0_y_xxxxxyz_y, g_0_y_xxxxxyz_z, g_0_y_xxxxyz_x, g_0_y_xxxxyz_xx, g_0_y_xxxxyz_xy, g_0_y_xxxxyz_xz, g_0_y_xxxxyz_y, g_0_y_xxxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxyz_x[k] = -g_0_y_xxxxyz_x[k] * ab_x + g_0_y_xxxxyz_xx[k];

                g_0_y_xxxxxyz_y[k] = -g_0_y_xxxxyz_y[k] * ab_x + g_0_y_xxxxyz_xy[k];

                g_0_y_xxxxxyz_z[k] = -g_0_y_xxxxyz_z[k] * ab_x + g_0_y_xxxxyz_xz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxzz_x = cbuffer.data(kp_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_y = cbuffer.data(kp_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_z = cbuffer.data(kp_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxzz_x, g_0_y_xxxxxzz_y, g_0_y_xxxxxzz_z, g_0_y_xxxxzz_x, g_0_y_xxxxzz_xx, g_0_y_xxxxzz_xy, g_0_y_xxxxzz_xz, g_0_y_xxxxzz_y, g_0_y_xxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxzz_x[k] = -g_0_y_xxxxzz_x[k] * ab_x + g_0_y_xxxxzz_xx[k];

                g_0_y_xxxxxzz_y[k] = -g_0_y_xxxxzz_y[k] * ab_x + g_0_y_xxxxzz_xy[k];

                g_0_y_xxxxxzz_z[k] = -g_0_y_xxxxzz_z[k] * ab_x + g_0_y_xxxxzz_xz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyyy_x = cbuffer.data(kp_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_y = cbuffer.data(kp_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_z = cbuffer.data(kp_geom_01_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyyy_x, g_0_y_xxxxyyy_y, g_0_y_xxxxyyy_z, g_0_y_xxxyyy_x, g_0_y_xxxyyy_xx, g_0_y_xxxyyy_xy, g_0_y_xxxyyy_xz, g_0_y_xxxyyy_y, g_0_y_xxxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyyy_x[k] = -g_0_y_xxxyyy_x[k] * ab_x + g_0_y_xxxyyy_xx[k];

                g_0_y_xxxxyyy_y[k] = -g_0_y_xxxyyy_y[k] * ab_x + g_0_y_xxxyyy_xy[k];

                g_0_y_xxxxyyy_z[k] = -g_0_y_xxxyyy_z[k] * ab_x + g_0_y_xxxyyy_xz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyyz_x = cbuffer.data(kp_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_y = cbuffer.data(kp_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_z = cbuffer.data(kp_geom_01_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyyz_x, g_0_y_xxxxyyz_y, g_0_y_xxxxyyz_z, g_0_y_xxxyyz_x, g_0_y_xxxyyz_xx, g_0_y_xxxyyz_xy, g_0_y_xxxyyz_xz, g_0_y_xxxyyz_y, g_0_y_xxxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyyz_x[k] = -g_0_y_xxxyyz_x[k] * ab_x + g_0_y_xxxyyz_xx[k];

                g_0_y_xxxxyyz_y[k] = -g_0_y_xxxyyz_y[k] * ab_x + g_0_y_xxxyyz_xy[k];

                g_0_y_xxxxyyz_z[k] = -g_0_y_xxxyyz_z[k] * ab_x + g_0_y_xxxyyz_xz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyzz_x = cbuffer.data(kp_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_y = cbuffer.data(kp_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_z = cbuffer.data(kp_geom_01_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyzz_x, g_0_y_xxxxyzz_y, g_0_y_xxxxyzz_z, g_0_y_xxxyzz_x, g_0_y_xxxyzz_xx, g_0_y_xxxyzz_xy, g_0_y_xxxyzz_xz, g_0_y_xxxyzz_y, g_0_y_xxxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyzz_x[k] = -g_0_y_xxxyzz_x[k] * ab_x + g_0_y_xxxyzz_xx[k];

                g_0_y_xxxxyzz_y[k] = -g_0_y_xxxyzz_y[k] * ab_x + g_0_y_xxxyzz_xy[k];

                g_0_y_xxxxyzz_z[k] = -g_0_y_xxxyzz_z[k] * ab_x + g_0_y_xxxyzz_xz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxzzz_x = cbuffer.data(kp_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_y = cbuffer.data(kp_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_z = cbuffer.data(kp_geom_01_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxzzz_x, g_0_y_xxxxzzz_y, g_0_y_xxxxzzz_z, g_0_y_xxxzzz_x, g_0_y_xxxzzz_xx, g_0_y_xxxzzz_xy, g_0_y_xxxzzz_xz, g_0_y_xxxzzz_y, g_0_y_xxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxzzz_x[k] = -g_0_y_xxxzzz_x[k] * ab_x + g_0_y_xxxzzz_xx[k];

                g_0_y_xxxxzzz_y[k] = -g_0_y_xxxzzz_y[k] * ab_x + g_0_y_xxxzzz_xy[k];

                g_0_y_xxxxzzz_z[k] = -g_0_y_xxxzzz_z[k] * ab_x + g_0_y_xxxzzz_xz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyyy_x = cbuffer.data(kp_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_y = cbuffer.data(kp_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_z = cbuffer.data(kp_geom_01_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyyy_x, g_0_y_xxxyyyy_y, g_0_y_xxxyyyy_z, g_0_y_xxyyyy_x, g_0_y_xxyyyy_xx, g_0_y_xxyyyy_xy, g_0_y_xxyyyy_xz, g_0_y_xxyyyy_y, g_0_y_xxyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyyy_x[k] = -g_0_y_xxyyyy_x[k] * ab_x + g_0_y_xxyyyy_xx[k];

                g_0_y_xxxyyyy_y[k] = -g_0_y_xxyyyy_y[k] * ab_x + g_0_y_xxyyyy_xy[k];

                g_0_y_xxxyyyy_z[k] = -g_0_y_xxyyyy_z[k] * ab_x + g_0_y_xxyyyy_xz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyyz_x = cbuffer.data(kp_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_y = cbuffer.data(kp_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_z = cbuffer.data(kp_geom_01_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyyz_x, g_0_y_xxxyyyz_y, g_0_y_xxxyyyz_z, g_0_y_xxyyyz_x, g_0_y_xxyyyz_xx, g_0_y_xxyyyz_xy, g_0_y_xxyyyz_xz, g_0_y_xxyyyz_y, g_0_y_xxyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyyz_x[k] = -g_0_y_xxyyyz_x[k] * ab_x + g_0_y_xxyyyz_xx[k];

                g_0_y_xxxyyyz_y[k] = -g_0_y_xxyyyz_y[k] * ab_x + g_0_y_xxyyyz_xy[k];

                g_0_y_xxxyyyz_z[k] = -g_0_y_xxyyyz_z[k] * ab_x + g_0_y_xxyyyz_xz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyzz_x = cbuffer.data(kp_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_y = cbuffer.data(kp_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_z = cbuffer.data(kp_geom_01_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyzz_x, g_0_y_xxxyyzz_y, g_0_y_xxxyyzz_z, g_0_y_xxyyzz_x, g_0_y_xxyyzz_xx, g_0_y_xxyyzz_xy, g_0_y_xxyyzz_xz, g_0_y_xxyyzz_y, g_0_y_xxyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyzz_x[k] = -g_0_y_xxyyzz_x[k] * ab_x + g_0_y_xxyyzz_xx[k];

                g_0_y_xxxyyzz_y[k] = -g_0_y_xxyyzz_y[k] * ab_x + g_0_y_xxyyzz_xy[k];

                g_0_y_xxxyyzz_z[k] = -g_0_y_xxyyzz_z[k] * ab_x + g_0_y_xxyyzz_xz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyzzz_x = cbuffer.data(kp_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_y = cbuffer.data(kp_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_z = cbuffer.data(kp_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyzzz_x, g_0_y_xxxyzzz_y, g_0_y_xxxyzzz_z, g_0_y_xxyzzz_x, g_0_y_xxyzzz_xx, g_0_y_xxyzzz_xy, g_0_y_xxyzzz_xz, g_0_y_xxyzzz_y, g_0_y_xxyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyzzz_x[k] = -g_0_y_xxyzzz_x[k] * ab_x + g_0_y_xxyzzz_xx[k];

                g_0_y_xxxyzzz_y[k] = -g_0_y_xxyzzz_y[k] * ab_x + g_0_y_xxyzzz_xy[k];

                g_0_y_xxxyzzz_z[k] = -g_0_y_xxyzzz_z[k] * ab_x + g_0_y_xxyzzz_xz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzzzz_x = cbuffer.data(kp_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_y = cbuffer.data(kp_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_z = cbuffer.data(kp_geom_01_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzzzz_x, g_0_y_xxxzzzz_y, g_0_y_xxxzzzz_z, g_0_y_xxzzzz_x, g_0_y_xxzzzz_xx, g_0_y_xxzzzz_xy, g_0_y_xxzzzz_xz, g_0_y_xxzzzz_y, g_0_y_xxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzzzz_x[k] = -g_0_y_xxzzzz_x[k] * ab_x + g_0_y_xxzzzz_xx[k];

                g_0_y_xxxzzzz_y[k] = -g_0_y_xxzzzz_y[k] * ab_x + g_0_y_xxzzzz_xy[k];

                g_0_y_xxxzzzz_z[k] = -g_0_y_xxzzzz_z[k] * ab_x + g_0_y_xxzzzz_xz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyyy_x = cbuffer.data(kp_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_y = cbuffer.data(kp_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_z = cbuffer.data(kp_geom_01_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyyy_x, g_0_y_xxyyyyy_y, g_0_y_xxyyyyy_z, g_0_y_xyyyyy_x, g_0_y_xyyyyy_xx, g_0_y_xyyyyy_xy, g_0_y_xyyyyy_xz, g_0_y_xyyyyy_y, g_0_y_xyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyyy_x[k] = -g_0_y_xyyyyy_x[k] * ab_x + g_0_y_xyyyyy_xx[k];

                g_0_y_xxyyyyy_y[k] = -g_0_y_xyyyyy_y[k] * ab_x + g_0_y_xyyyyy_xy[k];

                g_0_y_xxyyyyy_z[k] = -g_0_y_xyyyyy_z[k] * ab_x + g_0_y_xyyyyy_xz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyyz_x = cbuffer.data(kp_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_y = cbuffer.data(kp_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_z = cbuffer.data(kp_geom_01_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyyz_x, g_0_y_xxyyyyz_y, g_0_y_xxyyyyz_z, g_0_y_xyyyyz_x, g_0_y_xyyyyz_xx, g_0_y_xyyyyz_xy, g_0_y_xyyyyz_xz, g_0_y_xyyyyz_y, g_0_y_xyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyyz_x[k] = -g_0_y_xyyyyz_x[k] * ab_x + g_0_y_xyyyyz_xx[k];

                g_0_y_xxyyyyz_y[k] = -g_0_y_xyyyyz_y[k] * ab_x + g_0_y_xyyyyz_xy[k];

                g_0_y_xxyyyyz_z[k] = -g_0_y_xyyyyz_z[k] * ab_x + g_0_y_xyyyyz_xz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyzz_x = cbuffer.data(kp_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_y = cbuffer.data(kp_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_z = cbuffer.data(kp_geom_01_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyzz_x, g_0_y_xxyyyzz_y, g_0_y_xxyyyzz_z, g_0_y_xyyyzz_x, g_0_y_xyyyzz_xx, g_0_y_xyyyzz_xy, g_0_y_xyyyzz_xz, g_0_y_xyyyzz_y, g_0_y_xyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyzz_x[k] = -g_0_y_xyyyzz_x[k] * ab_x + g_0_y_xyyyzz_xx[k];

                g_0_y_xxyyyzz_y[k] = -g_0_y_xyyyzz_y[k] * ab_x + g_0_y_xyyyzz_xy[k];

                g_0_y_xxyyyzz_z[k] = -g_0_y_xyyyzz_z[k] * ab_x + g_0_y_xyyyzz_xz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyzzz_x = cbuffer.data(kp_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_y = cbuffer.data(kp_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_z = cbuffer.data(kp_geom_01_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyzzz_x, g_0_y_xxyyzzz_y, g_0_y_xxyyzzz_z, g_0_y_xyyzzz_x, g_0_y_xyyzzz_xx, g_0_y_xyyzzz_xy, g_0_y_xyyzzz_xz, g_0_y_xyyzzz_y, g_0_y_xyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyzzz_x[k] = -g_0_y_xyyzzz_x[k] * ab_x + g_0_y_xyyzzz_xx[k];

                g_0_y_xxyyzzz_y[k] = -g_0_y_xyyzzz_y[k] * ab_x + g_0_y_xyyzzz_xy[k];

                g_0_y_xxyyzzz_z[k] = -g_0_y_xyyzzz_z[k] * ab_x + g_0_y_xyyzzz_xz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzzzz_x = cbuffer.data(kp_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_y = cbuffer.data(kp_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_z = cbuffer.data(kp_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzzzz_x, g_0_y_xxyzzzz_y, g_0_y_xxyzzzz_z, g_0_y_xyzzzz_x, g_0_y_xyzzzz_xx, g_0_y_xyzzzz_xy, g_0_y_xyzzzz_xz, g_0_y_xyzzzz_y, g_0_y_xyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzzzz_x[k] = -g_0_y_xyzzzz_x[k] * ab_x + g_0_y_xyzzzz_xx[k];

                g_0_y_xxyzzzz_y[k] = -g_0_y_xyzzzz_y[k] * ab_x + g_0_y_xyzzzz_xy[k];

                g_0_y_xxyzzzz_z[k] = -g_0_y_xyzzzz_z[k] * ab_x + g_0_y_xyzzzz_xz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzzzz_x = cbuffer.data(kp_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_y = cbuffer.data(kp_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_z = cbuffer.data(kp_geom_01_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzzzz_x, g_0_y_xxzzzzz_y, g_0_y_xxzzzzz_z, g_0_y_xzzzzz_x, g_0_y_xzzzzz_xx, g_0_y_xzzzzz_xy, g_0_y_xzzzzz_xz, g_0_y_xzzzzz_y, g_0_y_xzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzzzz_x[k] = -g_0_y_xzzzzz_x[k] * ab_x + g_0_y_xzzzzz_xx[k];

                g_0_y_xxzzzzz_y[k] = -g_0_y_xzzzzz_y[k] * ab_x + g_0_y_xzzzzz_xy[k];

                g_0_y_xxzzzzz_z[k] = -g_0_y_xzzzzz_z[k] * ab_x + g_0_y_xzzzzz_xz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyyy_x = cbuffer.data(kp_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_y = cbuffer.data(kp_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_z = cbuffer.data(kp_geom_01_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyyy_x, g_0_y_xyyyyyy_y, g_0_y_xyyyyyy_z, g_0_y_yyyyyy_x, g_0_y_yyyyyy_xx, g_0_y_yyyyyy_xy, g_0_y_yyyyyy_xz, g_0_y_yyyyyy_y, g_0_y_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyyy_x[k] = -g_0_y_yyyyyy_x[k] * ab_x + g_0_y_yyyyyy_xx[k];

                g_0_y_xyyyyyy_y[k] = -g_0_y_yyyyyy_y[k] * ab_x + g_0_y_yyyyyy_xy[k];

                g_0_y_xyyyyyy_z[k] = -g_0_y_yyyyyy_z[k] * ab_x + g_0_y_yyyyyy_xz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyyz_x = cbuffer.data(kp_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_y = cbuffer.data(kp_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_z = cbuffer.data(kp_geom_01_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyyz_x, g_0_y_xyyyyyz_y, g_0_y_xyyyyyz_z, g_0_y_yyyyyz_x, g_0_y_yyyyyz_xx, g_0_y_yyyyyz_xy, g_0_y_yyyyyz_xz, g_0_y_yyyyyz_y, g_0_y_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyyz_x[k] = -g_0_y_yyyyyz_x[k] * ab_x + g_0_y_yyyyyz_xx[k];

                g_0_y_xyyyyyz_y[k] = -g_0_y_yyyyyz_y[k] * ab_x + g_0_y_yyyyyz_xy[k];

                g_0_y_xyyyyyz_z[k] = -g_0_y_yyyyyz_z[k] * ab_x + g_0_y_yyyyyz_xz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyzz_x = cbuffer.data(kp_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_y = cbuffer.data(kp_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_z = cbuffer.data(kp_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyzz_x, g_0_y_xyyyyzz_y, g_0_y_xyyyyzz_z, g_0_y_yyyyzz_x, g_0_y_yyyyzz_xx, g_0_y_yyyyzz_xy, g_0_y_yyyyzz_xz, g_0_y_yyyyzz_y, g_0_y_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyzz_x[k] = -g_0_y_yyyyzz_x[k] * ab_x + g_0_y_yyyyzz_xx[k];

                g_0_y_xyyyyzz_y[k] = -g_0_y_yyyyzz_y[k] * ab_x + g_0_y_yyyyzz_xy[k];

                g_0_y_xyyyyzz_z[k] = -g_0_y_yyyyzz_z[k] * ab_x + g_0_y_yyyyzz_xz[k];
            }

            /// Set up 180-183 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyzzz_x = cbuffer.data(kp_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_y = cbuffer.data(kp_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_z = cbuffer.data(kp_geom_01_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyzzz_x, g_0_y_xyyyzzz_y, g_0_y_xyyyzzz_z, g_0_y_yyyzzz_x, g_0_y_yyyzzz_xx, g_0_y_yyyzzz_xy, g_0_y_yyyzzz_xz, g_0_y_yyyzzz_y, g_0_y_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyzzz_x[k] = -g_0_y_yyyzzz_x[k] * ab_x + g_0_y_yyyzzz_xx[k];

                g_0_y_xyyyzzz_y[k] = -g_0_y_yyyzzz_y[k] * ab_x + g_0_y_yyyzzz_xy[k];

                g_0_y_xyyyzzz_z[k] = -g_0_y_yyyzzz_z[k] * ab_x + g_0_y_yyyzzz_xz[k];
            }

            /// Set up 183-186 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzzzz_x = cbuffer.data(kp_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_y = cbuffer.data(kp_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_z = cbuffer.data(kp_geom_01_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzzzz_x, g_0_y_xyyzzzz_y, g_0_y_xyyzzzz_z, g_0_y_yyzzzz_x, g_0_y_yyzzzz_xx, g_0_y_yyzzzz_xy, g_0_y_yyzzzz_xz, g_0_y_yyzzzz_y, g_0_y_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzzzz_x[k] = -g_0_y_yyzzzz_x[k] * ab_x + g_0_y_yyzzzz_xx[k];

                g_0_y_xyyzzzz_y[k] = -g_0_y_yyzzzz_y[k] * ab_x + g_0_y_yyzzzz_xy[k];

                g_0_y_xyyzzzz_z[k] = -g_0_y_yyzzzz_z[k] * ab_x + g_0_y_yyzzzz_xz[k];
            }

            /// Set up 186-189 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzzzz_x = cbuffer.data(kp_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_y = cbuffer.data(kp_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_z = cbuffer.data(kp_geom_01_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzzzz_x, g_0_y_xyzzzzz_y, g_0_y_xyzzzzz_z, g_0_y_yzzzzz_x, g_0_y_yzzzzz_xx, g_0_y_yzzzzz_xy, g_0_y_yzzzzz_xz, g_0_y_yzzzzz_y, g_0_y_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzzzz_x[k] = -g_0_y_yzzzzz_x[k] * ab_x + g_0_y_yzzzzz_xx[k];

                g_0_y_xyzzzzz_y[k] = -g_0_y_yzzzzz_y[k] * ab_x + g_0_y_yzzzzz_xy[k];

                g_0_y_xyzzzzz_z[k] = -g_0_y_yzzzzz_z[k] * ab_x + g_0_y_yzzzzz_xz[k];
            }

            /// Set up 189-192 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzzzz_x = cbuffer.data(kp_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_y = cbuffer.data(kp_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_z = cbuffer.data(kp_geom_01_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzzzz_x, g_0_y_xzzzzzz_y, g_0_y_xzzzzzz_z, g_0_y_zzzzzz_x, g_0_y_zzzzzz_xx, g_0_y_zzzzzz_xy, g_0_y_zzzzzz_xz, g_0_y_zzzzzz_y, g_0_y_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzzzz_x[k] = -g_0_y_zzzzzz_x[k] * ab_x + g_0_y_zzzzzz_xx[k];

                g_0_y_xzzzzzz_y[k] = -g_0_y_zzzzzz_y[k] * ab_x + g_0_y_zzzzzz_xy[k];

                g_0_y_xzzzzzz_z[k] = -g_0_y_zzzzzz_z[k] * ab_x + g_0_y_zzzzzz_xz[k];
            }

            /// Set up 192-195 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyyy_x = cbuffer.data(kp_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_y = cbuffer.data(kp_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_z = cbuffer.data(kp_geom_01_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyy_x, g_0_y_yyyyyy_xy, g_0_y_yyyyyy_y, g_0_y_yyyyyy_yy, g_0_y_yyyyyy_yz, g_0_y_yyyyyy_z, g_0_y_yyyyyyy_x, g_0_y_yyyyyyy_y, g_0_y_yyyyyyy_z, g_yyyyyy_x, g_yyyyyy_y, g_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyyy_x[k] = g_yyyyyy_x[k] - g_0_y_yyyyyy_x[k] * ab_y + g_0_y_yyyyyy_xy[k];

                g_0_y_yyyyyyy_y[k] = g_yyyyyy_y[k] - g_0_y_yyyyyy_y[k] * ab_y + g_0_y_yyyyyy_yy[k];

                g_0_y_yyyyyyy_z[k] = g_yyyyyy_z[k] - g_0_y_yyyyyy_z[k] * ab_y + g_0_y_yyyyyy_yz[k];
            }

            /// Set up 195-198 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyyz_x = cbuffer.data(kp_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_y = cbuffer.data(kp_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_z = cbuffer.data(kp_geom_01_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyy_x, g_0_y_yyyyyy_xz, g_0_y_yyyyyy_y, g_0_y_yyyyyy_yz, g_0_y_yyyyyy_z, g_0_y_yyyyyy_zz, g_0_y_yyyyyyz_x, g_0_y_yyyyyyz_y, g_0_y_yyyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyyz_x[k] = -g_0_y_yyyyyy_x[k] * ab_z + g_0_y_yyyyyy_xz[k];

                g_0_y_yyyyyyz_y[k] = -g_0_y_yyyyyy_y[k] * ab_z + g_0_y_yyyyyy_yz[k];

                g_0_y_yyyyyyz_z[k] = -g_0_y_yyyyyy_z[k] * ab_z + g_0_y_yyyyyy_zz[k];
            }

            /// Set up 198-201 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyzz_x = cbuffer.data(kp_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_y = cbuffer.data(kp_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_z = cbuffer.data(kp_geom_01_off + 200 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyz_x, g_0_y_yyyyyz_xz, g_0_y_yyyyyz_y, g_0_y_yyyyyz_yz, g_0_y_yyyyyz_z, g_0_y_yyyyyz_zz, g_0_y_yyyyyzz_x, g_0_y_yyyyyzz_y, g_0_y_yyyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyzz_x[k] = -g_0_y_yyyyyz_x[k] * ab_z + g_0_y_yyyyyz_xz[k];

                g_0_y_yyyyyzz_y[k] = -g_0_y_yyyyyz_y[k] * ab_z + g_0_y_yyyyyz_yz[k];

                g_0_y_yyyyyzz_z[k] = -g_0_y_yyyyyz_z[k] * ab_z + g_0_y_yyyyyz_zz[k];
            }

            /// Set up 201-204 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyzzz_x = cbuffer.data(kp_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_y = cbuffer.data(kp_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_z = cbuffer.data(kp_geom_01_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyzz_x, g_0_y_yyyyzz_xz, g_0_y_yyyyzz_y, g_0_y_yyyyzz_yz, g_0_y_yyyyzz_z, g_0_y_yyyyzz_zz, g_0_y_yyyyzzz_x, g_0_y_yyyyzzz_y, g_0_y_yyyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyzzz_x[k] = -g_0_y_yyyyzz_x[k] * ab_z + g_0_y_yyyyzz_xz[k];

                g_0_y_yyyyzzz_y[k] = -g_0_y_yyyyzz_y[k] * ab_z + g_0_y_yyyyzz_yz[k];

                g_0_y_yyyyzzz_z[k] = -g_0_y_yyyyzz_z[k] * ab_z + g_0_y_yyyyzz_zz[k];
            }

            /// Set up 204-207 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzzzz_x = cbuffer.data(kp_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_y = cbuffer.data(kp_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_z = cbuffer.data(kp_geom_01_off + 206 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyzzz_x, g_0_y_yyyzzz_xz, g_0_y_yyyzzz_y, g_0_y_yyyzzz_yz, g_0_y_yyyzzz_z, g_0_y_yyyzzz_zz, g_0_y_yyyzzzz_x, g_0_y_yyyzzzz_y, g_0_y_yyyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzzzz_x[k] = -g_0_y_yyyzzz_x[k] * ab_z + g_0_y_yyyzzz_xz[k];

                g_0_y_yyyzzzz_y[k] = -g_0_y_yyyzzz_y[k] * ab_z + g_0_y_yyyzzz_yz[k];

                g_0_y_yyyzzzz_z[k] = -g_0_y_yyyzzz_z[k] * ab_z + g_0_y_yyyzzz_zz[k];
            }

            /// Set up 207-210 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzzzz_x = cbuffer.data(kp_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_y = cbuffer.data(kp_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_z = cbuffer.data(kp_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzzzz_x, g_0_y_yyzzzz_xz, g_0_y_yyzzzz_y, g_0_y_yyzzzz_yz, g_0_y_yyzzzz_z, g_0_y_yyzzzz_zz, g_0_y_yyzzzzz_x, g_0_y_yyzzzzz_y, g_0_y_yyzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzzzz_x[k] = -g_0_y_yyzzzz_x[k] * ab_z + g_0_y_yyzzzz_xz[k];

                g_0_y_yyzzzzz_y[k] = -g_0_y_yyzzzz_y[k] * ab_z + g_0_y_yyzzzz_yz[k];

                g_0_y_yyzzzzz_z[k] = -g_0_y_yyzzzz_z[k] * ab_z + g_0_y_yyzzzz_zz[k];
            }

            /// Set up 210-213 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzzzz_x = cbuffer.data(kp_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_y = cbuffer.data(kp_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_z = cbuffer.data(kp_geom_01_off + 212 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzzzz_x, g_0_y_yzzzzz_xz, g_0_y_yzzzzz_y, g_0_y_yzzzzz_yz, g_0_y_yzzzzz_z, g_0_y_yzzzzz_zz, g_0_y_yzzzzzz_x, g_0_y_yzzzzzz_y, g_0_y_yzzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzzzz_x[k] = -g_0_y_yzzzzz_x[k] * ab_z + g_0_y_yzzzzz_xz[k];

                g_0_y_yzzzzzz_y[k] = -g_0_y_yzzzzz_y[k] * ab_z + g_0_y_yzzzzz_yz[k];

                g_0_y_yzzzzzz_z[k] = -g_0_y_yzzzzz_z[k] * ab_z + g_0_y_yzzzzz_zz[k];
            }

            /// Set up 213-216 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzzzz_x = cbuffer.data(kp_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_y = cbuffer.data(kp_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_z = cbuffer.data(kp_geom_01_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzzzz_x, g_0_y_zzzzzz_xz, g_0_y_zzzzzz_y, g_0_y_zzzzzz_yz, g_0_y_zzzzzz_z, g_0_y_zzzzzz_zz, g_0_y_zzzzzzz_x, g_0_y_zzzzzzz_y, g_0_y_zzzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzzzz_x[k] = -g_0_y_zzzzzz_x[k] * ab_z + g_0_y_zzzzzz_xz[k];

                g_0_y_zzzzzzz_y[k] = -g_0_y_zzzzzz_y[k] * ab_z + g_0_y_zzzzzz_yz[k];

                g_0_y_zzzzzzz_z[k] = -g_0_y_zzzzzz_z[k] * ab_z + g_0_y_zzzzzz_zz[k];
            }

            /// Set up 216-219 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxx_x = cbuffer.data(kp_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_y = cbuffer.data(kp_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_z = cbuffer.data(kp_geom_01_off + 218 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxx_x, g_0_z_xxxxxx_xx, g_0_z_xxxxxx_xy, g_0_z_xxxxxx_xz, g_0_z_xxxxxx_y, g_0_z_xxxxxx_z, g_0_z_xxxxxxx_x, g_0_z_xxxxxxx_y, g_0_z_xxxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxx_x[k] = -g_0_z_xxxxxx_x[k] * ab_x + g_0_z_xxxxxx_xx[k];

                g_0_z_xxxxxxx_y[k] = -g_0_z_xxxxxx_y[k] * ab_x + g_0_z_xxxxxx_xy[k];

                g_0_z_xxxxxxx_z[k] = -g_0_z_xxxxxx_z[k] * ab_x + g_0_z_xxxxxx_xz[k];
            }

            /// Set up 219-222 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxy_x = cbuffer.data(kp_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_y = cbuffer.data(kp_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_z = cbuffer.data(kp_geom_01_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxxy_x, g_0_z_xxxxxxy_y, g_0_z_xxxxxxy_z, g_0_z_xxxxxy_x, g_0_z_xxxxxy_xx, g_0_z_xxxxxy_xy, g_0_z_xxxxxy_xz, g_0_z_xxxxxy_y, g_0_z_xxxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxy_x[k] = -g_0_z_xxxxxy_x[k] * ab_x + g_0_z_xxxxxy_xx[k];

                g_0_z_xxxxxxy_y[k] = -g_0_z_xxxxxy_y[k] * ab_x + g_0_z_xxxxxy_xy[k];

                g_0_z_xxxxxxy_z[k] = -g_0_z_xxxxxy_z[k] * ab_x + g_0_z_xxxxxy_xz[k];
            }

            /// Set up 222-225 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxz_x = cbuffer.data(kp_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_y = cbuffer.data(kp_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_z = cbuffer.data(kp_geom_01_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxxz_x, g_0_z_xxxxxxz_y, g_0_z_xxxxxxz_z, g_0_z_xxxxxz_x, g_0_z_xxxxxz_xx, g_0_z_xxxxxz_xy, g_0_z_xxxxxz_xz, g_0_z_xxxxxz_y, g_0_z_xxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxz_x[k] = -g_0_z_xxxxxz_x[k] * ab_x + g_0_z_xxxxxz_xx[k];

                g_0_z_xxxxxxz_y[k] = -g_0_z_xxxxxz_y[k] * ab_x + g_0_z_xxxxxz_xy[k];

                g_0_z_xxxxxxz_z[k] = -g_0_z_xxxxxz_z[k] * ab_x + g_0_z_xxxxxz_xz[k];
            }

            /// Set up 225-228 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxyy_x = cbuffer.data(kp_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_y = cbuffer.data(kp_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_z = cbuffer.data(kp_geom_01_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxyy_x, g_0_z_xxxxxyy_y, g_0_z_xxxxxyy_z, g_0_z_xxxxyy_x, g_0_z_xxxxyy_xx, g_0_z_xxxxyy_xy, g_0_z_xxxxyy_xz, g_0_z_xxxxyy_y, g_0_z_xxxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxyy_x[k] = -g_0_z_xxxxyy_x[k] * ab_x + g_0_z_xxxxyy_xx[k];

                g_0_z_xxxxxyy_y[k] = -g_0_z_xxxxyy_y[k] * ab_x + g_0_z_xxxxyy_xy[k];

                g_0_z_xxxxxyy_z[k] = -g_0_z_xxxxyy_z[k] * ab_x + g_0_z_xxxxyy_xz[k];
            }

            /// Set up 228-231 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxyz_x = cbuffer.data(kp_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_y = cbuffer.data(kp_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_z = cbuffer.data(kp_geom_01_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxyz_x, g_0_z_xxxxxyz_y, g_0_z_xxxxxyz_z, g_0_z_xxxxyz_x, g_0_z_xxxxyz_xx, g_0_z_xxxxyz_xy, g_0_z_xxxxyz_xz, g_0_z_xxxxyz_y, g_0_z_xxxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxyz_x[k] = -g_0_z_xxxxyz_x[k] * ab_x + g_0_z_xxxxyz_xx[k];

                g_0_z_xxxxxyz_y[k] = -g_0_z_xxxxyz_y[k] * ab_x + g_0_z_xxxxyz_xy[k];

                g_0_z_xxxxxyz_z[k] = -g_0_z_xxxxyz_z[k] * ab_x + g_0_z_xxxxyz_xz[k];
            }

            /// Set up 231-234 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxzz_x = cbuffer.data(kp_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_y = cbuffer.data(kp_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_z = cbuffer.data(kp_geom_01_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxzz_x, g_0_z_xxxxxzz_y, g_0_z_xxxxxzz_z, g_0_z_xxxxzz_x, g_0_z_xxxxzz_xx, g_0_z_xxxxzz_xy, g_0_z_xxxxzz_xz, g_0_z_xxxxzz_y, g_0_z_xxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxzz_x[k] = -g_0_z_xxxxzz_x[k] * ab_x + g_0_z_xxxxzz_xx[k];

                g_0_z_xxxxxzz_y[k] = -g_0_z_xxxxzz_y[k] * ab_x + g_0_z_xxxxzz_xy[k];

                g_0_z_xxxxxzz_z[k] = -g_0_z_xxxxzz_z[k] * ab_x + g_0_z_xxxxzz_xz[k];
            }

            /// Set up 234-237 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyyy_x = cbuffer.data(kp_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_y = cbuffer.data(kp_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_z = cbuffer.data(kp_geom_01_off + 236 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyyy_x, g_0_z_xxxxyyy_y, g_0_z_xxxxyyy_z, g_0_z_xxxyyy_x, g_0_z_xxxyyy_xx, g_0_z_xxxyyy_xy, g_0_z_xxxyyy_xz, g_0_z_xxxyyy_y, g_0_z_xxxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyyy_x[k] = -g_0_z_xxxyyy_x[k] * ab_x + g_0_z_xxxyyy_xx[k];

                g_0_z_xxxxyyy_y[k] = -g_0_z_xxxyyy_y[k] * ab_x + g_0_z_xxxyyy_xy[k];

                g_0_z_xxxxyyy_z[k] = -g_0_z_xxxyyy_z[k] * ab_x + g_0_z_xxxyyy_xz[k];
            }

            /// Set up 237-240 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyyz_x = cbuffer.data(kp_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_y = cbuffer.data(kp_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_z = cbuffer.data(kp_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyyz_x, g_0_z_xxxxyyz_y, g_0_z_xxxxyyz_z, g_0_z_xxxyyz_x, g_0_z_xxxyyz_xx, g_0_z_xxxyyz_xy, g_0_z_xxxyyz_xz, g_0_z_xxxyyz_y, g_0_z_xxxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyyz_x[k] = -g_0_z_xxxyyz_x[k] * ab_x + g_0_z_xxxyyz_xx[k];

                g_0_z_xxxxyyz_y[k] = -g_0_z_xxxyyz_y[k] * ab_x + g_0_z_xxxyyz_xy[k];

                g_0_z_xxxxyyz_z[k] = -g_0_z_xxxyyz_z[k] * ab_x + g_0_z_xxxyyz_xz[k];
            }

            /// Set up 240-243 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyzz_x = cbuffer.data(kp_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_y = cbuffer.data(kp_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_z = cbuffer.data(kp_geom_01_off + 242 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyzz_x, g_0_z_xxxxyzz_y, g_0_z_xxxxyzz_z, g_0_z_xxxyzz_x, g_0_z_xxxyzz_xx, g_0_z_xxxyzz_xy, g_0_z_xxxyzz_xz, g_0_z_xxxyzz_y, g_0_z_xxxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyzz_x[k] = -g_0_z_xxxyzz_x[k] * ab_x + g_0_z_xxxyzz_xx[k];

                g_0_z_xxxxyzz_y[k] = -g_0_z_xxxyzz_y[k] * ab_x + g_0_z_xxxyzz_xy[k];

                g_0_z_xxxxyzz_z[k] = -g_0_z_xxxyzz_z[k] * ab_x + g_0_z_xxxyzz_xz[k];
            }

            /// Set up 243-246 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxzzz_x = cbuffer.data(kp_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_y = cbuffer.data(kp_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_z = cbuffer.data(kp_geom_01_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxzzz_x, g_0_z_xxxxzzz_y, g_0_z_xxxxzzz_z, g_0_z_xxxzzz_x, g_0_z_xxxzzz_xx, g_0_z_xxxzzz_xy, g_0_z_xxxzzz_xz, g_0_z_xxxzzz_y, g_0_z_xxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxzzz_x[k] = -g_0_z_xxxzzz_x[k] * ab_x + g_0_z_xxxzzz_xx[k];

                g_0_z_xxxxzzz_y[k] = -g_0_z_xxxzzz_y[k] * ab_x + g_0_z_xxxzzz_xy[k];

                g_0_z_xxxxzzz_z[k] = -g_0_z_xxxzzz_z[k] * ab_x + g_0_z_xxxzzz_xz[k];
            }

            /// Set up 246-249 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyyy_x = cbuffer.data(kp_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_y = cbuffer.data(kp_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_z = cbuffer.data(kp_geom_01_off + 248 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyyy_x, g_0_z_xxxyyyy_y, g_0_z_xxxyyyy_z, g_0_z_xxyyyy_x, g_0_z_xxyyyy_xx, g_0_z_xxyyyy_xy, g_0_z_xxyyyy_xz, g_0_z_xxyyyy_y, g_0_z_xxyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyyy_x[k] = -g_0_z_xxyyyy_x[k] * ab_x + g_0_z_xxyyyy_xx[k];

                g_0_z_xxxyyyy_y[k] = -g_0_z_xxyyyy_y[k] * ab_x + g_0_z_xxyyyy_xy[k];

                g_0_z_xxxyyyy_z[k] = -g_0_z_xxyyyy_z[k] * ab_x + g_0_z_xxyyyy_xz[k];
            }

            /// Set up 249-252 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyyz_x = cbuffer.data(kp_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_y = cbuffer.data(kp_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_z = cbuffer.data(kp_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyyz_x, g_0_z_xxxyyyz_y, g_0_z_xxxyyyz_z, g_0_z_xxyyyz_x, g_0_z_xxyyyz_xx, g_0_z_xxyyyz_xy, g_0_z_xxyyyz_xz, g_0_z_xxyyyz_y, g_0_z_xxyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyyz_x[k] = -g_0_z_xxyyyz_x[k] * ab_x + g_0_z_xxyyyz_xx[k];

                g_0_z_xxxyyyz_y[k] = -g_0_z_xxyyyz_y[k] * ab_x + g_0_z_xxyyyz_xy[k];

                g_0_z_xxxyyyz_z[k] = -g_0_z_xxyyyz_z[k] * ab_x + g_0_z_xxyyyz_xz[k];
            }

            /// Set up 252-255 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyzz_x = cbuffer.data(kp_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_y = cbuffer.data(kp_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_z = cbuffer.data(kp_geom_01_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyzz_x, g_0_z_xxxyyzz_y, g_0_z_xxxyyzz_z, g_0_z_xxyyzz_x, g_0_z_xxyyzz_xx, g_0_z_xxyyzz_xy, g_0_z_xxyyzz_xz, g_0_z_xxyyzz_y, g_0_z_xxyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyzz_x[k] = -g_0_z_xxyyzz_x[k] * ab_x + g_0_z_xxyyzz_xx[k];

                g_0_z_xxxyyzz_y[k] = -g_0_z_xxyyzz_y[k] * ab_x + g_0_z_xxyyzz_xy[k];

                g_0_z_xxxyyzz_z[k] = -g_0_z_xxyyzz_z[k] * ab_x + g_0_z_xxyyzz_xz[k];
            }

            /// Set up 255-258 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyzzz_x = cbuffer.data(kp_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_y = cbuffer.data(kp_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_z = cbuffer.data(kp_geom_01_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyzzz_x, g_0_z_xxxyzzz_y, g_0_z_xxxyzzz_z, g_0_z_xxyzzz_x, g_0_z_xxyzzz_xx, g_0_z_xxyzzz_xy, g_0_z_xxyzzz_xz, g_0_z_xxyzzz_y, g_0_z_xxyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyzzz_x[k] = -g_0_z_xxyzzz_x[k] * ab_x + g_0_z_xxyzzz_xx[k];

                g_0_z_xxxyzzz_y[k] = -g_0_z_xxyzzz_y[k] * ab_x + g_0_z_xxyzzz_xy[k];

                g_0_z_xxxyzzz_z[k] = -g_0_z_xxyzzz_z[k] * ab_x + g_0_z_xxyzzz_xz[k];
            }

            /// Set up 258-261 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzzzz_x = cbuffer.data(kp_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_y = cbuffer.data(kp_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_z = cbuffer.data(kp_geom_01_off + 260 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzzzz_x, g_0_z_xxxzzzz_y, g_0_z_xxxzzzz_z, g_0_z_xxzzzz_x, g_0_z_xxzzzz_xx, g_0_z_xxzzzz_xy, g_0_z_xxzzzz_xz, g_0_z_xxzzzz_y, g_0_z_xxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzzzz_x[k] = -g_0_z_xxzzzz_x[k] * ab_x + g_0_z_xxzzzz_xx[k];

                g_0_z_xxxzzzz_y[k] = -g_0_z_xxzzzz_y[k] * ab_x + g_0_z_xxzzzz_xy[k];

                g_0_z_xxxzzzz_z[k] = -g_0_z_xxzzzz_z[k] * ab_x + g_0_z_xxzzzz_xz[k];
            }

            /// Set up 261-264 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyyy_x = cbuffer.data(kp_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_y = cbuffer.data(kp_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_z = cbuffer.data(kp_geom_01_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyyy_x, g_0_z_xxyyyyy_y, g_0_z_xxyyyyy_z, g_0_z_xyyyyy_x, g_0_z_xyyyyy_xx, g_0_z_xyyyyy_xy, g_0_z_xyyyyy_xz, g_0_z_xyyyyy_y, g_0_z_xyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyyy_x[k] = -g_0_z_xyyyyy_x[k] * ab_x + g_0_z_xyyyyy_xx[k];

                g_0_z_xxyyyyy_y[k] = -g_0_z_xyyyyy_y[k] * ab_x + g_0_z_xyyyyy_xy[k];

                g_0_z_xxyyyyy_z[k] = -g_0_z_xyyyyy_z[k] * ab_x + g_0_z_xyyyyy_xz[k];
            }

            /// Set up 264-267 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyyz_x = cbuffer.data(kp_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_y = cbuffer.data(kp_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_z = cbuffer.data(kp_geom_01_off + 266 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyyz_x, g_0_z_xxyyyyz_y, g_0_z_xxyyyyz_z, g_0_z_xyyyyz_x, g_0_z_xyyyyz_xx, g_0_z_xyyyyz_xy, g_0_z_xyyyyz_xz, g_0_z_xyyyyz_y, g_0_z_xyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyyz_x[k] = -g_0_z_xyyyyz_x[k] * ab_x + g_0_z_xyyyyz_xx[k];

                g_0_z_xxyyyyz_y[k] = -g_0_z_xyyyyz_y[k] * ab_x + g_0_z_xyyyyz_xy[k];

                g_0_z_xxyyyyz_z[k] = -g_0_z_xyyyyz_z[k] * ab_x + g_0_z_xyyyyz_xz[k];
            }

            /// Set up 267-270 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyzz_x = cbuffer.data(kp_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_y = cbuffer.data(kp_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_z = cbuffer.data(kp_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyzz_x, g_0_z_xxyyyzz_y, g_0_z_xxyyyzz_z, g_0_z_xyyyzz_x, g_0_z_xyyyzz_xx, g_0_z_xyyyzz_xy, g_0_z_xyyyzz_xz, g_0_z_xyyyzz_y, g_0_z_xyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyzz_x[k] = -g_0_z_xyyyzz_x[k] * ab_x + g_0_z_xyyyzz_xx[k];

                g_0_z_xxyyyzz_y[k] = -g_0_z_xyyyzz_y[k] * ab_x + g_0_z_xyyyzz_xy[k];

                g_0_z_xxyyyzz_z[k] = -g_0_z_xyyyzz_z[k] * ab_x + g_0_z_xyyyzz_xz[k];
            }

            /// Set up 270-273 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyzzz_x = cbuffer.data(kp_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_y = cbuffer.data(kp_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_z = cbuffer.data(kp_geom_01_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyzzz_x, g_0_z_xxyyzzz_y, g_0_z_xxyyzzz_z, g_0_z_xyyzzz_x, g_0_z_xyyzzz_xx, g_0_z_xyyzzz_xy, g_0_z_xyyzzz_xz, g_0_z_xyyzzz_y, g_0_z_xyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyzzz_x[k] = -g_0_z_xyyzzz_x[k] * ab_x + g_0_z_xyyzzz_xx[k];

                g_0_z_xxyyzzz_y[k] = -g_0_z_xyyzzz_y[k] * ab_x + g_0_z_xyyzzz_xy[k];

                g_0_z_xxyyzzz_z[k] = -g_0_z_xyyzzz_z[k] * ab_x + g_0_z_xyyzzz_xz[k];
            }

            /// Set up 273-276 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzzzz_x = cbuffer.data(kp_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_y = cbuffer.data(kp_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_z = cbuffer.data(kp_geom_01_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzzzz_x, g_0_z_xxyzzzz_y, g_0_z_xxyzzzz_z, g_0_z_xyzzzz_x, g_0_z_xyzzzz_xx, g_0_z_xyzzzz_xy, g_0_z_xyzzzz_xz, g_0_z_xyzzzz_y, g_0_z_xyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzzzz_x[k] = -g_0_z_xyzzzz_x[k] * ab_x + g_0_z_xyzzzz_xx[k];

                g_0_z_xxyzzzz_y[k] = -g_0_z_xyzzzz_y[k] * ab_x + g_0_z_xyzzzz_xy[k];

                g_0_z_xxyzzzz_z[k] = -g_0_z_xyzzzz_z[k] * ab_x + g_0_z_xyzzzz_xz[k];
            }

            /// Set up 276-279 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzzzz_x = cbuffer.data(kp_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_y = cbuffer.data(kp_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_z = cbuffer.data(kp_geom_01_off + 278 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzzzz_x, g_0_z_xxzzzzz_y, g_0_z_xxzzzzz_z, g_0_z_xzzzzz_x, g_0_z_xzzzzz_xx, g_0_z_xzzzzz_xy, g_0_z_xzzzzz_xz, g_0_z_xzzzzz_y, g_0_z_xzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzzzz_x[k] = -g_0_z_xzzzzz_x[k] * ab_x + g_0_z_xzzzzz_xx[k];

                g_0_z_xxzzzzz_y[k] = -g_0_z_xzzzzz_y[k] * ab_x + g_0_z_xzzzzz_xy[k];

                g_0_z_xxzzzzz_z[k] = -g_0_z_xzzzzz_z[k] * ab_x + g_0_z_xzzzzz_xz[k];
            }

            /// Set up 279-282 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyyy_x = cbuffer.data(kp_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_y = cbuffer.data(kp_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_z = cbuffer.data(kp_geom_01_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyyy_x, g_0_z_xyyyyyy_y, g_0_z_xyyyyyy_z, g_0_z_yyyyyy_x, g_0_z_yyyyyy_xx, g_0_z_yyyyyy_xy, g_0_z_yyyyyy_xz, g_0_z_yyyyyy_y, g_0_z_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyyy_x[k] = -g_0_z_yyyyyy_x[k] * ab_x + g_0_z_yyyyyy_xx[k];

                g_0_z_xyyyyyy_y[k] = -g_0_z_yyyyyy_y[k] * ab_x + g_0_z_yyyyyy_xy[k];

                g_0_z_xyyyyyy_z[k] = -g_0_z_yyyyyy_z[k] * ab_x + g_0_z_yyyyyy_xz[k];
            }

            /// Set up 282-285 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyyz_x = cbuffer.data(kp_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_y = cbuffer.data(kp_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_z = cbuffer.data(kp_geom_01_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyyz_x, g_0_z_xyyyyyz_y, g_0_z_xyyyyyz_z, g_0_z_yyyyyz_x, g_0_z_yyyyyz_xx, g_0_z_yyyyyz_xy, g_0_z_yyyyyz_xz, g_0_z_yyyyyz_y, g_0_z_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyyz_x[k] = -g_0_z_yyyyyz_x[k] * ab_x + g_0_z_yyyyyz_xx[k];

                g_0_z_xyyyyyz_y[k] = -g_0_z_yyyyyz_y[k] * ab_x + g_0_z_yyyyyz_xy[k];

                g_0_z_xyyyyyz_z[k] = -g_0_z_yyyyyz_z[k] * ab_x + g_0_z_yyyyyz_xz[k];
            }

            /// Set up 285-288 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyzz_x = cbuffer.data(kp_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_y = cbuffer.data(kp_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_z = cbuffer.data(kp_geom_01_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyzz_x, g_0_z_xyyyyzz_y, g_0_z_xyyyyzz_z, g_0_z_yyyyzz_x, g_0_z_yyyyzz_xx, g_0_z_yyyyzz_xy, g_0_z_yyyyzz_xz, g_0_z_yyyyzz_y, g_0_z_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyzz_x[k] = -g_0_z_yyyyzz_x[k] * ab_x + g_0_z_yyyyzz_xx[k];

                g_0_z_xyyyyzz_y[k] = -g_0_z_yyyyzz_y[k] * ab_x + g_0_z_yyyyzz_xy[k];

                g_0_z_xyyyyzz_z[k] = -g_0_z_yyyyzz_z[k] * ab_x + g_0_z_yyyyzz_xz[k];
            }

            /// Set up 288-291 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyzzz_x = cbuffer.data(kp_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_y = cbuffer.data(kp_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_z = cbuffer.data(kp_geom_01_off + 290 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyzzz_x, g_0_z_xyyyzzz_y, g_0_z_xyyyzzz_z, g_0_z_yyyzzz_x, g_0_z_yyyzzz_xx, g_0_z_yyyzzz_xy, g_0_z_yyyzzz_xz, g_0_z_yyyzzz_y, g_0_z_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyzzz_x[k] = -g_0_z_yyyzzz_x[k] * ab_x + g_0_z_yyyzzz_xx[k];

                g_0_z_xyyyzzz_y[k] = -g_0_z_yyyzzz_y[k] * ab_x + g_0_z_yyyzzz_xy[k];

                g_0_z_xyyyzzz_z[k] = -g_0_z_yyyzzz_z[k] * ab_x + g_0_z_yyyzzz_xz[k];
            }

            /// Set up 291-294 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzzzz_x = cbuffer.data(kp_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_y = cbuffer.data(kp_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_z = cbuffer.data(kp_geom_01_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzzzz_x, g_0_z_xyyzzzz_y, g_0_z_xyyzzzz_z, g_0_z_yyzzzz_x, g_0_z_yyzzzz_xx, g_0_z_yyzzzz_xy, g_0_z_yyzzzz_xz, g_0_z_yyzzzz_y, g_0_z_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzzzz_x[k] = -g_0_z_yyzzzz_x[k] * ab_x + g_0_z_yyzzzz_xx[k];

                g_0_z_xyyzzzz_y[k] = -g_0_z_yyzzzz_y[k] * ab_x + g_0_z_yyzzzz_xy[k];

                g_0_z_xyyzzzz_z[k] = -g_0_z_yyzzzz_z[k] * ab_x + g_0_z_yyzzzz_xz[k];
            }

            /// Set up 294-297 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzzzz_x = cbuffer.data(kp_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_y = cbuffer.data(kp_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_z = cbuffer.data(kp_geom_01_off + 296 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzzzz_x, g_0_z_xyzzzzz_y, g_0_z_xyzzzzz_z, g_0_z_yzzzzz_x, g_0_z_yzzzzz_xx, g_0_z_yzzzzz_xy, g_0_z_yzzzzz_xz, g_0_z_yzzzzz_y, g_0_z_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzzzz_x[k] = -g_0_z_yzzzzz_x[k] * ab_x + g_0_z_yzzzzz_xx[k];

                g_0_z_xyzzzzz_y[k] = -g_0_z_yzzzzz_y[k] * ab_x + g_0_z_yzzzzz_xy[k];

                g_0_z_xyzzzzz_z[k] = -g_0_z_yzzzzz_z[k] * ab_x + g_0_z_yzzzzz_xz[k];
            }

            /// Set up 297-300 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzzzz_x = cbuffer.data(kp_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_y = cbuffer.data(kp_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_z = cbuffer.data(kp_geom_01_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzzzz_x, g_0_z_xzzzzzz_y, g_0_z_xzzzzzz_z, g_0_z_zzzzzz_x, g_0_z_zzzzzz_xx, g_0_z_zzzzzz_xy, g_0_z_zzzzzz_xz, g_0_z_zzzzzz_y, g_0_z_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzzzz_x[k] = -g_0_z_zzzzzz_x[k] * ab_x + g_0_z_zzzzzz_xx[k];

                g_0_z_xzzzzzz_y[k] = -g_0_z_zzzzzz_y[k] * ab_x + g_0_z_zzzzzz_xy[k];

                g_0_z_xzzzzzz_z[k] = -g_0_z_zzzzzz_z[k] * ab_x + g_0_z_zzzzzz_xz[k];
            }

            /// Set up 300-303 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyyy_x = cbuffer.data(kp_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_y = cbuffer.data(kp_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_z = cbuffer.data(kp_geom_01_off + 302 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyy_x, g_0_z_yyyyyy_xy, g_0_z_yyyyyy_y, g_0_z_yyyyyy_yy, g_0_z_yyyyyy_yz, g_0_z_yyyyyy_z, g_0_z_yyyyyyy_x, g_0_z_yyyyyyy_y, g_0_z_yyyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyyy_x[k] = -g_0_z_yyyyyy_x[k] * ab_y + g_0_z_yyyyyy_xy[k];

                g_0_z_yyyyyyy_y[k] = -g_0_z_yyyyyy_y[k] * ab_y + g_0_z_yyyyyy_yy[k];

                g_0_z_yyyyyyy_z[k] = -g_0_z_yyyyyy_z[k] * ab_y + g_0_z_yyyyyy_yz[k];
            }

            /// Set up 303-306 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyyz_x = cbuffer.data(kp_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_y = cbuffer.data(kp_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_z = cbuffer.data(kp_geom_01_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyyz_x, g_0_z_yyyyyyz_y, g_0_z_yyyyyyz_z, g_0_z_yyyyyz_x, g_0_z_yyyyyz_xy, g_0_z_yyyyyz_y, g_0_z_yyyyyz_yy, g_0_z_yyyyyz_yz, g_0_z_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyyz_x[k] = -g_0_z_yyyyyz_x[k] * ab_y + g_0_z_yyyyyz_xy[k];

                g_0_z_yyyyyyz_y[k] = -g_0_z_yyyyyz_y[k] * ab_y + g_0_z_yyyyyz_yy[k];

                g_0_z_yyyyyyz_z[k] = -g_0_z_yyyyyz_z[k] * ab_y + g_0_z_yyyyyz_yz[k];
            }

            /// Set up 306-309 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyzz_x = cbuffer.data(kp_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_y = cbuffer.data(kp_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_z = cbuffer.data(kp_geom_01_off + 308 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyzz_x, g_0_z_yyyyyzz_y, g_0_z_yyyyyzz_z, g_0_z_yyyyzz_x, g_0_z_yyyyzz_xy, g_0_z_yyyyzz_y, g_0_z_yyyyzz_yy, g_0_z_yyyyzz_yz, g_0_z_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyzz_x[k] = -g_0_z_yyyyzz_x[k] * ab_y + g_0_z_yyyyzz_xy[k];

                g_0_z_yyyyyzz_y[k] = -g_0_z_yyyyzz_y[k] * ab_y + g_0_z_yyyyzz_yy[k];

                g_0_z_yyyyyzz_z[k] = -g_0_z_yyyyzz_z[k] * ab_y + g_0_z_yyyyzz_yz[k];
            }

            /// Set up 309-312 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyzzz_x = cbuffer.data(kp_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_y = cbuffer.data(kp_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_z = cbuffer.data(kp_geom_01_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyzzz_x, g_0_z_yyyyzzz_y, g_0_z_yyyyzzz_z, g_0_z_yyyzzz_x, g_0_z_yyyzzz_xy, g_0_z_yyyzzz_y, g_0_z_yyyzzz_yy, g_0_z_yyyzzz_yz, g_0_z_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyzzz_x[k] = -g_0_z_yyyzzz_x[k] * ab_y + g_0_z_yyyzzz_xy[k];

                g_0_z_yyyyzzz_y[k] = -g_0_z_yyyzzz_y[k] * ab_y + g_0_z_yyyzzz_yy[k];

                g_0_z_yyyyzzz_z[k] = -g_0_z_yyyzzz_z[k] * ab_y + g_0_z_yyyzzz_yz[k];
            }

            /// Set up 312-315 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzzzz_x = cbuffer.data(kp_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_y = cbuffer.data(kp_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_z = cbuffer.data(kp_geom_01_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzzzz_x, g_0_z_yyyzzzz_y, g_0_z_yyyzzzz_z, g_0_z_yyzzzz_x, g_0_z_yyzzzz_xy, g_0_z_yyzzzz_y, g_0_z_yyzzzz_yy, g_0_z_yyzzzz_yz, g_0_z_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzzzz_x[k] = -g_0_z_yyzzzz_x[k] * ab_y + g_0_z_yyzzzz_xy[k];

                g_0_z_yyyzzzz_y[k] = -g_0_z_yyzzzz_y[k] * ab_y + g_0_z_yyzzzz_yy[k];

                g_0_z_yyyzzzz_z[k] = -g_0_z_yyzzzz_z[k] * ab_y + g_0_z_yyzzzz_yz[k];
            }

            /// Set up 315-318 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzzzz_x = cbuffer.data(kp_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_y = cbuffer.data(kp_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_z = cbuffer.data(kp_geom_01_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzzzz_x, g_0_z_yyzzzzz_y, g_0_z_yyzzzzz_z, g_0_z_yzzzzz_x, g_0_z_yzzzzz_xy, g_0_z_yzzzzz_y, g_0_z_yzzzzz_yy, g_0_z_yzzzzz_yz, g_0_z_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzzzz_x[k] = -g_0_z_yzzzzz_x[k] * ab_y + g_0_z_yzzzzz_xy[k];

                g_0_z_yyzzzzz_y[k] = -g_0_z_yzzzzz_y[k] * ab_y + g_0_z_yzzzzz_yy[k];

                g_0_z_yyzzzzz_z[k] = -g_0_z_yzzzzz_z[k] * ab_y + g_0_z_yzzzzz_yz[k];
            }

            /// Set up 318-321 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzzzz_x = cbuffer.data(kp_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_y = cbuffer.data(kp_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_z = cbuffer.data(kp_geom_01_off + 320 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzzzz_x, g_0_z_yzzzzzz_y, g_0_z_yzzzzzz_z, g_0_z_zzzzzz_x, g_0_z_zzzzzz_xy, g_0_z_zzzzzz_y, g_0_z_zzzzzz_yy, g_0_z_zzzzzz_yz, g_0_z_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzzzz_x[k] = -g_0_z_zzzzzz_x[k] * ab_y + g_0_z_zzzzzz_xy[k];

                g_0_z_yzzzzzz_y[k] = -g_0_z_zzzzzz_y[k] * ab_y + g_0_z_zzzzzz_yy[k];

                g_0_z_yzzzzzz_z[k] = -g_0_z_zzzzzz_z[k] * ab_y + g_0_z_zzzzzz_yz[k];
            }

            /// Set up 321-324 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzzzz_x = cbuffer.data(kp_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_y = cbuffer.data(kp_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_z = cbuffer.data(kp_geom_01_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzzzz_x, g_0_z_zzzzzz_xz, g_0_z_zzzzzz_y, g_0_z_zzzzzz_yz, g_0_z_zzzzzz_z, g_0_z_zzzzzz_zz, g_0_z_zzzzzzz_x, g_0_z_zzzzzzz_y, g_0_z_zzzzzzz_z, g_zzzzzz_x, g_zzzzzz_y, g_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzzzz_x[k] = g_zzzzzz_x[k] - g_0_z_zzzzzz_x[k] * ab_z + g_0_z_zzzzzz_xz[k];

                g_0_z_zzzzzzz_y[k] = g_zzzzzz_y[k] - g_0_z_zzzzzz_y[k] * ab_z + g_0_z_zzzzzz_yz[k];

                g_0_z_zzzzzzz_z[k] = g_zzzzzz_z[k] - g_0_z_zzzzzz_z[k] * ab_z + g_0_z_zzzzzz_zz[k];
            }
        }
    }
}

} // erirec namespace

