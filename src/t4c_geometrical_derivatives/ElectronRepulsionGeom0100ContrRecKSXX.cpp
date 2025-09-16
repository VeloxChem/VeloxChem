#include "ElectronRepulsionGeom0100ContrRecKSXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_ksxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_ksxx,
                                            const size_t idx_isxx,
                                            const size_t idx_geom_01_isxx,
                                            const size_t idx_geom_01_ipxx,
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
            /// Set up components of auxilary buffer : ISSS

            const auto is_off = idx_isxx + i * dcomps + j;

            auto g_xxxxxx_0 = cbuffer.data(is_off + 0 * ccomps * dcomps);

            auto g_xxxxxy_0 = cbuffer.data(is_off + 1 * ccomps * dcomps);

            auto g_xxxxxz_0 = cbuffer.data(is_off + 2 * ccomps * dcomps);

            auto g_xxxxyy_0 = cbuffer.data(is_off + 3 * ccomps * dcomps);

            auto g_xxxxyz_0 = cbuffer.data(is_off + 4 * ccomps * dcomps);

            auto g_xxxxzz_0 = cbuffer.data(is_off + 5 * ccomps * dcomps);

            auto g_xxxyyy_0 = cbuffer.data(is_off + 6 * ccomps * dcomps);

            auto g_xxxyyz_0 = cbuffer.data(is_off + 7 * ccomps * dcomps);

            auto g_xxxyzz_0 = cbuffer.data(is_off + 8 * ccomps * dcomps);

            auto g_xxxzzz_0 = cbuffer.data(is_off + 9 * ccomps * dcomps);

            auto g_xxyyyy_0 = cbuffer.data(is_off + 10 * ccomps * dcomps);

            auto g_xxyyyz_0 = cbuffer.data(is_off + 11 * ccomps * dcomps);

            auto g_xxyyzz_0 = cbuffer.data(is_off + 12 * ccomps * dcomps);

            auto g_xxyzzz_0 = cbuffer.data(is_off + 13 * ccomps * dcomps);

            auto g_xxzzzz_0 = cbuffer.data(is_off + 14 * ccomps * dcomps);

            auto g_xyyyyy_0 = cbuffer.data(is_off + 15 * ccomps * dcomps);

            auto g_xyyyyz_0 = cbuffer.data(is_off + 16 * ccomps * dcomps);

            auto g_xyyyzz_0 = cbuffer.data(is_off + 17 * ccomps * dcomps);

            auto g_xyyzzz_0 = cbuffer.data(is_off + 18 * ccomps * dcomps);

            auto g_xyzzzz_0 = cbuffer.data(is_off + 19 * ccomps * dcomps);

            auto g_xzzzzz_0 = cbuffer.data(is_off + 20 * ccomps * dcomps);

            auto g_yyyyyy_0 = cbuffer.data(is_off + 21 * ccomps * dcomps);

            auto g_yyyyyz_0 = cbuffer.data(is_off + 22 * ccomps * dcomps);

            auto g_yyyyzz_0 = cbuffer.data(is_off + 23 * ccomps * dcomps);

            auto g_yyyzzz_0 = cbuffer.data(is_off + 24 * ccomps * dcomps);

            auto g_yyzzzz_0 = cbuffer.data(is_off + 25 * ccomps * dcomps);

            auto g_yzzzzz_0 = cbuffer.data(is_off + 26 * ccomps * dcomps);

            auto g_zzzzzz_0 = cbuffer.data(is_off + 27 * ccomps * dcomps);

            /// Set up components of auxilary buffer : ISSS

            const auto is_geom_01_off = idx_geom_01_isxx + i * dcomps + j;

            auto g_0_x_xxxxxx_0 = cbuffer.data(is_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxy_0 = cbuffer.data(is_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxz_0 = cbuffer.data(is_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxyy_0 = cbuffer.data(is_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxyz_0 = cbuffer.data(is_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxzz_0 = cbuffer.data(is_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxyyy_0 = cbuffer.data(is_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxyyz_0 = cbuffer.data(is_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxyzz_0 = cbuffer.data(is_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxzzz_0 = cbuffer.data(is_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxyyyy_0 = cbuffer.data(is_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxyyyz_0 = cbuffer.data(is_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxyyzz_0 = cbuffer.data(is_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxyzzz_0 = cbuffer.data(is_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxzzzz_0 = cbuffer.data(is_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xyyyyy_0 = cbuffer.data(is_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xyyyyz_0 = cbuffer.data(is_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xyyyzz_0 = cbuffer.data(is_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xyyzzz_0 = cbuffer.data(is_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xyzzzz_0 = cbuffer.data(is_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xzzzzz_0 = cbuffer.data(is_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_yyyyyy_0 = cbuffer.data(is_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_yyyyyz_0 = cbuffer.data(is_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_yyyyzz_0 = cbuffer.data(is_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_yyyzzz_0 = cbuffer.data(is_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_yyzzzz_0 = cbuffer.data(is_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_yzzzzz_0 = cbuffer.data(is_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_zzzzzz_0 = cbuffer.data(is_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_y_xxxxxx_0 = cbuffer.data(is_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_y_xxxxxy_0 = cbuffer.data(is_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_y_xxxxxz_0 = cbuffer.data(is_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_y_xxxxyy_0 = cbuffer.data(is_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_y_xxxxyz_0 = cbuffer.data(is_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_y_xxxxzz_0 = cbuffer.data(is_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_y_xxxyyy_0 = cbuffer.data(is_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_y_xxxyyz_0 = cbuffer.data(is_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_y_xxxyzz_0 = cbuffer.data(is_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_y_xxxzzz_0 = cbuffer.data(is_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_y_xxyyyy_0 = cbuffer.data(is_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_y_xxyyyz_0 = cbuffer.data(is_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_y_xxyyzz_0 = cbuffer.data(is_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_y_xxyzzz_0 = cbuffer.data(is_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_y_xxzzzz_0 = cbuffer.data(is_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_y_xyyyyy_0 = cbuffer.data(is_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_y_xyyyyz_0 = cbuffer.data(is_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_y_xyyyzz_0 = cbuffer.data(is_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_y_xyyzzz_0 = cbuffer.data(is_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_y_xyzzzz_0 = cbuffer.data(is_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_y_xzzzzz_0 = cbuffer.data(is_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_y_yyyyyy_0 = cbuffer.data(is_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_y_yyyyyz_0 = cbuffer.data(is_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_y_yyyyzz_0 = cbuffer.data(is_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_y_yyyzzz_0 = cbuffer.data(is_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_y_yyzzzz_0 = cbuffer.data(is_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_y_yzzzzz_0 = cbuffer.data(is_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_y_zzzzzz_0 = cbuffer.data(is_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_z_xxxxxx_0 = cbuffer.data(is_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_z_xxxxxy_0 = cbuffer.data(is_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_z_xxxxxz_0 = cbuffer.data(is_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_z_xxxxyy_0 = cbuffer.data(is_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_z_xxxxyz_0 = cbuffer.data(is_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_z_xxxxzz_0 = cbuffer.data(is_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_z_xxxyyy_0 = cbuffer.data(is_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_z_xxxyyz_0 = cbuffer.data(is_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_z_xxxyzz_0 = cbuffer.data(is_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_z_xxxzzz_0 = cbuffer.data(is_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_z_xxyyyy_0 = cbuffer.data(is_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_z_xxyyyz_0 = cbuffer.data(is_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_z_xxyyzz_0 = cbuffer.data(is_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_z_xxyzzz_0 = cbuffer.data(is_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_z_xxzzzz_0 = cbuffer.data(is_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_z_xyyyyy_0 = cbuffer.data(is_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_z_xyyyyz_0 = cbuffer.data(is_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_z_xyyyzz_0 = cbuffer.data(is_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_z_xyyzzz_0 = cbuffer.data(is_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_z_xyzzzz_0 = cbuffer.data(is_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_z_xzzzzz_0 = cbuffer.data(is_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_z_yyyyyy_0 = cbuffer.data(is_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_z_yyyyyz_0 = cbuffer.data(is_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_z_yyyyzz_0 = cbuffer.data(is_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_z_yyyzzz_0 = cbuffer.data(is_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_z_yyzzzz_0 = cbuffer.data(is_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_z_yzzzzz_0 = cbuffer.data(is_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_z_zzzzzz_0 = cbuffer.data(is_geom_01_off + 83 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_ksxx

            const auto ks_geom_01_off = idx_geom_01_ksxx + i * dcomps + j;

            /// Set up 0-1 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxx_0 = cbuffer.data(ks_geom_01_off + 0 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_0, g_0_x_xxxxxx_x, g_0_x_xxxxxxx_0, g_xxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxx_0[k] = g_xxxxxx_0[k] - g_0_x_xxxxxx_0[k] * ab_x + g_0_x_xxxxxx_x[k];
            }

            /// Set up 1-2 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxy_0 = cbuffer.data(ks_geom_01_off + 1 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_0, g_0_x_xxxxxx_y, g_0_x_xxxxxxy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxy_0[k] = -g_0_x_xxxxxx_0[k] * ab_y + g_0_x_xxxxxx_y[k];
            }

            /// Set up 2-3 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxz_0 = cbuffer.data(ks_geom_01_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_0, g_0_x_xxxxxx_z, g_0_x_xxxxxxz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxz_0[k] = -g_0_x_xxxxxx_0[k] * ab_z + g_0_x_xxxxxx_z[k];
            }

            /// Set up 3-4 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxyy_0 = cbuffer.data(ks_geom_01_off + 3 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxy_0, g_0_x_xxxxxy_y, g_0_x_xxxxxyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxyy_0[k] = -g_0_x_xxxxxy_0[k] * ab_y + g_0_x_xxxxxy_y[k];
            }

            /// Set up 4-5 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxyz_0 = cbuffer.data(ks_geom_01_off + 4 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxyz_0, g_0_x_xxxxxz_0, g_0_x_xxxxxz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxyz_0[k] = -g_0_x_xxxxxz_0[k] * ab_y + g_0_x_xxxxxz_y[k];
            }

            /// Set up 5-6 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxzz_0 = cbuffer.data(ks_geom_01_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxz_0, g_0_x_xxxxxz_z, g_0_x_xxxxxzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxzz_0[k] = -g_0_x_xxxxxz_0[k] * ab_z + g_0_x_xxxxxz_z[k];
            }

            /// Set up 6-7 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyyy_0 = cbuffer.data(ks_geom_01_off + 6 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyy_0, g_0_x_xxxxyy_y, g_0_x_xxxxyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyyy_0[k] = -g_0_x_xxxxyy_0[k] * ab_y + g_0_x_xxxxyy_y[k];
            }

            /// Set up 7-8 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyyz_0 = cbuffer.data(ks_geom_01_off + 7 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyyz_0, g_0_x_xxxxyz_0, g_0_x_xxxxyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyyz_0[k] = -g_0_x_xxxxyz_0[k] * ab_y + g_0_x_xxxxyz_y[k];
            }

            /// Set up 8-9 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyzz_0 = cbuffer.data(ks_geom_01_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyzz_0, g_0_x_xxxxzz_0, g_0_x_xxxxzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyzz_0[k] = -g_0_x_xxxxzz_0[k] * ab_y + g_0_x_xxxxzz_y[k];
            }

            /// Set up 9-10 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxzzz_0 = cbuffer.data(ks_geom_01_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxzz_0, g_0_x_xxxxzz_z, g_0_x_xxxxzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxzzz_0[k] = -g_0_x_xxxxzz_0[k] * ab_z + g_0_x_xxxxzz_z[k];
            }

            /// Set up 10-11 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyyy_0 = cbuffer.data(ks_geom_01_off + 10 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyy_0, g_0_x_xxxyyy_y, g_0_x_xxxyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyyy_0[k] = -g_0_x_xxxyyy_0[k] * ab_y + g_0_x_xxxyyy_y[k];
            }

            /// Set up 11-12 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyyz_0 = cbuffer.data(ks_geom_01_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyyz_0, g_0_x_xxxyyz_0, g_0_x_xxxyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyyz_0[k] = -g_0_x_xxxyyz_0[k] * ab_y + g_0_x_xxxyyz_y[k];
            }

            /// Set up 12-13 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyzz_0 = cbuffer.data(ks_geom_01_off + 12 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyzz_0, g_0_x_xxxyzz_0, g_0_x_xxxyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyzz_0[k] = -g_0_x_xxxyzz_0[k] * ab_y + g_0_x_xxxyzz_y[k];
            }

            /// Set up 13-14 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyzzz_0 = cbuffer.data(ks_geom_01_off + 13 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyzzz_0, g_0_x_xxxzzz_0, g_0_x_xxxzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyzzz_0[k] = -g_0_x_xxxzzz_0[k] * ab_y + g_0_x_xxxzzz_y[k];
            }

            /// Set up 14-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzzzz_0 = cbuffer.data(ks_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxzzz_0, g_0_x_xxxzzz_z, g_0_x_xxxzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzzzz_0[k] = -g_0_x_xxxzzz_0[k] * ab_z + g_0_x_xxxzzz_z[k];
            }

            /// Set up 15-16 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyyy_0 = cbuffer.data(ks_geom_01_off + 15 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyy_0, g_0_x_xxyyyy_y, g_0_x_xxyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyyy_0[k] = -g_0_x_xxyyyy_0[k] * ab_y + g_0_x_xxyyyy_y[k];
            }

            /// Set up 16-17 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyyz_0 = cbuffer.data(ks_geom_01_off + 16 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyyz_0, g_0_x_xxyyyz_0, g_0_x_xxyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyyz_0[k] = -g_0_x_xxyyyz_0[k] * ab_y + g_0_x_xxyyyz_y[k];
            }

            /// Set up 17-18 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyzz_0 = cbuffer.data(ks_geom_01_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyzz_0, g_0_x_xxyyzz_0, g_0_x_xxyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyzz_0[k] = -g_0_x_xxyyzz_0[k] * ab_y + g_0_x_xxyyzz_y[k];
            }

            /// Set up 18-19 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyzzz_0 = cbuffer.data(ks_geom_01_off + 18 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyzzz_0, g_0_x_xxyzzz_0, g_0_x_xxyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyzzz_0[k] = -g_0_x_xxyzzz_0[k] * ab_y + g_0_x_xxyzzz_y[k];
            }

            /// Set up 19-20 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzzzz_0 = cbuffer.data(ks_geom_01_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzzzz_0, g_0_x_xxzzzz_0, g_0_x_xxzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzzzz_0[k] = -g_0_x_xxzzzz_0[k] * ab_y + g_0_x_xxzzzz_y[k];
            }

            /// Set up 20-21 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzzzz_0 = cbuffer.data(ks_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzzzz_0, g_0_x_xxzzzz_z, g_0_x_xxzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzzzz_0[k] = -g_0_x_xxzzzz_0[k] * ab_z + g_0_x_xxzzzz_z[k];
            }

            /// Set up 21-22 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyyy_0 = cbuffer.data(ks_geom_01_off + 21 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyy_0, g_0_x_xyyyyy_y, g_0_x_xyyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyyy_0[k] = -g_0_x_xyyyyy_0[k] * ab_y + g_0_x_xyyyyy_y[k];
            }

            /// Set up 22-23 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyyz_0 = cbuffer.data(ks_geom_01_off + 22 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyyz_0, g_0_x_xyyyyz_0, g_0_x_xyyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyyz_0[k] = -g_0_x_xyyyyz_0[k] * ab_y + g_0_x_xyyyyz_y[k];
            }

            /// Set up 23-24 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyzz_0 = cbuffer.data(ks_geom_01_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyzz_0, g_0_x_xyyyzz_0, g_0_x_xyyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyzz_0[k] = -g_0_x_xyyyzz_0[k] * ab_y + g_0_x_xyyyzz_y[k];
            }

            /// Set up 24-25 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyzzz_0 = cbuffer.data(ks_geom_01_off + 24 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyzzz_0, g_0_x_xyyzzz_0, g_0_x_xyyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyzzz_0[k] = -g_0_x_xyyzzz_0[k] * ab_y + g_0_x_xyyzzz_y[k];
            }

            /// Set up 25-26 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzzzz_0 = cbuffer.data(ks_geom_01_off + 25 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzzzz_0, g_0_x_xyzzzz_0, g_0_x_xyzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzzzz_0[k] = -g_0_x_xyzzzz_0[k] * ab_y + g_0_x_xyzzzz_y[k];
            }

            /// Set up 26-27 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzzzz_0 = cbuffer.data(ks_geom_01_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzzzz_0, g_0_x_xzzzzz_0, g_0_x_xzzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzzzz_0[k] = -g_0_x_xzzzzz_0[k] * ab_y + g_0_x_xzzzzz_y[k];
            }

            /// Set up 27-28 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzzzz_0 = cbuffer.data(ks_geom_01_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzzzz_0, g_0_x_xzzzzz_z, g_0_x_xzzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzzzz_0[k] = -g_0_x_xzzzzz_0[k] * ab_z + g_0_x_xzzzzz_z[k];
            }

            /// Set up 28-29 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyyy_0 = cbuffer.data(ks_geom_01_off + 28 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyy_0, g_0_x_yyyyyy_y, g_0_x_yyyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyyy_0[k] = -g_0_x_yyyyyy_0[k] * ab_y + g_0_x_yyyyyy_y[k];
            }

            /// Set up 29-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyyz_0 = cbuffer.data(ks_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyyz_0, g_0_x_yyyyyz_0, g_0_x_yyyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyyz_0[k] = -g_0_x_yyyyyz_0[k] * ab_y + g_0_x_yyyyyz_y[k];
            }

            /// Set up 30-31 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyzz_0 = cbuffer.data(ks_geom_01_off + 30 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyzz_0, g_0_x_yyyyzz_0, g_0_x_yyyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyzz_0[k] = -g_0_x_yyyyzz_0[k] * ab_y + g_0_x_yyyyzz_y[k];
            }

            /// Set up 31-32 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyzzz_0 = cbuffer.data(ks_geom_01_off + 31 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyzzz_0, g_0_x_yyyzzz_0, g_0_x_yyyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyzzz_0[k] = -g_0_x_yyyzzz_0[k] * ab_y + g_0_x_yyyzzz_y[k];
            }

            /// Set up 32-33 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzzzz_0 = cbuffer.data(ks_geom_01_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzzzz_0, g_0_x_yyzzzz_0, g_0_x_yyzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzzzz_0[k] = -g_0_x_yyzzzz_0[k] * ab_y + g_0_x_yyzzzz_y[k];
            }

            /// Set up 33-34 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzzzz_0 = cbuffer.data(ks_geom_01_off + 33 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzzzz_0, g_0_x_yzzzzz_0, g_0_x_yzzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzzzz_0[k] = -g_0_x_yzzzzz_0[k] * ab_y + g_0_x_yzzzzz_y[k];
            }

            /// Set up 34-35 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzzzz_0 = cbuffer.data(ks_geom_01_off + 34 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzzzz_0, g_0_x_zzzzzz_0, g_0_x_zzzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzzzz_0[k] = -g_0_x_zzzzzz_0[k] * ab_y + g_0_x_zzzzzz_y[k];
            }

            /// Set up 35-36 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzzzz_0 = cbuffer.data(ks_geom_01_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzzzz_0, g_0_x_zzzzzz_z, g_0_x_zzzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzzzz_0[k] = -g_0_x_zzzzzz_0[k] * ab_z + g_0_x_zzzzzz_z[k];
            }

            /// Set up 36-37 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxx_0 = cbuffer.data(ks_geom_01_off + 36 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxx_0, g_0_y_xxxxxx_x, g_0_y_xxxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxx_0[k] = -g_0_y_xxxxxx_0[k] * ab_x + g_0_y_xxxxxx_x[k];
            }

            /// Set up 37-38 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxy_0 = cbuffer.data(ks_geom_01_off + 37 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxxy_0, g_0_y_xxxxxy_0, g_0_y_xxxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxy_0[k] = -g_0_y_xxxxxy_0[k] * ab_x + g_0_y_xxxxxy_x[k];
            }

            /// Set up 38-39 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxz_0 = cbuffer.data(ks_geom_01_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxxz_0, g_0_y_xxxxxz_0, g_0_y_xxxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxz_0[k] = -g_0_y_xxxxxz_0[k] * ab_x + g_0_y_xxxxxz_x[k];
            }

            /// Set up 39-40 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxyy_0 = cbuffer.data(ks_geom_01_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxyy_0, g_0_y_xxxxyy_0, g_0_y_xxxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxyy_0[k] = -g_0_y_xxxxyy_0[k] * ab_x + g_0_y_xxxxyy_x[k];
            }

            /// Set up 40-41 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxyz_0 = cbuffer.data(ks_geom_01_off + 40 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxyz_0, g_0_y_xxxxyz_0, g_0_y_xxxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxyz_0[k] = -g_0_y_xxxxyz_0[k] * ab_x + g_0_y_xxxxyz_x[k];
            }

            /// Set up 41-42 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxzz_0 = cbuffer.data(ks_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxzz_0, g_0_y_xxxxzz_0, g_0_y_xxxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxzz_0[k] = -g_0_y_xxxxzz_0[k] * ab_x + g_0_y_xxxxzz_x[k];
            }

            /// Set up 42-43 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyyy_0 = cbuffer.data(ks_geom_01_off + 42 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyyy_0, g_0_y_xxxyyy_0, g_0_y_xxxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyyy_0[k] = -g_0_y_xxxyyy_0[k] * ab_x + g_0_y_xxxyyy_x[k];
            }

            /// Set up 43-44 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyyz_0 = cbuffer.data(ks_geom_01_off + 43 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyyz_0, g_0_y_xxxyyz_0, g_0_y_xxxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyyz_0[k] = -g_0_y_xxxyyz_0[k] * ab_x + g_0_y_xxxyyz_x[k];
            }

            /// Set up 44-45 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyzz_0 = cbuffer.data(ks_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyzz_0, g_0_y_xxxyzz_0, g_0_y_xxxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyzz_0[k] = -g_0_y_xxxyzz_0[k] * ab_x + g_0_y_xxxyzz_x[k];
            }

            /// Set up 45-46 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxzzz_0 = cbuffer.data(ks_geom_01_off + 45 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxzzz_0, g_0_y_xxxzzz_0, g_0_y_xxxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxzzz_0[k] = -g_0_y_xxxzzz_0[k] * ab_x + g_0_y_xxxzzz_x[k];
            }

            /// Set up 46-47 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyyy_0 = cbuffer.data(ks_geom_01_off + 46 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyyy_0, g_0_y_xxyyyy_0, g_0_y_xxyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyyy_0[k] = -g_0_y_xxyyyy_0[k] * ab_x + g_0_y_xxyyyy_x[k];
            }

            /// Set up 47-48 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyyz_0 = cbuffer.data(ks_geom_01_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyyz_0, g_0_y_xxyyyz_0, g_0_y_xxyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyyz_0[k] = -g_0_y_xxyyyz_0[k] * ab_x + g_0_y_xxyyyz_x[k];
            }

            /// Set up 48-49 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyzz_0 = cbuffer.data(ks_geom_01_off + 48 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyzz_0, g_0_y_xxyyzz_0, g_0_y_xxyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyzz_0[k] = -g_0_y_xxyyzz_0[k] * ab_x + g_0_y_xxyyzz_x[k];
            }

            /// Set up 49-50 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyzzz_0 = cbuffer.data(ks_geom_01_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyzzz_0, g_0_y_xxyzzz_0, g_0_y_xxyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyzzz_0[k] = -g_0_y_xxyzzz_0[k] * ab_x + g_0_y_xxyzzz_x[k];
            }

            /// Set up 50-51 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzzzz_0 = cbuffer.data(ks_geom_01_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzzzz_0, g_0_y_xxzzzz_0, g_0_y_xxzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzzzz_0[k] = -g_0_y_xxzzzz_0[k] * ab_x + g_0_y_xxzzzz_x[k];
            }

            /// Set up 51-52 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyyy_0 = cbuffer.data(ks_geom_01_off + 51 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyyy_0, g_0_y_xyyyyy_0, g_0_y_xyyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyyy_0[k] = -g_0_y_xyyyyy_0[k] * ab_x + g_0_y_xyyyyy_x[k];
            }

            /// Set up 52-53 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyyz_0 = cbuffer.data(ks_geom_01_off + 52 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyyz_0, g_0_y_xyyyyz_0, g_0_y_xyyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyyz_0[k] = -g_0_y_xyyyyz_0[k] * ab_x + g_0_y_xyyyyz_x[k];
            }

            /// Set up 53-54 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyzz_0 = cbuffer.data(ks_geom_01_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyzz_0, g_0_y_xyyyzz_0, g_0_y_xyyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyzz_0[k] = -g_0_y_xyyyzz_0[k] * ab_x + g_0_y_xyyyzz_x[k];
            }

            /// Set up 54-55 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyzzz_0 = cbuffer.data(ks_geom_01_off + 54 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyzzz_0, g_0_y_xyyzzz_0, g_0_y_xyyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyzzz_0[k] = -g_0_y_xyyzzz_0[k] * ab_x + g_0_y_xyyzzz_x[k];
            }

            /// Set up 55-56 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzzzz_0 = cbuffer.data(ks_geom_01_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzzzz_0, g_0_y_xyzzzz_0, g_0_y_xyzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzzzz_0[k] = -g_0_y_xyzzzz_0[k] * ab_x + g_0_y_xyzzzz_x[k];
            }

            /// Set up 56-57 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzzzz_0 = cbuffer.data(ks_geom_01_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzzzz_0, g_0_y_xzzzzz_0, g_0_y_xzzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzzzz_0[k] = -g_0_y_xzzzzz_0[k] * ab_x + g_0_y_xzzzzz_x[k];
            }

            /// Set up 57-58 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyyy_0 = cbuffer.data(ks_geom_01_off + 57 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyyy_0, g_0_y_yyyyyy_0, g_0_y_yyyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyyy_0[k] = -g_0_y_yyyyyy_0[k] * ab_x + g_0_y_yyyyyy_x[k];
            }

            /// Set up 58-59 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyyz_0 = cbuffer.data(ks_geom_01_off + 58 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyyz_0, g_0_y_yyyyyz_0, g_0_y_yyyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyyz_0[k] = -g_0_y_yyyyyz_0[k] * ab_x + g_0_y_yyyyyz_x[k];
            }

            /// Set up 59-60 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyzz_0 = cbuffer.data(ks_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyzz_0, g_0_y_yyyyzz_0, g_0_y_yyyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyzz_0[k] = -g_0_y_yyyyzz_0[k] * ab_x + g_0_y_yyyyzz_x[k];
            }

            /// Set up 60-61 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyzzz_0 = cbuffer.data(ks_geom_01_off + 60 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyzzz_0, g_0_y_yyyzzz_0, g_0_y_yyyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyzzz_0[k] = -g_0_y_yyyzzz_0[k] * ab_x + g_0_y_yyyzzz_x[k];
            }

            /// Set up 61-62 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzzzz_0 = cbuffer.data(ks_geom_01_off + 61 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzzzz_0, g_0_y_yyzzzz_0, g_0_y_yyzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzzzz_0[k] = -g_0_y_yyzzzz_0[k] * ab_x + g_0_y_yyzzzz_x[k];
            }

            /// Set up 62-63 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzzzz_0 = cbuffer.data(ks_geom_01_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzzzz_0, g_0_y_yzzzzz_0, g_0_y_yzzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzzzz_0[k] = -g_0_y_yzzzzz_0[k] * ab_x + g_0_y_yzzzzz_x[k];
            }

            /// Set up 63-64 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzzzz_0 = cbuffer.data(ks_geom_01_off + 63 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzzzz_0, g_0_y_zzzzzz_0, g_0_y_zzzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzzzz_0[k] = -g_0_y_zzzzzz_0[k] * ab_x + g_0_y_zzzzzz_x[k];
            }

            /// Set up 64-65 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyyy_0 = cbuffer.data(ks_geom_01_off + 64 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyy_0, g_0_y_yyyyyy_y, g_0_y_yyyyyyy_0, g_yyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyyy_0[k] = g_yyyyyy_0[k] - g_0_y_yyyyyy_0[k] * ab_y + g_0_y_yyyyyy_y[k];
            }

            /// Set up 65-66 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyyz_0 = cbuffer.data(ks_geom_01_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyy_0, g_0_y_yyyyyy_z, g_0_y_yyyyyyz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyyz_0[k] = -g_0_y_yyyyyy_0[k] * ab_z + g_0_y_yyyyyy_z[k];
            }

            /// Set up 66-67 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyzz_0 = cbuffer.data(ks_geom_01_off + 66 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyz_0, g_0_y_yyyyyz_z, g_0_y_yyyyyzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyzz_0[k] = -g_0_y_yyyyyz_0[k] * ab_z + g_0_y_yyyyyz_z[k];
            }

            /// Set up 67-68 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyzzz_0 = cbuffer.data(ks_geom_01_off + 67 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyzz_0, g_0_y_yyyyzz_z, g_0_y_yyyyzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyzzz_0[k] = -g_0_y_yyyyzz_0[k] * ab_z + g_0_y_yyyyzz_z[k];
            }

            /// Set up 68-69 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzzzz_0 = cbuffer.data(ks_geom_01_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyzzz_0, g_0_y_yyyzzz_z, g_0_y_yyyzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzzzz_0[k] = -g_0_y_yyyzzz_0[k] * ab_z + g_0_y_yyyzzz_z[k];
            }

            /// Set up 69-70 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzzzz_0 = cbuffer.data(ks_geom_01_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzzzz_0, g_0_y_yyzzzz_z, g_0_y_yyzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzzzz_0[k] = -g_0_y_yyzzzz_0[k] * ab_z + g_0_y_yyzzzz_z[k];
            }

            /// Set up 70-71 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzzzz_0 = cbuffer.data(ks_geom_01_off + 70 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzzzz_0, g_0_y_yzzzzz_z, g_0_y_yzzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzzzz_0[k] = -g_0_y_yzzzzz_0[k] * ab_z + g_0_y_yzzzzz_z[k];
            }

            /// Set up 71-72 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzzzz_0 = cbuffer.data(ks_geom_01_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzzzz_0, g_0_y_zzzzzz_z, g_0_y_zzzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzzzz_0[k] = -g_0_y_zzzzzz_0[k] * ab_z + g_0_y_zzzzzz_z[k];
            }

            /// Set up 72-73 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxx_0 = cbuffer.data(ks_geom_01_off + 72 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxx_0, g_0_z_xxxxxx_x, g_0_z_xxxxxxx_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxx_0[k] = -g_0_z_xxxxxx_0[k] * ab_x + g_0_z_xxxxxx_x[k];
            }

            /// Set up 73-74 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxy_0 = cbuffer.data(ks_geom_01_off + 73 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxxy_0, g_0_z_xxxxxy_0, g_0_z_xxxxxy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxy_0[k] = -g_0_z_xxxxxy_0[k] * ab_x + g_0_z_xxxxxy_x[k];
            }

            /// Set up 74-75 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxz_0 = cbuffer.data(ks_geom_01_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxxz_0, g_0_z_xxxxxz_0, g_0_z_xxxxxz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxz_0[k] = -g_0_z_xxxxxz_0[k] * ab_x + g_0_z_xxxxxz_x[k];
            }

            /// Set up 75-76 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxyy_0 = cbuffer.data(ks_geom_01_off + 75 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxyy_0, g_0_z_xxxxyy_0, g_0_z_xxxxyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxyy_0[k] = -g_0_z_xxxxyy_0[k] * ab_x + g_0_z_xxxxyy_x[k];
            }

            /// Set up 76-77 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxyz_0 = cbuffer.data(ks_geom_01_off + 76 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxyz_0, g_0_z_xxxxyz_0, g_0_z_xxxxyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxyz_0[k] = -g_0_z_xxxxyz_0[k] * ab_x + g_0_z_xxxxyz_x[k];
            }

            /// Set up 77-78 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxzz_0 = cbuffer.data(ks_geom_01_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxzz_0, g_0_z_xxxxzz_0, g_0_z_xxxxzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxzz_0[k] = -g_0_z_xxxxzz_0[k] * ab_x + g_0_z_xxxxzz_x[k];
            }

            /// Set up 78-79 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyyy_0 = cbuffer.data(ks_geom_01_off + 78 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyyy_0, g_0_z_xxxyyy_0, g_0_z_xxxyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyyy_0[k] = -g_0_z_xxxyyy_0[k] * ab_x + g_0_z_xxxyyy_x[k];
            }

            /// Set up 79-80 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyyz_0 = cbuffer.data(ks_geom_01_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyyz_0, g_0_z_xxxyyz_0, g_0_z_xxxyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyyz_0[k] = -g_0_z_xxxyyz_0[k] * ab_x + g_0_z_xxxyyz_x[k];
            }

            /// Set up 80-81 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyzz_0 = cbuffer.data(ks_geom_01_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyzz_0, g_0_z_xxxyzz_0, g_0_z_xxxyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyzz_0[k] = -g_0_z_xxxyzz_0[k] * ab_x + g_0_z_xxxyzz_x[k];
            }

            /// Set up 81-82 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxzzz_0 = cbuffer.data(ks_geom_01_off + 81 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxzzz_0, g_0_z_xxxzzz_0, g_0_z_xxxzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxzzz_0[k] = -g_0_z_xxxzzz_0[k] * ab_x + g_0_z_xxxzzz_x[k];
            }

            /// Set up 82-83 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyyy_0 = cbuffer.data(ks_geom_01_off + 82 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyyy_0, g_0_z_xxyyyy_0, g_0_z_xxyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyyy_0[k] = -g_0_z_xxyyyy_0[k] * ab_x + g_0_z_xxyyyy_x[k];
            }

            /// Set up 83-84 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyyz_0 = cbuffer.data(ks_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyyz_0, g_0_z_xxyyyz_0, g_0_z_xxyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyyz_0[k] = -g_0_z_xxyyyz_0[k] * ab_x + g_0_z_xxyyyz_x[k];
            }

            /// Set up 84-85 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyzz_0 = cbuffer.data(ks_geom_01_off + 84 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyzz_0, g_0_z_xxyyzz_0, g_0_z_xxyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyzz_0[k] = -g_0_z_xxyyzz_0[k] * ab_x + g_0_z_xxyyzz_x[k];
            }

            /// Set up 85-86 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyzzz_0 = cbuffer.data(ks_geom_01_off + 85 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyzzz_0, g_0_z_xxyzzz_0, g_0_z_xxyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyzzz_0[k] = -g_0_z_xxyzzz_0[k] * ab_x + g_0_z_xxyzzz_x[k];
            }

            /// Set up 86-87 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzzzz_0 = cbuffer.data(ks_geom_01_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzzzz_0, g_0_z_xxzzzz_0, g_0_z_xxzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzzzz_0[k] = -g_0_z_xxzzzz_0[k] * ab_x + g_0_z_xxzzzz_x[k];
            }

            /// Set up 87-88 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyyy_0 = cbuffer.data(ks_geom_01_off + 87 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyyy_0, g_0_z_xyyyyy_0, g_0_z_xyyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyyy_0[k] = -g_0_z_xyyyyy_0[k] * ab_x + g_0_z_xyyyyy_x[k];
            }

            /// Set up 88-89 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyyz_0 = cbuffer.data(ks_geom_01_off + 88 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyyz_0, g_0_z_xyyyyz_0, g_0_z_xyyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyyz_0[k] = -g_0_z_xyyyyz_0[k] * ab_x + g_0_z_xyyyyz_x[k];
            }

            /// Set up 89-90 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyzz_0 = cbuffer.data(ks_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyzz_0, g_0_z_xyyyzz_0, g_0_z_xyyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyzz_0[k] = -g_0_z_xyyyzz_0[k] * ab_x + g_0_z_xyyyzz_x[k];
            }

            /// Set up 90-91 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyzzz_0 = cbuffer.data(ks_geom_01_off + 90 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyzzz_0, g_0_z_xyyzzz_0, g_0_z_xyyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyzzz_0[k] = -g_0_z_xyyzzz_0[k] * ab_x + g_0_z_xyyzzz_x[k];
            }

            /// Set up 91-92 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzzzz_0 = cbuffer.data(ks_geom_01_off + 91 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzzzz_0, g_0_z_xyzzzz_0, g_0_z_xyzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzzzz_0[k] = -g_0_z_xyzzzz_0[k] * ab_x + g_0_z_xyzzzz_x[k];
            }

            /// Set up 92-93 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzzzz_0 = cbuffer.data(ks_geom_01_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzzzz_0, g_0_z_xzzzzz_0, g_0_z_xzzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzzzz_0[k] = -g_0_z_xzzzzz_0[k] * ab_x + g_0_z_xzzzzz_x[k];
            }

            /// Set up 93-94 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyyy_0 = cbuffer.data(ks_geom_01_off + 93 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyyy_0, g_0_z_yyyyyy_0, g_0_z_yyyyyy_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyyy_0[k] = -g_0_z_yyyyyy_0[k] * ab_x + g_0_z_yyyyyy_x[k];
            }

            /// Set up 94-95 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyyz_0 = cbuffer.data(ks_geom_01_off + 94 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyyz_0, g_0_z_yyyyyz_0, g_0_z_yyyyyz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyyz_0[k] = -g_0_z_yyyyyz_0[k] * ab_x + g_0_z_yyyyyz_x[k];
            }

            /// Set up 95-96 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyzz_0 = cbuffer.data(ks_geom_01_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyzz_0, g_0_z_yyyyzz_0, g_0_z_yyyyzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyzz_0[k] = -g_0_z_yyyyzz_0[k] * ab_x + g_0_z_yyyyzz_x[k];
            }

            /// Set up 96-97 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyzzz_0 = cbuffer.data(ks_geom_01_off + 96 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyzzz_0, g_0_z_yyyzzz_0, g_0_z_yyyzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyzzz_0[k] = -g_0_z_yyyzzz_0[k] * ab_x + g_0_z_yyyzzz_x[k];
            }

            /// Set up 97-98 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzzzz_0 = cbuffer.data(ks_geom_01_off + 97 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzzzz_0, g_0_z_yyzzzz_0, g_0_z_yyzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzzzz_0[k] = -g_0_z_yyzzzz_0[k] * ab_x + g_0_z_yyzzzz_x[k];
            }

            /// Set up 98-99 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzzzz_0 = cbuffer.data(ks_geom_01_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzzzz_0, g_0_z_yzzzzz_0, g_0_z_yzzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzzzz_0[k] = -g_0_z_yzzzzz_0[k] * ab_x + g_0_z_yzzzzz_x[k];
            }

            /// Set up 99-100 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzzzz_0 = cbuffer.data(ks_geom_01_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzzzz_0, g_0_z_zzzzzz_0, g_0_z_zzzzzz_x  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzzzz_0[k] = -g_0_z_zzzzzz_0[k] * ab_x + g_0_z_zzzzzz_x[k];
            }

            /// Set up 100-101 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyyy_0 = cbuffer.data(ks_geom_01_off + 100 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyy_0, g_0_z_yyyyyy_y, g_0_z_yyyyyyy_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyyy_0[k] = -g_0_z_yyyyyy_0[k] * ab_y + g_0_z_yyyyyy_y[k];
            }

            /// Set up 101-102 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyyz_0 = cbuffer.data(ks_geom_01_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyyz_0, g_0_z_yyyyyz_0, g_0_z_yyyyyz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyyz_0[k] = -g_0_z_yyyyyz_0[k] * ab_y + g_0_z_yyyyyz_y[k];
            }

            /// Set up 102-103 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyzz_0 = cbuffer.data(ks_geom_01_off + 102 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyzz_0, g_0_z_yyyyzz_0, g_0_z_yyyyzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyzz_0[k] = -g_0_z_yyyyzz_0[k] * ab_y + g_0_z_yyyyzz_y[k];
            }

            /// Set up 103-104 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyzzz_0 = cbuffer.data(ks_geom_01_off + 103 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyzzz_0, g_0_z_yyyzzz_0, g_0_z_yyyzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyzzz_0[k] = -g_0_z_yyyzzz_0[k] * ab_y + g_0_z_yyyzzz_y[k];
            }

            /// Set up 104-105 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzzzz_0 = cbuffer.data(ks_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzzzz_0, g_0_z_yyzzzz_0, g_0_z_yyzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzzzz_0[k] = -g_0_z_yyzzzz_0[k] * ab_y + g_0_z_yyzzzz_y[k];
            }

            /// Set up 105-106 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzzzz_0 = cbuffer.data(ks_geom_01_off + 105 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzzzz_0, g_0_z_yzzzzz_0, g_0_z_yzzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzzzz_0[k] = -g_0_z_yzzzzz_0[k] * ab_y + g_0_z_yzzzzz_y[k];
            }

            /// Set up 106-107 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzzzz_0 = cbuffer.data(ks_geom_01_off + 106 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzzzz_0, g_0_z_zzzzzz_0, g_0_z_zzzzzz_y  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzzzz_0[k] = -g_0_z_zzzzzz_0[k] * ab_y + g_0_z_zzzzzz_y[k];
            }

            /// Set up 107-108 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzzzz_0 = cbuffer.data(ks_geom_01_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzzzz_0, g_0_z_zzzzzz_z, g_0_z_zzzzzzz_0, g_zzzzzz_0  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzzzz_0[k] = g_zzzzzz_0[k] - g_0_z_zzzzzz_0[k] * ab_z + g_0_z_zzzzzz_z[k];
            }
        }
    }
}

} // erirec namespace

