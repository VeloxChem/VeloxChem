#include "ElectronRepulsionGeom1000ContrRecFIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_fixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_fixx,
                                            const size_t idx_dixx,
                                            const size_t idx_geom_10_dixx,
                                            const size_t idx_geom_10_dkxx,
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
            /// Set up components of auxilary buffer : DISS

            const auto di_off = idx_dixx + i * dcomps + j;

            auto g_xx_xxxxxx = cbuffer.data(di_off + 0 * ccomps * dcomps);

            auto g_xx_xxxxxy = cbuffer.data(di_off + 1 * ccomps * dcomps);

            auto g_xx_xxxxxz = cbuffer.data(di_off + 2 * ccomps * dcomps);

            auto g_xx_xxxxyy = cbuffer.data(di_off + 3 * ccomps * dcomps);

            auto g_xx_xxxxyz = cbuffer.data(di_off + 4 * ccomps * dcomps);

            auto g_xx_xxxxzz = cbuffer.data(di_off + 5 * ccomps * dcomps);

            auto g_xx_xxxyyy = cbuffer.data(di_off + 6 * ccomps * dcomps);

            auto g_xx_xxxyyz = cbuffer.data(di_off + 7 * ccomps * dcomps);

            auto g_xx_xxxyzz = cbuffer.data(di_off + 8 * ccomps * dcomps);

            auto g_xx_xxxzzz = cbuffer.data(di_off + 9 * ccomps * dcomps);

            auto g_xx_xxyyyy = cbuffer.data(di_off + 10 * ccomps * dcomps);

            auto g_xx_xxyyyz = cbuffer.data(di_off + 11 * ccomps * dcomps);

            auto g_xx_xxyyzz = cbuffer.data(di_off + 12 * ccomps * dcomps);

            auto g_xx_xxyzzz = cbuffer.data(di_off + 13 * ccomps * dcomps);

            auto g_xx_xxzzzz = cbuffer.data(di_off + 14 * ccomps * dcomps);

            auto g_xx_xyyyyy = cbuffer.data(di_off + 15 * ccomps * dcomps);

            auto g_xx_xyyyyz = cbuffer.data(di_off + 16 * ccomps * dcomps);

            auto g_xx_xyyyzz = cbuffer.data(di_off + 17 * ccomps * dcomps);

            auto g_xx_xyyzzz = cbuffer.data(di_off + 18 * ccomps * dcomps);

            auto g_xx_xyzzzz = cbuffer.data(di_off + 19 * ccomps * dcomps);

            auto g_xx_xzzzzz = cbuffer.data(di_off + 20 * ccomps * dcomps);

            auto g_xx_yyyyyy = cbuffer.data(di_off + 21 * ccomps * dcomps);

            auto g_xx_yyyyyz = cbuffer.data(di_off + 22 * ccomps * dcomps);

            auto g_xx_yyyyzz = cbuffer.data(di_off + 23 * ccomps * dcomps);

            auto g_xx_yyyzzz = cbuffer.data(di_off + 24 * ccomps * dcomps);

            auto g_xx_yyzzzz = cbuffer.data(di_off + 25 * ccomps * dcomps);

            auto g_xx_yzzzzz = cbuffer.data(di_off + 26 * ccomps * dcomps);

            auto g_xx_zzzzzz = cbuffer.data(di_off + 27 * ccomps * dcomps);

            auto g_yy_xxxxxx = cbuffer.data(di_off + 84 * ccomps * dcomps);

            auto g_yy_xxxxxy = cbuffer.data(di_off + 85 * ccomps * dcomps);

            auto g_yy_xxxxxz = cbuffer.data(di_off + 86 * ccomps * dcomps);

            auto g_yy_xxxxyy = cbuffer.data(di_off + 87 * ccomps * dcomps);

            auto g_yy_xxxxyz = cbuffer.data(di_off + 88 * ccomps * dcomps);

            auto g_yy_xxxxzz = cbuffer.data(di_off + 89 * ccomps * dcomps);

            auto g_yy_xxxyyy = cbuffer.data(di_off + 90 * ccomps * dcomps);

            auto g_yy_xxxyyz = cbuffer.data(di_off + 91 * ccomps * dcomps);

            auto g_yy_xxxyzz = cbuffer.data(di_off + 92 * ccomps * dcomps);

            auto g_yy_xxxzzz = cbuffer.data(di_off + 93 * ccomps * dcomps);

            auto g_yy_xxyyyy = cbuffer.data(di_off + 94 * ccomps * dcomps);

            auto g_yy_xxyyyz = cbuffer.data(di_off + 95 * ccomps * dcomps);

            auto g_yy_xxyyzz = cbuffer.data(di_off + 96 * ccomps * dcomps);

            auto g_yy_xxyzzz = cbuffer.data(di_off + 97 * ccomps * dcomps);

            auto g_yy_xxzzzz = cbuffer.data(di_off + 98 * ccomps * dcomps);

            auto g_yy_xyyyyy = cbuffer.data(di_off + 99 * ccomps * dcomps);

            auto g_yy_xyyyyz = cbuffer.data(di_off + 100 * ccomps * dcomps);

            auto g_yy_xyyyzz = cbuffer.data(di_off + 101 * ccomps * dcomps);

            auto g_yy_xyyzzz = cbuffer.data(di_off + 102 * ccomps * dcomps);

            auto g_yy_xyzzzz = cbuffer.data(di_off + 103 * ccomps * dcomps);

            auto g_yy_xzzzzz = cbuffer.data(di_off + 104 * ccomps * dcomps);

            auto g_yy_yyyyyy = cbuffer.data(di_off + 105 * ccomps * dcomps);

            auto g_yy_yyyyyz = cbuffer.data(di_off + 106 * ccomps * dcomps);

            auto g_yy_yyyyzz = cbuffer.data(di_off + 107 * ccomps * dcomps);

            auto g_yy_yyyzzz = cbuffer.data(di_off + 108 * ccomps * dcomps);

            auto g_yy_yyzzzz = cbuffer.data(di_off + 109 * ccomps * dcomps);

            auto g_yy_yzzzzz = cbuffer.data(di_off + 110 * ccomps * dcomps);

            auto g_yy_zzzzzz = cbuffer.data(di_off + 111 * ccomps * dcomps);

            auto g_zz_xxxxxx = cbuffer.data(di_off + 140 * ccomps * dcomps);

            auto g_zz_xxxxxy = cbuffer.data(di_off + 141 * ccomps * dcomps);

            auto g_zz_xxxxxz = cbuffer.data(di_off + 142 * ccomps * dcomps);

            auto g_zz_xxxxyy = cbuffer.data(di_off + 143 * ccomps * dcomps);

            auto g_zz_xxxxyz = cbuffer.data(di_off + 144 * ccomps * dcomps);

            auto g_zz_xxxxzz = cbuffer.data(di_off + 145 * ccomps * dcomps);

            auto g_zz_xxxyyy = cbuffer.data(di_off + 146 * ccomps * dcomps);

            auto g_zz_xxxyyz = cbuffer.data(di_off + 147 * ccomps * dcomps);

            auto g_zz_xxxyzz = cbuffer.data(di_off + 148 * ccomps * dcomps);

            auto g_zz_xxxzzz = cbuffer.data(di_off + 149 * ccomps * dcomps);

            auto g_zz_xxyyyy = cbuffer.data(di_off + 150 * ccomps * dcomps);

            auto g_zz_xxyyyz = cbuffer.data(di_off + 151 * ccomps * dcomps);

            auto g_zz_xxyyzz = cbuffer.data(di_off + 152 * ccomps * dcomps);

            auto g_zz_xxyzzz = cbuffer.data(di_off + 153 * ccomps * dcomps);

            auto g_zz_xxzzzz = cbuffer.data(di_off + 154 * ccomps * dcomps);

            auto g_zz_xyyyyy = cbuffer.data(di_off + 155 * ccomps * dcomps);

            auto g_zz_xyyyyz = cbuffer.data(di_off + 156 * ccomps * dcomps);

            auto g_zz_xyyyzz = cbuffer.data(di_off + 157 * ccomps * dcomps);

            auto g_zz_xyyzzz = cbuffer.data(di_off + 158 * ccomps * dcomps);

            auto g_zz_xyzzzz = cbuffer.data(di_off + 159 * ccomps * dcomps);

            auto g_zz_xzzzzz = cbuffer.data(di_off + 160 * ccomps * dcomps);

            auto g_zz_yyyyyy = cbuffer.data(di_off + 161 * ccomps * dcomps);

            auto g_zz_yyyyyz = cbuffer.data(di_off + 162 * ccomps * dcomps);

            auto g_zz_yyyyzz = cbuffer.data(di_off + 163 * ccomps * dcomps);

            auto g_zz_yyyzzz = cbuffer.data(di_off + 164 * ccomps * dcomps);

            auto g_zz_yyzzzz = cbuffer.data(di_off + 165 * ccomps * dcomps);

            auto g_zz_yzzzzz = cbuffer.data(di_off + 166 * ccomps * dcomps);

            auto g_zz_zzzzzz = cbuffer.data(di_off + 167 * ccomps * dcomps);

            /// Set up components of auxilary buffer : DISS

            const auto di_geom_10_off = idx_geom_10_dixx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : DKSS

            const auto dk_geom_10_off = idx_geom_10_dkxx + i * dcomps + j;

            auto g_x_0_xx_xxxxxxx = cbuffer.data(dk_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xx_xxxxxxy = cbuffer.data(dk_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xx_xxxxxxz = cbuffer.data(dk_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xx_xxxxxyy = cbuffer.data(dk_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xx_xxxxxyz = cbuffer.data(dk_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xx_xxxxxzz = cbuffer.data(dk_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xx_xxxxyyy = cbuffer.data(dk_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xx_xxxxyyz = cbuffer.data(dk_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xx_xxxxyzz = cbuffer.data(dk_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xx_xxxxzzz = cbuffer.data(dk_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xx_xxxyyyy = cbuffer.data(dk_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xx_xxxyyyz = cbuffer.data(dk_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xx_xxxyyzz = cbuffer.data(dk_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xx_xxxyzzz = cbuffer.data(dk_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xx_xxxzzzz = cbuffer.data(dk_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xx_xxyyyyy = cbuffer.data(dk_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xx_xxyyyyz = cbuffer.data(dk_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xx_xxyyyzz = cbuffer.data(dk_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xx_xxyyzzz = cbuffer.data(dk_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xx_xxyzzzz = cbuffer.data(dk_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xx_xxzzzzz = cbuffer.data(dk_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xx_xyyyyyy = cbuffer.data(dk_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xx_xyyyyyz = cbuffer.data(dk_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xx_xyyyyzz = cbuffer.data(dk_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xx_xyyyzzz = cbuffer.data(dk_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xx_xyyzzzz = cbuffer.data(dk_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xx_xyzzzzz = cbuffer.data(dk_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xx_xzzzzzz = cbuffer.data(dk_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xx_yyyyyyy = cbuffer.data(dk_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xx_yyyyyyz = cbuffer.data(dk_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xx_yyyyyzz = cbuffer.data(dk_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xx_yyyyzzz = cbuffer.data(dk_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xx_yyyzzzz = cbuffer.data(dk_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xx_yyzzzzz = cbuffer.data(dk_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xx_yzzzzzz = cbuffer.data(dk_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xx_zzzzzzz = cbuffer.data(dk_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xy_xxxxxxy = cbuffer.data(dk_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xy_xxxxxyy = cbuffer.data(dk_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xy_xxxxxyz = cbuffer.data(dk_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xy_xxxxyyy = cbuffer.data(dk_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xy_xxxxyyz = cbuffer.data(dk_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xy_xxxxyzz = cbuffer.data(dk_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xy_xxxyyyy = cbuffer.data(dk_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xy_xxxyyyz = cbuffer.data(dk_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xy_xxxyyzz = cbuffer.data(dk_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xy_xxxyzzz = cbuffer.data(dk_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xy_xxyyyyy = cbuffer.data(dk_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xy_xxyyyyz = cbuffer.data(dk_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xy_xxyyyzz = cbuffer.data(dk_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xy_xxyyzzz = cbuffer.data(dk_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xy_xxyzzzz = cbuffer.data(dk_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xy_xyyyyyy = cbuffer.data(dk_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xy_xyyyyyz = cbuffer.data(dk_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xy_xyyyyzz = cbuffer.data(dk_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xy_xyyyzzz = cbuffer.data(dk_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xy_xyyzzzz = cbuffer.data(dk_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xy_xyzzzzz = cbuffer.data(dk_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xy_yyyyyyy = cbuffer.data(dk_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xy_yyyyyyz = cbuffer.data(dk_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xy_yyyyyzz = cbuffer.data(dk_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xy_yyyyzzz = cbuffer.data(dk_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xy_yyyzzzz = cbuffer.data(dk_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xy_yyzzzzz = cbuffer.data(dk_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xy_yzzzzzz = cbuffer.data(dk_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xz_xxxxxxy = cbuffer.data(dk_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xz_xxxxxxz = cbuffer.data(dk_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xz_xxxxxyy = cbuffer.data(dk_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xz_xxxxxyz = cbuffer.data(dk_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xz_xxxxxzz = cbuffer.data(dk_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xz_xxxxyyy = cbuffer.data(dk_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xz_xxxxyyz = cbuffer.data(dk_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xz_xxxxyzz = cbuffer.data(dk_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xz_xxxxzzz = cbuffer.data(dk_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xz_xxxyyyy = cbuffer.data(dk_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xz_xxxyyyz = cbuffer.data(dk_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xz_xxxyyzz = cbuffer.data(dk_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xz_xxxyzzz = cbuffer.data(dk_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xz_xxxzzzz = cbuffer.data(dk_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xz_xxyyyyy = cbuffer.data(dk_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xz_xxyyyyz = cbuffer.data(dk_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xz_xxyyyzz = cbuffer.data(dk_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xz_xxyyzzz = cbuffer.data(dk_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xz_xxyzzzz = cbuffer.data(dk_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xz_xxzzzzz = cbuffer.data(dk_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xz_xyyyyyy = cbuffer.data(dk_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xz_xyyyyyz = cbuffer.data(dk_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xz_xyyyyzz = cbuffer.data(dk_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xz_xyyyzzz = cbuffer.data(dk_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xz_xyyzzzz = cbuffer.data(dk_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xz_xyzzzzz = cbuffer.data(dk_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xz_xzzzzzz = cbuffer.data(dk_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xz_yyyyyyy = cbuffer.data(dk_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xz_yyyyyyz = cbuffer.data(dk_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xz_yyyyyzz = cbuffer.data(dk_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xz_yyyyzzz = cbuffer.data(dk_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xz_yyyzzzz = cbuffer.data(dk_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xz_yyzzzzz = cbuffer.data(dk_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xz_yzzzzzz = cbuffer.data(dk_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xz_zzzzzzz = cbuffer.data(dk_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_yy_xxxxxxy = cbuffer.data(dk_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_yy_xxxxxyy = cbuffer.data(dk_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_yy_xxxxxyz = cbuffer.data(dk_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_yy_xxxxyyy = cbuffer.data(dk_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_yy_xxxxyyz = cbuffer.data(dk_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_yy_xxxxyzz = cbuffer.data(dk_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_yy_xxxyyyy = cbuffer.data(dk_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_yy_xxxyyyz = cbuffer.data(dk_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_yy_xxxyyzz = cbuffer.data(dk_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_yy_xxxyzzz = cbuffer.data(dk_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_yy_xxyyyyy = cbuffer.data(dk_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_yy_xxyyyyz = cbuffer.data(dk_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_yy_xxyyyzz = cbuffer.data(dk_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_yy_xxyyzzz = cbuffer.data(dk_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_yy_xxyzzzz = cbuffer.data(dk_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_yy_xyyyyyy = cbuffer.data(dk_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_yy_xyyyyyz = cbuffer.data(dk_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_yy_xyyyyzz = cbuffer.data(dk_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_yy_xyyyzzz = cbuffer.data(dk_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_yy_xyyzzzz = cbuffer.data(dk_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_yy_xyzzzzz = cbuffer.data(dk_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_yy_yyyyyyy = cbuffer.data(dk_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_yy_yyyyyyz = cbuffer.data(dk_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_yy_yyyyyzz = cbuffer.data(dk_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_yy_yyyyzzz = cbuffer.data(dk_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_yy_yyyzzzz = cbuffer.data(dk_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_yy_yyzzzzz = cbuffer.data(dk_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_yy_yzzzzzz = cbuffer.data(dk_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_yz_xxxxxxy = cbuffer.data(dk_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_yz_xxxxxyy = cbuffer.data(dk_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_yz_xxxxxyz = cbuffer.data(dk_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_yz_xxxxyyy = cbuffer.data(dk_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_yz_xxxxyyz = cbuffer.data(dk_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_yz_xxxxyzz = cbuffer.data(dk_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_yz_xxxyyyy = cbuffer.data(dk_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_yz_xxxyyyz = cbuffer.data(dk_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_yz_xxxyyzz = cbuffer.data(dk_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_yz_xxxyzzz = cbuffer.data(dk_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_yz_xxyyyyy = cbuffer.data(dk_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_yz_xxyyyyz = cbuffer.data(dk_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_yz_xxyyyzz = cbuffer.data(dk_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_yz_xxyyzzz = cbuffer.data(dk_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_yz_xxyzzzz = cbuffer.data(dk_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_yz_xyyyyyy = cbuffer.data(dk_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_yz_xyyyyyz = cbuffer.data(dk_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_yz_xyyyyzz = cbuffer.data(dk_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_yz_xyyyzzz = cbuffer.data(dk_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_yz_xyyzzzz = cbuffer.data(dk_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_yz_xyzzzzz = cbuffer.data(dk_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_yz_yyyyyyy = cbuffer.data(dk_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_yz_yyyyyyz = cbuffer.data(dk_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_yz_yyyyyzz = cbuffer.data(dk_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_yz_yyyyzzz = cbuffer.data(dk_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_yz_yyyzzzz = cbuffer.data(dk_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_yz_yyzzzzz = cbuffer.data(dk_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_yz_yzzzzzz = cbuffer.data(dk_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_zz_xxxxxxy = cbuffer.data(dk_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_zz_xxxxxxz = cbuffer.data(dk_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_zz_xxxxxyy = cbuffer.data(dk_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_zz_xxxxxyz = cbuffer.data(dk_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_zz_xxxxxzz = cbuffer.data(dk_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_zz_xxxxyyy = cbuffer.data(dk_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_zz_xxxxyyz = cbuffer.data(dk_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_zz_xxxxyzz = cbuffer.data(dk_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_zz_xxxxzzz = cbuffer.data(dk_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_zz_xxxyyyy = cbuffer.data(dk_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_zz_xxxyyyz = cbuffer.data(dk_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_zz_xxxyyzz = cbuffer.data(dk_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_zz_xxxyzzz = cbuffer.data(dk_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_zz_xxxzzzz = cbuffer.data(dk_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_zz_xxyyyyy = cbuffer.data(dk_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_zz_xxyyyyz = cbuffer.data(dk_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_zz_xxyyyzz = cbuffer.data(dk_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_zz_xxyyzzz = cbuffer.data(dk_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_zz_xxyzzzz = cbuffer.data(dk_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_zz_xxzzzzz = cbuffer.data(dk_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_zz_xyyyyyy = cbuffer.data(dk_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_zz_xyyyyyz = cbuffer.data(dk_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_zz_xyyyyzz = cbuffer.data(dk_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_zz_xyyyzzz = cbuffer.data(dk_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_zz_xyyzzzz = cbuffer.data(dk_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_zz_xyzzzzz = cbuffer.data(dk_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_zz_xzzzzzz = cbuffer.data(dk_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_zz_yyyyyyy = cbuffer.data(dk_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_zz_yyyyyyz = cbuffer.data(dk_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_zz_yyyyyzz = cbuffer.data(dk_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_zz_yyyyzzz = cbuffer.data(dk_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_zz_yyyzzzz = cbuffer.data(dk_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_zz_yyzzzzz = cbuffer.data(dk_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_zz_yzzzzzz = cbuffer.data(dk_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_zz_zzzzzzz = cbuffer.data(dk_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxxx = cbuffer.data(dk_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxxy = cbuffer.data(dk_geom_10_off + 217 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxxz = cbuffer.data(dk_geom_10_off + 218 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxyy = cbuffer.data(dk_geom_10_off + 219 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxyz = cbuffer.data(dk_geom_10_off + 220 * ccomps * dcomps);

            auto g_y_0_xx_xxxxxzz = cbuffer.data(dk_geom_10_off + 221 * ccomps * dcomps);

            auto g_y_0_xx_xxxxyyy = cbuffer.data(dk_geom_10_off + 222 * ccomps * dcomps);

            auto g_y_0_xx_xxxxyyz = cbuffer.data(dk_geom_10_off + 223 * ccomps * dcomps);

            auto g_y_0_xx_xxxxyzz = cbuffer.data(dk_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_xx_xxxxzzz = cbuffer.data(dk_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xx_xxxyyyy = cbuffer.data(dk_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xx_xxxyyyz = cbuffer.data(dk_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xx_xxxyyzz = cbuffer.data(dk_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xx_xxxyzzz = cbuffer.data(dk_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_xx_xxxzzzz = cbuffer.data(dk_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_xx_xxyyyyy = cbuffer.data(dk_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xx_xxyyyyz = cbuffer.data(dk_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xx_xxyyyzz = cbuffer.data(dk_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xx_xxyyzzz = cbuffer.data(dk_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xx_xxyzzzz = cbuffer.data(dk_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xx_xxzzzzz = cbuffer.data(dk_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xx_xyyyyyy = cbuffer.data(dk_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xx_xyyyyyz = cbuffer.data(dk_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xx_xyyyyzz = cbuffer.data(dk_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_xx_xyyyzzz = cbuffer.data(dk_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_xx_xyyzzzz = cbuffer.data(dk_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_xx_xyzzzzz = cbuffer.data(dk_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_xx_xzzzzzz = cbuffer.data(dk_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxxx = cbuffer.data(dk_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxxy = cbuffer.data(dk_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxxz = cbuffer.data(dk_geom_10_off + 254 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxyy = cbuffer.data(dk_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxyz = cbuffer.data(dk_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_xy_xxxxxzz = cbuffer.data(dk_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_xy_xxxxyyy = cbuffer.data(dk_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_xy_xxxxyyz = cbuffer.data(dk_geom_10_off + 259 * ccomps * dcomps);

            auto g_y_0_xy_xxxxyzz = cbuffer.data(dk_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_xy_xxxxzzz = cbuffer.data(dk_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_xy_xxxyyyy = cbuffer.data(dk_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_xy_xxxyyyz = cbuffer.data(dk_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_xy_xxxyyzz = cbuffer.data(dk_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_xy_xxxyzzz = cbuffer.data(dk_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_xy_xxxzzzz = cbuffer.data(dk_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_xy_xxyyyyy = cbuffer.data(dk_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_xy_xxyyyyz = cbuffer.data(dk_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_xy_xxyyyzz = cbuffer.data(dk_geom_10_off + 269 * ccomps * dcomps);

            auto g_y_0_xy_xxyyzzz = cbuffer.data(dk_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_xy_xxyzzzz = cbuffer.data(dk_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_xy_xxzzzzz = cbuffer.data(dk_geom_10_off + 272 * ccomps * dcomps);

            auto g_y_0_xy_xyyyyyy = cbuffer.data(dk_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_xy_xyyyyyz = cbuffer.data(dk_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_xy_xyyyyzz = cbuffer.data(dk_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_xy_xyyyzzz = cbuffer.data(dk_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_xy_xyyzzzz = cbuffer.data(dk_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_xy_xyzzzzz = cbuffer.data(dk_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_xy_xzzzzzz = cbuffer.data(dk_geom_10_off + 279 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxxx = cbuffer.data(dk_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxxy = cbuffer.data(dk_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxxz = cbuffer.data(dk_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxyy = cbuffer.data(dk_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxyz = cbuffer.data(dk_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_xz_xxxxxzz = cbuffer.data(dk_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_xz_xxxxyyy = cbuffer.data(dk_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_xz_xxxxyyz = cbuffer.data(dk_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_xz_xxxxyzz = cbuffer.data(dk_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_xz_xxxxzzz = cbuffer.data(dk_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_xz_xxxyyyy = cbuffer.data(dk_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_xz_xxxyyyz = cbuffer.data(dk_geom_10_off + 299 * ccomps * dcomps);

            auto g_y_0_xz_xxxyyzz = cbuffer.data(dk_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_xz_xxxyzzz = cbuffer.data(dk_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_xz_xxxzzzz = cbuffer.data(dk_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_xz_xxyyyyy = cbuffer.data(dk_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_xz_xxyyyyz = cbuffer.data(dk_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_xz_xxyyyzz = cbuffer.data(dk_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_xz_xxyyzzz = cbuffer.data(dk_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_xz_xxyzzzz = cbuffer.data(dk_geom_10_off + 307 * ccomps * dcomps);

            auto g_y_0_xz_xxzzzzz = cbuffer.data(dk_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_xz_xyyyyyy = cbuffer.data(dk_geom_10_off + 309 * ccomps * dcomps);

            auto g_y_0_xz_xyyyyyz = cbuffer.data(dk_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_xz_xyyyyzz = cbuffer.data(dk_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_xz_xyyyzzz = cbuffer.data(dk_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_xz_xyyzzzz = cbuffer.data(dk_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_xz_xyzzzzz = cbuffer.data(dk_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_xz_xzzzzzz = cbuffer.data(dk_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxxx = cbuffer.data(dk_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxxy = cbuffer.data(dk_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxxz = cbuffer.data(dk_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxyy = cbuffer.data(dk_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxyz = cbuffer.data(dk_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_yy_xxxxxzz = cbuffer.data(dk_geom_10_off + 329 * ccomps * dcomps);

            auto g_y_0_yy_xxxxyyy = cbuffer.data(dk_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_yy_xxxxyyz = cbuffer.data(dk_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_yy_xxxxyzz = cbuffer.data(dk_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_yy_xxxxzzz = cbuffer.data(dk_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_yy_xxxyyyy = cbuffer.data(dk_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_yy_xxxyyyz = cbuffer.data(dk_geom_10_off + 335 * ccomps * dcomps);

            auto g_y_0_yy_xxxyyzz = cbuffer.data(dk_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_yy_xxxyzzz = cbuffer.data(dk_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_yy_xxxzzzz = cbuffer.data(dk_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_yy_xxyyyyy = cbuffer.data(dk_geom_10_off + 339 * ccomps * dcomps);

            auto g_y_0_yy_xxyyyyz = cbuffer.data(dk_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_yy_xxyyyzz = cbuffer.data(dk_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_yy_xxyyzzz = cbuffer.data(dk_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_yy_xxyzzzz = cbuffer.data(dk_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_yy_xxzzzzz = cbuffer.data(dk_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_yy_xyyyyyy = cbuffer.data(dk_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_yy_xyyyyyz = cbuffer.data(dk_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_yy_xyyyyzz = cbuffer.data(dk_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_yy_xyyyzzz = cbuffer.data(dk_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_yy_xyyzzzz = cbuffer.data(dk_geom_10_off + 349 * ccomps * dcomps);

            auto g_y_0_yy_xyzzzzz = cbuffer.data(dk_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_yy_xzzzzzz = cbuffer.data(dk_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_yy_yyyyyyy = cbuffer.data(dk_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_yy_yyyyyyz = cbuffer.data(dk_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_yy_yyyyyzz = cbuffer.data(dk_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_yy_yyyyzzz = cbuffer.data(dk_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_yy_yyyzzzz = cbuffer.data(dk_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_yy_yyzzzzz = cbuffer.data(dk_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_yy_yzzzzzz = cbuffer.data(dk_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_yy_zzzzzzz = cbuffer.data(dk_geom_10_off + 359 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxxx = cbuffer.data(dk_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxxy = cbuffer.data(dk_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxxz = cbuffer.data(dk_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxyy = cbuffer.data(dk_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxyz = cbuffer.data(dk_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_yz_xxxxxzz = cbuffer.data(dk_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_yz_xxxxyyy = cbuffer.data(dk_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_yz_xxxxyyz = cbuffer.data(dk_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_yz_xxxxyzz = cbuffer.data(dk_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_yz_xxxxzzz = cbuffer.data(dk_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_yz_xxxyyyy = cbuffer.data(dk_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_yz_xxxyyyz = cbuffer.data(dk_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_yz_xxxyyzz = cbuffer.data(dk_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_yz_xxxyzzz = cbuffer.data(dk_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_yz_xxxzzzz = cbuffer.data(dk_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_yz_xxyyyyy = cbuffer.data(dk_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_yz_xxyyyyz = cbuffer.data(dk_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_yz_xxyyyzz = cbuffer.data(dk_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_yz_xxyyzzz = cbuffer.data(dk_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_yz_xxyzzzz = cbuffer.data(dk_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_yz_xxzzzzz = cbuffer.data(dk_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_yz_xyyyyyy = cbuffer.data(dk_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_yz_xyyyyyz = cbuffer.data(dk_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_yz_xyyyyzz = cbuffer.data(dk_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_yz_xyyyzzz = cbuffer.data(dk_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_yz_xyyzzzz = cbuffer.data(dk_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_yz_xyzzzzz = cbuffer.data(dk_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_yz_xzzzzzz = cbuffer.data(dk_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_yz_yyyyyyz = cbuffer.data(dk_geom_10_off + 389 * ccomps * dcomps);

            auto g_y_0_yz_yyyyyzz = cbuffer.data(dk_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_yz_yyyyzzz = cbuffer.data(dk_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_yz_yyyzzzz = cbuffer.data(dk_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_yz_yyzzzzz = cbuffer.data(dk_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_yz_yzzzzzz = cbuffer.data(dk_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_yz_zzzzzzz = cbuffer.data(dk_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxxx = cbuffer.data(dk_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxxy = cbuffer.data(dk_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxxz = cbuffer.data(dk_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxyy = cbuffer.data(dk_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxyz = cbuffer.data(dk_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_zz_xxxxxzz = cbuffer.data(dk_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_zz_xxxxyyy = cbuffer.data(dk_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_zz_xxxxyyz = cbuffer.data(dk_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_zz_xxxxyzz = cbuffer.data(dk_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_zz_xxxxzzz = cbuffer.data(dk_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_zz_xxxyyyy = cbuffer.data(dk_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_zz_xxxyyyz = cbuffer.data(dk_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_zz_xxxyyzz = cbuffer.data(dk_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_zz_xxxyzzz = cbuffer.data(dk_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_zz_xxxzzzz = cbuffer.data(dk_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_zz_xxyyyyy = cbuffer.data(dk_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_zz_xxyyyyz = cbuffer.data(dk_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_zz_xxyyyzz = cbuffer.data(dk_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_zz_xxyyzzz = cbuffer.data(dk_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_zz_xxyzzzz = cbuffer.data(dk_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_zz_xxzzzzz = cbuffer.data(dk_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_zz_xyyyyyy = cbuffer.data(dk_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_zz_xyyyyyz = cbuffer.data(dk_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_zz_xyyyyzz = cbuffer.data(dk_geom_10_off + 419 * ccomps * dcomps);

            auto g_y_0_zz_xyyyzzz = cbuffer.data(dk_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_zz_xyyzzzz = cbuffer.data(dk_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_zz_xyzzzzz = cbuffer.data(dk_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_zz_xzzzzzz = cbuffer.data(dk_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_zz_yyyyyyz = cbuffer.data(dk_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_zz_yyyyyzz = cbuffer.data(dk_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_zz_yyyyzzz = cbuffer.data(dk_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_zz_yyyzzzz = cbuffer.data(dk_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_zz_yyzzzzz = cbuffer.data(dk_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_zz_yzzzzzz = cbuffer.data(dk_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_zz_zzzzzzz = cbuffer.data(dk_geom_10_off + 431 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxxx = cbuffer.data(dk_geom_10_off + 432 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxxy = cbuffer.data(dk_geom_10_off + 433 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxxz = cbuffer.data(dk_geom_10_off + 434 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxyy = cbuffer.data(dk_geom_10_off + 435 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxyz = cbuffer.data(dk_geom_10_off + 436 * ccomps * dcomps);

            auto g_z_0_xx_xxxxxzz = cbuffer.data(dk_geom_10_off + 437 * ccomps * dcomps);

            auto g_z_0_xx_xxxxyyy = cbuffer.data(dk_geom_10_off + 438 * ccomps * dcomps);

            auto g_z_0_xx_xxxxyyz = cbuffer.data(dk_geom_10_off + 439 * ccomps * dcomps);

            auto g_z_0_xx_xxxxyzz = cbuffer.data(dk_geom_10_off + 440 * ccomps * dcomps);

            auto g_z_0_xx_xxxxzzz = cbuffer.data(dk_geom_10_off + 441 * ccomps * dcomps);

            auto g_z_0_xx_xxxyyyy = cbuffer.data(dk_geom_10_off + 442 * ccomps * dcomps);

            auto g_z_0_xx_xxxyyyz = cbuffer.data(dk_geom_10_off + 443 * ccomps * dcomps);

            auto g_z_0_xx_xxxyyzz = cbuffer.data(dk_geom_10_off + 444 * ccomps * dcomps);

            auto g_z_0_xx_xxxyzzz = cbuffer.data(dk_geom_10_off + 445 * ccomps * dcomps);

            auto g_z_0_xx_xxxzzzz = cbuffer.data(dk_geom_10_off + 446 * ccomps * dcomps);

            auto g_z_0_xx_xxyyyyy = cbuffer.data(dk_geom_10_off + 447 * ccomps * dcomps);

            auto g_z_0_xx_xxyyyyz = cbuffer.data(dk_geom_10_off + 448 * ccomps * dcomps);

            auto g_z_0_xx_xxyyyzz = cbuffer.data(dk_geom_10_off + 449 * ccomps * dcomps);

            auto g_z_0_xx_xxyyzzz = cbuffer.data(dk_geom_10_off + 450 * ccomps * dcomps);

            auto g_z_0_xx_xxyzzzz = cbuffer.data(dk_geom_10_off + 451 * ccomps * dcomps);

            auto g_z_0_xx_xxzzzzz = cbuffer.data(dk_geom_10_off + 452 * ccomps * dcomps);

            auto g_z_0_xx_xyyyyyy = cbuffer.data(dk_geom_10_off + 453 * ccomps * dcomps);

            auto g_z_0_xx_xyyyyyz = cbuffer.data(dk_geom_10_off + 454 * ccomps * dcomps);

            auto g_z_0_xx_xyyyyzz = cbuffer.data(dk_geom_10_off + 455 * ccomps * dcomps);

            auto g_z_0_xx_xyyyzzz = cbuffer.data(dk_geom_10_off + 456 * ccomps * dcomps);

            auto g_z_0_xx_xyyzzzz = cbuffer.data(dk_geom_10_off + 457 * ccomps * dcomps);

            auto g_z_0_xx_xyzzzzz = cbuffer.data(dk_geom_10_off + 458 * ccomps * dcomps);

            auto g_z_0_xx_xzzzzzz = cbuffer.data(dk_geom_10_off + 459 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxxx = cbuffer.data(dk_geom_10_off + 468 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxxy = cbuffer.data(dk_geom_10_off + 469 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxxz = cbuffer.data(dk_geom_10_off + 470 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxyy = cbuffer.data(dk_geom_10_off + 471 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxyz = cbuffer.data(dk_geom_10_off + 472 * ccomps * dcomps);

            auto g_z_0_xy_xxxxxzz = cbuffer.data(dk_geom_10_off + 473 * ccomps * dcomps);

            auto g_z_0_xy_xxxxyyy = cbuffer.data(dk_geom_10_off + 474 * ccomps * dcomps);

            auto g_z_0_xy_xxxxyyz = cbuffer.data(dk_geom_10_off + 475 * ccomps * dcomps);

            auto g_z_0_xy_xxxxyzz = cbuffer.data(dk_geom_10_off + 476 * ccomps * dcomps);

            auto g_z_0_xy_xxxxzzz = cbuffer.data(dk_geom_10_off + 477 * ccomps * dcomps);

            auto g_z_0_xy_xxxyyyy = cbuffer.data(dk_geom_10_off + 478 * ccomps * dcomps);

            auto g_z_0_xy_xxxyyyz = cbuffer.data(dk_geom_10_off + 479 * ccomps * dcomps);

            auto g_z_0_xy_xxxyyzz = cbuffer.data(dk_geom_10_off + 480 * ccomps * dcomps);

            auto g_z_0_xy_xxxyzzz = cbuffer.data(dk_geom_10_off + 481 * ccomps * dcomps);

            auto g_z_0_xy_xxxzzzz = cbuffer.data(dk_geom_10_off + 482 * ccomps * dcomps);

            auto g_z_0_xy_xxyyyyy = cbuffer.data(dk_geom_10_off + 483 * ccomps * dcomps);

            auto g_z_0_xy_xxyyyyz = cbuffer.data(dk_geom_10_off + 484 * ccomps * dcomps);

            auto g_z_0_xy_xxyyyzz = cbuffer.data(dk_geom_10_off + 485 * ccomps * dcomps);

            auto g_z_0_xy_xxyyzzz = cbuffer.data(dk_geom_10_off + 486 * ccomps * dcomps);

            auto g_z_0_xy_xxyzzzz = cbuffer.data(dk_geom_10_off + 487 * ccomps * dcomps);

            auto g_z_0_xy_xxzzzzz = cbuffer.data(dk_geom_10_off + 488 * ccomps * dcomps);

            auto g_z_0_xy_xyyyyyy = cbuffer.data(dk_geom_10_off + 489 * ccomps * dcomps);

            auto g_z_0_xy_xyyyyyz = cbuffer.data(dk_geom_10_off + 490 * ccomps * dcomps);

            auto g_z_0_xy_xyyyyzz = cbuffer.data(dk_geom_10_off + 491 * ccomps * dcomps);

            auto g_z_0_xy_xyyyzzz = cbuffer.data(dk_geom_10_off + 492 * ccomps * dcomps);

            auto g_z_0_xy_xyyzzzz = cbuffer.data(dk_geom_10_off + 493 * ccomps * dcomps);

            auto g_z_0_xy_xyzzzzz = cbuffer.data(dk_geom_10_off + 494 * ccomps * dcomps);

            auto g_z_0_xy_xzzzzzz = cbuffer.data(dk_geom_10_off + 495 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxxx = cbuffer.data(dk_geom_10_off + 504 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxxy = cbuffer.data(dk_geom_10_off + 505 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxxz = cbuffer.data(dk_geom_10_off + 506 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxyy = cbuffer.data(dk_geom_10_off + 507 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxyz = cbuffer.data(dk_geom_10_off + 508 * ccomps * dcomps);

            auto g_z_0_xz_xxxxxzz = cbuffer.data(dk_geom_10_off + 509 * ccomps * dcomps);

            auto g_z_0_xz_xxxxyyy = cbuffer.data(dk_geom_10_off + 510 * ccomps * dcomps);

            auto g_z_0_xz_xxxxyyz = cbuffer.data(dk_geom_10_off + 511 * ccomps * dcomps);

            auto g_z_0_xz_xxxxyzz = cbuffer.data(dk_geom_10_off + 512 * ccomps * dcomps);

            auto g_z_0_xz_xxxxzzz = cbuffer.data(dk_geom_10_off + 513 * ccomps * dcomps);

            auto g_z_0_xz_xxxyyyy = cbuffer.data(dk_geom_10_off + 514 * ccomps * dcomps);

            auto g_z_0_xz_xxxyyyz = cbuffer.data(dk_geom_10_off + 515 * ccomps * dcomps);

            auto g_z_0_xz_xxxyyzz = cbuffer.data(dk_geom_10_off + 516 * ccomps * dcomps);

            auto g_z_0_xz_xxxyzzz = cbuffer.data(dk_geom_10_off + 517 * ccomps * dcomps);

            auto g_z_0_xz_xxxzzzz = cbuffer.data(dk_geom_10_off + 518 * ccomps * dcomps);

            auto g_z_0_xz_xxyyyyy = cbuffer.data(dk_geom_10_off + 519 * ccomps * dcomps);

            auto g_z_0_xz_xxyyyyz = cbuffer.data(dk_geom_10_off + 520 * ccomps * dcomps);

            auto g_z_0_xz_xxyyyzz = cbuffer.data(dk_geom_10_off + 521 * ccomps * dcomps);

            auto g_z_0_xz_xxyyzzz = cbuffer.data(dk_geom_10_off + 522 * ccomps * dcomps);

            auto g_z_0_xz_xxyzzzz = cbuffer.data(dk_geom_10_off + 523 * ccomps * dcomps);

            auto g_z_0_xz_xxzzzzz = cbuffer.data(dk_geom_10_off + 524 * ccomps * dcomps);

            auto g_z_0_xz_xyyyyyy = cbuffer.data(dk_geom_10_off + 525 * ccomps * dcomps);

            auto g_z_0_xz_xyyyyyz = cbuffer.data(dk_geom_10_off + 526 * ccomps * dcomps);

            auto g_z_0_xz_xyyyyzz = cbuffer.data(dk_geom_10_off + 527 * ccomps * dcomps);

            auto g_z_0_xz_xyyyzzz = cbuffer.data(dk_geom_10_off + 528 * ccomps * dcomps);

            auto g_z_0_xz_xyyzzzz = cbuffer.data(dk_geom_10_off + 529 * ccomps * dcomps);

            auto g_z_0_xz_xyzzzzz = cbuffer.data(dk_geom_10_off + 530 * ccomps * dcomps);

            auto g_z_0_xz_xzzzzzz = cbuffer.data(dk_geom_10_off + 531 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxxx = cbuffer.data(dk_geom_10_off + 540 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxxy = cbuffer.data(dk_geom_10_off + 541 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxxz = cbuffer.data(dk_geom_10_off + 542 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxyy = cbuffer.data(dk_geom_10_off + 543 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxyz = cbuffer.data(dk_geom_10_off + 544 * ccomps * dcomps);

            auto g_z_0_yy_xxxxxzz = cbuffer.data(dk_geom_10_off + 545 * ccomps * dcomps);

            auto g_z_0_yy_xxxxyyy = cbuffer.data(dk_geom_10_off + 546 * ccomps * dcomps);

            auto g_z_0_yy_xxxxyyz = cbuffer.data(dk_geom_10_off + 547 * ccomps * dcomps);

            auto g_z_0_yy_xxxxyzz = cbuffer.data(dk_geom_10_off + 548 * ccomps * dcomps);

            auto g_z_0_yy_xxxxzzz = cbuffer.data(dk_geom_10_off + 549 * ccomps * dcomps);

            auto g_z_0_yy_xxxyyyy = cbuffer.data(dk_geom_10_off + 550 * ccomps * dcomps);

            auto g_z_0_yy_xxxyyyz = cbuffer.data(dk_geom_10_off + 551 * ccomps * dcomps);

            auto g_z_0_yy_xxxyyzz = cbuffer.data(dk_geom_10_off + 552 * ccomps * dcomps);

            auto g_z_0_yy_xxxyzzz = cbuffer.data(dk_geom_10_off + 553 * ccomps * dcomps);

            auto g_z_0_yy_xxxzzzz = cbuffer.data(dk_geom_10_off + 554 * ccomps * dcomps);

            auto g_z_0_yy_xxyyyyy = cbuffer.data(dk_geom_10_off + 555 * ccomps * dcomps);

            auto g_z_0_yy_xxyyyyz = cbuffer.data(dk_geom_10_off + 556 * ccomps * dcomps);

            auto g_z_0_yy_xxyyyzz = cbuffer.data(dk_geom_10_off + 557 * ccomps * dcomps);

            auto g_z_0_yy_xxyyzzz = cbuffer.data(dk_geom_10_off + 558 * ccomps * dcomps);

            auto g_z_0_yy_xxyzzzz = cbuffer.data(dk_geom_10_off + 559 * ccomps * dcomps);

            auto g_z_0_yy_xxzzzzz = cbuffer.data(dk_geom_10_off + 560 * ccomps * dcomps);

            auto g_z_0_yy_xyyyyyy = cbuffer.data(dk_geom_10_off + 561 * ccomps * dcomps);

            auto g_z_0_yy_xyyyyyz = cbuffer.data(dk_geom_10_off + 562 * ccomps * dcomps);

            auto g_z_0_yy_xyyyyzz = cbuffer.data(dk_geom_10_off + 563 * ccomps * dcomps);

            auto g_z_0_yy_xyyyzzz = cbuffer.data(dk_geom_10_off + 564 * ccomps * dcomps);

            auto g_z_0_yy_xyyzzzz = cbuffer.data(dk_geom_10_off + 565 * ccomps * dcomps);

            auto g_z_0_yy_xyzzzzz = cbuffer.data(dk_geom_10_off + 566 * ccomps * dcomps);

            auto g_z_0_yy_xzzzzzz = cbuffer.data(dk_geom_10_off + 567 * ccomps * dcomps);

            auto g_z_0_yy_yyyyyyy = cbuffer.data(dk_geom_10_off + 568 * ccomps * dcomps);

            auto g_z_0_yy_yyyyyyz = cbuffer.data(dk_geom_10_off + 569 * ccomps * dcomps);

            auto g_z_0_yy_yyyyyzz = cbuffer.data(dk_geom_10_off + 570 * ccomps * dcomps);

            auto g_z_0_yy_yyyyzzz = cbuffer.data(dk_geom_10_off + 571 * ccomps * dcomps);

            auto g_z_0_yy_yyyzzzz = cbuffer.data(dk_geom_10_off + 572 * ccomps * dcomps);

            auto g_z_0_yy_yyzzzzz = cbuffer.data(dk_geom_10_off + 573 * ccomps * dcomps);

            auto g_z_0_yy_yzzzzzz = cbuffer.data(dk_geom_10_off + 574 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxxx = cbuffer.data(dk_geom_10_off + 576 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxxy = cbuffer.data(dk_geom_10_off + 577 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxxz = cbuffer.data(dk_geom_10_off + 578 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxyy = cbuffer.data(dk_geom_10_off + 579 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxyz = cbuffer.data(dk_geom_10_off + 580 * ccomps * dcomps);

            auto g_z_0_yz_xxxxxzz = cbuffer.data(dk_geom_10_off + 581 * ccomps * dcomps);

            auto g_z_0_yz_xxxxyyy = cbuffer.data(dk_geom_10_off + 582 * ccomps * dcomps);

            auto g_z_0_yz_xxxxyyz = cbuffer.data(dk_geom_10_off + 583 * ccomps * dcomps);

            auto g_z_0_yz_xxxxyzz = cbuffer.data(dk_geom_10_off + 584 * ccomps * dcomps);

            auto g_z_0_yz_xxxxzzz = cbuffer.data(dk_geom_10_off + 585 * ccomps * dcomps);

            auto g_z_0_yz_xxxyyyy = cbuffer.data(dk_geom_10_off + 586 * ccomps * dcomps);

            auto g_z_0_yz_xxxyyyz = cbuffer.data(dk_geom_10_off + 587 * ccomps * dcomps);

            auto g_z_0_yz_xxxyyzz = cbuffer.data(dk_geom_10_off + 588 * ccomps * dcomps);

            auto g_z_0_yz_xxxyzzz = cbuffer.data(dk_geom_10_off + 589 * ccomps * dcomps);

            auto g_z_0_yz_xxxzzzz = cbuffer.data(dk_geom_10_off + 590 * ccomps * dcomps);

            auto g_z_0_yz_xxyyyyy = cbuffer.data(dk_geom_10_off + 591 * ccomps * dcomps);

            auto g_z_0_yz_xxyyyyz = cbuffer.data(dk_geom_10_off + 592 * ccomps * dcomps);

            auto g_z_0_yz_xxyyyzz = cbuffer.data(dk_geom_10_off + 593 * ccomps * dcomps);

            auto g_z_0_yz_xxyyzzz = cbuffer.data(dk_geom_10_off + 594 * ccomps * dcomps);

            auto g_z_0_yz_xxyzzzz = cbuffer.data(dk_geom_10_off + 595 * ccomps * dcomps);

            auto g_z_0_yz_xxzzzzz = cbuffer.data(dk_geom_10_off + 596 * ccomps * dcomps);

            auto g_z_0_yz_xyyyyyy = cbuffer.data(dk_geom_10_off + 597 * ccomps * dcomps);

            auto g_z_0_yz_xyyyyyz = cbuffer.data(dk_geom_10_off + 598 * ccomps * dcomps);

            auto g_z_0_yz_xyyyyzz = cbuffer.data(dk_geom_10_off + 599 * ccomps * dcomps);

            auto g_z_0_yz_xyyyzzz = cbuffer.data(dk_geom_10_off + 600 * ccomps * dcomps);

            auto g_z_0_yz_xyyzzzz = cbuffer.data(dk_geom_10_off + 601 * ccomps * dcomps);

            auto g_z_0_yz_xyzzzzz = cbuffer.data(dk_geom_10_off + 602 * ccomps * dcomps);

            auto g_z_0_yz_xzzzzzz = cbuffer.data(dk_geom_10_off + 603 * ccomps * dcomps);

            auto g_z_0_yz_yyyyyyy = cbuffer.data(dk_geom_10_off + 604 * ccomps * dcomps);

            auto g_z_0_yz_yyyyyyz = cbuffer.data(dk_geom_10_off + 605 * ccomps * dcomps);

            auto g_z_0_yz_yyyyyzz = cbuffer.data(dk_geom_10_off + 606 * ccomps * dcomps);

            auto g_z_0_yz_yyyyzzz = cbuffer.data(dk_geom_10_off + 607 * ccomps * dcomps);

            auto g_z_0_yz_yyyzzzz = cbuffer.data(dk_geom_10_off + 608 * ccomps * dcomps);

            auto g_z_0_yz_yyzzzzz = cbuffer.data(dk_geom_10_off + 609 * ccomps * dcomps);

            auto g_z_0_yz_yzzzzzz = cbuffer.data(dk_geom_10_off + 610 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxxx = cbuffer.data(dk_geom_10_off + 612 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxxy = cbuffer.data(dk_geom_10_off + 613 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxxz = cbuffer.data(dk_geom_10_off + 614 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxyy = cbuffer.data(dk_geom_10_off + 615 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxyz = cbuffer.data(dk_geom_10_off + 616 * ccomps * dcomps);

            auto g_z_0_zz_xxxxxzz = cbuffer.data(dk_geom_10_off + 617 * ccomps * dcomps);

            auto g_z_0_zz_xxxxyyy = cbuffer.data(dk_geom_10_off + 618 * ccomps * dcomps);

            auto g_z_0_zz_xxxxyyz = cbuffer.data(dk_geom_10_off + 619 * ccomps * dcomps);

            auto g_z_0_zz_xxxxyzz = cbuffer.data(dk_geom_10_off + 620 * ccomps * dcomps);

            auto g_z_0_zz_xxxxzzz = cbuffer.data(dk_geom_10_off + 621 * ccomps * dcomps);

            auto g_z_0_zz_xxxyyyy = cbuffer.data(dk_geom_10_off + 622 * ccomps * dcomps);

            auto g_z_0_zz_xxxyyyz = cbuffer.data(dk_geom_10_off + 623 * ccomps * dcomps);

            auto g_z_0_zz_xxxyyzz = cbuffer.data(dk_geom_10_off + 624 * ccomps * dcomps);

            auto g_z_0_zz_xxxyzzz = cbuffer.data(dk_geom_10_off + 625 * ccomps * dcomps);

            auto g_z_0_zz_xxxzzzz = cbuffer.data(dk_geom_10_off + 626 * ccomps * dcomps);

            auto g_z_0_zz_xxyyyyy = cbuffer.data(dk_geom_10_off + 627 * ccomps * dcomps);

            auto g_z_0_zz_xxyyyyz = cbuffer.data(dk_geom_10_off + 628 * ccomps * dcomps);

            auto g_z_0_zz_xxyyyzz = cbuffer.data(dk_geom_10_off + 629 * ccomps * dcomps);

            auto g_z_0_zz_xxyyzzz = cbuffer.data(dk_geom_10_off + 630 * ccomps * dcomps);

            auto g_z_0_zz_xxyzzzz = cbuffer.data(dk_geom_10_off + 631 * ccomps * dcomps);

            auto g_z_0_zz_xxzzzzz = cbuffer.data(dk_geom_10_off + 632 * ccomps * dcomps);

            auto g_z_0_zz_xyyyyyy = cbuffer.data(dk_geom_10_off + 633 * ccomps * dcomps);

            auto g_z_0_zz_xyyyyyz = cbuffer.data(dk_geom_10_off + 634 * ccomps * dcomps);

            auto g_z_0_zz_xyyyyzz = cbuffer.data(dk_geom_10_off + 635 * ccomps * dcomps);

            auto g_z_0_zz_xyyyzzz = cbuffer.data(dk_geom_10_off + 636 * ccomps * dcomps);

            auto g_z_0_zz_xyyzzzz = cbuffer.data(dk_geom_10_off + 637 * ccomps * dcomps);

            auto g_z_0_zz_xyzzzzz = cbuffer.data(dk_geom_10_off + 638 * ccomps * dcomps);

            auto g_z_0_zz_xzzzzzz = cbuffer.data(dk_geom_10_off + 639 * ccomps * dcomps);

            auto g_z_0_zz_yyyyyyy = cbuffer.data(dk_geom_10_off + 640 * ccomps * dcomps);

            auto g_z_0_zz_yyyyyyz = cbuffer.data(dk_geom_10_off + 641 * ccomps * dcomps);

            auto g_z_0_zz_yyyyyzz = cbuffer.data(dk_geom_10_off + 642 * ccomps * dcomps);

            auto g_z_0_zz_yyyyzzz = cbuffer.data(dk_geom_10_off + 643 * ccomps * dcomps);

            auto g_z_0_zz_yyyzzzz = cbuffer.data(dk_geom_10_off + 644 * ccomps * dcomps);

            auto g_z_0_zz_yyzzzzz = cbuffer.data(dk_geom_10_off + 645 * ccomps * dcomps);

            auto g_z_0_zz_yzzzzzz = cbuffer.data(dk_geom_10_off + 646 * ccomps * dcomps);

            auto g_z_0_zz_zzzzzzz = cbuffer.data(dk_geom_10_off + 647 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_fixx

            const auto fi_geom_10_off = idx_geom_10_fixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxx_xxxxxx = cbuffer.data(fi_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxxy = cbuffer.data(fi_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxxz = cbuffer.data(fi_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxyy = cbuffer.data(fi_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxyz = cbuffer.data(fi_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxx_xxxxzz = cbuffer.data(fi_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyyy = cbuffer.data(fi_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyyz = cbuffer.data(fi_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxx_xxxyzz = cbuffer.data(fi_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxx_xxxzzz = cbuffer.data(fi_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyyy = cbuffer.data(fi_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyyz = cbuffer.data(fi_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxx_xxyyzz = cbuffer.data(fi_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxx_xxyzzz = cbuffer.data(fi_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxx_xxzzzz = cbuffer.data(fi_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyyy = cbuffer.data(fi_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyyz = cbuffer.data(fi_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxx_xyyyzz = cbuffer.data(fi_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxx_xyyzzz = cbuffer.data(fi_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxx_xyzzzz = cbuffer.data(fi_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxx_xzzzzz = cbuffer.data(fi_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyyy = cbuffer.data(fi_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyyz = cbuffer.data(fi_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxx_yyyyzz = cbuffer.data(fi_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxx_yyyzzz = cbuffer.data(fi_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxx_yyzzzz = cbuffer.data(fi_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxx_yzzzzz = cbuffer.data(fi_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxx_zzzzzz = cbuffer.data(fi_geom_10_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_xxxxxx, g_x_0_xx_xxxxxxx, g_x_0_xx_xxxxxxy, g_x_0_xx_xxxxxxz, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxxyy, g_x_0_xx_xxxxxyz, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxxzz, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyyy, g_x_0_xx_xxxxyyz, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxyzz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxxzzz, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyyy, g_x_0_xx_xxxyyyz, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyyzz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxyzzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxxzzzz, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyyy, g_x_0_xx_xxyyyyz, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyyzz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyyzzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxyzzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xxzzzzz, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyyy, g_x_0_xx_xyyyyyz, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyyzz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyyzzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyyzzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xyzzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_xzzzzzz, g_x_0_xx_yyyyyy, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_zzzzzz, g_x_0_xxx_xxxxxx, g_x_0_xxx_xxxxxy, g_x_0_xxx_xxxxxz, g_x_0_xxx_xxxxyy, g_x_0_xxx_xxxxyz, g_x_0_xxx_xxxxzz, g_x_0_xxx_xxxyyy, g_x_0_xxx_xxxyyz, g_x_0_xxx_xxxyzz, g_x_0_xxx_xxxzzz, g_x_0_xxx_xxyyyy, g_x_0_xxx_xxyyyz, g_x_0_xxx_xxyyzz, g_x_0_xxx_xxyzzz, g_x_0_xxx_xxzzzz, g_x_0_xxx_xyyyyy, g_x_0_xxx_xyyyyz, g_x_0_xxx_xyyyzz, g_x_0_xxx_xyyzzz, g_x_0_xxx_xyzzzz, g_x_0_xxx_xzzzzz, g_x_0_xxx_yyyyyy, g_x_0_xxx_yyyyyz, g_x_0_xxx_yyyyzz, g_x_0_xxx_yyyzzz, g_x_0_xxx_yyzzzz, g_x_0_xxx_yzzzzz, g_x_0_xxx_zzzzzz, g_xx_xxxxxx, g_xx_xxxxxy, g_xx_xxxxxz, g_xx_xxxxyy, g_xx_xxxxyz, g_xx_xxxxzz, g_xx_xxxyyy, g_xx_xxxyyz, g_xx_xxxyzz, g_xx_xxxzzz, g_xx_xxyyyy, g_xx_xxyyyz, g_xx_xxyyzz, g_xx_xxyzzz, g_xx_xxzzzz, g_xx_xyyyyy, g_xx_xyyyyz, g_xx_xyyyzz, g_xx_xyyzzz, g_xx_xyzzzz, g_xx_xzzzzz, g_xx_yyyyyy, g_xx_yyyyyz, g_xx_yyyyzz, g_xx_yyyzzz, g_xx_yyzzzz, g_xx_yzzzzz, g_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxx_xxxxxx[k] = -g_xx_xxxxxx[k] - g_x_0_xx_xxxxxx[k] * ab_x + g_x_0_xx_xxxxxxx[k];

                g_x_0_xxx_xxxxxy[k] = -g_xx_xxxxxy[k] - g_x_0_xx_xxxxxy[k] * ab_x + g_x_0_xx_xxxxxxy[k];

                g_x_0_xxx_xxxxxz[k] = -g_xx_xxxxxz[k] - g_x_0_xx_xxxxxz[k] * ab_x + g_x_0_xx_xxxxxxz[k];

                g_x_0_xxx_xxxxyy[k] = -g_xx_xxxxyy[k] - g_x_0_xx_xxxxyy[k] * ab_x + g_x_0_xx_xxxxxyy[k];

                g_x_0_xxx_xxxxyz[k] = -g_xx_xxxxyz[k] - g_x_0_xx_xxxxyz[k] * ab_x + g_x_0_xx_xxxxxyz[k];

                g_x_0_xxx_xxxxzz[k] = -g_xx_xxxxzz[k] - g_x_0_xx_xxxxzz[k] * ab_x + g_x_0_xx_xxxxxzz[k];

                g_x_0_xxx_xxxyyy[k] = -g_xx_xxxyyy[k] - g_x_0_xx_xxxyyy[k] * ab_x + g_x_0_xx_xxxxyyy[k];

                g_x_0_xxx_xxxyyz[k] = -g_xx_xxxyyz[k] - g_x_0_xx_xxxyyz[k] * ab_x + g_x_0_xx_xxxxyyz[k];

                g_x_0_xxx_xxxyzz[k] = -g_xx_xxxyzz[k] - g_x_0_xx_xxxyzz[k] * ab_x + g_x_0_xx_xxxxyzz[k];

                g_x_0_xxx_xxxzzz[k] = -g_xx_xxxzzz[k] - g_x_0_xx_xxxzzz[k] * ab_x + g_x_0_xx_xxxxzzz[k];

                g_x_0_xxx_xxyyyy[k] = -g_xx_xxyyyy[k] - g_x_0_xx_xxyyyy[k] * ab_x + g_x_0_xx_xxxyyyy[k];

                g_x_0_xxx_xxyyyz[k] = -g_xx_xxyyyz[k] - g_x_0_xx_xxyyyz[k] * ab_x + g_x_0_xx_xxxyyyz[k];

                g_x_0_xxx_xxyyzz[k] = -g_xx_xxyyzz[k] - g_x_0_xx_xxyyzz[k] * ab_x + g_x_0_xx_xxxyyzz[k];

                g_x_0_xxx_xxyzzz[k] = -g_xx_xxyzzz[k] - g_x_0_xx_xxyzzz[k] * ab_x + g_x_0_xx_xxxyzzz[k];

                g_x_0_xxx_xxzzzz[k] = -g_xx_xxzzzz[k] - g_x_0_xx_xxzzzz[k] * ab_x + g_x_0_xx_xxxzzzz[k];

                g_x_0_xxx_xyyyyy[k] = -g_xx_xyyyyy[k] - g_x_0_xx_xyyyyy[k] * ab_x + g_x_0_xx_xxyyyyy[k];

                g_x_0_xxx_xyyyyz[k] = -g_xx_xyyyyz[k] - g_x_0_xx_xyyyyz[k] * ab_x + g_x_0_xx_xxyyyyz[k];

                g_x_0_xxx_xyyyzz[k] = -g_xx_xyyyzz[k] - g_x_0_xx_xyyyzz[k] * ab_x + g_x_0_xx_xxyyyzz[k];

                g_x_0_xxx_xyyzzz[k] = -g_xx_xyyzzz[k] - g_x_0_xx_xyyzzz[k] * ab_x + g_x_0_xx_xxyyzzz[k];

                g_x_0_xxx_xyzzzz[k] = -g_xx_xyzzzz[k] - g_x_0_xx_xyzzzz[k] * ab_x + g_x_0_xx_xxyzzzz[k];

                g_x_0_xxx_xzzzzz[k] = -g_xx_xzzzzz[k] - g_x_0_xx_xzzzzz[k] * ab_x + g_x_0_xx_xxzzzzz[k];

                g_x_0_xxx_yyyyyy[k] = -g_xx_yyyyyy[k] - g_x_0_xx_yyyyyy[k] * ab_x + g_x_0_xx_xyyyyyy[k];

                g_x_0_xxx_yyyyyz[k] = -g_xx_yyyyyz[k] - g_x_0_xx_yyyyyz[k] * ab_x + g_x_0_xx_xyyyyyz[k];

                g_x_0_xxx_yyyyzz[k] = -g_xx_yyyyzz[k] - g_x_0_xx_yyyyzz[k] * ab_x + g_x_0_xx_xyyyyzz[k];

                g_x_0_xxx_yyyzzz[k] = -g_xx_yyyzzz[k] - g_x_0_xx_yyyzzz[k] * ab_x + g_x_0_xx_xyyyzzz[k];

                g_x_0_xxx_yyzzzz[k] = -g_xx_yyzzzz[k] - g_x_0_xx_yyzzzz[k] * ab_x + g_x_0_xx_xyyzzzz[k];

                g_x_0_xxx_yzzzzz[k] = -g_xx_yzzzzz[k] - g_x_0_xx_yzzzzz[k] * ab_x + g_x_0_xx_xyzzzzz[k];

                g_x_0_xxx_zzzzzz[k] = -g_xx_zzzzzz[k] - g_x_0_xx_zzzzzz[k] * ab_x + g_x_0_xx_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxy_xxxxxx = cbuffer.data(fi_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxy = cbuffer.data(fi_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxxz = cbuffer.data(fi_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxyy = cbuffer.data(fi_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxyz = cbuffer.data(fi_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxy_xxxxzz = cbuffer.data(fi_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyyy = cbuffer.data(fi_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyyz = cbuffer.data(fi_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxy_xxxyzz = cbuffer.data(fi_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxy_xxxzzz = cbuffer.data(fi_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyyy = cbuffer.data(fi_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyyz = cbuffer.data(fi_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxy_xxyyzz = cbuffer.data(fi_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxy_xxyzzz = cbuffer.data(fi_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxy_xxzzzz = cbuffer.data(fi_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyyy = cbuffer.data(fi_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyyz = cbuffer.data(fi_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxy_xyyyzz = cbuffer.data(fi_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxy_xyyzzz = cbuffer.data(fi_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxy_xyzzzz = cbuffer.data(fi_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxy_xzzzzz = cbuffer.data(fi_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyyy = cbuffer.data(fi_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyyz = cbuffer.data(fi_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxy_yyyyzz = cbuffer.data(fi_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxy_yyyzzz = cbuffer.data(fi_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxy_yyzzzz = cbuffer.data(fi_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxy_yzzzzz = cbuffer.data(fi_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxy_zzzzzz = cbuffer.data(fi_geom_10_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_xxxxxx, g_x_0_xx_xxxxxxy, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxxyy, g_x_0_xx_xxxxxyz, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyyy, g_x_0_xx_xxxxyyz, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxyzz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyyy, g_x_0_xx_xxxyyyz, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyyzz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxyzzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyyy, g_x_0_xx_xxyyyyz, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyyzz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyyzzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxyzzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyyy, g_x_0_xx_xyyyyyz, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyyzz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyyzzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyyzzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xyzzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_yyyyyy, g_x_0_xx_yyyyyyy, g_x_0_xx_yyyyyyz, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyyzz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyyzzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyyzzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yyzzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_yzzzzzz, g_x_0_xx_zzzzzz, g_x_0_xxy_xxxxxx, g_x_0_xxy_xxxxxy, g_x_0_xxy_xxxxxz, g_x_0_xxy_xxxxyy, g_x_0_xxy_xxxxyz, g_x_0_xxy_xxxxzz, g_x_0_xxy_xxxyyy, g_x_0_xxy_xxxyyz, g_x_0_xxy_xxxyzz, g_x_0_xxy_xxxzzz, g_x_0_xxy_xxyyyy, g_x_0_xxy_xxyyyz, g_x_0_xxy_xxyyzz, g_x_0_xxy_xxyzzz, g_x_0_xxy_xxzzzz, g_x_0_xxy_xyyyyy, g_x_0_xxy_xyyyyz, g_x_0_xxy_xyyyzz, g_x_0_xxy_xyyzzz, g_x_0_xxy_xyzzzz, g_x_0_xxy_xzzzzz, g_x_0_xxy_yyyyyy, g_x_0_xxy_yyyyyz, g_x_0_xxy_yyyyzz, g_x_0_xxy_yyyzzz, g_x_0_xxy_yyzzzz, g_x_0_xxy_yzzzzz, g_x_0_xxy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxy_xxxxxx[k] = -g_x_0_xx_xxxxxx[k] * ab_y + g_x_0_xx_xxxxxxy[k];

                g_x_0_xxy_xxxxxy[k] = -g_x_0_xx_xxxxxy[k] * ab_y + g_x_0_xx_xxxxxyy[k];

                g_x_0_xxy_xxxxxz[k] = -g_x_0_xx_xxxxxz[k] * ab_y + g_x_0_xx_xxxxxyz[k];

                g_x_0_xxy_xxxxyy[k] = -g_x_0_xx_xxxxyy[k] * ab_y + g_x_0_xx_xxxxyyy[k];

                g_x_0_xxy_xxxxyz[k] = -g_x_0_xx_xxxxyz[k] * ab_y + g_x_0_xx_xxxxyyz[k];

                g_x_0_xxy_xxxxzz[k] = -g_x_0_xx_xxxxzz[k] * ab_y + g_x_0_xx_xxxxyzz[k];

                g_x_0_xxy_xxxyyy[k] = -g_x_0_xx_xxxyyy[k] * ab_y + g_x_0_xx_xxxyyyy[k];

                g_x_0_xxy_xxxyyz[k] = -g_x_0_xx_xxxyyz[k] * ab_y + g_x_0_xx_xxxyyyz[k];

                g_x_0_xxy_xxxyzz[k] = -g_x_0_xx_xxxyzz[k] * ab_y + g_x_0_xx_xxxyyzz[k];

                g_x_0_xxy_xxxzzz[k] = -g_x_0_xx_xxxzzz[k] * ab_y + g_x_0_xx_xxxyzzz[k];

                g_x_0_xxy_xxyyyy[k] = -g_x_0_xx_xxyyyy[k] * ab_y + g_x_0_xx_xxyyyyy[k];

                g_x_0_xxy_xxyyyz[k] = -g_x_0_xx_xxyyyz[k] * ab_y + g_x_0_xx_xxyyyyz[k];

                g_x_0_xxy_xxyyzz[k] = -g_x_0_xx_xxyyzz[k] * ab_y + g_x_0_xx_xxyyyzz[k];

                g_x_0_xxy_xxyzzz[k] = -g_x_0_xx_xxyzzz[k] * ab_y + g_x_0_xx_xxyyzzz[k];

                g_x_0_xxy_xxzzzz[k] = -g_x_0_xx_xxzzzz[k] * ab_y + g_x_0_xx_xxyzzzz[k];

                g_x_0_xxy_xyyyyy[k] = -g_x_0_xx_xyyyyy[k] * ab_y + g_x_0_xx_xyyyyyy[k];

                g_x_0_xxy_xyyyyz[k] = -g_x_0_xx_xyyyyz[k] * ab_y + g_x_0_xx_xyyyyyz[k];

                g_x_0_xxy_xyyyzz[k] = -g_x_0_xx_xyyyzz[k] * ab_y + g_x_0_xx_xyyyyzz[k];

                g_x_0_xxy_xyyzzz[k] = -g_x_0_xx_xyyzzz[k] * ab_y + g_x_0_xx_xyyyzzz[k];

                g_x_0_xxy_xyzzzz[k] = -g_x_0_xx_xyzzzz[k] * ab_y + g_x_0_xx_xyyzzzz[k];

                g_x_0_xxy_xzzzzz[k] = -g_x_0_xx_xzzzzz[k] * ab_y + g_x_0_xx_xyzzzzz[k];

                g_x_0_xxy_yyyyyy[k] = -g_x_0_xx_yyyyyy[k] * ab_y + g_x_0_xx_yyyyyyy[k];

                g_x_0_xxy_yyyyyz[k] = -g_x_0_xx_yyyyyz[k] * ab_y + g_x_0_xx_yyyyyyz[k];

                g_x_0_xxy_yyyyzz[k] = -g_x_0_xx_yyyyzz[k] * ab_y + g_x_0_xx_yyyyyzz[k];

                g_x_0_xxy_yyyzzz[k] = -g_x_0_xx_yyyzzz[k] * ab_y + g_x_0_xx_yyyyzzz[k];

                g_x_0_xxy_yyzzzz[k] = -g_x_0_xx_yyzzzz[k] * ab_y + g_x_0_xx_yyyzzzz[k];

                g_x_0_xxy_yzzzzz[k] = -g_x_0_xx_yzzzzz[k] * ab_y + g_x_0_xx_yyzzzzz[k];

                g_x_0_xxy_zzzzzz[k] = -g_x_0_xx_zzzzzz[k] * ab_y + g_x_0_xx_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxz_xxxxxx = cbuffer.data(fi_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxy = cbuffer.data(fi_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxxz = cbuffer.data(fi_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxyy = cbuffer.data(fi_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxyz = cbuffer.data(fi_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxz_xxxxzz = cbuffer.data(fi_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyyy = cbuffer.data(fi_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyyz = cbuffer.data(fi_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxz_xxxyzz = cbuffer.data(fi_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxz_xxxzzz = cbuffer.data(fi_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyyy = cbuffer.data(fi_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyyz = cbuffer.data(fi_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxz_xxyyzz = cbuffer.data(fi_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxz_xxyzzz = cbuffer.data(fi_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxz_xxzzzz = cbuffer.data(fi_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyyy = cbuffer.data(fi_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyyz = cbuffer.data(fi_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxz_xyyyzz = cbuffer.data(fi_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxz_xyyzzz = cbuffer.data(fi_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxz_xyzzzz = cbuffer.data(fi_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxz_xzzzzz = cbuffer.data(fi_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyyy = cbuffer.data(fi_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyyz = cbuffer.data(fi_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxz_yyyyzz = cbuffer.data(fi_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxz_yyyzzz = cbuffer.data(fi_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxz_yyzzzz = cbuffer.data(fi_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxz_yzzzzz = cbuffer.data(fi_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxz_zzzzzz = cbuffer.data(fi_geom_10_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xx_xxxxxx, g_x_0_xx_xxxxxxz, g_x_0_xx_xxxxxy, g_x_0_xx_xxxxxyz, g_x_0_xx_xxxxxz, g_x_0_xx_xxxxxzz, g_x_0_xx_xxxxyy, g_x_0_xx_xxxxyyz, g_x_0_xx_xxxxyz, g_x_0_xx_xxxxyzz, g_x_0_xx_xxxxzz, g_x_0_xx_xxxxzzz, g_x_0_xx_xxxyyy, g_x_0_xx_xxxyyyz, g_x_0_xx_xxxyyz, g_x_0_xx_xxxyyzz, g_x_0_xx_xxxyzz, g_x_0_xx_xxxyzzz, g_x_0_xx_xxxzzz, g_x_0_xx_xxxzzzz, g_x_0_xx_xxyyyy, g_x_0_xx_xxyyyyz, g_x_0_xx_xxyyyz, g_x_0_xx_xxyyyzz, g_x_0_xx_xxyyzz, g_x_0_xx_xxyyzzz, g_x_0_xx_xxyzzz, g_x_0_xx_xxyzzzz, g_x_0_xx_xxzzzz, g_x_0_xx_xxzzzzz, g_x_0_xx_xyyyyy, g_x_0_xx_xyyyyyz, g_x_0_xx_xyyyyz, g_x_0_xx_xyyyyzz, g_x_0_xx_xyyyzz, g_x_0_xx_xyyyzzz, g_x_0_xx_xyyzzz, g_x_0_xx_xyyzzzz, g_x_0_xx_xyzzzz, g_x_0_xx_xyzzzzz, g_x_0_xx_xzzzzz, g_x_0_xx_xzzzzzz, g_x_0_xx_yyyyyy, g_x_0_xx_yyyyyyz, g_x_0_xx_yyyyyz, g_x_0_xx_yyyyyzz, g_x_0_xx_yyyyzz, g_x_0_xx_yyyyzzz, g_x_0_xx_yyyzzz, g_x_0_xx_yyyzzzz, g_x_0_xx_yyzzzz, g_x_0_xx_yyzzzzz, g_x_0_xx_yzzzzz, g_x_0_xx_yzzzzzz, g_x_0_xx_zzzzzz, g_x_0_xx_zzzzzzz, g_x_0_xxz_xxxxxx, g_x_0_xxz_xxxxxy, g_x_0_xxz_xxxxxz, g_x_0_xxz_xxxxyy, g_x_0_xxz_xxxxyz, g_x_0_xxz_xxxxzz, g_x_0_xxz_xxxyyy, g_x_0_xxz_xxxyyz, g_x_0_xxz_xxxyzz, g_x_0_xxz_xxxzzz, g_x_0_xxz_xxyyyy, g_x_0_xxz_xxyyyz, g_x_0_xxz_xxyyzz, g_x_0_xxz_xxyzzz, g_x_0_xxz_xxzzzz, g_x_0_xxz_xyyyyy, g_x_0_xxz_xyyyyz, g_x_0_xxz_xyyyzz, g_x_0_xxz_xyyzzz, g_x_0_xxz_xyzzzz, g_x_0_xxz_xzzzzz, g_x_0_xxz_yyyyyy, g_x_0_xxz_yyyyyz, g_x_0_xxz_yyyyzz, g_x_0_xxz_yyyzzz, g_x_0_xxz_yyzzzz, g_x_0_xxz_yzzzzz, g_x_0_xxz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxz_xxxxxx[k] = -g_x_0_xx_xxxxxx[k] * ab_z + g_x_0_xx_xxxxxxz[k];

                g_x_0_xxz_xxxxxy[k] = -g_x_0_xx_xxxxxy[k] * ab_z + g_x_0_xx_xxxxxyz[k];

                g_x_0_xxz_xxxxxz[k] = -g_x_0_xx_xxxxxz[k] * ab_z + g_x_0_xx_xxxxxzz[k];

                g_x_0_xxz_xxxxyy[k] = -g_x_0_xx_xxxxyy[k] * ab_z + g_x_0_xx_xxxxyyz[k];

                g_x_0_xxz_xxxxyz[k] = -g_x_0_xx_xxxxyz[k] * ab_z + g_x_0_xx_xxxxyzz[k];

                g_x_0_xxz_xxxxzz[k] = -g_x_0_xx_xxxxzz[k] * ab_z + g_x_0_xx_xxxxzzz[k];

                g_x_0_xxz_xxxyyy[k] = -g_x_0_xx_xxxyyy[k] * ab_z + g_x_0_xx_xxxyyyz[k];

                g_x_0_xxz_xxxyyz[k] = -g_x_0_xx_xxxyyz[k] * ab_z + g_x_0_xx_xxxyyzz[k];

                g_x_0_xxz_xxxyzz[k] = -g_x_0_xx_xxxyzz[k] * ab_z + g_x_0_xx_xxxyzzz[k];

                g_x_0_xxz_xxxzzz[k] = -g_x_0_xx_xxxzzz[k] * ab_z + g_x_0_xx_xxxzzzz[k];

                g_x_0_xxz_xxyyyy[k] = -g_x_0_xx_xxyyyy[k] * ab_z + g_x_0_xx_xxyyyyz[k];

                g_x_0_xxz_xxyyyz[k] = -g_x_0_xx_xxyyyz[k] * ab_z + g_x_0_xx_xxyyyzz[k];

                g_x_0_xxz_xxyyzz[k] = -g_x_0_xx_xxyyzz[k] * ab_z + g_x_0_xx_xxyyzzz[k];

                g_x_0_xxz_xxyzzz[k] = -g_x_0_xx_xxyzzz[k] * ab_z + g_x_0_xx_xxyzzzz[k];

                g_x_0_xxz_xxzzzz[k] = -g_x_0_xx_xxzzzz[k] * ab_z + g_x_0_xx_xxzzzzz[k];

                g_x_0_xxz_xyyyyy[k] = -g_x_0_xx_xyyyyy[k] * ab_z + g_x_0_xx_xyyyyyz[k];

                g_x_0_xxz_xyyyyz[k] = -g_x_0_xx_xyyyyz[k] * ab_z + g_x_0_xx_xyyyyzz[k];

                g_x_0_xxz_xyyyzz[k] = -g_x_0_xx_xyyyzz[k] * ab_z + g_x_0_xx_xyyyzzz[k];

                g_x_0_xxz_xyyzzz[k] = -g_x_0_xx_xyyzzz[k] * ab_z + g_x_0_xx_xyyzzzz[k];

                g_x_0_xxz_xyzzzz[k] = -g_x_0_xx_xyzzzz[k] * ab_z + g_x_0_xx_xyzzzzz[k];

                g_x_0_xxz_xzzzzz[k] = -g_x_0_xx_xzzzzz[k] * ab_z + g_x_0_xx_xzzzzzz[k];

                g_x_0_xxz_yyyyyy[k] = -g_x_0_xx_yyyyyy[k] * ab_z + g_x_0_xx_yyyyyyz[k];

                g_x_0_xxz_yyyyyz[k] = -g_x_0_xx_yyyyyz[k] * ab_z + g_x_0_xx_yyyyyzz[k];

                g_x_0_xxz_yyyyzz[k] = -g_x_0_xx_yyyyzz[k] * ab_z + g_x_0_xx_yyyyzzz[k];

                g_x_0_xxz_yyyzzz[k] = -g_x_0_xx_yyyzzz[k] * ab_z + g_x_0_xx_yyyzzzz[k];

                g_x_0_xxz_yyzzzz[k] = -g_x_0_xx_yyzzzz[k] * ab_z + g_x_0_xx_yyzzzzz[k];

                g_x_0_xxz_yzzzzz[k] = -g_x_0_xx_yzzzzz[k] * ab_z + g_x_0_xx_yzzzzzz[k];

                g_x_0_xxz_zzzzzz[k] = -g_x_0_xx_zzzzzz[k] * ab_z + g_x_0_xx_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyy_xxxxxx = cbuffer.data(fi_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxy = cbuffer.data(fi_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxxz = cbuffer.data(fi_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxyy = cbuffer.data(fi_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxyz = cbuffer.data(fi_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xyy_xxxxzz = cbuffer.data(fi_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyyy = cbuffer.data(fi_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyyz = cbuffer.data(fi_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xyy_xxxyzz = cbuffer.data(fi_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xyy_xxxzzz = cbuffer.data(fi_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyyy = cbuffer.data(fi_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyyz = cbuffer.data(fi_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xyy_xxyyzz = cbuffer.data(fi_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xyy_xxyzzz = cbuffer.data(fi_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xyy_xxzzzz = cbuffer.data(fi_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyyy = cbuffer.data(fi_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyyz = cbuffer.data(fi_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xyy_xyyyzz = cbuffer.data(fi_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xyy_xyyzzz = cbuffer.data(fi_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xyy_xyzzzz = cbuffer.data(fi_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xyy_xzzzzz = cbuffer.data(fi_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyyy = cbuffer.data(fi_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyyz = cbuffer.data(fi_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xyy_yyyyzz = cbuffer.data(fi_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xyy_yyyzzz = cbuffer.data(fi_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xyy_yyzzzz = cbuffer.data(fi_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xyy_yzzzzz = cbuffer.data(fi_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xyy_zzzzzz = cbuffer.data(fi_geom_10_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xy_xxxxxx, g_x_0_xy_xxxxxxy, g_x_0_xy_xxxxxy, g_x_0_xy_xxxxxyy, g_x_0_xy_xxxxxyz, g_x_0_xy_xxxxxz, g_x_0_xy_xxxxyy, g_x_0_xy_xxxxyyy, g_x_0_xy_xxxxyyz, g_x_0_xy_xxxxyz, g_x_0_xy_xxxxyzz, g_x_0_xy_xxxxzz, g_x_0_xy_xxxyyy, g_x_0_xy_xxxyyyy, g_x_0_xy_xxxyyyz, g_x_0_xy_xxxyyz, g_x_0_xy_xxxyyzz, g_x_0_xy_xxxyzz, g_x_0_xy_xxxyzzz, g_x_0_xy_xxxzzz, g_x_0_xy_xxyyyy, g_x_0_xy_xxyyyyy, g_x_0_xy_xxyyyyz, g_x_0_xy_xxyyyz, g_x_0_xy_xxyyyzz, g_x_0_xy_xxyyzz, g_x_0_xy_xxyyzzz, g_x_0_xy_xxyzzz, g_x_0_xy_xxyzzzz, g_x_0_xy_xxzzzz, g_x_0_xy_xyyyyy, g_x_0_xy_xyyyyyy, g_x_0_xy_xyyyyyz, g_x_0_xy_xyyyyz, g_x_0_xy_xyyyyzz, g_x_0_xy_xyyyzz, g_x_0_xy_xyyyzzz, g_x_0_xy_xyyzzz, g_x_0_xy_xyyzzzz, g_x_0_xy_xyzzzz, g_x_0_xy_xyzzzzz, g_x_0_xy_xzzzzz, g_x_0_xy_yyyyyy, g_x_0_xy_yyyyyyy, g_x_0_xy_yyyyyyz, g_x_0_xy_yyyyyz, g_x_0_xy_yyyyyzz, g_x_0_xy_yyyyzz, g_x_0_xy_yyyyzzz, g_x_0_xy_yyyzzz, g_x_0_xy_yyyzzzz, g_x_0_xy_yyzzzz, g_x_0_xy_yyzzzzz, g_x_0_xy_yzzzzz, g_x_0_xy_yzzzzzz, g_x_0_xy_zzzzzz, g_x_0_xyy_xxxxxx, g_x_0_xyy_xxxxxy, g_x_0_xyy_xxxxxz, g_x_0_xyy_xxxxyy, g_x_0_xyy_xxxxyz, g_x_0_xyy_xxxxzz, g_x_0_xyy_xxxyyy, g_x_0_xyy_xxxyyz, g_x_0_xyy_xxxyzz, g_x_0_xyy_xxxzzz, g_x_0_xyy_xxyyyy, g_x_0_xyy_xxyyyz, g_x_0_xyy_xxyyzz, g_x_0_xyy_xxyzzz, g_x_0_xyy_xxzzzz, g_x_0_xyy_xyyyyy, g_x_0_xyy_xyyyyz, g_x_0_xyy_xyyyzz, g_x_0_xyy_xyyzzz, g_x_0_xyy_xyzzzz, g_x_0_xyy_xzzzzz, g_x_0_xyy_yyyyyy, g_x_0_xyy_yyyyyz, g_x_0_xyy_yyyyzz, g_x_0_xyy_yyyzzz, g_x_0_xyy_yyzzzz, g_x_0_xyy_yzzzzz, g_x_0_xyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyy_xxxxxx[k] = -g_x_0_xy_xxxxxx[k] * ab_y + g_x_0_xy_xxxxxxy[k];

                g_x_0_xyy_xxxxxy[k] = -g_x_0_xy_xxxxxy[k] * ab_y + g_x_0_xy_xxxxxyy[k];

                g_x_0_xyy_xxxxxz[k] = -g_x_0_xy_xxxxxz[k] * ab_y + g_x_0_xy_xxxxxyz[k];

                g_x_0_xyy_xxxxyy[k] = -g_x_0_xy_xxxxyy[k] * ab_y + g_x_0_xy_xxxxyyy[k];

                g_x_0_xyy_xxxxyz[k] = -g_x_0_xy_xxxxyz[k] * ab_y + g_x_0_xy_xxxxyyz[k];

                g_x_0_xyy_xxxxzz[k] = -g_x_0_xy_xxxxzz[k] * ab_y + g_x_0_xy_xxxxyzz[k];

                g_x_0_xyy_xxxyyy[k] = -g_x_0_xy_xxxyyy[k] * ab_y + g_x_0_xy_xxxyyyy[k];

                g_x_0_xyy_xxxyyz[k] = -g_x_0_xy_xxxyyz[k] * ab_y + g_x_0_xy_xxxyyyz[k];

                g_x_0_xyy_xxxyzz[k] = -g_x_0_xy_xxxyzz[k] * ab_y + g_x_0_xy_xxxyyzz[k];

                g_x_0_xyy_xxxzzz[k] = -g_x_0_xy_xxxzzz[k] * ab_y + g_x_0_xy_xxxyzzz[k];

                g_x_0_xyy_xxyyyy[k] = -g_x_0_xy_xxyyyy[k] * ab_y + g_x_0_xy_xxyyyyy[k];

                g_x_0_xyy_xxyyyz[k] = -g_x_0_xy_xxyyyz[k] * ab_y + g_x_0_xy_xxyyyyz[k];

                g_x_0_xyy_xxyyzz[k] = -g_x_0_xy_xxyyzz[k] * ab_y + g_x_0_xy_xxyyyzz[k];

                g_x_0_xyy_xxyzzz[k] = -g_x_0_xy_xxyzzz[k] * ab_y + g_x_0_xy_xxyyzzz[k];

                g_x_0_xyy_xxzzzz[k] = -g_x_0_xy_xxzzzz[k] * ab_y + g_x_0_xy_xxyzzzz[k];

                g_x_0_xyy_xyyyyy[k] = -g_x_0_xy_xyyyyy[k] * ab_y + g_x_0_xy_xyyyyyy[k];

                g_x_0_xyy_xyyyyz[k] = -g_x_0_xy_xyyyyz[k] * ab_y + g_x_0_xy_xyyyyyz[k];

                g_x_0_xyy_xyyyzz[k] = -g_x_0_xy_xyyyzz[k] * ab_y + g_x_0_xy_xyyyyzz[k];

                g_x_0_xyy_xyyzzz[k] = -g_x_0_xy_xyyzzz[k] * ab_y + g_x_0_xy_xyyyzzz[k];

                g_x_0_xyy_xyzzzz[k] = -g_x_0_xy_xyzzzz[k] * ab_y + g_x_0_xy_xyyzzzz[k];

                g_x_0_xyy_xzzzzz[k] = -g_x_0_xy_xzzzzz[k] * ab_y + g_x_0_xy_xyzzzzz[k];

                g_x_0_xyy_yyyyyy[k] = -g_x_0_xy_yyyyyy[k] * ab_y + g_x_0_xy_yyyyyyy[k];

                g_x_0_xyy_yyyyyz[k] = -g_x_0_xy_yyyyyz[k] * ab_y + g_x_0_xy_yyyyyyz[k];

                g_x_0_xyy_yyyyzz[k] = -g_x_0_xy_yyyyzz[k] * ab_y + g_x_0_xy_yyyyyzz[k];

                g_x_0_xyy_yyyzzz[k] = -g_x_0_xy_yyyzzz[k] * ab_y + g_x_0_xy_yyyyzzz[k];

                g_x_0_xyy_yyzzzz[k] = -g_x_0_xy_yyzzzz[k] * ab_y + g_x_0_xy_yyyzzzz[k];

                g_x_0_xyy_yzzzzz[k] = -g_x_0_xy_yzzzzz[k] * ab_y + g_x_0_xy_yyzzzzz[k];

                g_x_0_xyy_zzzzzz[k] = -g_x_0_xy_zzzzzz[k] * ab_y + g_x_0_xy_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyz_xxxxxx = cbuffer.data(fi_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxy = cbuffer.data(fi_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxxz = cbuffer.data(fi_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxyy = cbuffer.data(fi_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxyz = cbuffer.data(fi_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xyz_xxxxzz = cbuffer.data(fi_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyyy = cbuffer.data(fi_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyyz = cbuffer.data(fi_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xyz_xxxyzz = cbuffer.data(fi_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xyz_xxxzzz = cbuffer.data(fi_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyyy = cbuffer.data(fi_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyyz = cbuffer.data(fi_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xyz_xxyyzz = cbuffer.data(fi_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xyz_xxyzzz = cbuffer.data(fi_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xyz_xxzzzz = cbuffer.data(fi_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyyy = cbuffer.data(fi_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyyz = cbuffer.data(fi_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xyz_xyyyzz = cbuffer.data(fi_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xyz_xyyzzz = cbuffer.data(fi_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xyz_xyzzzz = cbuffer.data(fi_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xyz_xzzzzz = cbuffer.data(fi_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyyy = cbuffer.data(fi_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyyz = cbuffer.data(fi_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xyz_yyyyzz = cbuffer.data(fi_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xyz_yyyzzz = cbuffer.data(fi_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xyz_yyzzzz = cbuffer.data(fi_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xyz_yzzzzz = cbuffer.data(fi_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xyz_zzzzzz = cbuffer.data(fi_geom_10_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyz_xxxxxx, g_x_0_xyz_xxxxxy, g_x_0_xyz_xxxxxz, g_x_0_xyz_xxxxyy, g_x_0_xyz_xxxxyz, g_x_0_xyz_xxxxzz, g_x_0_xyz_xxxyyy, g_x_0_xyz_xxxyyz, g_x_0_xyz_xxxyzz, g_x_0_xyz_xxxzzz, g_x_0_xyz_xxyyyy, g_x_0_xyz_xxyyyz, g_x_0_xyz_xxyyzz, g_x_0_xyz_xxyzzz, g_x_0_xyz_xxzzzz, g_x_0_xyz_xyyyyy, g_x_0_xyz_xyyyyz, g_x_0_xyz_xyyyzz, g_x_0_xyz_xyyzzz, g_x_0_xyz_xyzzzz, g_x_0_xyz_xzzzzz, g_x_0_xyz_yyyyyy, g_x_0_xyz_yyyyyz, g_x_0_xyz_yyyyzz, g_x_0_xyz_yyyzzz, g_x_0_xyz_yyzzzz, g_x_0_xyz_yzzzzz, g_x_0_xyz_zzzzzz, g_x_0_xz_xxxxxx, g_x_0_xz_xxxxxxy, g_x_0_xz_xxxxxy, g_x_0_xz_xxxxxyy, g_x_0_xz_xxxxxyz, g_x_0_xz_xxxxxz, g_x_0_xz_xxxxyy, g_x_0_xz_xxxxyyy, g_x_0_xz_xxxxyyz, g_x_0_xz_xxxxyz, g_x_0_xz_xxxxyzz, g_x_0_xz_xxxxzz, g_x_0_xz_xxxyyy, g_x_0_xz_xxxyyyy, g_x_0_xz_xxxyyyz, g_x_0_xz_xxxyyz, g_x_0_xz_xxxyyzz, g_x_0_xz_xxxyzz, g_x_0_xz_xxxyzzz, g_x_0_xz_xxxzzz, g_x_0_xz_xxyyyy, g_x_0_xz_xxyyyyy, g_x_0_xz_xxyyyyz, g_x_0_xz_xxyyyz, g_x_0_xz_xxyyyzz, g_x_0_xz_xxyyzz, g_x_0_xz_xxyyzzz, g_x_0_xz_xxyzzz, g_x_0_xz_xxyzzzz, g_x_0_xz_xxzzzz, g_x_0_xz_xyyyyy, g_x_0_xz_xyyyyyy, g_x_0_xz_xyyyyyz, g_x_0_xz_xyyyyz, g_x_0_xz_xyyyyzz, g_x_0_xz_xyyyzz, g_x_0_xz_xyyyzzz, g_x_0_xz_xyyzzz, g_x_0_xz_xyyzzzz, g_x_0_xz_xyzzzz, g_x_0_xz_xyzzzzz, g_x_0_xz_xzzzzz, g_x_0_xz_yyyyyy, g_x_0_xz_yyyyyyy, g_x_0_xz_yyyyyyz, g_x_0_xz_yyyyyz, g_x_0_xz_yyyyyzz, g_x_0_xz_yyyyzz, g_x_0_xz_yyyyzzz, g_x_0_xz_yyyzzz, g_x_0_xz_yyyzzzz, g_x_0_xz_yyzzzz, g_x_0_xz_yyzzzzz, g_x_0_xz_yzzzzz, g_x_0_xz_yzzzzzz, g_x_0_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyz_xxxxxx[k] = -g_x_0_xz_xxxxxx[k] * ab_y + g_x_0_xz_xxxxxxy[k];

                g_x_0_xyz_xxxxxy[k] = -g_x_0_xz_xxxxxy[k] * ab_y + g_x_0_xz_xxxxxyy[k];

                g_x_0_xyz_xxxxxz[k] = -g_x_0_xz_xxxxxz[k] * ab_y + g_x_0_xz_xxxxxyz[k];

                g_x_0_xyz_xxxxyy[k] = -g_x_0_xz_xxxxyy[k] * ab_y + g_x_0_xz_xxxxyyy[k];

                g_x_0_xyz_xxxxyz[k] = -g_x_0_xz_xxxxyz[k] * ab_y + g_x_0_xz_xxxxyyz[k];

                g_x_0_xyz_xxxxzz[k] = -g_x_0_xz_xxxxzz[k] * ab_y + g_x_0_xz_xxxxyzz[k];

                g_x_0_xyz_xxxyyy[k] = -g_x_0_xz_xxxyyy[k] * ab_y + g_x_0_xz_xxxyyyy[k];

                g_x_0_xyz_xxxyyz[k] = -g_x_0_xz_xxxyyz[k] * ab_y + g_x_0_xz_xxxyyyz[k];

                g_x_0_xyz_xxxyzz[k] = -g_x_0_xz_xxxyzz[k] * ab_y + g_x_0_xz_xxxyyzz[k];

                g_x_0_xyz_xxxzzz[k] = -g_x_0_xz_xxxzzz[k] * ab_y + g_x_0_xz_xxxyzzz[k];

                g_x_0_xyz_xxyyyy[k] = -g_x_0_xz_xxyyyy[k] * ab_y + g_x_0_xz_xxyyyyy[k];

                g_x_0_xyz_xxyyyz[k] = -g_x_0_xz_xxyyyz[k] * ab_y + g_x_0_xz_xxyyyyz[k];

                g_x_0_xyz_xxyyzz[k] = -g_x_0_xz_xxyyzz[k] * ab_y + g_x_0_xz_xxyyyzz[k];

                g_x_0_xyz_xxyzzz[k] = -g_x_0_xz_xxyzzz[k] * ab_y + g_x_0_xz_xxyyzzz[k];

                g_x_0_xyz_xxzzzz[k] = -g_x_0_xz_xxzzzz[k] * ab_y + g_x_0_xz_xxyzzzz[k];

                g_x_0_xyz_xyyyyy[k] = -g_x_0_xz_xyyyyy[k] * ab_y + g_x_0_xz_xyyyyyy[k];

                g_x_0_xyz_xyyyyz[k] = -g_x_0_xz_xyyyyz[k] * ab_y + g_x_0_xz_xyyyyyz[k];

                g_x_0_xyz_xyyyzz[k] = -g_x_0_xz_xyyyzz[k] * ab_y + g_x_0_xz_xyyyyzz[k];

                g_x_0_xyz_xyyzzz[k] = -g_x_0_xz_xyyzzz[k] * ab_y + g_x_0_xz_xyyyzzz[k];

                g_x_0_xyz_xyzzzz[k] = -g_x_0_xz_xyzzzz[k] * ab_y + g_x_0_xz_xyyzzzz[k];

                g_x_0_xyz_xzzzzz[k] = -g_x_0_xz_xzzzzz[k] * ab_y + g_x_0_xz_xyzzzzz[k];

                g_x_0_xyz_yyyyyy[k] = -g_x_0_xz_yyyyyy[k] * ab_y + g_x_0_xz_yyyyyyy[k];

                g_x_0_xyz_yyyyyz[k] = -g_x_0_xz_yyyyyz[k] * ab_y + g_x_0_xz_yyyyyyz[k];

                g_x_0_xyz_yyyyzz[k] = -g_x_0_xz_yyyyzz[k] * ab_y + g_x_0_xz_yyyyyzz[k];

                g_x_0_xyz_yyyzzz[k] = -g_x_0_xz_yyyzzz[k] * ab_y + g_x_0_xz_yyyyzzz[k];

                g_x_0_xyz_yyzzzz[k] = -g_x_0_xz_yyzzzz[k] * ab_y + g_x_0_xz_yyyzzzz[k];

                g_x_0_xyz_yzzzzz[k] = -g_x_0_xz_yzzzzz[k] * ab_y + g_x_0_xz_yyzzzzz[k];

                g_x_0_xyz_zzzzzz[k] = -g_x_0_xz_zzzzzz[k] * ab_y + g_x_0_xz_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzz_xxxxxx = cbuffer.data(fi_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxy = cbuffer.data(fi_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxxz = cbuffer.data(fi_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxyy = cbuffer.data(fi_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxyz = cbuffer.data(fi_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xzz_xxxxzz = cbuffer.data(fi_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyyy = cbuffer.data(fi_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyyz = cbuffer.data(fi_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xzz_xxxyzz = cbuffer.data(fi_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xzz_xxxzzz = cbuffer.data(fi_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyyy = cbuffer.data(fi_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyyz = cbuffer.data(fi_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xzz_xxyyzz = cbuffer.data(fi_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xzz_xxyzzz = cbuffer.data(fi_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xzz_xxzzzz = cbuffer.data(fi_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyyy = cbuffer.data(fi_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyyz = cbuffer.data(fi_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xzz_xyyyzz = cbuffer.data(fi_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xzz_xyyzzz = cbuffer.data(fi_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xzz_xyzzzz = cbuffer.data(fi_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xzz_xzzzzz = cbuffer.data(fi_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyyy = cbuffer.data(fi_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyyz = cbuffer.data(fi_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xzz_yyyyzz = cbuffer.data(fi_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xzz_yyyzzz = cbuffer.data(fi_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xzz_yyzzzz = cbuffer.data(fi_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xzz_yzzzzz = cbuffer.data(fi_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xzz_zzzzzz = cbuffer.data(fi_geom_10_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xz_xxxxxx, g_x_0_xz_xxxxxxz, g_x_0_xz_xxxxxy, g_x_0_xz_xxxxxyz, g_x_0_xz_xxxxxz, g_x_0_xz_xxxxxzz, g_x_0_xz_xxxxyy, g_x_0_xz_xxxxyyz, g_x_0_xz_xxxxyz, g_x_0_xz_xxxxyzz, g_x_0_xz_xxxxzz, g_x_0_xz_xxxxzzz, g_x_0_xz_xxxyyy, g_x_0_xz_xxxyyyz, g_x_0_xz_xxxyyz, g_x_0_xz_xxxyyzz, g_x_0_xz_xxxyzz, g_x_0_xz_xxxyzzz, g_x_0_xz_xxxzzz, g_x_0_xz_xxxzzzz, g_x_0_xz_xxyyyy, g_x_0_xz_xxyyyyz, g_x_0_xz_xxyyyz, g_x_0_xz_xxyyyzz, g_x_0_xz_xxyyzz, g_x_0_xz_xxyyzzz, g_x_0_xz_xxyzzz, g_x_0_xz_xxyzzzz, g_x_0_xz_xxzzzz, g_x_0_xz_xxzzzzz, g_x_0_xz_xyyyyy, g_x_0_xz_xyyyyyz, g_x_0_xz_xyyyyz, g_x_0_xz_xyyyyzz, g_x_0_xz_xyyyzz, g_x_0_xz_xyyyzzz, g_x_0_xz_xyyzzz, g_x_0_xz_xyyzzzz, g_x_0_xz_xyzzzz, g_x_0_xz_xyzzzzz, g_x_0_xz_xzzzzz, g_x_0_xz_xzzzzzz, g_x_0_xz_yyyyyy, g_x_0_xz_yyyyyyz, g_x_0_xz_yyyyyz, g_x_0_xz_yyyyyzz, g_x_0_xz_yyyyzz, g_x_0_xz_yyyyzzz, g_x_0_xz_yyyzzz, g_x_0_xz_yyyzzzz, g_x_0_xz_yyzzzz, g_x_0_xz_yyzzzzz, g_x_0_xz_yzzzzz, g_x_0_xz_yzzzzzz, g_x_0_xz_zzzzzz, g_x_0_xz_zzzzzzz, g_x_0_xzz_xxxxxx, g_x_0_xzz_xxxxxy, g_x_0_xzz_xxxxxz, g_x_0_xzz_xxxxyy, g_x_0_xzz_xxxxyz, g_x_0_xzz_xxxxzz, g_x_0_xzz_xxxyyy, g_x_0_xzz_xxxyyz, g_x_0_xzz_xxxyzz, g_x_0_xzz_xxxzzz, g_x_0_xzz_xxyyyy, g_x_0_xzz_xxyyyz, g_x_0_xzz_xxyyzz, g_x_0_xzz_xxyzzz, g_x_0_xzz_xxzzzz, g_x_0_xzz_xyyyyy, g_x_0_xzz_xyyyyz, g_x_0_xzz_xyyyzz, g_x_0_xzz_xyyzzz, g_x_0_xzz_xyzzzz, g_x_0_xzz_xzzzzz, g_x_0_xzz_yyyyyy, g_x_0_xzz_yyyyyz, g_x_0_xzz_yyyyzz, g_x_0_xzz_yyyzzz, g_x_0_xzz_yyzzzz, g_x_0_xzz_yzzzzz, g_x_0_xzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzz_xxxxxx[k] = -g_x_0_xz_xxxxxx[k] * ab_z + g_x_0_xz_xxxxxxz[k];

                g_x_0_xzz_xxxxxy[k] = -g_x_0_xz_xxxxxy[k] * ab_z + g_x_0_xz_xxxxxyz[k];

                g_x_0_xzz_xxxxxz[k] = -g_x_0_xz_xxxxxz[k] * ab_z + g_x_0_xz_xxxxxzz[k];

                g_x_0_xzz_xxxxyy[k] = -g_x_0_xz_xxxxyy[k] * ab_z + g_x_0_xz_xxxxyyz[k];

                g_x_0_xzz_xxxxyz[k] = -g_x_0_xz_xxxxyz[k] * ab_z + g_x_0_xz_xxxxyzz[k];

                g_x_0_xzz_xxxxzz[k] = -g_x_0_xz_xxxxzz[k] * ab_z + g_x_0_xz_xxxxzzz[k];

                g_x_0_xzz_xxxyyy[k] = -g_x_0_xz_xxxyyy[k] * ab_z + g_x_0_xz_xxxyyyz[k];

                g_x_0_xzz_xxxyyz[k] = -g_x_0_xz_xxxyyz[k] * ab_z + g_x_0_xz_xxxyyzz[k];

                g_x_0_xzz_xxxyzz[k] = -g_x_0_xz_xxxyzz[k] * ab_z + g_x_0_xz_xxxyzzz[k];

                g_x_0_xzz_xxxzzz[k] = -g_x_0_xz_xxxzzz[k] * ab_z + g_x_0_xz_xxxzzzz[k];

                g_x_0_xzz_xxyyyy[k] = -g_x_0_xz_xxyyyy[k] * ab_z + g_x_0_xz_xxyyyyz[k];

                g_x_0_xzz_xxyyyz[k] = -g_x_0_xz_xxyyyz[k] * ab_z + g_x_0_xz_xxyyyzz[k];

                g_x_0_xzz_xxyyzz[k] = -g_x_0_xz_xxyyzz[k] * ab_z + g_x_0_xz_xxyyzzz[k];

                g_x_0_xzz_xxyzzz[k] = -g_x_0_xz_xxyzzz[k] * ab_z + g_x_0_xz_xxyzzzz[k];

                g_x_0_xzz_xxzzzz[k] = -g_x_0_xz_xxzzzz[k] * ab_z + g_x_0_xz_xxzzzzz[k];

                g_x_0_xzz_xyyyyy[k] = -g_x_0_xz_xyyyyy[k] * ab_z + g_x_0_xz_xyyyyyz[k];

                g_x_0_xzz_xyyyyz[k] = -g_x_0_xz_xyyyyz[k] * ab_z + g_x_0_xz_xyyyyzz[k];

                g_x_0_xzz_xyyyzz[k] = -g_x_0_xz_xyyyzz[k] * ab_z + g_x_0_xz_xyyyzzz[k];

                g_x_0_xzz_xyyzzz[k] = -g_x_0_xz_xyyzzz[k] * ab_z + g_x_0_xz_xyyzzzz[k];

                g_x_0_xzz_xyzzzz[k] = -g_x_0_xz_xyzzzz[k] * ab_z + g_x_0_xz_xyzzzzz[k];

                g_x_0_xzz_xzzzzz[k] = -g_x_0_xz_xzzzzz[k] * ab_z + g_x_0_xz_xzzzzzz[k];

                g_x_0_xzz_yyyyyy[k] = -g_x_0_xz_yyyyyy[k] * ab_z + g_x_0_xz_yyyyyyz[k];

                g_x_0_xzz_yyyyyz[k] = -g_x_0_xz_yyyyyz[k] * ab_z + g_x_0_xz_yyyyyzz[k];

                g_x_0_xzz_yyyyzz[k] = -g_x_0_xz_yyyyzz[k] * ab_z + g_x_0_xz_yyyyzzz[k];

                g_x_0_xzz_yyyzzz[k] = -g_x_0_xz_yyyzzz[k] * ab_z + g_x_0_xz_yyyzzzz[k];

                g_x_0_xzz_yyzzzz[k] = -g_x_0_xz_yyzzzz[k] * ab_z + g_x_0_xz_yyzzzzz[k];

                g_x_0_xzz_yzzzzz[k] = -g_x_0_xz_yzzzzz[k] * ab_z + g_x_0_xz_yzzzzzz[k];

                g_x_0_xzz_zzzzzz[k] = -g_x_0_xz_zzzzzz[k] * ab_z + g_x_0_xz_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyy_xxxxxx = cbuffer.data(fi_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxy = cbuffer.data(fi_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxxz = cbuffer.data(fi_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxyy = cbuffer.data(fi_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxyz = cbuffer.data(fi_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_yyy_xxxxzz = cbuffer.data(fi_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyyy = cbuffer.data(fi_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyyz = cbuffer.data(fi_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_yyy_xxxyzz = cbuffer.data(fi_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_yyy_xxxzzz = cbuffer.data(fi_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyyy = cbuffer.data(fi_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyyz = cbuffer.data(fi_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_yyy_xxyyzz = cbuffer.data(fi_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_yyy_xxyzzz = cbuffer.data(fi_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_yyy_xxzzzz = cbuffer.data(fi_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyyy = cbuffer.data(fi_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyyz = cbuffer.data(fi_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_yyy_xyyyzz = cbuffer.data(fi_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_yyy_xyyzzz = cbuffer.data(fi_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_yyy_xyzzzz = cbuffer.data(fi_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_yyy_xzzzzz = cbuffer.data(fi_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyyy = cbuffer.data(fi_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyyz = cbuffer.data(fi_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_yyy_yyyyzz = cbuffer.data(fi_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_yyy_yyyzzz = cbuffer.data(fi_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_yyy_yyzzzz = cbuffer.data(fi_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_yyy_yzzzzz = cbuffer.data(fi_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_yyy_zzzzzz = cbuffer.data(fi_geom_10_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yy_xxxxxx, g_x_0_yy_xxxxxxy, g_x_0_yy_xxxxxy, g_x_0_yy_xxxxxyy, g_x_0_yy_xxxxxyz, g_x_0_yy_xxxxxz, g_x_0_yy_xxxxyy, g_x_0_yy_xxxxyyy, g_x_0_yy_xxxxyyz, g_x_0_yy_xxxxyz, g_x_0_yy_xxxxyzz, g_x_0_yy_xxxxzz, g_x_0_yy_xxxyyy, g_x_0_yy_xxxyyyy, g_x_0_yy_xxxyyyz, g_x_0_yy_xxxyyz, g_x_0_yy_xxxyyzz, g_x_0_yy_xxxyzz, g_x_0_yy_xxxyzzz, g_x_0_yy_xxxzzz, g_x_0_yy_xxyyyy, g_x_0_yy_xxyyyyy, g_x_0_yy_xxyyyyz, g_x_0_yy_xxyyyz, g_x_0_yy_xxyyyzz, g_x_0_yy_xxyyzz, g_x_0_yy_xxyyzzz, g_x_0_yy_xxyzzz, g_x_0_yy_xxyzzzz, g_x_0_yy_xxzzzz, g_x_0_yy_xyyyyy, g_x_0_yy_xyyyyyy, g_x_0_yy_xyyyyyz, g_x_0_yy_xyyyyz, g_x_0_yy_xyyyyzz, g_x_0_yy_xyyyzz, g_x_0_yy_xyyyzzz, g_x_0_yy_xyyzzz, g_x_0_yy_xyyzzzz, g_x_0_yy_xyzzzz, g_x_0_yy_xyzzzzz, g_x_0_yy_xzzzzz, g_x_0_yy_yyyyyy, g_x_0_yy_yyyyyyy, g_x_0_yy_yyyyyyz, g_x_0_yy_yyyyyz, g_x_0_yy_yyyyyzz, g_x_0_yy_yyyyzz, g_x_0_yy_yyyyzzz, g_x_0_yy_yyyzzz, g_x_0_yy_yyyzzzz, g_x_0_yy_yyzzzz, g_x_0_yy_yyzzzzz, g_x_0_yy_yzzzzz, g_x_0_yy_yzzzzzz, g_x_0_yy_zzzzzz, g_x_0_yyy_xxxxxx, g_x_0_yyy_xxxxxy, g_x_0_yyy_xxxxxz, g_x_0_yyy_xxxxyy, g_x_0_yyy_xxxxyz, g_x_0_yyy_xxxxzz, g_x_0_yyy_xxxyyy, g_x_0_yyy_xxxyyz, g_x_0_yyy_xxxyzz, g_x_0_yyy_xxxzzz, g_x_0_yyy_xxyyyy, g_x_0_yyy_xxyyyz, g_x_0_yyy_xxyyzz, g_x_0_yyy_xxyzzz, g_x_0_yyy_xxzzzz, g_x_0_yyy_xyyyyy, g_x_0_yyy_xyyyyz, g_x_0_yyy_xyyyzz, g_x_0_yyy_xyyzzz, g_x_0_yyy_xyzzzz, g_x_0_yyy_xzzzzz, g_x_0_yyy_yyyyyy, g_x_0_yyy_yyyyyz, g_x_0_yyy_yyyyzz, g_x_0_yyy_yyyzzz, g_x_0_yyy_yyzzzz, g_x_0_yyy_yzzzzz, g_x_0_yyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyy_xxxxxx[k] = -g_x_0_yy_xxxxxx[k] * ab_y + g_x_0_yy_xxxxxxy[k];

                g_x_0_yyy_xxxxxy[k] = -g_x_0_yy_xxxxxy[k] * ab_y + g_x_0_yy_xxxxxyy[k];

                g_x_0_yyy_xxxxxz[k] = -g_x_0_yy_xxxxxz[k] * ab_y + g_x_0_yy_xxxxxyz[k];

                g_x_0_yyy_xxxxyy[k] = -g_x_0_yy_xxxxyy[k] * ab_y + g_x_0_yy_xxxxyyy[k];

                g_x_0_yyy_xxxxyz[k] = -g_x_0_yy_xxxxyz[k] * ab_y + g_x_0_yy_xxxxyyz[k];

                g_x_0_yyy_xxxxzz[k] = -g_x_0_yy_xxxxzz[k] * ab_y + g_x_0_yy_xxxxyzz[k];

                g_x_0_yyy_xxxyyy[k] = -g_x_0_yy_xxxyyy[k] * ab_y + g_x_0_yy_xxxyyyy[k];

                g_x_0_yyy_xxxyyz[k] = -g_x_0_yy_xxxyyz[k] * ab_y + g_x_0_yy_xxxyyyz[k];

                g_x_0_yyy_xxxyzz[k] = -g_x_0_yy_xxxyzz[k] * ab_y + g_x_0_yy_xxxyyzz[k];

                g_x_0_yyy_xxxzzz[k] = -g_x_0_yy_xxxzzz[k] * ab_y + g_x_0_yy_xxxyzzz[k];

                g_x_0_yyy_xxyyyy[k] = -g_x_0_yy_xxyyyy[k] * ab_y + g_x_0_yy_xxyyyyy[k];

                g_x_0_yyy_xxyyyz[k] = -g_x_0_yy_xxyyyz[k] * ab_y + g_x_0_yy_xxyyyyz[k];

                g_x_0_yyy_xxyyzz[k] = -g_x_0_yy_xxyyzz[k] * ab_y + g_x_0_yy_xxyyyzz[k];

                g_x_0_yyy_xxyzzz[k] = -g_x_0_yy_xxyzzz[k] * ab_y + g_x_0_yy_xxyyzzz[k];

                g_x_0_yyy_xxzzzz[k] = -g_x_0_yy_xxzzzz[k] * ab_y + g_x_0_yy_xxyzzzz[k];

                g_x_0_yyy_xyyyyy[k] = -g_x_0_yy_xyyyyy[k] * ab_y + g_x_0_yy_xyyyyyy[k];

                g_x_0_yyy_xyyyyz[k] = -g_x_0_yy_xyyyyz[k] * ab_y + g_x_0_yy_xyyyyyz[k];

                g_x_0_yyy_xyyyzz[k] = -g_x_0_yy_xyyyzz[k] * ab_y + g_x_0_yy_xyyyyzz[k];

                g_x_0_yyy_xyyzzz[k] = -g_x_0_yy_xyyzzz[k] * ab_y + g_x_0_yy_xyyyzzz[k];

                g_x_0_yyy_xyzzzz[k] = -g_x_0_yy_xyzzzz[k] * ab_y + g_x_0_yy_xyyzzzz[k];

                g_x_0_yyy_xzzzzz[k] = -g_x_0_yy_xzzzzz[k] * ab_y + g_x_0_yy_xyzzzzz[k];

                g_x_0_yyy_yyyyyy[k] = -g_x_0_yy_yyyyyy[k] * ab_y + g_x_0_yy_yyyyyyy[k];

                g_x_0_yyy_yyyyyz[k] = -g_x_0_yy_yyyyyz[k] * ab_y + g_x_0_yy_yyyyyyz[k];

                g_x_0_yyy_yyyyzz[k] = -g_x_0_yy_yyyyzz[k] * ab_y + g_x_0_yy_yyyyyzz[k];

                g_x_0_yyy_yyyzzz[k] = -g_x_0_yy_yyyzzz[k] * ab_y + g_x_0_yy_yyyyzzz[k];

                g_x_0_yyy_yyzzzz[k] = -g_x_0_yy_yyzzzz[k] * ab_y + g_x_0_yy_yyyzzzz[k];

                g_x_0_yyy_yzzzzz[k] = -g_x_0_yy_yzzzzz[k] * ab_y + g_x_0_yy_yyzzzzz[k];

                g_x_0_yyy_zzzzzz[k] = -g_x_0_yy_zzzzzz[k] * ab_y + g_x_0_yy_yzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyz_xxxxxx = cbuffer.data(fi_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxy = cbuffer.data(fi_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxxz = cbuffer.data(fi_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxyy = cbuffer.data(fi_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxyz = cbuffer.data(fi_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_yyz_xxxxzz = cbuffer.data(fi_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyyy = cbuffer.data(fi_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyyz = cbuffer.data(fi_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_yyz_xxxyzz = cbuffer.data(fi_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_yyz_xxxzzz = cbuffer.data(fi_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyyy = cbuffer.data(fi_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyyz = cbuffer.data(fi_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_yyz_xxyyzz = cbuffer.data(fi_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_yyz_xxyzzz = cbuffer.data(fi_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_yyz_xxzzzz = cbuffer.data(fi_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyyy = cbuffer.data(fi_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyyz = cbuffer.data(fi_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_yyz_xyyyzz = cbuffer.data(fi_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_yyz_xyyzzz = cbuffer.data(fi_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_yyz_xyzzzz = cbuffer.data(fi_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_yyz_xzzzzz = cbuffer.data(fi_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyyy = cbuffer.data(fi_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyyz = cbuffer.data(fi_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_yyz_yyyyzz = cbuffer.data(fi_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_yyz_yyyzzz = cbuffer.data(fi_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_yyz_yyzzzz = cbuffer.data(fi_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_yyz_yzzzzz = cbuffer.data(fi_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_yyz_zzzzzz = cbuffer.data(fi_geom_10_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyz_xxxxxx, g_x_0_yyz_xxxxxy, g_x_0_yyz_xxxxxz, g_x_0_yyz_xxxxyy, g_x_0_yyz_xxxxyz, g_x_0_yyz_xxxxzz, g_x_0_yyz_xxxyyy, g_x_0_yyz_xxxyyz, g_x_0_yyz_xxxyzz, g_x_0_yyz_xxxzzz, g_x_0_yyz_xxyyyy, g_x_0_yyz_xxyyyz, g_x_0_yyz_xxyyzz, g_x_0_yyz_xxyzzz, g_x_0_yyz_xxzzzz, g_x_0_yyz_xyyyyy, g_x_0_yyz_xyyyyz, g_x_0_yyz_xyyyzz, g_x_0_yyz_xyyzzz, g_x_0_yyz_xyzzzz, g_x_0_yyz_xzzzzz, g_x_0_yyz_yyyyyy, g_x_0_yyz_yyyyyz, g_x_0_yyz_yyyyzz, g_x_0_yyz_yyyzzz, g_x_0_yyz_yyzzzz, g_x_0_yyz_yzzzzz, g_x_0_yyz_zzzzzz, g_x_0_yz_xxxxxx, g_x_0_yz_xxxxxxy, g_x_0_yz_xxxxxy, g_x_0_yz_xxxxxyy, g_x_0_yz_xxxxxyz, g_x_0_yz_xxxxxz, g_x_0_yz_xxxxyy, g_x_0_yz_xxxxyyy, g_x_0_yz_xxxxyyz, g_x_0_yz_xxxxyz, g_x_0_yz_xxxxyzz, g_x_0_yz_xxxxzz, g_x_0_yz_xxxyyy, g_x_0_yz_xxxyyyy, g_x_0_yz_xxxyyyz, g_x_0_yz_xxxyyz, g_x_0_yz_xxxyyzz, g_x_0_yz_xxxyzz, g_x_0_yz_xxxyzzz, g_x_0_yz_xxxzzz, g_x_0_yz_xxyyyy, g_x_0_yz_xxyyyyy, g_x_0_yz_xxyyyyz, g_x_0_yz_xxyyyz, g_x_0_yz_xxyyyzz, g_x_0_yz_xxyyzz, g_x_0_yz_xxyyzzz, g_x_0_yz_xxyzzz, g_x_0_yz_xxyzzzz, g_x_0_yz_xxzzzz, g_x_0_yz_xyyyyy, g_x_0_yz_xyyyyyy, g_x_0_yz_xyyyyyz, g_x_0_yz_xyyyyz, g_x_0_yz_xyyyyzz, g_x_0_yz_xyyyzz, g_x_0_yz_xyyyzzz, g_x_0_yz_xyyzzz, g_x_0_yz_xyyzzzz, g_x_0_yz_xyzzzz, g_x_0_yz_xyzzzzz, g_x_0_yz_xzzzzz, g_x_0_yz_yyyyyy, g_x_0_yz_yyyyyyy, g_x_0_yz_yyyyyyz, g_x_0_yz_yyyyyz, g_x_0_yz_yyyyyzz, g_x_0_yz_yyyyzz, g_x_0_yz_yyyyzzz, g_x_0_yz_yyyzzz, g_x_0_yz_yyyzzzz, g_x_0_yz_yyzzzz, g_x_0_yz_yyzzzzz, g_x_0_yz_yzzzzz, g_x_0_yz_yzzzzzz, g_x_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyz_xxxxxx[k] = -g_x_0_yz_xxxxxx[k] * ab_y + g_x_0_yz_xxxxxxy[k];

                g_x_0_yyz_xxxxxy[k] = -g_x_0_yz_xxxxxy[k] * ab_y + g_x_0_yz_xxxxxyy[k];

                g_x_0_yyz_xxxxxz[k] = -g_x_0_yz_xxxxxz[k] * ab_y + g_x_0_yz_xxxxxyz[k];

                g_x_0_yyz_xxxxyy[k] = -g_x_0_yz_xxxxyy[k] * ab_y + g_x_0_yz_xxxxyyy[k];

                g_x_0_yyz_xxxxyz[k] = -g_x_0_yz_xxxxyz[k] * ab_y + g_x_0_yz_xxxxyyz[k];

                g_x_0_yyz_xxxxzz[k] = -g_x_0_yz_xxxxzz[k] * ab_y + g_x_0_yz_xxxxyzz[k];

                g_x_0_yyz_xxxyyy[k] = -g_x_0_yz_xxxyyy[k] * ab_y + g_x_0_yz_xxxyyyy[k];

                g_x_0_yyz_xxxyyz[k] = -g_x_0_yz_xxxyyz[k] * ab_y + g_x_0_yz_xxxyyyz[k];

                g_x_0_yyz_xxxyzz[k] = -g_x_0_yz_xxxyzz[k] * ab_y + g_x_0_yz_xxxyyzz[k];

                g_x_0_yyz_xxxzzz[k] = -g_x_0_yz_xxxzzz[k] * ab_y + g_x_0_yz_xxxyzzz[k];

                g_x_0_yyz_xxyyyy[k] = -g_x_0_yz_xxyyyy[k] * ab_y + g_x_0_yz_xxyyyyy[k];

                g_x_0_yyz_xxyyyz[k] = -g_x_0_yz_xxyyyz[k] * ab_y + g_x_0_yz_xxyyyyz[k];

                g_x_0_yyz_xxyyzz[k] = -g_x_0_yz_xxyyzz[k] * ab_y + g_x_0_yz_xxyyyzz[k];

                g_x_0_yyz_xxyzzz[k] = -g_x_0_yz_xxyzzz[k] * ab_y + g_x_0_yz_xxyyzzz[k];

                g_x_0_yyz_xxzzzz[k] = -g_x_0_yz_xxzzzz[k] * ab_y + g_x_0_yz_xxyzzzz[k];

                g_x_0_yyz_xyyyyy[k] = -g_x_0_yz_xyyyyy[k] * ab_y + g_x_0_yz_xyyyyyy[k];

                g_x_0_yyz_xyyyyz[k] = -g_x_0_yz_xyyyyz[k] * ab_y + g_x_0_yz_xyyyyyz[k];

                g_x_0_yyz_xyyyzz[k] = -g_x_0_yz_xyyyzz[k] * ab_y + g_x_0_yz_xyyyyzz[k];

                g_x_0_yyz_xyyzzz[k] = -g_x_0_yz_xyyzzz[k] * ab_y + g_x_0_yz_xyyyzzz[k];

                g_x_0_yyz_xyzzzz[k] = -g_x_0_yz_xyzzzz[k] * ab_y + g_x_0_yz_xyyzzzz[k];

                g_x_0_yyz_xzzzzz[k] = -g_x_0_yz_xzzzzz[k] * ab_y + g_x_0_yz_xyzzzzz[k];

                g_x_0_yyz_yyyyyy[k] = -g_x_0_yz_yyyyyy[k] * ab_y + g_x_0_yz_yyyyyyy[k];

                g_x_0_yyz_yyyyyz[k] = -g_x_0_yz_yyyyyz[k] * ab_y + g_x_0_yz_yyyyyyz[k];

                g_x_0_yyz_yyyyzz[k] = -g_x_0_yz_yyyyzz[k] * ab_y + g_x_0_yz_yyyyyzz[k];

                g_x_0_yyz_yyyzzz[k] = -g_x_0_yz_yyyzzz[k] * ab_y + g_x_0_yz_yyyyzzz[k];

                g_x_0_yyz_yyzzzz[k] = -g_x_0_yz_yyzzzz[k] * ab_y + g_x_0_yz_yyyzzzz[k];

                g_x_0_yyz_yzzzzz[k] = -g_x_0_yz_yzzzzz[k] * ab_y + g_x_0_yz_yyzzzzz[k];

                g_x_0_yyz_zzzzzz[k] = -g_x_0_yz_zzzzzz[k] * ab_y + g_x_0_yz_yzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzz_xxxxxx = cbuffer.data(fi_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxy = cbuffer.data(fi_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxxz = cbuffer.data(fi_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxyy = cbuffer.data(fi_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxyz = cbuffer.data(fi_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_yzz_xxxxzz = cbuffer.data(fi_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyyy = cbuffer.data(fi_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyyz = cbuffer.data(fi_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_yzz_xxxyzz = cbuffer.data(fi_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_yzz_xxxzzz = cbuffer.data(fi_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyyy = cbuffer.data(fi_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyyz = cbuffer.data(fi_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_yzz_xxyyzz = cbuffer.data(fi_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_yzz_xxyzzz = cbuffer.data(fi_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_yzz_xxzzzz = cbuffer.data(fi_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyyy = cbuffer.data(fi_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyyz = cbuffer.data(fi_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_yzz_xyyyzz = cbuffer.data(fi_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_yzz_xyyzzz = cbuffer.data(fi_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_yzz_xyzzzz = cbuffer.data(fi_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_yzz_xzzzzz = cbuffer.data(fi_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyyy = cbuffer.data(fi_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyyz = cbuffer.data(fi_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_yzz_yyyyzz = cbuffer.data(fi_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_yzz_yyyzzz = cbuffer.data(fi_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_yzz_yyzzzz = cbuffer.data(fi_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_yzz_yzzzzz = cbuffer.data(fi_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_yzz_zzzzzz = cbuffer.data(fi_geom_10_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzz_xxxxxx, g_x_0_yzz_xxxxxy, g_x_0_yzz_xxxxxz, g_x_0_yzz_xxxxyy, g_x_0_yzz_xxxxyz, g_x_0_yzz_xxxxzz, g_x_0_yzz_xxxyyy, g_x_0_yzz_xxxyyz, g_x_0_yzz_xxxyzz, g_x_0_yzz_xxxzzz, g_x_0_yzz_xxyyyy, g_x_0_yzz_xxyyyz, g_x_0_yzz_xxyyzz, g_x_0_yzz_xxyzzz, g_x_0_yzz_xxzzzz, g_x_0_yzz_xyyyyy, g_x_0_yzz_xyyyyz, g_x_0_yzz_xyyyzz, g_x_0_yzz_xyyzzz, g_x_0_yzz_xyzzzz, g_x_0_yzz_xzzzzz, g_x_0_yzz_yyyyyy, g_x_0_yzz_yyyyyz, g_x_0_yzz_yyyyzz, g_x_0_yzz_yyyzzz, g_x_0_yzz_yyzzzz, g_x_0_yzz_yzzzzz, g_x_0_yzz_zzzzzz, g_x_0_zz_xxxxxx, g_x_0_zz_xxxxxxy, g_x_0_zz_xxxxxy, g_x_0_zz_xxxxxyy, g_x_0_zz_xxxxxyz, g_x_0_zz_xxxxxz, g_x_0_zz_xxxxyy, g_x_0_zz_xxxxyyy, g_x_0_zz_xxxxyyz, g_x_0_zz_xxxxyz, g_x_0_zz_xxxxyzz, g_x_0_zz_xxxxzz, g_x_0_zz_xxxyyy, g_x_0_zz_xxxyyyy, g_x_0_zz_xxxyyyz, g_x_0_zz_xxxyyz, g_x_0_zz_xxxyyzz, g_x_0_zz_xxxyzz, g_x_0_zz_xxxyzzz, g_x_0_zz_xxxzzz, g_x_0_zz_xxyyyy, g_x_0_zz_xxyyyyy, g_x_0_zz_xxyyyyz, g_x_0_zz_xxyyyz, g_x_0_zz_xxyyyzz, g_x_0_zz_xxyyzz, g_x_0_zz_xxyyzzz, g_x_0_zz_xxyzzz, g_x_0_zz_xxyzzzz, g_x_0_zz_xxzzzz, g_x_0_zz_xyyyyy, g_x_0_zz_xyyyyyy, g_x_0_zz_xyyyyyz, g_x_0_zz_xyyyyz, g_x_0_zz_xyyyyzz, g_x_0_zz_xyyyzz, g_x_0_zz_xyyyzzz, g_x_0_zz_xyyzzz, g_x_0_zz_xyyzzzz, g_x_0_zz_xyzzzz, g_x_0_zz_xyzzzzz, g_x_0_zz_xzzzzz, g_x_0_zz_yyyyyy, g_x_0_zz_yyyyyyy, g_x_0_zz_yyyyyyz, g_x_0_zz_yyyyyz, g_x_0_zz_yyyyyzz, g_x_0_zz_yyyyzz, g_x_0_zz_yyyyzzz, g_x_0_zz_yyyzzz, g_x_0_zz_yyyzzzz, g_x_0_zz_yyzzzz, g_x_0_zz_yyzzzzz, g_x_0_zz_yzzzzz, g_x_0_zz_yzzzzzz, g_x_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzz_xxxxxx[k] = -g_x_0_zz_xxxxxx[k] * ab_y + g_x_0_zz_xxxxxxy[k];

                g_x_0_yzz_xxxxxy[k] = -g_x_0_zz_xxxxxy[k] * ab_y + g_x_0_zz_xxxxxyy[k];

                g_x_0_yzz_xxxxxz[k] = -g_x_0_zz_xxxxxz[k] * ab_y + g_x_0_zz_xxxxxyz[k];

                g_x_0_yzz_xxxxyy[k] = -g_x_0_zz_xxxxyy[k] * ab_y + g_x_0_zz_xxxxyyy[k];

                g_x_0_yzz_xxxxyz[k] = -g_x_0_zz_xxxxyz[k] * ab_y + g_x_0_zz_xxxxyyz[k];

                g_x_0_yzz_xxxxzz[k] = -g_x_0_zz_xxxxzz[k] * ab_y + g_x_0_zz_xxxxyzz[k];

                g_x_0_yzz_xxxyyy[k] = -g_x_0_zz_xxxyyy[k] * ab_y + g_x_0_zz_xxxyyyy[k];

                g_x_0_yzz_xxxyyz[k] = -g_x_0_zz_xxxyyz[k] * ab_y + g_x_0_zz_xxxyyyz[k];

                g_x_0_yzz_xxxyzz[k] = -g_x_0_zz_xxxyzz[k] * ab_y + g_x_0_zz_xxxyyzz[k];

                g_x_0_yzz_xxxzzz[k] = -g_x_0_zz_xxxzzz[k] * ab_y + g_x_0_zz_xxxyzzz[k];

                g_x_0_yzz_xxyyyy[k] = -g_x_0_zz_xxyyyy[k] * ab_y + g_x_0_zz_xxyyyyy[k];

                g_x_0_yzz_xxyyyz[k] = -g_x_0_zz_xxyyyz[k] * ab_y + g_x_0_zz_xxyyyyz[k];

                g_x_0_yzz_xxyyzz[k] = -g_x_0_zz_xxyyzz[k] * ab_y + g_x_0_zz_xxyyyzz[k];

                g_x_0_yzz_xxyzzz[k] = -g_x_0_zz_xxyzzz[k] * ab_y + g_x_0_zz_xxyyzzz[k];

                g_x_0_yzz_xxzzzz[k] = -g_x_0_zz_xxzzzz[k] * ab_y + g_x_0_zz_xxyzzzz[k];

                g_x_0_yzz_xyyyyy[k] = -g_x_0_zz_xyyyyy[k] * ab_y + g_x_0_zz_xyyyyyy[k];

                g_x_0_yzz_xyyyyz[k] = -g_x_0_zz_xyyyyz[k] * ab_y + g_x_0_zz_xyyyyyz[k];

                g_x_0_yzz_xyyyzz[k] = -g_x_0_zz_xyyyzz[k] * ab_y + g_x_0_zz_xyyyyzz[k];

                g_x_0_yzz_xyyzzz[k] = -g_x_0_zz_xyyzzz[k] * ab_y + g_x_0_zz_xyyyzzz[k];

                g_x_0_yzz_xyzzzz[k] = -g_x_0_zz_xyzzzz[k] * ab_y + g_x_0_zz_xyyzzzz[k];

                g_x_0_yzz_xzzzzz[k] = -g_x_0_zz_xzzzzz[k] * ab_y + g_x_0_zz_xyzzzzz[k];

                g_x_0_yzz_yyyyyy[k] = -g_x_0_zz_yyyyyy[k] * ab_y + g_x_0_zz_yyyyyyy[k];

                g_x_0_yzz_yyyyyz[k] = -g_x_0_zz_yyyyyz[k] * ab_y + g_x_0_zz_yyyyyyz[k];

                g_x_0_yzz_yyyyzz[k] = -g_x_0_zz_yyyyzz[k] * ab_y + g_x_0_zz_yyyyyzz[k];

                g_x_0_yzz_yyyzzz[k] = -g_x_0_zz_yyyzzz[k] * ab_y + g_x_0_zz_yyyyzzz[k];

                g_x_0_yzz_yyzzzz[k] = -g_x_0_zz_yyzzzz[k] * ab_y + g_x_0_zz_yyyzzzz[k];

                g_x_0_yzz_yzzzzz[k] = -g_x_0_zz_yzzzzz[k] * ab_y + g_x_0_zz_yyzzzzz[k];

                g_x_0_yzz_zzzzzz[k] = -g_x_0_zz_zzzzzz[k] * ab_y + g_x_0_zz_yzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzz_xxxxxx = cbuffer.data(fi_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxy = cbuffer.data(fi_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxxz = cbuffer.data(fi_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxyy = cbuffer.data(fi_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxyz = cbuffer.data(fi_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_zzz_xxxxzz = cbuffer.data(fi_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyyy = cbuffer.data(fi_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyyz = cbuffer.data(fi_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_zzz_xxxyzz = cbuffer.data(fi_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_zzz_xxxzzz = cbuffer.data(fi_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyyy = cbuffer.data(fi_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyyz = cbuffer.data(fi_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_zzz_xxyyzz = cbuffer.data(fi_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_zzz_xxyzzz = cbuffer.data(fi_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_zzz_xxzzzz = cbuffer.data(fi_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyyy = cbuffer.data(fi_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyyz = cbuffer.data(fi_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_zzz_xyyyzz = cbuffer.data(fi_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_zzz_xyyzzz = cbuffer.data(fi_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_zzz_xyzzzz = cbuffer.data(fi_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_zzz_xzzzzz = cbuffer.data(fi_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyyy = cbuffer.data(fi_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyyz = cbuffer.data(fi_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_zzz_yyyyzz = cbuffer.data(fi_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_zzz_yyyzzz = cbuffer.data(fi_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_zzz_yyzzzz = cbuffer.data(fi_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_zzz_yzzzzz = cbuffer.data(fi_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_zzz_zzzzzz = cbuffer.data(fi_geom_10_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zz_xxxxxx, g_x_0_zz_xxxxxxz, g_x_0_zz_xxxxxy, g_x_0_zz_xxxxxyz, g_x_0_zz_xxxxxz, g_x_0_zz_xxxxxzz, g_x_0_zz_xxxxyy, g_x_0_zz_xxxxyyz, g_x_0_zz_xxxxyz, g_x_0_zz_xxxxyzz, g_x_0_zz_xxxxzz, g_x_0_zz_xxxxzzz, g_x_0_zz_xxxyyy, g_x_0_zz_xxxyyyz, g_x_0_zz_xxxyyz, g_x_0_zz_xxxyyzz, g_x_0_zz_xxxyzz, g_x_0_zz_xxxyzzz, g_x_0_zz_xxxzzz, g_x_0_zz_xxxzzzz, g_x_0_zz_xxyyyy, g_x_0_zz_xxyyyyz, g_x_0_zz_xxyyyz, g_x_0_zz_xxyyyzz, g_x_0_zz_xxyyzz, g_x_0_zz_xxyyzzz, g_x_0_zz_xxyzzz, g_x_0_zz_xxyzzzz, g_x_0_zz_xxzzzz, g_x_0_zz_xxzzzzz, g_x_0_zz_xyyyyy, g_x_0_zz_xyyyyyz, g_x_0_zz_xyyyyz, g_x_0_zz_xyyyyzz, g_x_0_zz_xyyyzz, g_x_0_zz_xyyyzzz, g_x_0_zz_xyyzzz, g_x_0_zz_xyyzzzz, g_x_0_zz_xyzzzz, g_x_0_zz_xyzzzzz, g_x_0_zz_xzzzzz, g_x_0_zz_xzzzzzz, g_x_0_zz_yyyyyy, g_x_0_zz_yyyyyyz, g_x_0_zz_yyyyyz, g_x_0_zz_yyyyyzz, g_x_0_zz_yyyyzz, g_x_0_zz_yyyyzzz, g_x_0_zz_yyyzzz, g_x_0_zz_yyyzzzz, g_x_0_zz_yyzzzz, g_x_0_zz_yyzzzzz, g_x_0_zz_yzzzzz, g_x_0_zz_yzzzzzz, g_x_0_zz_zzzzzz, g_x_0_zz_zzzzzzz, g_x_0_zzz_xxxxxx, g_x_0_zzz_xxxxxy, g_x_0_zzz_xxxxxz, g_x_0_zzz_xxxxyy, g_x_0_zzz_xxxxyz, g_x_0_zzz_xxxxzz, g_x_0_zzz_xxxyyy, g_x_0_zzz_xxxyyz, g_x_0_zzz_xxxyzz, g_x_0_zzz_xxxzzz, g_x_0_zzz_xxyyyy, g_x_0_zzz_xxyyyz, g_x_0_zzz_xxyyzz, g_x_0_zzz_xxyzzz, g_x_0_zzz_xxzzzz, g_x_0_zzz_xyyyyy, g_x_0_zzz_xyyyyz, g_x_0_zzz_xyyyzz, g_x_0_zzz_xyyzzz, g_x_0_zzz_xyzzzz, g_x_0_zzz_xzzzzz, g_x_0_zzz_yyyyyy, g_x_0_zzz_yyyyyz, g_x_0_zzz_yyyyzz, g_x_0_zzz_yyyzzz, g_x_0_zzz_yyzzzz, g_x_0_zzz_yzzzzz, g_x_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzz_xxxxxx[k] = -g_x_0_zz_xxxxxx[k] * ab_z + g_x_0_zz_xxxxxxz[k];

                g_x_0_zzz_xxxxxy[k] = -g_x_0_zz_xxxxxy[k] * ab_z + g_x_0_zz_xxxxxyz[k];

                g_x_0_zzz_xxxxxz[k] = -g_x_0_zz_xxxxxz[k] * ab_z + g_x_0_zz_xxxxxzz[k];

                g_x_0_zzz_xxxxyy[k] = -g_x_0_zz_xxxxyy[k] * ab_z + g_x_0_zz_xxxxyyz[k];

                g_x_0_zzz_xxxxyz[k] = -g_x_0_zz_xxxxyz[k] * ab_z + g_x_0_zz_xxxxyzz[k];

                g_x_0_zzz_xxxxzz[k] = -g_x_0_zz_xxxxzz[k] * ab_z + g_x_0_zz_xxxxzzz[k];

                g_x_0_zzz_xxxyyy[k] = -g_x_0_zz_xxxyyy[k] * ab_z + g_x_0_zz_xxxyyyz[k];

                g_x_0_zzz_xxxyyz[k] = -g_x_0_zz_xxxyyz[k] * ab_z + g_x_0_zz_xxxyyzz[k];

                g_x_0_zzz_xxxyzz[k] = -g_x_0_zz_xxxyzz[k] * ab_z + g_x_0_zz_xxxyzzz[k];

                g_x_0_zzz_xxxzzz[k] = -g_x_0_zz_xxxzzz[k] * ab_z + g_x_0_zz_xxxzzzz[k];

                g_x_0_zzz_xxyyyy[k] = -g_x_0_zz_xxyyyy[k] * ab_z + g_x_0_zz_xxyyyyz[k];

                g_x_0_zzz_xxyyyz[k] = -g_x_0_zz_xxyyyz[k] * ab_z + g_x_0_zz_xxyyyzz[k];

                g_x_0_zzz_xxyyzz[k] = -g_x_0_zz_xxyyzz[k] * ab_z + g_x_0_zz_xxyyzzz[k];

                g_x_0_zzz_xxyzzz[k] = -g_x_0_zz_xxyzzz[k] * ab_z + g_x_0_zz_xxyzzzz[k];

                g_x_0_zzz_xxzzzz[k] = -g_x_0_zz_xxzzzz[k] * ab_z + g_x_0_zz_xxzzzzz[k];

                g_x_0_zzz_xyyyyy[k] = -g_x_0_zz_xyyyyy[k] * ab_z + g_x_0_zz_xyyyyyz[k];

                g_x_0_zzz_xyyyyz[k] = -g_x_0_zz_xyyyyz[k] * ab_z + g_x_0_zz_xyyyyzz[k];

                g_x_0_zzz_xyyyzz[k] = -g_x_0_zz_xyyyzz[k] * ab_z + g_x_0_zz_xyyyzzz[k];

                g_x_0_zzz_xyyzzz[k] = -g_x_0_zz_xyyzzz[k] * ab_z + g_x_0_zz_xyyzzzz[k];

                g_x_0_zzz_xyzzzz[k] = -g_x_0_zz_xyzzzz[k] * ab_z + g_x_0_zz_xyzzzzz[k];

                g_x_0_zzz_xzzzzz[k] = -g_x_0_zz_xzzzzz[k] * ab_z + g_x_0_zz_xzzzzzz[k];

                g_x_0_zzz_yyyyyy[k] = -g_x_0_zz_yyyyyy[k] * ab_z + g_x_0_zz_yyyyyyz[k];

                g_x_0_zzz_yyyyyz[k] = -g_x_0_zz_yyyyyz[k] * ab_z + g_x_0_zz_yyyyyzz[k];

                g_x_0_zzz_yyyyzz[k] = -g_x_0_zz_yyyyzz[k] * ab_z + g_x_0_zz_yyyyzzz[k];

                g_x_0_zzz_yyyzzz[k] = -g_x_0_zz_yyyzzz[k] * ab_z + g_x_0_zz_yyyzzzz[k];

                g_x_0_zzz_yyzzzz[k] = -g_x_0_zz_yyzzzz[k] * ab_z + g_x_0_zz_yyzzzzz[k];

                g_x_0_zzz_yzzzzz[k] = -g_x_0_zz_yzzzzz[k] * ab_z + g_x_0_zz_yzzzzzz[k];

                g_x_0_zzz_zzzzzz[k] = -g_x_0_zz_zzzzzz[k] * ab_z + g_x_0_zz_zzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxx_xxxxxx = cbuffer.data(fi_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxy = cbuffer.data(fi_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxxz = cbuffer.data(fi_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxyy = cbuffer.data(fi_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxyz = cbuffer.data(fi_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_xxx_xxxxzz = cbuffer.data(fi_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyyy = cbuffer.data(fi_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyyz = cbuffer.data(fi_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_xxx_xxxyzz = cbuffer.data(fi_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_xxx_xxxzzz = cbuffer.data(fi_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyyy = cbuffer.data(fi_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyyz = cbuffer.data(fi_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_xxx_xxyyzz = cbuffer.data(fi_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_xxx_xxyzzz = cbuffer.data(fi_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_xxx_xxzzzz = cbuffer.data(fi_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyyy = cbuffer.data(fi_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyyz = cbuffer.data(fi_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_xxx_xyyyzz = cbuffer.data(fi_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_xxx_xyyzzz = cbuffer.data(fi_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_xxx_xyzzzz = cbuffer.data(fi_geom_10_off + 299 * ccomps * dcomps);

            auto g_y_0_xxx_xzzzzz = cbuffer.data(fi_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyyy = cbuffer.data(fi_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyyz = cbuffer.data(fi_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_xxx_yyyyzz = cbuffer.data(fi_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_xxx_yyyzzz = cbuffer.data(fi_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_xxx_yyzzzz = cbuffer.data(fi_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_xxx_yzzzzz = cbuffer.data(fi_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_xxx_zzzzzz = cbuffer.data(fi_geom_10_off + 307 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xx_xxxxxx, g_y_0_xx_xxxxxxx, g_y_0_xx_xxxxxxy, g_y_0_xx_xxxxxxz, g_y_0_xx_xxxxxy, g_y_0_xx_xxxxxyy, g_y_0_xx_xxxxxyz, g_y_0_xx_xxxxxz, g_y_0_xx_xxxxxzz, g_y_0_xx_xxxxyy, g_y_0_xx_xxxxyyy, g_y_0_xx_xxxxyyz, g_y_0_xx_xxxxyz, g_y_0_xx_xxxxyzz, g_y_0_xx_xxxxzz, g_y_0_xx_xxxxzzz, g_y_0_xx_xxxyyy, g_y_0_xx_xxxyyyy, g_y_0_xx_xxxyyyz, g_y_0_xx_xxxyyz, g_y_0_xx_xxxyyzz, g_y_0_xx_xxxyzz, g_y_0_xx_xxxyzzz, g_y_0_xx_xxxzzz, g_y_0_xx_xxxzzzz, g_y_0_xx_xxyyyy, g_y_0_xx_xxyyyyy, g_y_0_xx_xxyyyyz, g_y_0_xx_xxyyyz, g_y_0_xx_xxyyyzz, g_y_0_xx_xxyyzz, g_y_0_xx_xxyyzzz, g_y_0_xx_xxyzzz, g_y_0_xx_xxyzzzz, g_y_0_xx_xxzzzz, g_y_0_xx_xxzzzzz, g_y_0_xx_xyyyyy, g_y_0_xx_xyyyyyy, g_y_0_xx_xyyyyyz, g_y_0_xx_xyyyyz, g_y_0_xx_xyyyyzz, g_y_0_xx_xyyyzz, g_y_0_xx_xyyyzzz, g_y_0_xx_xyyzzz, g_y_0_xx_xyyzzzz, g_y_0_xx_xyzzzz, g_y_0_xx_xyzzzzz, g_y_0_xx_xzzzzz, g_y_0_xx_xzzzzzz, g_y_0_xx_yyyyyy, g_y_0_xx_yyyyyz, g_y_0_xx_yyyyzz, g_y_0_xx_yyyzzz, g_y_0_xx_yyzzzz, g_y_0_xx_yzzzzz, g_y_0_xx_zzzzzz, g_y_0_xxx_xxxxxx, g_y_0_xxx_xxxxxy, g_y_0_xxx_xxxxxz, g_y_0_xxx_xxxxyy, g_y_0_xxx_xxxxyz, g_y_0_xxx_xxxxzz, g_y_0_xxx_xxxyyy, g_y_0_xxx_xxxyyz, g_y_0_xxx_xxxyzz, g_y_0_xxx_xxxzzz, g_y_0_xxx_xxyyyy, g_y_0_xxx_xxyyyz, g_y_0_xxx_xxyyzz, g_y_0_xxx_xxyzzz, g_y_0_xxx_xxzzzz, g_y_0_xxx_xyyyyy, g_y_0_xxx_xyyyyz, g_y_0_xxx_xyyyzz, g_y_0_xxx_xyyzzz, g_y_0_xxx_xyzzzz, g_y_0_xxx_xzzzzz, g_y_0_xxx_yyyyyy, g_y_0_xxx_yyyyyz, g_y_0_xxx_yyyyzz, g_y_0_xxx_yyyzzz, g_y_0_xxx_yyzzzz, g_y_0_xxx_yzzzzz, g_y_0_xxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxx_xxxxxx[k] = -g_y_0_xx_xxxxxx[k] * ab_x + g_y_0_xx_xxxxxxx[k];

                g_y_0_xxx_xxxxxy[k] = -g_y_0_xx_xxxxxy[k] * ab_x + g_y_0_xx_xxxxxxy[k];

                g_y_0_xxx_xxxxxz[k] = -g_y_0_xx_xxxxxz[k] * ab_x + g_y_0_xx_xxxxxxz[k];

                g_y_0_xxx_xxxxyy[k] = -g_y_0_xx_xxxxyy[k] * ab_x + g_y_0_xx_xxxxxyy[k];

                g_y_0_xxx_xxxxyz[k] = -g_y_0_xx_xxxxyz[k] * ab_x + g_y_0_xx_xxxxxyz[k];

                g_y_0_xxx_xxxxzz[k] = -g_y_0_xx_xxxxzz[k] * ab_x + g_y_0_xx_xxxxxzz[k];

                g_y_0_xxx_xxxyyy[k] = -g_y_0_xx_xxxyyy[k] * ab_x + g_y_0_xx_xxxxyyy[k];

                g_y_0_xxx_xxxyyz[k] = -g_y_0_xx_xxxyyz[k] * ab_x + g_y_0_xx_xxxxyyz[k];

                g_y_0_xxx_xxxyzz[k] = -g_y_0_xx_xxxyzz[k] * ab_x + g_y_0_xx_xxxxyzz[k];

                g_y_0_xxx_xxxzzz[k] = -g_y_0_xx_xxxzzz[k] * ab_x + g_y_0_xx_xxxxzzz[k];

                g_y_0_xxx_xxyyyy[k] = -g_y_0_xx_xxyyyy[k] * ab_x + g_y_0_xx_xxxyyyy[k];

                g_y_0_xxx_xxyyyz[k] = -g_y_0_xx_xxyyyz[k] * ab_x + g_y_0_xx_xxxyyyz[k];

                g_y_0_xxx_xxyyzz[k] = -g_y_0_xx_xxyyzz[k] * ab_x + g_y_0_xx_xxxyyzz[k];

                g_y_0_xxx_xxyzzz[k] = -g_y_0_xx_xxyzzz[k] * ab_x + g_y_0_xx_xxxyzzz[k];

                g_y_0_xxx_xxzzzz[k] = -g_y_0_xx_xxzzzz[k] * ab_x + g_y_0_xx_xxxzzzz[k];

                g_y_0_xxx_xyyyyy[k] = -g_y_0_xx_xyyyyy[k] * ab_x + g_y_0_xx_xxyyyyy[k];

                g_y_0_xxx_xyyyyz[k] = -g_y_0_xx_xyyyyz[k] * ab_x + g_y_0_xx_xxyyyyz[k];

                g_y_0_xxx_xyyyzz[k] = -g_y_0_xx_xyyyzz[k] * ab_x + g_y_0_xx_xxyyyzz[k];

                g_y_0_xxx_xyyzzz[k] = -g_y_0_xx_xyyzzz[k] * ab_x + g_y_0_xx_xxyyzzz[k];

                g_y_0_xxx_xyzzzz[k] = -g_y_0_xx_xyzzzz[k] * ab_x + g_y_0_xx_xxyzzzz[k];

                g_y_0_xxx_xzzzzz[k] = -g_y_0_xx_xzzzzz[k] * ab_x + g_y_0_xx_xxzzzzz[k];

                g_y_0_xxx_yyyyyy[k] = -g_y_0_xx_yyyyyy[k] * ab_x + g_y_0_xx_xyyyyyy[k];

                g_y_0_xxx_yyyyyz[k] = -g_y_0_xx_yyyyyz[k] * ab_x + g_y_0_xx_xyyyyyz[k];

                g_y_0_xxx_yyyyzz[k] = -g_y_0_xx_yyyyzz[k] * ab_x + g_y_0_xx_xyyyyzz[k];

                g_y_0_xxx_yyyzzz[k] = -g_y_0_xx_yyyzzz[k] * ab_x + g_y_0_xx_xyyyzzz[k];

                g_y_0_xxx_yyzzzz[k] = -g_y_0_xx_yyzzzz[k] * ab_x + g_y_0_xx_xyyzzzz[k];

                g_y_0_xxx_yzzzzz[k] = -g_y_0_xx_yzzzzz[k] * ab_x + g_y_0_xx_xyzzzzz[k];

                g_y_0_xxx_zzzzzz[k] = -g_y_0_xx_zzzzzz[k] * ab_x + g_y_0_xx_xzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxy_xxxxxx = cbuffer.data(fi_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxy = cbuffer.data(fi_geom_10_off + 309 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxxz = cbuffer.data(fi_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxyy = cbuffer.data(fi_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxyz = cbuffer.data(fi_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_xxy_xxxxzz = cbuffer.data(fi_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyyy = cbuffer.data(fi_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyyz = cbuffer.data(fi_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xxy_xxxyzz = cbuffer.data(fi_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xxy_xxxzzz = cbuffer.data(fi_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyyy = cbuffer.data(fi_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyyz = cbuffer.data(fi_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_xxy_xxyyzz = cbuffer.data(fi_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xxy_xxyzzz = cbuffer.data(fi_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xxy_xxzzzz = cbuffer.data(fi_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyyy = cbuffer.data(fi_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyyz = cbuffer.data(fi_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xxy_xyyyzz = cbuffer.data(fi_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xxy_xyyzzz = cbuffer.data(fi_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xxy_xyzzzz = cbuffer.data(fi_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xxy_xzzzzz = cbuffer.data(fi_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyyy = cbuffer.data(fi_geom_10_off + 329 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyyz = cbuffer.data(fi_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xxy_yyyyzz = cbuffer.data(fi_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xxy_yyyzzz = cbuffer.data(fi_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xxy_yyzzzz = cbuffer.data(fi_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xxy_yzzzzz = cbuffer.data(fi_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xxy_zzzzzz = cbuffer.data(fi_geom_10_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxy_xxxxxx, g_y_0_xxy_xxxxxy, g_y_0_xxy_xxxxxz, g_y_0_xxy_xxxxyy, g_y_0_xxy_xxxxyz, g_y_0_xxy_xxxxzz, g_y_0_xxy_xxxyyy, g_y_0_xxy_xxxyyz, g_y_0_xxy_xxxyzz, g_y_0_xxy_xxxzzz, g_y_0_xxy_xxyyyy, g_y_0_xxy_xxyyyz, g_y_0_xxy_xxyyzz, g_y_0_xxy_xxyzzz, g_y_0_xxy_xxzzzz, g_y_0_xxy_xyyyyy, g_y_0_xxy_xyyyyz, g_y_0_xxy_xyyyzz, g_y_0_xxy_xyyzzz, g_y_0_xxy_xyzzzz, g_y_0_xxy_xzzzzz, g_y_0_xxy_yyyyyy, g_y_0_xxy_yyyyyz, g_y_0_xxy_yyyyzz, g_y_0_xxy_yyyzzz, g_y_0_xxy_yyzzzz, g_y_0_xxy_yzzzzz, g_y_0_xxy_zzzzzz, g_y_0_xy_xxxxxx, g_y_0_xy_xxxxxxx, g_y_0_xy_xxxxxxy, g_y_0_xy_xxxxxxz, g_y_0_xy_xxxxxy, g_y_0_xy_xxxxxyy, g_y_0_xy_xxxxxyz, g_y_0_xy_xxxxxz, g_y_0_xy_xxxxxzz, g_y_0_xy_xxxxyy, g_y_0_xy_xxxxyyy, g_y_0_xy_xxxxyyz, g_y_0_xy_xxxxyz, g_y_0_xy_xxxxyzz, g_y_0_xy_xxxxzz, g_y_0_xy_xxxxzzz, g_y_0_xy_xxxyyy, g_y_0_xy_xxxyyyy, g_y_0_xy_xxxyyyz, g_y_0_xy_xxxyyz, g_y_0_xy_xxxyyzz, g_y_0_xy_xxxyzz, g_y_0_xy_xxxyzzz, g_y_0_xy_xxxzzz, g_y_0_xy_xxxzzzz, g_y_0_xy_xxyyyy, g_y_0_xy_xxyyyyy, g_y_0_xy_xxyyyyz, g_y_0_xy_xxyyyz, g_y_0_xy_xxyyyzz, g_y_0_xy_xxyyzz, g_y_0_xy_xxyyzzz, g_y_0_xy_xxyzzz, g_y_0_xy_xxyzzzz, g_y_0_xy_xxzzzz, g_y_0_xy_xxzzzzz, g_y_0_xy_xyyyyy, g_y_0_xy_xyyyyyy, g_y_0_xy_xyyyyyz, g_y_0_xy_xyyyyz, g_y_0_xy_xyyyyzz, g_y_0_xy_xyyyzz, g_y_0_xy_xyyyzzz, g_y_0_xy_xyyzzz, g_y_0_xy_xyyzzzz, g_y_0_xy_xyzzzz, g_y_0_xy_xyzzzzz, g_y_0_xy_xzzzzz, g_y_0_xy_xzzzzzz, g_y_0_xy_yyyyyy, g_y_0_xy_yyyyyz, g_y_0_xy_yyyyzz, g_y_0_xy_yyyzzz, g_y_0_xy_yyzzzz, g_y_0_xy_yzzzzz, g_y_0_xy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxy_xxxxxx[k] = -g_y_0_xy_xxxxxx[k] * ab_x + g_y_0_xy_xxxxxxx[k];

                g_y_0_xxy_xxxxxy[k] = -g_y_0_xy_xxxxxy[k] * ab_x + g_y_0_xy_xxxxxxy[k];

                g_y_0_xxy_xxxxxz[k] = -g_y_0_xy_xxxxxz[k] * ab_x + g_y_0_xy_xxxxxxz[k];

                g_y_0_xxy_xxxxyy[k] = -g_y_0_xy_xxxxyy[k] * ab_x + g_y_0_xy_xxxxxyy[k];

                g_y_0_xxy_xxxxyz[k] = -g_y_0_xy_xxxxyz[k] * ab_x + g_y_0_xy_xxxxxyz[k];

                g_y_0_xxy_xxxxzz[k] = -g_y_0_xy_xxxxzz[k] * ab_x + g_y_0_xy_xxxxxzz[k];

                g_y_0_xxy_xxxyyy[k] = -g_y_0_xy_xxxyyy[k] * ab_x + g_y_0_xy_xxxxyyy[k];

                g_y_0_xxy_xxxyyz[k] = -g_y_0_xy_xxxyyz[k] * ab_x + g_y_0_xy_xxxxyyz[k];

                g_y_0_xxy_xxxyzz[k] = -g_y_0_xy_xxxyzz[k] * ab_x + g_y_0_xy_xxxxyzz[k];

                g_y_0_xxy_xxxzzz[k] = -g_y_0_xy_xxxzzz[k] * ab_x + g_y_0_xy_xxxxzzz[k];

                g_y_0_xxy_xxyyyy[k] = -g_y_0_xy_xxyyyy[k] * ab_x + g_y_0_xy_xxxyyyy[k];

                g_y_0_xxy_xxyyyz[k] = -g_y_0_xy_xxyyyz[k] * ab_x + g_y_0_xy_xxxyyyz[k];

                g_y_0_xxy_xxyyzz[k] = -g_y_0_xy_xxyyzz[k] * ab_x + g_y_0_xy_xxxyyzz[k];

                g_y_0_xxy_xxyzzz[k] = -g_y_0_xy_xxyzzz[k] * ab_x + g_y_0_xy_xxxyzzz[k];

                g_y_0_xxy_xxzzzz[k] = -g_y_0_xy_xxzzzz[k] * ab_x + g_y_0_xy_xxxzzzz[k];

                g_y_0_xxy_xyyyyy[k] = -g_y_0_xy_xyyyyy[k] * ab_x + g_y_0_xy_xxyyyyy[k];

                g_y_0_xxy_xyyyyz[k] = -g_y_0_xy_xyyyyz[k] * ab_x + g_y_0_xy_xxyyyyz[k];

                g_y_0_xxy_xyyyzz[k] = -g_y_0_xy_xyyyzz[k] * ab_x + g_y_0_xy_xxyyyzz[k];

                g_y_0_xxy_xyyzzz[k] = -g_y_0_xy_xyyzzz[k] * ab_x + g_y_0_xy_xxyyzzz[k];

                g_y_0_xxy_xyzzzz[k] = -g_y_0_xy_xyzzzz[k] * ab_x + g_y_0_xy_xxyzzzz[k];

                g_y_0_xxy_xzzzzz[k] = -g_y_0_xy_xzzzzz[k] * ab_x + g_y_0_xy_xxzzzzz[k];

                g_y_0_xxy_yyyyyy[k] = -g_y_0_xy_yyyyyy[k] * ab_x + g_y_0_xy_xyyyyyy[k];

                g_y_0_xxy_yyyyyz[k] = -g_y_0_xy_yyyyyz[k] * ab_x + g_y_0_xy_xyyyyyz[k];

                g_y_0_xxy_yyyyzz[k] = -g_y_0_xy_yyyyzz[k] * ab_x + g_y_0_xy_xyyyyzz[k];

                g_y_0_xxy_yyyzzz[k] = -g_y_0_xy_yyyzzz[k] * ab_x + g_y_0_xy_xyyyzzz[k];

                g_y_0_xxy_yyzzzz[k] = -g_y_0_xy_yyzzzz[k] * ab_x + g_y_0_xy_xyyzzzz[k];

                g_y_0_xxy_yzzzzz[k] = -g_y_0_xy_yzzzzz[k] * ab_x + g_y_0_xy_xyzzzzz[k];

                g_y_0_xxy_zzzzzz[k] = -g_y_0_xy_zzzzzz[k] * ab_x + g_y_0_xy_xzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxz_xxxxxx = cbuffer.data(fi_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxy = cbuffer.data(fi_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxxz = cbuffer.data(fi_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxyy = cbuffer.data(fi_geom_10_off + 339 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxyz = cbuffer.data(fi_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_xxz_xxxxzz = cbuffer.data(fi_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyyy = cbuffer.data(fi_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyyz = cbuffer.data(fi_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_xxz_xxxyzz = cbuffer.data(fi_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_xxz_xxxzzz = cbuffer.data(fi_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyyy = cbuffer.data(fi_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyyz = cbuffer.data(fi_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_xxz_xxyyzz = cbuffer.data(fi_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_xxz_xxyzzz = cbuffer.data(fi_geom_10_off + 349 * ccomps * dcomps);

            auto g_y_0_xxz_xxzzzz = cbuffer.data(fi_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyyy = cbuffer.data(fi_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyyz = cbuffer.data(fi_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_xxz_xyyyzz = cbuffer.data(fi_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_xxz_xyyzzz = cbuffer.data(fi_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_xxz_xyzzzz = cbuffer.data(fi_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_xxz_xzzzzz = cbuffer.data(fi_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyyy = cbuffer.data(fi_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyyz = cbuffer.data(fi_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_xxz_yyyyzz = cbuffer.data(fi_geom_10_off + 359 * ccomps * dcomps);

            auto g_y_0_xxz_yyyzzz = cbuffer.data(fi_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_xxz_yyzzzz = cbuffer.data(fi_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_xxz_yzzzzz = cbuffer.data(fi_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_xxz_zzzzzz = cbuffer.data(fi_geom_10_off + 363 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxz_xxxxxx, g_y_0_xxz_xxxxxy, g_y_0_xxz_xxxxxz, g_y_0_xxz_xxxxyy, g_y_0_xxz_xxxxyz, g_y_0_xxz_xxxxzz, g_y_0_xxz_xxxyyy, g_y_0_xxz_xxxyyz, g_y_0_xxz_xxxyzz, g_y_0_xxz_xxxzzz, g_y_0_xxz_xxyyyy, g_y_0_xxz_xxyyyz, g_y_0_xxz_xxyyzz, g_y_0_xxz_xxyzzz, g_y_0_xxz_xxzzzz, g_y_0_xxz_xyyyyy, g_y_0_xxz_xyyyyz, g_y_0_xxz_xyyyzz, g_y_0_xxz_xyyzzz, g_y_0_xxz_xyzzzz, g_y_0_xxz_xzzzzz, g_y_0_xxz_yyyyyy, g_y_0_xxz_yyyyyz, g_y_0_xxz_yyyyzz, g_y_0_xxz_yyyzzz, g_y_0_xxz_yyzzzz, g_y_0_xxz_yzzzzz, g_y_0_xxz_zzzzzz, g_y_0_xz_xxxxxx, g_y_0_xz_xxxxxxx, g_y_0_xz_xxxxxxy, g_y_0_xz_xxxxxxz, g_y_0_xz_xxxxxy, g_y_0_xz_xxxxxyy, g_y_0_xz_xxxxxyz, g_y_0_xz_xxxxxz, g_y_0_xz_xxxxxzz, g_y_0_xz_xxxxyy, g_y_0_xz_xxxxyyy, g_y_0_xz_xxxxyyz, g_y_0_xz_xxxxyz, g_y_0_xz_xxxxyzz, g_y_0_xz_xxxxzz, g_y_0_xz_xxxxzzz, g_y_0_xz_xxxyyy, g_y_0_xz_xxxyyyy, g_y_0_xz_xxxyyyz, g_y_0_xz_xxxyyz, g_y_0_xz_xxxyyzz, g_y_0_xz_xxxyzz, g_y_0_xz_xxxyzzz, g_y_0_xz_xxxzzz, g_y_0_xz_xxxzzzz, g_y_0_xz_xxyyyy, g_y_0_xz_xxyyyyy, g_y_0_xz_xxyyyyz, g_y_0_xz_xxyyyz, g_y_0_xz_xxyyyzz, g_y_0_xz_xxyyzz, g_y_0_xz_xxyyzzz, g_y_0_xz_xxyzzz, g_y_0_xz_xxyzzzz, g_y_0_xz_xxzzzz, g_y_0_xz_xxzzzzz, g_y_0_xz_xyyyyy, g_y_0_xz_xyyyyyy, g_y_0_xz_xyyyyyz, g_y_0_xz_xyyyyz, g_y_0_xz_xyyyyzz, g_y_0_xz_xyyyzz, g_y_0_xz_xyyyzzz, g_y_0_xz_xyyzzz, g_y_0_xz_xyyzzzz, g_y_0_xz_xyzzzz, g_y_0_xz_xyzzzzz, g_y_0_xz_xzzzzz, g_y_0_xz_xzzzzzz, g_y_0_xz_yyyyyy, g_y_0_xz_yyyyyz, g_y_0_xz_yyyyzz, g_y_0_xz_yyyzzz, g_y_0_xz_yyzzzz, g_y_0_xz_yzzzzz, g_y_0_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxz_xxxxxx[k] = -g_y_0_xz_xxxxxx[k] * ab_x + g_y_0_xz_xxxxxxx[k];

                g_y_0_xxz_xxxxxy[k] = -g_y_0_xz_xxxxxy[k] * ab_x + g_y_0_xz_xxxxxxy[k];

                g_y_0_xxz_xxxxxz[k] = -g_y_0_xz_xxxxxz[k] * ab_x + g_y_0_xz_xxxxxxz[k];

                g_y_0_xxz_xxxxyy[k] = -g_y_0_xz_xxxxyy[k] * ab_x + g_y_0_xz_xxxxxyy[k];

                g_y_0_xxz_xxxxyz[k] = -g_y_0_xz_xxxxyz[k] * ab_x + g_y_0_xz_xxxxxyz[k];

                g_y_0_xxz_xxxxzz[k] = -g_y_0_xz_xxxxzz[k] * ab_x + g_y_0_xz_xxxxxzz[k];

                g_y_0_xxz_xxxyyy[k] = -g_y_0_xz_xxxyyy[k] * ab_x + g_y_0_xz_xxxxyyy[k];

                g_y_0_xxz_xxxyyz[k] = -g_y_0_xz_xxxyyz[k] * ab_x + g_y_0_xz_xxxxyyz[k];

                g_y_0_xxz_xxxyzz[k] = -g_y_0_xz_xxxyzz[k] * ab_x + g_y_0_xz_xxxxyzz[k];

                g_y_0_xxz_xxxzzz[k] = -g_y_0_xz_xxxzzz[k] * ab_x + g_y_0_xz_xxxxzzz[k];

                g_y_0_xxz_xxyyyy[k] = -g_y_0_xz_xxyyyy[k] * ab_x + g_y_0_xz_xxxyyyy[k];

                g_y_0_xxz_xxyyyz[k] = -g_y_0_xz_xxyyyz[k] * ab_x + g_y_0_xz_xxxyyyz[k];

                g_y_0_xxz_xxyyzz[k] = -g_y_0_xz_xxyyzz[k] * ab_x + g_y_0_xz_xxxyyzz[k];

                g_y_0_xxz_xxyzzz[k] = -g_y_0_xz_xxyzzz[k] * ab_x + g_y_0_xz_xxxyzzz[k];

                g_y_0_xxz_xxzzzz[k] = -g_y_0_xz_xxzzzz[k] * ab_x + g_y_0_xz_xxxzzzz[k];

                g_y_0_xxz_xyyyyy[k] = -g_y_0_xz_xyyyyy[k] * ab_x + g_y_0_xz_xxyyyyy[k];

                g_y_0_xxz_xyyyyz[k] = -g_y_0_xz_xyyyyz[k] * ab_x + g_y_0_xz_xxyyyyz[k];

                g_y_0_xxz_xyyyzz[k] = -g_y_0_xz_xyyyzz[k] * ab_x + g_y_0_xz_xxyyyzz[k];

                g_y_0_xxz_xyyzzz[k] = -g_y_0_xz_xyyzzz[k] * ab_x + g_y_0_xz_xxyyzzz[k];

                g_y_0_xxz_xyzzzz[k] = -g_y_0_xz_xyzzzz[k] * ab_x + g_y_0_xz_xxyzzzz[k];

                g_y_0_xxz_xzzzzz[k] = -g_y_0_xz_xzzzzz[k] * ab_x + g_y_0_xz_xxzzzzz[k];

                g_y_0_xxz_yyyyyy[k] = -g_y_0_xz_yyyyyy[k] * ab_x + g_y_0_xz_xyyyyyy[k];

                g_y_0_xxz_yyyyyz[k] = -g_y_0_xz_yyyyyz[k] * ab_x + g_y_0_xz_xyyyyyz[k];

                g_y_0_xxz_yyyyzz[k] = -g_y_0_xz_yyyyzz[k] * ab_x + g_y_0_xz_xyyyyzz[k];

                g_y_0_xxz_yyyzzz[k] = -g_y_0_xz_yyyzzz[k] * ab_x + g_y_0_xz_xyyyzzz[k];

                g_y_0_xxz_yyzzzz[k] = -g_y_0_xz_yyzzzz[k] * ab_x + g_y_0_xz_xyyzzzz[k];

                g_y_0_xxz_yzzzzz[k] = -g_y_0_xz_yzzzzz[k] * ab_x + g_y_0_xz_xyzzzzz[k];

                g_y_0_xxz_zzzzzz[k] = -g_y_0_xz_zzzzzz[k] * ab_x + g_y_0_xz_xzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyy_xxxxxx = cbuffer.data(fi_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxy = cbuffer.data(fi_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxxz = cbuffer.data(fi_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxyy = cbuffer.data(fi_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxyz = cbuffer.data(fi_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_xyy_xxxxzz = cbuffer.data(fi_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyyy = cbuffer.data(fi_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyyz = cbuffer.data(fi_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_xyy_xxxyzz = cbuffer.data(fi_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_xyy_xxxzzz = cbuffer.data(fi_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyyy = cbuffer.data(fi_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyyz = cbuffer.data(fi_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_xyy_xxyyzz = cbuffer.data(fi_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_xyy_xxyzzz = cbuffer.data(fi_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_xyy_xxzzzz = cbuffer.data(fi_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyyy = cbuffer.data(fi_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyyz = cbuffer.data(fi_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_xyy_xyyyzz = cbuffer.data(fi_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_xyy_xyyzzz = cbuffer.data(fi_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_xyy_xyzzzz = cbuffer.data(fi_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_xyy_xzzzzz = cbuffer.data(fi_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyyy = cbuffer.data(fi_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyyz = cbuffer.data(fi_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_xyy_yyyyzz = cbuffer.data(fi_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_xyy_yyyzzz = cbuffer.data(fi_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_xyy_yyzzzz = cbuffer.data(fi_geom_10_off + 389 * ccomps * dcomps);

            auto g_y_0_xyy_yzzzzz = cbuffer.data(fi_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_xyy_zzzzzz = cbuffer.data(fi_geom_10_off + 391 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyy_xxxxxx, g_y_0_xyy_xxxxxy, g_y_0_xyy_xxxxxz, g_y_0_xyy_xxxxyy, g_y_0_xyy_xxxxyz, g_y_0_xyy_xxxxzz, g_y_0_xyy_xxxyyy, g_y_0_xyy_xxxyyz, g_y_0_xyy_xxxyzz, g_y_0_xyy_xxxzzz, g_y_0_xyy_xxyyyy, g_y_0_xyy_xxyyyz, g_y_0_xyy_xxyyzz, g_y_0_xyy_xxyzzz, g_y_0_xyy_xxzzzz, g_y_0_xyy_xyyyyy, g_y_0_xyy_xyyyyz, g_y_0_xyy_xyyyzz, g_y_0_xyy_xyyzzz, g_y_0_xyy_xyzzzz, g_y_0_xyy_xzzzzz, g_y_0_xyy_yyyyyy, g_y_0_xyy_yyyyyz, g_y_0_xyy_yyyyzz, g_y_0_xyy_yyyzzz, g_y_0_xyy_yyzzzz, g_y_0_xyy_yzzzzz, g_y_0_xyy_zzzzzz, g_y_0_yy_xxxxxx, g_y_0_yy_xxxxxxx, g_y_0_yy_xxxxxxy, g_y_0_yy_xxxxxxz, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxxyy, g_y_0_yy_xxxxxyz, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxxzz, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyyy, g_y_0_yy_xxxxyyz, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxyzz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxxzzz, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyyy, g_y_0_yy_xxxyyyz, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyyzz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxyzzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxxzzzz, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyyy, g_y_0_yy_xxyyyyz, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyyzz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyyzzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxyzzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xxzzzzz, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyyy, g_y_0_yy_xyyyyyz, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyyzz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyyzzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyyzzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xyzzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_xzzzzzz, g_y_0_yy_yyyyyy, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyy_xxxxxx[k] = -g_y_0_yy_xxxxxx[k] * ab_x + g_y_0_yy_xxxxxxx[k];

                g_y_0_xyy_xxxxxy[k] = -g_y_0_yy_xxxxxy[k] * ab_x + g_y_0_yy_xxxxxxy[k];

                g_y_0_xyy_xxxxxz[k] = -g_y_0_yy_xxxxxz[k] * ab_x + g_y_0_yy_xxxxxxz[k];

                g_y_0_xyy_xxxxyy[k] = -g_y_0_yy_xxxxyy[k] * ab_x + g_y_0_yy_xxxxxyy[k];

                g_y_0_xyy_xxxxyz[k] = -g_y_0_yy_xxxxyz[k] * ab_x + g_y_0_yy_xxxxxyz[k];

                g_y_0_xyy_xxxxzz[k] = -g_y_0_yy_xxxxzz[k] * ab_x + g_y_0_yy_xxxxxzz[k];

                g_y_0_xyy_xxxyyy[k] = -g_y_0_yy_xxxyyy[k] * ab_x + g_y_0_yy_xxxxyyy[k];

                g_y_0_xyy_xxxyyz[k] = -g_y_0_yy_xxxyyz[k] * ab_x + g_y_0_yy_xxxxyyz[k];

                g_y_0_xyy_xxxyzz[k] = -g_y_0_yy_xxxyzz[k] * ab_x + g_y_0_yy_xxxxyzz[k];

                g_y_0_xyy_xxxzzz[k] = -g_y_0_yy_xxxzzz[k] * ab_x + g_y_0_yy_xxxxzzz[k];

                g_y_0_xyy_xxyyyy[k] = -g_y_0_yy_xxyyyy[k] * ab_x + g_y_0_yy_xxxyyyy[k];

                g_y_0_xyy_xxyyyz[k] = -g_y_0_yy_xxyyyz[k] * ab_x + g_y_0_yy_xxxyyyz[k];

                g_y_0_xyy_xxyyzz[k] = -g_y_0_yy_xxyyzz[k] * ab_x + g_y_0_yy_xxxyyzz[k];

                g_y_0_xyy_xxyzzz[k] = -g_y_0_yy_xxyzzz[k] * ab_x + g_y_0_yy_xxxyzzz[k];

                g_y_0_xyy_xxzzzz[k] = -g_y_0_yy_xxzzzz[k] * ab_x + g_y_0_yy_xxxzzzz[k];

                g_y_0_xyy_xyyyyy[k] = -g_y_0_yy_xyyyyy[k] * ab_x + g_y_0_yy_xxyyyyy[k];

                g_y_0_xyy_xyyyyz[k] = -g_y_0_yy_xyyyyz[k] * ab_x + g_y_0_yy_xxyyyyz[k];

                g_y_0_xyy_xyyyzz[k] = -g_y_0_yy_xyyyzz[k] * ab_x + g_y_0_yy_xxyyyzz[k];

                g_y_0_xyy_xyyzzz[k] = -g_y_0_yy_xyyzzz[k] * ab_x + g_y_0_yy_xxyyzzz[k];

                g_y_0_xyy_xyzzzz[k] = -g_y_0_yy_xyzzzz[k] * ab_x + g_y_0_yy_xxyzzzz[k];

                g_y_0_xyy_xzzzzz[k] = -g_y_0_yy_xzzzzz[k] * ab_x + g_y_0_yy_xxzzzzz[k];

                g_y_0_xyy_yyyyyy[k] = -g_y_0_yy_yyyyyy[k] * ab_x + g_y_0_yy_xyyyyyy[k];

                g_y_0_xyy_yyyyyz[k] = -g_y_0_yy_yyyyyz[k] * ab_x + g_y_0_yy_xyyyyyz[k];

                g_y_0_xyy_yyyyzz[k] = -g_y_0_yy_yyyyzz[k] * ab_x + g_y_0_yy_xyyyyzz[k];

                g_y_0_xyy_yyyzzz[k] = -g_y_0_yy_yyyzzz[k] * ab_x + g_y_0_yy_xyyyzzz[k];

                g_y_0_xyy_yyzzzz[k] = -g_y_0_yy_yyzzzz[k] * ab_x + g_y_0_yy_xyyzzzz[k];

                g_y_0_xyy_yzzzzz[k] = -g_y_0_yy_yzzzzz[k] * ab_x + g_y_0_yy_xyzzzzz[k];

                g_y_0_xyy_zzzzzz[k] = -g_y_0_yy_zzzzzz[k] * ab_x + g_y_0_yy_xzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyz_xxxxxx = cbuffer.data(fi_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxy = cbuffer.data(fi_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxxz = cbuffer.data(fi_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxyy = cbuffer.data(fi_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxyz = cbuffer.data(fi_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_xyz_xxxxzz = cbuffer.data(fi_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyyy = cbuffer.data(fi_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyyz = cbuffer.data(fi_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_xyz_xxxyzz = cbuffer.data(fi_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_xyz_xxxzzz = cbuffer.data(fi_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyyy = cbuffer.data(fi_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyyz = cbuffer.data(fi_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_xyz_xxyyzz = cbuffer.data(fi_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_xyz_xxyzzz = cbuffer.data(fi_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_xyz_xxzzzz = cbuffer.data(fi_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyyy = cbuffer.data(fi_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyyz = cbuffer.data(fi_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_xyz_xyyyzz = cbuffer.data(fi_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_xyz_xyyzzz = cbuffer.data(fi_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_xyz_xyzzzz = cbuffer.data(fi_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_xyz_xzzzzz = cbuffer.data(fi_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyyy = cbuffer.data(fi_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyyz = cbuffer.data(fi_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_xyz_yyyyzz = cbuffer.data(fi_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_xyz_yyyzzz = cbuffer.data(fi_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_xyz_yyzzzz = cbuffer.data(fi_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_xyz_yzzzzz = cbuffer.data(fi_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_xyz_zzzzzz = cbuffer.data(fi_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyz_xxxxxx, g_y_0_xyz_xxxxxy, g_y_0_xyz_xxxxxz, g_y_0_xyz_xxxxyy, g_y_0_xyz_xxxxyz, g_y_0_xyz_xxxxzz, g_y_0_xyz_xxxyyy, g_y_0_xyz_xxxyyz, g_y_0_xyz_xxxyzz, g_y_0_xyz_xxxzzz, g_y_0_xyz_xxyyyy, g_y_0_xyz_xxyyyz, g_y_0_xyz_xxyyzz, g_y_0_xyz_xxyzzz, g_y_0_xyz_xxzzzz, g_y_0_xyz_xyyyyy, g_y_0_xyz_xyyyyz, g_y_0_xyz_xyyyzz, g_y_0_xyz_xyyzzz, g_y_0_xyz_xyzzzz, g_y_0_xyz_xzzzzz, g_y_0_xyz_yyyyyy, g_y_0_xyz_yyyyyz, g_y_0_xyz_yyyyzz, g_y_0_xyz_yyyzzz, g_y_0_xyz_yyzzzz, g_y_0_xyz_yzzzzz, g_y_0_xyz_zzzzzz, g_y_0_yz_xxxxxx, g_y_0_yz_xxxxxxx, g_y_0_yz_xxxxxxy, g_y_0_yz_xxxxxxz, g_y_0_yz_xxxxxy, g_y_0_yz_xxxxxyy, g_y_0_yz_xxxxxyz, g_y_0_yz_xxxxxz, g_y_0_yz_xxxxxzz, g_y_0_yz_xxxxyy, g_y_0_yz_xxxxyyy, g_y_0_yz_xxxxyyz, g_y_0_yz_xxxxyz, g_y_0_yz_xxxxyzz, g_y_0_yz_xxxxzz, g_y_0_yz_xxxxzzz, g_y_0_yz_xxxyyy, g_y_0_yz_xxxyyyy, g_y_0_yz_xxxyyyz, g_y_0_yz_xxxyyz, g_y_0_yz_xxxyyzz, g_y_0_yz_xxxyzz, g_y_0_yz_xxxyzzz, g_y_0_yz_xxxzzz, g_y_0_yz_xxxzzzz, g_y_0_yz_xxyyyy, g_y_0_yz_xxyyyyy, g_y_0_yz_xxyyyyz, g_y_0_yz_xxyyyz, g_y_0_yz_xxyyyzz, g_y_0_yz_xxyyzz, g_y_0_yz_xxyyzzz, g_y_0_yz_xxyzzz, g_y_0_yz_xxyzzzz, g_y_0_yz_xxzzzz, g_y_0_yz_xxzzzzz, g_y_0_yz_xyyyyy, g_y_0_yz_xyyyyyy, g_y_0_yz_xyyyyyz, g_y_0_yz_xyyyyz, g_y_0_yz_xyyyyzz, g_y_0_yz_xyyyzz, g_y_0_yz_xyyyzzz, g_y_0_yz_xyyzzz, g_y_0_yz_xyyzzzz, g_y_0_yz_xyzzzz, g_y_0_yz_xyzzzzz, g_y_0_yz_xzzzzz, g_y_0_yz_xzzzzzz, g_y_0_yz_yyyyyy, g_y_0_yz_yyyyyz, g_y_0_yz_yyyyzz, g_y_0_yz_yyyzzz, g_y_0_yz_yyzzzz, g_y_0_yz_yzzzzz, g_y_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyz_xxxxxx[k] = -g_y_0_yz_xxxxxx[k] * ab_x + g_y_0_yz_xxxxxxx[k];

                g_y_0_xyz_xxxxxy[k] = -g_y_0_yz_xxxxxy[k] * ab_x + g_y_0_yz_xxxxxxy[k];

                g_y_0_xyz_xxxxxz[k] = -g_y_0_yz_xxxxxz[k] * ab_x + g_y_0_yz_xxxxxxz[k];

                g_y_0_xyz_xxxxyy[k] = -g_y_0_yz_xxxxyy[k] * ab_x + g_y_0_yz_xxxxxyy[k];

                g_y_0_xyz_xxxxyz[k] = -g_y_0_yz_xxxxyz[k] * ab_x + g_y_0_yz_xxxxxyz[k];

                g_y_0_xyz_xxxxzz[k] = -g_y_0_yz_xxxxzz[k] * ab_x + g_y_0_yz_xxxxxzz[k];

                g_y_0_xyz_xxxyyy[k] = -g_y_0_yz_xxxyyy[k] * ab_x + g_y_0_yz_xxxxyyy[k];

                g_y_0_xyz_xxxyyz[k] = -g_y_0_yz_xxxyyz[k] * ab_x + g_y_0_yz_xxxxyyz[k];

                g_y_0_xyz_xxxyzz[k] = -g_y_0_yz_xxxyzz[k] * ab_x + g_y_0_yz_xxxxyzz[k];

                g_y_0_xyz_xxxzzz[k] = -g_y_0_yz_xxxzzz[k] * ab_x + g_y_0_yz_xxxxzzz[k];

                g_y_0_xyz_xxyyyy[k] = -g_y_0_yz_xxyyyy[k] * ab_x + g_y_0_yz_xxxyyyy[k];

                g_y_0_xyz_xxyyyz[k] = -g_y_0_yz_xxyyyz[k] * ab_x + g_y_0_yz_xxxyyyz[k];

                g_y_0_xyz_xxyyzz[k] = -g_y_0_yz_xxyyzz[k] * ab_x + g_y_0_yz_xxxyyzz[k];

                g_y_0_xyz_xxyzzz[k] = -g_y_0_yz_xxyzzz[k] * ab_x + g_y_0_yz_xxxyzzz[k];

                g_y_0_xyz_xxzzzz[k] = -g_y_0_yz_xxzzzz[k] * ab_x + g_y_0_yz_xxxzzzz[k];

                g_y_0_xyz_xyyyyy[k] = -g_y_0_yz_xyyyyy[k] * ab_x + g_y_0_yz_xxyyyyy[k];

                g_y_0_xyz_xyyyyz[k] = -g_y_0_yz_xyyyyz[k] * ab_x + g_y_0_yz_xxyyyyz[k];

                g_y_0_xyz_xyyyzz[k] = -g_y_0_yz_xyyyzz[k] * ab_x + g_y_0_yz_xxyyyzz[k];

                g_y_0_xyz_xyyzzz[k] = -g_y_0_yz_xyyzzz[k] * ab_x + g_y_0_yz_xxyyzzz[k];

                g_y_0_xyz_xyzzzz[k] = -g_y_0_yz_xyzzzz[k] * ab_x + g_y_0_yz_xxyzzzz[k];

                g_y_0_xyz_xzzzzz[k] = -g_y_0_yz_xzzzzz[k] * ab_x + g_y_0_yz_xxzzzzz[k];

                g_y_0_xyz_yyyyyy[k] = -g_y_0_yz_yyyyyy[k] * ab_x + g_y_0_yz_xyyyyyy[k];

                g_y_0_xyz_yyyyyz[k] = -g_y_0_yz_yyyyyz[k] * ab_x + g_y_0_yz_xyyyyyz[k];

                g_y_0_xyz_yyyyzz[k] = -g_y_0_yz_yyyyzz[k] * ab_x + g_y_0_yz_xyyyyzz[k];

                g_y_0_xyz_yyyzzz[k] = -g_y_0_yz_yyyzzz[k] * ab_x + g_y_0_yz_xyyyzzz[k];

                g_y_0_xyz_yyzzzz[k] = -g_y_0_yz_yyzzzz[k] * ab_x + g_y_0_yz_xyyzzzz[k];

                g_y_0_xyz_yzzzzz[k] = -g_y_0_yz_yzzzzz[k] * ab_x + g_y_0_yz_xyzzzzz[k];

                g_y_0_xyz_zzzzzz[k] = -g_y_0_yz_zzzzzz[k] * ab_x + g_y_0_yz_xzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzz_xxxxxx = cbuffer.data(fi_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxy = cbuffer.data(fi_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxxz = cbuffer.data(fi_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxyy = cbuffer.data(fi_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxyz = cbuffer.data(fi_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_xzz_xxxxzz = cbuffer.data(fi_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyyy = cbuffer.data(fi_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyyz = cbuffer.data(fi_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_xzz_xxxyzz = cbuffer.data(fi_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_xzz_xxxzzz = cbuffer.data(fi_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyyy = cbuffer.data(fi_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyyz = cbuffer.data(fi_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_xzz_xxyyzz = cbuffer.data(fi_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_xzz_xxyzzz = cbuffer.data(fi_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_xzz_xxzzzz = cbuffer.data(fi_geom_10_off + 434 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyyy = cbuffer.data(fi_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyyz = cbuffer.data(fi_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_xzz_xyyyzz = cbuffer.data(fi_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_xzz_xyyzzz = cbuffer.data(fi_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_xzz_xyzzzz = cbuffer.data(fi_geom_10_off + 439 * ccomps * dcomps);

            auto g_y_0_xzz_xzzzzz = cbuffer.data(fi_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyyy = cbuffer.data(fi_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyyz = cbuffer.data(fi_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xzz_yyyyzz = cbuffer.data(fi_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xzz_yyyzzz = cbuffer.data(fi_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xzz_yyzzzz = cbuffer.data(fi_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xzz_yzzzzz = cbuffer.data(fi_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xzz_zzzzzz = cbuffer.data(fi_geom_10_off + 447 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzz_xxxxxx, g_y_0_xzz_xxxxxy, g_y_0_xzz_xxxxxz, g_y_0_xzz_xxxxyy, g_y_0_xzz_xxxxyz, g_y_0_xzz_xxxxzz, g_y_0_xzz_xxxyyy, g_y_0_xzz_xxxyyz, g_y_0_xzz_xxxyzz, g_y_0_xzz_xxxzzz, g_y_0_xzz_xxyyyy, g_y_0_xzz_xxyyyz, g_y_0_xzz_xxyyzz, g_y_0_xzz_xxyzzz, g_y_0_xzz_xxzzzz, g_y_0_xzz_xyyyyy, g_y_0_xzz_xyyyyz, g_y_0_xzz_xyyyzz, g_y_0_xzz_xyyzzz, g_y_0_xzz_xyzzzz, g_y_0_xzz_xzzzzz, g_y_0_xzz_yyyyyy, g_y_0_xzz_yyyyyz, g_y_0_xzz_yyyyzz, g_y_0_xzz_yyyzzz, g_y_0_xzz_yyzzzz, g_y_0_xzz_yzzzzz, g_y_0_xzz_zzzzzz, g_y_0_zz_xxxxxx, g_y_0_zz_xxxxxxx, g_y_0_zz_xxxxxxy, g_y_0_zz_xxxxxxz, g_y_0_zz_xxxxxy, g_y_0_zz_xxxxxyy, g_y_0_zz_xxxxxyz, g_y_0_zz_xxxxxz, g_y_0_zz_xxxxxzz, g_y_0_zz_xxxxyy, g_y_0_zz_xxxxyyy, g_y_0_zz_xxxxyyz, g_y_0_zz_xxxxyz, g_y_0_zz_xxxxyzz, g_y_0_zz_xxxxzz, g_y_0_zz_xxxxzzz, g_y_0_zz_xxxyyy, g_y_0_zz_xxxyyyy, g_y_0_zz_xxxyyyz, g_y_0_zz_xxxyyz, g_y_0_zz_xxxyyzz, g_y_0_zz_xxxyzz, g_y_0_zz_xxxyzzz, g_y_0_zz_xxxzzz, g_y_0_zz_xxxzzzz, g_y_0_zz_xxyyyy, g_y_0_zz_xxyyyyy, g_y_0_zz_xxyyyyz, g_y_0_zz_xxyyyz, g_y_0_zz_xxyyyzz, g_y_0_zz_xxyyzz, g_y_0_zz_xxyyzzz, g_y_0_zz_xxyzzz, g_y_0_zz_xxyzzzz, g_y_0_zz_xxzzzz, g_y_0_zz_xxzzzzz, g_y_0_zz_xyyyyy, g_y_0_zz_xyyyyyy, g_y_0_zz_xyyyyyz, g_y_0_zz_xyyyyz, g_y_0_zz_xyyyyzz, g_y_0_zz_xyyyzz, g_y_0_zz_xyyyzzz, g_y_0_zz_xyyzzz, g_y_0_zz_xyyzzzz, g_y_0_zz_xyzzzz, g_y_0_zz_xyzzzzz, g_y_0_zz_xzzzzz, g_y_0_zz_xzzzzzz, g_y_0_zz_yyyyyy, g_y_0_zz_yyyyyz, g_y_0_zz_yyyyzz, g_y_0_zz_yyyzzz, g_y_0_zz_yyzzzz, g_y_0_zz_yzzzzz, g_y_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzz_xxxxxx[k] = -g_y_0_zz_xxxxxx[k] * ab_x + g_y_0_zz_xxxxxxx[k];

                g_y_0_xzz_xxxxxy[k] = -g_y_0_zz_xxxxxy[k] * ab_x + g_y_0_zz_xxxxxxy[k];

                g_y_0_xzz_xxxxxz[k] = -g_y_0_zz_xxxxxz[k] * ab_x + g_y_0_zz_xxxxxxz[k];

                g_y_0_xzz_xxxxyy[k] = -g_y_0_zz_xxxxyy[k] * ab_x + g_y_0_zz_xxxxxyy[k];

                g_y_0_xzz_xxxxyz[k] = -g_y_0_zz_xxxxyz[k] * ab_x + g_y_0_zz_xxxxxyz[k];

                g_y_0_xzz_xxxxzz[k] = -g_y_0_zz_xxxxzz[k] * ab_x + g_y_0_zz_xxxxxzz[k];

                g_y_0_xzz_xxxyyy[k] = -g_y_0_zz_xxxyyy[k] * ab_x + g_y_0_zz_xxxxyyy[k];

                g_y_0_xzz_xxxyyz[k] = -g_y_0_zz_xxxyyz[k] * ab_x + g_y_0_zz_xxxxyyz[k];

                g_y_0_xzz_xxxyzz[k] = -g_y_0_zz_xxxyzz[k] * ab_x + g_y_0_zz_xxxxyzz[k];

                g_y_0_xzz_xxxzzz[k] = -g_y_0_zz_xxxzzz[k] * ab_x + g_y_0_zz_xxxxzzz[k];

                g_y_0_xzz_xxyyyy[k] = -g_y_0_zz_xxyyyy[k] * ab_x + g_y_0_zz_xxxyyyy[k];

                g_y_0_xzz_xxyyyz[k] = -g_y_0_zz_xxyyyz[k] * ab_x + g_y_0_zz_xxxyyyz[k];

                g_y_0_xzz_xxyyzz[k] = -g_y_0_zz_xxyyzz[k] * ab_x + g_y_0_zz_xxxyyzz[k];

                g_y_0_xzz_xxyzzz[k] = -g_y_0_zz_xxyzzz[k] * ab_x + g_y_0_zz_xxxyzzz[k];

                g_y_0_xzz_xxzzzz[k] = -g_y_0_zz_xxzzzz[k] * ab_x + g_y_0_zz_xxxzzzz[k];

                g_y_0_xzz_xyyyyy[k] = -g_y_0_zz_xyyyyy[k] * ab_x + g_y_0_zz_xxyyyyy[k];

                g_y_0_xzz_xyyyyz[k] = -g_y_0_zz_xyyyyz[k] * ab_x + g_y_0_zz_xxyyyyz[k];

                g_y_0_xzz_xyyyzz[k] = -g_y_0_zz_xyyyzz[k] * ab_x + g_y_0_zz_xxyyyzz[k];

                g_y_0_xzz_xyyzzz[k] = -g_y_0_zz_xyyzzz[k] * ab_x + g_y_0_zz_xxyyzzz[k];

                g_y_0_xzz_xyzzzz[k] = -g_y_0_zz_xyzzzz[k] * ab_x + g_y_0_zz_xxyzzzz[k];

                g_y_0_xzz_xzzzzz[k] = -g_y_0_zz_xzzzzz[k] * ab_x + g_y_0_zz_xxzzzzz[k];

                g_y_0_xzz_yyyyyy[k] = -g_y_0_zz_yyyyyy[k] * ab_x + g_y_0_zz_xyyyyyy[k];

                g_y_0_xzz_yyyyyz[k] = -g_y_0_zz_yyyyyz[k] * ab_x + g_y_0_zz_xyyyyyz[k];

                g_y_0_xzz_yyyyzz[k] = -g_y_0_zz_yyyyzz[k] * ab_x + g_y_0_zz_xyyyyzz[k];

                g_y_0_xzz_yyyzzz[k] = -g_y_0_zz_yyyzzz[k] * ab_x + g_y_0_zz_xyyyzzz[k];

                g_y_0_xzz_yyzzzz[k] = -g_y_0_zz_yyzzzz[k] * ab_x + g_y_0_zz_xyyzzzz[k];

                g_y_0_xzz_yzzzzz[k] = -g_y_0_zz_yzzzzz[k] * ab_x + g_y_0_zz_xyzzzzz[k];

                g_y_0_xzz_zzzzzz[k] = -g_y_0_zz_zzzzzz[k] * ab_x + g_y_0_zz_xzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyy_xxxxxx = cbuffer.data(fi_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxy = cbuffer.data(fi_geom_10_off + 449 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxxz = cbuffer.data(fi_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxyy = cbuffer.data(fi_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxyz = cbuffer.data(fi_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_yyy_xxxxzz = cbuffer.data(fi_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyyy = cbuffer.data(fi_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyyz = cbuffer.data(fi_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_yyy_xxxyzz = cbuffer.data(fi_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_yyy_xxxzzz = cbuffer.data(fi_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyyy = cbuffer.data(fi_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyyz = cbuffer.data(fi_geom_10_off + 459 * ccomps * dcomps);

            auto g_y_0_yyy_xxyyzz = cbuffer.data(fi_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_yyy_xxyzzz = cbuffer.data(fi_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_yyy_xxzzzz = cbuffer.data(fi_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyyy = cbuffer.data(fi_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyyz = cbuffer.data(fi_geom_10_off + 464 * ccomps * dcomps);

            auto g_y_0_yyy_xyyyzz = cbuffer.data(fi_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_yyy_xyyzzz = cbuffer.data(fi_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_yyy_xyzzzz = cbuffer.data(fi_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_yyy_xzzzzz = cbuffer.data(fi_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyyy = cbuffer.data(fi_geom_10_off + 469 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyyz = cbuffer.data(fi_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_yyy_yyyyzz = cbuffer.data(fi_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_yyy_yyyzzz = cbuffer.data(fi_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_yyy_yyzzzz = cbuffer.data(fi_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_yyy_yzzzzz = cbuffer.data(fi_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_yyy_zzzzzz = cbuffer.data(fi_geom_10_off + 475 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_xxxxxx, g_y_0_yy_xxxxxxy, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxxyy, g_y_0_yy_xxxxxyz, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyyy, g_y_0_yy_xxxxyyz, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxyzz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyyy, g_y_0_yy_xxxyyyz, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyyzz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxyzzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyyy, g_y_0_yy_xxyyyyz, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyyzz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyyzzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxyzzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyyy, g_y_0_yy_xyyyyyz, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyyzz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyyzzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyyzzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xyzzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_yyyyyy, g_y_0_yy_yyyyyyy, g_y_0_yy_yyyyyyz, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyyzz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyyzzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyyzzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yyzzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_yzzzzzz, g_y_0_yy_zzzzzz, g_y_0_yyy_xxxxxx, g_y_0_yyy_xxxxxy, g_y_0_yyy_xxxxxz, g_y_0_yyy_xxxxyy, g_y_0_yyy_xxxxyz, g_y_0_yyy_xxxxzz, g_y_0_yyy_xxxyyy, g_y_0_yyy_xxxyyz, g_y_0_yyy_xxxyzz, g_y_0_yyy_xxxzzz, g_y_0_yyy_xxyyyy, g_y_0_yyy_xxyyyz, g_y_0_yyy_xxyyzz, g_y_0_yyy_xxyzzz, g_y_0_yyy_xxzzzz, g_y_0_yyy_xyyyyy, g_y_0_yyy_xyyyyz, g_y_0_yyy_xyyyzz, g_y_0_yyy_xyyzzz, g_y_0_yyy_xyzzzz, g_y_0_yyy_xzzzzz, g_y_0_yyy_yyyyyy, g_y_0_yyy_yyyyyz, g_y_0_yyy_yyyyzz, g_y_0_yyy_yyyzzz, g_y_0_yyy_yyzzzz, g_y_0_yyy_yzzzzz, g_y_0_yyy_zzzzzz, g_yy_xxxxxx, g_yy_xxxxxy, g_yy_xxxxxz, g_yy_xxxxyy, g_yy_xxxxyz, g_yy_xxxxzz, g_yy_xxxyyy, g_yy_xxxyyz, g_yy_xxxyzz, g_yy_xxxzzz, g_yy_xxyyyy, g_yy_xxyyyz, g_yy_xxyyzz, g_yy_xxyzzz, g_yy_xxzzzz, g_yy_xyyyyy, g_yy_xyyyyz, g_yy_xyyyzz, g_yy_xyyzzz, g_yy_xyzzzz, g_yy_xzzzzz, g_yy_yyyyyy, g_yy_yyyyyz, g_yy_yyyyzz, g_yy_yyyzzz, g_yy_yyzzzz, g_yy_yzzzzz, g_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyy_xxxxxx[k] = -g_yy_xxxxxx[k] - g_y_0_yy_xxxxxx[k] * ab_y + g_y_0_yy_xxxxxxy[k];

                g_y_0_yyy_xxxxxy[k] = -g_yy_xxxxxy[k] - g_y_0_yy_xxxxxy[k] * ab_y + g_y_0_yy_xxxxxyy[k];

                g_y_0_yyy_xxxxxz[k] = -g_yy_xxxxxz[k] - g_y_0_yy_xxxxxz[k] * ab_y + g_y_0_yy_xxxxxyz[k];

                g_y_0_yyy_xxxxyy[k] = -g_yy_xxxxyy[k] - g_y_0_yy_xxxxyy[k] * ab_y + g_y_0_yy_xxxxyyy[k];

                g_y_0_yyy_xxxxyz[k] = -g_yy_xxxxyz[k] - g_y_0_yy_xxxxyz[k] * ab_y + g_y_0_yy_xxxxyyz[k];

                g_y_0_yyy_xxxxzz[k] = -g_yy_xxxxzz[k] - g_y_0_yy_xxxxzz[k] * ab_y + g_y_0_yy_xxxxyzz[k];

                g_y_0_yyy_xxxyyy[k] = -g_yy_xxxyyy[k] - g_y_0_yy_xxxyyy[k] * ab_y + g_y_0_yy_xxxyyyy[k];

                g_y_0_yyy_xxxyyz[k] = -g_yy_xxxyyz[k] - g_y_0_yy_xxxyyz[k] * ab_y + g_y_0_yy_xxxyyyz[k];

                g_y_0_yyy_xxxyzz[k] = -g_yy_xxxyzz[k] - g_y_0_yy_xxxyzz[k] * ab_y + g_y_0_yy_xxxyyzz[k];

                g_y_0_yyy_xxxzzz[k] = -g_yy_xxxzzz[k] - g_y_0_yy_xxxzzz[k] * ab_y + g_y_0_yy_xxxyzzz[k];

                g_y_0_yyy_xxyyyy[k] = -g_yy_xxyyyy[k] - g_y_0_yy_xxyyyy[k] * ab_y + g_y_0_yy_xxyyyyy[k];

                g_y_0_yyy_xxyyyz[k] = -g_yy_xxyyyz[k] - g_y_0_yy_xxyyyz[k] * ab_y + g_y_0_yy_xxyyyyz[k];

                g_y_0_yyy_xxyyzz[k] = -g_yy_xxyyzz[k] - g_y_0_yy_xxyyzz[k] * ab_y + g_y_0_yy_xxyyyzz[k];

                g_y_0_yyy_xxyzzz[k] = -g_yy_xxyzzz[k] - g_y_0_yy_xxyzzz[k] * ab_y + g_y_0_yy_xxyyzzz[k];

                g_y_0_yyy_xxzzzz[k] = -g_yy_xxzzzz[k] - g_y_0_yy_xxzzzz[k] * ab_y + g_y_0_yy_xxyzzzz[k];

                g_y_0_yyy_xyyyyy[k] = -g_yy_xyyyyy[k] - g_y_0_yy_xyyyyy[k] * ab_y + g_y_0_yy_xyyyyyy[k];

                g_y_0_yyy_xyyyyz[k] = -g_yy_xyyyyz[k] - g_y_0_yy_xyyyyz[k] * ab_y + g_y_0_yy_xyyyyyz[k];

                g_y_0_yyy_xyyyzz[k] = -g_yy_xyyyzz[k] - g_y_0_yy_xyyyzz[k] * ab_y + g_y_0_yy_xyyyyzz[k];

                g_y_0_yyy_xyyzzz[k] = -g_yy_xyyzzz[k] - g_y_0_yy_xyyzzz[k] * ab_y + g_y_0_yy_xyyyzzz[k];

                g_y_0_yyy_xyzzzz[k] = -g_yy_xyzzzz[k] - g_y_0_yy_xyzzzz[k] * ab_y + g_y_0_yy_xyyzzzz[k];

                g_y_0_yyy_xzzzzz[k] = -g_yy_xzzzzz[k] - g_y_0_yy_xzzzzz[k] * ab_y + g_y_0_yy_xyzzzzz[k];

                g_y_0_yyy_yyyyyy[k] = -g_yy_yyyyyy[k] - g_y_0_yy_yyyyyy[k] * ab_y + g_y_0_yy_yyyyyyy[k];

                g_y_0_yyy_yyyyyz[k] = -g_yy_yyyyyz[k] - g_y_0_yy_yyyyyz[k] * ab_y + g_y_0_yy_yyyyyyz[k];

                g_y_0_yyy_yyyyzz[k] = -g_yy_yyyyzz[k] - g_y_0_yy_yyyyzz[k] * ab_y + g_y_0_yy_yyyyyzz[k];

                g_y_0_yyy_yyyzzz[k] = -g_yy_yyyzzz[k] - g_y_0_yy_yyyzzz[k] * ab_y + g_y_0_yy_yyyyzzz[k];

                g_y_0_yyy_yyzzzz[k] = -g_yy_yyzzzz[k] - g_y_0_yy_yyzzzz[k] * ab_y + g_y_0_yy_yyyzzzz[k];

                g_y_0_yyy_yzzzzz[k] = -g_yy_yzzzzz[k] - g_y_0_yy_yzzzzz[k] * ab_y + g_y_0_yy_yyzzzzz[k];

                g_y_0_yyy_zzzzzz[k] = -g_yy_zzzzzz[k] - g_y_0_yy_zzzzzz[k] * ab_y + g_y_0_yy_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyz_xxxxxx = cbuffer.data(fi_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxy = cbuffer.data(fi_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxxz = cbuffer.data(fi_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxyy = cbuffer.data(fi_geom_10_off + 479 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxyz = cbuffer.data(fi_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_yyz_xxxxzz = cbuffer.data(fi_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyyy = cbuffer.data(fi_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyyz = cbuffer.data(fi_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_yyz_xxxyzz = cbuffer.data(fi_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_yyz_xxxzzz = cbuffer.data(fi_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyyy = cbuffer.data(fi_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyyz = cbuffer.data(fi_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_yyz_xxyyzz = cbuffer.data(fi_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_yyz_xxyzzz = cbuffer.data(fi_geom_10_off + 489 * ccomps * dcomps);

            auto g_y_0_yyz_xxzzzz = cbuffer.data(fi_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyyy = cbuffer.data(fi_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyyz = cbuffer.data(fi_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_yyz_xyyyzz = cbuffer.data(fi_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_yyz_xyyzzz = cbuffer.data(fi_geom_10_off + 494 * ccomps * dcomps);

            auto g_y_0_yyz_xyzzzz = cbuffer.data(fi_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_yyz_xzzzzz = cbuffer.data(fi_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyyy = cbuffer.data(fi_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyyz = cbuffer.data(fi_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_yyz_yyyyzz = cbuffer.data(fi_geom_10_off + 499 * ccomps * dcomps);

            auto g_y_0_yyz_yyyzzz = cbuffer.data(fi_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_yyz_yyzzzz = cbuffer.data(fi_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_yyz_yzzzzz = cbuffer.data(fi_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_yyz_zzzzzz = cbuffer.data(fi_geom_10_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yy_xxxxxx, g_y_0_yy_xxxxxxz, g_y_0_yy_xxxxxy, g_y_0_yy_xxxxxyz, g_y_0_yy_xxxxxz, g_y_0_yy_xxxxxzz, g_y_0_yy_xxxxyy, g_y_0_yy_xxxxyyz, g_y_0_yy_xxxxyz, g_y_0_yy_xxxxyzz, g_y_0_yy_xxxxzz, g_y_0_yy_xxxxzzz, g_y_0_yy_xxxyyy, g_y_0_yy_xxxyyyz, g_y_0_yy_xxxyyz, g_y_0_yy_xxxyyzz, g_y_0_yy_xxxyzz, g_y_0_yy_xxxyzzz, g_y_0_yy_xxxzzz, g_y_0_yy_xxxzzzz, g_y_0_yy_xxyyyy, g_y_0_yy_xxyyyyz, g_y_0_yy_xxyyyz, g_y_0_yy_xxyyyzz, g_y_0_yy_xxyyzz, g_y_0_yy_xxyyzzz, g_y_0_yy_xxyzzz, g_y_0_yy_xxyzzzz, g_y_0_yy_xxzzzz, g_y_0_yy_xxzzzzz, g_y_0_yy_xyyyyy, g_y_0_yy_xyyyyyz, g_y_0_yy_xyyyyz, g_y_0_yy_xyyyyzz, g_y_0_yy_xyyyzz, g_y_0_yy_xyyyzzz, g_y_0_yy_xyyzzz, g_y_0_yy_xyyzzzz, g_y_0_yy_xyzzzz, g_y_0_yy_xyzzzzz, g_y_0_yy_xzzzzz, g_y_0_yy_xzzzzzz, g_y_0_yy_yyyyyy, g_y_0_yy_yyyyyyz, g_y_0_yy_yyyyyz, g_y_0_yy_yyyyyzz, g_y_0_yy_yyyyzz, g_y_0_yy_yyyyzzz, g_y_0_yy_yyyzzz, g_y_0_yy_yyyzzzz, g_y_0_yy_yyzzzz, g_y_0_yy_yyzzzzz, g_y_0_yy_yzzzzz, g_y_0_yy_yzzzzzz, g_y_0_yy_zzzzzz, g_y_0_yy_zzzzzzz, g_y_0_yyz_xxxxxx, g_y_0_yyz_xxxxxy, g_y_0_yyz_xxxxxz, g_y_0_yyz_xxxxyy, g_y_0_yyz_xxxxyz, g_y_0_yyz_xxxxzz, g_y_0_yyz_xxxyyy, g_y_0_yyz_xxxyyz, g_y_0_yyz_xxxyzz, g_y_0_yyz_xxxzzz, g_y_0_yyz_xxyyyy, g_y_0_yyz_xxyyyz, g_y_0_yyz_xxyyzz, g_y_0_yyz_xxyzzz, g_y_0_yyz_xxzzzz, g_y_0_yyz_xyyyyy, g_y_0_yyz_xyyyyz, g_y_0_yyz_xyyyzz, g_y_0_yyz_xyyzzz, g_y_0_yyz_xyzzzz, g_y_0_yyz_xzzzzz, g_y_0_yyz_yyyyyy, g_y_0_yyz_yyyyyz, g_y_0_yyz_yyyyzz, g_y_0_yyz_yyyzzz, g_y_0_yyz_yyzzzz, g_y_0_yyz_yzzzzz, g_y_0_yyz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyz_xxxxxx[k] = -g_y_0_yy_xxxxxx[k] * ab_z + g_y_0_yy_xxxxxxz[k];

                g_y_0_yyz_xxxxxy[k] = -g_y_0_yy_xxxxxy[k] * ab_z + g_y_0_yy_xxxxxyz[k];

                g_y_0_yyz_xxxxxz[k] = -g_y_0_yy_xxxxxz[k] * ab_z + g_y_0_yy_xxxxxzz[k];

                g_y_0_yyz_xxxxyy[k] = -g_y_0_yy_xxxxyy[k] * ab_z + g_y_0_yy_xxxxyyz[k];

                g_y_0_yyz_xxxxyz[k] = -g_y_0_yy_xxxxyz[k] * ab_z + g_y_0_yy_xxxxyzz[k];

                g_y_0_yyz_xxxxzz[k] = -g_y_0_yy_xxxxzz[k] * ab_z + g_y_0_yy_xxxxzzz[k];

                g_y_0_yyz_xxxyyy[k] = -g_y_0_yy_xxxyyy[k] * ab_z + g_y_0_yy_xxxyyyz[k];

                g_y_0_yyz_xxxyyz[k] = -g_y_0_yy_xxxyyz[k] * ab_z + g_y_0_yy_xxxyyzz[k];

                g_y_0_yyz_xxxyzz[k] = -g_y_0_yy_xxxyzz[k] * ab_z + g_y_0_yy_xxxyzzz[k];

                g_y_0_yyz_xxxzzz[k] = -g_y_0_yy_xxxzzz[k] * ab_z + g_y_0_yy_xxxzzzz[k];

                g_y_0_yyz_xxyyyy[k] = -g_y_0_yy_xxyyyy[k] * ab_z + g_y_0_yy_xxyyyyz[k];

                g_y_0_yyz_xxyyyz[k] = -g_y_0_yy_xxyyyz[k] * ab_z + g_y_0_yy_xxyyyzz[k];

                g_y_0_yyz_xxyyzz[k] = -g_y_0_yy_xxyyzz[k] * ab_z + g_y_0_yy_xxyyzzz[k];

                g_y_0_yyz_xxyzzz[k] = -g_y_0_yy_xxyzzz[k] * ab_z + g_y_0_yy_xxyzzzz[k];

                g_y_0_yyz_xxzzzz[k] = -g_y_0_yy_xxzzzz[k] * ab_z + g_y_0_yy_xxzzzzz[k];

                g_y_0_yyz_xyyyyy[k] = -g_y_0_yy_xyyyyy[k] * ab_z + g_y_0_yy_xyyyyyz[k];

                g_y_0_yyz_xyyyyz[k] = -g_y_0_yy_xyyyyz[k] * ab_z + g_y_0_yy_xyyyyzz[k];

                g_y_0_yyz_xyyyzz[k] = -g_y_0_yy_xyyyzz[k] * ab_z + g_y_0_yy_xyyyzzz[k];

                g_y_0_yyz_xyyzzz[k] = -g_y_0_yy_xyyzzz[k] * ab_z + g_y_0_yy_xyyzzzz[k];

                g_y_0_yyz_xyzzzz[k] = -g_y_0_yy_xyzzzz[k] * ab_z + g_y_0_yy_xyzzzzz[k];

                g_y_0_yyz_xzzzzz[k] = -g_y_0_yy_xzzzzz[k] * ab_z + g_y_0_yy_xzzzzzz[k];

                g_y_0_yyz_yyyyyy[k] = -g_y_0_yy_yyyyyy[k] * ab_z + g_y_0_yy_yyyyyyz[k];

                g_y_0_yyz_yyyyyz[k] = -g_y_0_yy_yyyyyz[k] * ab_z + g_y_0_yy_yyyyyzz[k];

                g_y_0_yyz_yyyyzz[k] = -g_y_0_yy_yyyyzz[k] * ab_z + g_y_0_yy_yyyyzzz[k];

                g_y_0_yyz_yyyzzz[k] = -g_y_0_yy_yyyzzz[k] * ab_z + g_y_0_yy_yyyzzzz[k];

                g_y_0_yyz_yyzzzz[k] = -g_y_0_yy_yyzzzz[k] * ab_z + g_y_0_yy_yyzzzzz[k];

                g_y_0_yyz_yzzzzz[k] = -g_y_0_yy_yzzzzz[k] * ab_z + g_y_0_yy_yzzzzzz[k];

                g_y_0_yyz_zzzzzz[k] = -g_y_0_yy_zzzzzz[k] * ab_z + g_y_0_yy_zzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzz_xxxxxx = cbuffer.data(fi_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxy = cbuffer.data(fi_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxxz = cbuffer.data(fi_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxyy = cbuffer.data(fi_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxyz = cbuffer.data(fi_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_yzz_xxxxzz = cbuffer.data(fi_geom_10_off + 509 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyyy = cbuffer.data(fi_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyyz = cbuffer.data(fi_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_yzz_xxxyzz = cbuffer.data(fi_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_yzz_xxxzzz = cbuffer.data(fi_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyyy = cbuffer.data(fi_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyyz = cbuffer.data(fi_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_yzz_xxyyzz = cbuffer.data(fi_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_yzz_xxyzzz = cbuffer.data(fi_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_yzz_xxzzzz = cbuffer.data(fi_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyyy = cbuffer.data(fi_geom_10_off + 519 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyyz = cbuffer.data(fi_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_yzz_xyyyzz = cbuffer.data(fi_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_yzz_xyyzzz = cbuffer.data(fi_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_yzz_xyzzzz = cbuffer.data(fi_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_yzz_xzzzzz = cbuffer.data(fi_geom_10_off + 524 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyyy = cbuffer.data(fi_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyyz = cbuffer.data(fi_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_yzz_yyyyzz = cbuffer.data(fi_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_yzz_yyyzzz = cbuffer.data(fi_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_yzz_yyzzzz = cbuffer.data(fi_geom_10_off + 529 * ccomps * dcomps);

            auto g_y_0_yzz_yzzzzz = cbuffer.data(fi_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_yzz_zzzzzz = cbuffer.data(fi_geom_10_off + 531 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yz_xxxxxx, g_y_0_yz_xxxxxxz, g_y_0_yz_xxxxxy, g_y_0_yz_xxxxxyz, g_y_0_yz_xxxxxz, g_y_0_yz_xxxxxzz, g_y_0_yz_xxxxyy, g_y_0_yz_xxxxyyz, g_y_0_yz_xxxxyz, g_y_0_yz_xxxxyzz, g_y_0_yz_xxxxzz, g_y_0_yz_xxxxzzz, g_y_0_yz_xxxyyy, g_y_0_yz_xxxyyyz, g_y_0_yz_xxxyyz, g_y_0_yz_xxxyyzz, g_y_0_yz_xxxyzz, g_y_0_yz_xxxyzzz, g_y_0_yz_xxxzzz, g_y_0_yz_xxxzzzz, g_y_0_yz_xxyyyy, g_y_0_yz_xxyyyyz, g_y_0_yz_xxyyyz, g_y_0_yz_xxyyyzz, g_y_0_yz_xxyyzz, g_y_0_yz_xxyyzzz, g_y_0_yz_xxyzzz, g_y_0_yz_xxyzzzz, g_y_0_yz_xxzzzz, g_y_0_yz_xxzzzzz, g_y_0_yz_xyyyyy, g_y_0_yz_xyyyyyz, g_y_0_yz_xyyyyz, g_y_0_yz_xyyyyzz, g_y_0_yz_xyyyzz, g_y_0_yz_xyyyzzz, g_y_0_yz_xyyzzz, g_y_0_yz_xyyzzzz, g_y_0_yz_xyzzzz, g_y_0_yz_xyzzzzz, g_y_0_yz_xzzzzz, g_y_0_yz_xzzzzzz, g_y_0_yz_yyyyyy, g_y_0_yz_yyyyyyz, g_y_0_yz_yyyyyz, g_y_0_yz_yyyyyzz, g_y_0_yz_yyyyzz, g_y_0_yz_yyyyzzz, g_y_0_yz_yyyzzz, g_y_0_yz_yyyzzzz, g_y_0_yz_yyzzzz, g_y_0_yz_yyzzzzz, g_y_0_yz_yzzzzz, g_y_0_yz_yzzzzzz, g_y_0_yz_zzzzzz, g_y_0_yz_zzzzzzz, g_y_0_yzz_xxxxxx, g_y_0_yzz_xxxxxy, g_y_0_yzz_xxxxxz, g_y_0_yzz_xxxxyy, g_y_0_yzz_xxxxyz, g_y_0_yzz_xxxxzz, g_y_0_yzz_xxxyyy, g_y_0_yzz_xxxyyz, g_y_0_yzz_xxxyzz, g_y_0_yzz_xxxzzz, g_y_0_yzz_xxyyyy, g_y_0_yzz_xxyyyz, g_y_0_yzz_xxyyzz, g_y_0_yzz_xxyzzz, g_y_0_yzz_xxzzzz, g_y_0_yzz_xyyyyy, g_y_0_yzz_xyyyyz, g_y_0_yzz_xyyyzz, g_y_0_yzz_xyyzzz, g_y_0_yzz_xyzzzz, g_y_0_yzz_xzzzzz, g_y_0_yzz_yyyyyy, g_y_0_yzz_yyyyyz, g_y_0_yzz_yyyyzz, g_y_0_yzz_yyyzzz, g_y_0_yzz_yyzzzz, g_y_0_yzz_yzzzzz, g_y_0_yzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzz_xxxxxx[k] = -g_y_0_yz_xxxxxx[k] * ab_z + g_y_0_yz_xxxxxxz[k];

                g_y_0_yzz_xxxxxy[k] = -g_y_0_yz_xxxxxy[k] * ab_z + g_y_0_yz_xxxxxyz[k];

                g_y_0_yzz_xxxxxz[k] = -g_y_0_yz_xxxxxz[k] * ab_z + g_y_0_yz_xxxxxzz[k];

                g_y_0_yzz_xxxxyy[k] = -g_y_0_yz_xxxxyy[k] * ab_z + g_y_0_yz_xxxxyyz[k];

                g_y_0_yzz_xxxxyz[k] = -g_y_0_yz_xxxxyz[k] * ab_z + g_y_0_yz_xxxxyzz[k];

                g_y_0_yzz_xxxxzz[k] = -g_y_0_yz_xxxxzz[k] * ab_z + g_y_0_yz_xxxxzzz[k];

                g_y_0_yzz_xxxyyy[k] = -g_y_0_yz_xxxyyy[k] * ab_z + g_y_0_yz_xxxyyyz[k];

                g_y_0_yzz_xxxyyz[k] = -g_y_0_yz_xxxyyz[k] * ab_z + g_y_0_yz_xxxyyzz[k];

                g_y_0_yzz_xxxyzz[k] = -g_y_0_yz_xxxyzz[k] * ab_z + g_y_0_yz_xxxyzzz[k];

                g_y_0_yzz_xxxzzz[k] = -g_y_0_yz_xxxzzz[k] * ab_z + g_y_0_yz_xxxzzzz[k];

                g_y_0_yzz_xxyyyy[k] = -g_y_0_yz_xxyyyy[k] * ab_z + g_y_0_yz_xxyyyyz[k];

                g_y_0_yzz_xxyyyz[k] = -g_y_0_yz_xxyyyz[k] * ab_z + g_y_0_yz_xxyyyzz[k];

                g_y_0_yzz_xxyyzz[k] = -g_y_0_yz_xxyyzz[k] * ab_z + g_y_0_yz_xxyyzzz[k];

                g_y_0_yzz_xxyzzz[k] = -g_y_0_yz_xxyzzz[k] * ab_z + g_y_0_yz_xxyzzzz[k];

                g_y_0_yzz_xxzzzz[k] = -g_y_0_yz_xxzzzz[k] * ab_z + g_y_0_yz_xxzzzzz[k];

                g_y_0_yzz_xyyyyy[k] = -g_y_0_yz_xyyyyy[k] * ab_z + g_y_0_yz_xyyyyyz[k];

                g_y_0_yzz_xyyyyz[k] = -g_y_0_yz_xyyyyz[k] * ab_z + g_y_0_yz_xyyyyzz[k];

                g_y_0_yzz_xyyyzz[k] = -g_y_0_yz_xyyyzz[k] * ab_z + g_y_0_yz_xyyyzzz[k];

                g_y_0_yzz_xyyzzz[k] = -g_y_0_yz_xyyzzz[k] * ab_z + g_y_0_yz_xyyzzzz[k];

                g_y_0_yzz_xyzzzz[k] = -g_y_0_yz_xyzzzz[k] * ab_z + g_y_0_yz_xyzzzzz[k];

                g_y_0_yzz_xzzzzz[k] = -g_y_0_yz_xzzzzz[k] * ab_z + g_y_0_yz_xzzzzzz[k];

                g_y_0_yzz_yyyyyy[k] = -g_y_0_yz_yyyyyy[k] * ab_z + g_y_0_yz_yyyyyyz[k];

                g_y_0_yzz_yyyyyz[k] = -g_y_0_yz_yyyyyz[k] * ab_z + g_y_0_yz_yyyyyzz[k];

                g_y_0_yzz_yyyyzz[k] = -g_y_0_yz_yyyyzz[k] * ab_z + g_y_0_yz_yyyyzzz[k];

                g_y_0_yzz_yyyzzz[k] = -g_y_0_yz_yyyzzz[k] * ab_z + g_y_0_yz_yyyzzzz[k];

                g_y_0_yzz_yyzzzz[k] = -g_y_0_yz_yyzzzz[k] * ab_z + g_y_0_yz_yyzzzzz[k];

                g_y_0_yzz_yzzzzz[k] = -g_y_0_yz_yzzzzz[k] * ab_z + g_y_0_yz_yzzzzzz[k];

                g_y_0_yzz_zzzzzz[k] = -g_y_0_yz_zzzzzz[k] * ab_z + g_y_0_yz_zzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzz_xxxxxx = cbuffer.data(fi_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxy = cbuffer.data(fi_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxxz = cbuffer.data(fi_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxyy = cbuffer.data(fi_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxyz = cbuffer.data(fi_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_zzz_xxxxzz = cbuffer.data(fi_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyyy = cbuffer.data(fi_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyyz = cbuffer.data(fi_geom_10_off + 539 * ccomps * dcomps);

            auto g_y_0_zzz_xxxyzz = cbuffer.data(fi_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_zzz_xxxzzz = cbuffer.data(fi_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyyy = cbuffer.data(fi_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyyz = cbuffer.data(fi_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_zzz_xxyyzz = cbuffer.data(fi_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_zzz_xxyzzz = cbuffer.data(fi_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_zzz_xxzzzz = cbuffer.data(fi_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyyy = cbuffer.data(fi_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyyz = cbuffer.data(fi_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_zzz_xyyyzz = cbuffer.data(fi_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_zzz_xyyzzz = cbuffer.data(fi_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_zzz_xyzzzz = cbuffer.data(fi_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_zzz_xzzzzz = cbuffer.data(fi_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyyy = cbuffer.data(fi_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyyz = cbuffer.data(fi_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_zzz_yyyyzz = cbuffer.data(fi_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_zzz_yyyzzz = cbuffer.data(fi_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_zzz_yyzzzz = cbuffer.data(fi_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_zzz_yzzzzz = cbuffer.data(fi_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_zzz_zzzzzz = cbuffer.data(fi_geom_10_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zz_xxxxxx, g_y_0_zz_xxxxxxz, g_y_0_zz_xxxxxy, g_y_0_zz_xxxxxyz, g_y_0_zz_xxxxxz, g_y_0_zz_xxxxxzz, g_y_0_zz_xxxxyy, g_y_0_zz_xxxxyyz, g_y_0_zz_xxxxyz, g_y_0_zz_xxxxyzz, g_y_0_zz_xxxxzz, g_y_0_zz_xxxxzzz, g_y_0_zz_xxxyyy, g_y_0_zz_xxxyyyz, g_y_0_zz_xxxyyz, g_y_0_zz_xxxyyzz, g_y_0_zz_xxxyzz, g_y_0_zz_xxxyzzz, g_y_0_zz_xxxzzz, g_y_0_zz_xxxzzzz, g_y_0_zz_xxyyyy, g_y_0_zz_xxyyyyz, g_y_0_zz_xxyyyz, g_y_0_zz_xxyyyzz, g_y_0_zz_xxyyzz, g_y_0_zz_xxyyzzz, g_y_0_zz_xxyzzz, g_y_0_zz_xxyzzzz, g_y_0_zz_xxzzzz, g_y_0_zz_xxzzzzz, g_y_0_zz_xyyyyy, g_y_0_zz_xyyyyyz, g_y_0_zz_xyyyyz, g_y_0_zz_xyyyyzz, g_y_0_zz_xyyyzz, g_y_0_zz_xyyyzzz, g_y_0_zz_xyyzzz, g_y_0_zz_xyyzzzz, g_y_0_zz_xyzzzz, g_y_0_zz_xyzzzzz, g_y_0_zz_xzzzzz, g_y_0_zz_xzzzzzz, g_y_0_zz_yyyyyy, g_y_0_zz_yyyyyyz, g_y_0_zz_yyyyyz, g_y_0_zz_yyyyyzz, g_y_0_zz_yyyyzz, g_y_0_zz_yyyyzzz, g_y_0_zz_yyyzzz, g_y_0_zz_yyyzzzz, g_y_0_zz_yyzzzz, g_y_0_zz_yyzzzzz, g_y_0_zz_yzzzzz, g_y_0_zz_yzzzzzz, g_y_0_zz_zzzzzz, g_y_0_zz_zzzzzzz, g_y_0_zzz_xxxxxx, g_y_0_zzz_xxxxxy, g_y_0_zzz_xxxxxz, g_y_0_zzz_xxxxyy, g_y_0_zzz_xxxxyz, g_y_0_zzz_xxxxzz, g_y_0_zzz_xxxyyy, g_y_0_zzz_xxxyyz, g_y_0_zzz_xxxyzz, g_y_0_zzz_xxxzzz, g_y_0_zzz_xxyyyy, g_y_0_zzz_xxyyyz, g_y_0_zzz_xxyyzz, g_y_0_zzz_xxyzzz, g_y_0_zzz_xxzzzz, g_y_0_zzz_xyyyyy, g_y_0_zzz_xyyyyz, g_y_0_zzz_xyyyzz, g_y_0_zzz_xyyzzz, g_y_0_zzz_xyzzzz, g_y_0_zzz_xzzzzz, g_y_0_zzz_yyyyyy, g_y_0_zzz_yyyyyz, g_y_0_zzz_yyyyzz, g_y_0_zzz_yyyzzz, g_y_0_zzz_yyzzzz, g_y_0_zzz_yzzzzz, g_y_0_zzz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzz_xxxxxx[k] = -g_y_0_zz_xxxxxx[k] * ab_z + g_y_0_zz_xxxxxxz[k];

                g_y_0_zzz_xxxxxy[k] = -g_y_0_zz_xxxxxy[k] * ab_z + g_y_0_zz_xxxxxyz[k];

                g_y_0_zzz_xxxxxz[k] = -g_y_0_zz_xxxxxz[k] * ab_z + g_y_0_zz_xxxxxzz[k];

                g_y_0_zzz_xxxxyy[k] = -g_y_0_zz_xxxxyy[k] * ab_z + g_y_0_zz_xxxxyyz[k];

                g_y_0_zzz_xxxxyz[k] = -g_y_0_zz_xxxxyz[k] * ab_z + g_y_0_zz_xxxxyzz[k];

                g_y_0_zzz_xxxxzz[k] = -g_y_0_zz_xxxxzz[k] * ab_z + g_y_0_zz_xxxxzzz[k];

                g_y_0_zzz_xxxyyy[k] = -g_y_0_zz_xxxyyy[k] * ab_z + g_y_0_zz_xxxyyyz[k];

                g_y_0_zzz_xxxyyz[k] = -g_y_0_zz_xxxyyz[k] * ab_z + g_y_0_zz_xxxyyzz[k];

                g_y_0_zzz_xxxyzz[k] = -g_y_0_zz_xxxyzz[k] * ab_z + g_y_0_zz_xxxyzzz[k];

                g_y_0_zzz_xxxzzz[k] = -g_y_0_zz_xxxzzz[k] * ab_z + g_y_0_zz_xxxzzzz[k];

                g_y_0_zzz_xxyyyy[k] = -g_y_0_zz_xxyyyy[k] * ab_z + g_y_0_zz_xxyyyyz[k];

                g_y_0_zzz_xxyyyz[k] = -g_y_0_zz_xxyyyz[k] * ab_z + g_y_0_zz_xxyyyzz[k];

                g_y_0_zzz_xxyyzz[k] = -g_y_0_zz_xxyyzz[k] * ab_z + g_y_0_zz_xxyyzzz[k];

                g_y_0_zzz_xxyzzz[k] = -g_y_0_zz_xxyzzz[k] * ab_z + g_y_0_zz_xxyzzzz[k];

                g_y_0_zzz_xxzzzz[k] = -g_y_0_zz_xxzzzz[k] * ab_z + g_y_0_zz_xxzzzzz[k];

                g_y_0_zzz_xyyyyy[k] = -g_y_0_zz_xyyyyy[k] * ab_z + g_y_0_zz_xyyyyyz[k];

                g_y_0_zzz_xyyyyz[k] = -g_y_0_zz_xyyyyz[k] * ab_z + g_y_0_zz_xyyyyzz[k];

                g_y_0_zzz_xyyyzz[k] = -g_y_0_zz_xyyyzz[k] * ab_z + g_y_0_zz_xyyyzzz[k];

                g_y_0_zzz_xyyzzz[k] = -g_y_0_zz_xyyzzz[k] * ab_z + g_y_0_zz_xyyzzzz[k];

                g_y_0_zzz_xyzzzz[k] = -g_y_0_zz_xyzzzz[k] * ab_z + g_y_0_zz_xyzzzzz[k];

                g_y_0_zzz_xzzzzz[k] = -g_y_0_zz_xzzzzz[k] * ab_z + g_y_0_zz_xzzzzzz[k];

                g_y_0_zzz_yyyyyy[k] = -g_y_0_zz_yyyyyy[k] * ab_z + g_y_0_zz_yyyyyyz[k];

                g_y_0_zzz_yyyyyz[k] = -g_y_0_zz_yyyyyz[k] * ab_z + g_y_0_zz_yyyyyzz[k];

                g_y_0_zzz_yyyyzz[k] = -g_y_0_zz_yyyyzz[k] * ab_z + g_y_0_zz_yyyyzzz[k];

                g_y_0_zzz_yyyzzz[k] = -g_y_0_zz_yyyzzz[k] * ab_z + g_y_0_zz_yyyzzzz[k];

                g_y_0_zzz_yyzzzz[k] = -g_y_0_zz_yyzzzz[k] * ab_z + g_y_0_zz_yyzzzzz[k];

                g_y_0_zzz_yzzzzz[k] = -g_y_0_zz_yzzzzz[k] * ab_z + g_y_0_zz_yzzzzzz[k];

                g_y_0_zzz_zzzzzz[k] = -g_y_0_zz_zzzzzz[k] * ab_z + g_y_0_zz_zzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxx_xxxxxx = cbuffer.data(fi_geom_10_off + 560 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxy = cbuffer.data(fi_geom_10_off + 561 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxxz = cbuffer.data(fi_geom_10_off + 562 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxyy = cbuffer.data(fi_geom_10_off + 563 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxyz = cbuffer.data(fi_geom_10_off + 564 * ccomps * dcomps);

            auto g_z_0_xxx_xxxxzz = cbuffer.data(fi_geom_10_off + 565 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyyy = cbuffer.data(fi_geom_10_off + 566 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyyz = cbuffer.data(fi_geom_10_off + 567 * ccomps * dcomps);

            auto g_z_0_xxx_xxxyzz = cbuffer.data(fi_geom_10_off + 568 * ccomps * dcomps);

            auto g_z_0_xxx_xxxzzz = cbuffer.data(fi_geom_10_off + 569 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyyy = cbuffer.data(fi_geom_10_off + 570 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyyz = cbuffer.data(fi_geom_10_off + 571 * ccomps * dcomps);

            auto g_z_0_xxx_xxyyzz = cbuffer.data(fi_geom_10_off + 572 * ccomps * dcomps);

            auto g_z_0_xxx_xxyzzz = cbuffer.data(fi_geom_10_off + 573 * ccomps * dcomps);

            auto g_z_0_xxx_xxzzzz = cbuffer.data(fi_geom_10_off + 574 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyyy = cbuffer.data(fi_geom_10_off + 575 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyyz = cbuffer.data(fi_geom_10_off + 576 * ccomps * dcomps);

            auto g_z_0_xxx_xyyyzz = cbuffer.data(fi_geom_10_off + 577 * ccomps * dcomps);

            auto g_z_0_xxx_xyyzzz = cbuffer.data(fi_geom_10_off + 578 * ccomps * dcomps);

            auto g_z_0_xxx_xyzzzz = cbuffer.data(fi_geom_10_off + 579 * ccomps * dcomps);

            auto g_z_0_xxx_xzzzzz = cbuffer.data(fi_geom_10_off + 580 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyyy = cbuffer.data(fi_geom_10_off + 581 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyyz = cbuffer.data(fi_geom_10_off + 582 * ccomps * dcomps);

            auto g_z_0_xxx_yyyyzz = cbuffer.data(fi_geom_10_off + 583 * ccomps * dcomps);

            auto g_z_0_xxx_yyyzzz = cbuffer.data(fi_geom_10_off + 584 * ccomps * dcomps);

            auto g_z_0_xxx_yyzzzz = cbuffer.data(fi_geom_10_off + 585 * ccomps * dcomps);

            auto g_z_0_xxx_yzzzzz = cbuffer.data(fi_geom_10_off + 586 * ccomps * dcomps);

            auto g_z_0_xxx_zzzzzz = cbuffer.data(fi_geom_10_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xx_xxxxxx, g_z_0_xx_xxxxxxx, g_z_0_xx_xxxxxxy, g_z_0_xx_xxxxxxz, g_z_0_xx_xxxxxy, g_z_0_xx_xxxxxyy, g_z_0_xx_xxxxxyz, g_z_0_xx_xxxxxz, g_z_0_xx_xxxxxzz, g_z_0_xx_xxxxyy, g_z_0_xx_xxxxyyy, g_z_0_xx_xxxxyyz, g_z_0_xx_xxxxyz, g_z_0_xx_xxxxyzz, g_z_0_xx_xxxxzz, g_z_0_xx_xxxxzzz, g_z_0_xx_xxxyyy, g_z_0_xx_xxxyyyy, g_z_0_xx_xxxyyyz, g_z_0_xx_xxxyyz, g_z_0_xx_xxxyyzz, g_z_0_xx_xxxyzz, g_z_0_xx_xxxyzzz, g_z_0_xx_xxxzzz, g_z_0_xx_xxxzzzz, g_z_0_xx_xxyyyy, g_z_0_xx_xxyyyyy, g_z_0_xx_xxyyyyz, g_z_0_xx_xxyyyz, g_z_0_xx_xxyyyzz, g_z_0_xx_xxyyzz, g_z_0_xx_xxyyzzz, g_z_0_xx_xxyzzz, g_z_0_xx_xxyzzzz, g_z_0_xx_xxzzzz, g_z_0_xx_xxzzzzz, g_z_0_xx_xyyyyy, g_z_0_xx_xyyyyyy, g_z_0_xx_xyyyyyz, g_z_0_xx_xyyyyz, g_z_0_xx_xyyyyzz, g_z_0_xx_xyyyzz, g_z_0_xx_xyyyzzz, g_z_0_xx_xyyzzz, g_z_0_xx_xyyzzzz, g_z_0_xx_xyzzzz, g_z_0_xx_xyzzzzz, g_z_0_xx_xzzzzz, g_z_0_xx_xzzzzzz, g_z_0_xx_yyyyyy, g_z_0_xx_yyyyyz, g_z_0_xx_yyyyzz, g_z_0_xx_yyyzzz, g_z_0_xx_yyzzzz, g_z_0_xx_yzzzzz, g_z_0_xx_zzzzzz, g_z_0_xxx_xxxxxx, g_z_0_xxx_xxxxxy, g_z_0_xxx_xxxxxz, g_z_0_xxx_xxxxyy, g_z_0_xxx_xxxxyz, g_z_0_xxx_xxxxzz, g_z_0_xxx_xxxyyy, g_z_0_xxx_xxxyyz, g_z_0_xxx_xxxyzz, g_z_0_xxx_xxxzzz, g_z_0_xxx_xxyyyy, g_z_0_xxx_xxyyyz, g_z_0_xxx_xxyyzz, g_z_0_xxx_xxyzzz, g_z_0_xxx_xxzzzz, g_z_0_xxx_xyyyyy, g_z_0_xxx_xyyyyz, g_z_0_xxx_xyyyzz, g_z_0_xxx_xyyzzz, g_z_0_xxx_xyzzzz, g_z_0_xxx_xzzzzz, g_z_0_xxx_yyyyyy, g_z_0_xxx_yyyyyz, g_z_0_xxx_yyyyzz, g_z_0_xxx_yyyzzz, g_z_0_xxx_yyzzzz, g_z_0_xxx_yzzzzz, g_z_0_xxx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxx_xxxxxx[k] = -g_z_0_xx_xxxxxx[k] * ab_x + g_z_0_xx_xxxxxxx[k];

                g_z_0_xxx_xxxxxy[k] = -g_z_0_xx_xxxxxy[k] * ab_x + g_z_0_xx_xxxxxxy[k];

                g_z_0_xxx_xxxxxz[k] = -g_z_0_xx_xxxxxz[k] * ab_x + g_z_0_xx_xxxxxxz[k];

                g_z_0_xxx_xxxxyy[k] = -g_z_0_xx_xxxxyy[k] * ab_x + g_z_0_xx_xxxxxyy[k];

                g_z_0_xxx_xxxxyz[k] = -g_z_0_xx_xxxxyz[k] * ab_x + g_z_0_xx_xxxxxyz[k];

                g_z_0_xxx_xxxxzz[k] = -g_z_0_xx_xxxxzz[k] * ab_x + g_z_0_xx_xxxxxzz[k];

                g_z_0_xxx_xxxyyy[k] = -g_z_0_xx_xxxyyy[k] * ab_x + g_z_0_xx_xxxxyyy[k];

                g_z_0_xxx_xxxyyz[k] = -g_z_0_xx_xxxyyz[k] * ab_x + g_z_0_xx_xxxxyyz[k];

                g_z_0_xxx_xxxyzz[k] = -g_z_0_xx_xxxyzz[k] * ab_x + g_z_0_xx_xxxxyzz[k];

                g_z_0_xxx_xxxzzz[k] = -g_z_0_xx_xxxzzz[k] * ab_x + g_z_0_xx_xxxxzzz[k];

                g_z_0_xxx_xxyyyy[k] = -g_z_0_xx_xxyyyy[k] * ab_x + g_z_0_xx_xxxyyyy[k];

                g_z_0_xxx_xxyyyz[k] = -g_z_0_xx_xxyyyz[k] * ab_x + g_z_0_xx_xxxyyyz[k];

                g_z_0_xxx_xxyyzz[k] = -g_z_0_xx_xxyyzz[k] * ab_x + g_z_0_xx_xxxyyzz[k];

                g_z_0_xxx_xxyzzz[k] = -g_z_0_xx_xxyzzz[k] * ab_x + g_z_0_xx_xxxyzzz[k];

                g_z_0_xxx_xxzzzz[k] = -g_z_0_xx_xxzzzz[k] * ab_x + g_z_0_xx_xxxzzzz[k];

                g_z_0_xxx_xyyyyy[k] = -g_z_0_xx_xyyyyy[k] * ab_x + g_z_0_xx_xxyyyyy[k];

                g_z_0_xxx_xyyyyz[k] = -g_z_0_xx_xyyyyz[k] * ab_x + g_z_0_xx_xxyyyyz[k];

                g_z_0_xxx_xyyyzz[k] = -g_z_0_xx_xyyyzz[k] * ab_x + g_z_0_xx_xxyyyzz[k];

                g_z_0_xxx_xyyzzz[k] = -g_z_0_xx_xyyzzz[k] * ab_x + g_z_0_xx_xxyyzzz[k];

                g_z_0_xxx_xyzzzz[k] = -g_z_0_xx_xyzzzz[k] * ab_x + g_z_0_xx_xxyzzzz[k];

                g_z_0_xxx_xzzzzz[k] = -g_z_0_xx_xzzzzz[k] * ab_x + g_z_0_xx_xxzzzzz[k];

                g_z_0_xxx_yyyyyy[k] = -g_z_0_xx_yyyyyy[k] * ab_x + g_z_0_xx_xyyyyyy[k];

                g_z_0_xxx_yyyyyz[k] = -g_z_0_xx_yyyyyz[k] * ab_x + g_z_0_xx_xyyyyyz[k];

                g_z_0_xxx_yyyyzz[k] = -g_z_0_xx_yyyyzz[k] * ab_x + g_z_0_xx_xyyyyzz[k];

                g_z_0_xxx_yyyzzz[k] = -g_z_0_xx_yyyzzz[k] * ab_x + g_z_0_xx_xyyyzzz[k];

                g_z_0_xxx_yyzzzz[k] = -g_z_0_xx_yyzzzz[k] * ab_x + g_z_0_xx_xyyzzzz[k];

                g_z_0_xxx_yzzzzz[k] = -g_z_0_xx_yzzzzz[k] * ab_x + g_z_0_xx_xyzzzzz[k];

                g_z_0_xxx_zzzzzz[k] = -g_z_0_xx_zzzzzz[k] * ab_x + g_z_0_xx_xzzzzzz[k];
            }

            /// Set up 588-616 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxy_xxxxxx = cbuffer.data(fi_geom_10_off + 588 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxy = cbuffer.data(fi_geom_10_off + 589 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxxz = cbuffer.data(fi_geom_10_off + 590 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxyy = cbuffer.data(fi_geom_10_off + 591 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxyz = cbuffer.data(fi_geom_10_off + 592 * ccomps * dcomps);

            auto g_z_0_xxy_xxxxzz = cbuffer.data(fi_geom_10_off + 593 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyyy = cbuffer.data(fi_geom_10_off + 594 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyyz = cbuffer.data(fi_geom_10_off + 595 * ccomps * dcomps);

            auto g_z_0_xxy_xxxyzz = cbuffer.data(fi_geom_10_off + 596 * ccomps * dcomps);

            auto g_z_0_xxy_xxxzzz = cbuffer.data(fi_geom_10_off + 597 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyyy = cbuffer.data(fi_geom_10_off + 598 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyyz = cbuffer.data(fi_geom_10_off + 599 * ccomps * dcomps);

            auto g_z_0_xxy_xxyyzz = cbuffer.data(fi_geom_10_off + 600 * ccomps * dcomps);

            auto g_z_0_xxy_xxyzzz = cbuffer.data(fi_geom_10_off + 601 * ccomps * dcomps);

            auto g_z_0_xxy_xxzzzz = cbuffer.data(fi_geom_10_off + 602 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyyy = cbuffer.data(fi_geom_10_off + 603 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyyz = cbuffer.data(fi_geom_10_off + 604 * ccomps * dcomps);

            auto g_z_0_xxy_xyyyzz = cbuffer.data(fi_geom_10_off + 605 * ccomps * dcomps);

            auto g_z_0_xxy_xyyzzz = cbuffer.data(fi_geom_10_off + 606 * ccomps * dcomps);

            auto g_z_0_xxy_xyzzzz = cbuffer.data(fi_geom_10_off + 607 * ccomps * dcomps);

            auto g_z_0_xxy_xzzzzz = cbuffer.data(fi_geom_10_off + 608 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyyy = cbuffer.data(fi_geom_10_off + 609 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyyz = cbuffer.data(fi_geom_10_off + 610 * ccomps * dcomps);

            auto g_z_0_xxy_yyyyzz = cbuffer.data(fi_geom_10_off + 611 * ccomps * dcomps);

            auto g_z_0_xxy_yyyzzz = cbuffer.data(fi_geom_10_off + 612 * ccomps * dcomps);

            auto g_z_0_xxy_yyzzzz = cbuffer.data(fi_geom_10_off + 613 * ccomps * dcomps);

            auto g_z_0_xxy_yzzzzz = cbuffer.data(fi_geom_10_off + 614 * ccomps * dcomps);

            auto g_z_0_xxy_zzzzzz = cbuffer.data(fi_geom_10_off + 615 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxy_xxxxxx, g_z_0_xxy_xxxxxy, g_z_0_xxy_xxxxxz, g_z_0_xxy_xxxxyy, g_z_0_xxy_xxxxyz, g_z_0_xxy_xxxxzz, g_z_0_xxy_xxxyyy, g_z_0_xxy_xxxyyz, g_z_0_xxy_xxxyzz, g_z_0_xxy_xxxzzz, g_z_0_xxy_xxyyyy, g_z_0_xxy_xxyyyz, g_z_0_xxy_xxyyzz, g_z_0_xxy_xxyzzz, g_z_0_xxy_xxzzzz, g_z_0_xxy_xyyyyy, g_z_0_xxy_xyyyyz, g_z_0_xxy_xyyyzz, g_z_0_xxy_xyyzzz, g_z_0_xxy_xyzzzz, g_z_0_xxy_xzzzzz, g_z_0_xxy_yyyyyy, g_z_0_xxy_yyyyyz, g_z_0_xxy_yyyyzz, g_z_0_xxy_yyyzzz, g_z_0_xxy_yyzzzz, g_z_0_xxy_yzzzzz, g_z_0_xxy_zzzzzz, g_z_0_xy_xxxxxx, g_z_0_xy_xxxxxxx, g_z_0_xy_xxxxxxy, g_z_0_xy_xxxxxxz, g_z_0_xy_xxxxxy, g_z_0_xy_xxxxxyy, g_z_0_xy_xxxxxyz, g_z_0_xy_xxxxxz, g_z_0_xy_xxxxxzz, g_z_0_xy_xxxxyy, g_z_0_xy_xxxxyyy, g_z_0_xy_xxxxyyz, g_z_0_xy_xxxxyz, g_z_0_xy_xxxxyzz, g_z_0_xy_xxxxzz, g_z_0_xy_xxxxzzz, g_z_0_xy_xxxyyy, g_z_0_xy_xxxyyyy, g_z_0_xy_xxxyyyz, g_z_0_xy_xxxyyz, g_z_0_xy_xxxyyzz, g_z_0_xy_xxxyzz, g_z_0_xy_xxxyzzz, g_z_0_xy_xxxzzz, g_z_0_xy_xxxzzzz, g_z_0_xy_xxyyyy, g_z_0_xy_xxyyyyy, g_z_0_xy_xxyyyyz, g_z_0_xy_xxyyyz, g_z_0_xy_xxyyyzz, g_z_0_xy_xxyyzz, g_z_0_xy_xxyyzzz, g_z_0_xy_xxyzzz, g_z_0_xy_xxyzzzz, g_z_0_xy_xxzzzz, g_z_0_xy_xxzzzzz, g_z_0_xy_xyyyyy, g_z_0_xy_xyyyyyy, g_z_0_xy_xyyyyyz, g_z_0_xy_xyyyyz, g_z_0_xy_xyyyyzz, g_z_0_xy_xyyyzz, g_z_0_xy_xyyyzzz, g_z_0_xy_xyyzzz, g_z_0_xy_xyyzzzz, g_z_0_xy_xyzzzz, g_z_0_xy_xyzzzzz, g_z_0_xy_xzzzzz, g_z_0_xy_xzzzzzz, g_z_0_xy_yyyyyy, g_z_0_xy_yyyyyz, g_z_0_xy_yyyyzz, g_z_0_xy_yyyzzz, g_z_0_xy_yyzzzz, g_z_0_xy_yzzzzz, g_z_0_xy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxy_xxxxxx[k] = -g_z_0_xy_xxxxxx[k] * ab_x + g_z_0_xy_xxxxxxx[k];

                g_z_0_xxy_xxxxxy[k] = -g_z_0_xy_xxxxxy[k] * ab_x + g_z_0_xy_xxxxxxy[k];

                g_z_0_xxy_xxxxxz[k] = -g_z_0_xy_xxxxxz[k] * ab_x + g_z_0_xy_xxxxxxz[k];

                g_z_0_xxy_xxxxyy[k] = -g_z_0_xy_xxxxyy[k] * ab_x + g_z_0_xy_xxxxxyy[k];

                g_z_0_xxy_xxxxyz[k] = -g_z_0_xy_xxxxyz[k] * ab_x + g_z_0_xy_xxxxxyz[k];

                g_z_0_xxy_xxxxzz[k] = -g_z_0_xy_xxxxzz[k] * ab_x + g_z_0_xy_xxxxxzz[k];

                g_z_0_xxy_xxxyyy[k] = -g_z_0_xy_xxxyyy[k] * ab_x + g_z_0_xy_xxxxyyy[k];

                g_z_0_xxy_xxxyyz[k] = -g_z_0_xy_xxxyyz[k] * ab_x + g_z_0_xy_xxxxyyz[k];

                g_z_0_xxy_xxxyzz[k] = -g_z_0_xy_xxxyzz[k] * ab_x + g_z_0_xy_xxxxyzz[k];

                g_z_0_xxy_xxxzzz[k] = -g_z_0_xy_xxxzzz[k] * ab_x + g_z_0_xy_xxxxzzz[k];

                g_z_0_xxy_xxyyyy[k] = -g_z_0_xy_xxyyyy[k] * ab_x + g_z_0_xy_xxxyyyy[k];

                g_z_0_xxy_xxyyyz[k] = -g_z_0_xy_xxyyyz[k] * ab_x + g_z_0_xy_xxxyyyz[k];

                g_z_0_xxy_xxyyzz[k] = -g_z_0_xy_xxyyzz[k] * ab_x + g_z_0_xy_xxxyyzz[k];

                g_z_0_xxy_xxyzzz[k] = -g_z_0_xy_xxyzzz[k] * ab_x + g_z_0_xy_xxxyzzz[k];

                g_z_0_xxy_xxzzzz[k] = -g_z_0_xy_xxzzzz[k] * ab_x + g_z_0_xy_xxxzzzz[k];

                g_z_0_xxy_xyyyyy[k] = -g_z_0_xy_xyyyyy[k] * ab_x + g_z_0_xy_xxyyyyy[k];

                g_z_0_xxy_xyyyyz[k] = -g_z_0_xy_xyyyyz[k] * ab_x + g_z_0_xy_xxyyyyz[k];

                g_z_0_xxy_xyyyzz[k] = -g_z_0_xy_xyyyzz[k] * ab_x + g_z_0_xy_xxyyyzz[k];

                g_z_0_xxy_xyyzzz[k] = -g_z_0_xy_xyyzzz[k] * ab_x + g_z_0_xy_xxyyzzz[k];

                g_z_0_xxy_xyzzzz[k] = -g_z_0_xy_xyzzzz[k] * ab_x + g_z_0_xy_xxyzzzz[k];

                g_z_0_xxy_xzzzzz[k] = -g_z_0_xy_xzzzzz[k] * ab_x + g_z_0_xy_xxzzzzz[k];

                g_z_0_xxy_yyyyyy[k] = -g_z_0_xy_yyyyyy[k] * ab_x + g_z_0_xy_xyyyyyy[k];

                g_z_0_xxy_yyyyyz[k] = -g_z_0_xy_yyyyyz[k] * ab_x + g_z_0_xy_xyyyyyz[k];

                g_z_0_xxy_yyyyzz[k] = -g_z_0_xy_yyyyzz[k] * ab_x + g_z_0_xy_xyyyyzz[k];

                g_z_0_xxy_yyyzzz[k] = -g_z_0_xy_yyyzzz[k] * ab_x + g_z_0_xy_xyyyzzz[k];

                g_z_0_xxy_yyzzzz[k] = -g_z_0_xy_yyzzzz[k] * ab_x + g_z_0_xy_xyyzzzz[k];

                g_z_0_xxy_yzzzzz[k] = -g_z_0_xy_yzzzzz[k] * ab_x + g_z_0_xy_xyzzzzz[k];

                g_z_0_xxy_zzzzzz[k] = -g_z_0_xy_zzzzzz[k] * ab_x + g_z_0_xy_xzzzzzz[k];
            }

            /// Set up 616-644 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxz_xxxxxx = cbuffer.data(fi_geom_10_off + 616 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxy = cbuffer.data(fi_geom_10_off + 617 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxxz = cbuffer.data(fi_geom_10_off + 618 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxyy = cbuffer.data(fi_geom_10_off + 619 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxyz = cbuffer.data(fi_geom_10_off + 620 * ccomps * dcomps);

            auto g_z_0_xxz_xxxxzz = cbuffer.data(fi_geom_10_off + 621 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyyy = cbuffer.data(fi_geom_10_off + 622 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyyz = cbuffer.data(fi_geom_10_off + 623 * ccomps * dcomps);

            auto g_z_0_xxz_xxxyzz = cbuffer.data(fi_geom_10_off + 624 * ccomps * dcomps);

            auto g_z_0_xxz_xxxzzz = cbuffer.data(fi_geom_10_off + 625 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyyy = cbuffer.data(fi_geom_10_off + 626 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyyz = cbuffer.data(fi_geom_10_off + 627 * ccomps * dcomps);

            auto g_z_0_xxz_xxyyzz = cbuffer.data(fi_geom_10_off + 628 * ccomps * dcomps);

            auto g_z_0_xxz_xxyzzz = cbuffer.data(fi_geom_10_off + 629 * ccomps * dcomps);

            auto g_z_0_xxz_xxzzzz = cbuffer.data(fi_geom_10_off + 630 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyyy = cbuffer.data(fi_geom_10_off + 631 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyyz = cbuffer.data(fi_geom_10_off + 632 * ccomps * dcomps);

            auto g_z_0_xxz_xyyyzz = cbuffer.data(fi_geom_10_off + 633 * ccomps * dcomps);

            auto g_z_0_xxz_xyyzzz = cbuffer.data(fi_geom_10_off + 634 * ccomps * dcomps);

            auto g_z_0_xxz_xyzzzz = cbuffer.data(fi_geom_10_off + 635 * ccomps * dcomps);

            auto g_z_0_xxz_xzzzzz = cbuffer.data(fi_geom_10_off + 636 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyyy = cbuffer.data(fi_geom_10_off + 637 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyyz = cbuffer.data(fi_geom_10_off + 638 * ccomps * dcomps);

            auto g_z_0_xxz_yyyyzz = cbuffer.data(fi_geom_10_off + 639 * ccomps * dcomps);

            auto g_z_0_xxz_yyyzzz = cbuffer.data(fi_geom_10_off + 640 * ccomps * dcomps);

            auto g_z_0_xxz_yyzzzz = cbuffer.data(fi_geom_10_off + 641 * ccomps * dcomps);

            auto g_z_0_xxz_yzzzzz = cbuffer.data(fi_geom_10_off + 642 * ccomps * dcomps);

            auto g_z_0_xxz_zzzzzz = cbuffer.data(fi_geom_10_off + 643 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxz_xxxxxx, g_z_0_xxz_xxxxxy, g_z_0_xxz_xxxxxz, g_z_0_xxz_xxxxyy, g_z_0_xxz_xxxxyz, g_z_0_xxz_xxxxzz, g_z_0_xxz_xxxyyy, g_z_0_xxz_xxxyyz, g_z_0_xxz_xxxyzz, g_z_0_xxz_xxxzzz, g_z_0_xxz_xxyyyy, g_z_0_xxz_xxyyyz, g_z_0_xxz_xxyyzz, g_z_0_xxz_xxyzzz, g_z_0_xxz_xxzzzz, g_z_0_xxz_xyyyyy, g_z_0_xxz_xyyyyz, g_z_0_xxz_xyyyzz, g_z_0_xxz_xyyzzz, g_z_0_xxz_xyzzzz, g_z_0_xxz_xzzzzz, g_z_0_xxz_yyyyyy, g_z_0_xxz_yyyyyz, g_z_0_xxz_yyyyzz, g_z_0_xxz_yyyzzz, g_z_0_xxz_yyzzzz, g_z_0_xxz_yzzzzz, g_z_0_xxz_zzzzzz, g_z_0_xz_xxxxxx, g_z_0_xz_xxxxxxx, g_z_0_xz_xxxxxxy, g_z_0_xz_xxxxxxz, g_z_0_xz_xxxxxy, g_z_0_xz_xxxxxyy, g_z_0_xz_xxxxxyz, g_z_0_xz_xxxxxz, g_z_0_xz_xxxxxzz, g_z_0_xz_xxxxyy, g_z_0_xz_xxxxyyy, g_z_0_xz_xxxxyyz, g_z_0_xz_xxxxyz, g_z_0_xz_xxxxyzz, g_z_0_xz_xxxxzz, g_z_0_xz_xxxxzzz, g_z_0_xz_xxxyyy, g_z_0_xz_xxxyyyy, g_z_0_xz_xxxyyyz, g_z_0_xz_xxxyyz, g_z_0_xz_xxxyyzz, g_z_0_xz_xxxyzz, g_z_0_xz_xxxyzzz, g_z_0_xz_xxxzzz, g_z_0_xz_xxxzzzz, g_z_0_xz_xxyyyy, g_z_0_xz_xxyyyyy, g_z_0_xz_xxyyyyz, g_z_0_xz_xxyyyz, g_z_0_xz_xxyyyzz, g_z_0_xz_xxyyzz, g_z_0_xz_xxyyzzz, g_z_0_xz_xxyzzz, g_z_0_xz_xxyzzzz, g_z_0_xz_xxzzzz, g_z_0_xz_xxzzzzz, g_z_0_xz_xyyyyy, g_z_0_xz_xyyyyyy, g_z_0_xz_xyyyyyz, g_z_0_xz_xyyyyz, g_z_0_xz_xyyyyzz, g_z_0_xz_xyyyzz, g_z_0_xz_xyyyzzz, g_z_0_xz_xyyzzz, g_z_0_xz_xyyzzzz, g_z_0_xz_xyzzzz, g_z_0_xz_xyzzzzz, g_z_0_xz_xzzzzz, g_z_0_xz_xzzzzzz, g_z_0_xz_yyyyyy, g_z_0_xz_yyyyyz, g_z_0_xz_yyyyzz, g_z_0_xz_yyyzzz, g_z_0_xz_yyzzzz, g_z_0_xz_yzzzzz, g_z_0_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxz_xxxxxx[k] = -g_z_0_xz_xxxxxx[k] * ab_x + g_z_0_xz_xxxxxxx[k];

                g_z_0_xxz_xxxxxy[k] = -g_z_0_xz_xxxxxy[k] * ab_x + g_z_0_xz_xxxxxxy[k];

                g_z_0_xxz_xxxxxz[k] = -g_z_0_xz_xxxxxz[k] * ab_x + g_z_0_xz_xxxxxxz[k];

                g_z_0_xxz_xxxxyy[k] = -g_z_0_xz_xxxxyy[k] * ab_x + g_z_0_xz_xxxxxyy[k];

                g_z_0_xxz_xxxxyz[k] = -g_z_0_xz_xxxxyz[k] * ab_x + g_z_0_xz_xxxxxyz[k];

                g_z_0_xxz_xxxxzz[k] = -g_z_0_xz_xxxxzz[k] * ab_x + g_z_0_xz_xxxxxzz[k];

                g_z_0_xxz_xxxyyy[k] = -g_z_0_xz_xxxyyy[k] * ab_x + g_z_0_xz_xxxxyyy[k];

                g_z_0_xxz_xxxyyz[k] = -g_z_0_xz_xxxyyz[k] * ab_x + g_z_0_xz_xxxxyyz[k];

                g_z_0_xxz_xxxyzz[k] = -g_z_0_xz_xxxyzz[k] * ab_x + g_z_0_xz_xxxxyzz[k];

                g_z_0_xxz_xxxzzz[k] = -g_z_0_xz_xxxzzz[k] * ab_x + g_z_0_xz_xxxxzzz[k];

                g_z_0_xxz_xxyyyy[k] = -g_z_0_xz_xxyyyy[k] * ab_x + g_z_0_xz_xxxyyyy[k];

                g_z_0_xxz_xxyyyz[k] = -g_z_0_xz_xxyyyz[k] * ab_x + g_z_0_xz_xxxyyyz[k];

                g_z_0_xxz_xxyyzz[k] = -g_z_0_xz_xxyyzz[k] * ab_x + g_z_0_xz_xxxyyzz[k];

                g_z_0_xxz_xxyzzz[k] = -g_z_0_xz_xxyzzz[k] * ab_x + g_z_0_xz_xxxyzzz[k];

                g_z_0_xxz_xxzzzz[k] = -g_z_0_xz_xxzzzz[k] * ab_x + g_z_0_xz_xxxzzzz[k];

                g_z_0_xxz_xyyyyy[k] = -g_z_0_xz_xyyyyy[k] * ab_x + g_z_0_xz_xxyyyyy[k];

                g_z_0_xxz_xyyyyz[k] = -g_z_0_xz_xyyyyz[k] * ab_x + g_z_0_xz_xxyyyyz[k];

                g_z_0_xxz_xyyyzz[k] = -g_z_0_xz_xyyyzz[k] * ab_x + g_z_0_xz_xxyyyzz[k];

                g_z_0_xxz_xyyzzz[k] = -g_z_0_xz_xyyzzz[k] * ab_x + g_z_0_xz_xxyyzzz[k];

                g_z_0_xxz_xyzzzz[k] = -g_z_0_xz_xyzzzz[k] * ab_x + g_z_0_xz_xxyzzzz[k];

                g_z_0_xxz_xzzzzz[k] = -g_z_0_xz_xzzzzz[k] * ab_x + g_z_0_xz_xxzzzzz[k];

                g_z_0_xxz_yyyyyy[k] = -g_z_0_xz_yyyyyy[k] * ab_x + g_z_0_xz_xyyyyyy[k];

                g_z_0_xxz_yyyyyz[k] = -g_z_0_xz_yyyyyz[k] * ab_x + g_z_0_xz_xyyyyyz[k];

                g_z_0_xxz_yyyyzz[k] = -g_z_0_xz_yyyyzz[k] * ab_x + g_z_0_xz_xyyyyzz[k];

                g_z_0_xxz_yyyzzz[k] = -g_z_0_xz_yyyzzz[k] * ab_x + g_z_0_xz_xyyyzzz[k];

                g_z_0_xxz_yyzzzz[k] = -g_z_0_xz_yyzzzz[k] * ab_x + g_z_0_xz_xyyzzzz[k];

                g_z_0_xxz_yzzzzz[k] = -g_z_0_xz_yzzzzz[k] * ab_x + g_z_0_xz_xyzzzzz[k];

                g_z_0_xxz_zzzzzz[k] = -g_z_0_xz_zzzzzz[k] * ab_x + g_z_0_xz_xzzzzzz[k];
            }

            /// Set up 644-672 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyy_xxxxxx = cbuffer.data(fi_geom_10_off + 644 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxy = cbuffer.data(fi_geom_10_off + 645 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxxz = cbuffer.data(fi_geom_10_off + 646 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxyy = cbuffer.data(fi_geom_10_off + 647 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxyz = cbuffer.data(fi_geom_10_off + 648 * ccomps * dcomps);

            auto g_z_0_xyy_xxxxzz = cbuffer.data(fi_geom_10_off + 649 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyyy = cbuffer.data(fi_geom_10_off + 650 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyyz = cbuffer.data(fi_geom_10_off + 651 * ccomps * dcomps);

            auto g_z_0_xyy_xxxyzz = cbuffer.data(fi_geom_10_off + 652 * ccomps * dcomps);

            auto g_z_0_xyy_xxxzzz = cbuffer.data(fi_geom_10_off + 653 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyyy = cbuffer.data(fi_geom_10_off + 654 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyyz = cbuffer.data(fi_geom_10_off + 655 * ccomps * dcomps);

            auto g_z_0_xyy_xxyyzz = cbuffer.data(fi_geom_10_off + 656 * ccomps * dcomps);

            auto g_z_0_xyy_xxyzzz = cbuffer.data(fi_geom_10_off + 657 * ccomps * dcomps);

            auto g_z_0_xyy_xxzzzz = cbuffer.data(fi_geom_10_off + 658 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyyy = cbuffer.data(fi_geom_10_off + 659 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyyz = cbuffer.data(fi_geom_10_off + 660 * ccomps * dcomps);

            auto g_z_0_xyy_xyyyzz = cbuffer.data(fi_geom_10_off + 661 * ccomps * dcomps);

            auto g_z_0_xyy_xyyzzz = cbuffer.data(fi_geom_10_off + 662 * ccomps * dcomps);

            auto g_z_0_xyy_xyzzzz = cbuffer.data(fi_geom_10_off + 663 * ccomps * dcomps);

            auto g_z_0_xyy_xzzzzz = cbuffer.data(fi_geom_10_off + 664 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyyy = cbuffer.data(fi_geom_10_off + 665 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyyz = cbuffer.data(fi_geom_10_off + 666 * ccomps * dcomps);

            auto g_z_0_xyy_yyyyzz = cbuffer.data(fi_geom_10_off + 667 * ccomps * dcomps);

            auto g_z_0_xyy_yyyzzz = cbuffer.data(fi_geom_10_off + 668 * ccomps * dcomps);

            auto g_z_0_xyy_yyzzzz = cbuffer.data(fi_geom_10_off + 669 * ccomps * dcomps);

            auto g_z_0_xyy_yzzzzz = cbuffer.data(fi_geom_10_off + 670 * ccomps * dcomps);

            auto g_z_0_xyy_zzzzzz = cbuffer.data(fi_geom_10_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyy_xxxxxx, g_z_0_xyy_xxxxxy, g_z_0_xyy_xxxxxz, g_z_0_xyy_xxxxyy, g_z_0_xyy_xxxxyz, g_z_0_xyy_xxxxzz, g_z_0_xyy_xxxyyy, g_z_0_xyy_xxxyyz, g_z_0_xyy_xxxyzz, g_z_0_xyy_xxxzzz, g_z_0_xyy_xxyyyy, g_z_0_xyy_xxyyyz, g_z_0_xyy_xxyyzz, g_z_0_xyy_xxyzzz, g_z_0_xyy_xxzzzz, g_z_0_xyy_xyyyyy, g_z_0_xyy_xyyyyz, g_z_0_xyy_xyyyzz, g_z_0_xyy_xyyzzz, g_z_0_xyy_xyzzzz, g_z_0_xyy_xzzzzz, g_z_0_xyy_yyyyyy, g_z_0_xyy_yyyyyz, g_z_0_xyy_yyyyzz, g_z_0_xyy_yyyzzz, g_z_0_xyy_yyzzzz, g_z_0_xyy_yzzzzz, g_z_0_xyy_zzzzzz, g_z_0_yy_xxxxxx, g_z_0_yy_xxxxxxx, g_z_0_yy_xxxxxxy, g_z_0_yy_xxxxxxz, g_z_0_yy_xxxxxy, g_z_0_yy_xxxxxyy, g_z_0_yy_xxxxxyz, g_z_0_yy_xxxxxz, g_z_0_yy_xxxxxzz, g_z_0_yy_xxxxyy, g_z_0_yy_xxxxyyy, g_z_0_yy_xxxxyyz, g_z_0_yy_xxxxyz, g_z_0_yy_xxxxyzz, g_z_0_yy_xxxxzz, g_z_0_yy_xxxxzzz, g_z_0_yy_xxxyyy, g_z_0_yy_xxxyyyy, g_z_0_yy_xxxyyyz, g_z_0_yy_xxxyyz, g_z_0_yy_xxxyyzz, g_z_0_yy_xxxyzz, g_z_0_yy_xxxyzzz, g_z_0_yy_xxxzzz, g_z_0_yy_xxxzzzz, g_z_0_yy_xxyyyy, g_z_0_yy_xxyyyyy, g_z_0_yy_xxyyyyz, g_z_0_yy_xxyyyz, g_z_0_yy_xxyyyzz, g_z_0_yy_xxyyzz, g_z_0_yy_xxyyzzz, g_z_0_yy_xxyzzz, g_z_0_yy_xxyzzzz, g_z_0_yy_xxzzzz, g_z_0_yy_xxzzzzz, g_z_0_yy_xyyyyy, g_z_0_yy_xyyyyyy, g_z_0_yy_xyyyyyz, g_z_0_yy_xyyyyz, g_z_0_yy_xyyyyzz, g_z_0_yy_xyyyzz, g_z_0_yy_xyyyzzz, g_z_0_yy_xyyzzz, g_z_0_yy_xyyzzzz, g_z_0_yy_xyzzzz, g_z_0_yy_xyzzzzz, g_z_0_yy_xzzzzz, g_z_0_yy_xzzzzzz, g_z_0_yy_yyyyyy, g_z_0_yy_yyyyyz, g_z_0_yy_yyyyzz, g_z_0_yy_yyyzzz, g_z_0_yy_yyzzzz, g_z_0_yy_yzzzzz, g_z_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyy_xxxxxx[k] = -g_z_0_yy_xxxxxx[k] * ab_x + g_z_0_yy_xxxxxxx[k];

                g_z_0_xyy_xxxxxy[k] = -g_z_0_yy_xxxxxy[k] * ab_x + g_z_0_yy_xxxxxxy[k];

                g_z_0_xyy_xxxxxz[k] = -g_z_0_yy_xxxxxz[k] * ab_x + g_z_0_yy_xxxxxxz[k];

                g_z_0_xyy_xxxxyy[k] = -g_z_0_yy_xxxxyy[k] * ab_x + g_z_0_yy_xxxxxyy[k];

                g_z_0_xyy_xxxxyz[k] = -g_z_0_yy_xxxxyz[k] * ab_x + g_z_0_yy_xxxxxyz[k];

                g_z_0_xyy_xxxxzz[k] = -g_z_0_yy_xxxxzz[k] * ab_x + g_z_0_yy_xxxxxzz[k];

                g_z_0_xyy_xxxyyy[k] = -g_z_0_yy_xxxyyy[k] * ab_x + g_z_0_yy_xxxxyyy[k];

                g_z_0_xyy_xxxyyz[k] = -g_z_0_yy_xxxyyz[k] * ab_x + g_z_0_yy_xxxxyyz[k];

                g_z_0_xyy_xxxyzz[k] = -g_z_0_yy_xxxyzz[k] * ab_x + g_z_0_yy_xxxxyzz[k];

                g_z_0_xyy_xxxzzz[k] = -g_z_0_yy_xxxzzz[k] * ab_x + g_z_0_yy_xxxxzzz[k];

                g_z_0_xyy_xxyyyy[k] = -g_z_0_yy_xxyyyy[k] * ab_x + g_z_0_yy_xxxyyyy[k];

                g_z_0_xyy_xxyyyz[k] = -g_z_0_yy_xxyyyz[k] * ab_x + g_z_0_yy_xxxyyyz[k];

                g_z_0_xyy_xxyyzz[k] = -g_z_0_yy_xxyyzz[k] * ab_x + g_z_0_yy_xxxyyzz[k];

                g_z_0_xyy_xxyzzz[k] = -g_z_0_yy_xxyzzz[k] * ab_x + g_z_0_yy_xxxyzzz[k];

                g_z_0_xyy_xxzzzz[k] = -g_z_0_yy_xxzzzz[k] * ab_x + g_z_0_yy_xxxzzzz[k];

                g_z_0_xyy_xyyyyy[k] = -g_z_0_yy_xyyyyy[k] * ab_x + g_z_0_yy_xxyyyyy[k];

                g_z_0_xyy_xyyyyz[k] = -g_z_0_yy_xyyyyz[k] * ab_x + g_z_0_yy_xxyyyyz[k];

                g_z_0_xyy_xyyyzz[k] = -g_z_0_yy_xyyyzz[k] * ab_x + g_z_0_yy_xxyyyzz[k];

                g_z_0_xyy_xyyzzz[k] = -g_z_0_yy_xyyzzz[k] * ab_x + g_z_0_yy_xxyyzzz[k];

                g_z_0_xyy_xyzzzz[k] = -g_z_0_yy_xyzzzz[k] * ab_x + g_z_0_yy_xxyzzzz[k];

                g_z_0_xyy_xzzzzz[k] = -g_z_0_yy_xzzzzz[k] * ab_x + g_z_0_yy_xxzzzzz[k];

                g_z_0_xyy_yyyyyy[k] = -g_z_0_yy_yyyyyy[k] * ab_x + g_z_0_yy_xyyyyyy[k];

                g_z_0_xyy_yyyyyz[k] = -g_z_0_yy_yyyyyz[k] * ab_x + g_z_0_yy_xyyyyyz[k];

                g_z_0_xyy_yyyyzz[k] = -g_z_0_yy_yyyyzz[k] * ab_x + g_z_0_yy_xyyyyzz[k];

                g_z_0_xyy_yyyzzz[k] = -g_z_0_yy_yyyzzz[k] * ab_x + g_z_0_yy_xyyyzzz[k];

                g_z_0_xyy_yyzzzz[k] = -g_z_0_yy_yyzzzz[k] * ab_x + g_z_0_yy_xyyzzzz[k];

                g_z_0_xyy_yzzzzz[k] = -g_z_0_yy_yzzzzz[k] * ab_x + g_z_0_yy_xyzzzzz[k];

                g_z_0_xyy_zzzzzz[k] = -g_z_0_yy_zzzzzz[k] * ab_x + g_z_0_yy_xzzzzzz[k];
            }

            /// Set up 672-700 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyz_xxxxxx = cbuffer.data(fi_geom_10_off + 672 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxy = cbuffer.data(fi_geom_10_off + 673 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxxz = cbuffer.data(fi_geom_10_off + 674 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxyy = cbuffer.data(fi_geom_10_off + 675 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxyz = cbuffer.data(fi_geom_10_off + 676 * ccomps * dcomps);

            auto g_z_0_xyz_xxxxzz = cbuffer.data(fi_geom_10_off + 677 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyyy = cbuffer.data(fi_geom_10_off + 678 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyyz = cbuffer.data(fi_geom_10_off + 679 * ccomps * dcomps);

            auto g_z_0_xyz_xxxyzz = cbuffer.data(fi_geom_10_off + 680 * ccomps * dcomps);

            auto g_z_0_xyz_xxxzzz = cbuffer.data(fi_geom_10_off + 681 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyyy = cbuffer.data(fi_geom_10_off + 682 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyyz = cbuffer.data(fi_geom_10_off + 683 * ccomps * dcomps);

            auto g_z_0_xyz_xxyyzz = cbuffer.data(fi_geom_10_off + 684 * ccomps * dcomps);

            auto g_z_0_xyz_xxyzzz = cbuffer.data(fi_geom_10_off + 685 * ccomps * dcomps);

            auto g_z_0_xyz_xxzzzz = cbuffer.data(fi_geom_10_off + 686 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyyy = cbuffer.data(fi_geom_10_off + 687 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyyz = cbuffer.data(fi_geom_10_off + 688 * ccomps * dcomps);

            auto g_z_0_xyz_xyyyzz = cbuffer.data(fi_geom_10_off + 689 * ccomps * dcomps);

            auto g_z_0_xyz_xyyzzz = cbuffer.data(fi_geom_10_off + 690 * ccomps * dcomps);

            auto g_z_0_xyz_xyzzzz = cbuffer.data(fi_geom_10_off + 691 * ccomps * dcomps);

            auto g_z_0_xyz_xzzzzz = cbuffer.data(fi_geom_10_off + 692 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyyy = cbuffer.data(fi_geom_10_off + 693 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyyz = cbuffer.data(fi_geom_10_off + 694 * ccomps * dcomps);

            auto g_z_0_xyz_yyyyzz = cbuffer.data(fi_geom_10_off + 695 * ccomps * dcomps);

            auto g_z_0_xyz_yyyzzz = cbuffer.data(fi_geom_10_off + 696 * ccomps * dcomps);

            auto g_z_0_xyz_yyzzzz = cbuffer.data(fi_geom_10_off + 697 * ccomps * dcomps);

            auto g_z_0_xyz_yzzzzz = cbuffer.data(fi_geom_10_off + 698 * ccomps * dcomps);

            auto g_z_0_xyz_zzzzzz = cbuffer.data(fi_geom_10_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyz_xxxxxx, g_z_0_xyz_xxxxxy, g_z_0_xyz_xxxxxz, g_z_0_xyz_xxxxyy, g_z_0_xyz_xxxxyz, g_z_0_xyz_xxxxzz, g_z_0_xyz_xxxyyy, g_z_0_xyz_xxxyyz, g_z_0_xyz_xxxyzz, g_z_0_xyz_xxxzzz, g_z_0_xyz_xxyyyy, g_z_0_xyz_xxyyyz, g_z_0_xyz_xxyyzz, g_z_0_xyz_xxyzzz, g_z_0_xyz_xxzzzz, g_z_0_xyz_xyyyyy, g_z_0_xyz_xyyyyz, g_z_0_xyz_xyyyzz, g_z_0_xyz_xyyzzz, g_z_0_xyz_xyzzzz, g_z_0_xyz_xzzzzz, g_z_0_xyz_yyyyyy, g_z_0_xyz_yyyyyz, g_z_0_xyz_yyyyzz, g_z_0_xyz_yyyzzz, g_z_0_xyz_yyzzzz, g_z_0_xyz_yzzzzz, g_z_0_xyz_zzzzzz, g_z_0_yz_xxxxxx, g_z_0_yz_xxxxxxx, g_z_0_yz_xxxxxxy, g_z_0_yz_xxxxxxz, g_z_0_yz_xxxxxy, g_z_0_yz_xxxxxyy, g_z_0_yz_xxxxxyz, g_z_0_yz_xxxxxz, g_z_0_yz_xxxxxzz, g_z_0_yz_xxxxyy, g_z_0_yz_xxxxyyy, g_z_0_yz_xxxxyyz, g_z_0_yz_xxxxyz, g_z_0_yz_xxxxyzz, g_z_0_yz_xxxxzz, g_z_0_yz_xxxxzzz, g_z_0_yz_xxxyyy, g_z_0_yz_xxxyyyy, g_z_0_yz_xxxyyyz, g_z_0_yz_xxxyyz, g_z_0_yz_xxxyyzz, g_z_0_yz_xxxyzz, g_z_0_yz_xxxyzzz, g_z_0_yz_xxxzzz, g_z_0_yz_xxxzzzz, g_z_0_yz_xxyyyy, g_z_0_yz_xxyyyyy, g_z_0_yz_xxyyyyz, g_z_0_yz_xxyyyz, g_z_0_yz_xxyyyzz, g_z_0_yz_xxyyzz, g_z_0_yz_xxyyzzz, g_z_0_yz_xxyzzz, g_z_0_yz_xxyzzzz, g_z_0_yz_xxzzzz, g_z_0_yz_xxzzzzz, g_z_0_yz_xyyyyy, g_z_0_yz_xyyyyyy, g_z_0_yz_xyyyyyz, g_z_0_yz_xyyyyz, g_z_0_yz_xyyyyzz, g_z_0_yz_xyyyzz, g_z_0_yz_xyyyzzz, g_z_0_yz_xyyzzz, g_z_0_yz_xyyzzzz, g_z_0_yz_xyzzzz, g_z_0_yz_xyzzzzz, g_z_0_yz_xzzzzz, g_z_0_yz_xzzzzzz, g_z_0_yz_yyyyyy, g_z_0_yz_yyyyyz, g_z_0_yz_yyyyzz, g_z_0_yz_yyyzzz, g_z_0_yz_yyzzzz, g_z_0_yz_yzzzzz, g_z_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyz_xxxxxx[k] = -g_z_0_yz_xxxxxx[k] * ab_x + g_z_0_yz_xxxxxxx[k];

                g_z_0_xyz_xxxxxy[k] = -g_z_0_yz_xxxxxy[k] * ab_x + g_z_0_yz_xxxxxxy[k];

                g_z_0_xyz_xxxxxz[k] = -g_z_0_yz_xxxxxz[k] * ab_x + g_z_0_yz_xxxxxxz[k];

                g_z_0_xyz_xxxxyy[k] = -g_z_0_yz_xxxxyy[k] * ab_x + g_z_0_yz_xxxxxyy[k];

                g_z_0_xyz_xxxxyz[k] = -g_z_0_yz_xxxxyz[k] * ab_x + g_z_0_yz_xxxxxyz[k];

                g_z_0_xyz_xxxxzz[k] = -g_z_0_yz_xxxxzz[k] * ab_x + g_z_0_yz_xxxxxzz[k];

                g_z_0_xyz_xxxyyy[k] = -g_z_0_yz_xxxyyy[k] * ab_x + g_z_0_yz_xxxxyyy[k];

                g_z_0_xyz_xxxyyz[k] = -g_z_0_yz_xxxyyz[k] * ab_x + g_z_0_yz_xxxxyyz[k];

                g_z_0_xyz_xxxyzz[k] = -g_z_0_yz_xxxyzz[k] * ab_x + g_z_0_yz_xxxxyzz[k];

                g_z_0_xyz_xxxzzz[k] = -g_z_0_yz_xxxzzz[k] * ab_x + g_z_0_yz_xxxxzzz[k];

                g_z_0_xyz_xxyyyy[k] = -g_z_0_yz_xxyyyy[k] * ab_x + g_z_0_yz_xxxyyyy[k];

                g_z_0_xyz_xxyyyz[k] = -g_z_0_yz_xxyyyz[k] * ab_x + g_z_0_yz_xxxyyyz[k];

                g_z_0_xyz_xxyyzz[k] = -g_z_0_yz_xxyyzz[k] * ab_x + g_z_0_yz_xxxyyzz[k];

                g_z_0_xyz_xxyzzz[k] = -g_z_0_yz_xxyzzz[k] * ab_x + g_z_0_yz_xxxyzzz[k];

                g_z_0_xyz_xxzzzz[k] = -g_z_0_yz_xxzzzz[k] * ab_x + g_z_0_yz_xxxzzzz[k];

                g_z_0_xyz_xyyyyy[k] = -g_z_0_yz_xyyyyy[k] * ab_x + g_z_0_yz_xxyyyyy[k];

                g_z_0_xyz_xyyyyz[k] = -g_z_0_yz_xyyyyz[k] * ab_x + g_z_0_yz_xxyyyyz[k];

                g_z_0_xyz_xyyyzz[k] = -g_z_0_yz_xyyyzz[k] * ab_x + g_z_0_yz_xxyyyzz[k];

                g_z_0_xyz_xyyzzz[k] = -g_z_0_yz_xyyzzz[k] * ab_x + g_z_0_yz_xxyyzzz[k];

                g_z_0_xyz_xyzzzz[k] = -g_z_0_yz_xyzzzz[k] * ab_x + g_z_0_yz_xxyzzzz[k];

                g_z_0_xyz_xzzzzz[k] = -g_z_0_yz_xzzzzz[k] * ab_x + g_z_0_yz_xxzzzzz[k];

                g_z_0_xyz_yyyyyy[k] = -g_z_0_yz_yyyyyy[k] * ab_x + g_z_0_yz_xyyyyyy[k];

                g_z_0_xyz_yyyyyz[k] = -g_z_0_yz_yyyyyz[k] * ab_x + g_z_0_yz_xyyyyyz[k];

                g_z_0_xyz_yyyyzz[k] = -g_z_0_yz_yyyyzz[k] * ab_x + g_z_0_yz_xyyyyzz[k];

                g_z_0_xyz_yyyzzz[k] = -g_z_0_yz_yyyzzz[k] * ab_x + g_z_0_yz_xyyyzzz[k];

                g_z_0_xyz_yyzzzz[k] = -g_z_0_yz_yyzzzz[k] * ab_x + g_z_0_yz_xyyzzzz[k];

                g_z_0_xyz_yzzzzz[k] = -g_z_0_yz_yzzzzz[k] * ab_x + g_z_0_yz_xyzzzzz[k];

                g_z_0_xyz_zzzzzz[k] = -g_z_0_yz_zzzzzz[k] * ab_x + g_z_0_yz_xzzzzzz[k];
            }

            /// Set up 700-728 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzz_xxxxxx = cbuffer.data(fi_geom_10_off + 700 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxy = cbuffer.data(fi_geom_10_off + 701 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxxz = cbuffer.data(fi_geom_10_off + 702 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxyy = cbuffer.data(fi_geom_10_off + 703 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxyz = cbuffer.data(fi_geom_10_off + 704 * ccomps * dcomps);

            auto g_z_0_xzz_xxxxzz = cbuffer.data(fi_geom_10_off + 705 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyyy = cbuffer.data(fi_geom_10_off + 706 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyyz = cbuffer.data(fi_geom_10_off + 707 * ccomps * dcomps);

            auto g_z_0_xzz_xxxyzz = cbuffer.data(fi_geom_10_off + 708 * ccomps * dcomps);

            auto g_z_0_xzz_xxxzzz = cbuffer.data(fi_geom_10_off + 709 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyyy = cbuffer.data(fi_geom_10_off + 710 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyyz = cbuffer.data(fi_geom_10_off + 711 * ccomps * dcomps);

            auto g_z_0_xzz_xxyyzz = cbuffer.data(fi_geom_10_off + 712 * ccomps * dcomps);

            auto g_z_0_xzz_xxyzzz = cbuffer.data(fi_geom_10_off + 713 * ccomps * dcomps);

            auto g_z_0_xzz_xxzzzz = cbuffer.data(fi_geom_10_off + 714 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyyy = cbuffer.data(fi_geom_10_off + 715 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyyz = cbuffer.data(fi_geom_10_off + 716 * ccomps * dcomps);

            auto g_z_0_xzz_xyyyzz = cbuffer.data(fi_geom_10_off + 717 * ccomps * dcomps);

            auto g_z_0_xzz_xyyzzz = cbuffer.data(fi_geom_10_off + 718 * ccomps * dcomps);

            auto g_z_0_xzz_xyzzzz = cbuffer.data(fi_geom_10_off + 719 * ccomps * dcomps);

            auto g_z_0_xzz_xzzzzz = cbuffer.data(fi_geom_10_off + 720 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyyy = cbuffer.data(fi_geom_10_off + 721 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyyz = cbuffer.data(fi_geom_10_off + 722 * ccomps * dcomps);

            auto g_z_0_xzz_yyyyzz = cbuffer.data(fi_geom_10_off + 723 * ccomps * dcomps);

            auto g_z_0_xzz_yyyzzz = cbuffer.data(fi_geom_10_off + 724 * ccomps * dcomps);

            auto g_z_0_xzz_yyzzzz = cbuffer.data(fi_geom_10_off + 725 * ccomps * dcomps);

            auto g_z_0_xzz_yzzzzz = cbuffer.data(fi_geom_10_off + 726 * ccomps * dcomps);

            auto g_z_0_xzz_zzzzzz = cbuffer.data(fi_geom_10_off + 727 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzz_xxxxxx, g_z_0_xzz_xxxxxy, g_z_0_xzz_xxxxxz, g_z_0_xzz_xxxxyy, g_z_0_xzz_xxxxyz, g_z_0_xzz_xxxxzz, g_z_0_xzz_xxxyyy, g_z_0_xzz_xxxyyz, g_z_0_xzz_xxxyzz, g_z_0_xzz_xxxzzz, g_z_0_xzz_xxyyyy, g_z_0_xzz_xxyyyz, g_z_0_xzz_xxyyzz, g_z_0_xzz_xxyzzz, g_z_0_xzz_xxzzzz, g_z_0_xzz_xyyyyy, g_z_0_xzz_xyyyyz, g_z_0_xzz_xyyyzz, g_z_0_xzz_xyyzzz, g_z_0_xzz_xyzzzz, g_z_0_xzz_xzzzzz, g_z_0_xzz_yyyyyy, g_z_0_xzz_yyyyyz, g_z_0_xzz_yyyyzz, g_z_0_xzz_yyyzzz, g_z_0_xzz_yyzzzz, g_z_0_xzz_yzzzzz, g_z_0_xzz_zzzzzz, g_z_0_zz_xxxxxx, g_z_0_zz_xxxxxxx, g_z_0_zz_xxxxxxy, g_z_0_zz_xxxxxxz, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxxyy, g_z_0_zz_xxxxxyz, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxxzz, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyyy, g_z_0_zz_xxxxyyz, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxyzz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxxzzz, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyyy, g_z_0_zz_xxxyyyz, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyyzz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxyzzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxxzzzz, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyyy, g_z_0_zz_xxyyyyz, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyyzz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyyzzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxyzzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xxzzzzz, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyyy, g_z_0_zz_xyyyyyz, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyyzz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyyzzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyyzzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xyzzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_xzzzzzz, g_z_0_zz_yyyyyy, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzz_xxxxxx[k] = -g_z_0_zz_xxxxxx[k] * ab_x + g_z_0_zz_xxxxxxx[k];

                g_z_0_xzz_xxxxxy[k] = -g_z_0_zz_xxxxxy[k] * ab_x + g_z_0_zz_xxxxxxy[k];

                g_z_0_xzz_xxxxxz[k] = -g_z_0_zz_xxxxxz[k] * ab_x + g_z_0_zz_xxxxxxz[k];

                g_z_0_xzz_xxxxyy[k] = -g_z_0_zz_xxxxyy[k] * ab_x + g_z_0_zz_xxxxxyy[k];

                g_z_0_xzz_xxxxyz[k] = -g_z_0_zz_xxxxyz[k] * ab_x + g_z_0_zz_xxxxxyz[k];

                g_z_0_xzz_xxxxzz[k] = -g_z_0_zz_xxxxzz[k] * ab_x + g_z_0_zz_xxxxxzz[k];

                g_z_0_xzz_xxxyyy[k] = -g_z_0_zz_xxxyyy[k] * ab_x + g_z_0_zz_xxxxyyy[k];

                g_z_0_xzz_xxxyyz[k] = -g_z_0_zz_xxxyyz[k] * ab_x + g_z_0_zz_xxxxyyz[k];

                g_z_0_xzz_xxxyzz[k] = -g_z_0_zz_xxxyzz[k] * ab_x + g_z_0_zz_xxxxyzz[k];

                g_z_0_xzz_xxxzzz[k] = -g_z_0_zz_xxxzzz[k] * ab_x + g_z_0_zz_xxxxzzz[k];

                g_z_0_xzz_xxyyyy[k] = -g_z_0_zz_xxyyyy[k] * ab_x + g_z_0_zz_xxxyyyy[k];

                g_z_0_xzz_xxyyyz[k] = -g_z_0_zz_xxyyyz[k] * ab_x + g_z_0_zz_xxxyyyz[k];

                g_z_0_xzz_xxyyzz[k] = -g_z_0_zz_xxyyzz[k] * ab_x + g_z_0_zz_xxxyyzz[k];

                g_z_0_xzz_xxyzzz[k] = -g_z_0_zz_xxyzzz[k] * ab_x + g_z_0_zz_xxxyzzz[k];

                g_z_0_xzz_xxzzzz[k] = -g_z_0_zz_xxzzzz[k] * ab_x + g_z_0_zz_xxxzzzz[k];

                g_z_0_xzz_xyyyyy[k] = -g_z_0_zz_xyyyyy[k] * ab_x + g_z_0_zz_xxyyyyy[k];

                g_z_0_xzz_xyyyyz[k] = -g_z_0_zz_xyyyyz[k] * ab_x + g_z_0_zz_xxyyyyz[k];

                g_z_0_xzz_xyyyzz[k] = -g_z_0_zz_xyyyzz[k] * ab_x + g_z_0_zz_xxyyyzz[k];

                g_z_0_xzz_xyyzzz[k] = -g_z_0_zz_xyyzzz[k] * ab_x + g_z_0_zz_xxyyzzz[k];

                g_z_0_xzz_xyzzzz[k] = -g_z_0_zz_xyzzzz[k] * ab_x + g_z_0_zz_xxyzzzz[k];

                g_z_0_xzz_xzzzzz[k] = -g_z_0_zz_xzzzzz[k] * ab_x + g_z_0_zz_xxzzzzz[k];

                g_z_0_xzz_yyyyyy[k] = -g_z_0_zz_yyyyyy[k] * ab_x + g_z_0_zz_xyyyyyy[k];

                g_z_0_xzz_yyyyyz[k] = -g_z_0_zz_yyyyyz[k] * ab_x + g_z_0_zz_xyyyyyz[k];

                g_z_0_xzz_yyyyzz[k] = -g_z_0_zz_yyyyzz[k] * ab_x + g_z_0_zz_xyyyyzz[k];

                g_z_0_xzz_yyyzzz[k] = -g_z_0_zz_yyyzzz[k] * ab_x + g_z_0_zz_xyyyzzz[k];

                g_z_0_xzz_yyzzzz[k] = -g_z_0_zz_yyzzzz[k] * ab_x + g_z_0_zz_xyyzzzz[k];

                g_z_0_xzz_yzzzzz[k] = -g_z_0_zz_yzzzzz[k] * ab_x + g_z_0_zz_xyzzzzz[k];

                g_z_0_xzz_zzzzzz[k] = -g_z_0_zz_zzzzzz[k] * ab_x + g_z_0_zz_xzzzzzz[k];
            }

            /// Set up 728-756 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyy_xxxxxx = cbuffer.data(fi_geom_10_off + 728 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxy = cbuffer.data(fi_geom_10_off + 729 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxxz = cbuffer.data(fi_geom_10_off + 730 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxyy = cbuffer.data(fi_geom_10_off + 731 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxyz = cbuffer.data(fi_geom_10_off + 732 * ccomps * dcomps);

            auto g_z_0_yyy_xxxxzz = cbuffer.data(fi_geom_10_off + 733 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyyy = cbuffer.data(fi_geom_10_off + 734 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyyz = cbuffer.data(fi_geom_10_off + 735 * ccomps * dcomps);

            auto g_z_0_yyy_xxxyzz = cbuffer.data(fi_geom_10_off + 736 * ccomps * dcomps);

            auto g_z_0_yyy_xxxzzz = cbuffer.data(fi_geom_10_off + 737 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyyy = cbuffer.data(fi_geom_10_off + 738 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyyz = cbuffer.data(fi_geom_10_off + 739 * ccomps * dcomps);

            auto g_z_0_yyy_xxyyzz = cbuffer.data(fi_geom_10_off + 740 * ccomps * dcomps);

            auto g_z_0_yyy_xxyzzz = cbuffer.data(fi_geom_10_off + 741 * ccomps * dcomps);

            auto g_z_0_yyy_xxzzzz = cbuffer.data(fi_geom_10_off + 742 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyyy = cbuffer.data(fi_geom_10_off + 743 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyyz = cbuffer.data(fi_geom_10_off + 744 * ccomps * dcomps);

            auto g_z_0_yyy_xyyyzz = cbuffer.data(fi_geom_10_off + 745 * ccomps * dcomps);

            auto g_z_0_yyy_xyyzzz = cbuffer.data(fi_geom_10_off + 746 * ccomps * dcomps);

            auto g_z_0_yyy_xyzzzz = cbuffer.data(fi_geom_10_off + 747 * ccomps * dcomps);

            auto g_z_0_yyy_xzzzzz = cbuffer.data(fi_geom_10_off + 748 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyyy = cbuffer.data(fi_geom_10_off + 749 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyyz = cbuffer.data(fi_geom_10_off + 750 * ccomps * dcomps);

            auto g_z_0_yyy_yyyyzz = cbuffer.data(fi_geom_10_off + 751 * ccomps * dcomps);

            auto g_z_0_yyy_yyyzzz = cbuffer.data(fi_geom_10_off + 752 * ccomps * dcomps);

            auto g_z_0_yyy_yyzzzz = cbuffer.data(fi_geom_10_off + 753 * ccomps * dcomps);

            auto g_z_0_yyy_yzzzzz = cbuffer.data(fi_geom_10_off + 754 * ccomps * dcomps);

            auto g_z_0_yyy_zzzzzz = cbuffer.data(fi_geom_10_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yy_xxxxxx, g_z_0_yy_xxxxxxy, g_z_0_yy_xxxxxy, g_z_0_yy_xxxxxyy, g_z_0_yy_xxxxxyz, g_z_0_yy_xxxxxz, g_z_0_yy_xxxxyy, g_z_0_yy_xxxxyyy, g_z_0_yy_xxxxyyz, g_z_0_yy_xxxxyz, g_z_0_yy_xxxxyzz, g_z_0_yy_xxxxzz, g_z_0_yy_xxxyyy, g_z_0_yy_xxxyyyy, g_z_0_yy_xxxyyyz, g_z_0_yy_xxxyyz, g_z_0_yy_xxxyyzz, g_z_0_yy_xxxyzz, g_z_0_yy_xxxyzzz, g_z_0_yy_xxxzzz, g_z_0_yy_xxyyyy, g_z_0_yy_xxyyyyy, g_z_0_yy_xxyyyyz, g_z_0_yy_xxyyyz, g_z_0_yy_xxyyyzz, g_z_0_yy_xxyyzz, g_z_0_yy_xxyyzzz, g_z_0_yy_xxyzzz, g_z_0_yy_xxyzzzz, g_z_0_yy_xxzzzz, g_z_0_yy_xyyyyy, g_z_0_yy_xyyyyyy, g_z_0_yy_xyyyyyz, g_z_0_yy_xyyyyz, g_z_0_yy_xyyyyzz, g_z_0_yy_xyyyzz, g_z_0_yy_xyyyzzz, g_z_0_yy_xyyzzz, g_z_0_yy_xyyzzzz, g_z_0_yy_xyzzzz, g_z_0_yy_xyzzzzz, g_z_0_yy_xzzzzz, g_z_0_yy_yyyyyy, g_z_0_yy_yyyyyyy, g_z_0_yy_yyyyyyz, g_z_0_yy_yyyyyz, g_z_0_yy_yyyyyzz, g_z_0_yy_yyyyzz, g_z_0_yy_yyyyzzz, g_z_0_yy_yyyzzz, g_z_0_yy_yyyzzzz, g_z_0_yy_yyzzzz, g_z_0_yy_yyzzzzz, g_z_0_yy_yzzzzz, g_z_0_yy_yzzzzzz, g_z_0_yy_zzzzzz, g_z_0_yyy_xxxxxx, g_z_0_yyy_xxxxxy, g_z_0_yyy_xxxxxz, g_z_0_yyy_xxxxyy, g_z_0_yyy_xxxxyz, g_z_0_yyy_xxxxzz, g_z_0_yyy_xxxyyy, g_z_0_yyy_xxxyyz, g_z_0_yyy_xxxyzz, g_z_0_yyy_xxxzzz, g_z_0_yyy_xxyyyy, g_z_0_yyy_xxyyyz, g_z_0_yyy_xxyyzz, g_z_0_yyy_xxyzzz, g_z_0_yyy_xxzzzz, g_z_0_yyy_xyyyyy, g_z_0_yyy_xyyyyz, g_z_0_yyy_xyyyzz, g_z_0_yyy_xyyzzz, g_z_0_yyy_xyzzzz, g_z_0_yyy_xzzzzz, g_z_0_yyy_yyyyyy, g_z_0_yyy_yyyyyz, g_z_0_yyy_yyyyzz, g_z_0_yyy_yyyzzz, g_z_0_yyy_yyzzzz, g_z_0_yyy_yzzzzz, g_z_0_yyy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyy_xxxxxx[k] = -g_z_0_yy_xxxxxx[k] * ab_y + g_z_0_yy_xxxxxxy[k];

                g_z_0_yyy_xxxxxy[k] = -g_z_0_yy_xxxxxy[k] * ab_y + g_z_0_yy_xxxxxyy[k];

                g_z_0_yyy_xxxxxz[k] = -g_z_0_yy_xxxxxz[k] * ab_y + g_z_0_yy_xxxxxyz[k];

                g_z_0_yyy_xxxxyy[k] = -g_z_0_yy_xxxxyy[k] * ab_y + g_z_0_yy_xxxxyyy[k];

                g_z_0_yyy_xxxxyz[k] = -g_z_0_yy_xxxxyz[k] * ab_y + g_z_0_yy_xxxxyyz[k];

                g_z_0_yyy_xxxxzz[k] = -g_z_0_yy_xxxxzz[k] * ab_y + g_z_0_yy_xxxxyzz[k];

                g_z_0_yyy_xxxyyy[k] = -g_z_0_yy_xxxyyy[k] * ab_y + g_z_0_yy_xxxyyyy[k];

                g_z_0_yyy_xxxyyz[k] = -g_z_0_yy_xxxyyz[k] * ab_y + g_z_0_yy_xxxyyyz[k];

                g_z_0_yyy_xxxyzz[k] = -g_z_0_yy_xxxyzz[k] * ab_y + g_z_0_yy_xxxyyzz[k];

                g_z_0_yyy_xxxzzz[k] = -g_z_0_yy_xxxzzz[k] * ab_y + g_z_0_yy_xxxyzzz[k];

                g_z_0_yyy_xxyyyy[k] = -g_z_0_yy_xxyyyy[k] * ab_y + g_z_0_yy_xxyyyyy[k];

                g_z_0_yyy_xxyyyz[k] = -g_z_0_yy_xxyyyz[k] * ab_y + g_z_0_yy_xxyyyyz[k];

                g_z_0_yyy_xxyyzz[k] = -g_z_0_yy_xxyyzz[k] * ab_y + g_z_0_yy_xxyyyzz[k];

                g_z_0_yyy_xxyzzz[k] = -g_z_0_yy_xxyzzz[k] * ab_y + g_z_0_yy_xxyyzzz[k];

                g_z_0_yyy_xxzzzz[k] = -g_z_0_yy_xxzzzz[k] * ab_y + g_z_0_yy_xxyzzzz[k];

                g_z_0_yyy_xyyyyy[k] = -g_z_0_yy_xyyyyy[k] * ab_y + g_z_0_yy_xyyyyyy[k];

                g_z_0_yyy_xyyyyz[k] = -g_z_0_yy_xyyyyz[k] * ab_y + g_z_0_yy_xyyyyyz[k];

                g_z_0_yyy_xyyyzz[k] = -g_z_0_yy_xyyyzz[k] * ab_y + g_z_0_yy_xyyyyzz[k];

                g_z_0_yyy_xyyzzz[k] = -g_z_0_yy_xyyzzz[k] * ab_y + g_z_0_yy_xyyyzzz[k];

                g_z_0_yyy_xyzzzz[k] = -g_z_0_yy_xyzzzz[k] * ab_y + g_z_0_yy_xyyzzzz[k];

                g_z_0_yyy_xzzzzz[k] = -g_z_0_yy_xzzzzz[k] * ab_y + g_z_0_yy_xyzzzzz[k];

                g_z_0_yyy_yyyyyy[k] = -g_z_0_yy_yyyyyy[k] * ab_y + g_z_0_yy_yyyyyyy[k];

                g_z_0_yyy_yyyyyz[k] = -g_z_0_yy_yyyyyz[k] * ab_y + g_z_0_yy_yyyyyyz[k];

                g_z_0_yyy_yyyyzz[k] = -g_z_0_yy_yyyyzz[k] * ab_y + g_z_0_yy_yyyyyzz[k];

                g_z_0_yyy_yyyzzz[k] = -g_z_0_yy_yyyzzz[k] * ab_y + g_z_0_yy_yyyyzzz[k];

                g_z_0_yyy_yyzzzz[k] = -g_z_0_yy_yyzzzz[k] * ab_y + g_z_0_yy_yyyzzzz[k];

                g_z_0_yyy_yzzzzz[k] = -g_z_0_yy_yzzzzz[k] * ab_y + g_z_0_yy_yyzzzzz[k];

                g_z_0_yyy_zzzzzz[k] = -g_z_0_yy_zzzzzz[k] * ab_y + g_z_0_yy_yzzzzzz[k];
            }

            /// Set up 756-784 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyz_xxxxxx = cbuffer.data(fi_geom_10_off + 756 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxy = cbuffer.data(fi_geom_10_off + 757 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxxz = cbuffer.data(fi_geom_10_off + 758 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxyy = cbuffer.data(fi_geom_10_off + 759 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxyz = cbuffer.data(fi_geom_10_off + 760 * ccomps * dcomps);

            auto g_z_0_yyz_xxxxzz = cbuffer.data(fi_geom_10_off + 761 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyyy = cbuffer.data(fi_geom_10_off + 762 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyyz = cbuffer.data(fi_geom_10_off + 763 * ccomps * dcomps);

            auto g_z_0_yyz_xxxyzz = cbuffer.data(fi_geom_10_off + 764 * ccomps * dcomps);

            auto g_z_0_yyz_xxxzzz = cbuffer.data(fi_geom_10_off + 765 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyyy = cbuffer.data(fi_geom_10_off + 766 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyyz = cbuffer.data(fi_geom_10_off + 767 * ccomps * dcomps);

            auto g_z_0_yyz_xxyyzz = cbuffer.data(fi_geom_10_off + 768 * ccomps * dcomps);

            auto g_z_0_yyz_xxyzzz = cbuffer.data(fi_geom_10_off + 769 * ccomps * dcomps);

            auto g_z_0_yyz_xxzzzz = cbuffer.data(fi_geom_10_off + 770 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyyy = cbuffer.data(fi_geom_10_off + 771 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyyz = cbuffer.data(fi_geom_10_off + 772 * ccomps * dcomps);

            auto g_z_0_yyz_xyyyzz = cbuffer.data(fi_geom_10_off + 773 * ccomps * dcomps);

            auto g_z_0_yyz_xyyzzz = cbuffer.data(fi_geom_10_off + 774 * ccomps * dcomps);

            auto g_z_0_yyz_xyzzzz = cbuffer.data(fi_geom_10_off + 775 * ccomps * dcomps);

            auto g_z_0_yyz_xzzzzz = cbuffer.data(fi_geom_10_off + 776 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyyy = cbuffer.data(fi_geom_10_off + 777 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyyz = cbuffer.data(fi_geom_10_off + 778 * ccomps * dcomps);

            auto g_z_0_yyz_yyyyzz = cbuffer.data(fi_geom_10_off + 779 * ccomps * dcomps);

            auto g_z_0_yyz_yyyzzz = cbuffer.data(fi_geom_10_off + 780 * ccomps * dcomps);

            auto g_z_0_yyz_yyzzzz = cbuffer.data(fi_geom_10_off + 781 * ccomps * dcomps);

            auto g_z_0_yyz_yzzzzz = cbuffer.data(fi_geom_10_off + 782 * ccomps * dcomps);

            auto g_z_0_yyz_zzzzzz = cbuffer.data(fi_geom_10_off + 783 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyz_xxxxxx, g_z_0_yyz_xxxxxy, g_z_0_yyz_xxxxxz, g_z_0_yyz_xxxxyy, g_z_0_yyz_xxxxyz, g_z_0_yyz_xxxxzz, g_z_0_yyz_xxxyyy, g_z_0_yyz_xxxyyz, g_z_0_yyz_xxxyzz, g_z_0_yyz_xxxzzz, g_z_0_yyz_xxyyyy, g_z_0_yyz_xxyyyz, g_z_0_yyz_xxyyzz, g_z_0_yyz_xxyzzz, g_z_0_yyz_xxzzzz, g_z_0_yyz_xyyyyy, g_z_0_yyz_xyyyyz, g_z_0_yyz_xyyyzz, g_z_0_yyz_xyyzzz, g_z_0_yyz_xyzzzz, g_z_0_yyz_xzzzzz, g_z_0_yyz_yyyyyy, g_z_0_yyz_yyyyyz, g_z_0_yyz_yyyyzz, g_z_0_yyz_yyyzzz, g_z_0_yyz_yyzzzz, g_z_0_yyz_yzzzzz, g_z_0_yyz_zzzzzz, g_z_0_yz_xxxxxx, g_z_0_yz_xxxxxxy, g_z_0_yz_xxxxxy, g_z_0_yz_xxxxxyy, g_z_0_yz_xxxxxyz, g_z_0_yz_xxxxxz, g_z_0_yz_xxxxyy, g_z_0_yz_xxxxyyy, g_z_0_yz_xxxxyyz, g_z_0_yz_xxxxyz, g_z_0_yz_xxxxyzz, g_z_0_yz_xxxxzz, g_z_0_yz_xxxyyy, g_z_0_yz_xxxyyyy, g_z_0_yz_xxxyyyz, g_z_0_yz_xxxyyz, g_z_0_yz_xxxyyzz, g_z_0_yz_xxxyzz, g_z_0_yz_xxxyzzz, g_z_0_yz_xxxzzz, g_z_0_yz_xxyyyy, g_z_0_yz_xxyyyyy, g_z_0_yz_xxyyyyz, g_z_0_yz_xxyyyz, g_z_0_yz_xxyyyzz, g_z_0_yz_xxyyzz, g_z_0_yz_xxyyzzz, g_z_0_yz_xxyzzz, g_z_0_yz_xxyzzzz, g_z_0_yz_xxzzzz, g_z_0_yz_xyyyyy, g_z_0_yz_xyyyyyy, g_z_0_yz_xyyyyyz, g_z_0_yz_xyyyyz, g_z_0_yz_xyyyyzz, g_z_0_yz_xyyyzz, g_z_0_yz_xyyyzzz, g_z_0_yz_xyyzzz, g_z_0_yz_xyyzzzz, g_z_0_yz_xyzzzz, g_z_0_yz_xyzzzzz, g_z_0_yz_xzzzzz, g_z_0_yz_yyyyyy, g_z_0_yz_yyyyyyy, g_z_0_yz_yyyyyyz, g_z_0_yz_yyyyyz, g_z_0_yz_yyyyyzz, g_z_0_yz_yyyyzz, g_z_0_yz_yyyyzzz, g_z_0_yz_yyyzzz, g_z_0_yz_yyyzzzz, g_z_0_yz_yyzzzz, g_z_0_yz_yyzzzzz, g_z_0_yz_yzzzzz, g_z_0_yz_yzzzzzz, g_z_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyz_xxxxxx[k] = -g_z_0_yz_xxxxxx[k] * ab_y + g_z_0_yz_xxxxxxy[k];

                g_z_0_yyz_xxxxxy[k] = -g_z_0_yz_xxxxxy[k] * ab_y + g_z_0_yz_xxxxxyy[k];

                g_z_0_yyz_xxxxxz[k] = -g_z_0_yz_xxxxxz[k] * ab_y + g_z_0_yz_xxxxxyz[k];

                g_z_0_yyz_xxxxyy[k] = -g_z_0_yz_xxxxyy[k] * ab_y + g_z_0_yz_xxxxyyy[k];

                g_z_0_yyz_xxxxyz[k] = -g_z_0_yz_xxxxyz[k] * ab_y + g_z_0_yz_xxxxyyz[k];

                g_z_0_yyz_xxxxzz[k] = -g_z_0_yz_xxxxzz[k] * ab_y + g_z_0_yz_xxxxyzz[k];

                g_z_0_yyz_xxxyyy[k] = -g_z_0_yz_xxxyyy[k] * ab_y + g_z_0_yz_xxxyyyy[k];

                g_z_0_yyz_xxxyyz[k] = -g_z_0_yz_xxxyyz[k] * ab_y + g_z_0_yz_xxxyyyz[k];

                g_z_0_yyz_xxxyzz[k] = -g_z_0_yz_xxxyzz[k] * ab_y + g_z_0_yz_xxxyyzz[k];

                g_z_0_yyz_xxxzzz[k] = -g_z_0_yz_xxxzzz[k] * ab_y + g_z_0_yz_xxxyzzz[k];

                g_z_0_yyz_xxyyyy[k] = -g_z_0_yz_xxyyyy[k] * ab_y + g_z_0_yz_xxyyyyy[k];

                g_z_0_yyz_xxyyyz[k] = -g_z_0_yz_xxyyyz[k] * ab_y + g_z_0_yz_xxyyyyz[k];

                g_z_0_yyz_xxyyzz[k] = -g_z_0_yz_xxyyzz[k] * ab_y + g_z_0_yz_xxyyyzz[k];

                g_z_0_yyz_xxyzzz[k] = -g_z_0_yz_xxyzzz[k] * ab_y + g_z_0_yz_xxyyzzz[k];

                g_z_0_yyz_xxzzzz[k] = -g_z_0_yz_xxzzzz[k] * ab_y + g_z_0_yz_xxyzzzz[k];

                g_z_0_yyz_xyyyyy[k] = -g_z_0_yz_xyyyyy[k] * ab_y + g_z_0_yz_xyyyyyy[k];

                g_z_0_yyz_xyyyyz[k] = -g_z_0_yz_xyyyyz[k] * ab_y + g_z_0_yz_xyyyyyz[k];

                g_z_0_yyz_xyyyzz[k] = -g_z_0_yz_xyyyzz[k] * ab_y + g_z_0_yz_xyyyyzz[k];

                g_z_0_yyz_xyyzzz[k] = -g_z_0_yz_xyyzzz[k] * ab_y + g_z_0_yz_xyyyzzz[k];

                g_z_0_yyz_xyzzzz[k] = -g_z_0_yz_xyzzzz[k] * ab_y + g_z_0_yz_xyyzzzz[k];

                g_z_0_yyz_xzzzzz[k] = -g_z_0_yz_xzzzzz[k] * ab_y + g_z_0_yz_xyzzzzz[k];

                g_z_0_yyz_yyyyyy[k] = -g_z_0_yz_yyyyyy[k] * ab_y + g_z_0_yz_yyyyyyy[k];

                g_z_0_yyz_yyyyyz[k] = -g_z_0_yz_yyyyyz[k] * ab_y + g_z_0_yz_yyyyyyz[k];

                g_z_0_yyz_yyyyzz[k] = -g_z_0_yz_yyyyzz[k] * ab_y + g_z_0_yz_yyyyyzz[k];

                g_z_0_yyz_yyyzzz[k] = -g_z_0_yz_yyyzzz[k] * ab_y + g_z_0_yz_yyyyzzz[k];

                g_z_0_yyz_yyzzzz[k] = -g_z_0_yz_yyzzzz[k] * ab_y + g_z_0_yz_yyyzzzz[k];

                g_z_0_yyz_yzzzzz[k] = -g_z_0_yz_yzzzzz[k] * ab_y + g_z_0_yz_yyzzzzz[k];

                g_z_0_yyz_zzzzzz[k] = -g_z_0_yz_zzzzzz[k] * ab_y + g_z_0_yz_yzzzzzz[k];
            }

            /// Set up 784-812 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzz_xxxxxx = cbuffer.data(fi_geom_10_off + 784 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxy = cbuffer.data(fi_geom_10_off + 785 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxxz = cbuffer.data(fi_geom_10_off + 786 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxyy = cbuffer.data(fi_geom_10_off + 787 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxyz = cbuffer.data(fi_geom_10_off + 788 * ccomps * dcomps);

            auto g_z_0_yzz_xxxxzz = cbuffer.data(fi_geom_10_off + 789 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyyy = cbuffer.data(fi_geom_10_off + 790 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyyz = cbuffer.data(fi_geom_10_off + 791 * ccomps * dcomps);

            auto g_z_0_yzz_xxxyzz = cbuffer.data(fi_geom_10_off + 792 * ccomps * dcomps);

            auto g_z_0_yzz_xxxzzz = cbuffer.data(fi_geom_10_off + 793 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyyy = cbuffer.data(fi_geom_10_off + 794 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyyz = cbuffer.data(fi_geom_10_off + 795 * ccomps * dcomps);

            auto g_z_0_yzz_xxyyzz = cbuffer.data(fi_geom_10_off + 796 * ccomps * dcomps);

            auto g_z_0_yzz_xxyzzz = cbuffer.data(fi_geom_10_off + 797 * ccomps * dcomps);

            auto g_z_0_yzz_xxzzzz = cbuffer.data(fi_geom_10_off + 798 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyyy = cbuffer.data(fi_geom_10_off + 799 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyyz = cbuffer.data(fi_geom_10_off + 800 * ccomps * dcomps);

            auto g_z_0_yzz_xyyyzz = cbuffer.data(fi_geom_10_off + 801 * ccomps * dcomps);

            auto g_z_0_yzz_xyyzzz = cbuffer.data(fi_geom_10_off + 802 * ccomps * dcomps);

            auto g_z_0_yzz_xyzzzz = cbuffer.data(fi_geom_10_off + 803 * ccomps * dcomps);

            auto g_z_0_yzz_xzzzzz = cbuffer.data(fi_geom_10_off + 804 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyyy = cbuffer.data(fi_geom_10_off + 805 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyyz = cbuffer.data(fi_geom_10_off + 806 * ccomps * dcomps);

            auto g_z_0_yzz_yyyyzz = cbuffer.data(fi_geom_10_off + 807 * ccomps * dcomps);

            auto g_z_0_yzz_yyyzzz = cbuffer.data(fi_geom_10_off + 808 * ccomps * dcomps);

            auto g_z_0_yzz_yyzzzz = cbuffer.data(fi_geom_10_off + 809 * ccomps * dcomps);

            auto g_z_0_yzz_yzzzzz = cbuffer.data(fi_geom_10_off + 810 * ccomps * dcomps);

            auto g_z_0_yzz_zzzzzz = cbuffer.data(fi_geom_10_off + 811 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzz_xxxxxx, g_z_0_yzz_xxxxxy, g_z_0_yzz_xxxxxz, g_z_0_yzz_xxxxyy, g_z_0_yzz_xxxxyz, g_z_0_yzz_xxxxzz, g_z_0_yzz_xxxyyy, g_z_0_yzz_xxxyyz, g_z_0_yzz_xxxyzz, g_z_0_yzz_xxxzzz, g_z_0_yzz_xxyyyy, g_z_0_yzz_xxyyyz, g_z_0_yzz_xxyyzz, g_z_0_yzz_xxyzzz, g_z_0_yzz_xxzzzz, g_z_0_yzz_xyyyyy, g_z_0_yzz_xyyyyz, g_z_0_yzz_xyyyzz, g_z_0_yzz_xyyzzz, g_z_0_yzz_xyzzzz, g_z_0_yzz_xzzzzz, g_z_0_yzz_yyyyyy, g_z_0_yzz_yyyyyz, g_z_0_yzz_yyyyzz, g_z_0_yzz_yyyzzz, g_z_0_yzz_yyzzzz, g_z_0_yzz_yzzzzz, g_z_0_yzz_zzzzzz, g_z_0_zz_xxxxxx, g_z_0_zz_xxxxxxy, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxxyy, g_z_0_zz_xxxxxyz, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyyy, g_z_0_zz_xxxxyyz, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxyzz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyyy, g_z_0_zz_xxxyyyz, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyyzz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxyzzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyyy, g_z_0_zz_xxyyyyz, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyyzz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyyzzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxyzzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyyy, g_z_0_zz_xyyyyyz, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyyzz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyyzzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyyzzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xyzzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_yyyyyy, g_z_0_zz_yyyyyyy, g_z_0_zz_yyyyyyz, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyyzz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyyzzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyyzzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yyzzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_yzzzzzz, g_z_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzz_xxxxxx[k] = -g_z_0_zz_xxxxxx[k] * ab_y + g_z_0_zz_xxxxxxy[k];

                g_z_0_yzz_xxxxxy[k] = -g_z_0_zz_xxxxxy[k] * ab_y + g_z_0_zz_xxxxxyy[k];

                g_z_0_yzz_xxxxxz[k] = -g_z_0_zz_xxxxxz[k] * ab_y + g_z_0_zz_xxxxxyz[k];

                g_z_0_yzz_xxxxyy[k] = -g_z_0_zz_xxxxyy[k] * ab_y + g_z_0_zz_xxxxyyy[k];

                g_z_0_yzz_xxxxyz[k] = -g_z_0_zz_xxxxyz[k] * ab_y + g_z_0_zz_xxxxyyz[k];

                g_z_0_yzz_xxxxzz[k] = -g_z_0_zz_xxxxzz[k] * ab_y + g_z_0_zz_xxxxyzz[k];

                g_z_0_yzz_xxxyyy[k] = -g_z_0_zz_xxxyyy[k] * ab_y + g_z_0_zz_xxxyyyy[k];

                g_z_0_yzz_xxxyyz[k] = -g_z_0_zz_xxxyyz[k] * ab_y + g_z_0_zz_xxxyyyz[k];

                g_z_0_yzz_xxxyzz[k] = -g_z_0_zz_xxxyzz[k] * ab_y + g_z_0_zz_xxxyyzz[k];

                g_z_0_yzz_xxxzzz[k] = -g_z_0_zz_xxxzzz[k] * ab_y + g_z_0_zz_xxxyzzz[k];

                g_z_0_yzz_xxyyyy[k] = -g_z_0_zz_xxyyyy[k] * ab_y + g_z_0_zz_xxyyyyy[k];

                g_z_0_yzz_xxyyyz[k] = -g_z_0_zz_xxyyyz[k] * ab_y + g_z_0_zz_xxyyyyz[k];

                g_z_0_yzz_xxyyzz[k] = -g_z_0_zz_xxyyzz[k] * ab_y + g_z_0_zz_xxyyyzz[k];

                g_z_0_yzz_xxyzzz[k] = -g_z_0_zz_xxyzzz[k] * ab_y + g_z_0_zz_xxyyzzz[k];

                g_z_0_yzz_xxzzzz[k] = -g_z_0_zz_xxzzzz[k] * ab_y + g_z_0_zz_xxyzzzz[k];

                g_z_0_yzz_xyyyyy[k] = -g_z_0_zz_xyyyyy[k] * ab_y + g_z_0_zz_xyyyyyy[k];

                g_z_0_yzz_xyyyyz[k] = -g_z_0_zz_xyyyyz[k] * ab_y + g_z_0_zz_xyyyyyz[k];

                g_z_0_yzz_xyyyzz[k] = -g_z_0_zz_xyyyzz[k] * ab_y + g_z_0_zz_xyyyyzz[k];

                g_z_0_yzz_xyyzzz[k] = -g_z_0_zz_xyyzzz[k] * ab_y + g_z_0_zz_xyyyzzz[k];

                g_z_0_yzz_xyzzzz[k] = -g_z_0_zz_xyzzzz[k] * ab_y + g_z_0_zz_xyyzzzz[k];

                g_z_0_yzz_xzzzzz[k] = -g_z_0_zz_xzzzzz[k] * ab_y + g_z_0_zz_xyzzzzz[k];

                g_z_0_yzz_yyyyyy[k] = -g_z_0_zz_yyyyyy[k] * ab_y + g_z_0_zz_yyyyyyy[k];

                g_z_0_yzz_yyyyyz[k] = -g_z_0_zz_yyyyyz[k] * ab_y + g_z_0_zz_yyyyyyz[k];

                g_z_0_yzz_yyyyzz[k] = -g_z_0_zz_yyyyzz[k] * ab_y + g_z_0_zz_yyyyyzz[k];

                g_z_0_yzz_yyyzzz[k] = -g_z_0_zz_yyyzzz[k] * ab_y + g_z_0_zz_yyyyzzz[k];

                g_z_0_yzz_yyzzzz[k] = -g_z_0_zz_yyzzzz[k] * ab_y + g_z_0_zz_yyyzzzz[k];

                g_z_0_yzz_yzzzzz[k] = -g_z_0_zz_yzzzzz[k] * ab_y + g_z_0_zz_yyzzzzz[k];

                g_z_0_yzz_zzzzzz[k] = -g_z_0_zz_zzzzzz[k] * ab_y + g_z_0_zz_yzzzzzz[k];
            }

            /// Set up 812-840 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzz_xxxxxx = cbuffer.data(fi_geom_10_off + 812 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxy = cbuffer.data(fi_geom_10_off + 813 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxxz = cbuffer.data(fi_geom_10_off + 814 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxyy = cbuffer.data(fi_geom_10_off + 815 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxyz = cbuffer.data(fi_geom_10_off + 816 * ccomps * dcomps);

            auto g_z_0_zzz_xxxxzz = cbuffer.data(fi_geom_10_off + 817 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyyy = cbuffer.data(fi_geom_10_off + 818 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyyz = cbuffer.data(fi_geom_10_off + 819 * ccomps * dcomps);

            auto g_z_0_zzz_xxxyzz = cbuffer.data(fi_geom_10_off + 820 * ccomps * dcomps);

            auto g_z_0_zzz_xxxzzz = cbuffer.data(fi_geom_10_off + 821 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyyy = cbuffer.data(fi_geom_10_off + 822 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyyz = cbuffer.data(fi_geom_10_off + 823 * ccomps * dcomps);

            auto g_z_0_zzz_xxyyzz = cbuffer.data(fi_geom_10_off + 824 * ccomps * dcomps);

            auto g_z_0_zzz_xxyzzz = cbuffer.data(fi_geom_10_off + 825 * ccomps * dcomps);

            auto g_z_0_zzz_xxzzzz = cbuffer.data(fi_geom_10_off + 826 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyyy = cbuffer.data(fi_geom_10_off + 827 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyyz = cbuffer.data(fi_geom_10_off + 828 * ccomps * dcomps);

            auto g_z_0_zzz_xyyyzz = cbuffer.data(fi_geom_10_off + 829 * ccomps * dcomps);

            auto g_z_0_zzz_xyyzzz = cbuffer.data(fi_geom_10_off + 830 * ccomps * dcomps);

            auto g_z_0_zzz_xyzzzz = cbuffer.data(fi_geom_10_off + 831 * ccomps * dcomps);

            auto g_z_0_zzz_xzzzzz = cbuffer.data(fi_geom_10_off + 832 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyyy = cbuffer.data(fi_geom_10_off + 833 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyyz = cbuffer.data(fi_geom_10_off + 834 * ccomps * dcomps);

            auto g_z_0_zzz_yyyyzz = cbuffer.data(fi_geom_10_off + 835 * ccomps * dcomps);

            auto g_z_0_zzz_yyyzzz = cbuffer.data(fi_geom_10_off + 836 * ccomps * dcomps);

            auto g_z_0_zzz_yyzzzz = cbuffer.data(fi_geom_10_off + 837 * ccomps * dcomps);

            auto g_z_0_zzz_yzzzzz = cbuffer.data(fi_geom_10_off + 838 * ccomps * dcomps);

            auto g_z_0_zzz_zzzzzz = cbuffer.data(fi_geom_10_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zz_xxxxxx, g_z_0_zz_xxxxxxz, g_z_0_zz_xxxxxy, g_z_0_zz_xxxxxyz, g_z_0_zz_xxxxxz, g_z_0_zz_xxxxxzz, g_z_0_zz_xxxxyy, g_z_0_zz_xxxxyyz, g_z_0_zz_xxxxyz, g_z_0_zz_xxxxyzz, g_z_0_zz_xxxxzz, g_z_0_zz_xxxxzzz, g_z_0_zz_xxxyyy, g_z_0_zz_xxxyyyz, g_z_0_zz_xxxyyz, g_z_0_zz_xxxyyzz, g_z_0_zz_xxxyzz, g_z_0_zz_xxxyzzz, g_z_0_zz_xxxzzz, g_z_0_zz_xxxzzzz, g_z_0_zz_xxyyyy, g_z_0_zz_xxyyyyz, g_z_0_zz_xxyyyz, g_z_0_zz_xxyyyzz, g_z_0_zz_xxyyzz, g_z_0_zz_xxyyzzz, g_z_0_zz_xxyzzz, g_z_0_zz_xxyzzzz, g_z_0_zz_xxzzzz, g_z_0_zz_xxzzzzz, g_z_0_zz_xyyyyy, g_z_0_zz_xyyyyyz, g_z_0_zz_xyyyyz, g_z_0_zz_xyyyyzz, g_z_0_zz_xyyyzz, g_z_0_zz_xyyyzzz, g_z_0_zz_xyyzzz, g_z_0_zz_xyyzzzz, g_z_0_zz_xyzzzz, g_z_0_zz_xyzzzzz, g_z_0_zz_xzzzzz, g_z_0_zz_xzzzzzz, g_z_0_zz_yyyyyy, g_z_0_zz_yyyyyyz, g_z_0_zz_yyyyyz, g_z_0_zz_yyyyyzz, g_z_0_zz_yyyyzz, g_z_0_zz_yyyyzzz, g_z_0_zz_yyyzzz, g_z_0_zz_yyyzzzz, g_z_0_zz_yyzzzz, g_z_0_zz_yyzzzzz, g_z_0_zz_yzzzzz, g_z_0_zz_yzzzzzz, g_z_0_zz_zzzzzz, g_z_0_zz_zzzzzzz, g_z_0_zzz_xxxxxx, g_z_0_zzz_xxxxxy, g_z_0_zzz_xxxxxz, g_z_0_zzz_xxxxyy, g_z_0_zzz_xxxxyz, g_z_0_zzz_xxxxzz, g_z_0_zzz_xxxyyy, g_z_0_zzz_xxxyyz, g_z_0_zzz_xxxyzz, g_z_0_zzz_xxxzzz, g_z_0_zzz_xxyyyy, g_z_0_zzz_xxyyyz, g_z_0_zzz_xxyyzz, g_z_0_zzz_xxyzzz, g_z_0_zzz_xxzzzz, g_z_0_zzz_xyyyyy, g_z_0_zzz_xyyyyz, g_z_0_zzz_xyyyzz, g_z_0_zzz_xyyzzz, g_z_0_zzz_xyzzzz, g_z_0_zzz_xzzzzz, g_z_0_zzz_yyyyyy, g_z_0_zzz_yyyyyz, g_z_0_zzz_yyyyzz, g_z_0_zzz_yyyzzz, g_z_0_zzz_yyzzzz, g_z_0_zzz_yzzzzz, g_z_0_zzz_zzzzzz, g_zz_xxxxxx, g_zz_xxxxxy, g_zz_xxxxxz, g_zz_xxxxyy, g_zz_xxxxyz, g_zz_xxxxzz, g_zz_xxxyyy, g_zz_xxxyyz, g_zz_xxxyzz, g_zz_xxxzzz, g_zz_xxyyyy, g_zz_xxyyyz, g_zz_xxyyzz, g_zz_xxyzzz, g_zz_xxzzzz, g_zz_xyyyyy, g_zz_xyyyyz, g_zz_xyyyzz, g_zz_xyyzzz, g_zz_xyzzzz, g_zz_xzzzzz, g_zz_yyyyyy, g_zz_yyyyyz, g_zz_yyyyzz, g_zz_yyyzzz, g_zz_yyzzzz, g_zz_yzzzzz, g_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzz_xxxxxx[k] = -g_zz_xxxxxx[k] - g_z_0_zz_xxxxxx[k] * ab_z + g_z_0_zz_xxxxxxz[k];

                g_z_0_zzz_xxxxxy[k] = -g_zz_xxxxxy[k] - g_z_0_zz_xxxxxy[k] * ab_z + g_z_0_zz_xxxxxyz[k];

                g_z_0_zzz_xxxxxz[k] = -g_zz_xxxxxz[k] - g_z_0_zz_xxxxxz[k] * ab_z + g_z_0_zz_xxxxxzz[k];

                g_z_0_zzz_xxxxyy[k] = -g_zz_xxxxyy[k] - g_z_0_zz_xxxxyy[k] * ab_z + g_z_0_zz_xxxxyyz[k];

                g_z_0_zzz_xxxxyz[k] = -g_zz_xxxxyz[k] - g_z_0_zz_xxxxyz[k] * ab_z + g_z_0_zz_xxxxyzz[k];

                g_z_0_zzz_xxxxzz[k] = -g_zz_xxxxzz[k] - g_z_0_zz_xxxxzz[k] * ab_z + g_z_0_zz_xxxxzzz[k];

                g_z_0_zzz_xxxyyy[k] = -g_zz_xxxyyy[k] - g_z_0_zz_xxxyyy[k] * ab_z + g_z_0_zz_xxxyyyz[k];

                g_z_0_zzz_xxxyyz[k] = -g_zz_xxxyyz[k] - g_z_0_zz_xxxyyz[k] * ab_z + g_z_0_zz_xxxyyzz[k];

                g_z_0_zzz_xxxyzz[k] = -g_zz_xxxyzz[k] - g_z_0_zz_xxxyzz[k] * ab_z + g_z_0_zz_xxxyzzz[k];

                g_z_0_zzz_xxxzzz[k] = -g_zz_xxxzzz[k] - g_z_0_zz_xxxzzz[k] * ab_z + g_z_0_zz_xxxzzzz[k];

                g_z_0_zzz_xxyyyy[k] = -g_zz_xxyyyy[k] - g_z_0_zz_xxyyyy[k] * ab_z + g_z_0_zz_xxyyyyz[k];

                g_z_0_zzz_xxyyyz[k] = -g_zz_xxyyyz[k] - g_z_0_zz_xxyyyz[k] * ab_z + g_z_0_zz_xxyyyzz[k];

                g_z_0_zzz_xxyyzz[k] = -g_zz_xxyyzz[k] - g_z_0_zz_xxyyzz[k] * ab_z + g_z_0_zz_xxyyzzz[k];

                g_z_0_zzz_xxyzzz[k] = -g_zz_xxyzzz[k] - g_z_0_zz_xxyzzz[k] * ab_z + g_z_0_zz_xxyzzzz[k];

                g_z_0_zzz_xxzzzz[k] = -g_zz_xxzzzz[k] - g_z_0_zz_xxzzzz[k] * ab_z + g_z_0_zz_xxzzzzz[k];

                g_z_0_zzz_xyyyyy[k] = -g_zz_xyyyyy[k] - g_z_0_zz_xyyyyy[k] * ab_z + g_z_0_zz_xyyyyyz[k];

                g_z_0_zzz_xyyyyz[k] = -g_zz_xyyyyz[k] - g_z_0_zz_xyyyyz[k] * ab_z + g_z_0_zz_xyyyyzz[k];

                g_z_0_zzz_xyyyzz[k] = -g_zz_xyyyzz[k] - g_z_0_zz_xyyyzz[k] * ab_z + g_z_0_zz_xyyyzzz[k];

                g_z_0_zzz_xyyzzz[k] = -g_zz_xyyzzz[k] - g_z_0_zz_xyyzzz[k] * ab_z + g_z_0_zz_xyyzzzz[k];

                g_z_0_zzz_xyzzzz[k] = -g_zz_xyzzzz[k] - g_z_0_zz_xyzzzz[k] * ab_z + g_z_0_zz_xyzzzzz[k];

                g_z_0_zzz_xzzzzz[k] = -g_zz_xzzzzz[k] - g_z_0_zz_xzzzzz[k] * ab_z + g_z_0_zz_xzzzzzz[k];

                g_z_0_zzz_yyyyyy[k] = -g_zz_yyyyyy[k] - g_z_0_zz_yyyyyy[k] * ab_z + g_z_0_zz_yyyyyyz[k];

                g_z_0_zzz_yyyyyz[k] = -g_zz_yyyyyz[k] - g_z_0_zz_yyyyyz[k] * ab_z + g_z_0_zz_yyyyyzz[k];

                g_z_0_zzz_yyyyzz[k] = -g_zz_yyyyzz[k] - g_z_0_zz_yyyyzz[k] * ab_z + g_z_0_zz_yyyyzzz[k];

                g_z_0_zzz_yyyzzz[k] = -g_zz_yyyzzz[k] - g_z_0_zz_yyyzzz[k] * ab_z + g_z_0_zz_yyyzzzz[k];

                g_z_0_zzz_yyzzzz[k] = -g_zz_yyzzzz[k] - g_z_0_zz_yyzzzz[k] * ab_z + g_z_0_zz_yyzzzzz[k];

                g_z_0_zzz_yzzzzz[k] = -g_zz_yzzzzz[k] - g_z_0_zz_yzzzzz[k] * ab_z + g_z_0_zz_yzzzzzz[k];

                g_z_0_zzz_zzzzzz[k] = -g_zz_zzzzzz[k] - g_z_0_zz_zzzzzz[k] * ab_z + g_z_0_zz_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

