#include "ElectronRepulsionGeom2000ContrRecDIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_dixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_dixx,
                                            const size_t idx_geom_10_pixx,
                                            const size_t idx_geom_20_pixx,
                                            const size_t idx_geom_20_pkxx,
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

            /// Set up components of auxilary buffer : PISS

            const auto pi_geom_20_off = idx_geom_20_pixx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : PKSS

            const auto pk_geom_20_off = idx_geom_20_pkxx + i * dcomps + j;

            auto g_xx_0_x_xxxxxxx = cbuffer.data(pk_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_x_xxxxxxy = cbuffer.data(pk_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_x_xxxxxxz = cbuffer.data(pk_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_x_xxxxxyy = cbuffer.data(pk_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_x_xxxxxyz = cbuffer.data(pk_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_x_xxxxxzz = cbuffer.data(pk_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_x_xxxxyyy = cbuffer.data(pk_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_x_xxxxyyz = cbuffer.data(pk_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_x_xxxxyzz = cbuffer.data(pk_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_x_xxxxzzz = cbuffer.data(pk_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_x_xxxyyyy = cbuffer.data(pk_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_x_xxxyyyz = cbuffer.data(pk_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_x_xxxyyzz = cbuffer.data(pk_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_x_xxxyzzz = cbuffer.data(pk_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_x_xxxzzzz = cbuffer.data(pk_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_x_xxyyyyy = cbuffer.data(pk_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_x_xxyyyyz = cbuffer.data(pk_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_x_xxyyyzz = cbuffer.data(pk_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_x_xxyyzzz = cbuffer.data(pk_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_x_xxyzzzz = cbuffer.data(pk_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_x_xxzzzzz = cbuffer.data(pk_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_x_xyyyyyy = cbuffer.data(pk_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_x_xyyyyyz = cbuffer.data(pk_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_x_xyyyyzz = cbuffer.data(pk_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_x_xyyyzzz = cbuffer.data(pk_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_x_xyyzzzz = cbuffer.data(pk_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_x_xyzzzzz = cbuffer.data(pk_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_x_xzzzzzz = cbuffer.data(pk_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_x_yyyyyyy = cbuffer.data(pk_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_x_yyyyyyz = cbuffer.data(pk_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_x_yyyyyzz = cbuffer.data(pk_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_x_yyyyzzz = cbuffer.data(pk_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_x_yyyzzzz = cbuffer.data(pk_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_x_yyzzzzz = cbuffer.data(pk_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_x_yzzzzzz = cbuffer.data(pk_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_x_zzzzzzz = cbuffer.data(pk_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_y_xxxxxxx = cbuffer.data(pk_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_y_xxxxxxy = cbuffer.data(pk_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_y_xxxxxxz = cbuffer.data(pk_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_y_xxxxxyy = cbuffer.data(pk_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_y_xxxxxyz = cbuffer.data(pk_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_y_xxxxxzz = cbuffer.data(pk_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_y_xxxxyyy = cbuffer.data(pk_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_y_xxxxyyz = cbuffer.data(pk_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_y_xxxxyzz = cbuffer.data(pk_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_y_xxxxzzz = cbuffer.data(pk_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_y_xxxyyyy = cbuffer.data(pk_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_y_xxxyyyz = cbuffer.data(pk_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_y_xxxyyzz = cbuffer.data(pk_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_y_xxxyzzz = cbuffer.data(pk_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_y_xxxzzzz = cbuffer.data(pk_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_y_xxyyyyy = cbuffer.data(pk_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_y_xxyyyyz = cbuffer.data(pk_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_y_xxyyyzz = cbuffer.data(pk_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_y_xxyyzzz = cbuffer.data(pk_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_y_xxyzzzz = cbuffer.data(pk_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_y_xxzzzzz = cbuffer.data(pk_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_y_xyyyyyy = cbuffer.data(pk_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_y_xyyyyyz = cbuffer.data(pk_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_y_xyyyyzz = cbuffer.data(pk_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_y_xyyyzzz = cbuffer.data(pk_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_y_xyyzzzz = cbuffer.data(pk_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_y_xyzzzzz = cbuffer.data(pk_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_y_xzzzzzz = cbuffer.data(pk_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_y_yyyyyyy = cbuffer.data(pk_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_y_yyyyyyz = cbuffer.data(pk_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_y_yyyyyzz = cbuffer.data(pk_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_y_yyyyzzz = cbuffer.data(pk_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_y_yyyzzzz = cbuffer.data(pk_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_y_yyzzzzz = cbuffer.data(pk_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_y_yzzzzzz = cbuffer.data(pk_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_y_zzzzzzz = cbuffer.data(pk_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_z_xxxxxxx = cbuffer.data(pk_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_z_xxxxxxy = cbuffer.data(pk_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_z_xxxxxxz = cbuffer.data(pk_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_z_xxxxxyy = cbuffer.data(pk_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_z_xxxxxyz = cbuffer.data(pk_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_z_xxxxxzz = cbuffer.data(pk_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_z_xxxxyyy = cbuffer.data(pk_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_z_xxxxyyz = cbuffer.data(pk_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_z_xxxxyzz = cbuffer.data(pk_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_z_xxxxzzz = cbuffer.data(pk_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_z_xxxyyyy = cbuffer.data(pk_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_z_xxxyyyz = cbuffer.data(pk_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_z_xxxyyzz = cbuffer.data(pk_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_z_xxxyzzz = cbuffer.data(pk_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_z_xxxzzzz = cbuffer.data(pk_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_z_xxyyyyy = cbuffer.data(pk_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_z_xxyyyyz = cbuffer.data(pk_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_z_xxyyyzz = cbuffer.data(pk_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_z_xxyyzzz = cbuffer.data(pk_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_z_xxyzzzz = cbuffer.data(pk_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_z_xxzzzzz = cbuffer.data(pk_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_z_xyyyyyy = cbuffer.data(pk_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_z_xyyyyyz = cbuffer.data(pk_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_z_xyyyyzz = cbuffer.data(pk_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_z_xyyyzzz = cbuffer.data(pk_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_z_xyyzzzz = cbuffer.data(pk_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_z_xyzzzzz = cbuffer.data(pk_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_z_xzzzzzz = cbuffer.data(pk_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_z_yyyyyyy = cbuffer.data(pk_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_z_yyyyyyz = cbuffer.data(pk_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_z_yyyyyzz = cbuffer.data(pk_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_z_yyyyzzz = cbuffer.data(pk_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_z_yyyzzzz = cbuffer.data(pk_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_z_yyzzzzz = cbuffer.data(pk_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_z_yzzzzzz = cbuffer.data(pk_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_z_zzzzzzz = cbuffer.data(pk_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_x_xxxxxxx = cbuffer.data(pk_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_x_xxxxxxy = cbuffer.data(pk_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_x_xxxxxxz = cbuffer.data(pk_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_x_xxxxxyy = cbuffer.data(pk_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_x_xxxxxyz = cbuffer.data(pk_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_x_xxxxxzz = cbuffer.data(pk_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_x_xxxxyyy = cbuffer.data(pk_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_x_xxxxyyz = cbuffer.data(pk_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_x_xxxxyzz = cbuffer.data(pk_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_x_xxxxzzz = cbuffer.data(pk_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_x_xxxyyyy = cbuffer.data(pk_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_x_xxxyyyz = cbuffer.data(pk_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_x_xxxyyzz = cbuffer.data(pk_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_x_xxxyzzz = cbuffer.data(pk_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_x_xxxzzzz = cbuffer.data(pk_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_x_xxyyyyy = cbuffer.data(pk_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_x_xxyyyyz = cbuffer.data(pk_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_x_xxyyyzz = cbuffer.data(pk_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_x_xxyyzzz = cbuffer.data(pk_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_x_xxyzzzz = cbuffer.data(pk_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_x_xxzzzzz = cbuffer.data(pk_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_x_xyyyyyy = cbuffer.data(pk_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_x_xyyyyyz = cbuffer.data(pk_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_x_xyyyyzz = cbuffer.data(pk_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_x_xyyyzzz = cbuffer.data(pk_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_x_xyyzzzz = cbuffer.data(pk_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_x_xyzzzzz = cbuffer.data(pk_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_x_xzzzzzz = cbuffer.data(pk_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_x_yyyyyyy = cbuffer.data(pk_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_x_yyyyyyz = cbuffer.data(pk_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_x_yyyyyzz = cbuffer.data(pk_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_x_yyyyzzz = cbuffer.data(pk_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_x_yyyzzzz = cbuffer.data(pk_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_x_yyzzzzz = cbuffer.data(pk_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_x_yzzzzzz = cbuffer.data(pk_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_x_zzzzzzz = cbuffer.data(pk_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_y_xxxxxxx = cbuffer.data(pk_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_y_xxxxxxy = cbuffer.data(pk_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_y_xxxxxxz = cbuffer.data(pk_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_y_xxxxxyy = cbuffer.data(pk_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_y_xxxxxyz = cbuffer.data(pk_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_y_xxxxxzz = cbuffer.data(pk_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_y_xxxxyyy = cbuffer.data(pk_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_y_xxxxyyz = cbuffer.data(pk_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_y_xxxxyzz = cbuffer.data(pk_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_y_xxxxzzz = cbuffer.data(pk_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_y_xxxyyyy = cbuffer.data(pk_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_y_xxxyyyz = cbuffer.data(pk_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_y_xxxyyzz = cbuffer.data(pk_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_y_xxxyzzz = cbuffer.data(pk_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_y_xxxzzzz = cbuffer.data(pk_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_y_xxyyyyy = cbuffer.data(pk_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_y_xxyyyyz = cbuffer.data(pk_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_y_xxyyyzz = cbuffer.data(pk_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_y_xxyyzzz = cbuffer.data(pk_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_y_xxyzzzz = cbuffer.data(pk_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_y_xxzzzzz = cbuffer.data(pk_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_y_xyyyyyy = cbuffer.data(pk_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_y_xyyyyyz = cbuffer.data(pk_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_y_xyyyyzz = cbuffer.data(pk_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_y_xyyyzzz = cbuffer.data(pk_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_y_xyyzzzz = cbuffer.data(pk_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_y_xyzzzzz = cbuffer.data(pk_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_y_xzzzzzz = cbuffer.data(pk_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_y_yyyyyyy = cbuffer.data(pk_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_y_yyyyyyz = cbuffer.data(pk_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_y_yyyyyzz = cbuffer.data(pk_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_y_yyyyzzz = cbuffer.data(pk_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_y_yyyzzzz = cbuffer.data(pk_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_y_yyzzzzz = cbuffer.data(pk_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_y_yzzzzzz = cbuffer.data(pk_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_y_zzzzzzz = cbuffer.data(pk_geom_20_off + 179 * ccomps * dcomps);

            auto g_xy_0_z_xxxxxxx = cbuffer.data(pk_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_z_xxxxxxy = cbuffer.data(pk_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_z_xxxxxxz = cbuffer.data(pk_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_z_xxxxxyy = cbuffer.data(pk_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_z_xxxxxyz = cbuffer.data(pk_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_z_xxxxxzz = cbuffer.data(pk_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_z_xxxxyyy = cbuffer.data(pk_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_z_xxxxyyz = cbuffer.data(pk_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_z_xxxxyzz = cbuffer.data(pk_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_z_xxxxzzz = cbuffer.data(pk_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_z_xxxyyyy = cbuffer.data(pk_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_z_xxxyyyz = cbuffer.data(pk_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_z_xxxyyzz = cbuffer.data(pk_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_z_xxxyzzz = cbuffer.data(pk_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_z_xxxzzzz = cbuffer.data(pk_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_z_xxyyyyy = cbuffer.data(pk_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_z_xxyyyyz = cbuffer.data(pk_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_z_xxyyyzz = cbuffer.data(pk_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_z_xxyyzzz = cbuffer.data(pk_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_z_xxyzzzz = cbuffer.data(pk_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_z_xxzzzzz = cbuffer.data(pk_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_z_xyyyyyy = cbuffer.data(pk_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_z_xyyyyyz = cbuffer.data(pk_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_z_xyyyyzz = cbuffer.data(pk_geom_20_off + 203 * ccomps * dcomps);

            auto g_xy_0_z_xyyyzzz = cbuffer.data(pk_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_z_xyyzzzz = cbuffer.data(pk_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_z_xyzzzzz = cbuffer.data(pk_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_z_xzzzzzz = cbuffer.data(pk_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_z_yyyyyyy = cbuffer.data(pk_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_z_yyyyyyz = cbuffer.data(pk_geom_20_off + 209 * ccomps * dcomps);

            auto g_xy_0_z_yyyyyzz = cbuffer.data(pk_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_z_yyyyzzz = cbuffer.data(pk_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_z_yyyzzzz = cbuffer.data(pk_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_z_yyzzzzz = cbuffer.data(pk_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_z_yzzzzzz = cbuffer.data(pk_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_z_zzzzzzz = cbuffer.data(pk_geom_20_off + 215 * ccomps * dcomps);

            auto g_xz_0_x_xxxxxxx = cbuffer.data(pk_geom_20_off + 216 * ccomps * dcomps);

            auto g_xz_0_x_xxxxxxy = cbuffer.data(pk_geom_20_off + 217 * ccomps * dcomps);

            auto g_xz_0_x_xxxxxxz = cbuffer.data(pk_geom_20_off + 218 * ccomps * dcomps);

            auto g_xz_0_x_xxxxxyy = cbuffer.data(pk_geom_20_off + 219 * ccomps * dcomps);

            auto g_xz_0_x_xxxxxyz = cbuffer.data(pk_geom_20_off + 220 * ccomps * dcomps);

            auto g_xz_0_x_xxxxxzz = cbuffer.data(pk_geom_20_off + 221 * ccomps * dcomps);

            auto g_xz_0_x_xxxxyyy = cbuffer.data(pk_geom_20_off + 222 * ccomps * dcomps);

            auto g_xz_0_x_xxxxyyz = cbuffer.data(pk_geom_20_off + 223 * ccomps * dcomps);

            auto g_xz_0_x_xxxxyzz = cbuffer.data(pk_geom_20_off + 224 * ccomps * dcomps);

            auto g_xz_0_x_xxxxzzz = cbuffer.data(pk_geom_20_off + 225 * ccomps * dcomps);

            auto g_xz_0_x_xxxyyyy = cbuffer.data(pk_geom_20_off + 226 * ccomps * dcomps);

            auto g_xz_0_x_xxxyyyz = cbuffer.data(pk_geom_20_off + 227 * ccomps * dcomps);

            auto g_xz_0_x_xxxyyzz = cbuffer.data(pk_geom_20_off + 228 * ccomps * dcomps);

            auto g_xz_0_x_xxxyzzz = cbuffer.data(pk_geom_20_off + 229 * ccomps * dcomps);

            auto g_xz_0_x_xxxzzzz = cbuffer.data(pk_geom_20_off + 230 * ccomps * dcomps);

            auto g_xz_0_x_xxyyyyy = cbuffer.data(pk_geom_20_off + 231 * ccomps * dcomps);

            auto g_xz_0_x_xxyyyyz = cbuffer.data(pk_geom_20_off + 232 * ccomps * dcomps);

            auto g_xz_0_x_xxyyyzz = cbuffer.data(pk_geom_20_off + 233 * ccomps * dcomps);

            auto g_xz_0_x_xxyyzzz = cbuffer.data(pk_geom_20_off + 234 * ccomps * dcomps);

            auto g_xz_0_x_xxyzzzz = cbuffer.data(pk_geom_20_off + 235 * ccomps * dcomps);

            auto g_xz_0_x_xxzzzzz = cbuffer.data(pk_geom_20_off + 236 * ccomps * dcomps);

            auto g_xz_0_x_xyyyyyy = cbuffer.data(pk_geom_20_off + 237 * ccomps * dcomps);

            auto g_xz_0_x_xyyyyyz = cbuffer.data(pk_geom_20_off + 238 * ccomps * dcomps);

            auto g_xz_0_x_xyyyyzz = cbuffer.data(pk_geom_20_off + 239 * ccomps * dcomps);

            auto g_xz_0_x_xyyyzzz = cbuffer.data(pk_geom_20_off + 240 * ccomps * dcomps);

            auto g_xz_0_x_xyyzzzz = cbuffer.data(pk_geom_20_off + 241 * ccomps * dcomps);

            auto g_xz_0_x_xyzzzzz = cbuffer.data(pk_geom_20_off + 242 * ccomps * dcomps);

            auto g_xz_0_x_xzzzzzz = cbuffer.data(pk_geom_20_off + 243 * ccomps * dcomps);

            auto g_xz_0_x_yyyyyyy = cbuffer.data(pk_geom_20_off + 244 * ccomps * dcomps);

            auto g_xz_0_x_yyyyyyz = cbuffer.data(pk_geom_20_off + 245 * ccomps * dcomps);

            auto g_xz_0_x_yyyyyzz = cbuffer.data(pk_geom_20_off + 246 * ccomps * dcomps);

            auto g_xz_0_x_yyyyzzz = cbuffer.data(pk_geom_20_off + 247 * ccomps * dcomps);

            auto g_xz_0_x_yyyzzzz = cbuffer.data(pk_geom_20_off + 248 * ccomps * dcomps);

            auto g_xz_0_x_yyzzzzz = cbuffer.data(pk_geom_20_off + 249 * ccomps * dcomps);

            auto g_xz_0_x_yzzzzzz = cbuffer.data(pk_geom_20_off + 250 * ccomps * dcomps);

            auto g_xz_0_x_zzzzzzz = cbuffer.data(pk_geom_20_off + 251 * ccomps * dcomps);

            auto g_xz_0_y_xxxxxxx = cbuffer.data(pk_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_y_xxxxxxy = cbuffer.data(pk_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_y_xxxxxxz = cbuffer.data(pk_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_y_xxxxxyy = cbuffer.data(pk_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_y_xxxxxyz = cbuffer.data(pk_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_y_xxxxxzz = cbuffer.data(pk_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_y_xxxxyyy = cbuffer.data(pk_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_y_xxxxyyz = cbuffer.data(pk_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_y_xxxxyzz = cbuffer.data(pk_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_y_xxxxzzz = cbuffer.data(pk_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_y_xxxyyyy = cbuffer.data(pk_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_y_xxxyyyz = cbuffer.data(pk_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_y_xxxyyzz = cbuffer.data(pk_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_y_xxxyzzz = cbuffer.data(pk_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_y_xxxzzzz = cbuffer.data(pk_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_y_xxyyyyy = cbuffer.data(pk_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_y_xxyyyyz = cbuffer.data(pk_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_y_xxyyyzz = cbuffer.data(pk_geom_20_off + 269 * ccomps * dcomps);

            auto g_xz_0_y_xxyyzzz = cbuffer.data(pk_geom_20_off + 270 * ccomps * dcomps);

            auto g_xz_0_y_xxyzzzz = cbuffer.data(pk_geom_20_off + 271 * ccomps * dcomps);

            auto g_xz_0_y_xxzzzzz = cbuffer.data(pk_geom_20_off + 272 * ccomps * dcomps);

            auto g_xz_0_y_xyyyyyy = cbuffer.data(pk_geom_20_off + 273 * ccomps * dcomps);

            auto g_xz_0_y_xyyyyyz = cbuffer.data(pk_geom_20_off + 274 * ccomps * dcomps);

            auto g_xz_0_y_xyyyyzz = cbuffer.data(pk_geom_20_off + 275 * ccomps * dcomps);

            auto g_xz_0_y_xyyyzzz = cbuffer.data(pk_geom_20_off + 276 * ccomps * dcomps);

            auto g_xz_0_y_xyyzzzz = cbuffer.data(pk_geom_20_off + 277 * ccomps * dcomps);

            auto g_xz_0_y_xyzzzzz = cbuffer.data(pk_geom_20_off + 278 * ccomps * dcomps);

            auto g_xz_0_y_xzzzzzz = cbuffer.data(pk_geom_20_off + 279 * ccomps * dcomps);

            auto g_xz_0_y_yyyyyyy = cbuffer.data(pk_geom_20_off + 280 * ccomps * dcomps);

            auto g_xz_0_y_yyyyyyz = cbuffer.data(pk_geom_20_off + 281 * ccomps * dcomps);

            auto g_xz_0_y_yyyyyzz = cbuffer.data(pk_geom_20_off + 282 * ccomps * dcomps);

            auto g_xz_0_y_yyyyzzz = cbuffer.data(pk_geom_20_off + 283 * ccomps * dcomps);

            auto g_xz_0_y_yyyzzzz = cbuffer.data(pk_geom_20_off + 284 * ccomps * dcomps);

            auto g_xz_0_y_yyzzzzz = cbuffer.data(pk_geom_20_off + 285 * ccomps * dcomps);

            auto g_xz_0_y_yzzzzzz = cbuffer.data(pk_geom_20_off + 286 * ccomps * dcomps);

            auto g_xz_0_y_zzzzzzz = cbuffer.data(pk_geom_20_off + 287 * ccomps * dcomps);

            auto g_xz_0_z_xxxxxxx = cbuffer.data(pk_geom_20_off + 288 * ccomps * dcomps);

            auto g_xz_0_z_xxxxxxy = cbuffer.data(pk_geom_20_off + 289 * ccomps * dcomps);

            auto g_xz_0_z_xxxxxxz = cbuffer.data(pk_geom_20_off + 290 * ccomps * dcomps);

            auto g_xz_0_z_xxxxxyy = cbuffer.data(pk_geom_20_off + 291 * ccomps * dcomps);

            auto g_xz_0_z_xxxxxyz = cbuffer.data(pk_geom_20_off + 292 * ccomps * dcomps);

            auto g_xz_0_z_xxxxxzz = cbuffer.data(pk_geom_20_off + 293 * ccomps * dcomps);

            auto g_xz_0_z_xxxxyyy = cbuffer.data(pk_geom_20_off + 294 * ccomps * dcomps);

            auto g_xz_0_z_xxxxyyz = cbuffer.data(pk_geom_20_off + 295 * ccomps * dcomps);

            auto g_xz_0_z_xxxxyzz = cbuffer.data(pk_geom_20_off + 296 * ccomps * dcomps);

            auto g_xz_0_z_xxxxzzz = cbuffer.data(pk_geom_20_off + 297 * ccomps * dcomps);

            auto g_xz_0_z_xxxyyyy = cbuffer.data(pk_geom_20_off + 298 * ccomps * dcomps);

            auto g_xz_0_z_xxxyyyz = cbuffer.data(pk_geom_20_off + 299 * ccomps * dcomps);

            auto g_xz_0_z_xxxyyzz = cbuffer.data(pk_geom_20_off + 300 * ccomps * dcomps);

            auto g_xz_0_z_xxxyzzz = cbuffer.data(pk_geom_20_off + 301 * ccomps * dcomps);

            auto g_xz_0_z_xxxzzzz = cbuffer.data(pk_geom_20_off + 302 * ccomps * dcomps);

            auto g_xz_0_z_xxyyyyy = cbuffer.data(pk_geom_20_off + 303 * ccomps * dcomps);

            auto g_xz_0_z_xxyyyyz = cbuffer.data(pk_geom_20_off + 304 * ccomps * dcomps);

            auto g_xz_0_z_xxyyyzz = cbuffer.data(pk_geom_20_off + 305 * ccomps * dcomps);

            auto g_xz_0_z_xxyyzzz = cbuffer.data(pk_geom_20_off + 306 * ccomps * dcomps);

            auto g_xz_0_z_xxyzzzz = cbuffer.data(pk_geom_20_off + 307 * ccomps * dcomps);

            auto g_xz_0_z_xxzzzzz = cbuffer.data(pk_geom_20_off + 308 * ccomps * dcomps);

            auto g_xz_0_z_xyyyyyy = cbuffer.data(pk_geom_20_off + 309 * ccomps * dcomps);

            auto g_xz_0_z_xyyyyyz = cbuffer.data(pk_geom_20_off + 310 * ccomps * dcomps);

            auto g_xz_0_z_xyyyyzz = cbuffer.data(pk_geom_20_off + 311 * ccomps * dcomps);

            auto g_xz_0_z_xyyyzzz = cbuffer.data(pk_geom_20_off + 312 * ccomps * dcomps);

            auto g_xz_0_z_xyyzzzz = cbuffer.data(pk_geom_20_off + 313 * ccomps * dcomps);

            auto g_xz_0_z_xyzzzzz = cbuffer.data(pk_geom_20_off + 314 * ccomps * dcomps);

            auto g_xz_0_z_xzzzzzz = cbuffer.data(pk_geom_20_off + 315 * ccomps * dcomps);

            auto g_xz_0_z_yyyyyyy = cbuffer.data(pk_geom_20_off + 316 * ccomps * dcomps);

            auto g_xz_0_z_yyyyyyz = cbuffer.data(pk_geom_20_off + 317 * ccomps * dcomps);

            auto g_xz_0_z_yyyyyzz = cbuffer.data(pk_geom_20_off + 318 * ccomps * dcomps);

            auto g_xz_0_z_yyyyzzz = cbuffer.data(pk_geom_20_off + 319 * ccomps * dcomps);

            auto g_xz_0_z_yyyzzzz = cbuffer.data(pk_geom_20_off + 320 * ccomps * dcomps);

            auto g_xz_0_z_yyzzzzz = cbuffer.data(pk_geom_20_off + 321 * ccomps * dcomps);

            auto g_xz_0_z_yzzzzzz = cbuffer.data(pk_geom_20_off + 322 * ccomps * dcomps);

            auto g_xz_0_z_zzzzzzz = cbuffer.data(pk_geom_20_off + 323 * ccomps * dcomps);

            auto g_yy_0_x_xxxxxxx = cbuffer.data(pk_geom_20_off + 324 * ccomps * dcomps);

            auto g_yy_0_x_xxxxxxy = cbuffer.data(pk_geom_20_off + 325 * ccomps * dcomps);

            auto g_yy_0_x_xxxxxxz = cbuffer.data(pk_geom_20_off + 326 * ccomps * dcomps);

            auto g_yy_0_x_xxxxxyy = cbuffer.data(pk_geom_20_off + 327 * ccomps * dcomps);

            auto g_yy_0_x_xxxxxyz = cbuffer.data(pk_geom_20_off + 328 * ccomps * dcomps);

            auto g_yy_0_x_xxxxxzz = cbuffer.data(pk_geom_20_off + 329 * ccomps * dcomps);

            auto g_yy_0_x_xxxxyyy = cbuffer.data(pk_geom_20_off + 330 * ccomps * dcomps);

            auto g_yy_0_x_xxxxyyz = cbuffer.data(pk_geom_20_off + 331 * ccomps * dcomps);

            auto g_yy_0_x_xxxxyzz = cbuffer.data(pk_geom_20_off + 332 * ccomps * dcomps);

            auto g_yy_0_x_xxxxzzz = cbuffer.data(pk_geom_20_off + 333 * ccomps * dcomps);

            auto g_yy_0_x_xxxyyyy = cbuffer.data(pk_geom_20_off + 334 * ccomps * dcomps);

            auto g_yy_0_x_xxxyyyz = cbuffer.data(pk_geom_20_off + 335 * ccomps * dcomps);

            auto g_yy_0_x_xxxyyzz = cbuffer.data(pk_geom_20_off + 336 * ccomps * dcomps);

            auto g_yy_0_x_xxxyzzz = cbuffer.data(pk_geom_20_off + 337 * ccomps * dcomps);

            auto g_yy_0_x_xxxzzzz = cbuffer.data(pk_geom_20_off + 338 * ccomps * dcomps);

            auto g_yy_0_x_xxyyyyy = cbuffer.data(pk_geom_20_off + 339 * ccomps * dcomps);

            auto g_yy_0_x_xxyyyyz = cbuffer.data(pk_geom_20_off + 340 * ccomps * dcomps);

            auto g_yy_0_x_xxyyyzz = cbuffer.data(pk_geom_20_off + 341 * ccomps * dcomps);

            auto g_yy_0_x_xxyyzzz = cbuffer.data(pk_geom_20_off + 342 * ccomps * dcomps);

            auto g_yy_0_x_xxyzzzz = cbuffer.data(pk_geom_20_off + 343 * ccomps * dcomps);

            auto g_yy_0_x_xxzzzzz = cbuffer.data(pk_geom_20_off + 344 * ccomps * dcomps);

            auto g_yy_0_x_xyyyyyy = cbuffer.data(pk_geom_20_off + 345 * ccomps * dcomps);

            auto g_yy_0_x_xyyyyyz = cbuffer.data(pk_geom_20_off + 346 * ccomps * dcomps);

            auto g_yy_0_x_xyyyyzz = cbuffer.data(pk_geom_20_off + 347 * ccomps * dcomps);

            auto g_yy_0_x_xyyyzzz = cbuffer.data(pk_geom_20_off + 348 * ccomps * dcomps);

            auto g_yy_0_x_xyyzzzz = cbuffer.data(pk_geom_20_off + 349 * ccomps * dcomps);

            auto g_yy_0_x_xyzzzzz = cbuffer.data(pk_geom_20_off + 350 * ccomps * dcomps);

            auto g_yy_0_x_xzzzzzz = cbuffer.data(pk_geom_20_off + 351 * ccomps * dcomps);

            auto g_yy_0_x_yyyyyyy = cbuffer.data(pk_geom_20_off + 352 * ccomps * dcomps);

            auto g_yy_0_x_yyyyyyz = cbuffer.data(pk_geom_20_off + 353 * ccomps * dcomps);

            auto g_yy_0_x_yyyyyzz = cbuffer.data(pk_geom_20_off + 354 * ccomps * dcomps);

            auto g_yy_0_x_yyyyzzz = cbuffer.data(pk_geom_20_off + 355 * ccomps * dcomps);

            auto g_yy_0_x_yyyzzzz = cbuffer.data(pk_geom_20_off + 356 * ccomps * dcomps);

            auto g_yy_0_x_yyzzzzz = cbuffer.data(pk_geom_20_off + 357 * ccomps * dcomps);

            auto g_yy_0_x_yzzzzzz = cbuffer.data(pk_geom_20_off + 358 * ccomps * dcomps);

            auto g_yy_0_x_zzzzzzz = cbuffer.data(pk_geom_20_off + 359 * ccomps * dcomps);

            auto g_yy_0_y_xxxxxxx = cbuffer.data(pk_geom_20_off + 360 * ccomps * dcomps);

            auto g_yy_0_y_xxxxxxy = cbuffer.data(pk_geom_20_off + 361 * ccomps * dcomps);

            auto g_yy_0_y_xxxxxxz = cbuffer.data(pk_geom_20_off + 362 * ccomps * dcomps);

            auto g_yy_0_y_xxxxxyy = cbuffer.data(pk_geom_20_off + 363 * ccomps * dcomps);

            auto g_yy_0_y_xxxxxyz = cbuffer.data(pk_geom_20_off + 364 * ccomps * dcomps);

            auto g_yy_0_y_xxxxxzz = cbuffer.data(pk_geom_20_off + 365 * ccomps * dcomps);

            auto g_yy_0_y_xxxxyyy = cbuffer.data(pk_geom_20_off + 366 * ccomps * dcomps);

            auto g_yy_0_y_xxxxyyz = cbuffer.data(pk_geom_20_off + 367 * ccomps * dcomps);

            auto g_yy_0_y_xxxxyzz = cbuffer.data(pk_geom_20_off + 368 * ccomps * dcomps);

            auto g_yy_0_y_xxxxzzz = cbuffer.data(pk_geom_20_off + 369 * ccomps * dcomps);

            auto g_yy_0_y_xxxyyyy = cbuffer.data(pk_geom_20_off + 370 * ccomps * dcomps);

            auto g_yy_0_y_xxxyyyz = cbuffer.data(pk_geom_20_off + 371 * ccomps * dcomps);

            auto g_yy_0_y_xxxyyzz = cbuffer.data(pk_geom_20_off + 372 * ccomps * dcomps);

            auto g_yy_0_y_xxxyzzz = cbuffer.data(pk_geom_20_off + 373 * ccomps * dcomps);

            auto g_yy_0_y_xxxzzzz = cbuffer.data(pk_geom_20_off + 374 * ccomps * dcomps);

            auto g_yy_0_y_xxyyyyy = cbuffer.data(pk_geom_20_off + 375 * ccomps * dcomps);

            auto g_yy_0_y_xxyyyyz = cbuffer.data(pk_geom_20_off + 376 * ccomps * dcomps);

            auto g_yy_0_y_xxyyyzz = cbuffer.data(pk_geom_20_off + 377 * ccomps * dcomps);

            auto g_yy_0_y_xxyyzzz = cbuffer.data(pk_geom_20_off + 378 * ccomps * dcomps);

            auto g_yy_0_y_xxyzzzz = cbuffer.data(pk_geom_20_off + 379 * ccomps * dcomps);

            auto g_yy_0_y_xxzzzzz = cbuffer.data(pk_geom_20_off + 380 * ccomps * dcomps);

            auto g_yy_0_y_xyyyyyy = cbuffer.data(pk_geom_20_off + 381 * ccomps * dcomps);

            auto g_yy_0_y_xyyyyyz = cbuffer.data(pk_geom_20_off + 382 * ccomps * dcomps);

            auto g_yy_0_y_xyyyyzz = cbuffer.data(pk_geom_20_off + 383 * ccomps * dcomps);

            auto g_yy_0_y_xyyyzzz = cbuffer.data(pk_geom_20_off + 384 * ccomps * dcomps);

            auto g_yy_0_y_xyyzzzz = cbuffer.data(pk_geom_20_off + 385 * ccomps * dcomps);

            auto g_yy_0_y_xyzzzzz = cbuffer.data(pk_geom_20_off + 386 * ccomps * dcomps);

            auto g_yy_0_y_xzzzzzz = cbuffer.data(pk_geom_20_off + 387 * ccomps * dcomps);

            auto g_yy_0_y_yyyyyyy = cbuffer.data(pk_geom_20_off + 388 * ccomps * dcomps);

            auto g_yy_0_y_yyyyyyz = cbuffer.data(pk_geom_20_off + 389 * ccomps * dcomps);

            auto g_yy_0_y_yyyyyzz = cbuffer.data(pk_geom_20_off + 390 * ccomps * dcomps);

            auto g_yy_0_y_yyyyzzz = cbuffer.data(pk_geom_20_off + 391 * ccomps * dcomps);

            auto g_yy_0_y_yyyzzzz = cbuffer.data(pk_geom_20_off + 392 * ccomps * dcomps);

            auto g_yy_0_y_yyzzzzz = cbuffer.data(pk_geom_20_off + 393 * ccomps * dcomps);

            auto g_yy_0_y_yzzzzzz = cbuffer.data(pk_geom_20_off + 394 * ccomps * dcomps);

            auto g_yy_0_y_zzzzzzz = cbuffer.data(pk_geom_20_off + 395 * ccomps * dcomps);

            auto g_yy_0_z_xxxxxxx = cbuffer.data(pk_geom_20_off + 396 * ccomps * dcomps);

            auto g_yy_0_z_xxxxxxy = cbuffer.data(pk_geom_20_off + 397 * ccomps * dcomps);

            auto g_yy_0_z_xxxxxxz = cbuffer.data(pk_geom_20_off + 398 * ccomps * dcomps);

            auto g_yy_0_z_xxxxxyy = cbuffer.data(pk_geom_20_off + 399 * ccomps * dcomps);

            auto g_yy_0_z_xxxxxyz = cbuffer.data(pk_geom_20_off + 400 * ccomps * dcomps);

            auto g_yy_0_z_xxxxxzz = cbuffer.data(pk_geom_20_off + 401 * ccomps * dcomps);

            auto g_yy_0_z_xxxxyyy = cbuffer.data(pk_geom_20_off + 402 * ccomps * dcomps);

            auto g_yy_0_z_xxxxyyz = cbuffer.data(pk_geom_20_off + 403 * ccomps * dcomps);

            auto g_yy_0_z_xxxxyzz = cbuffer.data(pk_geom_20_off + 404 * ccomps * dcomps);

            auto g_yy_0_z_xxxxzzz = cbuffer.data(pk_geom_20_off + 405 * ccomps * dcomps);

            auto g_yy_0_z_xxxyyyy = cbuffer.data(pk_geom_20_off + 406 * ccomps * dcomps);

            auto g_yy_0_z_xxxyyyz = cbuffer.data(pk_geom_20_off + 407 * ccomps * dcomps);

            auto g_yy_0_z_xxxyyzz = cbuffer.data(pk_geom_20_off + 408 * ccomps * dcomps);

            auto g_yy_0_z_xxxyzzz = cbuffer.data(pk_geom_20_off + 409 * ccomps * dcomps);

            auto g_yy_0_z_xxxzzzz = cbuffer.data(pk_geom_20_off + 410 * ccomps * dcomps);

            auto g_yy_0_z_xxyyyyy = cbuffer.data(pk_geom_20_off + 411 * ccomps * dcomps);

            auto g_yy_0_z_xxyyyyz = cbuffer.data(pk_geom_20_off + 412 * ccomps * dcomps);

            auto g_yy_0_z_xxyyyzz = cbuffer.data(pk_geom_20_off + 413 * ccomps * dcomps);

            auto g_yy_0_z_xxyyzzz = cbuffer.data(pk_geom_20_off + 414 * ccomps * dcomps);

            auto g_yy_0_z_xxyzzzz = cbuffer.data(pk_geom_20_off + 415 * ccomps * dcomps);

            auto g_yy_0_z_xxzzzzz = cbuffer.data(pk_geom_20_off + 416 * ccomps * dcomps);

            auto g_yy_0_z_xyyyyyy = cbuffer.data(pk_geom_20_off + 417 * ccomps * dcomps);

            auto g_yy_0_z_xyyyyyz = cbuffer.data(pk_geom_20_off + 418 * ccomps * dcomps);

            auto g_yy_0_z_xyyyyzz = cbuffer.data(pk_geom_20_off + 419 * ccomps * dcomps);

            auto g_yy_0_z_xyyyzzz = cbuffer.data(pk_geom_20_off + 420 * ccomps * dcomps);

            auto g_yy_0_z_xyyzzzz = cbuffer.data(pk_geom_20_off + 421 * ccomps * dcomps);

            auto g_yy_0_z_xyzzzzz = cbuffer.data(pk_geom_20_off + 422 * ccomps * dcomps);

            auto g_yy_0_z_xzzzzzz = cbuffer.data(pk_geom_20_off + 423 * ccomps * dcomps);

            auto g_yy_0_z_yyyyyyy = cbuffer.data(pk_geom_20_off + 424 * ccomps * dcomps);

            auto g_yy_0_z_yyyyyyz = cbuffer.data(pk_geom_20_off + 425 * ccomps * dcomps);

            auto g_yy_0_z_yyyyyzz = cbuffer.data(pk_geom_20_off + 426 * ccomps * dcomps);

            auto g_yy_0_z_yyyyzzz = cbuffer.data(pk_geom_20_off + 427 * ccomps * dcomps);

            auto g_yy_0_z_yyyzzzz = cbuffer.data(pk_geom_20_off + 428 * ccomps * dcomps);

            auto g_yy_0_z_yyzzzzz = cbuffer.data(pk_geom_20_off + 429 * ccomps * dcomps);

            auto g_yy_0_z_yzzzzzz = cbuffer.data(pk_geom_20_off + 430 * ccomps * dcomps);

            auto g_yy_0_z_zzzzzzz = cbuffer.data(pk_geom_20_off + 431 * ccomps * dcomps);

            auto g_yz_0_x_xxxxxxx = cbuffer.data(pk_geom_20_off + 432 * ccomps * dcomps);

            auto g_yz_0_x_xxxxxxy = cbuffer.data(pk_geom_20_off + 433 * ccomps * dcomps);

            auto g_yz_0_x_xxxxxxz = cbuffer.data(pk_geom_20_off + 434 * ccomps * dcomps);

            auto g_yz_0_x_xxxxxyy = cbuffer.data(pk_geom_20_off + 435 * ccomps * dcomps);

            auto g_yz_0_x_xxxxxyz = cbuffer.data(pk_geom_20_off + 436 * ccomps * dcomps);

            auto g_yz_0_x_xxxxxzz = cbuffer.data(pk_geom_20_off + 437 * ccomps * dcomps);

            auto g_yz_0_x_xxxxyyy = cbuffer.data(pk_geom_20_off + 438 * ccomps * dcomps);

            auto g_yz_0_x_xxxxyyz = cbuffer.data(pk_geom_20_off + 439 * ccomps * dcomps);

            auto g_yz_0_x_xxxxyzz = cbuffer.data(pk_geom_20_off + 440 * ccomps * dcomps);

            auto g_yz_0_x_xxxxzzz = cbuffer.data(pk_geom_20_off + 441 * ccomps * dcomps);

            auto g_yz_0_x_xxxyyyy = cbuffer.data(pk_geom_20_off + 442 * ccomps * dcomps);

            auto g_yz_0_x_xxxyyyz = cbuffer.data(pk_geom_20_off + 443 * ccomps * dcomps);

            auto g_yz_0_x_xxxyyzz = cbuffer.data(pk_geom_20_off + 444 * ccomps * dcomps);

            auto g_yz_0_x_xxxyzzz = cbuffer.data(pk_geom_20_off + 445 * ccomps * dcomps);

            auto g_yz_0_x_xxxzzzz = cbuffer.data(pk_geom_20_off + 446 * ccomps * dcomps);

            auto g_yz_0_x_xxyyyyy = cbuffer.data(pk_geom_20_off + 447 * ccomps * dcomps);

            auto g_yz_0_x_xxyyyyz = cbuffer.data(pk_geom_20_off + 448 * ccomps * dcomps);

            auto g_yz_0_x_xxyyyzz = cbuffer.data(pk_geom_20_off + 449 * ccomps * dcomps);

            auto g_yz_0_x_xxyyzzz = cbuffer.data(pk_geom_20_off + 450 * ccomps * dcomps);

            auto g_yz_0_x_xxyzzzz = cbuffer.data(pk_geom_20_off + 451 * ccomps * dcomps);

            auto g_yz_0_x_xxzzzzz = cbuffer.data(pk_geom_20_off + 452 * ccomps * dcomps);

            auto g_yz_0_x_xyyyyyy = cbuffer.data(pk_geom_20_off + 453 * ccomps * dcomps);

            auto g_yz_0_x_xyyyyyz = cbuffer.data(pk_geom_20_off + 454 * ccomps * dcomps);

            auto g_yz_0_x_xyyyyzz = cbuffer.data(pk_geom_20_off + 455 * ccomps * dcomps);

            auto g_yz_0_x_xyyyzzz = cbuffer.data(pk_geom_20_off + 456 * ccomps * dcomps);

            auto g_yz_0_x_xyyzzzz = cbuffer.data(pk_geom_20_off + 457 * ccomps * dcomps);

            auto g_yz_0_x_xyzzzzz = cbuffer.data(pk_geom_20_off + 458 * ccomps * dcomps);

            auto g_yz_0_x_xzzzzzz = cbuffer.data(pk_geom_20_off + 459 * ccomps * dcomps);

            auto g_yz_0_x_yyyyyyy = cbuffer.data(pk_geom_20_off + 460 * ccomps * dcomps);

            auto g_yz_0_x_yyyyyyz = cbuffer.data(pk_geom_20_off + 461 * ccomps * dcomps);

            auto g_yz_0_x_yyyyyzz = cbuffer.data(pk_geom_20_off + 462 * ccomps * dcomps);

            auto g_yz_0_x_yyyyzzz = cbuffer.data(pk_geom_20_off + 463 * ccomps * dcomps);

            auto g_yz_0_x_yyyzzzz = cbuffer.data(pk_geom_20_off + 464 * ccomps * dcomps);

            auto g_yz_0_x_yyzzzzz = cbuffer.data(pk_geom_20_off + 465 * ccomps * dcomps);

            auto g_yz_0_x_yzzzzzz = cbuffer.data(pk_geom_20_off + 466 * ccomps * dcomps);

            auto g_yz_0_x_zzzzzzz = cbuffer.data(pk_geom_20_off + 467 * ccomps * dcomps);

            auto g_yz_0_y_xxxxxxx = cbuffer.data(pk_geom_20_off + 468 * ccomps * dcomps);

            auto g_yz_0_y_xxxxxxy = cbuffer.data(pk_geom_20_off + 469 * ccomps * dcomps);

            auto g_yz_0_y_xxxxxxz = cbuffer.data(pk_geom_20_off + 470 * ccomps * dcomps);

            auto g_yz_0_y_xxxxxyy = cbuffer.data(pk_geom_20_off + 471 * ccomps * dcomps);

            auto g_yz_0_y_xxxxxyz = cbuffer.data(pk_geom_20_off + 472 * ccomps * dcomps);

            auto g_yz_0_y_xxxxxzz = cbuffer.data(pk_geom_20_off + 473 * ccomps * dcomps);

            auto g_yz_0_y_xxxxyyy = cbuffer.data(pk_geom_20_off + 474 * ccomps * dcomps);

            auto g_yz_0_y_xxxxyyz = cbuffer.data(pk_geom_20_off + 475 * ccomps * dcomps);

            auto g_yz_0_y_xxxxyzz = cbuffer.data(pk_geom_20_off + 476 * ccomps * dcomps);

            auto g_yz_0_y_xxxxzzz = cbuffer.data(pk_geom_20_off + 477 * ccomps * dcomps);

            auto g_yz_0_y_xxxyyyy = cbuffer.data(pk_geom_20_off + 478 * ccomps * dcomps);

            auto g_yz_0_y_xxxyyyz = cbuffer.data(pk_geom_20_off + 479 * ccomps * dcomps);

            auto g_yz_0_y_xxxyyzz = cbuffer.data(pk_geom_20_off + 480 * ccomps * dcomps);

            auto g_yz_0_y_xxxyzzz = cbuffer.data(pk_geom_20_off + 481 * ccomps * dcomps);

            auto g_yz_0_y_xxxzzzz = cbuffer.data(pk_geom_20_off + 482 * ccomps * dcomps);

            auto g_yz_0_y_xxyyyyy = cbuffer.data(pk_geom_20_off + 483 * ccomps * dcomps);

            auto g_yz_0_y_xxyyyyz = cbuffer.data(pk_geom_20_off + 484 * ccomps * dcomps);

            auto g_yz_0_y_xxyyyzz = cbuffer.data(pk_geom_20_off + 485 * ccomps * dcomps);

            auto g_yz_0_y_xxyyzzz = cbuffer.data(pk_geom_20_off + 486 * ccomps * dcomps);

            auto g_yz_0_y_xxyzzzz = cbuffer.data(pk_geom_20_off + 487 * ccomps * dcomps);

            auto g_yz_0_y_xxzzzzz = cbuffer.data(pk_geom_20_off + 488 * ccomps * dcomps);

            auto g_yz_0_y_xyyyyyy = cbuffer.data(pk_geom_20_off + 489 * ccomps * dcomps);

            auto g_yz_0_y_xyyyyyz = cbuffer.data(pk_geom_20_off + 490 * ccomps * dcomps);

            auto g_yz_0_y_xyyyyzz = cbuffer.data(pk_geom_20_off + 491 * ccomps * dcomps);

            auto g_yz_0_y_xyyyzzz = cbuffer.data(pk_geom_20_off + 492 * ccomps * dcomps);

            auto g_yz_0_y_xyyzzzz = cbuffer.data(pk_geom_20_off + 493 * ccomps * dcomps);

            auto g_yz_0_y_xyzzzzz = cbuffer.data(pk_geom_20_off + 494 * ccomps * dcomps);

            auto g_yz_0_y_xzzzzzz = cbuffer.data(pk_geom_20_off + 495 * ccomps * dcomps);

            auto g_yz_0_y_yyyyyyy = cbuffer.data(pk_geom_20_off + 496 * ccomps * dcomps);

            auto g_yz_0_y_yyyyyyz = cbuffer.data(pk_geom_20_off + 497 * ccomps * dcomps);

            auto g_yz_0_y_yyyyyzz = cbuffer.data(pk_geom_20_off + 498 * ccomps * dcomps);

            auto g_yz_0_y_yyyyzzz = cbuffer.data(pk_geom_20_off + 499 * ccomps * dcomps);

            auto g_yz_0_y_yyyzzzz = cbuffer.data(pk_geom_20_off + 500 * ccomps * dcomps);

            auto g_yz_0_y_yyzzzzz = cbuffer.data(pk_geom_20_off + 501 * ccomps * dcomps);

            auto g_yz_0_y_yzzzzzz = cbuffer.data(pk_geom_20_off + 502 * ccomps * dcomps);

            auto g_yz_0_y_zzzzzzz = cbuffer.data(pk_geom_20_off + 503 * ccomps * dcomps);

            auto g_yz_0_z_xxxxxxx = cbuffer.data(pk_geom_20_off + 504 * ccomps * dcomps);

            auto g_yz_0_z_xxxxxxy = cbuffer.data(pk_geom_20_off + 505 * ccomps * dcomps);

            auto g_yz_0_z_xxxxxxz = cbuffer.data(pk_geom_20_off + 506 * ccomps * dcomps);

            auto g_yz_0_z_xxxxxyy = cbuffer.data(pk_geom_20_off + 507 * ccomps * dcomps);

            auto g_yz_0_z_xxxxxyz = cbuffer.data(pk_geom_20_off + 508 * ccomps * dcomps);

            auto g_yz_0_z_xxxxxzz = cbuffer.data(pk_geom_20_off + 509 * ccomps * dcomps);

            auto g_yz_0_z_xxxxyyy = cbuffer.data(pk_geom_20_off + 510 * ccomps * dcomps);

            auto g_yz_0_z_xxxxyyz = cbuffer.data(pk_geom_20_off + 511 * ccomps * dcomps);

            auto g_yz_0_z_xxxxyzz = cbuffer.data(pk_geom_20_off + 512 * ccomps * dcomps);

            auto g_yz_0_z_xxxxzzz = cbuffer.data(pk_geom_20_off + 513 * ccomps * dcomps);

            auto g_yz_0_z_xxxyyyy = cbuffer.data(pk_geom_20_off + 514 * ccomps * dcomps);

            auto g_yz_0_z_xxxyyyz = cbuffer.data(pk_geom_20_off + 515 * ccomps * dcomps);

            auto g_yz_0_z_xxxyyzz = cbuffer.data(pk_geom_20_off + 516 * ccomps * dcomps);

            auto g_yz_0_z_xxxyzzz = cbuffer.data(pk_geom_20_off + 517 * ccomps * dcomps);

            auto g_yz_0_z_xxxzzzz = cbuffer.data(pk_geom_20_off + 518 * ccomps * dcomps);

            auto g_yz_0_z_xxyyyyy = cbuffer.data(pk_geom_20_off + 519 * ccomps * dcomps);

            auto g_yz_0_z_xxyyyyz = cbuffer.data(pk_geom_20_off + 520 * ccomps * dcomps);

            auto g_yz_0_z_xxyyyzz = cbuffer.data(pk_geom_20_off + 521 * ccomps * dcomps);

            auto g_yz_0_z_xxyyzzz = cbuffer.data(pk_geom_20_off + 522 * ccomps * dcomps);

            auto g_yz_0_z_xxyzzzz = cbuffer.data(pk_geom_20_off + 523 * ccomps * dcomps);

            auto g_yz_0_z_xxzzzzz = cbuffer.data(pk_geom_20_off + 524 * ccomps * dcomps);

            auto g_yz_0_z_xyyyyyy = cbuffer.data(pk_geom_20_off + 525 * ccomps * dcomps);

            auto g_yz_0_z_xyyyyyz = cbuffer.data(pk_geom_20_off + 526 * ccomps * dcomps);

            auto g_yz_0_z_xyyyyzz = cbuffer.data(pk_geom_20_off + 527 * ccomps * dcomps);

            auto g_yz_0_z_xyyyzzz = cbuffer.data(pk_geom_20_off + 528 * ccomps * dcomps);

            auto g_yz_0_z_xyyzzzz = cbuffer.data(pk_geom_20_off + 529 * ccomps * dcomps);

            auto g_yz_0_z_xyzzzzz = cbuffer.data(pk_geom_20_off + 530 * ccomps * dcomps);

            auto g_yz_0_z_xzzzzzz = cbuffer.data(pk_geom_20_off + 531 * ccomps * dcomps);

            auto g_yz_0_z_yyyyyyy = cbuffer.data(pk_geom_20_off + 532 * ccomps * dcomps);

            auto g_yz_0_z_yyyyyyz = cbuffer.data(pk_geom_20_off + 533 * ccomps * dcomps);

            auto g_yz_0_z_yyyyyzz = cbuffer.data(pk_geom_20_off + 534 * ccomps * dcomps);

            auto g_yz_0_z_yyyyzzz = cbuffer.data(pk_geom_20_off + 535 * ccomps * dcomps);

            auto g_yz_0_z_yyyzzzz = cbuffer.data(pk_geom_20_off + 536 * ccomps * dcomps);

            auto g_yz_0_z_yyzzzzz = cbuffer.data(pk_geom_20_off + 537 * ccomps * dcomps);

            auto g_yz_0_z_yzzzzzz = cbuffer.data(pk_geom_20_off + 538 * ccomps * dcomps);

            auto g_yz_0_z_zzzzzzz = cbuffer.data(pk_geom_20_off + 539 * ccomps * dcomps);

            auto g_zz_0_x_xxxxxxx = cbuffer.data(pk_geom_20_off + 540 * ccomps * dcomps);

            auto g_zz_0_x_xxxxxxy = cbuffer.data(pk_geom_20_off + 541 * ccomps * dcomps);

            auto g_zz_0_x_xxxxxxz = cbuffer.data(pk_geom_20_off + 542 * ccomps * dcomps);

            auto g_zz_0_x_xxxxxyy = cbuffer.data(pk_geom_20_off + 543 * ccomps * dcomps);

            auto g_zz_0_x_xxxxxyz = cbuffer.data(pk_geom_20_off + 544 * ccomps * dcomps);

            auto g_zz_0_x_xxxxxzz = cbuffer.data(pk_geom_20_off + 545 * ccomps * dcomps);

            auto g_zz_0_x_xxxxyyy = cbuffer.data(pk_geom_20_off + 546 * ccomps * dcomps);

            auto g_zz_0_x_xxxxyyz = cbuffer.data(pk_geom_20_off + 547 * ccomps * dcomps);

            auto g_zz_0_x_xxxxyzz = cbuffer.data(pk_geom_20_off + 548 * ccomps * dcomps);

            auto g_zz_0_x_xxxxzzz = cbuffer.data(pk_geom_20_off + 549 * ccomps * dcomps);

            auto g_zz_0_x_xxxyyyy = cbuffer.data(pk_geom_20_off + 550 * ccomps * dcomps);

            auto g_zz_0_x_xxxyyyz = cbuffer.data(pk_geom_20_off + 551 * ccomps * dcomps);

            auto g_zz_0_x_xxxyyzz = cbuffer.data(pk_geom_20_off + 552 * ccomps * dcomps);

            auto g_zz_0_x_xxxyzzz = cbuffer.data(pk_geom_20_off + 553 * ccomps * dcomps);

            auto g_zz_0_x_xxxzzzz = cbuffer.data(pk_geom_20_off + 554 * ccomps * dcomps);

            auto g_zz_0_x_xxyyyyy = cbuffer.data(pk_geom_20_off + 555 * ccomps * dcomps);

            auto g_zz_0_x_xxyyyyz = cbuffer.data(pk_geom_20_off + 556 * ccomps * dcomps);

            auto g_zz_0_x_xxyyyzz = cbuffer.data(pk_geom_20_off + 557 * ccomps * dcomps);

            auto g_zz_0_x_xxyyzzz = cbuffer.data(pk_geom_20_off + 558 * ccomps * dcomps);

            auto g_zz_0_x_xxyzzzz = cbuffer.data(pk_geom_20_off + 559 * ccomps * dcomps);

            auto g_zz_0_x_xxzzzzz = cbuffer.data(pk_geom_20_off + 560 * ccomps * dcomps);

            auto g_zz_0_x_xyyyyyy = cbuffer.data(pk_geom_20_off + 561 * ccomps * dcomps);

            auto g_zz_0_x_xyyyyyz = cbuffer.data(pk_geom_20_off + 562 * ccomps * dcomps);

            auto g_zz_0_x_xyyyyzz = cbuffer.data(pk_geom_20_off + 563 * ccomps * dcomps);

            auto g_zz_0_x_xyyyzzz = cbuffer.data(pk_geom_20_off + 564 * ccomps * dcomps);

            auto g_zz_0_x_xyyzzzz = cbuffer.data(pk_geom_20_off + 565 * ccomps * dcomps);

            auto g_zz_0_x_xyzzzzz = cbuffer.data(pk_geom_20_off + 566 * ccomps * dcomps);

            auto g_zz_0_x_xzzzzzz = cbuffer.data(pk_geom_20_off + 567 * ccomps * dcomps);

            auto g_zz_0_x_yyyyyyy = cbuffer.data(pk_geom_20_off + 568 * ccomps * dcomps);

            auto g_zz_0_x_yyyyyyz = cbuffer.data(pk_geom_20_off + 569 * ccomps * dcomps);

            auto g_zz_0_x_yyyyyzz = cbuffer.data(pk_geom_20_off + 570 * ccomps * dcomps);

            auto g_zz_0_x_yyyyzzz = cbuffer.data(pk_geom_20_off + 571 * ccomps * dcomps);

            auto g_zz_0_x_yyyzzzz = cbuffer.data(pk_geom_20_off + 572 * ccomps * dcomps);

            auto g_zz_0_x_yyzzzzz = cbuffer.data(pk_geom_20_off + 573 * ccomps * dcomps);

            auto g_zz_0_x_yzzzzzz = cbuffer.data(pk_geom_20_off + 574 * ccomps * dcomps);

            auto g_zz_0_x_zzzzzzz = cbuffer.data(pk_geom_20_off + 575 * ccomps * dcomps);

            auto g_zz_0_y_xxxxxxx = cbuffer.data(pk_geom_20_off + 576 * ccomps * dcomps);

            auto g_zz_0_y_xxxxxxy = cbuffer.data(pk_geom_20_off + 577 * ccomps * dcomps);

            auto g_zz_0_y_xxxxxxz = cbuffer.data(pk_geom_20_off + 578 * ccomps * dcomps);

            auto g_zz_0_y_xxxxxyy = cbuffer.data(pk_geom_20_off + 579 * ccomps * dcomps);

            auto g_zz_0_y_xxxxxyz = cbuffer.data(pk_geom_20_off + 580 * ccomps * dcomps);

            auto g_zz_0_y_xxxxxzz = cbuffer.data(pk_geom_20_off + 581 * ccomps * dcomps);

            auto g_zz_0_y_xxxxyyy = cbuffer.data(pk_geom_20_off + 582 * ccomps * dcomps);

            auto g_zz_0_y_xxxxyyz = cbuffer.data(pk_geom_20_off + 583 * ccomps * dcomps);

            auto g_zz_0_y_xxxxyzz = cbuffer.data(pk_geom_20_off + 584 * ccomps * dcomps);

            auto g_zz_0_y_xxxxzzz = cbuffer.data(pk_geom_20_off + 585 * ccomps * dcomps);

            auto g_zz_0_y_xxxyyyy = cbuffer.data(pk_geom_20_off + 586 * ccomps * dcomps);

            auto g_zz_0_y_xxxyyyz = cbuffer.data(pk_geom_20_off + 587 * ccomps * dcomps);

            auto g_zz_0_y_xxxyyzz = cbuffer.data(pk_geom_20_off + 588 * ccomps * dcomps);

            auto g_zz_0_y_xxxyzzz = cbuffer.data(pk_geom_20_off + 589 * ccomps * dcomps);

            auto g_zz_0_y_xxxzzzz = cbuffer.data(pk_geom_20_off + 590 * ccomps * dcomps);

            auto g_zz_0_y_xxyyyyy = cbuffer.data(pk_geom_20_off + 591 * ccomps * dcomps);

            auto g_zz_0_y_xxyyyyz = cbuffer.data(pk_geom_20_off + 592 * ccomps * dcomps);

            auto g_zz_0_y_xxyyyzz = cbuffer.data(pk_geom_20_off + 593 * ccomps * dcomps);

            auto g_zz_0_y_xxyyzzz = cbuffer.data(pk_geom_20_off + 594 * ccomps * dcomps);

            auto g_zz_0_y_xxyzzzz = cbuffer.data(pk_geom_20_off + 595 * ccomps * dcomps);

            auto g_zz_0_y_xxzzzzz = cbuffer.data(pk_geom_20_off + 596 * ccomps * dcomps);

            auto g_zz_0_y_xyyyyyy = cbuffer.data(pk_geom_20_off + 597 * ccomps * dcomps);

            auto g_zz_0_y_xyyyyyz = cbuffer.data(pk_geom_20_off + 598 * ccomps * dcomps);

            auto g_zz_0_y_xyyyyzz = cbuffer.data(pk_geom_20_off + 599 * ccomps * dcomps);

            auto g_zz_0_y_xyyyzzz = cbuffer.data(pk_geom_20_off + 600 * ccomps * dcomps);

            auto g_zz_0_y_xyyzzzz = cbuffer.data(pk_geom_20_off + 601 * ccomps * dcomps);

            auto g_zz_0_y_xyzzzzz = cbuffer.data(pk_geom_20_off + 602 * ccomps * dcomps);

            auto g_zz_0_y_xzzzzzz = cbuffer.data(pk_geom_20_off + 603 * ccomps * dcomps);

            auto g_zz_0_y_yyyyyyy = cbuffer.data(pk_geom_20_off + 604 * ccomps * dcomps);

            auto g_zz_0_y_yyyyyyz = cbuffer.data(pk_geom_20_off + 605 * ccomps * dcomps);

            auto g_zz_0_y_yyyyyzz = cbuffer.data(pk_geom_20_off + 606 * ccomps * dcomps);

            auto g_zz_0_y_yyyyzzz = cbuffer.data(pk_geom_20_off + 607 * ccomps * dcomps);

            auto g_zz_0_y_yyyzzzz = cbuffer.data(pk_geom_20_off + 608 * ccomps * dcomps);

            auto g_zz_0_y_yyzzzzz = cbuffer.data(pk_geom_20_off + 609 * ccomps * dcomps);

            auto g_zz_0_y_yzzzzzz = cbuffer.data(pk_geom_20_off + 610 * ccomps * dcomps);

            auto g_zz_0_y_zzzzzzz = cbuffer.data(pk_geom_20_off + 611 * ccomps * dcomps);

            auto g_zz_0_z_xxxxxxx = cbuffer.data(pk_geom_20_off + 612 * ccomps * dcomps);

            auto g_zz_0_z_xxxxxxy = cbuffer.data(pk_geom_20_off + 613 * ccomps * dcomps);

            auto g_zz_0_z_xxxxxxz = cbuffer.data(pk_geom_20_off + 614 * ccomps * dcomps);

            auto g_zz_0_z_xxxxxyy = cbuffer.data(pk_geom_20_off + 615 * ccomps * dcomps);

            auto g_zz_0_z_xxxxxyz = cbuffer.data(pk_geom_20_off + 616 * ccomps * dcomps);

            auto g_zz_0_z_xxxxxzz = cbuffer.data(pk_geom_20_off + 617 * ccomps * dcomps);

            auto g_zz_0_z_xxxxyyy = cbuffer.data(pk_geom_20_off + 618 * ccomps * dcomps);

            auto g_zz_0_z_xxxxyyz = cbuffer.data(pk_geom_20_off + 619 * ccomps * dcomps);

            auto g_zz_0_z_xxxxyzz = cbuffer.data(pk_geom_20_off + 620 * ccomps * dcomps);

            auto g_zz_0_z_xxxxzzz = cbuffer.data(pk_geom_20_off + 621 * ccomps * dcomps);

            auto g_zz_0_z_xxxyyyy = cbuffer.data(pk_geom_20_off + 622 * ccomps * dcomps);

            auto g_zz_0_z_xxxyyyz = cbuffer.data(pk_geom_20_off + 623 * ccomps * dcomps);

            auto g_zz_0_z_xxxyyzz = cbuffer.data(pk_geom_20_off + 624 * ccomps * dcomps);

            auto g_zz_0_z_xxxyzzz = cbuffer.data(pk_geom_20_off + 625 * ccomps * dcomps);

            auto g_zz_0_z_xxxzzzz = cbuffer.data(pk_geom_20_off + 626 * ccomps * dcomps);

            auto g_zz_0_z_xxyyyyy = cbuffer.data(pk_geom_20_off + 627 * ccomps * dcomps);

            auto g_zz_0_z_xxyyyyz = cbuffer.data(pk_geom_20_off + 628 * ccomps * dcomps);

            auto g_zz_0_z_xxyyyzz = cbuffer.data(pk_geom_20_off + 629 * ccomps * dcomps);

            auto g_zz_0_z_xxyyzzz = cbuffer.data(pk_geom_20_off + 630 * ccomps * dcomps);

            auto g_zz_0_z_xxyzzzz = cbuffer.data(pk_geom_20_off + 631 * ccomps * dcomps);

            auto g_zz_0_z_xxzzzzz = cbuffer.data(pk_geom_20_off + 632 * ccomps * dcomps);

            auto g_zz_0_z_xyyyyyy = cbuffer.data(pk_geom_20_off + 633 * ccomps * dcomps);

            auto g_zz_0_z_xyyyyyz = cbuffer.data(pk_geom_20_off + 634 * ccomps * dcomps);

            auto g_zz_0_z_xyyyyzz = cbuffer.data(pk_geom_20_off + 635 * ccomps * dcomps);

            auto g_zz_0_z_xyyyzzz = cbuffer.data(pk_geom_20_off + 636 * ccomps * dcomps);

            auto g_zz_0_z_xyyzzzz = cbuffer.data(pk_geom_20_off + 637 * ccomps * dcomps);

            auto g_zz_0_z_xyzzzzz = cbuffer.data(pk_geom_20_off + 638 * ccomps * dcomps);

            auto g_zz_0_z_xzzzzzz = cbuffer.data(pk_geom_20_off + 639 * ccomps * dcomps);

            auto g_zz_0_z_yyyyyyy = cbuffer.data(pk_geom_20_off + 640 * ccomps * dcomps);

            auto g_zz_0_z_yyyyyyz = cbuffer.data(pk_geom_20_off + 641 * ccomps * dcomps);

            auto g_zz_0_z_yyyyyzz = cbuffer.data(pk_geom_20_off + 642 * ccomps * dcomps);

            auto g_zz_0_z_yyyyzzz = cbuffer.data(pk_geom_20_off + 643 * ccomps * dcomps);

            auto g_zz_0_z_yyyzzzz = cbuffer.data(pk_geom_20_off + 644 * ccomps * dcomps);

            auto g_zz_0_z_yyzzzzz = cbuffer.data(pk_geom_20_off + 645 * ccomps * dcomps);

            auto g_zz_0_z_yzzzzzz = cbuffer.data(pk_geom_20_off + 646 * ccomps * dcomps);

            auto g_zz_0_z_zzzzzzz = cbuffer.data(pk_geom_20_off + 647 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dixx

            const auto di_geom_20_off = idx_geom_20_dixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xx_xxxxxx = cbuffer.data(di_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xx_xxxxxy = cbuffer.data(di_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xx_xxxxxz = cbuffer.data(di_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xx_xxxxyy = cbuffer.data(di_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xx_xxxxyz = cbuffer.data(di_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xx_xxxxzz = cbuffer.data(di_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xx_xxxyyy = cbuffer.data(di_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xx_xxxyyz = cbuffer.data(di_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xx_xxxyzz = cbuffer.data(di_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xx_xxxzzz = cbuffer.data(di_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xx_xxyyyy = cbuffer.data(di_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xx_xxyyyz = cbuffer.data(di_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xx_xxyyzz = cbuffer.data(di_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xx_xxyzzz = cbuffer.data(di_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xx_xxzzzz = cbuffer.data(di_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xx_xyyyyy = cbuffer.data(di_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xx_xyyyyz = cbuffer.data(di_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xx_xyyyzz = cbuffer.data(di_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xx_xyyzzz = cbuffer.data(di_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xx_xyzzzz = cbuffer.data(di_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xx_xzzzzz = cbuffer.data(di_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xx_yyyyyy = cbuffer.data(di_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xx_yyyyyz = cbuffer.data(di_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xx_yyyyzz = cbuffer.data(di_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xx_yyyzzz = cbuffer.data(di_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xx_yyzzzz = cbuffer.data(di_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xx_yzzzzz = cbuffer.data(di_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xx_zzzzzz = cbuffer.data(di_geom_20_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxxxxx, g_x_0_x_xxxxxy, g_x_0_x_xxxxxz, g_x_0_x_xxxxyy, g_x_0_x_xxxxyz, g_x_0_x_xxxxzz, g_x_0_x_xxxyyy, g_x_0_x_xxxyyz, g_x_0_x_xxxyzz, g_x_0_x_xxxzzz, g_x_0_x_xxyyyy, g_x_0_x_xxyyyz, g_x_0_x_xxyyzz, g_x_0_x_xxyzzz, g_x_0_x_xxzzzz, g_x_0_x_xyyyyy, g_x_0_x_xyyyyz, g_x_0_x_xyyyzz, g_x_0_x_xyyzzz, g_x_0_x_xyzzzz, g_x_0_x_xzzzzz, g_x_0_x_yyyyyy, g_x_0_x_yyyyyz, g_x_0_x_yyyyzz, g_x_0_x_yyyzzz, g_x_0_x_yyzzzz, g_x_0_x_yzzzzz, g_x_0_x_zzzzzz, g_xx_0_x_xxxxxx, g_xx_0_x_xxxxxxx, g_xx_0_x_xxxxxxy, g_xx_0_x_xxxxxxz, g_xx_0_x_xxxxxy, g_xx_0_x_xxxxxyy, g_xx_0_x_xxxxxyz, g_xx_0_x_xxxxxz, g_xx_0_x_xxxxxzz, g_xx_0_x_xxxxyy, g_xx_0_x_xxxxyyy, g_xx_0_x_xxxxyyz, g_xx_0_x_xxxxyz, g_xx_0_x_xxxxyzz, g_xx_0_x_xxxxzz, g_xx_0_x_xxxxzzz, g_xx_0_x_xxxyyy, g_xx_0_x_xxxyyyy, g_xx_0_x_xxxyyyz, g_xx_0_x_xxxyyz, g_xx_0_x_xxxyyzz, g_xx_0_x_xxxyzz, g_xx_0_x_xxxyzzz, g_xx_0_x_xxxzzz, g_xx_0_x_xxxzzzz, g_xx_0_x_xxyyyy, g_xx_0_x_xxyyyyy, g_xx_0_x_xxyyyyz, g_xx_0_x_xxyyyz, g_xx_0_x_xxyyyzz, g_xx_0_x_xxyyzz, g_xx_0_x_xxyyzzz, g_xx_0_x_xxyzzz, g_xx_0_x_xxyzzzz, g_xx_0_x_xxzzzz, g_xx_0_x_xxzzzzz, g_xx_0_x_xyyyyy, g_xx_0_x_xyyyyyy, g_xx_0_x_xyyyyyz, g_xx_0_x_xyyyyz, g_xx_0_x_xyyyyzz, g_xx_0_x_xyyyzz, g_xx_0_x_xyyyzzz, g_xx_0_x_xyyzzz, g_xx_0_x_xyyzzzz, g_xx_0_x_xyzzzz, g_xx_0_x_xyzzzzz, g_xx_0_x_xzzzzz, g_xx_0_x_xzzzzzz, g_xx_0_x_yyyyyy, g_xx_0_x_yyyyyz, g_xx_0_x_yyyyzz, g_xx_0_x_yyyzzz, g_xx_0_x_yyzzzz, g_xx_0_x_yzzzzz, g_xx_0_x_zzzzzz, g_xx_0_xx_xxxxxx, g_xx_0_xx_xxxxxy, g_xx_0_xx_xxxxxz, g_xx_0_xx_xxxxyy, g_xx_0_xx_xxxxyz, g_xx_0_xx_xxxxzz, g_xx_0_xx_xxxyyy, g_xx_0_xx_xxxyyz, g_xx_0_xx_xxxyzz, g_xx_0_xx_xxxzzz, g_xx_0_xx_xxyyyy, g_xx_0_xx_xxyyyz, g_xx_0_xx_xxyyzz, g_xx_0_xx_xxyzzz, g_xx_0_xx_xxzzzz, g_xx_0_xx_xyyyyy, g_xx_0_xx_xyyyyz, g_xx_0_xx_xyyyzz, g_xx_0_xx_xyyzzz, g_xx_0_xx_xyzzzz, g_xx_0_xx_xzzzzz, g_xx_0_xx_yyyyyy, g_xx_0_xx_yyyyyz, g_xx_0_xx_yyyyzz, g_xx_0_xx_yyyzzz, g_xx_0_xx_yyzzzz, g_xx_0_xx_yzzzzz, g_xx_0_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xx_xxxxxx[k] = -2.0 * g_x_0_x_xxxxxx[k] - g_xx_0_x_xxxxxx[k] * ab_x + g_xx_0_x_xxxxxxx[k];

                g_xx_0_xx_xxxxxy[k] = -2.0 * g_x_0_x_xxxxxy[k] - g_xx_0_x_xxxxxy[k] * ab_x + g_xx_0_x_xxxxxxy[k];

                g_xx_0_xx_xxxxxz[k] = -2.0 * g_x_0_x_xxxxxz[k] - g_xx_0_x_xxxxxz[k] * ab_x + g_xx_0_x_xxxxxxz[k];

                g_xx_0_xx_xxxxyy[k] = -2.0 * g_x_0_x_xxxxyy[k] - g_xx_0_x_xxxxyy[k] * ab_x + g_xx_0_x_xxxxxyy[k];

                g_xx_0_xx_xxxxyz[k] = -2.0 * g_x_0_x_xxxxyz[k] - g_xx_0_x_xxxxyz[k] * ab_x + g_xx_0_x_xxxxxyz[k];

                g_xx_0_xx_xxxxzz[k] = -2.0 * g_x_0_x_xxxxzz[k] - g_xx_0_x_xxxxzz[k] * ab_x + g_xx_0_x_xxxxxzz[k];

                g_xx_0_xx_xxxyyy[k] = -2.0 * g_x_0_x_xxxyyy[k] - g_xx_0_x_xxxyyy[k] * ab_x + g_xx_0_x_xxxxyyy[k];

                g_xx_0_xx_xxxyyz[k] = -2.0 * g_x_0_x_xxxyyz[k] - g_xx_0_x_xxxyyz[k] * ab_x + g_xx_0_x_xxxxyyz[k];

                g_xx_0_xx_xxxyzz[k] = -2.0 * g_x_0_x_xxxyzz[k] - g_xx_0_x_xxxyzz[k] * ab_x + g_xx_0_x_xxxxyzz[k];

                g_xx_0_xx_xxxzzz[k] = -2.0 * g_x_0_x_xxxzzz[k] - g_xx_0_x_xxxzzz[k] * ab_x + g_xx_0_x_xxxxzzz[k];

                g_xx_0_xx_xxyyyy[k] = -2.0 * g_x_0_x_xxyyyy[k] - g_xx_0_x_xxyyyy[k] * ab_x + g_xx_0_x_xxxyyyy[k];

                g_xx_0_xx_xxyyyz[k] = -2.0 * g_x_0_x_xxyyyz[k] - g_xx_0_x_xxyyyz[k] * ab_x + g_xx_0_x_xxxyyyz[k];

                g_xx_0_xx_xxyyzz[k] = -2.0 * g_x_0_x_xxyyzz[k] - g_xx_0_x_xxyyzz[k] * ab_x + g_xx_0_x_xxxyyzz[k];

                g_xx_0_xx_xxyzzz[k] = -2.0 * g_x_0_x_xxyzzz[k] - g_xx_0_x_xxyzzz[k] * ab_x + g_xx_0_x_xxxyzzz[k];

                g_xx_0_xx_xxzzzz[k] = -2.0 * g_x_0_x_xxzzzz[k] - g_xx_0_x_xxzzzz[k] * ab_x + g_xx_0_x_xxxzzzz[k];

                g_xx_0_xx_xyyyyy[k] = -2.0 * g_x_0_x_xyyyyy[k] - g_xx_0_x_xyyyyy[k] * ab_x + g_xx_0_x_xxyyyyy[k];

                g_xx_0_xx_xyyyyz[k] = -2.0 * g_x_0_x_xyyyyz[k] - g_xx_0_x_xyyyyz[k] * ab_x + g_xx_0_x_xxyyyyz[k];

                g_xx_0_xx_xyyyzz[k] = -2.0 * g_x_0_x_xyyyzz[k] - g_xx_0_x_xyyyzz[k] * ab_x + g_xx_0_x_xxyyyzz[k];

                g_xx_0_xx_xyyzzz[k] = -2.0 * g_x_0_x_xyyzzz[k] - g_xx_0_x_xyyzzz[k] * ab_x + g_xx_0_x_xxyyzzz[k];

                g_xx_0_xx_xyzzzz[k] = -2.0 * g_x_0_x_xyzzzz[k] - g_xx_0_x_xyzzzz[k] * ab_x + g_xx_0_x_xxyzzzz[k];

                g_xx_0_xx_xzzzzz[k] = -2.0 * g_x_0_x_xzzzzz[k] - g_xx_0_x_xzzzzz[k] * ab_x + g_xx_0_x_xxzzzzz[k];

                g_xx_0_xx_yyyyyy[k] = -2.0 * g_x_0_x_yyyyyy[k] - g_xx_0_x_yyyyyy[k] * ab_x + g_xx_0_x_xyyyyyy[k];

                g_xx_0_xx_yyyyyz[k] = -2.0 * g_x_0_x_yyyyyz[k] - g_xx_0_x_yyyyyz[k] * ab_x + g_xx_0_x_xyyyyyz[k];

                g_xx_0_xx_yyyyzz[k] = -2.0 * g_x_0_x_yyyyzz[k] - g_xx_0_x_yyyyzz[k] * ab_x + g_xx_0_x_xyyyyzz[k];

                g_xx_0_xx_yyyzzz[k] = -2.0 * g_x_0_x_yyyzzz[k] - g_xx_0_x_yyyzzz[k] * ab_x + g_xx_0_x_xyyyzzz[k];

                g_xx_0_xx_yyzzzz[k] = -2.0 * g_x_0_x_yyzzzz[k] - g_xx_0_x_yyzzzz[k] * ab_x + g_xx_0_x_xyyzzzz[k];

                g_xx_0_xx_yzzzzz[k] = -2.0 * g_x_0_x_yzzzzz[k] - g_xx_0_x_yzzzzz[k] * ab_x + g_xx_0_x_xyzzzzz[k];

                g_xx_0_xx_zzzzzz[k] = -2.0 * g_x_0_x_zzzzzz[k] - g_xx_0_x_zzzzzz[k] * ab_x + g_xx_0_x_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xy_xxxxxx = cbuffer.data(di_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xy_xxxxxy = cbuffer.data(di_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xy_xxxxxz = cbuffer.data(di_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xy_xxxxyy = cbuffer.data(di_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xy_xxxxyz = cbuffer.data(di_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xy_xxxxzz = cbuffer.data(di_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xy_xxxyyy = cbuffer.data(di_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xy_xxxyyz = cbuffer.data(di_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xy_xxxyzz = cbuffer.data(di_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xy_xxxzzz = cbuffer.data(di_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xy_xxyyyy = cbuffer.data(di_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xy_xxyyyz = cbuffer.data(di_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xy_xxyyzz = cbuffer.data(di_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xy_xxyzzz = cbuffer.data(di_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xy_xxzzzz = cbuffer.data(di_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xy_xyyyyy = cbuffer.data(di_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xy_xyyyyz = cbuffer.data(di_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xy_xyyyzz = cbuffer.data(di_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xy_xyyzzz = cbuffer.data(di_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xy_xyzzzz = cbuffer.data(di_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xy_xzzzzz = cbuffer.data(di_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xy_yyyyyy = cbuffer.data(di_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xy_yyyyyz = cbuffer.data(di_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xy_yyyyzz = cbuffer.data(di_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xy_yyyzzz = cbuffer.data(di_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xy_yyzzzz = cbuffer.data(di_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xy_yzzzzz = cbuffer.data(di_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xy_zzzzzz = cbuffer.data(di_geom_20_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_x_xxxxxx, g_xx_0_x_xxxxxxy, g_xx_0_x_xxxxxy, g_xx_0_x_xxxxxyy, g_xx_0_x_xxxxxyz, g_xx_0_x_xxxxxz, g_xx_0_x_xxxxyy, g_xx_0_x_xxxxyyy, g_xx_0_x_xxxxyyz, g_xx_0_x_xxxxyz, g_xx_0_x_xxxxyzz, g_xx_0_x_xxxxzz, g_xx_0_x_xxxyyy, g_xx_0_x_xxxyyyy, g_xx_0_x_xxxyyyz, g_xx_0_x_xxxyyz, g_xx_0_x_xxxyyzz, g_xx_0_x_xxxyzz, g_xx_0_x_xxxyzzz, g_xx_0_x_xxxzzz, g_xx_0_x_xxyyyy, g_xx_0_x_xxyyyyy, g_xx_0_x_xxyyyyz, g_xx_0_x_xxyyyz, g_xx_0_x_xxyyyzz, g_xx_0_x_xxyyzz, g_xx_0_x_xxyyzzz, g_xx_0_x_xxyzzz, g_xx_0_x_xxyzzzz, g_xx_0_x_xxzzzz, g_xx_0_x_xyyyyy, g_xx_0_x_xyyyyyy, g_xx_0_x_xyyyyyz, g_xx_0_x_xyyyyz, g_xx_0_x_xyyyyzz, g_xx_0_x_xyyyzz, g_xx_0_x_xyyyzzz, g_xx_0_x_xyyzzz, g_xx_0_x_xyyzzzz, g_xx_0_x_xyzzzz, g_xx_0_x_xyzzzzz, g_xx_0_x_xzzzzz, g_xx_0_x_yyyyyy, g_xx_0_x_yyyyyyy, g_xx_0_x_yyyyyyz, g_xx_0_x_yyyyyz, g_xx_0_x_yyyyyzz, g_xx_0_x_yyyyzz, g_xx_0_x_yyyyzzz, g_xx_0_x_yyyzzz, g_xx_0_x_yyyzzzz, g_xx_0_x_yyzzzz, g_xx_0_x_yyzzzzz, g_xx_0_x_yzzzzz, g_xx_0_x_yzzzzzz, g_xx_0_x_zzzzzz, g_xx_0_xy_xxxxxx, g_xx_0_xy_xxxxxy, g_xx_0_xy_xxxxxz, g_xx_0_xy_xxxxyy, g_xx_0_xy_xxxxyz, g_xx_0_xy_xxxxzz, g_xx_0_xy_xxxyyy, g_xx_0_xy_xxxyyz, g_xx_0_xy_xxxyzz, g_xx_0_xy_xxxzzz, g_xx_0_xy_xxyyyy, g_xx_0_xy_xxyyyz, g_xx_0_xy_xxyyzz, g_xx_0_xy_xxyzzz, g_xx_0_xy_xxzzzz, g_xx_0_xy_xyyyyy, g_xx_0_xy_xyyyyz, g_xx_0_xy_xyyyzz, g_xx_0_xy_xyyzzz, g_xx_0_xy_xyzzzz, g_xx_0_xy_xzzzzz, g_xx_0_xy_yyyyyy, g_xx_0_xy_yyyyyz, g_xx_0_xy_yyyyzz, g_xx_0_xy_yyyzzz, g_xx_0_xy_yyzzzz, g_xx_0_xy_yzzzzz, g_xx_0_xy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xy_xxxxxx[k] = -g_xx_0_x_xxxxxx[k] * ab_y + g_xx_0_x_xxxxxxy[k];

                g_xx_0_xy_xxxxxy[k] = -g_xx_0_x_xxxxxy[k] * ab_y + g_xx_0_x_xxxxxyy[k];

                g_xx_0_xy_xxxxxz[k] = -g_xx_0_x_xxxxxz[k] * ab_y + g_xx_0_x_xxxxxyz[k];

                g_xx_0_xy_xxxxyy[k] = -g_xx_0_x_xxxxyy[k] * ab_y + g_xx_0_x_xxxxyyy[k];

                g_xx_0_xy_xxxxyz[k] = -g_xx_0_x_xxxxyz[k] * ab_y + g_xx_0_x_xxxxyyz[k];

                g_xx_0_xy_xxxxzz[k] = -g_xx_0_x_xxxxzz[k] * ab_y + g_xx_0_x_xxxxyzz[k];

                g_xx_0_xy_xxxyyy[k] = -g_xx_0_x_xxxyyy[k] * ab_y + g_xx_0_x_xxxyyyy[k];

                g_xx_0_xy_xxxyyz[k] = -g_xx_0_x_xxxyyz[k] * ab_y + g_xx_0_x_xxxyyyz[k];

                g_xx_0_xy_xxxyzz[k] = -g_xx_0_x_xxxyzz[k] * ab_y + g_xx_0_x_xxxyyzz[k];

                g_xx_0_xy_xxxzzz[k] = -g_xx_0_x_xxxzzz[k] * ab_y + g_xx_0_x_xxxyzzz[k];

                g_xx_0_xy_xxyyyy[k] = -g_xx_0_x_xxyyyy[k] * ab_y + g_xx_0_x_xxyyyyy[k];

                g_xx_0_xy_xxyyyz[k] = -g_xx_0_x_xxyyyz[k] * ab_y + g_xx_0_x_xxyyyyz[k];

                g_xx_0_xy_xxyyzz[k] = -g_xx_0_x_xxyyzz[k] * ab_y + g_xx_0_x_xxyyyzz[k];

                g_xx_0_xy_xxyzzz[k] = -g_xx_0_x_xxyzzz[k] * ab_y + g_xx_0_x_xxyyzzz[k];

                g_xx_0_xy_xxzzzz[k] = -g_xx_0_x_xxzzzz[k] * ab_y + g_xx_0_x_xxyzzzz[k];

                g_xx_0_xy_xyyyyy[k] = -g_xx_0_x_xyyyyy[k] * ab_y + g_xx_0_x_xyyyyyy[k];

                g_xx_0_xy_xyyyyz[k] = -g_xx_0_x_xyyyyz[k] * ab_y + g_xx_0_x_xyyyyyz[k];

                g_xx_0_xy_xyyyzz[k] = -g_xx_0_x_xyyyzz[k] * ab_y + g_xx_0_x_xyyyyzz[k];

                g_xx_0_xy_xyyzzz[k] = -g_xx_0_x_xyyzzz[k] * ab_y + g_xx_0_x_xyyyzzz[k];

                g_xx_0_xy_xyzzzz[k] = -g_xx_0_x_xyzzzz[k] * ab_y + g_xx_0_x_xyyzzzz[k];

                g_xx_0_xy_xzzzzz[k] = -g_xx_0_x_xzzzzz[k] * ab_y + g_xx_0_x_xyzzzzz[k];

                g_xx_0_xy_yyyyyy[k] = -g_xx_0_x_yyyyyy[k] * ab_y + g_xx_0_x_yyyyyyy[k];

                g_xx_0_xy_yyyyyz[k] = -g_xx_0_x_yyyyyz[k] * ab_y + g_xx_0_x_yyyyyyz[k];

                g_xx_0_xy_yyyyzz[k] = -g_xx_0_x_yyyyzz[k] * ab_y + g_xx_0_x_yyyyyzz[k];

                g_xx_0_xy_yyyzzz[k] = -g_xx_0_x_yyyzzz[k] * ab_y + g_xx_0_x_yyyyzzz[k];

                g_xx_0_xy_yyzzzz[k] = -g_xx_0_x_yyzzzz[k] * ab_y + g_xx_0_x_yyyzzzz[k];

                g_xx_0_xy_yzzzzz[k] = -g_xx_0_x_yzzzzz[k] * ab_y + g_xx_0_x_yyzzzzz[k];

                g_xx_0_xy_zzzzzz[k] = -g_xx_0_x_zzzzzz[k] * ab_y + g_xx_0_x_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xz_xxxxxx = cbuffer.data(di_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xz_xxxxxy = cbuffer.data(di_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xz_xxxxxz = cbuffer.data(di_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xz_xxxxyy = cbuffer.data(di_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_xz_xxxxyz = cbuffer.data(di_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xz_xxxxzz = cbuffer.data(di_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xz_xxxyyy = cbuffer.data(di_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xz_xxxyyz = cbuffer.data(di_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xz_xxxyzz = cbuffer.data(di_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xz_xxxzzz = cbuffer.data(di_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xz_xxyyyy = cbuffer.data(di_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xz_xxyyyz = cbuffer.data(di_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xz_xxyyzz = cbuffer.data(di_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xz_xxyzzz = cbuffer.data(di_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xz_xxzzzz = cbuffer.data(di_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xz_xyyyyy = cbuffer.data(di_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xz_xyyyyz = cbuffer.data(di_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xz_xyyyzz = cbuffer.data(di_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xz_xyyzzz = cbuffer.data(di_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xz_xyzzzz = cbuffer.data(di_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xz_xzzzzz = cbuffer.data(di_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xz_yyyyyy = cbuffer.data(di_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xz_yyyyyz = cbuffer.data(di_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xz_yyyyzz = cbuffer.data(di_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xz_yyyzzz = cbuffer.data(di_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xz_yyzzzz = cbuffer.data(di_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xz_yzzzzz = cbuffer.data(di_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xz_zzzzzz = cbuffer.data(di_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_x_xxxxxx, g_xx_0_x_xxxxxxz, g_xx_0_x_xxxxxy, g_xx_0_x_xxxxxyz, g_xx_0_x_xxxxxz, g_xx_0_x_xxxxxzz, g_xx_0_x_xxxxyy, g_xx_0_x_xxxxyyz, g_xx_0_x_xxxxyz, g_xx_0_x_xxxxyzz, g_xx_0_x_xxxxzz, g_xx_0_x_xxxxzzz, g_xx_0_x_xxxyyy, g_xx_0_x_xxxyyyz, g_xx_0_x_xxxyyz, g_xx_0_x_xxxyyzz, g_xx_0_x_xxxyzz, g_xx_0_x_xxxyzzz, g_xx_0_x_xxxzzz, g_xx_0_x_xxxzzzz, g_xx_0_x_xxyyyy, g_xx_0_x_xxyyyyz, g_xx_0_x_xxyyyz, g_xx_0_x_xxyyyzz, g_xx_0_x_xxyyzz, g_xx_0_x_xxyyzzz, g_xx_0_x_xxyzzz, g_xx_0_x_xxyzzzz, g_xx_0_x_xxzzzz, g_xx_0_x_xxzzzzz, g_xx_0_x_xyyyyy, g_xx_0_x_xyyyyyz, g_xx_0_x_xyyyyz, g_xx_0_x_xyyyyzz, g_xx_0_x_xyyyzz, g_xx_0_x_xyyyzzz, g_xx_0_x_xyyzzz, g_xx_0_x_xyyzzzz, g_xx_0_x_xyzzzz, g_xx_0_x_xyzzzzz, g_xx_0_x_xzzzzz, g_xx_0_x_xzzzzzz, g_xx_0_x_yyyyyy, g_xx_0_x_yyyyyyz, g_xx_0_x_yyyyyz, g_xx_0_x_yyyyyzz, g_xx_0_x_yyyyzz, g_xx_0_x_yyyyzzz, g_xx_0_x_yyyzzz, g_xx_0_x_yyyzzzz, g_xx_0_x_yyzzzz, g_xx_0_x_yyzzzzz, g_xx_0_x_yzzzzz, g_xx_0_x_yzzzzzz, g_xx_0_x_zzzzzz, g_xx_0_x_zzzzzzz, g_xx_0_xz_xxxxxx, g_xx_0_xz_xxxxxy, g_xx_0_xz_xxxxxz, g_xx_0_xz_xxxxyy, g_xx_0_xz_xxxxyz, g_xx_0_xz_xxxxzz, g_xx_0_xz_xxxyyy, g_xx_0_xz_xxxyyz, g_xx_0_xz_xxxyzz, g_xx_0_xz_xxxzzz, g_xx_0_xz_xxyyyy, g_xx_0_xz_xxyyyz, g_xx_0_xz_xxyyzz, g_xx_0_xz_xxyzzz, g_xx_0_xz_xxzzzz, g_xx_0_xz_xyyyyy, g_xx_0_xz_xyyyyz, g_xx_0_xz_xyyyzz, g_xx_0_xz_xyyzzz, g_xx_0_xz_xyzzzz, g_xx_0_xz_xzzzzz, g_xx_0_xz_yyyyyy, g_xx_0_xz_yyyyyz, g_xx_0_xz_yyyyzz, g_xx_0_xz_yyyzzz, g_xx_0_xz_yyzzzz, g_xx_0_xz_yzzzzz, g_xx_0_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xz_xxxxxx[k] = -g_xx_0_x_xxxxxx[k] * ab_z + g_xx_0_x_xxxxxxz[k];

                g_xx_0_xz_xxxxxy[k] = -g_xx_0_x_xxxxxy[k] * ab_z + g_xx_0_x_xxxxxyz[k];

                g_xx_0_xz_xxxxxz[k] = -g_xx_0_x_xxxxxz[k] * ab_z + g_xx_0_x_xxxxxzz[k];

                g_xx_0_xz_xxxxyy[k] = -g_xx_0_x_xxxxyy[k] * ab_z + g_xx_0_x_xxxxyyz[k];

                g_xx_0_xz_xxxxyz[k] = -g_xx_0_x_xxxxyz[k] * ab_z + g_xx_0_x_xxxxyzz[k];

                g_xx_0_xz_xxxxzz[k] = -g_xx_0_x_xxxxzz[k] * ab_z + g_xx_0_x_xxxxzzz[k];

                g_xx_0_xz_xxxyyy[k] = -g_xx_0_x_xxxyyy[k] * ab_z + g_xx_0_x_xxxyyyz[k];

                g_xx_0_xz_xxxyyz[k] = -g_xx_0_x_xxxyyz[k] * ab_z + g_xx_0_x_xxxyyzz[k];

                g_xx_0_xz_xxxyzz[k] = -g_xx_0_x_xxxyzz[k] * ab_z + g_xx_0_x_xxxyzzz[k];

                g_xx_0_xz_xxxzzz[k] = -g_xx_0_x_xxxzzz[k] * ab_z + g_xx_0_x_xxxzzzz[k];

                g_xx_0_xz_xxyyyy[k] = -g_xx_0_x_xxyyyy[k] * ab_z + g_xx_0_x_xxyyyyz[k];

                g_xx_0_xz_xxyyyz[k] = -g_xx_0_x_xxyyyz[k] * ab_z + g_xx_0_x_xxyyyzz[k];

                g_xx_0_xz_xxyyzz[k] = -g_xx_0_x_xxyyzz[k] * ab_z + g_xx_0_x_xxyyzzz[k];

                g_xx_0_xz_xxyzzz[k] = -g_xx_0_x_xxyzzz[k] * ab_z + g_xx_0_x_xxyzzzz[k];

                g_xx_0_xz_xxzzzz[k] = -g_xx_0_x_xxzzzz[k] * ab_z + g_xx_0_x_xxzzzzz[k];

                g_xx_0_xz_xyyyyy[k] = -g_xx_0_x_xyyyyy[k] * ab_z + g_xx_0_x_xyyyyyz[k];

                g_xx_0_xz_xyyyyz[k] = -g_xx_0_x_xyyyyz[k] * ab_z + g_xx_0_x_xyyyyzz[k];

                g_xx_0_xz_xyyyzz[k] = -g_xx_0_x_xyyyzz[k] * ab_z + g_xx_0_x_xyyyzzz[k];

                g_xx_0_xz_xyyzzz[k] = -g_xx_0_x_xyyzzz[k] * ab_z + g_xx_0_x_xyyzzzz[k];

                g_xx_0_xz_xyzzzz[k] = -g_xx_0_x_xyzzzz[k] * ab_z + g_xx_0_x_xyzzzzz[k];

                g_xx_0_xz_xzzzzz[k] = -g_xx_0_x_xzzzzz[k] * ab_z + g_xx_0_x_xzzzzzz[k];

                g_xx_0_xz_yyyyyy[k] = -g_xx_0_x_yyyyyy[k] * ab_z + g_xx_0_x_yyyyyyz[k];

                g_xx_0_xz_yyyyyz[k] = -g_xx_0_x_yyyyyz[k] * ab_z + g_xx_0_x_yyyyyzz[k];

                g_xx_0_xz_yyyyzz[k] = -g_xx_0_x_yyyyzz[k] * ab_z + g_xx_0_x_yyyyzzz[k];

                g_xx_0_xz_yyyzzz[k] = -g_xx_0_x_yyyzzz[k] * ab_z + g_xx_0_x_yyyzzzz[k];

                g_xx_0_xz_yyzzzz[k] = -g_xx_0_x_yyzzzz[k] * ab_z + g_xx_0_x_yyzzzzz[k];

                g_xx_0_xz_yzzzzz[k] = -g_xx_0_x_yzzzzz[k] * ab_z + g_xx_0_x_yzzzzzz[k];

                g_xx_0_xz_zzzzzz[k] = -g_xx_0_x_zzzzzz[k] * ab_z + g_xx_0_x_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yy_xxxxxx = cbuffer.data(di_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_yy_xxxxxy = cbuffer.data(di_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_yy_xxxxxz = cbuffer.data(di_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_yy_xxxxyy = cbuffer.data(di_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_yy_xxxxyz = cbuffer.data(di_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_yy_xxxxzz = cbuffer.data(di_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_yy_xxxyyy = cbuffer.data(di_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_yy_xxxyyz = cbuffer.data(di_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_yy_xxxyzz = cbuffer.data(di_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_yy_xxxzzz = cbuffer.data(di_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_yy_xxyyyy = cbuffer.data(di_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_yy_xxyyyz = cbuffer.data(di_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_yy_xxyyzz = cbuffer.data(di_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_yy_xxyzzz = cbuffer.data(di_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_yy_xxzzzz = cbuffer.data(di_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_yy_xyyyyy = cbuffer.data(di_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_yy_xyyyyz = cbuffer.data(di_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_yy_xyyyzz = cbuffer.data(di_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_yy_xyyzzz = cbuffer.data(di_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_yy_xyzzzz = cbuffer.data(di_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_yy_xzzzzz = cbuffer.data(di_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_yy_yyyyyy = cbuffer.data(di_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_yy_yyyyyz = cbuffer.data(di_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_yy_yyyyzz = cbuffer.data(di_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_yy_yyyzzz = cbuffer.data(di_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_yy_yyzzzz = cbuffer.data(di_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_yy_yzzzzz = cbuffer.data(di_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_yy_zzzzzz = cbuffer.data(di_geom_20_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_y_xxxxxx, g_xx_0_y_xxxxxxy, g_xx_0_y_xxxxxy, g_xx_0_y_xxxxxyy, g_xx_0_y_xxxxxyz, g_xx_0_y_xxxxxz, g_xx_0_y_xxxxyy, g_xx_0_y_xxxxyyy, g_xx_0_y_xxxxyyz, g_xx_0_y_xxxxyz, g_xx_0_y_xxxxyzz, g_xx_0_y_xxxxzz, g_xx_0_y_xxxyyy, g_xx_0_y_xxxyyyy, g_xx_0_y_xxxyyyz, g_xx_0_y_xxxyyz, g_xx_0_y_xxxyyzz, g_xx_0_y_xxxyzz, g_xx_0_y_xxxyzzz, g_xx_0_y_xxxzzz, g_xx_0_y_xxyyyy, g_xx_0_y_xxyyyyy, g_xx_0_y_xxyyyyz, g_xx_0_y_xxyyyz, g_xx_0_y_xxyyyzz, g_xx_0_y_xxyyzz, g_xx_0_y_xxyyzzz, g_xx_0_y_xxyzzz, g_xx_0_y_xxyzzzz, g_xx_0_y_xxzzzz, g_xx_0_y_xyyyyy, g_xx_0_y_xyyyyyy, g_xx_0_y_xyyyyyz, g_xx_0_y_xyyyyz, g_xx_0_y_xyyyyzz, g_xx_0_y_xyyyzz, g_xx_0_y_xyyyzzz, g_xx_0_y_xyyzzz, g_xx_0_y_xyyzzzz, g_xx_0_y_xyzzzz, g_xx_0_y_xyzzzzz, g_xx_0_y_xzzzzz, g_xx_0_y_yyyyyy, g_xx_0_y_yyyyyyy, g_xx_0_y_yyyyyyz, g_xx_0_y_yyyyyz, g_xx_0_y_yyyyyzz, g_xx_0_y_yyyyzz, g_xx_0_y_yyyyzzz, g_xx_0_y_yyyzzz, g_xx_0_y_yyyzzzz, g_xx_0_y_yyzzzz, g_xx_0_y_yyzzzzz, g_xx_0_y_yzzzzz, g_xx_0_y_yzzzzzz, g_xx_0_y_zzzzzz, g_xx_0_yy_xxxxxx, g_xx_0_yy_xxxxxy, g_xx_0_yy_xxxxxz, g_xx_0_yy_xxxxyy, g_xx_0_yy_xxxxyz, g_xx_0_yy_xxxxzz, g_xx_0_yy_xxxyyy, g_xx_0_yy_xxxyyz, g_xx_0_yy_xxxyzz, g_xx_0_yy_xxxzzz, g_xx_0_yy_xxyyyy, g_xx_0_yy_xxyyyz, g_xx_0_yy_xxyyzz, g_xx_0_yy_xxyzzz, g_xx_0_yy_xxzzzz, g_xx_0_yy_xyyyyy, g_xx_0_yy_xyyyyz, g_xx_0_yy_xyyyzz, g_xx_0_yy_xyyzzz, g_xx_0_yy_xyzzzz, g_xx_0_yy_xzzzzz, g_xx_0_yy_yyyyyy, g_xx_0_yy_yyyyyz, g_xx_0_yy_yyyyzz, g_xx_0_yy_yyyzzz, g_xx_0_yy_yyzzzz, g_xx_0_yy_yzzzzz, g_xx_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yy_xxxxxx[k] = -g_xx_0_y_xxxxxx[k] * ab_y + g_xx_0_y_xxxxxxy[k];

                g_xx_0_yy_xxxxxy[k] = -g_xx_0_y_xxxxxy[k] * ab_y + g_xx_0_y_xxxxxyy[k];

                g_xx_0_yy_xxxxxz[k] = -g_xx_0_y_xxxxxz[k] * ab_y + g_xx_0_y_xxxxxyz[k];

                g_xx_0_yy_xxxxyy[k] = -g_xx_0_y_xxxxyy[k] * ab_y + g_xx_0_y_xxxxyyy[k];

                g_xx_0_yy_xxxxyz[k] = -g_xx_0_y_xxxxyz[k] * ab_y + g_xx_0_y_xxxxyyz[k];

                g_xx_0_yy_xxxxzz[k] = -g_xx_0_y_xxxxzz[k] * ab_y + g_xx_0_y_xxxxyzz[k];

                g_xx_0_yy_xxxyyy[k] = -g_xx_0_y_xxxyyy[k] * ab_y + g_xx_0_y_xxxyyyy[k];

                g_xx_0_yy_xxxyyz[k] = -g_xx_0_y_xxxyyz[k] * ab_y + g_xx_0_y_xxxyyyz[k];

                g_xx_0_yy_xxxyzz[k] = -g_xx_0_y_xxxyzz[k] * ab_y + g_xx_0_y_xxxyyzz[k];

                g_xx_0_yy_xxxzzz[k] = -g_xx_0_y_xxxzzz[k] * ab_y + g_xx_0_y_xxxyzzz[k];

                g_xx_0_yy_xxyyyy[k] = -g_xx_0_y_xxyyyy[k] * ab_y + g_xx_0_y_xxyyyyy[k];

                g_xx_0_yy_xxyyyz[k] = -g_xx_0_y_xxyyyz[k] * ab_y + g_xx_0_y_xxyyyyz[k];

                g_xx_0_yy_xxyyzz[k] = -g_xx_0_y_xxyyzz[k] * ab_y + g_xx_0_y_xxyyyzz[k];

                g_xx_0_yy_xxyzzz[k] = -g_xx_0_y_xxyzzz[k] * ab_y + g_xx_0_y_xxyyzzz[k];

                g_xx_0_yy_xxzzzz[k] = -g_xx_0_y_xxzzzz[k] * ab_y + g_xx_0_y_xxyzzzz[k];

                g_xx_0_yy_xyyyyy[k] = -g_xx_0_y_xyyyyy[k] * ab_y + g_xx_0_y_xyyyyyy[k];

                g_xx_0_yy_xyyyyz[k] = -g_xx_0_y_xyyyyz[k] * ab_y + g_xx_0_y_xyyyyyz[k];

                g_xx_0_yy_xyyyzz[k] = -g_xx_0_y_xyyyzz[k] * ab_y + g_xx_0_y_xyyyyzz[k];

                g_xx_0_yy_xyyzzz[k] = -g_xx_0_y_xyyzzz[k] * ab_y + g_xx_0_y_xyyyzzz[k];

                g_xx_0_yy_xyzzzz[k] = -g_xx_0_y_xyzzzz[k] * ab_y + g_xx_0_y_xyyzzzz[k];

                g_xx_0_yy_xzzzzz[k] = -g_xx_0_y_xzzzzz[k] * ab_y + g_xx_0_y_xyzzzzz[k];

                g_xx_0_yy_yyyyyy[k] = -g_xx_0_y_yyyyyy[k] * ab_y + g_xx_0_y_yyyyyyy[k];

                g_xx_0_yy_yyyyyz[k] = -g_xx_0_y_yyyyyz[k] * ab_y + g_xx_0_y_yyyyyyz[k];

                g_xx_0_yy_yyyyzz[k] = -g_xx_0_y_yyyyzz[k] * ab_y + g_xx_0_y_yyyyyzz[k];

                g_xx_0_yy_yyyzzz[k] = -g_xx_0_y_yyyzzz[k] * ab_y + g_xx_0_y_yyyyzzz[k];

                g_xx_0_yy_yyzzzz[k] = -g_xx_0_y_yyzzzz[k] * ab_y + g_xx_0_y_yyyzzzz[k];

                g_xx_0_yy_yzzzzz[k] = -g_xx_0_y_yzzzzz[k] * ab_y + g_xx_0_y_yyzzzzz[k];

                g_xx_0_yy_zzzzzz[k] = -g_xx_0_y_zzzzzz[k] * ab_y + g_xx_0_y_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yz_xxxxxx = cbuffer.data(di_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_yz_xxxxxy = cbuffer.data(di_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_yz_xxxxxz = cbuffer.data(di_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_yz_xxxxyy = cbuffer.data(di_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_yz_xxxxyz = cbuffer.data(di_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_yz_xxxxzz = cbuffer.data(di_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_yz_xxxyyy = cbuffer.data(di_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_yz_xxxyyz = cbuffer.data(di_geom_20_off + 119 * ccomps * dcomps);

            auto g_xx_0_yz_xxxyzz = cbuffer.data(di_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_yz_xxxzzz = cbuffer.data(di_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_yz_xxyyyy = cbuffer.data(di_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_yz_xxyyyz = cbuffer.data(di_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_yz_xxyyzz = cbuffer.data(di_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_yz_xxyzzz = cbuffer.data(di_geom_20_off + 125 * ccomps * dcomps);

            auto g_xx_0_yz_xxzzzz = cbuffer.data(di_geom_20_off + 126 * ccomps * dcomps);

            auto g_xx_0_yz_xyyyyy = cbuffer.data(di_geom_20_off + 127 * ccomps * dcomps);

            auto g_xx_0_yz_xyyyyz = cbuffer.data(di_geom_20_off + 128 * ccomps * dcomps);

            auto g_xx_0_yz_xyyyzz = cbuffer.data(di_geom_20_off + 129 * ccomps * dcomps);

            auto g_xx_0_yz_xyyzzz = cbuffer.data(di_geom_20_off + 130 * ccomps * dcomps);

            auto g_xx_0_yz_xyzzzz = cbuffer.data(di_geom_20_off + 131 * ccomps * dcomps);

            auto g_xx_0_yz_xzzzzz = cbuffer.data(di_geom_20_off + 132 * ccomps * dcomps);

            auto g_xx_0_yz_yyyyyy = cbuffer.data(di_geom_20_off + 133 * ccomps * dcomps);

            auto g_xx_0_yz_yyyyyz = cbuffer.data(di_geom_20_off + 134 * ccomps * dcomps);

            auto g_xx_0_yz_yyyyzz = cbuffer.data(di_geom_20_off + 135 * ccomps * dcomps);

            auto g_xx_0_yz_yyyzzz = cbuffer.data(di_geom_20_off + 136 * ccomps * dcomps);

            auto g_xx_0_yz_yyzzzz = cbuffer.data(di_geom_20_off + 137 * ccomps * dcomps);

            auto g_xx_0_yz_yzzzzz = cbuffer.data(di_geom_20_off + 138 * ccomps * dcomps);

            auto g_xx_0_yz_zzzzzz = cbuffer.data(di_geom_20_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yz_xxxxxx, g_xx_0_yz_xxxxxy, g_xx_0_yz_xxxxxz, g_xx_0_yz_xxxxyy, g_xx_0_yz_xxxxyz, g_xx_0_yz_xxxxzz, g_xx_0_yz_xxxyyy, g_xx_0_yz_xxxyyz, g_xx_0_yz_xxxyzz, g_xx_0_yz_xxxzzz, g_xx_0_yz_xxyyyy, g_xx_0_yz_xxyyyz, g_xx_0_yz_xxyyzz, g_xx_0_yz_xxyzzz, g_xx_0_yz_xxzzzz, g_xx_0_yz_xyyyyy, g_xx_0_yz_xyyyyz, g_xx_0_yz_xyyyzz, g_xx_0_yz_xyyzzz, g_xx_0_yz_xyzzzz, g_xx_0_yz_xzzzzz, g_xx_0_yz_yyyyyy, g_xx_0_yz_yyyyyz, g_xx_0_yz_yyyyzz, g_xx_0_yz_yyyzzz, g_xx_0_yz_yyzzzz, g_xx_0_yz_yzzzzz, g_xx_0_yz_zzzzzz, g_xx_0_z_xxxxxx, g_xx_0_z_xxxxxxy, g_xx_0_z_xxxxxy, g_xx_0_z_xxxxxyy, g_xx_0_z_xxxxxyz, g_xx_0_z_xxxxxz, g_xx_0_z_xxxxyy, g_xx_0_z_xxxxyyy, g_xx_0_z_xxxxyyz, g_xx_0_z_xxxxyz, g_xx_0_z_xxxxyzz, g_xx_0_z_xxxxzz, g_xx_0_z_xxxyyy, g_xx_0_z_xxxyyyy, g_xx_0_z_xxxyyyz, g_xx_0_z_xxxyyz, g_xx_0_z_xxxyyzz, g_xx_0_z_xxxyzz, g_xx_0_z_xxxyzzz, g_xx_0_z_xxxzzz, g_xx_0_z_xxyyyy, g_xx_0_z_xxyyyyy, g_xx_0_z_xxyyyyz, g_xx_0_z_xxyyyz, g_xx_0_z_xxyyyzz, g_xx_0_z_xxyyzz, g_xx_0_z_xxyyzzz, g_xx_0_z_xxyzzz, g_xx_0_z_xxyzzzz, g_xx_0_z_xxzzzz, g_xx_0_z_xyyyyy, g_xx_0_z_xyyyyyy, g_xx_0_z_xyyyyyz, g_xx_0_z_xyyyyz, g_xx_0_z_xyyyyzz, g_xx_0_z_xyyyzz, g_xx_0_z_xyyyzzz, g_xx_0_z_xyyzzz, g_xx_0_z_xyyzzzz, g_xx_0_z_xyzzzz, g_xx_0_z_xyzzzzz, g_xx_0_z_xzzzzz, g_xx_0_z_yyyyyy, g_xx_0_z_yyyyyyy, g_xx_0_z_yyyyyyz, g_xx_0_z_yyyyyz, g_xx_0_z_yyyyyzz, g_xx_0_z_yyyyzz, g_xx_0_z_yyyyzzz, g_xx_0_z_yyyzzz, g_xx_0_z_yyyzzzz, g_xx_0_z_yyzzzz, g_xx_0_z_yyzzzzz, g_xx_0_z_yzzzzz, g_xx_0_z_yzzzzzz, g_xx_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yz_xxxxxx[k] = -g_xx_0_z_xxxxxx[k] * ab_y + g_xx_0_z_xxxxxxy[k];

                g_xx_0_yz_xxxxxy[k] = -g_xx_0_z_xxxxxy[k] * ab_y + g_xx_0_z_xxxxxyy[k];

                g_xx_0_yz_xxxxxz[k] = -g_xx_0_z_xxxxxz[k] * ab_y + g_xx_0_z_xxxxxyz[k];

                g_xx_0_yz_xxxxyy[k] = -g_xx_0_z_xxxxyy[k] * ab_y + g_xx_0_z_xxxxyyy[k];

                g_xx_0_yz_xxxxyz[k] = -g_xx_0_z_xxxxyz[k] * ab_y + g_xx_0_z_xxxxyyz[k];

                g_xx_0_yz_xxxxzz[k] = -g_xx_0_z_xxxxzz[k] * ab_y + g_xx_0_z_xxxxyzz[k];

                g_xx_0_yz_xxxyyy[k] = -g_xx_0_z_xxxyyy[k] * ab_y + g_xx_0_z_xxxyyyy[k];

                g_xx_0_yz_xxxyyz[k] = -g_xx_0_z_xxxyyz[k] * ab_y + g_xx_0_z_xxxyyyz[k];

                g_xx_0_yz_xxxyzz[k] = -g_xx_0_z_xxxyzz[k] * ab_y + g_xx_0_z_xxxyyzz[k];

                g_xx_0_yz_xxxzzz[k] = -g_xx_0_z_xxxzzz[k] * ab_y + g_xx_0_z_xxxyzzz[k];

                g_xx_0_yz_xxyyyy[k] = -g_xx_0_z_xxyyyy[k] * ab_y + g_xx_0_z_xxyyyyy[k];

                g_xx_0_yz_xxyyyz[k] = -g_xx_0_z_xxyyyz[k] * ab_y + g_xx_0_z_xxyyyyz[k];

                g_xx_0_yz_xxyyzz[k] = -g_xx_0_z_xxyyzz[k] * ab_y + g_xx_0_z_xxyyyzz[k];

                g_xx_0_yz_xxyzzz[k] = -g_xx_0_z_xxyzzz[k] * ab_y + g_xx_0_z_xxyyzzz[k];

                g_xx_0_yz_xxzzzz[k] = -g_xx_0_z_xxzzzz[k] * ab_y + g_xx_0_z_xxyzzzz[k];

                g_xx_0_yz_xyyyyy[k] = -g_xx_0_z_xyyyyy[k] * ab_y + g_xx_0_z_xyyyyyy[k];

                g_xx_0_yz_xyyyyz[k] = -g_xx_0_z_xyyyyz[k] * ab_y + g_xx_0_z_xyyyyyz[k];

                g_xx_0_yz_xyyyzz[k] = -g_xx_0_z_xyyyzz[k] * ab_y + g_xx_0_z_xyyyyzz[k];

                g_xx_0_yz_xyyzzz[k] = -g_xx_0_z_xyyzzz[k] * ab_y + g_xx_0_z_xyyyzzz[k];

                g_xx_0_yz_xyzzzz[k] = -g_xx_0_z_xyzzzz[k] * ab_y + g_xx_0_z_xyyzzzz[k];

                g_xx_0_yz_xzzzzz[k] = -g_xx_0_z_xzzzzz[k] * ab_y + g_xx_0_z_xyzzzzz[k];

                g_xx_0_yz_yyyyyy[k] = -g_xx_0_z_yyyyyy[k] * ab_y + g_xx_0_z_yyyyyyy[k];

                g_xx_0_yz_yyyyyz[k] = -g_xx_0_z_yyyyyz[k] * ab_y + g_xx_0_z_yyyyyyz[k];

                g_xx_0_yz_yyyyzz[k] = -g_xx_0_z_yyyyzz[k] * ab_y + g_xx_0_z_yyyyyzz[k];

                g_xx_0_yz_yyyzzz[k] = -g_xx_0_z_yyyzzz[k] * ab_y + g_xx_0_z_yyyyzzz[k];

                g_xx_0_yz_yyzzzz[k] = -g_xx_0_z_yyzzzz[k] * ab_y + g_xx_0_z_yyyzzzz[k];

                g_xx_0_yz_yzzzzz[k] = -g_xx_0_z_yzzzzz[k] * ab_y + g_xx_0_z_yyzzzzz[k];

                g_xx_0_yz_zzzzzz[k] = -g_xx_0_z_zzzzzz[k] * ab_y + g_xx_0_z_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zz_xxxxxx = cbuffer.data(di_geom_20_off + 140 * ccomps * dcomps);

            auto g_xx_0_zz_xxxxxy = cbuffer.data(di_geom_20_off + 141 * ccomps * dcomps);

            auto g_xx_0_zz_xxxxxz = cbuffer.data(di_geom_20_off + 142 * ccomps * dcomps);

            auto g_xx_0_zz_xxxxyy = cbuffer.data(di_geom_20_off + 143 * ccomps * dcomps);

            auto g_xx_0_zz_xxxxyz = cbuffer.data(di_geom_20_off + 144 * ccomps * dcomps);

            auto g_xx_0_zz_xxxxzz = cbuffer.data(di_geom_20_off + 145 * ccomps * dcomps);

            auto g_xx_0_zz_xxxyyy = cbuffer.data(di_geom_20_off + 146 * ccomps * dcomps);

            auto g_xx_0_zz_xxxyyz = cbuffer.data(di_geom_20_off + 147 * ccomps * dcomps);

            auto g_xx_0_zz_xxxyzz = cbuffer.data(di_geom_20_off + 148 * ccomps * dcomps);

            auto g_xx_0_zz_xxxzzz = cbuffer.data(di_geom_20_off + 149 * ccomps * dcomps);

            auto g_xx_0_zz_xxyyyy = cbuffer.data(di_geom_20_off + 150 * ccomps * dcomps);

            auto g_xx_0_zz_xxyyyz = cbuffer.data(di_geom_20_off + 151 * ccomps * dcomps);

            auto g_xx_0_zz_xxyyzz = cbuffer.data(di_geom_20_off + 152 * ccomps * dcomps);

            auto g_xx_0_zz_xxyzzz = cbuffer.data(di_geom_20_off + 153 * ccomps * dcomps);

            auto g_xx_0_zz_xxzzzz = cbuffer.data(di_geom_20_off + 154 * ccomps * dcomps);

            auto g_xx_0_zz_xyyyyy = cbuffer.data(di_geom_20_off + 155 * ccomps * dcomps);

            auto g_xx_0_zz_xyyyyz = cbuffer.data(di_geom_20_off + 156 * ccomps * dcomps);

            auto g_xx_0_zz_xyyyzz = cbuffer.data(di_geom_20_off + 157 * ccomps * dcomps);

            auto g_xx_0_zz_xyyzzz = cbuffer.data(di_geom_20_off + 158 * ccomps * dcomps);

            auto g_xx_0_zz_xyzzzz = cbuffer.data(di_geom_20_off + 159 * ccomps * dcomps);

            auto g_xx_0_zz_xzzzzz = cbuffer.data(di_geom_20_off + 160 * ccomps * dcomps);

            auto g_xx_0_zz_yyyyyy = cbuffer.data(di_geom_20_off + 161 * ccomps * dcomps);

            auto g_xx_0_zz_yyyyyz = cbuffer.data(di_geom_20_off + 162 * ccomps * dcomps);

            auto g_xx_0_zz_yyyyzz = cbuffer.data(di_geom_20_off + 163 * ccomps * dcomps);

            auto g_xx_0_zz_yyyzzz = cbuffer.data(di_geom_20_off + 164 * ccomps * dcomps);

            auto g_xx_0_zz_yyzzzz = cbuffer.data(di_geom_20_off + 165 * ccomps * dcomps);

            auto g_xx_0_zz_yzzzzz = cbuffer.data(di_geom_20_off + 166 * ccomps * dcomps);

            auto g_xx_0_zz_zzzzzz = cbuffer.data(di_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_z_xxxxxx, g_xx_0_z_xxxxxxz, g_xx_0_z_xxxxxy, g_xx_0_z_xxxxxyz, g_xx_0_z_xxxxxz, g_xx_0_z_xxxxxzz, g_xx_0_z_xxxxyy, g_xx_0_z_xxxxyyz, g_xx_0_z_xxxxyz, g_xx_0_z_xxxxyzz, g_xx_0_z_xxxxzz, g_xx_0_z_xxxxzzz, g_xx_0_z_xxxyyy, g_xx_0_z_xxxyyyz, g_xx_0_z_xxxyyz, g_xx_0_z_xxxyyzz, g_xx_0_z_xxxyzz, g_xx_0_z_xxxyzzz, g_xx_0_z_xxxzzz, g_xx_0_z_xxxzzzz, g_xx_0_z_xxyyyy, g_xx_0_z_xxyyyyz, g_xx_0_z_xxyyyz, g_xx_0_z_xxyyyzz, g_xx_0_z_xxyyzz, g_xx_0_z_xxyyzzz, g_xx_0_z_xxyzzz, g_xx_0_z_xxyzzzz, g_xx_0_z_xxzzzz, g_xx_0_z_xxzzzzz, g_xx_0_z_xyyyyy, g_xx_0_z_xyyyyyz, g_xx_0_z_xyyyyz, g_xx_0_z_xyyyyzz, g_xx_0_z_xyyyzz, g_xx_0_z_xyyyzzz, g_xx_0_z_xyyzzz, g_xx_0_z_xyyzzzz, g_xx_0_z_xyzzzz, g_xx_0_z_xyzzzzz, g_xx_0_z_xzzzzz, g_xx_0_z_xzzzzzz, g_xx_0_z_yyyyyy, g_xx_0_z_yyyyyyz, g_xx_0_z_yyyyyz, g_xx_0_z_yyyyyzz, g_xx_0_z_yyyyzz, g_xx_0_z_yyyyzzz, g_xx_0_z_yyyzzz, g_xx_0_z_yyyzzzz, g_xx_0_z_yyzzzz, g_xx_0_z_yyzzzzz, g_xx_0_z_yzzzzz, g_xx_0_z_yzzzzzz, g_xx_0_z_zzzzzz, g_xx_0_z_zzzzzzz, g_xx_0_zz_xxxxxx, g_xx_0_zz_xxxxxy, g_xx_0_zz_xxxxxz, g_xx_0_zz_xxxxyy, g_xx_0_zz_xxxxyz, g_xx_0_zz_xxxxzz, g_xx_0_zz_xxxyyy, g_xx_0_zz_xxxyyz, g_xx_0_zz_xxxyzz, g_xx_0_zz_xxxzzz, g_xx_0_zz_xxyyyy, g_xx_0_zz_xxyyyz, g_xx_0_zz_xxyyzz, g_xx_0_zz_xxyzzz, g_xx_0_zz_xxzzzz, g_xx_0_zz_xyyyyy, g_xx_0_zz_xyyyyz, g_xx_0_zz_xyyyzz, g_xx_0_zz_xyyzzz, g_xx_0_zz_xyzzzz, g_xx_0_zz_xzzzzz, g_xx_0_zz_yyyyyy, g_xx_0_zz_yyyyyz, g_xx_0_zz_yyyyzz, g_xx_0_zz_yyyzzz, g_xx_0_zz_yyzzzz, g_xx_0_zz_yzzzzz, g_xx_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zz_xxxxxx[k] = -g_xx_0_z_xxxxxx[k] * ab_z + g_xx_0_z_xxxxxxz[k];

                g_xx_0_zz_xxxxxy[k] = -g_xx_0_z_xxxxxy[k] * ab_z + g_xx_0_z_xxxxxyz[k];

                g_xx_0_zz_xxxxxz[k] = -g_xx_0_z_xxxxxz[k] * ab_z + g_xx_0_z_xxxxxzz[k];

                g_xx_0_zz_xxxxyy[k] = -g_xx_0_z_xxxxyy[k] * ab_z + g_xx_0_z_xxxxyyz[k];

                g_xx_0_zz_xxxxyz[k] = -g_xx_0_z_xxxxyz[k] * ab_z + g_xx_0_z_xxxxyzz[k];

                g_xx_0_zz_xxxxzz[k] = -g_xx_0_z_xxxxzz[k] * ab_z + g_xx_0_z_xxxxzzz[k];

                g_xx_0_zz_xxxyyy[k] = -g_xx_0_z_xxxyyy[k] * ab_z + g_xx_0_z_xxxyyyz[k];

                g_xx_0_zz_xxxyyz[k] = -g_xx_0_z_xxxyyz[k] * ab_z + g_xx_0_z_xxxyyzz[k];

                g_xx_0_zz_xxxyzz[k] = -g_xx_0_z_xxxyzz[k] * ab_z + g_xx_0_z_xxxyzzz[k];

                g_xx_0_zz_xxxzzz[k] = -g_xx_0_z_xxxzzz[k] * ab_z + g_xx_0_z_xxxzzzz[k];

                g_xx_0_zz_xxyyyy[k] = -g_xx_0_z_xxyyyy[k] * ab_z + g_xx_0_z_xxyyyyz[k];

                g_xx_0_zz_xxyyyz[k] = -g_xx_0_z_xxyyyz[k] * ab_z + g_xx_0_z_xxyyyzz[k];

                g_xx_0_zz_xxyyzz[k] = -g_xx_0_z_xxyyzz[k] * ab_z + g_xx_0_z_xxyyzzz[k];

                g_xx_0_zz_xxyzzz[k] = -g_xx_0_z_xxyzzz[k] * ab_z + g_xx_0_z_xxyzzzz[k];

                g_xx_0_zz_xxzzzz[k] = -g_xx_0_z_xxzzzz[k] * ab_z + g_xx_0_z_xxzzzzz[k];

                g_xx_0_zz_xyyyyy[k] = -g_xx_0_z_xyyyyy[k] * ab_z + g_xx_0_z_xyyyyyz[k];

                g_xx_0_zz_xyyyyz[k] = -g_xx_0_z_xyyyyz[k] * ab_z + g_xx_0_z_xyyyyzz[k];

                g_xx_0_zz_xyyyzz[k] = -g_xx_0_z_xyyyzz[k] * ab_z + g_xx_0_z_xyyyzzz[k];

                g_xx_0_zz_xyyzzz[k] = -g_xx_0_z_xyyzzz[k] * ab_z + g_xx_0_z_xyyzzzz[k];

                g_xx_0_zz_xyzzzz[k] = -g_xx_0_z_xyzzzz[k] * ab_z + g_xx_0_z_xyzzzzz[k];

                g_xx_0_zz_xzzzzz[k] = -g_xx_0_z_xzzzzz[k] * ab_z + g_xx_0_z_xzzzzzz[k];

                g_xx_0_zz_yyyyyy[k] = -g_xx_0_z_yyyyyy[k] * ab_z + g_xx_0_z_yyyyyyz[k];

                g_xx_0_zz_yyyyyz[k] = -g_xx_0_z_yyyyyz[k] * ab_z + g_xx_0_z_yyyyyzz[k];

                g_xx_0_zz_yyyyzz[k] = -g_xx_0_z_yyyyzz[k] * ab_z + g_xx_0_z_yyyyzzz[k];

                g_xx_0_zz_yyyzzz[k] = -g_xx_0_z_yyyzzz[k] * ab_z + g_xx_0_z_yyyzzzz[k];

                g_xx_0_zz_yyzzzz[k] = -g_xx_0_z_yyzzzz[k] * ab_z + g_xx_0_z_yyzzzzz[k];

                g_xx_0_zz_yzzzzz[k] = -g_xx_0_z_yzzzzz[k] * ab_z + g_xx_0_z_yzzzzzz[k];

                g_xx_0_zz_zzzzzz[k] = -g_xx_0_z_zzzzzz[k] * ab_z + g_xx_0_z_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xx_xxxxxx = cbuffer.data(di_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xx_xxxxxy = cbuffer.data(di_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_xx_xxxxxz = cbuffer.data(di_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xx_xxxxyy = cbuffer.data(di_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xx_xxxxyz = cbuffer.data(di_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xx_xxxxzz = cbuffer.data(di_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_xx_xxxyyy = cbuffer.data(di_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xx_xxxyyz = cbuffer.data(di_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xx_xxxyzz = cbuffer.data(di_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xx_xxxzzz = cbuffer.data(di_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xx_xxyyyy = cbuffer.data(di_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xx_xxyyyz = cbuffer.data(di_geom_20_off + 179 * ccomps * dcomps);

            auto g_xy_0_xx_xxyyzz = cbuffer.data(di_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xx_xxyzzz = cbuffer.data(di_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xx_xxzzzz = cbuffer.data(di_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xx_xyyyyy = cbuffer.data(di_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xx_xyyyyz = cbuffer.data(di_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xx_xyyyzz = cbuffer.data(di_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_xx_xyyzzz = cbuffer.data(di_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xx_xyzzzz = cbuffer.data(di_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xx_xzzzzz = cbuffer.data(di_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_xx_yyyyyy = cbuffer.data(di_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_xx_yyyyyz = cbuffer.data(di_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_xx_yyyyzz = cbuffer.data(di_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_xx_yyyzzz = cbuffer.data(di_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_xx_yyzzzz = cbuffer.data(di_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_xx_yzzzzz = cbuffer.data(di_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_xx_zzzzzz = cbuffer.data(di_geom_20_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_x_xxxxxx, g_xy_0_x_xxxxxxx, g_xy_0_x_xxxxxxy, g_xy_0_x_xxxxxxz, g_xy_0_x_xxxxxy, g_xy_0_x_xxxxxyy, g_xy_0_x_xxxxxyz, g_xy_0_x_xxxxxz, g_xy_0_x_xxxxxzz, g_xy_0_x_xxxxyy, g_xy_0_x_xxxxyyy, g_xy_0_x_xxxxyyz, g_xy_0_x_xxxxyz, g_xy_0_x_xxxxyzz, g_xy_0_x_xxxxzz, g_xy_0_x_xxxxzzz, g_xy_0_x_xxxyyy, g_xy_0_x_xxxyyyy, g_xy_0_x_xxxyyyz, g_xy_0_x_xxxyyz, g_xy_0_x_xxxyyzz, g_xy_0_x_xxxyzz, g_xy_0_x_xxxyzzz, g_xy_0_x_xxxzzz, g_xy_0_x_xxxzzzz, g_xy_0_x_xxyyyy, g_xy_0_x_xxyyyyy, g_xy_0_x_xxyyyyz, g_xy_0_x_xxyyyz, g_xy_0_x_xxyyyzz, g_xy_0_x_xxyyzz, g_xy_0_x_xxyyzzz, g_xy_0_x_xxyzzz, g_xy_0_x_xxyzzzz, g_xy_0_x_xxzzzz, g_xy_0_x_xxzzzzz, g_xy_0_x_xyyyyy, g_xy_0_x_xyyyyyy, g_xy_0_x_xyyyyyz, g_xy_0_x_xyyyyz, g_xy_0_x_xyyyyzz, g_xy_0_x_xyyyzz, g_xy_0_x_xyyyzzz, g_xy_0_x_xyyzzz, g_xy_0_x_xyyzzzz, g_xy_0_x_xyzzzz, g_xy_0_x_xyzzzzz, g_xy_0_x_xzzzzz, g_xy_0_x_xzzzzzz, g_xy_0_x_yyyyyy, g_xy_0_x_yyyyyz, g_xy_0_x_yyyyzz, g_xy_0_x_yyyzzz, g_xy_0_x_yyzzzz, g_xy_0_x_yzzzzz, g_xy_0_x_zzzzzz, g_xy_0_xx_xxxxxx, g_xy_0_xx_xxxxxy, g_xy_0_xx_xxxxxz, g_xy_0_xx_xxxxyy, g_xy_0_xx_xxxxyz, g_xy_0_xx_xxxxzz, g_xy_0_xx_xxxyyy, g_xy_0_xx_xxxyyz, g_xy_0_xx_xxxyzz, g_xy_0_xx_xxxzzz, g_xy_0_xx_xxyyyy, g_xy_0_xx_xxyyyz, g_xy_0_xx_xxyyzz, g_xy_0_xx_xxyzzz, g_xy_0_xx_xxzzzz, g_xy_0_xx_xyyyyy, g_xy_0_xx_xyyyyz, g_xy_0_xx_xyyyzz, g_xy_0_xx_xyyzzz, g_xy_0_xx_xyzzzz, g_xy_0_xx_xzzzzz, g_xy_0_xx_yyyyyy, g_xy_0_xx_yyyyyz, g_xy_0_xx_yyyyzz, g_xy_0_xx_yyyzzz, g_xy_0_xx_yyzzzz, g_xy_0_xx_yzzzzz, g_xy_0_xx_zzzzzz, g_y_0_x_xxxxxx, g_y_0_x_xxxxxy, g_y_0_x_xxxxxz, g_y_0_x_xxxxyy, g_y_0_x_xxxxyz, g_y_0_x_xxxxzz, g_y_0_x_xxxyyy, g_y_0_x_xxxyyz, g_y_0_x_xxxyzz, g_y_0_x_xxxzzz, g_y_0_x_xxyyyy, g_y_0_x_xxyyyz, g_y_0_x_xxyyzz, g_y_0_x_xxyzzz, g_y_0_x_xxzzzz, g_y_0_x_xyyyyy, g_y_0_x_xyyyyz, g_y_0_x_xyyyzz, g_y_0_x_xyyzzz, g_y_0_x_xyzzzz, g_y_0_x_xzzzzz, g_y_0_x_yyyyyy, g_y_0_x_yyyyyz, g_y_0_x_yyyyzz, g_y_0_x_yyyzzz, g_y_0_x_yyzzzz, g_y_0_x_yzzzzz, g_y_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xx_xxxxxx[k] = -g_y_0_x_xxxxxx[k] - g_xy_0_x_xxxxxx[k] * ab_x + g_xy_0_x_xxxxxxx[k];

                g_xy_0_xx_xxxxxy[k] = -g_y_0_x_xxxxxy[k] - g_xy_0_x_xxxxxy[k] * ab_x + g_xy_0_x_xxxxxxy[k];

                g_xy_0_xx_xxxxxz[k] = -g_y_0_x_xxxxxz[k] - g_xy_0_x_xxxxxz[k] * ab_x + g_xy_0_x_xxxxxxz[k];

                g_xy_0_xx_xxxxyy[k] = -g_y_0_x_xxxxyy[k] - g_xy_0_x_xxxxyy[k] * ab_x + g_xy_0_x_xxxxxyy[k];

                g_xy_0_xx_xxxxyz[k] = -g_y_0_x_xxxxyz[k] - g_xy_0_x_xxxxyz[k] * ab_x + g_xy_0_x_xxxxxyz[k];

                g_xy_0_xx_xxxxzz[k] = -g_y_0_x_xxxxzz[k] - g_xy_0_x_xxxxzz[k] * ab_x + g_xy_0_x_xxxxxzz[k];

                g_xy_0_xx_xxxyyy[k] = -g_y_0_x_xxxyyy[k] - g_xy_0_x_xxxyyy[k] * ab_x + g_xy_0_x_xxxxyyy[k];

                g_xy_0_xx_xxxyyz[k] = -g_y_0_x_xxxyyz[k] - g_xy_0_x_xxxyyz[k] * ab_x + g_xy_0_x_xxxxyyz[k];

                g_xy_0_xx_xxxyzz[k] = -g_y_0_x_xxxyzz[k] - g_xy_0_x_xxxyzz[k] * ab_x + g_xy_0_x_xxxxyzz[k];

                g_xy_0_xx_xxxzzz[k] = -g_y_0_x_xxxzzz[k] - g_xy_0_x_xxxzzz[k] * ab_x + g_xy_0_x_xxxxzzz[k];

                g_xy_0_xx_xxyyyy[k] = -g_y_0_x_xxyyyy[k] - g_xy_0_x_xxyyyy[k] * ab_x + g_xy_0_x_xxxyyyy[k];

                g_xy_0_xx_xxyyyz[k] = -g_y_0_x_xxyyyz[k] - g_xy_0_x_xxyyyz[k] * ab_x + g_xy_0_x_xxxyyyz[k];

                g_xy_0_xx_xxyyzz[k] = -g_y_0_x_xxyyzz[k] - g_xy_0_x_xxyyzz[k] * ab_x + g_xy_0_x_xxxyyzz[k];

                g_xy_0_xx_xxyzzz[k] = -g_y_0_x_xxyzzz[k] - g_xy_0_x_xxyzzz[k] * ab_x + g_xy_0_x_xxxyzzz[k];

                g_xy_0_xx_xxzzzz[k] = -g_y_0_x_xxzzzz[k] - g_xy_0_x_xxzzzz[k] * ab_x + g_xy_0_x_xxxzzzz[k];

                g_xy_0_xx_xyyyyy[k] = -g_y_0_x_xyyyyy[k] - g_xy_0_x_xyyyyy[k] * ab_x + g_xy_0_x_xxyyyyy[k];

                g_xy_0_xx_xyyyyz[k] = -g_y_0_x_xyyyyz[k] - g_xy_0_x_xyyyyz[k] * ab_x + g_xy_0_x_xxyyyyz[k];

                g_xy_0_xx_xyyyzz[k] = -g_y_0_x_xyyyzz[k] - g_xy_0_x_xyyyzz[k] * ab_x + g_xy_0_x_xxyyyzz[k];

                g_xy_0_xx_xyyzzz[k] = -g_y_0_x_xyyzzz[k] - g_xy_0_x_xyyzzz[k] * ab_x + g_xy_0_x_xxyyzzz[k];

                g_xy_0_xx_xyzzzz[k] = -g_y_0_x_xyzzzz[k] - g_xy_0_x_xyzzzz[k] * ab_x + g_xy_0_x_xxyzzzz[k];

                g_xy_0_xx_xzzzzz[k] = -g_y_0_x_xzzzzz[k] - g_xy_0_x_xzzzzz[k] * ab_x + g_xy_0_x_xxzzzzz[k];

                g_xy_0_xx_yyyyyy[k] = -g_y_0_x_yyyyyy[k] - g_xy_0_x_yyyyyy[k] * ab_x + g_xy_0_x_xyyyyyy[k];

                g_xy_0_xx_yyyyyz[k] = -g_y_0_x_yyyyyz[k] - g_xy_0_x_yyyyyz[k] * ab_x + g_xy_0_x_xyyyyyz[k];

                g_xy_0_xx_yyyyzz[k] = -g_y_0_x_yyyyzz[k] - g_xy_0_x_yyyyzz[k] * ab_x + g_xy_0_x_xyyyyzz[k];

                g_xy_0_xx_yyyzzz[k] = -g_y_0_x_yyyzzz[k] - g_xy_0_x_yyyzzz[k] * ab_x + g_xy_0_x_xyyyzzz[k];

                g_xy_0_xx_yyzzzz[k] = -g_y_0_x_yyzzzz[k] - g_xy_0_x_yyzzzz[k] * ab_x + g_xy_0_x_xyyzzzz[k];

                g_xy_0_xx_yzzzzz[k] = -g_y_0_x_yzzzzz[k] - g_xy_0_x_yzzzzz[k] * ab_x + g_xy_0_x_xyzzzzz[k];

                g_xy_0_xx_zzzzzz[k] = -g_y_0_x_zzzzzz[k] - g_xy_0_x_zzzzzz[k] * ab_x + g_xy_0_x_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xy_xxxxxx = cbuffer.data(di_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_xy_xxxxxy = cbuffer.data(di_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_xy_xxxxxz = cbuffer.data(di_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_xy_xxxxyy = cbuffer.data(di_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_xy_xxxxyz = cbuffer.data(di_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_xy_xxxxzz = cbuffer.data(di_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_xy_xxxyyy = cbuffer.data(di_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_xy_xxxyyz = cbuffer.data(di_geom_20_off + 203 * ccomps * dcomps);

            auto g_xy_0_xy_xxxyzz = cbuffer.data(di_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_xy_xxxzzz = cbuffer.data(di_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_xy_xxyyyy = cbuffer.data(di_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_xy_xxyyyz = cbuffer.data(di_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_xy_xxyyzz = cbuffer.data(di_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_xy_xxyzzz = cbuffer.data(di_geom_20_off + 209 * ccomps * dcomps);

            auto g_xy_0_xy_xxzzzz = cbuffer.data(di_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_xy_xyyyyy = cbuffer.data(di_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_xy_xyyyyz = cbuffer.data(di_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_xy_xyyyzz = cbuffer.data(di_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_xy_xyyzzz = cbuffer.data(di_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_xy_xyzzzz = cbuffer.data(di_geom_20_off + 215 * ccomps * dcomps);

            auto g_xy_0_xy_xzzzzz = cbuffer.data(di_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_xy_yyyyyy = cbuffer.data(di_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_xy_yyyyyz = cbuffer.data(di_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_xy_yyyyzz = cbuffer.data(di_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_xy_yyyzzz = cbuffer.data(di_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_xy_yyzzzz = cbuffer.data(di_geom_20_off + 221 * ccomps * dcomps);

            auto g_xy_0_xy_yzzzzz = cbuffer.data(di_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_xy_zzzzzz = cbuffer.data(di_geom_20_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xy_xxxxxx, g_xy_0_xy_xxxxxy, g_xy_0_xy_xxxxxz, g_xy_0_xy_xxxxyy, g_xy_0_xy_xxxxyz, g_xy_0_xy_xxxxzz, g_xy_0_xy_xxxyyy, g_xy_0_xy_xxxyyz, g_xy_0_xy_xxxyzz, g_xy_0_xy_xxxzzz, g_xy_0_xy_xxyyyy, g_xy_0_xy_xxyyyz, g_xy_0_xy_xxyyzz, g_xy_0_xy_xxyzzz, g_xy_0_xy_xxzzzz, g_xy_0_xy_xyyyyy, g_xy_0_xy_xyyyyz, g_xy_0_xy_xyyyzz, g_xy_0_xy_xyyzzz, g_xy_0_xy_xyzzzz, g_xy_0_xy_xzzzzz, g_xy_0_xy_yyyyyy, g_xy_0_xy_yyyyyz, g_xy_0_xy_yyyyzz, g_xy_0_xy_yyyzzz, g_xy_0_xy_yyzzzz, g_xy_0_xy_yzzzzz, g_xy_0_xy_zzzzzz, g_xy_0_y_xxxxxx, g_xy_0_y_xxxxxxx, g_xy_0_y_xxxxxxy, g_xy_0_y_xxxxxxz, g_xy_0_y_xxxxxy, g_xy_0_y_xxxxxyy, g_xy_0_y_xxxxxyz, g_xy_0_y_xxxxxz, g_xy_0_y_xxxxxzz, g_xy_0_y_xxxxyy, g_xy_0_y_xxxxyyy, g_xy_0_y_xxxxyyz, g_xy_0_y_xxxxyz, g_xy_0_y_xxxxyzz, g_xy_0_y_xxxxzz, g_xy_0_y_xxxxzzz, g_xy_0_y_xxxyyy, g_xy_0_y_xxxyyyy, g_xy_0_y_xxxyyyz, g_xy_0_y_xxxyyz, g_xy_0_y_xxxyyzz, g_xy_0_y_xxxyzz, g_xy_0_y_xxxyzzz, g_xy_0_y_xxxzzz, g_xy_0_y_xxxzzzz, g_xy_0_y_xxyyyy, g_xy_0_y_xxyyyyy, g_xy_0_y_xxyyyyz, g_xy_0_y_xxyyyz, g_xy_0_y_xxyyyzz, g_xy_0_y_xxyyzz, g_xy_0_y_xxyyzzz, g_xy_0_y_xxyzzz, g_xy_0_y_xxyzzzz, g_xy_0_y_xxzzzz, g_xy_0_y_xxzzzzz, g_xy_0_y_xyyyyy, g_xy_0_y_xyyyyyy, g_xy_0_y_xyyyyyz, g_xy_0_y_xyyyyz, g_xy_0_y_xyyyyzz, g_xy_0_y_xyyyzz, g_xy_0_y_xyyyzzz, g_xy_0_y_xyyzzz, g_xy_0_y_xyyzzzz, g_xy_0_y_xyzzzz, g_xy_0_y_xyzzzzz, g_xy_0_y_xzzzzz, g_xy_0_y_xzzzzzz, g_xy_0_y_yyyyyy, g_xy_0_y_yyyyyz, g_xy_0_y_yyyyzz, g_xy_0_y_yyyzzz, g_xy_0_y_yyzzzz, g_xy_0_y_yzzzzz, g_xy_0_y_zzzzzz, g_y_0_y_xxxxxx, g_y_0_y_xxxxxy, g_y_0_y_xxxxxz, g_y_0_y_xxxxyy, g_y_0_y_xxxxyz, g_y_0_y_xxxxzz, g_y_0_y_xxxyyy, g_y_0_y_xxxyyz, g_y_0_y_xxxyzz, g_y_0_y_xxxzzz, g_y_0_y_xxyyyy, g_y_0_y_xxyyyz, g_y_0_y_xxyyzz, g_y_0_y_xxyzzz, g_y_0_y_xxzzzz, g_y_0_y_xyyyyy, g_y_0_y_xyyyyz, g_y_0_y_xyyyzz, g_y_0_y_xyyzzz, g_y_0_y_xyzzzz, g_y_0_y_xzzzzz, g_y_0_y_yyyyyy, g_y_0_y_yyyyyz, g_y_0_y_yyyyzz, g_y_0_y_yyyzzz, g_y_0_y_yyzzzz, g_y_0_y_yzzzzz, g_y_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xy_xxxxxx[k] = -g_y_0_y_xxxxxx[k] - g_xy_0_y_xxxxxx[k] * ab_x + g_xy_0_y_xxxxxxx[k];

                g_xy_0_xy_xxxxxy[k] = -g_y_0_y_xxxxxy[k] - g_xy_0_y_xxxxxy[k] * ab_x + g_xy_0_y_xxxxxxy[k];

                g_xy_0_xy_xxxxxz[k] = -g_y_0_y_xxxxxz[k] - g_xy_0_y_xxxxxz[k] * ab_x + g_xy_0_y_xxxxxxz[k];

                g_xy_0_xy_xxxxyy[k] = -g_y_0_y_xxxxyy[k] - g_xy_0_y_xxxxyy[k] * ab_x + g_xy_0_y_xxxxxyy[k];

                g_xy_0_xy_xxxxyz[k] = -g_y_0_y_xxxxyz[k] - g_xy_0_y_xxxxyz[k] * ab_x + g_xy_0_y_xxxxxyz[k];

                g_xy_0_xy_xxxxzz[k] = -g_y_0_y_xxxxzz[k] - g_xy_0_y_xxxxzz[k] * ab_x + g_xy_0_y_xxxxxzz[k];

                g_xy_0_xy_xxxyyy[k] = -g_y_0_y_xxxyyy[k] - g_xy_0_y_xxxyyy[k] * ab_x + g_xy_0_y_xxxxyyy[k];

                g_xy_0_xy_xxxyyz[k] = -g_y_0_y_xxxyyz[k] - g_xy_0_y_xxxyyz[k] * ab_x + g_xy_0_y_xxxxyyz[k];

                g_xy_0_xy_xxxyzz[k] = -g_y_0_y_xxxyzz[k] - g_xy_0_y_xxxyzz[k] * ab_x + g_xy_0_y_xxxxyzz[k];

                g_xy_0_xy_xxxzzz[k] = -g_y_0_y_xxxzzz[k] - g_xy_0_y_xxxzzz[k] * ab_x + g_xy_0_y_xxxxzzz[k];

                g_xy_0_xy_xxyyyy[k] = -g_y_0_y_xxyyyy[k] - g_xy_0_y_xxyyyy[k] * ab_x + g_xy_0_y_xxxyyyy[k];

                g_xy_0_xy_xxyyyz[k] = -g_y_0_y_xxyyyz[k] - g_xy_0_y_xxyyyz[k] * ab_x + g_xy_0_y_xxxyyyz[k];

                g_xy_0_xy_xxyyzz[k] = -g_y_0_y_xxyyzz[k] - g_xy_0_y_xxyyzz[k] * ab_x + g_xy_0_y_xxxyyzz[k];

                g_xy_0_xy_xxyzzz[k] = -g_y_0_y_xxyzzz[k] - g_xy_0_y_xxyzzz[k] * ab_x + g_xy_0_y_xxxyzzz[k];

                g_xy_0_xy_xxzzzz[k] = -g_y_0_y_xxzzzz[k] - g_xy_0_y_xxzzzz[k] * ab_x + g_xy_0_y_xxxzzzz[k];

                g_xy_0_xy_xyyyyy[k] = -g_y_0_y_xyyyyy[k] - g_xy_0_y_xyyyyy[k] * ab_x + g_xy_0_y_xxyyyyy[k];

                g_xy_0_xy_xyyyyz[k] = -g_y_0_y_xyyyyz[k] - g_xy_0_y_xyyyyz[k] * ab_x + g_xy_0_y_xxyyyyz[k];

                g_xy_0_xy_xyyyzz[k] = -g_y_0_y_xyyyzz[k] - g_xy_0_y_xyyyzz[k] * ab_x + g_xy_0_y_xxyyyzz[k];

                g_xy_0_xy_xyyzzz[k] = -g_y_0_y_xyyzzz[k] - g_xy_0_y_xyyzzz[k] * ab_x + g_xy_0_y_xxyyzzz[k];

                g_xy_0_xy_xyzzzz[k] = -g_y_0_y_xyzzzz[k] - g_xy_0_y_xyzzzz[k] * ab_x + g_xy_0_y_xxyzzzz[k];

                g_xy_0_xy_xzzzzz[k] = -g_y_0_y_xzzzzz[k] - g_xy_0_y_xzzzzz[k] * ab_x + g_xy_0_y_xxzzzzz[k];

                g_xy_0_xy_yyyyyy[k] = -g_y_0_y_yyyyyy[k] - g_xy_0_y_yyyyyy[k] * ab_x + g_xy_0_y_xyyyyyy[k];

                g_xy_0_xy_yyyyyz[k] = -g_y_0_y_yyyyyz[k] - g_xy_0_y_yyyyyz[k] * ab_x + g_xy_0_y_xyyyyyz[k];

                g_xy_0_xy_yyyyzz[k] = -g_y_0_y_yyyyzz[k] - g_xy_0_y_yyyyzz[k] * ab_x + g_xy_0_y_xyyyyzz[k];

                g_xy_0_xy_yyyzzz[k] = -g_y_0_y_yyyzzz[k] - g_xy_0_y_yyyzzz[k] * ab_x + g_xy_0_y_xyyyzzz[k];

                g_xy_0_xy_yyzzzz[k] = -g_y_0_y_yyzzzz[k] - g_xy_0_y_yyzzzz[k] * ab_x + g_xy_0_y_xyyzzzz[k];

                g_xy_0_xy_yzzzzz[k] = -g_y_0_y_yzzzzz[k] - g_xy_0_y_yzzzzz[k] * ab_x + g_xy_0_y_xyzzzzz[k];

                g_xy_0_xy_zzzzzz[k] = -g_y_0_y_zzzzzz[k] - g_xy_0_y_zzzzzz[k] * ab_x + g_xy_0_y_xzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xz_xxxxxx = cbuffer.data(di_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_xz_xxxxxy = cbuffer.data(di_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_xz_xxxxxz = cbuffer.data(di_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_xz_xxxxyy = cbuffer.data(di_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_xz_xxxxyz = cbuffer.data(di_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_xz_xxxxzz = cbuffer.data(di_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_xz_xxxyyy = cbuffer.data(di_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_xz_xxxyyz = cbuffer.data(di_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_xz_xxxyzz = cbuffer.data(di_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_xz_xxxzzz = cbuffer.data(di_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_xz_xxyyyy = cbuffer.data(di_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_xz_xxyyyz = cbuffer.data(di_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_xz_xxyyzz = cbuffer.data(di_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_xz_xxyzzz = cbuffer.data(di_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_xz_xxzzzz = cbuffer.data(di_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_xz_xyyyyy = cbuffer.data(di_geom_20_off + 239 * ccomps * dcomps);

            auto g_xy_0_xz_xyyyyz = cbuffer.data(di_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_xz_xyyyzz = cbuffer.data(di_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_xz_xyyzzz = cbuffer.data(di_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_xz_xyzzzz = cbuffer.data(di_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_xz_xzzzzz = cbuffer.data(di_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_xz_yyyyyy = cbuffer.data(di_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_xz_yyyyyz = cbuffer.data(di_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_xz_yyyyzz = cbuffer.data(di_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_xz_yyyzzz = cbuffer.data(di_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_xz_yyzzzz = cbuffer.data(di_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_xz_yzzzzz = cbuffer.data(di_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_xz_zzzzzz = cbuffer.data(di_geom_20_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_x_xxxxxx, g_xy_0_x_xxxxxxz, g_xy_0_x_xxxxxy, g_xy_0_x_xxxxxyz, g_xy_0_x_xxxxxz, g_xy_0_x_xxxxxzz, g_xy_0_x_xxxxyy, g_xy_0_x_xxxxyyz, g_xy_0_x_xxxxyz, g_xy_0_x_xxxxyzz, g_xy_0_x_xxxxzz, g_xy_0_x_xxxxzzz, g_xy_0_x_xxxyyy, g_xy_0_x_xxxyyyz, g_xy_0_x_xxxyyz, g_xy_0_x_xxxyyzz, g_xy_0_x_xxxyzz, g_xy_0_x_xxxyzzz, g_xy_0_x_xxxzzz, g_xy_0_x_xxxzzzz, g_xy_0_x_xxyyyy, g_xy_0_x_xxyyyyz, g_xy_0_x_xxyyyz, g_xy_0_x_xxyyyzz, g_xy_0_x_xxyyzz, g_xy_0_x_xxyyzzz, g_xy_0_x_xxyzzz, g_xy_0_x_xxyzzzz, g_xy_0_x_xxzzzz, g_xy_0_x_xxzzzzz, g_xy_0_x_xyyyyy, g_xy_0_x_xyyyyyz, g_xy_0_x_xyyyyz, g_xy_0_x_xyyyyzz, g_xy_0_x_xyyyzz, g_xy_0_x_xyyyzzz, g_xy_0_x_xyyzzz, g_xy_0_x_xyyzzzz, g_xy_0_x_xyzzzz, g_xy_0_x_xyzzzzz, g_xy_0_x_xzzzzz, g_xy_0_x_xzzzzzz, g_xy_0_x_yyyyyy, g_xy_0_x_yyyyyyz, g_xy_0_x_yyyyyz, g_xy_0_x_yyyyyzz, g_xy_0_x_yyyyzz, g_xy_0_x_yyyyzzz, g_xy_0_x_yyyzzz, g_xy_0_x_yyyzzzz, g_xy_0_x_yyzzzz, g_xy_0_x_yyzzzzz, g_xy_0_x_yzzzzz, g_xy_0_x_yzzzzzz, g_xy_0_x_zzzzzz, g_xy_0_x_zzzzzzz, g_xy_0_xz_xxxxxx, g_xy_0_xz_xxxxxy, g_xy_0_xz_xxxxxz, g_xy_0_xz_xxxxyy, g_xy_0_xz_xxxxyz, g_xy_0_xz_xxxxzz, g_xy_0_xz_xxxyyy, g_xy_0_xz_xxxyyz, g_xy_0_xz_xxxyzz, g_xy_0_xz_xxxzzz, g_xy_0_xz_xxyyyy, g_xy_0_xz_xxyyyz, g_xy_0_xz_xxyyzz, g_xy_0_xz_xxyzzz, g_xy_0_xz_xxzzzz, g_xy_0_xz_xyyyyy, g_xy_0_xz_xyyyyz, g_xy_0_xz_xyyyzz, g_xy_0_xz_xyyzzz, g_xy_0_xz_xyzzzz, g_xy_0_xz_xzzzzz, g_xy_0_xz_yyyyyy, g_xy_0_xz_yyyyyz, g_xy_0_xz_yyyyzz, g_xy_0_xz_yyyzzz, g_xy_0_xz_yyzzzz, g_xy_0_xz_yzzzzz, g_xy_0_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xz_xxxxxx[k] = -g_xy_0_x_xxxxxx[k] * ab_z + g_xy_0_x_xxxxxxz[k];

                g_xy_0_xz_xxxxxy[k] = -g_xy_0_x_xxxxxy[k] * ab_z + g_xy_0_x_xxxxxyz[k];

                g_xy_0_xz_xxxxxz[k] = -g_xy_0_x_xxxxxz[k] * ab_z + g_xy_0_x_xxxxxzz[k];

                g_xy_0_xz_xxxxyy[k] = -g_xy_0_x_xxxxyy[k] * ab_z + g_xy_0_x_xxxxyyz[k];

                g_xy_0_xz_xxxxyz[k] = -g_xy_0_x_xxxxyz[k] * ab_z + g_xy_0_x_xxxxyzz[k];

                g_xy_0_xz_xxxxzz[k] = -g_xy_0_x_xxxxzz[k] * ab_z + g_xy_0_x_xxxxzzz[k];

                g_xy_0_xz_xxxyyy[k] = -g_xy_0_x_xxxyyy[k] * ab_z + g_xy_0_x_xxxyyyz[k];

                g_xy_0_xz_xxxyyz[k] = -g_xy_0_x_xxxyyz[k] * ab_z + g_xy_0_x_xxxyyzz[k];

                g_xy_0_xz_xxxyzz[k] = -g_xy_0_x_xxxyzz[k] * ab_z + g_xy_0_x_xxxyzzz[k];

                g_xy_0_xz_xxxzzz[k] = -g_xy_0_x_xxxzzz[k] * ab_z + g_xy_0_x_xxxzzzz[k];

                g_xy_0_xz_xxyyyy[k] = -g_xy_0_x_xxyyyy[k] * ab_z + g_xy_0_x_xxyyyyz[k];

                g_xy_0_xz_xxyyyz[k] = -g_xy_0_x_xxyyyz[k] * ab_z + g_xy_0_x_xxyyyzz[k];

                g_xy_0_xz_xxyyzz[k] = -g_xy_0_x_xxyyzz[k] * ab_z + g_xy_0_x_xxyyzzz[k];

                g_xy_0_xz_xxyzzz[k] = -g_xy_0_x_xxyzzz[k] * ab_z + g_xy_0_x_xxyzzzz[k];

                g_xy_0_xz_xxzzzz[k] = -g_xy_0_x_xxzzzz[k] * ab_z + g_xy_0_x_xxzzzzz[k];

                g_xy_0_xz_xyyyyy[k] = -g_xy_0_x_xyyyyy[k] * ab_z + g_xy_0_x_xyyyyyz[k];

                g_xy_0_xz_xyyyyz[k] = -g_xy_0_x_xyyyyz[k] * ab_z + g_xy_0_x_xyyyyzz[k];

                g_xy_0_xz_xyyyzz[k] = -g_xy_0_x_xyyyzz[k] * ab_z + g_xy_0_x_xyyyzzz[k];

                g_xy_0_xz_xyyzzz[k] = -g_xy_0_x_xyyzzz[k] * ab_z + g_xy_0_x_xyyzzzz[k];

                g_xy_0_xz_xyzzzz[k] = -g_xy_0_x_xyzzzz[k] * ab_z + g_xy_0_x_xyzzzzz[k];

                g_xy_0_xz_xzzzzz[k] = -g_xy_0_x_xzzzzz[k] * ab_z + g_xy_0_x_xzzzzzz[k];

                g_xy_0_xz_yyyyyy[k] = -g_xy_0_x_yyyyyy[k] * ab_z + g_xy_0_x_yyyyyyz[k];

                g_xy_0_xz_yyyyyz[k] = -g_xy_0_x_yyyyyz[k] * ab_z + g_xy_0_x_yyyyyzz[k];

                g_xy_0_xz_yyyyzz[k] = -g_xy_0_x_yyyyzz[k] * ab_z + g_xy_0_x_yyyyzzz[k];

                g_xy_0_xz_yyyzzz[k] = -g_xy_0_x_yyyzzz[k] * ab_z + g_xy_0_x_yyyzzzz[k];

                g_xy_0_xz_yyzzzz[k] = -g_xy_0_x_yyzzzz[k] * ab_z + g_xy_0_x_yyzzzzz[k];

                g_xy_0_xz_yzzzzz[k] = -g_xy_0_x_yzzzzz[k] * ab_z + g_xy_0_x_yzzzzzz[k];

                g_xy_0_xz_zzzzzz[k] = -g_xy_0_x_zzzzzz[k] * ab_z + g_xy_0_x_zzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yy_xxxxxx = cbuffer.data(di_geom_20_off + 252 * ccomps * dcomps);

            auto g_xy_0_yy_xxxxxy = cbuffer.data(di_geom_20_off + 253 * ccomps * dcomps);

            auto g_xy_0_yy_xxxxxz = cbuffer.data(di_geom_20_off + 254 * ccomps * dcomps);

            auto g_xy_0_yy_xxxxyy = cbuffer.data(di_geom_20_off + 255 * ccomps * dcomps);

            auto g_xy_0_yy_xxxxyz = cbuffer.data(di_geom_20_off + 256 * ccomps * dcomps);

            auto g_xy_0_yy_xxxxzz = cbuffer.data(di_geom_20_off + 257 * ccomps * dcomps);

            auto g_xy_0_yy_xxxyyy = cbuffer.data(di_geom_20_off + 258 * ccomps * dcomps);

            auto g_xy_0_yy_xxxyyz = cbuffer.data(di_geom_20_off + 259 * ccomps * dcomps);

            auto g_xy_0_yy_xxxyzz = cbuffer.data(di_geom_20_off + 260 * ccomps * dcomps);

            auto g_xy_0_yy_xxxzzz = cbuffer.data(di_geom_20_off + 261 * ccomps * dcomps);

            auto g_xy_0_yy_xxyyyy = cbuffer.data(di_geom_20_off + 262 * ccomps * dcomps);

            auto g_xy_0_yy_xxyyyz = cbuffer.data(di_geom_20_off + 263 * ccomps * dcomps);

            auto g_xy_0_yy_xxyyzz = cbuffer.data(di_geom_20_off + 264 * ccomps * dcomps);

            auto g_xy_0_yy_xxyzzz = cbuffer.data(di_geom_20_off + 265 * ccomps * dcomps);

            auto g_xy_0_yy_xxzzzz = cbuffer.data(di_geom_20_off + 266 * ccomps * dcomps);

            auto g_xy_0_yy_xyyyyy = cbuffer.data(di_geom_20_off + 267 * ccomps * dcomps);

            auto g_xy_0_yy_xyyyyz = cbuffer.data(di_geom_20_off + 268 * ccomps * dcomps);

            auto g_xy_0_yy_xyyyzz = cbuffer.data(di_geom_20_off + 269 * ccomps * dcomps);

            auto g_xy_0_yy_xyyzzz = cbuffer.data(di_geom_20_off + 270 * ccomps * dcomps);

            auto g_xy_0_yy_xyzzzz = cbuffer.data(di_geom_20_off + 271 * ccomps * dcomps);

            auto g_xy_0_yy_xzzzzz = cbuffer.data(di_geom_20_off + 272 * ccomps * dcomps);

            auto g_xy_0_yy_yyyyyy = cbuffer.data(di_geom_20_off + 273 * ccomps * dcomps);

            auto g_xy_0_yy_yyyyyz = cbuffer.data(di_geom_20_off + 274 * ccomps * dcomps);

            auto g_xy_0_yy_yyyyzz = cbuffer.data(di_geom_20_off + 275 * ccomps * dcomps);

            auto g_xy_0_yy_yyyzzz = cbuffer.data(di_geom_20_off + 276 * ccomps * dcomps);

            auto g_xy_0_yy_yyzzzz = cbuffer.data(di_geom_20_off + 277 * ccomps * dcomps);

            auto g_xy_0_yy_yzzzzz = cbuffer.data(di_geom_20_off + 278 * ccomps * dcomps);

            auto g_xy_0_yy_zzzzzz = cbuffer.data(di_geom_20_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxxxxx, g_x_0_y_xxxxxy, g_x_0_y_xxxxxz, g_x_0_y_xxxxyy, g_x_0_y_xxxxyz, g_x_0_y_xxxxzz, g_x_0_y_xxxyyy, g_x_0_y_xxxyyz, g_x_0_y_xxxyzz, g_x_0_y_xxxzzz, g_x_0_y_xxyyyy, g_x_0_y_xxyyyz, g_x_0_y_xxyyzz, g_x_0_y_xxyzzz, g_x_0_y_xxzzzz, g_x_0_y_xyyyyy, g_x_0_y_xyyyyz, g_x_0_y_xyyyzz, g_x_0_y_xyyzzz, g_x_0_y_xyzzzz, g_x_0_y_xzzzzz, g_x_0_y_yyyyyy, g_x_0_y_yyyyyz, g_x_0_y_yyyyzz, g_x_0_y_yyyzzz, g_x_0_y_yyzzzz, g_x_0_y_yzzzzz, g_x_0_y_zzzzzz, g_xy_0_y_xxxxxx, g_xy_0_y_xxxxxxy, g_xy_0_y_xxxxxy, g_xy_0_y_xxxxxyy, g_xy_0_y_xxxxxyz, g_xy_0_y_xxxxxz, g_xy_0_y_xxxxyy, g_xy_0_y_xxxxyyy, g_xy_0_y_xxxxyyz, g_xy_0_y_xxxxyz, g_xy_0_y_xxxxyzz, g_xy_0_y_xxxxzz, g_xy_0_y_xxxyyy, g_xy_0_y_xxxyyyy, g_xy_0_y_xxxyyyz, g_xy_0_y_xxxyyz, g_xy_0_y_xxxyyzz, g_xy_0_y_xxxyzz, g_xy_0_y_xxxyzzz, g_xy_0_y_xxxzzz, g_xy_0_y_xxyyyy, g_xy_0_y_xxyyyyy, g_xy_0_y_xxyyyyz, g_xy_0_y_xxyyyz, g_xy_0_y_xxyyyzz, g_xy_0_y_xxyyzz, g_xy_0_y_xxyyzzz, g_xy_0_y_xxyzzz, g_xy_0_y_xxyzzzz, g_xy_0_y_xxzzzz, g_xy_0_y_xyyyyy, g_xy_0_y_xyyyyyy, g_xy_0_y_xyyyyyz, g_xy_0_y_xyyyyz, g_xy_0_y_xyyyyzz, g_xy_0_y_xyyyzz, g_xy_0_y_xyyyzzz, g_xy_0_y_xyyzzz, g_xy_0_y_xyyzzzz, g_xy_0_y_xyzzzz, g_xy_0_y_xyzzzzz, g_xy_0_y_xzzzzz, g_xy_0_y_yyyyyy, g_xy_0_y_yyyyyyy, g_xy_0_y_yyyyyyz, g_xy_0_y_yyyyyz, g_xy_0_y_yyyyyzz, g_xy_0_y_yyyyzz, g_xy_0_y_yyyyzzz, g_xy_0_y_yyyzzz, g_xy_0_y_yyyzzzz, g_xy_0_y_yyzzzz, g_xy_0_y_yyzzzzz, g_xy_0_y_yzzzzz, g_xy_0_y_yzzzzzz, g_xy_0_y_zzzzzz, g_xy_0_yy_xxxxxx, g_xy_0_yy_xxxxxy, g_xy_0_yy_xxxxxz, g_xy_0_yy_xxxxyy, g_xy_0_yy_xxxxyz, g_xy_0_yy_xxxxzz, g_xy_0_yy_xxxyyy, g_xy_0_yy_xxxyyz, g_xy_0_yy_xxxyzz, g_xy_0_yy_xxxzzz, g_xy_0_yy_xxyyyy, g_xy_0_yy_xxyyyz, g_xy_0_yy_xxyyzz, g_xy_0_yy_xxyzzz, g_xy_0_yy_xxzzzz, g_xy_0_yy_xyyyyy, g_xy_0_yy_xyyyyz, g_xy_0_yy_xyyyzz, g_xy_0_yy_xyyzzz, g_xy_0_yy_xyzzzz, g_xy_0_yy_xzzzzz, g_xy_0_yy_yyyyyy, g_xy_0_yy_yyyyyz, g_xy_0_yy_yyyyzz, g_xy_0_yy_yyyzzz, g_xy_0_yy_yyzzzz, g_xy_0_yy_yzzzzz, g_xy_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yy_xxxxxx[k] = -g_x_0_y_xxxxxx[k] - g_xy_0_y_xxxxxx[k] * ab_y + g_xy_0_y_xxxxxxy[k];

                g_xy_0_yy_xxxxxy[k] = -g_x_0_y_xxxxxy[k] - g_xy_0_y_xxxxxy[k] * ab_y + g_xy_0_y_xxxxxyy[k];

                g_xy_0_yy_xxxxxz[k] = -g_x_0_y_xxxxxz[k] - g_xy_0_y_xxxxxz[k] * ab_y + g_xy_0_y_xxxxxyz[k];

                g_xy_0_yy_xxxxyy[k] = -g_x_0_y_xxxxyy[k] - g_xy_0_y_xxxxyy[k] * ab_y + g_xy_0_y_xxxxyyy[k];

                g_xy_0_yy_xxxxyz[k] = -g_x_0_y_xxxxyz[k] - g_xy_0_y_xxxxyz[k] * ab_y + g_xy_0_y_xxxxyyz[k];

                g_xy_0_yy_xxxxzz[k] = -g_x_0_y_xxxxzz[k] - g_xy_0_y_xxxxzz[k] * ab_y + g_xy_0_y_xxxxyzz[k];

                g_xy_0_yy_xxxyyy[k] = -g_x_0_y_xxxyyy[k] - g_xy_0_y_xxxyyy[k] * ab_y + g_xy_0_y_xxxyyyy[k];

                g_xy_0_yy_xxxyyz[k] = -g_x_0_y_xxxyyz[k] - g_xy_0_y_xxxyyz[k] * ab_y + g_xy_0_y_xxxyyyz[k];

                g_xy_0_yy_xxxyzz[k] = -g_x_0_y_xxxyzz[k] - g_xy_0_y_xxxyzz[k] * ab_y + g_xy_0_y_xxxyyzz[k];

                g_xy_0_yy_xxxzzz[k] = -g_x_0_y_xxxzzz[k] - g_xy_0_y_xxxzzz[k] * ab_y + g_xy_0_y_xxxyzzz[k];

                g_xy_0_yy_xxyyyy[k] = -g_x_0_y_xxyyyy[k] - g_xy_0_y_xxyyyy[k] * ab_y + g_xy_0_y_xxyyyyy[k];

                g_xy_0_yy_xxyyyz[k] = -g_x_0_y_xxyyyz[k] - g_xy_0_y_xxyyyz[k] * ab_y + g_xy_0_y_xxyyyyz[k];

                g_xy_0_yy_xxyyzz[k] = -g_x_0_y_xxyyzz[k] - g_xy_0_y_xxyyzz[k] * ab_y + g_xy_0_y_xxyyyzz[k];

                g_xy_0_yy_xxyzzz[k] = -g_x_0_y_xxyzzz[k] - g_xy_0_y_xxyzzz[k] * ab_y + g_xy_0_y_xxyyzzz[k];

                g_xy_0_yy_xxzzzz[k] = -g_x_0_y_xxzzzz[k] - g_xy_0_y_xxzzzz[k] * ab_y + g_xy_0_y_xxyzzzz[k];

                g_xy_0_yy_xyyyyy[k] = -g_x_0_y_xyyyyy[k] - g_xy_0_y_xyyyyy[k] * ab_y + g_xy_0_y_xyyyyyy[k];

                g_xy_0_yy_xyyyyz[k] = -g_x_0_y_xyyyyz[k] - g_xy_0_y_xyyyyz[k] * ab_y + g_xy_0_y_xyyyyyz[k];

                g_xy_0_yy_xyyyzz[k] = -g_x_0_y_xyyyzz[k] - g_xy_0_y_xyyyzz[k] * ab_y + g_xy_0_y_xyyyyzz[k];

                g_xy_0_yy_xyyzzz[k] = -g_x_0_y_xyyzzz[k] - g_xy_0_y_xyyzzz[k] * ab_y + g_xy_0_y_xyyyzzz[k];

                g_xy_0_yy_xyzzzz[k] = -g_x_0_y_xyzzzz[k] - g_xy_0_y_xyzzzz[k] * ab_y + g_xy_0_y_xyyzzzz[k];

                g_xy_0_yy_xzzzzz[k] = -g_x_0_y_xzzzzz[k] - g_xy_0_y_xzzzzz[k] * ab_y + g_xy_0_y_xyzzzzz[k];

                g_xy_0_yy_yyyyyy[k] = -g_x_0_y_yyyyyy[k] - g_xy_0_y_yyyyyy[k] * ab_y + g_xy_0_y_yyyyyyy[k];

                g_xy_0_yy_yyyyyz[k] = -g_x_0_y_yyyyyz[k] - g_xy_0_y_yyyyyz[k] * ab_y + g_xy_0_y_yyyyyyz[k];

                g_xy_0_yy_yyyyzz[k] = -g_x_0_y_yyyyzz[k] - g_xy_0_y_yyyyzz[k] * ab_y + g_xy_0_y_yyyyyzz[k];

                g_xy_0_yy_yyyzzz[k] = -g_x_0_y_yyyzzz[k] - g_xy_0_y_yyyzzz[k] * ab_y + g_xy_0_y_yyyyzzz[k];

                g_xy_0_yy_yyzzzz[k] = -g_x_0_y_yyzzzz[k] - g_xy_0_y_yyzzzz[k] * ab_y + g_xy_0_y_yyyzzzz[k];

                g_xy_0_yy_yzzzzz[k] = -g_x_0_y_yzzzzz[k] - g_xy_0_y_yzzzzz[k] * ab_y + g_xy_0_y_yyzzzzz[k];

                g_xy_0_yy_zzzzzz[k] = -g_x_0_y_zzzzzz[k] - g_xy_0_y_zzzzzz[k] * ab_y + g_xy_0_y_yzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yz_xxxxxx = cbuffer.data(di_geom_20_off + 280 * ccomps * dcomps);

            auto g_xy_0_yz_xxxxxy = cbuffer.data(di_geom_20_off + 281 * ccomps * dcomps);

            auto g_xy_0_yz_xxxxxz = cbuffer.data(di_geom_20_off + 282 * ccomps * dcomps);

            auto g_xy_0_yz_xxxxyy = cbuffer.data(di_geom_20_off + 283 * ccomps * dcomps);

            auto g_xy_0_yz_xxxxyz = cbuffer.data(di_geom_20_off + 284 * ccomps * dcomps);

            auto g_xy_0_yz_xxxxzz = cbuffer.data(di_geom_20_off + 285 * ccomps * dcomps);

            auto g_xy_0_yz_xxxyyy = cbuffer.data(di_geom_20_off + 286 * ccomps * dcomps);

            auto g_xy_0_yz_xxxyyz = cbuffer.data(di_geom_20_off + 287 * ccomps * dcomps);

            auto g_xy_0_yz_xxxyzz = cbuffer.data(di_geom_20_off + 288 * ccomps * dcomps);

            auto g_xy_0_yz_xxxzzz = cbuffer.data(di_geom_20_off + 289 * ccomps * dcomps);

            auto g_xy_0_yz_xxyyyy = cbuffer.data(di_geom_20_off + 290 * ccomps * dcomps);

            auto g_xy_0_yz_xxyyyz = cbuffer.data(di_geom_20_off + 291 * ccomps * dcomps);

            auto g_xy_0_yz_xxyyzz = cbuffer.data(di_geom_20_off + 292 * ccomps * dcomps);

            auto g_xy_0_yz_xxyzzz = cbuffer.data(di_geom_20_off + 293 * ccomps * dcomps);

            auto g_xy_0_yz_xxzzzz = cbuffer.data(di_geom_20_off + 294 * ccomps * dcomps);

            auto g_xy_0_yz_xyyyyy = cbuffer.data(di_geom_20_off + 295 * ccomps * dcomps);

            auto g_xy_0_yz_xyyyyz = cbuffer.data(di_geom_20_off + 296 * ccomps * dcomps);

            auto g_xy_0_yz_xyyyzz = cbuffer.data(di_geom_20_off + 297 * ccomps * dcomps);

            auto g_xy_0_yz_xyyzzz = cbuffer.data(di_geom_20_off + 298 * ccomps * dcomps);

            auto g_xy_0_yz_xyzzzz = cbuffer.data(di_geom_20_off + 299 * ccomps * dcomps);

            auto g_xy_0_yz_xzzzzz = cbuffer.data(di_geom_20_off + 300 * ccomps * dcomps);

            auto g_xy_0_yz_yyyyyy = cbuffer.data(di_geom_20_off + 301 * ccomps * dcomps);

            auto g_xy_0_yz_yyyyyz = cbuffer.data(di_geom_20_off + 302 * ccomps * dcomps);

            auto g_xy_0_yz_yyyyzz = cbuffer.data(di_geom_20_off + 303 * ccomps * dcomps);

            auto g_xy_0_yz_yyyzzz = cbuffer.data(di_geom_20_off + 304 * ccomps * dcomps);

            auto g_xy_0_yz_yyzzzz = cbuffer.data(di_geom_20_off + 305 * ccomps * dcomps);

            auto g_xy_0_yz_yzzzzz = cbuffer.data(di_geom_20_off + 306 * ccomps * dcomps);

            auto g_xy_0_yz_zzzzzz = cbuffer.data(di_geom_20_off + 307 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_y_xxxxxx, g_xy_0_y_xxxxxxz, g_xy_0_y_xxxxxy, g_xy_0_y_xxxxxyz, g_xy_0_y_xxxxxz, g_xy_0_y_xxxxxzz, g_xy_0_y_xxxxyy, g_xy_0_y_xxxxyyz, g_xy_0_y_xxxxyz, g_xy_0_y_xxxxyzz, g_xy_0_y_xxxxzz, g_xy_0_y_xxxxzzz, g_xy_0_y_xxxyyy, g_xy_0_y_xxxyyyz, g_xy_0_y_xxxyyz, g_xy_0_y_xxxyyzz, g_xy_0_y_xxxyzz, g_xy_0_y_xxxyzzz, g_xy_0_y_xxxzzz, g_xy_0_y_xxxzzzz, g_xy_0_y_xxyyyy, g_xy_0_y_xxyyyyz, g_xy_0_y_xxyyyz, g_xy_0_y_xxyyyzz, g_xy_0_y_xxyyzz, g_xy_0_y_xxyyzzz, g_xy_0_y_xxyzzz, g_xy_0_y_xxyzzzz, g_xy_0_y_xxzzzz, g_xy_0_y_xxzzzzz, g_xy_0_y_xyyyyy, g_xy_0_y_xyyyyyz, g_xy_0_y_xyyyyz, g_xy_0_y_xyyyyzz, g_xy_0_y_xyyyzz, g_xy_0_y_xyyyzzz, g_xy_0_y_xyyzzz, g_xy_0_y_xyyzzzz, g_xy_0_y_xyzzzz, g_xy_0_y_xyzzzzz, g_xy_0_y_xzzzzz, g_xy_0_y_xzzzzzz, g_xy_0_y_yyyyyy, g_xy_0_y_yyyyyyz, g_xy_0_y_yyyyyz, g_xy_0_y_yyyyyzz, g_xy_0_y_yyyyzz, g_xy_0_y_yyyyzzz, g_xy_0_y_yyyzzz, g_xy_0_y_yyyzzzz, g_xy_0_y_yyzzzz, g_xy_0_y_yyzzzzz, g_xy_0_y_yzzzzz, g_xy_0_y_yzzzzzz, g_xy_0_y_zzzzzz, g_xy_0_y_zzzzzzz, g_xy_0_yz_xxxxxx, g_xy_0_yz_xxxxxy, g_xy_0_yz_xxxxxz, g_xy_0_yz_xxxxyy, g_xy_0_yz_xxxxyz, g_xy_0_yz_xxxxzz, g_xy_0_yz_xxxyyy, g_xy_0_yz_xxxyyz, g_xy_0_yz_xxxyzz, g_xy_0_yz_xxxzzz, g_xy_0_yz_xxyyyy, g_xy_0_yz_xxyyyz, g_xy_0_yz_xxyyzz, g_xy_0_yz_xxyzzz, g_xy_0_yz_xxzzzz, g_xy_0_yz_xyyyyy, g_xy_0_yz_xyyyyz, g_xy_0_yz_xyyyzz, g_xy_0_yz_xyyzzz, g_xy_0_yz_xyzzzz, g_xy_0_yz_xzzzzz, g_xy_0_yz_yyyyyy, g_xy_0_yz_yyyyyz, g_xy_0_yz_yyyyzz, g_xy_0_yz_yyyzzz, g_xy_0_yz_yyzzzz, g_xy_0_yz_yzzzzz, g_xy_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yz_xxxxxx[k] = -g_xy_0_y_xxxxxx[k] * ab_z + g_xy_0_y_xxxxxxz[k];

                g_xy_0_yz_xxxxxy[k] = -g_xy_0_y_xxxxxy[k] * ab_z + g_xy_0_y_xxxxxyz[k];

                g_xy_0_yz_xxxxxz[k] = -g_xy_0_y_xxxxxz[k] * ab_z + g_xy_0_y_xxxxxzz[k];

                g_xy_0_yz_xxxxyy[k] = -g_xy_0_y_xxxxyy[k] * ab_z + g_xy_0_y_xxxxyyz[k];

                g_xy_0_yz_xxxxyz[k] = -g_xy_0_y_xxxxyz[k] * ab_z + g_xy_0_y_xxxxyzz[k];

                g_xy_0_yz_xxxxzz[k] = -g_xy_0_y_xxxxzz[k] * ab_z + g_xy_0_y_xxxxzzz[k];

                g_xy_0_yz_xxxyyy[k] = -g_xy_0_y_xxxyyy[k] * ab_z + g_xy_0_y_xxxyyyz[k];

                g_xy_0_yz_xxxyyz[k] = -g_xy_0_y_xxxyyz[k] * ab_z + g_xy_0_y_xxxyyzz[k];

                g_xy_0_yz_xxxyzz[k] = -g_xy_0_y_xxxyzz[k] * ab_z + g_xy_0_y_xxxyzzz[k];

                g_xy_0_yz_xxxzzz[k] = -g_xy_0_y_xxxzzz[k] * ab_z + g_xy_0_y_xxxzzzz[k];

                g_xy_0_yz_xxyyyy[k] = -g_xy_0_y_xxyyyy[k] * ab_z + g_xy_0_y_xxyyyyz[k];

                g_xy_0_yz_xxyyyz[k] = -g_xy_0_y_xxyyyz[k] * ab_z + g_xy_0_y_xxyyyzz[k];

                g_xy_0_yz_xxyyzz[k] = -g_xy_0_y_xxyyzz[k] * ab_z + g_xy_0_y_xxyyzzz[k];

                g_xy_0_yz_xxyzzz[k] = -g_xy_0_y_xxyzzz[k] * ab_z + g_xy_0_y_xxyzzzz[k];

                g_xy_0_yz_xxzzzz[k] = -g_xy_0_y_xxzzzz[k] * ab_z + g_xy_0_y_xxzzzzz[k];

                g_xy_0_yz_xyyyyy[k] = -g_xy_0_y_xyyyyy[k] * ab_z + g_xy_0_y_xyyyyyz[k];

                g_xy_0_yz_xyyyyz[k] = -g_xy_0_y_xyyyyz[k] * ab_z + g_xy_0_y_xyyyyzz[k];

                g_xy_0_yz_xyyyzz[k] = -g_xy_0_y_xyyyzz[k] * ab_z + g_xy_0_y_xyyyzzz[k];

                g_xy_0_yz_xyyzzz[k] = -g_xy_0_y_xyyzzz[k] * ab_z + g_xy_0_y_xyyzzzz[k];

                g_xy_0_yz_xyzzzz[k] = -g_xy_0_y_xyzzzz[k] * ab_z + g_xy_0_y_xyzzzzz[k];

                g_xy_0_yz_xzzzzz[k] = -g_xy_0_y_xzzzzz[k] * ab_z + g_xy_0_y_xzzzzzz[k];

                g_xy_0_yz_yyyyyy[k] = -g_xy_0_y_yyyyyy[k] * ab_z + g_xy_0_y_yyyyyyz[k];

                g_xy_0_yz_yyyyyz[k] = -g_xy_0_y_yyyyyz[k] * ab_z + g_xy_0_y_yyyyyzz[k];

                g_xy_0_yz_yyyyzz[k] = -g_xy_0_y_yyyyzz[k] * ab_z + g_xy_0_y_yyyyzzz[k];

                g_xy_0_yz_yyyzzz[k] = -g_xy_0_y_yyyzzz[k] * ab_z + g_xy_0_y_yyyzzzz[k];

                g_xy_0_yz_yyzzzz[k] = -g_xy_0_y_yyzzzz[k] * ab_z + g_xy_0_y_yyzzzzz[k];

                g_xy_0_yz_yzzzzz[k] = -g_xy_0_y_yzzzzz[k] * ab_z + g_xy_0_y_yzzzzzz[k];

                g_xy_0_yz_zzzzzz[k] = -g_xy_0_y_zzzzzz[k] * ab_z + g_xy_0_y_zzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zz_xxxxxx = cbuffer.data(di_geom_20_off + 308 * ccomps * dcomps);

            auto g_xy_0_zz_xxxxxy = cbuffer.data(di_geom_20_off + 309 * ccomps * dcomps);

            auto g_xy_0_zz_xxxxxz = cbuffer.data(di_geom_20_off + 310 * ccomps * dcomps);

            auto g_xy_0_zz_xxxxyy = cbuffer.data(di_geom_20_off + 311 * ccomps * dcomps);

            auto g_xy_0_zz_xxxxyz = cbuffer.data(di_geom_20_off + 312 * ccomps * dcomps);

            auto g_xy_0_zz_xxxxzz = cbuffer.data(di_geom_20_off + 313 * ccomps * dcomps);

            auto g_xy_0_zz_xxxyyy = cbuffer.data(di_geom_20_off + 314 * ccomps * dcomps);

            auto g_xy_0_zz_xxxyyz = cbuffer.data(di_geom_20_off + 315 * ccomps * dcomps);

            auto g_xy_0_zz_xxxyzz = cbuffer.data(di_geom_20_off + 316 * ccomps * dcomps);

            auto g_xy_0_zz_xxxzzz = cbuffer.data(di_geom_20_off + 317 * ccomps * dcomps);

            auto g_xy_0_zz_xxyyyy = cbuffer.data(di_geom_20_off + 318 * ccomps * dcomps);

            auto g_xy_0_zz_xxyyyz = cbuffer.data(di_geom_20_off + 319 * ccomps * dcomps);

            auto g_xy_0_zz_xxyyzz = cbuffer.data(di_geom_20_off + 320 * ccomps * dcomps);

            auto g_xy_0_zz_xxyzzz = cbuffer.data(di_geom_20_off + 321 * ccomps * dcomps);

            auto g_xy_0_zz_xxzzzz = cbuffer.data(di_geom_20_off + 322 * ccomps * dcomps);

            auto g_xy_0_zz_xyyyyy = cbuffer.data(di_geom_20_off + 323 * ccomps * dcomps);

            auto g_xy_0_zz_xyyyyz = cbuffer.data(di_geom_20_off + 324 * ccomps * dcomps);

            auto g_xy_0_zz_xyyyzz = cbuffer.data(di_geom_20_off + 325 * ccomps * dcomps);

            auto g_xy_0_zz_xyyzzz = cbuffer.data(di_geom_20_off + 326 * ccomps * dcomps);

            auto g_xy_0_zz_xyzzzz = cbuffer.data(di_geom_20_off + 327 * ccomps * dcomps);

            auto g_xy_0_zz_xzzzzz = cbuffer.data(di_geom_20_off + 328 * ccomps * dcomps);

            auto g_xy_0_zz_yyyyyy = cbuffer.data(di_geom_20_off + 329 * ccomps * dcomps);

            auto g_xy_0_zz_yyyyyz = cbuffer.data(di_geom_20_off + 330 * ccomps * dcomps);

            auto g_xy_0_zz_yyyyzz = cbuffer.data(di_geom_20_off + 331 * ccomps * dcomps);

            auto g_xy_0_zz_yyyzzz = cbuffer.data(di_geom_20_off + 332 * ccomps * dcomps);

            auto g_xy_0_zz_yyzzzz = cbuffer.data(di_geom_20_off + 333 * ccomps * dcomps);

            auto g_xy_0_zz_yzzzzz = cbuffer.data(di_geom_20_off + 334 * ccomps * dcomps);

            auto g_xy_0_zz_zzzzzz = cbuffer.data(di_geom_20_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_z_xxxxxx, g_xy_0_z_xxxxxxz, g_xy_0_z_xxxxxy, g_xy_0_z_xxxxxyz, g_xy_0_z_xxxxxz, g_xy_0_z_xxxxxzz, g_xy_0_z_xxxxyy, g_xy_0_z_xxxxyyz, g_xy_0_z_xxxxyz, g_xy_0_z_xxxxyzz, g_xy_0_z_xxxxzz, g_xy_0_z_xxxxzzz, g_xy_0_z_xxxyyy, g_xy_0_z_xxxyyyz, g_xy_0_z_xxxyyz, g_xy_0_z_xxxyyzz, g_xy_0_z_xxxyzz, g_xy_0_z_xxxyzzz, g_xy_0_z_xxxzzz, g_xy_0_z_xxxzzzz, g_xy_0_z_xxyyyy, g_xy_0_z_xxyyyyz, g_xy_0_z_xxyyyz, g_xy_0_z_xxyyyzz, g_xy_0_z_xxyyzz, g_xy_0_z_xxyyzzz, g_xy_0_z_xxyzzz, g_xy_0_z_xxyzzzz, g_xy_0_z_xxzzzz, g_xy_0_z_xxzzzzz, g_xy_0_z_xyyyyy, g_xy_0_z_xyyyyyz, g_xy_0_z_xyyyyz, g_xy_0_z_xyyyyzz, g_xy_0_z_xyyyzz, g_xy_0_z_xyyyzzz, g_xy_0_z_xyyzzz, g_xy_0_z_xyyzzzz, g_xy_0_z_xyzzzz, g_xy_0_z_xyzzzzz, g_xy_0_z_xzzzzz, g_xy_0_z_xzzzzzz, g_xy_0_z_yyyyyy, g_xy_0_z_yyyyyyz, g_xy_0_z_yyyyyz, g_xy_0_z_yyyyyzz, g_xy_0_z_yyyyzz, g_xy_0_z_yyyyzzz, g_xy_0_z_yyyzzz, g_xy_0_z_yyyzzzz, g_xy_0_z_yyzzzz, g_xy_0_z_yyzzzzz, g_xy_0_z_yzzzzz, g_xy_0_z_yzzzzzz, g_xy_0_z_zzzzzz, g_xy_0_z_zzzzzzz, g_xy_0_zz_xxxxxx, g_xy_0_zz_xxxxxy, g_xy_0_zz_xxxxxz, g_xy_0_zz_xxxxyy, g_xy_0_zz_xxxxyz, g_xy_0_zz_xxxxzz, g_xy_0_zz_xxxyyy, g_xy_0_zz_xxxyyz, g_xy_0_zz_xxxyzz, g_xy_0_zz_xxxzzz, g_xy_0_zz_xxyyyy, g_xy_0_zz_xxyyyz, g_xy_0_zz_xxyyzz, g_xy_0_zz_xxyzzz, g_xy_0_zz_xxzzzz, g_xy_0_zz_xyyyyy, g_xy_0_zz_xyyyyz, g_xy_0_zz_xyyyzz, g_xy_0_zz_xyyzzz, g_xy_0_zz_xyzzzz, g_xy_0_zz_xzzzzz, g_xy_0_zz_yyyyyy, g_xy_0_zz_yyyyyz, g_xy_0_zz_yyyyzz, g_xy_0_zz_yyyzzz, g_xy_0_zz_yyzzzz, g_xy_0_zz_yzzzzz, g_xy_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zz_xxxxxx[k] = -g_xy_0_z_xxxxxx[k] * ab_z + g_xy_0_z_xxxxxxz[k];

                g_xy_0_zz_xxxxxy[k] = -g_xy_0_z_xxxxxy[k] * ab_z + g_xy_0_z_xxxxxyz[k];

                g_xy_0_zz_xxxxxz[k] = -g_xy_0_z_xxxxxz[k] * ab_z + g_xy_0_z_xxxxxzz[k];

                g_xy_0_zz_xxxxyy[k] = -g_xy_0_z_xxxxyy[k] * ab_z + g_xy_0_z_xxxxyyz[k];

                g_xy_0_zz_xxxxyz[k] = -g_xy_0_z_xxxxyz[k] * ab_z + g_xy_0_z_xxxxyzz[k];

                g_xy_0_zz_xxxxzz[k] = -g_xy_0_z_xxxxzz[k] * ab_z + g_xy_0_z_xxxxzzz[k];

                g_xy_0_zz_xxxyyy[k] = -g_xy_0_z_xxxyyy[k] * ab_z + g_xy_0_z_xxxyyyz[k];

                g_xy_0_zz_xxxyyz[k] = -g_xy_0_z_xxxyyz[k] * ab_z + g_xy_0_z_xxxyyzz[k];

                g_xy_0_zz_xxxyzz[k] = -g_xy_0_z_xxxyzz[k] * ab_z + g_xy_0_z_xxxyzzz[k];

                g_xy_0_zz_xxxzzz[k] = -g_xy_0_z_xxxzzz[k] * ab_z + g_xy_0_z_xxxzzzz[k];

                g_xy_0_zz_xxyyyy[k] = -g_xy_0_z_xxyyyy[k] * ab_z + g_xy_0_z_xxyyyyz[k];

                g_xy_0_zz_xxyyyz[k] = -g_xy_0_z_xxyyyz[k] * ab_z + g_xy_0_z_xxyyyzz[k];

                g_xy_0_zz_xxyyzz[k] = -g_xy_0_z_xxyyzz[k] * ab_z + g_xy_0_z_xxyyzzz[k];

                g_xy_0_zz_xxyzzz[k] = -g_xy_0_z_xxyzzz[k] * ab_z + g_xy_0_z_xxyzzzz[k];

                g_xy_0_zz_xxzzzz[k] = -g_xy_0_z_xxzzzz[k] * ab_z + g_xy_0_z_xxzzzzz[k];

                g_xy_0_zz_xyyyyy[k] = -g_xy_0_z_xyyyyy[k] * ab_z + g_xy_0_z_xyyyyyz[k];

                g_xy_0_zz_xyyyyz[k] = -g_xy_0_z_xyyyyz[k] * ab_z + g_xy_0_z_xyyyyzz[k];

                g_xy_0_zz_xyyyzz[k] = -g_xy_0_z_xyyyzz[k] * ab_z + g_xy_0_z_xyyyzzz[k];

                g_xy_0_zz_xyyzzz[k] = -g_xy_0_z_xyyzzz[k] * ab_z + g_xy_0_z_xyyzzzz[k];

                g_xy_0_zz_xyzzzz[k] = -g_xy_0_z_xyzzzz[k] * ab_z + g_xy_0_z_xyzzzzz[k];

                g_xy_0_zz_xzzzzz[k] = -g_xy_0_z_xzzzzz[k] * ab_z + g_xy_0_z_xzzzzzz[k];

                g_xy_0_zz_yyyyyy[k] = -g_xy_0_z_yyyyyy[k] * ab_z + g_xy_0_z_yyyyyyz[k];

                g_xy_0_zz_yyyyyz[k] = -g_xy_0_z_yyyyyz[k] * ab_z + g_xy_0_z_yyyyyzz[k];

                g_xy_0_zz_yyyyzz[k] = -g_xy_0_z_yyyyzz[k] * ab_z + g_xy_0_z_yyyyzzz[k];

                g_xy_0_zz_yyyzzz[k] = -g_xy_0_z_yyyzzz[k] * ab_z + g_xy_0_z_yyyzzzz[k];

                g_xy_0_zz_yyzzzz[k] = -g_xy_0_z_yyzzzz[k] * ab_z + g_xy_0_z_yyzzzzz[k];

                g_xy_0_zz_yzzzzz[k] = -g_xy_0_z_yzzzzz[k] * ab_z + g_xy_0_z_yzzzzzz[k];

                g_xy_0_zz_zzzzzz[k] = -g_xy_0_z_zzzzzz[k] * ab_z + g_xy_0_z_zzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xx_xxxxxx = cbuffer.data(di_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_xx_xxxxxy = cbuffer.data(di_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_xx_xxxxxz = cbuffer.data(di_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_xx_xxxxyy = cbuffer.data(di_geom_20_off + 339 * ccomps * dcomps);

            auto g_xz_0_xx_xxxxyz = cbuffer.data(di_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_xx_xxxxzz = cbuffer.data(di_geom_20_off + 341 * ccomps * dcomps);

            auto g_xz_0_xx_xxxyyy = cbuffer.data(di_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_xx_xxxyyz = cbuffer.data(di_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_xx_xxxyzz = cbuffer.data(di_geom_20_off + 344 * ccomps * dcomps);

            auto g_xz_0_xx_xxxzzz = cbuffer.data(di_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_xx_xxyyyy = cbuffer.data(di_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_xx_xxyyyz = cbuffer.data(di_geom_20_off + 347 * ccomps * dcomps);

            auto g_xz_0_xx_xxyyzz = cbuffer.data(di_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_xx_xxyzzz = cbuffer.data(di_geom_20_off + 349 * ccomps * dcomps);

            auto g_xz_0_xx_xxzzzz = cbuffer.data(di_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_xx_xyyyyy = cbuffer.data(di_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_xx_xyyyyz = cbuffer.data(di_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_xx_xyyyzz = cbuffer.data(di_geom_20_off + 353 * ccomps * dcomps);

            auto g_xz_0_xx_xyyzzz = cbuffer.data(di_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_xx_xyzzzz = cbuffer.data(di_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_xx_xzzzzz = cbuffer.data(di_geom_20_off + 356 * ccomps * dcomps);

            auto g_xz_0_xx_yyyyyy = cbuffer.data(di_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_xx_yyyyyz = cbuffer.data(di_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_xx_yyyyzz = cbuffer.data(di_geom_20_off + 359 * ccomps * dcomps);

            auto g_xz_0_xx_yyyzzz = cbuffer.data(di_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_xx_yyzzzz = cbuffer.data(di_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_xx_yzzzzz = cbuffer.data(di_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_xx_zzzzzz = cbuffer.data(di_geom_20_off + 363 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_x_xxxxxx, g_xz_0_x_xxxxxxx, g_xz_0_x_xxxxxxy, g_xz_0_x_xxxxxxz, g_xz_0_x_xxxxxy, g_xz_0_x_xxxxxyy, g_xz_0_x_xxxxxyz, g_xz_0_x_xxxxxz, g_xz_0_x_xxxxxzz, g_xz_0_x_xxxxyy, g_xz_0_x_xxxxyyy, g_xz_0_x_xxxxyyz, g_xz_0_x_xxxxyz, g_xz_0_x_xxxxyzz, g_xz_0_x_xxxxzz, g_xz_0_x_xxxxzzz, g_xz_0_x_xxxyyy, g_xz_0_x_xxxyyyy, g_xz_0_x_xxxyyyz, g_xz_0_x_xxxyyz, g_xz_0_x_xxxyyzz, g_xz_0_x_xxxyzz, g_xz_0_x_xxxyzzz, g_xz_0_x_xxxzzz, g_xz_0_x_xxxzzzz, g_xz_0_x_xxyyyy, g_xz_0_x_xxyyyyy, g_xz_0_x_xxyyyyz, g_xz_0_x_xxyyyz, g_xz_0_x_xxyyyzz, g_xz_0_x_xxyyzz, g_xz_0_x_xxyyzzz, g_xz_0_x_xxyzzz, g_xz_0_x_xxyzzzz, g_xz_0_x_xxzzzz, g_xz_0_x_xxzzzzz, g_xz_0_x_xyyyyy, g_xz_0_x_xyyyyyy, g_xz_0_x_xyyyyyz, g_xz_0_x_xyyyyz, g_xz_0_x_xyyyyzz, g_xz_0_x_xyyyzz, g_xz_0_x_xyyyzzz, g_xz_0_x_xyyzzz, g_xz_0_x_xyyzzzz, g_xz_0_x_xyzzzz, g_xz_0_x_xyzzzzz, g_xz_0_x_xzzzzz, g_xz_0_x_xzzzzzz, g_xz_0_x_yyyyyy, g_xz_0_x_yyyyyz, g_xz_0_x_yyyyzz, g_xz_0_x_yyyzzz, g_xz_0_x_yyzzzz, g_xz_0_x_yzzzzz, g_xz_0_x_zzzzzz, g_xz_0_xx_xxxxxx, g_xz_0_xx_xxxxxy, g_xz_0_xx_xxxxxz, g_xz_0_xx_xxxxyy, g_xz_0_xx_xxxxyz, g_xz_0_xx_xxxxzz, g_xz_0_xx_xxxyyy, g_xz_0_xx_xxxyyz, g_xz_0_xx_xxxyzz, g_xz_0_xx_xxxzzz, g_xz_0_xx_xxyyyy, g_xz_0_xx_xxyyyz, g_xz_0_xx_xxyyzz, g_xz_0_xx_xxyzzz, g_xz_0_xx_xxzzzz, g_xz_0_xx_xyyyyy, g_xz_0_xx_xyyyyz, g_xz_0_xx_xyyyzz, g_xz_0_xx_xyyzzz, g_xz_0_xx_xyzzzz, g_xz_0_xx_xzzzzz, g_xz_0_xx_yyyyyy, g_xz_0_xx_yyyyyz, g_xz_0_xx_yyyyzz, g_xz_0_xx_yyyzzz, g_xz_0_xx_yyzzzz, g_xz_0_xx_yzzzzz, g_xz_0_xx_zzzzzz, g_z_0_x_xxxxxx, g_z_0_x_xxxxxy, g_z_0_x_xxxxxz, g_z_0_x_xxxxyy, g_z_0_x_xxxxyz, g_z_0_x_xxxxzz, g_z_0_x_xxxyyy, g_z_0_x_xxxyyz, g_z_0_x_xxxyzz, g_z_0_x_xxxzzz, g_z_0_x_xxyyyy, g_z_0_x_xxyyyz, g_z_0_x_xxyyzz, g_z_0_x_xxyzzz, g_z_0_x_xxzzzz, g_z_0_x_xyyyyy, g_z_0_x_xyyyyz, g_z_0_x_xyyyzz, g_z_0_x_xyyzzz, g_z_0_x_xyzzzz, g_z_0_x_xzzzzz, g_z_0_x_yyyyyy, g_z_0_x_yyyyyz, g_z_0_x_yyyyzz, g_z_0_x_yyyzzz, g_z_0_x_yyzzzz, g_z_0_x_yzzzzz, g_z_0_x_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xx_xxxxxx[k] = -g_z_0_x_xxxxxx[k] - g_xz_0_x_xxxxxx[k] * ab_x + g_xz_0_x_xxxxxxx[k];

                g_xz_0_xx_xxxxxy[k] = -g_z_0_x_xxxxxy[k] - g_xz_0_x_xxxxxy[k] * ab_x + g_xz_0_x_xxxxxxy[k];

                g_xz_0_xx_xxxxxz[k] = -g_z_0_x_xxxxxz[k] - g_xz_0_x_xxxxxz[k] * ab_x + g_xz_0_x_xxxxxxz[k];

                g_xz_0_xx_xxxxyy[k] = -g_z_0_x_xxxxyy[k] - g_xz_0_x_xxxxyy[k] * ab_x + g_xz_0_x_xxxxxyy[k];

                g_xz_0_xx_xxxxyz[k] = -g_z_0_x_xxxxyz[k] - g_xz_0_x_xxxxyz[k] * ab_x + g_xz_0_x_xxxxxyz[k];

                g_xz_0_xx_xxxxzz[k] = -g_z_0_x_xxxxzz[k] - g_xz_0_x_xxxxzz[k] * ab_x + g_xz_0_x_xxxxxzz[k];

                g_xz_0_xx_xxxyyy[k] = -g_z_0_x_xxxyyy[k] - g_xz_0_x_xxxyyy[k] * ab_x + g_xz_0_x_xxxxyyy[k];

                g_xz_0_xx_xxxyyz[k] = -g_z_0_x_xxxyyz[k] - g_xz_0_x_xxxyyz[k] * ab_x + g_xz_0_x_xxxxyyz[k];

                g_xz_0_xx_xxxyzz[k] = -g_z_0_x_xxxyzz[k] - g_xz_0_x_xxxyzz[k] * ab_x + g_xz_0_x_xxxxyzz[k];

                g_xz_0_xx_xxxzzz[k] = -g_z_0_x_xxxzzz[k] - g_xz_0_x_xxxzzz[k] * ab_x + g_xz_0_x_xxxxzzz[k];

                g_xz_0_xx_xxyyyy[k] = -g_z_0_x_xxyyyy[k] - g_xz_0_x_xxyyyy[k] * ab_x + g_xz_0_x_xxxyyyy[k];

                g_xz_0_xx_xxyyyz[k] = -g_z_0_x_xxyyyz[k] - g_xz_0_x_xxyyyz[k] * ab_x + g_xz_0_x_xxxyyyz[k];

                g_xz_0_xx_xxyyzz[k] = -g_z_0_x_xxyyzz[k] - g_xz_0_x_xxyyzz[k] * ab_x + g_xz_0_x_xxxyyzz[k];

                g_xz_0_xx_xxyzzz[k] = -g_z_0_x_xxyzzz[k] - g_xz_0_x_xxyzzz[k] * ab_x + g_xz_0_x_xxxyzzz[k];

                g_xz_0_xx_xxzzzz[k] = -g_z_0_x_xxzzzz[k] - g_xz_0_x_xxzzzz[k] * ab_x + g_xz_0_x_xxxzzzz[k];

                g_xz_0_xx_xyyyyy[k] = -g_z_0_x_xyyyyy[k] - g_xz_0_x_xyyyyy[k] * ab_x + g_xz_0_x_xxyyyyy[k];

                g_xz_0_xx_xyyyyz[k] = -g_z_0_x_xyyyyz[k] - g_xz_0_x_xyyyyz[k] * ab_x + g_xz_0_x_xxyyyyz[k];

                g_xz_0_xx_xyyyzz[k] = -g_z_0_x_xyyyzz[k] - g_xz_0_x_xyyyzz[k] * ab_x + g_xz_0_x_xxyyyzz[k];

                g_xz_0_xx_xyyzzz[k] = -g_z_0_x_xyyzzz[k] - g_xz_0_x_xyyzzz[k] * ab_x + g_xz_0_x_xxyyzzz[k];

                g_xz_0_xx_xyzzzz[k] = -g_z_0_x_xyzzzz[k] - g_xz_0_x_xyzzzz[k] * ab_x + g_xz_0_x_xxyzzzz[k];

                g_xz_0_xx_xzzzzz[k] = -g_z_0_x_xzzzzz[k] - g_xz_0_x_xzzzzz[k] * ab_x + g_xz_0_x_xxzzzzz[k];

                g_xz_0_xx_yyyyyy[k] = -g_z_0_x_yyyyyy[k] - g_xz_0_x_yyyyyy[k] * ab_x + g_xz_0_x_xyyyyyy[k];

                g_xz_0_xx_yyyyyz[k] = -g_z_0_x_yyyyyz[k] - g_xz_0_x_yyyyyz[k] * ab_x + g_xz_0_x_xyyyyyz[k];

                g_xz_0_xx_yyyyzz[k] = -g_z_0_x_yyyyzz[k] - g_xz_0_x_yyyyzz[k] * ab_x + g_xz_0_x_xyyyyzz[k];

                g_xz_0_xx_yyyzzz[k] = -g_z_0_x_yyyzzz[k] - g_xz_0_x_yyyzzz[k] * ab_x + g_xz_0_x_xyyyzzz[k];

                g_xz_0_xx_yyzzzz[k] = -g_z_0_x_yyzzzz[k] - g_xz_0_x_yyzzzz[k] * ab_x + g_xz_0_x_xyyzzzz[k];

                g_xz_0_xx_yzzzzz[k] = -g_z_0_x_yzzzzz[k] - g_xz_0_x_yzzzzz[k] * ab_x + g_xz_0_x_xyzzzzz[k];

                g_xz_0_xx_zzzzzz[k] = -g_z_0_x_zzzzzz[k] - g_xz_0_x_zzzzzz[k] * ab_x + g_xz_0_x_xzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xy_xxxxxx = cbuffer.data(di_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_xy_xxxxxy = cbuffer.data(di_geom_20_off + 365 * ccomps * dcomps);

            auto g_xz_0_xy_xxxxxz = cbuffer.data(di_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_xy_xxxxyy = cbuffer.data(di_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_xy_xxxxyz = cbuffer.data(di_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_xy_xxxxzz = cbuffer.data(di_geom_20_off + 369 * ccomps * dcomps);

            auto g_xz_0_xy_xxxyyy = cbuffer.data(di_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_xy_xxxyyz = cbuffer.data(di_geom_20_off + 371 * ccomps * dcomps);

            auto g_xz_0_xy_xxxyzz = cbuffer.data(di_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_xy_xxxzzz = cbuffer.data(di_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_xy_xxyyyy = cbuffer.data(di_geom_20_off + 374 * ccomps * dcomps);

            auto g_xz_0_xy_xxyyyz = cbuffer.data(di_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_xy_xxyyzz = cbuffer.data(di_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_xy_xxyzzz = cbuffer.data(di_geom_20_off + 377 * ccomps * dcomps);

            auto g_xz_0_xy_xxzzzz = cbuffer.data(di_geom_20_off + 378 * ccomps * dcomps);

            auto g_xz_0_xy_xyyyyy = cbuffer.data(di_geom_20_off + 379 * ccomps * dcomps);

            auto g_xz_0_xy_xyyyyz = cbuffer.data(di_geom_20_off + 380 * ccomps * dcomps);

            auto g_xz_0_xy_xyyyzz = cbuffer.data(di_geom_20_off + 381 * ccomps * dcomps);

            auto g_xz_0_xy_xyyzzz = cbuffer.data(di_geom_20_off + 382 * ccomps * dcomps);

            auto g_xz_0_xy_xyzzzz = cbuffer.data(di_geom_20_off + 383 * ccomps * dcomps);

            auto g_xz_0_xy_xzzzzz = cbuffer.data(di_geom_20_off + 384 * ccomps * dcomps);

            auto g_xz_0_xy_yyyyyy = cbuffer.data(di_geom_20_off + 385 * ccomps * dcomps);

            auto g_xz_0_xy_yyyyyz = cbuffer.data(di_geom_20_off + 386 * ccomps * dcomps);

            auto g_xz_0_xy_yyyyzz = cbuffer.data(di_geom_20_off + 387 * ccomps * dcomps);

            auto g_xz_0_xy_yyyzzz = cbuffer.data(di_geom_20_off + 388 * ccomps * dcomps);

            auto g_xz_0_xy_yyzzzz = cbuffer.data(di_geom_20_off + 389 * ccomps * dcomps);

            auto g_xz_0_xy_yzzzzz = cbuffer.data(di_geom_20_off + 390 * ccomps * dcomps);

            auto g_xz_0_xy_zzzzzz = cbuffer.data(di_geom_20_off + 391 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_x_xxxxxx, g_xz_0_x_xxxxxxy, g_xz_0_x_xxxxxy, g_xz_0_x_xxxxxyy, g_xz_0_x_xxxxxyz, g_xz_0_x_xxxxxz, g_xz_0_x_xxxxyy, g_xz_0_x_xxxxyyy, g_xz_0_x_xxxxyyz, g_xz_0_x_xxxxyz, g_xz_0_x_xxxxyzz, g_xz_0_x_xxxxzz, g_xz_0_x_xxxyyy, g_xz_0_x_xxxyyyy, g_xz_0_x_xxxyyyz, g_xz_0_x_xxxyyz, g_xz_0_x_xxxyyzz, g_xz_0_x_xxxyzz, g_xz_0_x_xxxyzzz, g_xz_0_x_xxxzzz, g_xz_0_x_xxyyyy, g_xz_0_x_xxyyyyy, g_xz_0_x_xxyyyyz, g_xz_0_x_xxyyyz, g_xz_0_x_xxyyyzz, g_xz_0_x_xxyyzz, g_xz_0_x_xxyyzzz, g_xz_0_x_xxyzzz, g_xz_0_x_xxyzzzz, g_xz_0_x_xxzzzz, g_xz_0_x_xyyyyy, g_xz_0_x_xyyyyyy, g_xz_0_x_xyyyyyz, g_xz_0_x_xyyyyz, g_xz_0_x_xyyyyzz, g_xz_0_x_xyyyzz, g_xz_0_x_xyyyzzz, g_xz_0_x_xyyzzz, g_xz_0_x_xyyzzzz, g_xz_0_x_xyzzzz, g_xz_0_x_xyzzzzz, g_xz_0_x_xzzzzz, g_xz_0_x_yyyyyy, g_xz_0_x_yyyyyyy, g_xz_0_x_yyyyyyz, g_xz_0_x_yyyyyz, g_xz_0_x_yyyyyzz, g_xz_0_x_yyyyzz, g_xz_0_x_yyyyzzz, g_xz_0_x_yyyzzz, g_xz_0_x_yyyzzzz, g_xz_0_x_yyzzzz, g_xz_0_x_yyzzzzz, g_xz_0_x_yzzzzz, g_xz_0_x_yzzzzzz, g_xz_0_x_zzzzzz, g_xz_0_xy_xxxxxx, g_xz_0_xy_xxxxxy, g_xz_0_xy_xxxxxz, g_xz_0_xy_xxxxyy, g_xz_0_xy_xxxxyz, g_xz_0_xy_xxxxzz, g_xz_0_xy_xxxyyy, g_xz_0_xy_xxxyyz, g_xz_0_xy_xxxyzz, g_xz_0_xy_xxxzzz, g_xz_0_xy_xxyyyy, g_xz_0_xy_xxyyyz, g_xz_0_xy_xxyyzz, g_xz_0_xy_xxyzzz, g_xz_0_xy_xxzzzz, g_xz_0_xy_xyyyyy, g_xz_0_xy_xyyyyz, g_xz_0_xy_xyyyzz, g_xz_0_xy_xyyzzz, g_xz_0_xy_xyzzzz, g_xz_0_xy_xzzzzz, g_xz_0_xy_yyyyyy, g_xz_0_xy_yyyyyz, g_xz_0_xy_yyyyzz, g_xz_0_xy_yyyzzz, g_xz_0_xy_yyzzzz, g_xz_0_xy_yzzzzz, g_xz_0_xy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xy_xxxxxx[k] = -g_xz_0_x_xxxxxx[k] * ab_y + g_xz_0_x_xxxxxxy[k];

                g_xz_0_xy_xxxxxy[k] = -g_xz_0_x_xxxxxy[k] * ab_y + g_xz_0_x_xxxxxyy[k];

                g_xz_0_xy_xxxxxz[k] = -g_xz_0_x_xxxxxz[k] * ab_y + g_xz_0_x_xxxxxyz[k];

                g_xz_0_xy_xxxxyy[k] = -g_xz_0_x_xxxxyy[k] * ab_y + g_xz_0_x_xxxxyyy[k];

                g_xz_0_xy_xxxxyz[k] = -g_xz_0_x_xxxxyz[k] * ab_y + g_xz_0_x_xxxxyyz[k];

                g_xz_0_xy_xxxxzz[k] = -g_xz_0_x_xxxxzz[k] * ab_y + g_xz_0_x_xxxxyzz[k];

                g_xz_0_xy_xxxyyy[k] = -g_xz_0_x_xxxyyy[k] * ab_y + g_xz_0_x_xxxyyyy[k];

                g_xz_0_xy_xxxyyz[k] = -g_xz_0_x_xxxyyz[k] * ab_y + g_xz_0_x_xxxyyyz[k];

                g_xz_0_xy_xxxyzz[k] = -g_xz_0_x_xxxyzz[k] * ab_y + g_xz_0_x_xxxyyzz[k];

                g_xz_0_xy_xxxzzz[k] = -g_xz_0_x_xxxzzz[k] * ab_y + g_xz_0_x_xxxyzzz[k];

                g_xz_0_xy_xxyyyy[k] = -g_xz_0_x_xxyyyy[k] * ab_y + g_xz_0_x_xxyyyyy[k];

                g_xz_0_xy_xxyyyz[k] = -g_xz_0_x_xxyyyz[k] * ab_y + g_xz_0_x_xxyyyyz[k];

                g_xz_0_xy_xxyyzz[k] = -g_xz_0_x_xxyyzz[k] * ab_y + g_xz_0_x_xxyyyzz[k];

                g_xz_0_xy_xxyzzz[k] = -g_xz_0_x_xxyzzz[k] * ab_y + g_xz_0_x_xxyyzzz[k];

                g_xz_0_xy_xxzzzz[k] = -g_xz_0_x_xxzzzz[k] * ab_y + g_xz_0_x_xxyzzzz[k];

                g_xz_0_xy_xyyyyy[k] = -g_xz_0_x_xyyyyy[k] * ab_y + g_xz_0_x_xyyyyyy[k];

                g_xz_0_xy_xyyyyz[k] = -g_xz_0_x_xyyyyz[k] * ab_y + g_xz_0_x_xyyyyyz[k];

                g_xz_0_xy_xyyyzz[k] = -g_xz_0_x_xyyyzz[k] * ab_y + g_xz_0_x_xyyyyzz[k];

                g_xz_0_xy_xyyzzz[k] = -g_xz_0_x_xyyzzz[k] * ab_y + g_xz_0_x_xyyyzzz[k];

                g_xz_0_xy_xyzzzz[k] = -g_xz_0_x_xyzzzz[k] * ab_y + g_xz_0_x_xyyzzzz[k];

                g_xz_0_xy_xzzzzz[k] = -g_xz_0_x_xzzzzz[k] * ab_y + g_xz_0_x_xyzzzzz[k];

                g_xz_0_xy_yyyyyy[k] = -g_xz_0_x_yyyyyy[k] * ab_y + g_xz_0_x_yyyyyyy[k];

                g_xz_0_xy_yyyyyz[k] = -g_xz_0_x_yyyyyz[k] * ab_y + g_xz_0_x_yyyyyyz[k];

                g_xz_0_xy_yyyyzz[k] = -g_xz_0_x_yyyyzz[k] * ab_y + g_xz_0_x_yyyyyzz[k];

                g_xz_0_xy_yyyzzz[k] = -g_xz_0_x_yyyzzz[k] * ab_y + g_xz_0_x_yyyyzzz[k];

                g_xz_0_xy_yyzzzz[k] = -g_xz_0_x_yyzzzz[k] * ab_y + g_xz_0_x_yyyzzzz[k];

                g_xz_0_xy_yzzzzz[k] = -g_xz_0_x_yzzzzz[k] * ab_y + g_xz_0_x_yyzzzzz[k];

                g_xz_0_xy_zzzzzz[k] = -g_xz_0_x_zzzzzz[k] * ab_y + g_xz_0_x_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xz_xxxxxx = cbuffer.data(di_geom_20_off + 392 * ccomps * dcomps);

            auto g_xz_0_xz_xxxxxy = cbuffer.data(di_geom_20_off + 393 * ccomps * dcomps);

            auto g_xz_0_xz_xxxxxz = cbuffer.data(di_geom_20_off + 394 * ccomps * dcomps);

            auto g_xz_0_xz_xxxxyy = cbuffer.data(di_geom_20_off + 395 * ccomps * dcomps);

            auto g_xz_0_xz_xxxxyz = cbuffer.data(di_geom_20_off + 396 * ccomps * dcomps);

            auto g_xz_0_xz_xxxxzz = cbuffer.data(di_geom_20_off + 397 * ccomps * dcomps);

            auto g_xz_0_xz_xxxyyy = cbuffer.data(di_geom_20_off + 398 * ccomps * dcomps);

            auto g_xz_0_xz_xxxyyz = cbuffer.data(di_geom_20_off + 399 * ccomps * dcomps);

            auto g_xz_0_xz_xxxyzz = cbuffer.data(di_geom_20_off + 400 * ccomps * dcomps);

            auto g_xz_0_xz_xxxzzz = cbuffer.data(di_geom_20_off + 401 * ccomps * dcomps);

            auto g_xz_0_xz_xxyyyy = cbuffer.data(di_geom_20_off + 402 * ccomps * dcomps);

            auto g_xz_0_xz_xxyyyz = cbuffer.data(di_geom_20_off + 403 * ccomps * dcomps);

            auto g_xz_0_xz_xxyyzz = cbuffer.data(di_geom_20_off + 404 * ccomps * dcomps);

            auto g_xz_0_xz_xxyzzz = cbuffer.data(di_geom_20_off + 405 * ccomps * dcomps);

            auto g_xz_0_xz_xxzzzz = cbuffer.data(di_geom_20_off + 406 * ccomps * dcomps);

            auto g_xz_0_xz_xyyyyy = cbuffer.data(di_geom_20_off + 407 * ccomps * dcomps);

            auto g_xz_0_xz_xyyyyz = cbuffer.data(di_geom_20_off + 408 * ccomps * dcomps);

            auto g_xz_0_xz_xyyyzz = cbuffer.data(di_geom_20_off + 409 * ccomps * dcomps);

            auto g_xz_0_xz_xyyzzz = cbuffer.data(di_geom_20_off + 410 * ccomps * dcomps);

            auto g_xz_0_xz_xyzzzz = cbuffer.data(di_geom_20_off + 411 * ccomps * dcomps);

            auto g_xz_0_xz_xzzzzz = cbuffer.data(di_geom_20_off + 412 * ccomps * dcomps);

            auto g_xz_0_xz_yyyyyy = cbuffer.data(di_geom_20_off + 413 * ccomps * dcomps);

            auto g_xz_0_xz_yyyyyz = cbuffer.data(di_geom_20_off + 414 * ccomps * dcomps);

            auto g_xz_0_xz_yyyyzz = cbuffer.data(di_geom_20_off + 415 * ccomps * dcomps);

            auto g_xz_0_xz_yyyzzz = cbuffer.data(di_geom_20_off + 416 * ccomps * dcomps);

            auto g_xz_0_xz_yyzzzz = cbuffer.data(di_geom_20_off + 417 * ccomps * dcomps);

            auto g_xz_0_xz_yzzzzz = cbuffer.data(di_geom_20_off + 418 * ccomps * dcomps);

            auto g_xz_0_xz_zzzzzz = cbuffer.data(di_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xz_xxxxxx, g_xz_0_xz_xxxxxy, g_xz_0_xz_xxxxxz, g_xz_0_xz_xxxxyy, g_xz_0_xz_xxxxyz, g_xz_0_xz_xxxxzz, g_xz_0_xz_xxxyyy, g_xz_0_xz_xxxyyz, g_xz_0_xz_xxxyzz, g_xz_0_xz_xxxzzz, g_xz_0_xz_xxyyyy, g_xz_0_xz_xxyyyz, g_xz_0_xz_xxyyzz, g_xz_0_xz_xxyzzz, g_xz_0_xz_xxzzzz, g_xz_0_xz_xyyyyy, g_xz_0_xz_xyyyyz, g_xz_0_xz_xyyyzz, g_xz_0_xz_xyyzzz, g_xz_0_xz_xyzzzz, g_xz_0_xz_xzzzzz, g_xz_0_xz_yyyyyy, g_xz_0_xz_yyyyyz, g_xz_0_xz_yyyyzz, g_xz_0_xz_yyyzzz, g_xz_0_xz_yyzzzz, g_xz_0_xz_yzzzzz, g_xz_0_xz_zzzzzz, g_xz_0_z_xxxxxx, g_xz_0_z_xxxxxxx, g_xz_0_z_xxxxxxy, g_xz_0_z_xxxxxxz, g_xz_0_z_xxxxxy, g_xz_0_z_xxxxxyy, g_xz_0_z_xxxxxyz, g_xz_0_z_xxxxxz, g_xz_0_z_xxxxxzz, g_xz_0_z_xxxxyy, g_xz_0_z_xxxxyyy, g_xz_0_z_xxxxyyz, g_xz_0_z_xxxxyz, g_xz_0_z_xxxxyzz, g_xz_0_z_xxxxzz, g_xz_0_z_xxxxzzz, g_xz_0_z_xxxyyy, g_xz_0_z_xxxyyyy, g_xz_0_z_xxxyyyz, g_xz_0_z_xxxyyz, g_xz_0_z_xxxyyzz, g_xz_0_z_xxxyzz, g_xz_0_z_xxxyzzz, g_xz_0_z_xxxzzz, g_xz_0_z_xxxzzzz, g_xz_0_z_xxyyyy, g_xz_0_z_xxyyyyy, g_xz_0_z_xxyyyyz, g_xz_0_z_xxyyyz, g_xz_0_z_xxyyyzz, g_xz_0_z_xxyyzz, g_xz_0_z_xxyyzzz, g_xz_0_z_xxyzzz, g_xz_0_z_xxyzzzz, g_xz_0_z_xxzzzz, g_xz_0_z_xxzzzzz, g_xz_0_z_xyyyyy, g_xz_0_z_xyyyyyy, g_xz_0_z_xyyyyyz, g_xz_0_z_xyyyyz, g_xz_0_z_xyyyyzz, g_xz_0_z_xyyyzz, g_xz_0_z_xyyyzzz, g_xz_0_z_xyyzzz, g_xz_0_z_xyyzzzz, g_xz_0_z_xyzzzz, g_xz_0_z_xyzzzzz, g_xz_0_z_xzzzzz, g_xz_0_z_xzzzzzz, g_xz_0_z_yyyyyy, g_xz_0_z_yyyyyz, g_xz_0_z_yyyyzz, g_xz_0_z_yyyzzz, g_xz_0_z_yyzzzz, g_xz_0_z_yzzzzz, g_xz_0_z_zzzzzz, g_z_0_z_xxxxxx, g_z_0_z_xxxxxy, g_z_0_z_xxxxxz, g_z_0_z_xxxxyy, g_z_0_z_xxxxyz, g_z_0_z_xxxxzz, g_z_0_z_xxxyyy, g_z_0_z_xxxyyz, g_z_0_z_xxxyzz, g_z_0_z_xxxzzz, g_z_0_z_xxyyyy, g_z_0_z_xxyyyz, g_z_0_z_xxyyzz, g_z_0_z_xxyzzz, g_z_0_z_xxzzzz, g_z_0_z_xyyyyy, g_z_0_z_xyyyyz, g_z_0_z_xyyyzz, g_z_0_z_xyyzzz, g_z_0_z_xyzzzz, g_z_0_z_xzzzzz, g_z_0_z_yyyyyy, g_z_0_z_yyyyyz, g_z_0_z_yyyyzz, g_z_0_z_yyyzzz, g_z_0_z_yyzzzz, g_z_0_z_yzzzzz, g_z_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xz_xxxxxx[k] = -g_z_0_z_xxxxxx[k] - g_xz_0_z_xxxxxx[k] * ab_x + g_xz_0_z_xxxxxxx[k];

                g_xz_0_xz_xxxxxy[k] = -g_z_0_z_xxxxxy[k] - g_xz_0_z_xxxxxy[k] * ab_x + g_xz_0_z_xxxxxxy[k];

                g_xz_0_xz_xxxxxz[k] = -g_z_0_z_xxxxxz[k] - g_xz_0_z_xxxxxz[k] * ab_x + g_xz_0_z_xxxxxxz[k];

                g_xz_0_xz_xxxxyy[k] = -g_z_0_z_xxxxyy[k] - g_xz_0_z_xxxxyy[k] * ab_x + g_xz_0_z_xxxxxyy[k];

                g_xz_0_xz_xxxxyz[k] = -g_z_0_z_xxxxyz[k] - g_xz_0_z_xxxxyz[k] * ab_x + g_xz_0_z_xxxxxyz[k];

                g_xz_0_xz_xxxxzz[k] = -g_z_0_z_xxxxzz[k] - g_xz_0_z_xxxxzz[k] * ab_x + g_xz_0_z_xxxxxzz[k];

                g_xz_0_xz_xxxyyy[k] = -g_z_0_z_xxxyyy[k] - g_xz_0_z_xxxyyy[k] * ab_x + g_xz_0_z_xxxxyyy[k];

                g_xz_0_xz_xxxyyz[k] = -g_z_0_z_xxxyyz[k] - g_xz_0_z_xxxyyz[k] * ab_x + g_xz_0_z_xxxxyyz[k];

                g_xz_0_xz_xxxyzz[k] = -g_z_0_z_xxxyzz[k] - g_xz_0_z_xxxyzz[k] * ab_x + g_xz_0_z_xxxxyzz[k];

                g_xz_0_xz_xxxzzz[k] = -g_z_0_z_xxxzzz[k] - g_xz_0_z_xxxzzz[k] * ab_x + g_xz_0_z_xxxxzzz[k];

                g_xz_0_xz_xxyyyy[k] = -g_z_0_z_xxyyyy[k] - g_xz_0_z_xxyyyy[k] * ab_x + g_xz_0_z_xxxyyyy[k];

                g_xz_0_xz_xxyyyz[k] = -g_z_0_z_xxyyyz[k] - g_xz_0_z_xxyyyz[k] * ab_x + g_xz_0_z_xxxyyyz[k];

                g_xz_0_xz_xxyyzz[k] = -g_z_0_z_xxyyzz[k] - g_xz_0_z_xxyyzz[k] * ab_x + g_xz_0_z_xxxyyzz[k];

                g_xz_0_xz_xxyzzz[k] = -g_z_0_z_xxyzzz[k] - g_xz_0_z_xxyzzz[k] * ab_x + g_xz_0_z_xxxyzzz[k];

                g_xz_0_xz_xxzzzz[k] = -g_z_0_z_xxzzzz[k] - g_xz_0_z_xxzzzz[k] * ab_x + g_xz_0_z_xxxzzzz[k];

                g_xz_0_xz_xyyyyy[k] = -g_z_0_z_xyyyyy[k] - g_xz_0_z_xyyyyy[k] * ab_x + g_xz_0_z_xxyyyyy[k];

                g_xz_0_xz_xyyyyz[k] = -g_z_0_z_xyyyyz[k] - g_xz_0_z_xyyyyz[k] * ab_x + g_xz_0_z_xxyyyyz[k];

                g_xz_0_xz_xyyyzz[k] = -g_z_0_z_xyyyzz[k] - g_xz_0_z_xyyyzz[k] * ab_x + g_xz_0_z_xxyyyzz[k];

                g_xz_0_xz_xyyzzz[k] = -g_z_0_z_xyyzzz[k] - g_xz_0_z_xyyzzz[k] * ab_x + g_xz_0_z_xxyyzzz[k];

                g_xz_0_xz_xyzzzz[k] = -g_z_0_z_xyzzzz[k] - g_xz_0_z_xyzzzz[k] * ab_x + g_xz_0_z_xxyzzzz[k];

                g_xz_0_xz_xzzzzz[k] = -g_z_0_z_xzzzzz[k] - g_xz_0_z_xzzzzz[k] * ab_x + g_xz_0_z_xxzzzzz[k];

                g_xz_0_xz_yyyyyy[k] = -g_z_0_z_yyyyyy[k] - g_xz_0_z_yyyyyy[k] * ab_x + g_xz_0_z_xyyyyyy[k];

                g_xz_0_xz_yyyyyz[k] = -g_z_0_z_yyyyyz[k] - g_xz_0_z_yyyyyz[k] * ab_x + g_xz_0_z_xyyyyyz[k];

                g_xz_0_xz_yyyyzz[k] = -g_z_0_z_yyyyzz[k] - g_xz_0_z_yyyyzz[k] * ab_x + g_xz_0_z_xyyyyzz[k];

                g_xz_0_xz_yyyzzz[k] = -g_z_0_z_yyyzzz[k] - g_xz_0_z_yyyzzz[k] * ab_x + g_xz_0_z_xyyyzzz[k];

                g_xz_0_xz_yyzzzz[k] = -g_z_0_z_yyzzzz[k] - g_xz_0_z_yyzzzz[k] * ab_x + g_xz_0_z_xyyzzzz[k];

                g_xz_0_xz_yzzzzz[k] = -g_z_0_z_yzzzzz[k] - g_xz_0_z_yzzzzz[k] * ab_x + g_xz_0_z_xyzzzzz[k];

                g_xz_0_xz_zzzzzz[k] = -g_z_0_z_zzzzzz[k] - g_xz_0_z_zzzzzz[k] * ab_x + g_xz_0_z_xzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yy_xxxxxx = cbuffer.data(di_geom_20_off + 420 * ccomps * dcomps);

            auto g_xz_0_yy_xxxxxy = cbuffer.data(di_geom_20_off + 421 * ccomps * dcomps);

            auto g_xz_0_yy_xxxxxz = cbuffer.data(di_geom_20_off + 422 * ccomps * dcomps);

            auto g_xz_0_yy_xxxxyy = cbuffer.data(di_geom_20_off + 423 * ccomps * dcomps);

            auto g_xz_0_yy_xxxxyz = cbuffer.data(di_geom_20_off + 424 * ccomps * dcomps);

            auto g_xz_0_yy_xxxxzz = cbuffer.data(di_geom_20_off + 425 * ccomps * dcomps);

            auto g_xz_0_yy_xxxyyy = cbuffer.data(di_geom_20_off + 426 * ccomps * dcomps);

            auto g_xz_0_yy_xxxyyz = cbuffer.data(di_geom_20_off + 427 * ccomps * dcomps);

            auto g_xz_0_yy_xxxyzz = cbuffer.data(di_geom_20_off + 428 * ccomps * dcomps);

            auto g_xz_0_yy_xxxzzz = cbuffer.data(di_geom_20_off + 429 * ccomps * dcomps);

            auto g_xz_0_yy_xxyyyy = cbuffer.data(di_geom_20_off + 430 * ccomps * dcomps);

            auto g_xz_0_yy_xxyyyz = cbuffer.data(di_geom_20_off + 431 * ccomps * dcomps);

            auto g_xz_0_yy_xxyyzz = cbuffer.data(di_geom_20_off + 432 * ccomps * dcomps);

            auto g_xz_0_yy_xxyzzz = cbuffer.data(di_geom_20_off + 433 * ccomps * dcomps);

            auto g_xz_0_yy_xxzzzz = cbuffer.data(di_geom_20_off + 434 * ccomps * dcomps);

            auto g_xz_0_yy_xyyyyy = cbuffer.data(di_geom_20_off + 435 * ccomps * dcomps);

            auto g_xz_0_yy_xyyyyz = cbuffer.data(di_geom_20_off + 436 * ccomps * dcomps);

            auto g_xz_0_yy_xyyyzz = cbuffer.data(di_geom_20_off + 437 * ccomps * dcomps);

            auto g_xz_0_yy_xyyzzz = cbuffer.data(di_geom_20_off + 438 * ccomps * dcomps);

            auto g_xz_0_yy_xyzzzz = cbuffer.data(di_geom_20_off + 439 * ccomps * dcomps);

            auto g_xz_0_yy_xzzzzz = cbuffer.data(di_geom_20_off + 440 * ccomps * dcomps);

            auto g_xz_0_yy_yyyyyy = cbuffer.data(di_geom_20_off + 441 * ccomps * dcomps);

            auto g_xz_0_yy_yyyyyz = cbuffer.data(di_geom_20_off + 442 * ccomps * dcomps);

            auto g_xz_0_yy_yyyyzz = cbuffer.data(di_geom_20_off + 443 * ccomps * dcomps);

            auto g_xz_0_yy_yyyzzz = cbuffer.data(di_geom_20_off + 444 * ccomps * dcomps);

            auto g_xz_0_yy_yyzzzz = cbuffer.data(di_geom_20_off + 445 * ccomps * dcomps);

            auto g_xz_0_yy_yzzzzz = cbuffer.data(di_geom_20_off + 446 * ccomps * dcomps);

            auto g_xz_0_yy_zzzzzz = cbuffer.data(di_geom_20_off + 447 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_y_xxxxxx, g_xz_0_y_xxxxxxy, g_xz_0_y_xxxxxy, g_xz_0_y_xxxxxyy, g_xz_0_y_xxxxxyz, g_xz_0_y_xxxxxz, g_xz_0_y_xxxxyy, g_xz_0_y_xxxxyyy, g_xz_0_y_xxxxyyz, g_xz_0_y_xxxxyz, g_xz_0_y_xxxxyzz, g_xz_0_y_xxxxzz, g_xz_0_y_xxxyyy, g_xz_0_y_xxxyyyy, g_xz_0_y_xxxyyyz, g_xz_0_y_xxxyyz, g_xz_0_y_xxxyyzz, g_xz_0_y_xxxyzz, g_xz_0_y_xxxyzzz, g_xz_0_y_xxxzzz, g_xz_0_y_xxyyyy, g_xz_0_y_xxyyyyy, g_xz_0_y_xxyyyyz, g_xz_0_y_xxyyyz, g_xz_0_y_xxyyyzz, g_xz_0_y_xxyyzz, g_xz_0_y_xxyyzzz, g_xz_0_y_xxyzzz, g_xz_0_y_xxyzzzz, g_xz_0_y_xxzzzz, g_xz_0_y_xyyyyy, g_xz_0_y_xyyyyyy, g_xz_0_y_xyyyyyz, g_xz_0_y_xyyyyz, g_xz_0_y_xyyyyzz, g_xz_0_y_xyyyzz, g_xz_0_y_xyyyzzz, g_xz_0_y_xyyzzz, g_xz_0_y_xyyzzzz, g_xz_0_y_xyzzzz, g_xz_0_y_xyzzzzz, g_xz_0_y_xzzzzz, g_xz_0_y_yyyyyy, g_xz_0_y_yyyyyyy, g_xz_0_y_yyyyyyz, g_xz_0_y_yyyyyz, g_xz_0_y_yyyyyzz, g_xz_0_y_yyyyzz, g_xz_0_y_yyyyzzz, g_xz_0_y_yyyzzz, g_xz_0_y_yyyzzzz, g_xz_0_y_yyzzzz, g_xz_0_y_yyzzzzz, g_xz_0_y_yzzzzz, g_xz_0_y_yzzzzzz, g_xz_0_y_zzzzzz, g_xz_0_yy_xxxxxx, g_xz_0_yy_xxxxxy, g_xz_0_yy_xxxxxz, g_xz_0_yy_xxxxyy, g_xz_0_yy_xxxxyz, g_xz_0_yy_xxxxzz, g_xz_0_yy_xxxyyy, g_xz_0_yy_xxxyyz, g_xz_0_yy_xxxyzz, g_xz_0_yy_xxxzzz, g_xz_0_yy_xxyyyy, g_xz_0_yy_xxyyyz, g_xz_0_yy_xxyyzz, g_xz_0_yy_xxyzzz, g_xz_0_yy_xxzzzz, g_xz_0_yy_xyyyyy, g_xz_0_yy_xyyyyz, g_xz_0_yy_xyyyzz, g_xz_0_yy_xyyzzz, g_xz_0_yy_xyzzzz, g_xz_0_yy_xzzzzz, g_xz_0_yy_yyyyyy, g_xz_0_yy_yyyyyz, g_xz_0_yy_yyyyzz, g_xz_0_yy_yyyzzz, g_xz_0_yy_yyzzzz, g_xz_0_yy_yzzzzz, g_xz_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yy_xxxxxx[k] = -g_xz_0_y_xxxxxx[k] * ab_y + g_xz_0_y_xxxxxxy[k];

                g_xz_0_yy_xxxxxy[k] = -g_xz_0_y_xxxxxy[k] * ab_y + g_xz_0_y_xxxxxyy[k];

                g_xz_0_yy_xxxxxz[k] = -g_xz_0_y_xxxxxz[k] * ab_y + g_xz_0_y_xxxxxyz[k];

                g_xz_0_yy_xxxxyy[k] = -g_xz_0_y_xxxxyy[k] * ab_y + g_xz_0_y_xxxxyyy[k];

                g_xz_0_yy_xxxxyz[k] = -g_xz_0_y_xxxxyz[k] * ab_y + g_xz_0_y_xxxxyyz[k];

                g_xz_0_yy_xxxxzz[k] = -g_xz_0_y_xxxxzz[k] * ab_y + g_xz_0_y_xxxxyzz[k];

                g_xz_0_yy_xxxyyy[k] = -g_xz_0_y_xxxyyy[k] * ab_y + g_xz_0_y_xxxyyyy[k];

                g_xz_0_yy_xxxyyz[k] = -g_xz_0_y_xxxyyz[k] * ab_y + g_xz_0_y_xxxyyyz[k];

                g_xz_0_yy_xxxyzz[k] = -g_xz_0_y_xxxyzz[k] * ab_y + g_xz_0_y_xxxyyzz[k];

                g_xz_0_yy_xxxzzz[k] = -g_xz_0_y_xxxzzz[k] * ab_y + g_xz_0_y_xxxyzzz[k];

                g_xz_0_yy_xxyyyy[k] = -g_xz_0_y_xxyyyy[k] * ab_y + g_xz_0_y_xxyyyyy[k];

                g_xz_0_yy_xxyyyz[k] = -g_xz_0_y_xxyyyz[k] * ab_y + g_xz_0_y_xxyyyyz[k];

                g_xz_0_yy_xxyyzz[k] = -g_xz_0_y_xxyyzz[k] * ab_y + g_xz_0_y_xxyyyzz[k];

                g_xz_0_yy_xxyzzz[k] = -g_xz_0_y_xxyzzz[k] * ab_y + g_xz_0_y_xxyyzzz[k];

                g_xz_0_yy_xxzzzz[k] = -g_xz_0_y_xxzzzz[k] * ab_y + g_xz_0_y_xxyzzzz[k];

                g_xz_0_yy_xyyyyy[k] = -g_xz_0_y_xyyyyy[k] * ab_y + g_xz_0_y_xyyyyyy[k];

                g_xz_0_yy_xyyyyz[k] = -g_xz_0_y_xyyyyz[k] * ab_y + g_xz_0_y_xyyyyyz[k];

                g_xz_0_yy_xyyyzz[k] = -g_xz_0_y_xyyyzz[k] * ab_y + g_xz_0_y_xyyyyzz[k];

                g_xz_0_yy_xyyzzz[k] = -g_xz_0_y_xyyzzz[k] * ab_y + g_xz_0_y_xyyyzzz[k];

                g_xz_0_yy_xyzzzz[k] = -g_xz_0_y_xyzzzz[k] * ab_y + g_xz_0_y_xyyzzzz[k];

                g_xz_0_yy_xzzzzz[k] = -g_xz_0_y_xzzzzz[k] * ab_y + g_xz_0_y_xyzzzzz[k];

                g_xz_0_yy_yyyyyy[k] = -g_xz_0_y_yyyyyy[k] * ab_y + g_xz_0_y_yyyyyyy[k];

                g_xz_0_yy_yyyyyz[k] = -g_xz_0_y_yyyyyz[k] * ab_y + g_xz_0_y_yyyyyyz[k];

                g_xz_0_yy_yyyyzz[k] = -g_xz_0_y_yyyyzz[k] * ab_y + g_xz_0_y_yyyyyzz[k];

                g_xz_0_yy_yyyzzz[k] = -g_xz_0_y_yyyzzz[k] * ab_y + g_xz_0_y_yyyyzzz[k];

                g_xz_0_yy_yyzzzz[k] = -g_xz_0_y_yyzzzz[k] * ab_y + g_xz_0_y_yyyzzzz[k];

                g_xz_0_yy_yzzzzz[k] = -g_xz_0_y_yzzzzz[k] * ab_y + g_xz_0_y_yyzzzzz[k];

                g_xz_0_yy_zzzzzz[k] = -g_xz_0_y_zzzzzz[k] * ab_y + g_xz_0_y_yzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yz_xxxxxx = cbuffer.data(di_geom_20_off + 448 * ccomps * dcomps);

            auto g_xz_0_yz_xxxxxy = cbuffer.data(di_geom_20_off + 449 * ccomps * dcomps);

            auto g_xz_0_yz_xxxxxz = cbuffer.data(di_geom_20_off + 450 * ccomps * dcomps);

            auto g_xz_0_yz_xxxxyy = cbuffer.data(di_geom_20_off + 451 * ccomps * dcomps);

            auto g_xz_0_yz_xxxxyz = cbuffer.data(di_geom_20_off + 452 * ccomps * dcomps);

            auto g_xz_0_yz_xxxxzz = cbuffer.data(di_geom_20_off + 453 * ccomps * dcomps);

            auto g_xz_0_yz_xxxyyy = cbuffer.data(di_geom_20_off + 454 * ccomps * dcomps);

            auto g_xz_0_yz_xxxyyz = cbuffer.data(di_geom_20_off + 455 * ccomps * dcomps);

            auto g_xz_0_yz_xxxyzz = cbuffer.data(di_geom_20_off + 456 * ccomps * dcomps);

            auto g_xz_0_yz_xxxzzz = cbuffer.data(di_geom_20_off + 457 * ccomps * dcomps);

            auto g_xz_0_yz_xxyyyy = cbuffer.data(di_geom_20_off + 458 * ccomps * dcomps);

            auto g_xz_0_yz_xxyyyz = cbuffer.data(di_geom_20_off + 459 * ccomps * dcomps);

            auto g_xz_0_yz_xxyyzz = cbuffer.data(di_geom_20_off + 460 * ccomps * dcomps);

            auto g_xz_0_yz_xxyzzz = cbuffer.data(di_geom_20_off + 461 * ccomps * dcomps);

            auto g_xz_0_yz_xxzzzz = cbuffer.data(di_geom_20_off + 462 * ccomps * dcomps);

            auto g_xz_0_yz_xyyyyy = cbuffer.data(di_geom_20_off + 463 * ccomps * dcomps);

            auto g_xz_0_yz_xyyyyz = cbuffer.data(di_geom_20_off + 464 * ccomps * dcomps);

            auto g_xz_0_yz_xyyyzz = cbuffer.data(di_geom_20_off + 465 * ccomps * dcomps);

            auto g_xz_0_yz_xyyzzz = cbuffer.data(di_geom_20_off + 466 * ccomps * dcomps);

            auto g_xz_0_yz_xyzzzz = cbuffer.data(di_geom_20_off + 467 * ccomps * dcomps);

            auto g_xz_0_yz_xzzzzz = cbuffer.data(di_geom_20_off + 468 * ccomps * dcomps);

            auto g_xz_0_yz_yyyyyy = cbuffer.data(di_geom_20_off + 469 * ccomps * dcomps);

            auto g_xz_0_yz_yyyyyz = cbuffer.data(di_geom_20_off + 470 * ccomps * dcomps);

            auto g_xz_0_yz_yyyyzz = cbuffer.data(di_geom_20_off + 471 * ccomps * dcomps);

            auto g_xz_0_yz_yyyzzz = cbuffer.data(di_geom_20_off + 472 * ccomps * dcomps);

            auto g_xz_0_yz_yyzzzz = cbuffer.data(di_geom_20_off + 473 * ccomps * dcomps);

            auto g_xz_0_yz_yzzzzz = cbuffer.data(di_geom_20_off + 474 * ccomps * dcomps);

            auto g_xz_0_yz_zzzzzz = cbuffer.data(di_geom_20_off + 475 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yz_xxxxxx, g_xz_0_yz_xxxxxy, g_xz_0_yz_xxxxxz, g_xz_0_yz_xxxxyy, g_xz_0_yz_xxxxyz, g_xz_0_yz_xxxxzz, g_xz_0_yz_xxxyyy, g_xz_0_yz_xxxyyz, g_xz_0_yz_xxxyzz, g_xz_0_yz_xxxzzz, g_xz_0_yz_xxyyyy, g_xz_0_yz_xxyyyz, g_xz_0_yz_xxyyzz, g_xz_0_yz_xxyzzz, g_xz_0_yz_xxzzzz, g_xz_0_yz_xyyyyy, g_xz_0_yz_xyyyyz, g_xz_0_yz_xyyyzz, g_xz_0_yz_xyyzzz, g_xz_0_yz_xyzzzz, g_xz_0_yz_xzzzzz, g_xz_0_yz_yyyyyy, g_xz_0_yz_yyyyyz, g_xz_0_yz_yyyyzz, g_xz_0_yz_yyyzzz, g_xz_0_yz_yyzzzz, g_xz_0_yz_yzzzzz, g_xz_0_yz_zzzzzz, g_xz_0_z_xxxxxx, g_xz_0_z_xxxxxxy, g_xz_0_z_xxxxxy, g_xz_0_z_xxxxxyy, g_xz_0_z_xxxxxyz, g_xz_0_z_xxxxxz, g_xz_0_z_xxxxyy, g_xz_0_z_xxxxyyy, g_xz_0_z_xxxxyyz, g_xz_0_z_xxxxyz, g_xz_0_z_xxxxyzz, g_xz_0_z_xxxxzz, g_xz_0_z_xxxyyy, g_xz_0_z_xxxyyyy, g_xz_0_z_xxxyyyz, g_xz_0_z_xxxyyz, g_xz_0_z_xxxyyzz, g_xz_0_z_xxxyzz, g_xz_0_z_xxxyzzz, g_xz_0_z_xxxzzz, g_xz_0_z_xxyyyy, g_xz_0_z_xxyyyyy, g_xz_0_z_xxyyyyz, g_xz_0_z_xxyyyz, g_xz_0_z_xxyyyzz, g_xz_0_z_xxyyzz, g_xz_0_z_xxyyzzz, g_xz_0_z_xxyzzz, g_xz_0_z_xxyzzzz, g_xz_0_z_xxzzzz, g_xz_0_z_xyyyyy, g_xz_0_z_xyyyyyy, g_xz_0_z_xyyyyyz, g_xz_0_z_xyyyyz, g_xz_0_z_xyyyyzz, g_xz_0_z_xyyyzz, g_xz_0_z_xyyyzzz, g_xz_0_z_xyyzzz, g_xz_0_z_xyyzzzz, g_xz_0_z_xyzzzz, g_xz_0_z_xyzzzzz, g_xz_0_z_xzzzzz, g_xz_0_z_yyyyyy, g_xz_0_z_yyyyyyy, g_xz_0_z_yyyyyyz, g_xz_0_z_yyyyyz, g_xz_0_z_yyyyyzz, g_xz_0_z_yyyyzz, g_xz_0_z_yyyyzzz, g_xz_0_z_yyyzzz, g_xz_0_z_yyyzzzz, g_xz_0_z_yyzzzz, g_xz_0_z_yyzzzzz, g_xz_0_z_yzzzzz, g_xz_0_z_yzzzzzz, g_xz_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yz_xxxxxx[k] = -g_xz_0_z_xxxxxx[k] * ab_y + g_xz_0_z_xxxxxxy[k];

                g_xz_0_yz_xxxxxy[k] = -g_xz_0_z_xxxxxy[k] * ab_y + g_xz_0_z_xxxxxyy[k];

                g_xz_0_yz_xxxxxz[k] = -g_xz_0_z_xxxxxz[k] * ab_y + g_xz_0_z_xxxxxyz[k];

                g_xz_0_yz_xxxxyy[k] = -g_xz_0_z_xxxxyy[k] * ab_y + g_xz_0_z_xxxxyyy[k];

                g_xz_0_yz_xxxxyz[k] = -g_xz_0_z_xxxxyz[k] * ab_y + g_xz_0_z_xxxxyyz[k];

                g_xz_0_yz_xxxxzz[k] = -g_xz_0_z_xxxxzz[k] * ab_y + g_xz_0_z_xxxxyzz[k];

                g_xz_0_yz_xxxyyy[k] = -g_xz_0_z_xxxyyy[k] * ab_y + g_xz_0_z_xxxyyyy[k];

                g_xz_0_yz_xxxyyz[k] = -g_xz_0_z_xxxyyz[k] * ab_y + g_xz_0_z_xxxyyyz[k];

                g_xz_0_yz_xxxyzz[k] = -g_xz_0_z_xxxyzz[k] * ab_y + g_xz_0_z_xxxyyzz[k];

                g_xz_0_yz_xxxzzz[k] = -g_xz_0_z_xxxzzz[k] * ab_y + g_xz_0_z_xxxyzzz[k];

                g_xz_0_yz_xxyyyy[k] = -g_xz_0_z_xxyyyy[k] * ab_y + g_xz_0_z_xxyyyyy[k];

                g_xz_0_yz_xxyyyz[k] = -g_xz_0_z_xxyyyz[k] * ab_y + g_xz_0_z_xxyyyyz[k];

                g_xz_0_yz_xxyyzz[k] = -g_xz_0_z_xxyyzz[k] * ab_y + g_xz_0_z_xxyyyzz[k];

                g_xz_0_yz_xxyzzz[k] = -g_xz_0_z_xxyzzz[k] * ab_y + g_xz_0_z_xxyyzzz[k];

                g_xz_0_yz_xxzzzz[k] = -g_xz_0_z_xxzzzz[k] * ab_y + g_xz_0_z_xxyzzzz[k];

                g_xz_0_yz_xyyyyy[k] = -g_xz_0_z_xyyyyy[k] * ab_y + g_xz_0_z_xyyyyyy[k];

                g_xz_0_yz_xyyyyz[k] = -g_xz_0_z_xyyyyz[k] * ab_y + g_xz_0_z_xyyyyyz[k];

                g_xz_0_yz_xyyyzz[k] = -g_xz_0_z_xyyyzz[k] * ab_y + g_xz_0_z_xyyyyzz[k];

                g_xz_0_yz_xyyzzz[k] = -g_xz_0_z_xyyzzz[k] * ab_y + g_xz_0_z_xyyyzzz[k];

                g_xz_0_yz_xyzzzz[k] = -g_xz_0_z_xyzzzz[k] * ab_y + g_xz_0_z_xyyzzzz[k];

                g_xz_0_yz_xzzzzz[k] = -g_xz_0_z_xzzzzz[k] * ab_y + g_xz_0_z_xyzzzzz[k];

                g_xz_0_yz_yyyyyy[k] = -g_xz_0_z_yyyyyy[k] * ab_y + g_xz_0_z_yyyyyyy[k];

                g_xz_0_yz_yyyyyz[k] = -g_xz_0_z_yyyyyz[k] * ab_y + g_xz_0_z_yyyyyyz[k];

                g_xz_0_yz_yyyyzz[k] = -g_xz_0_z_yyyyzz[k] * ab_y + g_xz_0_z_yyyyyzz[k];

                g_xz_0_yz_yyyzzz[k] = -g_xz_0_z_yyyzzz[k] * ab_y + g_xz_0_z_yyyyzzz[k];

                g_xz_0_yz_yyzzzz[k] = -g_xz_0_z_yyzzzz[k] * ab_y + g_xz_0_z_yyyzzzz[k];

                g_xz_0_yz_yzzzzz[k] = -g_xz_0_z_yzzzzz[k] * ab_y + g_xz_0_z_yyzzzzz[k];

                g_xz_0_yz_zzzzzz[k] = -g_xz_0_z_zzzzzz[k] * ab_y + g_xz_0_z_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zz_xxxxxx = cbuffer.data(di_geom_20_off + 476 * ccomps * dcomps);

            auto g_xz_0_zz_xxxxxy = cbuffer.data(di_geom_20_off + 477 * ccomps * dcomps);

            auto g_xz_0_zz_xxxxxz = cbuffer.data(di_geom_20_off + 478 * ccomps * dcomps);

            auto g_xz_0_zz_xxxxyy = cbuffer.data(di_geom_20_off + 479 * ccomps * dcomps);

            auto g_xz_0_zz_xxxxyz = cbuffer.data(di_geom_20_off + 480 * ccomps * dcomps);

            auto g_xz_0_zz_xxxxzz = cbuffer.data(di_geom_20_off + 481 * ccomps * dcomps);

            auto g_xz_0_zz_xxxyyy = cbuffer.data(di_geom_20_off + 482 * ccomps * dcomps);

            auto g_xz_0_zz_xxxyyz = cbuffer.data(di_geom_20_off + 483 * ccomps * dcomps);

            auto g_xz_0_zz_xxxyzz = cbuffer.data(di_geom_20_off + 484 * ccomps * dcomps);

            auto g_xz_0_zz_xxxzzz = cbuffer.data(di_geom_20_off + 485 * ccomps * dcomps);

            auto g_xz_0_zz_xxyyyy = cbuffer.data(di_geom_20_off + 486 * ccomps * dcomps);

            auto g_xz_0_zz_xxyyyz = cbuffer.data(di_geom_20_off + 487 * ccomps * dcomps);

            auto g_xz_0_zz_xxyyzz = cbuffer.data(di_geom_20_off + 488 * ccomps * dcomps);

            auto g_xz_0_zz_xxyzzz = cbuffer.data(di_geom_20_off + 489 * ccomps * dcomps);

            auto g_xz_0_zz_xxzzzz = cbuffer.data(di_geom_20_off + 490 * ccomps * dcomps);

            auto g_xz_0_zz_xyyyyy = cbuffer.data(di_geom_20_off + 491 * ccomps * dcomps);

            auto g_xz_0_zz_xyyyyz = cbuffer.data(di_geom_20_off + 492 * ccomps * dcomps);

            auto g_xz_0_zz_xyyyzz = cbuffer.data(di_geom_20_off + 493 * ccomps * dcomps);

            auto g_xz_0_zz_xyyzzz = cbuffer.data(di_geom_20_off + 494 * ccomps * dcomps);

            auto g_xz_0_zz_xyzzzz = cbuffer.data(di_geom_20_off + 495 * ccomps * dcomps);

            auto g_xz_0_zz_xzzzzz = cbuffer.data(di_geom_20_off + 496 * ccomps * dcomps);

            auto g_xz_0_zz_yyyyyy = cbuffer.data(di_geom_20_off + 497 * ccomps * dcomps);

            auto g_xz_0_zz_yyyyyz = cbuffer.data(di_geom_20_off + 498 * ccomps * dcomps);

            auto g_xz_0_zz_yyyyzz = cbuffer.data(di_geom_20_off + 499 * ccomps * dcomps);

            auto g_xz_0_zz_yyyzzz = cbuffer.data(di_geom_20_off + 500 * ccomps * dcomps);

            auto g_xz_0_zz_yyzzzz = cbuffer.data(di_geom_20_off + 501 * ccomps * dcomps);

            auto g_xz_0_zz_yzzzzz = cbuffer.data(di_geom_20_off + 502 * ccomps * dcomps);

            auto g_xz_0_zz_zzzzzz = cbuffer.data(di_geom_20_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxxxxx, g_x_0_z_xxxxxy, g_x_0_z_xxxxxz, g_x_0_z_xxxxyy, g_x_0_z_xxxxyz, g_x_0_z_xxxxzz, g_x_0_z_xxxyyy, g_x_0_z_xxxyyz, g_x_0_z_xxxyzz, g_x_0_z_xxxzzz, g_x_0_z_xxyyyy, g_x_0_z_xxyyyz, g_x_0_z_xxyyzz, g_x_0_z_xxyzzz, g_x_0_z_xxzzzz, g_x_0_z_xyyyyy, g_x_0_z_xyyyyz, g_x_0_z_xyyyzz, g_x_0_z_xyyzzz, g_x_0_z_xyzzzz, g_x_0_z_xzzzzz, g_x_0_z_yyyyyy, g_x_0_z_yyyyyz, g_x_0_z_yyyyzz, g_x_0_z_yyyzzz, g_x_0_z_yyzzzz, g_x_0_z_yzzzzz, g_x_0_z_zzzzzz, g_xz_0_z_xxxxxx, g_xz_0_z_xxxxxxz, g_xz_0_z_xxxxxy, g_xz_0_z_xxxxxyz, g_xz_0_z_xxxxxz, g_xz_0_z_xxxxxzz, g_xz_0_z_xxxxyy, g_xz_0_z_xxxxyyz, g_xz_0_z_xxxxyz, g_xz_0_z_xxxxyzz, g_xz_0_z_xxxxzz, g_xz_0_z_xxxxzzz, g_xz_0_z_xxxyyy, g_xz_0_z_xxxyyyz, g_xz_0_z_xxxyyz, g_xz_0_z_xxxyyzz, g_xz_0_z_xxxyzz, g_xz_0_z_xxxyzzz, g_xz_0_z_xxxzzz, g_xz_0_z_xxxzzzz, g_xz_0_z_xxyyyy, g_xz_0_z_xxyyyyz, g_xz_0_z_xxyyyz, g_xz_0_z_xxyyyzz, g_xz_0_z_xxyyzz, g_xz_0_z_xxyyzzz, g_xz_0_z_xxyzzz, g_xz_0_z_xxyzzzz, g_xz_0_z_xxzzzz, g_xz_0_z_xxzzzzz, g_xz_0_z_xyyyyy, g_xz_0_z_xyyyyyz, g_xz_0_z_xyyyyz, g_xz_0_z_xyyyyzz, g_xz_0_z_xyyyzz, g_xz_0_z_xyyyzzz, g_xz_0_z_xyyzzz, g_xz_0_z_xyyzzzz, g_xz_0_z_xyzzzz, g_xz_0_z_xyzzzzz, g_xz_0_z_xzzzzz, g_xz_0_z_xzzzzzz, g_xz_0_z_yyyyyy, g_xz_0_z_yyyyyyz, g_xz_0_z_yyyyyz, g_xz_0_z_yyyyyzz, g_xz_0_z_yyyyzz, g_xz_0_z_yyyyzzz, g_xz_0_z_yyyzzz, g_xz_0_z_yyyzzzz, g_xz_0_z_yyzzzz, g_xz_0_z_yyzzzzz, g_xz_0_z_yzzzzz, g_xz_0_z_yzzzzzz, g_xz_0_z_zzzzzz, g_xz_0_z_zzzzzzz, g_xz_0_zz_xxxxxx, g_xz_0_zz_xxxxxy, g_xz_0_zz_xxxxxz, g_xz_0_zz_xxxxyy, g_xz_0_zz_xxxxyz, g_xz_0_zz_xxxxzz, g_xz_0_zz_xxxyyy, g_xz_0_zz_xxxyyz, g_xz_0_zz_xxxyzz, g_xz_0_zz_xxxzzz, g_xz_0_zz_xxyyyy, g_xz_0_zz_xxyyyz, g_xz_0_zz_xxyyzz, g_xz_0_zz_xxyzzz, g_xz_0_zz_xxzzzz, g_xz_0_zz_xyyyyy, g_xz_0_zz_xyyyyz, g_xz_0_zz_xyyyzz, g_xz_0_zz_xyyzzz, g_xz_0_zz_xyzzzz, g_xz_0_zz_xzzzzz, g_xz_0_zz_yyyyyy, g_xz_0_zz_yyyyyz, g_xz_0_zz_yyyyzz, g_xz_0_zz_yyyzzz, g_xz_0_zz_yyzzzz, g_xz_0_zz_yzzzzz, g_xz_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zz_xxxxxx[k] = -g_x_0_z_xxxxxx[k] - g_xz_0_z_xxxxxx[k] * ab_z + g_xz_0_z_xxxxxxz[k];

                g_xz_0_zz_xxxxxy[k] = -g_x_0_z_xxxxxy[k] - g_xz_0_z_xxxxxy[k] * ab_z + g_xz_0_z_xxxxxyz[k];

                g_xz_0_zz_xxxxxz[k] = -g_x_0_z_xxxxxz[k] - g_xz_0_z_xxxxxz[k] * ab_z + g_xz_0_z_xxxxxzz[k];

                g_xz_0_zz_xxxxyy[k] = -g_x_0_z_xxxxyy[k] - g_xz_0_z_xxxxyy[k] * ab_z + g_xz_0_z_xxxxyyz[k];

                g_xz_0_zz_xxxxyz[k] = -g_x_0_z_xxxxyz[k] - g_xz_0_z_xxxxyz[k] * ab_z + g_xz_0_z_xxxxyzz[k];

                g_xz_0_zz_xxxxzz[k] = -g_x_0_z_xxxxzz[k] - g_xz_0_z_xxxxzz[k] * ab_z + g_xz_0_z_xxxxzzz[k];

                g_xz_0_zz_xxxyyy[k] = -g_x_0_z_xxxyyy[k] - g_xz_0_z_xxxyyy[k] * ab_z + g_xz_0_z_xxxyyyz[k];

                g_xz_0_zz_xxxyyz[k] = -g_x_0_z_xxxyyz[k] - g_xz_0_z_xxxyyz[k] * ab_z + g_xz_0_z_xxxyyzz[k];

                g_xz_0_zz_xxxyzz[k] = -g_x_0_z_xxxyzz[k] - g_xz_0_z_xxxyzz[k] * ab_z + g_xz_0_z_xxxyzzz[k];

                g_xz_0_zz_xxxzzz[k] = -g_x_0_z_xxxzzz[k] - g_xz_0_z_xxxzzz[k] * ab_z + g_xz_0_z_xxxzzzz[k];

                g_xz_0_zz_xxyyyy[k] = -g_x_0_z_xxyyyy[k] - g_xz_0_z_xxyyyy[k] * ab_z + g_xz_0_z_xxyyyyz[k];

                g_xz_0_zz_xxyyyz[k] = -g_x_0_z_xxyyyz[k] - g_xz_0_z_xxyyyz[k] * ab_z + g_xz_0_z_xxyyyzz[k];

                g_xz_0_zz_xxyyzz[k] = -g_x_0_z_xxyyzz[k] - g_xz_0_z_xxyyzz[k] * ab_z + g_xz_0_z_xxyyzzz[k];

                g_xz_0_zz_xxyzzz[k] = -g_x_0_z_xxyzzz[k] - g_xz_0_z_xxyzzz[k] * ab_z + g_xz_0_z_xxyzzzz[k];

                g_xz_0_zz_xxzzzz[k] = -g_x_0_z_xxzzzz[k] - g_xz_0_z_xxzzzz[k] * ab_z + g_xz_0_z_xxzzzzz[k];

                g_xz_0_zz_xyyyyy[k] = -g_x_0_z_xyyyyy[k] - g_xz_0_z_xyyyyy[k] * ab_z + g_xz_0_z_xyyyyyz[k];

                g_xz_0_zz_xyyyyz[k] = -g_x_0_z_xyyyyz[k] - g_xz_0_z_xyyyyz[k] * ab_z + g_xz_0_z_xyyyyzz[k];

                g_xz_0_zz_xyyyzz[k] = -g_x_0_z_xyyyzz[k] - g_xz_0_z_xyyyzz[k] * ab_z + g_xz_0_z_xyyyzzz[k];

                g_xz_0_zz_xyyzzz[k] = -g_x_0_z_xyyzzz[k] - g_xz_0_z_xyyzzz[k] * ab_z + g_xz_0_z_xyyzzzz[k];

                g_xz_0_zz_xyzzzz[k] = -g_x_0_z_xyzzzz[k] - g_xz_0_z_xyzzzz[k] * ab_z + g_xz_0_z_xyzzzzz[k];

                g_xz_0_zz_xzzzzz[k] = -g_x_0_z_xzzzzz[k] - g_xz_0_z_xzzzzz[k] * ab_z + g_xz_0_z_xzzzzzz[k];

                g_xz_0_zz_yyyyyy[k] = -g_x_0_z_yyyyyy[k] - g_xz_0_z_yyyyyy[k] * ab_z + g_xz_0_z_yyyyyyz[k];

                g_xz_0_zz_yyyyyz[k] = -g_x_0_z_yyyyyz[k] - g_xz_0_z_yyyyyz[k] * ab_z + g_xz_0_z_yyyyyzz[k];

                g_xz_0_zz_yyyyzz[k] = -g_x_0_z_yyyyzz[k] - g_xz_0_z_yyyyzz[k] * ab_z + g_xz_0_z_yyyyzzz[k];

                g_xz_0_zz_yyyzzz[k] = -g_x_0_z_yyyzzz[k] - g_xz_0_z_yyyzzz[k] * ab_z + g_xz_0_z_yyyzzzz[k];

                g_xz_0_zz_yyzzzz[k] = -g_x_0_z_yyzzzz[k] - g_xz_0_z_yyzzzz[k] * ab_z + g_xz_0_z_yyzzzzz[k];

                g_xz_0_zz_yzzzzz[k] = -g_x_0_z_yzzzzz[k] - g_xz_0_z_yzzzzz[k] * ab_z + g_xz_0_z_yzzzzzz[k];

                g_xz_0_zz_zzzzzz[k] = -g_x_0_z_zzzzzz[k] - g_xz_0_z_zzzzzz[k] * ab_z + g_xz_0_z_zzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xx_xxxxxx = cbuffer.data(di_geom_20_off + 504 * ccomps * dcomps);

            auto g_yy_0_xx_xxxxxy = cbuffer.data(di_geom_20_off + 505 * ccomps * dcomps);

            auto g_yy_0_xx_xxxxxz = cbuffer.data(di_geom_20_off + 506 * ccomps * dcomps);

            auto g_yy_0_xx_xxxxyy = cbuffer.data(di_geom_20_off + 507 * ccomps * dcomps);

            auto g_yy_0_xx_xxxxyz = cbuffer.data(di_geom_20_off + 508 * ccomps * dcomps);

            auto g_yy_0_xx_xxxxzz = cbuffer.data(di_geom_20_off + 509 * ccomps * dcomps);

            auto g_yy_0_xx_xxxyyy = cbuffer.data(di_geom_20_off + 510 * ccomps * dcomps);

            auto g_yy_0_xx_xxxyyz = cbuffer.data(di_geom_20_off + 511 * ccomps * dcomps);

            auto g_yy_0_xx_xxxyzz = cbuffer.data(di_geom_20_off + 512 * ccomps * dcomps);

            auto g_yy_0_xx_xxxzzz = cbuffer.data(di_geom_20_off + 513 * ccomps * dcomps);

            auto g_yy_0_xx_xxyyyy = cbuffer.data(di_geom_20_off + 514 * ccomps * dcomps);

            auto g_yy_0_xx_xxyyyz = cbuffer.data(di_geom_20_off + 515 * ccomps * dcomps);

            auto g_yy_0_xx_xxyyzz = cbuffer.data(di_geom_20_off + 516 * ccomps * dcomps);

            auto g_yy_0_xx_xxyzzz = cbuffer.data(di_geom_20_off + 517 * ccomps * dcomps);

            auto g_yy_0_xx_xxzzzz = cbuffer.data(di_geom_20_off + 518 * ccomps * dcomps);

            auto g_yy_0_xx_xyyyyy = cbuffer.data(di_geom_20_off + 519 * ccomps * dcomps);

            auto g_yy_0_xx_xyyyyz = cbuffer.data(di_geom_20_off + 520 * ccomps * dcomps);

            auto g_yy_0_xx_xyyyzz = cbuffer.data(di_geom_20_off + 521 * ccomps * dcomps);

            auto g_yy_0_xx_xyyzzz = cbuffer.data(di_geom_20_off + 522 * ccomps * dcomps);

            auto g_yy_0_xx_xyzzzz = cbuffer.data(di_geom_20_off + 523 * ccomps * dcomps);

            auto g_yy_0_xx_xzzzzz = cbuffer.data(di_geom_20_off + 524 * ccomps * dcomps);

            auto g_yy_0_xx_yyyyyy = cbuffer.data(di_geom_20_off + 525 * ccomps * dcomps);

            auto g_yy_0_xx_yyyyyz = cbuffer.data(di_geom_20_off + 526 * ccomps * dcomps);

            auto g_yy_0_xx_yyyyzz = cbuffer.data(di_geom_20_off + 527 * ccomps * dcomps);

            auto g_yy_0_xx_yyyzzz = cbuffer.data(di_geom_20_off + 528 * ccomps * dcomps);

            auto g_yy_0_xx_yyzzzz = cbuffer.data(di_geom_20_off + 529 * ccomps * dcomps);

            auto g_yy_0_xx_yzzzzz = cbuffer.data(di_geom_20_off + 530 * ccomps * dcomps);

            auto g_yy_0_xx_zzzzzz = cbuffer.data(di_geom_20_off + 531 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_x_xxxxxx, g_yy_0_x_xxxxxxx, g_yy_0_x_xxxxxxy, g_yy_0_x_xxxxxxz, g_yy_0_x_xxxxxy, g_yy_0_x_xxxxxyy, g_yy_0_x_xxxxxyz, g_yy_0_x_xxxxxz, g_yy_0_x_xxxxxzz, g_yy_0_x_xxxxyy, g_yy_0_x_xxxxyyy, g_yy_0_x_xxxxyyz, g_yy_0_x_xxxxyz, g_yy_0_x_xxxxyzz, g_yy_0_x_xxxxzz, g_yy_0_x_xxxxzzz, g_yy_0_x_xxxyyy, g_yy_0_x_xxxyyyy, g_yy_0_x_xxxyyyz, g_yy_0_x_xxxyyz, g_yy_0_x_xxxyyzz, g_yy_0_x_xxxyzz, g_yy_0_x_xxxyzzz, g_yy_0_x_xxxzzz, g_yy_0_x_xxxzzzz, g_yy_0_x_xxyyyy, g_yy_0_x_xxyyyyy, g_yy_0_x_xxyyyyz, g_yy_0_x_xxyyyz, g_yy_0_x_xxyyyzz, g_yy_0_x_xxyyzz, g_yy_0_x_xxyyzzz, g_yy_0_x_xxyzzz, g_yy_0_x_xxyzzzz, g_yy_0_x_xxzzzz, g_yy_0_x_xxzzzzz, g_yy_0_x_xyyyyy, g_yy_0_x_xyyyyyy, g_yy_0_x_xyyyyyz, g_yy_0_x_xyyyyz, g_yy_0_x_xyyyyzz, g_yy_0_x_xyyyzz, g_yy_0_x_xyyyzzz, g_yy_0_x_xyyzzz, g_yy_0_x_xyyzzzz, g_yy_0_x_xyzzzz, g_yy_0_x_xyzzzzz, g_yy_0_x_xzzzzz, g_yy_0_x_xzzzzzz, g_yy_0_x_yyyyyy, g_yy_0_x_yyyyyz, g_yy_0_x_yyyyzz, g_yy_0_x_yyyzzz, g_yy_0_x_yyzzzz, g_yy_0_x_yzzzzz, g_yy_0_x_zzzzzz, g_yy_0_xx_xxxxxx, g_yy_0_xx_xxxxxy, g_yy_0_xx_xxxxxz, g_yy_0_xx_xxxxyy, g_yy_0_xx_xxxxyz, g_yy_0_xx_xxxxzz, g_yy_0_xx_xxxyyy, g_yy_0_xx_xxxyyz, g_yy_0_xx_xxxyzz, g_yy_0_xx_xxxzzz, g_yy_0_xx_xxyyyy, g_yy_0_xx_xxyyyz, g_yy_0_xx_xxyyzz, g_yy_0_xx_xxyzzz, g_yy_0_xx_xxzzzz, g_yy_0_xx_xyyyyy, g_yy_0_xx_xyyyyz, g_yy_0_xx_xyyyzz, g_yy_0_xx_xyyzzz, g_yy_0_xx_xyzzzz, g_yy_0_xx_xzzzzz, g_yy_0_xx_yyyyyy, g_yy_0_xx_yyyyyz, g_yy_0_xx_yyyyzz, g_yy_0_xx_yyyzzz, g_yy_0_xx_yyzzzz, g_yy_0_xx_yzzzzz, g_yy_0_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xx_xxxxxx[k] = -g_yy_0_x_xxxxxx[k] * ab_x + g_yy_0_x_xxxxxxx[k];

                g_yy_0_xx_xxxxxy[k] = -g_yy_0_x_xxxxxy[k] * ab_x + g_yy_0_x_xxxxxxy[k];

                g_yy_0_xx_xxxxxz[k] = -g_yy_0_x_xxxxxz[k] * ab_x + g_yy_0_x_xxxxxxz[k];

                g_yy_0_xx_xxxxyy[k] = -g_yy_0_x_xxxxyy[k] * ab_x + g_yy_0_x_xxxxxyy[k];

                g_yy_0_xx_xxxxyz[k] = -g_yy_0_x_xxxxyz[k] * ab_x + g_yy_0_x_xxxxxyz[k];

                g_yy_0_xx_xxxxzz[k] = -g_yy_0_x_xxxxzz[k] * ab_x + g_yy_0_x_xxxxxzz[k];

                g_yy_0_xx_xxxyyy[k] = -g_yy_0_x_xxxyyy[k] * ab_x + g_yy_0_x_xxxxyyy[k];

                g_yy_0_xx_xxxyyz[k] = -g_yy_0_x_xxxyyz[k] * ab_x + g_yy_0_x_xxxxyyz[k];

                g_yy_0_xx_xxxyzz[k] = -g_yy_0_x_xxxyzz[k] * ab_x + g_yy_0_x_xxxxyzz[k];

                g_yy_0_xx_xxxzzz[k] = -g_yy_0_x_xxxzzz[k] * ab_x + g_yy_0_x_xxxxzzz[k];

                g_yy_0_xx_xxyyyy[k] = -g_yy_0_x_xxyyyy[k] * ab_x + g_yy_0_x_xxxyyyy[k];

                g_yy_0_xx_xxyyyz[k] = -g_yy_0_x_xxyyyz[k] * ab_x + g_yy_0_x_xxxyyyz[k];

                g_yy_0_xx_xxyyzz[k] = -g_yy_0_x_xxyyzz[k] * ab_x + g_yy_0_x_xxxyyzz[k];

                g_yy_0_xx_xxyzzz[k] = -g_yy_0_x_xxyzzz[k] * ab_x + g_yy_0_x_xxxyzzz[k];

                g_yy_0_xx_xxzzzz[k] = -g_yy_0_x_xxzzzz[k] * ab_x + g_yy_0_x_xxxzzzz[k];

                g_yy_0_xx_xyyyyy[k] = -g_yy_0_x_xyyyyy[k] * ab_x + g_yy_0_x_xxyyyyy[k];

                g_yy_0_xx_xyyyyz[k] = -g_yy_0_x_xyyyyz[k] * ab_x + g_yy_0_x_xxyyyyz[k];

                g_yy_0_xx_xyyyzz[k] = -g_yy_0_x_xyyyzz[k] * ab_x + g_yy_0_x_xxyyyzz[k];

                g_yy_0_xx_xyyzzz[k] = -g_yy_0_x_xyyzzz[k] * ab_x + g_yy_0_x_xxyyzzz[k];

                g_yy_0_xx_xyzzzz[k] = -g_yy_0_x_xyzzzz[k] * ab_x + g_yy_0_x_xxyzzzz[k];

                g_yy_0_xx_xzzzzz[k] = -g_yy_0_x_xzzzzz[k] * ab_x + g_yy_0_x_xxzzzzz[k];

                g_yy_0_xx_yyyyyy[k] = -g_yy_0_x_yyyyyy[k] * ab_x + g_yy_0_x_xyyyyyy[k];

                g_yy_0_xx_yyyyyz[k] = -g_yy_0_x_yyyyyz[k] * ab_x + g_yy_0_x_xyyyyyz[k];

                g_yy_0_xx_yyyyzz[k] = -g_yy_0_x_yyyyzz[k] * ab_x + g_yy_0_x_xyyyyzz[k];

                g_yy_0_xx_yyyzzz[k] = -g_yy_0_x_yyyzzz[k] * ab_x + g_yy_0_x_xyyyzzz[k];

                g_yy_0_xx_yyzzzz[k] = -g_yy_0_x_yyzzzz[k] * ab_x + g_yy_0_x_xyyzzzz[k];

                g_yy_0_xx_yzzzzz[k] = -g_yy_0_x_yzzzzz[k] * ab_x + g_yy_0_x_xyzzzzz[k];

                g_yy_0_xx_zzzzzz[k] = -g_yy_0_x_zzzzzz[k] * ab_x + g_yy_0_x_xzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xy_xxxxxx = cbuffer.data(di_geom_20_off + 532 * ccomps * dcomps);

            auto g_yy_0_xy_xxxxxy = cbuffer.data(di_geom_20_off + 533 * ccomps * dcomps);

            auto g_yy_0_xy_xxxxxz = cbuffer.data(di_geom_20_off + 534 * ccomps * dcomps);

            auto g_yy_0_xy_xxxxyy = cbuffer.data(di_geom_20_off + 535 * ccomps * dcomps);

            auto g_yy_0_xy_xxxxyz = cbuffer.data(di_geom_20_off + 536 * ccomps * dcomps);

            auto g_yy_0_xy_xxxxzz = cbuffer.data(di_geom_20_off + 537 * ccomps * dcomps);

            auto g_yy_0_xy_xxxyyy = cbuffer.data(di_geom_20_off + 538 * ccomps * dcomps);

            auto g_yy_0_xy_xxxyyz = cbuffer.data(di_geom_20_off + 539 * ccomps * dcomps);

            auto g_yy_0_xy_xxxyzz = cbuffer.data(di_geom_20_off + 540 * ccomps * dcomps);

            auto g_yy_0_xy_xxxzzz = cbuffer.data(di_geom_20_off + 541 * ccomps * dcomps);

            auto g_yy_0_xy_xxyyyy = cbuffer.data(di_geom_20_off + 542 * ccomps * dcomps);

            auto g_yy_0_xy_xxyyyz = cbuffer.data(di_geom_20_off + 543 * ccomps * dcomps);

            auto g_yy_0_xy_xxyyzz = cbuffer.data(di_geom_20_off + 544 * ccomps * dcomps);

            auto g_yy_0_xy_xxyzzz = cbuffer.data(di_geom_20_off + 545 * ccomps * dcomps);

            auto g_yy_0_xy_xxzzzz = cbuffer.data(di_geom_20_off + 546 * ccomps * dcomps);

            auto g_yy_0_xy_xyyyyy = cbuffer.data(di_geom_20_off + 547 * ccomps * dcomps);

            auto g_yy_0_xy_xyyyyz = cbuffer.data(di_geom_20_off + 548 * ccomps * dcomps);

            auto g_yy_0_xy_xyyyzz = cbuffer.data(di_geom_20_off + 549 * ccomps * dcomps);

            auto g_yy_0_xy_xyyzzz = cbuffer.data(di_geom_20_off + 550 * ccomps * dcomps);

            auto g_yy_0_xy_xyzzzz = cbuffer.data(di_geom_20_off + 551 * ccomps * dcomps);

            auto g_yy_0_xy_xzzzzz = cbuffer.data(di_geom_20_off + 552 * ccomps * dcomps);

            auto g_yy_0_xy_yyyyyy = cbuffer.data(di_geom_20_off + 553 * ccomps * dcomps);

            auto g_yy_0_xy_yyyyyz = cbuffer.data(di_geom_20_off + 554 * ccomps * dcomps);

            auto g_yy_0_xy_yyyyzz = cbuffer.data(di_geom_20_off + 555 * ccomps * dcomps);

            auto g_yy_0_xy_yyyzzz = cbuffer.data(di_geom_20_off + 556 * ccomps * dcomps);

            auto g_yy_0_xy_yyzzzz = cbuffer.data(di_geom_20_off + 557 * ccomps * dcomps);

            auto g_yy_0_xy_yzzzzz = cbuffer.data(di_geom_20_off + 558 * ccomps * dcomps);

            auto g_yy_0_xy_zzzzzz = cbuffer.data(di_geom_20_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xy_xxxxxx, g_yy_0_xy_xxxxxy, g_yy_0_xy_xxxxxz, g_yy_0_xy_xxxxyy, g_yy_0_xy_xxxxyz, g_yy_0_xy_xxxxzz, g_yy_0_xy_xxxyyy, g_yy_0_xy_xxxyyz, g_yy_0_xy_xxxyzz, g_yy_0_xy_xxxzzz, g_yy_0_xy_xxyyyy, g_yy_0_xy_xxyyyz, g_yy_0_xy_xxyyzz, g_yy_0_xy_xxyzzz, g_yy_0_xy_xxzzzz, g_yy_0_xy_xyyyyy, g_yy_0_xy_xyyyyz, g_yy_0_xy_xyyyzz, g_yy_0_xy_xyyzzz, g_yy_0_xy_xyzzzz, g_yy_0_xy_xzzzzz, g_yy_0_xy_yyyyyy, g_yy_0_xy_yyyyyz, g_yy_0_xy_yyyyzz, g_yy_0_xy_yyyzzz, g_yy_0_xy_yyzzzz, g_yy_0_xy_yzzzzz, g_yy_0_xy_zzzzzz, g_yy_0_y_xxxxxx, g_yy_0_y_xxxxxxx, g_yy_0_y_xxxxxxy, g_yy_0_y_xxxxxxz, g_yy_0_y_xxxxxy, g_yy_0_y_xxxxxyy, g_yy_0_y_xxxxxyz, g_yy_0_y_xxxxxz, g_yy_0_y_xxxxxzz, g_yy_0_y_xxxxyy, g_yy_0_y_xxxxyyy, g_yy_0_y_xxxxyyz, g_yy_0_y_xxxxyz, g_yy_0_y_xxxxyzz, g_yy_0_y_xxxxzz, g_yy_0_y_xxxxzzz, g_yy_0_y_xxxyyy, g_yy_0_y_xxxyyyy, g_yy_0_y_xxxyyyz, g_yy_0_y_xxxyyz, g_yy_0_y_xxxyyzz, g_yy_0_y_xxxyzz, g_yy_0_y_xxxyzzz, g_yy_0_y_xxxzzz, g_yy_0_y_xxxzzzz, g_yy_0_y_xxyyyy, g_yy_0_y_xxyyyyy, g_yy_0_y_xxyyyyz, g_yy_0_y_xxyyyz, g_yy_0_y_xxyyyzz, g_yy_0_y_xxyyzz, g_yy_0_y_xxyyzzz, g_yy_0_y_xxyzzz, g_yy_0_y_xxyzzzz, g_yy_0_y_xxzzzz, g_yy_0_y_xxzzzzz, g_yy_0_y_xyyyyy, g_yy_0_y_xyyyyyy, g_yy_0_y_xyyyyyz, g_yy_0_y_xyyyyz, g_yy_0_y_xyyyyzz, g_yy_0_y_xyyyzz, g_yy_0_y_xyyyzzz, g_yy_0_y_xyyzzz, g_yy_0_y_xyyzzzz, g_yy_0_y_xyzzzz, g_yy_0_y_xyzzzzz, g_yy_0_y_xzzzzz, g_yy_0_y_xzzzzzz, g_yy_0_y_yyyyyy, g_yy_0_y_yyyyyz, g_yy_0_y_yyyyzz, g_yy_0_y_yyyzzz, g_yy_0_y_yyzzzz, g_yy_0_y_yzzzzz, g_yy_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xy_xxxxxx[k] = -g_yy_0_y_xxxxxx[k] * ab_x + g_yy_0_y_xxxxxxx[k];

                g_yy_0_xy_xxxxxy[k] = -g_yy_0_y_xxxxxy[k] * ab_x + g_yy_0_y_xxxxxxy[k];

                g_yy_0_xy_xxxxxz[k] = -g_yy_0_y_xxxxxz[k] * ab_x + g_yy_0_y_xxxxxxz[k];

                g_yy_0_xy_xxxxyy[k] = -g_yy_0_y_xxxxyy[k] * ab_x + g_yy_0_y_xxxxxyy[k];

                g_yy_0_xy_xxxxyz[k] = -g_yy_0_y_xxxxyz[k] * ab_x + g_yy_0_y_xxxxxyz[k];

                g_yy_0_xy_xxxxzz[k] = -g_yy_0_y_xxxxzz[k] * ab_x + g_yy_0_y_xxxxxzz[k];

                g_yy_0_xy_xxxyyy[k] = -g_yy_0_y_xxxyyy[k] * ab_x + g_yy_0_y_xxxxyyy[k];

                g_yy_0_xy_xxxyyz[k] = -g_yy_0_y_xxxyyz[k] * ab_x + g_yy_0_y_xxxxyyz[k];

                g_yy_0_xy_xxxyzz[k] = -g_yy_0_y_xxxyzz[k] * ab_x + g_yy_0_y_xxxxyzz[k];

                g_yy_0_xy_xxxzzz[k] = -g_yy_0_y_xxxzzz[k] * ab_x + g_yy_0_y_xxxxzzz[k];

                g_yy_0_xy_xxyyyy[k] = -g_yy_0_y_xxyyyy[k] * ab_x + g_yy_0_y_xxxyyyy[k];

                g_yy_0_xy_xxyyyz[k] = -g_yy_0_y_xxyyyz[k] * ab_x + g_yy_0_y_xxxyyyz[k];

                g_yy_0_xy_xxyyzz[k] = -g_yy_0_y_xxyyzz[k] * ab_x + g_yy_0_y_xxxyyzz[k];

                g_yy_0_xy_xxyzzz[k] = -g_yy_0_y_xxyzzz[k] * ab_x + g_yy_0_y_xxxyzzz[k];

                g_yy_0_xy_xxzzzz[k] = -g_yy_0_y_xxzzzz[k] * ab_x + g_yy_0_y_xxxzzzz[k];

                g_yy_0_xy_xyyyyy[k] = -g_yy_0_y_xyyyyy[k] * ab_x + g_yy_0_y_xxyyyyy[k];

                g_yy_0_xy_xyyyyz[k] = -g_yy_0_y_xyyyyz[k] * ab_x + g_yy_0_y_xxyyyyz[k];

                g_yy_0_xy_xyyyzz[k] = -g_yy_0_y_xyyyzz[k] * ab_x + g_yy_0_y_xxyyyzz[k];

                g_yy_0_xy_xyyzzz[k] = -g_yy_0_y_xyyzzz[k] * ab_x + g_yy_0_y_xxyyzzz[k];

                g_yy_0_xy_xyzzzz[k] = -g_yy_0_y_xyzzzz[k] * ab_x + g_yy_0_y_xxyzzzz[k];

                g_yy_0_xy_xzzzzz[k] = -g_yy_0_y_xzzzzz[k] * ab_x + g_yy_0_y_xxzzzzz[k];

                g_yy_0_xy_yyyyyy[k] = -g_yy_0_y_yyyyyy[k] * ab_x + g_yy_0_y_xyyyyyy[k];

                g_yy_0_xy_yyyyyz[k] = -g_yy_0_y_yyyyyz[k] * ab_x + g_yy_0_y_xyyyyyz[k];

                g_yy_0_xy_yyyyzz[k] = -g_yy_0_y_yyyyzz[k] * ab_x + g_yy_0_y_xyyyyzz[k];

                g_yy_0_xy_yyyzzz[k] = -g_yy_0_y_yyyzzz[k] * ab_x + g_yy_0_y_xyyyzzz[k];

                g_yy_0_xy_yyzzzz[k] = -g_yy_0_y_yyzzzz[k] * ab_x + g_yy_0_y_xyyzzzz[k];

                g_yy_0_xy_yzzzzz[k] = -g_yy_0_y_yzzzzz[k] * ab_x + g_yy_0_y_xyzzzzz[k];

                g_yy_0_xy_zzzzzz[k] = -g_yy_0_y_zzzzzz[k] * ab_x + g_yy_0_y_xzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xz_xxxxxx = cbuffer.data(di_geom_20_off + 560 * ccomps * dcomps);

            auto g_yy_0_xz_xxxxxy = cbuffer.data(di_geom_20_off + 561 * ccomps * dcomps);

            auto g_yy_0_xz_xxxxxz = cbuffer.data(di_geom_20_off + 562 * ccomps * dcomps);

            auto g_yy_0_xz_xxxxyy = cbuffer.data(di_geom_20_off + 563 * ccomps * dcomps);

            auto g_yy_0_xz_xxxxyz = cbuffer.data(di_geom_20_off + 564 * ccomps * dcomps);

            auto g_yy_0_xz_xxxxzz = cbuffer.data(di_geom_20_off + 565 * ccomps * dcomps);

            auto g_yy_0_xz_xxxyyy = cbuffer.data(di_geom_20_off + 566 * ccomps * dcomps);

            auto g_yy_0_xz_xxxyyz = cbuffer.data(di_geom_20_off + 567 * ccomps * dcomps);

            auto g_yy_0_xz_xxxyzz = cbuffer.data(di_geom_20_off + 568 * ccomps * dcomps);

            auto g_yy_0_xz_xxxzzz = cbuffer.data(di_geom_20_off + 569 * ccomps * dcomps);

            auto g_yy_0_xz_xxyyyy = cbuffer.data(di_geom_20_off + 570 * ccomps * dcomps);

            auto g_yy_0_xz_xxyyyz = cbuffer.data(di_geom_20_off + 571 * ccomps * dcomps);

            auto g_yy_0_xz_xxyyzz = cbuffer.data(di_geom_20_off + 572 * ccomps * dcomps);

            auto g_yy_0_xz_xxyzzz = cbuffer.data(di_geom_20_off + 573 * ccomps * dcomps);

            auto g_yy_0_xz_xxzzzz = cbuffer.data(di_geom_20_off + 574 * ccomps * dcomps);

            auto g_yy_0_xz_xyyyyy = cbuffer.data(di_geom_20_off + 575 * ccomps * dcomps);

            auto g_yy_0_xz_xyyyyz = cbuffer.data(di_geom_20_off + 576 * ccomps * dcomps);

            auto g_yy_0_xz_xyyyzz = cbuffer.data(di_geom_20_off + 577 * ccomps * dcomps);

            auto g_yy_0_xz_xyyzzz = cbuffer.data(di_geom_20_off + 578 * ccomps * dcomps);

            auto g_yy_0_xz_xyzzzz = cbuffer.data(di_geom_20_off + 579 * ccomps * dcomps);

            auto g_yy_0_xz_xzzzzz = cbuffer.data(di_geom_20_off + 580 * ccomps * dcomps);

            auto g_yy_0_xz_yyyyyy = cbuffer.data(di_geom_20_off + 581 * ccomps * dcomps);

            auto g_yy_0_xz_yyyyyz = cbuffer.data(di_geom_20_off + 582 * ccomps * dcomps);

            auto g_yy_0_xz_yyyyzz = cbuffer.data(di_geom_20_off + 583 * ccomps * dcomps);

            auto g_yy_0_xz_yyyzzz = cbuffer.data(di_geom_20_off + 584 * ccomps * dcomps);

            auto g_yy_0_xz_yyzzzz = cbuffer.data(di_geom_20_off + 585 * ccomps * dcomps);

            auto g_yy_0_xz_yzzzzz = cbuffer.data(di_geom_20_off + 586 * ccomps * dcomps);

            auto g_yy_0_xz_zzzzzz = cbuffer.data(di_geom_20_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xz_xxxxxx, g_yy_0_xz_xxxxxy, g_yy_0_xz_xxxxxz, g_yy_0_xz_xxxxyy, g_yy_0_xz_xxxxyz, g_yy_0_xz_xxxxzz, g_yy_0_xz_xxxyyy, g_yy_0_xz_xxxyyz, g_yy_0_xz_xxxyzz, g_yy_0_xz_xxxzzz, g_yy_0_xz_xxyyyy, g_yy_0_xz_xxyyyz, g_yy_0_xz_xxyyzz, g_yy_0_xz_xxyzzz, g_yy_0_xz_xxzzzz, g_yy_0_xz_xyyyyy, g_yy_0_xz_xyyyyz, g_yy_0_xz_xyyyzz, g_yy_0_xz_xyyzzz, g_yy_0_xz_xyzzzz, g_yy_0_xz_xzzzzz, g_yy_0_xz_yyyyyy, g_yy_0_xz_yyyyyz, g_yy_0_xz_yyyyzz, g_yy_0_xz_yyyzzz, g_yy_0_xz_yyzzzz, g_yy_0_xz_yzzzzz, g_yy_0_xz_zzzzzz, g_yy_0_z_xxxxxx, g_yy_0_z_xxxxxxx, g_yy_0_z_xxxxxxy, g_yy_0_z_xxxxxxz, g_yy_0_z_xxxxxy, g_yy_0_z_xxxxxyy, g_yy_0_z_xxxxxyz, g_yy_0_z_xxxxxz, g_yy_0_z_xxxxxzz, g_yy_0_z_xxxxyy, g_yy_0_z_xxxxyyy, g_yy_0_z_xxxxyyz, g_yy_0_z_xxxxyz, g_yy_0_z_xxxxyzz, g_yy_0_z_xxxxzz, g_yy_0_z_xxxxzzz, g_yy_0_z_xxxyyy, g_yy_0_z_xxxyyyy, g_yy_0_z_xxxyyyz, g_yy_0_z_xxxyyz, g_yy_0_z_xxxyyzz, g_yy_0_z_xxxyzz, g_yy_0_z_xxxyzzz, g_yy_0_z_xxxzzz, g_yy_0_z_xxxzzzz, g_yy_0_z_xxyyyy, g_yy_0_z_xxyyyyy, g_yy_0_z_xxyyyyz, g_yy_0_z_xxyyyz, g_yy_0_z_xxyyyzz, g_yy_0_z_xxyyzz, g_yy_0_z_xxyyzzz, g_yy_0_z_xxyzzz, g_yy_0_z_xxyzzzz, g_yy_0_z_xxzzzz, g_yy_0_z_xxzzzzz, g_yy_0_z_xyyyyy, g_yy_0_z_xyyyyyy, g_yy_0_z_xyyyyyz, g_yy_0_z_xyyyyz, g_yy_0_z_xyyyyzz, g_yy_0_z_xyyyzz, g_yy_0_z_xyyyzzz, g_yy_0_z_xyyzzz, g_yy_0_z_xyyzzzz, g_yy_0_z_xyzzzz, g_yy_0_z_xyzzzzz, g_yy_0_z_xzzzzz, g_yy_0_z_xzzzzzz, g_yy_0_z_yyyyyy, g_yy_0_z_yyyyyz, g_yy_0_z_yyyyzz, g_yy_0_z_yyyzzz, g_yy_0_z_yyzzzz, g_yy_0_z_yzzzzz, g_yy_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xz_xxxxxx[k] = -g_yy_0_z_xxxxxx[k] * ab_x + g_yy_0_z_xxxxxxx[k];

                g_yy_0_xz_xxxxxy[k] = -g_yy_0_z_xxxxxy[k] * ab_x + g_yy_0_z_xxxxxxy[k];

                g_yy_0_xz_xxxxxz[k] = -g_yy_0_z_xxxxxz[k] * ab_x + g_yy_0_z_xxxxxxz[k];

                g_yy_0_xz_xxxxyy[k] = -g_yy_0_z_xxxxyy[k] * ab_x + g_yy_0_z_xxxxxyy[k];

                g_yy_0_xz_xxxxyz[k] = -g_yy_0_z_xxxxyz[k] * ab_x + g_yy_0_z_xxxxxyz[k];

                g_yy_0_xz_xxxxzz[k] = -g_yy_0_z_xxxxzz[k] * ab_x + g_yy_0_z_xxxxxzz[k];

                g_yy_0_xz_xxxyyy[k] = -g_yy_0_z_xxxyyy[k] * ab_x + g_yy_0_z_xxxxyyy[k];

                g_yy_0_xz_xxxyyz[k] = -g_yy_0_z_xxxyyz[k] * ab_x + g_yy_0_z_xxxxyyz[k];

                g_yy_0_xz_xxxyzz[k] = -g_yy_0_z_xxxyzz[k] * ab_x + g_yy_0_z_xxxxyzz[k];

                g_yy_0_xz_xxxzzz[k] = -g_yy_0_z_xxxzzz[k] * ab_x + g_yy_0_z_xxxxzzz[k];

                g_yy_0_xz_xxyyyy[k] = -g_yy_0_z_xxyyyy[k] * ab_x + g_yy_0_z_xxxyyyy[k];

                g_yy_0_xz_xxyyyz[k] = -g_yy_0_z_xxyyyz[k] * ab_x + g_yy_0_z_xxxyyyz[k];

                g_yy_0_xz_xxyyzz[k] = -g_yy_0_z_xxyyzz[k] * ab_x + g_yy_0_z_xxxyyzz[k];

                g_yy_0_xz_xxyzzz[k] = -g_yy_0_z_xxyzzz[k] * ab_x + g_yy_0_z_xxxyzzz[k];

                g_yy_0_xz_xxzzzz[k] = -g_yy_0_z_xxzzzz[k] * ab_x + g_yy_0_z_xxxzzzz[k];

                g_yy_0_xz_xyyyyy[k] = -g_yy_0_z_xyyyyy[k] * ab_x + g_yy_0_z_xxyyyyy[k];

                g_yy_0_xz_xyyyyz[k] = -g_yy_0_z_xyyyyz[k] * ab_x + g_yy_0_z_xxyyyyz[k];

                g_yy_0_xz_xyyyzz[k] = -g_yy_0_z_xyyyzz[k] * ab_x + g_yy_0_z_xxyyyzz[k];

                g_yy_0_xz_xyyzzz[k] = -g_yy_0_z_xyyzzz[k] * ab_x + g_yy_0_z_xxyyzzz[k];

                g_yy_0_xz_xyzzzz[k] = -g_yy_0_z_xyzzzz[k] * ab_x + g_yy_0_z_xxyzzzz[k];

                g_yy_0_xz_xzzzzz[k] = -g_yy_0_z_xzzzzz[k] * ab_x + g_yy_0_z_xxzzzzz[k];

                g_yy_0_xz_yyyyyy[k] = -g_yy_0_z_yyyyyy[k] * ab_x + g_yy_0_z_xyyyyyy[k];

                g_yy_0_xz_yyyyyz[k] = -g_yy_0_z_yyyyyz[k] * ab_x + g_yy_0_z_xyyyyyz[k];

                g_yy_0_xz_yyyyzz[k] = -g_yy_0_z_yyyyzz[k] * ab_x + g_yy_0_z_xyyyyzz[k];

                g_yy_0_xz_yyyzzz[k] = -g_yy_0_z_yyyzzz[k] * ab_x + g_yy_0_z_xyyyzzz[k];

                g_yy_0_xz_yyzzzz[k] = -g_yy_0_z_yyzzzz[k] * ab_x + g_yy_0_z_xyyzzzz[k];

                g_yy_0_xz_yzzzzz[k] = -g_yy_0_z_yzzzzz[k] * ab_x + g_yy_0_z_xyzzzzz[k];

                g_yy_0_xz_zzzzzz[k] = -g_yy_0_z_zzzzzz[k] * ab_x + g_yy_0_z_xzzzzzz[k];
            }

            /// Set up 588-616 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yy_xxxxxx = cbuffer.data(di_geom_20_off + 588 * ccomps * dcomps);

            auto g_yy_0_yy_xxxxxy = cbuffer.data(di_geom_20_off + 589 * ccomps * dcomps);

            auto g_yy_0_yy_xxxxxz = cbuffer.data(di_geom_20_off + 590 * ccomps * dcomps);

            auto g_yy_0_yy_xxxxyy = cbuffer.data(di_geom_20_off + 591 * ccomps * dcomps);

            auto g_yy_0_yy_xxxxyz = cbuffer.data(di_geom_20_off + 592 * ccomps * dcomps);

            auto g_yy_0_yy_xxxxzz = cbuffer.data(di_geom_20_off + 593 * ccomps * dcomps);

            auto g_yy_0_yy_xxxyyy = cbuffer.data(di_geom_20_off + 594 * ccomps * dcomps);

            auto g_yy_0_yy_xxxyyz = cbuffer.data(di_geom_20_off + 595 * ccomps * dcomps);

            auto g_yy_0_yy_xxxyzz = cbuffer.data(di_geom_20_off + 596 * ccomps * dcomps);

            auto g_yy_0_yy_xxxzzz = cbuffer.data(di_geom_20_off + 597 * ccomps * dcomps);

            auto g_yy_0_yy_xxyyyy = cbuffer.data(di_geom_20_off + 598 * ccomps * dcomps);

            auto g_yy_0_yy_xxyyyz = cbuffer.data(di_geom_20_off + 599 * ccomps * dcomps);

            auto g_yy_0_yy_xxyyzz = cbuffer.data(di_geom_20_off + 600 * ccomps * dcomps);

            auto g_yy_0_yy_xxyzzz = cbuffer.data(di_geom_20_off + 601 * ccomps * dcomps);

            auto g_yy_0_yy_xxzzzz = cbuffer.data(di_geom_20_off + 602 * ccomps * dcomps);

            auto g_yy_0_yy_xyyyyy = cbuffer.data(di_geom_20_off + 603 * ccomps * dcomps);

            auto g_yy_0_yy_xyyyyz = cbuffer.data(di_geom_20_off + 604 * ccomps * dcomps);

            auto g_yy_0_yy_xyyyzz = cbuffer.data(di_geom_20_off + 605 * ccomps * dcomps);

            auto g_yy_0_yy_xyyzzz = cbuffer.data(di_geom_20_off + 606 * ccomps * dcomps);

            auto g_yy_0_yy_xyzzzz = cbuffer.data(di_geom_20_off + 607 * ccomps * dcomps);

            auto g_yy_0_yy_xzzzzz = cbuffer.data(di_geom_20_off + 608 * ccomps * dcomps);

            auto g_yy_0_yy_yyyyyy = cbuffer.data(di_geom_20_off + 609 * ccomps * dcomps);

            auto g_yy_0_yy_yyyyyz = cbuffer.data(di_geom_20_off + 610 * ccomps * dcomps);

            auto g_yy_0_yy_yyyyzz = cbuffer.data(di_geom_20_off + 611 * ccomps * dcomps);

            auto g_yy_0_yy_yyyzzz = cbuffer.data(di_geom_20_off + 612 * ccomps * dcomps);

            auto g_yy_0_yy_yyzzzz = cbuffer.data(di_geom_20_off + 613 * ccomps * dcomps);

            auto g_yy_0_yy_yzzzzz = cbuffer.data(di_geom_20_off + 614 * ccomps * dcomps);

            auto g_yy_0_yy_zzzzzz = cbuffer.data(di_geom_20_off + 615 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxxxxx, g_y_0_y_xxxxxy, g_y_0_y_xxxxxz, g_y_0_y_xxxxyy, g_y_0_y_xxxxyz, g_y_0_y_xxxxzz, g_y_0_y_xxxyyy, g_y_0_y_xxxyyz, g_y_0_y_xxxyzz, g_y_0_y_xxxzzz, g_y_0_y_xxyyyy, g_y_0_y_xxyyyz, g_y_0_y_xxyyzz, g_y_0_y_xxyzzz, g_y_0_y_xxzzzz, g_y_0_y_xyyyyy, g_y_0_y_xyyyyz, g_y_0_y_xyyyzz, g_y_0_y_xyyzzz, g_y_0_y_xyzzzz, g_y_0_y_xzzzzz, g_y_0_y_yyyyyy, g_y_0_y_yyyyyz, g_y_0_y_yyyyzz, g_y_0_y_yyyzzz, g_y_0_y_yyzzzz, g_y_0_y_yzzzzz, g_y_0_y_zzzzzz, g_yy_0_y_xxxxxx, g_yy_0_y_xxxxxxy, g_yy_0_y_xxxxxy, g_yy_0_y_xxxxxyy, g_yy_0_y_xxxxxyz, g_yy_0_y_xxxxxz, g_yy_0_y_xxxxyy, g_yy_0_y_xxxxyyy, g_yy_0_y_xxxxyyz, g_yy_0_y_xxxxyz, g_yy_0_y_xxxxyzz, g_yy_0_y_xxxxzz, g_yy_0_y_xxxyyy, g_yy_0_y_xxxyyyy, g_yy_0_y_xxxyyyz, g_yy_0_y_xxxyyz, g_yy_0_y_xxxyyzz, g_yy_0_y_xxxyzz, g_yy_0_y_xxxyzzz, g_yy_0_y_xxxzzz, g_yy_0_y_xxyyyy, g_yy_0_y_xxyyyyy, g_yy_0_y_xxyyyyz, g_yy_0_y_xxyyyz, g_yy_0_y_xxyyyzz, g_yy_0_y_xxyyzz, g_yy_0_y_xxyyzzz, g_yy_0_y_xxyzzz, g_yy_0_y_xxyzzzz, g_yy_0_y_xxzzzz, g_yy_0_y_xyyyyy, g_yy_0_y_xyyyyyy, g_yy_0_y_xyyyyyz, g_yy_0_y_xyyyyz, g_yy_0_y_xyyyyzz, g_yy_0_y_xyyyzz, g_yy_0_y_xyyyzzz, g_yy_0_y_xyyzzz, g_yy_0_y_xyyzzzz, g_yy_0_y_xyzzzz, g_yy_0_y_xyzzzzz, g_yy_0_y_xzzzzz, g_yy_0_y_yyyyyy, g_yy_0_y_yyyyyyy, g_yy_0_y_yyyyyyz, g_yy_0_y_yyyyyz, g_yy_0_y_yyyyyzz, g_yy_0_y_yyyyzz, g_yy_0_y_yyyyzzz, g_yy_0_y_yyyzzz, g_yy_0_y_yyyzzzz, g_yy_0_y_yyzzzz, g_yy_0_y_yyzzzzz, g_yy_0_y_yzzzzz, g_yy_0_y_yzzzzzz, g_yy_0_y_zzzzzz, g_yy_0_yy_xxxxxx, g_yy_0_yy_xxxxxy, g_yy_0_yy_xxxxxz, g_yy_0_yy_xxxxyy, g_yy_0_yy_xxxxyz, g_yy_0_yy_xxxxzz, g_yy_0_yy_xxxyyy, g_yy_0_yy_xxxyyz, g_yy_0_yy_xxxyzz, g_yy_0_yy_xxxzzz, g_yy_0_yy_xxyyyy, g_yy_0_yy_xxyyyz, g_yy_0_yy_xxyyzz, g_yy_0_yy_xxyzzz, g_yy_0_yy_xxzzzz, g_yy_0_yy_xyyyyy, g_yy_0_yy_xyyyyz, g_yy_0_yy_xyyyzz, g_yy_0_yy_xyyzzz, g_yy_0_yy_xyzzzz, g_yy_0_yy_xzzzzz, g_yy_0_yy_yyyyyy, g_yy_0_yy_yyyyyz, g_yy_0_yy_yyyyzz, g_yy_0_yy_yyyzzz, g_yy_0_yy_yyzzzz, g_yy_0_yy_yzzzzz, g_yy_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yy_xxxxxx[k] = -2.0 * g_y_0_y_xxxxxx[k] - g_yy_0_y_xxxxxx[k] * ab_y + g_yy_0_y_xxxxxxy[k];

                g_yy_0_yy_xxxxxy[k] = -2.0 * g_y_0_y_xxxxxy[k] - g_yy_0_y_xxxxxy[k] * ab_y + g_yy_0_y_xxxxxyy[k];

                g_yy_0_yy_xxxxxz[k] = -2.0 * g_y_0_y_xxxxxz[k] - g_yy_0_y_xxxxxz[k] * ab_y + g_yy_0_y_xxxxxyz[k];

                g_yy_0_yy_xxxxyy[k] = -2.0 * g_y_0_y_xxxxyy[k] - g_yy_0_y_xxxxyy[k] * ab_y + g_yy_0_y_xxxxyyy[k];

                g_yy_0_yy_xxxxyz[k] = -2.0 * g_y_0_y_xxxxyz[k] - g_yy_0_y_xxxxyz[k] * ab_y + g_yy_0_y_xxxxyyz[k];

                g_yy_0_yy_xxxxzz[k] = -2.0 * g_y_0_y_xxxxzz[k] - g_yy_0_y_xxxxzz[k] * ab_y + g_yy_0_y_xxxxyzz[k];

                g_yy_0_yy_xxxyyy[k] = -2.0 * g_y_0_y_xxxyyy[k] - g_yy_0_y_xxxyyy[k] * ab_y + g_yy_0_y_xxxyyyy[k];

                g_yy_0_yy_xxxyyz[k] = -2.0 * g_y_0_y_xxxyyz[k] - g_yy_0_y_xxxyyz[k] * ab_y + g_yy_0_y_xxxyyyz[k];

                g_yy_0_yy_xxxyzz[k] = -2.0 * g_y_0_y_xxxyzz[k] - g_yy_0_y_xxxyzz[k] * ab_y + g_yy_0_y_xxxyyzz[k];

                g_yy_0_yy_xxxzzz[k] = -2.0 * g_y_0_y_xxxzzz[k] - g_yy_0_y_xxxzzz[k] * ab_y + g_yy_0_y_xxxyzzz[k];

                g_yy_0_yy_xxyyyy[k] = -2.0 * g_y_0_y_xxyyyy[k] - g_yy_0_y_xxyyyy[k] * ab_y + g_yy_0_y_xxyyyyy[k];

                g_yy_0_yy_xxyyyz[k] = -2.0 * g_y_0_y_xxyyyz[k] - g_yy_0_y_xxyyyz[k] * ab_y + g_yy_0_y_xxyyyyz[k];

                g_yy_0_yy_xxyyzz[k] = -2.0 * g_y_0_y_xxyyzz[k] - g_yy_0_y_xxyyzz[k] * ab_y + g_yy_0_y_xxyyyzz[k];

                g_yy_0_yy_xxyzzz[k] = -2.0 * g_y_0_y_xxyzzz[k] - g_yy_0_y_xxyzzz[k] * ab_y + g_yy_0_y_xxyyzzz[k];

                g_yy_0_yy_xxzzzz[k] = -2.0 * g_y_0_y_xxzzzz[k] - g_yy_0_y_xxzzzz[k] * ab_y + g_yy_0_y_xxyzzzz[k];

                g_yy_0_yy_xyyyyy[k] = -2.0 * g_y_0_y_xyyyyy[k] - g_yy_0_y_xyyyyy[k] * ab_y + g_yy_0_y_xyyyyyy[k];

                g_yy_0_yy_xyyyyz[k] = -2.0 * g_y_0_y_xyyyyz[k] - g_yy_0_y_xyyyyz[k] * ab_y + g_yy_0_y_xyyyyyz[k];

                g_yy_0_yy_xyyyzz[k] = -2.0 * g_y_0_y_xyyyzz[k] - g_yy_0_y_xyyyzz[k] * ab_y + g_yy_0_y_xyyyyzz[k];

                g_yy_0_yy_xyyzzz[k] = -2.0 * g_y_0_y_xyyzzz[k] - g_yy_0_y_xyyzzz[k] * ab_y + g_yy_0_y_xyyyzzz[k];

                g_yy_0_yy_xyzzzz[k] = -2.0 * g_y_0_y_xyzzzz[k] - g_yy_0_y_xyzzzz[k] * ab_y + g_yy_0_y_xyyzzzz[k];

                g_yy_0_yy_xzzzzz[k] = -2.0 * g_y_0_y_xzzzzz[k] - g_yy_0_y_xzzzzz[k] * ab_y + g_yy_0_y_xyzzzzz[k];

                g_yy_0_yy_yyyyyy[k] = -2.0 * g_y_0_y_yyyyyy[k] - g_yy_0_y_yyyyyy[k] * ab_y + g_yy_0_y_yyyyyyy[k];

                g_yy_0_yy_yyyyyz[k] = -2.0 * g_y_0_y_yyyyyz[k] - g_yy_0_y_yyyyyz[k] * ab_y + g_yy_0_y_yyyyyyz[k];

                g_yy_0_yy_yyyyzz[k] = -2.0 * g_y_0_y_yyyyzz[k] - g_yy_0_y_yyyyzz[k] * ab_y + g_yy_0_y_yyyyyzz[k];

                g_yy_0_yy_yyyzzz[k] = -2.0 * g_y_0_y_yyyzzz[k] - g_yy_0_y_yyyzzz[k] * ab_y + g_yy_0_y_yyyyzzz[k];

                g_yy_0_yy_yyzzzz[k] = -2.0 * g_y_0_y_yyzzzz[k] - g_yy_0_y_yyzzzz[k] * ab_y + g_yy_0_y_yyyzzzz[k];

                g_yy_0_yy_yzzzzz[k] = -2.0 * g_y_0_y_yzzzzz[k] - g_yy_0_y_yzzzzz[k] * ab_y + g_yy_0_y_yyzzzzz[k];

                g_yy_0_yy_zzzzzz[k] = -2.0 * g_y_0_y_zzzzzz[k] - g_yy_0_y_zzzzzz[k] * ab_y + g_yy_0_y_yzzzzzz[k];
            }

            /// Set up 616-644 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yz_xxxxxx = cbuffer.data(di_geom_20_off + 616 * ccomps * dcomps);

            auto g_yy_0_yz_xxxxxy = cbuffer.data(di_geom_20_off + 617 * ccomps * dcomps);

            auto g_yy_0_yz_xxxxxz = cbuffer.data(di_geom_20_off + 618 * ccomps * dcomps);

            auto g_yy_0_yz_xxxxyy = cbuffer.data(di_geom_20_off + 619 * ccomps * dcomps);

            auto g_yy_0_yz_xxxxyz = cbuffer.data(di_geom_20_off + 620 * ccomps * dcomps);

            auto g_yy_0_yz_xxxxzz = cbuffer.data(di_geom_20_off + 621 * ccomps * dcomps);

            auto g_yy_0_yz_xxxyyy = cbuffer.data(di_geom_20_off + 622 * ccomps * dcomps);

            auto g_yy_0_yz_xxxyyz = cbuffer.data(di_geom_20_off + 623 * ccomps * dcomps);

            auto g_yy_0_yz_xxxyzz = cbuffer.data(di_geom_20_off + 624 * ccomps * dcomps);

            auto g_yy_0_yz_xxxzzz = cbuffer.data(di_geom_20_off + 625 * ccomps * dcomps);

            auto g_yy_0_yz_xxyyyy = cbuffer.data(di_geom_20_off + 626 * ccomps * dcomps);

            auto g_yy_0_yz_xxyyyz = cbuffer.data(di_geom_20_off + 627 * ccomps * dcomps);

            auto g_yy_0_yz_xxyyzz = cbuffer.data(di_geom_20_off + 628 * ccomps * dcomps);

            auto g_yy_0_yz_xxyzzz = cbuffer.data(di_geom_20_off + 629 * ccomps * dcomps);

            auto g_yy_0_yz_xxzzzz = cbuffer.data(di_geom_20_off + 630 * ccomps * dcomps);

            auto g_yy_0_yz_xyyyyy = cbuffer.data(di_geom_20_off + 631 * ccomps * dcomps);

            auto g_yy_0_yz_xyyyyz = cbuffer.data(di_geom_20_off + 632 * ccomps * dcomps);

            auto g_yy_0_yz_xyyyzz = cbuffer.data(di_geom_20_off + 633 * ccomps * dcomps);

            auto g_yy_0_yz_xyyzzz = cbuffer.data(di_geom_20_off + 634 * ccomps * dcomps);

            auto g_yy_0_yz_xyzzzz = cbuffer.data(di_geom_20_off + 635 * ccomps * dcomps);

            auto g_yy_0_yz_xzzzzz = cbuffer.data(di_geom_20_off + 636 * ccomps * dcomps);

            auto g_yy_0_yz_yyyyyy = cbuffer.data(di_geom_20_off + 637 * ccomps * dcomps);

            auto g_yy_0_yz_yyyyyz = cbuffer.data(di_geom_20_off + 638 * ccomps * dcomps);

            auto g_yy_0_yz_yyyyzz = cbuffer.data(di_geom_20_off + 639 * ccomps * dcomps);

            auto g_yy_0_yz_yyyzzz = cbuffer.data(di_geom_20_off + 640 * ccomps * dcomps);

            auto g_yy_0_yz_yyzzzz = cbuffer.data(di_geom_20_off + 641 * ccomps * dcomps);

            auto g_yy_0_yz_yzzzzz = cbuffer.data(di_geom_20_off + 642 * ccomps * dcomps);

            auto g_yy_0_yz_zzzzzz = cbuffer.data(di_geom_20_off + 643 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_y_xxxxxx, g_yy_0_y_xxxxxxz, g_yy_0_y_xxxxxy, g_yy_0_y_xxxxxyz, g_yy_0_y_xxxxxz, g_yy_0_y_xxxxxzz, g_yy_0_y_xxxxyy, g_yy_0_y_xxxxyyz, g_yy_0_y_xxxxyz, g_yy_0_y_xxxxyzz, g_yy_0_y_xxxxzz, g_yy_0_y_xxxxzzz, g_yy_0_y_xxxyyy, g_yy_0_y_xxxyyyz, g_yy_0_y_xxxyyz, g_yy_0_y_xxxyyzz, g_yy_0_y_xxxyzz, g_yy_0_y_xxxyzzz, g_yy_0_y_xxxzzz, g_yy_0_y_xxxzzzz, g_yy_0_y_xxyyyy, g_yy_0_y_xxyyyyz, g_yy_0_y_xxyyyz, g_yy_0_y_xxyyyzz, g_yy_0_y_xxyyzz, g_yy_0_y_xxyyzzz, g_yy_0_y_xxyzzz, g_yy_0_y_xxyzzzz, g_yy_0_y_xxzzzz, g_yy_0_y_xxzzzzz, g_yy_0_y_xyyyyy, g_yy_0_y_xyyyyyz, g_yy_0_y_xyyyyz, g_yy_0_y_xyyyyzz, g_yy_0_y_xyyyzz, g_yy_0_y_xyyyzzz, g_yy_0_y_xyyzzz, g_yy_0_y_xyyzzzz, g_yy_0_y_xyzzzz, g_yy_0_y_xyzzzzz, g_yy_0_y_xzzzzz, g_yy_0_y_xzzzzzz, g_yy_0_y_yyyyyy, g_yy_0_y_yyyyyyz, g_yy_0_y_yyyyyz, g_yy_0_y_yyyyyzz, g_yy_0_y_yyyyzz, g_yy_0_y_yyyyzzz, g_yy_0_y_yyyzzz, g_yy_0_y_yyyzzzz, g_yy_0_y_yyzzzz, g_yy_0_y_yyzzzzz, g_yy_0_y_yzzzzz, g_yy_0_y_yzzzzzz, g_yy_0_y_zzzzzz, g_yy_0_y_zzzzzzz, g_yy_0_yz_xxxxxx, g_yy_0_yz_xxxxxy, g_yy_0_yz_xxxxxz, g_yy_0_yz_xxxxyy, g_yy_0_yz_xxxxyz, g_yy_0_yz_xxxxzz, g_yy_0_yz_xxxyyy, g_yy_0_yz_xxxyyz, g_yy_0_yz_xxxyzz, g_yy_0_yz_xxxzzz, g_yy_0_yz_xxyyyy, g_yy_0_yz_xxyyyz, g_yy_0_yz_xxyyzz, g_yy_0_yz_xxyzzz, g_yy_0_yz_xxzzzz, g_yy_0_yz_xyyyyy, g_yy_0_yz_xyyyyz, g_yy_0_yz_xyyyzz, g_yy_0_yz_xyyzzz, g_yy_0_yz_xyzzzz, g_yy_0_yz_xzzzzz, g_yy_0_yz_yyyyyy, g_yy_0_yz_yyyyyz, g_yy_0_yz_yyyyzz, g_yy_0_yz_yyyzzz, g_yy_0_yz_yyzzzz, g_yy_0_yz_yzzzzz, g_yy_0_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yz_xxxxxx[k] = -g_yy_0_y_xxxxxx[k] * ab_z + g_yy_0_y_xxxxxxz[k];

                g_yy_0_yz_xxxxxy[k] = -g_yy_0_y_xxxxxy[k] * ab_z + g_yy_0_y_xxxxxyz[k];

                g_yy_0_yz_xxxxxz[k] = -g_yy_0_y_xxxxxz[k] * ab_z + g_yy_0_y_xxxxxzz[k];

                g_yy_0_yz_xxxxyy[k] = -g_yy_0_y_xxxxyy[k] * ab_z + g_yy_0_y_xxxxyyz[k];

                g_yy_0_yz_xxxxyz[k] = -g_yy_0_y_xxxxyz[k] * ab_z + g_yy_0_y_xxxxyzz[k];

                g_yy_0_yz_xxxxzz[k] = -g_yy_0_y_xxxxzz[k] * ab_z + g_yy_0_y_xxxxzzz[k];

                g_yy_0_yz_xxxyyy[k] = -g_yy_0_y_xxxyyy[k] * ab_z + g_yy_0_y_xxxyyyz[k];

                g_yy_0_yz_xxxyyz[k] = -g_yy_0_y_xxxyyz[k] * ab_z + g_yy_0_y_xxxyyzz[k];

                g_yy_0_yz_xxxyzz[k] = -g_yy_0_y_xxxyzz[k] * ab_z + g_yy_0_y_xxxyzzz[k];

                g_yy_0_yz_xxxzzz[k] = -g_yy_0_y_xxxzzz[k] * ab_z + g_yy_0_y_xxxzzzz[k];

                g_yy_0_yz_xxyyyy[k] = -g_yy_0_y_xxyyyy[k] * ab_z + g_yy_0_y_xxyyyyz[k];

                g_yy_0_yz_xxyyyz[k] = -g_yy_0_y_xxyyyz[k] * ab_z + g_yy_0_y_xxyyyzz[k];

                g_yy_0_yz_xxyyzz[k] = -g_yy_0_y_xxyyzz[k] * ab_z + g_yy_0_y_xxyyzzz[k];

                g_yy_0_yz_xxyzzz[k] = -g_yy_0_y_xxyzzz[k] * ab_z + g_yy_0_y_xxyzzzz[k];

                g_yy_0_yz_xxzzzz[k] = -g_yy_0_y_xxzzzz[k] * ab_z + g_yy_0_y_xxzzzzz[k];

                g_yy_0_yz_xyyyyy[k] = -g_yy_0_y_xyyyyy[k] * ab_z + g_yy_0_y_xyyyyyz[k];

                g_yy_0_yz_xyyyyz[k] = -g_yy_0_y_xyyyyz[k] * ab_z + g_yy_0_y_xyyyyzz[k];

                g_yy_0_yz_xyyyzz[k] = -g_yy_0_y_xyyyzz[k] * ab_z + g_yy_0_y_xyyyzzz[k];

                g_yy_0_yz_xyyzzz[k] = -g_yy_0_y_xyyzzz[k] * ab_z + g_yy_0_y_xyyzzzz[k];

                g_yy_0_yz_xyzzzz[k] = -g_yy_0_y_xyzzzz[k] * ab_z + g_yy_0_y_xyzzzzz[k];

                g_yy_0_yz_xzzzzz[k] = -g_yy_0_y_xzzzzz[k] * ab_z + g_yy_0_y_xzzzzzz[k];

                g_yy_0_yz_yyyyyy[k] = -g_yy_0_y_yyyyyy[k] * ab_z + g_yy_0_y_yyyyyyz[k];

                g_yy_0_yz_yyyyyz[k] = -g_yy_0_y_yyyyyz[k] * ab_z + g_yy_0_y_yyyyyzz[k];

                g_yy_0_yz_yyyyzz[k] = -g_yy_0_y_yyyyzz[k] * ab_z + g_yy_0_y_yyyyzzz[k];

                g_yy_0_yz_yyyzzz[k] = -g_yy_0_y_yyyzzz[k] * ab_z + g_yy_0_y_yyyzzzz[k];

                g_yy_0_yz_yyzzzz[k] = -g_yy_0_y_yyzzzz[k] * ab_z + g_yy_0_y_yyzzzzz[k];

                g_yy_0_yz_yzzzzz[k] = -g_yy_0_y_yzzzzz[k] * ab_z + g_yy_0_y_yzzzzzz[k];

                g_yy_0_yz_zzzzzz[k] = -g_yy_0_y_zzzzzz[k] * ab_z + g_yy_0_y_zzzzzzz[k];
            }

            /// Set up 644-672 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zz_xxxxxx = cbuffer.data(di_geom_20_off + 644 * ccomps * dcomps);

            auto g_yy_0_zz_xxxxxy = cbuffer.data(di_geom_20_off + 645 * ccomps * dcomps);

            auto g_yy_0_zz_xxxxxz = cbuffer.data(di_geom_20_off + 646 * ccomps * dcomps);

            auto g_yy_0_zz_xxxxyy = cbuffer.data(di_geom_20_off + 647 * ccomps * dcomps);

            auto g_yy_0_zz_xxxxyz = cbuffer.data(di_geom_20_off + 648 * ccomps * dcomps);

            auto g_yy_0_zz_xxxxzz = cbuffer.data(di_geom_20_off + 649 * ccomps * dcomps);

            auto g_yy_0_zz_xxxyyy = cbuffer.data(di_geom_20_off + 650 * ccomps * dcomps);

            auto g_yy_0_zz_xxxyyz = cbuffer.data(di_geom_20_off + 651 * ccomps * dcomps);

            auto g_yy_0_zz_xxxyzz = cbuffer.data(di_geom_20_off + 652 * ccomps * dcomps);

            auto g_yy_0_zz_xxxzzz = cbuffer.data(di_geom_20_off + 653 * ccomps * dcomps);

            auto g_yy_0_zz_xxyyyy = cbuffer.data(di_geom_20_off + 654 * ccomps * dcomps);

            auto g_yy_0_zz_xxyyyz = cbuffer.data(di_geom_20_off + 655 * ccomps * dcomps);

            auto g_yy_0_zz_xxyyzz = cbuffer.data(di_geom_20_off + 656 * ccomps * dcomps);

            auto g_yy_0_zz_xxyzzz = cbuffer.data(di_geom_20_off + 657 * ccomps * dcomps);

            auto g_yy_0_zz_xxzzzz = cbuffer.data(di_geom_20_off + 658 * ccomps * dcomps);

            auto g_yy_0_zz_xyyyyy = cbuffer.data(di_geom_20_off + 659 * ccomps * dcomps);

            auto g_yy_0_zz_xyyyyz = cbuffer.data(di_geom_20_off + 660 * ccomps * dcomps);

            auto g_yy_0_zz_xyyyzz = cbuffer.data(di_geom_20_off + 661 * ccomps * dcomps);

            auto g_yy_0_zz_xyyzzz = cbuffer.data(di_geom_20_off + 662 * ccomps * dcomps);

            auto g_yy_0_zz_xyzzzz = cbuffer.data(di_geom_20_off + 663 * ccomps * dcomps);

            auto g_yy_0_zz_xzzzzz = cbuffer.data(di_geom_20_off + 664 * ccomps * dcomps);

            auto g_yy_0_zz_yyyyyy = cbuffer.data(di_geom_20_off + 665 * ccomps * dcomps);

            auto g_yy_0_zz_yyyyyz = cbuffer.data(di_geom_20_off + 666 * ccomps * dcomps);

            auto g_yy_0_zz_yyyyzz = cbuffer.data(di_geom_20_off + 667 * ccomps * dcomps);

            auto g_yy_0_zz_yyyzzz = cbuffer.data(di_geom_20_off + 668 * ccomps * dcomps);

            auto g_yy_0_zz_yyzzzz = cbuffer.data(di_geom_20_off + 669 * ccomps * dcomps);

            auto g_yy_0_zz_yzzzzz = cbuffer.data(di_geom_20_off + 670 * ccomps * dcomps);

            auto g_yy_0_zz_zzzzzz = cbuffer.data(di_geom_20_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_z_xxxxxx, g_yy_0_z_xxxxxxz, g_yy_0_z_xxxxxy, g_yy_0_z_xxxxxyz, g_yy_0_z_xxxxxz, g_yy_0_z_xxxxxzz, g_yy_0_z_xxxxyy, g_yy_0_z_xxxxyyz, g_yy_0_z_xxxxyz, g_yy_0_z_xxxxyzz, g_yy_0_z_xxxxzz, g_yy_0_z_xxxxzzz, g_yy_0_z_xxxyyy, g_yy_0_z_xxxyyyz, g_yy_0_z_xxxyyz, g_yy_0_z_xxxyyzz, g_yy_0_z_xxxyzz, g_yy_0_z_xxxyzzz, g_yy_0_z_xxxzzz, g_yy_0_z_xxxzzzz, g_yy_0_z_xxyyyy, g_yy_0_z_xxyyyyz, g_yy_0_z_xxyyyz, g_yy_0_z_xxyyyzz, g_yy_0_z_xxyyzz, g_yy_0_z_xxyyzzz, g_yy_0_z_xxyzzz, g_yy_0_z_xxyzzzz, g_yy_0_z_xxzzzz, g_yy_0_z_xxzzzzz, g_yy_0_z_xyyyyy, g_yy_0_z_xyyyyyz, g_yy_0_z_xyyyyz, g_yy_0_z_xyyyyzz, g_yy_0_z_xyyyzz, g_yy_0_z_xyyyzzz, g_yy_0_z_xyyzzz, g_yy_0_z_xyyzzzz, g_yy_0_z_xyzzzz, g_yy_0_z_xyzzzzz, g_yy_0_z_xzzzzz, g_yy_0_z_xzzzzzz, g_yy_0_z_yyyyyy, g_yy_0_z_yyyyyyz, g_yy_0_z_yyyyyz, g_yy_0_z_yyyyyzz, g_yy_0_z_yyyyzz, g_yy_0_z_yyyyzzz, g_yy_0_z_yyyzzz, g_yy_0_z_yyyzzzz, g_yy_0_z_yyzzzz, g_yy_0_z_yyzzzzz, g_yy_0_z_yzzzzz, g_yy_0_z_yzzzzzz, g_yy_0_z_zzzzzz, g_yy_0_z_zzzzzzz, g_yy_0_zz_xxxxxx, g_yy_0_zz_xxxxxy, g_yy_0_zz_xxxxxz, g_yy_0_zz_xxxxyy, g_yy_0_zz_xxxxyz, g_yy_0_zz_xxxxzz, g_yy_0_zz_xxxyyy, g_yy_0_zz_xxxyyz, g_yy_0_zz_xxxyzz, g_yy_0_zz_xxxzzz, g_yy_0_zz_xxyyyy, g_yy_0_zz_xxyyyz, g_yy_0_zz_xxyyzz, g_yy_0_zz_xxyzzz, g_yy_0_zz_xxzzzz, g_yy_0_zz_xyyyyy, g_yy_0_zz_xyyyyz, g_yy_0_zz_xyyyzz, g_yy_0_zz_xyyzzz, g_yy_0_zz_xyzzzz, g_yy_0_zz_xzzzzz, g_yy_0_zz_yyyyyy, g_yy_0_zz_yyyyyz, g_yy_0_zz_yyyyzz, g_yy_0_zz_yyyzzz, g_yy_0_zz_yyzzzz, g_yy_0_zz_yzzzzz, g_yy_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zz_xxxxxx[k] = -g_yy_0_z_xxxxxx[k] * ab_z + g_yy_0_z_xxxxxxz[k];

                g_yy_0_zz_xxxxxy[k] = -g_yy_0_z_xxxxxy[k] * ab_z + g_yy_0_z_xxxxxyz[k];

                g_yy_0_zz_xxxxxz[k] = -g_yy_0_z_xxxxxz[k] * ab_z + g_yy_0_z_xxxxxzz[k];

                g_yy_0_zz_xxxxyy[k] = -g_yy_0_z_xxxxyy[k] * ab_z + g_yy_0_z_xxxxyyz[k];

                g_yy_0_zz_xxxxyz[k] = -g_yy_0_z_xxxxyz[k] * ab_z + g_yy_0_z_xxxxyzz[k];

                g_yy_0_zz_xxxxzz[k] = -g_yy_0_z_xxxxzz[k] * ab_z + g_yy_0_z_xxxxzzz[k];

                g_yy_0_zz_xxxyyy[k] = -g_yy_0_z_xxxyyy[k] * ab_z + g_yy_0_z_xxxyyyz[k];

                g_yy_0_zz_xxxyyz[k] = -g_yy_0_z_xxxyyz[k] * ab_z + g_yy_0_z_xxxyyzz[k];

                g_yy_0_zz_xxxyzz[k] = -g_yy_0_z_xxxyzz[k] * ab_z + g_yy_0_z_xxxyzzz[k];

                g_yy_0_zz_xxxzzz[k] = -g_yy_0_z_xxxzzz[k] * ab_z + g_yy_0_z_xxxzzzz[k];

                g_yy_0_zz_xxyyyy[k] = -g_yy_0_z_xxyyyy[k] * ab_z + g_yy_0_z_xxyyyyz[k];

                g_yy_0_zz_xxyyyz[k] = -g_yy_0_z_xxyyyz[k] * ab_z + g_yy_0_z_xxyyyzz[k];

                g_yy_0_zz_xxyyzz[k] = -g_yy_0_z_xxyyzz[k] * ab_z + g_yy_0_z_xxyyzzz[k];

                g_yy_0_zz_xxyzzz[k] = -g_yy_0_z_xxyzzz[k] * ab_z + g_yy_0_z_xxyzzzz[k];

                g_yy_0_zz_xxzzzz[k] = -g_yy_0_z_xxzzzz[k] * ab_z + g_yy_0_z_xxzzzzz[k];

                g_yy_0_zz_xyyyyy[k] = -g_yy_0_z_xyyyyy[k] * ab_z + g_yy_0_z_xyyyyyz[k];

                g_yy_0_zz_xyyyyz[k] = -g_yy_0_z_xyyyyz[k] * ab_z + g_yy_0_z_xyyyyzz[k];

                g_yy_0_zz_xyyyzz[k] = -g_yy_0_z_xyyyzz[k] * ab_z + g_yy_0_z_xyyyzzz[k];

                g_yy_0_zz_xyyzzz[k] = -g_yy_0_z_xyyzzz[k] * ab_z + g_yy_0_z_xyyzzzz[k];

                g_yy_0_zz_xyzzzz[k] = -g_yy_0_z_xyzzzz[k] * ab_z + g_yy_0_z_xyzzzzz[k];

                g_yy_0_zz_xzzzzz[k] = -g_yy_0_z_xzzzzz[k] * ab_z + g_yy_0_z_xzzzzzz[k];

                g_yy_0_zz_yyyyyy[k] = -g_yy_0_z_yyyyyy[k] * ab_z + g_yy_0_z_yyyyyyz[k];

                g_yy_0_zz_yyyyyz[k] = -g_yy_0_z_yyyyyz[k] * ab_z + g_yy_0_z_yyyyyzz[k];

                g_yy_0_zz_yyyyzz[k] = -g_yy_0_z_yyyyzz[k] * ab_z + g_yy_0_z_yyyyzzz[k];

                g_yy_0_zz_yyyzzz[k] = -g_yy_0_z_yyyzzz[k] * ab_z + g_yy_0_z_yyyzzzz[k];

                g_yy_0_zz_yyzzzz[k] = -g_yy_0_z_yyzzzz[k] * ab_z + g_yy_0_z_yyzzzzz[k];

                g_yy_0_zz_yzzzzz[k] = -g_yy_0_z_yzzzzz[k] * ab_z + g_yy_0_z_yzzzzzz[k];

                g_yy_0_zz_zzzzzz[k] = -g_yy_0_z_zzzzzz[k] * ab_z + g_yy_0_z_zzzzzzz[k];
            }

            /// Set up 672-700 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xx_xxxxxx = cbuffer.data(di_geom_20_off + 672 * ccomps * dcomps);

            auto g_yz_0_xx_xxxxxy = cbuffer.data(di_geom_20_off + 673 * ccomps * dcomps);

            auto g_yz_0_xx_xxxxxz = cbuffer.data(di_geom_20_off + 674 * ccomps * dcomps);

            auto g_yz_0_xx_xxxxyy = cbuffer.data(di_geom_20_off + 675 * ccomps * dcomps);

            auto g_yz_0_xx_xxxxyz = cbuffer.data(di_geom_20_off + 676 * ccomps * dcomps);

            auto g_yz_0_xx_xxxxzz = cbuffer.data(di_geom_20_off + 677 * ccomps * dcomps);

            auto g_yz_0_xx_xxxyyy = cbuffer.data(di_geom_20_off + 678 * ccomps * dcomps);

            auto g_yz_0_xx_xxxyyz = cbuffer.data(di_geom_20_off + 679 * ccomps * dcomps);

            auto g_yz_0_xx_xxxyzz = cbuffer.data(di_geom_20_off + 680 * ccomps * dcomps);

            auto g_yz_0_xx_xxxzzz = cbuffer.data(di_geom_20_off + 681 * ccomps * dcomps);

            auto g_yz_0_xx_xxyyyy = cbuffer.data(di_geom_20_off + 682 * ccomps * dcomps);

            auto g_yz_0_xx_xxyyyz = cbuffer.data(di_geom_20_off + 683 * ccomps * dcomps);

            auto g_yz_0_xx_xxyyzz = cbuffer.data(di_geom_20_off + 684 * ccomps * dcomps);

            auto g_yz_0_xx_xxyzzz = cbuffer.data(di_geom_20_off + 685 * ccomps * dcomps);

            auto g_yz_0_xx_xxzzzz = cbuffer.data(di_geom_20_off + 686 * ccomps * dcomps);

            auto g_yz_0_xx_xyyyyy = cbuffer.data(di_geom_20_off + 687 * ccomps * dcomps);

            auto g_yz_0_xx_xyyyyz = cbuffer.data(di_geom_20_off + 688 * ccomps * dcomps);

            auto g_yz_0_xx_xyyyzz = cbuffer.data(di_geom_20_off + 689 * ccomps * dcomps);

            auto g_yz_0_xx_xyyzzz = cbuffer.data(di_geom_20_off + 690 * ccomps * dcomps);

            auto g_yz_0_xx_xyzzzz = cbuffer.data(di_geom_20_off + 691 * ccomps * dcomps);

            auto g_yz_0_xx_xzzzzz = cbuffer.data(di_geom_20_off + 692 * ccomps * dcomps);

            auto g_yz_0_xx_yyyyyy = cbuffer.data(di_geom_20_off + 693 * ccomps * dcomps);

            auto g_yz_0_xx_yyyyyz = cbuffer.data(di_geom_20_off + 694 * ccomps * dcomps);

            auto g_yz_0_xx_yyyyzz = cbuffer.data(di_geom_20_off + 695 * ccomps * dcomps);

            auto g_yz_0_xx_yyyzzz = cbuffer.data(di_geom_20_off + 696 * ccomps * dcomps);

            auto g_yz_0_xx_yyzzzz = cbuffer.data(di_geom_20_off + 697 * ccomps * dcomps);

            auto g_yz_0_xx_yzzzzz = cbuffer.data(di_geom_20_off + 698 * ccomps * dcomps);

            auto g_yz_0_xx_zzzzzz = cbuffer.data(di_geom_20_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_x_xxxxxx, g_yz_0_x_xxxxxxx, g_yz_0_x_xxxxxxy, g_yz_0_x_xxxxxxz, g_yz_0_x_xxxxxy, g_yz_0_x_xxxxxyy, g_yz_0_x_xxxxxyz, g_yz_0_x_xxxxxz, g_yz_0_x_xxxxxzz, g_yz_0_x_xxxxyy, g_yz_0_x_xxxxyyy, g_yz_0_x_xxxxyyz, g_yz_0_x_xxxxyz, g_yz_0_x_xxxxyzz, g_yz_0_x_xxxxzz, g_yz_0_x_xxxxzzz, g_yz_0_x_xxxyyy, g_yz_0_x_xxxyyyy, g_yz_0_x_xxxyyyz, g_yz_0_x_xxxyyz, g_yz_0_x_xxxyyzz, g_yz_0_x_xxxyzz, g_yz_0_x_xxxyzzz, g_yz_0_x_xxxzzz, g_yz_0_x_xxxzzzz, g_yz_0_x_xxyyyy, g_yz_0_x_xxyyyyy, g_yz_0_x_xxyyyyz, g_yz_0_x_xxyyyz, g_yz_0_x_xxyyyzz, g_yz_0_x_xxyyzz, g_yz_0_x_xxyyzzz, g_yz_0_x_xxyzzz, g_yz_0_x_xxyzzzz, g_yz_0_x_xxzzzz, g_yz_0_x_xxzzzzz, g_yz_0_x_xyyyyy, g_yz_0_x_xyyyyyy, g_yz_0_x_xyyyyyz, g_yz_0_x_xyyyyz, g_yz_0_x_xyyyyzz, g_yz_0_x_xyyyzz, g_yz_0_x_xyyyzzz, g_yz_0_x_xyyzzz, g_yz_0_x_xyyzzzz, g_yz_0_x_xyzzzz, g_yz_0_x_xyzzzzz, g_yz_0_x_xzzzzz, g_yz_0_x_xzzzzzz, g_yz_0_x_yyyyyy, g_yz_0_x_yyyyyz, g_yz_0_x_yyyyzz, g_yz_0_x_yyyzzz, g_yz_0_x_yyzzzz, g_yz_0_x_yzzzzz, g_yz_0_x_zzzzzz, g_yz_0_xx_xxxxxx, g_yz_0_xx_xxxxxy, g_yz_0_xx_xxxxxz, g_yz_0_xx_xxxxyy, g_yz_0_xx_xxxxyz, g_yz_0_xx_xxxxzz, g_yz_0_xx_xxxyyy, g_yz_0_xx_xxxyyz, g_yz_0_xx_xxxyzz, g_yz_0_xx_xxxzzz, g_yz_0_xx_xxyyyy, g_yz_0_xx_xxyyyz, g_yz_0_xx_xxyyzz, g_yz_0_xx_xxyzzz, g_yz_0_xx_xxzzzz, g_yz_0_xx_xyyyyy, g_yz_0_xx_xyyyyz, g_yz_0_xx_xyyyzz, g_yz_0_xx_xyyzzz, g_yz_0_xx_xyzzzz, g_yz_0_xx_xzzzzz, g_yz_0_xx_yyyyyy, g_yz_0_xx_yyyyyz, g_yz_0_xx_yyyyzz, g_yz_0_xx_yyyzzz, g_yz_0_xx_yyzzzz, g_yz_0_xx_yzzzzz, g_yz_0_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xx_xxxxxx[k] = -g_yz_0_x_xxxxxx[k] * ab_x + g_yz_0_x_xxxxxxx[k];

                g_yz_0_xx_xxxxxy[k] = -g_yz_0_x_xxxxxy[k] * ab_x + g_yz_0_x_xxxxxxy[k];

                g_yz_0_xx_xxxxxz[k] = -g_yz_0_x_xxxxxz[k] * ab_x + g_yz_0_x_xxxxxxz[k];

                g_yz_0_xx_xxxxyy[k] = -g_yz_0_x_xxxxyy[k] * ab_x + g_yz_0_x_xxxxxyy[k];

                g_yz_0_xx_xxxxyz[k] = -g_yz_0_x_xxxxyz[k] * ab_x + g_yz_0_x_xxxxxyz[k];

                g_yz_0_xx_xxxxzz[k] = -g_yz_0_x_xxxxzz[k] * ab_x + g_yz_0_x_xxxxxzz[k];

                g_yz_0_xx_xxxyyy[k] = -g_yz_0_x_xxxyyy[k] * ab_x + g_yz_0_x_xxxxyyy[k];

                g_yz_0_xx_xxxyyz[k] = -g_yz_0_x_xxxyyz[k] * ab_x + g_yz_0_x_xxxxyyz[k];

                g_yz_0_xx_xxxyzz[k] = -g_yz_0_x_xxxyzz[k] * ab_x + g_yz_0_x_xxxxyzz[k];

                g_yz_0_xx_xxxzzz[k] = -g_yz_0_x_xxxzzz[k] * ab_x + g_yz_0_x_xxxxzzz[k];

                g_yz_0_xx_xxyyyy[k] = -g_yz_0_x_xxyyyy[k] * ab_x + g_yz_0_x_xxxyyyy[k];

                g_yz_0_xx_xxyyyz[k] = -g_yz_0_x_xxyyyz[k] * ab_x + g_yz_0_x_xxxyyyz[k];

                g_yz_0_xx_xxyyzz[k] = -g_yz_0_x_xxyyzz[k] * ab_x + g_yz_0_x_xxxyyzz[k];

                g_yz_0_xx_xxyzzz[k] = -g_yz_0_x_xxyzzz[k] * ab_x + g_yz_0_x_xxxyzzz[k];

                g_yz_0_xx_xxzzzz[k] = -g_yz_0_x_xxzzzz[k] * ab_x + g_yz_0_x_xxxzzzz[k];

                g_yz_0_xx_xyyyyy[k] = -g_yz_0_x_xyyyyy[k] * ab_x + g_yz_0_x_xxyyyyy[k];

                g_yz_0_xx_xyyyyz[k] = -g_yz_0_x_xyyyyz[k] * ab_x + g_yz_0_x_xxyyyyz[k];

                g_yz_0_xx_xyyyzz[k] = -g_yz_0_x_xyyyzz[k] * ab_x + g_yz_0_x_xxyyyzz[k];

                g_yz_0_xx_xyyzzz[k] = -g_yz_0_x_xyyzzz[k] * ab_x + g_yz_0_x_xxyyzzz[k];

                g_yz_0_xx_xyzzzz[k] = -g_yz_0_x_xyzzzz[k] * ab_x + g_yz_0_x_xxyzzzz[k];

                g_yz_0_xx_xzzzzz[k] = -g_yz_0_x_xzzzzz[k] * ab_x + g_yz_0_x_xxzzzzz[k];

                g_yz_0_xx_yyyyyy[k] = -g_yz_0_x_yyyyyy[k] * ab_x + g_yz_0_x_xyyyyyy[k];

                g_yz_0_xx_yyyyyz[k] = -g_yz_0_x_yyyyyz[k] * ab_x + g_yz_0_x_xyyyyyz[k];

                g_yz_0_xx_yyyyzz[k] = -g_yz_0_x_yyyyzz[k] * ab_x + g_yz_0_x_xyyyyzz[k];

                g_yz_0_xx_yyyzzz[k] = -g_yz_0_x_yyyzzz[k] * ab_x + g_yz_0_x_xyyyzzz[k];

                g_yz_0_xx_yyzzzz[k] = -g_yz_0_x_yyzzzz[k] * ab_x + g_yz_0_x_xyyzzzz[k];

                g_yz_0_xx_yzzzzz[k] = -g_yz_0_x_yzzzzz[k] * ab_x + g_yz_0_x_xyzzzzz[k];

                g_yz_0_xx_zzzzzz[k] = -g_yz_0_x_zzzzzz[k] * ab_x + g_yz_0_x_xzzzzzz[k];
            }

            /// Set up 700-728 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xy_xxxxxx = cbuffer.data(di_geom_20_off + 700 * ccomps * dcomps);

            auto g_yz_0_xy_xxxxxy = cbuffer.data(di_geom_20_off + 701 * ccomps * dcomps);

            auto g_yz_0_xy_xxxxxz = cbuffer.data(di_geom_20_off + 702 * ccomps * dcomps);

            auto g_yz_0_xy_xxxxyy = cbuffer.data(di_geom_20_off + 703 * ccomps * dcomps);

            auto g_yz_0_xy_xxxxyz = cbuffer.data(di_geom_20_off + 704 * ccomps * dcomps);

            auto g_yz_0_xy_xxxxzz = cbuffer.data(di_geom_20_off + 705 * ccomps * dcomps);

            auto g_yz_0_xy_xxxyyy = cbuffer.data(di_geom_20_off + 706 * ccomps * dcomps);

            auto g_yz_0_xy_xxxyyz = cbuffer.data(di_geom_20_off + 707 * ccomps * dcomps);

            auto g_yz_0_xy_xxxyzz = cbuffer.data(di_geom_20_off + 708 * ccomps * dcomps);

            auto g_yz_0_xy_xxxzzz = cbuffer.data(di_geom_20_off + 709 * ccomps * dcomps);

            auto g_yz_0_xy_xxyyyy = cbuffer.data(di_geom_20_off + 710 * ccomps * dcomps);

            auto g_yz_0_xy_xxyyyz = cbuffer.data(di_geom_20_off + 711 * ccomps * dcomps);

            auto g_yz_0_xy_xxyyzz = cbuffer.data(di_geom_20_off + 712 * ccomps * dcomps);

            auto g_yz_0_xy_xxyzzz = cbuffer.data(di_geom_20_off + 713 * ccomps * dcomps);

            auto g_yz_0_xy_xxzzzz = cbuffer.data(di_geom_20_off + 714 * ccomps * dcomps);

            auto g_yz_0_xy_xyyyyy = cbuffer.data(di_geom_20_off + 715 * ccomps * dcomps);

            auto g_yz_0_xy_xyyyyz = cbuffer.data(di_geom_20_off + 716 * ccomps * dcomps);

            auto g_yz_0_xy_xyyyzz = cbuffer.data(di_geom_20_off + 717 * ccomps * dcomps);

            auto g_yz_0_xy_xyyzzz = cbuffer.data(di_geom_20_off + 718 * ccomps * dcomps);

            auto g_yz_0_xy_xyzzzz = cbuffer.data(di_geom_20_off + 719 * ccomps * dcomps);

            auto g_yz_0_xy_xzzzzz = cbuffer.data(di_geom_20_off + 720 * ccomps * dcomps);

            auto g_yz_0_xy_yyyyyy = cbuffer.data(di_geom_20_off + 721 * ccomps * dcomps);

            auto g_yz_0_xy_yyyyyz = cbuffer.data(di_geom_20_off + 722 * ccomps * dcomps);

            auto g_yz_0_xy_yyyyzz = cbuffer.data(di_geom_20_off + 723 * ccomps * dcomps);

            auto g_yz_0_xy_yyyzzz = cbuffer.data(di_geom_20_off + 724 * ccomps * dcomps);

            auto g_yz_0_xy_yyzzzz = cbuffer.data(di_geom_20_off + 725 * ccomps * dcomps);

            auto g_yz_0_xy_yzzzzz = cbuffer.data(di_geom_20_off + 726 * ccomps * dcomps);

            auto g_yz_0_xy_zzzzzz = cbuffer.data(di_geom_20_off + 727 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xy_xxxxxx, g_yz_0_xy_xxxxxy, g_yz_0_xy_xxxxxz, g_yz_0_xy_xxxxyy, g_yz_0_xy_xxxxyz, g_yz_0_xy_xxxxzz, g_yz_0_xy_xxxyyy, g_yz_0_xy_xxxyyz, g_yz_0_xy_xxxyzz, g_yz_0_xy_xxxzzz, g_yz_0_xy_xxyyyy, g_yz_0_xy_xxyyyz, g_yz_0_xy_xxyyzz, g_yz_0_xy_xxyzzz, g_yz_0_xy_xxzzzz, g_yz_0_xy_xyyyyy, g_yz_0_xy_xyyyyz, g_yz_0_xy_xyyyzz, g_yz_0_xy_xyyzzz, g_yz_0_xy_xyzzzz, g_yz_0_xy_xzzzzz, g_yz_0_xy_yyyyyy, g_yz_0_xy_yyyyyz, g_yz_0_xy_yyyyzz, g_yz_0_xy_yyyzzz, g_yz_0_xy_yyzzzz, g_yz_0_xy_yzzzzz, g_yz_0_xy_zzzzzz, g_yz_0_y_xxxxxx, g_yz_0_y_xxxxxxx, g_yz_0_y_xxxxxxy, g_yz_0_y_xxxxxxz, g_yz_0_y_xxxxxy, g_yz_0_y_xxxxxyy, g_yz_0_y_xxxxxyz, g_yz_0_y_xxxxxz, g_yz_0_y_xxxxxzz, g_yz_0_y_xxxxyy, g_yz_0_y_xxxxyyy, g_yz_0_y_xxxxyyz, g_yz_0_y_xxxxyz, g_yz_0_y_xxxxyzz, g_yz_0_y_xxxxzz, g_yz_0_y_xxxxzzz, g_yz_0_y_xxxyyy, g_yz_0_y_xxxyyyy, g_yz_0_y_xxxyyyz, g_yz_0_y_xxxyyz, g_yz_0_y_xxxyyzz, g_yz_0_y_xxxyzz, g_yz_0_y_xxxyzzz, g_yz_0_y_xxxzzz, g_yz_0_y_xxxzzzz, g_yz_0_y_xxyyyy, g_yz_0_y_xxyyyyy, g_yz_0_y_xxyyyyz, g_yz_0_y_xxyyyz, g_yz_0_y_xxyyyzz, g_yz_0_y_xxyyzz, g_yz_0_y_xxyyzzz, g_yz_0_y_xxyzzz, g_yz_0_y_xxyzzzz, g_yz_0_y_xxzzzz, g_yz_0_y_xxzzzzz, g_yz_0_y_xyyyyy, g_yz_0_y_xyyyyyy, g_yz_0_y_xyyyyyz, g_yz_0_y_xyyyyz, g_yz_0_y_xyyyyzz, g_yz_0_y_xyyyzz, g_yz_0_y_xyyyzzz, g_yz_0_y_xyyzzz, g_yz_0_y_xyyzzzz, g_yz_0_y_xyzzzz, g_yz_0_y_xyzzzzz, g_yz_0_y_xzzzzz, g_yz_0_y_xzzzzzz, g_yz_0_y_yyyyyy, g_yz_0_y_yyyyyz, g_yz_0_y_yyyyzz, g_yz_0_y_yyyzzz, g_yz_0_y_yyzzzz, g_yz_0_y_yzzzzz, g_yz_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xy_xxxxxx[k] = -g_yz_0_y_xxxxxx[k] * ab_x + g_yz_0_y_xxxxxxx[k];

                g_yz_0_xy_xxxxxy[k] = -g_yz_0_y_xxxxxy[k] * ab_x + g_yz_0_y_xxxxxxy[k];

                g_yz_0_xy_xxxxxz[k] = -g_yz_0_y_xxxxxz[k] * ab_x + g_yz_0_y_xxxxxxz[k];

                g_yz_0_xy_xxxxyy[k] = -g_yz_0_y_xxxxyy[k] * ab_x + g_yz_0_y_xxxxxyy[k];

                g_yz_0_xy_xxxxyz[k] = -g_yz_0_y_xxxxyz[k] * ab_x + g_yz_0_y_xxxxxyz[k];

                g_yz_0_xy_xxxxzz[k] = -g_yz_0_y_xxxxzz[k] * ab_x + g_yz_0_y_xxxxxzz[k];

                g_yz_0_xy_xxxyyy[k] = -g_yz_0_y_xxxyyy[k] * ab_x + g_yz_0_y_xxxxyyy[k];

                g_yz_0_xy_xxxyyz[k] = -g_yz_0_y_xxxyyz[k] * ab_x + g_yz_0_y_xxxxyyz[k];

                g_yz_0_xy_xxxyzz[k] = -g_yz_0_y_xxxyzz[k] * ab_x + g_yz_0_y_xxxxyzz[k];

                g_yz_0_xy_xxxzzz[k] = -g_yz_0_y_xxxzzz[k] * ab_x + g_yz_0_y_xxxxzzz[k];

                g_yz_0_xy_xxyyyy[k] = -g_yz_0_y_xxyyyy[k] * ab_x + g_yz_0_y_xxxyyyy[k];

                g_yz_0_xy_xxyyyz[k] = -g_yz_0_y_xxyyyz[k] * ab_x + g_yz_0_y_xxxyyyz[k];

                g_yz_0_xy_xxyyzz[k] = -g_yz_0_y_xxyyzz[k] * ab_x + g_yz_0_y_xxxyyzz[k];

                g_yz_0_xy_xxyzzz[k] = -g_yz_0_y_xxyzzz[k] * ab_x + g_yz_0_y_xxxyzzz[k];

                g_yz_0_xy_xxzzzz[k] = -g_yz_0_y_xxzzzz[k] * ab_x + g_yz_0_y_xxxzzzz[k];

                g_yz_0_xy_xyyyyy[k] = -g_yz_0_y_xyyyyy[k] * ab_x + g_yz_0_y_xxyyyyy[k];

                g_yz_0_xy_xyyyyz[k] = -g_yz_0_y_xyyyyz[k] * ab_x + g_yz_0_y_xxyyyyz[k];

                g_yz_0_xy_xyyyzz[k] = -g_yz_0_y_xyyyzz[k] * ab_x + g_yz_0_y_xxyyyzz[k];

                g_yz_0_xy_xyyzzz[k] = -g_yz_0_y_xyyzzz[k] * ab_x + g_yz_0_y_xxyyzzz[k];

                g_yz_0_xy_xyzzzz[k] = -g_yz_0_y_xyzzzz[k] * ab_x + g_yz_0_y_xxyzzzz[k];

                g_yz_0_xy_xzzzzz[k] = -g_yz_0_y_xzzzzz[k] * ab_x + g_yz_0_y_xxzzzzz[k];

                g_yz_0_xy_yyyyyy[k] = -g_yz_0_y_yyyyyy[k] * ab_x + g_yz_0_y_xyyyyyy[k];

                g_yz_0_xy_yyyyyz[k] = -g_yz_0_y_yyyyyz[k] * ab_x + g_yz_0_y_xyyyyyz[k];

                g_yz_0_xy_yyyyzz[k] = -g_yz_0_y_yyyyzz[k] * ab_x + g_yz_0_y_xyyyyzz[k];

                g_yz_0_xy_yyyzzz[k] = -g_yz_0_y_yyyzzz[k] * ab_x + g_yz_0_y_xyyyzzz[k];

                g_yz_0_xy_yyzzzz[k] = -g_yz_0_y_yyzzzz[k] * ab_x + g_yz_0_y_xyyzzzz[k];

                g_yz_0_xy_yzzzzz[k] = -g_yz_0_y_yzzzzz[k] * ab_x + g_yz_0_y_xyzzzzz[k];

                g_yz_0_xy_zzzzzz[k] = -g_yz_0_y_zzzzzz[k] * ab_x + g_yz_0_y_xzzzzzz[k];
            }

            /// Set up 728-756 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xz_xxxxxx = cbuffer.data(di_geom_20_off + 728 * ccomps * dcomps);

            auto g_yz_0_xz_xxxxxy = cbuffer.data(di_geom_20_off + 729 * ccomps * dcomps);

            auto g_yz_0_xz_xxxxxz = cbuffer.data(di_geom_20_off + 730 * ccomps * dcomps);

            auto g_yz_0_xz_xxxxyy = cbuffer.data(di_geom_20_off + 731 * ccomps * dcomps);

            auto g_yz_0_xz_xxxxyz = cbuffer.data(di_geom_20_off + 732 * ccomps * dcomps);

            auto g_yz_0_xz_xxxxzz = cbuffer.data(di_geom_20_off + 733 * ccomps * dcomps);

            auto g_yz_0_xz_xxxyyy = cbuffer.data(di_geom_20_off + 734 * ccomps * dcomps);

            auto g_yz_0_xz_xxxyyz = cbuffer.data(di_geom_20_off + 735 * ccomps * dcomps);

            auto g_yz_0_xz_xxxyzz = cbuffer.data(di_geom_20_off + 736 * ccomps * dcomps);

            auto g_yz_0_xz_xxxzzz = cbuffer.data(di_geom_20_off + 737 * ccomps * dcomps);

            auto g_yz_0_xz_xxyyyy = cbuffer.data(di_geom_20_off + 738 * ccomps * dcomps);

            auto g_yz_0_xz_xxyyyz = cbuffer.data(di_geom_20_off + 739 * ccomps * dcomps);

            auto g_yz_0_xz_xxyyzz = cbuffer.data(di_geom_20_off + 740 * ccomps * dcomps);

            auto g_yz_0_xz_xxyzzz = cbuffer.data(di_geom_20_off + 741 * ccomps * dcomps);

            auto g_yz_0_xz_xxzzzz = cbuffer.data(di_geom_20_off + 742 * ccomps * dcomps);

            auto g_yz_0_xz_xyyyyy = cbuffer.data(di_geom_20_off + 743 * ccomps * dcomps);

            auto g_yz_0_xz_xyyyyz = cbuffer.data(di_geom_20_off + 744 * ccomps * dcomps);

            auto g_yz_0_xz_xyyyzz = cbuffer.data(di_geom_20_off + 745 * ccomps * dcomps);

            auto g_yz_0_xz_xyyzzz = cbuffer.data(di_geom_20_off + 746 * ccomps * dcomps);

            auto g_yz_0_xz_xyzzzz = cbuffer.data(di_geom_20_off + 747 * ccomps * dcomps);

            auto g_yz_0_xz_xzzzzz = cbuffer.data(di_geom_20_off + 748 * ccomps * dcomps);

            auto g_yz_0_xz_yyyyyy = cbuffer.data(di_geom_20_off + 749 * ccomps * dcomps);

            auto g_yz_0_xz_yyyyyz = cbuffer.data(di_geom_20_off + 750 * ccomps * dcomps);

            auto g_yz_0_xz_yyyyzz = cbuffer.data(di_geom_20_off + 751 * ccomps * dcomps);

            auto g_yz_0_xz_yyyzzz = cbuffer.data(di_geom_20_off + 752 * ccomps * dcomps);

            auto g_yz_0_xz_yyzzzz = cbuffer.data(di_geom_20_off + 753 * ccomps * dcomps);

            auto g_yz_0_xz_yzzzzz = cbuffer.data(di_geom_20_off + 754 * ccomps * dcomps);

            auto g_yz_0_xz_zzzzzz = cbuffer.data(di_geom_20_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xz_xxxxxx, g_yz_0_xz_xxxxxy, g_yz_0_xz_xxxxxz, g_yz_0_xz_xxxxyy, g_yz_0_xz_xxxxyz, g_yz_0_xz_xxxxzz, g_yz_0_xz_xxxyyy, g_yz_0_xz_xxxyyz, g_yz_0_xz_xxxyzz, g_yz_0_xz_xxxzzz, g_yz_0_xz_xxyyyy, g_yz_0_xz_xxyyyz, g_yz_0_xz_xxyyzz, g_yz_0_xz_xxyzzz, g_yz_0_xz_xxzzzz, g_yz_0_xz_xyyyyy, g_yz_0_xz_xyyyyz, g_yz_0_xz_xyyyzz, g_yz_0_xz_xyyzzz, g_yz_0_xz_xyzzzz, g_yz_0_xz_xzzzzz, g_yz_0_xz_yyyyyy, g_yz_0_xz_yyyyyz, g_yz_0_xz_yyyyzz, g_yz_0_xz_yyyzzz, g_yz_0_xz_yyzzzz, g_yz_0_xz_yzzzzz, g_yz_0_xz_zzzzzz, g_yz_0_z_xxxxxx, g_yz_0_z_xxxxxxx, g_yz_0_z_xxxxxxy, g_yz_0_z_xxxxxxz, g_yz_0_z_xxxxxy, g_yz_0_z_xxxxxyy, g_yz_0_z_xxxxxyz, g_yz_0_z_xxxxxz, g_yz_0_z_xxxxxzz, g_yz_0_z_xxxxyy, g_yz_0_z_xxxxyyy, g_yz_0_z_xxxxyyz, g_yz_0_z_xxxxyz, g_yz_0_z_xxxxyzz, g_yz_0_z_xxxxzz, g_yz_0_z_xxxxzzz, g_yz_0_z_xxxyyy, g_yz_0_z_xxxyyyy, g_yz_0_z_xxxyyyz, g_yz_0_z_xxxyyz, g_yz_0_z_xxxyyzz, g_yz_0_z_xxxyzz, g_yz_0_z_xxxyzzz, g_yz_0_z_xxxzzz, g_yz_0_z_xxxzzzz, g_yz_0_z_xxyyyy, g_yz_0_z_xxyyyyy, g_yz_0_z_xxyyyyz, g_yz_0_z_xxyyyz, g_yz_0_z_xxyyyzz, g_yz_0_z_xxyyzz, g_yz_0_z_xxyyzzz, g_yz_0_z_xxyzzz, g_yz_0_z_xxyzzzz, g_yz_0_z_xxzzzz, g_yz_0_z_xxzzzzz, g_yz_0_z_xyyyyy, g_yz_0_z_xyyyyyy, g_yz_0_z_xyyyyyz, g_yz_0_z_xyyyyz, g_yz_0_z_xyyyyzz, g_yz_0_z_xyyyzz, g_yz_0_z_xyyyzzz, g_yz_0_z_xyyzzz, g_yz_0_z_xyyzzzz, g_yz_0_z_xyzzzz, g_yz_0_z_xyzzzzz, g_yz_0_z_xzzzzz, g_yz_0_z_xzzzzzz, g_yz_0_z_yyyyyy, g_yz_0_z_yyyyyz, g_yz_0_z_yyyyzz, g_yz_0_z_yyyzzz, g_yz_0_z_yyzzzz, g_yz_0_z_yzzzzz, g_yz_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xz_xxxxxx[k] = -g_yz_0_z_xxxxxx[k] * ab_x + g_yz_0_z_xxxxxxx[k];

                g_yz_0_xz_xxxxxy[k] = -g_yz_0_z_xxxxxy[k] * ab_x + g_yz_0_z_xxxxxxy[k];

                g_yz_0_xz_xxxxxz[k] = -g_yz_0_z_xxxxxz[k] * ab_x + g_yz_0_z_xxxxxxz[k];

                g_yz_0_xz_xxxxyy[k] = -g_yz_0_z_xxxxyy[k] * ab_x + g_yz_0_z_xxxxxyy[k];

                g_yz_0_xz_xxxxyz[k] = -g_yz_0_z_xxxxyz[k] * ab_x + g_yz_0_z_xxxxxyz[k];

                g_yz_0_xz_xxxxzz[k] = -g_yz_0_z_xxxxzz[k] * ab_x + g_yz_0_z_xxxxxzz[k];

                g_yz_0_xz_xxxyyy[k] = -g_yz_0_z_xxxyyy[k] * ab_x + g_yz_0_z_xxxxyyy[k];

                g_yz_0_xz_xxxyyz[k] = -g_yz_0_z_xxxyyz[k] * ab_x + g_yz_0_z_xxxxyyz[k];

                g_yz_0_xz_xxxyzz[k] = -g_yz_0_z_xxxyzz[k] * ab_x + g_yz_0_z_xxxxyzz[k];

                g_yz_0_xz_xxxzzz[k] = -g_yz_0_z_xxxzzz[k] * ab_x + g_yz_0_z_xxxxzzz[k];

                g_yz_0_xz_xxyyyy[k] = -g_yz_0_z_xxyyyy[k] * ab_x + g_yz_0_z_xxxyyyy[k];

                g_yz_0_xz_xxyyyz[k] = -g_yz_0_z_xxyyyz[k] * ab_x + g_yz_0_z_xxxyyyz[k];

                g_yz_0_xz_xxyyzz[k] = -g_yz_0_z_xxyyzz[k] * ab_x + g_yz_0_z_xxxyyzz[k];

                g_yz_0_xz_xxyzzz[k] = -g_yz_0_z_xxyzzz[k] * ab_x + g_yz_0_z_xxxyzzz[k];

                g_yz_0_xz_xxzzzz[k] = -g_yz_0_z_xxzzzz[k] * ab_x + g_yz_0_z_xxxzzzz[k];

                g_yz_0_xz_xyyyyy[k] = -g_yz_0_z_xyyyyy[k] * ab_x + g_yz_0_z_xxyyyyy[k];

                g_yz_0_xz_xyyyyz[k] = -g_yz_0_z_xyyyyz[k] * ab_x + g_yz_0_z_xxyyyyz[k];

                g_yz_0_xz_xyyyzz[k] = -g_yz_0_z_xyyyzz[k] * ab_x + g_yz_0_z_xxyyyzz[k];

                g_yz_0_xz_xyyzzz[k] = -g_yz_0_z_xyyzzz[k] * ab_x + g_yz_0_z_xxyyzzz[k];

                g_yz_0_xz_xyzzzz[k] = -g_yz_0_z_xyzzzz[k] * ab_x + g_yz_0_z_xxyzzzz[k];

                g_yz_0_xz_xzzzzz[k] = -g_yz_0_z_xzzzzz[k] * ab_x + g_yz_0_z_xxzzzzz[k];

                g_yz_0_xz_yyyyyy[k] = -g_yz_0_z_yyyyyy[k] * ab_x + g_yz_0_z_xyyyyyy[k];

                g_yz_0_xz_yyyyyz[k] = -g_yz_0_z_yyyyyz[k] * ab_x + g_yz_0_z_xyyyyyz[k];

                g_yz_0_xz_yyyyzz[k] = -g_yz_0_z_yyyyzz[k] * ab_x + g_yz_0_z_xyyyyzz[k];

                g_yz_0_xz_yyyzzz[k] = -g_yz_0_z_yyyzzz[k] * ab_x + g_yz_0_z_xyyyzzz[k];

                g_yz_0_xz_yyzzzz[k] = -g_yz_0_z_yyzzzz[k] * ab_x + g_yz_0_z_xyyzzzz[k];

                g_yz_0_xz_yzzzzz[k] = -g_yz_0_z_yzzzzz[k] * ab_x + g_yz_0_z_xyzzzzz[k];

                g_yz_0_xz_zzzzzz[k] = -g_yz_0_z_zzzzzz[k] * ab_x + g_yz_0_z_xzzzzzz[k];
            }

            /// Set up 756-784 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yy_xxxxxx = cbuffer.data(di_geom_20_off + 756 * ccomps * dcomps);

            auto g_yz_0_yy_xxxxxy = cbuffer.data(di_geom_20_off + 757 * ccomps * dcomps);

            auto g_yz_0_yy_xxxxxz = cbuffer.data(di_geom_20_off + 758 * ccomps * dcomps);

            auto g_yz_0_yy_xxxxyy = cbuffer.data(di_geom_20_off + 759 * ccomps * dcomps);

            auto g_yz_0_yy_xxxxyz = cbuffer.data(di_geom_20_off + 760 * ccomps * dcomps);

            auto g_yz_0_yy_xxxxzz = cbuffer.data(di_geom_20_off + 761 * ccomps * dcomps);

            auto g_yz_0_yy_xxxyyy = cbuffer.data(di_geom_20_off + 762 * ccomps * dcomps);

            auto g_yz_0_yy_xxxyyz = cbuffer.data(di_geom_20_off + 763 * ccomps * dcomps);

            auto g_yz_0_yy_xxxyzz = cbuffer.data(di_geom_20_off + 764 * ccomps * dcomps);

            auto g_yz_0_yy_xxxzzz = cbuffer.data(di_geom_20_off + 765 * ccomps * dcomps);

            auto g_yz_0_yy_xxyyyy = cbuffer.data(di_geom_20_off + 766 * ccomps * dcomps);

            auto g_yz_0_yy_xxyyyz = cbuffer.data(di_geom_20_off + 767 * ccomps * dcomps);

            auto g_yz_0_yy_xxyyzz = cbuffer.data(di_geom_20_off + 768 * ccomps * dcomps);

            auto g_yz_0_yy_xxyzzz = cbuffer.data(di_geom_20_off + 769 * ccomps * dcomps);

            auto g_yz_0_yy_xxzzzz = cbuffer.data(di_geom_20_off + 770 * ccomps * dcomps);

            auto g_yz_0_yy_xyyyyy = cbuffer.data(di_geom_20_off + 771 * ccomps * dcomps);

            auto g_yz_0_yy_xyyyyz = cbuffer.data(di_geom_20_off + 772 * ccomps * dcomps);

            auto g_yz_0_yy_xyyyzz = cbuffer.data(di_geom_20_off + 773 * ccomps * dcomps);

            auto g_yz_0_yy_xyyzzz = cbuffer.data(di_geom_20_off + 774 * ccomps * dcomps);

            auto g_yz_0_yy_xyzzzz = cbuffer.data(di_geom_20_off + 775 * ccomps * dcomps);

            auto g_yz_0_yy_xzzzzz = cbuffer.data(di_geom_20_off + 776 * ccomps * dcomps);

            auto g_yz_0_yy_yyyyyy = cbuffer.data(di_geom_20_off + 777 * ccomps * dcomps);

            auto g_yz_0_yy_yyyyyz = cbuffer.data(di_geom_20_off + 778 * ccomps * dcomps);

            auto g_yz_0_yy_yyyyzz = cbuffer.data(di_geom_20_off + 779 * ccomps * dcomps);

            auto g_yz_0_yy_yyyzzz = cbuffer.data(di_geom_20_off + 780 * ccomps * dcomps);

            auto g_yz_0_yy_yyzzzz = cbuffer.data(di_geom_20_off + 781 * ccomps * dcomps);

            auto g_yz_0_yy_yzzzzz = cbuffer.data(di_geom_20_off + 782 * ccomps * dcomps);

            auto g_yz_0_yy_zzzzzz = cbuffer.data(di_geom_20_off + 783 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_y_xxxxxx, g_yz_0_y_xxxxxxy, g_yz_0_y_xxxxxy, g_yz_0_y_xxxxxyy, g_yz_0_y_xxxxxyz, g_yz_0_y_xxxxxz, g_yz_0_y_xxxxyy, g_yz_0_y_xxxxyyy, g_yz_0_y_xxxxyyz, g_yz_0_y_xxxxyz, g_yz_0_y_xxxxyzz, g_yz_0_y_xxxxzz, g_yz_0_y_xxxyyy, g_yz_0_y_xxxyyyy, g_yz_0_y_xxxyyyz, g_yz_0_y_xxxyyz, g_yz_0_y_xxxyyzz, g_yz_0_y_xxxyzz, g_yz_0_y_xxxyzzz, g_yz_0_y_xxxzzz, g_yz_0_y_xxyyyy, g_yz_0_y_xxyyyyy, g_yz_0_y_xxyyyyz, g_yz_0_y_xxyyyz, g_yz_0_y_xxyyyzz, g_yz_0_y_xxyyzz, g_yz_0_y_xxyyzzz, g_yz_0_y_xxyzzz, g_yz_0_y_xxyzzzz, g_yz_0_y_xxzzzz, g_yz_0_y_xyyyyy, g_yz_0_y_xyyyyyy, g_yz_0_y_xyyyyyz, g_yz_0_y_xyyyyz, g_yz_0_y_xyyyyzz, g_yz_0_y_xyyyzz, g_yz_0_y_xyyyzzz, g_yz_0_y_xyyzzz, g_yz_0_y_xyyzzzz, g_yz_0_y_xyzzzz, g_yz_0_y_xyzzzzz, g_yz_0_y_xzzzzz, g_yz_0_y_yyyyyy, g_yz_0_y_yyyyyyy, g_yz_0_y_yyyyyyz, g_yz_0_y_yyyyyz, g_yz_0_y_yyyyyzz, g_yz_0_y_yyyyzz, g_yz_0_y_yyyyzzz, g_yz_0_y_yyyzzz, g_yz_0_y_yyyzzzz, g_yz_0_y_yyzzzz, g_yz_0_y_yyzzzzz, g_yz_0_y_yzzzzz, g_yz_0_y_yzzzzzz, g_yz_0_y_zzzzzz, g_yz_0_yy_xxxxxx, g_yz_0_yy_xxxxxy, g_yz_0_yy_xxxxxz, g_yz_0_yy_xxxxyy, g_yz_0_yy_xxxxyz, g_yz_0_yy_xxxxzz, g_yz_0_yy_xxxyyy, g_yz_0_yy_xxxyyz, g_yz_0_yy_xxxyzz, g_yz_0_yy_xxxzzz, g_yz_0_yy_xxyyyy, g_yz_0_yy_xxyyyz, g_yz_0_yy_xxyyzz, g_yz_0_yy_xxyzzz, g_yz_0_yy_xxzzzz, g_yz_0_yy_xyyyyy, g_yz_0_yy_xyyyyz, g_yz_0_yy_xyyyzz, g_yz_0_yy_xyyzzz, g_yz_0_yy_xyzzzz, g_yz_0_yy_xzzzzz, g_yz_0_yy_yyyyyy, g_yz_0_yy_yyyyyz, g_yz_0_yy_yyyyzz, g_yz_0_yy_yyyzzz, g_yz_0_yy_yyzzzz, g_yz_0_yy_yzzzzz, g_yz_0_yy_zzzzzz, g_z_0_y_xxxxxx, g_z_0_y_xxxxxy, g_z_0_y_xxxxxz, g_z_0_y_xxxxyy, g_z_0_y_xxxxyz, g_z_0_y_xxxxzz, g_z_0_y_xxxyyy, g_z_0_y_xxxyyz, g_z_0_y_xxxyzz, g_z_0_y_xxxzzz, g_z_0_y_xxyyyy, g_z_0_y_xxyyyz, g_z_0_y_xxyyzz, g_z_0_y_xxyzzz, g_z_0_y_xxzzzz, g_z_0_y_xyyyyy, g_z_0_y_xyyyyz, g_z_0_y_xyyyzz, g_z_0_y_xyyzzz, g_z_0_y_xyzzzz, g_z_0_y_xzzzzz, g_z_0_y_yyyyyy, g_z_0_y_yyyyyz, g_z_0_y_yyyyzz, g_z_0_y_yyyzzz, g_z_0_y_yyzzzz, g_z_0_y_yzzzzz, g_z_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yy_xxxxxx[k] = -g_z_0_y_xxxxxx[k] - g_yz_0_y_xxxxxx[k] * ab_y + g_yz_0_y_xxxxxxy[k];

                g_yz_0_yy_xxxxxy[k] = -g_z_0_y_xxxxxy[k] - g_yz_0_y_xxxxxy[k] * ab_y + g_yz_0_y_xxxxxyy[k];

                g_yz_0_yy_xxxxxz[k] = -g_z_0_y_xxxxxz[k] - g_yz_0_y_xxxxxz[k] * ab_y + g_yz_0_y_xxxxxyz[k];

                g_yz_0_yy_xxxxyy[k] = -g_z_0_y_xxxxyy[k] - g_yz_0_y_xxxxyy[k] * ab_y + g_yz_0_y_xxxxyyy[k];

                g_yz_0_yy_xxxxyz[k] = -g_z_0_y_xxxxyz[k] - g_yz_0_y_xxxxyz[k] * ab_y + g_yz_0_y_xxxxyyz[k];

                g_yz_0_yy_xxxxzz[k] = -g_z_0_y_xxxxzz[k] - g_yz_0_y_xxxxzz[k] * ab_y + g_yz_0_y_xxxxyzz[k];

                g_yz_0_yy_xxxyyy[k] = -g_z_0_y_xxxyyy[k] - g_yz_0_y_xxxyyy[k] * ab_y + g_yz_0_y_xxxyyyy[k];

                g_yz_0_yy_xxxyyz[k] = -g_z_0_y_xxxyyz[k] - g_yz_0_y_xxxyyz[k] * ab_y + g_yz_0_y_xxxyyyz[k];

                g_yz_0_yy_xxxyzz[k] = -g_z_0_y_xxxyzz[k] - g_yz_0_y_xxxyzz[k] * ab_y + g_yz_0_y_xxxyyzz[k];

                g_yz_0_yy_xxxzzz[k] = -g_z_0_y_xxxzzz[k] - g_yz_0_y_xxxzzz[k] * ab_y + g_yz_0_y_xxxyzzz[k];

                g_yz_0_yy_xxyyyy[k] = -g_z_0_y_xxyyyy[k] - g_yz_0_y_xxyyyy[k] * ab_y + g_yz_0_y_xxyyyyy[k];

                g_yz_0_yy_xxyyyz[k] = -g_z_0_y_xxyyyz[k] - g_yz_0_y_xxyyyz[k] * ab_y + g_yz_0_y_xxyyyyz[k];

                g_yz_0_yy_xxyyzz[k] = -g_z_0_y_xxyyzz[k] - g_yz_0_y_xxyyzz[k] * ab_y + g_yz_0_y_xxyyyzz[k];

                g_yz_0_yy_xxyzzz[k] = -g_z_0_y_xxyzzz[k] - g_yz_0_y_xxyzzz[k] * ab_y + g_yz_0_y_xxyyzzz[k];

                g_yz_0_yy_xxzzzz[k] = -g_z_0_y_xxzzzz[k] - g_yz_0_y_xxzzzz[k] * ab_y + g_yz_0_y_xxyzzzz[k];

                g_yz_0_yy_xyyyyy[k] = -g_z_0_y_xyyyyy[k] - g_yz_0_y_xyyyyy[k] * ab_y + g_yz_0_y_xyyyyyy[k];

                g_yz_0_yy_xyyyyz[k] = -g_z_0_y_xyyyyz[k] - g_yz_0_y_xyyyyz[k] * ab_y + g_yz_0_y_xyyyyyz[k];

                g_yz_0_yy_xyyyzz[k] = -g_z_0_y_xyyyzz[k] - g_yz_0_y_xyyyzz[k] * ab_y + g_yz_0_y_xyyyyzz[k];

                g_yz_0_yy_xyyzzz[k] = -g_z_0_y_xyyzzz[k] - g_yz_0_y_xyyzzz[k] * ab_y + g_yz_0_y_xyyyzzz[k];

                g_yz_0_yy_xyzzzz[k] = -g_z_0_y_xyzzzz[k] - g_yz_0_y_xyzzzz[k] * ab_y + g_yz_0_y_xyyzzzz[k];

                g_yz_0_yy_xzzzzz[k] = -g_z_0_y_xzzzzz[k] - g_yz_0_y_xzzzzz[k] * ab_y + g_yz_0_y_xyzzzzz[k];

                g_yz_0_yy_yyyyyy[k] = -g_z_0_y_yyyyyy[k] - g_yz_0_y_yyyyyy[k] * ab_y + g_yz_0_y_yyyyyyy[k];

                g_yz_0_yy_yyyyyz[k] = -g_z_0_y_yyyyyz[k] - g_yz_0_y_yyyyyz[k] * ab_y + g_yz_0_y_yyyyyyz[k];

                g_yz_0_yy_yyyyzz[k] = -g_z_0_y_yyyyzz[k] - g_yz_0_y_yyyyzz[k] * ab_y + g_yz_0_y_yyyyyzz[k];

                g_yz_0_yy_yyyzzz[k] = -g_z_0_y_yyyzzz[k] - g_yz_0_y_yyyzzz[k] * ab_y + g_yz_0_y_yyyyzzz[k];

                g_yz_0_yy_yyzzzz[k] = -g_z_0_y_yyzzzz[k] - g_yz_0_y_yyzzzz[k] * ab_y + g_yz_0_y_yyyzzzz[k];

                g_yz_0_yy_yzzzzz[k] = -g_z_0_y_yzzzzz[k] - g_yz_0_y_yzzzzz[k] * ab_y + g_yz_0_y_yyzzzzz[k];

                g_yz_0_yy_zzzzzz[k] = -g_z_0_y_zzzzzz[k] - g_yz_0_y_zzzzzz[k] * ab_y + g_yz_0_y_yzzzzzz[k];
            }

            /// Set up 784-812 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yz_xxxxxx = cbuffer.data(di_geom_20_off + 784 * ccomps * dcomps);

            auto g_yz_0_yz_xxxxxy = cbuffer.data(di_geom_20_off + 785 * ccomps * dcomps);

            auto g_yz_0_yz_xxxxxz = cbuffer.data(di_geom_20_off + 786 * ccomps * dcomps);

            auto g_yz_0_yz_xxxxyy = cbuffer.data(di_geom_20_off + 787 * ccomps * dcomps);

            auto g_yz_0_yz_xxxxyz = cbuffer.data(di_geom_20_off + 788 * ccomps * dcomps);

            auto g_yz_0_yz_xxxxzz = cbuffer.data(di_geom_20_off + 789 * ccomps * dcomps);

            auto g_yz_0_yz_xxxyyy = cbuffer.data(di_geom_20_off + 790 * ccomps * dcomps);

            auto g_yz_0_yz_xxxyyz = cbuffer.data(di_geom_20_off + 791 * ccomps * dcomps);

            auto g_yz_0_yz_xxxyzz = cbuffer.data(di_geom_20_off + 792 * ccomps * dcomps);

            auto g_yz_0_yz_xxxzzz = cbuffer.data(di_geom_20_off + 793 * ccomps * dcomps);

            auto g_yz_0_yz_xxyyyy = cbuffer.data(di_geom_20_off + 794 * ccomps * dcomps);

            auto g_yz_0_yz_xxyyyz = cbuffer.data(di_geom_20_off + 795 * ccomps * dcomps);

            auto g_yz_0_yz_xxyyzz = cbuffer.data(di_geom_20_off + 796 * ccomps * dcomps);

            auto g_yz_0_yz_xxyzzz = cbuffer.data(di_geom_20_off + 797 * ccomps * dcomps);

            auto g_yz_0_yz_xxzzzz = cbuffer.data(di_geom_20_off + 798 * ccomps * dcomps);

            auto g_yz_0_yz_xyyyyy = cbuffer.data(di_geom_20_off + 799 * ccomps * dcomps);

            auto g_yz_0_yz_xyyyyz = cbuffer.data(di_geom_20_off + 800 * ccomps * dcomps);

            auto g_yz_0_yz_xyyyzz = cbuffer.data(di_geom_20_off + 801 * ccomps * dcomps);

            auto g_yz_0_yz_xyyzzz = cbuffer.data(di_geom_20_off + 802 * ccomps * dcomps);

            auto g_yz_0_yz_xyzzzz = cbuffer.data(di_geom_20_off + 803 * ccomps * dcomps);

            auto g_yz_0_yz_xzzzzz = cbuffer.data(di_geom_20_off + 804 * ccomps * dcomps);

            auto g_yz_0_yz_yyyyyy = cbuffer.data(di_geom_20_off + 805 * ccomps * dcomps);

            auto g_yz_0_yz_yyyyyz = cbuffer.data(di_geom_20_off + 806 * ccomps * dcomps);

            auto g_yz_0_yz_yyyyzz = cbuffer.data(di_geom_20_off + 807 * ccomps * dcomps);

            auto g_yz_0_yz_yyyzzz = cbuffer.data(di_geom_20_off + 808 * ccomps * dcomps);

            auto g_yz_0_yz_yyzzzz = cbuffer.data(di_geom_20_off + 809 * ccomps * dcomps);

            auto g_yz_0_yz_yzzzzz = cbuffer.data(di_geom_20_off + 810 * ccomps * dcomps);

            auto g_yz_0_yz_zzzzzz = cbuffer.data(di_geom_20_off + 811 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yz_xxxxxx, g_yz_0_yz_xxxxxy, g_yz_0_yz_xxxxxz, g_yz_0_yz_xxxxyy, g_yz_0_yz_xxxxyz, g_yz_0_yz_xxxxzz, g_yz_0_yz_xxxyyy, g_yz_0_yz_xxxyyz, g_yz_0_yz_xxxyzz, g_yz_0_yz_xxxzzz, g_yz_0_yz_xxyyyy, g_yz_0_yz_xxyyyz, g_yz_0_yz_xxyyzz, g_yz_0_yz_xxyzzz, g_yz_0_yz_xxzzzz, g_yz_0_yz_xyyyyy, g_yz_0_yz_xyyyyz, g_yz_0_yz_xyyyzz, g_yz_0_yz_xyyzzz, g_yz_0_yz_xyzzzz, g_yz_0_yz_xzzzzz, g_yz_0_yz_yyyyyy, g_yz_0_yz_yyyyyz, g_yz_0_yz_yyyyzz, g_yz_0_yz_yyyzzz, g_yz_0_yz_yyzzzz, g_yz_0_yz_yzzzzz, g_yz_0_yz_zzzzzz, g_yz_0_z_xxxxxx, g_yz_0_z_xxxxxxy, g_yz_0_z_xxxxxy, g_yz_0_z_xxxxxyy, g_yz_0_z_xxxxxyz, g_yz_0_z_xxxxxz, g_yz_0_z_xxxxyy, g_yz_0_z_xxxxyyy, g_yz_0_z_xxxxyyz, g_yz_0_z_xxxxyz, g_yz_0_z_xxxxyzz, g_yz_0_z_xxxxzz, g_yz_0_z_xxxyyy, g_yz_0_z_xxxyyyy, g_yz_0_z_xxxyyyz, g_yz_0_z_xxxyyz, g_yz_0_z_xxxyyzz, g_yz_0_z_xxxyzz, g_yz_0_z_xxxyzzz, g_yz_0_z_xxxzzz, g_yz_0_z_xxyyyy, g_yz_0_z_xxyyyyy, g_yz_0_z_xxyyyyz, g_yz_0_z_xxyyyz, g_yz_0_z_xxyyyzz, g_yz_0_z_xxyyzz, g_yz_0_z_xxyyzzz, g_yz_0_z_xxyzzz, g_yz_0_z_xxyzzzz, g_yz_0_z_xxzzzz, g_yz_0_z_xyyyyy, g_yz_0_z_xyyyyyy, g_yz_0_z_xyyyyyz, g_yz_0_z_xyyyyz, g_yz_0_z_xyyyyzz, g_yz_0_z_xyyyzz, g_yz_0_z_xyyyzzz, g_yz_0_z_xyyzzz, g_yz_0_z_xyyzzzz, g_yz_0_z_xyzzzz, g_yz_0_z_xyzzzzz, g_yz_0_z_xzzzzz, g_yz_0_z_yyyyyy, g_yz_0_z_yyyyyyy, g_yz_0_z_yyyyyyz, g_yz_0_z_yyyyyz, g_yz_0_z_yyyyyzz, g_yz_0_z_yyyyzz, g_yz_0_z_yyyyzzz, g_yz_0_z_yyyzzz, g_yz_0_z_yyyzzzz, g_yz_0_z_yyzzzz, g_yz_0_z_yyzzzzz, g_yz_0_z_yzzzzz, g_yz_0_z_yzzzzzz, g_yz_0_z_zzzzzz, g_z_0_z_xxxxxx, g_z_0_z_xxxxxy, g_z_0_z_xxxxxz, g_z_0_z_xxxxyy, g_z_0_z_xxxxyz, g_z_0_z_xxxxzz, g_z_0_z_xxxyyy, g_z_0_z_xxxyyz, g_z_0_z_xxxyzz, g_z_0_z_xxxzzz, g_z_0_z_xxyyyy, g_z_0_z_xxyyyz, g_z_0_z_xxyyzz, g_z_0_z_xxyzzz, g_z_0_z_xxzzzz, g_z_0_z_xyyyyy, g_z_0_z_xyyyyz, g_z_0_z_xyyyzz, g_z_0_z_xyyzzz, g_z_0_z_xyzzzz, g_z_0_z_xzzzzz, g_z_0_z_yyyyyy, g_z_0_z_yyyyyz, g_z_0_z_yyyyzz, g_z_0_z_yyyzzz, g_z_0_z_yyzzzz, g_z_0_z_yzzzzz, g_z_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yz_xxxxxx[k] = -g_z_0_z_xxxxxx[k] - g_yz_0_z_xxxxxx[k] * ab_y + g_yz_0_z_xxxxxxy[k];

                g_yz_0_yz_xxxxxy[k] = -g_z_0_z_xxxxxy[k] - g_yz_0_z_xxxxxy[k] * ab_y + g_yz_0_z_xxxxxyy[k];

                g_yz_0_yz_xxxxxz[k] = -g_z_0_z_xxxxxz[k] - g_yz_0_z_xxxxxz[k] * ab_y + g_yz_0_z_xxxxxyz[k];

                g_yz_0_yz_xxxxyy[k] = -g_z_0_z_xxxxyy[k] - g_yz_0_z_xxxxyy[k] * ab_y + g_yz_0_z_xxxxyyy[k];

                g_yz_0_yz_xxxxyz[k] = -g_z_0_z_xxxxyz[k] - g_yz_0_z_xxxxyz[k] * ab_y + g_yz_0_z_xxxxyyz[k];

                g_yz_0_yz_xxxxzz[k] = -g_z_0_z_xxxxzz[k] - g_yz_0_z_xxxxzz[k] * ab_y + g_yz_0_z_xxxxyzz[k];

                g_yz_0_yz_xxxyyy[k] = -g_z_0_z_xxxyyy[k] - g_yz_0_z_xxxyyy[k] * ab_y + g_yz_0_z_xxxyyyy[k];

                g_yz_0_yz_xxxyyz[k] = -g_z_0_z_xxxyyz[k] - g_yz_0_z_xxxyyz[k] * ab_y + g_yz_0_z_xxxyyyz[k];

                g_yz_0_yz_xxxyzz[k] = -g_z_0_z_xxxyzz[k] - g_yz_0_z_xxxyzz[k] * ab_y + g_yz_0_z_xxxyyzz[k];

                g_yz_0_yz_xxxzzz[k] = -g_z_0_z_xxxzzz[k] - g_yz_0_z_xxxzzz[k] * ab_y + g_yz_0_z_xxxyzzz[k];

                g_yz_0_yz_xxyyyy[k] = -g_z_0_z_xxyyyy[k] - g_yz_0_z_xxyyyy[k] * ab_y + g_yz_0_z_xxyyyyy[k];

                g_yz_0_yz_xxyyyz[k] = -g_z_0_z_xxyyyz[k] - g_yz_0_z_xxyyyz[k] * ab_y + g_yz_0_z_xxyyyyz[k];

                g_yz_0_yz_xxyyzz[k] = -g_z_0_z_xxyyzz[k] - g_yz_0_z_xxyyzz[k] * ab_y + g_yz_0_z_xxyyyzz[k];

                g_yz_0_yz_xxyzzz[k] = -g_z_0_z_xxyzzz[k] - g_yz_0_z_xxyzzz[k] * ab_y + g_yz_0_z_xxyyzzz[k];

                g_yz_0_yz_xxzzzz[k] = -g_z_0_z_xxzzzz[k] - g_yz_0_z_xxzzzz[k] * ab_y + g_yz_0_z_xxyzzzz[k];

                g_yz_0_yz_xyyyyy[k] = -g_z_0_z_xyyyyy[k] - g_yz_0_z_xyyyyy[k] * ab_y + g_yz_0_z_xyyyyyy[k];

                g_yz_0_yz_xyyyyz[k] = -g_z_0_z_xyyyyz[k] - g_yz_0_z_xyyyyz[k] * ab_y + g_yz_0_z_xyyyyyz[k];

                g_yz_0_yz_xyyyzz[k] = -g_z_0_z_xyyyzz[k] - g_yz_0_z_xyyyzz[k] * ab_y + g_yz_0_z_xyyyyzz[k];

                g_yz_0_yz_xyyzzz[k] = -g_z_0_z_xyyzzz[k] - g_yz_0_z_xyyzzz[k] * ab_y + g_yz_0_z_xyyyzzz[k];

                g_yz_0_yz_xyzzzz[k] = -g_z_0_z_xyzzzz[k] - g_yz_0_z_xyzzzz[k] * ab_y + g_yz_0_z_xyyzzzz[k];

                g_yz_0_yz_xzzzzz[k] = -g_z_0_z_xzzzzz[k] - g_yz_0_z_xzzzzz[k] * ab_y + g_yz_0_z_xyzzzzz[k];

                g_yz_0_yz_yyyyyy[k] = -g_z_0_z_yyyyyy[k] - g_yz_0_z_yyyyyy[k] * ab_y + g_yz_0_z_yyyyyyy[k];

                g_yz_0_yz_yyyyyz[k] = -g_z_0_z_yyyyyz[k] - g_yz_0_z_yyyyyz[k] * ab_y + g_yz_0_z_yyyyyyz[k];

                g_yz_0_yz_yyyyzz[k] = -g_z_0_z_yyyyzz[k] - g_yz_0_z_yyyyzz[k] * ab_y + g_yz_0_z_yyyyyzz[k];

                g_yz_0_yz_yyyzzz[k] = -g_z_0_z_yyyzzz[k] - g_yz_0_z_yyyzzz[k] * ab_y + g_yz_0_z_yyyyzzz[k];

                g_yz_0_yz_yyzzzz[k] = -g_z_0_z_yyzzzz[k] - g_yz_0_z_yyzzzz[k] * ab_y + g_yz_0_z_yyyzzzz[k];

                g_yz_0_yz_yzzzzz[k] = -g_z_0_z_yzzzzz[k] - g_yz_0_z_yzzzzz[k] * ab_y + g_yz_0_z_yyzzzzz[k];

                g_yz_0_yz_zzzzzz[k] = -g_z_0_z_zzzzzz[k] - g_yz_0_z_zzzzzz[k] * ab_y + g_yz_0_z_yzzzzzz[k];
            }

            /// Set up 812-840 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zz_xxxxxx = cbuffer.data(di_geom_20_off + 812 * ccomps * dcomps);

            auto g_yz_0_zz_xxxxxy = cbuffer.data(di_geom_20_off + 813 * ccomps * dcomps);

            auto g_yz_0_zz_xxxxxz = cbuffer.data(di_geom_20_off + 814 * ccomps * dcomps);

            auto g_yz_0_zz_xxxxyy = cbuffer.data(di_geom_20_off + 815 * ccomps * dcomps);

            auto g_yz_0_zz_xxxxyz = cbuffer.data(di_geom_20_off + 816 * ccomps * dcomps);

            auto g_yz_0_zz_xxxxzz = cbuffer.data(di_geom_20_off + 817 * ccomps * dcomps);

            auto g_yz_0_zz_xxxyyy = cbuffer.data(di_geom_20_off + 818 * ccomps * dcomps);

            auto g_yz_0_zz_xxxyyz = cbuffer.data(di_geom_20_off + 819 * ccomps * dcomps);

            auto g_yz_0_zz_xxxyzz = cbuffer.data(di_geom_20_off + 820 * ccomps * dcomps);

            auto g_yz_0_zz_xxxzzz = cbuffer.data(di_geom_20_off + 821 * ccomps * dcomps);

            auto g_yz_0_zz_xxyyyy = cbuffer.data(di_geom_20_off + 822 * ccomps * dcomps);

            auto g_yz_0_zz_xxyyyz = cbuffer.data(di_geom_20_off + 823 * ccomps * dcomps);

            auto g_yz_0_zz_xxyyzz = cbuffer.data(di_geom_20_off + 824 * ccomps * dcomps);

            auto g_yz_0_zz_xxyzzz = cbuffer.data(di_geom_20_off + 825 * ccomps * dcomps);

            auto g_yz_0_zz_xxzzzz = cbuffer.data(di_geom_20_off + 826 * ccomps * dcomps);

            auto g_yz_0_zz_xyyyyy = cbuffer.data(di_geom_20_off + 827 * ccomps * dcomps);

            auto g_yz_0_zz_xyyyyz = cbuffer.data(di_geom_20_off + 828 * ccomps * dcomps);

            auto g_yz_0_zz_xyyyzz = cbuffer.data(di_geom_20_off + 829 * ccomps * dcomps);

            auto g_yz_0_zz_xyyzzz = cbuffer.data(di_geom_20_off + 830 * ccomps * dcomps);

            auto g_yz_0_zz_xyzzzz = cbuffer.data(di_geom_20_off + 831 * ccomps * dcomps);

            auto g_yz_0_zz_xzzzzz = cbuffer.data(di_geom_20_off + 832 * ccomps * dcomps);

            auto g_yz_0_zz_yyyyyy = cbuffer.data(di_geom_20_off + 833 * ccomps * dcomps);

            auto g_yz_0_zz_yyyyyz = cbuffer.data(di_geom_20_off + 834 * ccomps * dcomps);

            auto g_yz_0_zz_yyyyzz = cbuffer.data(di_geom_20_off + 835 * ccomps * dcomps);

            auto g_yz_0_zz_yyyzzz = cbuffer.data(di_geom_20_off + 836 * ccomps * dcomps);

            auto g_yz_0_zz_yyzzzz = cbuffer.data(di_geom_20_off + 837 * ccomps * dcomps);

            auto g_yz_0_zz_yzzzzz = cbuffer.data(di_geom_20_off + 838 * ccomps * dcomps);

            auto g_yz_0_zz_zzzzzz = cbuffer.data(di_geom_20_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxxxxx, g_y_0_z_xxxxxy, g_y_0_z_xxxxxz, g_y_0_z_xxxxyy, g_y_0_z_xxxxyz, g_y_0_z_xxxxzz, g_y_0_z_xxxyyy, g_y_0_z_xxxyyz, g_y_0_z_xxxyzz, g_y_0_z_xxxzzz, g_y_0_z_xxyyyy, g_y_0_z_xxyyyz, g_y_0_z_xxyyzz, g_y_0_z_xxyzzz, g_y_0_z_xxzzzz, g_y_0_z_xyyyyy, g_y_0_z_xyyyyz, g_y_0_z_xyyyzz, g_y_0_z_xyyzzz, g_y_0_z_xyzzzz, g_y_0_z_xzzzzz, g_y_0_z_yyyyyy, g_y_0_z_yyyyyz, g_y_0_z_yyyyzz, g_y_0_z_yyyzzz, g_y_0_z_yyzzzz, g_y_0_z_yzzzzz, g_y_0_z_zzzzzz, g_yz_0_z_xxxxxx, g_yz_0_z_xxxxxxz, g_yz_0_z_xxxxxy, g_yz_0_z_xxxxxyz, g_yz_0_z_xxxxxz, g_yz_0_z_xxxxxzz, g_yz_0_z_xxxxyy, g_yz_0_z_xxxxyyz, g_yz_0_z_xxxxyz, g_yz_0_z_xxxxyzz, g_yz_0_z_xxxxzz, g_yz_0_z_xxxxzzz, g_yz_0_z_xxxyyy, g_yz_0_z_xxxyyyz, g_yz_0_z_xxxyyz, g_yz_0_z_xxxyyzz, g_yz_0_z_xxxyzz, g_yz_0_z_xxxyzzz, g_yz_0_z_xxxzzz, g_yz_0_z_xxxzzzz, g_yz_0_z_xxyyyy, g_yz_0_z_xxyyyyz, g_yz_0_z_xxyyyz, g_yz_0_z_xxyyyzz, g_yz_0_z_xxyyzz, g_yz_0_z_xxyyzzz, g_yz_0_z_xxyzzz, g_yz_0_z_xxyzzzz, g_yz_0_z_xxzzzz, g_yz_0_z_xxzzzzz, g_yz_0_z_xyyyyy, g_yz_0_z_xyyyyyz, g_yz_0_z_xyyyyz, g_yz_0_z_xyyyyzz, g_yz_0_z_xyyyzz, g_yz_0_z_xyyyzzz, g_yz_0_z_xyyzzz, g_yz_0_z_xyyzzzz, g_yz_0_z_xyzzzz, g_yz_0_z_xyzzzzz, g_yz_0_z_xzzzzz, g_yz_0_z_xzzzzzz, g_yz_0_z_yyyyyy, g_yz_0_z_yyyyyyz, g_yz_0_z_yyyyyz, g_yz_0_z_yyyyyzz, g_yz_0_z_yyyyzz, g_yz_0_z_yyyyzzz, g_yz_0_z_yyyzzz, g_yz_0_z_yyyzzzz, g_yz_0_z_yyzzzz, g_yz_0_z_yyzzzzz, g_yz_0_z_yzzzzz, g_yz_0_z_yzzzzzz, g_yz_0_z_zzzzzz, g_yz_0_z_zzzzzzz, g_yz_0_zz_xxxxxx, g_yz_0_zz_xxxxxy, g_yz_0_zz_xxxxxz, g_yz_0_zz_xxxxyy, g_yz_0_zz_xxxxyz, g_yz_0_zz_xxxxzz, g_yz_0_zz_xxxyyy, g_yz_0_zz_xxxyyz, g_yz_0_zz_xxxyzz, g_yz_0_zz_xxxzzz, g_yz_0_zz_xxyyyy, g_yz_0_zz_xxyyyz, g_yz_0_zz_xxyyzz, g_yz_0_zz_xxyzzz, g_yz_0_zz_xxzzzz, g_yz_0_zz_xyyyyy, g_yz_0_zz_xyyyyz, g_yz_0_zz_xyyyzz, g_yz_0_zz_xyyzzz, g_yz_0_zz_xyzzzz, g_yz_0_zz_xzzzzz, g_yz_0_zz_yyyyyy, g_yz_0_zz_yyyyyz, g_yz_0_zz_yyyyzz, g_yz_0_zz_yyyzzz, g_yz_0_zz_yyzzzz, g_yz_0_zz_yzzzzz, g_yz_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zz_xxxxxx[k] = -g_y_0_z_xxxxxx[k] - g_yz_0_z_xxxxxx[k] * ab_z + g_yz_0_z_xxxxxxz[k];

                g_yz_0_zz_xxxxxy[k] = -g_y_0_z_xxxxxy[k] - g_yz_0_z_xxxxxy[k] * ab_z + g_yz_0_z_xxxxxyz[k];

                g_yz_0_zz_xxxxxz[k] = -g_y_0_z_xxxxxz[k] - g_yz_0_z_xxxxxz[k] * ab_z + g_yz_0_z_xxxxxzz[k];

                g_yz_0_zz_xxxxyy[k] = -g_y_0_z_xxxxyy[k] - g_yz_0_z_xxxxyy[k] * ab_z + g_yz_0_z_xxxxyyz[k];

                g_yz_0_zz_xxxxyz[k] = -g_y_0_z_xxxxyz[k] - g_yz_0_z_xxxxyz[k] * ab_z + g_yz_0_z_xxxxyzz[k];

                g_yz_0_zz_xxxxzz[k] = -g_y_0_z_xxxxzz[k] - g_yz_0_z_xxxxzz[k] * ab_z + g_yz_0_z_xxxxzzz[k];

                g_yz_0_zz_xxxyyy[k] = -g_y_0_z_xxxyyy[k] - g_yz_0_z_xxxyyy[k] * ab_z + g_yz_0_z_xxxyyyz[k];

                g_yz_0_zz_xxxyyz[k] = -g_y_0_z_xxxyyz[k] - g_yz_0_z_xxxyyz[k] * ab_z + g_yz_0_z_xxxyyzz[k];

                g_yz_0_zz_xxxyzz[k] = -g_y_0_z_xxxyzz[k] - g_yz_0_z_xxxyzz[k] * ab_z + g_yz_0_z_xxxyzzz[k];

                g_yz_0_zz_xxxzzz[k] = -g_y_0_z_xxxzzz[k] - g_yz_0_z_xxxzzz[k] * ab_z + g_yz_0_z_xxxzzzz[k];

                g_yz_0_zz_xxyyyy[k] = -g_y_0_z_xxyyyy[k] - g_yz_0_z_xxyyyy[k] * ab_z + g_yz_0_z_xxyyyyz[k];

                g_yz_0_zz_xxyyyz[k] = -g_y_0_z_xxyyyz[k] - g_yz_0_z_xxyyyz[k] * ab_z + g_yz_0_z_xxyyyzz[k];

                g_yz_0_zz_xxyyzz[k] = -g_y_0_z_xxyyzz[k] - g_yz_0_z_xxyyzz[k] * ab_z + g_yz_0_z_xxyyzzz[k];

                g_yz_0_zz_xxyzzz[k] = -g_y_0_z_xxyzzz[k] - g_yz_0_z_xxyzzz[k] * ab_z + g_yz_0_z_xxyzzzz[k];

                g_yz_0_zz_xxzzzz[k] = -g_y_0_z_xxzzzz[k] - g_yz_0_z_xxzzzz[k] * ab_z + g_yz_0_z_xxzzzzz[k];

                g_yz_0_zz_xyyyyy[k] = -g_y_0_z_xyyyyy[k] - g_yz_0_z_xyyyyy[k] * ab_z + g_yz_0_z_xyyyyyz[k];

                g_yz_0_zz_xyyyyz[k] = -g_y_0_z_xyyyyz[k] - g_yz_0_z_xyyyyz[k] * ab_z + g_yz_0_z_xyyyyzz[k];

                g_yz_0_zz_xyyyzz[k] = -g_y_0_z_xyyyzz[k] - g_yz_0_z_xyyyzz[k] * ab_z + g_yz_0_z_xyyyzzz[k];

                g_yz_0_zz_xyyzzz[k] = -g_y_0_z_xyyzzz[k] - g_yz_0_z_xyyzzz[k] * ab_z + g_yz_0_z_xyyzzzz[k];

                g_yz_0_zz_xyzzzz[k] = -g_y_0_z_xyzzzz[k] - g_yz_0_z_xyzzzz[k] * ab_z + g_yz_0_z_xyzzzzz[k];

                g_yz_0_zz_xzzzzz[k] = -g_y_0_z_xzzzzz[k] - g_yz_0_z_xzzzzz[k] * ab_z + g_yz_0_z_xzzzzzz[k];

                g_yz_0_zz_yyyyyy[k] = -g_y_0_z_yyyyyy[k] - g_yz_0_z_yyyyyy[k] * ab_z + g_yz_0_z_yyyyyyz[k];

                g_yz_0_zz_yyyyyz[k] = -g_y_0_z_yyyyyz[k] - g_yz_0_z_yyyyyz[k] * ab_z + g_yz_0_z_yyyyyzz[k];

                g_yz_0_zz_yyyyzz[k] = -g_y_0_z_yyyyzz[k] - g_yz_0_z_yyyyzz[k] * ab_z + g_yz_0_z_yyyyzzz[k];

                g_yz_0_zz_yyyzzz[k] = -g_y_0_z_yyyzzz[k] - g_yz_0_z_yyyzzz[k] * ab_z + g_yz_0_z_yyyzzzz[k];

                g_yz_0_zz_yyzzzz[k] = -g_y_0_z_yyzzzz[k] - g_yz_0_z_yyzzzz[k] * ab_z + g_yz_0_z_yyzzzzz[k];

                g_yz_0_zz_yzzzzz[k] = -g_y_0_z_yzzzzz[k] - g_yz_0_z_yzzzzz[k] * ab_z + g_yz_0_z_yzzzzzz[k];

                g_yz_0_zz_zzzzzz[k] = -g_y_0_z_zzzzzz[k] - g_yz_0_z_zzzzzz[k] * ab_z + g_yz_0_z_zzzzzzz[k];
            }

            /// Set up 840-868 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xx_xxxxxx = cbuffer.data(di_geom_20_off + 840 * ccomps * dcomps);

            auto g_zz_0_xx_xxxxxy = cbuffer.data(di_geom_20_off + 841 * ccomps * dcomps);

            auto g_zz_0_xx_xxxxxz = cbuffer.data(di_geom_20_off + 842 * ccomps * dcomps);

            auto g_zz_0_xx_xxxxyy = cbuffer.data(di_geom_20_off + 843 * ccomps * dcomps);

            auto g_zz_0_xx_xxxxyz = cbuffer.data(di_geom_20_off + 844 * ccomps * dcomps);

            auto g_zz_0_xx_xxxxzz = cbuffer.data(di_geom_20_off + 845 * ccomps * dcomps);

            auto g_zz_0_xx_xxxyyy = cbuffer.data(di_geom_20_off + 846 * ccomps * dcomps);

            auto g_zz_0_xx_xxxyyz = cbuffer.data(di_geom_20_off + 847 * ccomps * dcomps);

            auto g_zz_0_xx_xxxyzz = cbuffer.data(di_geom_20_off + 848 * ccomps * dcomps);

            auto g_zz_0_xx_xxxzzz = cbuffer.data(di_geom_20_off + 849 * ccomps * dcomps);

            auto g_zz_0_xx_xxyyyy = cbuffer.data(di_geom_20_off + 850 * ccomps * dcomps);

            auto g_zz_0_xx_xxyyyz = cbuffer.data(di_geom_20_off + 851 * ccomps * dcomps);

            auto g_zz_0_xx_xxyyzz = cbuffer.data(di_geom_20_off + 852 * ccomps * dcomps);

            auto g_zz_0_xx_xxyzzz = cbuffer.data(di_geom_20_off + 853 * ccomps * dcomps);

            auto g_zz_0_xx_xxzzzz = cbuffer.data(di_geom_20_off + 854 * ccomps * dcomps);

            auto g_zz_0_xx_xyyyyy = cbuffer.data(di_geom_20_off + 855 * ccomps * dcomps);

            auto g_zz_0_xx_xyyyyz = cbuffer.data(di_geom_20_off + 856 * ccomps * dcomps);

            auto g_zz_0_xx_xyyyzz = cbuffer.data(di_geom_20_off + 857 * ccomps * dcomps);

            auto g_zz_0_xx_xyyzzz = cbuffer.data(di_geom_20_off + 858 * ccomps * dcomps);

            auto g_zz_0_xx_xyzzzz = cbuffer.data(di_geom_20_off + 859 * ccomps * dcomps);

            auto g_zz_0_xx_xzzzzz = cbuffer.data(di_geom_20_off + 860 * ccomps * dcomps);

            auto g_zz_0_xx_yyyyyy = cbuffer.data(di_geom_20_off + 861 * ccomps * dcomps);

            auto g_zz_0_xx_yyyyyz = cbuffer.data(di_geom_20_off + 862 * ccomps * dcomps);

            auto g_zz_0_xx_yyyyzz = cbuffer.data(di_geom_20_off + 863 * ccomps * dcomps);

            auto g_zz_0_xx_yyyzzz = cbuffer.data(di_geom_20_off + 864 * ccomps * dcomps);

            auto g_zz_0_xx_yyzzzz = cbuffer.data(di_geom_20_off + 865 * ccomps * dcomps);

            auto g_zz_0_xx_yzzzzz = cbuffer.data(di_geom_20_off + 866 * ccomps * dcomps);

            auto g_zz_0_xx_zzzzzz = cbuffer.data(di_geom_20_off + 867 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_x_xxxxxx, g_zz_0_x_xxxxxxx, g_zz_0_x_xxxxxxy, g_zz_0_x_xxxxxxz, g_zz_0_x_xxxxxy, g_zz_0_x_xxxxxyy, g_zz_0_x_xxxxxyz, g_zz_0_x_xxxxxz, g_zz_0_x_xxxxxzz, g_zz_0_x_xxxxyy, g_zz_0_x_xxxxyyy, g_zz_0_x_xxxxyyz, g_zz_0_x_xxxxyz, g_zz_0_x_xxxxyzz, g_zz_0_x_xxxxzz, g_zz_0_x_xxxxzzz, g_zz_0_x_xxxyyy, g_zz_0_x_xxxyyyy, g_zz_0_x_xxxyyyz, g_zz_0_x_xxxyyz, g_zz_0_x_xxxyyzz, g_zz_0_x_xxxyzz, g_zz_0_x_xxxyzzz, g_zz_0_x_xxxzzz, g_zz_0_x_xxxzzzz, g_zz_0_x_xxyyyy, g_zz_0_x_xxyyyyy, g_zz_0_x_xxyyyyz, g_zz_0_x_xxyyyz, g_zz_0_x_xxyyyzz, g_zz_0_x_xxyyzz, g_zz_0_x_xxyyzzz, g_zz_0_x_xxyzzz, g_zz_0_x_xxyzzzz, g_zz_0_x_xxzzzz, g_zz_0_x_xxzzzzz, g_zz_0_x_xyyyyy, g_zz_0_x_xyyyyyy, g_zz_0_x_xyyyyyz, g_zz_0_x_xyyyyz, g_zz_0_x_xyyyyzz, g_zz_0_x_xyyyzz, g_zz_0_x_xyyyzzz, g_zz_0_x_xyyzzz, g_zz_0_x_xyyzzzz, g_zz_0_x_xyzzzz, g_zz_0_x_xyzzzzz, g_zz_0_x_xzzzzz, g_zz_0_x_xzzzzzz, g_zz_0_x_yyyyyy, g_zz_0_x_yyyyyz, g_zz_0_x_yyyyzz, g_zz_0_x_yyyzzz, g_zz_0_x_yyzzzz, g_zz_0_x_yzzzzz, g_zz_0_x_zzzzzz, g_zz_0_xx_xxxxxx, g_zz_0_xx_xxxxxy, g_zz_0_xx_xxxxxz, g_zz_0_xx_xxxxyy, g_zz_0_xx_xxxxyz, g_zz_0_xx_xxxxzz, g_zz_0_xx_xxxyyy, g_zz_0_xx_xxxyyz, g_zz_0_xx_xxxyzz, g_zz_0_xx_xxxzzz, g_zz_0_xx_xxyyyy, g_zz_0_xx_xxyyyz, g_zz_0_xx_xxyyzz, g_zz_0_xx_xxyzzz, g_zz_0_xx_xxzzzz, g_zz_0_xx_xyyyyy, g_zz_0_xx_xyyyyz, g_zz_0_xx_xyyyzz, g_zz_0_xx_xyyzzz, g_zz_0_xx_xyzzzz, g_zz_0_xx_xzzzzz, g_zz_0_xx_yyyyyy, g_zz_0_xx_yyyyyz, g_zz_0_xx_yyyyzz, g_zz_0_xx_yyyzzz, g_zz_0_xx_yyzzzz, g_zz_0_xx_yzzzzz, g_zz_0_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xx_xxxxxx[k] = -g_zz_0_x_xxxxxx[k] * ab_x + g_zz_0_x_xxxxxxx[k];

                g_zz_0_xx_xxxxxy[k] = -g_zz_0_x_xxxxxy[k] * ab_x + g_zz_0_x_xxxxxxy[k];

                g_zz_0_xx_xxxxxz[k] = -g_zz_0_x_xxxxxz[k] * ab_x + g_zz_0_x_xxxxxxz[k];

                g_zz_0_xx_xxxxyy[k] = -g_zz_0_x_xxxxyy[k] * ab_x + g_zz_0_x_xxxxxyy[k];

                g_zz_0_xx_xxxxyz[k] = -g_zz_0_x_xxxxyz[k] * ab_x + g_zz_0_x_xxxxxyz[k];

                g_zz_0_xx_xxxxzz[k] = -g_zz_0_x_xxxxzz[k] * ab_x + g_zz_0_x_xxxxxzz[k];

                g_zz_0_xx_xxxyyy[k] = -g_zz_0_x_xxxyyy[k] * ab_x + g_zz_0_x_xxxxyyy[k];

                g_zz_0_xx_xxxyyz[k] = -g_zz_0_x_xxxyyz[k] * ab_x + g_zz_0_x_xxxxyyz[k];

                g_zz_0_xx_xxxyzz[k] = -g_zz_0_x_xxxyzz[k] * ab_x + g_zz_0_x_xxxxyzz[k];

                g_zz_0_xx_xxxzzz[k] = -g_zz_0_x_xxxzzz[k] * ab_x + g_zz_0_x_xxxxzzz[k];

                g_zz_0_xx_xxyyyy[k] = -g_zz_0_x_xxyyyy[k] * ab_x + g_zz_0_x_xxxyyyy[k];

                g_zz_0_xx_xxyyyz[k] = -g_zz_0_x_xxyyyz[k] * ab_x + g_zz_0_x_xxxyyyz[k];

                g_zz_0_xx_xxyyzz[k] = -g_zz_0_x_xxyyzz[k] * ab_x + g_zz_0_x_xxxyyzz[k];

                g_zz_0_xx_xxyzzz[k] = -g_zz_0_x_xxyzzz[k] * ab_x + g_zz_0_x_xxxyzzz[k];

                g_zz_0_xx_xxzzzz[k] = -g_zz_0_x_xxzzzz[k] * ab_x + g_zz_0_x_xxxzzzz[k];

                g_zz_0_xx_xyyyyy[k] = -g_zz_0_x_xyyyyy[k] * ab_x + g_zz_0_x_xxyyyyy[k];

                g_zz_0_xx_xyyyyz[k] = -g_zz_0_x_xyyyyz[k] * ab_x + g_zz_0_x_xxyyyyz[k];

                g_zz_0_xx_xyyyzz[k] = -g_zz_0_x_xyyyzz[k] * ab_x + g_zz_0_x_xxyyyzz[k];

                g_zz_0_xx_xyyzzz[k] = -g_zz_0_x_xyyzzz[k] * ab_x + g_zz_0_x_xxyyzzz[k];

                g_zz_0_xx_xyzzzz[k] = -g_zz_0_x_xyzzzz[k] * ab_x + g_zz_0_x_xxyzzzz[k];

                g_zz_0_xx_xzzzzz[k] = -g_zz_0_x_xzzzzz[k] * ab_x + g_zz_0_x_xxzzzzz[k];

                g_zz_0_xx_yyyyyy[k] = -g_zz_0_x_yyyyyy[k] * ab_x + g_zz_0_x_xyyyyyy[k];

                g_zz_0_xx_yyyyyz[k] = -g_zz_0_x_yyyyyz[k] * ab_x + g_zz_0_x_xyyyyyz[k];

                g_zz_0_xx_yyyyzz[k] = -g_zz_0_x_yyyyzz[k] * ab_x + g_zz_0_x_xyyyyzz[k];

                g_zz_0_xx_yyyzzz[k] = -g_zz_0_x_yyyzzz[k] * ab_x + g_zz_0_x_xyyyzzz[k];

                g_zz_0_xx_yyzzzz[k] = -g_zz_0_x_yyzzzz[k] * ab_x + g_zz_0_x_xyyzzzz[k];

                g_zz_0_xx_yzzzzz[k] = -g_zz_0_x_yzzzzz[k] * ab_x + g_zz_0_x_xyzzzzz[k];

                g_zz_0_xx_zzzzzz[k] = -g_zz_0_x_zzzzzz[k] * ab_x + g_zz_0_x_xzzzzzz[k];
            }

            /// Set up 868-896 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xy_xxxxxx = cbuffer.data(di_geom_20_off + 868 * ccomps * dcomps);

            auto g_zz_0_xy_xxxxxy = cbuffer.data(di_geom_20_off + 869 * ccomps * dcomps);

            auto g_zz_0_xy_xxxxxz = cbuffer.data(di_geom_20_off + 870 * ccomps * dcomps);

            auto g_zz_0_xy_xxxxyy = cbuffer.data(di_geom_20_off + 871 * ccomps * dcomps);

            auto g_zz_0_xy_xxxxyz = cbuffer.data(di_geom_20_off + 872 * ccomps * dcomps);

            auto g_zz_0_xy_xxxxzz = cbuffer.data(di_geom_20_off + 873 * ccomps * dcomps);

            auto g_zz_0_xy_xxxyyy = cbuffer.data(di_geom_20_off + 874 * ccomps * dcomps);

            auto g_zz_0_xy_xxxyyz = cbuffer.data(di_geom_20_off + 875 * ccomps * dcomps);

            auto g_zz_0_xy_xxxyzz = cbuffer.data(di_geom_20_off + 876 * ccomps * dcomps);

            auto g_zz_0_xy_xxxzzz = cbuffer.data(di_geom_20_off + 877 * ccomps * dcomps);

            auto g_zz_0_xy_xxyyyy = cbuffer.data(di_geom_20_off + 878 * ccomps * dcomps);

            auto g_zz_0_xy_xxyyyz = cbuffer.data(di_geom_20_off + 879 * ccomps * dcomps);

            auto g_zz_0_xy_xxyyzz = cbuffer.data(di_geom_20_off + 880 * ccomps * dcomps);

            auto g_zz_0_xy_xxyzzz = cbuffer.data(di_geom_20_off + 881 * ccomps * dcomps);

            auto g_zz_0_xy_xxzzzz = cbuffer.data(di_geom_20_off + 882 * ccomps * dcomps);

            auto g_zz_0_xy_xyyyyy = cbuffer.data(di_geom_20_off + 883 * ccomps * dcomps);

            auto g_zz_0_xy_xyyyyz = cbuffer.data(di_geom_20_off + 884 * ccomps * dcomps);

            auto g_zz_0_xy_xyyyzz = cbuffer.data(di_geom_20_off + 885 * ccomps * dcomps);

            auto g_zz_0_xy_xyyzzz = cbuffer.data(di_geom_20_off + 886 * ccomps * dcomps);

            auto g_zz_0_xy_xyzzzz = cbuffer.data(di_geom_20_off + 887 * ccomps * dcomps);

            auto g_zz_0_xy_xzzzzz = cbuffer.data(di_geom_20_off + 888 * ccomps * dcomps);

            auto g_zz_0_xy_yyyyyy = cbuffer.data(di_geom_20_off + 889 * ccomps * dcomps);

            auto g_zz_0_xy_yyyyyz = cbuffer.data(di_geom_20_off + 890 * ccomps * dcomps);

            auto g_zz_0_xy_yyyyzz = cbuffer.data(di_geom_20_off + 891 * ccomps * dcomps);

            auto g_zz_0_xy_yyyzzz = cbuffer.data(di_geom_20_off + 892 * ccomps * dcomps);

            auto g_zz_0_xy_yyzzzz = cbuffer.data(di_geom_20_off + 893 * ccomps * dcomps);

            auto g_zz_0_xy_yzzzzz = cbuffer.data(di_geom_20_off + 894 * ccomps * dcomps);

            auto g_zz_0_xy_zzzzzz = cbuffer.data(di_geom_20_off + 895 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xy_xxxxxx, g_zz_0_xy_xxxxxy, g_zz_0_xy_xxxxxz, g_zz_0_xy_xxxxyy, g_zz_0_xy_xxxxyz, g_zz_0_xy_xxxxzz, g_zz_0_xy_xxxyyy, g_zz_0_xy_xxxyyz, g_zz_0_xy_xxxyzz, g_zz_0_xy_xxxzzz, g_zz_0_xy_xxyyyy, g_zz_0_xy_xxyyyz, g_zz_0_xy_xxyyzz, g_zz_0_xy_xxyzzz, g_zz_0_xy_xxzzzz, g_zz_0_xy_xyyyyy, g_zz_0_xy_xyyyyz, g_zz_0_xy_xyyyzz, g_zz_0_xy_xyyzzz, g_zz_0_xy_xyzzzz, g_zz_0_xy_xzzzzz, g_zz_0_xy_yyyyyy, g_zz_0_xy_yyyyyz, g_zz_0_xy_yyyyzz, g_zz_0_xy_yyyzzz, g_zz_0_xy_yyzzzz, g_zz_0_xy_yzzzzz, g_zz_0_xy_zzzzzz, g_zz_0_y_xxxxxx, g_zz_0_y_xxxxxxx, g_zz_0_y_xxxxxxy, g_zz_0_y_xxxxxxz, g_zz_0_y_xxxxxy, g_zz_0_y_xxxxxyy, g_zz_0_y_xxxxxyz, g_zz_0_y_xxxxxz, g_zz_0_y_xxxxxzz, g_zz_0_y_xxxxyy, g_zz_0_y_xxxxyyy, g_zz_0_y_xxxxyyz, g_zz_0_y_xxxxyz, g_zz_0_y_xxxxyzz, g_zz_0_y_xxxxzz, g_zz_0_y_xxxxzzz, g_zz_0_y_xxxyyy, g_zz_0_y_xxxyyyy, g_zz_0_y_xxxyyyz, g_zz_0_y_xxxyyz, g_zz_0_y_xxxyyzz, g_zz_0_y_xxxyzz, g_zz_0_y_xxxyzzz, g_zz_0_y_xxxzzz, g_zz_0_y_xxxzzzz, g_zz_0_y_xxyyyy, g_zz_0_y_xxyyyyy, g_zz_0_y_xxyyyyz, g_zz_0_y_xxyyyz, g_zz_0_y_xxyyyzz, g_zz_0_y_xxyyzz, g_zz_0_y_xxyyzzz, g_zz_0_y_xxyzzz, g_zz_0_y_xxyzzzz, g_zz_0_y_xxzzzz, g_zz_0_y_xxzzzzz, g_zz_0_y_xyyyyy, g_zz_0_y_xyyyyyy, g_zz_0_y_xyyyyyz, g_zz_0_y_xyyyyz, g_zz_0_y_xyyyyzz, g_zz_0_y_xyyyzz, g_zz_0_y_xyyyzzz, g_zz_0_y_xyyzzz, g_zz_0_y_xyyzzzz, g_zz_0_y_xyzzzz, g_zz_0_y_xyzzzzz, g_zz_0_y_xzzzzz, g_zz_0_y_xzzzzzz, g_zz_0_y_yyyyyy, g_zz_0_y_yyyyyz, g_zz_0_y_yyyyzz, g_zz_0_y_yyyzzz, g_zz_0_y_yyzzzz, g_zz_0_y_yzzzzz, g_zz_0_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xy_xxxxxx[k] = -g_zz_0_y_xxxxxx[k] * ab_x + g_zz_0_y_xxxxxxx[k];

                g_zz_0_xy_xxxxxy[k] = -g_zz_0_y_xxxxxy[k] * ab_x + g_zz_0_y_xxxxxxy[k];

                g_zz_0_xy_xxxxxz[k] = -g_zz_0_y_xxxxxz[k] * ab_x + g_zz_0_y_xxxxxxz[k];

                g_zz_0_xy_xxxxyy[k] = -g_zz_0_y_xxxxyy[k] * ab_x + g_zz_0_y_xxxxxyy[k];

                g_zz_0_xy_xxxxyz[k] = -g_zz_0_y_xxxxyz[k] * ab_x + g_zz_0_y_xxxxxyz[k];

                g_zz_0_xy_xxxxzz[k] = -g_zz_0_y_xxxxzz[k] * ab_x + g_zz_0_y_xxxxxzz[k];

                g_zz_0_xy_xxxyyy[k] = -g_zz_0_y_xxxyyy[k] * ab_x + g_zz_0_y_xxxxyyy[k];

                g_zz_0_xy_xxxyyz[k] = -g_zz_0_y_xxxyyz[k] * ab_x + g_zz_0_y_xxxxyyz[k];

                g_zz_0_xy_xxxyzz[k] = -g_zz_0_y_xxxyzz[k] * ab_x + g_zz_0_y_xxxxyzz[k];

                g_zz_0_xy_xxxzzz[k] = -g_zz_0_y_xxxzzz[k] * ab_x + g_zz_0_y_xxxxzzz[k];

                g_zz_0_xy_xxyyyy[k] = -g_zz_0_y_xxyyyy[k] * ab_x + g_zz_0_y_xxxyyyy[k];

                g_zz_0_xy_xxyyyz[k] = -g_zz_0_y_xxyyyz[k] * ab_x + g_zz_0_y_xxxyyyz[k];

                g_zz_0_xy_xxyyzz[k] = -g_zz_0_y_xxyyzz[k] * ab_x + g_zz_0_y_xxxyyzz[k];

                g_zz_0_xy_xxyzzz[k] = -g_zz_0_y_xxyzzz[k] * ab_x + g_zz_0_y_xxxyzzz[k];

                g_zz_0_xy_xxzzzz[k] = -g_zz_0_y_xxzzzz[k] * ab_x + g_zz_0_y_xxxzzzz[k];

                g_zz_0_xy_xyyyyy[k] = -g_zz_0_y_xyyyyy[k] * ab_x + g_zz_0_y_xxyyyyy[k];

                g_zz_0_xy_xyyyyz[k] = -g_zz_0_y_xyyyyz[k] * ab_x + g_zz_0_y_xxyyyyz[k];

                g_zz_0_xy_xyyyzz[k] = -g_zz_0_y_xyyyzz[k] * ab_x + g_zz_0_y_xxyyyzz[k];

                g_zz_0_xy_xyyzzz[k] = -g_zz_0_y_xyyzzz[k] * ab_x + g_zz_0_y_xxyyzzz[k];

                g_zz_0_xy_xyzzzz[k] = -g_zz_0_y_xyzzzz[k] * ab_x + g_zz_0_y_xxyzzzz[k];

                g_zz_0_xy_xzzzzz[k] = -g_zz_0_y_xzzzzz[k] * ab_x + g_zz_0_y_xxzzzzz[k];

                g_zz_0_xy_yyyyyy[k] = -g_zz_0_y_yyyyyy[k] * ab_x + g_zz_0_y_xyyyyyy[k];

                g_zz_0_xy_yyyyyz[k] = -g_zz_0_y_yyyyyz[k] * ab_x + g_zz_0_y_xyyyyyz[k];

                g_zz_0_xy_yyyyzz[k] = -g_zz_0_y_yyyyzz[k] * ab_x + g_zz_0_y_xyyyyzz[k];

                g_zz_0_xy_yyyzzz[k] = -g_zz_0_y_yyyzzz[k] * ab_x + g_zz_0_y_xyyyzzz[k];

                g_zz_0_xy_yyzzzz[k] = -g_zz_0_y_yyzzzz[k] * ab_x + g_zz_0_y_xyyzzzz[k];

                g_zz_0_xy_yzzzzz[k] = -g_zz_0_y_yzzzzz[k] * ab_x + g_zz_0_y_xyzzzzz[k];

                g_zz_0_xy_zzzzzz[k] = -g_zz_0_y_zzzzzz[k] * ab_x + g_zz_0_y_xzzzzzz[k];
            }

            /// Set up 896-924 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xz_xxxxxx = cbuffer.data(di_geom_20_off + 896 * ccomps * dcomps);

            auto g_zz_0_xz_xxxxxy = cbuffer.data(di_geom_20_off + 897 * ccomps * dcomps);

            auto g_zz_0_xz_xxxxxz = cbuffer.data(di_geom_20_off + 898 * ccomps * dcomps);

            auto g_zz_0_xz_xxxxyy = cbuffer.data(di_geom_20_off + 899 * ccomps * dcomps);

            auto g_zz_0_xz_xxxxyz = cbuffer.data(di_geom_20_off + 900 * ccomps * dcomps);

            auto g_zz_0_xz_xxxxzz = cbuffer.data(di_geom_20_off + 901 * ccomps * dcomps);

            auto g_zz_0_xz_xxxyyy = cbuffer.data(di_geom_20_off + 902 * ccomps * dcomps);

            auto g_zz_0_xz_xxxyyz = cbuffer.data(di_geom_20_off + 903 * ccomps * dcomps);

            auto g_zz_0_xz_xxxyzz = cbuffer.data(di_geom_20_off + 904 * ccomps * dcomps);

            auto g_zz_0_xz_xxxzzz = cbuffer.data(di_geom_20_off + 905 * ccomps * dcomps);

            auto g_zz_0_xz_xxyyyy = cbuffer.data(di_geom_20_off + 906 * ccomps * dcomps);

            auto g_zz_0_xz_xxyyyz = cbuffer.data(di_geom_20_off + 907 * ccomps * dcomps);

            auto g_zz_0_xz_xxyyzz = cbuffer.data(di_geom_20_off + 908 * ccomps * dcomps);

            auto g_zz_0_xz_xxyzzz = cbuffer.data(di_geom_20_off + 909 * ccomps * dcomps);

            auto g_zz_0_xz_xxzzzz = cbuffer.data(di_geom_20_off + 910 * ccomps * dcomps);

            auto g_zz_0_xz_xyyyyy = cbuffer.data(di_geom_20_off + 911 * ccomps * dcomps);

            auto g_zz_0_xz_xyyyyz = cbuffer.data(di_geom_20_off + 912 * ccomps * dcomps);

            auto g_zz_0_xz_xyyyzz = cbuffer.data(di_geom_20_off + 913 * ccomps * dcomps);

            auto g_zz_0_xz_xyyzzz = cbuffer.data(di_geom_20_off + 914 * ccomps * dcomps);

            auto g_zz_0_xz_xyzzzz = cbuffer.data(di_geom_20_off + 915 * ccomps * dcomps);

            auto g_zz_0_xz_xzzzzz = cbuffer.data(di_geom_20_off + 916 * ccomps * dcomps);

            auto g_zz_0_xz_yyyyyy = cbuffer.data(di_geom_20_off + 917 * ccomps * dcomps);

            auto g_zz_0_xz_yyyyyz = cbuffer.data(di_geom_20_off + 918 * ccomps * dcomps);

            auto g_zz_0_xz_yyyyzz = cbuffer.data(di_geom_20_off + 919 * ccomps * dcomps);

            auto g_zz_0_xz_yyyzzz = cbuffer.data(di_geom_20_off + 920 * ccomps * dcomps);

            auto g_zz_0_xz_yyzzzz = cbuffer.data(di_geom_20_off + 921 * ccomps * dcomps);

            auto g_zz_0_xz_yzzzzz = cbuffer.data(di_geom_20_off + 922 * ccomps * dcomps);

            auto g_zz_0_xz_zzzzzz = cbuffer.data(di_geom_20_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xz_xxxxxx, g_zz_0_xz_xxxxxy, g_zz_0_xz_xxxxxz, g_zz_0_xz_xxxxyy, g_zz_0_xz_xxxxyz, g_zz_0_xz_xxxxzz, g_zz_0_xz_xxxyyy, g_zz_0_xz_xxxyyz, g_zz_0_xz_xxxyzz, g_zz_0_xz_xxxzzz, g_zz_0_xz_xxyyyy, g_zz_0_xz_xxyyyz, g_zz_0_xz_xxyyzz, g_zz_0_xz_xxyzzz, g_zz_0_xz_xxzzzz, g_zz_0_xz_xyyyyy, g_zz_0_xz_xyyyyz, g_zz_0_xz_xyyyzz, g_zz_0_xz_xyyzzz, g_zz_0_xz_xyzzzz, g_zz_0_xz_xzzzzz, g_zz_0_xz_yyyyyy, g_zz_0_xz_yyyyyz, g_zz_0_xz_yyyyzz, g_zz_0_xz_yyyzzz, g_zz_0_xz_yyzzzz, g_zz_0_xz_yzzzzz, g_zz_0_xz_zzzzzz, g_zz_0_z_xxxxxx, g_zz_0_z_xxxxxxx, g_zz_0_z_xxxxxxy, g_zz_0_z_xxxxxxz, g_zz_0_z_xxxxxy, g_zz_0_z_xxxxxyy, g_zz_0_z_xxxxxyz, g_zz_0_z_xxxxxz, g_zz_0_z_xxxxxzz, g_zz_0_z_xxxxyy, g_zz_0_z_xxxxyyy, g_zz_0_z_xxxxyyz, g_zz_0_z_xxxxyz, g_zz_0_z_xxxxyzz, g_zz_0_z_xxxxzz, g_zz_0_z_xxxxzzz, g_zz_0_z_xxxyyy, g_zz_0_z_xxxyyyy, g_zz_0_z_xxxyyyz, g_zz_0_z_xxxyyz, g_zz_0_z_xxxyyzz, g_zz_0_z_xxxyzz, g_zz_0_z_xxxyzzz, g_zz_0_z_xxxzzz, g_zz_0_z_xxxzzzz, g_zz_0_z_xxyyyy, g_zz_0_z_xxyyyyy, g_zz_0_z_xxyyyyz, g_zz_0_z_xxyyyz, g_zz_0_z_xxyyyzz, g_zz_0_z_xxyyzz, g_zz_0_z_xxyyzzz, g_zz_0_z_xxyzzz, g_zz_0_z_xxyzzzz, g_zz_0_z_xxzzzz, g_zz_0_z_xxzzzzz, g_zz_0_z_xyyyyy, g_zz_0_z_xyyyyyy, g_zz_0_z_xyyyyyz, g_zz_0_z_xyyyyz, g_zz_0_z_xyyyyzz, g_zz_0_z_xyyyzz, g_zz_0_z_xyyyzzz, g_zz_0_z_xyyzzz, g_zz_0_z_xyyzzzz, g_zz_0_z_xyzzzz, g_zz_0_z_xyzzzzz, g_zz_0_z_xzzzzz, g_zz_0_z_xzzzzzz, g_zz_0_z_yyyyyy, g_zz_0_z_yyyyyz, g_zz_0_z_yyyyzz, g_zz_0_z_yyyzzz, g_zz_0_z_yyzzzz, g_zz_0_z_yzzzzz, g_zz_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xz_xxxxxx[k] = -g_zz_0_z_xxxxxx[k] * ab_x + g_zz_0_z_xxxxxxx[k];

                g_zz_0_xz_xxxxxy[k] = -g_zz_0_z_xxxxxy[k] * ab_x + g_zz_0_z_xxxxxxy[k];

                g_zz_0_xz_xxxxxz[k] = -g_zz_0_z_xxxxxz[k] * ab_x + g_zz_0_z_xxxxxxz[k];

                g_zz_0_xz_xxxxyy[k] = -g_zz_0_z_xxxxyy[k] * ab_x + g_zz_0_z_xxxxxyy[k];

                g_zz_0_xz_xxxxyz[k] = -g_zz_0_z_xxxxyz[k] * ab_x + g_zz_0_z_xxxxxyz[k];

                g_zz_0_xz_xxxxzz[k] = -g_zz_0_z_xxxxzz[k] * ab_x + g_zz_0_z_xxxxxzz[k];

                g_zz_0_xz_xxxyyy[k] = -g_zz_0_z_xxxyyy[k] * ab_x + g_zz_0_z_xxxxyyy[k];

                g_zz_0_xz_xxxyyz[k] = -g_zz_0_z_xxxyyz[k] * ab_x + g_zz_0_z_xxxxyyz[k];

                g_zz_0_xz_xxxyzz[k] = -g_zz_0_z_xxxyzz[k] * ab_x + g_zz_0_z_xxxxyzz[k];

                g_zz_0_xz_xxxzzz[k] = -g_zz_0_z_xxxzzz[k] * ab_x + g_zz_0_z_xxxxzzz[k];

                g_zz_0_xz_xxyyyy[k] = -g_zz_0_z_xxyyyy[k] * ab_x + g_zz_0_z_xxxyyyy[k];

                g_zz_0_xz_xxyyyz[k] = -g_zz_0_z_xxyyyz[k] * ab_x + g_zz_0_z_xxxyyyz[k];

                g_zz_0_xz_xxyyzz[k] = -g_zz_0_z_xxyyzz[k] * ab_x + g_zz_0_z_xxxyyzz[k];

                g_zz_0_xz_xxyzzz[k] = -g_zz_0_z_xxyzzz[k] * ab_x + g_zz_0_z_xxxyzzz[k];

                g_zz_0_xz_xxzzzz[k] = -g_zz_0_z_xxzzzz[k] * ab_x + g_zz_0_z_xxxzzzz[k];

                g_zz_0_xz_xyyyyy[k] = -g_zz_0_z_xyyyyy[k] * ab_x + g_zz_0_z_xxyyyyy[k];

                g_zz_0_xz_xyyyyz[k] = -g_zz_0_z_xyyyyz[k] * ab_x + g_zz_0_z_xxyyyyz[k];

                g_zz_0_xz_xyyyzz[k] = -g_zz_0_z_xyyyzz[k] * ab_x + g_zz_0_z_xxyyyzz[k];

                g_zz_0_xz_xyyzzz[k] = -g_zz_0_z_xyyzzz[k] * ab_x + g_zz_0_z_xxyyzzz[k];

                g_zz_0_xz_xyzzzz[k] = -g_zz_0_z_xyzzzz[k] * ab_x + g_zz_0_z_xxyzzzz[k];

                g_zz_0_xz_xzzzzz[k] = -g_zz_0_z_xzzzzz[k] * ab_x + g_zz_0_z_xxzzzzz[k];

                g_zz_0_xz_yyyyyy[k] = -g_zz_0_z_yyyyyy[k] * ab_x + g_zz_0_z_xyyyyyy[k];

                g_zz_0_xz_yyyyyz[k] = -g_zz_0_z_yyyyyz[k] * ab_x + g_zz_0_z_xyyyyyz[k];

                g_zz_0_xz_yyyyzz[k] = -g_zz_0_z_yyyyzz[k] * ab_x + g_zz_0_z_xyyyyzz[k];

                g_zz_0_xz_yyyzzz[k] = -g_zz_0_z_yyyzzz[k] * ab_x + g_zz_0_z_xyyyzzz[k];

                g_zz_0_xz_yyzzzz[k] = -g_zz_0_z_yyzzzz[k] * ab_x + g_zz_0_z_xyyzzzz[k];

                g_zz_0_xz_yzzzzz[k] = -g_zz_0_z_yzzzzz[k] * ab_x + g_zz_0_z_xyzzzzz[k];

                g_zz_0_xz_zzzzzz[k] = -g_zz_0_z_zzzzzz[k] * ab_x + g_zz_0_z_xzzzzzz[k];
            }

            /// Set up 924-952 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yy_xxxxxx = cbuffer.data(di_geom_20_off + 924 * ccomps * dcomps);

            auto g_zz_0_yy_xxxxxy = cbuffer.data(di_geom_20_off + 925 * ccomps * dcomps);

            auto g_zz_0_yy_xxxxxz = cbuffer.data(di_geom_20_off + 926 * ccomps * dcomps);

            auto g_zz_0_yy_xxxxyy = cbuffer.data(di_geom_20_off + 927 * ccomps * dcomps);

            auto g_zz_0_yy_xxxxyz = cbuffer.data(di_geom_20_off + 928 * ccomps * dcomps);

            auto g_zz_0_yy_xxxxzz = cbuffer.data(di_geom_20_off + 929 * ccomps * dcomps);

            auto g_zz_0_yy_xxxyyy = cbuffer.data(di_geom_20_off + 930 * ccomps * dcomps);

            auto g_zz_0_yy_xxxyyz = cbuffer.data(di_geom_20_off + 931 * ccomps * dcomps);

            auto g_zz_0_yy_xxxyzz = cbuffer.data(di_geom_20_off + 932 * ccomps * dcomps);

            auto g_zz_0_yy_xxxzzz = cbuffer.data(di_geom_20_off + 933 * ccomps * dcomps);

            auto g_zz_0_yy_xxyyyy = cbuffer.data(di_geom_20_off + 934 * ccomps * dcomps);

            auto g_zz_0_yy_xxyyyz = cbuffer.data(di_geom_20_off + 935 * ccomps * dcomps);

            auto g_zz_0_yy_xxyyzz = cbuffer.data(di_geom_20_off + 936 * ccomps * dcomps);

            auto g_zz_0_yy_xxyzzz = cbuffer.data(di_geom_20_off + 937 * ccomps * dcomps);

            auto g_zz_0_yy_xxzzzz = cbuffer.data(di_geom_20_off + 938 * ccomps * dcomps);

            auto g_zz_0_yy_xyyyyy = cbuffer.data(di_geom_20_off + 939 * ccomps * dcomps);

            auto g_zz_0_yy_xyyyyz = cbuffer.data(di_geom_20_off + 940 * ccomps * dcomps);

            auto g_zz_0_yy_xyyyzz = cbuffer.data(di_geom_20_off + 941 * ccomps * dcomps);

            auto g_zz_0_yy_xyyzzz = cbuffer.data(di_geom_20_off + 942 * ccomps * dcomps);

            auto g_zz_0_yy_xyzzzz = cbuffer.data(di_geom_20_off + 943 * ccomps * dcomps);

            auto g_zz_0_yy_xzzzzz = cbuffer.data(di_geom_20_off + 944 * ccomps * dcomps);

            auto g_zz_0_yy_yyyyyy = cbuffer.data(di_geom_20_off + 945 * ccomps * dcomps);

            auto g_zz_0_yy_yyyyyz = cbuffer.data(di_geom_20_off + 946 * ccomps * dcomps);

            auto g_zz_0_yy_yyyyzz = cbuffer.data(di_geom_20_off + 947 * ccomps * dcomps);

            auto g_zz_0_yy_yyyzzz = cbuffer.data(di_geom_20_off + 948 * ccomps * dcomps);

            auto g_zz_0_yy_yyzzzz = cbuffer.data(di_geom_20_off + 949 * ccomps * dcomps);

            auto g_zz_0_yy_yzzzzz = cbuffer.data(di_geom_20_off + 950 * ccomps * dcomps);

            auto g_zz_0_yy_zzzzzz = cbuffer.data(di_geom_20_off + 951 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_y_xxxxxx, g_zz_0_y_xxxxxxy, g_zz_0_y_xxxxxy, g_zz_0_y_xxxxxyy, g_zz_0_y_xxxxxyz, g_zz_0_y_xxxxxz, g_zz_0_y_xxxxyy, g_zz_0_y_xxxxyyy, g_zz_0_y_xxxxyyz, g_zz_0_y_xxxxyz, g_zz_0_y_xxxxyzz, g_zz_0_y_xxxxzz, g_zz_0_y_xxxyyy, g_zz_0_y_xxxyyyy, g_zz_0_y_xxxyyyz, g_zz_0_y_xxxyyz, g_zz_0_y_xxxyyzz, g_zz_0_y_xxxyzz, g_zz_0_y_xxxyzzz, g_zz_0_y_xxxzzz, g_zz_0_y_xxyyyy, g_zz_0_y_xxyyyyy, g_zz_0_y_xxyyyyz, g_zz_0_y_xxyyyz, g_zz_0_y_xxyyyzz, g_zz_0_y_xxyyzz, g_zz_0_y_xxyyzzz, g_zz_0_y_xxyzzz, g_zz_0_y_xxyzzzz, g_zz_0_y_xxzzzz, g_zz_0_y_xyyyyy, g_zz_0_y_xyyyyyy, g_zz_0_y_xyyyyyz, g_zz_0_y_xyyyyz, g_zz_0_y_xyyyyzz, g_zz_0_y_xyyyzz, g_zz_0_y_xyyyzzz, g_zz_0_y_xyyzzz, g_zz_0_y_xyyzzzz, g_zz_0_y_xyzzzz, g_zz_0_y_xyzzzzz, g_zz_0_y_xzzzzz, g_zz_0_y_yyyyyy, g_zz_0_y_yyyyyyy, g_zz_0_y_yyyyyyz, g_zz_0_y_yyyyyz, g_zz_0_y_yyyyyzz, g_zz_0_y_yyyyzz, g_zz_0_y_yyyyzzz, g_zz_0_y_yyyzzz, g_zz_0_y_yyyzzzz, g_zz_0_y_yyzzzz, g_zz_0_y_yyzzzzz, g_zz_0_y_yzzzzz, g_zz_0_y_yzzzzzz, g_zz_0_y_zzzzzz, g_zz_0_yy_xxxxxx, g_zz_0_yy_xxxxxy, g_zz_0_yy_xxxxxz, g_zz_0_yy_xxxxyy, g_zz_0_yy_xxxxyz, g_zz_0_yy_xxxxzz, g_zz_0_yy_xxxyyy, g_zz_0_yy_xxxyyz, g_zz_0_yy_xxxyzz, g_zz_0_yy_xxxzzz, g_zz_0_yy_xxyyyy, g_zz_0_yy_xxyyyz, g_zz_0_yy_xxyyzz, g_zz_0_yy_xxyzzz, g_zz_0_yy_xxzzzz, g_zz_0_yy_xyyyyy, g_zz_0_yy_xyyyyz, g_zz_0_yy_xyyyzz, g_zz_0_yy_xyyzzz, g_zz_0_yy_xyzzzz, g_zz_0_yy_xzzzzz, g_zz_0_yy_yyyyyy, g_zz_0_yy_yyyyyz, g_zz_0_yy_yyyyzz, g_zz_0_yy_yyyzzz, g_zz_0_yy_yyzzzz, g_zz_0_yy_yzzzzz, g_zz_0_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yy_xxxxxx[k] = -g_zz_0_y_xxxxxx[k] * ab_y + g_zz_0_y_xxxxxxy[k];

                g_zz_0_yy_xxxxxy[k] = -g_zz_0_y_xxxxxy[k] * ab_y + g_zz_0_y_xxxxxyy[k];

                g_zz_0_yy_xxxxxz[k] = -g_zz_0_y_xxxxxz[k] * ab_y + g_zz_0_y_xxxxxyz[k];

                g_zz_0_yy_xxxxyy[k] = -g_zz_0_y_xxxxyy[k] * ab_y + g_zz_0_y_xxxxyyy[k];

                g_zz_0_yy_xxxxyz[k] = -g_zz_0_y_xxxxyz[k] * ab_y + g_zz_0_y_xxxxyyz[k];

                g_zz_0_yy_xxxxzz[k] = -g_zz_0_y_xxxxzz[k] * ab_y + g_zz_0_y_xxxxyzz[k];

                g_zz_0_yy_xxxyyy[k] = -g_zz_0_y_xxxyyy[k] * ab_y + g_zz_0_y_xxxyyyy[k];

                g_zz_0_yy_xxxyyz[k] = -g_zz_0_y_xxxyyz[k] * ab_y + g_zz_0_y_xxxyyyz[k];

                g_zz_0_yy_xxxyzz[k] = -g_zz_0_y_xxxyzz[k] * ab_y + g_zz_0_y_xxxyyzz[k];

                g_zz_0_yy_xxxzzz[k] = -g_zz_0_y_xxxzzz[k] * ab_y + g_zz_0_y_xxxyzzz[k];

                g_zz_0_yy_xxyyyy[k] = -g_zz_0_y_xxyyyy[k] * ab_y + g_zz_0_y_xxyyyyy[k];

                g_zz_0_yy_xxyyyz[k] = -g_zz_0_y_xxyyyz[k] * ab_y + g_zz_0_y_xxyyyyz[k];

                g_zz_0_yy_xxyyzz[k] = -g_zz_0_y_xxyyzz[k] * ab_y + g_zz_0_y_xxyyyzz[k];

                g_zz_0_yy_xxyzzz[k] = -g_zz_0_y_xxyzzz[k] * ab_y + g_zz_0_y_xxyyzzz[k];

                g_zz_0_yy_xxzzzz[k] = -g_zz_0_y_xxzzzz[k] * ab_y + g_zz_0_y_xxyzzzz[k];

                g_zz_0_yy_xyyyyy[k] = -g_zz_0_y_xyyyyy[k] * ab_y + g_zz_0_y_xyyyyyy[k];

                g_zz_0_yy_xyyyyz[k] = -g_zz_0_y_xyyyyz[k] * ab_y + g_zz_0_y_xyyyyyz[k];

                g_zz_0_yy_xyyyzz[k] = -g_zz_0_y_xyyyzz[k] * ab_y + g_zz_0_y_xyyyyzz[k];

                g_zz_0_yy_xyyzzz[k] = -g_zz_0_y_xyyzzz[k] * ab_y + g_zz_0_y_xyyyzzz[k];

                g_zz_0_yy_xyzzzz[k] = -g_zz_0_y_xyzzzz[k] * ab_y + g_zz_0_y_xyyzzzz[k];

                g_zz_0_yy_xzzzzz[k] = -g_zz_0_y_xzzzzz[k] * ab_y + g_zz_0_y_xyzzzzz[k];

                g_zz_0_yy_yyyyyy[k] = -g_zz_0_y_yyyyyy[k] * ab_y + g_zz_0_y_yyyyyyy[k];

                g_zz_0_yy_yyyyyz[k] = -g_zz_0_y_yyyyyz[k] * ab_y + g_zz_0_y_yyyyyyz[k];

                g_zz_0_yy_yyyyzz[k] = -g_zz_0_y_yyyyzz[k] * ab_y + g_zz_0_y_yyyyyzz[k];

                g_zz_0_yy_yyyzzz[k] = -g_zz_0_y_yyyzzz[k] * ab_y + g_zz_0_y_yyyyzzz[k];

                g_zz_0_yy_yyzzzz[k] = -g_zz_0_y_yyzzzz[k] * ab_y + g_zz_0_y_yyyzzzz[k];

                g_zz_0_yy_yzzzzz[k] = -g_zz_0_y_yzzzzz[k] * ab_y + g_zz_0_y_yyzzzzz[k];

                g_zz_0_yy_zzzzzz[k] = -g_zz_0_y_zzzzzz[k] * ab_y + g_zz_0_y_yzzzzzz[k];
            }

            /// Set up 952-980 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yz_xxxxxx = cbuffer.data(di_geom_20_off + 952 * ccomps * dcomps);

            auto g_zz_0_yz_xxxxxy = cbuffer.data(di_geom_20_off + 953 * ccomps * dcomps);

            auto g_zz_0_yz_xxxxxz = cbuffer.data(di_geom_20_off + 954 * ccomps * dcomps);

            auto g_zz_0_yz_xxxxyy = cbuffer.data(di_geom_20_off + 955 * ccomps * dcomps);

            auto g_zz_0_yz_xxxxyz = cbuffer.data(di_geom_20_off + 956 * ccomps * dcomps);

            auto g_zz_0_yz_xxxxzz = cbuffer.data(di_geom_20_off + 957 * ccomps * dcomps);

            auto g_zz_0_yz_xxxyyy = cbuffer.data(di_geom_20_off + 958 * ccomps * dcomps);

            auto g_zz_0_yz_xxxyyz = cbuffer.data(di_geom_20_off + 959 * ccomps * dcomps);

            auto g_zz_0_yz_xxxyzz = cbuffer.data(di_geom_20_off + 960 * ccomps * dcomps);

            auto g_zz_0_yz_xxxzzz = cbuffer.data(di_geom_20_off + 961 * ccomps * dcomps);

            auto g_zz_0_yz_xxyyyy = cbuffer.data(di_geom_20_off + 962 * ccomps * dcomps);

            auto g_zz_0_yz_xxyyyz = cbuffer.data(di_geom_20_off + 963 * ccomps * dcomps);

            auto g_zz_0_yz_xxyyzz = cbuffer.data(di_geom_20_off + 964 * ccomps * dcomps);

            auto g_zz_0_yz_xxyzzz = cbuffer.data(di_geom_20_off + 965 * ccomps * dcomps);

            auto g_zz_0_yz_xxzzzz = cbuffer.data(di_geom_20_off + 966 * ccomps * dcomps);

            auto g_zz_0_yz_xyyyyy = cbuffer.data(di_geom_20_off + 967 * ccomps * dcomps);

            auto g_zz_0_yz_xyyyyz = cbuffer.data(di_geom_20_off + 968 * ccomps * dcomps);

            auto g_zz_0_yz_xyyyzz = cbuffer.data(di_geom_20_off + 969 * ccomps * dcomps);

            auto g_zz_0_yz_xyyzzz = cbuffer.data(di_geom_20_off + 970 * ccomps * dcomps);

            auto g_zz_0_yz_xyzzzz = cbuffer.data(di_geom_20_off + 971 * ccomps * dcomps);

            auto g_zz_0_yz_xzzzzz = cbuffer.data(di_geom_20_off + 972 * ccomps * dcomps);

            auto g_zz_0_yz_yyyyyy = cbuffer.data(di_geom_20_off + 973 * ccomps * dcomps);

            auto g_zz_0_yz_yyyyyz = cbuffer.data(di_geom_20_off + 974 * ccomps * dcomps);

            auto g_zz_0_yz_yyyyzz = cbuffer.data(di_geom_20_off + 975 * ccomps * dcomps);

            auto g_zz_0_yz_yyyzzz = cbuffer.data(di_geom_20_off + 976 * ccomps * dcomps);

            auto g_zz_0_yz_yyzzzz = cbuffer.data(di_geom_20_off + 977 * ccomps * dcomps);

            auto g_zz_0_yz_yzzzzz = cbuffer.data(di_geom_20_off + 978 * ccomps * dcomps);

            auto g_zz_0_yz_zzzzzz = cbuffer.data(di_geom_20_off + 979 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yz_xxxxxx, g_zz_0_yz_xxxxxy, g_zz_0_yz_xxxxxz, g_zz_0_yz_xxxxyy, g_zz_0_yz_xxxxyz, g_zz_0_yz_xxxxzz, g_zz_0_yz_xxxyyy, g_zz_0_yz_xxxyyz, g_zz_0_yz_xxxyzz, g_zz_0_yz_xxxzzz, g_zz_0_yz_xxyyyy, g_zz_0_yz_xxyyyz, g_zz_0_yz_xxyyzz, g_zz_0_yz_xxyzzz, g_zz_0_yz_xxzzzz, g_zz_0_yz_xyyyyy, g_zz_0_yz_xyyyyz, g_zz_0_yz_xyyyzz, g_zz_0_yz_xyyzzz, g_zz_0_yz_xyzzzz, g_zz_0_yz_xzzzzz, g_zz_0_yz_yyyyyy, g_zz_0_yz_yyyyyz, g_zz_0_yz_yyyyzz, g_zz_0_yz_yyyzzz, g_zz_0_yz_yyzzzz, g_zz_0_yz_yzzzzz, g_zz_0_yz_zzzzzz, g_zz_0_z_xxxxxx, g_zz_0_z_xxxxxxy, g_zz_0_z_xxxxxy, g_zz_0_z_xxxxxyy, g_zz_0_z_xxxxxyz, g_zz_0_z_xxxxxz, g_zz_0_z_xxxxyy, g_zz_0_z_xxxxyyy, g_zz_0_z_xxxxyyz, g_zz_0_z_xxxxyz, g_zz_0_z_xxxxyzz, g_zz_0_z_xxxxzz, g_zz_0_z_xxxyyy, g_zz_0_z_xxxyyyy, g_zz_0_z_xxxyyyz, g_zz_0_z_xxxyyz, g_zz_0_z_xxxyyzz, g_zz_0_z_xxxyzz, g_zz_0_z_xxxyzzz, g_zz_0_z_xxxzzz, g_zz_0_z_xxyyyy, g_zz_0_z_xxyyyyy, g_zz_0_z_xxyyyyz, g_zz_0_z_xxyyyz, g_zz_0_z_xxyyyzz, g_zz_0_z_xxyyzz, g_zz_0_z_xxyyzzz, g_zz_0_z_xxyzzz, g_zz_0_z_xxyzzzz, g_zz_0_z_xxzzzz, g_zz_0_z_xyyyyy, g_zz_0_z_xyyyyyy, g_zz_0_z_xyyyyyz, g_zz_0_z_xyyyyz, g_zz_0_z_xyyyyzz, g_zz_0_z_xyyyzz, g_zz_0_z_xyyyzzz, g_zz_0_z_xyyzzz, g_zz_0_z_xyyzzzz, g_zz_0_z_xyzzzz, g_zz_0_z_xyzzzzz, g_zz_0_z_xzzzzz, g_zz_0_z_yyyyyy, g_zz_0_z_yyyyyyy, g_zz_0_z_yyyyyyz, g_zz_0_z_yyyyyz, g_zz_0_z_yyyyyzz, g_zz_0_z_yyyyzz, g_zz_0_z_yyyyzzz, g_zz_0_z_yyyzzz, g_zz_0_z_yyyzzzz, g_zz_0_z_yyzzzz, g_zz_0_z_yyzzzzz, g_zz_0_z_yzzzzz, g_zz_0_z_yzzzzzz, g_zz_0_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yz_xxxxxx[k] = -g_zz_0_z_xxxxxx[k] * ab_y + g_zz_0_z_xxxxxxy[k];

                g_zz_0_yz_xxxxxy[k] = -g_zz_0_z_xxxxxy[k] * ab_y + g_zz_0_z_xxxxxyy[k];

                g_zz_0_yz_xxxxxz[k] = -g_zz_0_z_xxxxxz[k] * ab_y + g_zz_0_z_xxxxxyz[k];

                g_zz_0_yz_xxxxyy[k] = -g_zz_0_z_xxxxyy[k] * ab_y + g_zz_0_z_xxxxyyy[k];

                g_zz_0_yz_xxxxyz[k] = -g_zz_0_z_xxxxyz[k] * ab_y + g_zz_0_z_xxxxyyz[k];

                g_zz_0_yz_xxxxzz[k] = -g_zz_0_z_xxxxzz[k] * ab_y + g_zz_0_z_xxxxyzz[k];

                g_zz_0_yz_xxxyyy[k] = -g_zz_0_z_xxxyyy[k] * ab_y + g_zz_0_z_xxxyyyy[k];

                g_zz_0_yz_xxxyyz[k] = -g_zz_0_z_xxxyyz[k] * ab_y + g_zz_0_z_xxxyyyz[k];

                g_zz_0_yz_xxxyzz[k] = -g_zz_0_z_xxxyzz[k] * ab_y + g_zz_0_z_xxxyyzz[k];

                g_zz_0_yz_xxxzzz[k] = -g_zz_0_z_xxxzzz[k] * ab_y + g_zz_0_z_xxxyzzz[k];

                g_zz_0_yz_xxyyyy[k] = -g_zz_0_z_xxyyyy[k] * ab_y + g_zz_0_z_xxyyyyy[k];

                g_zz_0_yz_xxyyyz[k] = -g_zz_0_z_xxyyyz[k] * ab_y + g_zz_0_z_xxyyyyz[k];

                g_zz_0_yz_xxyyzz[k] = -g_zz_0_z_xxyyzz[k] * ab_y + g_zz_0_z_xxyyyzz[k];

                g_zz_0_yz_xxyzzz[k] = -g_zz_0_z_xxyzzz[k] * ab_y + g_zz_0_z_xxyyzzz[k];

                g_zz_0_yz_xxzzzz[k] = -g_zz_0_z_xxzzzz[k] * ab_y + g_zz_0_z_xxyzzzz[k];

                g_zz_0_yz_xyyyyy[k] = -g_zz_0_z_xyyyyy[k] * ab_y + g_zz_0_z_xyyyyyy[k];

                g_zz_0_yz_xyyyyz[k] = -g_zz_0_z_xyyyyz[k] * ab_y + g_zz_0_z_xyyyyyz[k];

                g_zz_0_yz_xyyyzz[k] = -g_zz_0_z_xyyyzz[k] * ab_y + g_zz_0_z_xyyyyzz[k];

                g_zz_0_yz_xyyzzz[k] = -g_zz_0_z_xyyzzz[k] * ab_y + g_zz_0_z_xyyyzzz[k];

                g_zz_0_yz_xyzzzz[k] = -g_zz_0_z_xyzzzz[k] * ab_y + g_zz_0_z_xyyzzzz[k];

                g_zz_0_yz_xzzzzz[k] = -g_zz_0_z_xzzzzz[k] * ab_y + g_zz_0_z_xyzzzzz[k];

                g_zz_0_yz_yyyyyy[k] = -g_zz_0_z_yyyyyy[k] * ab_y + g_zz_0_z_yyyyyyy[k];

                g_zz_0_yz_yyyyyz[k] = -g_zz_0_z_yyyyyz[k] * ab_y + g_zz_0_z_yyyyyyz[k];

                g_zz_0_yz_yyyyzz[k] = -g_zz_0_z_yyyyzz[k] * ab_y + g_zz_0_z_yyyyyzz[k];

                g_zz_0_yz_yyyzzz[k] = -g_zz_0_z_yyyzzz[k] * ab_y + g_zz_0_z_yyyyzzz[k];

                g_zz_0_yz_yyzzzz[k] = -g_zz_0_z_yyzzzz[k] * ab_y + g_zz_0_z_yyyzzzz[k];

                g_zz_0_yz_yzzzzz[k] = -g_zz_0_z_yzzzzz[k] * ab_y + g_zz_0_z_yyzzzzz[k];

                g_zz_0_yz_zzzzzz[k] = -g_zz_0_z_zzzzzz[k] * ab_y + g_zz_0_z_yzzzzzz[k];
            }

            /// Set up 980-1008 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zz_xxxxxx = cbuffer.data(di_geom_20_off + 980 * ccomps * dcomps);

            auto g_zz_0_zz_xxxxxy = cbuffer.data(di_geom_20_off + 981 * ccomps * dcomps);

            auto g_zz_0_zz_xxxxxz = cbuffer.data(di_geom_20_off + 982 * ccomps * dcomps);

            auto g_zz_0_zz_xxxxyy = cbuffer.data(di_geom_20_off + 983 * ccomps * dcomps);

            auto g_zz_0_zz_xxxxyz = cbuffer.data(di_geom_20_off + 984 * ccomps * dcomps);

            auto g_zz_0_zz_xxxxzz = cbuffer.data(di_geom_20_off + 985 * ccomps * dcomps);

            auto g_zz_0_zz_xxxyyy = cbuffer.data(di_geom_20_off + 986 * ccomps * dcomps);

            auto g_zz_0_zz_xxxyyz = cbuffer.data(di_geom_20_off + 987 * ccomps * dcomps);

            auto g_zz_0_zz_xxxyzz = cbuffer.data(di_geom_20_off + 988 * ccomps * dcomps);

            auto g_zz_0_zz_xxxzzz = cbuffer.data(di_geom_20_off + 989 * ccomps * dcomps);

            auto g_zz_0_zz_xxyyyy = cbuffer.data(di_geom_20_off + 990 * ccomps * dcomps);

            auto g_zz_0_zz_xxyyyz = cbuffer.data(di_geom_20_off + 991 * ccomps * dcomps);

            auto g_zz_0_zz_xxyyzz = cbuffer.data(di_geom_20_off + 992 * ccomps * dcomps);

            auto g_zz_0_zz_xxyzzz = cbuffer.data(di_geom_20_off + 993 * ccomps * dcomps);

            auto g_zz_0_zz_xxzzzz = cbuffer.data(di_geom_20_off + 994 * ccomps * dcomps);

            auto g_zz_0_zz_xyyyyy = cbuffer.data(di_geom_20_off + 995 * ccomps * dcomps);

            auto g_zz_0_zz_xyyyyz = cbuffer.data(di_geom_20_off + 996 * ccomps * dcomps);

            auto g_zz_0_zz_xyyyzz = cbuffer.data(di_geom_20_off + 997 * ccomps * dcomps);

            auto g_zz_0_zz_xyyzzz = cbuffer.data(di_geom_20_off + 998 * ccomps * dcomps);

            auto g_zz_0_zz_xyzzzz = cbuffer.data(di_geom_20_off + 999 * ccomps * dcomps);

            auto g_zz_0_zz_xzzzzz = cbuffer.data(di_geom_20_off + 1000 * ccomps * dcomps);

            auto g_zz_0_zz_yyyyyy = cbuffer.data(di_geom_20_off + 1001 * ccomps * dcomps);

            auto g_zz_0_zz_yyyyyz = cbuffer.data(di_geom_20_off + 1002 * ccomps * dcomps);

            auto g_zz_0_zz_yyyyzz = cbuffer.data(di_geom_20_off + 1003 * ccomps * dcomps);

            auto g_zz_0_zz_yyyzzz = cbuffer.data(di_geom_20_off + 1004 * ccomps * dcomps);

            auto g_zz_0_zz_yyzzzz = cbuffer.data(di_geom_20_off + 1005 * ccomps * dcomps);

            auto g_zz_0_zz_yzzzzz = cbuffer.data(di_geom_20_off + 1006 * ccomps * dcomps);

            auto g_zz_0_zz_zzzzzz = cbuffer.data(di_geom_20_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxxxx, g_z_0_z_xxxxxy, g_z_0_z_xxxxxz, g_z_0_z_xxxxyy, g_z_0_z_xxxxyz, g_z_0_z_xxxxzz, g_z_0_z_xxxyyy, g_z_0_z_xxxyyz, g_z_0_z_xxxyzz, g_z_0_z_xxxzzz, g_z_0_z_xxyyyy, g_z_0_z_xxyyyz, g_z_0_z_xxyyzz, g_z_0_z_xxyzzz, g_z_0_z_xxzzzz, g_z_0_z_xyyyyy, g_z_0_z_xyyyyz, g_z_0_z_xyyyzz, g_z_0_z_xyyzzz, g_z_0_z_xyzzzz, g_z_0_z_xzzzzz, g_z_0_z_yyyyyy, g_z_0_z_yyyyyz, g_z_0_z_yyyyzz, g_z_0_z_yyyzzz, g_z_0_z_yyzzzz, g_z_0_z_yzzzzz, g_z_0_z_zzzzzz, g_zz_0_z_xxxxxx, g_zz_0_z_xxxxxxz, g_zz_0_z_xxxxxy, g_zz_0_z_xxxxxyz, g_zz_0_z_xxxxxz, g_zz_0_z_xxxxxzz, g_zz_0_z_xxxxyy, g_zz_0_z_xxxxyyz, g_zz_0_z_xxxxyz, g_zz_0_z_xxxxyzz, g_zz_0_z_xxxxzz, g_zz_0_z_xxxxzzz, g_zz_0_z_xxxyyy, g_zz_0_z_xxxyyyz, g_zz_0_z_xxxyyz, g_zz_0_z_xxxyyzz, g_zz_0_z_xxxyzz, g_zz_0_z_xxxyzzz, g_zz_0_z_xxxzzz, g_zz_0_z_xxxzzzz, g_zz_0_z_xxyyyy, g_zz_0_z_xxyyyyz, g_zz_0_z_xxyyyz, g_zz_0_z_xxyyyzz, g_zz_0_z_xxyyzz, g_zz_0_z_xxyyzzz, g_zz_0_z_xxyzzz, g_zz_0_z_xxyzzzz, g_zz_0_z_xxzzzz, g_zz_0_z_xxzzzzz, g_zz_0_z_xyyyyy, g_zz_0_z_xyyyyyz, g_zz_0_z_xyyyyz, g_zz_0_z_xyyyyzz, g_zz_0_z_xyyyzz, g_zz_0_z_xyyyzzz, g_zz_0_z_xyyzzz, g_zz_0_z_xyyzzzz, g_zz_0_z_xyzzzz, g_zz_0_z_xyzzzzz, g_zz_0_z_xzzzzz, g_zz_0_z_xzzzzzz, g_zz_0_z_yyyyyy, g_zz_0_z_yyyyyyz, g_zz_0_z_yyyyyz, g_zz_0_z_yyyyyzz, g_zz_0_z_yyyyzz, g_zz_0_z_yyyyzzz, g_zz_0_z_yyyzzz, g_zz_0_z_yyyzzzz, g_zz_0_z_yyzzzz, g_zz_0_z_yyzzzzz, g_zz_0_z_yzzzzz, g_zz_0_z_yzzzzzz, g_zz_0_z_zzzzzz, g_zz_0_z_zzzzzzz, g_zz_0_zz_xxxxxx, g_zz_0_zz_xxxxxy, g_zz_0_zz_xxxxxz, g_zz_0_zz_xxxxyy, g_zz_0_zz_xxxxyz, g_zz_0_zz_xxxxzz, g_zz_0_zz_xxxyyy, g_zz_0_zz_xxxyyz, g_zz_0_zz_xxxyzz, g_zz_0_zz_xxxzzz, g_zz_0_zz_xxyyyy, g_zz_0_zz_xxyyyz, g_zz_0_zz_xxyyzz, g_zz_0_zz_xxyzzz, g_zz_0_zz_xxzzzz, g_zz_0_zz_xyyyyy, g_zz_0_zz_xyyyyz, g_zz_0_zz_xyyyzz, g_zz_0_zz_xyyzzz, g_zz_0_zz_xyzzzz, g_zz_0_zz_xzzzzz, g_zz_0_zz_yyyyyy, g_zz_0_zz_yyyyyz, g_zz_0_zz_yyyyzz, g_zz_0_zz_yyyzzz, g_zz_0_zz_yyzzzz, g_zz_0_zz_yzzzzz, g_zz_0_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zz_xxxxxx[k] = -2.0 * g_z_0_z_xxxxxx[k] - g_zz_0_z_xxxxxx[k] * ab_z + g_zz_0_z_xxxxxxz[k];

                g_zz_0_zz_xxxxxy[k] = -2.0 * g_z_0_z_xxxxxy[k] - g_zz_0_z_xxxxxy[k] * ab_z + g_zz_0_z_xxxxxyz[k];

                g_zz_0_zz_xxxxxz[k] = -2.0 * g_z_0_z_xxxxxz[k] - g_zz_0_z_xxxxxz[k] * ab_z + g_zz_0_z_xxxxxzz[k];

                g_zz_0_zz_xxxxyy[k] = -2.0 * g_z_0_z_xxxxyy[k] - g_zz_0_z_xxxxyy[k] * ab_z + g_zz_0_z_xxxxyyz[k];

                g_zz_0_zz_xxxxyz[k] = -2.0 * g_z_0_z_xxxxyz[k] - g_zz_0_z_xxxxyz[k] * ab_z + g_zz_0_z_xxxxyzz[k];

                g_zz_0_zz_xxxxzz[k] = -2.0 * g_z_0_z_xxxxzz[k] - g_zz_0_z_xxxxzz[k] * ab_z + g_zz_0_z_xxxxzzz[k];

                g_zz_0_zz_xxxyyy[k] = -2.0 * g_z_0_z_xxxyyy[k] - g_zz_0_z_xxxyyy[k] * ab_z + g_zz_0_z_xxxyyyz[k];

                g_zz_0_zz_xxxyyz[k] = -2.0 * g_z_0_z_xxxyyz[k] - g_zz_0_z_xxxyyz[k] * ab_z + g_zz_0_z_xxxyyzz[k];

                g_zz_0_zz_xxxyzz[k] = -2.0 * g_z_0_z_xxxyzz[k] - g_zz_0_z_xxxyzz[k] * ab_z + g_zz_0_z_xxxyzzz[k];

                g_zz_0_zz_xxxzzz[k] = -2.0 * g_z_0_z_xxxzzz[k] - g_zz_0_z_xxxzzz[k] * ab_z + g_zz_0_z_xxxzzzz[k];

                g_zz_0_zz_xxyyyy[k] = -2.0 * g_z_0_z_xxyyyy[k] - g_zz_0_z_xxyyyy[k] * ab_z + g_zz_0_z_xxyyyyz[k];

                g_zz_0_zz_xxyyyz[k] = -2.0 * g_z_0_z_xxyyyz[k] - g_zz_0_z_xxyyyz[k] * ab_z + g_zz_0_z_xxyyyzz[k];

                g_zz_0_zz_xxyyzz[k] = -2.0 * g_z_0_z_xxyyzz[k] - g_zz_0_z_xxyyzz[k] * ab_z + g_zz_0_z_xxyyzzz[k];

                g_zz_0_zz_xxyzzz[k] = -2.0 * g_z_0_z_xxyzzz[k] - g_zz_0_z_xxyzzz[k] * ab_z + g_zz_0_z_xxyzzzz[k];

                g_zz_0_zz_xxzzzz[k] = -2.0 * g_z_0_z_xxzzzz[k] - g_zz_0_z_xxzzzz[k] * ab_z + g_zz_0_z_xxzzzzz[k];

                g_zz_0_zz_xyyyyy[k] = -2.0 * g_z_0_z_xyyyyy[k] - g_zz_0_z_xyyyyy[k] * ab_z + g_zz_0_z_xyyyyyz[k];

                g_zz_0_zz_xyyyyz[k] = -2.0 * g_z_0_z_xyyyyz[k] - g_zz_0_z_xyyyyz[k] * ab_z + g_zz_0_z_xyyyyzz[k];

                g_zz_0_zz_xyyyzz[k] = -2.0 * g_z_0_z_xyyyzz[k] - g_zz_0_z_xyyyzz[k] * ab_z + g_zz_0_z_xyyyzzz[k];

                g_zz_0_zz_xyyzzz[k] = -2.0 * g_z_0_z_xyyzzz[k] - g_zz_0_z_xyyzzz[k] * ab_z + g_zz_0_z_xyyzzzz[k];

                g_zz_0_zz_xyzzzz[k] = -2.0 * g_z_0_z_xyzzzz[k] - g_zz_0_z_xyzzzz[k] * ab_z + g_zz_0_z_xyzzzzz[k];

                g_zz_0_zz_xzzzzz[k] = -2.0 * g_z_0_z_xzzzzz[k] - g_zz_0_z_xzzzzz[k] * ab_z + g_zz_0_z_xzzzzzz[k];

                g_zz_0_zz_yyyyyy[k] = -2.0 * g_z_0_z_yyyyyy[k] - g_zz_0_z_yyyyyy[k] * ab_z + g_zz_0_z_yyyyyyz[k];

                g_zz_0_zz_yyyyyz[k] = -2.0 * g_z_0_z_yyyyyz[k] - g_zz_0_z_yyyyyz[k] * ab_z + g_zz_0_z_yyyyyzz[k];

                g_zz_0_zz_yyyyzz[k] = -2.0 * g_z_0_z_yyyyzz[k] - g_zz_0_z_yyyyzz[k] * ab_z + g_zz_0_z_yyyyzzz[k];

                g_zz_0_zz_yyyzzz[k] = -2.0 * g_z_0_z_yyyzzz[k] - g_zz_0_z_yyyzzz[k] * ab_z + g_zz_0_z_yyyzzzz[k];

                g_zz_0_zz_yyzzzz[k] = -2.0 * g_z_0_z_yyzzzz[k] - g_zz_0_z_yyzzzz[k] * ab_z + g_zz_0_z_yyzzzzz[k];

                g_zz_0_zz_yzzzzz[k] = -2.0 * g_z_0_z_yzzzzz[k] - g_zz_0_z_yzzzzz[k] * ab_z + g_zz_0_z_yzzzzzz[k];

                g_zz_0_zz_zzzzzz[k] = -2.0 * g_z_0_z_zzzzzz[k] - g_zz_0_z_zzzzzz[k] * ab_z + g_zz_0_z_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

