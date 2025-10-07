#include "ElectronRepulsionGeom1100ContrRecDIXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_dixx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_dixx,
                                            const size_t idx_geom_01_pixx,
                                            const size_t idx_geom_10_pixx,
                                            const size_t idx_geom_11_pixx,
                                            const size_t idx_geom_11_pkxx,
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

            const auto pi_geom_01_off = idx_geom_01_pixx + i * dcomps + j;

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

            const auto pi_geom_11_off = idx_geom_11_pixx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : PKSS

            const auto pk_geom_11_off = idx_geom_11_pkxx + i * dcomps + j;

            auto g_x_x_x_xxxxxxx = cbuffer.data(pk_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xxxxxxy = cbuffer.data(pk_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xxxxxxz = cbuffer.data(pk_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_xxxxxyy = cbuffer.data(pk_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_xxxxxyz = cbuffer.data(pk_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_xxxxxzz = cbuffer.data(pk_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_x_xxxxyyy = cbuffer.data(pk_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_x_xxxxyyz = cbuffer.data(pk_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_x_xxxxyzz = cbuffer.data(pk_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_x_xxxxzzz = cbuffer.data(pk_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_x_xxxyyyy = cbuffer.data(pk_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_x_xxxyyyz = cbuffer.data(pk_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_x_xxxyyzz = cbuffer.data(pk_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_x_xxxyzzz = cbuffer.data(pk_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_x_xxxzzzz = cbuffer.data(pk_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_x_xxyyyyy = cbuffer.data(pk_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_x_xxyyyyz = cbuffer.data(pk_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_x_xxyyyzz = cbuffer.data(pk_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_x_xxyyzzz = cbuffer.data(pk_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_x_xxyzzzz = cbuffer.data(pk_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_x_xxzzzzz = cbuffer.data(pk_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_x_xyyyyyy = cbuffer.data(pk_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_x_xyyyyyz = cbuffer.data(pk_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_x_xyyyyzz = cbuffer.data(pk_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_x_xyyyzzz = cbuffer.data(pk_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_x_xyyzzzz = cbuffer.data(pk_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_x_xyzzzzz = cbuffer.data(pk_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_x_xzzzzzz = cbuffer.data(pk_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_x_yyyyyyy = cbuffer.data(pk_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_x_yyyyyyz = cbuffer.data(pk_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_x_yyyyyzz = cbuffer.data(pk_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_x_yyyyzzz = cbuffer.data(pk_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_x_yyyzzzz = cbuffer.data(pk_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_x_yyzzzzz = cbuffer.data(pk_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_x_yzzzzzz = cbuffer.data(pk_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_x_zzzzzzz = cbuffer.data(pk_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_y_xxxxxxx = cbuffer.data(pk_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_y_xxxxxxy = cbuffer.data(pk_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_y_xxxxxxz = cbuffer.data(pk_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_y_xxxxxyy = cbuffer.data(pk_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_y_xxxxxyz = cbuffer.data(pk_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_y_xxxxxzz = cbuffer.data(pk_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_y_xxxxyyy = cbuffer.data(pk_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_y_xxxxyyz = cbuffer.data(pk_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_y_xxxxyzz = cbuffer.data(pk_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_y_xxxxzzz = cbuffer.data(pk_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_y_xxxyyyy = cbuffer.data(pk_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_y_xxxyyyz = cbuffer.data(pk_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_y_xxxyyzz = cbuffer.data(pk_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_y_xxxyzzz = cbuffer.data(pk_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_y_xxxzzzz = cbuffer.data(pk_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_y_xxyyyyy = cbuffer.data(pk_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_y_xxyyyyz = cbuffer.data(pk_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_y_xxyyyzz = cbuffer.data(pk_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_y_xxyyzzz = cbuffer.data(pk_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_y_xxyzzzz = cbuffer.data(pk_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_y_xxzzzzz = cbuffer.data(pk_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_y_xyyyyyy = cbuffer.data(pk_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_y_xyyyyyz = cbuffer.data(pk_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_y_xyyyyzz = cbuffer.data(pk_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_y_xyyyzzz = cbuffer.data(pk_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_y_xyyzzzz = cbuffer.data(pk_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_y_xyzzzzz = cbuffer.data(pk_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_y_xzzzzzz = cbuffer.data(pk_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_y_yyyyyyy = cbuffer.data(pk_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_y_yyyyyyz = cbuffer.data(pk_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_y_yyyyyzz = cbuffer.data(pk_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_y_yyyyzzz = cbuffer.data(pk_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_y_yyyzzzz = cbuffer.data(pk_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_y_yyzzzzz = cbuffer.data(pk_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_y_yzzzzzz = cbuffer.data(pk_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_y_zzzzzzz = cbuffer.data(pk_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_z_xxxxxxx = cbuffer.data(pk_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_z_xxxxxxy = cbuffer.data(pk_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_z_xxxxxxz = cbuffer.data(pk_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_z_xxxxxyy = cbuffer.data(pk_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_z_xxxxxyz = cbuffer.data(pk_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_z_xxxxxzz = cbuffer.data(pk_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_z_xxxxyyy = cbuffer.data(pk_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_z_xxxxyyz = cbuffer.data(pk_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_z_xxxxyzz = cbuffer.data(pk_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_z_xxxxzzz = cbuffer.data(pk_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_z_xxxyyyy = cbuffer.data(pk_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_z_xxxyyyz = cbuffer.data(pk_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_x_z_xxxyyzz = cbuffer.data(pk_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_z_xxxyzzz = cbuffer.data(pk_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_z_xxxzzzz = cbuffer.data(pk_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_z_xxyyyyy = cbuffer.data(pk_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_z_xxyyyyz = cbuffer.data(pk_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_z_xxyyyzz = cbuffer.data(pk_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_x_z_xxyyzzz = cbuffer.data(pk_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_x_z_xxyzzzz = cbuffer.data(pk_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_x_z_xxzzzzz = cbuffer.data(pk_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_x_z_xyyyyyy = cbuffer.data(pk_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_x_z_xyyyyyz = cbuffer.data(pk_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_x_z_xyyyyzz = cbuffer.data(pk_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_x_z_xyyyzzz = cbuffer.data(pk_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_x_z_xyyzzzz = cbuffer.data(pk_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_x_z_xyzzzzz = cbuffer.data(pk_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_x_z_xzzzzzz = cbuffer.data(pk_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_x_z_yyyyyyy = cbuffer.data(pk_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_x_z_yyyyyyz = cbuffer.data(pk_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_x_z_yyyyyzz = cbuffer.data(pk_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_x_z_yyyyzzz = cbuffer.data(pk_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_x_z_yyyzzzz = cbuffer.data(pk_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_x_z_yyzzzzz = cbuffer.data(pk_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_x_z_yzzzzzz = cbuffer.data(pk_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_x_z_zzzzzzz = cbuffer.data(pk_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_x_xxxxxxx = cbuffer.data(pk_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_x_xxxxxxy = cbuffer.data(pk_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_x_xxxxxxz = cbuffer.data(pk_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_x_xxxxxyy = cbuffer.data(pk_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_x_xxxxxyz = cbuffer.data(pk_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_x_xxxxxzz = cbuffer.data(pk_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_x_xxxxyyy = cbuffer.data(pk_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_x_xxxxyyz = cbuffer.data(pk_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_x_xxxxyzz = cbuffer.data(pk_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_x_xxxxzzz = cbuffer.data(pk_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_x_xxxyyyy = cbuffer.data(pk_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_x_xxxyyyz = cbuffer.data(pk_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_y_x_xxxyyzz = cbuffer.data(pk_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_y_x_xxxyzzz = cbuffer.data(pk_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_y_x_xxxzzzz = cbuffer.data(pk_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_y_x_xxyyyyy = cbuffer.data(pk_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_y_x_xxyyyyz = cbuffer.data(pk_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_y_x_xxyyyzz = cbuffer.data(pk_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_y_x_xxyyzzz = cbuffer.data(pk_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_y_x_xxyzzzz = cbuffer.data(pk_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_y_x_xxzzzzz = cbuffer.data(pk_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_y_x_xyyyyyy = cbuffer.data(pk_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_y_x_xyyyyyz = cbuffer.data(pk_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_y_x_xyyyyzz = cbuffer.data(pk_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_y_x_xyyyzzz = cbuffer.data(pk_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_y_x_xyyzzzz = cbuffer.data(pk_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_y_x_xyzzzzz = cbuffer.data(pk_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_y_x_xzzzzzz = cbuffer.data(pk_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_y_x_yyyyyyy = cbuffer.data(pk_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_y_x_yyyyyyz = cbuffer.data(pk_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_y_x_yyyyyzz = cbuffer.data(pk_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_y_x_yyyyzzz = cbuffer.data(pk_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_y_x_yyyzzzz = cbuffer.data(pk_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_y_x_yyzzzzz = cbuffer.data(pk_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_y_x_yzzzzzz = cbuffer.data(pk_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_y_x_zzzzzzz = cbuffer.data(pk_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_y_y_xxxxxxx = cbuffer.data(pk_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_y_y_xxxxxxy = cbuffer.data(pk_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_y_y_xxxxxxz = cbuffer.data(pk_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_y_y_xxxxxyy = cbuffer.data(pk_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_y_y_xxxxxyz = cbuffer.data(pk_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_y_y_xxxxxzz = cbuffer.data(pk_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_y_y_xxxxyyy = cbuffer.data(pk_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_y_xxxxyyz = cbuffer.data(pk_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_y_xxxxyzz = cbuffer.data(pk_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_y_xxxxzzz = cbuffer.data(pk_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_y_xxxyyyy = cbuffer.data(pk_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_y_xxxyyyz = cbuffer.data(pk_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_y_y_xxxyyzz = cbuffer.data(pk_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_y_xxxyzzz = cbuffer.data(pk_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_y_xxxzzzz = cbuffer.data(pk_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_y_xxyyyyy = cbuffer.data(pk_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_y_xxyyyyz = cbuffer.data(pk_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_y_xxyyyzz = cbuffer.data(pk_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_y_y_xxyyzzz = cbuffer.data(pk_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_y_xxyzzzz = cbuffer.data(pk_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_y_xxzzzzz = cbuffer.data(pk_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_y_y_xyyyyyy = cbuffer.data(pk_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_y_xyyyyyz = cbuffer.data(pk_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_y_xyyyyzz = cbuffer.data(pk_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_y_y_xyyyzzz = cbuffer.data(pk_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_y_xyyzzzz = cbuffer.data(pk_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_y_xyzzzzz = cbuffer.data(pk_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_y_xzzzzzz = cbuffer.data(pk_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_y_yyyyyyy = cbuffer.data(pk_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_y_yyyyyyz = cbuffer.data(pk_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_y_y_yyyyyzz = cbuffer.data(pk_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_y_yyyyzzz = cbuffer.data(pk_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_y_yyyzzzz = cbuffer.data(pk_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_y_yyzzzzz = cbuffer.data(pk_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_y_yzzzzzz = cbuffer.data(pk_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_y_zzzzzzz = cbuffer.data(pk_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_y_z_xxxxxxx = cbuffer.data(pk_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_y_z_xxxxxxy = cbuffer.data(pk_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_y_z_xxxxxxz = cbuffer.data(pk_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_y_z_xxxxxyy = cbuffer.data(pk_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_y_z_xxxxxyz = cbuffer.data(pk_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_y_z_xxxxxzz = cbuffer.data(pk_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_y_z_xxxxyyy = cbuffer.data(pk_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_y_z_xxxxyyz = cbuffer.data(pk_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_y_z_xxxxyzz = cbuffer.data(pk_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_y_z_xxxxzzz = cbuffer.data(pk_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_y_z_xxxyyyy = cbuffer.data(pk_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_y_z_xxxyyyz = cbuffer.data(pk_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_y_z_xxxyyzz = cbuffer.data(pk_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_y_z_xxxyzzz = cbuffer.data(pk_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_y_z_xxxzzzz = cbuffer.data(pk_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_y_z_xxyyyyy = cbuffer.data(pk_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_y_z_xxyyyyz = cbuffer.data(pk_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_y_z_xxyyyzz = cbuffer.data(pk_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_y_z_xxyyzzz = cbuffer.data(pk_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_y_z_xxyzzzz = cbuffer.data(pk_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_y_z_xxzzzzz = cbuffer.data(pk_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_y_z_xyyyyyy = cbuffer.data(pk_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_y_z_xyyyyyz = cbuffer.data(pk_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_y_z_xyyyyzz = cbuffer.data(pk_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_y_z_xyyyzzz = cbuffer.data(pk_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_y_z_xyyzzzz = cbuffer.data(pk_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_y_z_xyzzzzz = cbuffer.data(pk_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_y_z_xzzzzzz = cbuffer.data(pk_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_y_z_yyyyyyy = cbuffer.data(pk_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_y_z_yyyyyyz = cbuffer.data(pk_geom_11_off + 209 * ccomps * dcomps);

            auto g_x_y_z_yyyyyzz = cbuffer.data(pk_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_y_z_yyyyzzz = cbuffer.data(pk_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_y_z_yyyzzzz = cbuffer.data(pk_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_y_z_yyzzzzz = cbuffer.data(pk_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_y_z_yzzzzzz = cbuffer.data(pk_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_y_z_zzzzzzz = cbuffer.data(pk_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_z_x_xxxxxxx = cbuffer.data(pk_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_z_x_xxxxxxy = cbuffer.data(pk_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_z_x_xxxxxxz = cbuffer.data(pk_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_z_x_xxxxxyy = cbuffer.data(pk_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_z_x_xxxxxyz = cbuffer.data(pk_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_z_x_xxxxxzz = cbuffer.data(pk_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_z_x_xxxxyyy = cbuffer.data(pk_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_z_x_xxxxyyz = cbuffer.data(pk_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_z_x_xxxxyzz = cbuffer.data(pk_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_z_x_xxxxzzz = cbuffer.data(pk_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_z_x_xxxyyyy = cbuffer.data(pk_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_z_x_xxxyyyz = cbuffer.data(pk_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_z_x_xxxyyzz = cbuffer.data(pk_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_z_x_xxxyzzz = cbuffer.data(pk_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_z_x_xxxzzzz = cbuffer.data(pk_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_z_x_xxyyyyy = cbuffer.data(pk_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_z_x_xxyyyyz = cbuffer.data(pk_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_z_x_xxyyyzz = cbuffer.data(pk_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_z_x_xxyyzzz = cbuffer.data(pk_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_z_x_xxyzzzz = cbuffer.data(pk_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_z_x_xxzzzzz = cbuffer.data(pk_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_z_x_xyyyyyy = cbuffer.data(pk_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_z_x_xyyyyyz = cbuffer.data(pk_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_z_x_xyyyyzz = cbuffer.data(pk_geom_11_off + 239 * ccomps * dcomps);

            auto g_x_z_x_xyyyzzz = cbuffer.data(pk_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_z_x_xyyzzzz = cbuffer.data(pk_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_z_x_xyzzzzz = cbuffer.data(pk_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_z_x_xzzzzzz = cbuffer.data(pk_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_z_x_yyyyyyy = cbuffer.data(pk_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_z_x_yyyyyyz = cbuffer.data(pk_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_z_x_yyyyyzz = cbuffer.data(pk_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_z_x_yyyyzzz = cbuffer.data(pk_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_z_x_yyyzzzz = cbuffer.data(pk_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_z_x_yyzzzzz = cbuffer.data(pk_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_z_x_yzzzzzz = cbuffer.data(pk_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_z_x_zzzzzzz = cbuffer.data(pk_geom_11_off + 251 * ccomps * dcomps);

            auto g_x_z_y_xxxxxxx = cbuffer.data(pk_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_z_y_xxxxxxy = cbuffer.data(pk_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_z_y_xxxxxxz = cbuffer.data(pk_geom_11_off + 254 * ccomps * dcomps);

            auto g_x_z_y_xxxxxyy = cbuffer.data(pk_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_z_y_xxxxxyz = cbuffer.data(pk_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_z_y_xxxxxzz = cbuffer.data(pk_geom_11_off + 257 * ccomps * dcomps);

            auto g_x_z_y_xxxxyyy = cbuffer.data(pk_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_z_y_xxxxyyz = cbuffer.data(pk_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_z_y_xxxxyzz = cbuffer.data(pk_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_z_y_xxxxzzz = cbuffer.data(pk_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_z_y_xxxyyyy = cbuffer.data(pk_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_z_y_xxxyyyz = cbuffer.data(pk_geom_11_off + 263 * ccomps * dcomps);

            auto g_x_z_y_xxxyyzz = cbuffer.data(pk_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_z_y_xxxyzzz = cbuffer.data(pk_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_z_y_xxxzzzz = cbuffer.data(pk_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_z_y_xxyyyyy = cbuffer.data(pk_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_z_y_xxyyyyz = cbuffer.data(pk_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_z_y_xxyyyzz = cbuffer.data(pk_geom_11_off + 269 * ccomps * dcomps);

            auto g_x_z_y_xxyyzzz = cbuffer.data(pk_geom_11_off + 270 * ccomps * dcomps);

            auto g_x_z_y_xxyzzzz = cbuffer.data(pk_geom_11_off + 271 * ccomps * dcomps);

            auto g_x_z_y_xxzzzzz = cbuffer.data(pk_geom_11_off + 272 * ccomps * dcomps);

            auto g_x_z_y_xyyyyyy = cbuffer.data(pk_geom_11_off + 273 * ccomps * dcomps);

            auto g_x_z_y_xyyyyyz = cbuffer.data(pk_geom_11_off + 274 * ccomps * dcomps);

            auto g_x_z_y_xyyyyzz = cbuffer.data(pk_geom_11_off + 275 * ccomps * dcomps);

            auto g_x_z_y_xyyyzzz = cbuffer.data(pk_geom_11_off + 276 * ccomps * dcomps);

            auto g_x_z_y_xyyzzzz = cbuffer.data(pk_geom_11_off + 277 * ccomps * dcomps);

            auto g_x_z_y_xyzzzzz = cbuffer.data(pk_geom_11_off + 278 * ccomps * dcomps);

            auto g_x_z_y_xzzzzzz = cbuffer.data(pk_geom_11_off + 279 * ccomps * dcomps);

            auto g_x_z_y_yyyyyyy = cbuffer.data(pk_geom_11_off + 280 * ccomps * dcomps);

            auto g_x_z_y_yyyyyyz = cbuffer.data(pk_geom_11_off + 281 * ccomps * dcomps);

            auto g_x_z_y_yyyyyzz = cbuffer.data(pk_geom_11_off + 282 * ccomps * dcomps);

            auto g_x_z_y_yyyyzzz = cbuffer.data(pk_geom_11_off + 283 * ccomps * dcomps);

            auto g_x_z_y_yyyzzzz = cbuffer.data(pk_geom_11_off + 284 * ccomps * dcomps);

            auto g_x_z_y_yyzzzzz = cbuffer.data(pk_geom_11_off + 285 * ccomps * dcomps);

            auto g_x_z_y_yzzzzzz = cbuffer.data(pk_geom_11_off + 286 * ccomps * dcomps);

            auto g_x_z_y_zzzzzzz = cbuffer.data(pk_geom_11_off + 287 * ccomps * dcomps);

            auto g_x_z_z_xxxxxxx = cbuffer.data(pk_geom_11_off + 288 * ccomps * dcomps);

            auto g_x_z_z_xxxxxxy = cbuffer.data(pk_geom_11_off + 289 * ccomps * dcomps);

            auto g_x_z_z_xxxxxxz = cbuffer.data(pk_geom_11_off + 290 * ccomps * dcomps);

            auto g_x_z_z_xxxxxyy = cbuffer.data(pk_geom_11_off + 291 * ccomps * dcomps);

            auto g_x_z_z_xxxxxyz = cbuffer.data(pk_geom_11_off + 292 * ccomps * dcomps);

            auto g_x_z_z_xxxxxzz = cbuffer.data(pk_geom_11_off + 293 * ccomps * dcomps);

            auto g_x_z_z_xxxxyyy = cbuffer.data(pk_geom_11_off + 294 * ccomps * dcomps);

            auto g_x_z_z_xxxxyyz = cbuffer.data(pk_geom_11_off + 295 * ccomps * dcomps);

            auto g_x_z_z_xxxxyzz = cbuffer.data(pk_geom_11_off + 296 * ccomps * dcomps);

            auto g_x_z_z_xxxxzzz = cbuffer.data(pk_geom_11_off + 297 * ccomps * dcomps);

            auto g_x_z_z_xxxyyyy = cbuffer.data(pk_geom_11_off + 298 * ccomps * dcomps);

            auto g_x_z_z_xxxyyyz = cbuffer.data(pk_geom_11_off + 299 * ccomps * dcomps);

            auto g_x_z_z_xxxyyzz = cbuffer.data(pk_geom_11_off + 300 * ccomps * dcomps);

            auto g_x_z_z_xxxyzzz = cbuffer.data(pk_geom_11_off + 301 * ccomps * dcomps);

            auto g_x_z_z_xxxzzzz = cbuffer.data(pk_geom_11_off + 302 * ccomps * dcomps);

            auto g_x_z_z_xxyyyyy = cbuffer.data(pk_geom_11_off + 303 * ccomps * dcomps);

            auto g_x_z_z_xxyyyyz = cbuffer.data(pk_geom_11_off + 304 * ccomps * dcomps);

            auto g_x_z_z_xxyyyzz = cbuffer.data(pk_geom_11_off + 305 * ccomps * dcomps);

            auto g_x_z_z_xxyyzzz = cbuffer.data(pk_geom_11_off + 306 * ccomps * dcomps);

            auto g_x_z_z_xxyzzzz = cbuffer.data(pk_geom_11_off + 307 * ccomps * dcomps);

            auto g_x_z_z_xxzzzzz = cbuffer.data(pk_geom_11_off + 308 * ccomps * dcomps);

            auto g_x_z_z_xyyyyyy = cbuffer.data(pk_geom_11_off + 309 * ccomps * dcomps);

            auto g_x_z_z_xyyyyyz = cbuffer.data(pk_geom_11_off + 310 * ccomps * dcomps);

            auto g_x_z_z_xyyyyzz = cbuffer.data(pk_geom_11_off + 311 * ccomps * dcomps);

            auto g_x_z_z_xyyyzzz = cbuffer.data(pk_geom_11_off + 312 * ccomps * dcomps);

            auto g_x_z_z_xyyzzzz = cbuffer.data(pk_geom_11_off + 313 * ccomps * dcomps);

            auto g_x_z_z_xyzzzzz = cbuffer.data(pk_geom_11_off + 314 * ccomps * dcomps);

            auto g_x_z_z_xzzzzzz = cbuffer.data(pk_geom_11_off + 315 * ccomps * dcomps);

            auto g_x_z_z_yyyyyyy = cbuffer.data(pk_geom_11_off + 316 * ccomps * dcomps);

            auto g_x_z_z_yyyyyyz = cbuffer.data(pk_geom_11_off + 317 * ccomps * dcomps);

            auto g_x_z_z_yyyyyzz = cbuffer.data(pk_geom_11_off + 318 * ccomps * dcomps);

            auto g_x_z_z_yyyyzzz = cbuffer.data(pk_geom_11_off + 319 * ccomps * dcomps);

            auto g_x_z_z_yyyzzzz = cbuffer.data(pk_geom_11_off + 320 * ccomps * dcomps);

            auto g_x_z_z_yyzzzzz = cbuffer.data(pk_geom_11_off + 321 * ccomps * dcomps);

            auto g_x_z_z_yzzzzzz = cbuffer.data(pk_geom_11_off + 322 * ccomps * dcomps);

            auto g_x_z_z_zzzzzzz = cbuffer.data(pk_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_x_x_xxxxxxx = cbuffer.data(pk_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_x_x_xxxxxxy = cbuffer.data(pk_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_x_x_xxxxxxz = cbuffer.data(pk_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_x_x_xxxxxyy = cbuffer.data(pk_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_x_x_xxxxxyz = cbuffer.data(pk_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_x_x_xxxxxzz = cbuffer.data(pk_geom_11_off + 329 * ccomps * dcomps);

            auto g_y_x_x_xxxxyyy = cbuffer.data(pk_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_x_x_xxxxyyz = cbuffer.data(pk_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_x_x_xxxxyzz = cbuffer.data(pk_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_x_x_xxxxzzz = cbuffer.data(pk_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_x_x_xxxyyyy = cbuffer.data(pk_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_x_x_xxxyyyz = cbuffer.data(pk_geom_11_off + 335 * ccomps * dcomps);

            auto g_y_x_x_xxxyyzz = cbuffer.data(pk_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_x_x_xxxyzzz = cbuffer.data(pk_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_x_x_xxxzzzz = cbuffer.data(pk_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_x_x_xxyyyyy = cbuffer.data(pk_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_x_x_xxyyyyz = cbuffer.data(pk_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_x_x_xxyyyzz = cbuffer.data(pk_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_x_x_xxyyzzz = cbuffer.data(pk_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_x_x_xxyzzzz = cbuffer.data(pk_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_x_x_xxzzzzz = cbuffer.data(pk_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_x_x_xyyyyyy = cbuffer.data(pk_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_x_x_xyyyyyz = cbuffer.data(pk_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_x_x_xyyyyzz = cbuffer.data(pk_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_x_x_xyyyzzz = cbuffer.data(pk_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_x_x_xyyzzzz = cbuffer.data(pk_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_x_x_xyzzzzz = cbuffer.data(pk_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_x_x_xzzzzzz = cbuffer.data(pk_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_x_x_yyyyyyy = cbuffer.data(pk_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_x_x_yyyyyyz = cbuffer.data(pk_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_x_x_yyyyyzz = cbuffer.data(pk_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_x_x_yyyyzzz = cbuffer.data(pk_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_x_x_yyyzzzz = cbuffer.data(pk_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_x_x_yyzzzzz = cbuffer.data(pk_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_x_x_yzzzzzz = cbuffer.data(pk_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_x_x_zzzzzzz = cbuffer.data(pk_geom_11_off + 359 * ccomps * dcomps);

            auto g_y_x_y_xxxxxxx = cbuffer.data(pk_geom_11_off + 360 * ccomps * dcomps);

            auto g_y_x_y_xxxxxxy = cbuffer.data(pk_geom_11_off + 361 * ccomps * dcomps);

            auto g_y_x_y_xxxxxxz = cbuffer.data(pk_geom_11_off + 362 * ccomps * dcomps);

            auto g_y_x_y_xxxxxyy = cbuffer.data(pk_geom_11_off + 363 * ccomps * dcomps);

            auto g_y_x_y_xxxxxyz = cbuffer.data(pk_geom_11_off + 364 * ccomps * dcomps);

            auto g_y_x_y_xxxxxzz = cbuffer.data(pk_geom_11_off + 365 * ccomps * dcomps);

            auto g_y_x_y_xxxxyyy = cbuffer.data(pk_geom_11_off + 366 * ccomps * dcomps);

            auto g_y_x_y_xxxxyyz = cbuffer.data(pk_geom_11_off + 367 * ccomps * dcomps);

            auto g_y_x_y_xxxxyzz = cbuffer.data(pk_geom_11_off + 368 * ccomps * dcomps);

            auto g_y_x_y_xxxxzzz = cbuffer.data(pk_geom_11_off + 369 * ccomps * dcomps);

            auto g_y_x_y_xxxyyyy = cbuffer.data(pk_geom_11_off + 370 * ccomps * dcomps);

            auto g_y_x_y_xxxyyyz = cbuffer.data(pk_geom_11_off + 371 * ccomps * dcomps);

            auto g_y_x_y_xxxyyzz = cbuffer.data(pk_geom_11_off + 372 * ccomps * dcomps);

            auto g_y_x_y_xxxyzzz = cbuffer.data(pk_geom_11_off + 373 * ccomps * dcomps);

            auto g_y_x_y_xxxzzzz = cbuffer.data(pk_geom_11_off + 374 * ccomps * dcomps);

            auto g_y_x_y_xxyyyyy = cbuffer.data(pk_geom_11_off + 375 * ccomps * dcomps);

            auto g_y_x_y_xxyyyyz = cbuffer.data(pk_geom_11_off + 376 * ccomps * dcomps);

            auto g_y_x_y_xxyyyzz = cbuffer.data(pk_geom_11_off + 377 * ccomps * dcomps);

            auto g_y_x_y_xxyyzzz = cbuffer.data(pk_geom_11_off + 378 * ccomps * dcomps);

            auto g_y_x_y_xxyzzzz = cbuffer.data(pk_geom_11_off + 379 * ccomps * dcomps);

            auto g_y_x_y_xxzzzzz = cbuffer.data(pk_geom_11_off + 380 * ccomps * dcomps);

            auto g_y_x_y_xyyyyyy = cbuffer.data(pk_geom_11_off + 381 * ccomps * dcomps);

            auto g_y_x_y_xyyyyyz = cbuffer.data(pk_geom_11_off + 382 * ccomps * dcomps);

            auto g_y_x_y_xyyyyzz = cbuffer.data(pk_geom_11_off + 383 * ccomps * dcomps);

            auto g_y_x_y_xyyyzzz = cbuffer.data(pk_geom_11_off + 384 * ccomps * dcomps);

            auto g_y_x_y_xyyzzzz = cbuffer.data(pk_geom_11_off + 385 * ccomps * dcomps);

            auto g_y_x_y_xyzzzzz = cbuffer.data(pk_geom_11_off + 386 * ccomps * dcomps);

            auto g_y_x_y_xzzzzzz = cbuffer.data(pk_geom_11_off + 387 * ccomps * dcomps);

            auto g_y_x_y_yyyyyyy = cbuffer.data(pk_geom_11_off + 388 * ccomps * dcomps);

            auto g_y_x_y_yyyyyyz = cbuffer.data(pk_geom_11_off + 389 * ccomps * dcomps);

            auto g_y_x_y_yyyyyzz = cbuffer.data(pk_geom_11_off + 390 * ccomps * dcomps);

            auto g_y_x_y_yyyyzzz = cbuffer.data(pk_geom_11_off + 391 * ccomps * dcomps);

            auto g_y_x_y_yyyzzzz = cbuffer.data(pk_geom_11_off + 392 * ccomps * dcomps);

            auto g_y_x_y_yyzzzzz = cbuffer.data(pk_geom_11_off + 393 * ccomps * dcomps);

            auto g_y_x_y_yzzzzzz = cbuffer.data(pk_geom_11_off + 394 * ccomps * dcomps);

            auto g_y_x_y_zzzzzzz = cbuffer.data(pk_geom_11_off + 395 * ccomps * dcomps);

            auto g_y_x_z_xxxxxxx = cbuffer.data(pk_geom_11_off + 396 * ccomps * dcomps);

            auto g_y_x_z_xxxxxxy = cbuffer.data(pk_geom_11_off + 397 * ccomps * dcomps);

            auto g_y_x_z_xxxxxxz = cbuffer.data(pk_geom_11_off + 398 * ccomps * dcomps);

            auto g_y_x_z_xxxxxyy = cbuffer.data(pk_geom_11_off + 399 * ccomps * dcomps);

            auto g_y_x_z_xxxxxyz = cbuffer.data(pk_geom_11_off + 400 * ccomps * dcomps);

            auto g_y_x_z_xxxxxzz = cbuffer.data(pk_geom_11_off + 401 * ccomps * dcomps);

            auto g_y_x_z_xxxxyyy = cbuffer.data(pk_geom_11_off + 402 * ccomps * dcomps);

            auto g_y_x_z_xxxxyyz = cbuffer.data(pk_geom_11_off + 403 * ccomps * dcomps);

            auto g_y_x_z_xxxxyzz = cbuffer.data(pk_geom_11_off + 404 * ccomps * dcomps);

            auto g_y_x_z_xxxxzzz = cbuffer.data(pk_geom_11_off + 405 * ccomps * dcomps);

            auto g_y_x_z_xxxyyyy = cbuffer.data(pk_geom_11_off + 406 * ccomps * dcomps);

            auto g_y_x_z_xxxyyyz = cbuffer.data(pk_geom_11_off + 407 * ccomps * dcomps);

            auto g_y_x_z_xxxyyzz = cbuffer.data(pk_geom_11_off + 408 * ccomps * dcomps);

            auto g_y_x_z_xxxyzzz = cbuffer.data(pk_geom_11_off + 409 * ccomps * dcomps);

            auto g_y_x_z_xxxzzzz = cbuffer.data(pk_geom_11_off + 410 * ccomps * dcomps);

            auto g_y_x_z_xxyyyyy = cbuffer.data(pk_geom_11_off + 411 * ccomps * dcomps);

            auto g_y_x_z_xxyyyyz = cbuffer.data(pk_geom_11_off + 412 * ccomps * dcomps);

            auto g_y_x_z_xxyyyzz = cbuffer.data(pk_geom_11_off + 413 * ccomps * dcomps);

            auto g_y_x_z_xxyyzzz = cbuffer.data(pk_geom_11_off + 414 * ccomps * dcomps);

            auto g_y_x_z_xxyzzzz = cbuffer.data(pk_geom_11_off + 415 * ccomps * dcomps);

            auto g_y_x_z_xxzzzzz = cbuffer.data(pk_geom_11_off + 416 * ccomps * dcomps);

            auto g_y_x_z_xyyyyyy = cbuffer.data(pk_geom_11_off + 417 * ccomps * dcomps);

            auto g_y_x_z_xyyyyyz = cbuffer.data(pk_geom_11_off + 418 * ccomps * dcomps);

            auto g_y_x_z_xyyyyzz = cbuffer.data(pk_geom_11_off + 419 * ccomps * dcomps);

            auto g_y_x_z_xyyyzzz = cbuffer.data(pk_geom_11_off + 420 * ccomps * dcomps);

            auto g_y_x_z_xyyzzzz = cbuffer.data(pk_geom_11_off + 421 * ccomps * dcomps);

            auto g_y_x_z_xyzzzzz = cbuffer.data(pk_geom_11_off + 422 * ccomps * dcomps);

            auto g_y_x_z_xzzzzzz = cbuffer.data(pk_geom_11_off + 423 * ccomps * dcomps);

            auto g_y_x_z_yyyyyyy = cbuffer.data(pk_geom_11_off + 424 * ccomps * dcomps);

            auto g_y_x_z_yyyyyyz = cbuffer.data(pk_geom_11_off + 425 * ccomps * dcomps);

            auto g_y_x_z_yyyyyzz = cbuffer.data(pk_geom_11_off + 426 * ccomps * dcomps);

            auto g_y_x_z_yyyyzzz = cbuffer.data(pk_geom_11_off + 427 * ccomps * dcomps);

            auto g_y_x_z_yyyzzzz = cbuffer.data(pk_geom_11_off + 428 * ccomps * dcomps);

            auto g_y_x_z_yyzzzzz = cbuffer.data(pk_geom_11_off + 429 * ccomps * dcomps);

            auto g_y_x_z_yzzzzzz = cbuffer.data(pk_geom_11_off + 430 * ccomps * dcomps);

            auto g_y_x_z_zzzzzzz = cbuffer.data(pk_geom_11_off + 431 * ccomps * dcomps);

            auto g_y_y_x_xxxxxxx = cbuffer.data(pk_geom_11_off + 432 * ccomps * dcomps);

            auto g_y_y_x_xxxxxxy = cbuffer.data(pk_geom_11_off + 433 * ccomps * dcomps);

            auto g_y_y_x_xxxxxxz = cbuffer.data(pk_geom_11_off + 434 * ccomps * dcomps);

            auto g_y_y_x_xxxxxyy = cbuffer.data(pk_geom_11_off + 435 * ccomps * dcomps);

            auto g_y_y_x_xxxxxyz = cbuffer.data(pk_geom_11_off + 436 * ccomps * dcomps);

            auto g_y_y_x_xxxxxzz = cbuffer.data(pk_geom_11_off + 437 * ccomps * dcomps);

            auto g_y_y_x_xxxxyyy = cbuffer.data(pk_geom_11_off + 438 * ccomps * dcomps);

            auto g_y_y_x_xxxxyyz = cbuffer.data(pk_geom_11_off + 439 * ccomps * dcomps);

            auto g_y_y_x_xxxxyzz = cbuffer.data(pk_geom_11_off + 440 * ccomps * dcomps);

            auto g_y_y_x_xxxxzzz = cbuffer.data(pk_geom_11_off + 441 * ccomps * dcomps);

            auto g_y_y_x_xxxyyyy = cbuffer.data(pk_geom_11_off + 442 * ccomps * dcomps);

            auto g_y_y_x_xxxyyyz = cbuffer.data(pk_geom_11_off + 443 * ccomps * dcomps);

            auto g_y_y_x_xxxyyzz = cbuffer.data(pk_geom_11_off + 444 * ccomps * dcomps);

            auto g_y_y_x_xxxyzzz = cbuffer.data(pk_geom_11_off + 445 * ccomps * dcomps);

            auto g_y_y_x_xxxzzzz = cbuffer.data(pk_geom_11_off + 446 * ccomps * dcomps);

            auto g_y_y_x_xxyyyyy = cbuffer.data(pk_geom_11_off + 447 * ccomps * dcomps);

            auto g_y_y_x_xxyyyyz = cbuffer.data(pk_geom_11_off + 448 * ccomps * dcomps);

            auto g_y_y_x_xxyyyzz = cbuffer.data(pk_geom_11_off + 449 * ccomps * dcomps);

            auto g_y_y_x_xxyyzzz = cbuffer.data(pk_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_y_x_xxyzzzz = cbuffer.data(pk_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_y_x_xxzzzzz = cbuffer.data(pk_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_y_x_xyyyyyy = cbuffer.data(pk_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_y_x_xyyyyyz = cbuffer.data(pk_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_y_x_xyyyyzz = cbuffer.data(pk_geom_11_off + 455 * ccomps * dcomps);

            auto g_y_y_x_xyyyzzz = cbuffer.data(pk_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_y_x_xyyzzzz = cbuffer.data(pk_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_y_x_xyzzzzz = cbuffer.data(pk_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_y_x_xzzzzzz = cbuffer.data(pk_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_y_x_yyyyyyy = cbuffer.data(pk_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_y_x_yyyyyyz = cbuffer.data(pk_geom_11_off + 461 * ccomps * dcomps);

            auto g_y_y_x_yyyyyzz = cbuffer.data(pk_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_y_x_yyyyzzz = cbuffer.data(pk_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_y_x_yyyzzzz = cbuffer.data(pk_geom_11_off + 464 * ccomps * dcomps);

            auto g_y_y_x_yyzzzzz = cbuffer.data(pk_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_y_x_yzzzzzz = cbuffer.data(pk_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_y_x_zzzzzzz = cbuffer.data(pk_geom_11_off + 467 * ccomps * dcomps);

            auto g_y_y_y_xxxxxxx = cbuffer.data(pk_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_y_y_xxxxxxy = cbuffer.data(pk_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_y_y_xxxxxxz = cbuffer.data(pk_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_y_y_xxxxxyy = cbuffer.data(pk_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_y_y_xxxxxyz = cbuffer.data(pk_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_y_y_xxxxxzz = cbuffer.data(pk_geom_11_off + 473 * ccomps * dcomps);

            auto g_y_y_y_xxxxyyy = cbuffer.data(pk_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_y_y_xxxxyyz = cbuffer.data(pk_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_y_y_xxxxyzz = cbuffer.data(pk_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_y_y_xxxxzzz = cbuffer.data(pk_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_y_y_xxxyyyy = cbuffer.data(pk_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_y_y_xxxyyyz = cbuffer.data(pk_geom_11_off + 479 * ccomps * dcomps);

            auto g_y_y_y_xxxyyzz = cbuffer.data(pk_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_y_y_xxxyzzz = cbuffer.data(pk_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_y_y_xxxzzzz = cbuffer.data(pk_geom_11_off + 482 * ccomps * dcomps);

            auto g_y_y_y_xxyyyyy = cbuffer.data(pk_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_y_y_xxyyyyz = cbuffer.data(pk_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_y_y_xxyyyzz = cbuffer.data(pk_geom_11_off + 485 * ccomps * dcomps);

            auto g_y_y_y_xxyyzzz = cbuffer.data(pk_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_y_y_xxyzzzz = cbuffer.data(pk_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_y_y_xxzzzzz = cbuffer.data(pk_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_y_y_xyyyyyy = cbuffer.data(pk_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_y_y_xyyyyyz = cbuffer.data(pk_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_y_y_xyyyyzz = cbuffer.data(pk_geom_11_off + 491 * ccomps * dcomps);

            auto g_y_y_y_xyyyzzz = cbuffer.data(pk_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_y_y_xyyzzzz = cbuffer.data(pk_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_y_y_xyzzzzz = cbuffer.data(pk_geom_11_off + 494 * ccomps * dcomps);

            auto g_y_y_y_xzzzzzz = cbuffer.data(pk_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_y_y_yyyyyyy = cbuffer.data(pk_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_y_y_yyyyyyz = cbuffer.data(pk_geom_11_off + 497 * ccomps * dcomps);

            auto g_y_y_y_yyyyyzz = cbuffer.data(pk_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_y_y_yyyyzzz = cbuffer.data(pk_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_y_y_yyyzzzz = cbuffer.data(pk_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_y_y_yyzzzzz = cbuffer.data(pk_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_y_y_yzzzzzz = cbuffer.data(pk_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_y_y_zzzzzzz = cbuffer.data(pk_geom_11_off + 503 * ccomps * dcomps);

            auto g_y_y_z_xxxxxxx = cbuffer.data(pk_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_y_z_xxxxxxy = cbuffer.data(pk_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_y_z_xxxxxxz = cbuffer.data(pk_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_y_z_xxxxxyy = cbuffer.data(pk_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_y_z_xxxxxyz = cbuffer.data(pk_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_y_z_xxxxxzz = cbuffer.data(pk_geom_11_off + 509 * ccomps * dcomps);

            auto g_y_y_z_xxxxyyy = cbuffer.data(pk_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_y_z_xxxxyyz = cbuffer.data(pk_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_y_z_xxxxyzz = cbuffer.data(pk_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_y_z_xxxxzzz = cbuffer.data(pk_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_y_z_xxxyyyy = cbuffer.data(pk_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_y_z_xxxyyyz = cbuffer.data(pk_geom_11_off + 515 * ccomps * dcomps);

            auto g_y_y_z_xxxyyzz = cbuffer.data(pk_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_y_z_xxxyzzz = cbuffer.data(pk_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_y_z_xxxzzzz = cbuffer.data(pk_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_y_z_xxyyyyy = cbuffer.data(pk_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_y_z_xxyyyyz = cbuffer.data(pk_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_y_z_xxyyyzz = cbuffer.data(pk_geom_11_off + 521 * ccomps * dcomps);

            auto g_y_y_z_xxyyzzz = cbuffer.data(pk_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_y_z_xxyzzzz = cbuffer.data(pk_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_y_z_xxzzzzz = cbuffer.data(pk_geom_11_off + 524 * ccomps * dcomps);

            auto g_y_y_z_xyyyyyy = cbuffer.data(pk_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_y_z_xyyyyyz = cbuffer.data(pk_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_y_z_xyyyyzz = cbuffer.data(pk_geom_11_off + 527 * ccomps * dcomps);

            auto g_y_y_z_xyyyzzz = cbuffer.data(pk_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_y_z_xyyzzzz = cbuffer.data(pk_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_y_z_xyzzzzz = cbuffer.data(pk_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_y_z_xzzzzzz = cbuffer.data(pk_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_y_z_yyyyyyy = cbuffer.data(pk_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_y_z_yyyyyyz = cbuffer.data(pk_geom_11_off + 533 * ccomps * dcomps);

            auto g_y_y_z_yyyyyzz = cbuffer.data(pk_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_y_z_yyyyzzz = cbuffer.data(pk_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_y_z_yyyzzzz = cbuffer.data(pk_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_y_z_yyzzzzz = cbuffer.data(pk_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_y_z_yzzzzzz = cbuffer.data(pk_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_y_z_zzzzzzz = cbuffer.data(pk_geom_11_off + 539 * ccomps * dcomps);

            auto g_y_z_x_xxxxxxx = cbuffer.data(pk_geom_11_off + 540 * ccomps * dcomps);

            auto g_y_z_x_xxxxxxy = cbuffer.data(pk_geom_11_off + 541 * ccomps * dcomps);

            auto g_y_z_x_xxxxxxz = cbuffer.data(pk_geom_11_off + 542 * ccomps * dcomps);

            auto g_y_z_x_xxxxxyy = cbuffer.data(pk_geom_11_off + 543 * ccomps * dcomps);

            auto g_y_z_x_xxxxxyz = cbuffer.data(pk_geom_11_off + 544 * ccomps * dcomps);

            auto g_y_z_x_xxxxxzz = cbuffer.data(pk_geom_11_off + 545 * ccomps * dcomps);

            auto g_y_z_x_xxxxyyy = cbuffer.data(pk_geom_11_off + 546 * ccomps * dcomps);

            auto g_y_z_x_xxxxyyz = cbuffer.data(pk_geom_11_off + 547 * ccomps * dcomps);

            auto g_y_z_x_xxxxyzz = cbuffer.data(pk_geom_11_off + 548 * ccomps * dcomps);

            auto g_y_z_x_xxxxzzz = cbuffer.data(pk_geom_11_off + 549 * ccomps * dcomps);

            auto g_y_z_x_xxxyyyy = cbuffer.data(pk_geom_11_off + 550 * ccomps * dcomps);

            auto g_y_z_x_xxxyyyz = cbuffer.data(pk_geom_11_off + 551 * ccomps * dcomps);

            auto g_y_z_x_xxxyyzz = cbuffer.data(pk_geom_11_off + 552 * ccomps * dcomps);

            auto g_y_z_x_xxxyzzz = cbuffer.data(pk_geom_11_off + 553 * ccomps * dcomps);

            auto g_y_z_x_xxxzzzz = cbuffer.data(pk_geom_11_off + 554 * ccomps * dcomps);

            auto g_y_z_x_xxyyyyy = cbuffer.data(pk_geom_11_off + 555 * ccomps * dcomps);

            auto g_y_z_x_xxyyyyz = cbuffer.data(pk_geom_11_off + 556 * ccomps * dcomps);

            auto g_y_z_x_xxyyyzz = cbuffer.data(pk_geom_11_off + 557 * ccomps * dcomps);

            auto g_y_z_x_xxyyzzz = cbuffer.data(pk_geom_11_off + 558 * ccomps * dcomps);

            auto g_y_z_x_xxyzzzz = cbuffer.data(pk_geom_11_off + 559 * ccomps * dcomps);

            auto g_y_z_x_xxzzzzz = cbuffer.data(pk_geom_11_off + 560 * ccomps * dcomps);

            auto g_y_z_x_xyyyyyy = cbuffer.data(pk_geom_11_off + 561 * ccomps * dcomps);

            auto g_y_z_x_xyyyyyz = cbuffer.data(pk_geom_11_off + 562 * ccomps * dcomps);

            auto g_y_z_x_xyyyyzz = cbuffer.data(pk_geom_11_off + 563 * ccomps * dcomps);

            auto g_y_z_x_xyyyzzz = cbuffer.data(pk_geom_11_off + 564 * ccomps * dcomps);

            auto g_y_z_x_xyyzzzz = cbuffer.data(pk_geom_11_off + 565 * ccomps * dcomps);

            auto g_y_z_x_xyzzzzz = cbuffer.data(pk_geom_11_off + 566 * ccomps * dcomps);

            auto g_y_z_x_xzzzzzz = cbuffer.data(pk_geom_11_off + 567 * ccomps * dcomps);

            auto g_y_z_x_yyyyyyy = cbuffer.data(pk_geom_11_off + 568 * ccomps * dcomps);

            auto g_y_z_x_yyyyyyz = cbuffer.data(pk_geom_11_off + 569 * ccomps * dcomps);

            auto g_y_z_x_yyyyyzz = cbuffer.data(pk_geom_11_off + 570 * ccomps * dcomps);

            auto g_y_z_x_yyyyzzz = cbuffer.data(pk_geom_11_off + 571 * ccomps * dcomps);

            auto g_y_z_x_yyyzzzz = cbuffer.data(pk_geom_11_off + 572 * ccomps * dcomps);

            auto g_y_z_x_yyzzzzz = cbuffer.data(pk_geom_11_off + 573 * ccomps * dcomps);

            auto g_y_z_x_yzzzzzz = cbuffer.data(pk_geom_11_off + 574 * ccomps * dcomps);

            auto g_y_z_x_zzzzzzz = cbuffer.data(pk_geom_11_off + 575 * ccomps * dcomps);

            auto g_y_z_y_xxxxxxx = cbuffer.data(pk_geom_11_off + 576 * ccomps * dcomps);

            auto g_y_z_y_xxxxxxy = cbuffer.data(pk_geom_11_off + 577 * ccomps * dcomps);

            auto g_y_z_y_xxxxxxz = cbuffer.data(pk_geom_11_off + 578 * ccomps * dcomps);

            auto g_y_z_y_xxxxxyy = cbuffer.data(pk_geom_11_off + 579 * ccomps * dcomps);

            auto g_y_z_y_xxxxxyz = cbuffer.data(pk_geom_11_off + 580 * ccomps * dcomps);

            auto g_y_z_y_xxxxxzz = cbuffer.data(pk_geom_11_off + 581 * ccomps * dcomps);

            auto g_y_z_y_xxxxyyy = cbuffer.data(pk_geom_11_off + 582 * ccomps * dcomps);

            auto g_y_z_y_xxxxyyz = cbuffer.data(pk_geom_11_off + 583 * ccomps * dcomps);

            auto g_y_z_y_xxxxyzz = cbuffer.data(pk_geom_11_off + 584 * ccomps * dcomps);

            auto g_y_z_y_xxxxzzz = cbuffer.data(pk_geom_11_off + 585 * ccomps * dcomps);

            auto g_y_z_y_xxxyyyy = cbuffer.data(pk_geom_11_off + 586 * ccomps * dcomps);

            auto g_y_z_y_xxxyyyz = cbuffer.data(pk_geom_11_off + 587 * ccomps * dcomps);

            auto g_y_z_y_xxxyyzz = cbuffer.data(pk_geom_11_off + 588 * ccomps * dcomps);

            auto g_y_z_y_xxxyzzz = cbuffer.data(pk_geom_11_off + 589 * ccomps * dcomps);

            auto g_y_z_y_xxxzzzz = cbuffer.data(pk_geom_11_off + 590 * ccomps * dcomps);

            auto g_y_z_y_xxyyyyy = cbuffer.data(pk_geom_11_off + 591 * ccomps * dcomps);

            auto g_y_z_y_xxyyyyz = cbuffer.data(pk_geom_11_off + 592 * ccomps * dcomps);

            auto g_y_z_y_xxyyyzz = cbuffer.data(pk_geom_11_off + 593 * ccomps * dcomps);

            auto g_y_z_y_xxyyzzz = cbuffer.data(pk_geom_11_off + 594 * ccomps * dcomps);

            auto g_y_z_y_xxyzzzz = cbuffer.data(pk_geom_11_off + 595 * ccomps * dcomps);

            auto g_y_z_y_xxzzzzz = cbuffer.data(pk_geom_11_off + 596 * ccomps * dcomps);

            auto g_y_z_y_xyyyyyy = cbuffer.data(pk_geom_11_off + 597 * ccomps * dcomps);

            auto g_y_z_y_xyyyyyz = cbuffer.data(pk_geom_11_off + 598 * ccomps * dcomps);

            auto g_y_z_y_xyyyyzz = cbuffer.data(pk_geom_11_off + 599 * ccomps * dcomps);

            auto g_y_z_y_xyyyzzz = cbuffer.data(pk_geom_11_off + 600 * ccomps * dcomps);

            auto g_y_z_y_xyyzzzz = cbuffer.data(pk_geom_11_off + 601 * ccomps * dcomps);

            auto g_y_z_y_xyzzzzz = cbuffer.data(pk_geom_11_off + 602 * ccomps * dcomps);

            auto g_y_z_y_xzzzzzz = cbuffer.data(pk_geom_11_off + 603 * ccomps * dcomps);

            auto g_y_z_y_yyyyyyy = cbuffer.data(pk_geom_11_off + 604 * ccomps * dcomps);

            auto g_y_z_y_yyyyyyz = cbuffer.data(pk_geom_11_off + 605 * ccomps * dcomps);

            auto g_y_z_y_yyyyyzz = cbuffer.data(pk_geom_11_off + 606 * ccomps * dcomps);

            auto g_y_z_y_yyyyzzz = cbuffer.data(pk_geom_11_off + 607 * ccomps * dcomps);

            auto g_y_z_y_yyyzzzz = cbuffer.data(pk_geom_11_off + 608 * ccomps * dcomps);

            auto g_y_z_y_yyzzzzz = cbuffer.data(pk_geom_11_off + 609 * ccomps * dcomps);

            auto g_y_z_y_yzzzzzz = cbuffer.data(pk_geom_11_off + 610 * ccomps * dcomps);

            auto g_y_z_y_zzzzzzz = cbuffer.data(pk_geom_11_off + 611 * ccomps * dcomps);

            auto g_y_z_z_xxxxxxx = cbuffer.data(pk_geom_11_off + 612 * ccomps * dcomps);

            auto g_y_z_z_xxxxxxy = cbuffer.data(pk_geom_11_off + 613 * ccomps * dcomps);

            auto g_y_z_z_xxxxxxz = cbuffer.data(pk_geom_11_off + 614 * ccomps * dcomps);

            auto g_y_z_z_xxxxxyy = cbuffer.data(pk_geom_11_off + 615 * ccomps * dcomps);

            auto g_y_z_z_xxxxxyz = cbuffer.data(pk_geom_11_off + 616 * ccomps * dcomps);

            auto g_y_z_z_xxxxxzz = cbuffer.data(pk_geom_11_off + 617 * ccomps * dcomps);

            auto g_y_z_z_xxxxyyy = cbuffer.data(pk_geom_11_off + 618 * ccomps * dcomps);

            auto g_y_z_z_xxxxyyz = cbuffer.data(pk_geom_11_off + 619 * ccomps * dcomps);

            auto g_y_z_z_xxxxyzz = cbuffer.data(pk_geom_11_off + 620 * ccomps * dcomps);

            auto g_y_z_z_xxxxzzz = cbuffer.data(pk_geom_11_off + 621 * ccomps * dcomps);

            auto g_y_z_z_xxxyyyy = cbuffer.data(pk_geom_11_off + 622 * ccomps * dcomps);

            auto g_y_z_z_xxxyyyz = cbuffer.data(pk_geom_11_off + 623 * ccomps * dcomps);

            auto g_y_z_z_xxxyyzz = cbuffer.data(pk_geom_11_off + 624 * ccomps * dcomps);

            auto g_y_z_z_xxxyzzz = cbuffer.data(pk_geom_11_off + 625 * ccomps * dcomps);

            auto g_y_z_z_xxxzzzz = cbuffer.data(pk_geom_11_off + 626 * ccomps * dcomps);

            auto g_y_z_z_xxyyyyy = cbuffer.data(pk_geom_11_off + 627 * ccomps * dcomps);

            auto g_y_z_z_xxyyyyz = cbuffer.data(pk_geom_11_off + 628 * ccomps * dcomps);

            auto g_y_z_z_xxyyyzz = cbuffer.data(pk_geom_11_off + 629 * ccomps * dcomps);

            auto g_y_z_z_xxyyzzz = cbuffer.data(pk_geom_11_off + 630 * ccomps * dcomps);

            auto g_y_z_z_xxyzzzz = cbuffer.data(pk_geom_11_off + 631 * ccomps * dcomps);

            auto g_y_z_z_xxzzzzz = cbuffer.data(pk_geom_11_off + 632 * ccomps * dcomps);

            auto g_y_z_z_xyyyyyy = cbuffer.data(pk_geom_11_off + 633 * ccomps * dcomps);

            auto g_y_z_z_xyyyyyz = cbuffer.data(pk_geom_11_off + 634 * ccomps * dcomps);

            auto g_y_z_z_xyyyyzz = cbuffer.data(pk_geom_11_off + 635 * ccomps * dcomps);

            auto g_y_z_z_xyyyzzz = cbuffer.data(pk_geom_11_off + 636 * ccomps * dcomps);

            auto g_y_z_z_xyyzzzz = cbuffer.data(pk_geom_11_off + 637 * ccomps * dcomps);

            auto g_y_z_z_xyzzzzz = cbuffer.data(pk_geom_11_off + 638 * ccomps * dcomps);

            auto g_y_z_z_xzzzzzz = cbuffer.data(pk_geom_11_off + 639 * ccomps * dcomps);

            auto g_y_z_z_yyyyyyy = cbuffer.data(pk_geom_11_off + 640 * ccomps * dcomps);

            auto g_y_z_z_yyyyyyz = cbuffer.data(pk_geom_11_off + 641 * ccomps * dcomps);

            auto g_y_z_z_yyyyyzz = cbuffer.data(pk_geom_11_off + 642 * ccomps * dcomps);

            auto g_y_z_z_yyyyzzz = cbuffer.data(pk_geom_11_off + 643 * ccomps * dcomps);

            auto g_y_z_z_yyyzzzz = cbuffer.data(pk_geom_11_off + 644 * ccomps * dcomps);

            auto g_y_z_z_yyzzzzz = cbuffer.data(pk_geom_11_off + 645 * ccomps * dcomps);

            auto g_y_z_z_yzzzzzz = cbuffer.data(pk_geom_11_off + 646 * ccomps * dcomps);

            auto g_y_z_z_zzzzzzz = cbuffer.data(pk_geom_11_off + 647 * ccomps * dcomps);

            auto g_z_x_x_xxxxxxx = cbuffer.data(pk_geom_11_off + 648 * ccomps * dcomps);

            auto g_z_x_x_xxxxxxy = cbuffer.data(pk_geom_11_off + 649 * ccomps * dcomps);

            auto g_z_x_x_xxxxxxz = cbuffer.data(pk_geom_11_off + 650 * ccomps * dcomps);

            auto g_z_x_x_xxxxxyy = cbuffer.data(pk_geom_11_off + 651 * ccomps * dcomps);

            auto g_z_x_x_xxxxxyz = cbuffer.data(pk_geom_11_off + 652 * ccomps * dcomps);

            auto g_z_x_x_xxxxxzz = cbuffer.data(pk_geom_11_off + 653 * ccomps * dcomps);

            auto g_z_x_x_xxxxyyy = cbuffer.data(pk_geom_11_off + 654 * ccomps * dcomps);

            auto g_z_x_x_xxxxyyz = cbuffer.data(pk_geom_11_off + 655 * ccomps * dcomps);

            auto g_z_x_x_xxxxyzz = cbuffer.data(pk_geom_11_off + 656 * ccomps * dcomps);

            auto g_z_x_x_xxxxzzz = cbuffer.data(pk_geom_11_off + 657 * ccomps * dcomps);

            auto g_z_x_x_xxxyyyy = cbuffer.data(pk_geom_11_off + 658 * ccomps * dcomps);

            auto g_z_x_x_xxxyyyz = cbuffer.data(pk_geom_11_off + 659 * ccomps * dcomps);

            auto g_z_x_x_xxxyyzz = cbuffer.data(pk_geom_11_off + 660 * ccomps * dcomps);

            auto g_z_x_x_xxxyzzz = cbuffer.data(pk_geom_11_off + 661 * ccomps * dcomps);

            auto g_z_x_x_xxxzzzz = cbuffer.data(pk_geom_11_off + 662 * ccomps * dcomps);

            auto g_z_x_x_xxyyyyy = cbuffer.data(pk_geom_11_off + 663 * ccomps * dcomps);

            auto g_z_x_x_xxyyyyz = cbuffer.data(pk_geom_11_off + 664 * ccomps * dcomps);

            auto g_z_x_x_xxyyyzz = cbuffer.data(pk_geom_11_off + 665 * ccomps * dcomps);

            auto g_z_x_x_xxyyzzz = cbuffer.data(pk_geom_11_off + 666 * ccomps * dcomps);

            auto g_z_x_x_xxyzzzz = cbuffer.data(pk_geom_11_off + 667 * ccomps * dcomps);

            auto g_z_x_x_xxzzzzz = cbuffer.data(pk_geom_11_off + 668 * ccomps * dcomps);

            auto g_z_x_x_xyyyyyy = cbuffer.data(pk_geom_11_off + 669 * ccomps * dcomps);

            auto g_z_x_x_xyyyyyz = cbuffer.data(pk_geom_11_off + 670 * ccomps * dcomps);

            auto g_z_x_x_xyyyyzz = cbuffer.data(pk_geom_11_off + 671 * ccomps * dcomps);

            auto g_z_x_x_xyyyzzz = cbuffer.data(pk_geom_11_off + 672 * ccomps * dcomps);

            auto g_z_x_x_xyyzzzz = cbuffer.data(pk_geom_11_off + 673 * ccomps * dcomps);

            auto g_z_x_x_xyzzzzz = cbuffer.data(pk_geom_11_off + 674 * ccomps * dcomps);

            auto g_z_x_x_xzzzzzz = cbuffer.data(pk_geom_11_off + 675 * ccomps * dcomps);

            auto g_z_x_x_yyyyyyy = cbuffer.data(pk_geom_11_off + 676 * ccomps * dcomps);

            auto g_z_x_x_yyyyyyz = cbuffer.data(pk_geom_11_off + 677 * ccomps * dcomps);

            auto g_z_x_x_yyyyyzz = cbuffer.data(pk_geom_11_off + 678 * ccomps * dcomps);

            auto g_z_x_x_yyyyzzz = cbuffer.data(pk_geom_11_off + 679 * ccomps * dcomps);

            auto g_z_x_x_yyyzzzz = cbuffer.data(pk_geom_11_off + 680 * ccomps * dcomps);

            auto g_z_x_x_yyzzzzz = cbuffer.data(pk_geom_11_off + 681 * ccomps * dcomps);

            auto g_z_x_x_yzzzzzz = cbuffer.data(pk_geom_11_off + 682 * ccomps * dcomps);

            auto g_z_x_x_zzzzzzz = cbuffer.data(pk_geom_11_off + 683 * ccomps * dcomps);

            auto g_z_x_y_xxxxxxx = cbuffer.data(pk_geom_11_off + 684 * ccomps * dcomps);

            auto g_z_x_y_xxxxxxy = cbuffer.data(pk_geom_11_off + 685 * ccomps * dcomps);

            auto g_z_x_y_xxxxxxz = cbuffer.data(pk_geom_11_off + 686 * ccomps * dcomps);

            auto g_z_x_y_xxxxxyy = cbuffer.data(pk_geom_11_off + 687 * ccomps * dcomps);

            auto g_z_x_y_xxxxxyz = cbuffer.data(pk_geom_11_off + 688 * ccomps * dcomps);

            auto g_z_x_y_xxxxxzz = cbuffer.data(pk_geom_11_off + 689 * ccomps * dcomps);

            auto g_z_x_y_xxxxyyy = cbuffer.data(pk_geom_11_off + 690 * ccomps * dcomps);

            auto g_z_x_y_xxxxyyz = cbuffer.data(pk_geom_11_off + 691 * ccomps * dcomps);

            auto g_z_x_y_xxxxyzz = cbuffer.data(pk_geom_11_off + 692 * ccomps * dcomps);

            auto g_z_x_y_xxxxzzz = cbuffer.data(pk_geom_11_off + 693 * ccomps * dcomps);

            auto g_z_x_y_xxxyyyy = cbuffer.data(pk_geom_11_off + 694 * ccomps * dcomps);

            auto g_z_x_y_xxxyyyz = cbuffer.data(pk_geom_11_off + 695 * ccomps * dcomps);

            auto g_z_x_y_xxxyyzz = cbuffer.data(pk_geom_11_off + 696 * ccomps * dcomps);

            auto g_z_x_y_xxxyzzz = cbuffer.data(pk_geom_11_off + 697 * ccomps * dcomps);

            auto g_z_x_y_xxxzzzz = cbuffer.data(pk_geom_11_off + 698 * ccomps * dcomps);

            auto g_z_x_y_xxyyyyy = cbuffer.data(pk_geom_11_off + 699 * ccomps * dcomps);

            auto g_z_x_y_xxyyyyz = cbuffer.data(pk_geom_11_off + 700 * ccomps * dcomps);

            auto g_z_x_y_xxyyyzz = cbuffer.data(pk_geom_11_off + 701 * ccomps * dcomps);

            auto g_z_x_y_xxyyzzz = cbuffer.data(pk_geom_11_off + 702 * ccomps * dcomps);

            auto g_z_x_y_xxyzzzz = cbuffer.data(pk_geom_11_off + 703 * ccomps * dcomps);

            auto g_z_x_y_xxzzzzz = cbuffer.data(pk_geom_11_off + 704 * ccomps * dcomps);

            auto g_z_x_y_xyyyyyy = cbuffer.data(pk_geom_11_off + 705 * ccomps * dcomps);

            auto g_z_x_y_xyyyyyz = cbuffer.data(pk_geom_11_off + 706 * ccomps * dcomps);

            auto g_z_x_y_xyyyyzz = cbuffer.data(pk_geom_11_off + 707 * ccomps * dcomps);

            auto g_z_x_y_xyyyzzz = cbuffer.data(pk_geom_11_off + 708 * ccomps * dcomps);

            auto g_z_x_y_xyyzzzz = cbuffer.data(pk_geom_11_off + 709 * ccomps * dcomps);

            auto g_z_x_y_xyzzzzz = cbuffer.data(pk_geom_11_off + 710 * ccomps * dcomps);

            auto g_z_x_y_xzzzzzz = cbuffer.data(pk_geom_11_off + 711 * ccomps * dcomps);

            auto g_z_x_y_yyyyyyy = cbuffer.data(pk_geom_11_off + 712 * ccomps * dcomps);

            auto g_z_x_y_yyyyyyz = cbuffer.data(pk_geom_11_off + 713 * ccomps * dcomps);

            auto g_z_x_y_yyyyyzz = cbuffer.data(pk_geom_11_off + 714 * ccomps * dcomps);

            auto g_z_x_y_yyyyzzz = cbuffer.data(pk_geom_11_off + 715 * ccomps * dcomps);

            auto g_z_x_y_yyyzzzz = cbuffer.data(pk_geom_11_off + 716 * ccomps * dcomps);

            auto g_z_x_y_yyzzzzz = cbuffer.data(pk_geom_11_off + 717 * ccomps * dcomps);

            auto g_z_x_y_yzzzzzz = cbuffer.data(pk_geom_11_off + 718 * ccomps * dcomps);

            auto g_z_x_y_zzzzzzz = cbuffer.data(pk_geom_11_off + 719 * ccomps * dcomps);

            auto g_z_x_z_xxxxxxx = cbuffer.data(pk_geom_11_off + 720 * ccomps * dcomps);

            auto g_z_x_z_xxxxxxy = cbuffer.data(pk_geom_11_off + 721 * ccomps * dcomps);

            auto g_z_x_z_xxxxxxz = cbuffer.data(pk_geom_11_off + 722 * ccomps * dcomps);

            auto g_z_x_z_xxxxxyy = cbuffer.data(pk_geom_11_off + 723 * ccomps * dcomps);

            auto g_z_x_z_xxxxxyz = cbuffer.data(pk_geom_11_off + 724 * ccomps * dcomps);

            auto g_z_x_z_xxxxxzz = cbuffer.data(pk_geom_11_off + 725 * ccomps * dcomps);

            auto g_z_x_z_xxxxyyy = cbuffer.data(pk_geom_11_off + 726 * ccomps * dcomps);

            auto g_z_x_z_xxxxyyz = cbuffer.data(pk_geom_11_off + 727 * ccomps * dcomps);

            auto g_z_x_z_xxxxyzz = cbuffer.data(pk_geom_11_off + 728 * ccomps * dcomps);

            auto g_z_x_z_xxxxzzz = cbuffer.data(pk_geom_11_off + 729 * ccomps * dcomps);

            auto g_z_x_z_xxxyyyy = cbuffer.data(pk_geom_11_off + 730 * ccomps * dcomps);

            auto g_z_x_z_xxxyyyz = cbuffer.data(pk_geom_11_off + 731 * ccomps * dcomps);

            auto g_z_x_z_xxxyyzz = cbuffer.data(pk_geom_11_off + 732 * ccomps * dcomps);

            auto g_z_x_z_xxxyzzz = cbuffer.data(pk_geom_11_off + 733 * ccomps * dcomps);

            auto g_z_x_z_xxxzzzz = cbuffer.data(pk_geom_11_off + 734 * ccomps * dcomps);

            auto g_z_x_z_xxyyyyy = cbuffer.data(pk_geom_11_off + 735 * ccomps * dcomps);

            auto g_z_x_z_xxyyyyz = cbuffer.data(pk_geom_11_off + 736 * ccomps * dcomps);

            auto g_z_x_z_xxyyyzz = cbuffer.data(pk_geom_11_off + 737 * ccomps * dcomps);

            auto g_z_x_z_xxyyzzz = cbuffer.data(pk_geom_11_off + 738 * ccomps * dcomps);

            auto g_z_x_z_xxyzzzz = cbuffer.data(pk_geom_11_off + 739 * ccomps * dcomps);

            auto g_z_x_z_xxzzzzz = cbuffer.data(pk_geom_11_off + 740 * ccomps * dcomps);

            auto g_z_x_z_xyyyyyy = cbuffer.data(pk_geom_11_off + 741 * ccomps * dcomps);

            auto g_z_x_z_xyyyyyz = cbuffer.data(pk_geom_11_off + 742 * ccomps * dcomps);

            auto g_z_x_z_xyyyyzz = cbuffer.data(pk_geom_11_off + 743 * ccomps * dcomps);

            auto g_z_x_z_xyyyzzz = cbuffer.data(pk_geom_11_off + 744 * ccomps * dcomps);

            auto g_z_x_z_xyyzzzz = cbuffer.data(pk_geom_11_off + 745 * ccomps * dcomps);

            auto g_z_x_z_xyzzzzz = cbuffer.data(pk_geom_11_off + 746 * ccomps * dcomps);

            auto g_z_x_z_xzzzzzz = cbuffer.data(pk_geom_11_off + 747 * ccomps * dcomps);

            auto g_z_x_z_yyyyyyy = cbuffer.data(pk_geom_11_off + 748 * ccomps * dcomps);

            auto g_z_x_z_yyyyyyz = cbuffer.data(pk_geom_11_off + 749 * ccomps * dcomps);

            auto g_z_x_z_yyyyyzz = cbuffer.data(pk_geom_11_off + 750 * ccomps * dcomps);

            auto g_z_x_z_yyyyzzz = cbuffer.data(pk_geom_11_off + 751 * ccomps * dcomps);

            auto g_z_x_z_yyyzzzz = cbuffer.data(pk_geom_11_off + 752 * ccomps * dcomps);

            auto g_z_x_z_yyzzzzz = cbuffer.data(pk_geom_11_off + 753 * ccomps * dcomps);

            auto g_z_x_z_yzzzzzz = cbuffer.data(pk_geom_11_off + 754 * ccomps * dcomps);

            auto g_z_x_z_zzzzzzz = cbuffer.data(pk_geom_11_off + 755 * ccomps * dcomps);

            auto g_z_y_x_xxxxxxx = cbuffer.data(pk_geom_11_off + 756 * ccomps * dcomps);

            auto g_z_y_x_xxxxxxy = cbuffer.data(pk_geom_11_off + 757 * ccomps * dcomps);

            auto g_z_y_x_xxxxxxz = cbuffer.data(pk_geom_11_off + 758 * ccomps * dcomps);

            auto g_z_y_x_xxxxxyy = cbuffer.data(pk_geom_11_off + 759 * ccomps * dcomps);

            auto g_z_y_x_xxxxxyz = cbuffer.data(pk_geom_11_off + 760 * ccomps * dcomps);

            auto g_z_y_x_xxxxxzz = cbuffer.data(pk_geom_11_off + 761 * ccomps * dcomps);

            auto g_z_y_x_xxxxyyy = cbuffer.data(pk_geom_11_off + 762 * ccomps * dcomps);

            auto g_z_y_x_xxxxyyz = cbuffer.data(pk_geom_11_off + 763 * ccomps * dcomps);

            auto g_z_y_x_xxxxyzz = cbuffer.data(pk_geom_11_off + 764 * ccomps * dcomps);

            auto g_z_y_x_xxxxzzz = cbuffer.data(pk_geom_11_off + 765 * ccomps * dcomps);

            auto g_z_y_x_xxxyyyy = cbuffer.data(pk_geom_11_off + 766 * ccomps * dcomps);

            auto g_z_y_x_xxxyyyz = cbuffer.data(pk_geom_11_off + 767 * ccomps * dcomps);

            auto g_z_y_x_xxxyyzz = cbuffer.data(pk_geom_11_off + 768 * ccomps * dcomps);

            auto g_z_y_x_xxxyzzz = cbuffer.data(pk_geom_11_off + 769 * ccomps * dcomps);

            auto g_z_y_x_xxxzzzz = cbuffer.data(pk_geom_11_off + 770 * ccomps * dcomps);

            auto g_z_y_x_xxyyyyy = cbuffer.data(pk_geom_11_off + 771 * ccomps * dcomps);

            auto g_z_y_x_xxyyyyz = cbuffer.data(pk_geom_11_off + 772 * ccomps * dcomps);

            auto g_z_y_x_xxyyyzz = cbuffer.data(pk_geom_11_off + 773 * ccomps * dcomps);

            auto g_z_y_x_xxyyzzz = cbuffer.data(pk_geom_11_off + 774 * ccomps * dcomps);

            auto g_z_y_x_xxyzzzz = cbuffer.data(pk_geom_11_off + 775 * ccomps * dcomps);

            auto g_z_y_x_xxzzzzz = cbuffer.data(pk_geom_11_off + 776 * ccomps * dcomps);

            auto g_z_y_x_xyyyyyy = cbuffer.data(pk_geom_11_off + 777 * ccomps * dcomps);

            auto g_z_y_x_xyyyyyz = cbuffer.data(pk_geom_11_off + 778 * ccomps * dcomps);

            auto g_z_y_x_xyyyyzz = cbuffer.data(pk_geom_11_off + 779 * ccomps * dcomps);

            auto g_z_y_x_xyyyzzz = cbuffer.data(pk_geom_11_off + 780 * ccomps * dcomps);

            auto g_z_y_x_xyyzzzz = cbuffer.data(pk_geom_11_off + 781 * ccomps * dcomps);

            auto g_z_y_x_xyzzzzz = cbuffer.data(pk_geom_11_off + 782 * ccomps * dcomps);

            auto g_z_y_x_xzzzzzz = cbuffer.data(pk_geom_11_off + 783 * ccomps * dcomps);

            auto g_z_y_x_yyyyyyy = cbuffer.data(pk_geom_11_off + 784 * ccomps * dcomps);

            auto g_z_y_x_yyyyyyz = cbuffer.data(pk_geom_11_off + 785 * ccomps * dcomps);

            auto g_z_y_x_yyyyyzz = cbuffer.data(pk_geom_11_off + 786 * ccomps * dcomps);

            auto g_z_y_x_yyyyzzz = cbuffer.data(pk_geom_11_off + 787 * ccomps * dcomps);

            auto g_z_y_x_yyyzzzz = cbuffer.data(pk_geom_11_off + 788 * ccomps * dcomps);

            auto g_z_y_x_yyzzzzz = cbuffer.data(pk_geom_11_off + 789 * ccomps * dcomps);

            auto g_z_y_x_yzzzzzz = cbuffer.data(pk_geom_11_off + 790 * ccomps * dcomps);

            auto g_z_y_x_zzzzzzz = cbuffer.data(pk_geom_11_off + 791 * ccomps * dcomps);

            auto g_z_y_y_xxxxxxx = cbuffer.data(pk_geom_11_off + 792 * ccomps * dcomps);

            auto g_z_y_y_xxxxxxy = cbuffer.data(pk_geom_11_off + 793 * ccomps * dcomps);

            auto g_z_y_y_xxxxxxz = cbuffer.data(pk_geom_11_off + 794 * ccomps * dcomps);

            auto g_z_y_y_xxxxxyy = cbuffer.data(pk_geom_11_off + 795 * ccomps * dcomps);

            auto g_z_y_y_xxxxxyz = cbuffer.data(pk_geom_11_off + 796 * ccomps * dcomps);

            auto g_z_y_y_xxxxxzz = cbuffer.data(pk_geom_11_off + 797 * ccomps * dcomps);

            auto g_z_y_y_xxxxyyy = cbuffer.data(pk_geom_11_off + 798 * ccomps * dcomps);

            auto g_z_y_y_xxxxyyz = cbuffer.data(pk_geom_11_off + 799 * ccomps * dcomps);

            auto g_z_y_y_xxxxyzz = cbuffer.data(pk_geom_11_off + 800 * ccomps * dcomps);

            auto g_z_y_y_xxxxzzz = cbuffer.data(pk_geom_11_off + 801 * ccomps * dcomps);

            auto g_z_y_y_xxxyyyy = cbuffer.data(pk_geom_11_off + 802 * ccomps * dcomps);

            auto g_z_y_y_xxxyyyz = cbuffer.data(pk_geom_11_off + 803 * ccomps * dcomps);

            auto g_z_y_y_xxxyyzz = cbuffer.data(pk_geom_11_off + 804 * ccomps * dcomps);

            auto g_z_y_y_xxxyzzz = cbuffer.data(pk_geom_11_off + 805 * ccomps * dcomps);

            auto g_z_y_y_xxxzzzz = cbuffer.data(pk_geom_11_off + 806 * ccomps * dcomps);

            auto g_z_y_y_xxyyyyy = cbuffer.data(pk_geom_11_off + 807 * ccomps * dcomps);

            auto g_z_y_y_xxyyyyz = cbuffer.data(pk_geom_11_off + 808 * ccomps * dcomps);

            auto g_z_y_y_xxyyyzz = cbuffer.data(pk_geom_11_off + 809 * ccomps * dcomps);

            auto g_z_y_y_xxyyzzz = cbuffer.data(pk_geom_11_off + 810 * ccomps * dcomps);

            auto g_z_y_y_xxyzzzz = cbuffer.data(pk_geom_11_off + 811 * ccomps * dcomps);

            auto g_z_y_y_xxzzzzz = cbuffer.data(pk_geom_11_off + 812 * ccomps * dcomps);

            auto g_z_y_y_xyyyyyy = cbuffer.data(pk_geom_11_off + 813 * ccomps * dcomps);

            auto g_z_y_y_xyyyyyz = cbuffer.data(pk_geom_11_off + 814 * ccomps * dcomps);

            auto g_z_y_y_xyyyyzz = cbuffer.data(pk_geom_11_off + 815 * ccomps * dcomps);

            auto g_z_y_y_xyyyzzz = cbuffer.data(pk_geom_11_off + 816 * ccomps * dcomps);

            auto g_z_y_y_xyyzzzz = cbuffer.data(pk_geom_11_off + 817 * ccomps * dcomps);

            auto g_z_y_y_xyzzzzz = cbuffer.data(pk_geom_11_off + 818 * ccomps * dcomps);

            auto g_z_y_y_xzzzzzz = cbuffer.data(pk_geom_11_off + 819 * ccomps * dcomps);

            auto g_z_y_y_yyyyyyy = cbuffer.data(pk_geom_11_off + 820 * ccomps * dcomps);

            auto g_z_y_y_yyyyyyz = cbuffer.data(pk_geom_11_off + 821 * ccomps * dcomps);

            auto g_z_y_y_yyyyyzz = cbuffer.data(pk_geom_11_off + 822 * ccomps * dcomps);

            auto g_z_y_y_yyyyzzz = cbuffer.data(pk_geom_11_off + 823 * ccomps * dcomps);

            auto g_z_y_y_yyyzzzz = cbuffer.data(pk_geom_11_off + 824 * ccomps * dcomps);

            auto g_z_y_y_yyzzzzz = cbuffer.data(pk_geom_11_off + 825 * ccomps * dcomps);

            auto g_z_y_y_yzzzzzz = cbuffer.data(pk_geom_11_off + 826 * ccomps * dcomps);

            auto g_z_y_y_zzzzzzz = cbuffer.data(pk_geom_11_off + 827 * ccomps * dcomps);

            auto g_z_y_z_xxxxxxx = cbuffer.data(pk_geom_11_off + 828 * ccomps * dcomps);

            auto g_z_y_z_xxxxxxy = cbuffer.data(pk_geom_11_off + 829 * ccomps * dcomps);

            auto g_z_y_z_xxxxxxz = cbuffer.data(pk_geom_11_off + 830 * ccomps * dcomps);

            auto g_z_y_z_xxxxxyy = cbuffer.data(pk_geom_11_off + 831 * ccomps * dcomps);

            auto g_z_y_z_xxxxxyz = cbuffer.data(pk_geom_11_off + 832 * ccomps * dcomps);

            auto g_z_y_z_xxxxxzz = cbuffer.data(pk_geom_11_off + 833 * ccomps * dcomps);

            auto g_z_y_z_xxxxyyy = cbuffer.data(pk_geom_11_off + 834 * ccomps * dcomps);

            auto g_z_y_z_xxxxyyz = cbuffer.data(pk_geom_11_off + 835 * ccomps * dcomps);

            auto g_z_y_z_xxxxyzz = cbuffer.data(pk_geom_11_off + 836 * ccomps * dcomps);

            auto g_z_y_z_xxxxzzz = cbuffer.data(pk_geom_11_off + 837 * ccomps * dcomps);

            auto g_z_y_z_xxxyyyy = cbuffer.data(pk_geom_11_off + 838 * ccomps * dcomps);

            auto g_z_y_z_xxxyyyz = cbuffer.data(pk_geom_11_off + 839 * ccomps * dcomps);

            auto g_z_y_z_xxxyyzz = cbuffer.data(pk_geom_11_off + 840 * ccomps * dcomps);

            auto g_z_y_z_xxxyzzz = cbuffer.data(pk_geom_11_off + 841 * ccomps * dcomps);

            auto g_z_y_z_xxxzzzz = cbuffer.data(pk_geom_11_off + 842 * ccomps * dcomps);

            auto g_z_y_z_xxyyyyy = cbuffer.data(pk_geom_11_off + 843 * ccomps * dcomps);

            auto g_z_y_z_xxyyyyz = cbuffer.data(pk_geom_11_off + 844 * ccomps * dcomps);

            auto g_z_y_z_xxyyyzz = cbuffer.data(pk_geom_11_off + 845 * ccomps * dcomps);

            auto g_z_y_z_xxyyzzz = cbuffer.data(pk_geom_11_off + 846 * ccomps * dcomps);

            auto g_z_y_z_xxyzzzz = cbuffer.data(pk_geom_11_off + 847 * ccomps * dcomps);

            auto g_z_y_z_xxzzzzz = cbuffer.data(pk_geom_11_off + 848 * ccomps * dcomps);

            auto g_z_y_z_xyyyyyy = cbuffer.data(pk_geom_11_off + 849 * ccomps * dcomps);

            auto g_z_y_z_xyyyyyz = cbuffer.data(pk_geom_11_off + 850 * ccomps * dcomps);

            auto g_z_y_z_xyyyyzz = cbuffer.data(pk_geom_11_off + 851 * ccomps * dcomps);

            auto g_z_y_z_xyyyzzz = cbuffer.data(pk_geom_11_off + 852 * ccomps * dcomps);

            auto g_z_y_z_xyyzzzz = cbuffer.data(pk_geom_11_off + 853 * ccomps * dcomps);

            auto g_z_y_z_xyzzzzz = cbuffer.data(pk_geom_11_off + 854 * ccomps * dcomps);

            auto g_z_y_z_xzzzzzz = cbuffer.data(pk_geom_11_off + 855 * ccomps * dcomps);

            auto g_z_y_z_yyyyyyy = cbuffer.data(pk_geom_11_off + 856 * ccomps * dcomps);

            auto g_z_y_z_yyyyyyz = cbuffer.data(pk_geom_11_off + 857 * ccomps * dcomps);

            auto g_z_y_z_yyyyyzz = cbuffer.data(pk_geom_11_off + 858 * ccomps * dcomps);

            auto g_z_y_z_yyyyzzz = cbuffer.data(pk_geom_11_off + 859 * ccomps * dcomps);

            auto g_z_y_z_yyyzzzz = cbuffer.data(pk_geom_11_off + 860 * ccomps * dcomps);

            auto g_z_y_z_yyzzzzz = cbuffer.data(pk_geom_11_off + 861 * ccomps * dcomps);

            auto g_z_y_z_yzzzzzz = cbuffer.data(pk_geom_11_off + 862 * ccomps * dcomps);

            auto g_z_y_z_zzzzzzz = cbuffer.data(pk_geom_11_off + 863 * ccomps * dcomps);

            auto g_z_z_x_xxxxxxx = cbuffer.data(pk_geom_11_off + 864 * ccomps * dcomps);

            auto g_z_z_x_xxxxxxy = cbuffer.data(pk_geom_11_off + 865 * ccomps * dcomps);

            auto g_z_z_x_xxxxxxz = cbuffer.data(pk_geom_11_off + 866 * ccomps * dcomps);

            auto g_z_z_x_xxxxxyy = cbuffer.data(pk_geom_11_off + 867 * ccomps * dcomps);

            auto g_z_z_x_xxxxxyz = cbuffer.data(pk_geom_11_off + 868 * ccomps * dcomps);

            auto g_z_z_x_xxxxxzz = cbuffer.data(pk_geom_11_off + 869 * ccomps * dcomps);

            auto g_z_z_x_xxxxyyy = cbuffer.data(pk_geom_11_off + 870 * ccomps * dcomps);

            auto g_z_z_x_xxxxyyz = cbuffer.data(pk_geom_11_off + 871 * ccomps * dcomps);

            auto g_z_z_x_xxxxyzz = cbuffer.data(pk_geom_11_off + 872 * ccomps * dcomps);

            auto g_z_z_x_xxxxzzz = cbuffer.data(pk_geom_11_off + 873 * ccomps * dcomps);

            auto g_z_z_x_xxxyyyy = cbuffer.data(pk_geom_11_off + 874 * ccomps * dcomps);

            auto g_z_z_x_xxxyyyz = cbuffer.data(pk_geom_11_off + 875 * ccomps * dcomps);

            auto g_z_z_x_xxxyyzz = cbuffer.data(pk_geom_11_off + 876 * ccomps * dcomps);

            auto g_z_z_x_xxxyzzz = cbuffer.data(pk_geom_11_off + 877 * ccomps * dcomps);

            auto g_z_z_x_xxxzzzz = cbuffer.data(pk_geom_11_off + 878 * ccomps * dcomps);

            auto g_z_z_x_xxyyyyy = cbuffer.data(pk_geom_11_off + 879 * ccomps * dcomps);

            auto g_z_z_x_xxyyyyz = cbuffer.data(pk_geom_11_off + 880 * ccomps * dcomps);

            auto g_z_z_x_xxyyyzz = cbuffer.data(pk_geom_11_off + 881 * ccomps * dcomps);

            auto g_z_z_x_xxyyzzz = cbuffer.data(pk_geom_11_off + 882 * ccomps * dcomps);

            auto g_z_z_x_xxyzzzz = cbuffer.data(pk_geom_11_off + 883 * ccomps * dcomps);

            auto g_z_z_x_xxzzzzz = cbuffer.data(pk_geom_11_off + 884 * ccomps * dcomps);

            auto g_z_z_x_xyyyyyy = cbuffer.data(pk_geom_11_off + 885 * ccomps * dcomps);

            auto g_z_z_x_xyyyyyz = cbuffer.data(pk_geom_11_off + 886 * ccomps * dcomps);

            auto g_z_z_x_xyyyyzz = cbuffer.data(pk_geom_11_off + 887 * ccomps * dcomps);

            auto g_z_z_x_xyyyzzz = cbuffer.data(pk_geom_11_off + 888 * ccomps * dcomps);

            auto g_z_z_x_xyyzzzz = cbuffer.data(pk_geom_11_off + 889 * ccomps * dcomps);

            auto g_z_z_x_xyzzzzz = cbuffer.data(pk_geom_11_off + 890 * ccomps * dcomps);

            auto g_z_z_x_xzzzzzz = cbuffer.data(pk_geom_11_off + 891 * ccomps * dcomps);

            auto g_z_z_x_yyyyyyy = cbuffer.data(pk_geom_11_off + 892 * ccomps * dcomps);

            auto g_z_z_x_yyyyyyz = cbuffer.data(pk_geom_11_off + 893 * ccomps * dcomps);

            auto g_z_z_x_yyyyyzz = cbuffer.data(pk_geom_11_off + 894 * ccomps * dcomps);

            auto g_z_z_x_yyyyzzz = cbuffer.data(pk_geom_11_off + 895 * ccomps * dcomps);

            auto g_z_z_x_yyyzzzz = cbuffer.data(pk_geom_11_off + 896 * ccomps * dcomps);

            auto g_z_z_x_yyzzzzz = cbuffer.data(pk_geom_11_off + 897 * ccomps * dcomps);

            auto g_z_z_x_yzzzzzz = cbuffer.data(pk_geom_11_off + 898 * ccomps * dcomps);

            auto g_z_z_x_zzzzzzz = cbuffer.data(pk_geom_11_off + 899 * ccomps * dcomps);

            auto g_z_z_y_xxxxxxx = cbuffer.data(pk_geom_11_off + 900 * ccomps * dcomps);

            auto g_z_z_y_xxxxxxy = cbuffer.data(pk_geom_11_off + 901 * ccomps * dcomps);

            auto g_z_z_y_xxxxxxz = cbuffer.data(pk_geom_11_off + 902 * ccomps * dcomps);

            auto g_z_z_y_xxxxxyy = cbuffer.data(pk_geom_11_off + 903 * ccomps * dcomps);

            auto g_z_z_y_xxxxxyz = cbuffer.data(pk_geom_11_off + 904 * ccomps * dcomps);

            auto g_z_z_y_xxxxxzz = cbuffer.data(pk_geom_11_off + 905 * ccomps * dcomps);

            auto g_z_z_y_xxxxyyy = cbuffer.data(pk_geom_11_off + 906 * ccomps * dcomps);

            auto g_z_z_y_xxxxyyz = cbuffer.data(pk_geom_11_off + 907 * ccomps * dcomps);

            auto g_z_z_y_xxxxyzz = cbuffer.data(pk_geom_11_off + 908 * ccomps * dcomps);

            auto g_z_z_y_xxxxzzz = cbuffer.data(pk_geom_11_off + 909 * ccomps * dcomps);

            auto g_z_z_y_xxxyyyy = cbuffer.data(pk_geom_11_off + 910 * ccomps * dcomps);

            auto g_z_z_y_xxxyyyz = cbuffer.data(pk_geom_11_off + 911 * ccomps * dcomps);

            auto g_z_z_y_xxxyyzz = cbuffer.data(pk_geom_11_off + 912 * ccomps * dcomps);

            auto g_z_z_y_xxxyzzz = cbuffer.data(pk_geom_11_off + 913 * ccomps * dcomps);

            auto g_z_z_y_xxxzzzz = cbuffer.data(pk_geom_11_off + 914 * ccomps * dcomps);

            auto g_z_z_y_xxyyyyy = cbuffer.data(pk_geom_11_off + 915 * ccomps * dcomps);

            auto g_z_z_y_xxyyyyz = cbuffer.data(pk_geom_11_off + 916 * ccomps * dcomps);

            auto g_z_z_y_xxyyyzz = cbuffer.data(pk_geom_11_off + 917 * ccomps * dcomps);

            auto g_z_z_y_xxyyzzz = cbuffer.data(pk_geom_11_off + 918 * ccomps * dcomps);

            auto g_z_z_y_xxyzzzz = cbuffer.data(pk_geom_11_off + 919 * ccomps * dcomps);

            auto g_z_z_y_xxzzzzz = cbuffer.data(pk_geom_11_off + 920 * ccomps * dcomps);

            auto g_z_z_y_xyyyyyy = cbuffer.data(pk_geom_11_off + 921 * ccomps * dcomps);

            auto g_z_z_y_xyyyyyz = cbuffer.data(pk_geom_11_off + 922 * ccomps * dcomps);

            auto g_z_z_y_xyyyyzz = cbuffer.data(pk_geom_11_off + 923 * ccomps * dcomps);

            auto g_z_z_y_xyyyzzz = cbuffer.data(pk_geom_11_off + 924 * ccomps * dcomps);

            auto g_z_z_y_xyyzzzz = cbuffer.data(pk_geom_11_off + 925 * ccomps * dcomps);

            auto g_z_z_y_xyzzzzz = cbuffer.data(pk_geom_11_off + 926 * ccomps * dcomps);

            auto g_z_z_y_xzzzzzz = cbuffer.data(pk_geom_11_off + 927 * ccomps * dcomps);

            auto g_z_z_y_yyyyyyy = cbuffer.data(pk_geom_11_off + 928 * ccomps * dcomps);

            auto g_z_z_y_yyyyyyz = cbuffer.data(pk_geom_11_off + 929 * ccomps * dcomps);

            auto g_z_z_y_yyyyyzz = cbuffer.data(pk_geom_11_off + 930 * ccomps * dcomps);

            auto g_z_z_y_yyyyzzz = cbuffer.data(pk_geom_11_off + 931 * ccomps * dcomps);

            auto g_z_z_y_yyyzzzz = cbuffer.data(pk_geom_11_off + 932 * ccomps * dcomps);

            auto g_z_z_y_yyzzzzz = cbuffer.data(pk_geom_11_off + 933 * ccomps * dcomps);

            auto g_z_z_y_yzzzzzz = cbuffer.data(pk_geom_11_off + 934 * ccomps * dcomps);

            auto g_z_z_y_zzzzzzz = cbuffer.data(pk_geom_11_off + 935 * ccomps * dcomps);

            auto g_z_z_z_xxxxxxx = cbuffer.data(pk_geom_11_off + 936 * ccomps * dcomps);

            auto g_z_z_z_xxxxxxy = cbuffer.data(pk_geom_11_off + 937 * ccomps * dcomps);

            auto g_z_z_z_xxxxxxz = cbuffer.data(pk_geom_11_off + 938 * ccomps * dcomps);

            auto g_z_z_z_xxxxxyy = cbuffer.data(pk_geom_11_off + 939 * ccomps * dcomps);

            auto g_z_z_z_xxxxxyz = cbuffer.data(pk_geom_11_off + 940 * ccomps * dcomps);

            auto g_z_z_z_xxxxxzz = cbuffer.data(pk_geom_11_off + 941 * ccomps * dcomps);

            auto g_z_z_z_xxxxyyy = cbuffer.data(pk_geom_11_off + 942 * ccomps * dcomps);

            auto g_z_z_z_xxxxyyz = cbuffer.data(pk_geom_11_off + 943 * ccomps * dcomps);

            auto g_z_z_z_xxxxyzz = cbuffer.data(pk_geom_11_off + 944 * ccomps * dcomps);

            auto g_z_z_z_xxxxzzz = cbuffer.data(pk_geom_11_off + 945 * ccomps * dcomps);

            auto g_z_z_z_xxxyyyy = cbuffer.data(pk_geom_11_off + 946 * ccomps * dcomps);

            auto g_z_z_z_xxxyyyz = cbuffer.data(pk_geom_11_off + 947 * ccomps * dcomps);

            auto g_z_z_z_xxxyyzz = cbuffer.data(pk_geom_11_off + 948 * ccomps * dcomps);

            auto g_z_z_z_xxxyzzz = cbuffer.data(pk_geom_11_off + 949 * ccomps * dcomps);

            auto g_z_z_z_xxxzzzz = cbuffer.data(pk_geom_11_off + 950 * ccomps * dcomps);

            auto g_z_z_z_xxyyyyy = cbuffer.data(pk_geom_11_off + 951 * ccomps * dcomps);

            auto g_z_z_z_xxyyyyz = cbuffer.data(pk_geom_11_off + 952 * ccomps * dcomps);

            auto g_z_z_z_xxyyyzz = cbuffer.data(pk_geom_11_off + 953 * ccomps * dcomps);

            auto g_z_z_z_xxyyzzz = cbuffer.data(pk_geom_11_off + 954 * ccomps * dcomps);

            auto g_z_z_z_xxyzzzz = cbuffer.data(pk_geom_11_off + 955 * ccomps * dcomps);

            auto g_z_z_z_xxzzzzz = cbuffer.data(pk_geom_11_off + 956 * ccomps * dcomps);

            auto g_z_z_z_xyyyyyy = cbuffer.data(pk_geom_11_off + 957 * ccomps * dcomps);

            auto g_z_z_z_xyyyyyz = cbuffer.data(pk_geom_11_off + 958 * ccomps * dcomps);

            auto g_z_z_z_xyyyyzz = cbuffer.data(pk_geom_11_off + 959 * ccomps * dcomps);

            auto g_z_z_z_xyyyzzz = cbuffer.data(pk_geom_11_off + 960 * ccomps * dcomps);

            auto g_z_z_z_xyyzzzz = cbuffer.data(pk_geom_11_off + 961 * ccomps * dcomps);

            auto g_z_z_z_xyzzzzz = cbuffer.data(pk_geom_11_off + 962 * ccomps * dcomps);

            auto g_z_z_z_xzzzzzz = cbuffer.data(pk_geom_11_off + 963 * ccomps * dcomps);

            auto g_z_z_z_yyyyyyy = cbuffer.data(pk_geom_11_off + 964 * ccomps * dcomps);

            auto g_z_z_z_yyyyyyz = cbuffer.data(pk_geom_11_off + 965 * ccomps * dcomps);

            auto g_z_z_z_yyyyyzz = cbuffer.data(pk_geom_11_off + 966 * ccomps * dcomps);

            auto g_z_z_z_yyyyzzz = cbuffer.data(pk_geom_11_off + 967 * ccomps * dcomps);

            auto g_z_z_z_yyyzzzz = cbuffer.data(pk_geom_11_off + 968 * ccomps * dcomps);

            auto g_z_z_z_yyzzzzz = cbuffer.data(pk_geom_11_off + 969 * ccomps * dcomps);

            auto g_z_z_z_yzzzzzz = cbuffer.data(pk_geom_11_off + 970 * ccomps * dcomps);

            auto g_z_z_z_zzzzzzz = cbuffer.data(pk_geom_11_off + 971 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_dixx

            const auto di_geom_11_off = idx_geom_11_dixx + i * dcomps + j;

            /// Set up 0-28 components of targeted buffer : cbuffer.data(

            auto g_x_x_xx_xxxxxx = cbuffer.data(di_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_xxxxxy = cbuffer.data(di_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_xxxxxz = cbuffer.data(di_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xx_xxxxyy = cbuffer.data(di_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xx_xxxxyz = cbuffer.data(di_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xx_xxxxzz = cbuffer.data(di_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xx_xxxyyy = cbuffer.data(di_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xx_xxxyyz = cbuffer.data(di_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xx_xxxyzz = cbuffer.data(di_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xx_xxxzzz = cbuffer.data(di_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xx_xxyyyy = cbuffer.data(di_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xx_xxyyyz = cbuffer.data(di_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xx_xxyyzz = cbuffer.data(di_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xx_xxyzzz = cbuffer.data(di_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xx_xxzzzz = cbuffer.data(di_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xx_xyyyyy = cbuffer.data(di_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xx_xyyyyz = cbuffer.data(di_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xx_xyyyzz = cbuffer.data(di_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xx_xyyzzz = cbuffer.data(di_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xx_xyzzzz = cbuffer.data(di_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xx_xzzzzz = cbuffer.data(di_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_xx_yyyyyy = cbuffer.data(di_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xx_yyyyyz = cbuffer.data(di_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xx_yyyyzz = cbuffer.data(di_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xx_yyyzzz = cbuffer.data(di_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xx_yyzzzz = cbuffer.data(di_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xx_yzzzzz = cbuffer.data(di_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xx_zzzzzz = cbuffer.data(di_geom_11_off + 27 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xxxxxx, g_0_x_x_xxxxxy, g_0_x_x_xxxxxz, g_0_x_x_xxxxyy, g_0_x_x_xxxxyz, g_0_x_x_xxxxzz, g_0_x_x_xxxyyy, g_0_x_x_xxxyyz, g_0_x_x_xxxyzz, g_0_x_x_xxxzzz, g_0_x_x_xxyyyy, g_0_x_x_xxyyyz, g_0_x_x_xxyyzz, g_0_x_x_xxyzzz, g_0_x_x_xxzzzz, g_0_x_x_xyyyyy, g_0_x_x_xyyyyz, g_0_x_x_xyyyzz, g_0_x_x_xyyzzz, g_0_x_x_xyzzzz, g_0_x_x_xzzzzz, g_0_x_x_yyyyyy, g_0_x_x_yyyyyz, g_0_x_x_yyyyzz, g_0_x_x_yyyzzz, g_0_x_x_yyzzzz, g_0_x_x_yzzzzz, g_0_x_x_zzzzzz, g_x_0_x_xxxxxx, g_x_0_x_xxxxxy, g_x_0_x_xxxxxz, g_x_0_x_xxxxyy, g_x_0_x_xxxxyz, g_x_0_x_xxxxzz, g_x_0_x_xxxyyy, g_x_0_x_xxxyyz, g_x_0_x_xxxyzz, g_x_0_x_xxxzzz, g_x_0_x_xxyyyy, g_x_0_x_xxyyyz, g_x_0_x_xxyyzz, g_x_0_x_xxyzzz, g_x_0_x_xxzzzz, g_x_0_x_xyyyyy, g_x_0_x_xyyyyz, g_x_0_x_xyyyzz, g_x_0_x_xyyzzz, g_x_0_x_xyzzzz, g_x_0_x_xzzzzz, g_x_0_x_yyyyyy, g_x_0_x_yyyyyz, g_x_0_x_yyyyzz, g_x_0_x_yyyzzz, g_x_0_x_yyzzzz, g_x_0_x_yzzzzz, g_x_0_x_zzzzzz, g_x_x_x_xxxxxx, g_x_x_x_xxxxxxx, g_x_x_x_xxxxxxy, g_x_x_x_xxxxxxz, g_x_x_x_xxxxxy, g_x_x_x_xxxxxyy, g_x_x_x_xxxxxyz, g_x_x_x_xxxxxz, g_x_x_x_xxxxxzz, g_x_x_x_xxxxyy, g_x_x_x_xxxxyyy, g_x_x_x_xxxxyyz, g_x_x_x_xxxxyz, g_x_x_x_xxxxyzz, g_x_x_x_xxxxzz, g_x_x_x_xxxxzzz, g_x_x_x_xxxyyy, g_x_x_x_xxxyyyy, g_x_x_x_xxxyyyz, g_x_x_x_xxxyyz, g_x_x_x_xxxyyzz, g_x_x_x_xxxyzz, g_x_x_x_xxxyzzz, g_x_x_x_xxxzzz, g_x_x_x_xxxzzzz, g_x_x_x_xxyyyy, g_x_x_x_xxyyyyy, g_x_x_x_xxyyyyz, g_x_x_x_xxyyyz, g_x_x_x_xxyyyzz, g_x_x_x_xxyyzz, g_x_x_x_xxyyzzz, g_x_x_x_xxyzzz, g_x_x_x_xxyzzzz, g_x_x_x_xxzzzz, g_x_x_x_xxzzzzz, g_x_x_x_xyyyyy, g_x_x_x_xyyyyyy, g_x_x_x_xyyyyyz, g_x_x_x_xyyyyz, g_x_x_x_xyyyyzz, g_x_x_x_xyyyzz, g_x_x_x_xyyyzzz, g_x_x_x_xyyzzz, g_x_x_x_xyyzzzz, g_x_x_x_xyzzzz, g_x_x_x_xyzzzzz, g_x_x_x_xzzzzz, g_x_x_x_xzzzzzz, g_x_x_x_yyyyyy, g_x_x_x_yyyyyz, g_x_x_x_yyyyzz, g_x_x_x_yyyzzz, g_x_x_x_yyzzzz, g_x_x_x_yzzzzz, g_x_x_x_zzzzzz, g_x_x_xx_xxxxxx, g_x_x_xx_xxxxxy, g_x_x_xx_xxxxxz, g_x_x_xx_xxxxyy, g_x_x_xx_xxxxyz, g_x_x_xx_xxxxzz, g_x_x_xx_xxxyyy, g_x_x_xx_xxxyyz, g_x_x_xx_xxxyzz, g_x_x_xx_xxxzzz, g_x_x_xx_xxyyyy, g_x_x_xx_xxyyyz, g_x_x_xx_xxyyzz, g_x_x_xx_xxyzzz, g_x_x_xx_xxzzzz, g_x_x_xx_xyyyyy, g_x_x_xx_xyyyyz, g_x_x_xx_xyyyzz, g_x_x_xx_xyyzzz, g_x_x_xx_xyzzzz, g_x_x_xx_xzzzzz, g_x_x_xx_yyyyyy, g_x_x_xx_yyyyyz, g_x_x_xx_yyyyzz, g_x_x_xx_yyyzzz, g_x_x_xx_yyzzzz, g_x_x_xx_yzzzzz, g_x_x_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xx_xxxxxx[k] = -g_0_x_x_xxxxxx[k] + g_x_0_x_xxxxxx[k] - g_x_x_x_xxxxxx[k] * ab_x + g_x_x_x_xxxxxxx[k];

                g_x_x_xx_xxxxxy[k] = -g_0_x_x_xxxxxy[k] + g_x_0_x_xxxxxy[k] - g_x_x_x_xxxxxy[k] * ab_x + g_x_x_x_xxxxxxy[k];

                g_x_x_xx_xxxxxz[k] = -g_0_x_x_xxxxxz[k] + g_x_0_x_xxxxxz[k] - g_x_x_x_xxxxxz[k] * ab_x + g_x_x_x_xxxxxxz[k];

                g_x_x_xx_xxxxyy[k] = -g_0_x_x_xxxxyy[k] + g_x_0_x_xxxxyy[k] - g_x_x_x_xxxxyy[k] * ab_x + g_x_x_x_xxxxxyy[k];

                g_x_x_xx_xxxxyz[k] = -g_0_x_x_xxxxyz[k] + g_x_0_x_xxxxyz[k] - g_x_x_x_xxxxyz[k] * ab_x + g_x_x_x_xxxxxyz[k];

                g_x_x_xx_xxxxzz[k] = -g_0_x_x_xxxxzz[k] + g_x_0_x_xxxxzz[k] - g_x_x_x_xxxxzz[k] * ab_x + g_x_x_x_xxxxxzz[k];

                g_x_x_xx_xxxyyy[k] = -g_0_x_x_xxxyyy[k] + g_x_0_x_xxxyyy[k] - g_x_x_x_xxxyyy[k] * ab_x + g_x_x_x_xxxxyyy[k];

                g_x_x_xx_xxxyyz[k] = -g_0_x_x_xxxyyz[k] + g_x_0_x_xxxyyz[k] - g_x_x_x_xxxyyz[k] * ab_x + g_x_x_x_xxxxyyz[k];

                g_x_x_xx_xxxyzz[k] = -g_0_x_x_xxxyzz[k] + g_x_0_x_xxxyzz[k] - g_x_x_x_xxxyzz[k] * ab_x + g_x_x_x_xxxxyzz[k];

                g_x_x_xx_xxxzzz[k] = -g_0_x_x_xxxzzz[k] + g_x_0_x_xxxzzz[k] - g_x_x_x_xxxzzz[k] * ab_x + g_x_x_x_xxxxzzz[k];

                g_x_x_xx_xxyyyy[k] = -g_0_x_x_xxyyyy[k] + g_x_0_x_xxyyyy[k] - g_x_x_x_xxyyyy[k] * ab_x + g_x_x_x_xxxyyyy[k];

                g_x_x_xx_xxyyyz[k] = -g_0_x_x_xxyyyz[k] + g_x_0_x_xxyyyz[k] - g_x_x_x_xxyyyz[k] * ab_x + g_x_x_x_xxxyyyz[k];

                g_x_x_xx_xxyyzz[k] = -g_0_x_x_xxyyzz[k] + g_x_0_x_xxyyzz[k] - g_x_x_x_xxyyzz[k] * ab_x + g_x_x_x_xxxyyzz[k];

                g_x_x_xx_xxyzzz[k] = -g_0_x_x_xxyzzz[k] + g_x_0_x_xxyzzz[k] - g_x_x_x_xxyzzz[k] * ab_x + g_x_x_x_xxxyzzz[k];

                g_x_x_xx_xxzzzz[k] = -g_0_x_x_xxzzzz[k] + g_x_0_x_xxzzzz[k] - g_x_x_x_xxzzzz[k] * ab_x + g_x_x_x_xxxzzzz[k];

                g_x_x_xx_xyyyyy[k] = -g_0_x_x_xyyyyy[k] + g_x_0_x_xyyyyy[k] - g_x_x_x_xyyyyy[k] * ab_x + g_x_x_x_xxyyyyy[k];

                g_x_x_xx_xyyyyz[k] = -g_0_x_x_xyyyyz[k] + g_x_0_x_xyyyyz[k] - g_x_x_x_xyyyyz[k] * ab_x + g_x_x_x_xxyyyyz[k];

                g_x_x_xx_xyyyzz[k] = -g_0_x_x_xyyyzz[k] + g_x_0_x_xyyyzz[k] - g_x_x_x_xyyyzz[k] * ab_x + g_x_x_x_xxyyyzz[k];

                g_x_x_xx_xyyzzz[k] = -g_0_x_x_xyyzzz[k] + g_x_0_x_xyyzzz[k] - g_x_x_x_xyyzzz[k] * ab_x + g_x_x_x_xxyyzzz[k];

                g_x_x_xx_xyzzzz[k] = -g_0_x_x_xyzzzz[k] + g_x_0_x_xyzzzz[k] - g_x_x_x_xyzzzz[k] * ab_x + g_x_x_x_xxyzzzz[k];

                g_x_x_xx_xzzzzz[k] = -g_0_x_x_xzzzzz[k] + g_x_0_x_xzzzzz[k] - g_x_x_x_xzzzzz[k] * ab_x + g_x_x_x_xxzzzzz[k];

                g_x_x_xx_yyyyyy[k] = -g_0_x_x_yyyyyy[k] + g_x_0_x_yyyyyy[k] - g_x_x_x_yyyyyy[k] * ab_x + g_x_x_x_xyyyyyy[k];

                g_x_x_xx_yyyyyz[k] = -g_0_x_x_yyyyyz[k] + g_x_0_x_yyyyyz[k] - g_x_x_x_yyyyyz[k] * ab_x + g_x_x_x_xyyyyyz[k];

                g_x_x_xx_yyyyzz[k] = -g_0_x_x_yyyyzz[k] + g_x_0_x_yyyyzz[k] - g_x_x_x_yyyyzz[k] * ab_x + g_x_x_x_xyyyyzz[k];

                g_x_x_xx_yyyzzz[k] = -g_0_x_x_yyyzzz[k] + g_x_0_x_yyyzzz[k] - g_x_x_x_yyyzzz[k] * ab_x + g_x_x_x_xyyyzzz[k];

                g_x_x_xx_yyzzzz[k] = -g_0_x_x_yyzzzz[k] + g_x_0_x_yyzzzz[k] - g_x_x_x_yyzzzz[k] * ab_x + g_x_x_x_xyyzzzz[k];

                g_x_x_xx_yzzzzz[k] = -g_0_x_x_yzzzzz[k] + g_x_0_x_yzzzzz[k] - g_x_x_x_yzzzzz[k] * ab_x + g_x_x_x_xyzzzzz[k];

                g_x_x_xx_zzzzzz[k] = -g_0_x_x_zzzzzz[k] + g_x_0_x_zzzzzz[k] - g_x_x_x_zzzzzz[k] * ab_x + g_x_x_x_xzzzzzz[k];
            }

            /// Set up 28-56 components of targeted buffer : cbuffer.data(

            auto g_x_x_xy_xxxxxx = cbuffer.data(di_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xy_xxxxxy = cbuffer.data(di_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_xy_xxxxxz = cbuffer.data(di_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xy_xxxxyy = cbuffer.data(di_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xy_xxxxyz = cbuffer.data(di_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xy_xxxxzz = cbuffer.data(di_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xy_xxxyyy = cbuffer.data(di_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xy_xxxyyz = cbuffer.data(di_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_xy_xxxyzz = cbuffer.data(di_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xy_xxxzzz = cbuffer.data(di_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xy_xxyyyy = cbuffer.data(di_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xy_xxyyyz = cbuffer.data(di_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_xy_xxyyzz = cbuffer.data(di_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xy_xxyzzz = cbuffer.data(di_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_xy_xxzzzz = cbuffer.data(di_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xy_xyyyyy = cbuffer.data(di_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xy_xyyyyz = cbuffer.data(di_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_xy_xyyyzz = cbuffer.data(di_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_xy_xyyzzz = cbuffer.data(di_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_xy_xyzzzz = cbuffer.data(di_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_xy_xzzzzz = cbuffer.data(di_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_xy_yyyyyy = cbuffer.data(di_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_xy_yyyyyz = cbuffer.data(di_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_xy_yyyyzz = cbuffer.data(di_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_xy_yyyzzz = cbuffer.data(di_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_xy_yyzzzz = cbuffer.data(di_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_xy_yzzzzz = cbuffer.data(di_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_xy_zzzzzz = cbuffer.data(di_geom_11_off + 55 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_xxxxxx, g_x_x_x_xxxxxxy, g_x_x_x_xxxxxy, g_x_x_x_xxxxxyy, g_x_x_x_xxxxxyz, g_x_x_x_xxxxxz, g_x_x_x_xxxxyy, g_x_x_x_xxxxyyy, g_x_x_x_xxxxyyz, g_x_x_x_xxxxyz, g_x_x_x_xxxxyzz, g_x_x_x_xxxxzz, g_x_x_x_xxxyyy, g_x_x_x_xxxyyyy, g_x_x_x_xxxyyyz, g_x_x_x_xxxyyz, g_x_x_x_xxxyyzz, g_x_x_x_xxxyzz, g_x_x_x_xxxyzzz, g_x_x_x_xxxzzz, g_x_x_x_xxyyyy, g_x_x_x_xxyyyyy, g_x_x_x_xxyyyyz, g_x_x_x_xxyyyz, g_x_x_x_xxyyyzz, g_x_x_x_xxyyzz, g_x_x_x_xxyyzzz, g_x_x_x_xxyzzz, g_x_x_x_xxyzzzz, g_x_x_x_xxzzzz, g_x_x_x_xyyyyy, g_x_x_x_xyyyyyy, g_x_x_x_xyyyyyz, g_x_x_x_xyyyyz, g_x_x_x_xyyyyzz, g_x_x_x_xyyyzz, g_x_x_x_xyyyzzz, g_x_x_x_xyyzzz, g_x_x_x_xyyzzzz, g_x_x_x_xyzzzz, g_x_x_x_xyzzzzz, g_x_x_x_xzzzzz, g_x_x_x_yyyyyy, g_x_x_x_yyyyyyy, g_x_x_x_yyyyyyz, g_x_x_x_yyyyyz, g_x_x_x_yyyyyzz, g_x_x_x_yyyyzz, g_x_x_x_yyyyzzz, g_x_x_x_yyyzzz, g_x_x_x_yyyzzzz, g_x_x_x_yyzzzz, g_x_x_x_yyzzzzz, g_x_x_x_yzzzzz, g_x_x_x_yzzzzzz, g_x_x_x_zzzzzz, g_x_x_xy_xxxxxx, g_x_x_xy_xxxxxy, g_x_x_xy_xxxxxz, g_x_x_xy_xxxxyy, g_x_x_xy_xxxxyz, g_x_x_xy_xxxxzz, g_x_x_xy_xxxyyy, g_x_x_xy_xxxyyz, g_x_x_xy_xxxyzz, g_x_x_xy_xxxzzz, g_x_x_xy_xxyyyy, g_x_x_xy_xxyyyz, g_x_x_xy_xxyyzz, g_x_x_xy_xxyzzz, g_x_x_xy_xxzzzz, g_x_x_xy_xyyyyy, g_x_x_xy_xyyyyz, g_x_x_xy_xyyyzz, g_x_x_xy_xyyzzz, g_x_x_xy_xyzzzz, g_x_x_xy_xzzzzz, g_x_x_xy_yyyyyy, g_x_x_xy_yyyyyz, g_x_x_xy_yyyyzz, g_x_x_xy_yyyzzz, g_x_x_xy_yyzzzz, g_x_x_xy_yzzzzz, g_x_x_xy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xy_xxxxxx[k] = -g_x_x_x_xxxxxx[k] * ab_y + g_x_x_x_xxxxxxy[k];

                g_x_x_xy_xxxxxy[k] = -g_x_x_x_xxxxxy[k] * ab_y + g_x_x_x_xxxxxyy[k];

                g_x_x_xy_xxxxxz[k] = -g_x_x_x_xxxxxz[k] * ab_y + g_x_x_x_xxxxxyz[k];

                g_x_x_xy_xxxxyy[k] = -g_x_x_x_xxxxyy[k] * ab_y + g_x_x_x_xxxxyyy[k];

                g_x_x_xy_xxxxyz[k] = -g_x_x_x_xxxxyz[k] * ab_y + g_x_x_x_xxxxyyz[k];

                g_x_x_xy_xxxxzz[k] = -g_x_x_x_xxxxzz[k] * ab_y + g_x_x_x_xxxxyzz[k];

                g_x_x_xy_xxxyyy[k] = -g_x_x_x_xxxyyy[k] * ab_y + g_x_x_x_xxxyyyy[k];

                g_x_x_xy_xxxyyz[k] = -g_x_x_x_xxxyyz[k] * ab_y + g_x_x_x_xxxyyyz[k];

                g_x_x_xy_xxxyzz[k] = -g_x_x_x_xxxyzz[k] * ab_y + g_x_x_x_xxxyyzz[k];

                g_x_x_xy_xxxzzz[k] = -g_x_x_x_xxxzzz[k] * ab_y + g_x_x_x_xxxyzzz[k];

                g_x_x_xy_xxyyyy[k] = -g_x_x_x_xxyyyy[k] * ab_y + g_x_x_x_xxyyyyy[k];

                g_x_x_xy_xxyyyz[k] = -g_x_x_x_xxyyyz[k] * ab_y + g_x_x_x_xxyyyyz[k];

                g_x_x_xy_xxyyzz[k] = -g_x_x_x_xxyyzz[k] * ab_y + g_x_x_x_xxyyyzz[k];

                g_x_x_xy_xxyzzz[k] = -g_x_x_x_xxyzzz[k] * ab_y + g_x_x_x_xxyyzzz[k];

                g_x_x_xy_xxzzzz[k] = -g_x_x_x_xxzzzz[k] * ab_y + g_x_x_x_xxyzzzz[k];

                g_x_x_xy_xyyyyy[k] = -g_x_x_x_xyyyyy[k] * ab_y + g_x_x_x_xyyyyyy[k];

                g_x_x_xy_xyyyyz[k] = -g_x_x_x_xyyyyz[k] * ab_y + g_x_x_x_xyyyyyz[k];

                g_x_x_xy_xyyyzz[k] = -g_x_x_x_xyyyzz[k] * ab_y + g_x_x_x_xyyyyzz[k];

                g_x_x_xy_xyyzzz[k] = -g_x_x_x_xyyzzz[k] * ab_y + g_x_x_x_xyyyzzz[k];

                g_x_x_xy_xyzzzz[k] = -g_x_x_x_xyzzzz[k] * ab_y + g_x_x_x_xyyzzzz[k];

                g_x_x_xy_xzzzzz[k] = -g_x_x_x_xzzzzz[k] * ab_y + g_x_x_x_xyzzzzz[k];

                g_x_x_xy_yyyyyy[k] = -g_x_x_x_yyyyyy[k] * ab_y + g_x_x_x_yyyyyyy[k];

                g_x_x_xy_yyyyyz[k] = -g_x_x_x_yyyyyz[k] * ab_y + g_x_x_x_yyyyyyz[k];

                g_x_x_xy_yyyyzz[k] = -g_x_x_x_yyyyzz[k] * ab_y + g_x_x_x_yyyyyzz[k];

                g_x_x_xy_yyyzzz[k] = -g_x_x_x_yyyzzz[k] * ab_y + g_x_x_x_yyyyzzz[k];

                g_x_x_xy_yyzzzz[k] = -g_x_x_x_yyzzzz[k] * ab_y + g_x_x_x_yyyzzzz[k];

                g_x_x_xy_yzzzzz[k] = -g_x_x_x_yzzzzz[k] * ab_y + g_x_x_x_yyzzzzz[k];

                g_x_x_xy_zzzzzz[k] = -g_x_x_x_zzzzzz[k] * ab_y + g_x_x_x_yzzzzzz[k];
            }

            /// Set up 56-84 components of targeted buffer : cbuffer.data(

            auto g_x_x_xz_xxxxxx = cbuffer.data(di_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_xz_xxxxxy = cbuffer.data(di_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_xz_xxxxxz = cbuffer.data(di_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_xz_xxxxyy = cbuffer.data(di_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_xz_xxxxyz = cbuffer.data(di_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_xz_xxxxzz = cbuffer.data(di_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_xz_xxxyyy = cbuffer.data(di_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_x_xz_xxxyyz = cbuffer.data(di_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_xz_xxxyzz = cbuffer.data(di_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_xz_xxxzzz = cbuffer.data(di_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_xz_xxyyyy = cbuffer.data(di_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_xz_xxyyyz = cbuffer.data(di_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_xz_xxyyzz = cbuffer.data(di_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_xz_xxyzzz = cbuffer.data(di_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_xz_xxzzzz = cbuffer.data(di_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_xz_xyyyyy = cbuffer.data(di_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_xz_xyyyyz = cbuffer.data(di_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_xz_xyyyzz = cbuffer.data(di_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_xz_xyyzzz = cbuffer.data(di_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_xz_xyzzzz = cbuffer.data(di_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_xz_xzzzzz = cbuffer.data(di_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_xz_yyyyyy = cbuffer.data(di_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_xz_yyyyyz = cbuffer.data(di_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_xz_yyyyzz = cbuffer.data(di_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_xz_yyyzzz = cbuffer.data(di_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_xz_yyzzzz = cbuffer.data(di_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_xz_yzzzzz = cbuffer.data(di_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_xz_zzzzzz = cbuffer.data(di_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_xxxxxx, g_x_x_x_xxxxxxz, g_x_x_x_xxxxxy, g_x_x_x_xxxxxyz, g_x_x_x_xxxxxz, g_x_x_x_xxxxxzz, g_x_x_x_xxxxyy, g_x_x_x_xxxxyyz, g_x_x_x_xxxxyz, g_x_x_x_xxxxyzz, g_x_x_x_xxxxzz, g_x_x_x_xxxxzzz, g_x_x_x_xxxyyy, g_x_x_x_xxxyyyz, g_x_x_x_xxxyyz, g_x_x_x_xxxyyzz, g_x_x_x_xxxyzz, g_x_x_x_xxxyzzz, g_x_x_x_xxxzzz, g_x_x_x_xxxzzzz, g_x_x_x_xxyyyy, g_x_x_x_xxyyyyz, g_x_x_x_xxyyyz, g_x_x_x_xxyyyzz, g_x_x_x_xxyyzz, g_x_x_x_xxyyzzz, g_x_x_x_xxyzzz, g_x_x_x_xxyzzzz, g_x_x_x_xxzzzz, g_x_x_x_xxzzzzz, g_x_x_x_xyyyyy, g_x_x_x_xyyyyyz, g_x_x_x_xyyyyz, g_x_x_x_xyyyyzz, g_x_x_x_xyyyzz, g_x_x_x_xyyyzzz, g_x_x_x_xyyzzz, g_x_x_x_xyyzzzz, g_x_x_x_xyzzzz, g_x_x_x_xyzzzzz, g_x_x_x_xzzzzz, g_x_x_x_xzzzzzz, g_x_x_x_yyyyyy, g_x_x_x_yyyyyyz, g_x_x_x_yyyyyz, g_x_x_x_yyyyyzz, g_x_x_x_yyyyzz, g_x_x_x_yyyyzzz, g_x_x_x_yyyzzz, g_x_x_x_yyyzzzz, g_x_x_x_yyzzzz, g_x_x_x_yyzzzzz, g_x_x_x_yzzzzz, g_x_x_x_yzzzzzz, g_x_x_x_zzzzzz, g_x_x_x_zzzzzzz, g_x_x_xz_xxxxxx, g_x_x_xz_xxxxxy, g_x_x_xz_xxxxxz, g_x_x_xz_xxxxyy, g_x_x_xz_xxxxyz, g_x_x_xz_xxxxzz, g_x_x_xz_xxxyyy, g_x_x_xz_xxxyyz, g_x_x_xz_xxxyzz, g_x_x_xz_xxxzzz, g_x_x_xz_xxyyyy, g_x_x_xz_xxyyyz, g_x_x_xz_xxyyzz, g_x_x_xz_xxyzzz, g_x_x_xz_xxzzzz, g_x_x_xz_xyyyyy, g_x_x_xz_xyyyyz, g_x_x_xz_xyyyzz, g_x_x_xz_xyyzzz, g_x_x_xz_xyzzzz, g_x_x_xz_xzzzzz, g_x_x_xz_yyyyyy, g_x_x_xz_yyyyyz, g_x_x_xz_yyyyzz, g_x_x_xz_yyyzzz, g_x_x_xz_yyzzzz, g_x_x_xz_yzzzzz, g_x_x_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xz_xxxxxx[k] = -g_x_x_x_xxxxxx[k] * ab_z + g_x_x_x_xxxxxxz[k];

                g_x_x_xz_xxxxxy[k] = -g_x_x_x_xxxxxy[k] * ab_z + g_x_x_x_xxxxxyz[k];

                g_x_x_xz_xxxxxz[k] = -g_x_x_x_xxxxxz[k] * ab_z + g_x_x_x_xxxxxzz[k];

                g_x_x_xz_xxxxyy[k] = -g_x_x_x_xxxxyy[k] * ab_z + g_x_x_x_xxxxyyz[k];

                g_x_x_xz_xxxxyz[k] = -g_x_x_x_xxxxyz[k] * ab_z + g_x_x_x_xxxxyzz[k];

                g_x_x_xz_xxxxzz[k] = -g_x_x_x_xxxxzz[k] * ab_z + g_x_x_x_xxxxzzz[k];

                g_x_x_xz_xxxyyy[k] = -g_x_x_x_xxxyyy[k] * ab_z + g_x_x_x_xxxyyyz[k];

                g_x_x_xz_xxxyyz[k] = -g_x_x_x_xxxyyz[k] * ab_z + g_x_x_x_xxxyyzz[k];

                g_x_x_xz_xxxyzz[k] = -g_x_x_x_xxxyzz[k] * ab_z + g_x_x_x_xxxyzzz[k];

                g_x_x_xz_xxxzzz[k] = -g_x_x_x_xxxzzz[k] * ab_z + g_x_x_x_xxxzzzz[k];

                g_x_x_xz_xxyyyy[k] = -g_x_x_x_xxyyyy[k] * ab_z + g_x_x_x_xxyyyyz[k];

                g_x_x_xz_xxyyyz[k] = -g_x_x_x_xxyyyz[k] * ab_z + g_x_x_x_xxyyyzz[k];

                g_x_x_xz_xxyyzz[k] = -g_x_x_x_xxyyzz[k] * ab_z + g_x_x_x_xxyyzzz[k];

                g_x_x_xz_xxyzzz[k] = -g_x_x_x_xxyzzz[k] * ab_z + g_x_x_x_xxyzzzz[k];

                g_x_x_xz_xxzzzz[k] = -g_x_x_x_xxzzzz[k] * ab_z + g_x_x_x_xxzzzzz[k];

                g_x_x_xz_xyyyyy[k] = -g_x_x_x_xyyyyy[k] * ab_z + g_x_x_x_xyyyyyz[k];

                g_x_x_xz_xyyyyz[k] = -g_x_x_x_xyyyyz[k] * ab_z + g_x_x_x_xyyyyzz[k];

                g_x_x_xz_xyyyzz[k] = -g_x_x_x_xyyyzz[k] * ab_z + g_x_x_x_xyyyzzz[k];

                g_x_x_xz_xyyzzz[k] = -g_x_x_x_xyyzzz[k] * ab_z + g_x_x_x_xyyzzzz[k];

                g_x_x_xz_xyzzzz[k] = -g_x_x_x_xyzzzz[k] * ab_z + g_x_x_x_xyzzzzz[k];

                g_x_x_xz_xzzzzz[k] = -g_x_x_x_xzzzzz[k] * ab_z + g_x_x_x_xzzzzzz[k];

                g_x_x_xz_yyyyyy[k] = -g_x_x_x_yyyyyy[k] * ab_z + g_x_x_x_yyyyyyz[k];

                g_x_x_xz_yyyyyz[k] = -g_x_x_x_yyyyyz[k] * ab_z + g_x_x_x_yyyyyzz[k];

                g_x_x_xz_yyyyzz[k] = -g_x_x_x_yyyyzz[k] * ab_z + g_x_x_x_yyyyzzz[k];

                g_x_x_xz_yyyzzz[k] = -g_x_x_x_yyyzzz[k] * ab_z + g_x_x_x_yyyzzzz[k];

                g_x_x_xz_yyzzzz[k] = -g_x_x_x_yyzzzz[k] * ab_z + g_x_x_x_yyzzzzz[k];

                g_x_x_xz_yzzzzz[k] = -g_x_x_x_yzzzzz[k] * ab_z + g_x_x_x_yzzzzzz[k];

                g_x_x_xz_zzzzzz[k] = -g_x_x_x_zzzzzz[k] * ab_z + g_x_x_x_zzzzzzz[k];
            }

            /// Set up 84-112 components of targeted buffer : cbuffer.data(

            auto g_x_x_yy_xxxxxx = cbuffer.data(di_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_yy_xxxxxy = cbuffer.data(di_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_yy_xxxxxz = cbuffer.data(di_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_yy_xxxxyy = cbuffer.data(di_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_yy_xxxxyz = cbuffer.data(di_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_yy_xxxxzz = cbuffer.data(di_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_x_yy_xxxyyy = cbuffer.data(di_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_x_yy_xxxyyz = cbuffer.data(di_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_x_yy_xxxyzz = cbuffer.data(di_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_x_yy_xxxzzz = cbuffer.data(di_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_x_yy_xxyyyy = cbuffer.data(di_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_x_yy_xxyyyz = cbuffer.data(di_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_x_yy_xxyyzz = cbuffer.data(di_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_x_yy_xxyzzz = cbuffer.data(di_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_x_yy_xxzzzz = cbuffer.data(di_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_x_yy_xyyyyy = cbuffer.data(di_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_x_yy_xyyyyz = cbuffer.data(di_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_x_yy_xyyyzz = cbuffer.data(di_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_x_yy_xyyzzz = cbuffer.data(di_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_x_yy_xyzzzz = cbuffer.data(di_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_x_yy_xzzzzz = cbuffer.data(di_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_x_yy_yyyyyy = cbuffer.data(di_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_x_yy_yyyyyz = cbuffer.data(di_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_x_yy_yyyyzz = cbuffer.data(di_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_x_yy_yyyzzz = cbuffer.data(di_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_x_yy_yyzzzz = cbuffer.data(di_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_x_yy_yzzzzz = cbuffer.data(di_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_x_yy_zzzzzz = cbuffer.data(di_geom_11_off + 111 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_y_xxxxxx, g_x_x_y_xxxxxxy, g_x_x_y_xxxxxy, g_x_x_y_xxxxxyy, g_x_x_y_xxxxxyz, g_x_x_y_xxxxxz, g_x_x_y_xxxxyy, g_x_x_y_xxxxyyy, g_x_x_y_xxxxyyz, g_x_x_y_xxxxyz, g_x_x_y_xxxxyzz, g_x_x_y_xxxxzz, g_x_x_y_xxxyyy, g_x_x_y_xxxyyyy, g_x_x_y_xxxyyyz, g_x_x_y_xxxyyz, g_x_x_y_xxxyyzz, g_x_x_y_xxxyzz, g_x_x_y_xxxyzzz, g_x_x_y_xxxzzz, g_x_x_y_xxyyyy, g_x_x_y_xxyyyyy, g_x_x_y_xxyyyyz, g_x_x_y_xxyyyz, g_x_x_y_xxyyyzz, g_x_x_y_xxyyzz, g_x_x_y_xxyyzzz, g_x_x_y_xxyzzz, g_x_x_y_xxyzzzz, g_x_x_y_xxzzzz, g_x_x_y_xyyyyy, g_x_x_y_xyyyyyy, g_x_x_y_xyyyyyz, g_x_x_y_xyyyyz, g_x_x_y_xyyyyzz, g_x_x_y_xyyyzz, g_x_x_y_xyyyzzz, g_x_x_y_xyyzzz, g_x_x_y_xyyzzzz, g_x_x_y_xyzzzz, g_x_x_y_xyzzzzz, g_x_x_y_xzzzzz, g_x_x_y_yyyyyy, g_x_x_y_yyyyyyy, g_x_x_y_yyyyyyz, g_x_x_y_yyyyyz, g_x_x_y_yyyyyzz, g_x_x_y_yyyyzz, g_x_x_y_yyyyzzz, g_x_x_y_yyyzzz, g_x_x_y_yyyzzzz, g_x_x_y_yyzzzz, g_x_x_y_yyzzzzz, g_x_x_y_yzzzzz, g_x_x_y_yzzzzzz, g_x_x_y_zzzzzz, g_x_x_yy_xxxxxx, g_x_x_yy_xxxxxy, g_x_x_yy_xxxxxz, g_x_x_yy_xxxxyy, g_x_x_yy_xxxxyz, g_x_x_yy_xxxxzz, g_x_x_yy_xxxyyy, g_x_x_yy_xxxyyz, g_x_x_yy_xxxyzz, g_x_x_yy_xxxzzz, g_x_x_yy_xxyyyy, g_x_x_yy_xxyyyz, g_x_x_yy_xxyyzz, g_x_x_yy_xxyzzz, g_x_x_yy_xxzzzz, g_x_x_yy_xyyyyy, g_x_x_yy_xyyyyz, g_x_x_yy_xyyyzz, g_x_x_yy_xyyzzz, g_x_x_yy_xyzzzz, g_x_x_yy_xzzzzz, g_x_x_yy_yyyyyy, g_x_x_yy_yyyyyz, g_x_x_yy_yyyyzz, g_x_x_yy_yyyzzz, g_x_x_yy_yyzzzz, g_x_x_yy_yzzzzz, g_x_x_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yy_xxxxxx[k] = -g_x_x_y_xxxxxx[k] * ab_y + g_x_x_y_xxxxxxy[k];

                g_x_x_yy_xxxxxy[k] = -g_x_x_y_xxxxxy[k] * ab_y + g_x_x_y_xxxxxyy[k];

                g_x_x_yy_xxxxxz[k] = -g_x_x_y_xxxxxz[k] * ab_y + g_x_x_y_xxxxxyz[k];

                g_x_x_yy_xxxxyy[k] = -g_x_x_y_xxxxyy[k] * ab_y + g_x_x_y_xxxxyyy[k];

                g_x_x_yy_xxxxyz[k] = -g_x_x_y_xxxxyz[k] * ab_y + g_x_x_y_xxxxyyz[k];

                g_x_x_yy_xxxxzz[k] = -g_x_x_y_xxxxzz[k] * ab_y + g_x_x_y_xxxxyzz[k];

                g_x_x_yy_xxxyyy[k] = -g_x_x_y_xxxyyy[k] * ab_y + g_x_x_y_xxxyyyy[k];

                g_x_x_yy_xxxyyz[k] = -g_x_x_y_xxxyyz[k] * ab_y + g_x_x_y_xxxyyyz[k];

                g_x_x_yy_xxxyzz[k] = -g_x_x_y_xxxyzz[k] * ab_y + g_x_x_y_xxxyyzz[k];

                g_x_x_yy_xxxzzz[k] = -g_x_x_y_xxxzzz[k] * ab_y + g_x_x_y_xxxyzzz[k];

                g_x_x_yy_xxyyyy[k] = -g_x_x_y_xxyyyy[k] * ab_y + g_x_x_y_xxyyyyy[k];

                g_x_x_yy_xxyyyz[k] = -g_x_x_y_xxyyyz[k] * ab_y + g_x_x_y_xxyyyyz[k];

                g_x_x_yy_xxyyzz[k] = -g_x_x_y_xxyyzz[k] * ab_y + g_x_x_y_xxyyyzz[k];

                g_x_x_yy_xxyzzz[k] = -g_x_x_y_xxyzzz[k] * ab_y + g_x_x_y_xxyyzzz[k];

                g_x_x_yy_xxzzzz[k] = -g_x_x_y_xxzzzz[k] * ab_y + g_x_x_y_xxyzzzz[k];

                g_x_x_yy_xyyyyy[k] = -g_x_x_y_xyyyyy[k] * ab_y + g_x_x_y_xyyyyyy[k];

                g_x_x_yy_xyyyyz[k] = -g_x_x_y_xyyyyz[k] * ab_y + g_x_x_y_xyyyyyz[k];

                g_x_x_yy_xyyyzz[k] = -g_x_x_y_xyyyzz[k] * ab_y + g_x_x_y_xyyyyzz[k];

                g_x_x_yy_xyyzzz[k] = -g_x_x_y_xyyzzz[k] * ab_y + g_x_x_y_xyyyzzz[k];

                g_x_x_yy_xyzzzz[k] = -g_x_x_y_xyzzzz[k] * ab_y + g_x_x_y_xyyzzzz[k];

                g_x_x_yy_xzzzzz[k] = -g_x_x_y_xzzzzz[k] * ab_y + g_x_x_y_xyzzzzz[k];

                g_x_x_yy_yyyyyy[k] = -g_x_x_y_yyyyyy[k] * ab_y + g_x_x_y_yyyyyyy[k];

                g_x_x_yy_yyyyyz[k] = -g_x_x_y_yyyyyz[k] * ab_y + g_x_x_y_yyyyyyz[k];

                g_x_x_yy_yyyyzz[k] = -g_x_x_y_yyyyzz[k] * ab_y + g_x_x_y_yyyyyzz[k];

                g_x_x_yy_yyyzzz[k] = -g_x_x_y_yyyzzz[k] * ab_y + g_x_x_y_yyyyzzz[k];

                g_x_x_yy_yyzzzz[k] = -g_x_x_y_yyzzzz[k] * ab_y + g_x_x_y_yyyzzzz[k];

                g_x_x_yy_yzzzzz[k] = -g_x_x_y_yzzzzz[k] * ab_y + g_x_x_y_yyzzzzz[k];

                g_x_x_yy_zzzzzz[k] = -g_x_x_y_zzzzzz[k] * ab_y + g_x_x_y_yzzzzzz[k];
            }

            /// Set up 112-140 components of targeted buffer : cbuffer.data(

            auto g_x_x_yz_xxxxxx = cbuffer.data(di_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_x_yz_xxxxxy = cbuffer.data(di_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_x_yz_xxxxxz = cbuffer.data(di_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_x_yz_xxxxyy = cbuffer.data(di_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_x_yz_xxxxyz = cbuffer.data(di_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_x_yz_xxxxzz = cbuffer.data(di_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_x_yz_xxxyyy = cbuffer.data(di_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_x_yz_xxxyyz = cbuffer.data(di_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_x_yz_xxxyzz = cbuffer.data(di_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_x_yz_xxxzzz = cbuffer.data(di_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_x_yz_xxyyyy = cbuffer.data(di_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_x_yz_xxyyyz = cbuffer.data(di_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_x_yz_xxyyzz = cbuffer.data(di_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_x_yz_xxyzzz = cbuffer.data(di_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_x_yz_xxzzzz = cbuffer.data(di_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_x_yz_xyyyyy = cbuffer.data(di_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_x_yz_xyyyyz = cbuffer.data(di_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_x_yz_xyyyzz = cbuffer.data(di_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_x_yz_xyyzzz = cbuffer.data(di_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_x_yz_xyzzzz = cbuffer.data(di_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_x_yz_xzzzzz = cbuffer.data(di_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_x_yz_yyyyyy = cbuffer.data(di_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_x_yz_yyyyyz = cbuffer.data(di_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_x_yz_yyyyzz = cbuffer.data(di_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_x_yz_yyyzzz = cbuffer.data(di_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_x_yz_yyzzzz = cbuffer.data(di_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_x_yz_yzzzzz = cbuffer.data(di_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_x_yz_zzzzzz = cbuffer.data(di_geom_11_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yz_xxxxxx, g_x_x_yz_xxxxxy, g_x_x_yz_xxxxxz, g_x_x_yz_xxxxyy, g_x_x_yz_xxxxyz, g_x_x_yz_xxxxzz, g_x_x_yz_xxxyyy, g_x_x_yz_xxxyyz, g_x_x_yz_xxxyzz, g_x_x_yz_xxxzzz, g_x_x_yz_xxyyyy, g_x_x_yz_xxyyyz, g_x_x_yz_xxyyzz, g_x_x_yz_xxyzzz, g_x_x_yz_xxzzzz, g_x_x_yz_xyyyyy, g_x_x_yz_xyyyyz, g_x_x_yz_xyyyzz, g_x_x_yz_xyyzzz, g_x_x_yz_xyzzzz, g_x_x_yz_xzzzzz, g_x_x_yz_yyyyyy, g_x_x_yz_yyyyyz, g_x_x_yz_yyyyzz, g_x_x_yz_yyyzzz, g_x_x_yz_yyzzzz, g_x_x_yz_yzzzzz, g_x_x_yz_zzzzzz, g_x_x_z_xxxxxx, g_x_x_z_xxxxxxy, g_x_x_z_xxxxxy, g_x_x_z_xxxxxyy, g_x_x_z_xxxxxyz, g_x_x_z_xxxxxz, g_x_x_z_xxxxyy, g_x_x_z_xxxxyyy, g_x_x_z_xxxxyyz, g_x_x_z_xxxxyz, g_x_x_z_xxxxyzz, g_x_x_z_xxxxzz, g_x_x_z_xxxyyy, g_x_x_z_xxxyyyy, g_x_x_z_xxxyyyz, g_x_x_z_xxxyyz, g_x_x_z_xxxyyzz, g_x_x_z_xxxyzz, g_x_x_z_xxxyzzz, g_x_x_z_xxxzzz, g_x_x_z_xxyyyy, g_x_x_z_xxyyyyy, g_x_x_z_xxyyyyz, g_x_x_z_xxyyyz, g_x_x_z_xxyyyzz, g_x_x_z_xxyyzz, g_x_x_z_xxyyzzz, g_x_x_z_xxyzzz, g_x_x_z_xxyzzzz, g_x_x_z_xxzzzz, g_x_x_z_xyyyyy, g_x_x_z_xyyyyyy, g_x_x_z_xyyyyyz, g_x_x_z_xyyyyz, g_x_x_z_xyyyyzz, g_x_x_z_xyyyzz, g_x_x_z_xyyyzzz, g_x_x_z_xyyzzz, g_x_x_z_xyyzzzz, g_x_x_z_xyzzzz, g_x_x_z_xyzzzzz, g_x_x_z_xzzzzz, g_x_x_z_yyyyyy, g_x_x_z_yyyyyyy, g_x_x_z_yyyyyyz, g_x_x_z_yyyyyz, g_x_x_z_yyyyyzz, g_x_x_z_yyyyzz, g_x_x_z_yyyyzzz, g_x_x_z_yyyzzz, g_x_x_z_yyyzzzz, g_x_x_z_yyzzzz, g_x_x_z_yyzzzzz, g_x_x_z_yzzzzz, g_x_x_z_yzzzzzz, g_x_x_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yz_xxxxxx[k] = -g_x_x_z_xxxxxx[k] * ab_y + g_x_x_z_xxxxxxy[k];

                g_x_x_yz_xxxxxy[k] = -g_x_x_z_xxxxxy[k] * ab_y + g_x_x_z_xxxxxyy[k];

                g_x_x_yz_xxxxxz[k] = -g_x_x_z_xxxxxz[k] * ab_y + g_x_x_z_xxxxxyz[k];

                g_x_x_yz_xxxxyy[k] = -g_x_x_z_xxxxyy[k] * ab_y + g_x_x_z_xxxxyyy[k];

                g_x_x_yz_xxxxyz[k] = -g_x_x_z_xxxxyz[k] * ab_y + g_x_x_z_xxxxyyz[k];

                g_x_x_yz_xxxxzz[k] = -g_x_x_z_xxxxzz[k] * ab_y + g_x_x_z_xxxxyzz[k];

                g_x_x_yz_xxxyyy[k] = -g_x_x_z_xxxyyy[k] * ab_y + g_x_x_z_xxxyyyy[k];

                g_x_x_yz_xxxyyz[k] = -g_x_x_z_xxxyyz[k] * ab_y + g_x_x_z_xxxyyyz[k];

                g_x_x_yz_xxxyzz[k] = -g_x_x_z_xxxyzz[k] * ab_y + g_x_x_z_xxxyyzz[k];

                g_x_x_yz_xxxzzz[k] = -g_x_x_z_xxxzzz[k] * ab_y + g_x_x_z_xxxyzzz[k];

                g_x_x_yz_xxyyyy[k] = -g_x_x_z_xxyyyy[k] * ab_y + g_x_x_z_xxyyyyy[k];

                g_x_x_yz_xxyyyz[k] = -g_x_x_z_xxyyyz[k] * ab_y + g_x_x_z_xxyyyyz[k];

                g_x_x_yz_xxyyzz[k] = -g_x_x_z_xxyyzz[k] * ab_y + g_x_x_z_xxyyyzz[k];

                g_x_x_yz_xxyzzz[k] = -g_x_x_z_xxyzzz[k] * ab_y + g_x_x_z_xxyyzzz[k];

                g_x_x_yz_xxzzzz[k] = -g_x_x_z_xxzzzz[k] * ab_y + g_x_x_z_xxyzzzz[k];

                g_x_x_yz_xyyyyy[k] = -g_x_x_z_xyyyyy[k] * ab_y + g_x_x_z_xyyyyyy[k];

                g_x_x_yz_xyyyyz[k] = -g_x_x_z_xyyyyz[k] * ab_y + g_x_x_z_xyyyyyz[k];

                g_x_x_yz_xyyyzz[k] = -g_x_x_z_xyyyzz[k] * ab_y + g_x_x_z_xyyyyzz[k];

                g_x_x_yz_xyyzzz[k] = -g_x_x_z_xyyzzz[k] * ab_y + g_x_x_z_xyyyzzz[k];

                g_x_x_yz_xyzzzz[k] = -g_x_x_z_xyzzzz[k] * ab_y + g_x_x_z_xyyzzzz[k];

                g_x_x_yz_xzzzzz[k] = -g_x_x_z_xzzzzz[k] * ab_y + g_x_x_z_xyzzzzz[k];

                g_x_x_yz_yyyyyy[k] = -g_x_x_z_yyyyyy[k] * ab_y + g_x_x_z_yyyyyyy[k];

                g_x_x_yz_yyyyyz[k] = -g_x_x_z_yyyyyz[k] * ab_y + g_x_x_z_yyyyyyz[k];

                g_x_x_yz_yyyyzz[k] = -g_x_x_z_yyyyzz[k] * ab_y + g_x_x_z_yyyyyzz[k];

                g_x_x_yz_yyyzzz[k] = -g_x_x_z_yyyzzz[k] * ab_y + g_x_x_z_yyyyzzz[k];

                g_x_x_yz_yyzzzz[k] = -g_x_x_z_yyzzzz[k] * ab_y + g_x_x_z_yyyzzzz[k];

                g_x_x_yz_yzzzzz[k] = -g_x_x_z_yzzzzz[k] * ab_y + g_x_x_z_yyzzzzz[k];

                g_x_x_yz_zzzzzz[k] = -g_x_x_z_zzzzzz[k] * ab_y + g_x_x_z_yzzzzzz[k];
            }

            /// Set up 140-168 components of targeted buffer : cbuffer.data(

            auto g_x_x_zz_xxxxxx = cbuffer.data(di_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_x_zz_xxxxxy = cbuffer.data(di_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_x_zz_xxxxxz = cbuffer.data(di_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_x_zz_xxxxyy = cbuffer.data(di_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_x_zz_xxxxyz = cbuffer.data(di_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_x_zz_xxxxzz = cbuffer.data(di_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_x_zz_xxxyyy = cbuffer.data(di_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_x_zz_xxxyyz = cbuffer.data(di_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_x_zz_xxxyzz = cbuffer.data(di_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_x_zz_xxxzzz = cbuffer.data(di_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_x_zz_xxyyyy = cbuffer.data(di_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_x_zz_xxyyyz = cbuffer.data(di_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_x_zz_xxyyzz = cbuffer.data(di_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_x_zz_xxyzzz = cbuffer.data(di_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_x_zz_xxzzzz = cbuffer.data(di_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_x_zz_xyyyyy = cbuffer.data(di_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_x_zz_xyyyyz = cbuffer.data(di_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_x_zz_xyyyzz = cbuffer.data(di_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_x_zz_xyyzzz = cbuffer.data(di_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_x_zz_xyzzzz = cbuffer.data(di_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_x_zz_xzzzzz = cbuffer.data(di_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_x_zz_yyyyyy = cbuffer.data(di_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_x_zz_yyyyyz = cbuffer.data(di_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_x_zz_yyyyzz = cbuffer.data(di_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_x_zz_yyyzzz = cbuffer.data(di_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_x_zz_yyzzzz = cbuffer.data(di_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_x_zz_yzzzzz = cbuffer.data(di_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_x_zz_zzzzzz = cbuffer.data(di_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_z_xxxxxx, g_x_x_z_xxxxxxz, g_x_x_z_xxxxxy, g_x_x_z_xxxxxyz, g_x_x_z_xxxxxz, g_x_x_z_xxxxxzz, g_x_x_z_xxxxyy, g_x_x_z_xxxxyyz, g_x_x_z_xxxxyz, g_x_x_z_xxxxyzz, g_x_x_z_xxxxzz, g_x_x_z_xxxxzzz, g_x_x_z_xxxyyy, g_x_x_z_xxxyyyz, g_x_x_z_xxxyyz, g_x_x_z_xxxyyzz, g_x_x_z_xxxyzz, g_x_x_z_xxxyzzz, g_x_x_z_xxxzzz, g_x_x_z_xxxzzzz, g_x_x_z_xxyyyy, g_x_x_z_xxyyyyz, g_x_x_z_xxyyyz, g_x_x_z_xxyyyzz, g_x_x_z_xxyyzz, g_x_x_z_xxyyzzz, g_x_x_z_xxyzzz, g_x_x_z_xxyzzzz, g_x_x_z_xxzzzz, g_x_x_z_xxzzzzz, g_x_x_z_xyyyyy, g_x_x_z_xyyyyyz, g_x_x_z_xyyyyz, g_x_x_z_xyyyyzz, g_x_x_z_xyyyzz, g_x_x_z_xyyyzzz, g_x_x_z_xyyzzz, g_x_x_z_xyyzzzz, g_x_x_z_xyzzzz, g_x_x_z_xyzzzzz, g_x_x_z_xzzzzz, g_x_x_z_xzzzzzz, g_x_x_z_yyyyyy, g_x_x_z_yyyyyyz, g_x_x_z_yyyyyz, g_x_x_z_yyyyyzz, g_x_x_z_yyyyzz, g_x_x_z_yyyyzzz, g_x_x_z_yyyzzz, g_x_x_z_yyyzzzz, g_x_x_z_yyzzzz, g_x_x_z_yyzzzzz, g_x_x_z_yzzzzz, g_x_x_z_yzzzzzz, g_x_x_z_zzzzzz, g_x_x_z_zzzzzzz, g_x_x_zz_xxxxxx, g_x_x_zz_xxxxxy, g_x_x_zz_xxxxxz, g_x_x_zz_xxxxyy, g_x_x_zz_xxxxyz, g_x_x_zz_xxxxzz, g_x_x_zz_xxxyyy, g_x_x_zz_xxxyyz, g_x_x_zz_xxxyzz, g_x_x_zz_xxxzzz, g_x_x_zz_xxyyyy, g_x_x_zz_xxyyyz, g_x_x_zz_xxyyzz, g_x_x_zz_xxyzzz, g_x_x_zz_xxzzzz, g_x_x_zz_xyyyyy, g_x_x_zz_xyyyyz, g_x_x_zz_xyyyzz, g_x_x_zz_xyyzzz, g_x_x_zz_xyzzzz, g_x_x_zz_xzzzzz, g_x_x_zz_yyyyyy, g_x_x_zz_yyyyyz, g_x_x_zz_yyyyzz, g_x_x_zz_yyyzzz, g_x_x_zz_yyzzzz, g_x_x_zz_yzzzzz, g_x_x_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zz_xxxxxx[k] = -g_x_x_z_xxxxxx[k] * ab_z + g_x_x_z_xxxxxxz[k];

                g_x_x_zz_xxxxxy[k] = -g_x_x_z_xxxxxy[k] * ab_z + g_x_x_z_xxxxxyz[k];

                g_x_x_zz_xxxxxz[k] = -g_x_x_z_xxxxxz[k] * ab_z + g_x_x_z_xxxxxzz[k];

                g_x_x_zz_xxxxyy[k] = -g_x_x_z_xxxxyy[k] * ab_z + g_x_x_z_xxxxyyz[k];

                g_x_x_zz_xxxxyz[k] = -g_x_x_z_xxxxyz[k] * ab_z + g_x_x_z_xxxxyzz[k];

                g_x_x_zz_xxxxzz[k] = -g_x_x_z_xxxxzz[k] * ab_z + g_x_x_z_xxxxzzz[k];

                g_x_x_zz_xxxyyy[k] = -g_x_x_z_xxxyyy[k] * ab_z + g_x_x_z_xxxyyyz[k];

                g_x_x_zz_xxxyyz[k] = -g_x_x_z_xxxyyz[k] * ab_z + g_x_x_z_xxxyyzz[k];

                g_x_x_zz_xxxyzz[k] = -g_x_x_z_xxxyzz[k] * ab_z + g_x_x_z_xxxyzzz[k];

                g_x_x_zz_xxxzzz[k] = -g_x_x_z_xxxzzz[k] * ab_z + g_x_x_z_xxxzzzz[k];

                g_x_x_zz_xxyyyy[k] = -g_x_x_z_xxyyyy[k] * ab_z + g_x_x_z_xxyyyyz[k];

                g_x_x_zz_xxyyyz[k] = -g_x_x_z_xxyyyz[k] * ab_z + g_x_x_z_xxyyyzz[k];

                g_x_x_zz_xxyyzz[k] = -g_x_x_z_xxyyzz[k] * ab_z + g_x_x_z_xxyyzzz[k];

                g_x_x_zz_xxyzzz[k] = -g_x_x_z_xxyzzz[k] * ab_z + g_x_x_z_xxyzzzz[k];

                g_x_x_zz_xxzzzz[k] = -g_x_x_z_xxzzzz[k] * ab_z + g_x_x_z_xxzzzzz[k];

                g_x_x_zz_xyyyyy[k] = -g_x_x_z_xyyyyy[k] * ab_z + g_x_x_z_xyyyyyz[k];

                g_x_x_zz_xyyyyz[k] = -g_x_x_z_xyyyyz[k] * ab_z + g_x_x_z_xyyyyzz[k];

                g_x_x_zz_xyyyzz[k] = -g_x_x_z_xyyyzz[k] * ab_z + g_x_x_z_xyyyzzz[k];

                g_x_x_zz_xyyzzz[k] = -g_x_x_z_xyyzzz[k] * ab_z + g_x_x_z_xyyzzzz[k];

                g_x_x_zz_xyzzzz[k] = -g_x_x_z_xyzzzz[k] * ab_z + g_x_x_z_xyzzzzz[k];

                g_x_x_zz_xzzzzz[k] = -g_x_x_z_xzzzzz[k] * ab_z + g_x_x_z_xzzzzzz[k];

                g_x_x_zz_yyyyyy[k] = -g_x_x_z_yyyyyy[k] * ab_z + g_x_x_z_yyyyyyz[k];

                g_x_x_zz_yyyyyz[k] = -g_x_x_z_yyyyyz[k] * ab_z + g_x_x_z_yyyyyzz[k];

                g_x_x_zz_yyyyzz[k] = -g_x_x_z_yyyyzz[k] * ab_z + g_x_x_z_yyyyzzz[k];

                g_x_x_zz_yyyzzz[k] = -g_x_x_z_yyyzzz[k] * ab_z + g_x_x_z_yyyzzzz[k];

                g_x_x_zz_yyzzzz[k] = -g_x_x_z_yyzzzz[k] * ab_z + g_x_x_z_yyzzzzz[k];

                g_x_x_zz_yzzzzz[k] = -g_x_x_z_yzzzzz[k] * ab_z + g_x_x_z_yzzzzzz[k];

                g_x_x_zz_zzzzzz[k] = -g_x_x_z_zzzzzz[k] * ab_z + g_x_x_z_zzzzzzz[k];
            }

            /// Set up 168-196 components of targeted buffer : cbuffer.data(

            auto g_x_y_xx_xxxxxx = cbuffer.data(di_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_xx_xxxxxy = cbuffer.data(di_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_xx_xxxxxz = cbuffer.data(di_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_xx_xxxxyy = cbuffer.data(di_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_xx_xxxxyz = cbuffer.data(di_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_xx_xxxxzz = cbuffer.data(di_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_y_xx_xxxyyy = cbuffer.data(di_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_xx_xxxyyz = cbuffer.data(di_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_xx_xxxyzz = cbuffer.data(di_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_xx_xxxzzz = cbuffer.data(di_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_xx_xxyyyy = cbuffer.data(di_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_xx_xxyyyz = cbuffer.data(di_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_y_xx_xxyyzz = cbuffer.data(di_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_y_xx_xxyzzz = cbuffer.data(di_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_y_xx_xxzzzz = cbuffer.data(di_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_y_xx_xyyyyy = cbuffer.data(di_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_y_xx_xyyyyz = cbuffer.data(di_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_y_xx_xyyyzz = cbuffer.data(di_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_y_xx_xyyzzz = cbuffer.data(di_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_y_xx_xyzzzz = cbuffer.data(di_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_y_xx_xzzzzz = cbuffer.data(di_geom_11_off + 188 * ccomps * dcomps);

            auto g_x_y_xx_yyyyyy = cbuffer.data(di_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_y_xx_yyyyyz = cbuffer.data(di_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_y_xx_yyyyzz = cbuffer.data(di_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_y_xx_yyyzzz = cbuffer.data(di_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_y_xx_yyzzzz = cbuffer.data(di_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_y_xx_yzzzzz = cbuffer.data(di_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_y_xx_zzzzzz = cbuffer.data(di_geom_11_off + 195 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_x_xxxxxx, g_0_y_x_xxxxxy, g_0_y_x_xxxxxz, g_0_y_x_xxxxyy, g_0_y_x_xxxxyz, g_0_y_x_xxxxzz, g_0_y_x_xxxyyy, g_0_y_x_xxxyyz, g_0_y_x_xxxyzz, g_0_y_x_xxxzzz, g_0_y_x_xxyyyy, g_0_y_x_xxyyyz, g_0_y_x_xxyyzz, g_0_y_x_xxyzzz, g_0_y_x_xxzzzz, g_0_y_x_xyyyyy, g_0_y_x_xyyyyz, g_0_y_x_xyyyzz, g_0_y_x_xyyzzz, g_0_y_x_xyzzzz, g_0_y_x_xzzzzz, g_0_y_x_yyyyyy, g_0_y_x_yyyyyz, g_0_y_x_yyyyzz, g_0_y_x_yyyzzz, g_0_y_x_yyzzzz, g_0_y_x_yzzzzz, g_0_y_x_zzzzzz, g_x_y_x_xxxxxx, g_x_y_x_xxxxxxx, g_x_y_x_xxxxxxy, g_x_y_x_xxxxxxz, g_x_y_x_xxxxxy, g_x_y_x_xxxxxyy, g_x_y_x_xxxxxyz, g_x_y_x_xxxxxz, g_x_y_x_xxxxxzz, g_x_y_x_xxxxyy, g_x_y_x_xxxxyyy, g_x_y_x_xxxxyyz, g_x_y_x_xxxxyz, g_x_y_x_xxxxyzz, g_x_y_x_xxxxzz, g_x_y_x_xxxxzzz, g_x_y_x_xxxyyy, g_x_y_x_xxxyyyy, g_x_y_x_xxxyyyz, g_x_y_x_xxxyyz, g_x_y_x_xxxyyzz, g_x_y_x_xxxyzz, g_x_y_x_xxxyzzz, g_x_y_x_xxxzzz, g_x_y_x_xxxzzzz, g_x_y_x_xxyyyy, g_x_y_x_xxyyyyy, g_x_y_x_xxyyyyz, g_x_y_x_xxyyyz, g_x_y_x_xxyyyzz, g_x_y_x_xxyyzz, g_x_y_x_xxyyzzz, g_x_y_x_xxyzzz, g_x_y_x_xxyzzzz, g_x_y_x_xxzzzz, g_x_y_x_xxzzzzz, g_x_y_x_xyyyyy, g_x_y_x_xyyyyyy, g_x_y_x_xyyyyyz, g_x_y_x_xyyyyz, g_x_y_x_xyyyyzz, g_x_y_x_xyyyzz, g_x_y_x_xyyyzzz, g_x_y_x_xyyzzz, g_x_y_x_xyyzzzz, g_x_y_x_xyzzzz, g_x_y_x_xyzzzzz, g_x_y_x_xzzzzz, g_x_y_x_xzzzzzz, g_x_y_x_yyyyyy, g_x_y_x_yyyyyz, g_x_y_x_yyyyzz, g_x_y_x_yyyzzz, g_x_y_x_yyzzzz, g_x_y_x_yzzzzz, g_x_y_x_zzzzzz, g_x_y_xx_xxxxxx, g_x_y_xx_xxxxxy, g_x_y_xx_xxxxxz, g_x_y_xx_xxxxyy, g_x_y_xx_xxxxyz, g_x_y_xx_xxxxzz, g_x_y_xx_xxxyyy, g_x_y_xx_xxxyyz, g_x_y_xx_xxxyzz, g_x_y_xx_xxxzzz, g_x_y_xx_xxyyyy, g_x_y_xx_xxyyyz, g_x_y_xx_xxyyzz, g_x_y_xx_xxyzzz, g_x_y_xx_xxzzzz, g_x_y_xx_xyyyyy, g_x_y_xx_xyyyyz, g_x_y_xx_xyyyzz, g_x_y_xx_xyyzzz, g_x_y_xx_xyzzzz, g_x_y_xx_xzzzzz, g_x_y_xx_yyyyyy, g_x_y_xx_yyyyyz, g_x_y_xx_yyyyzz, g_x_y_xx_yyyzzz, g_x_y_xx_yyzzzz, g_x_y_xx_yzzzzz, g_x_y_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xx_xxxxxx[k] = -g_0_y_x_xxxxxx[k] - g_x_y_x_xxxxxx[k] * ab_x + g_x_y_x_xxxxxxx[k];

                g_x_y_xx_xxxxxy[k] = -g_0_y_x_xxxxxy[k] - g_x_y_x_xxxxxy[k] * ab_x + g_x_y_x_xxxxxxy[k];

                g_x_y_xx_xxxxxz[k] = -g_0_y_x_xxxxxz[k] - g_x_y_x_xxxxxz[k] * ab_x + g_x_y_x_xxxxxxz[k];

                g_x_y_xx_xxxxyy[k] = -g_0_y_x_xxxxyy[k] - g_x_y_x_xxxxyy[k] * ab_x + g_x_y_x_xxxxxyy[k];

                g_x_y_xx_xxxxyz[k] = -g_0_y_x_xxxxyz[k] - g_x_y_x_xxxxyz[k] * ab_x + g_x_y_x_xxxxxyz[k];

                g_x_y_xx_xxxxzz[k] = -g_0_y_x_xxxxzz[k] - g_x_y_x_xxxxzz[k] * ab_x + g_x_y_x_xxxxxzz[k];

                g_x_y_xx_xxxyyy[k] = -g_0_y_x_xxxyyy[k] - g_x_y_x_xxxyyy[k] * ab_x + g_x_y_x_xxxxyyy[k];

                g_x_y_xx_xxxyyz[k] = -g_0_y_x_xxxyyz[k] - g_x_y_x_xxxyyz[k] * ab_x + g_x_y_x_xxxxyyz[k];

                g_x_y_xx_xxxyzz[k] = -g_0_y_x_xxxyzz[k] - g_x_y_x_xxxyzz[k] * ab_x + g_x_y_x_xxxxyzz[k];

                g_x_y_xx_xxxzzz[k] = -g_0_y_x_xxxzzz[k] - g_x_y_x_xxxzzz[k] * ab_x + g_x_y_x_xxxxzzz[k];

                g_x_y_xx_xxyyyy[k] = -g_0_y_x_xxyyyy[k] - g_x_y_x_xxyyyy[k] * ab_x + g_x_y_x_xxxyyyy[k];

                g_x_y_xx_xxyyyz[k] = -g_0_y_x_xxyyyz[k] - g_x_y_x_xxyyyz[k] * ab_x + g_x_y_x_xxxyyyz[k];

                g_x_y_xx_xxyyzz[k] = -g_0_y_x_xxyyzz[k] - g_x_y_x_xxyyzz[k] * ab_x + g_x_y_x_xxxyyzz[k];

                g_x_y_xx_xxyzzz[k] = -g_0_y_x_xxyzzz[k] - g_x_y_x_xxyzzz[k] * ab_x + g_x_y_x_xxxyzzz[k];

                g_x_y_xx_xxzzzz[k] = -g_0_y_x_xxzzzz[k] - g_x_y_x_xxzzzz[k] * ab_x + g_x_y_x_xxxzzzz[k];

                g_x_y_xx_xyyyyy[k] = -g_0_y_x_xyyyyy[k] - g_x_y_x_xyyyyy[k] * ab_x + g_x_y_x_xxyyyyy[k];

                g_x_y_xx_xyyyyz[k] = -g_0_y_x_xyyyyz[k] - g_x_y_x_xyyyyz[k] * ab_x + g_x_y_x_xxyyyyz[k];

                g_x_y_xx_xyyyzz[k] = -g_0_y_x_xyyyzz[k] - g_x_y_x_xyyyzz[k] * ab_x + g_x_y_x_xxyyyzz[k];

                g_x_y_xx_xyyzzz[k] = -g_0_y_x_xyyzzz[k] - g_x_y_x_xyyzzz[k] * ab_x + g_x_y_x_xxyyzzz[k];

                g_x_y_xx_xyzzzz[k] = -g_0_y_x_xyzzzz[k] - g_x_y_x_xyzzzz[k] * ab_x + g_x_y_x_xxyzzzz[k];

                g_x_y_xx_xzzzzz[k] = -g_0_y_x_xzzzzz[k] - g_x_y_x_xzzzzz[k] * ab_x + g_x_y_x_xxzzzzz[k];

                g_x_y_xx_yyyyyy[k] = -g_0_y_x_yyyyyy[k] - g_x_y_x_yyyyyy[k] * ab_x + g_x_y_x_xyyyyyy[k];

                g_x_y_xx_yyyyyz[k] = -g_0_y_x_yyyyyz[k] - g_x_y_x_yyyyyz[k] * ab_x + g_x_y_x_xyyyyyz[k];

                g_x_y_xx_yyyyzz[k] = -g_0_y_x_yyyyzz[k] - g_x_y_x_yyyyzz[k] * ab_x + g_x_y_x_xyyyyzz[k];

                g_x_y_xx_yyyzzz[k] = -g_0_y_x_yyyzzz[k] - g_x_y_x_yyyzzz[k] * ab_x + g_x_y_x_xyyyzzz[k];

                g_x_y_xx_yyzzzz[k] = -g_0_y_x_yyzzzz[k] - g_x_y_x_yyzzzz[k] * ab_x + g_x_y_x_xyyzzzz[k];

                g_x_y_xx_yzzzzz[k] = -g_0_y_x_yzzzzz[k] - g_x_y_x_yzzzzz[k] * ab_x + g_x_y_x_xyzzzzz[k];

                g_x_y_xx_zzzzzz[k] = -g_0_y_x_zzzzzz[k] - g_x_y_x_zzzzzz[k] * ab_x + g_x_y_x_xzzzzzz[k];
            }

            /// Set up 196-224 components of targeted buffer : cbuffer.data(

            auto g_x_y_xy_xxxxxx = cbuffer.data(di_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_y_xy_xxxxxy = cbuffer.data(di_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_y_xy_xxxxxz = cbuffer.data(di_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_y_xy_xxxxyy = cbuffer.data(di_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_y_xy_xxxxyz = cbuffer.data(di_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_y_xy_xxxxzz = cbuffer.data(di_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_y_xy_xxxyyy = cbuffer.data(di_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_y_xy_xxxyyz = cbuffer.data(di_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_y_xy_xxxyzz = cbuffer.data(di_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_y_xy_xxxzzz = cbuffer.data(di_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_y_xy_xxyyyy = cbuffer.data(di_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_y_xy_xxyyyz = cbuffer.data(di_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_y_xy_xxyyzz = cbuffer.data(di_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_y_xy_xxyzzz = cbuffer.data(di_geom_11_off + 209 * ccomps * dcomps);

            auto g_x_y_xy_xxzzzz = cbuffer.data(di_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_y_xy_xyyyyy = cbuffer.data(di_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_y_xy_xyyyyz = cbuffer.data(di_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_y_xy_xyyyzz = cbuffer.data(di_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_y_xy_xyyzzz = cbuffer.data(di_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_y_xy_xyzzzz = cbuffer.data(di_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_y_xy_xzzzzz = cbuffer.data(di_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_y_xy_yyyyyy = cbuffer.data(di_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_y_xy_yyyyyz = cbuffer.data(di_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_y_xy_yyyyzz = cbuffer.data(di_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_y_xy_yyyzzz = cbuffer.data(di_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_y_xy_yyzzzz = cbuffer.data(di_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_y_xy_yzzzzz = cbuffer.data(di_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_y_xy_zzzzzz = cbuffer.data(di_geom_11_off + 223 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxxxxx, g_0_y_y_xxxxxy, g_0_y_y_xxxxxz, g_0_y_y_xxxxyy, g_0_y_y_xxxxyz, g_0_y_y_xxxxzz, g_0_y_y_xxxyyy, g_0_y_y_xxxyyz, g_0_y_y_xxxyzz, g_0_y_y_xxxzzz, g_0_y_y_xxyyyy, g_0_y_y_xxyyyz, g_0_y_y_xxyyzz, g_0_y_y_xxyzzz, g_0_y_y_xxzzzz, g_0_y_y_xyyyyy, g_0_y_y_xyyyyz, g_0_y_y_xyyyzz, g_0_y_y_xyyzzz, g_0_y_y_xyzzzz, g_0_y_y_xzzzzz, g_0_y_y_yyyyyy, g_0_y_y_yyyyyz, g_0_y_y_yyyyzz, g_0_y_y_yyyzzz, g_0_y_y_yyzzzz, g_0_y_y_yzzzzz, g_0_y_y_zzzzzz, g_x_y_xy_xxxxxx, g_x_y_xy_xxxxxy, g_x_y_xy_xxxxxz, g_x_y_xy_xxxxyy, g_x_y_xy_xxxxyz, g_x_y_xy_xxxxzz, g_x_y_xy_xxxyyy, g_x_y_xy_xxxyyz, g_x_y_xy_xxxyzz, g_x_y_xy_xxxzzz, g_x_y_xy_xxyyyy, g_x_y_xy_xxyyyz, g_x_y_xy_xxyyzz, g_x_y_xy_xxyzzz, g_x_y_xy_xxzzzz, g_x_y_xy_xyyyyy, g_x_y_xy_xyyyyz, g_x_y_xy_xyyyzz, g_x_y_xy_xyyzzz, g_x_y_xy_xyzzzz, g_x_y_xy_xzzzzz, g_x_y_xy_yyyyyy, g_x_y_xy_yyyyyz, g_x_y_xy_yyyyzz, g_x_y_xy_yyyzzz, g_x_y_xy_yyzzzz, g_x_y_xy_yzzzzz, g_x_y_xy_zzzzzz, g_x_y_y_xxxxxx, g_x_y_y_xxxxxxx, g_x_y_y_xxxxxxy, g_x_y_y_xxxxxxz, g_x_y_y_xxxxxy, g_x_y_y_xxxxxyy, g_x_y_y_xxxxxyz, g_x_y_y_xxxxxz, g_x_y_y_xxxxxzz, g_x_y_y_xxxxyy, g_x_y_y_xxxxyyy, g_x_y_y_xxxxyyz, g_x_y_y_xxxxyz, g_x_y_y_xxxxyzz, g_x_y_y_xxxxzz, g_x_y_y_xxxxzzz, g_x_y_y_xxxyyy, g_x_y_y_xxxyyyy, g_x_y_y_xxxyyyz, g_x_y_y_xxxyyz, g_x_y_y_xxxyyzz, g_x_y_y_xxxyzz, g_x_y_y_xxxyzzz, g_x_y_y_xxxzzz, g_x_y_y_xxxzzzz, g_x_y_y_xxyyyy, g_x_y_y_xxyyyyy, g_x_y_y_xxyyyyz, g_x_y_y_xxyyyz, g_x_y_y_xxyyyzz, g_x_y_y_xxyyzz, g_x_y_y_xxyyzzz, g_x_y_y_xxyzzz, g_x_y_y_xxyzzzz, g_x_y_y_xxzzzz, g_x_y_y_xxzzzzz, g_x_y_y_xyyyyy, g_x_y_y_xyyyyyy, g_x_y_y_xyyyyyz, g_x_y_y_xyyyyz, g_x_y_y_xyyyyzz, g_x_y_y_xyyyzz, g_x_y_y_xyyyzzz, g_x_y_y_xyyzzz, g_x_y_y_xyyzzzz, g_x_y_y_xyzzzz, g_x_y_y_xyzzzzz, g_x_y_y_xzzzzz, g_x_y_y_xzzzzzz, g_x_y_y_yyyyyy, g_x_y_y_yyyyyz, g_x_y_y_yyyyzz, g_x_y_y_yyyzzz, g_x_y_y_yyzzzz, g_x_y_y_yzzzzz, g_x_y_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xy_xxxxxx[k] = -g_0_y_y_xxxxxx[k] - g_x_y_y_xxxxxx[k] * ab_x + g_x_y_y_xxxxxxx[k];

                g_x_y_xy_xxxxxy[k] = -g_0_y_y_xxxxxy[k] - g_x_y_y_xxxxxy[k] * ab_x + g_x_y_y_xxxxxxy[k];

                g_x_y_xy_xxxxxz[k] = -g_0_y_y_xxxxxz[k] - g_x_y_y_xxxxxz[k] * ab_x + g_x_y_y_xxxxxxz[k];

                g_x_y_xy_xxxxyy[k] = -g_0_y_y_xxxxyy[k] - g_x_y_y_xxxxyy[k] * ab_x + g_x_y_y_xxxxxyy[k];

                g_x_y_xy_xxxxyz[k] = -g_0_y_y_xxxxyz[k] - g_x_y_y_xxxxyz[k] * ab_x + g_x_y_y_xxxxxyz[k];

                g_x_y_xy_xxxxzz[k] = -g_0_y_y_xxxxzz[k] - g_x_y_y_xxxxzz[k] * ab_x + g_x_y_y_xxxxxzz[k];

                g_x_y_xy_xxxyyy[k] = -g_0_y_y_xxxyyy[k] - g_x_y_y_xxxyyy[k] * ab_x + g_x_y_y_xxxxyyy[k];

                g_x_y_xy_xxxyyz[k] = -g_0_y_y_xxxyyz[k] - g_x_y_y_xxxyyz[k] * ab_x + g_x_y_y_xxxxyyz[k];

                g_x_y_xy_xxxyzz[k] = -g_0_y_y_xxxyzz[k] - g_x_y_y_xxxyzz[k] * ab_x + g_x_y_y_xxxxyzz[k];

                g_x_y_xy_xxxzzz[k] = -g_0_y_y_xxxzzz[k] - g_x_y_y_xxxzzz[k] * ab_x + g_x_y_y_xxxxzzz[k];

                g_x_y_xy_xxyyyy[k] = -g_0_y_y_xxyyyy[k] - g_x_y_y_xxyyyy[k] * ab_x + g_x_y_y_xxxyyyy[k];

                g_x_y_xy_xxyyyz[k] = -g_0_y_y_xxyyyz[k] - g_x_y_y_xxyyyz[k] * ab_x + g_x_y_y_xxxyyyz[k];

                g_x_y_xy_xxyyzz[k] = -g_0_y_y_xxyyzz[k] - g_x_y_y_xxyyzz[k] * ab_x + g_x_y_y_xxxyyzz[k];

                g_x_y_xy_xxyzzz[k] = -g_0_y_y_xxyzzz[k] - g_x_y_y_xxyzzz[k] * ab_x + g_x_y_y_xxxyzzz[k];

                g_x_y_xy_xxzzzz[k] = -g_0_y_y_xxzzzz[k] - g_x_y_y_xxzzzz[k] * ab_x + g_x_y_y_xxxzzzz[k];

                g_x_y_xy_xyyyyy[k] = -g_0_y_y_xyyyyy[k] - g_x_y_y_xyyyyy[k] * ab_x + g_x_y_y_xxyyyyy[k];

                g_x_y_xy_xyyyyz[k] = -g_0_y_y_xyyyyz[k] - g_x_y_y_xyyyyz[k] * ab_x + g_x_y_y_xxyyyyz[k];

                g_x_y_xy_xyyyzz[k] = -g_0_y_y_xyyyzz[k] - g_x_y_y_xyyyzz[k] * ab_x + g_x_y_y_xxyyyzz[k];

                g_x_y_xy_xyyzzz[k] = -g_0_y_y_xyyzzz[k] - g_x_y_y_xyyzzz[k] * ab_x + g_x_y_y_xxyyzzz[k];

                g_x_y_xy_xyzzzz[k] = -g_0_y_y_xyzzzz[k] - g_x_y_y_xyzzzz[k] * ab_x + g_x_y_y_xxyzzzz[k];

                g_x_y_xy_xzzzzz[k] = -g_0_y_y_xzzzzz[k] - g_x_y_y_xzzzzz[k] * ab_x + g_x_y_y_xxzzzzz[k];

                g_x_y_xy_yyyyyy[k] = -g_0_y_y_yyyyyy[k] - g_x_y_y_yyyyyy[k] * ab_x + g_x_y_y_xyyyyyy[k];

                g_x_y_xy_yyyyyz[k] = -g_0_y_y_yyyyyz[k] - g_x_y_y_yyyyyz[k] * ab_x + g_x_y_y_xyyyyyz[k];

                g_x_y_xy_yyyyzz[k] = -g_0_y_y_yyyyzz[k] - g_x_y_y_yyyyzz[k] * ab_x + g_x_y_y_xyyyyzz[k];

                g_x_y_xy_yyyzzz[k] = -g_0_y_y_yyyzzz[k] - g_x_y_y_yyyzzz[k] * ab_x + g_x_y_y_xyyyzzz[k];

                g_x_y_xy_yyzzzz[k] = -g_0_y_y_yyzzzz[k] - g_x_y_y_yyzzzz[k] * ab_x + g_x_y_y_xyyzzzz[k];

                g_x_y_xy_yzzzzz[k] = -g_0_y_y_yzzzzz[k] - g_x_y_y_yzzzzz[k] * ab_x + g_x_y_y_xyzzzzz[k];

                g_x_y_xy_zzzzzz[k] = -g_0_y_y_zzzzzz[k] - g_x_y_y_zzzzzz[k] * ab_x + g_x_y_y_xzzzzzz[k];
            }

            /// Set up 224-252 components of targeted buffer : cbuffer.data(

            auto g_x_y_xz_xxxxxx = cbuffer.data(di_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_y_xz_xxxxxy = cbuffer.data(di_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_y_xz_xxxxxz = cbuffer.data(di_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_y_xz_xxxxyy = cbuffer.data(di_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_y_xz_xxxxyz = cbuffer.data(di_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_y_xz_xxxxzz = cbuffer.data(di_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_y_xz_xxxyyy = cbuffer.data(di_geom_11_off + 230 * ccomps * dcomps);

            auto g_x_y_xz_xxxyyz = cbuffer.data(di_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_y_xz_xxxyzz = cbuffer.data(di_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_y_xz_xxxzzz = cbuffer.data(di_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_y_xz_xxyyyy = cbuffer.data(di_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_y_xz_xxyyyz = cbuffer.data(di_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_y_xz_xxyyzz = cbuffer.data(di_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_y_xz_xxyzzz = cbuffer.data(di_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_y_xz_xxzzzz = cbuffer.data(di_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_y_xz_xyyyyy = cbuffer.data(di_geom_11_off + 239 * ccomps * dcomps);

            auto g_x_y_xz_xyyyyz = cbuffer.data(di_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_y_xz_xyyyzz = cbuffer.data(di_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_y_xz_xyyzzz = cbuffer.data(di_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_y_xz_xyzzzz = cbuffer.data(di_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_y_xz_xzzzzz = cbuffer.data(di_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_y_xz_yyyyyy = cbuffer.data(di_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_y_xz_yyyyyz = cbuffer.data(di_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_y_xz_yyyyzz = cbuffer.data(di_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_y_xz_yyyzzz = cbuffer.data(di_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_y_xz_yyzzzz = cbuffer.data(di_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_y_xz_yzzzzz = cbuffer.data(di_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_y_xz_zzzzzz = cbuffer.data(di_geom_11_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_x_xxxxxx, g_x_y_x_xxxxxxz, g_x_y_x_xxxxxy, g_x_y_x_xxxxxyz, g_x_y_x_xxxxxz, g_x_y_x_xxxxxzz, g_x_y_x_xxxxyy, g_x_y_x_xxxxyyz, g_x_y_x_xxxxyz, g_x_y_x_xxxxyzz, g_x_y_x_xxxxzz, g_x_y_x_xxxxzzz, g_x_y_x_xxxyyy, g_x_y_x_xxxyyyz, g_x_y_x_xxxyyz, g_x_y_x_xxxyyzz, g_x_y_x_xxxyzz, g_x_y_x_xxxyzzz, g_x_y_x_xxxzzz, g_x_y_x_xxxzzzz, g_x_y_x_xxyyyy, g_x_y_x_xxyyyyz, g_x_y_x_xxyyyz, g_x_y_x_xxyyyzz, g_x_y_x_xxyyzz, g_x_y_x_xxyyzzz, g_x_y_x_xxyzzz, g_x_y_x_xxyzzzz, g_x_y_x_xxzzzz, g_x_y_x_xxzzzzz, g_x_y_x_xyyyyy, g_x_y_x_xyyyyyz, g_x_y_x_xyyyyz, g_x_y_x_xyyyyzz, g_x_y_x_xyyyzz, g_x_y_x_xyyyzzz, g_x_y_x_xyyzzz, g_x_y_x_xyyzzzz, g_x_y_x_xyzzzz, g_x_y_x_xyzzzzz, g_x_y_x_xzzzzz, g_x_y_x_xzzzzzz, g_x_y_x_yyyyyy, g_x_y_x_yyyyyyz, g_x_y_x_yyyyyz, g_x_y_x_yyyyyzz, g_x_y_x_yyyyzz, g_x_y_x_yyyyzzz, g_x_y_x_yyyzzz, g_x_y_x_yyyzzzz, g_x_y_x_yyzzzz, g_x_y_x_yyzzzzz, g_x_y_x_yzzzzz, g_x_y_x_yzzzzzz, g_x_y_x_zzzzzz, g_x_y_x_zzzzzzz, g_x_y_xz_xxxxxx, g_x_y_xz_xxxxxy, g_x_y_xz_xxxxxz, g_x_y_xz_xxxxyy, g_x_y_xz_xxxxyz, g_x_y_xz_xxxxzz, g_x_y_xz_xxxyyy, g_x_y_xz_xxxyyz, g_x_y_xz_xxxyzz, g_x_y_xz_xxxzzz, g_x_y_xz_xxyyyy, g_x_y_xz_xxyyyz, g_x_y_xz_xxyyzz, g_x_y_xz_xxyzzz, g_x_y_xz_xxzzzz, g_x_y_xz_xyyyyy, g_x_y_xz_xyyyyz, g_x_y_xz_xyyyzz, g_x_y_xz_xyyzzz, g_x_y_xz_xyzzzz, g_x_y_xz_xzzzzz, g_x_y_xz_yyyyyy, g_x_y_xz_yyyyyz, g_x_y_xz_yyyyzz, g_x_y_xz_yyyzzz, g_x_y_xz_yyzzzz, g_x_y_xz_yzzzzz, g_x_y_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xz_xxxxxx[k] = -g_x_y_x_xxxxxx[k] * ab_z + g_x_y_x_xxxxxxz[k];

                g_x_y_xz_xxxxxy[k] = -g_x_y_x_xxxxxy[k] * ab_z + g_x_y_x_xxxxxyz[k];

                g_x_y_xz_xxxxxz[k] = -g_x_y_x_xxxxxz[k] * ab_z + g_x_y_x_xxxxxzz[k];

                g_x_y_xz_xxxxyy[k] = -g_x_y_x_xxxxyy[k] * ab_z + g_x_y_x_xxxxyyz[k];

                g_x_y_xz_xxxxyz[k] = -g_x_y_x_xxxxyz[k] * ab_z + g_x_y_x_xxxxyzz[k];

                g_x_y_xz_xxxxzz[k] = -g_x_y_x_xxxxzz[k] * ab_z + g_x_y_x_xxxxzzz[k];

                g_x_y_xz_xxxyyy[k] = -g_x_y_x_xxxyyy[k] * ab_z + g_x_y_x_xxxyyyz[k];

                g_x_y_xz_xxxyyz[k] = -g_x_y_x_xxxyyz[k] * ab_z + g_x_y_x_xxxyyzz[k];

                g_x_y_xz_xxxyzz[k] = -g_x_y_x_xxxyzz[k] * ab_z + g_x_y_x_xxxyzzz[k];

                g_x_y_xz_xxxzzz[k] = -g_x_y_x_xxxzzz[k] * ab_z + g_x_y_x_xxxzzzz[k];

                g_x_y_xz_xxyyyy[k] = -g_x_y_x_xxyyyy[k] * ab_z + g_x_y_x_xxyyyyz[k];

                g_x_y_xz_xxyyyz[k] = -g_x_y_x_xxyyyz[k] * ab_z + g_x_y_x_xxyyyzz[k];

                g_x_y_xz_xxyyzz[k] = -g_x_y_x_xxyyzz[k] * ab_z + g_x_y_x_xxyyzzz[k];

                g_x_y_xz_xxyzzz[k] = -g_x_y_x_xxyzzz[k] * ab_z + g_x_y_x_xxyzzzz[k];

                g_x_y_xz_xxzzzz[k] = -g_x_y_x_xxzzzz[k] * ab_z + g_x_y_x_xxzzzzz[k];

                g_x_y_xz_xyyyyy[k] = -g_x_y_x_xyyyyy[k] * ab_z + g_x_y_x_xyyyyyz[k];

                g_x_y_xz_xyyyyz[k] = -g_x_y_x_xyyyyz[k] * ab_z + g_x_y_x_xyyyyzz[k];

                g_x_y_xz_xyyyzz[k] = -g_x_y_x_xyyyzz[k] * ab_z + g_x_y_x_xyyyzzz[k];

                g_x_y_xz_xyyzzz[k] = -g_x_y_x_xyyzzz[k] * ab_z + g_x_y_x_xyyzzzz[k];

                g_x_y_xz_xyzzzz[k] = -g_x_y_x_xyzzzz[k] * ab_z + g_x_y_x_xyzzzzz[k];

                g_x_y_xz_xzzzzz[k] = -g_x_y_x_xzzzzz[k] * ab_z + g_x_y_x_xzzzzzz[k];

                g_x_y_xz_yyyyyy[k] = -g_x_y_x_yyyyyy[k] * ab_z + g_x_y_x_yyyyyyz[k];

                g_x_y_xz_yyyyyz[k] = -g_x_y_x_yyyyyz[k] * ab_z + g_x_y_x_yyyyyzz[k];

                g_x_y_xz_yyyyzz[k] = -g_x_y_x_yyyyzz[k] * ab_z + g_x_y_x_yyyyzzz[k];

                g_x_y_xz_yyyzzz[k] = -g_x_y_x_yyyzzz[k] * ab_z + g_x_y_x_yyyzzzz[k];

                g_x_y_xz_yyzzzz[k] = -g_x_y_x_yyzzzz[k] * ab_z + g_x_y_x_yyzzzzz[k];

                g_x_y_xz_yzzzzz[k] = -g_x_y_x_yzzzzz[k] * ab_z + g_x_y_x_yzzzzzz[k];

                g_x_y_xz_zzzzzz[k] = -g_x_y_x_zzzzzz[k] * ab_z + g_x_y_x_zzzzzzz[k];
            }

            /// Set up 252-280 components of targeted buffer : cbuffer.data(

            auto g_x_y_yy_xxxxxx = cbuffer.data(di_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_y_yy_xxxxxy = cbuffer.data(di_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_y_yy_xxxxxz = cbuffer.data(di_geom_11_off + 254 * ccomps * dcomps);

            auto g_x_y_yy_xxxxyy = cbuffer.data(di_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_y_yy_xxxxyz = cbuffer.data(di_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_y_yy_xxxxzz = cbuffer.data(di_geom_11_off + 257 * ccomps * dcomps);

            auto g_x_y_yy_xxxyyy = cbuffer.data(di_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_y_yy_xxxyyz = cbuffer.data(di_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_y_yy_xxxyzz = cbuffer.data(di_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_y_yy_xxxzzz = cbuffer.data(di_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_y_yy_xxyyyy = cbuffer.data(di_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_y_yy_xxyyyz = cbuffer.data(di_geom_11_off + 263 * ccomps * dcomps);

            auto g_x_y_yy_xxyyzz = cbuffer.data(di_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_y_yy_xxyzzz = cbuffer.data(di_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_y_yy_xxzzzz = cbuffer.data(di_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_y_yy_xyyyyy = cbuffer.data(di_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_y_yy_xyyyyz = cbuffer.data(di_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_y_yy_xyyyzz = cbuffer.data(di_geom_11_off + 269 * ccomps * dcomps);

            auto g_x_y_yy_xyyzzz = cbuffer.data(di_geom_11_off + 270 * ccomps * dcomps);

            auto g_x_y_yy_xyzzzz = cbuffer.data(di_geom_11_off + 271 * ccomps * dcomps);

            auto g_x_y_yy_xzzzzz = cbuffer.data(di_geom_11_off + 272 * ccomps * dcomps);

            auto g_x_y_yy_yyyyyy = cbuffer.data(di_geom_11_off + 273 * ccomps * dcomps);

            auto g_x_y_yy_yyyyyz = cbuffer.data(di_geom_11_off + 274 * ccomps * dcomps);

            auto g_x_y_yy_yyyyzz = cbuffer.data(di_geom_11_off + 275 * ccomps * dcomps);

            auto g_x_y_yy_yyyzzz = cbuffer.data(di_geom_11_off + 276 * ccomps * dcomps);

            auto g_x_y_yy_yyzzzz = cbuffer.data(di_geom_11_off + 277 * ccomps * dcomps);

            auto g_x_y_yy_yzzzzz = cbuffer.data(di_geom_11_off + 278 * ccomps * dcomps);

            auto g_x_y_yy_zzzzzz = cbuffer.data(di_geom_11_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxxxxx, g_x_0_y_xxxxxy, g_x_0_y_xxxxxz, g_x_0_y_xxxxyy, g_x_0_y_xxxxyz, g_x_0_y_xxxxzz, g_x_0_y_xxxyyy, g_x_0_y_xxxyyz, g_x_0_y_xxxyzz, g_x_0_y_xxxzzz, g_x_0_y_xxyyyy, g_x_0_y_xxyyyz, g_x_0_y_xxyyzz, g_x_0_y_xxyzzz, g_x_0_y_xxzzzz, g_x_0_y_xyyyyy, g_x_0_y_xyyyyz, g_x_0_y_xyyyzz, g_x_0_y_xyyzzz, g_x_0_y_xyzzzz, g_x_0_y_xzzzzz, g_x_0_y_yyyyyy, g_x_0_y_yyyyyz, g_x_0_y_yyyyzz, g_x_0_y_yyyzzz, g_x_0_y_yyzzzz, g_x_0_y_yzzzzz, g_x_0_y_zzzzzz, g_x_y_y_xxxxxx, g_x_y_y_xxxxxxy, g_x_y_y_xxxxxy, g_x_y_y_xxxxxyy, g_x_y_y_xxxxxyz, g_x_y_y_xxxxxz, g_x_y_y_xxxxyy, g_x_y_y_xxxxyyy, g_x_y_y_xxxxyyz, g_x_y_y_xxxxyz, g_x_y_y_xxxxyzz, g_x_y_y_xxxxzz, g_x_y_y_xxxyyy, g_x_y_y_xxxyyyy, g_x_y_y_xxxyyyz, g_x_y_y_xxxyyz, g_x_y_y_xxxyyzz, g_x_y_y_xxxyzz, g_x_y_y_xxxyzzz, g_x_y_y_xxxzzz, g_x_y_y_xxyyyy, g_x_y_y_xxyyyyy, g_x_y_y_xxyyyyz, g_x_y_y_xxyyyz, g_x_y_y_xxyyyzz, g_x_y_y_xxyyzz, g_x_y_y_xxyyzzz, g_x_y_y_xxyzzz, g_x_y_y_xxyzzzz, g_x_y_y_xxzzzz, g_x_y_y_xyyyyy, g_x_y_y_xyyyyyy, g_x_y_y_xyyyyyz, g_x_y_y_xyyyyz, g_x_y_y_xyyyyzz, g_x_y_y_xyyyzz, g_x_y_y_xyyyzzz, g_x_y_y_xyyzzz, g_x_y_y_xyyzzzz, g_x_y_y_xyzzzz, g_x_y_y_xyzzzzz, g_x_y_y_xzzzzz, g_x_y_y_yyyyyy, g_x_y_y_yyyyyyy, g_x_y_y_yyyyyyz, g_x_y_y_yyyyyz, g_x_y_y_yyyyyzz, g_x_y_y_yyyyzz, g_x_y_y_yyyyzzz, g_x_y_y_yyyzzz, g_x_y_y_yyyzzzz, g_x_y_y_yyzzzz, g_x_y_y_yyzzzzz, g_x_y_y_yzzzzz, g_x_y_y_yzzzzzz, g_x_y_y_zzzzzz, g_x_y_yy_xxxxxx, g_x_y_yy_xxxxxy, g_x_y_yy_xxxxxz, g_x_y_yy_xxxxyy, g_x_y_yy_xxxxyz, g_x_y_yy_xxxxzz, g_x_y_yy_xxxyyy, g_x_y_yy_xxxyyz, g_x_y_yy_xxxyzz, g_x_y_yy_xxxzzz, g_x_y_yy_xxyyyy, g_x_y_yy_xxyyyz, g_x_y_yy_xxyyzz, g_x_y_yy_xxyzzz, g_x_y_yy_xxzzzz, g_x_y_yy_xyyyyy, g_x_y_yy_xyyyyz, g_x_y_yy_xyyyzz, g_x_y_yy_xyyzzz, g_x_y_yy_xyzzzz, g_x_y_yy_xzzzzz, g_x_y_yy_yyyyyy, g_x_y_yy_yyyyyz, g_x_y_yy_yyyyzz, g_x_y_yy_yyyzzz, g_x_y_yy_yyzzzz, g_x_y_yy_yzzzzz, g_x_y_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yy_xxxxxx[k] = g_x_0_y_xxxxxx[k] - g_x_y_y_xxxxxx[k] * ab_y + g_x_y_y_xxxxxxy[k];

                g_x_y_yy_xxxxxy[k] = g_x_0_y_xxxxxy[k] - g_x_y_y_xxxxxy[k] * ab_y + g_x_y_y_xxxxxyy[k];

                g_x_y_yy_xxxxxz[k] = g_x_0_y_xxxxxz[k] - g_x_y_y_xxxxxz[k] * ab_y + g_x_y_y_xxxxxyz[k];

                g_x_y_yy_xxxxyy[k] = g_x_0_y_xxxxyy[k] - g_x_y_y_xxxxyy[k] * ab_y + g_x_y_y_xxxxyyy[k];

                g_x_y_yy_xxxxyz[k] = g_x_0_y_xxxxyz[k] - g_x_y_y_xxxxyz[k] * ab_y + g_x_y_y_xxxxyyz[k];

                g_x_y_yy_xxxxzz[k] = g_x_0_y_xxxxzz[k] - g_x_y_y_xxxxzz[k] * ab_y + g_x_y_y_xxxxyzz[k];

                g_x_y_yy_xxxyyy[k] = g_x_0_y_xxxyyy[k] - g_x_y_y_xxxyyy[k] * ab_y + g_x_y_y_xxxyyyy[k];

                g_x_y_yy_xxxyyz[k] = g_x_0_y_xxxyyz[k] - g_x_y_y_xxxyyz[k] * ab_y + g_x_y_y_xxxyyyz[k];

                g_x_y_yy_xxxyzz[k] = g_x_0_y_xxxyzz[k] - g_x_y_y_xxxyzz[k] * ab_y + g_x_y_y_xxxyyzz[k];

                g_x_y_yy_xxxzzz[k] = g_x_0_y_xxxzzz[k] - g_x_y_y_xxxzzz[k] * ab_y + g_x_y_y_xxxyzzz[k];

                g_x_y_yy_xxyyyy[k] = g_x_0_y_xxyyyy[k] - g_x_y_y_xxyyyy[k] * ab_y + g_x_y_y_xxyyyyy[k];

                g_x_y_yy_xxyyyz[k] = g_x_0_y_xxyyyz[k] - g_x_y_y_xxyyyz[k] * ab_y + g_x_y_y_xxyyyyz[k];

                g_x_y_yy_xxyyzz[k] = g_x_0_y_xxyyzz[k] - g_x_y_y_xxyyzz[k] * ab_y + g_x_y_y_xxyyyzz[k];

                g_x_y_yy_xxyzzz[k] = g_x_0_y_xxyzzz[k] - g_x_y_y_xxyzzz[k] * ab_y + g_x_y_y_xxyyzzz[k];

                g_x_y_yy_xxzzzz[k] = g_x_0_y_xxzzzz[k] - g_x_y_y_xxzzzz[k] * ab_y + g_x_y_y_xxyzzzz[k];

                g_x_y_yy_xyyyyy[k] = g_x_0_y_xyyyyy[k] - g_x_y_y_xyyyyy[k] * ab_y + g_x_y_y_xyyyyyy[k];

                g_x_y_yy_xyyyyz[k] = g_x_0_y_xyyyyz[k] - g_x_y_y_xyyyyz[k] * ab_y + g_x_y_y_xyyyyyz[k];

                g_x_y_yy_xyyyzz[k] = g_x_0_y_xyyyzz[k] - g_x_y_y_xyyyzz[k] * ab_y + g_x_y_y_xyyyyzz[k];

                g_x_y_yy_xyyzzz[k] = g_x_0_y_xyyzzz[k] - g_x_y_y_xyyzzz[k] * ab_y + g_x_y_y_xyyyzzz[k];

                g_x_y_yy_xyzzzz[k] = g_x_0_y_xyzzzz[k] - g_x_y_y_xyzzzz[k] * ab_y + g_x_y_y_xyyzzzz[k];

                g_x_y_yy_xzzzzz[k] = g_x_0_y_xzzzzz[k] - g_x_y_y_xzzzzz[k] * ab_y + g_x_y_y_xyzzzzz[k];

                g_x_y_yy_yyyyyy[k] = g_x_0_y_yyyyyy[k] - g_x_y_y_yyyyyy[k] * ab_y + g_x_y_y_yyyyyyy[k];

                g_x_y_yy_yyyyyz[k] = g_x_0_y_yyyyyz[k] - g_x_y_y_yyyyyz[k] * ab_y + g_x_y_y_yyyyyyz[k];

                g_x_y_yy_yyyyzz[k] = g_x_0_y_yyyyzz[k] - g_x_y_y_yyyyzz[k] * ab_y + g_x_y_y_yyyyyzz[k];

                g_x_y_yy_yyyzzz[k] = g_x_0_y_yyyzzz[k] - g_x_y_y_yyyzzz[k] * ab_y + g_x_y_y_yyyyzzz[k];

                g_x_y_yy_yyzzzz[k] = g_x_0_y_yyzzzz[k] - g_x_y_y_yyzzzz[k] * ab_y + g_x_y_y_yyyzzzz[k];

                g_x_y_yy_yzzzzz[k] = g_x_0_y_yzzzzz[k] - g_x_y_y_yzzzzz[k] * ab_y + g_x_y_y_yyzzzzz[k];

                g_x_y_yy_zzzzzz[k] = g_x_0_y_zzzzzz[k] - g_x_y_y_zzzzzz[k] * ab_y + g_x_y_y_yzzzzzz[k];
            }

            /// Set up 280-308 components of targeted buffer : cbuffer.data(

            auto g_x_y_yz_xxxxxx = cbuffer.data(di_geom_11_off + 280 * ccomps * dcomps);

            auto g_x_y_yz_xxxxxy = cbuffer.data(di_geom_11_off + 281 * ccomps * dcomps);

            auto g_x_y_yz_xxxxxz = cbuffer.data(di_geom_11_off + 282 * ccomps * dcomps);

            auto g_x_y_yz_xxxxyy = cbuffer.data(di_geom_11_off + 283 * ccomps * dcomps);

            auto g_x_y_yz_xxxxyz = cbuffer.data(di_geom_11_off + 284 * ccomps * dcomps);

            auto g_x_y_yz_xxxxzz = cbuffer.data(di_geom_11_off + 285 * ccomps * dcomps);

            auto g_x_y_yz_xxxyyy = cbuffer.data(di_geom_11_off + 286 * ccomps * dcomps);

            auto g_x_y_yz_xxxyyz = cbuffer.data(di_geom_11_off + 287 * ccomps * dcomps);

            auto g_x_y_yz_xxxyzz = cbuffer.data(di_geom_11_off + 288 * ccomps * dcomps);

            auto g_x_y_yz_xxxzzz = cbuffer.data(di_geom_11_off + 289 * ccomps * dcomps);

            auto g_x_y_yz_xxyyyy = cbuffer.data(di_geom_11_off + 290 * ccomps * dcomps);

            auto g_x_y_yz_xxyyyz = cbuffer.data(di_geom_11_off + 291 * ccomps * dcomps);

            auto g_x_y_yz_xxyyzz = cbuffer.data(di_geom_11_off + 292 * ccomps * dcomps);

            auto g_x_y_yz_xxyzzz = cbuffer.data(di_geom_11_off + 293 * ccomps * dcomps);

            auto g_x_y_yz_xxzzzz = cbuffer.data(di_geom_11_off + 294 * ccomps * dcomps);

            auto g_x_y_yz_xyyyyy = cbuffer.data(di_geom_11_off + 295 * ccomps * dcomps);

            auto g_x_y_yz_xyyyyz = cbuffer.data(di_geom_11_off + 296 * ccomps * dcomps);

            auto g_x_y_yz_xyyyzz = cbuffer.data(di_geom_11_off + 297 * ccomps * dcomps);

            auto g_x_y_yz_xyyzzz = cbuffer.data(di_geom_11_off + 298 * ccomps * dcomps);

            auto g_x_y_yz_xyzzzz = cbuffer.data(di_geom_11_off + 299 * ccomps * dcomps);

            auto g_x_y_yz_xzzzzz = cbuffer.data(di_geom_11_off + 300 * ccomps * dcomps);

            auto g_x_y_yz_yyyyyy = cbuffer.data(di_geom_11_off + 301 * ccomps * dcomps);

            auto g_x_y_yz_yyyyyz = cbuffer.data(di_geom_11_off + 302 * ccomps * dcomps);

            auto g_x_y_yz_yyyyzz = cbuffer.data(di_geom_11_off + 303 * ccomps * dcomps);

            auto g_x_y_yz_yyyzzz = cbuffer.data(di_geom_11_off + 304 * ccomps * dcomps);

            auto g_x_y_yz_yyzzzz = cbuffer.data(di_geom_11_off + 305 * ccomps * dcomps);

            auto g_x_y_yz_yzzzzz = cbuffer.data(di_geom_11_off + 306 * ccomps * dcomps);

            auto g_x_y_yz_zzzzzz = cbuffer.data(di_geom_11_off + 307 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_y_xxxxxx, g_x_y_y_xxxxxxz, g_x_y_y_xxxxxy, g_x_y_y_xxxxxyz, g_x_y_y_xxxxxz, g_x_y_y_xxxxxzz, g_x_y_y_xxxxyy, g_x_y_y_xxxxyyz, g_x_y_y_xxxxyz, g_x_y_y_xxxxyzz, g_x_y_y_xxxxzz, g_x_y_y_xxxxzzz, g_x_y_y_xxxyyy, g_x_y_y_xxxyyyz, g_x_y_y_xxxyyz, g_x_y_y_xxxyyzz, g_x_y_y_xxxyzz, g_x_y_y_xxxyzzz, g_x_y_y_xxxzzz, g_x_y_y_xxxzzzz, g_x_y_y_xxyyyy, g_x_y_y_xxyyyyz, g_x_y_y_xxyyyz, g_x_y_y_xxyyyzz, g_x_y_y_xxyyzz, g_x_y_y_xxyyzzz, g_x_y_y_xxyzzz, g_x_y_y_xxyzzzz, g_x_y_y_xxzzzz, g_x_y_y_xxzzzzz, g_x_y_y_xyyyyy, g_x_y_y_xyyyyyz, g_x_y_y_xyyyyz, g_x_y_y_xyyyyzz, g_x_y_y_xyyyzz, g_x_y_y_xyyyzzz, g_x_y_y_xyyzzz, g_x_y_y_xyyzzzz, g_x_y_y_xyzzzz, g_x_y_y_xyzzzzz, g_x_y_y_xzzzzz, g_x_y_y_xzzzzzz, g_x_y_y_yyyyyy, g_x_y_y_yyyyyyz, g_x_y_y_yyyyyz, g_x_y_y_yyyyyzz, g_x_y_y_yyyyzz, g_x_y_y_yyyyzzz, g_x_y_y_yyyzzz, g_x_y_y_yyyzzzz, g_x_y_y_yyzzzz, g_x_y_y_yyzzzzz, g_x_y_y_yzzzzz, g_x_y_y_yzzzzzz, g_x_y_y_zzzzzz, g_x_y_y_zzzzzzz, g_x_y_yz_xxxxxx, g_x_y_yz_xxxxxy, g_x_y_yz_xxxxxz, g_x_y_yz_xxxxyy, g_x_y_yz_xxxxyz, g_x_y_yz_xxxxzz, g_x_y_yz_xxxyyy, g_x_y_yz_xxxyyz, g_x_y_yz_xxxyzz, g_x_y_yz_xxxzzz, g_x_y_yz_xxyyyy, g_x_y_yz_xxyyyz, g_x_y_yz_xxyyzz, g_x_y_yz_xxyzzz, g_x_y_yz_xxzzzz, g_x_y_yz_xyyyyy, g_x_y_yz_xyyyyz, g_x_y_yz_xyyyzz, g_x_y_yz_xyyzzz, g_x_y_yz_xyzzzz, g_x_y_yz_xzzzzz, g_x_y_yz_yyyyyy, g_x_y_yz_yyyyyz, g_x_y_yz_yyyyzz, g_x_y_yz_yyyzzz, g_x_y_yz_yyzzzz, g_x_y_yz_yzzzzz, g_x_y_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yz_xxxxxx[k] = -g_x_y_y_xxxxxx[k] * ab_z + g_x_y_y_xxxxxxz[k];

                g_x_y_yz_xxxxxy[k] = -g_x_y_y_xxxxxy[k] * ab_z + g_x_y_y_xxxxxyz[k];

                g_x_y_yz_xxxxxz[k] = -g_x_y_y_xxxxxz[k] * ab_z + g_x_y_y_xxxxxzz[k];

                g_x_y_yz_xxxxyy[k] = -g_x_y_y_xxxxyy[k] * ab_z + g_x_y_y_xxxxyyz[k];

                g_x_y_yz_xxxxyz[k] = -g_x_y_y_xxxxyz[k] * ab_z + g_x_y_y_xxxxyzz[k];

                g_x_y_yz_xxxxzz[k] = -g_x_y_y_xxxxzz[k] * ab_z + g_x_y_y_xxxxzzz[k];

                g_x_y_yz_xxxyyy[k] = -g_x_y_y_xxxyyy[k] * ab_z + g_x_y_y_xxxyyyz[k];

                g_x_y_yz_xxxyyz[k] = -g_x_y_y_xxxyyz[k] * ab_z + g_x_y_y_xxxyyzz[k];

                g_x_y_yz_xxxyzz[k] = -g_x_y_y_xxxyzz[k] * ab_z + g_x_y_y_xxxyzzz[k];

                g_x_y_yz_xxxzzz[k] = -g_x_y_y_xxxzzz[k] * ab_z + g_x_y_y_xxxzzzz[k];

                g_x_y_yz_xxyyyy[k] = -g_x_y_y_xxyyyy[k] * ab_z + g_x_y_y_xxyyyyz[k];

                g_x_y_yz_xxyyyz[k] = -g_x_y_y_xxyyyz[k] * ab_z + g_x_y_y_xxyyyzz[k];

                g_x_y_yz_xxyyzz[k] = -g_x_y_y_xxyyzz[k] * ab_z + g_x_y_y_xxyyzzz[k];

                g_x_y_yz_xxyzzz[k] = -g_x_y_y_xxyzzz[k] * ab_z + g_x_y_y_xxyzzzz[k];

                g_x_y_yz_xxzzzz[k] = -g_x_y_y_xxzzzz[k] * ab_z + g_x_y_y_xxzzzzz[k];

                g_x_y_yz_xyyyyy[k] = -g_x_y_y_xyyyyy[k] * ab_z + g_x_y_y_xyyyyyz[k];

                g_x_y_yz_xyyyyz[k] = -g_x_y_y_xyyyyz[k] * ab_z + g_x_y_y_xyyyyzz[k];

                g_x_y_yz_xyyyzz[k] = -g_x_y_y_xyyyzz[k] * ab_z + g_x_y_y_xyyyzzz[k];

                g_x_y_yz_xyyzzz[k] = -g_x_y_y_xyyzzz[k] * ab_z + g_x_y_y_xyyzzzz[k];

                g_x_y_yz_xyzzzz[k] = -g_x_y_y_xyzzzz[k] * ab_z + g_x_y_y_xyzzzzz[k];

                g_x_y_yz_xzzzzz[k] = -g_x_y_y_xzzzzz[k] * ab_z + g_x_y_y_xzzzzzz[k];

                g_x_y_yz_yyyyyy[k] = -g_x_y_y_yyyyyy[k] * ab_z + g_x_y_y_yyyyyyz[k];

                g_x_y_yz_yyyyyz[k] = -g_x_y_y_yyyyyz[k] * ab_z + g_x_y_y_yyyyyzz[k];

                g_x_y_yz_yyyyzz[k] = -g_x_y_y_yyyyzz[k] * ab_z + g_x_y_y_yyyyzzz[k];

                g_x_y_yz_yyyzzz[k] = -g_x_y_y_yyyzzz[k] * ab_z + g_x_y_y_yyyzzzz[k];

                g_x_y_yz_yyzzzz[k] = -g_x_y_y_yyzzzz[k] * ab_z + g_x_y_y_yyzzzzz[k];

                g_x_y_yz_yzzzzz[k] = -g_x_y_y_yzzzzz[k] * ab_z + g_x_y_y_yzzzzzz[k];

                g_x_y_yz_zzzzzz[k] = -g_x_y_y_zzzzzz[k] * ab_z + g_x_y_y_zzzzzzz[k];
            }

            /// Set up 308-336 components of targeted buffer : cbuffer.data(

            auto g_x_y_zz_xxxxxx = cbuffer.data(di_geom_11_off + 308 * ccomps * dcomps);

            auto g_x_y_zz_xxxxxy = cbuffer.data(di_geom_11_off + 309 * ccomps * dcomps);

            auto g_x_y_zz_xxxxxz = cbuffer.data(di_geom_11_off + 310 * ccomps * dcomps);

            auto g_x_y_zz_xxxxyy = cbuffer.data(di_geom_11_off + 311 * ccomps * dcomps);

            auto g_x_y_zz_xxxxyz = cbuffer.data(di_geom_11_off + 312 * ccomps * dcomps);

            auto g_x_y_zz_xxxxzz = cbuffer.data(di_geom_11_off + 313 * ccomps * dcomps);

            auto g_x_y_zz_xxxyyy = cbuffer.data(di_geom_11_off + 314 * ccomps * dcomps);

            auto g_x_y_zz_xxxyyz = cbuffer.data(di_geom_11_off + 315 * ccomps * dcomps);

            auto g_x_y_zz_xxxyzz = cbuffer.data(di_geom_11_off + 316 * ccomps * dcomps);

            auto g_x_y_zz_xxxzzz = cbuffer.data(di_geom_11_off + 317 * ccomps * dcomps);

            auto g_x_y_zz_xxyyyy = cbuffer.data(di_geom_11_off + 318 * ccomps * dcomps);

            auto g_x_y_zz_xxyyyz = cbuffer.data(di_geom_11_off + 319 * ccomps * dcomps);

            auto g_x_y_zz_xxyyzz = cbuffer.data(di_geom_11_off + 320 * ccomps * dcomps);

            auto g_x_y_zz_xxyzzz = cbuffer.data(di_geom_11_off + 321 * ccomps * dcomps);

            auto g_x_y_zz_xxzzzz = cbuffer.data(di_geom_11_off + 322 * ccomps * dcomps);

            auto g_x_y_zz_xyyyyy = cbuffer.data(di_geom_11_off + 323 * ccomps * dcomps);

            auto g_x_y_zz_xyyyyz = cbuffer.data(di_geom_11_off + 324 * ccomps * dcomps);

            auto g_x_y_zz_xyyyzz = cbuffer.data(di_geom_11_off + 325 * ccomps * dcomps);

            auto g_x_y_zz_xyyzzz = cbuffer.data(di_geom_11_off + 326 * ccomps * dcomps);

            auto g_x_y_zz_xyzzzz = cbuffer.data(di_geom_11_off + 327 * ccomps * dcomps);

            auto g_x_y_zz_xzzzzz = cbuffer.data(di_geom_11_off + 328 * ccomps * dcomps);

            auto g_x_y_zz_yyyyyy = cbuffer.data(di_geom_11_off + 329 * ccomps * dcomps);

            auto g_x_y_zz_yyyyyz = cbuffer.data(di_geom_11_off + 330 * ccomps * dcomps);

            auto g_x_y_zz_yyyyzz = cbuffer.data(di_geom_11_off + 331 * ccomps * dcomps);

            auto g_x_y_zz_yyyzzz = cbuffer.data(di_geom_11_off + 332 * ccomps * dcomps);

            auto g_x_y_zz_yyzzzz = cbuffer.data(di_geom_11_off + 333 * ccomps * dcomps);

            auto g_x_y_zz_yzzzzz = cbuffer.data(di_geom_11_off + 334 * ccomps * dcomps);

            auto g_x_y_zz_zzzzzz = cbuffer.data(di_geom_11_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_z_xxxxxx, g_x_y_z_xxxxxxz, g_x_y_z_xxxxxy, g_x_y_z_xxxxxyz, g_x_y_z_xxxxxz, g_x_y_z_xxxxxzz, g_x_y_z_xxxxyy, g_x_y_z_xxxxyyz, g_x_y_z_xxxxyz, g_x_y_z_xxxxyzz, g_x_y_z_xxxxzz, g_x_y_z_xxxxzzz, g_x_y_z_xxxyyy, g_x_y_z_xxxyyyz, g_x_y_z_xxxyyz, g_x_y_z_xxxyyzz, g_x_y_z_xxxyzz, g_x_y_z_xxxyzzz, g_x_y_z_xxxzzz, g_x_y_z_xxxzzzz, g_x_y_z_xxyyyy, g_x_y_z_xxyyyyz, g_x_y_z_xxyyyz, g_x_y_z_xxyyyzz, g_x_y_z_xxyyzz, g_x_y_z_xxyyzzz, g_x_y_z_xxyzzz, g_x_y_z_xxyzzzz, g_x_y_z_xxzzzz, g_x_y_z_xxzzzzz, g_x_y_z_xyyyyy, g_x_y_z_xyyyyyz, g_x_y_z_xyyyyz, g_x_y_z_xyyyyzz, g_x_y_z_xyyyzz, g_x_y_z_xyyyzzz, g_x_y_z_xyyzzz, g_x_y_z_xyyzzzz, g_x_y_z_xyzzzz, g_x_y_z_xyzzzzz, g_x_y_z_xzzzzz, g_x_y_z_xzzzzzz, g_x_y_z_yyyyyy, g_x_y_z_yyyyyyz, g_x_y_z_yyyyyz, g_x_y_z_yyyyyzz, g_x_y_z_yyyyzz, g_x_y_z_yyyyzzz, g_x_y_z_yyyzzz, g_x_y_z_yyyzzzz, g_x_y_z_yyzzzz, g_x_y_z_yyzzzzz, g_x_y_z_yzzzzz, g_x_y_z_yzzzzzz, g_x_y_z_zzzzzz, g_x_y_z_zzzzzzz, g_x_y_zz_xxxxxx, g_x_y_zz_xxxxxy, g_x_y_zz_xxxxxz, g_x_y_zz_xxxxyy, g_x_y_zz_xxxxyz, g_x_y_zz_xxxxzz, g_x_y_zz_xxxyyy, g_x_y_zz_xxxyyz, g_x_y_zz_xxxyzz, g_x_y_zz_xxxzzz, g_x_y_zz_xxyyyy, g_x_y_zz_xxyyyz, g_x_y_zz_xxyyzz, g_x_y_zz_xxyzzz, g_x_y_zz_xxzzzz, g_x_y_zz_xyyyyy, g_x_y_zz_xyyyyz, g_x_y_zz_xyyyzz, g_x_y_zz_xyyzzz, g_x_y_zz_xyzzzz, g_x_y_zz_xzzzzz, g_x_y_zz_yyyyyy, g_x_y_zz_yyyyyz, g_x_y_zz_yyyyzz, g_x_y_zz_yyyzzz, g_x_y_zz_yyzzzz, g_x_y_zz_yzzzzz, g_x_y_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zz_xxxxxx[k] = -g_x_y_z_xxxxxx[k] * ab_z + g_x_y_z_xxxxxxz[k];

                g_x_y_zz_xxxxxy[k] = -g_x_y_z_xxxxxy[k] * ab_z + g_x_y_z_xxxxxyz[k];

                g_x_y_zz_xxxxxz[k] = -g_x_y_z_xxxxxz[k] * ab_z + g_x_y_z_xxxxxzz[k];

                g_x_y_zz_xxxxyy[k] = -g_x_y_z_xxxxyy[k] * ab_z + g_x_y_z_xxxxyyz[k];

                g_x_y_zz_xxxxyz[k] = -g_x_y_z_xxxxyz[k] * ab_z + g_x_y_z_xxxxyzz[k];

                g_x_y_zz_xxxxzz[k] = -g_x_y_z_xxxxzz[k] * ab_z + g_x_y_z_xxxxzzz[k];

                g_x_y_zz_xxxyyy[k] = -g_x_y_z_xxxyyy[k] * ab_z + g_x_y_z_xxxyyyz[k];

                g_x_y_zz_xxxyyz[k] = -g_x_y_z_xxxyyz[k] * ab_z + g_x_y_z_xxxyyzz[k];

                g_x_y_zz_xxxyzz[k] = -g_x_y_z_xxxyzz[k] * ab_z + g_x_y_z_xxxyzzz[k];

                g_x_y_zz_xxxzzz[k] = -g_x_y_z_xxxzzz[k] * ab_z + g_x_y_z_xxxzzzz[k];

                g_x_y_zz_xxyyyy[k] = -g_x_y_z_xxyyyy[k] * ab_z + g_x_y_z_xxyyyyz[k];

                g_x_y_zz_xxyyyz[k] = -g_x_y_z_xxyyyz[k] * ab_z + g_x_y_z_xxyyyzz[k];

                g_x_y_zz_xxyyzz[k] = -g_x_y_z_xxyyzz[k] * ab_z + g_x_y_z_xxyyzzz[k];

                g_x_y_zz_xxyzzz[k] = -g_x_y_z_xxyzzz[k] * ab_z + g_x_y_z_xxyzzzz[k];

                g_x_y_zz_xxzzzz[k] = -g_x_y_z_xxzzzz[k] * ab_z + g_x_y_z_xxzzzzz[k];

                g_x_y_zz_xyyyyy[k] = -g_x_y_z_xyyyyy[k] * ab_z + g_x_y_z_xyyyyyz[k];

                g_x_y_zz_xyyyyz[k] = -g_x_y_z_xyyyyz[k] * ab_z + g_x_y_z_xyyyyzz[k];

                g_x_y_zz_xyyyzz[k] = -g_x_y_z_xyyyzz[k] * ab_z + g_x_y_z_xyyyzzz[k];

                g_x_y_zz_xyyzzz[k] = -g_x_y_z_xyyzzz[k] * ab_z + g_x_y_z_xyyzzzz[k];

                g_x_y_zz_xyzzzz[k] = -g_x_y_z_xyzzzz[k] * ab_z + g_x_y_z_xyzzzzz[k];

                g_x_y_zz_xzzzzz[k] = -g_x_y_z_xzzzzz[k] * ab_z + g_x_y_z_xzzzzzz[k];

                g_x_y_zz_yyyyyy[k] = -g_x_y_z_yyyyyy[k] * ab_z + g_x_y_z_yyyyyyz[k];

                g_x_y_zz_yyyyyz[k] = -g_x_y_z_yyyyyz[k] * ab_z + g_x_y_z_yyyyyzz[k];

                g_x_y_zz_yyyyzz[k] = -g_x_y_z_yyyyzz[k] * ab_z + g_x_y_z_yyyyzzz[k];

                g_x_y_zz_yyyzzz[k] = -g_x_y_z_yyyzzz[k] * ab_z + g_x_y_z_yyyzzzz[k];

                g_x_y_zz_yyzzzz[k] = -g_x_y_z_yyzzzz[k] * ab_z + g_x_y_z_yyzzzzz[k];

                g_x_y_zz_yzzzzz[k] = -g_x_y_z_yzzzzz[k] * ab_z + g_x_y_z_yzzzzzz[k];

                g_x_y_zz_zzzzzz[k] = -g_x_y_z_zzzzzz[k] * ab_z + g_x_y_z_zzzzzzz[k];
            }

            /// Set up 336-364 components of targeted buffer : cbuffer.data(

            auto g_x_z_xx_xxxxxx = cbuffer.data(di_geom_11_off + 336 * ccomps * dcomps);

            auto g_x_z_xx_xxxxxy = cbuffer.data(di_geom_11_off + 337 * ccomps * dcomps);

            auto g_x_z_xx_xxxxxz = cbuffer.data(di_geom_11_off + 338 * ccomps * dcomps);

            auto g_x_z_xx_xxxxyy = cbuffer.data(di_geom_11_off + 339 * ccomps * dcomps);

            auto g_x_z_xx_xxxxyz = cbuffer.data(di_geom_11_off + 340 * ccomps * dcomps);

            auto g_x_z_xx_xxxxzz = cbuffer.data(di_geom_11_off + 341 * ccomps * dcomps);

            auto g_x_z_xx_xxxyyy = cbuffer.data(di_geom_11_off + 342 * ccomps * dcomps);

            auto g_x_z_xx_xxxyyz = cbuffer.data(di_geom_11_off + 343 * ccomps * dcomps);

            auto g_x_z_xx_xxxyzz = cbuffer.data(di_geom_11_off + 344 * ccomps * dcomps);

            auto g_x_z_xx_xxxzzz = cbuffer.data(di_geom_11_off + 345 * ccomps * dcomps);

            auto g_x_z_xx_xxyyyy = cbuffer.data(di_geom_11_off + 346 * ccomps * dcomps);

            auto g_x_z_xx_xxyyyz = cbuffer.data(di_geom_11_off + 347 * ccomps * dcomps);

            auto g_x_z_xx_xxyyzz = cbuffer.data(di_geom_11_off + 348 * ccomps * dcomps);

            auto g_x_z_xx_xxyzzz = cbuffer.data(di_geom_11_off + 349 * ccomps * dcomps);

            auto g_x_z_xx_xxzzzz = cbuffer.data(di_geom_11_off + 350 * ccomps * dcomps);

            auto g_x_z_xx_xyyyyy = cbuffer.data(di_geom_11_off + 351 * ccomps * dcomps);

            auto g_x_z_xx_xyyyyz = cbuffer.data(di_geom_11_off + 352 * ccomps * dcomps);

            auto g_x_z_xx_xyyyzz = cbuffer.data(di_geom_11_off + 353 * ccomps * dcomps);

            auto g_x_z_xx_xyyzzz = cbuffer.data(di_geom_11_off + 354 * ccomps * dcomps);

            auto g_x_z_xx_xyzzzz = cbuffer.data(di_geom_11_off + 355 * ccomps * dcomps);

            auto g_x_z_xx_xzzzzz = cbuffer.data(di_geom_11_off + 356 * ccomps * dcomps);

            auto g_x_z_xx_yyyyyy = cbuffer.data(di_geom_11_off + 357 * ccomps * dcomps);

            auto g_x_z_xx_yyyyyz = cbuffer.data(di_geom_11_off + 358 * ccomps * dcomps);

            auto g_x_z_xx_yyyyzz = cbuffer.data(di_geom_11_off + 359 * ccomps * dcomps);

            auto g_x_z_xx_yyyzzz = cbuffer.data(di_geom_11_off + 360 * ccomps * dcomps);

            auto g_x_z_xx_yyzzzz = cbuffer.data(di_geom_11_off + 361 * ccomps * dcomps);

            auto g_x_z_xx_yzzzzz = cbuffer.data(di_geom_11_off + 362 * ccomps * dcomps);

            auto g_x_z_xx_zzzzzz = cbuffer.data(di_geom_11_off + 363 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_x_xxxxxx, g_0_z_x_xxxxxy, g_0_z_x_xxxxxz, g_0_z_x_xxxxyy, g_0_z_x_xxxxyz, g_0_z_x_xxxxzz, g_0_z_x_xxxyyy, g_0_z_x_xxxyyz, g_0_z_x_xxxyzz, g_0_z_x_xxxzzz, g_0_z_x_xxyyyy, g_0_z_x_xxyyyz, g_0_z_x_xxyyzz, g_0_z_x_xxyzzz, g_0_z_x_xxzzzz, g_0_z_x_xyyyyy, g_0_z_x_xyyyyz, g_0_z_x_xyyyzz, g_0_z_x_xyyzzz, g_0_z_x_xyzzzz, g_0_z_x_xzzzzz, g_0_z_x_yyyyyy, g_0_z_x_yyyyyz, g_0_z_x_yyyyzz, g_0_z_x_yyyzzz, g_0_z_x_yyzzzz, g_0_z_x_yzzzzz, g_0_z_x_zzzzzz, g_x_z_x_xxxxxx, g_x_z_x_xxxxxxx, g_x_z_x_xxxxxxy, g_x_z_x_xxxxxxz, g_x_z_x_xxxxxy, g_x_z_x_xxxxxyy, g_x_z_x_xxxxxyz, g_x_z_x_xxxxxz, g_x_z_x_xxxxxzz, g_x_z_x_xxxxyy, g_x_z_x_xxxxyyy, g_x_z_x_xxxxyyz, g_x_z_x_xxxxyz, g_x_z_x_xxxxyzz, g_x_z_x_xxxxzz, g_x_z_x_xxxxzzz, g_x_z_x_xxxyyy, g_x_z_x_xxxyyyy, g_x_z_x_xxxyyyz, g_x_z_x_xxxyyz, g_x_z_x_xxxyyzz, g_x_z_x_xxxyzz, g_x_z_x_xxxyzzz, g_x_z_x_xxxzzz, g_x_z_x_xxxzzzz, g_x_z_x_xxyyyy, g_x_z_x_xxyyyyy, g_x_z_x_xxyyyyz, g_x_z_x_xxyyyz, g_x_z_x_xxyyyzz, g_x_z_x_xxyyzz, g_x_z_x_xxyyzzz, g_x_z_x_xxyzzz, g_x_z_x_xxyzzzz, g_x_z_x_xxzzzz, g_x_z_x_xxzzzzz, g_x_z_x_xyyyyy, g_x_z_x_xyyyyyy, g_x_z_x_xyyyyyz, g_x_z_x_xyyyyz, g_x_z_x_xyyyyzz, g_x_z_x_xyyyzz, g_x_z_x_xyyyzzz, g_x_z_x_xyyzzz, g_x_z_x_xyyzzzz, g_x_z_x_xyzzzz, g_x_z_x_xyzzzzz, g_x_z_x_xzzzzz, g_x_z_x_xzzzzzz, g_x_z_x_yyyyyy, g_x_z_x_yyyyyz, g_x_z_x_yyyyzz, g_x_z_x_yyyzzz, g_x_z_x_yyzzzz, g_x_z_x_yzzzzz, g_x_z_x_zzzzzz, g_x_z_xx_xxxxxx, g_x_z_xx_xxxxxy, g_x_z_xx_xxxxxz, g_x_z_xx_xxxxyy, g_x_z_xx_xxxxyz, g_x_z_xx_xxxxzz, g_x_z_xx_xxxyyy, g_x_z_xx_xxxyyz, g_x_z_xx_xxxyzz, g_x_z_xx_xxxzzz, g_x_z_xx_xxyyyy, g_x_z_xx_xxyyyz, g_x_z_xx_xxyyzz, g_x_z_xx_xxyzzz, g_x_z_xx_xxzzzz, g_x_z_xx_xyyyyy, g_x_z_xx_xyyyyz, g_x_z_xx_xyyyzz, g_x_z_xx_xyyzzz, g_x_z_xx_xyzzzz, g_x_z_xx_xzzzzz, g_x_z_xx_yyyyyy, g_x_z_xx_yyyyyz, g_x_z_xx_yyyyzz, g_x_z_xx_yyyzzz, g_x_z_xx_yyzzzz, g_x_z_xx_yzzzzz, g_x_z_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xx_xxxxxx[k] = -g_0_z_x_xxxxxx[k] - g_x_z_x_xxxxxx[k] * ab_x + g_x_z_x_xxxxxxx[k];

                g_x_z_xx_xxxxxy[k] = -g_0_z_x_xxxxxy[k] - g_x_z_x_xxxxxy[k] * ab_x + g_x_z_x_xxxxxxy[k];

                g_x_z_xx_xxxxxz[k] = -g_0_z_x_xxxxxz[k] - g_x_z_x_xxxxxz[k] * ab_x + g_x_z_x_xxxxxxz[k];

                g_x_z_xx_xxxxyy[k] = -g_0_z_x_xxxxyy[k] - g_x_z_x_xxxxyy[k] * ab_x + g_x_z_x_xxxxxyy[k];

                g_x_z_xx_xxxxyz[k] = -g_0_z_x_xxxxyz[k] - g_x_z_x_xxxxyz[k] * ab_x + g_x_z_x_xxxxxyz[k];

                g_x_z_xx_xxxxzz[k] = -g_0_z_x_xxxxzz[k] - g_x_z_x_xxxxzz[k] * ab_x + g_x_z_x_xxxxxzz[k];

                g_x_z_xx_xxxyyy[k] = -g_0_z_x_xxxyyy[k] - g_x_z_x_xxxyyy[k] * ab_x + g_x_z_x_xxxxyyy[k];

                g_x_z_xx_xxxyyz[k] = -g_0_z_x_xxxyyz[k] - g_x_z_x_xxxyyz[k] * ab_x + g_x_z_x_xxxxyyz[k];

                g_x_z_xx_xxxyzz[k] = -g_0_z_x_xxxyzz[k] - g_x_z_x_xxxyzz[k] * ab_x + g_x_z_x_xxxxyzz[k];

                g_x_z_xx_xxxzzz[k] = -g_0_z_x_xxxzzz[k] - g_x_z_x_xxxzzz[k] * ab_x + g_x_z_x_xxxxzzz[k];

                g_x_z_xx_xxyyyy[k] = -g_0_z_x_xxyyyy[k] - g_x_z_x_xxyyyy[k] * ab_x + g_x_z_x_xxxyyyy[k];

                g_x_z_xx_xxyyyz[k] = -g_0_z_x_xxyyyz[k] - g_x_z_x_xxyyyz[k] * ab_x + g_x_z_x_xxxyyyz[k];

                g_x_z_xx_xxyyzz[k] = -g_0_z_x_xxyyzz[k] - g_x_z_x_xxyyzz[k] * ab_x + g_x_z_x_xxxyyzz[k];

                g_x_z_xx_xxyzzz[k] = -g_0_z_x_xxyzzz[k] - g_x_z_x_xxyzzz[k] * ab_x + g_x_z_x_xxxyzzz[k];

                g_x_z_xx_xxzzzz[k] = -g_0_z_x_xxzzzz[k] - g_x_z_x_xxzzzz[k] * ab_x + g_x_z_x_xxxzzzz[k];

                g_x_z_xx_xyyyyy[k] = -g_0_z_x_xyyyyy[k] - g_x_z_x_xyyyyy[k] * ab_x + g_x_z_x_xxyyyyy[k];

                g_x_z_xx_xyyyyz[k] = -g_0_z_x_xyyyyz[k] - g_x_z_x_xyyyyz[k] * ab_x + g_x_z_x_xxyyyyz[k];

                g_x_z_xx_xyyyzz[k] = -g_0_z_x_xyyyzz[k] - g_x_z_x_xyyyzz[k] * ab_x + g_x_z_x_xxyyyzz[k];

                g_x_z_xx_xyyzzz[k] = -g_0_z_x_xyyzzz[k] - g_x_z_x_xyyzzz[k] * ab_x + g_x_z_x_xxyyzzz[k];

                g_x_z_xx_xyzzzz[k] = -g_0_z_x_xyzzzz[k] - g_x_z_x_xyzzzz[k] * ab_x + g_x_z_x_xxyzzzz[k];

                g_x_z_xx_xzzzzz[k] = -g_0_z_x_xzzzzz[k] - g_x_z_x_xzzzzz[k] * ab_x + g_x_z_x_xxzzzzz[k];

                g_x_z_xx_yyyyyy[k] = -g_0_z_x_yyyyyy[k] - g_x_z_x_yyyyyy[k] * ab_x + g_x_z_x_xyyyyyy[k];

                g_x_z_xx_yyyyyz[k] = -g_0_z_x_yyyyyz[k] - g_x_z_x_yyyyyz[k] * ab_x + g_x_z_x_xyyyyyz[k];

                g_x_z_xx_yyyyzz[k] = -g_0_z_x_yyyyzz[k] - g_x_z_x_yyyyzz[k] * ab_x + g_x_z_x_xyyyyzz[k];

                g_x_z_xx_yyyzzz[k] = -g_0_z_x_yyyzzz[k] - g_x_z_x_yyyzzz[k] * ab_x + g_x_z_x_xyyyzzz[k];

                g_x_z_xx_yyzzzz[k] = -g_0_z_x_yyzzzz[k] - g_x_z_x_yyzzzz[k] * ab_x + g_x_z_x_xyyzzzz[k];

                g_x_z_xx_yzzzzz[k] = -g_0_z_x_yzzzzz[k] - g_x_z_x_yzzzzz[k] * ab_x + g_x_z_x_xyzzzzz[k];

                g_x_z_xx_zzzzzz[k] = -g_0_z_x_zzzzzz[k] - g_x_z_x_zzzzzz[k] * ab_x + g_x_z_x_xzzzzzz[k];
            }

            /// Set up 364-392 components of targeted buffer : cbuffer.data(

            auto g_x_z_xy_xxxxxx = cbuffer.data(di_geom_11_off + 364 * ccomps * dcomps);

            auto g_x_z_xy_xxxxxy = cbuffer.data(di_geom_11_off + 365 * ccomps * dcomps);

            auto g_x_z_xy_xxxxxz = cbuffer.data(di_geom_11_off + 366 * ccomps * dcomps);

            auto g_x_z_xy_xxxxyy = cbuffer.data(di_geom_11_off + 367 * ccomps * dcomps);

            auto g_x_z_xy_xxxxyz = cbuffer.data(di_geom_11_off + 368 * ccomps * dcomps);

            auto g_x_z_xy_xxxxzz = cbuffer.data(di_geom_11_off + 369 * ccomps * dcomps);

            auto g_x_z_xy_xxxyyy = cbuffer.data(di_geom_11_off + 370 * ccomps * dcomps);

            auto g_x_z_xy_xxxyyz = cbuffer.data(di_geom_11_off + 371 * ccomps * dcomps);

            auto g_x_z_xy_xxxyzz = cbuffer.data(di_geom_11_off + 372 * ccomps * dcomps);

            auto g_x_z_xy_xxxzzz = cbuffer.data(di_geom_11_off + 373 * ccomps * dcomps);

            auto g_x_z_xy_xxyyyy = cbuffer.data(di_geom_11_off + 374 * ccomps * dcomps);

            auto g_x_z_xy_xxyyyz = cbuffer.data(di_geom_11_off + 375 * ccomps * dcomps);

            auto g_x_z_xy_xxyyzz = cbuffer.data(di_geom_11_off + 376 * ccomps * dcomps);

            auto g_x_z_xy_xxyzzz = cbuffer.data(di_geom_11_off + 377 * ccomps * dcomps);

            auto g_x_z_xy_xxzzzz = cbuffer.data(di_geom_11_off + 378 * ccomps * dcomps);

            auto g_x_z_xy_xyyyyy = cbuffer.data(di_geom_11_off + 379 * ccomps * dcomps);

            auto g_x_z_xy_xyyyyz = cbuffer.data(di_geom_11_off + 380 * ccomps * dcomps);

            auto g_x_z_xy_xyyyzz = cbuffer.data(di_geom_11_off + 381 * ccomps * dcomps);

            auto g_x_z_xy_xyyzzz = cbuffer.data(di_geom_11_off + 382 * ccomps * dcomps);

            auto g_x_z_xy_xyzzzz = cbuffer.data(di_geom_11_off + 383 * ccomps * dcomps);

            auto g_x_z_xy_xzzzzz = cbuffer.data(di_geom_11_off + 384 * ccomps * dcomps);

            auto g_x_z_xy_yyyyyy = cbuffer.data(di_geom_11_off + 385 * ccomps * dcomps);

            auto g_x_z_xy_yyyyyz = cbuffer.data(di_geom_11_off + 386 * ccomps * dcomps);

            auto g_x_z_xy_yyyyzz = cbuffer.data(di_geom_11_off + 387 * ccomps * dcomps);

            auto g_x_z_xy_yyyzzz = cbuffer.data(di_geom_11_off + 388 * ccomps * dcomps);

            auto g_x_z_xy_yyzzzz = cbuffer.data(di_geom_11_off + 389 * ccomps * dcomps);

            auto g_x_z_xy_yzzzzz = cbuffer.data(di_geom_11_off + 390 * ccomps * dcomps);

            auto g_x_z_xy_zzzzzz = cbuffer.data(di_geom_11_off + 391 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_x_xxxxxx, g_x_z_x_xxxxxxy, g_x_z_x_xxxxxy, g_x_z_x_xxxxxyy, g_x_z_x_xxxxxyz, g_x_z_x_xxxxxz, g_x_z_x_xxxxyy, g_x_z_x_xxxxyyy, g_x_z_x_xxxxyyz, g_x_z_x_xxxxyz, g_x_z_x_xxxxyzz, g_x_z_x_xxxxzz, g_x_z_x_xxxyyy, g_x_z_x_xxxyyyy, g_x_z_x_xxxyyyz, g_x_z_x_xxxyyz, g_x_z_x_xxxyyzz, g_x_z_x_xxxyzz, g_x_z_x_xxxyzzz, g_x_z_x_xxxzzz, g_x_z_x_xxyyyy, g_x_z_x_xxyyyyy, g_x_z_x_xxyyyyz, g_x_z_x_xxyyyz, g_x_z_x_xxyyyzz, g_x_z_x_xxyyzz, g_x_z_x_xxyyzzz, g_x_z_x_xxyzzz, g_x_z_x_xxyzzzz, g_x_z_x_xxzzzz, g_x_z_x_xyyyyy, g_x_z_x_xyyyyyy, g_x_z_x_xyyyyyz, g_x_z_x_xyyyyz, g_x_z_x_xyyyyzz, g_x_z_x_xyyyzz, g_x_z_x_xyyyzzz, g_x_z_x_xyyzzz, g_x_z_x_xyyzzzz, g_x_z_x_xyzzzz, g_x_z_x_xyzzzzz, g_x_z_x_xzzzzz, g_x_z_x_yyyyyy, g_x_z_x_yyyyyyy, g_x_z_x_yyyyyyz, g_x_z_x_yyyyyz, g_x_z_x_yyyyyzz, g_x_z_x_yyyyzz, g_x_z_x_yyyyzzz, g_x_z_x_yyyzzz, g_x_z_x_yyyzzzz, g_x_z_x_yyzzzz, g_x_z_x_yyzzzzz, g_x_z_x_yzzzzz, g_x_z_x_yzzzzzz, g_x_z_x_zzzzzz, g_x_z_xy_xxxxxx, g_x_z_xy_xxxxxy, g_x_z_xy_xxxxxz, g_x_z_xy_xxxxyy, g_x_z_xy_xxxxyz, g_x_z_xy_xxxxzz, g_x_z_xy_xxxyyy, g_x_z_xy_xxxyyz, g_x_z_xy_xxxyzz, g_x_z_xy_xxxzzz, g_x_z_xy_xxyyyy, g_x_z_xy_xxyyyz, g_x_z_xy_xxyyzz, g_x_z_xy_xxyzzz, g_x_z_xy_xxzzzz, g_x_z_xy_xyyyyy, g_x_z_xy_xyyyyz, g_x_z_xy_xyyyzz, g_x_z_xy_xyyzzz, g_x_z_xy_xyzzzz, g_x_z_xy_xzzzzz, g_x_z_xy_yyyyyy, g_x_z_xy_yyyyyz, g_x_z_xy_yyyyzz, g_x_z_xy_yyyzzz, g_x_z_xy_yyzzzz, g_x_z_xy_yzzzzz, g_x_z_xy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xy_xxxxxx[k] = -g_x_z_x_xxxxxx[k] * ab_y + g_x_z_x_xxxxxxy[k];

                g_x_z_xy_xxxxxy[k] = -g_x_z_x_xxxxxy[k] * ab_y + g_x_z_x_xxxxxyy[k];

                g_x_z_xy_xxxxxz[k] = -g_x_z_x_xxxxxz[k] * ab_y + g_x_z_x_xxxxxyz[k];

                g_x_z_xy_xxxxyy[k] = -g_x_z_x_xxxxyy[k] * ab_y + g_x_z_x_xxxxyyy[k];

                g_x_z_xy_xxxxyz[k] = -g_x_z_x_xxxxyz[k] * ab_y + g_x_z_x_xxxxyyz[k];

                g_x_z_xy_xxxxzz[k] = -g_x_z_x_xxxxzz[k] * ab_y + g_x_z_x_xxxxyzz[k];

                g_x_z_xy_xxxyyy[k] = -g_x_z_x_xxxyyy[k] * ab_y + g_x_z_x_xxxyyyy[k];

                g_x_z_xy_xxxyyz[k] = -g_x_z_x_xxxyyz[k] * ab_y + g_x_z_x_xxxyyyz[k];

                g_x_z_xy_xxxyzz[k] = -g_x_z_x_xxxyzz[k] * ab_y + g_x_z_x_xxxyyzz[k];

                g_x_z_xy_xxxzzz[k] = -g_x_z_x_xxxzzz[k] * ab_y + g_x_z_x_xxxyzzz[k];

                g_x_z_xy_xxyyyy[k] = -g_x_z_x_xxyyyy[k] * ab_y + g_x_z_x_xxyyyyy[k];

                g_x_z_xy_xxyyyz[k] = -g_x_z_x_xxyyyz[k] * ab_y + g_x_z_x_xxyyyyz[k];

                g_x_z_xy_xxyyzz[k] = -g_x_z_x_xxyyzz[k] * ab_y + g_x_z_x_xxyyyzz[k];

                g_x_z_xy_xxyzzz[k] = -g_x_z_x_xxyzzz[k] * ab_y + g_x_z_x_xxyyzzz[k];

                g_x_z_xy_xxzzzz[k] = -g_x_z_x_xxzzzz[k] * ab_y + g_x_z_x_xxyzzzz[k];

                g_x_z_xy_xyyyyy[k] = -g_x_z_x_xyyyyy[k] * ab_y + g_x_z_x_xyyyyyy[k];

                g_x_z_xy_xyyyyz[k] = -g_x_z_x_xyyyyz[k] * ab_y + g_x_z_x_xyyyyyz[k];

                g_x_z_xy_xyyyzz[k] = -g_x_z_x_xyyyzz[k] * ab_y + g_x_z_x_xyyyyzz[k];

                g_x_z_xy_xyyzzz[k] = -g_x_z_x_xyyzzz[k] * ab_y + g_x_z_x_xyyyzzz[k];

                g_x_z_xy_xyzzzz[k] = -g_x_z_x_xyzzzz[k] * ab_y + g_x_z_x_xyyzzzz[k];

                g_x_z_xy_xzzzzz[k] = -g_x_z_x_xzzzzz[k] * ab_y + g_x_z_x_xyzzzzz[k];

                g_x_z_xy_yyyyyy[k] = -g_x_z_x_yyyyyy[k] * ab_y + g_x_z_x_yyyyyyy[k];

                g_x_z_xy_yyyyyz[k] = -g_x_z_x_yyyyyz[k] * ab_y + g_x_z_x_yyyyyyz[k];

                g_x_z_xy_yyyyzz[k] = -g_x_z_x_yyyyzz[k] * ab_y + g_x_z_x_yyyyyzz[k];

                g_x_z_xy_yyyzzz[k] = -g_x_z_x_yyyzzz[k] * ab_y + g_x_z_x_yyyyzzz[k];

                g_x_z_xy_yyzzzz[k] = -g_x_z_x_yyzzzz[k] * ab_y + g_x_z_x_yyyzzzz[k];

                g_x_z_xy_yzzzzz[k] = -g_x_z_x_yzzzzz[k] * ab_y + g_x_z_x_yyzzzzz[k];

                g_x_z_xy_zzzzzz[k] = -g_x_z_x_zzzzzz[k] * ab_y + g_x_z_x_yzzzzzz[k];
            }

            /// Set up 392-420 components of targeted buffer : cbuffer.data(

            auto g_x_z_xz_xxxxxx = cbuffer.data(di_geom_11_off + 392 * ccomps * dcomps);

            auto g_x_z_xz_xxxxxy = cbuffer.data(di_geom_11_off + 393 * ccomps * dcomps);

            auto g_x_z_xz_xxxxxz = cbuffer.data(di_geom_11_off + 394 * ccomps * dcomps);

            auto g_x_z_xz_xxxxyy = cbuffer.data(di_geom_11_off + 395 * ccomps * dcomps);

            auto g_x_z_xz_xxxxyz = cbuffer.data(di_geom_11_off + 396 * ccomps * dcomps);

            auto g_x_z_xz_xxxxzz = cbuffer.data(di_geom_11_off + 397 * ccomps * dcomps);

            auto g_x_z_xz_xxxyyy = cbuffer.data(di_geom_11_off + 398 * ccomps * dcomps);

            auto g_x_z_xz_xxxyyz = cbuffer.data(di_geom_11_off + 399 * ccomps * dcomps);

            auto g_x_z_xz_xxxyzz = cbuffer.data(di_geom_11_off + 400 * ccomps * dcomps);

            auto g_x_z_xz_xxxzzz = cbuffer.data(di_geom_11_off + 401 * ccomps * dcomps);

            auto g_x_z_xz_xxyyyy = cbuffer.data(di_geom_11_off + 402 * ccomps * dcomps);

            auto g_x_z_xz_xxyyyz = cbuffer.data(di_geom_11_off + 403 * ccomps * dcomps);

            auto g_x_z_xz_xxyyzz = cbuffer.data(di_geom_11_off + 404 * ccomps * dcomps);

            auto g_x_z_xz_xxyzzz = cbuffer.data(di_geom_11_off + 405 * ccomps * dcomps);

            auto g_x_z_xz_xxzzzz = cbuffer.data(di_geom_11_off + 406 * ccomps * dcomps);

            auto g_x_z_xz_xyyyyy = cbuffer.data(di_geom_11_off + 407 * ccomps * dcomps);

            auto g_x_z_xz_xyyyyz = cbuffer.data(di_geom_11_off + 408 * ccomps * dcomps);

            auto g_x_z_xz_xyyyzz = cbuffer.data(di_geom_11_off + 409 * ccomps * dcomps);

            auto g_x_z_xz_xyyzzz = cbuffer.data(di_geom_11_off + 410 * ccomps * dcomps);

            auto g_x_z_xz_xyzzzz = cbuffer.data(di_geom_11_off + 411 * ccomps * dcomps);

            auto g_x_z_xz_xzzzzz = cbuffer.data(di_geom_11_off + 412 * ccomps * dcomps);

            auto g_x_z_xz_yyyyyy = cbuffer.data(di_geom_11_off + 413 * ccomps * dcomps);

            auto g_x_z_xz_yyyyyz = cbuffer.data(di_geom_11_off + 414 * ccomps * dcomps);

            auto g_x_z_xz_yyyyzz = cbuffer.data(di_geom_11_off + 415 * ccomps * dcomps);

            auto g_x_z_xz_yyyzzz = cbuffer.data(di_geom_11_off + 416 * ccomps * dcomps);

            auto g_x_z_xz_yyzzzz = cbuffer.data(di_geom_11_off + 417 * ccomps * dcomps);

            auto g_x_z_xz_yzzzzz = cbuffer.data(di_geom_11_off + 418 * ccomps * dcomps);

            auto g_x_z_xz_zzzzzz = cbuffer.data(di_geom_11_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxxxx, g_0_z_z_xxxxxy, g_0_z_z_xxxxxz, g_0_z_z_xxxxyy, g_0_z_z_xxxxyz, g_0_z_z_xxxxzz, g_0_z_z_xxxyyy, g_0_z_z_xxxyyz, g_0_z_z_xxxyzz, g_0_z_z_xxxzzz, g_0_z_z_xxyyyy, g_0_z_z_xxyyyz, g_0_z_z_xxyyzz, g_0_z_z_xxyzzz, g_0_z_z_xxzzzz, g_0_z_z_xyyyyy, g_0_z_z_xyyyyz, g_0_z_z_xyyyzz, g_0_z_z_xyyzzz, g_0_z_z_xyzzzz, g_0_z_z_xzzzzz, g_0_z_z_yyyyyy, g_0_z_z_yyyyyz, g_0_z_z_yyyyzz, g_0_z_z_yyyzzz, g_0_z_z_yyzzzz, g_0_z_z_yzzzzz, g_0_z_z_zzzzzz, g_x_z_xz_xxxxxx, g_x_z_xz_xxxxxy, g_x_z_xz_xxxxxz, g_x_z_xz_xxxxyy, g_x_z_xz_xxxxyz, g_x_z_xz_xxxxzz, g_x_z_xz_xxxyyy, g_x_z_xz_xxxyyz, g_x_z_xz_xxxyzz, g_x_z_xz_xxxzzz, g_x_z_xz_xxyyyy, g_x_z_xz_xxyyyz, g_x_z_xz_xxyyzz, g_x_z_xz_xxyzzz, g_x_z_xz_xxzzzz, g_x_z_xz_xyyyyy, g_x_z_xz_xyyyyz, g_x_z_xz_xyyyzz, g_x_z_xz_xyyzzz, g_x_z_xz_xyzzzz, g_x_z_xz_xzzzzz, g_x_z_xz_yyyyyy, g_x_z_xz_yyyyyz, g_x_z_xz_yyyyzz, g_x_z_xz_yyyzzz, g_x_z_xz_yyzzzz, g_x_z_xz_yzzzzz, g_x_z_xz_zzzzzz, g_x_z_z_xxxxxx, g_x_z_z_xxxxxxx, g_x_z_z_xxxxxxy, g_x_z_z_xxxxxxz, g_x_z_z_xxxxxy, g_x_z_z_xxxxxyy, g_x_z_z_xxxxxyz, g_x_z_z_xxxxxz, g_x_z_z_xxxxxzz, g_x_z_z_xxxxyy, g_x_z_z_xxxxyyy, g_x_z_z_xxxxyyz, g_x_z_z_xxxxyz, g_x_z_z_xxxxyzz, g_x_z_z_xxxxzz, g_x_z_z_xxxxzzz, g_x_z_z_xxxyyy, g_x_z_z_xxxyyyy, g_x_z_z_xxxyyyz, g_x_z_z_xxxyyz, g_x_z_z_xxxyyzz, g_x_z_z_xxxyzz, g_x_z_z_xxxyzzz, g_x_z_z_xxxzzz, g_x_z_z_xxxzzzz, g_x_z_z_xxyyyy, g_x_z_z_xxyyyyy, g_x_z_z_xxyyyyz, g_x_z_z_xxyyyz, g_x_z_z_xxyyyzz, g_x_z_z_xxyyzz, g_x_z_z_xxyyzzz, g_x_z_z_xxyzzz, g_x_z_z_xxyzzzz, g_x_z_z_xxzzzz, g_x_z_z_xxzzzzz, g_x_z_z_xyyyyy, g_x_z_z_xyyyyyy, g_x_z_z_xyyyyyz, g_x_z_z_xyyyyz, g_x_z_z_xyyyyzz, g_x_z_z_xyyyzz, g_x_z_z_xyyyzzz, g_x_z_z_xyyzzz, g_x_z_z_xyyzzzz, g_x_z_z_xyzzzz, g_x_z_z_xyzzzzz, g_x_z_z_xzzzzz, g_x_z_z_xzzzzzz, g_x_z_z_yyyyyy, g_x_z_z_yyyyyz, g_x_z_z_yyyyzz, g_x_z_z_yyyzzz, g_x_z_z_yyzzzz, g_x_z_z_yzzzzz, g_x_z_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xz_xxxxxx[k] = -g_0_z_z_xxxxxx[k] - g_x_z_z_xxxxxx[k] * ab_x + g_x_z_z_xxxxxxx[k];

                g_x_z_xz_xxxxxy[k] = -g_0_z_z_xxxxxy[k] - g_x_z_z_xxxxxy[k] * ab_x + g_x_z_z_xxxxxxy[k];

                g_x_z_xz_xxxxxz[k] = -g_0_z_z_xxxxxz[k] - g_x_z_z_xxxxxz[k] * ab_x + g_x_z_z_xxxxxxz[k];

                g_x_z_xz_xxxxyy[k] = -g_0_z_z_xxxxyy[k] - g_x_z_z_xxxxyy[k] * ab_x + g_x_z_z_xxxxxyy[k];

                g_x_z_xz_xxxxyz[k] = -g_0_z_z_xxxxyz[k] - g_x_z_z_xxxxyz[k] * ab_x + g_x_z_z_xxxxxyz[k];

                g_x_z_xz_xxxxzz[k] = -g_0_z_z_xxxxzz[k] - g_x_z_z_xxxxzz[k] * ab_x + g_x_z_z_xxxxxzz[k];

                g_x_z_xz_xxxyyy[k] = -g_0_z_z_xxxyyy[k] - g_x_z_z_xxxyyy[k] * ab_x + g_x_z_z_xxxxyyy[k];

                g_x_z_xz_xxxyyz[k] = -g_0_z_z_xxxyyz[k] - g_x_z_z_xxxyyz[k] * ab_x + g_x_z_z_xxxxyyz[k];

                g_x_z_xz_xxxyzz[k] = -g_0_z_z_xxxyzz[k] - g_x_z_z_xxxyzz[k] * ab_x + g_x_z_z_xxxxyzz[k];

                g_x_z_xz_xxxzzz[k] = -g_0_z_z_xxxzzz[k] - g_x_z_z_xxxzzz[k] * ab_x + g_x_z_z_xxxxzzz[k];

                g_x_z_xz_xxyyyy[k] = -g_0_z_z_xxyyyy[k] - g_x_z_z_xxyyyy[k] * ab_x + g_x_z_z_xxxyyyy[k];

                g_x_z_xz_xxyyyz[k] = -g_0_z_z_xxyyyz[k] - g_x_z_z_xxyyyz[k] * ab_x + g_x_z_z_xxxyyyz[k];

                g_x_z_xz_xxyyzz[k] = -g_0_z_z_xxyyzz[k] - g_x_z_z_xxyyzz[k] * ab_x + g_x_z_z_xxxyyzz[k];

                g_x_z_xz_xxyzzz[k] = -g_0_z_z_xxyzzz[k] - g_x_z_z_xxyzzz[k] * ab_x + g_x_z_z_xxxyzzz[k];

                g_x_z_xz_xxzzzz[k] = -g_0_z_z_xxzzzz[k] - g_x_z_z_xxzzzz[k] * ab_x + g_x_z_z_xxxzzzz[k];

                g_x_z_xz_xyyyyy[k] = -g_0_z_z_xyyyyy[k] - g_x_z_z_xyyyyy[k] * ab_x + g_x_z_z_xxyyyyy[k];

                g_x_z_xz_xyyyyz[k] = -g_0_z_z_xyyyyz[k] - g_x_z_z_xyyyyz[k] * ab_x + g_x_z_z_xxyyyyz[k];

                g_x_z_xz_xyyyzz[k] = -g_0_z_z_xyyyzz[k] - g_x_z_z_xyyyzz[k] * ab_x + g_x_z_z_xxyyyzz[k];

                g_x_z_xz_xyyzzz[k] = -g_0_z_z_xyyzzz[k] - g_x_z_z_xyyzzz[k] * ab_x + g_x_z_z_xxyyzzz[k];

                g_x_z_xz_xyzzzz[k] = -g_0_z_z_xyzzzz[k] - g_x_z_z_xyzzzz[k] * ab_x + g_x_z_z_xxyzzzz[k];

                g_x_z_xz_xzzzzz[k] = -g_0_z_z_xzzzzz[k] - g_x_z_z_xzzzzz[k] * ab_x + g_x_z_z_xxzzzzz[k];

                g_x_z_xz_yyyyyy[k] = -g_0_z_z_yyyyyy[k] - g_x_z_z_yyyyyy[k] * ab_x + g_x_z_z_xyyyyyy[k];

                g_x_z_xz_yyyyyz[k] = -g_0_z_z_yyyyyz[k] - g_x_z_z_yyyyyz[k] * ab_x + g_x_z_z_xyyyyyz[k];

                g_x_z_xz_yyyyzz[k] = -g_0_z_z_yyyyzz[k] - g_x_z_z_yyyyzz[k] * ab_x + g_x_z_z_xyyyyzz[k];

                g_x_z_xz_yyyzzz[k] = -g_0_z_z_yyyzzz[k] - g_x_z_z_yyyzzz[k] * ab_x + g_x_z_z_xyyyzzz[k];

                g_x_z_xz_yyzzzz[k] = -g_0_z_z_yyzzzz[k] - g_x_z_z_yyzzzz[k] * ab_x + g_x_z_z_xyyzzzz[k];

                g_x_z_xz_yzzzzz[k] = -g_0_z_z_yzzzzz[k] - g_x_z_z_yzzzzz[k] * ab_x + g_x_z_z_xyzzzzz[k];

                g_x_z_xz_zzzzzz[k] = -g_0_z_z_zzzzzz[k] - g_x_z_z_zzzzzz[k] * ab_x + g_x_z_z_xzzzzzz[k];
            }

            /// Set up 420-448 components of targeted buffer : cbuffer.data(

            auto g_x_z_yy_xxxxxx = cbuffer.data(di_geom_11_off + 420 * ccomps * dcomps);

            auto g_x_z_yy_xxxxxy = cbuffer.data(di_geom_11_off + 421 * ccomps * dcomps);

            auto g_x_z_yy_xxxxxz = cbuffer.data(di_geom_11_off + 422 * ccomps * dcomps);

            auto g_x_z_yy_xxxxyy = cbuffer.data(di_geom_11_off + 423 * ccomps * dcomps);

            auto g_x_z_yy_xxxxyz = cbuffer.data(di_geom_11_off + 424 * ccomps * dcomps);

            auto g_x_z_yy_xxxxzz = cbuffer.data(di_geom_11_off + 425 * ccomps * dcomps);

            auto g_x_z_yy_xxxyyy = cbuffer.data(di_geom_11_off + 426 * ccomps * dcomps);

            auto g_x_z_yy_xxxyyz = cbuffer.data(di_geom_11_off + 427 * ccomps * dcomps);

            auto g_x_z_yy_xxxyzz = cbuffer.data(di_geom_11_off + 428 * ccomps * dcomps);

            auto g_x_z_yy_xxxzzz = cbuffer.data(di_geom_11_off + 429 * ccomps * dcomps);

            auto g_x_z_yy_xxyyyy = cbuffer.data(di_geom_11_off + 430 * ccomps * dcomps);

            auto g_x_z_yy_xxyyyz = cbuffer.data(di_geom_11_off + 431 * ccomps * dcomps);

            auto g_x_z_yy_xxyyzz = cbuffer.data(di_geom_11_off + 432 * ccomps * dcomps);

            auto g_x_z_yy_xxyzzz = cbuffer.data(di_geom_11_off + 433 * ccomps * dcomps);

            auto g_x_z_yy_xxzzzz = cbuffer.data(di_geom_11_off + 434 * ccomps * dcomps);

            auto g_x_z_yy_xyyyyy = cbuffer.data(di_geom_11_off + 435 * ccomps * dcomps);

            auto g_x_z_yy_xyyyyz = cbuffer.data(di_geom_11_off + 436 * ccomps * dcomps);

            auto g_x_z_yy_xyyyzz = cbuffer.data(di_geom_11_off + 437 * ccomps * dcomps);

            auto g_x_z_yy_xyyzzz = cbuffer.data(di_geom_11_off + 438 * ccomps * dcomps);

            auto g_x_z_yy_xyzzzz = cbuffer.data(di_geom_11_off + 439 * ccomps * dcomps);

            auto g_x_z_yy_xzzzzz = cbuffer.data(di_geom_11_off + 440 * ccomps * dcomps);

            auto g_x_z_yy_yyyyyy = cbuffer.data(di_geom_11_off + 441 * ccomps * dcomps);

            auto g_x_z_yy_yyyyyz = cbuffer.data(di_geom_11_off + 442 * ccomps * dcomps);

            auto g_x_z_yy_yyyyzz = cbuffer.data(di_geom_11_off + 443 * ccomps * dcomps);

            auto g_x_z_yy_yyyzzz = cbuffer.data(di_geom_11_off + 444 * ccomps * dcomps);

            auto g_x_z_yy_yyzzzz = cbuffer.data(di_geom_11_off + 445 * ccomps * dcomps);

            auto g_x_z_yy_yzzzzz = cbuffer.data(di_geom_11_off + 446 * ccomps * dcomps);

            auto g_x_z_yy_zzzzzz = cbuffer.data(di_geom_11_off + 447 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_y_xxxxxx, g_x_z_y_xxxxxxy, g_x_z_y_xxxxxy, g_x_z_y_xxxxxyy, g_x_z_y_xxxxxyz, g_x_z_y_xxxxxz, g_x_z_y_xxxxyy, g_x_z_y_xxxxyyy, g_x_z_y_xxxxyyz, g_x_z_y_xxxxyz, g_x_z_y_xxxxyzz, g_x_z_y_xxxxzz, g_x_z_y_xxxyyy, g_x_z_y_xxxyyyy, g_x_z_y_xxxyyyz, g_x_z_y_xxxyyz, g_x_z_y_xxxyyzz, g_x_z_y_xxxyzz, g_x_z_y_xxxyzzz, g_x_z_y_xxxzzz, g_x_z_y_xxyyyy, g_x_z_y_xxyyyyy, g_x_z_y_xxyyyyz, g_x_z_y_xxyyyz, g_x_z_y_xxyyyzz, g_x_z_y_xxyyzz, g_x_z_y_xxyyzzz, g_x_z_y_xxyzzz, g_x_z_y_xxyzzzz, g_x_z_y_xxzzzz, g_x_z_y_xyyyyy, g_x_z_y_xyyyyyy, g_x_z_y_xyyyyyz, g_x_z_y_xyyyyz, g_x_z_y_xyyyyzz, g_x_z_y_xyyyzz, g_x_z_y_xyyyzzz, g_x_z_y_xyyzzz, g_x_z_y_xyyzzzz, g_x_z_y_xyzzzz, g_x_z_y_xyzzzzz, g_x_z_y_xzzzzz, g_x_z_y_yyyyyy, g_x_z_y_yyyyyyy, g_x_z_y_yyyyyyz, g_x_z_y_yyyyyz, g_x_z_y_yyyyyzz, g_x_z_y_yyyyzz, g_x_z_y_yyyyzzz, g_x_z_y_yyyzzz, g_x_z_y_yyyzzzz, g_x_z_y_yyzzzz, g_x_z_y_yyzzzzz, g_x_z_y_yzzzzz, g_x_z_y_yzzzzzz, g_x_z_y_zzzzzz, g_x_z_yy_xxxxxx, g_x_z_yy_xxxxxy, g_x_z_yy_xxxxxz, g_x_z_yy_xxxxyy, g_x_z_yy_xxxxyz, g_x_z_yy_xxxxzz, g_x_z_yy_xxxyyy, g_x_z_yy_xxxyyz, g_x_z_yy_xxxyzz, g_x_z_yy_xxxzzz, g_x_z_yy_xxyyyy, g_x_z_yy_xxyyyz, g_x_z_yy_xxyyzz, g_x_z_yy_xxyzzz, g_x_z_yy_xxzzzz, g_x_z_yy_xyyyyy, g_x_z_yy_xyyyyz, g_x_z_yy_xyyyzz, g_x_z_yy_xyyzzz, g_x_z_yy_xyzzzz, g_x_z_yy_xzzzzz, g_x_z_yy_yyyyyy, g_x_z_yy_yyyyyz, g_x_z_yy_yyyyzz, g_x_z_yy_yyyzzz, g_x_z_yy_yyzzzz, g_x_z_yy_yzzzzz, g_x_z_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yy_xxxxxx[k] = -g_x_z_y_xxxxxx[k] * ab_y + g_x_z_y_xxxxxxy[k];

                g_x_z_yy_xxxxxy[k] = -g_x_z_y_xxxxxy[k] * ab_y + g_x_z_y_xxxxxyy[k];

                g_x_z_yy_xxxxxz[k] = -g_x_z_y_xxxxxz[k] * ab_y + g_x_z_y_xxxxxyz[k];

                g_x_z_yy_xxxxyy[k] = -g_x_z_y_xxxxyy[k] * ab_y + g_x_z_y_xxxxyyy[k];

                g_x_z_yy_xxxxyz[k] = -g_x_z_y_xxxxyz[k] * ab_y + g_x_z_y_xxxxyyz[k];

                g_x_z_yy_xxxxzz[k] = -g_x_z_y_xxxxzz[k] * ab_y + g_x_z_y_xxxxyzz[k];

                g_x_z_yy_xxxyyy[k] = -g_x_z_y_xxxyyy[k] * ab_y + g_x_z_y_xxxyyyy[k];

                g_x_z_yy_xxxyyz[k] = -g_x_z_y_xxxyyz[k] * ab_y + g_x_z_y_xxxyyyz[k];

                g_x_z_yy_xxxyzz[k] = -g_x_z_y_xxxyzz[k] * ab_y + g_x_z_y_xxxyyzz[k];

                g_x_z_yy_xxxzzz[k] = -g_x_z_y_xxxzzz[k] * ab_y + g_x_z_y_xxxyzzz[k];

                g_x_z_yy_xxyyyy[k] = -g_x_z_y_xxyyyy[k] * ab_y + g_x_z_y_xxyyyyy[k];

                g_x_z_yy_xxyyyz[k] = -g_x_z_y_xxyyyz[k] * ab_y + g_x_z_y_xxyyyyz[k];

                g_x_z_yy_xxyyzz[k] = -g_x_z_y_xxyyzz[k] * ab_y + g_x_z_y_xxyyyzz[k];

                g_x_z_yy_xxyzzz[k] = -g_x_z_y_xxyzzz[k] * ab_y + g_x_z_y_xxyyzzz[k];

                g_x_z_yy_xxzzzz[k] = -g_x_z_y_xxzzzz[k] * ab_y + g_x_z_y_xxyzzzz[k];

                g_x_z_yy_xyyyyy[k] = -g_x_z_y_xyyyyy[k] * ab_y + g_x_z_y_xyyyyyy[k];

                g_x_z_yy_xyyyyz[k] = -g_x_z_y_xyyyyz[k] * ab_y + g_x_z_y_xyyyyyz[k];

                g_x_z_yy_xyyyzz[k] = -g_x_z_y_xyyyzz[k] * ab_y + g_x_z_y_xyyyyzz[k];

                g_x_z_yy_xyyzzz[k] = -g_x_z_y_xyyzzz[k] * ab_y + g_x_z_y_xyyyzzz[k];

                g_x_z_yy_xyzzzz[k] = -g_x_z_y_xyzzzz[k] * ab_y + g_x_z_y_xyyzzzz[k];

                g_x_z_yy_xzzzzz[k] = -g_x_z_y_xzzzzz[k] * ab_y + g_x_z_y_xyzzzzz[k];

                g_x_z_yy_yyyyyy[k] = -g_x_z_y_yyyyyy[k] * ab_y + g_x_z_y_yyyyyyy[k];

                g_x_z_yy_yyyyyz[k] = -g_x_z_y_yyyyyz[k] * ab_y + g_x_z_y_yyyyyyz[k];

                g_x_z_yy_yyyyzz[k] = -g_x_z_y_yyyyzz[k] * ab_y + g_x_z_y_yyyyyzz[k];

                g_x_z_yy_yyyzzz[k] = -g_x_z_y_yyyzzz[k] * ab_y + g_x_z_y_yyyyzzz[k];

                g_x_z_yy_yyzzzz[k] = -g_x_z_y_yyzzzz[k] * ab_y + g_x_z_y_yyyzzzz[k];

                g_x_z_yy_yzzzzz[k] = -g_x_z_y_yzzzzz[k] * ab_y + g_x_z_y_yyzzzzz[k];

                g_x_z_yy_zzzzzz[k] = -g_x_z_y_zzzzzz[k] * ab_y + g_x_z_y_yzzzzzz[k];
            }

            /// Set up 448-476 components of targeted buffer : cbuffer.data(

            auto g_x_z_yz_xxxxxx = cbuffer.data(di_geom_11_off + 448 * ccomps * dcomps);

            auto g_x_z_yz_xxxxxy = cbuffer.data(di_geom_11_off + 449 * ccomps * dcomps);

            auto g_x_z_yz_xxxxxz = cbuffer.data(di_geom_11_off + 450 * ccomps * dcomps);

            auto g_x_z_yz_xxxxyy = cbuffer.data(di_geom_11_off + 451 * ccomps * dcomps);

            auto g_x_z_yz_xxxxyz = cbuffer.data(di_geom_11_off + 452 * ccomps * dcomps);

            auto g_x_z_yz_xxxxzz = cbuffer.data(di_geom_11_off + 453 * ccomps * dcomps);

            auto g_x_z_yz_xxxyyy = cbuffer.data(di_geom_11_off + 454 * ccomps * dcomps);

            auto g_x_z_yz_xxxyyz = cbuffer.data(di_geom_11_off + 455 * ccomps * dcomps);

            auto g_x_z_yz_xxxyzz = cbuffer.data(di_geom_11_off + 456 * ccomps * dcomps);

            auto g_x_z_yz_xxxzzz = cbuffer.data(di_geom_11_off + 457 * ccomps * dcomps);

            auto g_x_z_yz_xxyyyy = cbuffer.data(di_geom_11_off + 458 * ccomps * dcomps);

            auto g_x_z_yz_xxyyyz = cbuffer.data(di_geom_11_off + 459 * ccomps * dcomps);

            auto g_x_z_yz_xxyyzz = cbuffer.data(di_geom_11_off + 460 * ccomps * dcomps);

            auto g_x_z_yz_xxyzzz = cbuffer.data(di_geom_11_off + 461 * ccomps * dcomps);

            auto g_x_z_yz_xxzzzz = cbuffer.data(di_geom_11_off + 462 * ccomps * dcomps);

            auto g_x_z_yz_xyyyyy = cbuffer.data(di_geom_11_off + 463 * ccomps * dcomps);

            auto g_x_z_yz_xyyyyz = cbuffer.data(di_geom_11_off + 464 * ccomps * dcomps);

            auto g_x_z_yz_xyyyzz = cbuffer.data(di_geom_11_off + 465 * ccomps * dcomps);

            auto g_x_z_yz_xyyzzz = cbuffer.data(di_geom_11_off + 466 * ccomps * dcomps);

            auto g_x_z_yz_xyzzzz = cbuffer.data(di_geom_11_off + 467 * ccomps * dcomps);

            auto g_x_z_yz_xzzzzz = cbuffer.data(di_geom_11_off + 468 * ccomps * dcomps);

            auto g_x_z_yz_yyyyyy = cbuffer.data(di_geom_11_off + 469 * ccomps * dcomps);

            auto g_x_z_yz_yyyyyz = cbuffer.data(di_geom_11_off + 470 * ccomps * dcomps);

            auto g_x_z_yz_yyyyzz = cbuffer.data(di_geom_11_off + 471 * ccomps * dcomps);

            auto g_x_z_yz_yyyzzz = cbuffer.data(di_geom_11_off + 472 * ccomps * dcomps);

            auto g_x_z_yz_yyzzzz = cbuffer.data(di_geom_11_off + 473 * ccomps * dcomps);

            auto g_x_z_yz_yzzzzz = cbuffer.data(di_geom_11_off + 474 * ccomps * dcomps);

            auto g_x_z_yz_zzzzzz = cbuffer.data(di_geom_11_off + 475 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yz_xxxxxx, g_x_z_yz_xxxxxy, g_x_z_yz_xxxxxz, g_x_z_yz_xxxxyy, g_x_z_yz_xxxxyz, g_x_z_yz_xxxxzz, g_x_z_yz_xxxyyy, g_x_z_yz_xxxyyz, g_x_z_yz_xxxyzz, g_x_z_yz_xxxzzz, g_x_z_yz_xxyyyy, g_x_z_yz_xxyyyz, g_x_z_yz_xxyyzz, g_x_z_yz_xxyzzz, g_x_z_yz_xxzzzz, g_x_z_yz_xyyyyy, g_x_z_yz_xyyyyz, g_x_z_yz_xyyyzz, g_x_z_yz_xyyzzz, g_x_z_yz_xyzzzz, g_x_z_yz_xzzzzz, g_x_z_yz_yyyyyy, g_x_z_yz_yyyyyz, g_x_z_yz_yyyyzz, g_x_z_yz_yyyzzz, g_x_z_yz_yyzzzz, g_x_z_yz_yzzzzz, g_x_z_yz_zzzzzz, g_x_z_z_xxxxxx, g_x_z_z_xxxxxxy, g_x_z_z_xxxxxy, g_x_z_z_xxxxxyy, g_x_z_z_xxxxxyz, g_x_z_z_xxxxxz, g_x_z_z_xxxxyy, g_x_z_z_xxxxyyy, g_x_z_z_xxxxyyz, g_x_z_z_xxxxyz, g_x_z_z_xxxxyzz, g_x_z_z_xxxxzz, g_x_z_z_xxxyyy, g_x_z_z_xxxyyyy, g_x_z_z_xxxyyyz, g_x_z_z_xxxyyz, g_x_z_z_xxxyyzz, g_x_z_z_xxxyzz, g_x_z_z_xxxyzzz, g_x_z_z_xxxzzz, g_x_z_z_xxyyyy, g_x_z_z_xxyyyyy, g_x_z_z_xxyyyyz, g_x_z_z_xxyyyz, g_x_z_z_xxyyyzz, g_x_z_z_xxyyzz, g_x_z_z_xxyyzzz, g_x_z_z_xxyzzz, g_x_z_z_xxyzzzz, g_x_z_z_xxzzzz, g_x_z_z_xyyyyy, g_x_z_z_xyyyyyy, g_x_z_z_xyyyyyz, g_x_z_z_xyyyyz, g_x_z_z_xyyyyzz, g_x_z_z_xyyyzz, g_x_z_z_xyyyzzz, g_x_z_z_xyyzzz, g_x_z_z_xyyzzzz, g_x_z_z_xyzzzz, g_x_z_z_xyzzzzz, g_x_z_z_xzzzzz, g_x_z_z_yyyyyy, g_x_z_z_yyyyyyy, g_x_z_z_yyyyyyz, g_x_z_z_yyyyyz, g_x_z_z_yyyyyzz, g_x_z_z_yyyyzz, g_x_z_z_yyyyzzz, g_x_z_z_yyyzzz, g_x_z_z_yyyzzzz, g_x_z_z_yyzzzz, g_x_z_z_yyzzzzz, g_x_z_z_yzzzzz, g_x_z_z_yzzzzzz, g_x_z_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yz_xxxxxx[k] = -g_x_z_z_xxxxxx[k] * ab_y + g_x_z_z_xxxxxxy[k];

                g_x_z_yz_xxxxxy[k] = -g_x_z_z_xxxxxy[k] * ab_y + g_x_z_z_xxxxxyy[k];

                g_x_z_yz_xxxxxz[k] = -g_x_z_z_xxxxxz[k] * ab_y + g_x_z_z_xxxxxyz[k];

                g_x_z_yz_xxxxyy[k] = -g_x_z_z_xxxxyy[k] * ab_y + g_x_z_z_xxxxyyy[k];

                g_x_z_yz_xxxxyz[k] = -g_x_z_z_xxxxyz[k] * ab_y + g_x_z_z_xxxxyyz[k];

                g_x_z_yz_xxxxzz[k] = -g_x_z_z_xxxxzz[k] * ab_y + g_x_z_z_xxxxyzz[k];

                g_x_z_yz_xxxyyy[k] = -g_x_z_z_xxxyyy[k] * ab_y + g_x_z_z_xxxyyyy[k];

                g_x_z_yz_xxxyyz[k] = -g_x_z_z_xxxyyz[k] * ab_y + g_x_z_z_xxxyyyz[k];

                g_x_z_yz_xxxyzz[k] = -g_x_z_z_xxxyzz[k] * ab_y + g_x_z_z_xxxyyzz[k];

                g_x_z_yz_xxxzzz[k] = -g_x_z_z_xxxzzz[k] * ab_y + g_x_z_z_xxxyzzz[k];

                g_x_z_yz_xxyyyy[k] = -g_x_z_z_xxyyyy[k] * ab_y + g_x_z_z_xxyyyyy[k];

                g_x_z_yz_xxyyyz[k] = -g_x_z_z_xxyyyz[k] * ab_y + g_x_z_z_xxyyyyz[k];

                g_x_z_yz_xxyyzz[k] = -g_x_z_z_xxyyzz[k] * ab_y + g_x_z_z_xxyyyzz[k];

                g_x_z_yz_xxyzzz[k] = -g_x_z_z_xxyzzz[k] * ab_y + g_x_z_z_xxyyzzz[k];

                g_x_z_yz_xxzzzz[k] = -g_x_z_z_xxzzzz[k] * ab_y + g_x_z_z_xxyzzzz[k];

                g_x_z_yz_xyyyyy[k] = -g_x_z_z_xyyyyy[k] * ab_y + g_x_z_z_xyyyyyy[k];

                g_x_z_yz_xyyyyz[k] = -g_x_z_z_xyyyyz[k] * ab_y + g_x_z_z_xyyyyyz[k];

                g_x_z_yz_xyyyzz[k] = -g_x_z_z_xyyyzz[k] * ab_y + g_x_z_z_xyyyyzz[k];

                g_x_z_yz_xyyzzz[k] = -g_x_z_z_xyyzzz[k] * ab_y + g_x_z_z_xyyyzzz[k];

                g_x_z_yz_xyzzzz[k] = -g_x_z_z_xyzzzz[k] * ab_y + g_x_z_z_xyyzzzz[k];

                g_x_z_yz_xzzzzz[k] = -g_x_z_z_xzzzzz[k] * ab_y + g_x_z_z_xyzzzzz[k];

                g_x_z_yz_yyyyyy[k] = -g_x_z_z_yyyyyy[k] * ab_y + g_x_z_z_yyyyyyy[k];

                g_x_z_yz_yyyyyz[k] = -g_x_z_z_yyyyyz[k] * ab_y + g_x_z_z_yyyyyyz[k];

                g_x_z_yz_yyyyzz[k] = -g_x_z_z_yyyyzz[k] * ab_y + g_x_z_z_yyyyyzz[k];

                g_x_z_yz_yyyzzz[k] = -g_x_z_z_yyyzzz[k] * ab_y + g_x_z_z_yyyyzzz[k];

                g_x_z_yz_yyzzzz[k] = -g_x_z_z_yyzzzz[k] * ab_y + g_x_z_z_yyyzzzz[k];

                g_x_z_yz_yzzzzz[k] = -g_x_z_z_yzzzzz[k] * ab_y + g_x_z_z_yyzzzzz[k];

                g_x_z_yz_zzzzzz[k] = -g_x_z_z_zzzzzz[k] * ab_y + g_x_z_z_yzzzzzz[k];
            }

            /// Set up 476-504 components of targeted buffer : cbuffer.data(

            auto g_x_z_zz_xxxxxx = cbuffer.data(di_geom_11_off + 476 * ccomps * dcomps);

            auto g_x_z_zz_xxxxxy = cbuffer.data(di_geom_11_off + 477 * ccomps * dcomps);

            auto g_x_z_zz_xxxxxz = cbuffer.data(di_geom_11_off + 478 * ccomps * dcomps);

            auto g_x_z_zz_xxxxyy = cbuffer.data(di_geom_11_off + 479 * ccomps * dcomps);

            auto g_x_z_zz_xxxxyz = cbuffer.data(di_geom_11_off + 480 * ccomps * dcomps);

            auto g_x_z_zz_xxxxzz = cbuffer.data(di_geom_11_off + 481 * ccomps * dcomps);

            auto g_x_z_zz_xxxyyy = cbuffer.data(di_geom_11_off + 482 * ccomps * dcomps);

            auto g_x_z_zz_xxxyyz = cbuffer.data(di_geom_11_off + 483 * ccomps * dcomps);

            auto g_x_z_zz_xxxyzz = cbuffer.data(di_geom_11_off + 484 * ccomps * dcomps);

            auto g_x_z_zz_xxxzzz = cbuffer.data(di_geom_11_off + 485 * ccomps * dcomps);

            auto g_x_z_zz_xxyyyy = cbuffer.data(di_geom_11_off + 486 * ccomps * dcomps);

            auto g_x_z_zz_xxyyyz = cbuffer.data(di_geom_11_off + 487 * ccomps * dcomps);

            auto g_x_z_zz_xxyyzz = cbuffer.data(di_geom_11_off + 488 * ccomps * dcomps);

            auto g_x_z_zz_xxyzzz = cbuffer.data(di_geom_11_off + 489 * ccomps * dcomps);

            auto g_x_z_zz_xxzzzz = cbuffer.data(di_geom_11_off + 490 * ccomps * dcomps);

            auto g_x_z_zz_xyyyyy = cbuffer.data(di_geom_11_off + 491 * ccomps * dcomps);

            auto g_x_z_zz_xyyyyz = cbuffer.data(di_geom_11_off + 492 * ccomps * dcomps);

            auto g_x_z_zz_xyyyzz = cbuffer.data(di_geom_11_off + 493 * ccomps * dcomps);

            auto g_x_z_zz_xyyzzz = cbuffer.data(di_geom_11_off + 494 * ccomps * dcomps);

            auto g_x_z_zz_xyzzzz = cbuffer.data(di_geom_11_off + 495 * ccomps * dcomps);

            auto g_x_z_zz_xzzzzz = cbuffer.data(di_geom_11_off + 496 * ccomps * dcomps);

            auto g_x_z_zz_yyyyyy = cbuffer.data(di_geom_11_off + 497 * ccomps * dcomps);

            auto g_x_z_zz_yyyyyz = cbuffer.data(di_geom_11_off + 498 * ccomps * dcomps);

            auto g_x_z_zz_yyyyzz = cbuffer.data(di_geom_11_off + 499 * ccomps * dcomps);

            auto g_x_z_zz_yyyzzz = cbuffer.data(di_geom_11_off + 500 * ccomps * dcomps);

            auto g_x_z_zz_yyzzzz = cbuffer.data(di_geom_11_off + 501 * ccomps * dcomps);

            auto g_x_z_zz_yzzzzz = cbuffer.data(di_geom_11_off + 502 * ccomps * dcomps);

            auto g_x_z_zz_zzzzzz = cbuffer.data(di_geom_11_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxxxxx, g_x_0_z_xxxxxy, g_x_0_z_xxxxxz, g_x_0_z_xxxxyy, g_x_0_z_xxxxyz, g_x_0_z_xxxxzz, g_x_0_z_xxxyyy, g_x_0_z_xxxyyz, g_x_0_z_xxxyzz, g_x_0_z_xxxzzz, g_x_0_z_xxyyyy, g_x_0_z_xxyyyz, g_x_0_z_xxyyzz, g_x_0_z_xxyzzz, g_x_0_z_xxzzzz, g_x_0_z_xyyyyy, g_x_0_z_xyyyyz, g_x_0_z_xyyyzz, g_x_0_z_xyyzzz, g_x_0_z_xyzzzz, g_x_0_z_xzzzzz, g_x_0_z_yyyyyy, g_x_0_z_yyyyyz, g_x_0_z_yyyyzz, g_x_0_z_yyyzzz, g_x_0_z_yyzzzz, g_x_0_z_yzzzzz, g_x_0_z_zzzzzz, g_x_z_z_xxxxxx, g_x_z_z_xxxxxxz, g_x_z_z_xxxxxy, g_x_z_z_xxxxxyz, g_x_z_z_xxxxxz, g_x_z_z_xxxxxzz, g_x_z_z_xxxxyy, g_x_z_z_xxxxyyz, g_x_z_z_xxxxyz, g_x_z_z_xxxxyzz, g_x_z_z_xxxxzz, g_x_z_z_xxxxzzz, g_x_z_z_xxxyyy, g_x_z_z_xxxyyyz, g_x_z_z_xxxyyz, g_x_z_z_xxxyyzz, g_x_z_z_xxxyzz, g_x_z_z_xxxyzzz, g_x_z_z_xxxzzz, g_x_z_z_xxxzzzz, g_x_z_z_xxyyyy, g_x_z_z_xxyyyyz, g_x_z_z_xxyyyz, g_x_z_z_xxyyyzz, g_x_z_z_xxyyzz, g_x_z_z_xxyyzzz, g_x_z_z_xxyzzz, g_x_z_z_xxyzzzz, g_x_z_z_xxzzzz, g_x_z_z_xxzzzzz, g_x_z_z_xyyyyy, g_x_z_z_xyyyyyz, g_x_z_z_xyyyyz, g_x_z_z_xyyyyzz, g_x_z_z_xyyyzz, g_x_z_z_xyyyzzz, g_x_z_z_xyyzzz, g_x_z_z_xyyzzzz, g_x_z_z_xyzzzz, g_x_z_z_xyzzzzz, g_x_z_z_xzzzzz, g_x_z_z_xzzzzzz, g_x_z_z_yyyyyy, g_x_z_z_yyyyyyz, g_x_z_z_yyyyyz, g_x_z_z_yyyyyzz, g_x_z_z_yyyyzz, g_x_z_z_yyyyzzz, g_x_z_z_yyyzzz, g_x_z_z_yyyzzzz, g_x_z_z_yyzzzz, g_x_z_z_yyzzzzz, g_x_z_z_yzzzzz, g_x_z_z_yzzzzzz, g_x_z_z_zzzzzz, g_x_z_z_zzzzzzz, g_x_z_zz_xxxxxx, g_x_z_zz_xxxxxy, g_x_z_zz_xxxxxz, g_x_z_zz_xxxxyy, g_x_z_zz_xxxxyz, g_x_z_zz_xxxxzz, g_x_z_zz_xxxyyy, g_x_z_zz_xxxyyz, g_x_z_zz_xxxyzz, g_x_z_zz_xxxzzz, g_x_z_zz_xxyyyy, g_x_z_zz_xxyyyz, g_x_z_zz_xxyyzz, g_x_z_zz_xxyzzz, g_x_z_zz_xxzzzz, g_x_z_zz_xyyyyy, g_x_z_zz_xyyyyz, g_x_z_zz_xyyyzz, g_x_z_zz_xyyzzz, g_x_z_zz_xyzzzz, g_x_z_zz_xzzzzz, g_x_z_zz_yyyyyy, g_x_z_zz_yyyyyz, g_x_z_zz_yyyyzz, g_x_z_zz_yyyzzz, g_x_z_zz_yyzzzz, g_x_z_zz_yzzzzz, g_x_z_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zz_xxxxxx[k] = g_x_0_z_xxxxxx[k] - g_x_z_z_xxxxxx[k] * ab_z + g_x_z_z_xxxxxxz[k];

                g_x_z_zz_xxxxxy[k] = g_x_0_z_xxxxxy[k] - g_x_z_z_xxxxxy[k] * ab_z + g_x_z_z_xxxxxyz[k];

                g_x_z_zz_xxxxxz[k] = g_x_0_z_xxxxxz[k] - g_x_z_z_xxxxxz[k] * ab_z + g_x_z_z_xxxxxzz[k];

                g_x_z_zz_xxxxyy[k] = g_x_0_z_xxxxyy[k] - g_x_z_z_xxxxyy[k] * ab_z + g_x_z_z_xxxxyyz[k];

                g_x_z_zz_xxxxyz[k] = g_x_0_z_xxxxyz[k] - g_x_z_z_xxxxyz[k] * ab_z + g_x_z_z_xxxxyzz[k];

                g_x_z_zz_xxxxzz[k] = g_x_0_z_xxxxzz[k] - g_x_z_z_xxxxzz[k] * ab_z + g_x_z_z_xxxxzzz[k];

                g_x_z_zz_xxxyyy[k] = g_x_0_z_xxxyyy[k] - g_x_z_z_xxxyyy[k] * ab_z + g_x_z_z_xxxyyyz[k];

                g_x_z_zz_xxxyyz[k] = g_x_0_z_xxxyyz[k] - g_x_z_z_xxxyyz[k] * ab_z + g_x_z_z_xxxyyzz[k];

                g_x_z_zz_xxxyzz[k] = g_x_0_z_xxxyzz[k] - g_x_z_z_xxxyzz[k] * ab_z + g_x_z_z_xxxyzzz[k];

                g_x_z_zz_xxxzzz[k] = g_x_0_z_xxxzzz[k] - g_x_z_z_xxxzzz[k] * ab_z + g_x_z_z_xxxzzzz[k];

                g_x_z_zz_xxyyyy[k] = g_x_0_z_xxyyyy[k] - g_x_z_z_xxyyyy[k] * ab_z + g_x_z_z_xxyyyyz[k];

                g_x_z_zz_xxyyyz[k] = g_x_0_z_xxyyyz[k] - g_x_z_z_xxyyyz[k] * ab_z + g_x_z_z_xxyyyzz[k];

                g_x_z_zz_xxyyzz[k] = g_x_0_z_xxyyzz[k] - g_x_z_z_xxyyzz[k] * ab_z + g_x_z_z_xxyyzzz[k];

                g_x_z_zz_xxyzzz[k] = g_x_0_z_xxyzzz[k] - g_x_z_z_xxyzzz[k] * ab_z + g_x_z_z_xxyzzzz[k];

                g_x_z_zz_xxzzzz[k] = g_x_0_z_xxzzzz[k] - g_x_z_z_xxzzzz[k] * ab_z + g_x_z_z_xxzzzzz[k];

                g_x_z_zz_xyyyyy[k] = g_x_0_z_xyyyyy[k] - g_x_z_z_xyyyyy[k] * ab_z + g_x_z_z_xyyyyyz[k];

                g_x_z_zz_xyyyyz[k] = g_x_0_z_xyyyyz[k] - g_x_z_z_xyyyyz[k] * ab_z + g_x_z_z_xyyyyzz[k];

                g_x_z_zz_xyyyzz[k] = g_x_0_z_xyyyzz[k] - g_x_z_z_xyyyzz[k] * ab_z + g_x_z_z_xyyyzzz[k];

                g_x_z_zz_xyyzzz[k] = g_x_0_z_xyyzzz[k] - g_x_z_z_xyyzzz[k] * ab_z + g_x_z_z_xyyzzzz[k];

                g_x_z_zz_xyzzzz[k] = g_x_0_z_xyzzzz[k] - g_x_z_z_xyzzzz[k] * ab_z + g_x_z_z_xyzzzzz[k];

                g_x_z_zz_xzzzzz[k] = g_x_0_z_xzzzzz[k] - g_x_z_z_xzzzzz[k] * ab_z + g_x_z_z_xzzzzzz[k];

                g_x_z_zz_yyyyyy[k] = g_x_0_z_yyyyyy[k] - g_x_z_z_yyyyyy[k] * ab_z + g_x_z_z_yyyyyyz[k];

                g_x_z_zz_yyyyyz[k] = g_x_0_z_yyyyyz[k] - g_x_z_z_yyyyyz[k] * ab_z + g_x_z_z_yyyyyzz[k];

                g_x_z_zz_yyyyzz[k] = g_x_0_z_yyyyzz[k] - g_x_z_z_yyyyzz[k] * ab_z + g_x_z_z_yyyyzzz[k];

                g_x_z_zz_yyyzzz[k] = g_x_0_z_yyyzzz[k] - g_x_z_z_yyyzzz[k] * ab_z + g_x_z_z_yyyzzzz[k];

                g_x_z_zz_yyzzzz[k] = g_x_0_z_yyzzzz[k] - g_x_z_z_yyzzzz[k] * ab_z + g_x_z_z_yyzzzzz[k];

                g_x_z_zz_yzzzzz[k] = g_x_0_z_yzzzzz[k] - g_x_z_z_yzzzzz[k] * ab_z + g_x_z_z_yzzzzzz[k];

                g_x_z_zz_zzzzzz[k] = g_x_0_z_zzzzzz[k] - g_x_z_z_zzzzzz[k] * ab_z + g_x_z_z_zzzzzzz[k];
            }

            /// Set up 504-532 components of targeted buffer : cbuffer.data(

            auto g_y_x_xx_xxxxxx = cbuffer.data(di_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_x_xx_xxxxxy = cbuffer.data(di_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_x_xx_xxxxxz = cbuffer.data(di_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_x_xx_xxxxyy = cbuffer.data(di_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_x_xx_xxxxyz = cbuffer.data(di_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_x_xx_xxxxzz = cbuffer.data(di_geom_11_off + 509 * ccomps * dcomps);

            auto g_y_x_xx_xxxyyy = cbuffer.data(di_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_x_xx_xxxyyz = cbuffer.data(di_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_x_xx_xxxyzz = cbuffer.data(di_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_x_xx_xxxzzz = cbuffer.data(di_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_x_xx_xxyyyy = cbuffer.data(di_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_x_xx_xxyyyz = cbuffer.data(di_geom_11_off + 515 * ccomps * dcomps);

            auto g_y_x_xx_xxyyzz = cbuffer.data(di_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_x_xx_xxyzzz = cbuffer.data(di_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_x_xx_xxzzzz = cbuffer.data(di_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_x_xx_xyyyyy = cbuffer.data(di_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_x_xx_xyyyyz = cbuffer.data(di_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_x_xx_xyyyzz = cbuffer.data(di_geom_11_off + 521 * ccomps * dcomps);

            auto g_y_x_xx_xyyzzz = cbuffer.data(di_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_x_xx_xyzzzz = cbuffer.data(di_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_x_xx_xzzzzz = cbuffer.data(di_geom_11_off + 524 * ccomps * dcomps);

            auto g_y_x_xx_yyyyyy = cbuffer.data(di_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_x_xx_yyyyyz = cbuffer.data(di_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_x_xx_yyyyzz = cbuffer.data(di_geom_11_off + 527 * ccomps * dcomps);

            auto g_y_x_xx_yyyzzz = cbuffer.data(di_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_x_xx_yyzzzz = cbuffer.data(di_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_x_xx_yzzzzz = cbuffer.data(di_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_x_xx_zzzzzz = cbuffer.data(di_geom_11_off + 531 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_xxxxxx, g_y_0_x_xxxxxy, g_y_0_x_xxxxxz, g_y_0_x_xxxxyy, g_y_0_x_xxxxyz, g_y_0_x_xxxxzz, g_y_0_x_xxxyyy, g_y_0_x_xxxyyz, g_y_0_x_xxxyzz, g_y_0_x_xxxzzz, g_y_0_x_xxyyyy, g_y_0_x_xxyyyz, g_y_0_x_xxyyzz, g_y_0_x_xxyzzz, g_y_0_x_xxzzzz, g_y_0_x_xyyyyy, g_y_0_x_xyyyyz, g_y_0_x_xyyyzz, g_y_0_x_xyyzzz, g_y_0_x_xyzzzz, g_y_0_x_xzzzzz, g_y_0_x_yyyyyy, g_y_0_x_yyyyyz, g_y_0_x_yyyyzz, g_y_0_x_yyyzzz, g_y_0_x_yyzzzz, g_y_0_x_yzzzzz, g_y_0_x_zzzzzz, g_y_x_x_xxxxxx, g_y_x_x_xxxxxxx, g_y_x_x_xxxxxxy, g_y_x_x_xxxxxxz, g_y_x_x_xxxxxy, g_y_x_x_xxxxxyy, g_y_x_x_xxxxxyz, g_y_x_x_xxxxxz, g_y_x_x_xxxxxzz, g_y_x_x_xxxxyy, g_y_x_x_xxxxyyy, g_y_x_x_xxxxyyz, g_y_x_x_xxxxyz, g_y_x_x_xxxxyzz, g_y_x_x_xxxxzz, g_y_x_x_xxxxzzz, g_y_x_x_xxxyyy, g_y_x_x_xxxyyyy, g_y_x_x_xxxyyyz, g_y_x_x_xxxyyz, g_y_x_x_xxxyyzz, g_y_x_x_xxxyzz, g_y_x_x_xxxyzzz, g_y_x_x_xxxzzz, g_y_x_x_xxxzzzz, g_y_x_x_xxyyyy, g_y_x_x_xxyyyyy, g_y_x_x_xxyyyyz, g_y_x_x_xxyyyz, g_y_x_x_xxyyyzz, g_y_x_x_xxyyzz, g_y_x_x_xxyyzzz, g_y_x_x_xxyzzz, g_y_x_x_xxyzzzz, g_y_x_x_xxzzzz, g_y_x_x_xxzzzzz, g_y_x_x_xyyyyy, g_y_x_x_xyyyyyy, g_y_x_x_xyyyyyz, g_y_x_x_xyyyyz, g_y_x_x_xyyyyzz, g_y_x_x_xyyyzz, g_y_x_x_xyyyzzz, g_y_x_x_xyyzzz, g_y_x_x_xyyzzzz, g_y_x_x_xyzzzz, g_y_x_x_xyzzzzz, g_y_x_x_xzzzzz, g_y_x_x_xzzzzzz, g_y_x_x_yyyyyy, g_y_x_x_yyyyyz, g_y_x_x_yyyyzz, g_y_x_x_yyyzzz, g_y_x_x_yyzzzz, g_y_x_x_yzzzzz, g_y_x_x_zzzzzz, g_y_x_xx_xxxxxx, g_y_x_xx_xxxxxy, g_y_x_xx_xxxxxz, g_y_x_xx_xxxxyy, g_y_x_xx_xxxxyz, g_y_x_xx_xxxxzz, g_y_x_xx_xxxyyy, g_y_x_xx_xxxyyz, g_y_x_xx_xxxyzz, g_y_x_xx_xxxzzz, g_y_x_xx_xxyyyy, g_y_x_xx_xxyyyz, g_y_x_xx_xxyyzz, g_y_x_xx_xxyzzz, g_y_x_xx_xxzzzz, g_y_x_xx_xyyyyy, g_y_x_xx_xyyyyz, g_y_x_xx_xyyyzz, g_y_x_xx_xyyzzz, g_y_x_xx_xyzzzz, g_y_x_xx_xzzzzz, g_y_x_xx_yyyyyy, g_y_x_xx_yyyyyz, g_y_x_xx_yyyyzz, g_y_x_xx_yyyzzz, g_y_x_xx_yyzzzz, g_y_x_xx_yzzzzz, g_y_x_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xx_xxxxxx[k] = g_y_0_x_xxxxxx[k] - g_y_x_x_xxxxxx[k] * ab_x + g_y_x_x_xxxxxxx[k];

                g_y_x_xx_xxxxxy[k] = g_y_0_x_xxxxxy[k] - g_y_x_x_xxxxxy[k] * ab_x + g_y_x_x_xxxxxxy[k];

                g_y_x_xx_xxxxxz[k] = g_y_0_x_xxxxxz[k] - g_y_x_x_xxxxxz[k] * ab_x + g_y_x_x_xxxxxxz[k];

                g_y_x_xx_xxxxyy[k] = g_y_0_x_xxxxyy[k] - g_y_x_x_xxxxyy[k] * ab_x + g_y_x_x_xxxxxyy[k];

                g_y_x_xx_xxxxyz[k] = g_y_0_x_xxxxyz[k] - g_y_x_x_xxxxyz[k] * ab_x + g_y_x_x_xxxxxyz[k];

                g_y_x_xx_xxxxzz[k] = g_y_0_x_xxxxzz[k] - g_y_x_x_xxxxzz[k] * ab_x + g_y_x_x_xxxxxzz[k];

                g_y_x_xx_xxxyyy[k] = g_y_0_x_xxxyyy[k] - g_y_x_x_xxxyyy[k] * ab_x + g_y_x_x_xxxxyyy[k];

                g_y_x_xx_xxxyyz[k] = g_y_0_x_xxxyyz[k] - g_y_x_x_xxxyyz[k] * ab_x + g_y_x_x_xxxxyyz[k];

                g_y_x_xx_xxxyzz[k] = g_y_0_x_xxxyzz[k] - g_y_x_x_xxxyzz[k] * ab_x + g_y_x_x_xxxxyzz[k];

                g_y_x_xx_xxxzzz[k] = g_y_0_x_xxxzzz[k] - g_y_x_x_xxxzzz[k] * ab_x + g_y_x_x_xxxxzzz[k];

                g_y_x_xx_xxyyyy[k] = g_y_0_x_xxyyyy[k] - g_y_x_x_xxyyyy[k] * ab_x + g_y_x_x_xxxyyyy[k];

                g_y_x_xx_xxyyyz[k] = g_y_0_x_xxyyyz[k] - g_y_x_x_xxyyyz[k] * ab_x + g_y_x_x_xxxyyyz[k];

                g_y_x_xx_xxyyzz[k] = g_y_0_x_xxyyzz[k] - g_y_x_x_xxyyzz[k] * ab_x + g_y_x_x_xxxyyzz[k];

                g_y_x_xx_xxyzzz[k] = g_y_0_x_xxyzzz[k] - g_y_x_x_xxyzzz[k] * ab_x + g_y_x_x_xxxyzzz[k];

                g_y_x_xx_xxzzzz[k] = g_y_0_x_xxzzzz[k] - g_y_x_x_xxzzzz[k] * ab_x + g_y_x_x_xxxzzzz[k];

                g_y_x_xx_xyyyyy[k] = g_y_0_x_xyyyyy[k] - g_y_x_x_xyyyyy[k] * ab_x + g_y_x_x_xxyyyyy[k];

                g_y_x_xx_xyyyyz[k] = g_y_0_x_xyyyyz[k] - g_y_x_x_xyyyyz[k] * ab_x + g_y_x_x_xxyyyyz[k];

                g_y_x_xx_xyyyzz[k] = g_y_0_x_xyyyzz[k] - g_y_x_x_xyyyzz[k] * ab_x + g_y_x_x_xxyyyzz[k];

                g_y_x_xx_xyyzzz[k] = g_y_0_x_xyyzzz[k] - g_y_x_x_xyyzzz[k] * ab_x + g_y_x_x_xxyyzzz[k];

                g_y_x_xx_xyzzzz[k] = g_y_0_x_xyzzzz[k] - g_y_x_x_xyzzzz[k] * ab_x + g_y_x_x_xxyzzzz[k];

                g_y_x_xx_xzzzzz[k] = g_y_0_x_xzzzzz[k] - g_y_x_x_xzzzzz[k] * ab_x + g_y_x_x_xxzzzzz[k];

                g_y_x_xx_yyyyyy[k] = g_y_0_x_yyyyyy[k] - g_y_x_x_yyyyyy[k] * ab_x + g_y_x_x_xyyyyyy[k];

                g_y_x_xx_yyyyyz[k] = g_y_0_x_yyyyyz[k] - g_y_x_x_yyyyyz[k] * ab_x + g_y_x_x_xyyyyyz[k];

                g_y_x_xx_yyyyzz[k] = g_y_0_x_yyyyzz[k] - g_y_x_x_yyyyzz[k] * ab_x + g_y_x_x_xyyyyzz[k];

                g_y_x_xx_yyyzzz[k] = g_y_0_x_yyyzzz[k] - g_y_x_x_yyyzzz[k] * ab_x + g_y_x_x_xyyyzzz[k];

                g_y_x_xx_yyzzzz[k] = g_y_0_x_yyzzzz[k] - g_y_x_x_yyzzzz[k] * ab_x + g_y_x_x_xyyzzzz[k];

                g_y_x_xx_yzzzzz[k] = g_y_0_x_yzzzzz[k] - g_y_x_x_yzzzzz[k] * ab_x + g_y_x_x_xyzzzzz[k];

                g_y_x_xx_zzzzzz[k] = g_y_0_x_zzzzzz[k] - g_y_x_x_zzzzzz[k] * ab_x + g_y_x_x_xzzzzzz[k];
            }

            /// Set up 532-560 components of targeted buffer : cbuffer.data(

            auto g_y_x_xy_xxxxxx = cbuffer.data(di_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_x_xy_xxxxxy = cbuffer.data(di_geom_11_off + 533 * ccomps * dcomps);

            auto g_y_x_xy_xxxxxz = cbuffer.data(di_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_x_xy_xxxxyy = cbuffer.data(di_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_x_xy_xxxxyz = cbuffer.data(di_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_x_xy_xxxxzz = cbuffer.data(di_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_x_xy_xxxyyy = cbuffer.data(di_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_x_xy_xxxyyz = cbuffer.data(di_geom_11_off + 539 * ccomps * dcomps);

            auto g_y_x_xy_xxxyzz = cbuffer.data(di_geom_11_off + 540 * ccomps * dcomps);

            auto g_y_x_xy_xxxzzz = cbuffer.data(di_geom_11_off + 541 * ccomps * dcomps);

            auto g_y_x_xy_xxyyyy = cbuffer.data(di_geom_11_off + 542 * ccomps * dcomps);

            auto g_y_x_xy_xxyyyz = cbuffer.data(di_geom_11_off + 543 * ccomps * dcomps);

            auto g_y_x_xy_xxyyzz = cbuffer.data(di_geom_11_off + 544 * ccomps * dcomps);

            auto g_y_x_xy_xxyzzz = cbuffer.data(di_geom_11_off + 545 * ccomps * dcomps);

            auto g_y_x_xy_xxzzzz = cbuffer.data(di_geom_11_off + 546 * ccomps * dcomps);

            auto g_y_x_xy_xyyyyy = cbuffer.data(di_geom_11_off + 547 * ccomps * dcomps);

            auto g_y_x_xy_xyyyyz = cbuffer.data(di_geom_11_off + 548 * ccomps * dcomps);

            auto g_y_x_xy_xyyyzz = cbuffer.data(di_geom_11_off + 549 * ccomps * dcomps);

            auto g_y_x_xy_xyyzzz = cbuffer.data(di_geom_11_off + 550 * ccomps * dcomps);

            auto g_y_x_xy_xyzzzz = cbuffer.data(di_geom_11_off + 551 * ccomps * dcomps);

            auto g_y_x_xy_xzzzzz = cbuffer.data(di_geom_11_off + 552 * ccomps * dcomps);

            auto g_y_x_xy_yyyyyy = cbuffer.data(di_geom_11_off + 553 * ccomps * dcomps);

            auto g_y_x_xy_yyyyyz = cbuffer.data(di_geom_11_off + 554 * ccomps * dcomps);

            auto g_y_x_xy_yyyyzz = cbuffer.data(di_geom_11_off + 555 * ccomps * dcomps);

            auto g_y_x_xy_yyyzzz = cbuffer.data(di_geom_11_off + 556 * ccomps * dcomps);

            auto g_y_x_xy_yyzzzz = cbuffer.data(di_geom_11_off + 557 * ccomps * dcomps);

            auto g_y_x_xy_yzzzzz = cbuffer.data(di_geom_11_off + 558 * ccomps * dcomps);

            auto g_y_x_xy_zzzzzz = cbuffer.data(di_geom_11_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxxxxx, g_y_0_y_xxxxxy, g_y_0_y_xxxxxz, g_y_0_y_xxxxyy, g_y_0_y_xxxxyz, g_y_0_y_xxxxzz, g_y_0_y_xxxyyy, g_y_0_y_xxxyyz, g_y_0_y_xxxyzz, g_y_0_y_xxxzzz, g_y_0_y_xxyyyy, g_y_0_y_xxyyyz, g_y_0_y_xxyyzz, g_y_0_y_xxyzzz, g_y_0_y_xxzzzz, g_y_0_y_xyyyyy, g_y_0_y_xyyyyz, g_y_0_y_xyyyzz, g_y_0_y_xyyzzz, g_y_0_y_xyzzzz, g_y_0_y_xzzzzz, g_y_0_y_yyyyyy, g_y_0_y_yyyyyz, g_y_0_y_yyyyzz, g_y_0_y_yyyzzz, g_y_0_y_yyzzzz, g_y_0_y_yzzzzz, g_y_0_y_zzzzzz, g_y_x_xy_xxxxxx, g_y_x_xy_xxxxxy, g_y_x_xy_xxxxxz, g_y_x_xy_xxxxyy, g_y_x_xy_xxxxyz, g_y_x_xy_xxxxzz, g_y_x_xy_xxxyyy, g_y_x_xy_xxxyyz, g_y_x_xy_xxxyzz, g_y_x_xy_xxxzzz, g_y_x_xy_xxyyyy, g_y_x_xy_xxyyyz, g_y_x_xy_xxyyzz, g_y_x_xy_xxyzzz, g_y_x_xy_xxzzzz, g_y_x_xy_xyyyyy, g_y_x_xy_xyyyyz, g_y_x_xy_xyyyzz, g_y_x_xy_xyyzzz, g_y_x_xy_xyzzzz, g_y_x_xy_xzzzzz, g_y_x_xy_yyyyyy, g_y_x_xy_yyyyyz, g_y_x_xy_yyyyzz, g_y_x_xy_yyyzzz, g_y_x_xy_yyzzzz, g_y_x_xy_yzzzzz, g_y_x_xy_zzzzzz, g_y_x_y_xxxxxx, g_y_x_y_xxxxxxx, g_y_x_y_xxxxxxy, g_y_x_y_xxxxxxz, g_y_x_y_xxxxxy, g_y_x_y_xxxxxyy, g_y_x_y_xxxxxyz, g_y_x_y_xxxxxz, g_y_x_y_xxxxxzz, g_y_x_y_xxxxyy, g_y_x_y_xxxxyyy, g_y_x_y_xxxxyyz, g_y_x_y_xxxxyz, g_y_x_y_xxxxyzz, g_y_x_y_xxxxzz, g_y_x_y_xxxxzzz, g_y_x_y_xxxyyy, g_y_x_y_xxxyyyy, g_y_x_y_xxxyyyz, g_y_x_y_xxxyyz, g_y_x_y_xxxyyzz, g_y_x_y_xxxyzz, g_y_x_y_xxxyzzz, g_y_x_y_xxxzzz, g_y_x_y_xxxzzzz, g_y_x_y_xxyyyy, g_y_x_y_xxyyyyy, g_y_x_y_xxyyyyz, g_y_x_y_xxyyyz, g_y_x_y_xxyyyzz, g_y_x_y_xxyyzz, g_y_x_y_xxyyzzz, g_y_x_y_xxyzzz, g_y_x_y_xxyzzzz, g_y_x_y_xxzzzz, g_y_x_y_xxzzzzz, g_y_x_y_xyyyyy, g_y_x_y_xyyyyyy, g_y_x_y_xyyyyyz, g_y_x_y_xyyyyz, g_y_x_y_xyyyyzz, g_y_x_y_xyyyzz, g_y_x_y_xyyyzzz, g_y_x_y_xyyzzz, g_y_x_y_xyyzzzz, g_y_x_y_xyzzzz, g_y_x_y_xyzzzzz, g_y_x_y_xzzzzz, g_y_x_y_xzzzzzz, g_y_x_y_yyyyyy, g_y_x_y_yyyyyz, g_y_x_y_yyyyzz, g_y_x_y_yyyzzz, g_y_x_y_yyzzzz, g_y_x_y_yzzzzz, g_y_x_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xy_xxxxxx[k] = g_y_0_y_xxxxxx[k] - g_y_x_y_xxxxxx[k] * ab_x + g_y_x_y_xxxxxxx[k];

                g_y_x_xy_xxxxxy[k] = g_y_0_y_xxxxxy[k] - g_y_x_y_xxxxxy[k] * ab_x + g_y_x_y_xxxxxxy[k];

                g_y_x_xy_xxxxxz[k] = g_y_0_y_xxxxxz[k] - g_y_x_y_xxxxxz[k] * ab_x + g_y_x_y_xxxxxxz[k];

                g_y_x_xy_xxxxyy[k] = g_y_0_y_xxxxyy[k] - g_y_x_y_xxxxyy[k] * ab_x + g_y_x_y_xxxxxyy[k];

                g_y_x_xy_xxxxyz[k] = g_y_0_y_xxxxyz[k] - g_y_x_y_xxxxyz[k] * ab_x + g_y_x_y_xxxxxyz[k];

                g_y_x_xy_xxxxzz[k] = g_y_0_y_xxxxzz[k] - g_y_x_y_xxxxzz[k] * ab_x + g_y_x_y_xxxxxzz[k];

                g_y_x_xy_xxxyyy[k] = g_y_0_y_xxxyyy[k] - g_y_x_y_xxxyyy[k] * ab_x + g_y_x_y_xxxxyyy[k];

                g_y_x_xy_xxxyyz[k] = g_y_0_y_xxxyyz[k] - g_y_x_y_xxxyyz[k] * ab_x + g_y_x_y_xxxxyyz[k];

                g_y_x_xy_xxxyzz[k] = g_y_0_y_xxxyzz[k] - g_y_x_y_xxxyzz[k] * ab_x + g_y_x_y_xxxxyzz[k];

                g_y_x_xy_xxxzzz[k] = g_y_0_y_xxxzzz[k] - g_y_x_y_xxxzzz[k] * ab_x + g_y_x_y_xxxxzzz[k];

                g_y_x_xy_xxyyyy[k] = g_y_0_y_xxyyyy[k] - g_y_x_y_xxyyyy[k] * ab_x + g_y_x_y_xxxyyyy[k];

                g_y_x_xy_xxyyyz[k] = g_y_0_y_xxyyyz[k] - g_y_x_y_xxyyyz[k] * ab_x + g_y_x_y_xxxyyyz[k];

                g_y_x_xy_xxyyzz[k] = g_y_0_y_xxyyzz[k] - g_y_x_y_xxyyzz[k] * ab_x + g_y_x_y_xxxyyzz[k];

                g_y_x_xy_xxyzzz[k] = g_y_0_y_xxyzzz[k] - g_y_x_y_xxyzzz[k] * ab_x + g_y_x_y_xxxyzzz[k];

                g_y_x_xy_xxzzzz[k] = g_y_0_y_xxzzzz[k] - g_y_x_y_xxzzzz[k] * ab_x + g_y_x_y_xxxzzzz[k];

                g_y_x_xy_xyyyyy[k] = g_y_0_y_xyyyyy[k] - g_y_x_y_xyyyyy[k] * ab_x + g_y_x_y_xxyyyyy[k];

                g_y_x_xy_xyyyyz[k] = g_y_0_y_xyyyyz[k] - g_y_x_y_xyyyyz[k] * ab_x + g_y_x_y_xxyyyyz[k];

                g_y_x_xy_xyyyzz[k] = g_y_0_y_xyyyzz[k] - g_y_x_y_xyyyzz[k] * ab_x + g_y_x_y_xxyyyzz[k];

                g_y_x_xy_xyyzzz[k] = g_y_0_y_xyyzzz[k] - g_y_x_y_xyyzzz[k] * ab_x + g_y_x_y_xxyyzzz[k];

                g_y_x_xy_xyzzzz[k] = g_y_0_y_xyzzzz[k] - g_y_x_y_xyzzzz[k] * ab_x + g_y_x_y_xxyzzzz[k];

                g_y_x_xy_xzzzzz[k] = g_y_0_y_xzzzzz[k] - g_y_x_y_xzzzzz[k] * ab_x + g_y_x_y_xxzzzzz[k];

                g_y_x_xy_yyyyyy[k] = g_y_0_y_yyyyyy[k] - g_y_x_y_yyyyyy[k] * ab_x + g_y_x_y_xyyyyyy[k];

                g_y_x_xy_yyyyyz[k] = g_y_0_y_yyyyyz[k] - g_y_x_y_yyyyyz[k] * ab_x + g_y_x_y_xyyyyyz[k];

                g_y_x_xy_yyyyzz[k] = g_y_0_y_yyyyzz[k] - g_y_x_y_yyyyzz[k] * ab_x + g_y_x_y_xyyyyzz[k];

                g_y_x_xy_yyyzzz[k] = g_y_0_y_yyyzzz[k] - g_y_x_y_yyyzzz[k] * ab_x + g_y_x_y_xyyyzzz[k];

                g_y_x_xy_yyzzzz[k] = g_y_0_y_yyzzzz[k] - g_y_x_y_yyzzzz[k] * ab_x + g_y_x_y_xyyzzzz[k];

                g_y_x_xy_yzzzzz[k] = g_y_0_y_yzzzzz[k] - g_y_x_y_yzzzzz[k] * ab_x + g_y_x_y_xyzzzzz[k];

                g_y_x_xy_zzzzzz[k] = g_y_0_y_zzzzzz[k] - g_y_x_y_zzzzzz[k] * ab_x + g_y_x_y_xzzzzzz[k];
            }

            /// Set up 560-588 components of targeted buffer : cbuffer.data(

            auto g_y_x_xz_xxxxxx = cbuffer.data(di_geom_11_off + 560 * ccomps * dcomps);

            auto g_y_x_xz_xxxxxy = cbuffer.data(di_geom_11_off + 561 * ccomps * dcomps);

            auto g_y_x_xz_xxxxxz = cbuffer.data(di_geom_11_off + 562 * ccomps * dcomps);

            auto g_y_x_xz_xxxxyy = cbuffer.data(di_geom_11_off + 563 * ccomps * dcomps);

            auto g_y_x_xz_xxxxyz = cbuffer.data(di_geom_11_off + 564 * ccomps * dcomps);

            auto g_y_x_xz_xxxxzz = cbuffer.data(di_geom_11_off + 565 * ccomps * dcomps);

            auto g_y_x_xz_xxxyyy = cbuffer.data(di_geom_11_off + 566 * ccomps * dcomps);

            auto g_y_x_xz_xxxyyz = cbuffer.data(di_geom_11_off + 567 * ccomps * dcomps);

            auto g_y_x_xz_xxxyzz = cbuffer.data(di_geom_11_off + 568 * ccomps * dcomps);

            auto g_y_x_xz_xxxzzz = cbuffer.data(di_geom_11_off + 569 * ccomps * dcomps);

            auto g_y_x_xz_xxyyyy = cbuffer.data(di_geom_11_off + 570 * ccomps * dcomps);

            auto g_y_x_xz_xxyyyz = cbuffer.data(di_geom_11_off + 571 * ccomps * dcomps);

            auto g_y_x_xz_xxyyzz = cbuffer.data(di_geom_11_off + 572 * ccomps * dcomps);

            auto g_y_x_xz_xxyzzz = cbuffer.data(di_geom_11_off + 573 * ccomps * dcomps);

            auto g_y_x_xz_xxzzzz = cbuffer.data(di_geom_11_off + 574 * ccomps * dcomps);

            auto g_y_x_xz_xyyyyy = cbuffer.data(di_geom_11_off + 575 * ccomps * dcomps);

            auto g_y_x_xz_xyyyyz = cbuffer.data(di_geom_11_off + 576 * ccomps * dcomps);

            auto g_y_x_xz_xyyyzz = cbuffer.data(di_geom_11_off + 577 * ccomps * dcomps);

            auto g_y_x_xz_xyyzzz = cbuffer.data(di_geom_11_off + 578 * ccomps * dcomps);

            auto g_y_x_xz_xyzzzz = cbuffer.data(di_geom_11_off + 579 * ccomps * dcomps);

            auto g_y_x_xz_xzzzzz = cbuffer.data(di_geom_11_off + 580 * ccomps * dcomps);

            auto g_y_x_xz_yyyyyy = cbuffer.data(di_geom_11_off + 581 * ccomps * dcomps);

            auto g_y_x_xz_yyyyyz = cbuffer.data(di_geom_11_off + 582 * ccomps * dcomps);

            auto g_y_x_xz_yyyyzz = cbuffer.data(di_geom_11_off + 583 * ccomps * dcomps);

            auto g_y_x_xz_yyyzzz = cbuffer.data(di_geom_11_off + 584 * ccomps * dcomps);

            auto g_y_x_xz_yyzzzz = cbuffer.data(di_geom_11_off + 585 * ccomps * dcomps);

            auto g_y_x_xz_yzzzzz = cbuffer.data(di_geom_11_off + 586 * ccomps * dcomps);

            auto g_y_x_xz_zzzzzz = cbuffer.data(di_geom_11_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_x_xxxxxx, g_y_x_x_xxxxxxz, g_y_x_x_xxxxxy, g_y_x_x_xxxxxyz, g_y_x_x_xxxxxz, g_y_x_x_xxxxxzz, g_y_x_x_xxxxyy, g_y_x_x_xxxxyyz, g_y_x_x_xxxxyz, g_y_x_x_xxxxyzz, g_y_x_x_xxxxzz, g_y_x_x_xxxxzzz, g_y_x_x_xxxyyy, g_y_x_x_xxxyyyz, g_y_x_x_xxxyyz, g_y_x_x_xxxyyzz, g_y_x_x_xxxyzz, g_y_x_x_xxxyzzz, g_y_x_x_xxxzzz, g_y_x_x_xxxzzzz, g_y_x_x_xxyyyy, g_y_x_x_xxyyyyz, g_y_x_x_xxyyyz, g_y_x_x_xxyyyzz, g_y_x_x_xxyyzz, g_y_x_x_xxyyzzz, g_y_x_x_xxyzzz, g_y_x_x_xxyzzzz, g_y_x_x_xxzzzz, g_y_x_x_xxzzzzz, g_y_x_x_xyyyyy, g_y_x_x_xyyyyyz, g_y_x_x_xyyyyz, g_y_x_x_xyyyyzz, g_y_x_x_xyyyzz, g_y_x_x_xyyyzzz, g_y_x_x_xyyzzz, g_y_x_x_xyyzzzz, g_y_x_x_xyzzzz, g_y_x_x_xyzzzzz, g_y_x_x_xzzzzz, g_y_x_x_xzzzzzz, g_y_x_x_yyyyyy, g_y_x_x_yyyyyyz, g_y_x_x_yyyyyz, g_y_x_x_yyyyyzz, g_y_x_x_yyyyzz, g_y_x_x_yyyyzzz, g_y_x_x_yyyzzz, g_y_x_x_yyyzzzz, g_y_x_x_yyzzzz, g_y_x_x_yyzzzzz, g_y_x_x_yzzzzz, g_y_x_x_yzzzzzz, g_y_x_x_zzzzzz, g_y_x_x_zzzzzzz, g_y_x_xz_xxxxxx, g_y_x_xz_xxxxxy, g_y_x_xz_xxxxxz, g_y_x_xz_xxxxyy, g_y_x_xz_xxxxyz, g_y_x_xz_xxxxzz, g_y_x_xz_xxxyyy, g_y_x_xz_xxxyyz, g_y_x_xz_xxxyzz, g_y_x_xz_xxxzzz, g_y_x_xz_xxyyyy, g_y_x_xz_xxyyyz, g_y_x_xz_xxyyzz, g_y_x_xz_xxyzzz, g_y_x_xz_xxzzzz, g_y_x_xz_xyyyyy, g_y_x_xz_xyyyyz, g_y_x_xz_xyyyzz, g_y_x_xz_xyyzzz, g_y_x_xz_xyzzzz, g_y_x_xz_xzzzzz, g_y_x_xz_yyyyyy, g_y_x_xz_yyyyyz, g_y_x_xz_yyyyzz, g_y_x_xz_yyyzzz, g_y_x_xz_yyzzzz, g_y_x_xz_yzzzzz, g_y_x_xz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xz_xxxxxx[k] = -g_y_x_x_xxxxxx[k] * ab_z + g_y_x_x_xxxxxxz[k];

                g_y_x_xz_xxxxxy[k] = -g_y_x_x_xxxxxy[k] * ab_z + g_y_x_x_xxxxxyz[k];

                g_y_x_xz_xxxxxz[k] = -g_y_x_x_xxxxxz[k] * ab_z + g_y_x_x_xxxxxzz[k];

                g_y_x_xz_xxxxyy[k] = -g_y_x_x_xxxxyy[k] * ab_z + g_y_x_x_xxxxyyz[k];

                g_y_x_xz_xxxxyz[k] = -g_y_x_x_xxxxyz[k] * ab_z + g_y_x_x_xxxxyzz[k];

                g_y_x_xz_xxxxzz[k] = -g_y_x_x_xxxxzz[k] * ab_z + g_y_x_x_xxxxzzz[k];

                g_y_x_xz_xxxyyy[k] = -g_y_x_x_xxxyyy[k] * ab_z + g_y_x_x_xxxyyyz[k];

                g_y_x_xz_xxxyyz[k] = -g_y_x_x_xxxyyz[k] * ab_z + g_y_x_x_xxxyyzz[k];

                g_y_x_xz_xxxyzz[k] = -g_y_x_x_xxxyzz[k] * ab_z + g_y_x_x_xxxyzzz[k];

                g_y_x_xz_xxxzzz[k] = -g_y_x_x_xxxzzz[k] * ab_z + g_y_x_x_xxxzzzz[k];

                g_y_x_xz_xxyyyy[k] = -g_y_x_x_xxyyyy[k] * ab_z + g_y_x_x_xxyyyyz[k];

                g_y_x_xz_xxyyyz[k] = -g_y_x_x_xxyyyz[k] * ab_z + g_y_x_x_xxyyyzz[k];

                g_y_x_xz_xxyyzz[k] = -g_y_x_x_xxyyzz[k] * ab_z + g_y_x_x_xxyyzzz[k];

                g_y_x_xz_xxyzzz[k] = -g_y_x_x_xxyzzz[k] * ab_z + g_y_x_x_xxyzzzz[k];

                g_y_x_xz_xxzzzz[k] = -g_y_x_x_xxzzzz[k] * ab_z + g_y_x_x_xxzzzzz[k];

                g_y_x_xz_xyyyyy[k] = -g_y_x_x_xyyyyy[k] * ab_z + g_y_x_x_xyyyyyz[k];

                g_y_x_xz_xyyyyz[k] = -g_y_x_x_xyyyyz[k] * ab_z + g_y_x_x_xyyyyzz[k];

                g_y_x_xz_xyyyzz[k] = -g_y_x_x_xyyyzz[k] * ab_z + g_y_x_x_xyyyzzz[k];

                g_y_x_xz_xyyzzz[k] = -g_y_x_x_xyyzzz[k] * ab_z + g_y_x_x_xyyzzzz[k];

                g_y_x_xz_xyzzzz[k] = -g_y_x_x_xyzzzz[k] * ab_z + g_y_x_x_xyzzzzz[k];

                g_y_x_xz_xzzzzz[k] = -g_y_x_x_xzzzzz[k] * ab_z + g_y_x_x_xzzzzzz[k];

                g_y_x_xz_yyyyyy[k] = -g_y_x_x_yyyyyy[k] * ab_z + g_y_x_x_yyyyyyz[k];

                g_y_x_xz_yyyyyz[k] = -g_y_x_x_yyyyyz[k] * ab_z + g_y_x_x_yyyyyzz[k];

                g_y_x_xz_yyyyzz[k] = -g_y_x_x_yyyyzz[k] * ab_z + g_y_x_x_yyyyzzz[k];

                g_y_x_xz_yyyzzz[k] = -g_y_x_x_yyyzzz[k] * ab_z + g_y_x_x_yyyzzzz[k];

                g_y_x_xz_yyzzzz[k] = -g_y_x_x_yyzzzz[k] * ab_z + g_y_x_x_yyzzzzz[k];

                g_y_x_xz_yzzzzz[k] = -g_y_x_x_yzzzzz[k] * ab_z + g_y_x_x_yzzzzzz[k];

                g_y_x_xz_zzzzzz[k] = -g_y_x_x_zzzzzz[k] * ab_z + g_y_x_x_zzzzzzz[k];
            }

            /// Set up 588-616 components of targeted buffer : cbuffer.data(

            auto g_y_x_yy_xxxxxx = cbuffer.data(di_geom_11_off + 588 * ccomps * dcomps);

            auto g_y_x_yy_xxxxxy = cbuffer.data(di_geom_11_off + 589 * ccomps * dcomps);

            auto g_y_x_yy_xxxxxz = cbuffer.data(di_geom_11_off + 590 * ccomps * dcomps);

            auto g_y_x_yy_xxxxyy = cbuffer.data(di_geom_11_off + 591 * ccomps * dcomps);

            auto g_y_x_yy_xxxxyz = cbuffer.data(di_geom_11_off + 592 * ccomps * dcomps);

            auto g_y_x_yy_xxxxzz = cbuffer.data(di_geom_11_off + 593 * ccomps * dcomps);

            auto g_y_x_yy_xxxyyy = cbuffer.data(di_geom_11_off + 594 * ccomps * dcomps);

            auto g_y_x_yy_xxxyyz = cbuffer.data(di_geom_11_off + 595 * ccomps * dcomps);

            auto g_y_x_yy_xxxyzz = cbuffer.data(di_geom_11_off + 596 * ccomps * dcomps);

            auto g_y_x_yy_xxxzzz = cbuffer.data(di_geom_11_off + 597 * ccomps * dcomps);

            auto g_y_x_yy_xxyyyy = cbuffer.data(di_geom_11_off + 598 * ccomps * dcomps);

            auto g_y_x_yy_xxyyyz = cbuffer.data(di_geom_11_off + 599 * ccomps * dcomps);

            auto g_y_x_yy_xxyyzz = cbuffer.data(di_geom_11_off + 600 * ccomps * dcomps);

            auto g_y_x_yy_xxyzzz = cbuffer.data(di_geom_11_off + 601 * ccomps * dcomps);

            auto g_y_x_yy_xxzzzz = cbuffer.data(di_geom_11_off + 602 * ccomps * dcomps);

            auto g_y_x_yy_xyyyyy = cbuffer.data(di_geom_11_off + 603 * ccomps * dcomps);

            auto g_y_x_yy_xyyyyz = cbuffer.data(di_geom_11_off + 604 * ccomps * dcomps);

            auto g_y_x_yy_xyyyzz = cbuffer.data(di_geom_11_off + 605 * ccomps * dcomps);

            auto g_y_x_yy_xyyzzz = cbuffer.data(di_geom_11_off + 606 * ccomps * dcomps);

            auto g_y_x_yy_xyzzzz = cbuffer.data(di_geom_11_off + 607 * ccomps * dcomps);

            auto g_y_x_yy_xzzzzz = cbuffer.data(di_geom_11_off + 608 * ccomps * dcomps);

            auto g_y_x_yy_yyyyyy = cbuffer.data(di_geom_11_off + 609 * ccomps * dcomps);

            auto g_y_x_yy_yyyyyz = cbuffer.data(di_geom_11_off + 610 * ccomps * dcomps);

            auto g_y_x_yy_yyyyzz = cbuffer.data(di_geom_11_off + 611 * ccomps * dcomps);

            auto g_y_x_yy_yyyzzz = cbuffer.data(di_geom_11_off + 612 * ccomps * dcomps);

            auto g_y_x_yy_yyzzzz = cbuffer.data(di_geom_11_off + 613 * ccomps * dcomps);

            auto g_y_x_yy_yzzzzz = cbuffer.data(di_geom_11_off + 614 * ccomps * dcomps);

            auto g_y_x_yy_zzzzzz = cbuffer.data(di_geom_11_off + 615 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_y_xxxxxx, g_0_x_y_xxxxxy, g_0_x_y_xxxxxz, g_0_x_y_xxxxyy, g_0_x_y_xxxxyz, g_0_x_y_xxxxzz, g_0_x_y_xxxyyy, g_0_x_y_xxxyyz, g_0_x_y_xxxyzz, g_0_x_y_xxxzzz, g_0_x_y_xxyyyy, g_0_x_y_xxyyyz, g_0_x_y_xxyyzz, g_0_x_y_xxyzzz, g_0_x_y_xxzzzz, g_0_x_y_xyyyyy, g_0_x_y_xyyyyz, g_0_x_y_xyyyzz, g_0_x_y_xyyzzz, g_0_x_y_xyzzzz, g_0_x_y_xzzzzz, g_0_x_y_yyyyyy, g_0_x_y_yyyyyz, g_0_x_y_yyyyzz, g_0_x_y_yyyzzz, g_0_x_y_yyzzzz, g_0_x_y_yzzzzz, g_0_x_y_zzzzzz, g_y_x_y_xxxxxx, g_y_x_y_xxxxxxy, g_y_x_y_xxxxxy, g_y_x_y_xxxxxyy, g_y_x_y_xxxxxyz, g_y_x_y_xxxxxz, g_y_x_y_xxxxyy, g_y_x_y_xxxxyyy, g_y_x_y_xxxxyyz, g_y_x_y_xxxxyz, g_y_x_y_xxxxyzz, g_y_x_y_xxxxzz, g_y_x_y_xxxyyy, g_y_x_y_xxxyyyy, g_y_x_y_xxxyyyz, g_y_x_y_xxxyyz, g_y_x_y_xxxyyzz, g_y_x_y_xxxyzz, g_y_x_y_xxxyzzz, g_y_x_y_xxxzzz, g_y_x_y_xxyyyy, g_y_x_y_xxyyyyy, g_y_x_y_xxyyyyz, g_y_x_y_xxyyyz, g_y_x_y_xxyyyzz, g_y_x_y_xxyyzz, g_y_x_y_xxyyzzz, g_y_x_y_xxyzzz, g_y_x_y_xxyzzzz, g_y_x_y_xxzzzz, g_y_x_y_xyyyyy, g_y_x_y_xyyyyyy, g_y_x_y_xyyyyyz, g_y_x_y_xyyyyz, g_y_x_y_xyyyyzz, g_y_x_y_xyyyzz, g_y_x_y_xyyyzzz, g_y_x_y_xyyzzz, g_y_x_y_xyyzzzz, g_y_x_y_xyzzzz, g_y_x_y_xyzzzzz, g_y_x_y_xzzzzz, g_y_x_y_yyyyyy, g_y_x_y_yyyyyyy, g_y_x_y_yyyyyyz, g_y_x_y_yyyyyz, g_y_x_y_yyyyyzz, g_y_x_y_yyyyzz, g_y_x_y_yyyyzzz, g_y_x_y_yyyzzz, g_y_x_y_yyyzzzz, g_y_x_y_yyzzzz, g_y_x_y_yyzzzzz, g_y_x_y_yzzzzz, g_y_x_y_yzzzzzz, g_y_x_y_zzzzzz, g_y_x_yy_xxxxxx, g_y_x_yy_xxxxxy, g_y_x_yy_xxxxxz, g_y_x_yy_xxxxyy, g_y_x_yy_xxxxyz, g_y_x_yy_xxxxzz, g_y_x_yy_xxxyyy, g_y_x_yy_xxxyyz, g_y_x_yy_xxxyzz, g_y_x_yy_xxxzzz, g_y_x_yy_xxyyyy, g_y_x_yy_xxyyyz, g_y_x_yy_xxyyzz, g_y_x_yy_xxyzzz, g_y_x_yy_xxzzzz, g_y_x_yy_xyyyyy, g_y_x_yy_xyyyyz, g_y_x_yy_xyyyzz, g_y_x_yy_xyyzzz, g_y_x_yy_xyzzzz, g_y_x_yy_xzzzzz, g_y_x_yy_yyyyyy, g_y_x_yy_yyyyyz, g_y_x_yy_yyyyzz, g_y_x_yy_yyyzzz, g_y_x_yy_yyzzzz, g_y_x_yy_yzzzzz, g_y_x_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yy_xxxxxx[k] = -g_0_x_y_xxxxxx[k] - g_y_x_y_xxxxxx[k] * ab_y + g_y_x_y_xxxxxxy[k];

                g_y_x_yy_xxxxxy[k] = -g_0_x_y_xxxxxy[k] - g_y_x_y_xxxxxy[k] * ab_y + g_y_x_y_xxxxxyy[k];

                g_y_x_yy_xxxxxz[k] = -g_0_x_y_xxxxxz[k] - g_y_x_y_xxxxxz[k] * ab_y + g_y_x_y_xxxxxyz[k];

                g_y_x_yy_xxxxyy[k] = -g_0_x_y_xxxxyy[k] - g_y_x_y_xxxxyy[k] * ab_y + g_y_x_y_xxxxyyy[k];

                g_y_x_yy_xxxxyz[k] = -g_0_x_y_xxxxyz[k] - g_y_x_y_xxxxyz[k] * ab_y + g_y_x_y_xxxxyyz[k];

                g_y_x_yy_xxxxzz[k] = -g_0_x_y_xxxxzz[k] - g_y_x_y_xxxxzz[k] * ab_y + g_y_x_y_xxxxyzz[k];

                g_y_x_yy_xxxyyy[k] = -g_0_x_y_xxxyyy[k] - g_y_x_y_xxxyyy[k] * ab_y + g_y_x_y_xxxyyyy[k];

                g_y_x_yy_xxxyyz[k] = -g_0_x_y_xxxyyz[k] - g_y_x_y_xxxyyz[k] * ab_y + g_y_x_y_xxxyyyz[k];

                g_y_x_yy_xxxyzz[k] = -g_0_x_y_xxxyzz[k] - g_y_x_y_xxxyzz[k] * ab_y + g_y_x_y_xxxyyzz[k];

                g_y_x_yy_xxxzzz[k] = -g_0_x_y_xxxzzz[k] - g_y_x_y_xxxzzz[k] * ab_y + g_y_x_y_xxxyzzz[k];

                g_y_x_yy_xxyyyy[k] = -g_0_x_y_xxyyyy[k] - g_y_x_y_xxyyyy[k] * ab_y + g_y_x_y_xxyyyyy[k];

                g_y_x_yy_xxyyyz[k] = -g_0_x_y_xxyyyz[k] - g_y_x_y_xxyyyz[k] * ab_y + g_y_x_y_xxyyyyz[k];

                g_y_x_yy_xxyyzz[k] = -g_0_x_y_xxyyzz[k] - g_y_x_y_xxyyzz[k] * ab_y + g_y_x_y_xxyyyzz[k];

                g_y_x_yy_xxyzzz[k] = -g_0_x_y_xxyzzz[k] - g_y_x_y_xxyzzz[k] * ab_y + g_y_x_y_xxyyzzz[k];

                g_y_x_yy_xxzzzz[k] = -g_0_x_y_xxzzzz[k] - g_y_x_y_xxzzzz[k] * ab_y + g_y_x_y_xxyzzzz[k];

                g_y_x_yy_xyyyyy[k] = -g_0_x_y_xyyyyy[k] - g_y_x_y_xyyyyy[k] * ab_y + g_y_x_y_xyyyyyy[k];

                g_y_x_yy_xyyyyz[k] = -g_0_x_y_xyyyyz[k] - g_y_x_y_xyyyyz[k] * ab_y + g_y_x_y_xyyyyyz[k];

                g_y_x_yy_xyyyzz[k] = -g_0_x_y_xyyyzz[k] - g_y_x_y_xyyyzz[k] * ab_y + g_y_x_y_xyyyyzz[k];

                g_y_x_yy_xyyzzz[k] = -g_0_x_y_xyyzzz[k] - g_y_x_y_xyyzzz[k] * ab_y + g_y_x_y_xyyyzzz[k];

                g_y_x_yy_xyzzzz[k] = -g_0_x_y_xyzzzz[k] - g_y_x_y_xyzzzz[k] * ab_y + g_y_x_y_xyyzzzz[k];

                g_y_x_yy_xzzzzz[k] = -g_0_x_y_xzzzzz[k] - g_y_x_y_xzzzzz[k] * ab_y + g_y_x_y_xyzzzzz[k];

                g_y_x_yy_yyyyyy[k] = -g_0_x_y_yyyyyy[k] - g_y_x_y_yyyyyy[k] * ab_y + g_y_x_y_yyyyyyy[k];

                g_y_x_yy_yyyyyz[k] = -g_0_x_y_yyyyyz[k] - g_y_x_y_yyyyyz[k] * ab_y + g_y_x_y_yyyyyyz[k];

                g_y_x_yy_yyyyzz[k] = -g_0_x_y_yyyyzz[k] - g_y_x_y_yyyyzz[k] * ab_y + g_y_x_y_yyyyyzz[k];

                g_y_x_yy_yyyzzz[k] = -g_0_x_y_yyyzzz[k] - g_y_x_y_yyyzzz[k] * ab_y + g_y_x_y_yyyyzzz[k];

                g_y_x_yy_yyzzzz[k] = -g_0_x_y_yyzzzz[k] - g_y_x_y_yyzzzz[k] * ab_y + g_y_x_y_yyyzzzz[k];

                g_y_x_yy_yzzzzz[k] = -g_0_x_y_yzzzzz[k] - g_y_x_y_yzzzzz[k] * ab_y + g_y_x_y_yyzzzzz[k];

                g_y_x_yy_zzzzzz[k] = -g_0_x_y_zzzzzz[k] - g_y_x_y_zzzzzz[k] * ab_y + g_y_x_y_yzzzzzz[k];
            }

            /// Set up 616-644 components of targeted buffer : cbuffer.data(

            auto g_y_x_yz_xxxxxx = cbuffer.data(di_geom_11_off + 616 * ccomps * dcomps);

            auto g_y_x_yz_xxxxxy = cbuffer.data(di_geom_11_off + 617 * ccomps * dcomps);

            auto g_y_x_yz_xxxxxz = cbuffer.data(di_geom_11_off + 618 * ccomps * dcomps);

            auto g_y_x_yz_xxxxyy = cbuffer.data(di_geom_11_off + 619 * ccomps * dcomps);

            auto g_y_x_yz_xxxxyz = cbuffer.data(di_geom_11_off + 620 * ccomps * dcomps);

            auto g_y_x_yz_xxxxzz = cbuffer.data(di_geom_11_off + 621 * ccomps * dcomps);

            auto g_y_x_yz_xxxyyy = cbuffer.data(di_geom_11_off + 622 * ccomps * dcomps);

            auto g_y_x_yz_xxxyyz = cbuffer.data(di_geom_11_off + 623 * ccomps * dcomps);

            auto g_y_x_yz_xxxyzz = cbuffer.data(di_geom_11_off + 624 * ccomps * dcomps);

            auto g_y_x_yz_xxxzzz = cbuffer.data(di_geom_11_off + 625 * ccomps * dcomps);

            auto g_y_x_yz_xxyyyy = cbuffer.data(di_geom_11_off + 626 * ccomps * dcomps);

            auto g_y_x_yz_xxyyyz = cbuffer.data(di_geom_11_off + 627 * ccomps * dcomps);

            auto g_y_x_yz_xxyyzz = cbuffer.data(di_geom_11_off + 628 * ccomps * dcomps);

            auto g_y_x_yz_xxyzzz = cbuffer.data(di_geom_11_off + 629 * ccomps * dcomps);

            auto g_y_x_yz_xxzzzz = cbuffer.data(di_geom_11_off + 630 * ccomps * dcomps);

            auto g_y_x_yz_xyyyyy = cbuffer.data(di_geom_11_off + 631 * ccomps * dcomps);

            auto g_y_x_yz_xyyyyz = cbuffer.data(di_geom_11_off + 632 * ccomps * dcomps);

            auto g_y_x_yz_xyyyzz = cbuffer.data(di_geom_11_off + 633 * ccomps * dcomps);

            auto g_y_x_yz_xyyzzz = cbuffer.data(di_geom_11_off + 634 * ccomps * dcomps);

            auto g_y_x_yz_xyzzzz = cbuffer.data(di_geom_11_off + 635 * ccomps * dcomps);

            auto g_y_x_yz_xzzzzz = cbuffer.data(di_geom_11_off + 636 * ccomps * dcomps);

            auto g_y_x_yz_yyyyyy = cbuffer.data(di_geom_11_off + 637 * ccomps * dcomps);

            auto g_y_x_yz_yyyyyz = cbuffer.data(di_geom_11_off + 638 * ccomps * dcomps);

            auto g_y_x_yz_yyyyzz = cbuffer.data(di_geom_11_off + 639 * ccomps * dcomps);

            auto g_y_x_yz_yyyzzz = cbuffer.data(di_geom_11_off + 640 * ccomps * dcomps);

            auto g_y_x_yz_yyzzzz = cbuffer.data(di_geom_11_off + 641 * ccomps * dcomps);

            auto g_y_x_yz_yzzzzz = cbuffer.data(di_geom_11_off + 642 * ccomps * dcomps);

            auto g_y_x_yz_zzzzzz = cbuffer.data(di_geom_11_off + 643 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_y_xxxxxx, g_y_x_y_xxxxxxz, g_y_x_y_xxxxxy, g_y_x_y_xxxxxyz, g_y_x_y_xxxxxz, g_y_x_y_xxxxxzz, g_y_x_y_xxxxyy, g_y_x_y_xxxxyyz, g_y_x_y_xxxxyz, g_y_x_y_xxxxyzz, g_y_x_y_xxxxzz, g_y_x_y_xxxxzzz, g_y_x_y_xxxyyy, g_y_x_y_xxxyyyz, g_y_x_y_xxxyyz, g_y_x_y_xxxyyzz, g_y_x_y_xxxyzz, g_y_x_y_xxxyzzz, g_y_x_y_xxxzzz, g_y_x_y_xxxzzzz, g_y_x_y_xxyyyy, g_y_x_y_xxyyyyz, g_y_x_y_xxyyyz, g_y_x_y_xxyyyzz, g_y_x_y_xxyyzz, g_y_x_y_xxyyzzz, g_y_x_y_xxyzzz, g_y_x_y_xxyzzzz, g_y_x_y_xxzzzz, g_y_x_y_xxzzzzz, g_y_x_y_xyyyyy, g_y_x_y_xyyyyyz, g_y_x_y_xyyyyz, g_y_x_y_xyyyyzz, g_y_x_y_xyyyzz, g_y_x_y_xyyyzzz, g_y_x_y_xyyzzz, g_y_x_y_xyyzzzz, g_y_x_y_xyzzzz, g_y_x_y_xyzzzzz, g_y_x_y_xzzzzz, g_y_x_y_xzzzzzz, g_y_x_y_yyyyyy, g_y_x_y_yyyyyyz, g_y_x_y_yyyyyz, g_y_x_y_yyyyyzz, g_y_x_y_yyyyzz, g_y_x_y_yyyyzzz, g_y_x_y_yyyzzz, g_y_x_y_yyyzzzz, g_y_x_y_yyzzzz, g_y_x_y_yyzzzzz, g_y_x_y_yzzzzz, g_y_x_y_yzzzzzz, g_y_x_y_zzzzzz, g_y_x_y_zzzzzzz, g_y_x_yz_xxxxxx, g_y_x_yz_xxxxxy, g_y_x_yz_xxxxxz, g_y_x_yz_xxxxyy, g_y_x_yz_xxxxyz, g_y_x_yz_xxxxzz, g_y_x_yz_xxxyyy, g_y_x_yz_xxxyyz, g_y_x_yz_xxxyzz, g_y_x_yz_xxxzzz, g_y_x_yz_xxyyyy, g_y_x_yz_xxyyyz, g_y_x_yz_xxyyzz, g_y_x_yz_xxyzzz, g_y_x_yz_xxzzzz, g_y_x_yz_xyyyyy, g_y_x_yz_xyyyyz, g_y_x_yz_xyyyzz, g_y_x_yz_xyyzzz, g_y_x_yz_xyzzzz, g_y_x_yz_xzzzzz, g_y_x_yz_yyyyyy, g_y_x_yz_yyyyyz, g_y_x_yz_yyyyzz, g_y_x_yz_yyyzzz, g_y_x_yz_yyzzzz, g_y_x_yz_yzzzzz, g_y_x_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yz_xxxxxx[k] = -g_y_x_y_xxxxxx[k] * ab_z + g_y_x_y_xxxxxxz[k];

                g_y_x_yz_xxxxxy[k] = -g_y_x_y_xxxxxy[k] * ab_z + g_y_x_y_xxxxxyz[k];

                g_y_x_yz_xxxxxz[k] = -g_y_x_y_xxxxxz[k] * ab_z + g_y_x_y_xxxxxzz[k];

                g_y_x_yz_xxxxyy[k] = -g_y_x_y_xxxxyy[k] * ab_z + g_y_x_y_xxxxyyz[k];

                g_y_x_yz_xxxxyz[k] = -g_y_x_y_xxxxyz[k] * ab_z + g_y_x_y_xxxxyzz[k];

                g_y_x_yz_xxxxzz[k] = -g_y_x_y_xxxxzz[k] * ab_z + g_y_x_y_xxxxzzz[k];

                g_y_x_yz_xxxyyy[k] = -g_y_x_y_xxxyyy[k] * ab_z + g_y_x_y_xxxyyyz[k];

                g_y_x_yz_xxxyyz[k] = -g_y_x_y_xxxyyz[k] * ab_z + g_y_x_y_xxxyyzz[k];

                g_y_x_yz_xxxyzz[k] = -g_y_x_y_xxxyzz[k] * ab_z + g_y_x_y_xxxyzzz[k];

                g_y_x_yz_xxxzzz[k] = -g_y_x_y_xxxzzz[k] * ab_z + g_y_x_y_xxxzzzz[k];

                g_y_x_yz_xxyyyy[k] = -g_y_x_y_xxyyyy[k] * ab_z + g_y_x_y_xxyyyyz[k];

                g_y_x_yz_xxyyyz[k] = -g_y_x_y_xxyyyz[k] * ab_z + g_y_x_y_xxyyyzz[k];

                g_y_x_yz_xxyyzz[k] = -g_y_x_y_xxyyzz[k] * ab_z + g_y_x_y_xxyyzzz[k];

                g_y_x_yz_xxyzzz[k] = -g_y_x_y_xxyzzz[k] * ab_z + g_y_x_y_xxyzzzz[k];

                g_y_x_yz_xxzzzz[k] = -g_y_x_y_xxzzzz[k] * ab_z + g_y_x_y_xxzzzzz[k];

                g_y_x_yz_xyyyyy[k] = -g_y_x_y_xyyyyy[k] * ab_z + g_y_x_y_xyyyyyz[k];

                g_y_x_yz_xyyyyz[k] = -g_y_x_y_xyyyyz[k] * ab_z + g_y_x_y_xyyyyzz[k];

                g_y_x_yz_xyyyzz[k] = -g_y_x_y_xyyyzz[k] * ab_z + g_y_x_y_xyyyzzz[k];

                g_y_x_yz_xyyzzz[k] = -g_y_x_y_xyyzzz[k] * ab_z + g_y_x_y_xyyzzzz[k];

                g_y_x_yz_xyzzzz[k] = -g_y_x_y_xyzzzz[k] * ab_z + g_y_x_y_xyzzzzz[k];

                g_y_x_yz_xzzzzz[k] = -g_y_x_y_xzzzzz[k] * ab_z + g_y_x_y_xzzzzzz[k];

                g_y_x_yz_yyyyyy[k] = -g_y_x_y_yyyyyy[k] * ab_z + g_y_x_y_yyyyyyz[k];

                g_y_x_yz_yyyyyz[k] = -g_y_x_y_yyyyyz[k] * ab_z + g_y_x_y_yyyyyzz[k];

                g_y_x_yz_yyyyzz[k] = -g_y_x_y_yyyyzz[k] * ab_z + g_y_x_y_yyyyzzz[k];

                g_y_x_yz_yyyzzz[k] = -g_y_x_y_yyyzzz[k] * ab_z + g_y_x_y_yyyzzzz[k];

                g_y_x_yz_yyzzzz[k] = -g_y_x_y_yyzzzz[k] * ab_z + g_y_x_y_yyzzzzz[k];

                g_y_x_yz_yzzzzz[k] = -g_y_x_y_yzzzzz[k] * ab_z + g_y_x_y_yzzzzzz[k];

                g_y_x_yz_zzzzzz[k] = -g_y_x_y_zzzzzz[k] * ab_z + g_y_x_y_zzzzzzz[k];
            }

            /// Set up 644-672 components of targeted buffer : cbuffer.data(

            auto g_y_x_zz_xxxxxx = cbuffer.data(di_geom_11_off + 644 * ccomps * dcomps);

            auto g_y_x_zz_xxxxxy = cbuffer.data(di_geom_11_off + 645 * ccomps * dcomps);

            auto g_y_x_zz_xxxxxz = cbuffer.data(di_geom_11_off + 646 * ccomps * dcomps);

            auto g_y_x_zz_xxxxyy = cbuffer.data(di_geom_11_off + 647 * ccomps * dcomps);

            auto g_y_x_zz_xxxxyz = cbuffer.data(di_geom_11_off + 648 * ccomps * dcomps);

            auto g_y_x_zz_xxxxzz = cbuffer.data(di_geom_11_off + 649 * ccomps * dcomps);

            auto g_y_x_zz_xxxyyy = cbuffer.data(di_geom_11_off + 650 * ccomps * dcomps);

            auto g_y_x_zz_xxxyyz = cbuffer.data(di_geom_11_off + 651 * ccomps * dcomps);

            auto g_y_x_zz_xxxyzz = cbuffer.data(di_geom_11_off + 652 * ccomps * dcomps);

            auto g_y_x_zz_xxxzzz = cbuffer.data(di_geom_11_off + 653 * ccomps * dcomps);

            auto g_y_x_zz_xxyyyy = cbuffer.data(di_geom_11_off + 654 * ccomps * dcomps);

            auto g_y_x_zz_xxyyyz = cbuffer.data(di_geom_11_off + 655 * ccomps * dcomps);

            auto g_y_x_zz_xxyyzz = cbuffer.data(di_geom_11_off + 656 * ccomps * dcomps);

            auto g_y_x_zz_xxyzzz = cbuffer.data(di_geom_11_off + 657 * ccomps * dcomps);

            auto g_y_x_zz_xxzzzz = cbuffer.data(di_geom_11_off + 658 * ccomps * dcomps);

            auto g_y_x_zz_xyyyyy = cbuffer.data(di_geom_11_off + 659 * ccomps * dcomps);

            auto g_y_x_zz_xyyyyz = cbuffer.data(di_geom_11_off + 660 * ccomps * dcomps);

            auto g_y_x_zz_xyyyzz = cbuffer.data(di_geom_11_off + 661 * ccomps * dcomps);

            auto g_y_x_zz_xyyzzz = cbuffer.data(di_geom_11_off + 662 * ccomps * dcomps);

            auto g_y_x_zz_xyzzzz = cbuffer.data(di_geom_11_off + 663 * ccomps * dcomps);

            auto g_y_x_zz_xzzzzz = cbuffer.data(di_geom_11_off + 664 * ccomps * dcomps);

            auto g_y_x_zz_yyyyyy = cbuffer.data(di_geom_11_off + 665 * ccomps * dcomps);

            auto g_y_x_zz_yyyyyz = cbuffer.data(di_geom_11_off + 666 * ccomps * dcomps);

            auto g_y_x_zz_yyyyzz = cbuffer.data(di_geom_11_off + 667 * ccomps * dcomps);

            auto g_y_x_zz_yyyzzz = cbuffer.data(di_geom_11_off + 668 * ccomps * dcomps);

            auto g_y_x_zz_yyzzzz = cbuffer.data(di_geom_11_off + 669 * ccomps * dcomps);

            auto g_y_x_zz_yzzzzz = cbuffer.data(di_geom_11_off + 670 * ccomps * dcomps);

            auto g_y_x_zz_zzzzzz = cbuffer.data(di_geom_11_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_z_xxxxxx, g_y_x_z_xxxxxxz, g_y_x_z_xxxxxy, g_y_x_z_xxxxxyz, g_y_x_z_xxxxxz, g_y_x_z_xxxxxzz, g_y_x_z_xxxxyy, g_y_x_z_xxxxyyz, g_y_x_z_xxxxyz, g_y_x_z_xxxxyzz, g_y_x_z_xxxxzz, g_y_x_z_xxxxzzz, g_y_x_z_xxxyyy, g_y_x_z_xxxyyyz, g_y_x_z_xxxyyz, g_y_x_z_xxxyyzz, g_y_x_z_xxxyzz, g_y_x_z_xxxyzzz, g_y_x_z_xxxzzz, g_y_x_z_xxxzzzz, g_y_x_z_xxyyyy, g_y_x_z_xxyyyyz, g_y_x_z_xxyyyz, g_y_x_z_xxyyyzz, g_y_x_z_xxyyzz, g_y_x_z_xxyyzzz, g_y_x_z_xxyzzz, g_y_x_z_xxyzzzz, g_y_x_z_xxzzzz, g_y_x_z_xxzzzzz, g_y_x_z_xyyyyy, g_y_x_z_xyyyyyz, g_y_x_z_xyyyyz, g_y_x_z_xyyyyzz, g_y_x_z_xyyyzz, g_y_x_z_xyyyzzz, g_y_x_z_xyyzzz, g_y_x_z_xyyzzzz, g_y_x_z_xyzzzz, g_y_x_z_xyzzzzz, g_y_x_z_xzzzzz, g_y_x_z_xzzzzzz, g_y_x_z_yyyyyy, g_y_x_z_yyyyyyz, g_y_x_z_yyyyyz, g_y_x_z_yyyyyzz, g_y_x_z_yyyyzz, g_y_x_z_yyyyzzz, g_y_x_z_yyyzzz, g_y_x_z_yyyzzzz, g_y_x_z_yyzzzz, g_y_x_z_yyzzzzz, g_y_x_z_yzzzzz, g_y_x_z_yzzzzzz, g_y_x_z_zzzzzz, g_y_x_z_zzzzzzz, g_y_x_zz_xxxxxx, g_y_x_zz_xxxxxy, g_y_x_zz_xxxxxz, g_y_x_zz_xxxxyy, g_y_x_zz_xxxxyz, g_y_x_zz_xxxxzz, g_y_x_zz_xxxyyy, g_y_x_zz_xxxyyz, g_y_x_zz_xxxyzz, g_y_x_zz_xxxzzz, g_y_x_zz_xxyyyy, g_y_x_zz_xxyyyz, g_y_x_zz_xxyyzz, g_y_x_zz_xxyzzz, g_y_x_zz_xxzzzz, g_y_x_zz_xyyyyy, g_y_x_zz_xyyyyz, g_y_x_zz_xyyyzz, g_y_x_zz_xyyzzz, g_y_x_zz_xyzzzz, g_y_x_zz_xzzzzz, g_y_x_zz_yyyyyy, g_y_x_zz_yyyyyz, g_y_x_zz_yyyyzz, g_y_x_zz_yyyzzz, g_y_x_zz_yyzzzz, g_y_x_zz_yzzzzz, g_y_x_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zz_xxxxxx[k] = -g_y_x_z_xxxxxx[k] * ab_z + g_y_x_z_xxxxxxz[k];

                g_y_x_zz_xxxxxy[k] = -g_y_x_z_xxxxxy[k] * ab_z + g_y_x_z_xxxxxyz[k];

                g_y_x_zz_xxxxxz[k] = -g_y_x_z_xxxxxz[k] * ab_z + g_y_x_z_xxxxxzz[k];

                g_y_x_zz_xxxxyy[k] = -g_y_x_z_xxxxyy[k] * ab_z + g_y_x_z_xxxxyyz[k];

                g_y_x_zz_xxxxyz[k] = -g_y_x_z_xxxxyz[k] * ab_z + g_y_x_z_xxxxyzz[k];

                g_y_x_zz_xxxxzz[k] = -g_y_x_z_xxxxzz[k] * ab_z + g_y_x_z_xxxxzzz[k];

                g_y_x_zz_xxxyyy[k] = -g_y_x_z_xxxyyy[k] * ab_z + g_y_x_z_xxxyyyz[k];

                g_y_x_zz_xxxyyz[k] = -g_y_x_z_xxxyyz[k] * ab_z + g_y_x_z_xxxyyzz[k];

                g_y_x_zz_xxxyzz[k] = -g_y_x_z_xxxyzz[k] * ab_z + g_y_x_z_xxxyzzz[k];

                g_y_x_zz_xxxzzz[k] = -g_y_x_z_xxxzzz[k] * ab_z + g_y_x_z_xxxzzzz[k];

                g_y_x_zz_xxyyyy[k] = -g_y_x_z_xxyyyy[k] * ab_z + g_y_x_z_xxyyyyz[k];

                g_y_x_zz_xxyyyz[k] = -g_y_x_z_xxyyyz[k] * ab_z + g_y_x_z_xxyyyzz[k];

                g_y_x_zz_xxyyzz[k] = -g_y_x_z_xxyyzz[k] * ab_z + g_y_x_z_xxyyzzz[k];

                g_y_x_zz_xxyzzz[k] = -g_y_x_z_xxyzzz[k] * ab_z + g_y_x_z_xxyzzzz[k];

                g_y_x_zz_xxzzzz[k] = -g_y_x_z_xxzzzz[k] * ab_z + g_y_x_z_xxzzzzz[k];

                g_y_x_zz_xyyyyy[k] = -g_y_x_z_xyyyyy[k] * ab_z + g_y_x_z_xyyyyyz[k];

                g_y_x_zz_xyyyyz[k] = -g_y_x_z_xyyyyz[k] * ab_z + g_y_x_z_xyyyyzz[k];

                g_y_x_zz_xyyyzz[k] = -g_y_x_z_xyyyzz[k] * ab_z + g_y_x_z_xyyyzzz[k];

                g_y_x_zz_xyyzzz[k] = -g_y_x_z_xyyzzz[k] * ab_z + g_y_x_z_xyyzzzz[k];

                g_y_x_zz_xyzzzz[k] = -g_y_x_z_xyzzzz[k] * ab_z + g_y_x_z_xyzzzzz[k];

                g_y_x_zz_xzzzzz[k] = -g_y_x_z_xzzzzz[k] * ab_z + g_y_x_z_xzzzzzz[k];

                g_y_x_zz_yyyyyy[k] = -g_y_x_z_yyyyyy[k] * ab_z + g_y_x_z_yyyyyyz[k];

                g_y_x_zz_yyyyyz[k] = -g_y_x_z_yyyyyz[k] * ab_z + g_y_x_z_yyyyyzz[k];

                g_y_x_zz_yyyyzz[k] = -g_y_x_z_yyyyzz[k] * ab_z + g_y_x_z_yyyyzzz[k];

                g_y_x_zz_yyyzzz[k] = -g_y_x_z_yyyzzz[k] * ab_z + g_y_x_z_yyyzzzz[k];

                g_y_x_zz_yyzzzz[k] = -g_y_x_z_yyzzzz[k] * ab_z + g_y_x_z_yyzzzzz[k];

                g_y_x_zz_yzzzzz[k] = -g_y_x_z_yzzzzz[k] * ab_z + g_y_x_z_yzzzzzz[k];

                g_y_x_zz_zzzzzz[k] = -g_y_x_z_zzzzzz[k] * ab_z + g_y_x_z_zzzzzzz[k];
            }

            /// Set up 672-700 components of targeted buffer : cbuffer.data(

            auto g_y_y_xx_xxxxxx = cbuffer.data(di_geom_11_off + 672 * ccomps * dcomps);

            auto g_y_y_xx_xxxxxy = cbuffer.data(di_geom_11_off + 673 * ccomps * dcomps);

            auto g_y_y_xx_xxxxxz = cbuffer.data(di_geom_11_off + 674 * ccomps * dcomps);

            auto g_y_y_xx_xxxxyy = cbuffer.data(di_geom_11_off + 675 * ccomps * dcomps);

            auto g_y_y_xx_xxxxyz = cbuffer.data(di_geom_11_off + 676 * ccomps * dcomps);

            auto g_y_y_xx_xxxxzz = cbuffer.data(di_geom_11_off + 677 * ccomps * dcomps);

            auto g_y_y_xx_xxxyyy = cbuffer.data(di_geom_11_off + 678 * ccomps * dcomps);

            auto g_y_y_xx_xxxyyz = cbuffer.data(di_geom_11_off + 679 * ccomps * dcomps);

            auto g_y_y_xx_xxxyzz = cbuffer.data(di_geom_11_off + 680 * ccomps * dcomps);

            auto g_y_y_xx_xxxzzz = cbuffer.data(di_geom_11_off + 681 * ccomps * dcomps);

            auto g_y_y_xx_xxyyyy = cbuffer.data(di_geom_11_off + 682 * ccomps * dcomps);

            auto g_y_y_xx_xxyyyz = cbuffer.data(di_geom_11_off + 683 * ccomps * dcomps);

            auto g_y_y_xx_xxyyzz = cbuffer.data(di_geom_11_off + 684 * ccomps * dcomps);

            auto g_y_y_xx_xxyzzz = cbuffer.data(di_geom_11_off + 685 * ccomps * dcomps);

            auto g_y_y_xx_xxzzzz = cbuffer.data(di_geom_11_off + 686 * ccomps * dcomps);

            auto g_y_y_xx_xyyyyy = cbuffer.data(di_geom_11_off + 687 * ccomps * dcomps);

            auto g_y_y_xx_xyyyyz = cbuffer.data(di_geom_11_off + 688 * ccomps * dcomps);

            auto g_y_y_xx_xyyyzz = cbuffer.data(di_geom_11_off + 689 * ccomps * dcomps);

            auto g_y_y_xx_xyyzzz = cbuffer.data(di_geom_11_off + 690 * ccomps * dcomps);

            auto g_y_y_xx_xyzzzz = cbuffer.data(di_geom_11_off + 691 * ccomps * dcomps);

            auto g_y_y_xx_xzzzzz = cbuffer.data(di_geom_11_off + 692 * ccomps * dcomps);

            auto g_y_y_xx_yyyyyy = cbuffer.data(di_geom_11_off + 693 * ccomps * dcomps);

            auto g_y_y_xx_yyyyyz = cbuffer.data(di_geom_11_off + 694 * ccomps * dcomps);

            auto g_y_y_xx_yyyyzz = cbuffer.data(di_geom_11_off + 695 * ccomps * dcomps);

            auto g_y_y_xx_yyyzzz = cbuffer.data(di_geom_11_off + 696 * ccomps * dcomps);

            auto g_y_y_xx_yyzzzz = cbuffer.data(di_geom_11_off + 697 * ccomps * dcomps);

            auto g_y_y_xx_yzzzzz = cbuffer.data(di_geom_11_off + 698 * ccomps * dcomps);

            auto g_y_y_xx_zzzzzz = cbuffer.data(di_geom_11_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_x_xxxxxx, g_y_y_x_xxxxxxx, g_y_y_x_xxxxxxy, g_y_y_x_xxxxxxz, g_y_y_x_xxxxxy, g_y_y_x_xxxxxyy, g_y_y_x_xxxxxyz, g_y_y_x_xxxxxz, g_y_y_x_xxxxxzz, g_y_y_x_xxxxyy, g_y_y_x_xxxxyyy, g_y_y_x_xxxxyyz, g_y_y_x_xxxxyz, g_y_y_x_xxxxyzz, g_y_y_x_xxxxzz, g_y_y_x_xxxxzzz, g_y_y_x_xxxyyy, g_y_y_x_xxxyyyy, g_y_y_x_xxxyyyz, g_y_y_x_xxxyyz, g_y_y_x_xxxyyzz, g_y_y_x_xxxyzz, g_y_y_x_xxxyzzz, g_y_y_x_xxxzzz, g_y_y_x_xxxzzzz, g_y_y_x_xxyyyy, g_y_y_x_xxyyyyy, g_y_y_x_xxyyyyz, g_y_y_x_xxyyyz, g_y_y_x_xxyyyzz, g_y_y_x_xxyyzz, g_y_y_x_xxyyzzz, g_y_y_x_xxyzzz, g_y_y_x_xxyzzzz, g_y_y_x_xxzzzz, g_y_y_x_xxzzzzz, g_y_y_x_xyyyyy, g_y_y_x_xyyyyyy, g_y_y_x_xyyyyyz, g_y_y_x_xyyyyz, g_y_y_x_xyyyyzz, g_y_y_x_xyyyzz, g_y_y_x_xyyyzzz, g_y_y_x_xyyzzz, g_y_y_x_xyyzzzz, g_y_y_x_xyzzzz, g_y_y_x_xyzzzzz, g_y_y_x_xzzzzz, g_y_y_x_xzzzzzz, g_y_y_x_yyyyyy, g_y_y_x_yyyyyz, g_y_y_x_yyyyzz, g_y_y_x_yyyzzz, g_y_y_x_yyzzzz, g_y_y_x_yzzzzz, g_y_y_x_zzzzzz, g_y_y_xx_xxxxxx, g_y_y_xx_xxxxxy, g_y_y_xx_xxxxxz, g_y_y_xx_xxxxyy, g_y_y_xx_xxxxyz, g_y_y_xx_xxxxzz, g_y_y_xx_xxxyyy, g_y_y_xx_xxxyyz, g_y_y_xx_xxxyzz, g_y_y_xx_xxxzzz, g_y_y_xx_xxyyyy, g_y_y_xx_xxyyyz, g_y_y_xx_xxyyzz, g_y_y_xx_xxyzzz, g_y_y_xx_xxzzzz, g_y_y_xx_xyyyyy, g_y_y_xx_xyyyyz, g_y_y_xx_xyyyzz, g_y_y_xx_xyyzzz, g_y_y_xx_xyzzzz, g_y_y_xx_xzzzzz, g_y_y_xx_yyyyyy, g_y_y_xx_yyyyyz, g_y_y_xx_yyyyzz, g_y_y_xx_yyyzzz, g_y_y_xx_yyzzzz, g_y_y_xx_yzzzzz, g_y_y_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xx_xxxxxx[k] = -g_y_y_x_xxxxxx[k] * ab_x + g_y_y_x_xxxxxxx[k];

                g_y_y_xx_xxxxxy[k] = -g_y_y_x_xxxxxy[k] * ab_x + g_y_y_x_xxxxxxy[k];

                g_y_y_xx_xxxxxz[k] = -g_y_y_x_xxxxxz[k] * ab_x + g_y_y_x_xxxxxxz[k];

                g_y_y_xx_xxxxyy[k] = -g_y_y_x_xxxxyy[k] * ab_x + g_y_y_x_xxxxxyy[k];

                g_y_y_xx_xxxxyz[k] = -g_y_y_x_xxxxyz[k] * ab_x + g_y_y_x_xxxxxyz[k];

                g_y_y_xx_xxxxzz[k] = -g_y_y_x_xxxxzz[k] * ab_x + g_y_y_x_xxxxxzz[k];

                g_y_y_xx_xxxyyy[k] = -g_y_y_x_xxxyyy[k] * ab_x + g_y_y_x_xxxxyyy[k];

                g_y_y_xx_xxxyyz[k] = -g_y_y_x_xxxyyz[k] * ab_x + g_y_y_x_xxxxyyz[k];

                g_y_y_xx_xxxyzz[k] = -g_y_y_x_xxxyzz[k] * ab_x + g_y_y_x_xxxxyzz[k];

                g_y_y_xx_xxxzzz[k] = -g_y_y_x_xxxzzz[k] * ab_x + g_y_y_x_xxxxzzz[k];

                g_y_y_xx_xxyyyy[k] = -g_y_y_x_xxyyyy[k] * ab_x + g_y_y_x_xxxyyyy[k];

                g_y_y_xx_xxyyyz[k] = -g_y_y_x_xxyyyz[k] * ab_x + g_y_y_x_xxxyyyz[k];

                g_y_y_xx_xxyyzz[k] = -g_y_y_x_xxyyzz[k] * ab_x + g_y_y_x_xxxyyzz[k];

                g_y_y_xx_xxyzzz[k] = -g_y_y_x_xxyzzz[k] * ab_x + g_y_y_x_xxxyzzz[k];

                g_y_y_xx_xxzzzz[k] = -g_y_y_x_xxzzzz[k] * ab_x + g_y_y_x_xxxzzzz[k];

                g_y_y_xx_xyyyyy[k] = -g_y_y_x_xyyyyy[k] * ab_x + g_y_y_x_xxyyyyy[k];

                g_y_y_xx_xyyyyz[k] = -g_y_y_x_xyyyyz[k] * ab_x + g_y_y_x_xxyyyyz[k];

                g_y_y_xx_xyyyzz[k] = -g_y_y_x_xyyyzz[k] * ab_x + g_y_y_x_xxyyyzz[k];

                g_y_y_xx_xyyzzz[k] = -g_y_y_x_xyyzzz[k] * ab_x + g_y_y_x_xxyyzzz[k];

                g_y_y_xx_xyzzzz[k] = -g_y_y_x_xyzzzz[k] * ab_x + g_y_y_x_xxyzzzz[k];

                g_y_y_xx_xzzzzz[k] = -g_y_y_x_xzzzzz[k] * ab_x + g_y_y_x_xxzzzzz[k];

                g_y_y_xx_yyyyyy[k] = -g_y_y_x_yyyyyy[k] * ab_x + g_y_y_x_xyyyyyy[k];

                g_y_y_xx_yyyyyz[k] = -g_y_y_x_yyyyyz[k] * ab_x + g_y_y_x_xyyyyyz[k];

                g_y_y_xx_yyyyzz[k] = -g_y_y_x_yyyyzz[k] * ab_x + g_y_y_x_xyyyyzz[k];

                g_y_y_xx_yyyzzz[k] = -g_y_y_x_yyyzzz[k] * ab_x + g_y_y_x_xyyyzzz[k];

                g_y_y_xx_yyzzzz[k] = -g_y_y_x_yyzzzz[k] * ab_x + g_y_y_x_xyyzzzz[k];

                g_y_y_xx_yzzzzz[k] = -g_y_y_x_yzzzzz[k] * ab_x + g_y_y_x_xyzzzzz[k];

                g_y_y_xx_zzzzzz[k] = -g_y_y_x_zzzzzz[k] * ab_x + g_y_y_x_xzzzzzz[k];
            }

            /// Set up 700-728 components of targeted buffer : cbuffer.data(

            auto g_y_y_xy_xxxxxx = cbuffer.data(di_geom_11_off + 700 * ccomps * dcomps);

            auto g_y_y_xy_xxxxxy = cbuffer.data(di_geom_11_off + 701 * ccomps * dcomps);

            auto g_y_y_xy_xxxxxz = cbuffer.data(di_geom_11_off + 702 * ccomps * dcomps);

            auto g_y_y_xy_xxxxyy = cbuffer.data(di_geom_11_off + 703 * ccomps * dcomps);

            auto g_y_y_xy_xxxxyz = cbuffer.data(di_geom_11_off + 704 * ccomps * dcomps);

            auto g_y_y_xy_xxxxzz = cbuffer.data(di_geom_11_off + 705 * ccomps * dcomps);

            auto g_y_y_xy_xxxyyy = cbuffer.data(di_geom_11_off + 706 * ccomps * dcomps);

            auto g_y_y_xy_xxxyyz = cbuffer.data(di_geom_11_off + 707 * ccomps * dcomps);

            auto g_y_y_xy_xxxyzz = cbuffer.data(di_geom_11_off + 708 * ccomps * dcomps);

            auto g_y_y_xy_xxxzzz = cbuffer.data(di_geom_11_off + 709 * ccomps * dcomps);

            auto g_y_y_xy_xxyyyy = cbuffer.data(di_geom_11_off + 710 * ccomps * dcomps);

            auto g_y_y_xy_xxyyyz = cbuffer.data(di_geom_11_off + 711 * ccomps * dcomps);

            auto g_y_y_xy_xxyyzz = cbuffer.data(di_geom_11_off + 712 * ccomps * dcomps);

            auto g_y_y_xy_xxyzzz = cbuffer.data(di_geom_11_off + 713 * ccomps * dcomps);

            auto g_y_y_xy_xxzzzz = cbuffer.data(di_geom_11_off + 714 * ccomps * dcomps);

            auto g_y_y_xy_xyyyyy = cbuffer.data(di_geom_11_off + 715 * ccomps * dcomps);

            auto g_y_y_xy_xyyyyz = cbuffer.data(di_geom_11_off + 716 * ccomps * dcomps);

            auto g_y_y_xy_xyyyzz = cbuffer.data(di_geom_11_off + 717 * ccomps * dcomps);

            auto g_y_y_xy_xyyzzz = cbuffer.data(di_geom_11_off + 718 * ccomps * dcomps);

            auto g_y_y_xy_xyzzzz = cbuffer.data(di_geom_11_off + 719 * ccomps * dcomps);

            auto g_y_y_xy_xzzzzz = cbuffer.data(di_geom_11_off + 720 * ccomps * dcomps);

            auto g_y_y_xy_yyyyyy = cbuffer.data(di_geom_11_off + 721 * ccomps * dcomps);

            auto g_y_y_xy_yyyyyz = cbuffer.data(di_geom_11_off + 722 * ccomps * dcomps);

            auto g_y_y_xy_yyyyzz = cbuffer.data(di_geom_11_off + 723 * ccomps * dcomps);

            auto g_y_y_xy_yyyzzz = cbuffer.data(di_geom_11_off + 724 * ccomps * dcomps);

            auto g_y_y_xy_yyzzzz = cbuffer.data(di_geom_11_off + 725 * ccomps * dcomps);

            auto g_y_y_xy_yzzzzz = cbuffer.data(di_geom_11_off + 726 * ccomps * dcomps);

            auto g_y_y_xy_zzzzzz = cbuffer.data(di_geom_11_off + 727 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xy_xxxxxx, g_y_y_xy_xxxxxy, g_y_y_xy_xxxxxz, g_y_y_xy_xxxxyy, g_y_y_xy_xxxxyz, g_y_y_xy_xxxxzz, g_y_y_xy_xxxyyy, g_y_y_xy_xxxyyz, g_y_y_xy_xxxyzz, g_y_y_xy_xxxzzz, g_y_y_xy_xxyyyy, g_y_y_xy_xxyyyz, g_y_y_xy_xxyyzz, g_y_y_xy_xxyzzz, g_y_y_xy_xxzzzz, g_y_y_xy_xyyyyy, g_y_y_xy_xyyyyz, g_y_y_xy_xyyyzz, g_y_y_xy_xyyzzz, g_y_y_xy_xyzzzz, g_y_y_xy_xzzzzz, g_y_y_xy_yyyyyy, g_y_y_xy_yyyyyz, g_y_y_xy_yyyyzz, g_y_y_xy_yyyzzz, g_y_y_xy_yyzzzz, g_y_y_xy_yzzzzz, g_y_y_xy_zzzzzz, g_y_y_y_xxxxxx, g_y_y_y_xxxxxxx, g_y_y_y_xxxxxxy, g_y_y_y_xxxxxxz, g_y_y_y_xxxxxy, g_y_y_y_xxxxxyy, g_y_y_y_xxxxxyz, g_y_y_y_xxxxxz, g_y_y_y_xxxxxzz, g_y_y_y_xxxxyy, g_y_y_y_xxxxyyy, g_y_y_y_xxxxyyz, g_y_y_y_xxxxyz, g_y_y_y_xxxxyzz, g_y_y_y_xxxxzz, g_y_y_y_xxxxzzz, g_y_y_y_xxxyyy, g_y_y_y_xxxyyyy, g_y_y_y_xxxyyyz, g_y_y_y_xxxyyz, g_y_y_y_xxxyyzz, g_y_y_y_xxxyzz, g_y_y_y_xxxyzzz, g_y_y_y_xxxzzz, g_y_y_y_xxxzzzz, g_y_y_y_xxyyyy, g_y_y_y_xxyyyyy, g_y_y_y_xxyyyyz, g_y_y_y_xxyyyz, g_y_y_y_xxyyyzz, g_y_y_y_xxyyzz, g_y_y_y_xxyyzzz, g_y_y_y_xxyzzz, g_y_y_y_xxyzzzz, g_y_y_y_xxzzzz, g_y_y_y_xxzzzzz, g_y_y_y_xyyyyy, g_y_y_y_xyyyyyy, g_y_y_y_xyyyyyz, g_y_y_y_xyyyyz, g_y_y_y_xyyyyzz, g_y_y_y_xyyyzz, g_y_y_y_xyyyzzz, g_y_y_y_xyyzzz, g_y_y_y_xyyzzzz, g_y_y_y_xyzzzz, g_y_y_y_xyzzzzz, g_y_y_y_xzzzzz, g_y_y_y_xzzzzzz, g_y_y_y_yyyyyy, g_y_y_y_yyyyyz, g_y_y_y_yyyyzz, g_y_y_y_yyyzzz, g_y_y_y_yyzzzz, g_y_y_y_yzzzzz, g_y_y_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xy_xxxxxx[k] = -g_y_y_y_xxxxxx[k] * ab_x + g_y_y_y_xxxxxxx[k];

                g_y_y_xy_xxxxxy[k] = -g_y_y_y_xxxxxy[k] * ab_x + g_y_y_y_xxxxxxy[k];

                g_y_y_xy_xxxxxz[k] = -g_y_y_y_xxxxxz[k] * ab_x + g_y_y_y_xxxxxxz[k];

                g_y_y_xy_xxxxyy[k] = -g_y_y_y_xxxxyy[k] * ab_x + g_y_y_y_xxxxxyy[k];

                g_y_y_xy_xxxxyz[k] = -g_y_y_y_xxxxyz[k] * ab_x + g_y_y_y_xxxxxyz[k];

                g_y_y_xy_xxxxzz[k] = -g_y_y_y_xxxxzz[k] * ab_x + g_y_y_y_xxxxxzz[k];

                g_y_y_xy_xxxyyy[k] = -g_y_y_y_xxxyyy[k] * ab_x + g_y_y_y_xxxxyyy[k];

                g_y_y_xy_xxxyyz[k] = -g_y_y_y_xxxyyz[k] * ab_x + g_y_y_y_xxxxyyz[k];

                g_y_y_xy_xxxyzz[k] = -g_y_y_y_xxxyzz[k] * ab_x + g_y_y_y_xxxxyzz[k];

                g_y_y_xy_xxxzzz[k] = -g_y_y_y_xxxzzz[k] * ab_x + g_y_y_y_xxxxzzz[k];

                g_y_y_xy_xxyyyy[k] = -g_y_y_y_xxyyyy[k] * ab_x + g_y_y_y_xxxyyyy[k];

                g_y_y_xy_xxyyyz[k] = -g_y_y_y_xxyyyz[k] * ab_x + g_y_y_y_xxxyyyz[k];

                g_y_y_xy_xxyyzz[k] = -g_y_y_y_xxyyzz[k] * ab_x + g_y_y_y_xxxyyzz[k];

                g_y_y_xy_xxyzzz[k] = -g_y_y_y_xxyzzz[k] * ab_x + g_y_y_y_xxxyzzz[k];

                g_y_y_xy_xxzzzz[k] = -g_y_y_y_xxzzzz[k] * ab_x + g_y_y_y_xxxzzzz[k];

                g_y_y_xy_xyyyyy[k] = -g_y_y_y_xyyyyy[k] * ab_x + g_y_y_y_xxyyyyy[k];

                g_y_y_xy_xyyyyz[k] = -g_y_y_y_xyyyyz[k] * ab_x + g_y_y_y_xxyyyyz[k];

                g_y_y_xy_xyyyzz[k] = -g_y_y_y_xyyyzz[k] * ab_x + g_y_y_y_xxyyyzz[k];

                g_y_y_xy_xyyzzz[k] = -g_y_y_y_xyyzzz[k] * ab_x + g_y_y_y_xxyyzzz[k];

                g_y_y_xy_xyzzzz[k] = -g_y_y_y_xyzzzz[k] * ab_x + g_y_y_y_xxyzzzz[k];

                g_y_y_xy_xzzzzz[k] = -g_y_y_y_xzzzzz[k] * ab_x + g_y_y_y_xxzzzzz[k];

                g_y_y_xy_yyyyyy[k] = -g_y_y_y_yyyyyy[k] * ab_x + g_y_y_y_xyyyyyy[k];

                g_y_y_xy_yyyyyz[k] = -g_y_y_y_yyyyyz[k] * ab_x + g_y_y_y_xyyyyyz[k];

                g_y_y_xy_yyyyzz[k] = -g_y_y_y_yyyyzz[k] * ab_x + g_y_y_y_xyyyyzz[k];

                g_y_y_xy_yyyzzz[k] = -g_y_y_y_yyyzzz[k] * ab_x + g_y_y_y_xyyyzzz[k];

                g_y_y_xy_yyzzzz[k] = -g_y_y_y_yyzzzz[k] * ab_x + g_y_y_y_xyyzzzz[k];

                g_y_y_xy_yzzzzz[k] = -g_y_y_y_yzzzzz[k] * ab_x + g_y_y_y_xyzzzzz[k];

                g_y_y_xy_zzzzzz[k] = -g_y_y_y_zzzzzz[k] * ab_x + g_y_y_y_xzzzzzz[k];
            }

            /// Set up 728-756 components of targeted buffer : cbuffer.data(

            auto g_y_y_xz_xxxxxx = cbuffer.data(di_geom_11_off + 728 * ccomps * dcomps);

            auto g_y_y_xz_xxxxxy = cbuffer.data(di_geom_11_off + 729 * ccomps * dcomps);

            auto g_y_y_xz_xxxxxz = cbuffer.data(di_geom_11_off + 730 * ccomps * dcomps);

            auto g_y_y_xz_xxxxyy = cbuffer.data(di_geom_11_off + 731 * ccomps * dcomps);

            auto g_y_y_xz_xxxxyz = cbuffer.data(di_geom_11_off + 732 * ccomps * dcomps);

            auto g_y_y_xz_xxxxzz = cbuffer.data(di_geom_11_off + 733 * ccomps * dcomps);

            auto g_y_y_xz_xxxyyy = cbuffer.data(di_geom_11_off + 734 * ccomps * dcomps);

            auto g_y_y_xz_xxxyyz = cbuffer.data(di_geom_11_off + 735 * ccomps * dcomps);

            auto g_y_y_xz_xxxyzz = cbuffer.data(di_geom_11_off + 736 * ccomps * dcomps);

            auto g_y_y_xz_xxxzzz = cbuffer.data(di_geom_11_off + 737 * ccomps * dcomps);

            auto g_y_y_xz_xxyyyy = cbuffer.data(di_geom_11_off + 738 * ccomps * dcomps);

            auto g_y_y_xz_xxyyyz = cbuffer.data(di_geom_11_off + 739 * ccomps * dcomps);

            auto g_y_y_xz_xxyyzz = cbuffer.data(di_geom_11_off + 740 * ccomps * dcomps);

            auto g_y_y_xz_xxyzzz = cbuffer.data(di_geom_11_off + 741 * ccomps * dcomps);

            auto g_y_y_xz_xxzzzz = cbuffer.data(di_geom_11_off + 742 * ccomps * dcomps);

            auto g_y_y_xz_xyyyyy = cbuffer.data(di_geom_11_off + 743 * ccomps * dcomps);

            auto g_y_y_xz_xyyyyz = cbuffer.data(di_geom_11_off + 744 * ccomps * dcomps);

            auto g_y_y_xz_xyyyzz = cbuffer.data(di_geom_11_off + 745 * ccomps * dcomps);

            auto g_y_y_xz_xyyzzz = cbuffer.data(di_geom_11_off + 746 * ccomps * dcomps);

            auto g_y_y_xz_xyzzzz = cbuffer.data(di_geom_11_off + 747 * ccomps * dcomps);

            auto g_y_y_xz_xzzzzz = cbuffer.data(di_geom_11_off + 748 * ccomps * dcomps);

            auto g_y_y_xz_yyyyyy = cbuffer.data(di_geom_11_off + 749 * ccomps * dcomps);

            auto g_y_y_xz_yyyyyz = cbuffer.data(di_geom_11_off + 750 * ccomps * dcomps);

            auto g_y_y_xz_yyyyzz = cbuffer.data(di_geom_11_off + 751 * ccomps * dcomps);

            auto g_y_y_xz_yyyzzz = cbuffer.data(di_geom_11_off + 752 * ccomps * dcomps);

            auto g_y_y_xz_yyzzzz = cbuffer.data(di_geom_11_off + 753 * ccomps * dcomps);

            auto g_y_y_xz_yzzzzz = cbuffer.data(di_geom_11_off + 754 * ccomps * dcomps);

            auto g_y_y_xz_zzzzzz = cbuffer.data(di_geom_11_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xz_xxxxxx, g_y_y_xz_xxxxxy, g_y_y_xz_xxxxxz, g_y_y_xz_xxxxyy, g_y_y_xz_xxxxyz, g_y_y_xz_xxxxzz, g_y_y_xz_xxxyyy, g_y_y_xz_xxxyyz, g_y_y_xz_xxxyzz, g_y_y_xz_xxxzzz, g_y_y_xz_xxyyyy, g_y_y_xz_xxyyyz, g_y_y_xz_xxyyzz, g_y_y_xz_xxyzzz, g_y_y_xz_xxzzzz, g_y_y_xz_xyyyyy, g_y_y_xz_xyyyyz, g_y_y_xz_xyyyzz, g_y_y_xz_xyyzzz, g_y_y_xz_xyzzzz, g_y_y_xz_xzzzzz, g_y_y_xz_yyyyyy, g_y_y_xz_yyyyyz, g_y_y_xz_yyyyzz, g_y_y_xz_yyyzzz, g_y_y_xz_yyzzzz, g_y_y_xz_yzzzzz, g_y_y_xz_zzzzzz, g_y_y_z_xxxxxx, g_y_y_z_xxxxxxx, g_y_y_z_xxxxxxy, g_y_y_z_xxxxxxz, g_y_y_z_xxxxxy, g_y_y_z_xxxxxyy, g_y_y_z_xxxxxyz, g_y_y_z_xxxxxz, g_y_y_z_xxxxxzz, g_y_y_z_xxxxyy, g_y_y_z_xxxxyyy, g_y_y_z_xxxxyyz, g_y_y_z_xxxxyz, g_y_y_z_xxxxyzz, g_y_y_z_xxxxzz, g_y_y_z_xxxxzzz, g_y_y_z_xxxyyy, g_y_y_z_xxxyyyy, g_y_y_z_xxxyyyz, g_y_y_z_xxxyyz, g_y_y_z_xxxyyzz, g_y_y_z_xxxyzz, g_y_y_z_xxxyzzz, g_y_y_z_xxxzzz, g_y_y_z_xxxzzzz, g_y_y_z_xxyyyy, g_y_y_z_xxyyyyy, g_y_y_z_xxyyyyz, g_y_y_z_xxyyyz, g_y_y_z_xxyyyzz, g_y_y_z_xxyyzz, g_y_y_z_xxyyzzz, g_y_y_z_xxyzzz, g_y_y_z_xxyzzzz, g_y_y_z_xxzzzz, g_y_y_z_xxzzzzz, g_y_y_z_xyyyyy, g_y_y_z_xyyyyyy, g_y_y_z_xyyyyyz, g_y_y_z_xyyyyz, g_y_y_z_xyyyyzz, g_y_y_z_xyyyzz, g_y_y_z_xyyyzzz, g_y_y_z_xyyzzz, g_y_y_z_xyyzzzz, g_y_y_z_xyzzzz, g_y_y_z_xyzzzzz, g_y_y_z_xzzzzz, g_y_y_z_xzzzzzz, g_y_y_z_yyyyyy, g_y_y_z_yyyyyz, g_y_y_z_yyyyzz, g_y_y_z_yyyzzz, g_y_y_z_yyzzzz, g_y_y_z_yzzzzz, g_y_y_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xz_xxxxxx[k] = -g_y_y_z_xxxxxx[k] * ab_x + g_y_y_z_xxxxxxx[k];

                g_y_y_xz_xxxxxy[k] = -g_y_y_z_xxxxxy[k] * ab_x + g_y_y_z_xxxxxxy[k];

                g_y_y_xz_xxxxxz[k] = -g_y_y_z_xxxxxz[k] * ab_x + g_y_y_z_xxxxxxz[k];

                g_y_y_xz_xxxxyy[k] = -g_y_y_z_xxxxyy[k] * ab_x + g_y_y_z_xxxxxyy[k];

                g_y_y_xz_xxxxyz[k] = -g_y_y_z_xxxxyz[k] * ab_x + g_y_y_z_xxxxxyz[k];

                g_y_y_xz_xxxxzz[k] = -g_y_y_z_xxxxzz[k] * ab_x + g_y_y_z_xxxxxzz[k];

                g_y_y_xz_xxxyyy[k] = -g_y_y_z_xxxyyy[k] * ab_x + g_y_y_z_xxxxyyy[k];

                g_y_y_xz_xxxyyz[k] = -g_y_y_z_xxxyyz[k] * ab_x + g_y_y_z_xxxxyyz[k];

                g_y_y_xz_xxxyzz[k] = -g_y_y_z_xxxyzz[k] * ab_x + g_y_y_z_xxxxyzz[k];

                g_y_y_xz_xxxzzz[k] = -g_y_y_z_xxxzzz[k] * ab_x + g_y_y_z_xxxxzzz[k];

                g_y_y_xz_xxyyyy[k] = -g_y_y_z_xxyyyy[k] * ab_x + g_y_y_z_xxxyyyy[k];

                g_y_y_xz_xxyyyz[k] = -g_y_y_z_xxyyyz[k] * ab_x + g_y_y_z_xxxyyyz[k];

                g_y_y_xz_xxyyzz[k] = -g_y_y_z_xxyyzz[k] * ab_x + g_y_y_z_xxxyyzz[k];

                g_y_y_xz_xxyzzz[k] = -g_y_y_z_xxyzzz[k] * ab_x + g_y_y_z_xxxyzzz[k];

                g_y_y_xz_xxzzzz[k] = -g_y_y_z_xxzzzz[k] * ab_x + g_y_y_z_xxxzzzz[k];

                g_y_y_xz_xyyyyy[k] = -g_y_y_z_xyyyyy[k] * ab_x + g_y_y_z_xxyyyyy[k];

                g_y_y_xz_xyyyyz[k] = -g_y_y_z_xyyyyz[k] * ab_x + g_y_y_z_xxyyyyz[k];

                g_y_y_xz_xyyyzz[k] = -g_y_y_z_xyyyzz[k] * ab_x + g_y_y_z_xxyyyzz[k];

                g_y_y_xz_xyyzzz[k] = -g_y_y_z_xyyzzz[k] * ab_x + g_y_y_z_xxyyzzz[k];

                g_y_y_xz_xyzzzz[k] = -g_y_y_z_xyzzzz[k] * ab_x + g_y_y_z_xxyzzzz[k];

                g_y_y_xz_xzzzzz[k] = -g_y_y_z_xzzzzz[k] * ab_x + g_y_y_z_xxzzzzz[k];

                g_y_y_xz_yyyyyy[k] = -g_y_y_z_yyyyyy[k] * ab_x + g_y_y_z_xyyyyyy[k];

                g_y_y_xz_yyyyyz[k] = -g_y_y_z_yyyyyz[k] * ab_x + g_y_y_z_xyyyyyz[k];

                g_y_y_xz_yyyyzz[k] = -g_y_y_z_yyyyzz[k] * ab_x + g_y_y_z_xyyyyzz[k];

                g_y_y_xz_yyyzzz[k] = -g_y_y_z_yyyzzz[k] * ab_x + g_y_y_z_xyyyzzz[k];

                g_y_y_xz_yyzzzz[k] = -g_y_y_z_yyzzzz[k] * ab_x + g_y_y_z_xyyzzzz[k];

                g_y_y_xz_yzzzzz[k] = -g_y_y_z_yzzzzz[k] * ab_x + g_y_y_z_xyzzzzz[k];

                g_y_y_xz_zzzzzz[k] = -g_y_y_z_zzzzzz[k] * ab_x + g_y_y_z_xzzzzzz[k];
            }

            /// Set up 756-784 components of targeted buffer : cbuffer.data(

            auto g_y_y_yy_xxxxxx = cbuffer.data(di_geom_11_off + 756 * ccomps * dcomps);

            auto g_y_y_yy_xxxxxy = cbuffer.data(di_geom_11_off + 757 * ccomps * dcomps);

            auto g_y_y_yy_xxxxxz = cbuffer.data(di_geom_11_off + 758 * ccomps * dcomps);

            auto g_y_y_yy_xxxxyy = cbuffer.data(di_geom_11_off + 759 * ccomps * dcomps);

            auto g_y_y_yy_xxxxyz = cbuffer.data(di_geom_11_off + 760 * ccomps * dcomps);

            auto g_y_y_yy_xxxxzz = cbuffer.data(di_geom_11_off + 761 * ccomps * dcomps);

            auto g_y_y_yy_xxxyyy = cbuffer.data(di_geom_11_off + 762 * ccomps * dcomps);

            auto g_y_y_yy_xxxyyz = cbuffer.data(di_geom_11_off + 763 * ccomps * dcomps);

            auto g_y_y_yy_xxxyzz = cbuffer.data(di_geom_11_off + 764 * ccomps * dcomps);

            auto g_y_y_yy_xxxzzz = cbuffer.data(di_geom_11_off + 765 * ccomps * dcomps);

            auto g_y_y_yy_xxyyyy = cbuffer.data(di_geom_11_off + 766 * ccomps * dcomps);

            auto g_y_y_yy_xxyyyz = cbuffer.data(di_geom_11_off + 767 * ccomps * dcomps);

            auto g_y_y_yy_xxyyzz = cbuffer.data(di_geom_11_off + 768 * ccomps * dcomps);

            auto g_y_y_yy_xxyzzz = cbuffer.data(di_geom_11_off + 769 * ccomps * dcomps);

            auto g_y_y_yy_xxzzzz = cbuffer.data(di_geom_11_off + 770 * ccomps * dcomps);

            auto g_y_y_yy_xyyyyy = cbuffer.data(di_geom_11_off + 771 * ccomps * dcomps);

            auto g_y_y_yy_xyyyyz = cbuffer.data(di_geom_11_off + 772 * ccomps * dcomps);

            auto g_y_y_yy_xyyyzz = cbuffer.data(di_geom_11_off + 773 * ccomps * dcomps);

            auto g_y_y_yy_xyyzzz = cbuffer.data(di_geom_11_off + 774 * ccomps * dcomps);

            auto g_y_y_yy_xyzzzz = cbuffer.data(di_geom_11_off + 775 * ccomps * dcomps);

            auto g_y_y_yy_xzzzzz = cbuffer.data(di_geom_11_off + 776 * ccomps * dcomps);

            auto g_y_y_yy_yyyyyy = cbuffer.data(di_geom_11_off + 777 * ccomps * dcomps);

            auto g_y_y_yy_yyyyyz = cbuffer.data(di_geom_11_off + 778 * ccomps * dcomps);

            auto g_y_y_yy_yyyyzz = cbuffer.data(di_geom_11_off + 779 * ccomps * dcomps);

            auto g_y_y_yy_yyyzzz = cbuffer.data(di_geom_11_off + 780 * ccomps * dcomps);

            auto g_y_y_yy_yyzzzz = cbuffer.data(di_geom_11_off + 781 * ccomps * dcomps);

            auto g_y_y_yy_yzzzzz = cbuffer.data(di_geom_11_off + 782 * ccomps * dcomps);

            auto g_y_y_yy_zzzzzz = cbuffer.data(di_geom_11_off + 783 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxxxxx, g_0_y_y_xxxxxy, g_0_y_y_xxxxxz, g_0_y_y_xxxxyy, g_0_y_y_xxxxyz, g_0_y_y_xxxxzz, g_0_y_y_xxxyyy, g_0_y_y_xxxyyz, g_0_y_y_xxxyzz, g_0_y_y_xxxzzz, g_0_y_y_xxyyyy, g_0_y_y_xxyyyz, g_0_y_y_xxyyzz, g_0_y_y_xxyzzz, g_0_y_y_xxzzzz, g_0_y_y_xyyyyy, g_0_y_y_xyyyyz, g_0_y_y_xyyyzz, g_0_y_y_xyyzzz, g_0_y_y_xyzzzz, g_0_y_y_xzzzzz, g_0_y_y_yyyyyy, g_0_y_y_yyyyyz, g_0_y_y_yyyyzz, g_0_y_y_yyyzzz, g_0_y_y_yyzzzz, g_0_y_y_yzzzzz, g_0_y_y_zzzzzz, g_y_0_y_xxxxxx, g_y_0_y_xxxxxy, g_y_0_y_xxxxxz, g_y_0_y_xxxxyy, g_y_0_y_xxxxyz, g_y_0_y_xxxxzz, g_y_0_y_xxxyyy, g_y_0_y_xxxyyz, g_y_0_y_xxxyzz, g_y_0_y_xxxzzz, g_y_0_y_xxyyyy, g_y_0_y_xxyyyz, g_y_0_y_xxyyzz, g_y_0_y_xxyzzz, g_y_0_y_xxzzzz, g_y_0_y_xyyyyy, g_y_0_y_xyyyyz, g_y_0_y_xyyyzz, g_y_0_y_xyyzzz, g_y_0_y_xyzzzz, g_y_0_y_xzzzzz, g_y_0_y_yyyyyy, g_y_0_y_yyyyyz, g_y_0_y_yyyyzz, g_y_0_y_yyyzzz, g_y_0_y_yyzzzz, g_y_0_y_yzzzzz, g_y_0_y_zzzzzz, g_y_y_y_xxxxxx, g_y_y_y_xxxxxxy, g_y_y_y_xxxxxy, g_y_y_y_xxxxxyy, g_y_y_y_xxxxxyz, g_y_y_y_xxxxxz, g_y_y_y_xxxxyy, g_y_y_y_xxxxyyy, g_y_y_y_xxxxyyz, g_y_y_y_xxxxyz, g_y_y_y_xxxxyzz, g_y_y_y_xxxxzz, g_y_y_y_xxxyyy, g_y_y_y_xxxyyyy, g_y_y_y_xxxyyyz, g_y_y_y_xxxyyz, g_y_y_y_xxxyyzz, g_y_y_y_xxxyzz, g_y_y_y_xxxyzzz, g_y_y_y_xxxzzz, g_y_y_y_xxyyyy, g_y_y_y_xxyyyyy, g_y_y_y_xxyyyyz, g_y_y_y_xxyyyz, g_y_y_y_xxyyyzz, g_y_y_y_xxyyzz, g_y_y_y_xxyyzzz, g_y_y_y_xxyzzz, g_y_y_y_xxyzzzz, g_y_y_y_xxzzzz, g_y_y_y_xyyyyy, g_y_y_y_xyyyyyy, g_y_y_y_xyyyyyz, g_y_y_y_xyyyyz, g_y_y_y_xyyyyzz, g_y_y_y_xyyyzz, g_y_y_y_xyyyzzz, g_y_y_y_xyyzzz, g_y_y_y_xyyzzzz, g_y_y_y_xyzzzz, g_y_y_y_xyzzzzz, g_y_y_y_xzzzzz, g_y_y_y_yyyyyy, g_y_y_y_yyyyyyy, g_y_y_y_yyyyyyz, g_y_y_y_yyyyyz, g_y_y_y_yyyyyzz, g_y_y_y_yyyyzz, g_y_y_y_yyyyzzz, g_y_y_y_yyyzzz, g_y_y_y_yyyzzzz, g_y_y_y_yyzzzz, g_y_y_y_yyzzzzz, g_y_y_y_yzzzzz, g_y_y_y_yzzzzzz, g_y_y_y_zzzzzz, g_y_y_yy_xxxxxx, g_y_y_yy_xxxxxy, g_y_y_yy_xxxxxz, g_y_y_yy_xxxxyy, g_y_y_yy_xxxxyz, g_y_y_yy_xxxxzz, g_y_y_yy_xxxyyy, g_y_y_yy_xxxyyz, g_y_y_yy_xxxyzz, g_y_y_yy_xxxzzz, g_y_y_yy_xxyyyy, g_y_y_yy_xxyyyz, g_y_y_yy_xxyyzz, g_y_y_yy_xxyzzz, g_y_y_yy_xxzzzz, g_y_y_yy_xyyyyy, g_y_y_yy_xyyyyz, g_y_y_yy_xyyyzz, g_y_y_yy_xyyzzz, g_y_y_yy_xyzzzz, g_y_y_yy_xzzzzz, g_y_y_yy_yyyyyy, g_y_y_yy_yyyyyz, g_y_y_yy_yyyyzz, g_y_y_yy_yyyzzz, g_y_y_yy_yyzzzz, g_y_y_yy_yzzzzz, g_y_y_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yy_xxxxxx[k] = -g_0_y_y_xxxxxx[k] + g_y_0_y_xxxxxx[k] - g_y_y_y_xxxxxx[k] * ab_y + g_y_y_y_xxxxxxy[k];

                g_y_y_yy_xxxxxy[k] = -g_0_y_y_xxxxxy[k] + g_y_0_y_xxxxxy[k] - g_y_y_y_xxxxxy[k] * ab_y + g_y_y_y_xxxxxyy[k];

                g_y_y_yy_xxxxxz[k] = -g_0_y_y_xxxxxz[k] + g_y_0_y_xxxxxz[k] - g_y_y_y_xxxxxz[k] * ab_y + g_y_y_y_xxxxxyz[k];

                g_y_y_yy_xxxxyy[k] = -g_0_y_y_xxxxyy[k] + g_y_0_y_xxxxyy[k] - g_y_y_y_xxxxyy[k] * ab_y + g_y_y_y_xxxxyyy[k];

                g_y_y_yy_xxxxyz[k] = -g_0_y_y_xxxxyz[k] + g_y_0_y_xxxxyz[k] - g_y_y_y_xxxxyz[k] * ab_y + g_y_y_y_xxxxyyz[k];

                g_y_y_yy_xxxxzz[k] = -g_0_y_y_xxxxzz[k] + g_y_0_y_xxxxzz[k] - g_y_y_y_xxxxzz[k] * ab_y + g_y_y_y_xxxxyzz[k];

                g_y_y_yy_xxxyyy[k] = -g_0_y_y_xxxyyy[k] + g_y_0_y_xxxyyy[k] - g_y_y_y_xxxyyy[k] * ab_y + g_y_y_y_xxxyyyy[k];

                g_y_y_yy_xxxyyz[k] = -g_0_y_y_xxxyyz[k] + g_y_0_y_xxxyyz[k] - g_y_y_y_xxxyyz[k] * ab_y + g_y_y_y_xxxyyyz[k];

                g_y_y_yy_xxxyzz[k] = -g_0_y_y_xxxyzz[k] + g_y_0_y_xxxyzz[k] - g_y_y_y_xxxyzz[k] * ab_y + g_y_y_y_xxxyyzz[k];

                g_y_y_yy_xxxzzz[k] = -g_0_y_y_xxxzzz[k] + g_y_0_y_xxxzzz[k] - g_y_y_y_xxxzzz[k] * ab_y + g_y_y_y_xxxyzzz[k];

                g_y_y_yy_xxyyyy[k] = -g_0_y_y_xxyyyy[k] + g_y_0_y_xxyyyy[k] - g_y_y_y_xxyyyy[k] * ab_y + g_y_y_y_xxyyyyy[k];

                g_y_y_yy_xxyyyz[k] = -g_0_y_y_xxyyyz[k] + g_y_0_y_xxyyyz[k] - g_y_y_y_xxyyyz[k] * ab_y + g_y_y_y_xxyyyyz[k];

                g_y_y_yy_xxyyzz[k] = -g_0_y_y_xxyyzz[k] + g_y_0_y_xxyyzz[k] - g_y_y_y_xxyyzz[k] * ab_y + g_y_y_y_xxyyyzz[k];

                g_y_y_yy_xxyzzz[k] = -g_0_y_y_xxyzzz[k] + g_y_0_y_xxyzzz[k] - g_y_y_y_xxyzzz[k] * ab_y + g_y_y_y_xxyyzzz[k];

                g_y_y_yy_xxzzzz[k] = -g_0_y_y_xxzzzz[k] + g_y_0_y_xxzzzz[k] - g_y_y_y_xxzzzz[k] * ab_y + g_y_y_y_xxyzzzz[k];

                g_y_y_yy_xyyyyy[k] = -g_0_y_y_xyyyyy[k] + g_y_0_y_xyyyyy[k] - g_y_y_y_xyyyyy[k] * ab_y + g_y_y_y_xyyyyyy[k];

                g_y_y_yy_xyyyyz[k] = -g_0_y_y_xyyyyz[k] + g_y_0_y_xyyyyz[k] - g_y_y_y_xyyyyz[k] * ab_y + g_y_y_y_xyyyyyz[k];

                g_y_y_yy_xyyyzz[k] = -g_0_y_y_xyyyzz[k] + g_y_0_y_xyyyzz[k] - g_y_y_y_xyyyzz[k] * ab_y + g_y_y_y_xyyyyzz[k];

                g_y_y_yy_xyyzzz[k] = -g_0_y_y_xyyzzz[k] + g_y_0_y_xyyzzz[k] - g_y_y_y_xyyzzz[k] * ab_y + g_y_y_y_xyyyzzz[k];

                g_y_y_yy_xyzzzz[k] = -g_0_y_y_xyzzzz[k] + g_y_0_y_xyzzzz[k] - g_y_y_y_xyzzzz[k] * ab_y + g_y_y_y_xyyzzzz[k];

                g_y_y_yy_xzzzzz[k] = -g_0_y_y_xzzzzz[k] + g_y_0_y_xzzzzz[k] - g_y_y_y_xzzzzz[k] * ab_y + g_y_y_y_xyzzzzz[k];

                g_y_y_yy_yyyyyy[k] = -g_0_y_y_yyyyyy[k] + g_y_0_y_yyyyyy[k] - g_y_y_y_yyyyyy[k] * ab_y + g_y_y_y_yyyyyyy[k];

                g_y_y_yy_yyyyyz[k] = -g_0_y_y_yyyyyz[k] + g_y_0_y_yyyyyz[k] - g_y_y_y_yyyyyz[k] * ab_y + g_y_y_y_yyyyyyz[k];

                g_y_y_yy_yyyyzz[k] = -g_0_y_y_yyyyzz[k] + g_y_0_y_yyyyzz[k] - g_y_y_y_yyyyzz[k] * ab_y + g_y_y_y_yyyyyzz[k];

                g_y_y_yy_yyyzzz[k] = -g_0_y_y_yyyzzz[k] + g_y_0_y_yyyzzz[k] - g_y_y_y_yyyzzz[k] * ab_y + g_y_y_y_yyyyzzz[k];

                g_y_y_yy_yyzzzz[k] = -g_0_y_y_yyzzzz[k] + g_y_0_y_yyzzzz[k] - g_y_y_y_yyzzzz[k] * ab_y + g_y_y_y_yyyzzzz[k];

                g_y_y_yy_yzzzzz[k] = -g_0_y_y_yzzzzz[k] + g_y_0_y_yzzzzz[k] - g_y_y_y_yzzzzz[k] * ab_y + g_y_y_y_yyzzzzz[k];

                g_y_y_yy_zzzzzz[k] = -g_0_y_y_zzzzzz[k] + g_y_0_y_zzzzzz[k] - g_y_y_y_zzzzzz[k] * ab_y + g_y_y_y_yzzzzzz[k];
            }

            /// Set up 784-812 components of targeted buffer : cbuffer.data(

            auto g_y_y_yz_xxxxxx = cbuffer.data(di_geom_11_off + 784 * ccomps * dcomps);

            auto g_y_y_yz_xxxxxy = cbuffer.data(di_geom_11_off + 785 * ccomps * dcomps);

            auto g_y_y_yz_xxxxxz = cbuffer.data(di_geom_11_off + 786 * ccomps * dcomps);

            auto g_y_y_yz_xxxxyy = cbuffer.data(di_geom_11_off + 787 * ccomps * dcomps);

            auto g_y_y_yz_xxxxyz = cbuffer.data(di_geom_11_off + 788 * ccomps * dcomps);

            auto g_y_y_yz_xxxxzz = cbuffer.data(di_geom_11_off + 789 * ccomps * dcomps);

            auto g_y_y_yz_xxxyyy = cbuffer.data(di_geom_11_off + 790 * ccomps * dcomps);

            auto g_y_y_yz_xxxyyz = cbuffer.data(di_geom_11_off + 791 * ccomps * dcomps);

            auto g_y_y_yz_xxxyzz = cbuffer.data(di_geom_11_off + 792 * ccomps * dcomps);

            auto g_y_y_yz_xxxzzz = cbuffer.data(di_geom_11_off + 793 * ccomps * dcomps);

            auto g_y_y_yz_xxyyyy = cbuffer.data(di_geom_11_off + 794 * ccomps * dcomps);

            auto g_y_y_yz_xxyyyz = cbuffer.data(di_geom_11_off + 795 * ccomps * dcomps);

            auto g_y_y_yz_xxyyzz = cbuffer.data(di_geom_11_off + 796 * ccomps * dcomps);

            auto g_y_y_yz_xxyzzz = cbuffer.data(di_geom_11_off + 797 * ccomps * dcomps);

            auto g_y_y_yz_xxzzzz = cbuffer.data(di_geom_11_off + 798 * ccomps * dcomps);

            auto g_y_y_yz_xyyyyy = cbuffer.data(di_geom_11_off + 799 * ccomps * dcomps);

            auto g_y_y_yz_xyyyyz = cbuffer.data(di_geom_11_off + 800 * ccomps * dcomps);

            auto g_y_y_yz_xyyyzz = cbuffer.data(di_geom_11_off + 801 * ccomps * dcomps);

            auto g_y_y_yz_xyyzzz = cbuffer.data(di_geom_11_off + 802 * ccomps * dcomps);

            auto g_y_y_yz_xyzzzz = cbuffer.data(di_geom_11_off + 803 * ccomps * dcomps);

            auto g_y_y_yz_xzzzzz = cbuffer.data(di_geom_11_off + 804 * ccomps * dcomps);

            auto g_y_y_yz_yyyyyy = cbuffer.data(di_geom_11_off + 805 * ccomps * dcomps);

            auto g_y_y_yz_yyyyyz = cbuffer.data(di_geom_11_off + 806 * ccomps * dcomps);

            auto g_y_y_yz_yyyyzz = cbuffer.data(di_geom_11_off + 807 * ccomps * dcomps);

            auto g_y_y_yz_yyyzzz = cbuffer.data(di_geom_11_off + 808 * ccomps * dcomps);

            auto g_y_y_yz_yyzzzz = cbuffer.data(di_geom_11_off + 809 * ccomps * dcomps);

            auto g_y_y_yz_yzzzzz = cbuffer.data(di_geom_11_off + 810 * ccomps * dcomps);

            auto g_y_y_yz_zzzzzz = cbuffer.data(di_geom_11_off + 811 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_y_xxxxxx, g_y_y_y_xxxxxxz, g_y_y_y_xxxxxy, g_y_y_y_xxxxxyz, g_y_y_y_xxxxxz, g_y_y_y_xxxxxzz, g_y_y_y_xxxxyy, g_y_y_y_xxxxyyz, g_y_y_y_xxxxyz, g_y_y_y_xxxxyzz, g_y_y_y_xxxxzz, g_y_y_y_xxxxzzz, g_y_y_y_xxxyyy, g_y_y_y_xxxyyyz, g_y_y_y_xxxyyz, g_y_y_y_xxxyyzz, g_y_y_y_xxxyzz, g_y_y_y_xxxyzzz, g_y_y_y_xxxzzz, g_y_y_y_xxxzzzz, g_y_y_y_xxyyyy, g_y_y_y_xxyyyyz, g_y_y_y_xxyyyz, g_y_y_y_xxyyyzz, g_y_y_y_xxyyzz, g_y_y_y_xxyyzzz, g_y_y_y_xxyzzz, g_y_y_y_xxyzzzz, g_y_y_y_xxzzzz, g_y_y_y_xxzzzzz, g_y_y_y_xyyyyy, g_y_y_y_xyyyyyz, g_y_y_y_xyyyyz, g_y_y_y_xyyyyzz, g_y_y_y_xyyyzz, g_y_y_y_xyyyzzz, g_y_y_y_xyyzzz, g_y_y_y_xyyzzzz, g_y_y_y_xyzzzz, g_y_y_y_xyzzzzz, g_y_y_y_xzzzzz, g_y_y_y_xzzzzzz, g_y_y_y_yyyyyy, g_y_y_y_yyyyyyz, g_y_y_y_yyyyyz, g_y_y_y_yyyyyzz, g_y_y_y_yyyyzz, g_y_y_y_yyyyzzz, g_y_y_y_yyyzzz, g_y_y_y_yyyzzzz, g_y_y_y_yyzzzz, g_y_y_y_yyzzzzz, g_y_y_y_yzzzzz, g_y_y_y_yzzzzzz, g_y_y_y_zzzzzz, g_y_y_y_zzzzzzz, g_y_y_yz_xxxxxx, g_y_y_yz_xxxxxy, g_y_y_yz_xxxxxz, g_y_y_yz_xxxxyy, g_y_y_yz_xxxxyz, g_y_y_yz_xxxxzz, g_y_y_yz_xxxyyy, g_y_y_yz_xxxyyz, g_y_y_yz_xxxyzz, g_y_y_yz_xxxzzz, g_y_y_yz_xxyyyy, g_y_y_yz_xxyyyz, g_y_y_yz_xxyyzz, g_y_y_yz_xxyzzz, g_y_y_yz_xxzzzz, g_y_y_yz_xyyyyy, g_y_y_yz_xyyyyz, g_y_y_yz_xyyyzz, g_y_y_yz_xyyzzz, g_y_y_yz_xyzzzz, g_y_y_yz_xzzzzz, g_y_y_yz_yyyyyy, g_y_y_yz_yyyyyz, g_y_y_yz_yyyyzz, g_y_y_yz_yyyzzz, g_y_y_yz_yyzzzz, g_y_y_yz_yzzzzz, g_y_y_yz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yz_xxxxxx[k] = -g_y_y_y_xxxxxx[k] * ab_z + g_y_y_y_xxxxxxz[k];

                g_y_y_yz_xxxxxy[k] = -g_y_y_y_xxxxxy[k] * ab_z + g_y_y_y_xxxxxyz[k];

                g_y_y_yz_xxxxxz[k] = -g_y_y_y_xxxxxz[k] * ab_z + g_y_y_y_xxxxxzz[k];

                g_y_y_yz_xxxxyy[k] = -g_y_y_y_xxxxyy[k] * ab_z + g_y_y_y_xxxxyyz[k];

                g_y_y_yz_xxxxyz[k] = -g_y_y_y_xxxxyz[k] * ab_z + g_y_y_y_xxxxyzz[k];

                g_y_y_yz_xxxxzz[k] = -g_y_y_y_xxxxzz[k] * ab_z + g_y_y_y_xxxxzzz[k];

                g_y_y_yz_xxxyyy[k] = -g_y_y_y_xxxyyy[k] * ab_z + g_y_y_y_xxxyyyz[k];

                g_y_y_yz_xxxyyz[k] = -g_y_y_y_xxxyyz[k] * ab_z + g_y_y_y_xxxyyzz[k];

                g_y_y_yz_xxxyzz[k] = -g_y_y_y_xxxyzz[k] * ab_z + g_y_y_y_xxxyzzz[k];

                g_y_y_yz_xxxzzz[k] = -g_y_y_y_xxxzzz[k] * ab_z + g_y_y_y_xxxzzzz[k];

                g_y_y_yz_xxyyyy[k] = -g_y_y_y_xxyyyy[k] * ab_z + g_y_y_y_xxyyyyz[k];

                g_y_y_yz_xxyyyz[k] = -g_y_y_y_xxyyyz[k] * ab_z + g_y_y_y_xxyyyzz[k];

                g_y_y_yz_xxyyzz[k] = -g_y_y_y_xxyyzz[k] * ab_z + g_y_y_y_xxyyzzz[k];

                g_y_y_yz_xxyzzz[k] = -g_y_y_y_xxyzzz[k] * ab_z + g_y_y_y_xxyzzzz[k];

                g_y_y_yz_xxzzzz[k] = -g_y_y_y_xxzzzz[k] * ab_z + g_y_y_y_xxzzzzz[k];

                g_y_y_yz_xyyyyy[k] = -g_y_y_y_xyyyyy[k] * ab_z + g_y_y_y_xyyyyyz[k];

                g_y_y_yz_xyyyyz[k] = -g_y_y_y_xyyyyz[k] * ab_z + g_y_y_y_xyyyyzz[k];

                g_y_y_yz_xyyyzz[k] = -g_y_y_y_xyyyzz[k] * ab_z + g_y_y_y_xyyyzzz[k];

                g_y_y_yz_xyyzzz[k] = -g_y_y_y_xyyzzz[k] * ab_z + g_y_y_y_xyyzzzz[k];

                g_y_y_yz_xyzzzz[k] = -g_y_y_y_xyzzzz[k] * ab_z + g_y_y_y_xyzzzzz[k];

                g_y_y_yz_xzzzzz[k] = -g_y_y_y_xzzzzz[k] * ab_z + g_y_y_y_xzzzzzz[k];

                g_y_y_yz_yyyyyy[k] = -g_y_y_y_yyyyyy[k] * ab_z + g_y_y_y_yyyyyyz[k];

                g_y_y_yz_yyyyyz[k] = -g_y_y_y_yyyyyz[k] * ab_z + g_y_y_y_yyyyyzz[k];

                g_y_y_yz_yyyyzz[k] = -g_y_y_y_yyyyzz[k] * ab_z + g_y_y_y_yyyyzzz[k];

                g_y_y_yz_yyyzzz[k] = -g_y_y_y_yyyzzz[k] * ab_z + g_y_y_y_yyyzzzz[k];

                g_y_y_yz_yyzzzz[k] = -g_y_y_y_yyzzzz[k] * ab_z + g_y_y_y_yyzzzzz[k];

                g_y_y_yz_yzzzzz[k] = -g_y_y_y_yzzzzz[k] * ab_z + g_y_y_y_yzzzzzz[k];

                g_y_y_yz_zzzzzz[k] = -g_y_y_y_zzzzzz[k] * ab_z + g_y_y_y_zzzzzzz[k];
            }

            /// Set up 812-840 components of targeted buffer : cbuffer.data(

            auto g_y_y_zz_xxxxxx = cbuffer.data(di_geom_11_off + 812 * ccomps * dcomps);

            auto g_y_y_zz_xxxxxy = cbuffer.data(di_geom_11_off + 813 * ccomps * dcomps);

            auto g_y_y_zz_xxxxxz = cbuffer.data(di_geom_11_off + 814 * ccomps * dcomps);

            auto g_y_y_zz_xxxxyy = cbuffer.data(di_geom_11_off + 815 * ccomps * dcomps);

            auto g_y_y_zz_xxxxyz = cbuffer.data(di_geom_11_off + 816 * ccomps * dcomps);

            auto g_y_y_zz_xxxxzz = cbuffer.data(di_geom_11_off + 817 * ccomps * dcomps);

            auto g_y_y_zz_xxxyyy = cbuffer.data(di_geom_11_off + 818 * ccomps * dcomps);

            auto g_y_y_zz_xxxyyz = cbuffer.data(di_geom_11_off + 819 * ccomps * dcomps);

            auto g_y_y_zz_xxxyzz = cbuffer.data(di_geom_11_off + 820 * ccomps * dcomps);

            auto g_y_y_zz_xxxzzz = cbuffer.data(di_geom_11_off + 821 * ccomps * dcomps);

            auto g_y_y_zz_xxyyyy = cbuffer.data(di_geom_11_off + 822 * ccomps * dcomps);

            auto g_y_y_zz_xxyyyz = cbuffer.data(di_geom_11_off + 823 * ccomps * dcomps);

            auto g_y_y_zz_xxyyzz = cbuffer.data(di_geom_11_off + 824 * ccomps * dcomps);

            auto g_y_y_zz_xxyzzz = cbuffer.data(di_geom_11_off + 825 * ccomps * dcomps);

            auto g_y_y_zz_xxzzzz = cbuffer.data(di_geom_11_off + 826 * ccomps * dcomps);

            auto g_y_y_zz_xyyyyy = cbuffer.data(di_geom_11_off + 827 * ccomps * dcomps);

            auto g_y_y_zz_xyyyyz = cbuffer.data(di_geom_11_off + 828 * ccomps * dcomps);

            auto g_y_y_zz_xyyyzz = cbuffer.data(di_geom_11_off + 829 * ccomps * dcomps);

            auto g_y_y_zz_xyyzzz = cbuffer.data(di_geom_11_off + 830 * ccomps * dcomps);

            auto g_y_y_zz_xyzzzz = cbuffer.data(di_geom_11_off + 831 * ccomps * dcomps);

            auto g_y_y_zz_xzzzzz = cbuffer.data(di_geom_11_off + 832 * ccomps * dcomps);

            auto g_y_y_zz_yyyyyy = cbuffer.data(di_geom_11_off + 833 * ccomps * dcomps);

            auto g_y_y_zz_yyyyyz = cbuffer.data(di_geom_11_off + 834 * ccomps * dcomps);

            auto g_y_y_zz_yyyyzz = cbuffer.data(di_geom_11_off + 835 * ccomps * dcomps);

            auto g_y_y_zz_yyyzzz = cbuffer.data(di_geom_11_off + 836 * ccomps * dcomps);

            auto g_y_y_zz_yyzzzz = cbuffer.data(di_geom_11_off + 837 * ccomps * dcomps);

            auto g_y_y_zz_yzzzzz = cbuffer.data(di_geom_11_off + 838 * ccomps * dcomps);

            auto g_y_y_zz_zzzzzz = cbuffer.data(di_geom_11_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_z_xxxxxx, g_y_y_z_xxxxxxz, g_y_y_z_xxxxxy, g_y_y_z_xxxxxyz, g_y_y_z_xxxxxz, g_y_y_z_xxxxxzz, g_y_y_z_xxxxyy, g_y_y_z_xxxxyyz, g_y_y_z_xxxxyz, g_y_y_z_xxxxyzz, g_y_y_z_xxxxzz, g_y_y_z_xxxxzzz, g_y_y_z_xxxyyy, g_y_y_z_xxxyyyz, g_y_y_z_xxxyyz, g_y_y_z_xxxyyzz, g_y_y_z_xxxyzz, g_y_y_z_xxxyzzz, g_y_y_z_xxxzzz, g_y_y_z_xxxzzzz, g_y_y_z_xxyyyy, g_y_y_z_xxyyyyz, g_y_y_z_xxyyyz, g_y_y_z_xxyyyzz, g_y_y_z_xxyyzz, g_y_y_z_xxyyzzz, g_y_y_z_xxyzzz, g_y_y_z_xxyzzzz, g_y_y_z_xxzzzz, g_y_y_z_xxzzzzz, g_y_y_z_xyyyyy, g_y_y_z_xyyyyyz, g_y_y_z_xyyyyz, g_y_y_z_xyyyyzz, g_y_y_z_xyyyzz, g_y_y_z_xyyyzzz, g_y_y_z_xyyzzz, g_y_y_z_xyyzzzz, g_y_y_z_xyzzzz, g_y_y_z_xyzzzzz, g_y_y_z_xzzzzz, g_y_y_z_xzzzzzz, g_y_y_z_yyyyyy, g_y_y_z_yyyyyyz, g_y_y_z_yyyyyz, g_y_y_z_yyyyyzz, g_y_y_z_yyyyzz, g_y_y_z_yyyyzzz, g_y_y_z_yyyzzz, g_y_y_z_yyyzzzz, g_y_y_z_yyzzzz, g_y_y_z_yyzzzzz, g_y_y_z_yzzzzz, g_y_y_z_yzzzzzz, g_y_y_z_zzzzzz, g_y_y_z_zzzzzzz, g_y_y_zz_xxxxxx, g_y_y_zz_xxxxxy, g_y_y_zz_xxxxxz, g_y_y_zz_xxxxyy, g_y_y_zz_xxxxyz, g_y_y_zz_xxxxzz, g_y_y_zz_xxxyyy, g_y_y_zz_xxxyyz, g_y_y_zz_xxxyzz, g_y_y_zz_xxxzzz, g_y_y_zz_xxyyyy, g_y_y_zz_xxyyyz, g_y_y_zz_xxyyzz, g_y_y_zz_xxyzzz, g_y_y_zz_xxzzzz, g_y_y_zz_xyyyyy, g_y_y_zz_xyyyyz, g_y_y_zz_xyyyzz, g_y_y_zz_xyyzzz, g_y_y_zz_xyzzzz, g_y_y_zz_xzzzzz, g_y_y_zz_yyyyyy, g_y_y_zz_yyyyyz, g_y_y_zz_yyyyzz, g_y_y_zz_yyyzzz, g_y_y_zz_yyzzzz, g_y_y_zz_yzzzzz, g_y_y_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zz_xxxxxx[k] = -g_y_y_z_xxxxxx[k] * ab_z + g_y_y_z_xxxxxxz[k];

                g_y_y_zz_xxxxxy[k] = -g_y_y_z_xxxxxy[k] * ab_z + g_y_y_z_xxxxxyz[k];

                g_y_y_zz_xxxxxz[k] = -g_y_y_z_xxxxxz[k] * ab_z + g_y_y_z_xxxxxzz[k];

                g_y_y_zz_xxxxyy[k] = -g_y_y_z_xxxxyy[k] * ab_z + g_y_y_z_xxxxyyz[k];

                g_y_y_zz_xxxxyz[k] = -g_y_y_z_xxxxyz[k] * ab_z + g_y_y_z_xxxxyzz[k];

                g_y_y_zz_xxxxzz[k] = -g_y_y_z_xxxxzz[k] * ab_z + g_y_y_z_xxxxzzz[k];

                g_y_y_zz_xxxyyy[k] = -g_y_y_z_xxxyyy[k] * ab_z + g_y_y_z_xxxyyyz[k];

                g_y_y_zz_xxxyyz[k] = -g_y_y_z_xxxyyz[k] * ab_z + g_y_y_z_xxxyyzz[k];

                g_y_y_zz_xxxyzz[k] = -g_y_y_z_xxxyzz[k] * ab_z + g_y_y_z_xxxyzzz[k];

                g_y_y_zz_xxxzzz[k] = -g_y_y_z_xxxzzz[k] * ab_z + g_y_y_z_xxxzzzz[k];

                g_y_y_zz_xxyyyy[k] = -g_y_y_z_xxyyyy[k] * ab_z + g_y_y_z_xxyyyyz[k];

                g_y_y_zz_xxyyyz[k] = -g_y_y_z_xxyyyz[k] * ab_z + g_y_y_z_xxyyyzz[k];

                g_y_y_zz_xxyyzz[k] = -g_y_y_z_xxyyzz[k] * ab_z + g_y_y_z_xxyyzzz[k];

                g_y_y_zz_xxyzzz[k] = -g_y_y_z_xxyzzz[k] * ab_z + g_y_y_z_xxyzzzz[k];

                g_y_y_zz_xxzzzz[k] = -g_y_y_z_xxzzzz[k] * ab_z + g_y_y_z_xxzzzzz[k];

                g_y_y_zz_xyyyyy[k] = -g_y_y_z_xyyyyy[k] * ab_z + g_y_y_z_xyyyyyz[k];

                g_y_y_zz_xyyyyz[k] = -g_y_y_z_xyyyyz[k] * ab_z + g_y_y_z_xyyyyzz[k];

                g_y_y_zz_xyyyzz[k] = -g_y_y_z_xyyyzz[k] * ab_z + g_y_y_z_xyyyzzz[k];

                g_y_y_zz_xyyzzz[k] = -g_y_y_z_xyyzzz[k] * ab_z + g_y_y_z_xyyzzzz[k];

                g_y_y_zz_xyzzzz[k] = -g_y_y_z_xyzzzz[k] * ab_z + g_y_y_z_xyzzzzz[k];

                g_y_y_zz_xzzzzz[k] = -g_y_y_z_xzzzzz[k] * ab_z + g_y_y_z_xzzzzzz[k];

                g_y_y_zz_yyyyyy[k] = -g_y_y_z_yyyyyy[k] * ab_z + g_y_y_z_yyyyyyz[k];

                g_y_y_zz_yyyyyz[k] = -g_y_y_z_yyyyyz[k] * ab_z + g_y_y_z_yyyyyzz[k];

                g_y_y_zz_yyyyzz[k] = -g_y_y_z_yyyyzz[k] * ab_z + g_y_y_z_yyyyzzz[k];

                g_y_y_zz_yyyzzz[k] = -g_y_y_z_yyyzzz[k] * ab_z + g_y_y_z_yyyzzzz[k];

                g_y_y_zz_yyzzzz[k] = -g_y_y_z_yyzzzz[k] * ab_z + g_y_y_z_yyzzzzz[k];

                g_y_y_zz_yzzzzz[k] = -g_y_y_z_yzzzzz[k] * ab_z + g_y_y_z_yzzzzzz[k];

                g_y_y_zz_zzzzzz[k] = -g_y_y_z_zzzzzz[k] * ab_z + g_y_y_z_zzzzzzz[k];
            }

            /// Set up 840-868 components of targeted buffer : cbuffer.data(

            auto g_y_z_xx_xxxxxx = cbuffer.data(di_geom_11_off + 840 * ccomps * dcomps);

            auto g_y_z_xx_xxxxxy = cbuffer.data(di_geom_11_off + 841 * ccomps * dcomps);

            auto g_y_z_xx_xxxxxz = cbuffer.data(di_geom_11_off + 842 * ccomps * dcomps);

            auto g_y_z_xx_xxxxyy = cbuffer.data(di_geom_11_off + 843 * ccomps * dcomps);

            auto g_y_z_xx_xxxxyz = cbuffer.data(di_geom_11_off + 844 * ccomps * dcomps);

            auto g_y_z_xx_xxxxzz = cbuffer.data(di_geom_11_off + 845 * ccomps * dcomps);

            auto g_y_z_xx_xxxyyy = cbuffer.data(di_geom_11_off + 846 * ccomps * dcomps);

            auto g_y_z_xx_xxxyyz = cbuffer.data(di_geom_11_off + 847 * ccomps * dcomps);

            auto g_y_z_xx_xxxyzz = cbuffer.data(di_geom_11_off + 848 * ccomps * dcomps);

            auto g_y_z_xx_xxxzzz = cbuffer.data(di_geom_11_off + 849 * ccomps * dcomps);

            auto g_y_z_xx_xxyyyy = cbuffer.data(di_geom_11_off + 850 * ccomps * dcomps);

            auto g_y_z_xx_xxyyyz = cbuffer.data(di_geom_11_off + 851 * ccomps * dcomps);

            auto g_y_z_xx_xxyyzz = cbuffer.data(di_geom_11_off + 852 * ccomps * dcomps);

            auto g_y_z_xx_xxyzzz = cbuffer.data(di_geom_11_off + 853 * ccomps * dcomps);

            auto g_y_z_xx_xxzzzz = cbuffer.data(di_geom_11_off + 854 * ccomps * dcomps);

            auto g_y_z_xx_xyyyyy = cbuffer.data(di_geom_11_off + 855 * ccomps * dcomps);

            auto g_y_z_xx_xyyyyz = cbuffer.data(di_geom_11_off + 856 * ccomps * dcomps);

            auto g_y_z_xx_xyyyzz = cbuffer.data(di_geom_11_off + 857 * ccomps * dcomps);

            auto g_y_z_xx_xyyzzz = cbuffer.data(di_geom_11_off + 858 * ccomps * dcomps);

            auto g_y_z_xx_xyzzzz = cbuffer.data(di_geom_11_off + 859 * ccomps * dcomps);

            auto g_y_z_xx_xzzzzz = cbuffer.data(di_geom_11_off + 860 * ccomps * dcomps);

            auto g_y_z_xx_yyyyyy = cbuffer.data(di_geom_11_off + 861 * ccomps * dcomps);

            auto g_y_z_xx_yyyyyz = cbuffer.data(di_geom_11_off + 862 * ccomps * dcomps);

            auto g_y_z_xx_yyyyzz = cbuffer.data(di_geom_11_off + 863 * ccomps * dcomps);

            auto g_y_z_xx_yyyzzz = cbuffer.data(di_geom_11_off + 864 * ccomps * dcomps);

            auto g_y_z_xx_yyzzzz = cbuffer.data(di_geom_11_off + 865 * ccomps * dcomps);

            auto g_y_z_xx_yzzzzz = cbuffer.data(di_geom_11_off + 866 * ccomps * dcomps);

            auto g_y_z_xx_zzzzzz = cbuffer.data(di_geom_11_off + 867 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_x_xxxxxx, g_y_z_x_xxxxxxx, g_y_z_x_xxxxxxy, g_y_z_x_xxxxxxz, g_y_z_x_xxxxxy, g_y_z_x_xxxxxyy, g_y_z_x_xxxxxyz, g_y_z_x_xxxxxz, g_y_z_x_xxxxxzz, g_y_z_x_xxxxyy, g_y_z_x_xxxxyyy, g_y_z_x_xxxxyyz, g_y_z_x_xxxxyz, g_y_z_x_xxxxyzz, g_y_z_x_xxxxzz, g_y_z_x_xxxxzzz, g_y_z_x_xxxyyy, g_y_z_x_xxxyyyy, g_y_z_x_xxxyyyz, g_y_z_x_xxxyyz, g_y_z_x_xxxyyzz, g_y_z_x_xxxyzz, g_y_z_x_xxxyzzz, g_y_z_x_xxxzzz, g_y_z_x_xxxzzzz, g_y_z_x_xxyyyy, g_y_z_x_xxyyyyy, g_y_z_x_xxyyyyz, g_y_z_x_xxyyyz, g_y_z_x_xxyyyzz, g_y_z_x_xxyyzz, g_y_z_x_xxyyzzz, g_y_z_x_xxyzzz, g_y_z_x_xxyzzzz, g_y_z_x_xxzzzz, g_y_z_x_xxzzzzz, g_y_z_x_xyyyyy, g_y_z_x_xyyyyyy, g_y_z_x_xyyyyyz, g_y_z_x_xyyyyz, g_y_z_x_xyyyyzz, g_y_z_x_xyyyzz, g_y_z_x_xyyyzzz, g_y_z_x_xyyzzz, g_y_z_x_xyyzzzz, g_y_z_x_xyzzzz, g_y_z_x_xyzzzzz, g_y_z_x_xzzzzz, g_y_z_x_xzzzzzz, g_y_z_x_yyyyyy, g_y_z_x_yyyyyz, g_y_z_x_yyyyzz, g_y_z_x_yyyzzz, g_y_z_x_yyzzzz, g_y_z_x_yzzzzz, g_y_z_x_zzzzzz, g_y_z_xx_xxxxxx, g_y_z_xx_xxxxxy, g_y_z_xx_xxxxxz, g_y_z_xx_xxxxyy, g_y_z_xx_xxxxyz, g_y_z_xx_xxxxzz, g_y_z_xx_xxxyyy, g_y_z_xx_xxxyyz, g_y_z_xx_xxxyzz, g_y_z_xx_xxxzzz, g_y_z_xx_xxyyyy, g_y_z_xx_xxyyyz, g_y_z_xx_xxyyzz, g_y_z_xx_xxyzzz, g_y_z_xx_xxzzzz, g_y_z_xx_xyyyyy, g_y_z_xx_xyyyyz, g_y_z_xx_xyyyzz, g_y_z_xx_xyyzzz, g_y_z_xx_xyzzzz, g_y_z_xx_xzzzzz, g_y_z_xx_yyyyyy, g_y_z_xx_yyyyyz, g_y_z_xx_yyyyzz, g_y_z_xx_yyyzzz, g_y_z_xx_yyzzzz, g_y_z_xx_yzzzzz, g_y_z_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xx_xxxxxx[k] = -g_y_z_x_xxxxxx[k] * ab_x + g_y_z_x_xxxxxxx[k];

                g_y_z_xx_xxxxxy[k] = -g_y_z_x_xxxxxy[k] * ab_x + g_y_z_x_xxxxxxy[k];

                g_y_z_xx_xxxxxz[k] = -g_y_z_x_xxxxxz[k] * ab_x + g_y_z_x_xxxxxxz[k];

                g_y_z_xx_xxxxyy[k] = -g_y_z_x_xxxxyy[k] * ab_x + g_y_z_x_xxxxxyy[k];

                g_y_z_xx_xxxxyz[k] = -g_y_z_x_xxxxyz[k] * ab_x + g_y_z_x_xxxxxyz[k];

                g_y_z_xx_xxxxzz[k] = -g_y_z_x_xxxxzz[k] * ab_x + g_y_z_x_xxxxxzz[k];

                g_y_z_xx_xxxyyy[k] = -g_y_z_x_xxxyyy[k] * ab_x + g_y_z_x_xxxxyyy[k];

                g_y_z_xx_xxxyyz[k] = -g_y_z_x_xxxyyz[k] * ab_x + g_y_z_x_xxxxyyz[k];

                g_y_z_xx_xxxyzz[k] = -g_y_z_x_xxxyzz[k] * ab_x + g_y_z_x_xxxxyzz[k];

                g_y_z_xx_xxxzzz[k] = -g_y_z_x_xxxzzz[k] * ab_x + g_y_z_x_xxxxzzz[k];

                g_y_z_xx_xxyyyy[k] = -g_y_z_x_xxyyyy[k] * ab_x + g_y_z_x_xxxyyyy[k];

                g_y_z_xx_xxyyyz[k] = -g_y_z_x_xxyyyz[k] * ab_x + g_y_z_x_xxxyyyz[k];

                g_y_z_xx_xxyyzz[k] = -g_y_z_x_xxyyzz[k] * ab_x + g_y_z_x_xxxyyzz[k];

                g_y_z_xx_xxyzzz[k] = -g_y_z_x_xxyzzz[k] * ab_x + g_y_z_x_xxxyzzz[k];

                g_y_z_xx_xxzzzz[k] = -g_y_z_x_xxzzzz[k] * ab_x + g_y_z_x_xxxzzzz[k];

                g_y_z_xx_xyyyyy[k] = -g_y_z_x_xyyyyy[k] * ab_x + g_y_z_x_xxyyyyy[k];

                g_y_z_xx_xyyyyz[k] = -g_y_z_x_xyyyyz[k] * ab_x + g_y_z_x_xxyyyyz[k];

                g_y_z_xx_xyyyzz[k] = -g_y_z_x_xyyyzz[k] * ab_x + g_y_z_x_xxyyyzz[k];

                g_y_z_xx_xyyzzz[k] = -g_y_z_x_xyyzzz[k] * ab_x + g_y_z_x_xxyyzzz[k];

                g_y_z_xx_xyzzzz[k] = -g_y_z_x_xyzzzz[k] * ab_x + g_y_z_x_xxyzzzz[k];

                g_y_z_xx_xzzzzz[k] = -g_y_z_x_xzzzzz[k] * ab_x + g_y_z_x_xxzzzzz[k];

                g_y_z_xx_yyyyyy[k] = -g_y_z_x_yyyyyy[k] * ab_x + g_y_z_x_xyyyyyy[k];

                g_y_z_xx_yyyyyz[k] = -g_y_z_x_yyyyyz[k] * ab_x + g_y_z_x_xyyyyyz[k];

                g_y_z_xx_yyyyzz[k] = -g_y_z_x_yyyyzz[k] * ab_x + g_y_z_x_xyyyyzz[k];

                g_y_z_xx_yyyzzz[k] = -g_y_z_x_yyyzzz[k] * ab_x + g_y_z_x_xyyyzzz[k];

                g_y_z_xx_yyzzzz[k] = -g_y_z_x_yyzzzz[k] * ab_x + g_y_z_x_xyyzzzz[k];

                g_y_z_xx_yzzzzz[k] = -g_y_z_x_yzzzzz[k] * ab_x + g_y_z_x_xyzzzzz[k];

                g_y_z_xx_zzzzzz[k] = -g_y_z_x_zzzzzz[k] * ab_x + g_y_z_x_xzzzzzz[k];
            }

            /// Set up 868-896 components of targeted buffer : cbuffer.data(

            auto g_y_z_xy_xxxxxx = cbuffer.data(di_geom_11_off + 868 * ccomps * dcomps);

            auto g_y_z_xy_xxxxxy = cbuffer.data(di_geom_11_off + 869 * ccomps * dcomps);

            auto g_y_z_xy_xxxxxz = cbuffer.data(di_geom_11_off + 870 * ccomps * dcomps);

            auto g_y_z_xy_xxxxyy = cbuffer.data(di_geom_11_off + 871 * ccomps * dcomps);

            auto g_y_z_xy_xxxxyz = cbuffer.data(di_geom_11_off + 872 * ccomps * dcomps);

            auto g_y_z_xy_xxxxzz = cbuffer.data(di_geom_11_off + 873 * ccomps * dcomps);

            auto g_y_z_xy_xxxyyy = cbuffer.data(di_geom_11_off + 874 * ccomps * dcomps);

            auto g_y_z_xy_xxxyyz = cbuffer.data(di_geom_11_off + 875 * ccomps * dcomps);

            auto g_y_z_xy_xxxyzz = cbuffer.data(di_geom_11_off + 876 * ccomps * dcomps);

            auto g_y_z_xy_xxxzzz = cbuffer.data(di_geom_11_off + 877 * ccomps * dcomps);

            auto g_y_z_xy_xxyyyy = cbuffer.data(di_geom_11_off + 878 * ccomps * dcomps);

            auto g_y_z_xy_xxyyyz = cbuffer.data(di_geom_11_off + 879 * ccomps * dcomps);

            auto g_y_z_xy_xxyyzz = cbuffer.data(di_geom_11_off + 880 * ccomps * dcomps);

            auto g_y_z_xy_xxyzzz = cbuffer.data(di_geom_11_off + 881 * ccomps * dcomps);

            auto g_y_z_xy_xxzzzz = cbuffer.data(di_geom_11_off + 882 * ccomps * dcomps);

            auto g_y_z_xy_xyyyyy = cbuffer.data(di_geom_11_off + 883 * ccomps * dcomps);

            auto g_y_z_xy_xyyyyz = cbuffer.data(di_geom_11_off + 884 * ccomps * dcomps);

            auto g_y_z_xy_xyyyzz = cbuffer.data(di_geom_11_off + 885 * ccomps * dcomps);

            auto g_y_z_xy_xyyzzz = cbuffer.data(di_geom_11_off + 886 * ccomps * dcomps);

            auto g_y_z_xy_xyzzzz = cbuffer.data(di_geom_11_off + 887 * ccomps * dcomps);

            auto g_y_z_xy_xzzzzz = cbuffer.data(di_geom_11_off + 888 * ccomps * dcomps);

            auto g_y_z_xy_yyyyyy = cbuffer.data(di_geom_11_off + 889 * ccomps * dcomps);

            auto g_y_z_xy_yyyyyz = cbuffer.data(di_geom_11_off + 890 * ccomps * dcomps);

            auto g_y_z_xy_yyyyzz = cbuffer.data(di_geom_11_off + 891 * ccomps * dcomps);

            auto g_y_z_xy_yyyzzz = cbuffer.data(di_geom_11_off + 892 * ccomps * dcomps);

            auto g_y_z_xy_yyzzzz = cbuffer.data(di_geom_11_off + 893 * ccomps * dcomps);

            auto g_y_z_xy_yzzzzz = cbuffer.data(di_geom_11_off + 894 * ccomps * dcomps);

            auto g_y_z_xy_zzzzzz = cbuffer.data(di_geom_11_off + 895 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xy_xxxxxx, g_y_z_xy_xxxxxy, g_y_z_xy_xxxxxz, g_y_z_xy_xxxxyy, g_y_z_xy_xxxxyz, g_y_z_xy_xxxxzz, g_y_z_xy_xxxyyy, g_y_z_xy_xxxyyz, g_y_z_xy_xxxyzz, g_y_z_xy_xxxzzz, g_y_z_xy_xxyyyy, g_y_z_xy_xxyyyz, g_y_z_xy_xxyyzz, g_y_z_xy_xxyzzz, g_y_z_xy_xxzzzz, g_y_z_xy_xyyyyy, g_y_z_xy_xyyyyz, g_y_z_xy_xyyyzz, g_y_z_xy_xyyzzz, g_y_z_xy_xyzzzz, g_y_z_xy_xzzzzz, g_y_z_xy_yyyyyy, g_y_z_xy_yyyyyz, g_y_z_xy_yyyyzz, g_y_z_xy_yyyzzz, g_y_z_xy_yyzzzz, g_y_z_xy_yzzzzz, g_y_z_xy_zzzzzz, g_y_z_y_xxxxxx, g_y_z_y_xxxxxxx, g_y_z_y_xxxxxxy, g_y_z_y_xxxxxxz, g_y_z_y_xxxxxy, g_y_z_y_xxxxxyy, g_y_z_y_xxxxxyz, g_y_z_y_xxxxxz, g_y_z_y_xxxxxzz, g_y_z_y_xxxxyy, g_y_z_y_xxxxyyy, g_y_z_y_xxxxyyz, g_y_z_y_xxxxyz, g_y_z_y_xxxxyzz, g_y_z_y_xxxxzz, g_y_z_y_xxxxzzz, g_y_z_y_xxxyyy, g_y_z_y_xxxyyyy, g_y_z_y_xxxyyyz, g_y_z_y_xxxyyz, g_y_z_y_xxxyyzz, g_y_z_y_xxxyzz, g_y_z_y_xxxyzzz, g_y_z_y_xxxzzz, g_y_z_y_xxxzzzz, g_y_z_y_xxyyyy, g_y_z_y_xxyyyyy, g_y_z_y_xxyyyyz, g_y_z_y_xxyyyz, g_y_z_y_xxyyyzz, g_y_z_y_xxyyzz, g_y_z_y_xxyyzzz, g_y_z_y_xxyzzz, g_y_z_y_xxyzzzz, g_y_z_y_xxzzzz, g_y_z_y_xxzzzzz, g_y_z_y_xyyyyy, g_y_z_y_xyyyyyy, g_y_z_y_xyyyyyz, g_y_z_y_xyyyyz, g_y_z_y_xyyyyzz, g_y_z_y_xyyyzz, g_y_z_y_xyyyzzz, g_y_z_y_xyyzzz, g_y_z_y_xyyzzzz, g_y_z_y_xyzzzz, g_y_z_y_xyzzzzz, g_y_z_y_xzzzzz, g_y_z_y_xzzzzzz, g_y_z_y_yyyyyy, g_y_z_y_yyyyyz, g_y_z_y_yyyyzz, g_y_z_y_yyyzzz, g_y_z_y_yyzzzz, g_y_z_y_yzzzzz, g_y_z_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xy_xxxxxx[k] = -g_y_z_y_xxxxxx[k] * ab_x + g_y_z_y_xxxxxxx[k];

                g_y_z_xy_xxxxxy[k] = -g_y_z_y_xxxxxy[k] * ab_x + g_y_z_y_xxxxxxy[k];

                g_y_z_xy_xxxxxz[k] = -g_y_z_y_xxxxxz[k] * ab_x + g_y_z_y_xxxxxxz[k];

                g_y_z_xy_xxxxyy[k] = -g_y_z_y_xxxxyy[k] * ab_x + g_y_z_y_xxxxxyy[k];

                g_y_z_xy_xxxxyz[k] = -g_y_z_y_xxxxyz[k] * ab_x + g_y_z_y_xxxxxyz[k];

                g_y_z_xy_xxxxzz[k] = -g_y_z_y_xxxxzz[k] * ab_x + g_y_z_y_xxxxxzz[k];

                g_y_z_xy_xxxyyy[k] = -g_y_z_y_xxxyyy[k] * ab_x + g_y_z_y_xxxxyyy[k];

                g_y_z_xy_xxxyyz[k] = -g_y_z_y_xxxyyz[k] * ab_x + g_y_z_y_xxxxyyz[k];

                g_y_z_xy_xxxyzz[k] = -g_y_z_y_xxxyzz[k] * ab_x + g_y_z_y_xxxxyzz[k];

                g_y_z_xy_xxxzzz[k] = -g_y_z_y_xxxzzz[k] * ab_x + g_y_z_y_xxxxzzz[k];

                g_y_z_xy_xxyyyy[k] = -g_y_z_y_xxyyyy[k] * ab_x + g_y_z_y_xxxyyyy[k];

                g_y_z_xy_xxyyyz[k] = -g_y_z_y_xxyyyz[k] * ab_x + g_y_z_y_xxxyyyz[k];

                g_y_z_xy_xxyyzz[k] = -g_y_z_y_xxyyzz[k] * ab_x + g_y_z_y_xxxyyzz[k];

                g_y_z_xy_xxyzzz[k] = -g_y_z_y_xxyzzz[k] * ab_x + g_y_z_y_xxxyzzz[k];

                g_y_z_xy_xxzzzz[k] = -g_y_z_y_xxzzzz[k] * ab_x + g_y_z_y_xxxzzzz[k];

                g_y_z_xy_xyyyyy[k] = -g_y_z_y_xyyyyy[k] * ab_x + g_y_z_y_xxyyyyy[k];

                g_y_z_xy_xyyyyz[k] = -g_y_z_y_xyyyyz[k] * ab_x + g_y_z_y_xxyyyyz[k];

                g_y_z_xy_xyyyzz[k] = -g_y_z_y_xyyyzz[k] * ab_x + g_y_z_y_xxyyyzz[k];

                g_y_z_xy_xyyzzz[k] = -g_y_z_y_xyyzzz[k] * ab_x + g_y_z_y_xxyyzzz[k];

                g_y_z_xy_xyzzzz[k] = -g_y_z_y_xyzzzz[k] * ab_x + g_y_z_y_xxyzzzz[k];

                g_y_z_xy_xzzzzz[k] = -g_y_z_y_xzzzzz[k] * ab_x + g_y_z_y_xxzzzzz[k];

                g_y_z_xy_yyyyyy[k] = -g_y_z_y_yyyyyy[k] * ab_x + g_y_z_y_xyyyyyy[k];

                g_y_z_xy_yyyyyz[k] = -g_y_z_y_yyyyyz[k] * ab_x + g_y_z_y_xyyyyyz[k];

                g_y_z_xy_yyyyzz[k] = -g_y_z_y_yyyyzz[k] * ab_x + g_y_z_y_xyyyyzz[k];

                g_y_z_xy_yyyzzz[k] = -g_y_z_y_yyyzzz[k] * ab_x + g_y_z_y_xyyyzzz[k];

                g_y_z_xy_yyzzzz[k] = -g_y_z_y_yyzzzz[k] * ab_x + g_y_z_y_xyyzzzz[k];

                g_y_z_xy_yzzzzz[k] = -g_y_z_y_yzzzzz[k] * ab_x + g_y_z_y_xyzzzzz[k];

                g_y_z_xy_zzzzzz[k] = -g_y_z_y_zzzzzz[k] * ab_x + g_y_z_y_xzzzzzz[k];
            }

            /// Set up 896-924 components of targeted buffer : cbuffer.data(

            auto g_y_z_xz_xxxxxx = cbuffer.data(di_geom_11_off + 896 * ccomps * dcomps);

            auto g_y_z_xz_xxxxxy = cbuffer.data(di_geom_11_off + 897 * ccomps * dcomps);

            auto g_y_z_xz_xxxxxz = cbuffer.data(di_geom_11_off + 898 * ccomps * dcomps);

            auto g_y_z_xz_xxxxyy = cbuffer.data(di_geom_11_off + 899 * ccomps * dcomps);

            auto g_y_z_xz_xxxxyz = cbuffer.data(di_geom_11_off + 900 * ccomps * dcomps);

            auto g_y_z_xz_xxxxzz = cbuffer.data(di_geom_11_off + 901 * ccomps * dcomps);

            auto g_y_z_xz_xxxyyy = cbuffer.data(di_geom_11_off + 902 * ccomps * dcomps);

            auto g_y_z_xz_xxxyyz = cbuffer.data(di_geom_11_off + 903 * ccomps * dcomps);

            auto g_y_z_xz_xxxyzz = cbuffer.data(di_geom_11_off + 904 * ccomps * dcomps);

            auto g_y_z_xz_xxxzzz = cbuffer.data(di_geom_11_off + 905 * ccomps * dcomps);

            auto g_y_z_xz_xxyyyy = cbuffer.data(di_geom_11_off + 906 * ccomps * dcomps);

            auto g_y_z_xz_xxyyyz = cbuffer.data(di_geom_11_off + 907 * ccomps * dcomps);

            auto g_y_z_xz_xxyyzz = cbuffer.data(di_geom_11_off + 908 * ccomps * dcomps);

            auto g_y_z_xz_xxyzzz = cbuffer.data(di_geom_11_off + 909 * ccomps * dcomps);

            auto g_y_z_xz_xxzzzz = cbuffer.data(di_geom_11_off + 910 * ccomps * dcomps);

            auto g_y_z_xz_xyyyyy = cbuffer.data(di_geom_11_off + 911 * ccomps * dcomps);

            auto g_y_z_xz_xyyyyz = cbuffer.data(di_geom_11_off + 912 * ccomps * dcomps);

            auto g_y_z_xz_xyyyzz = cbuffer.data(di_geom_11_off + 913 * ccomps * dcomps);

            auto g_y_z_xz_xyyzzz = cbuffer.data(di_geom_11_off + 914 * ccomps * dcomps);

            auto g_y_z_xz_xyzzzz = cbuffer.data(di_geom_11_off + 915 * ccomps * dcomps);

            auto g_y_z_xz_xzzzzz = cbuffer.data(di_geom_11_off + 916 * ccomps * dcomps);

            auto g_y_z_xz_yyyyyy = cbuffer.data(di_geom_11_off + 917 * ccomps * dcomps);

            auto g_y_z_xz_yyyyyz = cbuffer.data(di_geom_11_off + 918 * ccomps * dcomps);

            auto g_y_z_xz_yyyyzz = cbuffer.data(di_geom_11_off + 919 * ccomps * dcomps);

            auto g_y_z_xz_yyyzzz = cbuffer.data(di_geom_11_off + 920 * ccomps * dcomps);

            auto g_y_z_xz_yyzzzz = cbuffer.data(di_geom_11_off + 921 * ccomps * dcomps);

            auto g_y_z_xz_yzzzzz = cbuffer.data(di_geom_11_off + 922 * ccomps * dcomps);

            auto g_y_z_xz_zzzzzz = cbuffer.data(di_geom_11_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xz_xxxxxx, g_y_z_xz_xxxxxy, g_y_z_xz_xxxxxz, g_y_z_xz_xxxxyy, g_y_z_xz_xxxxyz, g_y_z_xz_xxxxzz, g_y_z_xz_xxxyyy, g_y_z_xz_xxxyyz, g_y_z_xz_xxxyzz, g_y_z_xz_xxxzzz, g_y_z_xz_xxyyyy, g_y_z_xz_xxyyyz, g_y_z_xz_xxyyzz, g_y_z_xz_xxyzzz, g_y_z_xz_xxzzzz, g_y_z_xz_xyyyyy, g_y_z_xz_xyyyyz, g_y_z_xz_xyyyzz, g_y_z_xz_xyyzzz, g_y_z_xz_xyzzzz, g_y_z_xz_xzzzzz, g_y_z_xz_yyyyyy, g_y_z_xz_yyyyyz, g_y_z_xz_yyyyzz, g_y_z_xz_yyyzzz, g_y_z_xz_yyzzzz, g_y_z_xz_yzzzzz, g_y_z_xz_zzzzzz, g_y_z_z_xxxxxx, g_y_z_z_xxxxxxx, g_y_z_z_xxxxxxy, g_y_z_z_xxxxxxz, g_y_z_z_xxxxxy, g_y_z_z_xxxxxyy, g_y_z_z_xxxxxyz, g_y_z_z_xxxxxz, g_y_z_z_xxxxxzz, g_y_z_z_xxxxyy, g_y_z_z_xxxxyyy, g_y_z_z_xxxxyyz, g_y_z_z_xxxxyz, g_y_z_z_xxxxyzz, g_y_z_z_xxxxzz, g_y_z_z_xxxxzzz, g_y_z_z_xxxyyy, g_y_z_z_xxxyyyy, g_y_z_z_xxxyyyz, g_y_z_z_xxxyyz, g_y_z_z_xxxyyzz, g_y_z_z_xxxyzz, g_y_z_z_xxxyzzz, g_y_z_z_xxxzzz, g_y_z_z_xxxzzzz, g_y_z_z_xxyyyy, g_y_z_z_xxyyyyy, g_y_z_z_xxyyyyz, g_y_z_z_xxyyyz, g_y_z_z_xxyyyzz, g_y_z_z_xxyyzz, g_y_z_z_xxyyzzz, g_y_z_z_xxyzzz, g_y_z_z_xxyzzzz, g_y_z_z_xxzzzz, g_y_z_z_xxzzzzz, g_y_z_z_xyyyyy, g_y_z_z_xyyyyyy, g_y_z_z_xyyyyyz, g_y_z_z_xyyyyz, g_y_z_z_xyyyyzz, g_y_z_z_xyyyzz, g_y_z_z_xyyyzzz, g_y_z_z_xyyzzz, g_y_z_z_xyyzzzz, g_y_z_z_xyzzzz, g_y_z_z_xyzzzzz, g_y_z_z_xzzzzz, g_y_z_z_xzzzzzz, g_y_z_z_yyyyyy, g_y_z_z_yyyyyz, g_y_z_z_yyyyzz, g_y_z_z_yyyzzz, g_y_z_z_yyzzzz, g_y_z_z_yzzzzz, g_y_z_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xz_xxxxxx[k] = -g_y_z_z_xxxxxx[k] * ab_x + g_y_z_z_xxxxxxx[k];

                g_y_z_xz_xxxxxy[k] = -g_y_z_z_xxxxxy[k] * ab_x + g_y_z_z_xxxxxxy[k];

                g_y_z_xz_xxxxxz[k] = -g_y_z_z_xxxxxz[k] * ab_x + g_y_z_z_xxxxxxz[k];

                g_y_z_xz_xxxxyy[k] = -g_y_z_z_xxxxyy[k] * ab_x + g_y_z_z_xxxxxyy[k];

                g_y_z_xz_xxxxyz[k] = -g_y_z_z_xxxxyz[k] * ab_x + g_y_z_z_xxxxxyz[k];

                g_y_z_xz_xxxxzz[k] = -g_y_z_z_xxxxzz[k] * ab_x + g_y_z_z_xxxxxzz[k];

                g_y_z_xz_xxxyyy[k] = -g_y_z_z_xxxyyy[k] * ab_x + g_y_z_z_xxxxyyy[k];

                g_y_z_xz_xxxyyz[k] = -g_y_z_z_xxxyyz[k] * ab_x + g_y_z_z_xxxxyyz[k];

                g_y_z_xz_xxxyzz[k] = -g_y_z_z_xxxyzz[k] * ab_x + g_y_z_z_xxxxyzz[k];

                g_y_z_xz_xxxzzz[k] = -g_y_z_z_xxxzzz[k] * ab_x + g_y_z_z_xxxxzzz[k];

                g_y_z_xz_xxyyyy[k] = -g_y_z_z_xxyyyy[k] * ab_x + g_y_z_z_xxxyyyy[k];

                g_y_z_xz_xxyyyz[k] = -g_y_z_z_xxyyyz[k] * ab_x + g_y_z_z_xxxyyyz[k];

                g_y_z_xz_xxyyzz[k] = -g_y_z_z_xxyyzz[k] * ab_x + g_y_z_z_xxxyyzz[k];

                g_y_z_xz_xxyzzz[k] = -g_y_z_z_xxyzzz[k] * ab_x + g_y_z_z_xxxyzzz[k];

                g_y_z_xz_xxzzzz[k] = -g_y_z_z_xxzzzz[k] * ab_x + g_y_z_z_xxxzzzz[k];

                g_y_z_xz_xyyyyy[k] = -g_y_z_z_xyyyyy[k] * ab_x + g_y_z_z_xxyyyyy[k];

                g_y_z_xz_xyyyyz[k] = -g_y_z_z_xyyyyz[k] * ab_x + g_y_z_z_xxyyyyz[k];

                g_y_z_xz_xyyyzz[k] = -g_y_z_z_xyyyzz[k] * ab_x + g_y_z_z_xxyyyzz[k];

                g_y_z_xz_xyyzzz[k] = -g_y_z_z_xyyzzz[k] * ab_x + g_y_z_z_xxyyzzz[k];

                g_y_z_xz_xyzzzz[k] = -g_y_z_z_xyzzzz[k] * ab_x + g_y_z_z_xxyzzzz[k];

                g_y_z_xz_xzzzzz[k] = -g_y_z_z_xzzzzz[k] * ab_x + g_y_z_z_xxzzzzz[k];

                g_y_z_xz_yyyyyy[k] = -g_y_z_z_yyyyyy[k] * ab_x + g_y_z_z_xyyyyyy[k];

                g_y_z_xz_yyyyyz[k] = -g_y_z_z_yyyyyz[k] * ab_x + g_y_z_z_xyyyyyz[k];

                g_y_z_xz_yyyyzz[k] = -g_y_z_z_yyyyzz[k] * ab_x + g_y_z_z_xyyyyzz[k];

                g_y_z_xz_yyyzzz[k] = -g_y_z_z_yyyzzz[k] * ab_x + g_y_z_z_xyyyzzz[k];

                g_y_z_xz_yyzzzz[k] = -g_y_z_z_yyzzzz[k] * ab_x + g_y_z_z_xyyzzzz[k];

                g_y_z_xz_yzzzzz[k] = -g_y_z_z_yzzzzz[k] * ab_x + g_y_z_z_xyzzzzz[k];

                g_y_z_xz_zzzzzz[k] = -g_y_z_z_zzzzzz[k] * ab_x + g_y_z_z_xzzzzzz[k];
            }

            /// Set up 924-952 components of targeted buffer : cbuffer.data(

            auto g_y_z_yy_xxxxxx = cbuffer.data(di_geom_11_off + 924 * ccomps * dcomps);

            auto g_y_z_yy_xxxxxy = cbuffer.data(di_geom_11_off + 925 * ccomps * dcomps);

            auto g_y_z_yy_xxxxxz = cbuffer.data(di_geom_11_off + 926 * ccomps * dcomps);

            auto g_y_z_yy_xxxxyy = cbuffer.data(di_geom_11_off + 927 * ccomps * dcomps);

            auto g_y_z_yy_xxxxyz = cbuffer.data(di_geom_11_off + 928 * ccomps * dcomps);

            auto g_y_z_yy_xxxxzz = cbuffer.data(di_geom_11_off + 929 * ccomps * dcomps);

            auto g_y_z_yy_xxxyyy = cbuffer.data(di_geom_11_off + 930 * ccomps * dcomps);

            auto g_y_z_yy_xxxyyz = cbuffer.data(di_geom_11_off + 931 * ccomps * dcomps);

            auto g_y_z_yy_xxxyzz = cbuffer.data(di_geom_11_off + 932 * ccomps * dcomps);

            auto g_y_z_yy_xxxzzz = cbuffer.data(di_geom_11_off + 933 * ccomps * dcomps);

            auto g_y_z_yy_xxyyyy = cbuffer.data(di_geom_11_off + 934 * ccomps * dcomps);

            auto g_y_z_yy_xxyyyz = cbuffer.data(di_geom_11_off + 935 * ccomps * dcomps);

            auto g_y_z_yy_xxyyzz = cbuffer.data(di_geom_11_off + 936 * ccomps * dcomps);

            auto g_y_z_yy_xxyzzz = cbuffer.data(di_geom_11_off + 937 * ccomps * dcomps);

            auto g_y_z_yy_xxzzzz = cbuffer.data(di_geom_11_off + 938 * ccomps * dcomps);

            auto g_y_z_yy_xyyyyy = cbuffer.data(di_geom_11_off + 939 * ccomps * dcomps);

            auto g_y_z_yy_xyyyyz = cbuffer.data(di_geom_11_off + 940 * ccomps * dcomps);

            auto g_y_z_yy_xyyyzz = cbuffer.data(di_geom_11_off + 941 * ccomps * dcomps);

            auto g_y_z_yy_xyyzzz = cbuffer.data(di_geom_11_off + 942 * ccomps * dcomps);

            auto g_y_z_yy_xyzzzz = cbuffer.data(di_geom_11_off + 943 * ccomps * dcomps);

            auto g_y_z_yy_xzzzzz = cbuffer.data(di_geom_11_off + 944 * ccomps * dcomps);

            auto g_y_z_yy_yyyyyy = cbuffer.data(di_geom_11_off + 945 * ccomps * dcomps);

            auto g_y_z_yy_yyyyyz = cbuffer.data(di_geom_11_off + 946 * ccomps * dcomps);

            auto g_y_z_yy_yyyyzz = cbuffer.data(di_geom_11_off + 947 * ccomps * dcomps);

            auto g_y_z_yy_yyyzzz = cbuffer.data(di_geom_11_off + 948 * ccomps * dcomps);

            auto g_y_z_yy_yyzzzz = cbuffer.data(di_geom_11_off + 949 * ccomps * dcomps);

            auto g_y_z_yy_yzzzzz = cbuffer.data(di_geom_11_off + 950 * ccomps * dcomps);

            auto g_y_z_yy_zzzzzz = cbuffer.data(di_geom_11_off + 951 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_y_xxxxxx, g_0_z_y_xxxxxy, g_0_z_y_xxxxxz, g_0_z_y_xxxxyy, g_0_z_y_xxxxyz, g_0_z_y_xxxxzz, g_0_z_y_xxxyyy, g_0_z_y_xxxyyz, g_0_z_y_xxxyzz, g_0_z_y_xxxzzz, g_0_z_y_xxyyyy, g_0_z_y_xxyyyz, g_0_z_y_xxyyzz, g_0_z_y_xxyzzz, g_0_z_y_xxzzzz, g_0_z_y_xyyyyy, g_0_z_y_xyyyyz, g_0_z_y_xyyyzz, g_0_z_y_xyyzzz, g_0_z_y_xyzzzz, g_0_z_y_xzzzzz, g_0_z_y_yyyyyy, g_0_z_y_yyyyyz, g_0_z_y_yyyyzz, g_0_z_y_yyyzzz, g_0_z_y_yyzzzz, g_0_z_y_yzzzzz, g_0_z_y_zzzzzz, g_y_z_y_xxxxxx, g_y_z_y_xxxxxxy, g_y_z_y_xxxxxy, g_y_z_y_xxxxxyy, g_y_z_y_xxxxxyz, g_y_z_y_xxxxxz, g_y_z_y_xxxxyy, g_y_z_y_xxxxyyy, g_y_z_y_xxxxyyz, g_y_z_y_xxxxyz, g_y_z_y_xxxxyzz, g_y_z_y_xxxxzz, g_y_z_y_xxxyyy, g_y_z_y_xxxyyyy, g_y_z_y_xxxyyyz, g_y_z_y_xxxyyz, g_y_z_y_xxxyyzz, g_y_z_y_xxxyzz, g_y_z_y_xxxyzzz, g_y_z_y_xxxzzz, g_y_z_y_xxyyyy, g_y_z_y_xxyyyyy, g_y_z_y_xxyyyyz, g_y_z_y_xxyyyz, g_y_z_y_xxyyyzz, g_y_z_y_xxyyzz, g_y_z_y_xxyyzzz, g_y_z_y_xxyzzz, g_y_z_y_xxyzzzz, g_y_z_y_xxzzzz, g_y_z_y_xyyyyy, g_y_z_y_xyyyyyy, g_y_z_y_xyyyyyz, g_y_z_y_xyyyyz, g_y_z_y_xyyyyzz, g_y_z_y_xyyyzz, g_y_z_y_xyyyzzz, g_y_z_y_xyyzzz, g_y_z_y_xyyzzzz, g_y_z_y_xyzzzz, g_y_z_y_xyzzzzz, g_y_z_y_xzzzzz, g_y_z_y_yyyyyy, g_y_z_y_yyyyyyy, g_y_z_y_yyyyyyz, g_y_z_y_yyyyyz, g_y_z_y_yyyyyzz, g_y_z_y_yyyyzz, g_y_z_y_yyyyzzz, g_y_z_y_yyyzzz, g_y_z_y_yyyzzzz, g_y_z_y_yyzzzz, g_y_z_y_yyzzzzz, g_y_z_y_yzzzzz, g_y_z_y_yzzzzzz, g_y_z_y_zzzzzz, g_y_z_yy_xxxxxx, g_y_z_yy_xxxxxy, g_y_z_yy_xxxxxz, g_y_z_yy_xxxxyy, g_y_z_yy_xxxxyz, g_y_z_yy_xxxxzz, g_y_z_yy_xxxyyy, g_y_z_yy_xxxyyz, g_y_z_yy_xxxyzz, g_y_z_yy_xxxzzz, g_y_z_yy_xxyyyy, g_y_z_yy_xxyyyz, g_y_z_yy_xxyyzz, g_y_z_yy_xxyzzz, g_y_z_yy_xxzzzz, g_y_z_yy_xyyyyy, g_y_z_yy_xyyyyz, g_y_z_yy_xyyyzz, g_y_z_yy_xyyzzz, g_y_z_yy_xyzzzz, g_y_z_yy_xzzzzz, g_y_z_yy_yyyyyy, g_y_z_yy_yyyyyz, g_y_z_yy_yyyyzz, g_y_z_yy_yyyzzz, g_y_z_yy_yyzzzz, g_y_z_yy_yzzzzz, g_y_z_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yy_xxxxxx[k] = -g_0_z_y_xxxxxx[k] - g_y_z_y_xxxxxx[k] * ab_y + g_y_z_y_xxxxxxy[k];

                g_y_z_yy_xxxxxy[k] = -g_0_z_y_xxxxxy[k] - g_y_z_y_xxxxxy[k] * ab_y + g_y_z_y_xxxxxyy[k];

                g_y_z_yy_xxxxxz[k] = -g_0_z_y_xxxxxz[k] - g_y_z_y_xxxxxz[k] * ab_y + g_y_z_y_xxxxxyz[k];

                g_y_z_yy_xxxxyy[k] = -g_0_z_y_xxxxyy[k] - g_y_z_y_xxxxyy[k] * ab_y + g_y_z_y_xxxxyyy[k];

                g_y_z_yy_xxxxyz[k] = -g_0_z_y_xxxxyz[k] - g_y_z_y_xxxxyz[k] * ab_y + g_y_z_y_xxxxyyz[k];

                g_y_z_yy_xxxxzz[k] = -g_0_z_y_xxxxzz[k] - g_y_z_y_xxxxzz[k] * ab_y + g_y_z_y_xxxxyzz[k];

                g_y_z_yy_xxxyyy[k] = -g_0_z_y_xxxyyy[k] - g_y_z_y_xxxyyy[k] * ab_y + g_y_z_y_xxxyyyy[k];

                g_y_z_yy_xxxyyz[k] = -g_0_z_y_xxxyyz[k] - g_y_z_y_xxxyyz[k] * ab_y + g_y_z_y_xxxyyyz[k];

                g_y_z_yy_xxxyzz[k] = -g_0_z_y_xxxyzz[k] - g_y_z_y_xxxyzz[k] * ab_y + g_y_z_y_xxxyyzz[k];

                g_y_z_yy_xxxzzz[k] = -g_0_z_y_xxxzzz[k] - g_y_z_y_xxxzzz[k] * ab_y + g_y_z_y_xxxyzzz[k];

                g_y_z_yy_xxyyyy[k] = -g_0_z_y_xxyyyy[k] - g_y_z_y_xxyyyy[k] * ab_y + g_y_z_y_xxyyyyy[k];

                g_y_z_yy_xxyyyz[k] = -g_0_z_y_xxyyyz[k] - g_y_z_y_xxyyyz[k] * ab_y + g_y_z_y_xxyyyyz[k];

                g_y_z_yy_xxyyzz[k] = -g_0_z_y_xxyyzz[k] - g_y_z_y_xxyyzz[k] * ab_y + g_y_z_y_xxyyyzz[k];

                g_y_z_yy_xxyzzz[k] = -g_0_z_y_xxyzzz[k] - g_y_z_y_xxyzzz[k] * ab_y + g_y_z_y_xxyyzzz[k];

                g_y_z_yy_xxzzzz[k] = -g_0_z_y_xxzzzz[k] - g_y_z_y_xxzzzz[k] * ab_y + g_y_z_y_xxyzzzz[k];

                g_y_z_yy_xyyyyy[k] = -g_0_z_y_xyyyyy[k] - g_y_z_y_xyyyyy[k] * ab_y + g_y_z_y_xyyyyyy[k];

                g_y_z_yy_xyyyyz[k] = -g_0_z_y_xyyyyz[k] - g_y_z_y_xyyyyz[k] * ab_y + g_y_z_y_xyyyyyz[k];

                g_y_z_yy_xyyyzz[k] = -g_0_z_y_xyyyzz[k] - g_y_z_y_xyyyzz[k] * ab_y + g_y_z_y_xyyyyzz[k];

                g_y_z_yy_xyyzzz[k] = -g_0_z_y_xyyzzz[k] - g_y_z_y_xyyzzz[k] * ab_y + g_y_z_y_xyyyzzz[k];

                g_y_z_yy_xyzzzz[k] = -g_0_z_y_xyzzzz[k] - g_y_z_y_xyzzzz[k] * ab_y + g_y_z_y_xyyzzzz[k];

                g_y_z_yy_xzzzzz[k] = -g_0_z_y_xzzzzz[k] - g_y_z_y_xzzzzz[k] * ab_y + g_y_z_y_xyzzzzz[k];

                g_y_z_yy_yyyyyy[k] = -g_0_z_y_yyyyyy[k] - g_y_z_y_yyyyyy[k] * ab_y + g_y_z_y_yyyyyyy[k];

                g_y_z_yy_yyyyyz[k] = -g_0_z_y_yyyyyz[k] - g_y_z_y_yyyyyz[k] * ab_y + g_y_z_y_yyyyyyz[k];

                g_y_z_yy_yyyyzz[k] = -g_0_z_y_yyyyzz[k] - g_y_z_y_yyyyzz[k] * ab_y + g_y_z_y_yyyyyzz[k];

                g_y_z_yy_yyyzzz[k] = -g_0_z_y_yyyzzz[k] - g_y_z_y_yyyzzz[k] * ab_y + g_y_z_y_yyyyzzz[k];

                g_y_z_yy_yyzzzz[k] = -g_0_z_y_yyzzzz[k] - g_y_z_y_yyzzzz[k] * ab_y + g_y_z_y_yyyzzzz[k];

                g_y_z_yy_yzzzzz[k] = -g_0_z_y_yzzzzz[k] - g_y_z_y_yzzzzz[k] * ab_y + g_y_z_y_yyzzzzz[k];

                g_y_z_yy_zzzzzz[k] = -g_0_z_y_zzzzzz[k] - g_y_z_y_zzzzzz[k] * ab_y + g_y_z_y_yzzzzzz[k];
            }

            /// Set up 952-980 components of targeted buffer : cbuffer.data(

            auto g_y_z_yz_xxxxxx = cbuffer.data(di_geom_11_off + 952 * ccomps * dcomps);

            auto g_y_z_yz_xxxxxy = cbuffer.data(di_geom_11_off + 953 * ccomps * dcomps);

            auto g_y_z_yz_xxxxxz = cbuffer.data(di_geom_11_off + 954 * ccomps * dcomps);

            auto g_y_z_yz_xxxxyy = cbuffer.data(di_geom_11_off + 955 * ccomps * dcomps);

            auto g_y_z_yz_xxxxyz = cbuffer.data(di_geom_11_off + 956 * ccomps * dcomps);

            auto g_y_z_yz_xxxxzz = cbuffer.data(di_geom_11_off + 957 * ccomps * dcomps);

            auto g_y_z_yz_xxxyyy = cbuffer.data(di_geom_11_off + 958 * ccomps * dcomps);

            auto g_y_z_yz_xxxyyz = cbuffer.data(di_geom_11_off + 959 * ccomps * dcomps);

            auto g_y_z_yz_xxxyzz = cbuffer.data(di_geom_11_off + 960 * ccomps * dcomps);

            auto g_y_z_yz_xxxzzz = cbuffer.data(di_geom_11_off + 961 * ccomps * dcomps);

            auto g_y_z_yz_xxyyyy = cbuffer.data(di_geom_11_off + 962 * ccomps * dcomps);

            auto g_y_z_yz_xxyyyz = cbuffer.data(di_geom_11_off + 963 * ccomps * dcomps);

            auto g_y_z_yz_xxyyzz = cbuffer.data(di_geom_11_off + 964 * ccomps * dcomps);

            auto g_y_z_yz_xxyzzz = cbuffer.data(di_geom_11_off + 965 * ccomps * dcomps);

            auto g_y_z_yz_xxzzzz = cbuffer.data(di_geom_11_off + 966 * ccomps * dcomps);

            auto g_y_z_yz_xyyyyy = cbuffer.data(di_geom_11_off + 967 * ccomps * dcomps);

            auto g_y_z_yz_xyyyyz = cbuffer.data(di_geom_11_off + 968 * ccomps * dcomps);

            auto g_y_z_yz_xyyyzz = cbuffer.data(di_geom_11_off + 969 * ccomps * dcomps);

            auto g_y_z_yz_xyyzzz = cbuffer.data(di_geom_11_off + 970 * ccomps * dcomps);

            auto g_y_z_yz_xyzzzz = cbuffer.data(di_geom_11_off + 971 * ccomps * dcomps);

            auto g_y_z_yz_xzzzzz = cbuffer.data(di_geom_11_off + 972 * ccomps * dcomps);

            auto g_y_z_yz_yyyyyy = cbuffer.data(di_geom_11_off + 973 * ccomps * dcomps);

            auto g_y_z_yz_yyyyyz = cbuffer.data(di_geom_11_off + 974 * ccomps * dcomps);

            auto g_y_z_yz_yyyyzz = cbuffer.data(di_geom_11_off + 975 * ccomps * dcomps);

            auto g_y_z_yz_yyyzzz = cbuffer.data(di_geom_11_off + 976 * ccomps * dcomps);

            auto g_y_z_yz_yyzzzz = cbuffer.data(di_geom_11_off + 977 * ccomps * dcomps);

            auto g_y_z_yz_yzzzzz = cbuffer.data(di_geom_11_off + 978 * ccomps * dcomps);

            auto g_y_z_yz_zzzzzz = cbuffer.data(di_geom_11_off + 979 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxxxx, g_0_z_z_xxxxxy, g_0_z_z_xxxxxz, g_0_z_z_xxxxyy, g_0_z_z_xxxxyz, g_0_z_z_xxxxzz, g_0_z_z_xxxyyy, g_0_z_z_xxxyyz, g_0_z_z_xxxyzz, g_0_z_z_xxxzzz, g_0_z_z_xxyyyy, g_0_z_z_xxyyyz, g_0_z_z_xxyyzz, g_0_z_z_xxyzzz, g_0_z_z_xxzzzz, g_0_z_z_xyyyyy, g_0_z_z_xyyyyz, g_0_z_z_xyyyzz, g_0_z_z_xyyzzz, g_0_z_z_xyzzzz, g_0_z_z_xzzzzz, g_0_z_z_yyyyyy, g_0_z_z_yyyyyz, g_0_z_z_yyyyzz, g_0_z_z_yyyzzz, g_0_z_z_yyzzzz, g_0_z_z_yzzzzz, g_0_z_z_zzzzzz, g_y_z_yz_xxxxxx, g_y_z_yz_xxxxxy, g_y_z_yz_xxxxxz, g_y_z_yz_xxxxyy, g_y_z_yz_xxxxyz, g_y_z_yz_xxxxzz, g_y_z_yz_xxxyyy, g_y_z_yz_xxxyyz, g_y_z_yz_xxxyzz, g_y_z_yz_xxxzzz, g_y_z_yz_xxyyyy, g_y_z_yz_xxyyyz, g_y_z_yz_xxyyzz, g_y_z_yz_xxyzzz, g_y_z_yz_xxzzzz, g_y_z_yz_xyyyyy, g_y_z_yz_xyyyyz, g_y_z_yz_xyyyzz, g_y_z_yz_xyyzzz, g_y_z_yz_xyzzzz, g_y_z_yz_xzzzzz, g_y_z_yz_yyyyyy, g_y_z_yz_yyyyyz, g_y_z_yz_yyyyzz, g_y_z_yz_yyyzzz, g_y_z_yz_yyzzzz, g_y_z_yz_yzzzzz, g_y_z_yz_zzzzzz, g_y_z_z_xxxxxx, g_y_z_z_xxxxxxy, g_y_z_z_xxxxxy, g_y_z_z_xxxxxyy, g_y_z_z_xxxxxyz, g_y_z_z_xxxxxz, g_y_z_z_xxxxyy, g_y_z_z_xxxxyyy, g_y_z_z_xxxxyyz, g_y_z_z_xxxxyz, g_y_z_z_xxxxyzz, g_y_z_z_xxxxzz, g_y_z_z_xxxyyy, g_y_z_z_xxxyyyy, g_y_z_z_xxxyyyz, g_y_z_z_xxxyyz, g_y_z_z_xxxyyzz, g_y_z_z_xxxyzz, g_y_z_z_xxxyzzz, g_y_z_z_xxxzzz, g_y_z_z_xxyyyy, g_y_z_z_xxyyyyy, g_y_z_z_xxyyyyz, g_y_z_z_xxyyyz, g_y_z_z_xxyyyzz, g_y_z_z_xxyyzz, g_y_z_z_xxyyzzz, g_y_z_z_xxyzzz, g_y_z_z_xxyzzzz, g_y_z_z_xxzzzz, g_y_z_z_xyyyyy, g_y_z_z_xyyyyyy, g_y_z_z_xyyyyyz, g_y_z_z_xyyyyz, g_y_z_z_xyyyyzz, g_y_z_z_xyyyzz, g_y_z_z_xyyyzzz, g_y_z_z_xyyzzz, g_y_z_z_xyyzzzz, g_y_z_z_xyzzzz, g_y_z_z_xyzzzzz, g_y_z_z_xzzzzz, g_y_z_z_yyyyyy, g_y_z_z_yyyyyyy, g_y_z_z_yyyyyyz, g_y_z_z_yyyyyz, g_y_z_z_yyyyyzz, g_y_z_z_yyyyzz, g_y_z_z_yyyyzzz, g_y_z_z_yyyzzz, g_y_z_z_yyyzzzz, g_y_z_z_yyzzzz, g_y_z_z_yyzzzzz, g_y_z_z_yzzzzz, g_y_z_z_yzzzzzz, g_y_z_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yz_xxxxxx[k] = -g_0_z_z_xxxxxx[k] - g_y_z_z_xxxxxx[k] * ab_y + g_y_z_z_xxxxxxy[k];

                g_y_z_yz_xxxxxy[k] = -g_0_z_z_xxxxxy[k] - g_y_z_z_xxxxxy[k] * ab_y + g_y_z_z_xxxxxyy[k];

                g_y_z_yz_xxxxxz[k] = -g_0_z_z_xxxxxz[k] - g_y_z_z_xxxxxz[k] * ab_y + g_y_z_z_xxxxxyz[k];

                g_y_z_yz_xxxxyy[k] = -g_0_z_z_xxxxyy[k] - g_y_z_z_xxxxyy[k] * ab_y + g_y_z_z_xxxxyyy[k];

                g_y_z_yz_xxxxyz[k] = -g_0_z_z_xxxxyz[k] - g_y_z_z_xxxxyz[k] * ab_y + g_y_z_z_xxxxyyz[k];

                g_y_z_yz_xxxxzz[k] = -g_0_z_z_xxxxzz[k] - g_y_z_z_xxxxzz[k] * ab_y + g_y_z_z_xxxxyzz[k];

                g_y_z_yz_xxxyyy[k] = -g_0_z_z_xxxyyy[k] - g_y_z_z_xxxyyy[k] * ab_y + g_y_z_z_xxxyyyy[k];

                g_y_z_yz_xxxyyz[k] = -g_0_z_z_xxxyyz[k] - g_y_z_z_xxxyyz[k] * ab_y + g_y_z_z_xxxyyyz[k];

                g_y_z_yz_xxxyzz[k] = -g_0_z_z_xxxyzz[k] - g_y_z_z_xxxyzz[k] * ab_y + g_y_z_z_xxxyyzz[k];

                g_y_z_yz_xxxzzz[k] = -g_0_z_z_xxxzzz[k] - g_y_z_z_xxxzzz[k] * ab_y + g_y_z_z_xxxyzzz[k];

                g_y_z_yz_xxyyyy[k] = -g_0_z_z_xxyyyy[k] - g_y_z_z_xxyyyy[k] * ab_y + g_y_z_z_xxyyyyy[k];

                g_y_z_yz_xxyyyz[k] = -g_0_z_z_xxyyyz[k] - g_y_z_z_xxyyyz[k] * ab_y + g_y_z_z_xxyyyyz[k];

                g_y_z_yz_xxyyzz[k] = -g_0_z_z_xxyyzz[k] - g_y_z_z_xxyyzz[k] * ab_y + g_y_z_z_xxyyyzz[k];

                g_y_z_yz_xxyzzz[k] = -g_0_z_z_xxyzzz[k] - g_y_z_z_xxyzzz[k] * ab_y + g_y_z_z_xxyyzzz[k];

                g_y_z_yz_xxzzzz[k] = -g_0_z_z_xxzzzz[k] - g_y_z_z_xxzzzz[k] * ab_y + g_y_z_z_xxyzzzz[k];

                g_y_z_yz_xyyyyy[k] = -g_0_z_z_xyyyyy[k] - g_y_z_z_xyyyyy[k] * ab_y + g_y_z_z_xyyyyyy[k];

                g_y_z_yz_xyyyyz[k] = -g_0_z_z_xyyyyz[k] - g_y_z_z_xyyyyz[k] * ab_y + g_y_z_z_xyyyyyz[k];

                g_y_z_yz_xyyyzz[k] = -g_0_z_z_xyyyzz[k] - g_y_z_z_xyyyzz[k] * ab_y + g_y_z_z_xyyyyzz[k];

                g_y_z_yz_xyyzzz[k] = -g_0_z_z_xyyzzz[k] - g_y_z_z_xyyzzz[k] * ab_y + g_y_z_z_xyyyzzz[k];

                g_y_z_yz_xyzzzz[k] = -g_0_z_z_xyzzzz[k] - g_y_z_z_xyzzzz[k] * ab_y + g_y_z_z_xyyzzzz[k];

                g_y_z_yz_xzzzzz[k] = -g_0_z_z_xzzzzz[k] - g_y_z_z_xzzzzz[k] * ab_y + g_y_z_z_xyzzzzz[k];

                g_y_z_yz_yyyyyy[k] = -g_0_z_z_yyyyyy[k] - g_y_z_z_yyyyyy[k] * ab_y + g_y_z_z_yyyyyyy[k];

                g_y_z_yz_yyyyyz[k] = -g_0_z_z_yyyyyz[k] - g_y_z_z_yyyyyz[k] * ab_y + g_y_z_z_yyyyyyz[k];

                g_y_z_yz_yyyyzz[k] = -g_0_z_z_yyyyzz[k] - g_y_z_z_yyyyzz[k] * ab_y + g_y_z_z_yyyyyzz[k];

                g_y_z_yz_yyyzzz[k] = -g_0_z_z_yyyzzz[k] - g_y_z_z_yyyzzz[k] * ab_y + g_y_z_z_yyyyzzz[k];

                g_y_z_yz_yyzzzz[k] = -g_0_z_z_yyzzzz[k] - g_y_z_z_yyzzzz[k] * ab_y + g_y_z_z_yyyzzzz[k];

                g_y_z_yz_yzzzzz[k] = -g_0_z_z_yzzzzz[k] - g_y_z_z_yzzzzz[k] * ab_y + g_y_z_z_yyzzzzz[k];

                g_y_z_yz_zzzzzz[k] = -g_0_z_z_zzzzzz[k] - g_y_z_z_zzzzzz[k] * ab_y + g_y_z_z_yzzzzzz[k];
            }

            /// Set up 980-1008 components of targeted buffer : cbuffer.data(

            auto g_y_z_zz_xxxxxx = cbuffer.data(di_geom_11_off + 980 * ccomps * dcomps);

            auto g_y_z_zz_xxxxxy = cbuffer.data(di_geom_11_off + 981 * ccomps * dcomps);

            auto g_y_z_zz_xxxxxz = cbuffer.data(di_geom_11_off + 982 * ccomps * dcomps);

            auto g_y_z_zz_xxxxyy = cbuffer.data(di_geom_11_off + 983 * ccomps * dcomps);

            auto g_y_z_zz_xxxxyz = cbuffer.data(di_geom_11_off + 984 * ccomps * dcomps);

            auto g_y_z_zz_xxxxzz = cbuffer.data(di_geom_11_off + 985 * ccomps * dcomps);

            auto g_y_z_zz_xxxyyy = cbuffer.data(di_geom_11_off + 986 * ccomps * dcomps);

            auto g_y_z_zz_xxxyyz = cbuffer.data(di_geom_11_off + 987 * ccomps * dcomps);

            auto g_y_z_zz_xxxyzz = cbuffer.data(di_geom_11_off + 988 * ccomps * dcomps);

            auto g_y_z_zz_xxxzzz = cbuffer.data(di_geom_11_off + 989 * ccomps * dcomps);

            auto g_y_z_zz_xxyyyy = cbuffer.data(di_geom_11_off + 990 * ccomps * dcomps);

            auto g_y_z_zz_xxyyyz = cbuffer.data(di_geom_11_off + 991 * ccomps * dcomps);

            auto g_y_z_zz_xxyyzz = cbuffer.data(di_geom_11_off + 992 * ccomps * dcomps);

            auto g_y_z_zz_xxyzzz = cbuffer.data(di_geom_11_off + 993 * ccomps * dcomps);

            auto g_y_z_zz_xxzzzz = cbuffer.data(di_geom_11_off + 994 * ccomps * dcomps);

            auto g_y_z_zz_xyyyyy = cbuffer.data(di_geom_11_off + 995 * ccomps * dcomps);

            auto g_y_z_zz_xyyyyz = cbuffer.data(di_geom_11_off + 996 * ccomps * dcomps);

            auto g_y_z_zz_xyyyzz = cbuffer.data(di_geom_11_off + 997 * ccomps * dcomps);

            auto g_y_z_zz_xyyzzz = cbuffer.data(di_geom_11_off + 998 * ccomps * dcomps);

            auto g_y_z_zz_xyzzzz = cbuffer.data(di_geom_11_off + 999 * ccomps * dcomps);

            auto g_y_z_zz_xzzzzz = cbuffer.data(di_geom_11_off + 1000 * ccomps * dcomps);

            auto g_y_z_zz_yyyyyy = cbuffer.data(di_geom_11_off + 1001 * ccomps * dcomps);

            auto g_y_z_zz_yyyyyz = cbuffer.data(di_geom_11_off + 1002 * ccomps * dcomps);

            auto g_y_z_zz_yyyyzz = cbuffer.data(di_geom_11_off + 1003 * ccomps * dcomps);

            auto g_y_z_zz_yyyzzz = cbuffer.data(di_geom_11_off + 1004 * ccomps * dcomps);

            auto g_y_z_zz_yyzzzz = cbuffer.data(di_geom_11_off + 1005 * ccomps * dcomps);

            auto g_y_z_zz_yzzzzz = cbuffer.data(di_geom_11_off + 1006 * ccomps * dcomps);

            auto g_y_z_zz_zzzzzz = cbuffer.data(di_geom_11_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxxxxx, g_y_0_z_xxxxxy, g_y_0_z_xxxxxz, g_y_0_z_xxxxyy, g_y_0_z_xxxxyz, g_y_0_z_xxxxzz, g_y_0_z_xxxyyy, g_y_0_z_xxxyyz, g_y_0_z_xxxyzz, g_y_0_z_xxxzzz, g_y_0_z_xxyyyy, g_y_0_z_xxyyyz, g_y_0_z_xxyyzz, g_y_0_z_xxyzzz, g_y_0_z_xxzzzz, g_y_0_z_xyyyyy, g_y_0_z_xyyyyz, g_y_0_z_xyyyzz, g_y_0_z_xyyzzz, g_y_0_z_xyzzzz, g_y_0_z_xzzzzz, g_y_0_z_yyyyyy, g_y_0_z_yyyyyz, g_y_0_z_yyyyzz, g_y_0_z_yyyzzz, g_y_0_z_yyzzzz, g_y_0_z_yzzzzz, g_y_0_z_zzzzzz, g_y_z_z_xxxxxx, g_y_z_z_xxxxxxz, g_y_z_z_xxxxxy, g_y_z_z_xxxxxyz, g_y_z_z_xxxxxz, g_y_z_z_xxxxxzz, g_y_z_z_xxxxyy, g_y_z_z_xxxxyyz, g_y_z_z_xxxxyz, g_y_z_z_xxxxyzz, g_y_z_z_xxxxzz, g_y_z_z_xxxxzzz, g_y_z_z_xxxyyy, g_y_z_z_xxxyyyz, g_y_z_z_xxxyyz, g_y_z_z_xxxyyzz, g_y_z_z_xxxyzz, g_y_z_z_xxxyzzz, g_y_z_z_xxxzzz, g_y_z_z_xxxzzzz, g_y_z_z_xxyyyy, g_y_z_z_xxyyyyz, g_y_z_z_xxyyyz, g_y_z_z_xxyyyzz, g_y_z_z_xxyyzz, g_y_z_z_xxyyzzz, g_y_z_z_xxyzzz, g_y_z_z_xxyzzzz, g_y_z_z_xxzzzz, g_y_z_z_xxzzzzz, g_y_z_z_xyyyyy, g_y_z_z_xyyyyyz, g_y_z_z_xyyyyz, g_y_z_z_xyyyyzz, g_y_z_z_xyyyzz, g_y_z_z_xyyyzzz, g_y_z_z_xyyzzz, g_y_z_z_xyyzzzz, g_y_z_z_xyzzzz, g_y_z_z_xyzzzzz, g_y_z_z_xzzzzz, g_y_z_z_xzzzzzz, g_y_z_z_yyyyyy, g_y_z_z_yyyyyyz, g_y_z_z_yyyyyz, g_y_z_z_yyyyyzz, g_y_z_z_yyyyzz, g_y_z_z_yyyyzzz, g_y_z_z_yyyzzz, g_y_z_z_yyyzzzz, g_y_z_z_yyzzzz, g_y_z_z_yyzzzzz, g_y_z_z_yzzzzz, g_y_z_z_yzzzzzz, g_y_z_z_zzzzzz, g_y_z_z_zzzzzzz, g_y_z_zz_xxxxxx, g_y_z_zz_xxxxxy, g_y_z_zz_xxxxxz, g_y_z_zz_xxxxyy, g_y_z_zz_xxxxyz, g_y_z_zz_xxxxzz, g_y_z_zz_xxxyyy, g_y_z_zz_xxxyyz, g_y_z_zz_xxxyzz, g_y_z_zz_xxxzzz, g_y_z_zz_xxyyyy, g_y_z_zz_xxyyyz, g_y_z_zz_xxyyzz, g_y_z_zz_xxyzzz, g_y_z_zz_xxzzzz, g_y_z_zz_xyyyyy, g_y_z_zz_xyyyyz, g_y_z_zz_xyyyzz, g_y_z_zz_xyyzzz, g_y_z_zz_xyzzzz, g_y_z_zz_xzzzzz, g_y_z_zz_yyyyyy, g_y_z_zz_yyyyyz, g_y_z_zz_yyyyzz, g_y_z_zz_yyyzzz, g_y_z_zz_yyzzzz, g_y_z_zz_yzzzzz, g_y_z_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zz_xxxxxx[k] = g_y_0_z_xxxxxx[k] - g_y_z_z_xxxxxx[k] * ab_z + g_y_z_z_xxxxxxz[k];

                g_y_z_zz_xxxxxy[k] = g_y_0_z_xxxxxy[k] - g_y_z_z_xxxxxy[k] * ab_z + g_y_z_z_xxxxxyz[k];

                g_y_z_zz_xxxxxz[k] = g_y_0_z_xxxxxz[k] - g_y_z_z_xxxxxz[k] * ab_z + g_y_z_z_xxxxxzz[k];

                g_y_z_zz_xxxxyy[k] = g_y_0_z_xxxxyy[k] - g_y_z_z_xxxxyy[k] * ab_z + g_y_z_z_xxxxyyz[k];

                g_y_z_zz_xxxxyz[k] = g_y_0_z_xxxxyz[k] - g_y_z_z_xxxxyz[k] * ab_z + g_y_z_z_xxxxyzz[k];

                g_y_z_zz_xxxxzz[k] = g_y_0_z_xxxxzz[k] - g_y_z_z_xxxxzz[k] * ab_z + g_y_z_z_xxxxzzz[k];

                g_y_z_zz_xxxyyy[k] = g_y_0_z_xxxyyy[k] - g_y_z_z_xxxyyy[k] * ab_z + g_y_z_z_xxxyyyz[k];

                g_y_z_zz_xxxyyz[k] = g_y_0_z_xxxyyz[k] - g_y_z_z_xxxyyz[k] * ab_z + g_y_z_z_xxxyyzz[k];

                g_y_z_zz_xxxyzz[k] = g_y_0_z_xxxyzz[k] - g_y_z_z_xxxyzz[k] * ab_z + g_y_z_z_xxxyzzz[k];

                g_y_z_zz_xxxzzz[k] = g_y_0_z_xxxzzz[k] - g_y_z_z_xxxzzz[k] * ab_z + g_y_z_z_xxxzzzz[k];

                g_y_z_zz_xxyyyy[k] = g_y_0_z_xxyyyy[k] - g_y_z_z_xxyyyy[k] * ab_z + g_y_z_z_xxyyyyz[k];

                g_y_z_zz_xxyyyz[k] = g_y_0_z_xxyyyz[k] - g_y_z_z_xxyyyz[k] * ab_z + g_y_z_z_xxyyyzz[k];

                g_y_z_zz_xxyyzz[k] = g_y_0_z_xxyyzz[k] - g_y_z_z_xxyyzz[k] * ab_z + g_y_z_z_xxyyzzz[k];

                g_y_z_zz_xxyzzz[k] = g_y_0_z_xxyzzz[k] - g_y_z_z_xxyzzz[k] * ab_z + g_y_z_z_xxyzzzz[k];

                g_y_z_zz_xxzzzz[k] = g_y_0_z_xxzzzz[k] - g_y_z_z_xxzzzz[k] * ab_z + g_y_z_z_xxzzzzz[k];

                g_y_z_zz_xyyyyy[k] = g_y_0_z_xyyyyy[k] - g_y_z_z_xyyyyy[k] * ab_z + g_y_z_z_xyyyyyz[k];

                g_y_z_zz_xyyyyz[k] = g_y_0_z_xyyyyz[k] - g_y_z_z_xyyyyz[k] * ab_z + g_y_z_z_xyyyyzz[k];

                g_y_z_zz_xyyyzz[k] = g_y_0_z_xyyyzz[k] - g_y_z_z_xyyyzz[k] * ab_z + g_y_z_z_xyyyzzz[k];

                g_y_z_zz_xyyzzz[k] = g_y_0_z_xyyzzz[k] - g_y_z_z_xyyzzz[k] * ab_z + g_y_z_z_xyyzzzz[k];

                g_y_z_zz_xyzzzz[k] = g_y_0_z_xyzzzz[k] - g_y_z_z_xyzzzz[k] * ab_z + g_y_z_z_xyzzzzz[k];

                g_y_z_zz_xzzzzz[k] = g_y_0_z_xzzzzz[k] - g_y_z_z_xzzzzz[k] * ab_z + g_y_z_z_xzzzzzz[k];

                g_y_z_zz_yyyyyy[k] = g_y_0_z_yyyyyy[k] - g_y_z_z_yyyyyy[k] * ab_z + g_y_z_z_yyyyyyz[k];

                g_y_z_zz_yyyyyz[k] = g_y_0_z_yyyyyz[k] - g_y_z_z_yyyyyz[k] * ab_z + g_y_z_z_yyyyyzz[k];

                g_y_z_zz_yyyyzz[k] = g_y_0_z_yyyyzz[k] - g_y_z_z_yyyyzz[k] * ab_z + g_y_z_z_yyyyzzz[k];

                g_y_z_zz_yyyzzz[k] = g_y_0_z_yyyzzz[k] - g_y_z_z_yyyzzz[k] * ab_z + g_y_z_z_yyyzzzz[k];

                g_y_z_zz_yyzzzz[k] = g_y_0_z_yyzzzz[k] - g_y_z_z_yyzzzz[k] * ab_z + g_y_z_z_yyzzzzz[k];

                g_y_z_zz_yzzzzz[k] = g_y_0_z_yzzzzz[k] - g_y_z_z_yzzzzz[k] * ab_z + g_y_z_z_yzzzzzz[k];

                g_y_z_zz_zzzzzz[k] = g_y_0_z_zzzzzz[k] - g_y_z_z_zzzzzz[k] * ab_z + g_y_z_z_zzzzzzz[k];
            }

            /// Set up 1008-1036 components of targeted buffer : cbuffer.data(

            auto g_z_x_xx_xxxxxx = cbuffer.data(di_geom_11_off + 1008 * ccomps * dcomps);

            auto g_z_x_xx_xxxxxy = cbuffer.data(di_geom_11_off + 1009 * ccomps * dcomps);

            auto g_z_x_xx_xxxxxz = cbuffer.data(di_geom_11_off + 1010 * ccomps * dcomps);

            auto g_z_x_xx_xxxxyy = cbuffer.data(di_geom_11_off + 1011 * ccomps * dcomps);

            auto g_z_x_xx_xxxxyz = cbuffer.data(di_geom_11_off + 1012 * ccomps * dcomps);

            auto g_z_x_xx_xxxxzz = cbuffer.data(di_geom_11_off + 1013 * ccomps * dcomps);

            auto g_z_x_xx_xxxyyy = cbuffer.data(di_geom_11_off + 1014 * ccomps * dcomps);

            auto g_z_x_xx_xxxyyz = cbuffer.data(di_geom_11_off + 1015 * ccomps * dcomps);

            auto g_z_x_xx_xxxyzz = cbuffer.data(di_geom_11_off + 1016 * ccomps * dcomps);

            auto g_z_x_xx_xxxzzz = cbuffer.data(di_geom_11_off + 1017 * ccomps * dcomps);

            auto g_z_x_xx_xxyyyy = cbuffer.data(di_geom_11_off + 1018 * ccomps * dcomps);

            auto g_z_x_xx_xxyyyz = cbuffer.data(di_geom_11_off + 1019 * ccomps * dcomps);

            auto g_z_x_xx_xxyyzz = cbuffer.data(di_geom_11_off + 1020 * ccomps * dcomps);

            auto g_z_x_xx_xxyzzz = cbuffer.data(di_geom_11_off + 1021 * ccomps * dcomps);

            auto g_z_x_xx_xxzzzz = cbuffer.data(di_geom_11_off + 1022 * ccomps * dcomps);

            auto g_z_x_xx_xyyyyy = cbuffer.data(di_geom_11_off + 1023 * ccomps * dcomps);

            auto g_z_x_xx_xyyyyz = cbuffer.data(di_geom_11_off + 1024 * ccomps * dcomps);

            auto g_z_x_xx_xyyyzz = cbuffer.data(di_geom_11_off + 1025 * ccomps * dcomps);

            auto g_z_x_xx_xyyzzz = cbuffer.data(di_geom_11_off + 1026 * ccomps * dcomps);

            auto g_z_x_xx_xyzzzz = cbuffer.data(di_geom_11_off + 1027 * ccomps * dcomps);

            auto g_z_x_xx_xzzzzz = cbuffer.data(di_geom_11_off + 1028 * ccomps * dcomps);

            auto g_z_x_xx_yyyyyy = cbuffer.data(di_geom_11_off + 1029 * ccomps * dcomps);

            auto g_z_x_xx_yyyyyz = cbuffer.data(di_geom_11_off + 1030 * ccomps * dcomps);

            auto g_z_x_xx_yyyyzz = cbuffer.data(di_geom_11_off + 1031 * ccomps * dcomps);

            auto g_z_x_xx_yyyzzz = cbuffer.data(di_geom_11_off + 1032 * ccomps * dcomps);

            auto g_z_x_xx_yyzzzz = cbuffer.data(di_geom_11_off + 1033 * ccomps * dcomps);

            auto g_z_x_xx_yzzzzz = cbuffer.data(di_geom_11_off + 1034 * ccomps * dcomps);

            auto g_z_x_xx_zzzzzz = cbuffer.data(di_geom_11_off + 1035 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_xxxxxx, g_z_0_x_xxxxxy, g_z_0_x_xxxxxz, g_z_0_x_xxxxyy, g_z_0_x_xxxxyz, g_z_0_x_xxxxzz, g_z_0_x_xxxyyy, g_z_0_x_xxxyyz, g_z_0_x_xxxyzz, g_z_0_x_xxxzzz, g_z_0_x_xxyyyy, g_z_0_x_xxyyyz, g_z_0_x_xxyyzz, g_z_0_x_xxyzzz, g_z_0_x_xxzzzz, g_z_0_x_xyyyyy, g_z_0_x_xyyyyz, g_z_0_x_xyyyzz, g_z_0_x_xyyzzz, g_z_0_x_xyzzzz, g_z_0_x_xzzzzz, g_z_0_x_yyyyyy, g_z_0_x_yyyyyz, g_z_0_x_yyyyzz, g_z_0_x_yyyzzz, g_z_0_x_yyzzzz, g_z_0_x_yzzzzz, g_z_0_x_zzzzzz, g_z_x_x_xxxxxx, g_z_x_x_xxxxxxx, g_z_x_x_xxxxxxy, g_z_x_x_xxxxxxz, g_z_x_x_xxxxxy, g_z_x_x_xxxxxyy, g_z_x_x_xxxxxyz, g_z_x_x_xxxxxz, g_z_x_x_xxxxxzz, g_z_x_x_xxxxyy, g_z_x_x_xxxxyyy, g_z_x_x_xxxxyyz, g_z_x_x_xxxxyz, g_z_x_x_xxxxyzz, g_z_x_x_xxxxzz, g_z_x_x_xxxxzzz, g_z_x_x_xxxyyy, g_z_x_x_xxxyyyy, g_z_x_x_xxxyyyz, g_z_x_x_xxxyyz, g_z_x_x_xxxyyzz, g_z_x_x_xxxyzz, g_z_x_x_xxxyzzz, g_z_x_x_xxxzzz, g_z_x_x_xxxzzzz, g_z_x_x_xxyyyy, g_z_x_x_xxyyyyy, g_z_x_x_xxyyyyz, g_z_x_x_xxyyyz, g_z_x_x_xxyyyzz, g_z_x_x_xxyyzz, g_z_x_x_xxyyzzz, g_z_x_x_xxyzzz, g_z_x_x_xxyzzzz, g_z_x_x_xxzzzz, g_z_x_x_xxzzzzz, g_z_x_x_xyyyyy, g_z_x_x_xyyyyyy, g_z_x_x_xyyyyyz, g_z_x_x_xyyyyz, g_z_x_x_xyyyyzz, g_z_x_x_xyyyzz, g_z_x_x_xyyyzzz, g_z_x_x_xyyzzz, g_z_x_x_xyyzzzz, g_z_x_x_xyzzzz, g_z_x_x_xyzzzzz, g_z_x_x_xzzzzz, g_z_x_x_xzzzzzz, g_z_x_x_yyyyyy, g_z_x_x_yyyyyz, g_z_x_x_yyyyzz, g_z_x_x_yyyzzz, g_z_x_x_yyzzzz, g_z_x_x_yzzzzz, g_z_x_x_zzzzzz, g_z_x_xx_xxxxxx, g_z_x_xx_xxxxxy, g_z_x_xx_xxxxxz, g_z_x_xx_xxxxyy, g_z_x_xx_xxxxyz, g_z_x_xx_xxxxzz, g_z_x_xx_xxxyyy, g_z_x_xx_xxxyyz, g_z_x_xx_xxxyzz, g_z_x_xx_xxxzzz, g_z_x_xx_xxyyyy, g_z_x_xx_xxyyyz, g_z_x_xx_xxyyzz, g_z_x_xx_xxyzzz, g_z_x_xx_xxzzzz, g_z_x_xx_xyyyyy, g_z_x_xx_xyyyyz, g_z_x_xx_xyyyzz, g_z_x_xx_xyyzzz, g_z_x_xx_xyzzzz, g_z_x_xx_xzzzzz, g_z_x_xx_yyyyyy, g_z_x_xx_yyyyyz, g_z_x_xx_yyyyzz, g_z_x_xx_yyyzzz, g_z_x_xx_yyzzzz, g_z_x_xx_yzzzzz, g_z_x_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xx_xxxxxx[k] = g_z_0_x_xxxxxx[k] - g_z_x_x_xxxxxx[k] * ab_x + g_z_x_x_xxxxxxx[k];

                g_z_x_xx_xxxxxy[k] = g_z_0_x_xxxxxy[k] - g_z_x_x_xxxxxy[k] * ab_x + g_z_x_x_xxxxxxy[k];

                g_z_x_xx_xxxxxz[k] = g_z_0_x_xxxxxz[k] - g_z_x_x_xxxxxz[k] * ab_x + g_z_x_x_xxxxxxz[k];

                g_z_x_xx_xxxxyy[k] = g_z_0_x_xxxxyy[k] - g_z_x_x_xxxxyy[k] * ab_x + g_z_x_x_xxxxxyy[k];

                g_z_x_xx_xxxxyz[k] = g_z_0_x_xxxxyz[k] - g_z_x_x_xxxxyz[k] * ab_x + g_z_x_x_xxxxxyz[k];

                g_z_x_xx_xxxxzz[k] = g_z_0_x_xxxxzz[k] - g_z_x_x_xxxxzz[k] * ab_x + g_z_x_x_xxxxxzz[k];

                g_z_x_xx_xxxyyy[k] = g_z_0_x_xxxyyy[k] - g_z_x_x_xxxyyy[k] * ab_x + g_z_x_x_xxxxyyy[k];

                g_z_x_xx_xxxyyz[k] = g_z_0_x_xxxyyz[k] - g_z_x_x_xxxyyz[k] * ab_x + g_z_x_x_xxxxyyz[k];

                g_z_x_xx_xxxyzz[k] = g_z_0_x_xxxyzz[k] - g_z_x_x_xxxyzz[k] * ab_x + g_z_x_x_xxxxyzz[k];

                g_z_x_xx_xxxzzz[k] = g_z_0_x_xxxzzz[k] - g_z_x_x_xxxzzz[k] * ab_x + g_z_x_x_xxxxzzz[k];

                g_z_x_xx_xxyyyy[k] = g_z_0_x_xxyyyy[k] - g_z_x_x_xxyyyy[k] * ab_x + g_z_x_x_xxxyyyy[k];

                g_z_x_xx_xxyyyz[k] = g_z_0_x_xxyyyz[k] - g_z_x_x_xxyyyz[k] * ab_x + g_z_x_x_xxxyyyz[k];

                g_z_x_xx_xxyyzz[k] = g_z_0_x_xxyyzz[k] - g_z_x_x_xxyyzz[k] * ab_x + g_z_x_x_xxxyyzz[k];

                g_z_x_xx_xxyzzz[k] = g_z_0_x_xxyzzz[k] - g_z_x_x_xxyzzz[k] * ab_x + g_z_x_x_xxxyzzz[k];

                g_z_x_xx_xxzzzz[k] = g_z_0_x_xxzzzz[k] - g_z_x_x_xxzzzz[k] * ab_x + g_z_x_x_xxxzzzz[k];

                g_z_x_xx_xyyyyy[k] = g_z_0_x_xyyyyy[k] - g_z_x_x_xyyyyy[k] * ab_x + g_z_x_x_xxyyyyy[k];

                g_z_x_xx_xyyyyz[k] = g_z_0_x_xyyyyz[k] - g_z_x_x_xyyyyz[k] * ab_x + g_z_x_x_xxyyyyz[k];

                g_z_x_xx_xyyyzz[k] = g_z_0_x_xyyyzz[k] - g_z_x_x_xyyyzz[k] * ab_x + g_z_x_x_xxyyyzz[k];

                g_z_x_xx_xyyzzz[k] = g_z_0_x_xyyzzz[k] - g_z_x_x_xyyzzz[k] * ab_x + g_z_x_x_xxyyzzz[k];

                g_z_x_xx_xyzzzz[k] = g_z_0_x_xyzzzz[k] - g_z_x_x_xyzzzz[k] * ab_x + g_z_x_x_xxyzzzz[k];

                g_z_x_xx_xzzzzz[k] = g_z_0_x_xzzzzz[k] - g_z_x_x_xzzzzz[k] * ab_x + g_z_x_x_xxzzzzz[k];

                g_z_x_xx_yyyyyy[k] = g_z_0_x_yyyyyy[k] - g_z_x_x_yyyyyy[k] * ab_x + g_z_x_x_xyyyyyy[k];

                g_z_x_xx_yyyyyz[k] = g_z_0_x_yyyyyz[k] - g_z_x_x_yyyyyz[k] * ab_x + g_z_x_x_xyyyyyz[k];

                g_z_x_xx_yyyyzz[k] = g_z_0_x_yyyyzz[k] - g_z_x_x_yyyyzz[k] * ab_x + g_z_x_x_xyyyyzz[k];

                g_z_x_xx_yyyzzz[k] = g_z_0_x_yyyzzz[k] - g_z_x_x_yyyzzz[k] * ab_x + g_z_x_x_xyyyzzz[k];

                g_z_x_xx_yyzzzz[k] = g_z_0_x_yyzzzz[k] - g_z_x_x_yyzzzz[k] * ab_x + g_z_x_x_xyyzzzz[k];

                g_z_x_xx_yzzzzz[k] = g_z_0_x_yzzzzz[k] - g_z_x_x_yzzzzz[k] * ab_x + g_z_x_x_xyzzzzz[k];

                g_z_x_xx_zzzzzz[k] = g_z_0_x_zzzzzz[k] - g_z_x_x_zzzzzz[k] * ab_x + g_z_x_x_xzzzzzz[k];
            }

            /// Set up 1036-1064 components of targeted buffer : cbuffer.data(

            auto g_z_x_xy_xxxxxx = cbuffer.data(di_geom_11_off + 1036 * ccomps * dcomps);

            auto g_z_x_xy_xxxxxy = cbuffer.data(di_geom_11_off + 1037 * ccomps * dcomps);

            auto g_z_x_xy_xxxxxz = cbuffer.data(di_geom_11_off + 1038 * ccomps * dcomps);

            auto g_z_x_xy_xxxxyy = cbuffer.data(di_geom_11_off + 1039 * ccomps * dcomps);

            auto g_z_x_xy_xxxxyz = cbuffer.data(di_geom_11_off + 1040 * ccomps * dcomps);

            auto g_z_x_xy_xxxxzz = cbuffer.data(di_geom_11_off + 1041 * ccomps * dcomps);

            auto g_z_x_xy_xxxyyy = cbuffer.data(di_geom_11_off + 1042 * ccomps * dcomps);

            auto g_z_x_xy_xxxyyz = cbuffer.data(di_geom_11_off + 1043 * ccomps * dcomps);

            auto g_z_x_xy_xxxyzz = cbuffer.data(di_geom_11_off + 1044 * ccomps * dcomps);

            auto g_z_x_xy_xxxzzz = cbuffer.data(di_geom_11_off + 1045 * ccomps * dcomps);

            auto g_z_x_xy_xxyyyy = cbuffer.data(di_geom_11_off + 1046 * ccomps * dcomps);

            auto g_z_x_xy_xxyyyz = cbuffer.data(di_geom_11_off + 1047 * ccomps * dcomps);

            auto g_z_x_xy_xxyyzz = cbuffer.data(di_geom_11_off + 1048 * ccomps * dcomps);

            auto g_z_x_xy_xxyzzz = cbuffer.data(di_geom_11_off + 1049 * ccomps * dcomps);

            auto g_z_x_xy_xxzzzz = cbuffer.data(di_geom_11_off + 1050 * ccomps * dcomps);

            auto g_z_x_xy_xyyyyy = cbuffer.data(di_geom_11_off + 1051 * ccomps * dcomps);

            auto g_z_x_xy_xyyyyz = cbuffer.data(di_geom_11_off + 1052 * ccomps * dcomps);

            auto g_z_x_xy_xyyyzz = cbuffer.data(di_geom_11_off + 1053 * ccomps * dcomps);

            auto g_z_x_xy_xyyzzz = cbuffer.data(di_geom_11_off + 1054 * ccomps * dcomps);

            auto g_z_x_xy_xyzzzz = cbuffer.data(di_geom_11_off + 1055 * ccomps * dcomps);

            auto g_z_x_xy_xzzzzz = cbuffer.data(di_geom_11_off + 1056 * ccomps * dcomps);

            auto g_z_x_xy_yyyyyy = cbuffer.data(di_geom_11_off + 1057 * ccomps * dcomps);

            auto g_z_x_xy_yyyyyz = cbuffer.data(di_geom_11_off + 1058 * ccomps * dcomps);

            auto g_z_x_xy_yyyyzz = cbuffer.data(di_geom_11_off + 1059 * ccomps * dcomps);

            auto g_z_x_xy_yyyzzz = cbuffer.data(di_geom_11_off + 1060 * ccomps * dcomps);

            auto g_z_x_xy_yyzzzz = cbuffer.data(di_geom_11_off + 1061 * ccomps * dcomps);

            auto g_z_x_xy_yzzzzz = cbuffer.data(di_geom_11_off + 1062 * ccomps * dcomps);

            auto g_z_x_xy_zzzzzz = cbuffer.data(di_geom_11_off + 1063 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_x_xxxxxx, g_z_x_x_xxxxxxy, g_z_x_x_xxxxxy, g_z_x_x_xxxxxyy, g_z_x_x_xxxxxyz, g_z_x_x_xxxxxz, g_z_x_x_xxxxyy, g_z_x_x_xxxxyyy, g_z_x_x_xxxxyyz, g_z_x_x_xxxxyz, g_z_x_x_xxxxyzz, g_z_x_x_xxxxzz, g_z_x_x_xxxyyy, g_z_x_x_xxxyyyy, g_z_x_x_xxxyyyz, g_z_x_x_xxxyyz, g_z_x_x_xxxyyzz, g_z_x_x_xxxyzz, g_z_x_x_xxxyzzz, g_z_x_x_xxxzzz, g_z_x_x_xxyyyy, g_z_x_x_xxyyyyy, g_z_x_x_xxyyyyz, g_z_x_x_xxyyyz, g_z_x_x_xxyyyzz, g_z_x_x_xxyyzz, g_z_x_x_xxyyzzz, g_z_x_x_xxyzzz, g_z_x_x_xxyzzzz, g_z_x_x_xxzzzz, g_z_x_x_xyyyyy, g_z_x_x_xyyyyyy, g_z_x_x_xyyyyyz, g_z_x_x_xyyyyz, g_z_x_x_xyyyyzz, g_z_x_x_xyyyzz, g_z_x_x_xyyyzzz, g_z_x_x_xyyzzz, g_z_x_x_xyyzzzz, g_z_x_x_xyzzzz, g_z_x_x_xyzzzzz, g_z_x_x_xzzzzz, g_z_x_x_yyyyyy, g_z_x_x_yyyyyyy, g_z_x_x_yyyyyyz, g_z_x_x_yyyyyz, g_z_x_x_yyyyyzz, g_z_x_x_yyyyzz, g_z_x_x_yyyyzzz, g_z_x_x_yyyzzz, g_z_x_x_yyyzzzz, g_z_x_x_yyzzzz, g_z_x_x_yyzzzzz, g_z_x_x_yzzzzz, g_z_x_x_yzzzzzz, g_z_x_x_zzzzzz, g_z_x_xy_xxxxxx, g_z_x_xy_xxxxxy, g_z_x_xy_xxxxxz, g_z_x_xy_xxxxyy, g_z_x_xy_xxxxyz, g_z_x_xy_xxxxzz, g_z_x_xy_xxxyyy, g_z_x_xy_xxxyyz, g_z_x_xy_xxxyzz, g_z_x_xy_xxxzzz, g_z_x_xy_xxyyyy, g_z_x_xy_xxyyyz, g_z_x_xy_xxyyzz, g_z_x_xy_xxyzzz, g_z_x_xy_xxzzzz, g_z_x_xy_xyyyyy, g_z_x_xy_xyyyyz, g_z_x_xy_xyyyzz, g_z_x_xy_xyyzzz, g_z_x_xy_xyzzzz, g_z_x_xy_xzzzzz, g_z_x_xy_yyyyyy, g_z_x_xy_yyyyyz, g_z_x_xy_yyyyzz, g_z_x_xy_yyyzzz, g_z_x_xy_yyzzzz, g_z_x_xy_yzzzzz, g_z_x_xy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xy_xxxxxx[k] = -g_z_x_x_xxxxxx[k] * ab_y + g_z_x_x_xxxxxxy[k];

                g_z_x_xy_xxxxxy[k] = -g_z_x_x_xxxxxy[k] * ab_y + g_z_x_x_xxxxxyy[k];

                g_z_x_xy_xxxxxz[k] = -g_z_x_x_xxxxxz[k] * ab_y + g_z_x_x_xxxxxyz[k];

                g_z_x_xy_xxxxyy[k] = -g_z_x_x_xxxxyy[k] * ab_y + g_z_x_x_xxxxyyy[k];

                g_z_x_xy_xxxxyz[k] = -g_z_x_x_xxxxyz[k] * ab_y + g_z_x_x_xxxxyyz[k];

                g_z_x_xy_xxxxzz[k] = -g_z_x_x_xxxxzz[k] * ab_y + g_z_x_x_xxxxyzz[k];

                g_z_x_xy_xxxyyy[k] = -g_z_x_x_xxxyyy[k] * ab_y + g_z_x_x_xxxyyyy[k];

                g_z_x_xy_xxxyyz[k] = -g_z_x_x_xxxyyz[k] * ab_y + g_z_x_x_xxxyyyz[k];

                g_z_x_xy_xxxyzz[k] = -g_z_x_x_xxxyzz[k] * ab_y + g_z_x_x_xxxyyzz[k];

                g_z_x_xy_xxxzzz[k] = -g_z_x_x_xxxzzz[k] * ab_y + g_z_x_x_xxxyzzz[k];

                g_z_x_xy_xxyyyy[k] = -g_z_x_x_xxyyyy[k] * ab_y + g_z_x_x_xxyyyyy[k];

                g_z_x_xy_xxyyyz[k] = -g_z_x_x_xxyyyz[k] * ab_y + g_z_x_x_xxyyyyz[k];

                g_z_x_xy_xxyyzz[k] = -g_z_x_x_xxyyzz[k] * ab_y + g_z_x_x_xxyyyzz[k];

                g_z_x_xy_xxyzzz[k] = -g_z_x_x_xxyzzz[k] * ab_y + g_z_x_x_xxyyzzz[k];

                g_z_x_xy_xxzzzz[k] = -g_z_x_x_xxzzzz[k] * ab_y + g_z_x_x_xxyzzzz[k];

                g_z_x_xy_xyyyyy[k] = -g_z_x_x_xyyyyy[k] * ab_y + g_z_x_x_xyyyyyy[k];

                g_z_x_xy_xyyyyz[k] = -g_z_x_x_xyyyyz[k] * ab_y + g_z_x_x_xyyyyyz[k];

                g_z_x_xy_xyyyzz[k] = -g_z_x_x_xyyyzz[k] * ab_y + g_z_x_x_xyyyyzz[k];

                g_z_x_xy_xyyzzz[k] = -g_z_x_x_xyyzzz[k] * ab_y + g_z_x_x_xyyyzzz[k];

                g_z_x_xy_xyzzzz[k] = -g_z_x_x_xyzzzz[k] * ab_y + g_z_x_x_xyyzzzz[k];

                g_z_x_xy_xzzzzz[k] = -g_z_x_x_xzzzzz[k] * ab_y + g_z_x_x_xyzzzzz[k];

                g_z_x_xy_yyyyyy[k] = -g_z_x_x_yyyyyy[k] * ab_y + g_z_x_x_yyyyyyy[k];

                g_z_x_xy_yyyyyz[k] = -g_z_x_x_yyyyyz[k] * ab_y + g_z_x_x_yyyyyyz[k];

                g_z_x_xy_yyyyzz[k] = -g_z_x_x_yyyyzz[k] * ab_y + g_z_x_x_yyyyyzz[k];

                g_z_x_xy_yyyzzz[k] = -g_z_x_x_yyyzzz[k] * ab_y + g_z_x_x_yyyyzzz[k];

                g_z_x_xy_yyzzzz[k] = -g_z_x_x_yyzzzz[k] * ab_y + g_z_x_x_yyyzzzz[k];

                g_z_x_xy_yzzzzz[k] = -g_z_x_x_yzzzzz[k] * ab_y + g_z_x_x_yyzzzzz[k];

                g_z_x_xy_zzzzzz[k] = -g_z_x_x_zzzzzz[k] * ab_y + g_z_x_x_yzzzzzz[k];
            }

            /// Set up 1064-1092 components of targeted buffer : cbuffer.data(

            auto g_z_x_xz_xxxxxx = cbuffer.data(di_geom_11_off + 1064 * ccomps * dcomps);

            auto g_z_x_xz_xxxxxy = cbuffer.data(di_geom_11_off + 1065 * ccomps * dcomps);

            auto g_z_x_xz_xxxxxz = cbuffer.data(di_geom_11_off + 1066 * ccomps * dcomps);

            auto g_z_x_xz_xxxxyy = cbuffer.data(di_geom_11_off + 1067 * ccomps * dcomps);

            auto g_z_x_xz_xxxxyz = cbuffer.data(di_geom_11_off + 1068 * ccomps * dcomps);

            auto g_z_x_xz_xxxxzz = cbuffer.data(di_geom_11_off + 1069 * ccomps * dcomps);

            auto g_z_x_xz_xxxyyy = cbuffer.data(di_geom_11_off + 1070 * ccomps * dcomps);

            auto g_z_x_xz_xxxyyz = cbuffer.data(di_geom_11_off + 1071 * ccomps * dcomps);

            auto g_z_x_xz_xxxyzz = cbuffer.data(di_geom_11_off + 1072 * ccomps * dcomps);

            auto g_z_x_xz_xxxzzz = cbuffer.data(di_geom_11_off + 1073 * ccomps * dcomps);

            auto g_z_x_xz_xxyyyy = cbuffer.data(di_geom_11_off + 1074 * ccomps * dcomps);

            auto g_z_x_xz_xxyyyz = cbuffer.data(di_geom_11_off + 1075 * ccomps * dcomps);

            auto g_z_x_xz_xxyyzz = cbuffer.data(di_geom_11_off + 1076 * ccomps * dcomps);

            auto g_z_x_xz_xxyzzz = cbuffer.data(di_geom_11_off + 1077 * ccomps * dcomps);

            auto g_z_x_xz_xxzzzz = cbuffer.data(di_geom_11_off + 1078 * ccomps * dcomps);

            auto g_z_x_xz_xyyyyy = cbuffer.data(di_geom_11_off + 1079 * ccomps * dcomps);

            auto g_z_x_xz_xyyyyz = cbuffer.data(di_geom_11_off + 1080 * ccomps * dcomps);

            auto g_z_x_xz_xyyyzz = cbuffer.data(di_geom_11_off + 1081 * ccomps * dcomps);

            auto g_z_x_xz_xyyzzz = cbuffer.data(di_geom_11_off + 1082 * ccomps * dcomps);

            auto g_z_x_xz_xyzzzz = cbuffer.data(di_geom_11_off + 1083 * ccomps * dcomps);

            auto g_z_x_xz_xzzzzz = cbuffer.data(di_geom_11_off + 1084 * ccomps * dcomps);

            auto g_z_x_xz_yyyyyy = cbuffer.data(di_geom_11_off + 1085 * ccomps * dcomps);

            auto g_z_x_xz_yyyyyz = cbuffer.data(di_geom_11_off + 1086 * ccomps * dcomps);

            auto g_z_x_xz_yyyyzz = cbuffer.data(di_geom_11_off + 1087 * ccomps * dcomps);

            auto g_z_x_xz_yyyzzz = cbuffer.data(di_geom_11_off + 1088 * ccomps * dcomps);

            auto g_z_x_xz_yyzzzz = cbuffer.data(di_geom_11_off + 1089 * ccomps * dcomps);

            auto g_z_x_xz_yzzzzz = cbuffer.data(di_geom_11_off + 1090 * ccomps * dcomps);

            auto g_z_x_xz_zzzzzz = cbuffer.data(di_geom_11_off + 1091 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxxxx, g_z_0_z_xxxxxy, g_z_0_z_xxxxxz, g_z_0_z_xxxxyy, g_z_0_z_xxxxyz, g_z_0_z_xxxxzz, g_z_0_z_xxxyyy, g_z_0_z_xxxyyz, g_z_0_z_xxxyzz, g_z_0_z_xxxzzz, g_z_0_z_xxyyyy, g_z_0_z_xxyyyz, g_z_0_z_xxyyzz, g_z_0_z_xxyzzz, g_z_0_z_xxzzzz, g_z_0_z_xyyyyy, g_z_0_z_xyyyyz, g_z_0_z_xyyyzz, g_z_0_z_xyyzzz, g_z_0_z_xyzzzz, g_z_0_z_xzzzzz, g_z_0_z_yyyyyy, g_z_0_z_yyyyyz, g_z_0_z_yyyyzz, g_z_0_z_yyyzzz, g_z_0_z_yyzzzz, g_z_0_z_yzzzzz, g_z_0_z_zzzzzz, g_z_x_xz_xxxxxx, g_z_x_xz_xxxxxy, g_z_x_xz_xxxxxz, g_z_x_xz_xxxxyy, g_z_x_xz_xxxxyz, g_z_x_xz_xxxxzz, g_z_x_xz_xxxyyy, g_z_x_xz_xxxyyz, g_z_x_xz_xxxyzz, g_z_x_xz_xxxzzz, g_z_x_xz_xxyyyy, g_z_x_xz_xxyyyz, g_z_x_xz_xxyyzz, g_z_x_xz_xxyzzz, g_z_x_xz_xxzzzz, g_z_x_xz_xyyyyy, g_z_x_xz_xyyyyz, g_z_x_xz_xyyyzz, g_z_x_xz_xyyzzz, g_z_x_xz_xyzzzz, g_z_x_xz_xzzzzz, g_z_x_xz_yyyyyy, g_z_x_xz_yyyyyz, g_z_x_xz_yyyyzz, g_z_x_xz_yyyzzz, g_z_x_xz_yyzzzz, g_z_x_xz_yzzzzz, g_z_x_xz_zzzzzz, g_z_x_z_xxxxxx, g_z_x_z_xxxxxxx, g_z_x_z_xxxxxxy, g_z_x_z_xxxxxxz, g_z_x_z_xxxxxy, g_z_x_z_xxxxxyy, g_z_x_z_xxxxxyz, g_z_x_z_xxxxxz, g_z_x_z_xxxxxzz, g_z_x_z_xxxxyy, g_z_x_z_xxxxyyy, g_z_x_z_xxxxyyz, g_z_x_z_xxxxyz, g_z_x_z_xxxxyzz, g_z_x_z_xxxxzz, g_z_x_z_xxxxzzz, g_z_x_z_xxxyyy, g_z_x_z_xxxyyyy, g_z_x_z_xxxyyyz, g_z_x_z_xxxyyz, g_z_x_z_xxxyyzz, g_z_x_z_xxxyzz, g_z_x_z_xxxyzzz, g_z_x_z_xxxzzz, g_z_x_z_xxxzzzz, g_z_x_z_xxyyyy, g_z_x_z_xxyyyyy, g_z_x_z_xxyyyyz, g_z_x_z_xxyyyz, g_z_x_z_xxyyyzz, g_z_x_z_xxyyzz, g_z_x_z_xxyyzzz, g_z_x_z_xxyzzz, g_z_x_z_xxyzzzz, g_z_x_z_xxzzzz, g_z_x_z_xxzzzzz, g_z_x_z_xyyyyy, g_z_x_z_xyyyyyy, g_z_x_z_xyyyyyz, g_z_x_z_xyyyyz, g_z_x_z_xyyyyzz, g_z_x_z_xyyyzz, g_z_x_z_xyyyzzz, g_z_x_z_xyyzzz, g_z_x_z_xyyzzzz, g_z_x_z_xyzzzz, g_z_x_z_xyzzzzz, g_z_x_z_xzzzzz, g_z_x_z_xzzzzzz, g_z_x_z_yyyyyy, g_z_x_z_yyyyyz, g_z_x_z_yyyyzz, g_z_x_z_yyyzzz, g_z_x_z_yyzzzz, g_z_x_z_yzzzzz, g_z_x_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xz_xxxxxx[k] = g_z_0_z_xxxxxx[k] - g_z_x_z_xxxxxx[k] * ab_x + g_z_x_z_xxxxxxx[k];

                g_z_x_xz_xxxxxy[k] = g_z_0_z_xxxxxy[k] - g_z_x_z_xxxxxy[k] * ab_x + g_z_x_z_xxxxxxy[k];

                g_z_x_xz_xxxxxz[k] = g_z_0_z_xxxxxz[k] - g_z_x_z_xxxxxz[k] * ab_x + g_z_x_z_xxxxxxz[k];

                g_z_x_xz_xxxxyy[k] = g_z_0_z_xxxxyy[k] - g_z_x_z_xxxxyy[k] * ab_x + g_z_x_z_xxxxxyy[k];

                g_z_x_xz_xxxxyz[k] = g_z_0_z_xxxxyz[k] - g_z_x_z_xxxxyz[k] * ab_x + g_z_x_z_xxxxxyz[k];

                g_z_x_xz_xxxxzz[k] = g_z_0_z_xxxxzz[k] - g_z_x_z_xxxxzz[k] * ab_x + g_z_x_z_xxxxxzz[k];

                g_z_x_xz_xxxyyy[k] = g_z_0_z_xxxyyy[k] - g_z_x_z_xxxyyy[k] * ab_x + g_z_x_z_xxxxyyy[k];

                g_z_x_xz_xxxyyz[k] = g_z_0_z_xxxyyz[k] - g_z_x_z_xxxyyz[k] * ab_x + g_z_x_z_xxxxyyz[k];

                g_z_x_xz_xxxyzz[k] = g_z_0_z_xxxyzz[k] - g_z_x_z_xxxyzz[k] * ab_x + g_z_x_z_xxxxyzz[k];

                g_z_x_xz_xxxzzz[k] = g_z_0_z_xxxzzz[k] - g_z_x_z_xxxzzz[k] * ab_x + g_z_x_z_xxxxzzz[k];

                g_z_x_xz_xxyyyy[k] = g_z_0_z_xxyyyy[k] - g_z_x_z_xxyyyy[k] * ab_x + g_z_x_z_xxxyyyy[k];

                g_z_x_xz_xxyyyz[k] = g_z_0_z_xxyyyz[k] - g_z_x_z_xxyyyz[k] * ab_x + g_z_x_z_xxxyyyz[k];

                g_z_x_xz_xxyyzz[k] = g_z_0_z_xxyyzz[k] - g_z_x_z_xxyyzz[k] * ab_x + g_z_x_z_xxxyyzz[k];

                g_z_x_xz_xxyzzz[k] = g_z_0_z_xxyzzz[k] - g_z_x_z_xxyzzz[k] * ab_x + g_z_x_z_xxxyzzz[k];

                g_z_x_xz_xxzzzz[k] = g_z_0_z_xxzzzz[k] - g_z_x_z_xxzzzz[k] * ab_x + g_z_x_z_xxxzzzz[k];

                g_z_x_xz_xyyyyy[k] = g_z_0_z_xyyyyy[k] - g_z_x_z_xyyyyy[k] * ab_x + g_z_x_z_xxyyyyy[k];

                g_z_x_xz_xyyyyz[k] = g_z_0_z_xyyyyz[k] - g_z_x_z_xyyyyz[k] * ab_x + g_z_x_z_xxyyyyz[k];

                g_z_x_xz_xyyyzz[k] = g_z_0_z_xyyyzz[k] - g_z_x_z_xyyyzz[k] * ab_x + g_z_x_z_xxyyyzz[k];

                g_z_x_xz_xyyzzz[k] = g_z_0_z_xyyzzz[k] - g_z_x_z_xyyzzz[k] * ab_x + g_z_x_z_xxyyzzz[k];

                g_z_x_xz_xyzzzz[k] = g_z_0_z_xyzzzz[k] - g_z_x_z_xyzzzz[k] * ab_x + g_z_x_z_xxyzzzz[k];

                g_z_x_xz_xzzzzz[k] = g_z_0_z_xzzzzz[k] - g_z_x_z_xzzzzz[k] * ab_x + g_z_x_z_xxzzzzz[k];

                g_z_x_xz_yyyyyy[k] = g_z_0_z_yyyyyy[k] - g_z_x_z_yyyyyy[k] * ab_x + g_z_x_z_xyyyyyy[k];

                g_z_x_xz_yyyyyz[k] = g_z_0_z_yyyyyz[k] - g_z_x_z_yyyyyz[k] * ab_x + g_z_x_z_xyyyyyz[k];

                g_z_x_xz_yyyyzz[k] = g_z_0_z_yyyyzz[k] - g_z_x_z_yyyyzz[k] * ab_x + g_z_x_z_xyyyyzz[k];

                g_z_x_xz_yyyzzz[k] = g_z_0_z_yyyzzz[k] - g_z_x_z_yyyzzz[k] * ab_x + g_z_x_z_xyyyzzz[k];

                g_z_x_xz_yyzzzz[k] = g_z_0_z_yyzzzz[k] - g_z_x_z_yyzzzz[k] * ab_x + g_z_x_z_xyyzzzz[k];

                g_z_x_xz_yzzzzz[k] = g_z_0_z_yzzzzz[k] - g_z_x_z_yzzzzz[k] * ab_x + g_z_x_z_xyzzzzz[k];

                g_z_x_xz_zzzzzz[k] = g_z_0_z_zzzzzz[k] - g_z_x_z_zzzzzz[k] * ab_x + g_z_x_z_xzzzzzz[k];
            }

            /// Set up 1092-1120 components of targeted buffer : cbuffer.data(

            auto g_z_x_yy_xxxxxx = cbuffer.data(di_geom_11_off + 1092 * ccomps * dcomps);

            auto g_z_x_yy_xxxxxy = cbuffer.data(di_geom_11_off + 1093 * ccomps * dcomps);

            auto g_z_x_yy_xxxxxz = cbuffer.data(di_geom_11_off + 1094 * ccomps * dcomps);

            auto g_z_x_yy_xxxxyy = cbuffer.data(di_geom_11_off + 1095 * ccomps * dcomps);

            auto g_z_x_yy_xxxxyz = cbuffer.data(di_geom_11_off + 1096 * ccomps * dcomps);

            auto g_z_x_yy_xxxxzz = cbuffer.data(di_geom_11_off + 1097 * ccomps * dcomps);

            auto g_z_x_yy_xxxyyy = cbuffer.data(di_geom_11_off + 1098 * ccomps * dcomps);

            auto g_z_x_yy_xxxyyz = cbuffer.data(di_geom_11_off + 1099 * ccomps * dcomps);

            auto g_z_x_yy_xxxyzz = cbuffer.data(di_geom_11_off + 1100 * ccomps * dcomps);

            auto g_z_x_yy_xxxzzz = cbuffer.data(di_geom_11_off + 1101 * ccomps * dcomps);

            auto g_z_x_yy_xxyyyy = cbuffer.data(di_geom_11_off + 1102 * ccomps * dcomps);

            auto g_z_x_yy_xxyyyz = cbuffer.data(di_geom_11_off + 1103 * ccomps * dcomps);

            auto g_z_x_yy_xxyyzz = cbuffer.data(di_geom_11_off + 1104 * ccomps * dcomps);

            auto g_z_x_yy_xxyzzz = cbuffer.data(di_geom_11_off + 1105 * ccomps * dcomps);

            auto g_z_x_yy_xxzzzz = cbuffer.data(di_geom_11_off + 1106 * ccomps * dcomps);

            auto g_z_x_yy_xyyyyy = cbuffer.data(di_geom_11_off + 1107 * ccomps * dcomps);

            auto g_z_x_yy_xyyyyz = cbuffer.data(di_geom_11_off + 1108 * ccomps * dcomps);

            auto g_z_x_yy_xyyyzz = cbuffer.data(di_geom_11_off + 1109 * ccomps * dcomps);

            auto g_z_x_yy_xyyzzz = cbuffer.data(di_geom_11_off + 1110 * ccomps * dcomps);

            auto g_z_x_yy_xyzzzz = cbuffer.data(di_geom_11_off + 1111 * ccomps * dcomps);

            auto g_z_x_yy_xzzzzz = cbuffer.data(di_geom_11_off + 1112 * ccomps * dcomps);

            auto g_z_x_yy_yyyyyy = cbuffer.data(di_geom_11_off + 1113 * ccomps * dcomps);

            auto g_z_x_yy_yyyyyz = cbuffer.data(di_geom_11_off + 1114 * ccomps * dcomps);

            auto g_z_x_yy_yyyyzz = cbuffer.data(di_geom_11_off + 1115 * ccomps * dcomps);

            auto g_z_x_yy_yyyzzz = cbuffer.data(di_geom_11_off + 1116 * ccomps * dcomps);

            auto g_z_x_yy_yyzzzz = cbuffer.data(di_geom_11_off + 1117 * ccomps * dcomps);

            auto g_z_x_yy_yzzzzz = cbuffer.data(di_geom_11_off + 1118 * ccomps * dcomps);

            auto g_z_x_yy_zzzzzz = cbuffer.data(di_geom_11_off + 1119 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_y_xxxxxx, g_z_x_y_xxxxxxy, g_z_x_y_xxxxxy, g_z_x_y_xxxxxyy, g_z_x_y_xxxxxyz, g_z_x_y_xxxxxz, g_z_x_y_xxxxyy, g_z_x_y_xxxxyyy, g_z_x_y_xxxxyyz, g_z_x_y_xxxxyz, g_z_x_y_xxxxyzz, g_z_x_y_xxxxzz, g_z_x_y_xxxyyy, g_z_x_y_xxxyyyy, g_z_x_y_xxxyyyz, g_z_x_y_xxxyyz, g_z_x_y_xxxyyzz, g_z_x_y_xxxyzz, g_z_x_y_xxxyzzz, g_z_x_y_xxxzzz, g_z_x_y_xxyyyy, g_z_x_y_xxyyyyy, g_z_x_y_xxyyyyz, g_z_x_y_xxyyyz, g_z_x_y_xxyyyzz, g_z_x_y_xxyyzz, g_z_x_y_xxyyzzz, g_z_x_y_xxyzzz, g_z_x_y_xxyzzzz, g_z_x_y_xxzzzz, g_z_x_y_xyyyyy, g_z_x_y_xyyyyyy, g_z_x_y_xyyyyyz, g_z_x_y_xyyyyz, g_z_x_y_xyyyyzz, g_z_x_y_xyyyzz, g_z_x_y_xyyyzzz, g_z_x_y_xyyzzz, g_z_x_y_xyyzzzz, g_z_x_y_xyzzzz, g_z_x_y_xyzzzzz, g_z_x_y_xzzzzz, g_z_x_y_yyyyyy, g_z_x_y_yyyyyyy, g_z_x_y_yyyyyyz, g_z_x_y_yyyyyz, g_z_x_y_yyyyyzz, g_z_x_y_yyyyzz, g_z_x_y_yyyyzzz, g_z_x_y_yyyzzz, g_z_x_y_yyyzzzz, g_z_x_y_yyzzzz, g_z_x_y_yyzzzzz, g_z_x_y_yzzzzz, g_z_x_y_yzzzzzz, g_z_x_y_zzzzzz, g_z_x_yy_xxxxxx, g_z_x_yy_xxxxxy, g_z_x_yy_xxxxxz, g_z_x_yy_xxxxyy, g_z_x_yy_xxxxyz, g_z_x_yy_xxxxzz, g_z_x_yy_xxxyyy, g_z_x_yy_xxxyyz, g_z_x_yy_xxxyzz, g_z_x_yy_xxxzzz, g_z_x_yy_xxyyyy, g_z_x_yy_xxyyyz, g_z_x_yy_xxyyzz, g_z_x_yy_xxyzzz, g_z_x_yy_xxzzzz, g_z_x_yy_xyyyyy, g_z_x_yy_xyyyyz, g_z_x_yy_xyyyzz, g_z_x_yy_xyyzzz, g_z_x_yy_xyzzzz, g_z_x_yy_xzzzzz, g_z_x_yy_yyyyyy, g_z_x_yy_yyyyyz, g_z_x_yy_yyyyzz, g_z_x_yy_yyyzzz, g_z_x_yy_yyzzzz, g_z_x_yy_yzzzzz, g_z_x_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yy_xxxxxx[k] = -g_z_x_y_xxxxxx[k] * ab_y + g_z_x_y_xxxxxxy[k];

                g_z_x_yy_xxxxxy[k] = -g_z_x_y_xxxxxy[k] * ab_y + g_z_x_y_xxxxxyy[k];

                g_z_x_yy_xxxxxz[k] = -g_z_x_y_xxxxxz[k] * ab_y + g_z_x_y_xxxxxyz[k];

                g_z_x_yy_xxxxyy[k] = -g_z_x_y_xxxxyy[k] * ab_y + g_z_x_y_xxxxyyy[k];

                g_z_x_yy_xxxxyz[k] = -g_z_x_y_xxxxyz[k] * ab_y + g_z_x_y_xxxxyyz[k];

                g_z_x_yy_xxxxzz[k] = -g_z_x_y_xxxxzz[k] * ab_y + g_z_x_y_xxxxyzz[k];

                g_z_x_yy_xxxyyy[k] = -g_z_x_y_xxxyyy[k] * ab_y + g_z_x_y_xxxyyyy[k];

                g_z_x_yy_xxxyyz[k] = -g_z_x_y_xxxyyz[k] * ab_y + g_z_x_y_xxxyyyz[k];

                g_z_x_yy_xxxyzz[k] = -g_z_x_y_xxxyzz[k] * ab_y + g_z_x_y_xxxyyzz[k];

                g_z_x_yy_xxxzzz[k] = -g_z_x_y_xxxzzz[k] * ab_y + g_z_x_y_xxxyzzz[k];

                g_z_x_yy_xxyyyy[k] = -g_z_x_y_xxyyyy[k] * ab_y + g_z_x_y_xxyyyyy[k];

                g_z_x_yy_xxyyyz[k] = -g_z_x_y_xxyyyz[k] * ab_y + g_z_x_y_xxyyyyz[k];

                g_z_x_yy_xxyyzz[k] = -g_z_x_y_xxyyzz[k] * ab_y + g_z_x_y_xxyyyzz[k];

                g_z_x_yy_xxyzzz[k] = -g_z_x_y_xxyzzz[k] * ab_y + g_z_x_y_xxyyzzz[k];

                g_z_x_yy_xxzzzz[k] = -g_z_x_y_xxzzzz[k] * ab_y + g_z_x_y_xxyzzzz[k];

                g_z_x_yy_xyyyyy[k] = -g_z_x_y_xyyyyy[k] * ab_y + g_z_x_y_xyyyyyy[k];

                g_z_x_yy_xyyyyz[k] = -g_z_x_y_xyyyyz[k] * ab_y + g_z_x_y_xyyyyyz[k];

                g_z_x_yy_xyyyzz[k] = -g_z_x_y_xyyyzz[k] * ab_y + g_z_x_y_xyyyyzz[k];

                g_z_x_yy_xyyzzz[k] = -g_z_x_y_xyyzzz[k] * ab_y + g_z_x_y_xyyyzzz[k];

                g_z_x_yy_xyzzzz[k] = -g_z_x_y_xyzzzz[k] * ab_y + g_z_x_y_xyyzzzz[k];

                g_z_x_yy_xzzzzz[k] = -g_z_x_y_xzzzzz[k] * ab_y + g_z_x_y_xyzzzzz[k];

                g_z_x_yy_yyyyyy[k] = -g_z_x_y_yyyyyy[k] * ab_y + g_z_x_y_yyyyyyy[k];

                g_z_x_yy_yyyyyz[k] = -g_z_x_y_yyyyyz[k] * ab_y + g_z_x_y_yyyyyyz[k];

                g_z_x_yy_yyyyzz[k] = -g_z_x_y_yyyyzz[k] * ab_y + g_z_x_y_yyyyyzz[k];

                g_z_x_yy_yyyzzz[k] = -g_z_x_y_yyyzzz[k] * ab_y + g_z_x_y_yyyyzzz[k];

                g_z_x_yy_yyzzzz[k] = -g_z_x_y_yyzzzz[k] * ab_y + g_z_x_y_yyyzzzz[k];

                g_z_x_yy_yzzzzz[k] = -g_z_x_y_yzzzzz[k] * ab_y + g_z_x_y_yyzzzzz[k];

                g_z_x_yy_zzzzzz[k] = -g_z_x_y_zzzzzz[k] * ab_y + g_z_x_y_yzzzzzz[k];
            }

            /// Set up 1120-1148 components of targeted buffer : cbuffer.data(

            auto g_z_x_yz_xxxxxx = cbuffer.data(di_geom_11_off + 1120 * ccomps * dcomps);

            auto g_z_x_yz_xxxxxy = cbuffer.data(di_geom_11_off + 1121 * ccomps * dcomps);

            auto g_z_x_yz_xxxxxz = cbuffer.data(di_geom_11_off + 1122 * ccomps * dcomps);

            auto g_z_x_yz_xxxxyy = cbuffer.data(di_geom_11_off + 1123 * ccomps * dcomps);

            auto g_z_x_yz_xxxxyz = cbuffer.data(di_geom_11_off + 1124 * ccomps * dcomps);

            auto g_z_x_yz_xxxxzz = cbuffer.data(di_geom_11_off + 1125 * ccomps * dcomps);

            auto g_z_x_yz_xxxyyy = cbuffer.data(di_geom_11_off + 1126 * ccomps * dcomps);

            auto g_z_x_yz_xxxyyz = cbuffer.data(di_geom_11_off + 1127 * ccomps * dcomps);

            auto g_z_x_yz_xxxyzz = cbuffer.data(di_geom_11_off + 1128 * ccomps * dcomps);

            auto g_z_x_yz_xxxzzz = cbuffer.data(di_geom_11_off + 1129 * ccomps * dcomps);

            auto g_z_x_yz_xxyyyy = cbuffer.data(di_geom_11_off + 1130 * ccomps * dcomps);

            auto g_z_x_yz_xxyyyz = cbuffer.data(di_geom_11_off + 1131 * ccomps * dcomps);

            auto g_z_x_yz_xxyyzz = cbuffer.data(di_geom_11_off + 1132 * ccomps * dcomps);

            auto g_z_x_yz_xxyzzz = cbuffer.data(di_geom_11_off + 1133 * ccomps * dcomps);

            auto g_z_x_yz_xxzzzz = cbuffer.data(di_geom_11_off + 1134 * ccomps * dcomps);

            auto g_z_x_yz_xyyyyy = cbuffer.data(di_geom_11_off + 1135 * ccomps * dcomps);

            auto g_z_x_yz_xyyyyz = cbuffer.data(di_geom_11_off + 1136 * ccomps * dcomps);

            auto g_z_x_yz_xyyyzz = cbuffer.data(di_geom_11_off + 1137 * ccomps * dcomps);

            auto g_z_x_yz_xyyzzz = cbuffer.data(di_geom_11_off + 1138 * ccomps * dcomps);

            auto g_z_x_yz_xyzzzz = cbuffer.data(di_geom_11_off + 1139 * ccomps * dcomps);

            auto g_z_x_yz_xzzzzz = cbuffer.data(di_geom_11_off + 1140 * ccomps * dcomps);

            auto g_z_x_yz_yyyyyy = cbuffer.data(di_geom_11_off + 1141 * ccomps * dcomps);

            auto g_z_x_yz_yyyyyz = cbuffer.data(di_geom_11_off + 1142 * ccomps * dcomps);

            auto g_z_x_yz_yyyyzz = cbuffer.data(di_geom_11_off + 1143 * ccomps * dcomps);

            auto g_z_x_yz_yyyzzz = cbuffer.data(di_geom_11_off + 1144 * ccomps * dcomps);

            auto g_z_x_yz_yyzzzz = cbuffer.data(di_geom_11_off + 1145 * ccomps * dcomps);

            auto g_z_x_yz_yzzzzz = cbuffer.data(di_geom_11_off + 1146 * ccomps * dcomps);

            auto g_z_x_yz_zzzzzz = cbuffer.data(di_geom_11_off + 1147 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yz_xxxxxx, g_z_x_yz_xxxxxy, g_z_x_yz_xxxxxz, g_z_x_yz_xxxxyy, g_z_x_yz_xxxxyz, g_z_x_yz_xxxxzz, g_z_x_yz_xxxyyy, g_z_x_yz_xxxyyz, g_z_x_yz_xxxyzz, g_z_x_yz_xxxzzz, g_z_x_yz_xxyyyy, g_z_x_yz_xxyyyz, g_z_x_yz_xxyyzz, g_z_x_yz_xxyzzz, g_z_x_yz_xxzzzz, g_z_x_yz_xyyyyy, g_z_x_yz_xyyyyz, g_z_x_yz_xyyyzz, g_z_x_yz_xyyzzz, g_z_x_yz_xyzzzz, g_z_x_yz_xzzzzz, g_z_x_yz_yyyyyy, g_z_x_yz_yyyyyz, g_z_x_yz_yyyyzz, g_z_x_yz_yyyzzz, g_z_x_yz_yyzzzz, g_z_x_yz_yzzzzz, g_z_x_yz_zzzzzz, g_z_x_z_xxxxxx, g_z_x_z_xxxxxxy, g_z_x_z_xxxxxy, g_z_x_z_xxxxxyy, g_z_x_z_xxxxxyz, g_z_x_z_xxxxxz, g_z_x_z_xxxxyy, g_z_x_z_xxxxyyy, g_z_x_z_xxxxyyz, g_z_x_z_xxxxyz, g_z_x_z_xxxxyzz, g_z_x_z_xxxxzz, g_z_x_z_xxxyyy, g_z_x_z_xxxyyyy, g_z_x_z_xxxyyyz, g_z_x_z_xxxyyz, g_z_x_z_xxxyyzz, g_z_x_z_xxxyzz, g_z_x_z_xxxyzzz, g_z_x_z_xxxzzz, g_z_x_z_xxyyyy, g_z_x_z_xxyyyyy, g_z_x_z_xxyyyyz, g_z_x_z_xxyyyz, g_z_x_z_xxyyyzz, g_z_x_z_xxyyzz, g_z_x_z_xxyyzzz, g_z_x_z_xxyzzz, g_z_x_z_xxyzzzz, g_z_x_z_xxzzzz, g_z_x_z_xyyyyy, g_z_x_z_xyyyyyy, g_z_x_z_xyyyyyz, g_z_x_z_xyyyyz, g_z_x_z_xyyyyzz, g_z_x_z_xyyyzz, g_z_x_z_xyyyzzz, g_z_x_z_xyyzzz, g_z_x_z_xyyzzzz, g_z_x_z_xyzzzz, g_z_x_z_xyzzzzz, g_z_x_z_xzzzzz, g_z_x_z_yyyyyy, g_z_x_z_yyyyyyy, g_z_x_z_yyyyyyz, g_z_x_z_yyyyyz, g_z_x_z_yyyyyzz, g_z_x_z_yyyyzz, g_z_x_z_yyyyzzz, g_z_x_z_yyyzzz, g_z_x_z_yyyzzzz, g_z_x_z_yyzzzz, g_z_x_z_yyzzzzz, g_z_x_z_yzzzzz, g_z_x_z_yzzzzzz, g_z_x_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yz_xxxxxx[k] = -g_z_x_z_xxxxxx[k] * ab_y + g_z_x_z_xxxxxxy[k];

                g_z_x_yz_xxxxxy[k] = -g_z_x_z_xxxxxy[k] * ab_y + g_z_x_z_xxxxxyy[k];

                g_z_x_yz_xxxxxz[k] = -g_z_x_z_xxxxxz[k] * ab_y + g_z_x_z_xxxxxyz[k];

                g_z_x_yz_xxxxyy[k] = -g_z_x_z_xxxxyy[k] * ab_y + g_z_x_z_xxxxyyy[k];

                g_z_x_yz_xxxxyz[k] = -g_z_x_z_xxxxyz[k] * ab_y + g_z_x_z_xxxxyyz[k];

                g_z_x_yz_xxxxzz[k] = -g_z_x_z_xxxxzz[k] * ab_y + g_z_x_z_xxxxyzz[k];

                g_z_x_yz_xxxyyy[k] = -g_z_x_z_xxxyyy[k] * ab_y + g_z_x_z_xxxyyyy[k];

                g_z_x_yz_xxxyyz[k] = -g_z_x_z_xxxyyz[k] * ab_y + g_z_x_z_xxxyyyz[k];

                g_z_x_yz_xxxyzz[k] = -g_z_x_z_xxxyzz[k] * ab_y + g_z_x_z_xxxyyzz[k];

                g_z_x_yz_xxxzzz[k] = -g_z_x_z_xxxzzz[k] * ab_y + g_z_x_z_xxxyzzz[k];

                g_z_x_yz_xxyyyy[k] = -g_z_x_z_xxyyyy[k] * ab_y + g_z_x_z_xxyyyyy[k];

                g_z_x_yz_xxyyyz[k] = -g_z_x_z_xxyyyz[k] * ab_y + g_z_x_z_xxyyyyz[k];

                g_z_x_yz_xxyyzz[k] = -g_z_x_z_xxyyzz[k] * ab_y + g_z_x_z_xxyyyzz[k];

                g_z_x_yz_xxyzzz[k] = -g_z_x_z_xxyzzz[k] * ab_y + g_z_x_z_xxyyzzz[k];

                g_z_x_yz_xxzzzz[k] = -g_z_x_z_xxzzzz[k] * ab_y + g_z_x_z_xxyzzzz[k];

                g_z_x_yz_xyyyyy[k] = -g_z_x_z_xyyyyy[k] * ab_y + g_z_x_z_xyyyyyy[k];

                g_z_x_yz_xyyyyz[k] = -g_z_x_z_xyyyyz[k] * ab_y + g_z_x_z_xyyyyyz[k];

                g_z_x_yz_xyyyzz[k] = -g_z_x_z_xyyyzz[k] * ab_y + g_z_x_z_xyyyyzz[k];

                g_z_x_yz_xyyzzz[k] = -g_z_x_z_xyyzzz[k] * ab_y + g_z_x_z_xyyyzzz[k];

                g_z_x_yz_xyzzzz[k] = -g_z_x_z_xyzzzz[k] * ab_y + g_z_x_z_xyyzzzz[k];

                g_z_x_yz_xzzzzz[k] = -g_z_x_z_xzzzzz[k] * ab_y + g_z_x_z_xyzzzzz[k];

                g_z_x_yz_yyyyyy[k] = -g_z_x_z_yyyyyy[k] * ab_y + g_z_x_z_yyyyyyy[k];

                g_z_x_yz_yyyyyz[k] = -g_z_x_z_yyyyyz[k] * ab_y + g_z_x_z_yyyyyyz[k];

                g_z_x_yz_yyyyzz[k] = -g_z_x_z_yyyyzz[k] * ab_y + g_z_x_z_yyyyyzz[k];

                g_z_x_yz_yyyzzz[k] = -g_z_x_z_yyyzzz[k] * ab_y + g_z_x_z_yyyyzzz[k];

                g_z_x_yz_yyzzzz[k] = -g_z_x_z_yyzzzz[k] * ab_y + g_z_x_z_yyyzzzz[k];

                g_z_x_yz_yzzzzz[k] = -g_z_x_z_yzzzzz[k] * ab_y + g_z_x_z_yyzzzzz[k];

                g_z_x_yz_zzzzzz[k] = -g_z_x_z_zzzzzz[k] * ab_y + g_z_x_z_yzzzzzz[k];
            }

            /// Set up 1148-1176 components of targeted buffer : cbuffer.data(

            auto g_z_x_zz_xxxxxx = cbuffer.data(di_geom_11_off + 1148 * ccomps * dcomps);

            auto g_z_x_zz_xxxxxy = cbuffer.data(di_geom_11_off + 1149 * ccomps * dcomps);

            auto g_z_x_zz_xxxxxz = cbuffer.data(di_geom_11_off + 1150 * ccomps * dcomps);

            auto g_z_x_zz_xxxxyy = cbuffer.data(di_geom_11_off + 1151 * ccomps * dcomps);

            auto g_z_x_zz_xxxxyz = cbuffer.data(di_geom_11_off + 1152 * ccomps * dcomps);

            auto g_z_x_zz_xxxxzz = cbuffer.data(di_geom_11_off + 1153 * ccomps * dcomps);

            auto g_z_x_zz_xxxyyy = cbuffer.data(di_geom_11_off + 1154 * ccomps * dcomps);

            auto g_z_x_zz_xxxyyz = cbuffer.data(di_geom_11_off + 1155 * ccomps * dcomps);

            auto g_z_x_zz_xxxyzz = cbuffer.data(di_geom_11_off + 1156 * ccomps * dcomps);

            auto g_z_x_zz_xxxzzz = cbuffer.data(di_geom_11_off + 1157 * ccomps * dcomps);

            auto g_z_x_zz_xxyyyy = cbuffer.data(di_geom_11_off + 1158 * ccomps * dcomps);

            auto g_z_x_zz_xxyyyz = cbuffer.data(di_geom_11_off + 1159 * ccomps * dcomps);

            auto g_z_x_zz_xxyyzz = cbuffer.data(di_geom_11_off + 1160 * ccomps * dcomps);

            auto g_z_x_zz_xxyzzz = cbuffer.data(di_geom_11_off + 1161 * ccomps * dcomps);

            auto g_z_x_zz_xxzzzz = cbuffer.data(di_geom_11_off + 1162 * ccomps * dcomps);

            auto g_z_x_zz_xyyyyy = cbuffer.data(di_geom_11_off + 1163 * ccomps * dcomps);

            auto g_z_x_zz_xyyyyz = cbuffer.data(di_geom_11_off + 1164 * ccomps * dcomps);

            auto g_z_x_zz_xyyyzz = cbuffer.data(di_geom_11_off + 1165 * ccomps * dcomps);

            auto g_z_x_zz_xyyzzz = cbuffer.data(di_geom_11_off + 1166 * ccomps * dcomps);

            auto g_z_x_zz_xyzzzz = cbuffer.data(di_geom_11_off + 1167 * ccomps * dcomps);

            auto g_z_x_zz_xzzzzz = cbuffer.data(di_geom_11_off + 1168 * ccomps * dcomps);

            auto g_z_x_zz_yyyyyy = cbuffer.data(di_geom_11_off + 1169 * ccomps * dcomps);

            auto g_z_x_zz_yyyyyz = cbuffer.data(di_geom_11_off + 1170 * ccomps * dcomps);

            auto g_z_x_zz_yyyyzz = cbuffer.data(di_geom_11_off + 1171 * ccomps * dcomps);

            auto g_z_x_zz_yyyzzz = cbuffer.data(di_geom_11_off + 1172 * ccomps * dcomps);

            auto g_z_x_zz_yyzzzz = cbuffer.data(di_geom_11_off + 1173 * ccomps * dcomps);

            auto g_z_x_zz_yzzzzz = cbuffer.data(di_geom_11_off + 1174 * ccomps * dcomps);

            auto g_z_x_zz_zzzzzz = cbuffer.data(di_geom_11_off + 1175 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_z_xxxxxx, g_0_x_z_xxxxxy, g_0_x_z_xxxxxz, g_0_x_z_xxxxyy, g_0_x_z_xxxxyz, g_0_x_z_xxxxzz, g_0_x_z_xxxyyy, g_0_x_z_xxxyyz, g_0_x_z_xxxyzz, g_0_x_z_xxxzzz, g_0_x_z_xxyyyy, g_0_x_z_xxyyyz, g_0_x_z_xxyyzz, g_0_x_z_xxyzzz, g_0_x_z_xxzzzz, g_0_x_z_xyyyyy, g_0_x_z_xyyyyz, g_0_x_z_xyyyzz, g_0_x_z_xyyzzz, g_0_x_z_xyzzzz, g_0_x_z_xzzzzz, g_0_x_z_yyyyyy, g_0_x_z_yyyyyz, g_0_x_z_yyyyzz, g_0_x_z_yyyzzz, g_0_x_z_yyzzzz, g_0_x_z_yzzzzz, g_0_x_z_zzzzzz, g_z_x_z_xxxxxx, g_z_x_z_xxxxxxz, g_z_x_z_xxxxxy, g_z_x_z_xxxxxyz, g_z_x_z_xxxxxz, g_z_x_z_xxxxxzz, g_z_x_z_xxxxyy, g_z_x_z_xxxxyyz, g_z_x_z_xxxxyz, g_z_x_z_xxxxyzz, g_z_x_z_xxxxzz, g_z_x_z_xxxxzzz, g_z_x_z_xxxyyy, g_z_x_z_xxxyyyz, g_z_x_z_xxxyyz, g_z_x_z_xxxyyzz, g_z_x_z_xxxyzz, g_z_x_z_xxxyzzz, g_z_x_z_xxxzzz, g_z_x_z_xxxzzzz, g_z_x_z_xxyyyy, g_z_x_z_xxyyyyz, g_z_x_z_xxyyyz, g_z_x_z_xxyyyzz, g_z_x_z_xxyyzz, g_z_x_z_xxyyzzz, g_z_x_z_xxyzzz, g_z_x_z_xxyzzzz, g_z_x_z_xxzzzz, g_z_x_z_xxzzzzz, g_z_x_z_xyyyyy, g_z_x_z_xyyyyyz, g_z_x_z_xyyyyz, g_z_x_z_xyyyyzz, g_z_x_z_xyyyzz, g_z_x_z_xyyyzzz, g_z_x_z_xyyzzz, g_z_x_z_xyyzzzz, g_z_x_z_xyzzzz, g_z_x_z_xyzzzzz, g_z_x_z_xzzzzz, g_z_x_z_xzzzzzz, g_z_x_z_yyyyyy, g_z_x_z_yyyyyyz, g_z_x_z_yyyyyz, g_z_x_z_yyyyyzz, g_z_x_z_yyyyzz, g_z_x_z_yyyyzzz, g_z_x_z_yyyzzz, g_z_x_z_yyyzzzz, g_z_x_z_yyzzzz, g_z_x_z_yyzzzzz, g_z_x_z_yzzzzz, g_z_x_z_yzzzzzz, g_z_x_z_zzzzzz, g_z_x_z_zzzzzzz, g_z_x_zz_xxxxxx, g_z_x_zz_xxxxxy, g_z_x_zz_xxxxxz, g_z_x_zz_xxxxyy, g_z_x_zz_xxxxyz, g_z_x_zz_xxxxzz, g_z_x_zz_xxxyyy, g_z_x_zz_xxxyyz, g_z_x_zz_xxxyzz, g_z_x_zz_xxxzzz, g_z_x_zz_xxyyyy, g_z_x_zz_xxyyyz, g_z_x_zz_xxyyzz, g_z_x_zz_xxyzzz, g_z_x_zz_xxzzzz, g_z_x_zz_xyyyyy, g_z_x_zz_xyyyyz, g_z_x_zz_xyyyzz, g_z_x_zz_xyyzzz, g_z_x_zz_xyzzzz, g_z_x_zz_xzzzzz, g_z_x_zz_yyyyyy, g_z_x_zz_yyyyyz, g_z_x_zz_yyyyzz, g_z_x_zz_yyyzzz, g_z_x_zz_yyzzzz, g_z_x_zz_yzzzzz, g_z_x_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zz_xxxxxx[k] = -g_0_x_z_xxxxxx[k] - g_z_x_z_xxxxxx[k] * ab_z + g_z_x_z_xxxxxxz[k];

                g_z_x_zz_xxxxxy[k] = -g_0_x_z_xxxxxy[k] - g_z_x_z_xxxxxy[k] * ab_z + g_z_x_z_xxxxxyz[k];

                g_z_x_zz_xxxxxz[k] = -g_0_x_z_xxxxxz[k] - g_z_x_z_xxxxxz[k] * ab_z + g_z_x_z_xxxxxzz[k];

                g_z_x_zz_xxxxyy[k] = -g_0_x_z_xxxxyy[k] - g_z_x_z_xxxxyy[k] * ab_z + g_z_x_z_xxxxyyz[k];

                g_z_x_zz_xxxxyz[k] = -g_0_x_z_xxxxyz[k] - g_z_x_z_xxxxyz[k] * ab_z + g_z_x_z_xxxxyzz[k];

                g_z_x_zz_xxxxzz[k] = -g_0_x_z_xxxxzz[k] - g_z_x_z_xxxxzz[k] * ab_z + g_z_x_z_xxxxzzz[k];

                g_z_x_zz_xxxyyy[k] = -g_0_x_z_xxxyyy[k] - g_z_x_z_xxxyyy[k] * ab_z + g_z_x_z_xxxyyyz[k];

                g_z_x_zz_xxxyyz[k] = -g_0_x_z_xxxyyz[k] - g_z_x_z_xxxyyz[k] * ab_z + g_z_x_z_xxxyyzz[k];

                g_z_x_zz_xxxyzz[k] = -g_0_x_z_xxxyzz[k] - g_z_x_z_xxxyzz[k] * ab_z + g_z_x_z_xxxyzzz[k];

                g_z_x_zz_xxxzzz[k] = -g_0_x_z_xxxzzz[k] - g_z_x_z_xxxzzz[k] * ab_z + g_z_x_z_xxxzzzz[k];

                g_z_x_zz_xxyyyy[k] = -g_0_x_z_xxyyyy[k] - g_z_x_z_xxyyyy[k] * ab_z + g_z_x_z_xxyyyyz[k];

                g_z_x_zz_xxyyyz[k] = -g_0_x_z_xxyyyz[k] - g_z_x_z_xxyyyz[k] * ab_z + g_z_x_z_xxyyyzz[k];

                g_z_x_zz_xxyyzz[k] = -g_0_x_z_xxyyzz[k] - g_z_x_z_xxyyzz[k] * ab_z + g_z_x_z_xxyyzzz[k];

                g_z_x_zz_xxyzzz[k] = -g_0_x_z_xxyzzz[k] - g_z_x_z_xxyzzz[k] * ab_z + g_z_x_z_xxyzzzz[k];

                g_z_x_zz_xxzzzz[k] = -g_0_x_z_xxzzzz[k] - g_z_x_z_xxzzzz[k] * ab_z + g_z_x_z_xxzzzzz[k];

                g_z_x_zz_xyyyyy[k] = -g_0_x_z_xyyyyy[k] - g_z_x_z_xyyyyy[k] * ab_z + g_z_x_z_xyyyyyz[k];

                g_z_x_zz_xyyyyz[k] = -g_0_x_z_xyyyyz[k] - g_z_x_z_xyyyyz[k] * ab_z + g_z_x_z_xyyyyzz[k];

                g_z_x_zz_xyyyzz[k] = -g_0_x_z_xyyyzz[k] - g_z_x_z_xyyyzz[k] * ab_z + g_z_x_z_xyyyzzz[k];

                g_z_x_zz_xyyzzz[k] = -g_0_x_z_xyyzzz[k] - g_z_x_z_xyyzzz[k] * ab_z + g_z_x_z_xyyzzzz[k];

                g_z_x_zz_xyzzzz[k] = -g_0_x_z_xyzzzz[k] - g_z_x_z_xyzzzz[k] * ab_z + g_z_x_z_xyzzzzz[k];

                g_z_x_zz_xzzzzz[k] = -g_0_x_z_xzzzzz[k] - g_z_x_z_xzzzzz[k] * ab_z + g_z_x_z_xzzzzzz[k];

                g_z_x_zz_yyyyyy[k] = -g_0_x_z_yyyyyy[k] - g_z_x_z_yyyyyy[k] * ab_z + g_z_x_z_yyyyyyz[k];

                g_z_x_zz_yyyyyz[k] = -g_0_x_z_yyyyyz[k] - g_z_x_z_yyyyyz[k] * ab_z + g_z_x_z_yyyyyzz[k];

                g_z_x_zz_yyyyzz[k] = -g_0_x_z_yyyyzz[k] - g_z_x_z_yyyyzz[k] * ab_z + g_z_x_z_yyyyzzz[k];

                g_z_x_zz_yyyzzz[k] = -g_0_x_z_yyyzzz[k] - g_z_x_z_yyyzzz[k] * ab_z + g_z_x_z_yyyzzzz[k];

                g_z_x_zz_yyzzzz[k] = -g_0_x_z_yyzzzz[k] - g_z_x_z_yyzzzz[k] * ab_z + g_z_x_z_yyzzzzz[k];

                g_z_x_zz_yzzzzz[k] = -g_0_x_z_yzzzzz[k] - g_z_x_z_yzzzzz[k] * ab_z + g_z_x_z_yzzzzzz[k];

                g_z_x_zz_zzzzzz[k] = -g_0_x_z_zzzzzz[k] - g_z_x_z_zzzzzz[k] * ab_z + g_z_x_z_zzzzzzz[k];
            }

            /// Set up 1176-1204 components of targeted buffer : cbuffer.data(

            auto g_z_y_xx_xxxxxx = cbuffer.data(di_geom_11_off + 1176 * ccomps * dcomps);

            auto g_z_y_xx_xxxxxy = cbuffer.data(di_geom_11_off + 1177 * ccomps * dcomps);

            auto g_z_y_xx_xxxxxz = cbuffer.data(di_geom_11_off + 1178 * ccomps * dcomps);

            auto g_z_y_xx_xxxxyy = cbuffer.data(di_geom_11_off + 1179 * ccomps * dcomps);

            auto g_z_y_xx_xxxxyz = cbuffer.data(di_geom_11_off + 1180 * ccomps * dcomps);

            auto g_z_y_xx_xxxxzz = cbuffer.data(di_geom_11_off + 1181 * ccomps * dcomps);

            auto g_z_y_xx_xxxyyy = cbuffer.data(di_geom_11_off + 1182 * ccomps * dcomps);

            auto g_z_y_xx_xxxyyz = cbuffer.data(di_geom_11_off + 1183 * ccomps * dcomps);

            auto g_z_y_xx_xxxyzz = cbuffer.data(di_geom_11_off + 1184 * ccomps * dcomps);

            auto g_z_y_xx_xxxzzz = cbuffer.data(di_geom_11_off + 1185 * ccomps * dcomps);

            auto g_z_y_xx_xxyyyy = cbuffer.data(di_geom_11_off + 1186 * ccomps * dcomps);

            auto g_z_y_xx_xxyyyz = cbuffer.data(di_geom_11_off + 1187 * ccomps * dcomps);

            auto g_z_y_xx_xxyyzz = cbuffer.data(di_geom_11_off + 1188 * ccomps * dcomps);

            auto g_z_y_xx_xxyzzz = cbuffer.data(di_geom_11_off + 1189 * ccomps * dcomps);

            auto g_z_y_xx_xxzzzz = cbuffer.data(di_geom_11_off + 1190 * ccomps * dcomps);

            auto g_z_y_xx_xyyyyy = cbuffer.data(di_geom_11_off + 1191 * ccomps * dcomps);

            auto g_z_y_xx_xyyyyz = cbuffer.data(di_geom_11_off + 1192 * ccomps * dcomps);

            auto g_z_y_xx_xyyyzz = cbuffer.data(di_geom_11_off + 1193 * ccomps * dcomps);

            auto g_z_y_xx_xyyzzz = cbuffer.data(di_geom_11_off + 1194 * ccomps * dcomps);

            auto g_z_y_xx_xyzzzz = cbuffer.data(di_geom_11_off + 1195 * ccomps * dcomps);

            auto g_z_y_xx_xzzzzz = cbuffer.data(di_geom_11_off + 1196 * ccomps * dcomps);

            auto g_z_y_xx_yyyyyy = cbuffer.data(di_geom_11_off + 1197 * ccomps * dcomps);

            auto g_z_y_xx_yyyyyz = cbuffer.data(di_geom_11_off + 1198 * ccomps * dcomps);

            auto g_z_y_xx_yyyyzz = cbuffer.data(di_geom_11_off + 1199 * ccomps * dcomps);

            auto g_z_y_xx_yyyzzz = cbuffer.data(di_geom_11_off + 1200 * ccomps * dcomps);

            auto g_z_y_xx_yyzzzz = cbuffer.data(di_geom_11_off + 1201 * ccomps * dcomps);

            auto g_z_y_xx_yzzzzz = cbuffer.data(di_geom_11_off + 1202 * ccomps * dcomps);

            auto g_z_y_xx_zzzzzz = cbuffer.data(di_geom_11_off + 1203 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_x_xxxxxx, g_z_y_x_xxxxxxx, g_z_y_x_xxxxxxy, g_z_y_x_xxxxxxz, g_z_y_x_xxxxxy, g_z_y_x_xxxxxyy, g_z_y_x_xxxxxyz, g_z_y_x_xxxxxz, g_z_y_x_xxxxxzz, g_z_y_x_xxxxyy, g_z_y_x_xxxxyyy, g_z_y_x_xxxxyyz, g_z_y_x_xxxxyz, g_z_y_x_xxxxyzz, g_z_y_x_xxxxzz, g_z_y_x_xxxxzzz, g_z_y_x_xxxyyy, g_z_y_x_xxxyyyy, g_z_y_x_xxxyyyz, g_z_y_x_xxxyyz, g_z_y_x_xxxyyzz, g_z_y_x_xxxyzz, g_z_y_x_xxxyzzz, g_z_y_x_xxxzzz, g_z_y_x_xxxzzzz, g_z_y_x_xxyyyy, g_z_y_x_xxyyyyy, g_z_y_x_xxyyyyz, g_z_y_x_xxyyyz, g_z_y_x_xxyyyzz, g_z_y_x_xxyyzz, g_z_y_x_xxyyzzz, g_z_y_x_xxyzzz, g_z_y_x_xxyzzzz, g_z_y_x_xxzzzz, g_z_y_x_xxzzzzz, g_z_y_x_xyyyyy, g_z_y_x_xyyyyyy, g_z_y_x_xyyyyyz, g_z_y_x_xyyyyz, g_z_y_x_xyyyyzz, g_z_y_x_xyyyzz, g_z_y_x_xyyyzzz, g_z_y_x_xyyzzz, g_z_y_x_xyyzzzz, g_z_y_x_xyzzzz, g_z_y_x_xyzzzzz, g_z_y_x_xzzzzz, g_z_y_x_xzzzzzz, g_z_y_x_yyyyyy, g_z_y_x_yyyyyz, g_z_y_x_yyyyzz, g_z_y_x_yyyzzz, g_z_y_x_yyzzzz, g_z_y_x_yzzzzz, g_z_y_x_zzzzzz, g_z_y_xx_xxxxxx, g_z_y_xx_xxxxxy, g_z_y_xx_xxxxxz, g_z_y_xx_xxxxyy, g_z_y_xx_xxxxyz, g_z_y_xx_xxxxzz, g_z_y_xx_xxxyyy, g_z_y_xx_xxxyyz, g_z_y_xx_xxxyzz, g_z_y_xx_xxxzzz, g_z_y_xx_xxyyyy, g_z_y_xx_xxyyyz, g_z_y_xx_xxyyzz, g_z_y_xx_xxyzzz, g_z_y_xx_xxzzzz, g_z_y_xx_xyyyyy, g_z_y_xx_xyyyyz, g_z_y_xx_xyyyzz, g_z_y_xx_xyyzzz, g_z_y_xx_xyzzzz, g_z_y_xx_xzzzzz, g_z_y_xx_yyyyyy, g_z_y_xx_yyyyyz, g_z_y_xx_yyyyzz, g_z_y_xx_yyyzzz, g_z_y_xx_yyzzzz, g_z_y_xx_yzzzzz, g_z_y_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xx_xxxxxx[k] = -g_z_y_x_xxxxxx[k] * ab_x + g_z_y_x_xxxxxxx[k];

                g_z_y_xx_xxxxxy[k] = -g_z_y_x_xxxxxy[k] * ab_x + g_z_y_x_xxxxxxy[k];

                g_z_y_xx_xxxxxz[k] = -g_z_y_x_xxxxxz[k] * ab_x + g_z_y_x_xxxxxxz[k];

                g_z_y_xx_xxxxyy[k] = -g_z_y_x_xxxxyy[k] * ab_x + g_z_y_x_xxxxxyy[k];

                g_z_y_xx_xxxxyz[k] = -g_z_y_x_xxxxyz[k] * ab_x + g_z_y_x_xxxxxyz[k];

                g_z_y_xx_xxxxzz[k] = -g_z_y_x_xxxxzz[k] * ab_x + g_z_y_x_xxxxxzz[k];

                g_z_y_xx_xxxyyy[k] = -g_z_y_x_xxxyyy[k] * ab_x + g_z_y_x_xxxxyyy[k];

                g_z_y_xx_xxxyyz[k] = -g_z_y_x_xxxyyz[k] * ab_x + g_z_y_x_xxxxyyz[k];

                g_z_y_xx_xxxyzz[k] = -g_z_y_x_xxxyzz[k] * ab_x + g_z_y_x_xxxxyzz[k];

                g_z_y_xx_xxxzzz[k] = -g_z_y_x_xxxzzz[k] * ab_x + g_z_y_x_xxxxzzz[k];

                g_z_y_xx_xxyyyy[k] = -g_z_y_x_xxyyyy[k] * ab_x + g_z_y_x_xxxyyyy[k];

                g_z_y_xx_xxyyyz[k] = -g_z_y_x_xxyyyz[k] * ab_x + g_z_y_x_xxxyyyz[k];

                g_z_y_xx_xxyyzz[k] = -g_z_y_x_xxyyzz[k] * ab_x + g_z_y_x_xxxyyzz[k];

                g_z_y_xx_xxyzzz[k] = -g_z_y_x_xxyzzz[k] * ab_x + g_z_y_x_xxxyzzz[k];

                g_z_y_xx_xxzzzz[k] = -g_z_y_x_xxzzzz[k] * ab_x + g_z_y_x_xxxzzzz[k];

                g_z_y_xx_xyyyyy[k] = -g_z_y_x_xyyyyy[k] * ab_x + g_z_y_x_xxyyyyy[k];

                g_z_y_xx_xyyyyz[k] = -g_z_y_x_xyyyyz[k] * ab_x + g_z_y_x_xxyyyyz[k];

                g_z_y_xx_xyyyzz[k] = -g_z_y_x_xyyyzz[k] * ab_x + g_z_y_x_xxyyyzz[k];

                g_z_y_xx_xyyzzz[k] = -g_z_y_x_xyyzzz[k] * ab_x + g_z_y_x_xxyyzzz[k];

                g_z_y_xx_xyzzzz[k] = -g_z_y_x_xyzzzz[k] * ab_x + g_z_y_x_xxyzzzz[k];

                g_z_y_xx_xzzzzz[k] = -g_z_y_x_xzzzzz[k] * ab_x + g_z_y_x_xxzzzzz[k];

                g_z_y_xx_yyyyyy[k] = -g_z_y_x_yyyyyy[k] * ab_x + g_z_y_x_xyyyyyy[k];

                g_z_y_xx_yyyyyz[k] = -g_z_y_x_yyyyyz[k] * ab_x + g_z_y_x_xyyyyyz[k];

                g_z_y_xx_yyyyzz[k] = -g_z_y_x_yyyyzz[k] * ab_x + g_z_y_x_xyyyyzz[k];

                g_z_y_xx_yyyzzz[k] = -g_z_y_x_yyyzzz[k] * ab_x + g_z_y_x_xyyyzzz[k];

                g_z_y_xx_yyzzzz[k] = -g_z_y_x_yyzzzz[k] * ab_x + g_z_y_x_xyyzzzz[k];

                g_z_y_xx_yzzzzz[k] = -g_z_y_x_yzzzzz[k] * ab_x + g_z_y_x_xyzzzzz[k];

                g_z_y_xx_zzzzzz[k] = -g_z_y_x_zzzzzz[k] * ab_x + g_z_y_x_xzzzzzz[k];
            }

            /// Set up 1204-1232 components of targeted buffer : cbuffer.data(

            auto g_z_y_xy_xxxxxx = cbuffer.data(di_geom_11_off + 1204 * ccomps * dcomps);

            auto g_z_y_xy_xxxxxy = cbuffer.data(di_geom_11_off + 1205 * ccomps * dcomps);

            auto g_z_y_xy_xxxxxz = cbuffer.data(di_geom_11_off + 1206 * ccomps * dcomps);

            auto g_z_y_xy_xxxxyy = cbuffer.data(di_geom_11_off + 1207 * ccomps * dcomps);

            auto g_z_y_xy_xxxxyz = cbuffer.data(di_geom_11_off + 1208 * ccomps * dcomps);

            auto g_z_y_xy_xxxxzz = cbuffer.data(di_geom_11_off + 1209 * ccomps * dcomps);

            auto g_z_y_xy_xxxyyy = cbuffer.data(di_geom_11_off + 1210 * ccomps * dcomps);

            auto g_z_y_xy_xxxyyz = cbuffer.data(di_geom_11_off + 1211 * ccomps * dcomps);

            auto g_z_y_xy_xxxyzz = cbuffer.data(di_geom_11_off + 1212 * ccomps * dcomps);

            auto g_z_y_xy_xxxzzz = cbuffer.data(di_geom_11_off + 1213 * ccomps * dcomps);

            auto g_z_y_xy_xxyyyy = cbuffer.data(di_geom_11_off + 1214 * ccomps * dcomps);

            auto g_z_y_xy_xxyyyz = cbuffer.data(di_geom_11_off + 1215 * ccomps * dcomps);

            auto g_z_y_xy_xxyyzz = cbuffer.data(di_geom_11_off + 1216 * ccomps * dcomps);

            auto g_z_y_xy_xxyzzz = cbuffer.data(di_geom_11_off + 1217 * ccomps * dcomps);

            auto g_z_y_xy_xxzzzz = cbuffer.data(di_geom_11_off + 1218 * ccomps * dcomps);

            auto g_z_y_xy_xyyyyy = cbuffer.data(di_geom_11_off + 1219 * ccomps * dcomps);

            auto g_z_y_xy_xyyyyz = cbuffer.data(di_geom_11_off + 1220 * ccomps * dcomps);

            auto g_z_y_xy_xyyyzz = cbuffer.data(di_geom_11_off + 1221 * ccomps * dcomps);

            auto g_z_y_xy_xyyzzz = cbuffer.data(di_geom_11_off + 1222 * ccomps * dcomps);

            auto g_z_y_xy_xyzzzz = cbuffer.data(di_geom_11_off + 1223 * ccomps * dcomps);

            auto g_z_y_xy_xzzzzz = cbuffer.data(di_geom_11_off + 1224 * ccomps * dcomps);

            auto g_z_y_xy_yyyyyy = cbuffer.data(di_geom_11_off + 1225 * ccomps * dcomps);

            auto g_z_y_xy_yyyyyz = cbuffer.data(di_geom_11_off + 1226 * ccomps * dcomps);

            auto g_z_y_xy_yyyyzz = cbuffer.data(di_geom_11_off + 1227 * ccomps * dcomps);

            auto g_z_y_xy_yyyzzz = cbuffer.data(di_geom_11_off + 1228 * ccomps * dcomps);

            auto g_z_y_xy_yyzzzz = cbuffer.data(di_geom_11_off + 1229 * ccomps * dcomps);

            auto g_z_y_xy_yzzzzz = cbuffer.data(di_geom_11_off + 1230 * ccomps * dcomps);

            auto g_z_y_xy_zzzzzz = cbuffer.data(di_geom_11_off + 1231 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xy_xxxxxx, g_z_y_xy_xxxxxy, g_z_y_xy_xxxxxz, g_z_y_xy_xxxxyy, g_z_y_xy_xxxxyz, g_z_y_xy_xxxxzz, g_z_y_xy_xxxyyy, g_z_y_xy_xxxyyz, g_z_y_xy_xxxyzz, g_z_y_xy_xxxzzz, g_z_y_xy_xxyyyy, g_z_y_xy_xxyyyz, g_z_y_xy_xxyyzz, g_z_y_xy_xxyzzz, g_z_y_xy_xxzzzz, g_z_y_xy_xyyyyy, g_z_y_xy_xyyyyz, g_z_y_xy_xyyyzz, g_z_y_xy_xyyzzz, g_z_y_xy_xyzzzz, g_z_y_xy_xzzzzz, g_z_y_xy_yyyyyy, g_z_y_xy_yyyyyz, g_z_y_xy_yyyyzz, g_z_y_xy_yyyzzz, g_z_y_xy_yyzzzz, g_z_y_xy_yzzzzz, g_z_y_xy_zzzzzz, g_z_y_y_xxxxxx, g_z_y_y_xxxxxxx, g_z_y_y_xxxxxxy, g_z_y_y_xxxxxxz, g_z_y_y_xxxxxy, g_z_y_y_xxxxxyy, g_z_y_y_xxxxxyz, g_z_y_y_xxxxxz, g_z_y_y_xxxxxzz, g_z_y_y_xxxxyy, g_z_y_y_xxxxyyy, g_z_y_y_xxxxyyz, g_z_y_y_xxxxyz, g_z_y_y_xxxxyzz, g_z_y_y_xxxxzz, g_z_y_y_xxxxzzz, g_z_y_y_xxxyyy, g_z_y_y_xxxyyyy, g_z_y_y_xxxyyyz, g_z_y_y_xxxyyz, g_z_y_y_xxxyyzz, g_z_y_y_xxxyzz, g_z_y_y_xxxyzzz, g_z_y_y_xxxzzz, g_z_y_y_xxxzzzz, g_z_y_y_xxyyyy, g_z_y_y_xxyyyyy, g_z_y_y_xxyyyyz, g_z_y_y_xxyyyz, g_z_y_y_xxyyyzz, g_z_y_y_xxyyzz, g_z_y_y_xxyyzzz, g_z_y_y_xxyzzz, g_z_y_y_xxyzzzz, g_z_y_y_xxzzzz, g_z_y_y_xxzzzzz, g_z_y_y_xyyyyy, g_z_y_y_xyyyyyy, g_z_y_y_xyyyyyz, g_z_y_y_xyyyyz, g_z_y_y_xyyyyzz, g_z_y_y_xyyyzz, g_z_y_y_xyyyzzz, g_z_y_y_xyyzzz, g_z_y_y_xyyzzzz, g_z_y_y_xyzzzz, g_z_y_y_xyzzzzz, g_z_y_y_xzzzzz, g_z_y_y_xzzzzzz, g_z_y_y_yyyyyy, g_z_y_y_yyyyyz, g_z_y_y_yyyyzz, g_z_y_y_yyyzzz, g_z_y_y_yyzzzz, g_z_y_y_yzzzzz, g_z_y_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xy_xxxxxx[k] = -g_z_y_y_xxxxxx[k] * ab_x + g_z_y_y_xxxxxxx[k];

                g_z_y_xy_xxxxxy[k] = -g_z_y_y_xxxxxy[k] * ab_x + g_z_y_y_xxxxxxy[k];

                g_z_y_xy_xxxxxz[k] = -g_z_y_y_xxxxxz[k] * ab_x + g_z_y_y_xxxxxxz[k];

                g_z_y_xy_xxxxyy[k] = -g_z_y_y_xxxxyy[k] * ab_x + g_z_y_y_xxxxxyy[k];

                g_z_y_xy_xxxxyz[k] = -g_z_y_y_xxxxyz[k] * ab_x + g_z_y_y_xxxxxyz[k];

                g_z_y_xy_xxxxzz[k] = -g_z_y_y_xxxxzz[k] * ab_x + g_z_y_y_xxxxxzz[k];

                g_z_y_xy_xxxyyy[k] = -g_z_y_y_xxxyyy[k] * ab_x + g_z_y_y_xxxxyyy[k];

                g_z_y_xy_xxxyyz[k] = -g_z_y_y_xxxyyz[k] * ab_x + g_z_y_y_xxxxyyz[k];

                g_z_y_xy_xxxyzz[k] = -g_z_y_y_xxxyzz[k] * ab_x + g_z_y_y_xxxxyzz[k];

                g_z_y_xy_xxxzzz[k] = -g_z_y_y_xxxzzz[k] * ab_x + g_z_y_y_xxxxzzz[k];

                g_z_y_xy_xxyyyy[k] = -g_z_y_y_xxyyyy[k] * ab_x + g_z_y_y_xxxyyyy[k];

                g_z_y_xy_xxyyyz[k] = -g_z_y_y_xxyyyz[k] * ab_x + g_z_y_y_xxxyyyz[k];

                g_z_y_xy_xxyyzz[k] = -g_z_y_y_xxyyzz[k] * ab_x + g_z_y_y_xxxyyzz[k];

                g_z_y_xy_xxyzzz[k] = -g_z_y_y_xxyzzz[k] * ab_x + g_z_y_y_xxxyzzz[k];

                g_z_y_xy_xxzzzz[k] = -g_z_y_y_xxzzzz[k] * ab_x + g_z_y_y_xxxzzzz[k];

                g_z_y_xy_xyyyyy[k] = -g_z_y_y_xyyyyy[k] * ab_x + g_z_y_y_xxyyyyy[k];

                g_z_y_xy_xyyyyz[k] = -g_z_y_y_xyyyyz[k] * ab_x + g_z_y_y_xxyyyyz[k];

                g_z_y_xy_xyyyzz[k] = -g_z_y_y_xyyyzz[k] * ab_x + g_z_y_y_xxyyyzz[k];

                g_z_y_xy_xyyzzz[k] = -g_z_y_y_xyyzzz[k] * ab_x + g_z_y_y_xxyyzzz[k];

                g_z_y_xy_xyzzzz[k] = -g_z_y_y_xyzzzz[k] * ab_x + g_z_y_y_xxyzzzz[k];

                g_z_y_xy_xzzzzz[k] = -g_z_y_y_xzzzzz[k] * ab_x + g_z_y_y_xxzzzzz[k];

                g_z_y_xy_yyyyyy[k] = -g_z_y_y_yyyyyy[k] * ab_x + g_z_y_y_xyyyyyy[k];

                g_z_y_xy_yyyyyz[k] = -g_z_y_y_yyyyyz[k] * ab_x + g_z_y_y_xyyyyyz[k];

                g_z_y_xy_yyyyzz[k] = -g_z_y_y_yyyyzz[k] * ab_x + g_z_y_y_xyyyyzz[k];

                g_z_y_xy_yyyzzz[k] = -g_z_y_y_yyyzzz[k] * ab_x + g_z_y_y_xyyyzzz[k];

                g_z_y_xy_yyzzzz[k] = -g_z_y_y_yyzzzz[k] * ab_x + g_z_y_y_xyyzzzz[k];

                g_z_y_xy_yzzzzz[k] = -g_z_y_y_yzzzzz[k] * ab_x + g_z_y_y_xyzzzzz[k];

                g_z_y_xy_zzzzzz[k] = -g_z_y_y_zzzzzz[k] * ab_x + g_z_y_y_xzzzzzz[k];
            }

            /// Set up 1232-1260 components of targeted buffer : cbuffer.data(

            auto g_z_y_xz_xxxxxx = cbuffer.data(di_geom_11_off + 1232 * ccomps * dcomps);

            auto g_z_y_xz_xxxxxy = cbuffer.data(di_geom_11_off + 1233 * ccomps * dcomps);

            auto g_z_y_xz_xxxxxz = cbuffer.data(di_geom_11_off + 1234 * ccomps * dcomps);

            auto g_z_y_xz_xxxxyy = cbuffer.data(di_geom_11_off + 1235 * ccomps * dcomps);

            auto g_z_y_xz_xxxxyz = cbuffer.data(di_geom_11_off + 1236 * ccomps * dcomps);

            auto g_z_y_xz_xxxxzz = cbuffer.data(di_geom_11_off + 1237 * ccomps * dcomps);

            auto g_z_y_xz_xxxyyy = cbuffer.data(di_geom_11_off + 1238 * ccomps * dcomps);

            auto g_z_y_xz_xxxyyz = cbuffer.data(di_geom_11_off + 1239 * ccomps * dcomps);

            auto g_z_y_xz_xxxyzz = cbuffer.data(di_geom_11_off + 1240 * ccomps * dcomps);

            auto g_z_y_xz_xxxzzz = cbuffer.data(di_geom_11_off + 1241 * ccomps * dcomps);

            auto g_z_y_xz_xxyyyy = cbuffer.data(di_geom_11_off + 1242 * ccomps * dcomps);

            auto g_z_y_xz_xxyyyz = cbuffer.data(di_geom_11_off + 1243 * ccomps * dcomps);

            auto g_z_y_xz_xxyyzz = cbuffer.data(di_geom_11_off + 1244 * ccomps * dcomps);

            auto g_z_y_xz_xxyzzz = cbuffer.data(di_geom_11_off + 1245 * ccomps * dcomps);

            auto g_z_y_xz_xxzzzz = cbuffer.data(di_geom_11_off + 1246 * ccomps * dcomps);

            auto g_z_y_xz_xyyyyy = cbuffer.data(di_geom_11_off + 1247 * ccomps * dcomps);

            auto g_z_y_xz_xyyyyz = cbuffer.data(di_geom_11_off + 1248 * ccomps * dcomps);

            auto g_z_y_xz_xyyyzz = cbuffer.data(di_geom_11_off + 1249 * ccomps * dcomps);

            auto g_z_y_xz_xyyzzz = cbuffer.data(di_geom_11_off + 1250 * ccomps * dcomps);

            auto g_z_y_xz_xyzzzz = cbuffer.data(di_geom_11_off + 1251 * ccomps * dcomps);

            auto g_z_y_xz_xzzzzz = cbuffer.data(di_geom_11_off + 1252 * ccomps * dcomps);

            auto g_z_y_xz_yyyyyy = cbuffer.data(di_geom_11_off + 1253 * ccomps * dcomps);

            auto g_z_y_xz_yyyyyz = cbuffer.data(di_geom_11_off + 1254 * ccomps * dcomps);

            auto g_z_y_xz_yyyyzz = cbuffer.data(di_geom_11_off + 1255 * ccomps * dcomps);

            auto g_z_y_xz_yyyzzz = cbuffer.data(di_geom_11_off + 1256 * ccomps * dcomps);

            auto g_z_y_xz_yyzzzz = cbuffer.data(di_geom_11_off + 1257 * ccomps * dcomps);

            auto g_z_y_xz_yzzzzz = cbuffer.data(di_geom_11_off + 1258 * ccomps * dcomps);

            auto g_z_y_xz_zzzzzz = cbuffer.data(di_geom_11_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xz_xxxxxx, g_z_y_xz_xxxxxy, g_z_y_xz_xxxxxz, g_z_y_xz_xxxxyy, g_z_y_xz_xxxxyz, g_z_y_xz_xxxxzz, g_z_y_xz_xxxyyy, g_z_y_xz_xxxyyz, g_z_y_xz_xxxyzz, g_z_y_xz_xxxzzz, g_z_y_xz_xxyyyy, g_z_y_xz_xxyyyz, g_z_y_xz_xxyyzz, g_z_y_xz_xxyzzz, g_z_y_xz_xxzzzz, g_z_y_xz_xyyyyy, g_z_y_xz_xyyyyz, g_z_y_xz_xyyyzz, g_z_y_xz_xyyzzz, g_z_y_xz_xyzzzz, g_z_y_xz_xzzzzz, g_z_y_xz_yyyyyy, g_z_y_xz_yyyyyz, g_z_y_xz_yyyyzz, g_z_y_xz_yyyzzz, g_z_y_xz_yyzzzz, g_z_y_xz_yzzzzz, g_z_y_xz_zzzzzz, g_z_y_z_xxxxxx, g_z_y_z_xxxxxxx, g_z_y_z_xxxxxxy, g_z_y_z_xxxxxxz, g_z_y_z_xxxxxy, g_z_y_z_xxxxxyy, g_z_y_z_xxxxxyz, g_z_y_z_xxxxxz, g_z_y_z_xxxxxzz, g_z_y_z_xxxxyy, g_z_y_z_xxxxyyy, g_z_y_z_xxxxyyz, g_z_y_z_xxxxyz, g_z_y_z_xxxxyzz, g_z_y_z_xxxxzz, g_z_y_z_xxxxzzz, g_z_y_z_xxxyyy, g_z_y_z_xxxyyyy, g_z_y_z_xxxyyyz, g_z_y_z_xxxyyz, g_z_y_z_xxxyyzz, g_z_y_z_xxxyzz, g_z_y_z_xxxyzzz, g_z_y_z_xxxzzz, g_z_y_z_xxxzzzz, g_z_y_z_xxyyyy, g_z_y_z_xxyyyyy, g_z_y_z_xxyyyyz, g_z_y_z_xxyyyz, g_z_y_z_xxyyyzz, g_z_y_z_xxyyzz, g_z_y_z_xxyyzzz, g_z_y_z_xxyzzz, g_z_y_z_xxyzzzz, g_z_y_z_xxzzzz, g_z_y_z_xxzzzzz, g_z_y_z_xyyyyy, g_z_y_z_xyyyyyy, g_z_y_z_xyyyyyz, g_z_y_z_xyyyyz, g_z_y_z_xyyyyzz, g_z_y_z_xyyyzz, g_z_y_z_xyyyzzz, g_z_y_z_xyyzzz, g_z_y_z_xyyzzzz, g_z_y_z_xyzzzz, g_z_y_z_xyzzzzz, g_z_y_z_xzzzzz, g_z_y_z_xzzzzzz, g_z_y_z_yyyyyy, g_z_y_z_yyyyyz, g_z_y_z_yyyyzz, g_z_y_z_yyyzzz, g_z_y_z_yyzzzz, g_z_y_z_yzzzzz, g_z_y_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xz_xxxxxx[k] = -g_z_y_z_xxxxxx[k] * ab_x + g_z_y_z_xxxxxxx[k];

                g_z_y_xz_xxxxxy[k] = -g_z_y_z_xxxxxy[k] * ab_x + g_z_y_z_xxxxxxy[k];

                g_z_y_xz_xxxxxz[k] = -g_z_y_z_xxxxxz[k] * ab_x + g_z_y_z_xxxxxxz[k];

                g_z_y_xz_xxxxyy[k] = -g_z_y_z_xxxxyy[k] * ab_x + g_z_y_z_xxxxxyy[k];

                g_z_y_xz_xxxxyz[k] = -g_z_y_z_xxxxyz[k] * ab_x + g_z_y_z_xxxxxyz[k];

                g_z_y_xz_xxxxzz[k] = -g_z_y_z_xxxxzz[k] * ab_x + g_z_y_z_xxxxxzz[k];

                g_z_y_xz_xxxyyy[k] = -g_z_y_z_xxxyyy[k] * ab_x + g_z_y_z_xxxxyyy[k];

                g_z_y_xz_xxxyyz[k] = -g_z_y_z_xxxyyz[k] * ab_x + g_z_y_z_xxxxyyz[k];

                g_z_y_xz_xxxyzz[k] = -g_z_y_z_xxxyzz[k] * ab_x + g_z_y_z_xxxxyzz[k];

                g_z_y_xz_xxxzzz[k] = -g_z_y_z_xxxzzz[k] * ab_x + g_z_y_z_xxxxzzz[k];

                g_z_y_xz_xxyyyy[k] = -g_z_y_z_xxyyyy[k] * ab_x + g_z_y_z_xxxyyyy[k];

                g_z_y_xz_xxyyyz[k] = -g_z_y_z_xxyyyz[k] * ab_x + g_z_y_z_xxxyyyz[k];

                g_z_y_xz_xxyyzz[k] = -g_z_y_z_xxyyzz[k] * ab_x + g_z_y_z_xxxyyzz[k];

                g_z_y_xz_xxyzzz[k] = -g_z_y_z_xxyzzz[k] * ab_x + g_z_y_z_xxxyzzz[k];

                g_z_y_xz_xxzzzz[k] = -g_z_y_z_xxzzzz[k] * ab_x + g_z_y_z_xxxzzzz[k];

                g_z_y_xz_xyyyyy[k] = -g_z_y_z_xyyyyy[k] * ab_x + g_z_y_z_xxyyyyy[k];

                g_z_y_xz_xyyyyz[k] = -g_z_y_z_xyyyyz[k] * ab_x + g_z_y_z_xxyyyyz[k];

                g_z_y_xz_xyyyzz[k] = -g_z_y_z_xyyyzz[k] * ab_x + g_z_y_z_xxyyyzz[k];

                g_z_y_xz_xyyzzz[k] = -g_z_y_z_xyyzzz[k] * ab_x + g_z_y_z_xxyyzzz[k];

                g_z_y_xz_xyzzzz[k] = -g_z_y_z_xyzzzz[k] * ab_x + g_z_y_z_xxyzzzz[k];

                g_z_y_xz_xzzzzz[k] = -g_z_y_z_xzzzzz[k] * ab_x + g_z_y_z_xxzzzzz[k];

                g_z_y_xz_yyyyyy[k] = -g_z_y_z_yyyyyy[k] * ab_x + g_z_y_z_xyyyyyy[k];

                g_z_y_xz_yyyyyz[k] = -g_z_y_z_yyyyyz[k] * ab_x + g_z_y_z_xyyyyyz[k];

                g_z_y_xz_yyyyzz[k] = -g_z_y_z_yyyyzz[k] * ab_x + g_z_y_z_xyyyyzz[k];

                g_z_y_xz_yyyzzz[k] = -g_z_y_z_yyyzzz[k] * ab_x + g_z_y_z_xyyyzzz[k];

                g_z_y_xz_yyzzzz[k] = -g_z_y_z_yyzzzz[k] * ab_x + g_z_y_z_xyyzzzz[k];

                g_z_y_xz_yzzzzz[k] = -g_z_y_z_yzzzzz[k] * ab_x + g_z_y_z_xyzzzzz[k];

                g_z_y_xz_zzzzzz[k] = -g_z_y_z_zzzzzz[k] * ab_x + g_z_y_z_xzzzzzz[k];
            }

            /// Set up 1260-1288 components of targeted buffer : cbuffer.data(

            auto g_z_y_yy_xxxxxx = cbuffer.data(di_geom_11_off + 1260 * ccomps * dcomps);

            auto g_z_y_yy_xxxxxy = cbuffer.data(di_geom_11_off + 1261 * ccomps * dcomps);

            auto g_z_y_yy_xxxxxz = cbuffer.data(di_geom_11_off + 1262 * ccomps * dcomps);

            auto g_z_y_yy_xxxxyy = cbuffer.data(di_geom_11_off + 1263 * ccomps * dcomps);

            auto g_z_y_yy_xxxxyz = cbuffer.data(di_geom_11_off + 1264 * ccomps * dcomps);

            auto g_z_y_yy_xxxxzz = cbuffer.data(di_geom_11_off + 1265 * ccomps * dcomps);

            auto g_z_y_yy_xxxyyy = cbuffer.data(di_geom_11_off + 1266 * ccomps * dcomps);

            auto g_z_y_yy_xxxyyz = cbuffer.data(di_geom_11_off + 1267 * ccomps * dcomps);

            auto g_z_y_yy_xxxyzz = cbuffer.data(di_geom_11_off + 1268 * ccomps * dcomps);

            auto g_z_y_yy_xxxzzz = cbuffer.data(di_geom_11_off + 1269 * ccomps * dcomps);

            auto g_z_y_yy_xxyyyy = cbuffer.data(di_geom_11_off + 1270 * ccomps * dcomps);

            auto g_z_y_yy_xxyyyz = cbuffer.data(di_geom_11_off + 1271 * ccomps * dcomps);

            auto g_z_y_yy_xxyyzz = cbuffer.data(di_geom_11_off + 1272 * ccomps * dcomps);

            auto g_z_y_yy_xxyzzz = cbuffer.data(di_geom_11_off + 1273 * ccomps * dcomps);

            auto g_z_y_yy_xxzzzz = cbuffer.data(di_geom_11_off + 1274 * ccomps * dcomps);

            auto g_z_y_yy_xyyyyy = cbuffer.data(di_geom_11_off + 1275 * ccomps * dcomps);

            auto g_z_y_yy_xyyyyz = cbuffer.data(di_geom_11_off + 1276 * ccomps * dcomps);

            auto g_z_y_yy_xyyyzz = cbuffer.data(di_geom_11_off + 1277 * ccomps * dcomps);

            auto g_z_y_yy_xyyzzz = cbuffer.data(di_geom_11_off + 1278 * ccomps * dcomps);

            auto g_z_y_yy_xyzzzz = cbuffer.data(di_geom_11_off + 1279 * ccomps * dcomps);

            auto g_z_y_yy_xzzzzz = cbuffer.data(di_geom_11_off + 1280 * ccomps * dcomps);

            auto g_z_y_yy_yyyyyy = cbuffer.data(di_geom_11_off + 1281 * ccomps * dcomps);

            auto g_z_y_yy_yyyyyz = cbuffer.data(di_geom_11_off + 1282 * ccomps * dcomps);

            auto g_z_y_yy_yyyyzz = cbuffer.data(di_geom_11_off + 1283 * ccomps * dcomps);

            auto g_z_y_yy_yyyzzz = cbuffer.data(di_geom_11_off + 1284 * ccomps * dcomps);

            auto g_z_y_yy_yyzzzz = cbuffer.data(di_geom_11_off + 1285 * ccomps * dcomps);

            auto g_z_y_yy_yzzzzz = cbuffer.data(di_geom_11_off + 1286 * ccomps * dcomps);

            auto g_z_y_yy_zzzzzz = cbuffer.data(di_geom_11_off + 1287 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_xxxxxx, g_z_0_y_xxxxxy, g_z_0_y_xxxxxz, g_z_0_y_xxxxyy, g_z_0_y_xxxxyz, g_z_0_y_xxxxzz, g_z_0_y_xxxyyy, g_z_0_y_xxxyyz, g_z_0_y_xxxyzz, g_z_0_y_xxxzzz, g_z_0_y_xxyyyy, g_z_0_y_xxyyyz, g_z_0_y_xxyyzz, g_z_0_y_xxyzzz, g_z_0_y_xxzzzz, g_z_0_y_xyyyyy, g_z_0_y_xyyyyz, g_z_0_y_xyyyzz, g_z_0_y_xyyzzz, g_z_0_y_xyzzzz, g_z_0_y_xzzzzz, g_z_0_y_yyyyyy, g_z_0_y_yyyyyz, g_z_0_y_yyyyzz, g_z_0_y_yyyzzz, g_z_0_y_yyzzzz, g_z_0_y_yzzzzz, g_z_0_y_zzzzzz, g_z_y_y_xxxxxx, g_z_y_y_xxxxxxy, g_z_y_y_xxxxxy, g_z_y_y_xxxxxyy, g_z_y_y_xxxxxyz, g_z_y_y_xxxxxz, g_z_y_y_xxxxyy, g_z_y_y_xxxxyyy, g_z_y_y_xxxxyyz, g_z_y_y_xxxxyz, g_z_y_y_xxxxyzz, g_z_y_y_xxxxzz, g_z_y_y_xxxyyy, g_z_y_y_xxxyyyy, g_z_y_y_xxxyyyz, g_z_y_y_xxxyyz, g_z_y_y_xxxyyzz, g_z_y_y_xxxyzz, g_z_y_y_xxxyzzz, g_z_y_y_xxxzzz, g_z_y_y_xxyyyy, g_z_y_y_xxyyyyy, g_z_y_y_xxyyyyz, g_z_y_y_xxyyyz, g_z_y_y_xxyyyzz, g_z_y_y_xxyyzz, g_z_y_y_xxyyzzz, g_z_y_y_xxyzzz, g_z_y_y_xxyzzzz, g_z_y_y_xxzzzz, g_z_y_y_xyyyyy, g_z_y_y_xyyyyyy, g_z_y_y_xyyyyyz, g_z_y_y_xyyyyz, g_z_y_y_xyyyyzz, g_z_y_y_xyyyzz, g_z_y_y_xyyyzzz, g_z_y_y_xyyzzz, g_z_y_y_xyyzzzz, g_z_y_y_xyzzzz, g_z_y_y_xyzzzzz, g_z_y_y_xzzzzz, g_z_y_y_yyyyyy, g_z_y_y_yyyyyyy, g_z_y_y_yyyyyyz, g_z_y_y_yyyyyz, g_z_y_y_yyyyyzz, g_z_y_y_yyyyzz, g_z_y_y_yyyyzzz, g_z_y_y_yyyzzz, g_z_y_y_yyyzzzz, g_z_y_y_yyzzzz, g_z_y_y_yyzzzzz, g_z_y_y_yzzzzz, g_z_y_y_yzzzzzz, g_z_y_y_zzzzzz, g_z_y_yy_xxxxxx, g_z_y_yy_xxxxxy, g_z_y_yy_xxxxxz, g_z_y_yy_xxxxyy, g_z_y_yy_xxxxyz, g_z_y_yy_xxxxzz, g_z_y_yy_xxxyyy, g_z_y_yy_xxxyyz, g_z_y_yy_xxxyzz, g_z_y_yy_xxxzzz, g_z_y_yy_xxyyyy, g_z_y_yy_xxyyyz, g_z_y_yy_xxyyzz, g_z_y_yy_xxyzzz, g_z_y_yy_xxzzzz, g_z_y_yy_xyyyyy, g_z_y_yy_xyyyyz, g_z_y_yy_xyyyzz, g_z_y_yy_xyyzzz, g_z_y_yy_xyzzzz, g_z_y_yy_xzzzzz, g_z_y_yy_yyyyyy, g_z_y_yy_yyyyyz, g_z_y_yy_yyyyzz, g_z_y_yy_yyyzzz, g_z_y_yy_yyzzzz, g_z_y_yy_yzzzzz, g_z_y_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yy_xxxxxx[k] = g_z_0_y_xxxxxx[k] - g_z_y_y_xxxxxx[k] * ab_y + g_z_y_y_xxxxxxy[k];

                g_z_y_yy_xxxxxy[k] = g_z_0_y_xxxxxy[k] - g_z_y_y_xxxxxy[k] * ab_y + g_z_y_y_xxxxxyy[k];

                g_z_y_yy_xxxxxz[k] = g_z_0_y_xxxxxz[k] - g_z_y_y_xxxxxz[k] * ab_y + g_z_y_y_xxxxxyz[k];

                g_z_y_yy_xxxxyy[k] = g_z_0_y_xxxxyy[k] - g_z_y_y_xxxxyy[k] * ab_y + g_z_y_y_xxxxyyy[k];

                g_z_y_yy_xxxxyz[k] = g_z_0_y_xxxxyz[k] - g_z_y_y_xxxxyz[k] * ab_y + g_z_y_y_xxxxyyz[k];

                g_z_y_yy_xxxxzz[k] = g_z_0_y_xxxxzz[k] - g_z_y_y_xxxxzz[k] * ab_y + g_z_y_y_xxxxyzz[k];

                g_z_y_yy_xxxyyy[k] = g_z_0_y_xxxyyy[k] - g_z_y_y_xxxyyy[k] * ab_y + g_z_y_y_xxxyyyy[k];

                g_z_y_yy_xxxyyz[k] = g_z_0_y_xxxyyz[k] - g_z_y_y_xxxyyz[k] * ab_y + g_z_y_y_xxxyyyz[k];

                g_z_y_yy_xxxyzz[k] = g_z_0_y_xxxyzz[k] - g_z_y_y_xxxyzz[k] * ab_y + g_z_y_y_xxxyyzz[k];

                g_z_y_yy_xxxzzz[k] = g_z_0_y_xxxzzz[k] - g_z_y_y_xxxzzz[k] * ab_y + g_z_y_y_xxxyzzz[k];

                g_z_y_yy_xxyyyy[k] = g_z_0_y_xxyyyy[k] - g_z_y_y_xxyyyy[k] * ab_y + g_z_y_y_xxyyyyy[k];

                g_z_y_yy_xxyyyz[k] = g_z_0_y_xxyyyz[k] - g_z_y_y_xxyyyz[k] * ab_y + g_z_y_y_xxyyyyz[k];

                g_z_y_yy_xxyyzz[k] = g_z_0_y_xxyyzz[k] - g_z_y_y_xxyyzz[k] * ab_y + g_z_y_y_xxyyyzz[k];

                g_z_y_yy_xxyzzz[k] = g_z_0_y_xxyzzz[k] - g_z_y_y_xxyzzz[k] * ab_y + g_z_y_y_xxyyzzz[k];

                g_z_y_yy_xxzzzz[k] = g_z_0_y_xxzzzz[k] - g_z_y_y_xxzzzz[k] * ab_y + g_z_y_y_xxyzzzz[k];

                g_z_y_yy_xyyyyy[k] = g_z_0_y_xyyyyy[k] - g_z_y_y_xyyyyy[k] * ab_y + g_z_y_y_xyyyyyy[k];

                g_z_y_yy_xyyyyz[k] = g_z_0_y_xyyyyz[k] - g_z_y_y_xyyyyz[k] * ab_y + g_z_y_y_xyyyyyz[k];

                g_z_y_yy_xyyyzz[k] = g_z_0_y_xyyyzz[k] - g_z_y_y_xyyyzz[k] * ab_y + g_z_y_y_xyyyyzz[k];

                g_z_y_yy_xyyzzz[k] = g_z_0_y_xyyzzz[k] - g_z_y_y_xyyzzz[k] * ab_y + g_z_y_y_xyyyzzz[k];

                g_z_y_yy_xyzzzz[k] = g_z_0_y_xyzzzz[k] - g_z_y_y_xyzzzz[k] * ab_y + g_z_y_y_xyyzzzz[k];

                g_z_y_yy_xzzzzz[k] = g_z_0_y_xzzzzz[k] - g_z_y_y_xzzzzz[k] * ab_y + g_z_y_y_xyzzzzz[k];

                g_z_y_yy_yyyyyy[k] = g_z_0_y_yyyyyy[k] - g_z_y_y_yyyyyy[k] * ab_y + g_z_y_y_yyyyyyy[k];

                g_z_y_yy_yyyyyz[k] = g_z_0_y_yyyyyz[k] - g_z_y_y_yyyyyz[k] * ab_y + g_z_y_y_yyyyyyz[k];

                g_z_y_yy_yyyyzz[k] = g_z_0_y_yyyyzz[k] - g_z_y_y_yyyyzz[k] * ab_y + g_z_y_y_yyyyyzz[k];

                g_z_y_yy_yyyzzz[k] = g_z_0_y_yyyzzz[k] - g_z_y_y_yyyzzz[k] * ab_y + g_z_y_y_yyyyzzz[k];

                g_z_y_yy_yyzzzz[k] = g_z_0_y_yyzzzz[k] - g_z_y_y_yyzzzz[k] * ab_y + g_z_y_y_yyyzzzz[k];

                g_z_y_yy_yzzzzz[k] = g_z_0_y_yzzzzz[k] - g_z_y_y_yzzzzz[k] * ab_y + g_z_y_y_yyzzzzz[k];

                g_z_y_yy_zzzzzz[k] = g_z_0_y_zzzzzz[k] - g_z_y_y_zzzzzz[k] * ab_y + g_z_y_y_yzzzzzz[k];
            }

            /// Set up 1288-1316 components of targeted buffer : cbuffer.data(

            auto g_z_y_yz_xxxxxx = cbuffer.data(di_geom_11_off + 1288 * ccomps * dcomps);

            auto g_z_y_yz_xxxxxy = cbuffer.data(di_geom_11_off + 1289 * ccomps * dcomps);

            auto g_z_y_yz_xxxxxz = cbuffer.data(di_geom_11_off + 1290 * ccomps * dcomps);

            auto g_z_y_yz_xxxxyy = cbuffer.data(di_geom_11_off + 1291 * ccomps * dcomps);

            auto g_z_y_yz_xxxxyz = cbuffer.data(di_geom_11_off + 1292 * ccomps * dcomps);

            auto g_z_y_yz_xxxxzz = cbuffer.data(di_geom_11_off + 1293 * ccomps * dcomps);

            auto g_z_y_yz_xxxyyy = cbuffer.data(di_geom_11_off + 1294 * ccomps * dcomps);

            auto g_z_y_yz_xxxyyz = cbuffer.data(di_geom_11_off + 1295 * ccomps * dcomps);

            auto g_z_y_yz_xxxyzz = cbuffer.data(di_geom_11_off + 1296 * ccomps * dcomps);

            auto g_z_y_yz_xxxzzz = cbuffer.data(di_geom_11_off + 1297 * ccomps * dcomps);

            auto g_z_y_yz_xxyyyy = cbuffer.data(di_geom_11_off + 1298 * ccomps * dcomps);

            auto g_z_y_yz_xxyyyz = cbuffer.data(di_geom_11_off + 1299 * ccomps * dcomps);

            auto g_z_y_yz_xxyyzz = cbuffer.data(di_geom_11_off + 1300 * ccomps * dcomps);

            auto g_z_y_yz_xxyzzz = cbuffer.data(di_geom_11_off + 1301 * ccomps * dcomps);

            auto g_z_y_yz_xxzzzz = cbuffer.data(di_geom_11_off + 1302 * ccomps * dcomps);

            auto g_z_y_yz_xyyyyy = cbuffer.data(di_geom_11_off + 1303 * ccomps * dcomps);

            auto g_z_y_yz_xyyyyz = cbuffer.data(di_geom_11_off + 1304 * ccomps * dcomps);

            auto g_z_y_yz_xyyyzz = cbuffer.data(di_geom_11_off + 1305 * ccomps * dcomps);

            auto g_z_y_yz_xyyzzz = cbuffer.data(di_geom_11_off + 1306 * ccomps * dcomps);

            auto g_z_y_yz_xyzzzz = cbuffer.data(di_geom_11_off + 1307 * ccomps * dcomps);

            auto g_z_y_yz_xzzzzz = cbuffer.data(di_geom_11_off + 1308 * ccomps * dcomps);

            auto g_z_y_yz_yyyyyy = cbuffer.data(di_geom_11_off + 1309 * ccomps * dcomps);

            auto g_z_y_yz_yyyyyz = cbuffer.data(di_geom_11_off + 1310 * ccomps * dcomps);

            auto g_z_y_yz_yyyyzz = cbuffer.data(di_geom_11_off + 1311 * ccomps * dcomps);

            auto g_z_y_yz_yyyzzz = cbuffer.data(di_geom_11_off + 1312 * ccomps * dcomps);

            auto g_z_y_yz_yyzzzz = cbuffer.data(di_geom_11_off + 1313 * ccomps * dcomps);

            auto g_z_y_yz_yzzzzz = cbuffer.data(di_geom_11_off + 1314 * ccomps * dcomps);

            auto g_z_y_yz_zzzzzz = cbuffer.data(di_geom_11_off + 1315 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxxxx, g_z_0_z_xxxxxy, g_z_0_z_xxxxxz, g_z_0_z_xxxxyy, g_z_0_z_xxxxyz, g_z_0_z_xxxxzz, g_z_0_z_xxxyyy, g_z_0_z_xxxyyz, g_z_0_z_xxxyzz, g_z_0_z_xxxzzz, g_z_0_z_xxyyyy, g_z_0_z_xxyyyz, g_z_0_z_xxyyzz, g_z_0_z_xxyzzz, g_z_0_z_xxzzzz, g_z_0_z_xyyyyy, g_z_0_z_xyyyyz, g_z_0_z_xyyyzz, g_z_0_z_xyyzzz, g_z_0_z_xyzzzz, g_z_0_z_xzzzzz, g_z_0_z_yyyyyy, g_z_0_z_yyyyyz, g_z_0_z_yyyyzz, g_z_0_z_yyyzzz, g_z_0_z_yyzzzz, g_z_0_z_yzzzzz, g_z_0_z_zzzzzz, g_z_y_yz_xxxxxx, g_z_y_yz_xxxxxy, g_z_y_yz_xxxxxz, g_z_y_yz_xxxxyy, g_z_y_yz_xxxxyz, g_z_y_yz_xxxxzz, g_z_y_yz_xxxyyy, g_z_y_yz_xxxyyz, g_z_y_yz_xxxyzz, g_z_y_yz_xxxzzz, g_z_y_yz_xxyyyy, g_z_y_yz_xxyyyz, g_z_y_yz_xxyyzz, g_z_y_yz_xxyzzz, g_z_y_yz_xxzzzz, g_z_y_yz_xyyyyy, g_z_y_yz_xyyyyz, g_z_y_yz_xyyyzz, g_z_y_yz_xyyzzz, g_z_y_yz_xyzzzz, g_z_y_yz_xzzzzz, g_z_y_yz_yyyyyy, g_z_y_yz_yyyyyz, g_z_y_yz_yyyyzz, g_z_y_yz_yyyzzz, g_z_y_yz_yyzzzz, g_z_y_yz_yzzzzz, g_z_y_yz_zzzzzz, g_z_y_z_xxxxxx, g_z_y_z_xxxxxxy, g_z_y_z_xxxxxy, g_z_y_z_xxxxxyy, g_z_y_z_xxxxxyz, g_z_y_z_xxxxxz, g_z_y_z_xxxxyy, g_z_y_z_xxxxyyy, g_z_y_z_xxxxyyz, g_z_y_z_xxxxyz, g_z_y_z_xxxxyzz, g_z_y_z_xxxxzz, g_z_y_z_xxxyyy, g_z_y_z_xxxyyyy, g_z_y_z_xxxyyyz, g_z_y_z_xxxyyz, g_z_y_z_xxxyyzz, g_z_y_z_xxxyzz, g_z_y_z_xxxyzzz, g_z_y_z_xxxzzz, g_z_y_z_xxyyyy, g_z_y_z_xxyyyyy, g_z_y_z_xxyyyyz, g_z_y_z_xxyyyz, g_z_y_z_xxyyyzz, g_z_y_z_xxyyzz, g_z_y_z_xxyyzzz, g_z_y_z_xxyzzz, g_z_y_z_xxyzzzz, g_z_y_z_xxzzzz, g_z_y_z_xyyyyy, g_z_y_z_xyyyyyy, g_z_y_z_xyyyyyz, g_z_y_z_xyyyyz, g_z_y_z_xyyyyzz, g_z_y_z_xyyyzz, g_z_y_z_xyyyzzz, g_z_y_z_xyyzzz, g_z_y_z_xyyzzzz, g_z_y_z_xyzzzz, g_z_y_z_xyzzzzz, g_z_y_z_xzzzzz, g_z_y_z_yyyyyy, g_z_y_z_yyyyyyy, g_z_y_z_yyyyyyz, g_z_y_z_yyyyyz, g_z_y_z_yyyyyzz, g_z_y_z_yyyyzz, g_z_y_z_yyyyzzz, g_z_y_z_yyyzzz, g_z_y_z_yyyzzzz, g_z_y_z_yyzzzz, g_z_y_z_yyzzzzz, g_z_y_z_yzzzzz, g_z_y_z_yzzzzzz, g_z_y_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yz_xxxxxx[k] = g_z_0_z_xxxxxx[k] - g_z_y_z_xxxxxx[k] * ab_y + g_z_y_z_xxxxxxy[k];

                g_z_y_yz_xxxxxy[k] = g_z_0_z_xxxxxy[k] - g_z_y_z_xxxxxy[k] * ab_y + g_z_y_z_xxxxxyy[k];

                g_z_y_yz_xxxxxz[k] = g_z_0_z_xxxxxz[k] - g_z_y_z_xxxxxz[k] * ab_y + g_z_y_z_xxxxxyz[k];

                g_z_y_yz_xxxxyy[k] = g_z_0_z_xxxxyy[k] - g_z_y_z_xxxxyy[k] * ab_y + g_z_y_z_xxxxyyy[k];

                g_z_y_yz_xxxxyz[k] = g_z_0_z_xxxxyz[k] - g_z_y_z_xxxxyz[k] * ab_y + g_z_y_z_xxxxyyz[k];

                g_z_y_yz_xxxxzz[k] = g_z_0_z_xxxxzz[k] - g_z_y_z_xxxxzz[k] * ab_y + g_z_y_z_xxxxyzz[k];

                g_z_y_yz_xxxyyy[k] = g_z_0_z_xxxyyy[k] - g_z_y_z_xxxyyy[k] * ab_y + g_z_y_z_xxxyyyy[k];

                g_z_y_yz_xxxyyz[k] = g_z_0_z_xxxyyz[k] - g_z_y_z_xxxyyz[k] * ab_y + g_z_y_z_xxxyyyz[k];

                g_z_y_yz_xxxyzz[k] = g_z_0_z_xxxyzz[k] - g_z_y_z_xxxyzz[k] * ab_y + g_z_y_z_xxxyyzz[k];

                g_z_y_yz_xxxzzz[k] = g_z_0_z_xxxzzz[k] - g_z_y_z_xxxzzz[k] * ab_y + g_z_y_z_xxxyzzz[k];

                g_z_y_yz_xxyyyy[k] = g_z_0_z_xxyyyy[k] - g_z_y_z_xxyyyy[k] * ab_y + g_z_y_z_xxyyyyy[k];

                g_z_y_yz_xxyyyz[k] = g_z_0_z_xxyyyz[k] - g_z_y_z_xxyyyz[k] * ab_y + g_z_y_z_xxyyyyz[k];

                g_z_y_yz_xxyyzz[k] = g_z_0_z_xxyyzz[k] - g_z_y_z_xxyyzz[k] * ab_y + g_z_y_z_xxyyyzz[k];

                g_z_y_yz_xxyzzz[k] = g_z_0_z_xxyzzz[k] - g_z_y_z_xxyzzz[k] * ab_y + g_z_y_z_xxyyzzz[k];

                g_z_y_yz_xxzzzz[k] = g_z_0_z_xxzzzz[k] - g_z_y_z_xxzzzz[k] * ab_y + g_z_y_z_xxyzzzz[k];

                g_z_y_yz_xyyyyy[k] = g_z_0_z_xyyyyy[k] - g_z_y_z_xyyyyy[k] * ab_y + g_z_y_z_xyyyyyy[k];

                g_z_y_yz_xyyyyz[k] = g_z_0_z_xyyyyz[k] - g_z_y_z_xyyyyz[k] * ab_y + g_z_y_z_xyyyyyz[k];

                g_z_y_yz_xyyyzz[k] = g_z_0_z_xyyyzz[k] - g_z_y_z_xyyyzz[k] * ab_y + g_z_y_z_xyyyyzz[k];

                g_z_y_yz_xyyzzz[k] = g_z_0_z_xyyzzz[k] - g_z_y_z_xyyzzz[k] * ab_y + g_z_y_z_xyyyzzz[k];

                g_z_y_yz_xyzzzz[k] = g_z_0_z_xyzzzz[k] - g_z_y_z_xyzzzz[k] * ab_y + g_z_y_z_xyyzzzz[k];

                g_z_y_yz_xzzzzz[k] = g_z_0_z_xzzzzz[k] - g_z_y_z_xzzzzz[k] * ab_y + g_z_y_z_xyzzzzz[k];

                g_z_y_yz_yyyyyy[k] = g_z_0_z_yyyyyy[k] - g_z_y_z_yyyyyy[k] * ab_y + g_z_y_z_yyyyyyy[k];

                g_z_y_yz_yyyyyz[k] = g_z_0_z_yyyyyz[k] - g_z_y_z_yyyyyz[k] * ab_y + g_z_y_z_yyyyyyz[k];

                g_z_y_yz_yyyyzz[k] = g_z_0_z_yyyyzz[k] - g_z_y_z_yyyyzz[k] * ab_y + g_z_y_z_yyyyyzz[k];

                g_z_y_yz_yyyzzz[k] = g_z_0_z_yyyzzz[k] - g_z_y_z_yyyzzz[k] * ab_y + g_z_y_z_yyyyzzz[k];

                g_z_y_yz_yyzzzz[k] = g_z_0_z_yyzzzz[k] - g_z_y_z_yyzzzz[k] * ab_y + g_z_y_z_yyyzzzz[k];

                g_z_y_yz_yzzzzz[k] = g_z_0_z_yzzzzz[k] - g_z_y_z_yzzzzz[k] * ab_y + g_z_y_z_yyzzzzz[k];

                g_z_y_yz_zzzzzz[k] = g_z_0_z_zzzzzz[k] - g_z_y_z_zzzzzz[k] * ab_y + g_z_y_z_yzzzzzz[k];
            }

            /// Set up 1316-1344 components of targeted buffer : cbuffer.data(

            auto g_z_y_zz_xxxxxx = cbuffer.data(di_geom_11_off + 1316 * ccomps * dcomps);

            auto g_z_y_zz_xxxxxy = cbuffer.data(di_geom_11_off + 1317 * ccomps * dcomps);

            auto g_z_y_zz_xxxxxz = cbuffer.data(di_geom_11_off + 1318 * ccomps * dcomps);

            auto g_z_y_zz_xxxxyy = cbuffer.data(di_geom_11_off + 1319 * ccomps * dcomps);

            auto g_z_y_zz_xxxxyz = cbuffer.data(di_geom_11_off + 1320 * ccomps * dcomps);

            auto g_z_y_zz_xxxxzz = cbuffer.data(di_geom_11_off + 1321 * ccomps * dcomps);

            auto g_z_y_zz_xxxyyy = cbuffer.data(di_geom_11_off + 1322 * ccomps * dcomps);

            auto g_z_y_zz_xxxyyz = cbuffer.data(di_geom_11_off + 1323 * ccomps * dcomps);

            auto g_z_y_zz_xxxyzz = cbuffer.data(di_geom_11_off + 1324 * ccomps * dcomps);

            auto g_z_y_zz_xxxzzz = cbuffer.data(di_geom_11_off + 1325 * ccomps * dcomps);

            auto g_z_y_zz_xxyyyy = cbuffer.data(di_geom_11_off + 1326 * ccomps * dcomps);

            auto g_z_y_zz_xxyyyz = cbuffer.data(di_geom_11_off + 1327 * ccomps * dcomps);

            auto g_z_y_zz_xxyyzz = cbuffer.data(di_geom_11_off + 1328 * ccomps * dcomps);

            auto g_z_y_zz_xxyzzz = cbuffer.data(di_geom_11_off + 1329 * ccomps * dcomps);

            auto g_z_y_zz_xxzzzz = cbuffer.data(di_geom_11_off + 1330 * ccomps * dcomps);

            auto g_z_y_zz_xyyyyy = cbuffer.data(di_geom_11_off + 1331 * ccomps * dcomps);

            auto g_z_y_zz_xyyyyz = cbuffer.data(di_geom_11_off + 1332 * ccomps * dcomps);

            auto g_z_y_zz_xyyyzz = cbuffer.data(di_geom_11_off + 1333 * ccomps * dcomps);

            auto g_z_y_zz_xyyzzz = cbuffer.data(di_geom_11_off + 1334 * ccomps * dcomps);

            auto g_z_y_zz_xyzzzz = cbuffer.data(di_geom_11_off + 1335 * ccomps * dcomps);

            auto g_z_y_zz_xzzzzz = cbuffer.data(di_geom_11_off + 1336 * ccomps * dcomps);

            auto g_z_y_zz_yyyyyy = cbuffer.data(di_geom_11_off + 1337 * ccomps * dcomps);

            auto g_z_y_zz_yyyyyz = cbuffer.data(di_geom_11_off + 1338 * ccomps * dcomps);

            auto g_z_y_zz_yyyyzz = cbuffer.data(di_geom_11_off + 1339 * ccomps * dcomps);

            auto g_z_y_zz_yyyzzz = cbuffer.data(di_geom_11_off + 1340 * ccomps * dcomps);

            auto g_z_y_zz_yyzzzz = cbuffer.data(di_geom_11_off + 1341 * ccomps * dcomps);

            auto g_z_y_zz_yzzzzz = cbuffer.data(di_geom_11_off + 1342 * ccomps * dcomps);

            auto g_z_y_zz_zzzzzz = cbuffer.data(di_geom_11_off + 1343 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_z_xxxxxx, g_0_y_z_xxxxxy, g_0_y_z_xxxxxz, g_0_y_z_xxxxyy, g_0_y_z_xxxxyz, g_0_y_z_xxxxzz, g_0_y_z_xxxyyy, g_0_y_z_xxxyyz, g_0_y_z_xxxyzz, g_0_y_z_xxxzzz, g_0_y_z_xxyyyy, g_0_y_z_xxyyyz, g_0_y_z_xxyyzz, g_0_y_z_xxyzzz, g_0_y_z_xxzzzz, g_0_y_z_xyyyyy, g_0_y_z_xyyyyz, g_0_y_z_xyyyzz, g_0_y_z_xyyzzz, g_0_y_z_xyzzzz, g_0_y_z_xzzzzz, g_0_y_z_yyyyyy, g_0_y_z_yyyyyz, g_0_y_z_yyyyzz, g_0_y_z_yyyzzz, g_0_y_z_yyzzzz, g_0_y_z_yzzzzz, g_0_y_z_zzzzzz, g_z_y_z_xxxxxx, g_z_y_z_xxxxxxz, g_z_y_z_xxxxxy, g_z_y_z_xxxxxyz, g_z_y_z_xxxxxz, g_z_y_z_xxxxxzz, g_z_y_z_xxxxyy, g_z_y_z_xxxxyyz, g_z_y_z_xxxxyz, g_z_y_z_xxxxyzz, g_z_y_z_xxxxzz, g_z_y_z_xxxxzzz, g_z_y_z_xxxyyy, g_z_y_z_xxxyyyz, g_z_y_z_xxxyyz, g_z_y_z_xxxyyzz, g_z_y_z_xxxyzz, g_z_y_z_xxxyzzz, g_z_y_z_xxxzzz, g_z_y_z_xxxzzzz, g_z_y_z_xxyyyy, g_z_y_z_xxyyyyz, g_z_y_z_xxyyyz, g_z_y_z_xxyyyzz, g_z_y_z_xxyyzz, g_z_y_z_xxyyzzz, g_z_y_z_xxyzzz, g_z_y_z_xxyzzzz, g_z_y_z_xxzzzz, g_z_y_z_xxzzzzz, g_z_y_z_xyyyyy, g_z_y_z_xyyyyyz, g_z_y_z_xyyyyz, g_z_y_z_xyyyyzz, g_z_y_z_xyyyzz, g_z_y_z_xyyyzzz, g_z_y_z_xyyzzz, g_z_y_z_xyyzzzz, g_z_y_z_xyzzzz, g_z_y_z_xyzzzzz, g_z_y_z_xzzzzz, g_z_y_z_xzzzzzz, g_z_y_z_yyyyyy, g_z_y_z_yyyyyyz, g_z_y_z_yyyyyz, g_z_y_z_yyyyyzz, g_z_y_z_yyyyzz, g_z_y_z_yyyyzzz, g_z_y_z_yyyzzz, g_z_y_z_yyyzzzz, g_z_y_z_yyzzzz, g_z_y_z_yyzzzzz, g_z_y_z_yzzzzz, g_z_y_z_yzzzzzz, g_z_y_z_zzzzzz, g_z_y_z_zzzzzzz, g_z_y_zz_xxxxxx, g_z_y_zz_xxxxxy, g_z_y_zz_xxxxxz, g_z_y_zz_xxxxyy, g_z_y_zz_xxxxyz, g_z_y_zz_xxxxzz, g_z_y_zz_xxxyyy, g_z_y_zz_xxxyyz, g_z_y_zz_xxxyzz, g_z_y_zz_xxxzzz, g_z_y_zz_xxyyyy, g_z_y_zz_xxyyyz, g_z_y_zz_xxyyzz, g_z_y_zz_xxyzzz, g_z_y_zz_xxzzzz, g_z_y_zz_xyyyyy, g_z_y_zz_xyyyyz, g_z_y_zz_xyyyzz, g_z_y_zz_xyyzzz, g_z_y_zz_xyzzzz, g_z_y_zz_xzzzzz, g_z_y_zz_yyyyyy, g_z_y_zz_yyyyyz, g_z_y_zz_yyyyzz, g_z_y_zz_yyyzzz, g_z_y_zz_yyzzzz, g_z_y_zz_yzzzzz, g_z_y_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zz_xxxxxx[k] = -g_0_y_z_xxxxxx[k] - g_z_y_z_xxxxxx[k] * ab_z + g_z_y_z_xxxxxxz[k];

                g_z_y_zz_xxxxxy[k] = -g_0_y_z_xxxxxy[k] - g_z_y_z_xxxxxy[k] * ab_z + g_z_y_z_xxxxxyz[k];

                g_z_y_zz_xxxxxz[k] = -g_0_y_z_xxxxxz[k] - g_z_y_z_xxxxxz[k] * ab_z + g_z_y_z_xxxxxzz[k];

                g_z_y_zz_xxxxyy[k] = -g_0_y_z_xxxxyy[k] - g_z_y_z_xxxxyy[k] * ab_z + g_z_y_z_xxxxyyz[k];

                g_z_y_zz_xxxxyz[k] = -g_0_y_z_xxxxyz[k] - g_z_y_z_xxxxyz[k] * ab_z + g_z_y_z_xxxxyzz[k];

                g_z_y_zz_xxxxzz[k] = -g_0_y_z_xxxxzz[k] - g_z_y_z_xxxxzz[k] * ab_z + g_z_y_z_xxxxzzz[k];

                g_z_y_zz_xxxyyy[k] = -g_0_y_z_xxxyyy[k] - g_z_y_z_xxxyyy[k] * ab_z + g_z_y_z_xxxyyyz[k];

                g_z_y_zz_xxxyyz[k] = -g_0_y_z_xxxyyz[k] - g_z_y_z_xxxyyz[k] * ab_z + g_z_y_z_xxxyyzz[k];

                g_z_y_zz_xxxyzz[k] = -g_0_y_z_xxxyzz[k] - g_z_y_z_xxxyzz[k] * ab_z + g_z_y_z_xxxyzzz[k];

                g_z_y_zz_xxxzzz[k] = -g_0_y_z_xxxzzz[k] - g_z_y_z_xxxzzz[k] * ab_z + g_z_y_z_xxxzzzz[k];

                g_z_y_zz_xxyyyy[k] = -g_0_y_z_xxyyyy[k] - g_z_y_z_xxyyyy[k] * ab_z + g_z_y_z_xxyyyyz[k];

                g_z_y_zz_xxyyyz[k] = -g_0_y_z_xxyyyz[k] - g_z_y_z_xxyyyz[k] * ab_z + g_z_y_z_xxyyyzz[k];

                g_z_y_zz_xxyyzz[k] = -g_0_y_z_xxyyzz[k] - g_z_y_z_xxyyzz[k] * ab_z + g_z_y_z_xxyyzzz[k];

                g_z_y_zz_xxyzzz[k] = -g_0_y_z_xxyzzz[k] - g_z_y_z_xxyzzz[k] * ab_z + g_z_y_z_xxyzzzz[k];

                g_z_y_zz_xxzzzz[k] = -g_0_y_z_xxzzzz[k] - g_z_y_z_xxzzzz[k] * ab_z + g_z_y_z_xxzzzzz[k];

                g_z_y_zz_xyyyyy[k] = -g_0_y_z_xyyyyy[k] - g_z_y_z_xyyyyy[k] * ab_z + g_z_y_z_xyyyyyz[k];

                g_z_y_zz_xyyyyz[k] = -g_0_y_z_xyyyyz[k] - g_z_y_z_xyyyyz[k] * ab_z + g_z_y_z_xyyyyzz[k];

                g_z_y_zz_xyyyzz[k] = -g_0_y_z_xyyyzz[k] - g_z_y_z_xyyyzz[k] * ab_z + g_z_y_z_xyyyzzz[k];

                g_z_y_zz_xyyzzz[k] = -g_0_y_z_xyyzzz[k] - g_z_y_z_xyyzzz[k] * ab_z + g_z_y_z_xyyzzzz[k];

                g_z_y_zz_xyzzzz[k] = -g_0_y_z_xyzzzz[k] - g_z_y_z_xyzzzz[k] * ab_z + g_z_y_z_xyzzzzz[k];

                g_z_y_zz_xzzzzz[k] = -g_0_y_z_xzzzzz[k] - g_z_y_z_xzzzzz[k] * ab_z + g_z_y_z_xzzzzzz[k];

                g_z_y_zz_yyyyyy[k] = -g_0_y_z_yyyyyy[k] - g_z_y_z_yyyyyy[k] * ab_z + g_z_y_z_yyyyyyz[k];

                g_z_y_zz_yyyyyz[k] = -g_0_y_z_yyyyyz[k] - g_z_y_z_yyyyyz[k] * ab_z + g_z_y_z_yyyyyzz[k];

                g_z_y_zz_yyyyzz[k] = -g_0_y_z_yyyyzz[k] - g_z_y_z_yyyyzz[k] * ab_z + g_z_y_z_yyyyzzz[k];

                g_z_y_zz_yyyzzz[k] = -g_0_y_z_yyyzzz[k] - g_z_y_z_yyyzzz[k] * ab_z + g_z_y_z_yyyzzzz[k];

                g_z_y_zz_yyzzzz[k] = -g_0_y_z_yyzzzz[k] - g_z_y_z_yyzzzz[k] * ab_z + g_z_y_z_yyzzzzz[k];

                g_z_y_zz_yzzzzz[k] = -g_0_y_z_yzzzzz[k] - g_z_y_z_yzzzzz[k] * ab_z + g_z_y_z_yzzzzzz[k];

                g_z_y_zz_zzzzzz[k] = -g_0_y_z_zzzzzz[k] - g_z_y_z_zzzzzz[k] * ab_z + g_z_y_z_zzzzzzz[k];
            }

            /// Set up 1344-1372 components of targeted buffer : cbuffer.data(

            auto g_z_z_xx_xxxxxx = cbuffer.data(di_geom_11_off + 1344 * ccomps * dcomps);

            auto g_z_z_xx_xxxxxy = cbuffer.data(di_geom_11_off + 1345 * ccomps * dcomps);

            auto g_z_z_xx_xxxxxz = cbuffer.data(di_geom_11_off + 1346 * ccomps * dcomps);

            auto g_z_z_xx_xxxxyy = cbuffer.data(di_geom_11_off + 1347 * ccomps * dcomps);

            auto g_z_z_xx_xxxxyz = cbuffer.data(di_geom_11_off + 1348 * ccomps * dcomps);

            auto g_z_z_xx_xxxxzz = cbuffer.data(di_geom_11_off + 1349 * ccomps * dcomps);

            auto g_z_z_xx_xxxyyy = cbuffer.data(di_geom_11_off + 1350 * ccomps * dcomps);

            auto g_z_z_xx_xxxyyz = cbuffer.data(di_geom_11_off + 1351 * ccomps * dcomps);

            auto g_z_z_xx_xxxyzz = cbuffer.data(di_geom_11_off + 1352 * ccomps * dcomps);

            auto g_z_z_xx_xxxzzz = cbuffer.data(di_geom_11_off + 1353 * ccomps * dcomps);

            auto g_z_z_xx_xxyyyy = cbuffer.data(di_geom_11_off + 1354 * ccomps * dcomps);

            auto g_z_z_xx_xxyyyz = cbuffer.data(di_geom_11_off + 1355 * ccomps * dcomps);

            auto g_z_z_xx_xxyyzz = cbuffer.data(di_geom_11_off + 1356 * ccomps * dcomps);

            auto g_z_z_xx_xxyzzz = cbuffer.data(di_geom_11_off + 1357 * ccomps * dcomps);

            auto g_z_z_xx_xxzzzz = cbuffer.data(di_geom_11_off + 1358 * ccomps * dcomps);

            auto g_z_z_xx_xyyyyy = cbuffer.data(di_geom_11_off + 1359 * ccomps * dcomps);

            auto g_z_z_xx_xyyyyz = cbuffer.data(di_geom_11_off + 1360 * ccomps * dcomps);

            auto g_z_z_xx_xyyyzz = cbuffer.data(di_geom_11_off + 1361 * ccomps * dcomps);

            auto g_z_z_xx_xyyzzz = cbuffer.data(di_geom_11_off + 1362 * ccomps * dcomps);

            auto g_z_z_xx_xyzzzz = cbuffer.data(di_geom_11_off + 1363 * ccomps * dcomps);

            auto g_z_z_xx_xzzzzz = cbuffer.data(di_geom_11_off + 1364 * ccomps * dcomps);

            auto g_z_z_xx_yyyyyy = cbuffer.data(di_geom_11_off + 1365 * ccomps * dcomps);

            auto g_z_z_xx_yyyyyz = cbuffer.data(di_geom_11_off + 1366 * ccomps * dcomps);

            auto g_z_z_xx_yyyyzz = cbuffer.data(di_geom_11_off + 1367 * ccomps * dcomps);

            auto g_z_z_xx_yyyzzz = cbuffer.data(di_geom_11_off + 1368 * ccomps * dcomps);

            auto g_z_z_xx_yyzzzz = cbuffer.data(di_geom_11_off + 1369 * ccomps * dcomps);

            auto g_z_z_xx_yzzzzz = cbuffer.data(di_geom_11_off + 1370 * ccomps * dcomps);

            auto g_z_z_xx_zzzzzz = cbuffer.data(di_geom_11_off + 1371 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_x_xxxxxx, g_z_z_x_xxxxxxx, g_z_z_x_xxxxxxy, g_z_z_x_xxxxxxz, g_z_z_x_xxxxxy, g_z_z_x_xxxxxyy, g_z_z_x_xxxxxyz, g_z_z_x_xxxxxz, g_z_z_x_xxxxxzz, g_z_z_x_xxxxyy, g_z_z_x_xxxxyyy, g_z_z_x_xxxxyyz, g_z_z_x_xxxxyz, g_z_z_x_xxxxyzz, g_z_z_x_xxxxzz, g_z_z_x_xxxxzzz, g_z_z_x_xxxyyy, g_z_z_x_xxxyyyy, g_z_z_x_xxxyyyz, g_z_z_x_xxxyyz, g_z_z_x_xxxyyzz, g_z_z_x_xxxyzz, g_z_z_x_xxxyzzz, g_z_z_x_xxxzzz, g_z_z_x_xxxzzzz, g_z_z_x_xxyyyy, g_z_z_x_xxyyyyy, g_z_z_x_xxyyyyz, g_z_z_x_xxyyyz, g_z_z_x_xxyyyzz, g_z_z_x_xxyyzz, g_z_z_x_xxyyzzz, g_z_z_x_xxyzzz, g_z_z_x_xxyzzzz, g_z_z_x_xxzzzz, g_z_z_x_xxzzzzz, g_z_z_x_xyyyyy, g_z_z_x_xyyyyyy, g_z_z_x_xyyyyyz, g_z_z_x_xyyyyz, g_z_z_x_xyyyyzz, g_z_z_x_xyyyzz, g_z_z_x_xyyyzzz, g_z_z_x_xyyzzz, g_z_z_x_xyyzzzz, g_z_z_x_xyzzzz, g_z_z_x_xyzzzzz, g_z_z_x_xzzzzz, g_z_z_x_xzzzzzz, g_z_z_x_yyyyyy, g_z_z_x_yyyyyz, g_z_z_x_yyyyzz, g_z_z_x_yyyzzz, g_z_z_x_yyzzzz, g_z_z_x_yzzzzz, g_z_z_x_zzzzzz, g_z_z_xx_xxxxxx, g_z_z_xx_xxxxxy, g_z_z_xx_xxxxxz, g_z_z_xx_xxxxyy, g_z_z_xx_xxxxyz, g_z_z_xx_xxxxzz, g_z_z_xx_xxxyyy, g_z_z_xx_xxxyyz, g_z_z_xx_xxxyzz, g_z_z_xx_xxxzzz, g_z_z_xx_xxyyyy, g_z_z_xx_xxyyyz, g_z_z_xx_xxyyzz, g_z_z_xx_xxyzzz, g_z_z_xx_xxzzzz, g_z_z_xx_xyyyyy, g_z_z_xx_xyyyyz, g_z_z_xx_xyyyzz, g_z_z_xx_xyyzzz, g_z_z_xx_xyzzzz, g_z_z_xx_xzzzzz, g_z_z_xx_yyyyyy, g_z_z_xx_yyyyyz, g_z_z_xx_yyyyzz, g_z_z_xx_yyyzzz, g_z_z_xx_yyzzzz, g_z_z_xx_yzzzzz, g_z_z_xx_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xx_xxxxxx[k] = -g_z_z_x_xxxxxx[k] * ab_x + g_z_z_x_xxxxxxx[k];

                g_z_z_xx_xxxxxy[k] = -g_z_z_x_xxxxxy[k] * ab_x + g_z_z_x_xxxxxxy[k];

                g_z_z_xx_xxxxxz[k] = -g_z_z_x_xxxxxz[k] * ab_x + g_z_z_x_xxxxxxz[k];

                g_z_z_xx_xxxxyy[k] = -g_z_z_x_xxxxyy[k] * ab_x + g_z_z_x_xxxxxyy[k];

                g_z_z_xx_xxxxyz[k] = -g_z_z_x_xxxxyz[k] * ab_x + g_z_z_x_xxxxxyz[k];

                g_z_z_xx_xxxxzz[k] = -g_z_z_x_xxxxzz[k] * ab_x + g_z_z_x_xxxxxzz[k];

                g_z_z_xx_xxxyyy[k] = -g_z_z_x_xxxyyy[k] * ab_x + g_z_z_x_xxxxyyy[k];

                g_z_z_xx_xxxyyz[k] = -g_z_z_x_xxxyyz[k] * ab_x + g_z_z_x_xxxxyyz[k];

                g_z_z_xx_xxxyzz[k] = -g_z_z_x_xxxyzz[k] * ab_x + g_z_z_x_xxxxyzz[k];

                g_z_z_xx_xxxzzz[k] = -g_z_z_x_xxxzzz[k] * ab_x + g_z_z_x_xxxxzzz[k];

                g_z_z_xx_xxyyyy[k] = -g_z_z_x_xxyyyy[k] * ab_x + g_z_z_x_xxxyyyy[k];

                g_z_z_xx_xxyyyz[k] = -g_z_z_x_xxyyyz[k] * ab_x + g_z_z_x_xxxyyyz[k];

                g_z_z_xx_xxyyzz[k] = -g_z_z_x_xxyyzz[k] * ab_x + g_z_z_x_xxxyyzz[k];

                g_z_z_xx_xxyzzz[k] = -g_z_z_x_xxyzzz[k] * ab_x + g_z_z_x_xxxyzzz[k];

                g_z_z_xx_xxzzzz[k] = -g_z_z_x_xxzzzz[k] * ab_x + g_z_z_x_xxxzzzz[k];

                g_z_z_xx_xyyyyy[k] = -g_z_z_x_xyyyyy[k] * ab_x + g_z_z_x_xxyyyyy[k];

                g_z_z_xx_xyyyyz[k] = -g_z_z_x_xyyyyz[k] * ab_x + g_z_z_x_xxyyyyz[k];

                g_z_z_xx_xyyyzz[k] = -g_z_z_x_xyyyzz[k] * ab_x + g_z_z_x_xxyyyzz[k];

                g_z_z_xx_xyyzzz[k] = -g_z_z_x_xyyzzz[k] * ab_x + g_z_z_x_xxyyzzz[k];

                g_z_z_xx_xyzzzz[k] = -g_z_z_x_xyzzzz[k] * ab_x + g_z_z_x_xxyzzzz[k];

                g_z_z_xx_xzzzzz[k] = -g_z_z_x_xzzzzz[k] * ab_x + g_z_z_x_xxzzzzz[k];

                g_z_z_xx_yyyyyy[k] = -g_z_z_x_yyyyyy[k] * ab_x + g_z_z_x_xyyyyyy[k];

                g_z_z_xx_yyyyyz[k] = -g_z_z_x_yyyyyz[k] * ab_x + g_z_z_x_xyyyyyz[k];

                g_z_z_xx_yyyyzz[k] = -g_z_z_x_yyyyzz[k] * ab_x + g_z_z_x_xyyyyzz[k];

                g_z_z_xx_yyyzzz[k] = -g_z_z_x_yyyzzz[k] * ab_x + g_z_z_x_xyyyzzz[k];

                g_z_z_xx_yyzzzz[k] = -g_z_z_x_yyzzzz[k] * ab_x + g_z_z_x_xyyzzzz[k];

                g_z_z_xx_yzzzzz[k] = -g_z_z_x_yzzzzz[k] * ab_x + g_z_z_x_xyzzzzz[k];

                g_z_z_xx_zzzzzz[k] = -g_z_z_x_zzzzzz[k] * ab_x + g_z_z_x_xzzzzzz[k];
            }

            /// Set up 1372-1400 components of targeted buffer : cbuffer.data(

            auto g_z_z_xy_xxxxxx = cbuffer.data(di_geom_11_off + 1372 * ccomps * dcomps);

            auto g_z_z_xy_xxxxxy = cbuffer.data(di_geom_11_off + 1373 * ccomps * dcomps);

            auto g_z_z_xy_xxxxxz = cbuffer.data(di_geom_11_off + 1374 * ccomps * dcomps);

            auto g_z_z_xy_xxxxyy = cbuffer.data(di_geom_11_off + 1375 * ccomps * dcomps);

            auto g_z_z_xy_xxxxyz = cbuffer.data(di_geom_11_off + 1376 * ccomps * dcomps);

            auto g_z_z_xy_xxxxzz = cbuffer.data(di_geom_11_off + 1377 * ccomps * dcomps);

            auto g_z_z_xy_xxxyyy = cbuffer.data(di_geom_11_off + 1378 * ccomps * dcomps);

            auto g_z_z_xy_xxxyyz = cbuffer.data(di_geom_11_off + 1379 * ccomps * dcomps);

            auto g_z_z_xy_xxxyzz = cbuffer.data(di_geom_11_off + 1380 * ccomps * dcomps);

            auto g_z_z_xy_xxxzzz = cbuffer.data(di_geom_11_off + 1381 * ccomps * dcomps);

            auto g_z_z_xy_xxyyyy = cbuffer.data(di_geom_11_off + 1382 * ccomps * dcomps);

            auto g_z_z_xy_xxyyyz = cbuffer.data(di_geom_11_off + 1383 * ccomps * dcomps);

            auto g_z_z_xy_xxyyzz = cbuffer.data(di_geom_11_off + 1384 * ccomps * dcomps);

            auto g_z_z_xy_xxyzzz = cbuffer.data(di_geom_11_off + 1385 * ccomps * dcomps);

            auto g_z_z_xy_xxzzzz = cbuffer.data(di_geom_11_off + 1386 * ccomps * dcomps);

            auto g_z_z_xy_xyyyyy = cbuffer.data(di_geom_11_off + 1387 * ccomps * dcomps);

            auto g_z_z_xy_xyyyyz = cbuffer.data(di_geom_11_off + 1388 * ccomps * dcomps);

            auto g_z_z_xy_xyyyzz = cbuffer.data(di_geom_11_off + 1389 * ccomps * dcomps);

            auto g_z_z_xy_xyyzzz = cbuffer.data(di_geom_11_off + 1390 * ccomps * dcomps);

            auto g_z_z_xy_xyzzzz = cbuffer.data(di_geom_11_off + 1391 * ccomps * dcomps);

            auto g_z_z_xy_xzzzzz = cbuffer.data(di_geom_11_off + 1392 * ccomps * dcomps);

            auto g_z_z_xy_yyyyyy = cbuffer.data(di_geom_11_off + 1393 * ccomps * dcomps);

            auto g_z_z_xy_yyyyyz = cbuffer.data(di_geom_11_off + 1394 * ccomps * dcomps);

            auto g_z_z_xy_yyyyzz = cbuffer.data(di_geom_11_off + 1395 * ccomps * dcomps);

            auto g_z_z_xy_yyyzzz = cbuffer.data(di_geom_11_off + 1396 * ccomps * dcomps);

            auto g_z_z_xy_yyzzzz = cbuffer.data(di_geom_11_off + 1397 * ccomps * dcomps);

            auto g_z_z_xy_yzzzzz = cbuffer.data(di_geom_11_off + 1398 * ccomps * dcomps);

            auto g_z_z_xy_zzzzzz = cbuffer.data(di_geom_11_off + 1399 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xy_xxxxxx, g_z_z_xy_xxxxxy, g_z_z_xy_xxxxxz, g_z_z_xy_xxxxyy, g_z_z_xy_xxxxyz, g_z_z_xy_xxxxzz, g_z_z_xy_xxxyyy, g_z_z_xy_xxxyyz, g_z_z_xy_xxxyzz, g_z_z_xy_xxxzzz, g_z_z_xy_xxyyyy, g_z_z_xy_xxyyyz, g_z_z_xy_xxyyzz, g_z_z_xy_xxyzzz, g_z_z_xy_xxzzzz, g_z_z_xy_xyyyyy, g_z_z_xy_xyyyyz, g_z_z_xy_xyyyzz, g_z_z_xy_xyyzzz, g_z_z_xy_xyzzzz, g_z_z_xy_xzzzzz, g_z_z_xy_yyyyyy, g_z_z_xy_yyyyyz, g_z_z_xy_yyyyzz, g_z_z_xy_yyyzzz, g_z_z_xy_yyzzzz, g_z_z_xy_yzzzzz, g_z_z_xy_zzzzzz, g_z_z_y_xxxxxx, g_z_z_y_xxxxxxx, g_z_z_y_xxxxxxy, g_z_z_y_xxxxxxz, g_z_z_y_xxxxxy, g_z_z_y_xxxxxyy, g_z_z_y_xxxxxyz, g_z_z_y_xxxxxz, g_z_z_y_xxxxxzz, g_z_z_y_xxxxyy, g_z_z_y_xxxxyyy, g_z_z_y_xxxxyyz, g_z_z_y_xxxxyz, g_z_z_y_xxxxyzz, g_z_z_y_xxxxzz, g_z_z_y_xxxxzzz, g_z_z_y_xxxyyy, g_z_z_y_xxxyyyy, g_z_z_y_xxxyyyz, g_z_z_y_xxxyyz, g_z_z_y_xxxyyzz, g_z_z_y_xxxyzz, g_z_z_y_xxxyzzz, g_z_z_y_xxxzzz, g_z_z_y_xxxzzzz, g_z_z_y_xxyyyy, g_z_z_y_xxyyyyy, g_z_z_y_xxyyyyz, g_z_z_y_xxyyyz, g_z_z_y_xxyyyzz, g_z_z_y_xxyyzz, g_z_z_y_xxyyzzz, g_z_z_y_xxyzzz, g_z_z_y_xxyzzzz, g_z_z_y_xxzzzz, g_z_z_y_xxzzzzz, g_z_z_y_xyyyyy, g_z_z_y_xyyyyyy, g_z_z_y_xyyyyyz, g_z_z_y_xyyyyz, g_z_z_y_xyyyyzz, g_z_z_y_xyyyzz, g_z_z_y_xyyyzzz, g_z_z_y_xyyzzz, g_z_z_y_xyyzzzz, g_z_z_y_xyzzzz, g_z_z_y_xyzzzzz, g_z_z_y_xzzzzz, g_z_z_y_xzzzzzz, g_z_z_y_yyyyyy, g_z_z_y_yyyyyz, g_z_z_y_yyyyzz, g_z_z_y_yyyzzz, g_z_z_y_yyzzzz, g_z_z_y_yzzzzz, g_z_z_y_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xy_xxxxxx[k] = -g_z_z_y_xxxxxx[k] * ab_x + g_z_z_y_xxxxxxx[k];

                g_z_z_xy_xxxxxy[k] = -g_z_z_y_xxxxxy[k] * ab_x + g_z_z_y_xxxxxxy[k];

                g_z_z_xy_xxxxxz[k] = -g_z_z_y_xxxxxz[k] * ab_x + g_z_z_y_xxxxxxz[k];

                g_z_z_xy_xxxxyy[k] = -g_z_z_y_xxxxyy[k] * ab_x + g_z_z_y_xxxxxyy[k];

                g_z_z_xy_xxxxyz[k] = -g_z_z_y_xxxxyz[k] * ab_x + g_z_z_y_xxxxxyz[k];

                g_z_z_xy_xxxxzz[k] = -g_z_z_y_xxxxzz[k] * ab_x + g_z_z_y_xxxxxzz[k];

                g_z_z_xy_xxxyyy[k] = -g_z_z_y_xxxyyy[k] * ab_x + g_z_z_y_xxxxyyy[k];

                g_z_z_xy_xxxyyz[k] = -g_z_z_y_xxxyyz[k] * ab_x + g_z_z_y_xxxxyyz[k];

                g_z_z_xy_xxxyzz[k] = -g_z_z_y_xxxyzz[k] * ab_x + g_z_z_y_xxxxyzz[k];

                g_z_z_xy_xxxzzz[k] = -g_z_z_y_xxxzzz[k] * ab_x + g_z_z_y_xxxxzzz[k];

                g_z_z_xy_xxyyyy[k] = -g_z_z_y_xxyyyy[k] * ab_x + g_z_z_y_xxxyyyy[k];

                g_z_z_xy_xxyyyz[k] = -g_z_z_y_xxyyyz[k] * ab_x + g_z_z_y_xxxyyyz[k];

                g_z_z_xy_xxyyzz[k] = -g_z_z_y_xxyyzz[k] * ab_x + g_z_z_y_xxxyyzz[k];

                g_z_z_xy_xxyzzz[k] = -g_z_z_y_xxyzzz[k] * ab_x + g_z_z_y_xxxyzzz[k];

                g_z_z_xy_xxzzzz[k] = -g_z_z_y_xxzzzz[k] * ab_x + g_z_z_y_xxxzzzz[k];

                g_z_z_xy_xyyyyy[k] = -g_z_z_y_xyyyyy[k] * ab_x + g_z_z_y_xxyyyyy[k];

                g_z_z_xy_xyyyyz[k] = -g_z_z_y_xyyyyz[k] * ab_x + g_z_z_y_xxyyyyz[k];

                g_z_z_xy_xyyyzz[k] = -g_z_z_y_xyyyzz[k] * ab_x + g_z_z_y_xxyyyzz[k];

                g_z_z_xy_xyyzzz[k] = -g_z_z_y_xyyzzz[k] * ab_x + g_z_z_y_xxyyzzz[k];

                g_z_z_xy_xyzzzz[k] = -g_z_z_y_xyzzzz[k] * ab_x + g_z_z_y_xxyzzzz[k];

                g_z_z_xy_xzzzzz[k] = -g_z_z_y_xzzzzz[k] * ab_x + g_z_z_y_xxzzzzz[k];

                g_z_z_xy_yyyyyy[k] = -g_z_z_y_yyyyyy[k] * ab_x + g_z_z_y_xyyyyyy[k];

                g_z_z_xy_yyyyyz[k] = -g_z_z_y_yyyyyz[k] * ab_x + g_z_z_y_xyyyyyz[k];

                g_z_z_xy_yyyyzz[k] = -g_z_z_y_yyyyzz[k] * ab_x + g_z_z_y_xyyyyzz[k];

                g_z_z_xy_yyyzzz[k] = -g_z_z_y_yyyzzz[k] * ab_x + g_z_z_y_xyyyzzz[k];

                g_z_z_xy_yyzzzz[k] = -g_z_z_y_yyzzzz[k] * ab_x + g_z_z_y_xyyzzzz[k];

                g_z_z_xy_yzzzzz[k] = -g_z_z_y_yzzzzz[k] * ab_x + g_z_z_y_xyzzzzz[k];

                g_z_z_xy_zzzzzz[k] = -g_z_z_y_zzzzzz[k] * ab_x + g_z_z_y_xzzzzzz[k];
            }

            /// Set up 1400-1428 components of targeted buffer : cbuffer.data(

            auto g_z_z_xz_xxxxxx = cbuffer.data(di_geom_11_off + 1400 * ccomps * dcomps);

            auto g_z_z_xz_xxxxxy = cbuffer.data(di_geom_11_off + 1401 * ccomps * dcomps);

            auto g_z_z_xz_xxxxxz = cbuffer.data(di_geom_11_off + 1402 * ccomps * dcomps);

            auto g_z_z_xz_xxxxyy = cbuffer.data(di_geom_11_off + 1403 * ccomps * dcomps);

            auto g_z_z_xz_xxxxyz = cbuffer.data(di_geom_11_off + 1404 * ccomps * dcomps);

            auto g_z_z_xz_xxxxzz = cbuffer.data(di_geom_11_off + 1405 * ccomps * dcomps);

            auto g_z_z_xz_xxxyyy = cbuffer.data(di_geom_11_off + 1406 * ccomps * dcomps);

            auto g_z_z_xz_xxxyyz = cbuffer.data(di_geom_11_off + 1407 * ccomps * dcomps);

            auto g_z_z_xz_xxxyzz = cbuffer.data(di_geom_11_off + 1408 * ccomps * dcomps);

            auto g_z_z_xz_xxxzzz = cbuffer.data(di_geom_11_off + 1409 * ccomps * dcomps);

            auto g_z_z_xz_xxyyyy = cbuffer.data(di_geom_11_off + 1410 * ccomps * dcomps);

            auto g_z_z_xz_xxyyyz = cbuffer.data(di_geom_11_off + 1411 * ccomps * dcomps);

            auto g_z_z_xz_xxyyzz = cbuffer.data(di_geom_11_off + 1412 * ccomps * dcomps);

            auto g_z_z_xz_xxyzzz = cbuffer.data(di_geom_11_off + 1413 * ccomps * dcomps);

            auto g_z_z_xz_xxzzzz = cbuffer.data(di_geom_11_off + 1414 * ccomps * dcomps);

            auto g_z_z_xz_xyyyyy = cbuffer.data(di_geom_11_off + 1415 * ccomps * dcomps);

            auto g_z_z_xz_xyyyyz = cbuffer.data(di_geom_11_off + 1416 * ccomps * dcomps);

            auto g_z_z_xz_xyyyzz = cbuffer.data(di_geom_11_off + 1417 * ccomps * dcomps);

            auto g_z_z_xz_xyyzzz = cbuffer.data(di_geom_11_off + 1418 * ccomps * dcomps);

            auto g_z_z_xz_xyzzzz = cbuffer.data(di_geom_11_off + 1419 * ccomps * dcomps);

            auto g_z_z_xz_xzzzzz = cbuffer.data(di_geom_11_off + 1420 * ccomps * dcomps);

            auto g_z_z_xz_yyyyyy = cbuffer.data(di_geom_11_off + 1421 * ccomps * dcomps);

            auto g_z_z_xz_yyyyyz = cbuffer.data(di_geom_11_off + 1422 * ccomps * dcomps);

            auto g_z_z_xz_yyyyzz = cbuffer.data(di_geom_11_off + 1423 * ccomps * dcomps);

            auto g_z_z_xz_yyyzzz = cbuffer.data(di_geom_11_off + 1424 * ccomps * dcomps);

            auto g_z_z_xz_yyzzzz = cbuffer.data(di_geom_11_off + 1425 * ccomps * dcomps);

            auto g_z_z_xz_yzzzzz = cbuffer.data(di_geom_11_off + 1426 * ccomps * dcomps);

            auto g_z_z_xz_zzzzzz = cbuffer.data(di_geom_11_off + 1427 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xz_xxxxxx, g_z_z_xz_xxxxxy, g_z_z_xz_xxxxxz, g_z_z_xz_xxxxyy, g_z_z_xz_xxxxyz, g_z_z_xz_xxxxzz, g_z_z_xz_xxxyyy, g_z_z_xz_xxxyyz, g_z_z_xz_xxxyzz, g_z_z_xz_xxxzzz, g_z_z_xz_xxyyyy, g_z_z_xz_xxyyyz, g_z_z_xz_xxyyzz, g_z_z_xz_xxyzzz, g_z_z_xz_xxzzzz, g_z_z_xz_xyyyyy, g_z_z_xz_xyyyyz, g_z_z_xz_xyyyzz, g_z_z_xz_xyyzzz, g_z_z_xz_xyzzzz, g_z_z_xz_xzzzzz, g_z_z_xz_yyyyyy, g_z_z_xz_yyyyyz, g_z_z_xz_yyyyzz, g_z_z_xz_yyyzzz, g_z_z_xz_yyzzzz, g_z_z_xz_yzzzzz, g_z_z_xz_zzzzzz, g_z_z_z_xxxxxx, g_z_z_z_xxxxxxx, g_z_z_z_xxxxxxy, g_z_z_z_xxxxxxz, g_z_z_z_xxxxxy, g_z_z_z_xxxxxyy, g_z_z_z_xxxxxyz, g_z_z_z_xxxxxz, g_z_z_z_xxxxxzz, g_z_z_z_xxxxyy, g_z_z_z_xxxxyyy, g_z_z_z_xxxxyyz, g_z_z_z_xxxxyz, g_z_z_z_xxxxyzz, g_z_z_z_xxxxzz, g_z_z_z_xxxxzzz, g_z_z_z_xxxyyy, g_z_z_z_xxxyyyy, g_z_z_z_xxxyyyz, g_z_z_z_xxxyyz, g_z_z_z_xxxyyzz, g_z_z_z_xxxyzz, g_z_z_z_xxxyzzz, g_z_z_z_xxxzzz, g_z_z_z_xxxzzzz, g_z_z_z_xxyyyy, g_z_z_z_xxyyyyy, g_z_z_z_xxyyyyz, g_z_z_z_xxyyyz, g_z_z_z_xxyyyzz, g_z_z_z_xxyyzz, g_z_z_z_xxyyzzz, g_z_z_z_xxyzzz, g_z_z_z_xxyzzzz, g_z_z_z_xxzzzz, g_z_z_z_xxzzzzz, g_z_z_z_xyyyyy, g_z_z_z_xyyyyyy, g_z_z_z_xyyyyyz, g_z_z_z_xyyyyz, g_z_z_z_xyyyyzz, g_z_z_z_xyyyzz, g_z_z_z_xyyyzzz, g_z_z_z_xyyzzz, g_z_z_z_xyyzzzz, g_z_z_z_xyzzzz, g_z_z_z_xyzzzzz, g_z_z_z_xzzzzz, g_z_z_z_xzzzzzz, g_z_z_z_yyyyyy, g_z_z_z_yyyyyz, g_z_z_z_yyyyzz, g_z_z_z_yyyzzz, g_z_z_z_yyzzzz, g_z_z_z_yzzzzz, g_z_z_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xz_xxxxxx[k] = -g_z_z_z_xxxxxx[k] * ab_x + g_z_z_z_xxxxxxx[k];

                g_z_z_xz_xxxxxy[k] = -g_z_z_z_xxxxxy[k] * ab_x + g_z_z_z_xxxxxxy[k];

                g_z_z_xz_xxxxxz[k] = -g_z_z_z_xxxxxz[k] * ab_x + g_z_z_z_xxxxxxz[k];

                g_z_z_xz_xxxxyy[k] = -g_z_z_z_xxxxyy[k] * ab_x + g_z_z_z_xxxxxyy[k];

                g_z_z_xz_xxxxyz[k] = -g_z_z_z_xxxxyz[k] * ab_x + g_z_z_z_xxxxxyz[k];

                g_z_z_xz_xxxxzz[k] = -g_z_z_z_xxxxzz[k] * ab_x + g_z_z_z_xxxxxzz[k];

                g_z_z_xz_xxxyyy[k] = -g_z_z_z_xxxyyy[k] * ab_x + g_z_z_z_xxxxyyy[k];

                g_z_z_xz_xxxyyz[k] = -g_z_z_z_xxxyyz[k] * ab_x + g_z_z_z_xxxxyyz[k];

                g_z_z_xz_xxxyzz[k] = -g_z_z_z_xxxyzz[k] * ab_x + g_z_z_z_xxxxyzz[k];

                g_z_z_xz_xxxzzz[k] = -g_z_z_z_xxxzzz[k] * ab_x + g_z_z_z_xxxxzzz[k];

                g_z_z_xz_xxyyyy[k] = -g_z_z_z_xxyyyy[k] * ab_x + g_z_z_z_xxxyyyy[k];

                g_z_z_xz_xxyyyz[k] = -g_z_z_z_xxyyyz[k] * ab_x + g_z_z_z_xxxyyyz[k];

                g_z_z_xz_xxyyzz[k] = -g_z_z_z_xxyyzz[k] * ab_x + g_z_z_z_xxxyyzz[k];

                g_z_z_xz_xxyzzz[k] = -g_z_z_z_xxyzzz[k] * ab_x + g_z_z_z_xxxyzzz[k];

                g_z_z_xz_xxzzzz[k] = -g_z_z_z_xxzzzz[k] * ab_x + g_z_z_z_xxxzzzz[k];

                g_z_z_xz_xyyyyy[k] = -g_z_z_z_xyyyyy[k] * ab_x + g_z_z_z_xxyyyyy[k];

                g_z_z_xz_xyyyyz[k] = -g_z_z_z_xyyyyz[k] * ab_x + g_z_z_z_xxyyyyz[k];

                g_z_z_xz_xyyyzz[k] = -g_z_z_z_xyyyzz[k] * ab_x + g_z_z_z_xxyyyzz[k];

                g_z_z_xz_xyyzzz[k] = -g_z_z_z_xyyzzz[k] * ab_x + g_z_z_z_xxyyzzz[k];

                g_z_z_xz_xyzzzz[k] = -g_z_z_z_xyzzzz[k] * ab_x + g_z_z_z_xxyzzzz[k];

                g_z_z_xz_xzzzzz[k] = -g_z_z_z_xzzzzz[k] * ab_x + g_z_z_z_xxzzzzz[k];

                g_z_z_xz_yyyyyy[k] = -g_z_z_z_yyyyyy[k] * ab_x + g_z_z_z_xyyyyyy[k];

                g_z_z_xz_yyyyyz[k] = -g_z_z_z_yyyyyz[k] * ab_x + g_z_z_z_xyyyyyz[k];

                g_z_z_xz_yyyyzz[k] = -g_z_z_z_yyyyzz[k] * ab_x + g_z_z_z_xyyyyzz[k];

                g_z_z_xz_yyyzzz[k] = -g_z_z_z_yyyzzz[k] * ab_x + g_z_z_z_xyyyzzz[k];

                g_z_z_xz_yyzzzz[k] = -g_z_z_z_yyzzzz[k] * ab_x + g_z_z_z_xyyzzzz[k];

                g_z_z_xz_yzzzzz[k] = -g_z_z_z_yzzzzz[k] * ab_x + g_z_z_z_xyzzzzz[k];

                g_z_z_xz_zzzzzz[k] = -g_z_z_z_zzzzzz[k] * ab_x + g_z_z_z_xzzzzzz[k];
            }

            /// Set up 1428-1456 components of targeted buffer : cbuffer.data(

            auto g_z_z_yy_xxxxxx = cbuffer.data(di_geom_11_off + 1428 * ccomps * dcomps);

            auto g_z_z_yy_xxxxxy = cbuffer.data(di_geom_11_off + 1429 * ccomps * dcomps);

            auto g_z_z_yy_xxxxxz = cbuffer.data(di_geom_11_off + 1430 * ccomps * dcomps);

            auto g_z_z_yy_xxxxyy = cbuffer.data(di_geom_11_off + 1431 * ccomps * dcomps);

            auto g_z_z_yy_xxxxyz = cbuffer.data(di_geom_11_off + 1432 * ccomps * dcomps);

            auto g_z_z_yy_xxxxzz = cbuffer.data(di_geom_11_off + 1433 * ccomps * dcomps);

            auto g_z_z_yy_xxxyyy = cbuffer.data(di_geom_11_off + 1434 * ccomps * dcomps);

            auto g_z_z_yy_xxxyyz = cbuffer.data(di_geom_11_off + 1435 * ccomps * dcomps);

            auto g_z_z_yy_xxxyzz = cbuffer.data(di_geom_11_off + 1436 * ccomps * dcomps);

            auto g_z_z_yy_xxxzzz = cbuffer.data(di_geom_11_off + 1437 * ccomps * dcomps);

            auto g_z_z_yy_xxyyyy = cbuffer.data(di_geom_11_off + 1438 * ccomps * dcomps);

            auto g_z_z_yy_xxyyyz = cbuffer.data(di_geom_11_off + 1439 * ccomps * dcomps);

            auto g_z_z_yy_xxyyzz = cbuffer.data(di_geom_11_off + 1440 * ccomps * dcomps);

            auto g_z_z_yy_xxyzzz = cbuffer.data(di_geom_11_off + 1441 * ccomps * dcomps);

            auto g_z_z_yy_xxzzzz = cbuffer.data(di_geom_11_off + 1442 * ccomps * dcomps);

            auto g_z_z_yy_xyyyyy = cbuffer.data(di_geom_11_off + 1443 * ccomps * dcomps);

            auto g_z_z_yy_xyyyyz = cbuffer.data(di_geom_11_off + 1444 * ccomps * dcomps);

            auto g_z_z_yy_xyyyzz = cbuffer.data(di_geom_11_off + 1445 * ccomps * dcomps);

            auto g_z_z_yy_xyyzzz = cbuffer.data(di_geom_11_off + 1446 * ccomps * dcomps);

            auto g_z_z_yy_xyzzzz = cbuffer.data(di_geom_11_off + 1447 * ccomps * dcomps);

            auto g_z_z_yy_xzzzzz = cbuffer.data(di_geom_11_off + 1448 * ccomps * dcomps);

            auto g_z_z_yy_yyyyyy = cbuffer.data(di_geom_11_off + 1449 * ccomps * dcomps);

            auto g_z_z_yy_yyyyyz = cbuffer.data(di_geom_11_off + 1450 * ccomps * dcomps);

            auto g_z_z_yy_yyyyzz = cbuffer.data(di_geom_11_off + 1451 * ccomps * dcomps);

            auto g_z_z_yy_yyyzzz = cbuffer.data(di_geom_11_off + 1452 * ccomps * dcomps);

            auto g_z_z_yy_yyzzzz = cbuffer.data(di_geom_11_off + 1453 * ccomps * dcomps);

            auto g_z_z_yy_yzzzzz = cbuffer.data(di_geom_11_off + 1454 * ccomps * dcomps);

            auto g_z_z_yy_zzzzzz = cbuffer.data(di_geom_11_off + 1455 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_y_xxxxxx, g_z_z_y_xxxxxxy, g_z_z_y_xxxxxy, g_z_z_y_xxxxxyy, g_z_z_y_xxxxxyz, g_z_z_y_xxxxxz, g_z_z_y_xxxxyy, g_z_z_y_xxxxyyy, g_z_z_y_xxxxyyz, g_z_z_y_xxxxyz, g_z_z_y_xxxxyzz, g_z_z_y_xxxxzz, g_z_z_y_xxxyyy, g_z_z_y_xxxyyyy, g_z_z_y_xxxyyyz, g_z_z_y_xxxyyz, g_z_z_y_xxxyyzz, g_z_z_y_xxxyzz, g_z_z_y_xxxyzzz, g_z_z_y_xxxzzz, g_z_z_y_xxyyyy, g_z_z_y_xxyyyyy, g_z_z_y_xxyyyyz, g_z_z_y_xxyyyz, g_z_z_y_xxyyyzz, g_z_z_y_xxyyzz, g_z_z_y_xxyyzzz, g_z_z_y_xxyzzz, g_z_z_y_xxyzzzz, g_z_z_y_xxzzzz, g_z_z_y_xyyyyy, g_z_z_y_xyyyyyy, g_z_z_y_xyyyyyz, g_z_z_y_xyyyyz, g_z_z_y_xyyyyzz, g_z_z_y_xyyyzz, g_z_z_y_xyyyzzz, g_z_z_y_xyyzzz, g_z_z_y_xyyzzzz, g_z_z_y_xyzzzz, g_z_z_y_xyzzzzz, g_z_z_y_xzzzzz, g_z_z_y_yyyyyy, g_z_z_y_yyyyyyy, g_z_z_y_yyyyyyz, g_z_z_y_yyyyyz, g_z_z_y_yyyyyzz, g_z_z_y_yyyyzz, g_z_z_y_yyyyzzz, g_z_z_y_yyyzzz, g_z_z_y_yyyzzzz, g_z_z_y_yyzzzz, g_z_z_y_yyzzzzz, g_z_z_y_yzzzzz, g_z_z_y_yzzzzzz, g_z_z_y_zzzzzz, g_z_z_yy_xxxxxx, g_z_z_yy_xxxxxy, g_z_z_yy_xxxxxz, g_z_z_yy_xxxxyy, g_z_z_yy_xxxxyz, g_z_z_yy_xxxxzz, g_z_z_yy_xxxyyy, g_z_z_yy_xxxyyz, g_z_z_yy_xxxyzz, g_z_z_yy_xxxzzz, g_z_z_yy_xxyyyy, g_z_z_yy_xxyyyz, g_z_z_yy_xxyyzz, g_z_z_yy_xxyzzz, g_z_z_yy_xxzzzz, g_z_z_yy_xyyyyy, g_z_z_yy_xyyyyz, g_z_z_yy_xyyyzz, g_z_z_yy_xyyzzz, g_z_z_yy_xyzzzz, g_z_z_yy_xzzzzz, g_z_z_yy_yyyyyy, g_z_z_yy_yyyyyz, g_z_z_yy_yyyyzz, g_z_z_yy_yyyzzz, g_z_z_yy_yyzzzz, g_z_z_yy_yzzzzz, g_z_z_yy_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yy_xxxxxx[k] = -g_z_z_y_xxxxxx[k] * ab_y + g_z_z_y_xxxxxxy[k];

                g_z_z_yy_xxxxxy[k] = -g_z_z_y_xxxxxy[k] * ab_y + g_z_z_y_xxxxxyy[k];

                g_z_z_yy_xxxxxz[k] = -g_z_z_y_xxxxxz[k] * ab_y + g_z_z_y_xxxxxyz[k];

                g_z_z_yy_xxxxyy[k] = -g_z_z_y_xxxxyy[k] * ab_y + g_z_z_y_xxxxyyy[k];

                g_z_z_yy_xxxxyz[k] = -g_z_z_y_xxxxyz[k] * ab_y + g_z_z_y_xxxxyyz[k];

                g_z_z_yy_xxxxzz[k] = -g_z_z_y_xxxxzz[k] * ab_y + g_z_z_y_xxxxyzz[k];

                g_z_z_yy_xxxyyy[k] = -g_z_z_y_xxxyyy[k] * ab_y + g_z_z_y_xxxyyyy[k];

                g_z_z_yy_xxxyyz[k] = -g_z_z_y_xxxyyz[k] * ab_y + g_z_z_y_xxxyyyz[k];

                g_z_z_yy_xxxyzz[k] = -g_z_z_y_xxxyzz[k] * ab_y + g_z_z_y_xxxyyzz[k];

                g_z_z_yy_xxxzzz[k] = -g_z_z_y_xxxzzz[k] * ab_y + g_z_z_y_xxxyzzz[k];

                g_z_z_yy_xxyyyy[k] = -g_z_z_y_xxyyyy[k] * ab_y + g_z_z_y_xxyyyyy[k];

                g_z_z_yy_xxyyyz[k] = -g_z_z_y_xxyyyz[k] * ab_y + g_z_z_y_xxyyyyz[k];

                g_z_z_yy_xxyyzz[k] = -g_z_z_y_xxyyzz[k] * ab_y + g_z_z_y_xxyyyzz[k];

                g_z_z_yy_xxyzzz[k] = -g_z_z_y_xxyzzz[k] * ab_y + g_z_z_y_xxyyzzz[k];

                g_z_z_yy_xxzzzz[k] = -g_z_z_y_xxzzzz[k] * ab_y + g_z_z_y_xxyzzzz[k];

                g_z_z_yy_xyyyyy[k] = -g_z_z_y_xyyyyy[k] * ab_y + g_z_z_y_xyyyyyy[k];

                g_z_z_yy_xyyyyz[k] = -g_z_z_y_xyyyyz[k] * ab_y + g_z_z_y_xyyyyyz[k];

                g_z_z_yy_xyyyzz[k] = -g_z_z_y_xyyyzz[k] * ab_y + g_z_z_y_xyyyyzz[k];

                g_z_z_yy_xyyzzz[k] = -g_z_z_y_xyyzzz[k] * ab_y + g_z_z_y_xyyyzzz[k];

                g_z_z_yy_xyzzzz[k] = -g_z_z_y_xyzzzz[k] * ab_y + g_z_z_y_xyyzzzz[k];

                g_z_z_yy_xzzzzz[k] = -g_z_z_y_xzzzzz[k] * ab_y + g_z_z_y_xyzzzzz[k];

                g_z_z_yy_yyyyyy[k] = -g_z_z_y_yyyyyy[k] * ab_y + g_z_z_y_yyyyyyy[k];

                g_z_z_yy_yyyyyz[k] = -g_z_z_y_yyyyyz[k] * ab_y + g_z_z_y_yyyyyyz[k];

                g_z_z_yy_yyyyzz[k] = -g_z_z_y_yyyyzz[k] * ab_y + g_z_z_y_yyyyyzz[k];

                g_z_z_yy_yyyzzz[k] = -g_z_z_y_yyyzzz[k] * ab_y + g_z_z_y_yyyyzzz[k];

                g_z_z_yy_yyzzzz[k] = -g_z_z_y_yyzzzz[k] * ab_y + g_z_z_y_yyyzzzz[k];

                g_z_z_yy_yzzzzz[k] = -g_z_z_y_yzzzzz[k] * ab_y + g_z_z_y_yyzzzzz[k];

                g_z_z_yy_zzzzzz[k] = -g_z_z_y_zzzzzz[k] * ab_y + g_z_z_y_yzzzzzz[k];
            }

            /// Set up 1456-1484 components of targeted buffer : cbuffer.data(

            auto g_z_z_yz_xxxxxx = cbuffer.data(di_geom_11_off + 1456 * ccomps * dcomps);

            auto g_z_z_yz_xxxxxy = cbuffer.data(di_geom_11_off + 1457 * ccomps * dcomps);

            auto g_z_z_yz_xxxxxz = cbuffer.data(di_geom_11_off + 1458 * ccomps * dcomps);

            auto g_z_z_yz_xxxxyy = cbuffer.data(di_geom_11_off + 1459 * ccomps * dcomps);

            auto g_z_z_yz_xxxxyz = cbuffer.data(di_geom_11_off + 1460 * ccomps * dcomps);

            auto g_z_z_yz_xxxxzz = cbuffer.data(di_geom_11_off + 1461 * ccomps * dcomps);

            auto g_z_z_yz_xxxyyy = cbuffer.data(di_geom_11_off + 1462 * ccomps * dcomps);

            auto g_z_z_yz_xxxyyz = cbuffer.data(di_geom_11_off + 1463 * ccomps * dcomps);

            auto g_z_z_yz_xxxyzz = cbuffer.data(di_geom_11_off + 1464 * ccomps * dcomps);

            auto g_z_z_yz_xxxzzz = cbuffer.data(di_geom_11_off + 1465 * ccomps * dcomps);

            auto g_z_z_yz_xxyyyy = cbuffer.data(di_geom_11_off + 1466 * ccomps * dcomps);

            auto g_z_z_yz_xxyyyz = cbuffer.data(di_geom_11_off + 1467 * ccomps * dcomps);

            auto g_z_z_yz_xxyyzz = cbuffer.data(di_geom_11_off + 1468 * ccomps * dcomps);

            auto g_z_z_yz_xxyzzz = cbuffer.data(di_geom_11_off + 1469 * ccomps * dcomps);

            auto g_z_z_yz_xxzzzz = cbuffer.data(di_geom_11_off + 1470 * ccomps * dcomps);

            auto g_z_z_yz_xyyyyy = cbuffer.data(di_geom_11_off + 1471 * ccomps * dcomps);

            auto g_z_z_yz_xyyyyz = cbuffer.data(di_geom_11_off + 1472 * ccomps * dcomps);

            auto g_z_z_yz_xyyyzz = cbuffer.data(di_geom_11_off + 1473 * ccomps * dcomps);

            auto g_z_z_yz_xyyzzz = cbuffer.data(di_geom_11_off + 1474 * ccomps * dcomps);

            auto g_z_z_yz_xyzzzz = cbuffer.data(di_geom_11_off + 1475 * ccomps * dcomps);

            auto g_z_z_yz_xzzzzz = cbuffer.data(di_geom_11_off + 1476 * ccomps * dcomps);

            auto g_z_z_yz_yyyyyy = cbuffer.data(di_geom_11_off + 1477 * ccomps * dcomps);

            auto g_z_z_yz_yyyyyz = cbuffer.data(di_geom_11_off + 1478 * ccomps * dcomps);

            auto g_z_z_yz_yyyyzz = cbuffer.data(di_geom_11_off + 1479 * ccomps * dcomps);

            auto g_z_z_yz_yyyzzz = cbuffer.data(di_geom_11_off + 1480 * ccomps * dcomps);

            auto g_z_z_yz_yyzzzz = cbuffer.data(di_geom_11_off + 1481 * ccomps * dcomps);

            auto g_z_z_yz_yzzzzz = cbuffer.data(di_geom_11_off + 1482 * ccomps * dcomps);

            auto g_z_z_yz_zzzzzz = cbuffer.data(di_geom_11_off + 1483 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yz_xxxxxx, g_z_z_yz_xxxxxy, g_z_z_yz_xxxxxz, g_z_z_yz_xxxxyy, g_z_z_yz_xxxxyz, g_z_z_yz_xxxxzz, g_z_z_yz_xxxyyy, g_z_z_yz_xxxyyz, g_z_z_yz_xxxyzz, g_z_z_yz_xxxzzz, g_z_z_yz_xxyyyy, g_z_z_yz_xxyyyz, g_z_z_yz_xxyyzz, g_z_z_yz_xxyzzz, g_z_z_yz_xxzzzz, g_z_z_yz_xyyyyy, g_z_z_yz_xyyyyz, g_z_z_yz_xyyyzz, g_z_z_yz_xyyzzz, g_z_z_yz_xyzzzz, g_z_z_yz_xzzzzz, g_z_z_yz_yyyyyy, g_z_z_yz_yyyyyz, g_z_z_yz_yyyyzz, g_z_z_yz_yyyzzz, g_z_z_yz_yyzzzz, g_z_z_yz_yzzzzz, g_z_z_yz_zzzzzz, g_z_z_z_xxxxxx, g_z_z_z_xxxxxxy, g_z_z_z_xxxxxy, g_z_z_z_xxxxxyy, g_z_z_z_xxxxxyz, g_z_z_z_xxxxxz, g_z_z_z_xxxxyy, g_z_z_z_xxxxyyy, g_z_z_z_xxxxyyz, g_z_z_z_xxxxyz, g_z_z_z_xxxxyzz, g_z_z_z_xxxxzz, g_z_z_z_xxxyyy, g_z_z_z_xxxyyyy, g_z_z_z_xxxyyyz, g_z_z_z_xxxyyz, g_z_z_z_xxxyyzz, g_z_z_z_xxxyzz, g_z_z_z_xxxyzzz, g_z_z_z_xxxzzz, g_z_z_z_xxyyyy, g_z_z_z_xxyyyyy, g_z_z_z_xxyyyyz, g_z_z_z_xxyyyz, g_z_z_z_xxyyyzz, g_z_z_z_xxyyzz, g_z_z_z_xxyyzzz, g_z_z_z_xxyzzz, g_z_z_z_xxyzzzz, g_z_z_z_xxzzzz, g_z_z_z_xyyyyy, g_z_z_z_xyyyyyy, g_z_z_z_xyyyyyz, g_z_z_z_xyyyyz, g_z_z_z_xyyyyzz, g_z_z_z_xyyyzz, g_z_z_z_xyyyzzz, g_z_z_z_xyyzzz, g_z_z_z_xyyzzzz, g_z_z_z_xyzzzz, g_z_z_z_xyzzzzz, g_z_z_z_xzzzzz, g_z_z_z_yyyyyy, g_z_z_z_yyyyyyy, g_z_z_z_yyyyyyz, g_z_z_z_yyyyyz, g_z_z_z_yyyyyzz, g_z_z_z_yyyyzz, g_z_z_z_yyyyzzz, g_z_z_z_yyyzzz, g_z_z_z_yyyzzzz, g_z_z_z_yyzzzz, g_z_z_z_yyzzzzz, g_z_z_z_yzzzzz, g_z_z_z_yzzzzzz, g_z_z_z_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yz_xxxxxx[k] = -g_z_z_z_xxxxxx[k] * ab_y + g_z_z_z_xxxxxxy[k];

                g_z_z_yz_xxxxxy[k] = -g_z_z_z_xxxxxy[k] * ab_y + g_z_z_z_xxxxxyy[k];

                g_z_z_yz_xxxxxz[k] = -g_z_z_z_xxxxxz[k] * ab_y + g_z_z_z_xxxxxyz[k];

                g_z_z_yz_xxxxyy[k] = -g_z_z_z_xxxxyy[k] * ab_y + g_z_z_z_xxxxyyy[k];

                g_z_z_yz_xxxxyz[k] = -g_z_z_z_xxxxyz[k] * ab_y + g_z_z_z_xxxxyyz[k];

                g_z_z_yz_xxxxzz[k] = -g_z_z_z_xxxxzz[k] * ab_y + g_z_z_z_xxxxyzz[k];

                g_z_z_yz_xxxyyy[k] = -g_z_z_z_xxxyyy[k] * ab_y + g_z_z_z_xxxyyyy[k];

                g_z_z_yz_xxxyyz[k] = -g_z_z_z_xxxyyz[k] * ab_y + g_z_z_z_xxxyyyz[k];

                g_z_z_yz_xxxyzz[k] = -g_z_z_z_xxxyzz[k] * ab_y + g_z_z_z_xxxyyzz[k];

                g_z_z_yz_xxxzzz[k] = -g_z_z_z_xxxzzz[k] * ab_y + g_z_z_z_xxxyzzz[k];

                g_z_z_yz_xxyyyy[k] = -g_z_z_z_xxyyyy[k] * ab_y + g_z_z_z_xxyyyyy[k];

                g_z_z_yz_xxyyyz[k] = -g_z_z_z_xxyyyz[k] * ab_y + g_z_z_z_xxyyyyz[k];

                g_z_z_yz_xxyyzz[k] = -g_z_z_z_xxyyzz[k] * ab_y + g_z_z_z_xxyyyzz[k];

                g_z_z_yz_xxyzzz[k] = -g_z_z_z_xxyzzz[k] * ab_y + g_z_z_z_xxyyzzz[k];

                g_z_z_yz_xxzzzz[k] = -g_z_z_z_xxzzzz[k] * ab_y + g_z_z_z_xxyzzzz[k];

                g_z_z_yz_xyyyyy[k] = -g_z_z_z_xyyyyy[k] * ab_y + g_z_z_z_xyyyyyy[k];

                g_z_z_yz_xyyyyz[k] = -g_z_z_z_xyyyyz[k] * ab_y + g_z_z_z_xyyyyyz[k];

                g_z_z_yz_xyyyzz[k] = -g_z_z_z_xyyyzz[k] * ab_y + g_z_z_z_xyyyyzz[k];

                g_z_z_yz_xyyzzz[k] = -g_z_z_z_xyyzzz[k] * ab_y + g_z_z_z_xyyyzzz[k];

                g_z_z_yz_xyzzzz[k] = -g_z_z_z_xyzzzz[k] * ab_y + g_z_z_z_xyyzzzz[k];

                g_z_z_yz_xzzzzz[k] = -g_z_z_z_xzzzzz[k] * ab_y + g_z_z_z_xyzzzzz[k];

                g_z_z_yz_yyyyyy[k] = -g_z_z_z_yyyyyy[k] * ab_y + g_z_z_z_yyyyyyy[k];

                g_z_z_yz_yyyyyz[k] = -g_z_z_z_yyyyyz[k] * ab_y + g_z_z_z_yyyyyyz[k];

                g_z_z_yz_yyyyzz[k] = -g_z_z_z_yyyyzz[k] * ab_y + g_z_z_z_yyyyyzz[k];

                g_z_z_yz_yyyzzz[k] = -g_z_z_z_yyyzzz[k] * ab_y + g_z_z_z_yyyyzzz[k];

                g_z_z_yz_yyzzzz[k] = -g_z_z_z_yyzzzz[k] * ab_y + g_z_z_z_yyyzzzz[k];

                g_z_z_yz_yzzzzz[k] = -g_z_z_z_yzzzzz[k] * ab_y + g_z_z_z_yyzzzzz[k];

                g_z_z_yz_zzzzzz[k] = -g_z_z_z_zzzzzz[k] * ab_y + g_z_z_z_yzzzzzz[k];
            }

            /// Set up 1484-1512 components of targeted buffer : cbuffer.data(

            auto g_z_z_zz_xxxxxx = cbuffer.data(di_geom_11_off + 1484 * ccomps * dcomps);

            auto g_z_z_zz_xxxxxy = cbuffer.data(di_geom_11_off + 1485 * ccomps * dcomps);

            auto g_z_z_zz_xxxxxz = cbuffer.data(di_geom_11_off + 1486 * ccomps * dcomps);

            auto g_z_z_zz_xxxxyy = cbuffer.data(di_geom_11_off + 1487 * ccomps * dcomps);

            auto g_z_z_zz_xxxxyz = cbuffer.data(di_geom_11_off + 1488 * ccomps * dcomps);

            auto g_z_z_zz_xxxxzz = cbuffer.data(di_geom_11_off + 1489 * ccomps * dcomps);

            auto g_z_z_zz_xxxyyy = cbuffer.data(di_geom_11_off + 1490 * ccomps * dcomps);

            auto g_z_z_zz_xxxyyz = cbuffer.data(di_geom_11_off + 1491 * ccomps * dcomps);

            auto g_z_z_zz_xxxyzz = cbuffer.data(di_geom_11_off + 1492 * ccomps * dcomps);

            auto g_z_z_zz_xxxzzz = cbuffer.data(di_geom_11_off + 1493 * ccomps * dcomps);

            auto g_z_z_zz_xxyyyy = cbuffer.data(di_geom_11_off + 1494 * ccomps * dcomps);

            auto g_z_z_zz_xxyyyz = cbuffer.data(di_geom_11_off + 1495 * ccomps * dcomps);

            auto g_z_z_zz_xxyyzz = cbuffer.data(di_geom_11_off + 1496 * ccomps * dcomps);

            auto g_z_z_zz_xxyzzz = cbuffer.data(di_geom_11_off + 1497 * ccomps * dcomps);

            auto g_z_z_zz_xxzzzz = cbuffer.data(di_geom_11_off + 1498 * ccomps * dcomps);

            auto g_z_z_zz_xyyyyy = cbuffer.data(di_geom_11_off + 1499 * ccomps * dcomps);

            auto g_z_z_zz_xyyyyz = cbuffer.data(di_geom_11_off + 1500 * ccomps * dcomps);

            auto g_z_z_zz_xyyyzz = cbuffer.data(di_geom_11_off + 1501 * ccomps * dcomps);

            auto g_z_z_zz_xyyzzz = cbuffer.data(di_geom_11_off + 1502 * ccomps * dcomps);

            auto g_z_z_zz_xyzzzz = cbuffer.data(di_geom_11_off + 1503 * ccomps * dcomps);

            auto g_z_z_zz_xzzzzz = cbuffer.data(di_geom_11_off + 1504 * ccomps * dcomps);

            auto g_z_z_zz_yyyyyy = cbuffer.data(di_geom_11_off + 1505 * ccomps * dcomps);

            auto g_z_z_zz_yyyyyz = cbuffer.data(di_geom_11_off + 1506 * ccomps * dcomps);

            auto g_z_z_zz_yyyyzz = cbuffer.data(di_geom_11_off + 1507 * ccomps * dcomps);

            auto g_z_z_zz_yyyzzz = cbuffer.data(di_geom_11_off + 1508 * ccomps * dcomps);

            auto g_z_z_zz_yyzzzz = cbuffer.data(di_geom_11_off + 1509 * ccomps * dcomps);

            auto g_z_z_zz_yzzzzz = cbuffer.data(di_geom_11_off + 1510 * ccomps * dcomps);

            auto g_z_z_zz_zzzzzz = cbuffer.data(di_geom_11_off + 1511 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxxxx, g_0_z_z_xxxxxy, g_0_z_z_xxxxxz, g_0_z_z_xxxxyy, g_0_z_z_xxxxyz, g_0_z_z_xxxxzz, g_0_z_z_xxxyyy, g_0_z_z_xxxyyz, g_0_z_z_xxxyzz, g_0_z_z_xxxzzz, g_0_z_z_xxyyyy, g_0_z_z_xxyyyz, g_0_z_z_xxyyzz, g_0_z_z_xxyzzz, g_0_z_z_xxzzzz, g_0_z_z_xyyyyy, g_0_z_z_xyyyyz, g_0_z_z_xyyyzz, g_0_z_z_xyyzzz, g_0_z_z_xyzzzz, g_0_z_z_xzzzzz, g_0_z_z_yyyyyy, g_0_z_z_yyyyyz, g_0_z_z_yyyyzz, g_0_z_z_yyyzzz, g_0_z_z_yyzzzz, g_0_z_z_yzzzzz, g_0_z_z_zzzzzz, g_z_0_z_xxxxxx, g_z_0_z_xxxxxy, g_z_0_z_xxxxxz, g_z_0_z_xxxxyy, g_z_0_z_xxxxyz, g_z_0_z_xxxxzz, g_z_0_z_xxxyyy, g_z_0_z_xxxyyz, g_z_0_z_xxxyzz, g_z_0_z_xxxzzz, g_z_0_z_xxyyyy, g_z_0_z_xxyyyz, g_z_0_z_xxyyzz, g_z_0_z_xxyzzz, g_z_0_z_xxzzzz, g_z_0_z_xyyyyy, g_z_0_z_xyyyyz, g_z_0_z_xyyyzz, g_z_0_z_xyyzzz, g_z_0_z_xyzzzz, g_z_0_z_xzzzzz, g_z_0_z_yyyyyy, g_z_0_z_yyyyyz, g_z_0_z_yyyyzz, g_z_0_z_yyyzzz, g_z_0_z_yyzzzz, g_z_0_z_yzzzzz, g_z_0_z_zzzzzz, g_z_z_z_xxxxxx, g_z_z_z_xxxxxxz, g_z_z_z_xxxxxy, g_z_z_z_xxxxxyz, g_z_z_z_xxxxxz, g_z_z_z_xxxxxzz, g_z_z_z_xxxxyy, g_z_z_z_xxxxyyz, g_z_z_z_xxxxyz, g_z_z_z_xxxxyzz, g_z_z_z_xxxxzz, g_z_z_z_xxxxzzz, g_z_z_z_xxxyyy, g_z_z_z_xxxyyyz, g_z_z_z_xxxyyz, g_z_z_z_xxxyyzz, g_z_z_z_xxxyzz, g_z_z_z_xxxyzzz, g_z_z_z_xxxzzz, g_z_z_z_xxxzzzz, g_z_z_z_xxyyyy, g_z_z_z_xxyyyyz, g_z_z_z_xxyyyz, g_z_z_z_xxyyyzz, g_z_z_z_xxyyzz, g_z_z_z_xxyyzzz, g_z_z_z_xxyzzz, g_z_z_z_xxyzzzz, g_z_z_z_xxzzzz, g_z_z_z_xxzzzzz, g_z_z_z_xyyyyy, g_z_z_z_xyyyyyz, g_z_z_z_xyyyyz, g_z_z_z_xyyyyzz, g_z_z_z_xyyyzz, g_z_z_z_xyyyzzz, g_z_z_z_xyyzzz, g_z_z_z_xyyzzzz, g_z_z_z_xyzzzz, g_z_z_z_xyzzzzz, g_z_z_z_xzzzzz, g_z_z_z_xzzzzzz, g_z_z_z_yyyyyy, g_z_z_z_yyyyyyz, g_z_z_z_yyyyyz, g_z_z_z_yyyyyzz, g_z_z_z_yyyyzz, g_z_z_z_yyyyzzz, g_z_z_z_yyyzzz, g_z_z_z_yyyzzzz, g_z_z_z_yyzzzz, g_z_z_z_yyzzzzz, g_z_z_z_yzzzzz, g_z_z_z_yzzzzzz, g_z_z_z_zzzzzz, g_z_z_z_zzzzzzz, g_z_z_zz_xxxxxx, g_z_z_zz_xxxxxy, g_z_z_zz_xxxxxz, g_z_z_zz_xxxxyy, g_z_z_zz_xxxxyz, g_z_z_zz_xxxxzz, g_z_z_zz_xxxyyy, g_z_z_zz_xxxyyz, g_z_z_zz_xxxyzz, g_z_z_zz_xxxzzz, g_z_z_zz_xxyyyy, g_z_z_zz_xxyyyz, g_z_z_zz_xxyyzz, g_z_z_zz_xxyzzz, g_z_z_zz_xxzzzz, g_z_z_zz_xyyyyy, g_z_z_zz_xyyyyz, g_z_z_zz_xyyyzz, g_z_z_zz_xyyzzz, g_z_z_zz_xyzzzz, g_z_z_zz_xzzzzz, g_z_z_zz_yyyyyy, g_z_z_zz_yyyyyz, g_z_z_zz_yyyyzz, g_z_z_zz_yyyzzz, g_z_z_zz_yyzzzz, g_z_z_zz_yzzzzz, g_z_z_zz_zzzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zz_xxxxxx[k] = -g_0_z_z_xxxxxx[k] + g_z_0_z_xxxxxx[k] - g_z_z_z_xxxxxx[k] * ab_z + g_z_z_z_xxxxxxz[k];

                g_z_z_zz_xxxxxy[k] = -g_0_z_z_xxxxxy[k] + g_z_0_z_xxxxxy[k] - g_z_z_z_xxxxxy[k] * ab_z + g_z_z_z_xxxxxyz[k];

                g_z_z_zz_xxxxxz[k] = -g_0_z_z_xxxxxz[k] + g_z_0_z_xxxxxz[k] - g_z_z_z_xxxxxz[k] * ab_z + g_z_z_z_xxxxxzz[k];

                g_z_z_zz_xxxxyy[k] = -g_0_z_z_xxxxyy[k] + g_z_0_z_xxxxyy[k] - g_z_z_z_xxxxyy[k] * ab_z + g_z_z_z_xxxxyyz[k];

                g_z_z_zz_xxxxyz[k] = -g_0_z_z_xxxxyz[k] + g_z_0_z_xxxxyz[k] - g_z_z_z_xxxxyz[k] * ab_z + g_z_z_z_xxxxyzz[k];

                g_z_z_zz_xxxxzz[k] = -g_0_z_z_xxxxzz[k] + g_z_0_z_xxxxzz[k] - g_z_z_z_xxxxzz[k] * ab_z + g_z_z_z_xxxxzzz[k];

                g_z_z_zz_xxxyyy[k] = -g_0_z_z_xxxyyy[k] + g_z_0_z_xxxyyy[k] - g_z_z_z_xxxyyy[k] * ab_z + g_z_z_z_xxxyyyz[k];

                g_z_z_zz_xxxyyz[k] = -g_0_z_z_xxxyyz[k] + g_z_0_z_xxxyyz[k] - g_z_z_z_xxxyyz[k] * ab_z + g_z_z_z_xxxyyzz[k];

                g_z_z_zz_xxxyzz[k] = -g_0_z_z_xxxyzz[k] + g_z_0_z_xxxyzz[k] - g_z_z_z_xxxyzz[k] * ab_z + g_z_z_z_xxxyzzz[k];

                g_z_z_zz_xxxzzz[k] = -g_0_z_z_xxxzzz[k] + g_z_0_z_xxxzzz[k] - g_z_z_z_xxxzzz[k] * ab_z + g_z_z_z_xxxzzzz[k];

                g_z_z_zz_xxyyyy[k] = -g_0_z_z_xxyyyy[k] + g_z_0_z_xxyyyy[k] - g_z_z_z_xxyyyy[k] * ab_z + g_z_z_z_xxyyyyz[k];

                g_z_z_zz_xxyyyz[k] = -g_0_z_z_xxyyyz[k] + g_z_0_z_xxyyyz[k] - g_z_z_z_xxyyyz[k] * ab_z + g_z_z_z_xxyyyzz[k];

                g_z_z_zz_xxyyzz[k] = -g_0_z_z_xxyyzz[k] + g_z_0_z_xxyyzz[k] - g_z_z_z_xxyyzz[k] * ab_z + g_z_z_z_xxyyzzz[k];

                g_z_z_zz_xxyzzz[k] = -g_0_z_z_xxyzzz[k] + g_z_0_z_xxyzzz[k] - g_z_z_z_xxyzzz[k] * ab_z + g_z_z_z_xxyzzzz[k];

                g_z_z_zz_xxzzzz[k] = -g_0_z_z_xxzzzz[k] + g_z_0_z_xxzzzz[k] - g_z_z_z_xxzzzz[k] * ab_z + g_z_z_z_xxzzzzz[k];

                g_z_z_zz_xyyyyy[k] = -g_0_z_z_xyyyyy[k] + g_z_0_z_xyyyyy[k] - g_z_z_z_xyyyyy[k] * ab_z + g_z_z_z_xyyyyyz[k];

                g_z_z_zz_xyyyyz[k] = -g_0_z_z_xyyyyz[k] + g_z_0_z_xyyyyz[k] - g_z_z_z_xyyyyz[k] * ab_z + g_z_z_z_xyyyyzz[k];

                g_z_z_zz_xyyyzz[k] = -g_0_z_z_xyyyzz[k] + g_z_0_z_xyyyzz[k] - g_z_z_z_xyyyzz[k] * ab_z + g_z_z_z_xyyyzzz[k];

                g_z_z_zz_xyyzzz[k] = -g_0_z_z_xyyzzz[k] + g_z_0_z_xyyzzz[k] - g_z_z_z_xyyzzz[k] * ab_z + g_z_z_z_xyyzzzz[k];

                g_z_z_zz_xyzzzz[k] = -g_0_z_z_xyzzzz[k] + g_z_0_z_xyzzzz[k] - g_z_z_z_xyzzzz[k] * ab_z + g_z_z_z_xyzzzzz[k];

                g_z_z_zz_xzzzzz[k] = -g_0_z_z_xzzzzz[k] + g_z_0_z_xzzzzz[k] - g_z_z_z_xzzzzz[k] * ab_z + g_z_z_z_xzzzzzz[k];

                g_z_z_zz_yyyyyy[k] = -g_0_z_z_yyyyyy[k] + g_z_0_z_yyyyyy[k] - g_z_z_z_yyyyyy[k] * ab_z + g_z_z_z_yyyyyyz[k];

                g_z_z_zz_yyyyyz[k] = -g_0_z_z_yyyyyz[k] + g_z_0_z_yyyyyz[k] - g_z_z_z_yyyyyz[k] * ab_z + g_z_z_z_yyyyyzz[k];

                g_z_z_zz_yyyyzz[k] = -g_0_z_z_yyyyzz[k] + g_z_0_z_yyyyzz[k] - g_z_z_z_yyyyzz[k] * ab_z + g_z_z_z_yyyyzzz[k];

                g_z_z_zz_yyyzzz[k] = -g_0_z_z_yyyzzz[k] + g_z_0_z_yyyzzz[k] - g_z_z_z_yyyzzz[k] * ab_z + g_z_z_z_yyyzzzz[k];

                g_z_z_zz_yyzzzz[k] = -g_0_z_z_yyzzzz[k] + g_z_0_z_yyzzzz[k] - g_z_z_z_yyzzzz[k] * ab_z + g_z_z_z_yyzzzzz[k];

                g_z_z_zz_yzzzzz[k] = -g_0_z_z_yzzzzz[k] + g_z_0_z_yzzzzz[k] - g_z_z_z_yzzzzz[k] * ab_z + g_z_z_z_yzzzzzz[k];

                g_z_z_zz_zzzzzz[k] = -g_0_z_z_zzzzzz[k] + g_z_0_z_zzzzzz[k] - g_z_z_z_zzzzzz[k] * ab_z + g_z_z_z_zzzzzzz[k];
            }
        }
    }
}

} // erirec namespace

