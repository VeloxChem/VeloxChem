#include "ElectronRepulsionGeom2000ContrRecDHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_dhxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_dhxx,
                                            const size_t idx_geom_10_phxx,
                                            const size_t idx_geom_20_phxx,
                                            const size_t idx_geom_20_pixx,
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
            /// Set up components of auxilary buffer : PHSS

            const auto ph_geom_10_off = idx_geom_10_phxx + i * dcomps + j;

            auto g_x_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 89 * ccomps * dcomps);

            auto g_y_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 119 * ccomps * dcomps);

            auto g_y_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_x_xxxxx = cbuffer.data(ph_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_x_xxxxy = cbuffer.data(ph_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_x_xxxxz = cbuffer.data(ph_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_x_xxxyy = cbuffer.data(ph_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_x_xxxyz = cbuffer.data(ph_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_x_xxxzz = cbuffer.data(ph_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_x_xxyyy = cbuffer.data(ph_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_x_xxyyz = cbuffer.data(ph_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_x_xxyzz = cbuffer.data(ph_geom_10_off + 134 * ccomps * dcomps);

            auto g_z_0_x_xxzzz = cbuffer.data(ph_geom_10_off + 135 * ccomps * dcomps);

            auto g_z_0_x_xyyyy = cbuffer.data(ph_geom_10_off + 136 * ccomps * dcomps);

            auto g_z_0_x_xyyyz = cbuffer.data(ph_geom_10_off + 137 * ccomps * dcomps);

            auto g_z_0_x_xyyzz = cbuffer.data(ph_geom_10_off + 138 * ccomps * dcomps);

            auto g_z_0_x_xyzzz = cbuffer.data(ph_geom_10_off + 139 * ccomps * dcomps);

            auto g_z_0_x_xzzzz = cbuffer.data(ph_geom_10_off + 140 * ccomps * dcomps);

            auto g_z_0_x_yyyyy = cbuffer.data(ph_geom_10_off + 141 * ccomps * dcomps);

            auto g_z_0_x_yyyyz = cbuffer.data(ph_geom_10_off + 142 * ccomps * dcomps);

            auto g_z_0_x_yyyzz = cbuffer.data(ph_geom_10_off + 143 * ccomps * dcomps);

            auto g_z_0_x_yyzzz = cbuffer.data(ph_geom_10_off + 144 * ccomps * dcomps);

            auto g_z_0_x_yzzzz = cbuffer.data(ph_geom_10_off + 145 * ccomps * dcomps);

            auto g_z_0_x_zzzzz = cbuffer.data(ph_geom_10_off + 146 * ccomps * dcomps);

            auto g_z_0_y_xxxxx = cbuffer.data(ph_geom_10_off + 147 * ccomps * dcomps);

            auto g_z_0_y_xxxxy = cbuffer.data(ph_geom_10_off + 148 * ccomps * dcomps);

            auto g_z_0_y_xxxxz = cbuffer.data(ph_geom_10_off + 149 * ccomps * dcomps);

            auto g_z_0_y_xxxyy = cbuffer.data(ph_geom_10_off + 150 * ccomps * dcomps);

            auto g_z_0_y_xxxyz = cbuffer.data(ph_geom_10_off + 151 * ccomps * dcomps);

            auto g_z_0_y_xxxzz = cbuffer.data(ph_geom_10_off + 152 * ccomps * dcomps);

            auto g_z_0_y_xxyyy = cbuffer.data(ph_geom_10_off + 153 * ccomps * dcomps);

            auto g_z_0_y_xxyyz = cbuffer.data(ph_geom_10_off + 154 * ccomps * dcomps);

            auto g_z_0_y_xxyzz = cbuffer.data(ph_geom_10_off + 155 * ccomps * dcomps);

            auto g_z_0_y_xxzzz = cbuffer.data(ph_geom_10_off + 156 * ccomps * dcomps);

            auto g_z_0_y_xyyyy = cbuffer.data(ph_geom_10_off + 157 * ccomps * dcomps);

            auto g_z_0_y_xyyyz = cbuffer.data(ph_geom_10_off + 158 * ccomps * dcomps);

            auto g_z_0_y_xyyzz = cbuffer.data(ph_geom_10_off + 159 * ccomps * dcomps);

            auto g_z_0_y_xyzzz = cbuffer.data(ph_geom_10_off + 160 * ccomps * dcomps);

            auto g_z_0_y_xzzzz = cbuffer.data(ph_geom_10_off + 161 * ccomps * dcomps);

            auto g_z_0_y_yyyyy = cbuffer.data(ph_geom_10_off + 162 * ccomps * dcomps);

            auto g_z_0_y_yyyyz = cbuffer.data(ph_geom_10_off + 163 * ccomps * dcomps);

            auto g_z_0_y_yyyzz = cbuffer.data(ph_geom_10_off + 164 * ccomps * dcomps);

            auto g_z_0_y_yyzzz = cbuffer.data(ph_geom_10_off + 165 * ccomps * dcomps);

            auto g_z_0_y_yzzzz = cbuffer.data(ph_geom_10_off + 166 * ccomps * dcomps);

            auto g_z_0_y_zzzzz = cbuffer.data(ph_geom_10_off + 167 * ccomps * dcomps);

            auto g_z_0_z_xxxxx = cbuffer.data(ph_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_z_xxxxy = cbuffer.data(ph_geom_10_off + 169 * ccomps * dcomps);

            auto g_z_0_z_xxxxz = cbuffer.data(ph_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_z_xxxyy = cbuffer.data(ph_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_z_xxxyz = cbuffer.data(ph_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_z_xxxzz = cbuffer.data(ph_geom_10_off + 173 * ccomps * dcomps);

            auto g_z_0_z_xxyyy = cbuffer.data(ph_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_z_xxyyz = cbuffer.data(ph_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_z_xxyzz = cbuffer.data(ph_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_z_xxzzz = cbuffer.data(ph_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_z_xyyyy = cbuffer.data(ph_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_z_xyyyz = cbuffer.data(ph_geom_10_off + 179 * ccomps * dcomps);

            auto g_z_0_z_xyyzz = cbuffer.data(ph_geom_10_off + 180 * ccomps * dcomps);

            auto g_z_0_z_xyzzz = cbuffer.data(ph_geom_10_off + 181 * ccomps * dcomps);

            auto g_z_0_z_xzzzz = cbuffer.data(ph_geom_10_off + 182 * ccomps * dcomps);

            auto g_z_0_z_yyyyy = cbuffer.data(ph_geom_10_off + 183 * ccomps * dcomps);

            auto g_z_0_z_yyyyz = cbuffer.data(ph_geom_10_off + 184 * ccomps * dcomps);

            auto g_z_0_z_yyyzz = cbuffer.data(ph_geom_10_off + 185 * ccomps * dcomps);

            auto g_z_0_z_yyzzz = cbuffer.data(ph_geom_10_off + 186 * ccomps * dcomps);

            auto g_z_0_z_yzzzz = cbuffer.data(ph_geom_10_off + 187 * ccomps * dcomps);

            auto g_z_0_z_zzzzz = cbuffer.data(ph_geom_10_off + 188 * ccomps * dcomps);

            /// Set up components of auxilary buffer : PHSS

            const auto ph_geom_20_off = idx_geom_20_phxx + i * dcomps + j;

            auto g_xx_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 71 * ccomps * dcomps);

            auto g_xy_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 74 * ccomps * dcomps);

            auto g_xy_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 77 * ccomps * dcomps);

            auto g_xy_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 79 * ccomps * dcomps);

            auto g_xy_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 80 * ccomps * dcomps);

            auto g_xy_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 83 * ccomps * dcomps);

            auto g_xy_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 86 * ccomps * dcomps);

            auto g_xy_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 89 * ccomps * dcomps);

            auto g_xy_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 125 * ccomps * dcomps);

            auto g_xz_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 128 * ccomps * dcomps);

            auto g_xz_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 129 * ccomps * dcomps);

            auto g_xz_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 131 * ccomps * dcomps);

            auto g_xz_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 134 * ccomps * dcomps);

            auto g_xz_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 135 * ccomps * dcomps);

            auto g_xz_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 136 * ccomps * dcomps);

            auto g_xz_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 137 * ccomps * dcomps);

            auto g_xz_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 138 * ccomps * dcomps);

            auto g_xz_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 139 * ccomps * dcomps);

            auto g_xz_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 140 * ccomps * dcomps);

            auto g_xz_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 141 * ccomps * dcomps);

            auto g_xz_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 142 * ccomps * dcomps);

            auto g_xz_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 143 * ccomps * dcomps);

            auto g_xz_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 144 * ccomps * dcomps);

            auto g_xz_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 145 * ccomps * dcomps);

            auto g_xz_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 146 * ccomps * dcomps);

            auto g_xz_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 147 * ccomps * dcomps);

            auto g_xz_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 148 * ccomps * dcomps);

            auto g_xz_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 149 * ccomps * dcomps);

            auto g_xz_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 150 * ccomps * dcomps);

            auto g_xz_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 151 * ccomps * dcomps);

            auto g_xz_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 152 * ccomps * dcomps);

            auto g_xz_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 153 * ccomps * dcomps);

            auto g_xz_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 154 * ccomps * dcomps);

            auto g_xz_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 155 * ccomps * dcomps);

            auto g_xz_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 156 * ccomps * dcomps);

            auto g_xz_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 157 * ccomps * dcomps);

            auto g_xz_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 158 * ccomps * dcomps);

            auto g_xz_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 159 * ccomps * dcomps);

            auto g_xz_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 160 * ccomps * dcomps);

            auto g_xz_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 161 * ccomps * dcomps);

            auto g_xz_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 162 * ccomps * dcomps);

            auto g_xz_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 163 * ccomps * dcomps);

            auto g_xz_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 164 * ccomps * dcomps);

            auto g_xz_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 165 * ccomps * dcomps);

            auto g_xz_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 166 * ccomps * dcomps);

            auto g_xz_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 167 * ccomps * dcomps);

            auto g_xz_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 168 * ccomps * dcomps);

            auto g_xz_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 169 * ccomps * dcomps);

            auto g_xz_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 170 * ccomps * dcomps);

            auto g_xz_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 171 * ccomps * dcomps);

            auto g_xz_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 172 * ccomps * dcomps);

            auto g_xz_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 173 * ccomps * dcomps);

            auto g_xz_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 174 * ccomps * dcomps);

            auto g_xz_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 175 * ccomps * dcomps);

            auto g_xz_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 176 * ccomps * dcomps);

            auto g_xz_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 177 * ccomps * dcomps);

            auto g_xz_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 178 * ccomps * dcomps);

            auto g_xz_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 179 * ccomps * dcomps);

            auto g_xz_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 182 * ccomps * dcomps);

            auto g_xz_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 185 * ccomps * dcomps);

            auto g_xz_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 188 * ccomps * dcomps);

            auto g_yy_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 189 * ccomps * dcomps);

            auto g_yy_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 190 * ccomps * dcomps);

            auto g_yy_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 191 * ccomps * dcomps);

            auto g_yy_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 192 * ccomps * dcomps);

            auto g_yy_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 193 * ccomps * dcomps);

            auto g_yy_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 194 * ccomps * dcomps);

            auto g_yy_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 195 * ccomps * dcomps);

            auto g_yy_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 196 * ccomps * dcomps);

            auto g_yy_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 197 * ccomps * dcomps);

            auto g_yy_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 198 * ccomps * dcomps);

            auto g_yy_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 199 * ccomps * dcomps);

            auto g_yy_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 200 * ccomps * dcomps);

            auto g_yy_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 201 * ccomps * dcomps);

            auto g_yy_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 202 * ccomps * dcomps);

            auto g_yy_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 203 * ccomps * dcomps);

            auto g_yy_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 204 * ccomps * dcomps);

            auto g_yy_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 205 * ccomps * dcomps);

            auto g_yy_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 206 * ccomps * dcomps);

            auto g_yy_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 207 * ccomps * dcomps);

            auto g_yy_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 208 * ccomps * dcomps);

            auto g_yy_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 209 * ccomps * dcomps);

            auto g_yy_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 210 * ccomps * dcomps);

            auto g_yy_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 211 * ccomps * dcomps);

            auto g_yy_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 212 * ccomps * dcomps);

            auto g_yy_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 213 * ccomps * dcomps);

            auto g_yy_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 214 * ccomps * dcomps);

            auto g_yy_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 215 * ccomps * dcomps);

            auto g_yy_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 216 * ccomps * dcomps);

            auto g_yy_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 217 * ccomps * dcomps);

            auto g_yy_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 218 * ccomps * dcomps);

            auto g_yy_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 219 * ccomps * dcomps);

            auto g_yy_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 220 * ccomps * dcomps);

            auto g_yy_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 221 * ccomps * dcomps);

            auto g_yy_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 222 * ccomps * dcomps);

            auto g_yy_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 223 * ccomps * dcomps);

            auto g_yy_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 224 * ccomps * dcomps);

            auto g_yy_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 225 * ccomps * dcomps);

            auto g_yy_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 226 * ccomps * dcomps);

            auto g_yy_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 227 * ccomps * dcomps);

            auto g_yy_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 228 * ccomps * dcomps);

            auto g_yy_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 229 * ccomps * dcomps);

            auto g_yy_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 230 * ccomps * dcomps);

            auto g_yy_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 231 * ccomps * dcomps);

            auto g_yy_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 232 * ccomps * dcomps);

            auto g_yy_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 233 * ccomps * dcomps);

            auto g_yy_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 234 * ccomps * dcomps);

            auto g_yy_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 235 * ccomps * dcomps);

            auto g_yy_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 236 * ccomps * dcomps);

            auto g_yy_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 237 * ccomps * dcomps);

            auto g_yy_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 238 * ccomps * dcomps);

            auto g_yy_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 239 * ccomps * dcomps);

            auto g_yy_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 240 * ccomps * dcomps);

            auto g_yy_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 241 * ccomps * dcomps);

            auto g_yy_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 242 * ccomps * dcomps);

            auto g_yy_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 243 * ccomps * dcomps);

            auto g_yy_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 244 * ccomps * dcomps);

            auto g_yy_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 245 * ccomps * dcomps);

            auto g_yy_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 246 * ccomps * dcomps);

            auto g_yy_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 247 * ccomps * dcomps);

            auto g_yy_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 248 * ccomps * dcomps);

            auto g_yy_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 249 * ccomps * dcomps);

            auto g_yy_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 250 * ccomps * dcomps);

            auto g_yy_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 251 * ccomps * dcomps);

            auto g_yz_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 252 * ccomps * dcomps);

            auto g_yz_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 253 * ccomps * dcomps);

            auto g_yz_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 254 * ccomps * dcomps);

            auto g_yz_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 255 * ccomps * dcomps);

            auto g_yz_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 256 * ccomps * dcomps);

            auto g_yz_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 257 * ccomps * dcomps);

            auto g_yz_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 258 * ccomps * dcomps);

            auto g_yz_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 259 * ccomps * dcomps);

            auto g_yz_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 260 * ccomps * dcomps);

            auto g_yz_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 261 * ccomps * dcomps);

            auto g_yz_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 262 * ccomps * dcomps);

            auto g_yz_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 263 * ccomps * dcomps);

            auto g_yz_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 264 * ccomps * dcomps);

            auto g_yz_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 265 * ccomps * dcomps);

            auto g_yz_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 266 * ccomps * dcomps);

            auto g_yz_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 267 * ccomps * dcomps);

            auto g_yz_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 268 * ccomps * dcomps);

            auto g_yz_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 269 * ccomps * dcomps);

            auto g_yz_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 270 * ccomps * dcomps);

            auto g_yz_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 271 * ccomps * dcomps);

            auto g_yz_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 272 * ccomps * dcomps);

            auto g_yz_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 273 * ccomps * dcomps);

            auto g_yz_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 274 * ccomps * dcomps);

            auto g_yz_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 275 * ccomps * dcomps);

            auto g_yz_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 276 * ccomps * dcomps);

            auto g_yz_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 277 * ccomps * dcomps);

            auto g_yz_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 278 * ccomps * dcomps);

            auto g_yz_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 279 * ccomps * dcomps);

            auto g_yz_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 280 * ccomps * dcomps);

            auto g_yz_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 281 * ccomps * dcomps);

            auto g_yz_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 282 * ccomps * dcomps);

            auto g_yz_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 283 * ccomps * dcomps);

            auto g_yz_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 284 * ccomps * dcomps);

            auto g_yz_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 285 * ccomps * dcomps);

            auto g_yz_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 286 * ccomps * dcomps);

            auto g_yz_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 287 * ccomps * dcomps);

            auto g_yz_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 288 * ccomps * dcomps);

            auto g_yz_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 289 * ccomps * dcomps);

            auto g_yz_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 290 * ccomps * dcomps);

            auto g_yz_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 291 * ccomps * dcomps);

            auto g_yz_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 292 * ccomps * dcomps);

            auto g_yz_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 293 * ccomps * dcomps);

            auto g_yz_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 294 * ccomps * dcomps);

            auto g_yz_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 295 * ccomps * dcomps);

            auto g_yz_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 296 * ccomps * dcomps);

            auto g_yz_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 297 * ccomps * dcomps);

            auto g_yz_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 298 * ccomps * dcomps);

            auto g_yz_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 299 * ccomps * dcomps);

            auto g_yz_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 300 * ccomps * dcomps);

            auto g_yz_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 301 * ccomps * dcomps);

            auto g_yz_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 302 * ccomps * dcomps);

            auto g_yz_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 303 * ccomps * dcomps);

            auto g_yz_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 304 * ccomps * dcomps);

            auto g_yz_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 305 * ccomps * dcomps);

            auto g_yz_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 306 * ccomps * dcomps);

            auto g_yz_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 307 * ccomps * dcomps);

            auto g_yz_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 308 * ccomps * dcomps);

            auto g_yz_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 309 * ccomps * dcomps);

            auto g_yz_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 310 * ccomps * dcomps);

            auto g_yz_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 311 * ccomps * dcomps);

            auto g_yz_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 312 * ccomps * dcomps);

            auto g_yz_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 313 * ccomps * dcomps);

            auto g_yz_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 314 * ccomps * dcomps);

            auto g_zz_0_x_xxxxx = cbuffer.data(ph_geom_20_off + 315 * ccomps * dcomps);

            auto g_zz_0_x_xxxxy = cbuffer.data(ph_geom_20_off + 316 * ccomps * dcomps);

            auto g_zz_0_x_xxxxz = cbuffer.data(ph_geom_20_off + 317 * ccomps * dcomps);

            auto g_zz_0_x_xxxyy = cbuffer.data(ph_geom_20_off + 318 * ccomps * dcomps);

            auto g_zz_0_x_xxxyz = cbuffer.data(ph_geom_20_off + 319 * ccomps * dcomps);

            auto g_zz_0_x_xxxzz = cbuffer.data(ph_geom_20_off + 320 * ccomps * dcomps);

            auto g_zz_0_x_xxyyy = cbuffer.data(ph_geom_20_off + 321 * ccomps * dcomps);

            auto g_zz_0_x_xxyyz = cbuffer.data(ph_geom_20_off + 322 * ccomps * dcomps);

            auto g_zz_0_x_xxyzz = cbuffer.data(ph_geom_20_off + 323 * ccomps * dcomps);

            auto g_zz_0_x_xxzzz = cbuffer.data(ph_geom_20_off + 324 * ccomps * dcomps);

            auto g_zz_0_x_xyyyy = cbuffer.data(ph_geom_20_off + 325 * ccomps * dcomps);

            auto g_zz_0_x_xyyyz = cbuffer.data(ph_geom_20_off + 326 * ccomps * dcomps);

            auto g_zz_0_x_xyyzz = cbuffer.data(ph_geom_20_off + 327 * ccomps * dcomps);

            auto g_zz_0_x_xyzzz = cbuffer.data(ph_geom_20_off + 328 * ccomps * dcomps);

            auto g_zz_0_x_xzzzz = cbuffer.data(ph_geom_20_off + 329 * ccomps * dcomps);

            auto g_zz_0_x_yyyyy = cbuffer.data(ph_geom_20_off + 330 * ccomps * dcomps);

            auto g_zz_0_x_yyyyz = cbuffer.data(ph_geom_20_off + 331 * ccomps * dcomps);

            auto g_zz_0_x_yyyzz = cbuffer.data(ph_geom_20_off + 332 * ccomps * dcomps);

            auto g_zz_0_x_yyzzz = cbuffer.data(ph_geom_20_off + 333 * ccomps * dcomps);

            auto g_zz_0_x_yzzzz = cbuffer.data(ph_geom_20_off + 334 * ccomps * dcomps);

            auto g_zz_0_x_zzzzz = cbuffer.data(ph_geom_20_off + 335 * ccomps * dcomps);

            auto g_zz_0_y_xxxxx = cbuffer.data(ph_geom_20_off + 336 * ccomps * dcomps);

            auto g_zz_0_y_xxxxy = cbuffer.data(ph_geom_20_off + 337 * ccomps * dcomps);

            auto g_zz_0_y_xxxxz = cbuffer.data(ph_geom_20_off + 338 * ccomps * dcomps);

            auto g_zz_0_y_xxxyy = cbuffer.data(ph_geom_20_off + 339 * ccomps * dcomps);

            auto g_zz_0_y_xxxyz = cbuffer.data(ph_geom_20_off + 340 * ccomps * dcomps);

            auto g_zz_0_y_xxxzz = cbuffer.data(ph_geom_20_off + 341 * ccomps * dcomps);

            auto g_zz_0_y_xxyyy = cbuffer.data(ph_geom_20_off + 342 * ccomps * dcomps);

            auto g_zz_0_y_xxyyz = cbuffer.data(ph_geom_20_off + 343 * ccomps * dcomps);

            auto g_zz_0_y_xxyzz = cbuffer.data(ph_geom_20_off + 344 * ccomps * dcomps);

            auto g_zz_0_y_xxzzz = cbuffer.data(ph_geom_20_off + 345 * ccomps * dcomps);

            auto g_zz_0_y_xyyyy = cbuffer.data(ph_geom_20_off + 346 * ccomps * dcomps);

            auto g_zz_0_y_xyyyz = cbuffer.data(ph_geom_20_off + 347 * ccomps * dcomps);

            auto g_zz_0_y_xyyzz = cbuffer.data(ph_geom_20_off + 348 * ccomps * dcomps);

            auto g_zz_0_y_xyzzz = cbuffer.data(ph_geom_20_off + 349 * ccomps * dcomps);

            auto g_zz_0_y_xzzzz = cbuffer.data(ph_geom_20_off + 350 * ccomps * dcomps);

            auto g_zz_0_y_yyyyy = cbuffer.data(ph_geom_20_off + 351 * ccomps * dcomps);

            auto g_zz_0_y_yyyyz = cbuffer.data(ph_geom_20_off + 352 * ccomps * dcomps);

            auto g_zz_0_y_yyyzz = cbuffer.data(ph_geom_20_off + 353 * ccomps * dcomps);

            auto g_zz_0_y_yyzzz = cbuffer.data(ph_geom_20_off + 354 * ccomps * dcomps);

            auto g_zz_0_y_yzzzz = cbuffer.data(ph_geom_20_off + 355 * ccomps * dcomps);

            auto g_zz_0_y_zzzzz = cbuffer.data(ph_geom_20_off + 356 * ccomps * dcomps);

            auto g_zz_0_z_xxxxx = cbuffer.data(ph_geom_20_off + 357 * ccomps * dcomps);

            auto g_zz_0_z_xxxxy = cbuffer.data(ph_geom_20_off + 358 * ccomps * dcomps);

            auto g_zz_0_z_xxxxz = cbuffer.data(ph_geom_20_off + 359 * ccomps * dcomps);

            auto g_zz_0_z_xxxyy = cbuffer.data(ph_geom_20_off + 360 * ccomps * dcomps);

            auto g_zz_0_z_xxxyz = cbuffer.data(ph_geom_20_off + 361 * ccomps * dcomps);

            auto g_zz_0_z_xxxzz = cbuffer.data(ph_geom_20_off + 362 * ccomps * dcomps);

            auto g_zz_0_z_xxyyy = cbuffer.data(ph_geom_20_off + 363 * ccomps * dcomps);

            auto g_zz_0_z_xxyyz = cbuffer.data(ph_geom_20_off + 364 * ccomps * dcomps);

            auto g_zz_0_z_xxyzz = cbuffer.data(ph_geom_20_off + 365 * ccomps * dcomps);

            auto g_zz_0_z_xxzzz = cbuffer.data(ph_geom_20_off + 366 * ccomps * dcomps);

            auto g_zz_0_z_xyyyy = cbuffer.data(ph_geom_20_off + 367 * ccomps * dcomps);

            auto g_zz_0_z_xyyyz = cbuffer.data(ph_geom_20_off + 368 * ccomps * dcomps);

            auto g_zz_0_z_xyyzz = cbuffer.data(ph_geom_20_off + 369 * ccomps * dcomps);

            auto g_zz_0_z_xyzzz = cbuffer.data(ph_geom_20_off + 370 * ccomps * dcomps);

            auto g_zz_0_z_xzzzz = cbuffer.data(ph_geom_20_off + 371 * ccomps * dcomps);

            auto g_zz_0_z_yyyyy = cbuffer.data(ph_geom_20_off + 372 * ccomps * dcomps);

            auto g_zz_0_z_yyyyz = cbuffer.data(ph_geom_20_off + 373 * ccomps * dcomps);

            auto g_zz_0_z_yyyzz = cbuffer.data(ph_geom_20_off + 374 * ccomps * dcomps);

            auto g_zz_0_z_yyzzz = cbuffer.data(ph_geom_20_off + 375 * ccomps * dcomps);

            auto g_zz_0_z_yzzzz = cbuffer.data(ph_geom_20_off + 376 * ccomps * dcomps);

            auto g_zz_0_z_zzzzz = cbuffer.data(ph_geom_20_off + 377 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_dhxx

            const auto dh_geom_20_off = idx_geom_20_dhxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_x_xxxxx, g_x_0_x_xxxxy, g_x_0_x_xxxxz, g_x_0_x_xxxyy, g_x_0_x_xxxyz, g_x_0_x_xxxzz, g_x_0_x_xxyyy, g_x_0_x_xxyyz, g_x_0_x_xxyzz, g_x_0_x_xxzzz, g_x_0_x_xyyyy, g_x_0_x_xyyyz, g_x_0_x_xyyzz, g_x_0_x_xyzzz, g_x_0_x_xzzzz, g_x_0_x_yyyyy, g_x_0_x_yyyyz, g_x_0_x_yyyzz, g_x_0_x_yyzzz, g_x_0_x_yzzzz, g_x_0_x_zzzzz, g_xx_0_x_xxxxx, g_xx_0_x_xxxxxx, g_xx_0_x_xxxxxy, g_xx_0_x_xxxxxz, g_xx_0_x_xxxxy, g_xx_0_x_xxxxyy, g_xx_0_x_xxxxyz, g_xx_0_x_xxxxz, g_xx_0_x_xxxxzz, g_xx_0_x_xxxyy, g_xx_0_x_xxxyyy, g_xx_0_x_xxxyyz, g_xx_0_x_xxxyz, g_xx_0_x_xxxyzz, g_xx_0_x_xxxzz, g_xx_0_x_xxxzzz, g_xx_0_x_xxyyy, g_xx_0_x_xxyyyy, g_xx_0_x_xxyyyz, g_xx_0_x_xxyyz, g_xx_0_x_xxyyzz, g_xx_0_x_xxyzz, g_xx_0_x_xxyzzz, g_xx_0_x_xxzzz, g_xx_0_x_xxzzzz, g_xx_0_x_xyyyy, g_xx_0_x_xyyyyy, g_xx_0_x_xyyyyz, g_xx_0_x_xyyyz, g_xx_0_x_xyyyzz, g_xx_0_x_xyyzz, g_xx_0_x_xyyzzz, g_xx_0_x_xyzzz, g_xx_0_x_xyzzzz, g_xx_0_x_xzzzz, g_xx_0_x_xzzzzz, g_xx_0_x_yyyyy, g_xx_0_x_yyyyz, g_xx_0_x_yyyzz, g_xx_0_x_yyzzz, g_xx_0_x_yzzzz, g_xx_0_x_zzzzz, g_xx_0_xx_xxxxx, g_xx_0_xx_xxxxy, g_xx_0_xx_xxxxz, g_xx_0_xx_xxxyy, g_xx_0_xx_xxxyz, g_xx_0_xx_xxxzz, g_xx_0_xx_xxyyy, g_xx_0_xx_xxyyz, g_xx_0_xx_xxyzz, g_xx_0_xx_xxzzz, g_xx_0_xx_xyyyy, g_xx_0_xx_xyyyz, g_xx_0_xx_xyyzz, g_xx_0_xx_xyzzz, g_xx_0_xx_xzzzz, g_xx_0_xx_yyyyy, g_xx_0_xx_yyyyz, g_xx_0_xx_yyyzz, g_xx_0_xx_yyzzz, g_xx_0_xx_yzzzz, g_xx_0_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xx_xxxxx[k] = -2.0 * g_x_0_x_xxxxx[k] - g_xx_0_x_xxxxx[k] * ab_x + g_xx_0_x_xxxxxx[k];

                g_xx_0_xx_xxxxy[k] = -2.0 * g_x_0_x_xxxxy[k] - g_xx_0_x_xxxxy[k] * ab_x + g_xx_0_x_xxxxxy[k];

                g_xx_0_xx_xxxxz[k] = -2.0 * g_x_0_x_xxxxz[k] - g_xx_0_x_xxxxz[k] * ab_x + g_xx_0_x_xxxxxz[k];

                g_xx_0_xx_xxxyy[k] = -2.0 * g_x_0_x_xxxyy[k] - g_xx_0_x_xxxyy[k] * ab_x + g_xx_0_x_xxxxyy[k];

                g_xx_0_xx_xxxyz[k] = -2.0 * g_x_0_x_xxxyz[k] - g_xx_0_x_xxxyz[k] * ab_x + g_xx_0_x_xxxxyz[k];

                g_xx_0_xx_xxxzz[k] = -2.0 * g_x_0_x_xxxzz[k] - g_xx_0_x_xxxzz[k] * ab_x + g_xx_0_x_xxxxzz[k];

                g_xx_0_xx_xxyyy[k] = -2.0 * g_x_0_x_xxyyy[k] - g_xx_0_x_xxyyy[k] * ab_x + g_xx_0_x_xxxyyy[k];

                g_xx_0_xx_xxyyz[k] = -2.0 * g_x_0_x_xxyyz[k] - g_xx_0_x_xxyyz[k] * ab_x + g_xx_0_x_xxxyyz[k];

                g_xx_0_xx_xxyzz[k] = -2.0 * g_x_0_x_xxyzz[k] - g_xx_0_x_xxyzz[k] * ab_x + g_xx_0_x_xxxyzz[k];

                g_xx_0_xx_xxzzz[k] = -2.0 * g_x_0_x_xxzzz[k] - g_xx_0_x_xxzzz[k] * ab_x + g_xx_0_x_xxxzzz[k];

                g_xx_0_xx_xyyyy[k] = -2.0 * g_x_0_x_xyyyy[k] - g_xx_0_x_xyyyy[k] * ab_x + g_xx_0_x_xxyyyy[k];

                g_xx_0_xx_xyyyz[k] = -2.0 * g_x_0_x_xyyyz[k] - g_xx_0_x_xyyyz[k] * ab_x + g_xx_0_x_xxyyyz[k];

                g_xx_0_xx_xyyzz[k] = -2.0 * g_x_0_x_xyyzz[k] - g_xx_0_x_xyyzz[k] * ab_x + g_xx_0_x_xxyyzz[k];

                g_xx_0_xx_xyzzz[k] = -2.0 * g_x_0_x_xyzzz[k] - g_xx_0_x_xyzzz[k] * ab_x + g_xx_0_x_xxyzzz[k];

                g_xx_0_xx_xzzzz[k] = -2.0 * g_x_0_x_xzzzz[k] - g_xx_0_x_xzzzz[k] * ab_x + g_xx_0_x_xxzzzz[k];

                g_xx_0_xx_yyyyy[k] = -2.0 * g_x_0_x_yyyyy[k] - g_xx_0_x_yyyyy[k] * ab_x + g_xx_0_x_xyyyyy[k];

                g_xx_0_xx_yyyyz[k] = -2.0 * g_x_0_x_yyyyz[k] - g_xx_0_x_yyyyz[k] * ab_x + g_xx_0_x_xyyyyz[k];

                g_xx_0_xx_yyyzz[k] = -2.0 * g_x_0_x_yyyzz[k] - g_xx_0_x_yyyzz[k] * ab_x + g_xx_0_x_xyyyzz[k];

                g_xx_0_xx_yyzzz[k] = -2.0 * g_x_0_x_yyzzz[k] - g_xx_0_x_yyzzz[k] * ab_x + g_xx_0_x_xyyzzz[k];

                g_xx_0_xx_yzzzz[k] = -2.0 * g_x_0_x_yzzzz[k] - g_xx_0_x_yzzzz[k] * ab_x + g_xx_0_x_xyzzzz[k];

                g_xx_0_xx_zzzzz[k] = -2.0 * g_x_0_x_zzzzz[k] - g_xx_0_x_zzzzz[k] * ab_x + g_xx_0_x_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_x_xxxxx, g_xx_0_x_xxxxxy, g_xx_0_x_xxxxy, g_xx_0_x_xxxxyy, g_xx_0_x_xxxxyz, g_xx_0_x_xxxxz, g_xx_0_x_xxxyy, g_xx_0_x_xxxyyy, g_xx_0_x_xxxyyz, g_xx_0_x_xxxyz, g_xx_0_x_xxxyzz, g_xx_0_x_xxxzz, g_xx_0_x_xxyyy, g_xx_0_x_xxyyyy, g_xx_0_x_xxyyyz, g_xx_0_x_xxyyz, g_xx_0_x_xxyyzz, g_xx_0_x_xxyzz, g_xx_0_x_xxyzzz, g_xx_0_x_xxzzz, g_xx_0_x_xyyyy, g_xx_0_x_xyyyyy, g_xx_0_x_xyyyyz, g_xx_0_x_xyyyz, g_xx_0_x_xyyyzz, g_xx_0_x_xyyzz, g_xx_0_x_xyyzzz, g_xx_0_x_xyzzz, g_xx_0_x_xyzzzz, g_xx_0_x_xzzzz, g_xx_0_x_yyyyy, g_xx_0_x_yyyyyy, g_xx_0_x_yyyyyz, g_xx_0_x_yyyyz, g_xx_0_x_yyyyzz, g_xx_0_x_yyyzz, g_xx_0_x_yyyzzz, g_xx_0_x_yyzzz, g_xx_0_x_yyzzzz, g_xx_0_x_yzzzz, g_xx_0_x_yzzzzz, g_xx_0_x_zzzzz, g_xx_0_xy_xxxxx, g_xx_0_xy_xxxxy, g_xx_0_xy_xxxxz, g_xx_0_xy_xxxyy, g_xx_0_xy_xxxyz, g_xx_0_xy_xxxzz, g_xx_0_xy_xxyyy, g_xx_0_xy_xxyyz, g_xx_0_xy_xxyzz, g_xx_0_xy_xxzzz, g_xx_0_xy_xyyyy, g_xx_0_xy_xyyyz, g_xx_0_xy_xyyzz, g_xx_0_xy_xyzzz, g_xx_0_xy_xzzzz, g_xx_0_xy_yyyyy, g_xx_0_xy_yyyyz, g_xx_0_xy_yyyzz, g_xx_0_xy_yyzzz, g_xx_0_xy_yzzzz, g_xx_0_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xy_xxxxx[k] = -g_xx_0_x_xxxxx[k] * ab_y + g_xx_0_x_xxxxxy[k];

                g_xx_0_xy_xxxxy[k] = -g_xx_0_x_xxxxy[k] * ab_y + g_xx_0_x_xxxxyy[k];

                g_xx_0_xy_xxxxz[k] = -g_xx_0_x_xxxxz[k] * ab_y + g_xx_0_x_xxxxyz[k];

                g_xx_0_xy_xxxyy[k] = -g_xx_0_x_xxxyy[k] * ab_y + g_xx_0_x_xxxyyy[k];

                g_xx_0_xy_xxxyz[k] = -g_xx_0_x_xxxyz[k] * ab_y + g_xx_0_x_xxxyyz[k];

                g_xx_0_xy_xxxzz[k] = -g_xx_0_x_xxxzz[k] * ab_y + g_xx_0_x_xxxyzz[k];

                g_xx_0_xy_xxyyy[k] = -g_xx_0_x_xxyyy[k] * ab_y + g_xx_0_x_xxyyyy[k];

                g_xx_0_xy_xxyyz[k] = -g_xx_0_x_xxyyz[k] * ab_y + g_xx_0_x_xxyyyz[k];

                g_xx_0_xy_xxyzz[k] = -g_xx_0_x_xxyzz[k] * ab_y + g_xx_0_x_xxyyzz[k];

                g_xx_0_xy_xxzzz[k] = -g_xx_0_x_xxzzz[k] * ab_y + g_xx_0_x_xxyzzz[k];

                g_xx_0_xy_xyyyy[k] = -g_xx_0_x_xyyyy[k] * ab_y + g_xx_0_x_xyyyyy[k];

                g_xx_0_xy_xyyyz[k] = -g_xx_0_x_xyyyz[k] * ab_y + g_xx_0_x_xyyyyz[k];

                g_xx_0_xy_xyyzz[k] = -g_xx_0_x_xyyzz[k] * ab_y + g_xx_0_x_xyyyzz[k];

                g_xx_0_xy_xyzzz[k] = -g_xx_0_x_xyzzz[k] * ab_y + g_xx_0_x_xyyzzz[k];

                g_xx_0_xy_xzzzz[k] = -g_xx_0_x_xzzzz[k] * ab_y + g_xx_0_x_xyzzzz[k];

                g_xx_0_xy_yyyyy[k] = -g_xx_0_x_yyyyy[k] * ab_y + g_xx_0_x_yyyyyy[k];

                g_xx_0_xy_yyyyz[k] = -g_xx_0_x_yyyyz[k] * ab_y + g_xx_0_x_yyyyyz[k];

                g_xx_0_xy_yyyzz[k] = -g_xx_0_x_yyyzz[k] * ab_y + g_xx_0_x_yyyyzz[k];

                g_xx_0_xy_yyzzz[k] = -g_xx_0_x_yyzzz[k] * ab_y + g_xx_0_x_yyyzzz[k];

                g_xx_0_xy_yzzzz[k] = -g_xx_0_x_yzzzz[k] * ab_y + g_xx_0_x_yyzzzz[k];

                g_xx_0_xy_zzzzz[k] = -g_xx_0_x_zzzzz[k] * ab_y + g_xx_0_x_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_x_xxxxx, g_xx_0_x_xxxxxz, g_xx_0_x_xxxxy, g_xx_0_x_xxxxyz, g_xx_0_x_xxxxz, g_xx_0_x_xxxxzz, g_xx_0_x_xxxyy, g_xx_0_x_xxxyyz, g_xx_0_x_xxxyz, g_xx_0_x_xxxyzz, g_xx_0_x_xxxzz, g_xx_0_x_xxxzzz, g_xx_0_x_xxyyy, g_xx_0_x_xxyyyz, g_xx_0_x_xxyyz, g_xx_0_x_xxyyzz, g_xx_0_x_xxyzz, g_xx_0_x_xxyzzz, g_xx_0_x_xxzzz, g_xx_0_x_xxzzzz, g_xx_0_x_xyyyy, g_xx_0_x_xyyyyz, g_xx_0_x_xyyyz, g_xx_0_x_xyyyzz, g_xx_0_x_xyyzz, g_xx_0_x_xyyzzz, g_xx_0_x_xyzzz, g_xx_0_x_xyzzzz, g_xx_0_x_xzzzz, g_xx_0_x_xzzzzz, g_xx_0_x_yyyyy, g_xx_0_x_yyyyyz, g_xx_0_x_yyyyz, g_xx_0_x_yyyyzz, g_xx_0_x_yyyzz, g_xx_0_x_yyyzzz, g_xx_0_x_yyzzz, g_xx_0_x_yyzzzz, g_xx_0_x_yzzzz, g_xx_0_x_yzzzzz, g_xx_0_x_zzzzz, g_xx_0_x_zzzzzz, g_xx_0_xz_xxxxx, g_xx_0_xz_xxxxy, g_xx_0_xz_xxxxz, g_xx_0_xz_xxxyy, g_xx_0_xz_xxxyz, g_xx_0_xz_xxxzz, g_xx_0_xz_xxyyy, g_xx_0_xz_xxyyz, g_xx_0_xz_xxyzz, g_xx_0_xz_xxzzz, g_xx_0_xz_xyyyy, g_xx_0_xz_xyyyz, g_xx_0_xz_xyyzz, g_xx_0_xz_xyzzz, g_xx_0_xz_xzzzz, g_xx_0_xz_yyyyy, g_xx_0_xz_yyyyz, g_xx_0_xz_yyyzz, g_xx_0_xz_yyzzz, g_xx_0_xz_yzzzz, g_xx_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xz_xxxxx[k] = -g_xx_0_x_xxxxx[k] * ab_z + g_xx_0_x_xxxxxz[k];

                g_xx_0_xz_xxxxy[k] = -g_xx_0_x_xxxxy[k] * ab_z + g_xx_0_x_xxxxyz[k];

                g_xx_0_xz_xxxxz[k] = -g_xx_0_x_xxxxz[k] * ab_z + g_xx_0_x_xxxxzz[k];

                g_xx_0_xz_xxxyy[k] = -g_xx_0_x_xxxyy[k] * ab_z + g_xx_0_x_xxxyyz[k];

                g_xx_0_xz_xxxyz[k] = -g_xx_0_x_xxxyz[k] * ab_z + g_xx_0_x_xxxyzz[k];

                g_xx_0_xz_xxxzz[k] = -g_xx_0_x_xxxzz[k] * ab_z + g_xx_0_x_xxxzzz[k];

                g_xx_0_xz_xxyyy[k] = -g_xx_0_x_xxyyy[k] * ab_z + g_xx_0_x_xxyyyz[k];

                g_xx_0_xz_xxyyz[k] = -g_xx_0_x_xxyyz[k] * ab_z + g_xx_0_x_xxyyzz[k];

                g_xx_0_xz_xxyzz[k] = -g_xx_0_x_xxyzz[k] * ab_z + g_xx_0_x_xxyzzz[k];

                g_xx_0_xz_xxzzz[k] = -g_xx_0_x_xxzzz[k] * ab_z + g_xx_0_x_xxzzzz[k];

                g_xx_0_xz_xyyyy[k] = -g_xx_0_x_xyyyy[k] * ab_z + g_xx_0_x_xyyyyz[k];

                g_xx_0_xz_xyyyz[k] = -g_xx_0_x_xyyyz[k] * ab_z + g_xx_0_x_xyyyzz[k];

                g_xx_0_xz_xyyzz[k] = -g_xx_0_x_xyyzz[k] * ab_z + g_xx_0_x_xyyzzz[k];

                g_xx_0_xz_xyzzz[k] = -g_xx_0_x_xyzzz[k] * ab_z + g_xx_0_x_xyzzzz[k];

                g_xx_0_xz_xzzzz[k] = -g_xx_0_x_xzzzz[k] * ab_z + g_xx_0_x_xzzzzz[k];

                g_xx_0_xz_yyyyy[k] = -g_xx_0_x_yyyyy[k] * ab_z + g_xx_0_x_yyyyyz[k];

                g_xx_0_xz_yyyyz[k] = -g_xx_0_x_yyyyz[k] * ab_z + g_xx_0_x_yyyyzz[k];

                g_xx_0_xz_yyyzz[k] = -g_xx_0_x_yyyzz[k] * ab_z + g_xx_0_x_yyyzzz[k];

                g_xx_0_xz_yyzzz[k] = -g_xx_0_x_yyzzz[k] * ab_z + g_xx_0_x_yyzzzz[k];

                g_xx_0_xz_yzzzz[k] = -g_xx_0_x_yzzzz[k] * ab_z + g_xx_0_x_yzzzzz[k];

                g_xx_0_xz_zzzzz[k] = -g_xx_0_x_zzzzz[k] * ab_z + g_xx_0_x_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_y_xxxxx, g_xx_0_y_xxxxxy, g_xx_0_y_xxxxy, g_xx_0_y_xxxxyy, g_xx_0_y_xxxxyz, g_xx_0_y_xxxxz, g_xx_0_y_xxxyy, g_xx_0_y_xxxyyy, g_xx_0_y_xxxyyz, g_xx_0_y_xxxyz, g_xx_0_y_xxxyzz, g_xx_0_y_xxxzz, g_xx_0_y_xxyyy, g_xx_0_y_xxyyyy, g_xx_0_y_xxyyyz, g_xx_0_y_xxyyz, g_xx_0_y_xxyyzz, g_xx_0_y_xxyzz, g_xx_0_y_xxyzzz, g_xx_0_y_xxzzz, g_xx_0_y_xyyyy, g_xx_0_y_xyyyyy, g_xx_0_y_xyyyyz, g_xx_0_y_xyyyz, g_xx_0_y_xyyyzz, g_xx_0_y_xyyzz, g_xx_0_y_xyyzzz, g_xx_0_y_xyzzz, g_xx_0_y_xyzzzz, g_xx_0_y_xzzzz, g_xx_0_y_yyyyy, g_xx_0_y_yyyyyy, g_xx_0_y_yyyyyz, g_xx_0_y_yyyyz, g_xx_0_y_yyyyzz, g_xx_0_y_yyyzz, g_xx_0_y_yyyzzz, g_xx_0_y_yyzzz, g_xx_0_y_yyzzzz, g_xx_0_y_yzzzz, g_xx_0_y_yzzzzz, g_xx_0_y_zzzzz, g_xx_0_yy_xxxxx, g_xx_0_yy_xxxxy, g_xx_0_yy_xxxxz, g_xx_0_yy_xxxyy, g_xx_0_yy_xxxyz, g_xx_0_yy_xxxzz, g_xx_0_yy_xxyyy, g_xx_0_yy_xxyyz, g_xx_0_yy_xxyzz, g_xx_0_yy_xxzzz, g_xx_0_yy_xyyyy, g_xx_0_yy_xyyyz, g_xx_0_yy_xyyzz, g_xx_0_yy_xyzzz, g_xx_0_yy_xzzzz, g_xx_0_yy_yyyyy, g_xx_0_yy_yyyyz, g_xx_0_yy_yyyzz, g_xx_0_yy_yyzzz, g_xx_0_yy_yzzzz, g_xx_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yy_xxxxx[k] = -g_xx_0_y_xxxxx[k] * ab_y + g_xx_0_y_xxxxxy[k];

                g_xx_0_yy_xxxxy[k] = -g_xx_0_y_xxxxy[k] * ab_y + g_xx_0_y_xxxxyy[k];

                g_xx_0_yy_xxxxz[k] = -g_xx_0_y_xxxxz[k] * ab_y + g_xx_0_y_xxxxyz[k];

                g_xx_0_yy_xxxyy[k] = -g_xx_0_y_xxxyy[k] * ab_y + g_xx_0_y_xxxyyy[k];

                g_xx_0_yy_xxxyz[k] = -g_xx_0_y_xxxyz[k] * ab_y + g_xx_0_y_xxxyyz[k];

                g_xx_0_yy_xxxzz[k] = -g_xx_0_y_xxxzz[k] * ab_y + g_xx_0_y_xxxyzz[k];

                g_xx_0_yy_xxyyy[k] = -g_xx_0_y_xxyyy[k] * ab_y + g_xx_0_y_xxyyyy[k];

                g_xx_0_yy_xxyyz[k] = -g_xx_0_y_xxyyz[k] * ab_y + g_xx_0_y_xxyyyz[k];

                g_xx_0_yy_xxyzz[k] = -g_xx_0_y_xxyzz[k] * ab_y + g_xx_0_y_xxyyzz[k];

                g_xx_0_yy_xxzzz[k] = -g_xx_0_y_xxzzz[k] * ab_y + g_xx_0_y_xxyzzz[k];

                g_xx_0_yy_xyyyy[k] = -g_xx_0_y_xyyyy[k] * ab_y + g_xx_0_y_xyyyyy[k];

                g_xx_0_yy_xyyyz[k] = -g_xx_0_y_xyyyz[k] * ab_y + g_xx_0_y_xyyyyz[k];

                g_xx_0_yy_xyyzz[k] = -g_xx_0_y_xyyzz[k] * ab_y + g_xx_0_y_xyyyzz[k];

                g_xx_0_yy_xyzzz[k] = -g_xx_0_y_xyzzz[k] * ab_y + g_xx_0_y_xyyzzz[k];

                g_xx_0_yy_xzzzz[k] = -g_xx_0_y_xzzzz[k] * ab_y + g_xx_0_y_xyzzzz[k];

                g_xx_0_yy_yyyyy[k] = -g_xx_0_y_yyyyy[k] * ab_y + g_xx_0_y_yyyyyy[k];

                g_xx_0_yy_yyyyz[k] = -g_xx_0_y_yyyyz[k] * ab_y + g_xx_0_y_yyyyyz[k];

                g_xx_0_yy_yyyzz[k] = -g_xx_0_y_yyyzz[k] * ab_y + g_xx_0_y_yyyyzz[k];

                g_xx_0_yy_yyzzz[k] = -g_xx_0_y_yyzzz[k] * ab_y + g_xx_0_y_yyyzzz[k];

                g_xx_0_yy_yzzzz[k] = -g_xx_0_y_yzzzz[k] * ab_y + g_xx_0_y_yyzzzz[k];

                g_xx_0_yy_zzzzz[k] = -g_xx_0_y_zzzzz[k] * ab_y + g_xx_0_y_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yz_xxxxx, g_xx_0_yz_xxxxy, g_xx_0_yz_xxxxz, g_xx_0_yz_xxxyy, g_xx_0_yz_xxxyz, g_xx_0_yz_xxxzz, g_xx_0_yz_xxyyy, g_xx_0_yz_xxyyz, g_xx_0_yz_xxyzz, g_xx_0_yz_xxzzz, g_xx_0_yz_xyyyy, g_xx_0_yz_xyyyz, g_xx_0_yz_xyyzz, g_xx_0_yz_xyzzz, g_xx_0_yz_xzzzz, g_xx_0_yz_yyyyy, g_xx_0_yz_yyyyz, g_xx_0_yz_yyyzz, g_xx_0_yz_yyzzz, g_xx_0_yz_yzzzz, g_xx_0_yz_zzzzz, g_xx_0_z_xxxxx, g_xx_0_z_xxxxxy, g_xx_0_z_xxxxy, g_xx_0_z_xxxxyy, g_xx_0_z_xxxxyz, g_xx_0_z_xxxxz, g_xx_0_z_xxxyy, g_xx_0_z_xxxyyy, g_xx_0_z_xxxyyz, g_xx_0_z_xxxyz, g_xx_0_z_xxxyzz, g_xx_0_z_xxxzz, g_xx_0_z_xxyyy, g_xx_0_z_xxyyyy, g_xx_0_z_xxyyyz, g_xx_0_z_xxyyz, g_xx_0_z_xxyyzz, g_xx_0_z_xxyzz, g_xx_0_z_xxyzzz, g_xx_0_z_xxzzz, g_xx_0_z_xyyyy, g_xx_0_z_xyyyyy, g_xx_0_z_xyyyyz, g_xx_0_z_xyyyz, g_xx_0_z_xyyyzz, g_xx_0_z_xyyzz, g_xx_0_z_xyyzzz, g_xx_0_z_xyzzz, g_xx_0_z_xyzzzz, g_xx_0_z_xzzzz, g_xx_0_z_yyyyy, g_xx_0_z_yyyyyy, g_xx_0_z_yyyyyz, g_xx_0_z_yyyyz, g_xx_0_z_yyyyzz, g_xx_0_z_yyyzz, g_xx_0_z_yyyzzz, g_xx_0_z_yyzzz, g_xx_0_z_yyzzzz, g_xx_0_z_yzzzz, g_xx_0_z_yzzzzz, g_xx_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yz_xxxxx[k] = -g_xx_0_z_xxxxx[k] * ab_y + g_xx_0_z_xxxxxy[k];

                g_xx_0_yz_xxxxy[k] = -g_xx_0_z_xxxxy[k] * ab_y + g_xx_0_z_xxxxyy[k];

                g_xx_0_yz_xxxxz[k] = -g_xx_0_z_xxxxz[k] * ab_y + g_xx_0_z_xxxxyz[k];

                g_xx_0_yz_xxxyy[k] = -g_xx_0_z_xxxyy[k] * ab_y + g_xx_0_z_xxxyyy[k];

                g_xx_0_yz_xxxyz[k] = -g_xx_0_z_xxxyz[k] * ab_y + g_xx_0_z_xxxyyz[k];

                g_xx_0_yz_xxxzz[k] = -g_xx_0_z_xxxzz[k] * ab_y + g_xx_0_z_xxxyzz[k];

                g_xx_0_yz_xxyyy[k] = -g_xx_0_z_xxyyy[k] * ab_y + g_xx_0_z_xxyyyy[k];

                g_xx_0_yz_xxyyz[k] = -g_xx_0_z_xxyyz[k] * ab_y + g_xx_0_z_xxyyyz[k];

                g_xx_0_yz_xxyzz[k] = -g_xx_0_z_xxyzz[k] * ab_y + g_xx_0_z_xxyyzz[k];

                g_xx_0_yz_xxzzz[k] = -g_xx_0_z_xxzzz[k] * ab_y + g_xx_0_z_xxyzzz[k];

                g_xx_0_yz_xyyyy[k] = -g_xx_0_z_xyyyy[k] * ab_y + g_xx_0_z_xyyyyy[k];

                g_xx_0_yz_xyyyz[k] = -g_xx_0_z_xyyyz[k] * ab_y + g_xx_0_z_xyyyyz[k];

                g_xx_0_yz_xyyzz[k] = -g_xx_0_z_xyyzz[k] * ab_y + g_xx_0_z_xyyyzz[k];

                g_xx_0_yz_xyzzz[k] = -g_xx_0_z_xyzzz[k] * ab_y + g_xx_0_z_xyyzzz[k];

                g_xx_0_yz_xzzzz[k] = -g_xx_0_z_xzzzz[k] * ab_y + g_xx_0_z_xyzzzz[k];

                g_xx_0_yz_yyyyy[k] = -g_xx_0_z_yyyyy[k] * ab_y + g_xx_0_z_yyyyyy[k];

                g_xx_0_yz_yyyyz[k] = -g_xx_0_z_yyyyz[k] * ab_y + g_xx_0_z_yyyyyz[k];

                g_xx_0_yz_yyyzz[k] = -g_xx_0_z_yyyzz[k] * ab_y + g_xx_0_z_yyyyzz[k];

                g_xx_0_yz_yyzzz[k] = -g_xx_0_z_yyzzz[k] * ab_y + g_xx_0_z_yyyzzz[k];

                g_xx_0_yz_yzzzz[k] = -g_xx_0_z_yzzzz[k] * ab_y + g_xx_0_z_yyzzzz[k];

                g_xx_0_yz_zzzzz[k] = -g_xx_0_z_zzzzz[k] * ab_y + g_xx_0_z_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 119 * ccomps * dcomps);

            auto g_xx_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_z_xxxxx, g_xx_0_z_xxxxxz, g_xx_0_z_xxxxy, g_xx_0_z_xxxxyz, g_xx_0_z_xxxxz, g_xx_0_z_xxxxzz, g_xx_0_z_xxxyy, g_xx_0_z_xxxyyz, g_xx_0_z_xxxyz, g_xx_0_z_xxxyzz, g_xx_0_z_xxxzz, g_xx_0_z_xxxzzz, g_xx_0_z_xxyyy, g_xx_0_z_xxyyyz, g_xx_0_z_xxyyz, g_xx_0_z_xxyyzz, g_xx_0_z_xxyzz, g_xx_0_z_xxyzzz, g_xx_0_z_xxzzz, g_xx_0_z_xxzzzz, g_xx_0_z_xyyyy, g_xx_0_z_xyyyyz, g_xx_0_z_xyyyz, g_xx_0_z_xyyyzz, g_xx_0_z_xyyzz, g_xx_0_z_xyyzzz, g_xx_0_z_xyzzz, g_xx_0_z_xyzzzz, g_xx_0_z_xzzzz, g_xx_0_z_xzzzzz, g_xx_0_z_yyyyy, g_xx_0_z_yyyyyz, g_xx_0_z_yyyyz, g_xx_0_z_yyyyzz, g_xx_0_z_yyyzz, g_xx_0_z_yyyzzz, g_xx_0_z_yyzzz, g_xx_0_z_yyzzzz, g_xx_0_z_yzzzz, g_xx_0_z_yzzzzz, g_xx_0_z_zzzzz, g_xx_0_z_zzzzzz, g_xx_0_zz_xxxxx, g_xx_0_zz_xxxxy, g_xx_0_zz_xxxxz, g_xx_0_zz_xxxyy, g_xx_0_zz_xxxyz, g_xx_0_zz_xxxzz, g_xx_0_zz_xxyyy, g_xx_0_zz_xxyyz, g_xx_0_zz_xxyzz, g_xx_0_zz_xxzzz, g_xx_0_zz_xyyyy, g_xx_0_zz_xyyyz, g_xx_0_zz_xyyzz, g_xx_0_zz_xyzzz, g_xx_0_zz_xzzzz, g_xx_0_zz_yyyyy, g_xx_0_zz_yyyyz, g_xx_0_zz_yyyzz, g_xx_0_zz_yyzzz, g_xx_0_zz_yzzzz, g_xx_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zz_xxxxx[k] = -g_xx_0_z_xxxxx[k] * ab_z + g_xx_0_z_xxxxxz[k];

                g_xx_0_zz_xxxxy[k] = -g_xx_0_z_xxxxy[k] * ab_z + g_xx_0_z_xxxxyz[k];

                g_xx_0_zz_xxxxz[k] = -g_xx_0_z_xxxxz[k] * ab_z + g_xx_0_z_xxxxzz[k];

                g_xx_0_zz_xxxyy[k] = -g_xx_0_z_xxxyy[k] * ab_z + g_xx_0_z_xxxyyz[k];

                g_xx_0_zz_xxxyz[k] = -g_xx_0_z_xxxyz[k] * ab_z + g_xx_0_z_xxxyzz[k];

                g_xx_0_zz_xxxzz[k] = -g_xx_0_z_xxxzz[k] * ab_z + g_xx_0_z_xxxzzz[k];

                g_xx_0_zz_xxyyy[k] = -g_xx_0_z_xxyyy[k] * ab_z + g_xx_0_z_xxyyyz[k];

                g_xx_0_zz_xxyyz[k] = -g_xx_0_z_xxyyz[k] * ab_z + g_xx_0_z_xxyyzz[k];

                g_xx_0_zz_xxyzz[k] = -g_xx_0_z_xxyzz[k] * ab_z + g_xx_0_z_xxyzzz[k];

                g_xx_0_zz_xxzzz[k] = -g_xx_0_z_xxzzz[k] * ab_z + g_xx_0_z_xxzzzz[k];

                g_xx_0_zz_xyyyy[k] = -g_xx_0_z_xyyyy[k] * ab_z + g_xx_0_z_xyyyyz[k];

                g_xx_0_zz_xyyyz[k] = -g_xx_0_z_xyyyz[k] * ab_z + g_xx_0_z_xyyyzz[k];

                g_xx_0_zz_xyyzz[k] = -g_xx_0_z_xyyzz[k] * ab_z + g_xx_0_z_xyyzzz[k];

                g_xx_0_zz_xyzzz[k] = -g_xx_0_z_xyzzz[k] * ab_z + g_xx_0_z_xyzzzz[k];

                g_xx_0_zz_xzzzz[k] = -g_xx_0_z_xzzzz[k] * ab_z + g_xx_0_z_xzzzzz[k];

                g_xx_0_zz_yyyyy[k] = -g_xx_0_z_yyyyy[k] * ab_z + g_xx_0_z_yyyyyz[k];

                g_xx_0_zz_yyyyz[k] = -g_xx_0_z_yyyyz[k] * ab_z + g_xx_0_z_yyyyzz[k];

                g_xx_0_zz_yyyzz[k] = -g_xx_0_z_yyyzz[k] * ab_z + g_xx_0_z_yyyzzz[k];

                g_xx_0_zz_yyzzz[k] = -g_xx_0_z_yyzzz[k] * ab_z + g_xx_0_z_yyzzzz[k];

                g_xx_0_zz_yzzzz[k] = -g_xx_0_z_yzzzz[k] * ab_z + g_xx_0_z_yzzzzz[k];

                g_xx_0_zz_zzzzz[k] = -g_xx_0_z_zzzzz[k] * ab_z + g_xx_0_z_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_x_xxxxx, g_xy_0_x_xxxxxx, g_xy_0_x_xxxxxy, g_xy_0_x_xxxxxz, g_xy_0_x_xxxxy, g_xy_0_x_xxxxyy, g_xy_0_x_xxxxyz, g_xy_0_x_xxxxz, g_xy_0_x_xxxxzz, g_xy_0_x_xxxyy, g_xy_0_x_xxxyyy, g_xy_0_x_xxxyyz, g_xy_0_x_xxxyz, g_xy_0_x_xxxyzz, g_xy_0_x_xxxzz, g_xy_0_x_xxxzzz, g_xy_0_x_xxyyy, g_xy_0_x_xxyyyy, g_xy_0_x_xxyyyz, g_xy_0_x_xxyyz, g_xy_0_x_xxyyzz, g_xy_0_x_xxyzz, g_xy_0_x_xxyzzz, g_xy_0_x_xxzzz, g_xy_0_x_xxzzzz, g_xy_0_x_xyyyy, g_xy_0_x_xyyyyy, g_xy_0_x_xyyyyz, g_xy_0_x_xyyyz, g_xy_0_x_xyyyzz, g_xy_0_x_xyyzz, g_xy_0_x_xyyzzz, g_xy_0_x_xyzzz, g_xy_0_x_xyzzzz, g_xy_0_x_xzzzz, g_xy_0_x_xzzzzz, g_xy_0_x_yyyyy, g_xy_0_x_yyyyz, g_xy_0_x_yyyzz, g_xy_0_x_yyzzz, g_xy_0_x_yzzzz, g_xy_0_x_zzzzz, g_xy_0_xx_xxxxx, g_xy_0_xx_xxxxy, g_xy_0_xx_xxxxz, g_xy_0_xx_xxxyy, g_xy_0_xx_xxxyz, g_xy_0_xx_xxxzz, g_xy_0_xx_xxyyy, g_xy_0_xx_xxyyz, g_xy_0_xx_xxyzz, g_xy_0_xx_xxzzz, g_xy_0_xx_xyyyy, g_xy_0_xx_xyyyz, g_xy_0_xx_xyyzz, g_xy_0_xx_xyzzz, g_xy_0_xx_xzzzz, g_xy_0_xx_yyyyy, g_xy_0_xx_yyyyz, g_xy_0_xx_yyyzz, g_xy_0_xx_yyzzz, g_xy_0_xx_yzzzz, g_xy_0_xx_zzzzz, g_y_0_x_xxxxx, g_y_0_x_xxxxy, g_y_0_x_xxxxz, g_y_0_x_xxxyy, g_y_0_x_xxxyz, g_y_0_x_xxxzz, g_y_0_x_xxyyy, g_y_0_x_xxyyz, g_y_0_x_xxyzz, g_y_0_x_xxzzz, g_y_0_x_xyyyy, g_y_0_x_xyyyz, g_y_0_x_xyyzz, g_y_0_x_xyzzz, g_y_0_x_xzzzz, g_y_0_x_yyyyy, g_y_0_x_yyyyz, g_y_0_x_yyyzz, g_y_0_x_yyzzz, g_y_0_x_yzzzz, g_y_0_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xx_xxxxx[k] = -g_y_0_x_xxxxx[k] - g_xy_0_x_xxxxx[k] * ab_x + g_xy_0_x_xxxxxx[k];

                g_xy_0_xx_xxxxy[k] = -g_y_0_x_xxxxy[k] - g_xy_0_x_xxxxy[k] * ab_x + g_xy_0_x_xxxxxy[k];

                g_xy_0_xx_xxxxz[k] = -g_y_0_x_xxxxz[k] - g_xy_0_x_xxxxz[k] * ab_x + g_xy_0_x_xxxxxz[k];

                g_xy_0_xx_xxxyy[k] = -g_y_0_x_xxxyy[k] - g_xy_0_x_xxxyy[k] * ab_x + g_xy_0_x_xxxxyy[k];

                g_xy_0_xx_xxxyz[k] = -g_y_0_x_xxxyz[k] - g_xy_0_x_xxxyz[k] * ab_x + g_xy_0_x_xxxxyz[k];

                g_xy_0_xx_xxxzz[k] = -g_y_0_x_xxxzz[k] - g_xy_0_x_xxxzz[k] * ab_x + g_xy_0_x_xxxxzz[k];

                g_xy_0_xx_xxyyy[k] = -g_y_0_x_xxyyy[k] - g_xy_0_x_xxyyy[k] * ab_x + g_xy_0_x_xxxyyy[k];

                g_xy_0_xx_xxyyz[k] = -g_y_0_x_xxyyz[k] - g_xy_0_x_xxyyz[k] * ab_x + g_xy_0_x_xxxyyz[k];

                g_xy_0_xx_xxyzz[k] = -g_y_0_x_xxyzz[k] - g_xy_0_x_xxyzz[k] * ab_x + g_xy_0_x_xxxyzz[k];

                g_xy_0_xx_xxzzz[k] = -g_y_0_x_xxzzz[k] - g_xy_0_x_xxzzz[k] * ab_x + g_xy_0_x_xxxzzz[k];

                g_xy_0_xx_xyyyy[k] = -g_y_0_x_xyyyy[k] - g_xy_0_x_xyyyy[k] * ab_x + g_xy_0_x_xxyyyy[k];

                g_xy_0_xx_xyyyz[k] = -g_y_0_x_xyyyz[k] - g_xy_0_x_xyyyz[k] * ab_x + g_xy_0_x_xxyyyz[k];

                g_xy_0_xx_xyyzz[k] = -g_y_0_x_xyyzz[k] - g_xy_0_x_xyyzz[k] * ab_x + g_xy_0_x_xxyyzz[k];

                g_xy_0_xx_xyzzz[k] = -g_y_0_x_xyzzz[k] - g_xy_0_x_xyzzz[k] * ab_x + g_xy_0_x_xxyzzz[k];

                g_xy_0_xx_xzzzz[k] = -g_y_0_x_xzzzz[k] - g_xy_0_x_xzzzz[k] * ab_x + g_xy_0_x_xxzzzz[k];

                g_xy_0_xx_yyyyy[k] = -g_y_0_x_yyyyy[k] - g_xy_0_x_yyyyy[k] * ab_x + g_xy_0_x_xyyyyy[k];

                g_xy_0_xx_yyyyz[k] = -g_y_0_x_yyyyz[k] - g_xy_0_x_yyyyz[k] * ab_x + g_xy_0_x_xyyyyz[k];

                g_xy_0_xx_yyyzz[k] = -g_y_0_x_yyyzz[k] - g_xy_0_x_yyyzz[k] * ab_x + g_xy_0_x_xyyyzz[k];

                g_xy_0_xx_yyzzz[k] = -g_y_0_x_yyzzz[k] - g_xy_0_x_yyzzz[k] * ab_x + g_xy_0_x_xyyzzz[k];

                g_xy_0_xx_yzzzz[k] = -g_y_0_x_yzzzz[k] - g_xy_0_x_yzzzz[k] * ab_x + g_xy_0_x_xyzzzz[k];

                g_xy_0_xx_zzzzz[k] = -g_y_0_x_zzzzz[k] - g_xy_0_x_zzzzz[k] * ab_x + g_xy_0_x_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xy_xxxxx, g_xy_0_xy_xxxxy, g_xy_0_xy_xxxxz, g_xy_0_xy_xxxyy, g_xy_0_xy_xxxyz, g_xy_0_xy_xxxzz, g_xy_0_xy_xxyyy, g_xy_0_xy_xxyyz, g_xy_0_xy_xxyzz, g_xy_0_xy_xxzzz, g_xy_0_xy_xyyyy, g_xy_0_xy_xyyyz, g_xy_0_xy_xyyzz, g_xy_0_xy_xyzzz, g_xy_0_xy_xzzzz, g_xy_0_xy_yyyyy, g_xy_0_xy_yyyyz, g_xy_0_xy_yyyzz, g_xy_0_xy_yyzzz, g_xy_0_xy_yzzzz, g_xy_0_xy_zzzzz, g_xy_0_y_xxxxx, g_xy_0_y_xxxxxx, g_xy_0_y_xxxxxy, g_xy_0_y_xxxxxz, g_xy_0_y_xxxxy, g_xy_0_y_xxxxyy, g_xy_0_y_xxxxyz, g_xy_0_y_xxxxz, g_xy_0_y_xxxxzz, g_xy_0_y_xxxyy, g_xy_0_y_xxxyyy, g_xy_0_y_xxxyyz, g_xy_0_y_xxxyz, g_xy_0_y_xxxyzz, g_xy_0_y_xxxzz, g_xy_0_y_xxxzzz, g_xy_0_y_xxyyy, g_xy_0_y_xxyyyy, g_xy_0_y_xxyyyz, g_xy_0_y_xxyyz, g_xy_0_y_xxyyzz, g_xy_0_y_xxyzz, g_xy_0_y_xxyzzz, g_xy_0_y_xxzzz, g_xy_0_y_xxzzzz, g_xy_0_y_xyyyy, g_xy_0_y_xyyyyy, g_xy_0_y_xyyyyz, g_xy_0_y_xyyyz, g_xy_0_y_xyyyzz, g_xy_0_y_xyyzz, g_xy_0_y_xyyzzz, g_xy_0_y_xyzzz, g_xy_0_y_xyzzzz, g_xy_0_y_xzzzz, g_xy_0_y_xzzzzz, g_xy_0_y_yyyyy, g_xy_0_y_yyyyz, g_xy_0_y_yyyzz, g_xy_0_y_yyzzz, g_xy_0_y_yzzzz, g_xy_0_y_zzzzz, g_y_0_y_xxxxx, g_y_0_y_xxxxy, g_y_0_y_xxxxz, g_y_0_y_xxxyy, g_y_0_y_xxxyz, g_y_0_y_xxxzz, g_y_0_y_xxyyy, g_y_0_y_xxyyz, g_y_0_y_xxyzz, g_y_0_y_xxzzz, g_y_0_y_xyyyy, g_y_0_y_xyyyz, g_y_0_y_xyyzz, g_y_0_y_xyzzz, g_y_0_y_xzzzz, g_y_0_y_yyyyy, g_y_0_y_yyyyz, g_y_0_y_yyyzz, g_y_0_y_yyzzz, g_y_0_y_yzzzz, g_y_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xy_xxxxx[k] = -g_y_0_y_xxxxx[k] - g_xy_0_y_xxxxx[k] * ab_x + g_xy_0_y_xxxxxx[k];

                g_xy_0_xy_xxxxy[k] = -g_y_0_y_xxxxy[k] - g_xy_0_y_xxxxy[k] * ab_x + g_xy_0_y_xxxxxy[k];

                g_xy_0_xy_xxxxz[k] = -g_y_0_y_xxxxz[k] - g_xy_0_y_xxxxz[k] * ab_x + g_xy_0_y_xxxxxz[k];

                g_xy_0_xy_xxxyy[k] = -g_y_0_y_xxxyy[k] - g_xy_0_y_xxxyy[k] * ab_x + g_xy_0_y_xxxxyy[k];

                g_xy_0_xy_xxxyz[k] = -g_y_0_y_xxxyz[k] - g_xy_0_y_xxxyz[k] * ab_x + g_xy_0_y_xxxxyz[k];

                g_xy_0_xy_xxxzz[k] = -g_y_0_y_xxxzz[k] - g_xy_0_y_xxxzz[k] * ab_x + g_xy_0_y_xxxxzz[k];

                g_xy_0_xy_xxyyy[k] = -g_y_0_y_xxyyy[k] - g_xy_0_y_xxyyy[k] * ab_x + g_xy_0_y_xxxyyy[k];

                g_xy_0_xy_xxyyz[k] = -g_y_0_y_xxyyz[k] - g_xy_0_y_xxyyz[k] * ab_x + g_xy_0_y_xxxyyz[k];

                g_xy_0_xy_xxyzz[k] = -g_y_0_y_xxyzz[k] - g_xy_0_y_xxyzz[k] * ab_x + g_xy_0_y_xxxyzz[k];

                g_xy_0_xy_xxzzz[k] = -g_y_0_y_xxzzz[k] - g_xy_0_y_xxzzz[k] * ab_x + g_xy_0_y_xxxzzz[k];

                g_xy_0_xy_xyyyy[k] = -g_y_0_y_xyyyy[k] - g_xy_0_y_xyyyy[k] * ab_x + g_xy_0_y_xxyyyy[k];

                g_xy_0_xy_xyyyz[k] = -g_y_0_y_xyyyz[k] - g_xy_0_y_xyyyz[k] * ab_x + g_xy_0_y_xxyyyz[k];

                g_xy_0_xy_xyyzz[k] = -g_y_0_y_xyyzz[k] - g_xy_0_y_xyyzz[k] * ab_x + g_xy_0_y_xxyyzz[k];

                g_xy_0_xy_xyzzz[k] = -g_y_0_y_xyzzz[k] - g_xy_0_y_xyzzz[k] * ab_x + g_xy_0_y_xxyzzz[k];

                g_xy_0_xy_xzzzz[k] = -g_y_0_y_xzzzz[k] - g_xy_0_y_xzzzz[k] * ab_x + g_xy_0_y_xxzzzz[k];

                g_xy_0_xy_yyyyy[k] = -g_y_0_y_yyyyy[k] - g_xy_0_y_yyyyy[k] * ab_x + g_xy_0_y_xyyyyy[k];

                g_xy_0_xy_yyyyz[k] = -g_y_0_y_yyyyz[k] - g_xy_0_y_yyyyz[k] * ab_x + g_xy_0_y_xyyyyz[k];

                g_xy_0_xy_yyyzz[k] = -g_y_0_y_yyyzz[k] - g_xy_0_y_yyyzz[k] * ab_x + g_xy_0_y_xyyyzz[k];

                g_xy_0_xy_yyzzz[k] = -g_y_0_y_yyzzz[k] - g_xy_0_y_yyzzz[k] * ab_x + g_xy_0_y_xyyzzz[k];

                g_xy_0_xy_yzzzz[k] = -g_y_0_y_yzzzz[k] - g_xy_0_y_yzzzz[k] * ab_x + g_xy_0_y_xyzzzz[k];

                g_xy_0_xy_zzzzz[k] = -g_y_0_y_zzzzz[k] - g_xy_0_y_zzzzz[k] * ab_x + g_xy_0_y_xzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 179 * ccomps * dcomps);

            auto g_xy_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_x_xxxxx, g_xy_0_x_xxxxxz, g_xy_0_x_xxxxy, g_xy_0_x_xxxxyz, g_xy_0_x_xxxxz, g_xy_0_x_xxxxzz, g_xy_0_x_xxxyy, g_xy_0_x_xxxyyz, g_xy_0_x_xxxyz, g_xy_0_x_xxxyzz, g_xy_0_x_xxxzz, g_xy_0_x_xxxzzz, g_xy_0_x_xxyyy, g_xy_0_x_xxyyyz, g_xy_0_x_xxyyz, g_xy_0_x_xxyyzz, g_xy_0_x_xxyzz, g_xy_0_x_xxyzzz, g_xy_0_x_xxzzz, g_xy_0_x_xxzzzz, g_xy_0_x_xyyyy, g_xy_0_x_xyyyyz, g_xy_0_x_xyyyz, g_xy_0_x_xyyyzz, g_xy_0_x_xyyzz, g_xy_0_x_xyyzzz, g_xy_0_x_xyzzz, g_xy_0_x_xyzzzz, g_xy_0_x_xzzzz, g_xy_0_x_xzzzzz, g_xy_0_x_yyyyy, g_xy_0_x_yyyyyz, g_xy_0_x_yyyyz, g_xy_0_x_yyyyzz, g_xy_0_x_yyyzz, g_xy_0_x_yyyzzz, g_xy_0_x_yyzzz, g_xy_0_x_yyzzzz, g_xy_0_x_yzzzz, g_xy_0_x_yzzzzz, g_xy_0_x_zzzzz, g_xy_0_x_zzzzzz, g_xy_0_xz_xxxxx, g_xy_0_xz_xxxxy, g_xy_0_xz_xxxxz, g_xy_0_xz_xxxyy, g_xy_0_xz_xxxyz, g_xy_0_xz_xxxzz, g_xy_0_xz_xxyyy, g_xy_0_xz_xxyyz, g_xy_0_xz_xxyzz, g_xy_0_xz_xxzzz, g_xy_0_xz_xyyyy, g_xy_0_xz_xyyyz, g_xy_0_xz_xyyzz, g_xy_0_xz_xyzzz, g_xy_0_xz_xzzzz, g_xy_0_xz_yyyyy, g_xy_0_xz_yyyyz, g_xy_0_xz_yyyzz, g_xy_0_xz_yyzzz, g_xy_0_xz_yzzzz, g_xy_0_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xz_xxxxx[k] = -g_xy_0_x_xxxxx[k] * ab_z + g_xy_0_x_xxxxxz[k];

                g_xy_0_xz_xxxxy[k] = -g_xy_0_x_xxxxy[k] * ab_z + g_xy_0_x_xxxxyz[k];

                g_xy_0_xz_xxxxz[k] = -g_xy_0_x_xxxxz[k] * ab_z + g_xy_0_x_xxxxzz[k];

                g_xy_0_xz_xxxyy[k] = -g_xy_0_x_xxxyy[k] * ab_z + g_xy_0_x_xxxyyz[k];

                g_xy_0_xz_xxxyz[k] = -g_xy_0_x_xxxyz[k] * ab_z + g_xy_0_x_xxxyzz[k];

                g_xy_0_xz_xxxzz[k] = -g_xy_0_x_xxxzz[k] * ab_z + g_xy_0_x_xxxzzz[k];

                g_xy_0_xz_xxyyy[k] = -g_xy_0_x_xxyyy[k] * ab_z + g_xy_0_x_xxyyyz[k];

                g_xy_0_xz_xxyyz[k] = -g_xy_0_x_xxyyz[k] * ab_z + g_xy_0_x_xxyyzz[k];

                g_xy_0_xz_xxyzz[k] = -g_xy_0_x_xxyzz[k] * ab_z + g_xy_0_x_xxyzzz[k];

                g_xy_0_xz_xxzzz[k] = -g_xy_0_x_xxzzz[k] * ab_z + g_xy_0_x_xxzzzz[k];

                g_xy_0_xz_xyyyy[k] = -g_xy_0_x_xyyyy[k] * ab_z + g_xy_0_x_xyyyyz[k];

                g_xy_0_xz_xyyyz[k] = -g_xy_0_x_xyyyz[k] * ab_z + g_xy_0_x_xyyyzz[k];

                g_xy_0_xz_xyyzz[k] = -g_xy_0_x_xyyzz[k] * ab_z + g_xy_0_x_xyyzzz[k];

                g_xy_0_xz_xyzzz[k] = -g_xy_0_x_xyzzz[k] * ab_z + g_xy_0_x_xyzzzz[k];

                g_xy_0_xz_xzzzz[k] = -g_xy_0_x_xzzzz[k] * ab_z + g_xy_0_x_xzzzzz[k];

                g_xy_0_xz_yyyyy[k] = -g_xy_0_x_yyyyy[k] * ab_z + g_xy_0_x_yyyyyz[k];

                g_xy_0_xz_yyyyz[k] = -g_xy_0_x_yyyyz[k] * ab_z + g_xy_0_x_yyyyzz[k];

                g_xy_0_xz_yyyzz[k] = -g_xy_0_x_yyyzz[k] * ab_z + g_xy_0_x_yyyzzz[k];

                g_xy_0_xz_yyzzz[k] = -g_xy_0_x_yyzzz[k] * ab_z + g_xy_0_x_yyzzzz[k];

                g_xy_0_xz_yzzzz[k] = -g_xy_0_x_yzzzz[k] * ab_z + g_xy_0_x_yzzzzz[k];

                g_xy_0_xz_zzzzz[k] = -g_xy_0_x_zzzzz[k] * ab_z + g_xy_0_x_zzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 203 * ccomps * dcomps);

            auto g_xy_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxxxx, g_x_0_y_xxxxy, g_x_0_y_xxxxz, g_x_0_y_xxxyy, g_x_0_y_xxxyz, g_x_0_y_xxxzz, g_x_0_y_xxyyy, g_x_0_y_xxyyz, g_x_0_y_xxyzz, g_x_0_y_xxzzz, g_x_0_y_xyyyy, g_x_0_y_xyyyz, g_x_0_y_xyyzz, g_x_0_y_xyzzz, g_x_0_y_xzzzz, g_x_0_y_yyyyy, g_x_0_y_yyyyz, g_x_0_y_yyyzz, g_x_0_y_yyzzz, g_x_0_y_yzzzz, g_x_0_y_zzzzz, g_xy_0_y_xxxxx, g_xy_0_y_xxxxxy, g_xy_0_y_xxxxy, g_xy_0_y_xxxxyy, g_xy_0_y_xxxxyz, g_xy_0_y_xxxxz, g_xy_0_y_xxxyy, g_xy_0_y_xxxyyy, g_xy_0_y_xxxyyz, g_xy_0_y_xxxyz, g_xy_0_y_xxxyzz, g_xy_0_y_xxxzz, g_xy_0_y_xxyyy, g_xy_0_y_xxyyyy, g_xy_0_y_xxyyyz, g_xy_0_y_xxyyz, g_xy_0_y_xxyyzz, g_xy_0_y_xxyzz, g_xy_0_y_xxyzzz, g_xy_0_y_xxzzz, g_xy_0_y_xyyyy, g_xy_0_y_xyyyyy, g_xy_0_y_xyyyyz, g_xy_0_y_xyyyz, g_xy_0_y_xyyyzz, g_xy_0_y_xyyzz, g_xy_0_y_xyyzzz, g_xy_0_y_xyzzz, g_xy_0_y_xyzzzz, g_xy_0_y_xzzzz, g_xy_0_y_yyyyy, g_xy_0_y_yyyyyy, g_xy_0_y_yyyyyz, g_xy_0_y_yyyyz, g_xy_0_y_yyyyzz, g_xy_0_y_yyyzz, g_xy_0_y_yyyzzz, g_xy_0_y_yyzzz, g_xy_0_y_yyzzzz, g_xy_0_y_yzzzz, g_xy_0_y_yzzzzz, g_xy_0_y_zzzzz, g_xy_0_yy_xxxxx, g_xy_0_yy_xxxxy, g_xy_0_yy_xxxxz, g_xy_0_yy_xxxyy, g_xy_0_yy_xxxyz, g_xy_0_yy_xxxzz, g_xy_0_yy_xxyyy, g_xy_0_yy_xxyyz, g_xy_0_yy_xxyzz, g_xy_0_yy_xxzzz, g_xy_0_yy_xyyyy, g_xy_0_yy_xyyyz, g_xy_0_yy_xyyzz, g_xy_0_yy_xyzzz, g_xy_0_yy_xzzzz, g_xy_0_yy_yyyyy, g_xy_0_yy_yyyyz, g_xy_0_yy_yyyzz, g_xy_0_yy_yyzzz, g_xy_0_yy_yzzzz, g_xy_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yy_xxxxx[k] = -g_x_0_y_xxxxx[k] - g_xy_0_y_xxxxx[k] * ab_y + g_xy_0_y_xxxxxy[k];

                g_xy_0_yy_xxxxy[k] = -g_x_0_y_xxxxy[k] - g_xy_0_y_xxxxy[k] * ab_y + g_xy_0_y_xxxxyy[k];

                g_xy_0_yy_xxxxz[k] = -g_x_0_y_xxxxz[k] - g_xy_0_y_xxxxz[k] * ab_y + g_xy_0_y_xxxxyz[k];

                g_xy_0_yy_xxxyy[k] = -g_x_0_y_xxxyy[k] - g_xy_0_y_xxxyy[k] * ab_y + g_xy_0_y_xxxyyy[k];

                g_xy_0_yy_xxxyz[k] = -g_x_0_y_xxxyz[k] - g_xy_0_y_xxxyz[k] * ab_y + g_xy_0_y_xxxyyz[k];

                g_xy_0_yy_xxxzz[k] = -g_x_0_y_xxxzz[k] - g_xy_0_y_xxxzz[k] * ab_y + g_xy_0_y_xxxyzz[k];

                g_xy_0_yy_xxyyy[k] = -g_x_0_y_xxyyy[k] - g_xy_0_y_xxyyy[k] * ab_y + g_xy_0_y_xxyyyy[k];

                g_xy_0_yy_xxyyz[k] = -g_x_0_y_xxyyz[k] - g_xy_0_y_xxyyz[k] * ab_y + g_xy_0_y_xxyyyz[k];

                g_xy_0_yy_xxyzz[k] = -g_x_0_y_xxyzz[k] - g_xy_0_y_xxyzz[k] * ab_y + g_xy_0_y_xxyyzz[k];

                g_xy_0_yy_xxzzz[k] = -g_x_0_y_xxzzz[k] - g_xy_0_y_xxzzz[k] * ab_y + g_xy_0_y_xxyzzz[k];

                g_xy_0_yy_xyyyy[k] = -g_x_0_y_xyyyy[k] - g_xy_0_y_xyyyy[k] * ab_y + g_xy_0_y_xyyyyy[k];

                g_xy_0_yy_xyyyz[k] = -g_x_0_y_xyyyz[k] - g_xy_0_y_xyyyz[k] * ab_y + g_xy_0_y_xyyyyz[k];

                g_xy_0_yy_xyyzz[k] = -g_x_0_y_xyyzz[k] - g_xy_0_y_xyyzz[k] * ab_y + g_xy_0_y_xyyyzz[k];

                g_xy_0_yy_xyzzz[k] = -g_x_0_y_xyzzz[k] - g_xy_0_y_xyzzz[k] * ab_y + g_xy_0_y_xyyzzz[k];

                g_xy_0_yy_xzzzz[k] = -g_x_0_y_xzzzz[k] - g_xy_0_y_xzzzz[k] * ab_y + g_xy_0_y_xyzzzz[k];

                g_xy_0_yy_yyyyy[k] = -g_x_0_y_yyyyy[k] - g_xy_0_y_yyyyy[k] * ab_y + g_xy_0_y_yyyyyy[k];

                g_xy_0_yy_yyyyz[k] = -g_x_0_y_yyyyz[k] - g_xy_0_y_yyyyz[k] * ab_y + g_xy_0_y_yyyyyz[k];

                g_xy_0_yy_yyyzz[k] = -g_x_0_y_yyyzz[k] - g_xy_0_y_yyyzz[k] * ab_y + g_xy_0_y_yyyyzz[k];

                g_xy_0_yy_yyzzz[k] = -g_x_0_y_yyzzz[k] - g_xy_0_y_yyzzz[k] * ab_y + g_xy_0_y_yyyzzz[k];

                g_xy_0_yy_yzzzz[k] = -g_x_0_y_yzzzz[k] - g_xy_0_y_yzzzz[k] * ab_y + g_xy_0_y_yyzzzz[k];

                g_xy_0_yy_zzzzz[k] = -g_x_0_y_zzzzz[k] - g_xy_0_y_zzzzz[k] * ab_y + g_xy_0_y_yzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 215 * ccomps * dcomps);

            auto g_xy_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 221 * ccomps * dcomps);

            auto g_xy_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_y_xxxxx, g_xy_0_y_xxxxxz, g_xy_0_y_xxxxy, g_xy_0_y_xxxxyz, g_xy_0_y_xxxxz, g_xy_0_y_xxxxzz, g_xy_0_y_xxxyy, g_xy_0_y_xxxyyz, g_xy_0_y_xxxyz, g_xy_0_y_xxxyzz, g_xy_0_y_xxxzz, g_xy_0_y_xxxzzz, g_xy_0_y_xxyyy, g_xy_0_y_xxyyyz, g_xy_0_y_xxyyz, g_xy_0_y_xxyyzz, g_xy_0_y_xxyzz, g_xy_0_y_xxyzzz, g_xy_0_y_xxzzz, g_xy_0_y_xxzzzz, g_xy_0_y_xyyyy, g_xy_0_y_xyyyyz, g_xy_0_y_xyyyz, g_xy_0_y_xyyyzz, g_xy_0_y_xyyzz, g_xy_0_y_xyyzzz, g_xy_0_y_xyzzz, g_xy_0_y_xyzzzz, g_xy_0_y_xzzzz, g_xy_0_y_xzzzzz, g_xy_0_y_yyyyy, g_xy_0_y_yyyyyz, g_xy_0_y_yyyyz, g_xy_0_y_yyyyzz, g_xy_0_y_yyyzz, g_xy_0_y_yyyzzz, g_xy_0_y_yyzzz, g_xy_0_y_yyzzzz, g_xy_0_y_yzzzz, g_xy_0_y_yzzzzz, g_xy_0_y_zzzzz, g_xy_0_y_zzzzzz, g_xy_0_yz_xxxxx, g_xy_0_yz_xxxxy, g_xy_0_yz_xxxxz, g_xy_0_yz_xxxyy, g_xy_0_yz_xxxyz, g_xy_0_yz_xxxzz, g_xy_0_yz_xxyyy, g_xy_0_yz_xxyyz, g_xy_0_yz_xxyzz, g_xy_0_yz_xxzzz, g_xy_0_yz_xyyyy, g_xy_0_yz_xyyyz, g_xy_0_yz_xyyzz, g_xy_0_yz_xyzzz, g_xy_0_yz_xzzzz, g_xy_0_yz_yyyyy, g_xy_0_yz_yyyyz, g_xy_0_yz_yyyzz, g_xy_0_yz_yyzzz, g_xy_0_yz_yzzzz, g_xy_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yz_xxxxx[k] = -g_xy_0_y_xxxxx[k] * ab_z + g_xy_0_y_xxxxxz[k];

                g_xy_0_yz_xxxxy[k] = -g_xy_0_y_xxxxy[k] * ab_z + g_xy_0_y_xxxxyz[k];

                g_xy_0_yz_xxxxz[k] = -g_xy_0_y_xxxxz[k] * ab_z + g_xy_0_y_xxxxzz[k];

                g_xy_0_yz_xxxyy[k] = -g_xy_0_y_xxxyy[k] * ab_z + g_xy_0_y_xxxyyz[k];

                g_xy_0_yz_xxxyz[k] = -g_xy_0_y_xxxyz[k] * ab_z + g_xy_0_y_xxxyzz[k];

                g_xy_0_yz_xxxzz[k] = -g_xy_0_y_xxxzz[k] * ab_z + g_xy_0_y_xxxzzz[k];

                g_xy_0_yz_xxyyy[k] = -g_xy_0_y_xxyyy[k] * ab_z + g_xy_0_y_xxyyyz[k];

                g_xy_0_yz_xxyyz[k] = -g_xy_0_y_xxyyz[k] * ab_z + g_xy_0_y_xxyyzz[k];

                g_xy_0_yz_xxyzz[k] = -g_xy_0_y_xxyzz[k] * ab_z + g_xy_0_y_xxyzzz[k];

                g_xy_0_yz_xxzzz[k] = -g_xy_0_y_xxzzz[k] * ab_z + g_xy_0_y_xxzzzz[k];

                g_xy_0_yz_xyyyy[k] = -g_xy_0_y_xyyyy[k] * ab_z + g_xy_0_y_xyyyyz[k];

                g_xy_0_yz_xyyyz[k] = -g_xy_0_y_xyyyz[k] * ab_z + g_xy_0_y_xyyyzz[k];

                g_xy_0_yz_xyyzz[k] = -g_xy_0_y_xyyzz[k] * ab_z + g_xy_0_y_xyyzzz[k];

                g_xy_0_yz_xyzzz[k] = -g_xy_0_y_xyzzz[k] * ab_z + g_xy_0_y_xyzzzz[k];

                g_xy_0_yz_xzzzz[k] = -g_xy_0_y_xzzzz[k] * ab_z + g_xy_0_y_xzzzzz[k];

                g_xy_0_yz_yyyyy[k] = -g_xy_0_y_yyyyy[k] * ab_z + g_xy_0_y_yyyyyz[k];

                g_xy_0_yz_yyyyz[k] = -g_xy_0_y_yyyyz[k] * ab_z + g_xy_0_y_yyyyzz[k];

                g_xy_0_yz_yyyzz[k] = -g_xy_0_y_yyyzz[k] * ab_z + g_xy_0_y_yyyzzz[k];

                g_xy_0_yz_yyzzz[k] = -g_xy_0_y_yyzzz[k] * ab_z + g_xy_0_y_yyzzzz[k];

                g_xy_0_yz_yzzzz[k] = -g_xy_0_y_yzzzz[k] * ab_z + g_xy_0_y_yzzzzz[k];

                g_xy_0_yz_zzzzz[k] = -g_xy_0_y_zzzzz[k] * ab_z + g_xy_0_y_zzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 239 * ccomps * dcomps);

            auto g_xy_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_z_xxxxx, g_xy_0_z_xxxxxz, g_xy_0_z_xxxxy, g_xy_0_z_xxxxyz, g_xy_0_z_xxxxz, g_xy_0_z_xxxxzz, g_xy_0_z_xxxyy, g_xy_0_z_xxxyyz, g_xy_0_z_xxxyz, g_xy_0_z_xxxyzz, g_xy_0_z_xxxzz, g_xy_0_z_xxxzzz, g_xy_0_z_xxyyy, g_xy_0_z_xxyyyz, g_xy_0_z_xxyyz, g_xy_0_z_xxyyzz, g_xy_0_z_xxyzz, g_xy_0_z_xxyzzz, g_xy_0_z_xxzzz, g_xy_0_z_xxzzzz, g_xy_0_z_xyyyy, g_xy_0_z_xyyyyz, g_xy_0_z_xyyyz, g_xy_0_z_xyyyzz, g_xy_0_z_xyyzz, g_xy_0_z_xyyzzz, g_xy_0_z_xyzzz, g_xy_0_z_xyzzzz, g_xy_0_z_xzzzz, g_xy_0_z_xzzzzz, g_xy_0_z_yyyyy, g_xy_0_z_yyyyyz, g_xy_0_z_yyyyz, g_xy_0_z_yyyyzz, g_xy_0_z_yyyzz, g_xy_0_z_yyyzzz, g_xy_0_z_yyzzz, g_xy_0_z_yyzzzz, g_xy_0_z_yzzzz, g_xy_0_z_yzzzzz, g_xy_0_z_zzzzz, g_xy_0_z_zzzzzz, g_xy_0_zz_xxxxx, g_xy_0_zz_xxxxy, g_xy_0_zz_xxxxz, g_xy_0_zz_xxxyy, g_xy_0_zz_xxxyz, g_xy_0_zz_xxxzz, g_xy_0_zz_xxyyy, g_xy_0_zz_xxyyz, g_xy_0_zz_xxyzz, g_xy_0_zz_xxzzz, g_xy_0_zz_xyyyy, g_xy_0_zz_xyyyz, g_xy_0_zz_xyyzz, g_xy_0_zz_xyzzz, g_xy_0_zz_xzzzz, g_xy_0_zz_yyyyy, g_xy_0_zz_yyyyz, g_xy_0_zz_yyyzz, g_xy_0_zz_yyzzz, g_xy_0_zz_yzzzz, g_xy_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zz_xxxxx[k] = -g_xy_0_z_xxxxx[k] * ab_z + g_xy_0_z_xxxxxz[k];

                g_xy_0_zz_xxxxy[k] = -g_xy_0_z_xxxxy[k] * ab_z + g_xy_0_z_xxxxyz[k];

                g_xy_0_zz_xxxxz[k] = -g_xy_0_z_xxxxz[k] * ab_z + g_xy_0_z_xxxxzz[k];

                g_xy_0_zz_xxxyy[k] = -g_xy_0_z_xxxyy[k] * ab_z + g_xy_0_z_xxxyyz[k];

                g_xy_0_zz_xxxyz[k] = -g_xy_0_z_xxxyz[k] * ab_z + g_xy_0_z_xxxyzz[k];

                g_xy_0_zz_xxxzz[k] = -g_xy_0_z_xxxzz[k] * ab_z + g_xy_0_z_xxxzzz[k];

                g_xy_0_zz_xxyyy[k] = -g_xy_0_z_xxyyy[k] * ab_z + g_xy_0_z_xxyyyz[k];

                g_xy_0_zz_xxyyz[k] = -g_xy_0_z_xxyyz[k] * ab_z + g_xy_0_z_xxyyzz[k];

                g_xy_0_zz_xxyzz[k] = -g_xy_0_z_xxyzz[k] * ab_z + g_xy_0_z_xxyzzz[k];

                g_xy_0_zz_xxzzz[k] = -g_xy_0_z_xxzzz[k] * ab_z + g_xy_0_z_xxzzzz[k];

                g_xy_0_zz_xyyyy[k] = -g_xy_0_z_xyyyy[k] * ab_z + g_xy_0_z_xyyyyz[k];

                g_xy_0_zz_xyyyz[k] = -g_xy_0_z_xyyyz[k] * ab_z + g_xy_0_z_xyyyzz[k];

                g_xy_0_zz_xyyzz[k] = -g_xy_0_z_xyyzz[k] * ab_z + g_xy_0_z_xyyzzz[k];

                g_xy_0_zz_xyzzz[k] = -g_xy_0_z_xyzzz[k] * ab_z + g_xy_0_z_xyzzzz[k];

                g_xy_0_zz_xzzzz[k] = -g_xy_0_z_xzzzz[k] * ab_z + g_xy_0_z_xzzzzz[k];

                g_xy_0_zz_yyyyy[k] = -g_xy_0_z_yyyyy[k] * ab_z + g_xy_0_z_yyyyyz[k];

                g_xy_0_zz_yyyyz[k] = -g_xy_0_z_yyyyz[k] * ab_z + g_xy_0_z_yyyyzz[k];

                g_xy_0_zz_yyyzz[k] = -g_xy_0_z_yyyzz[k] * ab_z + g_xy_0_z_yyyzzz[k];

                g_xy_0_zz_yyzzz[k] = -g_xy_0_z_yyzzz[k] * ab_z + g_xy_0_z_yyzzzz[k];

                g_xy_0_zz_yzzzz[k] = -g_xy_0_z_yzzzz[k] * ab_z + g_xy_0_z_yzzzzz[k];

                g_xy_0_zz_zzzzz[k] = -g_xy_0_z_zzzzz[k] * ab_z + g_xy_0_z_zzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 269 * ccomps * dcomps);

            auto g_xz_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 270 * ccomps * dcomps);

            auto g_xz_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 271 * ccomps * dcomps);

            auto g_xz_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_x_xxxxx, g_xz_0_x_xxxxxx, g_xz_0_x_xxxxxy, g_xz_0_x_xxxxxz, g_xz_0_x_xxxxy, g_xz_0_x_xxxxyy, g_xz_0_x_xxxxyz, g_xz_0_x_xxxxz, g_xz_0_x_xxxxzz, g_xz_0_x_xxxyy, g_xz_0_x_xxxyyy, g_xz_0_x_xxxyyz, g_xz_0_x_xxxyz, g_xz_0_x_xxxyzz, g_xz_0_x_xxxzz, g_xz_0_x_xxxzzz, g_xz_0_x_xxyyy, g_xz_0_x_xxyyyy, g_xz_0_x_xxyyyz, g_xz_0_x_xxyyz, g_xz_0_x_xxyyzz, g_xz_0_x_xxyzz, g_xz_0_x_xxyzzz, g_xz_0_x_xxzzz, g_xz_0_x_xxzzzz, g_xz_0_x_xyyyy, g_xz_0_x_xyyyyy, g_xz_0_x_xyyyyz, g_xz_0_x_xyyyz, g_xz_0_x_xyyyzz, g_xz_0_x_xyyzz, g_xz_0_x_xyyzzz, g_xz_0_x_xyzzz, g_xz_0_x_xyzzzz, g_xz_0_x_xzzzz, g_xz_0_x_xzzzzz, g_xz_0_x_yyyyy, g_xz_0_x_yyyyz, g_xz_0_x_yyyzz, g_xz_0_x_yyzzz, g_xz_0_x_yzzzz, g_xz_0_x_zzzzz, g_xz_0_xx_xxxxx, g_xz_0_xx_xxxxy, g_xz_0_xx_xxxxz, g_xz_0_xx_xxxyy, g_xz_0_xx_xxxyz, g_xz_0_xx_xxxzz, g_xz_0_xx_xxyyy, g_xz_0_xx_xxyyz, g_xz_0_xx_xxyzz, g_xz_0_xx_xxzzz, g_xz_0_xx_xyyyy, g_xz_0_xx_xyyyz, g_xz_0_xx_xyyzz, g_xz_0_xx_xyzzz, g_xz_0_xx_xzzzz, g_xz_0_xx_yyyyy, g_xz_0_xx_yyyyz, g_xz_0_xx_yyyzz, g_xz_0_xx_yyzzz, g_xz_0_xx_yzzzz, g_xz_0_xx_zzzzz, g_z_0_x_xxxxx, g_z_0_x_xxxxy, g_z_0_x_xxxxz, g_z_0_x_xxxyy, g_z_0_x_xxxyz, g_z_0_x_xxxzz, g_z_0_x_xxyyy, g_z_0_x_xxyyz, g_z_0_x_xxyzz, g_z_0_x_xxzzz, g_z_0_x_xyyyy, g_z_0_x_xyyyz, g_z_0_x_xyyzz, g_z_0_x_xyzzz, g_z_0_x_xzzzz, g_z_0_x_yyyyy, g_z_0_x_yyyyz, g_z_0_x_yyyzz, g_z_0_x_yyzzz, g_z_0_x_yzzzz, g_z_0_x_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xx_xxxxx[k] = -g_z_0_x_xxxxx[k] - g_xz_0_x_xxxxx[k] * ab_x + g_xz_0_x_xxxxxx[k];

                g_xz_0_xx_xxxxy[k] = -g_z_0_x_xxxxy[k] - g_xz_0_x_xxxxy[k] * ab_x + g_xz_0_x_xxxxxy[k];

                g_xz_0_xx_xxxxz[k] = -g_z_0_x_xxxxz[k] - g_xz_0_x_xxxxz[k] * ab_x + g_xz_0_x_xxxxxz[k];

                g_xz_0_xx_xxxyy[k] = -g_z_0_x_xxxyy[k] - g_xz_0_x_xxxyy[k] * ab_x + g_xz_0_x_xxxxyy[k];

                g_xz_0_xx_xxxyz[k] = -g_z_0_x_xxxyz[k] - g_xz_0_x_xxxyz[k] * ab_x + g_xz_0_x_xxxxyz[k];

                g_xz_0_xx_xxxzz[k] = -g_z_0_x_xxxzz[k] - g_xz_0_x_xxxzz[k] * ab_x + g_xz_0_x_xxxxzz[k];

                g_xz_0_xx_xxyyy[k] = -g_z_0_x_xxyyy[k] - g_xz_0_x_xxyyy[k] * ab_x + g_xz_0_x_xxxyyy[k];

                g_xz_0_xx_xxyyz[k] = -g_z_0_x_xxyyz[k] - g_xz_0_x_xxyyz[k] * ab_x + g_xz_0_x_xxxyyz[k];

                g_xz_0_xx_xxyzz[k] = -g_z_0_x_xxyzz[k] - g_xz_0_x_xxyzz[k] * ab_x + g_xz_0_x_xxxyzz[k];

                g_xz_0_xx_xxzzz[k] = -g_z_0_x_xxzzz[k] - g_xz_0_x_xxzzz[k] * ab_x + g_xz_0_x_xxxzzz[k];

                g_xz_0_xx_xyyyy[k] = -g_z_0_x_xyyyy[k] - g_xz_0_x_xyyyy[k] * ab_x + g_xz_0_x_xxyyyy[k];

                g_xz_0_xx_xyyyz[k] = -g_z_0_x_xyyyz[k] - g_xz_0_x_xyyyz[k] * ab_x + g_xz_0_x_xxyyyz[k];

                g_xz_0_xx_xyyzz[k] = -g_z_0_x_xyyzz[k] - g_xz_0_x_xyyzz[k] * ab_x + g_xz_0_x_xxyyzz[k];

                g_xz_0_xx_xyzzz[k] = -g_z_0_x_xyzzz[k] - g_xz_0_x_xyzzz[k] * ab_x + g_xz_0_x_xxyzzz[k];

                g_xz_0_xx_xzzzz[k] = -g_z_0_x_xzzzz[k] - g_xz_0_x_xzzzz[k] * ab_x + g_xz_0_x_xxzzzz[k];

                g_xz_0_xx_yyyyy[k] = -g_z_0_x_yyyyy[k] - g_xz_0_x_yyyyy[k] * ab_x + g_xz_0_x_xyyyyy[k];

                g_xz_0_xx_yyyyz[k] = -g_z_0_x_yyyyz[k] - g_xz_0_x_yyyyz[k] * ab_x + g_xz_0_x_xyyyyz[k];

                g_xz_0_xx_yyyzz[k] = -g_z_0_x_yyyzz[k] - g_xz_0_x_yyyzz[k] * ab_x + g_xz_0_x_xyyyzz[k];

                g_xz_0_xx_yyzzz[k] = -g_z_0_x_yyzzz[k] - g_xz_0_x_yyzzz[k] * ab_x + g_xz_0_x_xyyzzz[k];

                g_xz_0_xx_yzzzz[k] = -g_z_0_x_yzzzz[k] - g_xz_0_x_yzzzz[k] * ab_x + g_xz_0_x_xyzzzz[k];

                g_xz_0_xx_zzzzz[k] = -g_z_0_x_zzzzz[k] - g_xz_0_x_zzzzz[k] * ab_x + g_xz_0_x_xzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 273 * ccomps * dcomps);

            auto g_xz_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 274 * ccomps * dcomps);

            auto g_xz_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 275 * ccomps * dcomps);

            auto g_xz_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 276 * ccomps * dcomps);

            auto g_xz_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 277 * ccomps * dcomps);

            auto g_xz_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 278 * ccomps * dcomps);

            auto g_xz_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 279 * ccomps * dcomps);

            auto g_xz_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 280 * ccomps * dcomps);

            auto g_xz_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 281 * ccomps * dcomps);

            auto g_xz_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 282 * ccomps * dcomps);

            auto g_xz_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 283 * ccomps * dcomps);

            auto g_xz_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 284 * ccomps * dcomps);

            auto g_xz_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 285 * ccomps * dcomps);

            auto g_xz_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 286 * ccomps * dcomps);

            auto g_xz_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 287 * ccomps * dcomps);

            auto g_xz_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 288 * ccomps * dcomps);

            auto g_xz_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 289 * ccomps * dcomps);

            auto g_xz_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 290 * ccomps * dcomps);

            auto g_xz_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 291 * ccomps * dcomps);

            auto g_xz_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 292 * ccomps * dcomps);

            auto g_xz_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_x_xxxxx, g_xz_0_x_xxxxxy, g_xz_0_x_xxxxy, g_xz_0_x_xxxxyy, g_xz_0_x_xxxxyz, g_xz_0_x_xxxxz, g_xz_0_x_xxxyy, g_xz_0_x_xxxyyy, g_xz_0_x_xxxyyz, g_xz_0_x_xxxyz, g_xz_0_x_xxxyzz, g_xz_0_x_xxxzz, g_xz_0_x_xxyyy, g_xz_0_x_xxyyyy, g_xz_0_x_xxyyyz, g_xz_0_x_xxyyz, g_xz_0_x_xxyyzz, g_xz_0_x_xxyzz, g_xz_0_x_xxyzzz, g_xz_0_x_xxzzz, g_xz_0_x_xyyyy, g_xz_0_x_xyyyyy, g_xz_0_x_xyyyyz, g_xz_0_x_xyyyz, g_xz_0_x_xyyyzz, g_xz_0_x_xyyzz, g_xz_0_x_xyyzzz, g_xz_0_x_xyzzz, g_xz_0_x_xyzzzz, g_xz_0_x_xzzzz, g_xz_0_x_yyyyy, g_xz_0_x_yyyyyy, g_xz_0_x_yyyyyz, g_xz_0_x_yyyyz, g_xz_0_x_yyyyzz, g_xz_0_x_yyyzz, g_xz_0_x_yyyzzz, g_xz_0_x_yyzzz, g_xz_0_x_yyzzzz, g_xz_0_x_yzzzz, g_xz_0_x_yzzzzz, g_xz_0_x_zzzzz, g_xz_0_xy_xxxxx, g_xz_0_xy_xxxxy, g_xz_0_xy_xxxxz, g_xz_0_xy_xxxyy, g_xz_0_xy_xxxyz, g_xz_0_xy_xxxzz, g_xz_0_xy_xxyyy, g_xz_0_xy_xxyyz, g_xz_0_xy_xxyzz, g_xz_0_xy_xxzzz, g_xz_0_xy_xyyyy, g_xz_0_xy_xyyyz, g_xz_0_xy_xyyzz, g_xz_0_xy_xyzzz, g_xz_0_xy_xzzzz, g_xz_0_xy_yyyyy, g_xz_0_xy_yyyyz, g_xz_0_xy_yyyzz, g_xz_0_xy_yyzzz, g_xz_0_xy_yzzzz, g_xz_0_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xy_xxxxx[k] = -g_xz_0_x_xxxxx[k] * ab_y + g_xz_0_x_xxxxxy[k];

                g_xz_0_xy_xxxxy[k] = -g_xz_0_x_xxxxy[k] * ab_y + g_xz_0_x_xxxxyy[k];

                g_xz_0_xy_xxxxz[k] = -g_xz_0_x_xxxxz[k] * ab_y + g_xz_0_x_xxxxyz[k];

                g_xz_0_xy_xxxyy[k] = -g_xz_0_x_xxxyy[k] * ab_y + g_xz_0_x_xxxyyy[k];

                g_xz_0_xy_xxxyz[k] = -g_xz_0_x_xxxyz[k] * ab_y + g_xz_0_x_xxxyyz[k];

                g_xz_0_xy_xxxzz[k] = -g_xz_0_x_xxxzz[k] * ab_y + g_xz_0_x_xxxyzz[k];

                g_xz_0_xy_xxyyy[k] = -g_xz_0_x_xxyyy[k] * ab_y + g_xz_0_x_xxyyyy[k];

                g_xz_0_xy_xxyyz[k] = -g_xz_0_x_xxyyz[k] * ab_y + g_xz_0_x_xxyyyz[k];

                g_xz_0_xy_xxyzz[k] = -g_xz_0_x_xxyzz[k] * ab_y + g_xz_0_x_xxyyzz[k];

                g_xz_0_xy_xxzzz[k] = -g_xz_0_x_xxzzz[k] * ab_y + g_xz_0_x_xxyzzz[k];

                g_xz_0_xy_xyyyy[k] = -g_xz_0_x_xyyyy[k] * ab_y + g_xz_0_x_xyyyyy[k];

                g_xz_0_xy_xyyyz[k] = -g_xz_0_x_xyyyz[k] * ab_y + g_xz_0_x_xyyyyz[k];

                g_xz_0_xy_xyyzz[k] = -g_xz_0_x_xyyzz[k] * ab_y + g_xz_0_x_xyyyzz[k];

                g_xz_0_xy_xyzzz[k] = -g_xz_0_x_xyzzz[k] * ab_y + g_xz_0_x_xyyzzz[k];

                g_xz_0_xy_xzzzz[k] = -g_xz_0_x_xzzzz[k] * ab_y + g_xz_0_x_xyzzzz[k];

                g_xz_0_xy_yyyyy[k] = -g_xz_0_x_yyyyy[k] * ab_y + g_xz_0_x_yyyyyy[k];

                g_xz_0_xy_yyyyz[k] = -g_xz_0_x_yyyyz[k] * ab_y + g_xz_0_x_yyyyyz[k];

                g_xz_0_xy_yyyzz[k] = -g_xz_0_x_yyyzz[k] * ab_y + g_xz_0_x_yyyyzz[k];

                g_xz_0_xy_yyzzz[k] = -g_xz_0_x_yyzzz[k] * ab_y + g_xz_0_x_yyyzzz[k];

                g_xz_0_xy_yzzzz[k] = -g_xz_0_x_yzzzz[k] * ab_y + g_xz_0_x_yyzzzz[k];

                g_xz_0_xy_zzzzz[k] = -g_xz_0_x_zzzzz[k] * ab_y + g_xz_0_x_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 294 * ccomps * dcomps);

            auto g_xz_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 295 * ccomps * dcomps);

            auto g_xz_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 296 * ccomps * dcomps);

            auto g_xz_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 297 * ccomps * dcomps);

            auto g_xz_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 298 * ccomps * dcomps);

            auto g_xz_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 299 * ccomps * dcomps);

            auto g_xz_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 300 * ccomps * dcomps);

            auto g_xz_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 301 * ccomps * dcomps);

            auto g_xz_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 302 * ccomps * dcomps);

            auto g_xz_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 303 * ccomps * dcomps);

            auto g_xz_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 304 * ccomps * dcomps);

            auto g_xz_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 305 * ccomps * dcomps);

            auto g_xz_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 306 * ccomps * dcomps);

            auto g_xz_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 307 * ccomps * dcomps);

            auto g_xz_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 308 * ccomps * dcomps);

            auto g_xz_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 309 * ccomps * dcomps);

            auto g_xz_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 310 * ccomps * dcomps);

            auto g_xz_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 311 * ccomps * dcomps);

            auto g_xz_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 312 * ccomps * dcomps);

            auto g_xz_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 313 * ccomps * dcomps);

            auto g_xz_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xz_xxxxx, g_xz_0_xz_xxxxy, g_xz_0_xz_xxxxz, g_xz_0_xz_xxxyy, g_xz_0_xz_xxxyz, g_xz_0_xz_xxxzz, g_xz_0_xz_xxyyy, g_xz_0_xz_xxyyz, g_xz_0_xz_xxyzz, g_xz_0_xz_xxzzz, g_xz_0_xz_xyyyy, g_xz_0_xz_xyyyz, g_xz_0_xz_xyyzz, g_xz_0_xz_xyzzz, g_xz_0_xz_xzzzz, g_xz_0_xz_yyyyy, g_xz_0_xz_yyyyz, g_xz_0_xz_yyyzz, g_xz_0_xz_yyzzz, g_xz_0_xz_yzzzz, g_xz_0_xz_zzzzz, g_xz_0_z_xxxxx, g_xz_0_z_xxxxxx, g_xz_0_z_xxxxxy, g_xz_0_z_xxxxxz, g_xz_0_z_xxxxy, g_xz_0_z_xxxxyy, g_xz_0_z_xxxxyz, g_xz_0_z_xxxxz, g_xz_0_z_xxxxzz, g_xz_0_z_xxxyy, g_xz_0_z_xxxyyy, g_xz_0_z_xxxyyz, g_xz_0_z_xxxyz, g_xz_0_z_xxxyzz, g_xz_0_z_xxxzz, g_xz_0_z_xxxzzz, g_xz_0_z_xxyyy, g_xz_0_z_xxyyyy, g_xz_0_z_xxyyyz, g_xz_0_z_xxyyz, g_xz_0_z_xxyyzz, g_xz_0_z_xxyzz, g_xz_0_z_xxyzzz, g_xz_0_z_xxzzz, g_xz_0_z_xxzzzz, g_xz_0_z_xyyyy, g_xz_0_z_xyyyyy, g_xz_0_z_xyyyyz, g_xz_0_z_xyyyz, g_xz_0_z_xyyyzz, g_xz_0_z_xyyzz, g_xz_0_z_xyyzzz, g_xz_0_z_xyzzz, g_xz_0_z_xyzzzz, g_xz_0_z_xzzzz, g_xz_0_z_xzzzzz, g_xz_0_z_yyyyy, g_xz_0_z_yyyyz, g_xz_0_z_yyyzz, g_xz_0_z_yyzzz, g_xz_0_z_yzzzz, g_xz_0_z_zzzzz, g_z_0_z_xxxxx, g_z_0_z_xxxxy, g_z_0_z_xxxxz, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxzz, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyzz, g_z_0_z_xxzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyzz, g_z_0_z_xyzzz, g_z_0_z_xzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyzz, g_z_0_z_yyzzz, g_z_0_z_yzzzz, g_z_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xz_xxxxx[k] = -g_z_0_z_xxxxx[k] - g_xz_0_z_xxxxx[k] * ab_x + g_xz_0_z_xxxxxx[k];

                g_xz_0_xz_xxxxy[k] = -g_z_0_z_xxxxy[k] - g_xz_0_z_xxxxy[k] * ab_x + g_xz_0_z_xxxxxy[k];

                g_xz_0_xz_xxxxz[k] = -g_z_0_z_xxxxz[k] - g_xz_0_z_xxxxz[k] * ab_x + g_xz_0_z_xxxxxz[k];

                g_xz_0_xz_xxxyy[k] = -g_z_0_z_xxxyy[k] - g_xz_0_z_xxxyy[k] * ab_x + g_xz_0_z_xxxxyy[k];

                g_xz_0_xz_xxxyz[k] = -g_z_0_z_xxxyz[k] - g_xz_0_z_xxxyz[k] * ab_x + g_xz_0_z_xxxxyz[k];

                g_xz_0_xz_xxxzz[k] = -g_z_0_z_xxxzz[k] - g_xz_0_z_xxxzz[k] * ab_x + g_xz_0_z_xxxxzz[k];

                g_xz_0_xz_xxyyy[k] = -g_z_0_z_xxyyy[k] - g_xz_0_z_xxyyy[k] * ab_x + g_xz_0_z_xxxyyy[k];

                g_xz_0_xz_xxyyz[k] = -g_z_0_z_xxyyz[k] - g_xz_0_z_xxyyz[k] * ab_x + g_xz_0_z_xxxyyz[k];

                g_xz_0_xz_xxyzz[k] = -g_z_0_z_xxyzz[k] - g_xz_0_z_xxyzz[k] * ab_x + g_xz_0_z_xxxyzz[k];

                g_xz_0_xz_xxzzz[k] = -g_z_0_z_xxzzz[k] - g_xz_0_z_xxzzz[k] * ab_x + g_xz_0_z_xxxzzz[k];

                g_xz_0_xz_xyyyy[k] = -g_z_0_z_xyyyy[k] - g_xz_0_z_xyyyy[k] * ab_x + g_xz_0_z_xxyyyy[k];

                g_xz_0_xz_xyyyz[k] = -g_z_0_z_xyyyz[k] - g_xz_0_z_xyyyz[k] * ab_x + g_xz_0_z_xxyyyz[k];

                g_xz_0_xz_xyyzz[k] = -g_z_0_z_xyyzz[k] - g_xz_0_z_xyyzz[k] * ab_x + g_xz_0_z_xxyyzz[k];

                g_xz_0_xz_xyzzz[k] = -g_z_0_z_xyzzz[k] - g_xz_0_z_xyzzz[k] * ab_x + g_xz_0_z_xxyzzz[k];

                g_xz_0_xz_xzzzz[k] = -g_z_0_z_xzzzz[k] - g_xz_0_z_xzzzz[k] * ab_x + g_xz_0_z_xxzzzz[k];

                g_xz_0_xz_yyyyy[k] = -g_z_0_z_yyyyy[k] - g_xz_0_z_yyyyy[k] * ab_x + g_xz_0_z_xyyyyy[k];

                g_xz_0_xz_yyyyz[k] = -g_z_0_z_yyyyz[k] - g_xz_0_z_yyyyz[k] * ab_x + g_xz_0_z_xyyyyz[k];

                g_xz_0_xz_yyyzz[k] = -g_z_0_z_yyyzz[k] - g_xz_0_z_yyyzz[k] * ab_x + g_xz_0_z_xyyyzz[k];

                g_xz_0_xz_yyzzz[k] = -g_z_0_z_yyzzz[k] - g_xz_0_z_yyzzz[k] * ab_x + g_xz_0_z_xyyzzz[k];

                g_xz_0_xz_yzzzz[k] = -g_z_0_z_yzzzz[k] - g_xz_0_z_yzzzz[k] * ab_x + g_xz_0_z_xyzzzz[k];

                g_xz_0_xz_zzzzz[k] = -g_z_0_z_zzzzz[k] - g_xz_0_z_zzzzz[k] * ab_x + g_xz_0_z_xzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 315 * ccomps * dcomps);

            auto g_xz_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 316 * ccomps * dcomps);

            auto g_xz_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 317 * ccomps * dcomps);

            auto g_xz_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 318 * ccomps * dcomps);

            auto g_xz_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 319 * ccomps * dcomps);

            auto g_xz_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 320 * ccomps * dcomps);

            auto g_xz_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 321 * ccomps * dcomps);

            auto g_xz_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 322 * ccomps * dcomps);

            auto g_xz_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 323 * ccomps * dcomps);

            auto g_xz_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 324 * ccomps * dcomps);

            auto g_xz_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 325 * ccomps * dcomps);

            auto g_xz_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 326 * ccomps * dcomps);

            auto g_xz_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 327 * ccomps * dcomps);

            auto g_xz_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 328 * ccomps * dcomps);

            auto g_xz_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 329 * ccomps * dcomps);

            auto g_xz_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 330 * ccomps * dcomps);

            auto g_xz_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 331 * ccomps * dcomps);

            auto g_xz_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 332 * ccomps * dcomps);

            auto g_xz_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 333 * ccomps * dcomps);

            auto g_xz_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 334 * ccomps * dcomps);

            auto g_xz_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_y_xxxxx, g_xz_0_y_xxxxxy, g_xz_0_y_xxxxy, g_xz_0_y_xxxxyy, g_xz_0_y_xxxxyz, g_xz_0_y_xxxxz, g_xz_0_y_xxxyy, g_xz_0_y_xxxyyy, g_xz_0_y_xxxyyz, g_xz_0_y_xxxyz, g_xz_0_y_xxxyzz, g_xz_0_y_xxxzz, g_xz_0_y_xxyyy, g_xz_0_y_xxyyyy, g_xz_0_y_xxyyyz, g_xz_0_y_xxyyz, g_xz_0_y_xxyyzz, g_xz_0_y_xxyzz, g_xz_0_y_xxyzzz, g_xz_0_y_xxzzz, g_xz_0_y_xyyyy, g_xz_0_y_xyyyyy, g_xz_0_y_xyyyyz, g_xz_0_y_xyyyz, g_xz_0_y_xyyyzz, g_xz_0_y_xyyzz, g_xz_0_y_xyyzzz, g_xz_0_y_xyzzz, g_xz_0_y_xyzzzz, g_xz_0_y_xzzzz, g_xz_0_y_yyyyy, g_xz_0_y_yyyyyy, g_xz_0_y_yyyyyz, g_xz_0_y_yyyyz, g_xz_0_y_yyyyzz, g_xz_0_y_yyyzz, g_xz_0_y_yyyzzz, g_xz_0_y_yyzzz, g_xz_0_y_yyzzzz, g_xz_0_y_yzzzz, g_xz_0_y_yzzzzz, g_xz_0_y_zzzzz, g_xz_0_yy_xxxxx, g_xz_0_yy_xxxxy, g_xz_0_yy_xxxxz, g_xz_0_yy_xxxyy, g_xz_0_yy_xxxyz, g_xz_0_yy_xxxzz, g_xz_0_yy_xxyyy, g_xz_0_yy_xxyyz, g_xz_0_yy_xxyzz, g_xz_0_yy_xxzzz, g_xz_0_yy_xyyyy, g_xz_0_yy_xyyyz, g_xz_0_yy_xyyzz, g_xz_0_yy_xyzzz, g_xz_0_yy_xzzzz, g_xz_0_yy_yyyyy, g_xz_0_yy_yyyyz, g_xz_0_yy_yyyzz, g_xz_0_yy_yyzzz, g_xz_0_yy_yzzzz, g_xz_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yy_xxxxx[k] = -g_xz_0_y_xxxxx[k] * ab_y + g_xz_0_y_xxxxxy[k];

                g_xz_0_yy_xxxxy[k] = -g_xz_0_y_xxxxy[k] * ab_y + g_xz_0_y_xxxxyy[k];

                g_xz_0_yy_xxxxz[k] = -g_xz_0_y_xxxxz[k] * ab_y + g_xz_0_y_xxxxyz[k];

                g_xz_0_yy_xxxyy[k] = -g_xz_0_y_xxxyy[k] * ab_y + g_xz_0_y_xxxyyy[k];

                g_xz_0_yy_xxxyz[k] = -g_xz_0_y_xxxyz[k] * ab_y + g_xz_0_y_xxxyyz[k];

                g_xz_0_yy_xxxzz[k] = -g_xz_0_y_xxxzz[k] * ab_y + g_xz_0_y_xxxyzz[k];

                g_xz_0_yy_xxyyy[k] = -g_xz_0_y_xxyyy[k] * ab_y + g_xz_0_y_xxyyyy[k];

                g_xz_0_yy_xxyyz[k] = -g_xz_0_y_xxyyz[k] * ab_y + g_xz_0_y_xxyyyz[k];

                g_xz_0_yy_xxyzz[k] = -g_xz_0_y_xxyzz[k] * ab_y + g_xz_0_y_xxyyzz[k];

                g_xz_0_yy_xxzzz[k] = -g_xz_0_y_xxzzz[k] * ab_y + g_xz_0_y_xxyzzz[k];

                g_xz_0_yy_xyyyy[k] = -g_xz_0_y_xyyyy[k] * ab_y + g_xz_0_y_xyyyyy[k];

                g_xz_0_yy_xyyyz[k] = -g_xz_0_y_xyyyz[k] * ab_y + g_xz_0_y_xyyyyz[k];

                g_xz_0_yy_xyyzz[k] = -g_xz_0_y_xyyzz[k] * ab_y + g_xz_0_y_xyyyzz[k];

                g_xz_0_yy_xyzzz[k] = -g_xz_0_y_xyzzz[k] * ab_y + g_xz_0_y_xyyzzz[k];

                g_xz_0_yy_xzzzz[k] = -g_xz_0_y_xzzzz[k] * ab_y + g_xz_0_y_xyzzzz[k];

                g_xz_0_yy_yyyyy[k] = -g_xz_0_y_yyyyy[k] * ab_y + g_xz_0_y_yyyyyy[k];

                g_xz_0_yy_yyyyz[k] = -g_xz_0_y_yyyyz[k] * ab_y + g_xz_0_y_yyyyyz[k];

                g_xz_0_yy_yyyzz[k] = -g_xz_0_y_yyyzz[k] * ab_y + g_xz_0_y_yyyyzz[k];

                g_xz_0_yy_yyzzz[k] = -g_xz_0_y_yyzzz[k] * ab_y + g_xz_0_y_yyyzzz[k];

                g_xz_0_yy_yzzzz[k] = -g_xz_0_y_yzzzz[k] * ab_y + g_xz_0_y_yyzzzz[k];

                g_xz_0_yy_zzzzz[k] = -g_xz_0_y_zzzzz[k] * ab_y + g_xz_0_y_yzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 339 * ccomps * dcomps);

            auto g_xz_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 341 * ccomps * dcomps);

            auto g_xz_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 344 * ccomps * dcomps);

            auto g_xz_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 347 * ccomps * dcomps);

            auto g_xz_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 349 * ccomps * dcomps);

            auto g_xz_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 353 * ccomps * dcomps);

            auto g_xz_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yz_xxxxx, g_xz_0_yz_xxxxy, g_xz_0_yz_xxxxz, g_xz_0_yz_xxxyy, g_xz_0_yz_xxxyz, g_xz_0_yz_xxxzz, g_xz_0_yz_xxyyy, g_xz_0_yz_xxyyz, g_xz_0_yz_xxyzz, g_xz_0_yz_xxzzz, g_xz_0_yz_xyyyy, g_xz_0_yz_xyyyz, g_xz_0_yz_xyyzz, g_xz_0_yz_xyzzz, g_xz_0_yz_xzzzz, g_xz_0_yz_yyyyy, g_xz_0_yz_yyyyz, g_xz_0_yz_yyyzz, g_xz_0_yz_yyzzz, g_xz_0_yz_yzzzz, g_xz_0_yz_zzzzz, g_xz_0_z_xxxxx, g_xz_0_z_xxxxxy, g_xz_0_z_xxxxy, g_xz_0_z_xxxxyy, g_xz_0_z_xxxxyz, g_xz_0_z_xxxxz, g_xz_0_z_xxxyy, g_xz_0_z_xxxyyy, g_xz_0_z_xxxyyz, g_xz_0_z_xxxyz, g_xz_0_z_xxxyzz, g_xz_0_z_xxxzz, g_xz_0_z_xxyyy, g_xz_0_z_xxyyyy, g_xz_0_z_xxyyyz, g_xz_0_z_xxyyz, g_xz_0_z_xxyyzz, g_xz_0_z_xxyzz, g_xz_0_z_xxyzzz, g_xz_0_z_xxzzz, g_xz_0_z_xyyyy, g_xz_0_z_xyyyyy, g_xz_0_z_xyyyyz, g_xz_0_z_xyyyz, g_xz_0_z_xyyyzz, g_xz_0_z_xyyzz, g_xz_0_z_xyyzzz, g_xz_0_z_xyzzz, g_xz_0_z_xyzzzz, g_xz_0_z_xzzzz, g_xz_0_z_yyyyy, g_xz_0_z_yyyyyy, g_xz_0_z_yyyyyz, g_xz_0_z_yyyyz, g_xz_0_z_yyyyzz, g_xz_0_z_yyyzz, g_xz_0_z_yyyzzz, g_xz_0_z_yyzzz, g_xz_0_z_yyzzzz, g_xz_0_z_yzzzz, g_xz_0_z_yzzzzz, g_xz_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yz_xxxxx[k] = -g_xz_0_z_xxxxx[k] * ab_y + g_xz_0_z_xxxxxy[k];

                g_xz_0_yz_xxxxy[k] = -g_xz_0_z_xxxxy[k] * ab_y + g_xz_0_z_xxxxyy[k];

                g_xz_0_yz_xxxxz[k] = -g_xz_0_z_xxxxz[k] * ab_y + g_xz_0_z_xxxxyz[k];

                g_xz_0_yz_xxxyy[k] = -g_xz_0_z_xxxyy[k] * ab_y + g_xz_0_z_xxxyyy[k];

                g_xz_0_yz_xxxyz[k] = -g_xz_0_z_xxxyz[k] * ab_y + g_xz_0_z_xxxyyz[k];

                g_xz_0_yz_xxxzz[k] = -g_xz_0_z_xxxzz[k] * ab_y + g_xz_0_z_xxxyzz[k];

                g_xz_0_yz_xxyyy[k] = -g_xz_0_z_xxyyy[k] * ab_y + g_xz_0_z_xxyyyy[k];

                g_xz_0_yz_xxyyz[k] = -g_xz_0_z_xxyyz[k] * ab_y + g_xz_0_z_xxyyyz[k];

                g_xz_0_yz_xxyzz[k] = -g_xz_0_z_xxyzz[k] * ab_y + g_xz_0_z_xxyyzz[k];

                g_xz_0_yz_xxzzz[k] = -g_xz_0_z_xxzzz[k] * ab_y + g_xz_0_z_xxyzzz[k];

                g_xz_0_yz_xyyyy[k] = -g_xz_0_z_xyyyy[k] * ab_y + g_xz_0_z_xyyyyy[k];

                g_xz_0_yz_xyyyz[k] = -g_xz_0_z_xyyyz[k] * ab_y + g_xz_0_z_xyyyyz[k];

                g_xz_0_yz_xyyzz[k] = -g_xz_0_z_xyyzz[k] * ab_y + g_xz_0_z_xyyyzz[k];

                g_xz_0_yz_xyzzz[k] = -g_xz_0_z_xyzzz[k] * ab_y + g_xz_0_z_xyyzzz[k];

                g_xz_0_yz_xzzzz[k] = -g_xz_0_z_xzzzz[k] * ab_y + g_xz_0_z_xyzzzz[k];

                g_xz_0_yz_yyyyy[k] = -g_xz_0_z_yyyyy[k] * ab_y + g_xz_0_z_yyyyyy[k];

                g_xz_0_yz_yyyyz[k] = -g_xz_0_z_yyyyz[k] * ab_y + g_xz_0_z_yyyyyz[k];

                g_xz_0_yz_yyyzz[k] = -g_xz_0_z_yyyzz[k] * ab_y + g_xz_0_z_yyyyzz[k];

                g_xz_0_yz_yyzzz[k] = -g_xz_0_z_yyzzz[k] * ab_y + g_xz_0_z_yyyzzz[k];

                g_xz_0_yz_yzzzz[k] = -g_xz_0_z_yzzzz[k] * ab_y + g_xz_0_z_yyzzzz[k];

                g_xz_0_yz_zzzzz[k] = -g_xz_0_z_zzzzz[k] * ab_y + g_xz_0_z_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 359 * ccomps * dcomps);

            auto g_xz_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 363 * ccomps * dcomps);

            auto g_xz_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 365 * ccomps * dcomps);

            auto g_xz_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 369 * ccomps * dcomps);

            auto g_xz_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 371 * ccomps * dcomps);

            auto g_xz_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 374 * ccomps * dcomps);

            auto g_xz_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxxxx, g_x_0_z_xxxxy, g_x_0_z_xxxxz, g_x_0_z_xxxyy, g_x_0_z_xxxyz, g_x_0_z_xxxzz, g_x_0_z_xxyyy, g_x_0_z_xxyyz, g_x_0_z_xxyzz, g_x_0_z_xxzzz, g_x_0_z_xyyyy, g_x_0_z_xyyyz, g_x_0_z_xyyzz, g_x_0_z_xyzzz, g_x_0_z_xzzzz, g_x_0_z_yyyyy, g_x_0_z_yyyyz, g_x_0_z_yyyzz, g_x_0_z_yyzzz, g_x_0_z_yzzzz, g_x_0_z_zzzzz, g_xz_0_z_xxxxx, g_xz_0_z_xxxxxz, g_xz_0_z_xxxxy, g_xz_0_z_xxxxyz, g_xz_0_z_xxxxz, g_xz_0_z_xxxxzz, g_xz_0_z_xxxyy, g_xz_0_z_xxxyyz, g_xz_0_z_xxxyz, g_xz_0_z_xxxyzz, g_xz_0_z_xxxzz, g_xz_0_z_xxxzzz, g_xz_0_z_xxyyy, g_xz_0_z_xxyyyz, g_xz_0_z_xxyyz, g_xz_0_z_xxyyzz, g_xz_0_z_xxyzz, g_xz_0_z_xxyzzz, g_xz_0_z_xxzzz, g_xz_0_z_xxzzzz, g_xz_0_z_xyyyy, g_xz_0_z_xyyyyz, g_xz_0_z_xyyyz, g_xz_0_z_xyyyzz, g_xz_0_z_xyyzz, g_xz_0_z_xyyzzz, g_xz_0_z_xyzzz, g_xz_0_z_xyzzzz, g_xz_0_z_xzzzz, g_xz_0_z_xzzzzz, g_xz_0_z_yyyyy, g_xz_0_z_yyyyyz, g_xz_0_z_yyyyz, g_xz_0_z_yyyyzz, g_xz_0_z_yyyzz, g_xz_0_z_yyyzzz, g_xz_0_z_yyzzz, g_xz_0_z_yyzzzz, g_xz_0_z_yzzzz, g_xz_0_z_yzzzzz, g_xz_0_z_zzzzz, g_xz_0_z_zzzzzz, g_xz_0_zz_xxxxx, g_xz_0_zz_xxxxy, g_xz_0_zz_xxxxz, g_xz_0_zz_xxxyy, g_xz_0_zz_xxxyz, g_xz_0_zz_xxxzz, g_xz_0_zz_xxyyy, g_xz_0_zz_xxyyz, g_xz_0_zz_xxyzz, g_xz_0_zz_xxzzz, g_xz_0_zz_xyyyy, g_xz_0_zz_xyyyz, g_xz_0_zz_xyyzz, g_xz_0_zz_xyzzz, g_xz_0_zz_xzzzz, g_xz_0_zz_yyyyy, g_xz_0_zz_yyyyz, g_xz_0_zz_yyyzz, g_xz_0_zz_yyzzz, g_xz_0_zz_yzzzz, g_xz_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zz_xxxxx[k] = -g_x_0_z_xxxxx[k] - g_xz_0_z_xxxxx[k] * ab_z + g_xz_0_z_xxxxxz[k];

                g_xz_0_zz_xxxxy[k] = -g_x_0_z_xxxxy[k] - g_xz_0_z_xxxxy[k] * ab_z + g_xz_0_z_xxxxyz[k];

                g_xz_0_zz_xxxxz[k] = -g_x_0_z_xxxxz[k] - g_xz_0_z_xxxxz[k] * ab_z + g_xz_0_z_xxxxzz[k];

                g_xz_0_zz_xxxyy[k] = -g_x_0_z_xxxyy[k] - g_xz_0_z_xxxyy[k] * ab_z + g_xz_0_z_xxxyyz[k];

                g_xz_0_zz_xxxyz[k] = -g_x_0_z_xxxyz[k] - g_xz_0_z_xxxyz[k] * ab_z + g_xz_0_z_xxxyzz[k];

                g_xz_0_zz_xxxzz[k] = -g_x_0_z_xxxzz[k] - g_xz_0_z_xxxzz[k] * ab_z + g_xz_0_z_xxxzzz[k];

                g_xz_0_zz_xxyyy[k] = -g_x_0_z_xxyyy[k] - g_xz_0_z_xxyyy[k] * ab_z + g_xz_0_z_xxyyyz[k];

                g_xz_0_zz_xxyyz[k] = -g_x_0_z_xxyyz[k] - g_xz_0_z_xxyyz[k] * ab_z + g_xz_0_z_xxyyzz[k];

                g_xz_0_zz_xxyzz[k] = -g_x_0_z_xxyzz[k] - g_xz_0_z_xxyzz[k] * ab_z + g_xz_0_z_xxyzzz[k];

                g_xz_0_zz_xxzzz[k] = -g_x_0_z_xxzzz[k] - g_xz_0_z_xxzzz[k] * ab_z + g_xz_0_z_xxzzzz[k];

                g_xz_0_zz_xyyyy[k] = -g_x_0_z_xyyyy[k] - g_xz_0_z_xyyyy[k] * ab_z + g_xz_0_z_xyyyyz[k];

                g_xz_0_zz_xyyyz[k] = -g_x_0_z_xyyyz[k] - g_xz_0_z_xyyyz[k] * ab_z + g_xz_0_z_xyyyzz[k];

                g_xz_0_zz_xyyzz[k] = -g_x_0_z_xyyzz[k] - g_xz_0_z_xyyzz[k] * ab_z + g_xz_0_z_xyyzzz[k];

                g_xz_0_zz_xyzzz[k] = -g_x_0_z_xyzzz[k] - g_xz_0_z_xyzzz[k] * ab_z + g_xz_0_z_xyzzzz[k];

                g_xz_0_zz_xzzzz[k] = -g_x_0_z_xzzzz[k] - g_xz_0_z_xzzzz[k] * ab_z + g_xz_0_z_xzzzzz[k];

                g_xz_0_zz_yyyyy[k] = -g_x_0_z_yyyyy[k] - g_xz_0_z_yyyyy[k] * ab_z + g_xz_0_z_yyyyyz[k];

                g_xz_0_zz_yyyyz[k] = -g_x_0_z_yyyyz[k] - g_xz_0_z_yyyyz[k] * ab_z + g_xz_0_z_yyyyzz[k];

                g_xz_0_zz_yyyzz[k] = -g_x_0_z_yyyzz[k] - g_xz_0_z_yyyzz[k] * ab_z + g_xz_0_z_yyyzzz[k];

                g_xz_0_zz_yyzzz[k] = -g_x_0_z_yyzzz[k] - g_xz_0_z_yyzzz[k] * ab_z + g_xz_0_z_yyzzzz[k];

                g_xz_0_zz_yzzzz[k] = -g_x_0_z_yzzzz[k] - g_xz_0_z_yzzzz[k] * ab_z + g_xz_0_z_yzzzzz[k];

                g_xz_0_zz_zzzzz[k] = -g_x_0_z_zzzzz[k] - g_xz_0_z_zzzzz[k] * ab_z + g_xz_0_z_zzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 378 * ccomps * dcomps);

            auto g_yy_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 379 * ccomps * dcomps);

            auto g_yy_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 380 * ccomps * dcomps);

            auto g_yy_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 381 * ccomps * dcomps);

            auto g_yy_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 382 * ccomps * dcomps);

            auto g_yy_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 383 * ccomps * dcomps);

            auto g_yy_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 384 * ccomps * dcomps);

            auto g_yy_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 385 * ccomps * dcomps);

            auto g_yy_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 386 * ccomps * dcomps);

            auto g_yy_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 387 * ccomps * dcomps);

            auto g_yy_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 388 * ccomps * dcomps);

            auto g_yy_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 389 * ccomps * dcomps);

            auto g_yy_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 390 * ccomps * dcomps);

            auto g_yy_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 391 * ccomps * dcomps);

            auto g_yy_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 392 * ccomps * dcomps);

            auto g_yy_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 393 * ccomps * dcomps);

            auto g_yy_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 394 * ccomps * dcomps);

            auto g_yy_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 395 * ccomps * dcomps);

            auto g_yy_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 396 * ccomps * dcomps);

            auto g_yy_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 397 * ccomps * dcomps);

            auto g_yy_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_x_xxxxx, g_yy_0_x_xxxxxx, g_yy_0_x_xxxxxy, g_yy_0_x_xxxxxz, g_yy_0_x_xxxxy, g_yy_0_x_xxxxyy, g_yy_0_x_xxxxyz, g_yy_0_x_xxxxz, g_yy_0_x_xxxxzz, g_yy_0_x_xxxyy, g_yy_0_x_xxxyyy, g_yy_0_x_xxxyyz, g_yy_0_x_xxxyz, g_yy_0_x_xxxyzz, g_yy_0_x_xxxzz, g_yy_0_x_xxxzzz, g_yy_0_x_xxyyy, g_yy_0_x_xxyyyy, g_yy_0_x_xxyyyz, g_yy_0_x_xxyyz, g_yy_0_x_xxyyzz, g_yy_0_x_xxyzz, g_yy_0_x_xxyzzz, g_yy_0_x_xxzzz, g_yy_0_x_xxzzzz, g_yy_0_x_xyyyy, g_yy_0_x_xyyyyy, g_yy_0_x_xyyyyz, g_yy_0_x_xyyyz, g_yy_0_x_xyyyzz, g_yy_0_x_xyyzz, g_yy_0_x_xyyzzz, g_yy_0_x_xyzzz, g_yy_0_x_xyzzzz, g_yy_0_x_xzzzz, g_yy_0_x_xzzzzz, g_yy_0_x_yyyyy, g_yy_0_x_yyyyz, g_yy_0_x_yyyzz, g_yy_0_x_yyzzz, g_yy_0_x_yzzzz, g_yy_0_x_zzzzz, g_yy_0_xx_xxxxx, g_yy_0_xx_xxxxy, g_yy_0_xx_xxxxz, g_yy_0_xx_xxxyy, g_yy_0_xx_xxxyz, g_yy_0_xx_xxxzz, g_yy_0_xx_xxyyy, g_yy_0_xx_xxyyz, g_yy_0_xx_xxyzz, g_yy_0_xx_xxzzz, g_yy_0_xx_xyyyy, g_yy_0_xx_xyyyz, g_yy_0_xx_xyyzz, g_yy_0_xx_xyzzz, g_yy_0_xx_xzzzz, g_yy_0_xx_yyyyy, g_yy_0_xx_yyyyz, g_yy_0_xx_yyyzz, g_yy_0_xx_yyzzz, g_yy_0_xx_yzzzz, g_yy_0_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xx_xxxxx[k] = -g_yy_0_x_xxxxx[k] * ab_x + g_yy_0_x_xxxxxx[k];

                g_yy_0_xx_xxxxy[k] = -g_yy_0_x_xxxxy[k] * ab_x + g_yy_0_x_xxxxxy[k];

                g_yy_0_xx_xxxxz[k] = -g_yy_0_x_xxxxz[k] * ab_x + g_yy_0_x_xxxxxz[k];

                g_yy_0_xx_xxxyy[k] = -g_yy_0_x_xxxyy[k] * ab_x + g_yy_0_x_xxxxyy[k];

                g_yy_0_xx_xxxyz[k] = -g_yy_0_x_xxxyz[k] * ab_x + g_yy_0_x_xxxxyz[k];

                g_yy_0_xx_xxxzz[k] = -g_yy_0_x_xxxzz[k] * ab_x + g_yy_0_x_xxxxzz[k];

                g_yy_0_xx_xxyyy[k] = -g_yy_0_x_xxyyy[k] * ab_x + g_yy_0_x_xxxyyy[k];

                g_yy_0_xx_xxyyz[k] = -g_yy_0_x_xxyyz[k] * ab_x + g_yy_0_x_xxxyyz[k];

                g_yy_0_xx_xxyzz[k] = -g_yy_0_x_xxyzz[k] * ab_x + g_yy_0_x_xxxyzz[k];

                g_yy_0_xx_xxzzz[k] = -g_yy_0_x_xxzzz[k] * ab_x + g_yy_0_x_xxxzzz[k];

                g_yy_0_xx_xyyyy[k] = -g_yy_0_x_xyyyy[k] * ab_x + g_yy_0_x_xxyyyy[k];

                g_yy_0_xx_xyyyz[k] = -g_yy_0_x_xyyyz[k] * ab_x + g_yy_0_x_xxyyyz[k];

                g_yy_0_xx_xyyzz[k] = -g_yy_0_x_xyyzz[k] * ab_x + g_yy_0_x_xxyyzz[k];

                g_yy_0_xx_xyzzz[k] = -g_yy_0_x_xyzzz[k] * ab_x + g_yy_0_x_xxyzzz[k];

                g_yy_0_xx_xzzzz[k] = -g_yy_0_x_xzzzz[k] * ab_x + g_yy_0_x_xxzzzz[k];

                g_yy_0_xx_yyyyy[k] = -g_yy_0_x_yyyyy[k] * ab_x + g_yy_0_x_xyyyyy[k];

                g_yy_0_xx_yyyyz[k] = -g_yy_0_x_yyyyz[k] * ab_x + g_yy_0_x_xyyyyz[k];

                g_yy_0_xx_yyyzz[k] = -g_yy_0_x_yyyzz[k] * ab_x + g_yy_0_x_xyyyzz[k];

                g_yy_0_xx_yyzzz[k] = -g_yy_0_x_yyzzz[k] * ab_x + g_yy_0_x_xyyzzz[k];

                g_yy_0_xx_yzzzz[k] = -g_yy_0_x_yzzzz[k] * ab_x + g_yy_0_x_xyzzzz[k];

                g_yy_0_xx_zzzzz[k] = -g_yy_0_x_zzzzz[k] * ab_x + g_yy_0_x_xzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 399 * ccomps * dcomps);

            auto g_yy_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 400 * ccomps * dcomps);

            auto g_yy_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 401 * ccomps * dcomps);

            auto g_yy_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 402 * ccomps * dcomps);

            auto g_yy_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 403 * ccomps * dcomps);

            auto g_yy_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 404 * ccomps * dcomps);

            auto g_yy_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 405 * ccomps * dcomps);

            auto g_yy_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 406 * ccomps * dcomps);

            auto g_yy_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 407 * ccomps * dcomps);

            auto g_yy_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 408 * ccomps * dcomps);

            auto g_yy_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 409 * ccomps * dcomps);

            auto g_yy_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 410 * ccomps * dcomps);

            auto g_yy_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 411 * ccomps * dcomps);

            auto g_yy_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 412 * ccomps * dcomps);

            auto g_yy_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 413 * ccomps * dcomps);

            auto g_yy_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 414 * ccomps * dcomps);

            auto g_yy_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 415 * ccomps * dcomps);

            auto g_yy_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 416 * ccomps * dcomps);

            auto g_yy_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 417 * ccomps * dcomps);

            auto g_yy_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 418 * ccomps * dcomps);

            auto g_yy_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xy_xxxxx, g_yy_0_xy_xxxxy, g_yy_0_xy_xxxxz, g_yy_0_xy_xxxyy, g_yy_0_xy_xxxyz, g_yy_0_xy_xxxzz, g_yy_0_xy_xxyyy, g_yy_0_xy_xxyyz, g_yy_0_xy_xxyzz, g_yy_0_xy_xxzzz, g_yy_0_xy_xyyyy, g_yy_0_xy_xyyyz, g_yy_0_xy_xyyzz, g_yy_0_xy_xyzzz, g_yy_0_xy_xzzzz, g_yy_0_xy_yyyyy, g_yy_0_xy_yyyyz, g_yy_0_xy_yyyzz, g_yy_0_xy_yyzzz, g_yy_0_xy_yzzzz, g_yy_0_xy_zzzzz, g_yy_0_y_xxxxx, g_yy_0_y_xxxxxx, g_yy_0_y_xxxxxy, g_yy_0_y_xxxxxz, g_yy_0_y_xxxxy, g_yy_0_y_xxxxyy, g_yy_0_y_xxxxyz, g_yy_0_y_xxxxz, g_yy_0_y_xxxxzz, g_yy_0_y_xxxyy, g_yy_0_y_xxxyyy, g_yy_0_y_xxxyyz, g_yy_0_y_xxxyz, g_yy_0_y_xxxyzz, g_yy_0_y_xxxzz, g_yy_0_y_xxxzzz, g_yy_0_y_xxyyy, g_yy_0_y_xxyyyy, g_yy_0_y_xxyyyz, g_yy_0_y_xxyyz, g_yy_0_y_xxyyzz, g_yy_0_y_xxyzz, g_yy_0_y_xxyzzz, g_yy_0_y_xxzzz, g_yy_0_y_xxzzzz, g_yy_0_y_xyyyy, g_yy_0_y_xyyyyy, g_yy_0_y_xyyyyz, g_yy_0_y_xyyyz, g_yy_0_y_xyyyzz, g_yy_0_y_xyyzz, g_yy_0_y_xyyzzz, g_yy_0_y_xyzzz, g_yy_0_y_xyzzzz, g_yy_0_y_xzzzz, g_yy_0_y_xzzzzz, g_yy_0_y_yyyyy, g_yy_0_y_yyyyz, g_yy_0_y_yyyzz, g_yy_0_y_yyzzz, g_yy_0_y_yzzzz, g_yy_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xy_xxxxx[k] = -g_yy_0_y_xxxxx[k] * ab_x + g_yy_0_y_xxxxxx[k];

                g_yy_0_xy_xxxxy[k] = -g_yy_0_y_xxxxy[k] * ab_x + g_yy_0_y_xxxxxy[k];

                g_yy_0_xy_xxxxz[k] = -g_yy_0_y_xxxxz[k] * ab_x + g_yy_0_y_xxxxxz[k];

                g_yy_0_xy_xxxyy[k] = -g_yy_0_y_xxxyy[k] * ab_x + g_yy_0_y_xxxxyy[k];

                g_yy_0_xy_xxxyz[k] = -g_yy_0_y_xxxyz[k] * ab_x + g_yy_0_y_xxxxyz[k];

                g_yy_0_xy_xxxzz[k] = -g_yy_0_y_xxxzz[k] * ab_x + g_yy_0_y_xxxxzz[k];

                g_yy_0_xy_xxyyy[k] = -g_yy_0_y_xxyyy[k] * ab_x + g_yy_0_y_xxxyyy[k];

                g_yy_0_xy_xxyyz[k] = -g_yy_0_y_xxyyz[k] * ab_x + g_yy_0_y_xxxyyz[k];

                g_yy_0_xy_xxyzz[k] = -g_yy_0_y_xxyzz[k] * ab_x + g_yy_0_y_xxxyzz[k];

                g_yy_0_xy_xxzzz[k] = -g_yy_0_y_xxzzz[k] * ab_x + g_yy_0_y_xxxzzz[k];

                g_yy_0_xy_xyyyy[k] = -g_yy_0_y_xyyyy[k] * ab_x + g_yy_0_y_xxyyyy[k];

                g_yy_0_xy_xyyyz[k] = -g_yy_0_y_xyyyz[k] * ab_x + g_yy_0_y_xxyyyz[k];

                g_yy_0_xy_xyyzz[k] = -g_yy_0_y_xyyzz[k] * ab_x + g_yy_0_y_xxyyzz[k];

                g_yy_0_xy_xyzzz[k] = -g_yy_0_y_xyzzz[k] * ab_x + g_yy_0_y_xxyzzz[k];

                g_yy_0_xy_xzzzz[k] = -g_yy_0_y_xzzzz[k] * ab_x + g_yy_0_y_xxzzzz[k];

                g_yy_0_xy_yyyyy[k] = -g_yy_0_y_yyyyy[k] * ab_x + g_yy_0_y_xyyyyy[k];

                g_yy_0_xy_yyyyz[k] = -g_yy_0_y_yyyyz[k] * ab_x + g_yy_0_y_xyyyyz[k];

                g_yy_0_xy_yyyzz[k] = -g_yy_0_y_yyyzz[k] * ab_x + g_yy_0_y_xyyyzz[k];

                g_yy_0_xy_yyzzz[k] = -g_yy_0_y_yyzzz[k] * ab_x + g_yy_0_y_xyyzzz[k];

                g_yy_0_xy_yzzzz[k] = -g_yy_0_y_yzzzz[k] * ab_x + g_yy_0_y_xyzzzz[k];

                g_yy_0_xy_zzzzz[k] = -g_yy_0_y_zzzzz[k] * ab_x + g_yy_0_y_xzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 420 * ccomps * dcomps);

            auto g_yy_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 421 * ccomps * dcomps);

            auto g_yy_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 422 * ccomps * dcomps);

            auto g_yy_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 423 * ccomps * dcomps);

            auto g_yy_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 424 * ccomps * dcomps);

            auto g_yy_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 425 * ccomps * dcomps);

            auto g_yy_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 426 * ccomps * dcomps);

            auto g_yy_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 427 * ccomps * dcomps);

            auto g_yy_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 428 * ccomps * dcomps);

            auto g_yy_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 429 * ccomps * dcomps);

            auto g_yy_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 430 * ccomps * dcomps);

            auto g_yy_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 431 * ccomps * dcomps);

            auto g_yy_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 432 * ccomps * dcomps);

            auto g_yy_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 433 * ccomps * dcomps);

            auto g_yy_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 434 * ccomps * dcomps);

            auto g_yy_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 435 * ccomps * dcomps);

            auto g_yy_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 436 * ccomps * dcomps);

            auto g_yy_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 437 * ccomps * dcomps);

            auto g_yy_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 438 * ccomps * dcomps);

            auto g_yy_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 439 * ccomps * dcomps);

            auto g_yy_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 440 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xz_xxxxx, g_yy_0_xz_xxxxy, g_yy_0_xz_xxxxz, g_yy_0_xz_xxxyy, g_yy_0_xz_xxxyz, g_yy_0_xz_xxxzz, g_yy_0_xz_xxyyy, g_yy_0_xz_xxyyz, g_yy_0_xz_xxyzz, g_yy_0_xz_xxzzz, g_yy_0_xz_xyyyy, g_yy_0_xz_xyyyz, g_yy_0_xz_xyyzz, g_yy_0_xz_xyzzz, g_yy_0_xz_xzzzz, g_yy_0_xz_yyyyy, g_yy_0_xz_yyyyz, g_yy_0_xz_yyyzz, g_yy_0_xz_yyzzz, g_yy_0_xz_yzzzz, g_yy_0_xz_zzzzz, g_yy_0_z_xxxxx, g_yy_0_z_xxxxxx, g_yy_0_z_xxxxxy, g_yy_0_z_xxxxxz, g_yy_0_z_xxxxy, g_yy_0_z_xxxxyy, g_yy_0_z_xxxxyz, g_yy_0_z_xxxxz, g_yy_0_z_xxxxzz, g_yy_0_z_xxxyy, g_yy_0_z_xxxyyy, g_yy_0_z_xxxyyz, g_yy_0_z_xxxyz, g_yy_0_z_xxxyzz, g_yy_0_z_xxxzz, g_yy_0_z_xxxzzz, g_yy_0_z_xxyyy, g_yy_0_z_xxyyyy, g_yy_0_z_xxyyyz, g_yy_0_z_xxyyz, g_yy_0_z_xxyyzz, g_yy_0_z_xxyzz, g_yy_0_z_xxyzzz, g_yy_0_z_xxzzz, g_yy_0_z_xxzzzz, g_yy_0_z_xyyyy, g_yy_0_z_xyyyyy, g_yy_0_z_xyyyyz, g_yy_0_z_xyyyz, g_yy_0_z_xyyyzz, g_yy_0_z_xyyzz, g_yy_0_z_xyyzzz, g_yy_0_z_xyzzz, g_yy_0_z_xyzzzz, g_yy_0_z_xzzzz, g_yy_0_z_xzzzzz, g_yy_0_z_yyyyy, g_yy_0_z_yyyyz, g_yy_0_z_yyyzz, g_yy_0_z_yyzzz, g_yy_0_z_yzzzz, g_yy_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xz_xxxxx[k] = -g_yy_0_z_xxxxx[k] * ab_x + g_yy_0_z_xxxxxx[k];

                g_yy_0_xz_xxxxy[k] = -g_yy_0_z_xxxxy[k] * ab_x + g_yy_0_z_xxxxxy[k];

                g_yy_0_xz_xxxxz[k] = -g_yy_0_z_xxxxz[k] * ab_x + g_yy_0_z_xxxxxz[k];

                g_yy_0_xz_xxxyy[k] = -g_yy_0_z_xxxyy[k] * ab_x + g_yy_0_z_xxxxyy[k];

                g_yy_0_xz_xxxyz[k] = -g_yy_0_z_xxxyz[k] * ab_x + g_yy_0_z_xxxxyz[k];

                g_yy_0_xz_xxxzz[k] = -g_yy_0_z_xxxzz[k] * ab_x + g_yy_0_z_xxxxzz[k];

                g_yy_0_xz_xxyyy[k] = -g_yy_0_z_xxyyy[k] * ab_x + g_yy_0_z_xxxyyy[k];

                g_yy_0_xz_xxyyz[k] = -g_yy_0_z_xxyyz[k] * ab_x + g_yy_0_z_xxxyyz[k];

                g_yy_0_xz_xxyzz[k] = -g_yy_0_z_xxyzz[k] * ab_x + g_yy_0_z_xxxyzz[k];

                g_yy_0_xz_xxzzz[k] = -g_yy_0_z_xxzzz[k] * ab_x + g_yy_0_z_xxxzzz[k];

                g_yy_0_xz_xyyyy[k] = -g_yy_0_z_xyyyy[k] * ab_x + g_yy_0_z_xxyyyy[k];

                g_yy_0_xz_xyyyz[k] = -g_yy_0_z_xyyyz[k] * ab_x + g_yy_0_z_xxyyyz[k];

                g_yy_0_xz_xyyzz[k] = -g_yy_0_z_xyyzz[k] * ab_x + g_yy_0_z_xxyyzz[k];

                g_yy_0_xz_xyzzz[k] = -g_yy_0_z_xyzzz[k] * ab_x + g_yy_0_z_xxyzzz[k];

                g_yy_0_xz_xzzzz[k] = -g_yy_0_z_xzzzz[k] * ab_x + g_yy_0_z_xxzzzz[k];

                g_yy_0_xz_yyyyy[k] = -g_yy_0_z_yyyyy[k] * ab_x + g_yy_0_z_xyyyyy[k];

                g_yy_0_xz_yyyyz[k] = -g_yy_0_z_yyyyz[k] * ab_x + g_yy_0_z_xyyyyz[k];

                g_yy_0_xz_yyyzz[k] = -g_yy_0_z_yyyzz[k] * ab_x + g_yy_0_z_xyyyzz[k];

                g_yy_0_xz_yyzzz[k] = -g_yy_0_z_yyzzz[k] * ab_x + g_yy_0_z_xyyzzz[k];

                g_yy_0_xz_yzzzz[k] = -g_yy_0_z_yzzzz[k] * ab_x + g_yy_0_z_xyzzzz[k];

                g_yy_0_xz_zzzzz[k] = -g_yy_0_z_zzzzz[k] * ab_x + g_yy_0_z_xzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 441 * ccomps * dcomps);

            auto g_yy_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 442 * ccomps * dcomps);

            auto g_yy_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 443 * ccomps * dcomps);

            auto g_yy_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 444 * ccomps * dcomps);

            auto g_yy_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 445 * ccomps * dcomps);

            auto g_yy_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 446 * ccomps * dcomps);

            auto g_yy_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 447 * ccomps * dcomps);

            auto g_yy_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 448 * ccomps * dcomps);

            auto g_yy_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 449 * ccomps * dcomps);

            auto g_yy_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 450 * ccomps * dcomps);

            auto g_yy_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 451 * ccomps * dcomps);

            auto g_yy_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 452 * ccomps * dcomps);

            auto g_yy_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 453 * ccomps * dcomps);

            auto g_yy_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 454 * ccomps * dcomps);

            auto g_yy_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 455 * ccomps * dcomps);

            auto g_yy_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 456 * ccomps * dcomps);

            auto g_yy_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 457 * ccomps * dcomps);

            auto g_yy_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 458 * ccomps * dcomps);

            auto g_yy_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 459 * ccomps * dcomps);

            auto g_yy_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 460 * ccomps * dcomps);

            auto g_yy_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxxxx, g_y_0_y_xxxxy, g_y_0_y_xxxxz, g_y_0_y_xxxyy, g_y_0_y_xxxyz, g_y_0_y_xxxzz, g_y_0_y_xxyyy, g_y_0_y_xxyyz, g_y_0_y_xxyzz, g_y_0_y_xxzzz, g_y_0_y_xyyyy, g_y_0_y_xyyyz, g_y_0_y_xyyzz, g_y_0_y_xyzzz, g_y_0_y_xzzzz, g_y_0_y_yyyyy, g_y_0_y_yyyyz, g_y_0_y_yyyzz, g_y_0_y_yyzzz, g_y_0_y_yzzzz, g_y_0_y_zzzzz, g_yy_0_y_xxxxx, g_yy_0_y_xxxxxy, g_yy_0_y_xxxxy, g_yy_0_y_xxxxyy, g_yy_0_y_xxxxyz, g_yy_0_y_xxxxz, g_yy_0_y_xxxyy, g_yy_0_y_xxxyyy, g_yy_0_y_xxxyyz, g_yy_0_y_xxxyz, g_yy_0_y_xxxyzz, g_yy_0_y_xxxzz, g_yy_0_y_xxyyy, g_yy_0_y_xxyyyy, g_yy_0_y_xxyyyz, g_yy_0_y_xxyyz, g_yy_0_y_xxyyzz, g_yy_0_y_xxyzz, g_yy_0_y_xxyzzz, g_yy_0_y_xxzzz, g_yy_0_y_xyyyy, g_yy_0_y_xyyyyy, g_yy_0_y_xyyyyz, g_yy_0_y_xyyyz, g_yy_0_y_xyyyzz, g_yy_0_y_xyyzz, g_yy_0_y_xyyzzz, g_yy_0_y_xyzzz, g_yy_0_y_xyzzzz, g_yy_0_y_xzzzz, g_yy_0_y_yyyyy, g_yy_0_y_yyyyyy, g_yy_0_y_yyyyyz, g_yy_0_y_yyyyz, g_yy_0_y_yyyyzz, g_yy_0_y_yyyzz, g_yy_0_y_yyyzzz, g_yy_0_y_yyzzz, g_yy_0_y_yyzzzz, g_yy_0_y_yzzzz, g_yy_0_y_yzzzzz, g_yy_0_y_zzzzz, g_yy_0_yy_xxxxx, g_yy_0_yy_xxxxy, g_yy_0_yy_xxxxz, g_yy_0_yy_xxxyy, g_yy_0_yy_xxxyz, g_yy_0_yy_xxxzz, g_yy_0_yy_xxyyy, g_yy_0_yy_xxyyz, g_yy_0_yy_xxyzz, g_yy_0_yy_xxzzz, g_yy_0_yy_xyyyy, g_yy_0_yy_xyyyz, g_yy_0_yy_xyyzz, g_yy_0_yy_xyzzz, g_yy_0_yy_xzzzz, g_yy_0_yy_yyyyy, g_yy_0_yy_yyyyz, g_yy_0_yy_yyyzz, g_yy_0_yy_yyzzz, g_yy_0_yy_yzzzz, g_yy_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yy_xxxxx[k] = -2.0 * g_y_0_y_xxxxx[k] - g_yy_0_y_xxxxx[k] * ab_y + g_yy_0_y_xxxxxy[k];

                g_yy_0_yy_xxxxy[k] = -2.0 * g_y_0_y_xxxxy[k] - g_yy_0_y_xxxxy[k] * ab_y + g_yy_0_y_xxxxyy[k];

                g_yy_0_yy_xxxxz[k] = -2.0 * g_y_0_y_xxxxz[k] - g_yy_0_y_xxxxz[k] * ab_y + g_yy_0_y_xxxxyz[k];

                g_yy_0_yy_xxxyy[k] = -2.0 * g_y_0_y_xxxyy[k] - g_yy_0_y_xxxyy[k] * ab_y + g_yy_0_y_xxxyyy[k];

                g_yy_0_yy_xxxyz[k] = -2.0 * g_y_0_y_xxxyz[k] - g_yy_0_y_xxxyz[k] * ab_y + g_yy_0_y_xxxyyz[k];

                g_yy_0_yy_xxxzz[k] = -2.0 * g_y_0_y_xxxzz[k] - g_yy_0_y_xxxzz[k] * ab_y + g_yy_0_y_xxxyzz[k];

                g_yy_0_yy_xxyyy[k] = -2.0 * g_y_0_y_xxyyy[k] - g_yy_0_y_xxyyy[k] * ab_y + g_yy_0_y_xxyyyy[k];

                g_yy_0_yy_xxyyz[k] = -2.0 * g_y_0_y_xxyyz[k] - g_yy_0_y_xxyyz[k] * ab_y + g_yy_0_y_xxyyyz[k];

                g_yy_0_yy_xxyzz[k] = -2.0 * g_y_0_y_xxyzz[k] - g_yy_0_y_xxyzz[k] * ab_y + g_yy_0_y_xxyyzz[k];

                g_yy_0_yy_xxzzz[k] = -2.0 * g_y_0_y_xxzzz[k] - g_yy_0_y_xxzzz[k] * ab_y + g_yy_0_y_xxyzzz[k];

                g_yy_0_yy_xyyyy[k] = -2.0 * g_y_0_y_xyyyy[k] - g_yy_0_y_xyyyy[k] * ab_y + g_yy_0_y_xyyyyy[k];

                g_yy_0_yy_xyyyz[k] = -2.0 * g_y_0_y_xyyyz[k] - g_yy_0_y_xyyyz[k] * ab_y + g_yy_0_y_xyyyyz[k];

                g_yy_0_yy_xyyzz[k] = -2.0 * g_y_0_y_xyyzz[k] - g_yy_0_y_xyyzz[k] * ab_y + g_yy_0_y_xyyyzz[k];

                g_yy_0_yy_xyzzz[k] = -2.0 * g_y_0_y_xyzzz[k] - g_yy_0_y_xyzzz[k] * ab_y + g_yy_0_y_xyyzzz[k];

                g_yy_0_yy_xzzzz[k] = -2.0 * g_y_0_y_xzzzz[k] - g_yy_0_y_xzzzz[k] * ab_y + g_yy_0_y_xyzzzz[k];

                g_yy_0_yy_yyyyy[k] = -2.0 * g_y_0_y_yyyyy[k] - g_yy_0_y_yyyyy[k] * ab_y + g_yy_0_y_yyyyyy[k];

                g_yy_0_yy_yyyyz[k] = -2.0 * g_y_0_y_yyyyz[k] - g_yy_0_y_yyyyz[k] * ab_y + g_yy_0_y_yyyyyz[k];

                g_yy_0_yy_yyyzz[k] = -2.0 * g_y_0_y_yyyzz[k] - g_yy_0_y_yyyzz[k] * ab_y + g_yy_0_y_yyyyzz[k];

                g_yy_0_yy_yyzzz[k] = -2.0 * g_y_0_y_yyzzz[k] - g_yy_0_y_yyzzz[k] * ab_y + g_yy_0_y_yyyzzz[k];

                g_yy_0_yy_yzzzz[k] = -2.0 * g_y_0_y_yzzzz[k] - g_yy_0_y_yzzzz[k] * ab_y + g_yy_0_y_yyzzzz[k];

                g_yy_0_yy_zzzzz[k] = -2.0 * g_y_0_y_zzzzz[k] - g_yy_0_y_zzzzz[k] * ab_y + g_yy_0_y_yzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 462 * ccomps * dcomps);

            auto g_yy_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 463 * ccomps * dcomps);

            auto g_yy_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 464 * ccomps * dcomps);

            auto g_yy_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 465 * ccomps * dcomps);

            auto g_yy_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 466 * ccomps * dcomps);

            auto g_yy_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 467 * ccomps * dcomps);

            auto g_yy_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 468 * ccomps * dcomps);

            auto g_yy_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 469 * ccomps * dcomps);

            auto g_yy_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 470 * ccomps * dcomps);

            auto g_yy_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 471 * ccomps * dcomps);

            auto g_yy_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 472 * ccomps * dcomps);

            auto g_yy_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 473 * ccomps * dcomps);

            auto g_yy_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 474 * ccomps * dcomps);

            auto g_yy_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 475 * ccomps * dcomps);

            auto g_yy_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 476 * ccomps * dcomps);

            auto g_yy_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 477 * ccomps * dcomps);

            auto g_yy_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 478 * ccomps * dcomps);

            auto g_yy_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 479 * ccomps * dcomps);

            auto g_yy_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 480 * ccomps * dcomps);

            auto g_yy_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 481 * ccomps * dcomps);

            auto g_yy_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 482 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_y_xxxxx, g_yy_0_y_xxxxxz, g_yy_0_y_xxxxy, g_yy_0_y_xxxxyz, g_yy_0_y_xxxxz, g_yy_0_y_xxxxzz, g_yy_0_y_xxxyy, g_yy_0_y_xxxyyz, g_yy_0_y_xxxyz, g_yy_0_y_xxxyzz, g_yy_0_y_xxxzz, g_yy_0_y_xxxzzz, g_yy_0_y_xxyyy, g_yy_0_y_xxyyyz, g_yy_0_y_xxyyz, g_yy_0_y_xxyyzz, g_yy_0_y_xxyzz, g_yy_0_y_xxyzzz, g_yy_0_y_xxzzz, g_yy_0_y_xxzzzz, g_yy_0_y_xyyyy, g_yy_0_y_xyyyyz, g_yy_0_y_xyyyz, g_yy_0_y_xyyyzz, g_yy_0_y_xyyzz, g_yy_0_y_xyyzzz, g_yy_0_y_xyzzz, g_yy_0_y_xyzzzz, g_yy_0_y_xzzzz, g_yy_0_y_xzzzzz, g_yy_0_y_yyyyy, g_yy_0_y_yyyyyz, g_yy_0_y_yyyyz, g_yy_0_y_yyyyzz, g_yy_0_y_yyyzz, g_yy_0_y_yyyzzz, g_yy_0_y_yyzzz, g_yy_0_y_yyzzzz, g_yy_0_y_yzzzz, g_yy_0_y_yzzzzz, g_yy_0_y_zzzzz, g_yy_0_y_zzzzzz, g_yy_0_yz_xxxxx, g_yy_0_yz_xxxxy, g_yy_0_yz_xxxxz, g_yy_0_yz_xxxyy, g_yy_0_yz_xxxyz, g_yy_0_yz_xxxzz, g_yy_0_yz_xxyyy, g_yy_0_yz_xxyyz, g_yy_0_yz_xxyzz, g_yy_0_yz_xxzzz, g_yy_0_yz_xyyyy, g_yy_0_yz_xyyyz, g_yy_0_yz_xyyzz, g_yy_0_yz_xyzzz, g_yy_0_yz_xzzzz, g_yy_0_yz_yyyyy, g_yy_0_yz_yyyyz, g_yy_0_yz_yyyzz, g_yy_0_yz_yyzzz, g_yy_0_yz_yzzzz, g_yy_0_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yz_xxxxx[k] = -g_yy_0_y_xxxxx[k] * ab_z + g_yy_0_y_xxxxxz[k];

                g_yy_0_yz_xxxxy[k] = -g_yy_0_y_xxxxy[k] * ab_z + g_yy_0_y_xxxxyz[k];

                g_yy_0_yz_xxxxz[k] = -g_yy_0_y_xxxxz[k] * ab_z + g_yy_0_y_xxxxzz[k];

                g_yy_0_yz_xxxyy[k] = -g_yy_0_y_xxxyy[k] * ab_z + g_yy_0_y_xxxyyz[k];

                g_yy_0_yz_xxxyz[k] = -g_yy_0_y_xxxyz[k] * ab_z + g_yy_0_y_xxxyzz[k];

                g_yy_0_yz_xxxzz[k] = -g_yy_0_y_xxxzz[k] * ab_z + g_yy_0_y_xxxzzz[k];

                g_yy_0_yz_xxyyy[k] = -g_yy_0_y_xxyyy[k] * ab_z + g_yy_0_y_xxyyyz[k];

                g_yy_0_yz_xxyyz[k] = -g_yy_0_y_xxyyz[k] * ab_z + g_yy_0_y_xxyyzz[k];

                g_yy_0_yz_xxyzz[k] = -g_yy_0_y_xxyzz[k] * ab_z + g_yy_0_y_xxyzzz[k];

                g_yy_0_yz_xxzzz[k] = -g_yy_0_y_xxzzz[k] * ab_z + g_yy_0_y_xxzzzz[k];

                g_yy_0_yz_xyyyy[k] = -g_yy_0_y_xyyyy[k] * ab_z + g_yy_0_y_xyyyyz[k];

                g_yy_0_yz_xyyyz[k] = -g_yy_0_y_xyyyz[k] * ab_z + g_yy_0_y_xyyyzz[k];

                g_yy_0_yz_xyyzz[k] = -g_yy_0_y_xyyzz[k] * ab_z + g_yy_0_y_xyyzzz[k];

                g_yy_0_yz_xyzzz[k] = -g_yy_0_y_xyzzz[k] * ab_z + g_yy_0_y_xyzzzz[k];

                g_yy_0_yz_xzzzz[k] = -g_yy_0_y_xzzzz[k] * ab_z + g_yy_0_y_xzzzzz[k];

                g_yy_0_yz_yyyyy[k] = -g_yy_0_y_yyyyy[k] * ab_z + g_yy_0_y_yyyyyz[k];

                g_yy_0_yz_yyyyz[k] = -g_yy_0_y_yyyyz[k] * ab_z + g_yy_0_y_yyyyzz[k];

                g_yy_0_yz_yyyzz[k] = -g_yy_0_y_yyyzz[k] * ab_z + g_yy_0_y_yyyzzz[k];

                g_yy_0_yz_yyzzz[k] = -g_yy_0_y_yyzzz[k] * ab_z + g_yy_0_y_yyzzzz[k];

                g_yy_0_yz_yzzzz[k] = -g_yy_0_y_yzzzz[k] * ab_z + g_yy_0_y_yzzzzz[k];

                g_yy_0_yz_zzzzz[k] = -g_yy_0_y_zzzzz[k] * ab_z + g_yy_0_y_zzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 483 * ccomps * dcomps);

            auto g_yy_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 484 * ccomps * dcomps);

            auto g_yy_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 485 * ccomps * dcomps);

            auto g_yy_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 486 * ccomps * dcomps);

            auto g_yy_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 487 * ccomps * dcomps);

            auto g_yy_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 488 * ccomps * dcomps);

            auto g_yy_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 489 * ccomps * dcomps);

            auto g_yy_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 490 * ccomps * dcomps);

            auto g_yy_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 491 * ccomps * dcomps);

            auto g_yy_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 492 * ccomps * dcomps);

            auto g_yy_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 493 * ccomps * dcomps);

            auto g_yy_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 494 * ccomps * dcomps);

            auto g_yy_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 495 * ccomps * dcomps);

            auto g_yy_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 496 * ccomps * dcomps);

            auto g_yy_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 497 * ccomps * dcomps);

            auto g_yy_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 498 * ccomps * dcomps);

            auto g_yy_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 499 * ccomps * dcomps);

            auto g_yy_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 500 * ccomps * dcomps);

            auto g_yy_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 501 * ccomps * dcomps);

            auto g_yy_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 502 * ccomps * dcomps);

            auto g_yy_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_z_xxxxx, g_yy_0_z_xxxxxz, g_yy_0_z_xxxxy, g_yy_0_z_xxxxyz, g_yy_0_z_xxxxz, g_yy_0_z_xxxxzz, g_yy_0_z_xxxyy, g_yy_0_z_xxxyyz, g_yy_0_z_xxxyz, g_yy_0_z_xxxyzz, g_yy_0_z_xxxzz, g_yy_0_z_xxxzzz, g_yy_0_z_xxyyy, g_yy_0_z_xxyyyz, g_yy_0_z_xxyyz, g_yy_0_z_xxyyzz, g_yy_0_z_xxyzz, g_yy_0_z_xxyzzz, g_yy_0_z_xxzzz, g_yy_0_z_xxzzzz, g_yy_0_z_xyyyy, g_yy_0_z_xyyyyz, g_yy_0_z_xyyyz, g_yy_0_z_xyyyzz, g_yy_0_z_xyyzz, g_yy_0_z_xyyzzz, g_yy_0_z_xyzzz, g_yy_0_z_xyzzzz, g_yy_0_z_xzzzz, g_yy_0_z_xzzzzz, g_yy_0_z_yyyyy, g_yy_0_z_yyyyyz, g_yy_0_z_yyyyz, g_yy_0_z_yyyyzz, g_yy_0_z_yyyzz, g_yy_0_z_yyyzzz, g_yy_0_z_yyzzz, g_yy_0_z_yyzzzz, g_yy_0_z_yzzzz, g_yy_0_z_yzzzzz, g_yy_0_z_zzzzz, g_yy_0_z_zzzzzz, g_yy_0_zz_xxxxx, g_yy_0_zz_xxxxy, g_yy_0_zz_xxxxz, g_yy_0_zz_xxxyy, g_yy_0_zz_xxxyz, g_yy_0_zz_xxxzz, g_yy_0_zz_xxyyy, g_yy_0_zz_xxyyz, g_yy_0_zz_xxyzz, g_yy_0_zz_xxzzz, g_yy_0_zz_xyyyy, g_yy_0_zz_xyyyz, g_yy_0_zz_xyyzz, g_yy_0_zz_xyzzz, g_yy_0_zz_xzzzz, g_yy_0_zz_yyyyy, g_yy_0_zz_yyyyz, g_yy_0_zz_yyyzz, g_yy_0_zz_yyzzz, g_yy_0_zz_yzzzz, g_yy_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zz_xxxxx[k] = -g_yy_0_z_xxxxx[k] * ab_z + g_yy_0_z_xxxxxz[k];

                g_yy_0_zz_xxxxy[k] = -g_yy_0_z_xxxxy[k] * ab_z + g_yy_0_z_xxxxyz[k];

                g_yy_0_zz_xxxxz[k] = -g_yy_0_z_xxxxz[k] * ab_z + g_yy_0_z_xxxxzz[k];

                g_yy_0_zz_xxxyy[k] = -g_yy_0_z_xxxyy[k] * ab_z + g_yy_0_z_xxxyyz[k];

                g_yy_0_zz_xxxyz[k] = -g_yy_0_z_xxxyz[k] * ab_z + g_yy_0_z_xxxyzz[k];

                g_yy_0_zz_xxxzz[k] = -g_yy_0_z_xxxzz[k] * ab_z + g_yy_0_z_xxxzzz[k];

                g_yy_0_zz_xxyyy[k] = -g_yy_0_z_xxyyy[k] * ab_z + g_yy_0_z_xxyyyz[k];

                g_yy_0_zz_xxyyz[k] = -g_yy_0_z_xxyyz[k] * ab_z + g_yy_0_z_xxyyzz[k];

                g_yy_0_zz_xxyzz[k] = -g_yy_0_z_xxyzz[k] * ab_z + g_yy_0_z_xxyzzz[k];

                g_yy_0_zz_xxzzz[k] = -g_yy_0_z_xxzzz[k] * ab_z + g_yy_0_z_xxzzzz[k];

                g_yy_0_zz_xyyyy[k] = -g_yy_0_z_xyyyy[k] * ab_z + g_yy_0_z_xyyyyz[k];

                g_yy_0_zz_xyyyz[k] = -g_yy_0_z_xyyyz[k] * ab_z + g_yy_0_z_xyyyzz[k];

                g_yy_0_zz_xyyzz[k] = -g_yy_0_z_xyyzz[k] * ab_z + g_yy_0_z_xyyzzz[k];

                g_yy_0_zz_xyzzz[k] = -g_yy_0_z_xyzzz[k] * ab_z + g_yy_0_z_xyzzzz[k];

                g_yy_0_zz_xzzzz[k] = -g_yy_0_z_xzzzz[k] * ab_z + g_yy_0_z_xzzzzz[k];

                g_yy_0_zz_yyyyy[k] = -g_yy_0_z_yyyyy[k] * ab_z + g_yy_0_z_yyyyyz[k];

                g_yy_0_zz_yyyyz[k] = -g_yy_0_z_yyyyz[k] * ab_z + g_yy_0_z_yyyyzz[k];

                g_yy_0_zz_yyyzz[k] = -g_yy_0_z_yyyzz[k] * ab_z + g_yy_0_z_yyyzzz[k];

                g_yy_0_zz_yyzzz[k] = -g_yy_0_z_yyzzz[k] * ab_z + g_yy_0_z_yyzzzz[k];

                g_yy_0_zz_yzzzz[k] = -g_yy_0_z_yzzzz[k] * ab_z + g_yy_0_z_yzzzzz[k];

                g_yy_0_zz_zzzzz[k] = -g_yy_0_z_zzzzz[k] * ab_z + g_yy_0_z_zzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 504 * ccomps * dcomps);

            auto g_yz_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 505 * ccomps * dcomps);

            auto g_yz_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 506 * ccomps * dcomps);

            auto g_yz_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 507 * ccomps * dcomps);

            auto g_yz_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 508 * ccomps * dcomps);

            auto g_yz_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 509 * ccomps * dcomps);

            auto g_yz_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 510 * ccomps * dcomps);

            auto g_yz_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 511 * ccomps * dcomps);

            auto g_yz_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 512 * ccomps * dcomps);

            auto g_yz_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 513 * ccomps * dcomps);

            auto g_yz_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 514 * ccomps * dcomps);

            auto g_yz_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 515 * ccomps * dcomps);

            auto g_yz_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 516 * ccomps * dcomps);

            auto g_yz_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 517 * ccomps * dcomps);

            auto g_yz_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 518 * ccomps * dcomps);

            auto g_yz_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 519 * ccomps * dcomps);

            auto g_yz_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 520 * ccomps * dcomps);

            auto g_yz_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 521 * ccomps * dcomps);

            auto g_yz_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 522 * ccomps * dcomps);

            auto g_yz_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 523 * ccomps * dcomps);

            auto g_yz_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_x_xxxxx, g_yz_0_x_xxxxxx, g_yz_0_x_xxxxxy, g_yz_0_x_xxxxxz, g_yz_0_x_xxxxy, g_yz_0_x_xxxxyy, g_yz_0_x_xxxxyz, g_yz_0_x_xxxxz, g_yz_0_x_xxxxzz, g_yz_0_x_xxxyy, g_yz_0_x_xxxyyy, g_yz_0_x_xxxyyz, g_yz_0_x_xxxyz, g_yz_0_x_xxxyzz, g_yz_0_x_xxxzz, g_yz_0_x_xxxzzz, g_yz_0_x_xxyyy, g_yz_0_x_xxyyyy, g_yz_0_x_xxyyyz, g_yz_0_x_xxyyz, g_yz_0_x_xxyyzz, g_yz_0_x_xxyzz, g_yz_0_x_xxyzzz, g_yz_0_x_xxzzz, g_yz_0_x_xxzzzz, g_yz_0_x_xyyyy, g_yz_0_x_xyyyyy, g_yz_0_x_xyyyyz, g_yz_0_x_xyyyz, g_yz_0_x_xyyyzz, g_yz_0_x_xyyzz, g_yz_0_x_xyyzzz, g_yz_0_x_xyzzz, g_yz_0_x_xyzzzz, g_yz_0_x_xzzzz, g_yz_0_x_xzzzzz, g_yz_0_x_yyyyy, g_yz_0_x_yyyyz, g_yz_0_x_yyyzz, g_yz_0_x_yyzzz, g_yz_0_x_yzzzz, g_yz_0_x_zzzzz, g_yz_0_xx_xxxxx, g_yz_0_xx_xxxxy, g_yz_0_xx_xxxxz, g_yz_0_xx_xxxyy, g_yz_0_xx_xxxyz, g_yz_0_xx_xxxzz, g_yz_0_xx_xxyyy, g_yz_0_xx_xxyyz, g_yz_0_xx_xxyzz, g_yz_0_xx_xxzzz, g_yz_0_xx_xyyyy, g_yz_0_xx_xyyyz, g_yz_0_xx_xyyzz, g_yz_0_xx_xyzzz, g_yz_0_xx_xzzzz, g_yz_0_xx_yyyyy, g_yz_0_xx_yyyyz, g_yz_0_xx_yyyzz, g_yz_0_xx_yyzzz, g_yz_0_xx_yzzzz, g_yz_0_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xx_xxxxx[k] = -g_yz_0_x_xxxxx[k] * ab_x + g_yz_0_x_xxxxxx[k];

                g_yz_0_xx_xxxxy[k] = -g_yz_0_x_xxxxy[k] * ab_x + g_yz_0_x_xxxxxy[k];

                g_yz_0_xx_xxxxz[k] = -g_yz_0_x_xxxxz[k] * ab_x + g_yz_0_x_xxxxxz[k];

                g_yz_0_xx_xxxyy[k] = -g_yz_0_x_xxxyy[k] * ab_x + g_yz_0_x_xxxxyy[k];

                g_yz_0_xx_xxxyz[k] = -g_yz_0_x_xxxyz[k] * ab_x + g_yz_0_x_xxxxyz[k];

                g_yz_0_xx_xxxzz[k] = -g_yz_0_x_xxxzz[k] * ab_x + g_yz_0_x_xxxxzz[k];

                g_yz_0_xx_xxyyy[k] = -g_yz_0_x_xxyyy[k] * ab_x + g_yz_0_x_xxxyyy[k];

                g_yz_0_xx_xxyyz[k] = -g_yz_0_x_xxyyz[k] * ab_x + g_yz_0_x_xxxyyz[k];

                g_yz_0_xx_xxyzz[k] = -g_yz_0_x_xxyzz[k] * ab_x + g_yz_0_x_xxxyzz[k];

                g_yz_0_xx_xxzzz[k] = -g_yz_0_x_xxzzz[k] * ab_x + g_yz_0_x_xxxzzz[k];

                g_yz_0_xx_xyyyy[k] = -g_yz_0_x_xyyyy[k] * ab_x + g_yz_0_x_xxyyyy[k];

                g_yz_0_xx_xyyyz[k] = -g_yz_0_x_xyyyz[k] * ab_x + g_yz_0_x_xxyyyz[k];

                g_yz_0_xx_xyyzz[k] = -g_yz_0_x_xyyzz[k] * ab_x + g_yz_0_x_xxyyzz[k];

                g_yz_0_xx_xyzzz[k] = -g_yz_0_x_xyzzz[k] * ab_x + g_yz_0_x_xxyzzz[k];

                g_yz_0_xx_xzzzz[k] = -g_yz_0_x_xzzzz[k] * ab_x + g_yz_0_x_xxzzzz[k];

                g_yz_0_xx_yyyyy[k] = -g_yz_0_x_yyyyy[k] * ab_x + g_yz_0_x_xyyyyy[k];

                g_yz_0_xx_yyyyz[k] = -g_yz_0_x_yyyyz[k] * ab_x + g_yz_0_x_xyyyyz[k];

                g_yz_0_xx_yyyzz[k] = -g_yz_0_x_yyyzz[k] * ab_x + g_yz_0_x_xyyyzz[k];

                g_yz_0_xx_yyzzz[k] = -g_yz_0_x_yyzzz[k] * ab_x + g_yz_0_x_xyyzzz[k];

                g_yz_0_xx_yzzzz[k] = -g_yz_0_x_yzzzz[k] * ab_x + g_yz_0_x_xyzzzz[k];

                g_yz_0_xx_zzzzz[k] = -g_yz_0_x_zzzzz[k] * ab_x + g_yz_0_x_xzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 525 * ccomps * dcomps);

            auto g_yz_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 526 * ccomps * dcomps);

            auto g_yz_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 527 * ccomps * dcomps);

            auto g_yz_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 528 * ccomps * dcomps);

            auto g_yz_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 529 * ccomps * dcomps);

            auto g_yz_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 530 * ccomps * dcomps);

            auto g_yz_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 531 * ccomps * dcomps);

            auto g_yz_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 532 * ccomps * dcomps);

            auto g_yz_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 533 * ccomps * dcomps);

            auto g_yz_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 534 * ccomps * dcomps);

            auto g_yz_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 535 * ccomps * dcomps);

            auto g_yz_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 536 * ccomps * dcomps);

            auto g_yz_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 537 * ccomps * dcomps);

            auto g_yz_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 538 * ccomps * dcomps);

            auto g_yz_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 539 * ccomps * dcomps);

            auto g_yz_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 540 * ccomps * dcomps);

            auto g_yz_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 541 * ccomps * dcomps);

            auto g_yz_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 542 * ccomps * dcomps);

            auto g_yz_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 543 * ccomps * dcomps);

            auto g_yz_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 544 * ccomps * dcomps);

            auto g_yz_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xy_xxxxx, g_yz_0_xy_xxxxy, g_yz_0_xy_xxxxz, g_yz_0_xy_xxxyy, g_yz_0_xy_xxxyz, g_yz_0_xy_xxxzz, g_yz_0_xy_xxyyy, g_yz_0_xy_xxyyz, g_yz_0_xy_xxyzz, g_yz_0_xy_xxzzz, g_yz_0_xy_xyyyy, g_yz_0_xy_xyyyz, g_yz_0_xy_xyyzz, g_yz_0_xy_xyzzz, g_yz_0_xy_xzzzz, g_yz_0_xy_yyyyy, g_yz_0_xy_yyyyz, g_yz_0_xy_yyyzz, g_yz_0_xy_yyzzz, g_yz_0_xy_yzzzz, g_yz_0_xy_zzzzz, g_yz_0_y_xxxxx, g_yz_0_y_xxxxxx, g_yz_0_y_xxxxxy, g_yz_0_y_xxxxxz, g_yz_0_y_xxxxy, g_yz_0_y_xxxxyy, g_yz_0_y_xxxxyz, g_yz_0_y_xxxxz, g_yz_0_y_xxxxzz, g_yz_0_y_xxxyy, g_yz_0_y_xxxyyy, g_yz_0_y_xxxyyz, g_yz_0_y_xxxyz, g_yz_0_y_xxxyzz, g_yz_0_y_xxxzz, g_yz_0_y_xxxzzz, g_yz_0_y_xxyyy, g_yz_0_y_xxyyyy, g_yz_0_y_xxyyyz, g_yz_0_y_xxyyz, g_yz_0_y_xxyyzz, g_yz_0_y_xxyzz, g_yz_0_y_xxyzzz, g_yz_0_y_xxzzz, g_yz_0_y_xxzzzz, g_yz_0_y_xyyyy, g_yz_0_y_xyyyyy, g_yz_0_y_xyyyyz, g_yz_0_y_xyyyz, g_yz_0_y_xyyyzz, g_yz_0_y_xyyzz, g_yz_0_y_xyyzzz, g_yz_0_y_xyzzz, g_yz_0_y_xyzzzz, g_yz_0_y_xzzzz, g_yz_0_y_xzzzzz, g_yz_0_y_yyyyy, g_yz_0_y_yyyyz, g_yz_0_y_yyyzz, g_yz_0_y_yyzzz, g_yz_0_y_yzzzz, g_yz_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xy_xxxxx[k] = -g_yz_0_y_xxxxx[k] * ab_x + g_yz_0_y_xxxxxx[k];

                g_yz_0_xy_xxxxy[k] = -g_yz_0_y_xxxxy[k] * ab_x + g_yz_0_y_xxxxxy[k];

                g_yz_0_xy_xxxxz[k] = -g_yz_0_y_xxxxz[k] * ab_x + g_yz_0_y_xxxxxz[k];

                g_yz_0_xy_xxxyy[k] = -g_yz_0_y_xxxyy[k] * ab_x + g_yz_0_y_xxxxyy[k];

                g_yz_0_xy_xxxyz[k] = -g_yz_0_y_xxxyz[k] * ab_x + g_yz_0_y_xxxxyz[k];

                g_yz_0_xy_xxxzz[k] = -g_yz_0_y_xxxzz[k] * ab_x + g_yz_0_y_xxxxzz[k];

                g_yz_0_xy_xxyyy[k] = -g_yz_0_y_xxyyy[k] * ab_x + g_yz_0_y_xxxyyy[k];

                g_yz_0_xy_xxyyz[k] = -g_yz_0_y_xxyyz[k] * ab_x + g_yz_0_y_xxxyyz[k];

                g_yz_0_xy_xxyzz[k] = -g_yz_0_y_xxyzz[k] * ab_x + g_yz_0_y_xxxyzz[k];

                g_yz_0_xy_xxzzz[k] = -g_yz_0_y_xxzzz[k] * ab_x + g_yz_0_y_xxxzzz[k];

                g_yz_0_xy_xyyyy[k] = -g_yz_0_y_xyyyy[k] * ab_x + g_yz_0_y_xxyyyy[k];

                g_yz_0_xy_xyyyz[k] = -g_yz_0_y_xyyyz[k] * ab_x + g_yz_0_y_xxyyyz[k];

                g_yz_0_xy_xyyzz[k] = -g_yz_0_y_xyyzz[k] * ab_x + g_yz_0_y_xxyyzz[k];

                g_yz_0_xy_xyzzz[k] = -g_yz_0_y_xyzzz[k] * ab_x + g_yz_0_y_xxyzzz[k];

                g_yz_0_xy_xzzzz[k] = -g_yz_0_y_xzzzz[k] * ab_x + g_yz_0_y_xxzzzz[k];

                g_yz_0_xy_yyyyy[k] = -g_yz_0_y_yyyyy[k] * ab_x + g_yz_0_y_xyyyyy[k];

                g_yz_0_xy_yyyyz[k] = -g_yz_0_y_yyyyz[k] * ab_x + g_yz_0_y_xyyyyz[k];

                g_yz_0_xy_yyyzz[k] = -g_yz_0_y_yyyzz[k] * ab_x + g_yz_0_y_xyyyzz[k];

                g_yz_0_xy_yyzzz[k] = -g_yz_0_y_yyzzz[k] * ab_x + g_yz_0_y_xyyzzz[k];

                g_yz_0_xy_yzzzz[k] = -g_yz_0_y_yzzzz[k] * ab_x + g_yz_0_y_xyzzzz[k];

                g_yz_0_xy_zzzzz[k] = -g_yz_0_y_zzzzz[k] * ab_x + g_yz_0_y_xzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 546 * ccomps * dcomps);

            auto g_yz_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 547 * ccomps * dcomps);

            auto g_yz_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 548 * ccomps * dcomps);

            auto g_yz_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 549 * ccomps * dcomps);

            auto g_yz_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 550 * ccomps * dcomps);

            auto g_yz_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 551 * ccomps * dcomps);

            auto g_yz_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 552 * ccomps * dcomps);

            auto g_yz_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 553 * ccomps * dcomps);

            auto g_yz_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 554 * ccomps * dcomps);

            auto g_yz_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 555 * ccomps * dcomps);

            auto g_yz_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 556 * ccomps * dcomps);

            auto g_yz_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 557 * ccomps * dcomps);

            auto g_yz_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 558 * ccomps * dcomps);

            auto g_yz_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 559 * ccomps * dcomps);

            auto g_yz_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 560 * ccomps * dcomps);

            auto g_yz_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 561 * ccomps * dcomps);

            auto g_yz_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 562 * ccomps * dcomps);

            auto g_yz_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 563 * ccomps * dcomps);

            auto g_yz_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 564 * ccomps * dcomps);

            auto g_yz_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 565 * ccomps * dcomps);

            auto g_yz_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 566 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xz_xxxxx, g_yz_0_xz_xxxxy, g_yz_0_xz_xxxxz, g_yz_0_xz_xxxyy, g_yz_0_xz_xxxyz, g_yz_0_xz_xxxzz, g_yz_0_xz_xxyyy, g_yz_0_xz_xxyyz, g_yz_0_xz_xxyzz, g_yz_0_xz_xxzzz, g_yz_0_xz_xyyyy, g_yz_0_xz_xyyyz, g_yz_0_xz_xyyzz, g_yz_0_xz_xyzzz, g_yz_0_xz_xzzzz, g_yz_0_xz_yyyyy, g_yz_0_xz_yyyyz, g_yz_0_xz_yyyzz, g_yz_0_xz_yyzzz, g_yz_0_xz_yzzzz, g_yz_0_xz_zzzzz, g_yz_0_z_xxxxx, g_yz_0_z_xxxxxx, g_yz_0_z_xxxxxy, g_yz_0_z_xxxxxz, g_yz_0_z_xxxxy, g_yz_0_z_xxxxyy, g_yz_0_z_xxxxyz, g_yz_0_z_xxxxz, g_yz_0_z_xxxxzz, g_yz_0_z_xxxyy, g_yz_0_z_xxxyyy, g_yz_0_z_xxxyyz, g_yz_0_z_xxxyz, g_yz_0_z_xxxyzz, g_yz_0_z_xxxzz, g_yz_0_z_xxxzzz, g_yz_0_z_xxyyy, g_yz_0_z_xxyyyy, g_yz_0_z_xxyyyz, g_yz_0_z_xxyyz, g_yz_0_z_xxyyzz, g_yz_0_z_xxyzz, g_yz_0_z_xxyzzz, g_yz_0_z_xxzzz, g_yz_0_z_xxzzzz, g_yz_0_z_xyyyy, g_yz_0_z_xyyyyy, g_yz_0_z_xyyyyz, g_yz_0_z_xyyyz, g_yz_0_z_xyyyzz, g_yz_0_z_xyyzz, g_yz_0_z_xyyzzz, g_yz_0_z_xyzzz, g_yz_0_z_xyzzzz, g_yz_0_z_xzzzz, g_yz_0_z_xzzzzz, g_yz_0_z_yyyyy, g_yz_0_z_yyyyz, g_yz_0_z_yyyzz, g_yz_0_z_yyzzz, g_yz_0_z_yzzzz, g_yz_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xz_xxxxx[k] = -g_yz_0_z_xxxxx[k] * ab_x + g_yz_0_z_xxxxxx[k];

                g_yz_0_xz_xxxxy[k] = -g_yz_0_z_xxxxy[k] * ab_x + g_yz_0_z_xxxxxy[k];

                g_yz_0_xz_xxxxz[k] = -g_yz_0_z_xxxxz[k] * ab_x + g_yz_0_z_xxxxxz[k];

                g_yz_0_xz_xxxyy[k] = -g_yz_0_z_xxxyy[k] * ab_x + g_yz_0_z_xxxxyy[k];

                g_yz_0_xz_xxxyz[k] = -g_yz_0_z_xxxyz[k] * ab_x + g_yz_0_z_xxxxyz[k];

                g_yz_0_xz_xxxzz[k] = -g_yz_0_z_xxxzz[k] * ab_x + g_yz_0_z_xxxxzz[k];

                g_yz_0_xz_xxyyy[k] = -g_yz_0_z_xxyyy[k] * ab_x + g_yz_0_z_xxxyyy[k];

                g_yz_0_xz_xxyyz[k] = -g_yz_0_z_xxyyz[k] * ab_x + g_yz_0_z_xxxyyz[k];

                g_yz_0_xz_xxyzz[k] = -g_yz_0_z_xxyzz[k] * ab_x + g_yz_0_z_xxxyzz[k];

                g_yz_0_xz_xxzzz[k] = -g_yz_0_z_xxzzz[k] * ab_x + g_yz_0_z_xxxzzz[k];

                g_yz_0_xz_xyyyy[k] = -g_yz_0_z_xyyyy[k] * ab_x + g_yz_0_z_xxyyyy[k];

                g_yz_0_xz_xyyyz[k] = -g_yz_0_z_xyyyz[k] * ab_x + g_yz_0_z_xxyyyz[k];

                g_yz_0_xz_xyyzz[k] = -g_yz_0_z_xyyzz[k] * ab_x + g_yz_0_z_xxyyzz[k];

                g_yz_0_xz_xyzzz[k] = -g_yz_0_z_xyzzz[k] * ab_x + g_yz_0_z_xxyzzz[k];

                g_yz_0_xz_xzzzz[k] = -g_yz_0_z_xzzzz[k] * ab_x + g_yz_0_z_xxzzzz[k];

                g_yz_0_xz_yyyyy[k] = -g_yz_0_z_yyyyy[k] * ab_x + g_yz_0_z_xyyyyy[k];

                g_yz_0_xz_yyyyz[k] = -g_yz_0_z_yyyyz[k] * ab_x + g_yz_0_z_xyyyyz[k];

                g_yz_0_xz_yyyzz[k] = -g_yz_0_z_yyyzz[k] * ab_x + g_yz_0_z_xyyyzz[k];

                g_yz_0_xz_yyzzz[k] = -g_yz_0_z_yyzzz[k] * ab_x + g_yz_0_z_xyyzzz[k];

                g_yz_0_xz_yzzzz[k] = -g_yz_0_z_yzzzz[k] * ab_x + g_yz_0_z_xyzzzz[k];

                g_yz_0_xz_zzzzz[k] = -g_yz_0_z_zzzzz[k] * ab_x + g_yz_0_z_xzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 567 * ccomps * dcomps);

            auto g_yz_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 568 * ccomps * dcomps);

            auto g_yz_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 569 * ccomps * dcomps);

            auto g_yz_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 570 * ccomps * dcomps);

            auto g_yz_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 571 * ccomps * dcomps);

            auto g_yz_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 572 * ccomps * dcomps);

            auto g_yz_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 573 * ccomps * dcomps);

            auto g_yz_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 574 * ccomps * dcomps);

            auto g_yz_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 575 * ccomps * dcomps);

            auto g_yz_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 576 * ccomps * dcomps);

            auto g_yz_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 577 * ccomps * dcomps);

            auto g_yz_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 578 * ccomps * dcomps);

            auto g_yz_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 579 * ccomps * dcomps);

            auto g_yz_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 580 * ccomps * dcomps);

            auto g_yz_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 581 * ccomps * dcomps);

            auto g_yz_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 582 * ccomps * dcomps);

            auto g_yz_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 583 * ccomps * dcomps);

            auto g_yz_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 584 * ccomps * dcomps);

            auto g_yz_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 585 * ccomps * dcomps);

            auto g_yz_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 586 * ccomps * dcomps);

            auto g_yz_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_y_xxxxx, g_yz_0_y_xxxxxy, g_yz_0_y_xxxxy, g_yz_0_y_xxxxyy, g_yz_0_y_xxxxyz, g_yz_0_y_xxxxz, g_yz_0_y_xxxyy, g_yz_0_y_xxxyyy, g_yz_0_y_xxxyyz, g_yz_0_y_xxxyz, g_yz_0_y_xxxyzz, g_yz_0_y_xxxzz, g_yz_0_y_xxyyy, g_yz_0_y_xxyyyy, g_yz_0_y_xxyyyz, g_yz_0_y_xxyyz, g_yz_0_y_xxyyzz, g_yz_0_y_xxyzz, g_yz_0_y_xxyzzz, g_yz_0_y_xxzzz, g_yz_0_y_xyyyy, g_yz_0_y_xyyyyy, g_yz_0_y_xyyyyz, g_yz_0_y_xyyyz, g_yz_0_y_xyyyzz, g_yz_0_y_xyyzz, g_yz_0_y_xyyzzz, g_yz_0_y_xyzzz, g_yz_0_y_xyzzzz, g_yz_0_y_xzzzz, g_yz_0_y_yyyyy, g_yz_0_y_yyyyyy, g_yz_0_y_yyyyyz, g_yz_0_y_yyyyz, g_yz_0_y_yyyyzz, g_yz_0_y_yyyzz, g_yz_0_y_yyyzzz, g_yz_0_y_yyzzz, g_yz_0_y_yyzzzz, g_yz_0_y_yzzzz, g_yz_0_y_yzzzzz, g_yz_0_y_zzzzz, g_yz_0_yy_xxxxx, g_yz_0_yy_xxxxy, g_yz_0_yy_xxxxz, g_yz_0_yy_xxxyy, g_yz_0_yy_xxxyz, g_yz_0_yy_xxxzz, g_yz_0_yy_xxyyy, g_yz_0_yy_xxyyz, g_yz_0_yy_xxyzz, g_yz_0_yy_xxzzz, g_yz_0_yy_xyyyy, g_yz_0_yy_xyyyz, g_yz_0_yy_xyyzz, g_yz_0_yy_xyzzz, g_yz_0_yy_xzzzz, g_yz_0_yy_yyyyy, g_yz_0_yy_yyyyz, g_yz_0_yy_yyyzz, g_yz_0_yy_yyzzz, g_yz_0_yy_yzzzz, g_yz_0_yy_zzzzz, g_z_0_y_xxxxx, g_z_0_y_xxxxy, g_z_0_y_xxxxz, g_z_0_y_xxxyy, g_z_0_y_xxxyz, g_z_0_y_xxxzz, g_z_0_y_xxyyy, g_z_0_y_xxyyz, g_z_0_y_xxyzz, g_z_0_y_xxzzz, g_z_0_y_xyyyy, g_z_0_y_xyyyz, g_z_0_y_xyyzz, g_z_0_y_xyzzz, g_z_0_y_xzzzz, g_z_0_y_yyyyy, g_z_0_y_yyyyz, g_z_0_y_yyyzz, g_z_0_y_yyzzz, g_z_0_y_yzzzz, g_z_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yy_xxxxx[k] = -g_z_0_y_xxxxx[k] - g_yz_0_y_xxxxx[k] * ab_y + g_yz_0_y_xxxxxy[k];

                g_yz_0_yy_xxxxy[k] = -g_z_0_y_xxxxy[k] - g_yz_0_y_xxxxy[k] * ab_y + g_yz_0_y_xxxxyy[k];

                g_yz_0_yy_xxxxz[k] = -g_z_0_y_xxxxz[k] - g_yz_0_y_xxxxz[k] * ab_y + g_yz_0_y_xxxxyz[k];

                g_yz_0_yy_xxxyy[k] = -g_z_0_y_xxxyy[k] - g_yz_0_y_xxxyy[k] * ab_y + g_yz_0_y_xxxyyy[k];

                g_yz_0_yy_xxxyz[k] = -g_z_0_y_xxxyz[k] - g_yz_0_y_xxxyz[k] * ab_y + g_yz_0_y_xxxyyz[k];

                g_yz_0_yy_xxxzz[k] = -g_z_0_y_xxxzz[k] - g_yz_0_y_xxxzz[k] * ab_y + g_yz_0_y_xxxyzz[k];

                g_yz_0_yy_xxyyy[k] = -g_z_0_y_xxyyy[k] - g_yz_0_y_xxyyy[k] * ab_y + g_yz_0_y_xxyyyy[k];

                g_yz_0_yy_xxyyz[k] = -g_z_0_y_xxyyz[k] - g_yz_0_y_xxyyz[k] * ab_y + g_yz_0_y_xxyyyz[k];

                g_yz_0_yy_xxyzz[k] = -g_z_0_y_xxyzz[k] - g_yz_0_y_xxyzz[k] * ab_y + g_yz_0_y_xxyyzz[k];

                g_yz_0_yy_xxzzz[k] = -g_z_0_y_xxzzz[k] - g_yz_0_y_xxzzz[k] * ab_y + g_yz_0_y_xxyzzz[k];

                g_yz_0_yy_xyyyy[k] = -g_z_0_y_xyyyy[k] - g_yz_0_y_xyyyy[k] * ab_y + g_yz_0_y_xyyyyy[k];

                g_yz_0_yy_xyyyz[k] = -g_z_0_y_xyyyz[k] - g_yz_0_y_xyyyz[k] * ab_y + g_yz_0_y_xyyyyz[k];

                g_yz_0_yy_xyyzz[k] = -g_z_0_y_xyyzz[k] - g_yz_0_y_xyyzz[k] * ab_y + g_yz_0_y_xyyyzz[k];

                g_yz_0_yy_xyzzz[k] = -g_z_0_y_xyzzz[k] - g_yz_0_y_xyzzz[k] * ab_y + g_yz_0_y_xyyzzz[k];

                g_yz_0_yy_xzzzz[k] = -g_z_0_y_xzzzz[k] - g_yz_0_y_xzzzz[k] * ab_y + g_yz_0_y_xyzzzz[k];

                g_yz_0_yy_yyyyy[k] = -g_z_0_y_yyyyy[k] - g_yz_0_y_yyyyy[k] * ab_y + g_yz_0_y_yyyyyy[k];

                g_yz_0_yy_yyyyz[k] = -g_z_0_y_yyyyz[k] - g_yz_0_y_yyyyz[k] * ab_y + g_yz_0_y_yyyyyz[k];

                g_yz_0_yy_yyyzz[k] = -g_z_0_y_yyyzz[k] - g_yz_0_y_yyyzz[k] * ab_y + g_yz_0_y_yyyyzz[k];

                g_yz_0_yy_yyzzz[k] = -g_z_0_y_yyzzz[k] - g_yz_0_y_yyzzz[k] * ab_y + g_yz_0_y_yyyzzz[k];

                g_yz_0_yy_yzzzz[k] = -g_z_0_y_yzzzz[k] - g_yz_0_y_yzzzz[k] * ab_y + g_yz_0_y_yyzzzz[k];

                g_yz_0_yy_zzzzz[k] = -g_z_0_y_zzzzz[k] - g_yz_0_y_zzzzz[k] * ab_y + g_yz_0_y_yzzzzz[k];
            }

            /// Set up 588-609 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 588 * ccomps * dcomps);

            auto g_yz_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 589 * ccomps * dcomps);

            auto g_yz_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 590 * ccomps * dcomps);

            auto g_yz_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 591 * ccomps * dcomps);

            auto g_yz_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 592 * ccomps * dcomps);

            auto g_yz_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 593 * ccomps * dcomps);

            auto g_yz_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 594 * ccomps * dcomps);

            auto g_yz_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 595 * ccomps * dcomps);

            auto g_yz_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 596 * ccomps * dcomps);

            auto g_yz_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 597 * ccomps * dcomps);

            auto g_yz_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 598 * ccomps * dcomps);

            auto g_yz_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 599 * ccomps * dcomps);

            auto g_yz_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 600 * ccomps * dcomps);

            auto g_yz_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 601 * ccomps * dcomps);

            auto g_yz_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 602 * ccomps * dcomps);

            auto g_yz_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 603 * ccomps * dcomps);

            auto g_yz_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 604 * ccomps * dcomps);

            auto g_yz_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 605 * ccomps * dcomps);

            auto g_yz_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 606 * ccomps * dcomps);

            auto g_yz_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 607 * ccomps * dcomps);

            auto g_yz_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 608 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yz_xxxxx, g_yz_0_yz_xxxxy, g_yz_0_yz_xxxxz, g_yz_0_yz_xxxyy, g_yz_0_yz_xxxyz, g_yz_0_yz_xxxzz, g_yz_0_yz_xxyyy, g_yz_0_yz_xxyyz, g_yz_0_yz_xxyzz, g_yz_0_yz_xxzzz, g_yz_0_yz_xyyyy, g_yz_0_yz_xyyyz, g_yz_0_yz_xyyzz, g_yz_0_yz_xyzzz, g_yz_0_yz_xzzzz, g_yz_0_yz_yyyyy, g_yz_0_yz_yyyyz, g_yz_0_yz_yyyzz, g_yz_0_yz_yyzzz, g_yz_0_yz_yzzzz, g_yz_0_yz_zzzzz, g_yz_0_z_xxxxx, g_yz_0_z_xxxxxy, g_yz_0_z_xxxxy, g_yz_0_z_xxxxyy, g_yz_0_z_xxxxyz, g_yz_0_z_xxxxz, g_yz_0_z_xxxyy, g_yz_0_z_xxxyyy, g_yz_0_z_xxxyyz, g_yz_0_z_xxxyz, g_yz_0_z_xxxyzz, g_yz_0_z_xxxzz, g_yz_0_z_xxyyy, g_yz_0_z_xxyyyy, g_yz_0_z_xxyyyz, g_yz_0_z_xxyyz, g_yz_0_z_xxyyzz, g_yz_0_z_xxyzz, g_yz_0_z_xxyzzz, g_yz_0_z_xxzzz, g_yz_0_z_xyyyy, g_yz_0_z_xyyyyy, g_yz_0_z_xyyyyz, g_yz_0_z_xyyyz, g_yz_0_z_xyyyzz, g_yz_0_z_xyyzz, g_yz_0_z_xyyzzz, g_yz_0_z_xyzzz, g_yz_0_z_xyzzzz, g_yz_0_z_xzzzz, g_yz_0_z_yyyyy, g_yz_0_z_yyyyyy, g_yz_0_z_yyyyyz, g_yz_0_z_yyyyz, g_yz_0_z_yyyyzz, g_yz_0_z_yyyzz, g_yz_0_z_yyyzzz, g_yz_0_z_yyzzz, g_yz_0_z_yyzzzz, g_yz_0_z_yzzzz, g_yz_0_z_yzzzzz, g_yz_0_z_zzzzz, g_z_0_z_xxxxx, g_z_0_z_xxxxy, g_z_0_z_xxxxz, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxzz, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyzz, g_z_0_z_xxzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyzz, g_z_0_z_xyzzz, g_z_0_z_xzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyzz, g_z_0_z_yyzzz, g_z_0_z_yzzzz, g_z_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yz_xxxxx[k] = -g_z_0_z_xxxxx[k] - g_yz_0_z_xxxxx[k] * ab_y + g_yz_0_z_xxxxxy[k];

                g_yz_0_yz_xxxxy[k] = -g_z_0_z_xxxxy[k] - g_yz_0_z_xxxxy[k] * ab_y + g_yz_0_z_xxxxyy[k];

                g_yz_0_yz_xxxxz[k] = -g_z_0_z_xxxxz[k] - g_yz_0_z_xxxxz[k] * ab_y + g_yz_0_z_xxxxyz[k];

                g_yz_0_yz_xxxyy[k] = -g_z_0_z_xxxyy[k] - g_yz_0_z_xxxyy[k] * ab_y + g_yz_0_z_xxxyyy[k];

                g_yz_0_yz_xxxyz[k] = -g_z_0_z_xxxyz[k] - g_yz_0_z_xxxyz[k] * ab_y + g_yz_0_z_xxxyyz[k];

                g_yz_0_yz_xxxzz[k] = -g_z_0_z_xxxzz[k] - g_yz_0_z_xxxzz[k] * ab_y + g_yz_0_z_xxxyzz[k];

                g_yz_0_yz_xxyyy[k] = -g_z_0_z_xxyyy[k] - g_yz_0_z_xxyyy[k] * ab_y + g_yz_0_z_xxyyyy[k];

                g_yz_0_yz_xxyyz[k] = -g_z_0_z_xxyyz[k] - g_yz_0_z_xxyyz[k] * ab_y + g_yz_0_z_xxyyyz[k];

                g_yz_0_yz_xxyzz[k] = -g_z_0_z_xxyzz[k] - g_yz_0_z_xxyzz[k] * ab_y + g_yz_0_z_xxyyzz[k];

                g_yz_0_yz_xxzzz[k] = -g_z_0_z_xxzzz[k] - g_yz_0_z_xxzzz[k] * ab_y + g_yz_0_z_xxyzzz[k];

                g_yz_0_yz_xyyyy[k] = -g_z_0_z_xyyyy[k] - g_yz_0_z_xyyyy[k] * ab_y + g_yz_0_z_xyyyyy[k];

                g_yz_0_yz_xyyyz[k] = -g_z_0_z_xyyyz[k] - g_yz_0_z_xyyyz[k] * ab_y + g_yz_0_z_xyyyyz[k];

                g_yz_0_yz_xyyzz[k] = -g_z_0_z_xyyzz[k] - g_yz_0_z_xyyzz[k] * ab_y + g_yz_0_z_xyyyzz[k];

                g_yz_0_yz_xyzzz[k] = -g_z_0_z_xyzzz[k] - g_yz_0_z_xyzzz[k] * ab_y + g_yz_0_z_xyyzzz[k];

                g_yz_0_yz_xzzzz[k] = -g_z_0_z_xzzzz[k] - g_yz_0_z_xzzzz[k] * ab_y + g_yz_0_z_xyzzzz[k];

                g_yz_0_yz_yyyyy[k] = -g_z_0_z_yyyyy[k] - g_yz_0_z_yyyyy[k] * ab_y + g_yz_0_z_yyyyyy[k];

                g_yz_0_yz_yyyyz[k] = -g_z_0_z_yyyyz[k] - g_yz_0_z_yyyyz[k] * ab_y + g_yz_0_z_yyyyyz[k];

                g_yz_0_yz_yyyzz[k] = -g_z_0_z_yyyzz[k] - g_yz_0_z_yyyzz[k] * ab_y + g_yz_0_z_yyyyzz[k];

                g_yz_0_yz_yyzzz[k] = -g_z_0_z_yyzzz[k] - g_yz_0_z_yyzzz[k] * ab_y + g_yz_0_z_yyyzzz[k];

                g_yz_0_yz_yzzzz[k] = -g_z_0_z_yzzzz[k] - g_yz_0_z_yzzzz[k] * ab_y + g_yz_0_z_yyzzzz[k];

                g_yz_0_yz_zzzzz[k] = -g_z_0_z_zzzzz[k] - g_yz_0_z_zzzzz[k] * ab_y + g_yz_0_z_yzzzzz[k];
            }

            /// Set up 609-630 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 609 * ccomps * dcomps);

            auto g_yz_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 610 * ccomps * dcomps);

            auto g_yz_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 611 * ccomps * dcomps);

            auto g_yz_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 612 * ccomps * dcomps);

            auto g_yz_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 613 * ccomps * dcomps);

            auto g_yz_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 614 * ccomps * dcomps);

            auto g_yz_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 615 * ccomps * dcomps);

            auto g_yz_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 616 * ccomps * dcomps);

            auto g_yz_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 617 * ccomps * dcomps);

            auto g_yz_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 618 * ccomps * dcomps);

            auto g_yz_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 619 * ccomps * dcomps);

            auto g_yz_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 620 * ccomps * dcomps);

            auto g_yz_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 621 * ccomps * dcomps);

            auto g_yz_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 622 * ccomps * dcomps);

            auto g_yz_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 623 * ccomps * dcomps);

            auto g_yz_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 624 * ccomps * dcomps);

            auto g_yz_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 625 * ccomps * dcomps);

            auto g_yz_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 626 * ccomps * dcomps);

            auto g_yz_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 627 * ccomps * dcomps);

            auto g_yz_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 628 * ccomps * dcomps);

            auto g_yz_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxxxx, g_y_0_z_xxxxy, g_y_0_z_xxxxz, g_y_0_z_xxxyy, g_y_0_z_xxxyz, g_y_0_z_xxxzz, g_y_0_z_xxyyy, g_y_0_z_xxyyz, g_y_0_z_xxyzz, g_y_0_z_xxzzz, g_y_0_z_xyyyy, g_y_0_z_xyyyz, g_y_0_z_xyyzz, g_y_0_z_xyzzz, g_y_0_z_xzzzz, g_y_0_z_yyyyy, g_y_0_z_yyyyz, g_y_0_z_yyyzz, g_y_0_z_yyzzz, g_y_0_z_yzzzz, g_y_0_z_zzzzz, g_yz_0_z_xxxxx, g_yz_0_z_xxxxxz, g_yz_0_z_xxxxy, g_yz_0_z_xxxxyz, g_yz_0_z_xxxxz, g_yz_0_z_xxxxzz, g_yz_0_z_xxxyy, g_yz_0_z_xxxyyz, g_yz_0_z_xxxyz, g_yz_0_z_xxxyzz, g_yz_0_z_xxxzz, g_yz_0_z_xxxzzz, g_yz_0_z_xxyyy, g_yz_0_z_xxyyyz, g_yz_0_z_xxyyz, g_yz_0_z_xxyyzz, g_yz_0_z_xxyzz, g_yz_0_z_xxyzzz, g_yz_0_z_xxzzz, g_yz_0_z_xxzzzz, g_yz_0_z_xyyyy, g_yz_0_z_xyyyyz, g_yz_0_z_xyyyz, g_yz_0_z_xyyyzz, g_yz_0_z_xyyzz, g_yz_0_z_xyyzzz, g_yz_0_z_xyzzz, g_yz_0_z_xyzzzz, g_yz_0_z_xzzzz, g_yz_0_z_xzzzzz, g_yz_0_z_yyyyy, g_yz_0_z_yyyyyz, g_yz_0_z_yyyyz, g_yz_0_z_yyyyzz, g_yz_0_z_yyyzz, g_yz_0_z_yyyzzz, g_yz_0_z_yyzzz, g_yz_0_z_yyzzzz, g_yz_0_z_yzzzz, g_yz_0_z_yzzzzz, g_yz_0_z_zzzzz, g_yz_0_z_zzzzzz, g_yz_0_zz_xxxxx, g_yz_0_zz_xxxxy, g_yz_0_zz_xxxxz, g_yz_0_zz_xxxyy, g_yz_0_zz_xxxyz, g_yz_0_zz_xxxzz, g_yz_0_zz_xxyyy, g_yz_0_zz_xxyyz, g_yz_0_zz_xxyzz, g_yz_0_zz_xxzzz, g_yz_0_zz_xyyyy, g_yz_0_zz_xyyyz, g_yz_0_zz_xyyzz, g_yz_0_zz_xyzzz, g_yz_0_zz_xzzzz, g_yz_0_zz_yyyyy, g_yz_0_zz_yyyyz, g_yz_0_zz_yyyzz, g_yz_0_zz_yyzzz, g_yz_0_zz_yzzzz, g_yz_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zz_xxxxx[k] = -g_y_0_z_xxxxx[k] - g_yz_0_z_xxxxx[k] * ab_z + g_yz_0_z_xxxxxz[k];

                g_yz_0_zz_xxxxy[k] = -g_y_0_z_xxxxy[k] - g_yz_0_z_xxxxy[k] * ab_z + g_yz_0_z_xxxxyz[k];

                g_yz_0_zz_xxxxz[k] = -g_y_0_z_xxxxz[k] - g_yz_0_z_xxxxz[k] * ab_z + g_yz_0_z_xxxxzz[k];

                g_yz_0_zz_xxxyy[k] = -g_y_0_z_xxxyy[k] - g_yz_0_z_xxxyy[k] * ab_z + g_yz_0_z_xxxyyz[k];

                g_yz_0_zz_xxxyz[k] = -g_y_0_z_xxxyz[k] - g_yz_0_z_xxxyz[k] * ab_z + g_yz_0_z_xxxyzz[k];

                g_yz_0_zz_xxxzz[k] = -g_y_0_z_xxxzz[k] - g_yz_0_z_xxxzz[k] * ab_z + g_yz_0_z_xxxzzz[k];

                g_yz_0_zz_xxyyy[k] = -g_y_0_z_xxyyy[k] - g_yz_0_z_xxyyy[k] * ab_z + g_yz_0_z_xxyyyz[k];

                g_yz_0_zz_xxyyz[k] = -g_y_0_z_xxyyz[k] - g_yz_0_z_xxyyz[k] * ab_z + g_yz_0_z_xxyyzz[k];

                g_yz_0_zz_xxyzz[k] = -g_y_0_z_xxyzz[k] - g_yz_0_z_xxyzz[k] * ab_z + g_yz_0_z_xxyzzz[k];

                g_yz_0_zz_xxzzz[k] = -g_y_0_z_xxzzz[k] - g_yz_0_z_xxzzz[k] * ab_z + g_yz_0_z_xxzzzz[k];

                g_yz_0_zz_xyyyy[k] = -g_y_0_z_xyyyy[k] - g_yz_0_z_xyyyy[k] * ab_z + g_yz_0_z_xyyyyz[k];

                g_yz_0_zz_xyyyz[k] = -g_y_0_z_xyyyz[k] - g_yz_0_z_xyyyz[k] * ab_z + g_yz_0_z_xyyyzz[k];

                g_yz_0_zz_xyyzz[k] = -g_y_0_z_xyyzz[k] - g_yz_0_z_xyyzz[k] * ab_z + g_yz_0_z_xyyzzz[k];

                g_yz_0_zz_xyzzz[k] = -g_y_0_z_xyzzz[k] - g_yz_0_z_xyzzz[k] * ab_z + g_yz_0_z_xyzzzz[k];

                g_yz_0_zz_xzzzz[k] = -g_y_0_z_xzzzz[k] - g_yz_0_z_xzzzz[k] * ab_z + g_yz_0_z_xzzzzz[k];

                g_yz_0_zz_yyyyy[k] = -g_y_0_z_yyyyy[k] - g_yz_0_z_yyyyy[k] * ab_z + g_yz_0_z_yyyyyz[k];

                g_yz_0_zz_yyyyz[k] = -g_y_0_z_yyyyz[k] - g_yz_0_z_yyyyz[k] * ab_z + g_yz_0_z_yyyyzz[k];

                g_yz_0_zz_yyyzz[k] = -g_y_0_z_yyyzz[k] - g_yz_0_z_yyyzz[k] * ab_z + g_yz_0_z_yyyzzz[k];

                g_yz_0_zz_yyzzz[k] = -g_y_0_z_yyzzz[k] - g_yz_0_z_yyzzz[k] * ab_z + g_yz_0_z_yyzzzz[k];

                g_yz_0_zz_yzzzz[k] = -g_y_0_z_yzzzz[k] - g_yz_0_z_yzzzz[k] * ab_z + g_yz_0_z_yzzzzz[k];

                g_yz_0_zz_zzzzz[k] = -g_y_0_z_zzzzz[k] - g_yz_0_z_zzzzz[k] * ab_z + g_yz_0_z_zzzzzz[k];
            }

            /// Set up 630-651 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xx_xxxxx = cbuffer.data(dh_geom_20_off + 630 * ccomps * dcomps);

            auto g_zz_0_xx_xxxxy = cbuffer.data(dh_geom_20_off + 631 * ccomps * dcomps);

            auto g_zz_0_xx_xxxxz = cbuffer.data(dh_geom_20_off + 632 * ccomps * dcomps);

            auto g_zz_0_xx_xxxyy = cbuffer.data(dh_geom_20_off + 633 * ccomps * dcomps);

            auto g_zz_0_xx_xxxyz = cbuffer.data(dh_geom_20_off + 634 * ccomps * dcomps);

            auto g_zz_0_xx_xxxzz = cbuffer.data(dh_geom_20_off + 635 * ccomps * dcomps);

            auto g_zz_0_xx_xxyyy = cbuffer.data(dh_geom_20_off + 636 * ccomps * dcomps);

            auto g_zz_0_xx_xxyyz = cbuffer.data(dh_geom_20_off + 637 * ccomps * dcomps);

            auto g_zz_0_xx_xxyzz = cbuffer.data(dh_geom_20_off + 638 * ccomps * dcomps);

            auto g_zz_0_xx_xxzzz = cbuffer.data(dh_geom_20_off + 639 * ccomps * dcomps);

            auto g_zz_0_xx_xyyyy = cbuffer.data(dh_geom_20_off + 640 * ccomps * dcomps);

            auto g_zz_0_xx_xyyyz = cbuffer.data(dh_geom_20_off + 641 * ccomps * dcomps);

            auto g_zz_0_xx_xyyzz = cbuffer.data(dh_geom_20_off + 642 * ccomps * dcomps);

            auto g_zz_0_xx_xyzzz = cbuffer.data(dh_geom_20_off + 643 * ccomps * dcomps);

            auto g_zz_0_xx_xzzzz = cbuffer.data(dh_geom_20_off + 644 * ccomps * dcomps);

            auto g_zz_0_xx_yyyyy = cbuffer.data(dh_geom_20_off + 645 * ccomps * dcomps);

            auto g_zz_0_xx_yyyyz = cbuffer.data(dh_geom_20_off + 646 * ccomps * dcomps);

            auto g_zz_0_xx_yyyzz = cbuffer.data(dh_geom_20_off + 647 * ccomps * dcomps);

            auto g_zz_0_xx_yyzzz = cbuffer.data(dh_geom_20_off + 648 * ccomps * dcomps);

            auto g_zz_0_xx_yzzzz = cbuffer.data(dh_geom_20_off + 649 * ccomps * dcomps);

            auto g_zz_0_xx_zzzzz = cbuffer.data(dh_geom_20_off + 650 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_x_xxxxx, g_zz_0_x_xxxxxx, g_zz_0_x_xxxxxy, g_zz_0_x_xxxxxz, g_zz_0_x_xxxxy, g_zz_0_x_xxxxyy, g_zz_0_x_xxxxyz, g_zz_0_x_xxxxz, g_zz_0_x_xxxxzz, g_zz_0_x_xxxyy, g_zz_0_x_xxxyyy, g_zz_0_x_xxxyyz, g_zz_0_x_xxxyz, g_zz_0_x_xxxyzz, g_zz_0_x_xxxzz, g_zz_0_x_xxxzzz, g_zz_0_x_xxyyy, g_zz_0_x_xxyyyy, g_zz_0_x_xxyyyz, g_zz_0_x_xxyyz, g_zz_0_x_xxyyzz, g_zz_0_x_xxyzz, g_zz_0_x_xxyzzz, g_zz_0_x_xxzzz, g_zz_0_x_xxzzzz, g_zz_0_x_xyyyy, g_zz_0_x_xyyyyy, g_zz_0_x_xyyyyz, g_zz_0_x_xyyyz, g_zz_0_x_xyyyzz, g_zz_0_x_xyyzz, g_zz_0_x_xyyzzz, g_zz_0_x_xyzzz, g_zz_0_x_xyzzzz, g_zz_0_x_xzzzz, g_zz_0_x_xzzzzz, g_zz_0_x_yyyyy, g_zz_0_x_yyyyz, g_zz_0_x_yyyzz, g_zz_0_x_yyzzz, g_zz_0_x_yzzzz, g_zz_0_x_zzzzz, g_zz_0_xx_xxxxx, g_zz_0_xx_xxxxy, g_zz_0_xx_xxxxz, g_zz_0_xx_xxxyy, g_zz_0_xx_xxxyz, g_zz_0_xx_xxxzz, g_zz_0_xx_xxyyy, g_zz_0_xx_xxyyz, g_zz_0_xx_xxyzz, g_zz_0_xx_xxzzz, g_zz_0_xx_xyyyy, g_zz_0_xx_xyyyz, g_zz_0_xx_xyyzz, g_zz_0_xx_xyzzz, g_zz_0_xx_xzzzz, g_zz_0_xx_yyyyy, g_zz_0_xx_yyyyz, g_zz_0_xx_yyyzz, g_zz_0_xx_yyzzz, g_zz_0_xx_yzzzz, g_zz_0_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xx_xxxxx[k] = -g_zz_0_x_xxxxx[k] * ab_x + g_zz_0_x_xxxxxx[k];

                g_zz_0_xx_xxxxy[k] = -g_zz_0_x_xxxxy[k] * ab_x + g_zz_0_x_xxxxxy[k];

                g_zz_0_xx_xxxxz[k] = -g_zz_0_x_xxxxz[k] * ab_x + g_zz_0_x_xxxxxz[k];

                g_zz_0_xx_xxxyy[k] = -g_zz_0_x_xxxyy[k] * ab_x + g_zz_0_x_xxxxyy[k];

                g_zz_0_xx_xxxyz[k] = -g_zz_0_x_xxxyz[k] * ab_x + g_zz_0_x_xxxxyz[k];

                g_zz_0_xx_xxxzz[k] = -g_zz_0_x_xxxzz[k] * ab_x + g_zz_0_x_xxxxzz[k];

                g_zz_0_xx_xxyyy[k] = -g_zz_0_x_xxyyy[k] * ab_x + g_zz_0_x_xxxyyy[k];

                g_zz_0_xx_xxyyz[k] = -g_zz_0_x_xxyyz[k] * ab_x + g_zz_0_x_xxxyyz[k];

                g_zz_0_xx_xxyzz[k] = -g_zz_0_x_xxyzz[k] * ab_x + g_zz_0_x_xxxyzz[k];

                g_zz_0_xx_xxzzz[k] = -g_zz_0_x_xxzzz[k] * ab_x + g_zz_0_x_xxxzzz[k];

                g_zz_0_xx_xyyyy[k] = -g_zz_0_x_xyyyy[k] * ab_x + g_zz_0_x_xxyyyy[k];

                g_zz_0_xx_xyyyz[k] = -g_zz_0_x_xyyyz[k] * ab_x + g_zz_0_x_xxyyyz[k];

                g_zz_0_xx_xyyzz[k] = -g_zz_0_x_xyyzz[k] * ab_x + g_zz_0_x_xxyyzz[k];

                g_zz_0_xx_xyzzz[k] = -g_zz_0_x_xyzzz[k] * ab_x + g_zz_0_x_xxyzzz[k];

                g_zz_0_xx_xzzzz[k] = -g_zz_0_x_xzzzz[k] * ab_x + g_zz_0_x_xxzzzz[k];

                g_zz_0_xx_yyyyy[k] = -g_zz_0_x_yyyyy[k] * ab_x + g_zz_0_x_xyyyyy[k];

                g_zz_0_xx_yyyyz[k] = -g_zz_0_x_yyyyz[k] * ab_x + g_zz_0_x_xyyyyz[k];

                g_zz_0_xx_yyyzz[k] = -g_zz_0_x_yyyzz[k] * ab_x + g_zz_0_x_xyyyzz[k];

                g_zz_0_xx_yyzzz[k] = -g_zz_0_x_yyzzz[k] * ab_x + g_zz_0_x_xyyzzz[k];

                g_zz_0_xx_yzzzz[k] = -g_zz_0_x_yzzzz[k] * ab_x + g_zz_0_x_xyzzzz[k];

                g_zz_0_xx_zzzzz[k] = -g_zz_0_x_zzzzz[k] * ab_x + g_zz_0_x_xzzzzz[k];
            }

            /// Set up 651-672 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xy_xxxxx = cbuffer.data(dh_geom_20_off + 651 * ccomps * dcomps);

            auto g_zz_0_xy_xxxxy = cbuffer.data(dh_geom_20_off + 652 * ccomps * dcomps);

            auto g_zz_0_xy_xxxxz = cbuffer.data(dh_geom_20_off + 653 * ccomps * dcomps);

            auto g_zz_0_xy_xxxyy = cbuffer.data(dh_geom_20_off + 654 * ccomps * dcomps);

            auto g_zz_0_xy_xxxyz = cbuffer.data(dh_geom_20_off + 655 * ccomps * dcomps);

            auto g_zz_0_xy_xxxzz = cbuffer.data(dh_geom_20_off + 656 * ccomps * dcomps);

            auto g_zz_0_xy_xxyyy = cbuffer.data(dh_geom_20_off + 657 * ccomps * dcomps);

            auto g_zz_0_xy_xxyyz = cbuffer.data(dh_geom_20_off + 658 * ccomps * dcomps);

            auto g_zz_0_xy_xxyzz = cbuffer.data(dh_geom_20_off + 659 * ccomps * dcomps);

            auto g_zz_0_xy_xxzzz = cbuffer.data(dh_geom_20_off + 660 * ccomps * dcomps);

            auto g_zz_0_xy_xyyyy = cbuffer.data(dh_geom_20_off + 661 * ccomps * dcomps);

            auto g_zz_0_xy_xyyyz = cbuffer.data(dh_geom_20_off + 662 * ccomps * dcomps);

            auto g_zz_0_xy_xyyzz = cbuffer.data(dh_geom_20_off + 663 * ccomps * dcomps);

            auto g_zz_0_xy_xyzzz = cbuffer.data(dh_geom_20_off + 664 * ccomps * dcomps);

            auto g_zz_0_xy_xzzzz = cbuffer.data(dh_geom_20_off + 665 * ccomps * dcomps);

            auto g_zz_0_xy_yyyyy = cbuffer.data(dh_geom_20_off + 666 * ccomps * dcomps);

            auto g_zz_0_xy_yyyyz = cbuffer.data(dh_geom_20_off + 667 * ccomps * dcomps);

            auto g_zz_0_xy_yyyzz = cbuffer.data(dh_geom_20_off + 668 * ccomps * dcomps);

            auto g_zz_0_xy_yyzzz = cbuffer.data(dh_geom_20_off + 669 * ccomps * dcomps);

            auto g_zz_0_xy_yzzzz = cbuffer.data(dh_geom_20_off + 670 * ccomps * dcomps);

            auto g_zz_0_xy_zzzzz = cbuffer.data(dh_geom_20_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xy_xxxxx, g_zz_0_xy_xxxxy, g_zz_0_xy_xxxxz, g_zz_0_xy_xxxyy, g_zz_0_xy_xxxyz, g_zz_0_xy_xxxzz, g_zz_0_xy_xxyyy, g_zz_0_xy_xxyyz, g_zz_0_xy_xxyzz, g_zz_0_xy_xxzzz, g_zz_0_xy_xyyyy, g_zz_0_xy_xyyyz, g_zz_0_xy_xyyzz, g_zz_0_xy_xyzzz, g_zz_0_xy_xzzzz, g_zz_0_xy_yyyyy, g_zz_0_xy_yyyyz, g_zz_0_xy_yyyzz, g_zz_0_xy_yyzzz, g_zz_0_xy_yzzzz, g_zz_0_xy_zzzzz, g_zz_0_y_xxxxx, g_zz_0_y_xxxxxx, g_zz_0_y_xxxxxy, g_zz_0_y_xxxxxz, g_zz_0_y_xxxxy, g_zz_0_y_xxxxyy, g_zz_0_y_xxxxyz, g_zz_0_y_xxxxz, g_zz_0_y_xxxxzz, g_zz_0_y_xxxyy, g_zz_0_y_xxxyyy, g_zz_0_y_xxxyyz, g_zz_0_y_xxxyz, g_zz_0_y_xxxyzz, g_zz_0_y_xxxzz, g_zz_0_y_xxxzzz, g_zz_0_y_xxyyy, g_zz_0_y_xxyyyy, g_zz_0_y_xxyyyz, g_zz_0_y_xxyyz, g_zz_0_y_xxyyzz, g_zz_0_y_xxyzz, g_zz_0_y_xxyzzz, g_zz_0_y_xxzzz, g_zz_0_y_xxzzzz, g_zz_0_y_xyyyy, g_zz_0_y_xyyyyy, g_zz_0_y_xyyyyz, g_zz_0_y_xyyyz, g_zz_0_y_xyyyzz, g_zz_0_y_xyyzz, g_zz_0_y_xyyzzz, g_zz_0_y_xyzzz, g_zz_0_y_xyzzzz, g_zz_0_y_xzzzz, g_zz_0_y_xzzzzz, g_zz_0_y_yyyyy, g_zz_0_y_yyyyz, g_zz_0_y_yyyzz, g_zz_0_y_yyzzz, g_zz_0_y_yzzzz, g_zz_0_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xy_xxxxx[k] = -g_zz_0_y_xxxxx[k] * ab_x + g_zz_0_y_xxxxxx[k];

                g_zz_0_xy_xxxxy[k] = -g_zz_0_y_xxxxy[k] * ab_x + g_zz_0_y_xxxxxy[k];

                g_zz_0_xy_xxxxz[k] = -g_zz_0_y_xxxxz[k] * ab_x + g_zz_0_y_xxxxxz[k];

                g_zz_0_xy_xxxyy[k] = -g_zz_0_y_xxxyy[k] * ab_x + g_zz_0_y_xxxxyy[k];

                g_zz_0_xy_xxxyz[k] = -g_zz_0_y_xxxyz[k] * ab_x + g_zz_0_y_xxxxyz[k];

                g_zz_0_xy_xxxzz[k] = -g_zz_0_y_xxxzz[k] * ab_x + g_zz_0_y_xxxxzz[k];

                g_zz_0_xy_xxyyy[k] = -g_zz_0_y_xxyyy[k] * ab_x + g_zz_0_y_xxxyyy[k];

                g_zz_0_xy_xxyyz[k] = -g_zz_0_y_xxyyz[k] * ab_x + g_zz_0_y_xxxyyz[k];

                g_zz_0_xy_xxyzz[k] = -g_zz_0_y_xxyzz[k] * ab_x + g_zz_0_y_xxxyzz[k];

                g_zz_0_xy_xxzzz[k] = -g_zz_0_y_xxzzz[k] * ab_x + g_zz_0_y_xxxzzz[k];

                g_zz_0_xy_xyyyy[k] = -g_zz_0_y_xyyyy[k] * ab_x + g_zz_0_y_xxyyyy[k];

                g_zz_0_xy_xyyyz[k] = -g_zz_0_y_xyyyz[k] * ab_x + g_zz_0_y_xxyyyz[k];

                g_zz_0_xy_xyyzz[k] = -g_zz_0_y_xyyzz[k] * ab_x + g_zz_0_y_xxyyzz[k];

                g_zz_0_xy_xyzzz[k] = -g_zz_0_y_xyzzz[k] * ab_x + g_zz_0_y_xxyzzz[k];

                g_zz_0_xy_xzzzz[k] = -g_zz_0_y_xzzzz[k] * ab_x + g_zz_0_y_xxzzzz[k];

                g_zz_0_xy_yyyyy[k] = -g_zz_0_y_yyyyy[k] * ab_x + g_zz_0_y_xyyyyy[k];

                g_zz_0_xy_yyyyz[k] = -g_zz_0_y_yyyyz[k] * ab_x + g_zz_0_y_xyyyyz[k];

                g_zz_0_xy_yyyzz[k] = -g_zz_0_y_yyyzz[k] * ab_x + g_zz_0_y_xyyyzz[k];

                g_zz_0_xy_yyzzz[k] = -g_zz_0_y_yyzzz[k] * ab_x + g_zz_0_y_xyyzzz[k];

                g_zz_0_xy_yzzzz[k] = -g_zz_0_y_yzzzz[k] * ab_x + g_zz_0_y_xyzzzz[k];

                g_zz_0_xy_zzzzz[k] = -g_zz_0_y_zzzzz[k] * ab_x + g_zz_0_y_xzzzzz[k];
            }

            /// Set up 672-693 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xz_xxxxx = cbuffer.data(dh_geom_20_off + 672 * ccomps * dcomps);

            auto g_zz_0_xz_xxxxy = cbuffer.data(dh_geom_20_off + 673 * ccomps * dcomps);

            auto g_zz_0_xz_xxxxz = cbuffer.data(dh_geom_20_off + 674 * ccomps * dcomps);

            auto g_zz_0_xz_xxxyy = cbuffer.data(dh_geom_20_off + 675 * ccomps * dcomps);

            auto g_zz_0_xz_xxxyz = cbuffer.data(dh_geom_20_off + 676 * ccomps * dcomps);

            auto g_zz_0_xz_xxxzz = cbuffer.data(dh_geom_20_off + 677 * ccomps * dcomps);

            auto g_zz_0_xz_xxyyy = cbuffer.data(dh_geom_20_off + 678 * ccomps * dcomps);

            auto g_zz_0_xz_xxyyz = cbuffer.data(dh_geom_20_off + 679 * ccomps * dcomps);

            auto g_zz_0_xz_xxyzz = cbuffer.data(dh_geom_20_off + 680 * ccomps * dcomps);

            auto g_zz_0_xz_xxzzz = cbuffer.data(dh_geom_20_off + 681 * ccomps * dcomps);

            auto g_zz_0_xz_xyyyy = cbuffer.data(dh_geom_20_off + 682 * ccomps * dcomps);

            auto g_zz_0_xz_xyyyz = cbuffer.data(dh_geom_20_off + 683 * ccomps * dcomps);

            auto g_zz_0_xz_xyyzz = cbuffer.data(dh_geom_20_off + 684 * ccomps * dcomps);

            auto g_zz_0_xz_xyzzz = cbuffer.data(dh_geom_20_off + 685 * ccomps * dcomps);

            auto g_zz_0_xz_xzzzz = cbuffer.data(dh_geom_20_off + 686 * ccomps * dcomps);

            auto g_zz_0_xz_yyyyy = cbuffer.data(dh_geom_20_off + 687 * ccomps * dcomps);

            auto g_zz_0_xz_yyyyz = cbuffer.data(dh_geom_20_off + 688 * ccomps * dcomps);

            auto g_zz_0_xz_yyyzz = cbuffer.data(dh_geom_20_off + 689 * ccomps * dcomps);

            auto g_zz_0_xz_yyzzz = cbuffer.data(dh_geom_20_off + 690 * ccomps * dcomps);

            auto g_zz_0_xz_yzzzz = cbuffer.data(dh_geom_20_off + 691 * ccomps * dcomps);

            auto g_zz_0_xz_zzzzz = cbuffer.data(dh_geom_20_off + 692 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xz_xxxxx, g_zz_0_xz_xxxxy, g_zz_0_xz_xxxxz, g_zz_0_xz_xxxyy, g_zz_0_xz_xxxyz, g_zz_0_xz_xxxzz, g_zz_0_xz_xxyyy, g_zz_0_xz_xxyyz, g_zz_0_xz_xxyzz, g_zz_0_xz_xxzzz, g_zz_0_xz_xyyyy, g_zz_0_xz_xyyyz, g_zz_0_xz_xyyzz, g_zz_0_xz_xyzzz, g_zz_0_xz_xzzzz, g_zz_0_xz_yyyyy, g_zz_0_xz_yyyyz, g_zz_0_xz_yyyzz, g_zz_0_xz_yyzzz, g_zz_0_xz_yzzzz, g_zz_0_xz_zzzzz, g_zz_0_z_xxxxx, g_zz_0_z_xxxxxx, g_zz_0_z_xxxxxy, g_zz_0_z_xxxxxz, g_zz_0_z_xxxxy, g_zz_0_z_xxxxyy, g_zz_0_z_xxxxyz, g_zz_0_z_xxxxz, g_zz_0_z_xxxxzz, g_zz_0_z_xxxyy, g_zz_0_z_xxxyyy, g_zz_0_z_xxxyyz, g_zz_0_z_xxxyz, g_zz_0_z_xxxyzz, g_zz_0_z_xxxzz, g_zz_0_z_xxxzzz, g_zz_0_z_xxyyy, g_zz_0_z_xxyyyy, g_zz_0_z_xxyyyz, g_zz_0_z_xxyyz, g_zz_0_z_xxyyzz, g_zz_0_z_xxyzz, g_zz_0_z_xxyzzz, g_zz_0_z_xxzzz, g_zz_0_z_xxzzzz, g_zz_0_z_xyyyy, g_zz_0_z_xyyyyy, g_zz_0_z_xyyyyz, g_zz_0_z_xyyyz, g_zz_0_z_xyyyzz, g_zz_0_z_xyyzz, g_zz_0_z_xyyzzz, g_zz_0_z_xyzzz, g_zz_0_z_xyzzzz, g_zz_0_z_xzzzz, g_zz_0_z_xzzzzz, g_zz_0_z_yyyyy, g_zz_0_z_yyyyz, g_zz_0_z_yyyzz, g_zz_0_z_yyzzz, g_zz_0_z_yzzzz, g_zz_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xz_xxxxx[k] = -g_zz_0_z_xxxxx[k] * ab_x + g_zz_0_z_xxxxxx[k];

                g_zz_0_xz_xxxxy[k] = -g_zz_0_z_xxxxy[k] * ab_x + g_zz_0_z_xxxxxy[k];

                g_zz_0_xz_xxxxz[k] = -g_zz_0_z_xxxxz[k] * ab_x + g_zz_0_z_xxxxxz[k];

                g_zz_0_xz_xxxyy[k] = -g_zz_0_z_xxxyy[k] * ab_x + g_zz_0_z_xxxxyy[k];

                g_zz_0_xz_xxxyz[k] = -g_zz_0_z_xxxyz[k] * ab_x + g_zz_0_z_xxxxyz[k];

                g_zz_0_xz_xxxzz[k] = -g_zz_0_z_xxxzz[k] * ab_x + g_zz_0_z_xxxxzz[k];

                g_zz_0_xz_xxyyy[k] = -g_zz_0_z_xxyyy[k] * ab_x + g_zz_0_z_xxxyyy[k];

                g_zz_0_xz_xxyyz[k] = -g_zz_0_z_xxyyz[k] * ab_x + g_zz_0_z_xxxyyz[k];

                g_zz_0_xz_xxyzz[k] = -g_zz_0_z_xxyzz[k] * ab_x + g_zz_0_z_xxxyzz[k];

                g_zz_0_xz_xxzzz[k] = -g_zz_0_z_xxzzz[k] * ab_x + g_zz_0_z_xxxzzz[k];

                g_zz_0_xz_xyyyy[k] = -g_zz_0_z_xyyyy[k] * ab_x + g_zz_0_z_xxyyyy[k];

                g_zz_0_xz_xyyyz[k] = -g_zz_0_z_xyyyz[k] * ab_x + g_zz_0_z_xxyyyz[k];

                g_zz_0_xz_xyyzz[k] = -g_zz_0_z_xyyzz[k] * ab_x + g_zz_0_z_xxyyzz[k];

                g_zz_0_xz_xyzzz[k] = -g_zz_0_z_xyzzz[k] * ab_x + g_zz_0_z_xxyzzz[k];

                g_zz_0_xz_xzzzz[k] = -g_zz_0_z_xzzzz[k] * ab_x + g_zz_0_z_xxzzzz[k];

                g_zz_0_xz_yyyyy[k] = -g_zz_0_z_yyyyy[k] * ab_x + g_zz_0_z_xyyyyy[k];

                g_zz_0_xz_yyyyz[k] = -g_zz_0_z_yyyyz[k] * ab_x + g_zz_0_z_xyyyyz[k];

                g_zz_0_xz_yyyzz[k] = -g_zz_0_z_yyyzz[k] * ab_x + g_zz_0_z_xyyyzz[k];

                g_zz_0_xz_yyzzz[k] = -g_zz_0_z_yyzzz[k] * ab_x + g_zz_0_z_xyyzzz[k];

                g_zz_0_xz_yzzzz[k] = -g_zz_0_z_yzzzz[k] * ab_x + g_zz_0_z_xyzzzz[k];

                g_zz_0_xz_zzzzz[k] = -g_zz_0_z_zzzzz[k] * ab_x + g_zz_0_z_xzzzzz[k];
            }

            /// Set up 693-714 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yy_xxxxx = cbuffer.data(dh_geom_20_off + 693 * ccomps * dcomps);

            auto g_zz_0_yy_xxxxy = cbuffer.data(dh_geom_20_off + 694 * ccomps * dcomps);

            auto g_zz_0_yy_xxxxz = cbuffer.data(dh_geom_20_off + 695 * ccomps * dcomps);

            auto g_zz_0_yy_xxxyy = cbuffer.data(dh_geom_20_off + 696 * ccomps * dcomps);

            auto g_zz_0_yy_xxxyz = cbuffer.data(dh_geom_20_off + 697 * ccomps * dcomps);

            auto g_zz_0_yy_xxxzz = cbuffer.data(dh_geom_20_off + 698 * ccomps * dcomps);

            auto g_zz_0_yy_xxyyy = cbuffer.data(dh_geom_20_off + 699 * ccomps * dcomps);

            auto g_zz_0_yy_xxyyz = cbuffer.data(dh_geom_20_off + 700 * ccomps * dcomps);

            auto g_zz_0_yy_xxyzz = cbuffer.data(dh_geom_20_off + 701 * ccomps * dcomps);

            auto g_zz_0_yy_xxzzz = cbuffer.data(dh_geom_20_off + 702 * ccomps * dcomps);

            auto g_zz_0_yy_xyyyy = cbuffer.data(dh_geom_20_off + 703 * ccomps * dcomps);

            auto g_zz_0_yy_xyyyz = cbuffer.data(dh_geom_20_off + 704 * ccomps * dcomps);

            auto g_zz_0_yy_xyyzz = cbuffer.data(dh_geom_20_off + 705 * ccomps * dcomps);

            auto g_zz_0_yy_xyzzz = cbuffer.data(dh_geom_20_off + 706 * ccomps * dcomps);

            auto g_zz_0_yy_xzzzz = cbuffer.data(dh_geom_20_off + 707 * ccomps * dcomps);

            auto g_zz_0_yy_yyyyy = cbuffer.data(dh_geom_20_off + 708 * ccomps * dcomps);

            auto g_zz_0_yy_yyyyz = cbuffer.data(dh_geom_20_off + 709 * ccomps * dcomps);

            auto g_zz_0_yy_yyyzz = cbuffer.data(dh_geom_20_off + 710 * ccomps * dcomps);

            auto g_zz_0_yy_yyzzz = cbuffer.data(dh_geom_20_off + 711 * ccomps * dcomps);

            auto g_zz_0_yy_yzzzz = cbuffer.data(dh_geom_20_off + 712 * ccomps * dcomps);

            auto g_zz_0_yy_zzzzz = cbuffer.data(dh_geom_20_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_y_xxxxx, g_zz_0_y_xxxxxy, g_zz_0_y_xxxxy, g_zz_0_y_xxxxyy, g_zz_0_y_xxxxyz, g_zz_0_y_xxxxz, g_zz_0_y_xxxyy, g_zz_0_y_xxxyyy, g_zz_0_y_xxxyyz, g_zz_0_y_xxxyz, g_zz_0_y_xxxyzz, g_zz_0_y_xxxzz, g_zz_0_y_xxyyy, g_zz_0_y_xxyyyy, g_zz_0_y_xxyyyz, g_zz_0_y_xxyyz, g_zz_0_y_xxyyzz, g_zz_0_y_xxyzz, g_zz_0_y_xxyzzz, g_zz_0_y_xxzzz, g_zz_0_y_xyyyy, g_zz_0_y_xyyyyy, g_zz_0_y_xyyyyz, g_zz_0_y_xyyyz, g_zz_0_y_xyyyzz, g_zz_0_y_xyyzz, g_zz_0_y_xyyzzz, g_zz_0_y_xyzzz, g_zz_0_y_xyzzzz, g_zz_0_y_xzzzz, g_zz_0_y_yyyyy, g_zz_0_y_yyyyyy, g_zz_0_y_yyyyyz, g_zz_0_y_yyyyz, g_zz_0_y_yyyyzz, g_zz_0_y_yyyzz, g_zz_0_y_yyyzzz, g_zz_0_y_yyzzz, g_zz_0_y_yyzzzz, g_zz_0_y_yzzzz, g_zz_0_y_yzzzzz, g_zz_0_y_zzzzz, g_zz_0_yy_xxxxx, g_zz_0_yy_xxxxy, g_zz_0_yy_xxxxz, g_zz_0_yy_xxxyy, g_zz_0_yy_xxxyz, g_zz_0_yy_xxxzz, g_zz_0_yy_xxyyy, g_zz_0_yy_xxyyz, g_zz_0_yy_xxyzz, g_zz_0_yy_xxzzz, g_zz_0_yy_xyyyy, g_zz_0_yy_xyyyz, g_zz_0_yy_xyyzz, g_zz_0_yy_xyzzz, g_zz_0_yy_xzzzz, g_zz_0_yy_yyyyy, g_zz_0_yy_yyyyz, g_zz_0_yy_yyyzz, g_zz_0_yy_yyzzz, g_zz_0_yy_yzzzz, g_zz_0_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yy_xxxxx[k] = -g_zz_0_y_xxxxx[k] * ab_y + g_zz_0_y_xxxxxy[k];

                g_zz_0_yy_xxxxy[k] = -g_zz_0_y_xxxxy[k] * ab_y + g_zz_0_y_xxxxyy[k];

                g_zz_0_yy_xxxxz[k] = -g_zz_0_y_xxxxz[k] * ab_y + g_zz_0_y_xxxxyz[k];

                g_zz_0_yy_xxxyy[k] = -g_zz_0_y_xxxyy[k] * ab_y + g_zz_0_y_xxxyyy[k];

                g_zz_0_yy_xxxyz[k] = -g_zz_0_y_xxxyz[k] * ab_y + g_zz_0_y_xxxyyz[k];

                g_zz_0_yy_xxxzz[k] = -g_zz_0_y_xxxzz[k] * ab_y + g_zz_0_y_xxxyzz[k];

                g_zz_0_yy_xxyyy[k] = -g_zz_0_y_xxyyy[k] * ab_y + g_zz_0_y_xxyyyy[k];

                g_zz_0_yy_xxyyz[k] = -g_zz_0_y_xxyyz[k] * ab_y + g_zz_0_y_xxyyyz[k];

                g_zz_0_yy_xxyzz[k] = -g_zz_0_y_xxyzz[k] * ab_y + g_zz_0_y_xxyyzz[k];

                g_zz_0_yy_xxzzz[k] = -g_zz_0_y_xxzzz[k] * ab_y + g_zz_0_y_xxyzzz[k];

                g_zz_0_yy_xyyyy[k] = -g_zz_0_y_xyyyy[k] * ab_y + g_zz_0_y_xyyyyy[k];

                g_zz_0_yy_xyyyz[k] = -g_zz_0_y_xyyyz[k] * ab_y + g_zz_0_y_xyyyyz[k];

                g_zz_0_yy_xyyzz[k] = -g_zz_0_y_xyyzz[k] * ab_y + g_zz_0_y_xyyyzz[k];

                g_zz_0_yy_xyzzz[k] = -g_zz_0_y_xyzzz[k] * ab_y + g_zz_0_y_xyyzzz[k];

                g_zz_0_yy_xzzzz[k] = -g_zz_0_y_xzzzz[k] * ab_y + g_zz_0_y_xyzzzz[k];

                g_zz_0_yy_yyyyy[k] = -g_zz_0_y_yyyyy[k] * ab_y + g_zz_0_y_yyyyyy[k];

                g_zz_0_yy_yyyyz[k] = -g_zz_0_y_yyyyz[k] * ab_y + g_zz_0_y_yyyyyz[k];

                g_zz_0_yy_yyyzz[k] = -g_zz_0_y_yyyzz[k] * ab_y + g_zz_0_y_yyyyzz[k];

                g_zz_0_yy_yyzzz[k] = -g_zz_0_y_yyzzz[k] * ab_y + g_zz_0_y_yyyzzz[k];

                g_zz_0_yy_yzzzz[k] = -g_zz_0_y_yzzzz[k] * ab_y + g_zz_0_y_yyzzzz[k];

                g_zz_0_yy_zzzzz[k] = -g_zz_0_y_zzzzz[k] * ab_y + g_zz_0_y_yzzzzz[k];
            }

            /// Set up 714-735 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yz_xxxxx = cbuffer.data(dh_geom_20_off + 714 * ccomps * dcomps);

            auto g_zz_0_yz_xxxxy = cbuffer.data(dh_geom_20_off + 715 * ccomps * dcomps);

            auto g_zz_0_yz_xxxxz = cbuffer.data(dh_geom_20_off + 716 * ccomps * dcomps);

            auto g_zz_0_yz_xxxyy = cbuffer.data(dh_geom_20_off + 717 * ccomps * dcomps);

            auto g_zz_0_yz_xxxyz = cbuffer.data(dh_geom_20_off + 718 * ccomps * dcomps);

            auto g_zz_0_yz_xxxzz = cbuffer.data(dh_geom_20_off + 719 * ccomps * dcomps);

            auto g_zz_0_yz_xxyyy = cbuffer.data(dh_geom_20_off + 720 * ccomps * dcomps);

            auto g_zz_0_yz_xxyyz = cbuffer.data(dh_geom_20_off + 721 * ccomps * dcomps);

            auto g_zz_0_yz_xxyzz = cbuffer.data(dh_geom_20_off + 722 * ccomps * dcomps);

            auto g_zz_0_yz_xxzzz = cbuffer.data(dh_geom_20_off + 723 * ccomps * dcomps);

            auto g_zz_0_yz_xyyyy = cbuffer.data(dh_geom_20_off + 724 * ccomps * dcomps);

            auto g_zz_0_yz_xyyyz = cbuffer.data(dh_geom_20_off + 725 * ccomps * dcomps);

            auto g_zz_0_yz_xyyzz = cbuffer.data(dh_geom_20_off + 726 * ccomps * dcomps);

            auto g_zz_0_yz_xyzzz = cbuffer.data(dh_geom_20_off + 727 * ccomps * dcomps);

            auto g_zz_0_yz_xzzzz = cbuffer.data(dh_geom_20_off + 728 * ccomps * dcomps);

            auto g_zz_0_yz_yyyyy = cbuffer.data(dh_geom_20_off + 729 * ccomps * dcomps);

            auto g_zz_0_yz_yyyyz = cbuffer.data(dh_geom_20_off + 730 * ccomps * dcomps);

            auto g_zz_0_yz_yyyzz = cbuffer.data(dh_geom_20_off + 731 * ccomps * dcomps);

            auto g_zz_0_yz_yyzzz = cbuffer.data(dh_geom_20_off + 732 * ccomps * dcomps);

            auto g_zz_0_yz_yzzzz = cbuffer.data(dh_geom_20_off + 733 * ccomps * dcomps);

            auto g_zz_0_yz_zzzzz = cbuffer.data(dh_geom_20_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yz_xxxxx, g_zz_0_yz_xxxxy, g_zz_0_yz_xxxxz, g_zz_0_yz_xxxyy, g_zz_0_yz_xxxyz, g_zz_0_yz_xxxzz, g_zz_0_yz_xxyyy, g_zz_0_yz_xxyyz, g_zz_0_yz_xxyzz, g_zz_0_yz_xxzzz, g_zz_0_yz_xyyyy, g_zz_0_yz_xyyyz, g_zz_0_yz_xyyzz, g_zz_0_yz_xyzzz, g_zz_0_yz_xzzzz, g_zz_0_yz_yyyyy, g_zz_0_yz_yyyyz, g_zz_0_yz_yyyzz, g_zz_0_yz_yyzzz, g_zz_0_yz_yzzzz, g_zz_0_yz_zzzzz, g_zz_0_z_xxxxx, g_zz_0_z_xxxxxy, g_zz_0_z_xxxxy, g_zz_0_z_xxxxyy, g_zz_0_z_xxxxyz, g_zz_0_z_xxxxz, g_zz_0_z_xxxyy, g_zz_0_z_xxxyyy, g_zz_0_z_xxxyyz, g_zz_0_z_xxxyz, g_zz_0_z_xxxyzz, g_zz_0_z_xxxzz, g_zz_0_z_xxyyy, g_zz_0_z_xxyyyy, g_zz_0_z_xxyyyz, g_zz_0_z_xxyyz, g_zz_0_z_xxyyzz, g_zz_0_z_xxyzz, g_zz_0_z_xxyzzz, g_zz_0_z_xxzzz, g_zz_0_z_xyyyy, g_zz_0_z_xyyyyy, g_zz_0_z_xyyyyz, g_zz_0_z_xyyyz, g_zz_0_z_xyyyzz, g_zz_0_z_xyyzz, g_zz_0_z_xyyzzz, g_zz_0_z_xyzzz, g_zz_0_z_xyzzzz, g_zz_0_z_xzzzz, g_zz_0_z_yyyyy, g_zz_0_z_yyyyyy, g_zz_0_z_yyyyyz, g_zz_0_z_yyyyz, g_zz_0_z_yyyyzz, g_zz_0_z_yyyzz, g_zz_0_z_yyyzzz, g_zz_0_z_yyzzz, g_zz_0_z_yyzzzz, g_zz_0_z_yzzzz, g_zz_0_z_yzzzzz, g_zz_0_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yz_xxxxx[k] = -g_zz_0_z_xxxxx[k] * ab_y + g_zz_0_z_xxxxxy[k];

                g_zz_0_yz_xxxxy[k] = -g_zz_0_z_xxxxy[k] * ab_y + g_zz_0_z_xxxxyy[k];

                g_zz_0_yz_xxxxz[k] = -g_zz_0_z_xxxxz[k] * ab_y + g_zz_0_z_xxxxyz[k];

                g_zz_0_yz_xxxyy[k] = -g_zz_0_z_xxxyy[k] * ab_y + g_zz_0_z_xxxyyy[k];

                g_zz_0_yz_xxxyz[k] = -g_zz_0_z_xxxyz[k] * ab_y + g_zz_0_z_xxxyyz[k];

                g_zz_0_yz_xxxzz[k] = -g_zz_0_z_xxxzz[k] * ab_y + g_zz_0_z_xxxyzz[k];

                g_zz_0_yz_xxyyy[k] = -g_zz_0_z_xxyyy[k] * ab_y + g_zz_0_z_xxyyyy[k];

                g_zz_0_yz_xxyyz[k] = -g_zz_0_z_xxyyz[k] * ab_y + g_zz_0_z_xxyyyz[k];

                g_zz_0_yz_xxyzz[k] = -g_zz_0_z_xxyzz[k] * ab_y + g_zz_0_z_xxyyzz[k];

                g_zz_0_yz_xxzzz[k] = -g_zz_0_z_xxzzz[k] * ab_y + g_zz_0_z_xxyzzz[k];

                g_zz_0_yz_xyyyy[k] = -g_zz_0_z_xyyyy[k] * ab_y + g_zz_0_z_xyyyyy[k];

                g_zz_0_yz_xyyyz[k] = -g_zz_0_z_xyyyz[k] * ab_y + g_zz_0_z_xyyyyz[k];

                g_zz_0_yz_xyyzz[k] = -g_zz_0_z_xyyzz[k] * ab_y + g_zz_0_z_xyyyzz[k];

                g_zz_0_yz_xyzzz[k] = -g_zz_0_z_xyzzz[k] * ab_y + g_zz_0_z_xyyzzz[k];

                g_zz_0_yz_xzzzz[k] = -g_zz_0_z_xzzzz[k] * ab_y + g_zz_0_z_xyzzzz[k];

                g_zz_0_yz_yyyyy[k] = -g_zz_0_z_yyyyy[k] * ab_y + g_zz_0_z_yyyyyy[k];

                g_zz_0_yz_yyyyz[k] = -g_zz_0_z_yyyyz[k] * ab_y + g_zz_0_z_yyyyyz[k];

                g_zz_0_yz_yyyzz[k] = -g_zz_0_z_yyyzz[k] * ab_y + g_zz_0_z_yyyyzz[k];

                g_zz_0_yz_yyzzz[k] = -g_zz_0_z_yyzzz[k] * ab_y + g_zz_0_z_yyyzzz[k];

                g_zz_0_yz_yzzzz[k] = -g_zz_0_z_yzzzz[k] * ab_y + g_zz_0_z_yyzzzz[k];

                g_zz_0_yz_zzzzz[k] = -g_zz_0_z_zzzzz[k] * ab_y + g_zz_0_z_yzzzzz[k];
            }

            /// Set up 735-756 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zz_xxxxx = cbuffer.data(dh_geom_20_off + 735 * ccomps * dcomps);

            auto g_zz_0_zz_xxxxy = cbuffer.data(dh_geom_20_off + 736 * ccomps * dcomps);

            auto g_zz_0_zz_xxxxz = cbuffer.data(dh_geom_20_off + 737 * ccomps * dcomps);

            auto g_zz_0_zz_xxxyy = cbuffer.data(dh_geom_20_off + 738 * ccomps * dcomps);

            auto g_zz_0_zz_xxxyz = cbuffer.data(dh_geom_20_off + 739 * ccomps * dcomps);

            auto g_zz_0_zz_xxxzz = cbuffer.data(dh_geom_20_off + 740 * ccomps * dcomps);

            auto g_zz_0_zz_xxyyy = cbuffer.data(dh_geom_20_off + 741 * ccomps * dcomps);

            auto g_zz_0_zz_xxyyz = cbuffer.data(dh_geom_20_off + 742 * ccomps * dcomps);

            auto g_zz_0_zz_xxyzz = cbuffer.data(dh_geom_20_off + 743 * ccomps * dcomps);

            auto g_zz_0_zz_xxzzz = cbuffer.data(dh_geom_20_off + 744 * ccomps * dcomps);

            auto g_zz_0_zz_xyyyy = cbuffer.data(dh_geom_20_off + 745 * ccomps * dcomps);

            auto g_zz_0_zz_xyyyz = cbuffer.data(dh_geom_20_off + 746 * ccomps * dcomps);

            auto g_zz_0_zz_xyyzz = cbuffer.data(dh_geom_20_off + 747 * ccomps * dcomps);

            auto g_zz_0_zz_xyzzz = cbuffer.data(dh_geom_20_off + 748 * ccomps * dcomps);

            auto g_zz_0_zz_xzzzz = cbuffer.data(dh_geom_20_off + 749 * ccomps * dcomps);

            auto g_zz_0_zz_yyyyy = cbuffer.data(dh_geom_20_off + 750 * ccomps * dcomps);

            auto g_zz_0_zz_yyyyz = cbuffer.data(dh_geom_20_off + 751 * ccomps * dcomps);

            auto g_zz_0_zz_yyyzz = cbuffer.data(dh_geom_20_off + 752 * ccomps * dcomps);

            auto g_zz_0_zz_yyzzz = cbuffer.data(dh_geom_20_off + 753 * ccomps * dcomps);

            auto g_zz_0_zz_yzzzz = cbuffer.data(dh_geom_20_off + 754 * ccomps * dcomps);

            auto g_zz_0_zz_zzzzz = cbuffer.data(dh_geom_20_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxxx, g_z_0_z_xxxxy, g_z_0_z_xxxxz, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxzz, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyzz, g_z_0_z_xxzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyzz, g_z_0_z_xyzzz, g_z_0_z_xzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyzz, g_z_0_z_yyzzz, g_z_0_z_yzzzz, g_z_0_z_zzzzz, g_zz_0_z_xxxxx, g_zz_0_z_xxxxxz, g_zz_0_z_xxxxy, g_zz_0_z_xxxxyz, g_zz_0_z_xxxxz, g_zz_0_z_xxxxzz, g_zz_0_z_xxxyy, g_zz_0_z_xxxyyz, g_zz_0_z_xxxyz, g_zz_0_z_xxxyzz, g_zz_0_z_xxxzz, g_zz_0_z_xxxzzz, g_zz_0_z_xxyyy, g_zz_0_z_xxyyyz, g_zz_0_z_xxyyz, g_zz_0_z_xxyyzz, g_zz_0_z_xxyzz, g_zz_0_z_xxyzzz, g_zz_0_z_xxzzz, g_zz_0_z_xxzzzz, g_zz_0_z_xyyyy, g_zz_0_z_xyyyyz, g_zz_0_z_xyyyz, g_zz_0_z_xyyyzz, g_zz_0_z_xyyzz, g_zz_0_z_xyyzzz, g_zz_0_z_xyzzz, g_zz_0_z_xyzzzz, g_zz_0_z_xzzzz, g_zz_0_z_xzzzzz, g_zz_0_z_yyyyy, g_zz_0_z_yyyyyz, g_zz_0_z_yyyyz, g_zz_0_z_yyyyzz, g_zz_0_z_yyyzz, g_zz_0_z_yyyzzz, g_zz_0_z_yyzzz, g_zz_0_z_yyzzzz, g_zz_0_z_yzzzz, g_zz_0_z_yzzzzz, g_zz_0_z_zzzzz, g_zz_0_z_zzzzzz, g_zz_0_zz_xxxxx, g_zz_0_zz_xxxxy, g_zz_0_zz_xxxxz, g_zz_0_zz_xxxyy, g_zz_0_zz_xxxyz, g_zz_0_zz_xxxzz, g_zz_0_zz_xxyyy, g_zz_0_zz_xxyyz, g_zz_0_zz_xxyzz, g_zz_0_zz_xxzzz, g_zz_0_zz_xyyyy, g_zz_0_zz_xyyyz, g_zz_0_zz_xyyzz, g_zz_0_zz_xyzzz, g_zz_0_zz_xzzzz, g_zz_0_zz_yyyyy, g_zz_0_zz_yyyyz, g_zz_0_zz_yyyzz, g_zz_0_zz_yyzzz, g_zz_0_zz_yzzzz, g_zz_0_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zz_xxxxx[k] = -2.0 * g_z_0_z_xxxxx[k] - g_zz_0_z_xxxxx[k] * ab_z + g_zz_0_z_xxxxxz[k];

                g_zz_0_zz_xxxxy[k] = -2.0 * g_z_0_z_xxxxy[k] - g_zz_0_z_xxxxy[k] * ab_z + g_zz_0_z_xxxxyz[k];

                g_zz_0_zz_xxxxz[k] = -2.0 * g_z_0_z_xxxxz[k] - g_zz_0_z_xxxxz[k] * ab_z + g_zz_0_z_xxxxzz[k];

                g_zz_0_zz_xxxyy[k] = -2.0 * g_z_0_z_xxxyy[k] - g_zz_0_z_xxxyy[k] * ab_z + g_zz_0_z_xxxyyz[k];

                g_zz_0_zz_xxxyz[k] = -2.0 * g_z_0_z_xxxyz[k] - g_zz_0_z_xxxyz[k] * ab_z + g_zz_0_z_xxxyzz[k];

                g_zz_0_zz_xxxzz[k] = -2.0 * g_z_0_z_xxxzz[k] - g_zz_0_z_xxxzz[k] * ab_z + g_zz_0_z_xxxzzz[k];

                g_zz_0_zz_xxyyy[k] = -2.0 * g_z_0_z_xxyyy[k] - g_zz_0_z_xxyyy[k] * ab_z + g_zz_0_z_xxyyyz[k];

                g_zz_0_zz_xxyyz[k] = -2.0 * g_z_0_z_xxyyz[k] - g_zz_0_z_xxyyz[k] * ab_z + g_zz_0_z_xxyyzz[k];

                g_zz_0_zz_xxyzz[k] = -2.0 * g_z_0_z_xxyzz[k] - g_zz_0_z_xxyzz[k] * ab_z + g_zz_0_z_xxyzzz[k];

                g_zz_0_zz_xxzzz[k] = -2.0 * g_z_0_z_xxzzz[k] - g_zz_0_z_xxzzz[k] * ab_z + g_zz_0_z_xxzzzz[k];

                g_zz_0_zz_xyyyy[k] = -2.0 * g_z_0_z_xyyyy[k] - g_zz_0_z_xyyyy[k] * ab_z + g_zz_0_z_xyyyyz[k];

                g_zz_0_zz_xyyyz[k] = -2.0 * g_z_0_z_xyyyz[k] - g_zz_0_z_xyyyz[k] * ab_z + g_zz_0_z_xyyyzz[k];

                g_zz_0_zz_xyyzz[k] = -2.0 * g_z_0_z_xyyzz[k] - g_zz_0_z_xyyzz[k] * ab_z + g_zz_0_z_xyyzzz[k];

                g_zz_0_zz_xyzzz[k] = -2.0 * g_z_0_z_xyzzz[k] - g_zz_0_z_xyzzz[k] * ab_z + g_zz_0_z_xyzzzz[k];

                g_zz_0_zz_xzzzz[k] = -2.0 * g_z_0_z_xzzzz[k] - g_zz_0_z_xzzzz[k] * ab_z + g_zz_0_z_xzzzzz[k];

                g_zz_0_zz_yyyyy[k] = -2.0 * g_z_0_z_yyyyy[k] - g_zz_0_z_yyyyy[k] * ab_z + g_zz_0_z_yyyyyz[k];

                g_zz_0_zz_yyyyz[k] = -2.0 * g_z_0_z_yyyyz[k] - g_zz_0_z_yyyyz[k] * ab_z + g_zz_0_z_yyyyzz[k];

                g_zz_0_zz_yyyzz[k] = -2.0 * g_z_0_z_yyyzz[k] - g_zz_0_z_yyyzz[k] * ab_z + g_zz_0_z_yyyzzz[k];

                g_zz_0_zz_yyzzz[k] = -2.0 * g_z_0_z_yyzzz[k] - g_zz_0_z_yyzzz[k] * ab_z + g_zz_0_z_yyzzzz[k];

                g_zz_0_zz_yzzzz[k] = -2.0 * g_z_0_z_yzzzz[k] - g_zz_0_z_yzzzz[k] * ab_z + g_zz_0_z_yzzzzz[k];

                g_zz_0_zz_zzzzz[k] = -2.0 * g_z_0_z_zzzzz[k] - g_zz_0_z_zzzzz[k] * ab_z + g_zz_0_z_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

