#include "ElectronRepulsionGeom1100ContrRecDHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom11_hrr_electron_repulsion_dhxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_11_dhxx,
                                            const size_t idx_geom_01_phxx,
                                            const size_t idx_geom_10_phxx,
                                            const size_t idx_geom_11_phxx,
                                            const size_t idx_geom_11_pixx,
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

            const auto ph_geom_01_off = idx_geom_01_phxx + i * dcomps + j;

            auto g_0_x_x_xxxxx = cbuffer.data(ph_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_x_xxxxy = cbuffer.data(ph_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_x_xxxxz = cbuffer.data(ph_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_x_xxxyy = cbuffer.data(ph_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_x_xxxyz = cbuffer.data(ph_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_x_xxxzz = cbuffer.data(ph_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_x_xxyyy = cbuffer.data(ph_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_x_xxyyz = cbuffer.data(ph_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_x_xxyzz = cbuffer.data(ph_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_x_xxzzz = cbuffer.data(ph_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_x_xyyyy = cbuffer.data(ph_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_x_xyyyz = cbuffer.data(ph_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_x_xyyzz = cbuffer.data(ph_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_x_xyzzz = cbuffer.data(ph_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_x_xzzzz = cbuffer.data(ph_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_x_yyyyy = cbuffer.data(ph_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_x_yyyyz = cbuffer.data(ph_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_x_yyyzz = cbuffer.data(ph_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_x_yyzzz = cbuffer.data(ph_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_x_yzzzz = cbuffer.data(ph_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_x_zzzzz = cbuffer.data(ph_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_y_xxxxx = cbuffer.data(ph_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_y_xxxxy = cbuffer.data(ph_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_y_xxxxz = cbuffer.data(ph_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_y_xxxyy = cbuffer.data(ph_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_y_xxxyz = cbuffer.data(ph_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_y_xxxzz = cbuffer.data(ph_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_y_xxyyy = cbuffer.data(ph_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_y_xxyyz = cbuffer.data(ph_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_y_xxyzz = cbuffer.data(ph_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_y_xxzzz = cbuffer.data(ph_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_y_xyyyy = cbuffer.data(ph_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_y_xyyyz = cbuffer.data(ph_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_y_xyyzz = cbuffer.data(ph_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_y_xyzzz = cbuffer.data(ph_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_y_xzzzz = cbuffer.data(ph_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_y_yyyyy = cbuffer.data(ph_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_y_yyyyz = cbuffer.data(ph_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_y_yyyzz = cbuffer.data(ph_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_y_yyzzz = cbuffer.data(ph_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_y_yzzzz = cbuffer.data(ph_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_y_zzzzz = cbuffer.data(ph_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_z_xxxxx = cbuffer.data(ph_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_z_xxxxy = cbuffer.data(ph_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_z_xxxxz = cbuffer.data(ph_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_z_xxxyy = cbuffer.data(ph_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_z_xxxyz = cbuffer.data(ph_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_z_xxxzz = cbuffer.data(ph_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_z_xxyyy = cbuffer.data(ph_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_z_xxyyz = cbuffer.data(ph_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_z_xxyzz = cbuffer.data(ph_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_z_xxzzz = cbuffer.data(ph_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_z_xyyyy = cbuffer.data(ph_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_z_xyyyz = cbuffer.data(ph_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_z_xyyzz = cbuffer.data(ph_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_z_xyzzz = cbuffer.data(ph_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_z_xzzzz = cbuffer.data(ph_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_z_yyyyy = cbuffer.data(ph_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_z_yyyyz = cbuffer.data(ph_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_z_yyyzz = cbuffer.data(ph_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_z_yyzzz = cbuffer.data(ph_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_z_yzzzz = cbuffer.data(ph_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_z_zzzzz = cbuffer.data(ph_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_y_x_xxxxx = cbuffer.data(ph_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_y_x_xxxxy = cbuffer.data(ph_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_y_x_xxxxz = cbuffer.data(ph_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_y_x_xxxyy = cbuffer.data(ph_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_y_x_xxxyz = cbuffer.data(ph_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_y_x_xxxzz = cbuffer.data(ph_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_y_x_xxyyy = cbuffer.data(ph_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_y_x_xxyyz = cbuffer.data(ph_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_y_x_xxyzz = cbuffer.data(ph_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_y_x_xxzzz = cbuffer.data(ph_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_y_x_xyyyy = cbuffer.data(ph_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_y_x_xyyyz = cbuffer.data(ph_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_y_x_xyyzz = cbuffer.data(ph_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_y_x_xyzzz = cbuffer.data(ph_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_y_x_xzzzz = cbuffer.data(ph_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_y_x_yyyyy = cbuffer.data(ph_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_y_x_yyyyz = cbuffer.data(ph_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_y_x_yyyzz = cbuffer.data(ph_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_y_x_yyzzz = cbuffer.data(ph_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_y_x_yzzzz = cbuffer.data(ph_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_y_x_zzzzz = cbuffer.data(ph_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_y_y_xxxxx = cbuffer.data(ph_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_y_y_xxxxy = cbuffer.data(ph_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_y_y_xxxxz = cbuffer.data(ph_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_y_y_xxxyy = cbuffer.data(ph_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_y_y_xxxyz = cbuffer.data(ph_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_y_y_xxxzz = cbuffer.data(ph_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_y_y_xxyyy = cbuffer.data(ph_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_y_y_xxyyz = cbuffer.data(ph_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_y_y_xxyzz = cbuffer.data(ph_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_y_y_xxzzz = cbuffer.data(ph_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_y_y_xyyyy = cbuffer.data(ph_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_y_y_xyyyz = cbuffer.data(ph_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_y_y_xyyzz = cbuffer.data(ph_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_y_y_xyzzz = cbuffer.data(ph_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_y_y_xzzzz = cbuffer.data(ph_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_y_y_yyyyy = cbuffer.data(ph_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_y_y_yyyyz = cbuffer.data(ph_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_y_y_yyyzz = cbuffer.data(ph_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_y_y_yyzzz = cbuffer.data(ph_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_y_y_yzzzz = cbuffer.data(ph_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_y_y_zzzzz = cbuffer.data(ph_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_y_z_xxxxx = cbuffer.data(ph_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_y_z_xxxxy = cbuffer.data(ph_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_y_z_xxxxz = cbuffer.data(ph_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_y_z_xxxyy = cbuffer.data(ph_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_y_z_xxxyz = cbuffer.data(ph_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_y_z_xxxzz = cbuffer.data(ph_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_y_z_xxyyy = cbuffer.data(ph_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_y_z_xxyyz = cbuffer.data(ph_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_y_z_xxyzz = cbuffer.data(ph_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_y_z_xxzzz = cbuffer.data(ph_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_y_z_xyyyy = cbuffer.data(ph_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_y_z_xyyyz = cbuffer.data(ph_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_y_z_xyyzz = cbuffer.data(ph_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_y_z_xyzzz = cbuffer.data(ph_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_y_z_xzzzz = cbuffer.data(ph_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_y_z_yyyyy = cbuffer.data(ph_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_y_z_yyyyz = cbuffer.data(ph_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_y_z_yyyzz = cbuffer.data(ph_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_y_z_yyzzz = cbuffer.data(ph_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_y_z_yzzzz = cbuffer.data(ph_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_y_z_zzzzz = cbuffer.data(ph_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_z_x_xxxxx = cbuffer.data(ph_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_z_x_xxxxy = cbuffer.data(ph_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_z_x_xxxxz = cbuffer.data(ph_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_z_x_xxxyy = cbuffer.data(ph_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_z_x_xxxyz = cbuffer.data(ph_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_z_x_xxxzz = cbuffer.data(ph_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_z_x_xxyyy = cbuffer.data(ph_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_z_x_xxyyz = cbuffer.data(ph_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_z_x_xxyzz = cbuffer.data(ph_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_z_x_xxzzz = cbuffer.data(ph_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_z_x_xyyyy = cbuffer.data(ph_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_z_x_xyyyz = cbuffer.data(ph_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_z_x_xyyzz = cbuffer.data(ph_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_z_x_xyzzz = cbuffer.data(ph_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_z_x_xzzzz = cbuffer.data(ph_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_z_x_yyyyy = cbuffer.data(ph_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_z_x_yyyyz = cbuffer.data(ph_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_z_x_yyyzz = cbuffer.data(ph_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_z_x_yyzzz = cbuffer.data(ph_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_z_x_yzzzz = cbuffer.data(ph_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_z_x_zzzzz = cbuffer.data(ph_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_z_y_xxxxx = cbuffer.data(ph_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_z_y_xxxxy = cbuffer.data(ph_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_z_y_xxxxz = cbuffer.data(ph_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_z_y_xxxyy = cbuffer.data(ph_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_z_y_xxxyz = cbuffer.data(ph_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_z_y_xxxzz = cbuffer.data(ph_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_z_y_xxyyy = cbuffer.data(ph_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_z_y_xxyyz = cbuffer.data(ph_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_z_y_xxyzz = cbuffer.data(ph_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_z_y_xxzzz = cbuffer.data(ph_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_z_y_xyyyy = cbuffer.data(ph_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_z_y_xyyyz = cbuffer.data(ph_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_z_y_xyyzz = cbuffer.data(ph_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_z_y_xyzzz = cbuffer.data(ph_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_z_y_xzzzz = cbuffer.data(ph_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_z_y_yyyyy = cbuffer.data(ph_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_z_y_yyyyz = cbuffer.data(ph_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_z_y_yyyzz = cbuffer.data(ph_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_z_y_yyzzz = cbuffer.data(ph_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_z_y_yzzzz = cbuffer.data(ph_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_z_y_zzzzz = cbuffer.data(ph_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_z_z_xxxxx = cbuffer.data(ph_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_z_z_xxxxy = cbuffer.data(ph_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_z_z_xxxxz = cbuffer.data(ph_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_z_z_xxxyy = cbuffer.data(ph_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_z_z_xxxyz = cbuffer.data(ph_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_z_z_xxxzz = cbuffer.data(ph_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_z_z_xxyyy = cbuffer.data(ph_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_z_z_xxyyz = cbuffer.data(ph_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_z_z_xxyzz = cbuffer.data(ph_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_z_z_xxzzz = cbuffer.data(ph_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_z_z_xyyyy = cbuffer.data(ph_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_z_z_xyyyz = cbuffer.data(ph_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_z_z_xyyzz = cbuffer.data(ph_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_z_z_xyzzz = cbuffer.data(ph_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_z_z_xzzzz = cbuffer.data(ph_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_z_z_yyyyy = cbuffer.data(ph_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_z_z_yyyyz = cbuffer.data(ph_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_z_z_yyyzz = cbuffer.data(ph_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_z_z_yyzzz = cbuffer.data(ph_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_z_z_yzzzz = cbuffer.data(ph_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_z_z_zzzzz = cbuffer.data(ph_geom_01_off + 188 * ccomps * dcomps);

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

            const auto ph_geom_11_off = idx_geom_11_phxx + i * dcomps + j;

            auto g_x_x_x_xxxxx = cbuffer.data(ph_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_x_xxxxy = cbuffer.data(ph_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_x_xxxxz = cbuffer.data(ph_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_x_xxxyy = cbuffer.data(ph_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_x_xxxyz = cbuffer.data(ph_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_x_xxxzz = cbuffer.data(ph_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_x_xxyyy = cbuffer.data(ph_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_x_xxyyz = cbuffer.data(ph_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_x_xxyzz = cbuffer.data(ph_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_x_xxzzz = cbuffer.data(ph_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_x_xyyyy = cbuffer.data(ph_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_x_xyyyz = cbuffer.data(ph_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_x_xyyzz = cbuffer.data(ph_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_x_xyzzz = cbuffer.data(ph_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_x_xzzzz = cbuffer.data(ph_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_x_yyyyy = cbuffer.data(ph_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_x_yyyyz = cbuffer.data(ph_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_x_yyyzz = cbuffer.data(ph_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_x_yyzzz = cbuffer.data(ph_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_x_yzzzz = cbuffer.data(ph_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_x_zzzzz = cbuffer.data(ph_geom_11_off + 20 * ccomps * dcomps);

            auto g_x_x_y_xxxxx = cbuffer.data(ph_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_y_xxxxy = cbuffer.data(ph_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_y_xxxxz = cbuffer.data(ph_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_y_xxxyy = cbuffer.data(ph_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_y_xxxyz = cbuffer.data(ph_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_y_xxxzz = cbuffer.data(ph_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_y_xxyyy = cbuffer.data(ph_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_y_xxyyz = cbuffer.data(ph_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_y_xxyzz = cbuffer.data(ph_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_y_xxzzz = cbuffer.data(ph_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_y_xyyyy = cbuffer.data(ph_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_y_xyyyz = cbuffer.data(ph_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_y_xyyzz = cbuffer.data(ph_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_y_xyzzz = cbuffer.data(ph_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_y_xzzzz = cbuffer.data(ph_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_y_yyyyy = cbuffer.data(ph_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_y_yyyyz = cbuffer.data(ph_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_y_yyyzz = cbuffer.data(ph_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_y_yyzzz = cbuffer.data(ph_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_y_yzzzz = cbuffer.data(ph_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_y_zzzzz = cbuffer.data(ph_geom_11_off + 41 * ccomps * dcomps);

            auto g_x_x_z_xxxxx = cbuffer.data(ph_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_z_xxxxy = cbuffer.data(ph_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_z_xxxxz = cbuffer.data(ph_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_z_xxxyy = cbuffer.data(ph_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_z_xxxyz = cbuffer.data(ph_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_z_xxxzz = cbuffer.data(ph_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_z_xxyyy = cbuffer.data(ph_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_z_xxyyz = cbuffer.data(ph_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_z_xxyzz = cbuffer.data(ph_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_z_xxzzz = cbuffer.data(ph_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_z_xyyyy = cbuffer.data(ph_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_z_xyyyz = cbuffer.data(ph_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_z_xyyzz = cbuffer.data(ph_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_z_xyzzz = cbuffer.data(ph_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_z_xzzzz = cbuffer.data(ph_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_z_yyyyy = cbuffer.data(ph_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_z_yyyyz = cbuffer.data(ph_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_z_yyyzz = cbuffer.data(ph_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_z_yyzzz = cbuffer.data(ph_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_z_yzzzz = cbuffer.data(ph_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_z_zzzzz = cbuffer.data(ph_geom_11_off + 62 * ccomps * dcomps);

            auto g_x_y_x_xxxxx = cbuffer.data(ph_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_y_x_xxxxy = cbuffer.data(ph_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_y_x_xxxxz = cbuffer.data(ph_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_y_x_xxxyy = cbuffer.data(ph_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_y_x_xxxyz = cbuffer.data(ph_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_y_x_xxxzz = cbuffer.data(ph_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_y_x_xxyyy = cbuffer.data(ph_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_y_x_xxyyz = cbuffer.data(ph_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_y_x_xxyzz = cbuffer.data(ph_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_y_x_xxzzz = cbuffer.data(ph_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_y_x_xyyyy = cbuffer.data(ph_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_y_x_xyyyz = cbuffer.data(ph_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_y_x_xyyzz = cbuffer.data(ph_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_y_x_xyzzz = cbuffer.data(ph_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_y_x_xzzzz = cbuffer.data(ph_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_y_x_yyyyy = cbuffer.data(ph_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_y_x_yyyyz = cbuffer.data(ph_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_y_x_yyyzz = cbuffer.data(ph_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_y_x_yyzzz = cbuffer.data(ph_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_y_x_yzzzz = cbuffer.data(ph_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_y_x_zzzzz = cbuffer.data(ph_geom_11_off + 83 * ccomps * dcomps);

            auto g_x_y_y_xxxxx = cbuffer.data(ph_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_y_y_xxxxy = cbuffer.data(ph_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_y_y_xxxxz = cbuffer.data(ph_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_y_y_xxxyy = cbuffer.data(ph_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_y_y_xxxyz = cbuffer.data(ph_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_y_y_xxxzz = cbuffer.data(ph_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_y_y_xxyyy = cbuffer.data(ph_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_y_y_xxyyz = cbuffer.data(ph_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_y_y_xxyzz = cbuffer.data(ph_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_y_y_xxzzz = cbuffer.data(ph_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_y_y_xyyyy = cbuffer.data(ph_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_y_y_xyyyz = cbuffer.data(ph_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_y_y_xyyzz = cbuffer.data(ph_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_y_y_xyzzz = cbuffer.data(ph_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_y_y_xzzzz = cbuffer.data(ph_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_y_y_yyyyy = cbuffer.data(ph_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_y_y_yyyyz = cbuffer.data(ph_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_y_y_yyyzz = cbuffer.data(ph_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_y_y_yyzzz = cbuffer.data(ph_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_y_y_yzzzz = cbuffer.data(ph_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_y_y_zzzzz = cbuffer.data(ph_geom_11_off + 104 * ccomps * dcomps);

            auto g_x_y_z_xxxxx = cbuffer.data(ph_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_y_z_xxxxy = cbuffer.data(ph_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_y_z_xxxxz = cbuffer.data(ph_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_y_z_xxxyy = cbuffer.data(ph_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_y_z_xxxyz = cbuffer.data(ph_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_y_z_xxxzz = cbuffer.data(ph_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_y_z_xxyyy = cbuffer.data(ph_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_y_z_xxyyz = cbuffer.data(ph_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_y_z_xxyzz = cbuffer.data(ph_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_y_z_xxzzz = cbuffer.data(ph_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_y_z_xyyyy = cbuffer.data(ph_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_y_z_xyyyz = cbuffer.data(ph_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_y_z_xyyzz = cbuffer.data(ph_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_y_z_xyzzz = cbuffer.data(ph_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_y_z_xzzzz = cbuffer.data(ph_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_y_z_yyyyy = cbuffer.data(ph_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_y_z_yyyyz = cbuffer.data(ph_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_y_z_yyyzz = cbuffer.data(ph_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_y_z_yyzzz = cbuffer.data(ph_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_y_z_yzzzz = cbuffer.data(ph_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_y_z_zzzzz = cbuffer.data(ph_geom_11_off + 125 * ccomps * dcomps);

            auto g_x_z_x_xxxxx = cbuffer.data(ph_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_z_x_xxxxy = cbuffer.data(ph_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_z_x_xxxxz = cbuffer.data(ph_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_z_x_xxxyy = cbuffer.data(ph_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_z_x_xxxyz = cbuffer.data(ph_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_z_x_xxxzz = cbuffer.data(ph_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_z_x_xxyyy = cbuffer.data(ph_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_z_x_xxyyz = cbuffer.data(ph_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_z_x_xxyzz = cbuffer.data(ph_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_z_x_xxzzz = cbuffer.data(ph_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_z_x_xyyyy = cbuffer.data(ph_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_z_x_xyyyz = cbuffer.data(ph_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_z_x_xyyzz = cbuffer.data(ph_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_z_x_xyzzz = cbuffer.data(ph_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_z_x_xzzzz = cbuffer.data(ph_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_z_x_yyyyy = cbuffer.data(ph_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_z_x_yyyyz = cbuffer.data(ph_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_z_x_yyyzz = cbuffer.data(ph_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_z_x_yyzzz = cbuffer.data(ph_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_z_x_yzzzz = cbuffer.data(ph_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_z_x_zzzzz = cbuffer.data(ph_geom_11_off + 146 * ccomps * dcomps);

            auto g_x_z_y_xxxxx = cbuffer.data(ph_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_z_y_xxxxy = cbuffer.data(ph_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_z_y_xxxxz = cbuffer.data(ph_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_z_y_xxxyy = cbuffer.data(ph_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_z_y_xxxyz = cbuffer.data(ph_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_z_y_xxxzz = cbuffer.data(ph_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_z_y_xxyyy = cbuffer.data(ph_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_z_y_xxyyz = cbuffer.data(ph_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_z_y_xxyzz = cbuffer.data(ph_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_z_y_xxzzz = cbuffer.data(ph_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_z_y_xyyyy = cbuffer.data(ph_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_z_y_xyyyz = cbuffer.data(ph_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_z_y_xyyzz = cbuffer.data(ph_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_z_y_xyzzz = cbuffer.data(ph_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_z_y_xzzzz = cbuffer.data(ph_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_z_y_yyyyy = cbuffer.data(ph_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_z_y_yyyyz = cbuffer.data(ph_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_z_y_yyyzz = cbuffer.data(ph_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_z_y_yyzzz = cbuffer.data(ph_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_z_y_yzzzz = cbuffer.data(ph_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_z_y_zzzzz = cbuffer.data(ph_geom_11_off + 167 * ccomps * dcomps);

            auto g_x_z_z_xxxxx = cbuffer.data(ph_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_z_z_xxxxy = cbuffer.data(ph_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_z_z_xxxxz = cbuffer.data(ph_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_z_z_xxxyy = cbuffer.data(ph_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_z_z_xxxyz = cbuffer.data(ph_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_z_z_xxxzz = cbuffer.data(ph_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_z_z_xxyyy = cbuffer.data(ph_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_z_z_xxyyz = cbuffer.data(ph_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_z_z_xxyzz = cbuffer.data(ph_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_z_z_xxzzz = cbuffer.data(ph_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_z_z_xyyyy = cbuffer.data(ph_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_z_z_xyyyz = cbuffer.data(ph_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_z_z_xyyzz = cbuffer.data(ph_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_z_z_xyzzz = cbuffer.data(ph_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_z_z_xzzzz = cbuffer.data(ph_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_z_z_yyyyy = cbuffer.data(ph_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_z_z_yyyyz = cbuffer.data(ph_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_z_z_yyyzz = cbuffer.data(ph_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_z_z_yyzzz = cbuffer.data(ph_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_z_z_yzzzz = cbuffer.data(ph_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_z_z_zzzzz = cbuffer.data(ph_geom_11_off + 188 * ccomps * dcomps);

            auto g_y_x_x_xxxxx = cbuffer.data(ph_geom_11_off + 189 * ccomps * dcomps);

            auto g_y_x_x_xxxxy = cbuffer.data(ph_geom_11_off + 190 * ccomps * dcomps);

            auto g_y_x_x_xxxxz = cbuffer.data(ph_geom_11_off + 191 * ccomps * dcomps);

            auto g_y_x_x_xxxyy = cbuffer.data(ph_geom_11_off + 192 * ccomps * dcomps);

            auto g_y_x_x_xxxyz = cbuffer.data(ph_geom_11_off + 193 * ccomps * dcomps);

            auto g_y_x_x_xxxzz = cbuffer.data(ph_geom_11_off + 194 * ccomps * dcomps);

            auto g_y_x_x_xxyyy = cbuffer.data(ph_geom_11_off + 195 * ccomps * dcomps);

            auto g_y_x_x_xxyyz = cbuffer.data(ph_geom_11_off + 196 * ccomps * dcomps);

            auto g_y_x_x_xxyzz = cbuffer.data(ph_geom_11_off + 197 * ccomps * dcomps);

            auto g_y_x_x_xxzzz = cbuffer.data(ph_geom_11_off + 198 * ccomps * dcomps);

            auto g_y_x_x_xyyyy = cbuffer.data(ph_geom_11_off + 199 * ccomps * dcomps);

            auto g_y_x_x_xyyyz = cbuffer.data(ph_geom_11_off + 200 * ccomps * dcomps);

            auto g_y_x_x_xyyzz = cbuffer.data(ph_geom_11_off + 201 * ccomps * dcomps);

            auto g_y_x_x_xyzzz = cbuffer.data(ph_geom_11_off + 202 * ccomps * dcomps);

            auto g_y_x_x_xzzzz = cbuffer.data(ph_geom_11_off + 203 * ccomps * dcomps);

            auto g_y_x_x_yyyyy = cbuffer.data(ph_geom_11_off + 204 * ccomps * dcomps);

            auto g_y_x_x_yyyyz = cbuffer.data(ph_geom_11_off + 205 * ccomps * dcomps);

            auto g_y_x_x_yyyzz = cbuffer.data(ph_geom_11_off + 206 * ccomps * dcomps);

            auto g_y_x_x_yyzzz = cbuffer.data(ph_geom_11_off + 207 * ccomps * dcomps);

            auto g_y_x_x_yzzzz = cbuffer.data(ph_geom_11_off + 208 * ccomps * dcomps);

            auto g_y_x_x_zzzzz = cbuffer.data(ph_geom_11_off + 209 * ccomps * dcomps);

            auto g_y_x_y_xxxxx = cbuffer.data(ph_geom_11_off + 210 * ccomps * dcomps);

            auto g_y_x_y_xxxxy = cbuffer.data(ph_geom_11_off + 211 * ccomps * dcomps);

            auto g_y_x_y_xxxxz = cbuffer.data(ph_geom_11_off + 212 * ccomps * dcomps);

            auto g_y_x_y_xxxyy = cbuffer.data(ph_geom_11_off + 213 * ccomps * dcomps);

            auto g_y_x_y_xxxyz = cbuffer.data(ph_geom_11_off + 214 * ccomps * dcomps);

            auto g_y_x_y_xxxzz = cbuffer.data(ph_geom_11_off + 215 * ccomps * dcomps);

            auto g_y_x_y_xxyyy = cbuffer.data(ph_geom_11_off + 216 * ccomps * dcomps);

            auto g_y_x_y_xxyyz = cbuffer.data(ph_geom_11_off + 217 * ccomps * dcomps);

            auto g_y_x_y_xxyzz = cbuffer.data(ph_geom_11_off + 218 * ccomps * dcomps);

            auto g_y_x_y_xxzzz = cbuffer.data(ph_geom_11_off + 219 * ccomps * dcomps);

            auto g_y_x_y_xyyyy = cbuffer.data(ph_geom_11_off + 220 * ccomps * dcomps);

            auto g_y_x_y_xyyyz = cbuffer.data(ph_geom_11_off + 221 * ccomps * dcomps);

            auto g_y_x_y_xyyzz = cbuffer.data(ph_geom_11_off + 222 * ccomps * dcomps);

            auto g_y_x_y_xyzzz = cbuffer.data(ph_geom_11_off + 223 * ccomps * dcomps);

            auto g_y_x_y_xzzzz = cbuffer.data(ph_geom_11_off + 224 * ccomps * dcomps);

            auto g_y_x_y_yyyyy = cbuffer.data(ph_geom_11_off + 225 * ccomps * dcomps);

            auto g_y_x_y_yyyyz = cbuffer.data(ph_geom_11_off + 226 * ccomps * dcomps);

            auto g_y_x_y_yyyzz = cbuffer.data(ph_geom_11_off + 227 * ccomps * dcomps);

            auto g_y_x_y_yyzzz = cbuffer.data(ph_geom_11_off + 228 * ccomps * dcomps);

            auto g_y_x_y_yzzzz = cbuffer.data(ph_geom_11_off + 229 * ccomps * dcomps);

            auto g_y_x_y_zzzzz = cbuffer.data(ph_geom_11_off + 230 * ccomps * dcomps);

            auto g_y_x_z_xxxxx = cbuffer.data(ph_geom_11_off + 231 * ccomps * dcomps);

            auto g_y_x_z_xxxxy = cbuffer.data(ph_geom_11_off + 232 * ccomps * dcomps);

            auto g_y_x_z_xxxxz = cbuffer.data(ph_geom_11_off + 233 * ccomps * dcomps);

            auto g_y_x_z_xxxyy = cbuffer.data(ph_geom_11_off + 234 * ccomps * dcomps);

            auto g_y_x_z_xxxyz = cbuffer.data(ph_geom_11_off + 235 * ccomps * dcomps);

            auto g_y_x_z_xxxzz = cbuffer.data(ph_geom_11_off + 236 * ccomps * dcomps);

            auto g_y_x_z_xxyyy = cbuffer.data(ph_geom_11_off + 237 * ccomps * dcomps);

            auto g_y_x_z_xxyyz = cbuffer.data(ph_geom_11_off + 238 * ccomps * dcomps);

            auto g_y_x_z_xxyzz = cbuffer.data(ph_geom_11_off + 239 * ccomps * dcomps);

            auto g_y_x_z_xxzzz = cbuffer.data(ph_geom_11_off + 240 * ccomps * dcomps);

            auto g_y_x_z_xyyyy = cbuffer.data(ph_geom_11_off + 241 * ccomps * dcomps);

            auto g_y_x_z_xyyyz = cbuffer.data(ph_geom_11_off + 242 * ccomps * dcomps);

            auto g_y_x_z_xyyzz = cbuffer.data(ph_geom_11_off + 243 * ccomps * dcomps);

            auto g_y_x_z_xyzzz = cbuffer.data(ph_geom_11_off + 244 * ccomps * dcomps);

            auto g_y_x_z_xzzzz = cbuffer.data(ph_geom_11_off + 245 * ccomps * dcomps);

            auto g_y_x_z_yyyyy = cbuffer.data(ph_geom_11_off + 246 * ccomps * dcomps);

            auto g_y_x_z_yyyyz = cbuffer.data(ph_geom_11_off + 247 * ccomps * dcomps);

            auto g_y_x_z_yyyzz = cbuffer.data(ph_geom_11_off + 248 * ccomps * dcomps);

            auto g_y_x_z_yyzzz = cbuffer.data(ph_geom_11_off + 249 * ccomps * dcomps);

            auto g_y_x_z_yzzzz = cbuffer.data(ph_geom_11_off + 250 * ccomps * dcomps);

            auto g_y_x_z_zzzzz = cbuffer.data(ph_geom_11_off + 251 * ccomps * dcomps);

            auto g_y_y_x_xxxxx = cbuffer.data(ph_geom_11_off + 252 * ccomps * dcomps);

            auto g_y_y_x_xxxxy = cbuffer.data(ph_geom_11_off + 253 * ccomps * dcomps);

            auto g_y_y_x_xxxxz = cbuffer.data(ph_geom_11_off + 254 * ccomps * dcomps);

            auto g_y_y_x_xxxyy = cbuffer.data(ph_geom_11_off + 255 * ccomps * dcomps);

            auto g_y_y_x_xxxyz = cbuffer.data(ph_geom_11_off + 256 * ccomps * dcomps);

            auto g_y_y_x_xxxzz = cbuffer.data(ph_geom_11_off + 257 * ccomps * dcomps);

            auto g_y_y_x_xxyyy = cbuffer.data(ph_geom_11_off + 258 * ccomps * dcomps);

            auto g_y_y_x_xxyyz = cbuffer.data(ph_geom_11_off + 259 * ccomps * dcomps);

            auto g_y_y_x_xxyzz = cbuffer.data(ph_geom_11_off + 260 * ccomps * dcomps);

            auto g_y_y_x_xxzzz = cbuffer.data(ph_geom_11_off + 261 * ccomps * dcomps);

            auto g_y_y_x_xyyyy = cbuffer.data(ph_geom_11_off + 262 * ccomps * dcomps);

            auto g_y_y_x_xyyyz = cbuffer.data(ph_geom_11_off + 263 * ccomps * dcomps);

            auto g_y_y_x_xyyzz = cbuffer.data(ph_geom_11_off + 264 * ccomps * dcomps);

            auto g_y_y_x_xyzzz = cbuffer.data(ph_geom_11_off + 265 * ccomps * dcomps);

            auto g_y_y_x_xzzzz = cbuffer.data(ph_geom_11_off + 266 * ccomps * dcomps);

            auto g_y_y_x_yyyyy = cbuffer.data(ph_geom_11_off + 267 * ccomps * dcomps);

            auto g_y_y_x_yyyyz = cbuffer.data(ph_geom_11_off + 268 * ccomps * dcomps);

            auto g_y_y_x_yyyzz = cbuffer.data(ph_geom_11_off + 269 * ccomps * dcomps);

            auto g_y_y_x_yyzzz = cbuffer.data(ph_geom_11_off + 270 * ccomps * dcomps);

            auto g_y_y_x_yzzzz = cbuffer.data(ph_geom_11_off + 271 * ccomps * dcomps);

            auto g_y_y_x_zzzzz = cbuffer.data(ph_geom_11_off + 272 * ccomps * dcomps);

            auto g_y_y_y_xxxxx = cbuffer.data(ph_geom_11_off + 273 * ccomps * dcomps);

            auto g_y_y_y_xxxxy = cbuffer.data(ph_geom_11_off + 274 * ccomps * dcomps);

            auto g_y_y_y_xxxxz = cbuffer.data(ph_geom_11_off + 275 * ccomps * dcomps);

            auto g_y_y_y_xxxyy = cbuffer.data(ph_geom_11_off + 276 * ccomps * dcomps);

            auto g_y_y_y_xxxyz = cbuffer.data(ph_geom_11_off + 277 * ccomps * dcomps);

            auto g_y_y_y_xxxzz = cbuffer.data(ph_geom_11_off + 278 * ccomps * dcomps);

            auto g_y_y_y_xxyyy = cbuffer.data(ph_geom_11_off + 279 * ccomps * dcomps);

            auto g_y_y_y_xxyyz = cbuffer.data(ph_geom_11_off + 280 * ccomps * dcomps);

            auto g_y_y_y_xxyzz = cbuffer.data(ph_geom_11_off + 281 * ccomps * dcomps);

            auto g_y_y_y_xxzzz = cbuffer.data(ph_geom_11_off + 282 * ccomps * dcomps);

            auto g_y_y_y_xyyyy = cbuffer.data(ph_geom_11_off + 283 * ccomps * dcomps);

            auto g_y_y_y_xyyyz = cbuffer.data(ph_geom_11_off + 284 * ccomps * dcomps);

            auto g_y_y_y_xyyzz = cbuffer.data(ph_geom_11_off + 285 * ccomps * dcomps);

            auto g_y_y_y_xyzzz = cbuffer.data(ph_geom_11_off + 286 * ccomps * dcomps);

            auto g_y_y_y_xzzzz = cbuffer.data(ph_geom_11_off + 287 * ccomps * dcomps);

            auto g_y_y_y_yyyyy = cbuffer.data(ph_geom_11_off + 288 * ccomps * dcomps);

            auto g_y_y_y_yyyyz = cbuffer.data(ph_geom_11_off + 289 * ccomps * dcomps);

            auto g_y_y_y_yyyzz = cbuffer.data(ph_geom_11_off + 290 * ccomps * dcomps);

            auto g_y_y_y_yyzzz = cbuffer.data(ph_geom_11_off + 291 * ccomps * dcomps);

            auto g_y_y_y_yzzzz = cbuffer.data(ph_geom_11_off + 292 * ccomps * dcomps);

            auto g_y_y_y_zzzzz = cbuffer.data(ph_geom_11_off + 293 * ccomps * dcomps);

            auto g_y_y_z_xxxxx = cbuffer.data(ph_geom_11_off + 294 * ccomps * dcomps);

            auto g_y_y_z_xxxxy = cbuffer.data(ph_geom_11_off + 295 * ccomps * dcomps);

            auto g_y_y_z_xxxxz = cbuffer.data(ph_geom_11_off + 296 * ccomps * dcomps);

            auto g_y_y_z_xxxyy = cbuffer.data(ph_geom_11_off + 297 * ccomps * dcomps);

            auto g_y_y_z_xxxyz = cbuffer.data(ph_geom_11_off + 298 * ccomps * dcomps);

            auto g_y_y_z_xxxzz = cbuffer.data(ph_geom_11_off + 299 * ccomps * dcomps);

            auto g_y_y_z_xxyyy = cbuffer.data(ph_geom_11_off + 300 * ccomps * dcomps);

            auto g_y_y_z_xxyyz = cbuffer.data(ph_geom_11_off + 301 * ccomps * dcomps);

            auto g_y_y_z_xxyzz = cbuffer.data(ph_geom_11_off + 302 * ccomps * dcomps);

            auto g_y_y_z_xxzzz = cbuffer.data(ph_geom_11_off + 303 * ccomps * dcomps);

            auto g_y_y_z_xyyyy = cbuffer.data(ph_geom_11_off + 304 * ccomps * dcomps);

            auto g_y_y_z_xyyyz = cbuffer.data(ph_geom_11_off + 305 * ccomps * dcomps);

            auto g_y_y_z_xyyzz = cbuffer.data(ph_geom_11_off + 306 * ccomps * dcomps);

            auto g_y_y_z_xyzzz = cbuffer.data(ph_geom_11_off + 307 * ccomps * dcomps);

            auto g_y_y_z_xzzzz = cbuffer.data(ph_geom_11_off + 308 * ccomps * dcomps);

            auto g_y_y_z_yyyyy = cbuffer.data(ph_geom_11_off + 309 * ccomps * dcomps);

            auto g_y_y_z_yyyyz = cbuffer.data(ph_geom_11_off + 310 * ccomps * dcomps);

            auto g_y_y_z_yyyzz = cbuffer.data(ph_geom_11_off + 311 * ccomps * dcomps);

            auto g_y_y_z_yyzzz = cbuffer.data(ph_geom_11_off + 312 * ccomps * dcomps);

            auto g_y_y_z_yzzzz = cbuffer.data(ph_geom_11_off + 313 * ccomps * dcomps);

            auto g_y_y_z_zzzzz = cbuffer.data(ph_geom_11_off + 314 * ccomps * dcomps);

            auto g_y_z_x_xxxxx = cbuffer.data(ph_geom_11_off + 315 * ccomps * dcomps);

            auto g_y_z_x_xxxxy = cbuffer.data(ph_geom_11_off + 316 * ccomps * dcomps);

            auto g_y_z_x_xxxxz = cbuffer.data(ph_geom_11_off + 317 * ccomps * dcomps);

            auto g_y_z_x_xxxyy = cbuffer.data(ph_geom_11_off + 318 * ccomps * dcomps);

            auto g_y_z_x_xxxyz = cbuffer.data(ph_geom_11_off + 319 * ccomps * dcomps);

            auto g_y_z_x_xxxzz = cbuffer.data(ph_geom_11_off + 320 * ccomps * dcomps);

            auto g_y_z_x_xxyyy = cbuffer.data(ph_geom_11_off + 321 * ccomps * dcomps);

            auto g_y_z_x_xxyyz = cbuffer.data(ph_geom_11_off + 322 * ccomps * dcomps);

            auto g_y_z_x_xxyzz = cbuffer.data(ph_geom_11_off + 323 * ccomps * dcomps);

            auto g_y_z_x_xxzzz = cbuffer.data(ph_geom_11_off + 324 * ccomps * dcomps);

            auto g_y_z_x_xyyyy = cbuffer.data(ph_geom_11_off + 325 * ccomps * dcomps);

            auto g_y_z_x_xyyyz = cbuffer.data(ph_geom_11_off + 326 * ccomps * dcomps);

            auto g_y_z_x_xyyzz = cbuffer.data(ph_geom_11_off + 327 * ccomps * dcomps);

            auto g_y_z_x_xyzzz = cbuffer.data(ph_geom_11_off + 328 * ccomps * dcomps);

            auto g_y_z_x_xzzzz = cbuffer.data(ph_geom_11_off + 329 * ccomps * dcomps);

            auto g_y_z_x_yyyyy = cbuffer.data(ph_geom_11_off + 330 * ccomps * dcomps);

            auto g_y_z_x_yyyyz = cbuffer.data(ph_geom_11_off + 331 * ccomps * dcomps);

            auto g_y_z_x_yyyzz = cbuffer.data(ph_geom_11_off + 332 * ccomps * dcomps);

            auto g_y_z_x_yyzzz = cbuffer.data(ph_geom_11_off + 333 * ccomps * dcomps);

            auto g_y_z_x_yzzzz = cbuffer.data(ph_geom_11_off + 334 * ccomps * dcomps);

            auto g_y_z_x_zzzzz = cbuffer.data(ph_geom_11_off + 335 * ccomps * dcomps);

            auto g_y_z_y_xxxxx = cbuffer.data(ph_geom_11_off + 336 * ccomps * dcomps);

            auto g_y_z_y_xxxxy = cbuffer.data(ph_geom_11_off + 337 * ccomps * dcomps);

            auto g_y_z_y_xxxxz = cbuffer.data(ph_geom_11_off + 338 * ccomps * dcomps);

            auto g_y_z_y_xxxyy = cbuffer.data(ph_geom_11_off + 339 * ccomps * dcomps);

            auto g_y_z_y_xxxyz = cbuffer.data(ph_geom_11_off + 340 * ccomps * dcomps);

            auto g_y_z_y_xxxzz = cbuffer.data(ph_geom_11_off + 341 * ccomps * dcomps);

            auto g_y_z_y_xxyyy = cbuffer.data(ph_geom_11_off + 342 * ccomps * dcomps);

            auto g_y_z_y_xxyyz = cbuffer.data(ph_geom_11_off + 343 * ccomps * dcomps);

            auto g_y_z_y_xxyzz = cbuffer.data(ph_geom_11_off + 344 * ccomps * dcomps);

            auto g_y_z_y_xxzzz = cbuffer.data(ph_geom_11_off + 345 * ccomps * dcomps);

            auto g_y_z_y_xyyyy = cbuffer.data(ph_geom_11_off + 346 * ccomps * dcomps);

            auto g_y_z_y_xyyyz = cbuffer.data(ph_geom_11_off + 347 * ccomps * dcomps);

            auto g_y_z_y_xyyzz = cbuffer.data(ph_geom_11_off + 348 * ccomps * dcomps);

            auto g_y_z_y_xyzzz = cbuffer.data(ph_geom_11_off + 349 * ccomps * dcomps);

            auto g_y_z_y_xzzzz = cbuffer.data(ph_geom_11_off + 350 * ccomps * dcomps);

            auto g_y_z_y_yyyyy = cbuffer.data(ph_geom_11_off + 351 * ccomps * dcomps);

            auto g_y_z_y_yyyyz = cbuffer.data(ph_geom_11_off + 352 * ccomps * dcomps);

            auto g_y_z_y_yyyzz = cbuffer.data(ph_geom_11_off + 353 * ccomps * dcomps);

            auto g_y_z_y_yyzzz = cbuffer.data(ph_geom_11_off + 354 * ccomps * dcomps);

            auto g_y_z_y_yzzzz = cbuffer.data(ph_geom_11_off + 355 * ccomps * dcomps);

            auto g_y_z_y_zzzzz = cbuffer.data(ph_geom_11_off + 356 * ccomps * dcomps);

            auto g_y_z_z_xxxxx = cbuffer.data(ph_geom_11_off + 357 * ccomps * dcomps);

            auto g_y_z_z_xxxxy = cbuffer.data(ph_geom_11_off + 358 * ccomps * dcomps);

            auto g_y_z_z_xxxxz = cbuffer.data(ph_geom_11_off + 359 * ccomps * dcomps);

            auto g_y_z_z_xxxyy = cbuffer.data(ph_geom_11_off + 360 * ccomps * dcomps);

            auto g_y_z_z_xxxyz = cbuffer.data(ph_geom_11_off + 361 * ccomps * dcomps);

            auto g_y_z_z_xxxzz = cbuffer.data(ph_geom_11_off + 362 * ccomps * dcomps);

            auto g_y_z_z_xxyyy = cbuffer.data(ph_geom_11_off + 363 * ccomps * dcomps);

            auto g_y_z_z_xxyyz = cbuffer.data(ph_geom_11_off + 364 * ccomps * dcomps);

            auto g_y_z_z_xxyzz = cbuffer.data(ph_geom_11_off + 365 * ccomps * dcomps);

            auto g_y_z_z_xxzzz = cbuffer.data(ph_geom_11_off + 366 * ccomps * dcomps);

            auto g_y_z_z_xyyyy = cbuffer.data(ph_geom_11_off + 367 * ccomps * dcomps);

            auto g_y_z_z_xyyyz = cbuffer.data(ph_geom_11_off + 368 * ccomps * dcomps);

            auto g_y_z_z_xyyzz = cbuffer.data(ph_geom_11_off + 369 * ccomps * dcomps);

            auto g_y_z_z_xyzzz = cbuffer.data(ph_geom_11_off + 370 * ccomps * dcomps);

            auto g_y_z_z_xzzzz = cbuffer.data(ph_geom_11_off + 371 * ccomps * dcomps);

            auto g_y_z_z_yyyyy = cbuffer.data(ph_geom_11_off + 372 * ccomps * dcomps);

            auto g_y_z_z_yyyyz = cbuffer.data(ph_geom_11_off + 373 * ccomps * dcomps);

            auto g_y_z_z_yyyzz = cbuffer.data(ph_geom_11_off + 374 * ccomps * dcomps);

            auto g_y_z_z_yyzzz = cbuffer.data(ph_geom_11_off + 375 * ccomps * dcomps);

            auto g_y_z_z_yzzzz = cbuffer.data(ph_geom_11_off + 376 * ccomps * dcomps);

            auto g_y_z_z_zzzzz = cbuffer.data(ph_geom_11_off + 377 * ccomps * dcomps);

            auto g_z_x_x_xxxxx = cbuffer.data(ph_geom_11_off + 378 * ccomps * dcomps);

            auto g_z_x_x_xxxxy = cbuffer.data(ph_geom_11_off + 379 * ccomps * dcomps);

            auto g_z_x_x_xxxxz = cbuffer.data(ph_geom_11_off + 380 * ccomps * dcomps);

            auto g_z_x_x_xxxyy = cbuffer.data(ph_geom_11_off + 381 * ccomps * dcomps);

            auto g_z_x_x_xxxyz = cbuffer.data(ph_geom_11_off + 382 * ccomps * dcomps);

            auto g_z_x_x_xxxzz = cbuffer.data(ph_geom_11_off + 383 * ccomps * dcomps);

            auto g_z_x_x_xxyyy = cbuffer.data(ph_geom_11_off + 384 * ccomps * dcomps);

            auto g_z_x_x_xxyyz = cbuffer.data(ph_geom_11_off + 385 * ccomps * dcomps);

            auto g_z_x_x_xxyzz = cbuffer.data(ph_geom_11_off + 386 * ccomps * dcomps);

            auto g_z_x_x_xxzzz = cbuffer.data(ph_geom_11_off + 387 * ccomps * dcomps);

            auto g_z_x_x_xyyyy = cbuffer.data(ph_geom_11_off + 388 * ccomps * dcomps);

            auto g_z_x_x_xyyyz = cbuffer.data(ph_geom_11_off + 389 * ccomps * dcomps);

            auto g_z_x_x_xyyzz = cbuffer.data(ph_geom_11_off + 390 * ccomps * dcomps);

            auto g_z_x_x_xyzzz = cbuffer.data(ph_geom_11_off + 391 * ccomps * dcomps);

            auto g_z_x_x_xzzzz = cbuffer.data(ph_geom_11_off + 392 * ccomps * dcomps);

            auto g_z_x_x_yyyyy = cbuffer.data(ph_geom_11_off + 393 * ccomps * dcomps);

            auto g_z_x_x_yyyyz = cbuffer.data(ph_geom_11_off + 394 * ccomps * dcomps);

            auto g_z_x_x_yyyzz = cbuffer.data(ph_geom_11_off + 395 * ccomps * dcomps);

            auto g_z_x_x_yyzzz = cbuffer.data(ph_geom_11_off + 396 * ccomps * dcomps);

            auto g_z_x_x_yzzzz = cbuffer.data(ph_geom_11_off + 397 * ccomps * dcomps);

            auto g_z_x_x_zzzzz = cbuffer.data(ph_geom_11_off + 398 * ccomps * dcomps);

            auto g_z_x_y_xxxxx = cbuffer.data(ph_geom_11_off + 399 * ccomps * dcomps);

            auto g_z_x_y_xxxxy = cbuffer.data(ph_geom_11_off + 400 * ccomps * dcomps);

            auto g_z_x_y_xxxxz = cbuffer.data(ph_geom_11_off + 401 * ccomps * dcomps);

            auto g_z_x_y_xxxyy = cbuffer.data(ph_geom_11_off + 402 * ccomps * dcomps);

            auto g_z_x_y_xxxyz = cbuffer.data(ph_geom_11_off + 403 * ccomps * dcomps);

            auto g_z_x_y_xxxzz = cbuffer.data(ph_geom_11_off + 404 * ccomps * dcomps);

            auto g_z_x_y_xxyyy = cbuffer.data(ph_geom_11_off + 405 * ccomps * dcomps);

            auto g_z_x_y_xxyyz = cbuffer.data(ph_geom_11_off + 406 * ccomps * dcomps);

            auto g_z_x_y_xxyzz = cbuffer.data(ph_geom_11_off + 407 * ccomps * dcomps);

            auto g_z_x_y_xxzzz = cbuffer.data(ph_geom_11_off + 408 * ccomps * dcomps);

            auto g_z_x_y_xyyyy = cbuffer.data(ph_geom_11_off + 409 * ccomps * dcomps);

            auto g_z_x_y_xyyyz = cbuffer.data(ph_geom_11_off + 410 * ccomps * dcomps);

            auto g_z_x_y_xyyzz = cbuffer.data(ph_geom_11_off + 411 * ccomps * dcomps);

            auto g_z_x_y_xyzzz = cbuffer.data(ph_geom_11_off + 412 * ccomps * dcomps);

            auto g_z_x_y_xzzzz = cbuffer.data(ph_geom_11_off + 413 * ccomps * dcomps);

            auto g_z_x_y_yyyyy = cbuffer.data(ph_geom_11_off + 414 * ccomps * dcomps);

            auto g_z_x_y_yyyyz = cbuffer.data(ph_geom_11_off + 415 * ccomps * dcomps);

            auto g_z_x_y_yyyzz = cbuffer.data(ph_geom_11_off + 416 * ccomps * dcomps);

            auto g_z_x_y_yyzzz = cbuffer.data(ph_geom_11_off + 417 * ccomps * dcomps);

            auto g_z_x_y_yzzzz = cbuffer.data(ph_geom_11_off + 418 * ccomps * dcomps);

            auto g_z_x_y_zzzzz = cbuffer.data(ph_geom_11_off + 419 * ccomps * dcomps);

            auto g_z_x_z_xxxxx = cbuffer.data(ph_geom_11_off + 420 * ccomps * dcomps);

            auto g_z_x_z_xxxxy = cbuffer.data(ph_geom_11_off + 421 * ccomps * dcomps);

            auto g_z_x_z_xxxxz = cbuffer.data(ph_geom_11_off + 422 * ccomps * dcomps);

            auto g_z_x_z_xxxyy = cbuffer.data(ph_geom_11_off + 423 * ccomps * dcomps);

            auto g_z_x_z_xxxyz = cbuffer.data(ph_geom_11_off + 424 * ccomps * dcomps);

            auto g_z_x_z_xxxzz = cbuffer.data(ph_geom_11_off + 425 * ccomps * dcomps);

            auto g_z_x_z_xxyyy = cbuffer.data(ph_geom_11_off + 426 * ccomps * dcomps);

            auto g_z_x_z_xxyyz = cbuffer.data(ph_geom_11_off + 427 * ccomps * dcomps);

            auto g_z_x_z_xxyzz = cbuffer.data(ph_geom_11_off + 428 * ccomps * dcomps);

            auto g_z_x_z_xxzzz = cbuffer.data(ph_geom_11_off + 429 * ccomps * dcomps);

            auto g_z_x_z_xyyyy = cbuffer.data(ph_geom_11_off + 430 * ccomps * dcomps);

            auto g_z_x_z_xyyyz = cbuffer.data(ph_geom_11_off + 431 * ccomps * dcomps);

            auto g_z_x_z_xyyzz = cbuffer.data(ph_geom_11_off + 432 * ccomps * dcomps);

            auto g_z_x_z_xyzzz = cbuffer.data(ph_geom_11_off + 433 * ccomps * dcomps);

            auto g_z_x_z_xzzzz = cbuffer.data(ph_geom_11_off + 434 * ccomps * dcomps);

            auto g_z_x_z_yyyyy = cbuffer.data(ph_geom_11_off + 435 * ccomps * dcomps);

            auto g_z_x_z_yyyyz = cbuffer.data(ph_geom_11_off + 436 * ccomps * dcomps);

            auto g_z_x_z_yyyzz = cbuffer.data(ph_geom_11_off + 437 * ccomps * dcomps);

            auto g_z_x_z_yyzzz = cbuffer.data(ph_geom_11_off + 438 * ccomps * dcomps);

            auto g_z_x_z_yzzzz = cbuffer.data(ph_geom_11_off + 439 * ccomps * dcomps);

            auto g_z_x_z_zzzzz = cbuffer.data(ph_geom_11_off + 440 * ccomps * dcomps);

            auto g_z_y_x_xxxxx = cbuffer.data(ph_geom_11_off + 441 * ccomps * dcomps);

            auto g_z_y_x_xxxxy = cbuffer.data(ph_geom_11_off + 442 * ccomps * dcomps);

            auto g_z_y_x_xxxxz = cbuffer.data(ph_geom_11_off + 443 * ccomps * dcomps);

            auto g_z_y_x_xxxyy = cbuffer.data(ph_geom_11_off + 444 * ccomps * dcomps);

            auto g_z_y_x_xxxyz = cbuffer.data(ph_geom_11_off + 445 * ccomps * dcomps);

            auto g_z_y_x_xxxzz = cbuffer.data(ph_geom_11_off + 446 * ccomps * dcomps);

            auto g_z_y_x_xxyyy = cbuffer.data(ph_geom_11_off + 447 * ccomps * dcomps);

            auto g_z_y_x_xxyyz = cbuffer.data(ph_geom_11_off + 448 * ccomps * dcomps);

            auto g_z_y_x_xxyzz = cbuffer.data(ph_geom_11_off + 449 * ccomps * dcomps);

            auto g_z_y_x_xxzzz = cbuffer.data(ph_geom_11_off + 450 * ccomps * dcomps);

            auto g_z_y_x_xyyyy = cbuffer.data(ph_geom_11_off + 451 * ccomps * dcomps);

            auto g_z_y_x_xyyyz = cbuffer.data(ph_geom_11_off + 452 * ccomps * dcomps);

            auto g_z_y_x_xyyzz = cbuffer.data(ph_geom_11_off + 453 * ccomps * dcomps);

            auto g_z_y_x_xyzzz = cbuffer.data(ph_geom_11_off + 454 * ccomps * dcomps);

            auto g_z_y_x_xzzzz = cbuffer.data(ph_geom_11_off + 455 * ccomps * dcomps);

            auto g_z_y_x_yyyyy = cbuffer.data(ph_geom_11_off + 456 * ccomps * dcomps);

            auto g_z_y_x_yyyyz = cbuffer.data(ph_geom_11_off + 457 * ccomps * dcomps);

            auto g_z_y_x_yyyzz = cbuffer.data(ph_geom_11_off + 458 * ccomps * dcomps);

            auto g_z_y_x_yyzzz = cbuffer.data(ph_geom_11_off + 459 * ccomps * dcomps);

            auto g_z_y_x_yzzzz = cbuffer.data(ph_geom_11_off + 460 * ccomps * dcomps);

            auto g_z_y_x_zzzzz = cbuffer.data(ph_geom_11_off + 461 * ccomps * dcomps);

            auto g_z_y_y_xxxxx = cbuffer.data(ph_geom_11_off + 462 * ccomps * dcomps);

            auto g_z_y_y_xxxxy = cbuffer.data(ph_geom_11_off + 463 * ccomps * dcomps);

            auto g_z_y_y_xxxxz = cbuffer.data(ph_geom_11_off + 464 * ccomps * dcomps);

            auto g_z_y_y_xxxyy = cbuffer.data(ph_geom_11_off + 465 * ccomps * dcomps);

            auto g_z_y_y_xxxyz = cbuffer.data(ph_geom_11_off + 466 * ccomps * dcomps);

            auto g_z_y_y_xxxzz = cbuffer.data(ph_geom_11_off + 467 * ccomps * dcomps);

            auto g_z_y_y_xxyyy = cbuffer.data(ph_geom_11_off + 468 * ccomps * dcomps);

            auto g_z_y_y_xxyyz = cbuffer.data(ph_geom_11_off + 469 * ccomps * dcomps);

            auto g_z_y_y_xxyzz = cbuffer.data(ph_geom_11_off + 470 * ccomps * dcomps);

            auto g_z_y_y_xxzzz = cbuffer.data(ph_geom_11_off + 471 * ccomps * dcomps);

            auto g_z_y_y_xyyyy = cbuffer.data(ph_geom_11_off + 472 * ccomps * dcomps);

            auto g_z_y_y_xyyyz = cbuffer.data(ph_geom_11_off + 473 * ccomps * dcomps);

            auto g_z_y_y_xyyzz = cbuffer.data(ph_geom_11_off + 474 * ccomps * dcomps);

            auto g_z_y_y_xyzzz = cbuffer.data(ph_geom_11_off + 475 * ccomps * dcomps);

            auto g_z_y_y_xzzzz = cbuffer.data(ph_geom_11_off + 476 * ccomps * dcomps);

            auto g_z_y_y_yyyyy = cbuffer.data(ph_geom_11_off + 477 * ccomps * dcomps);

            auto g_z_y_y_yyyyz = cbuffer.data(ph_geom_11_off + 478 * ccomps * dcomps);

            auto g_z_y_y_yyyzz = cbuffer.data(ph_geom_11_off + 479 * ccomps * dcomps);

            auto g_z_y_y_yyzzz = cbuffer.data(ph_geom_11_off + 480 * ccomps * dcomps);

            auto g_z_y_y_yzzzz = cbuffer.data(ph_geom_11_off + 481 * ccomps * dcomps);

            auto g_z_y_y_zzzzz = cbuffer.data(ph_geom_11_off + 482 * ccomps * dcomps);

            auto g_z_y_z_xxxxx = cbuffer.data(ph_geom_11_off + 483 * ccomps * dcomps);

            auto g_z_y_z_xxxxy = cbuffer.data(ph_geom_11_off + 484 * ccomps * dcomps);

            auto g_z_y_z_xxxxz = cbuffer.data(ph_geom_11_off + 485 * ccomps * dcomps);

            auto g_z_y_z_xxxyy = cbuffer.data(ph_geom_11_off + 486 * ccomps * dcomps);

            auto g_z_y_z_xxxyz = cbuffer.data(ph_geom_11_off + 487 * ccomps * dcomps);

            auto g_z_y_z_xxxzz = cbuffer.data(ph_geom_11_off + 488 * ccomps * dcomps);

            auto g_z_y_z_xxyyy = cbuffer.data(ph_geom_11_off + 489 * ccomps * dcomps);

            auto g_z_y_z_xxyyz = cbuffer.data(ph_geom_11_off + 490 * ccomps * dcomps);

            auto g_z_y_z_xxyzz = cbuffer.data(ph_geom_11_off + 491 * ccomps * dcomps);

            auto g_z_y_z_xxzzz = cbuffer.data(ph_geom_11_off + 492 * ccomps * dcomps);

            auto g_z_y_z_xyyyy = cbuffer.data(ph_geom_11_off + 493 * ccomps * dcomps);

            auto g_z_y_z_xyyyz = cbuffer.data(ph_geom_11_off + 494 * ccomps * dcomps);

            auto g_z_y_z_xyyzz = cbuffer.data(ph_geom_11_off + 495 * ccomps * dcomps);

            auto g_z_y_z_xyzzz = cbuffer.data(ph_geom_11_off + 496 * ccomps * dcomps);

            auto g_z_y_z_xzzzz = cbuffer.data(ph_geom_11_off + 497 * ccomps * dcomps);

            auto g_z_y_z_yyyyy = cbuffer.data(ph_geom_11_off + 498 * ccomps * dcomps);

            auto g_z_y_z_yyyyz = cbuffer.data(ph_geom_11_off + 499 * ccomps * dcomps);

            auto g_z_y_z_yyyzz = cbuffer.data(ph_geom_11_off + 500 * ccomps * dcomps);

            auto g_z_y_z_yyzzz = cbuffer.data(ph_geom_11_off + 501 * ccomps * dcomps);

            auto g_z_y_z_yzzzz = cbuffer.data(ph_geom_11_off + 502 * ccomps * dcomps);

            auto g_z_y_z_zzzzz = cbuffer.data(ph_geom_11_off + 503 * ccomps * dcomps);

            auto g_z_z_x_xxxxx = cbuffer.data(ph_geom_11_off + 504 * ccomps * dcomps);

            auto g_z_z_x_xxxxy = cbuffer.data(ph_geom_11_off + 505 * ccomps * dcomps);

            auto g_z_z_x_xxxxz = cbuffer.data(ph_geom_11_off + 506 * ccomps * dcomps);

            auto g_z_z_x_xxxyy = cbuffer.data(ph_geom_11_off + 507 * ccomps * dcomps);

            auto g_z_z_x_xxxyz = cbuffer.data(ph_geom_11_off + 508 * ccomps * dcomps);

            auto g_z_z_x_xxxzz = cbuffer.data(ph_geom_11_off + 509 * ccomps * dcomps);

            auto g_z_z_x_xxyyy = cbuffer.data(ph_geom_11_off + 510 * ccomps * dcomps);

            auto g_z_z_x_xxyyz = cbuffer.data(ph_geom_11_off + 511 * ccomps * dcomps);

            auto g_z_z_x_xxyzz = cbuffer.data(ph_geom_11_off + 512 * ccomps * dcomps);

            auto g_z_z_x_xxzzz = cbuffer.data(ph_geom_11_off + 513 * ccomps * dcomps);

            auto g_z_z_x_xyyyy = cbuffer.data(ph_geom_11_off + 514 * ccomps * dcomps);

            auto g_z_z_x_xyyyz = cbuffer.data(ph_geom_11_off + 515 * ccomps * dcomps);

            auto g_z_z_x_xyyzz = cbuffer.data(ph_geom_11_off + 516 * ccomps * dcomps);

            auto g_z_z_x_xyzzz = cbuffer.data(ph_geom_11_off + 517 * ccomps * dcomps);

            auto g_z_z_x_xzzzz = cbuffer.data(ph_geom_11_off + 518 * ccomps * dcomps);

            auto g_z_z_x_yyyyy = cbuffer.data(ph_geom_11_off + 519 * ccomps * dcomps);

            auto g_z_z_x_yyyyz = cbuffer.data(ph_geom_11_off + 520 * ccomps * dcomps);

            auto g_z_z_x_yyyzz = cbuffer.data(ph_geom_11_off + 521 * ccomps * dcomps);

            auto g_z_z_x_yyzzz = cbuffer.data(ph_geom_11_off + 522 * ccomps * dcomps);

            auto g_z_z_x_yzzzz = cbuffer.data(ph_geom_11_off + 523 * ccomps * dcomps);

            auto g_z_z_x_zzzzz = cbuffer.data(ph_geom_11_off + 524 * ccomps * dcomps);

            auto g_z_z_y_xxxxx = cbuffer.data(ph_geom_11_off + 525 * ccomps * dcomps);

            auto g_z_z_y_xxxxy = cbuffer.data(ph_geom_11_off + 526 * ccomps * dcomps);

            auto g_z_z_y_xxxxz = cbuffer.data(ph_geom_11_off + 527 * ccomps * dcomps);

            auto g_z_z_y_xxxyy = cbuffer.data(ph_geom_11_off + 528 * ccomps * dcomps);

            auto g_z_z_y_xxxyz = cbuffer.data(ph_geom_11_off + 529 * ccomps * dcomps);

            auto g_z_z_y_xxxzz = cbuffer.data(ph_geom_11_off + 530 * ccomps * dcomps);

            auto g_z_z_y_xxyyy = cbuffer.data(ph_geom_11_off + 531 * ccomps * dcomps);

            auto g_z_z_y_xxyyz = cbuffer.data(ph_geom_11_off + 532 * ccomps * dcomps);

            auto g_z_z_y_xxyzz = cbuffer.data(ph_geom_11_off + 533 * ccomps * dcomps);

            auto g_z_z_y_xxzzz = cbuffer.data(ph_geom_11_off + 534 * ccomps * dcomps);

            auto g_z_z_y_xyyyy = cbuffer.data(ph_geom_11_off + 535 * ccomps * dcomps);

            auto g_z_z_y_xyyyz = cbuffer.data(ph_geom_11_off + 536 * ccomps * dcomps);

            auto g_z_z_y_xyyzz = cbuffer.data(ph_geom_11_off + 537 * ccomps * dcomps);

            auto g_z_z_y_xyzzz = cbuffer.data(ph_geom_11_off + 538 * ccomps * dcomps);

            auto g_z_z_y_xzzzz = cbuffer.data(ph_geom_11_off + 539 * ccomps * dcomps);

            auto g_z_z_y_yyyyy = cbuffer.data(ph_geom_11_off + 540 * ccomps * dcomps);

            auto g_z_z_y_yyyyz = cbuffer.data(ph_geom_11_off + 541 * ccomps * dcomps);

            auto g_z_z_y_yyyzz = cbuffer.data(ph_geom_11_off + 542 * ccomps * dcomps);

            auto g_z_z_y_yyzzz = cbuffer.data(ph_geom_11_off + 543 * ccomps * dcomps);

            auto g_z_z_y_yzzzz = cbuffer.data(ph_geom_11_off + 544 * ccomps * dcomps);

            auto g_z_z_y_zzzzz = cbuffer.data(ph_geom_11_off + 545 * ccomps * dcomps);

            auto g_z_z_z_xxxxx = cbuffer.data(ph_geom_11_off + 546 * ccomps * dcomps);

            auto g_z_z_z_xxxxy = cbuffer.data(ph_geom_11_off + 547 * ccomps * dcomps);

            auto g_z_z_z_xxxxz = cbuffer.data(ph_geom_11_off + 548 * ccomps * dcomps);

            auto g_z_z_z_xxxyy = cbuffer.data(ph_geom_11_off + 549 * ccomps * dcomps);

            auto g_z_z_z_xxxyz = cbuffer.data(ph_geom_11_off + 550 * ccomps * dcomps);

            auto g_z_z_z_xxxzz = cbuffer.data(ph_geom_11_off + 551 * ccomps * dcomps);

            auto g_z_z_z_xxyyy = cbuffer.data(ph_geom_11_off + 552 * ccomps * dcomps);

            auto g_z_z_z_xxyyz = cbuffer.data(ph_geom_11_off + 553 * ccomps * dcomps);

            auto g_z_z_z_xxyzz = cbuffer.data(ph_geom_11_off + 554 * ccomps * dcomps);

            auto g_z_z_z_xxzzz = cbuffer.data(ph_geom_11_off + 555 * ccomps * dcomps);

            auto g_z_z_z_xyyyy = cbuffer.data(ph_geom_11_off + 556 * ccomps * dcomps);

            auto g_z_z_z_xyyyz = cbuffer.data(ph_geom_11_off + 557 * ccomps * dcomps);

            auto g_z_z_z_xyyzz = cbuffer.data(ph_geom_11_off + 558 * ccomps * dcomps);

            auto g_z_z_z_xyzzz = cbuffer.data(ph_geom_11_off + 559 * ccomps * dcomps);

            auto g_z_z_z_xzzzz = cbuffer.data(ph_geom_11_off + 560 * ccomps * dcomps);

            auto g_z_z_z_yyyyy = cbuffer.data(ph_geom_11_off + 561 * ccomps * dcomps);

            auto g_z_z_z_yyyyz = cbuffer.data(ph_geom_11_off + 562 * ccomps * dcomps);

            auto g_z_z_z_yyyzz = cbuffer.data(ph_geom_11_off + 563 * ccomps * dcomps);

            auto g_z_z_z_yyzzz = cbuffer.data(ph_geom_11_off + 564 * ccomps * dcomps);

            auto g_z_z_z_yzzzz = cbuffer.data(ph_geom_11_off + 565 * ccomps * dcomps);

            auto g_z_z_z_zzzzz = cbuffer.data(ph_geom_11_off + 566 * ccomps * dcomps);

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

            /// set up bra offset for contr_buffer_dhxx

            const auto dh_geom_11_off = idx_geom_11_dhxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_x_x_xx_xxxxx = cbuffer.data(dh_geom_11_off + 0 * ccomps * dcomps);

            auto g_x_x_xx_xxxxy = cbuffer.data(dh_geom_11_off + 1 * ccomps * dcomps);

            auto g_x_x_xx_xxxxz = cbuffer.data(dh_geom_11_off + 2 * ccomps * dcomps);

            auto g_x_x_xx_xxxyy = cbuffer.data(dh_geom_11_off + 3 * ccomps * dcomps);

            auto g_x_x_xx_xxxyz = cbuffer.data(dh_geom_11_off + 4 * ccomps * dcomps);

            auto g_x_x_xx_xxxzz = cbuffer.data(dh_geom_11_off + 5 * ccomps * dcomps);

            auto g_x_x_xx_xxyyy = cbuffer.data(dh_geom_11_off + 6 * ccomps * dcomps);

            auto g_x_x_xx_xxyyz = cbuffer.data(dh_geom_11_off + 7 * ccomps * dcomps);

            auto g_x_x_xx_xxyzz = cbuffer.data(dh_geom_11_off + 8 * ccomps * dcomps);

            auto g_x_x_xx_xxzzz = cbuffer.data(dh_geom_11_off + 9 * ccomps * dcomps);

            auto g_x_x_xx_xyyyy = cbuffer.data(dh_geom_11_off + 10 * ccomps * dcomps);

            auto g_x_x_xx_xyyyz = cbuffer.data(dh_geom_11_off + 11 * ccomps * dcomps);

            auto g_x_x_xx_xyyzz = cbuffer.data(dh_geom_11_off + 12 * ccomps * dcomps);

            auto g_x_x_xx_xyzzz = cbuffer.data(dh_geom_11_off + 13 * ccomps * dcomps);

            auto g_x_x_xx_xzzzz = cbuffer.data(dh_geom_11_off + 14 * ccomps * dcomps);

            auto g_x_x_xx_yyyyy = cbuffer.data(dh_geom_11_off + 15 * ccomps * dcomps);

            auto g_x_x_xx_yyyyz = cbuffer.data(dh_geom_11_off + 16 * ccomps * dcomps);

            auto g_x_x_xx_yyyzz = cbuffer.data(dh_geom_11_off + 17 * ccomps * dcomps);

            auto g_x_x_xx_yyzzz = cbuffer.data(dh_geom_11_off + 18 * ccomps * dcomps);

            auto g_x_x_xx_yzzzz = cbuffer.data(dh_geom_11_off + 19 * ccomps * dcomps);

            auto g_x_x_xx_zzzzz = cbuffer.data(dh_geom_11_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_x_xxxxx, g_0_x_x_xxxxy, g_0_x_x_xxxxz, g_0_x_x_xxxyy, g_0_x_x_xxxyz, g_0_x_x_xxxzz, g_0_x_x_xxyyy, g_0_x_x_xxyyz, g_0_x_x_xxyzz, g_0_x_x_xxzzz, g_0_x_x_xyyyy, g_0_x_x_xyyyz, g_0_x_x_xyyzz, g_0_x_x_xyzzz, g_0_x_x_xzzzz, g_0_x_x_yyyyy, g_0_x_x_yyyyz, g_0_x_x_yyyzz, g_0_x_x_yyzzz, g_0_x_x_yzzzz, g_0_x_x_zzzzz, g_x_0_x_xxxxx, g_x_0_x_xxxxy, g_x_0_x_xxxxz, g_x_0_x_xxxyy, g_x_0_x_xxxyz, g_x_0_x_xxxzz, g_x_0_x_xxyyy, g_x_0_x_xxyyz, g_x_0_x_xxyzz, g_x_0_x_xxzzz, g_x_0_x_xyyyy, g_x_0_x_xyyyz, g_x_0_x_xyyzz, g_x_0_x_xyzzz, g_x_0_x_xzzzz, g_x_0_x_yyyyy, g_x_0_x_yyyyz, g_x_0_x_yyyzz, g_x_0_x_yyzzz, g_x_0_x_yzzzz, g_x_0_x_zzzzz, g_x_x_x_xxxxx, g_x_x_x_xxxxxx, g_x_x_x_xxxxxy, g_x_x_x_xxxxxz, g_x_x_x_xxxxy, g_x_x_x_xxxxyy, g_x_x_x_xxxxyz, g_x_x_x_xxxxz, g_x_x_x_xxxxzz, g_x_x_x_xxxyy, g_x_x_x_xxxyyy, g_x_x_x_xxxyyz, g_x_x_x_xxxyz, g_x_x_x_xxxyzz, g_x_x_x_xxxzz, g_x_x_x_xxxzzz, g_x_x_x_xxyyy, g_x_x_x_xxyyyy, g_x_x_x_xxyyyz, g_x_x_x_xxyyz, g_x_x_x_xxyyzz, g_x_x_x_xxyzz, g_x_x_x_xxyzzz, g_x_x_x_xxzzz, g_x_x_x_xxzzzz, g_x_x_x_xyyyy, g_x_x_x_xyyyyy, g_x_x_x_xyyyyz, g_x_x_x_xyyyz, g_x_x_x_xyyyzz, g_x_x_x_xyyzz, g_x_x_x_xyyzzz, g_x_x_x_xyzzz, g_x_x_x_xyzzzz, g_x_x_x_xzzzz, g_x_x_x_xzzzzz, g_x_x_x_yyyyy, g_x_x_x_yyyyz, g_x_x_x_yyyzz, g_x_x_x_yyzzz, g_x_x_x_yzzzz, g_x_x_x_zzzzz, g_x_x_xx_xxxxx, g_x_x_xx_xxxxy, g_x_x_xx_xxxxz, g_x_x_xx_xxxyy, g_x_x_xx_xxxyz, g_x_x_xx_xxxzz, g_x_x_xx_xxyyy, g_x_x_xx_xxyyz, g_x_x_xx_xxyzz, g_x_x_xx_xxzzz, g_x_x_xx_xyyyy, g_x_x_xx_xyyyz, g_x_x_xx_xyyzz, g_x_x_xx_xyzzz, g_x_x_xx_xzzzz, g_x_x_xx_yyyyy, g_x_x_xx_yyyyz, g_x_x_xx_yyyzz, g_x_x_xx_yyzzz, g_x_x_xx_yzzzz, g_x_x_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xx_xxxxx[k] = -g_0_x_x_xxxxx[k] + g_x_0_x_xxxxx[k] - g_x_x_x_xxxxx[k] * ab_x + g_x_x_x_xxxxxx[k];

                g_x_x_xx_xxxxy[k] = -g_0_x_x_xxxxy[k] + g_x_0_x_xxxxy[k] - g_x_x_x_xxxxy[k] * ab_x + g_x_x_x_xxxxxy[k];

                g_x_x_xx_xxxxz[k] = -g_0_x_x_xxxxz[k] + g_x_0_x_xxxxz[k] - g_x_x_x_xxxxz[k] * ab_x + g_x_x_x_xxxxxz[k];

                g_x_x_xx_xxxyy[k] = -g_0_x_x_xxxyy[k] + g_x_0_x_xxxyy[k] - g_x_x_x_xxxyy[k] * ab_x + g_x_x_x_xxxxyy[k];

                g_x_x_xx_xxxyz[k] = -g_0_x_x_xxxyz[k] + g_x_0_x_xxxyz[k] - g_x_x_x_xxxyz[k] * ab_x + g_x_x_x_xxxxyz[k];

                g_x_x_xx_xxxzz[k] = -g_0_x_x_xxxzz[k] + g_x_0_x_xxxzz[k] - g_x_x_x_xxxzz[k] * ab_x + g_x_x_x_xxxxzz[k];

                g_x_x_xx_xxyyy[k] = -g_0_x_x_xxyyy[k] + g_x_0_x_xxyyy[k] - g_x_x_x_xxyyy[k] * ab_x + g_x_x_x_xxxyyy[k];

                g_x_x_xx_xxyyz[k] = -g_0_x_x_xxyyz[k] + g_x_0_x_xxyyz[k] - g_x_x_x_xxyyz[k] * ab_x + g_x_x_x_xxxyyz[k];

                g_x_x_xx_xxyzz[k] = -g_0_x_x_xxyzz[k] + g_x_0_x_xxyzz[k] - g_x_x_x_xxyzz[k] * ab_x + g_x_x_x_xxxyzz[k];

                g_x_x_xx_xxzzz[k] = -g_0_x_x_xxzzz[k] + g_x_0_x_xxzzz[k] - g_x_x_x_xxzzz[k] * ab_x + g_x_x_x_xxxzzz[k];

                g_x_x_xx_xyyyy[k] = -g_0_x_x_xyyyy[k] + g_x_0_x_xyyyy[k] - g_x_x_x_xyyyy[k] * ab_x + g_x_x_x_xxyyyy[k];

                g_x_x_xx_xyyyz[k] = -g_0_x_x_xyyyz[k] + g_x_0_x_xyyyz[k] - g_x_x_x_xyyyz[k] * ab_x + g_x_x_x_xxyyyz[k];

                g_x_x_xx_xyyzz[k] = -g_0_x_x_xyyzz[k] + g_x_0_x_xyyzz[k] - g_x_x_x_xyyzz[k] * ab_x + g_x_x_x_xxyyzz[k];

                g_x_x_xx_xyzzz[k] = -g_0_x_x_xyzzz[k] + g_x_0_x_xyzzz[k] - g_x_x_x_xyzzz[k] * ab_x + g_x_x_x_xxyzzz[k];

                g_x_x_xx_xzzzz[k] = -g_0_x_x_xzzzz[k] + g_x_0_x_xzzzz[k] - g_x_x_x_xzzzz[k] * ab_x + g_x_x_x_xxzzzz[k];

                g_x_x_xx_yyyyy[k] = -g_0_x_x_yyyyy[k] + g_x_0_x_yyyyy[k] - g_x_x_x_yyyyy[k] * ab_x + g_x_x_x_xyyyyy[k];

                g_x_x_xx_yyyyz[k] = -g_0_x_x_yyyyz[k] + g_x_0_x_yyyyz[k] - g_x_x_x_yyyyz[k] * ab_x + g_x_x_x_xyyyyz[k];

                g_x_x_xx_yyyzz[k] = -g_0_x_x_yyyzz[k] + g_x_0_x_yyyzz[k] - g_x_x_x_yyyzz[k] * ab_x + g_x_x_x_xyyyzz[k];

                g_x_x_xx_yyzzz[k] = -g_0_x_x_yyzzz[k] + g_x_0_x_yyzzz[k] - g_x_x_x_yyzzz[k] * ab_x + g_x_x_x_xyyzzz[k];

                g_x_x_xx_yzzzz[k] = -g_0_x_x_yzzzz[k] + g_x_0_x_yzzzz[k] - g_x_x_x_yzzzz[k] * ab_x + g_x_x_x_xyzzzz[k];

                g_x_x_xx_zzzzz[k] = -g_0_x_x_zzzzz[k] + g_x_0_x_zzzzz[k] - g_x_x_x_zzzzz[k] * ab_x + g_x_x_x_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_x_x_xy_xxxxx = cbuffer.data(dh_geom_11_off + 21 * ccomps * dcomps);

            auto g_x_x_xy_xxxxy = cbuffer.data(dh_geom_11_off + 22 * ccomps * dcomps);

            auto g_x_x_xy_xxxxz = cbuffer.data(dh_geom_11_off + 23 * ccomps * dcomps);

            auto g_x_x_xy_xxxyy = cbuffer.data(dh_geom_11_off + 24 * ccomps * dcomps);

            auto g_x_x_xy_xxxyz = cbuffer.data(dh_geom_11_off + 25 * ccomps * dcomps);

            auto g_x_x_xy_xxxzz = cbuffer.data(dh_geom_11_off + 26 * ccomps * dcomps);

            auto g_x_x_xy_xxyyy = cbuffer.data(dh_geom_11_off + 27 * ccomps * dcomps);

            auto g_x_x_xy_xxyyz = cbuffer.data(dh_geom_11_off + 28 * ccomps * dcomps);

            auto g_x_x_xy_xxyzz = cbuffer.data(dh_geom_11_off + 29 * ccomps * dcomps);

            auto g_x_x_xy_xxzzz = cbuffer.data(dh_geom_11_off + 30 * ccomps * dcomps);

            auto g_x_x_xy_xyyyy = cbuffer.data(dh_geom_11_off + 31 * ccomps * dcomps);

            auto g_x_x_xy_xyyyz = cbuffer.data(dh_geom_11_off + 32 * ccomps * dcomps);

            auto g_x_x_xy_xyyzz = cbuffer.data(dh_geom_11_off + 33 * ccomps * dcomps);

            auto g_x_x_xy_xyzzz = cbuffer.data(dh_geom_11_off + 34 * ccomps * dcomps);

            auto g_x_x_xy_xzzzz = cbuffer.data(dh_geom_11_off + 35 * ccomps * dcomps);

            auto g_x_x_xy_yyyyy = cbuffer.data(dh_geom_11_off + 36 * ccomps * dcomps);

            auto g_x_x_xy_yyyyz = cbuffer.data(dh_geom_11_off + 37 * ccomps * dcomps);

            auto g_x_x_xy_yyyzz = cbuffer.data(dh_geom_11_off + 38 * ccomps * dcomps);

            auto g_x_x_xy_yyzzz = cbuffer.data(dh_geom_11_off + 39 * ccomps * dcomps);

            auto g_x_x_xy_yzzzz = cbuffer.data(dh_geom_11_off + 40 * ccomps * dcomps);

            auto g_x_x_xy_zzzzz = cbuffer.data(dh_geom_11_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_xxxxx, g_x_x_x_xxxxxy, g_x_x_x_xxxxy, g_x_x_x_xxxxyy, g_x_x_x_xxxxyz, g_x_x_x_xxxxz, g_x_x_x_xxxyy, g_x_x_x_xxxyyy, g_x_x_x_xxxyyz, g_x_x_x_xxxyz, g_x_x_x_xxxyzz, g_x_x_x_xxxzz, g_x_x_x_xxyyy, g_x_x_x_xxyyyy, g_x_x_x_xxyyyz, g_x_x_x_xxyyz, g_x_x_x_xxyyzz, g_x_x_x_xxyzz, g_x_x_x_xxyzzz, g_x_x_x_xxzzz, g_x_x_x_xyyyy, g_x_x_x_xyyyyy, g_x_x_x_xyyyyz, g_x_x_x_xyyyz, g_x_x_x_xyyyzz, g_x_x_x_xyyzz, g_x_x_x_xyyzzz, g_x_x_x_xyzzz, g_x_x_x_xyzzzz, g_x_x_x_xzzzz, g_x_x_x_yyyyy, g_x_x_x_yyyyyy, g_x_x_x_yyyyyz, g_x_x_x_yyyyz, g_x_x_x_yyyyzz, g_x_x_x_yyyzz, g_x_x_x_yyyzzz, g_x_x_x_yyzzz, g_x_x_x_yyzzzz, g_x_x_x_yzzzz, g_x_x_x_yzzzzz, g_x_x_x_zzzzz, g_x_x_xy_xxxxx, g_x_x_xy_xxxxy, g_x_x_xy_xxxxz, g_x_x_xy_xxxyy, g_x_x_xy_xxxyz, g_x_x_xy_xxxzz, g_x_x_xy_xxyyy, g_x_x_xy_xxyyz, g_x_x_xy_xxyzz, g_x_x_xy_xxzzz, g_x_x_xy_xyyyy, g_x_x_xy_xyyyz, g_x_x_xy_xyyzz, g_x_x_xy_xyzzz, g_x_x_xy_xzzzz, g_x_x_xy_yyyyy, g_x_x_xy_yyyyz, g_x_x_xy_yyyzz, g_x_x_xy_yyzzz, g_x_x_xy_yzzzz, g_x_x_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xy_xxxxx[k] = -g_x_x_x_xxxxx[k] * ab_y + g_x_x_x_xxxxxy[k];

                g_x_x_xy_xxxxy[k] = -g_x_x_x_xxxxy[k] * ab_y + g_x_x_x_xxxxyy[k];

                g_x_x_xy_xxxxz[k] = -g_x_x_x_xxxxz[k] * ab_y + g_x_x_x_xxxxyz[k];

                g_x_x_xy_xxxyy[k] = -g_x_x_x_xxxyy[k] * ab_y + g_x_x_x_xxxyyy[k];

                g_x_x_xy_xxxyz[k] = -g_x_x_x_xxxyz[k] * ab_y + g_x_x_x_xxxyyz[k];

                g_x_x_xy_xxxzz[k] = -g_x_x_x_xxxzz[k] * ab_y + g_x_x_x_xxxyzz[k];

                g_x_x_xy_xxyyy[k] = -g_x_x_x_xxyyy[k] * ab_y + g_x_x_x_xxyyyy[k];

                g_x_x_xy_xxyyz[k] = -g_x_x_x_xxyyz[k] * ab_y + g_x_x_x_xxyyyz[k];

                g_x_x_xy_xxyzz[k] = -g_x_x_x_xxyzz[k] * ab_y + g_x_x_x_xxyyzz[k];

                g_x_x_xy_xxzzz[k] = -g_x_x_x_xxzzz[k] * ab_y + g_x_x_x_xxyzzz[k];

                g_x_x_xy_xyyyy[k] = -g_x_x_x_xyyyy[k] * ab_y + g_x_x_x_xyyyyy[k];

                g_x_x_xy_xyyyz[k] = -g_x_x_x_xyyyz[k] * ab_y + g_x_x_x_xyyyyz[k];

                g_x_x_xy_xyyzz[k] = -g_x_x_x_xyyzz[k] * ab_y + g_x_x_x_xyyyzz[k];

                g_x_x_xy_xyzzz[k] = -g_x_x_x_xyzzz[k] * ab_y + g_x_x_x_xyyzzz[k];

                g_x_x_xy_xzzzz[k] = -g_x_x_x_xzzzz[k] * ab_y + g_x_x_x_xyzzzz[k];

                g_x_x_xy_yyyyy[k] = -g_x_x_x_yyyyy[k] * ab_y + g_x_x_x_yyyyyy[k];

                g_x_x_xy_yyyyz[k] = -g_x_x_x_yyyyz[k] * ab_y + g_x_x_x_yyyyyz[k];

                g_x_x_xy_yyyzz[k] = -g_x_x_x_yyyzz[k] * ab_y + g_x_x_x_yyyyzz[k];

                g_x_x_xy_yyzzz[k] = -g_x_x_x_yyzzz[k] * ab_y + g_x_x_x_yyyzzz[k];

                g_x_x_xy_yzzzz[k] = -g_x_x_x_yzzzz[k] * ab_y + g_x_x_x_yyzzzz[k];

                g_x_x_xy_zzzzz[k] = -g_x_x_x_zzzzz[k] * ab_y + g_x_x_x_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_x_x_xz_xxxxx = cbuffer.data(dh_geom_11_off + 42 * ccomps * dcomps);

            auto g_x_x_xz_xxxxy = cbuffer.data(dh_geom_11_off + 43 * ccomps * dcomps);

            auto g_x_x_xz_xxxxz = cbuffer.data(dh_geom_11_off + 44 * ccomps * dcomps);

            auto g_x_x_xz_xxxyy = cbuffer.data(dh_geom_11_off + 45 * ccomps * dcomps);

            auto g_x_x_xz_xxxyz = cbuffer.data(dh_geom_11_off + 46 * ccomps * dcomps);

            auto g_x_x_xz_xxxzz = cbuffer.data(dh_geom_11_off + 47 * ccomps * dcomps);

            auto g_x_x_xz_xxyyy = cbuffer.data(dh_geom_11_off + 48 * ccomps * dcomps);

            auto g_x_x_xz_xxyyz = cbuffer.data(dh_geom_11_off + 49 * ccomps * dcomps);

            auto g_x_x_xz_xxyzz = cbuffer.data(dh_geom_11_off + 50 * ccomps * dcomps);

            auto g_x_x_xz_xxzzz = cbuffer.data(dh_geom_11_off + 51 * ccomps * dcomps);

            auto g_x_x_xz_xyyyy = cbuffer.data(dh_geom_11_off + 52 * ccomps * dcomps);

            auto g_x_x_xz_xyyyz = cbuffer.data(dh_geom_11_off + 53 * ccomps * dcomps);

            auto g_x_x_xz_xyyzz = cbuffer.data(dh_geom_11_off + 54 * ccomps * dcomps);

            auto g_x_x_xz_xyzzz = cbuffer.data(dh_geom_11_off + 55 * ccomps * dcomps);

            auto g_x_x_xz_xzzzz = cbuffer.data(dh_geom_11_off + 56 * ccomps * dcomps);

            auto g_x_x_xz_yyyyy = cbuffer.data(dh_geom_11_off + 57 * ccomps * dcomps);

            auto g_x_x_xz_yyyyz = cbuffer.data(dh_geom_11_off + 58 * ccomps * dcomps);

            auto g_x_x_xz_yyyzz = cbuffer.data(dh_geom_11_off + 59 * ccomps * dcomps);

            auto g_x_x_xz_yyzzz = cbuffer.data(dh_geom_11_off + 60 * ccomps * dcomps);

            auto g_x_x_xz_yzzzz = cbuffer.data(dh_geom_11_off + 61 * ccomps * dcomps);

            auto g_x_x_xz_zzzzz = cbuffer.data(dh_geom_11_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_x_xxxxx, g_x_x_x_xxxxxz, g_x_x_x_xxxxy, g_x_x_x_xxxxyz, g_x_x_x_xxxxz, g_x_x_x_xxxxzz, g_x_x_x_xxxyy, g_x_x_x_xxxyyz, g_x_x_x_xxxyz, g_x_x_x_xxxyzz, g_x_x_x_xxxzz, g_x_x_x_xxxzzz, g_x_x_x_xxyyy, g_x_x_x_xxyyyz, g_x_x_x_xxyyz, g_x_x_x_xxyyzz, g_x_x_x_xxyzz, g_x_x_x_xxyzzz, g_x_x_x_xxzzz, g_x_x_x_xxzzzz, g_x_x_x_xyyyy, g_x_x_x_xyyyyz, g_x_x_x_xyyyz, g_x_x_x_xyyyzz, g_x_x_x_xyyzz, g_x_x_x_xyyzzz, g_x_x_x_xyzzz, g_x_x_x_xyzzzz, g_x_x_x_xzzzz, g_x_x_x_xzzzzz, g_x_x_x_yyyyy, g_x_x_x_yyyyyz, g_x_x_x_yyyyz, g_x_x_x_yyyyzz, g_x_x_x_yyyzz, g_x_x_x_yyyzzz, g_x_x_x_yyzzz, g_x_x_x_yyzzzz, g_x_x_x_yzzzz, g_x_x_x_yzzzzz, g_x_x_x_zzzzz, g_x_x_x_zzzzzz, g_x_x_xz_xxxxx, g_x_x_xz_xxxxy, g_x_x_xz_xxxxz, g_x_x_xz_xxxyy, g_x_x_xz_xxxyz, g_x_x_xz_xxxzz, g_x_x_xz_xxyyy, g_x_x_xz_xxyyz, g_x_x_xz_xxyzz, g_x_x_xz_xxzzz, g_x_x_xz_xyyyy, g_x_x_xz_xyyyz, g_x_x_xz_xyyzz, g_x_x_xz_xyzzz, g_x_x_xz_xzzzz, g_x_x_xz_yyyyy, g_x_x_xz_yyyyz, g_x_x_xz_yyyzz, g_x_x_xz_yyzzz, g_x_x_xz_yzzzz, g_x_x_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_xz_xxxxx[k] = -g_x_x_x_xxxxx[k] * ab_z + g_x_x_x_xxxxxz[k];

                g_x_x_xz_xxxxy[k] = -g_x_x_x_xxxxy[k] * ab_z + g_x_x_x_xxxxyz[k];

                g_x_x_xz_xxxxz[k] = -g_x_x_x_xxxxz[k] * ab_z + g_x_x_x_xxxxzz[k];

                g_x_x_xz_xxxyy[k] = -g_x_x_x_xxxyy[k] * ab_z + g_x_x_x_xxxyyz[k];

                g_x_x_xz_xxxyz[k] = -g_x_x_x_xxxyz[k] * ab_z + g_x_x_x_xxxyzz[k];

                g_x_x_xz_xxxzz[k] = -g_x_x_x_xxxzz[k] * ab_z + g_x_x_x_xxxzzz[k];

                g_x_x_xz_xxyyy[k] = -g_x_x_x_xxyyy[k] * ab_z + g_x_x_x_xxyyyz[k];

                g_x_x_xz_xxyyz[k] = -g_x_x_x_xxyyz[k] * ab_z + g_x_x_x_xxyyzz[k];

                g_x_x_xz_xxyzz[k] = -g_x_x_x_xxyzz[k] * ab_z + g_x_x_x_xxyzzz[k];

                g_x_x_xz_xxzzz[k] = -g_x_x_x_xxzzz[k] * ab_z + g_x_x_x_xxzzzz[k];

                g_x_x_xz_xyyyy[k] = -g_x_x_x_xyyyy[k] * ab_z + g_x_x_x_xyyyyz[k];

                g_x_x_xz_xyyyz[k] = -g_x_x_x_xyyyz[k] * ab_z + g_x_x_x_xyyyzz[k];

                g_x_x_xz_xyyzz[k] = -g_x_x_x_xyyzz[k] * ab_z + g_x_x_x_xyyzzz[k];

                g_x_x_xz_xyzzz[k] = -g_x_x_x_xyzzz[k] * ab_z + g_x_x_x_xyzzzz[k];

                g_x_x_xz_xzzzz[k] = -g_x_x_x_xzzzz[k] * ab_z + g_x_x_x_xzzzzz[k];

                g_x_x_xz_yyyyy[k] = -g_x_x_x_yyyyy[k] * ab_z + g_x_x_x_yyyyyz[k];

                g_x_x_xz_yyyyz[k] = -g_x_x_x_yyyyz[k] * ab_z + g_x_x_x_yyyyzz[k];

                g_x_x_xz_yyyzz[k] = -g_x_x_x_yyyzz[k] * ab_z + g_x_x_x_yyyzzz[k];

                g_x_x_xz_yyzzz[k] = -g_x_x_x_yyzzz[k] * ab_z + g_x_x_x_yyzzzz[k];

                g_x_x_xz_yzzzz[k] = -g_x_x_x_yzzzz[k] * ab_z + g_x_x_x_yzzzzz[k];

                g_x_x_xz_zzzzz[k] = -g_x_x_x_zzzzz[k] * ab_z + g_x_x_x_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_x_x_yy_xxxxx = cbuffer.data(dh_geom_11_off + 63 * ccomps * dcomps);

            auto g_x_x_yy_xxxxy = cbuffer.data(dh_geom_11_off + 64 * ccomps * dcomps);

            auto g_x_x_yy_xxxxz = cbuffer.data(dh_geom_11_off + 65 * ccomps * dcomps);

            auto g_x_x_yy_xxxyy = cbuffer.data(dh_geom_11_off + 66 * ccomps * dcomps);

            auto g_x_x_yy_xxxyz = cbuffer.data(dh_geom_11_off + 67 * ccomps * dcomps);

            auto g_x_x_yy_xxxzz = cbuffer.data(dh_geom_11_off + 68 * ccomps * dcomps);

            auto g_x_x_yy_xxyyy = cbuffer.data(dh_geom_11_off + 69 * ccomps * dcomps);

            auto g_x_x_yy_xxyyz = cbuffer.data(dh_geom_11_off + 70 * ccomps * dcomps);

            auto g_x_x_yy_xxyzz = cbuffer.data(dh_geom_11_off + 71 * ccomps * dcomps);

            auto g_x_x_yy_xxzzz = cbuffer.data(dh_geom_11_off + 72 * ccomps * dcomps);

            auto g_x_x_yy_xyyyy = cbuffer.data(dh_geom_11_off + 73 * ccomps * dcomps);

            auto g_x_x_yy_xyyyz = cbuffer.data(dh_geom_11_off + 74 * ccomps * dcomps);

            auto g_x_x_yy_xyyzz = cbuffer.data(dh_geom_11_off + 75 * ccomps * dcomps);

            auto g_x_x_yy_xyzzz = cbuffer.data(dh_geom_11_off + 76 * ccomps * dcomps);

            auto g_x_x_yy_xzzzz = cbuffer.data(dh_geom_11_off + 77 * ccomps * dcomps);

            auto g_x_x_yy_yyyyy = cbuffer.data(dh_geom_11_off + 78 * ccomps * dcomps);

            auto g_x_x_yy_yyyyz = cbuffer.data(dh_geom_11_off + 79 * ccomps * dcomps);

            auto g_x_x_yy_yyyzz = cbuffer.data(dh_geom_11_off + 80 * ccomps * dcomps);

            auto g_x_x_yy_yyzzz = cbuffer.data(dh_geom_11_off + 81 * ccomps * dcomps);

            auto g_x_x_yy_yzzzz = cbuffer.data(dh_geom_11_off + 82 * ccomps * dcomps);

            auto g_x_x_yy_zzzzz = cbuffer.data(dh_geom_11_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_y_xxxxx, g_x_x_y_xxxxxy, g_x_x_y_xxxxy, g_x_x_y_xxxxyy, g_x_x_y_xxxxyz, g_x_x_y_xxxxz, g_x_x_y_xxxyy, g_x_x_y_xxxyyy, g_x_x_y_xxxyyz, g_x_x_y_xxxyz, g_x_x_y_xxxyzz, g_x_x_y_xxxzz, g_x_x_y_xxyyy, g_x_x_y_xxyyyy, g_x_x_y_xxyyyz, g_x_x_y_xxyyz, g_x_x_y_xxyyzz, g_x_x_y_xxyzz, g_x_x_y_xxyzzz, g_x_x_y_xxzzz, g_x_x_y_xyyyy, g_x_x_y_xyyyyy, g_x_x_y_xyyyyz, g_x_x_y_xyyyz, g_x_x_y_xyyyzz, g_x_x_y_xyyzz, g_x_x_y_xyyzzz, g_x_x_y_xyzzz, g_x_x_y_xyzzzz, g_x_x_y_xzzzz, g_x_x_y_yyyyy, g_x_x_y_yyyyyy, g_x_x_y_yyyyyz, g_x_x_y_yyyyz, g_x_x_y_yyyyzz, g_x_x_y_yyyzz, g_x_x_y_yyyzzz, g_x_x_y_yyzzz, g_x_x_y_yyzzzz, g_x_x_y_yzzzz, g_x_x_y_yzzzzz, g_x_x_y_zzzzz, g_x_x_yy_xxxxx, g_x_x_yy_xxxxy, g_x_x_yy_xxxxz, g_x_x_yy_xxxyy, g_x_x_yy_xxxyz, g_x_x_yy_xxxzz, g_x_x_yy_xxyyy, g_x_x_yy_xxyyz, g_x_x_yy_xxyzz, g_x_x_yy_xxzzz, g_x_x_yy_xyyyy, g_x_x_yy_xyyyz, g_x_x_yy_xyyzz, g_x_x_yy_xyzzz, g_x_x_yy_xzzzz, g_x_x_yy_yyyyy, g_x_x_yy_yyyyz, g_x_x_yy_yyyzz, g_x_x_yy_yyzzz, g_x_x_yy_yzzzz, g_x_x_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yy_xxxxx[k] = -g_x_x_y_xxxxx[k] * ab_y + g_x_x_y_xxxxxy[k];

                g_x_x_yy_xxxxy[k] = -g_x_x_y_xxxxy[k] * ab_y + g_x_x_y_xxxxyy[k];

                g_x_x_yy_xxxxz[k] = -g_x_x_y_xxxxz[k] * ab_y + g_x_x_y_xxxxyz[k];

                g_x_x_yy_xxxyy[k] = -g_x_x_y_xxxyy[k] * ab_y + g_x_x_y_xxxyyy[k];

                g_x_x_yy_xxxyz[k] = -g_x_x_y_xxxyz[k] * ab_y + g_x_x_y_xxxyyz[k];

                g_x_x_yy_xxxzz[k] = -g_x_x_y_xxxzz[k] * ab_y + g_x_x_y_xxxyzz[k];

                g_x_x_yy_xxyyy[k] = -g_x_x_y_xxyyy[k] * ab_y + g_x_x_y_xxyyyy[k];

                g_x_x_yy_xxyyz[k] = -g_x_x_y_xxyyz[k] * ab_y + g_x_x_y_xxyyyz[k];

                g_x_x_yy_xxyzz[k] = -g_x_x_y_xxyzz[k] * ab_y + g_x_x_y_xxyyzz[k];

                g_x_x_yy_xxzzz[k] = -g_x_x_y_xxzzz[k] * ab_y + g_x_x_y_xxyzzz[k];

                g_x_x_yy_xyyyy[k] = -g_x_x_y_xyyyy[k] * ab_y + g_x_x_y_xyyyyy[k];

                g_x_x_yy_xyyyz[k] = -g_x_x_y_xyyyz[k] * ab_y + g_x_x_y_xyyyyz[k];

                g_x_x_yy_xyyzz[k] = -g_x_x_y_xyyzz[k] * ab_y + g_x_x_y_xyyyzz[k];

                g_x_x_yy_xyzzz[k] = -g_x_x_y_xyzzz[k] * ab_y + g_x_x_y_xyyzzz[k];

                g_x_x_yy_xzzzz[k] = -g_x_x_y_xzzzz[k] * ab_y + g_x_x_y_xyzzzz[k];

                g_x_x_yy_yyyyy[k] = -g_x_x_y_yyyyy[k] * ab_y + g_x_x_y_yyyyyy[k];

                g_x_x_yy_yyyyz[k] = -g_x_x_y_yyyyz[k] * ab_y + g_x_x_y_yyyyyz[k];

                g_x_x_yy_yyyzz[k] = -g_x_x_y_yyyzz[k] * ab_y + g_x_x_y_yyyyzz[k];

                g_x_x_yy_yyzzz[k] = -g_x_x_y_yyzzz[k] * ab_y + g_x_x_y_yyyzzz[k];

                g_x_x_yy_yzzzz[k] = -g_x_x_y_yzzzz[k] * ab_y + g_x_x_y_yyzzzz[k];

                g_x_x_yy_zzzzz[k] = -g_x_x_y_zzzzz[k] * ab_y + g_x_x_y_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_x_x_yz_xxxxx = cbuffer.data(dh_geom_11_off + 84 * ccomps * dcomps);

            auto g_x_x_yz_xxxxy = cbuffer.data(dh_geom_11_off + 85 * ccomps * dcomps);

            auto g_x_x_yz_xxxxz = cbuffer.data(dh_geom_11_off + 86 * ccomps * dcomps);

            auto g_x_x_yz_xxxyy = cbuffer.data(dh_geom_11_off + 87 * ccomps * dcomps);

            auto g_x_x_yz_xxxyz = cbuffer.data(dh_geom_11_off + 88 * ccomps * dcomps);

            auto g_x_x_yz_xxxzz = cbuffer.data(dh_geom_11_off + 89 * ccomps * dcomps);

            auto g_x_x_yz_xxyyy = cbuffer.data(dh_geom_11_off + 90 * ccomps * dcomps);

            auto g_x_x_yz_xxyyz = cbuffer.data(dh_geom_11_off + 91 * ccomps * dcomps);

            auto g_x_x_yz_xxyzz = cbuffer.data(dh_geom_11_off + 92 * ccomps * dcomps);

            auto g_x_x_yz_xxzzz = cbuffer.data(dh_geom_11_off + 93 * ccomps * dcomps);

            auto g_x_x_yz_xyyyy = cbuffer.data(dh_geom_11_off + 94 * ccomps * dcomps);

            auto g_x_x_yz_xyyyz = cbuffer.data(dh_geom_11_off + 95 * ccomps * dcomps);

            auto g_x_x_yz_xyyzz = cbuffer.data(dh_geom_11_off + 96 * ccomps * dcomps);

            auto g_x_x_yz_xyzzz = cbuffer.data(dh_geom_11_off + 97 * ccomps * dcomps);

            auto g_x_x_yz_xzzzz = cbuffer.data(dh_geom_11_off + 98 * ccomps * dcomps);

            auto g_x_x_yz_yyyyy = cbuffer.data(dh_geom_11_off + 99 * ccomps * dcomps);

            auto g_x_x_yz_yyyyz = cbuffer.data(dh_geom_11_off + 100 * ccomps * dcomps);

            auto g_x_x_yz_yyyzz = cbuffer.data(dh_geom_11_off + 101 * ccomps * dcomps);

            auto g_x_x_yz_yyzzz = cbuffer.data(dh_geom_11_off + 102 * ccomps * dcomps);

            auto g_x_x_yz_yzzzz = cbuffer.data(dh_geom_11_off + 103 * ccomps * dcomps);

            auto g_x_x_yz_zzzzz = cbuffer.data(dh_geom_11_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_yz_xxxxx, g_x_x_yz_xxxxy, g_x_x_yz_xxxxz, g_x_x_yz_xxxyy, g_x_x_yz_xxxyz, g_x_x_yz_xxxzz, g_x_x_yz_xxyyy, g_x_x_yz_xxyyz, g_x_x_yz_xxyzz, g_x_x_yz_xxzzz, g_x_x_yz_xyyyy, g_x_x_yz_xyyyz, g_x_x_yz_xyyzz, g_x_x_yz_xyzzz, g_x_x_yz_xzzzz, g_x_x_yz_yyyyy, g_x_x_yz_yyyyz, g_x_x_yz_yyyzz, g_x_x_yz_yyzzz, g_x_x_yz_yzzzz, g_x_x_yz_zzzzz, g_x_x_z_xxxxx, g_x_x_z_xxxxxy, g_x_x_z_xxxxy, g_x_x_z_xxxxyy, g_x_x_z_xxxxyz, g_x_x_z_xxxxz, g_x_x_z_xxxyy, g_x_x_z_xxxyyy, g_x_x_z_xxxyyz, g_x_x_z_xxxyz, g_x_x_z_xxxyzz, g_x_x_z_xxxzz, g_x_x_z_xxyyy, g_x_x_z_xxyyyy, g_x_x_z_xxyyyz, g_x_x_z_xxyyz, g_x_x_z_xxyyzz, g_x_x_z_xxyzz, g_x_x_z_xxyzzz, g_x_x_z_xxzzz, g_x_x_z_xyyyy, g_x_x_z_xyyyyy, g_x_x_z_xyyyyz, g_x_x_z_xyyyz, g_x_x_z_xyyyzz, g_x_x_z_xyyzz, g_x_x_z_xyyzzz, g_x_x_z_xyzzz, g_x_x_z_xyzzzz, g_x_x_z_xzzzz, g_x_x_z_yyyyy, g_x_x_z_yyyyyy, g_x_x_z_yyyyyz, g_x_x_z_yyyyz, g_x_x_z_yyyyzz, g_x_x_z_yyyzz, g_x_x_z_yyyzzz, g_x_x_z_yyzzz, g_x_x_z_yyzzzz, g_x_x_z_yzzzz, g_x_x_z_yzzzzz, g_x_x_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_yz_xxxxx[k] = -g_x_x_z_xxxxx[k] * ab_y + g_x_x_z_xxxxxy[k];

                g_x_x_yz_xxxxy[k] = -g_x_x_z_xxxxy[k] * ab_y + g_x_x_z_xxxxyy[k];

                g_x_x_yz_xxxxz[k] = -g_x_x_z_xxxxz[k] * ab_y + g_x_x_z_xxxxyz[k];

                g_x_x_yz_xxxyy[k] = -g_x_x_z_xxxyy[k] * ab_y + g_x_x_z_xxxyyy[k];

                g_x_x_yz_xxxyz[k] = -g_x_x_z_xxxyz[k] * ab_y + g_x_x_z_xxxyyz[k];

                g_x_x_yz_xxxzz[k] = -g_x_x_z_xxxzz[k] * ab_y + g_x_x_z_xxxyzz[k];

                g_x_x_yz_xxyyy[k] = -g_x_x_z_xxyyy[k] * ab_y + g_x_x_z_xxyyyy[k];

                g_x_x_yz_xxyyz[k] = -g_x_x_z_xxyyz[k] * ab_y + g_x_x_z_xxyyyz[k];

                g_x_x_yz_xxyzz[k] = -g_x_x_z_xxyzz[k] * ab_y + g_x_x_z_xxyyzz[k];

                g_x_x_yz_xxzzz[k] = -g_x_x_z_xxzzz[k] * ab_y + g_x_x_z_xxyzzz[k];

                g_x_x_yz_xyyyy[k] = -g_x_x_z_xyyyy[k] * ab_y + g_x_x_z_xyyyyy[k];

                g_x_x_yz_xyyyz[k] = -g_x_x_z_xyyyz[k] * ab_y + g_x_x_z_xyyyyz[k];

                g_x_x_yz_xyyzz[k] = -g_x_x_z_xyyzz[k] * ab_y + g_x_x_z_xyyyzz[k];

                g_x_x_yz_xyzzz[k] = -g_x_x_z_xyzzz[k] * ab_y + g_x_x_z_xyyzzz[k];

                g_x_x_yz_xzzzz[k] = -g_x_x_z_xzzzz[k] * ab_y + g_x_x_z_xyzzzz[k];

                g_x_x_yz_yyyyy[k] = -g_x_x_z_yyyyy[k] * ab_y + g_x_x_z_yyyyyy[k];

                g_x_x_yz_yyyyz[k] = -g_x_x_z_yyyyz[k] * ab_y + g_x_x_z_yyyyyz[k];

                g_x_x_yz_yyyzz[k] = -g_x_x_z_yyyzz[k] * ab_y + g_x_x_z_yyyyzz[k];

                g_x_x_yz_yyzzz[k] = -g_x_x_z_yyzzz[k] * ab_y + g_x_x_z_yyyzzz[k];

                g_x_x_yz_yzzzz[k] = -g_x_x_z_yzzzz[k] * ab_y + g_x_x_z_yyzzzz[k];

                g_x_x_yz_zzzzz[k] = -g_x_x_z_zzzzz[k] * ab_y + g_x_x_z_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_x_x_zz_xxxxx = cbuffer.data(dh_geom_11_off + 105 * ccomps * dcomps);

            auto g_x_x_zz_xxxxy = cbuffer.data(dh_geom_11_off + 106 * ccomps * dcomps);

            auto g_x_x_zz_xxxxz = cbuffer.data(dh_geom_11_off + 107 * ccomps * dcomps);

            auto g_x_x_zz_xxxyy = cbuffer.data(dh_geom_11_off + 108 * ccomps * dcomps);

            auto g_x_x_zz_xxxyz = cbuffer.data(dh_geom_11_off + 109 * ccomps * dcomps);

            auto g_x_x_zz_xxxzz = cbuffer.data(dh_geom_11_off + 110 * ccomps * dcomps);

            auto g_x_x_zz_xxyyy = cbuffer.data(dh_geom_11_off + 111 * ccomps * dcomps);

            auto g_x_x_zz_xxyyz = cbuffer.data(dh_geom_11_off + 112 * ccomps * dcomps);

            auto g_x_x_zz_xxyzz = cbuffer.data(dh_geom_11_off + 113 * ccomps * dcomps);

            auto g_x_x_zz_xxzzz = cbuffer.data(dh_geom_11_off + 114 * ccomps * dcomps);

            auto g_x_x_zz_xyyyy = cbuffer.data(dh_geom_11_off + 115 * ccomps * dcomps);

            auto g_x_x_zz_xyyyz = cbuffer.data(dh_geom_11_off + 116 * ccomps * dcomps);

            auto g_x_x_zz_xyyzz = cbuffer.data(dh_geom_11_off + 117 * ccomps * dcomps);

            auto g_x_x_zz_xyzzz = cbuffer.data(dh_geom_11_off + 118 * ccomps * dcomps);

            auto g_x_x_zz_xzzzz = cbuffer.data(dh_geom_11_off + 119 * ccomps * dcomps);

            auto g_x_x_zz_yyyyy = cbuffer.data(dh_geom_11_off + 120 * ccomps * dcomps);

            auto g_x_x_zz_yyyyz = cbuffer.data(dh_geom_11_off + 121 * ccomps * dcomps);

            auto g_x_x_zz_yyyzz = cbuffer.data(dh_geom_11_off + 122 * ccomps * dcomps);

            auto g_x_x_zz_yyzzz = cbuffer.data(dh_geom_11_off + 123 * ccomps * dcomps);

            auto g_x_x_zz_yzzzz = cbuffer.data(dh_geom_11_off + 124 * ccomps * dcomps);

            auto g_x_x_zz_zzzzz = cbuffer.data(dh_geom_11_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_x_z_xxxxx, g_x_x_z_xxxxxz, g_x_x_z_xxxxy, g_x_x_z_xxxxyz, g_x_x_z_xxxxz, g_x_x_z_xxxxzz, g_x_x_z_xxxyy, g_x_x_z_xxxyyz, g_x_x_z_xxxyz, g_x_x_z_xxxyzz, g_x_x_z_xxxzz, g_x_x_z_xxxzzz, g_x_x_z_xxyyy, g_x_x_z_xxyyyz, g_x_x_z_xxyyz, g_x_x_z_xxyyzz, g_x_x_z_xxyzz, g_x_x_z_xxyzzz, g_x_x_z_xxzzz, g_x_x_z_xxzzzz, g_x_x_z_xyyyy, g_x_x_z_xyyyyz, g_x_x_z_xyyyz, g_x_x_z_xyyyzz, g_x_x_z_xyyzz, g_x_x_z_xyyzzz, g_x_x_z_xyzzz, g_x_x_z_xyzzzz, g_x_x_z_xzzzz, g_x_x_z_xzzzzz, g_x_x_z_yyyyy, g_x_x_z_yyyyyz, g_x_x_z_yyyyz, g_x_x_z_yyyyzz, g_x_x_z_yyyzz, g_x_x_z_yyyzzz, g_x_x_z_yyzzz, g_x_x_z_yyzzzz, g_x_x_z_yzzzz, g_x_x_z_yzzzzz, g_x_x_z_zzzzz, g_x_x_z_zzzzzz, g_x_x_zz_xxxxx, g_x_x_zz_xxxxy, g_x_x_zz_xxxxz, g_x_x_zz_xxxyy, g_x_x_zz_xxxyz, g_x_x_zz_xxxzz, g_x_x_zz_xxyyy, g_x_x_zz_xxyyz, g_x_x_zz_xxyzz, g_x_x_zz_xxzzz, g_x_x_zz_xyyyy, g_x_x_zz_xyyyz, g_x_x_zz_xyyzz, g_x_x_zz_xyzzz, g_x_x_zz_xzzzz, g_x_x_zz_yyyyy, g_x_x_zz_yyyyz, g_x_x_zz_yyyzz, g_x_x_zz_yyzzz, g_x_x_zz_yzzzz, g_x_x_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_x_zz_xxxxx[k] = -g_x_x_z_xxxxx[k] * ab_z + g_x_x_z_xxxxxz[k];

                g_x_x_zz_xxxxy[k] = -g_x_x_z_xxxxy[k] * ab_z + g_x_x_z_xxxxyz[k];

                g_x_x_zz_xxxxz[k] = -g_x_x_z_xxxxz[k] * ab_z + g_x_x_z_xxxxzz[k];

                g_x_x_zz_xxxyy[k] = -g_x_x_z_xxxyy[k] * ab_z + g_x_x_z_xxxyyz[k];

                g_x_x_zz_xxxyz[k] = -g_x_x_z_xxxyz[k] * ab_z + g_x_x_z_xxxyzz[k];

                g_x_x_zz_xxxzz[k] = -g_x_x_z_xxxzz[k] * ab_z + g_x_x_z_xxxzzz[k];

                g_x_x_zz_xxyyy[k] = -g_x_x_z_xxyyy[k] * ab_z + g_x_x_z_xxyyyz[k];

                g_x_x_zz_xxyyz[k] = -g_x_x_z_xxyyz[k] * ab_z + g_x_x_z_xxyyzz[k];

                g_x_x_zz_xxyzz[k] = -g_x_x_z_xxyzz[k] * ab_z + g_x_x_z_xxyzzz[k];

                g_x_x_zz_xxzzz[k] = -g_x_x_z_xxzzz[k] * ab_z + g_x_x_z_xxzzzz[k];

                g_x_x_zz_xyyyy[k] = -g_x_x_z_xyyyy[k] * ab_z + g_x_x_z_xyyyyz[k];

                g_x_x_zz_xyyyz[k] = -g_x_x_z_xyyyz[k] * ab_z + g_x_x_z_xyyyzz[k];

                g_x_x_zz_xyyzz[k] = -g_x_x_z_xyyzz[k] * ab_z + g_x_x_z_xyyzzz[k];

                g_x_x_zz_xyzzz[k] = -g_x_x_z_xyzzz[k] * ab_z + g_x_x_z_xyzzzz[k];

                g_x_x_zz_xzzzz[k] = -g_x_x_z_xzzzz[k] * ab_z + g_x_x_z_xzzzzz[k];

                g_x_x_zz_yyyyy[k] = -g_x_x_z_yyyyy[k] * ab_z + g_x_x_z_yyyyyz[k];

                g_x_x_zz_yyyyz[k] = -g_x_x_z_yyyyz[k] * ab_z + g_x_x_z_yyyyzz[k];

                g_x_x_zz_yyyzz[k] = -g_x_x_z_yyyzz[k] * ab_z + g_x_x_z_yyyzzz[k];

                g_x_x_zz_yyzzz[k] = -g_x_x_z_yyzzz[k] * ab_z + g_x_x_z_yyzzzz[k];

                g_x_x_zz_yzzzz[k] = -g_x_x_z_yzzzz[k] * ab_z + g_x_x_z_yzzzzz[k];

                g_x_x_zz_zzzzz[k] = -g_x_x_z_zzzzz[k] * ab_z + g_x_x_z_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_x_y_xx_xxxxx = cbuffer.data(dh_geom_11_off + 126 * ccomps * dcomps);

            auto g_x_y_xx_xxxxy = cbuffer.data(dh_geom_11_off + 127 * ccomps * dcomps);

            auto g_x_y_xx_xxxxz = cbuffer.data(dh_geom_11_off + 128 * ccomps * dcomps);

            auto g_x_y_xx_xxxyy = cbuffer.data(dh_geom_11_off + 129 * ccomps * dcomps);

            auto g_x_y_xx_xxxyz = cbuffer.data(dh_geom_11_off + 130 * ccomps * dcomps);

            auto g_x_y_xx_xxxzz = cbuffer.data(dh_geom_11_off + 131 * ccomps * dcomps);

            auto g_x_y_xx_xxyyy = cbuffer.data(dh_geom_11_off + 132 * ccomps * dcomps);

            auto g_x_y_xx_xxyyz = cbuffer.data(dh_geom_11_off + 133 * ccomps * dcomps);

            auto g_x_y_xx_xxyzz = cbuffer.data(dh_geom_11_off + 134 * ccomps * dcomps);

            auto g_x_y_xx_xxzzz = cbuffer.data(dh_geom_11_off + 135 * ccomps * dcomps);

            auto g_x_y_xx_xyyyy = cbuffer.data(dh_geom_11_off + 136 * ccomps * dcomps);

            auto g_x_y_xx_xyyyz = cbuffer.data(dh_geom_11_off + 137 * ccomps * dcomps);

            auto g_x_y_xx_xyyzz = cbuffer.data(dh_geom_11_off + 138 * ccomps * dcomps);

            auto g_x_y_xx_xyzzz = cbuffer.data(dh_geom_11_off + 139 * ccomps * dcomps);

            auto g_x_y_xx_xzzzz = cbuffer.data(dh_geom_11_off + 140 * ccomps * dcomps);

            auto g_x_y_xx_yyyyy = cbuffer.data(dh_geom_11_off + 141 * ccomps * dcomps);

            auto g_x_y_xx_yyyyz = cbuffer.data(dh_geom_11_off + 142 * ccomps * dcomps);

            auto g_x_y_xx_yyyzz = cbuffer.data(dh_geom_11_off + 143 * ccomps * dcomps);

            auto g_x_y_xx_yyzzz = cbuffer.data(dh_geom_11_off + 144 * ccomps * dcomps);

            auto g_x_y_xx_yzzzz = cbuffer.data(dh_geom_11_off + 145 * ccomps * dcomps);

            auto g_x_y_xx_zzzzz = cbuffer.data(dh_geom_11_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_x_xxxxx, g_0_y_x_xxxxy, g_0_y_x_xxxxz, g_0_y_x_xxxyy, g_0_y_x_xxxyz, g_0_y_x_xxxzz, g_0_y_x_xxyyy, g_0_y_x_xxyyz, g_0_y_x_xxyzz, g_0_y_x_xxzzz, g_0_y_x_xyyyy, g_0_y_x_xyyyz, g_0_y_x_xyyzz, g_0_y_x_xyzzz, g_0_y_x_xzzzz, g_0_y_x_yyyyy, g_0_y_x_yyyyz, g_0_y_x_yyyzz, g_0_y_x_yyzzz, g_0_y_x_yzzzz, g_0_y_x_zzzzz, g_x_y_x_xxxxx, g_x_y_x_xxxxxx, g_x_y_x_xxxxxy, g_x_y_x_xxxxxz, g_x_y_x_xxxxy, g_x_y_x_xxxxyy, g_x_y_x_xxxxyz, g_x_y_x_xxxxz, g_x_y_x_xxxxzz, g_x_y_x_xxxyy, g_x_y_x_xxxyyy, g_x_y_x_xxxyyz, g_x_y_x_xxxyz, g_x_y_x_xxxyzz, g_x_y_x_xxxzz, g_x_y_x_xxxzzz, g_x_y_x_xxyyy, g_x_y_x_xxyyyy, g_x_y_x_xxyyyz, g_x_y_x_xxyyz, g_x_y_x_xxyyzz, g_x_y_x_xxyzz, g_x_y_x_xxyzzz, g_x_y_x_xxzzz, g_x_y_x_xxzzzz, g_x_y_x_xyyyy, g_x_y_x_xyyyyy, g_x_y_x_xyyyyz, g_x_y_x_xyyyz, g_x_y_x_xyyyzz, g_x_y_x_xyyzz, g_x_y_x_xyyzzz, g_x_y_x_xyzzz, g_x_y_x_xyzzzz, g_x_y_x_xzzzz, g_x_y_x_xzzzzz, g_x_y_x_yyyyy, g_x_y_x_yyyyz, g_x_y_x_yyyzz, g_x_y_x_yyzzz, g_x_y_x_yzzzz, g_x_y_x_zzzzz, g_x_y_xx_xxxxx, g_x_y_xx_xxxxy, g_x_y_xx_xxxxz, g_x_y_xx_xxxyy, g_x_y_xx_xxxyz, g_x_y_xx_xxxzz, g_x_y_xx_xxyyy, g_x_y_xx_xxyyz, g_x_y_xx_xxyzz, g_x_y_xx_xxzzz, g_x_y_xx_xyyyy, g_x_y_xx_xyyyz, g_x_y_xx_xyyzz, g_x_y_xx_xyzzz, g_x_y_xx_xzzzz, g_x_y_xx_yyyyy, g_x_y_xx_yyyyz, g_x_y_xx_yyyzz, g_x_y_xx_yyzzz, g_x_y_xx_yzzzz, g_x_y_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xx_xxxxx[k] = -g_0_y_x_xxxxx[k] - g_x_y_x_xxxxx[k] * ab_x + g_x_y_x_xxxxxx[k];

                g_x_y_xx_xxxxy[k] = -g_0_y_x_xxxxy[k] - g_x_y_x_xxxxy[k] * ab_x + g_x_y_x_xxxxxy[k];

                g_x_y_xx_xxxxz[k] = -g_0_y_x_xxxxz[k] - g_x_y_x_xxxxz[k] * ab_x + g_x_y_x_xxxxxz[k];

                g_x_y_xx_xxxyy[k] = -g_0_y_x_xxxyy[k] - g_x_y_x_xxxyy[k] * ab_x + g_x_y_x_xxxxyy[k];

                g_x_y_xx_xxxyz[k] = -g_0_y_x_xxxyz[k] - g_x_y_x_xxxyz[k] * ab_x + g_x_y_x_xxxxyz[k];

                g_x_y_xx_xxxzz[k] = -g_0_y_x_xxxzz[k] - g_x_y_x_xxxzz[k] * ab_x + g_x_y_x_xxxxzz[k];

                g_x_y_xx_xxyyy[k] = -g_0_y_x_xxyyy[k] - g_x_y_x_xxyyy[k] * ab_x + g_x_y_x_xxxyyy[k];

                g_x_y_xx_xxyyz[k] = -g_0_y_x_xxyyz[k] - g_x_y_x_xxyyz[k] * ab_x + g_x_y_x_xxxyyz[k];

                g_x_y_xx_xxyzz[k] = -g_0_y_x_xxyzz[k] - g_x_y_x_xxyzz[k] * ab_x + g_x_y_x_xxxyzz[k];

                g_x_y_xx_xxzzz[k] = -g_0_y_x_xxzzz[k] - g_x_y_x_xxzzz[k] * ab_x + g_x_y_x_xxxzzz[k];

                g_x_y_xx_xyyyy[k] = -g_0_y_x_xyyyy[k] - g_x_y_x_xyyyy[k] * ab_x + g_x_y_x_xxyyyy[k];

                g_x_y_xx_xyyyz[k] = -g_0_y_x_xyyyz[k] - g_x_y_x_xyyyz[k] * ab_x + g_x_y_x_xxyyyz[k];

                g_x_y_xx_xyyzz[k] = -g_0_y_x_xyyzz[k] - g_x_y_x_xyyzz[k] * ab_x + g_x_y_x_xxyyzz[k];

                g_x_y_xx_xyzzz[k] = -g_0_y_x_xyzzz[k] - g_x_y_x_xyzzz[k] * ab_x + g_x_y_x_xxyzzz[k];

                g_x_y_xx_xzzzz[k] = -g_0_y_x_xzzzz[k] - g_x_y_x_xzzzz[k] * ab_x + g_x_y_x_xxzzzz[k];

                g_x_y_xx_yyyyy[k] = -g_0_y_x_yyyyy[k] - g_x_y_x_yyyyy[k] * ab_x + g_x_y_x_xyyyyy[k];

                g_x_y_xx_yyyyz[k] = -g_0_y_x_yyyyz[k] - g_x_y_x_yyyyz[k] * ab_x + g_x_y_x_xyyyyz[k];

                g_x_y_xx_yyyzz[k] = -g_0_y_x_yyyzz[k] - g_x_y_x_yyyzz[k] * ab_x + g_x_y_x_xyyyzz[k];

                g_x_y_xx_yyzzz[k] = -g_0_y_x_yyzzz[k] - g_x_y_x_yyzzz[k] * ab_x + g_x_y_x_xyyzzz[k];

                g_x_y_xx_yzzzz[k] = -g_0_y_x_yzzzz[k] - g_x_y_x_yzzzz[k] * ab_x + g_x_y_x_xyzzzz[k];

                g_x_y_xx_zzzzz[k] = -g_0_y_x_zzzzz[k] - g_x_y_x_zzzzz[k] * ab_x + g_x_y_x_xzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_x_y_xy_xxxxx = cbuffer.data(dh_geom_11_off + 147 * ccomps * dcomps);

            auto g_x_y_xy_xxxxy = cbuffer.data(dh_geom_11_off + 148 * ccomps * dcomps);

            auto g_x_y_xy_xxxxz = cbuffer.data(dh_geom_11_off + 149 * ccomps * dcomps);

            auto g_x_y_xy_xxxyy = cbuffer.data(dh_geom_11_off + 150 * ccomps * dcomps);

            auto g_x_y_xy_xxxyz = cbuffer.data(dh_geom_11_off + 151 * ccomps * dcomps);

            auto g_x_y_xy_xxxzz = cbuffer.data(dh_geom_11_off + 152 * ccomps * dcomps);

            auto g_x_y_xy_xxyyy = cbuffer.data(dh_geom_11_off + 153 * ccomps * dcomps);

            auto g_x_y_xy_xxyyz = cbuffer.data(dh_geom_11_off + 154 * ccomps * dcomps);

            auto g_x_y_xy_xxyzz = cbuffer.data(dh_geom_11_off + 155 * ccomps * dcomps);

            auto g_x_y_xy_xxzzz = cbuffer.data(dh_geom_11_off + 156 * ccomps * dcomps);

            auto g_x_y_xy_xyyyy = cbuffer.data(dh_geom_11_off + 157 * ccomps * dcomps);

            auto g_x_y_xy_xyyyz = cbuffer.data(dh_geom_11_off + 158 * ccomps * dcomps);

            auto g_x_y_xy_xyyzz = cbuffer.data(dh_geom_11_off + 159 * ccomps * dcomps);

            auto g_x_y_xy_xyzzz = cbuffer.data(dh_geom_11_off + 160 * ccomps * dcomps);

            auto g_x_y_xy_xzzzz = cbuffer.data(dh_geom_11_off + 161 * ccomps * dcomps);

            auto g_x_y_xy_yyyyy = cbuffer.data(dh_geom_11_off + 162 * ccomps * dcomps);

            auto g_x_y_xy_yyyyz = cbuffer.data(dh_geom_11_off + 163 * ccomps * dcomps);

            auto g_x_y_xy_yyyzz = cbuffer.data(dh_geom_11_off + 164 * ccomps * dcomps);

            auto g_x_y_xy_yyzzz = cbuffer.data(dh_geom_11_off + 165 * ccomps * dcomps);

            auto g_x_y_xy_yzzzz = cbuffer.data(dh_geom_11_off + 166 * ccomps * dcomps);

            auto g_x_y_xy_zzzzz = cbuffer.data(dh_geom_11_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxxxx, g_0_y_y_xxxxy, g_0_y_y_xxxxz, g_0_y_y_xxxyy, g_0_y_y_xxxyz, g_0_y_y_xxxzz, g_0_y_y_xxyyy, g_0_y_y_xxyyz, g_0_y_y_xxyzz, g_0_y_y_xxzzz, g_0_y_y_xyyyy, g_0_y_y_xyyyz, g_0_y_y_xyyzz, g_0_y_y_xyzzz, g_0_y_y_xzzzz, g_0_y_y_yyyyy, g_0_y_y_yyyyz, g_0_y_y_yyyzz, g_0_y_y_yyzzz, g_0_y_y_yzzzz, g_0_y_y_zzzzz, g_x_y_xy_xxxxx, g_x_y_xy_xxxxy, g_x_y_xy_xxxxz, g_x_y_xy_xxxyy, g_x_y_xy_xxxyz, g_x_y_xy_xxxzz, g_x_y_xy_xxyyy, g_x_y_xy_xxyyz, g_x_y_xy_xxyzz, g_x_y_xy_xxzzz, g_x_y_xy_xyyyy, g_x_y_xy_xyyyz, g_x_y_xy_xyyzz, g_x_y_xy_xyzzz, g_x_y_xy_xzzzz, g_x_y_xy_yyyyy, g_x_y_xy_yyyyz, g_x_y_xy_yyyzz, g_x_y_xy_yyzzz, g_x_y_xy_yzzzz, g_x_y_xy_zzzzz, g_x_y_y_xxxxx, g_x_y_y_xxxxxx, g_x_y_y_xxxxxy, g_x_y_y_xxxxxz, g_x_y_y_xxxxy, g_x_y_y_xxxxyy, g_x_y_y_xxxxyz, g_x_y_y_xxxxz, g_x_y_y_xxxxzz, g_x_y_y_xxxyy, g_x_y_y_xxxyyy, g_x_y_y_xxxyyz, g_x_y_y_xxxyz, g_x_y_y_xxxyzz, g_x_y_y_xxxzz, g_x_y_y_xxxzzz, g_x_y_y_xxyyy, g_x_y_y_xxyyyy, g_x_y_y_xxyyyz, g_x_y_y_xxyyz, g_x_y_y_xxyyzz, g_x_y_y_xxyzz, g_x_y_y_xxyzzz, g_x_y_y_xxzzz, g_x_y_y_xxzzzz, g_x_y_y_xyyyy, g_x_y_y_xyyyyy, g_x_y_y_xyyyyz, g_x_y_y_xyyyz, g_x_y_y_xyyyzz, g_x_y_y_xyyzz, g_x_y_y_xyyzzz, g_x_y_y_xyzzz, g_x_y_y_xyzzzz, g_x_y_y_xzzzz, g_x_y_y_xzzzzz, g_x_y_y_yyyyy, g_x_y_y_yyyyz, g_x_y_y_yyyzz, g_x_y_y_yyzzz, g_x_y_y_yzzzz, g_x_y_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xy_xxxxx[k] = -g_0_y_y_xxxxx[k] - g_x_y_y_xxxxx[k] * ab_x + g_x_y_y_xxxxxx[k];

                g_x_y_xy_xxxxy[k] = -g_0_y_y_xxxxy[k] - g_x_y_y_xxxxy[k] * ab_x + g_x_y_y_xxxxxy[k];

                g_x_y_xy_xxxxz[k] = -g_0_y_y_xxxxz[k] - g_x_y_y_xxxxz[k] * ab_x + g_x_y_y_xxxxxz[k];

                g_x_y_xy_xxxyy[k] = -g_0_y_y_xxxyy[k] - g_x_y_y_xxxyy[k] * ab_x + g_x_y_y_xxxxyy[k];

                g_x_y_xy_xxxyz[k] = -g_0_y_y_xxxyz[k] - g_x_y_y_xxxyz[k] * ab_x + g_x_y_y_xxxxyz[k];

                g_x_y_xy_xxxzz[k] = -g_0_y_y_xxxzz[k] - g_x_y_y_xxxzz[k] * ab_x + g_x_y_y_xxxxzz[k];

                g_x_y_xy_xxyyy[k] = -g_0_y_y_xxyyy[k] - g_x_y_y_xxyyy[k] * ab_x + g_x_y_y_xxxyyy[k];

                g_x_y_xy_xxyyz[k] = -g_0_y_y_xxyyz[k] - g_x_y_y_xxyyz[k] * ab_x + g_x_y_y_xxxyyz[k];

                g_x_y_xy_xxyzz[k] = -g_0_y_y_xxyzz[k] - g_x_y_y_xxyzz[k] * ab_x + g_x_y_y_xxxyzz[k];

                g_x_y_xy_xxzzz[k] = -g_0_y_y_xxzzz[k] - g_x_y_y_xxzzz[k] * ab_x + g_x_y_y_xxxzzz[k];

                g_x_y_xy_xyyyy[k] = -g_0_y_y_xyyyy[k] - g_x_y_y_xyyyy[k] * ab_x + g_x_y_y_xxyyyy[k];

                g_x_y_xy_xyyyz[k] = -g_0_y_y_xyyyz[k] - g_x_y_y_xyyyz[k] * ab_x + g_x_y_y_xxyyyz[k];

                g_x_y_xy_xyyzz[k] = -g_0_y_y_xyyzz[k] - g_x_y_y_xyyzz[k] * ab_x + g_x_y_y_xxyyzz[k];

                g_x_y_xy_xyzzz[k] = -g_0_y_y_xyzzz[k] - g_x_y_y_xyzzz[k] * ab_x + g_x_y_y_xxyzzz[k];

                g_x_y_xy_xzzzz[k] = -g_0_y_y_xzzzz[k] - g_x_y_y_xzzzz[k] * ab_x + g_x_y_y_xxzzzz[k];

                g_x_y_xy_yyyyy[k] = -g_0_y_y_yyyyy[k] - g_x_y_y_yyyyy[k] * ab_x + g_x_y_y_xyyyyy[k];

                g_x_y_xy_yyyyz[k] = -g_0_y_y_yyyyz[k] - g_x_y_y_yyyyz[k] * ab_x + g_x_y_y_xyyyyz[k];

                g_x_y_xy_yyyzz[k] = -g_0_y_y_yyyzz[k] - g_x_y_y_yyyzz[k] * ab_x + g_x_y_y_xyyyzz[k];

                g_x_y_xy_yyzzz[k] = -g_0_y_y_yyzzz[k] - g_x_y_y_yyzzz[k] * ab_x + g_x_y_y_xyyzzz[k];

                g_x_y_xy_yzzzz[k] = -g_0_y_y_yzzzz[k] - g_x_y_y_yzzzz[k] * ab_x + g_x_y_y_xyzzzz[k];

                g_x_y_xy_zzzzz[k] = -g_0_y_y_zzzzz[k] - g_x_y_y_zzzzz[k] * ab_x + g_x_y_y_xzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_x_y_xz_xxxxx = cbuffer.data(dh_geom_11_off + 168 * ccomps * dcomps);

            auto g_x_y_xz_xxxxy = cbuffer.data(dh_geom_11_off + 169 * ccomps * dcomps);

            auto g_x_y_xz_xxxxz = cbuffer.data(dh_geom_11_off + 170 * ccomps * dcomps);

            auto g_x_y_xz_xxxyy = cbuffer.data(dh_geom_11_off + 171 * ccomps * dcomps);

            auto g_x_y_xz_xxxyz = cbuffer.data(dh_geom_11_off + 172 * ccomps * dcomps);

            auto g_x_y_xz_xxxzz = cbuffer.data(dh_geom_11_off + 173 * ccomps * dcomps);

            auto g_x_y_xz_xxyyy = cbuffer.data(dh_geom_11_off + 174 * ccomps * dcomps);

            auto g_x_y_xz_xxyyz = cbuffer.data(dh_geom_11_off + 175 * ccomps * dcomps);

            auto g_x_y_xz_xxyzz = cbuffer.data(dh_geom_11_off + 176 * ccomps * dcomps);

            auto g_x_y_xz_xxzzz = cbuffer.data(dh_geom_11_off + 177 * ccomps * dcomps);

            auto g_x_y_xz_xyyyy = cbuffer.data(dh_geom_11_off + 178 * ccomps * dcomps);

            auto g_x_y_xz_xyyyz = cbuffer.data(dh_geom_11_off + 179 * ccomps * dcomps);

            auto g_x_y_xz_xyyzz = cbuffer.data(dh_geom_11_off + 180 * ccomps * dcomps);

            auto g_x_y_xz_xyzzz = cbuffer.data(dh_geom_11_off + 181 * ccomps * dcomps);

            auto g_x_y_xz_xzzzz = cbuffer.data(dh_geom_11_off + 182 * ccomps * dcomps);

            auto g_x_y_xz_yyyyy = cbuffer.data(dh_geom_11_off + 183 * ccomps * dcomps);

            auto g_x_y_xz_yyyyz = cbuffer.data(dh_geom_11_off + 184 * ccomps * dcomps);

            auto g_x_y_xz_yyyzz = cbuffer.data(dh_geom_11_off + 185 * ccomps * dcomps);

            auto g_x_y_xz_yyzzz = cbuffer.data(dh_geom_11_off + 186 * ccomps * dcomps);

            auto g_x_y_xz_yzzzz = cbuffer.data(dh_geom_11_off + 187 * ccomps * dcomps);

            auto g_x_y_xz_zzzzz = cbuffer.data(dh_geom_11_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_x_xxxxx, g_x_y_x_xxxxxz, g_x_y_x_xxxxy, g_x_y_x_xxxxyz, g_x_y_x_xxxxz, g_x_y_x_xxxxzz, g_x_y_x_xxxyy, g_x_y_x_xxxyyz, g_x_y_x_xxxyz, g_x_y_x_xxxyzz, g_x_y_x_xxxzz, g_x_y_x_xxxzzz, g_x_y_x_xxyyy, g_x_y_x_xxyyyz, g_x_y_x_xxyyz, g_x_y_x_xxyyzz, g_x_y_x_xxyzz, g_x_y_x_xxyzzz, g_x_y_x_xxzzz, g_x_y_x_xxzzzz, g_x_y_x_xyyyy, g_x_y_x_xyyyyz, g_x_y_x_xyyyz, g_x_y_x_xyyyzz, g_x_y_x_xyyzz, g_x_y_x_xyyzzz, g_x_y_x_xyzzz, g_x_y_x_xyzzzz, g_x_y_x_xzzzz, g_x_y_x_xzzzzz, g_x_y_x_yyyyy, g_x_y_x_yyyyyz, g_x_y_x_yyyyz, g_x_y_x_yyyyzz, g_x_y_x_yyyzz, g_x_y_x_yyyzzz, g_x_y_x_yyzzz, g_x_y_x_yyzzzz, g_x_y_x_yzzzz, g_x_y_x_yzzzzz, g_x_y_x_zzzzz, g_x_y_x_zzzzzz, g_x_y_xz_xxxxx, g_x_y_xz_xxxxy, g_x_y_xz_xxxxz, g_x_y_xz_xxxyy, g_x_y_xz_xxxyz, g_x_y_xz_xxxzz, g_x_y_xz_xxyyy, g_x_y_xz_xxyyz, g_x_y_xz_xxyzz, g_x_y_xz_xxzzz, g_x_y_xz_xyyyy, g_x_y_xz_xyyyz, g_x_y_xz_xyyzz, g_x_y_xz_xyzzz, g_x_y_xz_xzzzz, g_x_y_xz_yyyyy, g_x_y_xz_yyyyz, g_x_y_xz_yyyzz, g_x_y_xz_yyzzz, g_x_y_xz_yzzzz, g_x_y_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_xz_xxxxx[k] = -g_x_y_x_xxxxx[k] * ab_z + g_x_y_x_xxxxxz[k];

                g_x_y_xz_xxxxy[k] = -g_x_y_x_xxxxy[k] * ab_z + g_x_y_x_xxxxyz[k];

                g_x_y_xz_xxxxz[k] = -g_x_y_x_xxxxz[k] * ab_z + g_x_y_x_xxxxzz[k];

                g_x_y_xz_xxxyy[k] = -g_x_y_x_xxxyy[k] * ab_z + g_x_y_x_xxxyyz[k];

                g_x_y_xz_xxxyz[k] = -g_x_y_x_xxxyz[k] * ab_z + g_x_y_x_xxxyzz[k];

                g_x_y_xz_xxxzz[k] = -g_x_y_x_xxxzz[k] * ab_z + g_x_y_x_xxxzzz[k];

                g_x_y_xz_xxyyy[k] = -g_x_y_x_xxyyy[k] * ab_z + g_x_y_x_xxyyyz[k];

                g_x_y_xz_xxyyz[k] = -g_x_y_x_xxyyz[k] * ab_z + g_x_y_x_xxyyzz[k];

                g_x_y_xz_xxyzz[k] = -g_x_y_x_xxyzz[k] * ab_z + g_x_y_x_xxyzzz[k];

                g_x_y_xz_xxzzz[k] = -g_x_y_x_xxzzz[k] * ab_z + g_x_y_x_xxzzzz[k];

                g_x_y_xz_xyyyy[k] = -g_x_y_x_xyyyy[k] * ab_z + g_x_y_x_xyyyyz[k];

                g_x_y_xz_xyyyz[k] = -g_x_y_x_xyyyz[k] * ab_z + g_x_y_x_xyyyzz[k];

                g_x_y_xz_xyyzz[k] = -g_x_y_x_xyyzz[k] * ab_z + g_x_y_x_xyyzzz[k];

                g_x_y_xz_xyzzz[k] = -g_x_y_x_xyzzz[k] * ab_z + g_x_y_x_xyzzzz[k];

                g_x_y_xz_xzzzz[k] = -g_x_y_x_xzzzz[k] * ab_z + g_x_y_x_xzzzzz[k];

                g_x_y_xz_yyyyy[k] = -g_x_y_x_yyyyy[k] * ab_z + g_x_y_x_yyyyyz[k];

                g_x_y_xz_yyyyz[k] = -g_x_y_x_yyyyz[k] * ab_z + g_x_y_x_yyyyzz[k];

                g_x_y_xz_yyyzz[k] = -g_x_y_x_yyyzz[k] * ab_z + g_x_y_x_yyyzzz[k];

                g_x_y_xz_yyzzz[k] = -g_x_y_x_yyzzz[k] * ab_z + g_x_y_x_yyzzzz[k];

                g_x_y_xz_yzzzz[k] = -g_x_y_x_yzzzz[k] * ab_z + g_x_y_x_yzzzzz[k];

                g_x_y_xz_zzzzz[k] = -g_x_y_x_zzzzz[k] * ab_z + g_x_y_x_zzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_x_y_yy_xxxxx = cbuffer.data(dh_geom_11_off + 189 * ccomps * dcomps);

            auto g_x_y_yy_xxxxy = cbuffer.data(dh_geom_11_off + 190 * ccomps * dcomps);

            auto g_x_y_yy_xxxxz = cbuffer.data(dh_geom_11_off + 191 * ccomps * dcomps);

            auto g_x_y_yy_xxxyy = cbuffer.data(dh_geom_11_off + 192 * ccomps * dcomps);

            auto g_x_y_yy_xxxyz = cbuffer.data(dh_geom_11_off + 193 * ccomps * dcomps);

            auto g_x_y_yy_xxxzz = cbuffer.data(dh_geom_11_off + 194 * ccomps * dcomps);

            auto g_x_y_yy_xxyyy = cbuffer.data(dh_geom_11_off + 195 * ccomps * dcomps);

            auto g_x_y_yy_xxyyz = cbuffer.data(dh_geom_11_off + 196 * ccomps * dcomps);

            auto g_x_y_yy_xxyzz = cbuffer.data(dh_geom_11_off + 197 * ccomps * dcomps);

            auto g_x_y_yy_xxzzz = cbuffer.data(dh_geom_11_off + 198 * ccomps * dcomps);

            auto g_x_y_yy_xyyyy = cbuffer.data(dh_geom_11_off + 199 * ccomps * dcomps);

            auto g_x_y_yy_xyyyz = cbuffer.data(dh_geom_11_off + 200 * ccomps * dcomps);

            auto g_x_y_yy_xyyzz = cbuffer.data(dh_geom_11_off + 201 * ccomps * dcomps);

            auto g_x_y_yy_xyzzz = cbuffer.data(dh_geom_11_off + 202 * ccomps * dcomps);

            auto g_x_y_yy_xzzzz = cbuffer.data(dh_geom_11_off + 203 * ccomps * dcomps);

            auto g_x_y_yy_yyyyy = cbuffer.data(dh_geom_11_off + 204 * ccomps * dcomps);

            auto g_x_y_yy_yyyyz = cbuffer.data(dh_geom_11_off + 205 * ccomps * dcomps);

            auto g_x_y_yy_yyyzz = cbuffer.data(dh_geom_11_off + 206 * ccomps * dcomps);

            auto g_x_y_yy_yyzzz = cbuffer.data(dh_geom_11_off + 207 * ccomps * dcomps);

            auto g_x_y_yy_yzzzz = cbuffer.data(dh_geom_11_off + 208 * ccomps * dcomps);

            auto g_x_y_yy_zzzzz = cbuffer.data(dh_geom_11_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_y_xxxxx, g_x_0_y_xxxxy, g_x_0_y_xxxxz, g_x_0_y_xxxyy, g_x_0_y_xxxyz, g_x_0_y_xxxzz, g_x_0_y_xxyyy, g_x_0_y_xxyyz, g_x_0_y_xxyzz, g_x_0_y_xxzzz, g_x_0_y_xyyyy, g_x_0_y_xyyyz, g_x_0_y_xyyzz, g_x_0_y_xyzzz, g_x_0_y_xzzzz, g_x_0_y_yyyyy, g_x_0_y_yyyyz, g_x_0_y_yyyzz, g_x_0_y_yyzzz, g_x_0_y_yzzzz, g_x_0_y_zzzzz, g_x_y_y_xxxxx, g_x_y_y_xxxxxy, g_x_y_y_xxxxy, g_x_y_y_xxxxyy, g_x_y_y_xxxxyz, g_x_y_y_xxxxz, g_x_y_y_xxxyy, g_x_y_y_xxxyyy, g_x_y_y_xxxyyz, g_x_y_y_xxxyz, g_x_y_y_xxxyzz, g_x_y_y_xxxzz, g_x_y_y_xxyyy, g_x_y_y_xxyyyy, g_x_y_y_xxyyyz, g_x_y_y_xxyyz, g_x_y_y_xxyyzz, g_x_y_y_xxyzz, g_x_y_y_xxyzzz, g_x_y_y_xxzzz, g_x_y_y_xyyyy, g_x_y_y_xyyyyy, g_x_y_y_xyyyyz, g_x_y_y_xyyyz, g_x_y_y_xyyyzz, g_x_y_y_xyyzz, g_x_y_y_xyyzzz, g_x_y_y_xyzzz, g_x_y_y_xyzzzz, g_x_y_y_xzzzz, g_x_y_y_yyyyy, g_x_y_y_yyyyyy, g_x_y_y_yyyyyz, g_x_y_y_yyyyz, g_x_y_y_yyyyzz, g_x_y_y_yyyzz, g_x_y_y_yyyzzz, g_x_y_y_yyzzz, g_x_y_y_yyzzzz, g_x_y_y_yzzzz, g_x_y_y_yzzzzz, g_x_y_y_zzzzz, g_x_y_yy_xxxxx, g_x_y_yy_xxxxy, g_x_y_yy_xxxxz, g_x_y_yy_xxxyy, g_x_y_yy_xxxyz, g_x_y_yy_xxxzz, g_x_y_yy_xxyyy, g_x_y_yy_xxyyz, g_x_y_yy_xxyzz, g_x_y_yy_xxzzz, g_x_y_yy_xyyyy, g_x_y_yy_xyyyz, g_x_y_yy_xyyzz, g_x_y_yy_xyzzz, g_x_y_yy_xzzzz, g_x_y_yy_yyyyy, g_x_y_yy_yyyyz, g_x_y_yy_yyyzz, g_x_y_yy_yyzzz, g_x_y_yy_yzzzz, g_x_y_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yy_xxxxx[k] = g_x_0_y_xxxxx[k] - g_x_y_y_xxxxx[k] * ab_y + g_x_y_y_xxxxxy[k];

                g_x_y_yy_xxxxy[k] = g_x_0_y_xxxxy[k] - g_x_y_y_xxxxy[k] * ab_y + g_x_y_y_xxxxyy[k];

                g_x_y_yy_xxxxz[k] = g_x_0_y_xxxxz[k] - g_x_y_y_xxxxz[k] * ab_y + g_x_y_y_xxxxyz[k];

                g_x_y_yy_xxxyy[k] = g_x_0_y_xxxyy[k] - g_x_y_y_xxxyy[k] * ab_y + g_x_y_y_xxxyyy[k];

                g_x_y_yy_xxxyz[k] = g_x_0_y_xxxyz[k] - g_x_y_y_xxxyz[k] * ab_y + g_x_y_y_xxxyyz[k];

                g_x_y_yy_xxxzz[k] = g_x_0_y_xxxzz[k] - g_x_y_y_xxxzz[k] * ab_y + g_x_y_y_xxxyzz[k];

                g_x_y_yy_xxyyy[k] = g_x_0_y_xxyyy[k] - g_x_y_y_xxyyy[k] * ab_y + g_x_y_y_xxyyyy[k];

                g_x_y_yy_xxyyz[k] = g_x_0_y_xxyyz[k] - g_x_y_y_xxyyz[k] * ab_y + g_x_y_y_xxyyyz[k];

                g_x_y_yy_xxyzz[k] = g_x_0_y_xxyzz[k] - g_x_y_y_xxyzz[k] * ab_y + g_x_y_y_xxyyzz[k];

                g_x_y_yy_xxzzz[k] = g_x_0_y_xxzzz[k] - g_x_y_y_xxzzz[k] * ab_y + g_x_y_y_xxyzzz[k];

                g_x_y_yy_xyyyy[k] = g_x_0_y_xyyyy[k] - g_x_y_y_xyyyy[k] * ab_y + g_x_y_y_xyyyyy[k];

                g_x_y_yy_xyyyz[k] = g_x_0_y_xyyyz[k] - g_x_y_y_xyyyz[k] * ab_y + g_x_y_y_xyyyyz[k];

                g_x_y_yy_xyyzz[k] = g_x_0_y_xyyzz[k] - g_x_y_y_xyyzz[k] * ab_y + g_x_y_y_xyyyzz[k];

                g_x_y_yy_xyzzz[k] = g_x_0_y_xyzzz[k] - g_x_y_y_xyzzz[k] * ab_y + g_x_y_y_xyyzzz[k];

                g_x_y_yy_xzzzz[k] = g_x_0_y_xzzzz[k] - g_x_y_y_xzzzz[k] * ab_y + g_x_y_y_xyzzzz[k];

                g_x_y_yy_yyyyy[k] = g_x_0_y_yyyyy[k] - g_x_y_y_yyyyy[k] * ab_y + g_x_y_y_yyyyyy[k];

                g_x_y_yy_yyyyz[k] = g_x_0_y_yyyyz[k] - g_x_y_y_yyyyz[k] * ab_y + g_x_y_y_yyyyyz[k];

                g_x_y_yy_yyyzz[k] = g_x_0_y_yyyzz[k] - g_x_y_y_yyyzz[k] * ab_y + g_x_y_y_yyyyzz[k];

                g_x_y_yy_yyzzz[k] = g_x_0_y_yyzzz[k] - g_x_y_y_yyzzz[k] * ab_y + g_x_y_y_yyyzzz[k];

                g_x_y_yy_yzzzz[k] = g_x_0_y_yzzzz[k] - g_x_y_y_yzzzz[k] * ab_y + g_x_y_y_yyzzzz[k];

                g_x_y_yy_zzzzz[k] = g_x_0_y_zzzzz[k] - g_x_y_y_zzzzz[k] * ab_y + g_x_y_y_yzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_x_y_yz_xxxxx = cbuffer.data(dh_geom_11_off + 210 * ccomps * dcomps);

            auto g_x_y_yz_xxxxy = cbuffer.data(dh_geom_11_off + 211 * ccomps * dcomps);

            auto g_x_y_yz_xxxxz = cbuffer.data(dh_geom_11_off + 212 * ccomps * dcomps);

            auto g_x_y_yz_xxxyy = cbuffer.data(dh_geom_11_off + 213 * ccomps * dcomps);

            auto g_x_y_yz_xxxyz = cbuffer.data(dh_geom_11_off + 214 * ccomps * dcomps);

            auto g_x_y_yz_xxxzz = cbuffer.data(dh_geom_11_off + 215 * ccomps * dcomps);

            auto g_x_y_yz_xxyyy = cbuffer.data(dh_geom_11_off + 216 * ccomps * dcomps);

            auto g_x_y_yz_xxyyz = cbuffer.data(dh_geom_11_off + 217 * ccomps * dcomps);

            auto g_x_y_yz_xxyzz = cbuffer.data(dh_geom_11_off + 218 * ccomps * dcomps);

            auto g_x_y_yz_xxzzz = cbuffer.data(dh_geom_11_off + 219 * ccomps * dcomps);

            auto g_x_y_yz_xyyyy = cbuffer.data(dh_geom_11_off + 220 * ccomps * dcomps);

            auto g_x_y_yz_xyyyz = cbuffer.data(dh_geom_11_off + 221 * ccomps * dcomps);

            auto g_x_y_yz_xyyzz = cbuffer.data(dh_geom_11_off + 222 * ccomps * dcomps);

            auto g_x_y_yz_xyzzz = cbuffer.data(dh_geom_11_off + 223 * ccomps * dcomps);

            auto g_x_y_yz_xzzzz = cbuffer.data(dh_geom_11_off + 224 * ccomps * dcomps);

            auto g_x_y_yz_yyyyy = cbuffer.data(dh_geom_11_off + 225 * ccomps * dcomps);

            auto g_x_y_yz_yyyyz = cbuffer.data(dh_geom_11_off + 226 * ccomps * dcomps);

            auto g_x_y_yz_yyyzz = cbuffer.data(dh_geom_11_off + 227 * ccomps * dcomps);

            auto g_x_y_yz_yyzzz = cbuffer.data(dh_geom_11_off + 228 * ccomps * dcomps);

            auto g_x_y_yz_yzzzz = cbuffer.data(dh_geom_11_off + 229 * ccomps * dcomps);

            auto g_x_y_yz_zzzzz = cbuffer.data(dh_geom_11_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_y_xxxxx, g_x_y_y_xxxxxz, g_x_y_y_xxxxy, g_x_y_y_xxxxyz, g_x_y_y_xxxxz, g_x_y_y_xxxxzz, g_x_y_y_xxxyy, g_x_y_y_xxxyyz, g_x_y_y_xxxyz, g_x_y_y_xxxyzz, g_x_y_y_xxxzz, g_x_y_y_xxxzzz, g_x_y_y_xxyyy, g_x_y_y_xxyyyz, g_x_y_y_xxyyz, g_x_y_y_xxyyzz, g_x_y_y_xxyzz, g_x_y_y_xxyzzz, g_x_y_y_xxzzz, g_x_y_y_xxzzzz, g_x_y_y_xyyyy, g_x_y_y_xyyyyz, g_x_y_y_xyyyz, g_x_y_y_xyyyzz, g_x_y_y_xyyzz, g_x_y_y_xyyzzz, g_x_y_y_xyzzz, g_x_y_y_xyzzzz, g_x_y_y_xzzzz, g_x_y_y_xzzzzz, g_x_y_y_yyyyy, g_x_y_y_yyyyyz, g_x_y_y_yyyyz, g_x_y_y_yyyyzz, g_x_y_y_yyyzz, g_x_y_y_yyyzzz, g_x_y_y_yyzzz, g_x_y_y_yyzzzz, g_x_y_y_yzzzz, g_x_y_y_yzzzzz, g_x_y_y_zzzzz, g_x_y_y_zzzzzz, g_x_y_yz_xxxxx, g_x_y_yz_xxxxy, g_x_y_yz_xxxxz, g_x_y_yz_xxxyy, g_x_y_yz_xxxyz, g_x_y_yz_xxxzz, g_x_y_yz_xxyyy, g_x_y_yz_xxyyz, g_x_y_yz_xxyzz, g_x_y_yz_xxzzz, g_x_y_yz_xyyyy, g_x_y_yz_xyyyz, g_x_y_yz_xyyzz, g_x_y_yz_xyzzz, g_x_y_yz_xzzzz, g_x_y_yz_yyyyy, g_x_y_yz_yyyyz, g_x_y_yz_yyyzz, g_x_y_yz_yyzzz, g_x_y_yz_yzzzz, g_x_y_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_yz_xxxxx[k] = -g_x_y_y_xxxxx[k] * ab_z + g_x_y_y_xxxxxz[k];

                g_x_y_yz_xxxxy[k] = -g_x_y_y_xxxxy[k] * ab_z + g_x_y_y_xxxxyz[k];

                g_x_y_yz_xxxxz[k] = -g_x_y_y_xxxxz[k] * ab_z + g_x_y_y_xxxxzz[k];

                g_x_y_yz_xxxyy[k] = -g_x_y_y_xxxyy[k] * ab_z + g_x_y_y_xxxyyz[k];

                g_x_y_yz_xxxyz[k] = -g_x_y_y_xxxyz[k] * ab_z + g_x_y_y_xxxyzz[k];

                g_x_y_yz_xxxzz[k] = -g_x_y_y_xxxzz[k] * ab_z + g_x_y_y_xxxzzz[k];

                g_x_y_yz_xxyyy[k] = -g_x_y_y_xxyyy[k] * ab_z + g_x_y_y_xxyyyz[k];

                g_x_y_yz_xxyyz[k] = -g_x_y_y_xxyyz[k] * ab_z + g_x_y_y_xxyyzz[k];

                g_x_y_yz_xxyzz[k] = -g_x_y_y_xxyzz[k] * ab_z + g_x_y_y_xxyzzz[k];

                g_x_y_yz_xxzzz[k] = -g_x_y_y_xxzzz[k] * ab_z + g_x_y_y_xxzzzz[k];

                g_x_y_yz_xyyyy[k] = -g_x_y_y_xyyyy[k] * ab_z + g_x_y_y_xyyyyz[k];

                g_x_y_yz_xyyyz[k] = -g_x_y_y_xyyyz[k] * ab_z + g_x_y_y_xyyyzz[k];

                g_x_y_yz_xyyzz[k] = -g_x_y_y_xyyzz[k] * ab_z + g_x_y_y_xyyzzz[k];

                g_x_y_yz_xyzzz[k] = -g_x_y_y_xyzzz[k] * ab_z + g_x_y_y_xyzzzz[k];

                g_x_y_yz_xzzzz[k] = -g_x_y_y_xzzzz[k] * ab_z + g_x_y_y_xzzzzz[k];

                g_x_y_yz_yyyyy[k] = -g_x_y_y_yyyyy[k] * ab_z + g_x_y_y_yyyyyz[k];

                g_x_y_yz_yyyyz[k] = -g_x_y_y_yyyyz[k] * ab_z + g_x_y_y_yyyyzz[k];

                g_x_y_yz_yyyzz[k] = -g_x_y_y_yyyzz[k] * ab_z + g_x_y_y_yyyzzz[k];

                g_x_y_yz_yyzzz[k] = -g_x_y_y_yyzzz[k] * ab_z + g_x_y_y_yyzzzz[k];

                g_x_y_yz_yzzzz[k] = -g_x_y_y_yzzzz[k] * ab_z + g_x_y_y_yzzzzz[k];

                g_x_y_yz_zzzzz[k] = -g_x_y_y_zzzzz[k] * ab_z + g_x_y_y_zzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_x_y_zz_xxxxx = cbuffer.data(dh_geom_11_off + 231 * ccomps * dcomps);

            auto g_x_y_zz_xxxxy = cbuffer.data(dh_geom_11_off + 232 * ccomps * dcomps);

            auto g_x_y_zz_xxxxz = cbuffer.data(dh_geom_11_off + 233 * ccomps * dcomps);

            auto g_x_y_zz_xxxyy = cbuffer.data(dh_geom_11_off + 234 * ccomps * dcomps);

            auto g_x_y_zz_xxxyz = cbuffer.data(dh_geom_11_off + 235 * ccomps * dcomps);

            auto g_x_y_zz_xxxzz = cbuffer.data(dh_geom_11_off + 236 * ccomps * dcomps);

            auto g_x_y_zz_xxyyy = cbuffer.data(dh_geom_11_off + 237 * ccomps * dcomps);

            auto g_x_y_zz_xxyyz = cbuffer.data(dh_geom_11_off + 238 * ccomps * dcomps);

            auto g_x_y_zz_xxyzz = cbuffer.data(dh_geom_11_off + 239 * ccomps * dcomps);

            auto g_x_y_zz_xxzzz = cbuffer.data(dh_geom_11_off + 240 * ccomps * dcomps);

            auto g_x_y_zz_xyyyy = cbuffer.data(dh_geom_11_off + 241 * ccomps * dcomps);

            auto g_x_y_zz_xyyyz = cbuffer.data(dh_geom_11_off + 242 * ccomps * dcomps);

            auto g_x_y_zz_xyyzz = cbuffer.data(dh_geom_11_off + 243 * ccomps * dcomps);

            auto g_x_y_zz_xyzzz = cbuffer.data(dh_geom_11_off + 244 * ccomps * dcomps);

            auto g_x_y_zz_xzzzz = cbuffer.data(dh_geom_11_off + 245 * ccomps * dcomps);

            auto g_x_y_zz_yyyyy = cbuffer.data(dh_geom_11_off + 246 * ccomps * dcomps);

            auto g_x_y_zz_yyyyz = cbuffer.data(dh_geom_11_off + 247 * ccomps * dcomps);

            auto g_x_y_zz_yyyzz = cbuffer.data(dh_geom_11_off + 248 * ccomps * dcomps);

            auto g_x_y_zz_yyzzz = cbuffer.data(dh_geom_11_off + 249 * ccomps * dcomps);

            auto g_x_y_zz_yzzzz = cbuffer.data(dh_geom_11_off + 250 * ccomps * dcomps);

            auto g_x_y_zz_zzzzz = cbuffer.data(dh_geom_11_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_y_z_xxxxx, g_x_y_z_xxxxxz, g_x_y_z_xxxxy, g_x_y_z_xxxxyz, g_x_y_z_xxxxz, g_x_y_z_xxxxzz, g_x_y_z_xxxyy, g_x_y_z_xxxyyz, g_x_y_z_xxxyz, g_x_y_z_xxxyzz, g_x_y_z_xxxzz, g_x_y_z_xxxzzz, g_x_y_z_xxyyy, g_x_y_z_xxyyyz, g_x_y_z_xxyyz, g_x_y_z_xxyyzz, g_x_y_z_xxyzz, g_x_y_z_xxyzzz, g_x_y_z_xxzzz, g_x_y_z_xxzzzz, g_x_y_z_xyyyy, g_x_y_z_xyyyyz, g_x_y_z_xyyyz, g_x_y_z_xyyyzz, g_x_y_z_xyyzz, g_x_y_z_xyyzzz, g_x_y_z_xyzzz, g_x_y_z_xyzzzz, g_x_y_z_xzzzz, g_x_y_z_xzzzzz, g_x_y_z_yyyyy, g_x_y_z_yyyyyz, g_x_y_z_yyyyz, g_x_y_z_yyyyzz, g_x_y_z_yyyzz, g_x_y_z_yyyzzz, g_x_y_z_yyzzz, g_x_y_z_yyzzzz, g_x_y_z_yzzzz, g_x_y_z_yzzzzz, g_x_y_z_zzzzz, g_x_y_z_zzzzzz, g_x_y_zz_xxxxx, g_x_y_zz_xxxxy, g_x_y_zz_xxxxz, g_x_y_zz_xxxyy, g_x_y_zz_xxxyz, g_x_y_zz_xxxzz, g_x_y_zz_xxyyy, g_x_y_zz_xxyyz, g_x_y_zz_xxyzz, g_x_y_zz_xxzzz, g_x_y_zz_xyyyy, g_x_y_zz_xyyyz, g_x_y_zz_xyyzz, g_x_y_zz_xyzzz, g_x_y_zz_xzzzz, g_x_y_zz_yyyyy, g_x_y_zz_yyyyz, g_x_y_zz_yyyzz, g_x_y_zz_yyzzz, g_x_y_zz_yzzzz, g_x_y_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_y_zz_xxxxx[k] = -g_x_y_z_xxxxx[k] * ab_z + g_x_y_z_xxxxxz[k];

                g_x_y_zz_xxxxy[k] = -g_x_y_z_xxxxy[k] * ab_z + g_x_y_z_xxxxyz[k];

                g_x_y_zz_xxxxz[k] = -g_x_y_z_xxxxz[k] * ab_z + g_x_y_z_xxxxzz[k];

                g_x_y_zz_xxxyy[k] = -g_x_y_z_xxxyy[k] * ab_z + g_x_y_z_xxxyyz[k];

                g_x_y_zz_xxxyz[k] = -g_x_y_z_xxxyz[k] * ab_z + g_x_y_z_xxxyzz[k];

                g_x_y_zz_xxxzz[k] = -g_x_y_z_xxxzz[k] * ab_z + g_x_y_z_xxxzzz[k];

                g_x_y_zz_xxyyy[k] = -g_x_y_z_xxyyy[k] * ab_z + g_x_y_z_xxyyyz[k];

                g_x_y_zz_xxyyz[k] = -g_x_y_z_xxyyz[k] * ab_z + g_x_y_z_xxyyzz[k];

                g_x_y_zz_xxyzz[k] = -g_x_y_z_xxyzz[k] * ab_z + g_x_y_z_xxyzzz[k];

                g_x_y_zz_xxzzz[k] = -g_x_y_z_xxzzz[k] * ab_z + g_x_y_z_xxzzzz[k];

                g_x_y_zz_xyyyy[k] = -g_x_y_z_xyyyy[k] * ab_z + g_x_y_z_xyyyyz[k];

                g_x_y_zz_xyyyz[k] = -g_x_y_z_xyyyz[k] * ab_z + g_x_y_z_xyyyzz[k];

                g_x_y_zz_xyyzz[k] = -g_x_y_z_xyyzz[k] * ab_z + g_x_y_z_xyyzzz[k];

                g_x_y_zz_xyzzz[k] = -g_x_y_z_xyzzz[k] * ab_z + g_x_y_z_xyzzzz[k];

                g_x_y_zz_xzzzz[k] = -g_x_y_z_xzzzz[k] * ab_z + g_x_y_z_xzzzzz[k];

                g_x_y_zz_yyyyy[k] = -g_x_y_z_yyyyy[k] * ab_z + g_x_y_z_yyyyyz[k];

                g_x_y_zz_yyyyz[k] = -g_x_y_z_yyyyz[k] * ab_z + g_x_y_z_yyyyzz[k];

                g_x_y_zz_yyyzz[k] = -g_x_y_z_yyyzz[k] * ab_z + g_x_y_z_yyyzzz[k];

                g_x_y_zz_yyzzz[k] = -g_x_y_z_yyzzz[k] * ab_z + g_x_y_z_yyzzzz[k];

                g_x_y_zz_yzzzz[k] = -g_x_y_z_yzzzz[k] * ab_z + g_x_y_z_yzzzzz[k];

                g_x_y_zz_zzzzz[k] = -g_x_y_z_zzzzz[k] * ab_z + g_x_y_z_zzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_x_z_xx_xxxxx = cbuffer.data(dh_geom_11_off + 252 * ccomps * dcomps);

            auto g_x_z_xx_xxxxy = cbuffer.data(dh_geom_11_off + 253 * ccomps * dcomps);

            auto g_x_z_xx_xxxxz = cbuffer.data(dh_geom_11_off + 254 * ccomps * dcomps);

            auto g_x_z_xx_xxxyy = cbuffer.data(dh_geom_11_off + 255 * ccomps * dcomps);

            auto g_x_z_xx_xxxyz = cbuffer.data(dh_geom_11_off + 256 * ccomps * dcomps);

            auto g_x_z_xx_xxxzz = cbuffer.data(dh_geom_11_off + 257 * ccomps * dcomps);

            auto g_x_z_xx_xxyyy = cbuffer.data(dh_geom_11_off + 258 * ccomps * dcomps);

            auto g_x_z_xx_xxyyz = cbuffer.data(dh_geom_11_off + 259 * ccomps * dcomps);

            auto g_x_z_xx_xxyzz = cbuffer.data(dh_geom_11_off + 260 * ccomps * dcomps);

            auto g_x_z_xx_xxzzz = cbuffer.data(dh_geom_11_off + 261 * ccomps * dcomps);

            auto g_x_z_xx_xyyyy = cbuffer.data(dh_geom_11_off + 262 * ccomps * dcomps);

            auto g_x_z_xx_xyyyz = cbuffer.data(dh_geom_11_off + 263 * ccomps * dcomps);

            auto g_x_z_xx_xyyzz = cbuffer.data(dh_geom_11_off + 264 * ccomps * dcomps);

            auto g_x_z_xx_xyzzz = cbuffer.data(dh_geom_11_off + 265 * ccomps * dcomps);

            auto g_x_z_xx_xzzzz = cbuffer.data(dh_geom_11_off + 266 * ccomps * dcomps);

            auto g_x_z_xx_yyyyy = cbuffer.data(dh_geom_11_off + 267 * ccomps * dcomps);

            auto g_x_z_xx_yyyyz = cbuffer.data(dh_geom_11_off + 268 * ccomps * dcomps);

            auto g_x_z_xx_yyyzz = cbuffer.data(dh_geom_11_off + 269 * ccomps * dcomps);

            auto g_x_z_xx_yyzzz = cbuffer.data(dh_geom_11_off + 270 * ccomps * dcomps);

            auto g_x_z_xx_yzzzz = cbuffer.data(dh_geom_11_off + 271 * ccomps * dcomps);

            auto g_x_z_xx_zzzzz = cbuffer.data(dh_geom_11_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_x_xxxxx, g_0_z_x_xxxxy, g_0_z_x_xxxxz, g_0_z_x_xxxyy, g_0_z_x_xxxyz, g_0_z_x_xxxzz, g_0_z_x_xxyyy, g_0_z_x_xxyyz, g_0_z_x_xxyzz, g_0_z_x_xxzzz, g_0_z_x_xyyyy, g_0_z_x_xyyyz, g_0_z_x_xyyzz, g_0_z_x_xyzzz, g_0_z_x_xzzzz, g_0_z_x_yyyyy, g_0_z_x_yyyyz, g_0_z_x_yyyzz, g_0_z_x_yyzzz, g_0_z_x_yzzzz, g_0_z_x_zzzzz, g_x_z_x_xxxxx, g_x_z_x_xxxxxx, g_x_z_x_xxxxxy, g_x_z_x_xxxxxz, g_x_z_x_xxxxy, g_x_z_x_xxxxyy, g_x_z_x_xxxxyz, g_x_z_x_xxxxz, g_x_z_x_xxxxzz, g_x_z_x_xxxyy, g_x_z_x_xxxyyy, g_x_z_x_xxxyyz, g_x_z_x_xxxyz, g_x_z_x_xxxyzz, g_x_z_x_xxxzz, g_x_z_x_xxxzzz, g_x_z_x_xxyyy, g_x_z_x_xxyyyy, g_x_z_x_xxyyyz, g_x_z_x_xxyyz, g_x_z_x_xxyyzz, g_x_z_x_xxyzz, g_x_z_x_xxyzzz, g_x_z_x_xxzzz, g_x_z_x_xxzzzz, g_x_z_x_xyyyy, g_x_z_x_xyyyyy, g_x_z_x_xyyyyz, g_x_z_x_xyyyz, g_x_z_x_xyyyzz, g_x_z_x_xyyzz, g_x_z_x_xyyzzz, g_x_z_x_xyzzz, g_x_z_x_xyzzzz, g_x_z_x_xzzzz, g_x_z_x_xzzzzz, g_x_z_x_yyyyy, g_x_z_x_yyyyz, g_x_z_x_yyyzz, g_x_z_x_yyzzz, g_x_z_x_yzzzz, g_x_z_x_zzzzz, g_x_z_xx_xxxxx, g_x_z_xx_xxxxy, g_x_z_xx_xxxxz, g_x_z_xx_xxxyy, g_x_z_xx_xxxyz, g_x_z_xx_xxxzz, g_x_z_xx_xxyyy, g_x_z_xx_xxyyz, g_x_z_xx_xxyzz, g_x_z_xx_xxzzz, g_x_z_xx_xyyyy, g_x_z_xx_xyyyz, g_x_z_xx_xyyzz, g_x_z_xx_xyzzz, g_x_z_xx_xzzzz, g_x_z_xx_yyyyy, g_x_z_xx_yyyyz, g_x_z_xx_yyyzz, g_x_z_xx_yyzzz, g_x_z_xx_yzzzz, g_x_z_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xx_xxxxx[k] = -g_0_z_x_xxxxx[k] - g_x_z_x_xxxxx[k] * ab_x + g_x_z_x_xxxxxx[k];

                g_x_z_xx_xxxxy[k] = -g_0_z_x_xxxxy[k] - g_x_z_x_xxxxy[k] * ab_x + g_x_z_x_xxxxxy[k];

                g_x_z_xx_xxxxz[k] = -g_0_z_x_xxxxz[k] - g_x_z_x_xxxxz[k] * ab_x + g_x_z_x_xxxxxz[k];

                g_x_z_xx_xxxyy[k] = -g_0_z_x_xxxyy[k] - g_x_z_x_xxxyy[k] * ab_x + g_x_z_x_xxxxyy[k];

                g_x_z_xx_xxxyz[k] = -g_0_z_x_xxxyz[k] - g_x_z_x_xxxyz[k] * ab_x + g_x_z_x_xxxxyz[k];

                g_x_z_xx_xxxzz[k] = -g_0_z_x_xxxzz[k] - g_x_z_x_xxxzz[k] * ab_x + g_x_z_x_xxxxzz[k];

                g_x_z_xx_xxyyy[k] = -g_0_z_x_xxyyy[k] - g_x_z_x_xxyyy[k] * ab_x + g_x_z_x_xxxyyy[k];

                g_x_z_xx_xxyyz[k] = -g_0_z_x_xxyyz[k] - g_x_z_x_xxyyz[k] * ab_x + g_x_z_x_xxxyyz[k];

                g_x_z_xx_xxyzz[k] = -g_0_z_x_xxyzz[k] - g_x_z_x_xxyzz[k] * ab_x + g_x_z_x_xxxyzz[k];

                g_x_z_xx_xxzzz[k] = -g_0_z_x_xxzzz[k] - g_x_z_x_xxzzz[k] * ab_x + g_x_z_x_xxxzzz[k];

                g_x_z_xx_xyyyy[k] = -g_0_z_x_xyyyy[k] - g_x_z_x_xyyyy[k] * ab_x + g_x_z_x_xxyyyy[k];

                g_x_z_xx_xyyyz[k] = -g_0_z_x_xyyyz[k] - g_x_z_x_xyyyz[k] * ab_x + g_x_z_x_xxyyyz[k];

                g_x_z_xx_xyyzz[k] = -g_0_z_x_xyyzz[k] - g_x_z_x_xyyzz[k] * ab_x + g_x_z_x_xxyyzz[k];

                g_x_z_xx_xyzzz[k] = -g_0_z_x_xyzzz[k] - g_x_z_x_xyzzz[k] * ab_x + g_x_z_x_xxyzzz[k];

                g_x_z_xx_xzzzz[k] = -g_0_z_x_xzzzz[k] - g_x_z_x_xzzzz[k] * ab_x + g_x_z_x_xxzzzz[k];

                g_x_z_xx_yyyyy[k] = -g_0_z_x_yyyyy[k] - g_x_z_x_yyyyy[k] * ab_x + g_x_z_x_xyyyyy[k];

                g_x_z_xx_yyyyz[k] = -g_0_z_x_yyyyz[k] - g_x_z_x_yyyyz[k] * ab_x + g_x_z_x_xyyyyz[k];

                g_x_z_xx_yyyzz[k] = -g_0_z_x_yyyzz[k] - g_x_z_x_yyyzz[k] * ab_x + g_x_z_x_xyyyzz[k];

                g_x_z_xx_yyzzz[k] = -g_0_z_x_yyzzz[k] - g_x_z_x_yyzzz[k] * ab_x + g_x_z_x_xyyzzz[k];

                g_x_z_xx_yzzzz[k] = -g_0_z_x_yzzzz[k] - g_x_z_x_yzzzz[k] * ab_x + g_x_z_x_xyzzzz[k];

                g_x_z_xx_zzzzz[k] = -g_0_z_x_zzzzz[k] - g_x_z_x_zzzzz[k] * ab_x + g_x_z_x_xzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_x_z_xy_xxxxx = cbuffer.data(dh_geom_11_off + 273 * ccomps * dcomps);

            auto g_x_z_xy_xxxxy = cbuffer.data(dh_geom_11_off + 274 * ccomps * dcomps);

            auto g_x_z_xy_xxxxz = cbuffer.data(dh_geom_11_off + 275 * ccomps * dcomps);

            auto g_x_z_xy_xxxyy = cbuffer.data(dh_geom_11_off + 276 * ccomps * dcomps);

            auto g_x_z_xy_xxxyz = cbuffer.data(dh_geom_11_off + 277 * ccomps * dcomps);

            auto g_x_z_xy_xxxzz = cbuffer.data(dh_geom_11_off + 278 * ccomps * dcomps);

            auto g_x_z_xy_xxyyy = cbuffer.data(dh_geom_11_off + 279 * ccomps * dcomps);

            auto g_x_z_xy_xxyyz = cbuffer.data(dh_geom_11_off + 280 * ccomps * dcomps);

            auto g_x_z_xy_xxyzz = cbuffer.data(dh_geom_11_off + 281 * ccomps * dcomps);

            auto g_x_z_xy_xxzzz = cbuffer.data(dh_geom_11_off + 282 * ccomps * dcomps);

            auto g_x_z_xy_xyyyy = cbuffer.data(dh_geom_11_off + 283 * ccomps * dcomps);

            auto g_x_z_xy_xyyyz = cbuffer.data(dh_geom_11_off + 284 * ccomps * dcomps);

            auto g_x_z_xy_xyyzz = cbuffer.data(dh_geom_11_off + 285 * ccomps * dcomps);

            auto g_x_z_xy_xyzzz = cbuffer.data(dh_geom_11_off + 286 * ccomps * dcomps);

            auto g_x_z_xy_xzzzz = cbuffer.data(dh_geom_11_off + 287 * ccomps * dcomps);

            auto g_x_z_xy_yyyyy = cbuffer.data(dh_geom_11_off + 288 * ccomps * dcomps);

            auto g_x_z_xy_yyyyz = cbuffer.data(dh_geom_11_off + 289 * ccomps * dcomps);

            auto g_x_z_xy_yyyzz = cbuffer.data(dh_geom_11_off + 290 * ccomps * dcomps);

            auto g_x_z_xy_yyzzz = cbuffer.data(dh_geom_11_off + 291 * ccomps * dcomps);

            auto g_x_z_xy_yzzzz = cbuffer.data(dh_geom_11_off + 292 * ccomps * dcomps);

            auto g_x_z_xy_zzzzz = cbuffer.data(dh_geom_11_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_x_xxxxx, g_x_z_x_xxxxxy, g_x_z_x_xxxxy, g_x_z_x_xxxxyy, g_x_z_x_xxxxyz, g_x_z_x_xxxxz, g_x_z_x_xxxyy, g_x_z_x_xxxyyy, g_x_z_x_xxxyyz, g_x_z_x_xxxyz, g_x_z_x_xxxyzz, g_x_z_x_xxxzz, g_x_z_x_xxyyy, g_x_z_x_xxyyyy, g_x_z_x_xxyyyz, g_x_z_x_xxyyz, g_x_z_x_xxyyzz, g_x_z_x_xxyzz, g_x_z_x_xxyzzz, g_x_z_x_xxzzz, g_x_z_x_xyyyy, g_x_z_x_xyyyyy, g_x_z_x_xyyyyz, g_x_z_x_xyyyz, g_x_z_x_xyyyzz, g_x_z_x_xyyzz, g_x_z_x_xyyzzz, g_x_z_x_xyzzz, g_x_z_x_xyzzzz, g_x_z_x_xzzzz, g_x_z_x_yyyyy, g_x_z_x_yyyyyy, g_x_z_x_yyyyyz, g_x_z_x_yyyyz, g_x_z_x_yyyyzz, g_x_z_x_yyyzz, g_x_z_x_yyyzzz, g_x_z_x_yyzzz, g_x_z_x_yyzzzz, g_x_z_x_yzzzz, g_x_z_x_yzzzzz, g_x_z_x_zzzzz, g_x_z_xy_xxxxx, g_x_z_xy_xxxxy, g_x_z_xy_xxxxz, g_x_z_xy_xxxyy, g_x_z_xy_xxxyz, g_x_z_xy_xxxzz, g_x_z_xy_xxyyy, g_x_z_xy_xxyyz, g_x_z_xy_xxyzz, g_x_z_xy_xxzzz, g_x_z_xy_xyyyy, g_x_z_xy_xyyyz, g_x_z_xy_xyyzz, g_x_z_xy_xyzzz, g_x_z_xy_xzzzz, g_x_z_xy_yyyyy, g_x_z_xy_yyyyz, g_x_z_xy_yyyzz, g_x_z_xy_yyzzz, g_x_z_xy_yzzzz, g_x_z_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xy_xxxxx[k] = -g_x_z_x_xxxxx[k] * ab_y + g_x_z_x_xxxxxy[k];

                g_x_z_xy_xxxxy[k] = -g_x_z_x_xxxxy[k] * ab_y + g_x_z_x_xxxxyy[k];

                g_x_z_xy_xxxxz[k] = -g_x_z_x_xxxxz[k] * ab_y + g_x_z_x_xxxxyz[k];

                g_x_z_xy_xxxyy[k] = -g_x_z_x_xxxyy[k] * ab_y + g_x_z_x_xxxyyy[k];

                g_x_z_xy_xxxyz[k] = -g_x_z_x_xxxyz[k] * ab_y + g_x_z_x_xxxyyz[k];

                g_x_z_xy_xxxzz[k] = -g_x_z_x_xxxzz[k] * ab_y + g_x_z_x_xxxyzz[k];

                g_x_z_xy_xxyyy[k] = -g_x_z_x_xxyyy[k] * ab_y + g_x_z_x_xxyyyy[k];

                g_x_z_xy_xxyyz[k] = -g_x_z_x_xxyyz[k] * ab_y + g_x_z_x_xxyyyz[k];

                g_x_z_xy_xxyzz[k] = -g_x_z_x_xxyzz[k] * ab_y + g_x_z_x_xxyyzz[k];

                g_x_z_xy_xxzzz[k] = -g_x_z_x_xxzzz[k] * ab_y + g_x_z_x_xxyzzz[k];

                g_x_z_xy_xyyyy[k] = -g_x_z_x_xyyyy[k] * ab_y + g_x_z_x_xyyyyy[k];

                g_x_z_xy_xyyyz[k] = -g_x_z_x_xyyyz[k] * ab_y + g_x_z_x_xyyyyz[k];

                g_x_z_xy_xyyzz[k] = -g_x_z_x_xyyzz[k] * ab_y + g_x_z_x_xyyyzz[k];

                g_x_z_xy_xyzzz[k] = -g_x_z_x_xyzzz[k] * ab_y + g_x_z_x_xyyzzz[k];

                g_x_z_xy_xzzzz[k] = -g_x_z_x_xzzzz[k] * ab_y + g_x_z_x_xyzzzz[k];

                g_x_z_xy_yyyyy[k] = -g_x_z_x_yyyyy[k] * ab_y + g_x_z_x_yyyyyy[k];

                g_x_z_xy_yyyyz[k] = -g_x_z_x_yyyyz[k] * ab_y + g_x_z_x_yyyyyz[k];

                g_x_z_xy_yyyzz[k] = -g_x_z_x_yyyzz[k] * ab_y + g_x_z_x_yyyyzz[k];

                g_x_z_xy_yyzzz[k] = -g_x_z_x_yyzzz[k] * ab_y + g_x_z_x_yyyzzz[k];

                g_x_z_xy_yzzzz[k] = -g_x_z_x_yzzzz[k] * ab_y + g_x_z_x_yyzzzz[k];

                g_x_z_xy_zzzzz[k] = -g_x_z_x_zzzzz[k] * ab_y + g_x_z_x_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_x_z_xz_xxxxx = cbuffer.data(dh_geom_11_off + 294 * ccomps * dcomps);

            auto g_x_z_xz_xxxxy = cbuffer.data(dh_geom_11_off + 295 * ccomps * dcomps);

            auto g_x_z_xz_xxxxz = cbuffer.data(dh_geom_11_off + 296 * ccomps * dcomps);

            auto g_x_z_xz_xxxyy = cbuffer.data(dh_geom_11_off + 297 * ccomps * dcomps);

            auto g_x_z_xz_xxxyz = cbuffer.data(dh_geom_11_off + 298 * ccomps * dcomps);

            auto g_x_z_xz_xxxzz = cbuffer.data(dh_geom_11_off + 299 * ccomps * dcomps);

            auto g_x_z_xz_xxyyy = cbuffer.data(dh_geom_11_off + 300 * ccomps * dcomps);

            auto g_x_z_xz_xxyyz = cbuffer.data(dh_geom_11_off + 301 * ccomps * dcomps);

            auto g_x_z_xz_xxyzz = cbuffer.data(dh_geom_11_off + 302 * ccomps * dcomps);

            auto g_x_z_xz_xxzzz = cbuffer.data(dh_geom_11_off + 303 * ccomps * dcomps);

            auto g_x_z_xz_xyyyy = cbuffer.data(dh_geom_11_off + 304 * ccomps * dcomps);

            auto g_x_z_xz_xyyyz = cbuffer.data(dh_geom_11_off + 305 * ccomps * dcomps);

            auto g_x_z_xz_xyyzz = cbuffer.data(dh_geom_11_off + 306 * ccomps * dcomps);

            auto g_x_z_xz_xyzzz = cbuffer.data(dh_geom_11_off + 307 * ccomps * dcomps);

            auto g_x_z_xz_xzzzz = cbuffer.data(dh_geom_11_off + 308 * ccomps * dcomps);

            auto g_x_z_xz_yyyyy = cbuffer.data(dh_geom_11_off + 309 * ccomps * dcomps);

            auto g_x_z_xz_yyyyz = cbuffer.data(dh_geom_11_off + 310 * ccomps * dcomps);

            auto g_x_z_xz_yyyzz = cbuffer.data(dh_geom_11_off + 311 * ccomps * dcomps);

            auto g_x_z_xz_yyzzz = cbuffer.data(dh_geom_11_off + 312 * ccomps * dcomps);

            auto g_x_z_xz_yzzzz = cbuffer.data(dh_geom_11_off + 313 * ccomps * dcomps);

            auto g_x_z_xz_zzzzz = cbuffer.data(dh_geom_11_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxxx, g_0_z_z_xxxxy, g_0_z_z_xxxxz, g_0_z_z_xxxyy, g_0_z_z_xxxyz, g_0_z_z_xxxzz, g_0_z_z_xxyyy, g_0_z_z_xxyyz, g_0_z_z_xxyzz, g_0_z_z_xxzzz, g_0_z_z_xyyyy, g_0_z_z_xyyyz, g_0_z_z_xyyzz, g_0_z_z_xyzzz, g_0_z_z_xzzzz, g_0_z_z_yyyyy, g_0_z_z_yyyyz, g_0_z_z_yyyzz, g_0_z_z_yyzzz, g_0_z_z_yzzzz, g_0_z_z_zzzzz, g_x_z_xz_xxxxx, g_x_z_xz_xxxxy, g_x_z_xz_xxxxz, g_x_z_xz_xxxyy, g_x_z_xz_xxxyz, g_x_z_xz_xxxzz, g_x_z_xz_xxyyy, g_x_z_xz_xxyyz, g_x_z_xz_xxyzz, g_x_z_xz_xxzzz, g_x_z_xz_xyyyy, g_x_z_xz_xyyyz, g_x_z_xz_xyyzz, g_x_z_xz_xyzzz, g_x_z_xz_xzzzz, g_x_z_xz_yyyyy, g_x_z_xz_yyyyz, g_x_z_xz_yyyzz, g_x_z_xz_yyzzz, g_x_z_xz_yzzzz, g_x_z_xz_zzzzz, g_x_z_z_xxxxx, g_x_z_z_xxxxxx, g_x_z_z_xxxxxy, g_x_z_z_xxxxxz, g_x_z_z_xxxxy, g_x_z_z_xxxxyy, g_x_z_z_xxxxyz, g_x_z_z_xxxxz, g_x_z_z_xxxxzz, g_x_z_z_xxxyy, g_x_z_z_xxxyyy, g_x_z_z_xxxyyz, g_x_z_z_xxxyz, g_x_z_z_xxxyzz, g_x_z_z_xxxzz, g_x_z_z_xxxzzz, g_x_z_z_xxyyy, g_x_z_z_xxyyyy, g_x_z_z_xxyyyz, g_x_z_z_xxyyz, g_x_z_z_xxyyzz, g_x_z_z_xxyzz, g_x_z_z_xxyzzz, g_x_z_z_xxzzz, g_x_z_z_xxzzzz, g_x_z_z_xyyyy, g_x_z_z_xyyyyy, g_x_z_z_xyyyyz, g_x_z_z_xyyyz, g_x_z_z_xyyyzz, g_x_z_z_xyyzz, g_x_z_z_xyyzzz, g_x_z_z_xyzzz, g_x_z_z_xyzzzz, g_x_z_z_xzzzz, g_x_z_z_xzzzzz, g_x_z_z_yyyyy, g_x_z_z_yyyyz, g_x_z_z_yyyzz, g_x_z_z_yyzzz, g_x_z_z_yzzzz, g_x_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_xz_xxxxx[k] = -g_0_z_z_xxxxx[k] - g_x_z_z_xxxxx[k] * ab_x + g_x_z_z_xxxxxx[k];

                g_x_z_xz_xxxxy[k] = -g_0_z_z_xxxxy[k] - g_x_z_z_xxxxy[k] * ab_x + g_x_z_z_xxxxxy[k];

                g_x_z_xz_xxxxz[k] = -g_0_z_z_xxxxz[k] - g_x_z_z_xxxxz[k] * ab_x + g_x_z_z_xxxxxz[k];

                g_x_z_xz_xxxyy[k] = -g_0_z_z_xxxyy[k] - g_x_z_z_xxxyy[k] * ab_x + g_x_z_z_xxxxyy[k];

                g_x_z_xz_xxxyz[k] = -g_0_z_z_xxxyz[k] - g_x_z_z_xxxyz[k] * ab_x + g_x_z_z_xxxxyz[k];

                g_x_z_xz_xxxzz[k] = -g_0_z_z_xxxzz[k] - g_x_z_z_xxxzz[k] * ab_x + g_x_z_z_xxxxzz[k];

                g_x_z_xz_xxyyy[k] = -g_0_z_z_xxyyy[k] - g_x_z_z_xxyyy[k] * ab_x + g_x_z_z_xxxyyy[k];

                g_x_z_xz_xxyyz[k] = -g_0_z_z_xxyyz[k] - g_x_z_z_xxyyz[k] * ab_x + g_x_z_z_xxxyyz[k];

                g_x_z_xz_xxyzz[k] = -g_0_z_z_xxyzz[k] - g_x_z_z_xxyzz[k] * ab_x + g_x_z_z_xxxyzz[k];

                g_x_z_xz_xxzzz[k] = -g_0_z_z_xxzzz[k] - g_x_z_z_xxzzz[k] * ab_x + g_x_z_z_xxxzzz[k];

                g_x_z_xz_xyyyy[k] = -g_0_z_z_xyyyy[k] - g_x_z_z_xyyyy[k] * ab_x + g_x_z_z_xxyyyy[k];

                g_x_z_xz_xyyyz[k] = -g_0_z_z_xyyyz[k] - g_x_z_z_xyyyz[k] * ab_x + g_x_z_z_xxyyyz[k];

                g_x_z_xz_xyyzz[k] = -g_0_z_z_xyyzz[k] - g_x_z_z_xyyzz[k] * ab_x + g_x_z_z_xxyyzz[k];

                g_x_z_xz_xyzzz[k] = -g_0_z_z_xyzzz[k] - g_x_z_z_xyzzz[k] * ab_x + g_x_z_z_xxyzzz[k];

                g_x_z_xz_xzzzz[k] = -g_0_z_z_xzzzz[k] - g_x_z_z_xzzzz[k] * ab_x + g_x_z_z_xxzzzz[k];

                g_x_z_xz_yyyyy[k] = -g_0_z_z_yyyyy[k] - g_x_z_z_yyyyy[k] * ab_x + g_x_z_z_xyyyyy[k];

                g_x_z_xz_yyyyz[k] = -g_0_z_z_yyyyz[k] - g_x_z_z_yyyyz[k] * ab_x + g_x_z_z_xyyyyz[k];

                g_x_z_xz_yyyzz[k] = -g_0_z_z_yyyzz[k] - g_x_z_z_yyyzz[k] * ab_x + g_x_z_z_xyyyzz[k];

                g_x_z_xz_yyzzz[k] = -g_0_z_z_yyzzz[k] - g_x_z_z_yyzzz[k] * ab_x + g_x_z_z_xyyzzz[k];

                g_x_z_xz_yzzzz[k] = -g_0_z_z_yzzzz[k] - g_x_z_z_yzzzz[k] * ab_x + g_x_z_z_xyzzzz[k];

                g_x_z_xz_zzzzz[k] = -g_0_z_z_zzzzz[k] - g_x_z_z_zzzzz[k] * ab_x + g_x_z_z_xzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_x_z_yy_xxxxx = cbuffer.data(dh_geom_11_off + 315 * ccomps * dcomps);

            auto g_x_z_yy_xxxxy = cbuffer.data(dh_geom_11_off + 316 * ccomps * dcomps);

            auto g_x_z_yy_xxxxz = cbuffer.data(dh_geom_11_off + 317 * ccomps * dcomps);

            auto g_x_z_yy_xxxyy = cbuffer.data(dh_geom_11_off + 318 * ccomps * dcomps);

            auto g_x_z_yy_xxxyz = cbuffer.data(dh_geom_11_off + 319 * ccomps * dcomps);

            auto g_x_z_yy_xxxzz = cbuffer.data(dh_geom_11_off + 320 * ccomps * dcomps);

            auto g_x_z_yy_xxyyy = cbuffer.data(dh_geom_11_off + 321 * ccomps * dcomps);

            auto g_x_z_yy_xxyyz = cbuffer.data(dh_geom_11_off + 322 * ccomps * dcomps);

            auto g_x_z_yy_xxyzz = cbuffer.data(dh_geom_11_off + 323 * ccomps * dcomps);

            auto g_x_z_yy_xxzzz = cbuffer.data(dh_geom_11_off + 324 * ccomps * dcomps);

            auto g_x_z_yy_xyyyy = cbuffer.data(dh_geom_11_off + 325 * ccomps * dcomps);

            auto g_x_z_yy_xyyyz = cbuffer.data(dh_geom_11_off + 326 * ccomps * dcomps);

            auto g_x_z_yy_xyyzz = cbuffer.data(dh_geom_11_off + 327 * ccomps * dcomps);

            auto g_x_z_yy_xyzzz = cbuffer.data(dh_geom_11_off + 328 * ccomps * dcomps);

            auto g_x_z_yy_xzzzz = cbuffer.data(dh_geom_11_off + 329 * ccomps * dcomps);

            auto g_x_z_yy_yyyyy = cbuffer.data(dh_geom_11_off + 330 * ccomps * dcomps);

            auto g_x_z_yy_yyyyz = cbuffer.data(dh_geom_11_off + 331 * ccomps * dcomps);

            auto g_x_z_yy_yyyzz = cbuffer.data(dh_geom_11_off + 332 * ccomps * dcomps);

            auto g_x_z_yy_yyzzz = cbuffer.data(dh_geom_11_off + 333 * ccomps * dcomps);

            auto g_x_z_yy_yzzzz = cbuffer.data(dh_geom_11_off + 334 * ccomps * dcomps);

            auto g_x_z_yy_zzzzz = cbuffer.data(dh_geom_11_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_y_xxxxx, g_x_z_y_xxxxxy, g_x_z_y_xxxxy, g_x_z_y_xxxxyy, g_x_z_y_xxxxyz, g_x_z_y_xxxxz, g_x_z_y_xxxyy, g_x_z_y_xxxyyy, g_x_z_y_xxxyyz, g_x_z_y_xxxyz, g_x_z_y_xxxyzz, g_x_z_y_xxxzz, g_x_z_y_xxyyy, g_x_z_y_xxyyyy, g_x_z_y_xxyyyz, g_x_z_y_xxyyz, g_x_z_y_xxyyzz, g_x_z_y_xxyzz, g_x_z_y_xxyzzz, g_x_z_y_xxzzz, g_x_z_y_xyyyy, g_x_z_y_xyyyyy, g_x_z_y_xyyyyz, g_x_z_y_xyyyz, g_x_z_y_xyyyzz, g_x_z_y_xyyzz, g_x_z_y_xyyzzz, g_x_z_y_xyzzz, g_x_z_y_xyzzzz, g_x_z_y_xzzzz, g_x_z_y_yyyyy, g_x_z_y_yyyyyy, g_x_z_y_yyyyyz, g_x_z_y_yyyyz, g_x_z_y_yyyyzz, g_x_z_y_yyyzz, g_x_z_y_yyyzzz, g_x_z_y_yyzzz, g_x_z_y_yyzzzz, g_x_z_y_yzzzz, g_x_z_y_yzzzzz, g_x_z_y_zzzzz, g_x_z_yy_xxxxx, g_x_z_yy_xxxxy, g_x_z_yy_xxxxz, g_x_z_yy_xxxyy, g_x_z_yy_xxxyz, g_x_z_yy_xxxzz, g_x_z_yy_xxyyy, g_x_z_yy_xxyyz, g_x_z_yy_xxyzz, g_x_z_yy_xxzzz, g_x_z_yy_xyyyy, g_x_z_yy_xyyyz, g_x_z_yy_xyyzz, g_x_z_yy_xyzzz, g_x_z_yy_xzzzz, g_x_z_yy_yyyyy, g_x_z_yy_yyyyz, g_x_z_yy_yyyzz, g_x_z_yy_yyzzz, g_x_z_yy_yzzzz, g_x_z_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yy_xxxxx[k] = -g_x_z_y_xxxxx[k] * ab_y + g_x_z_y_xxxxxy[k];

                g_x_z_yy_xxxxy[k] = -g_x_z_y_xxxxy[k] * ab_y + g_x_z_y_xxxxyy[k];

                g_x_z_yy_xxxxz[k] = -g_x_z_y_xxxxz[k] * ab_y + g_x_z_y_xxxxyz[k];

                g_x_z_yy_xxxyy[k] = -g_x_z_y_xxxyy[k] * ab_y + g_x_z_y_xxxyyy[k];

                g_x_z_yy_xxxyz[k] = -g_x_z_y_xxxyz[k] * ab_y + g_x_z_y_xxxyyz[k];

                g_x_z_yy_xxxzz[k] = -g_x_z_y_xxxzz[k] * ab_y + g_x_z_y_xxxyzz[k];

                g_x_z_yy_xxyyy[k] = -g_x_z_y_xxyyy[k] * ab_y + g_x_z_y_xxyyyy[k];

                g_x_z_yy_xxyyz[k] = -g_x_z_y_xxyyz[k] * ab_y + g_x_z_y_xxyyyz[k];

                g_x_z_yy_xxyzz[k] = -g_x_z_y_xxyzz[k] * ab_y + g_x_z_y_xxyyzz[k];

                g_x_z_yy_xxzzz[k] = -g_x_z_y_xxzzz[k] * ab_y + g_x_z_y_xxyzzz[k];

                g_x_z_yy_xyyyy[k] = -g_x_z_y_xyyyy[k] * ab_y + g_x_z_y_xyyyyy[k];

                g_x_z_yy_xyyyz[k] = -g_x_z_y_xyyyz[k] * ab_y + g_x_z_y_xyyyyz[k];

                g_x_z_yy_xyyzz[k] = -g_x_z_y_xyyzz[k] * ab_y + g_x_z_y_xyyyzz[k];

                g_x_z_yy_xyzzz[k] = -g_x_z_y_xyzzz[k] * ab_y + g_x_z_y_xyyzzz[k];

                g_x_z_yy_xzzzz[k] = -g_x_z_y_xzzzz[k] * ab_y + g_x_z_y_xyzzzz[k];

                g_x_z_yy_yyyyy[k] = -g_x_z_y_yyyyy[k] * ab_y + g_x_z_y_yyyyyy[k];

                g_x_z_yy_yyyyz[k] = -g_x_z_y_yyyyz[k] * ab_y + g_x_z_y_yyyyyz[k];

                g_x_z_yy_yyyzz[k] = -g_x_z_y_yyyzz[k] * ab_y + g_x_z_y_yyyyzz[k];

                g_x_z_yy_yyzzz[k] = -g_x_z_y_yyzzz[k] * ab_y + g_x_z_y_yyyzzz[k];

                g_x_z_yy_yzzzz[k] = -g_x_z_y_yzzzz[k] * ab_y + g_x_z_y_yyzzzz[k];

                g_x_z_yy_zzzzz[k] = -g_x_z_y_zzzzz[k] * ab_y + g_x_z_y_yzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_x_z_yz_xxxxx = cbuffer.data(dh_geom_11_off + 336 * ccomps * dcomps);

            auto g_x_z_yz_xxxxy = cbuffer.data(dh_geom_11_off + 337 * ccomps * dcomps);

            auto g_x_z_yz_xxxxz = cbuffer.data(dh_geom_11_off + 338 * ccomps * dcomps);

            auto g_x_z_yz_xxxyy = cbuffer.data(dh_geom_11_off + 339 * ccomps * dcomps);

            auto g_x_z_yz_xxxyz = cbuffer.data(dh_geom_11_off + 340 * ccomps * dcomps);

            auto g_x_z_yz_xxxzz = cbuffer.data(dh_geom_11_off + 341 * ccomps * dcomps);

            auto g_x_z_yz_xxyyy = cbuffer.data(dh_geom_11_off + 342 * ccomps * dcomps);

            auto g_x_z_yz_xxyyz = cbuffer.data(dh_geom_11_off + 343 * ccomps * dcomps);

            auto g_x_z_yz_xxyzz = cbuffer.data(dh_geom_11_off + 344 * ccomps * dcomps);

            auto g_x_z_yz_xxzzz = cbuffer.data(dh_geom_11_off + 345 * ccomps * dcomps);

            auto g_x_z_yz_xyyyy = cbuffer.data(dh_geom_11_off + 346 * ccomps * dcomps);

            auto g_x_z_yz_xyyyz = cbuffer.data(dh_geom_11_off + 347 * ccomps * dcomps);

            auto g_x_z_yz_xyyzz = cbuffer.data(dh_geom_11_off + 348 * ccomps * dcomps);

            auto g_x_z_yz_xyzzz = cbuffer.data(dh_geom_11_off + 349 * ccomps * dcomps);

            auto g_x_z_yz_xzzzz = cbuffer.data(dh_geom_11_off + 350 * ccomps * dcomps);

            auto g_x_z_yz_yyyyy = cbuffer.data(dh_geom_11_off + 351 * ccomps * dcomps);

            auto g_x_z_yz_yyyyz = cbuffer.data(dh_geom_11_off + 352 * ccomps * dcomps);

            auto g_x_z_yz_yyyzz = cbuffer.data(dh_geom_11_off + 353 * ccomps * dcomps);

            auto g_x_z_yz_yyzzz = cbuffer.data(dh_geom_11_off + 354 * ccomps * dcomps);

            auto g_x_z_yz_yzzzz = cbuffer.data(dh_geom_11_off + 355 * ccomps * dcomps);

            auto g_x_z_yz_zzzzz = cbuffer.data(dh_geom_11_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_z_yz_xxxxx, g_x_z_yz_xxxxy, g_x_z_yz_xxxxz, g_x_z_yz_xxxyy, g_x_z_yz_xxxyz, g_x_z_yz_xxxzz, g_x_z_yz_xxyyy, g_x_z_yz_xxyyz, g_x_z_yz_xxyzz, g_x_z_yz_xxzzz, g_x_z_yz_xyyyy, g_x_z_yz_xyyyz, g_x_z_yz_xyyzz, g_x_z_yz_xyzzz, g_x_z_yz_xzzzz, g_x_z_yz_yyyyy, g_x_z_yz_yyyyz, g_x_z_yz_yyyzz, g_x_z_yz_yyzzz, g_x_z_yz_yzzzz, g_x_z_yz_zzzzz, g_x_z_z_xxxxx, g_x_z_z_xxxxxy, g_x_z_z_xxxxy, g_x_z_z_xxxxyy, g_x_z_z_xxxxyz, g_x_z_z_xxxxz, g_x_z_z_xxxyy, g_x_z_z_xxxyyy, g_x_z_z_xxxyyz, g_x_z_z_xxxyz, g_x_z_z_xxxyzz, g_x_z_z_xxxzz, g_x_z_z_xxyyy, g_x_z_z_xxyyyy, g_x_z_z_xxyyyz, g_x_z_z_xxyyz, g_x_z_z_xxyyzz, g_x_z_z_xxyzz, g_x_z_z_xxyzzz, g_x_z_z_xxzzz, g_x_z_z_xyyyy, g_x_z_z_xyyyyy, g_x_z_z_xyyyyz, g_x_z_z_xyyyz, g_x_z_z_xyyyzz, g_x_z_z_xyyzz, g_x_z_z_xyyzzz, g_x_z_z_xyzzz, g_x_z_z_xyzzzz, g_x_z_z_xzzzz, g_x_z_z_yyyyy, g_x_z_z_yyyyyy, g_x_z_z_yyyyyz, g_x_z_z_yyyyz, g_x_z_z_yyyyzz, g_x_z_z_yyyzz, g_x_z_z_yyyzzz, g_x_z_z_yyzzz, g_x_z_z_yyzzzz, g_x_z_z_yzzzz, g_x_z_z_yzzzzz, g_x_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_yz_xxxxx[k] = -g_x_z_z_xxxxx[k] * ab_y + g_x_z_z_xxxxxy[k];

                g_x_z_yz_xxxxy[k] = -g_x_z_z_xxxxy[k] * ab_y + g_x_z_z_xxxxyy[k];

                g_x_z_yz_xxxxz[k] = -g_x_z_z_xxxxz[k] * ab_y + g_x_z_z_xxxxyz[k];

                g_x_z_yz_xxxyy[k] = -g_x_z_z_xxxyy[k] * ab_y + g_x_z_z_xxxyyy[k];

                g_x_z_yz_xxxyz[k] = -g_x_z_z_xxxyz[k] * ab_y + g_x_z_z_xxxyyz[k];

                g_x_z_yz_xxxzz[k] = -g_x_z_z_xxxzz[k] * ab_y + g_x_z_z_xxxyzz[k];

                g_x_z_yz_xxyyy[k] = -g_x_z_z_xxyyy[k] * ab_y + g_x_z_z_xxyyyy[k];

                g_x_z_yz_xxyyz[k] = -g_x_z_z_xxyyz[k] * ab_y + g_x_z_z_xxyyyz[k];

                g_x_z_yz_xxyzz[k] = -g_x_z_z_xxyzz[k] * ab_y + g_x_z_z_xxyyzz[k];

                g_x_z_yz_xxzzz[k] = -g_x_z_z_xxzzz[k] * ab_y + g_x_z_z_xxyzzz[k];

                g_x_z_yz_xyyyy[k] = -g_x_z_z_xyyyy[k] * ab_y + g_x_z_z_xyyyyy[k];

                g_x_z_yz_xyyyz[k] = -g_x_z_z_xyyyz[k] * ab_y + g_x_z_z_xyyyyz[k];

                g_x_z_yz_xyyzz[k] = -g_x_z_z_xyyzz[k] * ab_y + g_x_z_z_xyyyzz[k];

                g_x_z_yz_xyzzz[k] = -g_x_z_z_xyzzz[k] * ab_y + g_x_z_z_xyyzzz[k];

                g_x_z_yz_xzzzz[k] = -g_x_z_z_xzzzz[k] * ab_y + g_x_z_z_xyzzzz[k];

                g_x_z_yz_yyyyy[k] = -g_x_z_z_yyyyy[k] * ab_y + g_x_z_z_yyyyyy[k];

                g_x_z_yz_yyyyz[k] = -g_x_z_z_yyyyz[k] * ab_y + g_x_z_z_yyyyyz[k];

                g_x_z_yz_yyyzz[k] = -g_x_z_z_yyyzz[k] * ab_y + g_x_z_z_yyyyzz[k];

                g_x_z_yz_yyzzz[k] = -g_x_z_z_yyzzz[k] * ab_y + g_x_z_z_yyyzzz[k];

                g_x_z_yz_yzzzz[k] = -g_x_z_z_yzzzz[k] * ab_y + g_x_z_z_yyzzzz[k];

                g_x_z_yz_zzzzz[k] = -g_x_z_z_zzzzz[k] * ab_y + g_x_z_z_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_x_z_zz_xxxxx = cbuffer.data(dh_geom_11_off + 357 * ccomps * dcomps);

            auto g_x_z_zz_xxxxy = cbuffer.data(dh_geom_11_off + 358 * ccomps * dcomps);

            auto g_x_z_zz_xxxxz = cbuffer.data(dh_geom_11_off + 359 * ccomps * dcomps);

            auto g_x_z_zz_xxxyy = cbuffer.data(dh_geom_11_off + 360 * ccomps * dcomps);

            auto g_x_z_zz_xxxyz = cbuffer.data(dh_geom_11_off + 361 * ccomps * dcomps);

            auto g_x_z_zz_xxxzz = cbuffer.data(dh_geom_11_off + 362 * ccomps * dcomps);

            auto g_x_z_zz_xxyyy = cbuffer.data(dh_geom_11_off + 363 * ccomps * dcomps);

            auto g_x_z_zz_xxyyz = cbuffer.data(dh_geom_11_off + 364 * ccomps * dcomps);

            auto g_x_z_zz_xxyzz = cbuffer.data(dh_geom_11_off + 365 * ccomps * dcomps);

            auto g_x_z_zz_xxzzz = cbuffer.data(dh_geom_11_off + 366 * ccomps * dcomps);

            auto g_x_z_zz_xyyyy = cbuffer.data(dh_geom_11_off + 367 * ccomps * dcomps);

            auto g_x_z_zz_xyyyz = cbuffer.data(dh_geom_11_off + 368 * ccomps * dcomps);

            auto g_x_z_zz_xyyzz = cbuffer.data(dh_geom_11_off + 369 * ccomps * dcomps);

            auto g_x_z_zz_xyzzz = cbuffer.data(dh_geom_11_off + 370 * ccomps * dcomps);

            auto g_x_z_zz_xzzzz = cbuffer.data(dh_geom_11_off + 371 * ccomps * dcomps);

            auto g_x_z_zz_yyyyy = cbuffer.data(dh_geom_11_off + 372 * ccomps * dcomps);

            auto g_x_z_zz_yyyyz = cbuffer.data(dh_geom_11_off + 373 * ccomps * dcomps);

            auto g_x_z_zz_yyyzz = cbuffer.data(dh_geom_11_off + 374 * ccomps * dcomps);

            auto g_x_z_zz_yyzzz = cbuffer.data(dh_geom_11_off + 375 * ccomps * dcomps);

            auto g_x_z_zz_yzzzz = cbuffer.data(dh_geom_11_off + 376 * ccomps * dcomps);

            auto g_x_z_zz_zzzzz = cbuffer.data(dh_geom_11_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_z_xxxxx, g_x_0_z_xxxxy, g_x_0_z_xxxxz, g_x_0_z_xxxyy, g_x_0_z_xxxyz, g_x_0_z_xxxzz, g_x_0_z_xxyyy, g_x_0_z_xxyyz, g_x_0_z_xxyzz, g_x_0_z_xxzzz, g_x_0_z_xyyyy, g_x_0_z_xyyyz, g_x_0_z_xyyzz, g_x_0_z_xyzzz, g_x_0_z_xzzzz, g_x_0_z_yyyyy, g_x_0_z_yyyyz, g_x_0_z_yyyzz, g_x_0_z_yyzzz, g_x_0_z_yzzzz, g_x_0_z_zzzzz, g_x_z_z_xxxxx, g_x_z_z_xxxxxz, g_x_z_z_xxxxy, g_x_z_z_xxxxyz, g_x_z_z_xxxxz, g_x_z_z_xxxxzz, g_x_z_z_xxxyy, g_x_z_z_xxxyyz, g_x_z_z_xxxyz, g_x_z_z_xxxyzz, g_x_z_z_xxxzz, g_x_z_z_xxxzzz, g_x_z_z_xxyyy, g_x_z_z_xxyyyz, g_x_z_z_xxyyz, g_x_z_z_xxyyzz, g_x_z_z_xxyzz, g_x_z_z_xxyzzz, g_x_z_z_xxzzz, g_x_z_z_xxzzzz, g_x_z_z_xyyyy, g_x_z_z_xyyyyz, g_x_z_z_xyyyz, g_x_z_z_xyyyzz, g_x_z_z_xyyzz, g_x_z_z_xyyzzz, g_x_z_z_xyzzz, g_x_z_z_xyzzzz, g_x_z_z_xzzzz, g_x_z_z_xzzzzz, g_x_z_z_yyyyy, g_x_z_z_yyyyyz, g_x_z_z_yyyyz, g_x_z_z_yyyyzz, g_x_z_z_yyyzz, g_x_z_z_yyyzzz, g_x_z_z_yyzzz, g_x_z_z_yyzzzz, g_x_z_z_yzzzz, g_x_z_z_yzzzzz, g_x_z_z_zzzzz, g_x_z_z_zzzzzz, g_x_z_zz_xxxxx, g_x_z_zz_xxxxy, g_x_z_zz_xxxxz, g_x_z_zz_xxxyy, g_x_z_zz_xxxyz, g_x_z_zz_xxxzz, g_x_z_zz_xxyyy, g_x_z_zz_xxyyz, g_x_z_zz_xxyzz, g_x_z_zz_xxzzz, g_x_z_zz_xyyyy, g_x_z_zz_xyyyz, g_x_z_zz_xyyzz, g_x_z_zz_xyzzz, g_x_z_zz_xzzzz, g_x_z_zz_yyyyy, g_x_z_zz_yyyyz, g_x_z_zz_yyyzz, g_x_z_zz_yyzzz, g_x_z_zz_yzzzz, g_x_z_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_z_zz_xxxxx[k] = g_x_0_z_xxxxx[k] - g_x_z_z_xxxxx[k] * ab_z + g_x_z_z_xxxxxz[k];

                g_x_z_zz_xxxxy[k] = g_x_0_z_xxxxy[k] - g_x_z_z_xxxxy[k] * ab_z + g_x_z_z_xxxxyz[k];

                g_x_z_zz_xxxxz[k] = g_x_0_z_xxxxz[k] - g_x_z_z_xxxxz[k] * ab_z + g_x_z_z_xxxxzz[k];

                g_x_z_zz_xxxyy[k] = g_x_0_z_xxxyy[k] - g_x_z_z_xxxyy[k] * ab_z + g_x_z_z_xxxyyz[k];

                g_x_z_zz_xxxyz[k] = g_x_0_z_xxxyz[k] - g_x_z_z_xxxyz[k] * ab_z + g_x_z_z_xxxyzz[k];

                g_x_z_zz_xxxzz[k] = g_x_0_z_xxxzz[k] - g_x_z_z_xxxzz[k] * ab_z + g_x_z_z_xxxzzz[k];

                g_x_z_zz_xxyyy[k] = g_x_0_z_xxyyy[k] - g_x_z_z_xxyyy[k] * ab_z + g_x_z_z_xxyyyz[k];

                g_x_z_zz_xxyyz[k] = g_x_0_z_xxyyz[k] - g_x_z_z_xxyyz[k] * ab_z + g_x_z_z_xxyyzz[k];

                g_x_z_zz_xxyzz[k] = g_x_0_z_xxyzz[k] - g_x_z_z_xxyzz[k] * ab_z + g_x_z_z_xxyzzz[k];

                g_x_z_zz_xxzzz[k] = g_x_0_z_xxzzz[k] - g_x_z_z_xxzzz[k] * ab_z + g_x_z_z_xxzzzz[k];

                g_x_z_zz_xyyyy[k] = g_x_0_z_xyyyy[k] - g_x_z_z_xyyyy[k] * ab_z + g_x_z_z_xyyyyz[k];

                g_x_z_zz_xyyyz[k] = g_x_0_z_xyyyz[k] - g_x_z_z_xyyyz[k] * ab_z + g_x_z_z_xyyyzz[k];

                g_x_z_zz_xyyzz[k] = g_x_0_z_xyyzz[k] - g_x_z_z_xyyzz[k] * ab_z + g_x_z_z_xyyzzz[k];

                g_x_z_zz_xyzzz[k] = g_x_0_z_xyzzz[k] - g_x_z_z_xyzzz[k] * ab_z + g_x_z_z_xyzzzz[k];

                g_x_z_zz_xzzzz[k] = g_x_0_z_xzzzz[k] - g_x_z_z_xzzzz[k] * ab_z + g_x_z_z_xzzzzz[k];

                g_x_z_zz_yyyyy[k] = g_x_0_z_yyyyy[k] - g_x_z_z_yyyyy[k] * ab_z + g_x_z_z_yyyyyz[k];

                g_x_z_zz_yyyyz[k] = g_x_0_z_yyyyz[k] - g_x_z_z_yyyyz[k] * ab_z + g_x_z_z_yyyyzz[k];

                g_x_z_zz_yyyzz[k] = g_x_0_z_yyyzz[k] - g_x_z_z_yyyzz[k] * ab_z + g_x_z_z_yyyzzz[k];

                g_x_z_zz_yyzzz[k] = g_x_0_z_yyzzz[k] - g_x_z_z_yyzzz[k] * ab_z + g_x_z_z_yyzzzz[k];

                g_x_z_zz_yzzzz[k] = g_x_0_z_yzzzz[k] - g_x_z_z_yzzzz[k] * ab_z + g_x_z_z_yzzzzz[k];

                g_x_z_zz_zzzzz[k] = g_x_0_z_zzzzz[k] - g_x_z_z_zzzzz[k] * ab_z + g_x_z_z_zzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_y_x_xx_xxxxx = cbuffer.data(dh_geom_11_off + 378 * ccomps * dcomps);

            auto g_y_x_xx_xxxxy = cbuffer.data(dh_geom_11_off + 379 * ccomps * dcomps);

            auto g_y_x_xx_xxxxz = cbuffer.data(dh_geom_11_off + 380 * ccomps * dcomps);

            auto g_y_x_xx_xxxyy = cbuffer.data(dh_geom_11_off + 381 * ccomps * dcomps);

            auto g_y_x_xx_xxxyz = cbuffer.data(dh_geom_11_off + 382 * ccomps * dcomps);

            auto g_y_x_xx_xxxzz = cbuffer.data(dh_geom_11_off + 383 * ccomps * dcomps);

            auto g_y_x_xx_xxyyy = cbuffer.data(dh_geom_11_off + 384 * ccomps * dcomps);

            auto g_y_x_xx_xxyyz = cbuffer.data(dh_geom_11_off + 385 * ccomps * dcomps);

            auto g_y_x_xx_xxyzz = cbuffer.data(dh_geom_11_off + 386 * ccomps * dcomps);

            auto g_y_x_xx_xxzzz = cbuffer.data(dh_geom_11_off + 387 * ccomps * dcomps);

            auto g_y_x_xx_xyyyy = cbuffer.data(dh_geom_11_off + 388 * ccomps * dcomps);

            auto g_y_x_xx_xyyyz = cbuffer.data(dh_geom_11_off + 389 * ccomps * dcomps);

            auto g_y_x_xx_xyyzz = cbuffer.data(dh_geom_11_off + 390 * ccomps * dcomps);

            auto g_y_x_xx_xyzzz = cbuffer.data(dh_geom_11_off + 391 * ccomps * dcomps);

            auto g_y_x_xx_xzzzz = cbuffer.data(dh_geom_11_off + 392 * ccomps * dcomps);

            auto g_y_x_xx_yyyyy = cbuffer.data(dh_geom_11_off + 393 * ccomps * dcomps);

            auto g_y_x_xx_yyyyz = cbuffer.data(dh_geom_11_off + 394 * ccomps * dcomps);

            auto g_y_x_xx_yyyzz = cbuffer.data(dh_geom_11_off + 395 * ccomps * dcomps);

            auto g_y_x_xx_yyzzz = cbuffer.data(dh_geom_11_off + 396 * ccomps * dcomps);

            auto g_y_x_xx_yzzzz = cbuffer.data(dh_geom_11_off + 397 * ccomps * dcomps);

            auto g_y_x_xx_zzzzz = cbuffer.data(dh_geom_11_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_x_xxxxx, g_y_0_x_xxxxy, g_y_0_x_xxxxz, g_y_0_x_xxxyy, g_y_0_x_xxxyz, g_y_0_x_xxxzz, g_y_0_x_xxyyy, g_y_0_x_xxyyz, g_y_0_x_xxyzz, g_y_0_x_xxzzz, g_y_0_x_xyyyy, g_y_0_x_xyyyz, g_y_0_x_xyyzz, g_y_0_x_xyzzz, g_y_0_x_xzzzz, g_y_0_x_yyyyy, g_y_0_x_yyyyz, g_y_0_x_yyyzz, g_y_0_x_yyzzz, g_y_0_x_yzzzz, g_y_0_x_zzzzz, g_y_x_x_xxxxx, g_y_x_x_xxxxxx, g_y_x_x_xxxxxy, g_y_x_x_xxxxxz, g_y_x_x_xxxxy, g_y_x_x_xxxxyy, g_y_x_x_xxxxyz, g_y_x_x_xxxxz, g_y_x_x_xxxxzz, g_y_x_x_xxxyy, g_y_x_x_xxxyyy, g_y_x_x_xxxyyz, g_y_x_x_xxxyz, g_y_x_x_xxxyzz, g_y_x_x_xxxzz, g_y_x_x_xxxzzz, g_y_x_x_xxyyy, g_y_x_x_xxyyyy, g_y_x_x_xxyyyz, g_y_x_x_xxyyz, g_y_x_x_xxyyzz, g_y_x_x_xxyzz, g_y_x_x_xxyzzz, g_y_x_x_xxzzz, g_y_x_x_xxzzzz, g_y_x_x_xyyyy, g_y_x_x_xyyyyy, g_y_x_x_xyyyyz, g_y_x_x_xyyyz, g_y_x_x_xyyyzz, g_y_x_x_xyyzz, g_y_x_x_xyyzzz, g_y_x_x_xyzzz, g_y_x_x_xyzzzz, g_y_x_x_xzzzz, g_y_x_x_xzzzzz, g_y_x_x_yyyyy, g_y_x_x_yyyyz, g_y_x_x_yyyzz, g_y_x_x_yyzzz, g_y_x_x_yzzzz, g_y_x_x_zzzzz, g_y_x_xx_xxxxx, g_y_x_xx_xxxxy, g_y_x_xx_xxxxz, g_y_x_xx_xxxyy, g_y_x_xx_xxxyz, g_y_x_xx_xxxzz, g_y_x_xx_xxyyy, g_y_x_xx_xxyyz, g_y_x_xx_xxyzz, g_y_x_xx_xxzzz, g_y_x_xx_xyyyy, g_y_x_xx_xyyyz, g_y_x_xx_xyyzz, g_y_x_xx_xyzzz, g_y_x_xx_xzzzz, g_y_x_xx_yyyyy, g_y_x_xx_yyyyz, g_y_x_xx_yyyzz, g_y_x_xx_yyzzz, g_y_x_xx_yzzzz, g_y_x_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xx_xxxxx[k] = g_y_0_x_xxxxx[k] - g_y_x_x_xxxxx[k] * ab_x + g_y_x_x_xxxxxx[k];

                g_y_x_xx_xxxxy[k] = g_y_0_x_xxxxy[k] - g_y_x_x_xxxxy[k] * ab_x + g_y_x_x_xxxxxy[k];

                g_y_x_xx_xxxxz[k] = g_y_0_x_xxxxz[k] - g_y_x_x_xxxxz[k] * ab_x + g_y_x_x_xxxxxz[k];

                g_y_x_xx_xxxyy[k] = g_y_0_x_xxxyy[k] - g_y_x_x_xxxyy[k] * ab_x + g_y_x_x_xxxxyy[k];

                g_y_x_xx_xxxyz[k] = g_y_0_x_xxxyz[k] - g_y_x_x_xxxyz[k] * ab_x + g_y_x_x_xxxxyz[k];

                g_y_x_xx_xxxzz[k] = g_y_0_x_xxxzz[k] - g_y_x_x_xxxzz[k] * ab_x + g_y_x_x_xxxxzz[k];

                g_y_x_xx_xxyyy[k] = g_y_0_x_xxyyy[k] - g_y_x_x_xxyyy[k] * ab_x + g_y_x_x_xxxyyy[k];

                g_y_x_xx_xxyyz[k] = g_y_0_x_xxyyz[k] - g_y_x_x_xxyyz[k] * ab_x + g_y_x_x_xxxyyz[k];

                g_y_x_xx_xxyzz[k] = g_y_0_x_xxyzz[k] - g_y_x_x_xxyzz[k] * ab_x + g_y_x_x_xxxyzz[k];

                g_y_x_xx_xxzzz[k] = g_y_0_x_xxzzz[k] - g_y_x_x_xxzzz[k] * ab_x + g_y_x_x_xxxzzz[k];

                g_y_x_xx_xyyyy[k] = g_y_0_x_xyyyy[k] - g_y_x_x_xyyyy[k] * ab_x + g_y_x_x_xxyyyy[k];

                g_y_x_xx_xyyyz[k] = g_y_0_x_xyyyz[k] - g_y_x_x_xyyyz[k] * ab_x + g_y_x_x_xxyyyz[k];

                g_y_x_xx_xyyzz[k] = g_y_0_x_xyyzz[k] - g_y_x_x_xyyzz[k] * ab_x + g_y_x_x_xxyyzz[k];

                g_y_x_xx_xyzzz[k] = g_y_0_x_xyzzz[k] - g_y_x_x_xyzzz[k] * ab_x + g_y_x_x_xxyzzz[k];

                g_y_x_xx_xzzzz[k] = g_y_0_x_xzzzz[k] - g_y_x_x_xzzzz[k] * ab_x + g_y_x_x_xxzzzz[k];

                g_y_x_xx_yyyyy[k] = g_y_0_x_yyyyy[k] - g_y_x_x_yyyyy[k] * ab_x + g_y_x_x_xyyyyy[k];

                g_y_x_xx_yyyyz[k] = g_y_0_x_yyyyz[k] - g_y_x_x_yyyyz[k] * ab_x + g_y_x_x_xyyyyz[k];

                g_y_x_xx_yyyzz[k] = g_y_0_x_yyyzz[k] - g_y_x_x_yyyzz[k] * ab_x + g_y_x_x_xyyyzz[k];

                g_y_x_xx_yyzzz[k] = g_y_0_x_yyzzz[k] - g_y_x_x_yyzzz[k] * ab_x + g_y_x_x_xyyzzz[k];

                g_y_x_xx_yzzzz[k] = g_y_0_x_yzzzz[k] - g_y_x_x_yzzzz[k] * ab_x + g_y_x_x_xyzzzz[k];

                g_y_x_xx_zzzzz[k] = g_y_0_x_zzzzz[k] - g_y_x_x_zzzzz[k] * ab_x + g_y_x_x_xzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_y_x_xy_xxxxx = cbuffer.data(dh_geom_11_off + 399 * ccomps * dcomps);

            auto g_y_x_xy_xxxxy = cbuffer.data(dh_geom_11_off + 400 * ccomps * dcomps);

            auto g_y_x_xy_xxxxz = cbuffer.data(dh_geom_11_off + 401 * ccomps * dcomps);

            auto g_y_x_xy_xxxyy = cbuffer.data(dh_geom_11_off + 402 * ccomps * dcomps);

            auto g_y_x_xy_xxxyz = cbuffer.data(dh_geom_11_off + 403 * ccomps * dcomps);

            auto g_y_x_xy_xxxzz = cbuffer.data(dh_geom_11_off + 404 * ccomps * dcomps);

            auto g_y_x_xy_xxyyy = cbuffer.data(dh_geom_11_off + 405 * ccomps * dcomps);

            auto g_y_x_xy_xxyyz = cbuffer.data(dh_geom_11_off + 406 * ccomps * dcomps);

            auto g_y_x_xy_xxyzz = cbuffer.data(dh_geom_11_off + 407 * ccomps * dcomps);

            auto g_y_x_xy_xxzzz = cbuffer.data(dh_geom_11_off + 408 * ccomps * dcomps);

            auto g_y_x_xy_xyyyy = cbuffer.data(dh_geom_11_off + 409 * ccomps * dcomps);

            auto g_y_x_xy_xyyyz = cbuffer.data(dh_geom_11_off + 410 * ccomps * dcomps);

            auto g_y_x_xy_xyyzz = cbuffer.data(dh_geom_11_off + 411 * ccomps * dcomps);

            auto g_y_x_xy_xyzzz = cbuffer.data(dh_geom_11_off + 412 * ccomps * dcomps);

            auto g_y_x_xy_xzzzz = cbuffer.data(dh_geom_11_off + 413 * ccomps * dcomps);

            auto g_y_x_xy_yyyyy = cbuffer.data(dh_geom_11_off + 414 * ccomps * dcomps);

            auto g_y_x_xy_yyyyz = cbuffer.data(dh_geom_11_off + 415 * ccomps * dcomps);

            auto g_y_x_xy_yyyzz = cbuffer.data(dh_geom_11_off + 416 * ccomps * dcomps);

            auto g_y_x_xy_yyzzz = cbuffer.data(dh_geom_11_off + 417 * ccomps * dcomps);

            auto g_y_x_xy_yzzzz = cbuffer.data(dh_geom_11_off + 418 * ccomps * dcomps);

            auto g_y_x_xy_zzzzz = cbuffer.data(dh_geom_11_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_y_xxxxx, g_y_0_y_xxxxy, g_y_0_y_xxxxz, g_y_0_y_xxxyy, g_y_0_y_xxxyz, g_y_0_y_xxxzz, g_y_0_y_xxyyy, g_y_0_y_xxyyz, g_y_0_y_xxyzz, g_y_0_y_xxzzz, g_y_0_y_xyyyy, g_y_0_y_xyyyz, g_y_0_y_xyyzz, g_y_0_y_xyzzz, g_y_0_y_xzzzz, g_y_0_y_yyyyy, g_y_0_y_yyyyz, g_y_0_y_yyyzz, g_y_0_y_yyzzz, g_y_0_y_yzzzz, g_y_0_y_zzzzz, g_y_x_xy_xxxxx, g_y_x_xy_xxxxy, g_y_x_xy_xxxxz, g_y_x_xy_xxxyy, g_y_x_xy_xxxyz, g_y_x_xy_xxxzz, g_y_x_xy_xxyyy, g_y_x_xy_xxyyz, g_y_x_xy_xxyzz, g_y_x_xy_xxzzz, g_y_x_xy_xyyyy, g_y_x_xy_xyyyz, g_y_x_xy_xyyzz, g_y_x_xy_xyzzz, g_y_x_xy_xzzzz, g_y_x_xy_yyyyy, g_y_x_xy_yyyyz, g_y_x_xy_yyyzz, g_y_x_xy_yyzzz, g_y_x_xy_yzzzz, g_y_x_xy_zzzzz, g_y_x_y_xxxxx, g_y_x_y_xxxxxx, g_y_x_y_xxxxxy, g_y_x_y_xxxxxz, g_y_x_y_xxxxy, g_y_x_y_xxxxyy, g_y_x_y_xxxxyz, g_y_x_y_xxxxz, g_y_x_y_xxxxzz, g_y_x_y_xxxyy, g_y_x_y_xxxyyy, g_y_x_y_xxxyyz, g_y_x_y_xxxyz, g_y_x_y_xxxyzz, g_y_x_y_xxxzz, g_y_x_y_xxxzzz, g_y_x_y_xxyyy, g_y_x_y_xxyyyy, g_y_x_y_xxyyyz, g_y_x_y_xxyyz, g_y_x_y_xxyyzz, g_y_x_y_xxyzz, g_y_x_y_xxyzzz, g_y_x_y_xxzzz, g_y_x_y_xxzzzz, g_y_x_y_xyyyy, g_y_x_y_xyyyyy, g_y_x_y_xyyyyz, g_y_x_y_xyyyz, g_y_x_y_xyyyzz, g_y_x_y_xyyzz, g_y_x_y_xyyzzz, g_y_x_y_xyzzz, g_y_x_y_xyzzzz, g_y_x_y_xzzzz, g_y_x_y_xzzzzz, g_y_x_y_yyyyy, g_y_x_y_yyyyz, g_y_x_y_yyyzz, g_y_x_y_yyzzz, g_y_x_y_yzzzz, g_y_x_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xy_xxxxx[k] = g_y_0_y_xxxxx[k] - g_y_x_y_xxxxx[k] * ab_x + g_y_x_y_xxxxxx[k];

                g_y_x_xy_xxxxy[k] = g_y_0_y_xxxxy[k] - g_y_x_y_xxxxy[k] * ab_x + g_y_x_y_xxxxxy[k];

                g_y_x_xy_xxxxz[k] = g_y_0_y_xxxxz[k] - g_y_x_y_xxxxz[k] * ab_x + g_y_x_y_xxxxxz[k];

                g_y_x_xy_xxxyy[k] = g_y_0_y_xxxyy[k] - g_y_x_y_xxxyy[k] * ab_x + g_y_x_y_xxxxyy[k];

                g_y_x_xy_xxxyz[k] = g_y_0_y_xxxyz[k] - g_y_x_y_xxxyz[k] * ab_x + g_y_x_y_xxxxyz[k];

                g_y_x_xy_xxxzz[k] = g_y_0_y_xxxzz[k] - g_y_x_y_xxxzz[k] * ab_x + g_y_x_y_xxxxzz[k];

                g_y_x_xy_xxyyy[k] = g_y_0_y_xxyyy[k] - g_y_x_y_xxyyy[k] * ab_x + g_y_x_y_xxxyyy[k];

                g_y_x_xy_xxyyz[k] = g_y_0_y_xxyyz[k] - g_y_x_y_xxyyz[k] * ab_x + g_y_x_y_xxxyyz[k];

                g_y_x_xy_xxyzz[k] = g_y_0_y_xxyzz[k] - g_y_x_y_xxyzz[k] * ab_x + g_y_x_y_xxxyzz[k];

                g_y_x_xy_xxzzz[k] = g_y_0_y_xxzzz[k] - g_y_x_y_xxzzz[k] * ab_x + g_y_x_y_xxxzzz[k];

                g_y_x_xy_xyyyy[k] = g_y_0_y_xyyyy[k] - g_y_x_y_xyyyy[k] * ab_x + g_y_x_y_xxyyyy[k];

                g_y_x_xy_xyyyz[k] = g_y_0_y_xyyyz[k] - g_y_x_y_xyyyz[k] * ab_x + g_y_x_y_xxyyyz[k];

                g_y_x_xy_xyyzz[k] = g_y_0_y_xyyzz[k] - g_y_x_y_xyyzz[k] * ab_x + g_y_x_y_xxyyzz[k];

                g_y_x_xy_xyzzz[k] = g_y_0_y_xyzzz[k] - g_y_x_y_xyzzz[k] * ab_x + g_y_x_y_xxyzzz[k];

                g_y_x_xy_xzzzz[k] = g_y_0_y_xzzzz[k] - g_y_x_y_xzzzz[k] * ab_x + g_y_x_y_xxzzzz[k];

                g_y_x_xy_yyyyy[k] = g_y_0_y_yyyyy[k] - g_y_x_y_yyyyy[k] * ab_x + g_y_x_y_xyyyyy[k];

                g_y_x_xy_yyyyz[k] = g_y_0_y_yyyyz[k] - g_y_x_y_yyyyz[k] * ab_x + g_y_x_y_xyyyyz[k];

                g_y_x_xy_yyyzz[k] = g_y_0_y_yyyzz[k] - g_y_x_y_yyyzz[k] * ab_x + g_y_x_y_xyyyzz[k];

                g_y_x_xy_yyzzz[k] = g_y_0_y_yyzzz[k] - g_y_x_y_yyzzz[k] * ab_x + g_y_x_y_xyyzzz[k];

                g_y_x_xy_yzzzz[k] = g_y_0_y_yzzzz[k] - g_y_x_y_yzzzz[k] * ab_x + g_y_x_y_xyzzzz[k];

                g_y_x_xy_zzzzz[k] = g_y_0_y_zzzzz[k] - g_y_x_y_zzzzz[k] * ab_x + g_y_x_y_xzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_y_x_xz_xxxxx = cbuffer.data(dh_geom_11_off + 420 * ccomps * dcomps);

            auto g_y_x_xz_xxxxy = cbuffer.data(dh_geom_11_off + 421 * ccomps * dcomps);

            auto g_y_x_xz_xxxxz = cbuffer.data(dh_geom_11_off + 422 * ccomps * dcomps);

            auto g_y_x_xz_xxxyy = cbuffer.data(dh_geom_11_off + 423 * ccomps * dcomps);

            auto g_y_x_xz_xxxyz = cbuffer.data(dh_geom_11_off + 424 * ccomps * dcomps);

            auto g_y_x_xz_xxxzz = cbuffer.data(dh_geom_11_off + 425 * ccomps * dcomps);

            auto g_y_x_xz_xxyyy = cbuffer.data(dh_geom_11_off + 426 * ccomps * dcomps);

            auto g_y_x_xz_xxyyz = cbuffer.data(dh_geom_11_off + 427 * ccomps * dcomps);

            auto g_y_x_xz_xxyzz = cbuffer.data(dh_geom_11_off + 428 * ccomps * dcomps);

            auto g_y_x_xz_xxzzz = cbuffer.data(dh_geom_11_off + 429 * ccomps * dcomps);

            auto g_y_x_xz_xyyyy = cbuffer.data(dh_geom_11_off + 430 * ccomps * dcomps);

            auto g_y_x_xz_xyyyz = cbuffer.data(dh_geom_11_off + 431 * ccomps * dcomps);

            auto g_y_x_xz_xyyzz = cbuffer.data(dh_geom_11_off + 432 * ccomps * dcomps);

            auto g_y_x_xz_xyzzz = cbuffer.data(dh_geom_11_off + 433 * ccomps * dcomps);

            auto g_y_x_xz_xzzzz = cbuffer.data(dh_geom_11_off + 434 * ccomps * dcomps);

            auto g_y_x_xz_yyyyy = cbuffer.data(dh_geom_11_off + 435 * ccomps * dcomps);

            auto g_y_x_xz_yyyyz = cbuffer.data(dh_geom_11_off + 436 * ccomps * dcomps);

            auto g_y_x_xz_yyyzz = cbuffer.data(dh_geom_11_off + 437 * ccomps * dcomps);

            auto g_y_x_xz_yyzzz = cbuffer.data(dh_geom_11_off + 438 * ccomps * dcomps);

            auto g_y_x_xz_yzzzz = cbuffer.data(dh_geom_11_off + 439 * ccomps * dcomps);

            auto g_y_x_xz_zzzzz = cbuffer.data(dh_geom_11_off + 440 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_x_xxxxx, g_y_x_x_xxxxxz, g_y_x_x_xxxxy, g_y_x_x_xxxxyz, g_y_x_x_xxxxz, g_y_x_x_xxxxzz, g_y_x_x_xxxyy, g_y_x_x_xxxyyz, g_y_x_x_xxxyz, g_y_x_x_xxxyzz, g_y_x_x_xxxzz, g_y_x_x_xxxzzz, g_y_x_x_xxyyy, g_y_x_x_xxyyyz, g_y_x_x_xxyyz, g_y_x_x_xxyyzz, g_y_x_x_xxyzz, g_y_x_x_xxyzzz, g_y_x_x_xxzzz, g_y_x_x_xxzzzz, g_y_x_x_xyyyy, g_y_x_x_xyyyyz, g_y_x_x_xyyyz, g_y_x_x_xyyyzz, g_y_x_x_xyyzz, g_y_x_x_xyyzzz, g_y_x_x_xyzzz, g_y_x_x_xyzzzz, g_y_x_x_xzzzz, g_y_x_x_xzzzzz, g_y_x_x_yyyyy, g_y_x_x_yyyyyz, g_y_x_x_yyyyz, g_y_x_x_yyyyzz, g_y_x_x_yyyzz, g_y_x_x_yyyzzz, g_y_x_x_yyzzz, g_y_x_x_yyzzzz, g_y_x_x_yzzzz, g_y_x_x_yzzzzz, g_y_x_x_zzzzz, g_y_x_x_zzzzzz, g_y_x_xz_xxxxx, g_y_x_xz_xxxxy, g_y_x_xz_xxxxz, g_y_x_xz_xxxyy, g_y_x_xz_xxxyz, g_y_x_xz_xxxzz, g_y_x_xz_xxyyy, g_y_x_xz_xxyyz, g_y_x_xz_xxyzz, g_y_x_xz_xxzzz, g_y_x_xz_xyyyy, g_y_x_xz_xyyyz, g_y_x_xz_xyyzz, g_y_x_xz_xyzzz, g_y_x_xz_xzzzz, g_y_x_xz_yyyyy, g_y_x_xz_yyyyz, g_y_x_xz_yyyzz, g_y_x_xz_yyzzz, g_y_x_xz_yzzzz, g_y_x_xz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_xz_xxxxx[k] = -g_y_x_x_xxxxx[k] * ab_z + g_y_x_x_xxxxxz[k];

                g_y_x_xz_xxxxy[k] = -g_y_x_x_xxxxy[k] * ab_z + g_y_x_x_xxxxyz[k];

                g_y_x_xz_xxxxz[k] = -g_y_x_x_xxxxz[k] * ab_z + g_y_x_x_xxxxzz[k];

                g_y_x_xz_xxxyy[k] = -g_y_x_x_xxxyy[k] * ab_z + g_y_x_x_xxxyyz[k];

                g_y_x_xz_xxxyz[k] = -g_y_x_x_xxxyz[k] * ab_z + g_y_x_x_xxxyzz[k];

                g_y_x_xz_xxxzz[k] = -g_y_x_x_xxxzz[k] * ab_z + g_y_x_x_xxxzzz[k];

                g_y_x_xz_xxyyy[k] = -g_y_x_x_xxyyy[k] * ab_z + g_y_x_x_xxyyyz[k];

                g_y_x_xz_xxyyz[k] = -g_y_x_x_xxyyz[k] * ab_z + g_y_x_x_xxyyzz[k];

                g_y_x_xz_xxyzz[k] = -g_y_x_x_xxyzz[k] * ab_z + g_y_x_x_xxyzzz[k];

                g_y_x_xz_xxzzz[k] = -g_y_x_x_xxzzz[k] * ab_z + g_y_x_x_xxzzzz[k];

                g_y_x_xz_xyyyy[k] = -g_y_x_x_xyyyy[k] * ab_z + g_y_x_x_xyyyyz[k];

                g_y_x_xz_xyyyz[k] = -g_y_x_x_xyyyz[k] * ab_z + g_y_x_x_xyyyzz[k];

                g_y_x_xz_xyyzz[k] = -g_y_x_x_xyyzz[k] * ab_z + g_y_x_x_xyyzzz[k];

                g_y_x_xz_xyzzz[k] = -g_y_x_x_xyzzz[k] * ab_z + g_y_x_x_xyzzzz[k];

                g_y_x_xz_xzzzz[k] = -g_y_x_x_xzzzz[k] * ab_z + g_y_x_x_xzzzzz[k];

                g_y_x_xz_yyyyy[k] = -g_y_x_x_yyyyy[k] * ab_z + g_y_x_x_yyyyyz[k];

                g_y_x_xz_yyyyz[k] = -g_y_x_x_yyyyz[k] * ab_z + g_y_x_x_yyyyzz[k];

                g_y_x_xz_yyyzz[k] = -g_y_x_x_yyyzz[k] * ab_z + g_y_x_x_yyyzzz[k];

                g_y_x_xz_yyzzz[k] = -g_y_x_x_yyzzz[k] * ab_z + g_y_x_x_yyzzzz[k];

                g_y_x_xz_yzzzz[k] = -g_y_x_x_yzzzz[k] * ab_z + g_y_x_x_yzzzzz[k];

                g_y_x_xz_zzzzz[k] = -g_y_x_x_zzzzz[k] * ab_z + g_y_x_x_zzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

            auto g_y_x_yy_xxxxx = cbuffer.data(dh_geom_11_off + 441 * ccomps * dcomps);

            auto g_y_x_yy_xxxxy = cbuffer.data(dh_geom_11_off + 442 * ccomps * dcomps);

            auto g_y_x_yy_xxxxz = cbuffer.data(dh_geom_11_off + 443 * ccomps * dcomps);

            auto g_y_x_yy_xxxyy = cbuffer.data(dh_geom_11_off + 444 * ccomps * dcomps);

            auto g_y_x_yy_xxxyz = cbuffer.data(dh_geom_11_off + 445 * ccomps * dcomps);

            auto g_y_x_yy_xxxzz = cbuffer.data(dh_geom_11_off + 446 * ccomps * dcomps);

            auto g_y_x_yy_xxyyy = cbuffer.data(dh_geom_11_off + 447 * ccomps * dcomps);

            auto g_y_x_yy_xxyyz = cbuffer.data(dh_geom_11_off + 448 * ccomps * dcomps);

            auto g_y_x_yy_xxyzz = cbuffer.data(dh_geom_11_off + 449 * ccomps * dcomps);

            auto g_y_x_yy_xxzzz = cbuffer.data(dh_geom_11_off + 450 * ccomps * dcomps);

            auto g_y_x_yy_xyyyy = cbuffer.data(dh_geom_11_off + 451 * ccomps * dcomps);

            auto g_y_x_yy_xyyyz = cbuffer.data(dh_geom_11_off + 452 * ccomps * dcomps);

            auto g_y_x_yy_xyyzz = cbuffer.data(dh_geom_11_off + 453 * ccomps * dcomps);

            auto g_y_x_yy_xyzzz = cbuffer.data(dh_geom_11_off + 454 * ccomps * dcomps);

            auto g_y_x_yy_xzzzz = cbuffer.data(dh_geom_11_off + 455 * ccomps * dcomps);

            auto g_y_x_yy_yyyyy = cbuffer.data(dh_geom_11_off + 456 * ccomps * dcomps);

            auto g_y_x_yy_yyyyz = cbuffer.data(dh_geom_11_off + 457 * ccomps * dcomps);

            auto g_y_x_yy_yyyzz = cbuffer.data(dh_geom_11_off + 458 * ccomps * dcomps);

            auto g_y_x_yy_yyzzz = cbuffer.data(dh_geom_11_off + 459 * ccomps * dcomps);

            auto g_y_x_yy_yzzzz = cbuffer.data(dh_geom_11_off + 460 * ccomps * dcomps);

            auto g_y_x_yy_zzzzz = cbuffer.data(dh_geom_11_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_y_xxxxx, g_0_x_y_xxxxy, g_0_x_y_xxxxz, g_0_x_y_xxxyy, g_0_x_y_xxxyz, g_0_x_y_xxxzz, g_0_x_y_xxyyy, g_0_x_y_xxyyz, g_0_x_y_xxyzz, g_0_x_y_xxzzz, g_0_x_y_xyyyy, g_0_x_y_xyyyz, g_0_x_y_xyyzz, g_0_x_y_xyzzz, g_0_x_y_xzzzz, g_0_x_y_yyyyy, g_0_x_y_yyyyz, g_0_x_y_yyyzz, g_0_x_y_yyzzz, g_0_x_y_yzzzz, g_0_x_y_zzzzz, g_y_x_y_xxxxx, g_y_x_y_xxxxxy, g_y_x_y_xxxxy, g_y_x_y_xxxxyy, g_y_x_y_xxxxyz, g_y_x_y_xxxxz, g_y_x_y_xxxyy, g_y_x_y_xxxyyy, g_y_x_y_xxxyyz, g_y_x_y_xxxyz, g_y_x_y_xxxyzz, g_y_x_y_xxxzz, g_y_x_y_xxyyy, g_y_x_y_xxyyyy, g_y_x_y_xxyyyz, g_y_x_y_xxyyz, g_y_x_y_xxyyzz, g_y_x_y_xxyzz, g_y_x_y_xxyzzz, g_y_x_y_xxzzz, g_y_x_y_xyyyy, g_y_x_y_xyyyyy, g_y_x_y_xyyyyz, g_y_x_y_xyyyz, g_y_x_y_xyyyzz, g_y_x_y_xyyzz, g_y_x_y_xyyzzz, g_y_x_y_xyzzz, g_y_x_y_xyzzzz, g_y_x_y_xzzzz, g_y_x_y_yyyyy, g_y_x_y_yyyyyy, g_y_x_y_yyyyyz, g_y_x_y_yyyyz, g_y_x_y_yyyyzz, g_y_x_y_yyyzz, g_y_x_y_yyyzzz, g_y_x_y_yyzzz, g_y_x_y_yyzzzz, g_y_x_y_yzzzz, g_y_x_y_yzzzzz, g_y_x_y_zzzzz, g_y_x_yy_xxxxx, g_y_x_yy_xxxxy, g_y_x_yy_xxxxz, g_y_x_yy_xxxyy, g_y_x_yy_xxxyz, g_y_x_yy_xxxzz, g_y_x_yy_xxyyy, g_y_x_yy_xxyyz, g_y_x_yy_xxyzz, g_y_x_yy_xxzzz, g_y_x_yy_xyyyy, g_y_x_yy_xyyyz, g_y_x_yy_xyyzz, g_y_x_yy_xyzzz, g_y_x_yy_xzzzz, g_y_x_yy_yyyyy, g_y_x_yy_yyyyz, g_y_x_yy_yyyzz, g_y_x_yy_yyzzz, g_y_x_yy_yzzzz, g_y_x_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yy_xxxxx[k] = -g_0_x_y_xxxxx[k] - g_y_x_y_xxxxx[k] * ab_y + g_y_x_y_xxxxxy[k];

                g_y_x_yy_xxxxy[k] = -g_0_x_y_xxxxy[k] - g_y_x_y_xxxxy[k] * ab_y + g_y_x_y_xxxxyy[k];

                g_y_x_yy_xxxxz[k] = -g_0_x_y_xxxxz[k] - g_y_x_y_xxxxz[k] * ab_y + g_y_x_y_xxxxyz[k];

                g_y_x_yy_xxxyy[k] = -g_0_x_y_xxxyy[k] - g_y_x_y_xxxyy[k] * ab_y + g_y_x_y_xxxyyy[k];

                g_y_x_yy_xxxyz[k] = -g_0_x_y_xxxyz[k] - g_y_x_y_xxxyz[k] * ab_y + g_y_x_y_xxxyyz[k];

                g_y_x_yy_xxxzz[k] = -g_0_x_y_xxxzz[k] - g_y_x_y_xxxzz[k] * ab_y + g_y_x_y_xxxyzz[k];

                g_y_x_yy_xxyyy[k] = -g_0_x_y_xxyyy[k] - g_y_x_y_xxyyy[k] * ab_y + g_y_x_y_xxyyyy[k];

                g_y_x_yy_xxyyz[k] = -g_0_x_y_xxyyz[k] - g_y_x_y_xxyyz[k] * ab_y + g_y_x_y_xxyyyz[k];

                g_y_x_yy_xxyzz[k] = -g_0_x_y_xxyzz[k] - g_y_x_y_xxyzz[k] * ab_y + g_y_x_y_xxyyzz[k];

                g_y_x_yy_xxzzz[k] = -g_0_x_y_xxzzz[k] - g_y_x_y_xxzzz[k] * ab_y + g_y_x_y_xxyzzz[k];

                g_y_x_yy_xyyyy[k] = -g_0_x_y_xyyyy[k] - g_y_x_y_xyyyy[k] * ab_y + g_y_x_y_xyyyyy[k];

                g_y_x_yy_xyyyz[k] = -g_0_x_y_xyyyz[k] - g_y_x_y_xyyyz[k] * ab_y + g_y_x_y_xyyyyz[k];

                g_y_x_yy_xyyzz[k] = -g_0_x_y_xyyzz[k] - g_y_x_y_xyyzz[k] * ab_y + g_y_x_y_xyyyzz[k];

                g_y_x_yy_xyzzz[k] = -g_0_x_y_xyzzz[k] - g_y_x_y_xyzzz[k] * ab_y + g_y_x_y_xyyzzz[k];

                g_y_x_yy_xzzzz[k] = -g_0_x_y_xzzzz[k] - g_y_x_y_xzzzz[k] * ab_y + g_y_x_y_xyzzzz[k];

                g_y_x_yy_yyyyy[k] = -g_0_x_y_yyyyy[k] - g_y_x_y_yyyyy[k] * ab_y + g_y_x_y_yyyyyy[k];

                g_y_x_yy_yyyyz[k] = -g_0_x_y_yyyyz[k] - g_y_x_y_yyyyz[k] * ab_y + g_y_x_y_yyyyyz[k];

                g_y_x_yy_yyyzz[k] = -g_0_x_y_yyyzz[k] - g_y_x_y_yyyzz[k] * ab_y + g_y_x_y_yyyyzz[k];

                g_y_x_yy_yyzzz[k] = -g_0_x_y_yyzzz[k] - g_y_x_y_yyzzz[k] * ab_y + g_y_x_y_yyyzzz[k];

                g_y_x_yy_yzzzz[k] = -g_0_x_y_yzzzz[k] - g_y_x_y_yzzzz[k] * ab_y + g_y_x_y_yyzzzz[k];

                g_y_x_yy_zzzzz[k] = -g_0_x_y_zzzzz[k] - g_y_x_y_zzzzz[k] * ab_y + g_y_x_y_yzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

            auto g_y_x_yz_xxxxx = cbuffer.data(dh_geom_11_off + 462 * ccomps * dcomps);

            auto g_y_x_yz_xxxxy = cbuffer.data(dh_geom_11_off + 463 * ccomps * dcomps);

            auto g_y_x_yz_xxxxz = cbuffer.data(dh_geom_11_off + 464 * ccomps * dcomps);

            auto g_y_x_yz_xxxyy = cbuffer.data(dh_geom_11_off + 465 * ccomps * dcomps);

            auto g_y_x_yz_xxxyz = cbuffer.data(dh_geom_11_off + 466 * ccomps * dcomps);

            auto g_y_x_yz_xxxzz = cbuffer.data(dh_geom_11_off + 467 * ccomps * dcomps);

            auto g_y_x_yz_xxyyy = cbuffer.data(dh_geom_11_off + 468 * ccomps * dcomps);

            auto g_y_x_yz_xxyyz = cbuffer.data(dh_geom_11_off + 469 * ccomps * dcomps);

            auto g_y_x_yz_xxyzz = cbuffer.data(dh_geom_11_off + 470 * ccomps * dcomps);

            auto g_y_x_yz_xxzzz = cbuffer.data(dh_geom_11_off + 471 * ccomps * dcomps);

            auto g_y_x_yz_xyyyy = cbuffer.data(dh_geom_11_off + 472 * ccomps * dcomps);

            auto g_y_x_yz_xyyyz = cbuffer.data(dh_geom_11_off + 473 * ccomps * dcomps);

            auto g_y_x_yz_xyyzz = cbuffer.data(dh_geom_11_off + 474 * ccomps * dcomps);

            auto g_y_x_yz_xyzzz = cbuffer.data(dh_geom_11_off + 475 * ccomps * dcomps);

            auto g_y_x_yz_xzzzz = cbuffer.data(dh_geom_11_off + 476 * ccomps * dcomps);

            auto g_y_x_yz_yyyyy = cbuffer.data(dh_geom_11_off + 477 * ccomps * dcomps);

            auto g_y_x_yz_yyyyz = cbuffer.data(dh_geom_11_off + 478 * ccomps * dcomps);

            auto g_y_x_yz_yyyzz = cbuffer.data(dh_geom_11_off + 479 * ccomps * dcomps);

            auto g_y_x_yz_yyzzz = cbuffer.data(dh_geom_11_off + 480 * ccomps * dcomps);

            auto g_y_x_yz_yzzzz = cbuffer.data(dh_geom_11_off + 481 * ccomps * dcomps);

            auto g_y_x_yz_zzzzz = cbuffer.data(dh_geom_11_off + 482 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_y_xxxxx, g_y_x_y_xxxxxz, g_y_x_y_xxxxy, g_y_x_y_xxxxyz, g_y_x_y_xxxxz, g_y_x_y_xxxxzz, g_y_x_y_xxxyy, g_y_x_y_xxxyyz, g_y_x_y_xxxyz, g_y_x_y_xxxyzz, g_y_x_y_xxxzz, g_y_x_y_xxxzzz, g_y_x_y_xxyyy, g_y_x_y_xxyyyz, g_y_x_y_xxyyz, g_y_x_y_xxyyzz, g_y_x_y_xxyzz, g_y_x_y_xxyzzz, g_y_x_y_xxzzz, g_y_x_y_xxzzzz, g_y_x_y_xyyyy, g_y_x_y_xyyyyz, g_y_x_y_xyyyz, g_y_x_y_xyyyzz, g_y_x_y_xyyzz, g_y_x_y_xyyzzz, g_y_x_y_xyzzz, g_y_x_y_xyzzzz, g_y_x_y_xzzzz, g_y_x_y_xzzzzz, g_y_x_y_yyyyy, g_y_x_y_yyyyyz, g_y_x_y_yyyyz, g_y_x_y_yyyyzz, g_y_x_y_yyyzz, g_y_x_y_yyyzzz, g_y_x_y_yyzzz, g_y_x_y_yyzzzz, g_y_x_y_yzzzz, g_y_x_y_yzzzzz, g_y_x_y_zzzzz, g_y_x_y_zzzzzz, g_y_x_yz_xxxxx, g_y_x_yz_xxxxy, g_y_x_yz_xxxxz, g_y_x_yz_xxxyy, g_y_x_yz_xxxyz, g_y_x_yz_xxxzz, g_y_x_yz_xxyyy, g_y_x_yz_xxyyz, g_y_x_yz_xxyzz, g_y_x_yz_xxzzz, g_y_x_yz_xyyyy, g_y_x_yz_xyyyz, g_y_x_yz_xyyzz, g_y_x_yz_xyzzz, g_y_x_yz_xzzzz, g_y_x_yz_yyyyy, g_y_x_yz_yyyyz, g_y_x_yz_yyyzz, g_y_x_yz_yyzzz, g_y_x_yz_yzzzz, g_y_x_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_yz_xxxxx[k] = -g_y_x_y_xxxxx[k] * ab_z + g_y_x_y_xxxxxz[k];

                g_y_x_yz_xxxxy[k] = -g_y_x_y_xxxxy[k] * ab_z + g_y_x_y_xxxxyz[k];

                g_y_x_yz_xxxxz[k] = -g_y_x_y_xxxxz[k] * ab_z + g_y_x_y_xxxxzz[k];

                g_y_x_yz_xxxyy[k] = -g_y_x_y_xxxyy[k] * ab_z + g_y_x_y_xxxyyz[k];

                g_y_x_yz_xxxyz[k] = -g_y_x_y_xxxyz[k] * ab_z + g_y_x_y_xxxyzz[k];

                g_y_x_yz_xxxzz[k] = -g_y_x_y_xxxzz[k] * ab_z + g_y_x_y_xxxzzz[k];

                g_y_x_yz_xxyyy[k] = -g_y_x_y_xxyyy[k] * ab_z + g_y_x_y_xxyyyz[k];

                g_y_x_yz_xxyyz[k] = -g_y_x_y_xxyyz[k] * ab_z + g_y_x_y_xxyyzz[k];

                g_y_x_yz_xxyzz[k] = -g_y_x_y_xxyzz[k] * ab_z + g_y_x_y_xxyzzz[k];

                g_y_x_yz_xxzzz[k] = -g_y_x_y_xxzzz[k] * ab_z + g_y_x_y_xxzzzz[k];

                g_y_x_yz_xyyyy[k] = -g_y_x_y_xyyyy[k] * ab_z + g_y_x_y_xyyyyz[k];

                g_y_x_yz_xyyyz[k] = -g_y_x_y_xyyyz[k] * ab_z + g_y_x_y_xyyyzz[k];

                g_y_x_yz_xyyzz[k] = -g_y_x_y_xyyzz[k] * ab_z + g_y_x_y_xyyzzz[k];

                g_y_x_yz_xyzzz[k] = -g_y_x_y_xyzzz[k] * ab_z + g_y_x_y_xyzzzz[k];

                g_y_x_yz_xzzzz[k] = -g_y_x_y_xzzzz[k] * ab_z + g_y_x_y_xzzzzz[k];

                g_y_x_yz_yyyyy[k] = -g_y_x_y_yyyyy[k] * ab_z + g_y_x_y_yyyyyz[k];

                g_y_x_yz_yyyyz[k] = -g_y_x_y_yyyyz[k] * ab_z + g_y_x_y_yyyyzz[k];

                g_y_x_yz_yyyzz[k] = -g_y_x_y_yyyzz[k] * ab_z + g_y_x_y_yyyzzz[k];

                g_y_x_yz_yyzzz[k] = -g_y_x_y_yyzzz[k] * ab_z + g_y_x_y_yyzzzz[k];

                g_y_x_yz_yzzzz[k] = -g_y_x_y_yzzzz[k] * ab_z + g_y_x_y_yzzzzz[k];

                g_y_x_yz_zzzzz[k] = -g_y_x_y_zzzzz[k] * ab_z + g_y_x_y_zzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

            auto g_y_x_zz_xxxxx = cbuffer.data(dh_geom_11_off + 483 * ccomps * dcomps);

            auto g_y_x_zz_xxxxy = cbuffer.data(dh_geom_11_off + 484 * ccomps * dcomps);

            auto g_y_x_zz_xxxxz = cbuffer.data(dh_geom_11_off + 485 * ccomps * dcomps);

            auto g_y_x_zz_xxxyy = cbuffer.data(dh_geom_11_off + 486 * ccomps * dcomps);

            auto g_y_x_zz_xxxyz = cbuffer.data(dh_geom_11_off + 487 * ccomps * dcomps);

            auto g_y_x_zz_xxxzz = cbuffer.data(dh_geom_11_off + 488 * ccomps * dcomps);

            auto g_y_x_zz_xxyyy = cbuffer.data(dh_geom_11_off + 489 * ccomps * dcomps);

            auto g_y_x_zz_xxyyz = cbuffer.data(dh_geom_11_off + 490 * ccomps * dcomps);

            auto g_y_x_zz_xxyzz = cbuffer.data(dh_geom_11_off + 491 * ccomps * dcomps);

            auto g_y_x_zz_xxzzz = cbuffer.data(dh_geom_11_off + 492 * ccomps * dcomps);

            auto g_y_x_zz_xyyyy = cbuffer.data(dh_geom_11_off + 493 * ccomps * dcomps);

            auto g_y_x_zz_xyyyz = cbuffer.data(dh_geom_11_off + 494 * ccomps * dcomps);

            auto g_y_x_zz_xyyzz = cbuffer.data(dh_geom_11_off + 495 * ccomps * dcomps);

            auto g_y_x_zz_xyzzz = cbuffer.data(dh_geom_11_off + 496 * ccomps * dcomps);

            auto g_y_x_zz_xzzzz = cbuffer.data(dh_geom_11_off + 497 * ccomps * dcomps);

            auto g_y_x_zz_yyyyy = cbuffer.data(dh_geom_11_off + 498 * ccomps * dcomps);

            auto g_y_x_zz_yyyyz = cbuffer.data(dh_geom_11_off + 499 * ccomps * dcomps);

            auto g_y_x_zz_yyyzz = cbuffer.data(dh_geom_11_off + 500 * ccomps * dcomps);

            auto g_y_x_zz_yyzzz = cbuffer.data(dh_geom_11_off + 501 * ccomps * dcomps);

            auto g_y_x_zz_yzzzz = cbuffer.data(dh_geom_11_off + 502 * ccomps * dcomps);

            auto g_y_x_zz_zzzzz = cbuffer.data(dh_geom_11_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_x_z_xxxxx, g_y_x_z_xxxxxz, g_y_x_z_xxxxy, g_y_x_z_xxxxyz, g_y_x_z_xxxxz, g_y_x_z_xxxxzz, g_y_x_z_xxxyy, g_y_x_z_xxxyyz, g_y_x_z_xxxyz, g_y_x_z_xxxyzz, g_y_x_z_xxxzz, g_y_x_z_xxxzzz, g_y_x_z_xxyyy, g_y_x_z_xxyyyz, g_y_x_z_xxyyz, g_y_x_z_xxyyzz, g_y_x_z_xxyzz, g_y_x_z_xxyzzz, g_y_x_z_xxzzz, g_y_x_z_xxzzzz, g_y_x_z_xyyyy, g_y_x_z_xyyyyz, g_y_x_z_xyyyz, g_y_x_z_xyyyzz, g_y_x_z_xyyzz, g_y_x_z_xyyzzz, g_y_x_z_xyzzz, g_y_x_z_xyzzzz, g_y_x_z_xzzzz, g_y_x_z_xzzzzz, g_y_x_z_yyyyy, g_y_x_z_yyyyyz, g_y_x_z_yyyyz, g_y_x_z_yyyyzz, g_y_x_z_yyyzz, g_y_x_z_yyyzzz, g_y_x_z_yyzzz, g_y_x_z_yyzzzz, g_y_x_z_yzzzz, g_y_x_z_yzzzzz, g_y_x_z_zzzzz, g_y_x_z_zzzzzz, g_y_x_zz_xxxxx, g_y_x_zz_xxxxy, g_y_x_zz_xxxxz, g_y_x_zz_xxxyy, g_y_x_zz_xxxyz, g_y_x_zz_xxxzz, g_y_x_zz_xxyyy, g_y_x_zz_xxyyz, g_y_x_zz_xxyzz, g_y_x_zz_xxzzz, g_y_x_zz_xyyyy, g_y_x_zz_xyyyz, g_y_x_zz_xyyzz, g_y_x_zz_xyzzz, g_y_x_zz_xzzzz, g_y_x_zz_yyyyy, g_y_x_zz_yyyyz, g_y_x_zz_yyyzz, g_y_x_zz_yyzzz, g_y_x_zz_yzzzz, g_y_x_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_x_zz_xxxxx[k] = -g_y_x_z_xxxxx[k] * ab_z + g_y_x_z_xxxxxz[k];

                g_y_x_zz_xxxxy[k] = -g_y_x_z_xxxxy[k] * ab_z + g_y_x_z_xxxxyz[k];

                g_y_x_zz_xxxxz[k] = -g_y_x_z_xxxxz[k] * ab_z + g_y_x_z_xxxxzz[k];

                g_y_x_zz_xxxyy[k] = -g_y_x_z_xxxyy[k] * ab_z + g_y_x_z_xxxyyz[k];

                g_y_x_zz_xxxyz[k] = -g_y_x_z_xxxyz[k] * ab_z + g_y_x_z_xxxyzz[k];

                g_y_x_zz_xxxzz[k] = -g_y_x_z_xxxzz[k] * ab_z + g_y_x_z_xxxzzz[k];

                g_y_x_zz_xxyyy[k] = -g_y_x_z_xxyyy[k] * ab_z + g_y_x_z_xxyyyz[k];

                g_y_x_zz_xxyyz[k] = -g_y_x_z_xxyyz[k] * ab_z + g_y_x_z_xxyyzz[k];

                g_y_x_zz_xxyzz[k] = -g_y_x_z_xxyzz[k] * ab_z + g_y_x_z_xxyzzz[k];

                g_y_x_zz_xxzzz[k] = -g_y_x_z_xxzzz[k] * ab_z + g_y_x_z_xxzzzz[k];

                g_y_x_zz_xyyyy[k] = -g_y_x_z_xyyyy[k] * ab_z + g_y_x_z_xyyyyz[k];

                g_y_x_zz_xyyyz[k] = -g_y_x_z_xyyyz[k] * ab_z + g_y_x_z_xyyyzz[k];

                g_y_x_zz_xyyzz[k] = -g_y_x_z_xyyzz[k] * ab_z + g_y_x_z_xyyzzz[k];

                g_y_x_zz_xyzzz[k] = -g_y_x_z_xyzzz[k] * ab_z + g_y_x_z_xyzzzz[k];

                g_y_x_zz_xzzzz[k] = -g_y_x_z_xzzzz[k] * ab_z + g_y_x_z_xzzzzz[k];

                g_y_x_zz_yyyyy[k] = -g_y_x_z_yyyyy[k] * ab_z + g_y_x_z_yyyyyz[k];

                g_y_x_zz_yyyyz[k] = -g_y_x_z_yyyyz[k] * ab_z + g_y_x_z_yyyyzz[k];

                g_y_x_zz_yyyzz[k] = -g_y_x_z_yyyzz[k] * ab_z + g_y_x_z_yyyzzz[k];

                g_y_x_zz_yyzzz[k] = -g_y_x_z_yyzzz[k] * ab_z + g_y_x_z_yyzzzz[k];

                g_y_x_zz_yzzzz[k] = -g_y_x_z_yzzzz[k] * ab_z + g_y_x_z_yzzzzz[k];

                g_y_x_zz_zzzzz[k] = -g_y_x_z_zzzzz[k] * ab_z + g_y_x_z_zzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

            auto g_y_y_xx_xxxxx = cbuffer.data(dh_geom_11_off + 504 * ccomps * dcomps);

            auto g_y_y_xx_xxxxy = cbuffer.data(dh_geom_11_off + 505 * ccomps * dcomps);

            auto g_y_y_xx_xxxxz = cbuffer.data(dh_geom_11_off + 506 * ccomps * dcomps);

            auto g_y_y_xx_xxxyy = cbuffer.data(dh_geom_11_off + 507 * ccomps * dcomps);

            auto g_y_y_xx_xxxyz = cbuffer.data(dh_geom_11_off + 508 * ccomps * dcomps);

            auto g_y_y_xx_xxxzz = cbuffer.data(dh_geom_11_off + 509 * ccomps * dcomps);

            auto g_y_y_xx_xxyyy = cbuffer.data(dh_geom_11_off + 510 * ccomps * dcomps);

            auto g_y_y_xx_xxyyz = cbuffer.data(dh_geom_11_off + 511 * ccomps * dcomps);

            auto g_y_y_xx_xxyzz = cbuffer.data(dh_geom_11_off + 512 * ccomps * dcomps);

            auto g_y_y_xx_xxzzz = cbuffer.data(dh_geom_11_off + 513 * ccomps * dcomps);

            auto g_y_y_xx_xyyyy = cbuffer.data(dh_geom_11_off + 514 * ccomps * dcomps);

            auto g_y_y_xx_xyyyz = cbuffer.data(dh_geom_11_off + 515 * ccomps * dcomps);

            auto g_y_y_xx_xyyzz = cbuffer.data(dh_geom_11_off + 516 * ccomps * dcomps);

            auto g_y_y_xx_xyzzz = cbuffer.data(dh_geom_11_off + 517 * ccomps * dcomps);

            auto g_y_y_xx_xzzzz = cbuffer.data(dh_geom_11_off + 518 * ccomps * dcomps);

            auto g_y_y_xx_yyyyy = cbuffer.data(dh_geom_11_off + 519 * ccomps * dcomps);

            auto g_y_y_xx_yyyyz = cbuffer.data(dh_geom_11_off + 520 * ccomps * dcomps);

            auto g_y_y_xx_yyyzz = cbuffer.data(dh_geom_11_off + 521 * ccomps * dcomps);

            auto g_y_y_xx_yyzzz = cbuffer.data(dh_geom_11_off + 522 * ccomps * dcomps);

            auto g_y_y_xx_yzzzz = cbuffer.data(dh_geom_11_off + 523 * ccomps * dcomps);

            auto g_y_y_xx_zzzzz = cbuffer.data(dh_geom_11_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_x_xxxxx, g_y_y_x_xxxxxx, g_y_y_x_xxxxxy, g_y_y_x_xxxxxz, g_y_y_x_xxxxy, g_y_y_x_xxxxyy, g_y_y_x_xxxxyz, g_y_y_x_xxxxz, g_y_y_x_xxxxzz, g_y_y_x_xxxyy, g_y_y_x_xxxyyy, g_y_y_x_xxxyyz, g_y_y_x_xxxyz, g_y_y_x_xxxyzz, g_y_y_x_xxxzz, g_y_y_x_xxxzzz, g_y_y_x_xxyyy, g_y_y_x_xxyyyy, g_y_y_x_xxyyyz, g_y_y_x_xxyyz, g_y_y_x_xxyyzz, g_y_y_x_xxyzz, g_y_y_x_xxyzzz, g_y_y_x_xxzzz, g_y_y_x_xxzzzz, g_y_y_x_xyyyy, g_y_y_x_xyyyyy, g_y_y_x_xyyyyz, g_y_y_x_xyyyz, g_y_y_x_xyyyzz, g_y_y_x_xyyzz, g_y_y_x_xyyzzz, g_y_y_x_xyzzz, g_y_y_x_xyzzzz, g_y_y_x_xzzzz, g_y_y_x_xzzzzz, g_y_y_x_yyyyy, g_y_y_x_yyyyz, g_y_y_x_yyyzz, g_y_y_x_yyzzz, g_y_y_x_yzzzz, g_y_y_x_zzzzz, g_y_y_xx_xxxxx, g_y_y_xx_xxxxy, g_y_y_xx_xxxxz, g_y_y_xx_xxxyy, g_y_y_xx_xxxyz, g_y_y_xx_xxxzz, g_y_y_xx_xxyyy, g_y_y_xx_xxyyz, g_y_y_xx_xxyzz, g_y_y_xx_xxzzz, g_y_y_xx_xyyyy, g_y_y_xx_xyyyz, g_y_y_xx_xyyzz, g_y_y_xx_xyzzz, g_y_y_xx_xzzzz, g_y_y_xx_yyyyy, g_y_y_xx_yyyyz, g_y_y_xx_yyyzz, g_y_y_xx_yyzzz, g_y_y_xx_yzzzz, g_y_y_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xx_xxxxx[k] = -g_y_y_x_xxxxx[k] * ab_x + g_y_y_x_xxxxxx[k];

                g_y_y_xx_xxxxy[k] = -g_y_y_x_xxxxy[k] * ab_x + g_y_y_x_xxxxxy[k];

                g_y_y_xx_xxxxz[k] = -g_y_y_x_xxxxz[k] * ab_x + g_y_y_x_xxxxxz[k];

                g_y_y_xx_xxxyy[k] = -g_y_y_x_xxxyy[k] * ab_x + g_y_y_x_xxxxyy[k];

                g_y_y_xx_xxxyz[k] = -g_y_y_x_xxxyz[k] * ab_x + g_y_y_x_xxxxyz[k];

                g_y_y_xx_xxxzz[k] = -g_y_y_x_xxxzz[k] * ab_x + g_y_y_x_xxxxzz[k];

                g_y_y_xx_xxyyy[k] = -g_y_y_x_xxyyy[k] * ab_x + g_y_y_x_xxxyyy[k];

                g_y_y_xx_xxyyz[k] = -g_y_y_x_xxyyz[k] * ab_x + g_y_y_x_xxxyyz[k];

                g_y_y_xx_xxyzz[k] = -g_y_y_x_xxyzz[k] * ab_x + g_y_y_x_xxxyzz[k];

                g_y_y_xx_xxzzz[k] = -g_y_y_x_xxzzz[k] * ab_x + g_y_y_x_xxxzzz[k];

                g_y_y_xx_xyyyy[k] = -g_y_y_x_xyyyy[k] * ab_x + g_y_y_x_xxyyyy[k];

                g_y_y_xx_xyyyz[k] = -g_y_y_x_xyyyz[k] * ab_x + g_y_y_x_xxyyyz[k];

                g_y_y_xx_xyyzz[k] = -g_y_y_x_xyyzz[k] * ab_x + g_y_y_x_xxyyzz[k];

                g_y_y_xx_xyzzz[k] = -g_y_y_x_xyzzz[k] * ab_x + g_y_y_x_xxyzzz[k];

                g_y_y_xx_xzzzz[k] = -g_y_y_x_xzzzz[k] * ab_x + g_y_y_x_xxzzzz[k];

                g_y_y_xx_yyyyy[k] = -g_y_y_x_yyyyy[k] * ab_x + g_y_y_x_xyyyyy[k];

                g_y_y_xx_yyyyz[k] = -g_y_y_x_yyyyz[k] * ab_x + g_y_y_x_xyyyyz[k];

                g_y_y_xx_yyyzz[k] = -g_y_y_x_yyyzz[k] * ab_x + g_y_y_x_xyyyzz[k];

                g_y_y_xx_yyzzz[k] = -g_y_y_x_yyzzz[k] * ab_x + g_y_y_x_xyyzzz[k];

                g_y_y_xx_yzzzz[k] = -g_y_y_x_yzzzz[k] * ab_x + g_y_y_x_xyzzzz[k];

                g_y_y_xx_zzzzz[k] = -g_y_y_x_zzzzz[k] * ab_x + g_y_y_x_xzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

            auto g_y_y_xy_xxxxx = cbuffer.data(dh_geom_11_off + 525 * ccomps * dcomps);

            auto g_y_y_xy_xxxxy = cbuffer.data(dh_geom_11_off + 526 * ccomps * dcomps);

            auto g_y_y_xy_xxxxz = cbuffer.data(dh_geom_11_off + 527 * ccomps * dcomps);

            auto g_y_y_xy_xxxyy = cbuffer.data(dh_geom_11_off + 528 * ccomps * dcomps);

            auto g_y_y_xy_xxxyz = cbuffer.data(dh_geom_11_off + 529 * ccomps * dcomps);

            auto g_y_y_xy_xxxzz = cbuffer.data(dh_geom_11_off + 530 * ccomps * dcomps);

            auto g_y_y_xy_xxyyy = cbuffer.data(dh_geom_11_off + 531 * ccomps * dcomps);

            auto g_y_y_xy_xxyyz = cbuffer.data(dh_geom_11_off + 532 * ccomps * dcomps);

            auto g_y_y_xy_xxyzz = cbuffer.data(dh_geom_11_off + 533 * ccomps * dcomps);

            auto g_y_y_xy_xxzzz = cbuffer.data(dh_geom_11_off + 534 * ccomps * dcomps);

            auto g_y_y_xy_xyyyy = cbuffer.data(dh_geom_11_off + 535 * ccomps * dcomps);

            auto g_y_y_xy_xyyyz = cbuffer.data(dh_geom_11_off + 536 * ccomps * dcomps);

            auto g_y_y_xy_xyyzz = cbuffer.data(dh_geom_11_off + 537 * ccomps * dcomps);

            auto g_y_y_xy_xyzzz = cbuffer.data(dh_geom_11_off + 538 * ccomps * dcomps);

            auto g_y_y_xy_xzzzz = cbuffer.data(dh_geom_11_off + 539 * ccomps * dcomps);

            auto g_y_y_xy_yyyyy = cbuffer.data(dh_geom_11_off + 540 * ccomps * dcomps);

            auto g_y_y_xy_yyyyz = cbuffer.data(dh_geom_11_off + 541 * ccomps * dcomps);

            auto g_y_y_xy_yyyzz = cbuffer.data(dh_geom_11_off + 542 * ccomps * dcomps);

            auto g_y_y_xy_yyzzz = cbuffer.data(dh_geom_11_off + 543 * ccomps * dcomps);

            auto g_y_y_xy_yzzzz = cbuffer.data(dh_geom_11_off + 544 * ccomps * dcomps);

            auto g_y_y_xy_zzzzz = cbuffer.data(dh_geom_11_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xy_xxxxx, g_y_y_xy_xxxxy, g_y_y_xy_xxxxz, g_y_y_xy_xxxyy, g_y_y_xy_xxxyz, g_y_y_xy_xxxzz, g_y_y_xy_xxyyy, g_y_y_xy_xxyyz, g_y_y_xy_xxyzz, g_y_y_xy_xxzzz, g_y_y_xy_xyyyy, g_y_y_xy_xyyyz, g_y_y_xy_xyyzz, g_y_y_xy_xyzzz, g_y_y_xy_xzzzz, g_y_y_xy_yyyyy, g_y_y_xy_yyyyz, g_y_y_xy_yyyzz, g_y_y_xy_yyzzz, g_y_y_xy_yzzzz, g_y_y_xy_zzzzz, g_y_y_y_xxxxx, g_y_y_y_xxxxxx, g_y_y_y_xxxxxy, g_y_y_y_xxxxxz, g_y_y_y_xxxxy, g_y_y_y_xxxxyy, g_y_y_y_xxxxyz, g_y_y_y_xxxxz, g_y_y_y_xxxxzz, g_y_y_y_xxxyy, g_y_y_y_xxxyyy, g_y_y_y_xxxyyz, g_y_y_y_xxxyz, g_y_y_y_xxxyzz, g_y_y_y_xxxzz, g_y_y_y_xxxzzz, g_y_y_y_xxyyy, g_y_y_y_xxyyyy, g_y_y_y_xxyyyz, g_y_y_y_xxyyz, g_y_y_y_xxyyzz, g_y_y_y_xxyzz, g_y_y_y_xxyzzz, g_y_y_y_xxzzz, g_y_y_y_xxzzzz, g_y_y_y_xyyyy, g_y_y_y_xyyyyy, g_y_y_y_xyyyyz, g_y_y_y_xyyyz, g_y_y_y_xyyyzz, g_y_y_y_xyyzz, g_y_y_y_xyyzzz, g_y_y_y_xyzzz, g_y_y_y_xyzzzz, g_y_y_y_xzzzz, g_y_y_y_xzzzzz, g_y_y_y_yyyyy, g_y_y_y_yyyyz, g_y_y_y_yyyzz, g_y_y_y_yyzzz, g_y_y_y_yzzzz, g_y_y_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xy_xxxxx[k] = -g_y_y_y_xxxxx[k] * ab_x + g_y_y_y_xxxxxx[k];

                g_y_y_xy_xxxxy[k] = -g_y_y_y_xxxxy[k] * ab_x + g_y_y_y_xxxxxy[k];

                g_y_y_xy_xxxxz[k] = -g_y_y_y_xxxxz[k] * ab_x + g_y_y_y_xxxxxz[k];

                g_y_y_xy_xxxyy[k] = -g_y_y_y_xxxyy[k] * ab_x + g_y_y_y_xxxxyy[k];

                g_y_y_xy_xxxyz[k] = -g_y_y_y_xxxyz[k] * ab_x + g_y_y_y_xxxxyz[k];

                g_y_y_xy_xxxzz[k] = -g_y_y_y_xxxzz[k] * ab_x + g_y_y_y_xxxxzz[k];

                g_y_y_xy_xxyyy[k] = -g_y_y_y_xxyyy[k] * ab_x + g_y_y_y_xxxyyy[k];

                g_y_y_xy_xxyyz[k] = -g_y_y_y_xxyyz[k] * ab_x + g_y_y_y_xxxyyz[k];

                g_y_y_xy_xxyzz[k] = -g_y_y_y_xxyzz[k] * ab_x + g_y_y_y_xxxyzz[k];

                g_y_y_xy_xxzzz[k] = -g_y_y_y_xxzzz[k] * ab_x + g_y_y_y_xxxzzz[k];

                g_y_y_xy_xyyyy[k] = -g_y_y_y_xyyyy[k] * ab_x + g_y_y_y_xxyyyy[k];

                g_y_y_xy_xyyyz[k] = -g_y_y_y_xyyyz[k] * ab_x + g_y_y_y_xxyyyz[k];

                g_y_y_xy_xyyzz[k] = -g_y_y_y_xyyzz[k] * ab_x + g_y_y_y_xxyyzz[k];

                g_y_y_xy_xyzzz[k] = -g_y_y_y_xyzzz[k] * ab_x + g_y_y_y_xxyzzz[k];

                g_y_y_xy_xzzzz[k] = -g_y_y_y_xzzzz[k] * ab_x + g_y_y_y_xxzzzz[k];

                g_y_y_xy_yyyyy[k] = -g_y_y_y_yyyyy[k] * ab_x + g_y_y_y_xyyyyy[k];

                g_y_y_xy_yyyyz[k] = -g_y_y_y_yyyyz[k] * ab_x + g_y_y_y_xyyyyz[k];

                g_y_y_xy_yyyzz[k] = -g_y_y_y_yyyzz[k] * ab_x + g_y_y_y_xyyyzz[k];

                g_y_y_xy_yyzzz[k] = -g_y_y_y_yyzzz[k] * ab_x + g_y_y_y_xyyzzz[k];

                g_y_y_xy_yzzzz[k] = -g_y_y_y_yzzzz[k] * ab_x + g_y_y_y_xyzzzz[k];

                g_y_y_xy_zzzzz[k] = -g_y_y_y_zzzzz[k] * ab_x + g_y_y_y_xzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

            auto g_y_y_xz_xxxxx = cbuffer.data(dh_geom_11_off + 546 * ccomps * dcomps);

            auto g_y_y_xz_xxxxy = cbuffer.data(dh_geom_11_off + 547 * ccomps * dcomps);

            auto g_y_y_xz_xxxxz = cbuffer.data(dh_geom_11_off + 548 * ccomps * dcomps);

            auto g_y_y_xz_xxxyy = cbuffer.data(dh_geom_11_off + 549 * ccomps * dcomps);

            auto g_y_y_xz_xxxyz = cbuffer.data(dh_geom_11_off + 550 * ccomps * dcomps);

            auto g_y_y_xz_xxxzz = cbuffer.data(dh_geom_11_off + 551 * ccomps * dcomps);

            auto g_y_y_xz_xxyyy = cbuffer.data(dh_geom_11_off + 552 * ccomps * dcomps);

            auto g_y_y_xz_xxyyz = cbuffer.data(dh_geom_11_off + 553 * ccomps * dcomps);

            auto g_y_y_xz_xxyzz = cbuffer.data(dh_geom_11_off + 554 * ccomps * dcomps);

            auto g_y_y_xz_xxzzz = cbuffer.data(dh_geom_11_off + 555 * ccomps * dcomps);

            auto g_y_y_xz_xyyyy = cbuffer.data(dh_geom_11_off + 556 * ccomps * dcomps);

            auto g_y_y_xz_xyyyz = cbuffer.data(dh_geom_11_off + 557 * ccomps * dcomps);

            auto g_y_y_xz_xyyzz = cbuffer.data(dh_geom_11_off + 558 * ccomps * dcomps);

            auto g_y_y_xz_xyzzz = cbuffer.data(dh_geom_11_off + 559 * ccomps * dcomps);

            auto g_y_y_xz_xzzzz = cbuffer.data(dh_geom_11_off + 560 * ccomps * dcomps);

            auto g_y_y_xz_yyyyy = cbuffer.data(dh_geom_11_off + 561 * ccomps * dcomps);

            auto g_y_y_xz_yyyyz = cbuffer.data(dh_geom_11_off + 562 * ccomps * dcomps);

            auto g_y_y_xz_yyyzz = cbuffer.data(dh_geom_11_off + 563 * ccomps * dcomps);

            auto g_y_y_xz_yyzzz = cbuffer.data(dh_geom_11_off + 564 * ccomps * dcomps);

            auto g_y_y_xz_yzzzz = cbuffer.data(dh_geom_11_off + 565 * ccomps * dcomps);

            auto g_y_y_xz_zzzzz = cbuffer.data(dh_geom_11_off + 566 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_xz_xxxxx, g_y_y_xz_xxxxy, g_y_y_xz_xxxxz, g_y_y_xz_xxxyy, g_y_y_xz_xxxyz, g_y_y_xz_xxxzz, g_y_y_xz_xxyyy, g_y_y_xz_xxyyz, g_y_y_xz_xxyzz, g_y_y_xz_xxzzz, g_y_y_xz_xyyyy, g_y_y_xz_xyyyz, g_y_y_xz_xyyzz, g_y_y_xz_xyzzz, g_y_y_xz_xzzzz, g_y_y_xz_yyyyy, g_y_y_xz_yyyyz, g_y_y_xz_yyyzz, g_y_y_xz_yyzzz, g_y_y_xz_yzzzz, g_y_y_xz_zzzzz, g_y_y_z_xxxxx, g_y_y_z_xxxxxx, g_y_y_z_xxxxxy, g_y_y_z_xxxxxz, g_y_y_z_xxxxy, g_y_y_z_xxxxyy, g_y_y_z_xxxxyz, g_y_y_z_xxxxz, g_y_y_z_xxxxzz, g_y_y_z_xxxyy, g_y_y_z_xxxyyy, g_y_y_z_xxxyyz, g_y_y_z_xxxyz, g_y_y_z_xxxyzz, g_y_y_z_xxxzz, g_y_y_z_xxxzzz, g_y_y_z_xxyyy, g_y_y_z_xxyyyy, g_y_y_z_xxyyyz, g_y_y_z_xxyyz, g_y_y_z_xxyyzz, g_y_y_z_xxyzz, g_y_y_z_xxyzzz, g_y_y_z_xxzzz, g_y_y_z_xxzzzz, g_y_y_z_xyyyy, g_y_y_z_xyyyyy, g_y_y_z_xyyyyz, g_y_y_z_xyyyz, g_y_y_z_xyyyzz, g_y_y_z_xyyzz, g_y_y_z_xyyzzz, g_y_y_z_xyzzz, g_y_y_z_xyzzzz, g_y_y_z_xzzzz, g_y_y_z_xzzzzz, g_y_y_z_yyyyy, g_y_y_z_yyyyz, g_y_y_z_yyyzz, g_y_y_z_yyzzz, g_y_y_z_yzzzz, g_y_y_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_xz_xxxxx[k] = -g_y_y_z_xxxxx[k] * ab_x + g_y_y_z_xxxxxx[k];

                g_y_y_xz_xxxxy[k] = -g_y_y_z_xxxxy[k] * ab_x + g_y_y_z_xxxxxy[k];

                g_y_y_xz_xxxxz[k] = -g_y_y_z_xxxxz[k] * ab_x + g_y_y_z_xxxxxz[k];

                g_y_y_xz_xxxyy[k] = -g_y_y_z_xxxyy[k] * ab_x + g_y_y_z_xxxxyy[k];

                g_y_y_xz_xxxyz[k] = -g_y_y_z_xxxyz[k] * ab_x + g_y_y_z_xxxxyz[k];

                g_y_y_xz_xxxzz[k] = -g_y_y_z_xxxzz[k] * ab_x + g_y_y_z_xxxxzz[k];

                g_y_y_xz_xxyyy[k] = -g_y_y_z_xxyyy[k] * ab_x + g_y_y_z_xxxyyy[k];

                g_y_y_xz_xxyyz[k] = -g_y_y_z_xxyyz[k] * ab_x + g_y_y_z_xxxyyz[k];

                g_y_y_xz_xxyzz[k] = -g_y_y_z_xxyzz[k] * ab_x + g_y_y_z_xxxyzz[k];

                g_y_y_xz_xxzzz[k] = -g_y_y_z_xxzzz[k] * ab_x + g_y_y_z_xxxzzz[k];

                g_y_y_xz_xyyyy[k] = -g_y_y_z_xyyyy[k] * ab_x + g_y_y_z_xxyyyy[k];

                g_y_y_xz_xyyyz[k] = -g_y_y_z_xyyyz[k] * ab_x + g_y_y_z_xxyyyz[k];

                g_y_y_xz_xyyzz[k] = -g_y_y_z_xyyzz[k] * ab_x + g_y_y_z_xxyyzz[k];

                g_y_y_xz_xyzzz[k] = -g_y_y_z_xyzzz[k] * ab_x + g_y_y_z_xxyzzz[k];

                g_y_y_xz_xzzzz[k] = -g_y_y_z_xzzzz[k] * ab_x + g_y_y_z_xxzzzz[k];

                g_y_y_xz_yyyyy[k] = -g_y_y_z_yyyyy[k] * ab_x + g_y_y_z_xyyyyy[k];

                g_y_y_xz_yyyyz[k] = -g_y_y_z_yyyyz[k] * ab_x + g_y_y_z_xyyyyz[k];

                g_y_y_xz_yyyzz[k] = -g_y_y_z_yyyzz[k] * ab_x + g_y_y_z_xyyyzz[k];

                g_y_y_xz_yyzzz[k] = -g_y_y_z_yyzzz[k] * ab_x + g_y_y_z_xyyzzz[k];

                g_y_y_xz_yzzzz[k] = -g_y_y_z_yzzzz[k] * ab_x + g_y_y_z_xyzzzz[k];

                g_y_y_xz_zzzzz[k] = -g_y_y_z_zzzzz[k] * ab_x + g_y_y_z_xzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

            auto g_y_y_yy_xxxxx = cbuffer.data(dh_geom_11_off + 567 * ccomps * dcomps);

            auto g_y_y_yy_xxxxy = cbuffer.data(dh_geom_11_off + 568 * ccomps * dcomps);

            auto g_y_y_yy_xxxxz = cbuffer.data(dh_geom_11_off + 569 * ccomps * dcomps);

            auto g_y_y_yy_xxxyy = cbuffer.data(dh_geom_11_off + 570 * ccomps * dcomps);

            auto g_y_y_yy_xxxyz = cbuffer.data(dh_geom_11_off + 571 * ccomps * dcomps);

            auto g_y_y_yy_xxxzz = cbuffer.data(dh_geom_11_off + 572 * ccomps * dcomps);

            auto g_y_y_yy_xxyyy = cbuffer.data(dh_geom_11_off + 573 * ccomps * dcomps);

            auto g_y_y_yy_xxyyz = cbuffer.data(dh_geom_11_off + 574 * ccomps * dcomps);

            auto g_y_y_yy_xxyzz = cbuffer.data(dh_geom_11_off + 575 * ccomps * dcomps);

            auto g_y_y_yy_xxzzz = cbuffer.data(dh_geom_11_off + 576 * ccomps * dcomps);

            auto g_y_y_yy_xyyyy = cbuffer.data(dh_geom_11_off + 577 * ccomps * dcomps);

            auto g_y_y_yy_xyyyz = cbuffer.data(dh_geom_11_off + 578 * ccomps * dcomps);

            auto g_y_y_yy_xyyzz = cbuffer.data(dh_geom_11_off + 579 * ccomps * dcomps);

            auto g_y_y_yy_xyzzz = cbuffer.data(dh_geom_11_off + 580 * ccomps * dcomps);

            auto g_y_y_yy_xzzzz = cbuffer.data(dh_geom_11_off + 581 * ccomps * dcomps);

            auto g_y_y_yy_yyyyy = cbuffer.data(dh_geom_11_off + 582 * ccomps * dcomps);

            auto g_y_y_yy_yyyyz = cbuffer.data(dh_geom_11_off + 583 * ccomps * dcomps);

            auto g_y_y_yy_yyyzz = cbuffer.data(dh_geom_11_off + 584 * ccomps * dcomps);

            auto g_y_y_yy_yyzzz = cbuffer.data(dh_geom_11_off + 585 * ccomps * dcomps);

            auto g_y_y_yy_yzzzz = cbuffer.data(dh_geom_11_off + 586 * ccomps * dcomps);

            auto g_y_y_yy_zzzzz = cbuffer.data(dh_geom_11_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_y_xxxxx, g_0_y_y_xxxxy, g_0_y_y_xxxxz, g_0_y_y_xxxyy, g_0_y_y_xxxyz, g_0_y_y_xxxzz, g_0_y_y_xxyyy, g_0_y_y_xxyyz, g_0_y_y_xxyzz, g_0_y_y_xxzzz, g_0_y_y_xyyyy, g_0_y_y_xyyyz, g_0_y_y_xyyzz, g_0_y_y_xyzzz, g_0_y_y_xzzzz, g_0_y_y_yyyyy, g_0_y_y_yyyyz, g_0_y_y_yyyzz, g_0_y_y_yyzzz, g_0_y_y_yzzzz, g_0_y_y_zzzzz, g_y_0_y_xxxxx, g_y_0_y_xxxxy, g_y_0_y_xxxxz, g_y_0_y_xxxyy, g_y_0_y_xxxyz, g_y_0_y_xxxzz, g_y_0_y_xxyyy, g_y_0_y_xxyyz, g_y_0_y_xxyzz, g_y_0_y_xxzzz, g_y_0_y_xyyyy, g_y_0_y_xyyyz, g_y_0_y_xyyzz, g_y_0_y_xyzzz, g_y_0_y_xzzzz, g_y_0_y_yyyyy, g_y_0_y_yyyyz, g_y_0_y_yyyzz, g_y_0_y_yyzzz, g_y_0_y_yzzzz, g_y_0_y_zzzzz, g_y_y_y_xxxxx, g_y_y_y_xxxxxy, g_y_y_y_xxxxy, g_y_y_y_xxxxyy, g_y_y_y_xxxxyz, g_y_y_y_xxxxz, g_y_y_y_xxxyy, g_y_y_y_xxxyyy, g_y_y_y_xxxyyz, g_y_y_y_xxxyz, g_y_y_y_xxxyzz, g_y_y_y_xxxzz, g_y_y_y_xxyyy, g_y_y_y_xxyyyy, g_y_y_y_xxyyyz, g_y_y_y_xxyyz, g_y_y_y_xxyyzz, g_y_y_y_xxyzz, g_y_y_y_xxyzzz, g_y_y_y_xxzzz, g_y_y_y_xyyyy, g_y_y_y_xyyyyy, g_y_y_y_xyyyyz, g_y_y_y_xyyyz, g_y_y_y_xyyyzz, g_y_y_y_xyyzz, g_y_y_y_xyyzzz, g_y_y_y_xyzzz, g_y_y_y_xyzzzz, g_y_y_y_xzzzz, g_y_y_y_yyyyy, g_y_y_y_yyyyyy, g_y_y_y_yyyyyz, g_y_y_y_yyyyz, g_y_y_y_yyyyzz, g_y_y_y_yyyzz, g_y_y_y_yyyzzz, g_y_y_y_yyzzz, g_y_y_y_yyzzzz, g_y_y_y_yzzzz, g_y_y_y_yzzzzz, g_y_y_y_zzzzz, g_y_y_yy_xxxxx, g_y_y_yy_xxxxy, g_y_y_yy_xxxxz, g_y_y_yy_xxxyy, g_y_y_yy_xxxyz, g_y_y_yy_xxxzz, g_y_y_yy_xxyyy, g_y_y_yy_xxyyz, g_y_y_yy_xxyzz, g_y_y_yy_xxzzz, g_y_y_yy_xyyyy, g_y_y_yy_xyyyz, g_y_y_yy_xyyzz, g_y_y_yy_xyzzz, g_y_y_yy_xzzzz, g_y_y_yy_yyyyy, g_y_y_yy_yyyyz, g_y_y_yy_yyyzz, g_y_y_yy_yyzzz, g_y_y_yy_yzzzz, g_y_y_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yy_xxxxx[k] = -g_0_y_y_xxxxx[k] + g_y_0_y_xxxxx[k] - g_y_y_y_xxxxx[k] * ab_y + g_y_y_y_xxxxxy[k];

                g_y_y_yy_xxxxy[k] = -g_0_y_y_xxxxy[k] + g_y_0_y_xxxxy[k] - g_y_y_y_xxxxy[k] * ab_y + g_y_y_y_xxxxyy[k];

                g_y_y_yy_xxxxz[k] = -g_0_y_y_xxxxz[k] + g_y_0_y_xxxxz[k] - g_y_y_y_xxxxz[k] * ab_y + g_y_y_y_xxxxyz[k];

                g_y_y_yy_xxxyy[k] = -g_0_y_y_xxxyy[k] + g_y_0_y_xxxyy[k] - g_y_y_y_xxxyy[k] * ab_y + g_y_y_y_xxxyyy[k];

                g_y_y_yy_xxxyz[k] = -g_0_y_y_xxxyz[k] + g_y_0_y_xxxyz[k] - g_y_y_y_xxxyz[k] * ab_y + g_y_y_y_xxxyyz[k];

                g_y_y_yy_xxxzz[k] = -g_0_y_y_xxxzz[k] + g_y_0_y_xxxzz[k] - g_y_y_y_xxxzz[k] * ab_y + g_y_y_y_xxxyzz[k];

                g_y_y_yy_xxyyy[k] = -g_0_y_y_xxyyy[k] + g_y_0_y_xxyyy[k] - g_y_y_y_xxyyy[k] * ab_y + g_y_y_y_xxyyyy[k];

                g_y_y_yy_xxyyz[k] = -g_0_y_y_xxyyz[k] + g_y_0_y_xxyyz[k] - g_y_y_y_xxyyz[k] * ab_y + g_y_y_y_xxyyyz[k];

                g_y_y_yy_xxyzz[k] = -g_0_y_y_xxyzz[k] + g_y_0_y_xxyzz[k] - g_y_y_y_xxyzz[k] * ab_y + g_y_y_y_xxyyzz[k];

                g_y_y_yy_xxzzz[k] = -g_0_y_y_xxzzz[k] + g_y_0_y_xxzzz[k] - g_y_y_y_xxzzz[k] * ab_y + g_y_y_y_xxyzzz[k];

                g_y_y_yy_xyyyy[k] = -g_0_y_y_xyyyy[k] + g_y_0_y_xyyyy[k] - g_y_y_y_xyyyy[k] * ab_y + g_y_y_y_xyyyyy[k];

                g_y_y_yy_xyyyz[k] = -g_0_y_y_xyyyz[k] + g_y_0_y_xyyyz[k] - g_y_y_y_xyyyz[k] * ab_y + g_y_y_y_xyyyyz[k];

                g_y_y_yy_xyyzz[k] = -g_0_y_y_xyyzz[k] + g_y_0_y_xyyzz[k] - g_y_y_y_xyyzz[k] * ab_y + g_y_y_y_xyyyzz[k];

                g_y_y_yy_xyzzz[k] = -g_0_y_y_xyzzz[k] + g_y_0_y_xyzzz[k] - g_y_y_y_xyzzz[k] * ab_y + g_y_y_y_xyyzzz[k];

                g_y_y_yy_xzzzz[k] = -g_0_y_y_xzzzz[k] + g_y_0_y_xzzzz[k] - g_y_y_y_xzzzz[k] * ab_y + g_y_y_y_xyzzzz[k];

                g_y_y_yy_yyyyy[k] = -g_0_y_y_yyyyy[k] + g_y_0_y_yyyyy[k] - g_y_y_y_yyyyy[k] * ab_y + g_y_y_y_yyyyyy[k];

                g_y_y_yy_yyyyz[k] = -g_0_y_y_yyyyz[k] + g_y_0_y_yyyyz[k] - g_y_y_y_yyyyz[k] * ab_y + g_y_y_y_yyyyyz[k];

                g_y_y_yy_yyyzz[k] = -g_0_y_y_yyyzz[k] + g_y_0_y_yyyzz[k] - g_y_y_y_yyyzz[k] * ab_y + g_y_y_y_yyyyzz[k];

                g_y_y_yy_yyzzz[k] = -g_0_y_y_yyzzz[k] + g_y_0_y_yyzzz[k] - g_y_y_y_yyzzz[k] * ab_y + g_y_y_y_yyyzzz[k];

                g_y_y_yy_yzzzz[k] = -g_0_y_y_yzzzz[k] + g_y_0_y_yzzzz[k] - g_y_y_y_yzzzz[k] * ab_y + g_y_y_y_yyzzzz[k];

                g_y_y_yy_zzzzz[k] = -g_0_y_y_zzzzz[k] + g_y_0_y_zzzzz[k] - g_y_y_y_zzzzz[k] * ab_y + g_y_y_y_yzzzzz[k];
            }

            /// Set up 588-609 components of targeted buffer : cbuffer.data(

            auto g_y_y_yz_xxxxx = cbuffer.data(dh_geom_11_off + 588 * ccomps * dcomps);

            auto g_y_y_yz_xxxxy = cbuffer.data(dh_geom_11_off + 589 * ccomps * dcomps);

            auto g_y_y_yz_xxxxz = cbuffer.data(dh_geom_11_off + 590 * ccomps * dcomps);

            auto g_y_y_yz_xxxyy = cbuffer.data(dh_geom_11_off + 591 * ccomps * dcomps);

            auto g_y_y_yz_xxxyz = cbuffer.data(dh_geom_11_off + 592 * ccomps * dcomps);

            auto g_y_y_yz_xxxzz = cbuffer.data(dh_geom_11_off + 593 * ccomps * dcomps);

            auto g_y_y_yz_xxyyy = cbuffer.data(dh_geom_11_off + 594 * ccomps * dcomps);

            auto g_y_y_yz_xxyyz = cbuffer.data(dh_geom_11_off + 595 * ccomps * dcomps);

            auto g_y_y_yz_xxyzz = cbuffer.data(dh_geom_11_off + 596 * ccomps * dcomps);

            auto g_y_y_yz_xxzzz = cbuffer.data(dh_geom_11_off + 597 * ccomps * dcomps);

            auto g_y_y_yz_xyyyy = cbuffer.data(dh_geom_11_off + 598 * ccomps * dcomps);

            auto g_y_y_yz_xyyyz = cbuffer.data(dh_geom_11_off + 599 * ccomps * dcomps);

            auto g_y_y_yz_xyyzz = cbuffer.data(dh_geom_11_off + 600 * ccomps * dcomps);

            auto g_y_y_yz_xyzzz = cbuffer.data(dh_geom_11_off + 601 * ccomps * dcomps);

            auto g_y_y_yz_xzzzz = cbuffer.data(dh_geom_11_off + 602 * ccomps * dcomps);

            auto g_y_y_yz_yyyyy = cbuffer.data(dh_geom_11_off + 603 * ccomps * dcomps);

            auto g_y_y_yz_yyyyz = cbuffer.data(dh_geom_11_off + 604 * ccomps * dcomps);

            auto g_y_y_yz_yyyzz = cbuffer.data(dh_geom_11_off + 605 * ccomps * dcomps);

            auto g_y_y_yz_yyzzz = cbuffer.data(dh_geom_11_off + 606 * ccomps * dcomps);

            auto g_y_y_yz_yzzzz = cbuffer.data(dh_geom_11_off + 607 * ccomps * dcomps);

            auto g_y_y_yz_zzzzz = cbuffer.data(dh_geom_11_off + 608 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_y_xxxxx, g_y_y_y_xxxxxz, g_y_y_y_xxxxy, g_y_y_y_xxxxyz, g_y_y_y_xxxxz, g_y_y_y_xxxxzz, g_y_y_y_xxxyy, g_y_y_y_xxxyyz, g_y_y_y_xxxyz, g_y_y_y_xxxyzz, g_y_y_y_xxxzz, g_y_y_y_xxxzzz, g_y_y_y_xxyyy, g_y_y_y_xxyyyz, g_y_y_y_xxyyz, g_y_y_y_xxyyzz, g_y_y_y_xxyzz, g_y_y_y_xxyzzz, g_y_y_y_xxzzz, g_y_y_y_xxzzzz, g_y_y_y_xyyyy, g_y_y_y_xyyyyz, g_y_y_y_xyyyz, g_y_y_y_xyyyzz, g_y_y_y_xyyzz, g_y_y_y_xyyzzz, g_y_y_y_xyzzz, g_y_y_y_xyzzzz, g_y_y_y_xzzzz, g_y_y_y_xzzzzz, g_y_y_y_yyyyy, g_y_y_y_yyyyyz, g_y_y_y_yyyyz, g_y_y_y_yyyyzz, g_y_y_y_yyyzz, g_y_y_y_yyyzzz, g_y_y_y_yyzzz, g_y_y_y_yyzzzz, g_y_y_y_yzzzz, g_y_y_y_yzzzzz, g_y_y_y_zzzzz, g_y_y_y_zzzzzz, g_y_y_yz_xxxxx, g_y_y_yz_xxxxy, g_y_y_yz_xxxxz, g_y_y_yz_xxxyy, g_y_y_yz_xxxyz, g_y_y_yz_xxxzz, g_y_y_yz_xxyyy, g_y_y_yz_xxyyz, g_y_y_yz_xxyzz, g_y_y_yz_xxzzz, g_y_y_yz_xyyyy, g_y_y_yz_xyyyz, g_y_y_yz_xyyzz, g_y_y_yz_xyzzz, g_y_y_yz_xzzzz, g_y_y_yz_yyyyy, g_y_y_yz_yyyyz, g_y_y_yz_yyyzz, g_y_y_yz_yyzzz, g_y_y_yz_yzzzz, g_y_y_yz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_yz_xxxxx[k] = -g_y_y_y_xxxxx[k] * ab_z + g_y_y_y_xxxxxz[k];

                g_y_y_yz_xxxxy[k] = -g_y_y_y_xxxxy[k] * ab_z + g_y_y_y_xxxxyz[k];

                g_y_y_yz_xxxxz[k] = -g_y_y_y_xxxxz[k] * ab_z + g_y_y_y_xxxxzz[k];

                g_y_y_yz_xxxyy[k] = -g_y_y_y_xxxyy[k] * ab_z + g_y_y_y_xxxyyz[k];

                g_y_y_yz_xxxyz[k] = -g_y_y_y_xxxyz[k] * ab_z + g_y_y_y_xxxyzz[k];

                g_y_y_yz_xxxzz[k] = -g_y_y_y_xxxzz[k] * ab_z + g_y_y_y_xxxzzz[k];

                g_y_y_yz_xxyyy[k] = -g_y_y_y_xxyyy[k] * ab_z + g_y_y_y_xxyyyz[k];

                g_y_y_yz_xxyyz[k] = -g_y_y_y_xxyyz[k] * ab_z + g_y_y_y_xxyyzz[k];

                g_y_y_yz_xxyzz[k] = -g_y_y_y_xxyzz[k] * ab_z + g_y_y_y_xxyzzz[k];

                g_y_y_yz_xxzzz[k] = -g_y_y_y_xxzzz[k] * ab_z + g_y_y_y_xxzzzz[k];

                g_y_y_yz_xyyyy[k] = -g_y_y_y_xyyyy[k] * ab_z + g_y_y_y_xyyyyz[k];

                g_y_y_yz_xyyyz[k] = -g_y_y_y_xyyyz[k] * ab_z + g_y_y_y_xyyyzz[k];

                g_y_y_yz_xyyzz[k] = -g_y_y_y_xyyzz[k] * ab_z + g_y_y_y_xyyzzz[k];

                g_y_y_yz_xyzzz[k] = -g_y_y_y_xyzzz[k] * ab_z + g_y_y_y_xyzzzz[k];

                g_y_y_yz_xzzzz[k] = -g_y_y_y_xzzzz[k] * ab_z + g_y_y_y_xzzzzz[k];

                g_y_y_yz_yyyyy[k] = -g_y_y_y_yyyyy[k] * ab_z + g_y_y_y_yyyyyz[k];

                g_y_y_yz_yyyyz[k] = -g_y_y_y_yyyyz[k] * ab_z + g_y_y_y_yyyyzz[k];

                g_y_y_yz_yyyzz[k] = -g_y_y_y_yyyzz[k] * ab_z + g_y_y_y_yyyzzz[k];

                g_y_y_yz_yyzzz[k] = -g_y_y_y_yyzzz[k] * ab_z + g_y_y_y_yyzzzz[k];

                g_y_y_yz_yzzzz[k] = -g_y_y_y_yzzzz[k] * ab_z + g_y_y_y_yzzzzz[k];

                g_y_y_yz_zzzzz[k] = -g_y_y_y_zzzzz[k] * ab_z + g_y_y_y_zzzzzz[k];
            }

            /// Set up 609-630 components of targeted buffer : cbuffer.data(

            auto g_y_y_zz_xxxxx = cbuffer.data(dh_geom_11_off + 609 * ccomps * dcomps);

            auto g_y_y_zz_xxxxy = cbuffer.data(dh_geom_11_off + 610 * ccomps * dcomps);

            auto g_y_y_zz_xxxxz = cbuffer.data(dh_geom_11_off + 611 * ccomps * dcomps);

            auto g_y_y_zz_xxxyy = cbuffer.data(dh_geom_11_off + 612 * ccomps * dcomps);

            auto g_y_y_zz_xxxyz = cbuffer.data(dh_geom_11_off + 613 * ccomps * dcomps);

            auto g_y_y_zz_xxxzz = cbuffer.data(dh_geom_11_off + 614 * ccomps * dcomps);

            auto g_y_y_zz_xxyyy = cbuffer.data(dh_geom_11_off + 615 * ccomps * dcomps);

            auto g_y_y_zz_xxyyz = cbuffer.data(dh_geom_11_off + 616 * ccomps * dcomps);

            auto g_y_y_zz_xxyzz = cbuffer.data(dh_geom_11_off + 617 * ccomps * dcomps);

            auto g_y_y_zz_xxzzz = cbuffer.data(dh_geom_11_off + 618 * ccomps * dcomps);

            auto g_y_y_zz_xyyyy = cbuffer.data(dh_geom_11_off + 619 * ccomps * dcomps);

            auto g_y_y_zz_xyyyz = cbuffer.data(dh_geom_11_off + 620 * ccomps * dcomps);

            auto g_y_y_zz_xyyzz = cbuffer.data(dh_geom_11_off + 621 * ccomps * dcomps);

            auto g_y_y_zz_xyzzz = cbuffer.data(dh_geom_11_off + 622 * ccomps * dcomps);

            auto g_y_y_zz_xzzzz = cbuffer.data(dh_geom_11_off + 623 * ccomps * dcomps);

            auto g_y_y_zz_yyyyy = cbuffer.data(dh_geom_11_off + 624 * ccomps * dcomps);

            auto g_y_y_zz_yyyyz = cbuffer.data(dh_geom_11_off + 625 * ccomps * dcomps);

            auto g_y_y_zz_yyyzz = cbuffer.data(dh_geom_11_off + 626 * ccomps * dcomps);

            auto g_y_y_zz_yyzzz = cbuffer.data(dh_geom_11_off + 627 * ccomps * dcomps);

            auto g_y_y_zz_yzzzz = cbuffer.data(dh_geom_11_off + 628 * ccomps * dcomps);

            auto g_y_y_zz_zzzzz = cbuffer.data(dh_geom_11_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_y_z_xxxxx, g_y_y_z_xxxxxz, g_y_y_z_xxxxy, g_y_y_z_xxxxyz, g_y_y_z_xxxxz, g_y_y_z_xxxxzz, g_y_y_z_xxxyy, g_y_y_z_xxxyyz, g_y_y_z_xxxyz, g_y_y_z_xxxyzz, g_y_y_z_xxxzz, g_y_y_z_xxxzzz, g_y_y_z_xxyyy, g_y_y_z_xxyyyz, g_y_y_z_xxyyz, g_y_y_z_xxyyzz, g_y_y_z_xxyzz, g_y_y_z_xxyzzz, g_y_y_z_xxzzz, g_y_y_z_xxzzzz, g_y_y_z_xyyyy, g_y_y_z_xyyyyz, g_y_y_z_xyyyz, g_y_y_z_xyyyzz, g_y_y_z_xyyzz, g_y_y_z_xyyzzz, g_y_y_z_xyzzz, g_y_y_z_xyzzzz, g_y_y_z_xzzzz, g_y_y_z_xzzzzz, g_y_y_z_yyyyy, g_y_y_z_yyyyyz, g_y_y_z_yyyyz, g_y_y_z_yyyyzz, g_y_y_z_yyyzz, g_y_y_z_yyyzzz, g_y_y_z_yyzzz, g_y_y_z_yyzzzz, g_y_y_z_yzzzz, g_y_y_z_yzzzzz, g_y_y_z_zzzzz, g_y_y_z_zzzzzz, g_y_y_zz_xxxxx, g_y_y_zz_xxxxy, g_y_y_zz_xxxxz, g_y_y_zz_xxxyy, g_y_y_zz_xxxyz, g_y_y_zz_xxxzz, g_y_y_zz_xxyyy, g_y_y_zz_xxyyz, g_y_y_zz_xxyzz, g_y_y_zz_xxzzz, g_y_y_zz_xyyyy, g_y_y_zz_xyyyz, g_y_y_zz_xyyzz, g_y_y_zz_xyzzz, g_y_y_zz_xzzzz, g_y_y_zz_yyyyy, g_y_y_zz_yyyyz, g_y_y_zz_yyyzz, g_y_y_zz_yyzzz, g_y_y_zz_yzzzz, g_y_y_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_y_zz_xxxxx[k] = -g_y_y_z_xxxxx[k] * ab_z + g_y_y_z_xxxxxz[k];

                g_y_y_zz_xxxxy[k] = -g_y_y_z_xxxxy[k] * ab_z + g_y_y_z_xxxxyz[k];

                g_y_y_zz_xxxxz[k] = -g_y_y_z_xxxxz[k] * ab_z + g_y_y_z_xxxxzz[k];

                g_y_y_zz_xxxyy[k] = -g_y_y_z_xxxyy[k] * ab_z + g_y_y_z_xxxyyz[k];

                g_y_y_zz_xxxyz[k] = -g_y_y_z_xxxyz[k] * ab_z + g_y_y_z_xxxyzz[k];

                g_y_y_zz_xxxzz[k] = -g_y_y_z_xxxzz[k] * ab_z + g_y_y_z_xxxzzz[k];

                g_y_y_zz_xxyyy[k] = -g_y_y_z_xxyyy[k] * ab_z + g_y_y_z_xxyyyz[k];

                g_y_y_zz_xxyyz[k] = -g_y_y_z_xxyyz[k] * ab_z + g_y_y_z_xxyyzz[k];

                g_y_y_zz_xxyzz[k] = -g_y_y_z_xxyzz[k] * ab_z + g_y_y_z_xxyzzz[k];

                g_y_y_zz_xxzzz[k] = -g_y_y_z_xxzzz[k] * ab_z + g_y_y_z_xxzzzz[k];

                g_y_y_zz_xyyyy[k] = -g_y_y_z_xyyyy[k] * ab_z + g_y_y_z_xyyyyz[k];

                g_y_y_zz_xyyyz[k] = -g_y_y_z_xyyyz[k] * ab_z + g_y_y_z_xyyyzz[k];

                g_y_y_zz_xyyzz[k] = -g_y_y_z_xyyzz[k] * ab_z + g_y_y_z_xyyzzz[k];

                g_y_y_zz_xyzzz[k] = -g_y_y_z_xyzzz[k] * ab_z + g_y_y_z_xyzzzz[k];

                g_y_y_zz_xzzzz[k] = -g_y_y_z_xzzzz[k] * ab_z + g_y_y_z_xzzzzz[k];

                g_y_y_zz_yyyyy[k] = -g_y_y_z_yyyyy[k] * ab_z + g_y_y_z_yyyyyz[k];

                g_y_y_zz_yyyyz[k] = -g_y_y_z_yyyyz[k] * ab_z + g_y_y_z_yyyyzz[k];

                g_y_y_zz_yyyzz[k] = -g_y_y_z_yyyzz[k] * ab_z + g_y_y_z_yyyzzz[k];

                g_y_y_zz_yyzzz[k] = -g_y_y_z_yyzzz[k] * ab_z + g_y_y_z_yyzzzz[k];

                g_y_y_zz_yzzzz[k] = -g_y_y_z_yzzzz[k] * ab_z + g_y_y_z_yzzzzz[k];

                g_y_y_zz_zzzzz[k] = -g_y_y_z_zzzzz[k] * ab_z + g_y_y_z_zzzzzz[k];
            }

            /// Set up 630-651 components of targeted buffer : cbuffer.data(

            auto g_y_z_xx_xxxxx = cbuffer.data(dh_geom_11_off + 630 * ccomps * dcomps);

            auto g_y_z_xx_xxxxy = cbuffer.data(dh_geom_11_off + 631 * ccomps * dcomps);

            auto g_y_z_xx_xxxxz = cbuffer.data(dh_geom_11_off + 632 * ccomps * dcomps);

            auto g_y_z_xx_xxxyy = cbuffer.data(dh_geom_11_off + 633 * ccomps * dcomps);

            auto g_y_z_xx_xxxyz = cbuffer.data(dh_geom_11_off + 634 * ccomps * dcomps);

            auto g_y_z_xx_xxxzz = cbuffer.data(dh_geom_11_off + 635 * ccomps * dcomps);

            auto g_y_z_xx_xxyyy = cbuffer.data(dh_geom_11_off + 636 * ccomps * dcomps);

            auto g_y_z_xx_xxyyz = cbuffer.data(dh_geom_11_off + 637 * ccomps * dcomps);

            auto g_y_z_xx_xxyzz = cbuffer.data(dh_geom_11_off + 638 * ccomps * dcomps);

            auto g_y_z_xx_xxzzz = cbuffer.data(dh_geom_11_off + 639 * ccomps * dcomps);

            auto g_y_z_xx_xyyyy = cbuffer.data(dh_geom_11_off + 640 * ccomps * dcomps);

            auto g_y_z_xx_xyyyz = cbuffer.data(dh_geom_11_off + 641 * ccomps * dcomps);

            auto g_y_z_xx_xyyzz = cbuffer.data(dh_geom_11_off + 642 * ccomps * dcomps);

            auto g_y_z_xx_xyzzz = cbuffer.data(dh_geom_11_off + 643 * ccomps * dcomps);

            auto g_y_z_xx_xzzzz = cbuffer.data(dh_geom_11_off + 644 * ccomps * dcomps);

            auto g_y_z_xx_yyyyy = cbuffer.data(dh_geom_11_off + 645 * ccomps * dcomps);

            auto g_y_z_xx_yyyyz = cbuffer.data(dh_geom_11_off + 646 * ccomps * dcomps);

            auto g_y_z_xx_yyyzz = cbuffer.data(dh_geom_11_off + 647 * ccomps * dcomps);

            auto g_y_z_xx_yyzzz = cbuffer.data(dh_geom_11_off + 648 * ccomps * dcomps);

            auto g_y_z_xx_yzzzz = cbuffer.data(dh_geom_11_off + 649 * ccomps * dcomps);

            auto g_y_z_xx_zzzzz = cbuffer.data(dh_geom_11_off + 650 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_x_xxxxx, g_y_z_x_xxxxxx, g_y_z_x_xxxxxy, g_y_z_x_xxxxxz, g_y_z_x_xxxxy, g_y_z_x_xxxxyy, g_y_z_x_xxxxyz, g_y_z_x_xxxxz, g_y_z_x_xxxxzz, g_y_z_x_xxxyy, g_y_z_x_xxxyyy, g_y_z_x_xxxyyz, g_y_z_x_xxxyz, g_y_z_x_xxxyzz, g_y_z_x_xxxzz, g_y_z_x_xxxzzz, g_y_z_x_xxyyy, g_y_z_x_xxyyyy, g_y_z_x_xxyyyz, g_y_z_x_xxyyz, g_y_z_x_xxyyzz, g_y_z_x_xxyzz, g_y_z_x_xxyzzz, g_y_z_x_xxzzz, g_y_z_x_xxzzzz, g_y_z_x_xyyyy, g_y_z_x_xyyyyy, g_y_z_x_xyyyyz, g_y_z_x_xyyyz, g_y_z_x_xyyyzz, g_y_z_x_xyyzz, g_y_z_x_xyyzzz, g_y_z_x_xyzzz, g_y_z_x_xyzzzz, g_y_z_x_xzzzz, g_y_z_x_xzzzzz, g_y_z_x_yyyyy, g_y_z_x_yyyyz, g_y_z_x_yyyzz, g_y_z_x_yyzzz, g_y_z_x_yzzzz, g_y_z_x_zzzzz, g_y_z_xx_xxxxx, g_y_z_xx_xxxxy, g_y_z_xx_xxxxz, g_y_z_xx_xxxyy, g_y_z_xx_xxxyz, g_y_z_xx_xxxzz, g_y_z_xx_xxyyy, g_y_z_xx_xxyyz, g_y_z_xx_xxyzz, g_y_z_xx_xxzzz, g_y_z_xx_xyyyy, g_y_z_xx_xyyyz, g_y_z_xx_xyyzz, g_y_z_xx_xyzzz, g_y_z_xx_xzzzz, g_y_z_xx_yyyyy, g_y_z_xx_yyyyz, g_y_z_xx_yyyzz, g_y_z_xx_yyzzz, g_y_z_xx_yzzzz, g_y_z_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xx_xxxxx[k] = -g_y_z_x_xxxxx[k] * ab_x + g_y_z_x_xxxxxx[k];

                g_y_z_xx_xxxxy[k] = -g_y_z_x_xxxxy[k] * ab_x + g_y_z_x_xxxxxy[k];

                g_y_z_xx_xxxxz[k] = -g_y_z_x_xxxxz[k] * ab_x + g_y_z_x_xxxxxz[k];

                g_y_z_xx_xxxyy[k] = -g_y_z_x_xxxyy[k] * ab_x + g_y_z_x_xxxxyy[k];

                g_y_z_xx_xxxyz[k] = -g_y_z_x_xxxyz[k] * ab_x + g_y_z_x_xxxxyz[k];

                g_y_z_xx_xxxzz[k] = -g_y_z_x_xxxzz[k] * ab_x + g_y_z_x_xxxxzz[k];

                g_y_z_xx_xxyyy[k] = -g_y_z_x_xxyyy[k] * ab_x + g_y_z_x_xxxyyy[k];

                g_y_z_xx_xxyyz[k] = -g_y_z_x_xxyyz[k] * ab_x + g_y_z_x_xxxyyz[k];

                g_y_z_xx_xxyzz[k] = -g_y_z_x_xxyzz[k] * ab_x + g_y_z_x_xxxyzz[k];

                g_y_z_xx_xxzzz[k] = -g_y_z_x_xxzzz[k] * ab_x + g_y_z_x_xxxzzz[k];

                g_y_z_xx_xyyyy[k] = -g_y_z_x_xyyyy[k] * ab_x + g_y_z_x_xxyyyy[k];

                g_y_z_xx_xyyyz[k] = -g_y_z_x_xyyyz[k] * ab_x + g_y_z_x_xxyyyz[k];

                g_y_z_xx_xyyzz[k] = -g_y_z_x_xyyzz[k] * ab_x + g_y_z_x_xxyyzz[k];

                g_y_z_xx_xyzzz[k] = -g_y_z_x_xyzzz[k] * ab_x + g_y_z_x_xxyzzz[k];

                g_y_z_xx_xzzzz[k] = -g_y_z_x_xzzzz[k] * ab_x + g_y_z_x_xxzzzz[k];

                g_y_z_xx_yyyyy[k] = -g_y_z_x_yyyyy[k] * ab_x + g_y_z_x_xyyyyy[k];

                g_y_z_xx_yyyyz[k] = -g_y_z_x_yyyyz[k] * ab_x + g_y_z_x_xyyyyz[k];

                g_y_z_xx_yyyzz[k] = -g_y_z_x_yyyzz[k] * ab_x + g_y_z_x_xyyyzz[k];

                g_y_z_xx_yyzzz[k] = -g_y_z_x_yyzzz[k] * ab_x + g_y_z_x_xyyzzz[k];

                g_y_z_xx_yzzzz[k] = -g_y_z_x_yzzzz[k] * ab_x + g_y_z_x_xyzzzz[k];

                g_y_z_xx_zzzzz[k] = -g_y_z_x_zzzzz[k] * ab_x + g_y_z_x_xzzzzz[k];
            }

            /// Set up 651-672 components of targeted buffer : cbuffer.data(

            auto g_y_z_xy_xxxxx = cbuffer.data(dh_geom_11_off + 651 * ccomps * dcomps);

            auto g_y_z_xy_xxxxy = cbuffer.data(dh_geom_11_off + 652 * ccomps * dcomps);

            auto g_y_z_xy_xxxxz = cbuffer.data(dh_geom_11_off + 653 * ccomps * dcomps);

            auto g_y_z_xy_xxxyy = cbuffer.data(dh_geom_11_off + 654 * ccomps * dcomps);

            auto g_y_z_xy_xxxyz = cbuffer.data(dh_geom_11_off + 655 * ccomps * dcomps);

            auto g_y_z_xy_xxxzz = cbuffer.data(dh_geom_11_off + 656 * ccomps * dcomps);

            auto g_y_z_xy_xxyyy = cbuffer.data(dh_geom_11_off + 657 * ccomps * dcomps);

            auto g_y_z_xy_xxyyz = cbuffer.data(dh_geom_11_off + 658 * ccomps * dcomps);

            auto g_y_z_xy_xxyzz = cbuffer.data(dh_geom_11_off + 659 * ccomps * dcomps);

            auto g_y_z_xy_xxzzz = cbuffer.data(dh_geom_11_off + 660 * ccomps * dcomps);

            auto g_y_z_xy_xyyyy = cbuffer.data(dh_geom_11_off + 661 * ccomps * dcomps);

            auto g_y_z_xy_xyyyz = cbuffer.data(dh_geom_11_off + 662 * ccomps * dcomps);

            auto g_y_z_xy_xyyzz = cbuffer.data(dh_geom_11_off + 663 * ccomps * dcomps);

            auto g_y_z_xy_xyzzz = cbuffer.data(dh_geom_11_off + 664 * ccomps * dcomps);

            auto g_y_z_xy_xzzzz = cbuffer.data(dh_geom_11_off + 665 * ccomps * dcomps);

            auto g_y_z_xy_yyyyy = cbuffer.data(dh_geom_11_off + 666 * ccomps * dcomps);

            auto g_y_z_xy_yyyyz = cbuffer.data(dh_geom_11_off + 667 * ccomps * dcomps);

            auto g_y_z_xy_yyyzz = cbuffer.data(dh_geom_11_off + 668 * ccomps * dcomps);

            auto g_y_z_xy_yyzzz = cbuffer.data(dh_geom_11_off + 669 * ccomps * dcomps);

            auto g_y_z_xy_yzzzz = cbuffer.data(dh_geom_11_off + 670 * ccomps * dcomps);

            auto g_y_z_xy_zzzzz = cbuffer.data(dh_geom_11_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xy_xxxxx, g_y_z_xy_xxxxy, g_y_z_xy_xxxxz, g_y_z_xy_xxxyy, g_y_z_xy_xxxyz, g_y_z_xy_xxxzz, g_y_z_xy_xxyyy, g_y_z_xy_xxyyz, g_y_z_xy_xxyzz, g_y_z_xy_xxzzz, g_y_z_xy_xyyyy, g_y_z_xy_xyyyz, g_y_z_xy_xyyzz, g_y_z_xy_xyzzz, g_y_z_xy_xzzzz, g_y_z_xy_yyyyy, g_y_z_xy_yyyyz, g_y_z_xy_yyyzz, g_y_z_xy_yyzzz, g_y_z_xy_yzzzz, g_y_z_xy_zzzzz, g_y_z_y_xxxxx, g_y_z_y_xxxxxx, g_y_z_y_xxxxxy, g_y_z_y_xxxxxz, g_y_z_y_xxxxy, g_y_z_y_xxxxyy, g_y_z_y_xxxxyz, g_y_z_y_xxxxz, g_y_z_y_xxxxzz, g_y_z_y_xxxyy, g_y_z_y_xxxyyy, g_y_z_y_xxxyyz, g_y_z_y_xxxyz, g_y_z_y_xxxyzz, g_y_z_y_xxxzz, g_y_z_y_xxxzzz, g_y_z_y_xxyyy, g_y_z_y_xxyyyy, g_y_z_y_xxyyyz, g_y_z_y_xxyyz, g_y_z_y_xxyyzz, g_y_z_y_xxyzz, g_y_z_y_xxyzzz, g_y_z_y_xxzzz, g_y_z_y_xxzzzz, g_y_z_y_xyyyy, g_y_z_y_xyyyyy, g_y_z_y_xyyyyz, g_y_z_y_xyyyz, g_y_z_y_xyyyzz, g_y_z_y_xyyzz, g_y_z_y_xyyzzz, g_y_z_y_xyzzz, g_y_z_y_xyzzzz, g_y_z_y_xzzzz, g_y_z_y_xzzzzz, g_y_z_y_yyyyy, g_y_z_y_yyyyz, g_y_z_y_yyyzz, g_y_z_y_yyzzz, g_y_z_y_yzzzz, g_y_z_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xy_xxxxx[k] = -g_y_z_y_xxxxx[k] * ab_x + g_y_z_y_xxxxxx[k];

                g_y_z_xy_xxxxy[k] = -g_y_z_y_xxxxy[k] * ab_x + g_y_z_y_xxxxxy[k];

                g_y_z_xy_xxxxz[k] = -g_y_z_y_xxxxz[k] * ab_x + g_y_z_y_xxxxxz[k];

                g_y_z_xy_xxxyy[k] = -g_y_z_y_xxxyy[k] * ab_x + g_y_z_y_xxxxyy[k];

                g_y_z_xy_xxxyz[k] = -g_y_z_y_xxxyz[k] * ab_x + g_y_z_y_xxxxyz[k];

                g_y_z_xy_xxxzz[k] = -g_y_z_y_xxxzz[k] * ab_x + g_y_z_y_xxxxzz[k];

                g_y_z_xy_xxyyy[k] = -g_y_z_y_xxyyy[k] * ab_x + g_y_z_y_xxxyyy[k];

                g_y_z_xy_xxyyz[k] = -g_y_z_y_xxyyz[k] * ab_x + g_y_z_y_xxxyyz[k];

                g_y_z_xy_xxyzz[k] = -g_y_z_y_xxyzz[k] * ab_x + g_y_z_y_xxxyzz[k];

                g_y_z_xy_xxzzz[k] = -g_y_z_y_xxzzz[k] * ab_x + g_y_z_y_xxxzzz[k];

                g_y_z_xy_xyyyy[k] = -g_y_z_y_xyyyy[k] * ab_x + g_y_z_y_xxyyyy[k];

                g_y_z_xy_xyyyz[k] = -g_y_z_y_xyyyz[k] * ab_x + g_y_z_y_xxyyyz[k];

                g_y_z_xy_xyyzz[k] = -g_y_z_y_xyyzz[k] * ab_x + g_y_z_y_xxyyzz[k];

                g_y_z_xy_xyzzz[k] = -g_y_z_y_xyzzz[k] * ab_x + g_y_z_y_xxyzzz[k];

                g_y_z_xy_xzzzz[k] = -g_y_z_y_xzzzz[k] * ab_x + g_y_z_y_xxzzzz[k];

                g_y_z_xy_yyyyy[k] = -g_y_z_y_yyyyy[k] * ab_x + g_y_z_y_xyyyyy[k];

                g_y_z_xy_yyyyz[k] = -g_y_z_y_yyyyz[k] * ab_x + g_y_z_y_xyyyyz[k];

                g_y_z_xy_yyyzz[k] = -g_y_z_y_yyyzz[k] * ab_x + g_y_z_y_xyyyzz[k];

                g_y_z_xy_yyzzz[k] = -g_y_z_y_yyzzz[k] * ab_x + g_y_z_y_xyyzzz[k];

                g_y_z_xy_yzzzz[k] = -g_y_z_y_yzzzz[k] * ab_x + g_y_z_y_xyzzzz[k];

                g_y_z_xy_zzzzz[k] = -g_y_z_y_zzzzz[k] * ab_x + g_y_z_y_xzzzzz[k];
            }

            /// Set up 672-693 components of targeted buffer : cbuffer.data(

            auto g_y_z_xz_xxxxx = cbuffer.data(dh_geom_11_off + 672 * ccomps * dcomps);

            auto g_y_z_xz_xxxxy = cbuffer.data(dh_geom_11_off + 673 * ccomps * dcomps);

            auto g_y_z_xz_xxxxz = cbuffer.data(dh_geom_11_off + 674 * ccomps * dcomps);

            auto g_y_z_xz_xxxyy = cbuffer.data(dh_geom_11_off + 675 * ccomps * dcomps);

            auto g_y_z_xz_xxxyz = cbuffer.data(dh_geom_11_off + 676 * ccomps * dcomps);

            auto g_y_z_xz_xxxzz = cbuffer.data(dh_geom_11_off + 677 * ccomps * dcomps);

            auto g_y_z_xz_xxyyy = cbuffer.data(dh_geom_11_off + 678 * ccomps * dcomps);

            auto g_y_z_xz_xxyyz = cbuffer.data(dh_geom_11_off + 679 * ccomps * dcomps);

            auto g_y_z_xz_xxyzz = cbuffer.data(dh_geom_11_off + 680 * ccomps * dcomps);

            auto g_y_z_xz_xxzzz = cbuffer.data(dh_geom_11_off + 681 * ccomps * dcomps);

            auto g_y_z_xz_xyyyy = cbuffer.data(dh_geom_11_off + 682 * ccomps * dcomps);

            auto g_y_z_xz_xyyyz = cbuffer.data(dh_geom_11_off + 683 * ccomps * dcomps);

            auto g_y_z_xz_xyyzz = cbuffer.data(dh_geom_11_off + 684 * ccomps * dcomps);

            auto g_y_z_xz_xyzzz = cbuffer.data(dh_geom_11_off + 685 * ccomps * dcomps);

            auto g_y_z_xz_xzzzz = cbuffer.data(dh_geom_11_off + 686 * ccomps * dcomps);

            auto g_y_z_xz_yyyyy = cbuffer.data(dh_geom_11_off + 687 * ccomps * dcomps);

            auto g_y_z_xz_yyyyz = cbuffer.data(dh_geom_11_off + 688 * ccomps * dcomps);

            auto g_y_z_xz_yyyzz = cbuffer.data(dh_geom_11_off + 689 * ccomps * dcomps);

            auto g_y_z_xz_yyzzz = cbuffer.data(dh_geom_11_off + 690 * ccomps * dcomps);

            auto g_y_z_xz_yzzzz = cbuffer.data(dh_geom_11_off + 691 * ccomps * dcomps);

            auto g_y_z_xz_zzzzz = cbuffer.data(dh_geom_11_off + 692 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_z_xz_xxxxx, g_y_z_xz_xxxxy, g_y_z_xz_xxxxz, g_y_z_xz_xxxyy, g_y_z_xz_xxxyz, g_y_z_xz_xxxzz, g_y_z_xz_xxyyy, g_y_z_xz_xxyyz, g_y_z_xz_xxyzz, g_y_z_xz_xxzzz, g_y_z_xz_xyyyy, g_y_z_xz_xyyyz, g_y_z_xz_xyyzz, g_y_z_xz_xyzzz, g_y_z_xz_xzzzz, g_y_z_xz_yyyyy, g_y_z_xz_yyyyz, g_y_z_xz_yyyzz, g_y_z_xz_yyzzz, g_y_z_xz_yzzzz, g_y_z_xz_zzzzz, g_y_z_z_xxxxx, g_y_z_z_xxxxxx, g_y_z_z_xxxxxy, g_y_z_z_xxxxxz, g_y_z_z_xxxxy, g_y_z_z_xxxxyy, g_y_z_z_xxxxyz, g_y_z_z_xxxxz, g_y_z_z_xxxxzz, g_y_z_z_xxxyy, g_y_z_z_xxxyyy, g_y_z_z_xxxyyz, g_y_z_z_xxxyz, g_y_z_z_xxxyzz, g_y_z_z_xxxzz, g_y_z_z_xxxzzz, g_y_z_z_xxyyy, g_y_z_z_xxyyyy, g_y_z_z_xxyyyz, g_y_z_z_xxyyz, g_y_z_z_xxyyzz, g_y_z_z_xxyzz, g_y_z_z_xxyzzz, g_y_z_z_xxzzz, g_y_z_z_xxzzzz, g_y_z_z_xyyyy, g_y_z_z_xyyyyy, g_y_z_z_xyyyyz, g_y_z_z_xyyyz, g_y_z_z_xyyyzz, g_y_z_z_xyyzz, g_y_z_z_xyyzzz, g_y_z_z_xyzzz, g_y_z_z_xyzzzz, g_y_z_z_xzzzz, g_y_z_z_xzzzzz, g_y_z_z_yyyyy, g_y_z_z_yyyyz, g_y_z_z_yyyzz, g_y_z_z_yyzzz, g_y_z_z_yzzzz, g_y_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_xz_xxxxx[k] = -g_y_z_z_xxxxx[k] * ab_x + g_y_z_z_xxxxxx[k];

                g_y_z_xz_xxxxy[k] = -g_y_z_z_xxxxy[k] * ab_x + g_y_z_z_xxxxxy[k];

                g_y_z_xz_xxxxz[k] = -g_y_z_z_xxxxz[k] * ab_x + g_y_z_z_xxxxxz[k];

                g_y_z_xz_xxxyy[k] = -g_y_z_z_xxxyy[k] * ab_x + g_y_z_z_xxxxyy[k];

                g_y_z_xz_xxxyz[k] = -g_y_z_z_xxxyz[k] * ab_x + g_y_z_z_xxxxyz[k];

                g_y_z_xz_xxxzz[k] = -g_y_z_z_xxxzz[k] * ab_x + g_y_z_z_xxxxzz[k];

                g_y_z_xz_xxyyy[k] = -g_y_z_z_xxyyy[k] * ab_x + g_y_z_z_xxxyyy[k];

                g_y_z_xz_xxyyz[k] = -g_y_z_z_xxyyz[k] * ab_x + g_y_z_z_xxxyyz[k];

                g_y_z_xz_xxyzz[k] = -g_y_z_z_xxyzz[k] * ab_x + g_y_z_z_xxxyzz[k];

                g_y_z_xz_xxzzz[k] = -g_y_z_z_xxzzz[k] * ab_x + g_y_z_z_xxxzzz[k];

                g_y_z_xz_xyyyy[k] = -g_y_z_z_xyyyy[k] * ab_x + g_y_z_z_xxyyyy[k];

                g_y_z_xz_xyyyz[k] = -g_y_z_z_xyyyz[k] * ab_x + g_y_z_z_xxyyyz[k];

                g_y_z_xz_xyyzz[k] = -g_y_z_z_xyyzz[k] * ab_x + g_y_z_z_xxyyzz[k];

                g_y_z_xz_xyzzz[k] = -g_y_z_z_xyzzz[k] * ab_x + g_y_z_z_xxyzzz[k];

                g_y_z_xz_xzzzz[k] = -g_y_z_z_xzzzz[k] * ab_x + g_y_z_z_xxzzzz[k];

                g_y_z_xz_yyyyy[k] = -g_y_z_z_yyyyy[k] * ab_x + g_y_z_z_xyyyyy[k];

                g_y_z_xz_yyyyz[k] = -g_y_z_z_yyyyz[k] * ab_x + g_y_z_z_xyyyyz[k];

                g_y_z_xz_yyyzz[k] = -g_y_z_z_yyyzz[k] * ab_x + g_y_z_z_xyyyzz[k];

                g_y_z_xz_yyzzz[k] = -g_y_z_z_yyzzz[k] * ab_x + g_y_z_z_xyyzzz[k];

                g_y_z_xz_yzzzz[k] = -g_y_z_z_yzzzz[k] * ab_x + g_y_z_z_xyzzzz[k];

                g_y_z_xz_zzzzz[k] = -g_y_z_z_zzzzz[k] * ab_x + g_y_z_z_xzzzzz[k];
            }

            /// Set up 693-714 components of targeted buffer : cbuffer.data(

            auto g_y_z_yy_xxxxx = cbuffer.data(dh_geom_11_off + 693 * ccomps * dcomps);

            auto g_y_z_yy_xxxxy = cbuffer.data(dh_geom_11_off + 694 * ccomps * dcomps);

            auto g_y_z_yy_xxxxz = cbuffer.data(dh_geom_11_off + 695 * ccomps * dcomps);

            auto g_y_z_yy_xxxyy = cbuffer.data(dh_geom_11_off + 696 * ccomps * dcomps);

            auto g_y_z_yy_xxxyz = cbuffer.data(dh_geom_11_off + 697 * ccomps * dcomps);

            auto g_y_z_yy_xxxzz = cbuffer.data(dh_geom_11_off + 698 * ccomps * dcomps);

            auto g_y_z_yy_xxyyy = cbuffer.data(dh_geom_11_off + 699 * ccomps * dcomps);

            auto g_y_z_yy_xxyyz = cbuffer.data(dh_geom_11_off + 700 * ccomps * dcomps);

            auto g_y_z_yy_xxyzz = cbuffer.data(dh_geom_11_off + 701 * ccomps * dcomps);

            auto g_y_z_yy_xxzzz = cbuffer.data(dh_geom_11_off + 702 * ccomps * dcomps);

            auto g_y_z_yy_xyyyy = cbuffer.data(dh_geom_11_off + 703 * ccomps * dcomps);

            auto g_y_z_yy_xyyyz = cbuffer.data(dh_geom_11_off + 704 * ccomps * dcomps);

            auto g_y_z_yy_xyyzz = cbuffer.data(dh_geom_11_off + 705 * ccomps * dcomps);

            auto g_y_z_yy_xyzzz = cbuffer.data(dh_geom_11_off + 706 * ccomps * dcomps);

            auto g_y_z_yy_xzzzz = cbuffer.data(dh_geom_11_off + 707 * ccomps * dcomps);

            auto g_y_z_yy_yyyyy = cbuffer.data(dh_geom_11_off + 708 * ccomps * dcomps);

            auto g_y_z_yy_yyyyz = cbuffer.data(dh_geom_11_off + 709 * ccomps * dcomps);

            auto g_y_z_yy_yyyzz = cbuffer.data(dh_geom_11_off + 710 * ccomps * dcomps);

            auto g_y_z_yy_yyzzz = cbuffer.data(dh_geom_11_off + 711 * ccomps * dcomps);

            auto g_y_z_yy_yzzzz = cbuffer.data(dh_geom_11_off + 712 * ccomps * dcomps);

            auto g_y_z_yy_zzzzz = cbuffer.data(dh_geom_11_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_y_xxxxx, g_0_z_y_xxxxy, g_0_z_y_xxxxz, g_0_z_y_xxxyy, g_0_z_y_xxxyz, g_0_z_y_xxxzz, g_0_z_y_xxyyy, g_0_z_y_xxyyz, g_0_z_y_xxyzz, g_0_z_y_xxzzz, g_0_z_y_xyyyy, g_0_z_y_xyyyz, g_0_z_y_xyyzz, g_0_z_y_xyzzz, g_0_z_y_xzzzz, g_0_z_y_yyyyy, g_0_z_y_yyyyz, g_0_z_y_yyyzz, g_0_z_y_yyzzz, g_0_z_y_yzzzz, g_0_z_y_zzzzz, g_y_z_y_xxxxx, g_y_z_y_xxxxxy, g_y_z_y_xxxxy, g_y_z_y_xxxxyy, g_y_z_y_xxxxyz, g_y_z_y_xxxxz, g_y_z_y_xxxyy, g_y_z_y_xxxyyy, g_y_z_y_xxxyyz, g_y_z_y_xxxyz, g_y_z_y_xxxyzz, g_y_z_y_xxxzz, g_y_z_y_xxyyy, g_y_z_y_xxyyyy, g_y_z_y_xxyyyz, g_y_z_y_xxyyz, g_y_z_y_xxyyzz, g_y_z_y_xxyzz, g_y_z_y_xxyzzz, g_y_z_y_xxzzz, g_y_z_y_xyyyy, g_y_z_y_xyyyyy, g_y_z_y_xyyyyz, g_y_z_y_xyyyz, g_y_z_y_xyyyzz, g_y_z_y_xyyzz, g_y_z_y_xyyzzz, g_y_z_y_xyzzz, g_y_z_y_xyzzzz, g_y_z_y_xzzzz, g_y_z_y_yyyyy, g_y_z_y_yyyyyy, g_y_z_y_yyyyyz, g_y_z_y_yyyyz, g_y_z_y_yyyyzz, g_y_z_y_yyyzz, g_y_z_y_yyyzzz, g_y_z_y_yyzzz, g_y_z_y_yyzzzz, g_y_z_y_yzzzz, g_y_z_y_yzzzzz, g_y_z_y_zzzzz, g_y_z_yy_xxxxx, g_y_z_yy_xxxxy, g_y_z_yy_xxxxz, g_y_z_yy_xxxyy, g_y_z_yy_xxxyz, g_y_z_yy_xxxzz, g_y_z_yy_xxyyy, g_y_z_yy_xxyyz, g_y_z_yy_xxyzz, g_y_z_yy_xxzzz, g_y_z_yy_xyyyy, g_y_z_yy_xyyyz, g_y_z_yy_xyyzz, g_y_z_yy_xyzzz, g_y_z_yy_xzzzz, g_y_z_yy_yyyyy, g_y_z_yy_yyyyz, g_y_z_yy_yyyzz, g_y_z_yy_yyzzz, g_y_z_yy_yzzzz, g_y_z_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yy_xxxxx[k] = -g_0_z_y_xxxxx[k] - g_y_z_y_xxxxx[k] * ab_y + g_y_z_y_xxxxxy[k];

                g_y_z_yy_xxxxy[k] = -g_0_z_y_xxxxy[k] - g_y_z_y_xxxxy[k] * ab_y + g_y_z_y_xxxxyy[k];

                g_y_z_yy_xxxxz[k] = -g_0_z_y_xxxxz[k] - g_y_z_y_xxxxz[k] * ab_y + g_y_z_y_xxxxyz[k];

                g_y_z_yy_xxxyy[k] = -g_0_z_y_xxxyy[k] - g_y_z_y_xxxyy[k] * ab_y + g_y_z_y_xxxyyy[k];

                g_y_z_yy_xxxyz[k] = -g_0_z_y_xxxyz[k] - g_y_z_y_xxxyz[k] * ab_y + g_y_z_y_xxxyyz[k];

                g_y_z_yy_xxxzz[k] = -g_0_z_y_xxxzz[k] - g_y_z_y_xxxzz[k] * ab_y + g_y_z_y_xxxyzz[k];

                g_y_z_yy_xxyyy[k] = -g_0_z_y_xxyyy[k] - g_y_z_y_xxyyy[k] * ab_y + g_y_z_y_xxyyyy[k];

                g_y_z_yy_xxyyz[k] = -g_0_z_y_xxyyz[k] - g_y_z_y_xxyyz[k] * ab_y + g_y_z_y_xxyyyz[k];

                g_y_z_yy_xxyzz[k] = -g_0_z_y_xxyzz[k] - g_y_z_y_xxyzz[k] * ab_y + g_y_z_y_xxyyzz[k];

                g_y_z_yy_xxzzz[k] = -g_0_z_y_xxzzz[k] - g_y_z_y_xxzzz[k] * ab_y + g_y_z_y_xxyzzz[k];

                g_y_z_yy_xyyyy[k] = -g_0_z_y_xyyyy[k] - g_y_z_y_xyyyy[k] * ab_y + g_y_z_y_xyyyyy[k];

                g_y_z_yy_xyyyz[k] = -g_0_z_y_xyyyz[k] - g_y_z_y_xyyyz[k] * ab_y + g_y_z_y_xyyyyz[k];

                g_y_z_yy_xyyzz[k] = -g_0_z_y_xyyzz[k] - g_y_z_y_xyyzz[k] * ab_y + g_y_z_y_xyyyzz[k];

                g_y_z_yy_xyzzz[k] = -g_0_z_y_xyzzz[k] - g_y_z_y_xyzzz[k] * ab_y + g_y_z_y_xyyzzz[k];

                g_y_z_yy_xzzzz[k] = -g_0_z_y_xzzzz[k] - g_y_z_y_xzzzz[k] * ab_y + g_y_z_y_xyzzzz[k];

                g_y_z_yy_yyyyy[k] = -g_0_z_y_yyyyy[k] - g_y_z_y_yyyyy[k] * ab_y + g_y_z_y_yyyyyy[k];

                g_y_z_yy_yyyyz[k] = -g_0_z_y_yyyyz[k] - g_y_z_y_yyyyz[k] * ab_y + g_y_z_y_yyyyyz[k];

                g_y_z_yy_yyyzz[k] = -g_0_z_y_yyyzz[k] - g_y_z_y_yyyzz[k] * ab_y + g_y_z_y_yyyyzz[k];

                g_y_z_yy_yyzzz[k] = -g_0_z_y_yyzzz[k] - g_y_z_y_yyzzz[k] * ab_y + g_y_z_y_yyyzzz[k];

                g_y_z_yy_yzzzz[k] = -g_0_z_y_yzzzz[k] - g_y_z_y_yzzzz[k] * ab_y + g_y_z_y_yyzzzz[k];

                g_y_z_yy_zzzzz[k] = -g_0_z_y_zzzzz[k] - g_y_z_y_zzzzz[k] * ab_y + g_y_z_y_yzzzzz[k];
            }

            /// Set up 714-735 components of targeted buffer : cbuffer.data(

            auto g_y_z_yz_xxxxx = cbuffer.data(dh_geom_11_off + 714 * ccomps * dcomps);

            auto g_y_z_yz_xxxxy = cbuffer.data(dh_geom_11_off + 715 * ccomps * dcomps);

            auto g_y_z_yz_xxxxz = cbuffer.data(dh_geom_11_off + 716 * ccomps * dcomps);

            auto g_y_z_yz_xxxyy = cbuffer.data(dh_geom_11_off + 717 * ccomps * dcomps);

            auto g_y_z_yz_xxxyz = cbuffer.data(dh_geom_11_off + 718 * ccomps * dcomps);

            auto g_y_z_yz_xxxzz = cbuffer.data(dh_geom_11_off + 719 * ccomps * dcomps);

            auto g_y_z_yz_xxyyy = cbuffer.data(dh_geom_11_off + 720 * ccomps * dcomps);

            auto g_y_z_yz_xxyyz = cbuffer.data(dh_geom_11_off + 721 * ccomps * dcomps);

            auto g_y_z_yz_xxyzz = cbuffer.data(dh_geom_11_off + 722 * ccomps * dcomps);

            auto g_y_z_yz_xxzzz = cbuffer.data(dh_geom_11_off + 723 * ccomps * dcomps);

            auto g_y_z_yz_xyyyy = cbuffer.data(dh_geom_11_off + 724 * ccomps * dcomps);

            auto g_y_z_yz_xyyyz = cbuffer.data(dh_geom_11_off + 725 * ccomps * dcomps);

            auto g_y_z_yz_xyyzz = cbuffer.data(dh_geom_11_off + 726 * ccomps * dcomps);

            auto g_y_z_yz_xyzzz = cbuffer.data(dh_geom_11_off + 727 * ccomps * dcomps);

            auto g_y_z_yz_xzzzz = cbuffer.data(dh_geom_11_off + 728 * ccomps * dcomps);

            auto g_y_z_yz_yyyyy = cbuffer.data(dh_geom_11_off + 729 * ccomps * dcomps);

            auto g_y_z_yz_yyyyz = cbuffer.data(dh_geom_11_off + 730 * ccomps * dcomps);

            auto g_y_z_yz_yyyzz = cbuffer.data(dh_geom_11_off + 731 * ccomps * dcomps);

            auto g_y_z_yz_yyzzz = cbuffer.data(dh_geom_11_off + 732 * ccomps * dcomps);

            auto g_y_z_yz_yzzzz = cbuffer.data(dh_geom_11_off + 733 * ccomps * dcomps);

            auto g_y_z_yz_zzzzz = cbuffer.data(dh_geom_11_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxxx, g_0_z_z_xxxxy, g_0_z_z_xxxxz, g_0_z_z_xxxyy, g_0_z_z_xxxyz, g_0_z_z_xxxzz, g_0_z_z_xxyyy, g_0_z_z_xxyyz, g_0_z_z_xxyzz, g_0_z_z_xxzzz, g_0_z_z_xyyyy, g_0_z_z_xyyyz, g_0_z_z_xyyzz, g_0_z_z_xyzzz, g_0_z_z_xzzzz, g_0_z_z_yyyyy, g_0_z_z_yyyyz, g_0_z_z_yyyzz, g_0_z_z_yyzzz, g_0_z_z_yzzzz, g_0_z_z_zzzzz, g_y_z_yz_xxxxx, g_y_z_yz_xxxxy, g_y_z_yz_xxxxz, g_y_z_yz_xxxyy, g_y_z_yz_xxxyz, g_y_z_yz_xxxzz, g_y_z_yz_xxyyy, g_y_z_yz_xxyyz, g_y_z_yz_xxyzz, g_y_z_yz_xxzzz, g_y_z_yz_xyyyy, g_y_z_yz_xyyyz, g_y_z_yz_xyyzz, g_y_z_yz_xyzzz, g_y_z_yz_xzzzz, g_y_z_yz_yyyyy, g_y_z_yz_yyyyz, g_y_z_yz_yyyzz, g_y_z_yz_yyzzz, g_y_z_yz_yzzzz, g_y_z_yz_zzzzz, g_y_z_z_xxxxx, g_y_z_z_xxxxxy, g_y_z_z_xxxxy, g_y_z_z_xxxxyy, g_y_z_z_xxxxyz, g_y_z_z_xxxxz, g_y_z_z_xxxyy, g_y_z_z_xxxyyy, g_y_z_z_xxxyyz, g_y_z_z_xxxyz, g_y_z_z_xxxyzz, g_y_z_z_xxxzz, g_y_z_z_xxyyy, g_y_z_z_xxyyyy, g_y_z_z_xxyyyz, g_y_z_z_xxyyz, g_y_z_z_xxyyzz, g_y_z_z_xxyzz, g_y_z_z_xxyzzz, g_y_z_z_xxzzz, g_y_z_z_xyyyy, g_y_z_z_xyyyyy, g_y_z_z_xyyyyz, g_y_z_z_xyyyz, g_y_z_z_xyyyzz, g_y_z_z_xyyzz, g_y_z_z_xyyzzz, g_y_z_z_xyzzz, g_y_z_z_xyzzzz, g_y_z_z_xzzzz, g_y_z_z_yyyyy, g_y_z_z_yyyyyy, g_y_z_z_yyyyyz, g_y_z_z_yyyyz, g_y_z_z_yyyyzz, g_y_z_z_yyyzz, g_y_z_z_yyyzzz, g_y_z_z_yyzzz, g_y_z_z_yyzzzz, g_y_z_z_yzzzz, g_y_z_z_yzzzzz, g_y_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_yz_xxxxx[k] = -g_0_z_z_xxxxx[k] - g_y_z_z_xxxxx[k] * ab_y + g_y_z_z_xxxxxy[k];

                g_y_z_yz_xxxxy[k] = -g_0_z_z_xxxxy[k] - g_y_z_z_xxxxy[k] * ab_y + g_y_z_z_xxxxyy[k];

                g_y_z_yz_xxxxz[k] = -g_0_z_z_xxxxz[k] - g_y_z_z_xxxxz[k] * ab_y + g_y_z_z_xxxxyz[k];

                g_y_z_yz_xxxyy[k] = -g_0_z_z_xxxyy[k] - g_y_z_z_xxxyy[k] * ab_y + g_y_z_z_xxxyyy[k];

                g_y_z_yz_xxxyz[k] = -g_0_z_z_xxxyz[k] - g_y_z_z_xxxyz[k] * ab_y + g_y_z_z_xxxyyz[k];

                g_y_z_yz_xxxzz[k] = -g_0_z_z_xxxzz[k] - g_y_z_z_xxxzz[k] * ab_y + g_y_z_z_xxxyzz[k];

                g_y_z_yz_xxyyy[k] = -g_0_z_z_xxyyy[k] - g_y_z_z_xxyyy[k] * ab_y + g_y_z_z_xxyyyy[k];

                g_y_z_yz_xxyyz[k] = -g_0_z_z_xxyyz[k] - g_y_z_z_xxyyz[k] * ab_y + g_y_z_z_xxyyyz[k];

                g_y_z_yz_xxyzz[k] = -g_0_z_z_xxyzz[k] - g_y_z_z_xxyzz[k] * ab_y + g_y_z_z_xxyyzz[k];

                g_y_z_yz_xxzzz[k] = -g_0_z_z_xxzzz[k] - g_y_z_z_xxzzz[k] * ab_y + g_y_z_z_xxyzzz[k];

                g_y_z_yz_xyyyy[k] = -g_0_z_z_xyyyy[k] - g_y_z_z_xyyyy[k] * ab_y + g_y_z_z_xyyyyy[k];

                g_y_z_yz_xyyyz[k] = -g_0_z_z_xyyyz[k] - g_y_z_z_xyyyz[k] * ab_y + g_y_z_z_xyyyyz[k];

                g_y_z_yz_xyyzz[k] = -g_0_z_z_xyyzz[k] - g_y_z_z_xyyzz[k] * ab_y + g_y_z_z_xyyyzz[k];

                g_y_z_yz_xyzzz[k] = -g_0_z_z_xyzzz[k] - g_y_z_z_xyzzz[k] * ab_y + g_y_z_z_xyyzzz[k];

                g_y_z_yz_xzzzz[k] = -g_0_z_z_xzzzz[k] - g_y_z_z_xzzzz[k] * ab_y + g_y_z_z_xyzzzz[k];

                g_y_z_yz_yyyyy[k] = -g_0_z_z_yyyyy[k] - g_y_z_z_yyyyy[k] * ab_y + g_y_z_z_yyyyyy[k];

                g_y_z_yz_yyyyz[k] = -g_0_z_z_yyyyz[k] - g_y_z_z_yyyyz[k] * ab_y + g_y_z_z_yyyyyz[k];

                g_y_z_yz_yyyzz[k] = -g_0_z_z_yyyzz[k] - g_y_z_z_yyyzz[k] * ab_y + g_y_z_z_yyyyzz[k];

                g_y_z_yz_yyzzz[k] = -g_0_z_z_yyzzz[k] - g_y_z_z_yyzzz[k] * ab_y + g_y_z_z_yyyzzz[k];

                g_y_z_yz_yzzzz[k] = -g_0_z_z_yzzzz[k] - g_y_z_z_yzzzz[k] * ab_y + g_y_z_z_yyzzzz[k];

                g_y_z_yz_zzzzz[k] = -g_0_z_z_zzzzz[k] - g_y_z_z_zzzzz[k] * ab_y + g_y_z_z_yzzzzz[k];
            }

            /// Set up 735-756 components of targeted buffer : cbuffer.data(

            auto g_y_z_zz_xxxxx = cbuffer.data(dh_geom_11_off + 735 * ccomps * dcomps);

            auto g_y_z_zz_xxxxy = cbuffer.data(dh_geom_11_off + 736 * ccomps * dcomps);

            auto g_y_z_zz_xxxxz = cbuffer.data(dh_geom_11_off + 737 * ccomps * dcomps);

            auto g_y_z_zz_xxxyy = cbuffer.data(dh_geom_11_off + 738 * ccomps * dcomps);

            auto g_y_z_zz_xxxyz = cbuffer.data(dh_geom_11_off + 739 * ccomps * dcomps);

            auto g_y_z_zz_xxxzz = cbuffer.data(dh_geom_11_off + 740 * ccomps * dcomps);

            auto g_y_z_zz_xxyyy = cbuffer.data(dh_geom_11_off + 741 * ccomps * dcomps);

            auto g_y_z_zz_xxyyz = cbuffer.data(dh_geom_11_off + 742 * ccomps * dcomps);

            auto g_y_z_zz_xxyzz = cbuffer.data(dh_geom_11_off + 743 * ccomps * dcomps);

            auto g_y_z_zz_xxzzz = cbuffer.data(dh_geom_11_off + 744 * ccomps * dcomps);

            auto g_y_z_zz_xyyyy = cbuffer.data(dh_geom_11_off + 745 * ccomps * dcomps);

            auto g_y_z_zz_xyyyz = cbuffer.data(dh_geom_11_off + 746 * ccomps * dcomps);

            auto g_y_z_zz_xyyzz = cbuffer.data(dh_geom_11_off + 747 * ccomps * dcomps);

            auto g_y_z_zz_xyzzz = cbuffer.data(dh_geom_11_off + 748 * ccomps * dcomps);

            auto g_y_z_zz_xzzzz = cbuffer.data(dh_geom_11_off + 749 * ccomps * dcomps);

            auto g_y_z_zz_yyyyy = cbuffer.data(dh_geom_11_off + 750 * ccomps * dcomps);

            auto g_y_z_zz_yyyyz = cbuffer.data(dh_geom_11_off + 751 * ccomps * dcomps);

            auto g_y_z_zz_yyyzz = cbuffer.data(dh_geom_11_off + 752 * ccomps * dcomps);

            auto g_y_z_zz_yyzzz = cbuffer.data(dh_geom_11_off + 753 * ccomps * dcomps);

            auto g_y_z_zz_yzzzz = cbuffer.data(dh_geom_11_off + 754 * ccomps * dcomps);

            auto g_y_z_zz_zzzzz = cbuffer.data(dh_geom_11_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_z_xxxxx, g_y_0_z_xxxxy, g_y_0_z_xxxxz, g_y_0_z_xxxyy, g_y_0_z_xxxyz, g_y_0_z_xxxzz, g_y_0_z_xxyyy, g_y_0_z_xxyyz, g_y_0_z_xxyzz, g_y_0_z_xxzzz, g_y_0_z_xyyyy, g_y_0_z_xyyyz, g_y_0_z_xyyzz, g_y_0_z_xyzzz, g_y_0_z_xzzzz, g_y_0_z_yyyyy, g_y_0_z_yyyyz, g_y_0_z_yyyzz, g_y_0_z_yyzzz, g_y_0_z_yzzzz, g_y_0_z_zzzzz, g_y_z_z_xxxxx, g_y_z_z_xxxxxz, g_y_z_z_xxxxy, g_y_z_z_xxxxyz, g_y_z_z_xxxxz, g_y_z_z_xxxxzz, g_y_z_z_xxxyy, g_y_z_z_xxxyyz, g_y_z_z_xxxyz, g_y_z_z_xxxyzz, g_y_z_z_xxxzz, g_y_z_z_xxxzzz, g_y_z_z_xxyyy, g_y_z_z_xxyyyz, g_y_z_z_xxyyz, g_y_z_z_xxyyzz, g_y_z_z_xxyzz, g_y_z_z_xxyzzz, g_y_z_z_xxzzz, g_y_z_z_xxzzzz, g_y_z_z_xyyyy, g_y_z_z_xyyyyz, g_y_z_z_xyyyz, g_y_z_z_xyyyzz, g_y_z_z_xyyzz, g_y_z_z_xyyzzz, g_y_z_z_xyzzz, g_y_z_z_xyzzzz, g_y_z_z_xzzzz, g_y_z_z_xzzzzz, g_y_z_z_yyyyy, g_y_z_z_yyyyyz, g_y_z_z_yyyyz, g_y_z_z_yyyyzz, g_y_z_z_yyyzz, g_y_z_z_yyyzzz, g_y_z_z_yyzzz, g_y_z_z_yyzzzz, g_y_z_z_yzzzz, g_y_z_z_yzzzzz, g_y_z_z_zzzzz, g_y_z_z_zzzzzz, g_y_z_zz_xxxxx, g_y_z_zz_xxxxy, g_y_z_zz_xxxxz, g_y_z_zz_xxxyy, g_y_z_zz_xxxyz, g_y_z_zz_xxxzz, g_y_z_zz_xxyyy, g_y_z_zz_xxyyz, g_y_z_zz_xxyzz, g_y_z_zz_xxzzz, g_y_z_zz_xyyyy, g_y_z_zz_xyyyz, g_y_z_zz_xyyzz, g_y_z_zz_xyzzz, g_y_z_zz_xzzzz, g_y_z_zz_yyyyy, g_y_z_zz_yyyyz, g_y_z_zz_yyyzz, g_y_z_zz_yyzzz, g_y_z_zz_yzzzz, g_y_z_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_z_zz_xxxxx[k] = g_y_0_z_xxxxx[k] - g_y_z_z_xxxxx[k] * ab_z + g_y_z_z_xxxxxz[k];

                g_y_z_zz_xxxxy[k] = g_y_0_z_xxxxy[k] - g_y_z_z_xxxxy[k] * ab_z + g_y_z_z_xxxxyz[k];

                g_y_z_zz_xxxxz[k] = g_y_0_z_xxxxz[k] - g_y_z_z_xxxxz[k] * ab_z + g_y_z_z_xxxxzz[k];

                g_y_z_zz_xxxyy[k] = g_y_0_z_xxxyy[k] - g_y_z_z_xxxyy[k] * ab_z + g_y_z_z_xxxyyz[k];

                g_y_z_zz_xxxyz[k] = g_y_0_z_xxxyz[k] - g_y_z_z_xxxyz[k] * ab_z + g_y_z_z_xxxyzz[k];

                g_y_z_zz_xxxzz[k] = g_y_0_z_xxxzz[k] - g_y_z_z_xxxzz[k] * ab_z + g_y_z_z_xxxzzz[k];

                g_y_z_zz_xxyyy[k] = g_y_0_z_xxyyy[k] - g_y_z_z_xxyyy[k] * ab_z + g_y_z_z_xxyyyz[k];

                g_y_z_zz_xxyyz[k] = g_y_0_z_xxyyz[k] - g_y_z_z_xxyyz[k] * ab_z + g_y_z_z_xxyyzz[k];

                g_y_z_zz_xxyzz[k] = g_y_0_z_xxyzz[k] - g_y_z_z_xxyzz[k] * ab_z + g_y_z_z_xxyzzz[k];

                g_y_z_zz_xxzzz[k] = g_y_0_z_xxzzz[k] - g_y_z_z_xxzzz[k] * ab_z + g_y_z_z_xxzzzz[k];

                g_y_z_zz_xyyyy[k] = g_y_0_z_xyyyy[k] - g_y_z_z_xyyyy[k] * ab_z + g_y_z_z_xyyyyz[k];

                g_y_z_zz_xyyyz[k] = g_y_0_z_xyyyz[k] - g_y_z_z_xyyyz[k] * ab_z + g_y_z_z_xyyyzz[k];

                g_y_z_zz_xyyzz[k] = g_y_0_z_xyyzz[k] - g_y_z_z_xyyzz[k] * ab_z + g_y_z_z_xyyzzz[k];

                g_y_z_zz_xyzzz[k] = g_y_0_z_xyzzz[k] - g_y_z_z_xyzzz[k] * ab_z + g_y_z_z_xyzzzz[k];

                g_y_z_zz_xzzzz[k] = g_y_0_z_xzzzz[k] - g_y_z_z_xzzzz[k] * ab_z + g_y_z_z_xzzzzz[k];

                g_y_z_zz_yyyyy[k] = g_y_0_z_yyyyy[k] - g_y_z_z_yyyyy[k] * ab_z + g_y_z_z_yyyyyz[k];

                g_y_z_zz_yyyyz[k] = g_y_0_z_yyyyz[k] - g_y_z_z_yyyyz[k] * ab_z + g_y_z_z_yyyyzz[k];

                g_y_z_zz_yyyzz[k] = g_y_0_z_yyyzz[k] - g_y_z_z_yyyzz[k] * ab_z + g_y_z_z_yyyzzz[k];

                g_y_z_zz_yyzzz[k] = g_y_0_z_yyzzz[k] - g_y_z_z_yyzzz[k] * ab_z + g_y_z_z_yyzzzz[k];

                g_y_z_zz_yzzzz[k] = g_y_0_z_yzzzz[k] - g_y_z_z_yzzzz[k] * ab_z + g_y_z_z_yzzzzz[k];

                g_y_z_zz_zzzzz[k] = g_y_0_z_zzzzz[k] - g_y_z_z_zzzzz[k] * ab_z + g_y_z_z_zzzzzz[k];
            }

            /// Set up 756-777 components of targeted buffer : cbuffer.data(

            auto g_z_x_xx_xxxxx = cbuffer.data(dh_geom_11_off + 756 * ccomps * dcomps);

            auto g_z_x_xx_xxxxy = cbuffer.data(dh_geom_11_off + 757 * ccomps * dcomps);

            auto g_z_x_xx_xxxxz = cbuffer.data(dh_geom_11_off + 758 * ccomps * dcomps);

            auto g_z_x_xx_xxxyy = cbuffer.data(dh_geom_11_off + 759 * ccomps * dcomps);

            auto g_z_x_xx_xxxyz = cbuffer.data(dh_geom_11_off + 760 * ccomps * dcomps);

            auto g_z_x_xx_xxxzz = cbuffer.data(dh_geom_11_off + 761 * ccomps * dcomps);

            auto g_z_x_xx_xxyyy = cbuffer.data(dh_geom_11_off + 762 * ccomps * dcomps);

            auto g_z_x_xx_xxyyz = cbuffer.data(dh_geom_11_off + 763 * ccomps * dcomps);

            auto g_z_x_xx_xxyzz = cbuffer.data(dh_geom_11_off + 764 * ccomps * dcomps);

            auto g_z_x_xx_xxzzz = cbuffer.data(dh_geom_11_off + 765 * ccomps * dcomps);

            auto g_z_x_xx_xyyyy = cbuffer.data(dh_geom_11_off + 766 * ccomps * dcomps);

            auto g_z_x_xx_xyyyz = cbuffer.data(dh_geom_11_off + 767 * ccomps * dcomps);

            auto g_z_x_xx_xyyzz = cbuffer.data(dh_geom_11_off + 768 * ccomps * dcomps);

            auto g_z_x_xx_xyzzz = cbuffer.data(dh_geom_11_off + 769 * ccomps * dcomps);

            auto g_z_x_xx_xzzzz = cbuffer.data(dh_geom_11_off + 770 * ccomps * dcomps);

            auto g_z_x_xx_yyyyy = cbuffer.data(dh_geom_11_off + 771 * ccomps * dcomps);

            auto g_z_x_xx_yyyyz = cbuffer.data(dh_geom_11_off + 772 * ccomps * dcomps);

            auto g_z_x_xx_yyyzz = cbuffer.data(dh_geom_11_off + 773 * ccomps * dcomps);

            auto g_z_x_xx_yyzzz = cbuffer.data(dh_geom_11_off + 774 * ccomps * dcomps);

            auto g_z_x_xx_yzzzz = cbuffer.data(dh_geom_11_off + 775 * ccomps * dcomps);

            auto g_z_x_xx_zzzzz = cbuffer.data(dh_geom_11_off + 776 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_x_xxxxx, g_z_0_x_xxxxy, g_z_0_x_xxxxz, g_z_0_x_xxxyy, g_z_0_x_xxxyz, g_z_0_x_xxxzz, g_z_0_x_xxyyy, g_z_0_x_xxyyz, g_z_0_x_xxyzz, g_z_0_x_xxzzz, g_z_0_x_xyyyy, g_z_0_x_xyyyz, g_z_0_x_xyyzz, g_z_0_x_xyzzz, g_z_0_x_xzzzz, g_z_0_x_yyyyy, g_z_0_x_yyyyz, g_z_0_x_yyyzz, g_z_0_x_yyzzz, g_z_0_x_yzzzz, g_z_0_x_zzzzz, g_z_x_x_xxxxx, g_z_x_x_xxxxxx, g_z_x_x_xxxxxy, g_z_x_x_xxxxxz, g_z_x_x_xxxxy, g_z_x_x_xxxxyy, g_z_x_x_xxxxyz, g_z_x_x_xxxxz, g_z_x_x_xxxxzz, g_z_x_x_xxxyy, g_z_x_x_xxxyyy, g_z_x_x_xxxyyz, g_z_x_x_xxxyz, g_z_x_x_xxxyzz, g_z_x_x_xxxzz, g_z_x_x_xxxzzz, g_z_x_x_xxyyy, g_z_x_x_xxyyyy, g_z_x_x_xxyyyz, g_z_x_x_xxyyz, g_z_x_x_xxyyzz, g_z_x_x_xxyzz, g_z_x_x_xxyzzz, g_z_x_x_xxzzz, g_z_x_x_xxzzzz, g_z_x_x_xyyyy, g_z_x_x_xyyyyy, g_z_x_x_xyyyyz, g_z_x_x_xyyyz, g_z_x_x_xyyyzz, g_z_x_x_xyyzz, g_z_x_x_xyyzzz, g_z_x_x_xyzzz, g_z_x_x_xyzzzz, g_z_x_x_xzzzz, g_z_x_x_xzzzzz, g_z_x_x_yyyyy, g_z_x_x_yyyyz, g_z_x_x_yyyzz, g_z_x_x_yyzzz, g_z_x_x_yzzzz, g_z_x_x_zzzzz, g_z_x_xx_xxxxx, g_z_x_xx_xxxxy, g_z_x_xx_xxxxz, g_z_x_xx_xxxyy, g_z_x_xx_xxxyz, g_z_x_xx_xxxzz, g_z_x_xx_xxyyy, g_z_x_xx_xxyyz, g_z_x_xx_xxyzz, g_z_x_xx_xxzzz, g_z_x_xx_xyyyy, g_z_x_xx_xyyyz, g_z_x_xx_xyyzz, g_z_x_xx_xyzzz, g_z_x_xx_xzzzz, g_z_x_xx_yyyyy, g_z_x_xx_yyyyz, g_z_x_xx_yyyzz, g_z_x_xx_yyzzz, g_z_x_xx_yzzzz, g_z_x_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xx_xxxxx[k] = g_z_0_x_xxxxx[k] - g_z_x_x_xxxxx[k] * ab_x + g_z_x_x_xxxxxx[k];

                g_z_x_xx_xxxxy[k] = g_z_0_x_xxxxy[k] - g_z_x_x_xxxxy[k] * ab_x + g_z_x_x_xxxxxy[k];

                g_z_x_xx_xxxxz[k] = g_z_0_x_xxxxz[k] - g_z_x_x_xxxxz[k] * ab_x + g_z_x_x_xxxxxz[k];

                g_z_x_xx_xxxyy[k] = g_z_0_x_xxxyy[k] - g_z_x_x_xxxyy[k] * ab_x + g_z_x_x_xxxxyy[k];

                g_z_x_xx_xxxyz[k] = g_z_0_x_xxxyz[k] - g_z_x_x_xxxyz[k] * ab_x + g_z_x_x_xxxxyz[k];

                g_z_x_xx_xxxzz[k] = g_z_0_x_xxxzz[k] - g_z_x_x_xxxzz[k] * ab_x + g_z_x_x_xxxxzz[k];

                g_z_x_xx_xxyyy[k] = g_z_0_x_xxyyy[k] - g_z_x_x_xxyyy[k] * ab_x + g_z_x_x_xxxyyy[k];

                g_z_x_xx_xxyyz[k] = g_z_0_x_xxyyz[k] - g_z_x_x_xxyyz[k] * ab_x + g_z_x_x_xxxyyz[k];

                g_z_x_xx_xxyzz[k] = g_z_0_x_xxyzz[k] - g_z_x_x_xxyzz[k] * ab_x + g_z_x_x_xxxyzz[k];

                g_z_x_xx_xxzzz[k] = g_z_0_x_xxzzz[k] - g_z_x_x_xxzzz[k] * ab_x + g_z_x_x_xxxzzz[k];

                g_z_x_xx_xyyyy[k] = g_z_0_x_xyyyy[k] - g_z_x_x_xyyyy[k] * ab_x + g_z_x_x_xxyyyy[k];

                g_z_x_xx_xyyyz[k] = g_z_0_x_xyyyz[k] - g_z_x_x_xyyyz[k] * ab_x + g_z_x_x_xxyyyz[k];

                g_z_x_xx_xyyzz[k] = g_z_0_x_xyyzz[k] - g_z_x_x_xyyzz[k] * ab_x + g_z_x_x_xxyyzz[k];

                g_z_x_xx_xyzzz[k] = g_z_0_x_xyzzz[k] - g_z_x_x_xyzzz[k] * ab_x + g_z_x_x_xxyzzz[k];

                g_z_x_xx_xzzzz[k] = g_z_0_x_xzzzz[k] - g_z_x_x_xzzzz[k] * ab_x + g_z_x_x_xxzzzz[k];

                g_z_x_xx_yyyyy[k] = g_z_0_x_yyyyy[k] - g_z_x_x_yyyyy[k] * ab_x + g_z_x_x_xyyyyy[k];

                g_z_x_xx_yyyyz[k] = g_z_0_x_yyyyz[k] - g_z_x_x_yyyyz[k] * ab_x + g_z_x_x_xyyyyz[k];

                g_z_x_xx_yyyzz[k] = g_z_0_x_yyyzz[k] - g_z_x_x_yyyzz[k] * ab_x + g_z_x_x_xyyyzz[k];

                g_z_x_xx_yyzzz[k] = g_z_0_x_yyzzz[k] - g_z_x_x_yyzzz[k] * ab_x + g_z_x_x_xyyzzz[k];

                g_z_x_xx_yzzzz[k] = g_z_0_x_yzzzz[k] - g_z_x_x_yzzzz[k] * ab_x + g_z_x_x_xyzzzz[k];

                g_z_x_xx_zzzzz[k] = g_z_0_x_zzzzz[k] - g_z_x_x_zzzzz[k] * ab_x + g_z_x_x_xzzzzz[k];
            }

            /// Set up 777-798 components of targeted buffer : cbuffer.data(

            auto g_z_x_xy_xxxxx = cbuffer.data(dh_geom_11_off + 777 * ccomps * dcomps);

            auto g_z_x_xy_xxxxy = cbuffer.data(dh_geom_11_off + 778 * ccomps * dcomps);

            auto g_z_x_xy_xxxxz = cbuffer.data(dh_geom_11_off + 779 * ccomps * dcomps);

            auto g_z_x_xy_xxxyy = cbuffer.data(dh_geom_11_off + 780 * ccomps * dcomps);

            auto g_z_x_xy_xxxyz = cbuffer.data(dh_geom_11_off + 781 * ccomps * dcomps);

            auto g_z_x_xy_xxxzz = cbuffer.data(dh_geom_11_off + 782 * ccomps * dcomps);

            auto g_z_x_xy_xxyyy = cbuffer.data(dh_geom_11_off + 783 * ccomps * dcomps);

            auto g_z_x_xy_xxyyz = cbuffer.data(dh_geom_11_off + 784 * ccomps * dcomps);

            auto g_z_x_xy_xxyzz = cbuffer.data(dh_geom_11_off + 785 * ccomps * dcomps);

            auto g_z_x_xy_xxzzz = cbuffer.data(dh_geom_11_off + 786 * ccomps * dcomps);

            auto g_z_x_xy_xyyyy = cbuffer.data(dh_geom_11_off + 787 * ccomps * dcomps);

            auto g_z_x_xy_xyyyz = cbuffer.data(dh_geom_11_off + 788 * ccomps * dcomps);

            auto g_z_x_xy_xyyzz = cbuffer.data(dh_geom_11_off + 789 * ccomps * dcomps);

            auto g_z_x_xy_xyzzz = cbuffer.data(dh_geom_11_off + 790 * ccomps * dcomps);

            auto g_z_x_xy_xzzzz = cbuffer.data(dh_geom_11_off + 791 * ccomps * dcomps);

            auto g_z_x_xy_yyyyy = cbuffer.data(dh_geom_11_off + 792 * ccomps * dcomps);

            auto g_z_x_xy_yyyyz = cbuffer.data(dh_geom_11_off + 793 * ccomps * dcomps);

            auto g_z_x_xy_yyyzz = cbuffer.data(dh_geom_11_off + 794 * ccomps * dcomps);

            auto g_z_x_xy_yyzzz = cbuffer.data(dh_geom_11_off + 795 * ccomps * dcomps);

            auto g_z_x_xy_yzzzz = cbuffer.data(dh_geom_11_off + 796 * ccomps * dcomps);

            auto g_z_x_xy_zzzzz = cbuffer.data(dh_geom_11_off + 797 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_x_xxxxx, g_z_x_x_xxxxxy, g_z_x_x_xxxxy, g_z_x_x_xxxxyy, g_z_x_x_xxxxyz, g_z_x_x_xxxxz, g_z_x_x_xxxyy, g_z_x_x_xxxyyy, g_z_x_x_xxxyyz, g_z_x_x_xxxyz, g_z_x_x_xxxyzz, g_z_x_x_xxxzz, g_z_x_x_xxyyy, g_z_x_x_xxyyyy, g_z_x_x_xxyyyz, g_z_x_x_xxyyz, g_z_x_x_xxyyzz, g_z_x_x_xxyzz, g_z_x_x_xxyzzz, g_z_x_x_xxzzz, g_z_x_x_xyyyy, g_z_x_x_xyyyyy, g_z_x_x_xyyyyz, g_z_x_x_xyyyz, g_z_x_x_xyyyzz, g_z_x_x_xyyzz, g_z_x_x_xyyzzz, g_z_x_x_xyzzz, g_z_x_x_xyzzzz, g_z_x_x_xzzzz, g_z_x_x_yyyyy, g_z_x_x_yyyyyy, g_z_x_x_yyyyyz, g_z_x_x_yyyyz, g_z_x_x_yyyyzz, g_z_x_x_yyyzz, g_z_x_x_yyyzzz, g_z_x_x_yyzzz, g_z_x_x_yyzzzz, g_z_x_x_yzzzz, g_z_x_x_yzzzzz, g_z_x_x_zzzzz, g_z_x_xy_xxxxx, g_z_x_xy_xxxxy, g_z_x_xy_xxxxz, g_z_x_xy_xxxyy, g_z_x_xy_xxxyz, g_z_x_xy_xxxzz, g_z_x_xy_xxyyy, g_z_x_xy_xxyyz, g_z_x_xy_xxyzz, g_z_x_xy_xxzzz, g_z_x_xy_xyyyy, g_z_x_xy_xyyyz, g_z_x_xy_xyyzz, g_z_x_xy_xyzzz, g_z_x_xy_xzzzz, g_z_x_xy_yyyyy, g_z_x_xy_yyyyz, g_z_x_xy_yyyzz, g_z_x_xy_yyzzz, g_z_x_xy_yzzzz, g_z_x_xy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xy_xxxxx[k] = -g_z_x_x_xxxxx[k] * ab_y + g_z_x_x_xxxxxy[k];

                g_z_x_xy_xxxxy[k] = -g_z_x_x_xxxxy[k] * ab_y + g_z_x_x_xxxxyy[k];

                g_z_x_xy_xxxxz[k] = -g_z_x_x_xxxxz[k] * ab_y + g_z_x_x_xxxxyz[k];

                g_z_x_xy_xxxyy[k] = -g_z_x_x_xxxyy[k] * ab_y + g_z_x_x_xxxyyy[k];

                g_z_x_xy_xxxyz[k] = -g_z_x_x_xxxyz[k] * ab_y + g_z_x_x_xxxyyz[k];

                g_z_x_xy_xxxzz[k] = -g_z_x_x_xxxzz[k] * ab_y + g_z_x_x_xxxyzz[k];

                g_z_x_xy_xxyyy[k] = -g_z_x_x_xxyyy[k] * ab_y + g_z_x_x_xxyyyy[k];

                g_z_x_xy_xxyyz[k] = -g_z_x_x_xxyyz[k] * ab_y + g_z_x_x_xxyyyz[k];

                g_z_x_xy_xxyzz[k] = -g_z_x_x_xxyzz[k] * ab_y + g_z_x_x_xxyyzz[k];

                g_z_x_xy_xxzzz[k] = -g_z_x_x_xxzzz[k] * ab_y + g_z_x_x_xxyzzz[k];

                g_z_x_xy_xyyyy[k] = -g_z_x_x_xyyyy[k] * ab_y + g_z_x_x_xyyyyy[k];

                g_z_x_xy_xyyyz[k] = -g_z_x_x_xyyyz[k] * ab_y + g_z_x_x_xyyyyz[k];

                g_z_x_xy_xyyzz[k] = -g_z_x_x_xyyzz[k] * ab_y + g_z_x_x_xyyyzz[k];

                g_z_x_xy_xyzzz[k] = -g_z_x_x_xyzzz[k] * ab_y + g_z_x_x_xyyzzz[k];

                g_z_x_xy_xzzzz[k] = -g_z_x_x_xzzzz[k] * ab_y + g_z_x_x_xyzzzz[k];

                g_z_x_xy_yyyyy[k] = -g_z_x_x_yyyyy[k] * ab_y + g_z_x_x_yyyyyy[k];

                g_z_x_xy_yyyyz[k] = -g_z_x_x_yyyyz[k] * ab_y + g_z_x_x_yyyyyz[k];

                g_z_x_xy_yyyzz[k] = -g_z_x_x_yyyzz[k] * ab_y + g_z_x_x_yyyyzz[k];

                g_z_x_xy_yyzzz[k] = -g_z_x_x_yyzzz[k] * ab_y + g_z_x_x_yyyzzz[k];

                g_z_x_xy_yzzzz[k] = -g_z_x_x_yzzzz[k] * ab_y + g_z_x_x_yyzzzz[k];

                g_z_x_xy_zzzzz[k] = -g_z_x_x_zzzzz[k] * ab_y + g_z_x_x_yzzzzz[k];
            }

            /// Set up 798-819 components of targeted buffer : cbuffer.data(

            auto g_z_x_xz_xxxxx = cbuffer.data(dh_geom_11_off + 798 * ccomps * dcomps);

            auto g_z_x_xz_xxxxy = cbuffer.data(dh_geom_11_off + 799 * ccomps * dcomps);

            auto g_z_x_xz_xxxxz = cbuffer.data(dh_geom_11_off + 800 * ccomps * dcomps);

            auto g_z_x_xz_xxxyy = cbuffer.data(dh_geom_11_off + 801 * ccomps * dcomps);

            auto g_z_x_xz_xxxyz = cbuffer.data(dh_geom_11_off + 802 * ccomps * dcomps);

            auto g_z_x_xz_xxxzz = cbuffer.data(dh_geom_11_off + 803 * ccomps * dcomps);

            auto g_z_x_xz_xxyyy = cbuffer.data(dh_geom_11_off + 804 * ccomps * dcomps);

            auto g_z_x_xz_xxyyz = cbuffer.data(dh_geom_11_off + 805 * ccomps * dcomps);

            auto g_z_x_xz_xxyzz = cbuffer.data(dh_geom_11_off + 806 * ccomps * dcomps);

            auto g_z_x_xz_xxzzz = cbuffer.data(dh_geom_11_off + 807 * ccomps * dcomps);

            auto g_z_x_xz_xyyyy = cbuffer.data(dh_geom_11_off + 808 * ccomps * dcomps);

            auto g_z_x_xz_xyyyz = cbuffer.data(dh_geom_11_off + 809 * ccomps * dcomps);

            auto g_z_x_xz_xyyzz = cbuffer.data(dh_geom_11_off + 810 * ccomps * dcomps);

            auto g_z_x_xz_xyzzz = cbuffer.data(dh_geom_11_off + 811 * ccomps * dcomps);

            auto g_z_x_xz_xzzzz = cbuffer.data(dh_geom_11_off + 812 * ccomps * dcomps);

            auto g_z_x_xz_yyyyy = cbuffer.data(dh_geom_11_off + 813 * ccomps * dcomps);

            auto g_z_x_xz_yyyyz = cbuffer.data(dh_geom_11_off + 814 * ccomps * dcomps);

            auto g_z_x_xz_yyyzz = cbuffer.data(dh_geom_11_off + 815 * ccomps * dcomps);

            auto g_z_x_xz_yyzzz = cbuffer.data(dh_geom_11_off + 816 * ccomps * dcomps);

            auto g_z_x_xz_yzzzz = cbuffer.data(dh_geom_11_off + 817 * ccomps * dcomps);

            auto g_z_x_xz_zzzzz = cbuffer.data(dh_geom_11_off + 818 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxxx, g_z_0_z_xxxxy, g_z_0_z_xxxxz, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxzz, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyzz, g_z_0_z_xxzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyzz, g_z_0_z_xyzzz, g_z_0_z_xzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyzz, g_z_0_z_yyzzz, g_z_0_z_yzzzz, g_z_0_z_zzzzz, g_z_x_xz_xxxxx, g_z_x_xz_xxxxy, g_z_x_xz_xxxxz, g_z_x_xz_xxxyy, g_z_x_xz_xxxyz, g_z_x_xz_xxxzz, g_z_x_xz_xxyyy, g_z_x_xz_xxyyz, g_z_x_xz_xxyzz, g_z_x_xz_xxzzz, g_z_x_xz_xyyyy, g_z_x_xz_xyyyz, g_z_x_xz_xyyzz, g_z_x_xz_xyzzz, g_z_x_xz_xzzzz, g_z_x_xz_yyyyy, g_z_x_xz_yyyyz, g_z_x_xz_yyyzz, g_z_x_xz_yyzzz, g_z_x_xz_yzzzz, g_z_x_xz_zzzzz, g_z_x_z_xxxxx, g_z_x_z_xxxxxx, g_z_x_z_xxxxxy, g_z_x_z_xxxxxz, g_z_x_z_xxxxy, g_z_x_z_xxxxyy, g_z_x_z_xxxxyz, g_z_x_z_xxxxz, g_z_x_z_xxxxzz, g_z_x_z_xxxyy, g_z_x_z_xxxyyy, g_z_x_z_xxxyyz, g_z_x_z_xxxyz, g_z_x_z_xxxyzz, g_z_x_z_xxxzz, g_z_x_z_xxxzzz, g_z_x_z_xxyyy, g_z_x_z_xxyyyy, g_z_x_z_xxyyyz, g_z_x_z_xxyyz, g_z_x_z_xxyyzz, g_z_x_z_xxyzz, g_z_x_z_xxyzzz, g_z_x_z_xxzzz, g_z_x_z_xxzzzz, g_z_x_z_xyyyy, g_z_x_z_xyyyyy, g_z_x_z_xyyyyz, g_z_x_z_xyyyz, g_z_x_z_xyyyzz, g_z_x_z_xyyzz, g_z_x_z_xyyzzz, g_z_x_z_xyzzz, g_z_x_z_xyzzzz, g_z_x_z_xzzzz, g_z_x_z_xzzzzz, g_z_x_z_yyyyy, g_z_x_z_yyyyz, g_z_x_z_yyyzz, g_z_x_z_yyzzz, g_z_x_z_yzzzz, g_z_x_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_xz_xxxxx[k] = g_z_0_z_xxxxx[k] - g_z_x_z_xxxxx[k] * ab_x + g_z_x_z_xxxxxx[k];

                g_z_x_xz_xxxxy[k] = g_z_0_z_xxxxy[k] - g_z_x_z_xxxxy[k] * ab_x + g_z_x_z_xxxxxy[k];

                g_z_x_xz_xxxxz[k] = g_z_0_z_xxxxz[k] - g_z_x_z_xxxxz[k] * ab_x + g_z_x_z_xxxxxz[k];

                g_z_x_xz_xxxyy[k] = g_z_0_z_xxxyy[k] - g_z_x_z_xxxyy[k] * ab_x + g_z_x_z_xxxxyy[k];

                g_z_x_xz_xxxyz[k] = g_z_0_z_xxxyz[k] - g_z_x_z_xxxyz[k] * ab_x + g_z_x_z_xxxxyz[k];

                g_z_x_xz_xxxzz[k] = g_z_0_z_xxxzz[k] - g_z_x_z_xxxzz[k] * ab_x + g_z_x_z_xxxxzz[k];

                g_z_x_xz_xxyyy[k] = g_z_0_z_xxyyy[k] - g_z_x_z_xxyyy[k] * ab_x + g_z_x_z_xxxyyy[k];

                g_z_x_xz_xxyyz[k] = g_z_0_z_xxyyz[k] - g_z_x_z_xxyyz[k] * ab_x + g_z_x_z_xxxyyz[k];

                g_z_x_xz_xxyzz[k] = g_z_0_z_xxyzz[k] - g_z_x_z_xxyzz[k] * ab_x + g_z_x_z_xxxyzz[k];

                g_z_x_xz_xxzzz[k] = g_z_0_z_xxzzz[k] - g_z_x_z_xxzzz[k] * ab_x + g_z_x_z_xxxzzz[k];

                g_z_x_xz_xyyyy[k] = g_z_0_z_xyyyy[k] - g_z_x_z_xyyyy[k] * ab_x + g_z_x_z_xxyyyy[k];

                g_z_x_xz_xyyyz[k] = g_z_0_z_xyyyz[k] - g_z_x_z_xyyyz[k] * ab_x + g_z_x_z_xxyyyz[k];

                g_z_x_xz_xyyzz[k] = g_z_0_z_xyyzz[k] - g_z_x_z_xyyzz[k] * ab_x + g_z_x_z_xxyyzz[k];

                g_z_x_xz_xyzzz[k] = g_z_0_z_xyzzz[k] - g_z_x_z_xyzzz[k] * ab_x + g_z_x_z_xxyzzz[k];

                g_z_x_xz_xzzzz[k] = g_z_0_z_xzzzz[k] - g_z_x_z_xzzzz[k] * ab_x + g_z_x_z_xxzzzz[k];

                g_z_x_xz_yyyyy[k] = g_z_0_z_yyyyy[k] - g_z_x_z_yyyyy[k] * ab_x + g_z_x_z_xyyyyy[k];

                g_z_x_xz_yyyyz[k] = g_z_0_z_yyyyz[k] - g_z_x_z_yyyyz[k] * ab_x + g_z_x_z_xyyyyz[k];

                g_z_x_xz_yyyzz[k] = g_z_0_z_yyyzz[k] - g_z_x_z_yyyzz[k] * ab_x + g_z_x_z_xyyyzz[k];

                g_z_x_xz_yyzzz[k] = g_z_0_z_yyzzz[k] - g_z_x_z_yyzzz[k] * ab_x + g_z_x_z_xyyzzz[k];

                g_z_x_xz_yzzzz[k] = g_z_0_z_yzzzz[k] - g_z_x_z_yzzzz[k] * ab_x + g_z_x_z_xyzzzz[k];

                g_z_x_xz_zzzzz[k] = g_z_0_z_zzzzz[k] - g_z_x_z_zzzzz[k] * ab_x + g_z_x_z_xzzzzz[k];
            }

            /// Set up 819-840 components of targeted buffer : cbuffer.data(

            auto g_z_x_yy_xxxxx = cbuffer.data(dh_geom_11_off + 819 * ccomps * dcomps);

            auto g_z_x_yy_xxxxy = cbuffer.data(dh_geom_11_off + 820 * ccomps * dcomps);

            auto g_z_x_yy_xxxxz = cbuffer.data(dh_geom_11_off + 821 * ccomps * dcomps);

            auto g_z_x_yy_xxxyy = cbuffer.data(dh_geom_11_off + 822 * ccomps * dcomps);

            auto g_z_x_yy_xxxyz = cbuffer.data(dh_geom_11_off + 823 * ccomps * dcomps);

            auto g_z_x_yy_xxxzz = cbuffer.data(dh_geom_11_off + 824 * ccomps * dcomps);

            auto g_z_x_yy_xxyyy = cbuffer.data(dh_geom_11_off + 825 * ccomps * dcomps);

            auto g_z_x_yy_xxyyz = cbuffer.data(dh_geom_11_off + 826 * ccomps * dcomps);

            auto g_z_x_yy_xxyzz = cbuffer.data(dh_geom_11_off + 827 * ccomps * dcomps);

            auto g_z_x_yy_xxzzz = cbuffer.data(dh_geom_11_off + 828 * ccomps * dcomps);

            auto g_z_x_yy_xyyyy = cbuffer.data(dh_geom_11_off + 829 * ccomps * dcomps);

            auto g_z_x_yy_xyyyz = cbuffer.data(dh_geom_11_off + 830 * ccomps * dcomps);

            auto g_z_x_yy_xyyzz = cbuffer.data(dh_geom_11_off + 831 * ccomps * dcomps);

            auto g_z_x_yy_xyzzz = cbuffer.data(dh_geom_11_off + 832 * ccomps * dcomps);

            auto g_z_x_yy_xzzzz = cbuffer.data(dh_geom_11_off + 833 * ccomps * dcomps);

            auto g_z_x_yy_yyyyy = cbuffer.data(dh_geom_11_off + 834 * ccomps * dcomps);

            auto g_z_x_yy_yyyyz = cbuffer.data(dh_geom_11_off + 835 * ccomps * dcomps);

            auto g_z_x_yy_yyyzz = cbuffer.data(dh_geom_11_off + 836 * ccomps * dcomps);

            auto g_z_x_yy_yyzzz = cbuffer.data(dh_geom_11_off + 837 * ccomps * dcomps);

            auto g_z_x_yy_yzzzz = cbuffer.data(dh_geom_11_off + 838 * ccomps * dcomps);

            auto g_z_x_yy_zzzzz = cbuffer.data(dh_geom_11_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_y_xxxxx, g_z_x_y_xxxxxy, g_z_x_y_xxxxy, g_z_x_y_xxxxyy, g_z_x_y_xxxxyz, g_z_x_y_xxxxz, g_z_x_y_xxxyy, g_z_x_y_xxxyyy, g_z_x_y_xxxyyz, g_z_x_y_xxxyz, g_z_x_y_xxxyzz, g_z_x_y_xxxzz, g_z_x_y_xxyyy, g_z_x_y_xxyyyy, g_z_x_y_xxyyyz, g_z_x_y_xxyyz, g_z_x_y_xxyyzz, g_z_x_y_xxyzz, g_z_x_y_xxyzzz, g_z_x_y_xxzzz, g_z_x_y_xyyyy, g_z_x_y_xyyyyy, g_z_x_y_xyyyyz, g_z_x_y_xyyyz, g_z_x_y_xyyyzz, g_z_x_y_xyyzz, g_z_x_y_xyyzzz, g_z_x_y_xyzzz, g_z_x_y_xyzzzz, g_z_x_y_xzzzz, g_z_x_y_yyyyy, g_z_x_y_yyyyyy, g_z_x_y_yyyyyz, g_z_x_y_yyyyz, g_z_x_y_yyyyzz, g_z_x_y_yyyzz, g_z_x_y_yyyzzz, g_z_x_y_yyzzz, g_z_x_y_yyzzzz, g_z_x_y_yzzzz, g_z_x_y_yzzzzz, g_z_x_y_zzzzz, g_z_x_yy_xxxxx, g_z_x_yy_xxxxy, g_z_x_yy_xxxxz, g_z_x_yy_xxxyy, g_z_x_yy_xxxyz, g_z_x_yy_xxxzz, g_z_x_yy_xxyyy, g_z_x_yy_xxyyz, g_z_x_yy_xxyzz, g_z_x_yy_xxzzz, g_z_x_yy_xyyyy, g_z_x_yy_xyyyz, g_z_x_yy_xyyzz, g_z_x_yy_xyzzz, g_z_x_yy_xzzzz, g_z_x_yy_yyyyy, g_z_x_yy_yyyyz, g_z_x_yy_yyyzz, g_z_x_yy_yyzzz, g_z_x_yy_yzzzz, g_z_x_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yy_xxxxx[k] = -g_z_x_y_xxxxx[k] * ab_y + g_z_x_y_xxxxxy[k];

                g_z_x_yy_xxxxy[k] = -g_z_x_y_xxxxy[k] * ab_y + g_z_x_y_xxxxyy[k];

                g_z_x_yy_xxxxz[k] = -g_z_x_y_xxxxz[k] * ab_y + g_z_x_y_xxxxyz[k];

                g_z_x_yy_xxxyy[k] = -g_z_x_y_xxxyy[k] * ab_y + g_z_x_y_xxxyyy[k];

                g_z_x_yy_xxxyz[k] = -g_z_x_y_xxxyz[k] * ab_y + g_z_x_y_xxxyyz[k];

                g_z_x_yy_xxxzz[k] = -g_z_x_y_xxxzz[k] * ab_y + g_z_x_y_xxxyzz[k];

                g_z_x_yy_xxyyy[k] = -g_z_x_y_xxyyy[k] * ab_y + g_z_x_y_xxyyyy[k];

                g_z_x_yy_xxyyz[k] = -g_z_x_y_xxyyz[k] * ab_y + g_z_x_y_xxyyyz[k];

                g_z_x_yy_xxyzz[k] = -g_z_x_y_xxyzz[k] * ab_y + g_z_x_y_xxyyzz[k];

                g_z_x_yy_xxzzz[k] = -g_z_x_y_xxzzz[k] * ab_y + g_z_x_y_xxyzzz[k];

                g_z_x_yy_xyyyy[k] = -g_z_x_y_xyyyy[k] * ab_y + g_z_x_y_xyyyyy[k];

                g_z_x_yy_xyyyz[k] = -g_z_x_y_xyyyz[k] * ab_y + g_z_x_y_xyyyyz[k];

                g_z_x_yy_xyyzz[k] = -g_z_x_y_xyyzz[k] * ab_y + g_z_x_y_xyyyzz[k];

                g_z_x_yy_xyzzz[k] = -g_z_x_y_xyzzz[k] * ab_y + g_z_x_y_xyyzzz[k];

                g_z_x_yy_xzzzz[k] = -g_z_x_y_xzzzz[k] * ab_y + g_z_x_y_xyzzzz[k];

                g_z_x_yy_yyyyy[k] = -g_z_x_y_yyyyy[k] * ab_y + g_z_x_y_yyyyyy[k];

                g_z_x_yy_yyyyz[k] = -g_z_x_y_yyyyz[k] * ab_y + g_z_x_y_yyyyyz[k];

                g_z_x_yy_yyyzz[k] = -g_z_x_y_yyyzz[k] * ab_y + g_z_x_y_yyyyzz[k];

                g_z_x_yy_yyzzz[k] = -g_z_x_y_yyzzz[k] * ab_y + g_z_x_y_yyyzzz[k];

                g_z_x_yy_yzzzz[k] = -g_z_x_y_yzzzz[k] * ab_y + g_z_x_y_yyzzzz[k];

                g_z_x_yy_zzzzz[k] = -g_z_x_y_zzzzz[k] * ab_y + g_z_x_y_yzzzzz[k];
            }

            /// Set up 840-861 components of targeted buffer : cbuffer.data(

            auto g_z_x_yz_xxxxx = cbuffer.data(dh_geom_11_off + 840 * ccomps * dcomps);

            auto g_z_x_yz_xxxxy = cbuffer.data(dh_geom_11_off + 841 * ccomps * dcomps);

            auto g_z_x_yz_xxxxz = cbuffer.data(dh_geom_11_off + 842 * ccomps * dcomps);

            auto g_z_x_yz_xxxyy = cbuffer.data(dh_geom_11_off + 843 * ccomps * dcomps);

            auto g_z_x_yz_xxxyz = cbuffer.data(dh_geom_11_off + 844 * ccomps * dcomps);

            auto g_z_x_yz_xxxzz = cbuffer.data(dh_geom_11_off + 845 * ccomps * dcomps);

            auto g_z_x_yz_xxyyy = cbuffer.data(dh_geom_11_off + 846 * ccomps * dcomps);

            auto g_z_x_yz_xxyyz = cbuffer.data(dh_geom_11_off + 847 * ccomps * dcomps);

            auto g_z_x_yz_xxyzz = cbuffer.data(dh_geom_11_off + 848 * ccomps * dcomps);

            auto g_z_x_yz_xxzzz = cbuffer.data(dh_geom_11_off + 849 * ccomps * dcomps);

            auto g_z_x_yz_xyyyy = cbuffer.data(dh_geom_11_off + 850 * ccomps * dcomps);

            auto g_z_x_yz_xyyyz = cbuffer.data(dh_geom_11_off + 851 * ccomps * dcomps);

            auto g_z_x_yz_xyyzz = cbuffer.data(dh_geom_11_off + 852 * ccomps * dcomps);

            auto g_z_x_yz_xyzzz = cbuffer.data(dh_geom_11_off + 853 * ccomps * dcomps);

            auto g_z_x_yz_xzzzz = cbuffer.data(dh_geom_11_off + 854 * ccomps * dcomps);

            auto g_z_x_yz_yyyyy = cbuffer.data(dh_geom_11_off + 855 * ccomps * dcomps);

            auto g_z_x_yz_yyyyz = cbuffer.data(dh_geom_11_off + 856 * ccomps * dcomps);

            auto g_z_x_yz_yyyzz = cbuffer.data(dh_geom_11_off + 857 * ccomps * dcomps);

            auto g_z_x_yz_yyzzz = cbuffer.data(dh_geom_11_off + 858 * ccomps * dcomps);

            auto g_z_x_yz_yzzzz = cbuffer.data(dh_geom_11_off + 859 * ccomps * dcomps);

            auto g_z_x_yz_zzzzz = cbuffer.data(dh_geom_11_off + 860 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_x_yz_xxxxx, g_z_x_yz_xxxxy, g_z_x_yz_xxxxz, g_z_x_yz_xxxyy, g_z_x_yz_xxxyz, g_z_x_yz_xxxzz, g_z_x_yz_xxyyy, g_z_x_yz_xxyyz, g_z_x_yz_xxyzz, g_z_x_yz_xxzzz, g_z_x_yz_xyyyy, g_z_x_yz_xyyyz, g_z_x_yz_xyyzz, g_z_x_yz_xyzzz, g_z_x_yz_xzzzz, g_z_x_yz_yyyyy, g_z_x_yz_yyyyz, g_z_x_yz_yyyzz, g_z_x_yz_yyzzz, g_z_x_yz_yzzzz, g_z_x_yz_zzzzz, g_z_x_z_xxxxx, g_z_x_z_xxxxxy, g_z_x_z_xxxxy, g_z_x_z_xxxxyy, g_z_x_z_xxxxyz, g_z_x_z_xxxxz, g_z_x_z_xxxyy, g_z_x_z_xxxyyy, g_z_x_z_xxxyyz, g_z_x_z_xxxyz, g_z_x_z_xxxyzz, g_z_x_z_xxxzz, g_z_x_z_xxyyy, g_z_x_z_xxyyyy, g_z_x_z_xxyyyz, g_z_x_z_xxyyz, g_z_x_z_xxyyzz, g_z_x_z_xxyzz, g_z_x_z_xxyzzz, g_z_x_z_xxzzz, g_z_x_z_xyyyy, g_z_x_z_xyyyyy, g_z_x_z_xyyyyz, g_z_x_z_xyyyz, g_z_x_z_xyyyzz, g_z_x_z_xyyzz, g_z_x_z_xyyzzz, g_z_x_z_xyzzz, g_z_x_z_xyzzzz, g_z_x_z_xzzzz, g_z_x_z_yyyyy, g_z_x_z_yyyyyy, g_z_x_z_yyyyyz, g_z_x_z_yyyyz, g_z_x_z_yyyyzz, g_z_x_z_yyyzz, g_z_x_z_yyyzzz, g_z_x_z_yyzzz, g_z_x_z_yyzzzz, g_z_x_z_yzzzz, g_z_x_z_yzzzzz, g_z_x_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_yz_xxxxx[k] = -g_z_x_z_xxxxx[k] * ab_y + g_z_x_z_xxxxxy[k];

                g_z_x_yz_xxxxy[k] = -g_z_x_z_xxxxy[k] * ab_y + g_z_x_z_xxxxyy[k];

                g_z_x_yz_xxxxz[k] = -g_z_x_z_xxxxz[k] * ab_y + g_z_x_z_xxxxyz[k];

                g_z_x_yz_xxxyy[k] = -g_z_x_z_xxxyy[k] * ab_y + g_z_x_z_xxxyyy[k];

                g_z_x_yz_xxxyz[k] = -g_z_x_z_xxxyz[k] * ab_y + g_z_x_z_xxxyyz[k];

                g_z_x_yz_xxxzz[k] = -g_z_x_z_xxxzz[k] * ab_y + g_z_x_z_xxxyzz[k];

                g_z_x_yz_xxyyy[k] = -g_z_x_z_xxyyy[k] * ab_y + g_z_x_z_xxyyyy[k];

                g_z_x_yz_xxyyz[k] = -g_z_x_z_xxyyz[k] * ab_y + g_z_x_z_xxyyyz[k];

                g_z_x_yz_xxyzz[k] = -g_z_x_z_xxyzz[k] * ab_y + g_z_x_z_xxyyzz[k];

                g_z_x_yz_xxzzz[k] = -g_z_x_z_xxzzz[k] * ab_y + g_z_x_z_xxyzzz[k];

                g_z_x_yz_xyyyy[k] = -g_z_x_z_xyyyy[k] * ab_y + g_z_x_z_xyyyyy[k];

                g_z_x_yz_xyyyz[k] = -g_z_x_z_xyyyz[k] * ab_y + g_z_x_z_xyyyyz[k];

                g_z_x_yz_xyyzz[k] = -g_z_x_z_xyyzz[k] * ab_y + g_z_x_z_xyyyzz[k];

                g_z_x_yz_xyzzz[k] = -g_z_x_z_xyzzz[k] * ab_y + g_z_x_z_xyyzzz[k];

                g_z_x_yz_xzzzz[k] = -g_z_x_z_xzzzz[k] * ab_y + g_z_x_z_xyzzzz[k];

                g_z_x_yz_yyyyy[k] = -g_z_x_z_yyyyy[k] * ab_y + g_z_x_z_yyyyyy[k];

                g_z_x_yz_yyyyz[k] = -g_z_x_z_yyyyz[k] * ab_y + g_z_x_z_yyyyyz[k];

                g_z_x_yz_yyyzz[k] = -g_z_x_z_yyyzz[k] * ab_y + g_z_x_z_yyyyzz[k];

                g_z_x_yz_yyzzz[k] = -g_z_x_z_yyzzz[k] * ab_y + g_z_x_z_yyyzzz[k];

                g_z_x_yz_yzzzz[k] = -g_z_x_z_yzzzz[k] * ab_y + g_z_x_z_yyzzzz[k];

                g_z_x_yz_zzzzz[k] = -g_z_x_z_zzzzz[k] * ab_y + g_z_x_z_yzzzzz[k];
            }

            /// Set up 861-882 components of targeted buffer : cbuffer.data(

            auto g_z_x_zz_xxxxx = cbuffer.data(dh_geom_11_off + 861 * ccomps * dcomps);

            auto g_z_x_zz_xxxxy = cbuffer.data(dh_geom_11_off + 862 * ccomps * dcomps);

            auto g_z_x_zz_xxxxz = cbuffer.data(dh_geom_11_off + 863 * ccomps * dcomps);

            auto g_z_x_zz_xxxyy = cbuffer.data(dh_geom_11_off + 864 * ccomps * dcomps);

            auto g_z_x_zz_xxxyz = cbuffer.data(dh_geom_11_off + 865 * ccomps * dcomps);

            auto g_z_x_zz_xxxzz = cbuffer.data(dh_geom_11_off + 866 * ccomps * dcomps);

            auto g_z_x_zz_xxyyy = cbuffer.data(dh_geom_11_off + 867 * ccomps * dcomps);

            auto g_z_x_zz_xxyyz = cbuffer.data(dh_geom_11_off + 868 * ccomps * dcomps);

            auto g_z_x_zz_xxyzz = cbuffer.data(dh_geom_11_off + 869 * ccomps * dcomps);

            auto g_z_x_zz_xxzzz = cbuffer.data(dh_geom_11_off + 870 * ccomps * dcomps);

            auto g_z_x_zz_xyyyy = cbuffer.data(dh_geom_11_off + 871 * ccomps * dcomps);

            auto g_z_x_zz_xyyyz = cbuffer.data(dh_geom_11_off + 872 * ccomps * dcomps);

            auto g_z_x_zz_xyyzz = cbuffer.data(dh_geom_11_off + 873 * ccomps * dcomps);

            auto g_z_x_zz_xyzzz = cbuffer.data(dh_geom_11_off + 874 * ccomps * dcomps);

            auto g_z_x_zz_xzzzz = cbuffer.data(dh_geom_11_off + 875 * ccomps * dcomps);

            auto g_z_x_zz_yyyyy = cbuffer.data(dh_geom_11_off + 876 * ccomps * dcomps);

            auto g_z_x_zz_yyyyz = cbuffer.data(dh_geom_11_off + 877 * ccomps * dcomps);

            auto g_z_x_zz_yyyzz = cbuffer.data(dh_geom_11_off + 878 * ccomps * dcomps);

            auto g_z_x_zz_yyzzz = cbuffer.data(dh_geom_11_off + 879 * ccomps * dcomps);

            auto g_z_x_zz_yzzzz = cbuffer.data(dh_geom_11_off + 880 * ccomps * dcomps);

            auto g_z_x_zz_zzzzz = cbuffer.data(dh_geom_11_off + 881 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_z_xxxxx, g_0_x_z_xxxxy, g_0_x_z_xxxxz, g_0_x_z_xxxyy, g_0_x_z_xxxyz, g_0_x_z_xxxzz, g_0_x_z_xxyyy, g_0_x_z_xxyyz, g_0_x_z_xxyzz, g_0_x_z_xxzzz, g_0_x_z_xyyyy, g_0_x_z_xyyyz, g_0_x_z_xyyzz, g_0_x_z_xyzzz, g_0_x_z_xzzzz, g_0_x_z_yyyyy, g_0_x_z_yyyyz, g_0_x_z_yyyzz, g_0_x_z_yyzzz, g_0_x_z_yzzzz, g_0_x_z_zzzzz, g_z_x_z_xxxxx, g_z_x_z_xxxxxz, g_z_x_z_xxxxy, g_z_x_z_xxxxyz, g_z_x_z_xxxxz, g_z_x_z_xxxxzz, g_z_x_z_xxxyy, g_z_x_z_xxxyyz, g_z_x_z_xxxyz, g_z_x_z_xxxyzz, g_z_x_z_xxxzz, g_z_x_z_xxxzzz, g_z_x_z_xxyyy, g_z_x_z_xxyyyz, g_z_x_z_xxyyz, g_z_x_z_xxyyzz, g_z_x_z_xxyzz, g_z_x_z_xxyzzz, g_z_x_z_xxzzz, g_z_x_z_xxzzzz, g_z_x_z_xyyyy, g_z_x_z_xyyyyz, g_z_x_z_xyyyz, g_z_x_z_xyyyzz, g_z_x_z_xyyzz, g_z_x_z_xyyzzz, g_z_x_z_xyzzz, g_z_x_z_xyzzzz, g_z_x_z_xzzzz, g_z_x_z_xzzzzz, g_z_x_z_yyyyy, g_z_x_z_yyyyyz, g_z_x_z_yyyyz, g_z_x_z_yyyyzz, g_z_x_z_yyyzz, g_z_x_z_yyyzzz, g_z_x_z_yyzzz, g_z_x_z_yyzzzz, g_z_x_z_yzzzz, g_z_x_z_yzzzzz, g_z_x_z_zzzzz, g_z_x_z_zzzzzz, g_z_x_zz_xxxxx, g_z_x_zz_xxxxy, g_z_x_zz_xxxxz, g_z_x_zz_xxxyy, g_z_x_zz_xxxyz, g_z_x_zz_xxxzz, g_z_x_zz_xxyyy, g_z_x_zz_xxyyz, g_z_x_zz_xxyzz, g_z_x_zz_xxzzz, g_z_x_zz_xyyyy, g_z_x_zz_xyyyz, g_z_x_zz_xyyzz, g_z_x_zz_xyzzz, g_z_x_zz_xzzzz, g_z_x_zz_yyyyy, g_z_x_zz_yyyyz, g_z_x_zz_yyyzz, g_z_x_zz_yyzzz, g_z_x_zz_yzzzz, g_z_x_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_x_zz_xxxxx[k] = -g_0_x_z_xxxxx[k] - g_z_x_z_xxxxx[k] * ab_z + g_z_x_z_xxxxxz[k];

                g_z_x_zz_xxxxy[k] = -g_0_x_z_xxxxy[k] - g_z_x_z_xxxxy[k] * ab_z + g_z_x_z_xxxxyz[k];

                g_z_x_zz_xxxxz[k] = -g_0_x_z_xxxxz[k] - g_z_x_z_xxxxz[k] * ab_z + g_z_x_z_xxxxzz[k];

                g_z_x_zz_xxxyy[k] = -g_0_x_z_xxxyy[k] - g_z_x_z_xxxyy[k] * ab_z + g_z_x_z_xxxyyz[k];

                g_z_x_zz_xxxyz[k] = -g_0_x_z_xxxyz[k] - g_z_x_z_xxxyz[k] * ab_z + g_z_x_z_xxxyzz[k];

                g_z_x_zz_xxxzz[k] = -g_0_x_z_xxxzz[k] - g_z_x_z_xxxzz[k] * ab_z + g_z_x_z_xxxzzz[k];

                g_z_x_zz_xxyyy[k] = -g_0_x_z_xxyyy[k] - g_z_x_z_xxyyy[k] * ab_z + g_z_x_z_xxyyyz[k];

                g_z_x_zz_xxyyz[k] = -g_0_x_z_xxyyz[k] - g_z_x_z_xxyyz[k] * ab_z + g_z_x_z_xxyyzz[k];

                g_z_x_zz_xxyzz[k] = -g_0_x_z_xxyzz[k] - g_z_x_z_xxyzz[k] * ab_z + g_z_x_z_xxyzzz[k];

                g_z_x_zz_xxzzz[k] = -g_0_x_z_xxzzz[k] - g_z_x_z_xxzzz[k] * ab_z + g_z_x_z_xxzzzz[k];

                g_z_x_zz_xyyyy[k] = -g_0_x_z_xyyyy[k] - g_z_x_z_xyyyy[k] * ab_z + g_z_x_z_xyyyyz[k];

                g_z_x_zz_xyyyz[k] = -g_0_x_z_xyyyz[k] - g_z_x_z_xyyyz[k] * ab_z + g_z_x_z_xyyyzz[k];

                g_z_x_zz_xyyzz[k] = -g_0_x_z_xyyzz[k] - g_z_x_z_xyyzz[k] * ab_z + g_z_x_z_xyyzzz[k];

                g_z_x_zz_xyzzz[k] = -g_0_x_z_xyzzz[k] - g_z_x_z_xyzzz[k] * ab_z + g_z_x_z_xyzzzz[k];

                g_z_x_zz_xzzzz[k] = -g_0_x_z_xzzzz[k] - g_z_x_z_xzzzz[k] * ab_z + g_z_x_z_xzzzzz[k];

                g_z_x_zz_yyyyy[k] = -g_0_x_z_yyyyy[k] - g_z_x_z_yyyyy[k] * ab_z + g_z_x_z_yyyyyz[k];

                g_z_x_zz_yyyyz[k] = -g_0_x_z_yyyyz[k] - g_z_x_z_yyyyz[k] * ab_z + g_z_x_z_yyyyzz[k];

                g_z_x_zz_yyyzz[k] = -g_0_x_z_yyyzz[k] - g_z_x_z_yyyzz[k] * ab_z + g_z_x_z_yyyzzz[k];

                g_z_x_zz_yyzzz[k] = -g_0_x_z_yyzzz[k] - g_z_x_z_yyzzz[k] * ab_z + g_z_x_z_yyzzzz[k];

                g_z_x_zz_yzzzz[k] = -g_0_x_z_yzzzz[k] - g_z_x_z_yzzzz[k] * ab_z + g_z_x_z_yzzzzz[k];

                g_z_x_zz_zzzzz[k] = -g_0_x_z_zzzzz[k] - g_z_x_z_zzzzz[k] * ab_z + g_z_x_z_zzzzzz[k];
            }

            /// Set up 882-903 components of targeted buffer : cbuffer.data(

            auto g_z_y_xx_xxxxx = cbuffer.data(dh_geom_11_off + 882 * ccomps * dcomps);

            auto g_z_y_xx_xxxxy = cbuffer.data(dh_geom_11_off + 883 * ccomps * dcomps);

            auto g_z_y_xx_xxxxz = cbuffer.data(dh_geom_11_off + 884 * ccomps * dcomps);

            auto g_z_y_xx_xxxyy = cbuffer.data(dh_geom_11_off + 885 * ccomps * dcomps);

            auto g_z_y_xx_xxxyz = cbuffer.data(dh_geom_11_off + 886 * ccomps * dcomps);

            auto g_z_y_xx_xxxzz = cbuffer.data(dh_geom_11_off + 887 * ccomps * dcomps);

            auto g_z_y_xx_xxyyy = cbuffer.data(dh_geom_11_off + 888 * ccomps * dcomps);

            auto g_z_y_xx_xxyyz = cbuffer.data(dh_geom_11_off + 889 * ccomps * dcomps);

            auto g_z_y_xx_xxyzz = cbuffer.data(dh_geom_11_off + 890 * ccomps * dcomps);

            auto g_z_y_xx_xxzzz = cbuffer.data(dh_geom_11_off + 891 * ccomps * dcomps);

            auto g_z_y_xx_xyyyy = cbuffer.data(dh_geom_11_off + 892 * ccomps * dcomps);

            auto g_z_y_xx_xyyyz = cbuffer.data(dh_geom_11_off + 893 * ccomps * dcomps);

            auto g_z_y_xx_xyyzz = cbuffer.data(dh_geom_11_off + 894 * ccomps * dcomps);

            auto g_z_y_xx_xyzzz = cbuffer.data(dh_geom_11_off + 895 * ccomps * dcomps);

            auto g_z_y_xx_xzzzz = cbuffer.data(dh_geom_11_off + 896 * ccomps * dcomps);

            auto g_z_y_xx_yyyyy = cbuffer.data(dh_geom_11_off + 897 * ccomps * dcomps);

            auto g_z_y_xx_yyyyz = cbuffer.data(dh_geom_11_off + 898 * ccomps * dcomps);

            auto g_z_y_xx_yyyzz = cbuffer.data(dh_geom_11_off + 899 * ccomps * dcomps);

            auto g_z_y_xx_yyzzz = cbuffer.data(dh_geom_11_off + 900 * ccomps * dcomps);

            auto g_z_y_xx_yzzzz = cbuffer.data(dh_geom_11_off + 901 * ccomps * dcomps);

            auto g_z_y_xx_zzzzz = cbuffer.data(dh_geom_11_off + 902 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_x_xxxxx, g_z_y_x_xxxxxx, g_z_y_x_xxxxxy, g_z_y_x_xxxxxz, g_z_y_x_xxxxy, g_z_y_x_xxxxyy, g_z_y_x_xxxxyz, g_z_y_x_xxxxz, g_z_y_x_xxxxzz, g_z_y_x_xxxyy, g_z_y_x_xxxyyy, g_z_y_x_xxxyyz, g_z_y_x_xxxyz, g_z_y_x_xxxyzz, g_z_y_x_xxxzz, g_z_y_x_xxxzzz, g_z_y_x_xxyyy, g_z_y_x_xxyyyy, g_z_y_x_xxyyyz, g_z_y_x_xxyyz, g_z_y_x_xxyyzz, g_z_y_x_xxyzz, g_z_y_x_xxyzzz, g_z_y_x_xxzzz, g_z_y_x_xxzzzz, g_z_y_x_xyyyy, g_z_y_x_xyyyyy, g_z_y_x_xyyyyz, g_z_y_x_xyyyz, g_z_y_x_xyyyzz, g_z_y_x_xyyzz, g_z_y_x_xyyzzz, g_z_y_x_xyzzz, g_z_y_x_xyzzzz, g_z_y_x_xzzzz, g_z_y_x_xzzzzz, g_z_y_x_yyyyy, g_z_y_x_yyyyz, g_z_y_x_yyyzz, g_z_y_x_yyzzz, g_z_y_x_yzzzz, g_z_y_x_zzzzz, g_z_y_xx_xxxxx, g_z_y_xx_xxxxy, g_z_y_xx_xxxxz, g_z_y_xx_xxxyy, g_z_y_xx_xxxyz, g_z_y_xx_xxxzz, g_z_y_xx_xxyyy, g_z_y_xx_xxyyz, g_z_y_xx_xxyzz, g_z_y_xx_xxzzz, g_z_y_xx_xyyyy, g_z_y_xx_xyyyz, g_z_y_xx_xyyzz, g_z_y_xx_xyzzz, g_z_y_xx_xzzzz, g_z_y_xx_yyyyy, g_z_y_xx_yyyyz, g_z_y_xx_yyyzz, g_z_y_xx_yyzzz, g_z_y_xx_yzzzz, g_z_y_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xx_xxxxx[k] = -g_z_y_x_xxxxx[k] * ab_x + g_z_y_x_xxxxxx[k];

                g_z_y_xx_xxxxy[k] = -g_z_y_x_xxxxy[k] * ab_x + g_z_y_x_xxxxxy[k];

                g_z_y_xx_xxxxz[k] = -g_z_y_x_xxxxz[k] * ab_x + g_z_y_x_xxxxxz[k];

                g_z_y_xx_xxxyy[k] = -g_z_y_x_xxxyy[k] * ab_x + g_z_y_x_xxxxyy[k];

                g_z_y_xx_xxxyz[k] = -g_z_y_x_xxxyz[k] * ab_x + g_z_y_x_xxxxyz[k];

                g_z_y_xx_xxxzz[k] = -g_z_y_x_xxxzz[k] * ab_x + g_z_y_x_xxxxzz[k];

                g_z_y_xx_xxyyy[k] = -g_z_y_x_xxyyy[k] * ab_x + g_z_y_x_xxxyyy[k];

                g_z_y_xx_xxyyz[k] = -g_z_y_x_xxyyz[k] * ab_x + g_z_y_x_xxxyyz[k];

                g_z_y_xx_xxyzz[k] = -g_z_y_x_xxyzz[k] * ab_x + g_z_y_x_xxxyzz[k];

                g_z_y_xx_xxzzz[k] = -g_z_y_x_xxzzz[k] * ab_x + g_z_y_x_xxxzzz[k];

                g_z_y_xx_xyyyy[k] = -g_z_y_x_xyyyy[k] * ab_x + g_z_y_x_xxyyyy[k];

                g_z_y_xx_xyyyz[k] = -g_z_y_x_xyyyz[k] * ab_x + g_z_y_x_xxyyyz[k];

                g_z_y_xx_xyyzz[k] = -g_z_y_x_xyyzz[k] * ab_x + g_z_y_x_xxyyzz[k];

                g_z_y_xx_xyzzz[k] = -g_z_y_x_xyzzz[k] * ab_x + g_z_y_x_xxyzzz[k];

                g_z_y_xx_xzzzz[k] = -g_z_y_x_xzzzz[k] * ab_x + g_z_y_x_xxzzzz[k];

                g_z_y_xx_yyyyy[k] = -g_z_y_x_yyyyy[k] * ab_x + g_z_y_x_xyyyyy[k];

                g_z_y_xx_yyyyz[k] = -g_z_y_x_yyyyz[k] * ab_x + g_z_y_x_xyyyyz[k];

                g_z_y_xx_yyyzz[k] = -g_z_y_x_yyyzz[k] * ab_x + g_z_y_x_xyyyzz[k];

                g_z_y_xx_yyzzz[k] = -g_z_y_x_yyzzz[k] * ab_x + g_z_y_x_xyyzzz[k];

                g_z_y_xx_yzzzz[k] = -g_z_y_x_yzzzz[k] * ab_x + g_z_y_x_xyzzzz[k];

                g_z_y_xx_zzzzz[k] = -g_z_y_x_zzzzz[k] * ab_x + g_z_y_x_xzzzzz[k];
            }

            /// Set up 903-924 components of targeted buffer : cbuffer.data(

            auto g_z_y_xy_xxxxx = cbuffer.data(dh_geom_11_off + 903 * ccomps * dcomps);

            auto g_z_y_xy_xxxxy = cbuffer.data(dh_geom_11_off + 904 * ccomps * dcomps);

            auto g_z_y_xy_xxxxz = cbuffer.data(dh_geom_11_off + 905 * ccomps * dcomps);

            auto g_z_y_xy_xxxyy = cbuffer.data(dh_geom_11_off + 906 * ccomps * dcomps);

            auto g_z_y_xy_xxxyz = cbuffer.data(dh_geom_11_off + 907 * ccomps * dcomps);

            auto g_z_y_xy_xxxzz = cbuffer.data(dh_geom_11_off + 908 * ccomps * dcomps);

            auto g_z_y_xy_xxyyy = cbuffer.data(dh_geom_11_off + 909 * ccomps * dcomps);

            auto g_z_y_xy_xxyyz = cbuffer.data(dh_geom_11_off + 910 * ccomps * dcomps);

            auto g_z_y_xy_xxyzz = cbuffer.data(dh_geom_11_off + 911 * ccomps * dcomps);

            auto g_z_y_xy_xxzzz = cbuffer.data(dh_geom_11_off + 912 * ccomps * dcomps);

            auto g_z_y_xy_xyyyy = cbuffer.data(dh_geom_11_off + 913 * ccomps * dcomps);

            auto g_z_y_xy_xyyyz = cbuffer.data(dh_geom_11_off + 914 * ccomps * dcomps);

            auto g_z_y_xy_xyyzz = cbuffer.data(dh_geom_11_off + 915 * ccomps * dcomps);

            auto g_z_y_xy_xyzzz = cbuffer.data(dh_geom_11_off + 916 * ccomps * dcomps);

            auto g_z_y_xy_xzzzz = cbuffer.data(dh_geom_11_off + 917 * ccomps * dcomps);

            auto g_z_y_xy_yyyyy = cbuffer.data(dh_geom_11_off + 918 * ccomps * dcomps);

            auto g_z_y_xy_yyyyz = cbuffer.data(dh_geom_11_off + 919 * ccomps * dcomps);

            auto g_z_y_xy_yyyzz = cbuffer.data(dh_geom_11_off + 920 * ccomps * dcomps);

            auto g_z_y_xy_yyzzz = cbuffer.data(dh_geom_11_off + 921 * ccomps * dcomps);

            auto g_z_y_xy_yzzzz = cbuffer.data(dh_geom_11_off + 922 * ccomps * dcomps);

            auto g_z_y_xy_zzzzz = cbuffer.data(dh_geom_11_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xy_xxxxx, g_z_y_xy_xxxxy, g_z_y_xy_xxxxz, g_z_y_xy_xxxyy, g_z_y_xy_xxxyz, g_z_y_xy_xxxzz, g_z_y_xy_xxyyy, g_z_y_xy_xxyyz, g_z_y_xy_xxyzz, g_z_y_xy_xxzzz, g_z_y_xy_xyyyy, g_z_y_xy_xyyyz, g_z_y_xy_xyyzz, g_z_y_xy_xyzzz, g_z_y_xy_xzzzz, g_z_y_xy_yyyyy, g_z_y_xy_yyyyz, g_z_y_xy_yyyzz, g_z_y_xy_yyzzz, g_z_y_xy_yzzzz, g_z_y_xy_zzzzz, g_z_y_y_xxxxx, g_z_y_y_xxxxxx, g_z_y_y_xxxxxy, g_z_y_y_xxxxxz, g_z_y_y_xxxxy, g_z_y_y_xxxxyy, g_z_y_y_xxxxyz, g_z_y_y_xxxxz, g_z_y_y_xxxxzz, g_z_y_y_xxxyy, g_z_y_y_xxxyyy, g_z_y_y_xxxyyz, g_z_y_y_xxxyz, g_z_y_y_xxxyzz, g_z_y_y_xxxzz, g_z_y_y_xxxzzz, g_z_y_y_xxyyy, g_z_y_y_xxyyyy, g_z_y_y_xxyyyz, g_z_y_y_xxyyz, g_z_y_y_xxyyzz, g_z_y_y_xxyzz, g_z_y_y_xxyzzz, g_z_y_y_xxzzz, g_z_y_y_xxzzzz, g_z_y_y_xyyyy, g_z_y_y_xyyyyy, g_z_y_y_xyyyyz, g_z_y_y_xyyyz, g_z_y_y_xyyyzz, g_z_y_y_xyyzz, g_z_y_y_xyyzzz, g_z_y_y_xyzzz, g_z_y_y_xyzzzz, g_z_y_y_xzzzz, g_z_y_y_xzzzzz, g_z_y_y_yyyyy, g_z_y_y_yyyyz, g_z_y_y_yyyzz, g_z_y_y_yyzzz, g_z_y_y_yzzzz, g_z_y_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xy_xxxxx[k] = -g_z_y_y_xxxxx[k] * ab_x + g_z_y_y_xxxxxx[k];

                g_z_y_xy_xxxxy[k] = -g_z_y_y_xxxxy[k] * ab_x + g_z_y_y_xxxxxy[k];

                g_z_y_xy_xxxxz[k] = -g_z_y_y_xxxxz[k] * ab_x + g_z_y_y_xxxxxz[k];

                g_z_y_xy_xxxyy[k] = -g_z_y_y_xxxyy[k] * ab_x + g_z_y_y_xxxxyy[k];

                g_z_y_xy_xxxyz[k] = -g_z_y_y_xxxyz[k] * ab_x + g_z_y_y_xxxxyz[k];

                g_z_y_xy_xxxzz[k] = -g_z_y_y_xxxzz[k] * ab_x + g_z_y_y_xxxxzz[k];

                g_z_y_xy_xxyyy[k] = -g_z_y_y_xxyyy[k] * ab_x + g_z_y_y_xxxyyy[k];

                g_z_y_xy_xxyyz[k] = -g_z_y_y_xxyyz[k] * ab_x + g_z_y_y_xxxyyz[k];

                g_z_y_xy_xxyzz[k] = -g_z_y_y_xxyzz[k] * ab_x + g_z_y_y_xxxyzz[k];

                g_z_y_xy_xxzzz[k] = -g_z_y_y_xxzzz[k] * ab_x + g_z_y_y_xxxzzz[k];

                g_z_y_xy_xyyyy[k] = -g_z_y_y_xyyyy[k] * ab_x + g_z_y_y_xxyyyy[k];

                g_z_y_xy_xyyyz[k] = -g_z_y_y_xyyyz[k] * ab_x + g_z_y_y_xxyyyz[k];

                g_z_y_xy_xyyzz[k] = -g_z_y_y_xyyzz[k] * ab_x + g_z_y_y_xxyyzz[k];

                g_z_y_xy_xyzzz[k] = -g_z_y_y_xyzzz[k] * ab_x + g_z_y_y_xxyzzz[k];

                g_z_y_xy_xzzzz[k] = -g_z_y_y_xzzzz[k] * ab_x + g_z_y_y_xxzzzz[k];

                g_z_y_xy_yyyyy[k] = -g_z_y_y_yyyyy[k] * ab_x + g_z_y_y_xyyyyy[k];

                g_z_y_xy_yyyyz[k] = -g_z_y_y_yyyyz[k] * ab_x + g_z_y_y_xyyyyz[k];

                g_z_y_xy_yyyzz[k] = -g_z_y_y_yyyzz[k] * ab_x + g_z_y_y_xyyyzz[k];

                g_z_y_xy_yyzzz[k] = -g_z_y_y_yyzzz[k] * ab_x + g_z_y_y_xyyzzz[k];

                g_z_y_xy_yzzzz[k] = -g_z_y_y_yzzzz[k] * ab_x + g_z_y_y_xyzzzz[k];

                g_z_y_xy_zzzzz[k] = -g_z_y_y_zzzzz[k] * ab_x + g_z_y_y_xzzzzz[k];
            }

            /// Set up 924-945 components of targeted buffer : cbuffer.data(

            auto g_z_y_xz_xxxxx = cbuffer.data(dh_geom_11_off + 924 * ccomps * dcomps);

            auto g_z_y_xz_xxxxy = cbuffer.data(dh_geom_11_off + 925 * ccomps * dcomps);

            auto g_z_y_xz_xxxxz = cbuffer.data(dh_geom_11_off + 926 * ccomps * dcomps);

            auto g_z_y_xz_xxxyy = cbuffer.data(dh_geom_11_off + 927 * ccomps * dcomps);

            auto g_z_y_xz_xxxyz = cbuffer.data(dh_geom_11_off + 928 * ccomps * dcomps);

            auto g_z_y_xz_xxxzz = cbuffer.data(dh_geom_11_off + 929 * ccomps * dcomps);

            auto g_z_y_xz_xxyyy = cbuffer.data(dh_geom_11_off + 930 * ccomps * dcomps);

            auto g_z_y_xz_xxyyz = cbuffer.data(dh_geom_11_off + 931 * ccomps * dcomps);

            auto g_z_y_xz_xxyzz = cbuffer.data(dh_geom_11_off + 932 * ccomps * dcomps);

            auto g_z_y_xz_xxzzz = cbuffer.data(dh_geom_11_off + 933 * ccomps * dcomps);

            auto g_z_y_xz_xyyyy = cbuffer.data(dh_geom_11_off + 934 * ccomps * dcomps);

            auto g_z_y_xz_xyyyz = cbuffer.data(dh_geom_11_off + 935 * ccomps * dcomps);

            auto g_z_y_xz_xyyzz = cbuffer.data(dh_geom_11_off + 936 * ccomps * dcomps);

            auto g_z_y_xz_xyzzz = cbuffer.data(dh_geom_11_off + 937 * ccomps * dcomps);

            auto g_z_y_xz_xzzzz = cbuffer.data(dh_geom_11_off + 938 * ccomps * dcomps);

            auto g_z_y_xz_yyyyy = cbuffer.data(dh_geom_11_off + 939 * ccomps * dcomps);

            auto g_z_y_xz_yyyyz = cbuffer.data(dh_geom_11_off + 940 * ccomps * dcomps);

            auto g_z_y_xz_yyyzz = cbuffer.data(dh_geom_11_off + 941 * ccomps * dcomps);

            auto g_z_y_xz_yyzzz = cbuffer.data(dh_geom_11_off + 942 * ccomps * dcomps);

            auto g_z_y_xz_yzzzz = cbuffer.data(dh_geom_11_off + 943 * ccomps * dcomps);

            auto g_z_y_xz_zzzzz = cbuffer.data(dh_geom_11_off + 944 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_y_xz_xxxxx, g_z_y_xz_xxxxy, g_z_y_xz_xxxxz, g_z_y_xz_xxxyy, g_z_y_xz_xxxyz, g_z_y_xz_xxxzz, g_z_y_xz_xxyyy, g_z_y_xz_xxyyz, g_z_y_xz_xxyzz, g_z_y_xz_xxzzz, g_z_y_xz_xyyyy, g_z_y_xz_xyyyz, g_z_y_xz_xyyzz, g_z_y_xz_xyzzz, g_z_y_xz_xzzzz, g_z_y_xz_yyyyy, g_z_y_xz_yyyyz, g_z_y_xz_yyyzz, g_z_y_xz_yyzzz, g_z_y_xz_yzzzz, g_z_y_xz_zzzzz, g_z_y_z_xxxxx, g_z_y_z_xxxxxx, g_z_y_z_xxxxxy, g_z_y_z_xxxxxz, g_z_y_z_xxxxy, g_z_y_z_xxxxyy, g_z_y_z_xxxxyz, g_z_y_z_xxxxz, g_z_y_z_xxxxzz, g_z_y_z_xxxyy, g_z_y_z_xxxyyy, g_z_y_z_xxxyyz, g_z_y_z_xxxyz, g_z_y_z_xxxyzz, g_z_y_z_xxxzz, g_z_y_z_xxxzzz, g_z_y_z_xxyyy, g_z_y_z_xxyyyy, g_z_y_z_xxyyyz, g_z_y_z_xxyyz, g_z_y_z_xxyyzz, g_z_y_z_xxyzz, g_z_y_z_xxyzzz, g_z_y_z_xxzzz, g_z_y_z_xxzzzz, g_z_y_z_xyyyy, g_z_y_z_xyyyyy, g_z_y_z_xyyyyz, g_z_y_z_xyyyz, g_z_y_z_xyyyzz, g_z_y_z_xyyzz, g_z_y_z_xyyzzz, g_z_y_z_xyzzz, g_z_y_z_xyzzzz, g_z_y_z_xzzzz, g_z_y_z_xzzzzz, g_z_y_z_yyyyy, g_z_y_z_yyyyz, g_z_y_z_yyyzz, g_z_y_z_yyzzz, g_z_y_z_yzzzz, g_z_y_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_xz_xxxxx[k] = -g_z_y_z_xxxxx[k] * ab_x + g_z_y_z_xxxxxx[k];

                g_z_y_xz_xxxxy[k] = -g_z_y_z_xxxxy[k] * ab_x + g_z_y_z_xxxxxy[k];

                g_z_y_xz_xxxxz[k] = -g_z_y_z_xxxxz[k] * ab_x + g_z_y_z_xxxxxz[k];

                g_z_y_xz_xxxyy[k] = -g_z_y_z_xxxyy[k] * ab_x + g_z_y_z_xxxxyy[k];

                g_z_y_xz_xxxyz[k] = -g_z_y_z_xxxyz[k] * ab_x + g_z_y_z_xxxxyz[k];

                g_z_y_xz_xxxzz[k] = -g_z_y_z_xxxzz[k] * ab_x + g_z_y_z_xxxxzz[k];

                g_z_y_xz_xxyyy[k] = -g_z_y_z_xxyyy[k] * ab_x + g_z_y_z_xxxyyy[k];

                g_z_y_xz_xxyyz[k] = -g_z_y_z_xxyyz[k] * ab_x + g_z_y_z_xxxyyz[k];

                g_z_y_xz_xxyzz[k] = -g_z_y_z_xxyzz[k] * ab_x + g_z_y_z_xxxyzz[k];

                g_z_y_xz_xxzzz[k] = -g_z_y_z_xxzzz[k] * ab_x + g_z_y_z_xxxzzz[k];

                g_z_y_xz_xyyyy[k] = -g_z_y_z_xyyyy[k] * ab_x + g_z_y_z_xxyyyy[k];

                g_z_y_xz_xyyyz[k] = -g_z_y_z_xyyyz[k] * ab_x + g_z_y_z_xxyyyz[k];

                g_z_y_xz_xyyzz[k] = -g_z_y_z_xyyzz[k] * ab_x + g_z_y_z_xxyyzz[k];

                g_z_y_xz_xyzzz[k] = -g_z_y_z_xyzzz[k] * ab_x + g_z_y_z_xxyzzz[k];

                g_z_y_xz_xzzzz[k] = -g_z_y_z_xzzzz[k] * ab_x + g_z_y_z_xxzzzz[k];

                g_z_y_xz_yyyyy[k] = -g_z_y_z_yyyyy[k] * ab_x + g_z_y_z_xyyyyy[k];

                g_z_y_xz_yyyyz[k] = -g_z_y_z_yyyyz[k] * ab_x + g_z_y_z_xyyyyz[k];

                g_z_y_xz_yyyzz[k] = -g_z_y_z_yyyzz[k] * ab_x + g_z_y_z_xyyyzz[k];

                g_z_y_xz_yyzzz[k] = -g_z_y_z_yyzzz[k] * ab_x + g_z_y_z_xyyzzz[k];

                g_z_y_xz_yzzzz[k] = -g_z_y_z_yzzzz[k] * ab_x + g_z_y_z_xyzzzz[k];

                g_z_y_xz_zzzzz[k] = -g_z_y_z_zzzzz[k] * ab_x + g_z_y_z_xzzzzz[k];
            }

            /// Set up 945-966 components of targeted buffer : cbuffer.data(

            auto g_z_y_yy_xxxxx = cbuffer.data(dh_geom_11_off + 945 * ccomps * dcomps);

            auto g_z_y_yy_xxxxy = cbuffer.data(dh_geom_11_off + 946 * ccomps * dcomps);

            auto g_z_y_yy_xxxxz = cbuffer.data(dh_geom_11_off + 947 * ccomps * dcomps);

            auto g_z_y_yy_xxxyy = cbuffer.data(dh_geom_11_off + 948 * ccomps * dcomps);

            auto g_z_y_yy_xxxyz = cbuffer.data(dh_geom_11_off + 949 * ccomps * dcomps);

            auto g_z_y_yy_xxxzz = cbuffer.data(dh_geom_11_off + 950 * ccomps * dcomps);

            auto g_z_y_yy_xxyyy = cbuffer.data(dh_geom_11_off + 951 * ccomps * dcomps);

            auto g_z_y_yy_xxyyz = cbuffer.data(dh_geom_11_off + 952 * ccomps * dcomps);

            auto g_z_y_yy_xxyzz = cbuffer.data(dh_geom_11_off + 953 * ccomps * dcomps);

            auto g_z_y_yy_xxzzz = cbuffer.data(dh_geom_11_off + 954 * ccomps * dcomps);

            auto g_z_y_yy_xyyyy = cbuffer.data(dh_geom_11_off + 955 * ccomps * dcomps);

            auto g_z_y_yy_xyyyz = cbuffer.data(dh_geom_11_off + 956 * ccomps * dcomps);

            auto g_z_y_yy_xyyzz = cbuffer.data(dh_geom_11_off + 957 * ccomps * dcomps);

            auto g_z_y_yy_xyzzz = cbuffer.data(dh_geom_11_off + 958 * ccomps * dcomps);

            auto g_z_y_yy_xzzzz = cbuffer.data(dh_geom_11_off + 959 * ccomps * dcomps);

            auto g_z_y_yy_yyyyy = cbuffer.data(dh_geom_11_off + 960 * ccomps * dcomps);

            auto g_z_y_yy_yyyyz = cbuffer.data(dh_geom_11_off + 961 * ccomps * dcomps);

            auto g_z_y_yy_yyyzz = cbuffer.data(dh_geom_11_off + 962 * ccomps * dcomps);

            auto g_z_y_yy_yyzzz = cbuffer.data(dh_geom_11_off + 963 * ccomps * dcomps);

            auto g_z_y_yy_yzzzz = cbuffer.data(dh_geom_11_off + 964 * ccomps * dcomps);

            auto g_z_y_yy_zzzzz = cbuffer.data(dh_geom_11_off + 965 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_y_xxxxx, g_z_0_y_xxxxy, g_z_0_y_xxxxz, g_z_0_y_xxxyy, g_z_0_y_xxxyz, g_z_0_y_xxxzz, g_z_0_y_xxyyy, g_z_0_y_xxyyz, g_z_0_y_xxyzz, g_z_0_y_xxzzz, g_z_0_y_xyyyy, g_z_0_y_xyyyz, g_z_0_y_xyyzz, g_z_0_y_xyzzz, g_z_0_y_xzzzz, g_z_0_y_yyyyy, g_z_0_y_yyyyz, g_z_0_y_yyyzz, g_z_0_y_yyzzz, g_z_0_y_yzzzz, g_z_0_y_zzzzz, g_z_y_y_xxxxx, g_z_y_y_xxxxxy, g_z_y_y_xxxxy, g_z_y_y_xxxxyy, g_z_y_y_xxxxyz, g_z_y_y_xxxxz, g_z_y_y_xxxyy, g_z_y_y_xxxyyy, g_z_y_y_xxxyyz, g_z_y_y_xxxyz, g_z_y_y_xxxyzz, g_z_y_y_xxxzz, g_z_y_y_xxyyy, g_z_y_y_xxyyyy, g_z_y_y_xxyyyz, g_z_y_y_xxyyz, g_z_y_y_xxyyzz, g_z_y_y_xxyzz, g_z_y_y_xxyzzz, g_z_y_y_xxzzz, g_z_y_y_xyyyy, g_z_y_y_xyyyyy, g_z_y_y_xyyyyz, g_z_y_y_xyyyz, g_z_y_y_xyyyzz, g_z_y_y_xyyzz, g_z_y_y_xyyzzz, g_z_y_y_xyzzz, g_z_y_y_xyzzzz, g_z_y_y_xzzzz, g_z_y_y_yyyyy, g_z_y_y_yyyyyy, g_z_y_y_yyyyyz, g_z_y_y_yyyyz, g_z_y_y_yyyyzz, g_z_y_y_yyyzz, g_z_y_y_yyyzzz, g_z_y_y_yyzzz, g_z_y_y_yyzzzz, g_z_y_y_yzzzz, g_z_y_y_yzzzzz, g_z_y_y_zzzzz, g_z_y_yy_xxxxx, g_z_y_yy_xxxxy, g_z_y_yy_xxxxz, g_z_y_yy_xxxyy, g_z_y_yy_xxxyz, g_z_y_yy_xxxzz, g_z_y_yy_xxyyy, g_z_y_yy_xxyyz, g_z_y_yy_xxyzz, g_z_y_yy_xxzzz, g_z_y_yy_xyyyy, g_z_y_yy_xyyyz, g_z_y_yy_xyyzz, g_z_y_yy_xyzzz, g_z_y_yy_xzzzz, g_z_y_yy_yyyyy, g_z_y_yy_yyyyz, g_z_y_yy_yyyzz, g_z_y_yy_yyzzz, g_z_y_yy_yzzzz, g_z_y_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yy_xxxxx[k] = g_z_0_y_xxxxx[k] - g_z_y_y_xxxxx[k] * ab_y + g_z_y_y_xxxxxy[k];

                g_z_y_yy_xxxxy[k] = g_z_0_y_xxxxy[k] - g_z_y_y_xxxxy[k] * ab_y + g_z_y_y_xxxxyy[k];

                g_z_y_yy_xxxxz[k] = g_z_0_y_xxxxz[k] - g_z_y_y_xxxxz[k] * ab_y + g_z_y_y_xxxxyz[k];

                g_z_y_yy_xxxyy[k] = g_z_0_y_xxxyy[k] - g_z_y_y_xxxyy[k] * ab_y + g_z_y_y_xxxyyy[k];

                g_z_y_yy_xxxyz[k] = g_z_0_y_xxxyz[k] - g_z_y_y_xxxyz[k] * ab_y + g_z_y_y_xxxyyz[k];

                g_z_y_yy_xxxzz[k] = g_z_0_y_xxxzz[k] - g_z_y_y_xxxzz[k] * ab_y + g_z_y_y_xxxyzz[k];

                g_z_y_yy_xxyyy[k] = g_z_0_y_xxyyy[k] - g_z_y_y_xxyyy[k] * ab_y + g_z_y_y_xxyyyy[k];

                g_z_y_yy_xxyyz[k] = g_z_0_y_xxyyz[k] - g_z_y_y_xxyyz[k] * ab_y + g_z_y_y_xxyyyz[k];

                g_z_y_yy_xxyzz[k] = g_z_0_y_xxyzz[k] - g_z_y_y_xxyzz[k] * ab_y + g_z_y_y_xxyyzz[k];

                g_z_y_yy_xxzzz[k] = g_z_0_y_xxzzz[k] - g_z_y_y_xxzzz[k] * ab_y + g_z_y_y_xxyzzz[k];

                g_z_y_yy_xyyyy[k] = g_z_0_y_xyyyy[k] - g_z_y_y_xyyyy[k] * ab_y + g_z_y_y_xyyyyy[k];

                g_z_y_yy_xyyyz[k] = g_z_0_y_xyyyz[k] - g_z_y_y_xyyyz[k] * ab_y + g_z_y_y_xyyyyz[k];

                g_z_y_yy_xyyzz[k] = g_z_0_y_xyyzz[k] - g_z_y_y_xyyzz[k] * ab_y + g_z_y_y_xyyyzz[k];

                g_z_y_yy_xyzzz[k] = g_z_0_y_xyzzz[k] - g_z_y_y_xyzzz[k] * ab_y + g_z_y_y_xyyzzz[k];

                g_z_y_yy_xzzzz[k] = g_z_0_y_xzzzz[k] - g_z_y_y_xzzzz[k] * ab_y + g_z_y_y_xyzzzz[k];

                g_z_y_yy_yyyyy[k] = g_z_0_y_yyyyy[k] - g_z_y_y_yyyyy[k] * ab_y + g_z_y_y_yyyyyy[k];

                g_z_y_yy_yyyyz[k] = g_z_0_y_yyyyz[k] - g_z_y_y_yyyyz[k] * ab_y + g_z_y_y_yyyyyz[k];

                g_z_y_yy_yyyzz[k] = g_z_0_y_yyyzz[k] - g_z_y_y_yyyzz[k] * ab_y + g_z_y_y_yyyyzz[k];

                g_z_y_yy_yyzzz[k] = g_z_0_y_yyzzz[k] - g_z_y_y_yyzzz[k] * ab_y + g_z_y_y_yyyzzz[k];

                g_z_y_yy_yzzzz[k] = g_z_0_y_yzzzz[k] - g_z_y_y_yzzzz[k] * ab_y + g_z_y_y_yyzzzz[k];

                g_z_y_yy_zzzzz[k] = g_z_0_y_zzzzz[k] - g_z_y_y_zzzzz[k] * ab_y + g_z_y_y_yzzzzz[k];
            }

            /// Set up 966-987 components of targeted buffer : cbuffer.data(

            auto g_z_y_yz_xxxxx = cbuffer.data(dh_geom_11_off + 966 * ccomps * dcomps);

            auto g_z_y_yz_xxxxy = cbuffer.data(dh_geom_11_off + 967 * ccomps * dcomps);

            auto g_z_y_yz_xxxxz = cbuffer.data(dh_geom_11_off + 968 * ccomps * dcomps);

            auto g_z_y_yz_xxxyy = cbuffer.data(dh_geom_11_off + 969 * ccomps * dcomps);

            auto g_z_y_yz_xxxyz = cbuffer.data(dh_geom_11_off + 970 * ccomps * dcomps);

            auto g_z_y_yz_xxxzz = cbuffer.data(dh_geom_11_off + 971 * ccomps * dcomps);

            auto g_z_y_yz_xxyyy = cbuffer.data(dh_geom_11_off + 972 * ccomps * dcomps);

            auto g_z_y_yz_xxyyz = cbuffer.data(dh_geom_11_off + 973 * ccomps * dcomps);

            auto g_z_y_yz_xxyzz = cbuffer.data(dh_geom_11_off + 974 * ccomps * dcomps);

            auto g_z_y_yz_xxzzz = cbuffer.data(dh_geom_11_off + 975 * ccomps * dcomps);

            auto g_z_y_yz_xyyyy = cbuffer.data(dh_geom_11_off + 976 * ccomps * dcomps);

            auto g_z_y_yz_xyyyz = cbuffer.data(dh_geom_11_off + 977 * ccomps * dcomps);

            auto g_z_y_yz_xyyzz = cbuffer.data(dh_geom_11_off + 978 * ccomps * dcomps);

            auto g_z_y_yz_xyzzz = cbuffer.data(dh_geom_11_off + 979 * ccomps * dcomps);

            auto g_z_y_yz_xzzzz = cbuffer.data(dh_geom_11_off + 980 * ccomps * dcomps);

            auto g_z_y_yz_yyyyy = cbuffer.data(dh_geom_11_off + 981 * ccomps * dcomps);

            auto g_z_y_yz_yyyyz = cbuffer.data(dh_geom_11_off + 982 * ccomps * dcomps);

            auto g_z_y_yz_yyyzz = cbuffer.data(dh_geom_11_off + 983 * ccomps * dcomps);

            auto g_z_y_yz_yyzzz = cbuffer.data(dh_geom_11_off + 984 * ccomps * dcomps);

            auto g_z_y_yz_yzzzz = cbuffer.data(dh_geom_11_off + 985 * ccomps * dcomps);

            auto g_z_y_yz_zzzzz = cbuffer.data(dh_geom_11_off + 986 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_z_xxxxx, g_z_0_z_xxxxy, g_z_0_z_xxxxz, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxzz, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyzz, g_z_0_z_xxzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyzz, g_z_0_z_xyzzz, g_z_0_z_xzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyzz, g_z_0_z_yyzzz, g_z_0_z_yzzzz, g_z_0_z_zzzzz, g_z_y_yz_xxxxx, g_z_y_yz_xxxxy, g_z_y_yz_xxxxz, g_z_y_yz_xxxyy, g_z_y_yz_xxxyz, g_z_y_yz_xxxzz, g_z_y_yz_xxyyy, g_z_y_yz_xxyyz, g_z_y_yz_xxyzz, g_z_y_yz_xxzzz, g_z_y_yz_xyyyy, g_z_y_yz_xyyyz, g_z_y_yz_xyyzz, g_z_y_yz_xyzzz, g_z_y_yz_xzzzz, g_z_y_yz_yyyyy, g_z_y_yz_yyyyz, g_z_y_yz_yyyzz, g_z_y_yz_yyzzz, g_z_y_yz_yzzzz, g_z_y_yz_zzzzz, g_z_y_z_xxxxx, g_z_y_z_xxxxxy, g_z_y_z_xxxxy, g_z_y_z_xxxxyy, g_z_y_z_xxxxyz, g_z_y_z_xxxxz, g_z_y_z_xxxyy, g_z_y_z_xxxyyy, g_z_y_z_xxxyyz, g_z_y_z_xxxyz, g_z_y_z_xxxyzz, g_z_y_z_xxxzz, g_z_y_z_xxyyy, g_z_y_z_xxyyyy, g_z_y_z_xxyyyz, g_z_y_z_xxyyz, g_z_y_z_xxyyzz, g_z_y_z_xxyzz, g_z_y_z_xxyzzz, g_z_y_z_xxzzz, g_z_y_z_xyyyy, g_z_y_z_xyyyyy, g_z_y_z_xyyyyz, g_z_y_z_xyyyz, g_z_y_z_xyyyzz, g_z_y_z_xyyzz, g_z_y_z_xyyzzz, g_z_y_z_xyzzz, g_z_y_z_xyzzzz, g_z_y_z_xzzzz, g_z_y_z_yyyyy, g_z_y_z_yyyyyy, g_z_y_z_yyyyyz, g_z_y_z_yyyyz, g_z_y_z_yyyyzz, g_z_y_z_yyyzz, g_z_y_z_yyyzzz, g_z_y_z_yyzzz, g_z_y_z_yyzzzz, g_z_y_z_yzzzz, g_z_y_z_yzzzzz, g_z_y_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_yz_xxxxx[k] = g_z_0_z_xxxxx[k] - g_z_y_z_xxxxx[k] * ab_y + g_z_y_z_xxxxxy[k];

                g_z_y_yz_xxxxy[k] = g_z_0_z_xxxxy[k] - g_z_y_z_xxxxy[k] * ab_y + g_z_y_z_xxxxyy[k];

                g_z_y_yz_xxxxz[k] = g_z_0_z_xxxxz[k] - g_z_y_z_xxxxz[k] * ab_y + g_z_y_z_xxxxyz[k];

                g_z_y_yz_xxxyy[k] = g_z_0_z_xxxyy[k] - g_z_y_z_xxxyy[k] * ab_y + g_z_y_z_xxxyyy[k];

                g_z_y_yz_xxxyz[k] = g_z_0_z_xxxyz[k] - g_z_y_z_xxxyz[k] * ab_y + g_z_y_z_xxxyyz[k];

                g_z_y_yz_xxxzz[k] = g_z_0_z_xxxzz[k] - g_z_y_z_xxxzz[k] * ab_y + g_z_y_z_xxxyzz[k];

                g_z_y_yz_xxyyy[k] = g_z_0_z_xxyyy[k] - g_z_y_z_xxyyy[k] * ab_y + g_z_y_z_xxyyyy[k];

                g_z_y_yz_xxyyz[k] = g_z_0_z_xxyyz[k] - g_z_y_z_xxyyz[k] * ab_y + g_z_y_z_xxyyyz[k];

                g_z_y_yz_xxyzz[k] = g_z_0_z_xxyzz[k] - g_z_y_z_xxyzz[k] * ab_y + g_z_y_z_xxyyzz[k];

                g_z_y_yz_xxzzz[k] = g_z_0_z_xxzzz[k] - g_z_y_z_xxzzz[k] * ab_y + g_z_y_z_xxyzzz[k];

                g_z_y_yz_xyyyy[k] = g_z_0_z_xyyyy[k] - g_z_y_z_xyyyy[k] * ab_y + g_z_y_z_xyyyyy[k];

                g_z_y_yz_xyyyz[k] = g_z_0_z_xyyyz[k] - g_z_y_z_xyyyz[k] * ab_y + g_z_y_z_xyyyyz[k];

                g_z_y_yz_xyyzz[k] = g_z_0_z_xyyzz[k] - g_z_y_z_xyyzz[k] * ab_y + g_z_y_z_xyyyzz[k];

                g_z_y_yz_xyzzz[k] = g_z_0_z_xyzzz[k] - g_z_y_z_xyzzz[k] * ab_y + g_z_y_z_xyyzzz[k];

                g_z_y_yz_xzzzz[k] = g_z_0_z_xzzzz[k] - g_z_y_z_xzzzz[k] * ab_y + g_z_y_z_xyzzzz[k];

                g_z_y_yz_yyyyy[k] = g_z_0_z_yyyyy[k] - g_z_y_z_yyyyy[k] * ab_y + g_z_y_z_yyyyyy[k];

                g_z_y_yz_yyyyz[k] = g_z_0_z_yyyyz[k] - g_z_y_z_yyyyz[k] * ab_y + g_z_y_z_yyyyyz[k];

                g_z_y_yz_yyyzz[k] = g_z_0_z_yyyzz[k] - g_z_y_z_yyyzz[k] * ab_y + g_z_y_z_yyyyzz[k];

                g_z_y_yz_yyzzz[k] = g_z_0_z_yyzzz[k] - g_z_y_z_yyzzz[k] * ab_y + g_z_y_z_yyyzzz[k];

                g_z_y_yz_yzzzz[k] = g_z_0_z_yzzzz[k] - g_z_y_z_yzzzz[k] * ab_y + g_z_y_z_yyzzzz[k];

                g_z_y_yz_zzzzz[k] = g_z_0_z_zzzzz[k] - g_z_y_z_zzzzz[k] * ab_y + g_z_y_z_yzzzzz[k];
            }

            /// Set up 987-1008 components of targeted buffer : cbuffer.data(

            auto g_z_y_zz_xxxxx = cbuffer.data(dh_geom_11_off + 987 * ccomps * dcomps);

            auto g_z_y_zz_xxxxy = cbuffer.data(dh_geom_11_off + 988 * ccomps * dcomps);

            auto g_z_y_zz_xxxxz = cbuffer.data(dh_geom_11_off + 989 * ccomps * dcomps);

            auto g_z_y_zz_xxxyy = cbuffer.data(dh_geom_11_off + 990 * ccomps * dcomps);

            auto g_z_y_zz_xxxyz = cbuffer.data(dh_geom_11_off + 991 * ccomps * dcomps);

            auto g_z_y_zz_xxxzz = cbuffer.data(dh_geom_11_off + 992 * ccomps * dcomps);

            auto g_z_y_zz_xxyyy = cbuffer.data(dh_geom_11_off + 993 * ccomps * dcomps);

            auto g_z_y_zz_xxyyz = cbuffer.data(dh_geom_11_off + 994 * ccomps * dcomps);

            auto g_z_y_zz_xxyzz = cbuffer.data(dh_geom_11_off + 995 * ccomps * dcomps);

            auto g_z_y_zz_xxzzz = cbuffer.data(dh_geom_11_off + 996 * ccomps * dcomps);

            auto g_z_y_zz_xyyyy = cbuffer.data(dh_geom_11_off + 997 * ccomps * dcomps);

            auto g_z_y_zz_xyyyz = cbuffer.data(dh_geom_11_off + 998 * ccomps * dcomps);

            auto g_z_y_zz_xyyzz = cbuffer.data(dh_geom_11_off + 999 * ccomps * dcomps);

            auto g_z_y_zz_xyzzz = cbuffer.data(dh_geom_11_off + 1000 * ccomps * dcomps);

            auto g_z_y_zz_xzzzz = cbuffer.data(dh_geom_11_off + 1001 * ccomps * dcomps);

            auto g_z_y_zz_yyyyy = cbuffer.data(dh_geom_11_off + 1002 * ccomps * dcomps);

            auto g_z_y_zz_yyyyz = cbuffer.data(dh_geom_11_off + 1003 * ccomps * dcomps);

            auto g_z_y_zz_yyyzz = cbuffer.data(dh_geom_11_off + 1004 * ccomps * dcomps);

            auto g_z_y_zz_yyzzz = cbuffer.data(dh_geom_11_off + 1005 * ccomps * dcomps);

            auto g_z_y_zz_yzzzz = cbuffer.data(dh_geom_11_off + 1006 * ccomps * dcomps);

            auto g_z_y_zz_zzzzz = cbuffer.data(dh_geom_11_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_z_xxxxx, g_0_y_z_xxxxy, g_0_y_z_xxxxz, g_0_y_z_xxxyy, g_0_y_z_xxxyz, g_0_y_z_xxxzz, g_0_y_z_xxyyy, g_0_y_z_xxyyz, g_0_y_z_xxyzz, g_0_y_z_xxzzz, g_0_y_z_xyyyy, g_0_y_z_xyyyz, g_0_y_z_xyyzz, g_0_y_z_xyzzz, g_0_y_z_xzzzz, g_0_y_z_yyyyy, g_0_y_z_yyyyz, g_0_y_z_yyyzz, g_0_y_z_yyzzz, g_0_y_z_yzzzz, g_0_y_z_zzzzz, g_z_y_z_xxxxx, g_z_y_z_xxxxxz, g_z_y_z_xxxxy, g_z_y_z_xxxxyz, g_z_y_z_xxxxz, g_z_y_z_xxxxzz, g_z_y_z_xxxyy, g_z_y_z_xxxyyz, g_z_y_z_xxxyz, g_z_y_z_xxxyzz, g_z_y_z_xxxzz, g_z_y_z_xxxzzz, g_z_y_z_xxyyy, g_z_y_z_xxyyyz, g_z_y_z_xxyyz, g_z_y_z_xxyyzz, g_z_y_z_xxyzz, g_z_y_z_xxyzzz, g_z_y_z_xxzzz, g_z_y_z_xxzzzz, g_z_y_z_xyyyy, g_z_y_z_xyyyyz, g_z_y_z_xyyyz, g_z_y_z_xyyyzz, g_z_y_z_xyyzz, g_z_y_z_xyyzzz, g_z_y_z_xyzzz, g_z_y_z_xyzzzz, g_z_y_z_xzzzz, g_z_y_z_xzzzzz, g_z_y_z_yyyyy, g_z_y_z_yyyyyz, g_z_y_z_yyyyz, g_z_y_z_yyyyzz, g_z_y_z_yyyzz, g_z_y_z_yyyzzz, g_z_y_z_yyzzz, g_z_y_z_yyzzzz, g_z_y_z_yzzzz, g_z_y_z_yzzzzz, g_z_y_z_zzzzz, g_z_y_z_zzzzzz, g_z_y_zz_xxxxx, g_z_y_zz_xxxxy, g_z_y_zz_xxxxz, g_z_y_zz_xxxyy, g_z_y_zz_xxxyz, g_z_y_zz_xxxzz, g_z_y_zz_xxyyy, g_z_y_zz_xxyyz, g_z_y_zz_xxyzz, g_z_y_zz_xxzzz, g_z_y_zz_xyyyy, g_z_y_zz_xyyyz, g_z_y_zz_xyyzz, g_z_y_zz_xyzzz, g_z_y_zz_xzzzz, g_z_y_zz_yyyyy, g_z_y_zz_yyyyz, g_z_y_zz_yyyzz, g_z_y_zz_yyzzz, g_z_y_zz_yzzzz, g_z_y_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_y_zz_xxxxx[k] = -g_0_y_z_xxxxx[k] - g_z_y_z_xxxxx[k] * ab_z + g_z_y_z_xxxxxz[k];

                g_z_y_zz_xxxxy[k] = -g_0_y_z_xxxxy[k] - g_z_y_z_xxxxy[k] * ab_z + g_z_y_z_xxxxyz[k];

                g_z_y_zz_xxxxz[k] = -g_0_y_z_xxxxz[k] - g_z_y_z_xxxxz[k] * ab_z + g_z_y_z_xxxxzz[k];

                g_z_y_zz_xxxyy[k] = -g_0_y_z_xxxyy[k] - g_z_y_z_xxxyy[k] * ab_z + g_z_y_z_xxxyyz[k];

                g_z_y_zz_xxxyz[k] = -g_0_y_z_xxxyz[k] - g_z_y_z_xxxyz[k] * ab_z + g_z_y_z_xxxyzz[k];

                g_z_y_zz_xxxzz[k] = -g_0_y_z_xxxzz[k] - g_z_y_z_xxxzz[k] * ab_z + g_z_y_z_xxxzzz[k];

                g_z_y_zz_xxyyy[k] = -g_0_y_z_xxyyy[k] - g_z_y_z_xxyyy[k] * ab_z + g_z_y_z_xxyyyz[k];

                g_z_y_zz_xxyyz[k] = -g_0_y_z_xxyyz[k] - g_z_y_z_xxyyz[k] * ab_z + g_z_y_z_xxyyzz[k];

                g_z_y_zz_xxyzz[k] = -g_0_y_z_xxyzz[k] - g_z_y_z_xxyzz[k] * ab_z + g_z_y_z_xxyzzz[k];

                g_z_y_zz_xxzzz[k] = -g_0_y_z_xxzzz[k] - g_z_y_z_xxzzz[k] * ab_z + g_z_y_z_xxzzzz[k];

                g_z_y_zz_xyyyy[k] = -g_0_y_z_xyyyy[k] - g_z_y_z_xyyyy[k] * ab_z + g_z_y_z_xyyyyz[k];

                g_z_y_zz_xyyyz[k] = -g_0_y_z_xyyyz[k] - g_z_y_z_xyyyz[k] * ab_z + g_z_y_z_xyyyzz[k];

                g_z_y_zz_xyyzz[k] = -g_0_y_z_xyyzz[k] - g_z_y_z_xyyzz[k] * ab_z + g_z_y_z_xyyzzz[k];

                g_z_y_zz_xyzzz[k] = -g_0_y_z_xyzzz[k] - g_z_y_z_xyzzz[k] * ab_z + g_z_y_z_xyzzzz[k];

                g_z_y_zz_xzzzz[k] = -g_0_y_z_xzzzz[k] - g_z_y_z_xzzzz[k] * ab_z + g_z_y_z_xzzzzz[k];

                g_z_y_zz_yyyyy[k] = -g_0_y_z_yyyyy[k] - g_z_y_z_yyyyy[k] * ab_z + g_z_y_z_yyyyyz[k];

                g_z_y_zz_yyyyz[k] = -g_0_y_z_yyyyz[k] - g_z_y_z_yyyyz[k] * ab_z + g_z_y_z_yyyyzz[k];

                g_z_y_zz_yyyzz[k] = -g_0_y_z_yyyzz[k] - g_z_y_z_yyyzz[k] * ab_z + g_z_y_z_yyyzzz[k];

                g_z_y_zz_yyzzz[k] = -g_0_y_z_yyzzz[k] - g_z_y_z_yyzzz[k] * ab_z + g_z_y_z_yyzzzz[k];

                g_z_y_zz_yzzzz[k] = -g_0_y_z_yzzzz[k] - g_z_y_z_yzzzz[k] * ab_z + g_z_y_z_yzzzzz[k];

                g_z_y_zz_zzzzz[k] = -g_0_y_z_zzzzz[k] - g_z_y_z_zzzzz[k] * ab_z + g_z_y_z_zzzzzz[k];
            }

            /// Set up 1008-1029 components of targeted buffer : cbuffer.data(

            auto g_z_z_xx_xxxxx = cbuffer.data(dh_geom_11_off + 1008 * ccomps * dcomps);

            auto g_z_z_xx_xxxxy = cbuffer.data(dh_geom_11_off + 1009 * ccomps * dcomps);

            auto g_z_z_xx_xxxxz = cbuffer.data(dh_geom_11_off + 1010 * ccomps * dcomps);

            auto g_z_z_xx_xxxyy = cbuffer.data(dh_geom_11_off + 1011 * ccomps * dcomps);

            auto g_z_z_xx_xxxyz = cbuffer.data(dh_geom_11_off + 1012 * ccomps * dcomps);

            auto g_z_z_xx_xxxzz = cbuffer.data(dh_geom_11_off + 1013 * ccomps * dcomps);

            auto g_z_z_xx_xxyyy = cbuffer.data(dh_geom_11_off + 1014 * ccomps * dcomps);

            auto g_z_z_xx_xxyyz = cbuffer.data(dh_geom_11_off + 1015 * ccomps * dcomps);

            auto g_z_z_xx_xxyzz = cbuffer.data(dh_geom_11_off + 1016 * ccomps * dcomps);

            auto g_z_z_xx_xxzzz = cbuffer.data(dh_geom_11_off + 1017 * ccomps * dcomps);

            auto g_z_z_xx_xyyyy = cbuffer.data(dh_geom_11_off + 1018 * ccomps * dcomps);

            auto g_z_z_xx_xyyyz = cbuffer.data(dh_geom_11_off + 1019 * ccomps * dcomps);

            auto g_z_z_xx_xyyzz = cbuffer.data(dh_geom_11_off + 1020 * ccomps * dcomps);

            auto g_z_z_xx_xyzzz = cbuffer.data(dh_geom_11_off + 1021 * ccomps * dcomps);

            auto g_z_z_xx_xzzzz = cbuffer.data(dh_geom_11_off + 1022 * ccomps * dcomps);

            auto g_z_z_xx_yyyyy = cbuffer.data(dh_geom_11_off + 1023 * ccomps * dcomps);

            auto g_z_z_xx_yyyyz = cbuffer.data(dh_geom_11_off + 1024 * ccomps * dcomps);

            auto g_z_z_xx_yyyzz = cbuffer.data(dh_geom_11_off + 1025 * ccomps * dcomps);

            auto g_z_z_xx_yyzzz = cbuffer.data(dh_geom_11_off + 1026 * ccomps * dcomps);

            auto g_z_z_xx_yzzzz = cbuffer.data(dh_geom_11_off + 1027 * ccomps * dcomps);

            auto g_z_z_xx_zzzzz = cbuffer.data(dh_geom_11_off + 1028 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_x_xxxxx, g_z_z_x_xxxxxx, g_z_z_x_xxxxxy, g_z_z_x_xxxxxz, g_z_z_x_xxxxy, g_z_z_x_xxxxyy, g_z_z_x_xxxxyz, g_z_z_x_xxxxz, g_z_z_x_xxxxzz, g_z_z_x_xxxyy, g_z_z_x_xxxyyy, g_z_z_x_xxxyyz, g_z_z_x_xxxyz, g_z_z_x_xxxyzz, g_z_z_x_xxxzz, g_z_z_x_xxxzzz, g_z_z_x_xxyyy, g_z_z_x_xxyyyy, g_z_z_x_xxyyyz, g_z_z_x_xxyyz, g_z_z_x_xxyyzz, g_z_z_x_xxyzz, g_z_z_x_xxyzzz, g_z_z_x_xxzzz, g_z_z_x_xxzzzz, g_z_z_x_xyyyy, g_z_z_x_xyyyyy, g_z_z_x_xyyyyz, g_z_z_x_xyyyz, g_z_z_x_xyyyzz, g_z_z_x_xyyzz, g_z_z_x_xyyzzz, g_z_z_x_xyzzz, g_z_z_x_xyzzzz, g_z_z_x_xzzzz, g_z_z_x_xzzzzz, g_z_z_x_yyyyy, g_z_z_x_yyyyz, g_z_z_x_yyyzz, g_z_z_x_yyzzz, g_z_z_x_yzzzz, g_z_z_x_zzzzz, g_z_z_xx_xxxxx, g_z_z_xx_xxxxy, g_z_z_xx_xxxxz, g_z_z_xx_xxxyy, g_z_z_xx_xxxyz, g_z_z_xx_xxxzz, g_z_z_xx_xxyyy, g_z_z_xx_xxyyz, g_z_z_xx_xxyzz, g_z_z_xx_xxzzz, g_z_z_xx_xyyyy, g_z_z_xx_xyyyz, g_z_z_xx_xyyzz, g_z_z_xx_xyzzz, g_z_z_xx_xzzzz, g_z_z_xx_yyyyy, g_z_z_xx_yyyyz, g_z_z_xx_yyyzz, g_z_z_xx_yyzzz, g_z_z_xx_yzzzz, g_z_z_xx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xx_xxxxx[k] = -g_z_z_x_xxxxx[k] * ab_x + g_z_z_x_xxxxxx[k];

                g_z_z_xx_xxxxy[k] = -g_z_z_x_xxxxy[k] * ab_x + g_z_z_x_xxxxxy[k];

                g_z_z_xx_xxxxz[k] = -g_z_z_x_xxxxz[k] * ab_x + g_z_z_x_xxxxxz[k];

                g_z_z_xx_xxxyy[k] = -g_z_z_x_xxxyy[k] * ab_x + g_z_z_x_xxxxyy[k];

                g_z_z_xx_xxxyz[k] = -g_z_z_x_xxxyz[k] * ab_x + g_z_z_x_xxxxyz[k];

                g_z_z_xx_xxxzz[k] = -g_z_z_x_xxxzz[k] * ab_x + g_z_z_x_xxxxzz[k];

                g_z_z_xx_xxyyy[k] = -g_z_z_x_xxyyy[k] * ab_x + g_z_z_x_xxxyyy[k];

                g_z_z_xx_xxyyz[k] = -g_z_z_x_xxyyz[k] * ab_x + g_z_z_x_xxxyyz[k];

                g_z_z_xx_xxyzz[k] = -g_z_z_x_xxyzz[k] * ab_x + g_z_z_x_xxxyzz[k];

                g_z_z_xx_xxzzz[k] = -g_z_z_x_xxzzz[k] * ab_x + g_z_z_x_xxxzzz[k];

                g_z_z_xx_xyyyy[k] = -g_z_z_x_xyyyy[k] * ab_x + g_z_z_x_xxyyyy[k];

                g_z_z_xx_xyyyz[k] = -g_z_z_x_xyyyz[k] * ab_x + g_z_z_x_xxyyyz[k];

                g_z_z_xx_xyyzz[k] = -g_z_z_x_xyyzz[k] * ab_x + g_z_z_x_xxyyzz[k];

                g_z_z_xx_xyzzz[k] = -g_z_z_x_xyzzz[k] * ab_x + g_z_z_x_xxyzzz[k];

                g_z_z_xx_xzzzz[k] = -g_z_z_x_xzzzz[k] * ab_x + g_z_z_x_xxzzzz[k];

                g_z_z_xx_yyyyy[k] = -g_z_z_x_yyyyy[k] * ab_x + g_z_z_x_xyyyyy[k];

                g_z_z_xx_yyyyz[k] = -g_z_z_x_yyyyz[k] * ab_x + g_z_z_x_xyyyyz[k];

                g_z_z_xx_yyyzz[k] = -g_z_z_x_yyyzz[k] * ab_x + g_z_z_x_xyyyzz[k];

                g_z_z_xx_yyzzz[k] = -g_z_z_x_yyzzz[k] * ab_x + g_z_z_x_xyyzzz[k];

                g_z_z_xx_yzzzz[k] = -g_z_z_x_yzzzz[k] * ab_x + g_z_z_x_xyzzzz[k];

                g_z_z_xx_zzzzz[k] = -g_z_z_x_zzzzz[k] * ab_x + g_z_z_x_xzzzzz[k];
            }

            /// Set up 1029-1050 components of targeted buffer : cbuffer.data(

            auto g_z_z_xy_xxxxx = cbuffer.data(dh_geom_11_off + 1029 * ccomps * dcomps);

            auto g_z_z_xy_xxxxy = cbuffer.data(dh_geom_11_off + 1030 * ccomps * dcomps);

            auto g_z_z_xy_xxxxz = cbuffer.data(dh_geom_11_off + 1031 * ccomps * dcomps);

            auto g_z_z_xy_xxxyy = cbuffer.data(dh_geom_11_off + 1032 * ccomps * dcomps);

            auto g_z_z_xy_xxxyz = cbuffer.data(dh_geom_11_off + 1033 * ccomps * dcomps);

            auto g_z_z_xy_xxxzz = cbuffer.data(dh_geom_11_off + 1034 * ccomps * dcomps);

            auto g_z_z_xy_xxyyy = cbuffer.data(dh_geom_11_off + 1035 * ccomps * dcomps);

            auto g_z_z_xy_xxyyz = cbuffer.data(dh_geom_11_off + 1036 * ccomps * dcomps);

            auto g_z_z_xy_xxyzz = cbuffer.data(dh_geom_11_off + 1037 * ccomps * dcomps);

            auto g_z_z_xy_xxzzz = cbuffer.data(dh_geom_11_off + 1038 * ccomps * dcomps);

            auto g_z_z_xy_xyyyy = cbuffer.data(dh_geom_11_off + 1039 * ccomps * dcomps);

            auto g_z_z_xy_xyyyz = cbuffer.data(dh_geom_11_off + 1040 * ccomps * dcomps);

            auto g_z_z_xy_xyyzz = cbuffer.data(dh_geom_11_off + 1041 * ccomps * dcomps);

            auto g_z_z_xy_xyzzz = cbuffer.data(dh_geom_11_off + 1042 * ccomps * dcomps);

            auto g_z_z_xy_xzzzz = cbuffer.data(dh_geom_11_off + 1043 * ccomps * dcomps);

            auto g_z_z_xy_yyyyy = cbuffer.data(dh_geom_11_off + 1044 * ccomps * dcomps);

            auto g_z_z_xy_yyyyz = cbuffer.data(dh_geom_11_off + 1045 * ccomps * dcomps);

            auto g_z_z_xy_yyyzz = cbuffer.data(dh_geom_11_off + 1046 * ccomps * dcomps);

            auto g_z_z_xy_yyzzz = cbuffer.data(dh_geom_11_off + 1047 * ccomps * dcomps);

            auto g_z_z_xy_yzzzz = cbuffer.data(dh_geom_11_off + 1048 * ccomps * dcomps);

            auto g_z_z_xy_zzzzz = cbuffer.data(dh_geom_11_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xy_xxxxx, g_z_z_xy_xxxxy, g_z_z_xy_xxxxz, g_z_z_xy_xxxyy, g_z_z_xy_xxxyz, g_z_z_xy_xxxzz, g_z_z_xy_xxyyy, g_z_z_xy_xxyyz, g_z_z_xy_xxyzz, g_z_z_xy_xxzzz, g_z_z_xy_xyyyy, g_z_z_xy_xyyyz, g_z_z_xy_xyyzz, g_z_z_xy_xyzzz, g_z_z_xy_xzzzz, g_z_z_xy_yyyyy, g_z_z_xy_yyyyz, g_z_z_xy_yyyzz, g_z_z_xy_yyzzz, g_z_z_xy_yzzzz, g_z_z_xy_zzzzz, g_z_z_y_xxxxx, g_z_z_y_xxxxxx, g_z_z_y_xxxxxy, g_z_z_y_xxxxxz, g_z_z_y_xxxxy, g_z_z_y_xxxxyy, g_z_z_y_xxxxyz, g_z_z_y_xxxxz, g_z_z_y_xxxxzz, g_z_z_y_xxxyy, g_z_z_y_xxxyyy, g_z_z_y_xxxyyz, g_z_z_y_xxxyz, g_z_z_y_xxxyzz, g_z_z_y_xxxzz, g_z_z_y_xxxzzz, g_z_z_y_xxyyy, g_z_z_y_xxyyyy, g_z_z_y_xxyyyz, g_z_z_y_xxyyz, g_z_z_y_xxyyzz, g_z_z_y_xxyzz, g_z_z_y_xxyzzz, g_z_z_y_xxzzz, g_z_z_y_xxzzzz, g_z_z_y_xyyyy, g_z_z_y_xyyyyy, g_z_z_y_xyyyyz, g_z_z_y_xyyyz, g_z_z_y_xyyyzz, g_z_z_y_xyyzz, g_z_z_y_xyyzzz, g_z_z_y_xyzzz, g_z_z_y_xyzzzz, g_z_z_y_xzzzz, g_z_z_y_xzzzzz, g_z_z_y_yyyyy, g_z_z_y_yyyyz, g_z_z_y_yyyzz, g_z_z_y_yyzzz, g_z_z_y_yzzzz, g_z_z_y_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xy_xxxxx[k] = -g_z_z_y_xxxxx[k] * ab_x + g_z_z_y_xxxxxx[k];

                g_z_z_xy_xxxxy[k] = -g_z_z_y_xxxxy[k] * ab_x + g_z_z_y_xxxxxy[k];

                g_z_z_xy_xxxxz[k] = -g_z_z_y_xxxxz[k] * ab_x + g_z_z_y_xxxxxz[k];

                g_z_z_xy_xxxyy[k] = -g_z_z_y_xxxyy[k] * ab_x + g_z_z_y_xxxxyy[k];

                g_z_z_xy_xxxyz[k] = -g_z_z_y_xxxyz[k] * ab_x + g_z_z_y_xxxxyz[k];

                g_z_z_xy_xxxzz[k] = -g_z_z_y_xxxzz[k] * ab_x + g_z_z_y_xxxxzz[k];

                g_z_z_xy_xxyyy[k] = -g_z_z_y_xxyyy[k] * ab_x + g_z_z_y_xxxyyy[k];

                g_z_z_xy_xxyyz[k] = -g_z_z_y_xxyyz[k] * ab_x + g_z_z_y_xxxyyz[k];

                g_z_z_xy_xxyzz[k] = -g_z_z_y_xxyzz[k] * ab_x + g_z_z_y_xxxyzz[k];

                g_z_z_xy_xxzzz[k] = -g_z_z_y_xxzzz[k] * ab_x + g_z_z_y_xxxzzz[k];

                g_z_z_xy_xyyyy[k] = -g_z_z_y_xyyyy[k] * ab_x + g_z_z_y_xxyyyy[k];

                g_z_z_xy_xyyyz[k] = -g_z_z_y_xyyyz[k] * ab_x + g_z_z_y_xxyyyz[k];

                g_z_z_xy_xyyzz[k] = -g_z_z_y_xyyzz[k] * ab_x + g_z_z_y_xxyyzz[k];

                g_z_z_xy_xyzzz[k] = -g_z_z_y_xyzzz[k] * ab_x + g_z_z_y_xxyzzz[k];

                g_z_z_xy_xzzzz[k] = -g_z_z_y_xzzzz[k] * ab_x + g_z_z_y_xxzzzz[k];

                g_z_z_xy_yyyyy[k] = -g_z_z_y_yyyyy[k] * ab_x + g_z_z_y_xyyyyy[k];

                g_z_z_xy_yyyyz[k] = -g_z_z_y_yyyyz[k] * ab_x + g_z_z_y_xyyyyz[k];

                g_z_z_xy_yyyzz[k] = -g_z_z_y_yyyzz[k] * ab_x + g_z_z_y_xyyyzz[k];

                g_z_z_xy_yyzzz[k] = -g_z_z_y_yyzzz[k] * ab_x + g_z_z_y_xyyzzz[k];

                g_z_z_xy_yzzzz[k] = -g_z_z_y_yzzzz[k] * ab_x + g_z_z_y_xyzzzz[k];

                g_z_z_xy_zzzzz[k] = -g_z_z_y_zzzzz[k] * ab_x + g_z_z_y_xzzzzz[k];
            }

            /// Set up 1050-1071 components of targeted buffer : cbuffer.data(

            auto g_z_z_xz_xxxxx = cbuffer.data(dh_geom_11_off + 1050 * ccomps * dcomps);

            auto g_z_z_xz_xxxxy = cbuffer.data(dh_geom_11_off + 1051 * ccomps * dcomps);

            auto g_z_z_xz_xxxxz = cbuffer.data(dh_geom_11_off + 1052 * ccomps * dcomps);

            auto g_z_z_xz_xxxyy = cbuffer.data(dh_geom_11_off + 1053 * ccomps * dcomps);

            auto g_z_z_xz_xxxyz = cbuffer.data(dh_geom_11_off + 1054 * ccomps * dcomps);

            auto g_z_z_xz_xxxzz = cbuffer.data(dh_geom_11_off + 1055 * ccomps * dcomps);

            auto g_z_z_xz_xxyyy = cbuffer.data(dh_geom_11_off + 1056 * ccomps * dcomps);

            auto g_z_z_xz_xxyyz = cbuffer.data(dh_geom_11_off + 1057 * ccomps * dcomps);

            auto g_z_z_xz_xxyzz = cbuffer.data(dh_geom_11_off + 1058 * ccomps * dcomps);

            auto g_z_z_xz_xxzzz = cbuffer.data(dh_geom_11_off + 1059 * ccomps * dcomps);

            auto g_z_z_xz_xyyyy = cbuffer.data(dh_geom_11_off + 1060 * ccomps * dcomps);

            auto g_z_z_xz_xyyyz = cbuffer.data(dh_geom_11_off + 1061 * ccomps * dcomps);

            auto g_z_z_xz_xyyzz = cbuffer.data(dh_geom_11_off + 1062 * ccomps * dcomps);

            auto g_z_z_xz_xyzzz = cbuffer.data(dh_geom_11_off + 1063 * ccomps * dcomps);

            auto g_z_z_xz_xzzzz = cbuffer.data(dh_geom_11_off + 1064 * ccomps * dcomps);

            auto g_z_z_xz_yyyyy = cbuffer.data(dh_geom_11_off + 1065 * ccomps * dcomps);

            auto g_z_z_xz_yyyyz = cbuffer.data(dh_geom_11_off + 1066 * ccomps * dcomps);

            auto g_z_z_xz_yyyzz = cbuffer.data(dh_geom_11_off + 1067 * ccomps * dcomps);

            auto g_z_z_xz_yyzzz = cbuffer.data(dh_geom_11_off + 1068 * ccomps * dcomps);

            auto g_z_z_xz_yzzzz = cbuffer.data(dh_geom_11_off + 1069 * ccomps * dcomps);

            auto g_z_z_xz_zzzzz = cbuffer.data(dh_geom_11_off + 1070 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_xz_xxxxx, g_z_z_xz_xxxxy, g_z_z_xz_xxxxz, g_z_z_xz_xxxyy, g_z_z_xz_xxxyz, g_z_z_xz_xxxzz, g_z_z_xz_xxyyy, g_z_z_xz_xxyyz, g_z_z_xz_xxyzz, g_z_z_xz_xxzzz, g_z_z_xz_xyyyy, g_z_z_xz_xyyyz, g_z_z_xz_xyyzz, g_z_z_xz_xyzzz, g_z_z_xz_xzzzz, g_z_z_xz_yyyyy, g_z_z_xz_yyyyz, g_z_z_xz_yyyzz, g_z_z_xz_yyzzz, g_z_z_xz_yzzzz, g_z_z_xz_zzzzz, g_z_z_z_xxxxx, g_z_z_z_xxxxxx, g_z_z_z_xxxxxy, g_z_z_z_xxxxxz, g_z_z_z_xxxxy, g_z_z_z_xxxxyy, g_z_z_z_xxxxyz, g_z_z_z_xxxxz, g_z_z_z_xxxxzz, g_z_z_z_xxxyy, g_z_z_z_xxxyyy, g_z_z_z_xxxyyz, g_z_z_z_xxxyz, g_z_z_z_xxxyzz, g_z_z_z_xxxzz, g_z_z_z_xxxzzz, g_z_z_z_xxyyy, g_z_z_z_xxyyyy, g_z_z_z_xxyyyz, g_z_z_z_xxyyz, g_z_z_z_xxyyzz, g_z_z_z_xxyzz, g_z_z_z_xxyzzz, g_z_z_z_xxzzz, g_z_z_z_xxzzzz, g_z_z_z_xyyyy, g_z_z_z_xyyyyy, g_z_z_z_xyyyyz, g_z_z_z_xyyyz, g_z_z_z_xyyyzz, g_z_z_z_xyyzz, g_z_z_z_xyyzzz, g_z_z_z_xyzzz, g_z_z_z_xyzzzz, g_z_z_z_xzzzz, g_z_z_z_xzzzzz, g_z_z_z_yyyyy, g_z_z_z_yyyyz, g_z_z_z_yyyzz, g_z_z_z_yyzzz, g_z_z_z_yzzzz, g_z_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_xz_xxxxx[k] = -g_z_z_z_xxxxx[k] * ab_x + g_z_z_z_xxxxxx[k];

                g_z_z_xz_xxxxy[k] = -g_z_z_z_xxxxy[k] * ab_x + g_z_z_z_xxxxxy[k];

                g_z_z_xz_xxxxz[k] = -g_z_z_z_xxxxz[k] * ab_x + g_z_z_z_xxxxxz[k];

                g_z_z_xz_xxxyy[k] = -g_z_z_z_xxxyy[k] * ab_x + g_z_z_z_xxxxyy[k];

                g_z_z_xz_xxxyz[k] = -g_z_z_z_xxxyz[k] * ab_x + g_z_z_z_xxxxyz[k];

                g_z_z_xz_xxxzz[k] = -g_z_z_z_xxxzz[k] * ab_x + g_z_z_z_xxxxzz[k];

                g_z_z_xz_xxyyy[k] = -g_z_z_z_xxyyy[k] * ab_x + g_z_z_z_xxxyyy[k];

                g_z_z_xz_xxyyz[k] = -g_z_z_z_xxyyz[k] * ab_x + g_z_z_z_xxxyyz[k];

                g_z_z_xz_xxyzz[k] = -g_z_z_z_xxyzz[k] * ab_x + g_z_z_z_xxxyzz[k];

                g_z_z_xz_xxzzz[k] = -g_z_z_z_xxzzz[k] * ab_x + g_z_z_z_xxxzzz[k];

                g_z_z_xz_xyyyy[k] = -g_z_z_z_xyyyy[k] * ab_x + g_z_z_z_xxyyyy[k];

                g_z_z_xz_xyyyz[k] = -g_z_z_z_xyyyz[k] * ab_x + g_z_z_z_xxyyyz[k];

                g_z_z_xz_xyyzz[k] = -g_z_z_z_xyyzz[k] * ab_x + g_z_z_z_xxyyzz[k];

                g_z_z_xz_xyzzz[k] = -g_z_z_z_xyzzz[k] * ab_x + g_z_z_z_xxyzzz[k];

                g_z_z_xz_xzzzz[k] = -g_z_z_z_xzzzz[k] * ab_x + g_z_z_z_xxzzzz[k];

                g_z_z_xz_yyyyy[k] = -g_z_z_z_yyyyy[k] * ab_x + g_z_z_z_xyyyyy[k];

                g_z_z_xz_yyyyz[k] = -g_z_z_z_yyyyz[k] * ab_x + g_z_z_z_xyyyyz[k];

                g_z_z_xz_yyyzz[k] = -g_z_z_z_yyyzz[k] * ab_x + g_z_z_z_xyyyzz[k];

                g_z_z_xz_yyzzz[k] = -g_z_z_z_yyzzz[k] * ab_x + g_z_z_z_xyyzzz[k];

                g_z_z_xz_yzzzz[k] = -g_z_z_z_yzzzz[k] * ab_x + g_z_z_z_xyzzzz[k];

                g_z_z_xz_zzzzz[k] = -g_z_z_z_zzzzz[k] * ab_x + g_z_z_z_xzzzzz[k];
            }

            /// Set up 1071-1092 components of targeted buffer : cbuffer.data(

            auto g_z_z_yy_xxxxx = cbuffer.data(dh_geom_11_off + 1071 * ccomps * dcomps);

            auto g_z_z_yy_xxxxy = cbuffer.data(dh_geom_11_off + 1072 * ccomps * dcomps);

            auto g_z_z_yy_xxxxz = cbuffer.data(dh_geom_11_off + 1073 * ccomps * dcomps);

            auto g_z_z_yy_xxxyy = cbuffer.data(dh_geom_11_off + 1074 * ccomps * dcomps);

            auto g_z_z_yy_xxxyz = cbuffer.data(dh_geom_11_off + 1075 * ccomps * dcomps);

            auto g_z_z_yy_xxxzz = cbuffer.data(dh_geom_11_off + 1076 * ccomps * dcomps);

            auto g_z_z_yy_xxyyy = cbuffer.data(dh_geom_11_off + 1077 * ccomps * dcomps);

            auto g_z_z_yy_xxyyz = cbuffer.data(dh_geom_11_off + 1078 * ccomps * dcomps);

            auto g_z_z_yy_xxyzz = cbuffer.data(dh_geom_11_off + 1079 * ccomps * dcomps);

            auto g_z_z_yy_xxzzz = cbuffer.data(dh_geom_11_off + 1080 * ccomps * dcomps);

            auto g_z_z_yy_xyyyy = cbuffer.data(dh_geom_11_off + 1081 * ccomps * dcomps);

            auto g_z_z_yy_xyyyz = cbuffer.data(dh_geom_11_off + 1082 * ccomps * dcomps);

            auto g_z_z_yy_xyyzz = cbuffer.data(dh_geom_11_off + 1083 * ccomps * dcomps);

            auto g_z_z_yy_xyzzz = cbuffer.data(dh_geom_11_off + 1084 * ccomps * dcomps);

            auto g_z_z_yy_xzzzz = cbuffer.data(dh_geom_11_off + 1085 * ccomps * dcomps);

            auto g_z_z_yy_yyyyy = cbuffer.data(dh_geom_11_off + 1086 * ccomps * dcomps);

            auto g_z_z_yy_yyyyz = cbuffer.data(dh_geom_11_off + 1087 * ccomps * dcomps);

            auto g_z_z_yy_yyyzz = cbuffer.data(dh_geom_11_off + 1088 * ccomps * dcomps);

            auto g_z_z_yy_yyzzz = cbuffer.data(dh_geom_11_off + 1089 * ccomps * dcomps);

            auto g_z_z_yy_yzzzz = cbuffer.data(dh_geom_11_off + 1090 * ccomps * dcomps);

            auto g_z_z_yy_zzzzz = cbuffer.data(dh_geom_11_off + 1091 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_y_xxxxx, g_z_z_y_xxxxxy, g_z_z_y_xxxxy, g_z_z_y_xxxxyy, g_z_z_y_xxxxyz, g_z_z_y_xxxxz, g_z_z_y_xxxyy, g_z_z_y_xxxyyy, g_z_z_y_xxxyyz, g_z_z_y_xxxyz, g_z_z_y_xxxyzz, g_z_z_y_xxxzz, g_z_z_y_xxyyy, g_z_z_y_xxyyyy, g_z_z_y_xxyyyz, g_z_z_y_xxyyz, g_z_z_y_xxyyzz, g_z_z_y_xxyzz, g_z_z_y_xxyzzz, g_z_z_y_xxzzz, g_z_z_y_xyyyy, g_z_z_y_xyyyyy, g_z_z_y_xyyyyz, g_z_z_y_xyyyz, g_z_z_y_xyyyzz, g_z_z_y_xyyzz, g_z_z_y_xyyzzz, g_z_z_y_xyzzz, g_z_z_y_xyzzzz, g_z_z_y_xzzzz, g_z_z_y_yyyyy, g_z_z_y_yyyyyy, g_z_z_y_yyyyyz, g_z_z_y_yyyyz, g_z_z_y_yyyyzz, g_z_z_y_yyyzz, g_z_z_y_yyyzzz, g_z_z_y_yyzzz, g_z_z_y_yyzzzz, g_z_z_y_yzzzz, g_z_z_y_yzzzzz, g_z_z_y_zzzzz, g_z_z_yy_xxxxx, g_z_z_yy_xxxxy, g_z_z_yy_xxxxz, g_z_z_yy_xxxyy, g_z_z_yy_xxxyz, g_z_z_yy_xxxzz, g_z_z_yy_xxyyy, g_z_z_yy_xxyyz, g_z_z_yy_xxyzz, g_z_z_yy_xxzzz, g_z_z_yy_xyyyy, g_z_z_yy_xyyyz, g_z_z_yy_xyyzz, g_z_z_yy_xyzzz, g_z_z_yy_xzzzz, g_z_z_yy_yyyyy, g_z_z_yy_yyyyz, g_z_z_yy_yyyzz, g_z_z_yy_yyzzz, g_z_z_yy_yzzzz, g_z_z_yy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yy_xxxxx[k] = -g_z_z_y_xxxxx[k] * ab_y + g_z_z_y_xxxxxy[k];

                g_z_z_yy_xxxxy[k] = -g_z_z_y_xxxxy[k] * ab_y + g_z_z_y_xxxxyy[k];

                g_z_z_yy_xxxxz[k] = -g_z_z_y_xxxxz[k] * ab_y + g_z_z_y_xxxxyz[k];

                g_z_z_yy_xxxyy[k] = -g_z_z_y_xxxyy[k] * ab_y + g_z_z_y_xxxyyy[k];

                g_z_z_yy_xxxyz[k] = -g_z_z_y_xxxyz[k] * ab_y + g_z_z_y_xxxyyz[k];

                g_z_z_yy_xxxzz[k] = -g_z_z_y_xxxzz[k] * ab_y + g_z_z_y_xxxyzz[k];

                g_z_z_yy_xxyyy[k] = -g_z_z_y_xxyyy[k] * ab_y + g_z_z_y_xxyyyy[k];

                g_z_z_yy_xxyyz[k] = -g_z_z_y_xxyyz[k] * ab_y + g_z_z_y_xxyyyz[k];

                g_z_z_yy_xxyzz[k] = -g_z_z_y_xxyzz[k] * ab_y + g_z_z_y_xxyyzz[k];

                g_z_z_yy_xxzzz[k] = -g_z_z_y_xxzzz[k] * ab_y + g_z_z_y_xxyzzz[k];

                g_z_z_yy_xyyyy[k] = -g_z_z_y_xyyyy[k] * ab_y + g_z_z_y_xyyyyy[k];

                g_z_z_yy_xyyyz[k] = -g_z_z_y_xyyyz[k] * ab_y + g_z_z_y_xyyyyz[k];

                g_z_z_yy_xyyzz[k] = -g_z_z_y_xyyzz[k] * ab_y + g_z_z_y_xyyyzz[k];

                g_z_z_yy_xyzzz[k] = -g_z_z_y_xyzzz[k] * ab_y + g_z_z_y_xyyzzz[k];

                g_z_z_yy_xzzzz[k] = -g_z_z_y_xzzzz[k] * ab_y + g_z_z_y_xyzzzz[k];

                g_z_z_yy_yyyyy[k] = -g_z_z_y_yyyyy[k] * ab_y + g_z_z_y_yyyyyy[k];

                g_z_z_yy_yyyyz[k] = -g_z_z_y_yyyyz[k] * ab_y + g_z_z_y_yyyyyz[k];

                g_z_z_yy_yyyzz[k] = -g_z_z_y_yyyzz[k] * ab_y + g_z_z_y_yyyyzz[k];

                g_z_z_yy_yyzzz[k] = -g_z_z_y_yyzzz[k] * ab_y + g_z_z_y_yyyzzz[k];

                g_z_z_yy_yzzzz[k] = -g_z_z_y_yzzzz[k] * ab_y + g_z_z_y_yyzzzz[k];

                g_z_z_yy_zzzzz[k] = -g_z_z_y_zzzzz[k] * ab_y + g_z_z_y_yzzzzz[k];
            }

            /// Set up 1092-1113 components of targeted buffer : cbuffer.data(

            auto g_z_z_yz_xxxxx = cbuffer.data(dh_geom_11_off + 1092 * ccomps * dcomps);

            auto g_z_z_yz_xxxxy = cbuffer.data(dh_geom_11_off + 1093 * ccomps * dcomps);

            auto g_z_z_yz_xxxxz = cbuffer.data(dh_geom_11_off + 1094 * ccomps * dcomps);

            auto g_z_z_yz_xxxyy = cbuffer.data(dh_geom_11_off + 1095 * ccomps * dcomps);

            auto g_z_z_yz_xxxyz = cbuffer.data(dh_geom_11_off + 1096 * ccomps * dcomps);

            auto g_z_z_yz_xxxzz = cbuffer.data(dh_geom_11_off + 1097 * ccomps * dcomps);

            auto g_z_z_yz_xxyyy = cbuffer.data(dh_geom_11_off + 1098 * ccomps * dcomps);

            auto g_z_z_yz_xxyyz = cbuffer.data(dh_geom_11_off + 1099 * ccomps * dcomps);

            auto g_z_z_yz_xxyzz = cbuffer.data(dh_geom_11_off + 1100 * ccomps * dcomps);

            auto g_z_z_yz_xxzzz = cbuffer.data(dh_geom_11_off + 1101 * ccomps * dcomps);

            auto g_z_z_yz_xyyyy = cbuffer.data(dh_geom_11_off + 1102 * ccomps * dcomps);

            auto g_z_z_yz_xyyyz = cbuffer.data(dh_geom_11_off + 1103 * ccomps * dcomps);

            auto g_z_z_yz_xyyzz = cbuffer.data(dh_geom_11_off + 1104 * ccomps * dcomps);

            auto g_z_z_yz_xyzzz = cbuffer.data(dh_geom_11_off + 1105 * ccomps * dcomps);

            auto g_z_z_yz_xzzzz = cbuffer.data(dh_geom_11_off + 1106 * ccomps * dcomps);

            auto g_z_z_yz_yyyyy = cbuffer.data(dh_geom_11_off + 1107 * ccomps * dcomps);

            auto g_z_z_yz_yyyyz = cbuffer.data(dh_geom_11_off + 1108 * ccomps * dcomps);

            auto g_z_z_yz_yyyzz = cbuffer.data(dh_geom_11_off + 1109 * ccomps * dcomps);

            auto g_z_z_yz_yyzzz = cbuffer.data(dh_geom_11_off + 1110 * ccomps * dcomps);

            auto g_z_z_yz_yzzzz = cbuffer.data(dh_geom_11_off + 1111 * ccomps * dcomps);

            auto g_z_z_yz_zzzzz = cbuffer.data(dh_geom_11_off + 1112 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_z_yz_xxxxx, g_z_z_yz_xxxxy, g_z_z_yz_xxxxz, g_z_z_yz_xxxyy, g_z_z_yz_xxxyz, g_z_z_yz_xxxzz, g_z_z_yz_xxyyy, g_z_z_yz_xxyyz, g_z_z_yz_xxyzz, g_z_z_yz_xxzzz, g_z_z_yz_xyyyy, g_z_z_yz_xyyyz, g_z_z_yz_xyyzz, g_z_z_yz_xyzzz, g_z_z_yz_xzzzz, g_z_z_yz_yyyyy, g_z_z_yz_yyyyz, g_z_z_yz_yyyzz, g_z_z_yz_yyzzz, g_z_z_yz_yzzzz, g_z_z_yz_zzzzz, g_z_z_z_xxxxx, g_z_z_z_xxxxxy, g_z_z_z_xxxxy, g_z_z_z_xxxxyy, g_z_z_z_xxxxyz, g_z_z_z_xxxxz, g_z_z_z_xxxyy, g_z_z_z_xxxyyy, g_z_z_z_xxxyyz, g_z_z_z_xxxyz, g_z_z_z_xxxyzz, g_z_z_z_xxxzz, g_z_z_z_xxyyy, g_z_z_z_xxyyyy, g_z_z_z_xxyyyz, g_z_z_z_xxyyz, g_z_z_z_xxyyzz, g_z_z_z_xxyzz, g_z_z_z_xxyzzz, g_z_z_z_xxzzz, g_z_z_z_xyyyy, g_z_z_z_xyyyyy, g_z_z_z_xyyyyz, g_z_z_z_xyyyz, g_z_z_z_xyyyzz, g_z_z_z_xyyzz, g_z_z_z_xyyzzz, g_z_z_z_xyzzz, g_z_z_z_xyzzzz, g_z_z_z_xzzzz, g_z_z_z_yyyyy, g_z_z_z_yyyyyy, g_z_z_z_yyyyyz, g_z_z_z_yyyyz, g_z_z_z_yyyyzz, g_z_z_z_yyyzz, g_z_z_z_yyyzzz, g_z_z_z_yyzzz, g_z_z_z_yyzzzz, g_z_z_z_yzzzz, g_z_z_z_yzzzzz, g_z_z_z_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_yz_xxxxx[k] = -g_z_z_z_xxxxx[k] * ab_y + g_z_z_z_xxxxxy[k];

                g_z_z_yz_xxxxy[k] = -g_z_z_z_xxxxy[k] * ab_y + g_z_z_z_xxxxyy[k];

                g_z_z_yz_xxxxz[k] = -g_z_z_z_xxxxz[k] * ab_y + g_z_z_z_xxxxyz[k];

                g_z_z_yz_xxxyy[k] = -g_z_z_z_xxxyy[k] * ab_y + g_z_z_z_xxxyyy[k];

                g_z_z_yz_xxxyz[k] = -g_z_z_z_xxxyz[k] * ab_y + g_z_z_z_xxxyyz[k];

                g_z_z_yz_xxxzz[k] = -g_z_z_z_xxxzz[k] * ab_y + g_z_z_z_xxxyzz[k];

                g_z_z_yz_xxyyy[k] = -g_z_z_z_xxyyy[k] * ab_y + g_z_z_z_xxyyyy[k];

                g_z_z_yz_xxyyz[k] = -g_z_z_z_xxyyz[k] * ab_y + g_z_z_z_xxyyyz[k];

                g_z_z_yz_xxyzz[k] = -g_z_z_z_xxyzz[k] * ab_y + g_z_z_z_xxyyzz[k];

                g_z_z_yz_xxzzz[k] = -g_z_z_z_xxzzz[k] * ab_y + g_z_z_z_xxyzzz[k];

                g_z_z_yz_xyyyy[k] = -g_z_z_z_xyyyy[k] * ab_y + g_z_z_z_xyyyyy[k];

                g_z_z_yz_xyyyz[k] = -g_z_z_z_xyyyz[k] * ab_y + g_z_z_z_xyyyyz[k];

                g_z_z_yz_xyyzz[k] = -g_z_z_z_xyyzz[k] * ab_y + g_z_z_z_xyyyzz[k];

                g_z_z_yz_xyzzz[k] = -g_z_z_z_xyzzz[k] * ab_y + g_z_z_z_xyyzzz[k];

                g_z_z_yz_xzzzz[k] = -g_z_z_z_xzzzz[k] * ab_y + g_z_z_z_xyzzzz[k];

                g_z_z_yz_yyyyy[k] = -g_z_z_z_yyyyy[k] * ab_y + g_z_z_z_yyyyyy[k];

                g_z_z_yz_yyyyz[k] = -g_z_z_z_yyyyz[k] * ab_y + g_z_z_z_yyyyyz[k];

                g_z_z_yz_yyyzz[k] = -g_z_z_z_yyyzz[k] * ab_y + g_z_z_z_yyyyzz[k];

                g_z_z_yz_yyzzz[k] = -g_z_z_z_yyzzz[k] * ab_y + g_z_z_z_yyyzzz[k];

                g_z_z_yz_yzzzz[k] = -g_z_z_z_yzzzz[k] * ab_y + g_z_z_z_yyzzzz[k];

                g_z_z_yz_zzzzz[k] = -g_z_z_z_zzzzz[k] * ab_y + g_z_z_z_yzzzzz[k];
            }

            /// Set up 1113-1134 components of targeted buffer : cbuffer.data(

            auto g_z_z_zz_xxxxx = cbuffer.data(dh_geom_11_off + 1113 * ccomps * dcomps);

            auto g_z_z_zz_xxxxy = cbuffer.data(dh_geom_11_off + 1114 * ccomps * dcomps);

            auto g_z_z_zz_xxxxz = cbuffer.data(dh_geom_11_off + 1115 * ccomps * dcomps);

            auto g_z_z_zz_xxxyy = cbuffer.data(dh_geom_11_off + 1116 * ccomps * dcomps);

            auto g_z_z_zz_xxxyz = cbuffer.data(dh_geom_11_off + 1117 * ccomps * dcomps);

            auto g_z_z_zz_xxxzz = cbuffer.data(dh_geom_11_off + 1118 * ccomps * dcomps);

            auto g_z_z_zz_xxyyy = cbuffer.data(dh_geom_11_off + 1119 * ccomps * dcomps);

            auto g_z_z_zz_xxyyz = cbuffer.data(dh_geom_11_off + 1120 * ccomps * dcomps);

            auto g_z_z_zz_xxyzz = cbuffer.data(dh_geom_11_off + 1121 * ccomps * dcomps);

            auto g_z_z_zz_xxzzz = cbuffer.data(dh_geom_11_off + 1122 * ccomps * dcomps);

            auto g_z_z_zz_xyyyy = cbuffer.data(dh_geom_11_off + 1123 * ccomps * dcomps);

            auto g_z_z_zz_xyyyz = cbuffer.data(dh_geom_11_off + 1124 * ccomps * dcomps);

            auto g_z_z_zz_xyyzz = cbuffer.data(dh_geom_11_off + 1125 * ccomps * dcomps);

            auto g_z_z_zz_xyzzz = cbuffer.data(dh_geom_11_off + 1126 * ccomps * dcomps);

            auto g_z_z_zz_xzzzz = cbuffer.data(dh_geom_11_off + 1127 * ccomps * dcomps);

            auto g_z_z_zz_yyyyy = cbuffer.data(dh_geom_11_off + 1128 * ccomps * dcomps);

            auto g_z_z_zz_yyyyz = cbuffer.data(dh_geom_11_off + 1129 * ccomps * dcomps);

            auto g_z_z_zz_yyyzz = cbuffer.data(dh_geom_11_off + 1130 * ccomps * dcomps);

            auto g_z_z_zz_yyzzz = cbuffer.data(dh_geom_11_off + 1131 * ccomps * dcomps);

            auto g_z_z_zz_yzzzz = cbuffer.data(dh_geom_11_off + 1132 * ccomps * dcomps);

            auto g_z_z_zz_zzzzz = cbuffer.data(dh_geom_11_off + 1133 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_z_xxxxx, g_0_z_z_xxxxy, g_0_z_z_xxxxz, g_0_z_z_xxxyy, g_0_z_z_xxxyz, g_0_z_z_xxxzz, g_0_z_z_xxyyy, g_0_z_z_xxyyz, g_0_z_z_xxyzz, g_0_z_z_xxzzz, g_0_z_z_xyyyy, g_0_z_z_xyyyz, g_0_z_z_xyyzz, g_0_z_z_xyzzz, g_0_z_z_xzzzz, g_0_z_z_yyyyy, g_0_z_z_yyyyz, g_0_z_z_yyyzz, g_0_z_z_yyzzz, g_0_z_z_yzzzz, g_0_z_z_zzzzz, g_z_0_z_xxxxx, g_z_0_z_xxxxy, g_z_0_z_xxxxz, g_z_0_z_xxxyy, g_z_0_z_xxxyz, g_z_0_z_xxxzz, g_z_0_z_xxyyy, g_z_0_z_xxyyz, g_z_0_z_xxyzz, g_z_0_z_xxzzz, g_z_0_z_xyyyy, g_z_0_z_xyyyz, g_z_0_z_xyyzz, g_z_0_z_xyzzz, g_z_0_z_xzzzz, g_z_0_z_yyyyy, g_z_0_z_yyyyz, g_z_0_z_yyyzz, g_z_0_z_yyzzz, g_z_0_z_yzzzz, g_z_0_z_zzzzz, g_z_z_z_xxxxx, g_z_z_z_xxxxxz, g_z_z_z_xxxxy, g_z_z_z_xxxxyz, g_z_z_z_xxxxz, g_z_z_z_xxxxzz, g_z_z_z_xxxyy, g_z_z_z_xxxyyz, g_z_z_z_xxxyz, g_z_z_z_xxxyzz, g_z_z_z_xxxzz, g_z_z_z_xxxzzz, g_z_z_z_xxyyy, g_z_z_z_xxyyyz, g_z_z_z_xxyyz, g_z_z_z_xxyyzz, g_z_z_z_xxyzz, g_z_z_z_xxyzzz, g_z_z_z_xxzzz, g_z_z_z_xxzzzz, g_z_z_z_xyyyy, g_z_z_z_xyyyyz, g_z_z_z_xyyyz, g_z_z_z_xyyyzz, g_z_z_z_xyyzz, g_z_z_z_xyyzzz, g_z_z_z_xyzzz, g_z_z_z_xyzzzz, g_z_z_z_xzzzz, g_z_z_z_xzzzzz, g_z_z_z_yyyyy, g_z_z_z_yyyyyz, g_z_z_z_yyyyz, g_z_z_z_yyyyzz, g_z_z_z_yyyzz, g_z_z_z_yyyzzz, g_z_z_z_yyzzz, g_z_z_z_yyzzzz, g_z_z_z_yzzzz, g_z_z_z_yzzzzz, g_z_z_z_zzzzz, g_z_z_z_zzzzzz, g_z_z_zz_xxxxx, g_z_z_zz_xxxxy, g_z_z_zz_xxxxz, g_z_z_zz_xxxyy, g_z_z_zz_xxxyz, g_z_z_zz_xxxzz, g_z_z_zz_xxyyy, g_z_z_zz_xxyyz, g_z_z_zz_xxyzz, g_z_z_zz_xxzzz, g_z_z_zz_xyyyy, g_z_z_zz_xyyyz, g_z_z_zz_xyyzz, g_z_z_zz_xyzzz, g_z_z_zz_xzzzz, g_z_z_zz_yyyyy, g_z_z_zz_yyyyz, g_z_z_zz_yyyzz, g_z_z_zz_yyzzz, g_z_z_zz_yzzzz, g_z_z_zz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_z_zz_xxxxx[k] = -g_0_z_z_xxxxx[k] + g_z_0_z_xxxxx[k] - g_z_z_z_xxxxx[k] * ab_z + g_z_z_z_xxxxxz[k];

                g_z_z_zz_xxxxy[k] = -g_0_z_z_xxxxy[k] + g_z_0_z_xxxxy[k] - g_z_z_z_xxxxy[k] * ab_z + g_z_z_z_xxxxyz[k];

                g_z_z_zz_xxxxz[k] = -g_0_z_z_xxxxz[k] + g_z_0_z_xxxxz[k] - g_z_z_z_xxxxz[k] * ab_z + g_z_z_z_xxxxzz[k];

                g_z_z_zz_xxxyy[k] = -g_0_z_z_xxxyy[k] + g_z_0_z_xxxyy[k] - g_z_z_z_xxxyy[k] * ab_z + g_z_z_z_xxxyyz[k];

                g_z_z_zz_xxxyz[k] = -g_0_z_z_xxxyz[k] + g_z_0_z_xxxyz[k] - g_z_z_z_xxxyz[k] * ab_z + g_z_z_z_xxxyzz[k];

                g_z_z_zz_xxxzz[k] = -g_0_z_z_xxxzz[k] + g_z_0_z_xxxzz[k] - g_z_z_z_xxxzz[k] * ab_z + g_z_z_z_xxxzzz[k];

                g_z_z_zz_xxyyy[k] = -g_0_z_z_xxyyy[k] + g_z_0_z_xxyyy[k] - g_z_z_z_xxyyy[k] * ab_z + g_z_z_z_xxyyyz[k];

                g_z_z_zz_xxyyz[k] = -g_0_z_z_xxyyz[k] + g_z_0_z_xxyyz[k] - g_z_z_z_xxyyz[k] * ab_z + g_z_z_z_xxyyzz[k];

                g_z_z_zz_xxyzz[k] = -g_0_z_z_xxyzz[k] + g_z_0_z_xxyzz[k] - g_z_z_z_xxyzz[k] * ab_z + g_z_z_z_xxyzzz[k];

                g_z_z_zz_xxzzz[k] = -g_0_z_z_xxzzz[k] + g_z_0_z_xxzzz[k] - g_z_z_z_xxzzz[k] * ab_z + g_z_z_z_xxzzzz[k];

                g_z_z_zz_xyyyy[k] = -g_0_z_z_xyyyy[k] + g_z_0_z_xyyyy[k] - g_z_z_z_xyyyy[k] * ab_z + g_z_z_z_xyyyyz[k];

                g_z_z_zz_xyyyz[k] = -g_0_z_z_xyyyz[k] + g_z_0_z_xyyyz[k] - g_z_z_z_xyyyz[k] * ab_z + g_z_z_z_xyyyzz[k];

                g_z_z_zz_xyyzz[k] = -g_0_z_z_xyyzz[k] + g_z_0_z_xyyzz[k] - g_z_z_z_xyyzz[k] * ab_z + g_z_z_z_xyyzzz[k];

                g_z_z_zz_xyzzz[k] = -g_0_z_z_xyzzz[k] + g_z_0_z_xyzzz[k] - g_z_z_z_xyzzz[k] * ab_z + g_z_z_z_xyzzzz[k];

                g_z_z_zz_xzzzz[k] = -g_0_z_z_xzzzz[k] + g_z_0_z_xzzzz[k] - g_z_z_z_xzzzz[k] * ab_z + g_z_z_z_xzzzzz[k];

                g_z_z_zz_yyyyy[k] = -g_0_z_z_yyyyy[k] + g_z_0_z_yyyyy[k] - g_z_z_z_yyyyy[k] * ab_z + g_z_z_z_yyyyyz[k];

                g_z_z_zz_yyyyz[k] = -g_0_z_z_yyyyz[k] + g_z_0_z_yyyyz[k] - g_z_z_z_yyyyz[k] * ab_z + g_z_z_z_yyyyzz[k];

                g_z_z_zz_yyyzz[k] = -g_0_z_z_yyyzz[k] + g_z_0_z_yyyzz[k] - g_z_z_z_yyyzz[k] * ab_z + g_z_z_z_yyyzzz[k];

                g_z_z_zz_yyzzz[k] = -g_0_z_z_yyzzz[k] + g_z_0_z_yyzzz[k] - g_z_z_z_yyzzz[k] * ab_z + g_z_z_z_yyzzzz[k];

                g_z_z_zz_yzzzz[k] = -g_0_z_z_yzzzz[k] + g_z_0_z_yzzzz[k] - g_z_z_z_yzzzz[k] * ab_z + g_z_z_z_yzzzzz[k];

                g_z_z_zz_zzzzz[k] = -g_0_z_z_zzzzz[k] + g_z_0_z_zzzzz[k] - g_z_z_z_zzzzz[k] * ab_z + g_z_z_z_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

