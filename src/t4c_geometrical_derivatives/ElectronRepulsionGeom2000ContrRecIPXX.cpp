#include "ElectronRepulsionGeom2000ContrRecIPXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom20_hrr_electron_repulsion_ipxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_20_ipxx,
                                            const size_t idx_geom_10_hpxx,
                                            const size_t idx_geom_20_hpxx,
                                            const size_t idx_geom_20_hdxx,
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
            /// Set up components of auxilary buffer : HPSS

            const auto hp_geom_10_off = idx_geom_10_hpxx + i * dcomps + j;

            auto g_x_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 62 * ccomps * dcomps);

            auto g_y_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 63 * ccomps * dcomps);

            auto g_y_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 64 * ccomps * dcomps);

            auto g_y_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 65 * ccomps * dcomps);

            auto g_y_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 66 * ccomps * dcomps);

            auto g_y_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 67 * ccomps * dcomps);

            auto g_y_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 68 * ccomps * dcomps);

            auto g_y_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 69 * ccomps * dcomps);

            auto g_y_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 70 * ccomps * dcomps);

            auto g_y_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 71 * ccomps * dcomps);

            auto g_y_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 72 * ccomps * dcomps);

            auto g_y_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 73 * ccomps * dcomps);

            auto g_y_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 74 * ccomps * dcomps);

            auto g_y_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 75 * ccomps * dcomps);

            auto g_y_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 76 * ccomps * dcomps);

            auto g_y_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 77 * ccomps * dcomps);

            auto g_y_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 78 * ccomps * dcomps);

            auto g_y_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 79 * ccomps * dcomps);

            auto g_y_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 80 * ccomps * dcomps);

            auto g_y_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 81 * ccomps * dcomps);

            auto g_y_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 82 * ccomps * dcomps);

            auto g_y_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 83 * ccomps * dcomps);

            auto g_y_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 84 * ccomps * dcomps);

            auto g_y_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 85 * ccomps * dcomps);

            auto g_y_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 86 * ccomps * dcomps);

            auto g_y_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 87 * ccomps * dcomps);

            auto g_y_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 88 * ccomps * dcomps);

            auto g_y_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 89 * ccomps * dcomps);

            auto g_y_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 90 * ccomps * dcomps);

            auto g_y_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 91 * ccomps * dcomps);

            auto g_y_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 92 * ccomps * dcomps);

            auto g_y_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 93 * ccomps * dcomps);

            auto g_y_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 94 * ccomps * dcomps);

            auto g_y_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 95 * ccomps * dcomps);

            auto g_y_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 96 * ccomps * dcomps);

            auto g_y_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 97 * ccomps * dcomps);

            auto g_y_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 98 * ccomps * dcomps);

            auto g_y_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 99 * ccomps * dcomps);

            auto g_y_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 100 * ccomps * dcomps);

            auto g_y_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 101 * ccomps * dcomps);

            auto g_y_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 102 * ccomps * dcomps);

            auto g_y_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 103 * ccomps * dcomps);

            auto g_y_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 104 * ccomps * dcomps);

            auto g_y_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 105 * ccomps * dcomps);

            auto g_y_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 106 * ccomps * dcomps);

            auto g_y_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 107 * ccomps * dcomps);

            auto g_y_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 108 * ccomps * dcomps);

            auto g_y_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 109 * ccomps * dcomps);

            auto g_y_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 110 * ccomps * dcomps);

            auto g_y_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 111 * ccomps * dcomps);

            auto g_y_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 112 * ccomps * dcomps);

            auto g_y_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 113 * ccomps * dcomps);

            auto g_y_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 114 * ccomps * dcomps);

            auto g_y_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 115 * ccomps * dcomps);

            auto g_y_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 116 * ccomps * dcomps);

            auto g_y_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 117 * ccomps * dcomps);

            auto g_y_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 118 * ccomps * dcomps);

            auto g_y_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 119 * ccomps * dcomps);

            auto g_y_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 120 * ccomps * dcomps);

            auto g_y_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 121 * ccomps * dcomps);

            auto g_y_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 122 * ccomps * dcomps);

            auto g_y_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 123 * ccomps * dcomps);

            auto g_y_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 124 * ccomps * dcomps);

            auto g_y_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 125 * ccomps * dcomps);

            auto g_z_0_xxxxx_x = cbuffer.data(hp_geom_10_off + 126 * ccomps * dcomps);

            auto g_z_0_xxxxx_y = cbuffer.data(hp_geom_10_off + 127 * ccomps * dcomps);

            auto g_z_0_xxxxx_z = cbuffer.data(hp_geom_10_off + 128 * ccomps * dcomps);

            auto g_z_0_xxxxy_x = cbuffer.data(hp_geom_10_off + 129 * ccomps * dcomps);

            auto g_z_0_xxxxy_y = cbuffer.data(hp_geom_10_off + 130 * ccomps * dcomps);

            auto g_z_0_xxxxy_z = cbuffer.data(hp_geom_10_off + 131 * ccomps * dcomps);

            auto g_z_0_xxxxz_x = cbuffer.data(hp_geom_10_off + 132 * ccomps * dcomps);

            auto g_z_0_xxxxz_y = cbuffer.data(hp_geom_10_off + 133 * ccomps * dcomps);

            auto g_z_0_xxxxz_z = cbuffer.data(hp_geom_10_off + 134 * ccomps * dcomps);

            auto g_z_0_xxxyy_x = cbuffer.data(hp_geom_10_off + 135 * ccomps * dcomps);

            auto g_z_0_xxxyy_y = cbuffer.data(hp_geom_10_off + 136 * ccomps * dcomps);

            auto g_z_0_xxxyy_z = cbuffer.data(hp_geom_10_off + 137 * ccomps * dcomps);

            auto g_z_0_xxxyz_x = cbuffer.data(hp_geom_10_off + 138 * ccomps * dcomps);

            auto g_z_0_xxxyz_y = cbuffer.data(hp_geom_10_off + 139 * ccomps * dcomps);

            auto g_z_0_xxxyz_z = cbuffer.data(hp_geom_10_off + 140 * ccomps * dcomps);

            auto g_z_0_xxxzz_x = cbuffer.data(hp_geom_10_off + 141 * ccomps * dcomps);

            auto g_z_0_xxxzz_y = cbuffer.data(hp_geom_10_off + 142 * ccomps * dcomps);

            auto g_z_0_xxxzz_z = cbuffer.data(hp_geom_10_off + 143 * ccomps * dcomps);

            auto g_z_0_xxyyy_x = cbuffer.data(hp_geom_10_off + 144 * ccomps * dcomps);

            auto g_z_0_xxyyy_y = cbuffer.data(hp_geom_10_off + 145 * ccomps * dcomps);

            auto g_z_0_xxyyy_z = cbuffer.data(hp_geom_10_off + 146 * ccomps * dcomps);

            auto g_z_0_xxyyz_x = cbuffer.data(hp_geom_10_off + 147 * ccomps * dcomps);

            auto g_z_0_xxyyz_y = cbuffer.data(hp_geom_10_off + 148 * ccomps * dcomps);

            auto g_z_0_xxyyz_z = cbuffer.data(hp_geom_10_off + 149 * ccomps * dcomps);

            auto g_z_0_xxyzz_x = cbuffer.data(hp_geom_10_off + 150 * ccomps * dcomps);

            auto g_z_0_xxyzz_y = cbuffer.data(hp_geom_10_off + 151 * ccomps * dcomps);

            auto g_z_0_xxyzz_z = cbuffer.data(hp_geom_10_off + 152 * ccomps * dcomps);

            auto g_z_0_xxzzz_x = cbuffer.data(hp_geom_10_off + 153 * ccomps * dcomps);

            auto g_z_0_xxzzz_y = cbuffer.data(hp_geom_10_off + 154 * ccomps * dcomps);

            auto g_z_0_xxzzz_z = cbuffer.data(hp_geom_10_off + 155 * ccomps * dcomps);

            auto g_z_0_xyyyy_x = cbuffer.data(hp_geom_10_off + 156 * ccomps * dcomps);

            auto g_z_0_xyyyy_y = cbuffer.data(hp_geom_10_off + 157 * ccomps * dcomps);

            auto g_z_0_xyyyy_z = cbuffer.data(hp_geom_10_off + 158 * ccomps * dcomps);

            auto g_z_0_xyyyz_x = cbuffer.data(hp_geom_10_off + 159 * ccomps * dcomps);

            auto g_z_0_xyyyz_y = cbuffer.data(hp_geom_10_off + 160 * ccomps * dcomps);

            auto g_z_0_xyyyz_z = cbuffer.data(hp_geom_10_off + 161 * ccomps * dcomps);

            auto g_z_0_xyyzz_x = cbuffer.data(hp_geom_10_off + 162 * ccomps * dcomps);

            auto g_z_0_xyyzz_y = cbuffer.data(hp_geom_10_off + 163 * ccomps * dcomps);

            auto g_z_0_xyyzz_z = cbuffer.data(hp_geom_10_off + 164 * ccomps * dcomps);

            auto g_z_0_xyzzz_x = cbuffer.data(hp_geom_10_off + 165 * ccomps * dcomps);

            auto g_z_0_xyzzz_y = cbuffer.data(hp_geom_10_off + 166 * ccomps * dcomps);

            auto g_z_0_xyzzz_z = cbuffer.data(hp_geom_10_off + 167 * ccomps * dcomps);

            auto g_z_0_xzzzz_x = cbuffer.data(hp_geom_10_off + 168 * ccomps * dcomps);

            auto g_z_0_xzzzz_y = cbuffer.data(hp_geom_10_off + 169 * ccomps * dcomps);

            auto g_z_0_xzzzz_z = cbuffer.data(hp_geom_10_off + 170 * ccomps * dcomps);

            auto g_z_0_yyyyy_x = cbuffer.data(hp_geom_10_off + 171 * ccomps * dcomps);

            auto g_z_0_yyyyy_y = cbuffer.data(hp_geom_10_off + 172 * ccomps * dcomps);

            auto g_z_0_yyyyy_z = cbuffer.data(hp_geom_10_off + 173 * ccomps * dcomps);

            auto g_z_0_yyyyz_x = cbuffer.data(hp_geom_10_off + 174 * ccomps * dcomps);

            auto g_z_0_yyyyz_y = cbuffer.data(hp_geom_10_off + 175 * ccomps * dcomps);

            auto g_z_0_yyyyz_z = cbuffer.data(hp_geom_10_off + 176 * ccomps * dcomps);

            auto g_z_0_yyyzz_x = cbuffer.data(hp_geom_10_off + 177 * ccomps * dcomps);

            auto g_z_0_yyyzz_y = cbuffer.data(hp_geom_10_off + 178 * ccomps * dcomps);

            auto g_z_0_yyyzz_z = cbuffer.data(hp_geom_10_off + 179 * ccomps * dcomps);

            auto g_z_0_yyzzz_x = cbuffer.data(hp_geom_10_off + 180 * ccomps * dcomps);

            auto g_z_0_yyzzz_y = cbuffer.data(hp_geom_10_off + 181 * ccomps * dcomps);

            auto g_z_0_yyzzz_z = cbuffer.data(hp_geom_10_off + 182 * ccomps * dcomps);

            auto g_z_0_yzzzz_x = cbuffer.data(hp_geom_10_off + 183 * ccomps * dcomps);

            auto g_z_0_yzzzz_y = cbuffer.data(hp_geom_10_off + 184 * ccomps * dcomps);

            auto g_z_0_yzzzz_z = cbuffer.data(hp_geom_10_off + 185 * ccomps * dcomps);

            auto g_z_0_zzzzz_x = cbuffer.data(hp_geom_10_off + 186 * ccomps * dcomps);

            auto g_z_0_zzzzz_y = cbuffer.data(hp_geom_10_off + 187 * ccomps * dcomps);

            auto g_z_0_zzzzz_z = cbuffer.data(hp_geom_10_off + 188 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HPSS

            const auto hp_geom_20_off = idx_geom_20_hpxx + i * dcomps + j;

            auto g_xx_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 62 * ccomps * dcomps);

            auto g_xy_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 63 * ccomps * dcomps);

            auto g_xy_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 64 * ccomps * dcomps);

            auto g_xy_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 65 * ccomps * dcomps);

            auto g_xy_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 66 * ccomps * dcomps);

            auto g_xy_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 67 * ccomps * dcomps);

            auto g_xy_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 68 * ccomps * dcomps);

            auto g_xy_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 69 * ccomps * dcomps);

            auto g_xy_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 70 * ccomps * dcomps);

            auto g_xy_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 71 * ccomps * dcomps);

            auto g_xy_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 72 * ccomps * dcomps);

            auto g_xy_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 73 * ccomps * dcomps);

            auto g_xy_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 74 * ccomps * dcomps);

            auto g_xy_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 75 * ccomps * dcomps);

            auto g_xy_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 76 * ccomps * dcomps);

            auto g_xy_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 77 * ccomps * dcomps);

            auto g_xy_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 78 * ccomps * dcomps);

            auto g_xy_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 79 * ccomps * dcomps);

            auto g_xy_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 80 * ccomps * dcomps);

            auto g_xy_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 81 * ccomps * dcomps);

            auto g_xy_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 82 * ccomps * dcomps);

            auto g_xy_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 83 * ccomps * dcomps);

            auto g_xy_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 86 * ccomps * dcomps);

            auto g_xy_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 89 * ccomps * dcomps);

            auto g_xy_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 92 * ccomps * dcomps);

            auto g_xy_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 95 * ccomps * dcomps);

            auto g_xy_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 98 * ccomps * dcomps);

            auto g_xy_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 101 * ccomps * dcomps);

            auto g_xy_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 104 * ccomps * dcomps);

            auto g_xy_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 107 * ccomps * dcomps);

            auto g_xy_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 110 * ccomps * dcomps);

            auto g_xy_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 113 * ccomps * dcomps);

            auto g_xy_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 116 * ccomps * dcomps);

            auto g_xy_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 119 * ccomps * dcomps);

            auto g_xy_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 122 * ccomps * dcomps);

            auto g_xy_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 125 * ccomps * dcomps);

            auto g_xz_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 126 * ccomps * dcomps);

            auto g_xz_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 127 * ccomps * dcomps);

            auto g_xz_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 128 * ccomps * dcomps);

            auto g_xz_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 129 * ccomps * dcomps);

            auto g_xz_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 130 * ccomps * dcomps);

            auto g_xz_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 131 * ccomps * dcomps);

            auto g_xz_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 132 * ccomps * dcomps);

            auto g_xz_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 133 * ccomps * dcomps);

            auto g_xz_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 134 * ccomps * dcomps);

            auto g_xz_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 135 * ccomps * dcomps);

            auto g_xz_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 136 * ccomps * dcomps);

            auto g_xz_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 137 * ccomps * dcomps);

            auto g_xz_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 138 * ccomps * dcomps);

            auto g_xz_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 139 * ccomps * dcomps);

            auto g_xz_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 140 * ccomps * dcomps);

            auto g_xz_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 141 * ccomps * dcomps);

            auto g_xz_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 142 * ccomps * dcomps);

            auto g_xz_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 143 * ccomps * dcomps);

            auto g_xz_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 144 * ccomps * dcomps);

            auto g_xz_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 145 * ccomps * dcomps);

            auto g_xz_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 146 * ccomps * dcomps);

            auto g_xz_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 147 * ccomps * dcomps);

            auto g_xz_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 148 * ccomps * dcomps);

            auto g_xz_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 149 * ccomps * dcomps);

            auto g_xz_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 150 * ccomps * dcomps);

            auto g_xz_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 151 * ccomps * dcomps);

            auto g_xz_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 152 * ccomps * dcomps);

            auto g_xz_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 153 * ccomps * dcomps);

            auto g_xz_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 154 * ccomps * dcomps);

            auto g_xz_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 155 * ccomps * dcomps);

            auto g_xz_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 156 * ccomps * dcomps);

            auto g_xz_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 157 * ccomps * dcomps);

            auto g_xz_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 158 * ccomps * dcomps);

            auto g_xz_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 159 * ccomps * dcomps);

            auto g_xz_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 160 * ccomps * dcomps);

            auto g_xz_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 161 * ccomps * dcomps);

            auto g_xz_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 162 * ccomps * dcomps);

            auto g_xz_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 163 * ccomps * dcomps);

            auto g_xz_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 164 * ccomps * dcomps);

            auto g_xz_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 165 * ccomps * dcomps);

            auto g_xz_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 166 * ccomps * dcomps);

            auto g_xz_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 167 * ccomps * dcomps);

            auto g_xz_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 168 * ccomps * dcomps);

            auto g_xz_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 169 * ccomps * dcomps);

            auto g_xz_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 170 * ccomps * dcomps);

            auto g_xz_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 171 * ccomps * dcomps);

            auto g_xz_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 172 * ccomps * dcomps);

            auto g_xz_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 173 * ccomps * dcomps);

            auto g_xz_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 174 * ccomps * dcomps);

            auto g_xz_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 175 * ccomps * dcomps);

            auto g_xz_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 176 * ccomps * dcomps);

            auto g_xz_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 177 * ccomps * dcomps);

            auto g_xz_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 178 * ccomps * dcomps);

            auto g_xz_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 179 * ccomps * dcomps);

            auto g_xz_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 182 * ccomps * dcomps);

            auto g_xz_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 185 * ccomps * dcomps);

            auto g_xz_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 188 * ccomps * dcomps);

            auto g_yy_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 189 * ccomps * dcomps);

            auto g_yy_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 190 * ccomps * dcomps);

            auto g_yy_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 191 * ccomps * dcomps);

            auto g_yy_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 192 * ccomps * dcomps);

            auto g_yy_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 193 * ccomps * dcomps);

            auto g_yy_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 194 * ccomps * dcomps);

            auto g_yy_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 195 * ccomps * dcomps);

            auto g_yy_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 196 * ccomps * dcomps);

            auto g_yy_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 197 * ccomps * dcomps);

            auto g_yy_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 198 * ccomps * dcomps);

            auto g_yy_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 199 * ccomps * dcomps);

            auto g_yy_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 200 * ccomps * dcomps);

            auto g_yy_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 201 * ccomps * dcomps);

            auto g_yy_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 202 * ccomps * dcomps);

            auto g_yy_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 203 * ccomps * dcomps);

            auto g_yy_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 204 * ccomps * dcomps);

            auto g_yy_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 205 * ccomps * dcomps);

            auto g_yy_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 206 * ccomps * dcomps);

            auto g_yy_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 207 * ccomps * dcomps);

            auto g_yy_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 208 * ccomps * dcomps);

            auto g_yy_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 209 * ccomps * dcomps);

            auto g_yy_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 210 * ccomps * dcomps);

            auto g_yy_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 211 * ccomps * dcomps);

            auto g_yy_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 212 * ccomps * dcomps);

            auto g_yy_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 213 * ccomps * dcomps);

            auto g_yy_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 214 * ccomps * dcomps);

            auto g_yy_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 215 * ccomps * dcomps);

            auto g_yy_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 216 * ccomps * dcomps);

            auto g_yy_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 217 * ccomps * dcomps);

            auto g_yy_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 218 * ccomps * dcomps);

            auto g_yy_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 219 * ccomps * dcomps);

            auto g_yy_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 220 * ccomps * dcomps);

            auto g_yy_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 221 * ccomps * dcomps);

            auto g_yy_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 222 * ccomps * dcomps);

            auto g_yy_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 223 * ccomps * dcomps);

            auto g_yy_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 224 * ccomps * dcomps);

            auto g_yy_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 225 * ccomps * dcomps);

            auto g_yy_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 226 * ccomps * dcomps);

            auto g_yy_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 227 * ccomps * dcomps);

            auto g_yy_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 228 * ccomps * dcomps);

            auto g_yy_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 229 * ccomps * dcomps);

            auto g_yy_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 230 * ccomps * dcomps);

            auto g_yy_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 231 * ccomps * dcomps);

            auto g_yy_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 232 * ccomps * dcomps);

            auto g_yy_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 233 * ccomps * dcomps);

            auto g_yy_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 234 * ccomps * dcomps);

            auto g_yy_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 235 * ccomps * dcomps);

            auto g_yy_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 236 * ccomps * dcomps);

            auto g_yy_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 237 * ccomps * dcomps);

            auto g_yy_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 238 * ccomps * dcomps);

            auto g_yy_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 239 * ccomps * dcomps);

            auto g_yy_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 240 * ccomps * dcomps);

            auto g_yy_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 241 * ccomps * dcomps);

            auto g_yy_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 242 * ccomps * dcomps);

            auto g_yy_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 243 * ccomps * dcomps);

            auto g_yy_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 244 * ccomps * dcomps);

            auto g_yy_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 245 * ccomps * dcomps);

            auto g_yy_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 246 * ccomps * dcomps);

            auto g_yy_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 247 * ccomps * dcomps);

            auto g_yy_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 248 * ccomps * dcomps);

            auto g_yy_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 249 * ccomps * dcomps);

            auto g_yy_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 250 * ccomps * dcomps);

            auto g_yy_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 251 * ccomps * dcomps);

            auto g_yz_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 252 * ccomps * dcomps);

            auto g_yz_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 253 * ccomps * dcomps);

            auto g_yz_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 254 * ccomps * dcomps);

            auto g_yz_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 255 * ccomps * dcomps);

            auto g_yz_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 256 * ccomps * dcomps);

            auto g_yz_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 257 * ccomps * dcomps);

            auto g_yz_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 258 * ccomps * dcomps);

            auto g_yz_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 259 * ccomps * dcomps);

            auto g_yz_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 260 * ccomps * dcomps);

            auto g_yz_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 261 * ccomps * dcomps);

            auto g_yz_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 262 * ccomps * dcomps);

            auto g_yz_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 263 * ccomps * dcomps);

            auto g_yz_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 264 * ccomps * dcomps);

            auto g_yz_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 265 * ccomps * dcomps);

            auto g_yz_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 266 * ccomps * dcomps);

            auto g_yz_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 267 * ccomps * dcomps);

            auto g_yz_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 268 * ccomps * dcomps);

            auto g_yz_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 269 * ccomps * dcomps);

            auto g_yz_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 270 * ccomps * dcomps);

            auto g_yz_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 271 * ccomps * dcomps);

            auto g_yz_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 272 * ccomps * dcomps);

            auto g_yz_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 273 * ccomps * dcomps);

            auto g_yz_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 274 * ccomps * dcomps);

            auto g_yz_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 275 * ccomps * dcomps);

            auto g_yz_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 276 * ccomps * dcomps);

            auto g_yz_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 277 * ccomps * dcomps);

            auto g_yz_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 278 * ccomps * dcomps);

            auto g_yz_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 279 * ccomps * dcomps);

            auto g_yz_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 280 * ccomps * dcomps);

            auto g_yz_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 281 * ccomps * dcomps);

            auto g_yz_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 282 * ccomps * dcomps);

            auto g_yz_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 283 * ccomps * dcomps);

            auto g_yz_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 284 * ccomps * dcomps);

            auto g_yz_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 285 * ccomps * dcomps);

            auto g_yz_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 286 * ccomps * dcomps);

            auto g_yz_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 287 * ccomps * dcomps);

            auto g_yz_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 288 * ccomps * dcomps);

            auto g_yz_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 289 * ccomps * dcomps);

            auto g_yz_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 290 * ccomps * dcomps);

            auto g_yz_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 291 * ccomps * dcomps);

            auto g_yz_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 292 * ccomps * dcomps);

            auto g_yz_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 293 * ccomps * dcomps);

            auto g_yz_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 294 * ccomps * dcomps);

            auto g_yz_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 295 * ccomps * dcomps);

            auto g_yz_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 296 * ccomps * dcomps);

            auto g_yz_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 297 * ccomps * dcomps);

            auto g_yz_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 298 * ccomps * dcomps);

            auto g_yz_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 299 * ccomps * dcomps);

            auto g_yz_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 300 * ccomps * dcomps);

            auto g_yz_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 301 * ccomps * dcomps);

            auto g_yz_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 302 * ccomps * dcomps);

            auto g_yz_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 303 * ccomps * dcomps);

            auto g_yz_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 304 * ccomps * dcomps);

            auto g_yz_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 305 * ccomps * dcomps);

            auto g_yz_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 306 * ccomps * dcomps);

            auto g_yz_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 307 * ccomps * dcomps);

            auto g_yz_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 308 * ccomps * dcomps);

            auto g_yz_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 309 * ccomps * dcomps);

            auto g_yz_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 310 * ccomps * dcomps);

            auto g_yz_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 311 * ccomps * dcomps);

            auto g_yz_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 312 * ccomps * dcomps);

            auto g_yz_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 313 * ccomps * dcomps);

            auto g_yz_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 314 * ccomps * dcomps);

            auto g_zz_0_xxxxx_x = cbuffer.data(hp_geom_20_off + 315 * ccomps * dcomps);

            auto g_zz_0_xxxxx_y = cbuffer.data(hp_geom_20_off + 316 * ccomps * dcomps);

            auto g_zz_0_xxxxx_z = cbuffer.data(hp_geom_20_off + 317 * ccomps * dcomps);

            auto g_zz_0_xxxxy_x = cbuffer.data(hp_geom_20_off + 318 * ccomps * dcomps);

            auto g_zz_0_xxxxy_y = cbuffer.data(hp_geom_20_off + 319 * ccomps * dcomps);

            auto g_zz_0_xxxxy_z = cbuffer.data(hp_geom_20_off + 320 * ccomps * dcomps);

            auto g_zz_0_xxxxz_x = cbuffer.data(hp_geom_20_off + 321 * ccomps * dcomps);

            auto g_zz_0_xxxxz_y = cbuffer.data(hp_geom_20_off + 322 * ccomps * dcomps);

            auto g_zz_0_xxxxz_z = cbuffer.data(hp_geom_20_off + 323 * ccomps * dcomps);

            auto g_zz_0_xxxyy_x = cbuffer.data(hp_geom_20_off + 324 * ccomps * dcomps);

            auto g_zz_0_xxxyy_y = cbuffer.data(hp_geom_20_off + 325 * ccomps * dcomps);

            auto g_zz_0_xxxyy_z = cbuffer.data(hp_geom_20_off + 326 * ccomps * dcomps);

            auto g_zz_0_xxxyz_x = cbuffer.data(hp_geom_20_off + 327 * ccomps * dcomps);

            auto g_zz_0_xxxyz_y = cbuffer.data(hp_geom_20_off + 328 * ccomps * dcomps);

            auto g_zz_0_xxxyz_z = cbuffer.data(hp_geom_20_off + 329 * ccomps * dcomps);

            auto g_zz_0_xxxzz_x = cbuffer.data(hp_geom_20_off + 330 * ccomps * dcomps);

            auto g_zz_0_xxxzz_y = cbuffer.data(hp_geom_20_off + 331 * ccomps * dcomps);

            auto g_zz_0_xxxzz_z = cbuffer.data(hp_geom_20_off + 332 * ccomps * dcomps);

            auto g_zz_0_xxyyy_x = cbuffer.data(hp_geom_20_off + 333 * ccomps * dcomps);

            auto g_zz_0_xxyyy_y = cbuffer.data(hp_geom_20_off + 334 * ccomps * dcomps);

            auto g_zz_0_xxyyy_z = cbuffer.data(hp_geom_20_off + 335 * ccomps * dcomps);

            auto g_zz_0_xxyyz_x = cbuffer.data(hp_geom_20_off + 336 * ccomps * dcomps);

            auto g_zz_0_xxyyz_y = cbuffer.data(hp_geom_20_off + 337 * ccomps * dcomps);

            auto g_zz_0_xxyyz_z = cbuffer.data(hp_geom_20_off + 338 * ccomps * dcomps);

            auto g_zz_0_xxyzz_x = cbuffer.data(hp_geom_20_off + 339 * ccomps * dcomps);

            auto g_zz_0_xxyzz_y = cbuffer.data(hp_geom_20_off + 340 * ccomps * dcomps);

            auto g_zz_0_xxyzz_z = cbuffer.data(hp_geom_20_off + 341 * ccomps * dcomps);

            auto g_zz_0_xxzzz_x = cbuffer.data(hp_geom_20_off + 342 * ccomps * dcomps);

            auto g_zz_0_xxzzz_y = cbuffer.data(hp_geom_20_off + 343 * ccomps * dcomps);

            auto g_zz_0_xxzzz_z = cbuffer.data(hp_geom_20_off + 344 * ccomps * dcomps);

            auto g_zz_0_xyyyy_x = cbuffer.data(hp_geom_20_off + 345 * ccomps * dcomps);

            auto g_zz_0_xyyyy_y = cbuffer.data(hp_geom_20_off + 346 * ccomps * dcomps);

            auto g_zz_0_xyyyy_z = cbuffer.data(hp_geom_20_off + 347 * ccomps * dcomps);

            auto g_zz_0_xyyyz_x = cbuffer.data(hp_geom_20_off + 348 * ccomps * dcomps);

            auto g_zz_0_xyyyz_y = cbuffer.data(hp_geom_20_off + 349 * ccomps * dcomps);

            auto g_zz_0_xyyyz_z = cbuffer.data(hp_geom_20_off + 350 * ccomps * dcomps);

            auto g_zz_0_xyyzz_x = cbuffer.data(hp_geom_20_off + 351 * ccomps * dcomps);

            auto g_zz_0_xyyzz_y = cbuffer.data(hp_geom_20_off + 352 * ccomps * dcomps);

            auto g_zz_0_xyyzz_z = cbuffer.data(hp_geom_20_off + 353 * ccomps * dcomps);

            auto g_zz_0_xyzzz_x = cbuffer.data(hp_geom_20_off + 354 * ccomps * dcomps);

            auto g_zz_0_xyzzz_y = cbuffer.data(hp_geom_20_off + 355 * ccomps * dcomps);

            auto g_zz_0_xyzzz_z = cbuffer.data(hp_geom_20_off + 356 * ccomps * dcomps);

            auto g_zz_0_xzzzz_x = cbuffer.data(hp_geom_20_off + 357 * ccomps * dcomps);

            auto g_zz_0_xzzzz_y = cbuffer.data(hp_geom_20_off + 358 * ccomps * dcomps);

            auto g_zz_0_xzzzz_z = cbuffer.data(hp_geom_20_off + 359 * ccomps * dcomps);

            auto g_zz_0_yyyyy_x = cbuffer.data(hp_geom_20_off + 360 * ccomps * dcomps);

            auto g_zz_0_yyyyy_y = cbuffer.data(hp_geom_20_off + 361 * ccomps * dcomps);

            auto g_zz_0_yyyyy_z = cbuffer.data(hp_geom_20_off + 362 * ccomps * dcomps);

            auto g_zz_0_yyyyz_x = cbuffer.data(hp_geom_20_off + 363 * ccomps * dcomps);

            auto g_zz_0_yyyyz_y = cbuffer.data(hp_geom_20_off + 364 * ccomps * dcomps);

            auto g_zz_0_yyyyz_z = cbuffer.data(hp_geom_20_off + 365 * ccomps * dcomps);

            auto g_zz_0_yyyzz_x = cbuffer.data(hp_geom_20_off + 366 * ccomps * dcomps);

            auto g_zz_0_yyyzz_y = cbuffer.data(hp_geom_20_off + 367 * ccomps * dcomps);

            auto g_zz_0_yyyzz_z = cbuffer.data(hp_geom_20_off + 368 * ccomps * dcomps);

            auto g_zz_0_yyzzz_x = cbuffer.data(hp_geom_20_off + 369 * ccomps * dcomps);

            auto g_zz_0_yyzzz_y = cbuffer.data(hp_geom_20_off + 370 * ccomps * dcomps);

            auto g_zz_0_yyzzz_z = cbuffer.data(hp_geom_20_off + 371 * ccomps * dcomps);

            auto g_zz_0_yzzzz_x = cbuffer.data(hp_geom_20_off + 372 * ccomps * dcomps);

            auto g_zz_0_yzzzz_y = cbuffer.data(hp_geom_20_off + 373 * ccomps * dcomps);

            auto g_zz_0_yzzzz_z = cbuffer.data(hp_geom_20_off + 374 * ccomps * dcomps);

            auto g_zz_0_zzzzz_x = cbuffer.data(hp_geom_20_off + 375 * ccomps * dcomps);

            auto g_zz_0_zzzzz_y = cbuffer.data(hp_geom_20_off + 376 * ccomps * dcomps);

            auto g_zz_0_zzzzz_z = cbuffer.data(hp_geom_20_off + 377 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HDSS

            const auto hd_geom_20_off = idx_geom_20_hdxx + i * dcomps + j;

            auto g_xx_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 2 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 5 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 8 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 11 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 14 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 17 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 20 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 23 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 26 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 29 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 32 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 35 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 38 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 41 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 44 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 47 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 50 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 53 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 56 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 59 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 62 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 65 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 68 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 71 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 74 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 77 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 80 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 83 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 84 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 85 * ccomps * dcomps);

            auto g_xx_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 86 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 87 * ccomps * dcomps);

            auto g_xx_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 88 * ccomps * dcomps);

            auto g_xx_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 89 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 90 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 91 * ccomps * dcomps);

            auto g_xx_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 92 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 93 * ccomps * dcomps);

            auto g_xx_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 94 * ccomps * dcomps);

            auto g_xx_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 95 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 96 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 97 * ccomps * dcomps);

            auto g_xx_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 98 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 99 * ccomps * dcomps);

            auto g_xx_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 100 * ccomps * dcomps);

            auto g_xx_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 101 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 102 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 103 * ccomps * dcomps);

            auto g_xx_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 104 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 105 * ccomps * dcomps);

            auto g_xx_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 106 * ccomps * dcomps);

            auto g_xx_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 107 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 108 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 109 * ccomps * dcomps);

            auto g_xx_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 110 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 111 * ccomps * dcomps);

            auto g_xx_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 112 * ccomps * dcomps);

            auto g_xx_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 113 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 114 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 115 * ccomps * dcomps);

            auto g_xx_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 116 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 117 * ccomps * dcomps);

            auto g_xx_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 118 * ccomps * dcomps);

            auto g_xx_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 119 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 120 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 121 * ccomps * dcomps);

            auto g_xx_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 122 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 123 * ccomps * dcomps);

            auto g_xx_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 124 * ccomps * dcomps);

            auto g_xx_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 125 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 128 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 131 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 134 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 137 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 140 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 143 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 146 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 149 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 152 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 155 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 158 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 161 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 164 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 167 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 168 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 169 * ccomps * dcomps);

            auto g_xy_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 170 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 171 * ccomps * dcomps);

            auto g_xy_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 172 * ccomps * dcomps);

            auto g_xy_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 173 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 174 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 175 * ccomps * dcomps);

            auto g_xy_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 176 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 177 * ccomps * dcomps);

            auto g_xy_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 178 * ccomps * dcomps);

            auto g_xy_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 179 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 180 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 181 * ccomps * dcomps);

            auto g_xy_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 182 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 183 * ccomps * dcomps);

            auto g_xy_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 184 * ccomps * dcomps);

            auto g_xy_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 185 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 186 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 187 * ccomps * dcomps);

            auto g_xy_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 188 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 189 * ccomps * dcomps);

            auto g_xy_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 190 * ccomps * dcomps);

            auto g_xy_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 191 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 192 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 193 * ccomps * dcomps);

            auto g_xy_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 194 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 195 * ccomps * dcomps);

            auto g_xy_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 196 * ccomps * dcomps);

            auto g_xy_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 197 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 198 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 199 * ccomps * dcomps);

            auto g_xy_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 200 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 201 * ccomps * dcomps);

            auto g_xy_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 202 * ccomps * dcomps);

            auto g_xy_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 203 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 204 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 205 * ccomps * dcomps);

            auto g_xy_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 206 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 207 * ccomps * dcomps);

            auto g_xy_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 208 * ccomps * dcomps);

            auto g_xy_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 209 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 210 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 211 * ccomps * dcomps);

            auto g_xy_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 212 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 213 * ccomps * dcomps);

            auto g_xy_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 214 * ccomps * dcomps);

            auto g_xy_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 215 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 216 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 217 * ccomps * dcomps);

            auto g_xy_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 218 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 219 * ccomps * dcomps);

            auto g_xy_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 220 * ccomps * dcomps);

            auto g_xy_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 221 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 222 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 223 * ccomps * dcomps);

            auto g_xy_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 224 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 225 * ccomps * dcomps);

            auto g_xy_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 226 * ccomps * dcomps);

            auto g_xy_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 227 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 228 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 229 * ccomps * dcomps);

            auto g_xy_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 230 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 231 * ccomps * dcomps);

            auto g_xy_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 232 * ccomps * dcomps);

            auto g_xy_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 233 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 234 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 235 * ccomps * dcomps);

            auto g_xy_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 236 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 237 * ccomps * dcomps);

            auto g_xy_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 238 * ccomps * dcomps);

            auto g_xy_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 239 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 240 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 241 * ccomps * dcomps);

            auto g_xy_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 242 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 243 * ccomps * dcomps);

            auto g_xy_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 244 * ccomps * dcomps);

            auto g_xy_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 245 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 246 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 247 * ccomps * dcomps);

            auto g_xy_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 248 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 249 * ccomps * dcomps);

            auto g_xy_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 250 * ccomps * dcomps);

            auto g_xy_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 251 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 252 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 253 * ccomps * dcomps);

            auto g_xz_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 254 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 255 * ccomps * dcomps);

            auto g_xz_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 256 * ccomps * dcomps);

            auto g_xz_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 257 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 258 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 259 * ccomps * dcomps);

            auto g_xz_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 260 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 261 * ccomps * dcomps);

            auto g_xz_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 262 * ccomps * dcomps);

            auto g_xz_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 263 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 264 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 265 * ccomps * dcomps);

            auto g_xz_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 266 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 267 * ccomps * dcomps);

            auto g_xz_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 268 * ccomps * dcomps);

            auto g_xz_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 269 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 270 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 271 * ccomps * dcomps);

            auto g_xz_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 272 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 273 * ccomps * dcomps);

            auto g_xz_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 274 * ccomps * dcomps);

            auto g_xz_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 275 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 276 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 277 * ccomps * dcomps);

            auto g_xz_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 278 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 279 * ccomps * dcomps);

            auto g_xz_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 280 * ccomps * dcomps);

            auto g_xz_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 281 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 282 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 283 * ccomps * dcomps);

            auto g_xz_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 284 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 285 * ccomps * dcomps);

            auto g_xz_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 286 * ccomps * dcomps);

            auto g_xz_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 287 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 288 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 289 * ccomps * dcomps);

            auto g_xz_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 290 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 291 * ccomps * dcomps);

            auto g_xz_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 292 * ccomps * dcomps);

            auto g_xz_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 293 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 294 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 295 * ccomps * dcomps);

            auto g_xz_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 296 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 297 * ccomps * dcomps);

            auto g_xz_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 298 * ccomps * dcomps);

            auto g_xz_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 299 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 300 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 301 * ccomps * dcomps);

            auto g_xz_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 302 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 303 * ccomps * dcomps);

            auto g_xz_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 304 * ccomps * dcomps);

            auto g_xz_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 305 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 306 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 307 * ccomps * dcomps);

            auto g_xz_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 308 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 309 * ccomps * dcomps);

            auto g_xz_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 310 * ccomps * dcomps);

            auto g_xz_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 311 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 312 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 313 * ccomps * dcomps);

            auto g_xz_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 314 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 315 * ccomps * dcomps);

            auto g_xz_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 316 * ccomps * dcomps);

            auto g_xz_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 317 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 318 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 319 * ccomps * dcomps);

            auto g_xz_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 320 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 321 * ccomps * dcomps);

            auto g_xz_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 322 * ccomps * dcomps);

            auto g_xz_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 323 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 324 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 325 * ccomps * dcomps);

            auto g_xz_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 326 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 327 * ccomps * dcomps);

            auto g_xz_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 328 * ccomps * dcomps);

            auto g_xz_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 329 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 330 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 331 * ccomps * dcomps);

            auto g_xz_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 332 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 333 * ccomps * dcomps);

            auto g_xz_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 334 * ccomps * dcomps);

            auto g_xz_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 335 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 336 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 337 * ccomps * dcomps);

            auto g_xz_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 338 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 339 * ccomps * dcomps);

            auto g_xz_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 340 * ccomps * dcomps);

            auto g_xz_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 341 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 342 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 343 * ccomps * dcomps);

            auto g_xz_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 344 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 345 * ccomps * dcomps);

            auto g_xz_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 346 * ccomps * dcomps);

            auto g_xz_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 347 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 348 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 349 * ccomps * dcomps);

            auto g_xz_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 350 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 351 * ccomps * dcomps);

            auto g_xz_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 352 * ccomps * dcomps);

            auto g_xz_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 353 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 354 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 355 * ccomps * dcomps);

            auto g_xz_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 356 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 357 * ccomps * dcomps);

            auto g_xz_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 358 * ccomps * dcomps);

            auto g_xz_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 359 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 360 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 361 * ccomps * dcomps);

            auto g_xz_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 362 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 363 * ccomps * dcomps);

            auto g_xz_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 364 * ccomps * dcomps);

            auto g_xz_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 365 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 366 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 367 * ccomps * dcomps);

            auto g_xz_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 368 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 369 * ccomps * dcomps);

            auto g_xz_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 370 * ccomps * dcomps);

            auto g_xz_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 371 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 372 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 373 * ccomps * dcomps);

            auto g_xz_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 374 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 375 * ccomps * dcomps);

            auto g_xz_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 376 * ccomps * dcomps);

            auto g_xz_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 377 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 378 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 379 * ccomps * dcomps);

            auto g_yy_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 380 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 381 * ccomps * dcomps);

            auto g_yy_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 382 * ccomps * dcomps);

            auto g_yy_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 383 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 384 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 385 * ccomps * dcomps);

            auto g_yy_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 386 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 387 * ccomps * dcomps);

            auto g_yy_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 388 * ccomps * dcomps);

            auto g_yy_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 389 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 390 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 391 * ccomps * dcomps);

            auto g_yy_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 392 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 393 * ccomps * dcomps);

            auto g_yy_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 394 * ccomps * dcomps);

            auto g_yy_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 395 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 396 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 397 * ccomps * dcomps);

            auto g_yy_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 398 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 399 * ccomps * dcomps);

            auto g_yy_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 400 * ccomps * dcomps);

            auto g_yy_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 401 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 402 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 403 * ccomps * dcomps);

            auto g_yy_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 404 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 405 * ccomps * dcomps);

            auto g_yy_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 406 * ccomps * dcomps);

            auto g_yy_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 407 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 408 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 409 * ccomps * dcomps);

            auto g_yy_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 410 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 411 * ccomps * dcomps);

            auto g_yy_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 412 * ccomps * dcomps);

            auto g_yy_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 413 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 414 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 415 * ccomps * dcomps);

            auto g_yy_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 416 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 417 * ccomps * dcomps);

            auto g_yy_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 418 * ccomps * dcomps);

            auto g_yy_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 419 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 420 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 421 * ccomps * dcomps);

            auto g_yy_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 422 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 423 * ccomps * dcomps);

            auto g_yy_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 424 * ccomps * dcomps);

            auto g_yy_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 425 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 426 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 427 * ccomps * dcomps);

            auto g_yy_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 428 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 429 * ccomps * dcomps);

            auto g_yy_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 430 * ccomps * dcomps);

            auto g_yy_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 431 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 432 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 433 * ccomps * dcomps);

            auto g_yy_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 434 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 435 * ccomps * dcomps);

            auto g_yy_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 436 * ccomps * dcomps);

            auto g_yy_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 437 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 438 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 439 * ccomps * dcomps);

            auto g_yy_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 440 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 441 * ccomps * dcomps);

            auto g_yy_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 442 * ccomps * dcomps);

            auto g_yy_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 443 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 444 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 445 * ccomps * dcomps);

            auto g_yy_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 446 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 447 * ccomps * dcomps);

            auto g_yy_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 448 * ccomps * dcomps);

            auto g_yy_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 449 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 450 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 451 * ccomps * dcomps);

            auto g_yy_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 452 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 453 * ccomps * dcomps);

            auto g_yy_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 454 * ccomps * dcomps);

            auto g_yy_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 455 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 456 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 457 * ccomps * dcomps);

            auto g_yy_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 458 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 459 * ccomps * dcomps);

            auto g_yy_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 460 * ccomps * dcomps);

            auto g_yy_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 461 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 462 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 463 * ccomps * dcomps);

            auto g_yy_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 464 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 465 * ccomps * dcomps);

            auto g_yy_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 466 * ccomps * dcomps);

            auto g_yy_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 467 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 468 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 469 * ccomps * dcomps);

            auto g_yy_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 470 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 471 * ccomps * dcomps);

            auto g_yy_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 472 * ccomps * dcomps);

            auto g_yy_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 473 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 474 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 475 * ccomps * dcomps);

            auto g_yy_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 476 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 477 * ccomps * dcomps);

            auto g_yy_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 478 * ccomps * dcomps);

            auto g_yy_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 479 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 480 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 481 * ccomps * dcomps);

            auto g_yy_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 482 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 483 * ccomps * dcomps);

            auto g_yy_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 484 * ccomps * dcomps);

            auto g_yy_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 485 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 486 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 487 * ccomps * dcomps);

            auto g_yy_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 488 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 489 * ccomps * dcomps);

            auto g_yy_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 490 * ccomps * dcomps);

            auto g_yy_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 491 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 492 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 493 * ccomps * dcomps);

            auto g_yy_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 494 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 495 * ccomps * dcomps);

            auto g_yy_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 496 * ccomps * dcomps);

            auto g_yy_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 497 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 498 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 499 * ccomps * dcomps);

            auto g_yy_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 500 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 501 * ccomps * dcomps);

            auto g_yy_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 502 * ccomps * dcomps);

            auto g_yy_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 503 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 504 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 505 * ccomps * dcomps);

            auto g_yz_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 506 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 507 * ccomps * dcomps);

            auto g_yz_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 508 * ccomps * dcomps);

            auto g_yz_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 509 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 510 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 511 * ccomps * dcomps);

            auto g_yz_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 512 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 513 * ccomps * dcomps);

            auto g_yz_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 514 * ccomps * dcomps);

            auto g_yz_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 515 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 516 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 517 * ccomps * dcomps);

            auto g_yz_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 518 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 519 * ccomps * dcomps);

            auto g_yz_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 520 * ccomps * dcomps);

            auto g_yz_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 521 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 522 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 523 * ccomps * dcomps);

            auto g_yz_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 524 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 525 * ccomps * dcomps);

            auto g_yz_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 526 * ccomps * dcomps);

            auto g_yz_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 527 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 528 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 529 * ccomps * dcomps);

            auto g_yz_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 530 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 531 * ccomps * dcomps);

            auto g_yz_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 532 * ccomps * dcomps);

            auto g_yz_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 533 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 534 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 535 * ccomps * dcomps);

            auto g_yz_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 536 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 537 * ccomps * dcomps);

            auto g_yz_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 538 * ccomps * dcomps);

            auto g_yz_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 539 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 540 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 541 * ccomps * dcomps);

            auto g_yz_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 542 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 543 * ccomps * dcomps);

            auto g_yz_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 544 * ccomps * dcomps);

            auto g_yz_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 545 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 546 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 547 * ccomps * dcomps);

            auto g_yz_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 548 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 549 * ccomps * dcomps);

            auto g_yz_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 550 * ccomps * dcomps);

            auto g_yz_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 551 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 552 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 553 * ccomps * dcomps);

            auto g_yz_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 554 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 555 * ccomps * dcomps);

            auto g_yz_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 556 * ccomps * dcomps);

            auto g_yz_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 557 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 558 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 559 * ccomps * dcomps);

            auto g_yz_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 560 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 561 * ccomps * dcomps);

            auto g_yz_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 562 * ccomps * dcomps);

            auto g_yz_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 563 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 564 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 565 * ccomps * dcomps);

            auto g_yz_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 566 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 567 * ccomps * dcomps);

            auto g_yz_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 568 * ccomps * dcomps);

            auto g_yz_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 569 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 570 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 571 * ccomps * dcomps);

            auto g_yz_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 572 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 573 * ccomps * dcomps);

            auto g_yz_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 574 * ccomps * dcomps);

            auto g_yz_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 575 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 576 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 577 * ccomps * dcomps);

            auto g_yz_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 578 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 579 * ccomps * dcomps);

            auto g_yz_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 580 * ccomps * dcomps);

            auto g_yz_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 581 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 582 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 583 * ccomps * dcomps);

            auto g_yz_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 584 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 585 * ccomps * dcomps);

            auto g_yz_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 586 * ccomps * dcomps);

            auto g_yz_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 587 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 588 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 589 * ccomps * dcomps);

            auto g_yz_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 590 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 591 * ccomps * dcomps);

            auto g_yz_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 592 * ccomps * dcomps);

            auto g_yz_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 593 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 594 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 595 * ccomps * dcomps);

            auto g_yz_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 596 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 597 * ccomps * dcomps);

            auto g_yz_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 598 * ccomps * dcomps);

            auto g_yz_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 599 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 600 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 601 * ccomps * dcomps);

            auto g_yz_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 602 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 603 * ccomps * dcomps);

            auto g_yz_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 604 * ccomps * dcomps);

            auto g_yz_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 605 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 606 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 607 * ccomps * dcomps);

            auto g_yz_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 608 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 609 * ccomps * dcomps);

            auto g_yz_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 610 * ccomps * dcomps);

            auto g_yz_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 611 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 612 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 613 * ccomps * dcomps);

            auto g_yz_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 614 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 615 * ccomps * dcomps);

            auto g_yz_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 616 * ccomps * dcomps);

            auto g_yz_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 617 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 618 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 619 * ccomps * dcomps);

            auto g_yz_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 620 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 621 * ccomps * dcomps);

            auto g_yz_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 622 * ccomps * dcomps);

            auto g_yz_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 623 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 624 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 625 * ccomps * dcomps);

            auto g_yz_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 626 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 627 * ccomps * dcomps);

            auto g_yz_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 628 * ccomps * dcomps);

            auto g_yz_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 629 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xx = cbuffer.data(hd_geom_20_off + 630 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xy = cbuffer.data(hd_geom_20_off + 631 * ccomps * dcomps);

            auto g_zz_0_xxxxx_xz = cbuffer.data(hd_geom_20_off + 632 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yy = cbuffer.data(hd_geom_20_off + 633 * ccomps * dcomps);

            auto g_zz_0_xxxxx_yz = cbuffer.data(hd_geom_20_off + 634 * ccomps * dcomps);

            auto g_zz_0_xxxxx_zz = cbuffer.data(hd_geom_20_off + 635 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xx = cbuffer.data(hd_geom_20_off + 636 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xy = cbuffer.data(hd_geom_20_off + 637 * ccomps * dcomps);

            auto g_zz_0_xxxxy_xz = cbuffer.data(hd_geom_20_off + 638 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yy = cbuffer.data(hd_geom_20_off + 639 * ccomps * dcomps);

            auto g_zz_0_xxxxy_yz = cbuffer.data(hd_geom_20_off + 640 * ccomps * dcomps);

            auto g_zz_0_xxxxy_zz = cbuffer.data(hd_geom_20_off + 641 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xx = cbuffer.data(hd_geom_20_off + 642 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xy = cbuffer.data(hd_geom_20_off + 643 * ccomps * dcomps);

            auto g_zz_0_xxxxz_xz = cbuffer.data(hd_geom_20_off + 644 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yy = cbuffer.data(hd_geom_20_off + 645 * ccomps * dcomps);

            auto g_zz_0_xxxxz_yz = cbuffer.data(hd_geom_20_off + 646 * ccomps * dcomps);

            auto g_zz_0_xxxxz_zz = cbuffer.data(hd_geom_20_off + 647 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xx = cbuffer.data(hd_geom_20_off + 648 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xy = cbuffer.data(hd_geom_20_off + 649 * ccomps * dcomps);

            auto g_zz_0_xxxyy_xz = cbuffer.data(hd_geom_20_off + 650 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yy = cbuffer.data(hd_geom_20_off + 651 * ccomps * dcomps);

            auto g_zz_0_xxxyy_yz = cbuffer.data(hd_geom_20_off + 652 * ccomps * dcomps);

            auto g_zz_0_xxxyy_zz = cbuffer.data(hd_geom_20_off + 653 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xx = cbuffer.data(hd_geom_20_off + 654 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xy = cbuffer.data(hd_geom_20_off + 655 * ccomps * dcomps);

            auto g_zz_0_xxxyz_xz = cbuffer.data(hd_geom_20_off + 656 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yy = cbuffer.data(hd_geom_20_off + 657 * ccomps * dcomps);

            auto g_zz_0_xxxyz_yz = cbuffer.data(hd_geom_20_off + 658 * ccomps * dcomps);

            auto g_zz_0_xxxyz_zz = cbuffer.data(hd_geom_20_off + 659 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xx = cbuffer.data(hd_geom_20_off + 660 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xy = cbuffer.data(hd_geom_20_off + 661 * ccomps * dcomps);

            auto g_zz_0_xxxzz_xz = cbuffer.data(hd_geom_20_off + 662 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yy = cbuffer.data(hd_geom_20_off + 663 * ccomps * dcomps);

            auto g_zz_0_xxxzz_yz = cbuffer.data(hd_geom_20_off + 664 * ccomps * dcomps);

            auto g_zz_0_xxxzz_zz = cbuffer.data(hd_geom_20_off + 665 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xx = cbuffer.data(hd_geom_20_off + 666 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xy = cbuffer.data(hd_geom_20_off + 667 * ccomps * dcomps);

            auto g_zz_0_xxyyy_xz = cbuffer.data(hd_geom_20_off + 668 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yy = cbuffer.data(hd_geom_20_off + 669 * ccomps * dcomps);

            auto g_zz_0_xxyyy_yz = cbuffer.data(hd_geom_20_off + 670 * ccomps * dcomps);

            auto g_zz_0_xxyyy_zz = cbuffer.data(hd_geom_20_off + 671 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xx = cbuffer.data(hd_geom_20_off + 672 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xy = cbuffer.data(hd_geom_20_off + 673 * ccomps * dcomps);

            auto g_zz_0_xxyyz_xz = cbuffer.data(hd_geom_20_off + 674 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yy = cbuffer.data(hd_geom_20_off + 675 * ccomps * dcomps);

            auto g_zz_0_xxyyz_yz = cbuffer.data(hd_geom_20_off + 676 * ccomps * dcomps);

            auto g_zz_0_xxyyz_zz = cbuffer.data(hd_geom_20_off + 677 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xx = cbuffer.data(hd_geom_20_off + 678 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xy = cbuffer.data(hd_geom_20_off + 679 * ccomps * dcomps);

            auto g_zz_0_xxyzz_xz = cbuffer.data(hd_geom_20_off + 680 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yy = cbuffer.data(hd_geom_20_off + 681 * ccomps * dcomps);

            auto g_zz_0_xxyzz_yz = cbuffer.data(hd_geom_20_off + 682 * ccomps * dcomps);

            auto g_zz_0_xxyzz_zz = cbuffer.data(hd_geom_20_off + 683 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xx = cbuffer.data(hd_geom_20_off + 684 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xy = cbuffer.data(hd_geom_20_off + 685 * ccomps * dcomps);

            auto g_zz_0_xxzzz_xz = cbuffer.data(hd_geom_20_off + 686 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yy = cbuffer.data(hd_geom_20_off + 687 * ccomps * dcomps);

            auto g_zz_0_xxzzz_yz = cbuffer.data(hd_geom_20_off + 688 * ccomps * dcomps);

            auto g_zz_0_xxzzz_zz = cbuffer.data(hd_geom_20_off + 689 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xx = cbuffer.data(hd_geom_20_off + 690 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xy = cbuffer.data(hd_geom_20_off + 691 * ccomps * dcomps);

            auto g_zz_0_xyyyy_xz = cbuffer.data(hd_geom_20_off + 692 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yy = cbuffer.data(hd_geom_20_off + 693 * ccomps * dcomps);

            auto g_zz_0_xyyyy_yz = cbuffer.data(hd_geom_20_off + 694 * ccomps * dcomps);

            auto g_zz_0_xyyyy_zz = cbuffer.data(hd_geom_20_off + 695 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xx = cbuffer.data(hd_geom_20_off + 696 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xy = cbuffer.data(hd_geom_20_off + 697 * ccomps * dcomps);

            auto g_zz_0_xyyyz_xz = cbuffer.data(hd_geom_20_off + 698 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yy = cbuffer.data(hd_geom_20_off + 699 * ccomps * dcomps);

            auto g_zz_0_xyyyz_yz = cbuffer.data(hd_geom_20_off + 700 * ccomps * dcomps);

            auto g_zz_0_xyyyz_zz = cbuffer.data(hd_geom_20_off + 701 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xx = cbuffer.data(hd_geom_20_off + 702 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xy = cbuffer.data(hd_geom_20_off + 703 * ccomps * dcomps);

            auto g_zz_0_xyyzz_xz = cbuffer.data(hd_geom_20_off + 704 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yy = cbuffer.data(hd_geom_20_off + 705 * ccomps * dcomps);

            auto g_zz_0_xyyzz_yz = cbuffer.data(hd_geom_20_off + 706 * ccomps * dcomps);

            auto g_zz_0_xyyzz_zz = cbuffer.data(hd_geom_20_off + 707 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xx = cbuffer.data(hd_geom_20_off + 708 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xy = cbuffer.data(hd_geom_20_off + 709 * ccomps * dcomps);

            auto g_zz_0_xyzzz_xz = cbuffer.data(hd_geom_20_off + 710 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yy = cbuffer.data(hd_geom_20_off + 711 * ccomps * dcomps);

            auto g_zz_0_xyzzz_yz = cbuffer.data(hd_geom_20_off + 712 * ccomps * dcomps);

            auto g_zz_0_xyzzz_zz = cbuffer.data(hd_geom_20_off + 713 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xx = cbuffer.data(hd_geom_20_off + 714 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xy = cbuffer.data(hd_geom_20_off + 715 * ccomps * dcomps);

            auto g_zz_0_xzzzz_xz = cbuffer.data(hd_geom_20_off + 716 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yy = cbuffer.data(hd_geom_20_off + 717 * ccomps * dcomps);

            auto g_zz_0_xzzzz_yz = cbuffer.data(hd_geom_20_off + 718 * ccomps * dcomps);

            auto g_zz_0_xzzzz_zz = cbuffer.data(hd_geom_20_off + 719 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xx = cbuffer.data(hd_geom_20_off + 720 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xy = cbuffer.data(hd_geom_20_off + 721 * ccomps * dcomps);

            auto g_zz_0_yyyyy_xz = cbuffer.data(hd_geom_20_off + 722 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yy = cbuffer.data(hd_geom_20_off + 723 * ccomps * dcomps);

            auto g_zz_0_yyyyy_yz = cbuffer.data(hd_geom_20_off + 724 * ccomps * dcomps);

            auto g_zz_0_yyyyy_zz = cbuffer.data(hd_geom_20_off + 725 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xx = cbuffer.data(hd_geom_20_off + 726 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xy = cbuffer.data(hd_geom_20_off + 727 * ccomps * dcomps);

            auto g_zz_0_yyyyz_xz = cbuffer.data(hd_geom_20_off + 728 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yy = cbuffer.data(hd_geom_20_off + 729 * ccomps * dcomps);

            auto g_zz_0_yyyyz_yz = cbuffer.data(hd_geom_20_off + 730 * ccomps * dcomps);

            auto g_zz_0_yyyyz_zz = cbuffer.data(hd_geom_20_off + 731 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xx = cbuffer.data(hd_geom_20_off + 732 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xy = cbuffer.data(hd_geom_20_off + 733 * ccomps * dcomps);

            auto g_zz_0_yyyzz_xz = cbuffer.data(hd_geom_20_off + 734 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yy = cbuffer.data(hd_geom_20_off + 735 * ccomps * dcomps);

            auto g_zz_0_yyyzz_yz = cbuffer.data(hd_geom_20_off + 736 * ccomps * dcomps);

            auto g_zz_0_yyyzz_zz = cbuffer.data(hd_geom_20_off + 737 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xx = cbuffer.data(hd_geom_20_off + 738 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xy = cbuffer.data(hd_geom_20_off + 739 * ccomps * dcomps);

            auto g_zz_0_yyzzz_xz = cbuffer.data(hd_geom_20_off + 740 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yy = cbuffer.data(hd_geom_20_off + 741 * ccomps * dcomps);

            auto g_zz_0_yyzzz_yz = cbuffer.data(hd_geom_20_off + 742 * ccomps * dcomps);

            auto g_zz_0_yyzzz_zz = cbuffer.data(hd_geom_20_off + 743 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xx = cbuffer.data(hd_geom_20_off + 744 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xy = cbuffer.data(hd_geom_20_off + 745 * ccomps * dcomps);

            auto g_zz_0_yzzzz_xz = cbuffer.data(hd_geom_20_off + 746 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yy = cbuffer.data(hd_geom_20_off + 747 * ccomps * dcomps);

            auto g_zz_0_yzzzz_yz = cbuffer.data(hd_geom_20_off + 748 * ccomps * dcomps);

            auto g_zz_0_yzzzz_zz = cbuffer.data(hd_geom_20_off + 749 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xx = cbuffer.data(hd_geom_20_off + 750 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xy = cbuffer.data(hd_geom_20_off + 751 * ccomps * dcomps);

            auto g_zz_0_zzzzz_xz = cbuffer.data(hd_geom_20_off + 752 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yy = cbuffer.data(hd_geom_20_off + 753 * ccomps * dcomps);

            auto g_zz_0_zzzzz_yz = cbuffer.data(hd_geom_20_off + 754 * ccomps * dcomps);

            auto g_zz_0_zzzzz_zz = cbuffer.data(hd_geom_20_off + 755 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ipxx

            const auto ip_geom_20_off = idx_geom_20_ipxx + i * dcomps + j;

            /// Set up 0-3 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxxx_x = cbuffer.data(ip_geom_20_off + 0 * ccomps * dcomps);

            auto g_xx_0_xxxxxx_y = cbuffer.data(ip_geom_20_off + 1 * ccomps * dcomps);

            auto g_xx_0_xxxxxx_z = cbuffer.data(ip_geom_20_off + 2 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_x, g_x_0_xxxxx_y, g_x_0_xxxxx_z, g_xx_0_xxxxx_x, g_xx_0_xxxxx_xx, g_xx_0_xxxxx_xy, g_xx_0_xxxxx_xz, g_xx_0_xxxxx_y, g_xx_0_xxxxx_z, g_xx_0_xxxxxx_x, g_xx_0_xxxxxx_y, g_xx_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxxx_x[k] = -2.0 * g_x_0_xxxxx_x[k] - g_xx_0_xxxxx_x[k] * ab_x + g_xx_0_xxxxx_xx[k];

                g_xx_0_xxxxxx_y[k] = -2.0 * g_x_0_xxxxx_y[k] - g_xx_0_xxxxx_y[k] * ab_x + g_xx_0_xxxxx_xy[k];

                g_xx_0_xxxxxx_z[k] = -2.0 * g_x_0_xxxxx_z[k] - g_xx_0_xxxxx_z[k] * ab_x + g_xx_0_xxxxx_xz[k];
            }

            /// Set up 3-6 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxxy_x = cbuffer.data(ip_geom_20_off + 3 * ccomps * dcomps);

            auto g_xx_0_xxxxxy_y = cbuffer.data(ip_geom_20_off + 4 * ccomps * dcomps);

            auto g_xx_0_xxxxxy_z = cbuffer.data(ip_geom_20_off + 5 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxx_x, g_xx_0_xxxxx_xy, g_xx_0_xxxxx_y, g_xx_0_xxxxx_yy, g_xx_0_xxxxx_yz, g_xx_0_xxxxx_z, g_xx_0_xxxxxy_x, g_xx_0_xxxxxy_y, g_xx_0_xxxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxxy_x[k] = -g_xx_0_xxxxx_x[k] * ab_y + g_xx_0_xxxxx_xy[k];

                g_xx_0_xxxxxy_y[k] = -g_xx_0_xxxxx_y[k] * ab_y + g_xx_0_xxxxx_yy[k];

                g_xx_0_xxxxxy_z[k] = -g_xx_0_xxxxx_z[k] * ab_y + g_xx_0_xxxxx_yz[k];
            }

            /// Set up 6-9 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxxz_x = cbuffer.data(ip_geom_20_off + 6 * ccomps * dcomps);

            auto g_xx_0_xxxxxz_y = cbuffer.data(ip_geom_20_off + 7 * ccomps * dcomps);

            auto g_xx_0_xxxxxz_z = cbuffer.data(ip_geom_20_off + 8 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxx_x, g_xx_0_xxxxx_xz, g_xx_0_xxxxx_y, g_xx_0_xxxxx_yz, g_xx_0_xxxxx_z, g_xx_0_xxxxx_zz, g_xx_0_xxxxxz_x, g_xx_0_xxxxxz_y, g_xx_0_xxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxxz_x[k] = -g_xx_0_xxxxx_x[k] * ab_z + g_xx_0_xxxxx_xz[k];

                g_xx_0_xxxxxz_y[k] = -g_xx_0_xxxxx_y[k] * ab_z + g_xx_0_xxxxx_yz[k];

                g_xx_0_xxxxxz_z[k] = -g_xx_0_xxxxx_z[k] * ab_z + g_xx_0_xxxxx_zz[k];
            }

            /// Set up 9-12 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxyy_x = cbuffer.data(ip_geom_20_off + 9 * ccomps * dcomps);

            auto g_xx_0_xxxxyy_y = cbuffer.data(ip_geom_20_off + 10 * ccomps * dcomps);

            auto g_xx_0_xxxxyy_z = cbuffer.data(ip_geom_20_off + 11 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxy_x, g_xx_0_xxxxy_xy, g_xx_0_xxxxy_y, g_xx_0_xxxxy_yy, g_xx_0_xxxxy_yz, g_xx_0_xxxxy_z, g_xx_0_xxxxyy_x, g_xx_0_xxxxyy_y, g_xx_0_xxxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxyy_x[k] = -g_xx_0_xxxxy_x[k] * ab_y + g_xx_0_xxxxy_xy[k];

                g_xx_0_xxxxyy_y[k] = -g_xx_0_xxxxy_y[k] * ab_y + g_xx_0_xxxxy_yy[k];

                g_xx_0_xxxxyy_z[k] = -g_xx_0_xxxxy_z[k] * ab_y + g_xx_0_xxxxy_yz[k];
            }

            /// Set up 12-15 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxyz_x = cbuffer.data(ip_geom_20_off + 12 * ccomps * dcomps);

            auto g_xx_0_xxxxyz_y = cbuffer.data(ip_geom_20_off + 13 * ccomps * dcomps);

            auto g_xx_0_xxxxyz_z = cbuffer.data(ip_geom_20_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxyz_x, g_xx_0_xxxxyz_y, g_xx_0_xxxxyz_z, g_xx_0_xxxxz_x, g_xx_0_xxxxz_xy, g_xx_0_xxxxz_y, g_xx_0_xxxxz_yy, g_xx_0_xxxxz_yz, g_xx_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxyz_x[k] = -g_xx_0_xxxxz_x[k] * ab_y + g_xx_0_xxxxz_xy[k];

                g_xx_0_xxxxyz_y[k] = -g_xx_0_xxxxz_y[k] * ab_y + g_xx_0_xxxxz_yy[k];

                g_xx_0_xxxxyz_z[k] = -g_xx_0_xxxxz_z[k] * ab_y + g_xx_0_xxxxz_yz[k];
            }

            /// Set up 15-18 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxxzz_x = cbuffer.data(ip_geom_20_off + 15 * ccomps * dcomps);

            auto g_xx_0_xxxxzz_y = cbuffer.data(ip_geom_20_off + 16 * ccomps * dcomps);

            auto g_xx_0_xxxxzz_z = cbuffer.data(ip_geom_20_off + 17 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxxz_x, g_xx_0_xxxxz_xz, g_xx_0_xxxxz_y, g_xx_0_xxxxz_yz, g_xx_0_xxxxz_z, g_xx_0_xxxxz_zz, g_xx_0_xxxxzz_x, g_xx_0_xxxxzz_y, g_xx_0_xxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxxzz_x[k] = -g_xx_0_xxxxz_x[k] * ab_z + g_xx_0_xxxxz_xz[k];

                g_xx_0_xxxxzz_y[k] = -g_xx_0_xxxxz_y[k] * ab_z + g_xx_0_xxxxz_yz[k];

                g_xx_0_xxxxzz_z[k] = -g_xx_0_xxxxz_z[k] * ab_z + g_xx_0_xxxxz_zz[k];
            }

            /// Set up 18-21 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyyy_x = cbuffer.data(ip_geom_20_off + 18 * ccomps * dcomps);

            auto g_xx_0_xxxyyy_y = cbuffer.data(ip_geom_20_off + 19 * ccomps * dcomps);

            auto g_xx_0_xxxyyy_z = cbuffer.data(ip_geom_20_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyy_x, g_xx_0_xxxyy_xy, g_xx_0_xxxyy_y, g_xx_0_xxxyy_yy, g_xx_0_xxxyy_yz, g_xx_0_xxxyy_z, g_xx_0_xxxyyy_x, g_xx_0_xxxyyy_y, g_xx_0_xxxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyyy_x[k] = -g_xx_0_xxxyy_x[k] * ab_y + g_xx_0_xxxyy_xy[k];

                g_xx_0_xxxyyy_y[k] = -g_xx_0_xxxyy_y[k] * ab_y + g_xx_0_xxxyy_yy[k];

                g_xx_0_xxxyyy_z[k] = -g_xx_0_xxxyy_z[k] * ab_y + g_xx_0_xxxyy_yz[k];
            }

            /// Set up 21-24 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyyz_x = cbuffer.data(ip_geom_20_off + 21 * ccomps * dcomps);

            auto g_xx_0_xxxyyz_y = cbuffer.data(ip_geom_20_off + 22 * ccomps * dcomps);

            auto g_xx_0_xxxyyz_z = cbuffer.data(ip_geom_20_off + 23 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyyz_x, g_xx_0_xxxyyz_y, g_xx_0_xxxyyz_z, g_xx_0_xxxyz_x, g_xx_0_xxxyz_xy, g_xx_0_xxxyz_y, g_xx_0_xxxyz_yy, g_xx_0_xxxyz_yz, g_xx_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyyz_x[k] = -g_xx_0_xxxyz_x[k] * ab_y + g_xx_0_xxxyz_xy[k];

                g_xx_0_xxxyyz_y[k] = -g_xx_0_xxxyz_y[k] * ab_y + g_xx_0_xxxyz_yy[k];

                g_xx_0_xxxyyz_z[k] = -g_xx_0_xxxyz_z[k] * ab_y + g_xx_0_xxxyz_yz[k];
            }

            /// Set up 24-27 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxyzz_x = cbuffer.data(ip_geom_20_off + 24 * ccomps * dcomps);

            auto g_xx_0_xxxyzz_y = cbuffer.data(ip_geom_20_off + 25 * ccomps * dcomps);

            auto g_xx_0_xxxyzz_z = cbuffer.data(ip_geom_20_off + 26 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxyzz_x, g_xx_0_xxxyzz_y, g_xx_0_xxxyzz_z, g_xx_0_xxxzz_x, g_xx_0_xxxzz_xy, g_xx_0_xxxzz_y, g_xx_0_xxxzz_yy, g_xx_0_xxxzz_yz, g_xx_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxyzz_x[k] = -g_xx_0_xxxzz_x[k] * ab_y + g_xx_0_xxxzz_xy[k];

                g_xx_0_xxxyzz_y[k] = -g_xx_0_xxxzz_y[k] * ab_y + g_xx_0_xxxzz_yy[k];

                g_xx_0_xxxyzz_z[k] = -g_xx_0_xxxzz_z[k] * ab_y + g_xx_0_xxxzz_yz[k];
            }

            /// Set up 27-30 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxxzzz_x = cbuffer.data(ip_geom_20_off + 27 * ccomps * dcomps);

            auto g_xx_0_xxxzzz_y = cbuffer.data(ip_geom_20_off + 28 * ccomps * dcomps);

            auto g_xx_0_xxxzzz_z = cbuffer.data(ip_geom_20_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxxzz_x, g_xx_0_xxxzz_xz, g_xx_0_xxxzz_y, g_xx_0_xxxzz_yz, g_xx_0_xxxzz_z, g_xx_0_xxxzz_zz, g_xx_0_xxxzzz_x, g_xx_0_xxxzzz_y, g_xx_0_xxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxxzzz_x[k] = -g_xx_0_xxxzz_x[k] * ab_z + g_xx_0_xxxzz_xz[k];

                g_xx_0_xxxzzz_y[k] = -g_xx_0_xxxzz_y[k] * ab_z + g_xx_0_xxxzz_yz[k];

                g_xx_0_xxxzzz_z[k] = -g_xx_0_xxxzz_z[k] * ab_z + g_xx_0_xxxzz_zz[k];
            }

            /// Set up 30-33 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyyy_x = cbuffer.data(ip_geom_20_off + 30 * ccomps * dcomps);

            auto g_xx_0_xxyyyy_y = cbuffer.data(ip_geom_20_off + 31 * ccomps * dcomps);

            auto g_xx_0_xxyyyy_z = cbuffer.data(ip_geom_20_off + 32 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyy_x, g_xx_0_xxyyy_xy, g_xx_0_xxyyy_y, g_xx_0_xxyyy_yy, g_xx_0_xxyyy_yz, g_xx_0_xxyyy_z, g_xx_0_xxyyyy_x, g_xx_0_xxyyyy_y, g_xx_0_xxyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyyy_x[k] = -g_xx_0_xxyyy_x[k] * ab_y + g_xx_0_xxyyy_xy[k];

                g_xx_0_xxyyyy_y[k] = -g_xx_0_xxyyy_y[k] * ab_y + g_xx_0_xxyyy_yy[k];

                g_xx_0_xxyyyy_z[k] = -g_xx_0_xxyyy_z[k] * ab_y + g_xx_0_xxyyy_yz[k];
            }

            /// Set up 33-36 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyyz_x = cbuffer.data(ip_geom_20_off + 33 * ccomps * dcomps);

            auto g_xx_0_xxyyyz_y = cbuffer.data(ip_geom_20_off + 34 * ccomps * dcomps);

            auto g_xx_0_xxyyyz_z = cbuffer.data(ip_geom_20_off + 35 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyyz_x, g_xx_0_xxyyyz_y, g_xx_0_xxyyyz_z, g_xx_0_xxyyz_x, g_xx_0_xxyyz_xy, g_xx_0_xxyyz_y, g_xx_0_xxyyz_yy, g_xx_0_xxyyz_yz, g_xx_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyyz_x[k] = -g_xx_0_xxyyz_x[k] * ab_y + g_xx_0_xxyyz_xy[k];

                g_xx_0_xxyyyz_y[k] = -g_xx_0_xxyyz_y[k] * ab_y + g_xx_0_xxyyz_yy[k];

                g_xx_0_xxyyyz_z[k] = -g_xx_0_xxyyz_z[k] * ab_y + g_xx_0_xxyyz_yz[k];
            }

            /// Set up 36-39 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyyzz_x = cbuffer.data(ip_geom_20_off + 36 * ccomps * dcomps);

            auto g_xx_0_xxyyzz_y = cbuffer.data(ip_geom_20_off + 37 * ccomps * dcomps);

            auto g_xx_0_xxyyzz_z = cbuffer.data(ip_geom_20_off + 38 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyyzz_x, g_xx_0_xxyyzz_y, g_xx_0_xxyyzz_z, g_xx_0_xxyzz_x, g_xx_0_xxyzz_xy, g_xx_0_xxyzz_y, g_xx_0_xxyzz_yy, g_xx_0_xxyzz_yz, g_xx_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyyzz_x[k] = -g_xx_0_xxyzz_x[k] * ab_y + g_xx_0_xxyzz_xy[k];

                g_xx_0_xxyyzz_y[k] = -g_xx_0_xxyzz_y[k] * ab_y + g_xx_0_xxyzz_yy[k];

                g_xx_0_xxyyzz_z[k] = -g_xx_0_xxyzz_z[k] * ab_y + g_xx_0_xxyzz_yz[k];
            }

            /// Set up 39-42 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxyzzz_x = cbuffer.data(ip_geom_20_off + 39 * ccomps * dcomps);

            auto g_xx_0_xxyzzz_y = cbuffer.data(ip_geom_20_off + 40 * ccomps * dcomps);

            auto g_xx_0_xxyzzz_z = cbuffer.data(ip_geom_20_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxyzzz_x, g_xx_0_xxyzzz_y, g_xx_0_xxyzzz_z, g_xx_0_xxzzz_x, g_xx_0_xxzzz_xy, g_xx_0_xxzzz_y, g_xx_0_xxzzz_yy, g_xx_0_xxzzz_yz, g_xx_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxyzzz_x[k] = -g_xx_0_xxzzz_x[k] * ab_y + g_xx_0_xxzzz_xy[k];

                g_xx_0_xxyzzz_y[k] = -g_xx_0_xxzzz_y[k] * ab_y + g_xx_0_xxzzz_yy[k];

                g_xx_0_xxyzzz_z[k] = -g_xx_0_xxzzz_z[k] * ab_y + g_xx_0_xxzzz_yz[k];
            }

            /// Set up 42-45 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xxzzzz_x = cbuffer.data(ip_geom_20_off + 42 * ccomps * dcomps);

            auto g_xx_0_xxzzzz_y = cbuffer.data(ip_geom_20_off + 43 * ccomps * dcomps);

            auto g_xx_0_xxzzzz_z = cbuffer.data(ip_geom_20_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xxzzz_x, g_xx_0_xxzzz_xz, g_xx_0_xxzzz_y, g_xx_0_xxzzz_yz, g_xx_0_xxzzz_z, g_xx_0_xxzzz_zz, g_xx_0_xxzzzz_x, g_xx_0_xxzzzz_y, g_xx_0_xxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xxzzzz_x[k] = -g_xx_0_xxzzz_x[k] * ab_z + g_xx_0_xxzzz_xz[k];

                g_xx_0_xxzzzz_y[k] = -g_xx_0_xxzzz_y[k] * ab_z + g_xx_0_xxzzz_yz[k];

                g_xx_0_xxzzzz_z[k] = -g_xx_0_xxzzz_z[k] * ab_z + g_xx_0_xxzzz_zz[k];
            }

            /// Set up 45-48 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyyy_x = cbuffer.data(ip_geom_20_off + 45 * ccomps * dcomps);

            auto g_xx_0_xyyyyy_y = cbuffer.data(ip_geom_20_off + 46 * ccomps * dcomps);

            auto g_xx_0_xyyyyy_z = cbuffer.data(ip_geom_20_off + 47 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyy_x, g_xx_0_xyyyy_xy, g_xx_0_xyyyy_y, g_xx_0_xyyyy_yy, g_xx_0_xyyyy_yz, g_xx_0_xyyyy_z, g_xx_0_xyyyyy_x, g_xx_0_xyyyyy_y, g_xx_0_xyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyyy_x[k] = -g_xx_0_xyyyy_x[k] * ab_y + g_xx_0_xyyyy_xy[k];

                g_xx_0_xyyyyy_y[k] = -g_xx_0_xyyyy_y[k] * ab_y + g_xx_0_xyyyy_yy[k];

                g_xx_0_xyyyyy_z[k] = -g_xx_0_xyyyy_z[k] * ab_y + g_xx_0_xyyyy_yz[k];
            }

            /// Set up 48-51 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyyz_x = cbuffer.data(ip_geom_20_off + 48 * ccomps * dcomps);

            auto g_xx_0_xyyyyz_y = cbuffer.data(ip_geom_20_off + 49 * ccomps * dcomps);

            auto g_xx_0_xyyyyz_z = cbuffer.data(ip_geom_20_off + 50 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyyz_x, g_xx_0_xyyyyz_y, g_xx_0_xyyyyz_z, g_xx_0_xyyyz_x, g_xx_0_xyyyz_xy, g_xx_0_xyyyz_y, g_xx_0_xyyyz_yy, g_xx_0_xyyyz_yz, g_xx_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyyz_x[k] = -g_xx_0_xyyyz_x[k] * ab_y + g_xx_0_xyyyz_xy[k];

                g_xx_0_xyyyyz_y[k] = -g_xx_0_xyyyz_y[k] * ab_y + g_xx_0_xyyyz_yy[k];

                g_xx_0_xyyyyz_z[k] = -g_xx_0_xyyyz_z[k] * ab_y + g_xx_0_xyyyz_yz[k];
            }

            /// Set up 51-54 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyyzz_x = cbuffer.data(ip_geom_20_off + 51 * ccomps * dcomps);

            auto g_xx_0_xyyyzz_y = cbuffer.data(ip_geom_20_off + 52 * ccomps * dcomps);

            auto g_xx_0_xyyyzz_z = cbuffer.data(ip_geom_20_off + 53 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyyzz_x, g_xx_0_xyyyzz_y, g_xx_0_xyyyzz_z, g_xx_0_xyyzz_x, g_xx_0_xyyzz_xy, g_xx_0_xyyzz_y, g_xx_0_xyyzz_yy, g_xx_0_xyyzz_yz, g_xx_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyyzz_x[k] = -g_xx_0_xyyzz_x[k] * ab_y + g_xx_0_xyyzz_xy[k];

                g_xx_0_xyyyzz_y[k] = -g_xx_0_xyyzz_y[k] * ab_y + g_xx_0_xyyzz_yy[k];

                g_xx_0_xyyyzz_z[k] = -g_xx_0_xyyzz_z[k] * ab_y + g_xx_0_xyyzz_yz[k];
            }

            /// Set up 54-57 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyyzzz_x = cbuffer.data(ip_geom_20_off + 54 * ccomps * dcomps);

            auto g_xx_0_xyyzzz_y = cbuffer.data(ip_geom_20_off + 55 * ccomps * dcomps);

            auto g_xx_0_xyyzzz_z = cbuffer.data(ip_geom_20_off + 56 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyyzzz_x, g_xx_0_xyyzzz_y, g_xx_0_xyyzzz_z, g_xx_0_xyzzz_x, g_xx_0_xyzzz_xy, g_xx_0_xyzzz_y, g_xx_0_xyzzz_yy, g_xx_0_xyzzz_yz, g_xx_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyyzzz_x[k] = -g_xx_0_xyzzz_x[k] * ab_y + g_xx_0_xyzzz_xy[k];

                g_xx_0_xyyzzz_y[k] = -g_xx_0_xyzzz_y[k] * ab_y + g_xx_0_xyzzz_yy[k];

                g_xx_0_xyyzzz_z[k] = -g_xx_0_xyzzz_z[k] * ab_y + g_xx_0_xyzzz_yz[k];
            }

            /// Set up 57-60 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xyzzzz_x = cbuffer.data(ip_geom_20_off + 57 * ccomps * dcomps);

            auto g_xx_0_xyzzzz_y = cbuffer.data(ip_geom_20_off + 58 * ccomps * dcomps);

            auto g_xx_0_xyzzzz_z = cbuffer.data(ip_geom_20_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xyzzzz_x, g_xx_0_xyzzzz_y, g_xx_0_xyzzzz_z, g_xx_0_xzzzz_x, g_xx_0_xzzzz_xy, g_xx_0_xzzzz_y, g_xx_0_xzzzz_yy, g_xx_0_xzzzz_yz, g_xx_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xyzzzz_x[k] = -g_xx_0_xzzzz_x[k] * ab_y + g_xx_0_xzzzz_xy[k];

                g_xx_0_xyzzzz_y[k] = -g_xx_0_xzzzz_y[k] * ab_y + g_xx_0_xzzzz_yy[k];

                g_xx_0_xyzzzz_z[k] = -g_xx_0_xzzzz_z[k] * ab_y + g_xx_0_xzzzz_yz[k];
            }

            /// Set up 60-63 components of targeted buffer : cbuffer.data(

            auto g_xx_0_xzzzzz_x = cbuffer.data(ip_geom_20_off + 60 * ccomps * dcomps);

            auto g_xx_0_xzzzzz_y = cbuffer.data(ip_geom_20_off + 61 * ccomps * dcomps);

            auto g_xx_0_xzzzzz_z = cbuffer.data(ip_geom_20_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_xzzzz_x, g_xx_0_xzzzz_xz, g_xx_0_xzzzz_y, g_xx_0_xzzzz_yz, g_xx_0_xzzzz_z, g_xx_0_xzzzz_zz, g_xx_0_xzzzzz_x, g_xx_0_xzzzzz_y, g_xx_0_xzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_xzzzzz_x[k] = -g_xx_0_xzzzz_x[k] * ab_z + g_xx_0_xzzzz_xz[k];

                g_xx_0_xzzzzz_y[k] = -g_xx_0_xzzzz_y[k] * ab_z + g_xx_0_xzzzz_yz[k];

                g_xx_0_xzzzzz_z[k] = -g_xx_0_xzzzz_z[k] * ab_z + g_xx_0_xzzzz_zz[k];
            }

            /// Set up 63-66 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyyy_x = cbuffer.data(ip_geom_20_off + 63 * ccomps * dcomps);

            auto g_xx_0_yyyyyy_y = cbuffer.data(ip_geom_20_off + 64 * ccomps * dcomps);

            auto g_xx_0_yyyyyy_z = cbuffer.data(ip_geom_20_off + 65 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyy_x, g_xx_0_yyyyy_xy, g_xx_0_yyyyy_y, g_xx_0_yyyyy_yy, g_xx_0_yyyyy_yz, g_xx_0_yyyyy_z, g_xx_0_yyyyyy_x, g_xx_0_yyyyyy_y, g_xx_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyyy_x[k] = -g_xx_0_yyyyy_x[k] * ab_y + g_xx_0_yyyyy_xy[k];

                g_xx_0_yyyyyy_y[k] = -g_xx_0_yyyyy_y[k] * ab_y + g_xx_0_yyyyy_yy[k];

                g_xx_0_yyyyyy_z[k] = -g_xx_0_yyyyy_z[k] * ab_y + g_xx_0_yyyyy_yz[k];
            }

            /// Set up 66-69 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyyz_x = cbuffer.data(ip_geom_20_off + 66 * ccomps * dcomps);

            auto g_xx_0_yyyyyz_y = cbuffer.data(ip_geom_20_off + 67 * ccomps * dcomps);

            auto g_xx_0_yyyyyz_z = cbuffer.data(ip_geom_20_off + 68 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyyz_x, g_xx_0_yyyyyz_y, g_xx_0_yyyyyz_z, g_xx_0_yyyyz_x, g_xx_0_yyyyz_xy, g_xx_0_yyyyz_y, g_xx_0_yyyyz_yy, g_xx_0_yyyyz_yz, g_xx_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyyz_x[k] = -g_xx_0_yyyyz_x[k] * ab_y + g_xx_0_yyyyz_xy[k];

                g_xx_0_yyyyyz_y[k] = -g_xx_0_yyyyz_y[k] * ab_y + g_xx_0_yyyyz_yy[k];

                g_xx_0_yyyyyz_z[k] = -g_xx_0_yyyyz_z[k] * ab_y + g_xx_0_yyyyz_yz[k];
            }

            /// Set up 69-72 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyyzz_x = cbuffer.data(ip_geom_20_off + 69 * ccomps * dcomps);

            auto g_xx_0_yyyyzz_y = cbuffer.data(ip_geom_20_off + 70 * ccomps * dcomps);

            auto g_xx_0_yyyyzz_z = cbuffer.data(ip_geom_20_off + 71 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyyzz_x, g_xx_0_yyyyzz_y, g_xx_0_yyyyzz_z, g_xx_0_yyyzz_x, g_xx_0_yyyzz_xy, g_xx_0_yyyzz_y, g_xx_0_yyyzz_yy, g_xx_0_yyyzz_yz, g_xx_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyyzz_x[k] = -g_xx_0_yyyzz_x[k] * ab_y + g_xx_0_yyyzz_xy[k];

                g_xx_0_yyyyzz_y[k] = -g_xx_0_yyyzz_y[k] * ab_y + g_xx_0_yyyzz_yy[k];

                g_xx_0_yyyyzz_z[k] = -g_xx_0_yyyzz_z[k] * ab_y + g_xx_0_yyyzz_yz[k];
            }

            /// Set up 72-75 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyyzzz_x = cbuffer.data(ip_geom_20_off + 72 * ccomps * dcomps);

            auto g_xx_0_yyyzzz_y = cbuffer.data(ip_geom_20_off + 73 * ccomps * dcomps);

            auto g_xx_0_yyyzzz_z = cbuffer.data(ip_geom_20_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyyzzz_x, g_xx_0_yyyzzz_y, g_xx_0_yyyzzz_z, g_xx_0_yyzzz_x, g_xx_0_yyzzz_xy, g_xx_0_yyzzz_y, g_xx_0_yyzzz_yy, g_xx_0_yyzzz_yz, g_xx_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyyzzz_x[k] = -g_xx_0_yyzzz_x[k] * ab_y + g_xx_0_yyzzz_xy[k];

                g_xx_0_yyyzzz_y[k] = -g_xx_0_yyzzz_y[k] * ab_y + g_xx_0_yyzzz_yy[k];

                g_xx_0_yyyzzz_z[k] = -g_xx_0_yyzzz_z[k] * ab_y + g_xx_0_yyzzz_yz[k];
            }

            /// Set up 75-78 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yyzzzz_x = cbuffer.data(ip_geom_20_off + 75 * ccomps * dcomps);

            auto g_xx_0_yyzzzz_y = cbuffer.data(ip_geom_20_off + 76 * ccomps * dcomps);

            auto g_xx_0_yyzzzz_z = cbuffer.data(ip_geom_20_off + 77 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yyzzzz_x, g_xx_0_yyzzzz_y, g_xx_0_yyzzzz_z, g_xx_0_yzzzz_x, g_xx_0_yzzzz_xy, g_xx_0_yzzzz_y, g_xx_0_yzzzz_yy, g_xx_0_yzzzz_yz, g_xx_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yyzzzz_x[k] = -g_xx_0_yzzzz_x[k] * ab_y + g_xx_0_yzzzz_xy[k];

                g_xx_0_yyzzzz_y[k] = -g_xx_0_yzzzz_y[k] * ab_y + g_xx_0_yzzzz_yy[k];

                g_xx_0_yyzzzz_z[k] = -g_xx_0_yzzzz_z[k] * ab_y + g_xx_0_yzzzz_yz[k];
            }

            /// Set up 78-81 components of targeted buffer : cbuffer.data(

            auto g_xx_0_yzzzzz_x = cbuffer.data(ip_geom_20_off + 78 * ccomps * dcomps);

            auto g_xx_0_yzzzzz_y = cbuffer.data(ip_geom_20_off + 79 * ccomps * dcomps);

            auto g_xx_0_yzzzzz_z = cbuffer.data(ip_geom_20_off + 80 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_yzzzzz_x, g_xx_0_yzzzzz_y, g_xx_0_yzzzzz_z, g_xx_0_zzzzz_x, g_xx_0_zzzzz_xy, g_xx_0_zzzzz_y, g_xx_0_zzzzz_yy, g_xx_0_zzzzz_yz, g_xx_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_yzzzzz_x[k] = -g_xx_0_zzzzz_x[k] * ab_y + g_xx_0_zzzzz_xy[k];

                g_xx_0_yzzzzz_y[k] = -g_xx_0_zzzzz_y[k] * ab_y + g_xx_0_zzzzz_yy[k];

                g_xx_0_yzzzzz_z[k] = -g_xx_0_zzzzz_z[k] * ab_y + g_xx_0_zzzzz_yz[k];
            }

            /// Set up 81-84 components of targeted buffer : cbuffer.data(

            auto g_xx_0_zzzzzz_x = cbuffer.data(ip_geom_20_off + 81 * ccomps * dcomps);

            auto g_xx_0_zzzzzz_y = cbuffer.data(ip_geom_20_off + 82 * ccomps * dcomps);

            auto g_xx_0_zzzzzz_z = cbuffer.data(ip_geom_20_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_xx_0_zzzzz_x, g_xx_0_zzzzz_xz, g_xx_0_zzzzz_y, g_xx_0_zzzzz_yz, g_xx_0_zzzzz_z, g_xx_0_zzzzz_zz, g_xx_0_zzzzzz_x, g_xx_0_zzzzzz_y, g_xx_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xx_0_zzzzzz_x[k] = -g_xx_0_zzzzz_x[k] * ab_z + g_xx_0_zzzzz_xz[k];

                g_xx_0_zzzzzz_y[k] = -g_xx_0_zzzzz_y[k] * ab_z + g_xx_0_zzzzz_yz[k];

                g_xx_0_zzzzzz_z[k] = -g_xx_0_zzzzz_z[k] * ab_z + g_xx_0_zzzzz_zz[k];
            }

            /// Set up 84-87 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxxx_x = cbuffer.data(ip_geom_20_off + 84 * ccomps * dcomps);

            auto g_xy_0_xxxxxx_y = cbuffer.data(ip_geom_20_off + 85 * ccomps * dcomps);

            auto g_xy_0_xxxxxx_z = cbuffer.data(ip_geom_20_off + 86 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxx_x, g_xy_0_xxxxx_xx, g_xy_0_xxxxx_xy, g_xy_0_xxxxx_xz, g_xy_0_xxxxx_y, g_xy_0_xxxxx_z, g_xy_0_xxxxxx_x, g_xy_0_xxxxxx_y, g_xy_0_xxxxxx_z, g_y_0_xxxxx_x, g_y_0_xxxxx_y, g_y_0_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxxx_x[k] = -g_y_0_xxxxx_x[k] - g_xy_0_xxxxx_x[k] * ab_x + g_xy_0_xxxxx_xx[k];

                g_xy_0_xxxxxx_y[k] = -g_y_0_xxxxx_y[k] - g_xy_0_xxxxx_y[k] * ab_x + g_xy_0_xxxxx_xy[k];

                g_xy_0_xxxxxx_z[k] = -g_y_0_xxxxx_z[k] - g_xy_0_xxxxx_z[k] * ab_x + g_xy_0_xxxxx_xz[k];
            }

            /// Set up 87-90 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxxy_x = cbuffer.data(ip_geom_20_off + 87 * ccomps * dcomps);

            auto g_xy_0_xxxxxy_y = cbuffer.data(ip_geom_20_off + 88 * ccomps * dcomps);

            auto g_xy_0_xxxxxy_z = cbuffer.data(ip_geom_20_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxxy_x, g_xy_0_xxxxxy_y, g_xy_0_xxxxxy_z, g_xy_0_xxxxy_x, g_xy_0_xxxxy_xx, g_xy_0_xxxxy_xy, g_xy_0_xxxxy_xz, g_xy_0_xxxxy_y, g_xy_0_xxxxy_z, g_y_0_xxxxy_x, g_y_0_xxxxy_y, g_y_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxxy_x[k] = -g_y_0_xxxxy_x[k] - g_xy_0_xxxxy_x[k] * ab_x + g_xy_0_xxxxy_xx[k];

                g_xy_0_xxxxxy_y[k] = -g_y_0_xxxxy_y[k] - g_xy_0_xxxxy_y[k] * ab_x + g_xy_0_xxxxy_xy[k];

                g_xy_0_xxxxxy_z[k] = -g_y_0_xxxxy_z[k] - g_xy_0_xxxxy_z[k] * ab_x + g_xy_0_xxxxy_xz[k];
            }

            /// Set up 90-93 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxxz_x = cbuffer.data(ip_geom_20_off + 90 * ccomps * dcomps);

            auto g_xy_0_xxxxxz_y = cbuffer.data(ip_geom_20_off + 91 * ccomps * dcomps);

            auto g_xy_0_xxxxxz_z = cbuffer.data(ip_geom_20_off + 92 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxx_x, g_xy_0_xxxxx_xz, g_xy_0_xxxxx_y, g_xy_0_xxxxx_yz, g_xy_0_xxxxx_z, g_xy_0_xxxxx_zz, g_xy_0_xxxxxz_x, g_xy_0_xxxxxz_y, g_xy_0_xxxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxxz_x[k] = -g_xy_0_xxxxx_x[k] * ab_z + g_xy_0_xxxxx_xz[k];

                g_xy_0_xxxxxz_y[k] = -g_xy_0_xxxxx_y[k] * ab_z + g_xy_0_xxxxx_yz[k];

                g_xy_0_xxxxxz_z[k] = -g_xy_0_xxxxx_z[k] * ab_z + g_xy_0_xxxxx_zz[k];
            }

            /// Set up 93-96 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxyy_x = cbuffer.data(ip_geom_20_off + 93 * ccomps * dcomps);

            auto g_xy_0_xxxxyy_y = cbuffer.data(ip_geom_20_off + 94 * ccomps * dcomps);

            auto g_xy_0_xxxxyy_z = cbuffer.data(ip_geom_20_off + 95 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxyy_x, g_xy_0_xxxxyy_y, g_xy_0_xxxxyy_z, g_xy_0_xxxyy_x, g_xy_0_xxxyy_xx, g_xy_0_xxxyy_xy, g_xy_0_xxxyy_xz, g_xy_0_xxxyy_y, g_xy_0_xxxyy_z, g_y_0_xxxyy_x, g_y_0_xxxyy_y, g_y_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxyy_x[k] = -g_y_0_xxxyy_x[k] - g_xy_0_xxxyy_x[k] * ab_x + g_xy_0_xxxyy_xx[k];

                g_xy_0_xxxxyy_y[k] = -g_y_0_xxxyy_y[k] - g_xy_0_xxxyy_y[k] * ab_x + g_xy_0_xxxyy_xy[k];

                g_xy_0_xxxxyy_z[k] = -g_y_0_xxxyy_z[k] - g_xy_0_xxxyy_z[k] * ab_x + g_xy_0_xxxyy_xz[k];
            }

            /// Set up 96-99 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxyz_x = cbuffer.data(ip_geom_20_off + 96 * ccomps * dcomps);

            auto g_xy_0_xxxxyz_y = cbuffer.data(ip_geom_20_off + 97 * ccomps * dcomps);

            auto g_xy_0_xxxxyz_z = cbuffer.data(ip_geom_20_off + 98 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxy_x, g_xy_0_xxxxy_xz, g_xy_0_xxxxy_y, g_xy_0_xxxxy_yz, g_xy_0_xxxxy_z, g_xy_0_xxxxy_zz, g_xy_0_xxxxyz_x, g_xy_0_xxxxyz_y, g_xy_0_xxxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxyz_x[k] = -g_xy_0_xxxxy_x[k] * ab_z + g_xy_0_xxxxy_xz[k];

                g_xy_0_xxxxyz_y[k] = -g_xy_0_xxxxy_y[k] * ab_z + g_xy_0_xxxxy_yz[k];

                g_xy_0_xxxxyz_z[k] = -g_xy_0_xxxxy_z[k] * ab_z + g_xy_0_xxxxy_zz[k];
            }

            /// Set up 99-102 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxxzz_x = cbuffer.data(ip_geom_20_off + 99 * ccomps * dcomps);

            auto g_xy_0_xxxxzz_y = cbuffer.data(ip_geom_20_off + 100 * ccomps * dcomps);

            auto g_xy_0_xxxxzz_z = cbuffer.data(ip_geom_20_off + 101 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxxz_x, g_xy_0_xxxxz_xz, g_xy_0_xxxxz_y, g_xy_0_xxxxz_yz, g_xy_0_xxxxz_z, g_xy_0_xxxxz_zz, g_xy_0_xxxxzz_x, g_xy_0_xxxxzz_y, g_xy_0_xxxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxxzz_x[k] = -g_xy_0_xxxxz_x[k] * ab_z + g_xy_0_xxxxz_xz[k];

                g_xy_0_xxxxzz_y[k] = -g_xy_0_xxxxz_y[k] * ab_z + g_xy_0_xxxxz_yz[k];

                g_xy_0_xxxxzz_z[k] = -g_xy_0_xxxxz_z[k] * ab_z + g_xy_0_xxxxz_zz[k];
            }

            /// Set up 102-105 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyyy_x = cbuffer.data(ip_geom_20_off + 102 * ccomps * dcomps);

            auto g_xy_0_xxxyyy_y = cbuffer.data(ip_geom_20_off + 103 * ccomps * dcomps);

            auto g_xy_0_xxxyyy_z = cbuffer.data(ip_geom_20_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyyy_x, g_xy_0_xxxyyy_y, g_xy_0_xxxyyy_z, g_xy_0_xxyyy_x, g_xy_0_xxyyy_xx, g_xy_0_xxyyy_xy, g_xy_0_xxyyy_xz, g_xy_0_xxyyy_y, g_xy_0_xxyyy_z, g_y_0_xxyyy_x, g_y_0_xxyyy_y, g_y_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyyy_x[k] = -g_y_0_xxyyy_x[k] - g_xy_0_xxyyy_x[k] * ab_x + g_xy_0_xxyyy_xx[k];

                g_xy_0_xxxyyy_y[k] = -g_y_0_xxyyy_y[k] - g_xy_0_xxyyy_y[k] * ab_x + g_xy_0_xxyyy_xy[k];

                g_xy_0_xxxyyy_z[k] = -g_y_0_xxyyy_z[k] - g_xy_0_xxyyy_z[k] * ab_x + g_xy_0_xxyyy_xz[k];
            }

            /// Set up 105-108 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyyz_x = cbuffer.data(ip_geom_20_off + 105 * ccomps * dcomps);

            auto g_xy_0_xxxyyz_y = cbuffer.data(ip_geom_20_off + 106 * ccomps * dcomps);

            auto g_xy_0_xxxyyz_z = cbuffer.data(ip_geom_20_off + 107 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyy_x, g_xy_0_xxxyy_xz, g_xy_0_xxxyy_y, g_xy_0_xxxyy_yz, g_xy_0_xxxyy_z, g_xy_0_xxxyy_zz, g_xy_0_xxxyyz_x, g_xy_0_xxxyyz_y, g_xy_0_xxxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyyz_x[k] = -g_xy_0_xxxyy_x[k] * ab_z + g_xy_0_xxxyy_xz[k];

                g_xy_0_xxxyyz_y[k] = -g_xy_0_xxxyy_y[k] * ab_z + g_xy_0_xxxyy_yz[k];

                g_xy_0_xxxyyz_z[k] = -g_xy_0_xxxyy_z[k] * ab_z + g_xy_0_xxxyy_zz[k];
            }

            /// Set up 108-111 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxyzz_x = cbuffer.data(ip_geom_20_off + 108 * ccomps * dcomps);

            auto g_xy_0_xxxyzz_y = cbuffer.data(ip_geom_20_off + 109 * ccomps * dcomps);

            auto g_xy_0_xxxyzz_z = cbuffer.data(ip_geom_20_off + 110 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxyz_x, g_xy_0_xxxyz_xz, g_xy_0_xxxyz_y, g_xy_0_xxxyz_yz, g_xy_0_xxxyz_z, g_xy_0_xxxyz_zz, g_xy_0_xxxyzz_x, g_xy_0_xxxyzz_y, g_xy_0_xxxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxyzz_x[k] = -g_xy_0_xxxyz_x[k] * ab_z + g_xy_0_xxxyz_xz[k];

                g_xy_0_xxxyzz_y[k] = -g_xy_0_xxxyz_y[k] * ab_z + g_xy_0_xxxyz_yz[k];

                g_xy_0_xxxyzz_z[k] = -g_xy_0_xxxyz_z[k] * ab_z + g_xy_0_xxxyz_zz[k];
            }

            /// Set up 111-114 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxxzzz_x = cbuffer.data(ip_geom_20_off + 111 * ccomps * dcomps);

            auto g_xy_0_xxxzzz_y = cbuffer.data(ip_geom_20_off + 112 * ccomps * dcomps);

            auto g_xy_0_xxxzzz_z = cbuffer.data(ip_geom_20_off + 113 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxxzz_x, g_xy_0_xxxzz_xz, g_xy_0_xxxzz_y, g_xy_0_xxxzz_yz, g_xy_0_xxxzz_z, g_xy_0_xxxzz_zz, g_xy_0_xxxzzz_x, g_xy_0_xxxzzz_y, g_xy_0_xxxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxxzzz_x[k] = -g_xy_0_xxxzz_x[k] * ab_z + g_xy_0_xxxzz_xz[k];

                g_xy_0_xxxzzz_y[k] = -g_xy_0_xxxzz_y[k] * ab_z + g_xy_0_xxxzz_yz[k];

                g_xy_0_xxxzzz_z[k] = -g_xy_0_xxxzz_z[k] * ab_z + g_xy_0_xxxzz_zz[k];
            }

            /// Set up 114-117 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyyy_x = cbuffer.data(ip_geom_20_off + 114 * ccomps * dcomps);

            auto g_xy_0_xxyyyy_y = cbuffer.data(ip_geom_20_off + 115 * ccomps * dcomps);

            auto g_xy_0_xxyyyy_z = cbuffer.data(ip_geom_20_off + 116 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyyy_x, g_xy_0_xxyyyy_y, g_xy_0_xxyyyy_z, g_xy_0_xyyyy_x, g_xy_0_xyyyy_xx, g_xy_0_xyyyy_xy, g_xy_0_xyyyy_xz, g_xy_0_xyyyy_y, g_xy_0_xyyyy_z, g_y_0_xyyyy_x, g_y_0_xyyyy_y, g_y_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyyy_x[k] = -g_y_0_xyyyy_x[k] - g_xy_0_xyyyy_x[k] * ab_x + g_xy_0_xyyyy_xx[k];

                g_xy_0_xxyyyy_y[k] = -g_y_0_xyyyy_y[k] - g_xy_0_xyyyy_y[k] * ab_x + g_xy_0_xyyyy_xy[k];

                g_xy_0_xxyyyy_z[k] = -g_y_0_xyyyy_z[k] - g_xy_0_xyyyy_z[k] * ab_x + g_xy_0_xyyyy_xz[k];
            }

            /// Set up 117-120 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyyz_x = cbuffer.data(ip_geom_20_off + 117 * ccomps * dcomps);

            auto g_xy_0_xxyyyz_y = cbuffer.data(ip_geom_20_off + 118 * ccomps * dcomps);

            auto g_xy_0_xxyyyz_z = cbuffer.data(ip_geom_20_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyy_x, g_xy_0_xxyyy_xz, g_xy_0_xxyyy_y, g_xy_0_xxyyy_yz, g_xy_0_xxyyy_z, g_xy_0_xxyyy_zz, g_xy_0_xxyyyz_x, g_xy_0_xxyyyz_y, g_xy_0_xxyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyyz_x[k] = -g_xy_0_xxyyy_x[k] * ab_z + g_xy_0_xxyyy_xz[k];

                g_xy_0_xxyyyz_y[k] = -g_xy_0_xxyyy_y[k] * ab_z + g_xy_0_xxyyy_yz[k];

                g_xy_0_xxyyyz_z[k] = -g_xy_0_xxyyy_z[k] * ab_z + g_xy_0_xxyyy_zz[k];
            }

            /// Set up 120-123 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyyzz_x = cbuffer.data(ip_geom_20_off + 120 * ccomps * dcomps);

            auto g_xy_0_xxyyzz_y = cbuffer.data(ip_geom_20_off + 121 * ccomps * dcomps);

            auto g_xy_0_xxyyzz_z = cbuffer.data(ip_geom_20_off + 122 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyyz_x, g_xy_0_xxyyz_xz, g_xy_0_xxyyz_y, g_xy_0_xxyyz_yz, g_xy_0_xxyyz_z, g_xy_0_xxyyz_zz, g_xy_0_xxyyzz_x, g_xy_0_xxyyzz_y, g_xy_0_xxyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyyzz_x[k] = -g_xy_0_xxyyz_x[k] * ab_z + g_xy_0_xxyyz_xz[k];

                g_xy_0_xxyyzz_y[k] = -g_xy_0_xxyyz_y[k] * ab_z + g_xy_0_xxyyz_yz[k];

                g_xy_0_xxyyzz_z[k] = -g_xy_0_xxyyz_z[k] * ab_z + g_xy_0_xxyyz_zz[k];
            }

            /// Set up 123-126 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxyzzz_x = cbuffer.data(ip_geom_20_off + 123 * ccomps * dcomps);

            auto g_xy_0_xxyzzz_y = cbuffer.data(ip_geom_20_off + 124 * ccomps * dcomps);

            auto g_xy_0_xxyzzz_z = cbuffer.data(ip_geom_20_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxyzz_x, g_xy_0_xxyzz_xz, g_xy_0_xxyzz_y, g_xy_0_xxyzz_yz, g_xy_0_xxyzz_z, g_xy_0_xxyzz_zz, g_xy_0_xxyzzz_x, g_xy_0_xxyzzz_y, g_xy_0_xxyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxyzzz_x[k] = -g_xy_0_xxyzz_x[k] * ab_z + g_xy_0_xxyzz_xz[k];

                g_xy_0_xxyzzz_y[k] = -g_xy_0_xxyzz_y[k] * ab_z + g_xy_0_xxyzz_yz[k];

                g_xy_0_xxyzzz_z[k] = -g_xy_0_xxyzz_z[k] * ab_z + g_xy_0_xxyzz_zz[k];
            }

            /// Set up 126-129 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xxzzzz_x = cbuffer.data(ip_geom_20_off + 126 * ccomps * dcomps);

            auto g_xy_0_xxzzzz_y = cbuffer.data(ip_geom_20_off + 127 * ccomps * dcomps);

            auto g_xy_0_xxzzzz_z = cbuffer.data(ip_geom_20_off + 128 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xxzzz_x, g_xy_0_xxzzz_xz, g_xy_0_xxzzz_y, g_xy_0_xxzzz_yz, g_xy_0_xxzzz_z, g_xy_0_xxzzz_zz, g_xy_0_xxzzzz_x, g_xy_0_xxzzzz_y, g_xy_0_xxzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xxzzzz_x[k] = -g_xy_0_xxzzz_x[k] * ab_z + g_xy_0_xxzzz_xz[k];

                g_xy_0_xxzzzz_y[k] = -g_xy_0_xxzzz_y[k] * ab_z + g_xy_0_xxzzz_yz[k];

                g_xy_0_xxzzzz_z[k] = -g_xy_0_xxzzz_z[k] * ab_z + g_xy_0_xxzzz_zz[k];
            }

            /// Set up 129-132 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyyy_x = cbuffer.data(ip_geom_20_off + 129 * ccomps * dcomps);

            auto g_xy_0_xyyyyy_y = cbuffer.data(ip_geom_20_off + 130 * ccomps * dcomps);

            auto g_xy_0_xyyyyy_z = cbuffer.data(ip_geom_20_off + 131 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyyy_x, g_xy_0_xyyyyy_y, g_xy_0_xyyyyy_z, g_xy_0_yyyyy_x, g_xy_0_yyyyy_xx, g_xy_0_yyyyy_xy, g_xy_0_yyyyy_xz, g_xy_0_yyyyy_y, g_xy_0_yyyyy_z, g_y_0_yyyyy_x, g_y_0_yyyyy_y, g_y_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyyy_x[k] = -g_y_0_yyyyy_x[k] - g_xy_0_yyyyy_x[k] * ab_x + g_xy_0_yyyyy_xx[k];

                g_xy_0_xyyyyy_y[k] = -g_y_0_yyyyy_y[k] - g_xy_0_yyyyy_y[k] * ab_x + g_xy_0_yyyyy_xy[k];

                g_xy_0_xyyyyy_z[k] = -g_y_0_yyyyy_z[k] - g_xy_0_yyyyy_z[k] * ab_x + g_xy_0_yyyyy_xz[k];
            }

            /// Set up 132-135 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyyz_x = cbuffer.data(ip_geom_20_off + 132 * ccomps * dcomps);

            auto g_xy_0_xyyyyz_y = cbuffer.data(ip_geom_20_off + 133 * ccomps * dcomps);

            auto g_xy_0_xyyyyz_z = cbuffer.data(ip_geom_20_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyy_x, g_xy_0_xyyyy_xz, g_xy_0_xyyyy_y, g_xy_0_xyyyy_yz, g_xy_0_xyyyy_z, g_xy_0_xyyyy_zz, g_xy_0_xyyyyz_x, g_xy_0_xyyyyz_y, g_xy_0_xyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyyz_x[k] = -g_xy_0_xyyyy_x[k] * ab_z + g_xy_0_xyyyy_xz[k];

                g_xy_0_xyyyyz_y[k] = -g_xy_0_xyyyy_y[k] * ab_z + g_xy_0_xyyyy_yz[k];

                g_xy_0_xyyyyz_z[k] = -g_xy_0_xyyyy_z[k] * ab_z + g_xy_0_xyyyy_zz[k];
            }

            /// Set up 135-138 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyyzz_x = cbuffer.data(ip_geom_20_off + 135 * ccomps * dcomps);

            auto g_xy_0_xyyyzz_y = cbuffer.data(ip_geom_20_off + 136 * ccomps * dcomps);

            auto g_xy_0_xyyyzz_z = cbuffer.data(ip_geom_20_off + 137 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyyz_x, g_xy_0_xyyyz_xz, g_xy_0_xyyyz_y, g_xy_0_xyyyz_yz, g_xy_0_xyyyz_z, g_xy_0_xyyyz_zz, g_xy_0_xyyyzz_x, g_xy_0_xyyyzz_y, g_xy_0_xyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyyzz_x[k] = -g_xy_0_xyyyz_x[k] * ab_z + g_xy_0_xyyyz_xz[k];

                g_xy_0_xyyyzz_y[k] = -g_xy_0_xyyyz_y[k] * ab_z + g_xy_0_xyyyz_yz[k];

                g_xy_0_xyyyzz_z[k] = -g_xy_0_xyyyz_z[k] * ab_z + g_xy_0_xyyyz_zz[k];
            }

            /// Set up 138-141 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyyzzz_x = cbuffer.data(ip_geom_20_off + 138 * ccomps * dcomps);

            auto g_xy_0_xyyzzz_y = cbuffer.data(ip_geom_20_off + 139 * ccomps * dcomps);

            auto g_xy_0_xyyzzz_z = cbuffer.data(ip_geom_20_off + 140 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyyzz_x, g_xy_0_xyyzz_xz, g_xy_0_xyyzz_y, g_xy_0_xyyzz_yz, g_xy_0_xyyzz_z, g_xy_0_xyyzz_zz, g_xy_0_xyyzzz_x, g_xy_0_xyyzzz_y, g_xy_0_xyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyyzzz_x[k] = -g_xy_0_xyyzz_x[k] * ab_z + g_xy_0_xyyzz_xz[k];

                g_xy_0_xyyzzz_y[k] = -g_xy_0_xyyzz_y[k] * ab_z + g_xy_0_xyyzz_yz[k];

                g_xy_0_xyyzzz_z[k] = -g_xy_0_xyyzz_z[k] * ab_z + g_xy_0_xyyzz_zz[k];
            }

            /// Set up 141-144 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xyzzzz_x = cbuffer.data(ip_geom_20_off + 141 * ccomps * dcomps);

            auto g_xy_0_xyzzzz_y = cbuffer.data(ip_geom_20_off + 142 * ccomps * dcomps);

            auto g_xy_0_xyzzzz_z = cbuffer.data(ip_geom_20_off + 143 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xyzzz_x, g_xy_0_xyzzz_xz, g_xy_0_xyzzz_y, g_xy_0_xyzzz_yz, g_xy_0_xyzzz_z, g_xy_0_xyzzz_zz, g_xy_0_xyzzzz_x, g_xy_0_xyzzzz_y, g_xy_0_xyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xyzzzz_x[k] = -g_xy_0_xyzzz_x[k] * ab_z + g_xy_0_xyzzz_xz[k];

                g_xy_0_xyzzzz_y[k] = -g_xy_0_xyzzz_y[k] * ab_z + g_xy_0_xyzzz_yz[k];

                g_xy_0_xyzzzz_z[k] = -g_xy_0_xyzzz_z[k] * ab_z + g_xy_0_xyzzz_zz[k];
            }

            /// Set up 144-147 components of targeted buffer : cbuffer.data(

            auto g_xy_0_xzzzzz_x = cbuffer.data(ip_geom_20_off + 144 * ccomps * dcomps);

            auto g_xy_0_xzzzzz_y = cbuffer.data(ip_geom_20_off + 145 * ccomps * dcomps);

            auto g_xy_0_xzzzzz_z = cbuffer.data(ip_geom_20_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_xzzzz_x, g_xy_0_xzzzz_xz, g_xy_0_xzzzz_y, g_xy_0_xzzzz_yz, g_xy_0_xzzzz_z, g_xy_0_xzzzz_zz, g_xy_0_xzzzzz_x, g_xy_0_xzzzzz_y, g_xy_0_xzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_xzzzzz_x[k] = -g_xy_0_xzzzz_x[k] * ab_z + g_xy_0_xzzzz_xz[k];

                g_xy_0_xzzzzz_y[k] = -g_xy_0_xzzzz_y[k] * ab_z + g_xy_0_xzzzz_yz[k];

                g_xy_0_xzzzzz_z[k] = -g_xy_0_xzzzz_z[k] * ab_z + g_xy_0_xzzzz_zz[k];
            }

            /// Set up 147-150 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyyy_x = cbuffer.data(ip_geom_20_off + 147 * ccomps * dcomps);

            auto g_xy_0_yyyyyy_y = cbuffer.data(ip_geom_20_off + 148 * ccomps * dcomps);

            auto g_xy_0_yyyyyy_z = cbuffer.data(ip_geom_20_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyy_x, g_x_0_yyyyy_y, g_x_0_yyyyy_z, g_xy_0_yyyyy_x, g_xy_0_yyyyy_xy, g_xy_0_yyyyy_y, g_xy_0_yyyyy_yy, g_xy_0_yyyyy_yz, g_xy_0_yyyyy_z, g_xy_0_yyyyyy_x, g_xy_0_yyyyyy_y, g_xy_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyyy_x[k] = -g_x_0_yyyyy_x[k] - g_xy_0_yyyyy_x[k] * ab_y + g_xy_0_yyyyy_xy[k];

                g_xy_0_yyyyyy_y[k] = -g_x_0_yyyyy_y[k] - g_xy_0_yyyyy_y[k] * ab_y + g_xy_0_yyyyy_yy[k];

                g_xy_0_yyyyyy_z[k] = -g_x_0_yyyyy_z[k] - g_xy_0_yyyyy_z[k] * ab_y + g_xy_0_yyyyy_yz[k];
            }

            /// Set up 150-153 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyyz_x = cbuffer.data(ip_geom_20_off + 150 * ccomps * dcomps);

            auto g_xy_0_yyyyyz_y = cbuffer.data(ip_geom_20_off + 151 * ccomps * dcomps);

            auto g_xy_0_yyyyyz_z = cbuffer.data(ip_geom_20_off + 152 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyyy_x, g_xy_0_yyyyy_xz, g_xy_0_yyyyy_y, g_xy_0_yyyyy_yz, g_xy_0_yyyyy_z, g_xy_0_yyyyy_zz, g_xy_0_yyyyyz_x, g_xy_0_yyyyyz_y, g_xy_0_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyyz_x[k] = -g_xy_0_yyyyy_x[k] * ab_z + g_xy_0_yyyyy_xz[k];

                g_xy_0_yyyyyz_y[k] = -g_xy_0_yyyyy_y[k] * ab_z + g_xy_0_yyyyy_yz[k];

                g_xy_0_yyyyyz_z[k] = -g_xy_0_yyyyy_z[k] * ab_z + g_xy_0_yyyyy_zz[k];
            }

            /// Set up 153-156 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyyzz_x = cbuffer.data(ip_geom_20_off + 153 * ccomps * dcomps);

            auto g_xy_0_yyyyzz_y = cbuffer.data(ip_geom_20_off + 154 * ccomps * dcomps);

            auto g_xy_0_yyyyzz_z = cbuffer.data(ip_geom_20_off + 155 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyyz_x, g_xy_0_yyyyz_xz, g_xy_0_yyyyz_y, g_xy_0_yyyyz_yz, g_xy_0_yyyyz_z, g_xy_0_yyyyz_zz, g_xy_0_yyyyzz_x, g_xy_0_yyyyzz_y, g_xy_0_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyyzz_x[k] = -g_xy_0_yyyyz_x[k] * ab_z + g_xy_0_yyyyz_xz[k];

                g_xy_0_yyyyzz_y[k] = -g_xy_0_yyyyz_y[k] * ab_z + g_xy_0_yyyyz_yz[k];

                g_xy_0_yyyyzz_z[k] = -g_xy_0_yyyyz_z[k] * ab_z + g_xy_0_yyyyz_zz[k];
            }

            /// Set up 156-159 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyyzzz_x = cbuffer.data(ip_geom_20_off + 156 * ccomps * dcomps);

            auto g_xy_0_yyyzzz_y = cbuffer.data(ip_geom_20_off + 157 * ccomps * dcomps);

            auto g_xy_0_yyyzzz_z = cbuffer.data(ip_geom_20_off + 158 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyyzz_x, g_xy_0_yyyzz_xz, g_xy_0_yyyzz_y, g_xy_0_yyyzz_yz, g_xy_0_yyyzz_z, g_xy_0_yyyzz_zz, g_xy_0_yyyzzz_x, g_xy_0_yyyzzz_y, g_xy_0_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyyzzz_x[k] = -g_xy_0_yyyzz_x[k] * ab_z + g_xy_0_yyyzz_xz[k];

                g_xy_0_yyyzzz_y[k] = -g_xy_0_yyyzz_y[k] * ab_z + g_xy_0_yyyzz_yz[k];

                g_xy_0_yyyzzz_z[k] = -g_xy_0_yyyzz_z[k] * ab_z + g_xy_0_yyyzz_zz[k];
            }

            /// Set up 159-162 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yyzzzz_x = cbuffer.data(ip_geom_20_off + 159 * ccomps * dcomps);

            auto g_xy_0_yyzzzz_y = cbuffer.data(ip_geom_20_off + 160 * ccomps * dcomps);

            auto g_xy_0_yyzzzz_z = cbuffer.data(ip_geom_20_off + 161 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yyzzz_x, g_xy_0_yyzzz_xz, g_xy_0_yyzzz_y, g_xy_0_yyzzz_yz, g_xy_0_yyzzz_z, g_xy_0_yyzzz_zz, g_xy_0_yyzzzz_x, g_xy_0_yyzzzz_y, g_xy_0_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yyzzzz_x[k] = -g_xy_0_yyzzz_x[k] * ab_z + g_xy_0_yyzzz_xz[k];

                g_xy_0_yyzzzz_y[k] = -g_xy_0_yyzzz_y[k] * ab_z + g_xy_0_yyzzz_yz[k];

                g_xy_0_yyzzzz_z[k] = -g_xy_0_yyzzz_z[k] * ab_z + g_xy_0_yyzzz_zz[k];
            }

            /// Set up 162-165 components of targeted buffer : cbuffer.data(

            auto g_xy_0_yzzzzz_x = cbuffer.data(ip_geom_20_off + 162 * ccomps * dcomps);

            auto g_xy_0_yzzzzz_y = cbuffer.data(ip_geom_20_off + 163 * ccomps * dcomps);

            auto g_xy_0_yzzzzz_z = cbuffer.data(ip_geom_20_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_yzzzz_x, g_xy_0_yzzzz_xz, g_xy_0_yzzzz_y, g_xy_0_yzzzz_yz, g_xy_0_yzzzz_z, g_xy_0_yzzzz_zz, g_xy_0_yzzzzz_x, g_xy_0_yzzzzz_y, g_xy_0_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_yzzzzz_x[k] = -g_xy_0_yzzzz_x[k] * ab_z + g_xy_0_yzzzz_xz[k];

                g_xy_0_yzzzzz_y[k] = -g_xy_0_yzzzz_y[k] * ab_z + g_xy_0_yzzzz_yz[k];

                g_xy_0_yzzzzz_z[k] = -g_xy_0_yzzzz_z[k] * ab_z + g_xy_0_yzzzz_zz[k];
            }

            /// Set up 165-168 components of targeted buffer : cbuffer.data(

            auto g_xy_0_zzzzzz_x = cbuffer.data(ip_geom_20_off + 165 * ccomps * dcomps);

            auto g_xy_0_zzzzzz_y = cbuffer.data(ip_geom_20_off + 166 * ccomps * dcomps);

            auto g_xy_0_zzzzzz_z = cbuffer.data(ip_geom_20_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_xy_0_zzzzz_x, g_xy_0_zzzzz_xz, g_xy_0_zzzzz_y, g_xy_0_zzzzz_yz, g_xy_0_zzzzz_z, g_xy_0_zzzzz_zz, g_xy_0_zzzzzz_x, g_xy_0_zzzzzz_y, g_xy_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xy_0_zzzzzz_x[k] = -g_xy_0_zzzzz_x[k] * ab_z + g_xy_0_zzzzz_xz[k];

                g_xy_0_zzzzzz_y[k] = -g_xy_0_zzzzz_y[k] * ab_z + g_xy_0_zzzzz_yz[k];

                g_xy_0_zzzzzz_z[k] = -g_xy_0_zzzzz_z[k] * ab_z + g_xy_0_zzzzz_zz[k];
            }

            /// Set up 168-171 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxxx_x = cbuffer.data(ip_geom_20_off + 168 * ccomps * dcomps);

            auto g_xz_0_xxxxxx_y = cbuffer.data(ip_geom_20_off + 169 * ccomps * dcomps);

            auto g_xz_0_xxxxxx_z = cbuffer.data(ip_geom_20_off + 170 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxx_x, g_xz_0_xxxxx_xx, g_xz_0_xxxxx_xy, g_xz_0_xxxxx_xz, g_xz_0_xxxxx_y, g_xz_0_xxxxx_z, g_xz_0_xxxxxx_x, g_xz_0_xxxxxx_y, g_xz_0_xxxxxx_z, g_z_0_xxxxx_x, g_z_0_xxxxx_y, g_z_0_xxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxxx_x[k] = -g_z_0_xxxxx_x[k] - g_xz_0_xxxxx_x[k] * ab_x + g_xz_0_xxxxx_xx[k];

                g_xz_0_xxxxxx_y[k] = -g_z_0_xxxxx_y[k] - g_xz_0_xxxxx_y[k] * ab_x + g_xz_0_xxxxx_xy[k];

                g_xz_0_xxxxxx_z[k] = -g_z_0_xxxxx_z[k] - g_xz_0_xxxxx_z[k] * ab_x + g_xz_0_xxxxx_xz[k];
            }

            /// Set up 171-174 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxxy_x = cbuffer.data(ip_geom_20_off + 171 * ccomps * dcomps);

            auto g_xz_0_xxxxxy_y = cbuffer.data(ip_geom_20_off + 172 * ccomps * dcomps);

            auto g_xz_0_xxxxxy_z = cbuffer.data(ip_geom_20_off + 173 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxx_x, g_xz_0_xxxxx_xy, g_xz_0_xxxxx_y, g_xz_0_xxxxx_yy, g_xz_0_xxxxx_yz, g_xz_0_xxxxx_z, g_xz_0_xxxxxy_x, g_xz_0_xxxxxy_y, g_xz_0_xxxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxxy_x[k] = -g_xz_0_xxxxx_x[k] * ab_y + g_xz_0_xxxxx_xy[k];

                g_xz_0_xxxxxy_y[k] = -g_xz_0_xxxxx_y[k] * ab_y + g_xz_0_xxxxx_yy[k];

                g_xz_0_xxxxxy_z[k] = -g_xz_0_xxxxx_z[k] * ab_y + g_xz_0_xxxxx_yz[k];
            }

            /// Set up 174-177 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxxz_x = cbuffer.data(ip_geom_20_off + 174 * ccomps * dcomps);

            auto g_xz_0_xxxxxz_y = cbuffer.data(ip_geom_20_off + 175 * ccomps * dcomps);

            auto g_xz_0_xxxxxz_z = cbuffer.data(ip_geom_20_off + 176 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxxz_x, g_xz_0_xxxxxz_y, g_xz_0_xxxxxz_z, g_xz_0_xxxxz_x, g_xz_0_xxxxz_xx, g_xz_0_xxxxz_xy, g_xz_0_xxxxz_xz, g_xz_0_xxxxz_y, g_xz_0_xxxxz_z, g_z_0_xxxxz_x, g_z_0_xxxxz_y, g_z_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxxz_x[k] = -g_z_0_xxxxz_x[k] - g_xz_0_xxxxz_x[k] * ab_x + g_xz_0_xxxxz_xx[k];

                g_xz_0_xxxxxz_y[k] = -g_z_0_xxxxz_y[k] - g_xz_0_xxxxz_y[k] * ab_x + g_xz_0_xxxxz_xy[k];

                g_xz_0_xxxxxz_z[k] = -g_z_0_xxxxz_z[k] - g_xz_0_xxxxz_z[k] * ab_x + g_xz_0_xxxxz_xz[k];
            }

            /// Set up 177-180 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxyy_x = cbuffer.data(ip_geom_20_off + 177 * ccomps * dcomps);

            auto g_xz_0_xxxxyy_y = cbuffer.data(ip_geom_20_off + 178 * ccomps * dcomps);

            auto g_xz_0_xxxxyy_z = cbuffer.data(ip_geom_20_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxy_x, g_xz_0_xxxxy_xy, g_xz_0_xxxxy_y, g_xz_0_xxxxy_yy, g_xz_0_xxxxy_yz, g_xz_0_xxxxy_z, g_xz_0_xxxxyy_x, g_xz_0_xxxxyy_y, g_xz_0_xxxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxyy_x[k] = -g_xz_0_xxxxy_x[k] * ab_y + g_xz_0_xxxxy_xy[k];

                g_xz_0_xxxxyy_y[k] = -g_xz_0_xxxxy_y[k] * ab_y + g_xz_0_xxxxy_yy[k];

                g_xz_0_xxxxyy_z[k] = -g_xz_0_xxxxy_z[k] * ab_y + g_xz_0_xxxxy_yz[k];
            }

            /// Set up 180-183 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxyz_x = cbuffer.data(ip_geom_20_off + 180 * ccomps * dcomps);

            auto g_xz_0_xxxxyz_y = cbuffer.data(ip_geom_20_off + 181 * ccomps * dcomps);

            auto g_xz_0_xxxxyz_z = cbuffer.data(ip_geom_20_off + 182 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxyz_x, g_xz_0_xxxxyz_y, g_xz_0_xxxxyz_z, g_xz_0_xxxxz_x, g_xz_0_xxxxz_xy, g_xz_0_xxxxz_y, g_xz_0_xxxxz_yy, g_xz_0_xxxxz_yz, g_xz_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxyz_x[k] = -g_xz_0_xxxxz_x[k] * ab_y + g_xz_0_xxxxz_xy[k];

                g_xz_0_xxxxyz_y[k] = -g_xz_0_xxxxz_y[k] * ab_y + g_xz_0_xxxxz_yy[k];

                g_xz_0_xxxxyz_z[k] = -g_xz_0_xxxxz_z[k] * ab_y + g_xz_0_xxxxz_yz[k];
            }

            /// Set up 183-186 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxxzz_x = cbuffer.data(ip_geom_20_off + 183 * ccomps * dcomps);

            auto g_xz_0_xxxxzz_y = cbuffer.data(ip_geom_20_off + 184 * ccomps * dcomps);

            auto g_xz_0_xxxxzz_z = cbuffer.data(ip_geom_20_off + 185 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxxzz_x, g_xz_0_xxxxzz_y, g_xz_0_xxxxzz_z, g_xz_0_xxxzz_x, g_xz_0_xxxzz_xx, g_xz_0_xxxzz_xy, g_xz_0_xxxzz_xz, g_xz_0_xxxzz_y, g_xz_0_xxxzz_z, g_z_0_xxxzz_x, g_z_0_xxxzz_y, g_z_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxxzz_x[k] = -g_z_0_xxxzz_x[k] - g_xz_0_xxxzz_x[k] * ab_x + g_xz_0_xxxzz_xx[k];

                g_xz_0_xxxxzz_y[k] = -g_z_0_xxxzz_y[k] - g_xz_0_xxxzz_y[k] * ab_x + g_xz_0_xxxzz_xy[k];

                g_xz_0_xxxxzz_z[k] = -g_z_0_xxxzz_z[k] - g_xz_0_xxxzz_z[k] * ab_x + g_xz_0_xxxzz_xz[k];
            }

            /// Set up 186-189 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyyy_x = cbuffer.data(ip_geom_20_off + 186 * ccomps * dcomps);

            auto g_xz_0_xxxyyy_y = cbuffer.data(ip_geom_20_off + 187 * ccomps * dcomps);

            auto g_xz_0_xxxyyy_z = cbuffer.data(ip_geom_20_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyy_x, g_xz_0_xxxyy_xy, g_xz_0_xxxyy_y, g_xz_0_xxxyy_yy, g_xz_0_xxxyy_yz, g_xz_0_xxxyy_z, g_xz_0_xxxyyy_x, g_xz_0_xxxyyy_y, g_xz_0_xxxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyyy_x[k] = -g_xz_0_xxxyy_x[k] * ab_y + g_xz_0_xxxyy_xy[k];

                g_xz_0_xxxyyy_y[k] = -g_xz_0_xxxyy_y[k] * ab_y + g_xz_0_xxxyy_yy[k];

                g_xz_0_xxxyyy_z[k] = -g_xz_0_xxxyy_z[k] * ab_y + g_xz_0_xxxyy_yz[k];
            }

            /// Set up 189-192 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyyz_x = cbuffer.data(ip_geom_20_off + 189 * ccomps * dcomps);

            auto g_xz_0_xxxyyz_y = cbuffer.data(ip_geom_20_off + 190 * ccomps * dcomps);

            auto g_xz_0_xxxyyz_z = cbuffer.data(ip_geom_20_off + 191 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyyz_x, g_xz_0_xxxyyz_y, g_xz_0_xxxyyz_z, g_xz_0_xxxyz_x, g_xz_0_xxxyz_xy, g_xz_0_xxxyz_y, g_xz_0_xxxyz_yy, g_xz_0_xxxyz_yz, g_xz_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyyz_x[k] = -g_xz_0_xxxyz_x[k] * ab_y + g_xz_0_xxxyz_xy[k];

                g_xz_0_xxxyyz_y[k] = -g_xz_0_xxxyz_y[k] * ab_y + g_xz_0_xxxyz_yy[k];

                g_xz_0_xxxyyz_z[k] = -g_xz_0_xxxyz_z[k] * ab_y + g_xz_0_xxxyz_yz[k];
            }

            /// Set up 192-195 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxyzz_x = cbuffer.data(ip_geom_20_off + 192 * ccomps * dcomps);

            auto g_xz_0_xxxyzz_y = cbuffer.data(ip_geom_20_off + 193 * ccomps * dcomps);

            auto g_xz_0_xxxyzz_z = cbuffer.data(ip_geom_20_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxyzz_x, g_xz_0_xxxyzz_y, g_xz_0_xxxyzz_z, g_xz_0_xxxzz_x, g_xz_0_xxxzz_xy, g_xz_0_xxxzz_y, g_xz_0_xxxzz_yy, g_xz_0_xxxzz_yz, g_xz_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxyzz_x[k] = -g_xz_0_xxxzz_x[k] * ab_y + g_xz_0_xxxzz_xy[k];

                g_xz_0_xxxyzz_y[k] = -g_xz_0_xxxzz_y[k] * ab_y + g_xz_0_xxxzz_yy[k];

                g_xz_0_xxxyzz_z[k] = -g_xz_0_xxxzz_z[k] * ab_y + g_xz_0_xxxzz_yz[k];
            }

            /// Set up 195-198 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxxzzz_x = cbuffer.data(ip_geom_20_off + 195 * ccomps * dcomps);

            auto g_xz_0_xxxzzz_y = cbuffer.data(ip_geom_20_off + 196 * ccomps * dcomps);

            auto g_xz_0_xxxzzz_z = cbuffer.data(ip_geom_20_off + 197 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxxzzz_x, g_xz_0_xxxzzz_y, g_xz_0_xxxzzz_z, g_xz_0_xxzzz_x, g_xz_0_xxzzz_xx, g_xz_0_xxzzz_xy, g_xz_0_xxzzz_xz, g_xz_0_xxzzz_y, g_xz_0_xxzzz_z, g_z_0_xxzzz_x, g_z_0_xxzzz_y, g_z_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxxzzz_x[k] = -g_z_0_xxzzz_x[k] - g_xz_0_xxzzz_x[k] * ab_x + g_xz_0_xxzzz_xx[k];

                g_xz_0_xxxzzz_y[k] = -g_z_0_xxzzz_y[k] - g_xz_0_xxzzz_y[k] * ab_x + g_xz_0_xxzzz_xy[k];

                g_xz_0_xxxzzz_z[k] = -g_z_0_xxzzz_z[k] - g_xz_0_xxzzz_z[k] * ab_x + g_xz_0_xxzzz_xz[k];
            }

            /// Set up 198-201 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyyy_x = cbuffer.data(ip_geom_20_off + 198 * ccomps * dcomps);

            auto g_xz_0_xxyyyy_y = cbuffer.data(ip_geom_20_off + 199 * ccomps * dcomps);

            auto g_xz_0_xxyyyy_z = cbuffer.data(ip_geom_20_off + 200 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyy_x, g_xz_0_xxyyy_xy, g_xz_0_xxyyy_y, g_xz_0_xxyyy_yy, g_xz_0_xxyyy_yz, g_xz_0_xxyyy_z, g_xz_0_xxyyyy_x, g_xz_0_xxyyyy_y, g_xz_0_xxyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyyy_x[k] = -g_xz_0_xxyyy_x[k] * ab_y + g_xz_0_xxyyy_xy[k];

                g_xz_0_xxyyyy_y[k] = -g_xz_0_xxyyy_y[k] * ab_y + g_xz_0_xxyyy_yy[k];

                g_xz_0_xxyyyy_z[k] = -g_xz_0_xxyyy_z[k] * ab_y + g_xz_0_xxyyy_yz[k];
            }

            /// Set up 201-204 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyyz_x = cbuffer.data(ip_geom_20_off + 201 * ccomps * dcomps);

            auto g_xz_0_xxyyyz_y = cbuffer.data(ip_geom_20_off + 202 * ccomps * dcomps);

            auto g_xz_0_xxyyyz_z = cbuffer.data(ip_geom_20_off + 203 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyyz_x, g_xz_0_xxyyyz_y, g_xz_0_xxyyyz_z, g_xz_0_xxyyz_x, g_xz_0_xxyyz_xy, g_xz_0_xxyyz_y, g_xz_0_xxyyz_yy, g_xz_0_xxyyz_yz, g_xz_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyyz_x[k] = -g_xz_0_xxyyz_x[k] * ab_y + g_xz_0_xxyyz_xy[k];

                g_xz_0_xxyyyz_y[k] = -g_xz_0_xxyyz_y[k] * ab_y + g_xz_0_xxyyz_yy[k];

                g_xz_0_xxyyyz_z[k] = -g_xz_0_xxyyz_z[k] * ab_y + g_xz_0_xxyyz_yz[k];
            }

            /// Set up 204-207 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyyzz_x = cbuffer.data(ip_geom_20_off + 204 * ccomps * dcomps);

            auto g_xz_0_xxyyzz_y = cbuffer.data(ip_geom_20_off + 205 * ccomps * dcomps);

            auto g_xz_0_xxyyzz_z = cbuffer.data(ip_geom_20_off + 206 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyyzz_x, g_xz_0_xxyyzz_y, g_xz_0_xxyyzz_z, g_xz_0_xxyzz_x, g_xz_0_xxyzz_xy, g_xz_0_xxyzz_y, g_xz_0_xxyzz_yy, g_xz_0_xxyzz_yz, g_xz_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyyzz_x[k] = -g_xz_0_xxyzz_x[k] * ab_y + g_xz_0_xxyzz_xy[k];

                g_xz_0_xxyyzz_y[k] = -g_xz_0_xxyzz_y[k] * ab_y + g_xz_0_xxyzz_yy[k];

                g_xz_0_xxyyzz_z[k] = -g_xz_0_xxyzz_z[k] * ab_y + g_xz_0_xxyzz_yz[k];
            }

            /// Set up 207-210 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxyzzz_x = cbuffer.data(ip_geom_20_off + 207 * ccomps * dcomps);

            auto g_xz_0_xxyzzz_y = cbuffer.data(ip_geom_20_off + 208 * ccomps * dcomps);

            auto g_xz_0_xxyzzz_z = cbuffer.data(ip_geom_20_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxyzzz_x, g_xz_0_xxyzzz_y, g_xz_0_xxyzzz_z, g_xz_0_xxzzz_x, g_xz_0_xxzzz_xy, g_xz_0_xxzzz_y, g_xz_0_xxzzz_yy, g_xz_0_xxzzz_yz, g_xz_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxyzzz_x[k] = -g_xz_0_xxzzz_x[k] * ab_y + g_xz_0_xxzzz_xy[k];

                g_xz_0_xxyzzz_y[k] = -g_xz_0_xxzzz_y[k] * ab_y + g_xz_0_xxzzz_yy[k];

                g_xz_0_xxyzzz_z[k] = -g_xz_0_xxzzz_z[k] * ab_y + g_xz_0_xxzzz_yz[k];
            }

            /// Set up 210-213 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xxzzzz_x = cbuffer.data(ip_geom_20_off + 210 * ccomps * dcomps);

            auto g_xz_0_xxzzzz_y = cbuffer.data(ip_geom_20_off + 211 * ccomps * dcomps);

            auto g_xz_0_xxzzzz_z = cbuffer.data(ip_geom_20_off + 212 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xxzzzz_x, g_xz_0_xxzzzz_y, g_xz_0_xxzzzz_z, g_xz_0_xzzzz_x, g_xz_0_xzzzz_xx, g_xz_0_xzzzz_xy, g_xz_0_xzzzz_xz, g_xz_0_xzzzz_y, g_xz_0_xzzzz_z, g_z_0_xzzzz_x, g_z_0_xzzzz_y, g_z_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xxzzzz_x[k] = -g_z_0_xzzzz_x[k] - g_xz_0_xzzzz_x[k] * ab_x + g_xz_0_xzzzz_xx[k];

                g_xz_0_xxzzzz_y[k] = -g_z_0_xzzzz_y[k] - g_xz_0_xzzzz_y[k] * ab_x + g_xz_0_xzzzz_xy[k];

                g_xz_0_xxzzzz_z[k] = -g_z_0_xzzzz_z[k] - g_xz_0_xzzzz_z[k] * ab_x + g_xz_0_xzzzz_xz[k];
            }

            /// Set up 213-216 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyyy_x = cbuffer.data(ip_geom_20_off + 213 * ccomps * dcomps);

            auto g_xz_0_xyyyyy_y = cbuffer.data(ip_geom_20_off + 214 * ccomps * dcomps);

            auto g_xz_0_xyyyyy_z = cbuffer.data(ip_geom_20_off + 215 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyy_x, g_xz_0_xyyyy_xy, g_xz_0_xyyyy_y, g_xz_0_xyyyy_yy, g_xz_0_xyyyy_yz, g_xz_0_xyyyy_z, g_xz_0_xyyyyy_x, g_xz_0_xyyyyy_y, g_xz_0_xyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyyy_x[k] = -g_xz_0_xyyyy_x[k] * ab_y + g_xz_0_xyyyy_xy[k];

                g_xz_0_xyyyyy_y[k] = -g_xz_0_xyyyy_y[k] * ab_y + g_xz_0_xyyyy_yy[k];

                g_xz_0_xyyyyy_z[k] = -g_xz_0_xyyyy_z[k] * ab_y + g_xz_0_xyyyy_yz[k];
            }

            /// Set up 216-219 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyyz_x = cbuffer.data(ip_geom_20_off + 216 * ccomps * dcomps);

            auto g_xz_0_xyyyyz_y = cbuffer.data(ip_geom_20_off + 217 * ccomps * dcomps);

            auto g_xz_0_xyyyyz_z = cbuffer.data(ip_geom_20_off + 218 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyyz_x, g_xz_0_xyyyyz_y, g_xz_0_xyyyyz_z, g_xz_0_xyyyz_x, g_xz_0_xyyyz_xy, g_xz_0_xyyyz_y, g_xz_0_xyyyz_yy, g_xz_0_xyyyz_yz, g_xz_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyyz_x[k] = -g_xz_0_xyyyz_x[k] * ab_y + g_xz_0_xyyyz_xy[k];

                g_xz_0_xyyyyz_y[k] = -g_xz_0_xyyyz_y[k] * ab_y + g_xz_0_xyyyz_yy[k];

                g_xz_0_xyyyyz_z[k] = -g_xz_0_xyyyz_z[k] * ab_y + g_xz_0_xyyyz_yz[k];
            }

            /// Set up 219-222 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyyzz_x = cbuffer.data(ip_geom_20_off + 219 * ccomps * dcomps);

            auto g_xz_0_xyyyzz_y = cbuffer.data(ip_geom_20_off + 220 * ccomps * dcomps);

            auto g_xz_0_xyyyzz_z = cbuffer.data(ip_geom_20_off + 221 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyyzz_x, g_xz_0_xyyyzz_y, g_xz_0_xyyyzz_z, g_xz_0_xyyzz_x, g_xz_0_xyyzz_xy, g_xz_0_xyyzz_y, g_xz_0_xyyzz_yy, g_xz_0_xyyzz_yz, g_xz_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyyzz_x[k] = -g_xz_0_xyyzz_x[k] * ab_y + g_xz_0_xyyzz_xy[k];

                g_xz_0_xyyyzz_y[k] = -g_xz_0_xyyzz_y[k] * ab_y + g_xz_0_xyyzz_yy[k];

                g_xz_0_xyyyzz_z[k] = -g_xz_0_xyyzz_z[k] * ab_y + g_xz_0_xyyzz_yz[k];
            }

            /// Set up 222-225 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyyzzz_x = cbuffer.data(ip_geom_20_off + 222 * ccomps * dcomps);

            auto g_xz_0_xyyzzz_y = cbuffer.data(ip_geom_20_off + 223 * ccomps * dcomps);

            auto g_xz_0_xyyzzz_z = cbuffer.data(ip_geom_20_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyyzzz_x, g_xz_0_xyyzzz_y, g_xz_0_xyyzzz_z, g_xz_0_xyzzz_x, g_xz_0_xyzzz_xy, g_xz_0_xyzzz_y, g_xz_0_xyzzz_yy, g_xz_0_xyzzz_yz, g_xz_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyyzzz_x[k] = -g_xz_0_xyzzz_x[k] * ab_y + g_xz_0_xyzzz_xy[k];

                g_xz_0_xyyzzz_y[k] = -g_xz_0_xyzzz_y[k] * ab_y + g_xz_0_xyzzz_yy[k];

                g_xz_0_xyyzzz_z[k] = -g_xz_0_xyzzz_z[k] * ab_y + g_xz_0_xyzzz_yz[k];
            }

            /// Set up 225-228 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xyzzzz_x = cbuffer.data(ip_geom_20_off + 225 * ccomps * dcomps);

            auto g_xz_0_xyzzzz_y = cbuffer.data(ip_geom_20_off + 226 * ccomps * dcomps);

            auto g_xz_0_xyzzzz_z = cbuffer.data(ip_geom_20_off + 227 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xyzzzz_x, g_xz_0_xyzzzz_y, g_xz_0_xyzzzz_z, g_xz_0_xzzzz_x, g_xz_0_xzzzz_xy, g_xz_0_xzzzz_y, g_xz_0_xzzzz_yy, g_xz_0_xzzzz_yz, g_xz_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xyzzzz_x[k] = -g_xz_0_xzzzz_x[k] * ab_y + g_xz_0_xzzzz_xy[k];

                g_xz_0_xyzzzz_y[k] = -g_xz_0_xzzzz_y[k] * ab_y + g_xz_0_xzzzz_yy[k];

                g_xz_0_xyzzzz_z[k] = -g_xz_0_xzzzz_z[k] * ab_y + g_xz_0_xzzzz_yz[k];
            }

            /// Set up 228-231 components of targeted buffer : cbuffer.data(

            auto g_xz_0_xzzzzz_x = cbuffer.data(ip_geom_20_off + 228 * ccomps * dcomps);

            auto g_xz_0_xzzzzz_y = cbuffer.data(ip_geom_20_off + 229 * ccomps * dcomps);

            auto g_xz_0_xzzzzz_z = cbuffer.data(ip_geom_20_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_xzzzzz_x, g_xz_0_xzzzzz_y, g_xz_0_xzzzzz_z, g_xz_0_zzzzz_x, g_xz_0_zzzzz_xx, g_xz_0_zzzzz_xy, g_xz_0_zzzzz_xz, g_xz_0_zzzzz_y, g_xz_0_zzzzz_z, g_z_0_zzzzz_x, g_z_0_zzzzz_y, g_z_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_xzzzzz_x[k] = -g_z_0_zzzzz_x[k] - g_xz_0_zzzzz_x[k] * ab_x + g_xz_0_zzzzz_xx[k];

                g_xz_0_xzzzzz_y[k] = -g_z_0_zzzzz_y[k] - g_xz_0_zzzzz_y[k] * ab_x + g_xz_0_zzzzz_xy[k];

                g_xz_0_xzzzzz_z[k] = -g_z_0_zzzzz_z[k] - g_xz_0_zzzzz_z[k] * ab_x + g_xz_0_zzzzz_xz[k];
            }

            /// Set up 231-234 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyyy_x = cbuffer.data(ip_geom_20_off + 231 * ccomps * dcomps);

            auto g_xz_0_yyyyyy_y = cbuffer.data(ip_geom_20_off + 232 * ccomps * dcomps);

            auto g_xz_0_yyyyyy_z = cbuffer.data(ip_geom_20_off + 233 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyy_x, g_xz_0_yyyyy_xy, g_xz_0_yyyyy_y, g_xz_0_yyyyy_yy, g_xz_0_yyyyy_yz, g_xz_0_yyyyy_z, g_xz_0_yyyyyy_x, g_xz_0_yyyyyy_y, g_xz_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyyy_x[k] = -g_xz_0_yyyyy_x[k] * ab_y + g_xz_0_yyyyy_xy[k];

                g_xz_0_yyyyyy_y[k] = -g_xz_0_yyyyy_y[k] * ab_y + g_xz_0_yyyyy_yy[k];

                g_xz_0_yyyyyy_z[k] = -g_xz_0_yyyyy_z[k] * ab_y + g_xz_0_yyyyy_yz[k];
            }

            /// Set up 234-237 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyyz_x = cbuffer.data(ip_geom_20_off + 234 * ccomps * dcomps);

            auto g_xz_0_yyyyyz_y = cbuffer.data(ip_geom_20_off + 235 * ccomps * dcomps);

            auto g_xz_0_yyyyyz_z = cbuffer.data(ip_geom_20_off + 236 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyyz_x, g_xz_0_yyyyyz_y, g_xz_0_yyyyyz_z, g_xz_0_yyyyz_x, g_xz_0_yyyyz_xy, g_xz_0_yyyyz_y, g_xz_0_yyyyz_yy, g_xz_0_yyyyz_yz, g_xz_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyyz_x[k] = -g_xz_0_yyyyz_x[k] * ab_y + g_xz_0_yyyyz_xy[k];

                g_xz_0_yyyyyz_y[k] = -g_xz_0_yyyyz_y[k] * ab_y + g_xz_0_yyyyz_yy[k];

                g_xz_0_yyyyyz_z[k] = -g_xz_0_yyyyz_z[k] * ab_y + g_xz_0_yyyyz_yz[k];
            }

            /// Set up 237-240 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyyzz_x = cbuffer.data(ip_geom_20_off + 237 * ccomps * dcomps);

            auto g_xz_0_yyyyzz_y = cbuffer.data(ip_geom_20_off + 238 * ccomps * dcomps);

            auto g_xz_0_yyyyzz_z = cbuffer.data(ip_geom_20_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyyzz_x, g_xz_0_yyyyzz_y, g_xz_0_yyyyzz_z, g_xz_0_yyyzz_x, g_xz_0_yyyzz_xy, g_xz_0_yyyzz_y, g_xz_0_yyyzz_yy, g_xz_0_yyyzz_yz, g_xz_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyyzz_x[k] = -g_xz_0_yyyzz_x[k] * ab_y + g_xz_0_yyyzz_xy[k];

                g_xz_0_yyyyzz_y[k] = -g_xz_0_yyyzz_y[k] * ab_y + g_xz_0_yyyzz_yy[k];

                g_xz_0_yyyyzz_z[k] = -g_xz_0_yyyzz_z[k] * ab_y + g_xz_0_yyyzz_yz[k];
            }

            /// Set up 240-243 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyyzzz_x = cbuffer.data(ip_geom_20_off + 240 * ccomps * dcomps);

            auto g_xz_0_yyyzzz_y = cbuffer.data(ip_geom_20_off + 241 * ccomps * dcomps);

            auto g_xz_0_yyyzzz_z = cbuffer.data(ip_geom_20_off + 242 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyyzzz_x, g_xz_0_yyyzzz_y, g_xz_0_yyyzzz_z, g_xz_0_yyzzz_x, g_xz_0_yyzzz_xy, g_xz_0_yyzzz_y, g_xz_0_yyzzz_yy, g_xz_0_yyzzz_yz, g_xz_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyyzzz_x[k] = -g_xz_0_yyzzz_x[k] * ab_y + g_xz_0_yyzzz_xy[k];

                g_xz_0_yyyzzz_y[k] = -g_xz_0_yyzzz_y[k] * ab_y + g_xz_0_yyzzz_yy[k];

                g_xz_0_yyyzzz_z[k] = -g_xz_0_yyzzz_z[k] * ab_y + g_xz_0_yyzzz_yz[k];
            }

            /// Set up 243-246 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yyzzzz_x = cbuffer.data(ip_geom_20_off + 243 * ccomps * dcomps);

            auto g_xz_0_yyzzzz_y = cbuffer.data(ip_geom_20_off + 244 * ccomps * dcomps);

            auto g_xz_0_yyzzzz_z = cbuffer.data(ip_geom_20_off + 245 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yyzzzz_x, g_xz_0_yyzzzz_y, g_xz_0_yyzzzz_z, g_xz_0_yzzzz_x, g_xz_0_yzzzz_xy, g_xz_0_yzzzz_y, g_xz_0_yzzzz_yy, g_xz_0_yzzzz_yz, g_xz_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yyzzzz_x[k] = -g_xz_0_yzzzz_x[k] * ab_y + g_xz_0_yzzzz_xy[k];

                g_xz_0_yyzzzz_y[k] = -g_xz_0_yzzzz_y[k] * ab_y + g_xz_0_yzzzz_yy[k];

                g_xz_0_yyzzzz_z[k] = -g_xz_0_yzzzz_z[k] * ab_y + g_xz_0_yzzzz_yz[k];
            }

            /// Set up 246-249 components of targeted buffer : cbuffer.data(

            auto g_xz_0_yzzzzz_x = cbuffer.data(ip_geom_20_off + 246 * ccomps * dcomps);

            auto g_xz_0_yzzzzz_y = cbuffer.data(ip_geom_20_off + 247 * ccomps * dcomps);

            auto g_xz_0_yzzzzz_z = cbuffer.data(ip_geom_20_off + 248 * ccomps * dcomps);

            #pragma omp simd aligned(g_xz_0_yzzzzz_x, g_xz_0_yzzzzz_y, g_xz_0_yzzzzz_z, g_xz_0_zzzzz_x, g_xz_0_zzzzz_xy, g_xz_0_zzzzz_y, g_xz_0_zzzzz_yy, g_xz_0_zzzzz_yz, g_xz_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_yzzzzz_x[k] = -g_xz_0_zzzzz_x[k] * ab_y + g_xz_0_zzzzz_xy[k];

                g_xz_0_yzzzzz_y[k] = -g_xz_0_zzzzz_y[k] * ab_y + g_xz_0_zzzzz_yy[k];

                g_xz_0_yzzzzz_z[k] = -g_xz_0_zzzzz_z[k] * ab_y + g_xz_0_zzzzz_yz[k];
            }

            /// Set up 249-252 components of targeted buffer : cbuffer.data(

            auto g_xz_0_zzzzzz_x = cbuffer.data(ip_geom_20_off + 249 * ccomps * dcomps);

            auto g_xz_0_zzzzzz_y = cbuffer.data(ip_geom_20_off + 250 * ccomps * dcomps);

            auto g_xz_0_zzzzzz_z = cbuffer.data(ip_geom_20_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzzz_x, g_x_0_zzzzz_y, g_x_0_zzzzz_z, g_xz_0_zzzzz_x, g_xz_0_zzzzz_xz, g_xz_0_zzzzz_y, g_xz_0_zzzzz_yz, g_xz_0_zzzzz_z, g_xz_0_zzzzz_zz, g_xz_0_zzzzzz_x, g_xz_0_zzzzzz_y, g_xz_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_xz_0_zzzzzz_x[k] = -g_x_0_zzzzz_x[k] - g_xz_0_zzzzz_x[k] * ab_z + g_xz_0_zzzzz_xz[k];

                g_xz_0_zzzzzz_y[k] = -g_x_0_zzzzz_y[k] - g_xz_0_zzzzz_y[k] * ab_z + g_xz_0_zzzzz_yz[k];

                g_xz_0_zzzzzz_z[k] = -g_x_0_zzzzz_z[k] - g_xz_0_zzzzz_z[k] * ab_z + g_xz_0_zzzzz_zz[k];
            }

            /// Set up 252-255 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxxx_x = cbuffer.data(ip_geom_20_off + 252 * ccomps * dcomps);

            auto g_yy_0_xxxxxx_y = cbuffer.data(ip_geom_20_off + 253 * ccomps * dcomps);

            auto g_yy_0_xxxxxx_z = cbuffer.data(ip_geom_20_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxx_x, g_yy_0_xxxxx_xx, g_yy_0_xxxxx_xy, g_yy_0_xxxxx_xz, g_yy_0_xxxxx_y, g_yy_0_xxxxx_z, g_yy_0_xxxxxx_x, g_yy_0_xxxxxx_y, g_yy_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxxx_x[k] = -g_yy_0_xxxxx_x[k] * ab_x + g_yy_0_xxxxx_xx[k];

                g_yy_0_xxxxxx_y[k] = -g_yy_0_xxxxx_y[k] * ab_x + g_yy_0_xxxxx_xy[k];

                g_yy_0_xxxxxx_z[k] = -g_yy_0_xxxxx_z[k] * ab_x + g_yy_0_xxxxx_xz[k];
            }

            /// Set up 255-258 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxxy_x = cbuffer.data(ip_geom_20_off + 255 * ccomps * dcomps);

            auto g_yy_0_xxxxxy_y = cbuffer.data(ip_geom_20_off + 256 * ccomps * dcomps);

            auto g_yy_0_xxxxxy_z = cbuffer.data(ip_geom_20_off + 257 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxxy_x, g_yy_0_xxxxxy_y, g_yy_0_xxxxxy_z, g_yy_0_xxxxy_x, g_yy_0_xxxxy_xx, g_yy_0_xxxxy_xy, g_yy_0_xxxxy_xz, g_yy_0_xxxxy_y, g_yy_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxxy_x[k] = -g_yy_0_xxxxy_x[k] * ab_x + g_yy_0_xxxxy_xx[k];

                g_yy_0_xxxxxy_y[k] = -g_yy_0_xxxxy_y[k] * ab_x + g_yy_0_xxxxy_xy[k];

                g_yy_0_xxxxxy_z[k] = -g_yy_0_xxxxy_z[k] * ab_x + g_yy_0_xxxxy_xz[k];
            }

            /// Set up 258-261 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxxz_x = cbuffer.data(ip_geom_20_off + 258 * ccomps * dcomps);

            auto g_yy_0_xxxxxz_y = cbuffer.data(ip_geom_20_off + 259 * ccomps * dcomps);

            auto g_yy_0_xxxxxz_z = cbuffer.data(ip_geom_20_off + 260 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxxz_x, g_yy_0_xxxxxz_y, g_yy_0_xxxxxz_z, g_yy_0_xxxxz_x, g_yy_0_xxxxz_xx, g_yy_0_xxxxz_xy, g_yy_0_xxxxz_xz, g_yy_0_xxxxz_y, g_yy_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxxz_x[k] = -g_yy_0_xxxxz_x[k] * ab_x + g_yy_0_xxxxz_xx[k];

                g_yy_0_xxxxxz_y[k] = -g_yy_0_xxxxz_y[k] * ab_x + g_yy_0_xxxxz_xy[k];

                g_yy_0_xxxxxz_z[k] = -g_yy_0_xxxxz_z[k] * ab_x + g_yy_0_xxxxz_xz[k];
            }

            /// Set up 261-264 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxyy_x = cbuffer.data(ip_geom_20_off + 261 * ccomps * dcomps);

            auto g_yy_0_xxxxyy_y = cbuffer.data(ip_geom_20_off + 262 * ccomps * dcomps);

            auto g_yy_0_xxxxyy_z = cbuffer.data(ip_geom_20_off + 263 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxyy_x, g_yy_0_xxxxyy_y, g_yy_0_xxxxyy_z, g_yy_0_xxxyy_x, g_yy_0_xxxyy_xx, g_yy_0_xxxyy_xy, g_yy_0_xxxyy_xz, g_yy_0_xxxyy_y, g_yy_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxyy_x[k] = -g_yy_0_xxxyy_x[k] * ab_x + g_yy_0_xxxyy_xx[k];

                g_yy_0_xxxxyy_y[k] = -g_yy_0_xxxyy_y[k] * ab_x + g_yy_0_xxxyy_xy[k];

                g_yy_0_xxxxyy_z[k] = -g_yy_0_xxxyy_z[k] * ab_x + g_yy_0_xxxyy_xz[k];
            }

            /// Set up 264-267 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxyz_x = cbuffer.data(ip_geom_20_off + 264 * ccomps * dcomps);

            auto g_yy_0_xxxxyz_y = cbuffer.data(ip_geom_20_off + 265 * ccomps * dcomps);

            auto g_yy_0_xxxxyz_z = cbuffer.data(ip_geom_20_off + 266 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxyz_x, g_yy_0_xxxxyz_y, g_yy_0_xxxxyz_z, g_yy_0_xxxyz_x, g_yy_0_xxxyz_xx, g_yy_0_xxxyz_xy, g_yy_0_xxxyz_xz, g_yy_0_xxxyz_y, g_yy_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxyz_x[k] = -g_yy_0_xxxyz_x[k] * ab_x + g_yy_0_xxxyz_xx[k];

                g_yy_0_xxxxyz_y[k] = -g_yy_0_xxxyz_y[k] * ab_x + g_yy_0_xxxyz_xy[k];

                g_yy_0_xxxxyz_z[k] = -g_yy_0_xxxyz_z[k] * ab_x + g_yy_0_xxxyz_xz[k];
            }

            /// Set up 267-270 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxxzz_x = cbuffer.data(ip_geom_20_off + 267 * ccomps * dcomps);

            auto g_yy_0_xxxxzz_y = cbuffer.data(ip_geom_20_off + 268 * ccomps * dcomps);

            auto g_yy_0_xxxxzz_z = cbuffer.data(ip_geom_20_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxxzz_x, g_yy_0_xxxxzz_y, g_yy_0_xxxxzz_z, g_yy_0_xxxzz_x, g_yy_0_xxxzz_xx, g_yy_0_xxxzz_xy, g_yy_0_xxxzz_xz, g_yy_0_xxxzz_y, g_yy_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxxzz_x[k] = -g_yy_0_xxxzz_x[k] * ab_x + g_yy_0_xxxzz_xx[k];

                g_yy_0_xxxxzz_y[k] = -g_yy_0_xxxzz_y[k] * ab_x + g_yy_0_xxxzz_xy[k];

                g_yy_0_xxxxzz_z[k] = -g_yy_0_xxxzz_z[k] * ab_x + g_yy_0_xxxzz_xz[k];
            }

            /// Set up 270-273 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyyy_x = cbuffer.data(ip_geom_20_off + 270 * ccomps * dcomps);

            auto g_yy_0_xxxyyy_y = cbuffer.data(ip_geom_20_off + 271 * ccomps * dcomps);

            auto g_yy_0_xxxyyy_z = cbuffer.data(ip_geom_20_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyyy_x, g_yy_0_xxxyyy_y, g_yy_0_xxxyyy_z, g_yy_0_xxyyy_x, g_yy_0_xxyyy_xx, g_yy_0_xxyyy_xy, g_yy_0_xxyyy_xz, g_yy_0_xxyyy_y, g_yy_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyyy_x[k] = -g_yy_0_xxyyy_x[k] * ab_x + g_yy_0_xxyyy_xx[k];

                g_yy_0_xxxyyy_y[k] = -g_yy_0_xxyyy_y[k] * ab_x + g_yy_0_xxyyy_xy[k];

                g_yy_0_xxxyyy_z[k] = -g_yy_0_xxyyy_z[k] * ab_x + g_yy_0_xxyyy_xz[k];
            }

            /// Set up 273-276 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyyz_x = cbuffer.data(ip_geom_20_off + 273 * ccomps * dcomps);

            auto g_yy_0_xxxyyz_y = cbuffer.data(ip_geom_20_off + 274 * ccomps * dcomps);

            auto g_yy_0_xxxyyz_z = cbuffer.data(ip_geom_20_off + 275 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyyz_x, g_yy_0_xxxyyz_y, g_yy_0_xxxyyz_z, g_yy_0_xxyyz_x, g_yy_0_xxyyz_xx, g_yy_0_xxyyz_xy, g_yy_0_xxyyz_xz, g_yy_0_xxyyz_y, g_yy_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyyz_x[k] = -g_yy_0_xxyyz_x[k] * ab_x + g_yy_0_xxyyz_xx[k];

                g_yy_0_xxxyyz_y[k] = -g_yy_0_xxyyz_y[k] * ab_x + g_yy_0_xxyyz_xy[k];

                g_yy_0_xxxyyz_z[k] = -g_yy_0_xxyyz_z[k] * ab_x + g_yy_0_xxyyz_xz[k];
            }

            /// Set up 276-279 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxyzz_x = cbuffer.data(ip_geom_20_off + 276 * ccomps * dcomps);

            auto g_yy_0_xxxyzz_y = cbuffer.data(ip_geom_20_off + 277 * ccomps * dcomps);

            auto g_yy_0_xxxyzz_z = cbuffer.data(ip_geom_20_off + 278 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxyzz_x, g_yy_0_xxxyzz_y, g_yy_0_xxxyzz_z, g_yy_0_xxyzz_x, g_yy_0_xxyzz_xx, g_yy_0_xxyzz_xy, g_yy_0_xxyzz_xz, g_yy_0_xxyzz_y, g_yy_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxyzz_x[k] = -g_yy_0_xxyzz_x[k] * ab_x + g_yy_0_xxyzz_xx[k];

                g_yy_0_xxxyzz_y[k] = -g_yy_0_xxyzz_y[k] * ab_x + g_yy_0_xxyzz_xy[k];

                g_yy_0_xxxyzz_z[k] = -g_yy_0_xxyzz_z[k] * ab_x + g_yy_0_xxyzz_xz[k];
            }

            /// Set up 279-282 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxxzzz_x = cbuffer.data(ip_geom_20_off + 279 * ccomps * dcomps);

            auto g_yy_0_xxxzzz_y = cbuffer.data(ip_geom_20_off + 280 * ccomps * dcomps);

            auto g_yy_0_xxxzzz_z = cbuffer.data(ip_geom_20_off + 281 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxxzzz_x, g_yy_0_xxxzzz_y, g_yy_0_xxxzzz_z, g_yy_0_xxzzz_x, g_yy_0_xxzzz_xx, g_yy_0_xxzzz_xy, g_yy_0_xxzzz_xz, g_yy_0_xxzzz_y, g_yy_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxxzzz_x[k] = -g_yy_0_xxzzz_x[k] * ab_x + g_yy_0_xxzzz_xx[k];

                g_yy_0_xxxzzz_y[k] = -g_yy_0_xxzzz_y[k] * ab_x + g_yy_0_xxzzz_xy[k];

                g_yy_0_xxxzzz_z[k] = -g_yy_0_xxzzz_z[k] * ab_x + g_yy_0_xxzzz_xz[k];
            }

            /// Set up 282-285 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyyy_x = cbuffer.data(ip_geom_20_off + 282 * ccomps * dcomps);

            auto g_yy_0_xxyyyy_y = cbuffer.data(ip_geom_20_off + 283 * ccomps * dcomps);

            auto g_yy_0_xxyyyy_z = cbuffer.data(ip_geom_20_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyyy_x, g_yy_0_xxyyyy_y, g_yy_0_xxyyyy_z, g_yy_0_xyyyy_x, g_yy_0_xyyyy_xx, g_yy_0_xyyyy_xy, g_yy_0_xyyyy_xz, g_yy_0_xyyyy_y, g_yy_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyyy_x[k] = -g_yy_0_xyyyy_x[k] * ab_x + g_yy_0_xyyyy_xx[k];

                g_yy_0_xxyyyy_y[k] = -g_yy_0_xyyyy_y[k] * ab_x + g_yy_0_xyyyy_xy[k];

                g_yy_0_xxyyyy_z[k] = -g_yy_0_xyyyy_z[k] * ab_x + g_yy_0_xyyyy_xz[k];
            }

            /// Set up 285-288 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyyz_x = cbuffer.data(ip_geom_20_off + 285 * ccomps * dcomps);

            auto g_yy_0_xxyyyz_y = cbuffer.data(ip_geom_20_off + 286 * ccomps * dcomps);

            auto g_yy_0_xxyyyz_z = cbuffer.data(ip_geom_20_off + 287 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyyz_x, g_yy_0_xxyyyz_y, g_yy_0_xxyyyz_z, g_yy_0_xyyyz_x, g_yy_0_xyyyz_xx, g_yy_0_xyyyz_xy, g_yy_0_xyyyz_xz, g_yy_0_xyyyz_y, g_yy_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyyz_x[k] = -g_yy_0_xyyyz_x[k] * ab_x + g_yy_0_xyyyz_xx[k];

                g_yy_0_xxyyyz_y[k] = -g_yy_0_xyyyz_y[k] * ab_x + g_yy_0_xyyyz_xy[k];

                g_yy_0_xxyyyz_z[k] = -g_yy_0_xyyyz_z[k] * ab_x + g_yy_0_xyyyz_xz[k];
            }

            /// Set up 288-291 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyyzz_x = cbuffer.data(ip_geom_20_off + 288 * ccomps * dcomps);

            auto g_yy_0_xxyyzz_y = cbuffer.data(ip_geom_20_off + 289 * ccomps * dcomps);

            auto g_yy_0_xxyyzz_z = cbuffer.data(ip_geom_20_off + 290 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyyzz_x, g_yy_0_xxyyzz_y, g_yy_0_xxyyzz_z, g_yy_0_xyyzz_x, g_yy_0_xyyzz_xx, g_yy_0_xyyzz_xy, g_yy_0_xyyzz_xz, g_yy_0_xyyzz_y, g_yy_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyyzz_x[k] = -g_yy_0_xyyzz_x[k] * ab_x + g_yy_0_xyyzz_xx[k];

                g_yy_0_xxyyzz_y[k] = -g_yy_0_xyyzz_y[k] * ab_x + g_yy_0_xyyzz_xy[k];

                g_yy_0_xxyyzz_z[k] = -g_yy_0_xyyzz_z[k] * ab_x + g_yy_0_xyyzz_xz[k];
            }

            /// Set up 291-294 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxyzzz_x = cbuffer.data(ip_geom_20_off + 291 * ccomps * dcomps);

            auto g_yy_0_xxyzzz_y = cbuffer.data(ip_geom_20_off + 292 * ccomps * dcomps);

            auto g_yy_0_xxyzzz_z = cbuffer.data(ip_geom_20_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxyzzz_x, g_yy_0_xxyzzz_y, g_yy_0_xxyzzz_z, g_yy_0_xyzzz_x, g_yy_0_xyzzz_xx, g_yy_0_xyzzz_xy, g_yy_0_xyzzz_xz, g_yy_0_xyzzz_y, g_yy_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxyzzz_x[k] = -g_yy_0_xyzzz_x[k] * ab_x + g_yy_0_xyzzz_xx[k];

                g_yy_0_xxyzzz_y[k] = -g_yy_0_xyzzz_y[k] * ab_x + g_yy_0_xyzzz_xy[k];

                g_yy_0_xxyzzz_z[k] = -g_yy_0_xyzzz_z[k] * ab_x + g_yy_0_xyzzz_xz[k];
            }

            /// Set up 294-297 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xxzzzz_x = cbuffer.data(ip_geom_20_off + 294 * ccomps * dcomps);

            auto g_yy_0_xxzzzz_y = cbuffer.data(ip_geom_20_off + 295 * ccomps * dcomps);

            auto g_yy_0_xxzzzz_z = cbuffer.data(ip_geom_20_off + 296 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xxzzzz_x, g_yy_0_xxzzzz_y, g_yy_0_xxzzzz_z, g_yy_0_xzzzz_x, g_yy_0_xzzzz_xx, g_yy_0_xzzzz_xy, g_yy_0_xzzzz_xz, g_yy_0_xzzzz_y, g_yy_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xxzzzz_x[k] = -g_yy_0_xzzzz_x[k] * ab_x + g_yy_0_xzzzz_xx[k];

                g_yy_0_xxzzzz_y[k] = -g_yy_0_xzzzz_y[k] * ab_x + g_yy_0_xzzzz_xy[k];

                g_yy_0_xxzzzz_z[k] = -g_yy_0_xzzzz_z[k] * ab_x + g_yy_0_xzzzz_xz[k];
            }

            /// Set up 297-300 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyyy_x = cbuffer.data(ip_geom_20_off + 297 * ccomps * dcomps);

            auto g_yy_0_xyyyyy_y = cbuffer.data(ip_geom_20_off + 298 * ccomps * dcomps);

            auto g_yy_0_xyyyyy_z = cbuffer.data(ip_geom_20_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyyy_x, g_yy_0_xyyyyy_y, g_yy_0_xyyyyy_z, g_yy_0_yyyyy_x, g_yy_0_yyyyy_xx, g_yy_0_yyyyy_xy, g_yy_0_yyyyy_xz, g_yy_0_yyyyy_y, g_yy_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyyy_x[k] = -g_yy_0_yyyyy_x[k] * ab_x + g_yy_0_yyyyy_xx[k];

                g_yy_0_xyyyyy_y[k] = -g_yy_0_yyyyy_y[k] * ab_x + g_yy_0_yyyyy_xy[k];

                g_yy_0_xyyyyy_z[k] = -g_yy_0_yyyyy_z[k] * ab_x + g_yy_0_yyyyy_xz[k];
            }

            /// Set up 300-303 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyyz_x = cbuffer.data(ip_geom_20_off + 300 * ccomps * dcomps);

            auto g_yy_0_xyyyyz_y = cbuffer.data(ip_geom_20_off + 301 * ccomps * dcomps);

            auto g_yy_0_xyyyyz_z = cbuffer.data(ip_geom_20_off + 302 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyyz_x, g_yy_0_xyyyyz_y, g_yy_0_xyyyyz_z, g_yy_0_yyyyz_x, g_yy_0_yyyyz_xx, g_yy_0_yyyyz_xy, g_yy_0_yyyyz_xz, g_yy_0_yyyyz_y, g_yy_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyyz_x[k] = -g_yy_0_yyyyz_x[k] * ab_x + g_yy_0_yyyyz_xx[k];

                g_yy_0_xyyyyz_y[k] = -g_yy_0_yyyyz_y[k] * ab_x + g_yy_0_yyyyz_xy[k];

                g_yy_0_xyyyyz_z[k] = -g_yy_0_yyyyz_z[k] * ab_x + g_yy_0_yyyyz_xz[k];
            }

            /// Set up 303-306 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyyzz_x = cbuffer.data(ip_geom_20_off + 303 * ccomps * dcomps);

            auto g_yy_0_xyyyzz_y = cbuffer.data(ip_geom_20_off + 304 * ccomps * dcomps);

            auto g_yy_0_xyyyzz_z = cbuffer.data(ip_geom_20_off + 305 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyyzz_x, g_yy_0_xyyyzz_y, g_yy_0_xyyyzz_z, g_yy_0_yyyzz_x, g_yy_0_yyyzz_xx, g_yy_0_yyyzz_xy, g_yy_0_yyyzz_xz, g_yy_0_yyyzz_y, g_yy_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyyzz_x[k] = -g_yy_0_yyyzz_x[k] * ab_x + g_yy_0_yyyzz_xx[k];

                g_yy_0_xyyyzz_y[k] = -g_yy_0_yyyzz_y[k] * ab_x + g_yy_0_yyyzz_xy[k];

                g_yy_0_xyyyzz_z[k] = -g_yy_0_yyyzz_z[k] * ab_x + g_yy_0_yyyzz_xz[k];
            }

            /// Set up 306-309 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyyzzz_x = cbuffer.data(ip_geom_20_off + 306 * ccomps * dcomps);

            auto g_yy_0_xyyzzz_y = cbuffer.data(ip_geom_20_off + 307 * ccomps * dcomps);

            auto g_yy_0_xyyzzz_z = cbuffer.data(ip_geom_20_off + 308 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyyzzz_x, g_yy_0_xyyzzz_y, g_yy_0_xyyzzz_z, g_yy_0_yyzzz_x, g_yy_0_yyzzz_xx, g_yy_0_yyzzz_xy, g_yy_0_yyzzz_xz, g_yy_0_yyzzz_y, g_yy_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyyzzz_x[k] = -g_yy_0_yyzzz_x[k] * ab_x + g_yy_0_yyzzz_xx[k];

                g_yy_0_xyyzzz_y[k] = -g_yy_0_yyzzz_y[k] * ab_x + g_yy_0_yyzzz_xy[k];

                g_yy_0_xyyzzz_z[k] = -g_yy_0_yyzzz_z[k] * ab_x + g_yy_0_yyzzz_xz[k];
            }

            /// Set up 309-312 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xyzzzz_x = cbuffer.data(ip_geom_20_off + 309 * ccomps * dcomps);

            auto g_yy_0_xyzzzz_y = cbuffer.data(ip_geom_20_off + 310 * ccomps * dcomps);

            auto g_yy_0_xyzzzz_z = cbuffer.data(ip_geom_20_off + 311 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xyzzzz_x, g_yy_0_xyzzzz_y, g_yy_0_xyzzzz_z, g_yy_0_yzzzz_x, g_yy_0_yzzzz_xx, g_yy_0_yzzzz_xy, g_yy_0_yzzzz_xz, g_yy_0_yzzzz_y, g_yy_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xyzzzz_x[k] = -g_yy_0_yzzzz_x[k] * ab_x + g_yy_0_yzzzz_xx[k];

                g_yy_0_xyzzzz_y[k] = -g_yy_0_yzzzz_y[k] * ab_x + g_yy_0_yzzzz_xy[k];

                g_yy_0_xyzzzz_z[k] = -g_yy_0_yzzzz_z[k] * ab_x + g_yy_0_yzzzz_xz[k];
            }

            /// Set up 312-315 components of targeted buffer : cbuffer.data(

            auto g_yy_0_xzzzzz_x = cbuffer.data(ip_geom_20_off + 312 * ccomps * dcomps);

            auto g_yy_0_xzzzzz_y = cbuffer.data(ip_geom_20_off + 313 * ccomps * dcomps);

            auto g_yy_0_xzzzzz_z = cbuffer.data(ip_geom_20_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_xzzzzz_x, g_yy_0_xzzzzz_y, g_yy_0_xzzzzz_z, g_yy_0_zzzzz_x, g_yy_0_zzzzz_xx, g_yy_0_zzzzz_xy, g_yy_0_zzzzz_xz, g_yy_0_zzzzz_y, g_yy_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_xzzzzz_x[k] = -g_yy_0_zzzzz_x[k] * ab_x + g_yy_0_zzzzz_xx[k];

                g_yy_0_xzzzzz_y[k] = -g_yy_0_zzzzz_y[k] * ab_x + g_yy_0_zzzzz_xy[k];

                g_yy_0_xzzzzz_z[k] = -g_yy_0_zzzzz_z[k] * ab_x + g_yy_0_zzzzz_xz[k];
            }

            /// Set up 315-318 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyyy_x = cbuffer.data(ip_geom_20_off + 315 * ccomps * dcomps);

            auto g_yy_0_yyyyyy_y = cbuffer.data(ip_geom_20_off + 316 * ccomps * dcomps);

            auto g_yy_0_yyyyyy_z = cbuffer.data(ip_geom_20_off + 317 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_x, g_y_0_yyyyy_y, g_y_0_yyyyy_z, g_yy_0_yyyyy_x, g_yy_0_yyyyy_xy, g_yy_0_yyyyy_y, g_yy_0_yyyyy_yy, g_yy_0_yyyyy_yz, g_yy_0_yyyyy_z, g_yy_0_yyyyyy_x, g_yy_0_yyyyyy_y, g_yy_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyyy_x[k] = -2.0 * g_y_0_yyyyy_x[k] - g_yy_0_yyyyy_x[k] * ab_y + g_yy_0_yyyyy_xy[k];

                g_yy_0_yyyyyy_y[k] = -2.0 * g_y_0_yyyyy_y[k] - g_yy_0_yyyyy_y[k] * ab_y + g_yy_0_yyyyy_yy[k];

                g_yy_0_yyyyyy_z[k] = -2.0 * g_y_0_yyyyy_z[k] - g_yy_0_yyyyy_z[k] * ab_y + g_yy_0_yyyyy_yz[k];
            }

            /// Set up 318-321 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyyz_x = cbuffer.data(ip_geom_20_off + 318 * ccomps * dcomps);

            auto g_yy_0_yyyyyz_y = cbuffer.data(ip_geom_20_off + 319 * ccomps * dcomps);

            auto g_yy_0_yyyyyz_z = cbuffer.data(ip_geom_20_off + 320 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyyy_x, g_yy_0_yyyyy_xz, g_yy_0_yyyyy_y, g_yy_0_yyyyy_yz, g_yy_0_yyyyy_z, g_yy_0_yyyyy_zz, g_yy_0_yyyyyz_x, g_yy_0_yyyyyz_y, g_yy_0_yyyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyyz_x[k] = -g_yy_0_yyyyy_x[k] * ab_z + g_yy_0_yyyyy_xz[k];

                g_yy_0_yyyyyz_y[k] = -g_yy_0_yyyyy_y[k] * ab_z + g_yy_0_yyyyy_yz[k];

                g_yy_0_yyyyyz_z[k] = -g_yy_0_yyyyy_z[k] * ab_z + g_yy_0_yyyyy_zz[k];
            }

            /// Set up 321-324 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyyzz_x = cbuffer.data(ip_geom_20_off + 321 * ccomps * dcomps);

            auto g_yy_0_yyyyzz_y = cbuffer.data(ip_geom_20_off + 322 * ccomps * dcomps);

            auto g_yy_0_yyyyzz_z = cbuffer.data(ip_geom_20_off + 323 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyyz_x, g_yy_0_yyyyz_xz, g_yy_0_yyyyz_y, g_yy_0_yyyyz_yz, g_yy_0_yyyyz_z, g_yy_0_yyyyz_zz, g_yy_0_yyyyzz_x, g_yy_0_yyyyzz_y, g_yy_0_yyyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyyzz_x[k] = -g_yy_0_yyyyz_x[k] * ab_z + g_yy_0_yyyyz_xz[k];

                g_yy_0_yyyyzz_y[k] = -g_yy_0_yyyyz_y[k] * ab_z + g_yy_0_yyyyz_yz[k];

                g_yy_0_yyyyzz_z[k] = -g_yy_0_yyyyz_z[k] * ab_z + g_yy_0_yyyyz_zz[k];
            }

            /// Set up 324-327 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyyzzz_x = cbuffer.data(ip_geom_20_off + 324 * ccomps * dcomps);

            auto g_yy_0_yyyzzz_y = cbuffer.data(ip_geom_20_off + 325 * ccomps * dcomps);

            auto g_yy_0_yyyzzz_z = cbuffer.data(ip_geom_20_off + 326 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyyzz_x, g_yy_0_yyyzz_xz, g_yy_0_yyyzz_y, g_yy_0_yyyzz_yz, g_yy_0_yyyzz_z, g_yy_0_yyyzz_zz, g_yy_0_yyyzzz_x, g_yy_0_yyyzzz_y, g_yy_0_yyyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyyzzz_x[k] = -g_yy_0_yyyzz_x[k] * ab_z + g_yy_0_yyyzz_xz[k];

                g_yy_0_yyyzzz_y[k] = -g_yy_0_yyyzz_y[k] * ab_z + g_yy_0_yyyzz_yz[k];

                g_yy_0_yyyzzz_z[k] = -g_yy_0_yyyzz_z[k] * ab_z + g_yy_0_yyyzz_zz[k];
            }

            /// Set up 327-330 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yyzzzz_x = cbuffer.data(ip_geom_20_off + 327 * ccomps * dcomps);

            auto g_yy_0_yyzzzz_y = cbuffer.data(ip_geom_20_off + 328 * ccomps * dcomps);

            auto g_yy_0_yyzzzz_z = cbuffer.data(ip_geom_20_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yyzzz_x, g_yy_0_yyzzz_xz, g_yy_0_yyzzz_y, g_yy_0_yyzzz_yz, g_yy_0_yyzzz_z, g_yy_0_yyzzz_zz, g_yy_0_yyzzzz_x, g_yy_0_yyzzzz_y, g_yy_0_yyzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yyzzzz_x[k] = -g_yy_0_yyzzz_x[k] * ab_z + g_yy_0_yyzzz_xz[k];

                g_yy_0_yyzzzz_y[k] = -g_yy_0_yyzzz_y[k] * ab_z + g_yy_0_yyzzz_yz[k];

                g_yy_0_yyzzzz_z[k] = -g_yy_0_yyzzz_z[k] * ab_z + g_yy_0_yyzzz_zz[k];
            }

            /// Set up 330-333 components of targeted buffer : cbuffer.data(

            auto g_yy_0_yzzzzz_x = cbuffer.data(ip_geom_20_off + 330 * ccomps * dcomps);

            auto g_yy_0_yzzzzz_y = cbuffer.data(ip_geom_20_off + 331 * ccomps * dcomps);

            auto g_yy_0_yzzzzz_z = cbuffer.data(ip_geom_20_off + 332 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_yzzzz_x, g_yy_0_yzzzz_xz, g_yy_0_yzzzz_y, g_yy_0_yzzzz_yz, g_yy_0_yzzzz_z, g_yy_0_yzzzz_zz, g_yy_0_yzzzzz_x, g_yy_0_yzzzzz_y, g_yy_0_yzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_yzzzzz_x[k] = -g_yy_0_yzzzz_x[k] * ab_z + g_yy_0_yzzzz_xz[k];

                g_yy_0_yzzzzz_y[k] = -g_yy_0_yzzzz_y[k] * ab_z + g_yy_0_yzzzz_yz[k];

                g_yy_0_yzzzzz_z[k] = -g_yy_0_yzzzz_z[k] * ab_z + g_yy_0_yzzzz_zz[k];
            }

            /// Set up 333-336 components of targeted buffer : cbuffer.data(

            auto g_yy_0_zzzzzz_x = cbuffer.data(ip_geom_20_off + 333 * ccomps * dcomps);

            auto g_yy_0_zzzzzz_y = cbuffer.data(ip_geom_20_off + 334 * ccomps * dcomps);

            auto g_yy_0_zzzzzz_z = cbuffer.data(ip_geom_20_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_yy_0_zzzzz_x, g_yy_0_zzzzz_xz, g_yy_0_zzzzz_y, g_yy_0_zzzzz_yz, g_yy_0_zzzzz_z, g_yy_0_zzzzz_zz, g_yy_0_zzzzzz_x, g_yy_0_zzzzzz_y, g_yy_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yy_0_zzzzzz_x[k] = -g_yy_0_zzzzz_x[k] * ab_z + g_yy_0_zzzzz_xz[k];

                g_yy_0_zzzzzz_y[k] = -g_yy_0_zzzzz_y[k] * ab_z + g_yy_0_zzzzz_yz[k];

                g_yy_0_zzzzzz_z[k] = -g_yy_0_zzzzz_z[k] * ab_z + g_yy_0_zzzzz_zz[k];
            }

            /// Set up 336-339 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxxx_x = cbuffer.data(ip_geom_20_off + 336 * ccomps * dcomps);

            auto g_yz_0_xxxxxx_y = cbuffer.data(ip_geom_20_off + 337 * ccomps * dcomps);

            auto g_yz_0_xxxxxx_z = cbuffer.data(ip_geom_20_off + 338 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxx_x, g_yz_0_xxxxx_xx, g_yz_0_xxxxx_xy, g_yz_0_xxxxx_xz, g_yz_0_xxxxx_y, g_yz_0_xxxxx_z, g_yz_0_xxxxxx_x, g_yz_0_xxxxxx_y, g_yz_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxxx_x[k] = -g_yz_0_xxxxx_x[k] * ab_x + g_yz_0_xxxxx_xx[k];

                g_yz_0_xxxxxx_y[k] = -g_yz_0_xxxxx_y[k] * ab_x + g_yz_0_xxxxx_xy[k];

                g_yz_0_xxxxxx_z[k] = -g_yz_0_xxxxx_z[k] * ab_x + g_yz_0_xxxxx_xz[k];
            }

            /// Set up 339-342 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxxy_x = cbuffer.data(ip_geom_20_off + 339 * ccomps * dcomps);

            auto g_yz_0_xxxxxy_y = cbuffer.data(ip_geom_20_off + 340 * ccomps * dcomps);

            auto g_yz_0_xxxxxy_z = cbuffer.data(ip_geom_20_off + 341 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxxy_x, g_yz_0_xxxxxy_y, g_yz_0_xxxxxy_z, g_yz_0_xxxxy_x, g_yz_0_xxxxy_xx, g_yz_0_xxxxy_xy, g_yz_0_xxxxy_xz, g_yz_0_xxxxy_y, g_yz_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxxy_x[k] = -g_yz_0_xxxxy_x[k] * ab_x + g_yz_0_xxxxy_xx[k];

                g_yz_0_xxxxxy_y[k] = -g_yz_0_xxxxy_y[k] * ab_x + g_yz_0_xxxxy_xy[k];

                g_yz_0_xxxxxy_z[k] = -g_yz_0_xxxxy_z[k] * ab_x + g_yz_0_xxxxy_xz[k];
            }

            /// Set up 342-345 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxxz_x = cbuffer.data(ip_geom_20_off + 342 * ccomps * dcomps);

            auto g_yz_0_xxxxxz_y = cbuffer.data(ip_geom_20_off + 343 * ccomps * dcomps);

            auto g_yz_0_xxxxxz_z = cbuffer.data(ip_geom_20_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxxz_x, g_yz_0_xxxxxz_y, g_yz_0_xxxxxz_z, g_yz_0_xxxxz_x, g_yz_0_xxxxz_xx, g_yz_0_xxxxz_xy, g_yz_0_xxxxz_xz, g_yz_0_xxxxz_y, g_yz_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxxz_x[k] = -g_yz_0_xxxxz_x[k] * ab_x + g_yz_0_xxxxz_xx[k];

                g_yz_0_xxxxxz_y[k] = -g_yz_0_xxxxz_y[k] * ab_x + g_yz_0_xxxxz_xy[k];

                g_yz_0_xxxxxz_z[k] = -g_yz_0_xxxxz_z[k] * ab_x + g_yz_0_xxxxz_xz[k];
            }

            /// Set up 345-348 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxyy_x = cbuffer.data(ip_geom_20_off + 345 * ccomps * dcomps);

            auto g_yz_0_xxxxyy_y = cbuffer.data(ip_geom_20_off + 346 * ccomps * dcomps);

            auto g_yz_0_xxxxyy_z = cbuffer.data(ip_geom_20_off + 347 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxyy_x, g_yz_0_xxxxyy_y, g_yz_0_xxxxyy_z, g_yz_0_xxxyy_x, g_yz_0_xxxyy_xx, g_yz_0_xxxyy_xy, g_yz_0_xxxyy_xz, g_yz_0_xxxyy_y, g_yz_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxyy_x[k] = -g_yz_0_xxxyy_x[k] * ab_x + g_yz_0_xxxyy_xx[k];

                g_yz_0_xxxxyy_y[k] = -g_yz_0_xxxyy_y[k] * ab_x + g_yz_0_xxxyy_xy[k];

                g_yz_0_xxxxyy_z[k] = -g_yz_0_xxxyy_z[k] * ab_x + g_yz_0_xxxyy_xz[k];
            }

            /// Set up 348-351 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxyz_x = cbuffer.data(ip_geom_20_off + 348 * ccomps * dcomps);

            auto g_yz_0_xxxxyz_y = cbuffer.data(ip_geom_20_off + 349 * ccomps * dcomps);

            auto g_yz_0_xxxxyz_z = cbuffer.data(ip_geom_20_off + 350 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxyz_x, g_yz_0_xxxxyz_y, g_yz_0_xxxxyz_z, g_yz_0_xxxyz_x, g_yz_0_xxxyz_xx, g_yz_0_xxxyz_xy, g_yz_0_xxxyz_xz, g_yz_0_xxxyz_y, g_yz_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxyz_x[k] = -g_yz_0_xxxyz_x[k] * ab_x + g_yz_0_xxxyz_xx[k];

                g_yz_0_xxxxyz_y[k] = -g_yz_0_xxxyz_y[k] * ab_x + g_yz_0_xxxyz_xy[k];

                g_yz_0_xxxxyz_z[k] = -g_yz_0_xxxyz_z[k] * ab_x + g_yz_0_xxxyz_xz[k];
            }

            /// Set up 351-354 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxxzz_x = cbuffer.data(ip_geom_20_off + 351 * ccomps * dcomps);

            auto g_yz_0_xxxxzz_y = cbuffer.data(ip_geom_20_off + 352 * ccomps * dcomps);

            auto g_yz_0_xxxxzz_z = cbuffer.data(ip_geom_20_off + 353 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxxzz_x, g_yz_0_xxxxzz_y, g_yz_0_xxxxzz_z, g_yz_0_xxxzz_x, g_yz_0_xxxzz_xx, g_yz_0_xxxzz_xy, g_yz_0_xxxzz_xz, g_yz_0_xxxzz_y, g_yz_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxxzz_x[k] = -g_yz_0_xxxzz_x[k] * ab_x + g_yz_0_xxxzz_xx[k];

                g_yz_0_xxxxzz_y[k] = -g_yz_0_xxxzz_y[k] * ab_x + g_yz_0_xxxzz_xy[k];

                g_yz_0_xxxxzz_z[k] = -g_yz_0_xxxzz_z[k] * ab_x + g_yz_0_xxxzz_xz[k];
            }

            /// Set up 354-357 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyyy_x = cbuffer.data(ip_geom_20_off + 354 * ccomps * dcomps);

            auto g_yz_0_xxxyyy_y = cbuffer.data(ip_geom_20_off + 355 * ccomps * dcomps);

            auto g_yz_0_xxxyyy_z = cbuffer.data(ip_geom_20_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyyy_x, g_yz_0_xxxyyy_y, g_yz_0_xxxyyy_z, g_yz_0_xxyyy_x, g_yz_0_xxyyy_xx, g_yz_0_xxyyy_xy, g_yz_0_xxyyy_xz, g_yz_0_xxyyy_y, g_yz_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyyy_x[k] = -g_yz_0_xxyyy_x[k] * ab_x + g_yz_0_xxyyy_xx[k];

                g_yz_0_xxxyyy_y[k] = -g_yz_0_xxyyy_y[k] * ab_x + g_yz_0_xxyyy_xy[k];

                g_yz_0_xxxyyy_z[k] = -g_yz_0_xxyyy_z[k] * ab_x + g_yz_0_xxyyy_xz[k];
            }

            /// Set up 357-360 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyyz_x = cbuffer.data(ip_geom_20_off + 357 * ccomps * dcomps);

            auto g_yz_0_xxxyyz_y = cbuffer.data(ip_geom_20_off + 358 * ccomps * dcomps);

            auto g_yz_0_xxxyyz_z = cbuffer.data(ip_geom_20_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyyz_x, g_yz_0_xxxyyz_y, g_yz_0_xxxyyz_z, g_yz_0_xxyyz_x, g_yz_0_xxyyz_xx, g_yz_0_xxyyz_xy, g_yz_0_xxyyz_xz, g_yz_0_xxyyz_y, g_yz_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyyz_x[k] = -g_yz_0_xxyyz_x[k] * ab_x + g_yz_0_xxyyz_xx[k];

                g_yz_0_xxxyyz_y[k] = -g_yz_0_xxyyz_y[k] * ab_x + g_yz_0_xxyyz_xy[k];

                g_yz_0_xxxyyz_z[k] = -g_yz_0_xxyyz_z[k] * ab_x + g_yz_0_xxyyz_xz[k];
            }

            /// Set up 360-363 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxyzz_x = cbuffer.data(ip_geom_20_off + 360 * ccomps * dcomps);

            auto g_yz_0_xxxyzz_y = cbuffer.data(ip_geom_20_off + 361 * ccomps * dcomps);

            auto g_yz_0_xxxyzz_z = cbuffer.data(ip_geom_20_off + 362 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxyzz_x, g_yz_0_xxxyzz_y, g_yz_0_xxxyzz_z, g_yz_0_xxyzz_x, g_yz_0_xxyzz_xx, g_yz_0_xxyzz_xy, g_yz_0_xxyzz_xz, g_yz_0_xxyzz_y, g_yz_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxyzz_x[k] = -g_yz_0_xxyzz_x[k] * ab_x + g_yz_0_xxyzz_xx[k];

                g_yz_0_xxxyzz_y[k] = -g_yz_0_xxyzz_y[k] * ab_x + g_yz_0_xxyzz_xy[k];

                g_yz_0_xxxyzz_z[k] = -g_yz_0_xxyzz_z[k] * ab_x + g_yz_0_xxyzz_xz[k];
            }

            /// Set up 363-366 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxxzzz_x = cbuffer.data(ip_geom_20_off + 363 * ccomps * dcomps);

            auto g_yz_0_xxxzzz_y = cbuffer.data(ip_geom_20_off + 364 * ccomps * dcomps);

            auto g_yz_0_xxxzzz_z = cbuffer.data(ip_geom_20_off + 365 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxxzzz_x, g_yz_0_xxxzzz_y, g_yz_0_xxxzzz_z, g_yz_0_xxzzz_x, g_yz_0_xxzzz_xx, g_yz_0_xxzzz_xy, g_yz_0_xxzzz_xz, g_yz_0_xxzzz_y, g_yz_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxxzzz_x[k] = -g_yz_0_xxzzz_x[k] * ab_x + g_yz_0_xxzzz_xx[k];

                g_yz_0_xxxzzz_y[k] = -g_yz_0_xxzzz_y[k] * ab_x + g_yz_0_xxzzz_xy[k];

                g_yz_0_xxxzzz_z[k] = -g_yz_0_xxzzz_z[k] * ab_x + g_yz_0_xxzzz_xz[k];
            }

            /// Set up 366-369 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyyy_x = cbuffer.data(ip_geom_20_off + 366 * ccomps * dcomps);

            auto g_yz_0_xxyyyy_y = cbuffer.data(ip_geom_20_off + 367 * ccomps * dcomps);

            auto g_yz_0_xxyyyy_z = cbuffer.data(ip_geom_20_off + 368 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyyy_x, g_yz_0_xxyyyy_y, g_yz_0_xxyyyy_z, g_yz_0_xyyyy_x, g_yz_0_xyyyy_xx, g_yz_0_xyyyy_xy, g_yz_0_xyyyy_xz, g_yz_0_xyyyy_y, g_yz_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyyy_x[k] = -g_yz_0_xyyyy_x[k] * ab_x + g_yz_0_xyyyy_xx[k];

                g_yz_0_xxyyyy_y[k] = -g_yz_0_xyyyy_y[k] * ab_x + g_yz_0_xyyyy_xy[k];

                g_yz_0_xxyyyy_z[k] = -g_yz_0_xyyyy_z[k] * ab_x + g_yz_0_xyyyy_xz[k];
            }

            /// Set up 369-372 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyyz_x = cbuffer.data(ip_geom_20_off + 369 * ccomps * dcomps);

            auto g_yz_0_xxyyyz_y = cbuffer.data(ip_geom_20_off + 370 * ccomps * dcomps);

            auto g_yz_0_xxyyyz_z = cbuffer.data(ip_geom_20_off + 371 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyyz_x, g_yz_0_xxyyyz_y, g_yz_0_xxyyyz_z, g_yz_0_xyyyz_x, g_yz_0_xyyyz_xx, g_yz_0_xyyyz_xy, g_yz_0_xyyyz_xz, g_yz_0_xyyyz_y, g_yz_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyyz_x[k] = -g_yz_0_xyyyz_x[k] * ab_x + g_yz_0_xyyyz_xx[k];

                g_yz_0_xxyyyz_y[k] = -g_yz_0_xyyyz_y[k] * ab_x + g_yz_0_xyyyz_xy[k];

                g_yz_0_xxyyyz_z[k] = -g_yz_0_xyyyz_z[k] * ab_x + g_yz_0_xyyyz_xz[k];
            }

            /// Set up 372-375 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyyzz_x = cbuffer.data(ip_geom_20_off + 372 * ccomps * dcomps);

            auto g_yz_0_xxyyzz_y = cbuffer.data(ip_geom_20_off + 373 * ccomps * dcomps);

            auto g_yz_0_xxyyzz_z = cbuffer.data(ip_geom_20_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyyzz_x, g_yz_0_xxyyzz_y, g_yz_0_xxyyzz_z, g_yz_0_xyyzz_x, g_yz_0_xyyzz_xx, g_yz_0_xyyzz_xy, g_yz_0_xyyzz_xz, g_yz_0_xyyzz_y, g_yz_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyyzz_x[k] = -g_yz_0_xyyzz_x[k] * ab_x + g_yz_0_xyyzz_xx[k];

                g_yz_0_xxyyzz_y[k] = -g_yz_0_xyyzz_y[k] * ab_x + g_yz_0_xyyzz_xy[k];

                g_yz_0_xxyyzz_z[k] = -g_yz_0_xyyzz_z[k] * ab_x + g_yz_0_xyyzz_xz[k];
            }

            /// Set up 375-378 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxyzzz_x = cbuffer.data(ip_geom_20_off + 375 * ccomps * dcomps);

            auto g_yz_0_xxyzzz_y = cbuffer.data(ip_geom_20_off + 376 * ccomps * dcomps);

            auto g_yz_0_xxyzzz_z = cbuffer.data(ip_geom_20_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxyzzz_x, g_yz_0_xxyzzz_y, g_yz_0_xxyzzz_z, g_yz_0_xyzzz_x, g_yz_0_xyzzz_xx, g_yz_0_xyzzz_xy, g_yz_0_xyzzz_xz, g_yz_0_xyzzz_y, g_yz_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxyzzz_x[k] = -g_yz_0_xyzzz_x[k] * ab_x + g_yz_0_xyzzz_xx[k];

                g_yz_0_xxyzzz_y[k] = -g_yz_0_xyzzz_y[k] * ab_x + g_yz_0_xyzzz_xy[k];

                g_yz_0_xxyzzz_z[k] = -g_yz_0_xyzzz_z[k] * ab_x + g_yz_0_xyzzz_xz[k];
            }

            /// Set up 378-381 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xxzzzz_x = cbuffer.data(ip_geom_20_off + 378 * ccomps * dcomps);

            auto g_yz_0_xxzzzz_y = cbuffer.data(ip_geom_20_off + 379 * ccomps * dcomps);

            auto g_yz_0_xxzzzz_z = cbuffer.data(ip_geom_20_off + 380 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xxzzzz_x, g_yz_0_xxzzzz_y, g_yz_0_xxzzzz_z, g_yz_0_xzzzz_x, g_yz_0_xzzzz_xx, g_yz_0_xzzzz_xy, g_yz_0_xzzzz_xz, g_yz_0_xzzzz_y, g_yz_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xxzzzz_x[k] = -g_yz_0_xzzzz_x[k] * ab_x + g_yz_0_xzzzz_xx[k];

                g_yz_0_xxzzzz_y[k] = -g_yz_0_xzzzz_y[k] * ab_x + g_yz_0_xzzzz_xy[k];

                g_yz_0_xxzzzz_z[k] = -g_yz_0_xzzzz_z[k] * ab_x + g_yz_0_xzzzz_xz[k];
            }

            /// Set up 381-384 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyyy_x = cbuffer.data(ip_geom_20_off + 381 * ccomps * dcomps);

            auto g_yz_0_xyyyyy_y = cbuffer.data(ip_geom_20_off + 382 * ccomps * dcomps);

            auto g_yz_0_xyyyyy_z = cbuffer.data(ip_geom_20_off + 383 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyyy_x, g_yz_0_xyyyyy_y, g_yz_0_xyyyyy_z, g_yz_0_yyyyy_x, g_yz_0_yyyyy_xx, g_yz_0_yyyyy_xy, g_yz_0_yyyyy_xz, g_yz_0_yyyyy_y, g_yz_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyyy_x[k] = -g_yz_0_yyyyy_x[k] * ab_x + g_yz_0_yyyyy_xx[k];

                g_yz_0_xyyyyy_y[k] = -g_yz_0_yyyyy_y[k] * ab_x + g_yz_0_yyyyy_xy[k];

                g_yz_0_xyyyyy_z[k] = -g_yz_0_yyyyy_z[k] * ab_x + g_yz_0_yyyyy_xz[k];
            }

            /// Set up 384-387 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyyz_x = cbuffer.data(ip_geom_20_off + 384 * ccomps * dcomps);

            auto g_yz_0_xyyyyz_y = cbuffer.data(ip_geom_20_off + 385 * ccomps * dcomps);

            auto g_yz_0_xyyyyz_z = cbuffer.data(ip_geom_20_off + 386 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyyz_x, g_yz_0_xyyyyz_y, g_yz_0_xyyyyz_z, g_yz_0_yyyyz_x, g_yz_0_yyyyz_xx, g_yz_0_yyyyz_xy, g_yz_0_yyyyz_xz, g_yz_0_yyyyz_y, g_yz_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyyz_x[k] = -g_yz_0_yyyyz_x[k] * ab_x + g_yz_0_yyyyz_xx[k];

                g_yz_0_xyyyyz_y[k] = -g_yz_0_yyyyz_y[k] * ab_x + g_yz_0_yyyyz_xy[k];

                g_yz_0_xyyyyz_z[k] = -g_yz_0_yyyyz_z[k] * ab_x + g_yz_0_yyyyz_xz[k];
            }

            /// Set up 387-390 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyyzz_x = cbuffer.data(ip_geom_20_off + 387 * ccomps * dcomps);

            auto g_yz_0_xyyyzz_y = cbuffer.data(ip_geom_20_off + 388 * ccomps * dcomps);

            auto g_yz_0_xyyyzz_z = cbuffer.data(ip_geom_20_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyyzz_x, g_yz_0_xyyyzz_y, g_yz_0_xyyyzz_z, g_yz_0_yyyzz_x, g_yz_0_yyyzz_xx, g_yz_0_yyyzz_xy, g_yz_0_yyyzz_xz, g_yz_0_yyyzz_y, g_yz_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyyzz_x[k] = -g_yz_0_yyyzz_x[k] * ab_x + g_yz_0_yyyzz_xx[k];

                g_yz_0_xyyyzz_y[k] = -g_yz_0_yyyzz_y[k] * ab_x + g_yz_0_yyyzz_xy[k];

                g_yz_0_xyyyzz_z[k] = -g_yz_0_yyyzz_z[k] * ab_x + g_yz_0_yyyzz_xz[k];
            }

            /// Set up 390-393 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyyzzz_x = cbuffer.data(ip_geom_20_off + 390 * ccomps * dcomps);

            auto g_yz_0_xyyzzz_y = cbuffer.data(ip_geom_20_off + 391 * ccomps * dcomps);

            auto g_yz_0_xyyzzz_z = cbuffer.data(ip_geom_20_off + 392 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyyzzz_x, g_yz_0_xyyzzz_y, g_yz_0_xyyzzz_z, g_yz_0_yyzzz_x, g_yz_0_yyzzz_xx, g_yz_0_yyzzz_xy, g_yz_0_yyzzz_xz, g_yz_0_yyzzz_y, g_yz_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyyzzz_x[k] = -g_yz_0_yyzzz_x[k] * ab_x + g_yz_0_yyzzz_xx[k];

                g_yz_0_xyyzzz_y[k] = -g_yz_0_yyzzz_y[k] * ab_x + g_yz_0_yyzzz_xy[k];

                g_yz_0_xyyzzz_z[k] = -g_yz_0_yyzzz_z[k] * ab_x + g_yz_0_yyzzz_xz[k];
            }

            /// Set up 393-396 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xyzzzz_x = cbuffer.data(ip_geom_20_off + 393 * ccomps * dcomps);

            auto g_yz_0_xyzzzz_y = cbuffer.data(ip_geom_20_off + 394 * ccomps * dcomps);

            auto g_yz_0_xyzzzz_z = cbuffer.data(ip_geom_20_off + 395 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xyzzzz_x, g_yz_0_xyzzzz_y, g_yz_0_xyzzzz_z, g_yz_0_yzzzz_x, g_yz_0_yzzzz_xx, g_yz_0_yzzzz_xy, g_yz_0_yzzzz_xz, g_yz_0_yzzzz_y, g_yz_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xyzzzz_x[k] = -g_yz_0_yzzzz_x[k] * ab_x + g_yz_0_yzzzz_xx[k];

                g_yz_0_xyzzzz_y[k] = -g_yz_0_yzzzz_y[k] * ab_x + g_yz_0_yzzzz_xy[k];

                g_yz_0_xyzzzz_z[k] = -g_yz_0_yzzzz_z[k] * ab_x + g_yz_0_yzzzz_xz[k];
            }

            /// Set up 396-399 components of targeted buffer : cbuffer.data(

            auto g_yz_0_xzzzzz_x = cbuffer.data(ip_geom_20_off + 396 * ccomps * dcomps);

            auto g_yz_0_xzzzzz_y = cbuffer.data(ip_geom_20_off + 397 * ccomps * dcomps);

            auto g_yz_0_xzzzzz_z = cbuffer.data(ip_geom_20_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_xzzzzz_x, g_yz_0_xzzzzz_y, g_yz_0_xzzzzz_z, g_yz_0_zzzzz_x, g_yz_0_zzzzz_xx, g_yz_0_zzzzz_xy, g_yz_0_zzzzz_xz, g_yz_0_zzzzz_y, g_yz_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_xzzzzz_x[k] = -g_yz_0_zzzzz_x[k] * ab_x + g_yz_0_zzzzz_xx[k];

                g_yz_0_xzzzzz_y[k] = -g_yz_0_zzzzz_y[k] * ab_x + g_yz_0_zzzzz_xy[k];

                g_yz_0_xzzzzz_z[k] = -g_yz_0_zzzzz_z[k] * ab_x + g_yz_0_zzzzz_xz[k];
            }

            /// Set up 399-402 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyyy_x = cbuffer.data(ip_geom_20_off + 399 * ccomps * dcomps);

            auto g_yz_0_yyyyyy_y = cbuffer.data(ip_geom_20_off + 400 * ccomps * dcomps);

            auto g_yz_0_yyyyyy_z = cbuffer.data(ip_geom_20_off + 401 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyy_x, g_yz_0_yyyyy_xy, g_yz_0_yyyyy_y, g_yz_0_yyyyy_yy, g_yz_0_yyyyy_yz, g_yz_0_yyyyy_z, g_yz_0_yyyyyy_x, g_yz_0_yyyyyy_y, g_yz_0_yyyyyy_z, g_z_0_yyyyy_x, g_z_0_yyyyy_y, g_z_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyyy_x[k] = -g_z_0_yyyyy_x[k] - g_yz_0_yyyyy_x[k] * ab_y + g_yz_0_yyyyy_xy[k];

                g_yz_0_yyyyyy_y[k] = -g_z_0_yyyyy_y[k] - g_yz_0_yyyyy_y[k] * ab_y + g_yz_0_yyyyy_yy[k];

                g_yz_0_yyyyyy_z[k] = -g_z_0_yyyyy_z[k] - g_yz_0_yyyyy_z[k] * ab_y + g_yz_0_yyyyy_yz[k];
            }

            /// Set up 402-405 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyyz_x = cbuffer.data(ip_geom_20_off + 402 * ccomps * dcomps);

            auto g_yz_0_yyyyyz_y = cbuffer.data(ip_geom_20_off + 403 * ccomps * dcomps);

            auto g_yz_0_yyyyyz_z = cbuffer.data(ip_geom_20_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyyz_x, g_yz_0_yyyyyz_y, g_yz_0_yyyyyz_z, g_yz_0_yyyyz_x, g_yz_0_yyyyz_xy, g_yz_0_yyyyz_y, g_yz_0_yyyyz_yy, g_yz_0_yyyyz_yz, g_yz_0_yyyyz_z, g_z_0_yyyyz_x, g_z_0_yyyyz_y, g_z_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyyz_x[k] = -g_z_0_yyyyz_x[k] - g_yz_0_yyyyz_x[k] * ab_y + g_yz_0_yyyyz_xy[k];

                g_yz_0_yyyyyz_y[k] = -g_z_0_yyyyz_y[k] - g_yz_0_yyyyz_y[k] * ab_y + g_yz_0_yyyyz_yy[k];

                g_yz_0_yyyyyz_z[k] = -g_z_0_yyyyz_z[k] - g_yz_0_yyyyz_z[k] * ab_y + g_yz_0_yyyyz_yz[k];
            }

            /// Set up 405-408 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyyzz_x = cbuffer.data(ip_geom_20_off + 405 * ccomps * dcomps);

            auto g_yz_0_yyyyzz_y = cbuffer.data(ip_geom_20_off + 406 * ccomps * dcomps);

            auto g_yz_0_yyyyzz_z = cbuffer.data(ip_geom_20_off + 407 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyyzz_x, g_yz_0_yyyyzz_y, g_yz_0_yyyyzz_z, g_yz_0_yyyzz_x, g_yz_0_yyyzz_xy, g_yz_0_yyyzz_y, g_yz_0_yyyzz_yy, g_yz_0_yyyzz_yz, g_yz_0_yyyzz_z, g_z_0_yyyzz_x, g_z_0_yyyzz_y, g_z_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyyzz_x[k] = -g_z_0_yyyzz_x[k] - g_yz_0_yyyzz_x[k] * ab_y + g_yz_0_yyyzz_xy[k];

                g_yz_0_yyyyzz_y[k] = -g_z_0_yyyzz_y[k] - g_yz_0_yyyzz_y[k] * ab_y + g_yz_0_yyyzz_yy[k];

                g_yz_0_yyyyzz_z[k] = -g_z_0_yyyzz_z[k] - g_yz_0_yyyzz_z[k] * ab_y + g_yz_0_yyyzz_yz[k];
            }

            /// Set up 408-411 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyyzzz_x = cbuffer.data(ip_geom_20_off + 408 * ccomps * dcomps);

            auto g_yz_0_yyyzzz_y = cbuffer.data(ip_geom_20_off + 409 * ccomps * dcomps);

            auto g_yz_0_yyyzzz_z = cbuffer.data(ip_geom_20_off + 410 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyyzzz_x, g_yz_0_yyyzzz_y, g_yz_0_yyyzzz_z, g_yz_0_yyzzz_x, g_yz_0_yyzzz_xy, g_yz_0_yyzzz_y, g_yz_0_yyzzz_yy, g_yz_0_yyzzz_yz, g_yz_0_yyzzz_z, g_z_0_yyzzz_x, g_z_0_yyzzz_y, g_z_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyyzzz_x[k] = -g_z_0_yyzzz_x[k] - g_yz_0_yyzzz_x[k] * ab_y + g_yz_0_yyzzz_xy[k];

                g_yz_0_yyyzzz_y[k] = -g_z_0_yyzzz_y[k] - g_yz_0_yyzzz_y[k] * ab_y + g_yz_0_yyzzz_yy[k];

                g_yz_0_yyyzzz_z[k] = -g_z_0_yyzzz_z[k] - g_yz_0_yyzzz_z[k] * ab_y + g_yz_0_yyzzz_yz[k];
            }

            /// Set up 411-414 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yyzzzz_x = cbuffer.data(ip_geom_20_off + 411 * ccomps * dcomps);

            auto g_yz_0_yyzzzz_y = cbuffer.data(ip_geom_20_off + 412 * ccomps * dcomps);

            auto g_yz_0_yyzzzz_z = cbuffer.data(ip_geom_20_off + 413 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yyzzzz_x, g_yz_0_yyzzzz_y, g_yz_0_yyzzzz_z, g_yz_0_yzzzz_x, g_yz_0_yzzzz_xy, g_yz_0_yzzzz_y, g_yz_0_yzzzz_yy, g_yz_0_yzzzz_yz, g_yz_0_yzzzz_z, g_z_0_yzzzz_x, g_z_0_yzzzz_y, g_z_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yyzzzz_x[k] = -g_z_0_yzzzz_x[k] - g_yz_0_yzzzz_x[k] * ab_y + g_yz_0_yzzzz_xy[k];

                g_yz_0_yyzzzz_y[k] = -g_z_0_yzzzz_y[k] - g_yz_0_yzzzz_y[k] * ab_y + g_yz_0_yzzzz_yy[k];

                g_yz_0_yyzzzz_z[k] = -g_z_0_yzzzz_z[k] - g_yz_0_yzzzz_z[k] * ab_y + g_yz_0_yzzzz_yz[k];
            }

            /// Set up 414-417 components of targeted buffer : cbuffer.data(

            auto g_yz_0_yzzzzz_x = cbuffer.data(ip_geom_20_off + 414 * ccomps * dcomps);

            auto g_yz_0_yzzzzz_y = cbuffer.data(ip_geom_20_off + 415 * ccomps * dcomps);

            auto g_yz_0_yzzzzz_z = cbuffer.data(ip_geom_20_off + 416 * ccomps * dcomps);

            #pragma omp simd aligned(g_yz_0_yzzzzz_x, g_yz_0_yzzzzz_y, g_yz_0_yzzzzz_z, g_yz_0_zzzzz_x, g_yz_0_zzzzz_xy, g_yz_0_zzzzz_y, g_yz_0_zzzzz_yy, g_yz_0_zzzzz_yz, g_yz_0_zzzzz_z, g_z_0_zzzzz_x, g_z_0_zzzzz_y, g_z_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_yzzzzz_x[k] = -g_z_0_zzzzz_x[k] - g_yz_0_zzzzz_x[k] * ab_y + g_yz_0_zzzzz_xy[k];

                g_yz_0_yzzzzz_y[k] = -g_z_0_zzzzz_y[k] - g_yz_0_zzzzz_y[k] * ab_y + g_yz_0_zzzzz_yy[k];

                g_yz_0_yzzzzz_z[k] = -g_z_0_zzzzz_z[k] - g_yz_0_zzzzz_z[k] * ab_y + g_yz_0_zzzzz_yz[k];
            }

            /// Set up 417-420 components of targeted buffer : cbuffer.data(

            auto g_yz_0_zzzzzz_x = cbuffer.data(ip_geom_20_off + 417 * ccomps * dcomps);

            auto g_yz_0_zzzzzz_y = cbuffer.data(ip_geom_20_off + 418 * ccomps * dcomps);

            auto g_yz_0_zzzzzz_z = cbuffer.data(ip_geom_20_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzzz_x, g_y_0_zzzzz_y, g_y_0_zzzzz_z, g_yz_0_zzzzz_x, g_yz_0_zzzzz_xz, g_yz_0_zzzzz_y, g_yz_0_zzzzz_yz, g_yz_0_zzzzz_z, g_yz_0_zzzzz_zz, g_yz_0_zzzzzz_x, g_yz_0_zzzzzz_y, g_yz_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_yz_0_zzzzzz_x[k] = -g_y_0_zzzzz_x[k] - g_yz_0_zzzzz_x[k] * ab_z + g_yz_0_zzzzz_xz[k];

                g_yz_0_zzzzzz_y[k] = -g_y_0_zzzzz_y[k] - g_yz_0_zzzzz_y[k] * ab_z + g_yz_0_zzzzz_yz[k];

                g_yz_0_zzzzzz_z[k] = -g_y_0_zzzzz_z[k] - g_yz_0_zzzzz_z[k] * ab_z + g_yz_0_zzzzz_zz[k];
            }

            /// Set up 420-423 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxxx_x = cbuffer.data(ip_geom_20_off + 420 * ccomps * dcomps);

            auto g_zz_0_xxxxxx_y = cbuffer.data(ip_geom_20_off + 421 * ccomps * dcomps);

            auto g_zz_0_xxxxxx_z = cbuffer.data(ip_geom_20_off + 422 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxx_x, g_zz_0_xxxxx_xx, g_zz_0_xxxxx_xy, g_zz_0_xxxxx_xz, g_zz_0_xxxxx_y, g_zz_0_xxxxx_z, g_zz_0_xxxxxx_x, g_zz_0_xxxxxx_y, g_zz_0_xxxxxx_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxxx_x[k] = -g_zz_0_xxxxx_x[k] * ab_x + g_zz_0_xxxxx_xx[k];

                g_zz_0_xxxxxx_y[k] = -g_zz_0_xxxxx_y[k] * ab_x + g_zz_0_xxxxx_xy[k];

                g_zz_0_xxxxxx_z[k] = -g_zz_0_xxxxx_z[k] * ab_x + g_zz_0_xxxxx_xz[k];
            }

            /// Set up 423-426 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxxy_x = cbuffer.data(ip_geom_20_off + 423 * ccomps * dcomps);

            auto g_zz_0_xxxxxy_y = cbuffer.data(ip_geom_20_off + 424 * ccomps * dcomps);

            auto g_zz_0_xxxxxy_z = cbuffer.data(ip_geom_20_off + 425 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxxy_x, g_zz_0_xxxxxy_y, g_zz_0_xxxxxy_z, g_zz_0_xxxxy_x, g_zz_0_xxxxy_xx, g_zz_0_xxxxy_xy, g_zz_0_xxxxy_xz, g_zz_0_xxxxy_y, g_zz_0_xxxxy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxxy_x[k] = -g_zz_0_xxxxy_x[k] * ab_x + g_zz_0_xxxxy_xx[k];

                g_zz_0_xxxxxy_y[k] = -g_zz_0_xxxxy_y[k] * ab_x + g_zz_0_xxxxy_xy[k];

                g_zz_0_xxxxxy_z[k] = -g_zz_0_xxxxy_z[k] * ab_x + g_zz_0_xxxxy_xz[k];
            }

            /// Set up 426-429 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxxz_x = cbuffer.data(ip_geom_20_off + 426 * ccomps * dcomps);

            auto g_zz_0_xxxxxz_y = cbuffer.data(ip_geom_20_off + 427 * ccomps * dcomps);

            auto g_zz_0_xxxxxz_z = cbuffer.data(ip_geom_20_off + 428 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxxz_x, g_zz_0_xxxxxz_y, g_zz_0_xxxxxz_z, g_zz_0_xxxxz_x, g_zz_0_xxxxz_xx, g_zz_0_xxxxz_xy, g_zz_0_xxxxz_xz, g_zz_0_xxxxz_y, g_zz_0_xxxxz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxxz_x[k] = -g_zz_0_xxxxz_x[k] * ab_x + g_zz_0_xxxxz_xx[k];

                g_zz_0_xxxxxz_y[k] = -g_zz_0_xxxxz_y[k] * ab_x + g_zz_0_xxxxz_xy[k];

                g_zz_0_xxxxxz_z[k] = -g_zz_0_xxxxz_z[k] * ab_x + g_zz_0_xxxxz_xz[k];
            }

            /// Set up 429-432 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxyy_x = cbuffer.data(ip_geom_20_off + 429 * ccomps * dcomps);

            auto g_zz_0_xxxxyy_y = cbuffer.data(ip_geom_20_off + 430 * ccomps * dcomps);

            auto g_zz_0_xxxxyy_z = cbuffer.data(ip_geom_20_off + 431 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxyy_x, g_zz_0_xxxxyy_y, g_zz_0_xxxxyy_z, g_zz_0_xxxyy_x, g_zz_0_xxxyy_xx, g_zz_0_xxxyy_xy, g_zz_0_xxxyy_xz, g_zz_0_xxxyy_y, g_zz_0_xxxyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxyy_x[k] = -g_zz_0_xxxyy_x[k] * ab_x + g_zz_0_xxxyy_xx[k];

                g_zz_0_xxxxyy_y[k] = -g_zz_0_xxxyy_y[k] * ab_x + g_zz_0_xxxyy_xy[k];

                g_zz_0_xxxxyy_z[k] = -g_zz_0_xxxyy_z[k] * ab_x + g_zz_0_xxxyy_xz[k];
            }

            /// Set up 432-435 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxyz_x = cbuffer.data(ip_geom_20_off + 432 * ccomps * dcomps);

            auto g_zz_0_xxxxyz_y = cbuffer.data(ip_geom_20_off + 433 * ccomps * dcomps);

            auto g_zz_0_xxxxyz_z = cbuffer.data(ip_geom_20_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxyz_x, g_zz_0_xxxxyz_y, g_zz_0_xxxxyz_z, g_zz_0_xxxyz_x, g_zz_0_xxxyz_xx, g_zz_0_xxxyz_xy, g_zz_0_xxxyz_xz, g_zz_0_xxxyz_y, g_zz_0_xxxyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxyz_x[k] = -g_zz_0_xxxyz_x[k] * ab_x + g_zz_0_xxxyz_xx[k];

                g_zz_0_xxxxyz_y[k] = -g_zz_0_xxxyz_y[k] * ab_x + g_zz_0_xxxyz_xy[k];

                g_zz_0_xxxxyz_z[k] = -g_zz_0_xxxyz_z[k] * ab_x + g_zz_0_xxxyz_xz[k];
            }

            /// Set up 435-438 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxxzz_x = cbuffer.data(ip_geom_20_off + 435 * ccomps * dcomps);

            auto g_zz_0_xxxxzz_y = cbuffer.data(ip_geom_20_off + 436 * ccomps * dcomps);

            auto g_zz_0_xxxxzz_z = cbuffer.data(ip_geom_20_off + 437 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxxzz_x, g_zz_0_xxxxzz_y, g_zz_0_xxxxzz_z, g_zz_0_xxxzz_x, g_zz_0_xxxzz_xx, g_zz_0_xxxzz_xy, g_zz_0_xxxzz_xz, g_zz_0_xxxzz_y, g_zz_0_xxxzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxxzz_x[k] = -g_zz_0_xxxzz_x[k] * ab_x + g_zz_0_xxxzz_xx[k];

                g_zz_0_xxxxzz_y[k] = -g_zz_0_xxxzz_y[k] * ab_x + g_zz_0_xxxzz_xy[k];

                g_zz_0_xxxxzz_z[k] = -g_zz_0_xxxzz_z[k] * ab_x + g_zz_0_xxxzz_xz[k];
            }

            /// Set up 438-441 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyyy_x = cbuffer.data(ip_geom_20_off + 438 * ccomps * dcomps);

            auto g_zz_0_xxxyyy_y = cbuffer.data(ip_geom_20_off + 439 * ccomps * dcomps);

            auto g_zz_0_xxxyyy_z = cbuffer.data(ip_geom_20_off + 440 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyyy_x, g_zz_0_xxxyyy_y, g_zz_0_xxxyyy_z, g_zz_0_xxyyy_x, g_zz_0_xxyyy_xx, g_zz_0_xxyyy_xy, g_zz_0_xxyyy_xz, g_zz_0_xxyyy_y, g_zz_0_xxyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyyy_x[k] = -g_zz_0_xxyyy_x[k] * ab_x + g_zz_0_xxyyy_xx[k];

                g_zz_0_xxxyyy_y[k] = -g_zz_0_xxyyy_y[k] * ab_x + g_zz_0_xxyyy_xy[k];

                g_zz_0_xxxyyy_z[k] = -g_zz_0_xxyyy_z[k] * ab_x + g_zz_0_xxyyy_xz[k];
            }

            /// Set up 441-444 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyyz_x = cbuffer.data(ip_geom_20_off + 441 * ccomps * dcomps);

            auto g_zz_0_xxxyyz_y = cbuffer.data(ip_geom_20_off + 442 * ccomps * dcomps);

            auto g_zz_0_xxxyyz_z = cbuffer.data(ip_geom_20_off + 443 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyyz_x, g_zz_0_xxxyyz_y, g_zz_0_xxxyyz_z, g_zz_0_xxyyz_x, g_zz_0_xxyyz_xx, g_zz_0_xxyyz_xy, g_zz_0_xxyyz_xz, g_zz_0_xxyyz_y, g_zz_0_xxyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyyz_x[k] = -g_zz_0_xxyyz_x[k] * ab_x + g_zz_0_xxyyz_xx[k];

                g_zz_0_xxxyyz_y[k] = -g_zz_0_xxyyz_y[k] * ab_x + g_zz_0_xxyyz_xy[k];

                g_zz_0_xxxyyz_z[k] = -g_zz_0_xxyyz_z[k] * ab_x + g_zz_0_xxyyz_xz[k];
            }

            /// Set up 444-447 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxyzz_x = cbuffer.data(ip_geom_20_off + 444 * ccomps * dcomps);

            auto g_zz_0_xxxyzz_y = cbuffer.data(ip_geom_20_off + 445 * ccomps * dcomps);

            auto g_zz_0_xxxyzz_z = cbuffer.data(ip_geom_20_off + 446 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxyzz_x, g_zz_0_xxxyzz_y, g_zz_0_xxxyzz_z, g_zz_0_xxyzz_x, g_zz_0_xxyzz_xx, g_zz_0_xxyzz_xy, g_zz_0_xxyzz_xz, g_zz_0_xxyzz_y, g_zz_0_xxyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxyzz_x[k] = -g_zz_0_xxyzz_x[k] * ab_x + g_zz_0_xxyzz_xx[k];

                g_zz_0_xxxyzz_y[k] = -g_zz_0_xxyzz_y[k] * ab_x + g_zz_0_xxyzz_xy[k];

                g_zz_0_xxxyzz_z[k] = -g_zz_0_xxyzz_z[k] * ab_x + g_zz_0_xxyzz_xz[k];
            }

            /// Set up 447-450 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxxzzz_x = cbuffer.data(ip_geom_20_off + 447 * ccomps * dcomps);

            auto g_zz_0_xxxzzz_y = cbuffer.data(ip_geom_20_off + 448 * ccomps * dcomps);

            auto g_zz_0_xxxzzz_z = cbuffer.data(ip_geom_20_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxxzzz_x, g_zz_0_xxxzzz_y, g_zz_0_xxxzzz_z, g_zz_0_xxzzz_x, g_zz_0_xxzzz_xx, g_zz_0_xxzzz_xy, g_zz_0_xxzzz_xz, g_zz_0_xxzzz_y, g_zz_0_xxzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxxzzz_x[k] = -g_zz_0_xxzzz_x[k] * ab_x + g_zz_0_xxzzz_xx[k];

                g_zz_0_xxxzzz_y[k] = -g_zz_0_xxzzz_y[k] * ab_x + g_zz_0_xxzzz_xy[k];

                g_zz_0_xxxzzz_z[k] = -g_zz_0_xxzzz_z[k] * ab_x + g_zz_0_xxzzz_xz[k];
            }

            /// Set up 450-453 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyyy_x = cbuffer.data(ip_geom_20_off + 450 * ccomps * dcomps);

            auto g_zz_0_xxyyyy_y = cbuffer.data(ip_geom_20_off + 451 * ccomps * dcomps);

            auto g_zz_0_xxyyyy_z = cbuffer.data(ip_geom_20_off + 452 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyyy_x, g_zz_0_xxyyyy_y, g_zz_0_xxyyyy_z, g_zz_0_xyyyy_x, g_zz_0_xyyyy_xx, g_zz_0_xyyyy_xy, g_zz_0_xyyyy_xz, g_zz_0_xyyyy_y, g_zz_0_xyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyyy_x[k] = -g_zz_0_xyyyy_x[k] * ab_x + g_zz_0_xyyyy_xx[k];

                g_zz_0_xxyyyy_y[k] = -g_zz_0_xyyyy_y[k] * ab_x + g_zz_0_xyyyy_xy[k];

                g_zz_0_xxyyyy_z[k] = -g_zz_0_xyyyy_z[k] * ab_x + g_zz_0_xyyyy_xz[k];
            }

            /// Set up 453-456 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyyz_x = cbuffer.data(ip_geom_20_off + 453 * ccomps * dcomps);

            auto g_zz_0_xxyyyz_y = cbuffer.data(ip_geom_20_off + 454 * ccomps * dcomps);

            auto g_zz_0_xxyyyz_z = cbuffer.data(ip_geom_20_off + 455 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyyz_x, g_zz_0_xxyyyz_y, g_zz_0_xxyyyz_z, g_zz_0_xyyyz_x, g_zz_0_xyyyz_xx, g_zz_0_xyyyz_xy, g_zz_0_xyyyz_xz, g_zz_0_xyyyz_y, g_zz_0_xyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyyz_x[k] = -g_zz_0_xyyyz_x[k] * ab_x + g_zz_0_xyyyz_xx[k];

                g_zz_0_xxyyyz_y[k] = -g_zz_0_xyyyz_y[k] * ab_x + g_zz_0_xyyyz_xy[k];

                g_zz_0_xxyyyz_z[k] = -g_zz_0_xyyyz_z[k] * ab_x + g_zz_0_xyyyz_xz[k];
            }

            /// Set up 456-459 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyyzz_x = cbuffer.data(ip_geom_20_off + 456 * ccomps * dcomps);

            auto g_zz_0_xxyyzz_y = cbuffer.data(ip_geom_20_off + 457 * ccomps * dcomps);

            auto g_zz_0_xxyyzz_z = cbuffer.data(ip_geom_20_off + 458 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyyzz_x, g_zz_0_xxyyzz_y, g_zz_0_xxyyzz_z, g_zz_0_xyyzz_x, g_zz_0_xyyzz_xx, g_zz_0_xyyzz_xy, g_zz_0_xyyzz_xz, g_zz_0_xyyzz_y, g_zz_0_xyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyyzz_x[k] = -g_zz_0_xyyzz_x[k] * ab_x + g_zz_0_xyyzz_xx[k];

                g_zz_0_xxyyzz_y[k] = -g_zz_0_xyyzz_y[k] * ab_x + g_zz_0_xyyzz_xy[k];

                g_zz_0_xxyyzz_z[k] = -g_zz_0_xyyzz_z[k] * ab_x + g_zz_0_xyyzz_xz[k];
            }

            /// Set up 459-462 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxyzzz_x = cbuffer.data(ip_geom_20_off + 459 * ccomps * dcomps);

            auto g_zz_0_xxyzzz_y = cbuffer.data(ip_geom_20_off + 460 * ccomps * dcomps);

            auto g_zz_0_xxyzzz_z = cbuffer.data(ip_geom_20_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxyzzz_x, g_zz_0_xxyzzz_y, g_zz_0_xxyzzz_z, g_zz_0_xyzzz_x, g_zz_0_xyzzz_xx, g_zz_0_xyzzz_xy, g_zz_0_xyzzz_xz, g_zz_0_xyzzz_y, g_zz_0_xyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxyzzz_x[k] = -g_zz_0_xyzzz_x[k] * ab_x + g_zz_0_xyzzz_xx[k];

                g_zz_0_xxyzzz_y[k] = -g_zz_0_xyzzz_y[k] * ab_x + g_zz_0_xyzzz_xy[k];

                g_zz_0_xxyzzz_z[k] = -g_zz_0_xyzzz_z[k] * ab_x + g_zz_0_xyzzz_xz[k];
            }

            /// Set up 462-465 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xxzzzz_x = cbuffer.data(ip_geom_20_off + 462 * ccomps * dcomps);

            auto g_zz_0_xxzzzz_y = cbuffer.data(ip_geom_20_off + 463 * ccomps * dcomps);

            auto g_zz_0_xxzzzz_z = cbuffer.data(ip_geom_20_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xxzzzz_x, g_zz_0_xxzzzz_y, g_zz_0_xxzzzz_z, g_zz_0_xzzzz_x, g_zz_0_xzzzz_xx, g_zz_0_xzzzz_xy, g_zz_0_xzzzz_xz, g_zz_0_xzzzz_y, g_zz_0_xzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xxzzzz_x[k] = -g_zz_0_xzzzz_x[k] * ab_x + g_zz_0_xzzzz_xx[k];

                g_zz_0_xxzzzz_y[k] = -g_zz_0_xzzzz_y[k] * ab_x + g_zz_0_xzzzz_xy[k];

                g_zz_0_xxzzzz_z[k] = -g_zz_0_xzzzz_z[k] * ab_x + g_zz_0_xzzzz_xz[k];
            }

            /// Set up 465-468 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyyy_x = cbuffer.data(ip_geom_20_off + 465 * ccomps * dcomps);

            auto g_zz_0_xyyyyy_y = cbuffer.data(ip_geom_20_off + 466 * ccomps * dcomps);

            auto g_zz_0_xyyyyy_z = cbuffer.data(ip_geom_20_off + 467 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyyy_x, g_zz_0_xyyyyy_y, g_zz_0_xyyyyy_z, g_zz_0_yyyyy_x, g_zz_0_yyyyy_xx, g_zz_0_yyyyy_xy, g_zz_0_yyyyy_xz, g_zz_0_yyyyy_y, g_zz_0_yyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyyy_x[k] = -g_zz_0_yyyyy_x[k] * ab_x + g_zz_0_yyyyy_xx[k];

                g_zz_0_xyyyyy_y[k] = -g_zz_0_yyyyy_y[k] * ab_x + g_zz_0_yyyyy_xy[k];

                g_zz_0_xyyyyy_z[k] = -g_zz_0_yyyyy_z[k] * ab_x + g_zz_0_yyyyy_xz[k];
            }

            /// Set up 468-471 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyyz_x = cbuffer.data(ip_geom_20_off + 468 * ccomps * dcomps);

            auto g_zz_0_xyyyyz_y = cbuffer.data(ip_geom_20_off + 469 * ccomps * dcomps);

            auto g_zz_0_xyyyyz_z = cbuffer.data(ip_geom_20_off + 470 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyyz_x, g_zz_0_xyyyyz_y, g_zz_0_xyyyyz_z, g_zz_0_yyyyz_x, g_zz_0_yyyyz_xx, g_zz_0_yyyyz_xy, g_zz_0_yyyyz_xz, g_zz_0_yyyyz_y, g_zz_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyyz_x[k] = -g_zz_0_yyyyz_x[k] * ab_x + g_zz_0_yyyyz_xx[k];

                g_zz_0_xyyyyz_y[k] = -g_zz_0_yyyyz_y[k] * ab_x + g_zz_0_yyyyz_xy[k];

                g_zz_0_xyyyyz_z[k] = -g_zz_0_yyyyz_z[k] * ab_x + g_zz_0_yyyyz_xz[k];
            }

            /// Set up 471-474 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyyzz_x = cbuffer.data(ip_geom_20_off + 471 * ccomps * dcomps);

            auto g_zz_0_xyyyzz_y = cbuffer.data(ip_geom_20_off + 472 * ccomps * dcomps);

            auto g_zz_0_xyyyzz_z = cbuffer.data(ip_geom_20_off + 473 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyyzz_x, g_zz_0_xyyyzz_y, g_zz_0_xyyyzz_z, g_zz_0_yyyzz_x, g_zz_0_yyyzz_xx, g_zz_0_yyyzz_xy, g_zz_0_yyyzz_xz, g_zz_0_yyyzz_y, g_zz_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyyzz_x[k] = -g_zz_0_yyyzz_x[k] * ab_x + g_zz_0_yyyzz_xx[k];

                g_zz_0_xyyyzz_y[k] = -g_zz_0_yyyzz_y[k] * ab_x + g_zz_0_yyyzz_xy[k];

                g_zz_0_xyyyzz_z[k] = -g_zz_0_yyyzz_z[k] * ab_x + g_zz_0_yyyzz_xz[k];
            }

            /// Set up 474-477 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyyzzz_x = cbuffer.data(ip_geom_20_off + 474 * ccomps * dcomps);

            auto g_zz_0_xyyzzz_y = cbuffer.data(ip_geom_20_off + 475 * ccomps * dcomps);

            auto g_zz_0_xyyzzz_z = cbuffer.data(ip_geom_20_off + 476 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyyzzz_x, g_zz_0_xyyzzz_y, g_zz_0_xyyzzz_z, g_zz_0_yyzzz_x, g_zz_0_yyzzz_xx, g_zz_0_yyzzz_xy, g_zz_0_yyzzz_xz, g_zz_0_yyzzz_y, g_zz_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyyzzz_x[k] = -g_zz_0_yyzzz_x[k] * ab_x + g_zz_0_yyzzz_xx[k];

                g_zz_0_xyyzzz_y[k] = -g_zz_0_yyzzz_y[k] * ab_x + g_zz_0_yyzzz_xy[k];

                g_zz_0_xyyzzz_z[k] = -g_zz_0_yyzzz_z[k] * ab_x + g_zz_0_yyzzz_xz[k];
            }

            /// Set up 477-480 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xyzzzz_x = cbuffer.data(ip_geom_20_off + 477 * ccomps * dcomps);

            auto g_zz_0_xyzzzz_y = cbuffer.data(ip_geom_20_off + 478 * ccomps * dcomps);

            auto g_zz_0_xyzzzz_z = cbuffer.data(ip_geom_20_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xyzzzz_x, g_zz_0_xyzzzz_y, g_zz_0_xyzzzz_z, g_zz_0_yzzzz_x, g_zz_0_yzzzz_xx, g_zz_0_yzzzz_xy, g_zz_0_yzzzz_xz, g_zz_0_yzzzz_y, g_zz_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xyzzzz_x[k] = -g_zz_0_yzzzz_x[k] * ab_x + g_zz_0_yzzzz_xx[k];

                g_zz_0_xyzzzz_y[k] = -g_zz_0_yzzzz_y[k] * ab_x + g_zz_0_yzzzz_xy[k];

                g_zz_0_xyzzzz_z[k] = -g_zz_0_yzzzz_z[k] * ab_x + g_zz_0_yzzzz_xz[k];
            }

            /// Set up 480-483 components of targeted buffer : cbuffer.data(

            auto g_zz_0_xzzzzz_x = cbuffer.data(ip_geom_20_off + 480 * ccomps * dcomps);

            auto g_zz_0_xzzzzz_y = cbuffer.data(ip_geom_20_off + 481 * ccomps * dcomps);

            auto g_zz_0_xzzzzz_z = cbuffer.data(ip_geom_20_off + 482 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_xzzzzz_x, g_zz_0_xzzzzz_y, g_zz_0_xzzzzz_z, g_zz_0_zzzzz_x, g_zz_0_zzzzz_xx, g_zz_0_zzzzz_xy, g_zz_0_zzzzz_xz, g_zz_0_zzzzz_y, g_zz_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_xzzzzz_x[k] = -g_zz_0_zzzzz_x[k] * ab_x + g_zz_0_zzzzz_xx[k];

                g_zz_0_xzzzzz_y[k] = -g_zz_0_zzzzz_y[k] * ab_x + g_zz_0_zzzzz_xy[k];

                g_zz_0_xzzzzz_z[k] = -g_zz_0_zzzzz_z[k] * ab_x + g_zz_0_zzzzz_xz[k];
            }

            /// Set up 483-486 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyyy_x = cbuffer.data(ip_geom_20_off + 483 * ccomps * dcomps);

            auto g_zz_0_yyyyyy_y = cbuffer.data(ip_geom_20_off + 484 * ccomps * dcomps);

            auto g_zz_0_yyyyyy_z = cbuffer.data(ip_geom_20_off + 485 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyy_x, g_zz_0_yyyyy_xy, g_zz_0_yyyyy_y, g_zz_0_yyyyy_yy, g_zz_0_yyyyy_yz, g_zz_0_yyyyy_z, g_zz_0_yyyyyy_x, g_zz_0_yyyyyy_y, g_zz_0_yyyyyy_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyyy_x[k] = -g_zz_0_yyyyy_x[k] * ab_y + g_zz_0_yyyyy_xy[k];

                g_zz_0_yyyyyy_y[k] = -g_zz_0_yyyyy_y[k] * ab_y + g_zz_0_yyyyy_yy[k];

                g_zz_0_yyyyyy_z[k] = -g_zz_0_yyyyy_z[k] * ab_y + g_zz_0_yyyyy_yz[k];
            }

            /// Set up 486-489 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyyz_x = cbuffer.data(ip_geom_20_off + 486 * ccomps * dcomps);

            auto g_zz_0_yyyyyz_y = cbuffer.data(ip_geom_20_off + 487 * ccomps * dcomps);

            auto g_zz_0_yyyyyz_z = cbuffer.data(ip_geom_20_off + 488 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyyz_x, g_zz_0_yyyyyz_y, g_zz_0_yyyyyz_z, g_zz_0_yyyyz_x, g_zz_0_yyyyz_xy, g_zz_0_yyyyz_y, g_zz_0_yyyyz_yy, g_zz_0_yyyyz_yz, g_zz_0_yyyyz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyyz_x[k] = -g_zz_0_yyyyz_x[k] * ab_y + g_zz_0_yyyyz_xy[k];

                g_zz_0_yyyyyz_y[k] = -g_zz_0_yyyyz_y[k] * ab_y + g_zz_0_yyyyz_yy[k];

                g_zz_0_yyyyyz_z[k] = -g_zz_0_yyyyz_z[k] * ab_y + g_zz_0_yyyyz_yz[k];
            }

            /// Set up 489-492 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyyzz_x = cbuffer.data(ip_geom_20_off + 489 * ccomps * dcomps);

            auto g_zz_0_yyyyzz_y = cbuffer.data(ip_geom_20_off + 490 * ccomps * dcomps);

            auto g_zz_0_yyyyzz_z = cbuffer.data(ip_geom_20_off + 491 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyyzz_x, g_zz_0_yyyyzz_y, g_zz_0_yyyyzz_z, g_zz_0_yyyzz_x, g_zz_0_yyyzz_xy, g_zz_0_yyyzz_y, g_zz_0_yyyzz_yy, g_zz_0_yyyzz_yz, g_zz_0_yyyzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyyzz_x[k] = -g_zz_0_yyyzz_x[k] * ab_y + g_zz_0_yyyzz_xy[k];

                g_zz_0_yyyyzz_y[k] = -g_zz_0_yyyzz_y[k] * ab_y + g_zz_0_yyyzz_yy[k];

                g_zz_0_yyyyzz_z[k] = -g_zz_0_yyyzz_z[k] * ab_y + g_zz_0_yyyzz_yz[k];
            }

            /// Set up 492-495 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyyzzz_x = cbuffer.data(ip_geom_20_off + 492 * ccomps * dcomps);

            auto g_zz_0_yyyzzz_y = cbuffer.data(ip_geom_20_off + 493 * ccomps * dcomps);

            auto g_zz_0_yyyzzz_z = cbuffer.data(ip_geom_20_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyyzzz_x, g_zz_0_yyyzzz_y, g_zz_0_yyyzzz_z, g_zz_0_yyzzz_x, g_zz_0_yyzzz_xy, g_zz_0_yyzzz_y, g_zz_0_yyzzz_yy, g_zz_0_yyzzz_yz, g_zz_0_yyzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyyzzz_x[k] = -g_zz_0_yyzzz_x[k] * ab_y + g_zz_0_yyzzz_xy[k];

                g_zz_0_yyyzzz_y[k] = -g_zz_0_yyzzz_y[k] * ab_y + g_zz_0_yyzzz_yy[k];

                g_zz_0_yyyzzz_z[k] = -g_zz_0_yyzzz_z[k] * ab_y + g_zz_0_yyzzz_yz[k];
            }

            /// Set up 495-498 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yyzzzz_x = cbuffer.data(ip_geom_20_off + 495 * ccomps * dcomps);

            auto g_zz_0_yyzzzz_y = cbuffer.data(ip_geom_20_off + 496 * ccomps * dcomps);

            auto g_zz_0_yyzzzz_z = cbuffer.data(ip_geom_20_off + 497 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yyzzzz_x, g_zz_0_yyzzzz_y, g_zz_0_yyzzzz_z, g_zz_0_yzzzz_x, g_zz_0_yzzzz_xy, g_zz_0_yzzzz_y, g_zz_0_yzzzz_yy, g_zz_0_yzzzz_yz, g_zz_0_yzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yyzzzz_x[k] = -g_zz_0_yzzzz_x[k] * ab_y + g_zz_0_yzzzz_xy[k];

                g_zz_0_yyzzzz_y[k] = -g_zz_0_yzzzz_y[k] * ab_y + g_zz_0_yzzzz_yy[k];

                g_zz_0_yyzzzz_z[k] = -g_zz_0_yzzzz_z[k] * ab_y + g_zz_0_yzzzz_yz[k];
            }

            /// Set up 498-501 components of targeted buffer : cbuffer.data(

            auto g_zz_0_yzzzzz_x = cbuffer.data(ip_geom_20_off + 498 * ccomps * dcomps);

            auto g_zz_0_yzzzzz_y = cbuffer.data(ip_geom_20_off + 499 * ccomps * dcomps);

            auto g_zz_0_yzzzzz_z = cbuffer.data(ip_geom_20_off + 500 * ccomps * dcomps);

            #pragma omp simd aligned(g_zz_0_yzzzzz_x, g_zz_0_yzzzzz_y, g_zz_0_yzzzzz_z, g_zz_0_zzzzz_x, g_zz_0_zzzzz_xy, g_zz_0_zzzzz_y, g_zz_0_zzzzz_yy, g_zz_0_zzzzz_yz, g_zz_0_zzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_yzzzzz_x[k] = -g_zz_0_zzzzz_x[k] * ab_y + g_zz_0_zzzzz_xy[k];

                g_zz_0_yzzzzz_y[k] = -g_zz_0_zzzzz_y[k] * ab_y + g_zz_0_zzzzz_yy[k];

                g_zz_0_yzzzzz_z[k] = -g_zz_0_zzzzz_z[k] * ab_y + g_zz_0_zzzzz_yz[k];
            }

            /// Set up 501-504 components of targeted buffer : cbuffer.data(

            auto g_zz_0_zzzzzz_x = cbuffer.data(ip_geom_20_off + 501 * ccomps * dcomps);

            auto g_zz_0_zzzzzz_y = cbuffer.data(ip_geom_20_off + 502 * ccomps * dcomps);

            auto g_zz_0_zzzzzz_z = cbuffer.data(ip_geom_20_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzzz_x, g_z_0_zzzzz_y, g_z_0_zzzzz_z, g_zz_0_zzzzz_x, g_zz_0_zzzzz_xz, g_zz_0_zzzzz_y, g_zz_0_zzzzz_yz, g_zz_0_zzzzz_z, g_zz_0_zzzzz_zz, g_zz_0_zzzzzz_x, g_zz_0_zzzzzz_y, g_zz_0_zzzzzz_z  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_zz_0_zzzzzz_x[k] = -2.0 * g_z_0_zzzzz_x[k] - g_zz_0_zzzzz_x[k] * ab_z + g_zz_0_zzzzz_xz[k];

                g_zz_0_zzzzzz_y[k] = -2.0 * g_z_0_zzzzz_y[k] - g_zz_0_zzzzz_y[k] * ab_z + g_zz_0_zzzzz_yz[k];

                g_zz_0_zzzzzz_z[k] = -2.0 * g_z_0_zzzzz_z[k] - g_zz_0_zzzzz_z[k] * ab_z + g_zz_0_zzzzz_zz[k];
            }
        }
    }
}

} // erirec namespace

