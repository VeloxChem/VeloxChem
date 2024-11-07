#include "ElectronRepulsionGeom1000ContrRecIFXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_ifxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_ifxx,
                                            const size_t idx_hfxx,
                                            const size_t idx_geom_10_hfxx,
                                            const size_t idx_geom_10_hgxx,
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
            /// Set up components of auxilary buffer : HFSS

            const auto hf_off = idx_hfxx + i * dcomps + j;

            auto g_xxxxx_xxx = cbuffer.data(hf_off + 0 * ccomps * dcomps);

            auto g_xxxxx_xxy = cbuffer.data(hf_off + 1 * ccomps * dcomps);

            auto g_xxxxx_xxz = cbuffer.data(hf_off + 2 * ccomps * dcomps);

            auto g_xxxxx_xyy = cbuffer.data(hf_off + 3 * ccomps * dcomps);

            auto g_xxxxx_xyz = cbuffer.data(hf_off + 4 * ccomps * dcomps);

            auto g_xxxxx_xzz = cbuffer.data(hf_off + 5 * ccomps * dcomps);

            auto g_xxxxx_yyy = cbuffer.data(hf_off + 6 * ccomps * dcomps);

            auto g_xxxxx_yyz = cbuffer.data(hf_off + 7 * ccomps * dcomps);

            auto g_xxxxx_yzz = cbuffer.data(hf_off + 8 * ccomps * dcomps);

            auto g_xxxxx_zzz = cbuffer.data(hf_off + 9 * ccomps * dcomps);

            auto g_xxxxy_xxx = cbuffer.data(hf_off + 10 * ccomps * dcomps);

            auto g_xxxxy_xxy = cbuffer.data(hf_off + 11 * ccomps * dcomps);

            auto g_xxxxy_xxz = cbuffer.data(hf_off + 12 * ccomps * dcomps);

            auto g_xxxxy_xyy = cbuffer.data(hf_off + 13 * ccomps * dcomps);

            auto g_xxxxy_xyz = cbuffer.data(hf_off + 14 * ccomps * dcomps);

            auto g_xxxxy_xzz = cbuffer.data(hf_off + 15 * ccomps * dcomps);

            auto g_xxxxy_yyy = cbuffer.data(hf_off + 16 * ccomps * dcomps);

            auto g_xxxxy_yyz = cbuffer.data(hf_off + 17 * ccomps * dcomps);

            auto g_xxxxy_yzz = cbuffer.data(hf_off + 18 * ccomps * dcomps);

            auto g_xxxxy_zzz = cbuffer.data(hf_off + 19 * ccomps * dcomps);

            auto g_xxxxz_xxx = cbuffer.data(hf_off + 20 * ccomps * dcomps);

            auto g_xxxxz_xxy = cbuffer.data(hf_off + 21 * ccomps * dcomps);

            auto g_xxxxz_xxz = cbuffer.data(hf_off + 22 * ccomps * dcomps);

            auto g_xxxxz_xyy = cbuffer.data(hf_off + 23 * ccomps * dcomps);

            auto g_xxxxz_xyz = cbuffer.data(hf_off + 24 * ccomps * dcomps);

            auto g_xxxxz_xzz = cbuffer.data(hf_off + 25 * ccomps * dcomps);

            auto g_xxxxz_yyy = cbuffer.data(hf_off + 26 * ccomps * dcomps);

            auto g_xxxxz_yyz = cbuffer.data(hf_off + 27 * ccomps * dcomps);

            auto g_xxxxz_yzz = cbuffer.data(hf_off + 28 * ccomps * dcomps);

            auto g_xxxxz_zzz = cbuffer.data(hf_off + 29 * ccomps * dcomps);

            auto g_xxxyy_xxx = cbuffer.data(hf_off + 30 * ccomps * dcomps);

            auto g_xxxyy_xxy = cbuffer.data(hf_off + 31 * ccomps * dcomps);

            auto g_xxxyy_xxz = cbuffer.data(hf_off + 32 * ccomps * dcomps);

            auto g_xxxyy_xyy = cbuffer.data(hf_off + 33 * ccomps * dcomps);

            auto g_xxxyy_xyz = cbuffer.data(hf_off + 34 * ccomps * dcomps);

            auto g_xxxyy_xzz = cbuffer.data(hf_off + 35 * ccomps * dcomps);

            auto g_xxxyy_yyy = cbuffer.data(hf_off + 36 * ccomps * dcomps);

            auto g_xxxyy_yyz = cbuffer.data(hf_off + 37 * ccomps * dcomps);

            auto g_xxxyy_yzz = cbuffer.data(hf_off + 38 * ccomps * dcomps);

            auto g_xxxyy_zzz = cbuffer.data(hf_off + 39 * ccomps * dcomps);

            auto g_xxxyz_xxx = cbuffer.data(hf_off + 40 * ccomps * dcomps);

            auto g_xxxyz_xxy = cbuffer.data(hf_off + 41 * ccomps * dcomps);

            auto g_xxxyz_xxz = cbuffer.data(hf_off + 42 * ccomps * dcomps);

            auto g_xxxyz_xyy = cbuffer.data(hf_off + 43 * ccomps * dcomps);

            auto g_xxxyz_xyz = cbuffer.data(hf_off + 44 * ccomps * dcomps);

            auto g_xxxyz_xzz = cbuffer.data(hf_off + 45 * ccomps * dcomps);

            auto g_xxxyz_yyy = cbuffer.data(hf_off + 46 * ccomps * dcomps);

            auto g_xxxyz_yyz = cbuffer.data(hf_off + 47 * ccomps * dcomps);

            auto g_xxxyz_yzz = cbuffer.data(hf_off + 48 * ccomps * dcomps);

            auto g_xxxyz_zzz = cbuffer.data(hf_off + 49 * ccomps * dcomps);

            auto g_xxxzz_xxx = cbuffer.data(hf_off + 50 * ccomps * dcomps);

            auto g_xxxzz_xxy = cbuffer.data(hf_off + 51 * ccomps * dcomps);

            auto g_xxxzz_xxz = cbuffer.data(hf_off + 52 * ccomps * dcomps);

            auto g_xxxzz_xyy = cbuffer.data(hf_off + 53 * ccomps * dcomps);

            auto g_xxxzz_xyz = cbuffer.data(hf_off + 54 * ccomps * dcomps);

            auto g_xxxzz_xzz = cbuffer.data(hf_off + 55 * ccomps * dcomps);

            auto g_xxxzz_yyy = cbuffer.data(hf_off + 56 * ccomps * dcomps);

            auto g_xxxzz_yyz = cbuffer.data(hf_off + 57 * ccomps * dcomps);

            auto g_xxxzz_yzz = cbuffer.data(hf_off + 58 * ccomps * dcomps);

            auto g_xxxzz_zzz = cbuffer.data(hf_off + 59 * ccomps * dcomps);

            auto g_xxyyy_xxx = cbuffer.data(hf_off + 60 * ccomps * dcomps);

            auto g_xxyyy_xxy = cbuffer.data(hf_off + 61 * ccomps * dcomps);

            auto g_xxyyy_xxz = cbuffer.data(hf_off + 62 * ccomps * dcomps);

            auto g_xxyyy_xyy = cbuffer.data(hf_off + 63 * ccomps * dcomps);

            auto g_xxyyy_xyz = cbuffer.data(hf_off + 64 * ccomps * dcomps);

            auto g_xxyyy_xzz = cbuffer.data(hf_off + 65 * ccomps * dcomps);

            auto g_xxyyy_yyy = cbuffer.data(hf_off + 66 * ccomps * dcomps);

            auto g_xxyyy_yyz = cbuffer.data(hf_off + 67 * ccomps * dcomps);

            auto g_xxyyy_yzz = cbuffer.data(hf_off + 68 * ccomps * dcomps);

            auto g_xxyyy_zzz = cbuffer.data(hf_off + 69 * ccomps * dcomps);

            auto g_xxyyz_xxx = cbuffer.data(hf_off + 70 * ccomps * dcomps);

            auto g_xxyyz_xxy = cbuffer.data(hf_off + 71 * ccomps * dcomps);

            auto g_xxyyz_xxz = cbuffer.data(hf_off + 72 * ccomps * dcomps);

            auto g_xxyyz_xyy = cbuffer.data(hf_off + 73 * ccomps * dcomps);

            auto g_xxyyz_xyz = cbuffer.data(hf_off + 74 * ccomps * dcomps);

            auto g_xxyyz_xzz = cbuffer.data(hf_off + 75 * ccomps * dcomps);

            auto g_xxyyz_yyy = cbuffer.data(hf_off + 76 * ccomps * dcomps);

            auto g_xxyyz_yyz = cbuffer.data(hf_off + 77 * ccomps * dcomps);

            auto g_xxyyz_yzz = cbuffer.data(hf_off + 78 * ccomps * dcomps);

            auto g_xxyyz_zzz = cbuffer.data(hf_off + 79 * ccomps * dcomps);

            auto g_xxyzz_xxx = cbuffer.data(hf_off + 80 * ccomps * dcomps);

            auto g_xxyzz_xxy = cbuffer.data(hf_off + 81 * ccomps * dcomps);

            auto g_xxyzz_xxz = cbuffer.data(hf_off + 82 * ccomps * dcomps);

            auto g_xxyzz_xyy = cbuffer.data(hf_off + 83 * ccomps * dcomps);

            auto g_xxyzz_xyz = cbuffer.data(hf_off + 84 * ccomps * dcomps);

            auto g_xxyzz_xzz = cbuffer.data(hf_off + 85 * ccomps * dcomps);

            auto g_xxyzz_yyy = cbuffer.data(hf_off + 86 * ccomps * dcomps);

            auto g_xxyzz_yyz = cbuffer.data(hf_off + 87 * ccomps * dcomps);

            auto g_xxyzz_yzz = cbuffer.data(hf_off + 88 * ccomps * dcomps);

            auto g_xxyzz_zzz = cbuffer.data(hf_off + 89 * ccomps * dcomps);

            auto g_xxzzz_xxx = cbuffer.data(hf_off + 90 * ccomps * dcomps);

            auto g_xxzzz_xxy = cbuffer.data(hf_off + 91 * ccomps * dcomps);

            auto g_xxzzz_xxz = cbuffer.data(hf_off + 92 * ccomps * dcomps);

            auto g_xxzzz_xyy = cbuffer.data(hf_off + 93 * ccomps * dcomps);

            auto g_xxzzz_xyz = cbuffer.data(hf_off + 94 * ccomps * dcomps);

            auto g_xxzzz_xzz = cbuffer.data(hf_off + 95 * ccomps * dcomps);

            auto g_xxzzz_yyy = cbuffer.data(hf_off + 96 * ccomps * dcomps);

            auto g_xxzzz_yyz = cbuffer.data(hf_off + 97 * ccomps * dcomps);

            auto g_xxzzz_yzz = cbuffer.data(hf_off + 98 * ccomps * dcomps);

            auto g_xxzzz_zzz = cbuffer.data(hf_off + 99 * ccomps * dcomps);

            auto g_xyyyy_xxx = cbuffer.data(hf_off + 100 * ccomps * dcomps);

            auto g_xyyyy_xxy = cbuffer.data(hf_off + 101 * ccomps * dcomps);

            auto g_xyyyy_xxz = cbuffer.data(hf_off + 102 * ccomps * dcomps);

            auto g_xyyyy_xyy = cbuffer.data(hf_off + 103 * ccomps * dcomps);

            auto g_xyyyy_xyz = cbuffer.data(hf_off + 104 * ccomps * dcomps);

            auto g_xyyyy_xzz = cbuffer.data(hf_off + 105 * ccomps * dcomps);

            auto g_xyyyy_yyy = cbuffer.data(hf_off + 106 * ccomps * dcomps);

            auto g_xyyyy_yyz = cbuffer.data(hf_off + 107 * ccomps * dcomps);

            auto g_xyyyy_yzz = cbuffer.data(hf_off + 108 * ccomps * dcomps);

            auto g_xyyyy_zzz = cbuffer.data(hf_off + 109 * ccomps * dcomps);

            auto g_xyyyz_xxx = cbuffer.data(hf_off + 110 * ccomps * dcomps);

            auto g_xyyyz_xxy = cbuffer.data(hf_off + 111 * ccomps * dcomps);

            auto g_xyyyz_xxz = cbuffer.data(hf_off + 112 * ccomps * dcomps);

            auto g_xyyyz_xyy = cbuffer.data(hf_off + 113 * ccomps * dcomps);

            auto g_xyyyz_xyz = cbuffer.data(hf_off + 114 * ccomps * dcomps);

            auto g_xyyyz_xzz = cbuffer.data(hf_off + 115 * ccomps * dcomps);

            auto g_xyyyz_yyy = cbuffer.data(hf_off + 116 * ccomps * dcomps);

            auto g_xyyyz_yyz = cbuffer.data(hf_off + 117 * ccomps * dcomps);

            auto g_xyyyz_yzz = cbuffer.data(hf_off + 118 * ccomps * dcomps);

            auto g_xyyyz_zzz = cbuffer.data(hf_off + 119 * ccomps * dcomps);

            auto g_xyyzz_xxx = cbuffer.data(hf_off + 120 * ccomps * dcomps);

            auto g_xyyzz_xxy = cbuffer.data(hf_off + 121 * ccomps * dcomps);

            auto g_xyyzz_xxz = cbuffer.data(hf_off + 122 * ccomps * dcomps);

            auto g_xyyzz_xyy = cbuffer.data(hf_off + 123 * ccomps * dcomps);

            auto g_xyyzz_xyz = cbuffer.data(hf_off + 124 * ccomps * dcomps);

            auto g_xyyzz_xzz = cbuffer.data(hf_off + 125 * ccomps * dcomps);

            auto g_xyyzz_yyy = cbuffer.data(hf_off + 126 * ccomps * dcomps);

            auto g_xyyzz_yyz = cbuffer.data(hf_off + 127 * ccomps * dcomps);

            auto g_xyyzz_yzz = cbuffer.data(hf_off + 128 * ccomps * dcomps);

            auto g_xyyzz_zzz = cbuffer.data(hf_off + 129 * ccomps * dcomps);

            auto g_xyzzz_xxx = cbuffer.data(hf_off + 130 * ccomps * dcomps);

            auto g_xyzzz_xxy = cbuffer.data(hf_off + 131 * ccomps * dcomps);

            auto g_xyzzz_xxz = cbuffer.data(hf_off + 132 * ccomps * dcomps);

            auto g_xyzzz_xyy = cbuffer.data(hf_off + 133 * ccomps * dcomps);

            auto g_xyzzz_xyz = cbuffer.data(hf_off + 134 * ccomps * dcomps);

            auto g_xyzzz_xzz = cbuffer.data(hf_off + 135 * ccomps * dcomps);

            auto g_xyzzz_yyy = cbuffer.data(hf_off + 136 * ccomps * dcomps);

            auto g_xyzzz_yyz = cbuffer.data(hf_off + 137 * ccomps * dcomps);

            auto g_xyzzz_yzz = cbuffer.data(hf_off + 138 * ccomps * dcomps);

            auto g_xyzzz_zzz = cbuffer.data(hf_off + 139 * ccomps * dcomps);

            auto g_xzzzz_xxx = cbuffer.data(hf_off + 140 * ccomps * dcomps);

            auto g_xzzzz_xxy = cbuffer.data(hf_off + 141 * ccomps * dcomps);

            auto g_xzzzz_xxz = cbuffer.data(hf_off + 142 * ccomps * dcomps);

            auto g_xzzzz_xyy = cbuffer.data(hf_off + 143 * ccomps * dcomps);

            auto g_xzzzz_xyz = cbuffer.data(hf_off + 144 * ccomps * dcomps);

            auto g_xzzzz_xzz = cbuffer.data(hf_off + 145 * ccomps * dcomps);

            auto g_xzzzz_yyy = cbuffer.data(hf_off + 146 * ccomps * dcomps);

            auto g_xzzzz_yyz = cbuffer.data(hf_off + 147 * ccomps * dcomps);

            auto g_xzzzz_yzz = cbuffer.data(hf_off + 148 * ccomps * dcomps);

            auto g_xzzzz_zzz = cbuffer.data(hf_off + 149 * ccomps * dcomps);

            auto g_yyyyy_xxx = cbuffer.data(hf_off + 150 * ccomps * dcomps);

            auto g_yyyyy_xxy = cbuffer.data(hf_off + 151 * ccomps * dcomps);

            auto g_yyyyy_xxz = cbuffer.data(hf_off + 152 * ccomps * dcomps);

            auto g_yyyyy_xyy = cbuffer.data(hf_off + 153 * ccomps * dcomps);

            auto g_yyyyy_xyz = cbuffer.data(hf_off + 154 * ccomps * dcomps);

            auto g_yyyyy_xzz = cbuffer.data(hf_off + 155 * ccomps * dcomps);

            auto g_yyyyy_yyy = cbuffer.data(hf_off + 156 * ccomps * dcomps);

            auto g_yyyyy_yyz = cbuffer.data(hf_off + 157 * ccomps * dcomps);

            auto g_yyyyy_yzz = cbuffer.data(hf_off + 158 * ccomps * dcomps);

            auto g_yyyyy_zzz = cbuffer.data(hf_off + 159 * ccomps * dcomps);

            auto g_yyyyz_xxx = cbuffer.data(hf_off + 160 * ccomps * dcomps);

            auto g_yyyyz_xxy = cbuffer.data(hf_off + 161 * ccomps * dcomps);

            auto g_yyyyz_xxz = cbuffer.data(hf_off + 162 * ccomps * dcomps);

            auto g_yyyyz_xyy = cbuffer.data(hf_off + 163 * ccomps * dcomps);

            auto g_yyyyz_xyz = cbuffer.data(hf_off + 164 * ccomps * dcomps);

            auto g_yyyyz_xzz = cbuffer.data(hf_off + 165 * ccomps * dcomps);

            auto g_yyyyz_yyy = cbuffer.data(hf_off + 166 * ccomps * dcomps);

            auto g_yyyyz_yyz = cbuffer.data(hf_off + 167 * ccomps * dcomps);

            auto g_yyyyz_yzz = cbuffer.data(hf_off + 168 * ccomps * dcomps);

            auto g_yyyyz_zzz = cbuffer.data(hf_off + 169 * ccomps * dcomps);

            auto g_yyyzz_xxx = cbuffer.data(hf_off + 170 * ccomps * dcomps);

            auto g_yyyzz_xxy = cbuffer.data(hf_off + 171 * ccomps * dcomps);

            auto g_yyyzz_xxz = cbuffer.data(hf_off + 172 * ccomps * dcomps);

            auto g_yyyzz_xyy = cbuffer.data(hf_off + 173 * ccomps * dcomps);

            auto g_yyyzz_xyz = cbuffer.data(hf_off + 174 * ccomps * dcomps);

            auto g_yyyzz_xzz = cbuffer.data(hf_off + 175 * ccomps * dcomps);

            auto g_yyyzz_yyy = cbuffer.data(hf_off + 176 * ccomps * dcomps);

            auto g_yyyzz_yyz = cbuffer.data(hf_off + 177 * ccomps * dcomps);

            auto g_yyyzz_yzz = cbuffer.data(hf_off + 178 * ccomps * dcomps);

            auto g_yyyzz_zzz = cbuffer.data(hf_off + 179 * ccomps * dcomps);

            auto g_yyzzz_xxx = cbuffer.data(hf_off + 180 * ccomps * dcomps);

            auto g_yyzzz_xxy = cbuffer.data(hf_off + 181 * ccomps * dcomps);

            auto g_yyzzz_xxz = cbuffer.data(hf_off + 182 * ccomps * dcomps);

            auto g_yyzzz_xyy = cbuffer.data(hf_off + 183 * ccomps * dcomps);

            auto g_yyzzz_xyz = cbuffer.data(hf_off + 184 * ccomps * dcomps);

            auto g_yyzzz_xzz = cbuffer.data(hf_off + 185 * ccomps * dcomps);

            auto g_yyzzz_yyy = cbuffer.data(hf_off + 186 * ccomps * dcomps);

            auto g_yyzzz_yyz = cbuffer.data(hf_off + 187 * ccomps * dcomps);

            auto g_yyzzz_yzz = cbuffer.data(hf_off + 188 * ccomps * dcomps);

            auto g_yyzzz_zzz = cbuffer.data(hf_off + 189 * ccomps * dcomps);

            auto g_yzzzz_xxx = cbuffer.data(hf_off + 190 * ccomps * dcomps);

            auto g_yzzzz_xxy = cbuffer.data(hf_off + 191 * ccomps * dcomps);

            auto g_yzzzz_xxz = cbuffer.data(hf_off + 192 * ccomps * dcomps);

            auto g_yzzzz_xyy = cbuffer.data(hf_off + 193 * ccomps * dcomps);

            auto g_yzzzz_xyz = cbuffer.data(hf_off + 194 * ccomps * dcomps);

            auto g_yzzzz_xzz = cbuffer.data(hf_off + 195 * ccomps * dcomps);

            auto g_yzzzz_yyy = cbuffer.data(hf_off + 196 * ccomps * dcomps);

            auto g_yzzzz_yyz = cbuffer.data(hf_off + 197 * ccomps * dcomps);

            auto g_yzzzz_yzz = cbuffer.data(hf_off + 198 * ccomps * dcomps);

            auto g_yzzzz_zzz = cbuffer.data(hf_off + 199 * ccomps * dcomps);

            auto g_zzzzz_xxx = cbuffer.data(hf_off + 200 * ccomps * dcomps);

            auto g_zzzzz_xxy = cbuffer.data(hf_off + 201 * ccomps * dcomps);

            auto g_zzzzz_xxz = cbuffer.data(hf_off + 202 * ccomps * dcomps);

            auto g_zzzzz_xyy = cbuffer.data(hf_off + 203 * ccomps * dcomps);

            auto g_zzzzz_xyz = cbuffer.data(hf_off + 204 * ccomps * dcomps);

            auto g_zzzzz_xzz = cbuffer.data(hf_off + 205 * ccomps * dcomps);

            auto g_zzzzz_yyy = cbuffer.data(hf_off + 206 * ccomps * dcomps);

            auto g_zzzzz_yyz = cbuffer.data(hf_off + 207 * ccomps * dcomps);

            auto g_zzzzz_yzz = cbuffer.data(hf_off + 208 * ccomps * dcomps);

            auto g_zzzzz_zzz = cbuffer.data(hf_off + 209 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HFSS

            const auto hf_geom_10_off = idx_geom_10_hfxx + i * dcomps + j;

            auto g_x_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 209 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 210 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 211 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 212 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 213 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 214 * ccomps * dcomps);

            auto g_y_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 215 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 216 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 217 * ccomps * dcomps);

            auto g_y_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 218 * ccomps * dcomps);

            auto g_y_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 219 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 220 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 221 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 222 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 223 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 224 * ccomps * dcomps);

            auto g_y_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 225 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 226 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 227 * ccomps * dcomps);

            auto g_y_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 228 * ccomps * dcomps);

            auto g_y_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 229 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 230 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 231 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 232 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 233 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 234 * ccomps * dcomps);

            auto g_y_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 235 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 236 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 237 * ccomps * dcomps);

            auto g_y_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 238 * ccomps * dcomps);

            auto g_y_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 239 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 240 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 241 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 242 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 243 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 244 * ccomps * dcomps);

            auto g_y_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 245 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 246 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 247 * ccomps * dcomps);

            auto g_y_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 248 * ccomps * dcomps);

            auto g_y_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 249 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 250 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 251 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 252 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 253 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 254 * ccomps * dcomps);

            auto g_y_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 255 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 256 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 257 * ccomps * dcomps);

            auto g_y_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 258 * ccomps * dcomps);

            auto g_y_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 259 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 260 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 261 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 262 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 263 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 264 * ccomps * dcomps);

            auto g_y_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 265 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 266 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 267 * ccomps * dcomps);

            auto g_y_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 268 * ccomps * dcomps);

            auto g_y_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 269 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 270 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 271 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 272 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 273 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 274 * ccomps * dcomps);

            auto g_y_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 275 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 276 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 277 * ccomps * dcomps);

            auto g_y_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 278 * ccomps * dcomps);

            auto g_y_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 279 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 289 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 299 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 307 * ccomps * dcomps);

            auto g_y_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 309 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 329 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 335 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 339 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 349 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 359 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 389 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 419 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxx = cbuffer.data(hf_geom_10_off + 420 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxy = cbuffer.data(hf_geom_10_off + 421 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxz = cbuffer.data(hf_geom_10_off + 422 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyy = cbuffer.data(hf_geom_10_off + 423 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyz = cbuffer.data(hf_geom_10_off + 424 * ccomps * dcomps);

            auto g_z_0_xxxxx_xzz = cbuffer.data(hf_geom_10_off + 425 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyy = cbuffer.data(hf_geom_10_off + 426 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyz = cbuffer.data(hf_geom_10_off + 427 * ccomps * dcomps);

            auto g_z_0_xxxxx_yzz = cbuffer.data(hf_geom_10_off + 428 * ccomps * dcomps);

            auto g_z_0_xxxxx_zzz = cbuffer.data(hf_geom_10_off + 429 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxx = cbuffer.data(hf_geom_10_off + 430 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxy = cbuffer.data(hf_geom_10_off + 431 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxz = cbuffer.data(hf_geom_10_off + 432 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyy = cbuffer.data(hf_geom_10_off + 433 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyz = cbuffer.data(hf_geom_10_off + 434 * ccomps * dcomps);

            auto g_z_0_xxxxy_xzz = cbuffer.data(hf_geom_10_off + 435 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyy = cbuffer.data(hf_geom_10_off + 436 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyz = cbuffer.data(hf_geom_10_off + 437 * ccomps * dcomps);

            auto g_z_0_xxxxy_yzz = cbuffer.data(hf_geom_10_off + 438 * ccomps * dcomps);

            auto g_z_0_xxxxy_zzz = cbuffer.data(hf_geom_10_off + 439 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxx = cbuffer.data(hf_geom_10_off + 440 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxy = cbuffer.data(hf_geom_10_off + 441 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxz = cbuffer.data(hf_geom_10_off + 442 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyy = cbuffer.data(hf_geom_10_off + 443 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyz = cbuffer.data(hf_geom_10_off + 444 * ccomps * dcomps);

            auto g_z_0_xxxxz_xzz = cbuffer.data(hf_geom_10_off + 445 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyy = cbuffer.data(hf_geom_10_off + 446 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyz = cbuffer.data(hf_geom_10_off + 447 * ccomps * dcomps);

            auto g_z_0_xxxxz_yzz = cbuffer.data(hf_geom_10_off + 448 * ccomps * dcomps);

            auto g_z_0_xxxxz_zzz = cbuffer.data(hf_geom_10_off + 449 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxx = cbuffer.data(hf_geom_10_off + 450 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxy = cbuffer.data(hf_geom_10_off + 451 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxz = cbuffer.data(hf_geom_10_off + 452 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyy = cbuffer.data(hf_geom_10_off + 453 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyz = cbuffer.data(hf_geom_10_off + 454 * ccomps * dcomps);

            auto g_z_0_xxxyy_xzz = cbuffer.data(hf_geom_10_off + 455 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyy = cbuffer.data(hf_geom_10_off + 456 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyz = cbuffer.data(hf_geom_10_off + 457 * ccomps * dcomps);

            auto g_z_0_xxxyy_yzz = cbuffer.data(hf_geom_10_off + 458 * ccomps * dcomps);

            auto g_z_0_xxxyy_zzz = cbuffer.data(hf_geom_10_off + 459 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxx = cbuffer.data(hf_geom_10_off + 460 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxy = cbuffer.data(hf_geom_10_off + 461 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxz = cbuffer.data(hf_geom_10_off + 462 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyy = cbuffer.data(hf_geom_10_off + 463 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyz = cbuffer.data(hf_geom_10_off + 464 * ccomps * dcomps);

            auto g_z_0_xxxyz_xzz = cbuffer.data(hf_geom_10_off + 465 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyy = cbuffer.data(hf_geom_10_off + 466 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyz = cbuffer.data(hf_geom_10_off + 467 * ccomps * dcomps);

            auto g_z_0_xxxyz_yzz = cbuffer.data(hf_geom_10_off + 468 * ccomps * dcomps);

            auto g_z_0_xxxyz_zzz = cbuffer.data(hf_geom_10_off + 469 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxx = cbuffer.data(hf_geom_10_off + 470 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxy = cbuffer.data(hf_geom_10_off + 471 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxz = cbuffer.data(hf_geom_10_off + 472 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyy = cbuffer.data(hf_geom_10_off + 473 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyz = cbuffer.data(hf_geom_10_off + 474 * ccomps * dcomps);

            auto g_z_0_xxxzz_xzz = cbuffer.data(hf_geom_10_off + 475 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyy = cbuffer.data(hf_geom_10_off + 476 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyz = cbuffer.data(hf_geom_10_off + 477 * ccomps * dcomps);

            auto g_z_0_xxxzz_yzz = cbuffer.data(hf_geom_10_off + 478 * ccomps * dcomps);

            auto g_z_0_xxxzz_zzz = cbuffer.data(hf_geom_10_off + 479 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxx = cbuffer.data(hf_geom_10_off + 480 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxy = cbuffer.data(hf_geom_10_off + 481 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxz = cbuffer.data(hf_geom_10_off + 482 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyy = cbuffer.data(hf_geom_10_off + 483 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyz = cbuffer.data(hf_geom_10_off + 484 * ccomps * dcomps);

            auto g_z_0_xxyyy_xzz = cbuffer.data(hf_geom_10_off + 485 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyy = cbuffer.data(hf_geom_10_off + 486 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyz = cbuffer.data(hf_geom_10_off + 487 * ccomps * dcomps);

            auto g_z_0_xxyyy_yzz = cbuffer.data(hf_geom_10_off + 488 * ccomps * dcomps);

            auto g_z_0_xxyyy_zzz = cbuffer.data(hf_geom_10_off + 489 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxx = cbuffer.data(hf_geom_10_off + 490 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxy = cbuffer.data(hf_geom_10_off + 491 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxz = cbuffer.data(hf_geom_10_off + 492 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyy = cbuffer.data(hf_geom_10_off + 493 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyz = cbuffer.data(hf_geom_10_off + 494 * ccomps * dcomps);

            auto g_z_0_xxyyz_xzz = cbuffer.data(hf_geom_10_off + 495 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyy = cbuffer.data(hf_geom_10_off + 496 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyz = cbuffer.data(hf_geom_10_off + 497 * ccomps * dcomps);

            auto g_z_0_xxyyz_yzz = cbuffer.data(hf_geom_10_off + 498 * ccomps * dcomps);

            auto g_z_0_xxyyz_zzz = cbuffer.data(hf_geom_10_off + 499 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxx = cbuffer.data(hf_geom_10_off + 500 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxy = cbuffer.data(hf_geom_10_off + 501 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxz = cbuffer.data(hf_geom_10_off + 502 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyy = cbuffer.data(hf_geom_10_off + 503 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyz = cbuffer.data(hf_geom_10_off + 504 * ccomps * dcomps);

            auto g_z_0_xxyzz_xzz = cbuffer.data(hf_geom_10_off + 505 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyy = cbuffer.data(hf_geom_10_off + 506 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyz = cbuffer.data(hf_geom_10_off + 507 * ccomps * dcomps);

            auto g_z_0_xxyzz_yzz = cbuffer.data(hf_geom_10_off + 508 * ccomps * dcomps);

            auto g_z_0_xxyzz_zzz = cbuffer.data(hf_geom_10_off + 509 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxx = cbuffer.data(hf_geom_10_off + 510 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxy = cbuffer.data(hf_geom_10_off + 511 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxz = cbuffer.data(hf_geom_10_off + 512 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyy = cbuffer.data(hf_geom_10_off + 513 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyz = cbuffer.data(hf_geom_10_off + 514 * ccomps * dcomps);

            auto g_z_0_xxzzz_xzz = cbuffer.data(hf_geom_10_off + 515 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyy = cbuffer.data(hf_geom_10_off + 516 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyz = cbuffer.data(hf_geom_10_off + 517 * ccomps * dcomps);

            auto g_z_0_xxzzz_yzz = cbuffer.data(hf_geom_10_off + 518 * ccomps * dcomps);

            auto g_z_0_xxzzz_zzz = cbuffer.data(hf_geom_10_off + 519 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxx = cbuffer.data(hf_geom_10_off + 520 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxy = cbuffer.data(hf_geom_10_off + 521 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxz = cbuffer.data(hf_geom_10_off + 522 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyy = cbuffer.data(hf_geom_10_off + 523 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyz = cbuffer.data(hf_geom_10_off + 524 * ccomps * dcomps);

            auto g_z_0_xyyyy_xzz = cbuffer.data(hf_geom_10_off + 525 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyy = cbuffer.data(hf_geom_10_off + 526 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyz = cbuffer.data(hf_geom_10_off + 527 * ccomps * dcomps);

            auto g_z_0_xyyyy_yzz = cbuffer.data(hf_geom_10_off + 528 * ccomps * dcomps);

            auto g_z_0_xyyyy_zzz = cbuffer.data(hf_geom_10_off + 529 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxx = cbuffer.data(hf_geom_10_off + 530 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxy = cbuffer.data(hf_geom_10_off + 531 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxz = cbuffer.data(hf_geom_10_off + 532 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyy = cbuffer.data(hf_geom_10_off + 533 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyz = cbuffer.data(hf_geom_10_off + 534 * ccomps * dcomps);

            auto g_z_0_xyyyz_xzz = cbuffer.data(hf_geom_10_off + 535 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyy = cbuffer.data(hf_geom_10_off + 536 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyz = cbuffer.data(hf_geom_10_off + 537 * ccomps * dcomps);

            auto g_z_0_xyyyz_yzz = cbuffer.data(hf_geom_10_off + 538 * ccomps * dcomps);

            auto g_z_0_xyyyz_zzz = cbuffer.data(hf_geom_10_off + 539 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxx = cbuffer.data(hf_geom_10_off + 540 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxy = cbuffer.data(hf_geom_10_off + 541 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxz = cbuffer.data(hf_geom_10_off + 542 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyy = cbuffer.data(hf_geom_10_off + 543 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyz = cbuffer.data(hf_geom_10_off + 544 * ccomps * dcomps);

            auto g_z_0_xyyzz_xzz = cbuffer.data(hf_geom_10_off + 545 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyy = cbuffer.data(hf_geom_10_off + 546 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyz = cbuffer.data(hf_geom_10_off + 547 * ccomps * dcomps);

            auto g_z_0_xyyzz_yzz = cbuffer.data(hf_geom_10_off + 548 * ccomps * dcomps);

            auto g_z_0_xyyzz_zzz = cbuffer.data(hf_geom_10_off + 549 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxx = cbuffer.data(hf_geom_10_off + 550 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxy = cbuffer.data(hf_geom_10_off + 551 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxz = cbuffer.data(hf_geom_10_off + 552 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyy = cbuffer.data(hf_geom_10_off + 553 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyz = cbuffer.data(hf_geom_10_off + 554 * ccomps * dcomps);

            auto g_z_0_xyzzz_xzz = cbuffer.data(hf_geom_10_off + 555 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyy = cbuffer.data(hf_geom_10_off + 556 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyz = cbuffer.data(hf_geom_10_off + 557 * ccomps * dcomps);

            auto g_z_0_xyzzz_yzz = cbuffer.data(hf_geom_10_off + 558 * ccomps * dcomps);

            auto g_z_0_xyzzz_zzz = cbuffer.data(hf_geom_10_off + 559 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxx = cbuffer.data(hf_geom_10_off + 560 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxy = cbuffer.data(hf_geom_10_off + 561 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxz = cbuffer.data(hf_geom_10_off + 562 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyy = cbuffer.data(hf_geom_10_off + 563 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyz = cbuffer.data(hf_geom_10_off + 564 * ccomps * dcomps);

            auto g_z_0_xzzzz_xzz = cbuffer.data(hf_geom_10_off + 565 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyy = cbuffer.data(hf_geom_10_off + 566 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyz = cbuffer.data(hf_geom_10_off + 567 * ccomps * dcomps);

            auto g_z_0_xzzzz_yzz = cbuffer.data(hf_geom_10_off + 568 * ccomps * dcomps);

            auto g_z_0_xzzzz_zzz = cbuffer.data(hf_geom_10_off + 569 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxx = cbuffer.data(hf_geom_10_off + 570 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxy = cbuffer.data(hf_geom_10_off + 571 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxz = cbuffer.data(hf_geom_10_off + 572 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyy = cbuffer.data(hf_geom_10_off + 573 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyz = cbuffer.data(hf_geom_10_off + 574 * ccomps * dcomps);

            auto g_z_0_yyyyy_xzz = cbuffer.data(hf_geom_10_off + 575 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyy = cbuffer.data(hf_geom_10_off + 576 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyz = cbuffer.data(hf_geom_10_off + 577 * ccomps * dcomps);

            auto g_z_0_yyyyy_yzz = cbuffer.data(hf_geom_10_off + 578 * ccomps * dcomps);

            auto g_z_0_yyyyy_zzz = cbuffer.data(hf_geom_10_off + 579 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxx = cbuffer.data(hf_geom_10_off + 580 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxy = cbuffer.data(hf_geom_10_off + 581 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxz = cbuffer.data(hf_geom_10_off + 582 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyy = cbuffer.data(hf_geom_10_off + 583 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyz = cbuffer.data(hf_geom_10_off + 584 * ccomps * dcomps);

            auto g_z_0_yyyyz_xzz = cbuffer.data(hf_geom_10_off + 585 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyy = cbuffer.data(hf_geom_10_off + 586 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyz = cbuffer.data(hf_geom_10_off + 587 * ccomps * dcomps);

            auto g_z_0_yyyyz_yzz = cbuffer.data(hf_geom_10_off + 588 * ccomps * dcomps);

            auto g_z_0_yyyyz_zzz = cbuffer.data(hf_geom_10_off + 589 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxx = cbuffer.data(hf_geom_10_off + 590 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxy = cbuffer.data(hf_geom_10_off + 591 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxz = cbuffer.data(hf_geom_10_off + 592 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyy = cbuffer.data(hf_geom_10_off + 593 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyz = cbuffer.data(hf_geom_10_off + 594 * ccomps * dcomps);

            auto g_z_0_yyyzz_xzz = cbuffer.data(hf_geom_10_off + 595 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyy = cbuffer.data(hf_geom_10_off + 596 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyz = cbuffer.data(hf_geom_10_off + 597 * ccomps * dcomps);

            auto g_z_0_yyyzz_yzz = cbuffer.data(hf_geom_10_off + 598 * ccomps * dcomps);

            auto g_z_0_yyyzz_zzz = cbuffer.data(hf_geom_10_off + 599 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxx = cbuffer.data(hf_geom_10_off + 600 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxy = cbuffer.data(hf_geom_10_off + 601 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxz = cbuffer.data(hf_geom_10_off + 602 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyy = cbuffer.data(hf_geom_10_off + 603 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyz = cbuffer.data(hf_geom_10_off + 604 * ccomps * dcomps);

            auto g_z_0_yyzzz_xzz = cbuffer.data(hf_geom_10_off + 605 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyy = cbuffer.data(hf_geom_10_off + 606 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyz = cbuffer.data(hf_geom_10_off + 607 * ccomps * dcomps);

            auto g_z_0_yyzzz_yzz = cbuffer.data(hf_geom_10_off + 608 * ccomps * dcomps);

            auto g_z_0_yyzzz_zzz = cbuffer.data(hf_geom_10_off + 609 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxx = cbuffer.data(hf_geom_10_off + 610 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxy = cbuffer.data(hf_geom_10_off + 611 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxz = cbuffer.data(hf_geom_10_off + 612 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyy = cbuffer.data(hf_geom_10_off + 613 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyz = cbuffer.data(hf_geom_10_off + 614 * ccomps * dcomps);

            auto g_z_0_yzzzz_xzz = cbuffer.data(hf_geom_10_off + 615 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyy = cbuffer.data(hf_geom_10_off + 616 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyz = cbuffer.data(hf_geom_10_off + 617 * ccomps * dcomps);

            auto g_z_0_yzzzz_yzz = cbuffer.data(hf_geom_10_off + 618 * ccomps * dcomps);

            auto g_z_0_yzzzz_zzz = cbuffer.data(hf_geom_10_off + 619 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxx = cbuffer.data(hf_geom_10_off + 620 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxy = cbuffer.data(hf_geom_10_off + 621 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxz = cbuffer.data(hf_geom_10_off + 622 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyy = cbuffer.data(hf_geom_10_off + 623 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyz = cbuffer.data(hf_geom_10_off + 624 * ccomps * dcomps);

            auto g_z_0_zzzzz_xzz = cbuffer.data(hf_geom_10_off + 625 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyy = cbuffer.data(hf_geom_10_off + 626 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyz = cbuffer.data(hf_geom_10_off + 627 * ccomps * dcomps);

            auto g_z_0_zzzzz_yzz = cbuffer.data(hf_geom_10_off + 628 * ccomps * dcomps);

            auto g_z_0_zzzzz_zzz = cbuffer.data(hf_geom_10_off + 629 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HGSS

            const auto hg_geom_10_off = idx_geom_10_hgxx + i * dcomps + j;

            auto g_x_0_xxxxx_xxxx = cbuffer.data(hg_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxy = cbuffer.data(hg_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxz = cbuffer.data(hg_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyy = cbuffer.data(hg_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyz = cbuffer.data(hg_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxzz = cbuffer.data(hg_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyy = cbuffer.data(hg_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyz = cbuffer.data(hg_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyzz = cbuffer.data(hg_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxx_xzzz = cbuffer.data(hg_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyy = cbuffer.data(hg_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyz = cbuffer.data(hg_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyzz = cbuffer.data(hg_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxx_yzzz = cbuffer.data(hg_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxx_zzzz = cbuffer.data(hg_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxx = cbuffer.data(hg_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxy = cbuffer.data(hg_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxz = cbuffer.data(hg_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyy = cbuffer.data(hg_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyz = cbuffer.data(hg_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxzz = cbuffer.data(hg_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyy = cbuffer.data(hg_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyz = cbuffer.data(hg_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyzz = cbuffer.data(hg_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxy_xzzz = cbuffer.data(hg_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyy = cbuffer.data(hg_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyz = cbuffer.data(hg_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyzz = cbuffer.data(hg_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxxy_yzzz = cbuffer.data(hg_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxy_zzzz = cbuffer.data(hg_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxx = cbuffer.data(hg_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxy = cbuffer.data(hg_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxz = cbuffer.data(hg_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyy = cbuffer.data(hg_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyz = cbuffer.data(hg_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxzz = cbuffer.data(hg_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyy = cbuffer.data(hg_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyz = cbuffer.data(hg_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyzz = cbuffer.data(hg_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxxz_xzzz = cbuffer.data(hg_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyy = cbuffer.data(hg_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyz = cbuffer.data(hg_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyzz = cbuffer.data(hg_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxxz_yzzz = cbuffer.data(hg_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxxz_zzzz = cbuffer.data(hg_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxx = cbuffer.data(hg_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxy = cbuffer.data(hg_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxz = cbuffer.data(hg_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyy = cbuffer.data(hg_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyz = cbuffer.data(hg_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxzz = cbuffer.data(hg_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyy = cbuffer.data(hg_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyz = cbuffer.data(hg_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyzz = cbuffer.data(hg_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxyy_xzzz = cbuffer.data(hg_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyy = cbuffer.data(hg_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyz = cbuffer.data(hg_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyzz = cbuffer.data(hg_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxyy_yzzz = cbuffer.data(hg_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxyy_zzzz = cbuffer.data(hg_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxx = cbuffer.data(hg_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxy = cbuffer.data(hg_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxz = cbuffer.data(hg_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyy = cbuffer.data(hg_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyz = cbuffer.data(hg_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxzz = cbuffer.data(hg_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyy = cbuffer.data(hg_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyz = cbuffer.data(hg_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyzz = cbuffer.data(hg_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxyz_xzzz = cbuffer.data(hg_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyy = cbuffer.data(hg_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyz = cbuffer.data(hg_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyzz = cbuffer.data(hg_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxyz_yzzz = cbuffer.data(hg_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxyz_zzzz = cbuffer.data(hg_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxx = cbuffer.data(hg_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxy = cbuffer.data(hg_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxz = cbuffer.data(hg_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyy = cbuffer.data(hg_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyz = cbuffer.data(hg_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxzz = cbuffer.data(hg_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyy = cbuffer.data(hg_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyz = cbuffer.data(hg_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyzz = cbuffer.data(hg_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxxzz_xzzz = cbuffer.data(hg_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyy = cbuffer.data(hg_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyz = cbuffer.data(hg_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyzz = cbuffer.data(hg_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxzz_yzzz = cbuffer.data(hg_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxzz_zzzz = cbuffer.data(hg_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxx = cbuffer.data(hg_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxy = cbuffer.data(hg_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxz = cbuffer.data(hg_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyy = cbuffer.data(hg_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyz = cbuffer.data(hg_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxzz = cbuffer.data(hg_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyy = cbuffer.data(hg_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyz = cbuffer.data(hg_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyzz = cbuffer.data(hg_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxyyy_xzzz = cbuffer.data(hg_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyy = cbuffer.data(hg_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyz = cbuffer.data(hg_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyzz = cbuffer.data(hg_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxyyy_yzzz = cbuffer.data(hg_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxyyy_zzzz = cbuffer.data(hg_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxx = cbuffer.data(hg_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxy = cbuffer.data(hg_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxz = cbuffer.data(hg_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyy = cbuffer.data(hg_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyz = cbuffer.data(hg_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxzz = cbuffer.data(hg_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyy = cbuffer.data(hg_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyz = cbuffer.data(hg_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyzz = cbuffer.data(hg_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxyyz_xzzz = cbuffer.data(hg_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyy = cbuffer.data(hg_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyz = cbuffer.data(hg_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyzz = cbuffer.data(hg_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxyyz_yzzz = cbuffer.data(hg_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxyyz_zzzz = cbuffer.data(hg_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxx = cbuffer.data(hg_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxy = cbuffer.data(hg_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxz = cbuffer.data(hg_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyy = cbuffer.data(hg_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyz = cbuffer.data(hg_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxzz = cbuffer.data(hg_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyy = cbuffer.data(hg_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyz = cbuffer.data(hg_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyzz = cbuffer.data(hg_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxyzz_xzzz = cbuffer.data(hg_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyy = cbuffer.data(hg_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyz = cbuffer.data(hg_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyzz = cbuffer.data(hg_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxyzz_yzzz = cbuffer.data(hg_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxyzz_zzzz = cbuffer.data(hg_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxx = cbuffer.data(hg_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxy = cbuffer.data(hg_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxz = cbuffer.data(hg_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyy = cbuffer.data(hg_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyz = cbuffer.data(hg_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxzz = cbuffer.data(hg_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyy = cbuffer.data(hg_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyz = cbuffer.data(hg_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyzz = cbuffer.data(hg_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxzzz_xzzz = cbuffer.data(hg_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyy = cbuffer.data(hg_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyz = cbuffer.data(hg_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyzz = cbuffer.data(hg_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxzzz_yzzz = cbuffer.data(hg_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxzzz_zzzz = cbuffer.data(hg_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxx = cbuffer.data(hg_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxy = cbuffer.data(hg_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxz = cbuffer.data(hg_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyy = cbuffer.data(hg_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyz = cbuffer.data(hg_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxzz = cbuffer.data(hg_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyy = cbuffer.data(hg_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyz = cbuffer.data(hg_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyzz = cbuffer.data(hg_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xyyyy_xzzz = cbuffer.data(hg_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyy = cbuffer.data(hg_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyz = cbuffer.data(hg_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyzz = cbuffer.data(hg_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xyyyy_yzzz = cbuffer.data(hg_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xyyyy_zzzz = cbuffer.data(hg_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxx = cbuffer.data(hg_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxy = cbuffer.data(hg_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxz = cbuffer.data(hg_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyy = cbuffer.data(hg_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyz = cbuffer.data(hg_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxzz = cbuffer.data(hg_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyy = cbuffer.data(hg_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyz = cbuffer.data(hg_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyzz = cbuffer.data(hg_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xyyyz_xzzz = cbuffer.data(hg_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyy = cbuffer.data(hg_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyz = cbuffer.data(hg_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyzz = cbuffer.data(hg_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xyyyz_yzzz = cbuffer.data(hg_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xyyyz_zzzz = cbuffer.data(hg_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxx = cbuffer.data(hg_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxy = cbuffer.data(hg_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxz = cbuffer.data(hg_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyy = cbuffer.data(hg_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyz = cbuffer.data(hg_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxzz = cbuffer.data(hg_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyy = cbuffer.data(hg_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyz = cbuffer.data(hg_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyzz = cbuffer.data(hg_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xyyzz_xzzz = cbuffer.data(hg_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyy = cbuffer.data(hg_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyz = cbuffer.data(hg_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyzz = cbuffer.data(hg_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xyyzz_yzzz = cbuffer.data(hg_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xyyzz_zzzz = cbuffer.data(hg_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxx = cbuffer.data(hg_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxy = cbuffer.data(hg_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxz = cbuffer.data(hg_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyy = cbuffer.data(hg_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyz = cbuffer.data(hg_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxzz = cbuffer.data(hg_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyy = cbuffer.data(hg_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyz = cbuffer.data(hg_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyzz = cbuffer.data(hg_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xyzzz_xzzz = cbuffer.data(hg_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyy = cbuffer.data(hg_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyz = cbuffer.data(hg_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyzz = cbuffer.data(hg_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xyzzz_yzzz = cbuffer.data(hg_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xyzzz_zzzz = cbuffer.data(hg_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxx = cbuffer.data(hg_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxy = cbuffer.data(hg_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxz = cbuffer.data(hg_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyy = cbuffer.data(hg_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyz = cbuffer.data(hg_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxzz = cbuffer.data(hg_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyy = cbuffer.data(hg_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyz = cbuffer.data(hg_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyzz = cbuffer.data(hg_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xzzzz_xzzz = cbuffer.data(hg_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyy = cbuffer.data(hg_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyz = cbuffer.data(hg_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyzz = cbuffer.data(hg_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xzzzz_yzzz = cbuffer.data(hg_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_xzzzz_zzzz = cbuffer.data(hg_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxx = cbuffer.data(hg_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxy = cbuffer.data(hg_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxz = cbuffer.data(hg_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyy = cbuffer.data(hg_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyz = cbuffer.data(hg_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxzz = cbuffer.data(hg_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyy = cbuffer.data(hg_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyz = cbuffer.data(hg_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyzz = cbuffer.data(hg_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_yyyyy_xzzz = cbuffer.data(hg_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyy = cbuffer.data(hg_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyz = cbuffer.data(hg_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyzz = cbuffer.data(hg_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_yyyyy_yzzz = cbuffer.data(hg_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_yyyyy_zzzz = cbuffer.data(hg_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxx = cbuffer.data(hg_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxy = cbuffer.data(hg_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxz = cbuffer.data(hg_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyy = cbuffer.data(hg_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyz = cbuffer.data(hg_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxzz = cbuffer.data(hg_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyy = cbuffer.data(hg_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyz = cbuffer.data(hg_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyzz = cbuffer.data(hg_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_yyyyz_xzzz = cbuffer.data(hg_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyy = cbuffer.data(hg_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyz = cbuffer.data(hg_geom_10_off + 251 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyzz = cbuffer.data(hg_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_yyyyz_yzzz = cbuffer.data(hg_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_yyyyz_zzzz = cbuffer.data(hg_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxx = cbuffer.data(hg_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxy = cbuffer.data(hg_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxz = cbuffer.data(hg_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyy = cbuffer.data(hg_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyz = cbuffer.data(hg_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxzz = cbuffer.data(hg_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyy = cbuffer.data(hg_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyz = cbuffer.data(hg_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyzz = cbuffer.data(hg_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_yyyzz_xzzz = cbuffer.data(hg_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyy = cbuffer.data(hg_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyz = cbuffer.data(hg_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyzz = cbuffer.data(hg_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_yyyzz_yzzz = cbuffer.data(hg_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_yyyzz_zzzz = cbuffer.data(hg_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxx = cbuffer.data(hg_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxy = cbuffer.data(hg_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxz = cbuffer.data(hg_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyy = cbuffer.data(hg_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyz = cbuffer.data(hg_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxzz = cbuffer.data(hg_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyy = cbuffer.data(hg_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyz = cbuffer.data(hg_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyzz = cbuffer.data(hg_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_yyzzz_xzzz = cbuffer.data(hg_geom_10_off + 279 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyy = cbuffer.data(hg_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyz = cbuffer.data(hg_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyzz = cbuffer.data(hg_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_yyzzz_yzzz = cbuffer.data(hg_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_yyzzz_zzzz = cbuffer.data(hg_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxx = cbuffer.data(hg_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxy = cbuffer.data(hg_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxz = cbuffer.data(hg_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyy = cbuffer.data(hg_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyz = cbuffer.data(hg_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxzz = cbuffer.data(hg_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyy = cbuffer.data(hg_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyz = cbuffer.data(hg_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyzz = cbuffer.data(hg_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_yzzzz_xzzz = cbuffer.data(hg_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyy = cbuffer.data(hg_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyz = cbuffer.data(hg_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyzz = cbuffer.data(hg_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_yzzzz_yzzz = cbuffer.data(hg_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_yzzzz_zzzz = cbuffer.data(hg_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxx = cbuffer.data(hg_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxy = cbuffer.data(hg_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxz = cbuffer.data(hg_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyy = cbuffer.data(hg_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyz = cbuffer.data(hg_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxzz = cbuffer.data(hg_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyy = cbuffer.data(hg_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyz = cbuffer.data(hg_geom_10_off + 307 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyzz = cbuffer.data(hg_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_zzzzz_xzzz = cbuffer.data(hg_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyy = cbuffer.data(hg_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyz = cbuffer.data(hg_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyzz = cbuffer.data(hg_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_zzzzz_yzzz = cbuffer.data(hg_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_zzzzz_zzzz = cbuffer.data(hg_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxx = cbuffer.data(hg_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxy = cbuffer.data(hg_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxz = cbuffer.data(hg_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyy = cbuffer.data(hg_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyz = cbuffer.data(hg_geom_10_off + 319 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxzz = cbuffer.data(hg_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyy = cbuffer.data(hg_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyz = cbuffer.data(hg_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyzz = cbuffer.data(hg_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xxxxx_xzzz = cbuffer.data(hg_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyy = cbuffer.data(hg_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyz = cbuffer.data(hg_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyzz = cbuffer.data(hg_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xxxxx_yzzz = cbuffer.data(hg_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xxxxx_zzzz = cbuffer.data(hg_geom_10_off + 329 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxx = cbuffer.data(hg_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxy = cbuffer.data(hg_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxz = cbuffer.data(hg_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyy = cbuffer.data(hg_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyz = cbuffer.data(hg_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxzz = cbuffer.data(hg_geom_10_off + 335 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyy = cbuffer.data(hg_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyz = cbuffer.data(hg_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyzz = cbuffer.data(hg_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_xxxxy_xzzz = cbuffer.data(hg_geom_10_off + 339 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyy = cbuffer.data(hg_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyz = cbuffer.data(hg_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyzz = cbuffer.data(hg_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_xxxxy_yzzz = cbuffer.data(hg_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_xxxxy_zzzz = cbuffer.data(hg_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxx = cbuffer.data(hg_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxy = cbuffer.data(hg_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxz = cbuffer.data(hg_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyy = cbuffer.data(hg_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyz = cbuffer.data(hg_geom_10_off + 349 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxzz = cbuffer.data(hg_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyy = cbuffer.data(hg_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyz = cbuffer.data(hg_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyzz = cbuffer.data(hg_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_xxxxz_xzzz = cbuffer.data(hg_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyy = cbuffer.data(hg_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyz = cbuffer.data(hg_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyzz = cbuffer.data(hg_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_xxxxz_yzzz = cbuffer.data(hg_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_xxxxz_zzzz = cbuffer.data(hg_geom_10_off + 359 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxx = cbuffer.data(hg_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxy = cbuffer.data(hg_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxz = cbuffer.data(hg_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyy = cbuffer.data(hg_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyz = cbuffer.data(hg_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxzz = cbuffer.data(hg_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyy = cbuffer.data(hg_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyz = cbuffer.data(hg_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyzz = cbuffer.data(hg_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_xxxyy_xzzz = cbuffer.data(hg_geom_10_off + 369 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyy = cbuffer.data(hg_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyz = cbuffer.data(hg_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyzz = cbuffer.data(hg_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_xxxyy_yzzz = cbuffer.data(hg_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_xxxyy_zzzz = cbuffer.data(hg_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxx = cbuffer.data(hg_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxy = cbuffer.data(hg_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxz = cbuffer.data(hg_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyy = cbuffer.data(hg_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyz = cbuffer.data(hg_geom_10_off + 379 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxzz = cbuffer.data(hg_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyy = cbuffer.data(hg_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyz = cbuffer.data(hg_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyzz = cbuffer.data(hg_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_xxxyz_xzzz = cbuffer.data(hg_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyy = cbuffer.data(hg_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyz = cbuffer.data(hg_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyzz = cbuffer.data(hg_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_xxxyz_yzzz = cbuffer.data(hg_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_xxxyz_zzzz = cbuffer.data(hg_geom_10_off + 389 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxx = cbuffer.data(hg_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxy = cbuffer.data(hg_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxz = cbuffer.data(hg_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyy = cbuffer.data(hg_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyz = cbuffer.data(hg_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxzz = cbuffer.data(hg_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyy = cbuffer.data(hg_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyz = cbuffer.data(hg_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyzz = cbuffer.data(hg_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_xxxzz_xzzz = cbuffer.data(hg_geom_10_off + 399 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyy = cbuffer.data(hg_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyz = cbuffer.data(hg_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyzz = cbuffer.data(hg_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_xxxzz_yzzz = cbuffer.data(hg_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_xxxzz_zzzz = cbuffer.data(hg_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxx = cbuffer.data(hg_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxy = cbuffer.data(hg_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxz = cbuffer.data(hg_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyy = cbuffer.data(hg_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyz = cbuffer.data(hg_geom_10_off + 409 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxzz = cbuffer.data(hg_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyy = cbuffer.data(hg_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyz = cbuffer.data(hg_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyzz = cbuffer.data(hg_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_xxyyy_xzzz = cbuffer.data(hg_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyy = cbuffer.data(hg_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyz = cbuffer.data(hg_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyzz = cbuffer.data(hg_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_xxyyy_yzzz = cbuffer.data(hg_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_xxyyy_zzzz = cbuffer.data(hg_geom_10_off + 419 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxx = cbuffer.data(hg_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxy = cbuffer.data(hg_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxz = cbuffer.data(hg_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyy = cbuffer.data(hg_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyz = cbuffer.data(hg_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxzz = cbuffer.data(hg_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyy = cbuffer.data(hg_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyz = cbuffer.data(hg_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyzz = cbuffer.data(hg_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_xxyyz_xzzz = cbuffer.data(hg_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyy = cbuffer.data(hg_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyz = cbuffer.data(hg_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyzz = cbuffer.data(hg_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_xxyyz_yzzz = cbuffer.data(hg_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_xxyyz_zzzz = cbuffer.data(hg_geom_10_off + 434 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxx = cbuffer.data(hg_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxy = cbuffer.data(hg_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxz = cbuffer.data(hg_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyy = cbuffer.data(hg_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyz = cbuffer.data(hg_geom_10_off + 439 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxzz = cbuffer.data(hg_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyy = cbuffer.data(hg_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyz = cbuffer.data(hg_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyzz = cbuffer.data(hg_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xxyzz_xzzz = cbuffer.data(hg_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyy = cbuffer.data(hg_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyz = cbuffer.data(hg_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyzz = cbuffer.data(hg_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_xxyzz_yzzz = cbuffer.data(hg_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_xxyzz_zzzz = cbuffer.data(hg_geom_10_off + 449 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxx = cbuffer.data(hg_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxy = cbuffer.data(hg_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxz = cbuffer.data(hg_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyy = cbuffer.data(hg_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyz = cbuffer.data(hg_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxzz = cbuffer.data(hg_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyy = cbuffer.data(hg_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyz = cbuffer.data(hg_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyzz = cbuffer.data(hg_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_xxzzz_xzzz = cbuffer.data(hg_geom_10_off + 459 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyy = cbuffer.data(hg_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyz = cbuffer.data(hg_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyzz = cbuffer.data(hg_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_xxzzz_yzzz = cbuffer.data(hg_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_xxzzz_zzzz = cbuffer.data(hg_geom_10_off + 464 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxx = cbuffer.data(hg_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxy = cbuffer.data(hg_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxz = cbuffer.data(hg_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyy = cbuffer.data(hg_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyz = cbuffer.data(hg_geom_10_off + 469 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxzz = cbuffer.data(hg_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyy = cbuffer.data(hg_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyz = cbuffer.data(hg_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyzz = cbuffer.data(hg_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_xyyyy_xzzz = cbuffer.data(hg_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyy = cbuffer.data(hg_geom_10_off + 475 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyz = cbuffer.data(hg_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyzz = cbuffer.data(hg_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_xyyyy_yzzz = cbuffer.data(hg_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_xyyyy_zzzz = cbuffer.data(hg_geom_10_off + 479 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxx = cbuffer.data(hg_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxy = cbuffer.data(hg_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxz = cbuffer.data(hg_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyy = cbuffer.data(hg_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyz = cbuffer.data(hg_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxzz = cbuffer.data(hg_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyy = cbuffer.data(hg_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyz = cbuffer.data(hg_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyzz = cbuffer.data(hg_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_xyyyz_xzzz = cbuffer.data(hg_geom_10_off + 489 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyy = cbuffer.data(hg_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyz = cbuffer.data(hg_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyzz = cbuffer.data(hg_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_xyyyz_yzzz = cbuffer.data(hg_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_xyyyz_zzzz = cbuffer.data(hg_geom_10_off + 494 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxx = cbuffer.data(hg_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxy = cbuffer.data(hg_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxz = cbuffer.data(hg_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyy = cbuffer.data(hg_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyz = cbuffer.data(hg_geom_10_off + 499 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxzz = cbuffer.data(hg_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyy = cbuffer.data(hg_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyz = cbuffer.data(hg_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyzz = cbuffer.data(hg_geom_10_off + 503 * ccomps * dcomps);

            auto g_y_0_xyyzz_xzzz = cbuffer.data(hg_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyy = cbuffer.data(hg_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyz = cbuffer.data(hg_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyzz = cbuffer.data(hg_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_xyyzz_yzzz = cbuffer.data(hg_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_xyyzz_zzzz = cbuffer.data(hg_geom_10_off + 509 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxx = cbuffer.data(hg_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxy = cbuffer.data(hg_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxz = cbuffer.data(hg_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyy = cbuffer.data(hg_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyz = cbuffer.data(hg_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxzz = cbuffer.data(hg_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyy = cbuffer.data(hg_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyz = cbuffer.data(hg_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyzz = cbuffer.data(hg_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_xyzzz_xzzz = cbuffer.data(hg_geom_10_off + 519 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyy = cbuffer.data(hg_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyz = cbuffer.data(hg_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyzz = cbuffer.data(hg_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_xyzzz_yzzz = cbuffer.data(hg_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_xyzzz_zzzz = cbuffer.data(hg_geom_10_off + 524 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxx = cbuffer.data(hg_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxy = cbuffer.data(hg_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxz = cbuffer.data(hg_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyy = cbuffer.data(hg_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyz = cbuffer.data(hg_geom_10_off + 529 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxzz = cbuffer.data(hg_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyy = cbuffer.data(hg_geom_10_off + 531 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyz = cbuffer.data(hg_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyzz = cbuffer.data(hg_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_xzzzz_xzzz = cbuffer.data(hg_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyy = cbuffer.data(hg_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyz = cbuffer.data(hg_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyzz = cbuffer.data(hg_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_xzzzz_yzzz = cbuffer.data(hg_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_xzzzz_zzzz = cbuffer.data(hg_geom_10_off + 539 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxx = cbuffer.data(hg_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxy = cbuffer.data(hg_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxz = cbuffer.data(hg_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyy = cbuffer.data(hg_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyz = cbuffer.data(hg_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxzz = cbuffer.data(hg_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyy = cbuffer.data(hg_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyz = cbuffer.data(hg_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyzz = cbuffer.data(hg_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_yyyyy_xzzz = cbuffer.data(hg_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyy = cbuffer.data(hg_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyz = cbuffer.data(hg_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyzz = cbuffer.data(hg_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_yyyyy_yzzz = cbuffer.data(hg_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_yyyyy_zzzz = cbuffer.data(hg_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxx = cbuffer.data(hg_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxy = cbuffer.data(hg_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxz = cbuffer.data(hg_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyy = cbuffer.data(hg_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyz = cbuffer.data(hg_geom_10_off + 559 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxzz = cbuffer.data(hg_geom_10_off + 560 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyy = cbuffer.data(hg_geom_10_off + 561 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyz = cbuffer.data(hg_geom_10_off + 562 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyzz = cbuffer.data(hg_geom_10_off + 563 * ccomps * dcomps);

            auto g_y_0_yyyyz_xzzz = cbuffer.data(hg_geom_10_off + 564 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyy = cbuffer.data(hg_geom_10_off + 565 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyz = cbuffer.data(hg_geom_10_off + 566 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyzz = cbuffer.data(hg_geom_10_off + 567 * ccomps * dcomps);

            auto g_y_0_yyyyz_yzzz = cbuffer.data(hg_geom_10_off + 568 * ccomps * dcomps);

            auto g_y_0_yyyyz_zzzz = cbuffer.data(hg_geom_10_off + 569 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxx = cbuffer.data(hg_geom_10_off + 570 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxy = cbuffer.data(hg_geom_10_off + 571 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxz = cbuffer.data(hg_geom_10_off + 572 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyy = cbuffer.data(hg_geom_10_off + 573 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyz = cbuffer.data(hg_geom_10_off + 574 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxzz = cbuffer.data(hg_geom_10_off + 575 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyy = cbuffer.data(hg_geom_10_off + 576 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyz = cbuffer.data(hg_geom_10_off + 577 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyzz = cbuffer.data(hg_geom_10_off + 578 * ccomps * dcomps);

            auto g_y_0_yyyzz_xzzz = cbuffer.data(hg_geom_10_off + 579 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyy = cbuffer.data(hg_geom_10_off + 580 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyz = cbuffer.data(hg_geom_10_off + 581 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyzz = cbuffer.data(hg_geom_10_off + 582 * ccomps * dcomps);

            auto g_y_0_yyyzz_yzzz = cbuffer.data(hg_geom_10_off + 583 * ccomps * dcomps);

            auto g_y_0_yyyzz_zzzz = cbuffer.data(hg_geom_10_off + 584 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxx = cbuffer.data(hg_geom_10_off + 585 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxy = cbuffer.data(hg_geom_10_off + 586 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxz = cbuffer.data(hg_geom_10_off + 587 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyy = cbuffer.data(hg_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyz = cbuffer.data(hg_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxzz = cbuffer.data(hg_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyy = cbuffer.data(hg_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyz = cbuffer.data(hg_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyzz = cbuffer.data(hg_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_yyzzz_xzzz = cbuffer.data(hg_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyy = cbuffer.data(hg_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyz = cbuffer.data(hg_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyzz = cbuffer.data(hg_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_yyzzz_yzzz = cbuffer.data(hg_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_yyzzz_zzzz = cbuffer.data(hg_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxx = cbuffer.data(hg_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxy = cbuffer.data(hg_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxz = cbuffer.data(hg_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyy = cbuffer.data(hg_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyz = cbuffer.data(hg_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxzz = cbuffer.data(hg_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyy = cbuffer.data(hg_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyz = cbuffer.data(hg_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyzz = cbuffer.data(hg_geom_10_off + 608 * ccomps * dcomps);

            auto g_y_0_yzzzz_xzzz = cbuffer.data(hg_geom_10_off + 609 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyy = cbuffer.data(hg_geom_10_off + 610 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyz = cbuffer.data(hg_geom_10_off + 611 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyzz = cbuffer.data(hg_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_yzzzz_yzzz = cbuffer.data(hg_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_yzzzz_zzzz = cbuffer.data(hg_geom_10_off + 614 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxx = cbuffer.data(hg_geom_10_off + 615 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxy = cbuffer.data(hg_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxz = cbuffer.data(hg_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyy = cbuffer.data(hg_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyz = cbuffer.data(hg_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxzz = cbuffer.data(hg_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyy = cbuffer.data(hg_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyz = cbuffer.data(hg_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyzz = cbuffer.data(hg_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_zzzzz_xzzz = cbuffer.data(hg_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyy = cbuffer.data(hg_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyz = cbuffer.data(hg_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyzz = cbuffer.data(hg_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_zzzzz_yzzz = cbuffer.data(hg_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_zzzzz_zzzz = cbuffer.data(hg_geom_10_off + 629 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxx = cbuffer.data(hg_geom_10_off + 630 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxy = cbuffer.data(hg_geom_10_off + 631 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxz = cbuffer.data(hg_geom_10_off + 632 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyy = cbuffer.data(hg_geom_10_off + 633 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyz = cbuffer.data(hg_geom_10_off + 634 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxzz = cbuffer.data(hg_geom_10_off + 635 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyy = cbuffer.data(hg_geom_10_off + 636 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyz = cbuffer.data(hg_geom_10_off + 637 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyzz = cbuffer.data(hg_geom_10_off + 638 * ccomps * dcomps);

            auto g_z_0_xxxxx_xzzz = cbuffer.data(hg_geom_10_off + 639 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyy = cbuffer.data(hg_geom_10_off + 640 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyz = cbuffer.data(hg_geom_10_off + 641 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyzz = cbuffer.data(hg_geom_10_off + 642 * ccomps * dcomps);

            auto g_z_0_xxxxx_yzzz = cbuffer.data(hg_geom_10_off + 643 * ccomps * dcomps);

            auto g_z_0_xxxxx_zzzz = cbuffer.data(hg_geom_10_off + 644 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxx = cbuffer.data(hg_geom_10_off + 645 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxy = cbuffer.data(hg_geom_10_off + 646 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxz = cbuffer.data(hg_geom_10_off + 647 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyy = cbuffer.data(hg_geom_10_off + 648 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyz = cbuffer.data(hg_geom_10_off + 649 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxzz = cbuffer.data(hg_geom_10_off + 650 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyy = cbuffer.data(hg_geom_10_off + 651 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyz = cbuffer.data(hg_geom_10_off + 652 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyzz = cbuffer.data(hg_geom_10_off + 653 * ccomps * dcomps);

            auto g_z_0_xxxxy_xzzz = cbuffer.data(hg_geom_10_off + 654 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyy = cbuffer.data(hg_geom_10_off + 655 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyz = cbuffer.data(hg_geom_10_off + 656 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyzz = cbuffer.data(hg_geom_10_off + 657 * ccomps * dcomps);

            auto g_z_0_xxxxy_yzzz = cbuffer.data(hg_geom_10_off + 658 * ccomps * dcomps);

            auto g_z_0_xxxxy_zzzz = cbuffer.data(hg_geom_10_off + 659 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxx = cbuffer.data(hg_geom_10_off + 660 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxy = cbuffer.data(hg_geom_10_off + 661 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxz = cbuffer.data(hg_geom_10_off + 662 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyy = cbuffer.data(hg_geom_10_off + 663 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyz = cbuffer.data(hg_geom_10_off + 664 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxzz = cbuffer.data(hg_geom_10_off + 665 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyy = cbuffer.data(hg_geom_10_off + 666 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyz = cbuffer.data(hg_geom_10_off + 667 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyzz = cbuffer.data(hg_geom_10_off + 668 * ccomps * dcomps);

            auto g_z_0_xxxxz_xzzz = cbuffer.data(hg_geom_10_off + 669 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyy = cbuffer.data(hg_geom_10_off + 670 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyz = cbuffer.data(hg_geom_10_off + 671 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyzz = cbuffer.data(hg_geom_10_off + 672 * ccomps * dcomps);

            auto g_z_0_xxxxz_yzzz = cbuffer.data(hg_geom_10_off + 673 * ccomps * dcomps);

            auto g_z_0_xxxxz_zzzz = cbuffer.data(hg_geom_10_off + 674 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxx = cbuffer.data(hg_geom_10_off + 675 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxy = cbuffer.data(hg_geom_10_off + 676 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxz = cbuffer.data(hg_geom_10_off + 677 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyy = cbuffer.data(hg_geom_10_off + 678 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyz = cbuffer.data(hg_geom_10_off + 679 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxzz = cbuffer.data(hg_geom_10_off + 680 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyy = cbuffer.data(hg_geom_10_off + 681 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyz = cbuffer.data(hg_geom_10_off + 682 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyzz = cbuffer.data(hg_geom_10_off + 683 * ccomps * dcomps);

            auto g_z_0_xxxyy_xzzz = cbuffer.data(hg_geom_10_off + 684 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyy = cbuffer.data(hg_geom_10_off + 685 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyz = cbuffer.data(hg_geom_10_off + 686 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyzz = cbuffer.data(hg_geom_10_off + 687 * ccomps * dcomps);

            auto g_z_0_xxxyy_yzzz = cbuffer.data(hg_geom_10_off + 688 * ccomps * dcomps);

            auto g_z_0_xxxyy_zzzz = cbuffer.data(hg_geom_10_off + 689 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxx = cbuffer.data(hg_geom_10_off + 690 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxy = cbuffer.data(hg_geom_10_off + 691 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxz = cbuffer.data(hg_geom_10_off + 692 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyy = cbuffer.data(hg_geom_10_off + 693 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyz = cbuffer.data(hg_geom_10_off + 694 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxzz = cbuffer.data(hg_geom_10_off + 695 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyy = cbuffer.data(hg_geom_10_off + 696 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyz = cbuffer.data(hg_geom_10_off + 697 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyzz = cbuffer.data(hg_geom_10_off + 698 * ccomps * dcomps);

            auto g_z_0_xxxyz_xzzz = cbuffer.data(hg_geom_10_off + 699 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyy = cbuffer.data(hg_geom_10_off + 700 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyz = cbuffer.data(hg_geom_10_off + 701 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyzz = cbuffer.data(hg_geom_10_off + 702 * ccomps * dcomps);

            auto g_z_0_xxxyz_yzzz = cbuffer.data(hg_geom_10_off + 703 * ccomps * dcomps);

            auto g_z_0_xxxyz_zzzz = cbuffer.data(hg_geom_10_off + 704 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxx = cbuffer.data(hg_geom_10_off + 705 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxy = cbuffer.data(hg_geom_10_off + 706 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxz = cbuffer.data(hg_geom_10_off + 707 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyy = cbuffer.data(hg_geom_10_off + 708 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyz = cbuffer.data(hg_geom_10_off + 709 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxzz = cbuffer.data(hg_geom_10_off + 710 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyy = cbuffer.data(hg_geom_10_off + 711 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyz = cbuffer.data(hg_geom_10_off + 712 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyzz = cbuffer.data(hg_geom_10_off + 713 * ccomps * dcomps);

            auto g_z_0_xxxzz_xzzz = cbuffer.data(hg_geom_10_off + 714 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyy = cbuffer.data(hg_geom_10_off + 715 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyz = cbuffer.data(hg_geom_10_off + 716 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyzz = cbuffer.data(hg_geom_10_off + 717 * ccomps * dcomps);

            auto g_z_0_xxxzz_yzzz = cbuffer.data(hg_geom_10_off + 718 * ccomps * dcomps);

            auto g_z_0_xxxzz_zzzz = cbuffer.data(hg_geom_10_off + 719 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxx = cbuffer.data(hg_geom_10_off + 720 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxy = cbuffer.data(hg_geom_10_off + 721 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxz = cbuffer.data(hg_geom_10_off + 722 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyy = cbuffer.data(hg_geom_10_off + 723 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyz = cbuffer.data(hg_geom_10_off + 724 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxzz = cbuffer.data(hg_geom_10_off + 725 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyy = cbuffer.data(hg_geom_10_off + 726 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyz = cbuffer.data(hg_geom_10_off + 727 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyzz = cbuffer.data(hg_geom_10_off + 728 * ccomps * dcomps);

            auto g_z_0_xxyyy_xzzz = cbuffer.data(hg_geom_10_off + 729 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyy = cbuffer.data(hg_geom_10_off + 730 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyz = cbuffer.data(hg_geom_10_off + 731 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyzz = cbuffer.data(hg_geom_10_off + 732 * ccomps * dcomps);

            auto g_z_0_xxyyy_yzzz = cbuffer.data(hg_geom_10_off + 733 * ccomps * dcomps);

            auto g_z_0_xxyyy_zzzz = cbuffer.data(hg_geom_10_off + 734 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxx = cbuffer.data(hg_geom_10_off + 735 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxy = cbuffer.data(hg_geom_10_off + 736 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxz = cbuffer.data(hg_geom_10_off + 737 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyy = cbuffer.data(hg_geom_10_off + 738 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyz = cbuffer.data(hg_geom_10_off + 739 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxzz = cbuffer.data(hg_geom_10_off + 740 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyy = cbuffer.data(hg_geom_10_off + 741 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyz = cbuffer.data(hg_geom_10_off + 742 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyzz = cbuffer.data(hg_geom_10_off + 743 * ccomps * dcomps);

            auto g_z_0_xxyyz_xzzz = cbuffer.data(hg_geom_10_off + 744 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyy = cbuffer.data(hg_geom_10_off + 745 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyz = cbuffer.data(hg_geom_10_off + 746 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyzz = cbuffer.data(hg_geom_10_off + 747 * ccomps * dcomps);

            auto g_z_0_xxyyz_yzzz = cbuffer.data(hg_geom_10_off + 748 * ccomps * dcomps);

            auto g_z_0_xxyyz_zzzz = cbuffer.data(hg_geom_10_off + 749 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxx = cbuffer.data(hg_geom_10_off + 750 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxy = cbuffer.data(hg_geom_10_off + 751 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxz = cbuffer.data(hg_geom_10_off + 752 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyy = cbuffer.data(hg_geom_10_off + 753 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyz = cbuffer.data(hg_geom_10_off + 754 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxzz = cbuffer.data(hg_geom_10_off + 755 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyy = cbuffer.data(hg_geom_10_off + 756 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyz = cbuffer.data(hg_geom_10_off + 757 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyzz = cbuffer.data(hg_geom_10_off + 758 * ccomps * dcomps);

            auto g_z_0_xxyzz_xzzz = cbuffer.data(hg_geom_10_off + 759 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyy = cbuffer.data(hg_geom_10_off + 760 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyz = cbuffer.data(hg_geom_10_off + 761 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyzz = cbuffer.data(hg_geom_10_off + 762 * ccomps * dcomps);

            auto g_z_0_xxyzz_yzzz = cbuffer.data(hg_geom_10_off + 763 * ccomps * dcomps);

            auto g_z_0_xxyzz_zzzz = cbuffer.data(hg_geom_10_off + 764 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxx = cbuffer.data(hg_geom_10_off + 765 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxy = cbuffer.data(hg_geom_10_off + 766 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxz = cbuffer.data(hg_geom_10_off + 767 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyy = cbuffer.data(hg_geom_10_off + 768 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyz = cbuffer.data(hg_geom_10_off + 769 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxzz = cbuffer.data(hg_geom_10_off + 770 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyy = cbuffer.data(hg_geom_10_off + 771 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyz = cbuffer.data(hg_geom_10_off + 772 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyzz = cbuffer.data(hg_geom_10_off + 773 * ccomps * dcomps);

            auto g_z_0_xxzzz_xzzz = cbuffer.data(hg_geom_10_off + 774 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyy = cbuffer.data(hg_geom_10_off + 775 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyz = cbuffer.data(hg_geom_10_off + 776 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyzz = cbuffer.data(hg_geom_10_off + 777 * ccomps * dcomps);

            auto g_z_0_xxzzz_yzzz = cbuffer.data(hg_geom_10_off + 778 * ccomps * dcomps);

            auto g_z_0_xxzzz_zzzz = cbuffer.data(hg_geom_10_off + 779 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxx = cbuffer.data(hg_geom_10_off + 780 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxy = cbuffer.data(hg_geom_10_off + 781 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxz = cbuffer.data(hg_geom_10_off + 782 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyy = cbuffer.data(hg_geom_10_off + 783 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyz = cbuffer.data(hg_geom_10_off + 784 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxzz = cbuffer.data(hg_geom_10_off + 785 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyy = cbuffer.data(hg_geom_10_off + 786 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyz = cbuffer.data(hg_geom_10_off + 787 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyzz = cbuffer.data(hg_geom_10_off + 788 * ccomps * dcomps);

            auto g_z_0_xyyyy_xzzz = cbuffer.data(hg_geom_10_off + 789 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyy = cbuffer.data(hg_geom_10_off + 790 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyz = cbuffer.data(hg_geom_10_off + 791 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyzz = cbuffer.data(hg_geom_10_off + 792 * ccomps * dcomps);

            auto g_z_0_xyyyy_yzzz = cbuffer.data(hg_geom_10_off + 793 * ccomps * dcomps);

            auto g_z_0_xyyyy_zzzz = cbuffer.data(hg_geom_10_off + 794 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxx = cbuffer.data(hg_geom_10_off + 795 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxy = cbuffer.data(hg_geom_10_off + 796 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxz = cbuffer.data(hg_geom_10_off + 797 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyy = cbuffer.data(hg_geom_10_off + 798 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyz = cbuffer.data(hg_geom_10_off + 799 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxzz = cbuffer.data(hg_geom_10_off + 800 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyy = cbuffer.data(hg_geom_10_off + 801 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyz = cbuffer.data(hg_geom_10_off + 802 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyzz = cbuffer.data(hg_geom_10_off + 803 * ccomps * dcomps);

            auto g_z_0_xyyyz_xzzz = cbuffer.data(hg_geom_10_off + 804 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyy = cbuffer.data(hg_geom_10_off + 805 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyz = cbuffer.data(hg_geom_10_off + 806 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyzz = cbuffer.data(hg_geom_10_off + 807 * ccomps * dcomps);

            auto g_z_0_xyyyz_yzzz = cbuffer.data(hg_geom_10_off + 808 * ccomps * dcomps);

            auto g_z_0_xyyyz_zzzz = cbuffer.data(hg_geom_10_off + 809 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxx = cbuffer.data(hg_geom_10_off + 810 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxy = cbuffer.data(hg_geom_10_off + 811 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxz = cbuffer.data(hg_geom_10_off + 812 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyy = cbuffer.data(hg_geom_10_off + 813 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyz = cbuffer.data(hg_geom_10_off + 814 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxzz = cbuffer.data(hg_geom_10_off + 815 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyy = cbuffer.data(hg_geom_10_off + 816 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyz = cbuffer.data(hg_geom_10_off + 817 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyzz = cbuffer.data(hg_geom_10_off + 818 * ccomps * dcomps);

            auto g_z_0_xyyzz_xzzz = cbuffer.data(hg_geom_10_off + 819 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyy = cbuffer.data(hg_geom_10_off + 820 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyz = cbuffer.data(hg_geom_10_off + 821 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyzz = cbuffer.data(hg_geom_10_off + 822 * ccomps * dcomps);

            auto g_z_0_xyyzz_yzzz = cbuffer.data(hg_geom_10_off + 823 * ccomps * dcomps);

            auto g_z_0_xyyzz_zzzz = cbuffer.data(hg_geom_10_off + 824 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxx = cbuffer.data(hg_geom_10_off + 825 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxy = cbuffer.data(hg_geom_10_off + 826 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxz = cbuffer.data(hg_geom_10_off + 827 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyy = cbuffer.data(hg_geom_10_off + 828 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyz = cbuffer.data(hg_geom_10_off + 829 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxzz = cbuffer.data(hg_geom_10_off + 830 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyy = cbuffer.data(hg_geom_10_off + 831 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyz = cbuffer.data(hg_geom_10_off + 832 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyzz = cbuffer.data(hg_geom_10_off + 833 * ccomps * dcomps);

            auto g_z_0_xyzzz_xzzz = cbuffer.data(hg_geom_10_off + 834 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyy = cbuffer.data(hg_geom_10_off + 835 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyz = cbuffer.data(hg_geom_10_off + 836 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyzz = cbuffer.data(hg_geom_10_off + 837 * ccomps * dcomps);

            auto g_z_0_xyzzz_yzzz = cbuffer.data(hg_geom_10_off + 838 * ccomps * dcomps);

            auto g_z_0_xyzzz_zzzz = cbuffer.data(hg_geom_10_off + 839 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxx = cbuffer.data(hg_geom_10_off + 840 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxy = cbuffer.data(hg_geom_10_off + 841 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxz = cbuffer.data(hg_geom_10_off + 842 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyy = cbuffer.data(hg_geom_10_off + 843 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyz = cbuffer.data(hg_geom_10_off + 844 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxzz = cbuffer.data(hg_geom_10_off + 845 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyy = cbuffer.data(hg_geom_10_off + 846 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyz = cbuffer.data(hg_geom_10_off + 847 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyzz = cbuffer.data(hg_geom_10_off + 848 * ccomps * dcomps);

            auto g_z_0_xzzzz_xzzz = cbuffer.data(hg_geom_10_off + 849 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyy = cbuffer.data(hg_geom_10_off + 850 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyz = cbuffer.data(hg_geom_10_off + 851 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyzz = cbuffer.data(hg_geom_10_off + 852 * ccomps * dcomps);

            auto g_z_0_xzzzz_yzzz = cbuffer.data(hg_geom_10_off + 853 * ccomps * dcomps);

            auto g_z_0_xzzzz_zzzz = cbuffer.data(hg_geom_10_off + 854 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxx = cbuffer.data(hg_geom_10_off + 855 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxy = cbuffer.data(hg_geom_10_off + 856 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxz = cbuffer.data(hg_geom_10_off + 857 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyy = cbuffer.data(hg_geom_10_off + 858 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyz = cbuffer.data(hg_geom_10_off + 859 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxzz = cbuffer.data(hg_geom_10_off + 860 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyy = cbuffer.data(hg_geom_10_off + 861 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyz = cbuffer.data(hg_geom_10_off + 862 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyzz = cbuffer.data(hg_geom_10_off + 863 * ccomps * dcomps);

            auto g_z_0_yyyyy_xzzz = cbuffer.data(hg_geom_10_off + 864 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyy = cbuffer.data(hg_geom_10_off + 865 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyz = cbuffer.data(hg_geom_10_off + 866 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyzz = cbuffer.data(hg_geom_10_off + 867 * ccomps * dcomps);

            auto g_z_0_yyyyy_yzzz = cbuffer.data(hg_geom_10_off + 868 * ccomps * dcomps);

            auto g_z_0_yyyyy_zzzz = cbuffer.data(hg_geom_10_off + 869 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxx = cbuffer.data(hg_geom_10_off + 870 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxy = cbuffer.data(hg_geom_10_off + 871 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxz = cbuffer.data(hg_geom_10_off + 872 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyy = cbuffer.data(hg_geom_10_off + 873 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyz = cbuffer.data(hg_geom_10_off + 874 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxzz = cbuffer.data(hg_geom_10_off + 875 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyy = cbuffer.data(hg_geom_10_off + 876 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyz = cbuffer.data(hg_geom_10_off + 877 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyzz = cbuffer.data(hg_geom_10_off + 878 * ccomps * dcomps);

            auto g_z_0_yyyyz_xzzz = cbuffer.data(hg_geom_10_off + 879 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyy = cbuffer.data(hg_geom_10_off + 880 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyz = cbuffer.data(hg_geom_10_off + 881 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyzz = cbuffer.data(hg_geom_10_off + 882 * ccomps * dcomps);

            auto g_z_0_yyyyz_yzzz = cbuffer.data(hg_geom_10_off + 883 * ccomps * dcomps);

            auto g_z_0_yyyyz_zzzz = cbuffer.data(hg_geom_10_off + 884 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxx = cbuffer.data(hg_geom_10_off + 885 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxy = cbuffer.data(hg_geom_10_off + 886 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxz = cbuffer.data(hg_geom_10_off + 887 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyy = cbuffer.data(hg_geom_10_off + 888 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyz = cbuffer.data(hg_geom_10_off + 889 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxzz = cbuffer.data(hg_geom_10_off + 890 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyy = cbuffer.data(hg_geom_10_off + 891 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyz = cbuffer.data(hg_geom_10_off + 892 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyzz = cbuffer.data(hg_geom_10_off + 893 * ccomps * dcomps);

            auto g_z_0_yyyzz_xzzz = cbuffer.data(hg_geom_10_off + 894 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyy = cbuffer.data(hg_geom_10_off + 895 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyz = cbuffer.data(hg_geom_10_off + 896 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyzz = cbuffer.data(hg_geom_10_off + 897 * ccomps * dcomps);

            auto g_z_0_yyyzz_yzzz = cbuffer.data(hg_geom_10_off + 898 * ccomps * dcomps);

            auto g_z_0_yyyzz_zzzz = cbuffer.data(hg_geom_10_off + 899 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxx = cbuffer.data(hg_geom_10_off + 900 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxy = cbuffer.data(hg_geom_10_off + 901 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxz = cbuffer.data(hg_geom_10_off + 902 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyy = cbuffer.data(hg_geom_10_off + 903 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyz = cbuffer.data(hg_geom_10_off + 904 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxzz = cbuffer.data(hg_geom_10_off + 905 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyy = cbuffer.data(hg_geom_10_off + 906 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyz = cbuffer.data(hg_geom_10_off + 907 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyzz = cbuffer.data(hg_geom_10_off + 908 * ccomps * dcomps);

            auto g_z_0_yyzzz_xzzz = cbuffer.data(hg_geom_10_off + 909 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyy = cbuffer.data(hg_geom_10_off + 910 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyz = cbuffer.data(hg_geom_10_off + 911 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyzz = cbuffer.data(hg_geom_10_off + 912 * ccomps * dcomps);

            auto g_z_0_yyzzz_yzzz = cbuffer.data(hg_geom_10_off + 913 * ccomps * dcomps);

            auto g_z_0_yyzzz_zzzz = cbuffer.data(hg_geom_10_off + 914 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxx = cbuffer.data(hg_geom_10_off + 915 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxy = cbuffer.data(hg_geom_10_off + 916 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxz = cbuffer.data(hg_geom_10_off + 917 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyy = cbuffer.data(hg_geom_10_off + 918 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyz = cbuffer.data(hg_geom_10_off + 919 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxzz = cbuffer.data(hg_geom_10_off + 920 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyy = cbuffer.data(hg_geom_10_off + 921 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyz = cbuffer.data(hg_geom_10_off + 922 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyzz = cbuffer.data(hg_geom_10_off + 923 * ccomps * dcomps);

            auto g_z_0_yzzzz_xzzz = cbuffer.data(hg_geom_10_off + 924 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyy = cbuffer.data(hg_geom_10_off + 925 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyz = cbuffer.data(hg_geom_10_off + 926 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyzz = cbuffer.data(hg_geom_10_off + 927 * ccomps * dcomps);

            auto g_z_0_yzzzz_yzzz = cbuffer.data(hg_geom_10_off + 928 * ccomps * dcomps);

            auto g_z_0_yzzzz_zzzz = cbuffer.data(hg_geom_10_off + 929 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxx = cbuffer.data(hg_geom_10_off + 930 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxy = cbuffer.data(hg_geom_10_off + 931 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxz = cbuffer.data(hg_geom_10_off + 932 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyy = cbuffer.data(hg_geom_10_off + 933 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyz = cbuffer.data(hg_geom_10_off + 934 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxzz = cbuffer.data(hg_geom_10_off + 935 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyy = cbuffer.data(hg_geom_10_off + 936 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyz = cbuffer.data(hg_geom_10_off + 937 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyzz = cbuffer.data(hg_geom_10_off + 938 * ccomps * dcomps);

            auto g_z_0_zzzzz_xzzz = cbuffer.data(hg_geom_10_off + 939 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyy = cbuffer.data(hg_geom_10_off + 940 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyz = cbuffer.data(hg_geom_10_off + 941 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyzz = cbuffer.data(hg_geom_10_off + 942 * ccomps * dcomps);

            auto g_z_0_zzzzz_yzzz = cbuffer.data(hg_geom_10_off + 943 * ccomps * dcomps);

            auto g_z_0_zzzzz_zzzz = cbuffer.data(hg_geom_10_off + 944 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_ifxx

            const auto if_geom_10_off = idx_geom_10_ifxx + i * dcomps + j;

            /// Set up 0-10 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_xxx = cbuffer.data(if_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxy = cbuffer.data(if_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxz = cbuffer.data(if_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyy = cbuffer.data(if_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyz = cbuffer.data(if_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xzz = cbuffer.data(if_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyy = cbuffer.data(if_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyz = cbuffer.data(if_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yzz = cbuffer.data(if_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxxx_zzz = cbuffer.data(if_geom_10_off + 9 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxx, g_x_0_xxxxx_xxxx, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxy, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xyy, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_yyy, g_x_0_xxxxx_yyz, g_x_0_xxxxx_yzz, g_x_0_xxxxx_zzz, g_x_0_xxxxxx_xxx, g_x_0_xxxxxx_xxy, g_x_0_xxxxxx_xxz, g_x_0_xxxxxx_xyy, g_x_0_xxxxxx_xyz, g_x_0_xxxxxx_xzz, g_x_0_xxxxxx_yyy, g_x_0_xxxxxx_yyz, g_x_0_xxxxxx_yzz, g_x_0_xxxxxx_zzz, g_xxxxx_xxx, g_xxxxx_xxy, g_xxxxx_xxz, g_xxxxx_xyy, g_xxxxx_xyz, g_xxxxx_xzz, g_xxxxx_yyy, g_xxxxx_yyz, g_xxxxx_yzz, g_xxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_xxx[k] = -g_xxxxx_xxx[k] - g_x_0_xxxxx_xxx[k] * ab_x + g_x_0_xxxxx_xxxx[k];

                g_x_0_xxxxxx_xxy[k] = -g_xxxxx_xxy[k] - g_x_0_xxxxx_xxy[k] * ab_x + g_x_0_xxxxx_xxxy[k];

                g_x_0_xxxxxx_xxz[k] = -g_xxxxx_xxz[k] - g_x_0_xxxxx_xxz[k] * ab_x + g_x_0_xxxxx_xxxz[k];

                g_x_0_xxxxxx_xyy[k] = -g_xxxxx_xyy[k] - g_x_0_xxxxx_xyy[k] * ab_x + g_x_0_xxxxx_xxyy[k];

                g_x_0_xxxxxx_xyz[k] = -g_xxxxx_xyz[k] - g_x_0_xxxxx_xyz[k] * ab_x + g_x_0_xxxxx_xxyz[k];

                g_x_0_xxxxxx_xzz[k] = -g_xxxxx_xzz[k] - g_x_0_xxxxx_xzz[k] * ab_x + g_x_0_xxxxx_xxzz[k];

                g_x_0_xxxxxx_yyy[k] = -g_xxxxx_yyy[k] - g_x_0_xxxxx_yyy[k] * ab_x + g_x_0_xxxxx_xyyy[k];

                g_x_0_xxxxxx_yyz[k] = -g_xxxxx_yyz[k] - g_x_0_xxxxx_yyz[k] * ab_x + g_x_0_xxxxx_xyyz[k];

                g_x_0_xxxxxx_yzz[k] = -g_xxxxx_yzz[k] - g_x_0_xxxxx_yzz[k] * ab_x + g_x_0_xxxxx_xyzz[k];

                g_x_0_xxxxxx_zzz[k] = -g_xxxxx_zzz[k] - g_x_0_xxxxx_zzz[k] * ab_x + g_x_0_xxxxx_xzzz[k];
            }

            /// Set up 10-20 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_xxx = cbuffer.data(if_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxy = cbuffer.data(if_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxz = cbuffer.data(if_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyy = cbuffer.data(if_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyz = cbuffer.data(if_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xzz = cbuffer.data(if_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyy = cbuffer.data(if_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyz = cbuffer.data(if_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yzz = cbuffer.data(if_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxxy_zzz = cbuffer.data(if_geom_10_off + 19 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxx, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxy, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxz, g_x_0_xxxxx_xyy, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xzz, g_x_0_xxxxx_yyy, g_x_0_xxxxx_yyyy, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_zzz, g_x_0_xxxxxy_xxx, g_x_0_xxxxxy_xxy, g_x_0_xxxxxy_xxz, g_x_0_xxxxxy_xyy, g_x_0_xxxxxy_xyz, g_x_0_xxxxxy_xzz, g_x_0_xxxxxy_yyy, g_x_0_xxxxxy_yyz, g_x_0_xxxxxy_yzz, g_x_0_xxxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_xxx[k] = -g_x_0_xxxxx_xxx[k] * ab_y + g_x_0_xxxxx_xxxy[k];

                g_x_0_xxxxxy_xxy[k] = -g_x_0_xxxxx_xxy[k] * ab_y + g_x_0_xxxxx_xxyy[k];

                g_x_0_xxxxxy_xxz[k] = -g_x_0_xxxxx_xxz[k] * ab_y + g_x_0_xxxxx_xxyz[k];

                g_x_0_xxxxxy_xyy[k] = -g_x_0_xxxxx_xyy[k] * ab_y + g_x_0_xxxxx_xyyy[k];

                g_x_0_xxxxxy_xyz[k] = -g_x_0_xxxxx_xyz[k] * ab_y + g_x_0_xxxxx_xyyz[k];

                g_x_0_xxxxxy_xzz[k] = -g_x_0_xxxxx_xzz[k] * ab_y + g_x_0_xxxxx_xyzz[k];

                g_x_0_xxxxxy_yyy[k] = -g_x_0_xxxxx_yyy[k] * ab_y + g_x_0_xxxxx_yyyy[k];

                g_x_0_xxxxxy_yyz[k] = -g_x_0_xxxxx_yyz[k] * ab_y + g_x_0_xxxxx_yyyz[k];

                g_x_0_xxxxxy_yzz[k] = -g_x_0_xxxxx_yzz[k] * ab_y + g_x_0_xxxxx_yyzz[k];

                g_x_0_xxxxxy_zzz[k] = -g_x_0_xxxxx_zzz[k] * ab_y + g_x_0_xxxxx_yzzz[k];
            }

            /// Set up 20-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_xxx = cbuffer.data(if_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxy = cbuffer.data(if_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxz = cbuffer.data(if_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyy = cbuffer.data(if_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyz = cbuffer.data(if_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xzz = cbuffer.data(if_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyy = cbuffer.data(if_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyz = cbuffer.data(if_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yzz = cbuffer.data(if_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxxz_zzz = cbuffer.data(if_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxx, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxy, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xyy, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_yyy, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_zzz, g_x_0_xxxxx_zzzz, g_x_0_xxxxxz_xxx, g_x_0_xxxxxz_xxy, g_x_0_xxxxxz_xxz, g_x_0_xxxxxz_xyy, g_x_0_xxxxxz_xyz, g_x_0_xxxxxz_xzz, g_x_0_xxxxxz_yyy, g_x_0_xxxxxz_yyz, g_x_0_xxxxxz_yzz, g_x_0_xxxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_xxx[k] = -g_x_0_xxxxx_xxx[k] * ab_z + g_x_0_xxxxx_xxxz[k];

                g_x_0_xxxxxz_xxy[k] = -g_x_0_xxxxx_xxy[k] * ab_z + g_x_0_xxxxx_xxyz[k];

                g_x_0_xxxxxz_xxz[k] = -g_x_0_xxxxx_xxz[k] * ab_z + g_x_0_xxxxx_xxzz[k];

                g_x_0_xxxxxz_xyy[k] = -g_x_0_xxxxx_xyy[k] * ab_z + g_x_0_xxxxx_xyyz[k];

                g_x_0_xxxxxz_xyz[k] = -g_x_0_xxxxx_xyz[k] * ab_z + g_x_0_xxxxx_xyzz[k];

                g_x_0_xxxxxz_xzz[k] = -g_x_0_xxxxx_xzz[k] * ab_z + g_x_0_xxxxx_xzzz[k];

                g_x_0_xxxxxz_yyy[k] = -g_x_0_xxxxx_yyy[k] * ab_z + g_x_0_xxxxx_yyyz[k];

                g_x_0_xxxxxz_yyz[k] = -g_x_0_xxxxx_yyz[k] * ab_z + g_x_0_xxxxx_yyzz[k];

                g_x_0_xxxxxz_yzz[k] = -g_x_0_xxxxx_yzz[k] * ab_z + g_x_0_xxxxx_yzzz[k];

                g_x_0_xxxxxz_zzz[k] = -g_x_0_xxxxx_zzz[k] * ab_z + g_x_0_xxxxx_zzzz[k];
            }

            /// Set up 30-40 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_xxx = cbuffer.data(if_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxy = cbuffer.data(if_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxz = cbuffer.data(if_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyy = cbuffer.data(if_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyz = cbuffer.data(if_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xzz = cbuffer.data(if_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyy = cbuffer.data(if_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyz = cbuffer.data(if_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yzz = cbuffer.data(if_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxxyy_zzz = cbuffer.data(if_geom_10_off + 39 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxy_xxx, g_x_0_xxxxy_xxxy, g_x_0_xxxxy_xxy, g_x_0_xxxxy_xxyy, g_x_0_xxxxy_xxyz, g_x_0_xxxxy_xxz, g_x_0_xxxxy_xyy, g_x_0_xxxxy_xyyy, g_x_0_xxxxy_xyyz, g_x_0_xxxxy_xyz, g_x_0_xxxxy_xyzz, g_x_0_xxxxy_xzz, g_x_0_xxxxy_yyy, g_x_0_xxxxy_yyyy, g_x_0_xxxxy_yyyz, g_x_0_xxxxy_yyz, g_x_0_xxxxy_yyzz, g_x_0_xxxxy_yzz, g_x_0_xxxxy_yzzz, g_x_0_xxxxy_zzz, g_x_0_xxxxyy_xxx, g_x_0_xxxxyy_xxy, g_x_0_xxxxyy_xxz, g_x_0_xxxxyy_xyy, g_x_0_xxxxyy_xyz, g_x_0_xxxxyy_xzz, g_x_0_xxxxyy_yyy, g_x_0_xxxxyy_yyz, g_x_0_xxxxyy_yzz, g_x_0_xxxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_xxx[k] = -g_x_0_xxxxy_xxx[k] * ab_y + g_x_0_xxxxy_xxxy[k];

                g_x_0_xxxxyy_xxy[k] = -g_x_0_xxxxy_xxy[k] * ab_y + g_x_0_xxxxy_xxyy[k];

                g_x_0_xxxxyy_xxz[k] = -g_x_0_xxxxy_xxz[k] * ab_y + g_x_0_xxxxy_xxyz[k];

                g_x_0_xxxxyy_xyy[k] = -g_x_0_xxxxy_xyy[k] * ab_y + g_x_0_xxxxy_xyyy[k];

                g_x_0_xxxxyy_xyz[k] = -g_x_0_xxxxy_xyz[k] * ab_y + g_x_0_xxxxy_xyyz[k];

                g_x_0_xxxxyy_xzz[k] = -g_x_0_xxxxy_xzz[k] * ab_y + g_x_0_xxxxy_xyzz[k];

                g_x_0_xxxxyy_yyy[k] = -g_x_0_xxxxy_yyy[k] * ab_y + g_x_0_xxxxy_yyyy[k];

                g_x_0_xxxxyy_yyz[k] = -g_x_0_xxxxy_yyz[k] * ab_y + g_x_0_xxxxy_yyyz[k];

                g_x_0_xxxxyy_yzz[k] = -g_x_0_xxxxy_yzz[k] * ab_y + g_x_0_xxxxy_yyzz[k];

                g_x_0_xxxxyy_zzz[k] = -g_x_0_xxxxy_zzz[k] * ab_y + g_x_0_xxxxy_yzzz[k];
            }

            /// Set up 40-50 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_xxx = cbuffer.data(if_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxy = cbuffer.data(if_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxz = cbuffer.data(if_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyy = cbuffer.data(if_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyz = cbuffer.data(if_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xzz = cbuffer.data(if_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyy = cbuffer.data(if_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyz = cbuffer.data(if_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yzz = cbuffer.data(if_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxxyz_zzz = cbuffer.data(if_geom_10_off + 49 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxyz_xxx, g_x_0_xxxxyz_xxy, g_x_0_xxxxyz_xxz, g_x_0_xxxxyz_xyy, g_x_0_xxxxyz_xyz, g_x_0_xxxxyz_xzz, g_x_0_xxxxyz_yyy, g_x_0_xxxxyz_yyz, g_x_0_xxxxyz_yzz, g_x_0_xxxxyz_zzz, g_x_0_xxxxz_xxx, g_x_0_xxxxz_xxxy, g_x_0_xxxxz_xxy, g_x_0_xxxxz_xxyy, g_x_0_xxxxz_xxyz, g_x_0_xxxxz_xxz, g_x_0_xxxxz_xyy, g_x_0_xxxxz_xyyy, g_x_0_xxxxz_xyyz, g_x_0_xxxxz_xyz, g_x_0_xxxxz_xyzz, g_x_0_xxxxz_xzz, g_x_0_xxxxz_yyy, g_x_0_xxxxz_yyyy, g_x_0_xxxxz_yyyz, g_x_0_xxxxz_yyz, g_x_0_xxxxz_yyzz, g_x_0_xxxxz_yzz, g_x_0_xxxxz_yzzz, g_x_0_xxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_xxx[k] = -g_x_0_xxxxz_xxx[k] * ab_y + g_x_0_xxxxz_xxxy[k];

                g_x_0_xxxxyz_xxy[k] = -g_x_0_xxxxz_xxy[k] * ab_y + g_x_0_xxxxz_xxyy[k];

                g_x_0_xxxxyz_xxz[k] = -g_x_0_xxxxz_xxz[k] * ab_y + g_x_0_xxxxz_xxyz[k];

                g_x_0_xxxxyz_xyy[k] = -g_x_0_xxxxz_xyy[k] * ab_y + g_x_0_xxxxz_xyyy[k];

                g_x_0_xxxxyz_xyz[k] = -g_x_0_xxxxz_xyz[k] * ab_y + g_x_0_xxxxz_xyyz[k];

                g_x_0_xxxxyz_xzz[k] = -g_x_0_xxxxz_xzz[k] * ab_y + g_x_0_xxxxz_xyzz[k];

                g_x_0_xxxxyz_yyy[k] = -g_x_0_xxxxz_yyy[k] * ab_y + g_x_0_xxxxz_yyyy[k];

                g_x_0_xxxxyz_yyz[k] = -g_x_0_xxxxz_yyz[k] * ab_y + g_x_0_xxxxz_yyyz[k];

                g_x_0_xxxxyz_yzz[k] = -g_x_0_xxxxz_yzz[k] * ab_y + g_x_0_xxxxz_yyzz[k];

                g_x_0_xxxxyz_zzz[k] = -g_x_0_xxxxz_zzz[k] * ab_y + g_x_0_xxxxz_yzzz[k];
            }

            /// Set up 50-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_xxx = cbuffer.data(if_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxy = cbuffer.data(if_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxz = cbuffer.data(if_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyy = cbuffer.data(if_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyz = cbuffer.data(if_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xzz = cbuffer.data(if_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyy = cbuffer.data(if_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyz = cbuffer.data(if_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yzz = cbuffer.data(if_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxxzz_zzz = cbuffer.data(if_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxz_xxx, g_x_0_xxxxz_xxxz, g_x_0_xxxxz_xxy, g_x_0_xxxxz_xxyz, g_x_0_xxxxz_xxz, g_x_0_xxxxz_xxzz, g_x_0_xxxxz_xyy, g_x_0_xxxxz_xyyz, g_x_0_xxxxz_xyz, g_x_0_xxxxz_xyzz, g_x_0_xxxxz_xzz, g_x_0_xxxxz_xzzz, g_x_0_xxxxz_yyy, g_x_0_xxxxz_yyyz, g_x_0_xxxxz_yyz, g_x_0_xxxxz_yyzz, g_x_0_xxxxz_yzz, g_x_0_xxxxz_yzzz, g_x_0_xxxxz_zzz, g_x_0_xxxxz_zzzz, g_x_0_xxxxzz_xxx, g_x_0_xxxxzz_xxy, g_x_0_xxxxzz_xxz, g_x_0_xxxxzz_xyy, g_x_0_xxxxzz_xyz, g_x_0_xxxxzz_xzz, g_x_0_xxxxzz_yyy, g_x_0_xxxxzz_yyz, g_x_0_xxxxzz_yzz, g_x_0_xxxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_xxx[k] = -g_x_0_xxxxz_xxx[k] * ab_z + g_x_0_xxxxz_xxxz[k];

                g_x_0_xxxxzz_xxy[k] = -g_x_0_xxxxz_xxy[k] * ab_z + g_x_0_xxxxz_xxyz[k];

                g_x_0_xxxxzz_xxz[k] = -g_x_0_xxxxz_xxz[k] * ab_z + g_x_0_xxxxz_xxzz[k];

                g_x_0_xxxxzz_xyy[k] = -g_x_0_xxxxz_xyy[k] * ab_z + g_x_0_xxxxz_xyyz[k];

                g_x_0_xxxxzz_xyz[k] = -g_x_0_xxxxz_xyz[k] * ab_z + g_x_0_xxxxz_xyzz[k];

                g_x_0_xxxxzz_xzz[k] = -g_x_0_xxxxz_xzz[k] * ab_z + g_x_0_xxxxz_xzzz[k];

                g_x_0_xxxxzz_yyy[k] = -g_x_0_xxxxz_yyy[k] * ab_z + g_x_0_xxxxz_yyyz[k];

                g_x_0_xxxxzz_yyz[k] = -g_x_0_xxxxz_yyz[k] * ab_z + g_x_0_xxxxz_yyzz[k];

                g_x_0_xxxxzz_yzz[k] = -g_x_0_xxxxz_yzz[k] * ab_z + g_x_0_xxxxz_yzzz[k];

                g_x_0_xxxxzz_zzz[k] = -g_x_0_xxxxz_zzz[k] * ab_z + g_x_0_xxxxz_zzzz[k];
            }

            /// Set up 60-70 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_xxx = cbuffer.data(if_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxy = cbuffer.data(if_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxz = cbuffer.data(if_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyy = cbuffer.data(if_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyz = cbuffer.data(if_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xzz = cbuffer.data(if_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyy = cbuffer.data(if_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyz = cbuffer.data(if_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yzz = cbuffer.data(if_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxyyy_zzz = cbuffer.data(if_geom_10_off + 69 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyy_xxx, g_x_0_xxxyy_xxxy, g_x_0_xxxyy_xxy, g_x_0_xxxyy_xxyy, g_x_0_xxxyy_xxyz, g_x_0_xxxyy_xxz, g_x_0_xxxyy_xyy, g_x_0_xxxyy_xyyy, g_x_0_xxxyy_xyyz, g_x_0_xxxyy_xyz, g_x_0_xxxyy_xyzz, g_x_0_xxxyy_xzz, g_x_0_xxxyy_yyy, g_x_0_xxxyy_yyyy, g_x_0_xxxyy_yyyz, g_x_0_xxxyy_yyz, g_x_0_xxxyy_yyzz, g_x_0_xxxyy_yzz, g_x_0_xxxyy_yzzz, g_x_0_xxxyy_zzz, g_x_0_xxxyyy_xxx, g_x_0_xxxyyy_xxy, g_x_0_xxxyyy_xxz, g_x_0_xxxyyy_xyy, g_x_0_xxxyyy_xyz, g_x_0_xxxyyy_xzz, g_x_0_xxxyyy_yyy, g_x_0_xxxyyy_yyz, g_x_0_xxxyyy_yzz, g_x_0_xxxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_xxx[k] = -g_x_0_xxxyy_xxx[k] * ab_y + g_x_0_xxxyy_xxxy[k];

                g_x_0_xxxyyy_xxy[k] = -g_x_0_xxxyy_xxy[k] * ab_y + g_x_0_xxxyy_xxyy[k];

                g_x_0_xxxyyy_xxz[k] = -g_x_0_xxxyy_xxz[k] * ab_y + g_x_0_xxxyy_xxyz[k];

                g_x_0_xxxyyy_xyy[k] = -g_x_0_xxxyy_xyy[k] * ab_y + g_x_0_xxxyy_xyyy[k];

                g_x_0_xxxyyy_xyz[k] = -g_x_0_xxxyy_xyz[k] * ab_y + g_x_0_xxxyy_xyyz[k];

                g_x_0_xxxyyy_xzz[k] = -g_x_0_xxxyy_xzz[k] * ab_y + g_x_0_xxxyy_xyzz[k];

                g_x_0_xxxyyy_yyy[k] = -g_x_0_xxxyy_yyy[k] * ab_y + g_x_0_xxxyy_yyyy[k];

                g_x_0_xxxyyy_yyz[k] = -g_x_0_xxxyy_yyz[k] * ab_y + g_x_0_xxxyy_yyyz[k];

                g_x_0_xxxyyy_yzz[k] = -g_x_0_xxxyy_yzz[k] * ab_y + g_x_0_xxxyy_yyzz[k];

                g_x_0_xxxyyy_zzz[k] = -g_x_0_xxxyy_zzz[k] * ab_y + g_x_0_xxxyy_yzzz[k];
            }

            /// Set up 70-80 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_xxx = cbuffer.data(if_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxy = cbuffer.data(if_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxz = cbuffer.data(if_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyy = cbuffer.data(if_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyz = cbuffer.data(if_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xzz = cbuffer.data(if_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyy = cbuffer.data(if_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyz = cbuffer.data(if_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yzz = cbuffer.data(if_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxyyz_zzz = cbuffer.data(if_geom_10_off + 79 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyyz_xxx, g_x_0_xxxyyz_xxy, g_x_0_xxxyyz_xxz, g_x_0_xxxyyz_xyy, g_x_0_xxxyyz_xyz, g_x_0_xxxyyz_xzz, g_x_0_xxxyyz_yyy, g_x_0_xxxyyz_yyz, g_x_0_xxxyyz_yzz, g_x_0_xxxyyz_zzz, g_x_0_xxxyz_xxx, g_x_0_xxxyz_xxxy, g_x_0_xxxyz_xxy, g_x_0_xxxyz_xxyy, g_x_0_xxxyz_xxyz, g_x_0_xxxyz_xxz, g_x_0_xxxyz_xyy, g_x_0_xxxyz_xyyy, g_x_0_xxxyz_xyyz, g_x_0_xxxyz_xyz, g_x_0_xxxyz_xyzz, g_x_0_xxxyz_xzz, g_x_0_xxxyz_yyy, g_x_0_xxxyz_yyyy, g_x_0_xxxyz_yyyz, g_x_0_xxxyz_yyz, g_x_0_xxxyz_yyzz, g_x_0_xxxyz_yzz, g_x_0_xxxyz_yzzz, g_x_0_xxxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_xxx[k] = -g_x_0_xxxyz_xxx[k] * ab_y + g_x_0_xxxyz_xxxy[k];

                g_x_0_xxxyyz_xxy[k] = -g_x_0_xxxyz_xxy[k] * ab_y + g_x_0_xxxyz_xxyy[k];

                g_x_0_xxxyyz_xxz[k] = -g_x_0_xxxyz_xxz[k] * ab_y + g_x_0_xxxyz_xxyz[k];

                g_x_0_xxxyyz_xyy[k] = -g_x_0_xxxyz_xyy[k] * ab_y + g_x_0_xxxyz_xyyy[k];

                g_x_0_xxxyyz_xyz[k] = -g_x_0_xxxyz_xyz[k] * ab_y + g_x_0_xxxyz_xyyz[k];

                g_x_0_xxxyyz_xzz[k] = -g_x_0_xxxyz_xzz[k] * ab_y + g_x_0_xxxyz_xyzz[k];

                g_x_0_xxxyyz_yyy[k] = -g_x_0_xxxyz_yyy[k] * ab_y + g_x_0_xxxyz_yyyy[k];

                g_x_0_xxxyyz_yyz[k] = -g_x_0_xxxyz_yyz[k] * ab_y + g_x_0_xxxyz_yyyz[k];

                g_x_0_xxxyyz_yzz[k] = -g_x_0_xxxyz_yzz[k] * ab_y + g_x_0_xxxyz_yyzz[k];

                g_x_0_xxxyyz_zzz[k] = -g_x_0_xxxyz_zzz[k] * ab_y + g_x_0_xxxyz_yzzz[k];
            }

            /// Set up 80-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_xxx = cbuffer.data(if_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxy = cbuffer.data(if_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxz = cbuffer.data(if_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyy = cbuffer.data(if_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyz = cbuffer.data(if_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xzz = cbuffer.data(if_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyy = cbuffer.data(if_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyz = cbuffer.data(if_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yzz = cbuffer.data(if_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxyzz_zzz = cbuffer.data(if_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyzz_xxx, g_x_0_xxxyzz_xxy, g_x_0_xxxyzz_xxz, g_x_0_xxxyzz_xyy, g_x_0_xxxyzz_xyz, g_x_0_xxxyzz_xzz, g_x_0_xxxyzz_yyy, g_x_0_xxxyzz_yyz, g_x_0_xxxyzz_yzz, g_x_0_xxxyzz_zzz, g_x_0_xxxzz_xxx, g_x_0_xxxzz_xxxy, g_x_0_xxxzz_xxy, g_x_0_xxxzz_xxyy, g_x_0_xxxzz_xxyz, g_x_0_xxxzz_xxz, g_x_0_xxxzz_xyy, g_x_0_xxxzz_xyyy, g_x_0_xxxzz_xyyz, g_x_0_xxxzz_xyz, g_x_0_xxxzz_xyzz, g_x_0_xxxzz_xzz, g_x_0_xxxzz_yyy, g_x_0_xxxzz_yyyy, g_x_0_xxxzz_yyyz, g_x_0_xxxzz_yyz, g_x_0_xxxzz_yyzz, g_x_0_xxxzz_yzz, g_x_0_xxxzz_yzzz, g_x_0_xxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_xxx[k] = -g_x_0_xxxzz_xxx[k] * ab_y + g_x_0_xxxzz_xxxy[k];

                g_x_0_xxxyzz_xxy[k] = -g_x_0_xxxzz_xxy[k] * ab_y + g_x_0_xxxzz_xxyy[k];

                g_x_0_xxxyzz_xxz[k] = -g_x_0_xxxzz_xxz[k] * ab_y + g_x_0_xxxzz_xxyz[k];

                g_x_0_xxxyzz_xyy[k] = -g_x_0_xxxzz_xyy[k] * ab_y + g_x_0_xxxzz_xyyy[k];

                g_x_0_xxxyzz_xyz[k] = -g_x_0_xxxzz_xyz[k] * ab_y + g_x_0_xxxzz_xyyz[k];

                g_x_0_xxxyzz_xzz[k] = -g_x_0_xxxzz_xzz[k] * ab_y + g_x_0_xxxzz_xyzz[k];

                g_x_0_xxxyzz_yyy[k] = -g_x_0_xxxzz_yyy[k] * ab_y + g_x_0_xxxzz_yyyy[k];

                g_x_0_xxxyzz_yyz[k] = -g_x_0_xxxzz_yyz[k] * ab_y + g_x_0_xxxzz_yyyz[k];

                g_x_0_xxxyzz_yzz[k] = -g_x_0_xxxzz_yzz[k] * ab_y + g_x_0_xxxzz_yyzz[k];

                g_x_0_xxxyzz_zzz[k] = -g_x_0_xxxzz_zzz[k] * ab_y + g_x_0_xxxzz_yzzz[k];
            }

            /// Set up 90-100 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_xxx = cbuffer.data(if_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxy = cbuffer.data(if_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxz = cbuffer.data(if_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyy = cbuffer.data(if_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyz = cbuffer.data(if_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xzz = cbuffer.data(if_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyy = cbuffer.data(if_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyz = cbuffer.data(if_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yzz = cbuffer.data(if_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxxzzz_zzz = cbuffer.data(if_geom_10_off + 99 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxzz_xxx, g_x_0_xxxzz_xxxz, g_x_0_xxxzz_xxy, g_x_0_xxxzz_xxyz, g_x_0_xxxzz_xxz, g_x_0_xxxzz_xxzz, g_x_0_xxxzz_xyy, g_x_0_xxxzz_xyyz, g_x_0_xxxzz_xyz, g_x_0_xxxzz_xyzz, g_x_0_xxxzz_xzz, g_x_0_xxxzz_xzzz, g_x_0_xxxzz_yyy, g_x_0_xxxzz_yyyz, g_x_0_xxxzz_yyz, g_x_0_xxxzz_yyzz, g_x_0_xxxzz_yzz, g_x_0_xxxzz_yzzz, g_x_0_xxxzz_zzz, g_x_0_xxxzz_zzzz, g_x_0_xxxzzz_xxx, g_x_0_xxxzzz_xxy, g_x_0_xxxzzz_xxz, g_x_0_xxxzzz_xyy, g_x_0_xxxzzz_xyz, g_x_0_xxxzzz_xzz, g_x_0_xxxzzz_yyy, g_x_0_xxxzzz_yyz, g_x_0_xxxzzz_yzz, g_x_0_xxxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_xxx[k] = -g_x_0_xxxzz_xxx[k] * ab_z + g_x_0_xxxzz_xxxz[k];

                g_x_0_xxxzzz_xxy[k] = -g_x_0_xxxzz_xxy[k] * ab_z + g_x_0_xxxzz_xxyz[k];

                g_x_0_xxxzzz_xxz[k] = -g_x_0_xxxzz_xxz[k] * ab_z + g_x_0_xxxzz_xxzz[k];

                g_x_0_xxxzzz_xyy[k] = -g_x_0_xxxzz_xyy[k] * ab_z + g_x_0_xxxzz_xyyz[k];

                g_x_0_xxxzzz_xyz[k] = -g_x_0_xxxzz_xyz[k] * ab_z + g_x_0_xxxzz_xyzz[k];

                g_x_0_xxxzzz_xzz[k] = -g_x_0_xxxzz_xzz[k] * ab_z + g_x_0_xxxzz_xzzz[k];

                g_x_0_xxxzzz_yyy[k] = -g_x_0_xxxzz_yyy[k] * ab_z + g_x_0_xxxzz_yyyz[k];

                g_x_0_xxxzzz_yyz[k] = -g_x_0_xxxzz_yyz[k] * ab_z + g_x_0_xxxzz_yyzz[k];

                g_x_0_xxxzzz_yzz[k] = -g_x_0_xxxzz_yzz[k] * ab_z + g_x_0_xxxzz_yzzz[k];

                g_x_0_xxxzzz_zzz[k] = -g_x_0_xxxzz_zzz[k] * ab_z + g_x_0_xxxzz_zzzz[k];
            }

            /// Set up 100-110 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_xxx = cbuffer.data(if_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxy = cbuffer.data(if_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxz = cbuffer.data(if_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyy = cbuffer.data(if_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyz = cbuffer.data(if_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xzz = cbuffer.data(if_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyy = cbuffer.data(if_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyz = cbuffer.data(if_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yzz = cbuffer.data(if_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxyyyy_zzz = cbuffer.data(if_geom_10_off + 109 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyy_xxx, g_x_0_xxyyy_xxxy, g_x_0_xxyyy_xxy, g_x_0_xxyyy_xxyy, g_x_0_xxyyy_xxyz, g_x_0_xxyyy_xxz, g_x_0_xxyyy_xyy, g_x_0_xxyyy_xyyy, g_x_0_xxyyy_xyyz, g_x_0_xxyyy_xyz, g_x_0_xxyyy_xyzz, g_x_0_xxyyy_xzz, g_x_0_xxyyy_yyy, g_x_0_xxyyy_yyyy, g_x_0_xxyyy_yyyz, g_x_0_xxyyy_yyz, g_x_0_xxyyy_yyzz, g_x_0_xxyyy_yzz, g_x_0_xxyyy_yzzz, g_x_0_xxyyy_zzz, g_x_0_xxyyyy_xxx, g_x_0_xxyyyy_xxy, g_x_0_xxyyyy_xxz, g_x_0_xxyyyy_xyy, g_x_0_xxyyyy_xyz, g_x_0_xxyyyy_xzz, g_x_0_xxyyyy_yyy, g_x_0_xxyyyy_yyz, g_x_0_xxyyyy_yzz, g_x_0_xxyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_xxx[k] = -g_x_0_xxyyy_xxx[k] * ab_y + g_x_0_xxyyy_xxxy[k];

                g_x_0_xxyyyy_xxy[k] = -g_x_0_xxyyy_xxy[k] * ab_y + g_x_0_xxyyy_xxyy[k];

                g_x_0_xxyyyy_xxz[k] = -g_x_0_xxyyy_xxz[k] * ab_y + g_x_0_xxyyy_xxyz[k];

                g_x_0_xxyyyy_xyy[k] = -g_x_0_xxyyy_xyy[k] * ab_y + g_x_0_xxyyy_xyyy[k];

                g_x_0_xxyyyy_xyz[k] = -g_x_0_xxyyy_xyz[k] * ab_y + g_x_0_xxyyy_xyyz[k];

                g_x_0_xxyyyy_xzz[k] = -g_x_0_xxyyy_xzz[k] * ab_y + g_x_0_xxyyy_xyzz[k];

                g_x_0_xxyyyy_yyy[k] = -g_x_0_xxyyy_yyy[k] * ab_y + g_x_0_xxyyy_yyyy[k];

                g_x_0_xxyyyy_yyz[k] = -g_x_0_xxyyy_yyz[k] * ab_y + g_x_0_xxyyy_yyyz[k];

                g_x_0_xxyyyy_yzz[k] = -g_x_0_xxyyy_yzz[k] * ab_y + g_x_0_xxyyy_yyzz[k];

                g_x_0_xxyyyy_zzz[k] = -g_x_0_xxyyy_zzz[k] * ab_y + g_x_0_xxyyy_yzzz[k];
            }

            /// Set up 110-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_xxx = cbuffer.data(if_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxy = cbuffer.data(if_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxz = cbuffer.data(if_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyy = cbuffer.data(if_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyz = cbuffer.data(if_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xzz = cbuffer.data(if_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyy = cbuffer.data(if_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyz = cbuffer.data(if_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yzz = cbuffer.data(if_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxyyyz_zzz = cbuffer.data(if_geom_10_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyyz_xxx, g_x_0_xxyyyz_xxy, g_x_0_xxyyyz_xxz, g_x_0_xxyyyz_xyy, g_x_0_xxyyyz_xyz, g_x_0_xxyyyz_xzz, g_x_0_xxyyyz_yyy, g_x_0_xxyyyz_yyz, g_x_0_xxyyyz_yzz, g_x_0_xxyyyz_zzz, g_x_0_xxyyz_xxx, g_x_0_xxyyz_xxxy, g_x_0_xxyyz_xxy, g_x_0_xxyyz_xxyy, g_x_0_xxyyz_xxyz, g_x_0_xxyyz_xxz, g_x_0_xxyyz_xyy, g_x_0_xxyyz_xyyy, g_x_0_xxyyz_xyyz, g_x_0_xxyyz_xyz, g_x_0_xxyyz_xyzz, g_x_0_xxyyz_xzz, g_x_0_xxyyz_yyy, g_x_0_xxyyz_yyyy, g_x_0_xxyyz_yyyz, g_x_0_xxyyz_yyz, g_x_0_xxyyz_yyzz, g_x_0_xxyyz_yzz, g_x_0_xxyyz_yzzz, g_x_0_xxyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_xxx[k] = -g_x_0_xxyyz_xxx[k] * ab_y + g_x_0_xxyyz_xxxy[k];

                g_x_0_xxyyyz_xxy[k] = -g_x_0_xxyyz_xxy[k] * ab_y + g_x_0_xxyyz_xxyy[k];

                g_x_0_xxyyyz_xxz[k] = -g_x_0_xxyyz_xxz[k] * ab_y + g_x_0_xxyyz_xxyz[k];

                g_x_0_xxyyyz_xyy[k] = -g_x_0_xxyyz_xyy[k] * ab_y + g_x_0_xxyyz_xyyy[k];

                g_x_0_xxyyyz_xyz[k] = -g_x_0_xxyyz_xyz[k] * ab_y + g_x_0_xxyyz_xyyz[k];

                g_x_0_xxyyyz_xzz[k] = -g_x_0_xxyyz_xzz[k] * ab_y + g_x_0_xxyyz_xyzz[k];

                g_x_0_xxyyyz_yyy[k] = -g_x_0_xxyyz_yyy[k] * ab_y + g_x_0_xxyyz_yyyy[k];

                g_x_0_xxyyyz_yyz[k] = -g_x_0_xxyyz_yyz[k] * ab_y + g_x_0_xxyyz_yyyz[k];

                g_x_0_xxyyyz_yzz[k] = -g_x_0_xxyyz_yzz[k] * ab_y + g_x_0_xxyyz_yyzz[k];

                g_x_0_xxyyyz_zzz[k] = -g_x_0_xxyyz_zzz[k] * ab_y + g_x_0_xxyyz_yzzz[k];
            }

            /// Set up 120-130 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_xxx = cbuffer.data(if_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxy = cbuffer.data(if_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxz = cbuffer.data(if_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyy = cbuffer.data(if_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyz = cbuffer.data(if_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xzz = cbuffer.data(if_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyy = cbuffer.data(if_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyz = cbuffer.data(if_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yzz = cbuffer.data(if_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxyyzz_zzz = cbuffer.data(if_geom_10_off + 129 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyzz_xxx, g_x_0_xxyyzz_xxy, g_x_0_xxyyzz_xxz, g_x_0_xxyyzz_xyy, g_x_0_xxyyzz_xyz, g_x_0_xxyyzz_xzz, g_x_0_xxyyzz_yyy, g_x_0_xxyyzz_yyz, g_x_0_xxyyzz_yzz, g_x_0_xxyyzz_zzz, g_x_0_xxyzz_xxx, g_x_0_xxyzz_xxxy, g_x_0_xxyzz_xxy, g_x_0_xxyzz_xxyy, g_x_0_xxyzz_xxyz, g_x_0_xxyzz_xxz, g_x_0_xxyzz_xyy, g_x_0_xxyzz_xyyy, g_x_0_xxyzz_xyyz, g_x_0_xxyzz_xyz, g_x_0_xxyzz_xyzz, g_x_0_xxyzz_xzz, g_x_0_xxyzz_yyy, g_x_0_xxyzz_yyyy, g_x_0_xxyzz_yyyz, g_x_0_xxyzz_yyz, g_x_0_xxyzz_yyzz, g_x_0_xxyzz_yzz, g_x_0_xxyzz_yzzz, g_x_0_xxyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_xxx[k] = -g_x_0_xxyzz_xxx[k] * ab_y + g_x_0_xxyzz_xxxy[k];

                g_x_0_xxyyzz_xxy[k] = -g_x_0_xxyzz_xxy[k] * ab_y + g_x_0_xxyzz_xxyy[k];

                g_x_0_xxyyzz_xxz[k] = -g_x_0_xxyzz_xxz[k] * ab_y + g_x_0_xxyzz_xxyz[k];

                g_x_0_xxyyzz_xyy[k] = -g_x_0_xxyzz_xyy[k] * ab_y + g_x_0_xxyzz_xyyy[k];

                g_x_0_xxyyzz_xyz[k] = -g_x_0_xxyzz_xyz[k] * ab_y + g_x_0_xxyzz_xyyz[k];

                g_x_0_xxyyzz_xzz[k] = -g_x_0_xxyzz_xzz[k] * ab_y + g_x_0_xxyzz_xyzz[k];

                g_x_0_xxyyzz_yyy[k] = -g_x_0_xxyzz_yyy[k] * ab_y + g_x_0_xxyzz_yyyy[k];

                g_x_0_xxyyzz_yyz[k] = -g_x_0_xxyzz_yyz[k] * ab_y + g_x_0_xxyzz_yyyz[k];

                g_x_0_xxyyzz_yzz[k] = -g_x_0_xxyzz_yzz[k] * ab_y + g_x_0_xxyzz_yyzz[k];

                g_x_0_xxyyzz_zzz[k] = -g_x_0_xxyzz_zzz[k] * ab_y + g_x_0_xxyzz_yzzz[k];
            }

            /// Set up 130-140 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_xxx = cbuffer.data(if_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxy = cbuffer.data(if_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxz = cbuffer.data(if_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyy = cbuffer.data(if_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyz = cbuffer.data(if_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xzz = cbuffer.data(if_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyy = cbuffer.data(if_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyz = cbuffer.data(if_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yzz = cbuffer.data(if_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxyzzz_zzz = cbuffer.data(if_geom_10_off + 139 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyzzz_xxx, g_x_0_xxyzzz_xxy, g_x_0_xxyzzz_xxz, g_x_0_xxyzzz_xyy, g_x_0_xxyzzz_xyz, g_x_0_xxyzzz_xzz, g_x_0_xxyzzz_yyy, g_x_0_xxyzzz_yyz, g_x_0_xxyzzz_yzz, g_x_0_xxyzzz_zzz, g_x_0_xxzzz_xxx, g_x_0_xxzzz_xxxy, g_x_0_xxzzz_xxy, g_x_0_xxzzz_xxyy, g_x_0_xxzzz_xxyz, g_x_0_xxzzz_xxz, g_x_0_xxzzz_xyy, g_x_0_xxzzz_xyyy, g_x_0_xxzzz_xyyz, g_x_0_xxzzz_xyz, g_x_0_xxzzz_xyzz, g_x_0_xxzzz_xzz, g_x_0_xxzzz_yyy, g_x_0_xxzzz_yyyy, g_x_0_xxzzz_yyyz, g_x_0_xxzzz_yyz, g_x_0_xxzzz_yyzz, g_x_0_xxzzz_yzz, g_x_0_xxzzz_yzzz, g_x_0_xxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_xxx[k] = -g_x_0_xxzzz_xxx[k] * ab_y + g_x_0_xxzzz_xxxy[k];

                g_x_0_xxyzzz_xxy[k] = -g_x_0_xxzzz_xxy[k] * ab_y + g_x_0_xxzzz_xxyy[k];

                g_x_0_xxyzzz_xxz[k] = -g_x_0_xxzzz_xxz[k] * ab_y + g_x_0_xxzzz_xxyz[k];

                g_x_0_xxyzzz_xyy[k] = -g_x_0_xxzzz_xyy[k] * ab_y + g_x_0_xxzzz_xyyy[k];

                g_x_0_xxyzzz_xyz[k] = -g_x_0_xxzzz_xyz[k] * ab_y + g_x_0_xxzzz_xyyz[k];

                g_x_0_xxyzzz_xzz[k] = -g_x_0_xxzzz_xzz[k] * ab_y + g_x_0_xxzzz_xyzz[k];

                g_x_0_xxyzzz_yyy[k] = -g_x_0_xxzzz_yyy[k] * ab_y + g_x_0_xxzzz_yyyy[k];

                g_x_0_xxyzzz_yyz[k] = -g_x_0_xxzzz_yyz[k] * ab_y + g_x_0_xxzzz_yyyz[k];

                g_x_0_xxyzzz_yzz[k] = -g_x_0_xxzzz_yzz[k] * ab_y + g_x_0_xxzzz_yyzz[k];

                g_x_0_xxyzzz_zzz[k] = -g_x_0_xxzzz_zzz[k] * ab_y + g_x_0_xxzzz_yzzz[k];
            }

            /// Set up 140-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_xxx = cbuffer.data(if_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxy = cbuffer.data(if_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxz = cbuffer.data(if_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyy = cbuffer.data(if_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyz = cbuffer.data(if_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xzz = cbuffer.data(if_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyy = cbuffer.data(if_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyz = cbuffer.data(if_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yzz = cbuffer.data(if_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxzzzz_zzz = cbuffer.data(if_geom_10_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxzzz_xxx, g_x_0_xxzzz_xxxz, g_x_0_xxzzz_xxy, g_x_0_xxzzz_xxyz, g_x_0_xxzzz_xxz, g_x_0_xxzzz_xxzz, g_x_0_xxzzz_xyy, g_x_0_xxzzz_xyyz, g_x_0_xxzzz_xyz, g_x_0_xxzzz_xyzz, g_x_0_xxzzz_xzz, g_x_0_xxzzz_xzzz, g_x_0_xxzzz_yyy, g_x_0_xxzzz_yyyz, g_x_0_xxzzz_yyz, g_x_0_xxzzz_yyzz, g_x_0_xxzzz_yzz, g_x_0_xxzzz_yzzz, g_x_0_xxzzz_zzz, g_x_0_xxzzz_zzzz, g_x_0_xxzzzz_xxx, g_x_0_xxzzzz_xxy, g_x_0_xxzzzz_xxz, g_x_0_xxzzzz_xyy, g_x_0_xxzzzz_xyz, g_x_0_xxzzzz_xzz, g_x_0_xxzzzz_yyy, g_x_0_xxzzzz_yyz, g_x_0_xxzzzz_yzz, g_x_0_xxzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_xxx[k] = -g_x_0_xxzzz_xxx[k] * ab_z + g_x_0_xxzzz_xxxz[k];

                g_x_0_xxzzzz_xxy[k] = -g_x_0_xxzzz_xxy[k] * ab_z + g_x_0_xxzzz_xxyz[k];

                g_x_0_xxzzzz_xxz[k] = -g_x_0_xxzzz_xxz[k] * ab_z + g_x_0_xxzzz_xxzz[k];

                g_x_0_xxzzzz_xyy[k] = -g_x_0_xxzzz_xyy[k] * ab_z + g_x_0_xxzzz_xyyz[k];

                g_x_0_xxzzzz_xyz[k] = -g_x_0_xxzzz_xyz[k] * ab_z + g_x_0_xxzzz_xyzz[k];

                g_x_0_xxzzzz_xzz[k] = -g_x_0_xxzzz_xzz[k] * ab_z + g_x_0_xxzzz_xzzz[k];

                g_x_0_xxzzzz_yyy[k] = -g_x_0_xxzzz_yyy[k] * ab_z + g_x_0_xxzzz_yyyz[k];

                g_x_0_xxzzzz_yyz[k] = -g_x_0_xxzzz_yyz[k] * ab_z + g_x_0_xxzzz_yyzz[k];

                g_x_0_xxzzzz_yzz[k] = -g_x_0_xxzzz_yzz[k] * ab_z + g_x_0_xxzzz_yzzz[k];

                g_x_0_xxzzzz_zzz[k] = -g_x_0_xxzzz_zzz[k] * ab_z + g_x_0_xxzzz_zzzz[k];
            }

            /// Set up 150-160 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_xxx = cbuffer.data(if_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxy = cbuffer.data(if_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxz = cbuffer.data(if_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyy = cbuffer.data(if_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyz = cbuffer.data(if_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xzz = cbuffer.data(if_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyy = cbuffer.data(if_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyz = cbuffer.data(if_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yzz = cbuffer.data(if_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xyyyyy_zzz = cbuffer.data(if_geom_10_off + 159 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyy_xxx, g_x_0_xyyyy_xxxy, g_x_0_xyyyy_xxy, g_x_0_xyyyy_xxyy, g_x_0_xyyyy_xxyz, g_x_0_xyyyy_xxz, g_x_0_xyyyy_xyy, g_x_0_xyyyy_xyyy, g_x_0_xyyyy_xyyz, g_x_0_xyyyy_xyz, g_x_0_xyyyy_xyzz, g_x_0_xyyyy_xzz, g_x_0_xyyyy_yyy, g_x_0_xyyyy_yyyy, g_x_0_xyyyy_yyyz, g_x_0_xyyyy_yyz, g_x_0_xyyyy_yyzz, g_x_0_xyyyy_yzz, g_x_0_xyyyy_yzzz, g_x_0_xyyyy_zzz, g_x_0_xyyyyy_xxx, g_x_0_xyyyyy_xxy, g_x_0_xyyyyy_xxz, g_x_0_xyyyyy_xyy, g_x_0_xyyyyy_xyz, g_x_0_xyyyyy_xzz, g_x_0_xyyyyy_yyy, g_x_0_xyyyyy_yyz, g_x_0_xyyyyy_yzz, g_x_0_xyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_xxx[k] = -g_x_0_xyyyy_xxx[k] * ab_y + g_x_0_xyyyy_xxxy[k];

                g_x_0_xyyyyy_xxy[k] = -g_x_0_xyyyy_xxy[k] * ab_y + g_x_0_xyyyy_xxyy[k];

                g_x_0_xyyyyy_xxz[k] = -g_x_0_xyyyy_xxz[k] * ab_y + g_x_0_xyyyy_xxyz[k];

                g_x_0_xyyyyy_xyy[k] = -g_x_0_xyyyy_xyy[k] * ab_y + g_x_0_xyyyy_xyyy[k];

                g_x_0_xyyyyy_xyz[k] = -g_x_0_xyyyy_xyz[k] * ab_y + g_x_0_xyyyy_xyyz[k];

                g_x_0_xyyyyy_xzz[k] = -g_x_0_xyyyy_xzz[k] * ab_y + g_x_0_xyyyy_xyzz[k];

                g_x_0_xyyyyy_yyy[k] = -g_x_0_xyyyy_yyy[k] * ab_y + g_x_0_xyyyy_yyyy[k];

                g_x_0_xyyyyy_yyz[k] = -g_x_0_xyyyy_yyz[k] * ab_y + g_x_0_xyyyy_yyyz[k];

                g_x_0_xyyyyy_yzz[k] = -g_x_0_xyyyy_yzz[k] * ab_y + g_x_0_xyyyy_yyzz[k];

                g_x_0_xyyyyy_zzz[k] = -g_x_0_xyyyy_zzz[k] * ab_y + g_x_0_xyyyy_yzzz[k];
            }

            /// Set up 160-170 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_xxx = cbuffer.data(if_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxy = cbuffer.data(if_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxz = cbuffer.data(if_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyy = cbuffer.data(if_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyz = cbuffer.data(if_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xzz = cbuffer.data(if_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyy = cbuffer.data(if_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyz = cbuffer.data(if_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yzz = cbuffer.data(if_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xyyyyz_zzz = cbuffer.data(if_geom_10_off + 169 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyyz_xxx, g_x_0_xyyyyz_xxy, g_x_0_xyyyyz_xxz, g_x_0_xyyyyz_xyy, g_x_0_xyyyyz_xyz, g_x_0_xyyyyz_xzz, g_x_0_xyyyyz_yyy, g_x_0_xyyyyz_yyz, g_x_0_xyyyyz_yzz, g_x_0_xyyyyz_zzz, g_x_0_xyyyz_xxx, g_x_0_xyyyz_xxxy, g_x_0_xyyyz_xxy, g_x_0_xyyyz_xxyy, g_x_0_xyyyz_xxyz, g_x_0_xyyyz_xxz, g_x_0_xyyyz_xyy, g_x_0_xyyyz_xyyy, g_x_0_xyyyz_xyyz, g_x_0_xyyyz_xyz, g_x_0_xyyyz_xyzz, g_x_0_xyyyz_xzz, g_x_0_xyyyz_yyy, g_x_0_xyyyz_yyyy, g_x_0_xyyyz_yyyz, g_x_0_xyyyz_yyz, g_x_0_xyyyz_yyzz, g_x_0_xyyyz_yzz, g_x_0_xyyyz_yzzz, g_x_0_xyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_xxx[k] = -g_x_0_xyyyz_xxx[k] * ab_y + g_x_0_xyyyz_xxxy[k];

                g_x_0_xyyyyz_xxy[k] = -g_x_0_xyyyz_xxy[k] * ab_y + g_x_0_xyyyz_xxyy[k];

                g_x_0_xyyyyz_xxz[k] = -g_x_0_xyyyz_xxz[k] * ab_y + g_x_0_xyyyz_xxyz[k];

                g_x_0_xyyyyz_xyy[k] = -g_x_0_xyyyz_xyy[k] * ab_y + g_x_0_xyyyz_xyyy[k];

                g_x_0_xyyyyz_xyz[k] = -g_x_0_xyyyz_xyz[k] * ab_y + g_x_0_xyyyz_xyyz[k];

                g_x_0_xyyyyz_xzz[k] = -g_x_0_xyyyz_xzz[k] * ab_y + g_x_0_xyyyz_xyzz[k];

                g_x_0_xyyyyz_yyy[k] = -g_x_0_xyyyz_yyy[k] * ab_y + g_x_0_xyyyz_yyyy[k];

                g_x_0_xyyyyz_yyz[k] = -g_x_0_xyyyz_yyz[k] * ab_y + g_x_0_xyyyz_yyyz[k];

                g_x_0_xyyyyz_yzz[k] = -g_x_0_xyyyz_yzz[k] * ab_y + g_x_0_xyyyz_yyzz[k];

                g_x_0_xyyyyz_zzz[k] = -g_x_0_xyyyz_zzz[k] * ab_y + g_x_0_xyyyz_yzzz[k];
            }

            /// Set up 170-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_xxx = cbuffer.data(if_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxy = cbuffer.data(if_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxz = cbuffer.data(if_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyy = cbuffer.data(if_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyz = cbuffer.data(if_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xzz = cbuffer.data(if_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyy = cbuffer.data(if_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyz = cbuffer.data(if_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yzz = cbuffer.data(if_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xyyyzz_zzz = cbuffer.data(if_geom_10_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyzz_xxx, g_x_0_xyyyzz_xxy, g_x_0_xyyyzz_xxz, g_x_0_xyyyzz_xyy, g_x_0_xyyyzz_xyz, g_x_0_xyyyzz_xzz, g_x_0_xyyyzz_yyy, g_x_0_xyyyzz_yyz, g_x_0_xyyyzz_yzz, g_x_0_xyyyzz_zzz, g_x_0_xyyzz_xxx, g_x_0_xyyzz_xxxy, g_x_0_xyyzz_xxy, g_x_0_xyyzz_xxyy, g_x_0_xyyzz_xxyz, g_x_0_xyyzz_xxz, g_x_0_xyyzz_xyy, g_x_0_xyyzz_xyyy, g_x_0_xyyzz_xyyz, g_x_0_xyyzz_xyz, g_x_0_xyyzz_xyzz, g_x_0_xyyzz_xzz, g_x_0_xyyzz_yyy, g_x_0_xyyzz_yyyy, g_x_0_xyyzz_yyyz, g_x_0_xyyzz_yyz, g_x_0_xyyzz_yyzz, g_x_0_xyyzz_yzz, g_x_0_xyyzz_yzzz, g_x_0_xyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_xxx[k] = -g_x_0_xyyzz_xxx[k] * ab_y + g_x_0_xyyzz_xxxy[k];

                g_x_0_xyyyzz_xxy[k] = -g_x_0_xyyzz_xxy[k] * ab_y + g_x_0_xyyzz_xxyy[k];

                g_x_0_xyyyzz_xxz[k] = -g_x_0_xyyzz_xxz[k] * ab_y + g_x_0_xyyzz_xxyz[k];

                g_x_0_xyyyzz_xyy[k] = -g_x_0_xyyzz_xyy[k] * ab_y + g_x_0_xyyzz_xyyy[k];

                g_x_0_xyyyzz_xyz[k] = -g_x_0_xyyzz_xyz[k] * ab_y + g_x_0_xyyzz_xyyz[k];

                g_x_0_xyyyzz_xzz[k] = -g_x_0_xyyzz_xzz[k] * ab_y + g_x_0_xyyzz_xyzz[k];

                g_x_0_xyyyzz_yyy[k] = -g_x_0_xyyzz_yyy[k] * ab_y + g_x_0_xyyzz_yyyy[k];

                g_x_0_xyyyzz_yyz[k] = -g_x_0_xyyzz_yyz[k] * ab_y + g_x_0_xyyzz_yyyz[k];

                g_x_0_xyyyzz_yzz[k] = -g_x_0_xyyzz_yzz[k] * ab_y + g_x_0_xyyzz_yyzz[k];

                g_x_0_xyyyzz_zzz[k] = -g_x_0_xyyzz_zzz[k] * ab_y + g_x_0_xyyzz_yzzz[k];
            }

            /// Set up 180-190 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_xxx = cbuffer.data(if_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxy = cbuffer.data(if_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxz = cbuffer.data(if_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyy = cbuffer.data(if_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyz = cbuffer.data(if_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xzz = cbuffer.data(if_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyy = cbuffer.data(if_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyz = cbuffer.data(if_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yzz = cbuffer.data(if_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xyyzzz_zzz = cbuffer.data(if_geom_10_off + 189 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyzzz_xxx, g_x_0_xyyzzz_xxy, g_x_0_xyyzzz_xxz, g_x_0_xyyzzz_xyy, g_x_0_xyyzzz_xyz, g_x_0_xyyzzz_xzz, g_x_0_xyyzzz_yyy, g_x_0_xyyzzz_yyz, g_x_0_xyyzzz_yzz, g_x_0_xyyzzz_zzz, g_x_0_xyzzz_xxx, g_x_0_xyzzz_xxxy, g_x_0_xyzzz_xxy, g_x_0_xyzzz_xxyy, g_x_0_xyzzz_xxyz, g_x_0_xyzzz_xxz, g_x_0_xyzzz_xyy, g_x_0_xyzzz_xyyy, g_x_0_xyzzz_xyyz, g_x_0_xyzzz_xyz, g_x_0_xyzzz_xyzz, g_x_0_xyzzz_xzz, g_x_0_xyzzz_yyy, g_x_0_xyzzz_yyyy, g_x_0_xyzzz_yyyz, g_x_0_xyzzz_yyz, g_x_0_xyzzz_yyzz, g_x_0_xyzzz_yzz, g_x_0_xyzzz_yzzz, g_x_0_xyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_xxx[k] = -g_x_0_xyzzz_xxx[k] * ab_y + g_x_0_xyzzz_xxxy[k];

                g_x_0_xyyzzz_xxy[k] = -g_x_0_xyzzz_xxy[k] * ab_y + g_x_0_xyzzz_xxyy[k];

                g_x_0_xyyzzz_xxz[k] = -g_x_0_xyzzz_xxz[k] * ab_y + g_x_0_xyzzz_xxyz[k];

                g_x_0_xyyzzz_xyy[k] = -g_x_0_xyzzz_xyy[k] * ab_y + g_x_0_xyzzz_xyyy[k];

                g_x_0_xyyzzz_xyz[k] = -g_x_0_xyzzz_xyz[k] * ab_y + g_x_0_xyzzz_xyyz[k];

                g_x_0_xyyzzz_xzz[k] = -g_x_0_xyzzz_xzz[k] * ab_y + g_x_0_xyzzz_xyzz[k];

                g_x_0_xyyzzz_yyy[k] = -g_x_0_xyzzz_yyy[k] * ab_y + g_x_0_xyzzz_yyyy[k];

                g_x_0_xyyzzz_yyz[k] = -g_x_0_xyzzz_yyz[k] * ab_y + g_x_0_xyzzz_yyyz[k];

                g_x_0_xyyzzz_yzz[k] = -g_x_0_xyzzz_yzz[k] * ab_y + g_x_0_xyzzz_yyzz[k];

                g_x_0_xyyzzz_zzz[k] = -g_x_0_xyzzz_zzz[k] * ab_y + g_x_0_xyzzz_yzzz[k];
            }

            /// Set up 190-200 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_xxx = cbuffer.data(if_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxy = cbuffer.data(if_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxz = cbuffer.data(if_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyy = cbuffer.data(if_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyz = cbuffer.data(if_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xzz = cbuffer.data(if_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyy = cbuffer.data(if_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyz = cbuffer.data(if_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yzz = cbuffer.data(if_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xyzzzz_zzz = cbuffer.data(if_geom_10_off + 199 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzzzz_xxx, g_x_0_xyzzzz_xxy, g_x_0_xyzzzz_xxz, g_x_0_xyzzzz_xyy, g_x_0_xyzzzz_xyz, g_x_0_xyzzzz_xzz, g_x_0_xyzzzz_yyy, g_x_0_xyzzzz_yyz, g_x_0_xyzzzz_yzz, g_x_0_xyzzzz_zzz, g_x_0_xzzzz_xxx, g_x_0_xzzzz_xxxy, g_x_0_xzzzz_xxy, g_x_0_xzzzz_xxyy, g_x_0_xzzzz_xxyz, g_x_0_xzzzz_xxz, g_x_0_xzzzz_xyy, g_x_0_xzzzz_xyyy, g_x_0_xzzzz_xyyz, g_x_0_xzzzz_xyz, g_x_0_xzzzz_xyzz, g_x_0_xzzzz_xzz, g_x_0_xzzzz_yyy, g_x_0_xzzzz_yyyy, g_x_0_xzzzz_yyyz, g_x_0_xzzzz_yyz, g_x_0_xzzzz_yyzz, g_x_0_xzzzz_yzz, g_x_0_xzzzz_yzzz, g_x_0_xzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_xxx[k] = -g_x_0_xzzzz_xxx[k] * ab_y + g_x_0_xzzzz_xxxy[k];

                g_x_0_xyzzzz_xxy[k] = -g_x_0_xzzzz_xxy[k] * ab_y + g_x_0_xzzzz_xxyy[k];

                g_x_0_xyzzzz_xxz[k] = -g_x_0_xzzzz_xxz[k] * ab_y + g_x_0_xzzzz_xxyz[k];

                g_x_0_xyzzzz_xyy[k] = -g_x_0_xzzzz_xyy[k] * ab_y + g_x_0_xzzzz_xyyy[k];

                g_x_0_xyzzzz_xyz[k] = -g_x_0_xzzzz_xyz[k] * ab_y + g_x_0_xzzzz_xyyz[k];

                g_x_0_xyzzzz_xzz[k] = -g_x_0_xzzzz_xzz[k] * ab_y + g_x_0_xzzzz_xyzz[k];

                g_x_0_xyzzzz_yyy[k] = -g_x_0_xzzzz_yyy[k] * ab_y + g_x_0_xzzzz_yyyy[k];

                g_x_0_xyzzzz_yyz[k] = -g_x_0_xzzzz_yyz[k] * ab_y + g_x_0_xzzzz_yyyz[k];

                g_x_0_xyzzzz_yzz[k] = -g_x_0_xzzzz_yzz[k] * ab_y + g_x_0_xzzzz_yyzz[k];

                g_x_0_xyzzzz_zzz[k] = -g_x_0_xzzzz_zzz[k] * ab_y + g_x_0_xzzzz_yzzz[k];
            }

            /// Set up 200-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_xxx = cbuffer.data(if_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxy = cbuffer.data(if_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxz = cbuffer.data(if_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyy = cbuffer.data(if_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyz = cbuffer.data(if_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xzz = cbuffer.data(if_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyy = cbuffer.data(if_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyz = cbuffer.data(if_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yzz = cbuffer.data(if_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xzzzzz_zzz = cbuffer.data(if_geom_10_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzzzz_xxx, g_x_0_xzzzz_xxxz, g_x_0_xzzzz_xxy, g_x_0_xzzzz_xxyz, g_x_0_xzzzz_xxz, g_x_0_xzzzz_xxzz, g_x_0_xzzzz_xyy, g_x_0_xzzzz_xyyz, g_x_0_xzzzz_xyz, g_x_0_xzzzz_xyzz, g_x_0_xzzzz_xzz, g_x_0_xzzzz_xzzz, g_x_0_xzzzz_yyy, g_x_0_xzzzz_yyyz, g_x_0_xzzzz_yyz, g_x_0_xzzzz_yyzz, g_x_0_xzzzz_yzz, g_x_0_xzzzz_yzzz, g_x_0_xzzzz_zzz, g_x_0_xzzzz_zzzz, g_x_0_xzzzzz_xxx, g_x_0_xzzzzz_xxy, g_x_0_xzzzzz_xxz, g_x_0_xzzzzz_xyy, g_x_0_xzzzzz_xyz, g_x_0_xzzzzz_xzz, g_x_0_xzzzzz_yyy, g_x_0_xzzzzz_yyz, g_x_0_xzzzzz_yzz, g_x_0_xzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_xxx[k] = -g_x_0_xzzzz_xxx[k] * ab_z + g_x_0_xzzzz_xxxz[k];

                g_x_0_xzzzzz_xxy[k] = -g_x_0_xzzzz_xxy[k] * ab_z + g_x_0_xzzzz_xxyz[k];

                g_x_0_xzzzzz_xxz[k] = -g_x_0_xzzzz_xxz[k] * ab_z + g_x_0_xzzzz_xxzz[k];

                g_x_0_xzzzzz_xyy[k] = -g_x_0_xzzzz_xyy[k] * ab_z + g_x_0_xzzzz_xyyz[k];

                g_x_0_xzzzzz_xyz[k] = -g_x_0_xzzzz_xyz[k] * ab_z + g_x_0_xzzzz_xyzz[k];

                g_x_0_xzzzzz_xzz[k] = -g_x_0_xzzzz_xzz[k] * ab_z + g_x_0_xzzzz_xzzz[k];

                g_x_0_xzzzzz_yyy[k] = -g_x_0_xzzzz_yyy[k] * ab_z + g_x_0_xzzzz_yyyz[k];

                g_x_0_xzzzzz_yyz[k] = -g_x_0_xzzzz_yyz[k] * ab_z + g_x_0_xzzzz_yyzz[k];

                g_x_0_xzzzzz_yzz[k] = -g_x_0_xzzzz_yzz[k] * ab_z + g_x_0_xzzzz_yzzz[k];

                g_x_0_xzzzzz_zzz[k] = -g_x_0_xzzzz_zzz[k] * ab_z + g_x_0_xzzzz_zzzz[k];
            }

            /// Set up 210-220 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_xxx = cbuffer.data(if_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxy = cbuffer.data(if_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxz = cbuffer.data(if_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyy = cbuffer.data(if_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyz = cbuffer.data(if_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xzz = cbuffer.data(if_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyy = cbuffer.data(if_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyz = cbuffer.data(if_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yzz = cbuffer.data(if_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_yyyyyy_zzz = cbuffer.data(if_geom_10_off + 219 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyy_xxx, g_x_0_yyyyy_xxxy, g_x_0_yyyyy_xxy, g_x_0_yyyyy_xxyy, g_x_0_yyyyy_xxyz, g_x_0_yyyyy_xxz, g_x_0_yyyyy_xyy, g_x_0_yyyyy_xyyy, g_x_0_yyyyy_xyyz, g_x_0_yyyyy_xyz, g_x_0_yyyyy_xyzz, g_x_0_yyyyy_xzz, g_x_0_yyyyy_yyy, g_x_0_yyyyy_yyyy, g_x_0_yyyyy_yyyz, g_x_0_yyyyy_yyz, g_x_0_yyyyy_yyzz, g_x_0_yyyyy_yzz, g_x_0_yyyyy_yzzz, g_x_0_yyyyy_zzz, g_x_0_yyyyyy_xxx, g_x_0_yyyyyy_xxy, g_x_0_yyyyyy_xxz, g_x_0_yyyyyy_xyy, g_x_0_yyyyyy_xyz, g_x_0_yyyyyy_xzz, g_x_0_yyyyyy_yyy, g_x_0_yyyyyy_yyz, g_x_0_yyyyyy_yzz, g_x_0_yyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_xxx[k] = -g_x_0_yyyyy_xxx[k] * ab_y + g_x_0_yyyyy_xxxy[k];

                g_x_0_yyyyyy_xxy[k] = -g_x_0_yyyyy_xxy[k] * ab_y + g_x_0_yyyyy_xxyy[k];

                g_x_0_yyyyyy_xxz[k] = -g_x_0_yyyyy_xxz[k] * ab_y + g_x_0_yyyyy_xxyz[k];

                g_x_0_yyyyyy_xyy[k] = -g_x_0_yyyyy_xyy[k] * ab_y + g_x_0_yyyyy_xyyy[k];

                g_x_0_yyyyyy_xyz[k] = -g_x_0_yyyyy_xyz[k] * ab_y + g_x_0_yyyyy_xyyz[k];

                g_x_0_yyyyyy_xzz[k] = -g_x_0_yyyyy_xzz[k] * ab_y + g_x_0_yyyyy_xyzz[k];

                g_x_0_yyyyyy_yyy[k] = -g_x_0_yyyyy_yyy[k] * ab_y + g_x_0_yyyyy_yyyy[k];

                g_x_0_yyyyyy_yyz[k] = -g_x_0_yyyyy_yyz[k] * ab_y + g_x_0_yyyyy_yyyz[k];

                g_x_0_yyyyyy_yzz[k] = -g_x_0_yyyyy_yzz[k] * ab_y + g_x_0_yyyyy_yyzz[k];

                g_x_0_yyyyyy_zzz[k] = -g_x_0_yyyyy_zzz[k] * ab_y + g_x_0_yyyyy_yzzz[k];
            }

            /// Set up 220-230 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_xxx = cbuffer.data(if_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxy = cbuffer.data(if_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxz = cbuffer.data(if_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyy = cbuffer.data(if_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyz = cbuffer.data(if_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xzz = cbuffer.data(if_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyy = cbuffer.data(if_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyz = cbuffer.data(if_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yzz = cbuffer.data(if_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_yyyyyz_zzz = cbuffer.data(if_geom_10_off + 229 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyyz_xxx, g_x_0_yyyyyz_xxy, g_x_0_yyyyyz_xxz, g_x_0_yyyyyz_xyy, g_x_0_yyyyyz_xyz, g_x_0_yyyyyz_xzz, g_x_0_yyyyyz_yyy, g_x_0_yyyyyz_yyz, g_x_0_yyyyyz_yzz, g_x_0_yyyyyz_zzz, g_x_0_yyyyz_xxx, g_x_0_yyyyz_xxxy, g_x_0_yyyyz_xxy, g_x_0_yyyyz_xxyy, g_x_0_yyyyz_xxyz, g_x_0_yyyyz_xxz, g_x_0_yyyyz_xyy, g_x_0_yyyyz_xyyy, g_x_0_yyyyz_xyyz, g_x_0_yyyyz_xyz, g_x_0_yyyyz_xyzz, g_x_0_yyyyz_xzz, g_x_0_yyyyz_yyy, g_x_0_yyyyz_yyyy, g_x_0_yyyyz_yyyz, g_x_0_yyyyz_yyz, g_x_0_yyyyz_yyzz, g_x_0_yyyyz_yzz, g_x_0_yyyyz_yzzz, g_x_0_yyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_xxx[k] = -g_x_0_yyyyz_xxx[k] * ab_y + g_x_0_yyyyz_xxxy[k];

                g_x_0_yyyyyz_xxy[k] = -g_x_0_yyyyz_xxy[k] * ab_y + g_x_0_yyyyz_xxyy[k];

                g_x_0_yyyyyz_xxz[k] = -g_x_0_yyyyz_xxz[k] * ab_y + g_x_0_yyyyz_xxyz[k];

                g_x_0_yyyyyz_xyy[k] = -g_x_0_yyyyz_xyy[k] * ab_y + g_x_0_yyyyz_xyyy[k];

                g_x_0_yyyyyz_xyz[k] = -g_x_0_yyyyz_xyz[k] * ab_y + g_x_0_yyyyz_xyyz[k];

                g_x_0_yyyyyz_xzz[k] = -g_x_0_yyyyz_xzz[k] * ab_y + g_x_0_yyyyz_xyzz[k];

                g_x_0_yyyyyz_yyy[k] = -g_x_0_yyyyz_yyy[k] * ab_y + g_x_0_yyyyz_yyyy[k];

                g_x_0_yyyyyz_yyz[k] = -g_x_0_yyyyz_yyz[k] * ab_y + g_x_0_yyyyz_yyyz[k];

                g_x_0_yyyyyz_yzz[k] = -g_x_0_yyyyz_yzz[k] * ab_y + g_x_0_yyyyz_yyzz[k];

                g_x_0_yyyyyz_zzz[k] = -g_x_0_yyyyz_zzz[k] * ab_y + g_x_0_yyyyz_yzzz[k];
            }

            /// Set up 230-240 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_xxx = cbuffer.data(if_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxy = cbuffer.data(if_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxz = cbuffer.data(if_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyy = cbuffer.data(if_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyz = cbuffer.data(if_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xzz = cbuffer.data(if_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyy = cbuffer.data(if_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyz = cbuffer.data(if_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yzz = cbuffer.data(if_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_yyyyzz_zzz = cbuffer.data(if_geom_10_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyzz_xxx, g_x_0_yyyyzz_xxy, g_x_0_yyyyzz_xxz, g_x_0_yyyyzz_xyy, g_x_0_yyyyzz_xyz, g_x_0_yyyyzz_xzz, g_x_0_yyyyzz_yyy, g_x_0_yyyyzz_yyz, g_x_0_yyyyzz_yzz, g_x_0_yyyyzz_zzz, g_x_0_yyyzz_xxx, g_x_0_yyyzz_xxxy, g_x_0_yyyzz_xxy, g_x_0_yyyzz_xxyy, g_x_0_yyyzz_xxyz, g_x_0_yyyzz_xxz, g_x_0_yyyzz_xyy, g_x_0_yyyzz_xyyy, g_x_0_yyyzz_xyyz, g_x_0_yyyzz_xyz, g_x_0_yyyzz_xyzz, g_x_0_yyyzz_xzz, g_x_0_yyyzz_yyy, g_x_0_yyyzz_yyyy, g_x_0_yyyzz_yyyz, g_x_0_yyyzz_yyz, g_x_0_yyyzz_yyzz, g_x_0_yyyzz_yzz, g_x_0_yyyzz_yzzz, g_x_0_yyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_xxx[k] = -g_x_0_yyyzz_xxx[k] * ab_y + g_x_0_yyyzz_xxxy[k];

                g_x_0_yyyyzz_xxy[k] = -g_x_0_yyyzz_xxy[k] * ab_y + g_x_0_yyyzz_xxyy[k];

                g_x_0_yyyyzz_xxz[k] = -g_x_0_yyyzz_xxz[k] * ab_y + g_x_0_yyyzz_xxyz[k];

                g_x_0_yyyyzz_xyy[k] = -g_x_0_yyyzz_xyy[k] * ab_y + g_x_0_yyyzz_xyyy[k];

                g_x_0_yyyyzz_xyz[k] = -g_x_0_yyyzz_xyz[k] * ab_y + g_x_0_yyyzz_xyyz[k];

                g_x_0_yyyyzz_xzz[k] = -g_x_0_yyyzz_xzz[k] * ab_y + g_x_0_yyyzz_xyzz[k];

                g_x_0_yyyyzz_yyy[k] = -g_x_0_yyyzz_yyy[k] * ab_y + g_x_0_yyyzz_yyyy[k];

                g_x_0_yyyyzz_yyz[k] = -g_x_0_yyyzz_yyz[k] * ab_y + g_x_0_yyyzz_yyyz[k];

                g_x_0_yyyyzz_yzz[k] = -g_x_0_yyyzz_yzz[k] * ab_y + g_x_0_yyyzz_yyzz[k];

                g_x_0_yyyyzz_zzz[k] = -g_x_0_yyyzz_zzz[k] * ab_y + g_x_0_yyyzz_yzzz[k];
            }

            /// Set up 240-250 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_xxx = cbuffer.data(if_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxy = cbuffer.data(if_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxz = cbuffer.data(if_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyy = cbuffer.data(if_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyz = cbuffer.data(if_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xzz = cbuffer.data(if_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyy = cbuffer.data(if_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyz = cbuffer.data(if_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yzz = cbuffer.data(if_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_yyyzzz_zzz = cbuffer.data(if_geom_10_off + 249 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyzzz_xxx, g_x_0_yyyzzz_xxy, g_x_0_yyyzzz_xxz, g_x_0_yyyzzz_xyy, g_x_0_yyyzzz_xyz, g_x_0_yyyzzz_xzz, g_x_0_yyyzzz_yyy, g_x_0_yyyzzz_yyz, g_x_0_yyyzzz_yzz, g_x_0_yyyzzz_zzz, g_x_0_yyzzz_xxx, g_x_0_yyzzz_xxxy, g_x_0_yyzzz_xxy, g_x_0_yyzzz_xxyy, g_x_0_yyzzz_xxyz, g_x_0_yyzzz_xxz, g_x_0_yyzzz_xyy, g_x_0_yyzzz_xyyy, g_x_0_yyzzz_xyyz, g_x_0_yyzzz_xyz, g_x_0_yyzzz_xyzz, g_x_0_yyzzz_xzz, g_x_0_yyzzz_yyy, g_x_0_yyzzz_yyyy, g_x_0_yyzzz_yyyz, g_x_0_yyzzz_yyz, g_x_0_yyzzz_yyzz, g_x_0_yyzzz_yzz, g_x_0_yyzzz_yzzz, g_x_0_yyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_xxx[k] = -g_x_0_yyzzz_xxx[k] * ab_y + g_x_0_yyzzz_xxxy[k];

                g_x_0_yyyzzz_xxy[k] = -g_x_0_yyzzz_xxy[k] * ab_y + g_x_0_yyzzz_xxyy[k];

                g_x_0_yyyzzz_xxz[k] = -g_x_0_yyzzz_xxz[k] * ab_y + g_x_0_yyzzz_xxyz[k];

                g_x_0_yyyzzz_xyy[k] = -g_x_0_yyzzz_xyy[k] * ab_y + g_x_0_yyzzz_xyyy[k];

                g_x_0_yyyzzz_xyz[k] = -g_x_0_yyzzz_xyz[k] * ab_y + g_x_0_yyzzz_xyyz[k];

                g_x_0_yyyzzz_xzz[k] = -g_x_0_yyzzz_xzz[k] * ab_y + g_x_0_yyzzz_xyzz[k];

                g_x_0_yyyzzz_yyy[k] = -g_x_0_yyzzz_yyy[k] * ab_y + g_x_0_yyzzz_yyyy[k];

                g_x_0_yyyzzz_yyz[k] = -g_x_0_yyzzz_yyz[k] * ab_y + g_x_0_yyzzz_yyyz[k];

                g_x_0_yyyzzz_yzz[k] = -g_x_0_yyzzz_yzz[k] * ab_y + g_x_0_yyzzz_yyzz[k];

                g_x_0_yyyzzz_zzz[k] = -g_x_0_yyzzz_zzz[k] * ab_y + g_x_0_yyzzz_yzzz[k];
            }

            /// Set up 250-260 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_xxx = cbuffer.data(if_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxy = cbuffer.data(if_geom_10_off + 251 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxz = cbuffer.data(if_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyy = cbuffer.data(if_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyz = cbuffer.data(if_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xzz = cbuffer.data(if_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyy = cbuffer.data(if_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyz = cbuffer.data(if_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yzz = cbuffer.data(if_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_yyzzzz_zzz = cbuffer.data(if_geom_10_off + 259 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzzzz_xxx, g_x_0_yyzzzz_xxy, g_x_0_yyzzzz_xxz, g_x_0_yyzzzz_xyy, g_x_0_yyzzzz_xyz, g_x_0_yyzzzz_xzz, g_x_0_yyzzzz_yyy, g_x_0_yyzzzz_yyz, g_x_0_yyzzzz_yzz, g_x_0_yyzzzz_zzz, g_x_0_yzzzz_xxx, g_x_0_yzzzz_xxxy, g_x_0_yzzzz_xxy, g_x_0_yzzzz_xxyy, g_x_0_yzzzz_xxyz, g_x_0_yzzzz_xxz, g_x_0_yzzzz_xyy, g_x_0_yzzzz_xyyy, g_x_0_yzzzz_xyyz, g_x_0_yzzzz_xyz, g_x_0_yzzzz_xyzz, g_x_0_yzzzz_xzz, g_x_0_yzzzz_yyy, g_x_0_yzzzz_yyyy, g_x_0_yzzzz_yyyz, g_x_0_yzzzz_yyz, g_x_0_yzzzz_yyzz, g_x_0_yzzzz_yzz, g_x_0_yzzzz_yzzz, g_x_0_yzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_xxx[k] = -g_x_0_yzzzz_xxx[k] * ab_y + g_x_0_yzzzz_xxxy[k];

                g_x_0_yyzzzz_xxy[k] = -g_x_0_yzzzz_xxy[k] * ab_y + g_x_0_yzzzz_xxyy[k];

                g_x_0_yyzzzz_xxz[k] = -g_x_0_yzzzz_xxz[k] * ab_y + g_x_0_yzzzz_xxyz[k];

                g_x_0_yyzzzz_xyy[k] = -g_x_0_yzzzz_xyy[k] * ab_y + g_x_0_yzzzz_xyyy[k];

                g_x_0_yyzzzz_xyz[k] = -g_x_0_yzzzz_xyz[k] * ab_y + g_x_0_yzzzz_xyyz[k];

                g_x_0_yyzzzz_xzz[k] = -g_x_0_yzzzz_xzz[k] * ab_y + g_x_0_yzzzz_xyzz[k];

                g_x_0_yyzzzz_yyy[k] = -g_x_0_yzzzz_yyy[k] * ab_y + g_x_0_yzzzz_yyyy[k];

                g_x_0_yyzzzz_yyz[k] = -g_x_0_yzzzz_yyz[k] * ab_y + g_x_0_yzzzz_yyyz[k];

                g_x_0_yyzzzz_yzz[k] = -g_x_0_yzzzz_yzz[k] * ab_y + g_x_0_yzzzz_yyzz[k];

                g_x_0_yyzzzz_zzz[k] = -g_x_0_yzzzz_zzz[k] * ab_y + g_x_0_yzzzz_yzzz[k];
            }

            /// Set up 260-270 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_xxx = cbuffer.data(if_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxy = cbuffer.data(if_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxz = cbuffer.data(if_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyy = cbuffer.data(if_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyz = cbuffer.data(if_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xzz = cbuffer.data(if_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyy = cbuffer.data(if_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyz = cbuffer.data(if_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yzz = cbuffer.data(if_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_yzzzzz_zzz = cbuffer.data(if_geom_10_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzzzz_xxx, g_x_0_yzzzzz_xxy, g_x_0_yzzzzz_xxz, g_x_0_yzzzzz_xyy, g_x_0_yzzzzz_xyz, g_x_0_yzzzzz_xzz, g_x_0_yzzzzz_yyy, g_x_0_yzzzzz_yyz, g_x_0_yzzzzz_yzz, g_x_0_yzzzzz_zzz, g_x_0_zzzzz_xxx, g_x_0_zzzzz_xxxy, g_x_0_zzzzz_xxy, g_x_0_zzzzz_xxyy, g_x_0_zzzzz_xxyz, g_x_0_zzzzz_xxz, g_x_0_zzzzz_xyy, g_x_0_zzzzz_xyyy, g_x_0_zzzzz_xyyz, g_x_0_zzzzz_xyz, g_x_0_zzzzz_xyzz, g_x_0_zzzzz_xzz, g_x_0_zzzzz_yyy, g_x_0_zzzzz_yyyy, g_x_0_zzzzz_yyyz, g_x_0_zzzzz_yyz, g_x_0_zzzzz_yyzz, g_x_0_zzzzz_yzz, g_x_0_zzzzz_yzzz, g_x_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_xxx[k] = -g_x_0_zzzzz_xxx[k] * ab_y + g_x_0_zzzzz_xxxy[k];

                g_x_0_yzzzzz_xxy[k] = -g_x_0_zzzzz_xxy[k] * ab_y + g_x_0_zzzzz_xxyy[k];

                g_x_0_yzzzzz_xxz[k] = -g_x_0_zzzzz_xxz[k] * ab_y + g_x_0_zzzzz_xxyz[k];

                g_x_0_yzzzzz_xyy[k] = -g_x_0_zzzzz_xyy[k] * ab_y + g_x_0_zzzzz_xyyy[k];

                g_x_0_yzzzzz_xyz[k] = -g_x_0_zzzzz_xyz[k] * ab_y + g_x_0_zzzzz_xyyz[k];

                g_x_0_yzzzzz_xzz[k] = -g_x_0_zzzzz_xzz[k] * ab_y + g_x_0_zzzzz_xyzz[k];

                g_x_0_yzzzzz_yyy[k] = -g_x_0_zzzzz_yyy[k] * ab_y + g_x_0_zzzzz_yyyy[k];

                g_x_0_yzzzzz_yyz[k] = -g_x_0_zzzzz_yyz[k] * ab_y + g_x_0_zzzzz_yyyz[k];

                g_x_0_yzzzzz_yzz[k] = -g_x_0_zzzzz_yzz[k] * ab_y + g_x_0_zzzzz_yyzz[k];

                g_x_0_yzzzzz_zzz[k] = -g_x_0_zzzzz_zzz[k] * ab_y + g_x_0_zzzzz_yzzz[k];
            }

            /// Set up 270-280 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_xxx = cbuffer.data(if_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxy = cbuffer.data(if_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxz = cbuffer.data(if_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyy = cbuffer.data(if_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyz = cbuffer.data(if_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xzz = cbuffer.data(if_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyy = cbuffer.data(if_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyz = cbuffer.data(if_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yzz = cbuffer.data(if_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_zzzzzz_zzz = cbuffer.data(if_geom_10_off + 279 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzzz_xxx, g_x_0_zzzzz_xxxz, g_x_0_zzzzz_xxy, g_x_0_zzzzz_xxyz, g_x_0_zzzzz_xxz, g_x_0_zzzzz_xxzz, g_x_0_zzzzz_xyy, g_x_0_zzzzz_xyyz, g_x_0_zzzzz_xyz, g_x_0_zzzzz_xyzz, g_x_0_zzzzz_xzz, g_x_0_zzzzz_xzzz, g_x_0_zzzzz_yyy, g_x_0_zzzzz_yyyz, g_x_0_zzzzz_yyz, g_x_0_zzzzz_yyzz, g_x_0_zzzzz_yzz, g_x_0_zzzzz_yzzz, g_x_0_zzzzz_zzz, g_x_0_zzzzz_zzzz, g_x_0_zzzzzz_xxx, g_x_0_zzzzzz_xxy, g_x_0_zzzzzz_xxz, g_x_0_zzzzzz_xyy, g_x_0_zzzzzz_xyz, g_x_0_zzzzzz_xzz, g_x_0_zzzzzz_yyy, g_x_0_zzzzzz_yyz, g_x_0_zzzzzz_yzz, g_x_0_zzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_xxx[k] = -g_x_0_zzzzz_xxx[k] * ab_z + g_x_0_zzzzz_xxxz[k];

                g_x_0_zzzzzz_xxy[k] = -g_x_0_zzzzz_xxy[k] * ab_z + g_x_0_zzzzz_xxyz[k];

                g_x_0_zzzzzz_xxz[k] = -g_x_0_zzzzz_xxz[k] * ab_z + g_x_0_zzzzz_xxzz[k];

                g_x_0_zzzzzz_xyy[k] = -g_x_0_zzzzz_xyy[k] * ab_z + g_x_0_zzzzz_xyyz[k];

                g_x_0_zzzzzz_xyz[k] = -g_x_0_zzzzz_xyz[k] * ab_z + g_x_0_zzzzz_xyzz[k];

                g_x_0_zzzzzz_xzz[k] = -g_x_0_zzzzz_xzz[k] * ab_z + g_x_0_zzzzz_xzzz[k];

                g_x_0_zzzzzz_yyy[k] = -g_x_0_zzzzz_yyy[k] * ab_z + g_x_0_zzzzz_yyyz[k];

                g_x_0_zzzzzz_yyz[k] = -g_x_0_zzzzz_yyz[k] * ab_z + g_x_0_zzzzz_yyzz[k];

                g_x_0_zzzzzz_yzz[k] = -g_x_0_zzzzz_yzz[k] * ab_z + g_x_0_zzzzz_yzzz[k];

                g_x_0_zzzzzz_zzz[k] = -g_x_0_zzzzz_zzz[k] * ab_z + g_x_0_zzzzz_zzzz[k];
            }

            /// Set up 280-290 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_xxx = cbuffer.data(if_geom_10_off + 280 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxy = cbuffer.data(if_geom_10_off + 281 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxz = cbuffer.data(if_geom_10_off + 282 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyy = cbuffer.data(if_geom_10_off + 283 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyz = cbuffer.data(if_geom_10_off + 284 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xzz = cbuffer.data(if_geom_10_off + 285 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyy = cbuffer.data(if_geom_10_off + 286 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyz = cbuffer.data(if_geom_10_off + 287 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yzz = cbuffer.data(if_geom_10_off + 288 * ccomps * dcomps);

            auto g_y_0_xxxxxx_zzz = cbuffer.data(if_geom_10_off + 289 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxx_xxx, g_y_0_xxxxx_xxxx, g_y_0_xxxxx_xxxy, g_y_0_xxxxx_xxxz, g_y_0_xxxxx_xxy, g_y_0_xxxxx_xxyy, g_y_0_xxxxx_xxyz, g_y_0_xxxxx_xxz, g_y_0_xxxxx_xxzz, g_y_0_xxxxx_xyy, g_y_0_xxxxx_xyyy, g_y_0_xxxxx_xyyz, g_y_0_xxxxx_xyz, g_y_0_xxxxx_xyzz, g_y_0_xxxxx_xzz, g_y_0_xxxxx_xzzz, g_y_0_xxxxx_yyy, g_y_0_xxxxx_yyz, g_y_0_xxxxx_yzz, g_y_0_xxxxx_zzz, g_y_0_xxxxxx_xxx, g_y_0_xxxxxx_xxy, g_y_0_xxxxxx_xxz, g_y_0_xxxxxx_xyy, g_y_0_xxxxxx_xyz, g_y_0_xxxxxx_xzz, g_y_0_xxxxxx_yyy, g_y_0_xxxxxx_yyz, g_y_0_xxxxxx_yzz, g_y_0_xxxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_xxx[k] = -g_y_0_xxxxx_xxx[k] * ab_x + g_y_0_xxxxx_xxxx[k];

                g_y_0_xxxxxx_xxy[k] = -g_y_0_xxxxx_xxy[k] * ab_x + g_y_0_xxxxx_xxxy[k];

                g_y_0_xxxxxx_xxz[k] = -g_y_0_xxxxx_xxz[k] * ab_x + g_y_0_xxxxx_xxxz[k];

                g_y_0_xxxxxx_xyy[k] = -g_y_0_xxxxx_xyy[k] * ab_x + g_y_0_xxxxx_xxyy[k];

                g_y_0_xxxxxx_xyz[k] = -g_y_0_xxxxx_xyz[k] * ab_x + g_y_0_xxxxx_xxyz[k];

                g_y_0_xxxxxx_xzz[k] = -g_y_0_xxxxx_xzz[k] * ab_x + g_y_0_xxxxx_xxzz[k];

                g_y_0_xxxxxx_yyy[k] = -g_y_0_xxxxx_yyy[k] * ab_x + g_y_0_xxxxx_xyyy[k];

                g_y_0_xxxxxx_yyz[k] = -g_y_0_xxxxx_yyz[k] * ab_x + g_y_0_xxxxx_xyyz[k];

                g_y_0_xxxxxx_yzz[k] = -g_y_0_xxxxx_yzz[k] * ab_x + g_y_0_xxxxx_xyzz[k];

                g_y_0_xxxxxx_zzz[k] = -g_y_0_xxxxx_zzz[k] * ab_x + g_y_0_xxxxx_xzzz[k];
            }

            /// Set up 290-300 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_xxx = cbuffer.data(if_geom_10_off + 290 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxy = cbuffer.data(if_geom_10_off + 291 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxz = cbuffer.data(if_geom_10_off + 292 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyy = cbuffer.data(if_geom_10_off + 293 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyz = cbuffer.data(if_geom_10_off + 294 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xzz = cbuffer.data(if_geom_10_off + 295 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyy = cbuffer.data(if_geom_10_off + 296 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyz = cbuffer.data(if_geom_10_off + 297 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yzz = cbuffer.data(if_geom_10_off + 298 * ccomps * dcomps);

            auto g_y_0_xxxxxy_zzz = cbuffer.data(if_geom_10_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxxy_xxx, g_y_0_xxxxxy_xxy, g_y_0_xxxxxy_xxz, g_y_0_xxxxxy_xyy, g_y_0_xxxxxy_xyz, g_y_0_xxxxxy_xzz, g_y_0_xxxxxy_yyy, g_y_0_xxxxxy_yyz, g_y_0_xxxxxy_yzz, g_y_0_xxxxxy_zzz, g_y_0_xxxxy_xxx, g_y_0_xxxxy_xxxx, g_y_0_xxxxy_xxxy, g_y_0_xxxxy_xxxz, g_y_0_xxxxy_xxy, g_y_0_xxxxy_xxyy, g_y_0_xxxxy_xxyz, g_y_0_xxxxy_xxz, g_y_0_xxxxy_xxzz, g_y_0_xxxxy_xyy, g_y_0_xxxxy_xyyy, g_y_0_xxxxy_xyyz, g_y_0_xxxxy_xyz, g_y_0_xxxxy_xyzz, g_y_0_xxxxy_xzz, g_y_0_xxxxy_xzzz, g_y_0_xxxxy_yyy, g_y_0_xxxxy_yyz, g_y_0_xxxxy_yzz, g_y_0_xxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_xxx[k] = -g_y_0_xxxxy_xxx[k] * ab_x + g_y_0_xxxxy_xxxx[k];

                g_y_0_xxxxxy_xxy[k] = -g_y_0_xxxxy_xxy[k] * ab_x + g_y_0_xxxxy_xxxy[k];

                g_y_0_xxxxxy_xxz[k] = -g_y_0_xxxxy_xxz[k] * ab_x + g_y_0_xxxxy_xxxz[k];

                g_y_0_xxxxxy_xyy[k] = -g_y_0_xxxxy_xyy[k] * ab_x + g_y_0_xxxxy_xxyy[k];

                g_y_0_xxxxxy_xyz[k] = -g_y_0_xxxxy_xyz[k] * ab_x + g_y_0_xxxxy_xxyz[k];

                g_y_0_xxxxxy_xzz[k] = -g_y_0_xxxxy_xzz[k] * ab_x + g_y_0_xxxxy_xxzz[k];

                g_y_0_xxxxxy_yyy[k] = -g_y_0_xxxxy_yyy[k] * ab_x + g_y_0_xxxxy_xyyy[k];

                g_y_0_xxxxxy_yyz[k] = -g_y_0_xxxxy_yyz[k] * ab_x + g_y_0_xxxxy_xyyz[k];

                g_y_0_xxxxxy_yzz[k] = -g_y_0_xxxxy_yzz[k] * ab_x + g_y_0_xxxxy_xyzz[k];

                g_y_0_xxxxxy_zzz[k] = -g_y_0_xxxxy_zzz[k] * ab_x + g_y_0_xxxxy_xzzz[k];
            }

            /// Set up 300-310 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_xxx = cbuffer.data(if_geom_10_off + 300 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxy = cbuffer.data(if_geom_10_off + 301 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxz = cbuffer.data(if_geom_10_off + 302 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyy = cbuffer.data(if_geom_10_off + 303 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyz = cbuffer.data(if_geom_10_off + 304 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xzz = cbuffer.data(if_geom_10_off + 305 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyy = cbuffer.data(if_geom_10_off + 306 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyz = cbuffer.data(if_geom_10_off + 307 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yzz = cbuffer.data(if_geom_10_off + 308 * ccomps * dcomps);

            auto g_y_0_xxxxxz_zzz = cbuffer.data(if_geom_10_off + 309 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxxz_xxx, g_y_0_xxxxxz_xxy, g_y_0_xxxxxz_xxz, g_y_0_xxxxxz_xyy, g_y_0_xxxxxz_xyz, g_y_0_xxxxxz_xzz, g_y_0_xxxxxz_yyy, g_y_0_xxxxxz_yyz, g_y_0_xxxxxz_yzz, g_y_0_xxxxxz_zzz, g_y_0_xxxxz_xxx, g_y_0_xxxxz_xxxx, g_y_0_xxxxz_xxxy, g_y_0_xxxxz_xxxz, g_y_0_xxxxz_xxy, g_y_0_xxxxz_xxyy, g_y_0_xxxxz_xxyz, g_y_0_xxxxz_xxz, g_y_0_xxxxz_xxzz, g_y_0_xxxxz_xyy, g_y_0_xxxxz_xyyy, g_y_0_xxxxz_xyyz, g_y_0_xxxxz_xyz, g_y_0_xxxxz_xyzz, g_y_0_xxxxz_xzz, g_y_0_xxxxz_xzzz, g_y_0_xxxxz_yyy, g_y_0_xxxxz_yyz, g_y_0_xxxxz_yzz, g_y_0_xxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_xxx[k] = -g_y_0_xxxxz_xxx[k] * ab_x + g_y_0_xxxxz_xxxx[k];

                g_y_0_xxxxxz_xxy[k] = -g_y_0_xxxxz_xxy[k] * ab_x + g_y_0_xxxxz_xxxy[k];

                g_y_0_xxxxxz_xxz[k] = -g_y_0_xxxxz_xxz[k] * ab_x + g_y_0_xxxxz_xxxz[k];

                g_y_0_xxxxxz_xyy[k] = -g_y_0_xxxxz_xyy[k] * ab_x + g_y_0_xxxxz_xxyy[k];

                g_y_0_xxxxxz_xyz[k] = -g_y_0_xxxxz_xyz[k] * ab_x + g_y_0_xxxxz_xxyz[k];

                g_y_0_xxxxxz_xzz[k] = -g_y_0_xxxxz_xzz[k] * ab_x + g_y_0_xxxxz_xxzz[k];

                g_y_0_xxxxxz_yyy[k] = -g_y_0_xxxxz_yyy[k] * ab_x + g_y_0_xxxxz_xyyy[k];

                g_y_0_xxxxxz_yyz[k] = -g_y_0_xxxxz_yyz[k] * ab_x + g_y_0_xxxxz_xyyz[k];

                g_y_0_xxxxxz_yzz[k] = -g_y_0_xxxxz_yzz[k] * ab_x + g_y_0_xxxxz_xyzz[k];

                g_y_0_xxxxxz_zzz[k] = -g_y_0_xxxxz_zzz[k] * ab_x + g_y_0_xxxxz_xzzz[k];
            }

            /// Set up 310-320 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_xxx = cbuffer.data(if_geom_10_off + 310 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxy = cbuffer.data(if_geom_10_off + 311 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxz = cbuffer.data(if_geom_10_off + 312 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyy = cbuffer.data(if_geom_10_off + 313 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyz = cbuffer.data(if_geom_10_off + 314 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xzz = cbuffer.data(if_geom_10_off + 315 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyy = cbuffer.data(if_geom_10_off + 316 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyz = cbuffer.data(if_geom_10_off + 317 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yzz = cbuffer.data(if_geom_10_off + 318 * ccomps * dcomps);

            auto g_y_0_xxxxyy_zzz = cbuffer.data(if_geom_10_off + 319 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxyy_xxx, g_y_0_xxxxyy_xxy, g_y_0_xxxxyy_xxz, g_y_0_xxxxyy_xyy, g_y_0_xxxxyy_xyz, g_y_0_xxxxyy_xzz, g_y_0_xxxxyy_yyy, g_y_0_xxxxyy_yyz, g_y_0_xxxxyy_yzz, g_y_0_xxxxyy_zzz, g_y_0_xxxyy_xxx, g_y_0_xxxyy_xxxx, g_y_0_xxxyy_xxxy, g_y_0_xxxyy_xxxz, g_y_0_xxxyy_xxy, g_y_0_xxxyy_xxyy, g_y_0_xxxyy_xxyz, g_y_0_xxxyy_xxz, g_y_0_xxxyy_xxzz, g_y_0_xxxyy_xyy, g_y_0_xxxyy_xyyy, g_y_0_xxxyy_xyyz, g_y_0_xxxyy_xyz, g_y_0_xxxyy_xyzz, g_y_0_xxxyy_xzz, g_y_0_xxxyy_xzzz, g_y_0_xxxyy_yyy, g_y_0_xxxyy_yyz, g_y_0_xxxyy_yzz, g_y_0_xxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_xxx[k] = -g_y_0_xxxyy_xxx[k] * ab_x + g_y_0_xxxyy_xxxx[k];

                g_y_0_xxxxyy_xxy[k] = -g_y_0_xxxyy_xxy[k] * ab_x + g_y_0_xxxyy_xxxy[k];

                g_y_0_xxxxyy_xxz[k] = -g_y_0_xxxyy_xxz[k] * ab_x + g_y_0_xxxyy_xxxz[k];

                g_y_0_xxxxyy_xyy[k] = -g_y_0_xxxyy_xyy[k] * ab_x + g_y_0_xxxyy_xxyy[k];

                g_y_0_xxxxyy_xyz[k] = -g_y_0_xxxyy_xyz[k] * ab_x + g_y_0_xxxyy_xxyz[k];

                g_y_0_xxxxyy_xzz[k] = -g_y_0_xxxyy_xzz[k] * ab_x + g_y_0_xxxyy_xxzz[k];

                g_y_0_xxxxyy_yyy[k] = -g_y_0_xxxyy_yyy[k] * ab_x + g_y_0_xxxyy_xyyy[k];

                g_y_0_xxxxyy_yyz[k] = -g_y_0_xxxyy_yyz[k] * ab_x + g_y_0_xxxyy_xyyz[k];

                g_y_0_xxxxyy_yzz[k] = -g_y_0_xxxyy_yzz[k] * ab_x + g_y_0_xxxyy_xyzz[k];

                g_y_0_xxxxyy_zzz[k] = -g_y_0_xxxyy_zzz[k] * ab_x + g_y_0_xxxyy_xzzz[k];
            }

            /// Set up 320-330 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_xxx = cbuffer.data(if_geom_10_off + 320 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxy = cbuffer.data(if_geom_10_off + 321 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxz = cbuffer.data(if_geom_10_off + 322 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyy = cbuffer.data(if_geom_10_off + 323 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyz = cbuffer.data(if_geom_10_off + 324 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xzz = cbuffer.data(if_geom_10_off + 325 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyy = cbuffer.data(if_geom_10_off + 326 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyz = cbuffer.data(if_geom_10_off + 327 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yzz = cbuffer.data(if_geom_10_off + 328 * ccomps * dcomps);

            auto g_y_0_xxxxyz_zzz = cbuffer.data(if_geom_10_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxyz_xxx, g_y_0_xxxxyz_xxy, g_y_0_xxxxyz_xxz, g_y_0_xxxxyz_xyy, g_y_0_xxxxyz_xyz, g_y_0_xxxxyz_xzz, g_y_0_xxxxyz_yyy, g_y_0_xxxxyz_yyz, g_y_0_xxxxyz_yzz, g_y_0_xxxxyz_zzz, g_y_0_xxxyz_xxx, g_y_0_xxxyz_xxxx, g_y_0_xxxyz_xxxy, g_y_0_xxxyz_xxxz, g_y_0_xxxyz_xxy, g_y_0_xxxyz_xxyy, g_y_0_xxxyz_xxyz, g_y_0_xxxyz_xxz, g_y_0_xxxyz_xxzz, g_y_0_xxxyz_xyy, g_y_0_xxxyz_xyyy, g_y_0_xxxyz_xyyz, g_y_0_xxxyz_xyz, g_y_0_xxxyz_xyzz, g_y_0_xxxyz_xzz, g_y_0_xxxyz_xzzz, g_y_0_xxxyz_yyy, g_y_0_xxxyz_yyz, g_y_0_xxxyz_yzz, g_y_0_xxxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_xxx[k] = -g_y_0_xxxyz_xxx[k] * ab_x + g_y_0_xxxyz_xxxx[k];

                g_y_0_xxxxyz_xxy[k] = -g_y_0_xxxyz_xxy[k] * ab_x + g_y_0_xxxyz_xxxy[k];

                g_y_0_xxxxyz_xxz[k] = -g_y_0_xxxyz_xxz[k] * ab_x + g_y_0_xxxyz_xxxz[k];

                g_y_0_xxxxyz_xyy[k] = -g_y_0_xxxyz_xyy[k] * ab_x + g_y_0_xxxyz_xxyy[k];

                g_y_0_xxxxyz_xyz[k] = -g_y_0_xxxyz_xyz[k] * ab_x + g_y_0_xxxyz_xxyz[k];

                g_y_0_xxxxyz_xzz[k] = -g_y_0_xxxyz_xzz[k] * ab_x + g_y_0_xxxyz_xxzz[k];

                g_y_0_xxxxyz_yyy[k] = -g_y_0_xxxyz_yyy[k] * ab_x + g_y_0_xxxyz_xyyy[k];

                g_y_0_xxxxyz_yyz[k] = -g_y_0_xxxyz_yyz[k] * ab_x + g_y_0_xxxyz_xyyz[k];

                g_y_0_xxxxyz_yzz[k] = -g_y_0_xxxyz_yzz[k] * ab_x + g_y_0_xxxyz_xyzz[k];

                g_y_0_xxxxyz_zzz[k] = -g_y_0_xxxyz_zzz[k] * ab_x + g_y_0_xxxyz_xzzz[k];
            }

            /// Set up 330-340 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_xxx = cbuffer.data(if_geom_10_off + 330 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxy = cbuffer.data(if_geom_10_off + 331 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxz = cbuffer.data(if_geom_10_off + 332 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyy = cbuffer.data(if_geom_10_off + 333 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyz = cbuffer.data(if_geom_10_off + 334 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xzz = cbuffer.data(if_geom_10_off + 335 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyy = cbuffer.data(if_geom_10_off + 336 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyz = cbuffer.data(if_geom_10_off + 337 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yzz = cbuffer.data(if_geom_10_off + 338 * ccomps * dcomps);

            auto g_y_0_xxxxzz_zzz = cbuffer.data(if_geom_10_off + 339 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxzz_xxx, g_y_0_xxxxzz_xxy, g_y_0_xxxxzz_xxz, g_y_0_xxxxzz_xyy, g_y_0_xxxxzz_xyz, g_y_0_xxxxzz_xzz, g_y_0_xxxxzz_yyy, g_y_0_xxxxzz_yyz, g_y_0_xxxxzz_yzz, g_y_0_xxxxzz_zzz, g_y_0_xxxzz_xxx, g_y_0_xxxzz_xxxx, g_y_0_xxxzz_xxxy, g_y_0_xxxzz_xxxz, g_y_0_xxxzz_xxy, g_y_0_xxxzz_xxyy, g_y_0_xxxzz_xxyz, g_y_0_xxxzz_xxz, g_y_0_xxxzz_xxzz, g_y_0_xxxzz_xyy, g_y_0_xxxzz_xyyy, g_y_0_xxxzz_xyyz, g_y_0_xxxzz_xyz, g_y_0_xxxzz_xyzz, g_y_0_xxxzz_xzz, g_y_0_xxxzz_xzzz, g_y_0_xxxzz_yyy, g_y_0_xxxzz_yyz, g_y_0_xxxzz_yzz, g_y_0_xxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_xxx[k] = -g_y_0_xxxzz_xxx[k] * ab_x + g_y_0_xxxzz_xxxx[k];

                g_y_0_xxxxzz_xxy[k] = -g_y_0_xxxzz_xxy[k] * ab_x + g_y_0_xxxzz_xxxy[k];

                g_y_0_xxxxzz_xxz[k] = -g_y_0_xxxzz_xxz[k] * ab_x + g_y_0_xxxzz_xxxz[k];

                g_y_0_xxxxzz_xyy[k] = -g_y_0_xxxzz_xyy[k] * ab_x + g_y_0_xxxzz_xxyy[k];

                g_y_0_xxxxzz_xyz[k] = -g_y_0_xxxzz_xyz[k] * ab_x + g_y_0_xxxzz_xxyz[k];

                g_y_0_xxxxzz_xzz[k] = -g_y_0_xxxzz_xzz[k] * ab_x + g_y_0_xxxzz_xxzz[k];

                g_y_0_xxxxzz_yyy[k] = -g_y_0_xxxzz_yyy[k] * ab_x + g_y_0_xxxzz_xyyy[k];

                g_y_0_xxxxzz_yyz[k] = -g_y_0_xxxzz_yyz[k] * ab_x + g_y_0_xxxzz_xyyz[k];

                g_y_0_xxxxzz_yzz[k] = -g_y_0_xxxzz_yzz[k] * ab_x + g_y_0_xxxzz_xyzz[k];

                g_y_0_xxxxzz_zzz[k] = -g_y_0_xxxzz_zzz[k] * ab_x + g_y_0_xxxzz_xzzz[k];
            }

            /// Set up 340-350 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_xxx = cbuffer.data(if_geom_10_off + 340 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxy = cbuffer.data(if_geom_10_off + 341 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxz = cbuffer.data(if_geom_10_off + 342 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyy = cbuffer.data(if_geom_10_off + 343 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyz = cbuffer.data(if_geom_10_off + 344 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xzz = cbuffer.data(if_geom_10_off + 345 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyy = cbuffer.data(if_geom_10_off + 346 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyz = cbuffer.data(if_geom_10_off + 347 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yzz = cbuffer.data(if_geom_10_off + 348 * ccomps * dcomps);

            auto g_y_0_xxxyyy_zzz = cbuffer.data(if_geom_10_off + 349 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyyy_xxx, g_y_0_xxxyyy_xxy, g_y_0_xxxyyy_xxz, g_y_0_xxxyyy_xyy, g_y_0_xxxyyy_xyz, g_y_0_xxxyyy_xzz, g_y_0_xxxyyy_yyy, g_y_0_xxxyyy_yyz, g_y_0_xxxyyy_yzz, g_y_0_xxxyyy_zzz, g_y_0_xxyyy_xxx, g_y_0_xxyyy_xxxx, g_y_0_xxyyy_xxxy, g_y_0_xxyyy_xxxz, g_y_0_xxyyy_xxy, g_y_0_xxyyy_xxyy, g_y_0_xxyyy_xxyz, g_y_0_xxyyy_xxz, g_y_0_xxyyy_xxzz, g_y_0_xxyyy_xyy, g_y_0_xxyyy_xyyy, g_y_0_xxyyy_xyyz, g_y_0_xxyyy_xyz, g_y_0_xxyyy_xyzz, g_y_0_xxyyy_xzz, g_y_0_xxyyy_xzzz, g_y_0_xxyyy_yyy, g_y_0_xxyyy_yyz, g_y_0_xxyyy_yzz, g_y_0_xxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_xxx[k] = -g_y_0_xxyyy_xxx[k] * ab_x + g_y_0_xxyyy_xxxx[k];

                g_y_0_xxxyyy_xxy[k] = -g_y_0_xxyyy_xxy[k] * ab_x + g_y_0_xxyyy_xxxy[k];

                g_y_0_xxxyyy_xxz[k] = -g_y_0_xxyyy_xxz[k] * ab_x + g_y_0_xxyyy_xxxz[k];

                g_y_0_xxxyyy_xyy[k] = -g_y_0_xxyyy_xyy[k] * ab_x + g_y_0_xxyyy_xxyy[k];

                g_y_0_xxxyyy_xyz[k] = -g_y_0_xxyyy_xyz[k] * ab_x + g_y_0_xxyyy_xxyz[k];

                g_y_0_xxxyyy_xzz[k] = -g_y_0_xxyyy_xzz[k] * ab_x + g_y_0_xxyyy_xxzz[k];

                g_y_0_xxxyyy_yyy[k] = -g_y_0_xxyyy_yyy[k] * ab_x + g_y_0_xxyyy_xyyy[k];

                g_y_0_xxxyyy_yyz[k] = -g_y_0_xxyyy_yyz[k] * ab_x + g_y_0_xxyyy_xyyz[k];

                g_y_0_xxxyyy_yzz[k] = -g_y_0_xxyyy_yzz[k] * ab_x + g_y_0_xxyyy_xyzz[k];

                g_y_0_xxxyyy_zzz[k] = -g_y_0_xxyyy_zzz[k] * ab_x + g_y_0_xxyyy_xzzz[k];
            }

            /// Set up 350-360 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_xxx = cbuffer.data(if_geom_10_off + 350 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxy = cbuffer.data(if_geom_10_off + 351 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxz = cbuffer.data(if_geom_10_off + 352 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyy = cbuffer.data(if_geom_10_off + 353 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyz = cbuffer.data(if_geom_10_off + 354 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xzz = cbuffer.data(if_geom_10_off + 355 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyy = cbuffer.data(if_geom_10_off + 356 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyz = cbuffer.data(if_geom_10_off + 357 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yzz = cbuffer.data(if_geom_10_off + 358 * ccomps * dcomps);

            auto g_y_0_xxxyyz_zzz = cbuffer.data(if_geom_10_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyyz_xxx, g_y_0_xxxyyz_xxy, g_y_0_xxxyyz_xxz, g_y_0_xxxyyz_xyy, g_y_0_xxxyyz_xyz, g_y_0_xxxyyz_xzz, g_y_0_xxxyyz_yyy, g_y_0_xxxyyz_yyz, g_y_0_xxxyyz_yzz, g_y_0_xxxyyz_zzz, g_y_0_xxyyz_xxx, g_y_0_xxyyz_xxxx, g_y_0_xxyyz_xxxy, g_y_0_xxyyz_xxxz, g_y_0_xxyyz_xxy, g_y_0_xxyyz_xxyy, g_y_0_xxyyz_xxyz, g_y_0_xxyyz_xxz, g_y_0_xxyyz_xxzz, g_y_0_xxyyz_xyy, g_y_0_xxyyz_xyyy, g_y_0_xxyyz_xyyz, g_y_0_xxyyz_xyz, g_y_0_xxyyz_xyzz, g_y_0_xxyyz_xzz, g_y_0_xxyyz_xzzz, g_y_0_xxyyz_yyy, g_y_0_xxyyz_yyz, g_y_0_xxyyz_yzz, g_y_0_xxyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_xxx[k] = -g_y_0_xxyyz_xxx[k] * ab_x + g_y_0_xxyyz_xxxx[k];

                g_y_0_xxxyyz_xxy[k] = -g_y_0_xxyyz_xxy[k] * ab_x + g_y_0_xxyyz_xxxy[k];

                g_y_0_xxxyyz_xxz[k] = -g_y_0_xxyyz_xxz[k] * ab_x + g_y_0_xxyyz_xxxz[k];

                g_y_0_xxxyyz_xyy[k] = -g_y_0_xxyyz_xyy[k] * ab_x + g_y_0_xxyyz_xxyy[k];

                g_y_0_xxxyyz_xyz[k] = -g_y_0_xxyyz_xyz[k] * ab_x + g_y_0_xxyyz_xxyz[k];

                g_y_0_xxxyyz_xzz[k] = -g_y_0_xxyyz_xzz[k] * ab_x + g_y_0_xxyyz_xxzz[k];

                g_y_0_xxxyyz_yyy[k] = -g_y_0_xxyyz_yyy[k] * ab_x + g_y_0_xxyyz_xyyy[k];

                g_y_0_xxxyyz_yyz[k] = -g_y_0_xxyyz_yyz[k] * ab_x + g_y_0_xxyyz_xyyz[k];

                g_y_0_xxxyyz_yzz[k] = -g_y_0_xxyyz_yzz[k] * ab_x + g_y_0_xxyyz_xyzz[k];

                g_y_0_xxxyyz_zzz[k] = -g_y_0_xxyyz_zzz[k] * ab_x + g_y_0_xxyyz_xzzz[k];
            }

            /// Set up 360-370 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_xxx = cbuffer.data(if_geom_10_off + 360 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxy = cbuffer.data(if_geom_10_off + 361 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxz = cbuffer.data(if_geom_10_off + 362 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyy = cbuffer.data(if_geom_10_off + 363 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyz = cbuffer.data(if_geom_10_off + 364 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xzz = cbuffer.data(if_geom_10_off + 365 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyy = cbuffer.data(if_geom_10_off + 366 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyz = cbuffer.data(if_geom_10_off + 367 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yzz = cbuffer.data(if_geom_10_off + 368 * ccomps * dcomps);

            auto g_y_0_xxxyzz_zzz = cbuffer.data(if_geom_10_off + 369 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyzz_xxx, g_y_0_xxxyzz_xxy, g_y_0_xxxyzz_xxz, g_y_0_xxxyzz_xyy, g_y_0_xxxyzz_xyz, g_y_0_xxxyzz_xzz, g_y_0_xxxyzz_yyy, g_y_0_xxxyzz_yyz, g_y_0_xxxyzz_yzz, g_y_0_xxxyzz_zzz, g_y_0_xxyzz_xxx, g_y_0_xxyzz_xxxx, g_y_0_xxyzz_xxxy, g_y_0_xxyzz_xxxz, g_y_0_xxyzz_xxy, g_y_0_xxyzz_xxyy, g_y_0_xxyzz_xxyz, g_y_0_xxyzz_xxz, g_y_0_xxyzz_xxzz, g_y_0_xxyzz_xyy, g_y_0_xxyzz_xyyy, g_y_0_xxyzz_xyyz, g_y_0_xxyzz_xyz, g_y_0_xxyzz_xyzz, g_y_0_xxyzz_xzz, g_y_0_xxyzz_xzzz, g_y_0_xxyzz_yyy, g_y_0_xxyzz_yyz, g_y_0_xxyzz_yzz, g_y_0_xxyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_xxx[k] = -g_y_0_xxyzz_xxx[k] * ab_x + g_y_0_xxyzz_xxxx[k];

                g_y_0_xxxyzz_xxy[k] = -g_y_0_xxyzz_xxy[k] * ab_x + g_y_0_xxyzz_xxxy[k];

                g_y_0_xxxyzz_xxz[k] = -g_y_0_xxyzz_xxz[k] * ab_x + g_y_0_xxyzz_xxxz[k];

                g_y_0_xxxyzz_xyy[k] = -g_y_0_xxyzz_xyy[k] * ab_x + g_y_0_xxyzz_xxyy[k];

                g_y_0_xxxyzz_xyz[k] = -g_y_0_xxyzz_xyz[k] * ab_x + g_y_0_xxyzz_xxyz[k];

                g_y_0_xxxyzz_xzz[k] = -g_y_0_xxyzz_xzz[k] * ab_x + g_y_0_xxyzz_xxzz[k];

                g_y_0_xxxyzz_yyy[k] = -g_y_0_xxyzz_yyy[k] * ab_x + g_y_0_xxyzz_xyyy[k];

                g_y_0_xxxyzz_yyz[k] = -g_y_0_xxyzz_yyz[k] * ab_x + g_y_0_xxyzz_xyyz[k];

                g_y_0_xxxyzz_yzz[k] = -g_y_0_xxyzz_yzz[k] * ab_x + g_y_0_xxyzz_xyzz[k];

                g_y_0_xxxyzz_zzz[k] = -g_y_0_xxyzz_zzz[k] * ab_x + g_y_0_xxyzz_xzzz[k];
            }

            /// Set up 370-380 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_xxx = cbuffer.data(if_geom_10_off + 370 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxy = cbuffer.data(if_geom_10_off + 371 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxz = cbuffer.data(if_geom_10_off + 372 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyy = cbuffer.data(if_geom_10_off + 373 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyz = cbuffer.data(if_geom_10_off + 374 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xzz = cbuffer.data(if_geom_10_off + 375 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyy = cbuffer.data(if_geom_10_off + 376 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyz = cbuffer.data(if_geom_10_off + 377 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yzz = cbuffer.data(if_geom_10_off + 378 * ccomps * dcomps);

            auto g_y_0_xxxzzz_zzz = cbuffer.data(if_geom_10_off + 379 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxzzz_xxx, g_y_0_xxxzzz_xxy, g_y_0_xxxzzz_xxz, g_y_0_xxxzzz_xyy, g_y_0_xxxzzz_xyz, g_y_0_xxxzzz_xzz, g_y_0_xxxzzz_yyy, g_y_0_xxxzzz_yyz, g_y_0_xxxzzz_yzz, g_y_0_xxxzzz_zzz, g_y_0_xxzzz_xxx, g_y_0_xxzzz_xxxx, g_y_0_xxzzz_xxxy, g_y_0_xxzzz_xxxz, g_y_0_xxzzz_xxy, g_y_0_xxzzz_xxyy, g_y_0_xxzzz_xxyz, g_y_0_xxzzz_xxz, g_y_0_xxzzz_xxzz, g_y_0_xxzzz_xyy, g_y_0_xxzzz_xyyy, g_y_0_xxzzz_xyyz, g_y_0_xxzzz_xyz, g_y_0_xxzzz_xyzz, g_y_0_xxzzz_xzz, g_y_0_xxzzz_xzzz, g_y_0_xxzzz_yyy, g_y_0_xxzzz_yyz, g_y_0_xxzzz_yzz, g_y_0_xxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_xxx[k] = -g_y_0_xxzzz_xxx[k] * ab_x + g_y_0_xxzzz_xxxx[k];

                g_y_0_xxxzzz_xxy[k] = -g_y_0_xxzzz_xxy[k] * ab_x + g_y_0_xxzzz_xxxy[k];

                g_y_0_xxxzzz_xxz[k] = -g_y_0_xxzzz_xxz[k] * ab_x + g_y_0_xxzzz_xxxz[k];

                g_y_0_xxxzzz_xyy[k] = -g_y_0_xxzzz_xyy[k] * ab_x + g_y_0_xxzzz_xxyy[k];

                g_y_0_xxxzzz_xyz[k] = -g_y_0_xxzzz_xyz[k] * ab_x + g_y_0_xxzzz_xxyz[k];

                g_y_0_xxxzzz_xzz[k] = -g_y_0_xxzzz_xzz[k] * ab_x + g_y_0_xxzzz_xxzz[k];

                g_y_0_xxxzzz_yyy[k] = -g_y_0_xxzzz_yyy[k] * ab_x + g_y_0_xxzzz_xyyy[k];

                g_y_0_xxxzzz_yyz[k] = -g_y_0_xxzzz_yyz[k] * ab_x + g_y_0_xxzzz_xyyz[k];

                g_y_0_xxxzzz_yzz[k] = -g_y_0_xxzzz_yzz[k] * ab_x + g_y_0_xxzzz_xyzz[k];

                g_y_0_xxxzzz_zzz[k] = -g_y_0_xxzzz_zzz[k] * ab_x + g_y_0_xxzzz_xzzz[k];
            }

            /// Set up 380-390 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_xxx = cbuffer.data(if_geom_10_off + 380 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxy = cbuffer.data(if_geom_10_off + 381 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxz = cbuffer.data(if_geom_10_off + 382 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyy = cbuffer.data(if_geom_10_off + 383 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyz = cbuffer.data(if_geom_10_off + 384 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xzz = cbuffer.data(if_geom_10_off + 385 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyy = cbuffer.data(if_geom_10_off + 386 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyz = cbuffer.data(if_geom_10_off + 387 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yzz = cbuffer.data(if_geom_10_off + 388 * ccomps * dcomps);

            auto g_y_0_xxyyyy_zzz = cbuffer.data(if_geom_10_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyyy_xxx, g_y_0_xxyyyy_xxy, g_y_0_xxyyyy_xxz, g_y_0_xxyyyy_xyy, g_y_0_xxyyyy_xyz, g_y_0_xxyyyy_xzz, g_y_0_xxyyyy_yyy, g_y_0_xxyyyy_yyz, g_y_0_xxyyyy_yzz, g_y_0_xxyyyy_zzz, g_y_0_xyyyy_xxx, g_y_0_xyyyy_xxxx, g_y_0_xyyyy_xxxy, g_y_0_xyyyy_xxxz, g_y_0_xyyyy_xxy, g_y_0_xyyyy_xxyy, g_y_0_xyyyy_xxyz, g_y_0_xyyyy_xxz, g_y_0_xyyyy_xxzz, g_y_0_xyyyy_xyy, g_y_0_xyyyy_xyyy, g_y_0_xyyyy_xyyz, g_y_0_xyyyy_xyz, g_y_0_xyyyy_xyzz, g_y_0_xyyyy_xzz, g_y_0_xyyyy_xzzz, g_y_0_xyyyy_yyy, g_y_0_xyyyy_yyz, g_y_0_xyyyy_yzz, g_y_0_xyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_xxx[k] = -g_y_0_xyyyy_xxx[k] * ab_x + g_y_0_xyyyy_xxxx[k];

                g_y_0_xxyyyy_xxy[k] = -g_y_0_xyyyy_xxy[k] * ab_x + g_y_0_xyyyy_xxxy[k];

                g_y_0_xxyyyy_xxz[k] = -g_y_0_xyyyy_xxz[k] * ab_x + g_y_0_xyyyy_xxxz[k];

                g_y_0_xxyyyy_xyy[k] = -g_y_0_xyyyy_xyy[k] * ab_x + g_y_0_xyyyy_xxyy[k];

                g_y_0_xxyyyy_xyz[k] = -g_y_0_xyyyy_xyz[k] * ab_x + g_y_0_xyyyy_xxyz[k];

                g_y_0_xxyyyy_xzz[k] = -g_y_0_xyyyy_xzz[k] * ab_x + g_y_0_xyyyy_xxzz[k];

                g_y_0_xxyyyy_yyy[k] = -g_y_0_xyyyy_yyy[k] * ab_x + g_y_0_xyyyy_xyyy[k];

                g_y_0_xxyyyy_yyz[k] = -g_y_0_xyyyy_yyz[k] * ab_x + g_y_0_xyyyy_xyyz[k];

                g_y_0_xxyyyy_yzz[k] = -g_y_0_xyyyy_yzz[k] * ab_x + g_y_0_xyyyy_xyzz[k];

                g_y_0_xxyyyy_zzz[k] = -g_y_0_xyyyy_zzz[k] * ab_x + g_y_0_xyyyy_xzzz[k];
            }

            /// Set up 390-400 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_xxx = cbuffer.data(if_geom_10_off + 390 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxy = cbuffer.data(if_geom_10_off + 391 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxz = cbuffer.data(if_geom_10_off + 392 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyy = cbuffer.data(if_geom_10_off + 393 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyz = cbuffer.data(if_geom_10_off + 394 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xzz = cbuffer.data(if_geom_10_off + 395 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyy = cbuffer.data(if_geom_10_off + 396 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyz = cbuffer.data(if_geom_10_off + 397 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yzz = cbuffer.data(if_geom_10_off + 398 * ccomps * dcomps);

            auto g_y_0_xxyyyz_zzz = cbuffer.data(if_geom_10_off + 399 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyyz_xxx, g_y_0_xxyyyz_xxy, g_y_0_xxyyyz_xxz, g_y_0_xxyyyz_xyy, g_y_0_xxyyyz_xyz, g_y_0_xxyyyz_xzz, g_y_0_xxyyyz_yyy, g_y_0_xxyyyz_yyz, g_y_0_xxyyyz_yzz, g_y_0_xxyyyz_zzz, g_y_0_xyyyz_xxx, g_y_0_xyyyz_xxxx, g_y_0_xyyyz_xxxy, g_y_0_xyyyz_xxxz, g_y_0_xyyyz_xxy, g_y_0_xyyyz_xxyy, g_y_0_xyyyz_xxyz, g_y_0_xyyyz_xxz, g_y_0_xyyyz_xxzz, g_y_0_xyyyz_xyy, g_y_0_xyyyz_xyyy, g_y_0_xyyyz_xyyz, g_y_0_xyyyz_xyz, g_y_0_xyyyz_xyzz, g_y_0_xyyyz_xzz, g_y_0_xyyyz_xzzz, g_y_0_xyyyz_yyy, g_y_0_xyyyz_yyz, g_y_0_xyyyz_yzz, g_y_0_xyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_xxx[k] = -g_y_0_xyyyz_xxx[k] * ab_x + g_y_0_xyyyz_xxxx[k];

                g_y_0_xxyyyz_xxy[k] = -g_y_0_xyyyz_xxy[k] * ab_x + g_y_0_xyyyz_xxxy[k];

                g_y_0_xxyyyz_xxz[k] = -g_y_0_xyyyz_xxz[k] * ab_x + g_y_0_xyyyz_xxxz[k];

                g_y_0_xxyyyz_xyy[k] = -g_y_0_xyyyz_xyy[k] * ab_x + g_y_0_xyyyz_xxyy[k];

                g_y_0_xxyyyz_xyz[k] = -g_y_0_xyyyz_xyz[k] * ab_x + g_y_0_xyyyz_xxyz[k];

                g_y_0_xxyyyz_xzz[k] = -g_y_0_xyyyz_xzz[k] * ab_x + g_y_0_xyyyz_xxzz[k];

                g_y_0_xxyyyz_yyy[k] = -g_y_0_xyyyz_yyy[k] * ab_x + g_y_0_xyyyz_xyyy[k];

                g_y_0_xxyyyz_yyz[k] = -g_y_0_xyyyz_yyz[k] * ab_x + g_y_0_xyyyz_xyyz[k];

                g_y_0_xxyyyz_yzz[k] = -g_y_0_xyyyz_yzz[k] * ab_x + g_y_0_xyyyz_xyzz[k];

                g_y_0_xxyyyz_zzz[k] = -g_y_0_xyyyz_zzz[k] * ab_x + g_y_0_xyyyz_xzzz[k];
            }

            /// Set up 400-410 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_xxx = cbuffer.data(if_geom_10_off + 400 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxy = cbuffer.data(if_geom_10_off + 401 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxz = cbuffer.data(if_geom_10_off + 402 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyy = cbuffer.data(if_geom_10_off + 403 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyz = cbuffer.data(if_geom_10_off + 404 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xzz = cbuffer.data(if_geom_10_off + 405 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyy = cbuffer.data(if_geom_10_off + 406 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyz = cbuffer.data(if_geom_10_off + 407 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yzz = cbuffer.data(if_geom_10_off + 408 * ccomps * dcomps);

            auto g_y_0_xxyyzz_zzz = cbuffer.data(if_geom_10_off + 409 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyzz_xxx, g_y_0_xxyyzz_xxy, g_y_0_xxyyzz_xxz, g_y_0_xxyyzz_xyy, g_y_0_xxyyzz_xyz, g_y_0_xxyyzz_xzz, g_y_0_xxyyzz_yyy, g_y_0_xxyyzz_yyz, g_y_0_xxyyzz_yzz, g_y_0_xxyyzz_zzz, g_y_0_xyyzz_xxx, g_y_0_xyyzz_xxxx, g_y_0_xyyzz_xxxy, g_y_0_xyyzz_xxxz, g_y_0_xyyzz_xxy, g_y_0_xyyzz_xxyy, g_y_0_xyyzz_xxyz, g_y_0_xyyzz_xxz, g_y_0_xyyzz_xxzz, g_y_0_xyyzz_xyy, g_y_0_xyyzz_xyyy, g_y_0_xyyzz_xyyz, g_y_0_xyyzz_xyz, g_y_0_xyyzz_xyzz, g_y_0_xyyzz_xzz, g_y_0_xyyzz_xzzz, g_y_0_xyyzz_yyy, g_y_0_xyyzz_yyz, g_y_0_xyyzz_yzz, g_y_0_xyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_xxx[k] = -g_y_0_xyyzz_xxx[k] * ab_x + g_y_0_xyyzz_xxxx[k];

                g_y_0_xxyyzz_xxy[k] = -g_y_0_xyyzz_xxy[k] * ab_x + g_y_0_xyyzz_xxxy[k];

                g_y_0_xxyyzz_xxz[k] = -g_y_0_xyyzz_xxz[k] * ab_x + g_y_0_xyyzz_xxxz[k];

                g_y_0_xxyyzz_xyy[k] = -g_y_0_xyyzz_xyy[k] * ab_x + g_y_0_xyyzz_xxyy[k];

                g_y_0_xxyyzz_xyz[k] = -g_y_0_xyyzz_xyz[k] * ab_x + g_y_0_xyyzz_xxyz[k];

                g_y_0_xxyyzz_xzz[k] = -g_y_0_xyyzz_xzz[k] * ab_x + g_y_0_xyyzz_xxzz[k];

                g_y_0_xxyyzz_yyy[k] = -g_y_0_xyyzz_yyy[k] * ab_x + g_y_0_xyyzz_xyyy[k];

                g_y_0_xxyyzz_yyz[k] = -g_y_0_xyyzz_yyz[k] * ab_x + g_y_0_xyyzz_xyyz[k];

                g_y_0_xxyyzz_yzz[k] = -g_y_0_xyyzz_yzz[k] * ab_x + g_y_0_xyyzz_xyzz[k];

                g_y_0_xxyyzz_zzz[k] = -g_y_0_xyyzz_zzz[k] * ab_x + g_y_0_xyyzz_xzzz[k];
            }

            /// Set up 410-420 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_xxx = cbuffer.data(if_geom_10_off + 410 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxy = cbuffer.data(if_geom_10_off + 411 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxz = cbuffer.data(if_geom_10_off + 412 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyy = cbuffer.data(if_geom_10_off + 413 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyz = cbuffer.data(if_geom_10_off + 414 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xzz = cbuffer.data(if_geom_10_off + 415 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyy = cbuffer.data(if_geom_10_off + 416 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyz = cbuffer.data(if_geom_10_off + 417 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yzz = cbuffer.data(if_geom_10_off + 418 * ccomps * dcomps);

            auto g_y_0_xxyzzz_zzz = cbuffer.data(if_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyzzz_xxx, g_y_0_xxyzzz_xxy, g_y_0_xxyzzz_xxz, g_y_0_xxyzzz_xyy, g_y_0_xxyzzz_xyz, g_y_0_xxyzzz_xzz, g_y_0_xxyzzz_yyy, g_y_0_xxyzzz_yyz, g_y_0_xxyzzz_yzz, g_y_0_xxyzzz_zzz, g_y_0_xyzzz_xxx, g_y_0_xyzzz_xxxx, g_y_0_xyzzz_xxxy, g_y_0_xyzzz_xxxz, g_y_0_xyzzz_xxy, g_y_0_xyzzz_xxyy, g_y_0_xyzzz_xxyz, g_y_0_xyzzz_xxz, g_y_0_xyzzz_xxzz, g_y_0_xyzzz_xyy, g_y_0_xyzzz_xyyy, g_y_0_xyzzz_xyyz, g_y_0_xyzzz_xyz, g_y_0_xyzzz_xyzz, g_y_0_xyzzz_xzz, g_y_0_xyzzz_xzzz, g_y_0_xyzzz_yyy, g_y_0_xyzzz_yyz, g_y_0_xyzzz_yzz, g_y_0_xyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_xxx[k] = -g_y_0_xyzzz_xxx[k] * ab_x + g_y_0_xyzzz_xxxx[k];

                g_y_0_xxyzzz_xxy[k] = -g_y_0_xyzzz_xxy[k] * ab_x + g_y_0_xyzzz_xxxy[k];

                g_y_0_xxyzzz_xxz[k] = -g_y_0_xyzzz_xxz[k] * ab_x + g_y_0_xyzzz_xxxz[k];

                g_y_0_xxyzzz_xyy[k] = -g_y_0_xyzzz_xyy[k] * ab_x + g_y_0_xyzzz_xxyy[k];

                g_y_0_xxyzzz_xyz[k] = -g_y_0_xyzzz_xyz[k] * ab_x + g_y_0_xyzzz_xxyz[k];

                g_y_0_xxyzzz_xzz[k] = -g_y_0_xyzzz_xzz[k] * ab_x + g_y_0_xyzzz_xxzz[k];

                g_y_0_xxyzzz_yyy[k] = -g_y_0_xyzzz_yyy[k] * ab_x + g_y_0_xyzzz_xyyy[k];

                g_y_0_xxyzzz_yyz[k] = -g_y_0_xyzzz_yyz[k] * ab_x + g_y_0_xyzzz_xyyz[k];

                g_y_0_xxyzzz_yzz[k] = -g_y_0_xyzzz_yzz[k] * ab_x + g_y_0_xyzzz_xyzz[k];

                g_y_0_xxyzzz_zzz[k] = -g_y_0_xyzzz_zzz[k] * ab_x + g_y_0_xyzzz_xzzz[k];
            }

            /// Set up 420-430 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_xxx = cbuffer.data(if_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxy = cbuffer.data(if_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxz = cbuffer.data(if_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyy = cbuffer.data(if_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyz = cbuffer.data(if_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xzz = cbuffer.data(if_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyy = cbuffer.data(if_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyz = cbuffer.data(if_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yzz = cbuffer.data(if_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_xxzzzz_zzz = cbuffer.data(if_geom_10_off + 429 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzzzz_xxx, g_y_0_xxzzzz_xxy, g_y_0_xxzzzz_xxz, g_y_0_xxzzzz_xyy, g_y_0_xxzzzz_xyz, g_y_0_xxzzzz_xzz, g_y_0_xxzzzz_yyy, g_y_0_xxzzzz_yyz, g_y_0_xxzzzz_yzz, g_y_0_xxzzzz_zzz, g_y_0_xzzzz_xxx, g_y_0_xzzzz_xxxx, g_y_0_xzzzz_xxxy, g_y_0_xzzzz_xxxz, g_y_0_xzzzz_xxy, g_y_0_xzzzz_xxyy, g_y_0_xzzzz_xxyz, g_y_0_xzzzz_xxz, g_y_0_xzzzz_xxzz, g_y_0_xzzzz_xyy, g_y_0_xzzzz_xyyy, g_y_0_xzzzz_xyyz, g_y_0_xzzzz_xyz, g_y_0_xzzzz_xyzz, g_y_0_xzzzz_xzz, g_y_0_xzzzz_xzzz, g_y_0_xzzzz_yyy, g_y_0_xzzzz_yyz, g_y_0_xzzzz_yzz, g_y_0_xzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_xxx[k] = -g_y_0_xzzzz_xxx[k] * ab_x + g_y_0_xzzzz_xxxx[k];

                g_y_0_xxzzzz_xxy[k] = -g_y_0_xzzzz_xxy[k] * ab_x + g_y_0_xzzzz_xxxy[k];

                g_y_0_xxzzzz_xxz[k] = -g_y_0_xzzzz_xxz[k] * ab_x + g_y_0_xzzzz_xxxz[k];

                g_y_0_xxzzzz_xyy[k] = -g_y_0_xzzzz_xyy[k] * ab_x + g_y_0_xzzzz_xxyy[k];

                g_y_0_xxzzzz_xyz[k] = -g_y_0_xzzzz_xyz[k] * ab_x + g_y_0_xzzzz_xxyz[k];

                g_y_0_xxzzzz_xzz[k] = -g_y_0_xzzzz_xzz[k] * ab_x + g_y_0_xzzzz_xxzz[k];

                g_y_0_xxzzzz_yyy[k] = -g_y_0_xzzzz_yyy[k] * ab_x + g_y_0_xzzzz_xyyy[k];

                g_y_0_xxzzzz_yyz[k] = -g_y_0_xzzzz_yyz[k] * ab_x + g_y_0_xzzzz_xyyz[k];

                g_y_0_xxzzzz_yzz[k] = -g_y_0_xzzzz_yzz[k] * ab_x + g_y_0_xzzzz_xyzz[k];

                g_y_0_xxzzzz_zzz[k] = -g_y_0_xzzzz_zzz[k] * ab_x + g_y_0_xzzzz_xzzz[k];
            }

            /// Set up 430-440 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_xxx = cbuffer.data(if_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxy = cbuffer.data(if_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxz = cbuffer.data(if_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyy = cbuffer.data(if_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyz = cbuffer.data(if_geom_10_off + 434 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xzz = cbuffer.data(if_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyy = cbuffer.data(if_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyz = cbuffer.data(if_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yzz = cbuffer.data(if_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_xyyyyy_zzz = cbuffer.data(if_geom_10_off + 439 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyyy_xxx, g_y_0_xyyyyy_xxy, g_y_0_xyyyyy_xxz, g_y_0_xyyyyy_xyy, g_y_0_xyyyyy_xyz, g_y_0_xyyyyy_xzz, g_y_0_xyyyyy_yyy, g_y_0_xyyyyy_yyz, g_y_0_xyyyyy_yzz, g_y_0_xyyyyy_zzz, g_y_0_yyyyy_xxx, g_y_0_yyyyy_xxxx, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxy, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xyy, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_yyy, g_y_0_yyyyy_yyz, g_y_0_yyyyy_yzz, g_y_0_yyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_xxx[k] = -g_y_0_yyyyy_xxx[k] * ab_x + g_y_0_yyyyy_xxxx[k];

                g_y_0_xyyyyy_xxy[k] = -g_y_0_yyyyy_xxy[k] * ab_x + g_y_0_yyyyy_xxxy[k];

                g_y_0_xyyyyy_xxz[k] = -g_y_0_yyyyy_xxz[k] * ab_x + g_y_0_yyyyy_xxxz[k];

                g_y_0_xyyyyy_xyy[k] = -g_y_0_yyyyy_xyy[k] * ab_x + g_y_0_yyyyy_xxyy[k];

                g_y_0_xyyyyy_xyz[k] = -g_y_0_yyyyy_xyz[k] * ab_x + g_y_0_yyyyy_xxyz[k];

                g_y_0_xyyyyy_xzz[k] = -g_y_0_yyyyy_xzz[k] * ab_x + g_y_0_yyyyy_xxzz[k];

                g_y_0_xyyyyy_yyy[k] = -g_y_0_yyyyy_yyy[k] * ab_x + g_y_0_yyyyy_xyyy[k];

                g_y_0_xyyyyy_yyz[k] = -g_y_0_yyyyy_yyz[k] * ab_x + g_y_0_yyyyy_xyyz[k];

                g_y_0_xyyyyy_yzz[k] = -g_y_0_yyyyy_yzz[k] * ab_x + g_y_0_yyyyy_xyzz[k];

                g_y_0_xyyyyy_zzz[k] = -g_y_0_yyyyy_zzz[k] * ab_x + g_y_0_yyyyy_xzzz[k];
            }

            /// Set up 440-450 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_xxx = cbuffer.data(if_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxy = cbuffer.data(if_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxz = cbuffer.data(if_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyy = cbuffer.data(if_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyz = cbuffer.data(if_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xzz = cbuffer.data(if_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyy = cbuffer.data(if_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyz = cbuffer.data(if_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yzz = cbuffer.data(if_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_xyyyyz_zzz = cbuffer.data(if_geom_10_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyyz_xxx, g_y_0_xyyyyz_xxy, g_y_0_xyyyyz_xxz, g_y_0_xyyyyz_xyy, g_y_0_xyyyyz_xyz, g_y_0_xyyyyz_xzz, g_y_0_xyyyyz_yyy, g_y_0_xyyyyz_yyz, g_y_0_xyyyyz_yzz, g_y_0_xyyyyz_zzz, g_y_0_yyyyz_xxx, g_y_0_yyyyz_xxxx, g_y_0_yyyyz_xxxy, g_y_0_yyyyz_xxxz, g_y_0_yyyyz_xxy, g_y_0_yyyyz_xxyy, g_y_0_yyyyz_xxyz, g_y_0_yyyyz_xxz, g_y_0_yyyyz_xxzz, g_y_0_yyyyz_xyy, g_y_0_yyyyz_xyyy, g_y_0_yyyyz_xyyz, g_y_0_yyyyz_xyz, g_y_0_yyyyz_xyzz, g_y_0_yyyyz_xzz, g_y_0_yyyyz_xzzz, g_y_0_yyyyz_yyy, g_y_0_yyyyz_yyz, g_y_0_yyyyz_yzz, g_y_0_yyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_xxx[k] = -g_y_0_yyyyz_xxx[k] * ab_x + g_y_0_yyyyz_xxxx[k];

                g_y_0_xyyyyz_xxy[k] = -g_y_0_yyyyz_xxy[k] * ab_x + g_y_0_yyyyz_xxxy[k];

                g_y_0_xyyyyz_xxz[k] = -g_y_0_yyyyz_xxz[k] * ab_x + g_y_0_yyyyz_xxxz[k];

                g_y_0_xyyyyz_xyy[k] = -g_y_0_yyyyz_xyy[k] * ab_x + g_y_0_yyyyz_xxyy[k];

                g_y_0_xyyyyz_xyz[k] = -g_y_0_yyyyz_xyz[k] * ab_x + g_y_0_yyyyz_xxyz[k];

                g_y_0_xyyyyz_xzz[k] = -g_y_0_yyyyz_xzz[k] * ab_x + g_y_0_yyyyz_xxzz[k];

                g_y_0_xyyyyz_yyy[k] = -g_y_0_yyyyz_yyy[k] * ab_x + g_y_0_yyyyz_xyyy[k];

                g_y_0_xyyyyz_yyz[k] = -g_y_0_yyyyz_yyz[k] * ab_x + g_y_0_yyyyz_xyyz[k];

                g_y_0_xyyyyz_yzz[k] = -g_y_0_yyyyz_yzz[k] * ab_x + g_y_0_yyyyz_xyzz[k];

                g_y_0_xyyyyz_zzz[k] = -g_y_0_yyyyz_zzz[k] * ab_x + g_y_0_yyyyz_xzzz[k];
            }

            /// Set up 450-460 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_xxx = cbuffer.data(if_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxy = cbuffer.data(if_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxz = cbuffer.data(if_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyy = cbuffer.data(if_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyz = cbuffer.data(if_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xzz = cbuffer.data(if_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyy = cbuffer.data(if_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyz = cbuffer.data(if_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yzz = cbuffer.data(if_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_xyyyzz_zzz = cbuffer.data(if_geom_10_off + 459 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyzz_xxx, g_y_0_xyyyzz_xxy, g_y_0_xyyyzz_xxz, g_y_0_xyyyzz_xyy, g_y_0_xyyyzz_xyz, g_y_0_xyyyzz_xzz, g_y_0_xyyyzz_yyy, g_y_0_xyyyzz_yyz, g_y_0_xyyyzz_yzz, g_y_0_xyyyzz_zzz, g_y_0_yyyzz_xxx, g_y_0_yyyzz_xxxx, g_y_0_yyyzz_xxxy, g_y_0_yyyzz_xxxz, g_y_0_yyyzz_xxy, g_y_0_yyyzz_xxyy, g_y_0_yyyzz_xxyz, g_y_0_yyyzz_xxz, g_y_0_yyyzz_xxzz, g_y_0_yyyzz_xyy, g_y_0_yyyzz_xyyy, g_y_0_yyyzz_xyyz, g_y_0_yyyzz_xyz, g_y_0_yyyzz_xyzz, g_y_0_yyyzz_xzz, g_y_0_yyyzz_xzzz, g_y_0_yyyzz_yyy, g_y_0_yyyzz_yyz, g_y_0_yyyzz_yzz, g_y_0_yyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_xxx[k] = -g_y_0_yyyzz_xxx[k] * ab_x + g_y_0_yyyzz_xxxx[k];

                g_y_0_xyyyzz_xxy[k] = -g_y_0_yyyzz_xxy[k] * ab_x + g_y_0_yyyzz_xxxy[k];

                g_y_0_xyyyzz_xxz[k] = -g_y_0_yyyzz_xxz[k] * ab_x + g_y_0_yyyzz_xxxz[k];

                g_y_0_xyyyzz_xyy[k] = -g_y_0_yyyzz_xyy[k] * ab_x + g_y_0_yyyzz_xxyy[k];

                g_y_0_xyyyzz_xyz[k] = -g_y_0_yyyzz_xyz[k] * ab_x + g_y_0_yyyzz_xxyz[k];

                g_y_0_xyyyzz_xzz[k] = -g_y_0_yyyzz_xzz[k] * ab_x + g_y_0_yyyzz_xxzz[k];

                g_y_0_xyyyzz_yyy[k] = -g_y_0_yyyzz_yyy[k] * ab_x + g_y_0_yyyzz_xyyy[k];

                g_y_0_xyyyzz_yyz[k] = -g_y_0_yyyzz_yyz[k] * ab_x + g_y_0_yyyzz_xyyz[k];

                g_y_0_xyyyzz_yzz[k] = -g_y_0_yyyzz_yzz[k] * ab_x + g_y_0_yyyzz_xyzz[k];

                g_y_0_xyyyzz_zzz[k] = -g_y_0_yyyzz_zzz[k] * ab_x + g_y_0_yyyzz_xzzz[k];
            }

            /// Set up 460-470 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_xxx = cbuffer.data(if_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxy = cbuffer.data(if_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxz = cbuffer.data(if_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyy = cbuffer.data(if_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyz = cbuffer.data(if_geom_10_off + 464 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xzz = cbuffer.data(if_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyy = cbuffer.data(if_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyz = cbuffer.data(if_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yzz = cbuffer.data(if_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_xyyzzz_zzz = cbuffer.data(if_geom_10_off + 469 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyzzz_xxx, g_y_0_xyyzzz_xxy, g_y_0_xyyzzz_xxz, g_y_0_xyyzzz_xyy, g_y_0_xyyzzz_xyz, g_y_0_xyyzzz_xzz, g_y_0_xyyzzz_yyy, g_y_0_xyyzzz_yyz, g_y_0_xyyzzz_yzz, g_y_0_xyyzzz_zzz, g_y_0_yyzzz_xxx, g_y_0_yyzzz_xxxx, g_y_0_yyzzz_xxxy, g_y_0_yyzzz_xxxz, g_y_0_yyzzz_xxy, g_y_0_yyzzz_xxyy, g_y_0_yyzzz_xxyz, g_y_0_yyzzz_xxz, g_y_0_yyzzz_xxzz, g_y_0_yyzzz_xyy, g_y_0_yyzzz_xyyy, g_y_0_yyzzz_xyyz, g_y_0_yyzzz_xyz, g_y_0_yyzzz_xyzz, g_y_0_yyzzz_xzz, g_y_0_yyzzz_xzzz, g_y_0_yyzzz_yyy, g_y_0_yyzzz_yyz, g_y_0_yyzzz_yzz, g_y_0_yyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_xxx[k] = -g_y_0_yyzzz_xxx[k] * ab_x + g_y_0_yyzzz_xxxx[k];

                g_y_0_xyyzzz_xxy[k] = -g_y_0_yyzzz_xxy[k] * ab_x + g_y_0_yyzzz_xxxy[k];

                g_y_0_xyyzzz_xxz[k] = -g_y_0_yyzzz_xxz[k] * ab_x + g_y_0_yyzzz_xxxz[k];

                g_y_0_xyyzzz_xyy[k] = -g_y_0_yyzzz_xyy[k] * ab_x + g_y_0_yyzzz_xxyy[k];

                g_y_0_xyyzzz_xyz[k] = -g_y_0_yyzzz_xyz[k] * ab_x + g_y_0_yyzzz_xxyz[k];

                g_y_0_xyyzzz_xzz[k] = -g_y_0_yyzzz_xzz[k] * ab_x + g_y_0_yyzzz_xxzz[k];

                g_y_0_xyyzzz_yyy[k] = -g_y_0_yyzzz_yyy[k] * ab_x + g_y_0_yyzzz_xyyy[k];

                g_y_0_xyyzzz_yyz[k] = -g_y_0_yyzzz_yyz[k] * ab_x + g_y_0_yyzzz_xyyz[k];

                g_y_0_xyyzzz_yzz[k] = -g_y_0_yyzzz_yzz[k] * ab_x + g_y_0_yyzzz_xyzz[k];

                g_y_0_xyyzzz_zzz[k] = -g_y_0_yyzzz_zzz[k] * ab_x + g_y_0_yyzzz_xzzz[k];
            }

            /// Set up 470-480 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_xxx = cbuffer.data(if_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxy = cbuffer.data(if_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxz = cbuffer.data(if_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyy = cbuffer.data(if_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyz = cbuffer.data(if_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xzz = cbuffer.data(if_geom_10_off + 475 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyy = cbuffer.data(if_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyz = cbuffer.data(if_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yzz = cbuffer.data(if_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_xyzzzz_zzz = cbuffer.data(if_geom_10_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzzzz_xxx, g_y_0_xyzzzz_xxy, g_y_0_xyzzzz_xxz, g_y_0_xyzzzz_xyy, g_y_0_xyzzzz_xyz, g_y_0_xyzzzz_xzz, g_y_0_xyzzzz_yyy, g_y_0_xyzzzz_yyz, g_y_0_xyzzzz_yzz, g_y_0_xyzzzz_zzz, g_y_0_yzzzz_xxx, g_y_0_yzzzz_xxxx, g_y_0_yzzzz_xxxy, g_y_0_yzzzz_xxxz, g_y_0_yzzzz_xxy, g_y_0_yzzzz_xxyy, g_y_0_yzzzz_xxyz, g_y_0_yzzzz_xxz, g_y_0_yzzzz_xxzz, g_y_0_yzzzz_xyy, g_y_0_yzzzz_xyyy, g_y_0_yzzzz_xyyz, g_y_0_yzzzz_xyz, g_y_0_yzzzz_xyzz, g_y_0_yzzzz_xzz, g_y_0_yzzzz_xzzz, g_y_0_yzzzz_yyy, g_y_0_yzzzz_yyz, g_y_0_yzzzz_yzz, g_y_0_yzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_xxx[k] = -g_y_0_yzzzz_xxx[k] * ab_x + g_y_0_yzzzz_xxxx[k];

                g_y_0_xyzzzz_xxy[k] = -g_y_0_yzzzz_xxy[k] * ab_x + g_y_0_yzzzz_xxxy[k];

                g_y_0_xyzzzz_xxz[k] = -g_y_0_yzzzz_xxz[k] * ab_x + g_y_0_yzzzz_xxxz[k];

                g_y_0_xyzzzz_xyy[k] = -g_y_0_yzzzz_xyy[k] * ab_x + g_y_0_yzzzz_xxyy[k];

                g_y_0_xyzzzz_xyz[k] = -g_y_0_yzzzz_xyz[k] * ab_x + g_y_0_yzzzz_xxyz[k];

                g_y_0_xyzzzz_xzz[k] = -g_y_0_yzzzz_xzz[k] * ab_x + g_y_0_yzzzz_xxzz[k];

                g_y_0_xyzzzz_yyy[k] = -g_y_0_yzzzz_yyy[k] * ab_x + g_y_0_yzzzz_xyyy[k];

                g_y_0_xyzzzz_yyz[k] = -g_y_0_yzzzz_yyz[k] * ab_x + g_y_0_yzzzz_xyyz[k];

                g_y_0_xyzzzz_yzz[k] = -g_y_0_yzzzz_yzz[k] * ab_x + g_y_0_yzzzz_xyzz[k];

                g_y_0_xyzzzz_zzz[k] = -g_y_0_yzzzz_zzz[k] * ab_x + g_y_0_yzzzz_xzzz[k];
            }

            /// Set up 480-490 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_xxx = cbuffer.data(if_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxy = cbuffer.data(if_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxz = cbuffer.data(if_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyy = cbuffer.data(if_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyz = cbuffer.data(if_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xzz = cbuffer.data(if_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyy = cbuffer.data(if_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyz = cbuffer.data(if_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yzz = cbuffer.data(if_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_xzzzzz_zzz = cbuffer.data(if_geom_10_off + 489 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzzzz_xxx, g_y_0_xzzzzz_xxy, g_y_0_xzzzzz_xxz, g_y_0_xzzzzz_xyy, g_y_0_xzzzzz_xyz, g_y_0_xzzzzz_xzz, g_y_0_xzzzzz_yyy, g_y_0_xzzzzz_yyz, g_y_0_xzzzzz_yzz, g_y_0_xzzzzz_zzz, g_y_0_zzzzz_xxx, g_y_0_zzzzz_xxxx, g_y_0_zzzzz_xxxy, g_y_0_zzzzz_xxxz, g_y_0_zzzzz_xxy, g_y_0_zzzzz_xxyy, g_y_0_zzzzz_xxyz, g_y_0_zzzzz_xxz, g_y_0_zzzzz_xxzz, g_y_0_zzzzz_xyy, g_y_0_zzzzz_xyyy, g_y_0_zzzzz_xyyz, g_y_0_zzzzz_xyz, g_y_0_zzzzz_xyzz, g_y_0_zzzzz_xzz, g_y_0_zzzzz_xzzz, g_y_0_zzzzz_yyy, g_y_0_zzzzz_yyz, g_y_0_zzzzz_yzz, g_y_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_xxx[k] = -g_y_0_zzzzz_xxx[k] * ab_x + g_y_0_zzzzz_xxxx[k];

                g_y_0_xzzzzz_xxy[k] = -g_y_0_zzzzz_xxy[k] * ab_x + g_y_0_zzzzz_xxxy[k];

                g_y_0_xzzzzz_xxz[k] = -g_y_0_zzzzz_xxz[k] * ab_x + g_y_0_zzzzz_xxxz[k];

                g_y_0_xzzzzz_xyy[k] = -g_y_0_zzzzz_xyy[k] * ab_x + g_y_0_zzzzz_xxyy[k];

                g_y_0_xzzzzz_xyz[k] = -g_y_0_zzzzz_xyz[k] * ab_x + g_y_0_zzzzz_xxyz[k];

                g_y_0_xzzzzz_xzz[k] = -g_y_0_zzzzz_xzz[k] * ab_x + g_y_0_zzzzz_xxzz[k];

                g_y_0_xzzzzz_yyy[k] = -g_y_0_zzzzz_yyy[k] * ab_x + g_y_0_zzzzz_xyyy[k];

                g_y_0_xzzzzz_yyz[k] = -g_y_0_zzzzz_yyz[k] * ab_x + g_y_0_zzzzz_xyyz[k];

                g_y_0_xzzzzz_yzz[k] = -g_y_0_zzzzz_yzz[k] * ab_x + g_y_0_zzzzz_xyzz[k];

                g_y_0_xzzzzz_zzz[k] = -g_y_0_zzzzz_zzz[k] * ab_x + g_y_0_zzzzz_xzzz[k];
            }

            /// Set up 490-500 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_xxx = cbuffer.data(if_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxy = cbuffer.data(if_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxz = cbuffer.data(if_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyy = cbuffer.data(if_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyz = cbuffer.data(if_geom_10_off + 494 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xzz = cbuffer.data(if_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyy = cbuffer.data(if_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyz = cbuffer.data(if_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yzz = cbuffer.data(if_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_yyyyyy_zzz = cbuffer.data(if_geom_10_off + 499 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_xxx, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxy, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxz, g_y_0_yyyyy_xyy, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xzz, g_y_0_yyyyy_yyy, g_y_0_yyyyy_yyyy, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_zzz, g_y_0_yyyyyy_xxx, g_y_0_yyyyyy_xxy, g_y_0_yyyyyy_xxz, g_y_0_yyyyyy_xyy, g_y_0_yyyyyy_xyz, g_y_0_yyyyyy_xzz, g_y_0_yyyyyy_yyy, g_y_0_yyyyyy_yyz, g_y_0_yyyyyy_yzz, g_y_0_yyyyyy_zzz, g_yyyyy_xxx, g_yyyyy_xxy, g_yyyyy_xxz, g_yyyyy_xyy, g_yyyyy_xyz, g_yyyyy_xzz, g_yyyyy_yyy, g_yyyyy_yyz, g_yyyyy_yzz, g_yyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_xxx[k] = -g_yyyyy_xxx[k] - g_y_0_yyyyy_xxx[k] * ab_y + g_y_0_yyyyy_xxxy[k];

                g_y_0_yyyyyy_xxy[k] = -g_yyyyy_xxy[k] - g_y_0_yyyyy_xxy[k] * ab_y + g_y_0_yyyyy_xxyy[k];

                g_y_0_yyyyyy_xxz[k] = -g_yyyyy_xxz[k] - g_y_0_yyyyy_xxz[k] * ab_y + g_y_0_yyyyy_xxyz[k];

                g_y_0_yyyyyy_xyy[k] = -g_yyyyy_xyy[k] - g_y_0_yyyyy_xyy[k] * ab_y + g_y_0_yyyyy_xyyy[k];

                g_y_0_yyyyyy_xyz[k] = -g_yyyyy_xyz[k] - g_y_0_yyyyy_xyz[k] * ab_y + g_y_0_yyyyy_xyyz[k];

                g_y_0_yyyyyy_xzz[k] = -g_yyyyy_xzz[k] - g_y_0_yyyyy_xzz[k] * ab_y + g_y_0_yyyyy_xyzz[k];

                g_y_0_yyyyyy_yyy[k] = -g_yyyyy_yyy[k] - g_y_0_yyyyy_yyy[k] * ab_y + g_y_0_yyyyy_yyyy[k];

                g_y_0_yyyyyy_yyz[k] = -g_yyyyy_yyz[k] - g_y_0_yyyyy_yyz[k] * ab_y + g_y_0_yyyyy_yyyz[k];

                g_y_0_yyyyyy_yzz[k] = -g_yyyyy_yzz[k] - g_y_0_yyyyy_yzz[k] * ab_y + g_y_0_yyyyy_yyzz[k];

                g_y_0_yyyyyy_zzz[k] = -g_yyyyy_zzz[k] - g_y_0_yyyyy_zzz[k] * ab_y + g_y_0_yyyyy_yzzz[k];
            }

            /// Set up 500-510 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_xxx = cbuffer.data(if_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxy = cbuffer.data(if_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxz = cbuffer.data(if_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyy = cbuffer.data(if_geom_10_off + 503 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyz = cbuffer.data(if_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xzz = cbuffer.data(if_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyy = cbuffer.data(if_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyz = cbuffer.data(if_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yzz = cbuffer.data(if_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_yyyyyz_zzz = cbuffer.data(if_geom_10_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_xxx, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxy, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xyy, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_yyy, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_zzz, g_y_0_yyyyy_zzzz, g_y_0_yyyyyz_xxx, g_y_0_yyyyyz_xxy, g_y_0_yyyyyz_xxz, g_y_0_yyyyyz_xyy, g_y_0_yyyyyz_xyz, g_y_0_yyyyyz_xzz, g_y_0_yyyyyz_yyy, g_y_0_yyyyyz_yyz, g_y_0_yyyyyz_yzz, g_y_0_yyyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_xxx[k] = -g_y_0_yyyyy_xxx[k] * ab_z + g_y_0_yyyyy_xxxz[k];

                g_y_0_yyyyyz_xxy[k] = -g_y_0_yyyyy_xxy[k] * ab_z + g_y_0_yyyyy_xxyz[k];

                g_y_0_yyyyyz_xxz[k] = -g_y_0_yyyyy_xxz[k] * ab_z + g_y_0_yyyyy_xxzz[k];

                g_y_0_yyyyyz_xyy[k] = -g_y_0_yyyyy_xyy[k] * ab_z + g_y_0_yyyyy_xyyz[k];

                g_y_0_yyyyyz_xyz[k] = -g_y_0_yyyyy_xyz[k] * ab_z + g_y_0_yyyyy_xyzz[k];

                g_y_0_yyyyyz_xzz[k] = -g_y_0_yyyyy_xzz[k] * ab_z + g_y_0_yyyyy_xzzz[k];

                g_y_0_yyyyyz_yyy[k] = -g_y_0_yyyyy_yyy[k] * ab_z + g_y_0_yyyyy_yyyz[k];

                g_y_0_yyyyyz_yyz[k] = -g_y_0_yyyyy_yyz[k] * ab_z + g_y_0_yyyyy_yyzz[k];

                g_y_0_yyyyyz_yzz[k] = -g_y_0_yyyyy_yzz[k] * ab_z + g_y_0_yyyyy_yzzz[k];

                g_y_0_yyyyyz_zzz[k] = -g_y_0_yyyyy_zzz[k] * ab_z + g_y_0_yyyyy_zzzz[k];
            }

            /// Set up 510-520 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_xxx = cbuffer.data(if_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxy = cbuffer.data(if_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxz = cbuffer.data(if_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyy = cbuffer.data(if_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyz = cbuffer.data(if_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xzz = cbuffer.data(if_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyy = cbuffer.data(if_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyz = cbuffer.data(if_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yzz = cbuffer.data(if_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_yyyyzz_zzz = cbuffer.data(if_geom_10_off + 519 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyz_xxx, g_y_0_yyyyz_xxxz, g_y_0_yyyyz_xxy, g_y_0_yyyyz_xxyz, g_y_0_yyyyz_xxz, g_y_0_yyyyz_xxzz, g_y_0_yyyyz_xyy, g_y_0_yyyyz_xyyz, g_y_0_yyyyz_xyz, g_y_0_yyyyz_xyzz, g_y_0_yyyyz_xzz, g_y_0_yyyyz_xzzz, g_y_0_yyyyz_yyy, g_y_0_yyyyz_yyyz, g_y_0_yyyyz_yyz, g_y_0_yyyyz_yyzz, g_y_0_yyyyz_yzz, g_y_0_yyyyz_yzzz, g_y_0_yyyyz_zzz, g_y_0_yyyyz_zzzz, g_y_0_yyyyzz_xxx, g_y_0_yyyyzz_xxy, g_y_0_yyyyzz_xxz, g_y_0_yyyyzz_xyy, g_y_0_yyyyzz_xyz, g_y_0_yyyyzz_xzz, g_y_0_yyyyzz_yyy, g_y_0_yyyyzz_yyz, g_y_0_yyyyzz_yzz, g_y_0_yyyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_xxx[k] = -g_y_0_yyyyz_xxx[k] * ab_z + g_y_0_yyyyz_xxxz[k];

                g_y_0_yyyyzz_xxy[k] = -g_y_0_yyyyz_xxy[k] * ab_z + g_y_0_yyyyz_xxyz[k];

                g_y_0_yyyyzz_xxz[k] = -g_y_0_yyyyz_xxz[k] * ab_z + g_y_0_yyyyz_xxzz[k];

                g_y_0_yyyyzz_xyy[k] = -g_y_0_yyyyz_xyy[k] * ab_z + g_y_0_yyyyz_xyyz[k];

                g_y_0_yyyyzz_xyz[k] = -g_y_0_yyyyz_xyz[k] * ab_z + g_y_0_yyyyz_xyzz[k];

                g_y_0_yyyyzz_xzz[k] = -g_y_0_yyyyz_xzz[k] * ab_z + g_y_0_yyyyz_xzzz[k];

                g_y_0_yyyyzz_yyy[k] = -g_y_0_yyyyz_yyy[k] * ab_z + g_y_0_yyyyz_yyyz[k];

                g_y_0_yyyyzz_yyz[k] = -g_y_0_yyyyz_yyz[k] * ab_z + g_y_0_yyyyz_yyzz[k];

                g_y_0_yyyyzz_yzz[k] = -g_y_0_yyyyz_yzz[k] * ab_z + g_y_0_yyyyz_yzzz[k];

                g_y_0_yyyyzz_zzz[k] = -g_y_0_yyyyz_zzz[k] * ab_z + g_y_0_yyyyz_zzzz[k];
            }

            /// Set up 520-530 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_xxx = cbuffer.data(if_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxy = cbuffer.data(if_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxz = cbuffer.data(if_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyy = cbuffer.data(if_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyz = cbuffer.data(if_geom_10_off + 524 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xzz = cbuffer.data(if_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyy = cbuffer.data(if_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyz = cbuffer.data(if_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yzz = cbuffer.data(if_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_yyyzzz_zzz = cbuffer.data(if_geom_10_off + 529 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyzz_xxx, g_y_0_yyyzz_xxxz, g_y_0_yyyzz_xxy, g_y_0_yyyzz_xxyz, g_y_0_yyyzz_xxz, g_y_0_yyyzz_xxzz, g_y_0_yyyzz_xyy, g_y_0_yyyzz_xyyz, g_y_0_yyyzz_xyz, g_y_0_yyyzz_xyzz, g_y_0_yyyzz_xzz, g_y_0_yyyzz_xzzz, g_y_0_yyyzz_yyy, g_y_0_yyyzz_yyyz, g_y_0_yyyzz_yyz, g_y_0_yyyzz_yyzz, g_y_0_yyyzz_yzz, g_y_0_yyyzz_yzzz, g_y_0_yyyzz_zzz, g_y_0_yyyzz_zzzz, g_y_0_yyyzzz_xxx, g_y_0_yyyzzz_xxy, g_y_0_yyyzzz_xxz, g_y_0_yyyzzz_xyy, g_y_0_yyyzzz_xyz, g_y_0_yyyzzz_xzz, g_y_0_yyyzzz_yyy, g_y_0_yyyzzz_yyz, g_y_0_yyyzzz_yzz, g_y_0_yyyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_xxx[k] = -g_y_0_yyyzz_xxx[k] * ab_z + g_y_0_yyyzz_xxxz[k];

                g_y_0_yyyzzz_xxy[k] = -g_y_0_yyyzz_xxy[k] * ab_z + g_y_0_yyyzz_xxyz[k];

                g_y_0_yyyzzz_xxz[k] = -g_y_0_yyyzz_xxz[k] * ab_z + g_y_0_yyyzz_xxzz[k];

                g_y_0_yyyzzz_xyy[k] = -g_y_0_yyyzz_xyy[k] * ab_z + g_y_0_yyyzz_xyyz[k];

                g_y_0_yyyzzz_xyz[k] = -g_y_0_yyyzz_xyz[k] * ab_z + g_y_0_yyyzz_xyzz[k];

                g_y_0_yyyzzz_xzz[k] = -g_y_0_yyyzz_xzz[k] * ab_z + g_y_0_yyyzz_xzzz[k];

                g_y_0_yyyzzz_yyy[k] = -g_y_0_yyyzz_yyy[k] * ab_z + g_y_0_yyyzz_yyyz[k];

                g_y_0_yyyzzz_yyz[k] = -g_y_0_yyyzz_yyz[k] * ab_z + g_y_0_yyyzz_yyzz[k];

                g_y_0_yyyzzz_yzz[k] = -g_y_0_yyyzz_yzz[k] * ab_z + g_y_0_yyyzz_yzzz[k];

                g_y_0_yyyzzz_zzz[k] = -g_y_0_yyyzz_zzz[k] * ab_z + g_y_0_yyyzz_zzzz[k];
            }

            /// Set up 530-540 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_xxx = cbuffer.data(if_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxy = cbuffer.data(if_geom_10_off + 531 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxz = cbuffer.data(if_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyy = cbuffer.data(if_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyz = cbuffer.data(if_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xzz = cbuffer.data(if_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyy = cbuffer.data(if_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyz = cbuffer.data(if_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yzz = cbuffer.data(if_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_yyzzzz_zzz = cbuffer.data(if_geom_10_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyzzz_xxx, g_y_0_yyzzz_xxxz, g_y_0_yyzzz_xxy, g_y_0_yyzzz_xxyz, g_y_0_yyzzz_xxz, g_y_0_yyzzz_xxzz, g_y_0_yyzzz_xyy, g_y_0_yyzzz_xyyz, g_y_0_yyzzz_xyz, g_y_0_yyzzz_xyzz, g_y_0_yyzzz_xzz, g_y_0_yyzzz_xzzz, g_y_0_yyzzz_yyy, g_y_0_yyzzz_yyyz, g_y_0_yyzzz_yyz, g_y_0_yyzzz_yyzz, g_y_0_yyzzz_yzz, g_y_0_yyzzz_yzzz, g_y_0_yyzzz_zzz, g_y_0_yyzzz_zzzz, g_y_0_yyzzzz_xxx, g_y_0_yyzzzz_xxy, g_y_0_yyzzzz_xxz, g_y_0_yyzzzz_xyy, g_y_0_yyzzzz_xyz, g_y_0_yyzzzz_xzz, g_y_0_yyzzzz_yyy, g_y_0_yyzzzz_yyz, g_y_0_yyzzzz_yzz, g_y_0_yyzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_xxx[k] = -g_y_0_yyzzz_xxx[k] * ab_z + g_y_0_yyzzz_xxxz[k];

                g_y_0_yyzzzz_xxy[k] = -g_y_0_yyzzz_xxy[k] * ab_z + g_y_0_yyzzz_xxyz[k];

                g_y_0_yyzzzz_xxz[k] = -g_y_0_yyzzz_xxz[k] * ab_z + g_y_0_yyzzz_xxzz[k];

                g_y_0_yyzzzz_xyy[k] = -g_y_0_yyzzz_xyy[k] * ab_z + g_y_0_yyzzz_xyyz[k];

                g_y_0_yyzzzz_xyz[k] = -g_y_0_yyzzz_xyz[k] * ab_z + g_y_0_yyzzz_xyzz[k];

                g_y_0_yyzzzz_xzz[k] = -g_y_0_yyzzz_xzz[k] * ab_z + g_y_0_yyzzz_xzzz[k];

                g_y_0_yyzzzz_yyy[k] = -g_y_0_yyzzz_yyy[k] * ab_z + g_y_0_yyzzz_yyyz[k];

                g_y_0_yyzzzz_yyz[k] = -g_y_0_yyzzz_yyz[k] * ab_z + g_y_0_yyzzz_yyzz[k];

                g_y_0_yyzzzz_yzz[k] = -g_y_0_yyzzz_yzz[k] * ab_z + g_y_0_yyzzz_yzzz[k];

                g_y_0_yyzzzz_zzz[k] = -g_y_0_yyzzz_zzz[k] * ab_z + g_y_0_yyzzz_zzzz[k];
            }

            /// Set up 540-550 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_xxx = cbuffer.data(if_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxy = cbuffer.data(if_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxz = cbuffer.data(if_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyy = cbuffer.data(if_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyz = cbuffer.data(if_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xzz = cbuffer.data(if_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyy = cbuffer.data(if_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyz = cbuffer.data(if_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yzz = cbuffer.data(if_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_yzzzzz_zzz = cbuffer.data(if_geom_10_off + 549 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzzzz_xxx, g_y_0_yzzzz_xxxz, g_y_0_yzzzz_xxy, g_y_0_yzzzz_xxyz, g_y_0_yzzzz_xxz, g_y_0_yzzzz_xxzz, g_y_0_yzzzz_xyy, g_y_0_yzzzz_xyyz, g_y_0_yzzzz_xyz, g_y_0_yzzzz_xyzz, g_y_0_yzzzz_xzz, g_y_0_yzzzz_xzzz, g_y_0_yzzzz_yyy, g_y_0_yzzzz_yyyz, g_y_0_yzzzz_yyz, g_y_0_yzzzz_yyzz, g_y_0_yzzzz_yzz, g_y_0_yzzzz_yzzz, g_y_0_yzzzz_zzz, g_y_0_yzzzz_zzzz, g_y_0_yzzzzz_xxx, g_y_0_yzzzzz_xxy, g_y_0_yzzzzz_xxz, g_y_0_yzzzzz_xyy, g_y_0_yzzzzz_xyz, g_y_0_yzzzzz_xzz, g_y_0_yzzzzz_yyy, g_y_0_yzzzzz_yyz, g_y_0_yzzzzz_yzz, g_y_0_yzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_xxx[k] = -g_y_0_yzzzz_xxx[k] * ab_z + g_y_0_yzzzz_xxxz[k];

                g_y_0_yzzzzz_xxy[k] = -g_y_0_yzzzz_xxy[k] * ab_z + g_y_0_yzzzz_xxyz[k];

                g_y_0_yzzzzz_xxz[k] = -g_y_0_yzzzz_xxz[k] * ab_z + g_y_0_yzzzz_xxzz[k];

                g_y_0_yzzzzz_xyy[k] = -g_y_0_yzzzz_xyy[k] * ab_z + g_y_0_yzzzz_xyyz[k];

                g_y_0_yzzzzz_xyz[k] = -g_y_0_yzzzz_xyz[k] * ab_z + g_y_0_yzzzz_xyzz[k];

                g_y_0_yzzzzz_xzz[k] = -g_y_0_yzzzz_xzz[k] * ab_z + g_y_0_yzzzz_xzzz[k];

                g_y_0_yzzzzz_yyy[k] = -g_y_0_yzzzz_yyy[k] * ab_z + g_y_0_yzzzz_yyyz[k];

                g_y_0_yzzzzz_yyz[k] = -g_y_0_yzzzz_yyz[k] * ab_z + g_y_0_yzzzz_yyzz[k];

                g_y_0_yzzzzz_yzz[k] = -g_y_0_yzzzz_yzz[k] * ab_z + g_y_0_yzzzz_yzzz[k];

                g_y_0_yzzzzz_zzz[k] = -g_y_0_yzzzz_zzz[k] * ab_z + g_y_0_yzzzz_zzzz[k];
            }

            /// Set up 550-560 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_xxx = cbuffer.data(if_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxy = cbuffer.data(if_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxz = cbuffer.data(if_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyy = cbuffer.data(if_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyz = cbuffer.data(if_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xzz = cbuffer.data(if_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyy = cbuffer.data(if_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyz = cbuffer.data(if_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yzz = cbuffer.data(if_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_zzzzzz_zzz = cbuffer.data(if_geom_10_off + 559 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzzz_xxx, g_y_0_zzzzz_xxxz, g_y_0_zzzzz_xxy, g_y_0_zzzzz_xxyz, g_y_0_zzzzz_xxz, g_y_0_zzzzz_xxzz, g_y_0_zzzzz_xyy, g_y_0_zzzzz_xyyz, g_y_0_zzzzz_xyz, g_y_0_zzzzz_xyzz, g_y_0_zzzzz_xzz, g_y_0_zzzzz_xzzz, g_y_0_zzzzz_yyy, g_y_0_zzzzz_yyyz, g_y_0_zzzzz_yyz, g_y_0_zzzzz_yyzz, g_y_0_zzzzz_yzz, g_y_0_zzzzz_yzzz, g_y_0_zzzzz_zzz, g_y_0_zzzzz_zzzz, g_y_0_zzzzzz_xxx, g_y_0_zzzzzz_xxy, g_y_0_zzzzzz_xxz, g_y_0_zzzzzz_xyy, g_y_0_zzzzzz_xyz, g_y_0_zzzzzz_xzz, g_y_0_zzzzzz_yyy, g_y_0_zzzzzz_yyz, g_y_0_zzzzzz_yzz, g_y_0_zzzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_xxx[k] = -g_y_0_zzzzz_xxx[k] * ab_z + g_y_0_zzzzz_xxxz[k];

                g_y_0_zzzzzz_xxy[k] = -g_y_0_zzzzz_xxy[k] * ab_z + g_y_0_zzzzz_xxyz[k];

                g_y_0_zzzzzz_xxz[k] = -g_y_0_zzzzz_xxz[k] * ab_z + g_y_0_zzzzz_xxzz[k];

                g_y_0_zzzzzz_xyy[k] = -g_y_0_zzzzz_xyy[k] * ab_z + g_y_0_zzzzz_xyyz[k];

                g_y_0_zzzzzz_xyz[k] = -g_y_0_zzzzz_xyz[k] * ab_z + g_y_0_zzzzz_xyzz[k];

                g_y_0_zzzzzz_xzz[k] = -g_y_0_zzzzz_xzz[k] * ab_z + g_y_0_zzzzz_xzzz[k];

                g_y_0_zzzzzz_yyy[k] = -g_y_0_zzzzz_yyy[k] * ab_z + g_y_0_zzzzz_yyyz[k];

                g_y_0_zzzzzz_yyz[k] = -g_y_0_zzzzz_yyz[k] * ab_z + g_y_0_zzzzz_yyzz[k];

                g_y_0_zzzzzz_yzz[k] = -g_y_0_zzzzz_yzz[k] * ab_z + g_y_0_zzzzz_yzzz[k];

                g_y_0_zzzzzz_zzz[k] = -g_y_0_zzzzz_zzz[k] * ab_z + g_y_0_zzzzz_zzzz[k];
            }

            /// Set up 560-570 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_xxx = cbuffer.data(if_geom_10_off + 560 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxy = cbuffer.data(if_geom_10_off + 561 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxz = cbuffer.data(if_geom_10_off + 562 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyy = cbuffer.data(if_geom_10_off + 563 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyz = cbuffer.data(if_geom_10_off + 564 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xzz = cbuffer.data(if_geom_10_off + 565 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyy = cbuffer.data(if_geom_10_off + 566 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyz = cbuffer.data(if_geom_10_off + 567 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yzz = cbuffer.data(if_geom_10_off + 568 * ccomps * dcomps);

            auto g_z_0_xxxxxx_zzz = cbuffer.data(if_geom_10_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxx_xxx, g_z_0_xxxxx_xxxx, g_z_0_xxxxx_xxxy, g_z_0_xxxxx_xxxz, g_z_0_xxxxx_xxy, g_z_0_xxxxx_xxyy, g_z_0_xxxxx_xxyz, g_z_0_xxxxx_xxz, g_z_0_xxxxx_xxzz, g_z_0_xxxxx_xyy, g_z_0_xxxxx_xyyy, g_z_0_xxxxx_xyyz, g_z_0_xxxxx_xyz, g_z_0_xxxxx_xyzz, g_z_0_xxxxx_xzz, g_z_0_xxxxx_xzzz, g_z_0_xxxxx_yyy, g_z_0_xxxxx_yyz, g_z_0_xxxxx_yzz, g_z_0_xxxxx_zzz, g_z_0_xxxxxx_xxx, g_z_0_xxxxxx_xxy, g_z_0_xxxxxx_xxz, g_z_0_xxxxxx_xyy, g_z_0_xxxxxx_xyz, g_z_0_xxxxxx_xzz, g_z_0_xxxxxx_yyy, g_z_0_xxxxxx_yyz, g_z_0_xxxxxx_yzz, g_z_0_xxxxxx_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_xxx[k] = -g_z_0_xxxxx_xxx[k] * ab_x + g_z_0_xxxxx_xxxx[k];

                g_z_0_xxxxxx_xxy[k] = -g_z_0_xxxxx_xxy[k] * ab_x + g_z_0_xxxxx_xxxy[k];

                g_z_0_xxxxxx_xxz[k] = -g_z_0_xxxxx_xxz[k] * ab_x + g_z_0_xxxxx_xxxz[k];

                g_z_0_xxxxxx_xyy[k] = -g_z_0_xxxxx_xyy[k] * ab_x + g_z_0_xxxxx_xxyy[k];

                g_z_0_xxxxxx_xyz[k] = -g_z_0_xxxxx_xyz[k] * ab_x + g_z_0_xxxxx_xxyz[k];

                g_z_0_xxxxxx_xzz[k] = -g_z_0_xxxxx_xzz[k] * ab_x + g_z_0_xxxxx_xxzz[k];

                g_z_0_xxxxxx_yyy[k] = -g_z_0_xxxxx_yyy[k] * ab_x + g_z_0_xxxxx_xyyy[k];

                g_z_0_xxxxxx_yyz[k] = -g_z_0_xxxxx_yyz[k] * ab_x + g_z_0_xxxxx_xyyz[k];

                g_z_0_xxxxxx_yzz[k] = -g_z_0_xxxxx_yzz[k] * ab_x + g_z_0_xxxxx_xyzz[k];

                g_z_0_xxxxxx_zzz[k] = -g_z_0_xxxxx_zzz[k] * ab_x + g_z_0_xxxxx_xzzz[k];
            }

            /// Set up 570-580 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_xxx = cbuffer.data(if_geom_10_off + 570 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxy = cbuffer.data(if_geom_10_off + 571 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxz = cbuffer.data(if_geom_10_off + 572 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyy = cbuffer.data(if_geom_10_off + 573 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyz = cbuffer.data(if_geom_10_off + 574 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xzz = cbuffer.data(if_geom_10_off + 575 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyy = cbuffer.data(if_geom_10_off + 576 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyz = cbuffer.data(if_geom_10_off + 577 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yzz = cbuffer.data(if_geom_10_off + 578 * ccomps * dcomps);

            auto g_z_0_xxxxxy_zzz = cbuffer.data(if_geom_10_off + 579 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxxy_xxx, g_z_0_xxxxxy_xxy, g_z_0_xxxxxy_xxz, g_z_0_xxxxxy_xyy, g_z_0_xxxxxy_xyz, g_z_0_xxxxxy_xzz, g_z_0_xxxxxy_yyy, g_z_0_xxxxxy_yyz, g_z_0_xxxxxy_yzz, g_z_0_xxxxxy_zzz, g_z_0_xxxxy_xxx, g_z_0_xxxxy_xxxx, g_z_0_xxxxy_xxxy, g_z_0_xxxxy_xxxz, g_z_0_xxxxy_xxy, g_z_0_xxxxy_xxyy, g_z_0_xxxxy_xxyz, g_z_0_xxxxy_xxz, g_z_0_xxxxy_xxzz, g_z_0_xxxxy_xyy, g_z_0_xxxxy_xyyy, g_z_0_xxxxy_xyyz, g_z_0_xxxxy_xyz, g_z_0_xxxxy_xyzz, g_z_0_xxxxy_xzz, g_z_0_xxxxy_xzzz, g_z_0_xxxxy_yyy, g_z_0_xxxxy_yyz, g_z_0_xxxxy_yzz, g_z_0_xxxxy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_xxx[k] = -g_z_0_xxxxy_xxx[k] * ab_x + g_z_0_xxxxy_xxxx[k];

                g_z_0_xxxxxy_xxy[k] = -g_z_0_xxxxy_xxy[k] * ab_x + g_z_0_xxxxy_xxxy[k];

                g_z_0_xxxxxy_xxz[k] = -g_z_0_xxxxy_xxz[k] * ab_x + g_z_0_xxxxy_xxxz[k];

                g_z_0_xxxxxy_xyy[k] = -g_z_0_xxxxy_xyy[k] * ab_x + g_z_0_xxxxy_xxyy[k];

                g_z_0_xxxxxy_xyz[k] = -g_z_0_xxxxy_xyz[k] * ab_x + g_z_0_xxxxy_xxyz[k];

                g_z_0_xxxxxy_xzz[k] = -g_z_0_xxxxy_xzz[k] * ab_x + g_z_0_xxxxy_xxzz[k];

                g_z_0_xxxxxy_yyy[k] = -g_z_0_xxxxy_yyy[k] * ab_x + g_z_0_xxxxy_xyyy[k];

                g_z_0_xxxxxy_yyz[k] = -g_z_0_xxxxy_yyz[k] * ab_x + g_z_0_xxxxy_xyyz[k];

                g_z_0_xxxxxy_yzz[k] = -g_z_0_xxxxy_yzz[k] * ab_x + g_z_0_xxxxy_xyzz[k];

                g_z_0_xxxxxy_zzz[k] = -g_z_0_xxxxy_zzz[k] * ab_x + g_z_0_xxxxy_xzzz[k];
            }

            /// Set up 580-590 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_xxx = cbuffer.data(if_geom_10_off + 580 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxy = cbuffer.data(if_geom_10_off + 581 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxz = cbuffer.data(if_geom_10_off + 582 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyy = cbuffer.data(if_geom_10_off + 583 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyz = cbuffer.data(if_geom_10_off + 584 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xzz = cbuffer.data(if_geom_10_off + 585 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyy = cbuffer.data(if_geom_10_off + 586 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyz = cbuffer.data(if_geom_10_off + 587 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yzz = cbuffer.data(if_geom_10_off + 588 * ccomps * dcomps);

            auto g_z_0_xxxxxz_zzz = cbuffer.data(if_geom_10_off + 589 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxxz_xxx, g_z_0_xxxxxz_xxy, g_z_0_xxxxxz_xxz, g_z_0_xxxxxz_xyy, g_z_0_xxxxxz_xyz, g_z_0_xxxxxz_xzz, g_z_0_xxxxxz_yyy, g_z_0_xxxxxz_yyz, g_z_0_xxxxxz_yzz, g_z_0_xxxxxz_zzz, g_z_0_xxxxz_xxx, g_z_0_xxxxz_xxxx, g_z_0_xxxxz_xxxy, g_z_0_xxxxz_xxxz, g_z_0_xxxxz_xxy, g_z_0_xxxxz_xxyy, g_z_0_xxxxz_xxyz, g_z_0_xxxxz_xxz, g_z_0_xxxxz_xxzz, g_z_0_xxxxz_xyy, g_z_0_xxxxz_xyyy, g_z_0_xxxxz_xyyz, g_z_0_xxxxz_xyz, g_z_0_xxxxz_xyzz, g_z_0_xxxxz_xzz, g_z_0_xxxxz_xzzz, g_z_0_xxxxz_yyy, g_z_0_xxxxz_yyz, g_z_0_xxxxz_yzz, g_z_0_xxxxz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_xxx[k] = -g_z_0_xxxxz_xxx[k] * ab_x + g_z_0_xxxxz_xxxx[k];

                g_z_0_xxxxxz_xxy[k] = -g_z_0_xxxxz_xxy[k] * ab_x + g_z_0_xxxxz_xxxy[k];

                g_z_0_xxxxxz_xxz[k] = -g_z_0_xxxxz_xxz[k] * ab_x + g_z_0_xxxxz_xxxz[k];

                g_z_0_xxxxxz_xyy[k] = -g_z_0_xxxxz_xyy[k] * ab_x + g_z_0_xxxxz_xxyy[k];

                g_z_0_xxxxxz_xyz[k] = -g_z_0_xxxxz_xyz[k] * ab_x + g_z_0_xxxxz_xxyz[k];

                g_z_0_xxxxxz_xzz[k] = -g_z_0_xxxxz_xzz[k] * ab_x + g_z_0_xxxxz_xxzz[k];

                g_z_0_xxxxxz_yyy[k] = -g_z_0_xxxxz_yyy[k] * ab_x + g_z_0_xxxxz_xyyy[k];

                g_z_0_xxxxxz_yyz[k] = -g_z_0_xxxxz_yyz[k] * ab_x + g_z_0_xxxxz_xyyz[k];

                g_z_0_xxxxxz_yzz[k] = -g_z_0_xxxxz_yzz[k] * ab_x + g_z_0_xxxxz_xyzz[k];

                g_z_0_xxxxxz_zzz[k] = -g_z_0_xxxxz_zzz[k] * ab_x + g_z_0_xxxxz_xzzz[k];
            }

            /// Set up 590-600 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_xxx = cbuffer.data(if_geom_10_off + 590 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxy = cbuffer.data(if_geom_10_off + 591 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxz = cbuffer.data(if_geom_10_off + 592 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyy = cbuffer.data(if_geom_10_off + 593 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyz = cbuffer.data(if_geom_10_off + 594 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xzz = cbuffer.data(if_geom_10_off + 595 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyy = cbuffer.data(if_geom_10_off + 596 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyz = cbuffer.data(if_geom_10_off + 597 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yzz = cbuffer.data(if_geom_10_off + 598 * ccomps * dcomps);

            auto g_z_0_xxxxyy_zzz = cbuffer.data(if_geom_10_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxyy_xxx, g_z_0_xxxxyy_xxy, g_z_0_xxxxyy_xxz, g_z_0_xxxxyy_xyy, g_z_0_xxxxyy_xyz, g_z_0_xxxxyy_xzz, g_z_0_xxxxyy_yyy, g_z_0_xxxxyy_yyz, g_z_0_xxxxyy_yzz, g_z_0_xxxxyy_zzz, g_z_0_xxxyy_xxx, g_z_0_xxxyy_xxxx, g_z_0_xxxyy_xxxy, g_z_0_xxxyy_xxxz, g_z_0_xxxyy_xxy, g_z_0_xxxyy_xxyy, g_z_0_xxxyy_xxyz, g_z_0_xxxyy_xxz, g_z_0_xxxyy_xxzz, g_z_0_xxxyy_xyy, g_z_0_xxxyy_xyyy, g_z_0_xxxyy_xyyz, g_z_0_xxxyy_xyz, g_z_0_xxxyy_xyzz, g_z_0_xxxyy_xzz, g_z_0_xxxyy_xzzz, g_z_0_xxxyy_yyy, g_z_0_xxxyy_yyz, g_z_0_xxxyy_yzz, g_z_0_xxxyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_xxx[k] = -g_z_0_xxxyy_xxx[k] * ab_x + g_z_0_xxxyy_xxxx[k];

                g_z_0_xxxxyy_xxy[k] = -g_z_0_xxxyy_xxy[k] * ab_x + g_z_0_xxxyy_xxxy[k];

                g_z_0_xxxxyy_xxz[k] = -g_z_0_xxxyy_xxz[k] * ab_x + g_z_0_xxxyy_xxxz[k];

                g_z_0_xxxxyy_xyy[k] = -g_z_0_xxxyy_xyy[k] * ab_x + g_z_0_xxxyy_xxyy[k];

                g_z_0_xxxxyy_xyz[k] = -g_z_0_xxxyy_xyz[k] * ab_x + g_z_0_xxxyy_xxyz[k];

                g_z_0_xxxxyy_xzz[k] = -g_z_0_xxxyy_xzz[k] * ab_x + g_z_0_xxxyy_xxzz[k];

                g_z_0_xxxxyy_yyy[k] = -g_z_0_xxxyy_yyy[k] * ab_x + g_z_0_xxxyy_xyyy[k];

                g_z_0_xxxxyy_yyz[k] = -g_z_0_xxxyy_yyz[k] * ab_x + g_z_0_xxxyy_xyyz[k];

                g_z_0_xxxxyy_yzz[k] = -g_z_0_xxxyy_yzz[k] * ab_x + g_z_0_xxxyy_xyzz[k];

                g_z_0_xxxxyy_zzz[k] = -g_z_0_xxxyy_zzz[k] * ab_x + g_z_0_xxxyy_xzzz[k];
            }

            /// Set up 600-610 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_xxx = cbuffer.data(if_geom_10_off + 600 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxy = cbuffer.data(if_geom_10_off + 601 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxz = cbuffer.data(if_geom_10_off + 602 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyy = cbuffer.data(if_geom_10_off + 603 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyz = cbuffer.data(if_geom_10_off + 604 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xzz = cbuffer.data(if_geom_10_off + 605 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyy = cbuffer.data(if_geom_10_off + 606 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyz = cbuffer.data(if_geom_10_off + 607 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yzz = cbuffer.data(if_geom_10_off + 608 * ccomps * dcomps);

            auto g_z_0_xxxxyz_zzz = cbuffer.data(if_geom_10_off + 609 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxyz_xxx, g_z_0_xxxxyz_xxy, g_z_0_xxxxyz_xxz, g_z_0_xxxxyz_xyy, g_z_0_xxxxyz_xyz, g_z_0_xxxxyz_xzz, g_z_0_xxxxyz_yyy, g_z_0_xxxxyz_yyz, g_z_0_xxxxyz_yzz, g_z_0_xxxxyz_zzz, g_z_0_xxxyz_xxx, g_z_0_xxxyz_xxxx, g_z_0_xxxyz_xxxy, g_z_0_xxxyz_xxxz, g_z_0_xxxyz_xxy, g_z_0_xxxyz_xxyy, g_z_0_xxxyz_xxyz, g_z_0_xxxyz_xxz, g_z_0_xxxyz_xxzz, g_z_0_xxxyz_xyy, g_z_0_xxxyz_xyyy, g_z_0_xxxyz_xyyz, g_z_0_xxxyz_xyz, g_z_0_xxxyz_xyzz, g_z_0_xxxyz_xzz, g_z_0_xxxyz_xzzz, g_z_0_xxxyz_yyy, g_z_0_xxxyz_yyz, g_z_0_xxxyz_yzz, g_z_0_xxxyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_xxx[k] = -g_z_0_xxxyz_xxx[k] * ab_x + g_z_0_xxxyz_xxxx[k];

                g_z_0_xxxxyz_xxy[k] = -g_z_0_xxxyz_xxy[k] * ab_x + g_z_0_xxxyz_xxxy[k];

                g_z_0_xxxxyz_xxz[k] = -g_z_0_xxxyz_xxz[k] * ab_x + g_z_0_xxxyz_xxxz[k];

                g_z_0_xxxxyz_xyy[k] = -g_z_0_xxxyz_xyy[k] * ab_x + g_z_0_xxxyz_xxyy[k];

                g_z_0_xxxxyz_xyz[k] = -g_z_0_xxxyz_xyz[k] * ab_x + g_z_0_xxxyz_xxyz[k];

                g_z_0_xxxxyz_xzz[k] = -g_z_0_xxxyz_xzz[k] * ab_x + g_z_0_xxxyz_xxzz[k];

                g_z_0_xxxxyz_yyy[k] = -g_z_0_xxxyz_yyy[k] * ab_x + g_z_0_xxxyz_xyyy[k];

                g_z_0_xxxxyz_yyz[k] = -g_z_0_xxxyz_yyz[k] * ab_x + g_z_0_xxxyz_xyyz[k];

                g_z_0_xxxxyz_yzz[k] = -g_z_0_xxxyz_yzz[k] * ab_x + g_z_0_xxxyz_xyzz[k];

                g_z_0_xxxxyz_zzz[k] = -g_z_0_xxxyz_zzz[k] * ab_x + g_z_0_xxxyz_xzzz[k];
            }

            /// Set up 610-620 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_xxx = cbuffer.data(if_geom_10_off + 610 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxy = cbuffer.data(if_geom_10_off + 611 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxz = cbuffer.data(if_geom_10_off + 612 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyy = cbuffer.data(if_geom_10_off + 613 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyz = cbuffer.data(if_geom_10_off + 614 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xzz = cbuffer.data(if_geom_10_off + 615 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyy = cbuffer.data(if_geom_10_off + 616 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyz = cbuffer.data(if_geom_10_off + 617 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yzz = cbuffer.data(if_geom_10_off + 618 * ccomps * dcomps);

            auto g_z_0_xxxxzz_zzz = cbuffer.data(if_geom_10_off + 619 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxzz_xxx, g_z_0_xxxxzz_xxy, g_z_0_xxxxzz_xxz, g_z_0_xxxxzz_xyy, g_z_0_xxxxzz_xyz, g_z_0_xxxxzz_xzz, g_z_0_xxxxzz_yyy, g_z_0_xxxxzz_yyz, g_z_0_xxxxzz_yzz, g_z_0_xxxxzz_zzz, g_z_0_xxxzz_xxx, g_z_0_xxxzz_xxxx, g_z_0_xxxzz_xxxy, g_z_0_xxxzz_xxxz, g_z_0_xxxzz_xxy, g_z_0_xxxzz_xxyy, g_z_0_xxxzz_xxyz, g_z_0_xxxzz_xxz, g_z_0_xxxzz_xxzz, g_z_0_xxxzz_xyy, g_z_0_xxxzz_xyyy, g_z_0_xxxzz_xyyz, g_z_0_xxxzz_xyz, g_z_0_xxxzz_xyzz, g_z_0_xxxzz_xzz, g_z_0_xxxzz_xzzz, g_z_0_xxxzz_yyy, g_z_0_xxxzz_yyz, g_z_0_xxxzz_yzz, g_z_0_xxxzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_xxx[k] = -g_z_0_xxxzz_xxx[k] * ab_x + g_z_0_xxxzz_xxxx[k];

                g_z_0_xxxxzz_xxy[k] = -g_z_0_xxxzz_xxy[k] * ab_x + g_z_0_xxxzz_xxxy[k];

                g_z_0_xxxxzz_xxz[k] = -g_z_0_xxxzz_xxz[k] * ab_x + g_z_0_xxxzz_xxxz[k];

                g_z_0_xxxxzz_xyy[k] = -g_z_0_xxxzz_xyy[k] * ab_x + g_z_0_xxxzz_xxyy[k];

                g_z_0_xxxxzz_xyz[k] = -g_z_0_xxxzz_xyz[k] * ab_x + g_z_0_xxxzz_xxyz[k];

                g_z_0_xxxxzz_xzz[k] = -g_z_0_xxxzz_xzz[k] * ab_x + g_z_0_xxxzz_xxzz[k];

                g_z_0_xxxxzz_yyy[k] = -g_z_0_xxxzz_yyy[k] * ab_x + g_z_0_xxxzz_xyyy[k];

                g_z_0_xxxxzz_yyz[k] = -g_z_0_xxxzz_yyz[k] * ab_x + g_z_0_xxxzz_xyyz[k];

                g_z_0_xxxxzz_yzz[k] = -g_z_0_xxxzz_yzz[k] * ab_x + g_z_0_xxxzz_xyzz[k];

                g_z_0_xxxxzz_zzz[k] = -g_z_0_xxxzz_zzz[k] * ab_x + g_z_0_xxxzz_xzzz[k];
            }

            /// Set up 620-630 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_xxx = cbuffer.data(if_geom_10_off + 620 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxy = cbuffer.data(if_geom_10_off + 621 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxz = cbuffer.data(if_geom_10_off + 622 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyy = cbuffer.data(if_geom_10_off + 623 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyz = cbuffer.data(if_geom_10_off + 624 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xzz = cbuffer.data(if_geom_10_off + 625 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyy = cbuffer.data(if_geom_10_off + 626 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyz = cbuffer.data(if_geom_10_off + 627 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yzz = cbuffer.data(if_geom_10_off + 628 * ccomps * dcomps);

            auto g_z_0_xxxyyy_zzz = cbuffer.data(if_geom_10_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyyy_xxx, g_z_0_xxxyyy_xxy, g_z_0_xxxyyy_xxz, g_z_0_xxxyyy_xyy, g_z_0_xxxyyy_xyz, g_z_0_xxxyyy_xzz, g_z_0_xxxyyy_yyy, g_z_0_xxxyyy_yyz, g_z_0_xxxyyy_yzz, g_z_0_xxxyyy_zzz, g_z_0_xxyyy_xxx, g_z_0_xxyyy_xxxx, g_z_0_xxyyy_xxxy, g_z_0_xxyyy_xxxz, g_z_0_xxyyy_xxy, g_z_0_xxyyy_xxyy, g_z_0_xxyyy_xxyz, g_z_0_xxyyy_xxz, g_z_0_xxyyy_xxzz, g_z_0_xxyyy_xyy, g_z_0_xxyyy_xyyy, g_z_0_xxyyy_xyyz, g_z_0_xxyyy_xyz, g_z_0_xxyyy_xyzz, g_z_0_xxyyy_xzz, g_z_0_xxyyy_xzzz, g_z_0_xxyyy_yyy, g_z_0_xxyyy_yyz, g_z_0_xxyyy_yzz, g_z_0_xxyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_xxx[k] = -g_z_0_xxyyy_xxx[k] * ab_x + g_z_0_xxyyy_xxxx[k];

                g_z_0_xxxyyy_xxy[k] = -g_z_0_xxyyy_xxy[k] * ab_x + g_z_0_xxyyy_xxxy[k];

                g_z_0_xxxyyy_xxz[k] = -g_z_0_xxyyy_xxz[k] * ab_x + g_z_0_xxyyy_xxxz[k];

                g_z_0_xxxyyy_xyy[k] = -g_z_0_xxyyy_xyy[k] * ab_x + g_z_0_xxyyy_xxyy[k];

                g_z_0_xxxyyy_xyz[k] = -g_z_0_xxyyy_xyz[k] * ab_x + g_z_0_xxyyy_xxyz[k];

                g_z_0_xxxyyy_xzz[k] = -g_z_0_xxyyy_xzz[k] * ab_x + g_z_0_xxyyy_xxzz[k];

                g_z_0_xxxyyy_yyy[k] = -g_z_0_xxyyy_yyy[k] * ab_x + g_z_0_xxyyy_xyyy[k];

                g_z_0_xxxyyy_yyz[k] = -g_z_0_xxyyy_yyz[k] * ab_x + g_z_0_xxyyy_xyyz[k];

                g_z_0_xxxyyy_yzz[k] = -g_z_0_xxyyy_yzz[k] * ab_x + g_z_0_xxyyy_xyzz[k];

                g_z_0_xxxyyy_zzz[k] = -g_z_0_xxyyy_zzz[k] * ab_x + g_z_0_xxyyy_xzzz[k];
            }

            /// Set up 630-640 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_xxx = cbuffer.data(if_geom_10_off + 630 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxy = cbuffer.data(if_geom_10_off + 631 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxz = cbuffer.data(if_geom_10_off + 632 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyy = cbuffer.data(if_geom_10_off + 633 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyz = cbuffer.data(if_geom_10_off + 634 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xzz = cbuffer.data(if_geom_10_off + 635 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyy = cbuffer.data(if_geom_10_off + 636 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyz = cbuffer.data(if_geom_10_off + 637 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yzz = cbuffer.data(if_geom_10_off + 638 * ccomps * dcomps);

            auto g_z_0_xxxyyz_zzz = cbuffer.data(if_geom_10_off + 639 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyyz_xxx, g_z_0_xxxyyz_xxy, g_z_0_xxxyyz_xxz, g_z_0_xxxyyz_xyy, g_z_0_xxxyyz_xyz, g_z_0_xxxyyz_xzz, g_z_0_xxxyyz_yyy, g_z_0_xxxyyz_yyz, g_z_0_xxxyyz_yzz, g_z_0_xxxyyz_zzz, g_z_0_xxyyz_xxx, g_z_0_xxyyz_xxxx, g_z_0_xxyyz_xxxy, g_z_0_xxyyz_xxxz, g_z_0_xxyyz_xxy, g_z_0_xxyyz_xxyy, g_z_0_xxyyz_xxyz, g_z_0_xxyyz_xxz, g_z_0_xxyyz_xxzz, g_z_0_xxyyz_xyy, g_z_0_xxyyz_xyyy, g_z_0_xxyyz_xyyz, g_z_0_xxyyz_xyz, g_z_0_xxyyz_xyzz, g_z_0_xxyyz_xzz, g_z_0_xxyyz_xzzz, g_z_0_xxyyz_yyy, g_z_0_xxyyz_yyz, g_z_0_xxyyz_yzz, g_z_0_xxyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_xxx[k] = -g_z_0_xxyyz_xxx[k] * ab_x + g_z_0_xxyyz_xxxx[k];

                g_z_0_xxxyyz_xxy[k] = -g_z_0_xxyyz_xxy[k] * ab_x + g_z_0_xxyyz_xxxy[k];

                g_z_0_xxxyyz_xxz[k] = -g_z_0_xxyyz_xxz[k] * ab_x + g_z_0_xxyyz_xxxz[k];

                g_z_0_xxxyyz_xyy[k] = -g_z_0_xxyyz_xyy[k] * ab_x + g_z_0_xxyyz_xxyy[k];

                g_z_0_xxxyyz_xyz[k] = -g_z_0_xxyyz_xyz[k] * ab_x + g_z_0_xxyyz_xxyz[k];

                g_z_0_xxxyyz_xzz[k] = -g_z_0_xxyyz_xzz[k] * ab_x + g_z_0_xxyyz_xxzz[k];

                g_z_0_xxxyyz_yyy[k] = -g_z_0_xxyyz_yyy[k] * ab_x + g_z_0_xxyyz_xyyy[k];

                g_z_0_xxxyyz_yyz[k] = -g_z_0_xxyyz_yyz[k] * ab_x + g_z_0_xxyyz_xyyz[k];

                g_z_0_xxxyyz_yzz[k] = -g_z_0_xxyyz_yzz[k] * ab_x + g_z_0_xxyyz_xyzz[k];

                g_z_0_xxxyyz_zzz[k] = -g_z_0_xxyyz_zzz[k] * ab_x + g_z_0_xxyyz_xzzz[k];
            }

            /// Set up 640-650 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_xxx = cbuffer.data(if_geom_10_off + 640 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxy = cbuffer.data(if_geom_10_off + 641 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxz = cbuffer.data(if_geom_10_off + 642 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyy = cbuffer.data(if_geom_10_off + 643 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyz = cbuffer.data(if_geom_10_off + 644 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xzz = cbuffer.data(if_geom_10_off + 645 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyy = cbuffer.data(if_geom_10_off + 646 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyz = cbuffer.data(if_geom_10_off + 647 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yzz = cbuffer.data(if_geom_10_off + 648 * ccomps * dcomps);

            auto g_z_0_xxxyzz_zzz = cbuffer.data(if_geom_10_off + 649 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyzz_xxx, g_z_0_xxxyzz_xxy, g_z_0_xxxyzz_xxz, g_z_0_xxxyzz_xyy, g_z_0_xxxyzz_xyz, g_z_0_xxxyzz_xzz, g_z_0_xxxyzz_yyy, g_z_0_xxxyzz_yyz, g_z_0_xxxyzz_yzz, g_z_0_xxxyzz_zzz, g_z_0_xxyzz_xxx, g_z_0_xxyzz_xxxx, g_z_0_xxyzz_xxxy, g_z_0_xxyzz_xxxz, g_z_0_xxyzz_xxy, g_z_0_xxyzz_xxyy, g_z_0_xxyzz_xxyz, g_z_0_xxyzz_xxz, g_z_0_xxyzz_xxzz, g_z_0_xxyzz_xyy, g_z_0_xxyzz_xyyy, g_z_0_xxyzz_xyyz, g_z_0_xxyzz_xyz, g_z_0_xxyzz_xyzz, g_z_0_xxyzz_xzz, g_z_0_xxyzz_xzzz, g_z_0_xxyzz_yyy, g_z_0_xxyzz_yyz, g_z_0_xxyzz_yzz, g_z_0_xxyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_xxx[k] = -g_z_0_xxyzz_xxx[k] * ab_x + g_z_0_xxyzz_xxxx[k];

                g_z_0_xxxyzz_xxy[k] = -g_z_0_xxyzz_xxy[k] * ab_x + g_z_0_xxyzz_xxxy[k];

                g_z_0_xxxyzz_xxz[k] = -g_z_0_xxyzz_xxz[k] * ab_x + g_z_0_xxyzz_xxxz[k];

                g_z_0_xxxyzz_xyy[k] = -g_z_0_xxyzz_xyy[k] * ab_x + g_z_0_xxyzz_xxyy[k];

                g_z_0_xxxyzz_xyz[k] = -g_z_0_xxyzz_xyz[k] * ab_x + g_z_0_xxyzz_xxyz[k];

                g_z_0_xxxyzz_xzz[k] = -g_z_0_xxyzz_xzz[k] * ab_x + g_z_0_xxyzz_xxzz[k];

                g_z_0_xxxyzz_yyy[k] = -g_z_0_xxyzz_yyy[k] * ab_x + g_z_0_xxyzz_xyyy[k];

                g_z_0_xxxyzz_yyz[k] = -g_z_0_xxyzz_yyz[k] * ab_x + g_z_0_xxyzz_xyyz[k];

                g_z_0_xxxyzz_yzz[k] = -g_z_0_xxyzz_yzz[k] * ab_x + g_z_0_xxyzz_xyzz[k];

                g_z_0_xxxyzz_zzz[k] = -g_z_0_xxyzz_zzz[k] * ab_x + g_z_0_xxyzz_xzzz[k];
            }

            /// Set up 650-660 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_xxx = cbuffer.data(if_geom_10_off + 650 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxy = cbuffer.data(if_geom_10_off + 651 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxz = cbuffer.data(if_geom_10_off + 652 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyy = cbuffer.data(if_geom_10_off + 653 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyz = cbuffer.data(if_geom_10_off + 654 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xzz = cbuffer.data(if_geom_10_off + 655 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyy = cbuffer.data(if_geom_10_off + 656 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyz = cbuffer.data(if_geom_10_off + 657 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yzz = cbuffer.data(if_geom_10_off + 658 * ccomps * dcomps);

            auto g_z_0_xxxzzz_zzz = cbuffer.data(if_geom_10_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxzzz_xxx, g_z_0_xxxzzz_xxy, g_z_0_xxxzzz_xxz, g_z_0_xxxzzz_xyy, g_z_0_xxxzzz_xyz, g_z_0_xxxzzz_xzz, g_z_0_xxxzzz_yyy, g_z_0_xxxzzz_yyz, g_z_0_xxxzzz_yzz, g_z_0_xxxzzz_zzz, g_z_0_xxzzz_xxx, g_z_0_xxzzz_xxxx, g_z_0_xxzzz_xxxy, g_z_0_xxzzz_xxxz, g_z_0_xxzzz_xxy, g_z_0_xxzzz_xxyy, g_z_0_xxzzz_xxyz, g_z_0_xxzzz_xxz, g_z_0_xxzzz_xxzz, g_z_0_xxzzz_xyy, g_z_0_xxzzz_xyyy, g_z_0_xxzzz_xyyz, g_z_0_xxzzz_xyz, g_z_0_xxzzz_xyzz, g_z_0_xxzzz_xzz, g_z_0_xxzzz_xzzz, g_z_0_xxzzz_yyy, g_z_0_xxzzz_yyz, g_z_0_xxzzz_yzz, g_z_0_xxzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_xxx[k] = -g_z_0_xxzzz_xxx[k] * ab_x + g_z_0_xxzzz_xxxx[k];

                g_z_0_xxxzzz_xxy[k] = -g_z_0_xxzzz_xxy[k] * ab_x + g_z_0_xxzzz_xxxy[k];

                g_z_0_xxxzzz_xxz[k] = -g_z_0_xxzzz_xxz[k] * ab_x + g_z_0_xxzzz_xxxz[k];

                g_z_0_xxxzzz_xyy[k] = -g_z_0_xxzzz_xyy[k] * ab_x + g_z_0_xxzzz_xxyy[k];

                g_z_0_xxxzzz_xyz[k] = -g_z_0_xxzzz_xyz[k] * ab_x + g_z_0_xxzzz_xxyz[k];

                g_z_0_xxxzzz_xzz[k] = -g_z_0_xxzzz_xzz[k] * ab_x + g_z_0_xxzzz_xxzz[k];

                g_z_0_xxxzzz_yyy[k] = -g_z_0_xxzzz_yyy[k] * ab_x + g_z_0_xxzzz_xyyy[k];

                g_z_0_xxxzzz_yyz[k] = -g_z_0_xxzzz_yyz[k] * ab_x + g_z_0_xxzzz_xyyz[k];

                g_z_0_xxxzzz_yzz[k] = -g_z_0_xxzzz_yzz[k] * ab_x + g_z_0_xxzzz_xyzz[k];

                g_z_0_xxxzzz_zzz[k] = -g_z_0_xxzzz_zzz[k] * ab_x + g_z_0_xxzzz_xzzz[k];
            }

            /// Set up 660-670 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_xxx = cbuffer.data(if_geom_10_off + 660 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxy = cbuffer.data(if_geom_10_off + 661 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxz = cbuffer.data(if_geom_10_off + 662 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyy = cbuffer.data(if_geom_10_off + 663 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyz = cbuffer.data(if_geom_10_off + 664 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xzz = cbuffer.data(if_geom_10_off + 665 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyy = cbuffer.data(if_geom_10_off + 666 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyz = cbuffer.data(if_geom_10_off + 667 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yzz = cbuffer.data(if_geom_10_off + 668 * ccomps * dcomps);

            auto g_z_0_xxyyyy_zzz = cbuffer.data(if_geom_10_off + 669 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyyy_xxx, g_z_0_xxyyyy_xxy, g_z_0_xxyyyy_xxz, g_z_0_xxyyyy_xyy, g_z_0_xxyyyy_xyz, g_z_0_xxyyyy_xzz, g_z_0_xxyyyy_yyy, g_z_0_xxyyyy_yyz, g_z_0_xxyyyy_yzz, g_z_0_xxyyyy_zzz, g_z_0_xyyyy_xxx, g_z_0_xyyyy_xxxx, g_z_0_xyyyy_xxxy, g_z_0_xyyyy_xxxz, g_z_0_xyyyy_xxy, g_z_0_xyyyy_xxyy, g_z_0_xyyyy_xxyz, g_z_0_xyyyy_xxz, g_z_0_xyyyy_xxzz, g_z_0_xyyyy_xyy, g_z_0_xyyyy_xyyy, g_z_0_xyyyy_xyyz, g_z_0_xyyyy_xyz, g_z_0_xyyyy_xyzz, g_z_0_xyyyy_xzz, g_z_0_xyyyy_xzzz, g_z_0_xyyyy_yyy, g_z_0_xyyyy_yyz, g_z_0_xyyyy_yzz, g_z_0_xyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_xxx[k] = -g_z_0_xyyyy_xxx[k] * ab_x + g_z_0_xyyyy_xxxx[k];

                g_z_0_xxyyyy_xxy[k] = -g_z_0_xyyyy_xxy[k] * ab_x + g_z_0_xyyyy_xxxy[k];

                g_z_0_xxyyyy_xxz[k] = -g_z_0_xyyyy_xxz[k] * ab_x + g_z_0_xyyyy_xxxz[k];

                g_z_0_xxyyyy_xyy[k] = -g_z_0_xyyyy_xyy[k] * ab_x + g_z_0_xyyyy_xxyy[k];

                g_z_0_xxyyyy_xyz[k] = -g_z_0_xyyyy_xyz[k] * ab_x + g_z_0_xyyyy_xxyz[k];

                g_z_0_xxyyyy_xzz[k] = -g_z_0_xyyyy_xzz[k] * ab_x + g_z_0_xyyyy_xxzz[k];

                g_z_0_xxyyyy_yyy[k] = -g_z_0_xyyyy_yyy[k] * ab_x + g_z_0_xyyyy_xyyy[k];

                g_z_0_xxyyyy_yyz[k] = -g_z_0_xyyyy_yyz[k] * ab_x + g_z_0_xyyyy_xyyz[k];

                g_z_0_xxyyyy_yzz[k] = -g_z_0_xyyyy_yzz[k] * ab_x + g_z_0_xyyyy_xyzz[k];

                g_z_0_xxyyyy_zzz[k] = -g_z_0_xyyyy_zzz[k] * ab_x + g_z_0_xyyyy_xzzz[k];
            }

            /// Set up 670-680 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_xxx = cbuffer.data(if_geom_10_off + 670 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxy = cbuffer.data(if_geom_10_off + 671 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxz = cbuffer.data(if_geom_10_off + 672 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyy = cbuffer.data(if_geom_10_off + 673 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyz = cbuffer.data(if_geom_10_off + 674 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xzz = cbuffer.data(if_geom_10_off + 675 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyy = cbuffer.data(if_geom_10_off + 676 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyz = cbuffer.data(if_geom_10_off + 677 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yzz = cbuffer.data(if_geom_10_off + 678 * ccomps * dcomps);

            auto g_z_0_xxyyyz_zzz = cbuffer.data(if_geom_10_off + 679 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyyz_xxx, g_z_0_xxyyyz_xxy, g_z_0_xxyyyz_xxz, g_z_0_xxyyyz_xyy, g_z_0_xxyyyz_xyz, g_z_0_xxyyyz_xzz, g_z_0_xxyyyz_yyy, g_z_0_xxyyyz_yyz, g_z_0_xxyyyz_yzz, g_z_0_xxyyyz_zzz, g_z_0_xyyyz_xxx, g_z_0_xyyyz_xxxx, g_z_0_xyyyz_xxxy, g_z_0_xyyyz_xxxz, g_z_0_xyyyz_xxy, g_z_0_xyyyz_xxyy, g_z_0_xyyyz_xxyz, g_z_0_xyyyz_xxz, g_z_0_xyyyz_xxzz, g_z_0_xyyyz_xyy, g_z_0_xyyyz_xyyy, g_z_0_xyyyz_xyyz, g_z_0_xyyyz_xyz, g_z_0_xyyyz_xyzz, g_z_0_xyyyz_xzz, g_z_0_xyyyz_xzzz, g_z_0_xyyyz_yyy, g_z_0_xyyyz_yyz, g_z_0_xyyyz_yzz, g_z_0_xyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_xxx[k] = -g_z_0_xyyyz_xxx[k] * ab_x + g_z_0_xyyyz_xxxx[k];

                g_z_0_xxyyyz_xxy[k] = -g_z_0_xyyyz_xxy[k] * ab_x + g_z_0_xyyyz_xxxy[k];

                g_z_0_xxyyyz_xxz[k] = -g_z_0_xyyyz_xxz[k] * ab_x + g_z_0_xyyyz_xxxz[k];

                g_z_0_xxyyyz_xyy[k] = -g_z_0_xyyyz_xyy[k] * ab_x + g_z_0_xyyyz_xxyy[k];

                g_z_0_xxyyyz_xyz[k] = -g_z_0_xyyyz_xyz[k] * ab_x + g_z_0_xyyyz_xxyz[k];

                g_z_0_xxyyyz_xzz[k] = -g_z_0_xyyyz_xzz[k] * ab_x + g_z_0_xyyyz_xxzz[k];

                g_z_0_xxyyyz_yyy[k] = -g_z_0_xyyyz_yyy[k] * ab_x + g_z_0_xyyyz_xyyy[k];

                g_z_0_xxyyyz_yyz[k] = -g_z_0_xyyyz_yyz[k] * ab_x + g_z_0_xyyyz_xyyz[k];

                g_z_0_xxyyyz_yzz[k] = -g_z_0_xyyyz_yzz[k] * ab_x + g_z_0_xyyyz_xyzz[k];

                g_z_0_xxyyyz_zzz[k] = -g_z_0_xyyyz_zzz[k] * ab_x + g_z_0_xyyyz_xzzz[k];
            }

            /// Set up 680-690 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_xxx = cbuffer.data(if_geom_10_off + 680 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxy = cbuffer.data(if_geom_10_off + 681 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxz = cbuffer.data(if_geom_10_off + 682 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyy = cbuffer.data(if_geom_10_off + 683 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyz = cbuffer.data(if_geom_10_off + 684 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xzz = cbuffer.data(if_geom_10_off + 685 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyy = cbuffer.data(if_geom_10_off + 686 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyz = cbuffer.data(if_geom_10_off + 687 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yzz = cbuffer.data(if_geom_10_off + 688 * ccomps * dcomps);

            auto g_z_0_xxyyzz_zzz = cbuffer.data(if_geom_10_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyzz_xxx, g_z_0_xxyyzz_xxy, g_z_0_xxyyzz_xxz, g_z_0_xxyyzz_xyy, g_z_0_xxyyzz_xyz, g_z_0_xxyyzz_xzz, g_z_0_xxyyzz_yyy, g_z_0_xxyyzz_yyz, g_z_0_xxyyzz_yzz, g_z_0_xxyyzz_zzz, g_z_0_xyyzz_xxx, g_z_0_xyyzz_xxxx, g_z_0_xyyzz_xxxy, g_z_0_xyyzz_xxxz, g_z_0_xyyzz_xxy, g_z_0_xyyzz_xxyy, g_z_0_xyyzz_xxyz, g_z_0_xyyzz_xxz, g_z_0_xyyzz_xxzz, g_z_0_xyyzz_xyy, g_z_0_xyyzz_xyyy, g_z_0_xyyzz_xyyz, g_z_0_xyyzz_xyz, g_z_0_xyyzz_xyzz, g_z_0_xyyzz_xzz, g_z_0_xyyzz_xzzz, g_z_0_xyyzz_yyy, g_z_0_xyyzz_yyz, g_z_0_xyyzz_yzz, g_z_0_xyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_xxx[k] = -g_z_0_xyyzz_xxx[k] * ab_x + g_z_0_xyyzz_xxxx[k];

                g_z_0_xxyyzz_xxy[k] = -g_z_0_xyyzz_xxy[k] * ab_x + g_z_0_xyyzz_xxxy[k];

                g_z_0_xxyyzz_xxz[k] = -g_z_0_xyyzz_xxz[k] * ab_x + g_z_0_xyyzz_xxxz[k];

                g_z_0_xxyyzz_xyy[k] = -g_z_0_xyyzz_xyy[k] * ab_x + g_z_0_xyyzz_xxyy[k];

                g_z_0_xxyyzz_xyz[k] = -g_z_0_xyyzz_xyz[k] * ab_x + g_z_0_xyyzz_xxyz[k];

                g_z_0_xxyyzz_xzz[k] = -g_z_0_xyyzz_xzz[k] * ab_x + g_z_0_xyyzz_xxzz[k];

                g_z_0_xxyyzz_yyy[k] = -g_z_0_xyyzz_yyy[k] * ab_x + g_z_0_xyyzz_xyyy[k];

                g_z_0_xxyyzz_yyz[k] = -g_z_0_xyyzz_yyz[k] * ab_x + g_z_0_xyyzz_xyyz[k];

                g_z_0_xxyyzz_yzz[k] = -g_z_0_xyyzz_yzz[k] * ab_x + g_z_0_xyyzz_xyzz[k];

                g_z_0_xxyyzz_zzz[k] = -g_z_0_xyyzz_zzz[k] * ab_x + g_z_0_xyyzz_xzzz[k];
            }

            /// Set up 690-700 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_xxx = cbuffer.data(if_geom_10_off + 690 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxy = cbuffer.data(if_geom_10_off + 691 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxz = cbuffer.data(if_geom_10_off + 692 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyy = cbuffer.data(if_geom_10_off + 693 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyz = cbuffer.data(if_geom_10_off + 694 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xzz = cbuffer.data(if_geom_10_off + 695 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyy = cbuffer.data(if_geom_10_off + 696 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyz = cbuffer.data(if_geom_10_off + 697 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yzz = cbuffer.data(if_geom_10_off + 698 * ccomps * dcomps);

            auto g_z_0_xxyzzz_zzz = cbuffer.data(if_geom_10_off + 699 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyzzz_xxx, g_z_0_xxyzzz_xxy, g_z_0_xxyzzz_xxz, g_z_0_xxyzzz_xyy, g_z_0_xxyzzz_xyz, g_z_0_xxyzzz_xzz, g_z_0_xxyzzz_yyy, g_z_0_xxyzzz_yyz, g_z_0_xxyzzz_yzz, g_z_0_xxyzzz_zzz, g_z_0_xyzzz_xxx, g_z_0_xyzzz_xxxx, g_z_0_xyzzz_xxxy, g_z_0_xyzzz_xxxz, g_z_0_xyzzz_xxy, g_z_0_xyzzz_xxyy, g_z_0_xyzzz_xxyz, g_z_0_xyzzz_xxz, g_z_0_xyzzz_xxzz, g_z_0_xyzzz_xyy, g_z_0_xyzzz_xyyy, g_z_0_xyzzz_xyyz, g_z_0_xyzzz_xyz, g_z_0_xyzzz_xyzz, g_z_0_xyzzz_xzz, g_z_0_xyzzz_xzzz, g_z_0_xyzzz_yyy, g_z_0_xyzzz_yyz, g_z_0_xyzzz_yzz, g_z_0_xyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_xxx[k] = -g_z_0_xyzzz_xxx[k] * ab_x + g_z_0_xyzzz_xxxx[k];

                g_z_0_xxyzzz_xxy[k] = -g_z_0_xyzzz_xxy[k] * ab_x + g_z_0_xyzzz_xxxy[k];

                g_z_0_xxyzzz_xxz[k] = -g_z_0_xyzzz_xxz[k] * ab_x + g_z_0_xyzzz_xxxz[k];

                g_z_0_xxyzzz_xyy[k] = -g_z_0_xyzzz_xyy[k] * ab_x + g_z_0_xyzzz_xxyy[k];

                g_z_0_xxyzzz_xyz[k] = -g_z_0_xyzzz_xyz[k] * ab_x + g_z_0_xyzzz_xxyz[k];

                g_z_0_xxyzzz_xzz[k] = -g_z_0_xyzzz_xzz[k] * ab_x + g_z_0_xyzzz_xxzz[k];

                g_z_0_xxyzzz_yyy[k] = -g_z_0_xyzzz_yyy[k] * ab_x + g_z_0_xyzzz_xyyy[k];

                g_z_0_xxyzzz_yyz[k] = -g_z_0_xyzzz_yyz[k] * ab_x + g_z_0_xyzzz_xyyz[k];

                g_z_0_xxyzzz_yzz[k] = -g_z_0_xyzzz_yzz[k] * ab_x + g_z_0_xyzzz_xyzz[k];

                g_z_0_xxyzzz_zzz[k] = -g_z_0_xyzzz_zzz[k] * ab_x + g_z_0_xyzzz_xzzz[k];
            }

            /// Set up 700-710 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_xxx = cbuffer.data(if_geom_10_off + 700 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxy = cbuffer.data(if_geom_10_off + 701 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxz = cbuffer.data(if_geom_10_off + 702 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyy = cbuffer.data(if_geom_10_off + 703 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyz = cbuffer.data(if_geom_10_off + 704 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xzz = cbuffer.data(if_geom_10_off + 705 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyy = cbuffer.data(if_geom_10_off + 706 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyz = cbuffer.data(if_geom_10_off + 707 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yzz = cbuffer.data(if_geom_10_off + 708 * ccomps * dcomps);

            auto g_z_0_xxzzzz_zzz = cbuffer.data(if_geom_10_off + 709 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzzzz_xxx, g_z_0_xxzzzz_xxy, g_z_0_xxzzzz_xxz, g_z_0_xxzzzz_xyy, g_z_0_xxzzzz_xyz, g_z_0_xxzzzz_xzz, g_z_0_xxzzzz_yyy, g_z_0_xxzzzz_yyz, g_z_0_xxzzzz_yzz, g_z_0_xxzzzz_zzz, g_z_0_xzzzz_xxx, g_z_0_xzzzz_xxxx, g_z_0_xzzzz_xxxy, g_z_0_xzzzz_xxxz, g_z_0_xzzzz_xxy, g_z_0_xzzzz_xxyy, g_z_0_xzzzz_xxyz, g_z_0_xzzzz_xxz, g_z_0_xzzzz_xxzz, g_z_0_xzzzz_xyy, g_z_0_xzzzz_xyyy, g_z_0_xzzzz_xyyz, g_z_0_xzzzz_xyz, g_z_0_xzzzz_xyzz, g_z_0_xzzzz_xzz, g_z_0_xzzzz_xzzz, g_z_0_xzzzz_yyy, g_z_0_xzzzz_yyz, g_z_0_xzzzz_yzz, g_z_0_xzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_xxx[k] = -g_z_0_xzzzz_xxx[k] * ab_x + g_z_0_xzzzz_xxxx[k];

                g_z_0_xxzzzz_xxy[k] = -g_z_0_xzzzz_xxy[k] * ab_x + g_z_0_xzzzz_xxxy[k];

                g_z_0_xxzzzz_xxz[k] = -g_z_0_xzzzz_xxz[k] * ab_x + g_z_0_xzzzz_xxxz[k];

                g_z_0_xxzzzz_xyy[k] = -g_z_0_xzzzz_xyy[k] * ab_x + g_z_0_xzzzz_xxyy[k];

                g_z_0_xxzzzz_xyz[k] = -g_z_0_xzzzz_xyz[k] * ab_x + g_z_0_xzzzz_xxyz[k];

                g_z_0_xxzzzz_xzz[k] = -g_z_0_xzzzz_xzz[k] * ab_x + g_z_0_xzzzz_xxzz[k];

                g_z_0_xxzzzz_yyy[k] = -g_z_0_xzzzz_yyy[k] * ab_x + g_z_0_xzzzz_xyyy[k];

                g_z_0_xxzzzz_yyz[k] = -g_z_0_xzzzz_yyz[k] * ab_x + g_z_0_xzzzz_xyyz[k];

                g_z_0_xxzzzz_yzz[k] = -g_z_0_xzzzz_yzz[k] * ab_x + g_z_0_xzzzz_xyzz[k];

                g_z_0_xxzzzz_zzz[k] = -g_z_0_xzzzz_zzz[k] * ab_x + g_z_0_xzzzz_xzzz[k];
            }

            /// Set up 710-720 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_xxx = cbuffer.data(if_geom_10_off + 710 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxy = cbuffer.data(if_geom_10_off + 711 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxz = cbuffer.data(if_geom_10_off + 712 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyy = cbuffer.data(if_geom_10_off + 713 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyz = cbuffer.data(if_geom_10_off + 714 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xzz = cbuffer.data(if_geom_10_off + 715 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyy = cbuffer.data(if_geom_10_off + 716 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyz = cbuffer.data(if_geom_10_off + 717 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yzz = cbuffer.data(if_geom_10_off + 718 * ccomps * dcomps);

            auto g_z_0_xyyyyy_zzz = cbuffer.data(if_geom_10_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyyy_xxx, g_z_0_xyyyyy_xxy, g_z_0_xyyyyy_xxz, g_z_0_xyyyyy_xyy, g_z_0_xyyyyy_xyz, g_z_0_xyyyyy_xzz, g_z_0_xyyyyy_yyy, g_z_0_xyyyyy_yyz, g_z_0_xyyyyy_yzz, g_z_0_xyyyyy_zzz, g_z_0_yyyyy_xxx, g_z_0_yyyyy_xxxx, g_z_0_yyyyy_xxxy, g_z_0_yyyyy_xxxz, g_z_0_yyyyy_xxy, g_z_0_yyyyy_xxyy, g_z_0_yyyyy_xxyz, g_z_0_yyyyy_xxz, g_z_0_yyyyy_xxzz, g_z_0_yyyyy_xyy, g_z_0_yyyyy_xyyy, g_z_0_yyyyy_xyyz, g_z_0_yyyyy_xyz, g_z_0_yyyyy_xyzz, g_z_0_yyyyy_xzz, g_z_0_yyyyy_xzzz, g_z_0_yyyyy_yyy, g_z_0_yyyyy_yyz, g_z_0_yyyyy_yzz, g_z_0_yyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_xxx[k] = -g_z_0_yyyyy_xxx[k] * ab_x + g_z_0_yyyyy_xxxx[k];

                g_z_0_xyyyyy_xxy[k] = -g_z_0_yyyyy_xxy[k] * ab_x + g_z_0_yyyyy_xxxy[k];

                g_z_0_xyyyyy_xxz[k] = -g_z_0_yyyyy_xxz[k] * ab_x + g_z_0_yyyyy_xxxz[k];

                g_z_0_xyyyyy_xyy[k] = -g_z_0_yyyyy_xyy[k] * ab_x + g_z_0_yyyyy_xxyy[k];

                g_z_0_xyyyyy_xyz[k] = -g_z_0_yyyyy_xyz[k] * ab_x + g_z_0_yyyyy_xxyz[k];

                g_z_0_xyyyyy_xzz[k] = -g_z_0_yyyyy_xzz[k] * ab_x + g_z_0_yyyyy_xxzz[k];

                g_z_0_xyyyyy_yyy[k] = -g_z_0_yyyyy_yyy[k] * ab_x + g_z_0_yyyyy_xyyy[k];

                g_z_0_xyyyyy_yyz[k] = -g_z_0_yyyyy_yyz[k] * ab_x + g_z_0_yyyyy_xyyz[k];

                g_z_0_xyyyyy_yzz[k] = -g_z_0_yyyyy_yzz[k] * ab_x + g_z_0_yyyyy_xyzz[k];

                g_z_0_xyyyyy_zzz[k] = -g_z_0_yyyyy_zzz[k] * ab_x + g_z_0_yyyyy_xzzz[k];
            }

            /// Set up 720-730 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_xxx = cbuffer.data(if_geom_10_off + 720 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxy = cbuffer.data(if_geom_10_off + 721 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxz = cbuffer.data(if_geom_10_off + 722 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyy = cbuffer.data(if_geom_10_off + 723 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyz = cbuffer.data(if_geom_10_off + 724 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xzz = cbuffer.data(if_geom_10_off + 725 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyy = cbuffer.data(if_geom_10_off + 726 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyz = cbuffer.data(if_geom_10_off + 727 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yzz = cbuffer.data(if_geom_10_off + 728 * ccomps * dcomps);

            auto g_z_0_xyyyyz_zzz = cbuffer.data(if_geom_10_off + 729 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyyz_xxx, g_z_0_xyyyyz_xxy, g_z_0_xyyyyz_xxz, g_z_0_xyyyyz_xyy, g_z_0_xyyyyz_xyz, g_z_0_xyyyyz_xzz, g_z_0_xyyyyz_yyy, g_z_0_xyyyyz_yyz, g_z_0_xyyyyz_yzz, g_z_0_xyyyyz_zzz, g_z_0_yyyyz_xxx, g_z_0_yyyyz_xxxx, g_z_0_yyyyz_xxxy, g_z_0_yyyyz_xxxz, g_z_0_yyyyz_xxy, g_z_0_yyyyz_xxyy, g_z_0_yyyyz_xxyz, g_z_0_yyyyz_xxz, g_z_0_yyyyz_xxzz, g_z_0_yyyyz_xyy, g_z_0_yyyyz_xyyy, g_z_0_yyyyz_xyyz, g_z_0_yyyyz_xyz, g_z_0_yyyyz_xyzz, g_z_0_yyyyz_xzz, g_z_0_yyyyz_xzzz, g_z_0_yyyyz_yyy, g_z_0_yyyyz_yyz, g_z_0_yyyyz_yzz, g_z_0_yyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_xxx[k] = -g_z_0_yyyyz_xxx[k] * ab_x + g_z_0_yyyyz_xxxx[k];

                g_z_0_xyyyyz_xxy[k] = -g_z_0_yyyyz_xxy[k] * ab_x + g_z_0_yyyyz_xxxy[k];

                g_z_0_xyyyyz_xxz[k] = -g_z_0_yyyyz_xxz[k] * ab_x + g_z_0_yyyyz_xxxz[k];

                g_z_0_xyyyyz_xyy[k] = -g_z_0_yyyyz_xyy[k] * ab_x + g_z_0_yyyyz_xxyy[k];

                g_z_0_xyyyyz_xyz[k] = -g_z_0_yyyyz_xyz[k] * ab_x + g_z_0_yyyyz_xxyz[k];

                g_z_0_xyyyyz_xzz[k] = -g_z_0_yyyyz_xzz[k] * ab_x + g_z_0_yyyyz_xxzz[k];

                g_z_0_xyyyyz_yyy[k] = -g_z_0_yyyyz_yyy[k] * ab_x + g_z_0_yyyyz_xyyy[k];

                g_z_0_xyyyyz_yyz[k] = -g_z_0_yyyyz_yyz[k] * ab_x + g_z_0_yyyyz_xyyz[k];

                g_z_0_xyyyyz_yzz[k] = -g_z_0_yyyyz_yzz[k] * ab_x + g_z_0_yyyyz_xyzz[k];

                g_z_0_xyyyyz_zzz[k] = -g_z_0_yyyyz_zzz[k] * ab_x + g_z_0_yyyyz_xzzz[k];
            }

            /// Set up 730-740 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_xxx = cbuffer.data(if_geom_10_off + 730 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxy = cbuffer.data(if_geom_10_off + 731 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxz = cbuffer.data(if_geom_10_off + 732 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyy = cbuffer.data(if_geom_10_off + 733 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyz = cbuffer.data(if_geom_10_off + 734 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xzz = cbuffer.data(if_geom_10_off + 735 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyy = cbuffer.data(if_geom_10_off + 736 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyz = cbuffer.data(if_geom_10_off + 737 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yzz = cbuffer.data(if_geom_10_off + 738 * ccomps * dcomps);

            auto g_z_0_xyyyzz_zzz = cbuffer.data(if_geom_10_off + 739 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyzz_xxx, g_z_0_xyyyzz_xxy, g_z_0_xyyyzz_xxz, g_z_0_xyyyzz_xyy, g_z_0_xyyyzz_xyz, g_z_0_xyyyzz_xzz, g_z_0_xyyyzz_yyy, g_z_0_xyyyzz_yyz, g_z_0_xyyyzz_yzz, g_z_0_xyyyzz_zzz, g_z_0_yyyzz_xxx, g_z_0_yyyzz_xxxx, g_z_0_yyyzz_xxxy, g_z_0_yyyzz_xxxz, g_z_0_yyyzz_xxy, g_z_0_yyyzz_xxyy, g_z_0_yyyzz_xxyz, g_z_0_yyyzz_xxz, g_z_0_yyyzz_xxzz, g_z_0_yyyzz_xyy, g_z_0_yyyzz_xyyy, g_z_0_yyyzz_xyyz, g_z_0_yyyzz_xyz, g_z_0_yyyzz_xyzz, g_z_0_yyyzz_xzz, g_z_0_yyyzz_xzzz, g_z_0_yyyzz_yyy, g_z_0_yyyzz_yyz, g_z_0_yyyzz_yzz, g_z_0_yyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_xxx[k] = -g_z_0_yyyzz_xxx[k] * ab_x + g_z_0_yyyzz_xxxx[k];

                g_z_0_xyyyzz_xxy[k] = -g_z_0_yyyzz_xxy[k] * ab_x + g_z_0_yyyzz_xxxy[k];

                g_z_0_xyyyzz_xxz[k] = -g_z_0_yyyzz_xxz[k] * ab_x + g_z_0_yyyzz_xxxz[k];

                g_z_0_xyyyzz_xyy[k] = -g_z_0_yyyzz_xyy[k] * ab_x + g_z_0_yyyzz_xxyy[k];

                g_z_0_xyyyzz_xyz[k] = -g_z_0_yyyzz_xyz[k] * ab_x + g_z_0_yyyzz_xxyz[k];

                g_z_0_xyyyzz_xzz[k] = -g_z_0_yyyzz_xzz[k] * ab_x + g_z_0_yyyzz_xxzz[k];

                g_z_0_xyyyzz_yyy[k] = -g_z_0_yyyzz_yyy[k] * ab_x + g_z_0_yyyzz_xyyy[k];

                g_z_0_xyyyzz_yyz[k] = -g_z_0_yyyzz_yyz[k] * ab_x + g_z_0_yyyzz_xyyz[k];

                g_z_0_xyyyzz_yzz[k] = -g_z_0_yyyzz_yzz[k] * ab_x + g_z_0_yyyzz_xyzz[k];

                g_z_0_xyyyzz_zzz[k] = -g_z_0_yyyzz_zzz[k] * ab_x + g_z_0_yyyzz_xzzz[k];
            }

            /// Set up 740-750 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_xxx = cbuffer.data(if_geom_10_off + 740 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxy = cbuffer.data(if_geom_10_off + 741 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxz = cbuffer.data(if_geom_10_off + 742 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyy = cbuffer.data(if_geom_10_off + 743 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyz = cbuffer.data(if_geom_10_off + 744 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xzz = cbuffer.data(if_geom_10_off + 745 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyy = cbuffer.data(if_geom_10_off + 746 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyz = cbuffer.data(if_geom_10_off + 747 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yzz = cbuffer.data(if_geom_10_off + 748 * ccomps * dcomps);

            auto g_z_0_xyyzzz_zzz = cbuffer.data(if_geom_10_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyzzz_xxx, g_z_0_xyyzzz_xxy, g_z_0_xyyzzz_xxz, g_z_0_xyyzzz_xyy, g_z_0_xyyzzz_xyz, g_z_0_xyyzzz_xzz, g_z_0_xyyzzz_yyy, g_z_0_xyyzzz_yyz, g_z_0_xyyzzz_yzz, g_z_0_xyyzzz_zzz, g_z_0_yyzzz_xxx, g_z_0_yyzzz_xxxx, g_z_0_yyzzz_xxxy, g_z_0_yyzzz_xxxz, g_z_0_yyzzz_xxy, g_z_0_yyzzz_xxyy, g_z_0_yyzzz_xxyz, g_z_0_yyzzz_xxz, g_z_0_yyzzz_xxzz, g_z_0_yyzzz_xyy, g_z_0_yyzzz_xyyy, g_z_0_yyzzz_xyyz, g_z_0_yyzzz_xyz, g_z_0_yyzzz_xyzz, g_z_0_yyzzz_xzz, g_z_0_yyzzz_xzzz, g_z_0_yyzzz_yyy, g_z_0_yyzzz_yyz, g_z_0_yyzzz_yzz, g_z_0_yyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_xxx[k] = -g_z_0_yyzzz_xxx[k] * ab_x + g_z_0_yyzzz_xxxx[k];

                g_z_0_xyyzzz_xxy[k] = -g_z_0_yyzzz_xxy[k] * ab_x + g_z_0_yyzzz_xxxy[k];

                g_z_0_xyyzzz_xxz[k] = -g_z_0_yyzzz_xxz[k] * ab_x + g_z_0_yyzzz_xxxz[k];

                g_z_0_xyyzzz_xyy[k] = -g_z_0_yyzzz_xyy[k] * ab_x + g_z_0_yyzzz_xxyy[k];

                g_z_0_xyyzzz_xyz[k] = -g_z_0_yyzzz_xyz[k] * ab_x + g_z_0_yyzzz_xxyz[k];

                g_z_0_xyyzzz_xzz[k] = -g_z_0_yyzzz_xzz[k] * ab_x + g_z_0_yyzzz_xxzz[k];

                g_z_0_xyyzzz_yyy[k] = -g_z_0_yyzzz_yyy[k] * ab_x + g_z_0_yyzzz_xyyy[k];

                g_z_0_xyyzzz_yyz[k] = -g_z_0_yyzzz_yyz[k] * ab_x + g_z_0_yyzzz_xyyz[k];

                g_z_0_xyyzzz_yzz[k] = -g_z_0_yyzzz_yzz[k] * ab_x + g_z_0_yyzzz_xyzz[k];

                g_z_0_xyyzzz_zzz[k] = -g_z_0_yyzzz_zzz[k] * ab_x + g_z_0_yyzzz_xzzz[k];
            }

            /// Set up 750-760 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_xxx = cbuffer.data(if_geom_10_off + 750 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxy = cbuffer.data(if_geom_10_off + 751 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxz = cbuffer.data(if_geom_10_off + 752 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyy = cbuffer.data(if_geom_10_off + 753 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyz = cbuffer.data(if_geom_10_off + 754 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xzz = cbuffer.data(if_geom_10_off + 755 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyy = cbuffer.data(if_geom_10_off + 756 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyz = cbuffer.data(if_geom_10_off + 757 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yzz = cbuffer.data(if_geom_10_off + 758 * ccomps * dcomps);

            auto g_z_0_xyzzzz_zzz = cbuffer.data(if_geom_10_off + 759 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzzzz_xxx, g_z_0_xyzzzz_xxy, g_z_0_xyzzzz_xxz, g_z_0_xyzzzz_xyy, g_z_0_xyzzzz_xyz, g_z_0_xyzzzz_xzz, g_z_0_xyzzzz_yyy, g_z_0_xyzzzz_yyz, g_z_0_xyzzzz_yzz, g_z_0_xyzzzz_zzz, g_z_0_yzzzz_xxx, g_z_0_yzzzz_xxxx, g_z_0_yzzzz_xxxy, g_z_0_yzzzz_xxxz, g_z_0_yzzzz_xxy, g_z_0_yzzzz_xxyy, g_z_0_yzzzz_xxyz, g_z_0_yzzzz_xxz, g_z_0_yzzzz_xxzz, g_z_0_yzzzz_xyy, g_z_0_yzzzz_xyyy, g_z_0_yzzzz_xyyz, g_z_0_yzzzz_xyz, g_z_0_yzzzz_xyzz, g_z_0_yzzzz_xzz, g_z_0_yzzzz_xzzz, g_z_0_yzzzz_yyy, g_z_0_yzzzz_yyz, g_z_0_yzzzz_yzz, g_z_0_yzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_xxx[k] = -g_z_0_yzzzz_xxx[k] * ab_x + g_z_0_yzzzz_xxxx[k];

                g_z_0_xyzzzz_xxy[k] = -g_z_0_yzzzz_xxy[k] * ab_x + g_z_0_yzzzz_xxxy[k];

                g_z_0_xyzzzz_xxz[k] = -g_z_0_yzzzz_xxz[k] * ab_x + g_z_0_yzzzz_xxxz[k];

                g_z_0_xyzzzz_xyy[k] = -g_z_0_yzzzz_xyy[k] * ab_x + g_z_0_yzzzz_xxyy[k];

                g_z_0_xyzzzz_xyz[k] = -g_z_0_yzzzz_xyz[k] * ab_x + g_z_0_yzzzz_xxyz[k];

                g_z_0_xyzzzz_xzz[k] = -g_z_0_yzzzz_xzz[k] * ab_x + g_z_0_yzzzz_xxzz[k];

                g_z_0_xyzzzz_yyy[k] = -g_z_0_yzzzz_yyy[k] * ab_x + g_z_0_yzzzz_xyyy[k];

                g_z_0_xyzzzz_yyz[k] = -g_z_0_yzzzz_yyz[k] * ab_x + g_z_0_yzzzz_xyyz[k];

                g_z_0_xyzzzz_yzz[k] = -g_z_0_yzzzz_yzz[k] * ab_x + g_z_0_yzzzz_xyzz[k];

                g_z_0_xyzzzz_zzz[k] = -g_z_0_yzzzz_zzz[k] * ab_x + g_z_0_yzzzz_xzzz[k];
            }

            /// Set up 760-770 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_xxx = cbuffer.data(if_geom_10_off + 760 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxy = cbuffer.data(if_geom_10_off + 761 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxz = cbuffer.data(if_geom_10_off + 762 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyy = cbuffer.data(if_geom_10_off + 763 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyz = cbuffer.data(if_geom_10_off + 764 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xzz = cbuffer.data(if_geom_10_off + 765 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyy = cbuffer.data(if_geom_10_off + 766 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyz = cbuffer.data(if_geom_10_off + 767 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yzz = cbuffer.data(if_geom_10_off + 768 * ccomps * dcomps);

            auto g_z_0_xzzzzz_zzz = cbuffer.data(if_geom_10_off + 769 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzzzz_xxx, g_z_0_xzzzzz_xxy, g_z_0_xzzzzz_xxz, g_z_0_xzzzzz_xyy, g_z_0_xzzzzz_xyz, g_z_0_xzzzzz_xzz, g_z_0_xzzzzz_yyy, g_z_0_xzzzzz_yyz, g_z_0_xzzzzz_yzz, g_z_0_xzzzzz_zzz, g_z_0_zzzzz_xxx, g_z_0_zzzzz_xxxx, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxy, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xyy, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_yyy, g_z_0_zzzzz_yyz, g_z_0_zzzzz_yzz, g_z_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_xxx[k] = -g_z_0_zzzzz_xxx[k] * ab_x + g_z_0_zzzzz_xxxx[k];

                g_z_0_xzzzzz_xxy[k] = -g_z_0_zzzzz_xxy[k] * ab_x + g_z_0_zzzzz_xxxy[k];

                g_z_0_xzzzzz_xxz[k] = -g_z_0_zzzzz_xxz[k] * ab_x + g_z_0_zzzzz_xxxz[k];

                g_z_0_xzzzzz_xyy[k] = -g_z_0_zzzzz_xyy[k] * ab_x + g_z_0_zzzzz_xxyy[k];

                g_z_0_xzzzzz_xyz[k] = -g_z_0_zzzzz_xyz[k] * ab_x + g_z_0_zzzzz_xxyz[k];

                g_z_0_xzzzzz_xzz[k] = -g_z_0_zzzzz_xzz[k] * ab_x + g_z_0_zzzzz_xxzz[k];

                g_z_0_xzzzzz_yyy[k] = -g_z_0_zzzzz_yyy[k] * ab_x + g_z_0_zzzzz_xyyy[k];

                g_z_0_xzzzzz_yyz[k] = -g_z_0_zzzzz_yyz[k] * ab_x + g_z_0_zzzzz_xyyz[k];

                g_z_0_xzzzzz_yzz[k] = -g_z_0_zzzzz_yzz[k] * ab_x + g_z_0_zzzzz_xyzz[k];

                g_z_0_xzzzzz_zzz[k] = -g_z_0_zzzzz_zzz[k] * ab_x + g_z_0_zzzzz_xzzz[k];
            }

            /// Set up 770-780 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_xxx = cbuffer.data(if_geom_10_off + 770 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxy = cbuffer.data(if_geom_10_off + 771 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxz = cbuffer.data(if_geom_10_off + 772 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyy = cbuffer.data(if_geom_10_off + 773 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyz = cbuffer.data(if_geom_10_off + 774 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xzz = cbuffer.data(if_geom_10_off + 775 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyy = cbuffer.data(if_geom_10_off + 776 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyz = cbuffer.data(if_geom_10_off + 777 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yzz = cbuffer.data(if_geom_10_off + 778 * ccomps * dcomps);

            auto g_z_0_yyyyyy_zzz = cbuffer.data(if_geom_10_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyy_xxx, g_z_0_yyyyy_xxxy, g_z_0_yyyyy_xxy, g_z_0_yyyyy_xxyy, g_z_0_yyyyy_xxyz, g_z_0_yyyyy_xxz, g_z_0_yyyyy_xyy, g_z_0_yyyyy_xyyy, g_z_0_yyyyy_xyyz, g_z_0_yyyyy_xyz, g_z_0_yyyyy_xyzz, g_z_0_yyyyy_xzz, g_z_0_yyyyy_yyy, g_z_0_yyyyy_yyyy, g_z_0_yyyyy_yyyz, g_z_0_yyyyy_yyz, g_z_0_yyyyy_yyzz, g_z_0_yyyyy_yzz, g_z_0_yyyyy_yzzz, g_z_0_yyyyy_zzz, g_z_0_yyyyyy_xxx, g_z_0_yyyyyy_xxy, g_z_0_yyyyyy_xxz, g_z_0_yyyyyy_xyy, g_z_0_yyyyyy_xyz, g_z_0_yyyyyy_xzz, g_z_0_yyyyyy_yyy, g_z_0_yyyyyy_yyz, g_z_0_yyyyyy_yzz, g_z_0_yyyyyy_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_xxx[k] = -g_z_0_yyyyy_xxx[k] * ab_y + g_z_0_yyyyy_xxxy[k];

                g_z_0_yyyyyy_xxy[k] = -g_z_0_yyyyy_xxy[k] * ab_y + g_z_0_yyyyy_xxyy[k];

                g_z_0_yyyyyy_xxz[k] = -g_z_0_yyyyy_xxz[k] * ab_y + g_z_0_yyyyy_xxyz[k];

                g_z_0_yyyyyy_xyy[k] = -g_z_0_yyyyy_xyy[k] * ab_y + g_z_0_yyyyy_xyyy[k];

                g_z_0_yyyyyy_xyz[k] = -g_z_0_yyyyy_xyz[k] * ab_y + g_z_0_yyyyy_xyyz[k];

                g_z_0_yyyyyy_xzz[k] = -g_z_0_yyyyy_xzz[k] * ab_y + g_z_0_yyyyy_xyzz[k];

                g_z_0_yyyyyy_yyy[k] = -g_z_0_yyyyy_yyy[k] * ab_y + g_z_0_yyyyy_yyyy[k];

                g_z_0_yyyyyy_yyz[k] = -g_z_0_yyyyy_yyz[k] * ab_y + g_z_0_yyyyy_yyyz[k];

                g_z_0_yyyyyy_yzz[k] = -g_z_0_yyyyy_yzz[k] * ab_y + g_z_0_yyyyy_yyzz[k];

                g_z_0_yyyyyy_zzz[k] = -g_z_0_yyyyy_zzz[k] * ab_y + g_z_0_yyyyy_yzzz[k];
            }

            /// Set up 780-790 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_xxx = cbuffer.data(if_geom_10_off + 780 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxy = cbuffer.data(if_geom_10_off + 781 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxz = cbuffer.data(if_geom_10_off + 782 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyy = cbuffer.data(if_geom_10_off + 783 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyz = cbuffer.data(if_geom_10_off + 784 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xzz = cbuffer.data(if_geom_10_off + 785 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyy = cbuffer.data(if_geom_10_off + 786 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyz = cbuffer.data(if_geom_10_off + 787 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yzz = cbuffer.data(if_geom_10_off + 788 * ccomps * dcomps);

            auto g_z_0_yyyyyz_zzz = cbuffer.data(if_geom_10_off + 789 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyyz_xxx, g_z_0_yyyyyz_xxy, g_z_0_yyyyyz_xxz, g_z_0_yyyyyz_xyy, g_z_0_yyyyyz_xyz, g_z_0_yyyyyz_xzz, g_z_0_yyyyyz_yyy, g_z_0_yyyyyz_yyz, g_z_0_yyyyyz_yzz, g_z_0_yyyyyz_zzz, g_z_0_yyyyz_xxx, g_z_0_yyyyz_xxxy, g_z_0_yyyyz_xxy, g_z_0_yyyyz_xxyy, g_z_0_yyyyz_xxyz, g_z_0_yyyyz_xxz, g_z_0_yyyyz_xyy, g_z_0_yyyyz_xyyy, g_z_0_yyyyz_xyyz, g_z_0_yyyyz_xyz, g_z_0_yyyyz_xyzz, g_z_0_yyyyz_xzz, g_z_0_yyyyz_yyy, g_z_0_yyyyz_yyyy, g_z_0_yyyyz_yyyz, g_z_0_yyyyz_yyz, g_z_0_yyyyz_yyzz, g_z_0_yyyyz_yzz, g_z_0_yyyyz_yzzz, g_z_0_yyyyz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_xxx[k] = -g_z_0_yyyyz_xxx[k] * ab_y + g_z_0_yyyyz_xxxy[k];

                g_z_0_yyyyyz_xxy[k] = -g_z_0_yyyyz_xxy[k] * ab_y + g_z_0_yyyyz_xxyy[k];

                g_z_0_yyyyyz_xxz[k] = -g_z_0_yyyyz_xxz[k] * ab_y + g_z_0_yyyyz_xxyz[k];

                g_z_0_yyyyyz_xyy[k] = -g_z_0_yyyyz_xyy[k] * ab_y + g_z_0_yyyyz_xyyy[k];

                g_z_0_yyyyyz_xyz[k] = -g_z_0_yyyyz_xyz[k] * ab_y + g_z_0_yyyyz_xyyz[k];

                g_z_0_yyyyyz_xzz[k] = -g_z_0_yyyyz_xzz[k] * ab_y + g_z_0_yyyyz_xyzz[k];

                g_z_0_yyyyyz_yyy[k] = -g_z_0_yyyyz_yyy[k] * ab_y + g_z_0_yyyyz_yyyy[k];

                g_z_0_yyyyyz_yyz[k] = -g_z_0_yyyyz_yyz[k] * ab_y + g_z_0_yyyyz_yyyz[k];

                g_z_0_yyyyyz_yzz[k] = -g_z_0_yyyyz_yzz[k] * ab_y + g_z_0_yyyyz_yyzz[k];

                g_z_0_yyyyyz_zzz[k] = -g_z_0_yyyyz_zzz[k] * ab_y + g_z_0_yyyyz_yzzz[k];
            }

            /// Set up 790-800 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_xxx = cbuffer.data(if_geom_10_off + 790 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxy = cbuffer.data(if_geom_10_off + 791 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxz = cbuffer.data(if_geom_10_off + 792 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyy = cbuffer.data(if_geom_10_off + 793 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyz = cbuffer.data(if_geom_10_off + 794 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xzz = cbuffer.data(if_geom_10_off + 795 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyy = cbuffer.data(if_geom_10_off + 796 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyz = cbuffer.data(if_geom_10_off + 797 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yzz = cbuffer.data(if_geom_10_off + 798 * ccomps * dcomps);

            auto g_z_0_yyyyzz_zzz = cbuffer.data(if_geom_10_off + 799 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyzz_xxx, g_z_0_yyyyzz_xxy, g_z_0_yyyyzz_xxz, g_z_0_yyyyzz_xyy, g_z_0_yyyyzz_xyz, g_z_0_yyyyzz_xzz, g_z_0_yyyyzz_yyy, g_z_0_yyyyzz_yyz, g_z_0_yyyyzz_yzz, g_z_0_yyyyzz_zzz, g_z_0_yyyzz_xxx, g_z_0_yyyzz_xxxy, g_z_0_yyyzz_xxy, g_z_0_yyyzz_xxyy, g_z_0_yyyzz_xxyz, g_z_0_yyyzz_xxz, g_z_0_yyyzz_xyy, g_z_0_yyyzz_xyyy, g_z_0_yyyzz_xyyz, g_z_0_yyyzz_xyz, g_z_0_yyyzz_xyzz, g_z_0_yyyzz_xzz, g_z_0_yyyzz_yyy, g_z_0_yyyzz_yyyy, g_z_0_yyyzz_yyyz, g_z_0_yyyzz_yyz, g_z_0_yyyzz_yyzz, g_z_0_yyyzz_yzz, g_z_0_yyyzz_yzzz, g_z_0_yyyzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_xxx[k] = -g_z_0_yyyzz_xxx[k] * ab_y + g_z_0_yyyzz_xxxy[k];

                g_z_0_yyyyzz_xxy[k] = -g_z_0_yyyzz_xxy[k] * ab_y + g_z_0_yyyzz_xxyy[k];

                g_z_0_yyyyzz_xxz[k] = -g_z_0_yyyzz_xxz[k] * ab_y + g_z_0_yyyzz_xxyz[k];

                g_z_0_yyyyzz_xyy[k] = -g_z_0_yyyzz_xyy[k] * ab_y + g_z_0_yyyzz_xyyy[k];

                g_z_0_yyyyzz_xyz[k] = -g_z_0_yyyzz_xyz[k] * ab_y + g_z_0_yyyzz_xyyz[k];

                g_z_0_yyyyzz_xzz[k] = -g_z_0_yyyzz_xzz[k] * ab_y + g_z_0_yyyzz_xyzz[k];

                g_z_0_yyyyzz_yyy[k] = -g_z_0_yyyzz_yyy[k] * ab_y + g_z_0_yyyzz_yyyy[k];

                g_z_0_yyyyzz_yyz[k] = -g_z_0_yyyzz_yyz[k] * ab_y + g_z_0_yyyzz_yyyz[k];

                g_z_0_yyyyzz_yzz[k] = -g_z_0_yyyzz_yzz[k] * ab_y + g_z_0_yyyzz_yyzz[k];

                g_z_0_yyyyzz_zzz[k] = -g_z_0_yyyzz_zzz[k] * ab_y + g_z_0_yyyzz_yzzz[k];
            }

            /// Set up 800-810 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_xxx = cbuffer.data(if_geom_10_off + 800 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxy = cbuffer.data(if_geom_10_off + 801 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxz = cbuffer.data(if_geom_10_off + 802 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyy = cbuffer.data(if_geom_10_off + 803 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyz = cbuffer.data(if_geom_10_off + 804 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xzz = cbuffer.data(if_geom_10_off + 805 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyy = cbuffer.data(if_geom_10_off + 806 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyz = cbuffer.data(if_geom_10_off + 807 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yzz = cbuffer.data(if_geom_10_off + 808 * ccomps * dcomps);

            auto g_z_0_yyyzzz_zzz = cbuffer.data(if_geom_10_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyzzz_xxx, g_z_0_yyyzzz_xxy, g_z_0_yyyzzz_xxz, g_z_0_yyyzzz_xyy, g_z_0_yyyzzz_xyz, g_z_0_yyyzzz_xzz, g_z_0_yyyzzz_yyy, g_z_0_yyyzzz_yyz, g_z_0_yyyzzz_yzz, g_z_0_yyyzzz_zzz, g_z_0_yyzzz_xxx, g_z_0_yyzzz_xxxy, g_z_0_yyzzz_xxy, g_z_0_yyzzz_xxyy, g_z_0_yyzzz_xxyz, g_z_0_yyzzz_xxz, g_z_0_yyzzz_xyy, g_z_0_yyzzz_xyyy, g_z_0_yyzzz_xyyz, g_z_0_yyzzz_xyz, g_z_0_yyzzz_xyzz, g_z_0_yyzzz_xzz, g_z_0_yyzzz_yyy, g_z_0_yyzzz_yyyy, g_z_0_yyzzz_yyyz, g_z_0_yyzzz_yyz, g_z_0_yyzzz_yyzz, g_z_0_yyzzz_yzz, g_z_0_yyzzz_yzzz, g_z_0_yyzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_xxx[k] = -g_z_0_yyzzz_xxx[k] * ab_y + g_z_0_yyzzz_xxxy[k];

                g_z_0_yyyzzz_xxy[k] = -g_z_0_yyzzz_xxy[k] * ab_y + g_z_0_yyzzz_xxyy[k];

                g_z_0_yyyzzz_xxz[k] = -g_z_0_yyzzz_xxz[k] * ab_y + g_z_0_yyzzz_xxyz[k];

                g_z_0_yyyzzz_xyy[k] = -g_z_0_yyzzz_xyy[k] * ab_y + g_z_0_yyzzz_xyyy[k];

                g_z_0_yyyzzz_xyz[k] = -g_z_0_yyzzz_xyz[k] * ab_y + g_z_0_yyzzz_xyyz[k];

                g_z_0_yyyzzz_xzz[k] = -g_z_0_yyzzz_xzz[k] * ab_y + g_z_0_yyzzz_xyzz[k];

                g_z_0_yyyzzz_yyy[k] = -g_z_0_yyzzz_yyy[k] * ab_y + g_z_0_yyzzz_yyyy[k];

                g_z_0_yyyzzz_yyz[k] = -g_z_0_yyzzz_yyz[k] * ab_y + g_z_0_yyzzz_yyyz[k];

                g_z_0_yyyzzz_yzz[k] = -g_z_0_yyzzz_yzz[k] * ab_y + g_z_0_yyzzz_yyzz[k];

                g_z_0_yyyzzz_zzz[k] = -g_z_0_yyzzz_zzz[k] * ab_y + g_z_0_yyzzz_yzzz[k];
            }

            /// Set up 810-820 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_xxx = cbuffer.data(if_geom_10_off + 810 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxy = cbuffer.data(if_geom_10_off + 811 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxz = cbuffer.data(if_geom_10_off + 812 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyy = cbuffer.data(if_geom_10_off + 813 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyz = cbuffer.data(if_geom_10_off + 814 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xzz = cbuffer.data(if_geom_10_off + 815 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyy = cbuffer.data(if_geom_10_off + 816 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyz = cbuffer.data(if_geom_10_off + 817 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yzz = cbuffer.data(if_geom_10_off + 818 * ccomps * dcomps);

            auto g_z_0_yyzzzz_zzz = cbuffer.data(if_geom_10_off + 819 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzzzz_xxx, g_z_0_yyzzzz_xxy, g_z_0_yyzzzz_xxz, g_z_0_yyzzzz_xyy, g_z_0_yyzzzz_xyz, g_z_0_yyzzzz_xzz, g_z_0_yyzzzz_yyy, g_z_0_yyzzzz_yyz, g_z_0_yyzzzz_yzz, g_z_0_yyzzzz_zzz, g_z_0_yzzzz_xxx, g_z_0_yzzzz_xxxy, g_z_0_yzzzz_xxy, g_z_0_yzzzz_xxyy, g_z_0_yzzzz_xxyz, g_z_0_yzzzz_xxz, g_z_0_yzzzz_xyy, g_z_0_yzzzz_xyyy, g_z_0_yzzzz_xyyz, g_z_0_yzzzz_xyz, g_z_0_yzzzz_xyzz, g_z_0_yzzzz_xzz, g_z_0_yzzzz_yyy, g_z_0_yzzzz_yyyy, g_z_0_yzzzz_yyyz, g_z_0_yzzzz_yyz, g_z_0_yzzzz_yyzz, g_z_0_yzzzz_yzz, g_z_0_yzzzz_yzzz, g_z_0_yzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_xxx[k] = -g_z_0_yzzzz_xxx[k] * ab_y + g_z_0_yzzzz_xxxy[k];

                g_z_0_yyzzzz_xxy[k] = -g_z_0_yzzzz_xxy[k] * ab_y + g_z_0_yzzzz_xxyy[k];

                g_z_0_yyzzzz_xxz[k] = -g_z_0_yzzzz_xxz[k] * ab_y + g_z_0_yzzzz_xxyz[k];

                g_z_0_yyzzzz_xyy[k] = -g_z_0_yzzzz_xyy[k] * ab_y + g_z_0_yzzzz_xyyy[k];

                g_z_0_yyzzzz_xyz[k] = -g_z_0_yzzzz_xyz[k] * ab_y + g_z_0_yzzzz_xyyz[k];

                g_z_0_yyzzzz_xzz[k] = -g_z_0_yzzzz_xzz[k] * ab_y + g_z_0_yzzzz_xyzz[k];

                g_z_0_yyzzzz_yyy[k] = -g_z_0_yzzzz_yyy[k] * ab_y + g_z_0_yzzzz_yyyy[k];

                g_z_0_yyzzzz_yyz[k] = -g_z_0_yzzzz_yyz[k] * ab_y + g_z_0_yzzzz_yyyz[k];

                g_z_0_yyzzzz_yzz[k] = -g_z_0_yzzzz_yzz[k] * ab_y + g_z_0_yzzzz_yyzz[k];

                g_z_0_yyzzzz_zzz[k] = -g_z_0_yzzzz_zzz[k] * ab_y + g_z_0_yzzzz_yzzz[k];
            }

            /// Set up 820-830 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_xxx = cbuffer.data(if_geom_10_off + 820 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxy = cbuffer.data(if_geom_10_off + 821 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxz = cbuffer.data(if_geom_10_off + 822 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyy = cbuffer.data(if_geom_10_off + 823 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyz = cbuffer.data(if_geom_10_off + 824 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xzz = cbuffer.data(if_geom_10_off + 825 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyy = cbuffer.data(if_geom_10_off + 826 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyz = cbuffer.data(if_geom_10_off + 827 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yzz = cbuffer.data(if_geom_10_off + 828 * ccomps * dcomps);

            auto g_z_0_yzzzzz_zzz = cbuffer.data(if_geom_10_off + 829 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzzzz_xxx, g_z_0_yzzzzz_xxy, g_z_0_yzzzzz_xxz, g_z_0_yzzzzz_xyy, g_z_0_yzzzzz_xyz, g_z_0_yzzzzz_xzz, g_z_0_yzzzzz_yyy, g_z_0_yzzzzz_yyz, g_z_0_yzzzzz_yzz, g_z_0_yzzzzz_zzz, g_z_0_zzzzz_xxx, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxy, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxz, g_z_0_zzzzz_xyy, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xzz, g_z_0_zzzzz_yyy, g_z_0_zzzzz_yyyy, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_xxx[k] = -g_z_0_zzzzz_xxx[k] * ab_y + g_z_0_zzzzz_xxxy[k];

                g_z_0_yzzzzz_xxy[k] = -g_z_0_zzzzz_xxy[k] * ab_y + g_z_0_zzzzz_xxyy[k];

                g_z_0_yzzzzz_xxz[k] = -g_z_0_zzzzz_xxz[k] * ab_y + g_z_0_zzzzz_xxyz[k];

                g_z_0_yzzzzz_xyy[k] = -g_z_0_zzzzz_xyy[k] * ab_y + g_z_0_zzzzz_xyyy[k];

                g_z_0_yzzzzz_xyz[k] = -g_z_0_zzzzz_xyz[k] * ab_y + g_z_0_zzzzz_xyyz[k];

                g_z_0_yzzzzz_xzz[k] = -g_z_0_zzzzz_xzz[k] * ab_y + g_z_0_zzzzz_xyzz[k];

                g_z_0_yzzzzz_yyy[k] = -g_z_0_zzzzz_yyy[k] * ab_y + g_z_0_zzzzz_yyyy[k];

                g_z_0_yzzzzz_yyz[k] = -g_z_0_zzzzz_yyz[k] * ab_y + g_z_0_zzzzz_yyyz[k];

                g_z_0_yzzzzz_yzz[k] = -g_z_0_zzzzz_yzz[k] * ab_y + g_z_0_zzzzz_yyzz[k];

                g_z_0_yzzzzz_zzz[k] = -g_z_0_zzzzz_zzz[k] * ab_y + g_z_0_zzzzz_yzzz[k];
            }

            /// Set up 830-840 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_xxx = cbuffer.data(if_geom_10_off + 830 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxy = cbuffer.data(if_geom_10_off + 831 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxz = cbuffer.data(if_geom_10_off + 832 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyy = cbuffer.data(if_geom_10_off + 833 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyz = cbuffer.data(if_geom_10_off + 834 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xzz = cbuffer.data(if_geom_10_off + 835 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyy = cbuffer.data(if_geom_10_off + 836 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyz = cbuffer.data(if_geom_10_off + 837 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yzz = cbuffer.data(if_geom_10_off + 838 * ccomps * dcomps);

            auto g_z_0_zzzzzz_zzz = cbuffer.data(if_geom_10_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzzz_xxx, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxy, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xyy, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_yyy, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_zzz, g_z_0_zzzzz_zzzz, g_z_0_zzzzzz_xxx, g_z_0_zzzzzz_xxy, g_z_0_zzzzzz_xxz, g_z_0_zzzzzz_xyy, g_z_0_zzzzzz_xyz, g_z_0_zzzzzz_xzz, g_z_0_zzzzzz_yyy, g_z_0_zzzzzz_yyz, g_z_0_zzzzzz_yzz, g_z_0_zzzzzz_zzz, g_zzzzz_xxx, g_zzzzz_xxy, g_zzzzz_xxz, g_zzzzz_xyy, g_zzzzz_xyz, g_zzzzz_xzz, g_zzzzz_yyy, g_zzzzz_yyz, g_zzzzz_yzz, g_zzzzz_zzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_xxx[k] = -g_zzzzz_xxx[k] - g_z_0_zzzzz_xxx[k] * ab_z + g_z_0_zzzzz_xxxz[k];

                g_z_0_zzzzzz_xxy[k] = -g_zzzzz_xxy[k] - g_z_0_zzzzz_xxy[k] * ab_z + g_z_0_zzzzz_xxyz[k];

                g_z_0_zzzzzz_xxz[k] = -g_zzzzz_xxz[k] - g_z_0_zzzzz_xxz[k] * ab_z + g_z_0_zzzzz_xxzz[k];

                g_z_0_zzzzzz_xyy[k] = -g_zzzzz_xyy[k] - g_z_0_zzzzz_xyy[k] * ab_z + g_z_0_zzzzz_xyyz[k];

                g_z_0_zzzzzz_xyz[k] = -g_zzzzz_xyz[k] - g_z_0_zzzzz_xyz[k] * ab_z + g_z_0_zzzzz_xyzz[k];

                g_z_0_zzzzzz_xzz[k] = -g_zzzzz_xzz[k] - g_z_0_zzzzz_xzz[k] * ab_z + g_z_0_zzzzz_xzzz[k];

                g_z_0_zzzzzz_yyy[k] = -g_zzzzz_yyy[k] - g_z_0_zzzzz_yyy[k] * ab_z + g_z_0_zzzzz_yyyz[k];

                g_z_0_zzzzzz_yyz[k] = -g_zzzzz_yyz[k] - g_z_0_zzzzz_yyz[k] * ab_z + g_z_0_zzzzz_yyzz[k];

                g_z_0_zzzzzz_yzz[k] = -g_zzzzz_yzz[k] - g_z_0_zzzzz_yzz[k] * ab_z + g_z_0_zzzzz_yzzz[k];

                g_z_0_zzzzzz_zzz[k] = -g_zzzzz_zzz[k] - g_z_0_zzzzz_zzz[k] * ab_z + g_z_0_zzzzz_zzzz[k];
            }
        }
    }
}

} // erirec namespace

