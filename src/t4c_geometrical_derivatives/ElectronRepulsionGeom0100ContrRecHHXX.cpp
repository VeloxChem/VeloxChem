#include "ElectronRepulsionGeom0100ContrRecHHXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_hhxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_hhxx,
                                            const size_t idx_ghxx,
                                            const size_t idx_geom_01_ghxx,
                                            const size_t idx_geom_01_gixx,
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
            /// Set up components of auxilary buffer : GHSS

            const auto gh_off = idx_ghxx + i * dcomps + j;

            auto g_xxxx_xxxxx = cbuffer.data(gh_off + 0 * ccomps * dcomps);

            auto g_xxxx_xxxxy = cbuffer.data(gh_off + 1 * ccomps * dcomps);

            auto g_xxxx_xxxxz = cbuffer.data(gh_off + 2 * ccomps * dcomps);

            auto g_xxxx_xxxyy = cbuffer.data(gh_off + 3 * ccomps * dcomps);

            auto g_xxxx_xxxyz = cbuffer.data(gh_off + 4 * ccomps * dcomps);

            auto g_xxxx_xxxzz = cbuffer.data(gh_off + 5 * ccomps * dcomps);

            auto g_xxxx_xxyyy = cbuffer.data(gh_off + 6 * ccomps * dcomps);

            auto g_xxxx_xxyyz = cbuffer.data(gh_off + 7 * ccomps * dcomps);

            auto g_xxxx_xxyzz = cbuffer.data(gh_off + 8 * ccomps * dcomps);

            auto g_xxxx_xxzzz = cbuffer.data(gh_off + 9 * ccomps * dcomps);

            auto g_xxxx_xyyyy = cbuffer.data(gh_off + 10 * ccomps * dcomps);

            auto g_xxxx_xyyyz = cbuffer.data(gh_off + 11 * ccomps * dcomps);

            auto g_xxxx_xyyzz = cbuffer.data(gh_off + 12 * ccomps * dcomps);

            auto g_xxxx_xyzzz = cbuffer.data(gh_off + 13 * ccomps * dcomps);

            auto g_xxxx_xzzzz = cbuffer.data(gh_off + 14 * ccomps * dcomps);

            auto g_xxxx_yyyyy = cbuffer.data(gh_off + 15 * ccomps * dcomps);

            auto g_xxxx_yyyyz = cbuffer.data(gh_off + 16 * ccomps * dcomps);

            auto g_xxxx_yyyzz = cbuffer.data(gh_off + 17 * ccomps * dcomps);

            auto g_xxxx_yyzzz = cbuffer.data(gh_off + 18 * ccomps * dcomps);

            auto g_xxxx_yzzzz = cbuffer.data(gh_off + 19 * ccomps * dcomps);

            auto g_xxxx_zzzzz = cbuffer.data(gh_off + 20 * ccomps * dcomps);

            auto g_xxxy_xxxxx = cbuffer.data(gh_off + 21 * ccomps * dcomps);

            auto g_xxxy_xxxxy = cbuffer.data(gh_off + 22 * ccomps * dcomps);

            auto g_xxxy_xxxxz = cbuffer.data(gh_off + 23 * ccomps * dcomps);

            auto g_xxxy_xxxyy = cbuffer.data(gh_off + 24 * ccomps * dcomps);

            auto g_xxxy_xxxyz = cbuffer.data(gh_off + 25 * ccomps * dcomps);

            auto g_xxxy_xxxzz = cbuffer.data(gh_off + 26 * ccomps * dcomps);

            auto g_xxxy_xxyyy = cbuffer.data(gh_off + 27 * ccomps * dcomps);

            auto g_xxxy_xxyyz = cbuffer.data(gh_off + 28 * ccomps * dcomps);

            auto g_xxxy_xxyzz = cbuffer.data(gh_off + 29 * ccomps * dcomps);

            auto g_xxxy_xxzzz = cbuffer.data(gh_off + 30 * ccomps * dcomps);

            auto g_xxxy_xyyyy = cbuffer.data(gh_off + 31 * ccomps * dcomps);

            auto g_xxxy_xyyyz = cbuffer.data(gh_off + 32 * ccomps * dcomps);

            auto g_xxxy_xyyzz = cbuffer.data(gh_off + 33 * ccomps * dcomps);

            auto g_xxxy_xyzzz = cbuffer.data(gh_off + 34 * ccomps * dcomps);

            auto g_xxxy_xzzzz = cbuffer.data(gh_off + 35 * ccomps * dcomps);

            auto g_xxxy_yyyyy = cbuffer.data(gh_off + 36 * ccomps * dcomps);

            auto g_xxxy_yyyyz = cbuffer.data(gh_off + 37 * ccomps * dcomps);

            auto g_xxxy_yyyzz = cbuffer.data(gh_off + 38 * ccomps * dcomps);

            auto g_xxxy_yyzzz = cbuffer.data(gh_off + 39 * ccomps * dcomps);

            auto g_xxxy_yzzzz = cbuffer.data(gh_off + 40 * ccomps * dcomps);

            auto g_xxxy_zzzzz = cbuffer.data(gh_off + 41 * ccomps * dcomps);

            auto g_xxxz_xxxxx = cbuffer.data(gh_off + 42 * ccomps * dcomps);

            auto g_xxxz_xxxxy = cbuffer.data(gh_off + 43 * ccomps * dcomps);

            auto g_xxxz_xxxxz = cbuffer.data(gh_off + 44 * ccomps * dcomps);

            auto g_xxxz_xxxyy = cbuffer.data(gh_off + 45 * ccomps * dcomps);

            auto g_xxxz_xxxyz = cbuffer.data(gh_off + 46 * ccomps * dcomps);

            auto g_xxxz_xxxzz = cbuffer.data(gh_off + 47 * ccomps * dcomps);

            auto g_xxxz_xxyyy = cbuffer.data(gh_off + 48 * ccomps * dcomps);

            auto g_xxxz_xxyyz = cbuffer.data(gh_off + 49 * ccomps * dcomps);

            auto g_xxxz_xxyzz = cbuffer.data(gh_off + 50 * ccomps * dcomps);

            auto g_xxxz_xxzzz = cbuffer.data(gh_off + 51 * ccomps * dcomps);

            auto g_xxxz_xyyyy = cbuffer.data(gh_off + 52 * ccomps * dcomps);

            auto g_xxxz_xyyyz = cbuffer.data(gh_off + 53 * ccomps * dcomps);

            auto g_xxxz_xyyzz = cbuffer.data(gh_off + 54 * ccomps * dcomps);

            auto g_xxxz_xyzzz = cbuffer.data(gh_off + 55 * ccomps * dcomps);

            auto g_xxxz_xzzzz = cbuffer.data(gh_off + 56 * ccomps * dcomps);

            auto g_xxxz_yyyyy = cbuffer.data(gh_off + 57 * ccomps * dcomps);

            auto g_xxxz_yyyyz = cbuffer.data(gh_off + 58 * ccomps * dcomps);

            auto g_xxxz_yyyzz = cbuffer.data(gh_off + 59 * ccomps * dcomps);

            auto g_xxxz_yyzzz = cbuffer.data(gh_off + 60 * ccomps * dcomps);

            auto g_xxxz_yzzzz = cbuffer.data(gh_off + 61 * ccomps * dcomps);

            auto g_xxxz_zzzzz = cbuffer.data(gh_off + 62 * ccomps * dcomps);

            auto g_xxyy_xxxxx = cbuffer.data(gh_off + 63 * ccomps * dcomps);

            auto g_xxyy_xxxxy = cbuffer.data(gh_off + 64 * ccomps * dcomps);

            auto g_xxyy_xxxxz = cbuffer.data(gh_off + 65 * ccomps * dcomps);

            auto g_xxyy_xxxyy = cbuffer.data(gh_off + 66 * ccomps * dcomps);

            auto g_xxyy_xxxyz = cbuffer.data(gh_off + 67 * ccomps * dcomps);

            auto g_xxyy_xxxzz = cbuffer.data(gh_off + 68 * ccomps * dcomps);

            auto g_xxyy_xxyyy = cbuffer.data(gh_off + 69 * ccomps * dcomps);

            auto g_xxyy_xxyyz = cbuffer.data(gh_off + 70 * ccomps * dcomps);

            auto g_xxyy_xxyzz = cbuffer.data(gh_off + 71 * ccomps * dcomps);

            auto g_xxyy_xxzzz = cbuffer.data(gh_off + 72 * ccomps * dcomps);

            auto g_xxyy_xyyyy = cbuffer.data(gh_off + 73 * ccomps * dcomps);

            auto g_xxyy_xyyyz = cbuffer.data(gh_off + 74 * ccomps * dcomps);

            auto g_xxyy_xyyzz = cbuffer.data(gh_off + 75 * ccomps * dcomps);

            auto g_xxyy_xyzzz = cbuffer.data(gh_off + 76 * ccomps * dcomps);

            auto g_xxyy_xzzzz = cbuffer.data(gh_off + 77 * ccomps * dcomps);

            auto g_xxyy_yyyyy = cbuffer.data(gh_off + 78 * ccomps * dcomps);

            auto g_xxyy_yyyyz = cbuffer.data(gh_off + 79 * ccomps * dcomps);

            auto g_xxyy_yyyzz = cbuffer.data(gh_off + 80 * ccomps * dcomps);

            auto g_xxyy_yyzzz = cbuffer.data(gh_off + 81 * ccomps * dcomps);

            auto g_xxyy_yzzzz = cbuffer.data(gh_off + 82 * ccomps * dcomps);

            auto g_xxyy_zzzzz = cbuffer.data(gh_off + 83 * ccomps * dcomps);

            auto g_xxyz_xxxxx = cbuffer.data(gh_off + 84 * ccomps * dcomps);

            auto g_xxyz_xxxxy = cbuffer.data(gh_off + 85 * ccomps * dcomps);

            auto g_xxyz_xxxxz = cbuffer.data(gh_off + 86 * ccomps * dcomps);

            auto g_xxyz_xxxyy = cbuffer.data(gh_off + 87 * ccomps * dcomps);

            auto g_xxyz_xxxyz = cbuffer.data(gh_off + 88 * ccomps * dcomps);

            auto g_xxyz_xxxzz = cbuffer.data(gh_off + 89 * ccomps * dcomps);

            auto g_xxyz_xxyyy = cbuffer.data(gh_off + 90 * ccomps * dcomps);

            auto g_xxyz_xxyyz = cbuffer.data(gh_off + 91 * ccomps * dcomps);

            auto g_xxyz_xxyzz = cbuffer.data(gh_off + 92 * ccomps * dcomps);

            auto g_xxyz_xxzzz = cbuffer.data(gh_off + 93 * ccomps * dcomps);

            auto g_xxyz_xyyyy = cbuffer.data(gh_off + 94 * ccomps * dcomps);

            auto g_xxyz_xyyyz = cbuffer.data(gh_off + 95 * ccomps * dcomps);

            auto g_xxyz_xyyzz = cbuffer.data(gh_off + 96 * ccomps * dcomps);

            auto g_xxyz_xyzzz = cbuffer.data(gh_off + 97 * ccomps * dcomps);

            auto g_xxyz_xzzzz = cbuffer.data(gh_off + 98 * ccomps * dcomps);

            auto g_xxyz_yyyyy = cbuffer.data(gh_off + 99 * ccomps * dcomps);

            auto g_xxyz_yyyyz = cbuffer.data(gh_off + 100 * ccomps * dcomps);

            auto g_xxyz_yyyzz = cbuffer.data(gh_off + 101 * ccomps * dcomps);

            auto g_xxyz_yyzzz = cbuffer.data(gh_off + 102 * ccomps * dcomps);

            auto g_xxyz_yzzzz = cbuffer.data(gh_off + 103 * ccomps * dcomps);

            auto g_xxyz_zzzzz = cbuffer.data(gh_off + 104 * ccomps * dcomps);

            auto g_xxzz_xxxxx = cbuffer.data(gh_off + 105 * ccomps * dcomps);

            auto g_xxzz_xxxxy = cbuffer.data(gh_off + 106 * ccomps * dcomps);

            auto g_xxzz_xxxxz = cbuffer.data(gh_off + 107 * ccomps * dcomps);

            auto g_xxzz_xxxyy = cbuffer.data(gh_off + 108 * ccomps * dcomps);

            auto g_xxzz_xxxyz = cbuffer.data(gh_off + 109 * ccomps * dcomps);

            auto g_xxzz_xxxzz = cbuffer.data(gh_off + 110 * ccomps * dcomps);

            auto g_xxzz_xxyyy = cbuffer.data(gh_off + 111 * ccomps * dcomps);

            auto g_xxzz_xxyyz = cbuffer.data(gh_off + 112 * ccomps * dcomps);

            auto g_xxzz_xxyzz = cbuffer.data(gh_off + 113 * ccomps * dcomps);

            auto g_xxzz_xxzzz = cbuffer.data(gh_off + 114 * ccomps * dcomps);

            auto g_xxzz_xyyyy = cbuffer.data(gh_off + 115 * ccomps * dcomps);

            auto g_xxzz_xyyyz = cbuffer.data(gh_off + 116 * ccomps * dcomps);

            auto g_xxzz_xyyzz = cbuffer.data(gh_off + 117 * ccomps * dcomps);

            auto g_xxzz_xyzzz = cbuffer.data(gh_off + 118 * ccomps * dcomps);

            auto g_xxzz_xzzzz = cbuffer.data(gh_off + 119 * ccomps * dcomps);

            auto g_xxzz_yyyyy = cbuffer.data(gh_off + 120 * ccomps * dcomps);

            auto g_xxzz_yyyyz = cbuffer.data(gh_off + 121 * ccomps * dcomps);

            auto g_xxzz_yyyzz = cbuffer.data(gh_off + 122 * ccomps * dcomps);

            auto g_xxzz_yyzzz = cbuffer.data(gh_off + 123 * ccomps * dcomps);

            auto g_xxzz_yzzzz = cbuffer.data(gh_off + 124 * ccomps * dcomps);

            auto g_xxzz_zzzzz = cbuffer.data(gh_off + 125 * ccomps * dcomps);

            auto g_xyyy_xxxxx = cbuffer.data(gh_off + 126 * ccomps * dcomps);

            auto g_xyyy_xxxxy = cbuffer.data(gh_off + 127 * ccomps * dcomps);

            auto g_xyyy_xxxxz = cbuffer.data(gh_off + 128 * ccomps * dcomps);

            auto g_xyyy_xxxyy = cbuffer.data(gh_off + 129 * ccomps * dcomps);

            auto g_xyyy_xxxyz = cbuffer.data(gh_off + 130 * ccomps * dcomps);

            auto g_xyyy_xxxzz = cbuffer.data(gh_off + 131 * ccomps * dcomps);

            auto g_xyyy_xxyyy = cbuffer.data(gh_off + 132 * ccomps * dcomps);

            auto g_xyyy_xxyyz = cbuffer.data(gh_off + 133 * ccomps * dcomps);

            auto g_xyyy_xxyzz = cbuffer.data(gh_off + 134 * ccomps * dcomps);

            auto g_xyyy_xxzzz = cbuffer.data(gh_off + 135 * ccomps * dcomps);

            auto g_xyyy_xyyyy = cbuffer.data(gh_off + 136 * ccomps * dcomps);

            auto g_xyyy_xyyyz = cbuffer.data(gh_off + 137 * ccomps * dcomps);

            auto g_xyyy_xyyzz = cbuffer.data(gh_off + 138 * ccomps * dcomps);

            auto g_xyyy_xyzzz = cbuffer.data(gh_off + 139 * ccomps * dcomps);

            auto g_xyyy_xzzzz = cbuffer.data(gh_off + 140 * ccomps * dcomps);

            auto g_xyyy_yyyyy = cbuffer.data(gh_off + 141 * ccomps * dcomps);

            auto g_xyyy_yyyyz = cbuffer.data(gh_off + 142 * ccomps * dcomps);

            auto g_xyyy_yyyzz = cbuffer.data(gh_off + 143 * ccomps * dcomps);

            auto g_xyyy_yyzzz = cbuffer.data(gh_off + 144 * ccomps * dcomps);

            auto g_xyyy_yzzzz = cbuffer.data(gh_off + 145 * ccomps * dcomps);

            auto g_xyyy_zzzzz = cbuffer.data(gh_off + 146 * ccomps * dcomps);

            auto g_xyyz_xxxxx = cbuffer.data(gh_off + 147 * ccomps * dcomps);

            auto g_xyyz_xxxxy = cbuffer.data(gh_off + 148 * ccomps * dcomps);

            auto g_xyyz_xxxxz = cbuffer.data(gh_off + 149 * ccomps * dcomps);

            auto g_xyyz_xxxyy = cbuffer.data(gh_off + 150 * ccomps * dcomps);

            auto g_xyyz_xxxyz = cbuffer.data(gh_off + 151 * ccomps * dcomps);

            auto g_xyyz_xxxzz = cbuffer.data(gh_off + 152 * ccomps * dcomps);

            auto g_xyyz_xxyyy = cbuffer.data(gh_off + 153 * ccomps * dcomps);

            auto g_xyyz_xxyyz = cbuffer.data(gh_off + 154 * ccomps * dcomps);

            auto g_xyyz_xxyzz = cbuffer.data(gh_off + 155 * ccomps * dcomps);

            auto g_xyyz_xxzzz = cbuffer.data(gh_off + 156 * ccomps * dcomps);

            auto g_xyyz_xyyyy = cbuffer.data(gh_off + 157 * ccomps * dcomps);

            auto g_xyyz_xyyyz = cbuffer.data(gh_off + 158 * ccomps * dcomps);

            auto g_xyyz_xyyzz = cbuffer.data(gh_off + 159 * ccomps * dcomps);

            auto g_xyyz_xyzzz = cbuffer.data(gh_off + 160 * ccomps * dcomps);

            auto g_xyyz_xzzzz = cbuffer.data(gh_off + 161 * ccomps * dcomps);

            auto g_xyyz_yyyyy = cbuffer.data(gh_off + 162 * ccomps * dcomps);

            auto g_xyyz_yyyyz = cbuffer.data(gh_off + 163 * ccomps * dcomps);

            auto g_xyyz_yyyzz = cbuffer.data(gh_off + 164 * ccomps * dcomps);

            auto g_xyyz_yyzzz = cbuffer.data(gh_off + 165 * ccomps * dcomps);

            auto g_xyyz_yzzzz = cbuffer.data(gh_off + 166 * ccomps * dcomps);

            auto g_xyyz_zzzzz = cbuffer.data(gh_off + 167 * ccomps * dcomps);

            auto g_xyzz_xxxxx = cbuffer.data(gh_off + 168 * ccomps * dcomps);

            auto g_xyzz_xxxxy = cbuffer.data(gh_off + 169 * ccomps * dcomps);

            auto g_xyzz_xxxxz = cbuffer.data(gh_off + 170 * ccomps * dcomps);

            auto g_xyzz_xxxyy = cbuffer.data(gh_off + 171 * ccomps * dcomps);

            auto g_xyzz_xxxyz = cbuffer.data(gh_off + 172 * ccomps * dcomps);

            auto g_xyzz_xxxzz = cbuffer.data(gh_off + 173 * ccomps * dcomps);

            auto g_xyzz_xxyyy = cbuffer.data(gh_off + 174 * ccomps * dcomps);

            auto g_xyzz_xxyyz = cbuffer.data(gh_off + 175 * ccomps * dcomps);

            auto g_xyzz_xxyzz = cbuffer.data(gh_off + 176 * ccomps * dcomps);

            auto g_xyzz_xxzzz = cbuffer.data(gh_off + 177 * ccomps * dcomps);

            auto g_xyzz_xyyyy = cbuffer.data(gh_off + 178 * ccomps * dcomps);

            auto g_xyzz_xyyyz = cbuffer.data(gh_off + 179 * ccomps * dcomps);

            auto g_xyzz_xyyzz = cbuffer.data(gh_off + 180 * ccomps * dcomps);

            auto g_xyzz_xyzzz = cbuffer.data(gh_off + 181 * ccomps * dcomps);

            auto g_xyzz_xzzzz = cbuffer.data(gh_off + 182 * ccomps * dcomps);

            auto g_xyzz_yyyyy = cbuffer.data(gh_off + 183 * ccomps * dcomps);

            auto g_xyzz_yyyyz = cbuffer.data(gh_off + 184 * ccomps * dcomps);

            auto g_xyzz_yyyzz = cbuffer.data(gh_off + 185 * ccomps * dcomps);

            auto g_xyzz_yyzzz = cbuffer.data(gh_off + 186 * ccomps * dcomps);

            auto g_xyzz_yzzzz = cbuffer.data(gh_off + 187 * ccomps * dcomps);

            auto g_xyzz_zzzzz = cbuffer.data(gh_off + 188 * ccomps * dcomps);

            auto g_xzzz_xxxxx = cbuffer.data(gh_off + 189 * ccomps * dcomps);

            auto g_xzzz_xxxxy = cbuffer.data(gh_off + 190 * ccomps * dcomps);

            auto g_xzzz_xxxxz = cbuffer.data(gh_off + 191 * ccomps * dcomps);

            auto g_xzzz_xxxyy = cbuffer.data(gh_off + 192 * ccomps * dcomps);

            auto g_xzzz_xxxyz = cbuffer.data(gh_off + 193 * ccomps * dcomps);

            auto g_xzzz_xxxzz = cbuffer.data(gh_off + 194 * ccomps * dcomps);

            auto g_xzzz_xxyyy = cbuffer.data(gh_off + 195 * ccomps * dcomps);

            auto g_xzzz_xxyyz = cbuffer.data(gh_off + 196 * ccomps * dcomps);

            auto g_xzzz_xxyzz = cbuffer.data(gh_off + 197 * ccomps * dcomps);

            auto g_xzzz_xxzzz = cbuffer.data(gh_off + 198 * ccomps * dcomps);

            auto g_xzzz_xyyyy = cbuffer.data(gh_off + 199 * ccomps * dcomps);

            auto g_xzzz_xyyyz = cbuffer.data(gh_off + 200 * ccomps * dcomps);

            auto g_xzzz_xyyzz = cbuffer.data(gh_off + 201 * ccomps * dcomps);

            auto g_xzzz_xyzzz = cbuffer.data(gh_off + 202 * ccomps * dcomps);

            auto g_xzzz_xzzzz = cbuffer.data(gh_off + 203 * ccomps * dcomps);

            auto g_xzzz_yyyyy = cbuffer.data(gh_off + 204 * ccomps * dcomps);

            auto g_xzzz_yyyyz = cbuffer.data(gh_off + 205 * ccomps * dcomps);

            auto g_xzzz_yyyzz = cbuffer.data(gh_off + 206 * ccomps * dcomps);

            auto g_xzzz_yyzzz = cbuffer.data(gh_off + 207 * ccomps * dcomps);

            auto g_xzzz_yzzzz = cbuffer.data(gh_off + 208 * ccomps * dcomps);

            auto g_xzzz_zzzzz = cbuffer.data(gh_off + 209 * ccomps * dcomps);

            auto g_yyyy_xxxxx = cbuffer.data(gh_off + 210 * ccomps * dcomps);

            auto g_yyyy_xxxxy = cbuffer.data(gh_off + 211 * ccomps * dcomps);

            auto g_yyyy_xxxxz = cbuffer.data(gh_off + 212 * ccomps * dcomps);

            auto g_yyyy_xxxyy = cbuffer.data(gh_off + 213 * ccomps * dcomps);

            auto g_yyyy_xxxyz = cbuffer.data(gh_off + 214 * ccomps * dcomps);

            auto g_yyyy_xxxzz = cbuffer.data(gh_off + 215 * ccomps * dcomps);

            auto g_yyyy_xxyyy = cbuffer.data(gh_off + 216 * ccomps * dcomps);

            auto g_yyyy_xxyyz = cbuffer.data(gh_off + 217 * ccomps * dcomps);

            auto g_yyyy_xxyzz = cbuffer.data(gh_off + 218 * ccomps * dcomps);

            auto g_yyyy_xxzzz = cbuffer.data(gh_off + 219 * ccomps * dcomps);

            auto g_yyyy_xyyyy = cbuffer.data(gh_off + 220 * ccomps * dcomps);

            auto g_yyyy_xyyyz = cbuffer.data(gh_off + 221 * ccomps * dcomps);

            auto g_yyyy_xyyzz = cbuffer.data(gh_off + 222 * ccomps * dcomps);

            auto g_yyyy_xyzzz = cbuffer.data(gh_off + 223 * ccomps * dcomps);

            auto g_yyyy_xzzzz = cbuffer.data(gh_off + 224 * ccomps * dcomps);

            auto g_yyyy_yyyyy = cbuffer.data(gh_off + 225 * ccomps * dcomps);

            auto g_yyyy_yyyyz = cbuffer.data(gh_off + 226 * ccomps * dcomps);

            auto g_yyyy_yyyzz = cbuffer.data(gh_off + 227 * ccomps * dcomps);

            auto g_yyyy_yyzzz = cbuffer.data(gh_off + 228 * ccomps * dcomps);

            auto g_yyyy_yzzzz = cbuffer.data(gh_off + 229 * ccomps * dcomps);

            auto g_yyyy_zzzzz = cbuffer.data(gh_off + 230 * ccomps * dcomps);

            auto g_yyyz_xxxxx = cbuffer.data(gh_off + 231 * ccomps * dcomps);

            auto g_yyyz_xxxxy = cbuffer.data(gh_off + 232 * ccomps * dcomps);

            auto g_yyyz_xxxxz = cbuffer.data(gh_off + 233 * ccomps * dcomps);

            auto g_yyyz_xxxyy = cbuffer.data(gh_off + 234 * ccomps * dcomps);

            auto g_yyyz_xxxyz = cbuffer.data(gh_off + 235 * ccomps * dcomps);

            auto g_yyyz_xxxzz = cbuffer.data(gh_off + 236 * ccomps * dcomps);

            auto g_yyyz_xxyyy = cbuffer.data(gh_off + 237 * ccomps * dcomps);

            auto g_yyyz_xxyyz = cbuffer.data(gh_off + 238 * ccomps * dcomps);

            auto g_yyyz_xxyzz = cbuffer.data(gh_off + 239 * ccomps * dcomps);

            auto g_yyyz_xxzzz = cbuffer.data(gh_off + 240 * ccomps * dcomps);

            auto g_yyyz_xyyyy = cbuffer.data(gh_off + 241 * ccomps * dcomps);

            auto g_yyyz_xyyyz = cbuffer.data(gh_off + 242 * ccomps * dcomps);

            auto g_yyyz_xyyzz = cbuffer.data(gh_off + 243 * ccomps * dcomps);

            auto g_yyyz_xyzzz = cbuffer.data(gh_off + 244 * ccomps * dcomps);

            auto g_yyyz_xzzzz = cbuffer.data(gh_off + 245 * ccomps * dcomps);

            auto g_yyyz_yyyyy = cbuffer.data(gh_off + 246 * ccomps * dcomps);

            auto g_yyyz_yyyyz = cbuffer.data(gh_off + 247 * ccomps * dcomps);

            auto g_yyyz_yyyzz = cbuffer.data(gh_off + 248 * ccomps * dcomps);

            auto g_yyyz_yyzzz = cbuffer.data(gh_off + 249 * ccomps * dcomps);

            auto g_yyyz_yzzzz = cbuffer.data(gh_off + 250 * ccomps * dcomps);

            auto g_yyyz_zzzzz = cbuffer.data(gh_off + 251 * ccomps * dcomps);

            auto g_yyzz_xxxxx = cbuffer.data(gh_off + 252 * ccomps * dcomps);

            auto g_yyzz_xxxxy = cbuffer.data(gh_off + 253 * ccomps * dcomps);

            auto g_yyzz_xxxxz = cbuffer.data(gh_off + 254 * ccomps * dcomps);

            auto g_yyzz_xxxyy = cbuffer.data(gh_off + 255 * ccomps * dcomps);

            auto g_yyzz_xxxyz = cbuffer.data(gh_off + 256 * ccomps * dcomps);

            auto g_yyzz_xxxzz = cbuffer.data(gh_off + 257 * ccomps * dcomps);

            auto g_yyzz_xxyyy = cbuffer.data(gh_off + 258 * ccomps * dcomps);

            auto g_yyzz_xxyyz = cbuffer.data(gh_off + 259 * ccomps * dcomps);

            auto g_yyzz_xxyzz = cbuffer.data(gh_off + 260 * ccomps * dcomps);

            auto g_yyzz_xxzzz = cbuffer.data(gh_off + 261 * ccomps * dcomps);

            auto g_yyzz_xyyyy = cbuffer.data(gh_off + 262 * ccomps * dcomps);

            auto g_yyzz_xyyyz = cbuffer.data(gh_off + 263 * ccomps * dcomps);

            auto g_yyzz_xyyzz = cbuffer.data(gh_off + 264 * ccomps * dcomps);

            auto g_yyzz_xyzzz = cbuffer.data(gh_off + 265 * ccomps * dcomps);

            auto g_yyzz_xzzzz = cbuffer.data(gh_off + 266 * ccomps * dcomps);

            auto g_yyzz_yyyyy = cbuffer.data(gh_off + 267 * ccomps * dcomps);

            auto g_yyzz_yyyyz = cbuffer.data(gh_off + 268 * ccomps * dcomps);

            auto g_yyzz_yyyzz = cbuffer.data(gh_off + 269 * ccomps * dcomps);

            auto g_yyzz_yyzzz = cbuffer.data(gh_off + 270 * ccomps * dcomps);

            auto g_yyzz_yzzzz = cbuffer.data(gh_off + 271 * ccomps * dcomps);

            auto g_yyzz_zzzzz = cbuffer.data(gh_off + 272 * ccomps * dcomps);

            auto g_yzzz_xxxxx = cbuffer.data(gh_off + 273 * ccomps * dcomps);

            auto g_yzzz_xxxxy = cbuffer.data(gh_off + 274 * ccomps * dcomps);

            auto g_yzzz_xxxxz = cbuffer.data(gh_off + 275 * ccomps * dcomps);

            auto g_yzzz_xxxyy = cbuffer.data(gh_off + 276 * ccomps * dcomps);

            auto g_yzzz_xxxyz = cbuffer.data(gh_off + 277 * ccomps * dcomps);

            auto g_yzzz_xxxzz = cbuffer.data(gh_off + 278 * ccomps * dcomps);

            auto g_yzzz_xxyyy = cbuffer.data(gh_off + 279 * ccomps * dcomps);

            auto g_yzzz_xxyyz = cbuffer.data(gh_off + 280 * ccomps * dcomps);

            auto g_yzzz_xxyzz = cbuffer.data(gh_off + 281 * ccomps * dcomps);

            auto g_yzzz_xxzzz = cbuffer.data(gh_off + 282 * ccomps * dcomps);

            auto g_yzzz_xyyyy = cbuffer.data(gh_off + 283 * ccomps * dcomps);

            auto g_yzzz_xyyyz = cbuffer.data(gh_off + 284 * ccomps * dcomps);

            auto g_yzzz_xyyzz = cbuffer.data(gh_off + 285 * ccomps * dcomps);

            auto g_yzzz_xyzzz = cbuffer.data(gh_off + 286 * ccomps * dcomps);

            auto g_yzzz_xzzzz = cbuffer.data(gh_off + 287 * ccomps * dcomps);

            auto g_yzzz_yyyyy = cbuffer.data(gh_off + 288 * ccomps * dcomps);

            auto g_yzzz_yyyyz = cbuffer.data(gh_off + 289 * ccomps * dcomps);

            auto g_yzzz_yyyzz = cbuffer.data(gh_off + 290 * ccomps * dcomps);

            auto g_yzzz_yyzzz = cbuffer.data(gh_off + 291 * ccomps * dcomps);

            auto g_yzzz_yzzzz = cbuffer.data(gh_off + 292 * ccomps * dcomps);

            auto g_yzzz_zzzzz = cbuffer.data(gh_off + 293 * ccomps * dcomps);

            auto g_zzzz_xxxxx = cbuffer.data(gh_off + 294 * ccomps * dcomps);

            auto g_zzzz_xxxxy = cbuffer.data(gh_off + 295 * ccomps * dcomps);

            auto g_zzzz_xxxxz = cbuffer.data(gh_off + 296 * ccomps * dcomps);

            auto g_zzzz_xxxyy = cbuffer.data(gh_off + 297 * ccomps * dcomps);

            auto g_zzzz_xxxyz = cbuffer.data(gh_off + 298 * ccomps * dcomps);

            auto g_zzzz_xxxzz = cbuffer.data(gh_off + 299 * ccomps * dcomps);

            auto g_zzzz_xxyyy = cbuffer.data(gh_off + 300 * ccomps * dcomps);

            auto g_zzzz_xxyyz = cbuffer.data(gh_off + 301 * ccomps * dcomps);

            auto g_zzzz_xxyzz = cbuffer.data(gh_off + 302 * ccomps * dcomps);

            auto g_zzzz_xxzzz = cbuffer.data(gh_off + 303 * ccomps * dcomps);

            auto g_zzzz_xyyyy = cbuffer.data(gh_off + 304 * ccomps * dcomps);

            auto g_zzzz_xyyyz = cbuffer.data(gh_off + 305 * ccomps * dcomps);

            auto g_zzzz_xyyzz = cbuffer.data(gh_off + 306 * ccomps * dcomps);

            auto g_zzzz_xyzzz = cbuffer.data(gh_off + 307 * ccomps * dcomps);

            auto g_zzzz_xzzzz = cbuffer.data(gh_off + 308 * ccomps * dcomps);

            auto g_zzzz_yyyyy = cbuffer.data(gh_off + 309 * ccomps * dcomps);

            auto g_zzzz_yyyyz = cbuffer.data(gh_off + 310 * ccomps * dcomps);

            auto g_zzzz_yyyzz = cbuffer.data(gh_off + 311 * ccomps * dcomps);

            auto g_zzzz_yyzzz = cbuffer.data(gh_off + 312 * ccomps * dcomps);

            auto g_zzzz_yzzzz = cbuffer.data(gh_off + 313 * ccomps * dcomps);

            auto g_zzzz_zzzzz = cbuffer.data(gh_off + 314 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GHSS

            const auto gh_geom_01_off = idx_geom_01_ghxx + i * dcomps + j;

            auto g_0_x_xxxx_xxxxx = cbuffer.data(gh_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxy = cbuffer.data(gh_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxz = cbuffer.data(gh_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyy = cbuffer.data(gh_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyz = cbuffer.data(gh_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxzz = cbuffer.data(gh_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyy = cbuffer.data(gh_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyz = cbuffer.data(gh_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyzz = cbuffer.data(gh_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxx_xxzzz = cbuffer.data(gh_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyy = cbuffer.data(gh_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyz = cbuffer.data(gh_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyzz = cbuffer.data(gh_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxx_xyzzz = cbuffer.data(gh_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxx_xzzzz = cbuffer.data(gh_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyy = cbuffer.data(gh_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyz = cbuffer.data(gh_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyzz = cbuffer.data(gh_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxx_yyzzz = cbuffer.data(gh_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxx_yzzzz = cbuffer.data(gh_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxx_zzzzz = cbuffer.data(gh_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxx = cbuffer.data(gh_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxy = cbuffer.data(gh_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxz = cbuffer.data(gh_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyy = cbuffer.data(gh_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyz = cbuffer.data(gh_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxzz = cbuffer.data(gh_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyy = cbuffer.data(gh_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyz = cbuffer.data(gh_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyzz = cbuffer.data(gh_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxy_xxzzz = cbuffer.data(gh_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyy = cbuffer.data(gh_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyz = cbuffer.data(gh_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyzz = cbuffer.data(gh_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxy_xyzzz = cbuffer.data(gh_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxy_xzzzz = cbuffer.data(gh_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyy = cbuffer.data(gh_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyz = cbuffer.data(gh_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyzz = cbuffer.data(gh_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxy_yyzzz = cbuffer.data(gh_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxy_yzzzz = cbuffer.data(gh_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxy_zzzzz = cbuffer.data(gh_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxx = cbuffer.data(gh_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxy = cbuffer.data(gh_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxz = cbuffer.data(gh_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyy = cbuffer.data(gh_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyz = cbuffer.data(gh_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxzz = cbuffer.data(gh_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyy = cbuffer.data(gh_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyz = cbuffer.data(gh_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyzz = cbuffer.data(gh_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxz_xxzzz = cbuffer.data(gh_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyy = cbuffer.data(gh_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyz = cbuffer.data(gh_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyzz = cbuffer.data(gh_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxz_xyzzz = cbuffer.data(gh_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxz_xzzzz = cbuffer.data(gh_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyy = cbuffer.data(gh_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyz = cbuffer.data(gh_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyzz = cbuffer.data(gh_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxxz_yyzzz = cbuffer.data(gh_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxz_yzzzz = cbuffer.data(gh_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxz_zzzzz = cbuffer.data(gh_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxx = cbuffer.data(gh_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxy = cbuffer.data(gh_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxz = cbuffer.data(gh_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyy = cbuffer.data(gh_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyz = cbuffer.data(gh_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxzz = cbuffer.data(gh_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyy = cbuffer.data(gh_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyz = cbuffer.data(gh_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyzz = cbuffer.data(gh_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxyy_xxzzz = cbuffer.data(gh_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyy = cbuffer.data(gh_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyz = cbuffer.data(gh_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyzz = cbuffer.data(gh_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxyy_xyzzz = cbuffer.data(gh_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxyy_xzzzz = cbuffer.data(gh_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyy = cbuffer.data(gh_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyz = cbuffer.data(gh_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyzz = cbuffer.data(gh_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxyy_yyzzz = cbuffer.data(gh_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxyy_yzzzz = cbuffer.data(gh_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxyy_zzzzz = cbuffer.data(gh_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxx = cbuffer.data(gh_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxy = cbuffer.data(gh_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxz = cbuffer.data(gh_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyy = cbuffer.data(gh_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyz = cbuffer.data(gh_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxzz = cbuffer.data(gh_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyy = cbuffer.data(gh_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyz = cbuffer.data(gh_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyzz = cbuffer.data(gh_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxyz_xxzzz = cbuffer.data(gh_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyy = cbuffer.data(gh_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyz = cbuffer.data(gh_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyzz = cbuffer.data(gh_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxyz_xyzzz = cbuffer.data(gh_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxyz_xzzzz = cbuffer.data(gh_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyy = cbuffer.data(gh_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyz = cbuffer.data(gh_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyzz = cbuffer.data(gh_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxyz_yyzzz = cbuffer.data(gh_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxyz_yzzzz = cbuffer.data(gh_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxyz_zzzzz = cbuffer.data(gh_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxx = cbuffer.data(gh_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxy = cbuffer.data(gh_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxz = cbuffer.data(gh_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyy = cbuffer.data(gh_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyz = cbuffer.data(gh_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxzz = cbuffer.data(gh_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyy = cbuffer.data(gh_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyz = cbuffer.data(gh_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyzz = cbuffer.data(gh_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxzz_xxzzz = cbuffer.data(gh_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyy = cbuffer.data(gh_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyz = cbuffer.data(gh_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyzz = cbuffer.data(gh_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxzz_xyzzz = cbuffer.data(gh_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxzz_xzzzz = cbuffer.data(gh_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyy = cbuffer.data(gh_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyz = cbuffer.data(gh_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyzz = cbuffer.data(gh_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxzz_yyzzz = cbuffer.data(gh_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxzz_yzzzz = cbuffer.data(gh_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxzz_zzzzz = cbuffer.data(gh_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxx = cbuffer.data(gh_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxy = cbuffer.data(gh_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxz = cbuffer.data(gh_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyy = cbuffer.data(gh_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyz = cbuffer.data(gh_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxzz = cbuffer.data(gh_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyy = cbuffer.data(gh_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyz = cbuffer.data(gh_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyzz = cbuffer.data(gh_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xyyy_xxzzz = cbuffer.data(gh_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyy = cbuffer.data(gh_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyz = cbuffer.data(gh_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyzz = cbuffer.data(gh_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xyyy_xyzzz = cbuffer.data(gh_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xyyy_xzzzz = cbuffer.data(gh_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyy = cbuffer.data(gh_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyz = cbuffer.data(gh_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyzz = cbuffer.data(gh_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xyyy_yyzzz = cbuffer.data(gh_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xyyy_yzzzz = cbuffer.data(gh_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xyyy_zzzzz = cbuffer.data(gh_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxx = cbuffer.data(gh_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxy = cbuffer.data(gh_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxz = cbuffer.data(gh_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyy = cbuffer.data(gh_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyz = cbuffer.data(gh_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxzz = cbuffer.data(gh_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyy = cbuffer.data(gh_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyz = cbuffer.data(gh_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyzz = cbuffer.data(gh_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xyyz_xxzzz = cbuffer.data(gh_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyy = cbuffer.data(gh_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyz = cbuffer.data(gh_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyzz = cbuffer.data(gh_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xyyz_xyzzz = cbuffer.data(gh_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xyyz_xzzzz = cbuffer.data(gh_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyy = cbuffer.data(gh_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyz = cbuffer.data(gh_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyzz = cbuffer.data(gh_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xyyz_yyzzz = cbuffer.data(gh_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xyyz_yzzzz = cbuffer.data(gh_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xyyz_zzzzz = cbuffer.data(gh_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxx = cbuffer.data(gh_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxy = cbuffer.data(gh_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxz = cbuffer.data(gh_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyy = cbuffer.data(gh_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyz = cbuffer.data(gh_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxzz = cbuffer.data(gh_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyy = cbuffer.data(gh_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyz = cbuffer.data(gh_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyzz = cbuffer.data(gh_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xyzz_xxzzz = cbuffer.data(gh_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyy = cbuffer.data(gh_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyz = cbuffer.data(gh_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyzz = cbuffer.data(gh_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xyzz_xyzzz = cbuffer.data(gh_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xyzz_xzzzz = cbuffer.data(gh_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyy = cbuffer.data(gh_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyz = cbuffer.data(gh_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyzz = cbuffer.data(gh_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xyzz_yyzzz = cbuffer.data(gh_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xyzz_yzzzz = cbuffer.data(gh_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xyzz_zzzzz = cbuffer.data(gh_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxx = cbuffer.data(gh_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxy = cbuffer.data(gh_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxz = cbuffer.data(gh_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyy = cbuffer.data(gh_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyz = cbuffer.data(gh_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxzz = cbuffer.data(gh_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyy = cbuffer.data(gh_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyz = cbuffer.data(gh_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyzz = cbuffer.data(gh_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xzzz_xxzzz = cbuffer.data(gh_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyy = cbuffer.data(gh_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyz = cbuffer.data(gh_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyzz = cbuffer.data(gh_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xzzz_xyzzz = cbuffer.data(gh_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xzzz_xzzzz = cbuffer.data(gh_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyy = cbuffer.data(gh_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyz = cbuffer.data(gh_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyzz = cbuffer.data(gh_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xzzz_yyzzz = cbuffer.data(gh_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xzzz_yzzzz = cbuffer.data(gh_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xzzz_zzzzz = cbuffer.data(gh_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxx = cbuffer.data(gh_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxy = cbuffer.data(gh_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxz = cbuffer.data(gh_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyy = cbuffer.data(gh_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyz = cbuffer.data(gh_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxzz = cbuffer.data(gh_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyy = cbuffer.data(gh_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyz = cbuffer.data(gh_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyzz = cbuffer.data(gh_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_yyyy_xxzzz = cbuffer.data(gh_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyy = cbuffer.data(gh_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyz = cbuffer.data(gh_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyzz = cbuffer.data(gh_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_yyyy_xyzzz = cbuffer.data(gh_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_yyyy_xzzzz = cbuffer.data(gh_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyy = cbuffer.data(gh_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyz = cbuffer.data(gh_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyzz = cbuffer.data(gh_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_yyyy_yyzzz = cbuffer.data(gh_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_yyyy_yzzzz = cbuffer.data(gh_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_yyyy_zzzzz = cbuffer.data(gh_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxx = cbuffer.data(gh_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxy = cbuffer.data(gh_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxz = cbuffer.data(gh_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyy = cbuffer.data(gh_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyz = cbuffer.data(gh_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxzz = cbuffer.data(gh_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyy = cbuffer.data(gh_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyz = cbuffer.data(gh_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyzz = cbuffer.data(gh_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_yyyz_xxzzz = cbuffer.data(gh_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyy = cbuffer.data(gh_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyz = cbuffer.data(gh_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyzz = cbuffer.data(gh_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_yyyz_xyzzz = cbuffer.data(gh_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_yyyz_xzzzz = cbuffer.data(gh_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyy = cbuffer.data(gh_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyz = cbuffer.data(gh_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyzz = cbuffer.data(gh_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_yyyz_yyzzz = cbuffer.data(gh_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_yyyz_yzzzz = cbuffer.data(gh_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_yyyz_zzzzz = cbuffer.data(gh_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxx = cbuffer.data(gh_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxy = cbuffer.data(gh_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxz = cbuffer.data(gh_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyy = cbuffer.data(gh_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyz = cbuffer.data(gh_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxzz = cbuffer.data(gh_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyy = cbuffer.data(gh_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyz = cbuffer.data(gh_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyzz = cbuffer.data(gh_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_yyzz_xxzzz = cbuffer.data(gh_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyy = cbuffer.data(gh_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyz = cbuffer.data(gh_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyzz = cbuffer.data(gh_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_yyzz_xyzzz = cbuffer.data(gh_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_yyzz_xzzzz = cbuffer.data(gh_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyy = cbuffer.data(gh_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyz = cbuffer.data(gh_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyzz = cbuffer.data(gh_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_yyzz_yyzzz = cbuffer.data(gh_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_yyzz_yzzzz = cbuffer.data(gh_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_yyzz_zzzzz = cbuffer.data(gh_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxx = cbuffer.data(gh_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxy = cbuffer.data(gh_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxz = cbuffer.data(gh_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyy = cbuffer.data(gh_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyz = cbuffer.data(gh_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxzz = cbuffer.data(gh_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyy = cbuffer.data(gh_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyz = cbuffer.data(gh_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyzz = cbuffer.data(gh_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_yzzz_xxzzz = cbuffer.data(gh_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyy = cbuffer.data(gh_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyz = cbuffer.data(gh_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyzz = cbuffer.data(gh_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_yzzz_xyzzz = cbuffer.data(gh_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_yzzz_xzzzz = cbuffer.data(gh_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyy = cbuffer.data(gh_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyz = cbuffer.data(gh_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyzz = cbuffer.data(gh_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_yzzz_yyzzz = cbuffer.data(gh_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_yzzz_yzzzz = cbuffer.data(gh_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_yzzz_zzzzz = cbuffer.data(gh_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxx = cbuffer.data(gh_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxy = cbuffer.data(gh_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxz = cbuffer.data(gh_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyy = cbuffer.data(gh_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyz = cbuffer.data(gh_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxzz = cbuffer.data(gh_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyy = cbuffer.data(gh_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyz = cbuffer.data(gh_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyzz = cbuffer.data(gh_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_zzzz_xxzzz = cbuffer.data(gh_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyy = cbuffer.data(gh_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyz = cbuffer.data(gh_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyzz = cbuffer.data(gh_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_zzzz_xyzzz = cbuffer.data(gh_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_zzzz_xzzzz = cbuffer.data(gh_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyy = cbuffer.data(gh_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyz = cbuffer.data(gh_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyzz = cbuffer.data(gh_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_zzzz_yyzzz = cbuffer.data(gh_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_zzzz_yzzzz = cbuffer.data(gh_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_zzzz_zzzzz = cbuffer.data(gh_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxx = cbuffer.data(gh_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxy = cbuffer.data(gh_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxz = cbuffer.data(gh_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyy = cbuffer.data(gh_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyz = cbuffer.data(gh_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxzz = cbuffer.data(gh_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyy = cbuffer.data(gh_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyz = cbuffer.data(gh_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyzz = cbuffer.data(gh_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_xxxx_xxzzz = cbuffer.data(gh_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyy = cbuffer.data(gh_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyz = cbuffer.data(gh_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyzz = cbuffer.data(gh_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xxxx_xyzzz = cbuffer.data(gh_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xxxx_xzzzz = cbuffer.data(gh_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyy = cbuffer.data(gh_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyz = cbuffer.data(gh_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyzz = cbuffer.data(gh_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xxxx_yyzzz = cbuffer.data(gh_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xxxx_yzzzz = cbuffer.data(gh_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xxxx_zzzzz = cbuffer.data(gh_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxx = cbuffer.data(gh_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxy = cbuffer.data(gh_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxz = cbuffer.data(gh_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyy = cbuffer.data(gh_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyz = cbuffer.data(gh_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxzz = cbuffer.data(gh_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyy = cbuffer.data(gh_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyz = cbuffer.data(gh_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyzz = cbuffer.data(gh_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_xxxy_xxzzz = cbuffer.data(gh_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyy = cbuffer.data(gh_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyz = cbuffer.data(gh_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyzz = cbuffer.data(gh_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_xxxy_xyzzz = cbuffer.data(gh_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_xxxy_xzzzz = cbuffer.data(gh_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyy = cbuffer.data(gh_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyz = cbuffer.data(gh_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyzz = cbuffer.data(gh_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_xxxy_yyzzz = cbuffer.data(gh_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_xxxy_yzzzz = cbuffer.data(gh_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_xxxy_zzzzz = cbuffer.data(gh_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxx = cbuffer.data(gh_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxy = cbuffer.data(gh_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxz = cbuffer.data(gh_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyy = cbuffer.data(gh_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyz = cbuffer.data(gh_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxzz = cbuffer.data(gh_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyy = cbuffer.data(gh_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyz = cbuffer.data(gh_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyzz = cbuffer.data(gh_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_xxxz_xxzzz = cbuffer.data(gh_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyy = cbuffer.data(gh_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyz = cbuffer.data(gh_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyzz = cbuffer.data(gh_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_xxxz_xyzzz = cbuffer.data(gh_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_xxxz_xzzzz = cbuffer.data(gh_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyy = cbuffer.data(gh_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyz = cbuffer.data(gh_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyzz = cbuffer.data(gh_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_xxxz_yyzzz = cbuffer.data(gh_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_xxxz_yzzzz = cbuffer.data(gh_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_xxxz_zzzzz = cbuffer.data(gh_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxx = cbuffer.data(gh_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxy = cbuffer.data(gh_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxz = cbuffer.data(gh_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyy = cbuffer.data(gh_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyz = cbuffer.data(gh_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxzz = cbuffer.data(gh_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyy = cbuffer.data(gh_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyz = cbuffer.data(gh_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyzz = cbuffer.data(gh_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_xxyy_xxzzz = cbuffer.data(gh_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyy = cbuffer.data(gh_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyz = cbuffer.data(gh_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyzz = cbuffer.data(gh_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_xxyy_xyzzz = cbuffer.data(gh_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_xxyy_xzzzz = cbuffer.data(gh_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyy = cbuffer.data(gh_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyz = cbuffer.data(gh_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyzz = cbuffer.data(gh_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_xxyy_yyzzz = cbuffer.data(gh_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_xxyy_yzzzz = cbuffer.data(gh_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_xxyy_zzzzz = cbuffer.data(gh_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxx = cbuffer.data(gh_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxy = cbuffer.data(gh_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxz = cbuffer.data(gh_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyy = cbuffer.data(gh_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyz = cbuffer.data(gh_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxzz = cbuffer.data(gh_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyy = cbuffer.data(gh_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyz = cbuffer.data(gh_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyzz = cbuffer.data(gh_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_xxyz_xxzzz = cbuffer.data(gh_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyy = cbuffer.data(gh_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyz = cbuffer.data(gh_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyzz = cbuffer.data(gh_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_xxyz_xyzzz = cbuffer.data(gh_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_xxyz_xzzzz = cbuffer.data(gh_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyy = cbuffer.data(gh_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyz = cbuffer.data(gh_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyzz = cbuffer.data(gh_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_xxyz_yyzzz = cbuffer.data(gh_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_xxyz_yzzzz = cbuffer.data(gh_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_xxyz_zzzzz = cbuffer.data(gh_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxx = cbuffer.data(gh_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxy = cbuffer.data(gh_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxz = cbuffer.data(gh_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyy = cbuffer.data(gh_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyz = cbuffer.data(gh_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxzz = cbuffer.data(gh_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyy = cbuffer.data(gh_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyz = cbuffer.data(gh_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyzz = cbuffer.data(gh_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_xxzz_xxzzz = cbuffer.data(gh_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyy = cbuffer.data(gh_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyz = cbuffer.data(gh_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyzz = cbuffer.data(gh_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_xxzz_xyzzz = cbuffer.data(gh_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_xxzz_xzzzz = cbuffer.data(gh_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyy = cbuffer.data(gh_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyz = cbuffer.data(gh_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyzz = cbuffer.data(gh_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_xxzz_yyzzz = cbuffer.data(gh_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_xxzz_yzzzz = cbuffer.data(gh_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_y_xxzz_zzzzz = cbuffer.data(gh_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxx = cbuffer.data(gh_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxy = cbuffer.data(gh_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxz = cbuffer.data(gh_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyy = cbuffer.data(gh_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyz = cbuffer.data(gh_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxzz = cbuffer.data(gh_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyy = cbuffer.data(gh_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyz = cbuffer.data(gh_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyzz = cbuffer.data(gh_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_y_xyyy_xxzzz = cbuffer.data(gh_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyy = cbuffer.data(gh_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyz = cbuffer.data(gh_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyzz = cbuffer.data(gh_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_xyyy_xyzzz = cbuffer.data(gh_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_xyyy_xzzzz = cbuffer.data(gh_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyy = cbuffer.data(gh_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyz = cbuffer.data(gh_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyzz = cbuffer.data(gh_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_xyyy_yyzzz = cbuffer.data(gh_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_xyyy_yzzzz = cbuffer.data(gh_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_xyyy_zzzzz = cbuffer.data(gh_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxx = cbuffer.data(gh_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxy = cbuffer.data(gh_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxz = cbuffer.data(gh_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyy = cbuffer.data(gh_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyz = cbuffer.data(gh_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxzz = cbuffer.data(gh_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyy = cbuffer.data(gh_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyz = cbuffer.data(gh_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyzz = cbuffer.data(gh_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_xyyz_xxzzz = cbuffer.data(gh_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyy = cbuffer.data(gh_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyz = cbuffer.data(gh_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyzz = cbuffer.data(gh_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_xyyz_xyzzz = cbuffer.data(gh_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_xyyz_xzzzz = cbuffer.data(gh_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyy = cbuffer.data(gh_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyz = cbuffer.data(gh_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyzz = cbuffer.data(gh_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_y_xyyz_yyzzz = cbuffer.data(gh_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_xyyz_yzzzz = cbuffer.data(gh_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_xyyz_zzzzz = cbuffer.data(gh_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxx = cbuffer.data(gh_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxy = cbuffer.data(gh_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxz = cbuffer.data(gh_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyy = cbuffer.data(gh_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyz = cbuffer.data(gh_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxzz = cbuffer.data(gh_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyy = cbuffer.data(gh_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyz = cbuffer.data(gh_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyzz = cbuffer.data(gh_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_xyzz_xxzzz = cbuffer.data(gh_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyy = cbuffer.data(gh_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyz = cbuffer.data(gh_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyzz = cbuffer.data(gh_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_xyzz_xyzzz = cbuffer.data(gh_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_xyzz_xzzzz = cbuffer.data(gh_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyy = cbuffer.data(gh_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyz = cbuffer.data(gh_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyzz = cbuffer.data(gh_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_xyzz_yyzzz = cbuffer.data(gh_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_xyzz_yzzzz = cbuffer.data(gh_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_xyzz_zzzzz = cbuffer.data(gh_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxx = cbuffer.data(gh_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxy = cbuffer.data(gh_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxz = cbuffer.data(gh_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyy = cbuffer.data(gh_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyz = cbuffer.data(gh_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxzz = cbuffer.data(gh_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyy = cbuffer.data(gh_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyz = cbuffer.data(gh_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyzz = cbuffer.data(gh_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_xzzz_xxzzz = cbuffer.data(gh_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyy = cbuffer.data(gh_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyz = cbuffer.data(gh_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyzz = cbuffer.data(gh_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_xzzz_xyzzz = cbuffer.data(gh_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_xzzz_xzzzz = cbuffer.data(gh_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyy = cbuffer.data(gh_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyz = cbuffer.data(gh_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyzz = cbuffer.data(gh_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_xzzz_yyzzz = cbuffer.data(gh_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_xzzz_yzzzz = cbuffer.data(gh_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_xzzz_zzzzz = cbuffer.data(gh_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxx = cbuffer.data(gh_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxy = cbuffer.data(gh_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxz = cbuffer.data(gh_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyy = cbuffer.data(gh_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyz = cbuffer.data(gh_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxzz = cbuffer.data(gh_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyy = cbuffer.data(gh_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyz = cbuffer.data(gh_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyzz = cbuffer.data(gh_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_yyyy_xxzzz = cbuffer.data(gh_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyy = cbuffer.data(gh_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyz = cbuffer.data(gh_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyzz = cbuffer.data(gh_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_yyyy_xyzzz = cbuffer.data(gh_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_yyyy_xzzzz = cbuffer.data(gh_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyy = cbuffer.data(gh_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyz = cbuffer.data(gh_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyzz = cbuffer.data(gh_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_yyyy_yyzzz = cbuffer.data(gh_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_yyyy_yzzzz = cbuffer.data(gh_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_yyyy_zzzzz = cbuffer.data(gh_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxx = cbuffer.data(gh_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxy = cbuffer.data(gh_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxz = cbuffer.data(gh_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyy = cbuffer.data(gh_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyz = cbuffer.data(gh_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxzz = cbuffer.data(gh_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyy = cbuffer.data(gh_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyz = cbuffer.data(gh_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyzz = cbuffer.data(gh_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_yyyz_xxzzz = cbuffer.data(gh_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyy = cbuffer.data(gh_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyz = cbuffer.data(gh_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyzz = cbuffer.data(gh_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_yyyz_xyzzz = cbuffer.data(gh_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_y_yyyz_xzzzz = cbuffer.data(gh_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyy = cbuffer.data(gh_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyz = cbuffer.data(gh_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyzz = cbuffer.data(gh_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_y_yyyz_yyzzz = cbuffer.data(gh_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_y_yyyz_yzzzz = cbuffer.data(gh_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_y_yyyz_zzzzz = cbuffer.data(gh_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxx = cbuffer.data(gh_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxy = cbuffer.data(gh_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxz = cbuffer.data(gh_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyy = cbuffer.data(gh_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyz = cbuffer.data(gh_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxzz = cbuffer.data(gh_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyy = cbuffer.data(gh_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyz = cbuffer.data(gh_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyzz = cbuffer.data(gh_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_y_yyzz_xxzzz = cbuffer.data(gh_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyy = cbuffer.data(gh_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyz = cbuffer.data(gh_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyzz = cbuffer.data(gh_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_y_yyzz_xyzzz = cbuffer.data(gh_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_y_yyzz_xzzzz = cbuffer.data(gh_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyy = cbuffer.data(gh_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyz = cbuffer.data(gh_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyzz = cbuffer.data(gh_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_y_yyzz_yyzzz = cbuffer.data(gh_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_y_yyzz_yzzzz = cbuffer.data(gh_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_y_yyzz_zzzzz = cbuffer.data(gh_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxx = cbuffer.data(gh_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxy = cbuffer.data(gh_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxz = cbuffer.data(gh_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyy = cbuffer.data(gh_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyz = cbuffer.data(gh_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxzz = cbuffer.data(gh_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyy = cbuffer.data(gh_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyz = cbuffer.data(gh_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyzz = cbuffer.data(gh_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_yzzz_xxzzz = cbuffer.data(gh_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyy = cbuffer.data(gh_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyz = cbuffer.data(gh_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyzz = cbuffer.data(gh_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_yzzz_xyzzz = cbuffer.data(gh_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_yzzz_xzzzz = cbuffer.data(gh_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyy = cbuffer.data(gh_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyz = cbuffer.data(gh_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyzz = cbuffer.data(gh_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_yzzz_yyzzz = cbuffer.data(gh_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_yzzz_yzzzz = cbuffer.data(gh_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_yzzz_zzzzz = cbuffer.data(gh_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxx = cbuffer.data(gh_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxy = cbuffer.data(gh_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxz = cbuffer.data(gh_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyy = cbuffer.data(gh_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyz = cbuffer.data(gh_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxzz = cbuffer.data(gh_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyy = cbuffer.data(gh_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyz = cbuffer.data(gh_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyzz = cbuffer.data(gh_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_zzzz_xxzzz = cbuffer.data(gh_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyy = cbuffer.data(gh_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyz = cbuffer.data(gh_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyzz = cbuffer.data(gh_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_zzzz_xyzzz = cbuffer.data(gh_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_zzzz_xzzzz = cbuffer.data(gh_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyy = cbuffer.data(gh_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyz = cbuffer.data(gh_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyzz = cbuffer.data(gh_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_zzzz_yyzzz = cbuffer.data(gh_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_zzzz_yzzzz = cbuffer.data(gh_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_zzzz_zzzzz = cbuffer.data(gh_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxx = cbuffer.data(gh_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxy = cbuffer.data(gh_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxz = cbuffer.data(gh_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyy = cbuffer.data(gh_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyz = cbuffer.data(gh_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxzz = cbuffer.data(gh_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyy = cbuffer.data(gh_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyz = cbuffer.data(gh_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyzz = cbuffer.data(gh_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_z_xxxx_xxzzz = cbuffer.data(gh_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyy = cbuffer.data(gh_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyz = cbuffer.data(gh_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyzz = cbuffer.data(gh_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_z_xxxx_xyzzz = cbuffer.data(gh_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_z_xxxx_xzzzz = cbuffer.data(gh_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyy = cbuffer.data(gh_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyz = cbuffer.data(gh_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyzz = cbuffer.data(gh_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_z_xxxx_yyzzz = cbuffer.data(gh_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_z_xxxx_yzzzz = cbuffer.data(gh_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_z_xxxx_zzzzz = cbuffer.data(gh_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxx = cbuffer.data(gh_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxy = cbuffer.data(gh_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxz = cbuffer.data(gh_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyy = cbuffer.data(gh_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyz = cbuffer.data(gh_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxzz = cbuffer.data(gh_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyy = cbuffer.data(gh_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyz = cbuffer.data(gh_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyzz = cbuffer.data(gh_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_z_xxxy_xxzzz = cbuffer.data(gh_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyy = cbuffer.data(gh_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyz = cbuffer.data(gh_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyzz = cbuffer.data(gh_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_z_xxxy_xyzzz = cbuffer.data(gh_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_z_xxxy_xzzzz = cbuffer.data(gh_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyy = cbuffer.data(gh_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyz = cbuffer.data(gh_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyzz = cbuffer.data(gh_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_z_xxxy_yyzzz = cbuffer.data(gh_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_z_xxxy_yzzzz = cbuffer.data(gh_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_z_xxxy_zzzzz = cbuffer.data(gh_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxx = cbuffer.data(gh_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxy = cbuffer.data(gh_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxz = cbuffer.data(gh_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyy = cbuffer.data(gh_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyz = cbuffer.data(gh_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxzz = cbuffer.data(gh_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyy = cbuffer.data(gh_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyz = cbuffer.data(gh_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyzz = cbuffer.data(gh_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_z_xxxz_xxzzz = cbuffer.data(gh_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyy = cbuffer.data(gh_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyz = cbuffer.data(gh_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyzz = cbuffer.data(gh_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_z_xxxz_xyzzz = cbuffer.data(gh_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_z_xxxz_xzzzz = cbuffer.data(gh_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyy = cbuffer.data(gh_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyz = cbuffer.data(gh_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyzz = cbuffer.data(gh_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_z_xxxz_yyzzz = cbuffer.data(gh_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_z_xxxz_yzzzz = cbuffer.data(gh_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_z_xxxz_zzzzz = cbuffer.data(gh_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxx = cbuffer.data(gh_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxy = cbuffer.data(gh_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxz = cbuffer.data(gh_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyy = cbuffer.data(gh_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyz = cbuffer.data(gh_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxzz = cbuffer.data(gh_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyy = cbuffer.data(gh_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyz = cbuffer.data(gh_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyzz = cbuffer.data(gh_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_z_xxyy_xxzzz = cbuffer.data(gh_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyy = cbuffer.data(gh_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyz = cbuffer.data(gh_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyzz = cbuffer.data(gh_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_z_xxyy_xyzzz = cbuffer.data(gh_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_z_xxyy_xzzzz = cbuffer.data(gh_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyy = cbuffer.data(gh_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyz = cbuffer.data(gh_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyzz = cbuffer.data(gh_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_z_xxyy_yyzzz = cbuffer.data(gh_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_z_xxyy_yzzzz = cbuffer.data(gh_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_z_xxyy_zzzzz = cbuffer.data(gh_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxx = cbuffer.data(gh_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxy = cbuffer.data(gh_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxz = cbuffer.data(gh_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyy = cbuffer.data(gh_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyz = cbuffer.data(gh_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxzz = cbuffer.data(gh_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyy = cbuffer.data(gh_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyz = cbuffer.data(gh_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyzz = cbuffer.data(gh_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_z_xxyz_xxzzz = cbuffer.data(gh_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyy = cbuffer.data(gh_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyz = cbuffer.data(gh_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyzz = cbuffer.data(gh_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_z_xxyz_xyzzz = cbuffer.data(gh_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_z_xxyz_xzzzz = cbuffer.data(gh_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyy = cbuffer.data(gh_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyz = cbuffer.data(gh_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyzz = cbuffer.data(gh_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_z_xxyz_yyzzz = cbuffer.data(gh_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_z_xxyz_yzzzz = cbuffer.data(gh_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_z_xxyz_zzzzz = cbuffer.data(gh_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxx = cbuffer.data(gh_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxy = cbuffer.data(gh_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxz = cbuffer.data(gh_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyy = cbuffer.data(gh_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyz = cbuffer.data(gh_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxzz = cbuffer.data(gh_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyy = cbuffer.data(gh_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyz = cbuffer.data(gh_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyzz = cbuffer.data(gh_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_z_xxzz_xxzzz = cbuffer.data(gh_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyy = cbuffer.data(gh_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyz = cbuffer.data(gh_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyzz = cbuffer.data(gh_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_z_xxzz_xyzzz = cbuffer.data(gh_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_z_xxzz_xzzzz = cbuffer.data(gh_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyy = cbuffer.data(gh_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyz = cbuffer.data(gh_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyzz = cbuffer.data(gh_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_z_xxzz_yyzzz = cbuffer.data(gh_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_z_xxzz_yzzzz = cbuffer.data(gh_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_z_xxzz_zzzzz = cbuffer.data(gh_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxx = cbuffer.data(gh_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxy = cbuffer.data(gh_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxz = cbuffer.data(gh_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyy = cbuffer.data(gh_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyz = cbuffer.data(gh_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxzz = cbuffer.data(gh_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyy = cbuffer.data(gh_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyz = cbuffer.data(gh_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyzz = cbuffer.data(gh_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_z_xyyy_xxzzz = cbuffer.data(gh_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyy = cbuffer.data(gh_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyz = cbuffer.data(gh_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyzz = cbuffer.data(gh_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_z_xyyy_xyzzz = cbuffer.data(gh_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_z_xyyy_xzzzz = cbuffer.data(gh_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyy = cbuffer.data(gh_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyz = cbuffer.data(gh_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyzz = cbuffer.data(gh_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_z_xyyy_yyzzz = cbuffer.data(gh_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_z_xyyy_yzzzz = cbuffer.data(gh_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_z_xyyy_zzzzz = cbuffer.data(gh_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxx = cbuffer.data(gh_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxy = cbuffer.data(gh_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxz = cbuffer.data(gh_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyy = cbuffer.data(gh_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyz = cbuffer.data(gh_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxzz = cbuffer.data(gh_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyy = cbuffer.data(gh_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyz = cbuffer.data(gh_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyzz = cbuffer.data(gh_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_z_xyyz_xxzzz = cbuffer.data(gh_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyy = cbuffer.data(gh_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyz = cbuffer.data(gh_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyzz = cbuffer.data(gh_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_z_xyyz_xyzzz = cbuffer.data(gh_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_z_xyyz_xzzzz = cbuffer.data(gh_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyy = cbuffer.data(gh_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyz = cbuffer.data(gh_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyzz = cbuffer.data(gh_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_z_xyyz_yyzzz = cbuffer.data(gh_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_z_xyyz_yzzzz = cbuffer.data(gh_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_z_xyyz_zzzzz = cbuffer.data(gh_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxx = cbuffer.data(gh_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxy = cbuffer.data(gh_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxz = cbuffer.data(gh_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyy = cbuffer.data(gh_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyz = cbuffer.data(gh_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxzz = cbuffer.data(gh_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyy = cbuffer.data(gh_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyz = cbuffer.data(gh_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyzz = cbuffer.data(gh_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_z_xyzz_xxzzz = cbuffer.data(gh_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyy = cbuffer.data(gh_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyz = cbuffer.data(gh_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyzz = cbuffer.data(gh_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_z_xyzz_xyzzz = cbuffer.data(gh_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_z_xyzz_xzzzz = cbuffer.data(gh_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyy = cbuffer.data(gh_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyz = cbuffer.data(gh_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyzz = cbuffer.data(gh_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_z_xyzz_yyzzz = cbuffer.data(gh_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_z_xyzz_yzzzz = cbuffer.data(gh_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_z_xyzz_zzzzz = cbuffer.data(gh_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxx = cbuffer.data(gh_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxy = cbuffer.data(gh_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxz = cbuffer.data(gh_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyy = cbuffer.data(gh_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyz = cbuffer.data(gh_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxzz = cbuffer.data(gh_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyy = cbuffer.data(gh_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyz = cbuffer.data(gh_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyzz = cbuffer.data(gh_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_z_xzzz_xxzzz = cbuffer.data(gh_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyy = cbuffer.data(gh_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyz = cbuffer.data(gh_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyzz = cbuffer.data(gh_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_z_xzzz_xyzzz = cbuffer.data(gh_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_z_xzzz_xzzzz = cbuffer.data(gh_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyy = cbuffer.data(gh_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyz = cbuffer.data(gh_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyzz = cbuffer.data(gh_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_z_xzzz_yyzzz = cbuffer.data(gh_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_z_xzzz_yzzzz = cbuffer.data(gh_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_z_xzzz_zzzzz = cbuffer.data(gh_geom_01_off + 839 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxx = cbuffer.data(gh_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxy = cbuffer.data(gh_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxz = cbuffer.data(gh_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyy = cbuffer.data(gh_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyz = cbuffer.data(gh_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxzz = cbuffer.data(gh_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyy = cbuffer.data(gh_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyz = cbuffer.data(gh_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyzz = cbuffer.data(gh_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_z_yyyy_xxzzz = cbuffer.data(gh_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyy = cbuffer.data(gh_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyz = cbuffer.data(gh_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyzz = cbuffer.data(gh_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_z_yyyy_xyzzz = cbuffer.data(gh_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_z_yyyy_xzzzz = cbuffer.data(gh_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyy = cbuffer.data(gh_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyz = cbuffer.data(gh_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyzz = cbuffer.data(gh_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_z_yyyy_yyzzz = cbuffer.data(gh_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_z_yyyy_yzzzz = cbuffer.data(gh_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_z_yyyy_zzzzz = cbuffer.data(gh_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxx = cbuffer.data(gh_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxy = cbuffer.data(gh_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxz = cbuffer.data(gh_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyy = cbuffer.data(gh_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyz = cbuffer.data(gh_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxzz = cbuffer.data(gh_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyy = cbuffer.data(gh_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyz = cbuffer.data(gh_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyzz = cbuffer.data(gh_geom_01_off + 869 * ccomps * dcomps);

            auto g_0_z_yyyz_xxzzz = cbuffer.data(gh_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyy = cbuffer.data(gh_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyz = cbuffer.data(gh_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyzz = cbuffer.data(gh_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_z_yyyz_xyzzz = cbuffer.data(gh_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_z_yyyz_xzzzz = cbuffer.data(gh_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyy = cbuffer.data(gh_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyz = cbuffer.data(gh_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyzz = cbuffer.data(gh_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_z_yyyz_yyzzz = cbuffer.data(gh_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_z_yyyz_yzzzz = cbuffer.data(gh_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_z_yyyz_zzzzz = cbuffer.data(gh_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxx = cbuffer.data(gh_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxy = cbuffer.data(gh_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxz = cbuffer.data(gh_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyy = cbuffer.data(gh_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyz = cbuffer.data(gh_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxzz = cbuffer.data(gh_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyy = cbuffer.data(gh_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyz = cbuffer.data(gh_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyzz = cbuffer.data(gh_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_z_yyzz_xxzzz = cbuffer.data(gh_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyy = cbuffer.data(gh_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyz = cbuffer.data(gh_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyzz = cbuffer.data(gh_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_z_yyzz_xyzzz = cbuffer.data(gh_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_z_yyzz_xzzzz = cbuffer.data(gh_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyy = cbuffer.data(gh_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyz = cbuffer.data(gh_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyzz = cbuffer.data(gh_geom_01_off + 899 * ccomps * dcomps);

            auto g_0_z_yyzz_yyzzz = cbuffer.data(gh_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_z_yyzz_yzzzz = cbuffer.data(gh_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_z_yyzz_zzzzz = cbuffer.data(gh_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxx = cbuffer.data(gh_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxy = cbuffer.data(gh_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxz = cbuffer.data(gh_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyy = cbuffer.data(gh_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyz = cbuffer.data(gh_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxzz = cbuffer.data(gh_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyy = cbuffer.data(gh_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyz = cbuffer.data(gh_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyzz = cbuffer.data(gh_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_z_yzzz_xxzzz = cbuffer.data(gh_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyy = cbuffer.data(gh_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyz = cbuffer.data(gh_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyzz = cbuffer.data(gh_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_z_yzzz_xyzzz = cbuffer.data(gh_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_z_yzzz_xzzzz = cbuffer.data(gh_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyy = cbuffer.data(gh_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyz = cbuffer.data(gh_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyzz = cbuffer.data(gh_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_z_yzzz_yyzzz = cbuffer.data(gh_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_z_yzzz_yzzzz = cbuffer.data(gh_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_z_yzzz_zzzzz = cbuffer.data(gh_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxx = cbuffer.data(gh_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxy = cbuffer.data(gh_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxz = cbuffer.data(gh_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyy = cbuffer.data(gh_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyz = cbuffer.data(gh_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxzz = cbuffer.data(gh_geom_01_off + 929 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyy = cbuffer.data(gh_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyz = cbuffer.data(gh_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyzz = cbuffer.data(gh_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_z_zzzz_xxzzz = cbuffer.data(gh_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyy = cbuffer.data(gh_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyz = cbuffer.data(gh_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyzz = cbuffer.data(gh_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_z_zzzz_xyzzz = cbuffer.data(gh_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_z_zzzz_xzzzz = cbuffer.data(gh_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyy = cbuffer.data(gh_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyz = cbuffer.data(gh_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyzz = cbuffer.data(gh_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_z_zzzz_yyzzz = cbuffer.data(gh_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_z_zzzz_yzzzz = cbuffer.data(gh_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_z_zzzz_zzzzz = cbuffer.data(gh_geom_01_off + 944 * ccomps * dcomps);

            /// Set up components of auxilary buffer : GISS

            const auto gi_geom_01_off = idx_geom_01_gixx + i * dcomps + j;

            auto g_0_x_xxxx_xxxxxx = cbuffer.data(gi_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxy = cbuffer.data(gi_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxxz = cbuffer.data(gi_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxyy = cbuffer.data(gi_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxyz = cbuffer.data(gi_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxxzz = cbuffer.data(gi_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyyy = cbuffer.data(gi_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyyz = cbuffer.data(gi_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxyzz = cbuffer.data(gi_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxx_xxxzzz = cbuffer.data(gi_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyyy = cbuffer.data(gi_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyyz = cbuffer.data(gi_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyyzz = cbuffer.data(gi_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxx_xxyzzz = cbuffer.data(gi_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxx_xxzzzz = cbuffer.data(gi_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyyy = cbuffer.data(gi_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyyz = cbuffer.data(gi_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyyzz = cbuffer.data(gi_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxx_xyyzzz = cbuffer.data(gi_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxx_xyzzzz = cbuffer.data(gi_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxx_xzzzzz = cbuffer.data(gi_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyyy = cbuffer.data(gi_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyyz = cbuffer.data(gi_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyyzz = cbuffer.data(gi_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxx_yyyzzz = cbuffer.data(gi_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxx_yyzzzz = cbuffer.data(gi_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxx_yzzzzz = cbuffer.data(gi_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxx_zzzzzz = cbuffer.data(gi_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxx = cbuffer.data(gi_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxy = cbuffer.data(gi_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxxz = cbuffer.data(gi_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxyy = cbuffer.data(gi_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxyz = cbuffer.data(gi_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxxzz = cbuffer.data(gi_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyyy = cbuffer.data(gi_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyyz = cbuffer.data(gi_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxyzz = cbuffer.data(gi_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxy_xxxzzz = cbuffer.data(gi_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyyy = cbuffer.data(gi_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyyz = cbuffer.data(gi_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyyzz = cbuffer.data(gi_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxy_xxyzzz = cbuffer.data(gi_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxy_xxzzzz = cbuffer.data(gi_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyyy = cbuffer.data(gi_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyyz = cbuffer.data(gi_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyyzz = cbuffer.data(gi_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxy_xyyzzz = cbuffer.data(gi_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxy_xyzzzz = cbuffer.data(gi_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxy_xzzzzz = cbuffer.data(gi_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyyy = cbuffer.data(gi_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyyz = cbuffer.data(gi_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyyzz = cbuffer.data(gi_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxy_yyyzzz = cbuffer.data(gi_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxy_yyzzzz = cbuffer.data(gi_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxy_yzzzzz = cbuffer.data(gi_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxy_zzzzzz = cbuffer.data(gi_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxx = cbuffer.data(gi_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxy = cbuffer.data(gi_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxxz = cbuffer.data(gi_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxyy = cbuffer.data(gi_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxyz = cbuffer.data(gi_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxxzz = cbuffer.data(gi_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyyy = cbuffer.data(gi_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyyz = cbuffer.data(gi_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxyzz = cbuffer.data(gi_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxz_xxxzzz = cbuffer.data(gi_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyyy = cbuffer.data(gi_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyyz = cbuffer.data(gi_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyyzz = cbuffer.data(gi_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxz_xxyzzz = cbuffer.data(gi_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxz_xxzzzz = cbuffer.data(gi_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyyy = cbuffer.data(gi_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyyz = cbuffer.data(gi_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyyzz = cbuffer.data(gi_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxz_xyyzzz = cbuffer.data(gi_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxz_xyzzzz = cbuffer.data(gi_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxz_xzzzzz = cbuffer.data(gi_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyyy = cbuffer.data(gi_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyyz = cbuffer.data(gi_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyyzz = cbuffer.data(gi_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxz_yyyzzz = cbuffer.data(gi_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxz_yyzzzz = cbuffer.data(gi_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxz_yzzzzz = cbuffer.data(gi_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxz_zzzzzz = cbuffer.data(gi_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxx = cbuffer.data(gi_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxy = cbuffer.data(gi_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxxz = cbuffer.data(gi_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxyy = cbuffer.data(gi_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxyz = cbuffer.data(gi_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxxzz = cbuffer.data(gi_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyyy = cbuffer.data(gi_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyyz = cbuffer.data(gi_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxyzz = cbuffer.data(gi_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxyy_xxxzzz = cbuffer.data(gi_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyyy = cbuffer.data(gi_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyyz = cbuffer.data(gi_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyyzz = cbuffer.data(gi_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxyy_xxyzzz = cbuffer.data(gi_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxyy_xxzzzz = cbuffer.data(gi_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyyy = cbuffer.data(gi_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyyz = cbuffer.data(gi_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyyzz = cbuffer.data(gi_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxyy_xyyzzz = cbuffer.data(gi_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxyy_xyzzzz = cbuffer.data(gi_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxyy_xzzzzz = cbuffer.data(gi_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyyy = cbuffer.data(gi_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyyz = cbuffer.data(gi_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyyzz = cbuffer.data(gi_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxyy_yyyzzz = cbuffer.data(gi_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxyy_yyzzzz = cbuffer.data(gi_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxyy_yzzzzz = cbuffer.data(gi_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxyy_zzzzzz = cbuffer.data(gi_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxx = cbuffer.data(gi_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxy = cbuffer.data(gi_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxxz = cbuffer.data(gi_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxyy = cbuffer.data(gi_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxyz = cbuffer.data(gi_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxxzz = cbuffer.data(gi_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyyy = cbuffer.data(gi_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyyz = cbuffer.data(gi_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxyzz = cbuffer.data(gi_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxyz_xxxzzz = cbuffer.data(gi_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyyy = cbuffer.data(gi_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyyz = cbuffer.data(gi_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyyzz = cbuffer.data(gi_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxyz_xxyzzz = cbuffer.data(gi_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxyz_xxzzzz = cbuffer.data(gi_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyyy = cbuffer.data(gi_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyyz = cbuffer.data(gi_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyyzz = cbuffer.data(gi_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxyz_xyyzzz = cbuffer.data(gi_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxyz_xyzzzz = cbuffer.data(gi_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxyz_xzzzzz = cbuffer.data(gi_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyyy = cbuffer.data(gi_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyyz = cbuffer.data(gi_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyyzz = cbuffer.data(gi_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxyz_yyyzzz = cbuffer.data(gi_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxyz_yyzzzz = cbuffer.data(gi_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxyz_yzzzzz = cbuffer.data(gi_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxyz_zzzzzz = cbuffer.data(gi_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxx = cbuffer.data(gi_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxy = cbuffer.data(gi_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxxz = cbuffer.data(gi_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxyy = cbuffer.data(gi_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxyz = cbuffer.data(gi_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxxzz = cbuffer.data(gi_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyyy = cbuffer.data(gi_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyyz = cbuffer.data(gi_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxyzz = cbuffer.data(gi_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxzz_xxxzzz = cbuffer.data(gi_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyyy = cbuffer.data(gi_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyyz = cbuffer.data(gi_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyyzz = cbuffer.data(gi_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xxzz_xxyzzz = cbuffer.data(gi_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xxzz_xxzzzz = cbuffer.data(gi_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyyy = cbuffer.data(gi_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyyz = cbuffer.data(gi_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyyzz = cbuffer.data(gi_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xxzz_xyyzzz = cbuffer.data(gi_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xxzz_xyzzzz = cbuffer.data(gi_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xxzz_xzzzzz = cbuffer.data(gi_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyyy = cbuffer.data(gi_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyyz = cbuffer.data(gi_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyyzz = cbuffer.data(gi_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xxzz_yyyzzz = cbuffer.data(gi_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xxzz_yyzzzz = cbuffer.data(gi_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xxzz_yzzzzz = cbuffer.data(gi_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xxzz_zzzzzz = cbuffer.data(gi_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxx = cbuffer.data(gi_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxy = cbuffer.data(gi_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxxz = cbuffer.data(gi_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxyy = cbuffer.data(gi_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxyz = cbuffer.data(gi_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxxzz = cbuffer.data(gi_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyyy = cbuffer.data(gi_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyyz = cbuffer.data(gi_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxyzz = cbuffer.data(gi_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xyyy_xxxzzz = cbuffer.data(gi_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyyy = cbuffer.data(gi_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyyz = cbuffer.data(gi_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyyzz = cbuffer.data(gi_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xyyy_xxyzzz = cbuffer.data(gi_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xyyy_xxzzzz = cbuffer.data(gi_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyyy = cbuffer.data(gi_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyyz = cbuffer.data(gi_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyyzz = cbuffer.data(gi_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xyyy_xyyzzz = cbuffer.data(gi_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xyyy_xyzzzz = cbuffer.data(gi_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xyyy_xzzzzz = cbuffer.data(gi_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyyy = cbuffer.data(gi_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyyz = cbuffer.data(gi_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyyzz = cbuffer.data(gi_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xyyy_yyyzzz = cbuffer.data(gi_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xyyy_yyzzzz = cbuffer.data(gi_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xyyy_yzzzzz = cbuffer.data(gi_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xyyy_zzzzzz = cbuffer.data(gi_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxx = cbuffer.data(gi_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxy = cbuffer.data(gi_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxxz = cbuffer.data(gi_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxyy = cbuffer.data(gi_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxyz = cbuffer.data(gi_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxxzz = cbuffer.data(gi_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyyy = cbuffer.data(gi_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyyz = cbuffer.data(gi_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxyzz = cbuffer.data(gi_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xyyz_xxxzzz = cbuffer.data(gi_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyyy = cbuffer.data(gi_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyyz = cbuffer.data(gi_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyyzz = cbuffer.data(gi_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xyyz_xxyzzz = cbuffer.data(gi_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_xyyz_xxzzzz = cbuffer.data(gi_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyyy = cbuffer.data(gi_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyyz = cbuffer.data(gi_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyyzz = cbuffer.data(gi_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xyyz_xyyzzz = cbuffer.data(gi_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xyyz_xyzzzz = cbuffer.data(gi_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_xyyz_xzzzzz = cbuffer.data(gi_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyyy = cbuffer.data(gi_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyyz = cbuffer.data(gi_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyyzz = cbuffer.data(gi_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_xyyz_yyyzzz = cbuffer.data(gi_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xyyz_yyzzzz = cbuffer.data(gi_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xyyz_yzzzzz = cbuffer.data(gi_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xyyz_zzzzzz = cbuffer.data(gi_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxx = cbuffer.data(gi_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxy = cbuffer.data(gi_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxxz = cbuffer.data(gi_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxyy = cbuffer.data(gi_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxyz = cbuffer.data(gi_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxxzz = cbuffer.data(gi_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyyy = cbuffer.data(gi_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyyz = cbuffer.data(gi_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxyzz = cbuffer.data(gi_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_xyzz_xxxzzz = cbuffer.data(gi_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyyy = cbuffer.data(gi_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyyz = cbuffer.data(gi_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyyzz = cbuffer.data(gi_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_xyzz_xxyzzz = cbuffer.data(gi_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_xyzz_xxzzzz = cbuffer.data(gi_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyyy = cbuffer.data(gi_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyyz = cbuffer.data(gi_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyyzz = cbuffer.data(gi_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_xyzz_xyyzzz = cbuffer.data(gi_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_xyzz_xyzzzz = cbuffer.data(gi_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_xyzz_xzzzzz = cbuffer.data(gi_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyyy = cbuffer.data(gi_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyyz = cbuffer.data(gi_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyyzz = cbuffer.data(gi_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_xyzz_yyyzzz = cbuffer.data(gi_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_xyzz_yyzzzz = cbuffer.data(gi_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_xyzz_yzzzzz = cbuffer.data(gi_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_xyzz_zzzzzz = cbuffer.data(gi_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxx = cbuffer.data(gi_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxy = cbuffer.data(gi_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxxz = cbuffer.data(gi_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxyy = cbuffer.data(gi_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxyz = cbuffer.data(gi_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxxzz = cbuffer.data(gi_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyyy = cbuffer.data(gi_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyyz = cbuffer.data(gi_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxyzz = cbuffer.data(gi_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_xzzz_xxxzzz = cbuffer.data(gi_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyyy = cbuffer.data(gi_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyyz = cbuffer.data(gi_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyyzz = cbuffer.data(gi_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_xzzz_xxyzzz = cbuffer.data(gi_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_xzzz_xxzzzz = cbuffer.data(gi_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyyy = cbuffer.data(gi_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyyz = cbuffer.data(gi_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyyzz = cbuffer.data(gi_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_xzzz_xyyzzz = cbuffer.data(gi_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_xzzz_xyzzzz = cbuffer.data(gi_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_xzzz_xzzzzz = cbuffer.data(gi_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyyy = cbuffer.data(gi_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyyz = cbuffer.data(gi_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyyzz = cbuffer.data(gi_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_xzzz_yyyzzz = cbuffer.data(gi_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_xzzz_yyzzzz = cbuffer.data(gi_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_xzzz_yzzzzz = cbuffer.data(gi_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_xzzz_zzzzzz = cbuffer.data(gi_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxx = cbuffer.data(gi_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxy = cbuffer.data(gi_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxxz = cbuffer.data(gi_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxyy = cbuffer.data(gi_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxyz = cbuffer.data(gi_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxxzz = cbuffer.data(gi_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyyy = cbuffer.data(gi_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyyz = cbuffer.data(gi_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxyzz = cbuffer.data(gi_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_yyyy_xxxzzz = cbuffer.data(gi_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyyy = cbuffer.data(gi_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyyz = cbuffer.data(gi_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyyzz = cbuffer.data(gi_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_yyyy_xxyzzz = cbuffer.data(gi_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_yyyy_xxzzzz = cbuffer.data(gi_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyyy = cbuffer.data(gi_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyyz = cbuffer.data(gi_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyyzz = cbuffer.data(gi_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_yyyy_xyyzzz = cbuffer.data(gi_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_yyyy_xyzzzz = cbuffer.data(gi_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_x_yyyy_xzzzzz = cbuffer.data(gi_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyyy = cbuffer.data(gi_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyyz = cbuffer.data(gi_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyyzz = cbuffer.data(gi_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_yyyy_yyyzzz = cbuffer.data(gi_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_yyyy_yyzzzz = cbuffer.data(gi_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_yyyy_yzzzzz = cbuffer.data(gi_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_yyyy_zzzzzz = cbuffer.data(gi_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxx = cbuffer.data(gi_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxy = cbuffer.data(gi_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxxz = cbuffer.data(gi_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxyy = cbuffer.data(gi_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxyz = cbuffer.data(gi_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxxzz = cbuffer.data(gi_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyyy = cbuffer.data(gi_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyyz = cbuffer.data(gi_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxyzz = cbuffer.data(gi_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_yyyz_xxxzzz = cbuffer.data(gi_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyyy = cbuffer.data(gi_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyyz = cbuffer.data(gi_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyyzz = cbuffer.data(gi_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_yyyz_xxyzzz = cbuffer.data(gi_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_yyyz_xxzzzz = cbuffer.data(gi_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyyy = cbuffer.data(gi_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyyz = cbuffer.data(gi_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyyzz = cbuffer.data(gi_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_yyyz_xyyzzz = cbuffer.data(gi_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_yyyz_xyzzzz = cbuffer.data(gi_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_yyyz_xzzzzz = cbuffer.data(gi_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyyy = cbuffer.data(gi_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyyz = cbuffer.data(gi_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyyzz = cbuffer.data(gi_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_yyyz_yyyzzz = cbuffer.data(gi_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_yyyz_yyzzzz = cbuffer.data(gi_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_yyyz_yzzzzz = cbuffer.data(gi_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_yyyz_zzzzzz = cbuffer.data(gi_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxx = cbuffer.data(gi_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxy = cbuffer.data(gi_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxxz = cbuffer.data(gi_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxyy = cbuffer.data(gi_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxyz = cbuffer.data(gi_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxxzz = cbuffer.data(gi_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyyy = cbuffer.data(gi_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyyz = cbuffer.data(gi_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxyzz = cbuffer.data(gi_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_x_yyzz_xxxzzz = cbuffer.data(gi_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyyy = cbuffer.data(gi_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyyz = cbuffer.data(gi_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyyzz = cbuffer.data(gi_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_yyzz_xxyzzz = cbuffer.data(gi_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_x_yyzz_xxzzzz = cbuffer.data(gi_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyyy = cbuffer.data(gi_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyyz = cbuffer.data(gi_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyyzz = cbuffer.data(gi_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_yyzz_xyyzzz = cbuffer.data(gi_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_yyzz_xyzzzz = cbuffer.data(gi_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_yyzz_xzzzzz = cbuffer.data(gi_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyyy = cbuffer.data(gi_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyyz = cbuffer.data(gi_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyyzz = cbuffer.data(gi_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_x_yyzz_yyyzzz = cbuffer.data(gi_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_x_yyzz_yyzzzz = cbuffer.data(gi_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_x_yyzz_yzzzzz = cbuffer.data(gi_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_x_yyzz_zzzzzz = cbuffer.data(gi_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxx = cbuffer.data(gi_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxy = cbuffer.data(gi_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxxz = cbuffer.data(gi_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxyy = cbuffer.data(gi_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxyz = cbuffer.data(gi_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxxzz = cbuffer.data(gi_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyyy = cbuffer.data(gi_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyyz = cbuffer.data(gi_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxyzz = cbuffer.data(gi_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_x_yzzz_xxxzzz = cbuffer.data(gi_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyyy = cbuffer.data(gi_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyyz = cbuffer.data(gi_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyyzz = cbuffer.data(gi_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_x_yzzz_xxyzzz = cbuffer.data(gi_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_x_yzzz_xxzzzz = cbuffer.data(gi_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyyy = cbuffer.data(gi_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyyz = cbuffer.data(gi_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyyzz = cbuffer.data(gi_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_x_yzzz_xyyzzz = cbuffer.data(gi_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_x_yzzz_xyzzzz = cbuffer.data(gi_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_x_yzzz_xzzzzz = cbuffer.data(gi_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyyy = cbuffer.data(gi_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyyz = cbuffer.data(gi_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyyzz = cbuffer.data(gi_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_x_yzzz_yyyzzz = cbuffer.data(gi_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_x_yzzz_yyzzzz = cbuffer.data(gi_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_x_yzzz_yzzzzz = cbuffer.data(gi_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_x_yzzz_zzzzzz = cbuffer.data(gi_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxx = cbuffer.data(gi_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxy = cbuffer.data(gi_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxxz = cbuffer.data(gi_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxyy = cbuffer.data(gi_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxyz = cbuffer.data(gi_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxxzz = cbuffer.data(gi_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyyy = cbuffer.data(gi_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyyz = cbuffer.data(gi_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxyzz = cbuffer.data(gi_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_x_zzzz_xxxzzz = cbuffer.data(gi_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyyy = cbuffer.data(gi_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyyz = cbuffer.data(gi_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyyzz = cbuffer.data(gi_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_x_zzzz_xxyzzz = cbuffer.data(gi_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_x_zzzz_xxzzzz = cbuffer.data(gi_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyyy = cbuffer.data(gi_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyyz = cbuffer.data(gi_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyyzz = cbuffer.data(gi_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_x_zzzz_xyyzzz = cbuffer.data(gi_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_x_zzzz_xyzzzz = cbuffer.data(gi_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_x_zzzz_xzzzzz = cbuffer.data(gi_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyyy = cbuffer.data(gi_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyyz = cbuffer.data(gi_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyyzz = cbuffer.data(gi_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_x_zzzz_yyyzzz = cbuffer.data(gi_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_x_zzzz_yyzzzz = cbuffer.data(gi_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_x_zzzz_yzzzzz = cbuffer.data(gi_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_x_zzzz_zzzzzz = cbuffer.data(gi_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxx = cbuffer.data(gi_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxy = cbuffer.data(gi_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxxz = cbuffer.data(gi_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxyy = cbuffer.data(gi_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxyz = cbuffer.data(gi_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxxzz = cbuffer.data(gi_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyyy = cbuffer.data(gi_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyyz = cbuffer.data(gi_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxyzz = cbuffer.data(gi_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_xxxx_xxxzzz = cbuffer.data(gi_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyyy = cbuffer.data(gi_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyyz = cbuffer.data(gi_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyyzz = cbuffer.data(gi_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_xxxx_xxyzzz = cbuffer.data(gi_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_xxxx_xxzzzz = cbuffer.data(gi_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyyy = cbuffer.data(gi_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyyz = cbuffer.data(gi_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyyzz = cbuffer.data(gi_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_xxxx_xyyzzz = cbuffer.data(gi_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_xxxx_xyzzzz = cbuffer.data(gi_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_y_xxxx_xzzzzz = cbuffer.data(gi_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyyy = cbuffer.data(gi_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyyz = cbuffer.data(gi_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyyzz = cbuffer.data(gi_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_xxxx_yyyzzz = cbuffer.data(gi_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_xxxx_yyzzzz = cbuffer.data(gi_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_xxxx_yzzzzz = cbuffer.data(gi_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_xxxx_zzzzzz = cbuffer.data(gi_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxx = cbuffer.data(gi_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxy = cbuffer.data(gi_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxxz = cbuffer.data(gi_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxyy = cbuffer.data(gi_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxyz = cbuffer.data(gi_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxxzz = cbuffer.data(gi_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyyy = cbuffer.data(gi_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyyz = cbuffer.data(gi_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxyzz = cbuffer.data(gi_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_xxxy_xxxzzz = cbuffer.data(gi_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyyy = cbuffer.data(gi_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyyz = cbuffer.data(gi_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyyzz = cbuffer.data(gi_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_xxxy_xxyzzz = cbuffer.data(gi_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_xxxy_xxzzzz = cbuffer.data(gi_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyyy = cbuffer.data(gi_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyyz = cbuffer.data(gi_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyyzz = cbuffer.data(gi_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_xxxy_xyyzzz = cbuffer.data(gi_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_xxxy_xyzzzz = cbuffer.data(gi_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_xxxy_xzzzzz = cbuffer.data(gi_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyyy = cbuffer.data(gi_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyyz = cbuffer.data(gi_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyyzz = cbuffer.data(gi_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_xxxy_yyyzzz = cbuffer.data(gi_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_xxxy_yyzzzz = cbuffer.data(gi_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_xxxy_yzzzzz = cbuffer.data(gi_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_xxxy_zzzzzz = cbuffer.data(gi_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxx = cbuffer.data(gi_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxy = cbuffer.data(gi_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxxz = cbuffer.data(gi_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxyy = cbuffer.data(gi_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxyz = cbuffer.data(gi_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxxzz = cbuffer.data(gi_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyyy = cbuffer.data(gi_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyyz = cbuffer.data(gi_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxyzz = cbuffer.data(gi_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_xxxz_xxxzzz = cbuffer.data(gi_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyyy = cbuffer.data(gi_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyyz = cbuffer.data(gi_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyyzz = cbuffer.data(gi_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_xxxz_xxyzzz = cbuffer.data(gi_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_xxxz_xxzzzz = cbuffer.data(gi_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyyy = cbuffer.data(gi_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyyz = cbuffer.data(gi_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyyzz = cbuffer.data(gi_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_xxxz_xyyzzz = cbuffer.data(gi_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_xxxz_xyzzzz = cbuffer.data(gi_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_xxxz_xzzzzz = cbuffer.data(gi_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyyy = cbuffer.data(gi_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyyz = cbuffer.data(gi_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyyzz = cbuffer.data(gi_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_xxxz_yyyzzz = cbuffer.data(gi_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_xxxz_yyzzzz = cbuffer.data(gi_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_xxxz_yzzzzz = cbuffer.data(gi_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_xxxz_zzzzzz = cbuffer.data(gi_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxx = cbuffer.data(gi_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxy = cbuffer.data(gi_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxxz = cbuffer.data(gi_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxyy = cbuffer.data(gi_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxyz = cbuffer.data(gi_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxxzz = cbuffer.data(gi_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyyy = cbuffer.data(gi_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyyz = cbuffer.data(gi_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxyzz = cbuffer.data(gi_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_xxyy_xxxzzz = cbuffer.data(gi_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyyy = cbuffer.data(gi_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyyz = cbuffer.data(gi_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyyzz = cbuffer.data(gi_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_xxyy_xxyzzz = cbuffer.data(gi_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_xxyy_xxzzzz = cbuffer.data(gi_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyyy = cbuffer.data(gi_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyyz = cbuffer.data(gi_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyyzz = cbuffer.data(gi_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_xxyy_xyyzzz = cbuffer.data(gi_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_xxyy_xyzzzz = cbuffer.data(gi_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_xxyy_xzzzzz = cbuffer.data(gi_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyyy = cbuffer.data(gi_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyyz = cbuffer.data(gi_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyyzz = cbuffer.data(gi_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_xxyy_yyyzzz = cbuffer.data(gi_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_xxyy_yyzzzz = cbuffer.data(gi_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_xxyy_yzzzzz = cbuffer.data(gi_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_xxyy_zzzzzz = cbuffer.data(gi_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxx = cbuffer.data(gi_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxy = cbuffer.data(gi_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxxz = cbuffer.data(gi_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxyy = cbuffer.data(gi_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxyz = cbuffer.data(gi_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxxzz = cbuffer.data(gi_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyyy = cbuffer.data(gi_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyyz = cbuffer.data(gi_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxyzz = cbuffer.data(gi_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_xxyz_xxxzzz = cbuffer.data(gi_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyyy = cbuffer.data(gi_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyyz = cbuffer.data(gi_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyyzz = cbuffer.data(gi_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_xxyz_xxyzzz = cbuffer.data(gi_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_xxyz_xxzzzz = cbuffer.data(gi_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyyy = cbuffer.data(gi_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyyz = cbuffer.data(gi_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyyzz = cbuffer.data(gi_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_xxyz_xyyzzz = cbuffer.data(gi_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_xxyz_xyzzzz = cbuffer.data(gi_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_xxyz_xzzzzz = cbuffer.data(gi_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyyy = cbuffer.data(gi_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyyz = cbuffer.data(gi_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyyzz = cbuffer.data(gi_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_xxyz_yyyzzz = cbuffer.data(gi_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_xxyz_yyzzzz = cbuffer.data(gi_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_xxyz_yzzzzz = cbuffer.data(gi_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_xxyz_zzzzzz = cbuffer.data(gi_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxx = cbuffer.data(gi_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxy = cbuffer.data(gi_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxxz = cbuffer.data(gi_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxyy = cbuffer.data(gi_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxyz = cbuffer.data(gi_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxxzz = cbuffer.data(gi_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyyy = cbuffer.data(gi_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyyz = cbuffer.data(gi_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxyzz = cbuffer.data(gi_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_y_xxzz_xxxzzz = cbuffer.data(gi_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyyy = cbuffer.data(gi_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyyz = cbuffer.data(gi_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyyzz = cbuffer.data(gi_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_y_xxzz_xxyzzz = cbuffer.data(gi_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_y_xxzz_xxzzzz = cbuffer.data(gi_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyyy = cbuffer.data(gi_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyyz = cbuffer.data(gi_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyyzz = cbuffer.data(gi_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_y_xxzz_xyyzzz = cbuffer.data(gi_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_y_xxzz_xyzzzz = cbuffer.data(gi_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_y_xxzz_xzzzzz = cbuffer.data(gi_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyyy = cbuffer.data(gi_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyyz = cbuffer.data(gi_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyyzz = cbuffer.data(gi_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_y_xxzz_yyyzzz = cbuffer.data(gi_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_y_xxzz_yyzzzz = cbuffer.data(gi_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_y_xxzz_yzzzzz = cbuffer.data(gi_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_y_xxzz_zzzzzz = cbuffer.data(gi_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxx = cbuffer.data(gi_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxy = cbuffer.data(gi_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxxz = cbuffer.data(gi_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxyy = cbuffer.data(gi_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxyz = cbuffer.data(gi_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxxzz = cbuffer.data(gi_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyyy = cbuffer.data(gi_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyyz = cbuffer.data(gi_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxyzz = cbuffer.data(gi_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_xyyy_xxxzzz = cbuffer.data(gi_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyyy = cbuffer.data(gi_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyyz = cbuffer.data(gi_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyyzz = cbuffer.data(gi_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_xyyy_xxyzzz = cbuffer.data(gi_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_xyyy_xxzzzz = cbuffer.data(gi_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyyy = cbuffer.data(gi_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyyz = cbuffer.data(gi_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyyzz = cbuffer.data(gi_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_xyyy_xyyzzz = cbuffer.data(gi_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_xyyy_xyzzzz = cbuffer.data(gi_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_xyyy_xzzzzz = cbuffer.data(gi_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyyy = cbuffer.data(gi_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyyz = cbuffer.data(gi_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyyzz = cbuffer.data(gi_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_xyyy_yyyzzz = cbuffer.data(gi_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_xyyy_yyzzzz = cbuffer.data(gi_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_xyyy_yzzzzz = cbuffer.data(gi_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_y_xyyy_zzzzzz = cbuffer.data(gi_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxx = cbuffer.data(gi_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxy = cbuffer.data(gi_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxxz = cbuffer.data(gi_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxyy = cbuffer.data(gi_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxyz = cbuffer.data(gi_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxxzz = cbuffer.data(gi_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyyy = cbuffer.data(gi_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyyz = cbuffer.data(gi_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxyzz = cbuffer.data(gi_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_xyyz_xxxzzz = cbuffer.data(gi_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyyy = cbuffer.data(gi_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyyz = cbuffer.data(gi_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyyzz = cbuffer.data(gi_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_xyyz_xxyzzz = cbuffer.data(gi_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_y_xyyz_xxzzzz = cbuffer.data(gi_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyyy = cbuffer.data(gi_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyyz = cbuffer.data(gi_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyyzz = cbuffer.data(gi_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_y_xyyz_xyyzzz = cbuffer.data(gi_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_y_xyyz_xyzzzz = cbuffer.data(gi_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_y_xyyz_xzzzzz = cbuffer.data(gi_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyyy = cbuffer.data(gi_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyyz = cbuffer.data(gi_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyyzz = cbuffer.data(gi_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_y_xyyz_yyyzzz = cbuffer.data(gi_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_y_xyyz_yyzzzz = cbuffer.data(gi_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_y_xyyz_yzzzzz = cbuffer.data(gi_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_y_xyyz_zzzzzz = cbuffer.data(gi_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxx = cbuffer.data(gi_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxy = cbuffer.data(gi_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxxz = cbuffer.data(gi_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxyy = cbuffer.data(gi_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxyz = cbuffer.data(gi_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxxzz = cbuffer.data(gi_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyyy = cbuffer.data(gi_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyyz = cbuffer.data(gi_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxyzz = cbuffer.data(gi_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_y_xyzz_xxxzzz = cbuffer.data(gi_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyyy = cbuffer.data(gi_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyyz = cbuffer.data(gi_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyyzz = cbuffer.data(gi_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_y_xyzz_xxyzzz = cbuffer.data(gi_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_y_xyzz_xxzzzz = cbuffer.data(gi_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyyy = cbuffer.data(gi_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyyz = cbuffer.data(gi_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyyzz = cbuffer.data(gi_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_y_xyzz_xyyzzz = cbuffer.data(gi_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_y_xyzz_xyzzzz = cbuffer.data(gi_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_y_xyzz_xzzzzz = cbuffer.data(gi_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyyy = cbuffer.data(gi_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyyz = cbuffer.data(gi_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyyzz = cbuffer.data(gi_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_y_xyzz_yyyzzz = cbuffer.data(gi_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_y_xyzz_yyzzzz = cbuffer.data(gi_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_y_xyzz_yzzzzz = cbuffer.data(gi_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_y_xyzz_zzzzzz = cbuffer.data(gi_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxx = cbuffer.data(gi_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxy = cbuffer.data(gi_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxxz = cbuffer.data(gi_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxyy = cbuffer.data(gi_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxyz = cbuffer.data(gi_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxxzz = cbuffer.data(gi_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyyy = cbuffer.data(gi_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyyz = cbuffer.data(gi_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxyzz = cbuffer.data(gi_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_y_xzzz_xxxzzz = cbuffer.data(gi_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyyy = cbuffer.data(gi_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyyz = cbuffer.data(gi_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyyzz = cbuffer.data(gi_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_y_xzzz_xxyzzz = cbuffer.data(gi_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_y_xzzz_xxzzzz = cbuffer.data(gi_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyyy = cbuffer.data(gi_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyyz = cbuffer.data(gi_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyyzz = cbuffer.data(gi_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_y_xzzz_xyyzzz = cbuffer.data(gi_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_y_xzzz_xyzzzz = cbuffer.data(gi_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_y_xzzz_xzzzzz = cbuffer.data(gi_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyyy = cbuffer.data(gi_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyyz = cbuffer.data(gi_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyyzz = cbuffer.data(gi_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_y_xzzz_yyyzzz = cbuffer.data(gi_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_y_xzzz_yyzzzz = cbuffer.data(gi_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_y_xzzz_yzzzzz = cbuffer.data(gi_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_y_xzzz_zzzzzz = cbuffer.data(gi_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxx = cbuffer.data(gi_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxy = cbuffer.data(gi_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxxz = cbuffer.data(gi_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxyy = cbuffer.data(gi_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxyz = cbuffer.data(gi_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxxzz = cbuffer.data(gi_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyyy = cbuffer.data(gi_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyyz = cbuffer.data(gi_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxyzz = cbuffer.data(gi_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_y_yyyy_xxxzzz = cbuffer.data(gi_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyyy = cbuffer.data(gi_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyyz = cbuffer.data(gi_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyyzz = cbuffer.data(gi_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_y_yyyy_xxyzzz = cbuffer.data(gi_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_y_yyyy_xxzzzz = cbuffer.data(gi_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyyy = cbuffer.data(gi_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyyz = cbuffer.data(gi_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyyzz = cbuffer.data(gi_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_y_yyyy_xyyzzz = cbuffer.data(gi_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_y_yyyy_xyzzzz = cbuffer.data(gi_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_y_yyyy_xzzzzz = cbuffer.data(gi_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyyy = cbuffer.data(gi_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyyz = cbuffer.data(gi_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyyzz = cbuffer.data(gi_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_y_yyyy_yyyzzz = cbuffer.data(gi_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_y_yyyy_yyzzzz = cbuffer.data(gi_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_y_yyyy_yzzzzz = cbuffer.data(gi_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_y_yyyy_zzzzzz = cbuffer.data(gi_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxx = cbuffer.data(gi_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxy = cbuffer.data(gi_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxxz = cbuffer.data(gi_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxyy = cbuffer.data(gi_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxyz = cbuffer.data(gi_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxxzz = cbuffer.data(gi_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyyy = cbuffer.data(gi_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyyz = cbuffer.data(gi_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxyzz = cbuffer.data(gi_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_y_yyyz_xxxzzz = cbuffer.data(gi_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyyy = cbuffer.data(gi_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyyz = cbuffer.data(gi_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyyzz = cbuffer.data(gi_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_y_yyyz_xxyzzz = cbuffer.data(gi_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_y_yyyz_xxzzzz = cbuffer.data(gi_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyyy = cbuffer.data(gi_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyyz = cbuffer.data(gi_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyyzz = cbuffer.data(gi_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_y_yyyz_xyyzzz = cbuffer.data(gi_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_y_yyyz_xyzzzz = cbuffer.data(gi_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_y_yyyz_xzzzzz = cbuffer.data(gi_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyyy = cbuffer.data(gi_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyyz = cbuffer.data(gi_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyyzz = cbuffer.data(gi_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_y_yyyz_yyyzzz = cbuffer.data(gi_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_y_yyyz_yyzzzz = cbuffer.data(gi_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_y_yyyz_yzzzzz = cbuffer.data(gi_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_y_yyyz_zzzzzz = cbuffer.data(gi_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxx = cbuffer.data(gi_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxy = cbuffer.data(gi_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxxz = cbuffer.data(gi_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxyy = cbuffer.data(gi_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxyz = cbuffer.data(gi_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxxzz = cbuffer.data(gi_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyyy = cbuffer.data(gi_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyyz = cbuffer.data(gi_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxyzz = cbuffer.data(gi_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_y_yyzz_xxxzzz = cbuffer.data(gi_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyyy = cbuffer.data(gi_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyyz = cbuffer.data(gi_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyyzz = cbuffer.data(gi_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_y_yyzz_xxyzzz = cbuffer.data(gi_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_y_yyzz_xxzzzz = cbuffer.data(gi_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyyy = cbuffer.data(gi_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyyz = cbuffer.data(gi_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyyzz = cbuffer.data(gi_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_y_yyzz_xyyzzz = cbuffer.data(gi_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_y_yyzz_xyzzzz = cbuffer.data(gi_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_y_yyzz_xzzzzz = cbuffer.data(gi_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyyy = cbuffer.data(gi_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyyz = cbuffer.data(gi_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyyzz = cbuffer.data(gi_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_y_yyzz_yyyzzz = cbuffer.data(gi_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_y_yyzz_yyzzzz = cbuffer.data(gi_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_y_yyzz_yzzzzz = cbuffer.data(gi_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_y_yyzz_zzzzzz = cbuffer.data(gi_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxx = cbuffer.data(gi_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxy = cbuffer.data(gi_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxxz = cbuffer.data(gi_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxyy = cbuffer.data(gi_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxyz = cbuffer.data(gi_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxxzz = cbuffer.data(gi_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyyy = cbuffer.data(gi_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyyz = cbuffer.data(gi_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxyzz = cbuffer.data(gi_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_y_yzzz_xxxzzz = cbuffer.data(gi_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyyy = cbuffer.data(gi_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyyz = cbuffer.data(gi_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyyzz = cbuffer.data(gi_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_y_yzzz_xxyzzz = cbuffer.data(gi_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_y_yzzz_xxzzzz = cbuffer.data(gi_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyyy = cbuffer.data(gi_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyyz = cbuffer.data(gi_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyyzz = cbuffer.data(gi_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_y_yzzz_xyyzzz = cbuffer.data(gi_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_y_yzzz_xyzzzz = cbuffer.data(gi_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_y_yzzz_xzzzzz = cbuffer.data(gi_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyyy = cbuffer.data(gi_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyyz = cbuffer.data(gi_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyyzz = cbuffer.data(gi_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_y_yzzz_yyyzzz = cbuffer.data(gi_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_y_yzzz_yyzzzz = cbuffer.data(gi_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_y_yzzz_yzzzzz = cbuffer.data(gi_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_y_yzzz_zzzzzz = cbuffer.data(gi_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxx = cbuffer.data(gi_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxy = cbuffer.data(gi_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxxz = cbuffer.data(gi_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxyy = cbuffer.data(gi_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxyz = cbuffer.data(gi_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxxzz = cbuffer.data(gi_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyyy = cbuffer.data(gi_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyyz = cbuffer.data(gi_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxyzz = cbuffer.data(gi_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_y_zzzz_xxxzzz = cbuffer.data(gi_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyyy = cbuffer.data(gi_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyyz = cbuffer.data(gi_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyyzz = cbuffer.data(gi_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_y_zzzz_xxyzzz = cbuffer.data(gi_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_y_zzzz_xxzzzz = cbuffer.data(gi_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyyy = cbuffer.data(gi_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyyz = cbuffer.data(gi_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyyzz = cbuffer.data(gi_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_y_zzzz_xyyzzz = cbuffer.data(gi_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_y_zzzz_xyzzzz = cbuffer.data(gi_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_y_zzzz_xzzzzz = cbuffer.data(gi_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyyy = cbuffer.data(gi_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyyz = cbuffer.data(gi_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyyzz = cbuffer.data(gi_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_y_zzzz_yyyzzz = cbuffer.data(gi_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_y_zzzz_yyzzzz = cbuffer.data(gi_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_y_zzzz_yzzzzz = cbuffer.data(gi_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_y_zzzz_zzzzzz = cbuffer.data(gi_geom_01_off + 839 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxx = cbuffer.data(gi_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxy = cbuffer.data(gi_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxxz = cbuffer.data(gi_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxyy = cbuffer.data(gi_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxyz = cbuffer.data(gi_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxxzz = cbuffer.data(gi_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyyy = cbuffer.data(gi_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyyz = cbuffer.data(gi_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxyzz = cbuffer.data(gi_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_z_xxxx_xxxzzz = cbuffer.data(gi_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyyy = cbuffer.data(gi_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyyz = cbuffer.data(gi_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyyzz = cbuffer.data(gi_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_z_xxxx_xxyzzz = cbuffer.data(gi_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_z_xxxx_xxzzzz = cbuffer.data(gi_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyyy = cbuffer.data(gi_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyyz = cbuffer.data(gi_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyyzz = cbuffer.data(gi_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_z_xxxx_xyyzzz = cbuffer.data(gi_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_z_xxxx_xyzzzz = cbuffer.data(gi_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_z_xxxx_xzzzzz = cbuffer.data(gi_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyyy = cbuffer.data(gi_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyyz = cbuffer.data(gi_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyyzz = cbuffer.data(gi_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_z_xxxx_yyyzzz = cbuffer.data(gi_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_z_xxxx_yyzzzz = cbuffer.data(gi_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_z_xxxx_yzzzzz = cbuffer.data(gi_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_z_xxxx_zzzzzz = cbuffer.data(gi_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxx = cbuffer.data(gi_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxy = cbuffer.data(gi_geom_01_off + 869 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxxz = cbuffer.data(gi_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxyy = cbuffer.data(gi_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxyz = cbuffer.data(gi_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxxzz = cbuffer.data(gi_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyyy = cbuffer.data(gi_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyyz = cbuffer.data(gi_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxyzz = cbuffer.data(gi_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_z_xxxy_xxxzzz = cbuffer.data(gi_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyyy = cbuffer.data(gi_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyyz = cbuffer.data(gi_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyyzz = cbuffer.data(gi_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_z_xxxy_xxyzzz = cbuffer.data(gi_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_z_xxxy_xxzzzz = cbuffer.data(gi_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyyy = cbuffer.data(gi_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyyz = cbuffer.data(gi_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyyzz = cbuffer.data(gi_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_z_xxxy_xyyzzz = cbuffer.data(gi_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_z_xxxy_xyzzzz = cbuffer.data(gi_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_z_xxxy_xzzzzz = cbuffer.data(gi_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyyy = cbuffer.data(gi_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyyz = cbuffer.data(gi_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyyzz = cbuffer.data(gi_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_z_xxxy_yyyzzz = cbuffer.data(gi_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_z_xxxy_yyzzzz = cbuffer.data(gi_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_z_xxxy_yzzzzz = cbuffer.data(gi_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_z_xxxy_zzzzzz = cbuffer.data(gi_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxx = cbuffer.data(gi_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxy = cbuffer.data(gi_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxxz = cbuffer.data(gi_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxyy = cbuffer.data(gi_geom_01_off + 899 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxyz = cbuffer.data(gi_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxxzz = cbuffer.data(gi_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyyy = cbuffer.data(gi_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyyz = cbuffer.data(gi_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxyzz = cbuffer.data(gi_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_z_xxxz_xxxzzz = cbuffer.data(gi_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyyy = cbuffer.data(gi_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyyz = cbuffer.data(gi_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyyzz = cbuffer.data(gi_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_z_xxxz_xxyzzz = cbuffer.data(gi_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_z_xxxz_xxzzzz = cbuffer.data(gi_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyyy = cbuffer.data(gi_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyyz = cbuffer.data(gi_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyyzz = cbuffer.data(gi_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_z_xxxz_xyyzzz = cbuffer.data(gi_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_z_xxxz_xyzzzz = cbuffer.data(gi_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_z_xxxz_xzzzzz = cbuffer.data(gi_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyyy = cbuffer.data(gi_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyyz = cbuffer.data(gi_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyyzz = cbuffer.data(gi_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_z_xxxz_yyyzzz = cbuffer.data(gi_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_z_xxxz_yyzzzz = cbuffer.data(gi_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_z_xxxz_yzzzzz = cbuffer.data(gi_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_z_xxxz_zzzzzz = cbuffer.data(gi_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxx = cbuffer.data(gi_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxy = cbuffer.data(gi_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxxz = cbuffer.data(gi_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxyy = cbuffer.data(gi_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxyz = cbuffer.data(gi_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxxzz = cbuffer.data(gi_geom_01_off + 929 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyyy = cbuffer.data(gi_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyyz = cbuffer.data(gi_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxyzz = cbuffer.data(gi_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_z_xxyy_xxxzzz = cbuffer.data(gi_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyyy = cbuffer.data(gi_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyyz = cbuffer.data(gi_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyyzz = cbuffer.data(gi_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_z_xxyy_xxyzzz = cbuffer.data(gi_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_z_xxyy_xxzzzz = cbuffer.data(gi_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyyy = cbuffer.data(gi_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyyz = cbuffer.data(gi_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyyzz = cbuffer.data(gi_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_z_xxyy_xyyzzz = cbuffer.data(gi_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_z_xxyy_xyzzzz = cbuffer.data(gi_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_z_xxyy_xzzzzz = cbuffer.data(gi_geom_01_off + 944 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyyy = cbuffer.data(gi_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyyz = cbuffer.data(gi_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyyzz = cbuffer.data(gi_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_z_xxyy_yyyzzz = cbuffer.data(gi_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_z_xxyy_yyzzzz = cbuffer.data(gi_geom_01_off + 949 * ccomps * dcomps);

            auto g_0_z_xxyy_yzzzzz = cbuffer.data(gi_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_z_xxyy_zzzzzz = cbuffer.data(gi_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxx = cbuffer.data(gi_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxy = cbuffer.data(gi_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxxz = cbuffer.data(gi_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxyy = cbuffer.data(gi_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxyz = cbuffer.data(gi_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxxzz = cbuffer.data(gi_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyyy = cbuffer.data(gi_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyyz = cbuffer.data(gi_geom_01_off + 959 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxyzz = cbuffer.data(gi_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_z_xxyz_xxxzzz = cbuffer.data(gi_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyyy = cbuffer.data(gi_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyyz = cbuffer.data(gi_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyyzz = cbuffer.data(gi_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_z_xxyz_xxyzzz = cbuffer.data(gi_geom_01_off + 965 * ccomps * dcomps);

            auto g_0_z_xxyz_xxzzzz = cbuffer.data(gi_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyyy = cbuffer.data(gi_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyyz = cbuffer.data(gi_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyyzz = cbuffer.data(gi_geom_01_off + 969 * ccomps * dcomps);

            auto g_0_z_xxyz_xyyzzz = cbuffer.data(gi_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_z_xxyz_xyzzzz = cbuffer.data(gi_geom_01_off + 971 * ccomps * dcomps);

            auto g_0_z_xxyz_xzzzzz = cbuffer.data(gi_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyyy = cbuffer.data(gi_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyyz = cbuffer.data(gi_geom_01_off + 974 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyyzz = cbuffer.data(gi_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_z_xxyz_yyyzzz = cbuffer.data(gi_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_z_xxyz_yyzzzz = cbuffer.data(gi_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_z_xxyz_yzzzzz = cbuffer.data(gi_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_z_xxyz_zzzzzz = cbuffer.data(gi_geom_01_off + 979 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxx = cbuffer.data(gi_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxy = cbuffer.data(gi_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxxz = cbuffer.data(gi_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxyy = cbuffer.data(gi_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxyz = cbuffer.data(gi_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxxzz = cbuffer.data(gi_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyyy = cbuffer.data(gi_geom_01_off + 986 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyyz = cbuffer.data(gi_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxyzz = cbuffer.data(gi_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_z_xxzz_xxxzzz = cbuffer.data(gi_geom_01_off + 989 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyyy = cbuffer.data(gi_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyyz = cbuffer.data(gi_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyyzz = cbuffer.data(gi_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_z_xxzz_xxyzzz = cbuffer.data(gi_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_z_xxzz_xxzzzz = cbuffer.data(gi_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyyy = cbuffer.data(gi_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyyz = cbuffer.data(gi_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyyzz = cbuffer.data(gi_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_z_xxzz_xyyzzz = cbuffer.data(gi_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_z_xxzz_xyzzzz = cbuffer.data(gi_geom_01_off + 999 * ccomps * dcomps);

            auto g_0_z_xxzz_xzzzzz = cbuffer.data(gi_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyyy = cbuffer.data(gi_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyyz = cbuffer.data(gi_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyyzz = cbuffer.data(gi_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_z_xxzz_yyyzzz = cbuffer.data(gi_geom_01_off + 1004 * ccomps * dcomps);

            auto g_0_z_xxzz_yyzzzz = cbuffer.data(gi_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_z_xxzz_yzzzzz = cbuffer.data(gi_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_z_xxzz_zzzzzz = cbuffer.data(gi_geom_01_off + 1007 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxx = cbuffer.data(gi_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxy = cbuffer.data(gi_geom_01_off + 1009 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxxz = cbuffer.data(gi_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxyy = cbuffer.data(gi_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxyz = cbuffer.data(gi_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxxzz = cbuffer.data(gi_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyyy = cbuffer.data(gi_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyyz = cbuffer.data(gi_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxyzz = cbuffer.data(gi_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_z_xyyy_xxxzzz = cbuffer.data(gi_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyyy = cbuffer.data(gi_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyyz = cbuffer.data(gi_geom_01_off + 1019 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyyzz = cbuffer.data(gi_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_z_xyyy_xxyzzz = cbuffer.data(gi_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_z_xyyy_xxzzzz = cbuffer.data(gi_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyyy = cbuffer.data(gi_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyyz = cbuffer.data(gi_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyyzz = cbuffer.data(gi_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_z_xyyy_xyyzzz = cbuffer.data(gi_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_z_xyyy_xyzzzz = cbuffer.data(gi_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_z_xyyy_xzzzzz = cbuffer.data(gi_geom_01_off + 1028 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyyy = cbuffer.data(gi_geom_01_off + 1029 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyyz = cbuffer.data(gi_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyyzz = cbuffer.data(gi_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_z_xyyy_yyyzzz = cbuffer.data(gi_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_z_xyyy_yyzzzz = cbuffer.data(gi_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_z_xyyy_yzzzzz = cbuffer.data(gi_geom_01_off + 1034 * ccomps * dcomps);

            auto g_0_z_xyyy_zzzzzz = cbuffer.data(gi_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxx = cbuffer.data(gi_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxy = cbuffer.data(gi_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxxz = cbuffer.data(gi_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxyy = cbuffer.data(gi_geom_01_off + 1039 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxyz = cbuffer.data(gi_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxxzz = cbuffer.data(gi_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyyy = cbuffer.data(gi_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyyz = cbuffer.data(gi_geom_01_off + 1043 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxyzz = cbuffer.data(gi_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_z_xyyz_xxxzzz = cbuffer.data(gi_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyyy = cbuffer.data(gi_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyyz = cbuffer.data(gi_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyyzz = cbuffer.data(gi_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_z_xyyz_xxyzzz = cbuffer.data(gi_geom_01_off + 1049 * ccomps * dcomps);

            auto g_0_z_xyyz_xxzzzz = cbuffer.data(gi_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyyy = cbuffer.data(gi_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyyz = cbuffer.data(gi_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyyzz = cbuffer.data(gi_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_z_xyyz_xyyzzz = cbuffer.data(gi_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_z_xyyz_xyzzzz = cbuffer.data(gi_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_z_xyyz_xzzzzz = cbuffer.data(gi_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyyy = cbuffer.data(gi_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyyz = cbuffer.data(gi_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyyzz = cbuffer.data(gi_geom_01_off + 1059 * ccomps * dcomps);

            auto g_0_z_xyyz_yyyzzz = cbuffer.data(gi_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_z_xyyz_yyzzzz = cbuffer.data(gi_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_z_xyyz_yzzzzz = cbuffer.data(gi_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_z_xyyz_zzzzzz = cbuffer.data(gi_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxx = cbuffer.data(gi_geom_01_off + 1064 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxy = cbuffer.data(gi_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxxz = cbuffer.data(gi_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxyy = cbuffer.data(gi_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxyz = cbuffer.data(gi_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxxzz = cbuffer.data(gi_geom_01_off + 1069 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyyy = cbuffer.data(gi_geom_01_off + 1070 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyyz = cbuffer.data(gi_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxyzz = cbuffer.data(gi_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_z_xyzz_xxxzzz = cbuffer.data(gi_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyyy = cbuffer.data(gi_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyyz = cbuffer.data(gi_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyyzz = cbuffer.data(gi_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_z_xyzz_xxyzzz = cbuffer.data(gi_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_z_xyzz_xxzzzz = cbuffer.data(gi_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyyy = cbuffer.data(gi_geom_01_off + 1079 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyyz = cbuffer.data(gi_geom_01_off + 1080 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyyzz = cbuffer.data(gi_geom_01_off + 1081 * ccomps * dcomps);

            auto g_0_z_xyzz_xyyzzz = cbuffer.data(gi_geom_01_off + 1082 * ccomps * dcomps);

            auto g_0_z_xyzz_xyzzzz = cbuffer.data(gi_geom_01_off + 1083 * ccomps * dcomps);

            auto g_0_z_xyzz_xzzzzz = cbuffer.data(gi_geom_01_off + 1084 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyyy = cbuffer.data(gi_geom_01_off + 1085 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyyz = cbuffer.data(gi_geom_01_off + 1086 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyyzz = cbuffer.data(gi_geom_01_off + 1087 * ccomps * dcomps);

            auto g_0_z_xyzz_yyyzzz = cbuffer.data(gi_geom_01_off + 1088 * ccomps * dcomps);

            auto g_0_z_xyzz_yyzzzz = cbuffer.data(gi_geom_01_off + 1089 * ccomps * dcomps);

            auto g_0_z_xyzz_yzzzzz = cbuffer.data(gi_geom_01_off + 1090 * ccomps * dcomps);

            auto g_0_z_xyzz_zzzzzz = cbuffer.data(gi_geom_01_off + 1091 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxx = cbuffer.data(gi_geom_01_off + 1092 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxy = cbuffer.data(gi_geom_01_off + 1093 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxxz = cbuffer.data(gi_geom_01_off + 1094 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxyy = cbuffer.data(gi_geom_01_off + 1095 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxyz = cbuffer.data(gi_geom_01_off + 1096 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxxzz = cbuffer.data(gi_geom_01_off + 1097 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyyy = cbuffer.data(gi_geom_01_off + 1098 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyyz = cbuffer.data(gi_geom_01_off + 1099 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxyzz = cbuffer.data(gi_geom_01_off + 1100 * ccomps * dcomps);

            auto g_0_z_xzzz_xxxzzz = cbuffer.data(gi_geom_01_off + 1101 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyyy = cbuffer.data(gi_geom_01_off + 1102 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyyz = cbuffer.data(gi_geom_01_off + 1103 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyyzz = cbuffer.data(gi_geom_01_off + 1104 * ccomps * dcomps);

            auto g_0_z_xzzz_xxyzzz = cbuffer.data(gi_geom_01_off + 1105 * ccomps * dcomps);

            auto g_0_z_xzzz_xxzzzz = cbuffer.data(gi_geom_01_off + 1106 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyyy = cbuffer.data(gi_geom_01_off + 1107 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyyz = cbuffer.data(gi_geom_01_off + 1108 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyyzz = cbuffer.data(gi_geom_01_off + 1109 * ccomps * dcomps);

            auto g_0_z_xzzz_xyyzzz = cbuffer.data(gi_geom_01_off + 1110 * ccomps * dcomps);

            auto g_0_z_xzzz_xyzzzz = cbuffer.data(gi_geom_01_off + 1111 * ccomps * dcomps);

            auto g_0_z_xzzz_xzzzzz = cbuffer.data(gi_geom_01_off + 1112 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyyy = cbuffer.data(gi_geom_01_off + 1113 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyyz = cbuffer.data(gi_geom_01_off + 1114 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyyzz = cbuffer.data(gi_geom_01_off + 1115 * ccomps * dcomps);

            auto g_0_z_xzzz_yyyzzz = cbuffer.data(gi_geom_01_off + 1116 * ccomps * dcomps);

            auto g_0_z_xzzz_yyzzzz = cbuffer.data(gi_geom_01_off + 1117 * ccomps * dcomps);

            auto g_0_z_xzzz_yzzzzz = cbuffer.data(gi_geom_01_off + 1118 * ccomps * dcomps);

            auto g_0_z_xzzz_zzzzzz = cbuffer.data(gi_geom_01_off + 1119 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxx = cbuffer.data(gi_geom_01_off + 1120 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxy = cbuffer.data(gi_geom_01_off + 1121 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxxz = cbuffer.data(gi_geom_01_off + 1122 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxyy = cbuffer.data(gi_geom_01_off + 1123 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxyz = cbuffer.data(gi_geom_01_off + 1124 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxxzz = cbuffer.data(gi_geom_01_off + 1125 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyyy = cbuffer.data(gi_geom_01_off + 1126 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyyz = cbuffer.data(gi_geom_01_off + 1127 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxyzz = cbuffer.data(gi_geom_01_off + 1128 * ccomps * dcomps);

            auto g_0_z_yyyy_xxxzzz = cbuffer.data(gi_geom_01_off + 1129 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyyy = cbuffer.data(gi_geom_01_off + 1130 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyyz = cbuffer.data(gi_geom_01_off + 1131 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyyzz = cbuffer.data(gi_geom_01_off + 1132 * ccomps * dcomps);

            auto g_0_z_yyyy_xxyzzz = cbuffer.data(gi_geom_01_off + 1133 * ccomps * dcomps);

            auto g_0_z_yyyy_xxzzzz = cbuffer.data(gi_geom_01_off + 1134 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyyy = cbuffer.data(gi_geom_01_off + 1135 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyyz = cbuffer.data(gi_geom_01_off + 1136 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyyzz = cbuffer.data(gi_geom_01_off + 1137 * ccomps * dcomps);

            auto g_0_z_yyyy_xyyzzz = cbuffer.data(gi_geom_01_off + 1138 * ccomps * dcomps);

            auto g_0_z_yyyy_xyzzzz = cbuffer.data(gi_geom_01_off + 1139 * ccomps * dcomps);

            auto g_0_z_yyyy_xzzzzz = cbuffer.data(gi_geom_01_off + 1140 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyyy = cbuffer.data(gi_geom_01_off + 1141 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyyz = cbuffer.data(gi_geom_01_off + 1142 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyyzz = cbuffer.data(gi_geom_01_off + 1143 * ccomps * dcomps);

            auto g_0_z_yyyy_yyyzzz = cbuffer.data(gi_geom_01_off + 1144 * ccomps * dcomps);

            auto g_0_z_yyyy_yyzzzz = cbuffer.data(gi_geom_01_off + 1145 * ccomps * dcomps);

            auto g_0_z_yyyy_yzzzzz = cbuffer.data(gi_geom_01_off + 1146 * ccomps * dcomps);

            auto g_0_z_yyyy_zzzzzz = cbuffer.data(gi_geom_01_off + 1147 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxx = cbuffer.data(gi_geom_01_off + 1148 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxy = cbuffer.data(gi_geom_01_off + 1149 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxxz = cbuffer.data(gi_geom_01_off + 1150 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxyy = cbuffer.data(gi_geom_01_off + 1151 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxyz = cbuffer.data(gi_geom_01_off + 1152 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxxzz = cbuffer.data(gi_geom_01_off + 1153 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyyy = cbuffer.data(gi_geom_01_off + 1154 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyyz = cbuffer.data(gi_geom_01_off + 1155 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxyzz = cbuffer.data(gi_geom_01_off + 1156 * ccomps * dcomps);

            auto g_0_z_yyyz_xxxzzz = cbuffer.data(gi_geom_01_off + 1157 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyyy = cbuffer.data(gi_geom_01_off + 1158 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyyz = cbuffer.data(gi_geom_01_off + 1159 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyyzz = cbuffer.data(gi_geom_01_off + 1160 * ccomps * dcomps);

            auto g_0_z_yyyz_xxyzzz = cbuffer.data(gi_geom_01_off + 1161 * ccomps * dcomps);

            auto g_0_z_yyyz_xxzzzz = cbuffer.data(gi_geom_01_off + 1162 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyyy = cbuffer.data(gi_geom_01_off + 1163 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyyz = cbuffer.data(gi_geom_01_off + 1164 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyyzz = cbuffer.data(gi_geom_01_off + 1165 * ccomps * dcomps);

            auto g_0_z_yyyz_xyyzzz = cbuffer.data(gi_geom_01_off + 1166 * ccomps * dcomps);

            auto g_0_z_yyyz_xyzzzz = cbuffer.data(gi_geom_01_off + 1167 * ccomps * dcomps);

            auto g_0_z_yyyz_xzzzzz = cbuffer.data(gi_geom_01_off + 1168 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyyy = cbuffer.data(gi_geom_01_off + 1169 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyyz = cbuffer.data(gi_geom_01_off + 1170 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyyzz = cbuffer.data(gi_geom_01_off + 1171 * ccomps * dcomps);

            auto g_0_z_yyyz_yyyzzz = cbuffer.data(gi_geom_01_off + 1172 * ccomps * dcomps);

            auto g_0_z_yyyz_yyzzzz = cbuffer.data(gi_geom_01_off + 1173 * ccomps * dcomps);

            auto g_0_z_yyyz_yzzzzz = cbuffer.data(gi_geom_01_off + 1174 * ccomps * dcomps);

            auto g_0_z_yyyz_zzzzzz = cbuffer.data(gi_geom_01_off + 1175 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxx = cbuffer.data(gi_geom_01_off + 1176 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxy = cbuffer.data(gi_geom_01_off + 1177 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxxz = cbuffer.data(gi_geom_01_off + 1178 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxyy = cbuffer.data(gi_geom_01_off + 1179 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxyz = cbuffer.data(gi_geom_01_off + 1180 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxxzz = cbuffer.data(gi_geom_01_off + 1181 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyyy = cbuffer.data(gi_geom_01_off + 1182 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyyz = cbuffer.data(gi_geom_01_off + 1183 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxyzz = cbuffer.data(gi_geom_01_off + 1184 * ccomps * dcomps);

            auto g_0_z_yyzz_xxxzzz = cbuffer.data(gi_geom_01_off + 1185 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyyy = cbuffer.data(gi_geom_01_off + 1186 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyyz = cbuffer.data(gi_geom_01_off + 1187 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyyzz = cbuffer.data(gi_geom_01_off + 1188 * ccomps * dcomps);

            auto g_0_z_yyzz_xxyzzz = cbuffer.data(gi_geom_01_off + 1189 * ccomps * dcomps);

            auto g_0_z_yyzz_xxzzzz = cbuffer.data(gi_geom_01_off + 1190 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyyy = cbuffer.data(gi_geom_01_off + 1191 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyyz = cbuffer.data(gi_geom_01_off + 1192 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyyzz = cbuffer.data(gi_geom_01_off + 1193 * ccomps * dcomps);

            auto g_0_z_yyzz_xyyzzz = cbuffer.data(gi_geom_01_off + 1194 * ccomps * dcomps);

            auto g_0_z_yyzz_xyzzzz = cbuffer.data(gi_geom_01_off + 1195 * ccomps * dcomps);

            auto g_0_z_yyzz_xzzzzz = cbuffer.data(gi_geom_01_off + 1196 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyyy = cbuffer.data(gi_geom_01_off + 1197 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyyz = cbuffer.data(gi_geom_01_off + 1198 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyyzz = cbuffer.data(gi_geom_01_off + 1199 * ccomps * dcomps);

            auto g_0_z_yyzz_yyyzzz = cbuffer.data(gi_geom_01_off + 1200 * ccomps * dcomps);

            auto g_0_z_yyzz_yyzzzz = cbuffer.data(gi_geom_01_off + 1201 * ccomps * dcomps);

            auto g_0_z_yyzz_yzzzzz = cbuffer.data(gi_geom_01_off + 1202 * ccomps * dcomps);

            auto g_0_z_yyzz_zzzzzz = cbuffer.data(gi_geom_01_off + 1203 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxx = cbuffer.data(gi_geom_01_off + 1204 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxy = cbuffer.data(gi_geom_01_off + 1205 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxxz = cbuffer.data(gi_geom_01_off + 1206 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxyy = cbuffer.data(gi_geom_01_off + 1207 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxyz = cbuffer.data(gi_geom_01_off + 1208 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxxzz = cbuffer.data(gi_geom_01_off + 1209 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyyy = cbuffer.data(gi_geom_01_off + 1210 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyyz = cbuffer.data(gi_geom_01_off + 1211 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxyzz = cbuffer.data(gi_geom_01_off + 1212 * ccomps * dcomps);

            auto g_0_z_yzzz_xxxzzz = cbuffer.data(gi_geom_01_off + 1213 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyyy = cbuffer.data(gi_geom_01_off + 1214 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyyz = cbuffer.data(gi_geom_01_off + 1215 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyyzz = cbuffer.data(gi_geom_01_off + 1216 * ccomps * dcomps);

            auto g_0_z_yzzz_xxyzzz = cbuffer.data(gi_geom_01_off + 1217 * ccomps * dcomps);

            auto g_0_z_yzzz_xxzzzz = cbuffer.data(gi_geom_01_off + 1218 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyyy = cbuffer.data(gi_geom_01_off + 1219 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyyz = cbuffer.data(gi_geom_01_off + 1220 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyyzz = cbuffer.data(gi_geom_01_off + 1221 * ccomps * dcomps);

            auto g_0_z_yzzz_xyyzzz = cbuffer.data(gi_geom_01_off + 1222 * ccomps * dcomps);

            auto g_0_z_yzzz_xyzzzz = cbuffer.data(gi_geom_01_off + 1223 * ccomps * dcomps);

            auto g_0_z_yzzz_xzzzzz = cbuffer.data(gi_geom_01_off + 1224 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyyy = cbuffer.data(gi_geom_01_off + 1225 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyyz = cbuffer.data(gi_geom_01_off + 1226 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyyzz = cbuffer.data(gi_geom_01_off + 1227 * ccomps * dcomps);

            auto g_0_z_yzzz_yyyzzz = cbuffer.data(gi_geom_01_off + 1228 * ccomps * dcomps);

            auto g_0_z_yzzz_yyzzzz = cbuffer.data(gi_geom_01_off + 1229 * ccomps * dcomps);

            auto g_0_z_yzzz_yzzzzz = cbuffer.data(gi_geom_01_off + 1230 * ccomps * dcomps);

            auto g_0_z_yzzz_zzzzzz = cbuffer.data(gi_geom_01_off + 1231 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxx = cbuffer.data(gi_geom_01_off + 1232 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxy = cbuffer.data(gi_geom_01_off + 1233 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxxz = cbuffer.data(gi_geom_01_off + 1234 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxyy = cbuffer.data(gi_geom_01_off + 1235 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxyz = cbuffer.data(gi_geom_01_off + 1236 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxxzz = cbuffer.data(gi_geom_01_off + 1237 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyyy = cbuffer.data(gi_geom_01_off + 1238 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyyz = cbuffer.data(gi_geom_01_off + 1239 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxyzz = cbuffer.data(gi_geom_01_off + 1240 * ccomps * dcomps);

            auto g_0_z_zzzz_xxxzzz = cbuffer.data(gi_geom_01_off + 1241 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyyy = cbuffer.data(gi_geom_01_off + 1242 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyyz = cbuffer.data(gi_geom_01_off + 1243 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyyzz = cbuffer.data(gi_geom_01_off + 1244 * ccomps * dcomps);

            auto g_0_z_zzzz_xxyzzz = cbuffer.data(gi_geom_01_off + 1245 * ccomps * dcomps);

            auto g_0_z_zzzz_xxzzzz = cbuffer.data(gi_geom_01_off + 1246 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyyy = cbuffer.data(gi_geom_01_off + 1247 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyyz = cbuffer.data(gi_geom_01_off + 1248 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyyzz = cbuffer.data(gi_geom_01_off + 1249 * ccomps * dcomps);

            auto g_0_z_zzzz_xyyzzz = cbuffer.data(gi_geom_01_off + 1250 * ccomps * dcomps);

            auto g_0_z_zzzz_xyzzzz = cbuffer.data(gi_geom_01_off + 1251 * ccomps * dcomps);

            auto g_0_z_zzzz_xzzzzz = cbuffer.data(gi_geom_01_off + 1252 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyyy = cbuffer.data(gi_geom_01_off + 1253 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyyz = cbuffer.data(gi_geom_01_off + 1254 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyyzz = cbuffer.data(gi_geom_01_off + 1255 * ccomps * dcomps);

            auto g_0_z_zzzz_yyyzzz = cbuffer.data(gi_geom_01_off + 1256 * ccomps * dcomps);

            auto g_0_z_zzzz_yyzzzz = cbuffer.data(gi_geom_01_off + 1257 * ccomps * dcomps);

            auto g_0_z_zzzz_yzzzzz = cbuffer.data(gi_geom_01_off + 1258 * ccomps * dcomps);

            auto g_0_z_zzzz_zzzzzz = cbuffer.data(gi_geom_01_off + 1259 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_hhxx

            const auto hh_geom_01_off = idx_geom_01_hhxx + i * dcomps + j;

            /// Set up 0-21 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxx_xxxxx = cbuffer.data(hh_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxy = cbuffer.data(hh_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxxz = cbuffer.data(hh_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxyy = cbuffer.data(hh_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxyz = cbuffer.data(hh_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxzz = cbuffer.data(hh_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxyyy = cbuffer.data(hh_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxyyz = cbuffer.data(hh_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxyzz = cbuffer.data(hh_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxzzz = cbuffer.data(hh_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyyyy = cbuffer.data(hh_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyyyz = cbuffer.data(hh_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyyzz = cbuffer.data(hh_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyzzz = cbuffer.data(hh_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxx_xzzzz = cbuffer.data(hh_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyyyy = cbuffer.data(hh_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyyyz = cbuffer.data(hh_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyyzz = cbuffer.data(hh_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyzzz = cbuffer.data(hh_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxx_yzzzz = cbuffer.data(hh_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxx_zzzzz = cbuffer.data(hh_geom_01_off + 20 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_xxxxx, g_0_x_xxxx_xxxxxx, g_0_x_xxxx_xxxxxy, g_0_x_xxxx_xxxxxz, g_0_x_xxxx_xxxxy, g_0_x_xxxx_xxxxyy, g_0_x_xxxx_xxxxyz, g_0_x_xxxx_xxxxz, g_0_x_xxxx_xxxxzz, g_0_x_xxxx_xxxyy, g_0_x_xxxx_xxxyyy, g_0_x_xxxx_xxxyyz, g_0_x_xxxx_xxxyz, g_0_x_xxxx_xxxyzz, g_0_x_xxxx_xxxzz, g_0_x_xxxx_xxxzzz, g_0_x_xxxx_xxyyy, g_0_x_xxxx_xxyyyy, g_0_x_xxxx_xxyyyz, g_0_x_xxxx_xxyyz, g_0_x_xxxx_xxyyzz, g_0_x_xxxx_xxyzz, g_0_x_xxxx_xxyzzz, g_0_x_xxxx_xxzzz, g_0_x_xxxx_xxzzzz, g_0_x_xxxx_xyyyy, g_0_x_xxxx_xyyyyy, g_0_x_xxxx_xyyyyz, g_0_x_xxxx_xyyyz, g_0_x_xxxx_xyyyzz, g_0_x_xxxx_xyyzz, g_0_x_xxxx_xyyzzz, g_0_x_xxxx_xyzzz, g_0_x_xxxx_xyzzzz, g_0_x_xxxx_xzzzz, g_0_x_xxxx_xzzzzz, g_0_x_xxxx_yyyyy, g_0_x_xxxx_yyyyz, g_0_x_xxxx_yyyzz, g_0_x_xxxx_yyzzz, g_0_x_xxxx_yzzzz, g_0_x_xxxx_zzzzz, g_0_x_xxxxx_xxxxx, g_0_x_xxxxx_xxxxy, g_0_x_xxxxx_xxxxz, g_0_x_xxxxx_xxxyy, g_0_x_xxxxx_xxxyz, g_0_x_xxxxx_xxxzz, g_0_x_xxxxx_xxyyy, g_0_x_xxxxx_xxyyz, g_0_x_xxxxx_xxyzz, g_0_x_xxxxx_xxzzz, g_0_x_xxxxx_xyyyy, g_0_x_xxxxx_xyyyz, g_0_x_xxxxx_xyyzz, g_0_x_xxxxx_xyzzz, g_0_x_xxxxx_xzzzz, g_0_x_xxxxx_yyyyy, g_0_x_xxxxx_yyyyz, g_0_x_xxxxx_yyyzz, g_0_x_xxxxx_yyzzz, g_0_x_xxxxx_yzzzz, g_0_x_xxxxx_zzzzz, g_xxxx_xxxxx, g_xxxx_xxxxy, g_xxxx_xxxxz, g_xxxx_xxxyy, g_xxxx_xxxyz, g_xxxx_xxxzz, g_xxxx_xxyyy, g_xxxx_xxyyz, g_xxxx_xxyzz, g_xxxx_xxzzz, g_xxxx_xyyyy, g_xxxx_xyyyz, g_xxxx_xyyzz, g_xxxx_xyzzz, g_xxxx_xzzzz, g_xxxx_yyyyy, g_xxxx_yyyyz, g_xxxx_yyyzz, g_xxxx_yyzzz, g_xxxx_yzzzz, g_xxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxx_xxxxx[k] = g_xxxx_xxxxx[k] - g_0_x_xxxx_xxxxx[k] * ab_x + g_0_x_xxxx_xxxxxx[k];

                g_0_x_xxxxx_xxxxy[k] = g_xxxx_xxxxy[k] - g_0_x_xxxx_xxxxy[k] * ab_x + g_0_x_xxxx_xxxxxy[k];

                g_0_x_xxxxx_xxxxz[k] = g_xxxx_xxxxz[k] - g_0_x_xxxx_xxxxz[k] * ab_x + g_0_x_xxxx_xxxxxz[k];

                g_0_x_xxxxx_xxxyy[k] = g_xxxx_xxxyy[k] - g_0_x_xxxx_xxxyy[k] * ab_x + g_0_x_xxxx_xxxxyy[k];

                g_0_x_xxxxx_xxxyz[k] = g_xxxx_xxxyz[k] - g_0_x_xxxx_xxxyz[k] * ab_x + g_0_x_xxxx_xxxxyz[k];

                g_0_x_xxxxx_xxxzz[k] = g_xxxx_xxxzz[k] - g_0_x_xxxx_xxxzz[k] * ab_x + g_0_x_xxxx_xxxxzz[k];

                g_0_x_xxxxx_xxyyy[k] = g_xxxx_xxyyy[k] - g_0_x_xxxx_xxyyy[k] * ab_x + g_0_x_xxxx_xxxyyy[k];

                g_0_x_xxxxx_xxyyz[k] = g_xxxx_xxyyz[k] - g_0_x_xxxx_xxyyz[k] * ab_x + g_0_x_xxxx_xxxyyz[k];

                g_0_x_xxxxx_xxyzz[k] = g_xxxx_xxyzz[k] - g_0_x_xxxx_xxyzz[k] * ab_x + g_0_x_xxxx_xxxyzz[k];

                g_0_x_xxxxx_xxzzz[k] = g_xxxx_xxzzz[k] - g_0_x_xxxx_xxzzz[k] * ab_x + g_0_x_xxxx_xxxzzz[k];

                g_0_x_xxxxx_xyyyy[k] = g_xxxx_xyyyy[k] - g_0_x_xxxx_xyyyy[k] * ab_x + g_0_x_xxxx_xxyyyy[k];

                g_0_x_xxxxx_xyyyz[k] = g_xxxx_xyyyz[k] - g_0_x_xxxx_xyyyz[k] * ab_x + g_0_x_xxxx_xxyyyz[k];

                g_0_x_xxxxx_xyyzz[k] = g_xxxx_xyyzz[k] - g_0_x_xxxx_xyyzz[k] * ab_x + g_0_x_xxxx_xxyyzz[k];

                g_0_x_xxxxx_xyzzz[k] = g_xxxx_xyzzz[k] - g_0_x_xxxx_xyzzz[k] * ab_x + g_0_x_xxxx_xxyzzz[k];

                g_0_x_xxxxx_xzzzz[k] = g_xxxx_xzzzz[k] - g_0_x_xxxx_xzzzz[k] * ab_x + g_0_x_xxxx_xxzzzz[k];

                g_0_x_xxxxx_yyyyy[k] = g_xxxx_yyyyy[k] - g_0_x_xxxx_yyyyy[k] * ab_x + g_0_x_xxxx_xyyyyy[k];

                g_0_x_xxxxx_yyyyz[k] = g_xxxx_yyyyz[k] - g_0_x_xxxx_yyyyz[k] * ab_x + g_0_x_xxxx_xyyyyz[k];

                g_0_x_xxxxx_yyyzz[k] = g_xxxx_yyyzz[k] - g_0_x_xxxx_yyyzz[k] * ab_x + g_0_x_xxxx_xyyyzz[k];

                g_0_x_xxxxx_yyzzz[k] = g_xxxx_yyzzz[k] - g_0_x_xxxx_yyzzz[k] * ab_x + g_0_x_xxxx_xyyzzz[k];

                g_0_x_xxxxx_yzzzz[k] = g_xxxx_yzzzz[k] - g_0_x_xxxx_yzzzz[k] * ab_x + g_0_x_xxxx_xyzzzz[k];

                g_0_x_xxxxx_zzzzz[k] = g_xxxx_zzzzz[k] - g_0_x_xxxx_zzzzz[k] * ab_x + g_0_x_xxxx_xzzzzz[k];
            }

            /// Set up 21-42 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxy_xxxxx = cbuffer.data(hh_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxy = cbuffer.data(hh_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxxz = cbuffer.data(hh_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxyy = cbuffer.data(hh_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxyz = cbuffer.data(hh_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxzz = cbuffer.data(hh_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxyyy = cbuffer.data(hh_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxyyz = cbuffer.data(hh_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxyzz = cbuffer.data(hh_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxzzz = cbuffer.data(hh_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyyyy = cbuffer.data(hh_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyyyz = cbuffer.data(hh_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyyzz = cbuffer.data(hh_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyzzz = cbuffer.data(hh_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxy_xzzzz = cbuffer.data(hh_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyyyy = cbuffer.data(hh_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyyyz = cbuffer.data(hh_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyyzz = cbuffer.data(hh_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyzzz = cbuffer.data(hh_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxxy_yzzzz = cbuffer.data(hh_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxxy_zzzzz = cbuffer.data(hh_geom_01_off + 41 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_xxxxx, g_0_x_xxxx_xxxxxy, g_0_x_xxxx_xxxxy, g_0_x_xxxx_xxxxyy, g_0_x_xxxx_xxxxyz, g_0_x_xxxx_xxxxz, g_0_x_xxxx_xxxyy, g_0_x_xxxx_xxxyyy, g_0_x_xxxx_xxxyyz, g_0_x_xxxx_xxxyz, g_0_x_xxxx_xxxyzz, g_0_x_xxxx_xxxzz, g_0_x_xxxx_xxyyy, g_0_x_xxxx_xxyyyy, g_0_x_xxxx_xxyyyz, g_0_x_xxxx_xxyyz, g_0_x_xxxx_xxyyzz, g_0_x_xxxx_xxyzz, g_0_x_xxxx_xxyzzz, g_0_x_xxxx_xxzzz, g_0_x_xxxx_xyyyy, g_0_x_xxxx_xyyyyy, g_0_x_xxxx_xyyyyz, g_0_x_xxxx_xyyyz, g_0_x_xxxx_xyyyzz, g_0_x_xxxx_xyyzz, g_0_x_xxxx_xyyzzz, g_0_x_xxxx_xyzzz, g_0_x_xxxx_xyzzzz, g_0_x_xxxx_xzzzz, g_0_x_xxxx_yyyyy, g_0_x_xxxx_yyyyyy, g_0_x_xxxx_yyyyyz, g_0_x_xxxx_yyyyz, g_0_x_xxxx_yyyyzz, g_0_x_xxxx_yyyzz, g_0_x_xxxx_yyyzzz, g_0_x_xxxx_yyzzz, g_0_x_xxxx_yyzzzz, g_0_x_xxxx_yzzzz, g_0_x_xxxx_yzzzzz, g_0_x_xxxx_zzzzz, g_0_x_xxxxy_xxxxx, g_0_x_xxxxy_xxxxy, g_0_x_xxxxy_xxxxz, g_0_x_xxxxy_xxxyy, g_0_x_xxxxy_xxxyz, g_0_x_xxxxy_xxxzz, g_0_x_xxxxy_xxyyy, g_0_x_xxxxy_xxyyz, g_0_x_xxxxy_xxyzz, g_0_x_xxxxy_xxzzz, g_0_x_xxxxy_xyyyy, g_0_x_xxxxy_xyyyz, g_0_x_xxxxy_xyyzz, g_0_x_xxxxy_xyzzz, g_0_x_xxxxy_xzzzz, g_0_x_xxxxy_yyyyy, g_0_x_xxxxy_yyyyz, g_0_x_xxxxy_yyyzz, g_0_x_xxxxy_yyzzz, g_0_x_xxxxy_yzzzz, g_0_x_xxxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxy_xxxxx[k] = -g_0_x_xxxx_xxxxx[k] * ab_y + g_0_x_xxxx_xxxxxy[k];

                g_0_x_xxxxy_xxxxy[k] = -g_0_x_xxxx_xxxxy[k] * ab_y + g_0_x_xxxx_xxxxyy[k];

                g_0_x_xxxxy_xxxxz[k] = -g_0_x_xxxx_xxxxz[k] * ab_y + g_0_x_xxxx_xxxxyz[k];

                g_0_x_xxxxy_xxxyy[k] = -g_0_x_xxxx_xxxyy[k] * ab_y + g_0_x_xxxx_xxxyyy[k];

                g_0_x_xxxxy_xxxyz[k] = -g_0_x_xxxx_xxxyz[k] * ab_y + g_0_x_xxxx_xxxyyz[k];

                g_0_x_xxxxy_xxxzz[k] = -g_0_x_xxxx_xxxzz[k] * ab_y + g_0_x_xxxx_xxxyzz[k];

                g_0_x_xxxxy_xxyyy[k] = -g_0_x_xxxx_xxyyy[k] * ab_y + g_0_x_xxxx_xxyyyy[k];

                g_0_x_xxxxy_xxyyz[k] = -g_0_x_xxxx_xxyyz[k] * ab_y + g_0_x_xxxx_xxyyyz[k];

                g_0_x_xxxxy_xxyzz[k] = -g_0_x_xxxx_xxyzz[k] * ab_y + g_0_x_xxxx_xxyyzz[k];

                g_0_x_xxxxy_xxzzz[k] = -g_0_x_xxxx_xxzzz[k] * ab_y + g_0_x_xxxx_xxyzzz[k];

                g_0_x_xxxxy_xyyyy[k] = -g_0_x_xxxx_xyyyy[k] * ab_y + g_0_x_xxxx_xyyyyy[k];

                g_0_x_xxxxy_xyyyz[k] = -g_0_x_xxxx_xyyyz[k] * ab_y + g_0_x_xxxx_xyyyyz[k];

                g_0_x_xxxxy_xyyzz[k] = -g_0_x_xxxx_xyyzz[k] * ab_y + g_0_x_xxxx_xyyyzz[k];

                g_0_x_xxxxy_xyzzz[k] = -g_0_x_xxxx_xyzzz[k] * ab_y + g_0_x_xxxx_xyyzzz[k];

                g_0_x_xxxxy_xzzzz[k] = -g_0_x_xxxx_xzzzz[k] * ab_y + g_0_x_xxxx_xyzzzz[k];

                g_0_x_xxxxy_yyyyy[k] = -g_0_x_xxxx_yyyyy[k] * ab_y + g_0_x_xxxx_yyyyyy[k];

                g_0_x_xxxxy_yyyyz[k] = -g_0_x_xxxx_yyyyz[k] * ab_y + g_0_x_xxxx_yyyyyz[k];

                g_0_x_xxxxy_yyyzz[k] = -g_0_x_xxxx_yyyzz[k] * ab_y + g_0_x_xxxx_yyyyzz[k];

                g_0_x_xxxxy_yyzzz[k] = -g_0_x_xxxx_yyzzz[k] * ab_y + g_0_x_xxxx_yyyzzz[k];

                g_0_x_xxxxy_yzzzz[k] = -g_0_x_xxxx_yzzzz[k] * ab_y + g_0_x_xxxx_yyzzzz[k];

                g_0_x_xxxxy_zzzzz[k] = -g_0_x_xxxx_zzzzz[k] * ab_y + g_0_x_xxxx_yzzzzz[k];
            }

            /// Set up 42-63 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxz_xxxxx = cbuffer.data(hh_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxy = cbuffer.data(hh_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxxz = cbuffer.data(hh_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxyy = cbuffer.data(hh_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxyz = cbuffer.data(hh_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxzz = cbuffer.data(hh_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxyyy = cbuffer.data(hh_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxyyz = cbuffer.data(hh_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxyzz = cbuffer.data(hh_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxzzz = cbuffer.data(hh_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyyyy = cbuffer.data(hh_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyyyz = cbuffer.data(hh_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyyzz = cbuffer.data(hh_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyzzz = cbuffer.data(hh_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxxz_xzzzz = cbuffer.data(hh_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyyyy = cbuffer.data(hh_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyyyz = cbuffer.data(hh_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyyzz = cbuffer.data(hh_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyzzz = cbuffer.data(hh_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxxz_yzzzz = cbuffer.data(hh_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxxz_zzzzz = cbuffer.data(hh_geom_01_off + 62 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxx_xxxxx, g_0_x_xxxx_xxxxxz, g_0_x_xxxx_xxxxy, g_0_x_xxxx_xxxxyz, g_0_x_xxxx_xxxxz, g_0_x_xxxx_xxxxzz, g_0_x_xxxx_xxxyy, g_0_x_xxxx_xxxyyz, g_0_x_xxxx_xxxyz, g_0_x_xxxx_xxxyzz, g_0_x_xxxx_xxxzz, g_0_x_xxxx_xxxzzz, g_0_x_xxxx_xxyyy, g_0_x_xxxx_xxyyyz, g_0_x_xxxx_xxyyz, g_0_x_xxxx_xxyyzz, g_0_x_xxxx_xxyzz, g_0_x_xxxx_xxyzzz, g_0_x_xxxx_xxzzz, g_0_x_xxxx_xxzzzz, g_0_x_xxxx_xyyyy, g_0_x_xxxx_xyyyyz, g_0_x_xxxx_xyyyz, g_0_x_xxxx_xyyyzz, g_0_x_xxxx_xyyzz, g_0_x_xxxx_xyyzzz, g_0_x_xxxx_xyzzz, g_0_x_xxxx_xyzzzz, g_0_x_xxxx_xzzzz, g_0_x_xxxx_xzzzzz, g_0_x_xxxx_yyyyy, g_0_x_xxxx_yyyyyz, g_0_x_xxxx_yyyyz, g_0_x_xxxx_yyyyzz, g_0_x_xxxx_yyyzz, g_0_x_xxxx_yyyzzz, g_0_x_xxxx_yyzzz, g_0_x_xxxx_yyzzzz, g_0_x_xxxx_yzzzz, g_0_x_xxxx_yzzzzz, g_0_x_xxxx_zzzzz, g_0_x_xxxx_zzzzzz, g_0_x_xxxxz_xxxxx, g_0_x_xxxxz_xxxxy, g_0_x_xxxxz_xxxxz, g_0_x_xxxxz_xxxyy, g_0_x_xxxxz_xxxyz, g_0_x_xxxxz_xxxzz, g_0_x_xxxxz_xxyyy, g_0_x_xxxxz_xxyyz, g_0_x_xxxxz_xxyzz, g_0_x_xxxxz_xxzzz, g_0_x_xxxxz_xyyyy, g_0_x_xxxxz_xyyyz, g_0_x_xxxxz_xyyzz, g_0_x_xxxxz_xyzzz, g_0_x_xxxxz_xzzzz, g_0_x_xxxxz_yyyyy, g_0_x_xxxxz_yyyyz, g_0_x_xxxxz_yyyzz, g_0_x_xxxxz_yyzzz, g_0_x_xxxxz_yzzzz, g_0_x_xxxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxz_xxxxx[k] = -g_0_x_xxxx_xxxxx[k] * ab_z + g_0_x_xxxx_xxxxxz[k];

                g_0_x_xxxxz_xxxxy[k] = -g_0_x_xxxx_xxxxy[k] * ab_z + g_0_x_xxxx_xxxxyz[k];

                g_0_x_xxxxz_xxxxz[k] = -g_0_x_xxxx_xxxxz[k] * ab_z + g_0_x_xxxx_xxxxzz[k];

                g_0_x_xxxxz_xxxyy[k] = -g_0_x_xxxx_xxxyy[k] * ab_z + g_0_x_xxxx_xxxyyz[k];

                g_0_x_xxxxz_xxxyz[k] = -g_0_x_xxxx_xxxyz[k] * ab_z + g_0_x_xxxx_xxxyzz[k];

                g_0_x_xxxxz_xxxzz[k] = -g_0_x_xxxx_xxxzz[k] * ab_z + g_0_x_xxxx_xxxzzz[k];

                g_0_x_xxxxz_xxyyy[k] = -g_0_x_xxxx_xxyyy[k] * ab_z + g_0_x_xxxx_xxyyyz[k];

                g_0_x_xxxxz_xxyyz[k] = -g_0_x_xxxx_xxyyz[k] * ab_z + g_0_x_xxxx_xxyyzz[k];

                g_0_x_xxxxz_xxyzz[k] = -g_0_x_xxxx_xxyzz[k] * ab_z + g_0_x_xxxx_xxyzzz[k];

                g_0_x_xxxxz_xxzzz[k] = -g_0_x_xxxx_xxzzz[k] * ab_z + g_0_x_xxxx_xxzzzz[k];

                g_0_x_xxxxz_xyyyy[k] = -g_0_x_xxxx_xyyyy[k] * ab_z + g_0_x_xxxx_xyyyyz[k];

                g_0_x_xxxxz_xyyyz[k] = -g_0_x_xxxx_xyyyz[k] * ab_z + g_0_x_xxxx_xyyyzz[k];

                g_0_x_xxxxz_xyyzz[k] = -g_0_x_xxxx_xyyzz[k] * ab_z + g_0_x_xxxx_xyyzzz[k];

                g_0_x_xxxxz_xyzzz[k] = -g_0_x_xxxx_xyzzz[k] * ab_z + g_0_x_xxxx_xyzzzz[k];

                g_0_x_xxxxz_xzzzz[k] = -g_0_x_xxxx_xzzzz[k] * ab_z + g_0_x_xxxx_xzzzzz[k];

                g_0_x_xxxxz_yyyyy[k] = -g_0_x_xxxx_yyyyy[k] * ab_z + g_0_x_xxxx_yyyyyz[k];

                g_0_x_xxxxz_yyyyz[k] = -g_0_x_xxxx_yyyyz[k] * ab_z + g_0_x_xxxx_yyyyzz[k];

                g_0_x_xxxxz_yyyzz[k] = -g_0_x_xxxx_yyyzz[k] * ab_z + g_0_x_xxxx_yyyzzz[k];

                g_0_x_xxxxz_yyzzz[k] = -g_0_x_xxxx_yyzzz[k] * ab_z + g_0_x_xxxx_yyzzzz[k];

                g_0_x_xxxxz_yzzzz[k] = -g_0_x_xxxx_yzzzz[k] * ab_z + g_0_x_xxxx_yzzzzz[k];

                g_0_x_xxxxz_zzzzz[k] = -g_0_x_xxxx_zzzzz[k] * ab_z + g_0_x_xxxx_zzzzzz[k];
            }

            /// Set up 63-84 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyy_xxxxx = cbuffer.data(hh_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxy = cbuffer.data(hh_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxxz = cbuffer.data(hh_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxyy = cbuffer.data(hh_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxyz = cbuffer.data(hh_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxzz = cbuffer.data(hh_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxyyy = cbuffer.data(hh_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxyyz = cbuffer.data(hh_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxyzz = cbuffer.data(hh_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxzzz = cbuffer.data(hh_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyyyy = cbuffer.data(hh_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyyyz = cbuffer.data(hh_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyyzz = cbuffer.data(hh_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyzzz = cbuffer.data(hh_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxyy_xzzzz = cbuffer.data(hh_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyyyy = cbuffer.data(hh_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyyyz = cbuffer.data(hh_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyyzz = cbuffer.data(hh_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyzzz = cbuffer.data(hh_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxyy_yzzzz = cbuffer.data(hh_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxyy_zzzzz = cbuffer.data(hh_geom_01_off + 83 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxy_xxxxx, g_0_x_xxxy_xxxxxy, g_0_x_xxxy_xxxxy, g_0_x_xxxy_xxxxyy, g_0_x_xxxy_xxxxyz, g_0_x_xxxy_xxxxz, g_0_x_xxxy_xxxyy, g_0_x_xxxy_xxxyyy, g_0_x_xxxy_xxxyyz, g_0_x_xxxy_xxxyz, g_0_x_xxxy_xxxyzz, g_0_x_xxxy_xxxzz, g_0_x_xxxy_xxyyy, g_0_x_xxxy_xxyyyy, g_0_x_xxxy_xxyyyz, g_0_x_xxxy_xxyyz, g_0_x_xxxy_xxyyzz, g_0_x_xxxy_xxyzz, g_0_x_xxxy_xxyzzz, g_0_x_xxxy_xxzzz, g_0_x_xxxy_xyyyy, g_0_x_xxxy_xyyyyy, g_0_x_xxxy_xyyyyz, g_0_x_xxxy_xyyyz, g_0_x_xxxy_xyyyzz, g_0_x_xxxy_xyyzz, g_0_x_xxxy_xyyzzz, g_0_x_xxxy_xyzzz, g_0_x_xxxy_xyzzzz, g_0_x_xxxy_xzzzz, g_0_x_xxxy_yyyyy, g_0_x_xxxy_yyyyyy, g_0_x_xxxy_yyyyyz, g_0_x_xxxy_yyyyz, g_0_x_xxxy_yyyyzz, g_0_x_xxxy_yyyzz, g_0_x_xxxy_yyyzzz, g_0_x_xxxy_yyzzz, g_0_x_xxxy_yyzzzz, g_0_x_xxxy_yzzzz, g_0_x_xxxy_yzzzzz, g_0_x_xxxy_zzzzz, g_0_x_xxxyy_xxxxx, g_0_x_xxxyy_xxxxy, g_0_x_xxxyy_xxxxz, g_0_x_xxxyy_xxxyy, g_0_x_xxxyy_xxxyz, g_0_x_xxxyy_xxxzz, g_0_x_xxxyy_xxyyy, g_0_x_xxxyy_xxyyz, g_0_x_xxxyy_xxyzz, g_0_x_xxxyy_xxzzz, g_0_x_xxxyy_xyyyy, g_0_x_xxxyy_xyyyz, g_0_x_xxxyy_xyyzz, g_0_x_xxxyy_xyzzz, g_0_x_xxxyy_xzzzz, g_0_x_xxxyy_yyyyy, g_0_x_xxxyy_yyyyz, g_0_x_xxxyy_yyyzz, g_0_x_xxxyy_yyzzz, g_0_x_xxxyy_yzzzz, g_0_x_xxxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyy_xxxxx[k] = -g_0_x_xxxy_xxxxx[k] * ab_y + g_0_x_xxxy_xxxxxy[k];

                g_0_x_xxxyy_xxxxy[k] = -g_0_x_xxxy_xxxxy[k] * ab_y + g_0_x_xxxy_xxxxyy[k];

                g_0_x_xxxyy_xxxxz[k] = -g_0_x_xxxy_xxxxz[k] * ab_y + g_0_x_xxxy_xxxxyz[k];

                g_0_x_xxxyy_xxxyy[k] = -g_0_x_xxxy_xxxyy[k] * ab_y + g_0_x_xxxy_xxxyyy[k];

                g_0_x_xxxyy_xxxyz[k] = -g_0_x_xxxy_xxxyz[k] * ab_y + g_0_x_xxxy_xxxyyz[k];

                g_0_x_xxxyy_xxxzz[k] = -g_0_x_xxxy_xxxzz[k] * ab_y + g_0_x_xxxy_xxxyzz[k];

                g_0_x_xxxyy_xxyyy[k] = -g_0_x_xxxy_xxyyy[k] * ab_y + g_0_x_xxxy_xxyyyy[k];

                g_0_x_xxxyy_xxyyz[k] = -g_0_x_xxxy_xxyyz[k] * ab_y + g_0_x_xxxy_xxyyyz[k];

                g_0_x_xxxyy_xxyzz[k] = -g_0_x_xxxy_xxyzz[k] * ab_y + g_0_x_xxxy_xxyyzz[k];

                g_0_x_xxxyy_xxzzz[k] = -g_0_x_xxxy_xxzzz[k] * ab_y + g_0_x_xxxy_xxyzzz[k];

                g_0_x_xxxyy_xyyyy[k] = -g_0_x_xxxy_xyyyy[k] * ab_y + g_0_x_xxxy_xyyyyy[k];

                g_0_x_xxxyy_xyyyz[k] = -g_0_x_xxxy_xyyyz[k] * ab_y + g_0_x_xxxy_xyyyyz[k];

                g_0_x_xxxyy_xyyzz[k] = -g_0_x_xxxy_xyyzz[k] * ab_y + g_0_x_xxxy_xyyyzz[k];

                g_0_x_xxxyy_xyzzz[k] = -g_0_x_xxxy_xyzzz[k] * ab_y + g_0_x_xxxy_xyyzzz[k];

                g_0_x_xxxyy_xzzzz[k] = -g_0_x_xxxy_xzzzz[k] * ab_y + g_0_x_xxxy_xyzzzz[k];

                g_0_x_xxxyy_yyyyy[k] = -g_0_x_xxxy_yyyyy[k] * ab_y + g_0_x_xxxy_yyyyyy[k];

                g_0_x_xxxyy_yyyyz[k] = -g_0_x_xxxy_yyyyz[k] * ab_y + g_0_x_xxxy_yyyyyz[k];

                g_0_x_xxxyy_yyyzz[k] = -g_0_x_xxxy_yyyzz[k] * ab_y + g_0_x_xxxy_yyyyzz[k];

                g_0_x_xxxyy_yyzzz[k] = -g_0_x_xxxy_yyzzz[k] * ab_y + g_0_x_xxxy_yyyzzz[k];

                g_0_x_xxxyy_yzzzz[k] = -g_0_x_xxxy_yzzzz[k] * ab_y + g_0_x_xxxy_yyzzzz[k];

                g_0_x_xxxyy_zzzzz[k] = -g_0_x_xxxy_zzzzz[k] * ab_y + g_0_x_xxxy_yzzzzz[k];
            }

            /// Set up 84-105 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyz_xxxxx = cbuffer.data(hh_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxy = cbuffer.data(hh_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxxz = cbuffer.data(hh_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxyy = cbuffer.data(hh_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxyz = cbuffer.data(hh_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxzz = cbuffer.data(hh_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxyyy = cbuffer.data(hh_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxyyz = cbuffer.data(hh_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxyzz = cbuffer.data(hh_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxzzz = cbuffer.data(hh_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyyyy = cbuffer.data(hh_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyyyz = cbuffer.data(hh_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyyzz = cbuffer.data(hh_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyzzz = cbuffer.data(hh_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxxyz_xzzzz = cbuffer.data(hh_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyyyy = cbuffer.data(hh_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyyyz = cbuffer.data(hh_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyyzz = cbuffer.data(hh_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyzzz = cbuffer.data(hh_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxxyz_yzzzz = cbuffer.data(hh_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxxyz_zzzzz = cbuffer.data(hh_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyz_xxxxx, g_0_x_xxxyz_xxxxy, g_0_x_xxxyz_xxxxz, g_0_x_xxxyz_xxxyy, g_0_x_xxxyz_xxxyz, g_0_x_xxxyz_xxxzz, g_0_x_xxxyz_xxyyy, g_0_x_xxxyz_xxyyz, g_0_x_xxxyz_xxyzz, g_0_x_xxxyz_xxzzz, g_0_x_xxxyz_xyyyy, g_0_x_xxxyz_xyyyz, g_0_x_xxxyz_xyyzz, g_0_x_xxxyz_xyzzz, g_0_x_xxxyz_xzzzz, g_0_x_xxxyz_yyyyy, g_0_x_xxxyz_yyyyz, g_0_x_xxxyz_yyyzz, g_0_x_xxxyz_yyzzz, g_0_x_xxxyz_yzzzz, g_0_x_xxxyz_zzzzz, g_0_x_xxxz_xxxxx, g_0_x_xxxz_xxxxxy, g_0_x_xxxz_xxxxy, g_0_x_xxxz_xxxxyy, g_0_x_xxxz_xxxxyz, g_0_x_xxxz_xxxxz, g_0_x_xxxz_xxxyy, g_0_x_xxxz_xxxyyy, g_0_x_xxxz_xxxyyz, g_0_x_xxxz_xxxyz, g_0_x_xxxz_xxxyzz, g_0_x_xxxz_xxxzz, g_0_x_xxxz_xxyyy, g_0_x_xxxz_xxyyyy, g_0_x_xxxz_xxyyyz, g_0_x_xxxz_xxyyz, g_0_x_xxxz_xxyyzz, g_0_x_xxxz_xxyzz, g_0_x_xxxz_xxyzzz, g_0_x_xxxz_xxzzz, g_0_x_xxxz_xyyyy, g_0_x_xxxz_xyyyyy, g_0_x_xxxz_xyyyyz, g_0_x_xxxz_xyyyz, g_0_x_xxxz_xyyyzz, g_0_x_xxxz_xyyzz, g_0_x_xxxz_xyyzzz, g_0_x_xxxz_xyzzz, g_0_x_xxxz_xyzzzz, g_0_x_xxxz_xzzzz, g_0_x_xxxz_yyyyy, g_0_x_xxxz_yyyyyy, g_0_x_xxxz_yyyyyz, g_0_x_xxxz_yyyyz, g_0_x_xxxz_yyyyzz, g_0_x_xxxz_yyyzz, g_0_x_xxxz_yyyzzz, g_0_x_xxxz_yyzzz, g_0_x_xxxz_yyzzzz, g_0_x_xxxz_yzzzz, g_0_x_xxxz_yzzzzz, g_0_x_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyz_xxxxx[k] = -g_0_x_xxxz_xxxxx[k] * ab_y + g_0_x_xxxz_xxxxxy[k];

                g_0_x_xxxyz_xxxxy[k] = -g_0_x_xxxz_xxxxy[k] * ab_y + g_0_x_xxxz_xxxxyy[k];

                g_0_x_xxxyz_xxxxz[k] = -g_0_x_xxxz_xxxxz[k] * ab_y + g_0_x_xxxz_xxxxyz[k];

                g_0_x_xxxyz_xxxyy[k] = -g_0_x_xxxz_xxxyy[k] * ab_y + g_0_x_xxxz_xxxyyy[k];

                g_0_x_xxxyz_xxxyz[k] = -g_0_x_xxxz_xxxyz[k] * ab_y + g_0_x_xxxz_xxxyyz[k];

                g_0_x_xxxyz_xxxzz[k] = -g_0_x_xxxz_xxxzz[k] * ab_y + g_0_x_xxxz_xxxyzz[k];

                g_0_x_xxxyz_xxyyy[k] = -g_0_x_xxxz_xxyyy[k] * ab_y + g_0_x_xxxz_xxyyyy[k];

                g_0_x_xxxyz_xxyyz[k] = -g_0_x_xxxz_xxyyz[k] * ab_y + g_0_x_xxxz_xxyyyz[k];

                g_0_x_xxxyz_xxyzz[k] = -g_0_x_xxxz_xxyzz[k] * ab_y + g_0_x_xxxz_xxyyzz[k];

                g_0_x_xxxyz_xxzzz[k] = -g_0_x_xxxz_xxzzz[k] * ab_y + g_0_x_xxxz_xxyzzz[k];

                g_0_x_xxxyz_xyyyy[k] = -g_0_x_xxxz_xyyyy[k] * ab_y + g_0_x_xxxz_xyyyyy[k];

                g_0_x_xxxyz_xyyyz[k] = -g_0_x_xxxz_xyyyz[k] * ab_y + g_0_x_xxxz_xyyyyz[k];

                g_0_x_xxxyz_xyyzz[k] = -g_0_x_xxxz_xyyzz[k] * ab_y + g_0_x_xxxz_xyyyzz[k];

                g_0_x_xxxyz_xyzzz[k] = -g_0_x_xxxz_xyzzz[k] * ab_y + g_0_x_xxxz_xyyzzz[k];

                g_0_x_xxxyz_xzzzz[k] = -g_0_x_xxxz_xzzzz[k] * ab_y + g_0_x_xxxz_xyzzzz[k];

                g_0_x_xxxyz_yyyyy[k] = -g_0_x_xxxz_yyyyy[k] * ab_y + g_0_x_xxxz_yyyyyy[k];

                g_0_x_xxxyz_yyyyz[k] = -g_0_x_xxxz_yyyyz[k] * ab_y + g_0_x_xxxz_yyyyyz[k];

                g_0_x_xxxyz_yyyzz[k] = -g_0_x_xxxz_yyyzz[k] * ab_y + g_0_x_xxxz_yyyyzz[k];

                g_0_x_xxxyz_yyzzz[k] = -g_0_x_xxxz_yyzzz[k] * ab_y + g_0_x_xxxz_yyyzzz[k];

                g_0_x_xxxyz_yzzzz[k] = -g_0_x_xxxz_yzzzz[k] * ab_y + g_0_x_xxxz_yyzzzz[k];

                g_0_x_xxxyz_zzzzz[k] = -g_0_x_xxxz_zzzzz[k] * ab_y + g_0_x_xxxz_yzzzzz[k];
            }

            /// Set up 105-126 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzz_xxxxx = cbuffer.data(hh_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxy = cbuffer.data(hh_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxxz = cbuffer.data(hh_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxyy = cbuffer.data(hh_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxyz = cbuffer.data(hh_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxzz = cbuffer.data(hh_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxyyy = cbuffer.data(hh_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxyyz = cbuffer.data(hh_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxyzz = cbuffer.data(hh_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxzzz = cbuffer.data(hh_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyyyy = cbuffer.data(hh_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyyyz = cbuffer.data(hh_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyyzz = cbuffer.data(hh_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyzzz = cbuffer.data(hh_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxxzz_xzzzz = cbuffer.data(hh_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyyyy = cbuffer.data(hh_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyyyz = cbuffer.data(hh_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyyzz = cbuffer.data(hh_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyzzz = cbuffer.data(hh_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxxzz_yzzzz = cbuffer.data(hh_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxxzz_zzzzz = cbuffer.data(hh_geom_01_off + 125 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxz_xxxxx, g_0_x_xxxz_xxxxxz, g_0_x_xxxz_xxxxy, g_0_x_xxxz_xxxxyz, g_0_x_xxxz_xxxxz, g_0_x_xxxz_xxxxzz, g_0_x_xxxz_xxxyy, g_0_x_xxxz_xxxyyz, g_0_x_xxxz_xxxyz, g_0_x_xxxz_xxxyzz, g_0_x_xxxz_xxxzz, g_0_x_xxxz_xxxzzz, g_0_x_xxxz_xxyyy, g_0_x_xxxz_xxyyyz, g_0_x_xxxz_xxyyz, g_0_x_xxxz_xxyyzz, g_0_x_xxxz_xxyzz, g_0_x_xxxz_xxyzzz, g_0_x_xxxz_xxzzz, g_0_x_xxxz_xxzzzz, g_0_x_xxxz_xyyyy, g_0_x_xxxz_xyyyyz, g_0_x_xxxz_xyyyz, g_0_x_xxxz_xyyyzz, g_0_x_xxxz_xyyzz, g_0_x_xxxz_xyyzzz, g_0_x_xxxz_xyzzz, g_0_x_xxxz_xyzzzz, g_0_x_xxxz_xzzzz, g_0_x_xxxz_xzzzzz, g_0_x_xxxz_yyyyy, g_0_x_xxxz_yyyyyz, g_0_x_xxxz_yyyyz, g_0_x_xxxz_yyyyzz, g_0_x_xxxz_yyyzz, g_0_x_xxxz_yyyzzz, g_0_x_xxxz_yyzzz, g_0_x_xxxz_yyzzzz, g_0_x_xxxz_yzzzz, g_0_x_xxxz_yzzzzz, g_0_x_xxxz_zzzzz, g_0_x_xxxz_zzzzzz, g_0_x_xxxzz_xxxxx, g_0_x_xxxzz_xxxxy, g_0_x_xxxzz_xxxxz, g_0_x_xxxzz_xxxyy, g_0_x_xxxzz_xxxyz, g_0_x_xxxzz_xxxzz, g_0_x_xxxzz_xxyyy, g_0_x_xxxzz_xxyyz, g_0_x_xxxzz_xxyzz, g_0_x_xxxzz_xxzzz, g_0_x_xxxzz_xyyyy, g_0_x_xxxzz_xyyyz, g_0_x_xxxzz_xyyzz, g_0_x_xxxzz_xyzzz, g_0_x_xxxzz_xzzzz, g_0_x_xxxzz_yyyyy, g_0_x_xxxzz_yyyyz, g_0_x_xxxzz_yyyzz, g_0_x_xxxzz_yyzzz, g_0_x_xxxzz_yzzzz, g_0_x_xxxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzz_xxxxx[k] = -g_0_x_xxxz_xxxxx[k] * ab_z + g_0_x_xxxz_xxxxxz[k];

                g_0_x_xxxzz_xxxxy[k] = -g_0_x_xxxz_xxxxy[k] * ab_z + g_0_x_xxxz_xxxxyz[k];

                g_0_x_xxxzz_xxxxz[k] = -g_0_x_xxxz_xxxxz[k] * ab_z + g_0_x_xxxz_xxxxzz[k];

                g_0_x_xxxzz_xxxyy[k] = -g_0_x_xxxz_xxxyy[k] * ab_z + g_0_x_xxxz_xxxyyz[k];

                g_0_x_xxxzz_xxxyz[k] = -g_0_x_xxxz_xxxyz[k] * ab_z + g_0_x_xxxz_xxxyzz[k];

                g_0_x_xxxzz_xxxzz[k] = -g_0_x_xxxz_xxxzz[k] * ab_z + g_0_x_xxxz_xxxzzz[k];

                g_0_x_xxxzz_xxyyy[k] = -g_0_x_xxxz_xxyyy[k] * ab_z + g_0_x_xxxz_xxyyyz[k];

                g_0_x_xxxzz_xxyyz[k] = -g_0_x_xxxz_xxyyz[k] * ab_z + g_0_x_xxxz_xxyyzz[k];

                g_0_x_xxxzz_xxyzz[k] = -g_0_x_xxxz_xxyzz[k] * ab_z + g_0_x_xxxz_xxyzzz[k];

                g_0_x_xxxzz_xxzzz[k] = -g_0_x_xxxz_xxzzz[k] * ab_z + g_0_x_xxxz_xxzzzz[k];

                g_0_x_xxxzz_xyyyy[k] = -g_0_x_xxxz_xyyyy[k] * ab_z + g_0_x_xxxz_xyyyyz[k];

                g_0_x_xxxzz_xyyyz[k] = -g_0_x_xxxz_xyyyz[k] * ab_z + g_0_x_xxxz_xyyyzz[k];

                g_0_x_xxxzz_xyyzz[k] = -g_0_x_xxxz_xyyzz[k] * ab_z + g_0_x_xxxz_xyyzzz[k];

                g_0_x_xxxzz_xyzzz[k] = -g_0_x_xxxz_xyzzz[k] * ab_z + g_0_x_xxxz_xyzzzz[k];

                g_0_x_xxxzz_xzzzz[k] = -g_0_x_xxxz_xzzzz[k] * ab_z + g_0_x_xxxz_xzzzzz[k];

                g_0_x_xxxzz_yyyyy[k] = -g_0_x_xxxz_yyyyy[k] * ab_z + g_0_x_xxxz_yyyyyz[k];

                g_0_x_xxxzz_yyyyz[k] = -g_0_x_xxxz_yyyyz[k] * ab_z + g_0_x_xxxz_yyyyzz[k];

                g_0_x_xxxzz_yyyzz[k] = -g_0_x_xxxz_yyyzz[k] * ab_z + g_0_x_xxxz_yyyzzz[k];

                g_0_x_xxxzz_yyzzz[k] = -g_0_x_xxxz_yyzzz[k] * ab_z + g_0_x_xxxz_yyzzzz[k];

                g_0_x_xxxzz_yzzzz[k] = -g_0_x_xxxz_yzzzz[k] * ab_z + g_0_x_xxxz_yzzzzz[k];

                g_0_x_xxxzz_zzzzz[k] = -g_0_x_xxxz_zzzzz[k] * ab_z + g_0_x_xxxz_zzzzzz[k];
            }

            /// Set up 126-147 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyy_xxxxx = cbuffer.data(hh_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxy = cbuffer.data(hh_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxxz = cbuffer.data(hh_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxyy = cbuffer.data(hh_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxyz = cbuffer.data(hh_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxzz = cbuffer.data(hh_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxyyy = cbuffer.data(hh_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxyyz = cbuffer.data(hh_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxyzz = cbuffer.data(hh_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxzzz = cbuffer.data(hh_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyyyy = cbuffer.data(hh_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyyyz = cbuffer.data(hh_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyyzz = cbuffer.data(hh_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyzzz = cbuffer.data(hh_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxyyy_xzzzz = cbuffer.data(hh_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyyyy = cbuffer.data(hh_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyyyz = cbuffer.data(hh_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyyzz = cbuffer.data(hh_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyzzz = cbuffer.data(hh_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxyyy_yzzzz = cbuffer.data(hh_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxyyy_zzzzz = cbuffer.data(hh_geom_01_off + 146 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyy_xxxxx, g_0_x_xxyy_xxxxxy, g_0_x_xxyy_xxxxy, g_0_x_xxyy_xxxxyy, g_0_x_xxyy_xxxxyz, g_0_x_xxyy_xxxxz, g_0_x_xxyy_xxxyy, g_0_x_xxyy_xxxyyy, g_0_x_xxyy_xxxyyz, g_0_x_xxyy_xxxyz, g_0_x_xxyy_xxxyzz, g_0_x_xxyy_xxxzz, g_0_x_xxyy_xxyyy, g_0_x_xxyy_xxyyyy, g_0_x_xxyy_xxyyyz, g_0_x_xxyy_xxyyz, g_0_x_xxyy_xxyyzz, g_0_x_xxyy_xxyzz, g_0_x_xxyy_xxyzzz, g_0_x_xxyy_xxzzz, g_0_x_xxyy_xyyyy, g_0_x_xxyy_xyyyyy, g_0_x_xxyy_xyyyyz, g_0_x_xxyy_xyyyz, g_0_x_xxyy_xyyyzz, g_0_x_xxyy_xyyzz, g_0_x_xxyy_xyyzzz, g_0_x_xxyy_xyzzz, g_0_x_xxyy_xyzzzz, g_0_x_xxyy_xzzzz, g_0_x_xxyy_yyyyy, g_0_x_xxyy_yyyyyy, g_0_x_xxyy_yyyyyz, g_0_x_xxyy_yyyyz, g_0_x_xxyy_yyyyzz, g_0_x_xxyy_yyyzz, g_0_x_xxyy_yyyzzz, g_0_x_xxyy_yyzzz, g_0_x_xxyy_yyzzzz, g_0_x_xxyy_yzzzz, g_0_x_xxyy_yzzzzz, g_0_x_xxyy_zzzzz, g_0_x_xxyyy_xxxxx, g_0_x_xxyyy_xxxxy, g_0_x_xxyyy_xxxxz, g_0_x_xxyyy_xxxyy, g_0_x_xxyyy_xxxyz, g_0_x_xxyyy_xxxzz, g_0_x_xxyyy_xxyyy, g_0_x_xxyyy_xxyyz, g_0_x_xxyyy_xxyzz, g_0_x_xxyyy_xxzzz, g_0_x_xxyyy_xyyyy, g_0_x_xxyyy_xyyyz, g_0_x_xxyyy_xyyzz, g_0_x_xxyyy_xyzzz, g_0_x_xxyyy_xzzzz, g_0_x_xxyyy_yyyyy, g_0_x_xxyyy_yyyyz, g_0_x_xxyyy_yyyzz, g_0_x_xxyyy_yyzzz, g_0_x_xxyyy_yzzzz, g_0_x_xxyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyy_xxxxx[k] = -g_0_x_xxyy_xxxxx[k] * ab_y + g_0_x_xxyy_xxxxxy[k];

                g_0_x_xxyyy_xxxxy[k] = -g_0_x_xxyy_xxxxy[k] * ab_y + g_0_x_xxyy_xxxxyy[k];

                g_0_x_xxyyy_xxxxz[k] = -g_0_x_xxyy_xxxxz[k] * ab_y + g_0_x_xxyy_xxxxyz[k];

                g_0_x_xxyyy_xxxyy[k] = -g_0_x_xxyy_xxxyy[k] * ab_y + g_0_x_xxyy_xxxyyy[k];

                g_0_x_xxyyy_xxxyz[k] = -g_0_x_xxyy_xxxyz[k] * ab_y + g_0_x_xxyy_xxxyyz[k];

                g_0_x_xxyyy_xxxzz[k] = -g_0_x_xxyy_xxxzz[k] * ab_y + g_0_x_xxyy_xxxyzz[k];

                g_0_x_xxyyy_xxyyy[k] = -g_0_x_xxyy_xxyyy[k] * ab_y + g_0_x_xxyy_xxyyyy[k];

                g_0_x_xxyyy_xxyyz[k] = -g_0_x_xxyy_xxyyz[k] * ab_y + g_0_x_xxyy_xxyyyz[k];

                g_0_x_xxyyy_xxyzz[k] = -g_0_x_xxyy_xxyzz[k] * ab_y + g_0_x_xxyy_xxyyzz[k];

                g_0_x_xxyyy_xxzzz[k] = -g_0_x_xxyy_xxzzz[k] * ab_y + g_0_x_xxyy_xxyzzz[k];

                g_0_x_xxyyy_xyyyy[k] = -g_0_x_xxyy_xyyyy[k] * ab_y + g_0_x_xxyy_xyyyyy[k];

                g_0_x_xxyyy_xyyyz[k] = -g_0_x_xxyy_xyyyz[k] * ab_y + g_0_x_xxyy_xyyyyz[k];

                g_0_x_xxyyy_xyyzz[k] = -g_0_x_xxyy_xyyzz[k] * ab_y + g_0_x_xxyy_xyyyzz[k];

                g_0_x_xxyyy_xyzzz[k] = -g_0_x_xxyy_xyzzz[k] * ab_y + g_0_x_xxyy_xyyzzz[k];

                g_0_x_xxyyy_xzzzz[k] = -g_0_x_xxyy_xzzzz[k] * ab_y + g_0_x_xxyy_xyzzzz[k];

                g_0_x_xxyyy_yyyyy[k] = -g_0_x_xxyy_yyyyy[k] * ab_y + g_0_x_xxyy_yyyyyy[k];

                g_0_x_xxyyy_yyyyz[k] = -g_0_x_xxyy_yyyyz[k] * ab_y + g_0_x_xxyy_yyyyyz[k];

                g_0_x_xxyyy_yyyzz[k] = -g_0_x_xxyy_yyyzz[k] * ab_y + g_0_x_xxyy_yyyyzz[k];

                g_0_x_xxyyy_yyzzz[k] = -g_0_x_xxyy_yyzzz[k] * ab_y + g_0_x_xxyy_yyyzzz[k];

                g_0_x_xxyyy_yzzzz[k] = -g_0_x_xxyy_yzzzz[k] * ab_y + g_0_x_xxyy_yyzzzz[k];

                g_0_x_xxyyy_zzzzz[k] = -g_0_x_xxyy_zzzzz[k] * ab_y + g_0_x_xxyy_yzzzzz[k];
            }

            /// Set up 147-168 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyz_xxxxx = cbuffer.data(hh_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxy = cbuffer.data(hh_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxxz = cbuffer.data(hh_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxyy = cbuffer.data(hh_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxyz = cbuffer.data(hh_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxzz = cbuffer.data(hh_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxyyy = cbuffer.data(hh_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxyyz = cbuffer.data(hh_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxyzz = cbuffer.data(hh_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxzzz = cbuffer.data(hh_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyyyy = cbuffer.data(hh_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyyyz = cbuffer.data(hh_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyyzz = cbuffer.data(hh_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyzzz = cbuffer.data(hh_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xxyyz_xzzzz = cbuffer.data(hh_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyyyy = cbuffer.data(hh_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyyyz = cbuffer.data(hh_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyyzz = cbuffer.data(hh_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyzzz = cbuffer.data(hh_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xxyyz_yzzzz = cbuffer.data(hh_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xxyyz_zzzzz = cbuffer.data(hh_geom_01_off + 167 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyz_xxxxx, g_0_x_xxyyz_xxxxy, g_0_x_xxyyz_xxxxz, g_0_x_xxyyz_xxxyy, g_0_x_xxyyz_xxxyz, g_0_x_xxyyz_xxxzz, g_0_x_xxyyz_xxyyy, g_0_x_xxyyz_xxyyz, g_0_x_xxyyz_xxyzz, g_0_x_xxyyz_xxzzz, g_0_x_xxyyz_xyyyy, g_0_x_xxyyz_xyyyz, g_0_x_xxyyz_xyyzz, g_0_x_xxyyz_xyzzz, g_0_x_xxyyz_xzzzz, g_0_x_xxyyz_yyyyy, g_0_x_xxyyz_yyyyz, g_0_x_xxyyz_yyyzz, g_0_x_xxyyz_yyzzz, g_0_x_xxyyz_yzzzz, g_0_x_xxyyz_zzzzz, g_0_x_xxyz_xxxxx, g_0_x_xxyz_xxxxxy, g_0_x_xxyz_xxxxy, g_0_x_xxyz_xxxxyy, g_0_x_xxyz_xxxxyz, g_0_x_xxyz_xxxxz, g_0_x_xxyz_xxxyy, g_0_x_xxyz_xxxyyy, g_0_x_xxyz_xxxyyz, g_0_x_xxyz_xxxyz, g_0_x_xxyz_xxxyzz, g_0_x_xxyz_xxxzz, g_0_x_xxyz_xxyyy, g_0_x_xxyz_xxyyyy, g_0_x_xxyz_xxyyyz, g_0_x_xxyz_xxyyz, g_0_x_xxyz_xxyyzz, g_0_x_xxyz_xxyzz, g_0_x_xxyz_xxyzzz, g_0_x_xxyz_xxzzz, g_0_x_xxyz_xyyyy, g_0_x_xxyz_xyyyyy, g_0_x_xxyz_xyyyyz, g_0_x_xxyz_xyyyz, g_0_x_xxyz_xyyyzz, g_0_x_xxyz_xyyzz, g_0_x_xxyz_xyyzzz, g_0_x_xxyz_xyzzz, g_0_x_xxyz_xyzzzz, g_0_x_xxyz_xzzzz, g_0_x_xxyz_yyyyy, g_0_x_xxyz_yyyyyy, g_0_x_xxyz_yyyyyz, g_0_x_xxyz_yyyyz, g_0_x_xxyz_yyyyzz, g_0_x_xxyz_yyyzz, g_0_x_xxyz_yyyzzz, g_0_x_xxyz_yyzzz, g_0_x_xxyz_yyzzzz, g_0_x_xxyz_yzzzz, g_0_x_xxyz_yzzzzz, g_0_x_xxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyz_xxxxx[k] = -g_0_x_xxyz_xxxxx[k] * ab_y + g_0_x_xxyz_xxxxxy[k];

                g_0_x_xxyyz_xxxxy[k] = -g_0_x_xxyz_xxxxy[k] * ab_y + g_0_x_xxyz_xxxxyy[k];

                g_0_x_xxyyz_xxxxz[k] = -g_0_x_xxyz_xxxxz[k] * ab_y + g_0_x_xxyz_xxxxyz[k];

                g_0_x_xxyyz_xxxyy[k] = -g_0_x_xxyz_xxxyy[k] * ab_y + g_0_x_xxyz_xxxyyy[k];

                g_0_x_xxyyz_xxxyz[k] = -g_0_x_xxyz_xxxyz[k] * ab_y + g_0_x_xxyz_xxxyyz[k];

                g_0_x_xxyyz_xxxzz[k] = -g_0_x_xxyz_xxxzz[k] * ab_y + g_0_x_xxyz_xxxyzz[k];

                g_0_x_xxyyz_xxyyy[k] = -g_0_x_xxyz_xxyyy[k] * ab_y + g_0_x_xxyz_xxyyyy[k];

                g_0_x_xxyyz_xxyyz[k] = -g_0_x_xxyz_xxyyz[k] * ab_y + g_0_x_xxyz_xxyyyz[k];

                g_0_x_xxyyz_xxyzz[k] = -g_0_x_xxyz_xxyzz[k] * ab_y + g_0_x_xxyz_xxyyzz[k];

                g_0_x_xxyyz_xxzzz[k] = -g_0_x_xxyz_xxzzz[k] * ab_y + g_0_x_xxyz_xxyzzz[k];

                g_0_x_xxyyz_xyyyy[k] = -g_0_x_xxyz_xyyyy[k] * ab_y + g_0_x_xxyz_xyyyyy[k];

                g_0_x_xxyyz_xyyyz[k] = -g_0_x_xxyz_xyyyz[k] * ab_y + g_0_x_xxyz_xyyyyz[k];

                g_0_x_xxyyz_xyyzz[k] = -g_0_x_xxyz_xyyzz[k] * ab_y + g_0_x_xxyz_xyyyzz[k];

                g_0_x_xxyyz_xyzzz[k] = -g_0_x_xxyz_xyzzz[k] * ab_y + g_0_x_xxyz_xyyzzz[k];

                g_0_x_xxyyz_xzzzz[k] = -g_0_x_xxyz_xzzzz[k] * ab_y + g_0_x_xxyz_xyzzzz[k];

                g_0_x_xxyyz_yyyyy[k] = -g_0_x_xxyz_yyyyy[k] * ab_y + g_0_x_xxyz_yyyyyy[k];

                g_0_x_xxyyz_yyyyz[k] = -g_0_x_xxyz_yyyyz[k] * ab_y + g_0_x_xxyz_yyyyyz[k];

                g_0_x_xxyyz_yyyzz[k] = -g_0_x_xxyz_yyyzz[k] * ab_y + g_0_x_xxyz_yyyyzz[k];

                g_0_x_xxyyz_yyzzz[k] = -g_0_x_xxyz_yyzzz[k] * ab_y + g_0_x_xxyz_yyyzzz[k];

                g_0_x_xxyyz_yzzzz[k] = -g_0_x_xxyz_yzzzz[k] * ab_y + g_0_x_xxyz_yyzzzz[k];

                g_0_x_xxyyz_zzzzz[k] = -g_0_x_xxyz_zzzzz[k] * ab_y + g_0_x_xxyz_yzzzzz[k];
            }

            /// Set up 168-189 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzz_xxxxx = cbuffer.data(hh_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxy = cbuffer.data(hh_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxxz = cbuffer.data(hh_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxyy = cbuffer.data(hh_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxyz = cbuffer.data(hh_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxzz = cbuffer.data(hh_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxyyy = cbuffer.data(hh_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxyyz = cbuffer.data(hh_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxyzz = cbuffer.data(hh_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxzzz = cbuffer.data(hh_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyyyy = cbuffer.data(hh_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyyyz = cbuffer.data(hh_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyyzz = cbuffer.data(hh_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyzzz = cbuffer.data(hh_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xxyzz_xzzzz = cbuffer.data(hh_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyyyy = cbuffer.data(hh_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyyyz = cbuffer.data(hh_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyyzz = cbuffer.data(hh_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyzzz = cbuffer.data(hh_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xxyzz_yzzzz = cbuffer.data(hh_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xxyzz_zzzzz = cbuffer.data(hh_geom_01_off + 188 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzz_xxxxx, g_0_x_xxyzz_xxxxy, g_0_x_xxyzz_xxxxz, g_0_x_xxyzz_xxxyy, g_0_x_xxyzz_xxxyz, g_0_x_xxyzz_xxxzz, g_0_x_xxyzz_xxyyy, g_0_x_xxyzz_xxyyz, g_0_x_xxyzz_xxyzz, g_0_x_xxyzz_xxzzz, g_0_x_xxyzz_xyyyy, g_0_x_xxyzz_xyyyz, g_0_x_xxyzz_xyyzz, g_0_x_xxyzz_xyzzz, g_0_x_xxyzz_xzzzz, g_0_x_xxyzz_yyyyy, g_0_x_xxyzz_yyyyz, g_0_x_xxyzz_yyyzz, g_0_x_xxyzz_yyzzz, g_0_x_xxyzz_yzzzz, g_0_x_xxyzz_zzzzz, g_0_x_xxzz_xxxxx, g_0_x_xxzz_xxxxxy, g_0_x_xxzz_xxxxy, g_0_x_xxzz_xxxxyy, g_0_x_xxzz_xxxxyz, g_0_x_xxzz_xxxxz, g_0_x_xxzz_xxxyy, g_0_x_xxzz_xxxyyy, g_0_x_xxzz_xxxyyz, g_0_x_xxzz_xxxyz, g_0_x_xxzz_xxxyzz, g_0_x_xxzz_xxxzz, g_0_x_xxzz_xxyyy, g_0_x_xxzz_xxyyyy, g_0_x_xxzz_xxyyyz, g_0_x_xxzz_xxyyz, g_0_x_xxzz_xxyyzz, g_0_x_xxzz_xxyzz, g_0_x_xxzz_xxyzzz, g_0_x_xxzz_xxzzz, g_0_x_xxzz_xyyyy, g_0_x_xxzz_xyyyyy, g_0_x_xxzz_xyyyyz, g_0_x_xxzz_xyyyz, g_0_x_xxzz_xyyyzz, g_0_x_xxzz_xyyzz, g_0_x_xxzz_xyyzzz, g_0_x_xxzz_xyzzz, g_0_x_xxzz_xyzzzz, g_0_x_xxzz_xzzzz, g_0_x_xxzz_yyyyy, g_0_x_xxzz_yyyyyy, g_0_x_xxzz_yyyyyz, g_0_x_xxzz_yyyyz, g_0_x_xxzz_yyyyzz, g_0_x_xxzz_yyyzz, g_0_x_xxzz_yyyzzz, g_0_x_xxzz_yyzzz, g_0_x_xxzz_yyzzzz, g_0_x_xxzz_yzzzz, g_0_x_xxzz_yzzzzz, g_0_x_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzz_xxxxx[k] = -g_0_x_xxzz_xxxxx[k] * ab_y + g_0_x_xxzz_xxxxxy[k];

                g_0_x_xxyzz_xxxxy[k] = -g_0_x_xxzz_xxxxy[k] * ab_y + g_0_x_xxzz_xxxxyy[k];

                g_0_x_xxyzz_xxxxz[k] = -g_0_x_xxzz_xxxxz[k] * ab_y + g_0_x_xxzz_xxxxyz[k];

                g_0_x_xxyzz_xxxyy[k] = -g_0_x_xxzz_xxxyy[k] * ab_y + g_0_x_xxzz_xxxyyy[k];

                g_0_x_xxyzz_xxxyz[k] = -g_0_x_xxzz_xxxyz[k] * ab_y + g_0_x_xxzz_xxxyyz[k];

                g_0_x_xxyzz_xxxzz[k] = -g_0_x_xxzz_xxxzz[k] * ab_y + g_0_x_xxzz_xxxyzz[k];

                g_0_x_xxyzz_xxyyy[k] = -g_0_x_xxzz_xxyyy[k] * ab_y + g_0_x_xxzz_xxyyyy[k];

                g_0_x_xxyzz_xxyyz[k] = -g_0_x_xxzz_xxyyz[k] * ab_y + g_0_x_xxzz_xxyyyz[k];

                g_0_x_xxyzz_xxyzz[k] = -g_0_x_xxzz_xxyzz[k] * ab_y + g_0_x_xxzz_xxyyzz[k];

                g_0_x_xxyzz_xxzzz[k] = -g_0_x_xxzz_xxzzz[k] * ab_y + g_0_x_xxzz_xxyzzz[k];

                g_0_x_xxyzz_xyyyy[k] = -g_0_x_xxzz_xyyyy[k] * ab_y + g_0_x_xxzz_xyyyyy[k];

                g_0_x_xxyzz_xyyyz[k] = -g_0_x_xxzz_xyyyz[k] * ab_y + g_0_x_xxzz_xyyyyz[k];

                g_0_x_xxyzz_xyyzz[k] = -g_0_x_xxzz_xyyzz[k] * ab_y + g_0_x_xxzz_xyyyzz[k];

                g_0_x_xxyzz_xyzzz[k] = -g_0_x_xxzz_xyzzz[k] * ab_y + g_0_x_xxzz_xyyzzz[k];

                g_0_x_xxyzz_xzzzz[k] = -g_0_x_xxzz_xzzzz[k] * ab_y + g_0_x_xxzz_xyzzzz[k];

                g_0_x_xxyzz_yyyyy[k] = -g_0_x_xxzz_yyyyy[k] * ab_y + g_0_x_xxzz_yyyyyy[k];

                g_0_x_xxyzz_yyyyz[k] = -g_0_x_xxzz_yyyyz[k] * ab_y + g_0_x_xxzz_yyyyyz[k];

                g_0_x_xxyzz_yyyzz[k] = -g_0_x_xxzz_yyyzz[k] * ab_y + g_0_x_xxzz_yyyyzz[k];

                g_0_x_xxyzz_yyzzz[k] = -g_0_x_xxzz_yyzzz[k] * ab_y + g_0_x_xxzz_yyyzzz[k];

                g_0_x_xxyzz_yzzzz[k] = -g_0_x_xxzz_yzzzz[k] * ab_y + g_0_x_xxzz_yyzzzz[k];

                g_0_x_xxyzz_zzzzz[k] = -g_0_x_xxzz_zzzzz[k] * ab_y + g_0_x_xxzz_yzzzzz[k];
            }

            /// Set up 189-210 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzz_xxxxx = cbuffer.data(hh_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxy = cbuffer.data(hh_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxxz = cbuffer.data(hh_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxyy = cbuffer.data(hh_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxyz = cbuffer.data(hh_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxzz = cbuffer.data(hh_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxyyy = cbuffer.data(hh_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxyyz = cbuffer.data(hh_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxyzz = cbuffer.data(hh_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxzzz = cbuffer.data(hh_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyyyy = cbuffer.data(hh_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyyyz = cbuffer.data(hh_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyyzz = cbuffer.data(hh_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyzzz = cbuffer.data(hh_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xxzzz_xzzzz = cbuffer.data(hh_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyyyy = cbuffer.data(hh_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyyyz = cbuffer.data(hh_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyyzz = cbuffer.data(hh_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyzzz = cbuffer.data(hh_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xxzzz_yzzzz = cbuffer.data(hh_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xxzzz_zzzzz = cbuffer.data(hh_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzz_xxxxx, g_0_x_xxzz_xxxxxz, g_0_x_xxzz_xxxxy, g_0_x_xxzz_xxxxyz, g_0_x_xxzz_xxxxz, g_0_x_xxzz_xxxxzz, g_0_x_xxzz_xxxyy, g_0_x_xxzz_xxxyyz, g_0_x_xxzz_xxxyz, g_0_x_xxzz_xxxyzz, g_0_x_xxzz_xxxzz, g_0_x_xxzz_xxxzzz, g_0_x_xxzz_xxyyy, g_0_x_xxzz_xxyyyz, g_0_x_xxzz_xxyyz, g_0_x_xxzz_xxyyzz, g_0_x_xxzz_xxyzz, g_0_x_xxzz_xxyzzz, g_0_x_xxzz_xxzzz, g_0_x_xxzz_xxzzzz, g_0_x_xxzz_xyyyy, g_0_x_xxzz_xyyyyz, g_0_x_xxzz_xyyyz, g_0_x_xxzz_xyyyzz, g_0_x_xxzz_xyyzz, g_0_x_xxzz_xyyzzz, g_0_x_xxzz_xyzzz, g_0_x_xxzz_xyzzzz, g_0_x_xxzz_xzzzz, g_0_x_xxzz_xzzzzz, g_0_x_xxzz_yyyyy, g_0_x_xxzz_yyyyyz, g_0_x_xxzz_yyyyz, g_0_x_xxzz_yyyyzz, g_0_x_xxzz_yyyzz, g_0_x_xxzz_yyyzzz, g_0_x_xxzz_yyzzz, g_0_x_xxzz_yyzzzz, g_0_x_xxzz_yzzzz, g_0_x_xxzz_yzzzzz, g_0_x_xxzz_zzzzz, g_0_x_xxzz_zzzzzz, g_0_x_xxzzz_xxxxx, g_0_x_xxzzz_xxxxy, g_0_x_xxzzz_xxxxz, g_0_x_xxzzz_xxxyy, g_0_x_xxzzz_xxxyz, g_0_x_xxzzz_xxxzz, g_0_x_xxzzz_xxyyy, g_0_x_xxzzz_xxyyz, g_0_x_xxzzz_xxyzz, g_0_x_xxzzz_xxzzz, g_0_x_xxzzz_xyyyy, g_0_x_xxzzz_xyyyz, g_0_x_xxzzz_xyyzz, g_0_x_xxzzz_xyzzz, g_0_x_xxzzz_xzzzz, g_0_x_xxzzz_yyyyy, g_0_x_xxzzz_yyyyz, g_0_x_xxzzz_yyyzz, g_0_x_xxzzz_yyzzz, g_0_x_xxzzz_yzzzz, g_0_x_xxzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzz_xxxxx[k] = -g_0_x_xxzz_xxxxx[k] * ab_z + g_0_x_xxzz_xxxxxz[k];

                g_0_x_xxzzz_xxxxy[k] = -g_0_x_xxzz_xxxxy[k] * ab_z + g_0_x_xxzz_xxxxyz[k];

                g_0_x_xxzzz_xxxxz[k] = -g_0_x_xxzz_xxxxz[k] * ab_z + g_0_x_xxzz_xxxxzz[k];

                g_0_x_xxzzz_xxxyy[k] = -g_0_x_xxzz_xxxyy[k] * ab_z + g_0_x_xxzz_xxxyyz[k];

                g_0_x_xxzzz_xxxyz[k] = -g_0_x_xxzz_xxxyz[k] * ab_z + g_0_x_xxzz_xxxyzz[k];

                g_0_x_xxzzz_xxxzz[k] = -g_0_x_xxzz_xxxzz[k] * ab_z + g_0_x_xxzz_xxxzzz[k];

                g_0_x_xxzzz_xxyyy[k] = -g_0_x_xxzz_xxyyy[k] * ab_z + g_0_x_xxzz_xxyyyz[k];

                g_0_x_xxzzz_xxyyz[k] = -g_0_x_xxzz_xxyyz[k] * ab_z + g_0_x_xxzz_xxyyzz[k];

                g_0_x_xxzzz_xxyzz[k] = -g_0_x_xxzz_xxyzz[k] * ab_z + g_0_x_xxzz_xxyzzz[k];

                g_0_x_xxzzz_xxzzz[k] = -g_0_x_xxzz_xxzzz[k] * ab_z + g_0_x_xxzz_xxzzzz[k];

                g_0_x_xxzzz_xyyyy[k] = -g_0_x_xxzz_xyyyy[k] * ab_z + g_0_x_xxzz_xyyyyz[k];

                g_0_x_xxzzz_xyyyz[k] = -g_0_x_xxzz_xyyyz[k] * ab_z + g_0_x_xxzz_xyyyzz[k];

                g_0_x_xxzzz_xyyzz[k] = -g_0_x_xxzz_xyyzz[k] * ab_z + g_0_x_xxzz_xyyzzz[k];

                g_0_x_xxzzz_xyzzz[k] = -g_0_x_xxzz_xyzzz[k] * ab_z + g_0_x_xxzz_xyzzzz[k];

                g_0_x_xxzzz_xzzzz[k] = -g_0_x_xxzz_xzzzz[k] * ab_z + g_0_x_xxzz_xzzzzz[k];

                g_0_x_xxzzz_yyyyy[k] = -g_0_x_xxzz_yyyyy[k] * ab_z + g_0_x_xxzz_yyyyyz[k];

                g_0_x_xxzzz_yyyyz[k] = -g_0_x_xxzz_yyyyz[k] * ab_z + g_0_x_xxzz_yyyyzz[k];

                g_0_x_xxzzz_yyyzz[k] = -g_0_x_xxzz_yyyzz[k] * ab_z + g_0_x_xxzz_yyyzzz[k];

                g_0_x_xxzzz_yyzzz[k] = -g_0_x_xxzz_yyzzz[k] * ab_z + g_0_x_xxzz_yyzzzz[k];

                g_0_x_xxzzz_yzzzz[k] = -g_0_x_xxzz_yzzzz[k] * ab_z + g_0_x_xxzz_yzzzzz[k];

                g_0_x_xxzzz_zzzzz[k] = -g_0_x_xxzz_zzzzz[k] * ab_z + g_0_x_xxzz_zzzzzz[k];
            }

            /// Set up 210-231 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyy_xxxxx = cbuffer.data(hh_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxy = cbuffer.data(hh_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxxz = cbuffer.data(hh_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxyy = cbuffer.data(hh_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxyz = cbuffer.data(hh_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxzz = cbuffer.data(hh_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxyyy = cbuffer.data(hh_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxyyz = cbuffer.data(hh_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxyzz = cbuffer.data(hh_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxzzz = cbuffer.data(hh_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyyyy = cbuffer.data(hh_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyyyz = cbuffer.data(hh_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyyzz = cbuffer.data(hh_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyzzz = cbuffer.data(hh_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xyyyy_xzzzz = cbuffer.data(hh_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyyyy = cbuffer.data(hh_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyyyz = cbuffer.data(hh_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyyzz = cbuffer.data(hh_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyzzz = cbuffer.data(hh_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_xyyyy_yzzzz = cbuffer.data(hh_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_xyyyy_zzzzz = cbuffer.data(hh_geom_01_off + 230 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyy_xxxxx, g_0_x_xyyy_xxxxxy, g_0_x_xyyy_xxxxy, g_0_x_xyyy_xxxxyy, g_0_x_xyyy_xxxxyz, g_0_x_xyyy_xxxxz, g_0_x_xyyy_xxxyy, g_0_x_xyyy_xxxyyy, g_0_x_xyyy_xxxyyz, g_0_x_xyyy_xxxyz, g_0_x_xyyy_xxxyzz, g_0_x_xyyy_xxxzz, g_0_x_xyyy_xxyyy, g_0_x_xyyy_xxyyyy, g_0_x_xyyy_xxyyyz, g_0_x_xyyy_xxyyz, g_0_x_xyyy_xxyyzz, g_0_x_xyyy_xxyzz, g_0_x_xyyy_xxyzzz, g_0_x_xyyy_xxzzz, g_0_x_xyyy_xyyyy, g_0_x_xyyy_xyyyyy, g_0_x_xyyy_xyyyyz, g_0_x_xyyy_xyyyz, g_0_x_xyyy_xyyyzz, g_0_x_xyyy_xyyzz, g_0_x_xyyy_xyyzzz, g_0_x_xyyy_xyzzz, g_0_x_xyyy_xyzzzz, g_0_x_xyyy_xzzzz, g_0_x_xyyy_yyyyy, g_0_x_xyyy_yyyyyy, g_0_x_xyyy_yyyyyz, g_0_x_xyyy_yyyyz, g_0_x_xyyy_yyyyzz, g_0_x_xyyy_yyyzz, g_0_x_xyyy_yyyzzz, g_0_x_xyyy_yyzzz, g_0_x_xyyy_yyzzzz, g_0_x_xyyy_yzzzz, g_0_x_xyyy_yzzzzz, g_0_x_xyyy_zzzzz, g_0_x_xyyyy_xxxxx, g_0_x_xyyyy_xxxxy, g_0_x_xyyyy_xxxxz, g_0_x_xyyyy_xxxyy, g_0_x_xyyyy_xxxyz, g_0_x_xyyyy_xxxzz, g_0_x_xyyyy_xxyyy, g_0_x_xyyyy_xxyyz, g_0_x_xyyyy_xxyzz, g_0_x_xyyyy_xxzzz, g_0_x_xyyyy_xyyyy, g_0_x_xyyyy_xyyyz, g_0_x_xyyyy_xyyzz, g_0_x_xyyyy_xyzzz, g_0_x_xyyyy_xzzzz, g_0_x_xyyyy_yyyyy, g_0_x_xyyyy_yyyyz, g_0_x_xyyyy_yyyzz, g_0_x_xyyyy_yyzzz, g_0_x_xyyyy_yzzzz, g_0_x_xyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyy_xxxxx[k] = -g_0_x_xyyy_xxxxx[k] * ab_y + g_0_x_xyyy_xxxxxy[k];

                g_0_x_xyyyy_xxxxy[k] = -g_0_x_xyyy_xxxxy[k] * ab_y + g_0_x_xyyy_xxxxyy[k];

                g_0_x_xyyyy_xxxxz[k] = -g_0_x_xyyy_xxxxz[k] * ab_y + g_0_x_xyyy_xxxxyz[k];

                g_0_x_xyyyy_xxxyy[k] = -g_0_x_xyyy_xxxyy[k] * ab_y + g_0_x_xyyy_xxxyyy[k];

                g_0_x_xyyyy_xxxyz[k] = -g_0_x_xyyy_xxxyz[k] * ab_y + g_0_x_xyyy_xxxyyz[k];

                g_0_x_xyyyy_xxxzz[k] = -g_0_x_xyyy_xxxzz[k] * ab_y + g_0_x_xyyy_xxxyzz[k];

                g_0_x_xyyyy_xxyyy[k] = -g_0_x_xyyy_xxyyy[k] * ab_y + g_0_x_xyyy_xxyyyy[k];

                g_0_x_xyyyy_xxyyz[k] = -g_0_x_xyyy_xxyyz[k] * ab_y + g_0_x_xyyy_xxyyyz[k];

                g_0_x_xyyyy_xxyzz[k] = -g_0_x_xyyy_xxyzz[k] * ab_y + g_0_x_xyyy_xxyyzz[k];

                g_0_x_xyyyy_xxzzz[k] = -g_0_x_xyyy_xxzzz[k] * ab_y + g_0_x_xyyy_xxyzzz[k];

                g_0_x_xyyyy_xyyyy[k] = -g_0_x_xyyy_xyyyy[k] * ab_y + g_0_x_xyyy_xyyyyy[k];

                g_0_x_xyyyy_xyyyz[k] = -g_0_x_xyyy_xyyyz[k] * ab_y + g_0_x_xyyy_xyyyyz[k];

                g_0_x_xyyyy_xyyzz[k] = -g_0_x_xyyy_xyyzz[k] * ab_y + g_0_x_xyyy_xyyyzz[k];

                g_0_x_xyyyy_xyzzz[k] = -g_0_x_xyyy_xyzzz[k] * ab_y + g_0_x_xyyy_xyyzzz[k];

                g_0_x_xyyyy_xzzzz[k] = -g_0_x_xyyy_xzzzz[k] * ab_y + g_0_x_xyyy_xyzzzz[k];

                g_0_x_xyyyy_yyyyy[k] = -g_0_x_xyyy_yyyyy[k] * ab_y + g_0_x_xyyy_yyyyyy[k];

                g_0_x_xyyyy_yyyyz[k] = -g_0_x_xyyy_yyyyz[k] * ab_y + g_0_x_xyyy_yyyyyz[k];

                g_0_x_xyyyy_yyyzz[k] = -g_0_x_xyyy_yyyzz[k] * ab_y + g_0_x_xyyy_yyyyzz[k];

                g_0_x_xyyyy_yyzzz[k] = -g_0_x_xyyy_yyzzz[k] * ab_y + g_0_x_xyyy_yyyzzz[k];

                g_0_x_xyyyy_yzzzz[k] = -g_0_x_xyyy_yzzzz[k] * ab_y + g_0_x_xyyy_yyzzzz[k];

                g_0_x_xyyyy_zzzzz[k] = -g_0_x_xyyy_zzzzz[k] * ab_y + g_0_x_xyyy_yzzzzz[k];
            }

            /// Set up 231-252 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyz_xxxxx = cbuffer.data(hh_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxy = cbuffer.data(hh_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxxz = cbuffer.data(hh_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxyy = cbuffer.data(hh_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxyz = cbuffer.data(hh_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxzz = cbuffer.data(hh_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxyyy = cbuffer.data(hh_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxyyz = cbuffer.data(hh_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxyzz = cbuffer.data(hh_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxzzz = cbuffer.data(hh_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyyyy = cbuffer.data(hh_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyyyz = cbuffer.data(hh_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyyzz = cbuffer.data(hh_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyzzz = cbuffer.data(hh_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_xyyyz_xzzzz = cbuffer.data(hh_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyyyy = cbuffer.data(hh_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyyyz = cbuffer.data(hh_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyyzz = cbuffer.data(hh_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyzzz = cbuffer.data(hh_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_xyyyz_yzzzz = cbuffer.data(hh_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_xyyyz_zzzzz = cbuffer.data(hh_geom_01_off + 251 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyz_xxxxx, g_0_x_xyyyz_xxxxy, g_0_x_xyyyz_xxxxz, g_0_x_xyyyz_xxxyy, g_0_x_xyyyz_xxxyz, g_0_x_xyyyz_xxxzz, g_0_x_xyyyz_xxyyy, g_0_x_xyyyz_xxyyz, g_0_x_xyyyz_xxyzz, g_0_x_xyyyz_xxzzz, g_0_x_xyyyz_xyyyy, g_0_x_xyyyz_xyyyz, g_0_x_xyyyz_xyyzz, g_0_x_xyyyz_xyzzz, g_0_x_xyyyz_xzzzz, g_0_x_xyyyz_yyyyy, g_0_x_xyyyz_yyyyz, g_0_x_xyyyz_yyyzz, g_0_x_xyyyz_yyzzz, g_0_x_xyyyz_yzzzz, g_0_x_xyyyz_zzzzz, g_0_x_xyyz_xxxxx, g_0_x_xyyz_xxxxxy, g_0_x_xyyz_xxxxy, g_0_x_xyyz_xxxxyy, g_0_x_xyyz_xxxxyz, g_0_x_xyyz_xxxxz, g_0_x_xyyz_xxxyy, g_0_x_xyyz_xxxyyy, g_0_x_xyyz_xxxyyz, g_0_x_xyyz_xxxyz, g_0_x_xyyz_xxxyzz, g_0_x_xyyz_xxxzz, g_0_x_xyyz_xxyyy, g_0_x_xyyz_xxyyyy, g_0_x_xyyz_xxyyyz, g_0_x_xyyz_xxyyz, g_0_x_xyyz_xxyyzz, g_0_x_xyyz_xxyzz, g_0_x_xyyz_xxyzzz, g_0_x_xyyz_xxzzz, g_0_x_xyyz_xyyyy, g_0_x_xyyz_xyyyyy, g_0_x_xyyz_xyyyyz, g_0_x_xyyz_xyyyz, g_0_x_xyyz_xyyyzz, g_0_x_xyyz_xyyzz, g_0_x_xyyz_xyyzzz, g_0_x_xyyz_xyzzz, g_0_x_xyyz_xyzzzz, g_0_x_xyyz_xzzzz, g_0_x_xyyz_yyyyy, g_0_x_xyyz_yyyyyy, g_0_x_xyyz_yyyyyz, g_0_x_xyyz_yyyyz, g_0_x_xyyz_yyyyzz, g_0_x_xyyz_yyyzz, g_0_x_xyyz_yyyzzz, g_0_x_xyyz_yyzzz, g_0_x_xyyz_yyzzzz, g_0_x_xyyz_yzzzz, g_0_x_xyyz_yzzzzz, g_0_x_xyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyz_xxxxx[k] = -g_0_x_xyyz_xxxxx[k] * ab_y + g_0_x_xyyz_xxxxxy[k];

                g_0_x_xyyyz_xxxxy[k] = -g_0_x_xyyz_xxxxy[k] * ab_y + g_0_x_xyyz_xxxxyy[k];

                g_0_x_xyyyz_xxxxz[k] = -g_0_x_xyyz_xxxxz[k] * ab_y + g_0_x_xyyz_xxxxyz[k];

                g_0_x_xyyyz_xxxyy[k] = -g_0_x_xyyz_xxxyy[k] * ab_y + g_0_x_xyyz_xxxyyy[k];

                g_0_x_xyyyz_xxxyz[k] = -g_0_x_xyyz_xxxyz[k] * ab_y + g_0_x_xyyz_xxxyyz[k];

                g_0_x_xyyyz_xxxzz[k] = -g_0_x_xyyz_xxxzz[k] * ab_y + g_0_x_xyyz_xxxyzz[k];

                g_0_x_xyyyz_xxyyy[k] = -g_0_x_xyyz_xxyyy[k] * ab_y + g_0_x_xyyz_xxyyyy[k];

                g_0_x_xyyyz_xxyyz[k] = -g_0_x_xyyz_xxyyz[k] * ab_y + g_0_x_xyyz_xxyyyz[k];

                g_0_x_xyyyz_xxyzz[k] = -g_0_x_xyyz_xxyzz[k] * ab_y + g_0_x_xyyz_xxyyzz[k];

                g_0_x_xyyyz_xxzzz[k] = -g_0_x_xyyz_xxzzz[k] * ab_y + g_0_x_xyyz_xxyzzz[k];

                g_0_x_xyyyz_xyyyy[k] = -g_0_x_xyyz_xyyyy[k] * ab_y + g_0_x_xyyz_xyyyyy[k];

                g_0_x_xyyyz_xyyyz[k] = -g_0_x_xyyz_xyyyz[k] * ab_y + g_0_x_xyyz_xyyyyz[k];

                g_0_x_xyyyz_xyyzz[k] = -g_0_x_xyyz_xyyzz[k] * ab_y + g_0_x_xyyz_xyyyzz[k];

                g_0_x_xyyyz_xyzzz[k] = -g_0_x_xyyz_xyzzz[k] * ab_y + g_0_x_xyyz_xyyzzz[k];

                g_0_x_xyyyz_xzzzz[k] = -g_0_x_xyyz_xzzzz[k] * ab_y + g_0_x_xyyz_xyzzzz[k];

                g_0_x_xyyyz_yyyyy[k] = -g_0_x_xyyz_yyyyy[k] * ab_y + g_0_x_xyyz_yyyyyy[k];

                g_0_x_xyyyz_yyyyz[k] = -g_0_x_xyyz_yyyyz[k] * ab_y + g_0_x_xyyz_yyyyyz[k];

                g_0_x_xyyyz_yyyzz[k] = -g_0_x_xyyz_yyyzz[k] * ab_y + g_0_x_xyyz_yyyyzz[k];

                g_0_x_xyyyz_yyzzz[k] = -g_0_x_xyyz_yyzzz[k] * ab_y + g_0_x_xyyz_yyyzzz[k];

                g_0_x_xyyyz_yzzzz[k] = -g_0_x_xyyz_yzzzz[k] * ab_y + g_0_x_xyyz_yyzzzz[k];

                g_0_x_xyyyz_zzzzz[k] = -g_0_x_xyyz_zzzzz[k] * ab_y + g_0_x_xyyz_yzzzzz[k];
            }

            /// Set up 252-273 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzz_xxxxx = cbuffer.data(hh_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxy = cbuffer.data(hh_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxxz = cbuffer.data(hh_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxyy = cbuffer.data(hh_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxyz = cbuffer.data(hh_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxzz = cbuffer.data(hh_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxyyy = cbuffer.data(hh_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxyyz = cbuffer.data(hh_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxyzz = cbuffer.data(hh_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxzzz = cbuffer.data(hh_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyyyy = cbuffer.data(hh_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyyyz = cbuffer.data(hh_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyyzz = cbuffer.data(hh_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyzzz = cbuffer.data(hh_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_xyyzz_xzzzz = cbuffer.data(hh_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyyyy = cbuffer.data(hh_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyyyz = cbuffer.data(hh_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyyzz = cbuffer.data(hh_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyzzz = cbuffer.data(hh_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_xyyzz_yzzzz = cbuffer.data(hh_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_xyyzz_zzzzz = cbuffer.data(hh_geom_01_off + 272 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzz_xxxxx, g_0_x_xyyzz_xxxxy, g_0_x_xyyzz_xxxxz, g_0_x_xyyzz_xxxyy, g_0_x_xyyzz_xxxyz, g_0_x_xyyzz_xxxzz, g_0_x_xyyzz_xxyyy, g_0_x_xyyzz_xxyyz, g_0_x_xyyzz_xxyzz, g_0_x_xyyzz_xxzzz, g_0_x_xyyzz_xyyyy, g_0_x_xyyzz_xyyyz, g_0_x_xyyzz_xyyzz, g_0_x_xyyzz_xyzzz, g_0_x_xyyzz_xzzzz, g_0_x_xyyzz_yyyyy, g_0_x_xyyzz_yyyyz, g_0_x_xyyzz_yyyzz, g_0_x_xyyzz_yyzzz, g_0_x_xyyzz_yzzzz, g_0_x_xyyzz_zzzzz, g_0_x_xyzz_xxxxx, g_0_x_xyzz_xxxxxy, g_0_x_xyzz_xxxxy, g_0_x_xyzz_xxxxyy, g_0_x_xyzz_xxxxyz, g_0_x_xyzz_xxxxz, g_0_x_xyzz_xxxyy, g_0_x_xyzz_xxxyyy, g_0_x_xyzz_xxxyyz, g_0_x_xyzz_xxxyz, g_0_x_xyzz_xxxyzz, g_0_x_xyzz_xxxzz, g_0_x_xyzz_xxyyy, g_0_x_xyzz_xxyyyy, g_0_x_xyzz_xxyyyz, g_0_x_xyzz_xxyyz, g_0_x_xyzz_xxyyzz, g_0_x_xyzz_xxyzz, g_0_x_xyzz_xxyzzz, g_0_x_xyzz_xxzzz, g_0_x_xyzz_xyyyy, g_0_x_xyzz_xyyyyy, g_0_x_xyzz_xyyyyz, g_0_x_xyzz_xyyyz, g_0_x_xyzz_xyyyzz, g_0_x_xyzz_xyyzz, g_0_x_xyzz_xyyzzz, g_0_x_xyzz_xyzzz, g_0_x_xyzz_xyzzzz, g_0_x_xyzz_xzzzz, g_0_x_xyzz_yyyyy, g_0_x_xyzz_yyyyyy, g_0_x_xyzz_yyyyyz, g_0_x_xyzz_yyyyz, g_0_x_xyzz_yyyyzz, g_0_x_xyzz_yyyzz, g_0_x_xyzz_yyyzzz, g_0_x_xyzz_yyzzz, g_0_x_xyzz_yyzzzz, g_0_x_xyzz_yzzzz, g_0_x_xyzz_yzzzzz, g_0_x_xyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzz_xxxxx[k] = -g_0_x_xyzz_xxxxx[k] * ab_y + g_0_x_xyzz_xxxxxy[k];

                g_0_x_xyyzz_xxxxy[k] = -g_0_x_xyzz_xxxxy[k] * ab_y + g_0_x_xyzz_xxxxyy[k];

                g_0_x_xyyzz_xxxxz[k] = -g_0_x_xyzz_xxxxz[k] * ab_y + g_0_x_xyzz_xxxxyz[k];

                g_0_x_xyyzz_xxxyy[k] = -g_0_x_xyzz_xxxyy[k] * ab_y + g_0_x_xyzz_xxxyyy[k];

                g_0_x_xyyzz_xxxyz[k] = -g_0_x_xyzz_xxxyz[k] * ab_y + g_0_x_xyzz_xxxyyz[k];

                g_0_x_xyyzz_xxxzz[k] = -g_0_x_xyzz_xxxzz[k] * ab_y + g_0_x_xyzz_xxxyzz[k];

                g_0_x_xyyzz_xxyyy[k] = -g_0_x_xyzz_xxyyy[k] * ab_y + g_0_x_xyzz_xxyyyy[k];

                g_0_x_xyyzz_xxyyz[k] = -g_0_x_xyzz_xxyyz[k] * ab_y + g_0_x_xyzz_xxyyyz[k];

                g_0_x_xyyzz_xxyzz[k] = -g_0_x_xyzz_xxyzz[k] * ab_y + g_0_x_xyzz_xxyyzz[k];

                g_0_x_xyyzz_xxzzz[k] = -g_0_x_xyzz_xxzzz[k] * ab_y + g_0_x_xyzz_xxyzzz[k];

                g_0_x_xyyzz_xyyyy[k] = -g_0_x_xyzz_xyyyy[k] * ab_y + g_0_x_xyzz_xyyyyy[k];

                g_0_x_xyyzz_xyyyz[k] = -g_0_x_xyzz_xyyyz[k] * ab_y + g_0_x_xyzz_xyyyyz[k];

                g_0_x_xyyzz_xyyzz[k] = -g_0_x_xyzz_xyyzz[k] * ab_y + g_0_x_xyzz_xyyyzz[k];

                g_0_x_xyyzz_xyzzz[k] = -g_0_x_xyzz_xyzzz[k] * ab_y + g_0_x_xyzz_xyyzzz[k];

                g_0_x_xyyzz_xzzzz[k] = -g_0_x_xyzz_xzzzz[k] * ab_y + g_0_x_xyzz_xyzzzz[k];

                g_0_x_xyyzz_yyyyy[k] = -g_0_x_xyzz_yyyyy[k] * ab_y + g_0_x_xyzz_yyyyyy[k];

                g_0_x_xyyzz_yyyyz[k] = -g_0_x_xyzz_yyyyz[k] * ab_y + g_0_x_xyzz_yyyyyz[k];

                g_0_x_xyyzz_yyyzz[k] = -g_0_x_xyzz_yyyzz[k] * ab_y + g_0_x_xyzz_yyyyzz[k];

                g_0_x_xyyzz_yyzzz[k] = -g_0_x_xyzz_yyzzz[k] * ab_y + g_0_x_xyzz_yyyzzz[k];

                g_0_x_xyyzz_yzzzz[k] = -g_0_x_xyzz_yzzzz[k] * ab_y + g_0_x_xyzz_yyzzzz[k];

                g_0_x_xyyzz_zzzzz[k] = -g_0_x_xyzz_zzzzz[k] * ab_y + g_0_x_xyzz_yzzzzz[k];
            }

            /// Set up 273-294 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzz_xxxxx = cbuffer.data(hh_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxy = cbuffer.data(hh_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxxz = cbuffer.data(hh_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxyy = cbuffer.data(hh_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxyz = cbuffer.data(hh_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxzz = cbuffer.data(hh_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxyyy = cbuffer.data(hh_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxyyz = cbuffer.data(hh_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxyzz = cbuffer.data(hh_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxzzz = cbuffer.data(hh_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyyyy = cbuffer.data(hh_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyyyz = cbuffer.data(hh_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyyzz = cbuffer.data(hh_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyzzz = cbuffer.data(hh_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_xyzzz_xzzzz = cbuffer.data(hh_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyyyy = cbuffer.data(hh_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyyyz = cbuffer.data(hh_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyyzz = cbuffer.data(hh_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyzzz = cbuffer.data(hh_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_xyzzz_yzzzz = cbuffer.data(hh_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_xyzzz_zzzzz = cbuffer.data(hh_geom_01_off + 293 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzz_xxxxx, g_0_x_xyzzz_xxxxy, g_0_x_xyzzz_xxxxz, g_0_x_xyzzz_xxxyy, g_0_x_xyzzz_xxxyz, g_0_x_xyzzz_xxxzz, g_0_x_xyzzz_xxyyy, g_0_x_xyzzz_xxyyz, g_0_x_xyzzz_xxyzz, g_0_x_xyzzz_xxzzz, g_0_x_xyzzz_xyyyy, g_0_x_xyzzz_xyyyz, g_0_x_xyzzz_xyyzz, g_0_x_xyzzz_xyzzz, g_0_x_xyzzz_xzzzz, g_0_x_xyzzz_yyyyy, g_0_x_xyzzz_yyyyz, g_0_x_xyzzz_yyyzz, g_0_x_xyzzz_yyzzz, g_0_x_xyzzz_yzzzz, g_0_x_xyzzz_zzzzz, g_0_x_xzzz_xxxxx, g_0_x_xzzz_xxxxxy, g_0_x_xzzz_xxxxy, g_0_x_xzzz_xxxxyy, g_0_x_xzzz_xxxxyz, g_0_x_xzzz_xxxxz, g_0_x_xzzz_xxxyy, g_0_x_xzzz_xxxyyy, g_0_x_xzzz_xxxyyz, g_0_x_xzzz_xxxyz, g_0_x_xzzz_xxxyzz, g_0_x_xzzz_xxxzz, g_0_x_xzzz_xxyyy, g_0_x_xzzz_xxyyyy, g_0_x_xzzz_xxyyyz, g_0_x_xzzz_xxyyz, g_0_x_xzzz_xxyyzz, g_0_x_xzzz_xxyzz, g_0_x_xzzz_xxyzzz, g_0_x_xzzz_xxzzz, g_0_x_xzzz_xyyyy, g_0_x_xzzz_xyyyyy, g_0_x_xzzz_xyyyyz, g_0_x_xzzz_xyyyz, g_0_x_xzzz_xyyyzz, g_0_x_xzzz_xyyzz, g_0_x_xzzz_xyyzzz, g_0_x_xzzz_xyzzz, g_0_x_xzzz_xyzzzz, g_0_x_xzzz_xzzzz, g_0_x_xzzz_yyyyy, g_0_x_xzzz_yyyyyy, g_0_x_xzzz_yyyyyz, g_0_x_xzzz_yyyyz, g_0_x_xzzz_yyyyzz, g_0_x_xzzz_yyyzz, g_0_x_xzzz_yyyzzz, g_0_x_xzzz_yyzzz, g_0_x_xzzz_yyzzzz, g_0_x_xzzz_yzzzz, g_0_x_xzzz_yzzzzz, g_0_x_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzz_xxxxx[k] = -g_0_x_xzzz_xxxxx[k] * ab_y + g_0_x_xzzz_xxxxxy[k];

                g_0_x_xyzzz_xxxxy[k] = -g_0_x_xzzz_xxxxy[k] * ab_y + g_0_x_xzzz_xxxxyy[k];

                g_0_x_xyzzz_xxxxz[k] = -g_0_x_xzzz_xxxxz[k] * ab_y + g_0_x_xzzz_xxxxyz[k];

                g_0_x_xyzzz_xxxyy[k] = -g_0_x_xzzz_xxxyy[k] * ab_y + g_0_x_xzzz_xxxyyy[k];

                g_0_x_xyzzz_xxxyz[k] = -g_0_x_xzzz_xxxyz[k] * ab_y + g_0_x_xzzz_xxxyyz[k];

                g_0_x_xyzzz_xxxzz[k] = -g_0_x_xzzz_xxxzz[k] * ab_y + g_0_x_xzzz_xxxyzz[k];

                g_0_x_xyzzz_xxyyy[k] = -g_0_x_xzzz_xxyyy[k] * ab_y + g_0_x_xzzz_xxyyyy[k];

                g_0_x_xyzzz_xxyyz[k] = -g_0_x_xzzz_xxyyz[k] * ab_y + g_0_x_xzzz_xxyyyz[k];

                g_0_x_xyzzz_xxyzz[k] = -g_0_x_xzzz_xxyzz[k] * ab_y + g_0_x_xzzz_xxyyzz[k];

                g_0_x_xyzzz_xxzzz[k] = -g_0_x_xzzz_xxzzz[k] * ab_y + g_0_x_xzzz_xxyzzz[k];

                g_0_x_xyzzz_xyyyy[k] = -g_0_x_xzzz_xyyyy[k] * ab_y + g_0_x_xzzz_xyyyyy[k];

                g_0_x_xyzzz_xyyyz[k] = -g_0_x_xzzz_xyyyz[k] * ab_y + g_0_x_xzzz_xyyyyz[k];

                g_0_x_xyzzz_xyyzz[k] = -g_0_x_xzzz_xyyzz[k] * ab_y + g_0_x_xzzz_xyyyzz[k];

                g_0_x_xyzzz_xyzzz[k] = -g_0_x_xzzz_xyzzz[k] * ab_y + g_0_x_xzzz_xyyzzz[k];

                g_0_x_xyzzz_xzzzz[k] = -g_0_x_xzzz_xzzzz[k] * ab_y + g_0_x_xzzz_xyzzzz[k];

                g_0_x_xyzzz_yyyyy[k] = -g_0_x_xzzz_yyyyy[k] * ab_y + g_0_x_xzzz_yyyyyy[k];

                g_0_x_xyzzz_yyyyz[k] = -g_0_x_xzzz_yyyyz[k] * ab_y + g_0_x_xzzz_yyyyyz[k];

                g_0_x_xyzzz_yyyzz[k] = -g_0_x_xzzz_yyyzz[k] * ab_y + g_0_x_xzzz_yyyyzz[k];

                g_0_x_xyzzz_yyzzz[k] = -g_0_x_xzzz_yyzzz[k] * ab_y + g_0_x_xzzz_yyyzzz[k];

                g_0_x_xyzzz_yzzzz[k] = -g_0_x_xzzz_yzzzz[k] * ab_y + g_0_x_xzzz_yyzzzz[k];

                g_0_x_xyzzz_zzzzz[k] = -g_0_x_xzzz_zzzzz[k] * ab_y + g_0_x_xzzz_yzzzzz[k];
            }

            /// Set up 294-315 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzz_xxxxx = cbuffer.data(hh_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxy = cbuffer.data(hh_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxxz = cbuffer.data(hh_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxyy = cbuffer.data(hh_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxyz = cbuffer.data(hh_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxzz = cbuffer.data(hh_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxyyy = cbuffer.data(hh_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxyyz = cbuffer.data(hh_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxyzz = cbuffer.data(hh_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxzzz = cbuffer.data(hh_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyyyy = cbuffer.data(hh_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyyyz = cbuffer.data(hh_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyyzz = cbuffer.data(hh_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyzzz = cbuffer.data(hh_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_xzzzz_xzzzz = cbuffer.data(hh_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyyyy = cbuffer.data(hh_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyyyz = cbuffer.data(hh_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyyzz = cbuffer.data(hh_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyzzz = cbuffer.data(hh_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_xzzzz_yzzzz = cbuffer.data(hh_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_xzzzz_zzzzz = cbuffer.data(hh_geom_01_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzz_xxxxx, g_0_x_xzzz_xxxxxz, g_0_x_xzzz_xxxxy, g_0_x_xzzz_xxxxyz, g_0_x_xzzz_xxxxz, g_0_x_xzzz_xxxxzz, g_0_x_xzzz_xxxyy, g_0_x_xzzz_xxxyyz, g_0_x_xzzz_xxxyz, g_0_x_xzzz_xxxyzz, g_0_x_xzzz_xxxzz, g_0_x_xzzz_xxxzzz, g_0_x_xzzz_xxyyy, g_0_x_xzzz_xxyyyz, g_0_x_xzzz_xxyyz, g_0_x_xzzz_xxyyzz, g_0_x_xzzz_xxyzz, g_0_x_xzzz_xxyzzz, g_0_x_xzzz_xxzzz, g_0_x_xzzz_xxzzzz, g_0_x_xzzz_xyyyy, g_0_x_xzzz_xyyyyz, g_0_x_xzzz_xyyyz, g_0_x_xzzz_xyyyzz, g_0_x_xzzz_xyyzz, g_0_x_xzzz_xyyzzz, g_0_x_xzzz_xyzzz, g_0_x_xzzz_xyzzzz, g_0_x_xzzz_xzzzz, g_0_x_xzzz_xzzzzz, g_0_x_xzzz_yyyyy, g_0_x_xzzz_yyyyyz, g_0_x_xzzz_yyyyz, g_0_x_xzzz_yyyyzz, g_0_x_xzzz_yyyzz, g_0_x_xzzz_yyyzzz, g_0_x_xzzz_yyzzz, g_0_x_xzzz_yyzzzz, g_0_x_xzzz_yzzzz, g_0_x_xzzz_yzzzzz, g_0_x_xzzz_zzzzz, g_0_x_xzzz_zzzzzz, g_0_x_xzzzz_xxxxx, g_0_x_xzzzz_xxxxy, g_0_x_xzzzz_xxxxz, g_0_x_xzzzz_xxxyy, g_0_x_xzzzz_xxxyz, g_0_x_xzzzz_xxxzz, g_0_x_xzzzz_xxyyy, g_0_x_xzzzz_xxyyz, g_0_x_xzzzz_xxyzz, g_0_x_xzzzz_xxzzz, g_0_x_xzzzz_xyyyy, g_0_x_xzzzz_xyyyz, g_0_x_xzzzz_xyyzz, g_0_x_xzzzz_xyzzz, g_0_x_xzzzz_xzzzz, g_0_x_xzzzz_yyyyy, g_0_x_xzzzz_yyyyz, g_0_x_xzzzz_yyyzz, g_0_x_xzzzz_yyzzz, g_0_x_xzzzz_yzzzz, g_0_x_xzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzz_xxxxx[k] = -g_0_x_xzzz_xxxxx[k] * ab_z + g_0_x_xzzz_xxxxxz[k];

                g_0_x_xzzzz_xxxxy[k] = -g_0_x_xzzz_xxxxy[k] * ab_z + g_0_x_xzzz_xxxxyz[k];

                g_0_x_xzzzz_xxxxz[k] = -g_0_x_xzzz_xxxxz[k] * ab_z + g_0_x_xzzz_xxxxzz[k];

                g_0_x_xzzzz_xxxyy[k] = -g_0_x_xzzz_xxxyy[k] * ab_z + g_0_x_xzzz_xxxyyz[k];

                g_0_x_xzzzz_xxxyz[k] = -g_0_x_xzzz_xxxyz[k] * ab_z + g_0_x_xzzz_xxxyzz[k];

                g_0_x_xzzzz_xxxzz[k] = -g_0_x_xzzz_xxxzz[k] * ab_z + g_0_x_xzzz_xxxzzz[k];

                g_0_x_xzzzz_xxyyy[k] = -g_0_x_xzzz_xxyyy[k] * ab_z + g_0_x_xzzz_xxyyyz[k];

                g_0_x_xzzzz_xxyyz[k] = -g_0_x_xzzz_xxyyz[k] * ab_z + g_0_x_xzzz_xxyyzz[k];

                g_0_x_xzzzz_xxyzz[k] = -g_0_x_xzzz_xxyzz[k] * ab_z + g_0_x_xzzz_xxyzzz[k];

                g_0_x_xzzzz_xxzzz[k] = -g_0_x_xzzz_xxzzz[k] * ab_z + g_0_x_xzzz_xxzzzz[k];

                g_0_x_xzzzz_xyyyy[k] = -g_0_x_xzzz_xyyyy[k] * ab_z + g_0_x_xzzz_xyyyyz[k];

                g_0_x_xzzzz_xyyyz[k] = -g_0_x_xzzz_xyyyz[k] * ab_z + g_0_x_xzzz_xyyyzz[k];

                g_0_x_xzzzz_xyyzz[k] = -g_0_x_xzzz_xyyzz[k] * ab_z + g_0_x_xzzz_xyyzzz[k];

                g_0_x_xzzzz_xyzzz[k] = -g_0_x_xzzz_xyzzz[k] * ab_z + g_0_x_xzzz_xyzzzz[k];

                g_0_x_xzzzz_xzzzz[k] = -g_0_x_xzzz_xzzzz[k] * ab_z + g_0_x_xzzz_xzzzzz[k];

                g_0_x_xzzzz_yyyyy[k] = -g_0_x_xzzz_yyyyy[k] * ab_z + g_0_x_xzzz_yyyyyz[k];

                g_0_x_xzzzz_yyyyz[k] = -g_0_x_xzzz_yyyyz[k] * ab_z + g_0_x_xzzz_yyyyzz[k];

                g_0_x_xzzzz_yyyzz[k] = -g_0_x_xzzz_yyyzz[k] * ab_z + g_0_x_xzzz_yyyzzz[k];

                g_0_x_xzzzz_yyzzz[k] = -g_0_x_xzzz_yyzzz[k] * ab_z + g_0_x_xzzz_yyzzzz[k];

                g_0_x_xzzzz_yzzzz[k] = -g_0_x_xzzz_yzzzz[k] * ab_z + g_0_x_xzzz_yzzzzz[k];

                g_0_x_xzzzz_zzzzz[k] = -g_0_x_xzzz_zzzzz[k] * ab_z + g_0_x_xzzz_zzzzzz[k];
            }

            /// Set up 315-336 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyy_xxxxx = cbuffer.data(hh_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxy = cbuffer.data(hh_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxxz = cbuffer.data(hh_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxyy = cbuffer.data(hh_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxyz = cbuffer.data(hh_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxzz = cbuffer.data(hh_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxyyy = cbuffer.data(hh_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxyyz = cbuffer.data(hh_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxyzz = cbuffer.data(hh_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxzzz = cbuffer.data(hh_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyyyy = cbuffer.data(hh_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyyyz = cbuffer.data(hh_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyyzz = cbuffer.data(hh_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyzzz = cbuffer.data(hh_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_yyyyy_xzzzz = cbuffer.data(hh_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyyyy = cbuffer.data(hh_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyyyz = cbuffer.data(hh_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyyzz = cbuffer.data(hh_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyzzz = cbuffer.data(hh_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_yyyyy_yzzzz = cbuffer.data(hh_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_yyyyy_zzzzz = cbuffer.data(hh_geom_01_off + 335 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyy_xxxxx, g_0_x_yyyy_xxxxxy, g_0_x_yyyy_xxxxy, g_0_x_yyyy_xxxxyy, g_0_x_yyyy_xxxxyz, g_0_x_yyyy_xxxxz, g_0_x_yyyy_xxxyy, g_0_x_yyyy_xxxyyy, g_0_x_yyyy_xxxyyz, g_0_x_yyyy_xxxyz, g_0_x_yyyy_xxxyzz, g_0_x_yyyy_xxxzz, g_0_x_yyyy_xxyyy, g_0_x_yyyy_xxyyyy, g_0_x_yyyy_xxyyyz, g_0_x_yyyy_xxyyz, g_0_x_yyyy_xxyyzz, g_0_x_yyyy_xxyzz, g_0_x_yyyy_xxyzzz, g_0_x_yyyy_xxzzz, g_0_x_yyyy_xyyyy, g_0_x_yyyy_xyyyyy, g_0_x_yyyy_xyyyyz, g_0_x_yyyy_xyyyz, g_0_x_yyyy_xyyyzz, g_0_x_yyyy_xyyzz, g_0_x_yyyy_xyyzzz, g_0_x_yyyy_xyzzz, g_0_x_yyyy_xyzzzz, g_0_x_yyyy_xzzzz, g_0_x_yyyy_yyyyy, g_0_x_yyyy_yyyyyy, g_0_x_yyyy_yyyyyz, g_0_x_yyyy_yyyyz, g_0_x_yyyy_yyyyzz, g_0_x_yyyy_yyyzz, g_0_x_yyyy_yyyzzz, g_0_x_yyyy_yyzzz, g_0_x_yyyy_yyzzzz, g_0_x_yyyy_yzzzz, g_0_x_yyyy_yzzzzz, g_0_x_yyyy_zzzzz, g_0_x_yyyyy_xxxxx, g_0_x_yyyyy_xxxxy, g_0_x_yyyyy_xxxxz, g_0_x_yyyyy_xxxyy, g_0_x_yyyyy_xxxyz, g_0_x_yyyyy_xxxzz, g_0_x_yyyyy_xxyyy, g_0_x_yyyyy_xxyyz, g_0_x_yyyyy_xxyzz, g_0_x_yyyyy_xxzzz, g_0_x_yyyyy_xyyyy, g_0_x_yyyyy_xyyyz, g_0_x_yyyyy_xyyzz, g_0_x_yyyyy_xyzzz, g_0_x_yyyyy_xzzzz, g_0_x_yyyyy_yyyyy, g_0_x_yyyyy_yyyyz, g_0_x_yyyyy_yyyzz, g_0_x_yyyyy_yyzzz, g_0_x_yyyyy_yzzzz, g_0_x_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyy_xxxxx[k] = -g_0_x_yyyy_xxxxx[k] * ab_y + g_0_x_yyyy_xxxxxy[k];

                g_0_x_yyyyy_xxxxy[k] = -g_0_x_yyyy_xxxxy[k] * ab_y + g_0_x_yyyy_xxxxyy[k];

                g_0_x_yyyyy_xxxxz[k] = -g_0_x_yyyy_xxxxz[k] * ab_y + g_0_x_yyyy_xxxxyz[k];

                g_0_x_yyyyy_xxxyy[k] = -g_0_x_yyyy_xxxyy[k] * ab_y + g_0_x_yyyy_xxxyyy[k];

                g_0_x_yyyyy_xxxyz[k] = -g_0_x_yyyy_xxxyz[k] * ab_y + g_0_x_yyyy_xxxyyz[k];

                g_0_x_yyyyy_xxxzz[k] = -g_0_x_yyyy_xxxzz[k] * ab_y + g_0_x_yyyy_xxxyzz[k];

                g_0_x_yyyyy_xxyyy[k] = -g_0_x_yyyy_xxyyy[k] * ab_y + g_0_x_yyyy_xxyyyy[k];

                g_0_x_yyyyy_xxyyz[k] = -g_0_x_yyyy_xxyyz[k] * ab_y + g_0_x_yyyy_xxyyyz[k];

                g_0_x_yyyyy_xxyzz[k] = -g_0_x_yyyy_xxyzz[k] * ab_y + g_0_x_yyyy_xxyyzz[k];

                g_0_x_yyyyy_xxzzz[k] = -g_0_x_yyyy_xxzzz[k] * ab_y + g_0_x_yyyy_xxyzzz[k];

                g_0_x_yyyyy_xyyyy[k] = -g_0_x_yyyy_xyyyy[k] * ab_y + g_0_x_yyyy_xyyyyy[k];

                g_0_x_yyyyy_xyyyz[k] = -g_0_x_yyyy_xyyyz[k] * ab_y + g_0_x_yyyy_xyyyyz[k];

                g_0_x_yyyyy_xyyzz[k] = -g_0_x_yyyy_xyyzz[k] * ab_y + g_0_x_yyyy_xyyyzz[k];

                g_0_x_yyyyy_xyzzz[k] = -g_0_x_yyyy_xyzzz[k] * ab_y + g_0_x_yyyy_xyyzzz[k];

                g_0_x_yyyyy_xzzzz[k] = -g_0_x_yyyy_xzzzz[k] * ab_y + g_0_x_yyyy_xyzzzz[k];

                g_0_x_yyyyy_yyyyy[k] = -g_0_x_yyyy_yyyyy[k] * ab_y + g_0_x_yyyy_yyyyyy[k];

                g_0_x_yyyyy_yyyyz[k] = -g_0_x_yyyy_yyyyz[k] * ab_y + g_0_x_yyyy_yyyyyz[k];

                g_0_x_yyyyy_yyyzz[k] = -g_0_x_yyyy_yyyzz[k] * ab_y + g_0_x_yyyy_yyyyzz[k];

                g_0_x_yyyyy_yyzzz[k] = -g_0_x_yyyy_yyzzz[k] * ab_y + g_0_x_yyyy_yyyzzz[k];

                g_0_x_yyyyy_yzzzz[k] = -g_0_x_yyyy_yzzzz[k] * ab_y + g_0_x_yyyy_yyzzzz[k];

                g_0_x_yyyyy_zzzzz[k] = -g_0_x_yyyy_zzzzz[k] * ab_y + g_0_x_yyyy_yzzzzz[k];
            }

            /// Set up 336-357 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyz_xxxxx = cbuffer.data(hh_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxy = cbuffer.data(hh_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxxz = cbuffer.data(hh_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxyy = cbuffer.data(hh_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxyz = cbuffer.data(hh_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxzz = cbuffer.data(hh_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxyyy = cbuffer.data(hh_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxyyz = cbuffer.data(hh_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxyzz = cbuffer.data(hh_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxzzz = cbuffer.data(hh_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyyyy = cbuffer.data(hh_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyyyz = cbuffer.data(hh_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyyzz = cbuffer.data(hh_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyzzz = cbuffer.data(hh_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_x_yyyyz_xzzzz = cbuffer.data(hh_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyyyy = cbuffer.data(hh_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyyyz = cbuffer.data(hh_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyyzz = cbuffer.data(hh_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyzzz = cbuffer.data(hh_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_yyyyz_yzzzz = cbuffer.data(hh_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_yyyyz_zzzzz = cbuffer.data(hh_geom_01_off + 356 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyz_xxxxx, g_0_x_yyyyz_xxxxy, g_0_x_yyyyz_xxxxz, g_0_x_yyyyz_xxxyy, g_0_x_yyyyz_xxxyz, g_0_x_yyyyz_xxxzz, g_0_x_yyyyz_xxyyy, g_0_x_yyyyz_xxyyz, g_0_x_yyyyz_xxyzz, g_0_x_yyyyz_xxzzz, g_0_x_yyyyz_xyyyy, g_0_x_yyyyz_xyyyz, g_0_x_yyyyz_xyyzz, g_0_x_yyyyz_xyzzz, g_0_x_yyyyz_xzzzz, g_0_x_yyyyz_yyyyy, g_0_x_yyyyz_yyyyz, g_0_x_yyyyz_yyyzz, g_0_x_yyyyz_yyzzz, g_0_x_yyyyz_yzzzz, g_0_x_yyyyz_zzzzz, g_0_x_yyyz_xxxxx, g_0_x_yyyz_xxxxxy, g_0_x_yyyz_xxxxy, g_0_x_yyyz_xxxxyy, g_0_x_yyyz_xxxxyz, g_0_x_yyyz_xxxxz, g_0_x_yyyz_xxxyy, g_0_x_yyyz_xxxyyy, g_0_x_yyyz_xxxyyz, g_0_x_yyyz_xxxyz, g_0_x_yyyz_xxxyzz, g_0_x_yyyz_xxxzz, g_0_x_yyyz_xxyyy, g_0_x_yyyz_xxyyyy, g_0_x_yyyz_xxyyyz, g_0_x_yyyz_xxyyz, g_0_x_yyyz_xxyyzz, g_0_x_yyyz_xxyzz, g_0_x_yyyz_xxyzzz, g_0_x_yyyz_xxzzz, g_0_x_yyyz_xyyyy, g_0_x_yyyz_xyyyyy, g_0_x_yyyz_xyyyyz, g_0_x_yyyz_xyyyz, g_0_x_yyyz_xyyyzz, g_0_x_yyyz_xyyzz, g_0_x_yyyz_xyyzzz, g_0_x_yyyz_xyzzz, g_0_x_yyyz_xyzzzz, g_0_x_yyyz_xzzzz, g_0_x_yyyz_yyyyy, g_0_x_yyyz_yyyyyy, g_0_x_yyyz_yyyyyz, g_0_x_yyyz_yyyyz, g_0_x_yyyz_yyyyzz, g_0_x_yyyz_yyyzz, g_0_x_yyyz_yyyzzz, g_0_x_yyyz_yyzzz, g_0_x_yyyz_yyzzzz, g_0_x_yyyz_yzzzz, g_0_x_yyyz_yzzzzz, g_0_x_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyz_xxxxx[k] = -g_0_x_yyyz_xxxxx[k] * ab_y + g_0_x_yyyz_xxxxxy[k];

                g_0_x_yyyyz_xxxxy[k] = -g_0_x_yyyz_xxxxy[k] * ab_y + g_0_x_yyyz_xxxxyy[k];

                g_0_x_yyyyz_xxxxz[k] = -g_0_x_yyyz_xxxxz[k] * ab_y + g_0_x_yyyz_xxxxyz[k];

                g_0_x_yyyyz_xxxyy[k] = -g_0_x_yyyz_xxxyy[k] * ab_y + g_0_x_yyyz_xxxyyy[k];

                g_0_x_yyyyz_xxxyz[k] = -g_0_x_yyyz_xxxyz[k] * ab_y + g_0_x_yyyz_xxxyyz[k];

                g_0_x_yyyyz_xxxzz[k] = -g_0_x_yyyz_xxxzz[k] * ab_y + g_0_x_yyyz_xxxyzz[k];

                g_0_x_yyyyz_xxyyy[k] = -g_0_x_yyyz_xxyyy[k] * ab_y + g_0_x_yyyz_xxyyyy[k];

                g_0_x_yyyyz_xxyyz[k] = -g_0_x_yyyz_xxyyz[k] * ab_y + g_0_x_yyyz_xxyyyz[k];

                g_0_x_yyyyz_xxyzz[k] = -g_0_x_yyyz_xxyzz[k] * ab_y + g_0_x_yyyz_xxyyzz[k];

                g_0_x_yyyyz_xxzzz[k] = -g_0_x_yyyz_xxzzz[k] * ab_y + g_0_x_yyyz_xxyzzz[k];

                g_0_x_yyyyz_xyyyy[k] = -g_0_x_yyyz_xyyyy[k] * ab_y + g_0_x_yyyz_xyyyyy[k];

                g_0_x_yyyyz_xyyyz[k] = -g_0_x_yyyz_xyyyz[k] * ab_y + g_0_x_yyyz_xyyyyz[k];

                g_0_x_yyyyz_xyyzz[k] = -g_0_x_yyyz_xyyzz[k] * ab_y + g_0_x_yyyz_xyyyzz[k];

                g_0_x_yyyyz_xyzzz[k] = -g_0_x_yyyz_xyzzz[k] * ab_y + g_0_x_yyyz_xyyzzz[k];

                g_0_x_yyyyz_xzzzz[k] = -g_0_x_yyyz_xzzzz[k] * ab_y + g_0_x_yyyz_xyzzzz[k];

                g_0_x_yyyyz_yyyyy[k] = -g_0_x_yyyz_yyyyy[k] * ab_y + g_0_x_yyyz_yyyyyy[k];

                g_0_x_yyyyz_yyyyz[k] = -g_0_x_yyyz_yyyyz[k] * ab_y + g_0_x_yyyz_yyyyyz[k];

                g_0_x_yyyyz_yyyzz[k] = -g_0_x_yyyz_yyyzz[k] * ab_y + g_0_x_yyyz_yyyyzz[k];

                g_0_x_yyyyz_yyzzz[k] = -g_0_x_yyyz_yyzzz[k] * ab_y + g_0_x_yyyz_yyyzzz[k];

                g_0_x_yyyyz_yzzzz[k] = -g_0_x_yyyz_yzzzz[k] * ab_y + g_0_x_yyyz_yyzzzz[k];

                g_0_x_yyyyz_zzzzz[k] = -g_0_x_yyyz_zzzzz[k] * ab_y + g_0_x_yyyz_yzzzzz[k];
            }

            /// Set up 357-378 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzz_xxxxx = cbuffer.data(hh_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxy = cbuffer.data(hh_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxxz = cbuffer.data(hh_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxyy = cbuffer.data(hh_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxyz = cbuffer.data(hh_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxzz = cbuffer.data(hh_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxyyy = cbuffer.data(hh_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxyyz = cbuffer.data(hh_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxyzz = cbuffer.data(hh_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxzzz = cbuffer.data(hh_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyyyy = cbuffer.data(hh_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyyyz = cbuffer.data(hh_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyyzz = cbuffer.data(hh_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyzzz = cbuffer.data(hh_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_x_yyyzz_xzzzz = cbuffer.data(hh_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyyyy = cbuffer.data(hh_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyyyz = cbuffer.data(hh_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyyzz = cbuffer.data(hh_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyzzz = cbuffer.data(hh_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_x_yyyzz_yzzzz = cbuffer.data(hh_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_x_yyyzz_zzzzz = cbuffer.data(hh_geom_01_off + 377 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzz_xxxxx, g_0_x_yyyzz_xxxxy, g_0_x_yyyzz_xxxxz, g_0_x_yyyzz_xxxyy, g_0_x_yyyzz_xxxyz, g_0_x_yyyzz_xxxzz, g_0_x_yyyzz_xxyyy, g_0_x_yyyzz_xxyyz, g_0_x_yyyzz_xxyzz, g_0_x_yyyzz_xxzzz, g_0_x_yyyzz_xyyyy, g_0_x_yyyzz_xyyyz, g_0_x_yyyzz_xyyzz, g_0_x_yyyzz_xyzzz, g_0_x_yyyzz_xzzzz, g_0_x_yyyzz_yyyyy, g_0_x_yyyzz_yyyyz, g_0_x_yyyzz_yyyzz, g_0_x_yyyzz_yyzzz, g_0_x_yyyzz_yzzzz, g_0_x_yyyzz_zzzzz, g_0_x_yyzz_xxxxx, g_0_x_yyzz_xxxxxy, g_0_x_yyzz_xxxxy, g_0_x_yyzz_xxxxyy, g_0_x_yyzz_xxxxyz, g_0_x_yyzz_xxxxz, g_0_x_yyzz_xxxyy, g_0_x_yyzz_xxxyyy, g_0_x_yyzz_xxxyyz, g_0_x_yyzz_xxxyz, g_0_x_yyzz_xxxyzz, g_0_x_yyzz_xxxzz, g_0_x_yyzz_xxyyy, g_0_x_yyzz_xxyyyy, g_0_x_yyzz_xxyyyz, g_0_x_yyzz_xxyyz, g_0_x_yyzz_xxyyzz, g_0_x_yyzz_xxyzz, g_0_x_yyzz_xxyzzz, g_0_x_yyzz_xxzzz, g_0_x_yyzz_xyyyy, g_0_x_yyzz_xyyyyy, g_0_x_yyzz_xyyyyz, g_0_x_yyzz_xyyyz, g_0_x_yyzz_xyyyzz, g_0_x_yyzz_xyyzz, g_0_x_yyzz_xyyzzz, g_0_x_yyzz_xyzzz, g_0_x_yyzz_xyzzzz, g_0_x_yyzz_xzzzz, g_0_x_yyzz_yyyyy, g_0_x_yyzz_yyyyyy, g_0_x_yyzz_yyyyyz, g_0_x_yyzz_yyyyz, g_0_x_yyzz_yyyyzz, g_0_x_yyzz_yyyzz, g_0_x_yyzz_yyyzzz, g_0_x_yyzz_yyzzz, g_0_x_yyzz_yyzzzz, g_0_x_yyzz_yzzzz, g_0_x_yyzz_yzzzzz, g_0_x_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzz_xxxxx[k] = -g_0_x_yyzz_xxxxx[k] * ab_y + g_0_x_yyzz_xxxxxy[k];

                g_0_x_yyyzz_xxxxy[k] = -g_0_x_yyzz_xxxxy[k] * ab_y + g_0_x_yyzz_xxxxyy[k];

                g_0_x_yyyzz_xxxxz[k] = -g_0_x_yyzz_xxxxz[k] * ab_y + g_0_x_yyzz_xxxxyz[k];

                g_0_x_yyyzz_xxxyy[k] = -g_0_x_yyzz_xxxyy[k] * ab_y + g_0_x_yyzz_xxxyyy[k];

                g_0_x_yyyzz_xxxyz[k] = -g_0_x_yyzz_xxxyz[k] * ab_y + g_0_x_yyzz_xxxyyz[k];

                g_0_x_yyyzz_xxxzz[k] = -g_0_x_yyzz_xxxzz[k] * ab_y + g_0_x_yyzz_xxxyzz[k];

                g_0_x_yyyzz_xxyyy[k] = -g_0_x_yyzz_xxyyy[k] * ab_y + g_0_x_yyzz_xxyyyy[k];

                g_0_x_yyyzz_xxyyz[k] = -g_0_x_yyzz_xxyyz[k] * ab_y + g_0_x_yyzz_xxyyyz[k];

                g_0_x_yyyzz_xxyzz[k] = -g_0_x_yyzz_xxyzz[k] * ab_y + g_0_x_yyzz_xxyyzz[k];

                g_0_x_yyyzz_xxzzz[k] = -g_0_x_yyzz_xxzzz[k] * ab_y + g_0_x_yyzz_xxyzzz[k];

                g_0_x_yyyzz_xyyyy[k] = -g_0_x_yyzz_xyyyy[k] * ab_y + g_0_x_yyzz_xyyyyy[k];

                g_0_x_yyyzz_xyyyz[k] = -g_0_x_yyzz_xyyyz[k] * ab_y + g_0_x_yyzz_xyyyyz[k];

                g_0_x_yyyzz_xyyzz[k] = -g_0_x_yyzz_xyyzz[k] * ab_y + g_0_x_yyzz_xyyyzz[k];

                g_0_x_yyyzz_xyzzz[k] = -g_0_x_yyzz_xyzzz[k] * ab_y + g_0_x_yyzz_xyyzzz[k];

                g_0_x_yyyzz_xzzzz[k] = -g_0_x_yyzz_xzzzz[k] * ab_y + g_0_x_yyzz_xyzzzz[k];

                g_0_x_yyyzz_yyyyy[k] = -g_0_x_yyzz_yyyyy[k] * ab_y + g_0_x_yyzz_yyyyyy[k];

                g_0_x_yyyzz_yyyyz[k] = -g_0_x_yyzz_yyyyz[k] * ab_y + g_0_x_yyzz_yyyyyz[k];

                g_0_x_yyyzz_yyyzz[k] = -g_0_x_yyzz_yyyzz[k] * ab_y + g_0_x_yyzz_yyyyzz[k];

                g_0_x_yyyzz_yyzzz[k] = -g_0_x_yyzz_yyzzz[k] * ab_y + g_0_x_yyzz_yyyzzz[k];

                g_0_x_yyyzz_yzzzz[k] = -g_0_x_yyzz_yzzzz[k] * ab_y + g_0_x_yyzz_yyzzzz[k];

                g_0_x_yyyzz_zzzzz[k] = -g_0_x_yyzz_zzzzz[k] * ab_y + g_0_x_yyzz_yzzzzz[k];
            }

            /// Set up 378-399 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzz_xxxxx = cbuffer.data(hh_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxy = cbuffer.data(hh_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxxz = cbuffer.data(hh_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxyy = cbuffer.data(hh_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxyz = cbuffer.data(hh_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxzz = cbuffer.data(hh_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxyyy = cbuffer.data(hh_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxyyz = cbuffer.data(hh_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxyzz = cbuffer.data(hh_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxzzz = cbuffer.data(hh_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyyyy = cbuffer.data(hh_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyyyz = cbuffer.data(hh_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyyzz = cbuffer.data(hh_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyzzz = cbuffer.data(hh_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_x_yyzzz_xzzzz = cbuffer.data(hh_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyyyy = cbuffer.data(hh_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyyyz = cbuffer.data(hh_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyyzz = cbuffer.data(hh_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyzzz = cbuffer.data(hh_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_x_yyzzz_yzzzz = cbuffer.data(hh_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_x_yyzzz_zzzzz = cbuffer.data(hh_geom_01_off + 398 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzz_xxxxx, g_0_x_yyzzz_xxxxy, g_0_x_yyzzz_xxxxz, g_0_x_yyzzz_xxxyy, g_0_x_yyzzz_xxxyz, g_0_x_yyzzz_xxxzz, g_0_x_yyzzz_xxyyy, g_0_x_yyzzz_xxyyz, g_0_x_yyzzz_xxyzz, g_0_x_yyzzz_xxzzz, g_0_x_yyzzz_xyyyy, g_0_x_yyzzz_xyyyz, g_0_x_yyzzz_xyyzz, g_0_x_yyzzz_xyzzz, g_0_x_yyzzz_xzzzz, g_0_x_yyzzz_yyyyy, g_0_x_yyzzz_yyyyz, g_0_x_yyzzz_yyyzz, g_0_x_yyzzz_yyzzz, g_0_x_yyzzz_yzzzz, g_0_x_yyzzz_zzzzz, g_0_x_yzzz_xxxxx, g_0_x_yzzz_xxxxxy, g_0_x_yzzz_xxxxy, g_0_x_yzzz_xxxxyy, g_0_x_yzzz_xxxxyz, g_0_x_yzzz_xxxxz, g_0_x_yzzz_xxxyy, g_0_x_yzzz_xxxyyy, g_0_x_yzzz_xxxyyz, g_0_x_yzzz_xxxyz, g_0_x_yzzz_xxxyzz, g_0_x_yzzz_xxxzz, g_0_x_yzzz_xxyyy, g_0_x_yzzz_xxyyyy, g_0_x_yzzz_xxyyyz, g_0_x_yzzz_xxyyz, g_0_x_yzzz_xxyyzz, g_0_x_yzzz_xxyzz, g_0_x_yzzz_xxyzzz, g_0_x_yzzz_xxzzz, g_0_x_yzzz_xyyyy, g_0_x_yzzz_xyyyyy, g_0_x_yzzz_xyyyyz, g_0_x_yzzz_xyyyz, g_0_x_yzzz_xyyyzz, g_0_x_yzzz_xyyzz, g_0_x_yzzz_xyyzzz, g_0_x_yzzz_xyzzz, g_0_x_yzzz_xyzzzz, g_0_x_yzzz_xzzzz, g_0_x_yzzz_yyyyy, g_0_x_yzzz_yyyyyy, g_0_x_yzzz_yyyyyz, g_0_x_yzzz_yyyyz, g_0_x_yzzz_yyyyzz, g_0_x_yzzz_yyyzz, g_0_x_yzzz_yyyzzz, g_0_x_yzzz_yyzzz, g_0_x_yzzz_yyzzzz, g_0_x_yzzz_yzzzz, g_0_x_yzzz_yzzzzz, g_0_x_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzz_xxxxx[k] = -g_0_x_yzzz_xxxxx[k] * ab_y + g_0_x_yzzz_xxxxxy[k];

                g_0_x_yyzzz_xxxxy[k] = -g_0_x_yzzz_xxxxy[k] * ab_y + g_0_x_yzzz_xxxxyy[k];

                g_0_x_yyzzz_xxxxz[k] = -g_0_x_yzzz_xxxxz[k] * ab_y + g_0_x_yzzz_xxxxyz[k];

                g_0_x_yyzzz_xxxyy[k] = -g_0_x_yzzz_xxxyy[k] * ab_y + g_0_x_yzzz_xxxyyy[k];

                g_0_x_yyzzz_xxxyz[k] = -g_0_x_yzzz_xxxyz[k] * ab_y + g_0_x_yzzz_xxxyyz[k];

                g_0_x_yyzzz_xxxzz[k] = -g_0_x_yzzz_xxxzz[k] * ab_y + g_0_x_yzzz_xxxyzz[k];

                g_0_x_yyzzz_xxyyy[k] = -g_0_x_yzzz_xxyyy[k] * ab_y + g_0_x_yzzz_xxyyyy[k];

                g_0_x_yyzzz_xxyyz[k] = -g_0_x_yzzz_xxyyz[k] * ab_y + g_0_x_yzzz_xxyyyz[k];

                g_0_x_yyzzz_xxyzz[k] = -g_0_x_yzzz_xxyzz[k] * ab_y + g_0_x_yzzz_xxyyzz[k];

                g_0_x_yyzzz_xxzzz[k] = -g_0_x_yzzz_xxzzz[k] * ab_y + g_0_x_yzzz_xxyzzz[k];

                g_0_x_yyzzz_xyyyy[k] = -g_0_x_yzzz_xyyyy[k] * ab_y + g_0_x_yzzz_xyyyyy[k];

                g_0_x_yyzzz_xyyyz[k] = -g_0_x_yzzz_xyyyz[k] * ab_y + g_0_x_yzzz_xyyyyz[k];

                g_0_x_yyzzz_xyyzz[k] = -g_0_x_yzzz_xyyzz[k] * ab_y + g_0_x_yzzz_xyyyzz[k];

                g_0_x_yyzzz_xyzzz[k] = -g_0_x_yzzz_xyzzz[k] * ab_y + g_0_x_yzzz_xyyzzz[k];

                g_0_x_yyzzz_xzzzz[k] = -g_0_x_yzzz_xzzzz[k] * ab_y + g_0_x_yzzz_xyzzzz[k];

                g_0_x_yyzzz_yyyyy[k] = -g_0_x_yzzz_yyyyy[k] * ab_y + g_0_x_yzzz_yyyyyy[k];

                g_0_x_yyzzz_yyyyz[k] = -g_0_x_yzzz_yyyyz[k] * ab_y + g_0_x_yzzz_yyyyyz[k];

                g_0_x_yyzzz_yyyzz[k] = -g_0_x_yzzz_yyyzz[k] * ab_y + g_0_x_yzzz_yyyyzz[k];

                g_0_x_yyzzz_yyzzz[k] = -g_0_x_yzzz_yyzzz[k] * ab_y + g_0_x_yzzz_yyyzzz[k];

                g_0_x_yyzzz_yzzzz[k] = -g_0_x_yzzz_yzzzz[k] * ab_y + g_0_x_yzzz_yyzzzz[k];

                g_0_x_yyzzz_zzzzz[k] = -g_0_x_yzzz_zzzzz[k] * ab_y + g_0_x_yzzz_yzzzzz[k];
            }

            /// Set up 399-420 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzz_xxxxx = cbuffer.data(hh_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxy = cbuffer.data(hh_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxxz = cbuffer.data(hh_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxyy = cbuffer.data(hh_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxyz = cbuffer.data(hh_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxzz = cbuffer.data(hh_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxyyy = cbuffer.data(hh_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxyyz = cbuffer.data(hh_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxyzz = cbuffer.data(hh_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxzzz = cbuffer.data(hh_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyyyy = cbuffer.data(hh_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyyyz = cbuffer.data(hh_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyyzz = cbuffer.data(hh_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyzzz = cbuffer.data(hh_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_x_yzzzz_xzzzz = cbuffer.data(hh_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyyyy = cbuffer.data(hh_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyyyz = cbuffer.data(hh_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyyzz = cbuffer.data(hh_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyzzz = cbuffer.data(hh_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_x_yzzzz_yzzzz = cbuffer.data(hh_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_x_yzzzz_zzzzz = cbuffer.data(hh_geom_01_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzz_xxxxx, g_0_x_yzzzz_xxxxy, g_0_x_yzzzz_xxxxz, g_0_x_yzzzz_xxxyy, g_0_x_yzzzz_xxxyz, g_0_x_yzzzz_xxxzz, g_0_x_yzzzz_xxyyy, g_0_x_yzzzz_xxyyz, g_0_x_yzzzz_xxyzz, g_0_x_yzzzz_xxzzz, g_0_x_yzzzz_xyyyy, g_0_x_yzzzz_xyyyz, g_0_x_yzzzz_xyyzz, g_0_x_yzzzz_xyzzz, g_0_x_yzzzz_xzzzz, g_0_x_yzzzz_yyyyy, g_0_x_yzzzz_yyyyz, g_0_x_yzzzz_yyyzz, g_0_x_yzzzz_yyzzz, g_0_x_yzzzz_yzzzz, g_0_x_yzzzz_zzzzz, g_0_x_zzzz_xxxxx, g_0_x_zzzz_xxxxxy, g_0_x_zzzz_xxxxy, g_0_x_zzzz_xxxxyy, g_0_x_zzzz_xxxxyz, g_0_x_zzzz_xxxxz, g_0_x_zzzz_xxxyy, g_0_x_zzzz_xxxyyy, g_0_x_zzzz_xxxyyz, g_0_x_zzzz_xxxyz, g_0_x_zzzz_xxxyzz, g_0_x_zzzz_xxxzz, g_0_x_zzzz_xxyyy, g_0_x_zzzz_xxyyyy, g_0_x_zzzz_xxyyyz, g_0_x_zzzz_xxyyz, g_0_x_zzzz_xxyyzz, g_0_x_zzzz_xxyzz, g_0_x_zzzz_xxyzzz, g_0_x_zzzz_xxzzz, g_0_x_zzzz_xyyyy, g_0_x_zzzz_xyyyyy, g_0_x_zzzz_xyyyyz, g_0_x_zzzz_xyyyz, g_0_x_zzzz_xyyyzz, g_0_x_zzzz_xyyzz, g_0_x_zzzz_xyyzzz, g_0_x_zzzz_xyzzz, g_0_x_zzzz_xyzzzz, g_0_x_zzzz_xzzzz, g_0_x_zzzz_yyyyy, g_0_x_zzzz_yyyyyy, g_0_x_zzzz_yyyyyz, g_0_x_zzzz_yyyyz, g_0_x_zzzz_yyyyzz, g_0_x_zzzz_yyyzz, g_0_x_zzzz_yyyzzz, g_0_x_zzzz_yyzzz, g_0_x_zzzz_yyzzzz, g_0_x_zzzz_yzzzz, g_0_x_zzzz_yzzzzz, g_0_x_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzz_xxxxx[k] = -g_0_x_zzzz_xxxxx[k] * ab_y + g_0_x_zzzz_xxxxxy[k];

                g_0_x_yzzzz_xxxxy[k] = -g_0_x_zzzz_xxxxy[k] * ab_y + g_0_x_zzzz_xxxxyy[k];

                g_0_x_yzzzz_xxxxz[k] = -g_0_x_zzzz_xxxxz[k] * ab_y + g_0_x_zzzz_xxxxyz[k];

                g_0_x_yzzzz_xxxyy[k] = -g_0_x_zzzz_xxxyy[k] * ab_y + g_0_x_zzzz_xxxyyy[k];

                g_0_x_yzzzz_xxxyz[k] = -g_0_x_zzzz_xxxyz[k] * ab_y + g_0_x_zzzz_xxxyyz[k];

                g_0_x_yzzzz_xxxzz[k] = -g_0_x_zzzz_xxxzz[k] * ab_y + g_0_x_zzzz_xxxyzz[k];

                g_0_x_yzzzz_xxyyy[k] = -g_0_x_zzzz_xxyyy[k] * ab_y + g_0_x_zzzz_xxyyyy[k];

                g_0_x_yzzzz_xxyyz[k] = -g_0_x_zzzz_xxyyz[k] * ab_y + g_0_x_zzzz_xxyyyz[k];

                g_0_x_yzzzz_xxyzz[k] = -g_0_x_zzzz_xxyzz[k] * ab_y + g_0_x_zzzz_xxyyzz[k];

                g_0_x_yzzzz_xxzzz[k] = -g_0_x_zzzz_xxzzz[k] * ab_y + g_0_x_zzzz_xxyzzz[k];

                g_0_x_yzzzz_xyyyy[k] = -g_0_x_zzzz_xyyyy[k] * ab_y + g_0_x_zzzz_xyyyyy[k];

                g_0_x_yzzzz_xyyyz[k] = -g_0_x_zzzz_xyyyz[k] * ab_y + g_0_x_zzzz_xyyyyz[k];

                g_0_x_yzzzz_xyyzz[k] = -g_0_x_zzzz_xyyzz[k] * ab_y + g_0_x_zzzz_xyyyzz[k];

                g_0_x_yzzzz_xyzzz[k] = -g_0_x_zzzz_xyzzz[k] * ab_y + g_0_x_zzzz_xyyzzz[k];

                g_0_x_yzzzz_xzzzz[k] = -g_0_x_zzzz_xzzzz[k] * ab_y + g_0_x_zzzz_xyzzzz[k];

                g_0_x_yzzzz_yyyyy[k] = -g_0_x_zzzz_yyyyy[k] * ab_y + g_0_x_zzzz_yyyyyy[k];

                g_0_x_yzzzz_yyyyz[k] = -g_0_x_zzzz_yyyyz[k] * ab_y + g_0_x_zzzz_yyyyyz[k];

                g_0_x_yzzzz_yyyzz[k] = -g_0_x_zzzz_yyyzz[k] * ab_y + g_0_x_zzzz_yyyyzz[k];

                g_0_x_yzzzz_yyzzz[k] = -g_0_x_zzzz_yyzzz[k] * ab_y + g_0_x_zzzz_yyyzzz[k];

                g_0_x_yzzzz_yzzzz[k] = -g_0_x_zzzz_yzzzz[k] * ab_y + g_0_x_zzzz_yyzzzz[k];

                g_0_x_yzzzz_zzzzz[k] = -g_0_x_zzzz_zzzzz[k] * ab_y + g_0_x_zzzz_yzzzzz[k];
            }

            /// Set up 420-441 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzz_xxxxx = cbuffer.data(hh_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxy = cbuffer.data(hh_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxxz = cbuffer.data(hh_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxyy = cbuffer.data(hh_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxyz = cbuffer.data(hh_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxzz = cbuffer.data(hh_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxyyy = cbuffer.data(hh_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxyyz = cbuffer.data(hh_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxyzz = cbuffer.data(hh_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxzzz = cbuffer.data(hh_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyyyy = cbuffer.data(hh_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyyyz = cbuffer.data(hh_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyyzz = cbuffer.data(hh_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyzzz = cbuffer.data(hh_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_x_zzzzz_xzzzz = cbuffer.data(hh_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyyyy = cbuffer.data(hh_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyyyz = cbuffer.data(hh_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyyzz = cbuffer.data(hh_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyzzz = cbuffer.data(hh_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_x_zzzzz_yzzzz = cbuffer.data(hh_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_x_zzzzz_zzzzz = cbuffer.data(hh_geom_01_off + 440 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzz_xxxxx, g_0_x_zzzz_xxxxxz, g_0_x_zzzz_xxxxy, g_0_x_zzzz_xxxxyz, g_0_x_zzzz_xxxxz, g_0_x_zzzz_xxxxzz, g_0_x_zzzz_xxxyy, g_0_x_zzzz_xxxyyz, g_0_x_zzzz_xxxyz, g_0_x_zzzz_xxxyzz, g_0_x_zzzz_xxxzz, g_0_x_zzzz_xxxzzz, g_0_x_zzzz_xxyyy, g_0_x_zzzz_xxyyyz, g_0_x_zzzz_xxyyz, g_0_x_zzzz_xxyyzz, g_0_x_zzzz_xxyzz, g_0_x_zzzz_xxyzzz, g_0_x_zzzz_xxzzz, g_0_x_zzzz_xxzzzz, g_0_x_zzzz_xyyyy, g_0_x_zzzz_xyyyyz, g_0_x_zzzz_xyyyz, g_0_x_zzzz_xyyyzz, g_0_x_zzzz_xyyzz, g_0_x_zzzz_xyyzzz, g_0_x_zzzz_xyzzz, g_0_x_zzzz_xyzzzz, g_0_x_zzzz_xzzzz, g_0_x_zzzz_xzzzzz, g_0_x_zzzz_yyyyy, g_0_x_zzzz_yyyyyz, g_0_x_zzzz_yyyyz, g_0_x_zzzz_yyyyzz, g_0_x_zzzz_yyyzz, g_0_x_zzzz_yyyzzz, g_0_x_zzzz_yyzzz, g_0_x_zzzz_yyzzzz, g_0_x_zzzz_yzzzz, g_0_x_zzzz_yzzzzz, g_0_x_zzzz_zzzzz, g_0_x_zzzz_zzzzzz, g_0_x_zzzzz_xxxxx, g_0_x_zzzzz_xxxxy, g_0_x_zzzzz_xxxxz, g_0_x_zzzzz_xxxyy, g_0_x_zzzzz_xxxyz, g_0_x_zzzzz_xxxzz, g_0_x_zzzzz_xxyyy, g_0_x_zzzzz_xxyyz, g_0_x_zzzzz_xxyzz, g_0_x_zzzzz_xxzzz, g_0_x_zzzzz_xyyyy, g_0_x_zzzzz_xyyyz, g_0_x_zzzzz_xyyzz, g_0_x_zzzzz_xyzzz, g_0_x_zzzzz_xzzzz, g_0_x_zzzzz_yyyyy, g_0_x_zzzzz_yyyyz, g_0_x_zzzzz_yyyzz, g_0_x_zzzzz_yyzzz, g_0_x_zzzzz_yzzzz, g_0_x_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzz_xxxxx[k] = -g_0_x_zzzz_xxxxx[k] * ab_z + g_0_x_zzzz_xxxxxz[k];

                g_0_x_zzzzz_xxxxy[k] = -g_0_x_zzzz_xxxxy[k] * ab_z + g_0_x_zzzz_xxxxyz[k];

                g_0_x_zzzzz_xxxxz[k] = -g_0_x_zzzz_xxxxz[k] * ab_z + g_0_x_zzzz_xxxxzz[k];

                g_0_x_zzzzz_xxxyy[k] = -g_0_x_zzzz_xxxyy[k] * ab_z + g_0_x_zzzz_xxxyyz[k];

                g_0_x_zzzzz_xxxyz[k] = -g_0_x_zzzz_xxxyz[k] * ab_z + g_0_x_zzzz_xxxyzz[k];

                g_0_x_zzzzz_xxxzz[k] = -g_0_x_zzzz_xxxzz[k] * ab_z + g_0_x_zzzz_xxxzzz[k];

                g_0_x_zzzzz_xxyyy[k] = -g_0_x_zzzz_xxyyy[k] * ab_z + g_0_x_zzzz_xxyyyz[k];

                g_0_x_zzzzz_xxyyz[k] = -g_0_x_zzzz_xxyyz[k] * ab_z + g_0_x_zzzz_xxyyzz[k];

                g_0_x_zzzzz_xxyzz[k] = -g_0_x_zzzz_xxyzz[k] * ab_z + g_0_x_zzzz_xxyzzz[k];

                g_0_x_zzzzz_xxzzz[k] = -g_0_x_zzzz_xxzzz[k] * ab_z + g_0_x_zzzz_xxzzzz[k];

                g_0_x_zzzzz_xyyyy[k] = -g_0_x_zzzz_xyyyy[k] * ab_z + g_0_x_zzzz_xyyyyz[k];

                g_0_x_zzzzz_xyyyz[k] = -g_0_x_zzzz_xyyyz[k] * ab_z + g_0_x_zzzz_xyyyzz[k];

                g_0_x_zzzzz_xyyzz[k] = -g_0_x_zzzz_xyyzz[k] * ab_z + g_0_x_zzzz_xyyzzz[k];

                g_0_x_zzzzz_xyzzz[k] = -g_0_x_zzzz_xyzzz[k] * ab_z + g_0_x_zzzz_xyzzzz[k];

                g_0_x_zzzzz_xzzzz[k] = -g_0_x_zzzz_xzzzz[k] * ab_z + g_0_x_zzzz_xzzzzz[k];

                g_0_x_zzzzz_yyyyy[k] = -g_0_x_zzzz_yyyyy[k] * ab_z + g_0_x_zzzz_yyyyyz[k];

                g_0_x_zzzzz_yyyyz[k] = -g_0_x_zzzz_yyyyz[k] * ab_z + g_0_x_zzzz_yyyyzz[k];

                g_0_x_zzzzz_yyyzz[k] = -g_0_x_zzzz_yyyzz[k] * ab_z + g_0_x_zzzz_yyyzzz[k];

                g_0_x_zzzzz_yyzzz[k] = -g_0_x_zzzz_yyzzz[k] * ab_z + g_0_x_zzzz_yyzzzz[k];

                g_0_x_zzzzz_yzzzz[k] = -g_0_x_zzzz_yzzzz[k] * ab_z + g_0_x_zzzz_yzzzzz[k];

                g_0_x_zzzzz_zzzzz[k] = -g_0_x_zzzz_zzzzz[k] * ab_z + g_0_x_zzzz_zzzzzz[k];
            }

            /// Set up 441-462 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxx_xxxxx = cbuffer.data(hh_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxy = cbuffer.data(hh_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxxz = cbuffer.data(hh_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxyy = cbuffer.data(hh_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxyz = cbuffer.data(hh_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxzz = cbuffer.data(hh_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxyyy = cbuffer.data(hh_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxyyz = cbuffer.data(hh_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxyzz = cbuffer.data(hh_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxzzz = cbuffer.data(hh_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyyyy = cbuffer.data(hh_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyyyz = cbuffer.data(hh_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyyzz = cbuffer.data(hh_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyzzz = cbuffer.data(hh_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_xxxxx_xzzzz = cbuffer.data(hh_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyyyy = cbuffer.data(hh_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyyyz = cbuffer.data(hh_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyyzz = cbuffer.data(hh_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyzzz = cbuffer.data(hh_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_xxxxx_yzzzz = cbuffer.data(hh_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_xxxxx_zzzzz = cbuffer.data(hh_geom_01_off + 461 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxx_xxxxx, g_0_y_xxxx_xxxxxx, g_0_y_xxxx_xxxxxy, g_0_y_xxxx_xxxxxz, g_0_y_xxxx_xxxxy, g_0_y_xxxx_xxxxyy, g_0_y_xxxx_xxxxyz, g_0_y_xxxx_xxxxz, g_0_y_xxxx_xxxxzz, g_0_y_xxxx_xxxyy, g_0_y_xxxx_xxxyyy, g_0_y_xxxx_xxxyyz, g_0_y_xxxx_xxxyz, g_0_y_xxxx_xxxyzz, g_0_y_xxxx_xxxzz, g_0_y_xxxx_xxxzzz, g_0_y_xxxx_xxyyy, g_0_y_xxxx_xxyyyy, g_0_y_xxxx_xxyyyz, g_0_y_xxxx_xxyyz, g_0_y_xxxx_xxyyzz, g_0_y_xxxx_xxyzz, g_0_y_xxxx_xxyzzz, g_0_y_xxxx_xxzzz, g_0_y_xxxx_xxzzzz, g_0_y_xxxx_xyyyy, g_0_y_xxxx_xyyyyy, g_0_y_xxxx_xyyyyz, g_0_y_xxxx_xyyyz, g_0_y_xxxx_xyyyzz, g_0_y_xxxx_xyyzz, g_0_y_xxxx_xyyzzz, g_0_y_xxxx_xyzzz, g_0_y_xxxx_xyzzzz, g_0_y_xxxx_xzzzz, g_0_y_xxxx_xzzzzz, g_0_y_xxxx_yyyyy, g_0_y_xxxx_yyyyz, g_0_y_xxxx_yyyzz, g_0_y_xxxx_yyzzz, g_0_y_xxxx_yzzzz, g_0_y_xxxx_zzzzz, g_0_y_xxxxx_xxxxx, g_0_y_xxxxx_xxxxy, g_0_y_xxxxx_xxxxz, g_0_y_xxxxx_xxxyy, g_0_y_xxxxx_xxxyz, g_0_y_xxxxx_xxxzz, g_0_y_xxxxx_xxyyy, g_0_y_xxxxx_xxyyz, g_0_y_xxxxx_xxyzz, g_0_y_xxxxx_xxzzz, g_0_y_xxxxx_xyyyy, g_0_y_xxxxx_xyyyz, g_0_y_xxxxx_xyyzz, g_0_y_xxxxx_xyzzz, g_0_y_xxxxx_xzzzz, g_0_y_xxxxx_yyyyy, g_0_y_xxxxx_yyyyz, g_0_y_xxxxx_yyyzz, g_0_y_xxxxx_yyzzz, g_0_y_xxxxx_yzzzz, g_0_y_xxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxx_xxxxx[k] = -g_0_y_xxxx_xxxxx[k] * ab_x + g_0_y_xxxx_xxxxxx[k];

                g_0_y_xxxxx_xxxxy[k] = -g_0_y_xxxx_xxxxy[k] * ab_x + g_0_y_xxxx_xxxxxy[k];

                g_0_y_xxxxx_xxxxz[k] = -g_0_y_xxxx_xxxxz[k] * ab_x + g_0_y_xxxx_xxxxxz[k];

                g_0_y_xxxxx_xxxyy[k] = -g_0_y_xxxx_xxxyy[k] * ab_x + g_0_y_xxxx_xxxxyy[k];

                g_0_y_xxxxx_xxxyz[k] = -g_0_y_xxxx_xxxyz[k] * ab_x + g_0_y_xxxx_xxxxyz[k];

                g_0_y_xxxxx_xxxzz[k] = -g_0_y_xxxx_xxxzz[k] * ab_x + g_0_y_xxxx_xxxxzz[k];

                g_0_y_xxxxx_xxyyy[k] = -g_0_y_xxxx_xxyyy[k] * ab_x + g_0_y_xxxx_xxxyyy[k];

                g_0_y_xxxxx_xxyyz[k] = -g_0_y_xxxx_xxyyz[k] * ab_x + g_0_y_xxxx_xxxyyz[k];

                g_0_y_xxxxx_xxyzz[k] = -g_0_y_xxxx_xxyzz[k] * ab_x + g_0_y_xxxx_xxxyzz[k];

                g_0_y_xxxxx_xxzzz[k] = -g_0_y_xxxx_xxzzz[k] * ab_x + g_0_y_xxxx_xxxzzz[k];

                g_0_y_xxxxx_xyyyy[k] = -g_0_y_xxxx_xyyyy[k] * ab_x + g_0_y_xxxx_xxyyyy[k];

                g_0_y_xxxxx_xyyyz[k] = -g_0_y_xxxx_xyyyz[k] * ab_x + g_0_y_xxxx_xxyyyz[k];

                g_0_y_xxxxx_xyyzz[k] = -g_0_y_xxxx_xyyzz[k] * ab_x + g_0_y_xxxx_xxyyzz[k];

                g_0_y_xxxxx_xyzzz[k] = -g_0_y_xxxx_xyzzz[k] * ab_x + g_0_y_xxxx_xxyzzz[k];

                g_0_y_xxxxx_xzzzz[k] = -g_0_y_xxxx_xzzzz[k] * ab_x + g_0_y_xxxx_xxzzzz[k];

                g_0_y_xxxxx_yyyyy[k] = -g_0_y_xxxx_yyyyy[k] * ab_x + g_0_y_xxxx_xyyyyy[k];

                g_0_y_xxxxx_yyyyz[k] = -g_0_y_xxxx_yyyyz[k] * ab_x + g_0_y_xxxx_xyyyyz[k];

                g_0_y_xxxxx_yyyzz[k] = -g_0_y_xxxx_yyyzz[k] * ab_x + g_0_y_xxxx_xyyyzz[k];

                g_0_y_xxxxx_yyzzz[k] = -g_0_y_xxxx_yyzzz[k] * ab_x + g_0_y_xxxx_xyyzzz[k];

                g_0_y_xxxxx_yzzzz[k] = -g_0_y_xxxx_yzzzz[k] * ab_x + g_0_y_xxxx_xyzzzz[k];

                g_0_y_xxxxx_zzzzz[k] = -g_0_y_xxxx_zzzzz[k] * ab_x + g_0_y_xxxx_xzzzzz[k];
            }

            /// Set up 462-483 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxy_xxxxx = cbuffer.data(hh_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxy = cbuffer.data(hh_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxxz = cbuffer.data(hh_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxyy = cbuffer.data(hh_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxyz = cbuffer.data(hh_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxzz = cbuffer.data(hh_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxyyy = cbuffer.data(hh_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxyyz = cbuffer.data(hh_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxyzz = cbuffer.data(hh_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxzzz = cbuffer.data(hh_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyyyy = cbuffer.data(hh_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyyyz = cbuffer.data(hh_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyyzz = cbuffer.data(hh_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyzzz = cbuffer.data(hh_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_xxxxy_xzzzz = cbuffer.data(hh_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyyyy = cbuffer.data(hh_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyyyz = cbuffer.data(hh_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyyzz = cbuffer.data(hh_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyzzz = cbuffer.data(hh_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_xxxxy_yzzzz = cbuffer.data(hh_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_xxxxy_zzzzz = cbuffer.data(hh_geom_01_off + 482 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxy_xxxxx, g_0_y_xxxxy_xxxxy, g_0_y_xxxxy_xxxxz, g_0_y_xxxxy_xxxyy, g_0_y_xxxxy_xxxyz, g_0_y_xxxxy_xxxzz, g_0_y_xxxxy_xxyyy, g_0_y_xxxxy_xxyyz, g_0_y_xxxxy_xxyzz, g_0_y_xxxxy_xxzzz, g_0_y_xxxxy_xyyyy, g_0_y_xxxxy_xyyyz, g_0_y_xxxxy_xyyzz, g_0_y_xxxxy_xyzzz, g_0_y_xxxxy_xzzzz, g_0_y_xxxxy_yyyyy, g_0_y_xxxxy_yyyyz, g_0_y_xxxxy_yyyzz, g_0_y_xxxxy_yyzzz, g_0_y_xxxxy_yzzzz, g_0_y_xxxxy_zzzzz, g_0_y_xxxy_xxxxx, g_0_y_xxxy_xxxxxx, g_0_y_xxxy_xxxxxy, g_0_y_xxxy_xxxxxz, g_0_y_xxxy_xxxxy, g_0_y_xxxy_xxxxyy, g_0_y_xxxy_xxxxyz, g_0_y_xxxy_xxxxz, g_0_y_xxxy_xxxxzz, g_0_y_xxxy_xxxyy, g_0_y_xxxy_xxxyyy, g_0_y_xxxy_xxxyyz, g_0_y_xxxy_xxxyz, g_0_y_xxxy_xxxyzz, g_0_y_xxxy_xxxzz, g_0_y_xxxy_xxxzzz, g_0_y_xxxy_xxyyy, g_0_y_xxxy_xxyyyy, g_0_y_xxxy_xxyyyz, g_0_y_xxxy_xxyyz, g_0_y_xxxy_xxyyzz, g_0_y_xxxy_xxyzz, g_0_y_xxxy_xxyzzz, g_0_y_xxxy_xxzzz, g_0_y_xxxy_xxzzzz, g_0_y_xxxy_xyyyy, g_0_y_xxxy_xyyyyy, g_0_y_xxxy_xyyyyz, g_0_y_xxxy_xyyyz, g_0_y_xxxy_xyyyzz, g_0_y_xxxy_xyyzz, g_0_y_xxxy_xyyzzz, g_0_y_xxxy_xyzzz, g_0_y_xxxy_xyzzzz, g_0_y_xxxy_xzzzz, g_0_y_xxxy_xzzzzz, g_0_y_xxxy_yyyyy, g_0_y_xxxy_yyyyz, g_0_y_xxxy_yyyzz, g_0_y_xxxy_yyzzz, g_0_y_xxxy_yzzzz, g_0_y_xxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxy_xxxxx[k] = -g_0_y_xxxy_xxxxx[k] * ab_x + g_0_y_xxxy_xxxxxx[k];

                g_0_y_xxxxy_xxxxy[k] = -g_0_y_xxxy_xxxxy[k] * ab_x + g_0_y_xxxy_xxxxxy[k];

                g_0_y_xxxxy_xxxxz[k] = -g_0_y_xxxy_xxxxz[k] * ab_x + g_0_y_xxxy_xxxxxz[k];

                g_0_y_xxxxy_xxxyy[k] = -g_0_y_xxxy_xxxyy[k] * ab_x + g_0_y_xxxy_xxxxyy[k];

                g_0_y_xxxxy_xxxyz[k] = -g_0_y_xxxy_xxxyz[k] * ab_x + g_0_y_xxxy_xxxxyz[k];

                g_0_y_xxxxy_xxxzz[k] = -g_0_y_xxxy_xxxzz[k] * ab_x + g_0_y_xxxy_xxxxzz[k];

                g_0_y_xxxxy_xxyyy[k] = -g_0_y_xxxy_xxyyy[k] * ab_x + g_0_y_xxxy_xxxyyy[k];

                g_0_y_xxxxy_xxyyz[k] = -g_0_y_xxxy_xxyyz[k] * ab_x + g_0_y_xxxy_xxxyyz[k];

                g_0_y_xxxxy_xxyzz[k] = -g_0_y_xxxy_xxyzz[k] * ab_x + g_0_y_xxxy_xxxyzz[k];

                g_0_y_xxxxy_xxzzz[k] = -g_0_y_xxxy_xxzzz[k] * ab_x + g_0_y_xxxy_xxxzzz[k];

                g_0_y_xxxxy_xyyyy[k] = -g_0_y_xxxy_xyyyy[k] * ab_x + g_0_y_xxxy_xxyyyy[k];

                g_0_y_xxxxy_xyyyz[k] = -g_0_y_xxxy_xyyyz[k] * ab_x + g_0_y_xxxy_xxyyyz[k];

                g_0_y_xxxxy_xyyzz[k] = -g_0_y_xxxy_xyyzz[k] * ab_x + g_0_y_xxxy_xxyyzz[k];

                g_0_y_xxxxy_xyzzz[k] = -g_0_y_xxxy_xyzzz[k] * ab_x + g_0_y_xxxy_xxyzzz[k];

                g_0_y_xxxxy_xzzzz[k] = -g_0_y_xxxy_xzzzz[k] * ab_x + g_0_y_xxxy_xxzzzz[k];

                g_0_y_xxxxy_yyyyy[k] = -g_0_y_xxxy_yyyyy[k] * ab_x + g_0_y_xxxy_xyyyyy[k];

                g_0_y_xxxxy_yyyyz[k] = -g_0_y_xxxy_yyyyz[k] * ab_x + g_0_y_xxxy_xyyyyz[k];

                g_0_y_xxxxy_yyyzz[k] = -g_0_y_xxxy_yyyzz[k] * ab_x + g_0_y_xxxy_xyyyzz[k];

                g_0_y_xxxxy_yyzzz[k] = -g_0_y_xxxy_yyzzz[k] * ab_x + g_0_y_xxxy_xyyzzz[k];

                g_0_y_xxxxy_yzzzz[k] = -g_0_y_xxxy_yzzzz[k] * ab_x + g_0_y_xxxy_xyzzzz[k];

                g_0_y_xxxxy_zzzzz[k] = -g_0_y_xxxy_zzzzz[k] * ab_x + g_0_y_xxxy_xzzzzz[k];
            }

            /// Set up 483-504 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxz_xxxxx = cbuffer.data(hh_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxy = cbuffer.data(hh_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxxz = cbuffer.data(hh_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxyy = cbuffer.data(hh_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxyz = cbuffer.data(hh_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxzz = cbuffer.data(hh_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxyyy = cbuffer.data(hh_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxyyz = cbuffer.data(hh_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxyzz = cbuffer.data(hh_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxzzz = cbuffer.data(hh_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyyyy = cbuffer.data(hh_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyyyz = cbuffer.data(hh_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyyzz = cbuffer.data(hh_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyzzz = cbuffer.data(hh_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_xxxxz_xzzzz = cbuffer.data(hh_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyyyy = cbuffer.data(hh_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyyyz = cbuffer.data(hh_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyyzz = cbuffer.data(hh_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyzzz = cbuffer.data(hh_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_xxxxz_yzzzz = cbuffer.data(hh_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_xxxxz_zzzzz = cbuffer.data(hh_geom_01_off + 503 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxz_xxxxx, g_0_y_xxxxz_xxxxy, g_0_y_xxxxz_xxxxz, g_0_y_xxxxz_xxxyy, g_0_y_xxxxz_xxxyz, g_0_y_xxxxz_xxxzz, g_0_y_xxxxz_xxyyy, g_0_y_xxxxz_xxyyz, g_0_y_xxxxz_xxyzz, g_0_y_xxxxz_xxzzz, g_0_y_xxxxz_xyyyy, g_0_y_xxxxz_xyyyz, g_0_y_xxxxz_xyyzz, g_0_y_xxxxz_xyzzz, g_0_y_xxxxz_xzzzz, g_0_y_xxxxz_yyyyy, g_0_y_xxxxz_yyyyz, g_0_y_xxxxz_yyyzz, g_0_y_xxxxz_yyzzz, g_0_y_xxxxz_yzzzz, g_0_y_xxxxz_zzzzz, g_0_y_xxxz_xxxxx, g_0_y_xxxz_xxxxxx, g_0_y_xxxz_xxxxxy, g_0_y_xxxz_xxxxxz, g_0_y_xxxz_xxxxy, g_0_y_xxxz_xxxxyy, g_0_y_xxxz_xxxxyz, g_0_y_xxxz_xxxxz, g_0_y_xxxz_xxxxzz, g_0_y_xxxz_xxxyy, g_0_y_xxxz_xxxyyy, g_0_y_xxxz_xxxyyz, g_0_y_xxxz_xxxyz, g_0_y_xxxz_xxxyzz, g_0_y_xxxz_xxxzz, g_0_y_xxxz_xxxzzz, g_0_y_xxxz_xxyyy, g_0_y_xxxz_xxyyyy, g_0_y_xxxz_xxyyyz, g_0_y_xxxz_xxyyz, g_0_y_xxxz_xxyyzz, g_0_y_xxxz_xxyzz, g_0_y_xxxz_xxyzzz, g_0_y_xxxz_xxzzz, g_0_y_xxxz_xxzzzz, g_0_y_xxxz_xyyyy, g_0_y_xxxz_xyyyyy, g_0_y_xxxz_xyyyyz, g_0_y_xxxz_xyyyz, g_0_y_xxxz_xyyyzz, g_0_y_xxxz_xyyzz, g_0_y_xxxz_xyyzzz, g_0_y_xxxz_xyzzz, g_0_y_xxxz_xyzzzz, g_0_y_xxxz_xzzzz, g_0_y_xxxz_xzzzzz, g_0_y_xxxz_yyyyy, g_0_y_xxxz_yyyyz, g_0_y_xxxz_yyyzz, g_0_y_xxxz_yyzzz, g_0_y_xxxz_yzzzz, g_0_y_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxz_xxxxx[k] = -g_0_y_xxxz_xxxxx[k] * ab_x + g_0_y_xxxz_xxxxxx[k];

                g_0_y_xxxxz_xxxxy[k] = -g_0_y_xxxz_xxxxy[k] * ab_x + g_0_y_xxxz_xxxxxy[k];

                g_0_y_xxxxz_xxxxz[k] = -g_0_y_xxxz_xxxxz[k] * ab_x + g_0_y_xxxz_xxxxxz[k];

                g_0_y_xxxxz_xxxyy[k] = -g_0_y_xxxz_xxxyy[k] * ab_x + g_0_y_xxxz_xxxxyy[k];

                g_0_y_xxxxz_xxxyz[k] = -g_0_y_xxxz_xxxyz[k] * ab_x + g_0_y_xxxz_xxxxyz[k];

                g_0_y_xxxxz_xxxzz[k] = -g_0_y_xxxz_xxxzz[k] * ab_x + g_0_y_xxxz_xxxxzz[k];

                g_0_y_xxxxz_xxyyy[k] = -g_0_y_xxxz_xxyyy[k] * ab_x + g_0_y_xxxz_xxxyyy[k];

                g_0_y_xxxxz_xxyyz[k] = -g_0_y_xxxz_xxyyz[k] * ab_x + g_0_y_xxxz_xxxyyz[k];

                g_0_y_xxxxz_xxyzz[k] = -g_0_y_xxxz_xxyzz[k] * ab_x + g_0_y_xxxz_xxxyzz[k];

                g_0_y_xxxxz_xxzzz[k] = -g_0_y_xxxz_xxzzz[k] * ab_x + g_0_y_xxxz_xxxzzz[k];

                g_0_y_xxxxz_xyyyy[k] = -g_0_y_xxxz_xyyyy[k] * ab_x + g_0_y_xxxz_xxyyyy[k];

                g_0_y_xxxxz_xyyyz[k] = -g_0_y_xxxz_xyyyz[k] * ab_x + g_0_y_xxxz_xxyyyz[k];

                g_0_y_xxxxz_xyyzz[k] = -g_0_y_xxxz_xyyzz[k] * ab_x + g_0_y_xxxz_xxyyzz[k];

                g_0_y_xxxxz_xyzzz[k] = -g_0_y_xxxz_xyzzz[k] * ab_x + g_0_y_xxxz_xxyzzz[k];

                g_0_y_xxxxz_xzzzz[k] = -g_0_y_xxxz_xzzzz[k] * ab_x + g_0_y_xxxz_xxzzzz[k];

                g_0_y_xxxxz_yyyyy[k] = -g_0_y_xxxz_yyyyy[k] * ab_x + g_0_y_xxxz_xyyyyy[k];

                g_0_y_xxxxz_yyyyz[k] = -g_0_y_xxxz_yyyyz[k] * ab_x + g_0_y_xxxz_xyyyyz[k];

                g_0_y_xxxxz_yyyzz[k] = -g_0_y_xxxz_yyyzz[k] * ab_x + g_0_y_xxxz_xyyyzz[k];

                g_0_y_xxxxz_yyzzz[k] = -g_0_y_xxxz_yyzzz[k] * ab_x + g_0_y_xxxz_xyyzzz[k];

                g_0_y_xxxxz_yzzzz[k] = -g_0_y_xxxz_yzzzz[k] * ab_x + g_0_y_xxxz_xyzzzz[k];

                g_0_y_xxxxz_zzzzz[k] = -g_0_y_xxxz_zzzzz[k] * ab_x + g_0_y_xxxz_xzzzzz[k];
            }

            /// Set up 504-525 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyy_xxxxx = cbuffer.data(hh_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxy = cbuffer.data(hh_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxxz = cbuffer.data(hh_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxyy = cbuffer.data(hh_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxyz = cbuffer.data(hh_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxzz = cbuffer.data(hh_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxyyy = cbuffer.data(hh_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxyyz = cbuffer.data(hh_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxyzz = cbuffer.data(hh_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxzzz = cbuffer.data(hh_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyyyy = cbuffer.data(hh_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyyyz = cbuffer.data(hh_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyyzz = cbuffer.data(hh_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyzzz = cbuffer.data(hh_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_xxxyy_xzzzz = cbuffer.data(hh_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyyyy = cbuffer.data(hh_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyyyz = cbuffer.data(hh_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyyzz = cbuffer.data(hh_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyzzz = cbuffer.data(hh_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_xxxyy_yzzzz = cbuffer.data(hh_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_xxxyy_zzzzz = cbuffer.data(hh_geom_01_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyy_xxxxx, g_0_y_xxxyy_xxxxy, g_0_y_xxxyy_xxxxz, g_0_y_xxxyy_xxxyy, g_0_y_xxxyy_xxxyz, g_0_y_xxxyy_xxxzz, g_0_y_xxxyy_xxyyy, g_0_y_xxxyy_xxyyz, g_0_y_xxxyy_xxyzz, g_0_y_xxxyy_xxzzz, g_0_y_xxxyy_xyyyy, g_0_y_xxxyy_xyyyz, g_0_y_xxxyy_xyyzz, g_0_y_xxxyy_xyzzz, g_0_y_xxxyy_xzzzz, g_0_y_xxxyy_yyyyy, g_0_y_xxxyy_yyyyz, g_0_y_xxxyy_yyyzz, g_0_y_xxxyy_yyzzz, g_0_y_xxxyy_yzzzz, g_0_y_xxxyy_zzzzz, g_0_y_xxyy_xxxxx, g_0_y_xxyy_xxxxxx, g_0_y_xxyy_xxxxxy, g_0_y_xxyy_xxxxxz, g_0_y_xxyy_xxxxy, g_0_y_xxyy_xxxxyy, g_0_y_xxyy_xxxxyz, g_0_y_xxyy_xxxxz, g_0_y_xxyy_xxxxzz, g_0_y_xxyy_xxxyy, g_0_y_xxyy_xxxyyy, g_0_y_xxyy_xxxyyz, g_0_y_xxyy_xxxyz, g_0_y_xxyy_xxxyzz, g_0_y_xxyy_xxxzz, g_0_y_xxyy_xxxzzz, g_0_y_xxyy_xxyyy, g_0_y_xxyy_xxyyyy, g_0_y_xxyy_xxyyyz, g_0_y_xxyy_xxyyz, g_0_y_xxyy_xxyyzz, g_0_y_xxyy_xxyzz, g_0_y_xxyy_xxyzzz, g_0_y_xxyy_xxzzz, g_0_y_xxyy_xxzzzz, g_0_y_xxyy_xyyyy, g_0_y_xxyy_xyyyyy, g_0_y_xxyy_xyyyyz, g_0_y_xxyy_xyyyz, g_0_y_xxyy_xyyyzz, g_0_y_xxyy_xyyzz, g_0_y_xxyy_xyyzzz, g_0_y_xxyy_xyzzz, g_0_y_xxyy_xyzzzz, g_0_y_xxyy_xzzzz, g_0_y_xxyy_xzzzzz, g_0_y_xxyy_yyyyy, g_0_y_xxyy_yyyyz, g_0_y_xxyy_yyyzz, g_0_y_xxyy_yyzzz, g_0_y_xxyy_yzzzz, g_0_y_xxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyy_xxxxx[k] = -g_0_y_xxyy_xxxxx[k] * ab_x + g_0_y_xxyy_xxxxxx[k];

                g_0_y_xxxyy_xxxxy[k] = -g_0_y_xxyy_xxxxy[k] * ab_x + g_0_y_xxyy_xxxxxy[k];

                g_0_y_xxxyy_xxxxz[k] = -g_0_y_xxyy_xxxxz[k] * ab_x + g_0_y_xxyy_xxxxxz[k];

                g_0_y_xxxyy_xxxyy[k] = -g_0_y_xxyy_xxxyy[k] * ab_x + g_0_y_xxyy_xxxxyy[k];

                g_0_y_xxxyy_xxxyz[k] = -g_0_y_xxyy_xxxyz[k] * ab_x + g_0_y_xxyy_xxxxyz[k];

                g_0_y_xxxyy_xxxzz[k] = -g_0_y_xxyy_xxxzz[k] * ab_x + g_0_y_xxyy_xxxxzz[k];

                g_0_y_xxxyy_xxyyy[k] = -g_0_y_xxyy_xxyyy[k] * ab_x + g_0_y_xxyy_xxxyyy[k];

                g_0_y_xxxyy_xxyyz[k] = -g_0_y_xxyy_xxyyz[k] * ab_x + g_0_y_xxyy_xxxyyz[k];

                g_0_y_xxxyy_xxyzz[k] = -g_0_y_xxyy_xxyzz[k] * ab_x + g_0_y_xxyy_xxxyzz[k];

                g_0_y_xxxyy_xxzzz[k] = -g_0_y_xxyy_xxzzz[k] * ab_x + g_0_y_xxyy_xxxzzz[k];

                g_0_y_xxxyy_xyyyy[k] = -g_0_y_xxyy_xyyyy[k] * ab_x + g_0_y_xxyy_xxyyyy[k];

                g_0_y_xxxyy_xyyyz[k] = -g_0_y_xxyy_xyyyz[k] * ab_x + g_0_y_xxyy_xxyyyz[k];

                g_0_y_xxxyy_xyyzz[k] = -g_0_y_xxyy_xyyzz[k] * ab_x + g_0_y_xxyy_xxyyzz[k];

                g_0_y_xxxyy_xyzzz[k] = -g_0_y_xxyy_xyzzz[k] * ab_x + g_0_y_xxyy_xxyzzz[k];

                g_0_y_xxxyy_xzzzz[k] = -g_0_y_xxyy_xzzzz[k] * ab_x + g_0_y_xxyy_xxzzzz[k];

                g_0_y_xxxyy_yyyyy[k] = -g_0_y_xxyy_yyyyy[k] * ab_x + g_0_y_xxyy_xyyyyy[k];

                g_0_y_xxxyy_yyyyz[k] = -g_0_y_xxyy_yyyyz[k] * ab_x + g_0_y_xxyy_xyyyyz[k];

                g_0_y_xxxyy_yyyzz[k] = -g_0_y_xxyy_yyyzz[k] * ab_x + g_0_y_xxyy_xyyyzz[k];

                g_0_y_xxxyy_yyzzz[k] = -g_0_y_xxyy_yyzzz[k] * ab_x + g_0_y_xxyy_xyyzzz[k];

                g_0_y_xxxyy_yzzzz[k] = -g_0_y_xxyy_yzzzz[k] * ab_x + g_0_y_xxyy_xyzzzz[k];

                g_0_y_xxxyy_zzzzz[k] = -g_0_y_xxyy_zzzzz[k] * ab_x + g_0_y_xxyy_xzzzzz[k];
            }

            /// Set up 525-546 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyz_xxxxx = cbuffer.data(hh_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxy = cbuffer.data(hh_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxxz = cbuffer.data(hh_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxyy = cbuffer.data(hh_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxyz = cbuffer.data(hh_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxzz = cbuffer.data(hh_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxyyy = cbuffer.data(hh_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxyyz = cbuffer.data(hh_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxyzz = cbuffer.data(hh_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxzzz = cbuffer.data(hh_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyyyy = cbuffer.data(hh_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyyyz = cbuffer.data(hh_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyyzz = cbuffer.data(hh_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyzzz = cbuffer.data(hh_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_xxxyz_xzzzz = cbuffer.data(hh_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyyyy = cbuffer.data(hh_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyyyz = cbuffer.data(hh_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyyzz = cbuffer.data(hh_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyzzz = cbuffer.data(hh_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_xxxyz_yzzzz = cbuffer.data(hh_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_xxxyz_zzzzz = cbuffer.data(hh_geom_01_off + 545 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyz_xxxxx, g_0_y_xxxyz_xxxxy, g_0_y_xxxyz_xxxxz, g_0_y_xxxyz_xxxyy, g_0_y_xxxyz_xxxyz, g_0_y_xxxyz_xxxzz, g_0_y_xxxyz_xxyyy, g_0_y_xxxyz_xxyyz, g_0_y_xxxyz_xxyzz, g_0_y_xxxyz_xxzzz, g_0_y_xxxyz_xyyyy, g_0_y_xxxyz_xyyyz, g_0_y_xxxyz_xyyzz, g_0_y_xxxyz_xyzzz, g_0_y_xxxyz_xzzzz, g_0_y_xxxyz_yyyyy, g_0_y_xxxyz_yyyyz, g_0_y_xxxyz_yyyzz, g_0_y_xxxyz_yyzzz, g_0_y_xxxyz_yzzzz, g_0_y_xxxyz_zzzzz, g_0_y_xxyz_xxxxx, g_0_y_xxyz_xxxxxx, g_0_y_xxyz_xxxxxy, g_0_y_xxyz_xxxxxz, g_0_y_xxyz_xxxxy, g_0_y_xxyz_xxxxyy, g_0_y_xxyz_xxxxyz, g_0_y_xxyz_xxxxz, g_0_y_xxyz_xxxxzz, g_0_y_xxyz_xxxyy, g_0_y_xxyz_xxxyyy, g_0_y_xxyz_xxxyyz, g_0_y_xxyz_xxxyz, g_0_y_xxyz_xxxyzz, g_0_y_xxyz_xxxzz, g_0_y_xxyz_xxxzzz, g_0_y_xxyz_xxyyy, g_0_y_xxyz_xxyyyy, g_0_y_xxyz_xxyyyz, g_0_y_xxyz_xxyyz, g_0_y_xxyz_xxyyzz, g_0_y_xxyz_xxyzz, g_0_y_xxyz_xxyzzz, g_0_y_xxyz_xxzzz, g_0_y_xxyz_xxzzzz, g_0_y_xxyz_xyyyy, g_0_y_xxyz_xyyyyy, g_0_y_xxyz_xyyyyz, g_0_y_xxyz_xyyyz, g_0_y_xxyz_xyyyzz, g_0_y_xxyz_xyyzz, g_0_y_xxyz_xyyzzz, g_0_y_xxyz_xyzzz, g_0_y_xxyz_xyzzzz, g_0_y_xxyz_xzzzz, g_0_y_xxyz_xzzzzz, g_0_y_xxyz_yyyyy, g_0_y_xxyz_yyyyz, g_0_y_xxyz_yyyzz, g_0_y_xxyz_yyzzz, g_0_y_xxyz_yzzzz, g_0_y_xxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyz_xxxxx[k] = -g_0_y_xxyz_xxxxx[k] * ab_x + g_0_y_xxyz_xxxxxx[k];

                g_0_y_xxxyz_xxxxy[k] = -g_0_y_xxyz_xxxxy[k] * ab_x + g_0_y_xxyz_xxxxxy[k];

                g_0_y_xxxyz_xxxxz[k] = -g_0_y_xxyz_xxxxz[k] * ab_x + g_0_y_xxyz_xxxxxz[k];

                g_0_y_xxxyz_xxxyy[k] = -g_0_y_xxyz_xxxyy[k] * ab_x + g_0_y_xxyz_xxxxyy[k];

                g_0_y_xxxyz_xxxyz[k] = -g_0_y_xxyz_xxxyz[k] * ab_x + g_0_y_xxyz_xxxxyz[k];

                g_0_y_xxxyz_xxxzz[k] = -g_0_y_xxyz_xxxzz[k] * ab_x + g_0_y_xxyz_xxxxzz[k];

                g_0_y_xxxyz_xxyyy[k] = -g_0_y_xxyz_xxyyy[k] * ab_x + g_0_y_xxyz_xxxyyy[k];

                g_0_y_xxxyz_xxyyz[k] = -g_0_y_xxyz_xxyyz[k] * ab_x + g_0_y_xxyz_xxxyyz[k];

                g_0_y_xxxyz_xxyzz[k] = -g_0_y_xxyz_xxyzz[k] * ab_x + g_0_y_xxyz_xxxyzz[k];

                g_0_y_xxxyz_xxzzz[k] = -g_0_y_xxyz_xxzzz[k] * ab_x + g_0_y_xxyz_xxxzzz[k];

                g_0_y_xxxyz_xyyyy[k] = -g_0_y_xxyz_xyyyy[k] * ab_x + g_0_y_xxyz_xxyyyy[k];

                g_0_y_xxxyz_xyyyz[k] = -g_0_y_xxyz_xyyyz[k] * ab_x + g_0_y_xxyz_xxyyyz[k];

                g_0_y_xxxyz_xyyzz[k] = -g_0_y_xxyz_xyyzz[k] * ab_x + g_0_y_xxyz_xxyyzz[k];

                g_0_y_xxxyz_xyzzz[k] = -g_0_y_xxyz_xyzzz[k] * ab_x + g_0_y_xxyz_xxyzzz[k];

                g_0_y_xxxyz_xzzzz[k] = -g_0_y_xxyz_xzzzz[k] * ab_x + g_0_y_xxyz_xxzzzz[k];

                g_0_y_xxxyz_yyyyy[k] = -g_0_y_xxyz_yyyyy[k] * ab_x + g_0_y_xxyz_xyyyyy[k];

                g_0_y_xxxyz_yyyyz[k] = -g_0_y_xxyz_yyyyz[k] * ab_x + g_0_y_xxyz_xyyyyz[k];

                g_0_y_xxxyz_yyyzz[k] = -g_0_y_xxyz_yyyzz[k] * ab_x + g_0_y_xxyz_xyyyzz[k];

                g_0_y_xxxyz_yyzzz[k] = -g_0_y_xxyz_yyzzz[k] * ab_x + g_0_y_xxyz_xyyzzz[k];

                g_0_y_xxxyz_yzzzz[k] = -g_0_y_xxyz_yzzzz[k] * ab_x + g_0_y_xxyz_xyzzzz[k];

                g_0_y_xxxyz_zzzzz[k] = -g_0_y_xxyz_zzzzz[k] * ab_x + g_0_y_xxyz_xzzzzz[k];
            }

            /// Set up 546-567 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzz_xxxxx = cbuffer.data(hh_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxy = cbuffer.data(hh_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxxz = cbuffer.data(hh_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxyy = cbuffer.data(hh_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxyz = cbuffer.data(hh_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxzz = cbuffer.data(hh_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxyyy = cbuffer.data(hh_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxyyz = cbuffer.data(hh_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxyzz = cbuffer.data(hh_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxzzz = cbuffer.data(hh_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyyyy = cbuffer.data(hh_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyyyz = cbuffer.data(hh_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyyzz = cbuffer.data(hh_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyzzz = cbuffer.data(hh_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_y_xxxzz_xzzzz = cbuffer.data(hh_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyyyy = cbuffer.data(hh_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyyyz = cbuffer.data(hh_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyyzz = cbuffer.data(hh_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyzzz = cbuffer.data(hh_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_y_xxxzz_yzzzz = cbuffer.data(hh_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_y_xxxzz_zzzzz = cbuffer.data(hh_geom_01_off + 566 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzz_xxxxx, g_0_y_xxxzz_xxxxy, g_0_y_xxxzz_xxxxz, g_0_y_xxxzz_xxxyy, g_0_y_xxxzz_xxxyz, g_0_y_xxxzz_xxxzz, g_0_y_xxxzz_xxyyy, g_0_y_xxxzz_xxyyz, g_0_y_xxxzz_xxyzz, g_0_y_xxxzz_xxzzz, g_0_y_xxxzz_xyyyy, g_0_y_xxxzz_xyyyz, g_0_y_xxxzz_xyyzz, g_0_y_xxxzz_xyzzz, g_0_y_xxxzz_xzzzz, g_0_y_xxxzz_yyyyy, g_0_y_xxxzz_yyyyz, g_0_y_xxxzz_yyyzz, g_0_y_xxxzz_yyzzz, g_0_y_xxxzz_yzzzz, g_0_y_xxxzz_zzzzz, g_0_y_xxzz_xxxxx, g_0_y_xxzz_xxxxxx, g_0_y_xxzz_xxxxxy, g_0_y_xxzz_xxxxxz, g_0_y_xxzz_xxxxy, g_0_y_xxzz_xxxxyy, g_0_y_xxzz_xxxxyz, g_0_y_xxzz_xxxxz, g_0_y_xxzz_xxxxzz, g_0_y_xxzz_xxxyy, g_0_y_xxzz_xxxyyy, g_0_y_xxzz_xxxyyz, g_0_y_xxzz_xxxyz, g_0_y_xxzz_xxxyzz, g_0_y_xxzz_xxxzz, g_0_y_xxzz_xxxzzz, g_0_y_xxzz_xxyyy, g_0_y_xxzz_xxyyyy, g_0_y_xxzz_xxyyyz, g_0_y_xxzz_xxyyz, g_0_y_xxzz_xxyyzz, g_0_y_xxzz_xxyzz, g_0_y_xxzz_xxyzzz, g_0_y_xxzz_xxzzz, g_0_y_xxzz_xxzzzz, g_0_y_xxzz_xyyyy, g_0_y_xxzz_xyyyyy, g_0_y_xxzz_xyyyyz, g_0_y_xxzz_xyyyz, g_0_y_xxzz_xyyyzz, g_0_y_xxzz_xyyzz, g_0_y_xxzz_xyyzzz, g_0_y_xxzz_xyzzz, g_0_y_xxzz_xyzzzz, g_0_y_xxzz_xzzzz, g_0_y_xxzz_xzzzzz, g_0_y_xxzz_yyyyy, g_0_y_xxzz_yyyyz, g_0_y_xxzz_yyyzz, g_0_y_xxzz_yyzzz, g_0_y_xxzz_yzzzz, g_0_y_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzz_xxxxx[k] = -g_0_y_xxzz_xxxxx[k] * ab_x + g_0_y_xxzz_xxxxxx[k];

                g_0_y_xxxzz_xxxxy[k] = -g_0_y_xxzz_xxxxy[k] * ab_x + g_0_y_xxzz_xxxxxy[k];

                g_0_y_xxxzz_xxxxz[k] = -g_0_y_xxzz_xxxxz[k] * ab_x + g_0_y_xxzz_xxxxxz[k];

                g_0_y_xxxzz_xxxyy[k] = -g_0_y_xxzz_xxxyy[k] * ab_x + g_0_y_xxzz_xxxxyy[k];

                g_0_y_xxxzz_xxxyz[k] = -g_0_y_xxzz_xxxyz[k] * ab_x + g_0_y_xxzz_xxxxyz[k];

                g_0_y_xxxzz_xxxzz[k] = -g_0_y_xxzz_xxxzz[k] * ab_x + g_0_y_xxzz_xxxxzz[k];

                g_0_y_xxxzz_xxyyy[k] = -g_0_y_xxzz_xxyyy[k] * ab_x + g_0_y_xxzz_xxxyyy[k];

                g_0_y_xxxzz_xxyyz[k] = -g_0_y_xxzz_xxyyz[k] * ab_x + g_0_y_xxzz_xxxyyz[k];

                g_0_y_xxxzz_xxyzz[k] = -g_0_y_xxzz_xxyzz[k] * ab_x + g_0_y_xxzz_xxxyzz[k];

                g_0_y_xxxzz_xxzzz[k] = -g_0_y_xxzz_xxzzz[k] * ab_x + g_0_y_xxzz_xxxzzz[k];

                g_0_y_xxxzz_xyyyy[k] = -g_0_y_xxzz_xyyyy[k] * ab_x + g_0_y_xxzz_xxyyyy[k];

                g_0_y_xxxzz_xyyyz[k] = -g_0_y_xxzz_xyyyz[k] * ab_x + g_0_y_xxzz_xxyyyz[k];

                g_0_y_xxxzz_xyyzz[k] = -g_0_y_xxzz_xyyzz[k] * ab_x + g_0_y_xxzz_xxyyzz[k];

                g_0_y_xxxzz_xyzzz[k] = -g_0_y_xxzz_xyzzz[k] * ab_x + g_0_y_xxzz_xxyzzz[k];

                g_0_y_xxxzz_xzzzz[k] = -g_0_y_xxzz_xzzzz[k] * ab_x + g_0_y_xxzz_xxzzzz[k];

                g_0_y_xxxzz_yyyyy[k] = -g_0_y_xxzz_yyyyy[k] * ab_x + g_0_y_xxzz_xyyyyy[k];

                g_0_y_xxxzz_yyyyz[k] = -g_0_y_xxzz_yyyyz[k] * ab_x + g_0_y_xxzz_xyyyyz[k];

                g_0_y_xxxzz_yyyzz[k] = -g_0_y_xxzz_yyyzz[k] * ab_x + g_0_y_xxzz_xyyyzz[k];

                g_0_y_xxxzz_yyzzz[k] = -g_0_y_xxzz_yyzzz[k] * ab_x + g_0_y_xxzz_xyyzzz[k];

                g_0_y_xxxzz_yzzzz[k] = -g_0_y_xxzz_yzzzz[k] * ab_x + g_0_y_xxzz_xyzzzz[k];

                g_0_y_xxxzz_zzzzz[k] = -g_0_y_xxzz_zzzzz[k] * ab_x + g_0_y_xxzz_xzzzzz[k];
            }

            /// Set up 567-588 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyy_xxxxx = cbuffer.data(hh_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxy = cbuffer.data(hh_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxxz = cbuffer.data(hh_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxyy = cbuffer.data(hh_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxyz = cbuffer.data(hh_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxzz = cbuffer.data(hh_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxyyy = cbuffer.data(hh_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxyyz = cbuffer.data(hh_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxyzz = cbuffer.data(hh_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxzzz = cbuffer.data(hh_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyyyy = cbuffer.data(hh_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyyyz = cbuffer.data(hh_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyyzz = cbuffer.data(hh_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyzzz = cbuffer.data(hh_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_y_xxyyy_xzzzz = cbuffer.data(hh_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyyyy = cbuffer.data(hh_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyyyz = cbuffer.data(hh_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyyzz = cbuffer.data(hh_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyzzz = cbuffer.data(hh_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_y_xxyyy_yzzzz = cbuffer.data(hh_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_y_xxyyy_zzzzz = cbuffer.data(hh_geom_01_off + 587 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyy_xxxxx, g_0_y_xxyyy_xxxxy, g_0_y_xxyyy_xxxxz, g_0_y_xxyyy_xxxyy, g_0_y_xxyyy_xxxyz, g_0_y_xxyyy_xxxzz, g_0_y_xxyyy_xxyyy, g_0_y_xxyyy_xxyyz, g_0_y_xxyyy_xxyzz, g_0_y_xxyyy_xxzzz, g_0_y_xxyyy_xyyyy, g_0_y_xxyyy_xyyyz, g_0_y_xxyyy_xyyzz, g_0_y_xxyyy_xyzzz, g_0_y_xxyyy_xzzzz, g_0_y_xxyyy_yyyyy, g_0_y_xxyyy_yyyyz, g_0_y_xxyyy_yyyzz, g_0_y_xxyyy_yyzzz, g_0_y_xxyyy_yzzzz, g_0_y_xxyyy_zzzzz, g_0_y_xyyy_xxxxx, g_0_y_xyyy_xxxxxx, g_0_y_xyyy_xxxxxy, g_0_y_xyyy_xxxxxz, g_0_y_xyyy_xxxxy, g_0_y_xyyy_xxxxyy, g_0_y_xyyy_xxxxyz, g_0_y_xyyy_xxxxz, g_0_y_xyyy_xxxxzz, g_0_y_xyyy_xxxyy, g_0_y_xyyy_xxxyyy, g_0_y_xyyy_xxxyyz, g_0_y_xyyy_xxxyz, g_0_y_xyyy_xxxyzz, g_0_y_xyyy_xxxzz, g_0_y_xyyy_xxxzzz, g_0_y_xyyy_xxyyy, g_0_y_xyyy_xxyyyy, g_0_y_xyyy_xxyyyz, g_0_y_xyyy_xxyyz, g_0_y_xyyy_xxyyzz, g_0_y_xyyy_xxyzz, g_0_y_xyyy_xxyzzz, g_0_y_xyyy_xxzzz, g_0_y_xyyy_xxzzzz, g_0_y_xyyy_xyyyy, g_0_y_xyyy_xyyyyy, g_0_y_xyyy_xyyyyz, g_0_y_xyyy_xyyyz, g_0_y_xyyy_xyyyzz, g_0_y_xyyy_xyyzz, g_0_y_xyyy_xyyzzz, g_0_y_xyyy_xyzzz, g_0_y_xyyy_xyzzzz, g_0_y_xyyy_xzzzz, g_0_y_xyyy_xzzzzz, g_0_y_xyyy_yyyyy, g_0_y_xyyy_yyyyz, g_0_y_xyyy_yyyzz, g_0_y_xyyy_yyzzz, g_0_y_xyyy_yzzzz, g_0_y_xyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyy_xxxxx[k] = -g_0_y_xyyy_xxxxx[k] * ab_x + g_0_y_xyyy_xxxxxx[k];

                g_0_y_xxyyy_xxxxy[k] = -g_0_y_xyyy_xxxxy[k] * ab_x + g_0_y_xyyy_xxxxxy[k];

                g_0_y_xxyyy_xxxxz[k] = -g_0_y_xyyy_xxxxz[k] * ab_x + g_0_y_xyyy_xxxxxz[k];

                g_0_y_xxyyy_xxxyy[k] = -g_0_y_xyyy_xxxyy[k] * ab_x + g_0_y_xyyy_xxxxyy[k];

                g_0_y_xxyyy_xxxyz[k] = -g_0_y_xyyy_xxxyz[k] * ab_x + g_0_y_xyyy_xxxxyz[k];

                g_0_y_xxyyy_xxxzz[k] = -g_0_y_xyyy_xxxzz[k] * ab_x + g_0_y_xyyy_xxxxzz[k];

                g_0_y_xxyyy_xxyyy[k] = -g_0_y_xyyy_xxyyy[k] * ab_x + g_0_y_xyyy_xxxyyy[k];

                g_0_y_xxyyy_xxyyz[k] = -g_0_y_xyyy_xxyyz[k] * ab_x + g_0_y_xyyy_xxxyyz[k];

                g_0_y_xxyyy_xxyzz[k] = -g_0_y_xyyy_xxyzz[k] * ab_x + g_0_y_xyyy_xxxyzz[k];

                g_0_y_xxyyy_xxzzz[k] = -g_0_y_xyyy_xxzzz[k] * ab_x + g_0_y_xyyy_xxxzzz[k];

                g_0_y_xxyyy_xyyyy[k] = -g_0_y_xyyy_xyyyy[k] * ab_x + g_0_y_xyyy_xxyyyy[k];

                g_0_y_xxyyy_xyyyz[k] = -g_0_y_xyyy_xyyyz[k] * ab_x + g_0_y_xyyy_xxyyyz[k];

                g_0_y_xxyyy_xyyzz[k] = -g_0_y_xyyy_xyyzz[k] * ab_x + g_0_y_xyyy_xxyyzz[k];

                g_0_y_xxyyy_xyzzz[k] = -g_0_y_xyyy_xyzzz[k] * ab_x + g_0_y_xyyy_xxyzzz[k];

                g_0_y_xxyyy_xzzzz[k] = -g_0_y_xyyy_xzzzz[k] * ab_x + g_0_y_xyyy_xxzzzz[k];

                g_0_y_xxyyy_yyyyy[k] = -g_0_y_xyyy_yyyyy[k] * ab_x + g_0_y_xyyy_xyyyyy[k];

                g_0_y_xxyyy_yyyyz[k] = -g_0_y_xyyy_yyyyz[k] * ab_x + g_0_y_xyyy_xyyyyz[k];

                g_0_y_xxyyy_yyyzz[k] = -g_0_y_xyyy_yyyzz[k] * ab_x + g_0_y_xyyy_xyyyzz[k];

                g_0_y_xxyyy_yyzzz[k] = -g_0_y_xyyy_yyzzz[k] * ab_x + g_0_y_xyyy_xyyzzz[k];

                g_0_y_xxyyy_yzzzz[k] = -g_0_y_xyyy_yzzzz[k] * ab_x + g_0_y_xyyy_xyzzzz[k];

                g_0_y_xxyyy_zzzzz[k] = -g_0_y_xyyy_zzzzz[k] * ab_x + g_0_y_xyyy_xzzzzz[k];
            }

            /// Set up 588-609 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyz_xxxxx = cbuffer.data(hh_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxy = cbuffer.data(hh_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxxz = cbuffer.data(hh_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxyy = cbuffer.data(hh_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxyz = cbuffer.data(hh_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxzz = cbuffer.data(hh_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxyyy = cbuffer.data(hh_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxyyz = cbuffer.data(hh_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxyzz = cbuffer.data(hh_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxzzz = cbuffer.data(hh_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyyyy = cbuffer.data(hh_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyyyz = cbuffer.data(hh_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyyzz = cbuffer.data(hh_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyzzz = cbuffer.data(hh_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_xxyyz_xzzzz = cbuffer.data(hh_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyyyy = cbuffer.data(hh_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyyyz = cbuffer.data(hh_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyyzz = cbuffer.data(hh_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyzzz = cbuffer.data(hh_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_xxyyz_yzzzz = cbuffer.data(hh_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_xxyyz_zzzzz = cbuffer.data(hh_geom_01_off + 608 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyz_xxxxx, g_0_y_xxyyz_xxxxy, g_0_y_xxyyz_xxxxz, g_0_y_xxyyz_xxxyy, g_0_y_xxyyz_xxxyz, g_0_y_xxyyz_xxxzz, g_0_y_xxyyz_xxyyy, g_0_y_xxyyz_xxyyz, g_0_y_xxyyz_xxyzz, g_0_y_xxyyz_xxzzz, g_0_y_xxyyz_xyyyy, g_0_y_xxyyz_xyyyz, g_0_y_xxyyz_xyyzz, g_0_y_xxyyz_xyzzz, g_0_y_xxyyz_xzzzz, g_0_y_xxyyz_yyyyy, g_0_y_xxyyz_yyyyz, g_0_y_xxyyz_yyyzz, g_0_y_xxyyz_yyzzz, g_0_y_xxyyz_yzzzz, g_0_y_xxyyz_zzzzz, g_0_y_xyyz_xxxxx, g_0_y_xyyz_xxxxxx, g_0_y_xyyz_xxxxxy, g_0_y_xyyz_xxxxxz, g_0_y_xyyz_xxxxy, g_0_y_xyyz_xxxxyy, g_0_y_xyyz_xxxxyz, g_0_y_xyyz_xxxxz, g_0_y_xyyz_xxxxzz, g_0_y_xyyz_xxxyy, g_0_y_xyyz_xxxyyy, g_0_y_xyyz_xxxyyz, g_0_y_xyyz_xxxyz, g_0_y_xyyz_xxxyzz, g_0_y_xyyz_xxxzz, g_0_y_xyyz_xxxzzz, g_0_y_xyyz_xxyyy, g_0_y_xyyz_xxyyyy, g_0_y_xyyz_xxyyyz, g_0_y_xyyz_xxyyz, g_0_y_xyyz_xxyyzz, g_0_y_xyyz_xxyzz, g_0_y_xyyz_xxyzzz, g_0_y_xyyz_xxzzz, g_0_y_xyyz_xxzzzz, g_0_y_xyyz_xyyyy, g_0_y_xyyz_xyyyyy, g_0_y_xyyz_xyyyyz, g_0_y_xyyz_xyyyz, g_0_y_xyyz_xyyyzz, g_0_y_xyyz_xyyzz, g_0_y_xyyz_xyyzzz, g_0_y_xyyz_xyzzz, g_0_y_xyyz_xyzzzz, g_0_y_xyyz_xzzzz, g_0_y_xyyz_xzzzzz, g_0_y_xyyz_yyyyy, g_0_y_xyyz_yyyyz, g_0_y_xyyz_yyyzz, g_0_y_xyyz_yyzzz, g_0_y_xyyz_yzzzz, g_0_y_xyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyz_xxxxx[k] = -g_0_y_xyyz_xxxxx[k] * ab_x + g_0_y_xyyz_xxxxxx[k];

                g_0_y_xxyyz_xxxxy[k] = -g_0_y_xyyz_xxxxy[k] * ab_x + g_0_y_xyyz_xxxxxy[k];

                g_0_y_xxyyz_xxxxz[k] = -g_0_y_xyyz_xxxxz[k] * ab_x + g_0_y_xyyz_xxxxxz[k];

                g_0_y_xxyyz_xxxyy[k] = -g_0_y_xyyz_xxxyy[k] * ab_x + g_0_y_xyyz_xxxxyy[k];

                g_0_y_xxyyz_xxxyz[k] = -g_0_y_xyyz_xxxyz[k] * ab_x + g_0_y_xyyz_xxxxyz[k];

                g_0_y_xxyyz_xxxzz[k] = -g_0_y_xyyz_xxxzz[k] * ab_x + g_0_y_xyyz_xxxxzz[k];

                g_0_y_xxyyz_xxyyy[k] = -g_0_y_xyyz_xxyyy[k] * ab_x + g_0_y_xyyz_xxxyyy[k];

                g_0_y_xxyyz_xxyyz[k] = -g_0_y_xyyz_xxyyz[k] * ab_x + g_0_y_xyyz_xxxyyz[k];

                g_0_y_xxyyz_xxyzz[k] = -g_0_y_xyyz_xxyzz[k] * ab_x + g_0_y_xyyz_xxxyzz[k];

                g_0_y_xxyyz_xxzzz[k] = -g_0_y_xyyz_xxzzz[k] * ab_x + g_0_y_xyyz_xxxzzz[k];

                g_0_y_xxyyz_xyyyy[k] = -g_0_y_xyyz_xyyyy[k] * ab_x + g_0_y_xyyz_xxyyyy[k];

                g_0_y_xxyyz_xyyyz[k] = -g_0_y_xyyz_xyyyz[k] * ab_x + g_0_y_xyyz_xxyyyz[k];

                g_0_y_xxyyz_xyyzz[k] = -g_0_y_xyyz_xyyzz[k] * ab_x + g_0_y_xyyz_xxyyzz[k];

                g_0_y_xxyyz_xyzzz[k] = -g_0_y_xyyz_xyzzz[k] * ab_x + g_0_y_xyyz_xxyzzz[k];

                g_0_y_xxyyz_xzzzz[k] = -g_0_y_xyyz_xzzzz[k] * ab_x + g_0_y_xyyz_xxzzzz[k];

                g_0_y_xxyyz_yyyyy[k] = -g_0_y_xyyz_yyyyy[k] * ab_x + g_0_y_xyyz_xyyyyy[k];

                g_0_y_xxyyz_yyyyz[k] = -g_0_y_xyyz_yyyyz[k] * ab_x + g_0_y_xyyz_xyyyyz[k];

                g_0_y_xxyyz_yyyzz[k] = -g_0_y_xyyz_yyyzz[k] * ab_x + g_0_y_xyyz_xyyyzz[k];

                g_0_y_xxyyz_yyzzz[k] = -g_0_y_xyyz_yyzzz[k] * ab_x + g_0_y_xyyz_xyyzzz[k];

                g_0_y_xxyyz_yzzzz[k] = -g_0_y_xyyz_yzzzz[k] * ab_x + g_0_y_xyyz_xyzzzz[k];

                g_0_y_xxyyz_zzzzz[k] = -g_0_y_xyyz_zzzzz[k] * ab_x + g_0_y_xyyz_xzzzzz[k];
            }

            /// Set up 609-630 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzz_xxxxx = cbuffer.data(hh_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxy = cbuffer.data(hh_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxxz = cbuffer.data(hh_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxyy = cbuffer.data(hh_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxyz = cbuffer.data(hh_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxzz = cbuffer.data(hh_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxyyy = cbuffer.data(hh_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxyyz = cbuffer.data(hh_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxyzz = cbuffer.data(hh_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxzzz = cbuffer.data(hh_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyyyy = cbuffer.data(hh_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyyyz = cbuffer.data(hh_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyyzz = cbuffer.data(hh_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyzzz = cbuffer.data(hh_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_xxyzz_xzzzz = cbuffer.data(hh_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyyyy = cbuffer.data(hh_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyyyz = cbuffer.data(hh_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyyzz = cbuffer.data(hh_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyzzz = cbuffer.data(hh_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_xxyzz_yzzzz = cbuffer.data(hh_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_xxyzz_zzzzz = cbuffer.data(hh_geom_01_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzz_xxxxx, g_0_y_xxyzz_xxxxy, g_0_y_xxyzz_xxxxz, g_0_y_xxyzz_xxxyy, g_0_y_xxyzz_xxxyz, g_0_y_xxyzz_xxxzz, g_0_y_xxyzz_xxyyy, g_0_y_xxyzz_xxyyz, g_0_y_xxyzz_xxyzz, g_0_y_xxyzz_xxzzz, g_0_y_xxyzz_xyyyy, g_0_y_xxyzz_xyyyz, g_0_y_xxyzz_xyyzz, g_0_y_xxyzz_xyzzz, g_0_y_xxyzz_xzzzz, g_0_y_xxyzz_yyyyy, g_0_y_xxyzz_yyyyz, g_0_y_xxyzz_yyyzz, g_0_y_xxyzz_yyzzz, g_0_y_xxyzz_yzzzz, g_0_y_xxyzz_zzzzz, g_0_y_xyzz_xxxxx, g_0_y_xyzz_xxxxxx, g_0_y_xyzz_xxxxxy, g_0_y_xyzz_xxxxxz, g_0_y_xyzz_xxxxy, g_0_y_xyzz_xxxxyy, g_0_y_xyzz_xxxxyz, g_0_y_xyzz_xxxxz, g_0_y_xyzz_xxxxzz, g_0_y_xyzz_xxxyy, g_0_y_xyzz_xxxyyy, g_0_y_xyzz_xxxyyz, g_0_y_xyzz_xxxyz, g_0_y_xyzz_xxxyzz, g_0_y_xyzz_xxxzz, g_0_y_xyzz_xxxzzz, g_0_y_xyzz_xxyyy, g_0_y_xyzz_xxyyyy, g_0_y_xyzz_xxyyyz, g_0_y_xyzz_xxyyz, g_0_y_xyzz_xxyyzz, g_0_y_xyzz_xxyzz, g_0_y_xyzz_xxyzzz, g_0_y_xyzz_xxzzz, g_0_y_xyzz_xxzzzz, g_0_y_xyzz_xyyyy, g_0_y_xyzz_xyyyyy, g_0_y_xyzz_xyyyyz, g_0_y_xyzz_xyyyz, g_0_y_xyzz_xyyyzz, g_0_y_xyzz_xyyzz, g_0_y_xyzz_xyyzzz, g_0_y_xyzz_xyzzz, g_0_y_xyzz_xyzzzz, g_0_y_xyzz_xzzzz, g_0_y_xyzz_xzzzzz, g_0_y_xyzz_yyyyy, g_0_y_xyzz_yyyyz, g_0_y_xyzz_yyyzz, g_0_y_xyzz_yyzzz, g_0_y_xyzz_yzzzz, g_0_y_xyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzz_xxxxx[k] = -g_0_y_xyzz_xxxxx[k] * ab_x + g_0_y_xyzz_xxxxxx[k];

                g_0_y_xxyzz_xxxxy[k] = -g_0_y_xyzz_xxxxy[k] * ab_x + g_0_y_xyzz_xxxxxy[k];

                g_0_y_xxyzz_xxxxz[k] = -g_0_y_xyzz_xxxxz[k] * ab_x + g_0_y_xyzz_xxxxxz[k];

                g_0_y_xxyzz_xxxyy[k] = -g_0_y_xyzz_xxxyy[k] * ab_x + g_0_y_xyzz_xxxxyy[k];

                g_0_y_xxyzz_xxxyz[k] = -g_0_y_xyzz_xxxyz[k] * ab_x + g_0_y_xyzz_xxxxyz[k];

                g_0_y_xxyzz_xxxzz[k] = -g_0_y_xyzz_xxxzz[k] * ab_x + g_0_y_xyzz_xxxxzz[k];

                g_0_y_xxyzz_xxyyy[k] = -g_0_y_xyzz_xxyyy[k] * ab_x + g_0_y_xyzz_xxxyyy[k];

                g_0_y_xxyzz_xxyyz[k] = -g_0_y_xyzz_xxyyz[k] * ab_x + g_0_y_xyzz_xxxyyz[k];

                g_0_y_xxyzz_xxyzz[k] = -g_0_y_xyzz_xxyzz[k] * ab_x + g_0_y_xyzz_xxxyzz[k];

                g_0_y_xxyzz_xxzzz[k] = -g_0_y_xyzz_xxzzz[k] * ab_x + g_0_y_xyzz_xxxzzz[k];

                g_0_y_xxyzz_xyyyy[k] = -g_0_y_xyzz_xyyyy[k] * ab_x + g_0_y_xyzz_xxyyyy[k];

                g_0_y_xxyzz_xyyyz[k] = -g_0_y_xyzz_xyyyz[k] * ab_x + g_0_y_xyzz_xxyyyz[k];

                g_0_y_xxyzz_xyyzz[k] = -g_0_y_xyzz_xyyzz[k] * ab_x + g_0_y_xyzz_xxyyzz[k];

                g_0_y_xxyzz_xyzzz[k] = -g_0_y_xyzz_xyzzz[k] * ab_x + g_0_y_xyzz_xxyzzz[k];

                g_0_y_xxyzz_xzzzz[k] = -g_0_y_xyzz_xzzzz[k] * ab_x + g_0_y_xyzz_xxzzzz[k];

                g_0_y_xxyzz_yyyyy[k] = -g_0_y_xyzz_yyyyy[k] * ab_x + g_0_y_xyzz_xyyyyy[k];

                g_0_y_xxyzz_yyyyz[k] = -g_0_y_xyzz_yyyyz[k] * ab_x + g_0_y_xyzz_xyyyyz[k];

                g_0_y_xxyzz_yyyzz[k] = -g_0_y_xyzz_yyyzz[k] * ab_x + g_0_y_xyzz_xyyyzz[k];

                g_0_y_xxyzz_yyzzz[k] = -g_0_y_xyzz_yyzzz[k] * ab_x + g_0_y_xyzz_xyyzzz[k];

                g_0_y_xxyzz_yzzzz[k] = -g_0_y_xyzz_yzzzz[k] * ab_x + g_0_y_xyzz_xyzzzz[k];

                g_0_y_xxyzz_zzzzz[k] = -g_0_y_xyzz_zzzzz[k] * ab_x + g_0_y_xyzz_xzzzzz[k];
            }

            /// Set up 630-651 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzz_xxxxx = cbuffer.data(hh_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxy = cbuffer.data(hh_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxxz = cbuffer.data(hh_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxyy = cbuffer.data(hh_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxyz = cbuffer.data(hh_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxzz = cbuffer.data(hh_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxyyy = cbuffer.data(hh_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxyyz = cbuffer.data(hh_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxyzz = cbuffer.data(hh_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxzzz = cbuffer.data(hh_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyyyy = cbuffer.data(hh_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyyyz = cbuffer.data(hh_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyyzz = cbuffer.data(hh_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyzzz = cbuffer.data(hh_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_y_xxzzz_xzzzz = cbuffer.data(hh_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyyyy = cbuffer.data(hh_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyyyz = cbuffer.data(hh_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyyzz = cbuffer.data(hh_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyzzz = cbuffer.data(hh_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_y_xxzzz_yzzzz = cbuffer.data(hh_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_y_xxzzz_zzzzz = cbuffer.data(hh_geom_01_off + 650 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzz_xxxxx, g_0_y_xxzzz_xxxxy, g_0_y_xxzzz_xxxxz, g_0_y_xxzzz_xxxyy, g_0_y_xxzzz_xxxyz, g_0_y_xxzzz_xxxzz, g_0_y_xxzzz_xxyyy, g_0_y_xxzzz_xxyyz, g_0_y_xxzzz_xxyzz, g_0_y_xxzzz_xxzzz, g_0_y_xxzzz_xyyyy, g_0_y_xxzzz_xyyyz, g_0_y_xxzzz_xyyzz, g_0_y_xxzzz_xyzzz, g_0_y_xxzzz_xzzzz, g_0_y_xxzzz_yyyyy, g_0_y_xxzzz_yyyyz, g_0_y_xxzzz_yyyzz, g_0_y_xxzzz_yyzzz, g_0_y_xxzzz_yzzzz, g_0_y_xxzzz_zzzzz, g_0_y_xzzz_xxxxx, g_0_y_xzzz_xxxxxx, g_0_y_xzzz_xxxxxy, g_0_y_xzzz_xxxxxz, g_0_y_xzzz_xxxxy, g_0_y_xzzz_xxxxyy, g_0_y_xzzz_xxxxyz, g_0_y_xzzz_xxxxz, g_0_y_xzzz_xxxxzz, g_0_y_xzzz_xxxyy, g_0_y_xzzz_xxxyyy, g_0_y_xzzz_xxxyyz, g_0_y_xzzz_xxxyz, g_0_y_xzzz_xxxyzz, g_0_y_xzzz_xxxzz, g_0_y_xzzz_xxxzzz, g_0_y_xzzz_xxyyy, g_0_y_xzzz_xxyyyy, g_0_y_xzzz_xxyyyz, g_0_y_xzzz_xxyyz, g_0_y_xzzz_xxyyzz, g_0_y_xzzz_xxyzz, g_0_y_xzzz_xxyzzz, g_0_y_xzzz_xxzzz, g_0_y_xzzz_xxzzzz, g_0_y_xzzz_xyyyy, g_0_y_xzzz_xyyyyy, g_0_y_xzzz_xyyyyz, g_0_y_xzzz_xyyyz, g_0_y_xzzz_xyyyzz, g_0_y_xzzz_xyyzz, g_0_y_xzzz_xyyzzz, g_0_y_xzzz_xyzzz, g_0_y_xzzz_xyzzzz, g_0_y_xzzz_xzzzz, g_0_y_xzzz_xzzzzz, g_0_y_xzzz_yyyyy, g_0_y_xzzz_yyyyz, g_0_y_xzzz_yyyzz, g_0_y_xzzz_yyzzz, g_0_y_xzzz_yzzzz, g_0_y_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzz_xxxxx[k] = -g_0_y_xzzz_xxxxx[k] * ab_x + g_0_y_xzzz_xxxxxx[k];

                g_0_y_xxzzz_xxxxy[k] = -g_0_y_xzzz_xxxxy[k] * ab_x + g_0_y_xzzz_xxxxxy[k];

                g_0_y_xxzzz_xxxxz[k] = -g_0_y_xzzz_xxxxz[k] * ab_x + g_0_y_xzzz_xxxxxz[k];

                g_0_y_xxzzz_xxxyy[k] = -g_0_y_xzzz_xxxyy[k] * ab_x + g_0_y_xzzz_xxxxyy[k];

                g_0_y_xxzzz_xxxyz[k] = -g_0_y_xzzz_xxxyz[k] * ab_x + g_0_y_xzzz_xxxxyz[k];

                g_0_y_xxzzz_xxxzz[k] = -g_0_y_xzzz_xxxzz[k] * ab_x + g_0_y_xzzz_xxxxzz[k];

                g_0_y_xxzzz_xxyyy[k] = -g_0_y_xzzz_xxyyy[k] * ab_x + g_0_y_xzzz_xxxyyy[k];

                g_0_y_xxzzz_xxyyz[k] = -g_0_y_xzzz_xxyyz[k] * ab_x + g_0_y_xzzz_xxxyyz[k];

                g_0_y_xxzzz_xxyzz[k] = -g_0_y_xzzz_xxyzz[k] * ab_x + g_0_y_xzzz_xxxyzz[k];

                g_0_y_xxzzz_xxzzz[k] = -g_0_y_xzzz_xxzzz[k] * ab_x + g_0_y_xzzz_xxxzzz[k];

                g_0_y_xxzzz_xyyyy[k] = -g_0_y_xzzz_xyyyy[k] * ab_x + g_0_y_xzzz_xxyyyy[k];

                g_0_y_xxzzz_xyyyz[k] = -g_0_y_xzzz_xyyyz[k] * ab_x + g_0_y_xzzz_xxyyyz[k];

                g_0_y_xxzzz_xyyzz[k] = -g_0_y_xzzz_xyyzz[k] * ab_x + g_0_y_xzzz_xxyyzz[k];

                g_0_y_xxzzz_xyzzz[k] = -g_0_y_xzzz_xyzzz[k] * ab_x + g_0_y_xzzz_xxyzzz[k];

                g_0_y_xxzzz_xzzzz[k] = -g_0_y_xzzz_xzzzz[k] * ab_x + g_0_y_xzzz_xxzzzz[k];

                g_0_y_xxzzz_yyyyy[k] = -g_0_y_xzzz_yyyyy[k] * ab_x + g_0_y_xzzz_xyyyyy[k];

                g_0_y_xxzzz_yyyyz[k] = -g_0_y_xzzz_yyyyz[k] * ab_x + g_0_y_xzzz_xyyyyz[k];

                g_0_y_xxzzz_yyyzz[k] = -g_0_y_xzzz_yyyzz[k] * ab_x + g_0_y_xzzz_xyyyzz[k];

                g_0_y_xxzzz_yyzzz[k] = -g_0_y_xzzz_yyzzz[k] * ab_x + g_0_y_xzzz_xyyzzz[k];

                g_0_y_xxzzz_yzzzz[k] = -g_0_y_xzzz_yzzzz[k] * ab_x + g_0_y_xzzz_xyzzzz[k];

                g_0_y_xxzzz_zzzzz[k] = -g_0_y_xzzz_zzzzz[k] * ab_x + g_0_y_xzzz_xzzzzz[k];
            }

            /// Set up 651-672 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyy_xxxxx = cbuffer.data(hh_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxy = cbuffer.data(hh_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxxz = cbuffer.data(hh_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxyy = cbuffer.data(hh_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxyz = cbuffer.data(hh_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxzz = cbuffer.data(hh_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxyyy = cbuffer.data(hh_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxyyz = cbuffer.data(hh_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxyzz = cbuffer.data(hh_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxzzz = cbuffer.data(hh_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyyyy = cbuffer.data(hh_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyyyz = cbuffer.data(hh_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyyzz = cbuffer.data(hh_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyzzz = cbuffer.data(hh_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_y_xyyyy_xzzzz = cbuffer.data(hh_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyyyy = cbuffer.data(hh_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyyyz = cbuffer.data(hh_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyyzz = cbuffer.data(hh_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyzzz = cbuffer.data(hh_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_y_xyyyy_yzzzz = cbuffer.data(hh_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_y_xyyyy_zzzzz = cbuffer.data(hh_geom_01_off + 671 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyy_xxxxx, g_0_y_xyyyy_xxxxy, g_0_y_xyyyy_xxxxz, g_0_y_xyyyy_xxxyy, g_0_y_xyyyy_xxxyz, g_0_y_xyyyy_xxxzz, g_0_y_xyyyy_xxyyy, g_0_y_xyyyy_xxyyz, g_0_y_xyyyy_xxyzz, g_0_y_xyyyy_xxzzz, g_0_y_xyyyy_xyyyy, g_0_y_xyyyy_xyyyz, g_0_y_xyyyy_xyyzz, g_0_y_xyyyy_xyzzz, g_0_y_xyyyy_xzzzz, g_0_y_xyyyy_yyyyy, g_0_y_xyyyy_yyyyz, g_0_y_xyyyy_yyyzz, g_0_y_xyyyy_yyzzz, g_0_y_xyyyy_yzzzz, g_0_y_xyyyy_zzzzz, g_0_y_yyyy_xxxxx, g_0_y_yyyy_xxxxxx, g_0_y_yyyy_xxxxxy, g_0_y_yyyy_xxxxxz, g_0_y_yyyy_xxxxy, g_0_y_yyyy_xxxxyy, g_0_y_yyyy_xxxxyz, g_0_y_yyyy_xxxxz, g_0_y_yyyy_xxxxzz, g_0_y_yyyy_xxxyy, g_0_y_yyyy_xxxyyy, g_0_y_yyyy_xxxyyz, g_0_y_yyyy_xxxyz, g_0_y_yyyy_xxxyzz, g_0_y_yyyy_xxxzz, g_0_y_yyyy_xxxzzz, g_0_y_yyyy_xxyyy, g_0_y_yyyy_xxyyyy, g_0_y_yyyy_xxyyyz, g_0_y_yyyy_xxyyz, g_0_y_yyyy_xxyyzz, g_0_y_yyyy_xxyzz, g_0_y_yyyy_xxyzzz, g_0_y_yyyy_xxzzz, g_0_y_yyyy_xxzzzz, g_0_y_yyyy_xyyyy, g_0_y_yyyy_xyyyyy, g_0_y_yyyy_xyyyyz, g_0_y_yyyy_xyyyz, g_0_y_yyyy_xyyyzz, g_0_y_yyyy_xyyzz, g_0_y_yyyy_xyyzzz, g_0_y_yyyy_xyzzz, g_0_y_yyyy_xyzzzz, g_0_y_yyyy_xzzzz, g_0_y_yyyy_xzzzzz, g_0_y_yyyy_yyyyy, g_0_y_yyyy_yyyyz, g_0_y_yyyy_yyyzz, g_0_y_yyyy_yyzzz, g_0_y_yyyy_yzzzz, g_0_y_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyy_xxxxx[k] = -g_0_y_yyyy_xxxxx[k] * ab_x + g_0_y_yyyy_xxxxxx[k];

                g_0_y_xyyyy_xxxxy[k] = -g_0_y_yyyy_xxxxy[k] * ab_x + g_0_y_yyyy_xxxxxy[k];

                g_0_y_xyyyy_xxxxz[k] = -g_0_y_yyyy_xxxxz[k] * ab_x + g_0_y_yyyy_xxxxxz[k];

                g_0_y_xyyyy_xxxyy[k] = -g_0_y_yyyy_xxxyy[k] * ab_x + g_0_y_yyyy_xxxxyy[k];

                g_0_y_xyyyy_xxxyz[k] = -g_0_y_yyyy_xxxyz[k] * ab_x + g_0_y_yyyy_xxxxyz[k];

                g_0_y_xyyyy_xxxzz[k] = -g_0_y_yyyy_xxxzz[k] * ab_x + g_0_y_yyyy_xxxxzz[k];

                g_0_y_xyyyy_xxyyy[k] = -g_0_y_yyyy_xxyyy[k] * ab_x + g_0_y_yyyy_xxxyyy[k];

                g_0_y_xyyyy_xxyyz[k] = -g_0_y_yyyy_xxyyz[k] * ab_x + g_0_y_yyyy_xxxyyz[k];

                g_0_y_xyyyy_xxyzz[k] = -g_0_y_yyyy_xxyzz[k] * ab_x + g_0_y_yyyy_xxxyzz[k];

                g_0_y_xyyyy_xxzzz[k] = -g_0_y_yyyy_xxzzz[k] * ab_x + g_0_y_yyyy_xxxzzz[k];

                g_0_y_xyyyy_xyyyy[k] = -g_0_y_yyyy_xyyyy[k] * ab_x + g_0_y_yyyy_xxyyyy[k];

                g_0_y_xyyyy_xyyyz[k] = -g_0_y_yyyy_xyyyz[k] * ab_x + g_0_y_yyyy_xxyyyz[k];

                g_0_y_xyyyy_xyyzz[k] = -g_0_y_yyyy_xyyzz[k] * ab_x + g_0_y_yyyy_xxyyzz[k];

                g_0_y_xyyyy_xyzzz[k] = -g_0_y_yyyy_xyzzz[k] * ab_x + g_0_y_yyyy_xxyzzz[k];

                g_0_y_xyyyy_xzzzz[k] = -g_0_y_yyyy_xzzzz[k] * ab_x + g_0_y_yyyy_xxzzzz[k];

                g_0_y_xyyyy_yyyyy[k] = -g_0_y_yyyy_yyyyy[k] * ab_x + g_0_y_yyyy_xyyyyy[k];

                g_0_y_xyyyy_yyyyz[k] = -g_0_y_yyyy_yyyyz[k] * ab_x + g_0_y_yyyy_xyyyyz[k];

                g_0_y_xyyyy_yyyzz[k] = -g_0_y_yyyy_yyyzz[k] * ab_x + g_0_y_yyyy_xyyyzz[k];

                g_0_y_xyyyy_yyzzz[k] = -g_0_y_yyyy_yyzzz[k] * ab_x + g_0_y_yyyy_xyyzzz[k];

                g_0_y_xyyyy_yzzzz[k] = -g_0_y_yyyy_yzzzz[k] * ab_x + g_0_y_yyyy_xyzzzz[k];

                g_0_y_xyyyy_zzzzz[k] = -g_0_y_yyyy_zzzzz[k] * ab_x + g_0_y_yyyy_xzzzzz[k];
            }

            /// Set up 672-693 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyz_xxxxx = cbuffer.data(hh_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxy = cbuffer.data(hh_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxxz = cbuffer.data(hh_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxyy = cbuffer.data(hh_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxyz = cbuffer.data(hh_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxzz = cbuffer.data(hh_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxyyy = cbuffer.data(hh_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxyyz = cbuffer.data(hh_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxyzz = cbuffer.data(hh_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxzzz = cbuffer.data(hh_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyyyy = cbuffer.data(hh_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyyyz = cbuffer.data(hh_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyyzz = cbuffer.data(hh_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyzzz = cbuffer.data(hh_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_y_xyyyz_xzzzz = cbuffer.data(hh_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyyyy = cbuffer.data(hh_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyyyz = cbuffer.data(hh_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyyzz = cbuffer.data(hh_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyzzz = cbuffer.data(hh_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_y_xyyyz_yzzzz = cbuffer.data(hh_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_y_xyyyz_zzzzz = cbuffer.data(hh_geom_01_off + 692 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyz_xxxxx, g_0_y_xyyyz_xxxxy, g_0_y_xyyyz_xxxxz, g_0_y_xyyyz_xxxyy, g_0_y_xyyyz_xxxyz, g_0_y_xyyyz_xxxzz, g_0_y_xyyyz_xxyyy, g_0_y_xyyyz_xxyyz, g_0_y_xyyyz_xxyzz, g_0_y_xyyyz_xxzzz, g_0_y_xyyyz_xyyyy, g_0_y_xyyyz_xyyyz, g_0_y_xyyyz_xyyzz, g_0_y_xyyyz_xyzzz, g_0_y_xyyyz_xzzzz, g_0_y_xyyyz_yyyyy, g_0_y_xyyyz_yyyyz, g_0_y_xyyyz_yyyzz, g_0_y_xyyyz_yyzzz, g_0_y_xyyyz_yzzzz, g_0_y_xyyyz_zzzzz, g_0_y_yyyz_xxxxx, g_0_y_yyyz_xxxxxx, g_0_y_yyyz_xxxxxy, g_0_y_yyyz_xxxxxz, g_0_y_yyyz_xxxxy, g_0_y_yyyz_xxxxyy, g_0_y_yyyz_xxxxyz, g_0_y_yyyz_xxxxz, g_0_y_yyyz_xxxxzz, g_0_y_yyyz_xxxyy, g_0_y_yyyz_xxxyyy, g_0_y_yyyz_xxxyyz, g_0_y_yyyz_xxxyz, g_0_y_yyyz_xxxyzz, g_0_y_yyyz_xxxzz, g_0_y_yyyz_xxxzzz, g_0_y_yyyz_xxyyy, g_0_y_yyyz_xxyyyy, g_0_y_yyyz_xxyyyz, g_0_y_yyyz_xxyyz, g_0_y_yyyz_xxyyzz, g_0_y_yyyz_xxyzz, g_0_y_yyyz_xxyzzz, g_0_y_yyyz_xxzzz, g_0_y_yyyz_xxzzzz, g_0_y_yyyz_xyyyy, g_0_y_yyyz_xyyyyy, g_0_y_yyyz_xyyyyz, g_0_y_yyyz_xyyyz, g_0_y_yyyz_xyyyzz, g_0_y_yyyz_xyyzz, g_0_y_yyyz_xyyzzz, g_0_y_yyyz_xyzzz, g_0_y_yyyz_xyzzzz, g_0_y_yyyz_xzzzz, g_0_y_yyyz_xzzzzz, g_0_y_yyyz_yyyyy, g_0_y_yyyz_yyyyz, g_0_y_yyyz_yyyzz, g_0_y_yyyz_yyzzz, g_0_y_yyyz_yzzzz, g_0_y_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyz_xxxxx[k] = -g_0_y_yyyz_xxxxx[k] * ab_x + g_0_y_yyyz_xxxxxx[k];

                g_0_y_xyyyz_xxxxy[k] = -g_0_y_yyyz_xxxxy[k] * ab_x + g_0_y_yyyz_xxxxxy[k];

                g_0_y_xyyyz_xxxxz[k] = -g_0_y_yyyz_xxxxz[k] * ab_x + g_0_y_yyyz_xxxxxz[k];

                g_0_y_xyyyz_xxxyy[k] = -g_0_y_yyyz_xxxyy[k] * ab_x + g_0_y_yyyz_xxxxyy[k];

                g_0_y_xyyyz_xxxyz[k] = -g_0_y_yyyz_xxxyz[k] * ab_x + g_0_y_yyyz_xxxxyz[k];

                g_0_y_xyyyz_xxxzz[k] = -g_0_y_yyyz_xxxzz[k] * ab_x + g_0_y_yyyz_xxxxzz[k];

                g_0_y_xyyyz_xxyyy[k] = -g_0_y_yyyz_xxyyy[k] * ab_x + g_0_y_yyyz_xxxyyy[k];

                g_0_y_xyyyz_xxyyz[k] = -g_0_y_yyyz_xxyyz[k] * ab_x + g_0_y_yyyz_xxxyyz[k];

                g_0_y_xyyyz_xxyzz[k] = -g_0_y_yyyz_xxyzz[k] * ab_x + g_0_y_yyyz_xxxyzz[k];

                g_0_y_xyyyz_xxzzz[k] = -g_0_y_yyyz_xxzzz[k] * ab_x + g_0_y_yyyz_xxxzzz[k];

                g_0_y_xyyyz_xyyyy[k] = -g_0_y_yyyz_xyyyy[k] * ab_x + g_0_y_yyyz_xxyyyy[k];

                g_0_y_xyyyz_xyyyz[k] = -g_0_y_yyyz_xyyyz[k] * ab_x + g_0_y_yyyz_xxyyyz[k];

                g_0_y_xyyyz_xyyzz[k] = -g_0_y_yyyz_xyyzz[k] * ab_x + g_0_y_yyyz_xxyyzz[k];

                g_0_y_xyyyz_xyzzz[k] = -g_0_y_yyyz_xyzzz[k] * ab_x + g_0_y_yyyz_xxyzzz[k];

                g_0_y_xyyyz_xzzzz[k] = -g_0_y_yyyz_xzzzz[k] * ab_x + g_0_y_yyyz_xxzzzz[k];

                g_0_y_xyyyz_yyyyy[k] = -g_0_y_yyyz_yyyyy[k] * ab_x + g_0_y_yyyz_xyyyyy[k];

                g_0_y_xyyyz_yyyyz[k] = -g_0_y_yyyz_yyyyz[k] * ab_x + g_0_y_yyyz_xyyyyz[k];

                g_0_y_xyyyz_yyyzz[k] = -g_0_y_yyyz_yyyzz[k] * ab_x + g_0_y_yyyz_xyyyzz[k];

                g_0_y_xyyyz_yyzzz[k] = -g_0_y_yyyz_yyzzz[k] * ab_x + g_0_y_yyyz_xyyzzz[k];

                g_0_y_xyyyz_yzzzz[k] = -g_0_y_yyyz_yzzzz[k] * ab_x + g_0_y_yyyz_xyzzzz[k];

                g_0_y_xyyyz_zzzzz[k] = -g_0_y_yyyz_zzzzz[k] * ab_x + g_0_y_yyyz_xzzzzz[k];
            }

            /// Set up 693-714 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzz_xxxxx = cbuffer.data(hh_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxy = cbuffer.data(hh_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxxz = cbuffer.data(hh_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxyy = cbuffer.data(hh_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxyz = cbuffer.data(hh_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxzz = cbuffer.data(hh_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxyyy = cbuffer.data(hh_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxyyz = cbuffer.data(hh_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxyzz = cbuffer.data(hh_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxzzz = cbuffer.data(hh_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyyyy = cbuffer.data(hh_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyyyz = cbuffer.data(hh_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyyzz = cbuffer.data(hh_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyzzz = cbuffer.data(hh_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_y_xyyzz_xzzzz = cbuffer.data(hh_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyyyy = cbuffer.data(hh_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyyyz = cbuffer.data(hh_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyyzz = cbuffer.data(hh_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyzzz = cbuffer.data(hh_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_y_xyyzz_yzzzz = cbuffer.data(hh_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_y_xyyzz_zzzzz = cbuffer.data(hh_geom_01_off + 713 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzz_xxxxx, g_0_y_xyyzz_xxxxy, g_0_y_xyyzz_xxxxz, g_0_y_xyyzz_xxxyy, g_0_y_xyyzz_xxxyz, g_0_y_xyyzz_xxxzz, g_0_y_xyyzz_xxyyy, g_0_y_xyyzz_xxyyz, g_0_y_xyyzz_xxyzz, g_0_y_xyyzz_xxzzz, g_0_y_xyyzz_xyyyy, g_0_y_xyyzz_xyyyz, g_0_y_xyyzz_xyyzz, g_0_y_xyyzz_xyzzz, g_0_y_xyyzz_xzzzz, g_0_y_xyyzz_yyyyy, g_0_y_xyyzz_yyyyz, g_0_y_xyyzz_yyyzz, g_0_y_xyyzz_yyzzz, g_0_y_xyyzz_yzzzz, g_0_y_xyyzz_zzzzz, g_0_y_yyzz_xxxxx, g_0_y_yyzz_xxxxxx, g_0_y_yyzz_xxxxxy, g_0_y_yyzz_xxxxxz, g_0_y_yyzz_xxxxy, g_0_y_yyzz_xxxxyy, g_0_y_yyzz_xxxxyz, g_0_y_yyzz_xxxxz, g_0_y_yyzz_xxxxzz, g_0_y_yyzz_xxxyy, g_0_y_yyzz_xxxyyy, g_0_y_yyzz_xxxyyz, g_0_y_yyzz_xxxyz, g_0_y_yyzz_xxxyzz, g_0_y_yyzz_xxxzz, g_0_y_yyzz_xxxzzz, g_0_y_yyzz_xxyyy, g_0_y_yyzz_xxyyyy, g_0_y_yyzz_xxyyyz, g_0_y_yyzz_xxyyz, g_0_y_yyzz_xxyyzz, g_0_y_yyzz_xxyzz, g_0_y_yyzz_xxyzzz, g_0_y_yyzz_xxzzz, g_0_y_yyzz_xxzzzz, g_0_y_yyzz_xyyyy, g_0_y_yyzz_xyyyyy, g_0_y_yyzz_xyyyyz, g_0_y_yyzz_xyyyz, g_0_y_yyzz_xyyyzz, g_0_y_yyzz_xyyzz, g_0_y_yyzz_xyyzzz, g_0_y_yyzz_xyzzz, g_0_y_yyzz_xyzzzz, g_0_y_yyzz_xzzzz, g_0_y_yyzz_xzzzzz, g_0_y_yyzz_yyyyy, g_0_y_yyzz_yyyyz, g_0_y_yyzz_yyyzz, g_0_y_yyzz_yyzzz, g_0_y_yyzz_yzzzz, g_0_y_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzz_xxxxx[k] = -g_0_y_yyzz_xxxxx[k] * ab_x + g_0_y_yyzz_xxxxxx[k];

                g_0_y_xyyzz_xxxxy[k] = -g_0_y_yyzz_xxxxy[k] * ab_x + g_0_y_yyzz_xxxxxy[k];

                g_0_y_xyyzz_xxxxz[k] = -g_0_y_yyzz_xxxxz[k] * ab_x + g_0_y_yyzz_xxxxxz[k];

                g_0_y_xyyzz_xxxyy[k] = -g_0_y_yyzz_xxxyy[k] * ab_x + g_0_y_yyzz_xxxxyy[k];

                g_0_y_xyyzz_xxxyz[k] = -g_0_y_yyzz_xxxyz[k] * ab_x + g_0_y_yyzz_xxxxyz[k];

                g_0_y_xyyzz_xxxzz[k] = -g_0_y_yyzz_xxxzz[k] * ab_x + g_0_y_yyzz_xxxxzz[k];

                g_0_y_xyyzz_xxyyy[k] = -g_0_y_yyzz_xxyyy[k] * ab_x + g_0_y_yyzz_xxxyyy[k];

                g_0_y_xyyzz_xxyyz[k] = -g_0_y_yyzz_xxyyz[k] * ab_x + g_0_y_yyzz_xxxyyz[k];

                g_0_y_xyyzz_xxyzz[k] = -g_0_y_yyzz_xxyzz[k] * ab_x + g_0_y_yyzz_xxxyzz[k];

                g_0_y_xyyzz_xxzzz[k] = -g_0_y_yyzz_xxzzz[k] * ab_x + g_0_y_yyzz_xxxzzz[k];

                g_0_y_xyyzz_xyyyy[k] = -g_0_y_yyzz_xyyyy[k] * ab_x + g_0_y_yyzz_xxyyyy[k];

                g_0_y_xyyzz_xyyyz[k] = -g_0_y_yyzz_xyyyz[k] * ab_x + g_0_y_yyzz_xxyyyz[k];

                g_0_y_xyyzz_xyyzz[k] = -g_0_y_yyzz_xyyzz[k] * ab_x + g_0_y_yyzz_xxyyzz[k];

                g_0_y_xyyzz_xyzzz[k] = -g_0_y_yyzz_xyzzz[k] * ab_x + g_0_y_yyzz_xxyzzz[k];

                g_0_y_xyyzz_xzzzz[k] = -g_0_y_yyzz_xzzzz[k] * ab_x + g_0_y_yyzz_xxzzzz[k];

                g_0_y_xyyzz_yyyyy[k] = -g_0_y_yyzz_yyyyy[k] * ab_x + g_0_y_yyzz_xyyyyy[k];

                g_0_y_xyyzz_yyyyz[k] = -g_0_y_yyzz_yyyyz[k] * ab_x + g_0_y_yyzz_xyyyyz[k];

                g_0_y_xyyzz_yyyzz[k] = -g_0_y_yyzz_yyyzz[k] * ab_x + g_0_y_yyzz_xyyyzz[k];

                g_0_y_xyyzz_yyzzz[k] = -g_0_y_yyzz_yyzzz[k] * ab_x + g_0_y_yyzz_xyyzzz[k];

                g_0_y_xyyzz_yzzzz[k] = -g_0_y_yyzz_yzzzz[k] * ab_x + g_0_y_yyzz_xyzzzz[k];

                g_0_y_xyyzz_zzzzz[k] = -g_0_y_yyzz_zzzzz[k] * ab_x + g_0_y_yyzz_xzzzzz[k];
            }

            /// Set up 714-735 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzz_xxxxx = cbuffer.data(hh_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxy = cbuffer.data(hh_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxxz = cbuffer.data(hh_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxyy = cbuffer.data(hh_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxyz = cbuffer.data(hh_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxzz = cbuffer.data(hh_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxyyy = cbuffer.data(hh_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxyyz = cbuffer.data(hh_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxyzz = cbuffer.data(hh_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxzzz = cbuffer.data(hh_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyyyy = cbuffer.data(hh_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyyyz = cbuffer.data(hh_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyyzz = cbuffer.data(hh_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyzzz = cbuffer.data(hh_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_y_xyzzz_xzzzz = cbuffer.data(hh_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyyyy = cbuffer.data(hh_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyyyz = cbuffer.data(hh_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyyzz = cbuffer.data(hh_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyzzz = cbuffer.data(hh_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_y_xyzzz_yzzzz = cbuffer.data(hh_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_y_xyzzz_zzzzz = cbuffer.data(hh_geom_01_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzz_xxxxx, g_0_y_xyzzz_xxxxy, g_0_y_xyzzz_xxxxz, g_0_y_xyzzz_xxxyy, g_0_y_xyzzz_xxxyz, g_0_y_xyzzz_xxxzz, g_0_y_xyzzz_xxyyy, g_0_y_xyzzz_xxyyz, g_0_y_xyzzz_xxyzz, g_0_y_xyzzz_xxzzz, g_0_y_xyzzz_xyyyy, g_0_y_xyzzz_xyyyz, g_0_y_xyzzz_xyyzz, g_0_y_xyzzz_xyzzz, g_0_y_xyzzz_xzzzz, g_0_y_xyzzz_yyyyy, g_0_y_xyzzz_yyyyz, g_0_y_xyzzz_yyyzz, g_0_y_xyzzz_yyzzz, g_0_y_xyzzz_yzzzz, g_0_y_xyzzz_zzzzz, g_0_y_yzzz_xxxxx, g_0_y_yzzz_xxxxxx, g_0_y_yzzz_xxxxxy, g_0_y_yzzz_xxxxxz, g_0_y_yzzz_xxxxy, g_0_y_yzzz_xxxxyy, g_0_y_yzzz_xxxxyz, g_0_y_yzzz_xxxxz, g_0_y_yzzz_xxxxzz, g_0_y_yzzz_xxxyy, g_0_y_yzzz_xxxyyy, g_0_y_yzzz_xxxyyz, g_0_y_yzzz_xxxyz, g_0_y_yzzz_xxxyzz, g_0_y_yzzz_xxxzz, g_0_y_yzzz_xxxzzz, g_0_y_yzzz_xxyyy, g_0_y_yzzz_xxyyyy, g_0_y_yzzz_xxyyyz, g_0_y_yzzz_xxyyz, g_0_y_yzzz_xxyyzz, g_0_y_yzzz_xxyzz, g_0_y_yzzz_xxyzzz, g_0_y_yzzz_xxzzz, g_0_y_yzzz_xxzzzz, g_0_y_yzzz_xyyyy, g_0_y_yzzz_xyyyyy, g_0_y_yzzz_xyyyyz, g_0_y_yzzz_xyyyz, g_0_y_yzzz_xyyyzz, g_0_y_yzzz_xyyzz, g_0_y_yzzz_xyyzzz, g_0_y_yzzz_xyzzz, g_0_y_yzzz_xyzzzz, g_0_y_yzzz_xzzzz, g_0_y_yzzz_xzzzzz, g_0_y_yzzz_yyyyy, g_0_y_yzzz_yyyyz, g_0_y_yzzz_yyyzz, g_0_y_yzzz_yyzzz, g_0_y_yzzz_yzzzz, g_0_y_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzz_xxxxx[k] = -g_0_y_yzzz_xxxxx[k] * ab_x + g_0_y_yzzz_xxxxxx[k];

                g_0_y_xyzzz_xxxxy[k] = -g_0_y_yzzz_xxxxy[k] * ab_x + g_0_y_yzzz_xxxxxy[k];

                g_0_y_xyzzz_xxxxz[k] = -g_0_y_yzzz_xxxxz[k] * ab_x + g_0_y_yzzz_xxxxxz[k];

                g_0_y_xyzzz_xxxyy[k] = -g_0_y_yzzz_xxxyy[k] * ab_x + g_0_y_yzzz_xxxxyy[k];

                g_0_y_xyzzz_xxxyz[k] = -g_0_y_yzzz_xxxyz[k] * ab_x + g_0_y_yzzz_xxxxyz[k];

                g_0_y_xyzzz_xxxzz[k] = -g_0_y_yzzz_xxxzz[k] * ab_x + g_0_y_yzzz_xxxxzz[k];

                g_0_y_xyzzz_xxyyy[k] = -g_0_y_yzzz_xxyyy[k] * ab_x + g_0_y_yzzz_xxxyyy[k];

                g_0_y_xyzzz_xxyyz[k] = -g_0_y_yzzz_xxyyz[k] * ab_x + g_0_y_yzzz_xxxyyz[k];

                g_0_y_xyzzz_xxyzz[k] = -g_0_y_yzzz_xxyzz[k] * ab_x + g_0_y_yzzz_xxxyzz[k];

                g_0_y_xyzzz_xxzzz[k] = -g_0_y_yzzz_xxzzz[k] * ab_x + g_0_y_yzzz_xxxzzz[k];

                g_0_y_xyzzz_xyyyy[k] = -g_0_y_yzzz_xyyyy[k] * ab_x + g_0_y_yzzz_xxyyyy[k];

                g_0_y_xyzzz_xyyyz[k] = -g_0_y_yzzz_xyyyz[k] * ab_x + g_0_y_yzzz_xxyyyz[k];

                g_0_y_xyzzz_xyyzz[k] = -g_0_y_yzzz_xyyzz[k] * ab_x + g_0_y_yzzz_xxyyzz[k];

                g_0_y_xyzzz_xyzzz[k] = -g_0_y_yzzz_xyzzz[k] * ab_x + g_0_y_yzzz_xxyzzz[k];

                g_0_y_xyzzz_xzzzz[k] = -g_0_y_yzzz_xzzzz[k] * ab_x + g_0_y_yzzz_xxzzzz[k];

                g_0_y_xyzzz_yyyyy[k] = -g_0_y_yzzz_yyyyy[k] * ab_x + g_0_y_yzzz_xyyyyy[k];

                g_0_y_xyzzz_yyyyz[k] = -g_0_y_yzzz_yyyyz[k] * ab_x + g_0_y_yzzz_xyyyyz[k];

                g_0_y_xyzzz_yyyzz[k] = -g_0_y_yzzz_yyyzz[k] * ab_x + g_0_y_yzzz_xyyyzz[k];

                g_0_y_xyzzz_yyzzz[k] = -g_0_y_yzzz_yyzzz[k] * ab_x + g_0_y_yzzz_xyyzzz[k];

                g_0_y_xyzzz_yzzzz[k] = -g_0_y_yzzz_yzzzz[k] * ab_x + g_0_y_yzzz_xyzzzz[k];

                g_0_y_xyzzz_zzzzz[k] = -g_0_y_yzzz_zzzzz[k] * ab_x + g_0_y_yzzz_xzzzzz[k];
            }

            /// Set up 735-756 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzz_xxxxx = cbuffer.data(hh_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxy = cbuffer.data(hh_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxxz = cbuffer.data(hh_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxyy = cbuffer.data(hh_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxyz = cbuffer.data(hh_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxzz = cbuffer.data(hh_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxyyy = cbuffer.data(hh_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxyyz = cbuffer.data(hh_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxyzz = cbuffer.data(hh_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxzzz = cbuffer.data(hh_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyyyy = cbuffer.data(hh_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyyyz = cbuffer.data(hh_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyyzz = cbuffer.data(hh_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyzzz = cbuffer.data(hh_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_y_xzzzz_xzzzz = cbuffer.data(hh_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyyyy = cbuffer.data(hh_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyyyz = cbuffer.data(hh_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyyzz = cbuffer.data(hh_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyzzz = cbuffer.data(hh_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_y_xzzzz_yzzzz = cbuffer.data(hh_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_y_xzzzz_zzzzz = cbuffer.data(hh_geom_01_off + 755 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzz_xxxxx, g_0_y_xzzzz_xxxxy, g_0_y_xzzzz_xxxxz, g_0_y_xzzzz_xxxyy, g_0_y_xzzzz_xxxyz, g_0_y_xzzzz_xxxzz, g_0_y_xzzzz_xxyyy, g_0_y_xzzzz_xxyyz, g_0_y_xzzzz_xxyzz, g_0_y_xzzzz_xxzzz, g_0_y_xzzzz_xyyyy, g_0_y_xzzzz_xyyyz, g_0_y_xzzzz_xyyzz, g_0_y_xzzzz_xyzzz, g_0_y_xzzzz_xzzzz, g_0_y_xzzzz_yyyyy, g_0_y_xzzzz_yyyyz, g_0_y_xzzzz_yyyzz, g_0_y_xzzzz_yyzzz, g_0_y_xzzzz_yzzzz, g_0_y_xzzzz_zzzzz, g_0_y_zzzz_xxxxx, g_0_y_zzzz_xxxxxx, g_0_y_zzzz_xxxxxy, g_0_y_zzzz_xxxxxz, g_0_y_zzzz_xxxxy, g_0_y_zzzz_xxxxyy, g_0_y_zzzz_xxxxyz, g_0_y_zzzz_xxxxz, g_0_y_zzzz_xxxxzz, g_0_y_zzzz_xxxyy, g_0_y_zzzz_xxxyyy, g_0_y_zzzz_xxxyyz, g_0_y_zzzz_xxxyz, g_0_y_zzzz_xxxyzz, g_0_y_zzzz_xxxzz, g_0_y_zzzz_xxxzzz, g_0_y_zzzz_xxyyy, g_0_y_zzzz_xxyyyy, g_0_y_zzzz_xxyyyz, g_0_y_zzzz_xxyyz, g_0_y_zzzz_xxyyzz, g_0_y_zzzz_xxyzz, g_0_y_zzzz_xxyzzz, g_0_y_zzzz_xxzzz, g_0_y_zzzz_xxzzzz, g_0_y_zzzz_xyyyy, g_0_y_zzzz_xyyyyy, g_0_y_zzzz_xyyyyz, g_0_y_zzzz_xyyyz, g_0_y_zzzz_xyyyzz, g_0_y_zzzz_xyyzz, g_0_y_zzzz_xyyzzz, g_0_y_zzzz_xyzzz, g_0_y_zzzz_xyzzzz, g_0_y_zzzz_xzzzz, g_0_y_zzzz_xzzzzz, g_0_y_zzzz_yyyyy, g_0_y_zzzz_yyyyz, g_0_y_zzzz_yyyzz, g_0_y_zzzz_yyzzz, g_0_y_zzzz_yzzzz, g_0_y_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzz_xxxxx[k] = -g_0_y_zzzz_xxxxx[k] * ab_x + g_0_y_zzzz_xxxxxx[k];

                g_0_y_xzzzz_xxxxy[k] = -g_0_y_zzzz_xxxxy[k] * ab_x + g_0_y_zzzz_xxxxxy[k];

                g_0_y_xzzzz_xxxxz[k] = -g_0_y_zzzz_xxxxz[k] * ab_x + g_0_y_zzzz_xxxxxz[k];

                g_0_y_xzzzz_xxxyy[k] = -g_0_y_zzzz_xxxyy[k] * ab_x + g_0_y_zzzz_xxxxyy[k];

                g_0_y_xzzzz_xxxyz[k] = -g_0_y_zzzz_xxxyz[k] * ab_x + g_0_y_zzzz_xxxxyz[k];

                g_0_y_xzzzz_xxxzz[k] = -g_0_y_zzzz_xxxzz[k] * ab_x + g_0_y_zzzz_xxxxzz[k];

                g_0_y_xzzzz_xxyyy[k] = -g_0_y_zzzz_xxyyy[k] * ab_x + g_0_y_zzzz_xxxyyy[k];

                g_0_y_xzzzz_xxyyz[k] = -g_0_y_zzzz_xxyyz[k] * ab_x + g_0_y_zzzz_xxxyyz[k];

                g_0_y_xzzzz_xxyzz[k] = -g_0_y_zzzz_xxyzz[k] * ab_x + g_0_y_zzzz_xxxyzz[k];

                g_0_y_xzzzz_xxzzz[k] = -g_0_y_zzzz_xxzzz[k] * ab_x + g_0_y_zzzz_xxxzzz[k];

                g_0_y_xzzzz_xyyyy[k] = -g_0_y_zzzz_xyyyy[k] * ab_x + g_0_y_zzzz_xxyyyy[k];

                g_0_y_xzzzz_xyyyz[k] = -g_0_y_zzzz_xyyyz[k] * ab_x + g_0_y_zzzz_xxyyyz[k];

                g_0_y_xzzzz_xyyzz[k] = -g_0_y_zzzz_xyyzz[k] * ab_x + g_0_y_zzzz_xxyyzz[k];

                g_0_y_xzzzz_xyzzz[k] = -g_0_y_zzzz_xyzzz[k] * ab_x + g_0_y_zzzz_xxyzzz[k];

                g_0_y_xzzzz_xzzzz[k] = -g_0_y_zzzz_xzzzz[k] * ab_x + g_0_y_zzzz_xxzzzz[k];

                g_0_y_xzzzz_yyyyy[k] = -g_0_y_zzzz_yyyyy[k] * ab_x + g_0_y_zzzz_xyyyyy[k];

                g_0_y_xzzzz_yyyyz[k] = -g_0_y_zzzz_yyyyz[k] * ab_x + g_0_y_zzzz_xyyyyz[k];

                g_0_y_xzzzz_yyyzz[k] = -g_0_y_zzzz_yyyzz[k] * ab_x + g_0_y_zzzz_xyyyzz[k];

                g_0_y_xzzzz_yyzzz[k] = -g_0_y_zzzz_yyzzz[k] * ab_x + g_0_y_zzzz_xyyzzz[k];

                g_0_y_xzzzz_yzzzz[k] = -g_0_y_zzzz_yzzzz[k] * ab_x + g_0_y_zzzz_xyzzzz[k];

                g_0_y_xzzzz_zzzzz[k] = -g_0_y_zzzz_zzzzz[k] * ab_x + g_0_y_zzzz_xzzzzz[k];
            }

            /// Set up 756-777 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyy_xxxxx = cbuffer.data(hh_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxy = cbuffer.data(hh_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxxz = cbuffer.data(hh_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxyy = cbuffer.data(hh_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxyz = cbuffer.data(hh_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxzz = cbuffer.data(hh_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxyyy = cbuffer.data(hh_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxyyz = cbuffer.data(hh_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxyzz = cbuffer.data(hh_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxzzz = cbuffer.data(hh_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyyyy = cbuffer.data(hh_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyyyz = cbuffer.data(hh_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyyzz = cbuffer.data(hh_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyzzz = cbuffer.data(hh_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_y_yyyyy_xzzzz = cbuffer.data(hh_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyyyy = cbuffer.data(hh_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyyyz = cbuffer.data(hh_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyyzz = cbuffer.data(hh_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyzzz = cbuffer.data(hh_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_y_yyyyy_yzzzz = cbuffer.data(hh_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_y_yyyyy_zzzzz = cbuffer.data(hh_geom_01_off + 776 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_xxxxx, g_0_y_yyyy_xxxxxy, g_0_y_yyyy_xxxxy, g_0_y_yyyy_xxxxyy, g_0_y_yyyy_xxxxyz, g_0_y_yyyy_xxxxz, g_0_y_yyyy_xxxyy, g_0_y_yyyy_xxxyyy, g_0_y_yyyy_xxxyyz, g_0_y_yyyy_xxxyz, g_0_y_yyyy_xxxyzz, g_0_y_yyyy_xxxzz, g_0_y_yyyy_xxyyy, g_0_y_yyyy_xxyyyy, g_0_y_yyyy_xxyyyz, g_0_y_yyyy_xxyyz, g_0_y_yyyy_xxyyzz, g_0_y_yyyy_xxyzz, g_0_y_yyyy_xxyzzz, g_0_y_yyyy_xxzzz, g_0_y_yyyy_xyyyy, g_0_y_yyyy_xyyyyy, g_0_y_yyyy_xyyyyz, g_0_y_yyyy_xyyyz, g_0_y_yyyy_xyyyzz, g_0_y_yyyy_xyyzz, g_0_y_yyyy_xyyzzz, g_0_y_yyyy_xyzzz, g_0_y_yyyy_xyzzzz, g_0_y_yyyy_xzzzz, g_0_y_yyyy_yyyyy, g_0_y_yyyy_yyyyyy, g_0_y_yyyy_yyyyyz, g_0_y_yyyy_yyyyz, g_0_y_yyyy_yyyyzz, g_0_y_yyyy_yyyzz, g_0_y_yyyy_yyyzzz, g_0_y_yyyy_yyzzz, g_0_y_yyyy_yyzzzz, g_0_y_yyyy_yzzzz, g_0_y_yyyy_yzzzzz, g_0_y_yyyy_zzzzz, g_0_y_yyyyy_xxxxx, g_0_y_yyyyy_xxxxy, g_0_y_yyyyy_xxxxz, g_0_y_yyyyy_xxxyy, g_0_y_yyyyy_xxxyz, g_0_y_yyyyy_xxxzz, g_0_y_yyyyy_xxyyy, g_0_y_yyyyy_xxyyz, g_0_y_yyyyy_xxyzz, g_0_y_yyyyy_xxzzz, g_0_y_yyyyy_xyyyy, g_0_y_yyyyy_xyyyz, g_0_y_yyyyy_xyyzz, g_0_y_yyyyy_xyzzz, g_0_y_yyyyy_xzzzz, g_0_y_yyyyy_yyyyy, g_0_y_yyyyy_yyyyz, g_0_y_yyyyy_yyyzz, g_0_y_yyyyy_yyzzz, g_0_y_yyyyy_yzzzz, g_0_y_yyyyy_zzzzz, g_yyyy_xxxxx, g_yyyy_xxxxy, g_yyyy_xxxxz, g_yyyy_xxxyy, g_yyyy_xxxyz, g_yyyy_xxxzz, g_yyyy_xxyyy, g_yyyy_xxyyz, g_yyyy_xxyzz, g_yyyy_xxzzz, g_yyyy_xyyyy, g_yyyy_xyyyz, g_yyyy_xyyzz, g_yyyy_xyzzz, g_yyyy_xzzzz, g_yyyy_yyyyy, g_yyyy_yyyyz, g_yyyy_yyyzz, g_yyyy_yyzzz, g_yyyy_yzzzz, g_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyy_xxxxx[k] = g_yyyy_xxxxx[k] - g_0_y_yyyy_xxxxx[k] * ab_y + g_0_y_yyyy_xxxxxy[k];

                g_0_y_yyyyy_xxxxy[k] = g_yyyy_xxxxy[k] - g_0_y_yyyy_xxxxy[k] * ab_y + g_0_y_yyyy_xxxxyy[k];

                g_0_y_yyyyy_xxxxz[k] = g_yyyy_xxxxz[k] - g_0_y_yyyy_xxxxz[k] * ab_y + g_0_y_yyyy_xxxxyz[k];

                g_0_y_yyyyy_xxxyy[k] = g_yyyy_xxxyy[k] - g_0_y_yyyy_xxxyy[k] * ab_y + g_0_y_yyyy_xxxyyy[k];

                g_0_y_yyyyy_xxxyz[k] = g_yyyy_xxxyz[k] - g_0_y_yyyy_xxxyz[k] * ab_y + g_0_y_yyyy_xxxyyz[k];

                g_0_y_yyyyy_xxxzz[k] = g_yyyy_xxxzz[k] - g_0_y_yyyy_xxxzz[k] * ab_y + g_0_y_yyyy_xxxyzz[k];

                g_0_y_yyyyy_xxyyy[k] = g_yyyy_xxyyy[k] - g_0_y_yyyy_xxyyy[k] * ab_y + g_0_y_yyyy_xxyyyy[k];

                g_0_y_yyyyy_xxyyz[k] = g_yyyy_xxyyz[k] - g_0_y_yyyy_xxyyz[k] * ab_y + g_0_y_yyyy_xxyyyz[k];

                g_0_y_yyyyy_xxyzz[k] = g_yyyy_xxyzz[k] - g_0_y_yyyy_xxyzz[k] * ab_y + g_0_y_yyyy_xxyyzz[k];

                g_0_y_yyyyy_xxzzz[k] = g_yyyy_xxzzz[k] - g_0_y_yyyy_xxzzz[k] * ab_y + g_0_y_yyyy_xxyzzz[k];

                g_0_y_yyyyy_xyyyy[k] = g_yyyy_xyyyy[k] - g_0_y_yyyy_xyyyy[k] * ab_y + g_0_y_yyyy_xyyyyy[k];

                g_0_y_yyyyy_xyyyz[k] = g_yyyy_xyyyz[k] - g_0_y_yyyy_xyyyz[k] * ab_y + g_0_y_yyyy_xyyyyz[k];

                g_0_y_yyyyy_xyyzz[k] = g_yyyy_xyyzz[k] - g_0_y_yyyy_xyyzz[k] * ab_y + g_0_y_yyyy_xyyyzz[k];

                g_0_y_yyyyy_xyzzz[k] = g_yyyy_xyzzz[k] - g_0_y_yyyy_xyzzz[k] * ab_y + g_0_y_yyyy_xyyzzz[k];

                g_0_y_yyyyy_xzzzz[k] = g_yyyy_xzzzz[k] - g_0_y_yyyy_xzzzz[k] * ab_y + g_0_y_yyyy_xyzzzz[k];

                g_0_y_yyyyy_yyyyy[k] = g_yyyy_yyyyy[k] - g_0_y_yyyy_yyyyy[k] * ab_y + g_0_y_yyyy_yyyyyy[k];

                g_0_y_yyyyy_yyyyz[k] = g_yyyy_yyyyz[k] - g_0_y_yyyy_yyyyz[k] * ab_y + g_0_y_yyyy_yyyyyz[k];

                g_0_y_yyyyy_yyyzz[k] = g_yyyy_yyyzz[k] - g_0_y_yyyy_yyyzz[k] * ab_y + g_0_y_yyyy_yyyyzz[k];

                g_0_y_yyyyy_yyzzz[k] = g_yyyy_yyzzz[k] - g_0_y_yyyy_yyzzz[k] * ab_y + g_0_y_yyyy_yyyzzz[k];

                g_0_y_yyyyy_yzzzz[k] = g_yyyy_yzzzz[k] - g_0_y_yyyy_yzzzz[k] * ab_y + g_0_y_yyyy_yyzzzz[k];

                g_0_y_yyyyy_zzzzz[k] = g_yyyy_zzzzz[k] - g_0_y_yyyy_zzzzz[k] * ab_y + g_0_y_yyyy_yzzzzz[k];
            }

            /// Set up 777-798 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyz_xxxxx = cbuffer.data(hh_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxy = cbuffer.data(hh_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxxz = cbuffer.data(hh_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxyy = cbuffer.data(hh_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxyz = cbuffer.data(hh_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxzz = cbuffer.data(hh_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxyyy = cbuffer.data(hh_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxyyz = cbuffer.data(hh_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxyzz = cbuffer.data(hh_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxzzz = cbuffer.data(hh_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyyyy = cbuffer.data(hh_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyyyz = cbuffer.data(hh_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyyzz = cbuffer.data(hh_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyzzz = cbuffer.data(hh_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_y_yyyyz_xzzzz = cbuffer.data(hh_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyyyy = cbuffer.data(hh_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyyyz = cbuffer.data(hh_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyyzz = cbuffer.data(hh_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyzzz = cbuffer.data(hh_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_y_yyyyz_yzzzz = cbuffer.data(hh_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_y_yyyyz_zzzzz = cbuffer.data(hh_geom_01_off + 797 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyy_xxxxx, g_0_y_yyyy_xxxxxz, g_0_y_yyyy_xxxxy, g_0_y_yyyy_xxxxyz, g_0_y_yyyy_xxxxz, g_0_y_yyyy_xxxxzz, g_0_y_yyyy_xxxyy, g_0_y_yyyy_xxxyyz, g_0_y_yyyy_xxxyz, g_0_y_yyyy_xxxyzz, g_0_y_yyyy_xxxzz, g_0_y_yyyy_xxxzzz, g_0_y_yyyy_xxyyy, g_0_y_yyyy_xxyyyz, g_0_y_yyyy_xxyyz, g_0_y_yyyy_xxyyzz, g_0_y_yyyy_xxyzz, g_0_y_yyyy_xxyzzz, g_0_y_yyyy_xxzzz, g_0_y_yyyy_xxzzzz, g_0_y_yyyy_xyyyy, g_0_y_yyyy_xyyyyz, g_0_y_yyyy_xyyyz, g_0_y_yyyy_xyyyzz, g_0_y_yyyy_xyyzz, g_0_y_yyyy_xyyzzz, g_0_y_yyyy_xyzzz, g_0_y_yyyy_xyzzzz, g_0_y_yyyy_xzzzz, g_0_y_yyyy_xzzzzz, g_0_y_yyyy_yyyyy, g_0_y_yyyy_yyyyyz, g_0_y_yyyy_yyyyz, g_0_y_yyyy_yyyyzz, g_0_y_yyyy_yyyzz, g_0_y_yyyy_yyyzzz, g_0_y_yyyy_yyzzz, g_0_y_yyyy_yyzzzz, g_0_y_yyyy_yzzzz, g_0_y_yyyy_yzzzzz, g_0_y_yyyy_zzzzz, g_0_y_yyyy_zzzzzz, g_0_y_yyyyz_xxxxx, g_0_y_yyyyz_xxxxy, g_0_y_yyyyz_xxxxz, g_0_y_yyyyz_xxxyy, g_0_y_yyyyz_xxxyz, g_0_y_yyyyz_xxxzz, g_0_y_yyyyz_xxyyy, g_0_y_yyyyz_xxyyz, g_0_y_yyyyz_xxyzz, g_0_y_yyyyz_xxzzz, g_0_y_yyyyz_xyyyy, g_0_y_yyyyz_xyyyz, g_0_y_yyyyz_xyyzz, g_0_y_yyyyz_xyzzz, g_0_y_yyyyz_xzzzz, g_0_y_yyyyz_yyyyy, g_0_y_yyyyz_yyyyz, g_0_y_yyyyz_yyyzz, g_0_y_yyyyz_yyzzz, g_0_y_yyyyz_yzzzz, g_0_y_yyyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyz_xxxxx[k] = -g_0_y_yyyy_xxxxx[k] * ab_z + g_0_y_yyyy_xxxxxz[k];

                g_0_y_yyyyz_xxxxy[k] = -g_0_y_yyyy_xxxxy[k] * ab_z + g_0_y_yyyy_xxxxyz[k];

                g_0_y_yyyyz_xxxxz[k] = -g_0_y_yyyy_xxxxz[k] * ab_z + g_0_y_yyyy_xxxxzz[k];

                g_0_y_yyyyz_xxxyy[k] = -g_0_y_yyyy_xxxyy[k] * ab_z + g_0_y_yyyy_xxxyyz[k];

                g_0_y_yyyyz_xxxyz[k] = -g_0_y_yyyy_xxxyz[k] * ab_z + g_0_y_yyyy_xxxyzz[k];

                g_0_y_yyyyz_xxxzz[k] = -g_0_y_yyyy_xxxzz[k] * ab_z + g_0_y_yyyy_xxxzzz[k];

                g_0_y_yyyyz_xxyyy[k] = -g_0_y_yyyy_xxyyy[k] * ab_z + g_0_y_yyyy_xxyyyz[k];

                g_0_y_yyyyz_xxyyz[k] = -g_0_y_yyyy_xxyyz[k] * ab_z + g_0_y_yyyy_xxyyzz[k];

                g_0_y_yyyyz_xxyzz[k] = -g_0_y_yyyy_xxyzz[k] * ab_z + g_0_y_yyyy_xxyzzz[k];

                g_0_y_yyyyz_xxzzz[k] = -g_0_y_yyyy_xxzzz[k] * ab_z + g_0_y_yyyy_xxzzzz[k];

                g_0_y_yyyyz_xyyyy[k] = -g_0_y_yyyy_xyyyy[k] * ab_z + g_0_y_yyyy_xyyyyz[k];

                g_0_y_yyyyz_xyyyz[k] = -g_0_y_yyyy_xyyyz[k] * ab_z + g_0_y_yyyy_xyyyzz[k];

                g_0_y_yyyyz_xyyzz[k] = -g_0_y_yyyy_xyyzz[k] * ab_z + g_0_y_yyyy_xyyzzz[k];

                g_0_y_yyyyz_xyzzz[k] = -g_0_y_yyyy_xyzzz[k] * ab_z + g_0_y_yyyy_xyzzzz[k];

                g_0_y_yyyyz_xzzzz[k] = -g_0_y_yyyy_xzzzz[k] * ab_z + g_0_y_yyyy_xzzzzz[k];

                g_0_y_yyyyz_yyyyy[k] = -g_0_y_yyyy_yyyyy[k] * ab_z + g_0_y_yyyy_yyyyyz[k];

                g_0_y_yyyyz_yyyyz[k] = -g_0_y_yyyy_yyyyz[k] * ab_z + g_0_y_yyyy_yyyyzz[k];

                g_0_y_yyyyz_yyyzz[k] = -g_0_y_yyyy_yyyzz[k] * ab_z + g_0_y_yyyy_yyyzzz[k];

                g_0_y_yyyyz_yyzzz[k] = -g_0_y_yyyy_yyzzz[k] * ab_z + g_0_y_yyyy_yyzzzz[k];

                g_0_y_yyyyz_yzzzz[k] = -g_0_y_yyyy_yzzzz[k] * ab_z + g_0_y_yyyy_yzzzzz[k];

                g_0_y_yyyyz_zzzzz[k] = -g_0_y_yyyy_zzzzz[k] * ab_z + g_0_y_yyyy_zzzzzz[k];
            }

            /// Set up 798-819 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzz_xxxxx = cbuffer.data(hh_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxy = cbuffer.data(hh_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxxz = cbuffer.data(hh_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxyy = cbuffer.data(hh_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxyz = cbuffer.data(hh_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxzz = cbuffer.data(hh_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxyyy = cbuffer.data(hh_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxyyz = cbuffer.data(hh_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxyzz = cbuffer.data(hh_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxzzz = cbuffer.data(hh_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyyyy = cbuffer.data(hh_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyyyz = cbuffer.data(hh_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyyzz = cbuffer.data(hh_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyzzz = cbuffer.data(hh_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_y_yyyzz_xzzzz = cbuffer.data(hh_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyyyy = cbuffer.data(hh_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyyyz = cbuffer.data(hh_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyyzz = cbuffer.data(hh_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyzzz = cbuffer.data(hh_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_y_yyyzz_yzzzz = cbuffer.data(hh_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_y_yyyzz_zzzzz = cbuffer.data(hh_geom_01_off + 818 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyz_xxxxx, g_0_y_yyyz_xxxxxz, g_0_y_yyyz_xxxxy, g_0_y_yyyz_xxxxyz, g_0_y_yyyz_xxxxz, g_0_y_yyyz_xxxxzz, g_0_y_yyyz_xxxyy, g_0_y_yyyz_xxxyyz, g_0_y_yyyz_xxxyz, g_0_y_yyyz_xxxyzz, g_0_y_yyyz_xxxzz, g_0_y_yyyz_xxxzzz, g_0_y_yyyz_xxyyy, g_0_y_yyyz_xxyyyz, g_0_y_yyyz_xxyyz, g_0_y_yyyz_xxyyzz, g_0_y_yyyz_xxyzz, g_0_y_yyyz_xxyzzz, g_0_y_yyyz_xxzzz, g_0_y_yyyz_xxzzzz, g_0_y_yyyz_xyyyy, g_0_y_yyyz_xyyyyz, g_0_y_yyyz_xyyyz, g_0_y_yyyz_xyyyzz, g_0_y_yyyz_xyyzz, g_0_y_yyyz_xyyzzz, g_0_y_yyyz_xyzzz, g_0_y_yyyz_xyzzzz, g_0_y_yyyz_xzzzz, g_0_y_yyyz_xzzzzz, g_0_y_yyyz_yyyyy, g_0_y_yyyz_yyyyyz, g_0_y_yyyz_yyyyz, g_0_y_yyyz_yyyyzz, g_0_y_yyyz_yyyzz, g_0_y_yyyz_yyyzzz, g_0_y_yyyz_yyzzz, g_0_y_yyyz_yyzzzz, g_0_y_yyyz_yzzzz, g_0_y_yyyz_yzzzzz, g_0_y_yyyz_zzzzz, g_0_y_yyyz_zzzzzz, g_0_y_yyyzz_xxxxx, g_0_y_yyyzz_xxxxy, g_0_y_yyyzz_xxxxz, g_0_y_yyyzz_xxxyy, g_0_y_yyyzz_xxxyz, g_0_y_yyyzz_xxxzz, g_0_y_yyyzz_xxyyy, g_0_y_yyyzz_xxyyz, g_0_y_yyyzz_xxyzz, g_0_y_yyyzz_xxzzz, g_0_y_yyyzz_xyyyy, g_0_y_yyyzz_xyyyz, g_0_y_yyyzz_xyyzz, g_0_y_yyyzz_xyzzz, g_0_y_yyyzz_xzzzz, g_0_y_yyyzz_yyyyy, g_0_y_yyyzz_yyyyz, g_0_y_yyyzz_yyyzz, g_0_y_yyyzz_yyzzz, g_0_y_yyyzz_yzzzz, g_0_y_yyyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzz_xxxxx[k] = -g_0_y_yyyz_xxxxx[k] * ab_z + g_0_y_yyyz_xxxxxz[k];

                g_0_y_yyyzz_xxxxy[k] = -g_0_y_yyyz_xxxxy[k] * ab_z + g_0_y_yyyz_xxxxyz[k];

                g_0_y_yyyzz_xxxxz[k] = -g_0_y_yyyz_xxxxz[k] * ab_z + g_0_y_yyyz_xxxxzz[k];

                g_0_y_yyyzz_xxxyy[k] = -g_0_y_yyyz_xxxyy[k] * ab_z + g_0_y_yyyz_xxxyyz[k];

                g_0_y_yyyzz_xxxyz[k] = -g_0_y_yyyz_xxxyz[k] * ab_z + g_0_y_yyyz_xxxyzz[k];

                g_0_y_yyyzz_xxxzz[k] = -g_0_y_yyyz_xxxzz[k] * ab_z + g_0_y_yyyz_xxxzzz[k];

                g_0_y_yyyzz_xxyyy[k] = -g_0_y_yyyz_xxyyy[k] * ab_z + g_0_y_yyyz_xxyyyz[k];

                g_0_y_yyyzz_xxyyz[k] = -g_0_y_yyyz_xxyyz[k] * ab_z + g_0_y_yyyz_xxyyzz[k];

                g_0_y_yyyzz_xxyzz[k] = -g_0_y_yyyz_xxyzz[k] * ab_z + g_0_y_yyyz_xxyzzz[k];

                g_0_y_yyyzz_xxzzz[k] = -g_0_y_yyyz_xxzzz[k] * ab_z + g_0_y_yyyz_xxzzzz[k];

                g_0_y_yyyzz_xyyyy[k] = -g_0_y_yyyz_xyyyy[k] * ab_z + g_0_y_yyyz_xyyyyz[k];

                g_0_y_yyyzz_xyyyz[k] = -g_0_y_yyyz_xyyyz[k] * ab_z + g_0_y_yyyz_xyyyzz[k];

                g_0_y_yyyzz_xyyzz[k] = -g_0_y_yyyz_xyyzz[k] * ab_z + g_0_y_yyyz_xyyzzz[k];

                g_0_y_yyyzz_xyzzz[k] = -g_0_y_yyyz_xyzzz[k] * ab_z + g_0_y_yyyz_xyzzzz[k];

                g_0_y_yyyzz_xzzzz[k] = -g_0_y_yyyz_xzzzz[k] * ab_z + g_0_y_yyyz_xzzzzz[k];

                g_0_y_yyyzz_yyyyy[k] = -g_0_y_yyyz_yyyyy[k] * ab_z + g_0_y_yyyz_yyyyyz[k];

                g_0_y_yyyzz_yyyyz[k] = -g_0_y_yyyz_yyyyz[k] * ab_z + g_0_y_yyyz_yyyyzz[k];

                g_0_y_yyyzz_yyyzz[k] = -g_0_y_yyyz_yyyzz[k] * ab_z + g_0_y_yyyz_yyyzzz[k];

                g_0_y_yyyzz_yyzzz[k] = -g_0_y_yyyz_yyzzz[k] * ab_z + g_0_y_yyyz_yyzzzz[k];

                g_0_y_yyyzz_yzzzz[k] = -g_0_y_yyyz_yzzzz[k] * ab_z + g_0_y_yyyz_yzzzzz[k];

                g_0_y_yyyzz_zzzzz[k] = -g_0_y_yyyz_zzzzz[k] * ab_z + g_0_y_yyyz_zzzzzz[k];
            }

            /// Set up 819-840 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzz_xxxxx = cbuffer.data(hh_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxy = cbuffer.data(hh_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxxz = cbuffer.data(hh_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxyy = cbuffer.data(hh_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxyz = cbuffer.data(hh_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxzz = cbuffer.data(hh_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxyyy = cbuffer.data(hh_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxyyz = cbuffer.data(hh_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxyzz = cbuffer.data(hh_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxzzz = cbuffer.data(hh_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyyyy = cbuffer.data(hh_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyyyz = cbuffer.data(hh_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyyzz = cbuffer.data(hh_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyzzz = cbuffer.data(hh_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_y_yyzzz_xzzzz = cbuffer.data(hh_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyyyy = cbuffer.data(hh_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyyyz = cbuffer.data(hh_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyyzz = cbuffer.data(hh_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyzzz = cbuffer.data(hh_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_y_yyzzz_yzzzz = cbuffer.data(hh_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_y_yyzzz_zzzzz = cbuffer.data(hh_geom_01_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzz_xxxxx, g_0_y_yyzz_xxxxxz, g_0_y_yyzz_xxxxy, g_0_y_yyzz_xxxxyz, g_0_y_yyzz_xxxxz, g_0_y_yyzz_xxxxzz, g_0_y_yyzz_xxxyy, g_0_y_yyzz_xxxyyz, g_0_y_yyzz_xxxyz, g_0_y_yyzz_xxxyzz, g_0_y_yyzz_xxxzz, g_0_y_yyzz_xxxzzz, g_0_y_yyzz_xxyyy, g_0_y_yyzz_xxyyyz, g_0_y_yyzz_xxyyz, g_0_y_yyzz_xxyyzz, g_0_y_yyzz_xxyzz, g_0_y_yyzz_xxyzzz, g_0_y_yyzz_xxzzz, g_0_y_yyzz_xxzzzz, g_0_y_yyzz_xyyyy, g_0_y_yyzz_xyyyyz, g_0_y_yyzz_xyyyz, g_0_y_yyzz_xyyyzz, g_0_y_yyzz_xyyzz, g_0_y_yyzz_xyyzzz, g_0_y_yyzz_xyzzz, g_0_y_yyzz_xyzzzz, g_0_y_yyzz_xzzzz, g_0_y_yyzz_xzzzzz, g_0_y_yyzz_yyyyy, g_0_y_yyzz_yyyyyz, g_0_y_yyzz_yyyyz, g_0_y_yyzz_yyyyzz, g_0_y_yyzz_yyyzz, g_0_y_yyzz_yyyzzz, g_0_y_yyzz_yyzzz, g_0_y_yyzz_yyzzzz, g_0_y_yyzz_yzzzz, g_0_y_yyzz_yzzzzz, g_0_y_yyzz_zzzzz, g_0_y_yyzz_zzzzzz, g_0_y_yyzzz_xxxxx, g_0_y_yyzzz_xxxxy, g_0_y_yyzzz_xxxxz, g_0_y_yyzzz_xxxyy, g_0_y_yyzzz_xxxyz, g_0_y_yyzzz_xxxzz, g_0_y_yyzzz_xxyyy, g_0_y_yyzzz_xxyyz, g_0_y_yyzzz_xxyzz, g_0_y_yyzzz_xxzzz, g_0_y_yyzzz_xyyyy, g_0_y_yyzzz_xyyyz, g_0_y_yyzzz_xyyzz, g_0_y_yyzzz_xyzzz, g_0_y_yyzzz_xzzzz, g_0_y_yyzzz_yyyyy, g_0_y_yyzzz_yyyyz, g_0_y_yyzzz_yyyzz, g_0_y_yyzzz_yyzzz, g_0_y_yyzzz_yzzzz, g_0_y_yyzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzz_xxxxx[k] = -g_0_y_yyzz_xxxxx[k] * ab_z + g_0_y_yyzz_xxxxxz[k];

                g_0_y_yyzzz_xxxxy[k] = -g_0_y_yyzz_xxxxy[k] * ab_z + g_0_y_yyzz_xxxxyz[k];

                g_0_y_yyzzz_xxxxz[k] = -g_0_y_yyzz_xxxxz[k] * ab_z + g_0_y_yyzz_xxxxzz[k];

                g_0_y_yyzzz_xxxyy[k] = -g_0_y_yyzz_xxxyy[k] * ab_z + g_0_y_yyzz_xxxyyz[k];

                g_0_y_yyzzz_xxxyz[k] = -g_0_y_yyzz_xxxyz[k] * ab_z + g_0_y_yyzz_xxxyzz[k];

                g_0_y_yyzzz_xxxzz[k] = -g_0_y_yyzz_xxxzz[k] * ab_z + g_0_y_yyzz_xxxzzz[k];

                g_0_y_yyzzz_xxyyy[k] = -g_0_y_yyzz_xxyyy[k] * ab_z + g_0_y_yyzz_xxyyyz[k];

                g_0_y_yyzzz_xxyyz[k] = -g_0_y_yyzz_xxyyz[k] * ab_z + g_0_y_yyzz_xxyyzz[k];

                g_0_y_yyzzz_xxyzz[k] = -g_0_y_yyzz_xxyzz[k] * ab_z + g_0_y_yyzz_xxyzzz[k];

                g_0_y_yyzzz_xxzzz[k] = -g_0_y_yyzz_xxzzz[k] * ab_z + g_0_y_yyzz_xxzzzz[k];

                g_0_y_yyzzz_xyyyy[k] = -g_0_y_yyzz_xyyyy[k] * ab_z + g_0_y_yyzz_xyyyyz[k];

                g_0_y_yyzzz_xyyyz[k] = -g_0_y_yyzz_xyyyz[k] * ab_z + g_0_y_yyzz_xyyyzz[k];

                g_0_y_yyzzz_xyyzz[k] = -g_0_y_yyzz_xyyzz[k] * ab_z + g_0_y_yyzz_xyyzzz[k];

                g_0_y_yyzzz_xyzzz[k] = -g_0_y_yyzz_xyzzz[k] * ab_z + g_0_y_yyzz_xyzzzz[k];

                g_0_y_yyzzz_xzzzz[k] = -g_0_y_yyzz_xzzzz[k] * ab_z + g_0_y_yyzz_xzzzzz[k];

                g_0_y_yyzzz_yyyyy[k] = -g_0_y_yyzz_yyyyy[k] * ab_z + g_0_y_yyzz_yyyyyz[k];

                g_0_y_yyzzz_yyyyz[k] = -g_0_y_yyzz_yyyyz[k] * ab_z + g_0_y_yyzz_yyyyzz[k];

                g_0_y_yyzzz_yyyzz[k] = -g_0_y_yyzz_yyyzz[k] * ab_z + g_0_y_yyzz_yyyzzz[k];

                g_0_y_yyzzz_yyzzz[k] = -g_0_y_yyzz_yyzzz[k] * ab_z + g_0_y_yyzz_yyzzzz[k];

                g_0_y_yyzzz_yzzzz[k] = -g_0_y_yyzz_yzzzz[k] * ab_z + g_0_y_yyzz_yzzzzz[k];

                g_0_y_yyzzz_zzzzz[k] = -g_0_y_yyzz_zzzzz[k] * ab_z + g_0_y_yyzz_zzzzzz[k];
            }

            /// Set up 840-861 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzz_xxxxx = cbuffer.data(hh_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxy = cbuffer.data(hh_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxxz = cbuffer.data(hh_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxyy = cbuffer.data(hh_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxyz = cbuffer.data(hh_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxzz = cbuffer.data(hh_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxyyy = cbuffer.data(hh_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxyyz = cbuffer.data(hh_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxyzz = cbuffer.data(hh_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxzzz = cbuffer.data(hh_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyyyy = cbuffer.data(hh_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyyyz = cbuffer.data(hh_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyyzz = cbuffer.data(hh_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyzzz = cbuffer.data(hh_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_y_yzzzz_xzzzz = cbuffer.data(hh_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyyyy = cbuffer.data(hh_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyyyz = cbuffer.data(hh_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyyzz = cbuffer.data(hh_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyzzz = cbuffer.data(hh_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_y_yzzzz_yzzzz = cbuffer.data(hh_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_y_yzzzz_zzzzz = cbuffer.data(hh_geom_01_off + 860 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzz_xxxxx, g_0_y_yzzz_xxxxxz, g_0_y_yzzz_xxxxy, g_0_y_yzzz_xxxxyz, g_0_y_yzzz_xxxxz, g_0_y_yzzz_xxxxzz, g_0_y_yzzz_xxxyy, g_0_y_yzzz_xxxyyz, g_0_y_yzzz_xxxyz, g_0_y_yzzz_xxxyzz, g_0_y_yzzz_xxxzz, g_0_y_yzzz_xxxzzz, g_0_y_yzzz_xxyyy, g_0_y_yzzz_xxyyyz, g_0_y_yzzz_xxyyz, g_0_y_yzzz_xxyyzz, g_0_y_yzzz_xxyzz, g_0_y_yzzz_xxyzzz, g_0_y_yzzz_xxzzz, g_0_y_yzzz_xxzzzz, g_0_y_yzzz_xyyyy, g_0_y_yzzz_xyyyyz, g_0_y_yzzz_xyyyz, g_0_y_yzzz_xyyyzz, g_0_y_yzzz_xyyzz, g_0_y_yzzz_xyyzzz, g_0_y_yzzz_xyzzz, g_0_y_yzzz_xyzzzz, g_0_y_yzzz_xzzzz, g_0_y_yzzz_xzzzzz, g_0_y_yzzz_yyyyy, g_0_y_yzzz_yyyyyz, g_0_y_yzzz_yyyyz, g_0_y_yzzz_yyyyzz, g_0_y_yzzz_yyyzz, g_0_y_yzzz_yyyzzz, g_0_y_yzzz_yyzzz, g_0_y_yzzz_yyzzzz, g_0_y_yzzz_yzzzz, g_0_y_yzzz_yzzzzz, g_0_y_yzzz_zzzzz, g_0_y_yzzz_zzzzzz, g_0_y_yzzzz_xxxxx, g_0_y_yzzzz_xxxxy, g_0_y_yzzzz_xxxxz, g_0_y_yzzzz_xxxyy, g_0_y_yzzzz_xxxyz, g_0_y_yzzzz_xxxzz, g_0_y_yzzzz_xxyyy, g_0_y_yzzzz_xxyyz, g_0_y_yzzzz_xxyzz, g_0_y_yzzzz_xxzzz, g_0_y_yzzzz_xyyyy, g_0_y_yzzzz_xyyyz, g_0_y_yzzzz_xyyzz, g_0_y_yzzzz_xyzzz, g_0_y_yzzzz_xzzzz, g_0_y_yzzzz_yyyyy, g_0_y_yzzzz_yyyyz, g_0_y_yzzzz_yyyzz, g_0_y_yzzzz_yyzzz, g_0_y_yzzzz_yzzzz, g_0_y_yzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzz_xxxxx[k] = -g_0_y_yzzz_xxxxx[k] * ab_z + g_0_y_yzzz_xxxxxz[k];

                g_0_y_yzzzz_xxxxy[k] = -g_0_y_yzzz_xxxxy[k] * ab_z + g_0_y_yzzz_xxxxyz[k];

                g_0_y_yzzzz_xxxxz[k] = -g_0_y_yzzz_xxxxz[k] * ab_z + g_0_y_yzzz_xxxxzz[k];

                g_0_y_yzzzz_xxxyy[k] = -g_0_y_yzzz_xxxyy[k] * ab_z + g_0_y_yzzz_xxxyyz[k];

                g_0_y_yzzzz_xxxyz[k] = -g_0_y_yzzz_xxxyz[k] * ab_z + g_0_y_yzzz_xxxyzz[k];

                g_0_y_yzzzz_xxxzz[k] = -g_0_y_yzzz_xxxzz[k] * ab_z + g_0_y_yzzz_xxxzzz[k];

                g_0_y_yzzzz_xxyyy[k] = -g_0_y_yzzz_xxyyy[k] * ab_z + g_0_y_yzzz_xxyyyz[k];

                g_0_y_yzzzz_xxyyz[k] = -g_0_y_yzzz_xxyyz[k] * ab_z + g_0_y_yzzz_xxyyzz[k];

                g_0_y_yzzzz_xxyzz[k] = -g_0_y_yzzz_xxyzz[k] * ab_z + g_0_y_yzzz_xxyzzz[k];

                g_0_y_yzzzz_xxzzz[k] = -g_0_y_yzzz_xxzzz[k] * ab_z + g_0_y_yzzz_xxzzzz[k];

                g_0_y_yzzzz_xyyyy[k] = -g_0_y_yzzz_xyyyy[k] * ab_z + g_0_y_yzzz_xyyyyz[k];

                g_0_y_yzzzz_xyyyz[k] = -g_0_y_yzzz_xyyyz[k] * ab_z + g_0_y_yzzz_xyyyzz[k];

                g_0_y_yzzzz_xyyzz[k] = -g_0_y_yzzz_xyyzz[k] * ab_z + g_0_y_yzzz_xyyzzz[k];

                g_0_y_yzzzz_xyzzz[k] = -g_0_y_yzzz_xyzzz[k] * ab_z + g_0_y_yzzz_xyzzzz[k];

                g_0_y_yzzzz_xzzzz[k] = -g_0_y_yzzz_xzzzz[k] * ab_z + g_0_y_yzzz_xzzzzz[k];

                g_0_y_yzzzz_yyyyy[k] = -g_0_y_yzzz_yyyyy[k] * ab_z + g_0_y_yzzz_yyyyyz[k];

                g_0_y_yzzzz_yyyyz[k] = -g_0_y_yzzz_yyyyz[k] * ab_z + g_0_y_yzzz_yyyyzz[k];

                g_0_y_yzzzz_yyyzz[k] = -g_0_y_yzzz_yyyzz[k] * ab_z + g_0_y_yzzz_yyyzzz[k];

                g_0_y_yzzzz_yyzzz[k] = -g_0_y_yzzz_yyzzz[k] * ab_z + g_0_y_yzzz_yyzzzz[k];

                g_0_y_yzzzz_yzzzz[k] = -g_0_y_yzzz_yzzzz[k] * ab_z + g_0_y_yzzz_yzzzzz[k];

                g_0_y_yzzzz_zzzzz[k] = -g_0_y_yzzz_zzzzz[k] * ab_z + g_0_y_yzzz_zzzzzz[k];
            }

            /// Set up 861-882 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzz_xxxxx = cbuffer.data(hh_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxy = cbuffer.data(hh_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxxz = cbuffer.data(hh_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxyy = cbuffer.data(hh_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxyz = cbuffer.data(hh_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxzz = cbuffer.data(hh_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxyyy = cbuffer.data(hh_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxyyz = cbuffer.data(hh_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxyzz = cbuffer.data(hh_geom_01_off + 869 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxzzz = cbuffer.data(hh_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyyyy = cbuffer.data(hh_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyyyz = cbuffer.data(hh_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyyzz = cbuffer.data(hh_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyzzz = cbuffer.data(hh_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_y_zzzzz_xzzzz = cbuffer.data(hh_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyyyy = cbuffer.data(hh_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyyyz = cbuffer.data(hh_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyyzz = cbuffer.data(hh_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyzzz = cbuffer.data(hh_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_y_zzzzz_yzzzz = cbuffer.data(hh_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_y_zzzzz_zzzzz = cbuffer.data(hh_geom_01_off + 881 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzz_xxxxx, g_0_y_zzzz_xxxxxz, g_0_y_zzzz_xxxxy, g_0_y_zzzz_xxxxyz, g_0_y_zzzz_xxxxz, g_0_y_zzzz_xxxxzz, g_0_y_zzzz_xxxyy, g_0_y_zzzz_xxxyyz, g_0_y_zzzz_xxxyz, g_0_y_zzzz_xxxyzz, g_0_y_zzzz_xxxzz, g_0_y_zzzz_xxxzzz, g_0_y_zzzz_xxyyy, g_0_y_zzzz_xxyyyz, g_0_y_zzzz_xxyyz, g_0_y_zzzz_xxyyzz, g_0_y_zzzz_xxyzz, g_0_y_zzzz_xxyzzz, g_0_y_zzzz_xxzzz, g_0_y_zzzz_xxzzzz, g_0_y_zzzz_xyyyy, g_0_y_zzzz_xyyyyz, g_0_y_zzzz_xyyyz, g_0_y_zzzz_xyyyzz, g_0_y_zzzz_xyyzz, g_0_y_zzzz_xyyzzz, g_0_y_zzzz_xyzzz, g_0_y_zzzz_xyzzzz, g_0_y_zzzz_xzzzz, g_0_y_zzzz_xzzzzz, g_0_y_zzzz_yyyyy, g_0_y_zzzz_yyyyyz, g_0_y_zzzz_yyyyz, g_0_y_zzzz_yyyyzz, g_0_y_zzzz_yyyzz, g_0_y_zzzz_yyyzzz, g_0_y_zzzz_yyzzz, g_0_y_zzzz_yyzzzz, g_0_y_zzzz_yzzzz, g_0_y_zzzz_yzzzzz, g_0_y_zzzz_zzzzz, g_0_y_zzzz_zzzzzz, g_0_y_zzzzz_xxxxx, g_0_y_zzzzz_xxxxy, g_0_y_zzzzz_xxxxz, g_0_y_zzzzz_xxxyy, g_0_y_zzzzz_xxxyz, g_0_y_zzzzz_xxxzz, g_0_y_zzzzz_xxyyy, g_0_y_zzzzz_xxyyz, g_0_y_zzzzz_xxyzz, g_0_y_zzzzz_xxzzz, g_0_y_zzzzz_xyyyy, g_0_y_zzzzz_xyyyz, g_0_y_zzzzz_xyyzz, g_0_y_zzzzz_xyzzz, g_0_y_zzzzz_xzzzz, g_0_y_zzzzz_yyyyy, g_0_y_zzzzz_yyyyz, g_0_y_zzzzz_yyyzz, g_0_y_zzzzz_yyzzz, g_0_y_zzzzz_yzzzz, g_0_y_zzzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzz_xxxxx[k] = -g_0_y_zzzz_xxxxx[k] * ab_z + g_0_y_zzzz_xxxxxz[k];

                g_0_y_zzzzz_xxxxy[k] = -g_0_y_zzzz_xxxxy[k] * ab_z + g_0_y_zzzz_xxxxyz[k];

                g_0_y_zzzzz_xxxxz[k] = -g_0_y_zzzz_xxxxz[k] * ab_z + g_0_y_zzzz_xxxxzz[k];

                g_0_y_zzzzz_xxxyy[k] = -g_0_y_zzzz_xxxyy[k] * ab_z + g_0_y_zzzz_xxxyyz[k];

                g_0_y_zzzzz_xxxyz[k] = -g_0_y_zzzz_xxxyz[k] * ab_z + g_0_y_zzzz_xxxyzz[k];

                g_0_y_zzzzz_xxxzz[k] = -g_0_y_zzzz_xxxzz[k] * ab_z + g_0_y_zzzz_xxxzzz[k];

                g_0_y_zzzzz_xxyyy[k] = -g_0_y_zzzz_xxyyy[k] * ab_z + g_0_y_zzzz_xxyyyz[k];

                g_0_y_zzzzz_xxyyz[k] = -g_0_y_zzzz_xxyyz[k] * ab_z + g_0_y_zzzz_xxyyzz[k];

                g_0_y_zzzzz_xxyzz[k] = -g_0_y_zzzz_xxyzz[k] * ab_z + g_0_y_zzzz_xxyzzz[k];

                g_0_y_zzzzz_xxzzz[k] = -g_0_y_zzzz_xxzzz[k] * ab_z + g_0_y_zzzz_xxzzzz[k];

                g_0_y_zzzzz_xyyyy[k] = -g_0_y_zzzz_xyyyy[k] * ab_z + g_0_y_zzzz_xyyyyz[k];

                g_0_y_zzzzz_xyyyz[k] = -g_0_y_zzzz_xyyyz[k] * ab_z + g_0_y_zzzz_xyyyzz[k];

                g_0_y_zzzzz_xyyzz[k] = -g_0_y_zzzz_xyyzz[k] * ab_z + g_0_y_zzzz_xyyzzz[k];

                g_0_y_zzzzz_xyzzz[k] = -g_0_y_zzzz_xyzzz[k] * ab_z + g_0_y_zzzz_xyzzzz[k];

                g_0_y_zzzzz_xzzzz[k] = -g_0_y_zzzz_xzzzz[k] * ab_z + g_0_y_zzzz_xzzzzz[k];

                g_0_y_zzzzz_yyyyy[k] = -g_0_y_zzzz_yyyyy[k] * ab_z + g_0_y_zzzz_yyyyyz[k];

                g_0_y_zzzzz_yyyyz[k] = -g_0_y_zzzz_yyyyz[k] * ab_z + g_0_y_zzzz_yyyyzz[k];

                g_0_y_zzzzz_yyyzz[k] = -g_0_y_zzzz_yyyzz[k] * ab_z + g_0_y_zzzz_yyyzzz[k];

                g_0_y_zzzzz_yyzzz[k] = -g_0_y_zzzz_yyzzz[k] * ab_z + g_0_y_zzzz_yyzzzz[k];

                g_0_y_zzzzz_yzzzz[k] = -g_0_y_zzzz_yzzzz[k] * ab_z + g_0_y_zzzz_yzzzzz[k];

                g_0_y_zzzzz_zzzzz[k] = -g_0_y_zzzz_zzzzz[k] * ab_z + g_0_y_zzzz_zzzzzz[k];
            }

            /// Set up 882-903 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxx_xxxxx = cbuffer.data(hh_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxy = cbuffer.data(hh_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxxz = cbuffer.data(hh_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxyy = cbuffer.data(hh_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxyz = cbuffer.data(hh_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxzz = cbuffer.data(hh_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxyyy = cbuffer.data(hh_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxyyz = cbuffer.data(hh_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxyzz = cbuffer.data(hh_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxzzz = cbuffer.data(hh_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyyyy = cbuffer.data(hh_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyyyz = cbuffer.data(hh_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyyzz = cbuffer.data(hh_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyzzz = cbuffer.data(hh_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_z_xxxxx_xzzzz = cbuffer.data(hh_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyyyy = cbuffer.data(hh_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyyyz = cbuffer.data(hh_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyyzz = cbuffer.data(hh_geom_01_off + 899 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyzzz = cbuffer.data(hh_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_z_xxxxx_yzzzz = cbuffer.data(hh_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_z_xxxxx_zzzzz = cbuffer.data(hh_geom_01_off + 902 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxx_xxxxx, g_0_z_xxxx_xxxxxx, g_0_z_xxxx_xxxxxy, g_0_z_xxxx_xxxxxz, g_0_z_xxxx_xxxxy, g_0_z_xxxx_xxxxyy, g_0_z_xxxx_xxxxyz, g_0_z_xxxx_xxxxz, g_0_z_xxxx_xxxxzz, g_0_z_xxxx_xxxyy, g_0_z_xxxx_xxxyyy, g_0_z_xxxx_xxxyyz, g_0_z_xxxx_xxxyz, g_0_z_xxxx_xxxyzz, g_0_z_xxxx_xxxzz, g_0_z_xxxx_xxxzzz, g_0_z_xxxx_xxyyy, g_0_z_xxxx_xxyyyy, g_0_z_xxxx_xxyyyz, g_0_z_xxxx_xxyyz, g_0_z_xxxx_xxyyzz, g_0_z_xxxx_xxyzz, g_0_z_xxxx_xxyzzz, g_0_z_xxxx_xxzzz, g_0_z_xxxx_xxzzzz, g_0_z_xxxx_xyyyy, g_0_z_xxxx_xyyyyy, g_0_z_xxxx_xyyyyz, g_0_z_xxxx_xyyyz, g_0_z_xxxx_xyyyzz, g_0_z_xxxx_xyyzz, g_0_z_xxxx_xyyzzz, g_0_z_xxxx_xyzzz, g_0_z_xxxx_xyzzzz, g_0_z_xxxx_xzzzz, g_0_z_xxxx_xzzzzz, g_0_z_xxxx_yyyyy, g_0_z_xxxx_yyyyz, g_0_z_xxxx_yyyzz, g_0_z_xxxx_yyzzz, g_0_z_xxxx_yzzzz, g_0_z_xxxx_zzzzz, g_0_z_xxxxx_xxxxx, g_0_z_xxxxx_xxxxy, g_0_z_xxxxx_xxxxz, g_0_z_xxxxx_xxxyy, g_0_z_xxxxx_xxxyz, g_0_z_xxxxx_xxxzz, g_0_z_xxxxx_xxyyy, g_0_z_xxxxx_xxyyz, g_0_z_xxxxx_xxyzz, g_0_z_xxxxx_xxzzz, g_0_z_xxxxx_xyyyy, g_0_z_xxxxx_xyyyz, g_0_z_xxxxx_xyyzz, g_0_z_xxxxx_xyzzz, g_0_z_xxxxx_xzzzz, g_0_z_xxxxx_yyyyy, g_0_z_xxxxx_yyyyz, g_0_z_xxxxx_yyyzz, g_0_z_xxxxx_yyzzz, g_0_z_xxxxx_yzzzz, g_0_z_xxxxx_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxx_xxxxx[k] = -g_0_z_xxxx_xxxxx[k] * ab_x + g_0_z_xxxx_xxxxxx[k];

                g_0_z_xxxxx_xxxxy[k] = -g_0_z_xxxx_xxxxy[k] * ab_x + g_0_z_xxxx_xxxxxy[k];

                g_0_z_xxxxx_xxxxz[k] = -g_0_z_xxxx_xxxxz[k] * ab_x + g_0_z_xxxx_xxxxxz[k];

                g_0_z_xxxxx_xxxyy[k] = -g_0_z_xxxx_xxxyy[k] * ab_x + g_0_z_xxxx_xxxxyy[k];

                g_0_z_xxxxx_xxxyz[k] = -g_0_z_xxxx_xxxyz[k] * ab_x + g_0_z_xxxx_xxxxyz[k];

                g_0_z_xxxxx_xxxzz[k] = -g_0_z_xxxx_xxxzz[k] * ab_x + g_0_z_xxxx_xxxxzz[k];

                g_0_z_xxxxx_xxyyy[k] = -g_0_z_xxxx_xxyyy[k] * ab_x + g_0_z_xxxx_xxxyyy[k];

                g_0_z_xxxxx_xxyyz[k] = -g_0_z_xxxx_xxyyz[k] * ab_x + g_0_z_xxxx_xxxyyz[k];

                g_0_z_xxxxx_xxyzz[k] = -g_0_z_xxxx_xxyzz[k] * ab_x + g_0_z_xxxx_xxxyzz[k];

                g_0_z_xxxxx_xxzzz[k] = -g_0_z_xxxx_xxzzz[k] * ab_x + g_0_z_xxxx_xxxzzz[k];

                g_0_z_xxxxx_xyyyy[k] = -g_0_z_xxxx_xyyyy[k] * ab_x + g_0_z_xxxx_xxyyyy[k];

                g_0_z_xxxxx_xyyyz[k] = -g_0_z_xxxx_xyyyz[k] * ab_x + g_0_z_xxxx_xxyyyz[k];

                g_0_z_xxxxx_xyyzz[k] = -g_0_z_xxxx_xyyzz[k] * ab_x + g_0_z_xxxx_xxyyzz[k];

                g_0_z_xxxxx_xyzzz[k] = -g_0_z_xxxx_xyzzz[k] * ab_x + g_0_z_xxxx_xxyzzz[k];

                g_0_z_xxxxx_xzzzz[k] = -g_0_z_xxxx_xzzzz[k] * ab_x + g_0_z_xxxx_xxzzzz[k];

                g_0_z_xxxxx_yyyyy[k] = -g_0_z_xxxx_yyyyy[k] * ab_x + g_0_z_xxxx_xyyyyy[k];

                g_0_z_xxxxx_yyyyz[k] = -g_0_z_xxxx_yyyyz[k] * ab_x + g_0_z_xxxx_xyyyyz[k];

                g_0_z_xxxxx_yyyzz[k] = -g_0_z_xxxx_yyyzz[k] * ab_x + g_0_z_xxxx_xyyyzz[k];

                g_0_z_xxxxx_yyzzz[k] = -g_0_z_xxxx_yyzzz[k] * ab_x + g_0_z_xxxx_xyyzzz[k];

                g_0_z_xxxxx_yzzzz[k] = -g_0_z_xxxx_yzzzz[k] * ab_x + g_0_z_xxxx_xyzzzz[k];

                g_0_z_xxxxx_zzzzz[k] = -g_0_z_xxxx_zzzzz[k] * ab_x + g_0_z_xxxx_xzzzzz[k];
            }

            /// Set up 903-924 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxy_xxxxx = cbuffer.data(hh_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxy = cbuffer.data(hh_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxxz = cbuffer.data(hh_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxyy = cbuffer.data(hh_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxyz = cbuffer.data(hh_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxzz = cbuffer.data(hh_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxyyy = cbuffer.data(hh_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxyyz = cbuffer.data(hh_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxyzz = cbuffer.data(hh_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxzzz = cbuffer.data(hh_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyyyy = cbuffer.data(hh_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyyyz = cbuffer.data(hh_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyyzz = cbuffer.data(hh_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyzzz = cbuffer.data(hh_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_z_xxxxy_xzzzz = cbuffer.data(hh_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyyyy = cbuffer.data(hh_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyyyz = cbuffer.data(hh_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyyzz = cbuffer.data(hh_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyzzz = cbuffer.data(hh_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_z_xxxxy_yzzzz = cbuffer.data(hh_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_z_xxxxy_zzzzz = cbuffer.data(hh_geom_01_off + 923 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxy_xxxxx, g_0_z_xxxxy_xxxxy, g_0_z_xxxxy_xxxxz, g_0_z_xxxxy_xxxyy, g_0_z_xxxxy_xxxyz, g_0_z_xxxxy_xxxzz, g_0_z_xxxxy_xxyyy, g_0_z_xxxxy_xxyyz, g_0_z_xxxxy_xxyzz, g_0_z_xxxxy_xxzzz, g_0_z_xxxxy_xyyyy, g_0_z_xxxxy_xyyyz, g_0_z_xxxxy_xyyzz, g_0_z_xxxxy_xyzzz, g_0_z_xxxxy_xzzzz, g_0_z_xxxxy_yyyyy, g_0_z_xxxxy_yyyyz, g_0_z_xxxxy_yyyzz, g_0_z_xxxxy_yyzzz, g_0_z_xxxxy_yzzzz, g_0_z_xxxxy_zzzzz, g_0_z_xxxy_xxxxx, g_0_z_xxxy_xxxxxx, g_0_z_xxxy_xxxxxy, g_0_z_xxxy_xxxxxz, g_0_z_xxxy_xxxxy, g_0_z_xxxy_xxxxyy, g_0_z_xxxy_xxxxyz, g_0_z_xxxy_xxxxz, g_0_z_xxxy_xxxxzz, g_0_z_xxxy_xxxyy, g_0_z_xxxy_xxxyyy, g_0_z_xxxy_xxxyyz, g_0_z_xxxy_xxxyz, g_0_z_xxxy_xxxyzz, g_0_z_xxxy_xxxzz, g_0_z_xxxy_xxxzzz, g_0_z_xxxy_xxyyy, g_0_z_xxxy_xxyyyy, g_0_z_xxxy_xxyyyz, g_0_z_xxxy_xxyyz, g_0_z_xxxy_xxyyzz, g_0_z_xxxy_xxyzz, g_0_z_xxxy_xxyzzz, g_0_z_xxxy_xxzzz, g_0_z_xxxy_xxzzzz, g_0_z_xxxy_xyyyy, g_0_z_xxxy_xyyyyy, g_0_z_xxxy_xyyyyz, g_0_z_xxxy_xyyyz, g_0_z_xxxy_xyyyzz, g_0_z_xxxy_xyyzz, g_0_z_xxxy_xyyzzz, g_0_z_xxxy_xyzzz, g_0_z_xxxy_xyzzzz, g_0_z_xxxy_xzzzz, g_0_z_xxxy_xzzzzz, g_0_z_xxxy_yyyyy, g_0_z_xxxy_yyyyz, g_0_z_xxxy_yyyzz, g_0_z_xxxy_yyzzz, g_0_z_xxxy_yzzzz, g_0_z_xxxy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxy_xxxxx[k] = -g_0_z_xxxy_xxxxx[k] * ab_x + g_0_z_xxxy_xxxxxx[k];

                g_0_z_xxxxy_xxxxy[k] = -g_0_z_xxxy_xxxxy[k] * ab_x + g_0_z_xxxy_xxxxxy[k];

                g_0_z_xxxxy_xxxxz[k] = -g_0_z_xxxy_xxxxz[k] * ab_x + g_0_z_xxxy_xxxxxz[k];

                g_0_z_xxxxy_xxxyy[k] = -g_0_z_xxxy_xxxyy[k] * ab_x + g_0_z_xxxy_xxxxyy[k];

                g_0_z_xxxxy_xxxyz[k] = -g_0_z_xxxy_xxxyz[k] * ab_x + g_0_z_xxxy_xxxxyz[k];

                g_0_z_xxxxy_xxxzz[k] = -g_0_z_xxxy_xxxzz[k] * ab_x + g_0_z_xxxy_xxxxzz[k];

                g_0_z_xxxxy_xxyyy[k] = -g_0_z_xxxy_xxyyy[k] * ab_x + g_0_z_xxxy_xxxyyy[k];

                g_0_z_xxxxy_xxyyz[k] = -g_0_z_xxxy_xxyyz[k] * ab_x + g_0_z_xxxy_xxxyyz[k];

                g_0_z_xxxxy_xxyzz[k] = -g_0_z_xxxy_xxyzz[k] * ab_x + g_0_z_xxxy_xxxyzz[k];

                g_0_z_xxxxy_xxzzz[k] = -g_0_z_xxxy_xxzzz[k] * ab_x + g_0_z_xxxy_xxxzzz[k];

                g_0_z_xxxxy_xyyyy[k] = -g_0_z_xxxy_xyyyy[k] * ab_x + g_0_z_xxxy_xxyyyy[k];

                g_0_z_xxxxy_xyyyz[k] = -g_0_z_xxxy_xyyyz[k] * ab_x + g_0_z_xxxy_xxyyyz[k];

                g_0_z_xxxxy_xyyzz[k] = -g_0_z_xxxy_xyyzz[k] * ab_x + g_0_z_xxxy_xxyyzz[k];

                g_0_z_xxxxy_xyzzz[k] = -g_0_z_xxxy_xyzzz[k] * ab_x + g_0_z_xxxy_xxyzzz[k];

                g_0_z_xxxxy_xzzzz[k] = -g_0_z_xxxy_xzzzz[k] * ab_x + g_0_z_xxxy_xxzzzz[k];

                g_0_z_xxxxy_yyyyy[k] = -g_0_z_xxxy_yyyyy[k] * ab_x + g_0_z_xxxy_xyyyyy[k];

                g_0_z_xxxxy_yyyyz[k] = -g_0_z_xxxy_yyyyz[k] * ab_x + g_0_z_xxxy_xyyyyz[k];

                g_0_z_xxxxy_yyyzz[k] = -g_0_z_xxxy_yyyzz[k] * ab_x + g_0_z_xxxy_xyyyzz[k];

                g_0_z_xxxxy_yyzzz[k] = -g_0_z_xxxy_yyzzz[k] * ab_x + g_0_z_xxxy_xyyzzz[k];

                g_0_z_xxxxy_yzzzz[k] = -g_0_z_xxxy_yzzzz[k] * ab_x + g_0_z_xxxy_xyzzzz[k];

                g_0_z_xxxxy_zzzzz[k] = -g_0_z_xxxy_zzzzz[k] * ab_x + g_0_z_xxxy_xzzzzz[k];
            }

            /// Set up 924-945 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxz_xxxxx = cbuffer.data(hh_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxy = cbuffer.data(hh_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxxz = cbuffer.data(hh_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxyy = cbuffer.data(hh_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxyz = cbuffer.data(hh_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxzz = cbuffer.data(hh_geom_01_off + 929 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxyyy = cbuffer.data(hh_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxyyz = cbuffer.data(hh_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxyzz = cbuffer.data(hh_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxzzz = cbuffer.data(hh_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyyyy = cbuffer.data(hh_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyyyz = cbuffer.data(hh_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyyzz = cbuffer.data(hh_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyzzz = cbuffer.data(hh_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_z_xxxxz_xzzzz = cbuffer.data(hh_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyyyy = cbuffer.data(hh_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyyyz = cbuffer.data(hh_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyyzz = cbuffer.data(hh_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyzzz = cbuffer.data(hh_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_z_xxxxz_yzzzz = cbuffer.data(hh_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_z_xxxxz_zzzzz = cbuffer.data(hh_geom_01_off + 944 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxz_xxxxx, g_0_z_xxxxz_xxxxy, g_0_z_xxxxz_xxxxz, g_0_z_xxxxz_xxxyy, g_0_z_xxxxz_xxxyz, g_0_z_xxxxz_xxxzz, g_0_z_xxxxz_xxyyy, g_0_z_xxxxz_xxyyz, g_0_z_xxxxz_xxyzz, g_0_z_xxxxz_xxzzz, g_0_z_xxxxz_xyyyy, g_0_z_xxxxz_xyyyz, g_0_z_xxxxz_xyyzz, g_0_z_xxxxz_xyzzz, g_0_z_xxxxz_xzzzz, g_0_z_xxxxz_yyyyy, g_0_z_xxxxz_yyyyz, g_0_z_xxxxz_yyyzz, g_0_z_xxxxz_yyzzz, g_0_z_xxxxz_yzzzz, g_0_z_xxxxz_zzzzz, g_0_z_xxxz_xxxxx, g_0_z_xxxz_xxxxxx, g_0_z_xxxz_xxxxxy, g_0_z_xxxz_xxxxxz, g_0_z_xxxz_xxxxy, g_0_z_xxxz_xxxxyy, g_0_z_xxxz_xxxxyz, g_0_z_xxxz_xxxxz, g_0_z_xxxz_xxxxzz, g_0_z_xxxz_xxxyy, g_0_z_xxxz_xxxyyy, g_0_z_xxxz_xxxyyz, g_0_z_xxxz_xxxyz, g_0_z_xxxz_xxxyzz, g_0_z_xxxz_xxxzz, g_0_z_xxxz_xxxzzz, g_0_z_xxxz_xxyyy, g_0_z_xxxz_xxyyyy, g_0_z_xxxz_xxyyyz, g_0_z_xxxz_xxyyz, g_0_z_xxxz_xxyyzz, g_0_z_xxxz_xxyzz, g_0_z_xxxz_xxyzzz, g_0_z_xxxz_xxzzz, g_0_z_xxxz_xxzzzz, g_0_z_xxxz_xyyyy, g_0_z_xxxz_xyyyyy, g_0_z_xxxz_xyyyyz, g_0_z_xxxz_xyyyz, g_0_z_xxxz_xyyyzz, g_0_z_xxxz_xyyzz, g_0_z_xxxz_xyyzzz, g_0_z_xxxz_xyzzz, g_0_z_xxxz_xyzzzz, g_0_z_xxxz_xzzzz, g_0_z_xxxz_xzzzzz, g_0_z_xxxz_yyyyy, g_0_z_xxxz_yyyyz, g_0_z_xxxz_yyyzz, g_0_z_xxxz_yyzzz, g_0_z_xxxz_yzzzz, g_0_z_xxxz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxz_xxxxx[k] = -g_0_z_xxxz_xxxxx[k] * ab_x + g_0_z_xxxz_xxxxxx[k];

                g_0_z_xxxxz_xxxxy[k] = -g_0_z_xxxz_xxxxy[k] * ab_x + g_0_z_xxxz_xxxxxy[k];

                g_0_z_xxxxz_xxxxz[k] = -g_0_z_xxxz_xxxxz[k] * ab_x + g_0_z_xxxz_xxxxxz[k];

                g_0_z_xxxxz_xxxyy[k] = -g_0_z_xxxz_xxxyy[k] * ab_x + g_0_z_xxxz_xxxxyy[k];

                g_0_z_xxxxz_xxxyz[k] = -g_0_z_xxxz_xxxyz[k] * ab_x + g_0_z_xxxz_xxxxyz[k];

                g_0_z_xxxxz_xxxzz[k] = -g_0_z_xxxz_xxxzz[k] * ab_x + g_0_z_xxxz_xxxxzz[k];

                g_0_z_xxxxz_xxyyy[k] = -g_0_z_xxxz_xxyyy[k] * ab_x + g_0_z_xxxz_xxxyyy[k];

                g_0_z_xxxxz_xxyyz[k] = -g_0_z_xxxz_xxyyz[k] * ab_x + g_0_z_xxxz_xxxyyz[k];

                g_0_z_xxxxz_xxyzz[k] = -g_0_z_xxxz_xxyzz[k] * ab_x + g_0_z_xxxz_xxxyzz[k];

                g_0_z_xxxxz_xxzzz[k] = -g_0_z_xxxz_xxzzz[k] * ab_x + g_0_z_xxxz_xxxzzz[k];

                g_0_z_xxxxz_xyyyy[k] = -g_0_z_xxxz_xyyyy[k] * ab_x + g_0_z_xxxz_xxyyyy[k];

                g_0_z_xxxxz_xyyyz[k] = -g_0_z_xxxz_xyyyz[k] * ab_x + g_0_z_xxxz_xxyyyz[k];

                g_0_z_xxxxz_xyyzz[k] = -g_0_z_xxxz_xyyzz[k] * ab_x + g_0_z_xxxz_xxyyzz[k];

                g_0_z_xxxxz_xyzzz[k] = -g_0_z_xxxz_xyzzz[k] * ab_x + g_0_z_xxxz_xxyzzz[k];

                g_0_z_xxxxz_xzzzz[k] = -g_0_z_xxxz_xzzzz[k] * ab_x + g_0_z_xxxz_xxzzzz[k];

                g_0_z_xxxxz_yyyyy[k] = -g_0_z_xxxz_yyyyy[k] * ab_x + g_0_z_xxxz_xyyyyy[k];

                g_0_z_xxxxz_yyyyz[k] = -g_0_z_xxxz_yyyyz[k] * ab_x + g_0_z_xxxz_xyyyyz[k];

                g_0_z_xxxxz_yyyzz[k] = -g_0_z_xxxz_yyyzz[k] * ab_x + g_0_z_xxxz_xyyyzz[k];

                g_0_z_xxxxz_yyzzz[k] = -g_0_z_xxxz_yyzzz[k] * ab_x + g_0_z_xxxz_xyyzzz[k];

                g_0_z_xxxxz_yzzzz[k] = -g_0_z_xxxz_yzzzz[k] * ab_x + g_0_z_xxxz_xyzzzz[k];

                g_0_z_xxxxz_zzzzz[k] = -g_0_z_xxxz_zzzzz[k] * ab_x + g_0_z_xxxz_xzzzzz[k];
            }

            /// Set up 945-966 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyy_xxxxx = cbuffer.data(hh_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxy = cbuffer.data(hh_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxxz = cbuffer.data(hh_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxyy = cbuffer.data(hh_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxyz = cbuffer.data(hh_geom_01_off + 949 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxzz = cbuffer.data(hh_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxyyy = cbuffer.data(hh_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxyyz = cbuffer.data(hh_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxyzz = cbuffer.data(hh_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxzzz = cbuffer.data(hh_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyyyy = cbuffer.data(hh_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyyyz = cbuffer.data(hh_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyyzz = cbuffer.data(hh_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyzzz = cbuffer.data(hh_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_z_xxxyy_xzzzz = cbuffer.data(hh_geom_01_off + 959 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyyyy = cbuffer.data(hh_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyyyz = cbuffer.data(hh_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyyzz = cbuffer.data(hh_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyzzz = cbuffer.data(hh_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_z_xxxyy_yzzzz = cbuffer.data(hh_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_z_xxxyy_zzzzz = cbuffer.data(hh_geom_01_off + 965 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyy_xxxxx, g_0_z_xxxyy_xxxxy, g_0_z_xxxyy_xxxxz, g_0_z_xxxyy_xxxyy, g_0_z_xxxyy_xxxyz, g_0_z_xxxyy_xxxzz, g_0_z_xxxyy_xxyyy, g_0_z_xxxyy_xxyyz, g_0_z_xxxyy_xxyzz, g_0_z_xxxyy_xxzzz, g_0_z_xxxyy_xyyyy, g_0_z_xxxyy_xyyyz, g_0_z_xxxyy_xyyzz, g_0_z_xxxyy_xyzzz, g_0_z_xxxyy_xzzzz, g_0_z_xxxyy_yyyyy, g_0_z_xxxyy_yyyyz, g_0_z_xxxyy_yyyzz, g_0_z_xxxyy_yyzzz, g_0_z_xxxyy_yzzzz, g_0_z_xxxyy_zzzzz, g_0_z_xxyy_xxxxx, g_0_z_xxyy_xxxxxx, g_0_z_xxyy_xxxxxy, g_0_z_xxyy_xxxxxz, g_0_z_xxyy_xxxxy, g_0_z_xxyy_xxxxyy, g_0_z_xxyy_xxxxyz, g_0_z_xxyy_xxxxz, g_0_z_xxyy_xxxxzz, g_0_z_xxyy_xxxyy, g_0_z_xxyy_xxxyyy, g_0_z_xxyy_xxxyyz, g_0_z_xxyy_xxxyz, g_0_z_xxyy_xxxyzz, g_0_z_xxyy_xxxzz, g_0_z_xxyy_xxxzzz, g_0_z_xxyy_xxyyy, g_0_z_xxyy_xxyyyy, g_0_z_xxyy_xxyyyz, g_0_z_xxyy_xxyyz, g_0_z_xxyy_xxyyzz, g_0_z_xxyy_xxyzz, g_0_z_xxyy_xxyzzz, g_0_z_xxyy_xxzzz, g_0_z_xxyy_xxzzzz, g_0_z_xxyy_xyyyy, g_0_z_xxyy_xyyyyy, g_0_z_xxyy_xyyyyz, g_0_z_xxyy_xyyyz, g_0_z_xxyy_xyyyzz, g_0_z_xxyy_xyyzz, g_0_z_xxyy_xyyzzz, g_0_z_xxyy_xyzzz, g_0_z_xxyy_xyzzzz, g_0_z_xxyy_xzzzz, g_0_z_xxyy_xzzzzz, g_0_z_xxyy_yyyyy, g_0_z_xxyy_yyyyz, g_0_z_xxyy_yyyzz, g_0_z_xxyy_yyzzz, g_0_z_xxyy_yzzzz, g_0_z_xxyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyy_xxxxx[k] = -g_0_z_xxyy_xxxxx[k] * ab_x + g_0_z_xxyy_xxxxxx[k];

                g_0_z_xxxyy_xxxxy[k] = -g_0_z_xxyy_xxxxy[k] * ab_x + g_0_z_xxyy_xxxxxy[k];

                g_0_z_xxxyy_xxxxz[k] = -g_0_z_xxyy_xxxxz[k] * ab_x + g_0_z_xxyy_xxxxxz[k];

                g_0_z_xxxyy_xxxyy[k] = -g_0_z_xxyy_xxxyy[k] * ab_x + g_0_z_xxyy_xxxxyy[k];

                g_0_z_xxxyy_xxxyz[k] = -g_0_z_xxyy_xxxyz[k] * ab_x + g_0_z_xxyy_xxxxyz[k];

                g_0_z_xxxyy_xxxzz[k] = -g_0_z_xxyy_xxxzz[k] * ab_x + g_0_z_xxyy_xxxxzz[k];

                g_0_z_xxxyy_xxyyy[k] = -g_0_z_xxyy_xxyyy[k] * ab_x + g_0_z_xxyy_xxxyyy[k];

                g_0_z_xxxyy_xxyyz[k] = -g_0_z_xxyy_xxyyz[k] * ab_x + g_0_z_xxyy_xxxyyz[k];

                g_0_z_xxxyy_xxyzz[k] = -g_0_z_xxyy_xxyzz[k] * ab_x + g_0_z_xxyy_xxxyzz[k];

                g_0_z_xxxyy_xxzzz[k] = -g_0_z_xxyy_xxzzz[k] * ab_x + g_0_z_xxyy_xxxzzz[k];

                g_0_z_xxxyy_xyyyy[k] = -g_0_z_xxyy_xyyyy[k] * ab_x + g_0_z_xxyy_xxyyyy[k];

                g_0_z_xxxyy_xyyyz[k] = -g_0_z_xxyy_xyyyz[k] * ab_x + g_0_z_xxyy_xxyyyz[k];

                g_0_z_xxxyy_xyyzz[k] = -g_0_z_xxyy_xyyzz[k] * ab_x + g_0_z_xxyy_xxyyzz[k];

                g_0_z_xxxyy_xyzzz[k] = -g_0_z_xxyy_xyzzz[k] * ab_x + g_0_z_xxyy_xxyzzz[k];

                g_0_z_xxxyy_xzzzz[k] = -g_0_z_xxyy_xzzzz[k] * ab_x + g_0_z_xxyy_xxzzzz[k];

                g_0_z_xxxyy_yyyyy[k] = -g_0_z_xxyy_yyyyy[k] * ab_x + g_0_z_xxyy_xyyyyy[k];

                g_0_z_xxxyy_yyyyz[k] = -g_0_z_xxyy_yyyyz[k] * ab_x + g_0_z_xxyy_xyyyyz[k];

                g_0_z_xxxyy_yyyzz[k] = -g_0_z_xxyy_yyyzz[k] * ab_x + g_0_z_xxyy_xyyyzz[k];

                g_0_z_xxxyy_yyzzz[k] = -g_0_z_xxyy_yyzzz[k] * ab_x + g_0_z_xxyy_xyyzzz[k];

                g_0_z_xxxyy_yzzzz[k] = -g_0_z_xxyy_yzzzz[k] * ab_x + g_0_z_xxyy_xyzzzz[k];

                g_0_z_xxxyy_zzzzz[k] = -g_0_z_xxyy_zzzzz[k] * ab_x + g_0_z_xxyy_xzzzzz[k];
            }

            /// Set up 966-987 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyz_xxxxx = cbuffer.data(hh_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxy = cbuffer.data(hh_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxxz = cbuffer.data(hh_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxyy = cbuffer.data(hh_geom_01_off + 969 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxyz = cbuffer.data(hh_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxzz = cbuffer.data(hh_geom_01_off + 971 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxyyy = cbuffer.data(hh_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxyyz = cbuffer.data(hh_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxyzz = cbuffer.data(hh_geom_01_off + 974 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxzzz = cbuffer.data(hh_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyyyy = cbuffer.data(hh_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyyyz = cbuffer.data(hh_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyyzz = cbuffer.data(hh_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyzzz = cbuffer.data(hh_geom_01_off + 979 * ccomps * dcomps);

            auto g_0_z_xxxyz_xzzzz = cbuffer.data(hh_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyyyy = cbuffer.data(hh_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyyyz = cbuffer.data(hh_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyyzz = cbuffer.data(hh_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyzzz = cbuffer.data(hh_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_z_xxxyz_yzzzz = cbuffer.data(hh_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_z_xxxyz_zzzzz = cbuffer.data(hh_geom_01_off + 986 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyz_xxxxx, g_0_z_xxxyz_xxxxy, g_0_z_xxxyz_xxxxz, g_0_z_xxxyz_xxxyy, g_0_z_xxxyz_xxxyz, g_0_z_xxxyz_xxxzz, g_0_z_xxxyz_xxyyy, g_0_z_xxxyz_xxyyz, g_0_z_xxxyz_xxyzz, g_0_z_xxxyz_xxzzz, g_0_z_xxxyz_xyyyy, g_0_z_xxxyz_xyyyz, g_0_z_xxxyz_xyyzz, g_0_z_xxxyz_xyzzz, g_0_z_xxxyz_xzzzz, g_0_z_xxxyz_yyyyy, g_0_z_xxxyz_yyyyz, g_0_z_xxxyz_yyyzz, g_0_z_xxxyz_yyzzz, g_0_z_xxxyz_yzzzz, g_0_z_xxxyz_zzzzz, g_0_z_xxyz_xxxxx, g_0_z_xxyz_xxxxxx, g_0_z_xxyz_xxxxxy, g_0_z_xxyz_xxxxxz, g_0_z_xxyz_xxxxy, g_0_z_xxyz_xxxxyy, g_0_z_xxyz_xxxxyz, g_0_z_xxyz_xxxxz, g_0_z_xxyz_xxxxzz, g_0_z_xxyz_xxxyy, g_0_z_xxyz_xxxyyy, g_0_z_xxyz_xxxyyz, g_0_z_xxyz_xxxyz, g_0_z_xxyz_xxxyzz, g_0_z_xxyz_xxxzz, g_0_z_xxyz_xxxzzz, g_0_z_xxyz_xxyyy, g_0_z_xxyz_xxyyyy, g_0_z_xxyz_xxyyyz, g_0_z_xxyz_xxyyz, g_0_z_xxyz_xxyyzz, g_0_z_xxyz_xxyzz, g_0_z_xxyz_xxyzzz, g_0_z_xxyz_xxzzz, g_0_z_xxyz_xxzzzz, g_0_z_xxyz_xyyyy, g_0_z_xxyz_xyyyyy, g_0_z_xxyz_xyyyyz, g_0_z_xxyz_xyyyz, g_0_z_xxyz_xyyyzz, g_0_z_xxyz_xyyzz, g_0_z_xxyz_xyyzzz, g_0_z_xxyz_xyzzz, g_0_z_xxyz_xyzzzz, g_0_z_xxyz_xzzzz, g_0_z_xxyz_xzzzzz, g_0_z_xxyz_yyyyy, g_0_z_xxyz_yyyyz, g_0_z_xxyz_yyyzz, g_0_z_xxyz_yyzzz, g_0_z_xxyz_yzzzz, g_0_z_xxyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyz_xxxxx[k] = -g_0_z_xxyz_xxxxx[k] * ab_x + g_0_z_xxyz_xxxxxx[k];

                g_0_z_xxxyz_xxxxy[k] = -g_0_z_xxyz_xxxxy[k] * ab_x + g_0_z_xxyz_xxxxxy[k];

                g_0_z_xxxyz_xxxxz[k] = -g_0_z_xxyz_xxxxz[k] * ab_x + g_0_z_xxyz_xxxxxz[k];

                g_0_z_xxxyz_xxxyy[k] = -g_0_z_xxyz_xxxyy[k] * ab_x + g_0_z_xxyz_xxxxyy[k];

                g_0_z_xxxyz_xxxyz[k] = -g_0_z_xxyz_xxxyz[k] * ab_x + g_0_z_xxyz_xxxxyz[k];

                g_0_z_xxxyz_xxxzz[k] = -g_0_z_xxyz_xxxzz[k] * ab_x + g_0_z_xxyz_xxxxzz[k];

                g_0_z_xxxyz_xxyyy[k] = -g_0_z_xxyz_xxyyy[k] * ab_x + g_0_z_xxyz_xxxyyy[k];

                g_0_z_xxxyz_xxyyz[k] = -g_0_z_xxyz_xxyyz[k] * ab_x + g_0_z_xxyz_xxxyyz[k];

                g_0_z_xxxyz_xxyzz[k] = -g_0_z_xxyz_xxyzz[k] * ab_x + g_0_z_xxyz_xxxyzz[k];

                g_0_z_xxxyz_xxzzz[k] = -g_0_z_xxyz_xxzzz[k] * ab_x + g_0_z_xxyz_xxxzzz[k];

                g_0_z_xxxyz_xyyyy[k] = -g_0_z_xxyz_xyyyy[k] * ab_x + g_0_z_xxyz_xxyyyy[k];

                g_0_z_xxxyz_xyyyz[k] = -g_0_z_xxyz_xyyyz[k] * ab_x + g_0_z_xxyz_xxyyyz[k];

                g_0_z_xxxyz_xyyzz[k] = -g_0_z_xxyz_xyyzz[k] * ab_x + g_0_z_xxyz_xxyyzz[k];

                g_0_z_xxxyz_xyzzz[k] = -g_0_z_xxyz_xyzzz[k] * ab_x + g_0_z_xxyz_xxyzzz[k];

                g_0_z_xxxyz_xzzzz[k] = -g_0_z_xxyz_xzzzz[k] * ab_x + g_0_z_xxyz_xxzzzz[k];

                g_0_z_xxxyz_yyyyy[k] = -g_0_z_xxyz_yyyyy[k] * ab_x + g_0_z_xxyz_xyyyyy[k];

                g_0_z_xxxyz_yyyyz[k] = -g_0_z_xxyz_yyyyz[k] * ab_x + g_0_z_xxyz_xyyyyz[k];

                g_0_z_xxxyz_yyyzz[k] = -g_0_z_xxyz_yyyzz[k] * ab_x + g_0_z_xxyz_xyyyzz[k];

                g_0_z_xxxyz_yyzzz[k] = -g_0_z_xxyz_yyzzz[k] * ab_x + g_0_z_xxyz_xyyzzz[k];

                g_0_z_xxxyz_yzzzz[k] = -g_0_z_xxyz_yzzzz[k] * ab_x + g_0_z_xxyz_xyzzzz[k];

                g_0_z_xxxyz_zzzzz[k] = -g_0_z_xxyz_zzzzz[k] * ab_x + g_0_z_xxyz_xzzzzz[k];
            }

            /// Set up 987-1008 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzz_xxxxx = cbuffer.data(hh_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxy = cbuffer.data(hh_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxxz = cbuffer.data(hh_geom_01_off + 989 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxyy = cbuffer.data(hh_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxyz = cbuffer.data(hh_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxzz = cbuffer.data(hh_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxyyy = cbuffer.data(hh_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxyyz = cbuffer.data(hh_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxyzz = cbuffer.data(hh_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxzzz = cbuffer.data(hh_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyyyy = cbuffer.data(hh_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyyyz = cbuffer.data(hh_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyyzz = cbuffer.data(hh_geom_01_off + 999 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyzzz = cbuffer.data(hh_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_z_xxxzz_xzzzz = cbuffer.data(hh_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyyyy = cbuffer.data(hh_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyyyz = cbuffer.data(hh_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyyzz = cbuffer.data(hh_geom_01_off + 1004 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyzzz = cbuffer.data(hh_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_z_xxxzz_yzzzz = cbuffer.data(hh_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_z_xxxzz_zzzzz = cbuffer.data(hh_geom_01_off + 1007 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzz_xxxxx, g_0_z_xxxzz_xxxxy, g_0_z_xxxzz_xxxxz, g_0_z_xxxzz_xxxyy, g_0_z_xxxzz_xxxyz, g_0_z_xxxzz_xxxzz, g_0_z_xxxzz_xxyyy, g_0_z_xxxzz_xxyyz, g_0_z_xxxzz_xxyzz, g_0_z_xxxzz_xxzzz, g_0_z_xxxzz_xyyyy, g_0_z_xxxzz_xyyyz, g_0_z_xxxzz_xyyzz, g_0_z_xxxzz_xyzzz, g_0_z_xxxzz_xzzzz, g_0_z_xxxzz_yyyyy, g_0_z_xxxzz_yyyyz, g_0_z_xxxzz_yyyzz, g_0_z_xxxzz_yyzzz, g_0_z_xxxzz_yzzzz, g_0_z_xxxzz_zzzzz, g_0_z_xxzz_xxxxx, g_0_z_xxzz_xxxxxx, g_0_z_xxzz_xxxxxy, g_0_z_xxzz_xxxxxz, g_0_z_xxzz_xxxxy, g_0_z_xxzz_xxxxyy, g_0_z_xxzz_xxxxyz, g_0_z_xxzz_xxxxz, g_0_z_xxzz_xxxxzz, g_0_z_xxzz_xxxyy, g_0_z_xxzz_xxxyyy, g_0_z_xxzz_xxxyyz, g_0_z_xxzz_xxxyz, g_0_z_xxzz_xxxyzz, g_0_z_xxzz_xxxzz, g_0_z_xxzz_xxxzzz, g_0_z_xxzz_xxyyy, g_0_z_xxzz_xxyyyy, g_0_z_xxzz_xxyyyz, g_0_z_xxzz_xxyyz, g_0_z_xxzz_xxyyzz, g_0_z_xxzz_xxyzz, g_0_z_xxzz_xxyzzz, g_0_z_xxzz_xxzzz, g_0_z_xxzz_xxzzzz, g_0_z_xxzz_xyyyy, g_0_z_xxzz_xyyyyy, g_0_z_xxzz_xyyyyz, g_0_z_xxzz_xyyyz, g_0_z_xxzz_xyyyzz, g_0_z_xxzz_xyyzz, g_0_z_xxzz_xyyzzz, g_0_z_xxzz_xyzzz, g_0_z_xxzz_xyzzzz, g_0_z_xxzz_xzzzz, g_0_z_xxzz_xzzzzz, g_0_z_xxzz_yyyyy, g_0_z_xxzz_yyyyz, g_0_z_xxzz_yyyzz, g_0_z_xxzz_yyzzz, g_0_z_xxzz_yzzzz, g_0_z_xxzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzz_xxxxx[k] = -g_0_z_xxzz_xxxxx[k] * ab_x + g_0_z_xxzz_xxxxxx[k];

                g_0_z_xxxzz_xxxxy[k] = -g_0_z_xxzz_xxxxy[k] * ab_x + g_0_z_xxzz_xxxxxy[k];

                g_0_z_xxxzz_xxxxz[k] = -g_0_z_xxzz_xxxxz[k] * ab_x + g_0_z_xxzz_xxxxxz[k];

                g_0_z_xxxzz_xxxyy[k] = -g_0_z_xxzz_xxxyy[k] * ab_x + g_0_z_xxzz_xxxxyy[k];

                g_0_z_xxxzz_xxxyz[k] = -g_0_z_xxzz_xxxyz[k] * ab_x + g_0_z_xxzz_xxxxyz[k];

                g_0_z_xxxzz_xxxzz[k] = -g_0_z_xxzz_xxxzz[k] * ab_x + g_0_z_xxzz_xxxxzz[k];

                g_0_z_xxxzz_xxyyy[k] = -g_0_z_xxzz_xxyyy[k] * ab_x + g_0_z_xxzz_xxxyyy[k];

                g_0_z_xxxzz_xxyyz[k] = -g_0_z_xxzz_xxyyz[k] * ab_x + g_0_z_xxzz_xxxyyz[k];

                g_0_z_xxxzz_xxyzz[k] = -g_0_z_xxzz_xxyzz[k] * ab_x + g_0_z_xxzz_xxxyzz[k];

                g_0_z_xxxzz_xxzzz[k] = -g_0_z_xxzz_xxzzz[k] * ab_x + g_0_z_xxzz_xxxzzz[k];

                g_0_z_xxxzz_xyyyy[k] = -g_0_z_xxzz_xyyyy[k] * ab_x + g_0_z_xxzz_xxyyyy[k];

                g_0_z_xxxzz_xyyyz[k] = -g_0_z_xxzz_xyyyz[k] * ab_x + g_0_z_xxzz_xxyyyz[k];

                g_0_z_xxxzz_xyyzz[k] = -g_0_z_xxzz_xyyzz[k] * ab_x + g_0_z_xxzz_xxyyzz[k];

                g_0_z_xxxzz_xyzzz[k] = -g_0_z_xxzz_xyzzz[k] * ab_x + g_0_z_xxzz_xxyzzz[k];

                g_0_z_xxxzz_xzzzz[k] = -g_0_z_xxzz_xzzzz[k] * ab_x + g_0_z_xxzz_xxzzzz[k];

                g_0_z_xxxzz_yyyyy[k] = -g_0_z_xxzz_yyyyy[k] * ab_x + g_0_z_xxzz_xyyyyy[k];

                g_0_z_xxxzz_yyyyz[k] = -g_0_z_xxzz_yyyyz[k] * ab_x + g_0_z_xxzz_xyyyyz[k];

                g_0_z_xxxzz_yyyzz[k] = -g_0_z_xxzz_yyyzz[k] * ab_x + g_0_z_xxzz_xyyyzz[k];

                g_0_z_xxxzz_yyzzz[k] = -g_0_z_xxzz_yyzzz[k] * ab_x + g_0_z_xxzz_xyyzzz[k];

                g_0_z_xxxzz_yzzzz[k] = -g_0_z_xxzz_yzzzz[k] * ab_x + g_0_z_xxzz_xyzzzz[k];

                g_0_z_xxxzz_zzzzz[k] = -g_0_z_xxzz_zzzzz[k] * ab_x + g_0_z_xxzz_xzzzzz[k];
            }

            /// Set up 1008-1029 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyy_xxxxx = cbuffer.data(hh_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxy = cbuffer.data(hh_geom_01_off + 1009 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxxz = cbuffer.data(hh_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxyy = cbuffer.data(hh_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxyz = cbuffer.data(hh_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxzz = cbuffer.data(hh_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxyyy = cbuffer.data(hh_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxyyz = cbuffer.data(hh_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxyzz = cbuffer.data(hh_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxzzz = cbuffer.data(hh_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyyyy = cbuffer.data(hh_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyyyz = cbuffer.data(hh_geom_01_off + 1019 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyyzz = cbuffer.data(hh_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyzzz = cbuffer.data(hh_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_z_xxyyy_xzzzz = cbuffer.data(hh_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyyyy = cbuffer.data(hh_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyyyz = cbuffer.data(hh_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyyzz = cbuffer.data(hh_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyzzz = cbuffer.data(hh_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_z_xxyyy_yzzzz = cbuffer.data(hh_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_z_xxyyy_zzzzz = cbuffer.data(hh_geom_01_off + 1028 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyy_xxxxx, g_0_z_xxyyy_xxxxy, g_0_z_xxyyy_xxxxz, g_0_z_xxyyy_xxxyy, g_0_z_xxyyy_xxxyz, g_0_z_xxyyy_xxxzz, g_0_z_xxyyy_xxyyy, g_0_z_xxyyy_xxyyz, g_0_z_xxyyy_xxyzz, g_0_z_xxyyy_xxzzz, g_0_z_xxyyy_xyyyy, g_0_z_xxyyy_xyyyz, g_0_z_xxyyy_xyyzz, g_0_z_xxyyy_xyzzz, g_0_z_xxyyy_xzzzz, g_0_z_xxyyy_yyyyy, g_0_z_xxyyy_yyyyz, g_0_z_xxyyy_yyyzz, g_0_z_xxyyy_yyzzz, g_0_z_xxyyy_yzzzz, g_0_z_xxyyy_zzzzz, g_0_z_xyyy_xxxxx, g_0_z_xyyy_xxxxxx, g_0_z_xyyy_xxxxxy, g_0_z_xyyy_xxxxxz, g_0_z_xyyy_xxxxy, g_0_z_xyyy_xxxxyy, g_0_z_xyyy_xxxxyz, g_0_z_xyyy_xxxxz, g_0_z_xyyy_xxxxzz, g_0_z_xyyy_xxxyy, g_0_z_xyyy_xxxyyy, g_0_z_xyyy_xxxyyz, g_0_z_xyyy_xxxyz, g_0_z_xyyy_xxxyzz, g_0_z_xyyy_xxxzz, g_0_z_xyyy_xxxzzz, g_0_z_xyyy_xxyyy, g_0_z_xyyy_xxyyyy, g_0_z_xyyy_xxyyyz, g_0_z_xyyy_xxyyz, g_0_z_xyyy_xxyyzz, g_0_z_xyyy_xxyzz, g_0_z_xyyy_xxyzzz, g_0_z_xyyy_xxzzz, g_0_z_xyyy_xxzzzz, g_0_z_xyyy_xyyyy, g_0_z_xyyy_xyyyyy, g_0_z_xyyy_xyyyyz, g_0_z_xyyy_xyyyz, g_0_z_xyyy_xyyyzz, g_0_z_xyyy_xyyzz, g_0_z_xyyy_xyyzzz, g_0_z_xyyy_xyzzz, g_0_z_xyyy_xyzzzz, g_0_z_xyyy_xzzzz, g_0_z_xyyy_xzzzzz, g_0_z_xyyy_yyyyy, g_0_z_xyyy_yyyyz, g_0_z_xyyy_yyyzz, g_0_z_xyyy_yyzzz, g_0_z_xyyy_yzzzz, g_0_z_xyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyy_xxxxx[k] = -g_0_z_xyyy_xxxxx[k] * ab_x + g_0_z_xyyy_xxxxxx[k];

                g_0_z_xxyyy_xxxxy[k] = -g_0_z_xyyy_xxxxy[k] * ab_x + g_0_z_xyyy_xxxxxy[k];

                g_0_z_xxyyy_xxxxz[k] = -g_0_z_xyyy_xxxxz[k] * ab_x + g_0_z_xyyy_xxxxxz[k];

                g_0_z_xxyyy_xxxyy[k] = -g_0_z_xyyy_xxxyy[k] * ab_x + g_0_z_xyyy_xxxxyy[k];

                g_0_z_xxyyy_xxxyz[k] = -g_0_z_xyyy_xxxyz[k] * ab_x + g_0_z_xyyy_xxxxyz[k];

                g_0_z_xxyyy_xxxzz[k] = -g_0_z_xyyy_xxxzz[k] * ab_x + g_0_z_xyyy_xxxxzz[k];

                g_0_z_xxyyy_xxyyy[k] = -g_0_z_xyyy_xxyyy[k] * ab_x + g_0_z_xyyy_xxxyyy[k];

                g_0_z_xxyyy_xxyyz[k] = -g_0_z_xyyy_xxyyz[k] * ab_x + g_0_z_xyyy_xxxyyz[k];

                g_0_z_xxyyy_xxyzz[k] = -g_0_z_xyyy_xxyzz[k] * ab_x + g_0_z_xyyy_xxxyzz[k];

                g_0_z_xxyyy_xxzzz[k] = -g_0_z_xyyy_xxzzz[k] * ab_x + g_0_z_xyyy_xxxzzz[k];

                g_0_z_xxyyy_xyyyy[k] = -g_0_z_xyyy_xyyyy[k] * ab_x + g_0_z_xyyy_xxyyyy[k];

                g_0_z_xxyyy_xyyyz[k] = -g_0_z_xyyy_xyyyz[k] * ab_x + g_0_z_xyyy_xxyyyz[k];

                g_0_z_xxyyy_xyyzz[k] = -g_0_z_xyyy_xyyzz[k] * ab_x + g_0_z_xyyy_xxyyzz[k];

                g_0_z_xxyyy_xyzzz[k] = -g_0_z_xyyy_xyzzz[k] * ab_x + g_0_z_xyyy_xxyzzz[k];

                g_0_z_xxyyy_xzzzz[k] = -g_0_z_xyyy_xzzzz[k] * ab_x + g_0_z_xyyy_xxzzzz[k];

                g_0_z_xxyyy_yyyyy[k] = -g_0_z_xyyy_yyyyy[k] * ab_x + g_0_z_xyyy_xyyyyy[k];

                g_0_z_xxyyy_yyyyz[k] = -g_0_z_xyyy_yyyyz[k] * ab_x + g_0_z_xyyy_xyyyyz[k];

                g_0_z_xxyyy_yyyzz[k] = -g_0_z_xyyy_yyyzz[k] * ab_x + g_0_z_xyyy_xyyyzz[k];

                g_0_z_xxyyy_yyzzz[k] = -g_0_z_xyyy_yyzzz[k] * ab_x + g_0_z_xyyy_xyyzzz[k];

                g_0_z_xxyyy_yzzzz[k] = -g_0_z_xyyy_yzzzz[k] * ab_x + g_0_z_xyyy_xyzzzz[k];

                g_0_z_xxyyy_zzzzz[k] = -g_0_z_xyyy_zzzzz[k] * ab_x + g_0_z_xyyy_xzzzzz[k];
            }

            /// Set up 1029-1050 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyz_xxxxx = cbuffer.data(hh_geom_01_off + 1029 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxy = cbuffer.data(hh_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxxz = cbuffer.data(hh_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxyy = cbuffer.data(hh_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxyz = cbuffer.data(hh_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxzz = cbuffer.data(hh_geom_01_off + 1034 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxyyy = cbuffer.data(hh_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxyyz = cbuffer.data(hh_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxyzz = cbuffer.data(hh_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxzzz = cbuffer.data(hh_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyyyy = cbuffer.data(hh_geom_01_off + 1039 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyyyz = cbuffer.data(hh_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyyzz = cbuffer.data(hh_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyzzz = cbuffer.data(hh_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_z_xxyyz_xzzzz = cbuffer.data(hh_geom_01_off + 1043 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyyyy = cbuffer.data(hh_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyyyz = cbuffer.data(hh_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyyzz = cbuffer.data(hh_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyzzz = cbuffer.data(hh_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_z_xxyyz_yzzzz = cbuffer.data(hh_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_z_xxyyz_zzzzz = cbuffer.data(hh_geom_01_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyz_xxxxx, g_0_z_xxyyz_xxxxy, g_0_z_xxyyz_xxxxz, g_0_z_xxyyz_xxxyy, g_0_z_xxyyz_xxxyz, g_0_z_xxyyz_xxxzz, g_0_z_xxyyz_xxyyy, g_0_z_xxyyz_xxyyz, g_0_z_xxyyz_xxyzz, g_0_z_xxyyz_xxzzz, g_0_z_xxyyz_xyyyy, g_0_z_xxyyz_xyyyz, g_0_z_xxyyz_xyyzz, g_0_z_xxyyz_xyzzz, g_0_z_xxyyz_xzzzz, g_0_z_xxyyz_yyyyy, g_0_z_xxyyz_yyyyz, g_0_z_xxyyz_yyyzz, g_0_z_xxyyz_yyzzz, g_0_z_xxyyz_yzzzz, g_0_z_xxyyz_zzzzz, g_0_z_xyyz_xxxxx, g_0_z_xyyz_xxxxxx, g_0_z_xyyz_xxxxxy, g_0_z_xyyz_xxxxxz, g_0_z_xyyz_xxxxy, g_0_z_xyyz_xxxxyy, g_0_z_xyyz_xxxxyz, g_0_z_xyyz_xxxxz, g_0_z_xyyz_xxxxzz, g_0_z_xyyz_xxxyy, g_0_z_xyyz_xxxyyy, g_0_z_xyyz_xxxyyz, g_0_z_xyyz_xxxyz, g_0_z_xyyz_xxxyzz, g_0_z_xyyz_xxxzz, g_0_z_xyyz_xxxzzz, g_0_z_xyyz_xxyyy, g_0_z_xyyz_xxyyyy, g_0_z_xyyz_xxyyyz, g_0_z_xyyz_xxyyz, g_0_z_xyyz_xxyyzz, g_0_z_xyyz_xxyzz, g_0_z_xyyz_xxyzzz, g_0_z_xyyz_xxzzz, g_0_z_xyyz_xxzzzz, g_0_z_xyyz_xyyyy, g_0_z_xyyz_xyyyyy, g_0_z_xyyz_xyyyyz, g_0_z_xyyz_xyyyz, g_0_z_xyyz_xyyyzz, g_0_z_xyyz_xyyzz, g_0_z_xyyz_xyyzzz, g_0_z_xyyz_xyzzz, g_0_z_xyyz_xyzzzz, g_0_z_xyyz_xzzzz, g_0_z_xyyz_xzzzzz, g_0_z_xyyz_yyyyy, g_0_z_xyyz_yyyyz, g_0_z_xyyz_yyyzz, g_0_z_xyyz_yyzzz, g_0_z_xyyz_yzzzz, g_0_z_xyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyz_xxxxx[k] = -g_0_z_xyyz_xxxxx[k] * ab_x + g_0_z_xyyz_xxxxxx[k];

                g_0_z_xxyyz_xxxxy[k] = -g_0_z_xyyz_xxxxy[k] * ab_x + g_0_z_xyyz_xxxxxy[k];

                g_0_z_xxyyz_xxxxz[k] = -g_0_z_xyyz_xxxxz[k] * ab_x + g_0_z_xyyz_xxxxxz[k];

                g_0_z_xxyyz_xxxyy[k] = -g_0_z_xyyz_xxxyy[k] * ab_x + g_0_z_xyyz_xxxxyy[k];

                g_0_z_xxyyz_xxxyz[k] = -g_0_z_xyyz_xxxyz[k] * ab_x + g_0_z_xyyz_xxxxyz[k];

                g_0_z_xxyyz_xxxzz[k] = -g_0_z_xyyz_xxxzz[k] * ab_x + g_0_z_xyyz_xxxxzz[k];

                g_0_z_xxyyz_xxyyy[k] = -g_0_z_xyyz_xxyyy[k] * ab_x + g_0_z_xyyz_xxxyyy[k];

                g_0_z_xxyyz_xxyyz[k] = -g_0_z_xyyz_xxyyz[k] * ab_x + g_0_z_xyyz_xxxyyz[k];

                g_0_z_xxyyz_xxyzz[k] = -g_0_z_xyyz_xxyzz[k] * ab_x + g_0_z_xyyz_xxxyzz[k];

                g_0_z_xxyyz_xxzzz[k] = -g_0_z_xyyz_xxzzz[k] * ab_x + g_0_z_xyyz_xxxzzz[k];

                g_0_z_xxyyz_xyyyy[k] = -g_0_z_xyyz_xyyyy[k] * ab_x + g_0_z_xyyz_xxyyyy[k];

                g_0_z_xxyyz_xyyyz[k] = -g_0_z_xyyz_xyyyz[k] * ab_x + g_0_z_xyyz_xxyyyz[k];

                g_0_z_xxyyz_xyyzz[k] = -g_0_z_xyyz_xyyzz[k] * ab_x + g_0_z_xyyz_xxyyzz[k];

                g_0_z_xxyyz_xyzzz[k] = -g_0_z_xyyz_xyzzz[k] * ab_x + g_0_z_xyyz_xxyzzz[k];

                g_0_z_xxyyz_xzzzz[k] = -g_0_z_xyyz_xzzzz[k] * ab_x + g_0_z_xyyz_xxzzzz[k];

                g_0_z_xxyyz_yyyyy[k] = -g_0_z_xyyz_yyyyy[k] * ab_x + g_0_z_xyyz_xyyyyy[k];

                g_0_z_xxyyz_yyyyz[k] = -g_0_z_xyyz_yyyyz[k] * ab_x + g_0_z_xyyz_xyyyyz[k];

                g_0_z_xxyyz_yyyzz[k] = -g_0_z_xyyz_yyyzz[k] * ab_x + g_0_z_xyyz_xyyyzz[k];

                g_0_z_xxyyz_yyzzz[k] = -g_0_z_xyyz_yyzzz[k] * ab_x + g_0_z_xyyz_xyyzzz[k];

                g_0_z_xxyyz_yzzzz[k] = -g_0_z_xyyz_yzzzz[k] * ab_x + g_0_z_xyyz_xyzzzz[k];

                g_0_z_xxyyz_zzzzz[k] = -g_0_z_xyyz_zzzzz[k] * ab_x + g_0_z_xyyz_xzzzzz[k];
            }

            /// Set up 1050-1071 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzz_xxxxx = cbuffer.data(hh_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxy = cbuffer.data(hh_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxxz = cbuffer.data(hh_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxyy = cbuffer.data(hh_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxyz = cbuffer.data(hh_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxzz = cbuffer.data(hh_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxyyy = cbuffer.data(hh_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxyyz = cbuffer.data(hh_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxyzz = cbuffer.data(hh_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxzzz = cbuffer.data(hh_geom_01_off + 1059 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyyyy = cbuffer.data(hh_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyyyz = cbuffer.data(hh_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyyzz = cbuffer.data(hh_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyzzz = cbuffer.data(hh_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_z_xxyzz_xzzzz = cbuffer.data(hh_geom_01_off + 1064 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyyyy = cbuffer.data(hh_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyyyz = cbuffer.data(hh_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyyzz = cbuffer.data(hh_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyzzz = cbuffer.data(hh_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_z_xxyzz_yzzzz = cbuffer.data(hh_geom_01_off + 1069 * ccomps * dcomps);

            auto g_0_z_xxyzz_zzzzz = cbuffer.data(hh_geom_01_off + 1070 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzz_xxxxx, g_0_z_xxyzz_xxxxy, g_0_z_xxyzz_xxxxz, g_0_z_xxyzz_xxxyy, g_0_z_xxyzz_xxxyz, g_0_z_xxyzz_xxxzz, g_0_z_xxyzz_xxyyy, g_0_z_xxyzz_xxyyz, g_0_z_xxyzz_xxyzz, g_0_z_xxyzz_xxzzz, g_0_z_xxyzz_xyyyy, g_0_z_xxyzz_xyyyz, g_0_z_xxyzz_xyyzz, g_0_z_xxyzz_xyzzz, g_0_z_xxyzz_xzzzz, g_0_z_xxyzz_yyyyy, g_0_z_xxyzz_yyyyz, g_0_z_xxyzz_yyyzz, g_0_z_xxyzz_yyzzz, g_0_z_xxyzz_yzzzz, g_0_z_xxyzz_zzzzz, g_0_z_xyzz_xxxxx, g_0_z_xyzz_xxxxxx, g_0_z_xyzz_xxxxxy, g_0_z_xyzz_xxxxxz, g_0_z_xyzz_xxxxy, g_0_z_xyzz_xxxxyy, g_0_z_xyzz_xxxxyz, g_0_z_xyzz_xxxxz, g_0_z_xyzz_xxxxzz, g_0_z_xyzz_xxxyy, g_0_z_xyzz_xxxyyy, g_0_z_xyzz_xxxyyz, g_0_z_xyzz_xxxyz, g_0_z_xyzz_xxxyzz, g_0_z_xyzz_xxxzz, g_0_z_xyzz_xxxzzz, g_0_z_xyzz_xxyyy, g_0_z_xyzz_xxyyyy, g_0_z_xyzz_xxyyyz, g_0_z_xyzz_xxyyz, g_0_z_xyzz_xxyyzz, g_0_z_xyzz_xxyzz, g_0_z_xyzz_xxyzzz, g_0_z_xyzz_xxzzz, g_0_z_xyzz_xxzzzz, g_0_z_xyzz_xyyyy, g_0_z_xyzz_xyyyyy, g_0_z_xyzz_xyyyyz, g_0_z_xyzz_xyyyz, g_0_z_xyzz_xyyyzz, g_0_z_xyzz_xyyzz, g_0_z_xyzz_xyyzzz, g_0_z_xyzz_xyzzz, g_0_z_xyzz_xyzzzz, g_0_z_xyzz_xzzzz, g_0_z_xyzz_xzzzzz, g_0_z_xyzz_yyyyy, g_0_z_xyzz_yyyyz, g_0_z_xyzz_yyyzz, g_0_z_xyzz_yyzzz, g_0_z_xyzz_yzzzz, g_0_z_xyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzz_xxxxx[k] = -g_0_z_xyzz_xxxxx[k] * ab_x + g_0_z_xyzz_xxxxxx[k];

                g_0_z_xxyzz_xxxxy[k] = -g_0_z_xyzz_xxxxy[k] * ab_x + g_0_z_xyzz_xxxxxy[k];

                g_0_z_xxyzz_xxxxz[k] = -g_0_z_xyzz_xxxxz[k] * ab_x + g_0_z_xyzz_xxxxxz[k];

                g_0_z_xxyzz_xxxyy[k] = -g_0_z_xyzz_xxxyy[k] * ab_x + g_0_z_xyzz_xxxxyy[k];

                g_0_z_xxyzz_xxxyz[k] = -g_0_z_xyzz_xxxyz[k] * ab_x + g_0_z_xyzz_xxxxyz[k];

                g_0_z_xxyzz_xxxzz[k] = -g_0_z_xyzz_xxxzz[k] * ab_x + g_0_z_xyzz_xxxxzz[k];

                g_0_z_xxyzz_xxyyy[k] = -g_0_z_xyzz_xxyyy[k] * ab_x + g_0_z_xyzz_xxxyyy[k];

                g_0_z_xxyzz_xxyyz[k] = -g_0_z_xyzz_xxyyz[k] * ab_x + g_0_z_xyzz_xxxyyz[k];

                g_0_z_xxyzz_xxyzz[k] = -g_0_z_xyzz_xxyzz[k] * ab_x + g_0_z_xyzz_xxxyzz[k];

                g_0_z_xxyzz_xxzzz[k] = -g_0_z_xyzz_xxzzz[k] * ab_x + g_0_z_xyzz_xxxzzz[k];

                g_0_z_xxyzz_xyyyy[k] = -g_0_z_xyzz_xyyyy[k] * ab_x + g_0_z_xyzz_xxyyyy[k];

                g_0_z_xxyzz_xyyyz[k] = -g_0_z_xyzz_xyyyz[k] * ab_x + g_0_z_xyzz_xxyyyz[k];

                g_0_z_xxyzz_xyyzz[k] = -g_0_z_xyzz_xyyzz[k] * ab_x + g_0_z_xyzz_xxyyzz[k];

                g_0_z_xxyzz_xyzzz[k] = -g_0_z_xyzz_xyzzz[k] * ab_x + g_0_z_xyzz_xxyzzz[k];

                g_0_z_xxyzz_xzzzz[k] = -g_0_z_xyzz_xzzzz[k] * ab_x + g_0_z_xyzz_xxzzzz[k];

                g_0_z_xxyzz_yyyyy[k] = -g_0_z_xyzz_yyyyy[k] * ab_x + g_0_z_xyzz_xyyyyy[k];

                g_0_z_xxyzz_yyyyz[k] = -g_0_z_xyzz_yyyyz[k] * ab_x + g_0_z_xyzz_xyyyyz[k];

                g_0_z_xxyzz_yyyzz[k] = -g_0_z_xyzz_yyyzz[k] * ab_x + g_0_z_xyzz_xyyyzz[k];

                g_0_z_xxyzz_yyzzz[k] = -g_0_z_xyzz_yyzzz[k] * ab_x + g_0_z_xyzz_xyyzzz[k];

                g_0_z_xxyzz_yzzzz[k] = -g_0_z_xyzz_yzzzz[k] * ab_x + g_0_z_xyzz_xyzzzz[k];

                g_0_z_xxyzz_zzzzz[k] = -g_0_z_xyzz_zzzzz[k] * ab_x + g_0_z_xyzz_xzzzzz[k];
            }

            /// Set up 1071-1092 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzz_xxxxx = cbuffer.data(hh_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxy = cbuffer.data(hh_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxxz = cbuffer.data(hh_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxyy = cbuffer.data(hh_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxyz = cbuffer.data(hh_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxzz = cbuffer.data(hh_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxyyy = cbuffer.data(hh_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxyyz = cbuffer.data(hh_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxyzz = cbuffer.data(hh_geom_01_off + 1079 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxzzz = cbuffer.data(hh_geom_01_off + 1080 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyyyy = cbuffer.data(hh_geom_01_off + 1081 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyyyz = cbuffer.data(hh_geom_01_off + 1082 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyyzz = cbuffer.data(hh_geom_01_off + 1083 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyzzz = cbuffer.data(hh_geom_01_off + 1084 * ccomps * dcomps);

            auto g_0_z_xxzzz_xzzzz = cbuffer.data(hh_geom_01_off + 1085 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyyyy = cbuffer.data(hh_geom_01_off + 1086 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyyyz = cbuffer.data(hh_geom_01_off + 1087 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyyzz = cbuffer.data(hh_geom_01_off + 1088 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyzzz = cbuffer.data(hh_geom_01_off + 1089 * ccomps * dcomps);

            auto g_0_z_xxzzz_yzzzz = cbuffer.data(hh_geom_01_off + 1090 * ccomps * dcomps);

            auto g_0_z_xxzzz_zzzzz = cbuffer.data(hh_geom_01_off + 1091 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzz_xxxxx, g_0_z_xxzzz_xxxxy, g_0_z_xxzzz_xxxxz, g_0_z_xxzzz_xxxyy, g_0_z_xxzzz_xxxyz, g_0_z_xxzzz_xxxzz, g_0_z_xxzzz_xxyyy, g_0_z_xxzzz_xxyyz, g_0_z_xxzzz_xxyzz, g_0_z_xxzzz_xxzzz, g_0_z_xxzzz_xyyyy, g_0_z_xxzzz_xyyyz, g_0_z_xxzzz_xyyzz, g_0_z_xxzzz_xyzzz, g_0_z_xxzzz_xzzzz, g_0_z_xxzzz_yyyyy, g_0_z_xxzzz_yyyyz, g_0_z_xxzzz_yyyzz, g_0_z_xxzzz_yyzzz, g_0_z_xxzzz_yzzzz, g_0_z_xxzzz_zzzzz, g_0_z_xzzz_xxxxx, g_0_z_xzzz_xxxxxx, g_0_z_xzzz_xxxxxy, g_0_z_xzzz_xxxxxz, g_0_z_xzzz_xxxxy, g_0_z_xzzz_xxxxyy, g_0_z_xzzz_xxxxyz, g_0_z_xzzz_xxxxz, g_0_z_xzzz_xxxxzz, g_0_z_xzzz_xxxyy, g_0_z_xzzz_xxxyyy, g_0_z_xzzz_xxxyyz, g_0_z_xzzz_xxxyz, g_0_z_xzzz_xxxyzz, g_0_z_xzzz_xxxzz, g_0_z_xzzz_xxxzzz, g_0_z_xzzz_xxyyy, g_0_z_xzzz_xxyyyy, g_0_z_xzzz_xxyyyz, g_0_z_xzzz_xxyyz, g_0_z_xzzz_xxyyzz, g_0_z_xzzz_xxyzz, g_0_z_xzzz_xxyzzz, g_0_z_xzzz_xxzzz, g_0_z_xzzz_xxzzzz, g_0_z_xzzz_xyyyy, g_0_z_xzzz_xyyyyy, g_0_z_xzzz_xyyyyz, g_0_z_xzzz_xyyyz, g_0_z_xzzz_xyyyzz, g_0_z_xzzz_xyyzz, g_0_z_xzzz_xyyzzz, g_0_z_xzzz_xyzzz, g_0_z_xzzz_xyzzzz, g_0_z_xzzz_xzzzz, g_0_z_xzzz_xzzzzz, g_0_z_xzzz_yyyyy, g_0_z_xzzz_yyyyz, g_0_z_xzzz_yyyzz, g_0_z_xzzz_yyzzz, g_0_z_xzzz_yzzzz, g_0_z_xzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzz_xxxxx[k] = -g_0_z_xzzz_xxxxx[k] * ab_x + g_0_z_xzzz_xxxxxx[k];

                g_0_z_xxzzz_xxxxy[k] = -g_0_z_xzzz_xxxxy[k] * ab_x + g_0_z_xzzz_xxxxxy[k];

                g_0_z_xxzzz_xxxxz[k] = -g_0_z_xzzz_xxxxz[k] * ab_x + g_0_z_xzzz_xxxxxz[k];

                g_0_z_xxzzz_xxxyy[k] = -g_0_z_xzzz_xxxyy[k] * ab_x + g_0_z_xzzz_xxxxyy[k];

                g_0_z_xxzzz_xxxyz[k] = -g_0_z_xzzz_xxxyz[k] * ab_x + g_0_z_xzzz_xxxxyz[k];

                g_0_z_xxzzz_xxxzz[k] = -g_0_z_xzzz_xxxzz[k] * ab_x + g_0_z_xzzz_xxxxzz[k];

                g_0_z_xxzzz_xxyyy[k] = -g_0_z_xzzz_xxyyy[k] * ab_x + g_0_z_xzzz_xxxyyy[k];

                g_0_z_xxzzz_xxyyz[k] = -g_0_z_xzzz_xxyyz[k] * ab_x + g_0_z_xzzz_xxxyyz[k];

                g_0_z_xxzzz_xxyzz[k] = -g_0_z_xzzz_xxyzz[k] * ab_x + g_0_z_xzzz_xxxyzz[k];

                g_0_z_xxzzz_xxzzz[k] = -g_0_z_xzzz_xxzzz[k] * ab_x + g_0_z_xzzz_xxxzzz[k];

                g_0_z_xxzzz_xyyyy[k] = -g_0_z_xzzz_xyyyy[k] * ab_x + g_0_z_xzzz_xxyyyy[k];

                g_0_z_xxzzz_xyyyz[k] = -g_0_z_xzzz_xyyyz[k] * ab_x + g_0_z_xzzz_xxyyyz[k];

                g_0_z_xxzzz_xyyzz[k] = -g_0_z_xzzz_xyyzz[k] * ab_x + g_0_z_xzzz_xxyyzz[k];

                g_0_z_xxzzz_xyzzz[k] = -g_0_z_xzzz_xyzzz[k] * ab_x + g_0_z_xzzz_xxyzzz[k];

                g_0_z_xxzzz_xzzzz[k] = -g_0_z_xzzz_xzzzz[k] * ab_x + g_0_z_xzzz_xxzzzz[k];

                g_0_z_xxzzz_yyyyy[k] = -g_0_z_xzzz_yyyyy[k] * ab_x + g_0_z_xzzz_xyyyyy[k];

                g_0_z_xxzzz_yyyyz[k] = -g_0_z_xzzz_yyyyz[k] * ab_x + g_0_z_xzzz_xyyyyz[k];

                g_0_z_xxzzz_yyyzz[k] = -g_0_z_xzzz_yyyzz[k] * ab_x + g_0_z_xzzz_xyyyzz[k];

                g_0_z_xxzzz_yyzzz[k] = -g_0_z_xzzz_yyzzz[k] * ab_x + g_0_z_xzzz_xyyzzz[k];

                g_0_z_xxzzz_yzzzz[k] = -g_0_z_xzzz_yzzzz[k] * ab_x + g_0_z_xzzz_xyzzzz[k];

                g_0_z_xxzzz_zzzzz[k] = -g_0_z_xzzz_zzzzz[k] * ab_x + g_0_z_xzzz_xzzzzz[k];
            }

            /// Set up 1092-1113 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyy_xxxxx = cbuffer.data(hh_geom_01_off + 1092 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxy = cbuffer.data(hh_geom_01_off + 1093 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxxz = cbuffer.data(hh_geom_01_off + 1094 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxyy = cbuffer.data(hh_geom_01_off + 1095 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxyz = cbuffer.data(hh_geom_01_off + 1096 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxzz = cbuffer.data(hh_geom_01_off + 1097 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxyyy = cbuffer.data(hh_geom_01_off + 1098 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxyyz = cbuffer.data(hh_geom_01_off + 1099 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxyzz = cbuffer.data(hh_geom_01_off + 1100 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxzzz = cbuffer.data(hh_geom_01_off + 1101 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyyyy = cbuffer.data(hh_geom_01_off + 1102 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyyyz = cbuffer.data(hh_geom_01_off + 1103 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyyzz = cbuffer.data(hh_geom_01_off + 1104 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyzzz = cbuffer.data(hh_geom_01_off + 1105 * ccomps * dcomps);

            auto g_0_z_xyyyy_xzzzz = cbuffer.data(hh_geom_01_off + 1106 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyyyy = cbuffer.data(hh_geom_01_off + 1107 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyyyz = cbuffer.data(hh_geom_01_off + 1108 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyyzz = cbuffer.data(hh_geom_01_off + 1109 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyzzz = cbuffer.data(hh_geom_01_off + 1110 * ccomps * dcomps);

            auto g_0_z_xyyyy_yzzzz = cbuffer.data(hh_geom_01_off + 1111 * ccomps * dcomps);

            auto g_0_z_xyyyy_zzzzz = cbuffer.data(hh_geom_01_off + 1112 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyy_xxxxx, g_0_z_xyyyy_xxxxy, g_0_z_xyyyy_xxxxz, g_0_z_xyyyy_xxxyy, g_0_z_xyyyy_xxxyz, g_0_z_xyyyy_xxxzz, g_0_z_xyyyy_xxyyy, g_0_z_xyyyy_xxyyz, g_0_z_xyyyy_xxyzz, g_0_z_xyyyy_xxzzz, g_0_z_xyyyy_xyyyy, g_0_z_xyyyy_xyyyz, g_0_z_xyyyy_xyyzz, g_0_z_xyyyy_xyzzz, g_0_z_xyyyy_xzzzz, g_0_z_xyyyy_yyyyy, g_0_z_xyyyy_yyyyz, g_0_z_xyyyy_yyyzz, g_0_z_xyyyy_yyzzz, g_0_z_xyyyy_yzzzz, g_0_z_xyyyy_zzzzz, g_0_z_yyyy_xxxxx, g_0_z_yyyy_xxxxxx, g_0_z_yyyy_xxxxxy, g_0_z_yyyy_xxxxxz, g_0_z_yyyy_xxxxy, g_0_z_yyyy_xxxxyy, g_0_z_yyyy_xxxxyz, g_0_z_yyyy_xxxxz, g_0_z_yyyy_xxxxzz, g_0_z_yyyy_xxxyy, g_0_z_yyyy_xxxyyy, g_0_z_yyyy_xxxyyz, g_0_z_yyyy_xxxyz, g_0_z_yyyy_xxxyzz, g_0_z_yyyy_xxxzz, g_0_z_yyyy_xxxzzz, g_0_z_yyyy_xxyyy, g_0_z_yyyy_xxyyyy, g_0_z_yyyy_xxyyyz, g_0_z_yyyy_xxyyz, g_0_z_yyyy_xxyyzz, g_0_z_yyyy_xxyzz, g_0_z_yyyy_xxyzzz, g_0_z_yyyy_xxzzz, g_0_z_yyyy_xxzzzz, g_0_z_yyyy_xyyyy, g_0_z_yyyy_xyyyyy, g_0_z_yyyy_xyyyyz, g_0_z_yyyy_xyyyz, g_0_z_yyyy_xyyyzz, g_0_z_yyyy_xyyzz, g_0_z_yyyy_xyyzzz, g_0_z_yyyy_xyzzz, g_0_z_yyyy_xyzzzz, g_0_z_yyyy_xzzzz, g_0_z_yyyy_xzzzzz, g_0_z_yyyy_yyyyy, g_0_z_yyyy_yyyyz, g_0_z_yyyy_yyyzz, g_0_z_yyyy_yyzzz, g_0_z_yyyy_yzzzz, g_0_z_yyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyy_xxxxx[k] = -g_0_z_yyyy_xxxxx[k] * ab_x + g_0_z_yyyy_xxxxxx[k];

                g_0_z_xyyyy_xxxxy[k] = -g_0_z_yyyy_xxxxy[k] * ab_x + g_0_z_yyyy_xxxxxy[k];

                g_0_z_xyyyy_xxxxz[k] = -g_0_z_yyyy_xxxxz[k] * ab_x + g_0_z_yyyy_xxxxxz[k];

                g_0_z_xyyyy_xxxyy[k] = -g_0_z_yyyy_xxxyy[k] * ab_x + g_0_z_yyyy_xxxxyy[k];

                g_0_z_xyyyy_xxxyz[k] = -g_0_z_yyyy_xxxyz[k] * ab_x + g_0_z_yyyy_xxxxyz[k];

                g_0_z_xyyyy_xxxzz[k] = -g_0_z_yyyy_xxxzz[k] * ab_x + g_0_z_yyyy_xxxxzz[k];

                g_0_z_xyyyy_xxyyy[k] = -g_0_z_yyyy_xxyyy[k] * ab_x + g_0_z_yyyy_xxxyyy[k];

                g_0_z_xyyyy_xxyyz[k] = -g_0_z_yyyy_xxyyz[k] * ab_x + g_0_z_yyyy_xxxyyz[k];

                g_0_z_xyyyy_xxyzz[k] = -g_0_z_yyyy_xxyzz[k] * ab_x + g_0_z_yyyy_xxxyzz[k];

                g_0_z_xyyyy_xxzzz[k] = -g_0_z_yyyy_xxzzz[k] * ab_x + g_0_z_yyyy_xxxzzz[k];

                g_0_z_xyyyy_xyyyy[k] = -g_0_z_yyyy_xyyyy[k] * ab_x + g_0_z_yyyy_xxyyyy[k];

                g_0_z_xyyyy_xyyyz[k] = -g_0_z_yyyy_xyyyz[k] * ab_x + g_0_z_yyyy_xxyyyz[k];

                g_0_z_xyyyy_xyyzz[k] = -g_0_z_yyyy_xyyzz[k] * ab_x + g_0_z_yyyy_xxyyzz[k];

                g_0_z_xyyyy_xyzzz[k] = -g_0_z_yyyy_xyzzz[k] * ab_x + g_0_z_yyyy_xxyzzz[k];

                g_0_z_xyyyy_xzzzz[k] = -g_0_z_yyyy_xzzzz[k] * ab_x + g_0_z_yyyy_xxzzzz[k];

                g_0_z_xyyyy_yyyyy[k] = -g_0_z_yyyy_yyyyy[k] * ab_x + g_0_z_yyyy_xyyyyy[k];

                g_0_z_xyyyy_yyyyz[k] = -g_0_z_yyyy_yyyyz[k] * ab_x + g_0_z_yyyy_xyyyyz[k];

                g_0_z_xyyyy_yyyzz[k] = -g_0_z_yyyy_yyyzz[k] * ab_x + g_0_z_yyyy_xyyyzz[k];

                g_0_z_xyyyy_yyzzz[k] = -g_0_z_yyyy_yyzzz[k] * ab_x + g_0_z_yyyy_xyyzzz[k];

                g_0_z_xyyyy_yzzzz[k] = -g_0_z_yyyy_yzzzz[k] * ab_x + g_0_z_yyyy_xyzzzz[k];

                g_0_z_xyyyy_zzzzz[k] = -g_0_z_yyyy_zzzzz[k] * ab_x + g_0_z_yyyy_xzzzzz[k];
            }

            /// Set up 1113-1134 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyz_xxxxx = cbuffer.data(hh_geom_01_off + 1113 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxy = cbuffer.data(hh_geom_01_off + 1114 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxxz = cbuffer.data(hh_geom_01_off + 1115 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxyy = cbuffer.data(hh_geom_01_off + 1116 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxyz = cbuffer.data(hh_geom_01_off + 1117 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxzz = cbuffer.data(hh_geom_01_off + 1118 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxyyy = cbuffer.data(hh_geom_01_off + 1119 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxyyz = cbuffer.data(hh_geom_01_off + 1120 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxyzz = cbuffer.data(hh_geom_01_off + 1121 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxzzz = cbuffer.data(hh_geom_01_off + 1122 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyyyy = cbuffer.data(hh_geom_01_off + 1123 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyyyz = cbuffer.data(hh_geom_01_off + 1124 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyyzz = cbuffer.data(hh_geom_01_off + 1125 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyzzz = cbuffer.data(hh_geom_01_off + 1126 * ccomps * dcomps);

            auto g_0_z_xyyyz_xzzzz = cbuffer.data(hh_geom_01_off + 1127 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyyyy = cbuffer.data(hh_geom_01_off + 1128 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyyyz = cbuffer.data(hh_geom_01_off + 1129 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyyzz = cbuffer.data(hh_geom_01_off + 1130 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyzzz = cbuffer.data(hh_geom_01_off + 1131 * ccomps * dcomps);

            auto g_0_z_xyyyz_yzzzz = cbuffer.data(hh_geom_01_off + 1132 * ccomps * dcomps);

            auto g_0_z_xyyyz_zzzzz = cbuffer.data(hh_geom_01_off + 1133 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyz_xxxxx, g_0_z_xyyyz_xxxxy, g_0_z_xyyyz_xxxxz, g_0_z_xyyyz_xxxyy, g_0_z_xyyyz_xxxyz, g_0_z_xyyyz_xxxzz, g_0_z_xyyyz_xxyyy, g_0_z_xyyyz_xxyyz, g_0_z_xyyyz_xxyzz, g_0_z_xyyyz_xxzzz, g_0_z_xyyyz_xyyyy, g_0_z_xyyyz_xyyyz, g_0_z_xyyyz_xyyzz, g_0_z_xyyyz_xyzzz, g_0_z_xyyyz_xzzzz, g_0_z_xyyyz_yyyyy, g_0_z_xyyyz_yyyyz, g_0_z_xyyyz_yyyzz, g_0_z_xyyyz_yyzzz, g_0_z_xyyyz_yzzzz, g_0_z_xyyyz_zzzzz, g_0_z_yyyz_xxxxx, g_0_z_yyyz_xxxxxx, g_0_z_yyyz_xxxxxy, g_0_z_yyyz_xxxxxz, g_0_z_yyyz_xxxxy, g_0_z_yyyz_xxxxyy, g_0_z_yyyz_xxxxyz, g_0_z_yyyz_xxxxz, g_0_z_yyyz_xxxxzz, g_0_z_yyyz_xxxyy, g_0_z_yyyz_xxxyyy, g_0_z_yyyz_xxxyyz, g_0_z_yyyz_xxxyz, g_0_z_yyyz_xxxyzz, g_0_z_yyyz_xxxzz, g_0_z_yyyz_xxxzzz, g_0_z_yyyz_xxyyy, g_0_z_yyyz_xxyyyy, g_0_z_yyyz_xxyyyz, g_0_z_yyyz_xxyyz, g_0_z_yyyz_xxyyzz, g_0_z_yyyz_xxyzz, g_0_z_yyyz_xxyzzz, g_0_z_yyyz_xxzzz, g_0_z_yyyz_xxzzzz, g_0_z_yyyz_xyyyy, g_0_z_yyyz_xyyyyy, g_0_z_yyyz_xyyyyz, g_0_z_yyyz_xyyyz, g_0_z_yyyz_xyyyzz, g_0_z_yyyz_xyyzz, g_0_z_yyyz_xyyzzz, g_0_z_yyyz_xyzzz, g_0_z_yyyz_xyzzzz, g_0_z_yyyz_xzzzz, g_0_z_yyyz_xzzzzz, g_0_z_yyyz_yyyyy, g_0_z_yyyz_yyyyz, g_0_z_yyyz_yyyzz, g_0_z_yyyz_yyzzz, g_0_z_yyyz_yzzzz, g_0_z_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyz_xxxxx[k] = -g_0_z_yyyz_xxxxx[k] * ab_x + g_0_z_yyyz_xxxxxx[k];

                g_0_z_xyyyz_xxxxy[k] = -g_0_z_yyyz_xxxxy[k] * ab_x + g_0_z_yyyz_xxxxxy[k];

                g_0_z_xyyyz_xxxxz[k] = -g_0_z_yyyz_xxxxz[k] * ab_x + g_0_z_yyyz_xxxxxz[k];

                g_0_z_xyyyz_xxxyy[k] = -g_0_z_yyyz_xxxyy[k] * ab_x + g_0_z_yyyz_xxxxyy[k];

                g_0_z_xyyyz_xxxyz[k] = -g_0_z_yyyz_xxxyz[k] * ab_x + g_0_z_yyyz_xxxxyz[k];

                g_0_z_xyyyz_xxxzz[k] = -g_0_z_yyyz_xxxzz[k] * ab_x + g_0_z_yyyz_xxxxzz[k];

                g_0_z_xyyyz_xxyyy[k] = -g_0_z_yyyz_xxyyy[k] * ab_x + g_0_z_yyyz_xxxyyy[k];

                g_0_z_xyyyz_xxyyz[k] = -g_0_z_yyyz_xxyyz[k] * ab_x + g_0_z_yyyz_xxxyyz[k];

                g_0_z_xyyyz_xxyzz[k] = -g_0_z_yyyz_xxyzz[k] * ab_x + g_0_z_yyyz_xxxyzz[k];

                g_0_z_xyyyz_xxzzz[k] = -g_0_z_yyyz_xxzzz[k] * ab_x + g_0_z_yyyz_xxxzzz[k];

                g_0_z_xyyyz_xyyyy[k] = -g_0_z_yyyz_xyyyy[k] * ab_x + g_0_z_yyyz_xxyyyy[k];

                g_0_z_xyyyz_xyyyz[k] = -g_0_z_yyyz_xyyyz[k] * ab_x + g_0_z_yyyz_xxyyyz[k];

                g_0_z_xyyyz_xyyzz[k] = -g_0_z_yyyz_xyyzz[k] * ab_x + g_0_z_yyyz_xxyyzz[k];

                g_0_z_xyyyz_xyzzz[k] = -g_0_z_yyyz_xyzzz[k] * ab_x + g_0_z_yyyz_xxyzzz[k];

                g_0_z_xyyyz_xzzzz[k] = -g_0_z_yyyz_xzzzz[k] * ab_x + g_0_z_yyyz_xxzzzz[k];

                g_0_z_xyyyz_yyyyy[k] = -g_0_z_yyyz_yyyyy[k] * ab_x + g_0_z_yyyz_xyyyyy[k];

                g_0_z_xyyyz_yyyyz[k] = -g_0_z_yyyz_yyyyz[k] * ab_x + g_0_z_yyyz_xyyyyz[k];

                g_0_z_xyyyz_yyyzz[k] = -g_0_z_yyyz_yyyzz[k] * ab_x + g_0_z_yyyz_xyyyzz[k];

                g_0_z_xyyyz_yyzzz[k] = -g_0_z_yyyz_yyzzz[k] * ab_x + g_0_z_yyyz_xyyzzz[k];

                g_0_z_xyyyz_yzzzz[k] = -g_0_z_yyyz_yzzzz[k] * ab_x + g_0_z_yyyz_xyzzzz[k];

                g_0_z_xyyyz_zzzzz[k] = -g_0_z_yyyz_zzzzz[k] * ab_x + g_0_z_yyyz_xzzzzz[k];
            }

            /// Set up 1134-1155 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzz_xxxxx = cbuffer.data(hh_geom_01_off + 1134 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxy = cbuffer.data(hh_geom_01_off + 1135 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxxz = cbuffer.data(hh_geom_01_off + 1136 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxyy = cbuffer.data(hh_geom_01_off + 1137 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxyz = cbuffer.data(hh_geom_01_off + 1138 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxzz = cbuffer.data(hh_geom_01_off + 1139 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxyyy = cbuffer.data(hh_geom_01_off + 1140 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxyyz = cbuffer.data(hh_geom_01_off + 1141 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxyzz = cbuffer.data(hh_geom_01_off + 1142 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxzzz = cbuffer.data(hh_geom_01_off + 1143 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyyyy = cbuffer.data(hh_geom_01_off + 1144 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyyyz = cbuffer.data(hh_geom_01_off + 1145 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyyzz = cbuffer.data(hh_geom_01_off + 1146 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyzzz = cbuffer.data(hh_geom_01_off + 1147 * ccomps * dcomps);

            auto g_0_z_xyyzz_xzzzz = cbuffer.data(hh_geom_01_off + 1148 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyyyy = cbuffer.data(hh_geom_01_off + 1149 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyyyz = cbuffer.data(hh_geom_01_off + 1150 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyyzz = cbuffer.data(hh_geom_01_off + 1151 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyzzz = cbuffer.data(hh_geom_01_off + 1152 * ccomps * dcomps);

            auto g_0_z_xyyzz_yzzzz = cbuffer.data(hh_geom_01_off + 1153 * ccomps * dcomps);

            auto g_0_z_xyyzz_zzzzz = cbuffer.data(hh_geom_01_off + 1154 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzz_xxxxx, g_0_z_xyyzz_xxxxy, g_0_z_xyyzz_xxxxz, g_0_z_xyyzz_xxxyy, g_0_z_xyyzz_xxxyz, g_0_z_xyyzz_xxxzz, g_0_z_xyyzz_xxyyy, g_0_z_xyyzz_xxyyz, g_0_z_xyyzz_xxyzz, g_0_z_xyyzz_xxzzz, g_0_z_xyyzz_xyyyy, g_0_z_xyyzz_xyyyz, g_0_z_xyyzz_xyyzz, g_0_z_xyyzz_xyzzz, g_0_z_xyyzz_xzzzz, g_0_z_xyyzz_yyyyy, g_0_z_xyyzz_yyyyz, g_0_z_xyyzz_yyyzz, g_0_z_xyyzz_yyzzz, g_0_z_xyyzz_yzzzz, g_0_z_xyyzz_zzzzz, g_0_z_yyzz_xxxxx, g_0_z_yyzz_xxxxxx, g_0_z_yyzz_xxxxxy, g_0_z_yyzz_xxxxxz, g_0_z_yyzz_xxxxy, g_0_z_yyzz_xxxxyy, g_0_z_yyzz_xxxxyz, g_0_z_yyzz_xxxxz, g_0_z_yyzz_xxxxzz, g_0_z_yyzz_xxxyy, g_0_z_yyzz_xxxyyy, g_0_z_yyzz_xxxyyz, g_0_z_yyzz_xxxyz, g_0_z_yyzz_xxxyzz, g_0_z_yyzz_xxxzz, g_0_z_yyzz_xxxzzz, g_0_z_yyzz_xxyyy, g_0_z_yyzz_xxyyyy, g_0_z_yyzz_xxyyyz, g_0_z_yyzz_xxyyz, g_0_z_yyzz_xxyyzz, g_0_z_yyzz_xxyzz, g_0_z_yyzz_xxyzzz, g_0_z_yyzz_xxzzz, g_0_z_yyzz_xxzzzz, g_0_z_yyzz_xyyyy, g_0_z_yyzz_xyyyyy, g_0_z_yyzz_xyyyyz, g_0_z_yyzz_xyyyz, g_0_z_yyzz_xyyyzz, g_0_z_yyzz_xyyzz, g_0_z_yyzz_xyyzzz, g_0_z_yyzz_xyzzz, g_0_z_yyzz_xyzzzz, g_0_z_yyzz_xzzzz, g_0_z_yyzz_xzzzzz, g_0_z_yyzz_yyyyy, g_0_z_yyzz_yyyyz, g_0_z_yyzz_yyyzz, g_0_z_yyzz_yyzzz, g_0_z_yyzz_yzzzz, g_0_z_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzz_xxxxx[k] = -g_0_z_yyzz_xxxxx[k] * ab_x + g_0_z_yyzz_xxxxxx[k];

                g_0_z_xyyzz_xxxxy[k] = -g_0_z_yyzz_xxxxy[k] * ab_x + g_0_z_yyzz_xxxxxy[k];

                g_0_z_xyyzz_xxxxz[k] = -g_0_z_yyzz_xxxxz[k] * ab_x + g_0_z_yyzz_xxxxxz[k];

                g_0_z_xyyzz_xxxyy[k] = -g_0_z_yyzz_xxxyy[k] * ab_x + g_0_z_yyzz_xxxxyy[k];

                g_0_z_xyyzz_xxxyz[k] = -g_0_z_yyzz_xxxyz[k] * ab_x + g_0_z_yyzz_xxxxyz[k];

                g_0_z_xyyzz_xxxzz[k] = -g_0_z_yyzz_xxxzz[k] * ab_x + g_0_z_yyzz_xxxxzz[k];

                g_0_z_xyyzz_xxyyy[k] = -g_0_z_yyzz_xxyyy[k] * ab_x + g_0_z_yyzz_xxxyyy[k];

                g_0_z_xyyzz_xxyyz[k] = -g_0_z_yyzz_xxyyz[k] * ab_x + g_0_z_yyzz_xxxyyz[k];

                g_0_z_xyyzz_xxyzz[k] = -g_0_z_yyzz_xxyzz[k] * ab_x + g_0_z_yyzz_xxxyzz[k];

                g_0_z_xyyzz_xxzzz[k] = -g_0_z_yyzz_xxzzz[k] * ab_x + g_0_z_yyzz_xxxzzz[k];

                g_0_z_xyyzz_xyyyy[k] = -g_0_z_yyzz_xyyyy[k] * ab_x + g_0_z_yyzz_xxyyyy[k];

                g_0_z_xyyzz_xyyyz[k] = -g_0_z_yyzz_xyyyz[k] * ab_x + g_0_z_yyzz_xxyyyz[k];

                g_0_z_xyyzz_xyyzz[k] = -g_0_z_yyzz_xyyzz[k] * ab_x + g_0_z_yyzz_xxyyzz[k];

                g_0_z_xyyzz_xyzzz[k] = -g_0_z_yyzz_xyzzz[k] * ab_x + g_0_z_yyzz_xxyzzz[k];

                g_0_z_xyyzz_xzzzz[k] = -g_0_z_yyzz_xzzzz[k] * ab_x + g_0_z_yyzz_xxzzzz[k];

                g_0_z_xyyzz_yyyyy[k] = -g_0_z_yyzz_yyyyy[k] * ab_x + g_0_z_yyzz_xyyyyy[k];

                g_0_z_xyyzz_yyyyz[k] = -g_0_z_yyzz_yyyyz[k] * ab_x + g_0_z_yyzz_xyyyyz[k];

                g_0_z_xyyzz_yyyzz[k] = -g_0_z_yyzz_yyyzz[k] * ab_x + g_0_z_yyzz_xyyyzz[k];

                g_0_z_xyyzz_yyzzz[k] = -g_0_z_yyzz_yyzzz[k] * ab_x + g_0_z_yyzz_xyyzzz[k];

                g_0_z_xyyzz_yzzzz[k] = -g_0_z_yyzz_yzzzz[k] * ab_x + g_0_z_yyzz_xyzzzz[k];

                g_0_z_xyyzz_zzzzz[k] = -g_0_z_yyzz_zzzzz[k] * ab_x + g_0_z_yyzz_xzzzzz[k];
            }

            /// Set up 1155-1176 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzz_xxxxx = cbuffer.data(hh_geom_01_off + 1155 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxy = cbuffer.data(hh_geom_01_off + 1156 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxxz = cbuffer.data(hh_geom_01_off + 1157 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxyy = cbuffer.data(hh_geom_01_off + 1158 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxyz = cbuffer.data(hh_geom_01_off + 1159 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxzz = cbuffer.data(hh_geom_01_off + 1160 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxyyy = cbuffer.data(hh_geom_01_off + 1161 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxyyz = cbuffer.data(hh_geom_01_off + 1162 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxyzz = cbuffer.data(hh_geom_01_off + 1163 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxzzz = cbuffer.data(hh_geom_01_off + 1164 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyyyy = cbuffer.data(hh_geom_01_off + 1165 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyyyz = cbuffer.data(hh_geom_01_off + 1166 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyyzz = cbuffer.data(hh_geom_01_off + 1167 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyzzz = cbuffer.data(hh_geom_01_off + 1168 * ccomps * dcomps);

            auto g_0_z_xyzzz_xzzzz = cbuffer.data(hh_geom_01_off + 1169 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyyyy = cbuffer.data(hh_geom_01_off + 1170 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyyyz = cbuffer.data(hh_geom_01_off + 1171 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyyzz = cbuffer.data(hh_geom_01_off + 1172 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyzzz = cbuffer.data(hh_geom_01_off + 1173 * ccomps * dcomps);

            auto g_0_z_xyzzz_yzzzz = cbuffer.data(hh_geom_01_off + 1174 * ccomps * dcomps);

            auto g_0_z_xyzzz_zzzzz = cbuffer.data(hh_geom_01_off + 1175 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzz_xxxxx, g_0_z_xyzzz_xxxxy, g_0_z_xyzzz_xxxxz, g_0_z_xyzzz_xxxyy, g_0_z_xyzzz_xxxyz, g_0_z_xyzzz_xxxzz, g_0_z_xyzzz_xxyyy, g_0_z_xyzzz_xxyyz, g_0_z_xyzzz_xxyzz, g_0_z_xyzzz_xxzzz, g_0_z_xyzzz_xyyyy, g_0_z_xyzzz_xyyyz, g_0_z_xyzzz_xyyzz, g_0_z_xyzzz_xyzzz, g_0_z_xyzzz_xzzzz, g_0_z_xyzzz_yyyyy, g_0_z_xyzzz_yyyyz, g_0_z_xyzzz_yyyzz, g_0_z_xyzzz_yyzzz, g_0_z_xyzzz_yzzzz, g_0_z_xyzzz_zzzzz, g_0_z_yzzz_xxxxx, g_0_z_yzzz_xxxxxx, g_0_z_yzzz_xxxxxy, g_0_z_yzzz_xxxxxz, g_0_z_yzzz_xxxxy, g_0_z_yzzz_xxxxyy, g_0_z_yzzz_xxxxyz, g_0_z_yzzz_xxxxz, g_0_z_yzzz_xxxxzz, g_0_z_yzzz_xxxyy, g_0_z_yzzz_xxxyyy, g_0_z_yzzz_xxxyyz, g_0_z_yzzz_xxxyz, g_0_z_yzzz_xxxyzz, g_0_z_yzzz_xxxzz, g_0_z_yzzz_xxxzzz, g_0_z_yzzz_xxyyy, g_0_z_yzzz_xxyyyy, g_0_z_yzzz_xxyyyz, g_0_z_yzzz_xxyyz, g_0_z_yzzz_xxyyzz, g_0_z_yzzz_xxyzz, g_0_z_yzzz_xxyzzz, g_0_z_yzzz_xxzzz, g_0_z_yzzz_xxzzzz, g_0_z_yzzz_xyyyy, g_0_z_yzzz_xyyyyy, g_0_z_yzzz_xyyyyz, g_0_z_yzzz_xyyyz, g_0_z_yzzz_xyyyzz, g_0_z_yzzz_xyyzz, g_0_z_yzzz_xyyzzz, g_0_z_yzzz_xyzzz, g_0_z_yzzz_xyzzzz, g_0_z_yzzz_xzzzz, g_0_z_yzzz_xzzzzz, g_0_z_yzzz_yyyyy, g_0_z_yzzz_yyyyz, g_0_z_yzzz_yyyzz, g_0_z_yzzz_yyzzz, g_0_z_yzzz_yzzzz, g_0_z_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzz_xxxxx[k] = -g_0_z_yzzz_xxxxx[k] * ab_x + g_0_z_yzzz_xxxxxx[k];

                g_0_z_xyzzz_xxxxy[k] = -g_0_z_yzzz_xxxxy[k] * ab_x + g_0_z_yzzz_xxxxxy[k];

                g_0_z_xyzzz_xxxxz[k] = -g_0_z_yzzz_xxxxz[k] * ab_x + g_0_z_yzzz_xxxxxz[k];

                g_0_z_xyzzz_xxxyy[k] = -g_0_z_yzzz_xxxyy[k] * ab_x + g_0_z_yzzz_xxxxyy[k];

                g_0_z_xyzzz_xxxyz[k] = -g_0_z_yzzz_xxxyz[k] * ab_x + g_0_z_yzzz_xxxxyz[k];

                g_0_z_xyzzz_xxxzz[k] = -g_0_z_yzzz_xxxzz[k] * ab_x + g_0_z_yzzz_xxxxzz[k];

                g_0_z_xyzzz_xxyyy[k] = -g_0_z_yzzz_xxyyy[k] * ab_x + g_0_z_yzzz_xxxyyy[k];

                g_0_z_xyzzz_xxyyz[k] = -g_0_z_yzzz_xxyyz[k] * ab_x + g_0_z_yzzz_xxxyyz[k];

                g_0_z_xyzzz_xxyzz[k] = -g_0_z_yzzz_xxyzz[k] * ab_x + g_0_z_yzzz_xxxyzz[k];

                g_0_z_xyzzz_xxzzz[k] = -g_0_z_yzzz_xxzzz[k] * ab_x + g_0_z_yzzz_xxxzzz[k];

                g_0_z_xyzzz_xyyyy[k] = -g_0_z_yzzz_xyyyy[k] * ab_x + g_0_z_yzzz_xxyyyy[k];

                g_0_z_xyzzz_xyyyz[k] = -g_0_z_yzzz_xyyyz[k] * ab_x + g_0_z_yzzz_xxyyyz[k];

                g_0_z_xyzzz_xyyzz[k] = -g_0_z_yzzz_xyyzz[k] * ab_x + g_0_z_yzzz_xxyyzz[k];

                g_0_z_xyzzz_xyzzz[k] = -g_0_z_yzzz_xyzzz[k] * ab_x + g_0_z_yzzz_xxyzzz[k];

                g_0_z_xyzzz_xzzzz[k] = -g_0_z_yzzz_xzzzz[k] * ab_x + g_0_z_yzzz_xxzzzz[k];

                g_0_z_xyzzz_yyyyy[k] = -g_0_z_yzzz_yyyyy[k] * ab_x + g_0_z_yzzz_xyyyyy[k];

                g_0_z_xyzzz_yyyyz[k] = -g_0_z_yzzz_yyyyz[k] * ab_x + g_0_z_yzzz_xyyyyz[k];

                g_0_z_xyzzz_yyyzz[k] = -g_0_z_yzzz_yyyzz[k] * ab_x + g_0_z_yzzz_xyyyzz[k];

                g_0_z_xyzzz_yyzzz[k] = -g_0_z_yzzz_yyzzz[k] * ab_x + g_0_z_yzzz_xyyzzz[k];

                g_0_z_xyzzz_yzzzz[k] = -g_0_z_yzzz_yzzzz[k] * ab_x + g_0_z_yzzz_xyzzzz[k];

                g_0_z_xyzzz_zzzzz[k] = -g_0_z_yzzz_zzzzz[k] * ab_x + g_0_z_yzzz_xzzzzz[k];
            }

            /// Set up 1176-1197 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzz_xxxxx = cbuffer.data(hh_geom_01_off + 1176 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxy = cbuffer.data(hh_geom_01_off + 1177 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxxz = cbuffer.data(hh_geom_01_off + 1178 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxyy = cbuffer.data(hh_geom_01_off + 1179 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxyz = cbuffer.data(hh_geom_01_off + 1180 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxzz = cbuffer.data(hh_geom_01_off + 1181 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxyyy = cbuffer.data(hh_geom_01_off + 1182 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxyyz = cbuffer.data(hh_geom_01_off + 1183 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxyzz = cbuffer.data(hh_geom_01_off + 1184 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxzzz = cbuffer.data(hh_geom_01_off + 1185 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyyyy = cbuffer.data(hh_geom_01_off + 1186 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyyyz = cbuffer.data(hh_geom_01_off + 1187 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyyzz = cbuffer.data(hh_geom_01_off + 1188 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyzzz = cbuffer.data(hh_geom_01_off + 1189 * ccomps * dcomps);

            auto g_0_z_xzzzz_xzzzz = cbuffer.data(hh_geom_01_off + 1190 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyyyy = cbuffer.data(hh_geom_01_off + 1191 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyyyz = cbuffer.data(hh_geom_01_off + 1192 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyyzz = cbuffer.data(hh_geom_01_off + 1193 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyzzz = cbuffer.data(hh_geom_01_off + 1194 * ccomps * dcomps);

            auto g_0_z_xzzzz_yzzzz = cbuffer.data(hh_geom_01_off + 1195 * ccomps * dcomps);

            auto g_0_z_xzzzz_zzzzz = cbuffer.data(hh_geom_01_off + 1196 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzz_xxxxx, g_0_z_xzzzz_xxxxy, g_0_z_xzzzz_xxxxz, g_0_z_xzzzz_xxxyy, g_0_z_xzzzz_xxxyz, g_0_z_xzzzz_xxxzz, g_0_z_xzzzz_xxyyy, g_0_z_xzzzz_xxyyz, g_0_z_xzzzz_xxyzz, g_0_z_xzzzz_xxzzz, g_0_z_xzzzz_xyyyy, g_0_z_xzzzz_xyyyz, g_0_z_xzzzz_xyyzz, g_0_z_xzzzz_xyzzz, g_0_z_xzzzz_xzzzz, g_0_z_xzzzz_yyyyy, g_0_z_xzzzz_yyyyz, g_0_z_xzzzz_yyyzz, g_0_z_xzzzz_yyzzz, g_0_z_xzzzz_yzzzz, g_0_z_xzzzz_zzzzz, g_0_z_zzzz_xxxxx, g_0_z_zzzz_xxxxxx, g_0_z_zzzz_xxxxxy, g_0_z_zzzz_xxxxxz, g_0_z_zzzz_xxxxy, g_0_z_zzzz_xxxxyy, g_0_z_zzzz_xxxxyz, g_0_z_zzzz_xxxxz, g_0_z_zzzz_xxxxzz, g_0_z_zzzz_xxxyy, g_0_z_zzzz_xxxyyy, g_0_z_zzzz_xxxyyz, g_0_z_zzzz_xxxyz, g_0_z_zzzz_xxxyzz, g_0_z_zzzz_xxxzz, g_0_z_zzzz_xxxzzz, g_0_z_zzzz_xxyyy, g_0_z_zzzz_xxyyyy, g_0_z_zzzz_xxyyyz, g_0_z_zzzz_xxyyz, g_0_z_zzzz_xxyyzz, g_0_z_zzzz_xxyzz, g_0_z_zzzz_xxyzzz, g_0_z_zzzz_xxzzz, g_0_z_zzzz_xxzzzz, g_0_z_zzzz_xyyyy, g_0_z_zzzz_xyyyyy, g_0_z_zzzz_xyyyyz, g_0_z_zzzz_xyyyz, g_0_z_zzzz_xyyyzz, g_0_z_zzzz_xyyzz, g_0_z_zzzz_xyyzzz, g_0_z_zzzz_xyzzz, g_0_z_zzzz_xyzzzz, g_0_z_zzzz_xzzzz, g_0_z_zzzz_xzzzzz, g_0_z_zzzz_yyyyy, g_0_z_zzzz_yyyyz, g_0_z_zzzz_yyyzz, g_0_z_zzzz_yyzzz, g_0_z_zzzz_yzzzz, g_0_z_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzz_xxxxx[k] = -g_0_z_zzzz_xxxxx[k] * ab_x + g_0_z_zzzz_xxxxxx[k];

                g_0_z_xzzzz_xxxxy[k] = -g_0_z_zzzz_xxxxy[k] * ab_x + g_0_z_zzzz_xxxxxy[k];

                g_0_z_xzzzz_xxxxz[k] = -g_0_z_zzzz_xxxxz[k] * ab_x + g_0_z_zzzz_xxxxxz[k];

                g_0_z_xzzzz_xxxyy[k] = -g_0_z_zzzz_xxxyy[k] * ab_x + g_0_z_zzzz_xxxxyy[k];

                g_0_z_xzzzz_xxxyz[k] = -g_0_z_zzzz_xxxyz[k] * ab_x + g_0_z_zzzz_xxxxyz[k];

                g_0_z_xzzzz_xxxzz[k] = -g_0_z_zzzz_xxxzz[k] * ab_x + g_0_z_zzzz_xxxxzz[k];

                g_0_z_xzzzz_xxyyy[k] = -g_0_z_zzzz_xxyyy[k] * ab_x + g_0_z_zzzz_xxxyyy[k];

                g_0_z_xzzzz_xxyyz[k] = -g_0_z_zzzz_xxyyz[k] * ab_x + g_0_z_zzzz_xxxyyz[k];

                g_0_z_xzzzz_xxyzz[k] = -g_0_z_zzzz_xxyzz[k] * ab_x + g_0_z_zzzz_xxxyzz[k];

                g_0_z_xzzzz_xxzzz[k] = -g_0_z_zzzz_xxzzz[k] * ab_x + g_0_z_zzzz_xxxzzz[k];

                g_0_z_xzzzz_xyyyy[k] = -g_0_z_zzzz_xyyyy[k] * ab_x + g_0_z_zzzz_xxyyyy[k];

                g_0_z_xzzzz_xyyyz[k] = -g_0_z_zzzz_xyyyz[k] * ab_x + g_0_z_zzzz_xxyyyz[k];

                g_0_z_xzzzz_xyyzz[k] = -g_0_z_zzzz_xyyzz[k] * ab_x + g_0_z_zzzz_xxyyzz[k];

                g_0_z_xzzzz_xyzzz[k] = -g_0_z_zzzz_xyzzz[k] * ab_x + g_0_z_zzzz_xxyzzz[k];

                g_0_z_xzzzz_xzzzz[k] = -g_0_z_zzzz_xzzzz[k] * ab_x + g_0_z_zzzz_xxzzzz[k];

                g_0_z_xzzzz_yyyyy[k] = -g_0_z_zzzz_yyyyy[k] * ab_x + g_0_z_zzzz_xyyyyy[k];

                g_0_z_xzzzz_yyyyz[k] = -g_0_z_zzzz_yyyyz[k] * ab_x + g_0_z_zzzz_xyyyyz[k];

                g_0_z_xzzzz_yyyzz[k] = -g_0_z_zzzz_yyyzz[k] * ab_x + g_0_z_zzzz_xyyyzz[k];

                g_0_z_xzzzz_yyzzz[k] = -g_0_z_zzzz_yyzzz[k] * ab_x + g_0_z_zzzz_xyyzzz[k];

                g_0_z_xzzzz_yzzzz[k] = -g_0_z_zzzz_yzzzz[k] * ab_x + g_0_z_zzzz_xyzzzz[k];

                g_0_z_xzzzz_zzzzz[k] = -g_0_z_zzzz_zzzzz[k] * ab_x + g_0_z_zzzz_xzzzzz[k];
            }

            /// Set up 1197-1218 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyy_xxxxx = cbuffer.data(hh_geom_01_off + 1197 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxy = cbuffer.data(hh_geom_01_off + 1198 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxxz = cbuffer.data(hh_geom_01_off + 1199 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxyy = cbuffer.data(hh_geom_01_off + 1200 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxyz = cbuffer.data(hh_geom_01_off + 1201 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxzz = cbuffer.data(hh_geom_01_off + 1202 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxyyy = cbuffer.data(hh_geom_01_off + 1203 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxyyz = cbuffer.data(hh_geom_01_off + 1204 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxyzz = cbuffer.data(hh_geom_01_off + 1205 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxzzz = cbuffer.data(hh_geom_01_off + 1206 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyyyy = cbuffer.data(hh_geom_01_off + 1207 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyyyz = cbuffer.data(hh_geom_01_off + 1208 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyyzz = cbuffer.data(hh_geom_01_off + 1209 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyzzz = cbuffer.data(hh_geom_01_off + 1210 * ccomps * dcomps);

            auto g_0_z_yyyyy_xzzzz = cbuffer.data(hh_geom_01_off + 1211 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyyyy = cbuffer.data(hh_geom_01_off + 1212 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyyyz = cbuffer.data(hh_geom_01_off + 1213 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyyzz = cbuffer.data(hh_geom_01_off + 1214 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyzzz = cbuffer.data(hh_geom_01_off + 1215 * ccomps * dcomps);

            auto g_0_z_yyyyy_yzzzz = cbuffer.data(hh_geom_01_off + 1216 * ccomps * dcomps);

            auto g_0_z_yyyyy_zzzzz = cbuffer.data(hh_geom_01_off + 1217 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyy_xxxxx, g_0_z_yyyy_xxxxxy, g_0_z_yyyy_xxxxy, g_0_z_yyyy_xxxxyy, g_0_z_yyyy_xxxxyz, g_0_z_yyyy_xxxxz, g_0_z_yyyy_xxxyy, g_0_z_yyyy_xxxyyy, g_0_z_yyyy_xxxyyz, g_0_z_yyyy_xxxyz, g_0_z_yyyy_xxxyzz, g_0_z_yyyy_xxxzz, g_0_z_yyyy_xxyyy, g_0_z_yyyy_xxyyyy, g_0_z_yyyy_xxyyyz, g_0_z_yyyy_xxyyz, g_0_z_yyyy_xxyyzz, g_0_z_yyyy_xxyzz, g_0_z_yyyy_xxyzzz, g_0_z_yyyy_xxzzz, g_0_z_yyyy_xyyyy, g_0_z_yyyy_xyyyyy, g_0_z_yyyy_xyyyyz, g_0_z_yyyy_xyyyz, g_0_z_yyyy_xyyyzz, g_0_z_yyyy_xyyzz, g_0_z_yyyy_xyyzzz, g_0_z_yyyy_xyzzz, g_0_z_yyyy_xyzzzz, g_0_z_yyyy_xzzzz, g_0_z_yyyy_yyyyy, g_0_z_yyyy_yyyyyy, g_0_z_yyyy_yyyyyz, g_0_z_yyyy_yyyyz, g_0_z_yyyy_yyyyzz, g_0_z_yyyy_yyyzz, g_0_z_yyyy_yyyzzz, g_0_z_yyyy_yyzzz, g_0_z_yyyy_yyzzzz, g_0_z_yyyy_yzzzz, g_0_z_yyyy_yzzzzz, g_0_z_yyyy_zzzzz, g_0_z_yyyyy_xxxxx, g_0_z_yyyyy_xxxxy, g_0_z_yyyyy_xxxxz, g_0_z_yyyyy_xxxyy, g_0_z_yyyyy_xxxyz, g_0_z_yyyyy_xxxzz, g_0_z_yyyyy_xxyyy, g_0_z_yyyyy_xxyyz, g_0_z_yyyyy_xxyzz, g_0_z_yyyyy_xxzzz, g_0_z_yyyyy_xyyyy, g_0_z_yyyyy_xyyyz, g_0_z_yyyyy_xyyzz, g_0_z_yyyyy_xyzzz, g_0_z_yyyyy_xzzzz, g_0_z_yyyyy_yyyyy, g_0_z_yyyyy_yyyyz, g_0_z_yyyyy_yyyzz, g_0_z_yyyyy_yyzzz, g_0_z_yyyyy_yzzzz, g_0_z_yyyyy_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyy_xxxxx[k] = -g_0_z_yyyy_xxxxx[k] * ab_y + g_0_z_yyyy_xxxxxy[k];

                g_0_z_yyyyy_xxxxy[k] = -g_0_z_yyyy_xxxxy[k] * ab_y + g_0_z_yyyy_xxxxyy[k];

                g_0_z_yyyyy_xxxxz[k] = -g_0_z_yyyy_xxxxz[k] * ab_y + g_0_z_yyyy_xxxxyz[k];

                g_0_z_yyyyy_xxxyy[k] = -g_0_z_yyyy_xxxyy[k] * ab_y + g_0_z_yyyy_xxxyyy[k];

                g_0_z_yyyyy_xxxyz[k] = -g_0_z_yyyy_xxxyz[k] * ab_y + g_0_z_yyyy_xxxyyz[k];

                g_0_z_yyyyy_xxxzz[k] = -g_0_z_yyyy_xxxzz[k] * ab_y + g_0_z_yyyy_xxxyzz[k];

                g_0_z_yyyyy_xxyyy[k] = -g_0_z_yyyy_xxyyy[k] * ab_y + g_0_z_yyyy_xxyyyy[k];

                g_0_z_yyyyy_xxyyz[k] = -g_0_z_yyyy_xxyyz[k] * ab_y + g_0_z_yyyy_xxyyyz[k];

                g_0_z_yyyyy_xxyzz[k] = -g_0_z_yyyy_xxyzz[k] * ab_y + g_0_z_yyyy_xxyyzz[k];

                g_0_z_yyyyy_xxzzz[k] = -g_0_z_yyyy_xxzzz[k] * ab_y + g_0_z_yyyy_xxyzzz[k];

                g_0_z_yyyyy_xyyyy[k] = -g_0_z_yyyy_xyyyy[k] * ab_y + g_0_z_yyyy_xyyyyy[k];

                g_0_z_yyyyy_xyyyz[k] = -g_0_z_yyyy_xyyyz[k] * ab_y + g_0_z_yyyy_xyyyyz[k];

                g_0_z_yyyyy_xyyzz[k] = -g_0_z_yyyy_xyyzz[k] * ab_y + g_0_z_yyyy_xyyyzz[k];

                g_0_z_yyyyy_xyzzz[k] = -g_0_z_yyyy_xyzzz[k] * ab_y + g_0_z_yyyy_xyyzzz[k];

                g_0_z_yyyyy_xzzzz[k] = -g_0_z_yyyy_xzzzz[k] * ab_y + g_0_z_yyyy_xyzzzz[k];

                g_0_z_yyyyy_yyyyy[k] = -g_0_z_yyyy_yyyyy[k] * ab_y + g_0_z_yyyy_yyyyyy[k];

                g_0_z_yyyyy_yyyyz[k] = -g_0_z_yyyy_yyyyz[k] * ab_y + g_0_z_yyyy_yyyyyz[k];

                g_0_z_yyyyy_yyyzz[k] = -g_0_z_yyyy_yyyzz[k] * ab_y + g_0_z_yyyy_yyyyzz[k];

                g_0_z_yyyyy_yyzzz[k] = -g_0_z_yyyy_yyzzz[k] * ab_y + g_0_z_yyyy_yyyzzz[k];

                g_0_z_yyyyy_yzzzz[k] = -g_0_z_yyyy_yzzzz[k] * ab_y + g_0_z_yyyy_yyzzzz[k];

                g_0_z_yyyyy_zzzzz[k] = -g_0_z_yyyy_zzzzz[k] * ab_y + g_0_z_yyyy_yzzzzz[k];
            }

            /// Set up 1218-1239 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyz_xxxxx = cbuffer.data(hh_geom_01_off + 1218 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxy = cbuffer.data(hh_geom_01_off + 1219 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxxz = cbuffer.data(hh_geom_01_off + 1220 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxyy = cbuffer.data(hh_geom_01_off + 1221 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxyz = cbuffer.data(hh_geom_01_off + 1222 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxzz = cbuffer.data(hh_geom_01_off + 1223 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxyyy = cbuffer.data(hh_geom_01_off + 1224 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxyyz = cbuffer.data(hh_geom_01_off + 1225 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxyzz = cbuffer.data(hh_geom_01_off + 1226 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxzzz = cbuffer.data(hh_geom_01_off + 1227 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyyyy = cbuffer.data(hh_geom_01_off + 1228 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyyyz = cbuffer.data(hh_geom_01_off + 1229 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyyzz = cbuffer.data(hh_geom_01_off + 1230 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyzzz = cbuffer.data(hh_geom_01_off + 1231 * ccomps * dcomps);

            auto g_0_z_yyyyz_xzzzz = cbuffer.data(hh_geom_01_off + 1232 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyyyy = cbuffer.data(hh_geom_01_off + 1233 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyyyz = cbuffer.data(hh_geom_01_off + 1234 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyyzz = cbuffer.data(hh_geom_01_off + 1235 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyzzz = cbuffer.data(hh_geom_01_off + 1236 * ccomps * dcomps);

            auto g_0_z_yyyyz_yzzzz = cbuffer.data(hh_geom_01_off + 1237 * ccomps * dcomps);

            auto g_0_z_yyyyz_zzzzz = cbuffer.data(hh_geom_01_off + 1238 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyz_xxxxx, g_0_z_yyyyz_xxxxy, g_0_z_yyyyz_xxxxz, g_0_z_yyyyz_xxxyy, g_0_z_yyyyz_xxxyz, g_0_z_yyyyz_xxxzz, g_0_z_yyyyz_xxyyy, g_0_z_yyyyz_xxyyz, g_0_z_yyyyz_xxyzz, g_0_z_yyyyz_xxzzz, g_0_z_yyyyz_xyyyy, g_0_z_yyyyz_xyyyz, g_0_z_yyyyz_xyyzz, g_0_z_yyyyz_xyzzz, g_0_z_yyyyz_xzzzz, g_0_z_yyyyz_yyyyy, g_0_z_yyyyz_yyyyz, g_0_z_yyyyz_yyyzz, g_0_z_yyyyz_yyzzz, g_0_z_yyyyz_yzzzz, g_0_z_yyyyz_zzzzz, g_0_z_yyyz_xxxxx, g_0_z_yyyz_xxxxxy, g_0_z_yyyz_xxxxy, g_0_z_yyyz_xxxxyy, g_0_z_yyyz_xxxxyz, g_0_z_yyyz_xxxxz, g_0_z_yyyz_xxxyy, g_0_z_yyyz_xxxyyy, g_0_z_yyyz_xxxyyz, g_0_z_yyyz_xxxyz, g_0_z_yyyz_xxxyzz, g_0_z_yyyz_xxxzz, g_0_z_yyyz_xxyyy, g_0_z_yyyz_xxyyyy, g_0_z_yyyz_xxyyyz, g_0_z_yyyz_xxyyz, g_0_z_yyyz_xxyyzz, g_0_z_yyyz_xxyzz, g_0_z_yyyz_xxyzzz, g_0_z_yyyz_xxzzz, g_0_z_yyyz_xyyyy, g_0_z_yyyz_xyyyyy, g_0_z_yyyz_xyyyyz, g_0_z_yyyz_xyyyz, g_0_z_yyyz_xyyyzz, g_0_z_yyyz_xyyzz, g_0_z_yyyz_xyyzzz, g_0_z_yyyz_xyzzz, g_0_z_yyyz_xyzzzz, g_0_z_yyyz_xzzzz, g_0_z_yyyz_yyyyy, g_0_z_yyyz_yyyyyy, g_0_z_yyyz_yyyyyz, g_0_z_yyyz_yyyyz, g_0_z_yyyz_yyyyzz, g_0_z_yyyz_yyyzz, g_0_z_yyyz_yyyzzz, g_0_z_yyyz_yyzzz, g_0_z_yyyz_yyzzzz, g_0_z_yyyz_yzzzz, g_0_z_yyyz_yzzzzz, g_0_z_yyyz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyz_xxxxx[k] = -g_0_z_yyyz_xxxxx[k] * ab_y + g_0_z_yyyz_xxxxxy[k];

                g_0_z_yyyyz_xxxxy[k] = -g_0_z_yyyz_xxxxy[k] * ab_y + g_0_z_yyyz_xxxxyy[k];

                g_0_z_yyyyz_xxxxz[k] = -g_0_z_yyyz_xxxxz[k] * ab_y + g_0_z_yyyz_xxxxyz[k];

                g_0_z_yyyyz_xxxyy[k] = -g_0_z_yyyz_xxxyy[k] * ab_y + g_0_z_yyyz_xxxyyy[k];

                g_0_z_yyyyz_xxxyz[k] = -g_0_z_yyyz_xxxyz[k] * ab_y + g_0_z_yyyz_xxxyyz[k];

                g_0_z_yyyyz_xxxzz[k] = -g_0_z_yyyz_xxxzz[k] * ab_y + g_0_z_yyyz_xxxyzz[k];

                g_0_z_yyyyz_xxyyy[k] = -g_0_z_yyyz_xxyyy[k] * ab_y + g_0_z_yyyz_xxyyyy[k];

                g_0_z_yyyyz_xxyyz[k] = -g_0_z_yyyz_xxyyz[k] * ab_y + g_0_z_yyyz_xxyyyz[k];

                g_0_z_yyyyz_xxyzz[k] = -g_0_z_yyyz_xxyzz[k] * ab_y + g_0_z_yyyz_xxyyzz[k];

                g_0_z_yyyyz_xxzzz[k] = -g_0_z_yyyz_xxzzz[k] * ab_y + g_0_z_yyyz_xxyzzz[k];

                g_0_z_yyyyz_xyyyy[k] = -g_0_z_yyyz_xyyyy[k] * ab_y + g_0_z_yyyz_xyyyyy[k];

                g_0_z_yyyyz_xyyyz[k] = -g_0_z_yyyz_xyyyz[k] * ab_y + g_0_z_yyyz_xyyyyz[k];

                g_0_z_yyyyz_xyyzz[k] = -g_0_z_yyyz_xyyzz[k] * ab_y + g_0_z_yyyz_xyyyzz[k];

                g_0_z_yyyyz_xyzzz[k] = -g_0_z_yyyz_xyzzz[k] * ab_y + g_0_z_yyyz_xyyzzz[k];

                g_0_z_yyyyz_xzzzz[k] = -g_0_z_yyyz_xzzzz[k] * ab_y + g_0_z_yyyz_xyzzzz[k];

                g_0_z_yyyyz_yyyyy[k] = -g_0_z_yyyz_yyyyy[k] * ab_y + g_0_z_yyyz_yyyyyy[k];

                g_0_z_yyyyz_yyyyz[k] = -g_0_z_yyyz_yyyyz[k] * ab_y + g_0_z_yyyz_yyyyyz[k];

                g_0_z_yyyyz_yyyzz[k] = -g_0_z_yyyz_yyyzz[k] * ab_y + g_0_z_yyyz_yyyyzz[k];

                g_0_z_yyyyz_yyzzz[k] = -g_0_z_yyyz_yyzzz[k] * ab_y + g_0_z_yyyz_yyyzzz[k];

                g_0_z_yyyyz_yzzzz[k] = -g_0_z_yyyz_yzzzz[k] * ab_y + g_0_z_yyyz_yyzzzz[k];

                g_0_z_yyyyz_zzzzz[k] = -g_0_z_yyyz_zzzzz[k] * ab_y + g_0_z_yyyz_yzzzzz[k];
            }

            /// Set up 1239-1260 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzz_xxxxx = cbuffer.data(hh_geom_01_off + 1239 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxy = cbuffer.data(hh_geom_01_off + 1240 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxxz = cbuffer.data(hh_geom_01_off + 1241 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxyy = cbuffer.data(hh_geom_01_off + 1242 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxyz = cbuffer.data(hh_geom_01_off + 1243 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxzz = cbuffer.data(hh_geom_01_off + 1244 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxyyy = cbuffer.data(hh_geom_01_off + 1245 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxyyz = cbuffer.data(hh_geom_01_off + 1246 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxyzz = cbuffer.data(hh_geom_01_off + 1247 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxzzz = cbuffer.data(hh_geom_01_off + 1248 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyyyy = cbuffer.data(hh_geom_01_off + 1249 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyyyz = cbuffer.data(hh_geom_01_off + 1250 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyyzz = cbuffer.data(hh_geom_01_off + 1251 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyzzz = cbuffer.data(hh_geom_01_off + 1252 * ccomps * dcomps);

            auto g_0_z_yyyzz_xzzzz = cbuffer.data(hh_geom_01_off + 1253 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyyyy = cbuffer.data(hh_geom_01_off + 1254 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyyyz = cbuffer.data(hh_geom_01_off + 1255 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyyzz = cbuffer.data(hh_geom_01_off + 1256 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyzzz = cbuffer.data(hh_geom_01_off + 1257 * ccomps * dcomps);

            auto g_0_z_yyyzz_yzzzz = cbuffer.data(hh_geom_01_off + 1258 * ccomps * dcomps);

            auto g_0_z_yyyzz_zzzzz = cbuffer.data(hh_geom_01_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzz_xxxxx, g_0_z_yyyzz_xxxxy, g_0_z_yyyzz_xxxxz, g_0_z_yyyzz_xxxyy, g_0_z_yyyzz_xxxyz, g_0_z_yyyzz_xxxzz, g_0_z_yyyzz_xxyyy, g_0_z_yyyzz_xxyyz, g_0_z_yyyzz_xxyzz, g_0_z_yyyzz_xxzzz, g_0_z_yyyzz_xyyyy, g_0_z_yyyzz_xyyyz, g_0_z_yyyzz_xyyzz, g_0_z_yyyzz_xyzzz, g_0_z_yyyzz_xzzzz, g_0_z_yyyzz_yyyyy, g_0_z_yyyzz_yyyyz, g_0_z_yyyzz_yyyzz, g_0_z_yyyzz_yyzzz, g_0_z_yyyzz_yzzzz, g_0_z_yyyzz_zzzzz, g_0_z_yyzz_xxxxx, g_0_z_yyzz_xxxxxy, g_0_z_yyzz_xxxxy, g_0_z_yyzz_xxxxyy, g_0_z_yyzz_xxxxyz, g_0_z_yyzz_xxxxz, g_0_z_yyzz_xxxyy, g_0_z_yyzz_xxxyyy, g_0_z_yyzz_xxxyyz, g_0_z_yyzz_xxxyz, g_0_z_yyzz_xxxyzz, g_0_z_yyzz_xxxzz, g_0_z_yyzz_xxyyy, g_0_z_yyzz_xxyyyy, g_0_z_yyzz_xxyyyz, g_0_z_yyzz_xxyyz, g_0_z_yyzz_xxyyzz, g_0_z_yyzz_xxyzz, g_0_z_yyzz_xxyzzz, g_0_z_yyzz_xxzzz, g_0_z_yyzz_xyyyy, g_0_z_yyzz_xyyyyy, g_0_z_yyzz_xyyyyz, g_0_z_yyzz_xyyyz, g_0_z_yyzz_xyyyzz, g_0_z_yyzz_xyyzz, g_0_z_yyzz_xyyzzz, g_0_z_yyzz_xyzzz, g_0_z_yyzz_xyzzzz, g_0_z_yyzz_xzzzz, g_0_z_yyzz_yyyyy, g_0_z_yyzz_yyyyyy, g_0_z_yyzz_yyyyyz, g_0_z_yyzz_yyyyz, g_0_z_yyzz_yyyyzz, g_0_z_yyzz_yyyzz, g_0_z_yyzz_yyyzzz, g_0_z_yyzz_yyzzz, g_0_z_yyzz_yyzzzz, g_0_z_yyzz_yzzzz, g_0_z_yyzz_yzzzzz, g_0_z_yyzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzz_xxxxx[k] = -g_0_z_yyzz_xxxxx[k] * ab_y + g_0_z_yyzz_xxxxxy[k];

                g_0_z_yyyzz_xxxxy[k] = -g_0_z_yyzz_xxxxy[k] * ab_y + g_0_z_yyzz_xxxxyy[k];

                g_0_z_yyyzz_xxxxz[k] = -g_0_z_yyzz_xxxxz[k] * ab_y + g_0_z_yyzz_xxxxyz[k];

                g_0_z_yyyzz_xxxyy[k] = -g_0_z_yyzz_xxxyy[k] * ab_y + g_0_z_yyzz_xxxyyy[k];

                g_0_z_yyyzz_xxxyz[k] = -g_0_z_yyzz_xxxyz[k] * ab_y + g_0_z_yyzz_xxxyyz[k];

                g_0_z_yyyzz_xxxzz[k] = -g_0_z_yyzz_xxxzz[k] * ab_y + g_0_z_yyzz_xxxyzz[k];

                g_0_z_yyyzz_xxyyy[k] = -g_0_z_yyzz_xxyyy[k] * ab_y + g_0_z_yyzz_xxyyyy[k];

                g_0_z_yyyzz_xxyyz[k] = -g_0_z_yyzz_xxyyz[k] * ab_y + g_0_z_yyzz_xxyyyz[k];

                g_0_z_yyyzz_xxyzz[k] = -g_0_z_yyzz_xxyzz[k] * ab_y + g_0_z_yyzz_xxyyzz[k];

                g_0_z_yyyzz_xxzzz[k] = -g_0_z_yyzz_xxzzz[k] * ab_y + g_0_z_yyzz_xxyzzz[k];

                g_0_z_yyyzz_xyyyy[k] = -g_0_z_yyzz_xyyyy[k] * ab_y + g_0_z_yyzz_xyyyyy[k];

                g_0_z_yyyzz_xyyyz[k] = -g_0_z_yyzz_xyyyz[k] * ab_y + g_0_z_yyzz_xyyyyz[k];

                g_0_z_yyyzz_xyyzz[k] = -g_0_z_yyzz_xyyzz[k] * ab_y + g_0_z_yyzz_xyyyzz[k];

                g_0_z_yyyzz_xyzzz[k] = -g_0_z_yyzz_xyzzz[k] * ab_y + g_0_z_yyzz_xyyzzz[k];

                g_0_z_yyyzz_xzzzz[k] = -g_0_z_yyzz_xzzzz[k] * ab_y + g_0_z_yyzz_xyzzzz[k];

                g_0_z_yyyzz_yyyyy[k] = -g_0_z_yyzz_yyyyy[k] * ab_y + g_0_z_yyzz_yyyyyy[k];

                g_0_z_yyyzz_yyyyz[k] = -g_0_z_yyzz_yyyyz[k] * ab_y + g_0_z_yyzz_yyyyyz[k];

                g_0_z_yyyzz_yyyzz[k] = -g_0_z_yyzz_yyyzz[k] * ab_y + g_0_z_yyzz_yyyyzz[k];

                g_0_z_yyyzz_yyzzz[k] = -g_0_z_yyzz_yyzzz[k] * ab_y + g_0_z_yyzz_yyyzzz[k];

                g_0_z_yyyzz_yzzzz[k] = -g_0_z_yyzz_yzzzz[k] * ab_y + g_0_z_yyzz_yyzzzz[k];

                g_0_z_yyyzz_zzzzz[k] = -g_0_z_yyzz_zzzzz[k] * ab_y + g_0_z_yyzz_yzzzzz[k];
            }

            /// Set up 1260-1281 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzz_xxxxx = cbuffer.data(hh_geom_01_off + 1260 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxy = cbuffer.data(hh_geom_01_off + 1261 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxxz = cbuffer.data(hh_geom_01_off + 1262 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxyy = cbuffer.data(hh_geom_01_off + 1263 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxyz = cbuffer.data(hh_geom_01_off + 1264 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxzz = cbuffer.data(hh_geom_01_off + 1265 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxyyy = cbuffer.data(hh_geom_01_off + 1266 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxyyz = cbuffer.data(hh_geom_01_off + 1267 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxyzz = cbuffer.data(hh_geom_01_off + 1268 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxzzz = cbuffer.data(hh_geom_01_off + 1269 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyyyy = cbuffer.data(hh_geom_01_off + 1270 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyyyz = cbuffer.data(hh_geom_01_off + 1271 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyyzz = cbuffer.data(hh_geom_01_off + 1272 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyzzz = cbuffer.data(hh_geom_01_off + 1273 * ccomps * dcomps);

            auto g_0_z_yyzzz_xzzzz = cbuffer.data(hh_geom_01_off + 1274 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyyyy = cbuffer.data(hh_geom_01_off + 1275 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyyyz = cbuffer.data(hh_geom_01_off + 1276 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyyzz = cbuffer.data(hh_geom_01_off + 1277 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyzzz = cbuffer.data(hh_geom_01_off + 1278 * ccomps * dcomps);

            auto g_0_z_yyzzz_yzzzz = cbuffer.data(hh_geom_01_off + 1279 * ccomps * dcomps);

            auto g_0_z_yyzzz_zzzzz = cbuffer.data(hh_geom_01_off + 1280 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzz_xxxxx, g_0_z_yyzzz_xxxxy, g_0_z_yyzzz_xxxxz, g_0_z_yyzzz_xxxyy, g_0_z_yyzzz_xxxyz, g_0_z_yyzzz_xxxzz, g_0_z_yyzzz_xxyyy, g_0_z_yyzzz_xxyyz, g_0_z_yyzzz_xxyzz, g_0_z_yyzzz_xxzzz, g_0_z_yyzzz_xyyyy, g_0_z_yyzzz_xyyyz, g_0_z_yyzzz_xyyzz, g_0_z_yyzzz_xyzzz, g_0_z_yyzzz_xzzzz, g_0_z_yyzzz_yyyyy, g_0_z_yyzzz_yyyyz, g_0_z_yyzzz_yyyzz, g_0_z_yyzzz_yyzzz, g_0_z_yyzzz_yzzzz, g_0_z_yyzzz_zzzzz, g_0_z_yzzz_xxxxx, g_0_z_yzzz_xxxxxy, g_0_z_yzzz_xxxxy, g_0_z_yzzz_xxxxyy, g_0_z_yzzz_xxxxyz, g_0_z_yzzz_xxxxz, g_0_z_yzzz_xxxyy, g_0_z_yzzz_xxxyyy, g_0_z_yzzz_xxxyyz, g_0_z_yzzz_xxxyz, g_0_z_yzzz_xxxyzz, g_0_z_yzzz_xxxzz, g_0_z_yzzz_xxyyy, g_0_z_yzzz_xxyyyy, g_0_z_yzzz_xxyyyz, g_0_z_yzzz_xxyyz, g_0_z_yzzz_xxyyzz, g_0_z_yzzz_xxyzz, g_0_z_yzzz_xxyzzz, g_0_z_yzzz_xxzzz, g_0_z_yzzz_xyyyy, g_0_z_yzzz_xyyyyy, g_0_z_yzzz_xyyyyz, g_0_z_yzzz_xyyyz, g_0_z_yzzz_xyyyzz, g_0_z_yzzz_xyyzz, g_0_z_yzzz_xyyzzz, g_0_z_yzzz_xyzzz, g_0_z_yzzz_xyzzzz, g_0_z_yzzz_xzzzz, g_0_z_yzzz_yyyyy, g_0_z_yzzz_yyyyyy, g_0_z_yzzz_yyyyyz, g_0_z_yzzz_yyyyz, g_0_z_yzzz_yyyyzz, g_0_z_yzzz_yyyzz, g_0_z_yzzz_yyyzzz, g_0_z_yzzz_yyzzz, g_0_z_yzzz_yyzzzz, g_0_z_yzzz_yzzzz, g_0_z_yzzz_yzzzzz, g_0_z_yzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzz_xxxxx[k] = -g_0_z_yzzz_xxxxx[k] * ab_y + g_0_z_yzzz_xxxxxy[k];

                g_0_z_yyzzz_xxxxy[k] = -g_0_z_yzzz_xxxxy[k] * ab_y + g_0_z_yzzz_xxxxyy[k];

                g_0_z_yyzzz_xxxxz[k] = -g_0_z_yzzz_xxxxz[k] * ab_y + g_0_z_yzzz_xxxxyz[k];

                g_0_z_yyzzz_xxxyy[k] = -g_0_z_yzzz_xxxyy[k] * ab_y + g_0_z_yzzz_xxxyyy[k];

                g_0_z_yyzzz_xxxyz[k] = -g_0_z_yzzz_xxxyz[k] * ab_y + g_0_z_yzzz_xxxyyz[k];

                g_0_z_yyzzz_xxxzz[k] = -g_0_z_yzzz_xxxzz[k] * ab_y + g_0_z_yzzz_xxxyzz[k];

                g_0_z_yyzzz_xxyyy[k] = -g_0_z_yzzz_xxyyy[k] * ab_y + g_0_z_yzzz_xxyyyy[k];

                g_0_z_yyzzz_xxyyz[k] = -g_0_z_yzzz_xxyyz[k] * ab_y + g_0_z_yzzz_xxyyyz[k];

                g_0_z_yyzzz_xxyzz[k] = -g_0_z_yzzz_xxyzz[k] * ab_y + g_0_z_yzzz_xxyyzz[k];

                g_0_z_yyzzz_xxzzz[k] = -g_0_z_yzzz_xxzzz[k] * ab_y + g_0_z_yzzz_xxyzzz[k];

                g_0_z_yyzzz_xyyyy[k] = -g_0_z_yzzz_xyyyy[k] * ab_y + g_0_z_yzzz_xyyyyy[k];

                g_0_z_yyzzz_xyyyz[k] = -g_0_z_yzzz_xyyyz[k] * ab_y + g_0_z_yzzz_xyyyyz[k];

                g_0_z_yyzzz_xyyzz[k] = -g_0_z_yzzz_xyyzz[k] * ab_y + g_0_z_yzzz_xyyyzz[k];

                g_0_z_yyzzz_xyzzz[k] = -g_0_z_yzzz_xyzzz[k] * ab_y + g_0_z_yzzz_xyyzzz[k];

                g_0_z_yyzzz_xzzzz[k] = -g_0_z_yzzz_xzzzz[k] * ab_y + g_0_z_yzzz_xyzzzz[k];

                g_0_z_yyzzz_yyyyy[k] = -g_0_z_yzzz_yyyyy[k] * ab_y + g_0_z_yzzz_yyyyyy[k];

                g_0_z_yyzzz_yyyyz[k] = -g_0_z_yzzz_yyyyz[k] * ab_y + g_0_z_yzzz_yyyyyz[k];

                g_0_z_yyzzz_yyyzz[k] = -g_0_z_yzzz_yyyzz[k] * ab_y + g_0_z_yzzz_yyyyzz[k];

                g_0_z_yyzzz_yyzzz[k] = -g_0_z_yzzz_yyzzz[k] * ab_y + g_0_z_yzzz_yyyzzz[k];

                g_0_z_yyzzz_yzzzz[k] = -g_0_z_yzzz_yzzzz[k] * ab_y + g_0_z_yzzz_yyzzzz[k];

                g_0_z_yyzzz_zzzzz[k] = -g_0_z_yzzz_zzzzz[k] * ab_y + g_0_z_yzzz_yzzzzz[k];
            }

            /// Set up 1281-1302 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzz_xxxxx = cbuffer.data(hh_geom_01_off + 1281 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxy = cbuffer.data(hh_geom_01_off + 1282 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxxz = cbuffer.data(hh_geom_01_off + 1283 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxyy = cbuffer.data(hh_geom_01_off + 1284 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxyz = cbuffer.data(hh_geom_01_off + 1285 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxzz = cbuffer.data(hh_geom_01_off + 1286 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxyyy = cbuffer.data(hh_geom_01_off + 1287 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxyyz = cbuffer.data(hh_geom_01_off + 1288 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxyzz = cbuffer.data(hh_geom_01_off + 1289 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxzzz = cbuffer.data(hh_geom_01_off + 1290 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyyyy = cbuffer.data(hh_geom_01_off + 1291 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyyyz = cbuffer.data(hh_geom_01_off + 1292 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyyzz = cbuffer.data(hh_geom_01_off + 1293 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyzzz = cbuffer.data(hh_geom_01_off + 1294 * ccomps * dcomps);

            auto g_0_z_yzzzz_xzzzz = cbuffer.data(hh_geom_01_off + 1295 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyyyy = cbuffer.data(hh_geom_01_off + 1296 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyyyz = cbuffer.data(hh_geom_01_off + 1297 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyyzz = cbuffer.data(hh_geom_01_off + 1298 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyzzz = cbuffer.data(hh_geom_01_off + 1299 * ccomps * dcomps);

            auto g_0_z_yzzzz_yzzzz = cbuffer.data(hh_geom_01_off + 1300 * ccomps * dcomps);

            auto g_0_z_yzzzz_zzzzz = cbuffer.data(hh_geom_01_off + 1301 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzz_xxxxx, g_0_z_yzzzz_xxxxy, g_0_z_yzzzz_xxxxz, g_0_z_yzzzz_xxxyy, g_0_z_yzzzz_xxxyz, g_0_z_yzzzz_xxxzz, g_0_z_yzzzz_xxyyy, g_0_z_yzzzz_xxyyz, g_0_z_yzzzz_xxyzz, g_0_z_yzzzz_xxzzz, g_0_z_yzzzz_xyyyy, g_0_z_yzzzz_xyyyz, g_0_z_yzzzz_xyyzz, g_0_z_yzzzz_xyzzz, g_0_z_yzzzz_xzzzz, g_0_z_yzzzz_yyyyy, g_0_z_yzzzz_yyyyz, g_0_z_yzzzz_yyyzz, g_0_z_yzzzz_yyzzz, g_0_z_yzzzz_yzzzz, g_0_z_yzzzz_zzzzz, g_0_z_zzzz_xxxxx, g_0_z_zzzz_xxxxxy, g_0_z_zzzz_xxxxy, g_0_z_zzzz_xxxxyy, g_0_z_zzzz_xxxxyz, g_0_z_zzzz_xxxxz, g_0_z_zzzz_xxxyy, g_0_z_zzzz_xxxyyy, g_0_z_zzzz_xxxyyz, g_0_z_zzzz_xxxyz, g_0_z_zzzz_xxxyzz, g_0_z_zzzz_xxxzz, g_0_z_zzzz_xxyyy, g_0_z_zzzz_xxyyyy, g_0_z_zzzz_xxyyyz, g_0_z_zzzz_xxyyz, g_0_z_zzzz_xxyyzz, g_0_z_zzzz_xxyzz, g_0_z_zzzz_xxyzzz, g_0_z_zzzz_xxzzz, g_0_z_zzzz_xyyyy, g_0_z_zzzz_xyyyyy, g_0_z_zzzz_xyyyyz, g_0_z_zzzz_xyyyz, g_0_z_zzzz_xyyyzz, g_0_z_zzzz_xyyzz, g_0_z_zzzz_xyyzzz, g_0_z_zzzz_xyzzz, g_0_z_zzzz_xyzzzz, g_0_z_zzzz_xzzzz, g_0_z_zzzz_yyyyy, g_0_z_zzzz_yyyyyy, g_0_z_zzzz_yyyyyz, g_0_z_zzzz_yyyyz, g_0_z_zzzz_yyyyzz, g_0_z_zzzz_yyyzz, g_0_z_zzzz_yyyzzz, g_0_z_zzzz_yyzzz, g_0_z_zzzz_yyzzzz, g_0_z_zzzz_yzzzz, g_0_z_zzzz_yzzzzz, g_0_z_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzz_xxxxx[k] = -g_0_z_zzzz_xxxxx[k] * ab_y + g_0_z_zzzz_xxxxxy[k];

                g_0_z_yzzzz_xxxxy[k] = -g_0_z_zzzz_xxxxy[k] * ab_y + g_0_z_zzzz_xxxxyy[k];

                g_0_z_yzzzz_xxxxz[k] = -g_0_z_zzzz_xxxxz[k] * ab_y + g_0_z_zzzz_xxxxyz[k];

                g_0_z_yzzzz_xxxyy[k] = -g_0_z_zzzz_xxxyy[k] * ab_y + g_0_z_zzzz_xxxyyy[k];

                g_0_z_yzzzz_xxxyz[k] = -g_0_z_zzzz_xxxyz[k] * ab_y + g_0_z_zzzz_xxxyyz[k];

                g_0_z_yzzzz_xxxzz[k] = -g_0_z_zzzz_xxxzz[k] * ab_y + g_0_z_zzzz_xxxyzz[k];

                g_0_z_yzzzz_xxyyy[k] = -g_0_z_zzzz_xxyyy[k] * ab_y + g_0_z_zzzz_xxyyyy[k];

                g_0_z_yzzzz_xxyyz[k] = -g_0_z_zzzz_xxyyz[k] * ab_y + g_0_z_zzzz_xxyyyz[k];

                g_0_z_yzzzz_xxyzz[k] = -g_0_z_zzzz_xxyzz[k] * ab_y + g_0_z_zzzz_xxyyzz[k];

                g_0_z_yzzzz_xxzzz[k] = -g_0_z_zzzz_xxzzz[k] * ab_y + g_0_z_zzzz_xxyzzz[k];

                g_0_z_yzzzz_xyyyy[k] = -g_0_z_zzzz_xyyyy[k] * ab_y + g_0_z_zzzz_xyyyyy[k];

                g_0_z_yzzzz_xyyyz[k] = -g_0_z_zzzz_xyyyz[k] * ab_y + g_0_z_zzzz_xyyyyz[k];

                g_0_z_yzzzz_xyyzz[k] = -g_0_z_zzzz_xyyzz[k] * ab_y + g_0_z_zzzz_xyyyzz[k];

                g_0_z_yzzzz_xyzzz[k] = -g_0_z_zzzz_xyzzz[k] * ab_y + g_0_z_zzzz_xyyzzz[k];

                g_0_z_yzzzz_xzzzz[k] = -g_0_z_zzzz_xzzzz[k] * ab_y + g_0_z_zzzz_xyzzzz[k];

                g_0_z_yzzzz_yyyyy[k] = -g_0_z_zzzz_yyyyy[k] * ab_y + g_0_z_zzzz_yyyyyy[k];

                g_0_z_yzzzz_yyyyz[k] = -g_0_z_zzzz_yyyyz[k] * ab_y + g_0_z_zzzz_yyyyyz[k];

                g_0_z_yzzzz_yyyzz[k] = -g_0_z_zzzz_yyyzz[k] * ab_y + g_0_z_zzzz_yyyyzz[k];

                g_0_z_yzzzz_yyzzz[k] = -g_0_z_zzzz_yyzzz[k] * ab_y + g_0_z_zzzz_yyyzzz[k];

                g_0_z_yzzzz_yzzzz[k] = -g_0_z_zzzz_yzzzz[k] * ab_y + g_0_z_zzzz_yyzzzz[k];

                g_0_z_yzzzz_zzzzz[k] = -g_0_z_zzzz_zzzzz[k] * ab_y + g_0_z_zzzz_yzzzzz[k];
            }

            /// Set up 1302-1323 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzz_xxxxx = cbuffer.data(hh_geom_01_off + 1302 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxy = cbuffer.data(hh_geom_01_off + 1303 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxxz = cbuffer.data(hh_geom_01_off + 1304 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxyy = cbuffer.data(hh_geom_01_off + 1305 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxyz = cbuffer.data(hh_geom_01_off + 1306 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxzz = cbuffer.data(hh_geom_01_off + 1307 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxyyy = cbuffer.data(hh_geom_01_off + 1308 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxyyz = cbuffer.data(hh_geom_01_off + 1309 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxyzz = cbuffer.data(hh_geom_01_off + 1310 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxzzz = cbuffer.data(hh_geom_01_off + 1311 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyyyy = cbuffer.data(hh_geom_01_off + 1312 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyyyz = cbuffer.data(hh_geom_01_off + 1313 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyyzz = cbuffer.data(hh_geom_01_off + 1314 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyzzz = cbuffer.data(hh_geom_01_off + 1315 * ccomps * dcomps);

            auto g_0_z_zzzzz_xzzzz = cbuffer.data(hh_geom_01_off + 1316 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyyyy = cbuffer.data(hh_geom_01_off + 1317 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyyyz = cbuffer.data(hh_geom_01_off + 1318 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyyzz = cbuffer.data(hh_geom_01_off + 1319 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyzzz = cbuffer.data(hh_geom_01_off + 1320 * ccomps * dcomps);

            auto g_0_z_zzzzz_yzzzz = cbuffer.data(hh_geom_01_off + 1321 * ccomps * dcomps);

            auto g_0_z_zzzzz_zzzzz = cbuffer.data(hh_geom_01_off + 1322 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzz_xxxxx, g_0_z_zzzz_xxxxxz, g_0_z_zzzz_xxxxy, g_0_z_zzzz_xxxxyz, g_0_z_zzzz_xxxxz, g_0_z_zzzz_xxxxzz, g_0_z_zzzz_xxxyy, g_0_z_zzzz_xxxyyz, g_0_z_zzzz_xxxyz, g_0_z_zzzz_xxxyzz, g_0_z_zzzz_xxxzz, g_0_z_zzzz_xxxzzz, g_0_z_zzzz_xxyyy, g_0_z_zzzz_xxyyyz, g_0_z_zzzz_xxyyz, g_0_z_zzzz_xxyyzz, g_0_z_zzzz_xxyzz, g_0_z_zzzz_xxyzzz, g_0_z_zzzz_xxzzz, g_0_z_zzzz_xxzzzz, g_0_z_zzzz_xyyyy, g_0_z_zzzz_xyyyyz, g_0_z_zzzz_xyyyz, g_0_z_zzzz_xyyyzz, g_0_z_zzzz_xyyzz, g_0_z_zzzz_xyyzzz, g_0_z_zzzz_xyzzz, g_0_z_zzzz_xyzzzz, g_0_z_zzzz_xzzzz, g_0_z_zzzz_xzzzzz, g_0_z_zzzz_yyyyy, g_0_z_zzzz_yyyyyz, g_0_z_zzzz_yyyyz, g_0_z_zzzz_yyyyzz, g_0_z_zzzz_yyyzz, g_0_z_zzzz_yyyzzz, g_0_z_zzzz_yyzzz, g_0_z_zzzz_yyzzzz, g_0_z_zzzz_yzzzz, g_0_z_zzzz_yzzzzz, g_0_z_zzzz_zzzzz, g_0_z_zzzz_zzzzzz, g_0_z_zzzzz_xxxxx, g_0_z_zzzzz_xxxxy, g_0_z_zzzzz_xxxxz, g_0_z_zzzzz_xxxyy, g_0_z_zzzzz_xxxyz, g_0_z_zzzzz_xxxzz, g_0_z_zzzzz_xxyyy, g_0_z_zzzzz_xxyyz, g_0_z_zzzzz_xxyzz, g_0_z_zzzzz_xxzzz, g_0_z_zzzzz_xyyyy, g_0_z_zzzzz_xyyyz, g_0_z_zzzzz_xyyzz, g_0_z_zzzzz_xyzzz, g_0_z_zzzzz_xzzzz, g_0_z_zzzzz_yyyyy, g_0_z_zzzzz_yyyyz, g_0_z_zzzzz_yyyzz, g_0_z_zzzzz_yyzzz, g_0_z_zzzzz_yzzzz, g_0_z_zzzzz_zzzzz, g_zzzz_xxxxx, g_zzzz_xxxxy, g_zzzz_xxxxz, g_zzzz_xxxyy, g_zzzz_xxxyz, g_zzzz_xxxzz, g_zzzz_xxyyy, g_zzzz_xxyyz, g_zzzz_xxyzz, g_zzzz_xxzzz, g_zzzz_xyyyy, g_zzzz_xyyyz, g_zzzz_xyyzz, g_zzzz_xyzzz, g_zzzz_xzzzz, g_zzzz_yyyyy, g_zzzz_yyyyz, g_zzzz_yyyzz, g_zzzz_yyzzz, g_zzzz_yzzzz, g_zzzz_zzzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzz_xxxxx[k] = g_zzzz_xxxxx[k] - g_0_z_zzzz_xxxxx[k] * ab_z + g_0_z_zzzz_xxxxxz[k];

                g_0_z_zzzzz_xxxxy[k] = g_zzzz_xxxxy[k] - g_0_z_zzzz_xxxxy[k] * ab_z + g_0_z_zzzz_xxxxyz[k];

                g_0_z_zzzzz_xxxxz[k] = g_zzzz_xxxxz[k] - g_0_z_zzzz_xxxxz[k] * ab_z + g_0_z_zzzz_xxxxzz[k];

                g_0_z_zzzzz_xxxyy[k] = g_zzzz_xxxyy[k] - g_0_z_zzzz_xxxyy[k] * ab_z + g_0_z_zzzz_xxxyyz[k];

                g_0_z_zzzzz_xxxyz[k] = g_zzzz_xxxyz[k] - g_0_z_zzzz_xxxyz[k] * ab_z + g_0_z_zzzz_xxxyzz[k];

                g_0_z_zzzzz_xxxzz[k] = g_zzzz_xxxzz[k] - g_0_z_zzzz_xxxzz[k] * ab_z + g_0_z_zzzz_xxxzzz[k];

                g_0_z_zzzzz_xxyyy[k] = g_zzzz_xxyyy[k] - g_0_z_zzzz_xxyyy[k] * ab_z + g_0_z_zzzz_xxyyyz[k];

                g_0_z_zzzzz_xxyyz[k] = g_zzzz_xxyyz[k] - g_0_z_zzzz_xxyyz[k] * ab_z + g_0_z_zzzz_xxyyzz[k];

                g_0_z_zzzzz_xxyzz[k] = g_zzzz_xxyzz[k] - g_0_z_zzzz_xxyzz[k] * ab_z + g_0_z_zzzz_xxyzzz[k];

                g_0_z_zzzzz_xxzzz[k] = g_zzzz_xxzzz[k] - g_0_z_zzzz_xxzzz[k] * ab_z + g_0_z_zzzz_xxzzzz[k];

                g_0_z_zzzzz_xyyyy[k] = g_zzzz_xyyyy[k] - g_0_z_zzzz_xyyyy[k] * ab_z + g_0_z_zzzz_xyyyyz[k];

                g_0_z_zzzzz_xyyyz[k] = g_zzzz_xyyyz[k] - g_0_z_zzzz_xyyyz[k] * ab_z + g_0_z_zzzz_xyyyzz[k];

                g_0_z_zzzzz_xyyzz[k] = g_zzzz_xyyzz[k] - g_0_z_zzzz_xyyzz[k] * ab_z + g_0_z_zzzz_xyyzzz[k];

                g_0_z_zzzzz_xyzzz[k] = g_zzzz_xyzzz[k] - g_0_z_zzzz_xyzzz[k] * ab_z + g_0_z_zzzz_xyzzzz[k];

                g_0_z_zzzzz_xzzzz[k] = g_zzzz_xzzzz[k] - g_0_z_zzzz_xzzzz[k] * ab_z + g_0_z_zzzz_xzzzzz[k];

                g_0_z_zzzzz_yyyyy[k] = g_zzzz_yyyyy[k] - g_0_z_zzzz_yyyyy[k] * ab_z + g_0_z_zzzz_yyyyyz[k];

                g_0_z_zzzzz_yyyyz[k] = g_zzzz_yyyyz[k] - g_0_z_zzzz_yyyyz[k] * ab_z + g_0_z_zzzz_yyyyzz[k];

                g_0_z_zzzzz_yyyzz[k] = g_zzzz_yyyzz[k] - g_0_z_zzzz_yyyzz[k] * ab_z + g_0_z_zzzz_yyyzzz[k];

                g_0_z_zzzzz_yyzzz[k] = g_zzzz_yyzzz[k] - g_0_z_zzzz_yyzzz[k] * ab_z + g_0_z_zzzz_yyzzzz[k];

                g_0_z_zzzzz_yzzzz[k] = g_zzzz_yzzzz[k] - g_0_z_zzzz_yzzzz[k] * ab_z + g_0_z_zzzz_yzzzzz[k];

                g_0_z_zzzzz_zzzzz[k] = g_zzzz_zzzzz[k] - g_0_z_zzzz_zzzzz[k] * ab_z + g_0_z_zzzz_zzzzzz[k];
            }
        }
    }
}

} // erirec namespace

