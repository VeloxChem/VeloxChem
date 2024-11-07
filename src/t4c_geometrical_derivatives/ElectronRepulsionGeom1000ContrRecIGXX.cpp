#include "ElectronRepulsionGeom1000ContrRecIGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom10_hrr_electron_repulsion_igxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_10_igxx,
                                            const size_t idx_hgxx,
                                            const size_t idx_geom_10_hgxx,
                                            const size_t idx_geom_10_hhxx,
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
            /// Set up components of auxilary buffer : HGSS

            const auto hg_off = idx_hgxx + i * dcomps + j;

            auto g_xxxxx_xxxx = cbuffer.data(hg_off + 0 * ccomps * dcomps);

            auto g_xxxxx_xxxy = cbuffer.data(hg_off + 1 * ccomps * dcomps);

            auto g_xxxxx_xxxz = cbuffer.data(hg_off + 2 * ccomps * dcomps);

            auto g_xxxxx_xxyy = cbuffer.data(hg_off + 3 * ccomps * dcomps);

            auto g_xxxxx_xxyz = cbuffer.data(hg_off + 4 * ccomps * dcomps);

            auto g_xxxxx_xxzz = cbuffer.data(hg_off + 5 * ccomps * dcomps);

            auto g_xxxxx_xyyy = cbuffer.data(hg_off + 6 * ccomps * dcomps);

            auto g_xxxxx_xyyz = cbuffer.data(hg_off + 7 * ccomps * dcomps);

            auto g_xxxxx_xyzz = cbuffer.data(hg_off + 8 * ccomps * dcomps);

            auto g_xxxxx_xzzz = cbuffer.data(hg_off + 9 * ccomps * dcomps);

            auto g_xxxxx_yyyy = cbuffer.data(hg_off + 10 * ccomps * dcomps);

            auto g_xxxxx_yyyz = cbuffer.data(hg_off + 11 * ccomps * dcomps);

            auto g_xxxxx_yyzz = cbuffer.data(hg_off + 12 * ccomps * dcomps);

            auto g_xxxxx_yzzz = cbuffer.data(hg_off + 13 * ccomps * dcomps);

            auto g_xxxxx_zzzz = cbuffer.data(hg_off + 14 * ccomps * dcomps);

            auto g_xxxxy_xxxx = cbuffer.data(hg_off + 15 * ccomps * dcomps);

            auto g_xxxxy_xxxy = cbuffer.data(hg_off + 16 * ccomps * dcomps);

            auto g_xxxxy_xxxz = cbuffer.data(hg_off + 17 * ccomps * dcomps);

            auto g_xxxxy_xxyy = cbuffer.data(hg_off + 18 * ccomps * dcomps);

            auto g_xxxxy_xxyz = cbuffer.data(hg_off + 19 * ccomps * dcomps);

            auto g_xxxxy_xxzz = cbuffer.data(hg_off + 20 * ccomps * dcomps);

            auto g_xxxxy_xyyy = cbuffer.data(hg_off + 21 * ccomps * dcomps);

            auto g_xxxxy_xyyz = cbuffer.data(hg_off + 22 * ccomps * dcomps);

            auto g_xxxxy_xyzz = cbuffer.data(hg_off + 23 * ccomps * dcomps);

            auto g_xxxxy_xzzz = cbuffer.data(hg_off + 24 * ccomps * dcomps);

            auto g_xxxxy_yyyy = cbuffer.data(hg_off + 25 * ccomps * dcomps);

            auto g_xxxxy_yyyz = cbuffer.data(hg_off + 26 * ccomps * dcomps);

            auto g_xxxxy_yyzz = cbuffer.data(hg_off + 27 * ccomps * dcomps);

            auto g_xxxxy_yzzz = cbuffer.data(hg_off + 28 * ccomps * dcomps);

            auto g_xxxxy_zzzz = cbuffer.data(hg_off + 29 * ccomps * dcomps);

            auto g_xxxxz_xxxx = cbuffer.data(hg_off + 30 * ccomps * dcomps);

            auto g_xxxxz_xxxy = cbuffer.data(hg_off + 31 * ccomps * dcomps);

            auto g_xxxxz_xxxz = cbuffer.data(hg_off + 32 * ccomps * dcomps);

            auto g_xxxxz_xxyy = cbuffer.data(hg_off + 33 * ccomps * dcomps);

            auto g_xxxxz_xxyz = cbuffer.data(hg_off + 34 * ccomps * dcomps);

            auto g_xxxxz_xxzz = cbuffer.data(hg_off + 35 * ccomps * dcomps);

            auto g_xxxxz_xyyy = cbuffer.data(hg_off + 36 * ccomps * dcomps);

            auto g_xxxxz_xyyz = cbuffer.data(hg_off + 37 * ccomps * dcomps);

            auto g_xxxxz_xyzz = cbuffer.data(hg_off + 38 * ccomps * dcomps);

            auto g_xxxxz_xzzz = cbuffer.data(hg_off + 39 * ccomps * dcomps);

            auto g_xxxxz_yyyy = cbuffer.data(hg_off + 40 * ccomps * dcomps);

            auto g_xxxxz_yyyz = cbuffer.data(hg_off + 41 * ccomps * dcomps);

            auto g_xxxxz_yyzz = cbuffer.data(hg_off + 42 * ccomps * dcomps);

            auto g_xxxxz_yzzz = cbuffer.data(hg_off + 43 * ccomps * dcomps);

            auto g_xxxxz_zzzz = cbuffer.data(hg_off + 44 * ccomps * dcomps);

            auto g_xxxyy_xxxx = cbuffer.data(hg_off + 45 * ccomps * dcomps);

            auto g_xxxyy_xxxy = cbuffer.data(hg_off + 46 * ccomps * dcomps);

            auto g_xxxyy_xxxz = cbuffer.data(hg_off + 47 * ccomps * dcomps);

            auto g_xxxyy_xxyy = cbuffer.data(hg_off + 48 * ccomps * dcomps);

            auto g_xxxyy_xxyz = cbuffer.data(hg_off + 49 * ccomps * dcomps);

            auto g_xxxyy_xxzz = cbuffer.data(hg_off + 50 * ccomps * dcomps);

            auto g_xxxyy_xyyy = cbuffer.data(hg_off + 51 * ccomps * dcomps);

            auto g_xxxyy_xyyz = cbuffer.data(hg_off + 52 * ccomps * dcomps);

            auto g_xxxyy_xyzz = cbuffer.data(hg_off + 53 * ccomps * dcomps);

            auto g_xxxyy_xzzz = cbuffer.data(hg_off + 54 * ccomps * dcomps);

            auto g_xxxyy_yyyy = cbuffer.data(hg_off + 55 * ccomps * dcomps);

            auto g_xxxyy_yyyz = cbuffer.data(hg_off + 56 * ccomps * dcomps);

            auto g_xxxyy_yyzz = cbuffer.data(hg_off + 57 * ccomps * dcomps);

            auto g_xxxyy_yzzz = cbuffer.data(hg_off + 58 * ccomps * dcomps);

            auto g_xxxyy_zzzz = cbuffer.data(hg_off + 59 * ccomps * dcomps);

            auto g_xxxyz_xxxx = cbuffer.data(hg_off + 60 * ccomps * dcomps);

            auto g_xxxyz_xxxy = cbuffer.data(hg_off + 61 * ccomps * dcomps);

            auto g_xxxyz_xxxz = cbuffer.data(hg_off + 62 * ccomps * dcomps);

            auto g_xxxyz_xxyy = cbuffer.data(hg_off + 63 * ccomps * dcomps);

            auto g_xxxyz_xxyz = cbuffer.data(hg_off + 64 * ccomps * dcomps);

            auto g_xxxyz_xxzz = cbuffer.data(hg_off + 65 * ccomps * dcomps);

            auto g_xxxyz_xyyy = cbuffer.data(hg_off + 66 * ccomps * dcomps);

            auto g_xxxyz_xyyz = cbuffer.data(hg_off + 67 * ccomps * dcomps);

            auto g_xxxyz_xyzz = cbuffer.data(hg_off + 68 * ccomps * dcomps);

            auto g_xxxyz_xzzz = cbuffer.data(hg_off + 69 * ccomps * dcomps);

            auto g_xxxyz_yyyy = cbuffer.data(hg_off + 70 * ccomps * dcomps);

            auto g_xxxyz_yyyz = cbuffer.data(hg_off + 71 * ccomps * dcomps);

            auto g_xxxyz_yyzz = cbuffer.data(hg_off + 72 * ccomps * dcomps);

            auto g_xxxyz_yzzz = cbuffer.data(hg_off + 73 * ccomps * dcomps);

            auto g_xxxyz_zzzz = cbuffer.data(hg_off + 74 * ccomps * dcomps);

            auto g_xxxzz_xxxx = cbuffer.data(hg_off + 75 * ccomps * dcomps);

            auto g_xxxzz_xxxy = cbuffer.data(hg_off + 76 * ccomps * dcomps);

            auto g_xxxzz_xxxz = cbuffer.data(hg_off + 77 * ccomps * dcomps);

            auto g_xxxzz_xxyy = cbuffer.data(hg_off + 78 * ccomps * dcomps);

            auto g_xxxzz_xxyz = cbuffer.data(hg_off + 79 * ccomps * dcomps);

            auto g_xxxzz_xxzz = cbuffer.data(hg_off + 80 * ccomps * dcomps);

            auto g_xxxzz_xyyy = cbuffer.data(hg_off + 81 * ccomps * dcomps);

            auto g_xxxzz_xyyz = cbuffer.data(hg_off + 82 * ccomps * dcomps);

            auto g_xxxzz_xyzz = cbuffer.data(hg_off + 83 * ccomps * dcomps);

            auto g_xxxzz_xzzz = cbuffer.data(hg_off + 84 * ccomps * dcomps);

            auto g_xxxzz_yyyy = cbuffer.data(hg_off + 85 * ccomps * dcomps);

            auto g_xxxzz_yyyz = cbuffer.data(hg_off + 86 * ccomps * dcomps);

            auto g_xxxzz_yyzz = cbuffer.data(hg_off + 87 * ccomps * dcomps);

            auto g_xxxzz_yzzz = cbuffer.data(hg_off + 88 * ccomps * dcomps);

            auto g_xxxzz_zzzz = cbuffer.data(hg_off + 89 * ccomps * dcomps);

            auto g_xxyyy_xxxx = cbuffer.data(hg_off + 90 * ccomps * dcomps);

            auto g_xxyyy_xxxy = cbuffer.data(hg_off + 91 * ccomps * dcomps);

            auto g_xxyyy_xxxz = cbuffer.data(hg_off + 92 * ccomps * dcomps);

            auto g_xxyyy_xxyy = cbuffer.data(hg_off + 93 * ccomps * dcomps);

            auto g_xxyyy_xxyz = cbuffer.data(hg_off + 94 * ccomps * dcomps);

            auto g_xxyyy_xxzz = cbuffer.data(hg_off + 95 * ccomps * dcomps);

            auto g_xxyyy_xyyy = cbuffer.data(hg_off + 96 * ccomps * dcomps);

            auto g_xxyyy_xyyz = cbuffer.data(hg_off + 97 * ccomps * dcomps);

            auto g_xxyyy_xyzz = cbuffer.data(hg_off + 98 * ccomps * dcomps);

            auto g_xxyyy_xzzz = cbuffer.data(hg_off + 99 * ccomps * dcomps);

            auto g_xxyyy_yyyy = cbuffer.data(hg_off + 100 * ccomps * dcomps);

            auto g_xxyyy_yyyz = cbuffer.data(hg_off + 101 * ccomps * dcomps);

            auto g_xxyyy_yyzz = cbuffer.data(hg_off + 102 * ccomps * dcomps);

            auto g_xxyyy_yzzz = cbuffer.data(hg_off + 103 * ccomps * dcomps);

            auto g_xxyyy_zzzz = cbuffer.data(hg_off + 104 * ccomps * dcomps);

            auto g_xxyyz_xxxx = cbuffer.data(hg_off + 105 * ccomps * dcomps);

            auto g_xxyyz_xxxy = cbuffer.data(hg_off + 106 * ccomps * dcomps);

            auto g_xxyyz_xxxz = cbuffer.data(hg_off + 107 * ccomps * dcomps);

            auto g_xxyyz_xxyy = cbuffer.data(hg_off + 108 * ccomps * dcomps);

            auto g_xxyyz_xxyz = cbuffer.data(hg_off + 109 * ccomps * dcomps);

            auto g_xxyyz_xxzz = cbuffer.data(hg_off + 110 * ccomps * dcomps);

            auto g_xxyyz_xyyy = cbuffer.data(hg_off + 111 * ccomps * dcomps);

            auto g_xxyyz_xyyz = cbuffer.data(hg_off + 112 * ccomps * dcomps);

            auto g_xxyyz_xyzz = cbuffer.data(hg_off + 113 * ccomps * dcomps);

            auto g_xxyyz_xzzz = cbuffer.data(hg_off + 114 * ccomps * dcomps);

            auto g_xxyyz_yyyy = cbuffer.data(hg_off + 115 * ccomps * dcomps);

            auto g_xxyyz_yyyz = cbuffer.data(hg_off + 116 * ccomps * dcomps);

            auto g_xxyyz_yyzz = cbuffer.data(hg_off + 117 * ccomps * dcomps);

            auto g_xxyyz_yzzz = cbuffer.data(hg_off + 118 * ccomps * dcomps);

            auto g_xxyyz_zzzz = cbuffer.data(hg_off + 119 * ccomps * dcomps);

            auto g_xxyzz_xxxx = cbuffer.data(hg_off + 120 * ccomps * dcomps);

            auto g_xxyzz_xxxy = cbuffer.data(hg_off + 121 * ccomps * dcomps);

            auto g_xxyzz_xxxz = cbuffer.data(hg_off + 122 * ccomps * dcomps);

            auto g_xxyzz_xxyy = cbuffer.data(hg_off + 123 * ccomps * dcomps);

            auto g_xxyzz_xxyz = cbuffer.data(hg_off + 124 * ccomps * dcomps);

            auto g_xxyzz_xxzz = cbuffer.data(hg_off + 125 * ccomps * dcomps);

            auto g_xxyzz_xyyy = cbuffer.data(hg_off + 126 * ccomps * dcomps);

            auto g_xxyzz_xyyz = cbuffer.data(hg_off + 127 * ccomps * dcomps);

            auto g_xxyzz_xyzz = cbuffer.data(hg_off + 128 * ccomps * dcomps);

            auto g_xxyzz_xzzz = cbuffer.data(hg_off + 129 * ccomps * dcomps);

            auto g_xxyzz_yyyy = cbuffer.data(hg_off + 130 * ccomps * dcomps);

            auto g_xxyzz_yyyz = cbuffer.data(hg_off + 131 * ccomps * dcomps);

            auto g_xxyzz_yyzz = cbuffer.data(hg_off + 132 * ccomps * dcomps);

            auto g_xxyzz_yzzz = cbuffer.data(hg_off + 133 * ccomps * dcomps);

            auto g_xxyzz_zzzz = cbuffer.data(hg_off + 134 * ccomps * dcomps);

            auto g_xxzzz_xxxx = cbuffer.data(hg_off + 135 * ccomps * dcomps);

            auto g_xxzzz_xxxy = cbuffer.data(hg_off + 136 * ccomps * dcomps);

            auto g_xxzzz_xxxz = cbuffer.data(hg_off + 137 * ccomps * dcomps);

            auto g_xxzzz_xxyy = cbuffer.data(hg_off + 138 * ccomps * dcomps);

            auto g_xxzzz_xxyz = cbuffer.data(hg_off + 139 * ccomps * dcomps);

            auto g_xxzzz_xxzz = cbuffer.data(hg_off + 140 * ccomps * dcomps);

            auto g_xxzzz_xyyy = cbuffer.data(hg_off + 141 * ccomps * dcomps);

            auto g_xxzzz_xyyz = cbuffer.data(hg_off + 142 * ccomps * dcomps);

            auto g_xxzzz_xyzz = cbuffer.data(hg_off + 143 * ccomps * dcomps);

            auto g_xxzzz_xzzz = cbuffer.data(hg_off + 144 * ccomps * dcomps);

            auto g_xxzzz_yyyy = cbuffer.data(hg_off + 145 * ccomps * dcomps);

            auto g_xxzzz_yyyz = cbuffer.data(hg_off + 146 * ccomps * dcomps);

            auto g_xxzzz_yyzz = cbuffer.data(hg_off + 147 * ccomps * dcomps);

            auto g_xxzzz_yzzz = cbuffer.data(hg_off + 148 * ccomps * dcomps);

            auto g_xxzzz_zzzz = cbuffer.data(hg_off + 149 * ccomps * dcomps);

            auto g_xyyyy_xxxx = cbuffer.data(hg_off + 150 * ccomps * dcomps);

            auto g_xyyyy_xxxy = cbuffer.data(hg_off + 151 * ccomps * dcomps);

            auto g_xyyyy_xxxz = cbuffer.data(hg_off + 152 * ccomps * dcomps);

            auto g_xyyyy_xxyy = cbuffer.data(hg_off + 153 * ccomps * dcomps);

            auto g_xyyyy_xxyz = cbuffer.data(hg_off + 154 * ccomps * dcomps);

            auto g_xyyyy_xxzz = cbuffer.data(hg_off + 155 * ccomps * dcomps);

            auto g_xyyyy_xyyy = cbuffer.data(hg_off + 156 * ccomps * dcomps);

            auto g_xyyyy_xyyz = cbuffer.data(hg_off + 157 * ccomps * dcomps);

            auto g_xyyyy_xyzz = cbuffer.data(hg_off + 158 * ccomps * dcomps);

            auto g_xyyyy_xzzz = cbuffer.data(hg_off + 159 * ccomps * dcomps);

            auto g_xyyyy_yyyy = cbuffer.data(hg_off + 160 * ccomps * dcomps);

            auto g_xyyyy_yyyz = cbuffer.data(hg_off + 161 * ccomps * dcomps);

            auto g_xyyyy_yyzz = cbuffer.data(hg_off + 162 * ccomps * dcomps);

            auto g_xyyyy_yzzz = cbuffer.data(hg_off + 163 * ccomps * dcomps);

            auto g_xyyyy_zzzz = cbuffer.data(hg_off + 164 * ccomps * dcomps);

            auto g_xyyyz_xxxx = cbuffer.data(hg_off + 165 * ccomps * dcomps);

            auto g_xyyyz_xxxy = cbuffer.data(hg_off + 166 * ccomps * dcomps);

            auto g_xyyyz_xxxz = cbuffer.data(hg_off + 167 * ccomps * dcomps);

            auto g_xyyyz_xxyy = cbuffer.data(hg_off + 168 * ccomps * dcomps);

            auto g_xyyyz_xxyz = cbuffer.data(hg_off + 169 * ccomps * dcomps);

            auto g_xyyyz_xxzz = cbuffer.data(hg_off + 170 * ccomps * dcomps);

            auto g_xyyyz_xyyy = cbuffer.data(hg_off + 171 * ccomps * dcomps);

            auto g_xyyyz_xyyz = cbuffer.data(hg_off + 172 * ccomps * dcomps);

            auto g_xyyyz_xyzz = cbuffer.data(hg_off + 173 * ccomps * dcomps);

            auto g_xyyyz_xzzz = cbuffer.data(hg_off + 174 * ccomps * dcomps);

            auto g_xyyyz_yyyy = cbuffer.data(hg_off + 175 * ccomps * dcomps);

            auto g_xyyyz_yyyz = cbuffer.data(hg_off + 176 * ccomps * dcomps);

            auto g_xyyyz_yyzz = cbuffer.data(hg_off + 177 * ccomps * dcomps);

            auto g_xyyyz_yzzz = cbuffer.data(hg_off + 178 * ccomps * dcomps);

            auto g_xyyyz_zzzz = cbuffer.data(hg_off + 179 * ccomps * dcomps);

            auto g_xyyzz_xxxx = cbuffer.data(hg_off + 180 * ccomps * dcomps);

            auto g_xyyzz_xxxy = cbuffer.data(hg_off + 181 * ccomps * dcomps);

            auto g_xyyzz_xxxz = cbuffer.data(hg_off + 182 * ccomps * dcomps);

            auto g_xyyzz_xxyy = cbuffer.data(hg_off + 183 * ccomps * dcomps);

            auto g_xyyzz_xxyz = cbuffer.data(hg_off + 184 * ccomps * dcomps);

            auto g_xyyzz_xxzz = cbuffer.data(hg_off + 185 * ccomps * dcomps);

            auto g_xyyzz_xyyy = cbuffer.data(hg_off + 186 * ccomps * dcomps);

            auto g_xyyzz_xyyz = cbuffer.data(hg_off + 187 * ccomps * dcomps);

            auto g_xyyzz_xyzz = cbuffer.data(hg_off + 188 * ccomps * dcomps);

            auto g_xyyzz_xzzz = cbuffer.data(hg_off + 189 * ccomps * dcomps);

            auto g_xyyzz_yyyy = cbuffer.data(hg_off + 190 * ccomps * dcomps);

            auto g_xyyzz_yyyz = cbuffer.data(hg_off + 191 * ccomps * dcomps);

            auto g_xyyzz_yyzz = cbuffer.data(hg_off + 192 * ccomps * dcomps);

            auto g_xyyzz_yzzz = cbuffer.data(hg_off + 193 * ccomps * dcomps);

            auto g_xyyzz_zzzz = cbuffer.data(hg_off + 194 * ccomps * dcomps);

            auto g_xyzzz_xxxx = cbuffer.data(hg_off + 195 * ccomps * dcomps);

            auto g_xyzzz_xxxy = cbuffer.data(hg_off + 196 * ccomps * dcomps);

            auto g_xyzzz_xxxz = cbuffer.data(hg_off + 197 * ccomps * dcomps);

            auto g_xyzzz_xxyy = cbuffer.data(hg_off + 198 * ccomps * dcomps);

            auto g_xyzzz_xxyz = cbuffer.data(hg_off + 199 * ccomps * dcomps);

            auto g_xyzzz_xxzz = cbuffer.data(hg_off + 200 * ccomps * dcomps);

            auto g_xyzzz_xyyy = cbuffer.data(hg_off + 201 * ccomps * dcomps);

            auto g_xyzzz_xyyz = cbuffer.data(hg_off + 202 * ccomps * dcomps);

            auto g_xyzzz_xyzz = cbuffer.data(hg_off + 203 * ccomps * dcomps);

            auto g_xyzzz_xzzz = cbuffer.data(hg_off + 204 * ccomps * dcomps);

            auto g_xyzzz_yyyy = cbuffer.data(hg_off + 205 * ccomps * dcomps);

            auto g_xyzzz_yyyz = cbuffer.data(hg_off + 206 * ccomps * dcomps);

            auto g_xyzzz_yyzz = cbuffer.data(hg_off + 207 * ccomps * dcomps);

            auto g_xyzzz_yzzz = cbuffer.data(hg_off + 208 * ccomps * dcomps);

            auto g_xyzzz_zzzz = cbuffer.data(hg_off + 209 * ccomps * dcomps);

            auto g_xzzzz_xxxx = cbuffer.data(hg_off + 210 * ccomps * dcomps);

            auto g_xzzzz_xxxy = cbuffer.data(hg_off + 211 * ccomps * dcomps);

            auto g_xzzzz_xxxz = cbuffer.data(hg_off + 212 * ccomps * dcomps);

            auto g_xzzzz_xxyy = cbuffer.data(hg_off + 213 * ccomps * dcomps);

            auto g_xzzzz_xxyz = cbuffer.data(hg_off + 214 * ccomps * dcomps);

            auto g_xzzzz_xxzz = cbuffer.data(hg_off + 215 * ccomps * dcomps);

            auto g_xzzzz_xyyy = cbuffer.data(hg_off + 216 * ccomps * dcomps);

            auto g_xzzzz_xyyz = cbuffer.data(hg_off + 217 * ccomps * dcomps);

            auto g_xzzzz_xyzz = cbuffer.data(hg_off + 218 * ccomps * dcomps);

            auto g_xzzzz_xzzz = cbuffer.data(hg_off + 219 * ccomps * dcomps);

            auto g_xzzzz_yyyy = cbuffer.data(hg_off + 220 * ccomps * dcomps);

            auto g_xzzzz_yyyz = cbuffer.data(hg_off + 221 * ccomps * dcomps);

            auto g_xzzzz_yyzz = cbuffer.data(hg_off + 222 * ccomps * dcomps);

            auto g_xzzzz_yzzz = cbuffer.data(hg_off + 223 * ccomps * dcomps);

            auto g_xzzzz_zzzz = cbuffer.data(hg_off + 224 * ccomps * dcomps);

            auto g_yyyyy_xxxx = cbuffer.data(hg_off + 225 * ccomps * dcomps);

            auto g_yyyyy_xxxy = cbuffer.data(hg_off + 226 * ccomps * dcomps);

            auto g_yyyyy_xxxz = cbuffer.data(hg_off + 227 * ccomps * dcomps);

            auto g_yyyyy_xxyy = cbuffer.data(hg_off + 228 * ccomps * dcomps);

            auto g_yyyyy_xxyz = cbuffer.data(hg_off + 229 * ccomps * dcomps);

            auto g_yyyyy_xxzz = cbuffer.data(hg_off + 230 * ccomps * dcomps);

            auto g_yyyyy_xyyy = cbuffer.data(hg_off + 231 * ccomps * dcomps);

            auto g_yyyyy_xyyz = cbuffer.data(hg_off + 232 * ccomps * dcomps);

            auto g_yyyyy_xyzz = cbuffer.data(hg_off + 233 * ccomps * dcomps);

            auto g_yyyyy_xzzz = cbuffer.data(hg_off + 234 * ccomps * dcomps);

            auto g_yyyyy_yyyy = cbuffer.data(hg_off + 235 * ccomps * dcomps);

            auto g_yyyyy_yyyz = cbuffer.data(hg_off + 236 * ccomps * dcomps);

            auto g_yyyyy_yyzz = cbuffer.data(hg_off + 237 * ccomps * dcomps);

            auto g_yyyyy_yzzz = cbuffer.data(hg_off + 238 * ccomps * dcomps);

            auto g_yyyyy_zzzz = cbuffer.data(hg_off + 239 * ccomps * dcomps);

            auto g_yyyyz_xxxx = cbuffer.data(hg_off + 240 * ccomps * dcomps);

            auto g_yyyyz_xxxy = cbuffer.data(hg_off + 241 * ccomps * dcomps);

            auto g_yyyyz_xxxz = cbuffer.data(hg_off + 242 * ccomps * dcomps);

            auto g_yyyyz_xxyy = cbuffer.data(hg_off + 243 * ccomps * dcomps);

            auto g_yyyyz_xxyz = cbuffer.data(hg_off + 244 * ccomps * dcomps);

            auto g_yyyyz_xxzz = cbuffer.data(hg_off + 245 * ccomps * dcomps);

            auto g_yyyyz_xyyy = cbuffer.data(hg_off + 246 * ccomps * dcomps);

            auto g_yyyyz_xyyz = cbuffer.data(hg_off + 247 * ccomps * dcomps);

            auto g_yyyyz_xyzz = cbuffer.data(hg_off + 248 * ccomps * dcomps);

            auto g_yyyyz_xzzz = cbuffer.data(hg_off + 249 * ccomps * dcomps);

            auto g_yyyyz_yyyy = cbuffer.data(hg_off + 250 * ccomps * dcomps);

            auto g_yyyyz_yyyz = cbuffer.data(hg_off + 251 * ccomps * dcomps);

            auto g_yyyyz_yyzz = cbuffer.data(hg_off + 252 * ccomps * dcomps);

            auto g_yyyyz_yzzz = cbuffer.data(hg_off + 253 * ccomps * dcomps);

            auto g_yyyyz_zzzz = cbuffer.data(hg_off + 254 * ccomps * dcomps);

            auto g_yyyzz_xxxx = cbuffer.data(hg_off + 255 * ccomps * dcomps);

            auto g_yyyzz_xxxy = cbuffer.data(hg_off + 256 * ccomps * dcomps);

            auto g_yyyzz_xxxz = cbuffer.data(hg_off + 257 * ccomps * dcomps);

            auto g_yyyzz_xxyy = cbuffer.data(hg_off + 258 * ccomps * dcomps);

            auto g_yyyzz_xxyz = cbuffer.data(hg_off + 259 * ccomps * dcomps);

            auto g_yyyzz_xxzz = cbuffer.data(hg_off + 260 * ccomps * dcomps);

            auto g_yyyzz_xyyy = cbuffer.data(hg_off + 261 * ccomps * dcomps);

            auto g_yyyzz_xyyz = cbuffer.data(hg_off + 262 * ccomps * dcomps);

            auto g_yyyzz_xyzz = cbuffer.data(hg_off + 263 * ccomps * dcomps);

            auto g_yyyzz_xzzz = cbuffer.data(hg_off + 264 * ccomps * dcomps);

            auto g_yyyzz_yyyy = cbuffer.data(hg_off + 265 * ccomps * dcomps);

            auto g_yyyzz_yyyz = cbuffer.data(hg_off + 266 * ccomps * dcomps);

            auto g_yyyzz_yyzz = cbuffer.data(hg_off + 267 * ccomps * dcomps);

            auto g_yyyzz_yzzz = cbuffer.data(hg_off + 268 * ccomps * dcomps);

            auto g_yyyzz_zzzz = cbuffer.data(hg_off + 269 * ccomps * dcomps);

            auto g_yyzzz_xxxx = cbuffer.data(hg_off + 270 * ccomps * dcomps);

            auto g_yyzzz_xxxy = cbuffer.data(hg_off + 271 * ccomps * dcomps);

            auto g_yyzzz_xxxz = cbuffer.data(hg_off + 272 * ccomps * dcomps);

            auto g_yyzzz_xxyy = cbuffer.data(hg_off + 273 * ccomps * dcomps);

            auto g_yyzzz_xxyz = cbuffer.data(hg_off + 274 * ccomps * dcomps);

            auto g_yyzzz_xxzz = cbuffer.data(hg_off + 275 * ccomps * dcomps);

            auto g_yyzzz_xyyy = cbuffer.data(hg_off + 276 * ccomps * dcomps);

            auto g_yyzzz_xyyz = cbuffer.data(hg_off + 277 * ccomps * dcomps);

            auto g_yyzzz_xyzz = cbuffer.data(hg_off + 278 * ccomps * dcomps);

            auto g_yyzzz_xzzz = cbuffer.data(hg_off + 279 * ccomps * dcomps);

            auto g_yyzzz_yyyy = cbuffer.data(hg_off + 280 * ccomps * dcomps);

            auto g_yyzzz_yyyz = cbuffer.data(hg_off + 281 * ccomps * dcomps);

            auto g_yyzzz_yyzz = cbuffer.data(hg_off + 282 * ccomps * dcomps);

            auto g_yyzzz_yzzz = cbuffer.data(hg_off + 283 * ccomps * dcomps);

            auto g_yyzzz_zzzz = cbuffer.data(hg_off + 284 * ccomps * dcomps);

            auto g_yzzzz_xxxx = cbuffer.data(hg_off + 285 * ccomps * dcomps);

            auto g_yzzzz_xxxy = cbuffer.data(hg_off + 286 * ccomps * dcomps);

            auto g_yzzzz_xxxz = cbuffer.data(hg_off + 287 * ccomps * dcomps);

            auto g_yzzzz_xxyy = cbuffer.data(hg_off + 288 * ccomps * dcomps);

            auto g_yzzzz_xxyz = cbuffer.data(hg_off + 289 * ccomps * dcomps);

            auto g_yzzzz_xxzz = cbuffer.data(hg_off + 290 * ccomps * dcomps);

            auto g_yzzzz_xyyy = cbuffer.data(hg_off + 291 * ccomps * dcomps);

            auto g_yzzzz_xyyz = cbuffer.data(hg_off + 292 * ccomps * dcomps);

            auto g_yzzzz_xyzz = cbuffer.data(hg_off + 293 * ccomps * dcomps);

            auto g_yzzzz_xzzz = cbuffer.data(hg_off + 294 * ccomps * dcomps);

            auto g_yzzzz_yyyy = cbuffer.data(hg_off + 295 * ccomps * dcomps);

            auto g_yzzzz_yyyz = cbuffer.data(hg_off + 296 * ccomps * dcomps);

            auto g_yzzzz_yyzz = cbuffer.data(hg_off + 297 * ccomps * dcomps);

            auto g_yzzzz_yzzz = cbuffer.data(hg_off + 298 * ccomps * dcomps);

            auto g_yzzzz_zzzz = cbuffer.data(hg_off + 299 * ccomps * dcomps);

            auto g_zzzzz_xxxx = cbuffer.data(hg_off + 300 * ccomps * dcomps);

            auto g_zzzzz_xxxy = cbuffer.data(hg_off + 301 * ccomps * dcomps);

            auto g_zzzzz_xxxz = cbuffer.data(hg_off + 302 * ccomps * dcomps);

            auto g_zzzzz_xxyy = cbuffer.data(hg_off + 303 * ccomps * dcomps);

            auto g_zzzzz_xxyz = cbuffer.data(hg_off + 304 * ccomps * dcomps);

            auto g_zzzzz_xxzz = cbuffer.data(hg_off + 305 * ccomps * dcomps);

            auto g_zzzzz_xyyy = cbuffer.data(hg_off + 306 * ccomps * dcomps);

            auto g_zzzzz_xyyz = cbuffer.data(hg_off + 307 * ccomps * dcomps);

            auto g_zzzzz_xyzz = cbuffer.data(hg_off + 308 * ccomps * dcomps);

            auto g_zzzzz_xzzz = cbuffer.data(hg_off + 309 * ccomps * dcomps);

            auto g_zzzzz_yyyy = cbuffer.data(hg_off + 310 * ccomps * dcomps);

            auto g_zzzzz_yyyz = cbuffer.data(hg_off + 311 * ccomps * dcomps);

            auto g_zzzzz_yyzz = cbuffer.data(hg_off + 312 * ccomps * dcomps);

            auto g_zzzzz_yzzz = cbuffer.data(hg_off + 313 * ccomps * dcomps);

            auto g_zzzzz_zzzz = cbuffer.data(hg_off + 314 * ccomps * dcomps);

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

            /// Set up components of auxilary buffer : HHSS

            const auto hh_geom_10_off = idx_geom_10_hhxx + i * dcomps + j;

            auto g_x_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 14 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 29 * ccomps * dcomps);

            auto g_x_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 44 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 59 * ccomps * dcomps);

            auto g_x_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 74 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 89 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 104 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 119 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 134 * ccomps * dcomps);

            auto g_x_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 149 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 164 * ccomps * dcomps);

            auto g_x_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 179 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 194 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 209 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 224 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 239 * ccomps * dcomps);

            auto g_x_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 251 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 254 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 269 * ccomps * dcomps);

            auto g_x_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 279 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 284 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 299 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 307 * ccomps * dcomps);

            auto g_x_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 314 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 329 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 335 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 344 * ccomps * dcomps);

            auto g_x_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 356 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 359 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 360 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 362 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 363 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 364 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 365 * ccomps * dcomps);

            auto g_x_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 366 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 369 * ccomps * dcomps);

            auto g_x_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 373 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 374 * ccomps * dcomps);

            auto g_x_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 378 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 380 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 384 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 387 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 389 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 391 * ccomps * dcomps);

            auto g_x_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 392 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 395 * ccomps * dcomps);

            auto g_x_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 396 * ccomps * dcomps);

            auto g_x_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 398 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 401 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 404 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 405 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 410 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 416 * ccomps * dcomps);

            auto g_x_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 419 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 420 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 421 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 422 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 423 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 424 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 425 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 426 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 427 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 428 * ccomps * dcomps);

            auto g_x_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 429 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 430 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 431 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 432 * ccomps * dcomps);

            auto g_x_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 433 * ccomps * dcomps);

            auto g_x_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 434 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 435 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 436 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 437 * ccomps * dcomps);

            auto g_x_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 438 * ccomps * dcomps);

            auto g_x_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 439 * ccomps * dcomps);

            auto g_x_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 449 * ccomps * dcomps);

            auto g_y_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 459 * ccomps * dcomps);

            auto g_y_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 464 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 469 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 475 * ccomps * dcomps);

            auto g_y_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 479 * ccomps * dcomps);

            auto g_y_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 489 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 494 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 499 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 503 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 509 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 519 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 524 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 529 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 531 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 539 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 554 * ccomps * dcomps);

            auto g_y_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 559 * ccomps * dcomps);

            auto g_y_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 560 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 561 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 562 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 563 * ccomps * dcomps);

            auto g_y_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 564 * ccomps * dcomps);

            auto g_y_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 565 * ccomps * dcomps);

            auto g_y_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 566 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 567 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 568 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 569 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 570 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 571 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 572 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 573 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 574 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 575 * ccomps * dcomps);

            auto g_y_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 576 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 577 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 578 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 579 * ccomps * dcomps);

            auto g_y_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 580 * ccomps * dcomps);

            auto g_y_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 581 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 582 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 583 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 584 * ccomps * dcomps);

            auto g_y_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 585 * ccomps * dcomps);

            auto g_y_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 586 * ccomps * dcomps);

            auto g_y_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 587 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 599 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 608 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 609 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 610 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 611 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 614 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 615 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 629 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 630 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 631 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 632 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 633 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 634 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 635 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 636 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 637 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 638 * ccomps * dcomps);

            auto g_y_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 639 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 640 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 641 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 642 * ccomps * dcomps);

            auto g_y_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 643 * ccomps * dcomps);

            auto g_y_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 644 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 645 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 646 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 647 * ccomps * dcomps);

            auto g_y_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 648 * ccomps * dcomps);

            auto g_y_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 649 * ccomps * dcomps);

            auto g_y_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 650 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 651 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 652 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 653 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 654 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 655 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 656 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 657 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 658 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 659 * ccomps * dcomps);

            auto g_y_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 660 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 661 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 662 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 663 * ccomps * dcomps);

            auto g_y_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 664 * ccomps * dcomps);

            auto g_y_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 665 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 666 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 667 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 668 * ccomps * dcomps);

            auto g_y_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 669 * ccomps * dcomps);

            auto g_y_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 670 * ccomps * dcomps);

            auto g_y_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 671 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 672 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 673 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 674 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 675 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 676 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 677 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 678 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 679 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 680 * ccomps * dcomps);

            auto g_y_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 681 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 682 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 683 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 684 * ccomps * dcomps);

            auto g_y_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 685 * ccomps * dcomps);

            auto g_y_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 686 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 687 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 688 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 689 * ccomps * dcomps);

            auto g_y_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 690 * ccomps * dcomps);

            auto g_y_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 691 * ccomps * dcomps);

            auto g_y_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 692 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 693 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 694 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 695 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 696 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 697 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 698 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 699 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 700 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 701 * ccomps * dcomps);

            auto g_y_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 702 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 703 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 704 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 705 * ccomps * dcomps);

            auto g_y_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 706 * ccomps * dcomps);

            auto g_y_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 707 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 708 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 709 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 710 * ccomps * dcomps);

            auto g_y_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 711 * ccomps * dcomps);

            auto g_y_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 712 * ccomps * dcomps);

            auto g_y_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 713 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 714 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 715 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 716 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 717 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 718 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 719 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 720 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 721 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 722 * ccomps * dcomps);

            auto g_y_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 723 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 724 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 725 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 726 * ccomps * dcomps);

            auto g_y_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 727 * ccomps * dcomps);

            auto g_y_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 728 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 729 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 730 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 731 * ccomps * dcomps);

            auto g_y_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 732 * ccomps * dcomps);

            auto g_y_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 733 * ccomps * dcomps);

            auto g_y_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 734 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 735 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 736 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 737 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 738 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 739 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 740 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 741 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 742 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 743 * ccomps * dcomps);

            auto g_y_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 744 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 745 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 746 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 747 * ccomps * dcomps);

            auto g_y_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 748 * ccomps * dcomps);

            auto g_y_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 749 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 750 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 751 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 752 * ccomps * dcomps);

            auto g_y_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 753 * ccomps * dcomps);

            auto g_y_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 754 * ccomps * dcomps);

            auto g_y_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 755 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 756 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 757 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 758 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 759 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 760 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 761 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 762 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 763 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 764 * ccomps * dcomps);

            auto g_y_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 765 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 766 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 767 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 768 * ccomps * dcomps);

            auto g_y_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 769 * ccomps * dcomps);

            auto g_y_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 770 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 771 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 772 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 773 * ccomps * dcomps);

            auto g_y_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 774 * ccomps * dcomps);

            auto g_y_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 775 * ccomps * dcomps);

            auto g_y_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 776 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 777 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 778 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 779 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 780 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 781 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 782 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 783 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 784 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 785 * ccomps * dcomps);

            auto g_y_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 786 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 787 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 788 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 789 * ccomps * dcomps);

            auto g_y_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 790 * ccomps * dcomps);

            auto g_y_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 791 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 794 * ccomps * dcomps);

            auto g_y_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 797 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 805 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 809 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 811 * ccomps * dcomps);

            auto g_y_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 818 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 820 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 821 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 822 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 823 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 824 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 825 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 826 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 827 * ccomps * dcomps);

            auto g_y_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 833 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 839 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 840 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 841 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 842 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 843 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 844 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 845 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 846 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 847 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 848 * ccomps * dcomps);

            auto g_y_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 849 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 850 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 851 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 852 * ccomps * dcomps);

            auto g_y_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 853 * ccomps * dcomps);

            auto g_y_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 854 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 855 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 856 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 857 * ccomps * dcomps);

            auto g_y_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 858 * ccomps * dcomps);

            auto g_y_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 859 * ccomps * dcomps);

            auto g_y_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 860 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 861 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 862 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 863 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 864 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 865 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 866 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 867 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 868 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 869 * ccomps * dcomps);

            auto g_y_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 870 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 871 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 872 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 873 * ccomps * dcomps);

            auto g_y_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 874 * ccomps * dcomps);

            auto g_y_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 875 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 876 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 877 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 878 * ccomps * dcomps);

            auto g_y_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 879 * ccomps * dcomps);

            auto g_y_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 880 * ccomps * dcomps);

            auto g_y_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 881 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxx = cbuffer.data(hh_geom_10_off + 882 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxy = cbuffer.data(hh_geom_10_off + 883 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxxz = cbuffer.data(hh_geom_10_off + 884 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyy = cbuffer.data(hh_geom_10_off + 885 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxyz = cbuffer.data(hh_geom_10_off + 886 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxxzz = cbuffer.data(hh_geom_10_off + 887 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyy = cbuffer.data(hh_geom_10_off + 888 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyyz = cbuffer.data(hh_geom_10_off + 889 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxyzz = cbuffer.data(hh_geom_10_off + 890 * ccomps * dcomps);

            auto g_z_0_xxxxx_xxzzz = cbuffer.data(hh_geom_10_off + 891 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyy = cbuffer.data(hh_geom_10_off + 892 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyyz = cbuffer.data(hh_geom_10_off + 893 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyyzz = cbuffer.data(hh_geom_10_off + 894 * ccomps * dcomps);

            auto g_z_0_xxxxx_xyzzz = cbuffer.data(hh_geom_10_off + 895 * ccomps * dcomps);

            auto g_z_0_xxxxx_xzzzz = cbuffer.data(hh_geom_10_off + 896 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyyy = cbuffer.data(hh_geom_10_off + 897 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyyz = cbuffer.data(hh_geom_10_off + 898 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyyzz = cbuffer.data(hh_geom_10_off + 899 * ccomps * dcomps);

            auto g_z_0_xxxxx_yyzzz = cbuffer.data(hh_geom_10_off + 900 * ccomps * dcomps);

            auto g_z_0_xxxxx_yzzzz = cbuffer.data(hh_geom_10_off + 901 * ccomps * dcomps);

            auto g_z_0_xxxxx_zzzzz = cbuffer.data(hh_geom_10_off + 902 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxx = cbuffer.data(hh_geom_10_off + 903 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxy = cbuffer.data(hh_geom_10_off + 904 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxxz = cbuffer.data(hh_geom_10_off + 905 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyy = cbuffer.data(hh_geom_10_off + 906 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxyz = cbuffer.data(hh_geom_10_off + 907 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxxzz = cbuffer.data(hh_geom_10_off + 908 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyy = cbuffer.data(hh_geom_10_off + 909 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyyz = cbuffer.data(hh_geom_10_off + 910 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxyzz = cbuffer.data(hh_geom_10_off + 911 * ccomps * dcomps);

            auto g_z_0_xxxxy_xxzzz = cbuffer.data(hh_geom_10_off + 912 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyy = cbuffer.data(hh_geom_10_off + 913 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyyz = cbuffer.data(hh_geom_10_off + 914 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyyzz = cbuffer.data(hh_geom_10_off + 915 * ccomps * dcomps);

            auto g_z_0_xxxxy_xyzzz = cbuffer.data(hh_geom_10_off + 916 * ccomps * dcomps);

            auto g_z_0_xxxxy_xzzzz = cbuffer.data(hh_geom_10_off + 917 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyyy = cbuffer.data(hh_geom_10_off + 918 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyyz = cbuffer.data(hh_geom_10_off + 919 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyyzz = cbuffer.data(hh_geom_10_off + 920 * ccomps * dcomps);

            auto g_z_0_xxxxy_yyzzz = cbuffer.data(hh_geom_10_off + 921 * ccomps * dcomps);

            auto g_z_0_xxxxy_yzzzz = cbuffer.data(hh_geom_10_off + 922 * ccomps * dcomps);

            auto g_z_0_xxxxy_zzzzz = cbuffer.data(hh_geom_10_off + 923 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxx = cbuffer.data(hh_geom_10_off + 924 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxy = cbuffer.data(hh_geom_10_off + 925 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxxz = cbuffer.data(hh_geom_10_off + 926 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyy = cbuffer.data(hh_geom_10_off + 927 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxyz = cbuffer.data(hh_geom_10_off + 928 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxxzz = cbuffer.data(hh_geom_10_off + 929 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyy = cbuffer.data(hh_geom_10_off + 930 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyyz = cbuffer.data(hh_geom_10_off + 931 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxyzz = cbuffer.data(hh_geom_10_off + 932 * ccomps * dcomps);

            auto g_z_0_xxxxz_xxzzz = cbuffer.data(hh_geom_10_off + 933 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyy = cbuffer.data(hh_geom_10_off + 934 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyyz = cbuffer.data(hh_geom_10_off + 935 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyyzz = cbuffer.data(hh_geom_10_off + 936 * ccomps * dcomps);

            auto g_z_0_xxxxz_xyzzz = cbuffer.data(hh_geom_10_off + 937 * ccomps * dcomps);

            auto g_z_0_xxxxz_xzzzz = cbuffer.data(hh_geom_10_off + 938 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyyy = cbuffer.data(hh_geom_10_off + 939 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyyz = cbuffer.data(hh_geom_10_off + 940 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyyzz = cbuffer.data(hh_geom_10_off + 941 * ccomps * dcomps);

            auto g_z_0_xxxxz_yyzzz = cbuffer.data(hh_geom_10_off + 942 * ccomps * dcomps);

            auto g_z_0_xxxxz_yzzzz = cbuffer.data(hh_geom_10_off + 943 * ccomps * dcomps);

            auto g_z_0_xxxxz_zzzzz = cbuffer.data(hh_geom_10_off + 944 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxx = cbuffer.data(hh_geom_10_off + 945 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxy = cbuffer.data(hh_geom_10_off + 946 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxxz = cbuffer.data(hh_geom_10_off + 947 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyy = cbuffer.data(hh_geom_10_off + 948 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxyz = cbuffer.data(hh_geom_10_off + 949 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxxzz = cbuffer.data(hh_geom_10_off + 950 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyy = cbuffer.data(hh_geom_10_off + 951 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyyz = cbuffer.data(hh_geom_10_off + 952 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxyzz = cbuffer.data(hh_geom_10_off + 953 * ccomps * dcomps);

            auto g_z_0_xxxyy_xxzzz = cbuffer.data(hh_geom_10_off + 954 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyy = cbuffer.data(hh_geom_10_off + 955 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyyz = cbuffer.data(hh_geom_10_off + 956 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyyzz = cbuffer.data(hh_geom_10_off + 957 * ccomps * dcomps);

            auto g_z_0_xxxyy_xyzzz = cbuffer.data(hh_geom_10_off + 958 * ccomps * dcomps);

            auto g_z_0_xxxyy_xzzzz = cbuffer.data(hh_geom_10_off + 959 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyyy = cbuffer.data(hh_geom_10_off + 960 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyyz = cbuffer.data(hh_geom_10_off + 961 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyyzz = cbuffer.data(hh_geom_10_off + 962 * ccomps * dcomps);

            auto g_z_0_xxxyy_yyzzz = cbuffer.data(hh_geom_10_off + 963 * ccomps * dcomps);

            auto g_z_0_xxxyy_yzzzz = cbuffer.data(hh_geom_10_off + 964 * ccomps * dcomps);

            auto g_z_0_xxxyy_zzzzz = cbuffer.data(hh_geom_10_off + 965 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxx = cbuffer.data(hh_geom_10_off + 966 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxy = cbuffer.data(hh_geom_10_off + 967 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxxz = cbuffer.data(hh_geom_10_off + 968 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyy = cbuffer.data(hh_geom_10_off + 969 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxyz = cbuffer.data(hh_geom_10_off + 970 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxxzz = cbuffer.data(hh_geom_10_off + 971 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyy = cbuffer.data(hh_geom_10_off + 972 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyyz = cbuffer.data(hh_geom_10_off + 973 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxyzz = cbuffer.data(hh_geom_10_off + 974 * ccomps * dcomps);

            auto g_z_0_xxxyz_xxzzz = cbuffer.data(hh_geom_10_off + 975 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyy = cbuffer.data(hh_geom_10_off + 976 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyyz = cbuffer.data(hh_geom_10_off + 977 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyyzz = cbuffer.data(hh_geom_10_off + 978 * ccomps * dcomps);

            auto g_z_0_xxxyz_xyzzz = cbuffer.data(hh_geom_10_off + 979 * ccomps * dcomps);

            auto g_z_0_xxxyz_xzzzz = cbuffer.data(hh_geom_10_off + 980 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyyy = cbuffer.data(hh_geom_10_off + 981 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyyz = cbuffer.data(hh_geom_10_off + 982 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyyzz = cbuffer.data(hh_geom_10_off + 983 * ccomps * dcomps);

            auto g_z_0_xxxyz_yyzzz = cbuffer.data(hh_geom_10_off + 984 * ccomps * dcomps);

            auto g_z_0_xxxyz_yzzzz = cbuffer.data(hh_geom_10_off + 985 * ccomps * dcomps);

            auto g_z_0_xxxyz_zzzzz = cbuffer.data(hh_geom_10_off + 986 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxx = cbuffer.data(hh_geom_10_off + 987 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxy = cbuffer.data(hh_geom_10_off + 988 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxxz = cbuffer.data(hh_geom_10_off + 989 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyy = cbuffer.data(hh_geom_10_off + 990 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxyz = cbuffer.data(hh_geom_10_off + 991 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxxzz = cbuffer.data(hh_geom_10_off + 992 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyy = cbuffer.data(hh_geom_10_off + 993 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyyz = cbuffer.data(hh_geom_10_off + 994 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxyzz = cbuffer.data(hh_geom_10_off + 995 * ccomps * dcomps);

            auto g_z_0_xxxzz_xxzzz = cbuffer.data(hh_geom_10_off + 996 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyy = cbuffer.data(hh_geom_10_off + 997 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyyz = cbuffer.data(hh_geom_10_off + 998 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyyzz = cbuffer.data(hh_geom_10_off + 999 * ccomps * dcomps);

            auto g_z_0_xxxzz_xyzzz = cbuffer.data(hh_geom_10_off + 1000 * ccomps * dcomps);

            auto g_z_0_xxxzz_xzzzz = cbuffer.data(hh_geom_10_off + 1001 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyyy = cbuffer.data(hh_geom_10_off + 1002 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyyz = cbuffer.data(hh_geom_10_off + 1003 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyyzz = cbuffer.data(hh_geom_10_off + 1004 * ccomps * dcomps);

            auto g_z_0_xxxzz_yyzzz = cbuffer.data(hh_geom_10_off + 1005 * ccomps * dcomps);

            auto g_z_0_xxxzz_yzzzz = cbuffer.data(hh_geom_10_off + 1006 * ccomps * dcomps);

            auto g_z_0_xxxzz_zzzzz = cbuffer.data(hh_geom_10_off + 1007 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxx = cbuffer.data(hh_geom_10_off + 1008 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxy = cbuffer.data(hh_geom_10_off + 1009 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxxz = cbuffer.data(hh_geom_10_off + 1010 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyy = cbuffer.data(hh_geom_10_off + 1011 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxyz = cbuffer.data(hh_geom_10_off + 1012 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxxzz = cbuffer.data(hh_geom_10_off + 1013 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyy = cbuffer.data(hh_geom_10_off + 1014 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyyz = cbuffer.data(hh_geom_10_off + 1015 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxyzz = cbuffer.data(hh_geom_10_off + 1016 * ccomps * dcomps);

            auto g_z_0_xxyyy_xxzzz = cbuffer.data(hh_geom_10_off + 1017 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyy = cbuffer.data(hh_geom_10_off + 1018 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyyz = cbuffer.data(hh_geom_10_off + 1019 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyyzz = cbuffer.data(hh_geom_10_off + 1020 * ccomps * dcomps);

            auto g_z_0_xxyyy_xyzzz = cbuffer.data(hh_geom_10_off + 1021 * ccomps * dcomps);

            auto g_z_0_xxyyy_xzzzz = cbuffer.data(hh_geom_10_off + 1022 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyyy = cbuffer.data(hh_geom_10_off + 1023 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyyz = cbuffer.data(hh_geom_10_off + 1024 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyyzz = cbuffer.data(hh_geom_10_off + 1025 * ccomps * dcomps);

            auto g_z_0_xxyyy_yyzzz = cbuffer.data(hh_geom_10_off + 1026 * ccomps * dcomps);

            auto g_z_0_xxyyy_yzzzz = cbuffer.data(hh_geom_10_off + 1027 * ccomps * dcomps);

            auto g_z_0_xxyyy_zzzzz = cbuffer.data(hh_geom_10_off + 1028 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxx = cbuffer.data(hh_geom_10_off + 1029 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxy = cbuffer.data(hh_geom_10_off + 1030 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxxz = cbuffer.data(hh_geom_10_off + 1031 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyy = cbuffer.data(hh_geom_10_off + 1032 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxyz = cbuffer.data(hh_geom_10_off + 1033 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxxzz = cbuffer.data(hh_geom_10_off + 1034 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyy = cbuffer.data(hh_geom_10_off + 1035 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyyz = cbuffer.data(hh_geom_10_off + 1036 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxyzz = cbuffer.data(hh_geom_10_off + 1037 * ccomps * dcomps);

            auto g_z_0_xxyyz_xxzzz = cbuffer.data(hh_geom_10_off + 1038 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyy = cbuffer.data(hh_geom_10_off + 1039 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyyz = cbuffer.data(hh_geom_10_off + 1040 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyyzz = cbuffer.data(hh_geom_10_off + 1041 * ccomps * dcomps);

            auto g_z_0_xxyyz_xyzzz = cbuffer.data(hh_geom_10_off + 1042 * ccomps * dcomps);

            auto g_z_0_xxyyz_xzzzz = cbuffer.data(hh_geom_10_off + 1043 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyyy = cbuffer.data(hh_geom_10_off + 1044 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyyz = cbuffer.data(hh_geom_10_off + 1045 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyyzz = cbuffer.data(hh_geom_10_off + 1046 * ccomps * dcomps);

            auto g_z_0_xxyyz_yyzzz = cbuffer.data(hh_geom_10_off + 1047 * ccomps * dcomps);

            auto g_z_0_xxyyz_yzzzz = cbuffer.data(hh_geom_10_off + 1048 * ccomps * dcomps);

            auto g_z_0_xxyyz_zzzzz = cbuffer.data(hh_geom_10_off + 1049 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxx = cbuffer.data(hh_geom_10_off + 1050 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxy = cbuffer.data(hh_geom_10_off + 1051 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxxz = cbuffer.data(hh_geom_10_off + 1052 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyy = cbuffer.data(hh_geom_10_off + 1053 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxyz = cbuffer.data(hh_geom_10_off + 1054 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxxzz = cbuffer.data(hh_geom_10_off + 1055 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyy = cbuffer.data(hh_geom_10_off + 1056 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyyz = cbuffer.data(hh_geom_10_off + 1057 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxyzz = cbuffer.data(hh_geom_10_off + 1058 * ccomps * dcomps);

            auto g_z_0_xxyzz_xxzzz = cbuffer.data(hh_geom_10_off + 1059 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyy = cbuffer.data(hh_geom_10_off + 1060 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyyz = cbuffer.data(hh_geom_10_off + 1061 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyyzz = cbuffer.data(hh_geom_10_off + 1062 * ccomps * dcomps);

            auto g_z_0_xxyzz_xyzzz = cbuffer.data(hh_geom_10_off + 1063 * ccomps * dcomps);

            auto g_z_0_xxyzz_xzzzz = cbuffer.data(hh_geom_10_off + 1064 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyyy = cbuffer.data(hh_geom_10_off + 1065 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyyz = cbuffer.data(hh_geom_10_off + 1066 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyyzz = cbuffer.data(hh_geom_10_off + 1067 * ccomps * dcomps);

            auto g_z_0_xxyzz_yyzzz = cbuffer.data(hh_geom_10_off + 1068 * ccomps * dcomps);

            auto g_z_0_xxyzz_yzzzz = cbuffer.data(hh_geom_10_off + 1069 * ccomps * dcomps);

            auto g_z_0_xxyzz_zzzzz = cbuffer.data(hh_geom_10_off + 1070 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1071 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1072 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1073 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1074 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1075 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1076 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1077 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1078 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1079 * ccomps * dcomps);

            auto g_z_0_xxzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1080 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1081 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1082 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1083 * ccomps * dcomps);

            auto g_z_0_xxzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1084 * ccomps * dcomps);

            auto g_z_0_xxzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1085 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1086 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1087 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1088 * ccomps * dcomps);

            auto g_z_0_xxzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1089 * ccomps * dcomps);

            auto g_z_0_xxzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1090 * ccomps * dcomps);

            auto g_z_0_xxzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1091 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 1092 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 1093 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 1094 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 1095 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 1096 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 1097 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 1098 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 1099 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 1100 * ccomps * dcomps);

            auto g_z_0_xyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 1101 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 1102 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 1103 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 1104 * ccomps * dcomps);

            auto g_z_0_xyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 1105 * ccomps * dcomps);

            auto g_z_0_xyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 1106 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 1107 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 1108 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 1109 * ccomps * dcomps);

            auto g_z_0_xyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 1110 * ccomps * dcomps);

            auto g_z_0_xyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 1111 * ccomps * dcomps);

            auto g_z_0_xyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 1112 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 1113 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 1114 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 1115 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 1116 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 1117 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 1118 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 1119 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 1120 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 1121 * ccomps * dcomps);

            auto g_z_0_xyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 1122 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 1123 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 1124 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 1125 * ccomps * dcomps);

            auto g_z_0_xyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 1126 * ccomps * dcomps);

            auto g_z_0_xyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 1127 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 1128 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 1129 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 1130 * ccomps * dcomps);

            auto g_z_0_xyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 1131 * ccomps * dcomps);

            auto g_z_0_xyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 1132 * ccomps * dcomps);

            auto g_z_0_xyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 1133 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 1134 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 1135 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 1136 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 1137 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 1138 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 1139 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 1140 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 1141 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 1142 * ccomps * dcomps);

            auto g_z_0_xyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 1143 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 1144 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 1145 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 1146 * ccomps * dcomps);

            auto g_z_0_xyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 1147 * ccomps * dcomps);

            auto g_z_0_xyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 1148 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 1149 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 1150 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 1151 * ccomps * dcomps);

            auto g_z_0_xyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 1152 * ccomps * dcomps);

            auto g_z_0_xyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 1153 * ccomps * dcomps);

            auto g_z_0_xyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 1154 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1155 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1156 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1157 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1158 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1159 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1160 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1161 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1162 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1163 * ccomps * dcomps);

            auto g_z_0_xyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1164 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1165 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1166 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1167 * ccomps * dcomps);

            auto g_z_0_xyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1168 * ccomps * dcomps);

            auto g_z_0_xyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1169 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1170 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1171 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1172 * ccomps * dcomps);

            auto g_z_0_xyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1173 * ccomps * dcomps);

            auto g_z_0_xyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1174 * ccomps * dcomps);

            auto g_z_0_xyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1175 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1176 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1177 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1178 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1179 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1180 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1181 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1182 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1183 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1184 * ccomps * dcomps);

            auto g_z_0_xzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1185 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1186 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1187 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1188 * ccomps * dcomps);

            auto g_z_0_xzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1189 * ccomps * dcomps);

            auto g_z_0_xzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1190 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1191 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1192 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1193 * ccomps * dcomps);

            auto g_z_0_xzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1194 * ccomps * dcomps);

            auto g_z_0_xzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1195 * ccomps * dcomps);

            auto g_z_0_xzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1196 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxx = cbuffer.data(hh_geom_10_off + 1197 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxy = cbuffer.data(hh_geom_10_off + 1198 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxxz = cbuffer.data(hh_geom_10_off + 1199 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyy = cbuffer.data(hh_geom_10_off + 1200 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxyz = cbuffer.data(hh_geom_10_off + 1201 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxxzz = cbuffer.data(hh_geom_10_off + 1202 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyy = cbuffer.data(hh_geom_10_off + 1203 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyyz = cbuffer.data(hh_geom_10_off + 1204 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxyzz = cbuffer.data(hh_geom_10_off + 1205 * ccomps * dcomps);

            auto g_z_0_yyyyy_xxzzz = cbuffer.data(hh_geom_10_off + 1206 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyy = cbuffer.data(hh_geom_10_off + 1207 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyyz = cbuffer.data(hh_geom_10_off + 1208 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyyzz = cbuffer.data(hh_geom_10_off + 1209 * ccomps * dcomps);

            auto g_z_0_yyyyy_xyzzz = cbuffer.data(hh_geom_10_off + 1210 * ccomps * dcomps);

            auto g_z_0_yyyyy_xzzzz = cbuffer.data(hh_geom_10_off + 1211 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyy = cbuffer.data(hh_geom_10_off + 1212 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyyz = cbuffer.data(hh_geom_10_off + 1213 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyyzz = cbuffer.data(hh_geom_10_off + 1214 * ccomps * dcomps);

            auto g_z_0_yyyyy_yyzzz = cbuffer.data(hh_geom_10_off + 1215 * ccomps * dcomps);

            auto g_z_0_yyyyy_yzzzz = cbuffer.data(hh_geom_10_off + 1216 * ccomps * dcomps);

            auto g_z_0_yyyyy_zzzzz = cbuffer.data(hh_geom_10_off + 1217 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxx = cbuffer.data(hh_geom_10_off + 1218 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxy = cbuffer.data(hh_geom_10_off + 1219 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxxz = cbuffer.data(hh_geom_10_off + 1220 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyy = cbuffer.data(hh_geom_10_off + 1221 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxyz = cbuffer.data(hh_geom_10_off + 1222 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxxzz = cbuffer.data(hh_geom_10_off + 1223 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyy = cbuffer.data(hh_geom_10_off + 1224 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyyz = cbuffer.data(hh_geom_10_off + 1225 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxyzz = cbuffer.data(hh_geom_10_off + 1226 * ccomps * dcomps);

            auto g_z_0_yyyyz_xxzzz = cbuffer.data(hh_geom_10_off + 1227 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyy = cbuffer.data(hh_geom_10_off + 1228 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyyz = cbuffer.data(hh_geom_10_off + 1229 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyyzz = cbuffer.data(hh_geom_10_off + 1230 * ccomps * dcomps);

            auto g_z_0_yyyyz_xyzzz = cbuffer.data(hh_geom_10_off + 1231 * ccomps * dcomps);

            auto g_z_0_yyyyz_xzzzz = cbuffer.data(hh_geom_10_off + 1232 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyy = cbuffer.data(hh_geom_10_off + 1233 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyyz = cbuffer.data(hh_geom_10_off + 1234 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyyzz = cbuffer.data(hh_geom_10_off + 1235 * ccomps * dcomps);

            auto g_z_0_yyyyz_yyzzz = cbuffer.data(hh_geom_10_off + 1236 * ccomps * dcomps);

            auto g_z_0_yyyyz_yzzzz = cbuffer.data(hh_geom_10_off + 1237 * ccomps * dcomps);

            auto g_z_0_yyyyz_zzzzz = cbuffer.data(hh_geom_10_off + 1238 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxx = cbuffer.data(hh_geom_10_off + 1239 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxy = cbuffer.data(hh_geom_10_off + 1240 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxxz = cbuffer.data(hh_geom_10_off + 1241 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyy = cbuffer.data(hh_geom_10_off + 1242 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxyz = cbuffer.data(hh_geom_10_off + 1243 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxxzz = cbuffer.data(hh_geom_10_off + 1244 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyy = cbuffer.data(hh_geom_10_off + 1245 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyyz = cbuffer.data(hh_geom_10_off + 1246 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxyzz = cbuffer.data(hh_geom_10_off + 1247 * ccomps * dcomps);

            auto g_z_0_yyyzz_xxzzz = cbuffer.data(hh_geom_10_off + 1248 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyy = cbuffer.data(hh_geom_10_off + 1249 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyyz = cbuffer.data(hh_geom_10_off + 1250 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyyzz = cbuffer.data(hh_geom_10_off + 1251 * ccomps * dcomps);

            auto g_z_0_yyyzz_xyzzz = cbuffer.data(hh_geom_10_off + 1252 * ccomps * dcomps);

            auto g_z_0_yyyzz_xzzzz = cbuffer.data(hh_geom_10_off + 1253 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyy = cbuffer.data(hh_geom_10_off + 1254 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyyz = cbuffer.data(hh_geom_10_off + 1255 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyyzz = cbuffer.data(hh_geom_10_off + 1256 * ccomps * dcomps);

            auto g_z_0_yyyzz_yyzzz = cbuffer.data(hh_geom_10_off + 1257 * ccomps * dcomps);

            auto g_z_0_yyyzz_yzzzz = cbuffer.data(hh_geom_10_off + 1258 * ccomps * dcomps);

            auto g_z_0_yyyzz_zzzzz = cbuffer.data(hh_geom_10_off + 1259 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1260 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1261 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1262 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1263 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1264 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1265 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1266 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1267 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1268 * ccomps * dcomps);

            auto g_z_0_yyzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1269 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1270 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1271 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1272 * ccomps * dcomps);

            auto g_z_0_yyzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1273 * ccomps * dcomps);

            auto g_z_0_yyzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1274 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1275 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1276 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1277 * ccomps * dcomps);

            auto g_z_0_yyzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1278 * ccomps * dcomps);

            auto g_z_0_yyzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1279 * ccomps * dcomps);

            auto g_z_0_yyzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1280 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1281 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1282 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1283 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1284 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1285 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1286 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1287 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1288 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1289 * ccomps * dcomps);

            auto g_z_0_yzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1290 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1291 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1292 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1293 * ccomps * dcomps);

            auto g_z_0_yzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1294 * ccomps * dcomps);

            auto g_z_0_yzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1295 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1296 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1297 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1298 * ccomps * dcomps);

            auto g_z_0_yzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1299 * ccomps * dcomps);

            auto g_z_0_yzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1300 * ccomps * dcomps);

            auto g_z_0_yzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1301 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxx = cbuffer.data(hh_geom_10_off + 1302 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxy = cbuffer.data(hh_geom_10_off + 1303 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxxz = cbuffer.data(hh_geom_10_off + 1304 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyy = cbuffer.data(hh_geom_10_off + 1305 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxyz = cbuffer.data(hh_geom_10_off + 1306 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxxzz = cbuffer.data(hh_geom_10_off + 1307 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyy = cbuffer.data(hh_geom_10_off + 1308 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyyz = cbuffer.data(hh_geom_10_off + 1309 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxyzz = cbuffer.data(hh_geom_10_off + 1310 * ccomps * dcomps);

            auto g_z_0_zzzzz_xxzzz = cbuffer.data(hh_geom_10_off + 1311 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyy = cbuffer.data(hh_geom_10_off + 1312 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyyz = cbuffer.data(hh_geom_10_off + 1313 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyyzz = cbuffer.data(hh_geom_10_off + 1314 * ccomps * dcomps);

            auto g_z_0_zzzzz_xyzzz = cbuffer.data(hh_geom_10_off + 1315 * ccomps * dcomps);

            auto g_z_0_zzzzz_xzzzz = cbuffer.data(hh_geom_10_off + 1316 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyy = cbuffer.data(hh_geom_10_off + 1317 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyyz = cbuffer.data(hh_geom_10_off + 1318 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyyzz = cbuffer.data(hh_geom_10_off + 1319 * ccomps * dcomps);

            auto g_z_0_zzzzz_yyzzz = cbuffer.data(hh_geom_10_off + 1320 * ccomps * dcomps);

            auto g_z_0_zzzzz_yzzzz = cbuffer.data(hh_geom_10_off + 1321 * ccomps * dcomps);

            auto g_z_0_zzzzz_zzzzz = cbuffer.data(hh_geom_10_off + 1322 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_igxx

            const auto ig_geom_10_off = idx_geom_10_igxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxx_xxxx = cbuffer.data(ig_geom_10_off + 0 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxy = cbuffer.data(ig_geom_10_off + 1 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxxz = cbuffer.data(ig_geom_10_off + 2 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxyy = cbuffer.data(ig_geom_10_off + 3 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxyz = cbuffer.data(ig_geom_10_off + 4 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xxzz = cbuffer.data(ig_geom_10_off + 5 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyyy = cbuffer.data(ig_geom_10_off + 6 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyyz = cbuffer.data(ig_geom_10_off + 7 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xyzz = cbuffer.data(ig_geom_10_off + 8 * ccomps * dcomps);

            auto g_x_0_xxxxxx_xzzz = cbuffer.data(ig_geom_10_off + 9 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyyy = cbuffer.data(ig_geom_10_off + 10 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyyz = cbuffer.data(ig_geom_10_off + 11 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yyzz = cbuffer.data(ig_geom_10_off + 12 * ccomps * dcomps);

            auto g_x_0_xxxxxx_yzzz = cbuffer.data(ig_geom_10_off + 13 * ccomps * dcomps);

            auto g_x_0_xxxxxx_zzzz = cbuffer.data(ig_geom_10_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxxx, g_x_0_xxxxx_xxxxx, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_yyyy, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_zzzz, g_x_0_xxxxxx_xxxx, g_x_0_xxxxxx_xxxy, g_x_0_xxxxxx_xxxz, g_x_0_xxxxxx_xxyy, g_x_0_xxxxxx_xxyz, g_x_0_xxxxxx_xxzz, g_x_0_xxxxxx_xyyy, g_x_0_xxxxxx_xyyz, g_x_0_xxxxxx_xyzz, g_x_0_xxxxxx_xzzz, g_x_0_xxxxxx_yyyy, g_x_0_xxxxxx_yyyz, g_x_0_xxxxxx_yyzz, g_x_0_xxxxxx_yzzz, g_x_0_xxxxxx_zzzz, g_xxxxx_xxxx, g_xxxxx_xxxy, g_xxxxx_xxxz, g_xxxxx_xxyy, g_xxxxx_xxyz, g_xxxxx_xxzz, g_xxxxx_xyyy, g_xxxxx_xyyz, g_xxxxx_xyzz, g_xxxxx_xzzz, g_xxxxx_yyyy, g_xxxxx_yyyz, g_xxxxx_yyzz, g_xxxxx_yzzz, g_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxx_xxxx[k] = -g_xxxxx_xxxx[k] - g_x_0_xxxxx_xxxx[k] * ab_x + g_x_0_xxxxx_xxxxx[k];

                g_x_0_xxxxxx_xxxy[k] = -g_xxxxx_xxxy[k] - g_x_0_xxxxx_xxxy[k] * ab_x + g_x_0_xxxxx_xxxxy[k];

                g_x_0_xxxxxx_xxxz[k] = -g_xxxxx_xxxz[k] - g_x_0_xxxxx_xxxz[k] * ab_x + g_x_0_xxxxx_xxxxz[k];

                g_x_0_xxxxxx_xxyy[k] = -g_xxxxx_xxyy[k] - g_x_0_xxxxx_xxyy[k] * ab_x + g_x_0_xxxxx_xxxyy[k];

                g_x_0_xxxxxx_xxyz[k] = -g_xxxxx_xxyz[k] - g_x_0_xxxxx_xxyz[k] * ab_x + g_x_0_xxxxx_xxxyz[k];

                g_x_0_xxxxxx_xxzz[k] = -g_xxxxx_xxzz[k] - g_x_0_xxxxx_xxzz[k] * ab_x + g_x_0_xxxxx_xxxzz[k];

                g_x_0_xxxxxx_xyyy[k] = -g_xxxxx_xyyy[k] - g_x_0_xxxxx_xyyy[k] * ab_x + g_x_0_xxxxx_xxyyy[k];

                g_x_0_xxxxxx_xyyz[k] = -g_xxxxx_xyyz[k] - g_x_0_xxxxx_xyyz[k] * ab_x + g_x_0_xxxxx_xxyyz[k];

                g_x_0_xxxxxx_xyzz[k] = -g_xxxxx_xyzz[k] - g_x_0_xxxxx_xyzz[k] * ab_x + g_x_0_xxxxx_xxyzz[k];

                g_x_0_xxxxxx_xzzz[k] = -g_xxxxx_xzzz[k] - g_x_0_xxxxx_xzzz[k] * ab_x + g_x_0_xxxxx_xxzzz[k];

                g_x_0_xxxxxx_yyyy[k] = -g_xxxxx_yyyy[k] - g_x_0_xxxxx_yyyy[k] * ab_x + g_x_0_xxxxx_xyyyy[k];

                g_x_0_xxxxxx_yyyz[k] = -g_xxxxx_yyyz[k] - g_x_0_xxxxx_yyyz[k] * ab_x + g_x_0_xxxxx_xyyyz[k];

                g_x_0_xxxxxx_yyzz[k] = -g_xxxxx_yyzz[k] - g_x_0_xxxxx_yyzz[k] * ab_x + g_x_0_xxxxx_xyyzz[k];

                g_x_0_xxxxxx_yzzz[k] = -g_xxxxx_yzzz[k] - g_x_0_xxxxx_yzzz[k] * ab_x + g_x_0_xxxxx_xyzzz[k];

                g_x_0_xxxxxx_zzzz[k] = -g_xxxxx_zzzz[k] - g_x_0_xxxxx_zzzz[k] * ab_x + g_x_0_xxxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxy_xxxx = cbuffer.data(ig_geom_10_off + 15 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxy = cbuffer.data(ig_geom_10_off + 16 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxxz = cbuffer.data(ig_geom_10_off + 17 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxyy = cbuffer.data(ig_geom_10_off + 18 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxyz = cbuffer.data(ig_geom_10_off + 19 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xxzz = cbuffer.data(ig_geom_10_off + 20 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyyy = cbuffer.data(ig_geom_10_off + 21 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyyz = cbuffer.data(ig_geom_10_off + 22 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xyzz = cbuffer.data(ig_geom_10_off + 23 * ccomps * dcomps);

            auto g_x_0_xxxxxy_xzzz = cbuffer.data(ig_geom_10_off + 24 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyyy = cbuffer.data(ig_geom_10_off + 25 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyyz = cbuffer.data(ig_geom_10_off + 26 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yyzz = cbuffer.data(ig_geom_10_off + 27 * ccomps * dcomps);

            auto g_x_0_xxxxxy_yzzz = cbuffer.data(ig_geom_10_off + 28 * ccomps * dcomps);

            auto g_x_0_xxxxxy_zzzz = cbuffer.data(ig_geom_10_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxxx, g_x_0_xxxxx_xxxxy, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxxyy, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyyy, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyyy, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_yyyy, g_x_0_xxxxx_yyyyy, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_zzzz, g_x_0_xxxxxy_xxxx, g_x_0_xxxxxy_xxxy, g_x_0_xxxxxy_xxxz, g_x_0_xxxxxy_xxyy, g_x_0_xxxxxy_xxyz, g_x_0_xxxxxy_xxzz, g_x_0_xxxxxy_xyyy, g_x_0_xxxxxy_xyyz, g_x_0_xxxxxy_xyzz, g_x_0_xxxxxy_xzzz, g_x_0_xxxxxy_yyyy, g_x_0_xxxxxy_yyyz, g_x_0_xxxxxy_yyzz, g_x_0_xxxxxy_yzzz, g_x_0_xxxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxy_xxxx[k] = -g_x_0_xxxxx_xxxx[k] * ab_y + g_x_0_xxxxx_xxxxy[k];

                g_x_0_xxxxxy_xxxy[k] = -g_x_0_xxxxx_xxxy[k] * ab_y + g_x_0_xxxxx_xxxyy[k];

                g_x_0_xxxxxy_xxxz[k] = -g_x_0_xxxxx_xxxz[k] * ab_y + g_x_0_xxxxx_xxxyz[k];

                g_x_0_xxxxxy_xxyy[k] = -g_x_0_xxxxx_xxyy[k] * ab_y + g_x_0_xxxxx_xxyyy[k];

                g_x_0_xxxxxy_xxyz[k] = -g_x_0_xxxxx_xxyz[k] * ab_y + g_x_0_xxxxx_xxyyz[k];

                g_x_0_xxxxxy_xxzz[k] = -g_x_0_xxxxx_xxzz[k] * ab_y + g_x_0_xxxxx_xxyzz[k];

                g_x_0_xxxxxy_xyyy[k] = -g_x_0_xxxxx_xyyy[k] * ab_y + g_x_0_xxxxx_xyyyy[k];

                g_x_0_xxxxxy_xyyz[k] = -g_x_0_xxxxx_xyyz[k] * ab_y + g_x_0_xxxxx_xyyyz[k];

                g_x_0_xxxxxy_xyzz[k] = -g_x_0_xxxxx_xyzz[k] * ab_y + g_x_0_xxxxx_xyyzz[k];

                g_x_0_xxxxxy_xzzz[k] = -g_x_0_xxxxx_xzzz[k] * ab_y + g_x_0_xxxxx_xyzzz[k];

                g_x_0_xxxxxy_yyyy[k] = -g_x_0_xxxxx_yyyy[k] * ab_y + g_x_0_xxxxx_yyyyy[k];

                g_x_0_xxxxxy_yyyz[k] = -g_x_0_xxxxx_yyyz[k] * ab_y + g_x_0_xxxxx_yyyyz[k];

                g_x_0_xxxxxy_yyzz[k] = -g_x_0_xxxxx_yyzz[k] * ab_y + g_x_0_xxxxx_yyyzz[k];

                g_x_0_xxxxxy_yzzz[k] = -g_x_0_xxxxx_yzzz[k] * ab_y + g_x_0_xxxxx_yyzzz[k];

                g_x_0_xxxxxy_zzzz[k] = -g_x_0_xxxxx_zzzz[k] * ab_y + g_x_0_xxxxx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxxz_xxxx = cbuffer.data(ig_geom_10_off + 30 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxy = cbuffer.data(ig_geom_10_off + 31 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxxz = cbuffer.data(ig_geom_10_off + 32 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxyy = cbuffer.data(ig_geom_10_off + 33 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxyz = cbuffer.data(ig_geom_10_off + 34 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xxzz = cbuffer.data(ig_geom_10_off + 35 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyyy = cbuffer.data(ig_geom_10_off + 36 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyyz = cbuffer.data(ig_geom_10_off + 37 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xyzz = cbuffer.data(ig_geom_10_off + 38 * ccomps * dcomps);

            auto g_x_0_xxxxxz_xzzz = cbuffer.data(ig_geom_10_off + 39 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyyy = cbuffer.data(ig_geom_10_off + 40 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyyz = cbuffer.data(ig_geom_10_off + 41 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yyzz = cbuffer.data(ig_geom_10_off + 42 * ccomps * dcomps);

            auto g_x_0_xxxxxz_yzzz = cbuffer.data(ig_geom_10_off + 43 * ccomps * dcomps);

            auto g_x_0_xxxxxz_zzzz = cbuffer.data(ig_geom_10_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxx_xxxx, g_x_0_xxxxx_xxxxz, g_x_0_xxxxx_xxxy, g_x_0_xxxxx_xxxyz, g_x_0_xxxxx_xxxz, g_x_0_xxxxx_xxxzz, g_x_0_xxxxx_xxyy, g_x_0_xxxxx_xxyyz, g_x_0_xxxxx_xxyz, g_x_0_xxxxx_xxyzz, g_x_0_xxxxx_xxzz, g_x_0_xxxxx_xxzzz, g_x_0_xxxxx_xyyy, g_x_0_xxxxx_xyyyz, g_x_0_xxxxx_xyyz, g_x_0_xxxxx_xyyzz, g_x_0_xxxxx_xyzz, g_x_0_xxxxx_xyzzz, g_x_0_xxxxx_xzzz, g_x_0_xxxxx_xzzzz, g_x_0_xxxxx_yyyy, g_x_0_xxxxx_yyyyz, g_x_0_xxxxx_yyyz, g_x_0_xxxxx_yyyzz, g_x_0_xxxxx_yyzz, g_x_0_xxxxx_yyzzz, g_x_0_xxxxx_yzzz, g_x_0_xxxxx_yzzzz, g_x_0_xxxxx_zzzz, g_x_0_xxxxx_zzzzz, g_x_0_xxxxxz_xxxx, g_x_0_xxxxxz_xxxy, g_x_0_xxxxxz_xxxz, g_x_0_xxxxxz_xxyy, g_x_0_xxxxxz_xxyz, g_x_0_xxxxxz_xxzz, g_x_0_xxxxxz_xyyy, g_x_0_xxxxxz_xyyz, g_x_0_xxxxxz_xyzz, g_x_0_xxxxxz_xzzz, g_x_0_xxxxxz_yyyy, g_x_0_xxxxxz_yyyz, g_x_0_xxxxxz_yyzz, g_x_0_xxxxxz_yzzz, g_x_0_xxxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxxz_xxxx[k] = -g_x_0_xxxxx_xxxx[k] * ab_z + g_x_0_xxxxx_xxxxz[k];

                g_x_0_xxxxxz_xxxy[k] = -g_x_0_xxxxx_xxxy[k] * ab_z + g_x_0_xxxxx_xxxyz[k];

                g_x_0_xxxxxz_xxxz[k] = -g_x_0_xxxxx_xxxz[k] * ab_z + g_x_0_xxxxx_xxxzz[k];

                g_x_0_xxxxxz_xxyy[k] = -g_x_0_xxxxx_xxyy[k] * ab_z + g_x_0_xxxxx_xxyyz[k];

                g_x_0_xxxxxz_xxyz[k] = -g_x_0_xxxxx_xxyz[k] * ab_z + g_x_0_xxxxx_xxyzz[k];

                g_x_0_xxxxxz_xxzz[k] = -g_x_0_xxxxx_xxzz[k] * ab_z + g_x_0_xxxxx_xxzzz[k];

                g_x_0_xxxxxz_xyyy[k] = -g_x_0_xxxxx_xyyy[k] * ab_z + g_x_0_xxxxx_xyyyz[k];

                g_x_0_xxxxxz_xyyz[k] = -g_x_0_xxxxx_xyyz[k] * ab_z + g_x_0_xxxxx_xyyzz[k];

                g_x_0_xxxxxz_xyzz[k] = -g_x_0_xxxxx_xyzz[k] * ab_z + g_x_0_xxxxx_xyzzz[k];

                g_x_0_xxxxxz_xzzz[k] = -g_x_0_xxxxx_xzzz[k] * ab_z + g_x_0_xxxxx_xzzzz[k];

                g_x_0_xxxxxz_yyyy[k] = -g_x_0_xxxxx_yyyy[k] * ab_z + g_x_0_xxxxx_yyyyz[k];

                g_x_0_xxxxxz_yyyz[k] = -g_x_0_xxxxx_yyyz[k] * ab_z + g_x_0_xxxxx_yyyzz[k];

                g_x_0_xxxxxz_yyzz[k] = -g_x_0_xxxxx_yyzz[k] * ab_z + g_x_0_xxxxx_yyzzz[k];

                g_x_0_xxxxxz_yzzz[k] = -g_x_0_xxxxx_yzzz[k] * ab_z + g_x_0_xxxxx_yzzzz[k];

                g_x_0_xxxxxz_zzzz[k] = -g_x_0_xxxxx_zzzz[k] * ab_z + g_x_0_xxxxx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyy_xxxx = cbuffer.data(ig_geom_10_off + 45 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxy = cbuffer.data(ig_geom_10_off + 46 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxxz = cbuffer.data(ig_geom_10_off + 47 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxyy = cbuffer.data(ig_geom_10_off + 48 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxyz = cbuffer.data(ig_geom_10_off + 49 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xxzz = cbuffer.data(ig_geom_10_off + 50 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyyy = cbuffer.data(ig_geom_10_off + 51 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyyz = cbuffer.data(ig_geom_10_off + 52 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xyzz = cbuffer.data(ig_geom_10_off + 53 * ccomps * dcomps);

            auto g_x_0_xxxxyy_xzzz = cbuffer.data(ig_geom_10_off + 54 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyyy = cbuffer.data(ig_geom_10_off + 55 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyyz = cbuffer.data(ig_geom_10_off + 56 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yyzz = cbuffer.data(ig_geom_10_off + 57 * ccomps * dcomps);

            auto g_x_0_xxxxyy_yzzz = cbuffer.data(ig_geom_10_off + 58 * ccomps * dcomps);

            auto g_x_0_xxxxyy_zzzz = cbuffer.data(ig_geom_10_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxy_xxxx, g_x_0_xxxxy_xxxxy, g_x_0_xxxxy_xxxy, g_x_0_xxxxy_xxxyy, g_x_0_xxxxy_xxxyz, g_x_0_xxxxy_xxxz, g_x_0_xxxxy_xxyy, g_x_0_xxxxy_xxyyy, g_x_0_xxxxy_xxyyz, g_x_0_xxxxy_xxyz, g_x_0_xxxxy_xxyzz, g_x_0_xxxxy_xxzz, g_x_0_xxxxy_xyyy, g_x_0_xxxxy_xyyyy, g_x_0_xxxxy_xyyyz, g_x_0_xxxxy_xyyz, g_x_0_xxxxy_xyyzz, g_x_0_xxxxy_xyzz, g_x_0_xxxxy_xyzzz, g_x_0_xxxxy_xzzz, g_x_0_xxxxy_yyyy, g_x_0_xxxxy_yyyyy, g_x_0_xxxxy_yyyyz, g_x_0_xxxxy_yyyz, g_x_0_xxxxy_yyyzz, g_x_0_xxxxy_yyzz, g_x_0_xxxxy_yyzzz, g_x_0_xxxxy_yzzz, g_x_0_xxxxy_yzzzz, g_x_0_xxxxy_zzzz, g_x_0_xxxxyy_xxxx, g_x_0_xxxxyy_xxxy, g_x_0_xxxxyy_xxxz, g_x_0_xxxxyy_xxyy, g_x_0_xxxxyy_xxyz, g_x_0_xxxxyy_xxzz, g_x_0_xxxxyy_xyyy, g_x_0_xxxxyy_xyyz, g_x_0_xxxxyy_xyzz, g_x_0_xxxxyy_xzzz, g_x_0_xxxxyy_yyyy, g_x_0_xxxxyy_yyyz, g_x_0_xxxxyy_yyzz, g_x_0_xxxxyy_yzzz, g_x_0_xxxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyy_xxxx[k] = -g_x_0_xxxxy_xxxx[k] * ab_y + g_x_0_xxxxy_xxxxy[k];

                g_x_0_xxxxyy_xxxy[k] = -g_x_0_xxxxy_xxxy[k] * ab_y + g_x_0_xxxxy_xxxyy[k];

                g_x_0_xxxxyy_xxxz[k] = -g_x_0_xxxxy_xxxz[k] * ab_y + g_x_0_xxxxy_xxxyz[k];

                g_x_0_xxxxyy_xxyy[k] = -g_x_0_xxxxy_xxyy[k] * ab_y + g_x_0_xxxxy_xxyyy[k];

                g_x_0_xxxxyy_xxyz[k] = -g_x_0_xxxxy_xxyz[k] * ab_y + g_x_0_xxxxy_xxyyz[k];

                g_x_0_xxxxyy_xxzz[k] = -g_x_0_xxxxy_xxzz[k] * ab_y + g_x_0_xxxxy_xxyzz[k];

                g_x_0_xxxxyy_xyyy[k] = -g_x_0_xxxxy_xyyy[k] * ab_y + g_x_0_xxxxy_xyyyy[k];

                g_x_0_xxxxyy_xyyz[k] = -g_x_0_xxxxy_xyyz[k] * ab_y + g_x_0_xxxxy_xyyyz[k];

                g_x_0_xxxxyy_xyzz[k] = -g_x_0_xxxxy_xyzz[k] * ab_y + g_x_0_xxxxy_xyyzz[k];

                g_x_0_xxxxyy_xzzz[k] = -g_x_0_xxxxy_xzzz[k] * ab_y + g_x_0_xxxxy_xyzzz[k];

                g_x_0_xxxxyy_yyyy[k] = -g_x_0_xxxxy_yyyy[k] * ab_y + g_x_0_xxxxy_yyyyy[k];

                g_x_0_xxxxyy_yyyz[k] = -g_x_0_xxxxy_yyyz[k] * ab_y + g_x_0_xxxxy_yyyyz[k];

                g_x_0_xxxxyy_yyzz[k] = -g_x_0_xxxxy_yyzz[k] * ab_y + g_x_0_xxxxy_yyyzz[k];

                g_x_0_xxxxyy_yzzz[k] = -g_x_0_xxxxy_yzzz[k] * ab_y + g_x_0_xxxxy_yyzzz[k];

                g_x_0_xxxxyy_zzzz[k] = -g_x_0_xxxxy_zzzz[k] * ab_y + g_x_0_xxxxy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxyz_xxxx = cbuffer.data(ig_geom_10_off + 60 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxy = cbuffer.data(ig_geom_10_off + 61 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxxz = cbuffer.data(ig_geom_10_off + 62 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxyy = cbuffer.data(ig_geom_10_off + 63 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxyz = cbuffer.data(ig_geom_10_off + 64 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xxzz = cbuffer.data(ig_geom_10_off + 65 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyyy = cbuffer.data(ig_geom_10_off + 66 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyyz = cbuffer.data(ig_geom_10_off + 67 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xyzz = cbuffer.data(ig_geom_10_off + 68 * ccomps * dcomps);

            auto g_x_0_xxxxyz_xzzz = cbuffer.data(ig_geom_10_off + 69 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyyy = cbuffer.data(ig_geom_10_off + 70 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyyz = cbuffer.data(ig_geom_10_off + 71 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yyzz = cbuffer.data(ig_geom_10_off + 72 * ccomps * dcomps);

            auto g_x_0_xxxxyz_yzzz = cbuffer.data(ig_geom_10_off + 73 * ccomps * dcomps);

            auto g_x_0_xxxxyz_zzzz = cbuffer.data(ig_geom_10_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxyz_xxxx, g_x_0_xxxxyz_xxxy, g_x_0_xxxxyz_xxxz, g_x_0_xxxxyz_xxyy, g_x_0_xxxxyz_xxyz, g_x_0_xxxxyz_xxzz, g_x_0_xxxxyz_xyyy, g_x_0_xxxxyz_xyyz, g_x_0_xxxxyz_xyzz, g_x_0_xxxxyz_xzzz, g_x_0_xxxxyz_yyyy, g_x_0_xxxxyz_yyyz, g_x_0_xxxxyz_yyzz, g_x_0_xxxxyz_yzzz, g_x_0_xxxxyz_zzzz, g_x_0_xxxxz_xxxx, g_x_0_xxxxz_xxxxy, g_x_0_xxxxz_xxxy, g_x_0_xxxxz_xxxyy, g_x_0_xxxxz_xxxyz, g_x_0_xxxxz_xxxz, g_x_0_xxxxz_xxyy, g_x_0_xxxxz_xxyyy, g_x_0_xxxxz_xxyyz, g_x_0_xxxxz_xxyz, g_x_0_xxxxz_xxyzz, g_x_0_xxxxz_xxzz, g_x_0_xxxxz_xyyy, g_x_0_xxxxz_xyyyy, g_x_0_xxxxz_xyyyz, g_x_0_xxxxz_xyyz, g_x_0_xxxxz_xyyzz, g_x_0_xxxxz_xyzz, g_x_0_xxxxz_xyzzz, g_x_0_xxxxz_xzzz, g_x_0_xxxxz_yyyy, g_x_0_xxxxz_yyyyy, g_x_0_xxxxz_yyyyz, g_x_0_xxxxz_yyyz, g_x_0_xxxxz_yyyzz, g_x_0_xxxxz_yyzz, g_x_0_xxxxz_yyzzz, g_x_0_xxxxz_yzzz, g_x_0_xxxxz_yzzzz, g_x_0_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxyz_xxxx[k] = -g_x_0_xxxxz_xxxx[k] * ab_y + g_x_0_xxxxz_xxxxy[k];

                g_x_0_xxxxyz_xxxy[k] = -g_x_0_xxxxz_xxxy[k] * ab_y + g_x_0_xxxxz_xxxyy[k];

                g_x_0_xxxxyz_xxxz[k] = -g_x_0_xxxxz_xxxz[k] * ab_y + g_x_0_xxxxz_xxxyz[k];

                g_x_0_xxxxyz_xxyy[k] = -g_x_0_xxxxz_xxyy[k] * ab_y + g_x_0_xxxxz_xxyyy[k];

                g_x_0_xxxxyz_xxyz[k] = -g_x_0_xxxxz_xxyz[k] * ab_y + g_x_0_xxxxz_xxyyz[k];

                g_x_0_xxxxyz_xxzz[k] = -g_x_0_xxxxz_xxzz[k] * ab_y + g_x_0_xxxxz_xxyzz[k];

                g_x_0_xxxxyz_xyyy[k] = -g_x_0_xxxxz_xyyy[k] * ab_y + g_x_0_xxxxz_xyyyy[k];

                g_x_0_xxxxyz_xyyz[k] = -g_x_0_xxxxz_xyyz[k] * ab_y + g_x_0_xxxxz_xyyyz[k];

                g_x_0_xxxxyz_xyzz[k] = -g_x_0_xxxxz_xyzz[k] * ab_y + g_x_0_xxxxz_xyyzz[k];

                g_x_0_xxxxyz_xzzz[k] = -g_x_0_xxxxz_xzzz[k] * ab_y + g_x_0_xxxxz_xyzzz[k];

                g_x_0_xxxxyz_yyyy[k] = -g_x_0_xxxxz_yyyy[k] * ab_y + g_x_0_xxxxz_yyyyy[k];

                g_x_0_xxxxyz_yyyz[k] = -g_x_0_xxxxz_yyyz[k] * ab_y + g_x_0_xxxxz_yyyyz[k];

                g_x_0_xxxxyz_yyzz[k] = -g_x_0_xxxxz_yyzz[k] * ab_y + g_x_0_xxxxz_yyyzz[k];

                g_x_0_xxxxyz_yzzz[k] = -g_x_0_xxxxz_yzzz[k] * ab_y + g_x_0_xxxxz_yyzzz[k];

                g_x_0_xxxxyz_zzzz[k] = -g_x_0_xxxxz_zzzz[k] * ab_y + g_x_0_xxxxz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxxzz_xxxx = cbuffer.data(ig_geom_10_off + 75 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxy = cbuffer.data(ig_geom_10_off + 76 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxxz = cbuffer.data(ig_geom_10_off + 77 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxyy = cbuffer.data(ig_geom_10_off + 78 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxyz = cbuffer.data(ig_geom_10_off + 79 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xxzz = cbuffer.data(ig_geom_10_off + 80 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyyy = cbuffer.data(ig_geom_10_off + 81 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyyz = cbuffer.data(ig_geom_10_off + 82 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xyzz = cbuffer.data(ig_geom_10_off + 83 * ccomps * dcomps);

            auto g_x_0_xxxxzz_xzzz = cbuffer.data(ig_geom_10_off + 84 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyyy = cbuffer.data(ig_geom_10_off + 85 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyyz = cbuffer.data(ig_geom_10_off + 86 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yyzz = cbuffer.data(ig_geom_10_off + 87 * ccomps * dcomps);

            auto g_x_0_xxxxzz_yzzz = cbuffer.data(ig_geom_10_off + 88 * ccomps * dcomps);

            auto g_x_0_xxxxzz_zzzz = cbuffer.data(ig_geom_10_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxxz_xxxx, g_x_0_xxxxz_xxxxz, g_x_0_xxxxz_xxxy, g_x_0_xxxxz_xxxyz, g_x_0_xxxxz_xxxz, g_x_0_xxxxz_xxxzz, g_x_0_xxxxz_xxyy, g_x_0_xxxxz_xxyyz, g_x_0_xxxxz_xxyz, g_x_0_xxxxz_xxyzz, g_x_0_xxxxz_xxzz, g_x_0_xxxxz_xxzzz, g_x_0_xxxxz_xyyy, g_x_0_xxxxz_xyyyz, g_x_0_xxxxz_xyyz, g_x_0_xxxxz_xyyzz, g_x_0_xxxxz_xyzz, g_x_0_xxxxz_xyzzz, g_x_0_xxxxz_xzzz, g_x_0_xxxxz_xzzzz, g_x_0_xxxxz_yyyy, g_x_0_xxxxz_yyyyz, g_x_0_xxxxz_yyyz, g_x_0_xxxxz_yyyzz, g_x_0_xxxxz_yyzz, g_x_0_xxxxz_yyzzz, g_x_0_xxxxz_yzzz, g_x_0_xxxxz_yzzzz, g_x_0_xxxxz_zzzz, g_x_0_xxxxz_zzzzz, g_x_0_xxxxzz_xxxx, g_x_0_xxxxzz_xxxy, g_x_0_xxxxzz_xxxz, g_x_0_xxxxzz_xxyy, g_x_0_xxxxzz_xxyz, g_x_0_xxxxzz_xxzz, g_x_0_xxxxzz_xyyy, g_x_0_xxxxzz_xyyz, g_x_0_xxxxzz_xyzz, g_x_0_xxxxzz_xzzz, g_x_0_xxxxzz_yyyy, g_x_0_xxxxzz_yyyz, g_x_0_xxxxzz_yyzz, g_x_0_xxxxzz_yzzz, g_x_0_xxxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxxzz_xxxx[k] = -g_x_0_xxxxz_xxxx[k] * ab_z + g_x_0_xxxxz_xxxxz[k];

                g_x_0_xxxxzz_xxxy[k] = -g_x_0_xxxxz_xxxy[k] * ab_z + g_x_0_xxxxz_xxxyz[k];

                g_x_0_xxxxzz_xxxz[k] = -g_x_0_xxxxz_xxxz[k] * ab_z + g_x_0_xxxxz_xxxzz[k];

                g_x_0_xxxxzz_xxyy[k] = -g_x_0_xxxxz_xxyy[k] * ab_z + g_x_0_xxxxz_xxyyz[k];

                g_x_0_xxxxzz_xxyz[k] = -g_x_0_xxxxz_xxyz[k] * ab_z + g_x_0_xxxxz_xxyzz[k];

                g_x_0_xxxxzz_xxzz[k] = -g_x_0_xxxxz_xxzz[k] * ab_z + g_x_0_xxxxz_xxzzz[k];

                g_x_0_xxxxzz_xyyy[k] = -g_x_0_xxxxz_xyyy[k] * ab_z + g_x_0_xxxxz_xyyyz[k];

                g_x_0_xxxxzz_xyyz[k] = -g_x_0_xxxxz_xyyz[k] * ab_z + g_x_0_xxxxz_xyyzz[k];

                g_x_0_xxxxzz_xyzz[k] = -g_x_0_xxxxz_xyzz[k] * ab_z + g_x_0_xxxxz_xyzzz[k];

                g_x_0_xxxxzz_xzzz[k] = -g_x_0_xxxxz_xzzz[k] * ab_z + g_x_0_xxxxz_xzzzz[k];

                g_x_0_xxxxzz_yyyy[k] = -g_x_0_xxxxz_yyyy[k] * ab_z + g_x_0_xxxxz_yyyyz[k];

                g_x_0_xxxxzz_yyyz[k] = -g_x_0_xxxxz_yyyz[k] * ab_z + g_x_0_xxxxz_yyyzz[k];

                g_x_0_xxxxzz_yyzz[k] = -g_x_0_xxxxz_yyzz[k] * ab_z + g_x_0_xxxxz_yyzzz[k];

                g_x_0_xxxxzz_yzzz[k] = -g_x_0_xxxxz_yzzz[k] * ab_z + g_x_0_xxxxz_yzzzz[k];

                g_x_0_xxxxzz_zzzz[k] = -g_x_0_xxxxz_zzzz[k] * ab_z + g_x_0_xxxxz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyy_xxxx = cbuffer.data(ig_geom_10_off + 90 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxy = cbuffer.data(ig_geom_10_off + 91 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxxz = cbuffer.data(ig_geom_10_off + 92 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxyy = cbuffer.data(ig_geom_10_off + 93 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxyz = cbuffer.data(ig_geom_10_off + 94 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xxzz = cbuffer.data(ig_geom_10_off + 95 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyyy = cbuffer.data(ig_geom_10_off + 96 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyyz = cbuffer.data(ig_geom_10_off + 97 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xyzz = cbuffer.data(ig_geom_10_off + 98 * ccomps * dcomps);

            auto g_x_0_xxxyyy_xzzz = cbuffer.data(ig_geom_10_off + 99 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyyy = cbuffer.data(ig_geom_10_off + 100 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyyz = cbuffer.data(ig_geom_10_off + 101 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yyzz = cbuffer.data(ig_geom_10_off + 102 * ccomps * dcomps);

            auto g_x_0_xxxyyy_yzzz = cbuffer.data(ig_geom_10_off + 103 * ccomps * dcomps);

            auto g_x_0_xxxyyy_zzzz = cbuffer.data(ig_geom_10_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyy_xxxx, g_x_0_xxxyy_xxxxy, g_x_0_xxxyy_xxxy, g_x_0_xxxyy_xxxyy, g_x_0_xxxyy_xxxyz, g_x_0_xxxyy_xxxz, g_x_0_xxxyy_xxyy, g_x_0_xxxyy_xxyyy, g_x_0_xxxyy_xxyyz, g_x_0_xxxyy_xxyz, g_x_0_xxxyy_xxyzz, g_x_0_xxxyy_xxzz, g_x_0_xxxyy_xyyy, g_x_0_xxxyy_xyyyy, g_x_0_xxxyy_xyyyz, g_x_0_xxxyy_xyyz, g_x_0_xxxyy_xyyzz, g_x_0_xxxyy_xyzz, g_x_0_xxxyy_xyzzz, g_x_0_xxxyy_xzzz, g_x_0_xxxyy_yyyy, g_x_0_xxxyy_yyyyy, g_x_0_xxxyy_yyyyz, g_x_0_xxxyy_yyyz, g_x_0_xxxyy_yyyzz, g_x_0_xxxyy_yyzz, g_x_0_xxxyy_yyzzz, g_x_0_xxxyy_yzzz, g_x_0_xxxyy_yzzzz, g_x_0_xxxyy_zzzz, g_x_0_xxxyyy_xxxx, g_x_0_xxxyyy_xxxy, g_x_0_xxxyyy_xxxz, g_x_0_xxxyyy_xxyy, g_x_0_xxxyyy_xxyz, g_x_0_xxxyyy_xxzz, g_x_0_xxxyyy_xyyy, g_x_0_xxxyyy_xyyz, g_x_0_xxxyyy_xyzz, g_x_0_xxxyyy_xzzz, g_x_0_xxxyyy_yyyy, g_x_0_xxxyyy_yyyz, g_x_0_xxxyyy_yyzz, g_x_0_xxxyyy_yzzz, g_x_0_xxxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyy_xxxx[k] = -g_x_0_xxxyy_xxxx[k] * ab_y + g_x_0_xxxyy_xxxxy[k];

                g_x_0_xxxyyy_xxxy[k] = -g_x_0_xxxyy_xxxy[k] * ab_y + g_x_0_xxxyy_xxxyy[k];

                g_x_0_xxxyyy_xxxz[k] = -g_x_0_xxxyy_xxxz[k] * ab_y + g_x_0_xxxyy_xxxyz[k];

                g_x_0_xxxyyy_xxyy[k] = -g_x_0_xxxyy_xxyy[k] * ab_y + g_x_0_xxxyy_xxyyy[k];

                g_x_0_xxxyyy_xxyz[k] = -g_x_0_xxxyy_xxyz[k] * ab_y + g_x_0_xxxyy_xxyyz[k];

                g_x_0_xxxyyy_xxzz[k] = -g_x_0_xxxyy_xxzz[k] * ab_y + g_x_0_xxxyy_xxyzz[k];

                g_x_0_xxxyyy_xyyy[k] = -g_x_0_xxxyy_xyyy[k] * ab_y + g_x_0_xxxyy_xyyyy[k];

                g_x_0_xxxyyy_xyyz[k] = -g_x_0_xxxyy_xyyz[k] * ab_y + g_x_0_xxxyy_xyyyz[k];

                g_x_0_xxxyyy_xyzz[k] = -g_x_0_xxxyy_xyzz[k] * ab_y + g_x_0_xxxyy_xyyzz[k];

                g_x_0_xxxyyy_xzzz[k] = -g_x_0_xxxyy_xzzz[k] * ab_y + g_x_0_xxxyy_xyzzz[k];

                g_x_0_xxxyyy_yyyy[k] = -g_x_0_xxxyy_yyyy[k] * ab_y + g_x_0_xxxyy_yyyyy[k];

                g_x_0_xxxyyy_yyyz[k] = -g_x_0_xxxyy_yyyz[k] * ab_y + g_x_0_xxxyy_yyyyz[k];

                g_x_0_xxxyyy_yyzz[k] = -g_x_0_xxxyy_yyzz[k] * ab_y + g_x_0_xxxyy_yyyzz[k];

                g_x_0_xxxyyy_yzzz[k] = -g_x_0_xxxyy_yzzz[k] * ab_y + g_x_0_xxxyy_yyzzz[k];

                g_x_0_xxxyyy_zzzz[k] = -g_x_0_xxxyy_zzzz[k] * ab_y + g_x_0_xxxyy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyyz_xxxx = cbuffer.data(ig_geom_10_off + 105 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxy = cbuffer.data(ig_geom_10_off + 106 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxxz = cbuffer.data(ig_geom_10_off + 107 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxyy = cbuffer.data(ig_geom_10_off + 108 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxyz = cbuffer.data(ig_geom_10_off + 109 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xxzz = cbuffer.data(ig_geom_10_off + 110 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyyy = cbuffer.data(ig_geom_10_off + 111 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyyz = cbuffer.data(ig_geom_10_off + 112 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xyzz = cbuffer.data(ig_geom_10_off + 113 * ccomps * dcomps);

            auto g_x_0_xxxyyz_xzzz = cbuffer.data(ig_geom_10_off + 114 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyyy = cbuffer.data(ig_geom_10_off + 115 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyyz = cbuffer.data(ig_geom_10_off + 116 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yyzz = cbuffer.data(ig_geom_10_off + 117 * ccomps * dcomps);

            auto g_x_0_xxxyyz_yzzz = cbuffer.data(ig_geom_10_off + 118 * ccomps * dcomps);

            auto g_x_0_xxxyyz_zzzz = cbuffer.data(ig_geom_10_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyyz_xxxx, g_x_0_xxxyyz_xxxy, g_x_0_xxxyyz_xxxz, g_x_0_xxxyyz_xxyy, g_x_0_xxxyyz_xxyz, g_x_0_xxxyyz_xxzz, g_x_0_xxxyyz_xyyy, g_x_0_xxxyyz_xyyz, g_x_0_xxxyyz_xyzz, g_x_0_xxxyyz_xzzz, g_x_0_xxxyyz_yyyy, g_x_0_xxxyyz_yyyz, g_x_0_xxxyyz_yyzz, g_x_0_xxxyyz_yzzz, g_x_0_xxxyyz_zzzz, g_x_0_xxxyz_xxxx, g_x_0_xxxyz_xxxxy, g_x_0_xxxyz_xxxy, g_x_0_xxxyz_xxxyy, g_x_0_xxxyz_xxxyz, g_x_0_xxxyz_xxxz, g_x_0_xxxyz_xxyy, g_x_0_xxxyz_xxyyy, g_x_0_xxxyz_xxyyz, g_x_0_xxxyz_xxyz, g_x_0_xxxyz_xxyzz, g_x_0_xxxyz_xxzz, g_x_0_xxxyz_xyyy, g_x_0_xxxyz_xyyyy, g_x_0_xxxyz_xyyyz, g_x_0_xxxyz_xyyz, g_x_0_xxxyz_xyyzz, g_x_0_xxxyz_xyzz, g_x_0_xxxyz_xyzzz, g_x_0_xxxyz_xzzz, g_x_0_xxxyz_yyyy, g_x_0_xxxyz_yyyyy, g_x_0_xxxyz_yyyyz, g_x_0_xxxyz_yyyz, g_x_0_xxxyz_yyyzz, g_x_0_xxxyz_yyzz, g_x_0_xxxyz_yyzzz, g_x_0_xxxyz_yzzz, g_x_0_xxxyz_yzzzz, g_x_0_xxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyyz_xxxx[k] = -g_x_0_xxxyz_xxxx[k] * ab_y + g_x_0_xxxyz_xxxxy[k];

                g_x_0_xxxyyz_xxxy[k] = -g_x_0_xxxyz_xxxy[k] * ab_y + g_x_0_xxxyz_xxxyy[k];

                g_x_0_xxxyyz_xxxz[k] = -g_x_0_xxxyz_xxxz[k] * ab_y + g_x_0_xxxyz_xxxyz[k];

                g_x_0_xxxyyz_xxyy[k] = -g_x_0_xxxyz_xxyy[k] * ab_y + g_x_0_xxxyz_xxyyy[k];

                g_x_0_xxxyyz_xxyz[k] = -g_x_0_xxxyz_xxyz[k] * ab_y + g_x_0_xxxyz_xxyyz[k];

                g_x_0_xxxyyz_xxzz[k] = -g_x_0_xxxyz_xxzz[k] * ab_y + g_x_0_xxxyz_xxyzz[k];

                g_x_0_xxxyyz_xyyy[k] = -g_x_0_xxxyz_xyyy[k] * ab_y + g_x_0_xxxyz_xyyyy[k];

                g_x_0_xxxyyz_xyyz[k] = -g_x_0_xxxyz_xyyz[k] * ab_y + g_x_0_xxxyz_xyyyz[k];

                g_x_0_xxxyyz_xyzz[k] = -g_x_0_xxxyz_xyzz[k] * ab_y + g_x_0_xxxyz_xyyzz[k];

                g_x_0_xxxyyz_xzzz[k] = -g_x_0_xxxyz_xzzz[k] * ab_y + g_x_0_xxxyz_xyzzz[k];

                g_x_0_xxxyyz_yyyy[k] = -g_x_0_xxxyz_yyyy[k] * ab_y + g_x_0_xxxyz_yyyyy[k];

                g_x_0_xxxyyz_yyyz[k] = -g_x_0_xxxyz_yyyz[k] * ab_y + g_x_0_xxxyz_yyyyz[k];

                g_x_0_xxxyyz_yyzz[k] = -g_x_0_xxxyz_yyzz[k] * ab_y + g_x_0_xxxyz_yyyzz[k];

                g_x_0_xxxyyz_yzzz[k] = -g_x_0_xxxyz_yzzz[k] * ab_y + g_x_0_xxxyz_yyzzz[k];

                g_x_0_xxxyyz_zzzz[k] = -g_x_0_xxxyz_zzzz[k] * ab_y + g_x_0_xxxyz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxyzz_xxxx = cbuffer.data(ig_geom_10_off + 120 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxy = cbuffer.data(ig_geom_10_off + 121 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxxz = cbuffer.data(ig_geom_10_off + 122 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxyy = cbuffer.data(ig_geom_10_off + 123 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxyz = cbuffer.data(ig_geom_10_off + 124 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xxzz = cbuffer.data(ig_geom_10_off + 125 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyyy = cbuffer.data(ig_geom_10_off + 126 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyyz = cbuffer.data(ig_geom_10_off + 127 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xyzz = cbuffer.data(ig_geom_10_off + 128 * ccomps * dcomps);

            auto g_x_0_xxxyzz_xzzz = cbuffer.data(ig_geom_10_off + 129 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyyy = cbuffer.data(ig_geom_10_off + 130 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyyz = cbuffer.data(ig_geom_10_off + 131 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yyzz = cbuffer.data(ig_geom_10_off + 132 * ccomps * dcomps);

            auto g_x_0_xxxyzz_yzzz = cbuffer.data(ig_geom_10_off + 133 * ccomps * dcomps);

            auto g_x_0_xxxyzz_zzzz = cbuffer.data(ig_geom_10_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxyzz_xxxx, g_x_0_xxxyzz_xxxy, g_x_0_xxxyzz_xxxz, g_x_0_xxxyzz_xxyy, g_x_0_xxxyzz_xxyz, g_x_0_xxxyzz_xxzz, g_x_0_xxxyzz_xyyy, g_x_0_xxxyzz_xyyz, g_x_0_xxxyzz_xyzz, g_x_0_xxxyzz_xzzz, g_x_0_xxxyzz_yyyy, g_x_0_xxxyzz_yyyz, g_x_0_xxxyzz_yyzz, g_x_0_xxxyzz_yzzz, g_x_0_xxxyzz_zzzz, g_x_0_xxxzz_xxxx, g_x_0_xxxzz_xxxxy, g_x_0_xxxzz_xxxy, g_x_0_xxxzz_xxxyy, g_x_0_xxxzz_xxxyz, g_x_0_xxxzz_xxxz, g_x_0_xxxzz_xxyy, g_x_0_xxxzz_xxyyy, g_x_0_xxxzz_xxyyz, g_x_0_xxxzz_xxyz, g_x_0_xxxzz_xxyzz, g_x_0_xxxzz_xxzz, g_x_0_xxxzz_xyyy, g_x_0_xxxzz_xyyyy, g_x_0_xxxzz_xyyyz, g_x_0_xxxzz_xyyz, g_x_0_xxxzz_xyyzz, g_x_0_xxxzz_xyzz, g_x_0_xxxzz_xyzzz, g_x_0_xxxzz_xzzz, g_x_0_xxxzz_yyyy, g_x_0_xxxzz_yyyyy, g_x_0_xxxzz_yyyyz, g_x_0_xxxzz_yyyz, g_x_0_xxxzz_yyyzz, g_x_0_xxxzz_yyzz, g_x_0_xxxzz_yyzzz, g_x_0_xxxzz_yzzz, g_x_0_xxxzz_yzzzz, g_x_0_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxyzz_xxxx[k] = -g_x_0_xxxzz_xxxx[k] * ab_y + g_x_0_xxxzz_xxxxy[k];

                g_x_0_xxxyzz_xxxy[k] = -g_x_0_xxxzz_xxxy[k] * ab_y + g_x_0_xxxzz_xxxyy[k];

                g_x_0_xxxyzz_xxxz[k] = -g_x_0_xxxzz_xxxz[k] * ab_y + g_x_0_xxxzz_xxxyz[k];

                g_x_0_xxxyzz_xxyy[k] = -g_x_0_xxxzz_xxyy[k] * ab_y + g_x_0_xxxzz_xxyyy[k];

                g_x_0_xxxyzz_xxyz[k] = -g_x_0_xxxzz_xxyz[k] * ab_y + g_x_0_xxxzz_xxyyz[k];

                g_x_0_xxxyzz_xxzz[k] = -g_x_0_xxxzz_xxzz[k] * ab_y + g_x_0_xxxzz_xxyzz[k];

                g_x_0_xxxyzz_xyyy[k] = -g_x_0_xxxzz_xyyy[k] * ab_y + g_x_0_xxxzz_xyyyy[k];

                g_x_0_xxxyzz_xyyz[k] = -g_x_0_xxxzz_xyyz[k] * ab_y + g_x_0_xxxzz_xyyyz[k];

                g_x_0_xxxyzz_xyzz[k] = -g_x_0_xxxzz_xyzz[k] * ab_y + g_x_0_xxxzz_xyyzz[k];

                g_x_0_xxxyzz_xzzz[k] = -g_x_0_xxxzz_xzzz[k] * ab_y + g_x_0_xxxzz_xyzzz[k];

                g_x_0_xxxyzz_yyyy[k] = -g_x_0_xxxzz_yyyy[k] * ab_y + g_x_0_xxxzz_yyyyy[k];

                g_x_0_xxxyzz_yyyz[k] = -g_x_0_xxxzz_yyyz[k] * ab_y + g_x_0_xxxzz_yyyyz[k];

                g_x_0_xxxyzz_yyzz[k] = -g_x_0_xxxzz_yyzz[k] * ab_y + g_x_0_xxxzz_yyyzz[k];

                g_x_0_xxxyzz_yzzz[k] = -g_x_0_xxxzz_yzzz[k] * ab_y + g_x_0_xxxzz_yyzzz[k];

                g_x_0_xxxyzz_zzzz[k] = -g_x_0_xxxzz_zzzz[k] * ab_y + g_x_0_xxxzz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxxzzz_xxxx = cbuffer.data(ig_geom_10_off + 135 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxy = cbuffer.data(ig_geom_10_off + 136 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxxz = cbuffer.data(ig_geom_10_off + 137 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxyy = cbuffer.data(ig_geom_10_off + 138 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxyz = cbuffer.data(ig_geom_10_off + 139 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xxzz = cbuffer.data(ig_geom_10_off + 140 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyyy = cbuffer.data(ig_geom_10_off + 141 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyyz = cbuffer.data(ig_geom_10_off + 142 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xyzz = cbuffer.data(ig_geom_10_off + 143 * ccomps * dcomps);

            auto g_x_0_xxxzzz_xzzz = cbuffer.data(ig_geom_10_off + 144 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyyy = cbuffer.data(ig_geom_10_off + 145 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyyz = cbuffer.data(ig_geom_10_off + 146 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yyzz = cbuffer.data(ig_geom_10_off + 147 * ccomps * dcomps);

            auto g_x_0_xxxzzz_yzzz = cbuffer.data(ig_geom_10_off + 148 * ccomps * dcomps);

            auto g_x_0_xxxzzz_zzzz = cbuffer.data(ig_geom_10_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxxzz_xxxx, g_x_0_xxxzz_xxxxz, g_x_0_xxxzz_xxxy, g_x_0_xxxzz_xxxyz, g_x_0_xxxzz_xxxz, g_x_0_xxxzz_xxxzz, g_x_0_xxxzz_xxyy, g_x_0_xxxzz_xxyyz, g_x_0_xxxzz_xxyz, g_x_0_xxxzz_xxyzz, g_x_0_xxxzz_xxzz, g_x_0_xxxzz_xxzzz, g_x_0_xxxzz_xyyy, g_x_0_xxxzz_xyyyz, g_x_0_xxxzz_xyyz, g_x_0_xxxzz_xyyzz, g_x_0_xxxzz_xyzz, g_x_0_xxxzz_xyzzz, g_x_0_xxxzz_xzzz, g_x_0_xxxzz_xzzzz, g_x_0_xxxzz_yyyy, g_x_0_xxxzz_yyyyz, g_x_0_xxxzz_yyyz, g_x_0_xxxzz_yyyzz, g_x_0_xxxzz_yyzz, g_x_0_xxxzz_yyzzz, g_x_0_xxxzz_yzzz, g_x_0_xxxzz_yzzzz, g_x_0_xxxzz_zzzz, g_x_0_xxxzz_zzzzz, g_x_0_xxxzzz_xxxx, g_x_0_xxxzzz_xxxy, g_x_0_xxxzzz_xxxz, g_x_0_xxxzzz_xxyy, g_x_0_xxxzzz_xxyz, g_x_0_xxxzzz_xxzz, g_x_0_xxxzzz_xyyy, g_x_0_xxxzzz_xyyz, g_x_0_xxxzzz_xyzz, g_x_0_xxxzzz_xzzz, g_x_0_xxxzzz_yyyy, g_x_0_xxxzzz_yyyz, g_x_0_xxxzzz_yyzz, g_x_0_xxxzzz_yzzz, g_x_0_xxxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxxzzz_xxxx[k] = -g_x_0_xxxzz_xxxx[k] * ab_z + g_x_0_xxxzz_xxxxz[k];

                g_x_0_xxxzzz_xxxy[k] = -g_x_0_xxxzz_xxxy[k] * ab_z + g_x_0_xxxzz_xxxyz[k];

                g_x_0_xxxzzz_xxxz[k] = -g_x_0_xxxzz_xxxz[k] * ab_z + g_x_0_xxxzz_xxxzz[k];

                g_x_0_xxxzzz_xxyy[k] = -g_x_0_xxxzz_xxyy[k] * ab_z + g_x_0_xxxzz_xxyyz[k];

                g_x_0_xxxzzz_xxyz[k] = -g_x_0_xxxzz_xxyz[k] * ab_z + g_x_0_xxxzz_xxyzz[k];

                g_x_0_xxxzzz_xxzz[k] = -g_x_0_xxxzz_xxzz[k] * ab_z + g_x_0_xxxzz_xxzzz[k];

                g_x_0_xxxzzz_xyyy[k] = -g_x_0_xxxzz_xyyy[k] * ab_z + g_x_0_xxxzz_xyyyz[k];

                g_x_0_xxxzzz_xyyz[k] = -g_x_0_xxxzz_xyyz[k] * ab_z + g_x_0_xxxzz_xyyzz[k];

                g_x_0_xxxzzz_xyzz[k] = -g_x_0_xxxzz_xyzz[k] * ab_z + g_x_0_xxxzz_xyzzz[k];

                g_x_0_xxxzzz_xzzz[k] = -g_x_0_xxxzz_xzzz[k] * ab_z + g_x_0_xxxzz_xzzzz[k];

                g_x_0_xxxzzz_yyyy[k] = -g_x_0_xxxzz_yyyy[k] * ab_z + g_x_0_xxxzz_yyyyz[k];

                g_x_0_xxxzzz_yyyz[k] = -g_x_0_xxxzz_yyyz[k] * ab_z + g_x_0_xxxzz_yyyzz[k];

                g_x_0_xxxzzz_yyzz[k] = -g_x_0_xxxzz_yyzz[k] * ab_z + g_x_0_xxxzz_yyzzz[k];

                g_x_0_xxxzzz_yzzz[k] = -g_x_0_xxxzz_yzzz[k] * ab_z + g_x_0_xxxzz_yzzzz[k];

                g_x_0_xxxzzz_zzzz[k] = -g_x_0_xxxzz_zzzz[k] * ab_z + g_x_0_xxxzz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyy_xxxx = cbuffer.data(ig_geom_10_off + 150 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxy = cbuffer.data(ig_geom_10_off + 151 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxxz = cbuffer.data(ig_geom_10_off + 152 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxyy = cbuffer.data(ig_geom_10_off + 153 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxyz = cbuffer.data(ig_geom_10_off + 154 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xxzz = cbuffer.data(ig_geom_10_off + 155 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyyy = cbuffer.data(ig_geom_10_off + 156 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyyz = cbuffer.data(ig_geom_10_off + 157 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xyzz = cbuffer.data(ig_geom_10_off + 158 * ccomps * dcomps);

            auto g_x_0_xxyyyy_xzzz = cbuffer.data(ig_geom_10_off + 159 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyyy = cbuffer.data(ig_geom_10_off + 160 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyyz = cbuffer.data(ig_geom_10_off + 161 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yyzz = cbuffer.data(ig_geom_10_off + 162 * ccomps * dcomps);

            auto g_x_0_xxyyyy_yzzz = cbuffer.data(ig_geom_10_off + 163 * ccomps * dcomps);

            auto g_x_0_xxyyyy_zzzz = cbuffer.data(ig_geom_10_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyy_xxxx, g_x_0_xxyyy_xxxxy, g_x_0_xxyyy_xxxy, g_x_0_xxyyy_xxxyy, g_x_0_xxyyy_xxxyz, g_x_0_xxyyy_xxxz, g_x_0_xxyyy_xxyy, g_x_0_xxyyy_xxyyy, g_x_0_xxyyy_xxyyz, g_x_0_xxyyy_xxyz, g_x_0_xxyyy_xxyzz, g_x_0_xxyyy_xxzz, g_x_0_xxyyy_xyyy, g_x_0_xxyyy_xyyyy, g_x_0_xxyyy_xyyyz, g_x_0_xxyyy_xyyz, g_x_0_xxyyy_xyyzz, g_x_0_xxyyy_xyzz, g_x_0_xxyyy_xyzzz, g_x_0_xxyyy_xzzz, g_x_0_xxyyy_yyyy, g_x_0_xxyyy_yyyyy, g_x_0_xxyyy_yyyyz, g_x_0_xxyyy_yyyz, g_x_0_xxyyy_yyyzz, g_x_0_xxyyy_yyzz, g_x_0_xxyyy_yyzzz, g_x_0_xxyyy_yzzz, g_x_0_xxyyy_yzzzz, g_x_0_xxyyy_zzzz, g_x_0_xxyyyy_xxxx, g_x_0_xxyyyy_xxxy, g_x_0_xxyyyy_xxxz, g_x_0_xxyyyy_xxyy, g_x_0_xxyyyy_xxyz, g_x_0_xxyyyy_xxzz, g_x_0_xxyyyy_xyyy, g_x_0_xxyyyy_xyyz, g_x_0_xxyyyy_xyzz, g_x_0_xxyyyy_xzzz, g_x_0_xxyyyy_yyyy, g_x_0_xxyyyy_yyyz, g_x_0_xxyyyy_yyzz, g_x_0_xxyyyy_yzzz, g_x_0_xxyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyy_xxxx[k] = -g_x_0_xxyyy_xxxx[k] * ab_y + g_x_0_xxyyy_xxxxy[k];

                g_x_0_xxyyyy_xxxy[k] = -g_x_0_xxyyy_xxxy[k] * ab_y + g_x_0_xxyyy_xxxyy[k];

                g_x_0_xxyyyy_xxxz[k] = -g_x_0_xxyyy_xxxz[k] * ab_y + g_x_0_xxyyy_xxxyz[k];

                g_x_0_xxyyyy_xxyy[k] = -g_x_0_xxyyy_xxyy[k] * ab_y + g_x_0_xxyyy_xxyyy[k];

                g_x_0_xxyyyy_xxyz[k] = -g_x_0_xxyyy_xxyz[k] * ab_y + g_x_0_xxyyy_xxyyz[k];

                g_x_0_xxyyyy_xxzz[k] = -g_x_0_xxyyy_xxzz[k] * ab_y + g_x_0_xxyyy_xxyzz[k];

                g_x_0_xxyyyy_xyyy[k] = -g_x_0_xxyyy_xyyy[k] * ab_y + g_x_0_xxyyy_xyyyy[k];

                g_x_0_xxyyyy_xyyz[k] = -g_x_0_xxyyy_xyyz[k] * ab_y + g_x_0_xxyyy_xyyyz[k];

                g_x_0_xxyyyy_xyzz[k] = -g_x_0_xxyyy_xyzz[k] * ab_y + g_x_0_xxyyy_xyyzz[k];

                g_x_0_xxyyyy_xzzz[k] = -g_x_0_xxyyy_xzzz[k] * ab_y + g_x_0_xxyyy_xyzzz[k];

                g_x_0_xxyyyy_yyyy[k] = -g_x_0_xxyyy_yyyy[k] * ab_y + g_x_0_xxyyy_yyyyy[k];

                g_x_0_xxyyyy_yyyz[k] = -g_x_0_xxyyy_yyyz[k] * ab_y + g_x_0_xxyyy_yyyyz[k];

                g_x_0_xxyyyy_yyzz[k] = -g_x_0_xxyyy_yyzz[k] * ab_y + g_x_0_xxyyy_yyyzz[k];

                g_x_0_xxyyyy_yzzz[k] = -g_x_0_xxyyy_yzzz[k] * ab_y + g_x_0_xxyyy_yyzzz[k];

                g_x_0_xxyyyy_zzzz[k] = -g_x_0_xxyyy_zzzz[k] * ab_y + g_x_0_xxyyy_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyyz_xxxx = cbuffer.data(ig_geom_10_off + 165 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxy = cbuffer.data(ig_geom_10_off + 166 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxxz = cbuffer.data(ig_geom_10_off + 167 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxyy = cbuffer.data(ig_geom_10_off + 168 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxyz = cbuffer.data(ig_geom_10_off + 169 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xxzz = cbuffer.data(ig_geom_10_off + 170 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyyy = cbuffer.data(ig_geom_10_off + 171 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyyz = cbuffer.data(ig_geom_10_off + 172 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xyzz = cbuffer.data(ig_geom_10_off + 173 * ccomps * dcomps);

            auto g_x_0_xxyyyz_xzzz = cbuffer.data(ig_geom_10_off + 174 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyyy = cbuffer.data(ig_geom_10_off + 175 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyyz = cbuffer.data(ig_geom_10_off + 176 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yyzz = cbuffer.data(ig_geom_10_off + 177 * ccomps * dcomps);

            auto g_x_0_xxyyyz_yzzz = cbuffer.data(ig_geom_10_off + 178 * ccomps * dcomps);

            auto g_x_0_xxyyyz_zzzz = cbuffer.data(ig_geom_10_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyyz_xxxx, g_x_0_xxyyyz_xxxy, g_x_0_xxyyyz_xxxz, g_x_0_xxyyyz_xxyy, g_x_0_xxyyyz_xxyz, g_x_0_xxyyyz_xxzz, g_x_0_xxyyyz_xyyy, g_x_0_xxyyyz_xyyz, g_x_0_xxyyyz_xyzz, g_x_0_xxyyyz_xzzz, g_x_0_xxyyyz_yyyy, g_x_0_xxyyyz_yyyz, g_x_0_xxyyyz_yyzz, g_x_0_xxyyyz_yzzz, g_x_0_xxyyyz_zzzz, g_x_0_xxyyz_xxxx, g_x_0_xxyyz_xxxxy, g_x_0_xxyyz_xxxy, g_x_0_xxyyz_xxxyy, g_x_0_xxyyz_xxxyz, g_x_0_xxyyz_xxxz, g_x_0_xxyyz_xxyy, g_x_0_xxyyz_xxyyy, g_x_0_xxyyz_xxyyz, g_x_0_xxyyz_xxyz, g_x_0_xxyyz_xxyzz, g_x_0_xxyyz_xxzz, g_x_0_xxyyz_xyyy, g_x_0_xxyyz_xyyyy, g_x_0_xxyyz_xyyyz, g_x_0_xxyyz_xyyz, g_x_0_xxyyz_xyyzz, g_x_0_xxyyz_xyzz, g_x_0_xxyyz_xyzzz, g_x_0_xxyyz_xzzz, g_x_0_xxyyz_yyyy, g_x_0_xxyyz_yyyyy, g_x_0_xxyyz_yyyyz, g_x_0_xxyyz_yyyz, g_x_0_xxyyz_yyyzz, g_x_0_xxyyz_yyzz, g_x_0_xxyyz_yyzzz, g_x_0_xxyyz_yzzz, g_x_0_xxyyz_yzzzz, g_x_0_xxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyyz_xxxx[k] = -g_x_0_xxyyz_xxxx[k] * ab_y + g_x_0_xxyyz_xxxxy[k];

                g_x_0_xxyyyz_xxxy[k] = -g_x_0_xxyyz_xxxy[k] * ab_y + g_x_0_xxyyz_xxxyy[k];

                g_x_0_xxyyyz_xxxz[k] = -g_x_0_xxyyz_xxxz[k] * ab_y + g_x_0_xxyyz_xxxyz[k];

                g_x_0_xxyyyz_xxyy[k] = -g_x_0_xxyyz_xxyy[k] * ab_y + g_x_0_xxyyz_xxyyy[k];

                g_x_0_xxyyyz_xxyz[k] = -g_x_0_xxyyz_xxyz[k] * ab_y + g_x_0_xxyyz_xxyyz[k];

                g_x_0_xxyyyz_xxzz[k] = -g_x_0_xxyyz_xxzz[k] * ab_y + g_x_0_xxyyz_xxyzz[k];

                g_x_0_xxyyyz_xyyy[k] = -g_x_0_xxyyz_xyyy[k] * ab_y + g_x_0_xxyyz_xyyyy[k];

                g_x_0_xxyyyz_xyyz[k] = -g_x_0_xxyyz_xyyz[k] * ab_y + g_x_0_xxyyz_xyyyz[k];

                g_x_0_xxyyyz_xyzz[k] = -g_x_0_xxyyz_xyzz[k] * ab_y + g_x_0_xxyyz_xyyzz[k];

                g_x_0_xxyyyz_xzzz[k] = -g_x_0_xxyyz_xzzz[k] * ab_y + g_x_0_xxyyz_xyzzz[k];

                g_x_0_xxyyyz_yyyy[k] = -g_x_0_xxyyz_yyyy[k] * ab_y + g_x_0_xxyyz_yyyyy[k];

                g_x_0_xxyyyz_yyyz[k] = -g_x_0_xxyyz_yyyz[k] * ab_y + g_x_0_xxyyz_yyyyz[k];

                g_x_0_xxyyyz_yyzz[k] = -g_x_0_xxyyz_yyzz[k] * ab_y + g_x_0_xxyyz_yyyzz[k];

                g_x_0_xxyyyz_yzzz[k] = -g_x_0_xxyyz_yzzz[k] * ab_y + g_x_0_xxyyz_yyzzz[k];

                g_x_0_xxyyyz_zzzz[k] = -g_x_0_xxyyz_zzzz[k] * ab_y + g_x_0_xxyyz_yzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyyzz_xxxx = cbuffer.data(ig_geom_10_off + 180 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxy = cbuffer.data(ig_geom_10_off + 181 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxxz = cbuffer.data(ig_geom_10_off + 182 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxyy = cbuffer.data(ig_geom_10_off + 183 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxyz = cbuffer.data(ig_geom_10_off + 184 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xxzz = cbuffer.data(ig_geom_10_off + 185 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyyy = cbuffer.data(ig_geom_10_off + 186 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyyz = cbuffer.data(ig_geom_10_off + 187 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xyzz = cbuffer.data(ig_geom_10_off + 188 * ccomps * dcomps);

            auto g_x_0_xxyyzz_xzzz = cbuffer.data(ig_geom_10_off + 189 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyyy = cbuffer.data(ig_geom_10_off + 190 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyyz = cbuffer.data(ig_geom_10_off + 191 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yyzz = cbuffer.data(ig_geom_10_off + 192 * ccomps * dcomps);

            auto g_x_0_xxyyzz_yzzz = cbuffer.data(ig_geom_10_off + 193 * ccomps * dcomps);

            auto g_x_0_xxyyzz_zzzz = cbuffer.data(ig_geom_10_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyyzz_xxxx, g_x_0_xxyyzz_xxxy, g_x_0_xxyyzz_xxxz, g_x_0_xxyyzz_xxyy, g_x_0_xxyyzz_xxyz, g_x_0_xxyyzz_xxzz, g_x_0_xxyyzz_xyyy, g_x_0_xxyyzz_xyyz, g_x_0_xxyyzz_xyzz, g_x_0_xxyyzz_xzzz, g_x_0_xxyyzz_yyyy, g_x_0_xxyyzz_yyyz, g_x_0_xxyyzz_yyzz, g_x_0_xxyyzz_yzzz, g_x_0_xxyyzz_zzzz, g_x_0_xxyzz_xxxx, g_x_0_xxyzz_xxxxy, g_x_0_xxyzz_xxxy, g_x_0_xxyzz_xxxyy, g_x_0_xxyzz_xxxyz, g_x_0_xxyzz_xxxz, g_x_0_xxyzz_xxyy, g_x_0_xxyzz_xxyyy, g_x_0_xxyzz_xxyyz, g_x_0_xxyzz_xxyz, g_x_0_xxyzz_xxyzz, g_x_0_xxyzz_xxzz, g_x_0_xxyzz_xyyy, g_x_0_xxyzz_xyyyy, g_x_0_xxyzz_xyyyz, g_x_0_xxyzz_xyyz, g_x_0_xxyzz_xyyzz, g_x_0_xxyzz_xyzz, g_x_0_xxyzz_xyzzz, g_x_0_xxyzz_xzzz, g_x_0_xxyzz_yyyy, g_x_0_xxyzz_yyyyy, g_x_0_xxyzz_yyyyz, g_x_0_xxyzz_yyyz, g_x_0_xxyzz_yyyzz, g_x_0_xxyzz_yyzz, g_x_0_xxyzz_yyzzz, g_x_0_xxyzz_yzzz, g_x_0_xxyzz_yzzzz, g_x_0_xxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyyzz_xxxx[k] = -g_x_0_xxyzz_xxxx[k] * ab_y + g_x_0_xxyzz_xxxxy[k];

                g_x_0_xxyyzz_xxxy[k] = -g_x_0_xxyzz_xxxy[k] * ab_y + g_x_0_xxyzz_xxxyy[k];

                g_x_0_xxyyzz_xxxz[k] = -g_x_0_xxyzz_xxxz[k] * ab_y + g_x_0_xxyzz_xxxyz[k];

                g_x_0_xxyyzz_xxyy[k] = -g_x_0_xxyzz_xxyy[k] * ab_y + g_x_0_xxyzz_xxyyy[k];

                g_x_0_xxyyzz_xxyz[k] = -g_x_0_xxyzz_xxyz[k] * ab_y + g_x_0_xxyzz_xxyyz[k];

                g_x_0_xxyyzz_xxzz[k] = -g_x_0_xxyzz_xxzz[k] * ab_y + g_x_0_xxyzz_xxyzz[k];

                g_x_0_xxyyzz_xyyy[k] = -g_x_0_xxyzz_xyyy[k] * ab_y + g_x_0_xxyzz_xyyyy[k];

                g_x_0_xxyyzz_xyyz[k] = -g_x_0_xxyzz_xyyz[k] * ab_y + g_x_0_xxyzz_xyyyz[k];

                g_x_0_xxyyzz_xyzz[k] = -g_x_0_xxyzz_xyzz[k] * ab_y + g_x_0_xxyzz_xyyzz[k];

                g_x_0_xxyyzz_xzzz[k] = -g_x_0_xxyzz_xzzz[k] * ab_y + g_x_0_xxyzz_xyzzz[k];

                g_x_0_xxyyzz_yyyy[k] = -g_x_0_xxyzz_yyyy[k] * ab_y + g_x_0_xxyzz_yyyyy[k];

                g_x_0_xxyyzz_yyyz[k] = -g_x_0_xxyzz_yyyz[k] * ab_y + g_x_0_xxyzz_yyyyz[k];

                g_x_0_xxyyzz_yyzz[k] = -g_x_0_xxyzz_yyzz[k] * ab_y + g_x_0_xxyzz_yyyzz[k];

                g_x_0_xxyyzz_yzzz[k] = -g_x_0_xxyzz_yzzz[k] * ab_y + g_x_0_xxyzz_yyzzz[k];

                g_x_0_xxyyzz_zzzz[k] = -g_x_0_xxyzz_zzzz[k] * ab_y + g_x_0_xxyzz_yzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxyzzz_xxxx = cbuffer.data(ig_geom_10_off + 195 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxy = cbuffer.data(ig_geom_10_off + 196 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxxz = cbuffer.data(ig_geom_10_off + 197 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxyy = cbuffer.data(ig_geom_10_off + 198 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxyz = cbuffer.data(ig_geom_10_off + 199 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xxzz = cbuffer.data(ig_geom_10_off + 200 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyyy = cbuffer.data(ig_geom_10_off + 201 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyyz = cbuffer.data(ig_geom_10_off + 202 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xyzz = cbuffer.data(ig_geom_10_off + 203 * ccomps * dcomps);

            auto g_x_0_xxyzzz_xzzz = cbuffer.data(ig_geom_10_off + 204 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyyy = cbuffer.data(ig_geom_10_off + 205 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyyz = cbuffer.data(ig_geom_10_off + 206 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yyzz = cbuffer.data(ig_geom_10_off + 207 * ccomps * dcomps);

            auto g_x_0_xxyzzz_yzzz = cbuffer.data(ig_geom_10_off + 208 * ccomps * dcomps);

            auto g_x_0_xxyzzz_zzzz = cbuffer.data(ig_geom_10_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxyzzz_xxxx, g_x_0_xxyzzz_xxxy, g_x_0_xxyzzz_xxxz, g_x_0_xxyzzz_xxyy, g_x_0_xxyzzz_xxyz, g_x_0_xxyzzz_xxzz, g_x_0_xxyzzz_xyyy, g_x_0_xxyzzz_xyyz, g_x_0_xxyzzz_xyzz, g_x_0_xxyzzz_xzzz, g_x_0_xxyzzz_yyyy, g_x_0_xxyzzz_yyyz, g_x_0_xxyzzz_yyzz, g_x_0_xxyzzz_yzzz, g_x_0_xxyzzz_zzzz, g_x_0_xxzzz_xxxx, g_x_0_xxzzz_xxxxy, g_x_0_xxzzz_xxxy, g_x_0_xxzzz_xxxyy, g_x_0_xxzzz_xxxyz, g_x_0_xxzzz_xxxz, g_x_0_xxzzz_xxyy, g_x_0_xxzzz_xxyyy, g_x_0_xxzzz_xxyyz, g_x_0_xxzzz_xxyz, g_x_0_xxzzz_xxyzz, g_x_0_xxzzz_xxzz, g_x_0_xxzzz_xyyy, g_x_0_xxzzz_xyyyy, g_x_0_xxzzz_xyyyz, g_x_0_xxzzz_xyyz, g_x_0_xxzzz_xyyzz, g_x_0_xxzzz_xyzz, g_x_0_xxzzz_xyzzz, g_x_0_xxzzz_xzzz, g_x_0_xxzzz_yyyy, g_x_0_xxzzz_yyyyy, g_x_0_xxzzz_yyyyz, g_x_0_xxzzz_yyyz, g_x_0_xxzzz_yyyzz, g_x_0_xxzzz_yyzz, g_x_0_xxzzz_yyzzz, g_x_0_xxzzz_yzzz, g_x_0_xxzzz_yzzzz, g_x_0_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxyzzz_xxxx[k] = -g_x_0_xxzzz_xxxx[k] * ab_y + g_x_0_xxzzz_xxxxy[k];

                g_x_0_xxyzzz_xxxy[k] = -g_x_0_xxzzz_xxxy[k] * ab_y + g_x_0_xxzzz_xxxyy[k];

                g_x_0_xxyzzz_xxxz[k] = -g_x_0_xxzzz_xxxz[k] * ab_y + g_x_0_xxzzz_xxxyz[k];

                g_x_0_xxyzzz_xxyy[k] = -g_x_0_xxzzz_xxyy[k] * ab_y + g_x_0_xxzzz_xxyyy[k];

                g_x_0_xxyzzz_xxyz[k] = -g_x_0_xxzzz_xxyz[k] * ab_y + g_x_0_xxzzz_xxyyz[k];

                g_x_0_xxyzzz_xxzz[k] = -g_x_0_xxzzz_xxzz[k] * ab_y + g_x_0_xxzzz_xxyzz[k];

                g_x_0_xxyzzz_xyyy[k] = -g_x_0_xxzzz_xyyy[k] * ab_y + g_x_0_xxzzz_xyyyy[k];

                g_x_0_xxyzzz_xyyz[k] = -g_x_0_xxzzz_xyyz[k] * ab_y + g_x_0_xxzzz_xyyyz[k];

                g_x_0_xxyzzz_xyzz[k] = -g_x_0_xxzzz_xyzz[k] * ab_y + g_x_0_xxzzz_xyyzz[k];

                g_x_0_xxyzzz_xzzz[k] = -g_x_0_xxzzz_xzzz[k] * ab_y + g_x_0_xxzzz_xyzzz[k];

                g_x_0_xxyzzz_yyyy[k] = -g_x_0_xxzzz_yyyy[k] * ab_y + g_x_0_xxzzz_yyyyy[k];

                g_x_0_xxyzzz_yyyz[k] = -g_x_0_xxzzz_yyyz[k] * ab_y + g_x_0_xxzzz_yyyyz[k];

                g_x_0_xxyzzz_yyzz[k] = -g_x_0_xxzzz_yyzz[k] * ab_y + g_x_0_xxzzz_yyyzz[k];

                g_x_0_xxyzzz_yzzz[k] = -g_x_0_xxzzz_yzzz[k] * ab_y + g_x_0_xxzzz_yyzzz[k];

                g_x_0_xxyzzz_zzzz[k] = -g_x_0_xxzzz_zzzz[k] * ab_y + g_x_0_xxzzz_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_x_0_xxzzzz_xxxx = cbuffer.data(ig_geom_10_off + 210 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxy = cbuffer.data(ig_geom_10_off + 211 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxxz = cbuffer.data(ig_geom_10_off + 212 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxyy = cbuffer.data(ig_geom_10_off + 213 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxyz = cbuffer.data(ig_geom_10_off + 214 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xxzz = cbuffer.data(ig_geom_10_off + 215 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyyy = cbuffer.data(ig_geom_10_off + 216 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyyz = cbuffer.data(ig_geom_10_off + 217 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xyzz = cbuffer.data(ig_geom_10_off + 218 * ccomps * dcomps);

            auto g_x_0_xxzzzz_xzzz = cbuffer.data(ig_geom_10_off + 219 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyyy = cbuffer.data(ig_geom_10_off + 220 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyyz = cbuffer.data(ig_geom_10_off + 221 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yyzz = cbuffer.data(ig_geom_10_off + 222 * ccomps * dcomps);

            auto g_x_0_xxzzzz_yzzz = cbuffer.data(ig_geom_10_off + 223 * ccomps * dcomps);

            auto g_x_0_xxzzzz_zzzz = cbuffer.data(ig_geom_10_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xxzzz_xxxx, g_x_0_xxzzz_xxxxz, g_x_0_xxzzz_xxxy, g_x_0_xxzzz_xxxyz, g_x_0_xxzzz_xxxz, g_x_0_xxzzz_xxxzz, g_x_0_xxzzz_xxyy, g_x_0_xxzzz_xxyyz, g_x_0_xxzzz_xxyz, g_x_0_xxzzz_xxyzz, g_x_0_xxzzz_xxzz, g_x_0_xxzzz_xxzzz, g_x_0_xxzzz_xyyy, g_x_0_xxzzz_xyyyz, g_x_0_xxzzz_xyyz, g_x_0_xxzzz_xyyzz, g_x_0_xxzzz_xyzz, g_x_0_xxzzz_xyzzz, g_x_0_xxzzz_xzzz, g_x_0_xxzzz_xzzzz, g_x_0_xxzzz_yyyy, g_x_0_xxzzz_yyyyz, g_x_0_xxzzz_yyyz, g_x_0_xxzzz_yyyzz, g_x_0_xxzzz_yyzz, g_x_0_xxzzz_yyzzz, g_x_0_xxzzz_yzzz, g_x_0_xxzzz_yzzzz, g_x_0_xxzzz_zzzz, g_x_0_xxzzz_zzzzz, g_x_0_xxzzzz_xxxx, g_x_0_xxzzzz_xxxy, g_x_0_xxzzzz_xxxz, g_x_0_xxzzzz_xxyy, g_x_0_xxzzzz_xxyz, g_x_0_xxzzzz_xxzz, g_x_0_xxzzzz_xyyy, g_x_0_xxzzzz_xyyz, g_x_0_xxzzzz_xyzz, g_x_0_xxzzzz_xzzz, g_x_0_xxzzzz_yyyy, g_x_0_xxzzzz_yyyz, g_x_0_xxzzzz_yyzz, g_x_0_xxzzzz_yzzz, g_x_0_xxzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xxzzzz_xxxx[k] = -g_x_0_xxzzz_xxxx[k] * ab_z + g_x_0_xxzzz_xxxxz[k];

                g_x_0_xxzzzz_xxxy[k] = -g_x_0_xxzzz_xxxy[k] * ab_z + g_x_0_xxzzz_xxxyz[k];

                g_x_0_xxzzzz_xxxz[k] = -g_x_0_xxzzz_xxxz[k] * ab_z + g_x_0_xxzzz_xxxzz[k];

                g_x_0_xxzzzz_xxyy[k] = -g_x_0_xxzzz_xxyy[k] * ab_z + g_x_0_xxzzz_xxyyz[k];

                g_x_0_xxzzzz_xxyz[k] = -g_x_0_xxzzz_xxyz[k] * ab_z + g_x_0_xxzzz_xxyzz[k];

                g_x_0_xxzzzz_xxzz[k] = -g_x_0_xxzzz_xxzz[k] * ab_z + g_x_0_xxzzz_xxzzz[k];

                g_x_0_xxzzzz_xyyy[k] = -g_x_0_xxzzz_xyyy[k] * ab_z + g_x_0_xxzzz_xyyyz[k];

                g_x_0_xxzzzz_xyyz[k] = -g_x_0_xxzzz_xyyz[k] * ab_z + g_x_0_xxzzz_xyyzz[k];

                g_x_0_xxzzzz_xyzz[k] = -g_x_0_xxzzz_xyzz[k] * ab_z + g_x_0_xxzzz_xyzzz[k];

                g_x_0_xxzzzz_xzzz[k] = -g_x_0_xxzzz_xzzz[k] * ab_z + g_x_0_xxzzz_xzzzz[k];

                g_x_0_xxzzzz_yyyy[k] = -g_x_0_xxzzz_yyyy[k] * ab_z + g_x_0_xxzzz_yyyyz[k];

                g_x_0_xxzzzz_yyyz[k] = -g_x_0_xxzzz_yyyz[k] * ab_z + g_x_0_xxzzz_yyyzz[k];

                g_x_0_xxzzzz_yyzz[k] = -g_x_0_xxzzz_yyzz[k] * ab_z + g_x_0_xxzzz_yyzzz[k];

                g_x_0_xxzzzz_yzzz[k] = -g_x_0_xxzzz_yzzz[k] * ab_z + g_x_0_xxzzz_yzzzz[k];

                g_x_0_xxzzzz_zzzz[k] = -g_x_0_xxzzz_zzzz[k] * ab_z + g_x_0_xxzzz_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 225 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 226 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 227 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 228 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 229 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 230 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 231 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 232 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 233 * ccomps * dcomps);

            auto g_x_0_xyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 234 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 235 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 236 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 237 * ccomps * dcomps);

            auto g_x_0_xyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 238 * ccomps * dcomps);

            auto g_x_0_xyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyy_xxxx, g_x_0_xyyyy_xxxxy, g_x_0_xyyyy_xxxy, g_x_0_xyyyy_xxxyy, g_x_0_xyyyy_xxxyz, g_x_0_xyyyy_xxxz, g_x_0_xyyyy_xxyy, g_x_0_xyyyy_xxyyy, g_x_0_xyyyy_xxyyz, g_x_0_xyyyy_xxyz, g_x_0_xyyyy_xxyzz, g_x_0_xyyyy_xxzz, g_x_0_xyyyy_xyyy, g_x_0_xyyyy_xyyyy, g_x_0_xyyyy_xyyyz, g_x_0_xyyyy_xyyz, g_x_0_xyyyy_xyyzz, g_x_0_xyyyy_xyzz, g_x_0_xyyyy_xyzzz, g_x_0_xyyyy_xzzz, g_x_0_xyyyy_yyyy, g_x_0_xyyyy_yyyyy, g_x_0_xyyyy_yyyyz, g_x_0_xyyyy_yyyz, g_x_0_xyyyy_yyyzz, g_x_0_xyyyy_yyzz, g_x_0_xyyyy_yyzzz, g_x_0_xyyyy_yzzz, g_x_0_xyyyy_yzzzz, g_x_0_xyyyy_zzzz, g_x_0_xyyyyy_xxxx, g_x_0_xyyyyy_xxxy, g_x_0_xyyyyy_xxxz, g_x_0_xyyyyy_xxyy, g_x_0_xyyyyy_xxyz, g_x_0_xyyyyy_xxzz, g_x_0_xyyyyy_xyyy, g_x_0_xyyyyy_xyyz, g_x_0_xyyyyy_xyzz, g_x_0_xyyyyy_xzzz, g_x_0_xyyyyy_yyyy, g_x_0_xyyyyy_yyyz, g_x_0_xyyyyy_yyzz, g_x_0_xyyyyy_yzzz, g_x_0_xyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyy_xxxx[k] = -g_x_0_xyyyy_xxxx[k] * ab_y + g_x_0_xyyyy_xxxxy[k];

                g_x_0_xyyyyy_xxxy[k] = -g_x_0_xyyyy_xxxy[k] * ab_y + g_x_0_xyyyy_xxxyy[k];

                g_x_0_xyyyyy_xxxz[k] = -g_x_0_xyyyy_xxxz[k] * ab_y + g_x_0_xyyyy_xxxyz[k];

                g_x_0_xyyyyy_xxyy[k] = -g_x_0_xyyyy_xxyy[k] * ab_y + g_x_0_xyyyy_xxyyy[k];

                g_x_0_xyyyyy_xxyz[k] = -g_x_0_xyyyy_xxyz[k] * ab_y + g_x_0_xyyyy_xxyyz[k];

                g_x_0_xyyyyy_xxzz[k] = -g_x_0_xyyyy_xxzz[k] * ab_y + g_x_0_xyyyy_xxyzz[k];

                g_x_0_xyyyyy_xyyy[k] = -g_x_0_xyyyy_xyyy[k] * ab_y + g_x_0_xyyyy_xyyyy[k];

                g_x_0_xyyyyy_xyyz[k] = -g_x_0_xyyyy_xyyz[k] * ab_y + g_x_0_xyyyy_xyyyz[k];

                g_x_0_xyyyyy_xyzz[k] = -g_x_0_xyyyy_xyzz[k] * ab_y + g_x_0_xyyyy_xyyzz[k];

                g_x_0_xyyyyy_xzzz[k] = -g_x_0_xyyyy_xzzz[k] * ab_y + g_x_0_xyyyy_xyzzz[k];

                g_x_0_xyyyyy_yyyy[k] = -g_x_0_xyyyy_yyyy[k] * ab_y + g_x_0_xyyyy_yyyyy[k];

                g_x_0_xyyyyy_yyyz[k] = -g_x_0_xyyyy_yyyz[k] * ab_y + g_x_0_xyyyy_yyyyz[k];

                g_x_0_xyyyyy_yyzz[k] = -g_x_0_xyyyy_yyzz[k] * ab_y + g_x_0_xyyyy_yyyzz[k];

                g_x_0_xyyyyy_yzzz[k] = -g_x_0_xyyyy_yzzz[k] * ab_y + g_x_0_xyyyy_yyzzz[k];

                g_x_0_xyyyyy_zzzz[k] = -g_x_0_xyyyy_zzzz[k] * ab_y + g_x_0_xyyyy_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 240 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 241 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 242 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 243 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 244 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 245 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 246 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 247 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 248 * ccomps * dcomps);

            auto g_x_0_xyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 249 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 250 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 251 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 252 * ccomps * dcomps);

            auto g_x_0_xyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 253 * ccomps * dcomps);

            auto g_x_0_xyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyyz_xxxx, g_x_0_xyyyyz_xxxy, g_x_0_xyyyyz_xxxz, g_x_0_xyyyyz_xxyy, g_x_0_xyyyyz_xxyz, g_x_0_xyyyyz_xxzz, g_x_0_xyyyyz_xyyy, g_x_0_xyyyyz_xyyz, g_x_0_xyyyyz_xyzz, g_x_0_xyyyyz_xzzz, g_x_0_xyyyyz_yyyy, g_x_0_xyyyyz_yyyz, g_x_0_xyyyyz_yyzz, g_x_0_xyyyyz_yzzz, g_x_0_xyyyyz_zzzz, g_x_0_xyyyz_xxxx, g_x_0_xyyyz_xxxxy, g_x_0_xyyyz_xxxy, g_x_0_xyyyz_xxxyy, g_x_0_xyyyz_xxxyz, g_x_0_xyyyz_xxxz, g_x_0_xyyyz_xxyy, g_x_0_xyyyz_xxyyy, g_x_0_xyyyz_xxyyz, g_x_0_xyyyz_xxyz, g_x_0_xyyyz_xxyzz, g_x_0_xyyyz_xxzz, g_x_0_xyyyz_xyyy, g_x_0_xyyyz_xyyyy, g_x_0_xyyyz_xyyyz, g_x_0_xyyyz_xyyz, g_x_0_xyyyz_xyyzz, g_x_0_xyyyz_xyzz, g_x_0_xyyyz_xyzzz, g_x_0_xyyyz_xzzz, g_x_0_xyyyz_yyyy, g_x_0_xyyyz_yyyyy, g_x_0_xyyyz_yyyyz, g_x_0_xyyyz_yyyz, g_x_0_xyyyz_yyyzz, g_x_0_xyyyz_yyzz, g_x_0_xyyyz_yyzzz, g_x_0_xyyyz_yzzz, g_x_0_xyyyz_yzzzz, g_x_0_xyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyyz_xxxx[k] = -g_x_0_xyyyz_xxxx[k] * ab_y + g_x_0_xyyyz_xxxxy[k];

                g_x_0_xyyyyz_xxxy[k] = -g_x_0_xyyyz_xxxy[k] * ab_y + g_x_0_xyyyz_xxxyy[k];

                g_x_0_xyyyyz_xxxz[k] = -g_x_0_xyyyz_xxxz[k] * ab_y + g_x_0_xyyyz_xxxyz[k];

                g_x_0_xyyyyz_xxyy[k] = -g_x_0_xyyyz_xxyy[k] * ab_y + g_x_0_xyyyz_xxyyy[k];

                g_x_0_xyyyyz_xxyz[k] = -g_x_0_xyyyz_xxyz[k] * ab_y + g_x_0_xyyyz_xxyyz[k];

                g_x_0_xyyyyz_xxzz[k] = -g_x_0_xyyyz_xxzz[k] * ab_y + g_x_0_xyyyz_xxyzz[k];

                g_x_0_xyyyyz_xyyy[k] = -g_x_0_xyyyz_xyyy[k] * ab_y + g_x_0_xyyyz_xyyyy[k];

                g_x_0_xyyyyz_xyyz[k] = -g_x_0_xyyyz_xyyz[k] * ab_y + g_x_0_xyyyz_xyyyz[k];

                g_x_0_xyyyyz_xyzz[k] = -g_x_0_xyyyz_xyzz[k] * ab_y + g_x_0_xyyyz_xyyzz[k];

                g_x_0_xyyyyz_xzzz[k] = -g_x_0_xyyyz_xzzz[k] * ab_y + g_x_0_xyyyz_xyzzz[k];

                g_x_0_xyyyyz_yyyy[k] = -g_x_0_xyyyz_yyyy[k] * ab_y + g_x_0_xyyyz_yyyyy[k];

                g_x_0_xyyyyz_yyyz[k] = -g_x_0_xyyyz_yyyz[k] * ab_y + g_x_0_xyyyz_yyyyz[k];

                g_x_0_xyyyyz_yyzz[k] = -g_x_0_xyyyz_yyzz[k] * ab_y + g_x_0_xyyyz_yyyzz[k];

                g_x_0_xyyyyz_yzzz[k] = -g_x_0_xyyyz_yzzz[k] * ab_y + g_x_0_xyyyz_yyzzz[k];

                g_x_0_xyyyyz_zzzz[k] = -g_x_0_xyyyz_zzzz[k] * ab_y + g_x_0_xyyyz_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 255 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 256 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 257 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 258 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 259 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 260 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 261 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 262 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 263 * ccomps * dcomps);

            auto g_x_0_xyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 264 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 265 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 266 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 267 * ccomps * dcomps);

            auto g_x_0_xyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 268 * ccomps * dcomps);

            auto g_x_0_xyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyyzz_xxxx, g_x_0_xyyyzz_xxxy, g_x_0_xyyyzz_xxxz, g_x_0_xyyyzz_xxyy, g_x_0_xyyyzz_xxyz, g_x_0_xyyyzz_xxzz, g_x_0_xyyyzz_xyyy, g_x_0_xyyyzz_xyyz, g_x_0_xyyyzz_xyzz, g_x_0_xyyyzz_xzzz, g_x_0_xyyyzz_yyyy, g_x_0_xyyyzz_yyyz, g_x_0_xyyyzz_yyzz, g_x_0_xyyyzz_yzzz, g_x_0_xyyyzz_zzzz, g_x_0_xyyzz_xxxx, g_x_0_xyyzz_xxxxy, g_x_0_xyyzz_xxxy, g_x_0_xyyzz_xxxyy, g_x_0_xyyzz_xxxyz, g_x_0_xyyzz_xxxz, g_x_0_xyyzz_xxyy, g_x_0_xyyzz_xxyyy, g_x_0_xyyzz_xxyyz, g_x_0_xyyzz_xxyz, g_x_0_xyyzz_xxyzz, g_x_0_xyyzz_xxzz, g_x_0_xyyzz_xyyy, g_x_0_xyyzz_xyyyy, g_x_0_xyyzz_xyyyz, g_x_0_xyyzz_xyyz, g_x_0_xyyzz_xyyzz, g_x_0_xyyzz_xyzz, g_x_0_xyyzz_xyzzz, g_x_0_xyyzz_xzzz, g_x_0_xyyzz_yyyy, g_x_0_xyyzz_yyyyy, g_x_0_xyyzz_yyyyz, g_x_0_xyyzz_yyyz, g_x_0_xyyzz_yyyzz, g_x_0_xyyzz_yyzz, g_x_0_xyyzz_yyzzz, g_x_0_xyyzz_yzzz, g_x_0_xyyzz_yzzzz, g_x_0_xyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyyzz_xxxx[k] = -g_x_0_xyyzz_xxxx[k] * ab_y + g_x_0_xyyzz_xxxxy[k];

                g_x_0_xyyyzz_xxxy[k] = -g_x_0_xyyzz_xxxy[k] * ab_y + g_x_0_xyyzz_xxxyy[k];

                g_x_0_xyyyzz_xxxz[k] = -g_x_0_xyyzz_xxxz[k] * ab_y + g_x_0_xyyzz_xxxyz[k];

                g_x_0_xyyyzz_xxyy[k] = -g_x_0_xyyzz_xxyy[k] * ab_y + g_x_0_xyyzz_xxyyy[k];

                g_x_0_xyyyzz_xxyz[k] = -g_x_0_xyyzz_xxyz[k] * ab_y + g_x_0_xyyzz_xxyyz[k];

                g_x_0_xyyyzz_xxzz[k] = -g_x_0_xyyzz_xxzz[k] * ab_y + g_x_0_xyyzz_xxyzz[k];

                g_x_0_xyyyzz_xyyy[k] = -g_x_0_xyyzz_xyyy[k] * ab_y + g_x_0_xyyzz_xyyyy[k];

                g_x_0_xyyyzz_xyyz[k] = -g_x_0_xyyzz_xyyz[k] * ab_y + g_x_0_xyyzz_xyyyz[k];

                g_x_0_xyyyzz_xyzz[k] = -g_x_0_xyyzz_xyzz[k] * ab_y + g_x_0_xyyzz_xyyzz[k];

                g_x_0_xyyyzz_xzzz[k] = -g_x_0_xyyzz_xzzz[k] * ab_y + g_x_0_xyyzz_xyzzz[k];

                g_x_0_xyyyzz_yyyy[k] = -g_x_0_xyyzz_yyyy[k] * ab_y + g_x_0_xyyzz_yyyyy[k];

                g_x_0_xyyyzz_yyyz[k] = -g_x_0_xyyzz_yyyz[k] * ab_y + g_x_0_xyyzz_yyyyz[k];

                g_x_0_xyyyzz_yyzz[k] = -g_x_0_xyyzz_yyzz[k] * ab_y + g_x_0_xyyzz_yyyzz[k];

                g_x_0_xyyyzz_yzzz[k] = -g_x_0_xyyzz_yzzz[k] * ab_y + g_x_0_xyyzz_yyzzz[k];

                g_x_0_xyyyzz_zzzz[k] = -g_x_0_xyyzz_zzzz[k] * ab_y + g_x_0_xyyzz_yzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 270 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 271 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 272 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 273 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 274 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 275 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 276 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 277 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 278 * ccomps * dcomps);

            auto g_x_0_xyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 279 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 280 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 281 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 282 * ccomps * dcomps);

            auto g_x_0_xyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 283 * ccomps * dcomps);

            auto g_x_0_xyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyyzzz_xxxx, g_x_0_xyyzzz_xxxy, g_x_0_xyyzzz_xxxz, g_x_0_xyyzzz_xxyy, g_x_0_xyyzzz_xxyz, g_x_0_xyyzzz_xxzz, g_x_0_xyyzzz_xyyy, g_x_0_xyyzzz_xyyz, g_x_0_xyyzzz_xyzz, g_x_0_xyyzzz_xzzz, g_x_0_xyyzzz_yyyy, g_x_0_xyyzzz_yyyz, g_x_0_xyyzzz_yyzz, g_x_0_xyyzzz_yzzz, g_x_0_xyyzzz_zzzz, g_x_0_xyzzz_xxxx, g_x_0_xyzzz_xxxxy, g_x_0_xyzzz_xxxy, g_x_0_xyzzz_xxxyy, g_x_0_xyzzz_xxxyz, g_x_0_xyzzz_xxxz, g_x_0_xyzzz_xxyy, g_x_0_xyzzz_xxyyy, g_x_0_xyzzz_xxyyz, g_x_0_xyzzz_xxyz, g_x_0_xyzzz_xxyzz, g_x_0_xyzzz_xxzz, g_x_0_xyzzz_xyyy, g_x_0_xyzzz_xyyyy, g_x_0_xyzzz_xyyyz, g_x_0_xyzzz_xyyz, g_x_0_xyzzz_xyyzz, g_x_0_xyzzz_xyzz, g_x_0_xyzzz_xyzzz, g_x_0_xyzzz_xzzz, g_x_0_xyzzz_yyyy, g_x_0_xyzzz_yyyyy, g_x_0_xyzzz_yyyyz, g_x_0_xyzzz_yyyz, g_x_0_xyzzz_yyyzz, g_x_0_xyzzz_yyzz, g_x_0_xyzzz_yyzzz, g_x_0_xyzzz_yzzz, g_x_0_xyzzz_yzzzz, g_x_0_xyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyyzzz_xxxx[k] = -g_x_0_xyzzz_xxxx[k] * ab_y + g_x_0_xyzzz_xxxxy[k];

                g_x_0_xyyzzz_xxxy[k] = -g_x_0_xyzzz_xxxy[k] * ab_y + g_x_0_xyzzz_xxxyy[k];

                g_x_0_xyyzzz_xxxz[k] = -g_x_0_xyzzz_xxxz[k] * ab_y + g_x_0_xyzzz_xxxyz[k];

                g_x_0_xyyzzz_xxyy[k] = -g_x_0_xyzzz_xxyy[k] * ab_y + g_x_0_xyzzz_xxyyy[k];

                g_x_0_xyyzzz_xxyz[k] = -g_x_0_xyzzz_xxyz[k] * ab_y + g_x_0_xyzzz_xxyyz[k];

                g_x_0_xyyzzz_xxzz[k] = -g_x_0_xyzzz_xxzz[k] * ab_y + g_x_0_xyzzz_xxyzz[k];

                g_x_0_xyyzzz_xyyy[k] = -g_x_0_xyzzz_xyyy[k] * ab_y + g_x_0_xyzzz_xyyyy[k];

                g_x_0_xyyzzz_xyyz[k] = -g_x_0_xyzzz_xyyz[k] * ab_y + g_x_0_xyzzz_xyyyz[k];

                g_x_0_xyyzzz_xyzz[k] = -g_x_0_xyzzz_xyzz[k] * ab_y + g_x_0_xyzzz_xyyzz[k];

                g_x_0_xyyzzz_xzzz[k] = -g_x_0_xyzzz_xzzz[k] * ab_y + g_x_0_xyzzz_xyzzz[k];

                g_x_0_xyyzzz_yyyy[k] = -g_x_0_xyzzz_yyyy[k] * ab_y + g_x_0_xyzzz_yyyyy[k];

                g_x_0_xyyzzz_yyyz[k] = -g_x_0_xyzzz_yyyz[k] * ab_y + g_x_0_xyzzz_yyyyz[k];

                g_x_0_xyyzzz_yyzz[k] = -g_x_0_xyzzz_yyzz[k] * ab_y + g_x_0_xyzzz_yyyzz[k];

                g_x_0_xyyzzz_yzzz[k] = -g_x_0_xyzzz_yzzz[k] * ab_y + g_x_0_xyzzz_yyzzz[k];

                g_x_0_xyyzzz_zzzz[k] = -g_x_0_xyzzz_zzzz[k] * ab_y + g_x_0_xyzzz_yzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_x_0_xyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 285 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 286 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 287 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 288 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 289 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 290 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 291 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 292 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 293 * ccomps * dcomps);

            auto g_x_0_xyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 294 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 295 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 296 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 297 * ccomps * dcomps);

            auto g_x_0_xyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 298 * ccomps * dcomps);

            auto g_x_0_xyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xyzzzz_xxxx, g_x_0_xyzzzz_xxxy, g_x_0_xyzzzz_xxxz, g_x_0_xyzzzz_xxyy, g_x_0_xyzzzz_xxyz, g_x_0_xyzzzz_xxzz, g_x_0_xyzzzz_xyyy, g_x_0_xyzzzz_xyyz, g_x_0_xyzzzz_xyzz, g_x_0_xyzzzz_xzzz, g_x_0_xyzzzz_yyyy, g_x_0_xyzzzz_yyyz, g_x_0_xyzzzz_yyzz, g_x_0_xyzzzz_yzzz, g_x_0_xyzzzz_zzzz, g_x_0_xzzzz_xxxx, g_x_0_xzzzz_xxxxy, g_x_0_xzzzz_xxxy, g_x_0_xzzzz_xxxyy, g_x_0_xzzzz_xxxyz, g_x_0_xzzzz_xxxz, g_x_0_xzzzz_xxyy, g_x_0_xzzzz_xxyyy, g_x_0_xzzzz_xxyyz, g_x_0_xzzzz_xxyz, g_x_0_xzzzz_xxyzz, g_x_0_xzzzz_xxzz, g_x_0_xzzzz_xyyy, g_x_0_xzzzz_xyyyy, g_x_0_xzzzz_xyyyz, g_x_0_xzzzz_xyyz, g_x_0_xzzzz_xyyzz, g_x_0_xzzzz_xyzz, g_x_0_xzzzz_xyzzz, g_x_0_xzzzz_xzzz, g_x_0_xzzzz_yyyy, g_x_0_xzzzz_yyyyy, g_x_0_xzzzz_yyyyz, g_x_0_xzzzz_yyyz, g_x_0_xzzzz_yyyzz, g_x_0_xzzzz_yyzz, g_x_0_xzzzz_yyzzz, g_x_0_xzzzz_yzzz, g_x_0_xzzzz_yzzzz, g_x_0_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xyzzzz_xxxx[k] = -g_x_0_xzzzz_xxxx[k] * ab_y + g_x_0_xzzzz_xxxxy[k];

                g_x_0_xyzzzz_xxxy[k] = -g_x_0_xzzzz_xxxy[k] * ab_y + g_x_0_xzzzz_xxxyy[k];

                g_x_0_xyzzzz_xxxz[k] = -g_x_0_xzzzz_xxxz[k] * ab_y + g_x_0_xzzzz_xxxyz[k];

                g_x_0_xyzzzz_xxyy[k] = -g_x_0_xzzzz_xxyy[k] * ab_y + g_x_0_xzzzz_xxyyy[k];

                g_x_0_xyzzzz_xxyz[k] = -g_x_0_xzzzz_xxyz[k] * ab_y + g_x_0_xzzzz_xxyyz[k];

                g_x_0_xyzzzz_xxzz[k] = -g_x_0_xzzzz_xxzz[k] * ab_y + g_x_0_xzzzz_xxyzz[k];

                g_x_0_xyzzzz_xyyy[k] = -g_x_0_xzzzz_xyyy[k] * ab_y + g_x_0_xzzzz_xyyyy[k];

                g_x_0_xyzzzz_xyyz[k] = -g_x_0_xzzzz_xyyz[k] * ab_y + g_x_0_xzzzz_xyyyz[k];

                g_x_0_xyzzzz_xyzz[k] = -g_x_0_xzzzz_xyzz[k] * ab_y + g_x_0_xzzzz_xyyzz[k];

                g_x_0_xyzzzz_xzzz[k] = -g_x_0_xzzzz_xzzz[k] * ab_y + g_x_0_xzzzz_xyzzz[k];

                g_x_0_xyzzzz_yyyy[k] = -g_x_0_xzzzz_yyyy[k] * ab_y + g_x_0_xzzzz_yyyyy[k];

                g_x_0_xyzzzz_yyyz[k] = -g_x_0_xzzzz_yyyz[k] * ab_y + g_x_0_xzzzz_yyyyz[k];

                g_x_0_xyzzzz_yyzz[k] = -g_x_0_xzzzz_yyzz[k] * ab_y + g_x_0_xzzzz_yyyzz[k];

                g_x_0_xyzzzz_yzzz[k] = -g_x_0_xzzzz_yzzz[k] * ab_y + g_x_0_xzzzz_yyzzz[k];

                g_x_0_xyzzzz_zzzz[k] = -g_x_0_xzzzz_zzzz[k] * ab_y + g_x_0_xzzzz_yzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_x_0_xzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 300 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 301 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 302 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 303 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 304 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 305 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 306 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 307 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 308 * ccomps * dcomps);

            auto g_x_0_xzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 309 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 310 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 311 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 312 * ccomps * dcomps);

            auto g_x_0_xzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 313 * ccomps * dcomps);

            auto g_x_0_xzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_xzzzz_xxxx, g_x_0_xzzzz_xxxxz, g_x_0_xzzzz_xxxy, g_x_0_xzzzz_xxxyz, g_x_0_xzzzz_xxxz, g_x_0_xzzzz_xxxzz, g_x_0_xzzzz_xxyy, g_x_0_xzzzz_xxyyz, g_x_0_xzzzz_xxyz, g_x_0_xzzzz_xxyzz, g_x_0_xzzzz_xxzz, g_x_0_xzzzz_xxzzz, g_x_0_xzzzz_xyyy, g_x_0_xzzzz_xyyyz, g_x_0_xzzzz_xyyz, g_x_0_xzzzz_xyyzz, g_x_0_xzzzz_xyzz, g_x_0_xzzzz_xyzzz, g_x_0_xzzzz_xzzz, g_x_0_xzzzz_xzzzz, g_x_0_xzzzz_yyyy, g_x_0_xzzzz_yyyyz, g_x_0_xzzzz_yyyz, g_x_0_xzzzz_yyyzz, g_x_0_xzzzz_yyzz, g_x_0_xzzzz_yyzzz, g_x_0_xzzzz_yzzz, g_x_0_xzzzz_yzzzz, g_x_0_xzzzz_zzzz, g_x_0_xzzzz_zzzzz, g_x_0_xzzzzz_xxxx, g_x_0_xzzzzz_xxxy, g_x_0_xzzzzz_xxxz, g_x_0_xzzzzz_xxyy, g_x_0_xzzzzz_xxyz, g_x_0_xzzzzz_xxzz, g_x_0_xzzzzz_xyyy, g_x_0_xzzzzz_xyyz, g_x_0_xzzzzz_xyzz, g_x_0_xzzzzz_xzzz, g_x_0_xzzzzz_yyyy, g_x_0_xzzzzz_yyyz, g_x_0_xzzzzz_yyzz, g_x_0_xzzzzz_yzzz, g_x_0_xzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_xzzzzz_xxxx[k] = -g_x_0_xzzzz_xxxx[k] * ab_z + g_x_0_xzzzz_xxxxz[k];

                g_x_0_xzzzzz_xxxy[k] = -g_x_0_xzzzz_xxxy[k] * ab_z + g_x_0_xzzzz_xxxyz[k];

                g_x_0_xzzzzz_xxxz[k] = -g_x_0_xzzzz_xxxz[k] * ab_z + g_x_0_xzzzz_xxxzz[k];

                g_x_0_xzzzzz_xxyy[k] = -g_x_0_xzzzz_xxyy[k] * ab_z + g_x_0_xzzzz_xxyyz[k];

                g_x_0_xzzzzz_xxyz[k] = -g_x_0_xzzzz_xxyz[k] * ab_z + g_x_0_xzzzz_xxyzz[k];

                g_x_0_xzzzzz_xxzz[k] = -g_x_0_xzzzz_xxzz[k] * ab_z + g_x_0_xzzzz_xxzzz[k];

                g_x_0_xzzzzz_xyyy[k] = -g_x_0_xzzzz_xyyy[k] * ab_z + g_x_0_xzzzz_xyyyz[k];

                g_x_0_xzzzzz_xyyz[k] = -g_x_0_xzzzz_xyyz[k] * ab_z + g_x_0_xzzzz_xyyzz[k];

                g_x_0_xzzzzz_xyzz[k] = -g_x_0_xzzzz_xyzz[k] * ab_z + g_x_0_xzzzz_xyzzz[k];

                g_x_0_xzzzzz_xzzz[k] = -g_x_0_xzzzz_xzzz[k] * ab_z + g_x_0_xzzzz_xzzzz[k];

                g_x_0_xzzzzz_yyyy[k] = -g_x_0_xzzzz_yyyy[k] * ab_z + g_x_0_xzzzz_yyyyz[k];

                g_x_0_xzzzzz_yyyz[k] = -g_x_0_xzzzz_yyyz[k] * ab_z + g_x_0_xzzzz_yyyzz[k];

                g_x_0_xzzzzz_yyzz[k] = -g_x_0_xzzzz_yyzz[k] * ab_z + g_x_0_xzzzz_yyzzz[k];

                g_x_0_xzzzzz_yzzz[k] = -g_x_0_xzzzz_yzzz[k] * ab_z + g_x_0_xzzzz_yzzzz[k];

                g_x_0_xzzzzz_zzzz[k] = -g_x_0_xzzzz_zzzz[k] * ab_z + g_x_0_xzzzz_zzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 315 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 316 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 317 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 318 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 319 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 320 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 321 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 322 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 323 * ccomps * dcomps);

            auto g_x_0_yyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 324 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 325 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 326 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 327 * ccomps * dcomps);

            auto g_x_0_yyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 328 * ccomps * dcomps);

            auto g_x_0_yyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyy_xxxx, g_x_0_yyyyy_xxxxy, g_x_0_yyyyy_xxxy, g_x_0_yyyyy_xxxyy, g_x_0_yyyyy_xxxyz, g_x_0_yyyyy_xxxz, g_x_0_yyyyy_xxyy, g_x_0_yyyyy_xxyyy, g_x_0_yyyyy_xxyyz, g_x_0_yyyyy_xxyz, g_x_0_yyyyy_xxyzz, g_x_0_yyyyy_xxzz, g_x_0_yyyyy_xyyy, g_x_0_yyyyy_xyyyy, g_x_0_yyyyy_xyyyz, g_x_0_yyyyy_xyyz, g_x_0_yyyyy_xyyzz, g_x_0_yyyyy_xyzz, g_x_0_yyyyy_xyzzz, g_x_0_yyyyy_xzzz, g_x_0_yyyyy_yyyy, g_x_0_yyyyy_yyyyy, g_x_0_yyyyy_yyyyz, g_x_0_yyyyy_yyyz, g_x_0_yyyyy_yyyzz, g_x_0_yyyyy_yyzz, g_x_0_yyyyy_yyzzz, g_x_0_yyyyy_yzzz, g_x_0_yyyyy_yzzzz, g_x_0_yyyyy_zzzz, g_x_0_yyyyyy_xxxx, g_x_0_yyyyyy_xxxy, g_x_0_yyyyyy_xxxz, g_x_0_yyyyyy_xxyy, g_x_0_yyyyyy_xxyz, g_x_0_yyyyyy_xxzz, g_x_0_yyyyyy_xyyy, g_x_0_yyyyyy_xyyz, g_x_0_yyyyyy_xyzz, g_x_0_yyyyyy_xzzz, g_x_0_yyyyyy_yyyy, g_x_0_yyyyyy_yyyz, g_x_0_yyyyyy_yyzz, g_x_0_yyyyyy_yzzz, g_x_0_yyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyy_xxxx[k] = -g_x_0_yyyyy_xxxx[k] * ab_y + g_x_0_yyyyy_xxxxy[k];

                g_x_0_yyyyyy_xxxy[k] = -g_x_0_yyyyy_xxxy[k] * ab_y + g_x_0_yyyyy_xxxyy[k];

                g_x_0_yyyyyy_xxxz[k] = -g_x_0_yyyyy_xxxz[k] * ab_y + g_x_0_yyyyy_xxxyz[k];

                g_x_0_yyyyyy_xxyy[k] = -g_x_0_yyyyy_xxyy[k] * ab_y + g_x_0_yyyyy_xxyyy[k];

                g_x_0_yyyyyy_xxyz[k] = -g_x_0_yyyyy_xxyz[k] * ab_y + g_x_0_yyyyy_xxyyz[k];

                g_x_0_yyyyyy_xxzz[k] = -g_x_0_yyyyy_xxzz[k] * ab_y + g_x_0_yyyyy_xxyzz[k];

                g_x_0_yyyyyy_xyyy[k] = -g_x_0_yyyyy_xyyy[k] * ab_y + g_x_0_yyyyy_xyyyy[k];

                g_x_0_yyyyyy_xyyz[k] = -g_x_0_yyyyy_xyyz[k] * ab_y + g_x_0_yyyyy_xyyyz[k];

                g_x_0_yyyyyy_xyzz[k] = -g_x_0_yyyyy_xyzz[k] * ab_y + g_x_0_yyyyy_xyyzz[k];

                g_x_0_yyyyyy_xzzz[k] = -g_x_0_yyyyy_xzzz[k] * ab_y + g_x_0_yyyyy_xyzzz[k];

                g_x_0_yyyyyy_yyyy[k] = -g_x_0_yyyyy_yyyy[k] * ab_y + g_x_0_yyyyy_yyyyy[k];

                g_x_0_yyyyyy_yyyz[k] = -g_x_0_yyyyy_yyyz[k] * ab_y + g_x_0_yyyyy_yyyyz[k];

                g_x_0_yyyyyy_yyzz[k] = -g_x_0_yyyyy_yyzz[k] * ab_y + g_x_0_yyyyy_yyyzz[k];

                g_x_0_yyyyyy_yzzz[k] = -g_x_0_yyyyy_yzzz[k] * ab_y + g_x_0_yyyyy_yyzzz[k];

                g_x_0_yyyyyy_zzzz[k] = -g_x_0_yyyyy_zzzz[k] * ab_y + g_x_0_yyyyy_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 330 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 331 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 332 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 333 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 334 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 335 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 336 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 337 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 338 * ccomps * dcomps);

            auto g_x_0_yyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 339 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 340 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 341 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 342 * ccomps * dcomps);

            auto g_x_0_yyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 343 * ccomps * dcomps);

            auto g_x_0_yyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyyz_xxxx, g_x_0_yyyyyz_xxxy, g_x_0_yyyyyz_xxxz, g_x_0_yyyyyz_xxyy, g_x_0_yyyyyz_xxyz, g_x_0_yyyyyz_xxzz, g_x_0_yyyyyz_xyyy, g_x_0_yyyyyz_xyyz, g_x_0_yyyyyz_xyzz, g_x_0_yyyyyz_xzzz, g_x_0_yyyyyz_yyyy, g_x_0_yyyyyz_yyyz, g_x_0_yyyyyz_yyzz, g_x_0_yyyyyz_yzzz, g_x_0_yyyyyz_zzzz, g_x_0_yyyyz_xxxx, g_x_0_yyyyz_xxxxy, g_x_0_yyyyz_xxxy, g_x_0_yyyyz_xxxyy, g_x_0_yyyyz_xxxyz, g_x_0_yyyyz_xxxz, g_x_0_yyyyz_xxyy, g_x_0_yyyyz_xxyyy, g_x_0_yyyyz_xxyyz, g_x_0_yyyyz_xxyz, g_x_0_yyyyz_xxyzz, g_x_0_yyyyz_xxzz, g_x_0_yyyyz_xyyy, g_x_0_yyyyz_xyyyy, g_x_0_yyyyz_xyyyz, g_x_0_yyyyz_xyyz, g_x_0_yyyyz_xyyzz, g_x_0_yyyyz_xyzz, g_x_0_yyyyz_xyzzz, g_x_0_yyyyz_xzzz, g_x_0_yyyyz_yyyy, g_x_0_yyyyz_yyyyy, g_x_0_yyyyz_yyyyz, g_x_0_yyyyz_yyyz, g_x_0_yyyyz_yyyzz, g_x_0_yyyyz_yyzz, g_x_0_yyyyz_yyzzz, g_x_0_yyyyz_yzzz, g_x_0_yyyyz_yzzzz, g_x_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyyz_xxxx[k] = -g_x_0_yyyyz_xxxx[k] * ab_y + g_x_0_yyyyz_xxxxy[k];

                g_x_0_yyyyyz_xxxy[k] = -g_x_0_yyyyz_xxxy[k] * ab_y + g_x_0_yyyyz_xxxyy[k];

                g_x_0_yyyyyz_xxxz[k] = -g_x_0_yyyyz_xxxz[k] * ab_y + g_x_0_yyyyz_xxxyz[k];

                g_x_0_yyyyyz_xxyy[k] = -g_x_0_yyyyz_xxyy[k] * ab_y + g_x_0_yyyyz_xxyyy[k];

                g_x_0_yyyyyz_xxyz[k] = -g_x_0_yyyyz_xxyz[k] * ab_y + g_x_0_yyyyz_xxyyz[k];

                g_x_0_yyyyyz_xxzz[k] = -g_x_0_yyyyz_xxzz[k] * ab_y + g_x_0_yyyyz_xxyzz[k];

                g_x_0_yyyyyz_xyyy[k] = -g_x_0_yyyyz_xyyy[k] * ab_y + g_x_0_yyyyz_xyyyy[k];

                g_x_0_yyyyyz_xyyz[k] = -g_x_0_yyyyz_xyyz[k] * ab_y + g_x_0_yyyyz_xyyyz[k];

                g_x_0_yyyyyz_xyzz[k] = -g_x_0_yyyyz_xyzz[k] * ab_y + g_x_0_yyyyz_xyyzz[k];

                g_x_0_yyyyyz_xzzz[k] = -g_x_0_yyyyz_xzzz[k] * ab_y + g_x_0_yyyyz_xyzzz[k];

                g_x_0_yyyyyz_yyyy[k] = -g_x_0_yyyyz_yyyy[k] * ab_y + g_x_0_yyyyz_yyyyy[k];

                g_x_0_yyyyyz_yyyz[k] = -g_x_0_yyyyz_yyyz[k] * ab_y + g_x_0_yyyyz_yyyyz[k];

                g_x_0_yyyyyz_yyzz[k] = -g_x_0_yyyyz_yyzz[k] * ab_y + g_x_0_yyyyz_yyyzz[k];

                g_x_0_yyyyyz_yzzz[k] = -g_x_0_yyyyz_yzzz[k] * ab_y + g_x_0_yyyyz_yyzzz[k];

                g_x_0_yyyyyz_zzzz[k] = -g_x_0_yyyyz_zzzz[k] * ab_y + g_x_0_yyyyz_yzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 345 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 346 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 347 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 348 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 349 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 350 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 351 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 352 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 353 * ccomps * dcomps);

            auto g_x_0_yyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 354 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 355 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 356 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 357 * ccomps * dcomps);

            auto g_x_0_yyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 358 * ccomps * dcomps);

            auto g_x_0_yyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyyzz_xxxx, g_x_0_yyyyzz_xxxy, g_x_0_yyyyzz_xxxz, g_x_0_yyyyzz_xxyy, g_x_0_yyyyzz_xxyz, g_x_0_yyyyzz_xxzz, g_x_0_yyyyzz_xyyy, g_x_0_yyyyzz_xyyz, g_x_0_yyyyzz_xyzz, g_x_0_yyyyzz_xzzz, g_x_0_yyyyzz_yyyy, g_x_0_yyyyzz_yyyz, g_x_0_yyyyzz_yyzz, g_x_0_yyyyzz_yzzz, g_x_0_yyyyzz_zzzz, g_x_0_yyyzz_xxxx, g_x_0_yyyzz_xxxxy, g_x_0_yyyzz_xxxy, g_x_0_yyyzz_xxxyy, g_x_0_yyyzz_xxxyz, g_x_0_yyyzz_xxxz, g_x_0_yyyzz_xxyy, g_x_0_yyyzz_xxyyy, g_x_0_yyyzz_xxyyz, g_x_0_yyyzz_xxyz, g_x_0_yyyzz_xxyzz, g_x_0_yyyzz_xxzz, g_x_0_yyyzz_xyyy, g_x_0_yyyzz_xyyyy, g_x_0_yyyzz_xyyyz, g_x_0_yyyzz_xyyz, g_x_0_yyyzz_xyyzz, g_x_0_yyyzz_xyzz, g_x_0_yyyzz_xyzzz, g_x_0_yyyzz_xzzz, g_x_0_yyyzz_yyyy, g_x_0_yyyzz_yyyyy, g_x_0_yyyzz_yyyyz, g_x_0_yyyzz_yyyz, g_x_0_yyyzz_yyyzz, g_x_0_yyyzz_yyzz, g_x_0_yyyzz_yyzzz, g_x_0_yyyzz_yzzz, g_x_0_yyyzz_yzzzz, g_x_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyyzz_xxxx[k] = -g_x_0_yyyzz_xxxx[k] * ab_y + g_x_0_yyyzz_xxxxy[k];

                g_x_0_yyyyzz_xxxy[k] = -g_x_0_yyyzz_xxxy[k] * ab_y + g_x_0_yyyzz_xxxyy[k];

                g_x_0_yyyyzz_xxxz[k] = -g_x_0_yyyzz_xxxz[k] * ab_y + g_x_0_yyyzz_xxxyz[k];

                g_x_0_yyyyzz_xxyy[k] = -g_x_0_yyyzz_xxyy[k] * ab_y + g_x_0_yyyzz_xxyyy[k];

                g_x_0_yyyyzz_xxyz[k] = -g_x_0_yyyzz_xxyz[k] * ab_y + g_x_0_yyyzz_xxyyz[k];

                g_x_0_yyyyzz_xxzz[k] = -g_x_0_yyyzz_xxzz[k] * ab_y + g_x_0_yyyzz_xxyzz[k];

                g_x_0_yyyyzz_xyyy[k] = -g_x_0_yyyzz_xyyy[k] * ab_y + g_x_0_yyyzz_xyyyy[k];

                g_x_0_yyyyzz_xyyz[k] = -g_x_0_yyyzz_xyyz[k] * ab_y + g_x_0_yyyzz_xyyyz[k];

                g_x_0_yyyyzz_xyzz[k] = -g_x_0_yyyzz_xyzz[k] * ab_y + g_x_0_yyyzz_xyyzz[k];

                g_x_0_yyyyzz_xzzz[k] = -g_x_0_yyyzz_xzzz[k] * ab_y + g_x_0_yyyzz_xyzzz[k];

                g_x_0_yyyyzz_yyyy[k] = -g_x_0_yyyzz_yyyy[k] * ab_y + g_x_0_yyyzz_yyyyy[k];

                g_x_0_yyyyzz_yyyz[k] = -g_x_0_yyyzz_yyyz[k] * ab_y + g_x_0_yyyzz_yyyyz[k];

                g_x_0_yyyyzz_yyzz[k] = -g_x_0_yyyzz_yyzz[k] * ab_y + g_x_0_yyyzz_yyyzz[k];

                g_x_0_yyyyzz_yzzz[k] = -g_x_0_yyyzz_yzzz[k] * ab_y + g_x_0_yyyzz_yyzzz[k];

                g_x_0_yyyyzz_zzzz[k] = -g_x_0_yyyzz_zzzz[k] * ab_y + g_x_0_yyyzz_yzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 360 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 361 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 362 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 363 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 364 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 365 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 366 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 367 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 368 * ccomps * dcomps);

            auto g_x_0_yyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 369 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 370 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 371 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 372 * ccomps * dcomps);

            auto g_x_0_yyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 373 * ccomps * dcomps);

            auto g_x_0_yyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyyzzz_xxxx, g_x_0_yyyzzz_xxxy, g_x_0_yyyzzz_xxxz, g_x_0_yyyzzz_xxyy, g_x_0_yyyzzz_xxyz, g_x_0_yyyzzz_xxzz, g_x_0_yyyzzz_xyyy, g_x_0_yyyzzz_xyyz, g_x_0_yyyzzz_xyzz, g_x_0_yyyzzz_xzzz, g_x_0_yyyzzz_yyyy, g_x_0_yyyzzz_yyyz, g_x_0_yyyzzz_yyzz, g_x_0_yyyzzz_yzzz, g_x_0_yyyzzz_zzzz, g_x_0_yyzzz_xxxx, g_x_0_yyzzz_xxxxy, g_x_0_yyzzz_xxxy, g_x_0_yyzzz_xxxyy, g_x_0_yyzzz_xxxyz, g_x_0_yyzzz_xxxz, g_x_0_yyzzz_xxyy, g_x_0_yyzzz_xxyyy, g_x_0_yyzzz_xxyyz, g_x_0_yyzzz_xxyz, g_x_0_yyzzz_xxyzz, g_x_0_yyzzz_xxzz, g_x_0_yyzzz_xyyy, g_x_0_yyzzz_xyyyy, g_x_0_yyzzz_xyyyz, g_x_0_yyzzz_xyyz, g_x_0_yyzzz_xyyzz, g_x_0_yyzzz_xyzz, g_x_0_yyzzz_xyzzz, g_x_0_yyzzz_xzzz, g_x_0_yyzzz_yyyy, g_x_0_yyzzz_yyyyy, g_x_0_yyzzz_yyyyz, g_x_0_yyzzz_yyyz, g_x_0_yyzzz_yyyzz, g_x_0_yyzzz_yyzz, g_x_0_yyzzz_yyzzz, g_x_0_yyzzz_yzzz, g_x_0_yyzzz_yzzzz, g_x_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyyzzz_xxxx[k] = -g_x_0_yyzzz_xxxx[k] * ab_y + g_x_0_yyzzz_xxxxy[k];

                g_x_0_yyyzzz_xxxy[k] = -g_x_0_yyzzz_xxxy[k] * ab_y + g_x_0_yyzzz_xxxyy[k];

                g_x_0_yyyzzz_xxxz[k] = -g_x_0_yyzzz_xxxz[k] * ab_y + g_x_0_yyzzz_xxxyz[k];

                g_x_0_yyyzzz_xxyy[k] = -g_x_0_yyzzz_xxyy[k] * ab_y + g_x_0_yyzzz_xxyyy[k];

                g_x_0_yyyzzz_xxyz[k] = -g_x_0_yyzzz_xxyz[k] * ab_y + g_x_0_yyzzz_xxyyz[k];

                g_x_0_yyyzzz_xxzz[k] = -g_x_0_yyzzz_xxzz[k] * ab_y + g_x_0_yyzzz_xxyzz[k];

                g_x_0_yyyzzz_xyyy[k] = -g_x_0_yyzzz_xyyy[k] * ab_y + g_x_0_yyzzz_xyyyy[k];

                g_x_0_yyyzzz_xyyz[k] = -g_x_0_yyzzz_xyyz[k] * ab_y + g_x_0_yyzzz_xyyyz[k];

                g_x_0_yyyzzz_xyzz[k] = -g_x_0_yyzzz_xyzz[k] * ab_y + g_x_0_yyzzz_xyyzz[k];

                g_x_0_yyyzzz_xzzz[k] = -g_x_0_yyzzz_xzzz[k] * ab_y + g_x_0_yyzzz_xyzzz[k];

                g_x_0_yyyzzz_yyyy[k] = -g_x_0_yyzzz_yyyy[k] * ab_y + g_x_0_yyzzz_yyyyy[k];

                g_x_0_yyyzzz_yyyz[k] = -g_x_0_yyzzz_yyyz[k] * ab_y + g_x_0_yyzzz_yyyyz[k];

                g_x_0_yyyzzz_yyzz[k] = -g_x_0_yyzzz_yyzz[k] * ab_y + g_x_0_yyzzz_yyyzz[k];

                g_x_0_yyyzzz_yzzz[k] = -g_x_0_yyzzz_yzzz[k] * ab_y + g_x_0_yyzzz_yyzzz[k];

                g_x_0_yyyzzz_zzzz[k] = -g_x_0_yyzzz_zzzz[k] * ab_y + g_x_0_yyzzz_yzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_x_0_yyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 375 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 376 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 377 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 378 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 379 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 380 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 381 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 382 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 383 * ccomps * dcomps);

            auto g_x_0_yyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 384 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 385 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 386 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 387 * ccomps * dcomps);

            auto g_x_0_yyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 388 * ccomps * dcomps);

            auto g_x_0_yyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yyzzzz_xxxx, g_x_0_yyzzzz_xxxy, g_x_0_yyzzzz_xxxz, g_x_0_yyzzzz_xxyy, g_x_0_yyzzzz_xxyz, g_x_0_yyzzzz_xxzz, g_x_0_yyzzzz_xyyy, g_x_0_yyzzzz_xyyz, g_x_0_yyzzzz_xyzz, g_x_0_yyzzzz_xzzz, g_x_0_yyzzzz_yyyy, g_x_0_yyzzzz_yyyz, g_x_0_yyzzzz_yyzz, g_x_0_yyzzzz_yzzz, g_x_0_yyzzzz_zzzz, g_x_0_yzzzz_xxxx, g_x_0_yzzzz_xxxxy, g_x_0_yzzzz_xxxy, g_x_0_yzzzz_xxxyy, g_x_0_yzzzz_xxxyz, g_x_0_yzzzz_xxxz, g_x_0_yzzzz_xxyy, g_x_0_yzzzz_xxyyy, g_x_0_yzzzz_xxyyz, g_x_0_yzzzz_xxyz, g_x_0_yzzzz_xxyzz, g_x_0_yzzzz_xxzz, g_x_0_yzzzz_xyyy, g_x_0_yzzzz_xyyyy, g_x_0_yzzzz_xyyyz, g_x_0_yzzzz_xyyz, g_x_0_yzzzz_xyyzz, g_x_0_yzzzz_xyzz, g_x_0_yzzzz_xyzzz, g_x_0_yzzzz_xzzz, g_x_0_yzzzz_yyyy, g_x_0_yzzzz_yyyyy, g_x_0_yzzzz_yyyyz, g_x_0_yzzzz_yyyz, g_x_0_yzzzz_yyyzz, g_x_0_yzzzz_yyzz, g_x_0_yzzzz_yyzzz, g_x_0_yzzzz_yzzz, g_x_0_yzzzz_yzzzz, g_x_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yyzzzz_xxxx[k] = -g_x_0_yzzzz_xxxx[k] * ab_y + g_x_0_yzzzz_xxxxy[k];

                g_x_0_yyzzzz_xxxy[k] = -g_x_0_yzzzz_xxxy[k] * ab_y + g_x_0_yzzzz_xxxyy[k];

                g_x_0_yyzzzz_xxxz[k] = -g_x_0_yzzzz_xxxz[k] * ab_y + g_x_0_yzzzz_xxxyz[k];

                g_x_0_yyzzzz_xxyy[k] = -g_x_0_yzzzz_xxyy[k] * ab_y + g_x_0_yzzzz_xxyyy[k];

                g_x_0_yyzzzz_xxyz[k] = -g_x_0_yzzzz_xxyz[k] * ab_y + g_x_0_yzzzz_xxyyz[k];

                g_x_0_yyzzzz_xxzz[k] = -g_x_0_yzzzz_xxzz[k] * ab_y + g_x_0_yzzzz_xxyzz[k];

                g_x_0_yyzzzz_xyyy[k] = -g_x_0_yzzzz_xyyy[k] * ab_y + g_x_0_yzzzz_xyyyy[k];

                g_x_0_yyzzzz_xyyz[k] = -g_x_0_yzzzz_xyyz[k] * ab_y + g_x_0_yzzzz_xyyyz[k];

                g_x_0_yyzzzz_xyzz[k] = -g_x_0_yzzzz_xyzz[k] * ab_y + g_x_0_yzzzz_xyyzz[k];

                g_x_0_yyzzzz_xzzz[k] = -g_x_0_yzzzz_xzzz[k] * ab_y + g_x_0_yzzzz_xyzzz[k];

                g_x_0_yyzzzz_yyyy[k] = -g_x_0_yzzzz_yyyy[k] * ab_y + g_x_0_yzzzz_yyyyy[k];

                g_x_0_yyzzzz_yyyz[k] = -g_x_0_yzzzz_yyyz[k] * ab_y + g_x_0_yzzzz_yyyyz[k];

                g_x_0_yyzzzz_yyzz[k] = -g_x_0_yzzzz_yyzz[k] * ab_y + g_x_0_yzzzz_yyyzz[k];

                g_x_0_yyzzzz_yzzz[k] = -g_x_0_yzzzz_yzzz[k] * ab_y + g_x_0_yzzzz_yyzzz[k];

                g_x_0_yyzzzz_zzzz[k] = -g_x_0_yzzzz_zzzz[k] * ab_y + g_x_0_yzzzz_yzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_x_0_yzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 390 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 391 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 392 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 393 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 394 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 395 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 396 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 397 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 398 * ccomps * dcomps);

            auto g_x_0_yzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 399 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 400 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 401 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 402 * ccomps * dcomps);

            auto g_x_0_yzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 403 * ccomps * dcomps);

            auto g_x_0_yzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_yzzzzz_xxxx, g_x_0_yzzzzz_xxxy, g_x_0_yzzzzz_xxxz, g_x_0_yzzzzz_xxyy, g_x_0_yzzzzz_xxyz, g_x_0_yzzzzz_xxzz, g_x_0_yzzzzz_xyyy, g_x_0_yzzzzz_xyyz, g_x_0_yzzzzz_xyzz, g_x_0_yzzzzz_xzzz, g_x_0_yzzzzz_yyyy, g_x_0_yzzzzz_yyyz, g_x_0_yzzzzz_yyzz, g_x_0_yzzzzz_yzzz, g_x_0_yzzzzz_zzzz, g_x_0_zzzzz_xxxx, g_x_0_zzzzz_xxxxy, g_x_0_zzzzz_xxxy, g_x_0_zzzzz_xxxyy, g_x_0_zzzzz_xxxyz, g_x_0_zzzzz_xxxz, g_x_0_zzzzz_xxyy, g_x_0_zzzzz_xxyyy, g_x_0_zzzzz_xxyyz, g_x_0_zzzzz_xxyz, g_x_0_zzzzz_xxyzz, g_x_0_zzzzz_xxzz, g_x_0_zzzzz_xyyy, g_x_0_zzzzz_xyyyy, g_x_0_zzzzz_xyyyz, g_x_0_zzzzz_xyyz, g_x_0_zzzzz_xyyzz, g_x_0_zzzzz_xyzz, g_x_0_zzzzz_xyzzz, g_x_0_zzzzz_xzzz, g_x_0_zzzzz_yyyy, g_x_0_zzzzz_yyyyy, g_x_0_zzzzz_yyyyz, g_x_0_zzzzz_yyyz, g_x_0_zzzzz_yyyzz, g_x_0_zzzzz_yyzz, g_x_0_zzzzz_yyzzz, g_x_0_zzzzz_yzzz, g_x_0_zzzzz_yzzzz, g_x_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_yzzzzz_xxxx[k] = -g_x_0_zzzzz_xxxx[k] * ab_y + g_x_0_zzzzz_xxxxy[k];

                g_x_0_yzzzzz_xxxy[k] = -g_x_0_zzzzz_xxxy[k] * ab_y + g_x_0_zzzzz_xxxyy[k];

                g_x_0_yzzzzz_xxxz[k] = -g_x_0_zzzzz_xxxz[k] * ab_y + g_x_0_zzzzz_xxxyz[k];

                g_x_0_yzzzzz_xxyy[k] = -g_x_0_zzzzz_xxyy[k] * ab_y + g_x_0_zzzzz_xxyyy[k];

                g_x_0_yzzzzz_xxyz[k] = -g_x_0_zzzzz_xxyz[k] * ab_y + g_x_0_zzzzz_xxyyz[k];

                g_x_0_yzzzzz_xxzz[k] = -g_x_0_zzzzz_xxzz[k] * ab_y + g_x_0_zzzzz_xxyzz[k];

                g_x_0_yzzzzz_xyyy[k] = -g_x_0_zzzzz_xyyy[k] * ab_y + g_x_0_zzzzz_xyyyy[k];

                g_x_0_yzzzzz_xyyz[k] = -g_x_0_zzzzz_xyyz[k] * ab_y + g_x_0_zzzzz_xyyyz[k];

                g_x_0_yzzzzz_xyzz[k] = -g_x_0_zzzzz_xyzz[k] * ab_y + g_x_0_zzzzz_xyyzz[k];

                g_x_0_yzzzzz_xzzz[k] = -g_x_0_zzzzz_xzzz[k] * ab_y + g_x_0_zzzzz_xyzzz[k];

                g_x_0_yzzzzz_yyyy[k] = -g_x_0_zzzzz_yyyy[k] * ab_y + g_x_0_zzzzz_yyyyy[k];

                g_x_0_yzzzzz_yyyz[k] = -g_x_0_zzzzz_yyyz[k] * ab_y + g_x_0_zzzzz_yyyyz[k];

                g_x_0_yzzzzz_yyzz[k] = -g_x_0_zzzzz_yyzz[k] * ab_y + g_x_0_zzzzz_yyyzz[k];

                g_x_0_yzzzzz_yzzz[k] = -g_x_0_zzzzz_yzzz[k] * ab_y + g_x_0_zzzzz_yyzzz[k];

                g_x_0_yzzzzz_zzzz[k] = -g_x_0_zzzzz_zzzz[k] * ab_y + g_x_0_zzzzz_yzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_x_0_zzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 405 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 406 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 407 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 408 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 409 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 410 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 411 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 412 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 413 * ccomps * dcomps);

            auto g_x_0_zzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 414 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 415 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 416 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 417 * ccomps * dcomps);

            auto g_x_0_zzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 418 * ccomps * dcomps);

            auto g_x_0_zzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_x_0_zzzzz_xxxx, g_x_0_zzzzz_xxxxz, g_x_0_zzzzz_xxxy, g_x_0_zzzzz_xxxyz, g_x_0_zzzzz_xxxz, g_x_0_zzzzz_xxxzz, g_x_0_zzzzz_xxyy, g_x_0_zzzzz_xxyyz, g_x_0_zzzzz_xxyz, g_x_0_zzzzz_xxyzz, g_x_0_zzzzz_xxzz, g_x_0_zzzzz_xxzzz, g_x_0_zzzzz_xyyy, g_x_0_zzzzz_xyyyz, g_x_0_zzzzz_xyyz, g_x_0_zzzzz_xyyzz, g_x_0_zzzzz_xyzz, g_x_0_zzzzz_xyzzz, g_x_0_zzzzz_xzzz, g_x_0_zzzzz_xzzzz, g_x_0_zzzzz_yyyy, g_x_0_zzzzz_yyyyz, g_x_0_zzzzz_yyyz, g_x_0_zzzzz_yyyzz, g_x_0_zzzzz_yyzz, g_x_0_zzzzz_yyzzz, g_x_0_zzzzz_yzzz, g_x_0_zzzzz_yzzzz, g_x_0_zzzzz_zzzz, g_x_0_zzzzz_zzzzz, g_x_0_zzzzzz_xxxx, g_x_0_zzzzzz_xxxy, g_x_0_zzzzzz_xxxz, g_x_0_zzzzzz_xxyy, g_x_0_zzzzzz_xxyz, g_x_0_zzzzzz_xxzz, g_x_0_zzzzzz_xyyy, g_x_0_zzzzzz_xyyz, g_x_0_zzzzzz_xyzz, g_x_0_zzzzzz_xzzz, g_x_0_zzzzzz_yyyy, g_x_0_zzzzzz_yyyz, g_x_0_zzzzzz_yyzz, g_x_0_zzzzzz_yzzz, g_x_0_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_x_0_zzzzzz_xxxx[k] = -g_x_0_zzzzz_xxxx[k] * ab_z + g_x_0_zzzzz_xxxxz[k];

                g_x_0_zzzzzz_xxxy[k] = -g_x_0_zzzzz_xxxy[k] * ab_z + g_x_0_zzzzz_xxxyz[k];

                g_x_0_zzzzzz_xxxz[k] = -g_x_0_zzzzz_xxxz[k] * ab_z + g_x_0_zzzzz_xxxzz[k];

                g_x_0_zzzzzz_xxyy[k] = -g_x_0_zzzzz_xxyy[k] * ab_z + g_x_0_zzzzz_xxyyz[k];

                g_x_0_zzzzzz_xxyz[k] = -g_x_0_zzzzz_xxyz[k] * ab_z + g_x_0_zzzzz_xxyzz[k];

                g_x_0_zzzzzz_xxzz[k] = -g_x_0_zzzzz_xxzz[k] * ab_z + g_x_0_zzzzz_xxzzz[k];

                g_x_0_zzzzzz_xyyy[k] = -g_x_0_zzzzz_xyyy[k] * ab_z + g_x_0_zzzzz_xyyyz[k];

                g_x_0_zzzzzz_xyyz[k] = -g_x_0_zzzzz_xyyz[k] * ab_z + g_x_0_zzzzz_xyyzz[k];

                g_x_0_zzzzzz_xyzz[k] = -g_x_0_zzzzz_xyzz[k] * ab_z + g_x_0_zzzzz_xyzzz[k];

                g_x_0_zzzzzz_xzzz[k] = -g_x_0_zzzzz_xzzz[k] * ab_z + g_x_0_zzzzz_xzzzz[k];

                g_x_0_zzzzzz_yyyy[k] = -g_x_0_zzzzz_yyyy[k] * ab_z + g_x_0_zzzzz_yyyyz[k];

                g_x_0_zzzzzz_yyyz[k] = -g_x_0_zzzzz_yyyz[k] * ab_z + g_x_0_zzzzz_yyyzz[k];

                g_x_0_zzzzzz_yyzz[k] = -g_x_0_zzzzz_yyzz[k] * ab_z + g_x_0_zzzzz_yyzzz[k];

                g_x_0_zzzzzz_yzzz[k] = -g_x_0_zzzzz_yzzz[k] * ab_z + g_x_0_zzzzz_yzzzz[k];

                g_x_0_zzzzzz_zzzz[k] = -g_x_0_zzzzz_zzzz[k] * ab_z + g_x_0_zzzzz_zzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxx_xxxx = cbuffer.data(ig_geom_10_off + 420 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxy = cbuffer.data(ig_geom_10_off + 421 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxxz = cbuffer.data(ig_geom_10_off + 422 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxyy = cbuffer.data(ig_geom_10_off + 423 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxyz = cbuffer.data(ig_geom_10_off + 424 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xxzz = cbuffer.data(ig_geom_10_off + 425 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyyy = cbuffer.data(ig_geom_10_off + 426 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyyz = cbuffer.data(ig_geom_10_off + 427 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xyzz = cbuffer.data(ig_geom_10_off + 428 * ccomps * dcomps);

            auto g_y_0_xxxxxx_xzzz = cbuffer.data(ig_geom_10_off + 429 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyyy = cbuffer.data(ig_geom_10_off + 430 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyyz = cbuffer.data(ig_geom_10_off + 431 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yyzz = cbuffer.data(ig_geom_10_off + 432 * ccomps * dcomps);

            auto g_y_0_xxxxxx_yzzz = cbuffer.data(ig_geom_10_off + 433 * ccomps * dcomps);

            auto g_y_0_xxxxxx_zzzz = cbuffer.data(ig_geom_10_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxx_xxxx, g_y_0_xxxxx_xxxxx, g_y_0_xxxxx_xxxxy, g_y_0_xxxxx_xxxxz, g_y_0_xxxxx_xxxy, g_y_0_xxxxx_xxxyy, g_y_0_xxxxx_xxxyz, g_y_0_xxxxx_xxxz, g_y_0_xxxxx_xxxzz, g_y_0_xxxxx_xxyy, g_y_0_xxxxx_xxyyy, g_y_0_xxxxx_xxyyz, g_y_0_xxxxx_xxyz, g_y_0_xxxxx_xxyzz, g_y_0_xxxxx_xxzz, g_y_0_xxxxx_xxzzz, g_y_0_xxxxx_xyyy, g_y_0_xxxxx_xyyyy, g_y_0_xxxxx_xyyyz, g_y_0_xxxxx_xyyz, g_y_0_xxxxx_xyyzz, g_y_0_xxxxx_xyzz, g_y_0_xxxxx_xyzzz, g_y_0_xxxxx_xzzz, g_y_0_xxxxx_xzzzz, g_y_0_xxxxx_yyyy, g_y_0_xxxxx_yyyz, g_y_0_xxxxx_yyzz, g_y_0_xxxxx_yzzz, g_y_0_xxxxx_zzzz, g_y_0_xxxxxx_xxxx, g_y_0_xxxxxx_xxxy, g_y_0_xxxxxx_xxxz, g_y_0_xxxxxx_xxyy, g_y_0_xxxxxx_xxyz, g_y_0_xxxxxx_xxzz, g_y_0_xxxxxx_xyyy, g_y_0_xxxxxx_xyyz, g_y_0_xxxxxx_xyzz, g_y_0_xxxxxx_xzzz, g_y_0_xxxxxx_yyyy, g_y_0_xxxxxx_yyyz, g_y_0_xxxxxx_yyzz, g_y_0_xxxxxx_yzzz, g_y_0_xxxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxx_xxxx[k] = -g_y_0_xxxxx_xxxx[k] * ab_x + g_y_0_xxxxx_xxxxx[k];

                g_y_0_xxxxxx_xxxy[k] = -g_y_0_xxxxx_xxxy[k] * ab_x + g_y_0_xxxxx_xxxxy[k];

                g_y_0_xxxxxx_xxxz[k] = -g_y_0_xxxxx_xxxz[k] * ab_x + g_y_0_xxxxx_xxxxz[k];

                g_y_0_xxxxxx_xxyy[k] = -g_y_0_xxxxx_xxyy[k] * ab_x + g_y_0_xxxxx_xxxyy[k];

                g_y_0_xxxxxx_xxyz[k] = -g_y_0_xxxxx_xxyz[k] * ab_x + g_y_0_xxxxx_xxxyz[k];

                g_y_0_xxxxxx_xxzz[k] = -g_y_0_xxxxx_xxzz[k] * ab_x + g_y_0_xxxxx_xxxzz[k];

                g_y_0_xxxxxx_xyyy[k] = -g_y_0_xxxxx_xyyy[k] * ab_x + g_y_0_xxxxx_xxyyy[k];

                g_y_0_xxxxxx_xyyz[k] = -g_y_0_xxxxx_xyyz[k] * ab_x + g_y_0_xxxxx_xxyyz[k];

                g_y_0_xxxxxx_xyzz[k] = -g_y_0_xxxxx_xyzz[k] * ab_x + g_y_0_xxxxx_xxyzz[k];

                g_y_0_xxxxxx_xzzz[k] = -g_y_0_xxxxx_xzzz[k] * ab_x + g_y_0_xxxxx_xxzzz[k];

                g_y_0_xxxxxx_yyyy[k] = -g_y_0_xxxxx_yyyy[k] * ab_x + g_y_0_xxxxx_xyyyy[k];

                g_y_0_xxxxxx_yyyz[k] = -g_y_0_xxxxx_yyyz[k] * ab_x + g_y_0_xxxxx_xyyyz[k];

                g_y_0_xxxxxx_yyzz[k] = -g_y_0_xxxxx_yyzz[k] * ab_x + g_y_0_xxxxx_xyyzz[k];

                g_y_0_xxxxxx_yzzz[k] = -g_y_0_xxxxx_yzzz[k] * ab_x + g_y_0_xxxxx_xyzzz[k];

                g_y_0_xxxxxx_zzzz[k] = -g_y_0_xxxxx_zzzz[k] * ab_x + g_y_0_xxxxx_xzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxy_xxxx = cbuffer.data(ig_geom_10_off + 435 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxy = cbuffer.data(ig_geom_10_off + 436 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxxz = cbuffer.data(ig_geom_10_off + 437 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxyy = cbuffer.data(ig_geom_10_off + 438 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxyz = cbuffer.data(ig_geom_10_off + 439 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xxzz = cbuffer.data(ig_geom_10_off + 440 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyyy = cbuffer.data(ig_geom_10_off + 441 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyyz = cbuffer.data(ig_geom_10_off + 442 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xyzz = cbuffer.data(ig_geom_10_off + 443 * ccomps * dcomps);

            auto g_y_0_xxxxxy_xzzz = cbuffer.data(ig_geom_10_off + 444 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyyy = cbuffer.data(ig_geom_10_off + 445 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyyz = cbuffer.data(ig_geom_10_off + 446 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yyzz = cbuffer.data(ig_geom_10_off + 447 * ccomps * dcomps);

            auto g_y_0_xxxxxy_yzzz = cbuffer.data(ig_geom_10_off + 448 * ccomps * dcomps);

            auto g_y_0_xxxxxy_zzzz = cbuffer.data(ig_geom_10_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxxy_xxxx, g_y_0_xxxxxy_xxxy, g_y_0_xxxxxy_xxxz, g_y_0_xxxxxy_xxyy, g_y_0_xxxxxy_xxyz, g_y_0_xxxxxy_xxzz, g_y_0_xxxxxy_xyyy, g_y_0_xxxxxy_xyyz, g_y_0_xxxxxy_xyzz, g_y_0_xxxxxy_xzzz, g_y_0_xxxxxy_yyyy, g_y_0_xxxxxy_yyyz, g_y_0_xxxxxy_yyzz, g_y_0_xxxxxy_yzzz, g_y_0_xxxxxy_zzzz, g_y_0_xxxxy_xxxx, g_y_0_xxxxy_xxxxx, g_y_0_xxxxy_xxxxy, g_y_0_xxxxy_xxxxz, g_y_0_xxxxy_xxxy, g_y_0_xxxxy_xxxyy, g_y_0_xxxxy_xxxyz, g_y_0_xxxxy_xxxz, g_y_0_xxxxy_xxxzz, g_y_0_xxxxy_xxyy, g_y_0_xxxxy_xxyyy, g_y_0_xxxxy_xxyyz, g_y_0_xxxxy_xxyz, g_y_0_xxxxy_xxyzz, g_y_0_xxxxy_xxzz, g_y_0_xxxxy_xxzzz, g_y_0_xxxxy_xyyy, g_y_0_xxxxy_xyyyy, g_y_0_xxxxy_xyyyz, g_y_0_xxxxy_xyyz, g_y_0_xxxxy_xyyzz, g_y_0_xxxxy_xyzz, g_y_0_xxxxy_xyzzz, g_y_0_xxxxy_xzzz, g_y_0_xxxxy_xzzzz, g_y_0_xxxxy_yyyy, g_y_0_xxxxy_yyyz, g_y_0_xxxxy_yyzz, g_y_0_xxxxy_yzzz, g_y_0_xxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxy_xxxx[k] = -g_y_0_xxxxy_xxxx[k] * ab_x + g_y_0_xxxxy_xxxxx[k];

                g_y_0_xxxxxy_xxxy[k] = -g_y_0_xxxxy_xxxy[k] * ab_x + g_y_0_xxxxy_xxxxy[k];

                g_y_0_xxxxxy_xxxz[k] = -g_y_0_xxxxy_xxxz[k] * ab_x + g_y_0_xxxxy_xxxxz[k];

                g_y_0_xxxxxy_xxyy[k] = -g_y_0_xxxxy_xxyy[k] * ab_x + g_y_0_xxxxy_xxxyy[k];

                g_y_0_xxxxxy_xxyz[k] = -g_y_0_xxxxy_xxyz[k] * ab_x + g_y_0_xxxxy_xxxyz[k];

                g_y_0_xxxxxy_xxzz[k] = -g_y_0_xxxxy_xxzz[k] * ab_x + g_y_0_xxxxy_xxxzz[k];

                g_y_0_xxxxxy_xyyy[k] = -g_y_0_xxxxy_xyyy[k] * ab_x + g_y_0_xxxxy_xxyyy[k];

                g_y_0_xxxxxy_xyyz[k] = -g_y_0_xxxxy_xyyz[k] * ab_x + g_y_0_xxxxy_xxyyz[k];

                g_y_0_xxxxxy_xyzz[k] = -g_y_0_xxxxy_xyzz[k] * ab_x + g_y_0_xxxxy_xxyzz[k];

                g_y_0_xxxxxy_xzzz[k] = -g_y_0_xxxxy_xzzz[k] * ab_x + g_y_0_xxxxy_xxzzz[k];

                g_y_0_xxxxxy_yyyy[k] = -g_y_0_xxxxy_yyyy[k] * ab_x + g_y_0_xxxxy_xyyyy[k];

                g_y_0_xxxxxy_yyyz[k] = -g_y_0_xxxxy_yyyz[k] * ab_x + g_y_0_xxxxy_xyyyz[k];

                g_y_0_xxxxxy_yyzz[k] = -g_y_0_xxxxy_yyzz[k] * ab_x + g_y_0_xxxxy_xyyzz[k];

                g_y_0_xxxxxy_yzzz[k] = -g_y_0_xxxxy_yzzz[k] * ab_x + g_y_0_xxxxy_xyzzz[k];

                g_y_0_xxxxxy_zzzz[k] = -g_y_0_xxxxy_zzzz[k] * ab_x + g_y_0_xxxxy_xzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxxz_xxxx = cbuffer.data(ig_geom_10_off + 450 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxy = cbuffer.data(ig_geom_10_off + 451 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxxz = cbuffer.data(ig_geom_10_off + 452 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxyy = cbuffer.data(ig_geom_10_off + 453 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxyz = cbuffer.data(ig_geom_10_off + 454 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xxzz = cbuffer.data(ig_geom_10_off + 455 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyyy = cbuffer.data(ig_geom_10_off + 456 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyyz = cbuffer.data(ig_geom_10_off + 457 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xyzz = cbuffer.data(ig_geom_10_off + 458 * ccomps * dcomps);

            auto g_y_0_xxxxxz_xzzz = cbuffer.data(ig_geom_10_off + 459 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyyy = cbuffer.data(ig_geom_10_off + 460 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyyz = cbuffer.data(ig_geom_10_off + 461 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yyzz = cbuffer.data(ig_geom_10_off + 462 * ccomps * dcomps);

            auto g_y_0_xxxxxz_yzzz = cbuffer.data(ig_geom_10_off + 463 * ccomps * dcomps);

            auto g_y_0_xxxxxz_zzzz = cbuffer.data(ig_geom_10_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxxz_xxxx, g_y_0_xxxxxz_xxxy, g_y_0_xxxxxz_xxxz, g_y_0_xxxxxz_xxyy, g_y_0_xxxxxz_xxyz, g_y_0_xxxxxz_xxzz, g_y_0_xxxxxz_xyyy, g_y_0_xxxxxz_xyyz, g_y_0_xxxxxz_xyzz, g_y_0_xxxxxz_xzzz, g_y_0_xxxxxz_yyyy, g_y_0_xxxxxz_yyyz, g_y_0_xxxxxz_yyzz, g_y_0_xxxxxz_yzzz, g_y_0_xxxxxz_zzzz, g_y_0_xxxxz_xxxx, g_y_0_xxxxz_xxxxx, g_y_0_xxxxz_xxxxy, g_y_0_xxxxz_xxxxz, g_y_0_xxxxz_xxxy, g_y_0_xxxxz_xxxyy, g_y_0_xxxxz_xxxyz, g_y_0_xxxxz_xxxz, g_y_0_xxxxz_xxxzz, g_y_0_xxxxz_xxyy, g_y_0_xxxxz_xxyyy, g_y_0_xxxxz_xxyyz, g_y_0_xxxxz_xxyz, g_y_0_xxxxz_xxyzz, g_y_0_xxxxz_xxzz, g_y_0_xxxxz_xxzzz, g_y_0_xxxxz_xyyy, g_y_0_xxxxz_xyyyy, g_y_0_xxxxz_xyyyz, g_y_0_xxxxz_xyyz, g_y_0_xxxxz_xyyzz, g_y_0_xxxxz_xyzz, g_y_0_xxxxz_xyzzz, g_y_0_xxxxz_xzzz, g_y_0_xxxxz_xzzzz, g_y_0_xxxxz_yyyy, g_y_0_xxxxz_yyyz, g_y_0_xxxxz_yyzz, g_y_0_xxxxz_yzzz, g_y_0_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxxz_xxxx[k] = -g_y_0_xxxxz_xxxx[k] * ab_x + g_y_0_xxxxz_xxxxx[k];

                g_y_0_xxxxxz_xxxy[k] = -g_y_0_xxxxz_xxxy[k] * ab_x + g_y_0_xxxxz_xxxxy[k];

                g_y_0_xxxxxz_xxxz[k] = -g_y_0_xxxxz_xxxz[k] * ab_x + g_y_0_xxxxz_xxxxz[k];

                g_y_0_xxxxxz_xxyy[k] = -g_y_0_xxxxz_xxyy[k] * ab_x + g_y_0_xxxxz_xxxyy[k];

                g_y_0_xxxxxz_xxyz[k] = -g_y_0_xxxxz_xxyz[k] * ab_x + g_y_0_xxxxz_xxxyz[k];

                g_y_0_xxxxxz_xxzz[k] = -g_y_0_xxxxz_xxzz[k] * ab_x + g_y_0_xxxxz_xxxzz[k];

                g_y_0_xxxxxz_xyyy[k] = -g_y_0_xxxxz_xyyy[k] * ab_x + g_y_0_xxxxz_xxyyy[k];

                g_y_0_xxxxxz_xyyz[k] = -g_y_0_xxxxz_xyyz[k] * ab_x + g_y_0_xxxxz_xxyyz[k];

                g_y_0_xxxxxz_xyzz[k] = -g_y_0_xxxxz_xyzz[k] * ab_x + g_y_0_xxxxz_xxyzz[k];

                g_y_0_xxxxxz_xzzz[k] = -g_y_0_xxxxz_xzzz[k] * ab_x + g_y_0_xxxxz_xxzzz[k];

                g_y_0_xxxxxz_yyyy[k] = -g_y_0_xxxxz_yyyy[k] * ab_x + g_y_0_xxxxz_xyyyy[k];

                g_y_0_xxxxxz_yyyz[k] = -g_y_0_xxxxz_yyyz[k] * ab_x + g_y_0_xxxxz_xyyyz[k];

                g_y_0_xxxxxz_yyzz[k] = -g_y_0_xxxxz_yyzz[k] * ab_x + g_y_0_xxxxz_xyyzz[k];

                g_y_0_xxxxxz_yzzz[k] = -g_y_0_xxxxz_yzzz[k] * ab_x + g_y_0_xxxxz_xyzzz[k];

                g_y_0_xxxxxz_zzzz[k] = -g_y_0_xxxxz_zzzz[k] * ab_x + g_y_0_xxxxz_xzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyy_xxxx = cbuffer.data(ig_geom_10_off + 465 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxy = cbuffer.data(ig_geom_10_off + 466 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxxz = cbuffer.data(ig_geom_10_off + 467 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxyy = cbuffer.data(ig_geom_10_off + 468 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxyz = cbuffer.data(ig_geom_10_off + 469 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xxzz = cbuffer.data(ig_geom_10_off + 470 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyyy = cbuffer.data(ig_geom_10_off + 471 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyyz = cbuffer.data(ig_geom_10_off + 472 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xyzz = cbuffer.data(ig_geom_10_off + 473 * ccomps * dcomps);

            auto g_y_0_xxxxyy_xzzz = cbuffer.data(ig_geom_10_off + 474 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyyy = cbuffer.data(ig_geom_10_off + 475 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyyz = cbuffer.data(ig_geom_10_off + 476 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yyzz = cbuffer.data(ig_geom_10_off + 477 * ccomps * dcomps);

            auto g_y_0_xxxxyy_yzzz = cbuffer.data(ig_geom_10_off + 478 * ccomps * dcomps);

            auto g_y_0_xxxxyy_zzzz = cbuffer.data(ig_geom_10_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxyy_xxxx, g_y_0_xxxxyy_xxxy, g_y_0_xxxxyy_xxxz, g_y_0_xxxxyy_xxyy, g_y_0_xxxxyy_xxyz, g_y_0_xxxxyy_xxzz, g_y_0_xxxxyy_xyyy, g_y_0_xxxxyy_xyyz, g_y_0_xxxxyy_xyzz, g_y_0_xxxxyy_xzzz, g_y_0_xxxxyy_yyyy, g_y_0_xxxxyy_yyyz, g_y_0_xxxxyy_yyzz, g_y_0_xxxxyy_yzzz, g_y_0_xxxxyy_zzzz, g_y_0_xxxyy_xxxx, g_y_0_xxxyy_xxxxx, g_y_0_xxxyy_xxxxy, g_y_0_xxxyy_xxxxz, g_y_0_xxxyy_xxxy, g_y_0_xxxyy_xxxyy, g_y_0_xxxyy_xxxyz, g_y_0_xxxyy_xxxz, g_y_0_xxxyy_xxxzz, g_y_0_xxxyy_xxyy, g_y_0_xxxyy_xxyyy, g_y_0_xxxyy_xxyyz, g_y_0_xxxyy_xxyz, g_y_0_xxxyy_xxyzz, g_y_0_xxxyy_xxzz, g_y_0_xxxyy_xxzzz, g_y_0_xxxyy_xyyy, g_y_0_xxxyy_xyyyy, g_y_0_xxxyy_xyyyz, g_y_0_xxxyy_xyyz, g_y_0_xxxyy_xyyzz, g_y_0_xxxyy_xyzz, g_y_0_xxxyy_xyzzz, g_y_0_xxxyy_xzzz, g_y_0_xxxyy_xzzzz, g_y_0_xxxyy_yyyy, g_y_0_xxxyy_yyyz, g_y_0_xxxyy_yyzz, g_y_0_xxxyy_yzzz, g_y_0_xxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyy_xxxx[k] = -g_y_0_xxxyy_xxxx[k] * ab_x + g_y_0_xxxyy_xxxxx[k];

                g_y_0_xxxxyy_xxxy[k] = -g_y_0_xxxyy_xxxy[k] * ab_x + g_y_0_xxxyy_xxxxy[k];

                g_y_0_xxxxyy_xxxz[k] = -g_y_0_xxxyy_xxxz[k] * ab_x + g_y_0_xxxyy_xxxxz[k];

                g_y_0_xxxxyy_xxyy[k] = -g_y_0_xxxyy_xxyy[k] * ab_x + g_y_0_xxxyy_xxxyy[k];

                g_y_0_xxxxyy_xxyz[k] = -g_y_0_xxxyy_xxyz[k] * ab_x + g_y_0_xxxyy_xxxyz[k];

                g_y_0_xxxxyy_xxzz[k] = -g_y_0_xxxyy_xxzz[k] * ab_x + g_y_0_xxxyy_xxxzz[k];

                g_y_0_xxxxyy_xyyy[k] = -g_y_0_xxxyy_xyyy[k] * ab_x + g_y_0_xxxyy_xxyyy[k];

                g_y_0_xxxxyy_xyyz[k] = -g_y_0_xxxyy_xyyz[k] * ab_x + g_y_0_xxxyy_xxyyz[k];

                g_y_0_xxxxyy_xyzz[k] = -g_y_0_xxxyy_xyzz[k] * ab_x + g_y_0_xxxyy_xxyzz[k];

                g_y_0_xxxxyy_xzzz[k] = -g_y_0_xxxyy_xzzz[k] * ab_x + g_y_0_xxxyy_xxzzz[k];

                g_y_0_xxxxyy_yyyy[k] = -g_y_0_xxxyy_yyyy[k] * ab_x + g_y_0_xxxyy_xyyyy[k];

                g_y_0_xxxxyy_yyyz[k] = -g_y_0_xxxyy_yyyz[k] * ab_x + g_y_0_xxxyy_xyyyz[k];

                g_y_0_xxxxyy_yyzz[k] = -g_y_0_xxxyy_yyzz[k] * ab_x + g_y_0_xxxyy_xyyzz[k];

                g_y_0_xxxxyy_yzzz[k] = -g_y_0_xxxyy_yzzz[k] * ab_x + g_y_0_xxxyy_xyzzz[k];

                g_y_0_xxxxyy_zzzz[k] = -g_y_0_xxxyy_zzzz[k] * ab_x + g_y_0_xxxyy_xzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxyz_xxxx = cbuffer.data(ig_geom_10_off + 480 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxy = cbuffer.data(ig_geom_10_off + 481 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxxz = cbuffer.data(ig_geom_10_off + 482 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxyy = cbuffer.data(ig_geom_10_off + 483 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxyz = cbuffer.data(ig_geom_10_off + 484 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xxzz = cbuffer.data(ig_geom_10_off + 485 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyyy = cbuffer.data(ig_geom_10_off + 486 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyyz = cbuffer.data(ig_geom_10_off + 487 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xyzz = cbuffer.data(ig_geom_10_off + 488 * ccomps * dcomps);

            auto g_y_0_xxxxyz_xzzz = cbuffer.data(ig_geom_10_off + 489 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyyy = cbuffer.data(ig_geom_10_off + 490 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyyz = cbuffer.data(ig_geom_10_off + 491 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yyzz = cbuffer.data(ig_geom_10_off + 492 * ccomps * dcomps);

            auto g_y_0_xxxxyz_yzzz = cbuffer.data(ig_geom_10_off + 493 * ccomps * dcomps);

            auto g_y_0_xxxxyz_zzzz = cbuffer.data(ig_geom_10_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxyz_xxxx, g_y_0_xxxxyz_xxxy, g_y_0_xxxxyz_xxxz, g_y_0_xxxxyz_xxyy, g_y_0_xxxxyz_xxyz, g_y_0_xxxxyz_xxzz, g_y_0_xxxxyz_xyyy, g_y_0_xxxxyz_xyyz, g_y_0_xxxxyz_xyzz, g_y_0_xxxxyz_xzzz, g_y_0_xxxxyz_yyyy, g_y_0_xxxxyz_yyyz, g_y_0_xxxxyz_yyzz, g_y_0_xxxxyz_yzzz, g_y_0_xxxxyz_zzzz, g_y_0_xxxyz_xxxx, g_y_0_xxxyz_xxxxx, g_y_0_xxxyz_xxxxy, g_y_0_xxxyz_xxxxz, g_y_0_xxxyz_xxxy, g_y_0_xxxyz_xxxyy, g_y_0_xxxyz_xxxyz, g_y_0_xxxyz_xxxz, g_y_0_xxxyz_xxxzz, g_y_0_xxxyz_xxyy, g_y_0_xxxyz_xxyyy, g_y_0_xxxyz_xxyyz, g_y_0_xxxyz_xxyz, g_y_0_xxxyz_xxyzz, g_y_0_xxxyz_xxzz, g_y_0_xxxyz_xxzzz, g_y_0_xxxyz_xyyy, g_y_0_xxxyz_xyyyy, g_y_0_xxxyz_xyyyz, g_y_0_xxxyz_xyyz, g_y_0_xxxyz_xyyzz, g_y_0_xxxyz_xyzz, g_y_0_xxxyz_xyzzz, g_y_0_xxxyz_xzzz, g_y_0_xxxyz_xzzzz, g_y_0_xxxyz_yyyy, g_y_0_xxxyz_yyyz, g_y_0_xxxyz_yyzz, g_y_0_xxxyz_yzzz, g_y_0_xxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxyz_xxxx[k] = -g_y_0_xxxyz_xxxx[k] * ab_x + g_y_0_xxxyz_xxxxx[k];

                g_y_0_xxxxyz_xxxy[k] = -g_y_0_xxxyz_xxxy[k] * ab_x + g_y_0_xxxyz_xxxxy[k];

                g_y_0_xxxxyz_xxxz[k] = -g_y_0_xxxyz_xxxz[k] * ab_x + g_y_0_xxxyz_xxxxz[k];

                g_y_0_xxxxyz_xxyy[k] = -g_y_0_xxxyz_xxyy[k] * ab_x + g_y_0_xxxyz_xxxyy[k];

                g_y_0_xxxxyz_xxyz[k] = -g_y_0_xxxyz_xxyz[k] * ab_x + g_y_0_xxxyz_xxxyz[k];

                g_y_0_xxxxyz_xxzz[k] = -g_y_0_xxxyz_xxzz[k] * ab_x + g_y_0_xxxyz_xxxzz[k];

                g_y_0_xxxxyz_xyyy[k] = -g_y_0_xxxyz_xyyy[k] * ab_x + g_y_0_xxxyz_xxyyy[k];

                g_y_0_xxxxyz_xyyz[k] = -g_y_0_xxxyz_xyyz[k] * ab_x + g_y_0_xxxyz_xxyyz[k];

                g_y_0_xxxxyz_xyzz[k] = -g_y_0_xxxyz_xyzz[k] * ab_x + g_y_0_xxxyz_xxyzz[k];

                g_y_0_xxxxyz_xzzz[k] = -g_y_0_xxxyz_xzzz[k] * ab_x + g_y_0_xxxyz_xxzzz[k];

                g_y_0_xxxxyz_yyyy[k] = -g_y_0_xxxyz_yyyy[k] * ab_x + g_y_0_xxxyz_xyyyy[k];

                g_y_0_xxxxyz_yyyz[k] = -g_y_0_xxxyz_yyyz[k] * ab_x + g_y_0_xxxyz_xyyyz[k];

                g_y_0_xxxxyz_yyzz[k] = -g_y_0_xxxyz_yyzz[k] * ab_x + g_y_0_xxxyz_xyyzz[k];

                g_y_0_xxxxyz_yzzz[k] = -g_y_0_xxxyz_yzzz[k] * ab_x + g_y_0_xxxyz_xyzzz[k];

                g_y_0_xxxxyz_zzzz[k] = -g_y_0_xxxyz_zzzz[k] * ab_x + g_y_0_xxxyz_xzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxxzz_xxxx = cbuffer.data(ig_geom_10_off + 495 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxy = cbuffer.data(ig_geom_10_off + 496 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxxz = cbuffer.data(ig_geom_10_off + 497 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxyy = cbuffer.data(ig_geom_10_off + 498 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxyz = cbuffer.data(ig_geom_10_off + 499 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xxzz = cbuffer.data(ig_geom_10_off + 500 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyyy = cbuffer.data(ig_geom_10_off + 501 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyyz = cbuffer.data(ig_geom_10_off + 502 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xyzz = cbuffer.data(ig_geom_10_off + 503 * ccomps * dcomps);

            auto g_y_0_xxxxzz_xzzz = cbuffer.data(ig_geom_10_off + 504 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyyy = cbuffer.data(ig_geom_10_off + 505 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyyz = cbuffer.data(ig_geom_10_off + 506 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yyzz = cbuffer.data(ig_geom_10_off + 507 * ccomps * dcomps);

            auto g_y_0_xxxxzz_yzzz = cbuffer.data(ig_geom_10_off + 508 * ccomps * dcomps);

            auto g_y_0_xxxxzz_zzzz = cbuffer.data(ig_geom_10_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxxzz_xxxx, g_y_0_xxxxzz_xxxy, g_y_0_xxxxzz_xxxz, g_y_0_xxxxzz_xxyy, g_y_0_xxxxzz_xxyz, g_y_0_xxxxzz_xxzz, g_y_0_xxxxzz_xyyy, g_y_0_xxxxzz_xyyz, g_y_0_xxxxzz_xyzz, g_y_0_xxxxzz_xzzz, g_y_0_xxxxzz_yyyy, g_y_0_xxxxzz_yyyz, g_y_0_xxxxzz_yyzz, g_y_0_xxxxzz_yzzz, g_y_0_xxxxzz_zzzz, g_y_0_xxxzz_xxxx, g_y_0_xxxzz_xxxxx, g_y_0_xxxzz_xxxxy, g_y_0_xxxzz_xxxxz, g_y_0_xxxzz_xxxy, g_y_0_xxxzz_xxxyy, g_y_0_xxxzz_xxxyz, g_y_0_xxxzz_xxxz, g_y_0_xxxzz_xxxzz, g_y_0_xxxzz_xxyy, g_y_0_xxxzz_xxyyy, g_y_0_xxxzz_xxyyz, g_y_0_xxxzz_xxyz, g_y_0_xxxzz_xxyzz, g_y_0_xxxzz_xxzz, g_y_0_xxxzz_xxzzz, g_y_0_xxxzz_xyyy, g_y_0_xxxzz_xyyyy, g_y_0_xxxzz_xyyyz, g_y_0_xxxzz_xyyz, g_y_0_xxxzz_xyyzz, g_y_0_xxxzz_xyzz, g_y_0_xxxzz_xyzzz, g_y_0_xxxzz_xzzz, g_y_0_xxxzz_xzzzz, g_y_0_xxxzz_yyyy, g_y_0_xxxzz_yyyz, g_y_0_xxxzz_yyzz, g_y_0_xxxzz_yzzz, g_y_0_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxxzz_xxxx[k] = -g_y_0_xxxzz_xxxx[k] * ab_x + g_y_0_xxxzz_xxxxx[k];

                g_y_0_xxxxzz_xxxy[k] = -g_y_0_xxxzz_xxxy[k] * ab_x + g_y_0_xxxzz_xxxxy[k];

                g_y_0_xxxxzz_xxxz[k] = -g_y_0_xxxzz_xxxz[k] * ab_x + g_y_0_xxxzz_xxxxz[k];

                g_y_0_xxxxzz_xxyy[k] = -g_y_0_xxxzz_xxyy[k] * ab_x + g_y_0_xxxzz_xxxyy[k];

                g_y_0_xxxxzz_xxyz[k] = -g_y_0_xxxzz_xxyz[k] * ab_x + g_y_0_xxxzz_xxxyz[k];

                g_y_0_xxxxzz_xxzz[k] = -g_y_0_xxxzz_xxzz[k] * ab_x + g_y_0_xxxzz_xxxzz[k];

                g_y_0_xxxxzz_xyyy[k] = -g_y_0_xxxzz_xyyy[k] * ab_x + g_y_0_xxxzz_xxyyy[k];

                g_y_0_xxxxzz_xyyz[k] = -g_y_0_xxxzz_xyyz[k] * ab_x + g_y_0_xxxzz_xxyyz[k];

                g_y_0_xxxxzz_xyzz[k] = -g_y_0_xxxzz_xyzz[k] * ab_x + g_y_0_xxxzz_xxyzz[k];

                g_y_0_xxxxzz_xzzz[k] = -g_y_0_xxxzz_xzzz[k] * ab_x + g_y_0_xxxzz_xxzzz[k];

                g_y_0_xxxxzz_yyyy[k] = -g_y_0_xxxzz_yyyy[k] * ab_x + g_y_0_xxxzz_xyyyy[k];

                g_y_0_xxxxzz_yyyz[k] = -g_y_0_xxxzz_yyyz[k] * ab_x + g_y_0_xxxzz_xyyyz[k];

                g_y_0_xxxxzz_yyzz[k] = -g_y_0_xxxzz_yyzz[k] * ab_x + g_y_0_xxxzz_xyyzz[k];

                g_y_0_xxxxzz_yzzz[k] = -g_y_0_xxxzz_yzzz[k] * ab_x + g_y_0_xxxzz_xyzzz[k];

                g_y_0_xxxxzz_zzzz[k] = -g_y_0_xxxzz_zzzz[k] * ab_x + g_y_0_xxxzz_xzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyy_xxxx = cbuffer.data(ig_geom_10_off + 510 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxy = cbuffer.data(ig_geom_10_off + 511 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxxz = cbuffer.data(ig_geom_10_off + 512 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxyy = cbuffer.data(ig_geom_10_off + 513 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxyz = cbuffer.data(ig_geom_10_off + 514 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xxzz = cbuffer.data(ig_geom_10_off + 515 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyyy = cbuffer.data(ig_geom_10_off + 516 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyyz = cbuffer.data(ig_geom_10_off + 517 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xyzz = cbuffer.data(ig_geom_10_off + 518 * ccomps * dcomps);

            auto g_y_0_xxxyyy_xzzz = cbuffer.data(ig_geom_10_off + 519 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyyy = cbuffer.data(ig_geom_10_off + 520 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyyz = cbuffer.data(ig_geom_10_off + 521 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yyzz = cbuffer.data(ig_geom_10_off + 522 * ccomps * dcomps);

            auto g_y_0_xxxyyy_yzzz = cbuffer.data(ig_geom_10_off + 523 * ccomps * dcomps);

            auto g_y_0_xxxyyy_zzzz = cbuffer.data(ig_geom_10_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyyy_xxxx, g_y_0_xxxyyy_xxxy, g_y_0_xxxyyy_xxxz, g_y_0_xxxyyy_xxyy, g_y_0_xxxyyy_xxyz, g_y_0_xxxyyy_xxzz, g_y_0_xxxyyy_xyyy, g_y_0_xxxyyy_xyyz, g_y_0_xxxyyy_xyzz, g_y_0_xxxyyy_xzzz, g_y_0_xxxyyy_yyyy, g_y_0_xxxyyy_yyyz, g_y_0_xxxyyy_yyzz, g_y_0_xxxyyy_yzzz, g_y_0_xxxyyy_zzzz, g_y_0_xxyyy_xxxx, g_y_0_xxyyy_xxxxx, g_y_0_xxyyy_xxxxy, g_y_0_xxyyy_xxxxz, g_y_0_xxyyy_xxxy, g_y_0_xxyyy_xxxyy, g_y_0_xxyyy_xxxyz, g_y_0_xxyyy_xxxz, g_y_0_xxyyy_xxxzz, g_y_0_xxyyy_xxyy, g_y_0_xxyyy_xxyyy, g_y_0_xxyyy_xxyyz, g_y_0_xxyyy_xxyz, g_y_0_xxyyy_xxyzz, g_y_0_xxyyy_xxzz, g_y_0_xxyyy_xxzzz, g_y_0_xxyyy_xyyy, g_y_0_xxyyy_xyyyy, g_y_0_xxyyy_xyyyz, g_y_0_xxyyy_xyyz, g_y_0_xxyyy_xyyzz, g_y_0_xxyyy_xyzz, g_y_0_xxyyy_xyzzz, g_y_0_xxyyy_xzzz, g_y_0_xxyyy_xzzzz, g_y_0_xxyyy_yyyy, g_y_0_xxyyy_yyyz, g_y_0_xxyyy_yyzz, g_y_0_xxyyy_yzzz, g_y_0_xxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyy_xxxx[k] = -g_y_0_xxyyy_xxxx[k] * ab_x + g_y_0_xxyyy_xxxxx[k];

                g_y_0_xxxyyy_xxxy[k] = -g_y_0_xxyyy_xxxy[k] * ab_x + g_y_0_xxyyy_xxxxy[k];

                g_y_0_xxxyyy_xxxz[k] = -g_y_0_xxyyy_xxxz[k] * ab_x + g_y_0_xxyyy_xxxxz[k];

                g_y_0_xxxyyy_xxyy[k] = -g_y_0_xxyyy_xxyy[k] * ab_x + g_y_0_xxyyy_xxxyy[k];

                g_y_0_xxxyyy_xxyz[k] = -g_y_0_xxyyy_xxyz[k] * ab_x + g_y_0_xxyyy_xxxyz[k];

                g_y_0_xxxyyy_xxzz[k] = -g_y_0_xxyyy_xxzz[k] * ab_x + g_y_0_xxyyy_xxxzz[k];

                g_y_0_xxxyyy_xyyy[k] = -g_y_0_xxyyy_xyyy[k] * ab_x + g_y_0_xxyyy_xxyyy[k];

                g_y_0_xxxyyy_xyyz[k] = -g_y_0_xxyyy_xyyz[k] * ab_x + g_y_0_xxyyy_xxyyz[k];

                g_y_0_xxxyyy_xyzz[k] = -g_y_0_xxyyy_xyzz[k] * ab_x + g_y_0_xxyyy_xxyzz[k];

                g_y_0_xxxyyy_xzzz[k] = -g_y_0_xxyyy_xzzz[k] * ab_x + g_y_0_xxyyy_xxzzz[k];

                g_y_0_xxxyyy_yyyy[k] = -g_y_0_xxyyy_yyyy[k] * ab_x + g_y_0_xxyyy_xyyyy[k];

                g_y_0_xxxyyy_yyyz[k] = -g_y_0_xxyyy_yyyz[k] * ab_x + g_y_0_xxyyy_xyyyz[k];

                g_y_0_xxxyyy_yyzz[k] = -g_y_0_xxyyy_yyzz[k] * ab_x + g_y_0_xxyyy_xyyzz[k];

                g_y_0_xxxyyy_yzzz[k] = -g_y_0_xxyyy_yzzz[k] * ab_x + g_y_0_xxyyy_xyzzz[k];

                g_y_0_xxxyyy_zzzz[k] = -g_y_0_xxyyy_zzzz[k] * ab_x + g_y_0_xxyyy_xzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyyz_xxxx = cbuffer.data(ig_geom_10_off + 525 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxy = cbuffer.data(ig_geom_10_off + 526 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxxz = cbuffer.data(ig_geom_10_off + 527 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxyy = cbuffer.data(ig_geom_10_off + 528 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxyz = cbuffer.data(ig_geom_10_off + 529 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xxzz = cbuffer.data(ig_geom_10_off + 530 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyyy = cbuffer.data(ig_geom_10_off + 531 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyyz = cbuffer.data(ig_geom_10_off + 532 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xyzz = cbuffer.data(ig_geom_10_off + 533 * ccomps * dcomps);

            auto g_y_0_xxxyyz_xzzz = cbuffer.data(ig_geom_10_off + 534 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyyy = cbuffer.data(ig_geom_10_off + 535 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyyz = cbuffer.data(ig_geom_10_off + 536 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yyzz = cbuffer.data(ig_geom_10_off + 537 * ccomps * dcomps);

            auto g_y_0_xxxyyz_yzzz = cbuffer.data(ig_geom_10_off + 538 * ccomps * dcomps);

            auto g_y_0_xxxyyz_zzzz = cbuffer.data(ig_geom_10_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyyz_xxxx, g_y_0_xxxyyz_xxxy, g_y_0_xxxyyz_xxxz, g_y_0_xxxyyz_xxyy, g_y_0_xxxyyz_xxyz, g_y_0_xxxyyz_xxzz, g_y_0_xxxyyz_xyyy, g_y_0_xxxyyz_xyyz, g_y_0_xxxyyz_xyzz, g_y_0_xxxyyz_xzzz, g_y_0_xxxyyz_yyyy, g_y_0_xxxyyz_yyyz, g_y_0_xxxyyz_yyzz, g_y_0_xxxyyz_yzzz, g_y_0_xxxyyz_zzzz, g_y_0_xxyyz_xxxx, g_y_0_xxyyz_xxxxx, g_y_0_xxyyz_xxxxy, g_y_0_xxyyz_xxxxz, g_y_0_xxyyz_xxxy, g_y_0_xxyyz_xxxyy, g_y_0_xxyyz_xxxyz, g_y_0_xxyyz_xxxz, g_y_0_xxyyz_xxxzz, g_y_0_xxyyz_xxyy, g_y_0_xxyyz_xxyyy, g_y_0_xxyyz_xxyyz, g_y_0_xxyyz_xxyz, g_y_0_xxyyz_xxyzz, g_y_0_xxyyz_xxzz, g_y_0_xxyyz_xxzzz, g_y_0_xxyyz_xyyy, g_y_0_xxyyz_xyyyy, g_y_0_xxyyz_xyyyz, g_y_0_xxyyz_xyyz, g_y_0_xxyyz_xyyzz, g_y_0_xxyyz_xyzz, g_y_0_xxyyz_xyzzz, g_y_0_xxyyz_xzzz, g_y_0_xxyyz_xzzzz, g_y_0_xxyyz_yyyy, g_y_0_xxyyz_yyyz, g_y_0_xxyyz_yyzz, g_y_0_xxyyz_yzzz, g_y_0_xxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyyz_xxxx[k] = -g_y_0_xxyyz_xxxx[k] * ab_x + g_y_0_xxyyz_xxxxx[k];

                g_y_0_xxxyyz_xxxy[k] = -g_y_0_xxyyz_xxxy[k] * ab_x + g_y_0_xxyyz_xxxxy[k];

                g_y_0_xxxyyz_xxxz[k] = -g_y_0_xxyyz_xxxz[k] * ab_x + g_y_0_xxyyz_xxxxz[k];

                g_y_0_xxxyyz_xxyy[k] = -g_y_0_xxyyz_xxyy[k] * ab_x + g_y_0_xxyyz_xxxyy[k];

                g_y_0_xxxyyz_xxyz[k] = -g_y_0_xxyyz_xxyz[k] * ab_x + g_y_0_xxyyz_xxxyz[k];

                g_y_0_xxxyyz_xxzz[k] = -g_y_0_xxyyz_xxzz[k] * ab_x + g_y_0_xxyyz_xxxzz[k];

                g_y_0_xxxyyz_xyyy[k] = -g_y_0_xxyyz_xyyy[k] * ab_x + g_y_0_xxyyz_xxyyy[k];

                g_y_0_xxxyyz_xyyz[k] = -g_y_0_xxyyz_xyyz[k] * ab_x + g_y_0_xxyyz_xxyyz[k];

                g_y_0_xxxyyz_xyzz[k] = -g_y_0_xxyyz_xyzz[k] * ab_x + g_y_0_xxyyz_xxyzz[k];

                g_y_0_xxxyyz_xzzz[k] = -g_y_0_xxyyz_xzzz[k] * ab_x + g_y_0_xxyyz_xxzzz[k];

                g_y_0_xxxyyz_yyyy[k] = -g_y_0_xxyyz_yyyy[k] * ab_x + g_y_0_xxyyz_xyyyy[k];

                g_y_0_xxxyyz_yyyz[k] = -g_y_0_xxyyz_yyyz[k] * ab_x + g_y_0_xxyyz_xyyyz[k];

                g_y_0_xxxyyz_yyzz[k] = -g_y_0_xxyyz_yyzz[k] * ab_x + g_y_0_xxyyz_xyyzz[k];

                g_y_0_xxxyyz_yzzz[k] = -g_y_0_xxyyz_yzzz[k] * ab_x + g_y_0_xxyyz_xyzzz[k];

                g_y_0_xxxyyz_zzzz[k] = -g_y_0_xxyyz_zzzz[k] * ab_x + g_y_0_xxyyz_xzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxyzz_xxxx = cbuffer.data(ig_geom_10_off + 540 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxy = cbuffer.data(ig_geom_10_off + 541 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxxz = cbuffer.data(ig_geom_10_off + 542 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxyy = cbuffer.data(ig_geom_10_off + 543 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxyz = cbuffer.data(ig_geom_10_off + 544 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xxzz = cbuffer.data(ig_geom_10_off + 545 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyyy = cbuffer.data(ig_geom_10_off + 546 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyyz = cbuffer.data(ig_geom_10_off + 547 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xyzz = cbuffer.data(ig_geom_10_off + 548 * ccomps * dcomps);

            auto g_y_0_xxxyzz_xzzz = cbuffer.data(ig_geom_10_off + 549 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyyy = cbuffer.data(ig_geom_10_off + 550 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyyz = cbuffer.data(ig_geom_10_off + 551 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yyzz = cbuffer.data(ig_geom_10_off + 552 * ccomps * dcomps);

            auto g_y_0_xxxyzz_yzzz = cbuffer.data(ig_geom_10_off + 553 * ccomps * dcomps);

            auto g_y_0_xxxyzz_zzzz = cbuffer.data(ig_geom_10_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxyzz_xxxx, g_y_0_xxxyzz_xxxy, g_y_0_xxxyzz_xxxz, g_y_0_xxxyzz_xxyy, g_y_0_xxxyzz_xxyz, g_y_0_xxxyzz_xxzz, g_y_0_xxxyzz_xyyy, g_y_0_xxxyzz_xyyz, g_y_0_xxxyzz_xyzz, g_y_0_xxxyzz_xzzz, g_y_0_xxxyzz_yyyy, g_y_0_xxxyzz_yyyz, g_y_0_xxxyzz_yyzz, g_y_0_xxxyzz_yzzz, g_y_0_xxxyzz_zzzz, g_y_0_xxyzz_xxxx, g_y_0_xxyzz_xxxxx, g_y_0_xxyzz_xxxxy, g_y_0_xxyzz_xxxxz, g_y_0_xxyzz_xxxy, g_y_0_xxyzz_xxxyy, g_y_0_xxyzz_xxxyz, g_y_0_xxyzz_xxxz, g_y_0_xxyzz_xxxzz, g_y_0_xxyzz_xxyy, g_y_0_xxyzz_xxyyy, g_y_0_xxyzz_xxyyz, g_y_0_xxyzz_xxyz, g_y_0_xxyzz_xxyzz, g_y_0_xxyzz_xxzz, g_y_0_xxyzz_xxzzz, g_y_0_xxyzz_xyyy, g_y_0_xxyzz_xyyyy, g_y_0_xxyzz_xyyyz, g_y_0_xxyzz_xyyz, g_y_0_xxyzz_xyyzz, g_y_0_xxyzz_xyzz, g_y_0_xxyzz_xyzzz, g_y_0_xxyzz_xzzz, g_y_0_xxyzz_xzzzz, g_y_0_xxyzz_yyyy, g_y_0_xxyzz_yyyz, g_y_0_xxyzz_yyzz, g_y_0_xxyzz_yzzz, g_y_0_xxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxyzz_xxxx[k] = -g_y_0_xxyzz_xxxx[k] * ab_x + g_y_0_xxyzz_xxxxx[k];

                g_y_0_xxxyzz_xxxy[k] = -g_y_0_xxyzz_xxxy[k] * ab_x + g_y_0_xxyzz_xxxxy[k];

                g_y_0_xxxyzz_xxxz[k] = -g_y_0_xxyzz_xxxz[k] * ab_x + g_y_0_xxyzz_xxxxz[k];

                g_y_0_xxxyzz_xxyy[k] = -g_y_0_xxyzz_xxyy[k] * ab_x + g_y_0_xxyzz_xxxyy[k];

                g_y_0_xxxyzz_xxyz[k] = -g_y_0_xxyzz_xxyz[k] * ab_x + g_y_0_xxyzz_xxxyz[k];

                g_y_0_xxxyzz_xxzz[k] = -g_y_0_xxyzz_xxzz[k] * ab_x + g_y_0_xxyzz_xxxzz[k];

                g_y_0_xxxyzz_xyyy[k] = -g_y_0_xxyzz_xyyy[k] * ab_x + g_y_0_xxyzz_xxyyy[k];

                g_y_0_xxxyzz_xyyz[k] = -g_y_0_xxyzz_xyyz[k] * ab_x + g_y_0_xxyzz_xxyyz[k];

                g_y_0_xxxyzz_xyzz[k] = -g_y_0_xxyzz_xyzz[k] * ab_x + g_y_0_xxyzz_xxyzz[k];

                g_y_0_xxxyzz_xzzz[k] = -g_y_0_xxyzz_xzzz[k] * ab_x + g_y_0_xxyzz_xxzzz[k];

                g_y_0_xxxyzz_yyyy[k] = -g_y_0_xxyzz_yyyy[k] * ab_x + g_y_0_xxyzz_xyyyy[k];

                g_y_0_xxxyzz_yyyz[k] = -g_y_0_xxyzz_yyyz[k] * ab_x + g_y_0_xxyzz_xyyyz[k];

                g_y_0_xxxyzz_yyzz[k] = -g_y_0_xxyzz_yyzz[k] * ab_x + g_y_0_xxyzz_xyyzz[k];

                g_y_0_xxxyzz_yzzz[k] = -g_y_0_xxyzz_yzzz[k] * ab_x + g_y_0_xxyzz_xyzzz[k];

                g_y_0_xxxyzz_zzzz[k] = -g_y_0_xxyzz_zzzz[k] * ab_x + g_y_0_xxyzz_xzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxxzzz_xxxx = cbuffer.data(ig_geom_10_off + 555 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxy = cbuffer.data(ig_geom_10_off + 556 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxxz = cbuffer.data(ig_geom_10_off + 557 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxyy = cbuffer.data(ig_geom_10_off + 558 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxyz = cbuffer.data(ig_geom_10_off + 559 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xxzz = cbuffer.data(ig_geom_10_off + 560 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyyy = cbuffer.data(ig_geom_10_off + 561 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyyz = cbuffer.data(ig_geom_10_off + 562 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xyzz = cbuffer.data(ig_geom_10_off + 563 * ccomps * dcomps);

            auto g_y_0_xxxzzz_xzzz = cbuffer.data(ig_geom_10_off + 564 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyyy = cbuffer.data(ig_geom_10_off + 565 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyyz = cbuffer.data(ig_geom_10_off + 566 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yyzz = cbuffer.data(ig_geom_10_off + 567 * ccomps * dcomps);

            auto g_y_0_xxxzzz_yzzz = cbuffer.data(ig_geom_10_off + 568 * ccomps * dcomps);

            auto g_y_0_xxxzzz_zzzz = cbuffer.data(ig_geom_10_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxxzzz_xxxx, g_y_0_xxxzzz_xxxy, g_y_0_xxxzzz_xxxz, g_y_0_xxxzzz_xxyy, g_y_0_xxxzzz_xxyz, g_y_0_xxxzzz_xxzz, g_y_0_xxxzzz_xyyy, g_y_0_xxxzzz_xyyz, g_y_0_xxxzzz_xyzz, g_y_0_xxxzzz_xzzz, g_y_0_xxxzzz_yyyy, g_y_0_xxxzzz_yyyz, g_y_0_xxxzzz_yyzz, g_y_0_xxxzzz_yzzz, g_y_0_xxxzzz_zzzz, g_y_0_xxzzz_xxxx, g_y_0_xxzzz_xxxxx, g_y_0_xxzzz_xxxxy, g_y_0_xxzzz_xxxxz, g_y_0_xxzzz_xxxy, g_y_0_xxzzz_xxxyy, g_y_0_xxzzz_xxxyz, g_y_0_xxzzz_xxxz, g_y_0_xxzzz_xxxzz, g_y_0_xxzzz_xxyy, g_y_0_xxzzz_xxyyy, g_y_0_xxzzz_xxyyz, g_y_0_xxzzz_xxyz, g_y_0_xxzzz_xxyzz, g_y_0_xxzzz_xxzz, g_y_0_xxzzz_xxzzz, g_y_0_xxzzz_xyyy, g_y_0_xxzzz_xyyyy, g_y_0_xxzzz_xyyyz, g_y_0_xxzzz_xyyz, g_y_0_xxzzz_xyyzz, g_y_0_xxzzz_xyzz, g_y_0_xxzzz_xyzzz, g_y_0_xxzzz_xzzz, g_y_0_xxzzz_xzzzz, g_y_0_xxzzz_yyyy, g_y_0_xxzzz_yyyz, g_y_0_xxzzz_yyzz, g_y_0_xxzzz_yzzz, g_y_0_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxxzzz_xxxx[k] = -g_y_0_xxzzz_xxxx[k] * ab_x + g_y_0_xxzzz_xxxxx[k];

                g_y_0_xxxzzz_xxxy[k] = -g_y_0_xxzzz_xxxy[k] * ab_x + g_y_0_xxzzz_xxxxy[k];

                g_y_0_xxxzzz_xxxz[k] = -g_y_0_xxzzz_xxxz[k] * ab_x + g_y_0_xxzzz_xxxxz[k];

                g_y_0_xxxzzz_xxyy[k] = -g_y_0_xxzzz_xxyy[k] * ab_x + g_y_0_xxzzz_xxxyy[k];

                g_y_0_xxxzzz_xxyz[k] = -g_y_0_xxzzz_xxyz[k] * ab_x + g_y_0_xxzzz_xxxyz[k];

                g_y_0_xxxzzz_xxzz[k] = -g_y_0_xxzzz_xxzz[k] * ab_x + g_y_0_xxzzz_xxxzz[k];

                g_y_0_xxxzzz_xyyy[k] = -g_y_0_xxzzz_xyyy[k] * ab_x + g_y_0_xxzzz_xxyyy[k];

                g_y_0_xxxzzz_xyyz[k] = -g_y_0_xxzzz_xyyz[k] * ab_x + g_y_0_xxzzz_xxyyz[k];

                g_y_0_xxxzzz_xyzz[k] = -g_y_0_xxzzz_xyzz[k] * ab_x + g_y_0_xxzzz_xxyzz[k];

                g_y_0_xxxzzz_xzzz[k] = -g_y_0_xxzzz_xzzz[k] * ab_x + g_y_0_xxzzz_xxzzz[k];

                g_y_0_xxxzzz_yyyy[k] = -g_y_0_xxzzz_yyyy[k] * ab_x + g_y_0_xxzzz_xyyyy[k];

                g_y_0_xxxzzz_yyyz[k] = -g_y_0_xxzzz_yyyz[k] * ab_x + g_y_0_xxzzz_xyyyz[k];

                g_y_0_xxxzzz_yyzz[k] = -g_y_0_xxzzz_yyzz[k] * ab_x + g_y_0_xxzzz_xyyzz[k];

                g_y_0_xxxzzz_yzzz[k] = -g_y_0_xxzzz_yzzz[k] * ab_x + g_y_0_xxzzz_xyzzz[k];

                g_y_0_xxxzzz_zzzz[k] = -g_y_0_xxzzz_zzzz[k] * ab_x + g_y_0_xxzzz_xzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyy_xxxx = cbuffer.data(ig_geom_10_off + 570 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxy = cbuffer.data(ig_geom_10_off + 571 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxxz = cbuffer.data(ig_geom_10_off + 572 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxyy = cbuffer.data(ig_geom_10_off + 573 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxyz = cbuffer.data(ig_geom_10_off + 574 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xxzz = cbuffer.data(ig_geom_10_off + 575 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyyy = cbuffer.data(ig_geom_10_off + 576 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyyz = cbuffer.data(ig_geom_10_off + 577 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xyzz = cbuffer.data(ig_geom_10_off + 578 * ccomps * dcomps);

            auto g_y_0_xxyyyy_xzzz = cbuffer.data(ig_geom_10_off + 579 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyyy = cbuffer.data(ig_geom_10_off + 580 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyyz = cbuffer.data(ig_geom_10_off + 581 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yyzz = cbuffer.data(ig_geom_10_off + 582 * ccomps * dcomps);

            auto g_y_0_xxyyyy_yzzz = cbuffer.data(ig_geom_10_off + 583 * ccomps * dcomps);

            auto g_y_0_xxyyyy_zzzz = cbuffer.data(ig_geom_10_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyyy_xxxx, g_y_0_xxyyyy_xxxy, g_y_0_xxyyyy_xxxz, g_y_0_xxyyyy_xxyy, g_y_0_xxyyyy_xxyz, g_y_0_xxyyyy_xxzz, g_y_0_xxyyyy_xyyy, g_y_0_xxyyyy_xyyz, g_y_0_xxyyyy_xyzz, g_y_0_xxyyyy_xzzz, g_y_0_xxyyyy_yyyy, g_y_0_xxyyyy_yyyz, g_y_0_xxyyyy_yyzz, g_y_0_xxyyyy_yzzz, g_y_0_xxyyyy_zzzz, g_y_0_xyyyy_xxxx, g_y_0_xyyyy_xxxxx, g_y_0_xyyyy_xxxxy, g_y_0_xyyyy_xxxxz, g_y_0_xyyyy_xxxy, g_y_0_xyyyy_xxxyy, g_y_0_xyyyy_xxxyz, g_y_0_xyyyy_xxxz, g_y_0_xyyyy_xxxzz, g_y_0_xyyyy_xxyy, g_y_0_xyyyy_xxyyy, g_y_0_xyyyy_xxyyz, g_y_0_xyyyy_xxyz, g_y_0_xyyyy_xxyzz, g_y_0_xyyyy_xxzz, g_y_0_xyyyy_xxzzz, g_y_0_xyyyy_xyyy, g_y_0_xyyyy_xyyyy, g_y_0_xyyyy_xyyyz, g_y_0_xyyyy_xyyz, g_y_0_xyyyy_xyyzz, g_y_0_xyyyy_xyzz, g_y_0_xyyyy_xyzzz, g_y_0_xyyyy_xzzz, g_y_0_xyyyy_xzzzz, g_y_0_xyyyy_yyyy, g_y_0_xyyyy_yyyz, g_y_0_xyyyy_yyzz, g_y_0_xyyyy_yzzz, g_y_0_xyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyy_xxxx[k] = -g_y_0_xyyyy_xxxx[k] * ab_x + g_y_0_xyyyy_xxxxx[k];

                g_y_0_xxyyyy_xxxy[k] = -g_y_0_xyyyy_xxxy[k] * ab_x + g_y_0_xyyyy_xxxxy[k];

                g_y_0_xxyyyy_xxxz[k] = -g_y_0_xyyyy_xxxz[k] * ab_x + g_y_0_xyyyy_xxxxz[k];

                g_y_0_xxyyyy_xxyy[k] = -g_y_0_xyyyy_xxyy[k] * ab_x + g_y_0_xyyyy_xxxyy[k];

                g_y_0_xxyyyy_xxyz[k] = -g_y_0_xyyyy_xxyz[k] * ab_x + g_y_0_xyyyy_xxxyz[k];

                g_y_0_xxyyyy_xxzz[k] = -g_y_0_xyyyy_xxzz[k] * ab_x + g_y_0_xyyyy_xxxzz[k];

                g_y_0_xxyyyy_xyyy[k] = -g_y_0_xyyyy_xyyy[k] * ab_x + g_y_0_xyyyy_xxyyy[k];

                g_y_0_xxyyyy_xyyz[k] = -g_y_0_xyyyy_xyyz[k] * ab_x + g_y_0_xyyyy_xxyyz[k];

                g_y_0_xxyyyy_xyzz[k] = -g_y_0_xyyyy_xyzz[k] * ab_x + g_y_0_xyyyy_xxyzz[k];

                g_y_0_xxyyyy_xzzz[k] = -g_y_0_xyyyy_xzzz[k] * ab_x + g_y_0_xyyyy_xxzzz[k];

                g_y_0_xxyyyy_yyyy[k] = -g_y_0_xyyyy_yyyy[k] * ab_x + g_y_0_xyyyy_xyyyy[k];

                g_y_0_xxyyyy_yyyz[k] = -g_y_0_xyyyy_yyyz[k] * ab_x + g_y_0_xyyyy_xyyyz[k];

                g_y_0_xxyyyy_yyzz[k] = -g_y_0_xyyyy_yyzz[k] * ab_x + g_y_0_xyyyy_xyyzz[k];

                g_y_0_xxyyyy_yzzz[k] = -g_y_0_xyyyy_yzzz[k] * ab_x + g_y_0_xyyyy_xyzzz[k];

                g_y_0_xxyyyy_zzzz[k] = -g_y_0_xyyyy_zzzz[k] * ab_x + g_y_0_xyyyy_xzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyyz_xxxx = cbuffer.data(ig_geom_10_off + 585 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxy = cbuffer.data(ig_geom_10_off + 586 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxxz = cbuffer.data(ig_geom_10_off + 587 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxyy = cbuffer.data(ig_geom_10_off + 588 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxyz = cbuffer.data(ig_geom_10_off + 589 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xxzz = cbuffer.data(ig_geom_10_off + 590 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyyy = cbuffer.data(ig_geom_10_off + 591 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyyz = cbuffer.data(ig_geom_10_off + 592 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xyzz = cbuffer.data(ig_geom_10_off + 593 * ccomps * dcomps);

            auto g_y_0_xxyyyz_xzzz = cbuffer.data(ig_geom_10_off + 594 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyyy = cbuffer.data(ig_geom_10_off + 595 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyyz = cbuffer.data(ig_geom_10_off + 596 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yyzz = cbuffer.data(ig_geom_10_off + 597 * ccomps * dcomps);

            auto g_y_0_xxyyyz_yzzz = cbuffer.data(ig_geom_10_off + 598 * ccomps * dcomps);

            auto g_y_0_xxyyyz_zzzz = cbuffer.data(ig_geom_10_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyyz_xxxx, g_y_0_xxyyyz_xxxy, g_y_0_xxyyyz_xxxz, g_y_0_xxyyyz_xxyy, g_y_0_xxyyyz_xxyz, g_y_0_xxyyyz_xxzz, g_y_0_xxyyyz_xyyy, g_y_0_xxyyyz_xyyz, g_y_0_xxyyyz_xyzz, g_y_0_xxyyyz_xzzz, g_y_0_xxyyyz_yyyy, g_y_0_xxyyyz_yyyz, g_y_0_xxyyyz_yyzz, g_y_0_xxyyyz_yzzz, g_y_0_xxyyyz_zzzz, g_y_0_xyyyz_xxxx, g_y_0_xyyyz_xxxxx, g_y_0_xyyyz_xxxxy, g_y_0_xyyyz_xxxxz, g_y_0_xyyyz_xxxy, g_y_0_xyyyz_xxxyy, g_y_0_xyyyz_xxxyz, g_y_0_xyyyz_xxxz, g_y_0_xyyyz_xxxzz, g_y_0_xyyyz_xxyy, g_y_0_xyyyz_xxyyy, g_y_0_xyyyz_xxyyz, g_y_0_xyyyz_xxyz, g_y_0_xyyyz_xxyzz, g_y_0_xyyyz_xxzz, g_y_0_xyyyz_xxzzz, g_y_0_xyyyz_xyyy, g_y_0_xyyyz_xyyyy, g_y_0_xyyyz_xyyyz, g_y_0_xyyyz_xyyz, g_y_0_xyyyz_xyyzz, g_y_0_xyyyz_xyzz, g_y_0_xyyyz_xyzzz, g_y_0_xyyyz_xzzz, g_y_0_xyyyz_xzzzz, g_y_0_xyyyz_yyyy, g_y_0_xyyyz_yyyz, g_y_0_xyyyz_yyzz, g_y_0_xyyyz_yzzz, g_y_0_xyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyyz_xxxx[k] = -g_y_0_xyyyz_xxxx[k] * ab_x + g_y_0_xyyyz_xxxxx[k];

                g_y_0_xxyyyz_xxxy[k] = -g_y_0_xyyyz_xxxy[k] * ab_x + g_y_0_xyyyz_xxxxy[k];

                g_y_0_xxyyyz_xxxz[k] = -g_y_0_xyyyz_xxxz[k] * ab_x + g_y_0_xyyyz_xxxxz[k];

                g_y_0_xxyyyz_xxyy[k] = -g_y_0_xyyyz_xxyy[k] * ab_x + g_y_0_xyyyz_xxxyy[k];

                g_y_0_xxyyyz_xxyz[k] = -g_y_0_xyyyz_xxyz[k] * ab_x + g_y_0_xyyyz_xxxyz[k];

                g_y_0_xxyyyz_xxzz[k] = -g_y_0_xyyyz_xxzz[k] * ab_x + g_y_0_xyyyz_xxxzz[k];

                g_y_0_xxyyyz_xyyy[k] = -g_y_0_xyyyz_xyyy[k] * ab_x + g_y_0_xyyyz_xxyyy[k];

                g_y_0_xxyyyz_xyyz[k] = -g_y_0_xyyyz_xyyz[k] * ab_x + g_y_0_xyyyz_xxyyz[k];

                g_y_0_xxyyyz_xyzz[k] = -g_y_0_xyyyz_xyzz[k] * ab_x + g_y_0_xyyyz_xxyzz[k];

                g_y_0_xxyyyz_xzzz[k] = -g_y_0_xyyyz_xzzz[k] * ab_x + g_y_0_xyyyz_xxzzz[k];

                g_y_0_xxyyyz_yyyy[k] = -g_y_0_xyyyz_yyyy[k] * ab_x + g_y_0_xyyyz_xyyyy[k];

                g_y_0_xxyyyz_yyyz[k] = -g_y_0_xyyyz_yyyz[k] * ab_x + g_y_0_xyyyz_xyyyz[k];

                g_y_0_xxyyyz_yyzz[k] = -g_y_0_xyyyz_yyzz[k] * ab_x + g_y_0_xyyyz_xyyzz[k];

                g_y_0_xxyyyz_yzzz[k] = -g_y_0_xyyyz_yzzz[k] * ab_x + g_y_0_xyyyz_xyzzz[k];

                g_y_0_xxyyyz_zzzz[k] = -g_y_0_xyyyz_zzzz[k] * ab_x + g_y_0_xyyyz_xzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyyzz_xxxx = cbuffer.data(ig_geom_10_off + 600 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxy = cbuffer.data(ig_geom_10_off + 601 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxxz = cbuffer.data(ig_geom_10_off + 602 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxyy = cbuffer.data(ig_geom_10_off + 603 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxyz = cbuffer.data(ig_geom_10_off + 604 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xxzz = cbuffer.data(ig_geom_10_off + 605 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyyy = cbuffer.data(ig_geom_10_off + 606 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyyz = cbuffer.data(ig_geom_10_off + 607 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xyzz = cbuffer.data(ig_geom_10_off + 608 * ccomps * dcomps);

            auto g_y_0_xxyyzz_xzzz = cbuffer.data(ig_geom_10_off + 609 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyyy = cbuffer.data(ig_geom_10_off + 610 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyyz = cbuffer.data(ig_geom_10_off + 611 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yyzz = cbuffer.data(ig_geom_10_off + 612 * ccomps * dcomps);

            auto g_y_0_xxyyzz_yzzz = cbuffer.data(ig_geom_10_off + 613 * ccomps * dcomps);

            auto g_y_0_xxyyzz_zzzz = cbuffer.data(ig_geom_10_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyyzz_xxxx, g_y_0_xxyyzz_xxxy, g_y_0_xxyyzz_xxxz, g_y_0_xxyyzz_xxyy, g_y_0_xxyyzz_xxyz, g_y_0_xxyyzz_xxzz, g_y_0_xxyyzz_xyyy, g_y_0_xxyyzz_xyyz, g_y_0_xxyyzz_xyzz, g_y_0_xxyyzz_xzzz, g_y_0_xxyyzz_yyyy, g_y_0_xxyyzz_yyyz, g_y_0_xxyyzz_yyzz, g_y_0_xxyyzz_yzzz, g_y_0_xxyyzz_zzzz, g_y_0_xyyzz_xxxx, g_y_0_xyyzz_xxxxx, g_y_0_xyyzz_xxxxy, g_y_0_xyyzz_xxxxz, g_y_0_xyyzz_xxxy, g_y_0_xyyzz_xxxyy, g_y_0_xyyzz_xxxyz, g_y_0_xyyzz_xxxz, g_y_0_xyyzz_xxxzz, g_y_0_xyyzz_xxyy, g_y_0_xyyzz_xxyyy, g_y_0_xyyzz_xxyyz, g_y_0_xyyzz_xxyz, g_y_0_xyyzz_xxyzz, g_y_0_xyyzz_xxzz, g_y_0_xyyzz_xxzzz, g_y_0_xyyzz_xyyy, g_y_0_xyyzz_xyyyy, g_y_0_xyyzz_xyyyz, g_y_0_xyyzz_xyyz, g_y_0_xyyzz_xyyzz, g_y_0_xyyzz_xyzz, g_y_0_xyyzz_xyzzz, g_y_0_xyyzz_xzzz, g_y_0_xyyzz_xzzzz, g_y_0_xyyzz_yyyy, g_y_0_xyyzz_yyyz, g_y_0_xyyzz_yyzz, g_y_0_xyyzz_yzzz, g_y_0_xyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyyzz_xxxx[k] = -g_y_0_xyyzz_xxxx[k] * ab_x + g_y_0_xyyzz_xxxxx[k];

                g_y_0_xxyyzz_xxxy[k] = -g_y_0_xyyzz_xxxy[k] * ab_x + g_y_0_xyyzz_xxxxy[k];

                g_y_0_xxyyzz_xxxz[k] = -g_y_0_xyyzz_xxxz[k] * ab_x + g_y_0_xyyzz_xxxxz[k];

                g_y_0_xxyyzz_xxyy[k] = -g_y_0_xyyzz_xxyy[k] * ab_x + g_y_0_xyyzz_xxxyy[k];

                g_y_0_xxyyzz_xxyz[k] = -g_y_0_xyyzz_xxyz[k] * ab_x + g_y_0_xyyzz_xxxyz[k];

                g_y_0_xxyyzz_xxzz[k] = -g_y_0_xyyzz_xxzz[k] * ab_x + g_y_0_xyyzz_xxxzz[k];

                g_y_0_xxyyzz_xyyy[k] = -g_y_0_xyyzz_xyyy[k] * ab_x + g_y_0_xyyzz_xxyyy[k];

                g_y_0_xxyyzz_xyyz[k] = -g_y_0_xyyzz_xyyz[k] * ab_x + g_y_0_xyyzz_xxyyz[k];

                g_y_0_xxyyzz_xyzz[k] = -g_y_0_xyyzz_xyzz[k] * ab_x + g_y_0_xyyzz_xxyzz[k];

                g_y_0_xxyyzz_xzzz[k] = -g_y_0_xyyzz_xzzz[k] * ab_x + g_y_0_xyyzz_xxzzz[k];

                g_y_0_xxyyzz_yyyy[k] = -g_y_0_xyyzz_yyyy[k] * ab_x + g_y_0_xyyzz_xyyyy[k];

                g_y_0_xxyyzz_yyyz[k] = -g_y_0_xyyzz_yyyz[k] * ab_x + g_y_0_xyyzz_xyyyz[k];

                g_y_0_xxyyzz_yyzz[k] = -g_y_0_xyyzz_yyzz[k] * ab_x + g_y_0_xyyzz_xyyzz[k];

                g_y_0_xxyyzz_yzzz[k] = -g_y_0_xyyzz_yzzz[k] * ab_x + g_y_0_xyyzz_xyzzz[k];

                g_y_0_xxyyzz_zzzz[k] = -g_y_0_xyyzz_zzzz[k] * ab_x + g_y_0_xyyzz_xzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxyzzz_xxxx = cbuffer.data(ig_geom_10_off + 615 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxy = cbuffer.data(ig_geom_10_off + 616 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxxz = cbuffer.data(ig_geom_10_off + 617 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxyy = cbuffer.data(ig_geom_10_off + 618 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxyz = cbuffer.data(ig_geom_10_off + 619 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xxzz = cbuffer.data(ig_geom_10_off + 620 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyyy = cbuffer.data(ig_geom_10_off + 621 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyyz = cbuffer.data(ig_geom_10_off + 622 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xyzz = cbuffer.data(ig_geom_10_off + 623 * ccomps * dcomps);

            auto g_y_0_xxyzzz_xzzz = cbuffer.data(ig_geom_10_off + 624 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyyy = cbuffer.data(ig_geom_10_off + 625 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyyz = cbuffer.data(ig_geom_10_off + 626 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yyzz = cbuffer.data(ig_geom_10_off + 627 * ccomps * dcomps);

            auto g_y_0_xxyzzz_yzzz = cbuffer.data(ig_geom_10_off + 628 * ccomps * dcomps);

            auto g_y_0_xxyzzz_zzzz = cbuffer.data(ig_geom_10_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxyzzz_xxxx, g_y_0_xxyzzz_xxxy, g_y_0_xxyzzz_xxxz, g_y_0_xxyzzz_xxyy, g_y_0_xxyzzz_xxyz, g_y_0_xxyzzz_xxzz, g_y_0_xxyzzz_xyyy, g_y_0_xxyzzz_xyyz, g_y_0_xxyzzz_xyzz, g_y_0_xxyzzz_xzzz, g_y_0_xxyzzz_yyyy, g_y_0_xxyzzz_yyyz, g_y_0_xxyzzz_yyzz, g_y_0_xxyzzz_yzzz, g_y_0_xxyzzz_zzzz, g_y_0_xyzzz_xxxx, g_y_0_xyzzz_xxxxx, g_y_0_xyzzz_xxxxy, g_y_0_xyzzz_xxxxz, g_y_0_xyzzz_xxxy, g_y_0_xyzzz_xxxyy, g_y_0_xyzzz_xxxyz, g_y_0_xyzzz_xxxz, g_y_0_xyzzz_xxxzz, g_y_0_xyzzz_xxyy, g_y_0_xyzzz_xxyyy, g_y_0_xyzzz_xxyyz, g_y_0_xyzzz_xxyz, g_y_0_xyzzz_xxyzz, g_y_0_xyzzz_xxzz, g_y_0_xyzzz_xxzzz, g_y_0_xyzzz_xyyy, g_y_0_xyzzz_xyyyy, g_y_0_xyzzz_xyyyz, g_y_0_xyzzz_xyyz, g_y_0_xyzzz_xyyzz, g_y_0_xyzzz_xyzz, g_y_0_xyzzz_xyzzz, g_y_0_xyzzz_xzzz, g_y_0_xyzzz_xzzzz, g_y_0_xyzzz_yyyy, g_y_0_xyzzz_yyyz, g_y_0_xyzzz_yyzz, g_y_0_xyzzz_yzzz, g_y_0_xyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxyzzz_xxxx[k] = -g_y_0_xyzzz_xxxx[k] * ab_x + g_y_0_xyzzz_xxxxx[k];

                g_y_0_xxyzzz_xxxy[k] = -g_y_0_xyzzz_xxxy[k] * ab_x + g_y_0_xyzzz_xxxxy[k];

                g_y_0_xxyzzz_xxxz[k] = -g_y_0_xyzzz_xxxz[k] * ab_x + g_y_0_xyzzz_xxxxz[k];

                g_y_0_xxyzzz_xxyy[k] = -g_y_0_xyzzz_xxyy[k] * ab_x + g_y_0_xyzzz_xxxyy[k];

                g_y_0_xxyzzz_xxyz[k] = -g_y_0_xyzzz_xxyz[k] * ab_x + g_y_0_xyzzz_xxxyz[k];

                g_y_0_xxyzzz_xxzz[k] = -g_y_0_xyzzz_xxzz[k] * ab_x + g_y_0_xyzzz_xxxzz[k];

                g_y_0_xxyzzz_xyyy[k] = -g_y_0_xyzzz_xyyy[k] * ab_x + g_y_0_xyzzz_xxyyy[k];

                g_y_0_xxyzzz_xyyz[k] = -g_y_0_xyzzz_xyyz[k] * ab_x + g_y_0_xyzzz_xxyyz[k];

                g_y_0_xxyzzz_xyzz[k] = -g_y_0_xyzzz_xyzz[k] * ab_x + g_y_0_xyzzz_xxyzz[k];

                g_y_0_xxyzzz_xzzz[k] = -g_y_0_xyzzz_xzzz[k] * ab_x + g_y_0_xyzzz_xxzzz[k];

                g_y_0_xxyzzz_yyyy[k] = -g_y_0_xyzzz_yyyy[k] * ab_x + g_y_0_xyzzz_xyyyy[k];

                g_y_0_xxyzzz_yyyz[k] = -g_y_0_xyzzz_yyyz[k] * ab_x + g_y_0_xyzzz_xyyyz[k];

                g_y_0_xxyzzz_yyzz[k] = -g_y_0_xyzzz_yyzz[k] * ab_x + g_y_0_xyzzz_xyyzz[k];

                g_y_0_xxyzzz_yzzz[k] = -g_y_0_xyzzz_yzzz[k] * ab_x + g_y_0_xyzzz_xyzzz[k];

                g_y_0_xxyzzz_zzzz[k] = -g_y_0_xyzzz_zzzz[k] * ab_x + g_y_0_xyzzz_xzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

            auto g_y_0_xxzzzz_xxxx = cbuffer.data(ig_geom_10_off + 630 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxy = cbuffer.data(ig_geom_10_off + 631 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxxz = cbuffer.data(ig_geom_10_off + 632 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxyy = cbuffer.data(ig_geom_10_off + 633 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxyz = cbuffer.data(ig_geom_10_off + 634 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xxzz = cbuffer.data(ig_geom_10_off + 635 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyyy = cbuffer.data(ig_geom_10_off + 636 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyyz = cbuffer.data(ig_geom_10_off + 637 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xyzz = cbuffer.data(ig_geom_10_off + 638 * ccomps * dcomps);

            auto g_y_0_xxzzzz_xzzz = cbuffer.data(ig_geom_10_off + 639 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyyy = cbuffer.data(ig_geom_10_off + 640 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyyz = cbuffer.data(ig_geom_10_off + 641 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yyzz = cbuffer.data(ig_geom_10_off + 642 * ccomps * dcomps);

            auto g_y_0_xxzzzz_yzzz = cbuffer.data(ig_geom_10_off + 643 * ccomps * dcomps);

            auto g_y_0_xxzzzz_zzzz = cbuffer.data(ig_geom_10_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xxzzzz_xxxx, g_y_0_xxzzzz_xxxy, g_y_0_xxzzzz_xxxz, g_y_0_xxzzzz_xxyy, g_y_0_xxzzzz_xxyz, g_y_0_xxzzzz_xxzz, g_y_0_xxzzzz_xyyy, g_y_0_xxzzzz_xyyz, g_y_0_xxzzzz_xyzz, g_y_0_xxzzzz_xzzz, g_y_0_xxzzzz_yyyy, g_y_0_xxzzzz_yyyz, g_y_0_xxzzzz_yyzz, g_y_0_xxzzzz_yzzz, g_y_0_xxzzzz_zzzz, g_y_0_xzzzz_xxxx, g_y_0_xzzzz_xxxxx, g_y_0_xzzzz_xxxxy, g_y_0_xzzzz_xxxxz, g_y_0_xzzzz_xxxy, g_y_0_xzzzz_xxxyy, g_y_0_xzzzz_xxxyz, g_y_0_xzzzz_xxxz, g_y_0_xzzzz_xxxzz, g_y_0_xzzzz_xxyy, g_y_0_xzzzz_xxyyy, g_y_0_xzzzz_xxyyz, g_y_0_xzzzz_xxyz, g_y_0_xzzzz_xxyzz, g_y_0_xzzzz_xxzz, g_y_0_xzzzz_xxzzz, g_y_0_xzzzz_xyyy, g_y_0_xzzzz_xyyyy, g_y_0_xzzzz_xyyyz, g_y_0_xzzzz_xyyz, g_y_0_xzzzz_xyyzz, g_y_0_xzzzz_xyzz, g_y_0_xzzzz_xyzzz, g_y_0_xzzzz_xzzz, g_y_0_xzzzz_xzzzz, g_y_0_xzzzz_yyyy, g_y_0_xzzzz_yyyz, g_y_0_xzzzz_yyzz, g_y_0_xzzzz_yzzz, g_y_0_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xxzzzz_xxxx[k] = -g_y_0_xzzzz_xxxx[k] * ab_x + g_y_0_xzzzz_xxxxx[k];

                g_y_0_xxzzzz_xxxy[k] = -g_y_0_xzzzz_xxxy[k] * ab_x + g_y_0_xzzzz_xxxxy[k];

                g_y_0_xxzzzz_xxxz[k] = -g_y_0_xzzzz_xxxz[k] * ab_x + g_y_0_xzzzz_xxxxz[k];

                g_y_0_xxzzzz_xxyy[k] = -g_y_0_xzzzz_xxyy[k] * ab_x + g_y_0_xzzzz_xxxyy[k];

                g_y_0_xxzzzz_xxyz[k] = -g_y_0_xzzzz_xxyz[k] * ab_x + g_y_0_xzzzz_xxxyz[k];

                g_y_0_xxzzzz_xxzz[k] = -g_y_0_xzzzz_xxzz[k] * ab_x + g_y_0_xzzzz_xxxzz[k];

                g_y_0_xxzzzz_xyyy[k] = -g_y_0_xzzzz_xyyy[k] * ab_x + g_y_0_xzzzz_xxyyy[k];

                g_y_0_xxzzzz_xyyz[k] = -g_y_0_xzzzz_xyyz[k] * ab_x + g_y_0_xzzzz_xxyyz[k];

                g_y_0_xxzzzz_xyzz[k] = -g_y_0_xzzzz_xyzz[k] * ab_x + g_y_0_xzzzz_xxyzz[k];

                g_y_0_xxzzzz_xzzz[k] = -g_y_0_xzzzz_xzzz[k] * ab_x + g_y_0_xzzzz_xxzzz[k];

                g_y_0_xxzzzz_yyyy[k] = -g_y_0_xzzzz_yyyy[k] * ab_x + g_y_0_xzzzz_xyyyy[k];

                g_y_0_xxzzzz_yyyz[k] = -g_y_0_xzzzz_yyyz[k] * ab_x + g_y_0_xzzzz_xyyyz[k];

                g_y_0_xxzzzz_yyzz[k] = -g_y_0_xzzzz_yyzz[k] * ab_x + g_y_0_xzzzz_xyyzz[k];

                g_y_0_xxzzzz_yzzz[k] = -g_y_0_xzzzz_yzzz[k] * ab_x + g_y_0_xzzzz_xyzzz[k];

                g_y_0_xxzzzz_zzzz[k] = -g_y_0_xzzzz_zzzz[k] * ab_x + g_y_0_xzzzz_xzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 645 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 646 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 647 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 648 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 649 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 650 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 651 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 652 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 653 * ccomps * dcomps);

            auto g_y_0_xyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 654 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 655 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 656 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 657 * ccomps * dcomps);

            auto g_y_0_xyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 658 * ccomps * dcomps);

            auto g_y_0_xyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyyy_xxxx, g_y_0_xyyyyy_xxxy, g_y_0_xyyyyy_xxxz, g_y_0_xyyyyy_xxyy, g_y_0_xyyyyy_xxyz, g_y_0_xyyyyy_xxzz, g_y_0_xyyyyy_xyyy, g_y_0_xyyyyy_xyyz, g_y_0_xyyyyy_xyzz, g_y_0_xyyyyy_xzzz, g_y_0_xyyyyy_yyyy, g_y_0_xyyyyy_yyyz, g_y_0_xyyyyy_yyzz, g_y_0_xyyyyy_yzzz, g_y_0_xyyyyy_zzzz, g_y_0_yyyyy_xxxx, g_y_0_yyyyy_xxxxx, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_yyyy, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyy_xxxx[k] = -g_y_0_yyyyy_xxxx[k] * ab_x + g_y_0_yyyyy_xxxxx[k];

                g_y_0_xyyyyy_xxxy[k] = -g_y_0_yyyyy_xxxy[k] * ab_x + g_y_0_yyyyy_xxxxy[k];

                g_y_0_xyyyyy_xxxz[k] = -g_y_0_yyyyy_xxxz[k] * ab_x + g_y_0_yyyyy_xxxxz[k];

                g_y_0_xyyyyy_xxyy[k] = -g_y_0_yyyyy_xxyy[k] * ab_x + g_y_0_yyyyy_xxxyy[k];

                g_y_0_xyyyyy_xxyz[k] = -g_y_0_yyyyy_xxyz[k] * ab_x + g_y_0_yyyyy_xxxyz[k];

                g_y_0_xyyyyy_xxzz[k] = -g_y_0_yyyyy_xxzz[k] * ab_x + g_y_0_yyyyy_xxxzz[k];

                g_y_0_xyyyyy_xyyy[k] = -g_y_0_yyyyy_xyyy[k] * ab_x + g_y_0_yyyyy_xxyyy[k];

                g_y_0_xyyyyy_xyyz[k] = -g_y_0_yyyyy_xyyz[k] * ab_x + g_y_0_yyyyy_xxyyz[k];

                g_y_0_xyyyyy_xyzz[k] = -g_y_0_yyyyy_xyzz[k] * ab_x + g_y_0_yyyyy_xxyzz[k];

                g_y_0_xyyyyy_xzzz[k] = -g_y_0_yyyyy_xzzz[k] * ab_x + g_y_0_yyyyy_xxzzz[k];

                g_y_0_xyyyyy_yyyy[k] = -g_y_0_yyyyy_yyyy[k] * ab_x + g_y_0_yyyyy_xyyyy[k];

                g_y_0_xyyyyy_yyyz[k] = -g_y_0_yyyyy_yyyz[k] * ab_x + g_y_0_yyyyy_xyyyz[k];

                g_y_0_xyyyyy_yyzz[k] = -g_y_0_yyyyy_yyzz[k] * ab_x + g_y_0_yyyyy_xyyzz[k];

                g_y_0_xyyyyy_yzzz[k] = -g_y_0_yyyyy_yzzz[k] * ab_x + g_y_0_yyyyy_xyzzz[k];

                g_y_0_xyyyyy_zzzz[k] = -g_y_0_yyyyy_zzzz[k] * ab_x + g_y_0_yyyyy_xzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 660 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 661 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 662 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 663 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 664 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 665 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 666 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 667 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 668 * ccomps * dcomps);

            auto g_y_0_xyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 669 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 670 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 671 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 672 * ccomps * dcomps);

            auto g_y_0_xyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 673 * ccomps * dcomps);

            auto g_y_0_xyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyyz_xxxx, g_y_0_xyyyyz_xxxy, g_y_0_xyyyyz_xxxz, g_y_0_xyyyyz_xxyy, g_y_0_xyyyyz_xxyz, g_y_0_xyyyyz_xxzz, g_y_0_xyyyyz_xyyy, g_y_0_xyyyyz_xyyz, g_y_0_xyyyyz_xyzz, g_y_0_xyyyyz_xzzz, g_y_0_xyyyyz_yyyy, g_y_0_xyyyyz_yyyz, g_y_0_xyyyyz_yyzz, g_y_0_xyyyyz_yzzz, g_y_0_xyyyyz_zzzz, g_y_0_yyyyz_xxxx, g_y_0_yyyyz_xxxxx, g_y_0_yyyyz_xxxxy, g_y_0_yyyyz_xxxxz, g_y_0_yyyyz_xxxy, g_y_0_yyyyz_xxxyy, g_y_0_yyyyz_xxxyz, g_y_0_yyyyz_xxxz, g_y_0_yyyyz_xxxzz, g_y_0_yyyyz_xxyy, g_y_0_yyyyz_xxyyy, g_y_0_yyyyz_xxyyz, g_y_0_yyyyz_xxyz, g_y_0_yyyyz_xxyzz, g_y_0_yyyyz_xxzz, g_y_0_yyyyz_xxzzz, g_y_0_yyyyz_xyyy, g_y_0_yyyyz_xyyyy, g_y_0_yyyyz_xyyyz, g_y_0_yyyyz_xyyz, g_y_0_yyyyz_xyyzz, g_y_0_yyyyz_xyzz, g_y_0_yyyyz_xyzzz, g_y_0_yyyyz_xzzz, g_y_0_yyyyz_xzzzz, g_y_0_yyyyz_yyyy, g_y_0_yyyyz_yyyz, g_y_0_yyyyz_yyzz, g_y_0_yyyyz_yzzz, g_y_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyyz_xxxx[k] = -g_y_0_yyyyz_xxxx[k] * ab_x + g_y_0_yyyyz_xxxxx[k];

                g_y_0_xyyyyz_xxxy[k] = -g_y_0_yyyyz_xxxy[k] * ab_x + g_y_0_yyyyz_xxxxy[k];

                g_y_0_xyyyyz_xxxz[k] = -g_y_0_yyyyz_xxxz[k] * ab_x + g_y_0_yyyyz_xxxxz[k];

                g_y_0_xyyyyz_xxyy[k] = -g_y_0_yyyyz_xxyy[k] * ab_x + g_y_0_yyyyz_xxxyy[k];

                g_y_0_xyyyyz_xxyz[k] = -g_y_0_yyyyz_xxyz[k] * ab_x + g_y_0_yyyyz_xxxyz[k];

                g_y_0_xyyyyz_xxzz[k] = -g_y_0_yyyyz_xxzz[k] * ab_x + g_y_0_yyyyz_xxxzz[k];

                g_y_0_xyyyyz_xyyy[k] = -g_y_0_yyyyz_xyyy[k] * ab_x + g_y_0_yyyyz_xxyyy[k];

                g_y_0_xyyyyz_xyyz[k] = -g_y_0_yyyyz_xyyz[k] * ab_x + g_y_0_yyyyz_xxyyz[k];

                g_y_0_xyyyyz_xyzz[k] = -g_y_0_yyyyz_xyzz[k] * ab_x + g_y_0_yyyyz_xxyzz[k];

                g_y_0_xyyyyz_xzzz[k] = -g_y_0_yyyyz_xzzz[k] * ab_x + g_y_0_yyyyz_xxzzz[k];

                g_y_0_xyyyyz_yyyy[k] = -g_y_0_yyyyz_yyyy[k] * ab_x + g_y_0_yyyyz_xyyyy[k];

                g_y_0_xyyyyz_yyyz[k] = -g_y_0_yyyyz_yyyz[k] * ab_x + g_y_0_yyyyz_xyyyz[k];

                g_y_0_xyyyyz_yyzz[k] = -g_y_0_yyyyz_yyzz[k] * ab_x + g_y_0_yyyyz_xyyzz[k];

                g_y_0_xyyyyz_yzzz[k] = -g_y_0_yyyyz_yzzz[k] * ab_x + g_y_0_yyyyz_xyzzz[k];

                g_y_0_xyyyyz_zzzz[k] = -g_y_0_yyyyz_zzzz[k] * ab_x + g_y_0_yyyyz_xzzzz[k];
            }

            /// Set up 675-690 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 675 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 676 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 677 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 678 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 679 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 680 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 681 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 682 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 683 * ccomps * dcomps);

            auto g_y_0_xyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 684 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 685 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 686 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 687 * ccomps * dcomps);

            auto g_y_0_xyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 688 * ccomps * dcomps);

            auto g_y_0_xyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyyzz_xxxx, g_y_0_xyyyzz_xxxy, g_y_0_xyyyzz_xxxz, g_y_0_xyyyzz_xxyy, g_y_0_xyyyzz_xxyz, g_y_0_xyyyzz_xxzz, g_y_0_xyyyzz_xyyy, g_y_0_xyyyzz_xyyz, g_y_0_xyyyzz_xyzz, g_y_0_xyyyzz_xzzz, g_y_0_xyyyzz_yyyy, g_y_0_xyyyzz_yyyz, g_y_0_xyyyzz_yyzz, g_y_0_xyyyzz_yzzz, g_y_0_xyyyzz_zzzz, g_y_0_yyyzz_xxxx, g_y_0_yyyzz_xxxxx, g_y_0_yyyzz_xxxxy, g_y_0_yyyzz_xxxxz, g_y_0_yyyzz_xxxy, g_y_0_yyyzz_xxxyy, g_y_0_yyyzz_xxxyz, g_y_0_yyyzz_xxxz, g_y_0_yyyzz_xxxzz, g_y_0_yyyzz_xxyy, g_y_0_yyyzz_xxyyy, g_y_0_yyyzz_xxyyz, g_y_0_yyyzz_xxyz, g_y_0_yyyzz_xxyzz, g_y_0_yyyzz_xxzz, g_y_0_yyyzz_xxzzz, g_y_0_yyyzz_xyyy, g_y_0_yyyzz_xyyyy, g_y_0_yyyzz_xyyyz, g_y_0_yyyzz_xyyz, g_y_0_yyyzz_xyyzz, g_y_0_yyyzz_xyzz, g_y_0_yyyzz_xyzzz, g_y_0_yyyzz_xzzz, g_y_0_yyyzz_xzzzz, g_y_0_yyyzz_yyyy, g_y_0_yyyzz_yyyz, g_y_0_yyyzz_yyzz, g_y_0_yyyzz_yzzz, g_y_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyyzz_xxxx[k] = -g_y_0_yyyzz_xxxx[k] * ab_x + g_y_0_yyyzz_xxxxx[k];

                g_y_0_xyyyzz_xxxy[k] = -g_y_0_yyyzz_xxxy[k] * ab_x + g_y_0_yyyzz_xxxxy[k];

                g_y_0_xyyyzz_xxxz[k] = -g_y_0_yyyzz_xxxz[k] * ab_x + g_y_0_yyyzz_xxxxz[k];

                g_y_0_xyyyzz_xxyy[k] = -g_y_0_yyyzz_xxyy[k] * ab_x + g_y_0_yyyzz_xxxyy[k];

                g_y_0_xyyyzz_xxyz[k] = -g_y_0_yyyzz_xxyz[k] * ab_x + g_y_0_yyyzz_xxxyz[k];

                g_y_0_xyyyzz_xxzz[k] = -g_y_0_yyyzz_xxzz[k] * ab_x + g_y_0_yyyzz_xxxzz[k];

                g_y_0_xyyyzz_xyyy[k] = -g_y_0_yyyzz_xyyy[k] * ab_x + g_y_0_yyyzz_xxyyy[k];

                g_y_0_xyyyzz_xyyz[k] = -g_y_0_yyyzz_xyyz[k] * ab_x + g_y_0_yyyzz_xxyyz[k];

                g_y_0_xyyyzz_xyzz[k] = -g_y_0_yyyzz_xyzz[k] * ab_x + g_y_0_yyyzz_xxyzz[k];

                g_y_0_xyyyzz_xzzz[k] = -g_y_0_yyyzz_xzzz[k] * ab_x + g_y_0_yyyzz_xxzzz[k];

                g_y_0_xyyyzz_yyyy[k] = -g_y_0_yyyzz_yyyy[k] * ab_x + g_y_0_yyyzz_xyyyy[k];

                g_y_0_xyyyzz_yyyz[k] = -g_y_0_yyyzz_yyyz[k] * ab_x + g_y_0_yyyzz_xyyyz[k];

                g_y_0_xyyyzz_yyzz[k] = -g_y_0_yyyzz_yyzz[k] * ab_x + g_y_0_yyyzz_xyyzz[k];

                g_y_0_xyyyzz_yzzz[k] = -g_y_0_yyyzz_yzzz[k] * ab_x + g_y_0_yyyzz_xyzzz[k];

                g_y_0_xyyyzz_zzzz[k] = -g_y_0_yyyzz_zzzz[k] * ab_x + g_y_0_yyyzz_xzzzz[k];
            }

            /// Set up 690-705 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 690 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 691 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 692 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 693 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 694 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 695 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 696 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 697 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 698 * ccomps * dcomps);

            auto g_y_0_xyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 699 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 700 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 701 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 702 * ccomps * dcomps);

            auto g_y_0_xyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 703 * ccomps * dcomps);

            auto g_y_0_xyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 704 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyyzzz_xxxx, g_y_0_xyyzzz_xxxy, g_y_0_xyyzzz_xxxz, g_y_0_xyyzzz_xxyy, g_y_0_xyyzzz_xxyz, g_y_0_xyyzzz_xxzz, g_y_0_xyyzzz_xyyy, g_y_0_xyyzzz_xyyz, g_y_0_xyyzzz_xyzz, g_y_0_xyyzzz_xzzz, g_y_0_xyyzzz_yyyy, g_y_0_xyyzzz_yyyz, g_y_0_xyyzzz_yyzz, g_y_0_xyyzzz_yzzz, g_y_0_xyyzzz_zzzz, g_y_0_yyzzz_xxxx, g_y_0_yyzzz_xxxxx, g_y_0_yyzzz_xxxxy, g_y_0_yyzzz_xxxxz, g_y_0_yyzzz_xxxy, g_y_0_yyzzz_xxxyy, g_y_0_yyzzz_xxxyz, g_y_0_yyzzz_xxxz, g_y_0_yyzzz_xxxzz, g_y_0_yyzzz_xxyy, g_y_0_yyzzz_xxyyy, g_y_0_yyzzz_xxyyz, g_y_0_yyzzz_xxyz, g_y_0_yyzzz_xxyzz, g_y_0_yyzzz_xxzz, g_y_0_yyzzz_xxzzz, g_y_0_yyzzz_xyyy, g_y_0_yyzzz_xyyyy, g_y_0_yyzzz_xyyyz, g_y_0_yyzzz_xyyz, g_y_0_yyzzz_xyyzz, g_y_0_yyzzz_xyzz, g_y_0_yyzzz_xyzzz, g_y_0_yyzzz_xzzz, g_y_0_yyzzz_xzzzz, g_y_0_yyzzz_yyyy, g_y_0_yyzzz_yyyz, g_y_0_yyzzz_yyzz, g_y_0_yyzzz_yzzz, g_y_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyyzzz_xxxx[k] = -g_y_0_yyzzz_xxxx[k] * ab_x + g_y_0_yyzzz_xxxxx[k];

                g_y_0_xyyzzz_xxxy[k] = -g_y_0_yyzzz_xxxy[k] * ab_x + g_y_0_yyzzz_xxxxy[k];

                g_y_0_xyyzzz_xxxz[k] = -g_y_0_yyzzz_xxxz[k] * ab_x + g_y_0_yyzzz_xxxxz[k];

                g_y_0_xyyzzz_xxyy[k] = -g_y_0_yyzzz_xxyy[k] * ab_x + g_y_0_yyzzz_xxxyy[k];

                g_y_0_xyyzzz_xxyz[k] = -g_y_0_yyzzz_xxyz[k] * ab_x + g_y_0_yyzzz_xxxyz[k];

                g_y_0_xyyzzz_xxzz[k] = -g_y_0_yyzzz_xxzz[k] * ab_x + g_y_0_yyzzz_xxxzz[k];

                g_y_0_xyyzzz_xyyy[k] = -g_y_0_yyzzz_xyyy[k] * ab_x + g_y_0_yyzzz_xxyyy[k];

                g_y_0_xyyzzz_xyyz[k] = -g_y_0_yyzzz_xyyz[k] * ab_x + g_y_0_yyzzz_xxyyz[k];

                g_y_0_xyyzzz_xyzz[k] = -g_y_0_yyzzz_xyzz[k] * ab_x + g_y_0_yyzzz_xxyzz[k];

                g_y_0_xyyzzz_xzzz[k] = -g_y_0_yyzzz_xzzz[k] * ab_x + g_y_0_yyzzz_xxzzz[k];

                g_y_0_xyyzzz_yyyy[k] = -g_y_0_yyzzz_yyyy[k] * ab_x + g_y_0_yyzzz_xyyyy[k];

                g_y_0_xyyzzz_yyyz[k] = -g_y_0_yyzzz_yyyz[k] * ab_x + g_y_0_yyzzz_xyyyz[k];

                g_y_0_xyyzzz_yyzz[k] = -g_y_0_yyzzz_yyzz[k] * ab_x + g_y_0_yyzzz_xyyzz[k];

                g_y_0_xyyzzz_yzzz[k] = -g_y_0_yyzzz_yzzz[k] * ab_x + g_y_0_yyzzz_xyzzz[k];

                g_y_0_xyyzzz_zzzz[k] = -g_y_0_yyzzz_zzzz[k] * ab_x + g_y_0_yyzzz_xzzzz[k];
            }

            /// Set up 705-720 components of targeted buffer : cbuffer.data(

            auto g_y_0_xyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 705 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 706 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 707 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 708 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 709 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 710 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 711 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 712 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 713 * ccomps * dcomps);

            auto g_y_0_xyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 714 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 715 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 716 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 717 * ccomps * dcomps);

            auto g_y_0_xyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 718 * ccomps * dcomps);

            auto g_y_0_xyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xyzzzz_xxxx, g_y_0_xyzzzz_xxxy, g_y_0_xyzzzz_xxxz, g_y_0_xyzzzz_xxyy, g_y_0_xyzzzz_xxyz, g_y_0_xyzzzz_xxzz, g_y_0_xyzzzz_xyyy, g_y_0_xyzzzz_xyyz, g_y_0_xyzzzz_xyzz, g_y_0_xyzzzz_xzzz, g_y_0_xyzzzz_yyyy, g_y_0_xyzzzz_yyyz, g_y_0_xyzzzz_yyzz, g_y_0_xyzzzz_yzzz, g_y_0_xyzzzz_zzzz, g_y_0_yzzzz_xxxx, g_y_0_yzzzz_xxxxx, g_y_0_yzzzz_xxxxy, g_y_0_yzzzz_xxxxz, g_y_0_yzzzz_xxxy, g_y_0_yzzzz_xxxyy, g_y_0_yzzzz_xxxyz, g_y_0_yzzzz_xxxz, g_y_0_yzzzz_xxxzz, g_y_0_yzzzz_xxyy, g_y_0_yzzzz_xxyyy, g_y_0_yzzzz_xxyyz, g_y_0_yzzzz_xxyz, g_y_0_yzzzz_xxyzz, g_y_0_yzzzz_xxzz, g_y_0_yzzzz_xxzzz, g_y_0_yzzzz_xyyy, g_y_0_yzzzz_xyyyy, g_y_0_yzzzz_xyyyz, g_y_0_yzzzz_xyyz, g_y_0_yzzzz_xyyzz, g_y_0_yzzzz_xyzz, g_y_0_yzzzz_xyzzz, g_y_0_yzzzz_xzzz, g_y_0_yzzzz_xzzzz, g_y_0_yzzzz_yyyy, g_y_0_yzzzz_yyyz, g_y_0_yzzzz_yyzz, g_y_0_yzzzz_yzzz, g_y_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xyzzzz_xxxx[k] = -g_y_0_yzzzz_xxxx[k] * ab_x + g_y_0_yzzzz_xxxxx[k];

                g_y_0_xyzzzz_xxxy[k] = -g_y_0_yzzzz_xxxy[k] * ab_x + g_y_0_yzzzz_xxxxy[k];

                g_y_0_xyzzzz_xxxz[k] = -g_y_0_yzzzz_xxxz[k] * ab_x + g_y_0_yzzzz_xxxxz[k];

                g_y_0_xyzzzz_xxyy[k] = -g_y_0_yzzzz_xxyy[k] * ab_x + g_y_0_yzzzz_xxxyy[k];

                g_y_0_xyzzzz_xxyz[k] = -g_y_0_yzzzz_xxyz[k] * ab_x + g_y_0_yzzzz_xxxyz[k];

                g_y_0_xyzzzz_xxzz[k] = -g_y_0_yzzzz_xxzz[k] * ab_x + g_y_0_yzzzz_xxxzz[k];

                g_y_0_xyzzzz_xyyy[k] = -g_y_0_yzzzz_xyyy[k] * ab_x + g_y_0_yzzzz_xxyyy[k];

                g_y_0_xyzzzz_xyyz[k] = -g_y_0_yzzzz_xyyz[k] * ab_x + g_y_0_yzzzz_xxyyz[k];

                g_y_0_xyzzzz_xyzz[k] = -g_y_0_yzzzz_xyzz[k] * ab_x + g_y_0_yzzzz_xxyzz[k];

                g_y_0_xyzzzz_xzzz[k] = -g_y_0_yzzzz_xzzz[k] * ab_x + g_y_0_yzzzz_xxzzz[k];

                g_y_0_xyzzzz_yyyy[k] = -g_y_0_yzzzz_yyyy[k] * ab_x + g_y_0_yzzzz_xyyyy[k];

                g_y_0_xyzzzz_yyyz[k] = -g_y_0_yzzzz_yyyz[k] * ab_x + g_y_0_yzzzz_xyyyz[k];

                g_y_0_xyzzzz_yyzz[k] = -g_y_0_yzzzz_yyzz[k] * ab_x + g_y_0_yzzzz_xyyzz[k];

                g_y_0_xyzzzz_yzzz[k] = -g_y_0_yzzzz_yzzz[k] * ab_x + g_y_0_yzzzz_xyzzz[k];

                g_y_0_xyzzzz_zzzz[k] = -g_y_0_yzzzz_zzzz[k] * ab_x + g_y_0_yzzzz_xzzzz[k];
            }

            /// Set up 720-735 components of targeted buffer : cbuffer.data(

            auto g_y_0_xzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 720 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 721 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 722 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 723 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 724 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 725 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 726 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 727 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 728 * ccomps * dcomps);

            auto g_y_0_xzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 729 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 730 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 731 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 732 * ccomps * dcomps);

            auto g_y_0_xzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 733 * ccomps * dcomps);

            auto g_y_0_xzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_xzzzzz_xxxx, g_y_0_xzzzzz_xxxy, g_y_0_xzzzzz_xxxz, g_y_0_xzzzzz_xxyy, g_y_0_xzzzzz_xxyz, g_y_0_xzzzzz_xxzz, g_y_0_xzzzzz_xyyy, g_y_0_xzzzzz_xyyz, g_y_0_xzzzzz_xyzz, g_y_0_xzzzzz_xzzz, g_y_0_xzzzzz_yyyy, g_y_0_xzzzzz_yyyz, g_y_0_xzzzzz_yyzz, g_y_0_xzzzzz_yzzz, g_y_0_xzzzzz_zzzz, g_y_0_zzzzz_xxxx, g_y_0_zzzzz_xxxxx, g_y_0_zzzzz_xxxxy, g_y_0_zzzzz_xxxxz, g_y_0_zzzzz_xxxy, g_y_0_zzzzz_xxxyy, g_y_0_zzzzz_xxxyz, g_y_0_zzzzz_xxxz, g_y_0_zzzzz_xxxzz, g_y_0_zzzzz_xxyy, g_y_0_zzzzz_xxyyy, g_y_0_zzzzz_xxyyz, g_y_0_zzzzz_xxyz, g_y_0_zzzzz_xxyzz, g_y_0_zzzzz_xxzz, g_y_0_zzzzz_xxzzz, g_y_0_zzzzz_xyyy, g_y_0_zzzzz_xyyyy, g_y_0_zzzzz_xyyyz, g_y_0_zzzzz_xyyz, g_y_0_zzzzz_xyyzz, g_y_0_zzzzz_xyzz, g_y_0_zzzzz_xyzzz, g_y_0_zzzzz_xzzz, g_y_0_zzzzz_xzzzz, g_y_0_zzzzz_yyyy, g_y_0_zzzzz_yyyz, g_y_0_zzzzz_yyzz, g_y_0_zzzzz_yzzz, g_y_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_xzzzzz_xxxx[k] = -g_y_0_zzzzz_xxxx[k] * ab_x + g_y_0_zzzzz_xxxxx[k];

                g_y_0_xzzzzz_xxxy[k] = -g_y_0_zzzzz_xxxy[k] * ab_x + g_y_0_zzzzz_xxxxy[k];

                g_y_0_xzzzzz_xxxz[k] = -g_y_0_zzzzz_xxxz[k] * ab_x + g_y_0_zzzzz_xxxxz[k];

                g_y_0_xzzzzz_xxyy[k] = -g_y_0_zzzzz_xxyy[k] * ab_x + g_y_0_zzzzz_xxxyy[k];

                g_y_0_xzzzzz_xxyz[k] = -g_y_0_zzzzz_xxyz[k] * ab_x + g_y_0_zzzzz_xxxyz[k];

                g_y_0_xzzzzz_xxzz[k] = -g_y_0_zzzzz_xxzz[k] * ab_x + g_y_0_zzzzz_xxxzz[k];

                g_y_0_xzzzzz_xyyy[k] = -g_y_0_zzzzz_xyyy[k] * ab_x + g_y_0_zzzzz_xxyyy[k];

                g_y_0_xzzzzz_xyyz[k] = -g_y_0_zzzzz_xyyz[k] * ab_x + g_y_0_zzzzz_xxyyz[k];

                g_y_0_xzzzzz_xyzz[k] = -g_y_0_zzzzz_xyzz[k] * ab_x + g_y_0_zzzzz_xxyzz[k];

                g_y_0_xzzzzz_xzzz[k] = -g_y_0_zzzzz_xzzz[k] * ab_x + g_y_0_zzzzz_xxzzz[k];

                g_y_0_xzzzzz_yyyy[k] = -g_y_0_zzzzz_yyyy[k] * ab_x + g_y_0_zzzzz_xyyyy[k];

                g_y_0_xzzzzz_yyyz[k] = -g_y_0_zzzzz_yyyz[k] * ab_x + g_y_0_zzzzz_xyyyz[k];

                g_y_0_xzzzzz_yyzz[k] = -g_y_0_zzzzz_yyzz[k] * ab_x + g_y_0_zzzzz_xyyzz[k];

                g_y_0_xzzzzz_yzzz[k] = -g_y_0_zzzzz_yzzz[k] * ab_x + g_y_0_zzzzz_xyzzz[k];

                g_y_0_xzzzzz_zzzz[k] = -g_y_0_zzzzz_zzzz[k] * ab_x + g_y_0_zzzzz_xzzzz[k];
            }

            /// Set up 735-750 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 735 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 736 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 737 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 738 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 739 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 740 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 741 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 742 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 743 * ccomps * dcomps);

            auto g_y_0_yyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 744 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 745 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 746 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 747 * ccomps * dcomps);

            auto g_y_0_yyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 748 * ccomps * dcomps);

            auto g_y_0_yyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_xxxx, g_y_0_yyyyy_xxxxy, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxxyy, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyyy, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyyy, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_yyyy, g_y_0_yyyyy_yyyyy, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_zzzz, g_y_0_yyyyyy_xxxx, g_y_0_yyyyyy_xxxy, g_y_0_yyyyyy_xxxz, g_y_0_yyyyyy_xxyy, g_y_0_yyyyyy_xxyz, g_y_0_yyyyyy_xxzz, g_y_0_yyyyyy_xyyy, g_y_0_yyyyyy_xyyz, g_y_0_yyyyyy_xyzz, g_y_0_yyyyyy_xzzz, g_y_0_yyyyyy_yyyy, g_y_0_yyyyyy_yyyz, g_y_0_yyyyyy_yyzz, g_y_0_yyyyyy_yzzz, g_y_0_yyyyyy_zzzz, g_yyyyy_xxxx, g_yyyyy_xxxy, g_yyyyy_xxxz, g_yyyyy_xxyy, g_yyyyy_xxyz, g_yyyyy_xxzz, g_yyyyy_xyyy, g_yyyyy_xyyz, g_yyyyy_xyzz, g_yyyyy_xzzz, g_yyyyy_yyyy, g_yyyyy_yyyz, g_yyyyy_yyzz, g_yyyyy_yzzz, g_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyy_xxxx[k] = -g_yyyyy_xxxx[k] - g_y_0_yyyyy_xxxx[k] * ab_y + g_y_0_yyyyy_xxxxy[k];

                g_y_0_yyyyyy_xxxy[k] = -g_yyyyy_xxxy[k] - g_y_0_yyyyy_xxxy[k] * ab_y + g_y_0_yyyyy_xxxyy[k];

                g_y_0_yyyyyy_xxxz[k] = -g_yyyyy_xxxz[k] - g_y_0_yyyyy_xxxz[k] * ab_y + g_y_0_yyyyy_xxxyz[k];

                g_y_0_yyyyyy_xxyy[k] = -g_yyyyy_xxyy[k] - g_y_0_yyyyy_xxyy[k] * ab_y + g_y_0_yyyyy_xxyyy[k];

                g_y_0_yyyyyy_xxyz[k] = -g_yyyyy_xxyz[k] - g_y_0_yyyyy_xxyz[k] * ab_y + g_y_0_yyyyy_xxyyz[k];

                g_y_0_yyyyyy_xxzz[k] = -g_yyyyy_xxzz[k] - g_y_0_yyyyy_xxzz[k] * ab_y + g_y_0_yyyyy_xxyzz[k];

                g_y_0_yyyyyy_xyyy[k] = -g_yyyyy_xyyy[k] - g_y_0_yyyyy_xyyy[k] * ab_y + g_y_0_yyyyy_xyyyy[k];

                g_y_0_yyyyyy_xyyz[k] = -g_yyyyy_xyyz[k] - g_y_0_yyyyy_xyyz[k] * ab_y + g_y_0_yyyyy_xyyyz[k];

                g_y_0_yyyyyy_xyzz[k] = -g_yyyyy_xyzz[k] - g_y_0_yyyyy_xyzz[k] * ab_y + g_y_0_yyyyy_xyyzz[k];

                g_y_0_yyyyyy_xzzz[k] = -g_yyyyy_xzzz[k] - g_y_0_yyyyy_xzzz[k] * ab_y + g_y_0_yyyyy_xyzzz[k];

                g_y_0_yyyyyy_yyyy[k] = -g_yyyyy_yyyy[k] - g_y_0_yyyyy_yyyy[k] * ab_y + g_y_0_yyyyy_yyyyy[k];

                g_y_0_yyyyyy_yyyz[k] = -g_yyyyy_yyyz[k] - g_y_0_yyyyy_yyyz[k] * ab_y + g_y_0_yyyyy_yyyyz[k];

                g_y_0_yyyyyy_yyzz[k] = -g_yyyyy_yyzz[k] - g_y_0_yyyyy_yyzz[k] * ab_y + g_y_0_yyyyy_yyyzz[k];

                g_y_0_yyyyyy_yzzz[k] = -g_yyyyy_yzzz[k] - g_y_0_yyyyy_yzzz[k] * ab_y + g_y_0_yyyyy_yyzzz[k];

                g_y_0_yyyyyy_zzzz[k] = -g_yyyyy_zzzz[k] - g_y_0_yyyyy_zzzz[k] * ab_y + g_y_0_yyyyy_yzzzz[k];
            }

            /// Set up 750-765 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 750 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 751 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 752 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 753 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 754 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 755 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 756 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 757 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 758 * ccomps * dcomps);

            auto g_y_0_yyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 759 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 760 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 761 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 762 * ccomps * dcomps);

            auto g_y_0_yyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 763 * ccomps * dcomps);

            auto g_y_0_yyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 764 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyy_xxxx, g_y_0_yyyyy_xxxxz, g_y_0_yyyyy_xxxy, g_y_0_yyyyy_xxxyz, g_y_0_yyyyy_xxxz, g_y_0_yyyyy_xxxzz, g_y_0_yyyyy_xxyy, g_y_0_yyyyy_xxyyz, g_y_0_yyyyy_xxyz, g_y_0_yyyyy_xxyzz, g_y_0_yyyyy_xxzz, g_y_0_yyyyy_xxzzz, g_y_0_yyyyy_xyyy, g_y_0_yyyyy_xyyyz, g_y_0_yyyyy_xyyz, g_y_0_yyyyy_xyyzz, g_y_0_yyyyy_xyzz, g_y_0_yyyyy_xyzzz, g_y_0_yyyyy_xzzz, g_y_0_yyyyy_xzzzz, g_y_0_yyyyy_yyyy, g_y_0_yyyyy_yyyyz, g_y_0_yyyyy_yyyz, g_y_0_yyyyy_yyyzz, g_y_0_yyyyy_yyzz, g_y_0_yyyyy_yyzzz, g_y_0_yyyyy_yzzz, g_y_0_yyyyy_yzzzz, g_y_0_yyyyy_zzzz, g_y_0_yyyyy_zzzzz, g_y_0_yyyyyz_xxxx, g_y_0_yyyyyz_xxxy, g_y_0_yyyyyz_xxxz, g_y_0_yyyyyz_xxyy, g_y_0_yyyyyz_xxyz, g_y_0_yyyyyz_xxzz, g_y_0_yyyyyz_xyyy, g_y_0_yyyyyz_xyyz, g_y_0_yyyyyz_xyzz, g_y_0_yyyyyz_xzzz, g_y_0_yyyyyz_yyyy, g_y_0_yyyyyz_yyyz, g_y_0_yyyyyz_yyzz, g_y_0_yyyyyz_yzzz, g_y_0_yyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyyz_xxxx[k] = -g_y_0_yyyyy_xxxx[k] * ab_z + g_y_0_yyyyy_xxxxz[k];

                g_y_0_yyyyyz_xxxy[k] = -g_y_0_yyyyy_xxxy[k] * ab_z + g_y_0_yyyyy_xxxyz[k];

                g_y_0_yyyyyz_xxxz[k] = -g_y_0_yyyyy_xxxz[k] * ab_z + g_y_0_yyyyy_xxxzz[k];

                g_y_0_yyyyyz_xxyy[k] = -g_y_0_yyyyy_xxyy[k] * ab_z + g_y_0_yyyyy_xxyyz[k];

                g_y_0_yyyyyz_xxyz[k] = -g_y_0_yyyyy_xxyz[k] * ab_z + g_y_0_yyyyy_xxyzz[k];

                g_y_0_yyyyyz_xxzz[k] = -g_y_0_yyyyy_xxzz[k] * ab_z + g_y_0_yyyyy_xxzzz[k];

                g_y_0_yyyyyz_xyyy[k] = -g_y_0_yyyyy_xyyy[k] * ab_z + g_y_0_yyyyy_xyyyz[k];

                g_y_0_yyyyyz_xyyz[k] = -g_y_0_yyyyy_xyyz[k] * ab_z + g_y_0_yyyyy_xyyzz[k];

                g_y_0_yyyyyz_xyzz[k] = -g_y_0_yyyyy_xyzz[k] * ab_z + g_y_0_yyyyy_xyzzz[k];

                g_y_0_yyyyyz_xzzz[k] = -g_y_0_yyyyy_xzzz[k] * ab_z + g_y_0_yyyyy_xzzzz[k];

                g_y_0_yyyyyz_yyyy[k] = -g_y_0_yyyyy_yyyy[k] * ab_z + g_y_0_yyyyy_yyyyz[k];

                g_y_0_yyyyyz_yyyz[k] = -g_y_0_yyyyy_yyyz[k] * ab_z + g_y_0_yyyyy_yyyzz[k];

                g_y_0_yyyyyz_yyzz[k] = -g_y_0_yyyyy_yyzz[k] * ab_z + g_y_0_yyyyy_yyzzz[k];

                g_y_0_yyyyyz_yzzz[k] = -g_y_0_yyyyy_yzzz[k] * ab_z + g_y_0_yyyyy_yzzzz[k];

                g_y_0_yyyyyz_zzzz[k] = -g_y_0_yyyyy_zzzz[k] * ab_z + g_y_0_yyyyy_zzzzz[k];
            }

            /// Set up 765-780 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 765 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 766 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 767 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 768 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 769 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 770 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 771 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 772 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 773 * ccomps * dcomps);

            auto g_y_0_yyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 774 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 775 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 776 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 777 * ccomps * dcomps);

            auto g_y_0_yyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 778 * ccomps * dcomps);

            auto g_y_0_yyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyyz_xxxx, g_y_0_yyyyz_xxxxz, g_y_0_yyyyz_xxxy, g_y_0_yyyyz_xxxyz, g_y_0_yyyyz_xxxz, g_y_0_yyyyz_xxxzz, g_y_0_yyyyz_xxyy, g_y_0_yyyyz_xxyyz, g_y_0_yyyyz_xxyz, g_y_0_yyyyz_xxyzz, g_y_0_yyyyz_xxzz, g_y_0_yyyyz_xxzzz, g_y_0_yyyyz_xyyy, g_y_0_yyyyz_xyyyz, g_y_0_yyyyz_xyyz, g_y_0_yyyyz_xyyzz, g_y_0_yyyyz_xyzz, g_y_0_yyyyz_xyzzz, g_y_0_yyyyz_xzzz, g_y_0_yyyyz_xzzzz, g_y_0_yyyyz_yyyy, g_y_0_yyyyz_yyyyz, g_y_0_yyyyz_yyyz, g_y_0_yyyyz_yyyzz, g_y_0_yyyyz_yyzz, g_y_0_yyyyz_yyzzz, g_y_0_yyyyz_yzzz, g_y_0_yyyyz_yzzzz, g_y_0_yyyyz_zzzz, g_y_0_yyyyz_zzzzz, g_y_0_yyyyzz_xxxx, g_y_0_yyyyzz_xxxy, g_y_0_yyyyzz_xxxz, g_y_0_yyyyzz_xxyy, g_y_0_yyyyzz_xxyz, g_y_0_yyyyzz_xxzz, g_y_0_yyyyzz_xyyy, g_y_0_yyyyzz_xyyz, g_y_0_yyyyzz_xyzz, g_y_0_yyyyzz_xzzz, g_y_0_yyyyzz_yyyy, g_y_0_yyyyzz_yyyz, g_y_0_yyyyzz_yyzz, g_y_0_yyyyzz_yzzz, g_y_0_yyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyyzz_xxxx[k] = -g_y_0_yyyyz_xxxx[k] * ab_z + g_y_0_yyyyz_xxxxz[k];

                g_y_0_yyyyzz_xxxy[k] = -g_y_0_yyyyz_xxxy[k] * ab_z + g_y_0_yyyyz_xxxyz[k];

                g_y_0_yyyyzz_xxxz[k] = -g_y_0_yyyyz_xxxz[k] * ab_z + g_y_0_yyyyz_xxxzz[k];

                g_y_0_yyyyzz_xxyy[k] = -g_y_0_yyyyz_xxyy[k] * ab_z + g_y_0_yyyyz_xxyyz[k];

                g_y_0_yyyyzz_xxyz[k] = -g_y_0_yyyyz_xxyz[k] * ab_z + g_y_0_yyyyz_xxyzz[k];

                g_y_0_yyyyzz_xxzz[k] = -g_y_0_yyyyz_xxzz[k] * ab_z + g_y_0_yyyyz_xxzzz[k];

                g_y_0_yyyyzz_xyyy[k] = -g_y_0_yyyyz_xyyy[k] * ab_z + g_y_0_yyyyz_xyyyz[k];

                g_y_0_yyyyzz_xyyz[k] = -g_y_0_yyyyz_xyyz[k] * ab_z + g_y_0_yyyyz_xyyzz[k];

                g_y_0_yyyyzz_xyzz[k] = -g_y_0_yyyyz_xyzz[k] * ab_z + g_y_0_yyyyz_xyzzz[k];

                g_y_0_yyyyzz_xzzz[k] = -g_y_0_yyyyz_xzzz[k] * ab_z + g_y_0_yyyyz_xzzzz[k];

                g_y_0_yyyyzz_yyyy[k] = -g_y_0_yyyyz_yyyy[k] * ab_z + g_y_0_yyyyz_yyyyz[k];

                g_y_0_yyyyzz_yyyz[k] = -g_y_0_yyyyz_yyyz[k] * ab_z + g_y_0_yyyyz_yyyzz[k];

                g_y_0_yyyyzz_yyzz[k] = -g_y_0_yyyyz_yyzz[k] * ab_z + g_y_0_yyyyz_yyzzz[k];

                g_y_0_yyyyzz_yzzz[k] = -g_y_0_yyyyz_yzzz[k] * ab_z + g_y_0_yyyyz_yzzzz[k];

                g_y_0_yyyyzz_zzzz[k] = -g_y_0_yyyyz_zzzz[k] * ab_z + g_y_0_yyyyz_zzzzz[k];
            }

            /// Set up 780-795 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 780 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 781 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 782 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 783 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 784 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 785 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 786 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 787 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 788 * ccomps * dcomps);

            auto g_y_0_yyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 789 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 790 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 791 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 792 * ccomps * dcomps);

            auto g_y_0_yyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 793 * ccomps * dcomps);

            auto g_y_0_yyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 794 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyyzz_xxxx, g_y_0_yyyzz_xxxxz, g_y_0_yyyzz_xxxy, g_y_0_yyyzz_xxxyz, g_y_0_yyyzz_xxxz, g_y_0_yyyzz_xxxzz, g_y_0_yyyzz_xxyy, g_y_0_yyyzz_xxyyz, g_y_0_yyyzz_xxyz, g_y_0_yyyzz_xxyzz, g_y_0_yyyzz_xxzz, g_y_0_yyyzz_xxzzz, g_y_0_yyyzz_xyyy, g_y_0_yyyzz_xyyyz, g_y_0_yyyzz_xyyz, g_y_0_yyyzz_xyyzz, g_y_0_yyyzz_xyzz, g_y_0_yyyzz_xyzzz, g_y_0_yyyzz_xzzz, g_y_0_yyyzz_xzzzz, g_y_0_yyyzz_yyyy, g_y_0_yyyzz_yyyyz, g_y_0_yyyzz_yyyz, g_y_0_yyyzz_yyyzz, g_y_0_yyyzz_yyzz, g_y_0_yyyzz_yyzzz, g_y_0_yyyzz_yzzz, g_y_0_yyyzz_yzzzz, g_y_0_yyyzz_zzzz, g_y_0_yyyzz_zzzzz, g_y_0_yyyzzz_xxxx, g_y_0_yyyzzz_xxxy, g_y_0_yyyzzz_xxxz, g_y_0_yyyzzz_xxyy, g_y_0_yyyzzz_xxyz, g_y_0_yyyzzz_xxzz, g_y_0_yyyzzz_xyyy, g_y_0_yyyzzz_xyyz, g_y_0_yyyzzz_xyzz, g_y_0_yyyzzz_xzzz, g_y_0_yyyzzz_yyyy, g_y_0_yyyzzz_yyyz, g_y_0_yyyzzz_yyzz, g_y_0_yyyzzz_yzzz, g_y_0_yyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyyzzz_xxxx[k] = -g_y_0_yyyzz_xxxx[k] * ab_z + g_y_0_yyyzz_xxxxz[k];

                g_y_0_yyyzzz_xxxy[k] = -g_y_0_yyyzz_xxxy[k] * ab_z + g_y_0_yyyzz_xxxyz[k];

                g_y_0_yyyzzz_xxxz[k] = -g_y_0_yyyzz_xxxz[k] * ab_z + g_y_0_yyyzz_xxxzz[k];

                g_y_0_yyyzzz_xxyy[k] = -g_y_0_yyyzz_xxyy[k] * ab_z + g_y_0_yyyzz_xxyyz[k];

                g_y_0_yyyzzz_xxyz[k] = -g_y_0_yyyzz_xxyz[k] * ab_z + g_y_0_yyyzz_xxyzz[k];

                g_y_0_yyyzzz_xxzz[k] = -g_y_0_yyyzz_xxzz[k] * ab_z + g_y_0_yyyzz_xxzzz[k];

                g_y_0_yyyzzz_xyyy[k] = -g_y_0_yyyzz_xyyy[k] * ab_z + g_y_0_yyyzz_xyyyz[k];

                g_y_0_yyyzzz_xyyz[k] = -g_y_0_yyyzz_xyyz[k] * ab_z + g_y_0_yyyzz_xyyzz[k];

                g_y_0_yyyzzz_xyzz[k] = -g_y_0_yyyzz_xyzz[k] * ab_z + g_y_0_yyyzz_xyzzz[k];

                g_y_0_yyyzzz_xzzz[k] = -g_y_0_yyyzz_xzzz[k] * ab_z + g_y_0_yyyzz_xzzzz[k];

                g_y_0_yyyzzz_yyyy[k] = -g_y_0_yyyzz_yyyy[k] * ab_z + g_y_0_yyyzz_yyyyz[k];

                g_y_0_yyyzzz_yyyz[k] = -g_y_0_yyyzz_yyyz[k] * ab_z + g_y_0_yyyzz_yyyzz[k];

                g_y_0_yyyzzz_yyzz[k] = -g_y_0_yyyzz_yyzz[k] * ab_z + g_y_0_yyyzz_yyzzz[k];

                g_y_0_yyyzzz_yzzz[k] = -g_y_0_yyyzz_yzzz[k] * ab_z + g_y_0_yyyzz_yzzzz[k];

                g_y_0_yyyzzz_zzzz[k] = -g_y_0_yyyzz_zzzz[k] * ab_z + g_y_0_yyyzz_zzzzz[k];
            }

            /// Set up 795-810 components of targeted buffer : cbuffer.data(

            auto g_y_0_yyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 795 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 796 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 797 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 798 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 799 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 800 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 801 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 802 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 803 * ccomps * dcomps);

            auto g_y_0_yyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 804 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 805 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 806 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 807 * ccomps * dcomps);

            auto g_y_0_yyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 808 * ccomps * dcomps);

            auto g_y_0_yyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yyzzz_xxxx, g_y_0_yyzzz_xxxxz, g_y_0_yyzzz_xxxy, g_y_0_yyzzz_xxxyz, g_y_0_yyzzz_xxxz, g_y_0_yyzzz_xxxzz, g_y_0_yyzzz_xxyy, g_y_0_yyzzz_xxyyz, g_y_0_yyzzz_xxyz, g_y_0_yyzzz_xxyzz, g_y_0_yyzzz_xxzz, g_y_0_yyzzz_xxzzz, g_y_0_yyzzz_xyyy, g_y_0_yyzzz_xyyyz, g_y_0_yyzzz_xyyz, g_y_0_yyzzz_xyyzz, g_y_0_yyzzz_xyzz, g_y_0_yyzzz_xyzzz, g_y_0_yyzzz_xzzz, g_y_0_yyzzz_xzzzz, g_y_0_yyzzz_yyyy, g_y_0_yyzzz_yyyyz, g_y_0_yyzzz_yyyz, g_y_0_yyzzz_yyyzz, g_y_0_yyzzz_yyzz, g_y_0_yyzzz_yyzzz, g_y_0_yyzzz_yzzz, g_y_0_yyzzz_yzzzz, g_y_0_yyzzz_zzzz, g_y_0_yyzzz_zzzzz, g_y_0_yyzzzz_xxxx, g_y_0_yyzzzz_xxxy, g_y_0_yyzzzz_xxxz, g_y_0_yyzzzz_xxyy, g_y_0_yyzzzz_xxyz, g_y_0_yyzzzz_xxzz, g_y_0_yyzzzz_xyyy, g_y_0_yyzzzz_xyyz, g_y_0_yyzzzz_xyzz, g_y_0_yyzzzz_xzzz, g_y_0_yyzzzz_yyyy, g_y_0_yyzzzz_yyyz, g_y_0_yyzzzz_yyzz, g_y_0_yyzzzz_yzzz, g_y_0_yyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yyzzzz_xxxx[k] = -g_y_0_yyzzz_xxxx[k] * ab_z + g_y_0_yyzzz_xxxxz[k];

                g_y_0_yyzzzz_xxxy[k] = -g_y_0_yyzzz_xxxy[k] * ab_z + g_y_0_yyzzz_xxxyz[k];

                g_y_0_yyzzzz_xxxz[k] = -g_y_0_yyzzz_xxxz[k] * ab_z + g_y_0_yyzzz_xxxzz[k];

                g_y_0_yyzzzz_xxyy[k] = -g_y_0_yyzzz_xxyy[k] * ab_z + g_y_0_yyzzz_xxyyz[k];

                g_y_0_yyzzzz_xxyz[k] = -g_y_0_yyzzz_xxyz[k] * ab_z + g_y_0_yyzzz_xxyzz[k];

                g_y_0_yyzzzz_xxzz[k] = -g_y_0_yyzzz_xxzz[k] * ab_z + g_y_0_yyzzz_xxzzz[k];

                g_y_0_yyzzzz_xyyy[k] = -g_y_0_yyzzz_xyyy[k] * ab_z + g_y_0_yyzzz_xyyyz[k];

                g_y_0_yyzzzz_xyyz[k] = -g_y_0_yyzzz_xyyz[k] * ab_z + g_y_0_yyzzz_xyyzz[k];

                g_y_0_yyzzzz_xyzz[k] = -g_y_0_yyzzz_xyzz[k] * ab_z + g_y_0_yyzzz_xyzzz[k];

                g_y_0_yyzzzz_xzzz[k] = -g_y_0_yyzzz_xzzz[k] * ab_z + g_y_0_yyzzz_xzzzz[k];

                g_y_0_yyzzzz_yyyy[k] = -g_y_0_yyzzz_yyyy[k] * ab_z + g_y_0_yyzzz_yyyyz[k];

                g_y_0_yyzzzz_yyyz[k] = -g_y_0_yyzzz_yyyz[k] * ab_z + g_y_0_yyzzz_yyyzz[k];

                g_y_0_yyzzzz_yyzz[k] = -g_y_0_yyzzz_yyzz[k] * ab_z + g_y_0_yyzzz_yyzzz[k];

                g_y_0_yyzzzz_yzzz[k] = -g_y_0_yyzzz_yzzz[k] * ab_z + g_y_0_yyzzz_yzzzz[k];

                g_y_0_yyzzzz_zzzz[k] = -g_y_0_yyzzz_zzzz[k] * ab_z + g_y_0_yyzzz_zzzzz[k];
            }

            /// Set up 810-825 components of targeted buffer : cbuffer.data(

            auto g_y_0_yzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 810 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 811 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 812 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 813 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 814 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 815 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 816 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 817 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 818 * ccomps * dcomps);

            auto g_y_0_yzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 819 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 820 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 821 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 822 * ccomps * dcomps);

            auto g_y_0_yzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 823 * ccomps * dcomps);

            auto g_y_0_yzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 824 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_yzzzz_xxxx, g_y_0_yzzzz_xxxxz, g_y_0_yzzzz_xxxy, g_y_0_yzzzz_xxxyz, g_y_0_yzzzz_xxxz, g_y_0_yzzzz_xxxzz, g_y_0_yzzzz_xxyy, g_y_0_yzzzz_xxyyz, g_y_0_yzzzz_xxyz, g_y_0_yzzzz_xxyzz, g_y_0_yzzzz_xxzz, g_y_0_yzzzz_xxzzz, g_y_0_yzzzz_xyyy, g_y_0_yzzzz_xyyyz, g_y_0_yzzzz_xyyz, g_y_0_yzzzz_xyyzz, g_y_0_yzzzz_xyzz, g_y_0_yzzzz_xyzzz, g_y_0_yzzzz_xzzz, g_y_0_yzzzz_xzzzz, g_y_0_yzzzz_yyyy, g_y_0_yzzzz_yyyyz, g_y_0_yzzzz_yyyz, g_y_0_yzzzz_yyyzz, g_y_0_yzzzz_yyzz, g_y_0_yzzzz_yyzzz, g_y_0_yzzzz_yzzz, g_y_0_yzzzz_yzzzz, g_y_0_yzzzz_zzzz, g_y_0_yzzzz_zzzzz, g_y_0_yzzzzz_xxxx, g_y_0_yzzzzz_xxxy, g_y_0_yzzzzz_xxxz, g_y_0_yzzzzz_xxyy, g_y_0_yzzzzz_xxyz, g_y_0_yzzzzz_xxzz, g_y_0_yzzzzz_xyyy, g_y_0_yzzzzz_xyyz, g_y_0_yzzzzz_xyzz, g_y_0_yzzzzz_xzzz, g_y_0_yzzzzz_yyyy, g_y_0_yzzzzz_yyyz, g_y_0_yzzzzz_yyzz, g_y_0_yzzzzz_yzzz, g_y_0_yzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_yzzzzz_xxxx[k] = -g_y_0_yzzzz_xxxx[k] * ab_z + g_y_0_yzzzz_xxxxz[k];

                g_y_0_yzzzzz_xxxy[k] = -g_y_0_yzzzz_xxxy[k] * ab_z + g_y_0_yzzzz_xxxyz[k];

                g_y_0_yzzzzz_xxxz[k] = -g_y_0_yzzzz_xxxz[k] * ab_z + g_y_0_yzzzz_xxxzz[k];

                g_y_0_yzzzzz_xxyy[k] = -g_y_0_yzzzz_xxyy[k] * ab_z + g_y_0_yzzzz_xxyyz[k];

                g_y_0_yzzzzz_xxyz[k] = -g_y_0_yzzzz_xxyz[k] * ab_z + g_y_0_yzzzz_xxyzz[k];

                g_y_0_yzzzzz_xxzz[k] = -g_y_0_yzzzz_xxzz[k] * ab_z + g_y_0_yzzzz_xxzzz[k];

                g_y_0_yzzzzz_xyyy[k] = -g_y_0_yzzzz_xyyy[k] * ab_z + g_y_0_yzzzz_xyyyz[k];

                g_y_0_yzzzzz_xyyz[k] = -g_y_0_yzzzz_xyyz[k] * ab_z + g_y_0_yzzzz_xyyzz[k];

                g_y_0_yzzzzz_xyzz[k] = -g_y_0_yzzzz_xyzz[k] * ab_z + g_y_0_yzzzz_xyzzz[k];

                g_y_0_yzzzzz_xzzz[k] = -g_y_0_yzzzz_xzzz[k] * ab_z + g_y_0_yzzzz_xzzzz[k];

                g_y_0_yzzzzz_yyyy[k] = -g_y_0_yzzzz_yyyy[k] * ab_z + g_y_0_yzzzz_yyyyz[k];

                g_y_0_yzzzzz_yyyz[k] = -g_y_0_yzzzz_yyyz[k] * ab_z + g_y_0_yzzzz_yyyzz[k];

                g_y_0_yzzzzz_yyzz[k] = -g_y_0_yzzzz_yyzz[k] * ab_z + g_y_0_yzzzz_yyzzz[k];

                g_y_0_yzzzzz_yzzz[k] = -g_y_0_yzzzz_yzzz[k] * ab_z + g_y_0_yzzzz_yzzzz[k];

                g_y_0_yzzzzz_zzzz[k] = -g_y_0_yzzzz_zzzz[k] * ab_z + g_y_0_yzzzz_zzzzz[k];
            }

            /// Set up 825-840 components of targeted buffer : cbuffer.data(

            auto g_y_0_zzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 825 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 826 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 827 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 828 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 829 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 830 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 831 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 832 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 833 * ccomps * dcomps);

            auto g_y_0_zzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 834 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 835 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 836 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 837 * ccomps * dcomps);

            auto g_y_0_zzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 838 * ccomps * dcomps);

            auto g_y_0_zzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_y_0_zzzzz_xxxx, g_y_0_zzzzz_xxxxz, g_y_0_zzzzz_xxxy, g_y_0_zzzzz_xxxyz, g_y_0_zzzzz_xxxz, g_y_0_zzzzz_xxxzz, g_y_0_zzzzz_xxyy, g_y_0_zzzzz_xxyyz, g_y_0_zzzzz_xxyz, g_y_0_zzzzz_xxyzz, g_y_0_zzzzz_xxzz, g_y_0_zzzzz_xxzzz, g_y_0_zzzzz_xyyy, g_y_0_zzzzz_xyyyz, g_y_0_zzzzz_xyyz, g_y_0_zzzzz_xyyzz, g_y_0_zzzzz_xyzz, g_y_0_zzzzz_xyzzz, g_y_0_zzzzz_xzzz, g_y_0_zzzzz_xzzzz, g_y_0_zzzzz_yyyy, g_y_0_zzzzz_yyyyz, g_y_0_zzzzz_yyyz, g_y_0_zzzzz_yyyzz, g_y_0_zzzzz_yyzz, g_y_0_zzzzz_yyzzz, g_y_0_zzzzz_yzzz, g_y_0_zzzzz_yzzzz, g_y_0_zzzzz_zzzz, g_y_0_zzzzz_zzzzz, g_y_0_zzzzzz_xxxx, g_y_0_zzzzzz_xxxy, g_y_0_zzzzzz_xxxz, g_y_0_zzzzzz_xxyy, g_y_0_zzzzzz_xxyz, g_y_0_zzzzzz_xxzz, g_y_0_zzzzzz_xyyy, g_y_0_zzzzzz_xyyz, g_y_0_zzzzzz_xyzz, g_y_0_zzzzzz_xzzz, g_y_0_zzzzzz_yyyy, g_y_0_zzzzzz_yyyz, g_y_0_zzzzzz_yyzz, g_y_0_zzzzzz_yzzz, g_y_0_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_y_0_zzzzzz_xxxx[k] = -g_y_0_zzzzz_xxxx[k] * ab_z + g_y_0_zzzzz_xxxxz[k];

                g_y_0_zzzzzz_xxxy[k] = -g_y_0_zzzzz_xxxy[k] * ab_z + g_y_0_zzzzz_xxxyz[k];

                g_y_0_zzzzzz_xxxz[k] = -g_y_0_zzzzz_xxxz[k] * ab_z + g_y_0_zzzzz_xxxzz[k];

                g_y_0_zzzzzz_xxyy[k] = -g_y_0_zzzzz_xxyy[k] * ab_z + g_y_0_zzzzz_xxyyz[k];

                g_y_0_zzzzzz_xxyz[k] = -g_y_0_zzzzz_xxyz[k] * ab_z + g_y_0_zzzzz_xxyzz[k];

                g_y_0_zzzzzz_xxzz[k] = -g_y_0_zzzzz_xxzz[k] * ab_z + g_y_0_zzzzz_xxzzz[k];

                g_y_0_zzzzzz_xyyy[k] = -g_y_0_zzzzz_xyyy[k] * ab_z + g_y_0_zzzzz_xyyyz[k];

                g_y_0_zzzzzz_xyyz[k] = -g_y_0_zzzzz_xyyz[k] * ab_z + g_y_0_zzzzz_xyyzz[k];

                g_y_0_zzzzzz_xyzz[k] = -g_y_0_zzzzz_xyzz[k] * ab_z + g_y_0_zzzzz_xyzzz[k];

                g_y_0_zzzzzz_xzzz[k] = -g_y_0_zzzzz_xzzz[k] * ab_z + g_y_0_zzzzz_xzzzz[k];

                g_y_0_zzzzzz_yyyy[k] = -g_y_0_zzzzz_yyyy[k] * ab_z + g_y_0_zzzzz_yyyyz[k];

                g_y_0_zzzzzz_yyyz[k] = -g_y_0_zzzzz_yyyz[k] * ab_z + g_y_0_zzzzz_yyyzz[k];

                g_y_0_zzzzzz_yyzz[k] = -g_y_0_zzzzz_yyzz[k] * ab_z + g_y_0_zzzzz_yyzzz[k];

                g_y_0_zzzzzz_yzzz[k] = -g_y_0_zzzzz_yzzz[k] * ab_z + g_y_0_zzzzz_yzzzz[k];

                g_y_0_zzzzzz_zzzz[k] = -g_y_0_zzzzz_zzzz[k] * ab_z + g_y_0_zzzzz_zzzzz[k];
            }

            /// Set up 840-855 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxx_xxxx = cbuffer.data(ig_geom_10_off + 840 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxy = cbuffer.data(ig_geom_10_off + 841 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxxz = cbuffer.data(ig_geom_10_off + 842 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxyy = cbuffer.data(ig_geom_10_off + 843 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxyz = cbuffer.data(ig_geom_10_off + 844 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xxzz = cbuffer.data(ig_geom_10_off + 845 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyyy = cbuffer.data(ig_geom_10_off + 846 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyyz = cbuffer.data(ig_geom_10_off + 847 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xyzz = cbuffer.data(ig_geom_10_off + 848 * ccomps * dcomps);

            auto g_z_0_xxxxxx_xzzz = cbuffer.data(ig_geom_10_off + 849 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyyy = cbuffer.data(ig_geom_10_off + 850 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyyz = cbuffer.data(ig_geom_10_off + 851 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yyzz = cbuffer.data(ig_geom_10_off + 852 * ccomps * dcomps);

            auto g_z_0_xxxxxx_yzzz = cbuffer.data(ig_geom_10_off + 853 * ccomps * dcomps);

            auto g_z_0_xxxxxx_zzzz = cbuffer.data(ig_geom_10_off + 854 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxx_xxxx, g_z_0_xxxxx_xxxxx, g_z_0_xxxxx_xxxxy, g_z_0_xxxxx_xxxxz, g_z_0_xxxxx_xxxy, g_z_0_xxxxx_xxxyy, g_z_0_xxxxx_xxxyz, g_z_0_xxxxx_xxxz, g_z_0_xxxxx_xxxzz, g_z_0_xxxxx_xxyy, g_z_0_xxxxx_xxyyy, g_z_0_xxxxx_xxyyz, g_z_0_xxxxx_xxyz, g_z_0_xxxxx_xxyzz, g_z_0_xxxxx_xxzz, g_z_0_xxxxx_xxzzz, g_z_0_xxxxx_xyyy, g_z_0_xxxxx_xyyyy, g_z_0_xxxxx_xyyyz, g_z_0_xxxxx_xyyz, g_z_0_xxxxx_xyyzz, g_z_0_xxxxx_xyzz, g_z_0_xxxxx_xyzzz, g_z_0_xxxxx_xzzz, g_z_0_xxxxx_xzzzz, g_z_0_xxxxx_yyyy, g_z_0_xxxxx_yyyz, g_z_0_xxxxx_yyzz, g_z_0_xxxxx_yzzz, g_z_0_xxxxx_zzzz, g_z_0_xxxxxx_xxxx, g_z_0_xxxxxx_xxxy, g_z_0_xxxxxx_xxxz, g_z_0_xxxxxx_xxyy, g_z_0_xxxxxx_xxyz, g_z_0_xxxxxx_xxzz, g_z_0_xxxxxx_xyyy, g_z_0_xxxxxx_xyyz, g_z_0_xxxxxx_xyzz, g_z_0_xxxxxx_xzzz, g_z_0_xxxxxx_yyyy, g_z_0_xxxxxx_yyyz, g_z_0_xxxxxx_yyzz, g_z_0_xxxxxx_yzzz, g_z_0_xxxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxx_xxxx[k] = -g_z_0_xxxxx_xxxx[k] * ab_x + g_z_0_xxxxx_xxxxx[k];

                g_z_0_xxxxxx_xxxy[k] = -g_z_0_xxxxx_xxxy[k] * ab_x + g_z_0_xxxxx_xxxxy[k];

                g_z_0_xxxxxx_xxxz[k] = -g_z_0_xxxxx_xxxz[k] * ab_x + g_z_0_xxxxx_xxxxz[k];

                g_z_0_xxxxxx_xxyy[k] = -g_z_0_xxxxx_xxyy[k] * ab_x + g_z_0_xxxxx_xxxyy[k];

                g_z_0_xxxxxx_xxyz[k] = -g_z_0_xxxxx_xxyz[k] * ab_x + g_z_0_xxxxx_xxxyz[k];

                g_z_0_xxxxxx_xxzz[k] = -g_z_0_xxxxx_xxzz[k] * ab_x + g_z_0_xxxxx_xxxzz[k];

                g_z_0_xxxxxx_xyyy[k] = -g_z_0_xxxxx_xyyy[k] * ab_x + g_z_0_xxxxx_xxyyy[k];

                g_z_0_xxxxxx_xyyz[k] = -g_z_0_xxxxx_xyyz[k] * ab_x + g_z_0_xxxxx_xxyyz[k];

                g_z_0_xxxxxx_xyzz[k] = -g_z_0_xxxxx_xyzz[k] * ab_x + g_z_0_xxxxx_xxyzz[k];

                g_z_0_xxxxxx_xzzz[k] = -g_z_0_xxxxx_xzzz[k] * ab_x + g_z_0_xxxxx_xxzzz[k];

                g_z_0_xxxxxx_yyyy[k] = -g_z_0_xxxxx_yyyy[k] * ab_x + g_z_0_xxxxx_xyyyy[k];

                g_z_0_xxxxxx_yyyz[k] = -g_z_0_xxxxx_yyyz[k] * ab_x + g_z_0_xxxxx_xyyyz[k];

                g_z_0_xxxxxx_yyzz[k] = -g_z_0_xxxxx_yyzz[k] * ab_x + g_z_0_xxxxx_xyyzz[k];

                g_z_0_xxxxxx_yzzz[k] = -g_z_0_xxxxx_yzzz[k] * ab_x + g_z_0_xxxxx_xyzzz[k];

                g_z_0_xxxxxx_zzzz[k] = -g_z_0_xxxxx_zzzz[k] * ab_x + g_z_0_xxxxx_xzzzz[k];
            }

            /// Set up 855-870 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxy_xxxx = cbuffer.data(ig_geom_10_off + 855 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxy = cbuffer.data(ig_geom_10_off + 856 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxxz = cbuffer.data(ig_geom_10_off + 857 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxyy = cbuffer.data(ig_geom_10_off + 858 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxyz = cbuffer.data(ig_geom_10_off + 859 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xxzz = cbuffer.data(ig_geom_10_off + 860 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyyy = cbuffer.data(ig_geom_10_off + 861 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyyz = cbuffer.data(ig_geom_10_off + 862 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xyzz = cbuffer.data(ig_geom_10_off + 863 * ccomps * dcomps);

            auto g_z_0_xxxxxy_xzzz = cbuffer.data(ig_geom_10_off + 864 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyyy = cbuffer.data(ig_geom_10_off + 865 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyyz = cbuffer.data(ig_geom_10_off + 866 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yyzz = cbuffer.data(ig_geom_10_off + 867 * ccomps * dcomps);

            auto g_z_0_xxxxxy_yzzz = cbuffer.data(ig_geom_10_off + 868 * ccomps * dcomps);

            auto g_z_0_xxxxxy_zzzz = cbuffer.data(ig_geom_10_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxxy_xxxx, g_z_0_xxxxxy_xxxy, g_z_0_xxxxxy_xxxz, g_z_0_xxxxxy_xxyy, g_z_0_xxxxxy_xxyz, g_z_0_xxxxxy_xxzz, g_z_0_xxxxxy_xyyy, g_z_0_xxxxxy_xyyz, g_z_0_xxxxxy_xyzz, g_z_0_xxxxxy_xzzz, g_z_0_xxxxxy_yyyy, g_z_0_xxxxxy_yyyz, g_z_0_xxxxxy_yyzz, g_z_0_xxxxxy_yzzz, g_z_0_xxxxxy_zzzz, g_z_0_xxxxy_xxxx, g_z_0_xxxxy_xxxxx, g_z_0_xxxxy_xxxxy, g_z_0_xxxxy_xxxxz, g_z_0_xxxxy_xxxy, g_z_0_xxxxy_xxxyy, g_z_0_xxxxy_xxxyz, g_z_0_xxxxy_xxxz, g_z_0_xxxxy_xxxzz, g_z_0_xxxxy_xxyy, g_z_0_xxxxy_xxyyy, g_z_0_xxxxy_xxyyz, g_z_0_xxxxy_xxyz, g_z_0_xxxxy_xxyzz, g_z_0_xxxxy_xxzz, g_z_0_xxxxy_xxzzz, g_z_0_xxxxy_xyyy, g_z_0_xxxxy_xyyyy, g_z_0_xxxxy_xyyyz, g_z_0_xxxxy_xyyz, g_z_0_xxxxy_xyyzz, g_z_0_xxxxy_xyzz, g_z_0_xxxxy_xyzzz, g_z_0_xxxxy_xzzz, g_z_0_xxxxy_xzzzz, g_z_0_xxxxy_yyyy, g_z_0_xxxxy_yyyz, g_z_0_xxxxy_yyzz, g_z_0_xxxxy_yzzz, g_z_0_xxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxy_xxxx[k] = -g_z_0_xxxxy_xxxx[k] * ab_x + g_z_0_xxxxy_xxxxx[k];

                g_z_0_xxxxxy_xxxy[k] = -g_z_0_xxxxy_xxxy[k] * ab_x + g_z_0_xxxxy_xxxxy[k];

                g_z_0_xxxxxy_xxxz[k] = -g_z_0_xxxxy_xxxz[k] * ab_x + g_z_0_xxxxy_xxxxz[k];

                g_z_0_xxxxxy_xxyy[k] = -g_z_0_xxxxy_xxyy[k] * ab_x + g_z_0_xxxxy_xxxyy[k];

                g_z_0_xxxxxy_xxyz[k] = -g_z_0_xxxxy_xxyz[k] * ab_x + g_z_0_xxxxy_xxxyz[k];

                g_z_0_xxxxxy_xxzz[k] = -g_z_0_xxxxy_xxzz[k] * ab_x + g_z_0_xxxxy_xxxzz[k];

                g_z_0_xxxxxy_xyyy[k] = -g_z_0_xxxxy_xyyy[k] * ab_x + g_z_0_xxxxy_xxyyy[k];

                g_z_0_xxxxxy_xyyz[k] = -g_z_0_xxxxy_xyyz[k] * ab_x + g_z_0_xxxxy_xxyyz[k];

                g_z_0_xxxxxy_xyzz[k] = -g_z_0_xxxxy_xyzz[k] * ab_x + g_z_0_xxxxy_xxyzz[k];

                g_z_0_xxxxxy_xzzz[k] = -g_z_0_xxxxy_xzzz[k] * ab_x + g_z_0_xxxxy_xxzzz[k];

                g_z_0_xxxxxy_yyyy[k] = -g_z_0_xxxxy_yyyy[k] * ab_x + g_z_0_xxxxy_xyyyy[k];

                g_z_0_xxxxxy_yyyz[k] = -g_z_0_xxxxy_yyyz[k] * ab_x + g_z_0_xxxxy_xyyyz[k];

                g_z_0_xxxxxy_yyzz[k] = -g_z_0_xxxxy_yyzz[k] * ab_x + g_z_0_xxxxy_xyyzz[k];

                g_z_0_xxxxxy_yzzz[k] = -g_z_0_xxxxy_yzzz[k] * ab_x + g_z_0_xxxxy_xyzzz[k];

                g_z_0_xxxxxy_zzzz[k] = -g_z_0_xxxxy_zzzz[k] * ab_x + g_z_0_xxxxy_xzzzz[k];
            }

            /// Set up 870-885 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxxz_xxxx = cbuffer.data(ig_geom_10_off + 870 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxy = cbuffer.data(ig_geom_10_off + 871 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxxz = cbuffer.data(ig_geom_10_off + 872 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxyy = cbuffer.data(ig_geom_10_off + 873 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxyz = cbuffer.data(ig_geom_10_off + 874 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xxzz = cbuffer.data(ig_geom_10_off + 875 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyyy = cbuffer.data(ig_geom_10_off + 876 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyyz = cbuffer.data(ig_geom_10_off + 877 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xyzz = cbuffer.data(ig_geom_10_off + 878 * ccomps * dcomps);

            auto g_z_0_xxxxxz_xzzz = cbuffer.data(ig_geom_10_off + 879 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyyy = cbuffer.data(ig_geom_10_off + 880 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyyz = cbuffer.data(ig_geom_10_off + 881 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yyzz = cbuffer.data(ig_geom_10_off + 882 * ccomps * dcomps);

            auto g_z_0_xxxxxz_yzzz = cbuffer.data(ig_geom_10_off + 883 * ccomps * dcomps);

            auto g_z_0_xxxxxz_zzzz = cbuffer.data(ig_geom_10_off + 884 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxxz_xxxx, g_z_0_xxxxxz_xxxy, g_z_0_xxxxxz_xxxz, g_z_0_xxxxxz_xxyy, g_z_0_xxxxxz_xxyz, g_z_0_xxxxxz_xxzz, g_z_0_xxxxxz_xyyy, g_z_0_xxxxxz_xyyz, g_z_0_xxxxxz_xyzz, g_z_0_xxxxxz_xzzz, g_z_0_xxxxxz_yyyy, g_z_0_xxxxxz_yyyz, g_z_0_xxxxxz_yyzz, g_z_0_xxxxxz_yzzz, g_z_0_xxxxxz_zzzz, g_z_0_xxxxz_xxxx, g_z_0_xxxxz_xxxxx, g_z_0_xxxxz_xxxxy, g_z_0_xxxxz_xxxxz, g_z_0_xxxxz_xxxy, g_z_0_xxxxz_xxxyy, g_z_0_xxxxz_xxxyz, g_z_0_xxxxz_xxxz, g_z_0_xxxxz_xxxzz, g_z_0_xxxxz_xxyy, g_z_0_xxxxz_xxyyy, g_z_0_xxxxz_xxyyz, g_z_0_xxxxz_xxyz, g_z_0_xxxxz_xxyzz, g_z_0_xxxxz_xxzz, g_z_0_xxxxz_xxzzz, g_z_0_xxxxz_xyyy, g_z_0_xxxxz_xyyyy, g_z_0_xxxxz_xyyyz, g_z_0_xxxxz_xyyz, g_z_0_xxxxz_xyyzz, g_z_0_xxxxz_xyzz, g_z_0_xxxxz_xyzzz, g_z_0_xxxxz_xzzz, g_z_0_xxxxz_xzzzz, g_z_0_xxxxz_yyyy, g_z_0_xxxxz_yyyz, g_z_0_xxxxz_yyzz, g_z_0_xxxxz_yzzz, g_z_0_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxxz_xxxx[k] = -g_z_0_xxxxz_xxxx[k] * ab_x + g_z_0_xxxxz_xxxxx[k];

                g_z_0_xxxxxz_xxxy[k] = -g_z_0_xxxxz_xxxy[k] * ab_x + g_z_0_xxxxz_xxxxy[k];

                g_z_0_xxxxxz_xxxz[k] = -g_z_0_xxxxz_xxxz[k] * ab_x + g_z_0_xxxxz_xxxxz[k];

                g_z_0_xxxxxz_xxyy[k] = -g_z_0_xxxxz_xxyy[k] * ab_x + g_z_0_xxxxz_xxxyy[k];

                g_z_0_xxxxxz_xxyz[k] = -g_z_0_xxxxz_xxyz[k] * ab_x + g_z_0_xxxxz_xxxyz[k];

                g_z_0_xxxxxz_xxzz[k] = -g_z_0_xxxxz_xxzz[k] * ab_x + g_z_0_xxxxz_xxxzz[k];

                g_z_0_xxxxxz_xyyy[k] = -g_z_0_xxxxz_xyyy[k] * ab_x + g_z_0_xxxxz_xxyyy[k];

                g_z_0_xxxxxz_xyyz[k] = -g_z_0_xxxxz_xyyz[k] * ab_x + g_z_0_xxxxz_xxyyz[k];

                g_z_0_xxxxxz_xyzz[k] = -g_z_0_xxxxz_xyzz[k] * ab_x + g_z_0_xxxxz_xxyzz[k];

                g_z_0_xxxxxz_xzzz[k] = -g_z_0_xxxxz_xzzz[k] * ab_x + g_z_0_xxxxz_xxzzz[k];

                g_z_0_xxxxxz_yyyy[k] = -g_z_0_xxxxz_yyyy[k] * ab_x + g_z_0_xxxxz_xyyyy[k];

                g_z_0_xxxxxz_yyyz[k] = -g_z_0_xxxxz_yyyz[k] * ab_x + g_z_0_xxxxz_xyyyz[k];

                g_z_0_xxxxxz_yyzz[k] = -g_z_0_xxxxz_yyzz[k] * ab_x + g_z_0_xxxxz_xyyzz[k];

                g_z_0_xxxxxz_yzzz[k] = -g_z_0_xxxxz_yzzz[k] * ab_x + g_z_0_xxxxz_xyzzz[k];

                g_z_0_xxxxxz_zzzz[k] = -g_z_0_xxxxz_zzzz[k] * ab_x + g_z_0_xxxxz_xzzzz[k];
            }

            /// Set up 885-900 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyy_xxxx = cbuffer.data(ig_geom_10_off + 885 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxy = cbuffer.data(ig_geom_10_off + 886 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxxz = cbuffer.data(ig_geom_10_off + 887 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxyy = cbuffer.data(ig_geom_10_off + 888 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxyz = cbuffer.data(ig_geom_10_off + 889 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xxzz = cbuffer.data(ig_geom_10_off + 890 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyyy = cbuffer.data(ig_geom_10_off + 891 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyyz = cbuffer.data(ig_geom_10_off + 892 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xyzz = cbuffer.data(ig_geom_10_off + 893 * ccomps * dcomps);

            auto g_z_0_xxxxyy_xzzz = cbuffer.data(ig_geom_10_off + 894 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyyy = cbuffer.data(ig_geom_10_off + 895 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyyz = cbuffer.data(ig_geom_10_off + 896 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yyzz = cbuffer.data(ig_geom_10_off + 897 * ccomps * dcomps);

            auto g_z_0_xxxxyy_yzzz = cbuffer.data(ig_geom_10_off + 898 * ccomps * dcomps);

            auto g_z_0_xxxxyy_zzzz = cbuffer.data(ig_geom_10_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxyy_xxxx, g_z_0_xxxxyy_xxxy, g_z_0_xxxxyy_xxxz, g_z_0_xxxxyy_xxyy, g_z_0_xxxxyy_xxyz, g_z_0_xxxxyy_xxzz, g_z_0_xxxxyy_xyyy, g_z_0_xxxxyy_xyyz, g_z_0_xxxxyy_xyzz, g_z_0_xxxxyy_xzzz, g_z_0_xxxxyy_yyyy, g_z_0_xxxxyy_yyyz, g_z_0_xxxxyy_yyzz, g_z_0_xxxxyy_yzzz, g_z_0_xxxxyy_zzzz, g_z_0_xxxyy_xxxx, g_z_0_xxxyy_xxxxx, g_z_0_xxxyy_xxxxy, g_z_0_xxxyy_xxxxz, g_z_0_xxxyy_xxxy, g_z_0_xxxyy_xxxyy, g_z_0_xxxyy_xxxyz, g_z_0_xxxyy_xxxz, g_z_0_xxxyy_xxxzz, g_z_0_xxxyy_xxyy, g_z_0_xxxyy_xxyyy, g_z_0_xxxyy_xxyyz, g_z_0_xxxyy_xxyz, g_z_0_xxxyy_xxyzz, g_z_0_xxxyy_xxzz, g_z_0_xxxyy_xxzzz, g_z_0_xxxyy_xyyy, g_z_0_xxxyy_xyyyy, g_z_0_xxxyy_xyyyz, g_z_0_xxxyy_xyyz, g_z_0_xxxyy_xyyzz, g_z_0_xxxyy_xyzz, g_z_0_xxxyy_xyzzz, g_z_0_xxxyy_xzzz, g_z_0_xxxyy_xzzzz, g_z_0_xxxyy_yyyy, g_z_0_xxxyy_yyyz, g_z_0_xxxyy_yyzz, g_z_0_xxxyy_yzzz, g_z_0_xxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyy_xxxx[k] = -g_z_0_xxxyy_xxxx[k] * ab_x + g_z_0_xxxyy_xxxxx[k];

                g_z_0_xxxxyy_xxxy[k] = -g_z_0_xxxyy_xxxy[k] * ab_x + g_z_0_xxxyy_xxxxy[k];

                g_z_0_xxxxyy_xxxz[k] = -g_z_0_xxxyy_xxxz[k] * ab_x + g_z_0_xxxyy_xxxxz[k];

                g_z_0_xxxxyy_xxyy[k] = -g_z_0_xxxyy_xxyy[k] * ab_x + g_z_0_xxxyy_xxxyy[k];

                g_z_0_xxxxyy_xxyz[k] = -g_z_0_xxxyy_xxyz[k] * ab_x + g_z_0_xxxyy_xxxyz[k];

                g_z_0_xxxxyy_xxzz[k] = -g_z_0_xxxyy_xxzz[k] * ab_x + g_z_0_xxxyy_xxxzz[k];

                g_z_0_xxxxyy_xyyy[k] = -g_z_0_xxxyy_xyyy[k] * ab_x + g_z_0_xxxyy_xxyyy[k];

                g_z_0_xxxxyy_xyyz[k] = -g_z_0_xxxyy_xyyz[k] * ab_x + g_z_0_xxxyy_xxyyz[k];

                g_z_0_xxxxyy_xyzz[k] = -g_z_0_xxxyy_xyzz[k] * ab_x + g_z_0_xxxyy_xxyzz[k];

                g_z_0_xxxxyy_xzzz[k] = -g_z_0_xxxyy_xzzz[k] * ab_x + g_z_0_xxxyy_xxzzz[k];

                g_z_0_xxxxyy_yyyy[k] = -g_z_0_xxxyy_yyyy[k] * ab_x + g_z_0_xxxyy_xyyyy[k];

                g_z_0_xxxxyy_yyyz[k] = -g_z_0_xxxyy_yyyz[k] * ab_x + g_z_0_xxxyy_xyyyz[k];

                g_z_0_xxxxyy_yyzz[k] = -g_z_0_xxxyy_yyzz[k] * ab_x + g_z_0_xxxyy_xyyzz[k];

                g_z_0_xxxxyy_yzzz[k] = -g_z_0_xxxyy_yzzz[k] * ab_x + g_z_0_xxxyy_xyzzz[k];

                g_z_0_xxxxyy_zzzz[k] = -g_z_0_xxxyy_zzzz[k] * ab_x + g_z_0_xxxyy_xzzzz[k];
            }

            /// Set up 900-915 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxyz_xxxx = cbuffer.data(ig_geom_10_off + 900 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxy = cbuffer.data(ig_geom_10_off + 901 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxxz = cbuffer.data(ig_geom_10_off + 902 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxyy = cbuffer.data(ig_geom_10_off + 903 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxyz = cbuffer.data(ig_geom_10_off + 904 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xxzz = cbuffer.data(ig_geom_10_off + 905 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyyy = cbuffer.data(ig_geom_10_off + 906 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyyz = cbuffer.data(ig_geom_10_off + 907 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xyzz = cbuffer.data(ig_geom_10_off + 908 * ccomps * dcomps);

            auto g_z_0_xxxxyz_xzzz = cbuffer.data(ig_geom_10_off + 909 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyyy = cbuffer.data(ig_geom_10_off + 910 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyyz = cbuffer.data(ig_geom_10_off + 911 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yyzz = cbuffer.data(ig_geom_10_off + 912 * ccomps * dcomps);

            auto g_z_0_xxxxyz_yzzz = cbuffer.data(ig_geom_10_off + 913 * ccomps * dcomps);

            auto g_z_0_xxxxyz_zzzz = cbuffer.data(ig_geom_10_off + 914 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxyz_xxxx, g_z_0_xxxxyz_xxxy, g_z_0_xxxxyz_xxxz, g_z_0_xxxxyz_xxyy, g_z_0_xxxxyz_xxyz, g_z_0_xxxxyz_xxzz, g_z_0_xxxxyz_xyyy, g_z_0_xxxxyz_xyyz, g_z_0_xxxxyz_xyzz, g_z_0_xxxxyz_xzzz, g_z_0_xxxxyz_yyyy, g_z_0_xxxxyz_yyyz, g_z_0_xxxxyz_yyzz, g_z_0_xxxxyz_yzzz, g_z_0_xxxxyz_zzzz, g_z_0_xxxyz_xxxx, g_z_0_xxxyz_xxxxx, g_z_0_xxxyz_xxxxy, g_z_0_xxxyz_xxxxz, g_z_0_xxxyz_xxxy, g_z_0_xxxyz_xxxyy, g_z_0_xxxyz_xxxyz, g_z_0_xxxyz_xxxz, g_z_0_xxxyz_xxxzz, g_z_0_xxxyz_xxyy, g_z_0_xxxyz_xxyyy, g_z_0_xxxyz_xxyyz, g_z_0_xxxyz_xxyz, g_z_0_xxxyz_xxyzz, g_z_0_xxxyz_xxzz, g_z_0_xxxyz_xxzzz, g_z_0_xxxyz_xyyy, g_z_0_xxxyz_xyyyy, g_z_0_xxxyz_xyyyz, g_z_0_xxxyz_xyyz, g_z_0_xxxyz_xyyzz, g_z_0_xxxyz_xyzz, g_z_0_xxxyz_xyzzz, g_z_0_xxxyz_xzzz, g_z_0_xxxyz_xzzzz, g_z_0_xxxyz_yyyy, g_z_0_xxxyz_yyyz, g_z_0_xxxyz_yyzz, g_z_0_xxxyz_yzzz, g_z_0_xxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxyz_xxxx[k] = -g_z_0_xxxyz_xxxx[k] * ab_x + g_z_0_xxxyz_xxxxx[k];

                g_z_0_xxxxyz_xxxy[k] = -g_z_0_xxxyz_xxxy[k] * ab_x + g_z_0_xxxyz_xxxxy[k];

                g_z_0_xxxxyz_xxxz[k] = -g_z_0_xxxyz_xxxz[k] * ab_x + g_z_0_xxxyz_xxxxz[k];

                g_z_0_xxxxyz_xxyy[k] = -g_z_0_xxxyz_xxyy[k] * ab_x + g_z_0_xxxyz_xxxyy[k];

                g_z_0_xxxxyz_xxyz[k] = -g_z_0_xxxyz_xxyz[k] * ab_x + g_z_0_xxxyz_xxxyz[k];

                g_z_0_xxxxyz_xxzz[k] = -g_z_0_xxxyz_xxzz[k] * ab_x + g_z_0_xxxyz_xxxzz[k];

                g_z_0_xxxxyz_xyyy[k] = -g_z_0_xxxyz_xyyy[k] * ab_x + g_z_0_xxxyz_xxyyy[k];

                g_z_0_xxxxyz_xyyz[k] = -g_z_0_xxxyz_xyyz[k] * ab_x + g_z_0_xxxyz_xxyyz[k];

                g_z_0_xxxxyz_xyzz[k] = -g_z_0_xxxyz_xyzz[k] * ab_x + g_z_0_xxxyz_xxyzz[k];

                g_z_0_xxxxyz_xzzz[k] = -g_z_0_xxxyz_xzzz[k] * ab_x + g_z_0_xxxyz_xxzzz[k];

                g_z_0_xxxxyz_yyyy[k] = -g_z_0_xxxyz_yyyy[k] * ab_x + g_z_0_xxxyz_xyyyy[k];

                g_z_0_xxxxyz_yyyz[k] = -g_z_0_xxxyz_yyyz[k] * ab_x + g_z_0_xxxyz_xyyyz[k];

                g_z_0_xxxxyz_yyzz[k] = -g_z_0_xxxyz_yyzz[k] * ab_x + g_z_0_xxxyz_xyyzz[k];

                g_z_0_xxxxyz_yzzz[k] = -g_z_0_xxxyz_yzzz[k] * ab_x + g_z_0_xxxyz_xyzzz[k];

                g_z_0_xxxxyz_zzzz[k] = -g_z_0_xxxyz_zzzz[k] * ab_x + g_z_0_xxxyz_xzzzz[k];
            }

            /// Set up 915-930 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxxzz_xxxx = cbuffer.data(ig_geom_10_off + 915 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxy = cbuffer.data(ig_geom_10_off + 916 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxxz = cbuffer.data(ig_geom_10_off + 917 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxyy = cbuffer.data(ig_geom_10_off + 918 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxyz = cbuffer.data(ig_geom_10_off + 919 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xxzz = cbuffer.data(ig_geom_10_off + 920 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyyy = cbuffer.data(ig_geom_10_off + 921 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyyz = cbuffer.data(ig_geom_10_off + 922 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xyzz = cbuffer.data(ig_geom_10_off + 923 * ccomps * dcomps);

            auto g_z_0_xxxxzz_xzzz = cbuffer.data(ig_geom_10_off + 924 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyyy = cbuffer.data(ig_geom_10_off + 925 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyyz = cbuffer.data(ig_geom_10_off + 926 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yyzz = cbuffer.data(ig_geom_10_off + 927 * ccomps * dcomps);

            auto g_z_0_xxxxzz_yzzz = cbuffer.data(ig_geom_10_off + 928 * ccomps * dcomps);

            auto g_z_0_xxxxzz_zzzz = cbuffer.data(ig_geom_10_off + 929 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxxzz_xxxx, g_z_0_xxxxzz_xxxy, g_z_0_xxxxzz_xxxz, g_z_0_xxxxzz_xxyy, g_z_0_xxxxzz_xxyz, g_z_0_xxxxzz_xxzz, g_z_0_xxxxzz_xyyy, g_z_0_xxxxzz_xyyz, g_z_0_xxxxzz_xyzz, g_z_0_xxxxzz_xzzz, g_z_0_xxxxzz_yyyy, g_z_0_xxxxzz_yyyz, g_z_0_xxxxzz_yyzz, g_z_0_xxxxzz_yzzz, g_z_0_xxxxzz_zzzz, g_z_0_xxxzz_xxxx, g_z_0_xxxzz_xxxxx, g_z_0_xxxzz_xxxxy, g_z_0_xxxzz_xxxxz, g_z_0_xxxzz_xxxy, g_z_0_xxxzz_xxxyy, g_z_0_xxxzz_xxxyz, g_z_0_xxxzz_xxxz, g_z_0_xxxzz_xxxzz, g_z_0_xxxzz_xxyy, g_z_0_xxxzz_xxyyy, g_z_0_xxxzz_xxyyz, g_z_0_xxxzz_xxyz, g_z_0_xxxzz_xxyzz, g_z_0_xxxzz_xxzz, g_z_0_xxxzz_xxzzz, g_z_0_xxxzz_xyyy, g_z_0_xxxzz_xyyyy, g_z_0_xxxzz_xyyyz, g_z_0_xxxzz_xyyz, g_z_0_xxxzz_xyyzz, g_z_0_xxxzz_xyzz, g_z_0_xxxzz_xyzzz, g_z_0_xxxzz_xzzz, g_z_0_xxxzz_xzzzz, g_z_0_xxxzz_yyyy, g_z_0_xxxzz_yyyz, g_z_0_xxxzz_yyzz, g_z_0_xxxzz_yzzz, g_z_0_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxxzz_xxxx[k] = -g_z_0_xxxzz_xxxx[k] * ab_x + g_z_0_xxxzz_xxxxx[k];

                g_z_0_xxxxzz_xxxy[k] = -g_z_0_xxxzz_xxxy[k] * ab_x + g_z_0_xxxzz_xxxxy[k];

                g_z_0_xxxxzz_xxxz[k] = -g_z_0_xxxzz_xxxz[k] * ab_x + g_z_0_xxxzz_xxxxz[k];

                g_z_0_xxxxzz_xxyy[k] = -g_z_0_xxxzz_xxyy[k] * ab_x + g_z_0_xxxzz_xxxyy[k];

                g_z_0_xxxxzz_xxyz[k] = -g_z_0_xxxzz_xxyz[k] * ab_x + g_z_0_xxxzz_xxxyz[k];

                g_z_0_xxxxzz_xxzz[k] = -g_z_0_xxxzz_xxzz[k] * ab_x + g_z_0_xxxzz_xxxzz[k];

                g_z_0_xxxxzz_xyyy[k] = -g_z_0_xxxzz_xyyy[k] * ab_x + g_z_0_xxxzz_xxyyy[k];

                g_z_0_xxxxzz_xyyz[k] = -g_z_0_xxxzz_xyyz[k] * ab_x + g_z_0_xxxzz_xxyyz[k];

                g_z_0_xxxxzz_xyzz[k] = -g_z_0_xxxzz_xyzz[k] * ab_x + g_z_0_xxxzz_xxyzz[k];

                g_z_0_xxxxzz_xzzz[k] = -g_z_0_xxxzz_xzzz[k] * ab_x + g_z_0_xxxzz_xxzzz[k];

                g_z_0_xxxxzz_yyyy[k] = -g_z_0_xxxzz_yyyy[k] * ab_x + g_z_0_xxxzz_xyyyy[k];

                g_z_0_xxxxzz_yyyz[k] = -g_z_0_xxxzz_yyyz[k] * ab_x + g_z_0_xxxzz_xyyyz[k];

                g_z_0_xxxxzz_yyzz[k] = -g_z_0_xxxzz_yyzz[k] * ab_x + g_z_0_xxxzz_xyyzz[k];

                g_z_0_xxxxzz_yzzz[k] = -g_z_0_xxxzz_yzzz[k] * ab_x + g_z_0_xxxzz_xyzzz[k];

                g_z_0_xxxxzz_zzzz[k] = -g_z_0_xxxzz_zzzz[k] * ab_x + g_z_0_xxxzz_xzzzz[k];
            }

            /// Set up 930-945 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyy_xxxx = cbuffer.data(ig_geom_10_off + 930 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxy = cbuffer.data(ig_geom_10_off + 931 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxxz = cbuffer.data(ig_geom_10_off + 932 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxyy = cbuffer.data(ig_geom_10_off + 933 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxyz = cbuffer.data(ig_geom_10_off + 934 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xxzz = cbuffer.data(ig_geom_10_off + 935 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyyy = cbuffer.data(ig_geom_10_off + 936 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyyz = cbuffer.data(ig_geom_10_off + 937 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xyzz = cbuffer.data(ig_geom_10_off + 938 * ccomps * dcomps);

            auto g_z_0_xxxyyy_xzzz = cbuffer.data(ig_geom_10_off + 939 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyyy = cbuffer.data(ig_geom_10_off + 940 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyyz = cbuffer.data(ig_geom_10_off + 941 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yyzz = cbuffer.data(ig_geom_10_off + 942 * ccomps * dcomps);

            auto g_z_0_xxxyyy_yzzz = cbuffer.data(ig_geom_10_off + 943 * ccomps * dcomps);

            auto g_z_0_xxxyyy_zzzz = cbuffer.data(ig_geom_10_off + 944 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyyy_xxxx, g_z_0_xxxyyy_xxxy, g_z_0_xxxyyy_xxxz, g_z_0_xxxyyy_xxyy, g_z_0_xxxyyy_xxyz, g_z_0_xxxyyy_xxzz, g_z_0_xxxyyy_xyyy, g_z_0_xxxyyy_xyyz, g_z_0_xxxyyy_xyzz, g_z_0_xxxyyy_xzzz, g_z_0_xxxyyy_yyyy, g_z_0_xxxyyy_yyyz, g_z_0_xxxyyy_yyzz, g_z_0_xxxyyy_yzzz, g_z_0_xxxyyy_zzzz, g_z_0_xxyyy_xxxx, g_z_0_xxyyy_xxxxx, g_z_0_xxyyy_xxxxy, g_z_0_xxyyy_xxxxz, g_z_0_xxyyy_xxxy, g_z_0_xxyyy_xxxyy, g_z_0_xxyyy_xxxyz, g_z_0_xxyyy_xxxz, g_z_0_xxyyy_xxxzz, g_z_0_xxyyy_xxyy, g_z_0_xxyyy_xxyyy, g_z_0_xxyyy_xxyyz, g_z_0_xxyyy_xxyz, g_z_0_xxyyy_xxyzz, g_z_0_xxyyy_xxzz, g_z_0_xxyyy_xxzzz, g_z_0_xxyyy_xyyy, g_z_0_xxyyy_xyyyy, g_z_0_xxyyy_xyyyz, g_z_0_xxyyy_xyyz, g_z_0_xxyyy_xyyzz, g_z_0_xxyyy_xyzz, g_z_0_xxyyy_xyzzz, g_z_0_xxyyy_xzzz, g_z_0_xxyyy_xzzzz, g_z_0_xxyyy_yyyy, g_z_0_xxyyy_yyyz, g_z_0_xxyyy_yyzz, g_z_0_xxyyy_yzzz, g_z_0_xxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyy_xxxx[k] = -g_z_0_xxyyy_xxxx[k] * ab_x + g_z_0_xxyyy_xxxxx[k];

                g_z_0_xxxyyy_xxxy[k] = -g_z_0_xxyyy_xxxy[k] * ab_x + g_z_0_xxyyy_xxxxy[k];

                g_z_0_xxxyyy_xxxz[k] = -g_z_0_xxyyy_xxxz[k] * ab_x + g_z_0_xxyyy_xxxxz[k];

                g_z_0_xxxyyy_xxyy[k] = -g_z_0_xxyyy_xxyy[k] * ab_x + g_z_0_xxyyy_xxxyy[k];

                g_z_0_xxxyyy_xxyz[k] = -g_z_0_xxyyy_xxyz[k] * ab_x + g_z_0_xxyyy_xxxyz[k];

                g_z_0_xxxyyy_xxzz[k] = -g_z_0_xxyyy_xxzz[k] * ab_x + g_z_0_xxyyy_xxxzz[k];

                g_z_0_xxxyyy_xyyy[k] = -g_z_0_xxyyy_xyyy[k] * ab_x + g_z_0_xxyyy_xxyyy[k];

                g_z_0_xxxyyy_xyyz[k] = -g_z_0_xxyyy_xyyz[k] * ab_x + g_z_0_xxyyy_xxyyz[k];

                g_z_0_xxxyyy_xyzz[k] = -g_z_0_xxyyy_xyzz[k] * ab_x + g_z_0_xxyyy_xxyzz[k];

                g_z_0_xxxyyy_xzzz[k] = -g_z_0_xxyyy_xzzz[k] * ab_x + g_z_0_xxyyy_xxzzz[k];

                g_z_0_xxxyyy_yyyy[k] = -g_z_0_xxyyy_yyyy[k] * ab_x + g_z_0_xxyyy_xyyyy[k];

                g_z_0_xxxyyy_yyyz[k] = -g_z_0_xxyyy_yyyz[k] * ab_x + g_z_0_xxyyy_xyyyz[k];

                g_z_0_xxxyyy_yyzz[k] = -g_z_0_xxyyy_yyzz[k] * ab_x + g_z_0_xxyyy_xyyzz[k];

                g_z_0_xxxyyy_yzzz[k] = -g_z_0_xxyyy_yzzz[k] * ab_x + g_z_0_xxyyy_xyzzz[k];

                g_z_0_xxxyyy_zzzz[k] = -g_z_0_xxyyy_zzzz[k] * ab_x + g_z_0_xxyyy_xzzzz[k];
            }

            /// Set up 945-960 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyyz_xxxx = cbuffer.data(ig_geom_10_off + 945 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxy = cbuffer.data(ig_geom_10_off + 946 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxxz = cbuffer.data(ig_geom_10_off + 947 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxyy = cbuffer.data(ig_geom_10_off + 948 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxyz = cbuffer.data(ig_geom_10_off + 949 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xxzz = cbuffer.data(ig_geom_10_off + 950 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyyy = cbuffer.data(ig_geom_10_off + 951 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyyz = cbuffer.data(ig_geom_10_off + 952 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xyzz = cbuffer.data(ig_geom_10_off + 953 * ccomps * dcomps);

            auto g_z_0_xxxyyz_xzzz = cbuffer.data(ig_geom_10_off + 954 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyyy = cbuffer.data(ig_geom_10_off + 955 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyyz = cbuffer.data(ig_geom_10_off + 956 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yyzz = cbuffer.data(ig_geom_10_off + 957 * ccomps * dcomps);

            auto g_z_0_xxxyyz_yzzz = cbuffer.data(ig_geom_10_off + 958 * ccomps * dcomps);

            auto g_z_0_xxxyyz_zzzz = cbuffer.data(ig_geom_10_off + 959 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyyz_xxxx, g_z_0_xxxyyz_xxxy, g_z_0_xxxyyz_xxxz, g_z_0_xxxyyz_xxyy, g_z_0_xxxyyz_xxyz, g_z_0_xxxyyz_xxzz, g_z_0_xxxyyz_xyyy, g_z_0_xxxyyz_xyyz, g_z_0_xxxyyz_xyzz, g_z_0_xxxyyz_xzzz, g_z_0_xxxyyz_yyyy, g_z_0_xxxyyz_yyyz, g_z_0_xxxyyz_yyzz, g_z_0_xxxyyz_yzzz, g_z_0_xxxyyz_zzzz, g_z_0_xxyyz_xxxx, g_z_0_xxyyz_xxxxx, g_z_0_xxyyz_xxxxy, g_z_0_xxyyz_xxxxz, g_z_0_xxyyz_xxxy, g_z_0_xxyyz_xxxyy, g_z_0_xxyyz_xxxyz, g_z_0_xxyyz_xxxz, g_z_0_xxyyz_xxxzz, g_z_0_xxyyz_xxyy, g_z_0_xxyyz_xxyyy, g_z_0_xxyyz_xxyyz, g_z_0_xxyyz_xxyz, g_z_0_xxyyz_xxyzz, g_z_0_xxyyz_xxzz, g_z_0_xxyyz_xxzzz, g_z_0_xxyyz_xyyy, g_z_0_xxyyz_xyyyy, g_z_0_xxyyz_xyyyz, g_z_0_xxyyz_xyyz, g_z_0_xxyyz_xyyzz, g_z_0_xxyyz_xyzz, g_z_0_xxyyz_xyzzz, g_z_0_xxyyz_xzzz, g_z_0_xxyyz_xzzzz, g_z_0_xxyyz_yyyy, g_z_0_xxyyz_yyyz, g_z_0_xxyyz_yyzz, g_z_0_xxyyz_yzzz, g_z_0_xxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyyz_xxxx[k] = -g_z_0_xxyyz_xxxx[k] * ab_x + g_z_0_xxyyz_xxxxx[k];

                g_z_0_xxxyyz_xxxy[k] = -g_z_0_xxyyz_xxxy[k] * ab_x + g_z_0_xxyyz_xxxxy[k];

                g_z_0_xxxyyz_xxxz[k] = -g_z_0_xxyyz_xxxz[k] * ab_x + g_z_0_xxyyz_xxxxz[k];

                g_z_0_xxxyyz_xxyy[k] = -g_z_0_xxyyz_xxyy[k] * ab_x + g_z_0_xxyyz_xxxyy[k];

                g_z_0_xxxyyz_xxyz[k] = -g_z_0_xxyyz_xxyz[k] * ab_x + g_z_0_xxyyz_xxxyz[k];

                g_z_0_xxxyyz_xxzz[k] = -g_z_0_xxyyz_xxzz[k] * ab_x + g_z_0_xxyyz_xxxzz[k];

                g_z_0_xxxyyz_xyyy[k] = -g_z_0_xxyyz_xyyy[k] * ab_x + g_z_0_xxyyz_xxyyy[k];

                g_z_0_xxxyyz_xyyz[k] = -g_z_0_xxyyz_xyyz[k] * ab_x + g_z_0_xxyyz_xxyyz[k];

                g_z_0_xxxyyz_xyzz[k] = -g_z_0_xxyyz_xyzz[k] * ab_x + g_z_0_xxyyz_xxyzz[k];

                g_z_0_xxxyyz_xzzz[k] = -g_z_0_xxyyz_xzzz[k] * ab_x + g_z_0_xxyyz_xxzzz[k];

                g_z_0_xxxyyz_yyyy[k] = -g_z_0_xxyyz_yyyy[k] * ab_x + g_z_0_xxyyz_xyyyy[k];

                g_z_0_xxxyyz_yyyz[k] = -g_z_0_xxyyz_yyyz[k] * ab_x + g_z_0_xxyyz_xyyyz[k];

                g_z_0_xxxyyz_yyzz[k] = -g_z_0_xxyyz_yyzz[k] * ab_x + g_z_0_xxyyz_xyyzz[k];

                g_z_0_xxxyyz_yzzz[k] = -g_z_0_xxyyz_yzzz[k] * ab_x + g_z_0_xxyyz_xyzzz[k];

                g_z_0_xxxyyz_zzzz[k] = -g_z_0_xxyyz_zzzz[k] * ab_x + g_z_0_xxyyz_xzzzz[k];
            }

            /// Set up 960-975 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxyzz_xxxx = cbuffer.data(ig_geom_10_off + 960 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxy = cbuffer.data(ig_geom_10_off + 961 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxxz = cbuffer.data(ig_geom_10_off + 962 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxyy = cbuffer.data(ig_geom_10_off + 963 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxyz = cbuffer.data(ig_geom_10_off + 964 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xxzz = cbuffer.data(ig_geom_10_off + 965 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyyy = cbuffer.data(ig_geom_10_off + 966 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyyz = cbuffer.data(ig_geom_10_off + 967 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xyzz = cbuffer.data(ig_geom_10_off + 968 * ccomps * dcomps);

            auto g_z_0_xxxyzz_xzzz = cbuffer.data(ig_geom_10_off + 969 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyyy = cbuffer.data(ig_geom_10_off + 970 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyyz = cbuffer.data(ig_geom_10_off + 971 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yyzz = cbuffer.data(ig_geom_10_off + 972 * ccomps * dcomps);

            auto g_z_0_xxxyzz_yzzz = cbuffer.data(ig_geom_10_off + 973 * ccomps * dcomps);

            auto g_z_0_xxxyzz_zzzz = cbuffer.data(ig_geom_10_off + 974 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxyzz_xxxx, g_z_0_xxxyzz_xxxy, g_z_0_xxxyzz_xxxz, g_z_0_xxxyzz_xxyy, g_z_0_xxxyzz_xxyz, g_z_0_xxxyzz_xxzz, g_z_0_xxxyzz_xyyy, g_z_0_xxxyzz_xyyz, g_z_0_xxxyzz_xyzz, g_z_0_xxxyzz_xzzz, g_z_0_xxxyzz_yyyy, g_z_0_xxxyzz_yyyz, g_z_0_xxxyzz_yyzz, g_z_0_xxxyzz_yzzz, g_z_0_xxxyzz_zzzz, g_z_0_xxyzz_xxxx, g_z_0_xxyzz_xxxxx, g_z_0_xxyzz_xxxxy, g_z_0_xxyzz_xxxxz, g_z_0_xxyzz_xxxy, g_z_0_xxyzz_xxxyy, g_z_0_xxyzz_xxxyz, g_z_0_xxyzz_xxxz, g_z_0_xxyzz_xxxzz, g_z_0_xxyzz_xxyy, g_z_0_xxyzz_xxyyy, g_z_0_xxyzz_xxyyz, g_z_0_xxyzz_xxyz, g_z_0_xxyzz_xxyzz, g_z_0_xxyzz_xxzz, g_z_0_xxyzz_xxzzz, g_z_0_xxyzz_xyyy, g_z_0_xxyzz_xyyyy, g_z_0_xxyzz_xyyyz, g_z_0_xxyzz_xyyz, g_z_0_xxyzz_xyyzz, g_z_0_xxyzz_xyzz, g_z_0_xxyzz_xyzzz, g_z_0_xxyzz_xzzz, g_z_0_xxyzz_xzzzz, g_z_0_xxyzz_yyyy, g_z_0_xxyzz_yyyz, g_z_0_xxyzz_yyzz, g_z_0_xxyzz_yzzz, g_z_0_xxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxyzz_xxxx[k] = -g_z_0_xxyzz_xxxx[k] * ab_x + g_z_0_xxyzz_xxxxx[k];

                g_z_0_xxxyzz_xxxy[k] = -g_z_0_xxyzz_xxxy[k] * ab_x + g_z_0_xxyzz_xxxxy[k];

                g_z_0_xxxyzz_xxxz[k] = -g_z_0_xxyzz_xxxz[k] * ab_x + g_z_0_xxyzz_xxxxz[k];

                g_z_0_xxxyzz_xxyy[k] = -g_z_0_xxyzz_xxyy[k] * ab_x + g_z_0_xxyzz_xxxyy[k];

                g_z_0_xxxyzz_xxyz[k] = -g_z_0_xxyzz_xxyz[k] * ab_x + g_z_0_xxyzz_xxxyz[k];

                g_z_0_xxxyzz_xxzz[k] = -g_z_0_xxyzz_xxzz[k] * ab_x + g_z_0_xxyzz_xxxzz[k];

                g_z_0_xxxyzz_xyyy[k] = -g_z_0_xxyzz_xyyy[k] * ab_x + g_z_0_xxyzz_xxyyy[k];

                g_z_0_xxxyzz_xyyz[k] = -g_z_0_xxyzz_xyyz[k] * ab_x + g_z_0_xxyzz_xxyyz[k];

                g_z_0_xxxyzz_xyzz[k] = -g_z_0_xxyzz_xyzz[k] * ab_x + g_z_0_xxyzz_xxyzz[k];

                g_z_0_xxxyzz_xzzz[k] = -g_z_0_xxyzz_xzzz[k] * ab_x + g_z_0_xxyzz_xxzzz[k];

                g_z_0_xxxyzz_yyyy[k] = -g_z_0_xxyzz_yyyy[k] * ab_x + g_z_0_xxyzz_xyyyy[k];

                g_z_0_xxxyzz_yyyz[k] = -g_z_0_xxyzz_yyyz[k] * ab_x + g_z_0_xxyzz_xyyyz[k];

                g_z_0_xxxyzz_yyzz[k] = -g_z_0_xxyzz_yyzz[k] * ab_x + g_z_0_xxyzz_xyyzz[k];

                g_z_0_xxxyzz_yzzz[k] = -g_z_0_xxyzz_yzzz[k] * ab_x + g_z_0_xxyzz_xyzzz[k];

                g_z_0_xxxyzz_zzzz[k] = -g_z_0_xxyzz_zzzz[k] * ab_x + g_z_0_xxyzz_xzzzz[k];
            }

            /// Set up 975-990 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxxzzz_xxxx = cbuffer.data(ig_geom_10_off + 975 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxy = cbuffer.data(ig_geom_10_off + 976 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxxz = cbuffer.data(ig_geom_10_off + 977 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxyy = cbuffer.data(ig_geom_10_off + 978 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxyz = cbuffer.data(ig_geom_10_off + 979 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xxzz = cbuffer.data(ig_geom_10_off + 980 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyyy = cbuffer.data(ig_geom_10_off + 981 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyyz = cbuffer.data(ig_geom_10_off + 982 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xyzz = cbuffer.data(ig_geom_10_off + 983 * ccomps * dcomps);

            auto g_z_0_xxxzzz_xzzz = cbuffer.data(ig_geom_10_off + 984 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyyy = cbuffer.data(ig_geom_10_off + 985 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyyz = cbuffer.data(ig_geom_10_off + 986 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yyzz = cbuffer.data(ig_geom_10_off + 987 * ccomps * dcomps);

            auto g_z_0_xxxzzz_yzzz = cbuffer.data(ig_geom_10_off + 988 * ccomps * dcomps);

            auto g_z_0_xxxzzz_zzzz = cbuffer.data(ig_geom_10_off + 989 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxxzzz_xxxx, g_z_0_xxxzzz_xxxy, g_z_0_xxxzzz_xxxz, g_z_0_xxxzzz_xxyy, g_z_0_xxxzzz_xxyz, g_z_0_xxxzzz_xxzz, g_z_0_xxxzzz_xyyy, g_z_0_xxxzzz_xyyz, g_z_0_xxxzzz_xyzz, g_z_0_xxxzzz_xzzz, g_z_0_xxxzzz_yyyy, g_z_0_xxxzzz_yyyz, g_z_0_xxxzzz_yyzz, g_z_0_xxxzzz_yzzz, g_z_0_xxxzzz_zzzz, g_z_0_xxzzz_xxxx, g_z_0_xxzzz_xxxxx, g_z_0_xxzzz_xxxxy, g_z_0_xxzzz_xxxxz, g_z_0_xxzzz_xxxy, g_z_0_xxzzz_xxxyy, g_z_0_xxzzz_xxxyz, g_z_0_xxzzz_xxxz, g_z_0_xxzzz_xxxzz, g_z_0_xxzzz_xxyy, g_z_0_xxzzz_xxyyy, g_z_0_xxzzz_xxyyz, g_z_0_xxzzz_xxyz, g_z_0_xxzzz_xxyzz, g_z_0_xxzzz_xxzz, g_z_0_xxzzz_xxzzz, g_z_0_xxzzz_xyyy, g_z_0_xxzzz_xyyyy, g_z_0_xxzzz_xyyyz, g_z_0_xxzzz_xyyz, g_z_0_xxzzz_xyyzz, g_z_0_xxzzz_xyzz, g_z_0_xxzzz_xyzzz, g_z_0_xxzzz_xzzz, g_z_0_xxzzz_xzzzz, g_z_0_xxzzz_yyyy, g_z_0_xxzzz_yyyz, g_z_0_xxzzz_yyzz, g_z_0_xxzzz_yzzz, g_z_0_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxxzzz_xxxx[k] = -g_z_0_xxzzz_xxxx[k] * ab_x + g_z_0_xxzzz_xxxxx[k];

                g_z_0_xxxzzz_xxxy[k] = -g_z_0_xxzzz_xxxy[k] * ab_x + g_z_0_xxzzz_xxxxy[k];

                g_z_0_xxxzzz_xxxz[k] = -g_z_0_xxzzz_xxxz[k] * ab_x + g_z_0_xxzzz_xxxxz[k];

                g_z_0_xxxzzz_xxyy[k] = -g_z_0_xxzzz_xxyy[k] * ab_x + g_z_0_xxzzz_xxxyy[k];

                g_z_0_xxxzzz_xxyz[k] = -g_z_0_xxzzz_xxyz[k] * ab_x + g_z_0_xxzzz_xxxyz[k];

                g_z_0_xxxzzz_xxzz[k] = -g_z_0_xxzzz_xxzz[k] * ab_x + g_z_0_xxzzz_xxxzz[k];

                g_z_0_xxxzzz_xyyy[k] = -g_z_0_xxzzz_xyyy[k] * ab_x + g_z_0_xxzzz_xxyyy[k];

                g_z_0_xxxzzz_xyyz[k] = -g_z_0_xxzzz_xyyz[k] * ab_x + g_z_0_xxzzz_xxyyz[k];

                g_z_0_xxxzzz_xyzz[k] = -g_z_0_xxzzz_xyzz[k] * ab_x + g_z_0_xxzzz_xxyzz[k];

                g_z_0_xxxzzz_xzzz[k] = -g_z_0_xxzzz_xzzz[k] * ab_x + g_z_0_xxzzz_xxzzz[k];

                g_z_0_xxxzzz_yyyy[k] = -g_z_0_xxzzz_yyyy[k] * ab_x + g_z_0_xxzzz_xyyyy[k];

                g_z_0_xxxzzz_yyyz[k] = -g_z_0_xxzzz_yyyz[k] * ab_x + g_z_0_xxzzz_xyyyz[k];

                g_z_0_xxxzzz_yyzz[k] = -g_z_0_xxzzz_yyzz[k] * ab_x + g_z_0_xxzzz_xyyzz[k];

                g_z_0_xxxzzz_yzzz[k] = -g_z_0_xxzzz_yzzz[k] * ab_x + g_z_0_xxzzz_xyzzz[k];

                g_z_0_xxxzzz_zzzz[k] = -g_z_0_xxzzz_zzzz[k] * ab_x + g_z_0_xxzzz_xzzzz[k];
            }

            /// Set up 990-1005 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyy_xxxx = cbuffer.data(ig_geom_10_off + 990 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxy = cbuffer.data(ig_geom_10_off + 991 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxxz = cbuffer.data(ig_geom_10_off + 992 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxyy = cbuffer.data(ig_geom_10_off + 993 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxyz = cbuffer.data(ig_geom_10_off + 994 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xxzz = cbuffer.data(ig_geom_10_off + 995 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyyy = cbuffer.data(ig_geom_10_off + 996 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyyz = cbuffer.data(ig_geom_10_off + 997 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xyzz = cbuffer.data(ig_geom_10_off + 998 * ccomps * dcomps);

            auto g_z_0_xxyyyy_xzzz = cbuffer.data(ig_geom_10_off + 999 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyyy = cbuffer.data(ig_geom_10_off + 1000 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyyz = cbuffer.data(ig_geom_10_off + 1001 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yyzz = cbuffer.data(ig_geom_10_off + 1002 * ccomps * dcomps);

            auto g_z_0_xxyyyy_yzzz = cbuffer.data(ig_geom_10_off + 1003 * ccomps * dcomps);

            auto g_z_0_xxyyyy_zzzz = cbuffer.data(ig_geom_10_off + 1004 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyyy_xxxx, g_z_0_xxyyyy_xxxy, g_z_0_xxyyyy_xxxz, g_z_0_xxyyyy_xxyy, g_z_0_xxyyyy_xxyz, g_z_0_xxyyyy_xxzz, g_z_0_xxyyyy_xyyy, g_z_0_xxyyyy_xyyz, g_z_0_xxyyyy_xyzz, g_z_0_xxyyyy_xzzz, g_z_0_xxyyyy_yyyy, g_z_0_xxyyyy_yyyz, g_z_0_xxyyyy_yyzz, g_z_0_xxyyyy_yzzz, g_z_0_xxyyyy_zzzz, g_z_0_xyyyy_xxxx, g_z_0_xyyyy_xxxxx, g_z_0_xyyyy_xxxxy, g_z_0_xyyyy_xxxxz, g_z_0_xyyyy_xxxy, g_z_0_xyyyy_xxxyy, g_z_0_xyyyy_xxxyz, g_z_0_xyyyy_xxxz, g_z_0_xyyyy_xxxzz, g_z_0_xyyyy_xxyy, g_z_0_xyyyy_xxyyy, g_z_0_xyyyy_xxyyz, g_z_0_xyyyy_xxyz, g_z_0_xyyyy_xxyzz, g_z_0_xyyyy_xxzz, g_z_0_xyyyy_xxzzz, g_z_0_xyyyy_xyyy, g_z_0_xyyyy_xyyyy, g_z_0_xyyyy_xyyyz, g_z_0_xyyyy_xyyz, g_z_0_xyyyy_xyyzz, g_z_0_xyyyy_xyzz, g_z_0_xyyyy_xyzzz, g_z_0_xyyyy_xzzz, g_z_0_xyyyy_xzzzz, g_z_0_xyyyy_yyyy, g_z_0_xyyyy_yyyz, g_z_0_xyyyy_yyzz, g_z_0_xyyyy_yzzz, g_z_0_xyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyy_xxxx[k] = -g_z_0_xyyyy_xxxx[k] * ab_x + g_z_0_xyyyy_xxxxx[k];

                g_z_0_xxyyyy_xxxy[k] = -g_z_0_xyyyy_xxxy[k] * ab_x + g_z_0_xyyyy_xxxxy[k];

                g_z_0_xxyyyy_xxxz[k] = -g_z_0_xyyyy_xxxz[k] * ab_x + g_z_0_xyyyy_xxxxz[k];

                g_z_0_xxyyyy_xxyy[k] = -g_z_0_xyyyy_xxyy[k] * ab_x + g_z_0_xyyyy_xxxyy[k];

                g_z_0_xxyyyy_xxyz[k] = -g_z_0_xyyyy_xxyz[k] * ab_x + g_z_0_xyyyy_xxxyz[k];

                g_z_0_xxyyyy_xxzz[k] = -g_z_0_xyyyy_xxzz[k] * ab_x + g_z_0_xyyyy_xxxzz[k];

                g_z_0_xxyyyy_xyyy[k] = -g_z_0_xyyyy_xyyy[k] * ab_x + g_z_0_xyyyy_xxyyy[k];

                g_z_0_xxyyyy_xyyz[k] = -g_z_0_xyyyy_xyyz[k] * ab_x + g_z_0_xyyyy_xxyyz[k];

                g_z_0_xxyyyy_xyzz[k] = -g_z_0_xyyyy_xyzz[k] * ab_x + g_z_0_xyyyy_xxyzz[k];

                g_z_0_xxyyyy_xzzz[k] = -g_z_0_xyyyy_xzzz[k] * ab_x + g_z_0_xyyyy_xxzzz[k];

                g_z_0_xxyyyy_yyyy[k] = -g_z_0_xyyyy_yyyy[k] * ab_x + g_z_0_xyyyy_xyyyy[k];

                g_z_0_xxyyyy_yyyz[k] = -g_z_0_xyyyy_yyyz[k] * ab_x + g_z_0_xyyyy_xyyyz[k];

                g_z_0_xxyyyy_yyzz[k] = -g_z_0_xyyyy_yyzz[k] * ab_x + g_z_0_xyyyy_xyyzz[k];

                g_z_0_xxyyyy_yzzz[k] = -g_z_0_xyyyy_yzzz[k] * ab_x + g_z_0_xyyyy_xyzzz[k];

                g_z_0_xxyyyy_zzzz[k] = -g_z_0_xyyyy_zzzz[k] * ab_x + g_z_0_xyyyy_xzzzz[k];
            }

            /// Set up 1005-1020 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyyz_xxxx = cbuffer.data(ig_geom_10_off + 1005 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxy = cbuffer.data(ig_geom_10_off + 1006 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxxz = cbuffer.data(ig_geom_10_off + 1007 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxyy = cbuffer.data(ig_geom_10_off + 1008 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxyz = cbuffer.data(ig_geom_10_off + 1009 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xxzz = cbuffer.data(ig_geom_10_off + 1010 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyyy = cbuffer.data(ig_geom_10_off + 1011 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyyz = cbuffer.data(ig_geom_10_off + 1012 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xyzz = cbuffer.data(ig_geom_10_off + 1013 * ccomps * dcomps);

            auto g_z_0_xxyyyz_xzzz = cbuffer.data(ig_geom_10_off + 1014 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyyy = cbuffer.data(ig_geom_10_off + 1015 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyyz = cbuffer.data(ig_geom_10_off + 1016 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yyzz = cbuffer.data(ig_geom_10_off + 1017 * ccomps * dcomps);

            auto g_z_0_xxyyyz_yzzz = cbuffer.data(ig_geom_10_off + 1018 * ccomps * dcomps);

            auto g_z_0_xxyyyz_zzzz = cbuffer.data(ig_geom_10_off + 1019 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyyz_xxxx, g_z_0_xxyyyz_xxxy, g_z_0_xxyyyz_xxxz, g_z_0_xxyyyz_xxyy, g_z_0_xxyyyz_xxyz, g_z_0_xxyyyz_xxzz, g_z_0_xxyyyz_xyyy, g_z_0_xxyyyz_xyyz, g_z_0_xxyyyz_xyzz, g_z_0_xxyyyz_xzzz, g_z_0_xxyyyz_yyyy, g_z_0_xxyyyz_yyyz, g_z_0_xxyyyz_yyzz, g_z_0_xxyyyz_yzzz, g_z_0_xxyyyz_zzzz, g_z_0_xyyyz_xxxx, g_z_0_xyyyz_xxxxx, g_z_0_xyyyz_xxxxy, g_z_0_xyyyz_xxxxz, g_z_0_xyyyz_xxxy, g_z_0_xyyyz_xxxyy, g_z_0_xyyyz_xxxyz, g_z_0_xyyyz_xxxz, g_z_0_xyyyz_xxxzz, g_z_0_xyyyz_xxyy, g_z_0_xyyyz_xxyyy, g_z_0_xyyyz_xxyyz, g_z_0_xyyyz_xxyz, g_z_0_xyyyz_xxyzz, g_z_0_xyyyz_xxzz, g_z_0_xyyyz_xxzzz, g_z_0_xyyyz_xyyy, g_z_0_xyyyz_xyyyy, g_z_0_xyyyz_xyyyz, g_z_0_xyyyz_xyyz, g_z_0_xyyyz_xyyzz, g_z_0_xyyyz_xyzz, g_z_0_xyyyz_xyzzz, g_z_0_xyyyz_xzzz, g_z_0_xyyyz_xzzzz, g_z_0_xyyyz_yyyy, g_z_0_xyyyz_yyyz, g_z_0_xyyyz_yyzz, g_z_0_xyyyz_yzzz, g_z_0_xyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyyz_xxxx[k] = -g_z_0_xyyyz_xxxx[k] * ab_x + g_z_0_xyyyz_xxxxx[k];

                g_z_0_xxyyyz_xxxy[k] = -g_z_0_xyyyz_xxxy[k] * ab_x + g_z_0_xyyyz_xxxxy[k];

                g_z_0_xxyyyz_xxxz[k] = -g_z_0_xyyyz_xxxz[k] * ab_x + g_z_0_xyyyz_xxxxz[k];

                g_z_0_xxyyyz_xxyy[k] = -g_z_0_xyyyz_xxyy[k] * ab_x + g_z_0_xyyyz_xxxyy[k];

                g_z_0_xxyyyz_xxyz[k] = -g_z_0_xyyyz_xxyz[k] * ab_x + g_z_0_xyyyz_xxxyz[k];

                g_z_0_xxyyyz_xxzz[k] = -g_z_0_xyyyz_xxzz[k] * ab_x + g_z_0_xyyyz_xxxzz[k];

                g_z_0_xxyyyz_xyyy[k] = -g_z_0_xyyyz_xyyy[k] * ab_x + g_z_0_xyyyz_xxyyy[k];

                g_z_0_xxyyyz_xyyz[k] = -g_z_0_xyyyz_xyyz[k] * ab_x + g_z_0_xyyyz_xxyyz[k];

                g_z_0_xxyyyz_xyzz[k] = -g_z_0_xyyyz_xyzz[k] * ab_x + g_z_0_xyyyz_xxyzz[k];

                g_z_0_xxyyyz_xzzz[k] = -g_z_0_xyyyz_xzzz[k] * ab_x + g_z_0_xyyyz_xxzzz[k];

                g_z_0_xxyyyz_yyyy[k] = -g_z_0_xyyyz_yyyy[k] * ab_x + g_z_0_xyyyz_xyyyy[k];

                g_z_0_xxyyyz_yyyz[k] = -g_z_0_xyyyz_yyyz[k] * ab_x + g_z_0_xyyyz_xyyyz[k];

                g_z_0_xxyyyz_yyzz[k] = -g_z_0_xyyyz_yyzz[k] * ab_x + g_z_0_xyyyz_xyyzz[k];

                g_z_0_xxyyyz_yzzz[k] = -g_z_0_xyyyz_yzzz[k] * ab_x + g_z_0_xyyyz_xyzzz[k];

                g_z_0_xxyyyz_zzzz[k] = -g_z_0_xyyyz_zzzz[k] * ab_x + g_z_0_xyyyz_xzzzz[k];
            }

            /// Set up 1020-1035 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyyzz_xxxx = cbuffer.data(ig_geom_10_off + 1020 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxy = cbuffer.data(ig_geom_10_off + 1021 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxxz = cbuffer.data(ig_geom_10_off + 1022 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxyy = cbuffer.data(ig_geom_10_off + 1023 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxyz = cbuffer.data(ig_geom_10_off + 1024 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xxzz = cbuffer.data(ig_geom_10_off + 1025 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyyy = cbuffer.data(ig_geom_10_off + 1026 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyyz = cbuffer.data(ig_geom_10_off + 1027 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xyzz = cbuffer.data(ig_geom_10_off + 1028 * ccomps * dcomps);

            auto g_z_0_xxyyzz_xzzz = cbuffer.data(ig_geom_10_off + 1029 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyyy = cbuffer.data(ig_geom_10_off + 1030 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyyz = cbuffer.data(ig_geom_10_off + 1031 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yyzz = cbuffer.data(ig_geom_10_off + 1032 * ccomps * dcomps);

            auto g_z_0_xxyyzz_yzzz = cbuffer.data(ig_geom_10_off + 1033 * ccomps * dcomps);

            auto g_z_0_xxyyzz_zzzz = cbuffer.data(ig_geom_10_off + 1034 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyyzz_xxxx, g_z_0_xxyyzz_xxxy, g_z_0_xxyyzz_xxxz, g_z_0_xxyyzz_xxyy, g_z_0_xxyyzz_xxyz, g_z_0_xxyyzz_xxzz, g_z_0_xxyyzz_xyyy, g_z_0_xxyyzz_xyyz, g_z_0_xxyyzz_xyzz, g_z_0_xxyyzz_xzzz, g_z_0_xxyyzz_yyyy, g_z_0_xxyyzz_yyyz, g_z_0_xxyyzz_yyzz, g_z_0_xxyyzz_yzzz, g_z_0_xxyyzz_zzzz, g_z_0_xyyzz_xxxx, g_z_0_xyyzz_xxxxx, g_z_0_xyyzz_xxxxy, g_z_0_xyyzz_xxxxz, g_z_0_xyyzz_xxxy, g_z_0_xyyzz_xxxyy, g_z_0_xyyzz_xxxyz, g_z_0_xyyzz_xxxz, g_z_0_xyyzz_xxxzz, g_z_0_xyyzz_xxyy, g_z_0_xyyzz_xxyyy, g_z_0_xyyzz_xxyyz, g_z_0_xyyzz_xxyz, g_z_0_xyyzz_xxyzz, g_z_0_xyyzz_xxzz, g_z_0_xyyzz_xxzzz, g_z_0_xyyzz_xyyy, g_z_0_xyyzz_xyyyy, g_z_0_xyyzz_xyyyz, g_z_0_xyyzz_xyyz, g_z_0_xyyzz_xyyzz, g_z_0_xyyzz_xyzz, g_z_0_xyyzz_xyzzz, g_z_0_xyyzz_xzzz, g_z_0_xyyzz_xzzzz, g_z_0_xyyzz_yyyy, g_z_0_xyyzz_yyyz, g_z_0_xyyzz_yyzz, g_z_0_xyyzz_yzzz, g_z_0_xyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyyzz_xxxx[k] = -g_z_0_xyyzz_xxxx[k] * ab_x + g_z_0_xyyzz_xxxxx[k];

                g_z_0_xxyyzz_xxxy[k] = -g_z_0_xyyzz_xxxy[k] * ab_x + g_z_0_xyyzz_xxxxy[k];

                g_z_0_xxyyzz_xxxz[k] = -g_z_0_xyyzz_xxxz[k] * ab_x + g_z_0_xyyzz_xxxxz[k];

                g_z_0_xxyyzz_xxyy[k] = -g_z_0_xyyzz_xxyy[k] * ab_x + g_z_0_xyyzz_xxxyy[k];

                g_z_0_xxyyzz_xxyz[k] = -g_z_0_xyyzz_xxyz[k] * ab_x + g_z_0_xyyzz_xxxyz[k];

                g_z_0_xxyyzz_xxzz[k] = -g_z_0_xyyzz_xxzz[k] * ab_x + g_z_0_xyyzz_xxxzz[k];

                g_z_0_xxyyzz_xyyy[k] = -g_z_0_xyyzz_xyyy[k] * ab_x + g_z_0_xyyzz_xxyyy[k];

                g_z_0_xxyyzz_xyyz[k] = -g_z_0_xyyzz_xyyz[k] * ab_x + g_z_0_xyyzz_xxyyz[k];

                g_z_0_xxyyzz_xyzz[k] = -g_z_0_xyyzz_xyzz[k] * ab_x + g_z_0_xyyzz_xxyzz[k];

                g_z_0_xxyyzz_xzzz[k] = -g_z_0_xyyzz_xzzz[k] * ab_x + g_z_0_xyyzz_xxzzz[k];

                g_z_0_xxyyzz_yyyy[k] = -g_z_0_xyyzz_yyyy[k] * ab_x + g_z_0_xyyzz_xyyyy[k];

                g_z_0_xxyyzz_yyyz[k] = -g_z_0_xyyzz_yyyz[k] * ab_x + g_z_0_xyyzz_xyyyz[k];

                g_z_0_xxyyzz_yyzz[k] = -g_z_0_xyyzz_yyzz[k] * ab_x + g_z_0_xyyzz_xyyzz[k];

                g_z_0_xxyyzz_yzzz[k] = -g_z_0_xyyzz_yzzz[k] * ab_x + g_z_0_xyyzz_xyzzz[k];

                g_z_0_xxyyzz_zzzz[k] = -g_z_0_xyyzz_zzzz[k] * ab_x + g_z_0_xyyzz_xzzzz[k];
            }

            /// Set up 1035-1050 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxyzzz_xxxx = cbuffer.data(ig_geom_10_off + 1035 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxy = cbuffer.data(ig_geom_10_off + 1036 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxxz = cbuffer.data(ig_geom_10_off + 1037 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxyy = cbuffer.data(ig_geom_10_off + 1038 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxyz = cbuffer.data(ig_geom_10_off + 1039 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xxzz = cbuffer.data(ig_geom_10_off + 1040 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyyy = cbuffer.data(ig_geom_10_off + 1041 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyyz = cbuffer.data(ig_geom_10_off + 1042 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xyzz = cbuffer.data(ig_geom_10_off + 1043 * ccomps * dcomps);

            auto g_z_0_xxyzzz_xzzz = cbuffer.data(ig_geom_10_off + 1044 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyyy = cbuffer.data(ig_geom_10_off + 1045 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyyz = cbuffer.data(ig_geom_10_off + 1046 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yyzz = cbuffer.data(ig_geom_10_off + 1047 * ccomps * dcomps);

            auto g_z_0_xxyzzz_yzzz = cbuffer.data(ig_geom_10_off + 1048 * ccomps * dcomps);

            auto g_z_0_xxyzzz_zzzz = cbuffer.data(ig_geom_10_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxyzzz_xxxx, g_z_0_xxyzzz_xxxy, g_z_0_xxyzzz_xxxz, g_z_0_xxyzzz_xxyy, g_z_0_xxyzzz_xxyz, g_z_0_xxyzzz_xxzz, g_z_0_xxyzzz_xyyy, g_z_0_xxyzzz_xyyz, g_z_0_xxyzzz_xyzz, g_z_0_xxyzzz_xzzz, g_z_0_xxyzzz_yyyy, g_z_0_xxyzzz_yyyz, g_z_0_xxyzzz_yyzz, g_z_0_xxyzzz_yzzz, g_z_0_xxyzzz_zzzz, g_z_0_xyzzz_xxxx, g_z_0_xyzzz_xxxxx, g_z_0_xyzzz_xxxxy, g_z_0_xyzzz_xxxxz, g_z_0_xyzzz_xxxy, g_z_0_xyzzz_xxxyy, g_z_0_xyzzz_xxxyz, g_z_0_xyzzz_xxxz, g_z_0_xyzzz_xxxzz, g_z_0_xyzzz_xxyy, g_z_0_xyzzz_xxyyy, g_z_0_xyzzz_xxyyz, g_z_0_xyzzz_xxyz, g_z_0_xyzzz_xxyzz, g_z_0_xyzzz_xxzz, g_z_0_xyzzz_xxzzz, g_z_0_xyzzz_xyyy, g_z_0_xyzzz_xyyyy, g_z_0_xyzzz_xyyyz, g_z_0_xyzzz_xyyz, g_z_0_xyzzz_xyyzz, g_z_0_xyzzz_xyzz, g_z_0_xyzzz_xyzzz, g_z_0_xyzzz_xzzz, g_z_0_xyzzz_xzzzz, g_z_0_xyzzz_yyyy, g_z_0_xyzzz_yyyz, g_z_0_xyzzz_yyzz, g_z_0_xyzzz_yzzz, g_z_0_xyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxyzzz_xxxx[k] = -g_z_0_xyzzz_xxxx[k] * ab_x + g_z_0_xyzzz_xxxxx[k];

                g_z_0_xxyzzz_xxxy[k] = -g_z_0_xyzzz_xxxy[k] * ab_x + g_z_0_xyzzz_xxxxy[k];

                g_z_0_xxyzzz_xxxz[k] = -g_z_0_xyzzz_xxxz[k] * ab_x + g_z_0_xyzzz_xxxxz[k];

                g_z_0_xxyzzz_xxyy[k] = -g_z_0_xyzzz_xxyy[k] * ab_x + g_z_0_xyzzz_xxxyy[k];

                g_z_0_xxyzzz_xxyz[k] = -g_z_0_xyzzz_xxyz[k] * ab_x + g_z_0_xyzzz_xxxyz[k];

                g_z_0_xxyzzz_xxzz[k] = -g_z_0_xyzzz_xxzz[k] * ab_x + g_z_0_xyzzz_xxxzz[k];

                g_z_0_xxyzzz_xyyy[k] = -g_z_0_xyzzz_xyyy[k] * ab_x + g_z_0_xyzzz_xxyyy[k];

                g_z_0_xxyzzz_xyyz[k] = -g_z_0_xyzzz_xyyz[k] * ab_x + g_z_0_xyzzz_xxyyz[k];

                g_z_0_xxyzzz_xyzz[k] = -g_z_0_xyzzz_xyzz[k] * ab_x + g_z_0_xyzzz_xxyzz[k];

                g_z_0_xxyzzz_xzzz[k] = -g_z_0_xyzzz_xzzz[k] * ab_x + g_z_0_xyzzz_xxzzz[k];

                g_z_0_xxyzzz_yyyy[k] = -g_z_0_xyzzz_yyyy[k] * ab_x + g_z_0_xyzzz_xyyyy[k];

                g_z_0_xxyzzz_yyyz[k] = -g_z_0_xyzzz_yyyz[k] * ab_x + g_z_0_xyzzz_xyyyz[k];

                g_z_0_xxyzzz_yyzz[k] = -g_z_0_xyzzz_yyzz[k] * ab_x + g_z_0_xyzzz_xyyzz[k];

                g_z_0_xxyzzz_yzzz[k] = -g_z_0_xyzzz_yzzz[k] * ab_x + g_z_0_xyzzz_xyzzz[k];

                g_z_0_xxyzzz_zzzz[k] = -g_z_0_xyzzz_zzzz[k] * ab_x + g_z_0_xyzzz_xzzzz[k];
            }

            /// Set up 1050-1065 components of targeted buffer : cbuffer.data(

            auto g_z_0_xxzzzz_xxxx = cbuffer.data(ig_geom_10_off + 1050 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxy = cbuffer.data(ig_geom_10_off + 1051 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxxz = cbuffer.data(ig_geom_10_off + 1052 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxyy = cbuffer.data(ig_geom_10_off + 1053 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxyz = cbuffer.data(ig_geom_10_off + 1054 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xxzz = cbuffer.data(ig_geom_10_off + 1055 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyyy = cbuffer.data(ig_geom_10_off + 1056 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyyz = cbuffer.data(ig_geom_10_off + 1057 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xyzz = cbuffer.data(ig_geom_10_off + 1058 * ccomps * dcomps);

            auto g_z_0_xxzzzz_xzzz = cbuffer.data(ig_geom_10_off + 1059 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyyy = cbuffer.data(ig_geom_10_off + 1060 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyyz = cbuffer.data(ig_geom_10_off + 1061 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yyzz = cbuffer.data(ig_geom_10_off + 1062 * ccomps * dcomps);

            auto g_z_0_xxzzzz_yzzz = cbuffer.data(ig_geom_10_off + 1063 * ccomps * dcomps);

            auto g_z_0_xxzzzz_zzzz = cbuffer.data(ig_geom_10_off + 1064 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xxzzzz_xxxx, g_z_0_xxzzzz_xxxy, g_z_0_xxzzzz_xxxz, g_z_0_xxzzzz_xxyy, g_z_0_xxzzzz_xxyz, g_z_0_xxzzzz_xxzz, g_z_0_xxzzzz_xyyy, g_z_0_xxzzzz_xyyz, g_z_0_xxzzzz_xyzz, g_z_0_xxzzzz_xzzz, g_z_0_xxzzzz_yyyy, g_z_0_xxzzzz_yyyz, g_z_0_xxzzzz_yyzz, g_z_0_xxzzzz_yzzz, g_z_0_xxzzzz_zzzz, g_z_0_xzzzz_xxxx, g_z_0_xzzzz_xxxxx, g_z_0_xzzzz_xxxxy, g_z_0_xzzzz_xxxxz, g_z_0_xzzzz_xxxy, g_z_0_xzzzz_xxxyy, g_z_0_xzzzz_xxxyz, g_z_0_xzzzz_xxxz, g_z_0_xzzzz_xxxzz, g_z_0_xzzzz_xxyy, g_z_0_xzzzz_xxyyy, g_z_0_xzzzz_xxyyz, g_z_0_xzzzz_xxyz, g_z_0_xzzzz_xxyzz, g_z_0_xzzzz_xxzz, g_z_0_xzzzz_xxzzz, g_z_0_xzzzz_xyyy, g_z_0_xzzzz_xyyyy, g_z_0_xzzzz_xyyyz, g_z_0_xzzzz_xyyz, g_z_0_xzzzz_xyyzz, g_z_0_xzzzz_xyzz, g_z_0_xzzzz_xyzzz, g_z_0_xzzzz_xzzz, g_z_0_xzzzz_xzzzz, g_z_0_xzzzz_yyyy, g_z_0_xzzzz_yyyz, g_z_0_xzzzz_yyzz, g_z_0_xzzzz_yzzz, g_z_0_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xxzzzz_xxxx[k] = -g_z_0_xzzzz_xxxx[k] * ab_x + g_z_0_xzzzz_xxxxx[k];

                g_z_0_xxzzzz_xxxy[k] = -g_z_0_xzzzz_xxxy[k] * ab_x + g_z_0_xzzzz_xxxxy[k];

                g_z_0_xxzzzz_xxxz[k] = -g_z_0_xzzzz_xxxz[k] * ab_x + g_z_0_xzzzz_xxxxz[k];

                g_z_0_xxzzzz_xxyy[k] = -g_z_0_xzzzz_xxyy[k] * ab_x + g_z_0_xzzzz_xxxyy[k];

                g_z_0_xxzzzz_xxyz[k] = -g_z_0_xzzzz_xxyz[k] * ab_x + g_z_0_xzzzz_xxxyz[k];

                g_z_0_xxzzzz_xxzz[k] = -g_z_0_xzzzz_xxzz[k] * ab_x + g_z_0_xzzzz_xxxzz[k];

                g_z_0_xxzzzz_xyyy[k] = -g_z_0_xzzzz_xyyy[k] * ab_x + g_z_0_xzzzz_xxyyy[k];

                g_z_0_xxzzzz_xyyz[k] = -g_z_0_xzzzz_xyyz[k] * ab_x + g_z_0_xzzzz_xxyyz[k];

                g_z_0_xxzzzz_xyzz[k] = -g_z_0_xzzzz_xyzz[k] * ab_x + g_z_0_xzzzz_xxyzz[k];

                g_z_0_xxzzzz_xzzz[k] = -g_z_0_xzzzz_xzzz[k] * ab_x + g_z_0_xzzzz_xxzzz[k];

                g_z_0_xxzzzz_yyyy[k] = -g_z_0_xzzzz_yyyy[k] * ab_x + g_z_0_xzzzz_xyyyy[k];

                g_z_0_xxzzzz_yyyz[k] = -g_z_0_xzzzz_yyyz[k] * ab_x + g_z_0_xzzzz_xyyyz[k];

                g_z_0_xxzzzz_yyzz[k] = -g_z_0_xzzzz_yyzz[k] * ab_x + g_z_0_xzzzz_xyyzz[k];

                g_z_0_xxzzzz_yzzz[k] = -g_z_0_xzzzz_yzzz[k] * ab_x + g_z_0_xzzzz_xyzzz[k];

                g_z_0_xxzzzz_zzzz[k] = -g_z_0_xzzzz_zzzz[k] * ab_x + g_z_0_xzzzz_xzzzz[k];
            }

            /// Set up 1065-1080 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 1065 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 1066 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 1067 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 1068 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 1069 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 1070 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 1071 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 1072 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 1073 * ccomps * dcomps);

            auto g_z_0_xyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 1074 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 1075 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 1076 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 1077 * ccomps * dcomps);

            auto g_z_0_xyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 1078 * ccomps * dcomps);

            auto g_z_0_xyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 1079 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyyy_xxxx, g_z_0_xyyyyy_xxxy, g_z_0_xyyyyy_xxxz, g_z_0_xyyyyy_xxyy, g_z_0_xyyyyy_xxyz, g_z_0_xyyyyy_xxzz, g_z_0_xyyyyy_xyyy, g_z_0_xyyyyy_xyyz, g_z_0_xyyyyy_xyzz, g_z_0_xyyyyy_xzzz, g_z_0_xyyyyy_yyyy, g_z_0_xyyyyy_yyyz, g_z_0_xyyyyy_yyzz, g_z_0_xyyyyy_yzzz, g_z_0_xyyyyy_zzzz, g_z_0_yyyyy_xxxx, g_z_0_yyyyy_xxxxx, g_z_0_yyyyy_xxxxy, g_z_0_yyyyy_xxxxz, g_z_0_yyyyy_xxxy, g_z_0_yyyyy_xxxyy, g_z_0_yyyyy_xxxyz, g_z_0_yyyyy_xxxz, g_z_0_yyyyy_xxxzz, g_z_0_yyyyy_xxyy, g_z_0_yyyyy_xxyyy, g_z_0_yyyyy_xxyyz, g_z_0_yyyyy_xxyz, g_z_0_yyyyy_xxyzz, g_z_0_yyyyy_xxzz, g_z_0_yyyyy_xxzzz, g_z_0_yyyyy_xyyy, g_z_0_yyyyy_xyyyy, g_z_0_yyyyy_xyyyz, g_z_0_yyyyy_xyyz, g_z_0_yyyyy_xyyzz, g_z_0_yyyyy_xyzz, g_z_0_yyyyy_xyzzz, g_z_0_yyyyy_xzzz, g_z_0_yyyyy_xzzzz, g_z_0_yyyyy_yyyy, g_z_0_yyyyy_yyyz, g_z_0_yyyyy_yyzz, g_z_0_yyyyy_yzzz, g_z_0_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyy_xxxx[k] = -g_z_0_yyyyy_xxxx[k] * ab_x + g_z_0_yyyyy_xxxxx[k];

                g_z_0_xyyyyy_xxxy[k] = -g_z_0_yyyyy_xxxy[k] * ab_x + g_z_0_yyyyy_xxxxy[k];

                g_z_0_xyyyyy_xxxz[k] = -g_z_0_yyyyy_xxxz[k] * ab_x + g_z_0_yyyyy_xxxxz[k];

                g_z_0_xyyyyy_xxyy[k] = -g_z_0_yyyyy_xxyy[k] * ab_x + g_z_0_yyyyy_xxxyy[k];

                g_z_0_xyyyyy_xxyz[k] = -g_z_0_yyyyy_xxyz[k] * ab_x + g_z_0_yyyyy_xxxyz[k];

                g_z_0_xyyyyy_xxzz[k] = -g_z_0_yyyyy_xxzz[k] * ab_x + g_z_0_yyyyy_xxxzz[k];

                g_z_0_xyyyyy_xyyy[k] = -g_z_0_yyyyy_xyyy[k] * ab_x + g_z_0_yyyyy_xxyyy[k];

                g_z_0_xyyyyy_xyyz[k] = -g_z_0_yyyyy_xyyz[k] * ab_x + g_z_0_yyyyy_xxyyz[k];

                g_z_0_xyyyyy_xyzz[k] = -g_z_0_yyyyy_xyzz[k] * ab_x + g_z_0_yyyyy_xxyzz[k];

                g_z_0_xyyyyy_xzzz[k] = -g_z_0_yyyyy_xzzz[k] * ab_x + g_z_0_yyyyy_xxzzz[k];

                g_z_0_xyyyyy_yyyy[k] = -g_z_0_yyyyy_yyyy[k] * ab_x + g_z_0_yyyyy_xyyyy[k];

                g_z_0_xyyyyy_yyyz[k] = -g_z_0_yyyyy_yyyz[k] * ab_x + g_z_0_yyyyy_xyyyz[k];

                g_z_0_xyyyyy_yyzz[k] = -g_z_0_yyyyy_yyzz[k] * ab_x + g_z_0_yyyyy_xyyzz[k];

                g_z_0_xyyyyy_yzzz[k] = -g_z_0_yyyyy_yzzz[k] * ab_x + g_z_0_yyyyy_xyzzz[k];

                g_z_0_xyyyyy_zzzz[k] = -g_z_0_yyyyy_zzzz[k] * ab_x + g_z_0_yyyyy_xzzzz[k];
            }

            /// Set up 1080-1095 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 1080 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 1081 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 1082 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 1083 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 1084 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 1085 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 1086 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 1087 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 1088 * ccomps * dcomps);

            auto g_z_0_xyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 1089 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 1090 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 1091 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 1092 * ccomps * dcomps);

            auto g_z_0_xyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 1093 * ccomps * dcomps);

            auto g_z_0_xyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 1094 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyyz_xxxx, g_z_0_xyyyyz_xxxy, g_z_0_xyyyyz_xxxz, g_z_0_xyyyyz_xxyy, g_z_0_xyyyyz_xxyz, g_z_0_xyyyyz_xxzz, g_z_0_xyyyyz_xyyy, g_z_0_xyyyyz_xyyz, g_z_0_xyyyyz_xyzz, g_z_0_xyyyyz_xzzz, g_z_0_xyyyyz_yyyy, g_z_0_xyyyyz_yyyz, g_z_0_xyyyyz_yyzz, g_z_0_xyyyyz_yzzz, g_z_0_xyyyyz_zzzz, g_z_0_yyyyz_xxxx, g_z_0_yyyyz_xxxxx, g_z_0_yyyyz_xxxxy, g_z_0_yyyyz_xxxxz, g_z_0_yyyyz_xxxy, g_z_0_yyyyz_xxxyy, g_z_0_yyyyz_xxxyz, g_z_0_yyyyz_xxxz, g_z_0_yyyyz_xxxzz, g_z_0_yyyyz_xxyy, g_z_0_yyyyz_xxyyy, g_z_0_yyyyz_xxyyz, g_z_0_yyyyz_xxyz, g_z_0_yyyyz_xxyzz, g_z_0_yyyyz_xxzz, g_z_0_yyyyz_xxzzz, g_z_0_yyyyz_xyyy, g_z_0_yyyyz_xyyyy, g_z_0_yyyyz_xyyyz, g_z_0_yyyyz_xyyz, g_z_0_yyyyz_xyyzz, g_z_0_yyyyz_xyzz, g_z_0_yyyyz_xyzzz, g_z_0_yyyyz_xzzz, g_z_0_yyyyz_xzzzz, g_z_0_yyyyz_yyyy, g_z_0_yyyyz_yyyz, g_z_0_yyyyz_yyzz, g_z_0_yyyyz_yzzz, g_z_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyyz_xxxx[k] = -g_z_0_yyyyz_xxxx[k] * ab_x + g_z_0_yyyyz_xxxxx[k];

                g_z_0_xyyyyz_xxxy[k] = -g_z_0_yyyyz_xxxy[k] * ab_x + g_z_0_yyyyz_xxxxy[k];

                g_z_0_xyyyyz_xxxz[k] = -g_z_0_yyyyz_xxxz[k] * ab_x + g_z_0_yyyyz_xxxxz[k];

                g_z_0_xyyyyz_xxyy[k] = -g_z_0_yyyyz_xxyy[k] * ab_x + g_z_0_yyyyz_xxxyy[k];

                g_z_0_xyyyyz_xxyz[k] = -g_z_0_yyyyz_xxyz[k] * ab_x + g_z_0_yyyyz_xxxyz[k];

                g_z_0_xyyyyz_xxzz[k] = -g_z_0_yyyyz_xxzz[k] * ab_x + g_z_0_yyyyz_xxxzz[k];

                g_z_0_xyyyyz_xyyy[k] = -g_z_0_yyyyz_xyyy[k] * ab_x + g_z_0_yyyyz_xxyyy[k];

                g_z_0_xyyyyz_xyyz[k] = -g_z_0_yyyyz_xyyz[k] * ab_x + g_z_0_yyyyz_xxyyz[k];

                g_z_0_xyyyyz_xyzz[k] = -g_z_0_yyyyz_xyzz[k] * ab_x + g_z_0_yyyyz_xxyzz[k];

                g_z_0_xyyyyz_xzzz[k] = -g_z_0_yyyyz_xzzz[k] * ab_x + g_z_0_yyyyz_xxzzz[k];

                g_z_0_xyyyyz_yyyy[k] = -g_z_0_yyyyz_yyyy[k] * ab_x + g_z_0_yyyyz_xyyyy[k];

                g_z_0_xyyyyz_yyyz[k] = -g_z_0_yyyyz_yyyz[k] * ab_x + g_z_0_yyyyz_xyyyz[k];

                g_z_0_xyyyyz_yyzz[k] = -g_z_0_yyyyz_yyzz[k] * ab_x + g_z_0_yyyyz_xyyzz[k];

                g_z_0_xyyyyz_yzzz[k] = -g_z_0_yyyyz_yzzz[k] * ab_x + g_z_0_yyyyz_xyzzz[k];

                g_z_0_xyyyyz_zzzz[k] = -g_z_0_yyyyz_zzzz[k] * ab_x + g_z_0_yyyyz_xzzzz[k];
            }

            /// Set up 1095-1110 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 1095 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 1096 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 1097 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 1098 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 1099 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 1100 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 1101 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 1102 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 1103 * ccomps * dcomps);

            auto g_z_0_xyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 1104 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 1105 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 1106 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 1107 * ccomps * dcomps);

            auto g_z_0_xyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 1108 * ccomps * dcomps);

            auto g_z_0_xyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 1109 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyyzz_xxxx, g_z_0_xyyyzz_xxxy, g_z_0_xyyyzz_xxxz, g_z_0_xyyyzz_xxyy, g_z_0_xyyyzz_xxyz, g_z_0_xyyyzz_xxzz, g_z_0_xyyyzz_xyyy, g_z_0_xyyyzz_xyyz, g_z_0_xyyyzz_xyzz, g_z_0_xyyyzz_xzzz, g_z_0_xyyyzz_yyyy, g_z_0_xyyyzz_yyyz, g_z_0_xyyyzz_yyzz, g_z_0_xyyyzz_yzzz, g_z_0_xyyyzz_zzzz, g_z_0_yyyzz_xxxx, g_z_0_yyyzz_xxxxx, g_z_0_yyyzz_xxxxy, g_z_0_yyyzz_xxxxz, g_z_0_yyyzz_xxxy, g_z_0_yyyzz_xxxyy, g_z_0_yyyzz_xxxyz, g_z_0_yyyzz_xxxz, g_z_0_yyyzz_xxxzz, g_z_0_yyyzz_xxyy, g_z_0_yyyzz_xxyyy, g_z_0_yyyzz_xxyyz, g_z_0_yyyzz_xxyz, g_z_0_yyyzz_xxyzz, g_z_0_yyyzz_xxzz, g_z_0_yyyzz_xxzzz, g_z_0_yyyzz_xyyy, g_z_0_yyyzz_xyyyy, g_z_0_yyyzz_xyyyz, g_z_0_yyyzz_xyyz, g_z_0_yyyzz_xyyzz, g_z_0_yyyzz_xyzz, g_z_0_yyyzz_xyzzz, g_z_0_yyyzz_xzzz, g_z_0_yyyzz_xzzzz, g_z_0_yyyzz_yyyy, g_z_0_yyyzz_yyyz, g_z_0_yyyzz_yyzz, g_z_0_yyyzz_yzzz, g_z_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyyzz_xxxx[k] = -g_z_0_yyyzz_xxxx[k] * ab_x + g_z_0_yyyzz_xxxxx[k];

                g_z_0_xyyyzz_xxxy[k] = -g_z_0_yyyzz_xxxy[k] * ab_x + g_z_0_yyyzz_xxxxy[k];

                g_z_0_xyyyzz_xxxz[k] = -g_z_0_yyyzz_xxxz[k] * ab_x + g_z_0_yyyzz_xxxxz[k];

                g_z_0_xyyyzz_xxyy[k] = -g_z_0_yyyzz_xxyy[k] * ab_x + g_z_0_yyyzz_xxxyy[k];

                g_z_0_xyyyzz_xxyz[k] = -g_z_0_yyyzz_xxyz[k] * ab_x + g_z_0_yyyzz_xxxyz[k];

                g_z_0_xyyyzz_xxzz[k] = -g_z_0_yyyzz_xxzz[k] * ab_x + g_z_0_yyyzz_xxxzz[k];

                g_z_0_xyyyzz_xyyy[k] = -g_z_0_yyyzz_xyyy[k] * ab_x + g_z_0_yyyzz_xxyyy[k];

                g_z_0_xyyyzz_xyyz[k] = -g_z_0_yyyzz_xyyz[k] * ab_x + g_z_0_yyyzz_xxyyz[k];

                g_z_0_xyyyzz_xyzz[k] = -g_z_0_yyyzz_xyzz[k] * ab_x + g_z_0_yyyzz_xxyzz[k];

                g_z_0_xyyyzz_xzzz[k] = -g_z_0_yyyzz_xzzz[k] * ab_x + g_z_0_yyyzz_xxzzz[k];

                g_z_0_xyyyzz_yyyy[k] = -g_z_0_yyyzz_yyyy[k] * ab_x + g_z_0_yyyzz_xyyyy[k];

                g_z_0_xyyyzz_yyyz[k] = -g_z_0_yyyzz_yyyz[k] * ab_x + g_z_0_yyyzz_xyyyz[k];

                g_z_0_xyyyzz_yyzz[k] = -g_z_0_yyyzz_yyzz[k] * ab_x + g_z_0_yyyzz_xyyzz[k];

                g_z_0_xyyyzz_yzzz[k] = -g_z_0_yyyzz_yzzz[k] * ab_x + g_z_0_yyyzz_xyzzz[k];

                g_z_0_xyyyzz_zzzz[k] = -g_z_0_yyyzz_zzzz[k] * ab_x + g_z_0_yyyzz_xzzzz[k];
            }

            /// Set up 1110-1125 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 1110 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 1111 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 1112 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 1113 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 1114 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 1115 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 1116 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 1117 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 1118 * ccomps * dcomps);

            auto g_z_0_xyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 1119 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 1120 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 1121 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 1122 * ccomps * dcomps);

            auto g_z_0_xyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 1123 * ccomps * dcomps);

            auto g_z_0_xyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 1124 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyyzzz_xxxx, g_z_0_xyyzzz_xxxy, g_z_0_xyyzzz_xxxz, g_z_0_xyyzzz_xxyy, g_z_0_xyyzzz_xxyz, g_z_0_xyyzzz_xxzz, g_z_0_xyyzzz_xyyy, g_z_0_xyyzzz_xyyz, g_z_0_xyyzzz_xyzz, g_z_0_xyyzzz_xzzz, g_z_0_xyyzzz_yyyy, g_z_0_xyyzzz_yyyz, g_z_0_xyyzzz_yyzz, g_z_0_xyyzzz_yzzz, g_z_0_xyyzzz_zzzz, g_z_0_yyzzz_xxxx, g_z_0_yyzzz_xxxxx, g_z_0_yyzzz_xxxxy, g_z_0_yyzzz_xxxxz, g_z_0_yyzzz_xxxy, g_z_0_yyzzz_xxxyy, g_z_0_yyzzz_xxxyz, g_z_0_yyzzz_xxxz, g_z_0_yyzzz_xxxzz, g_z_0_yyzzz_xxyy, g_z_0_yyzzz_xxyyy, g_z_0_yyzzz_xxyyz, g_z_0_yyzzz_xxyz, g_z_0_yyzzz_xxyzz, g_z_0_yyzzz_xxzz, g_z_0_yyzzz_xxzzz, g_z_0_yyzzz_xyyy, g_z_0_yyzzz_xyyyy, g_z_0_yyzzz_xyyyz, g_z_0_yyzzz_xyyz, g_z_0_yyzzz_xyyzz, g_z_0_yyzzz_xyzz, g_z_0_yyzzz_xyzzz, g_z_0_yyzzz_xzzz, g_z_0_yyzzz_xzzzz, g_z_0_yyzzz_yyyy, g_z_0_yyzzz_yyyz, g_z_0_yyzzz_yyzz, g_z_0_yyzzz_yzzz, g_z_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyyzzz_xxxx[k] = -g_z_0_yyzzz_xxxx[k] * ab_x + g_z_0_yyzzz_xxxxx[k];

                g_z_0_xyyzzz_xxxy[k] = -g_z_0_yyzzz_xxxy[k] * ab_x + g_z_0_yyzzz_xxxxy[k];

                g_z_0_xyyzzz_xxxz[k] = -g_z_0_yyzzz_xxxz[k] * ab_x + g_z_0_yyzzz_xxxxz[k];

                g_z_0_xyyzzz_xxyy[k] = -g_z_0_yyzzz_xxyy[k] * ab_x + g_z_0_yyzzz_xxxyy[k];

                g_z_0_xyyzzz_xxyz[k] = -g_z_0_yyzzz_xxyz[k] * ab_x + g_z_0_yyzzz_xxxyz[k];

                g_z_0_xyyzzz_xxzz[k] = -g_z_0_yyzzz_xxzz[k] * ab_x + g_z_0_yyzzz_xxxzz[k];

                g_z_0_xyyzzz_xyyy[k] = -g_z_0_yyzzz_xyyy[k] * ab_x + g_z_0_yyzzz_xxyyy[k];

                g_z_0_xyyzzz_xyyz[k] = -g_z_0_yyzzz_xyyz[k] * ab_x + g_z_0_yyzzz_xxyyz[k];

                g_z_0_xyyzzz_xyzz[k] = -g_z_0_yyzzz_xyzz[k] * ab_x + g_z_0_yyzzz_xxyzz[k];

                g_z_0_xyyzzz_xzzz[k] = -g_z_0_yyzzz_xzzz[k] * ab_x + g_z_0_yyzzz_xxzzz[k];

                g_z_0_xyyzzz_yyyy[k] = -g_z_0_yyzzz_yyyy[k] * ab_x + g_z_0_yyzzz_xyyyy[k];

                g_z_0_xyyzzz_yyyz[k] = -g_z_0_yyzzz_yyyz[k] * ab_x + g_z_0_yyzzz_xyyyz[k];

                g_z_0_xyyzzz_yyzz[k] = -g_z_0_yyzzz_yyzz[k] * ab_x + g_z_0_yyzzz_xyyzz[k];

                g_z_0_xyyzzz_yzzz[k] = -g_z_0_yyzzz_yzzz[k] * ab_x + g_z_0_yyzzz_xyzzz[k];

                g_z_0_xyyzzz_zzzz[k] = -g_z_0_yyzzz_zzzz[k] * ab_x + g_z_0_yyzzz_xzzzz[k];
            }

            /// Set up 1125-1140 components of targeted buffer : cbuffer.data(

            auto g_z_0_xyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 1125 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 1126 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 1127 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 1128 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 1129 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 1130 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 1131 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 1132 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 1133 * ccomps * dcomps);

            auto g_z_0_xyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 1134 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 1135 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 1136 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 1137 * ccomps * dcomps);

            auto g_z_0_xyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 1138 * ccomps * dcomps);

            auto g_z_0_xyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 1139 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xyzzzz_xxxx, g_z_0_xyzzzz_xxxy, g_z_0_xyzzzz_xxxz, g_z_0_xyzzzz_xxyy, g_z_0_xyzzzz_xxyz, g_z_0_xyzzzz_xxzz, g_z_0_xyzzzz_xyyy, g_z_0_xyzzzz_xyyz, g_z_0_xyzzzz_xyzz, g_z_0_xyzzzz_xzzz, g_z_0_xyzzzz_yyyy, g_z_0_xyzzzz_yyyz, g_z_0_xyzzzz_yyzz, g_z_0_xyzzzz_yzzz, g_z_0_xyzzzz_zzzz, g_z_0_yzzzz_xxxx, g_z_0_yzzzz_xxxxx, g_z_0_yzzzz_xxxxy, g_z_0_yzzzz_xxxxz, g_z_0_yzzzz_xxxy, g_z_0_yzzzz_xxxyy, g_z_0_yzzzz_xxxyz, g_z_0_yzzzz_xxxz, g_z_0_yzzzz_xxxzz, g_z_0_yzzzz_xxyy, g_z_0_yzzzz_xxyyy, g_z_0_yzzzz_xxyyz, g_z_0_yzzzz_xxyz, g_z_0_yzzzz_xxyzz, g_z_0_yzzzz_xxzz, g_z_0_yzzzz_xxzzz, g_z_0_yzzzz_xyyy, g_z_0_yzzzz_xyyyy, g_z_0_yzzzz_xyyyz, g_z_0_yzzzz_xyyz, g_z_0_yzzzz_xyyzz, g_z_0_yzzzz_xyzz, g_z_0_yzzzz_xyzzz, g_z_0_yzzzz_xzzz, g_z_0_yzzzz_xzzzz, g_z_0_yzzzz_yyyy, g_z_0_yzzzz_yyyz, g_z_0_yzzzz_yyzz, g_z_0_yzzzz_yzzz, g_z_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xyzzzz_xxxx[k] = -g_z_0_yzzzz_xxxx[k] * ab_x + g_z_0_yzzzz_xxxxx[k];

                g_z_0_xyzzzz_xxxy[k] = -g_z_0_yzzzz_xxxy[k] * ab_x + g_z_0_yzzzz_xxxxy[k];

                g_z_0_xyzzzz_xxxz[k] = -g_z_0_yzzzz_xxxz[k] * ab_x + g_z_0_yzzzz_xxxxz[k];

                g_z_0_xyzzzz_xxyy[k] = -g_z_0_yzzzz_xxyy[k] * ab_x + g_z_0_yzzzz_xxxyy[k];

                g_z_0_xyzzzz_xxyz[k] = -g_z_0_yzzzz_xxyz[k] * ab_x + g_z_0_yzzzz_xxxyz[k];

                g_z_0_xyzzzz_xxzz[k] = -g_z_0_yzzzz_xxzz[k] * ab_x + g_z_0_yzzzz_xxxzz[k];

                g_z_0_xyzzzz_xyyy[k] = -g_z_0_yzzzz_xyyy[k] * ab_x + g_z_0_yzzzz_xxyyy[k];

                g_z_0_xyzzzz_xyyz[k] = -g_z_0_yzzzz_xyyz[k] * ab_x + g_z_0_yzzzz_xxyyz[k];

                g_z_0_xyzzzz_xyzz[k] = -g_z_0_yzzzz_xyzz[k] * ab_x + g_z_0_yzzzz_xxyzz[k];

                g_z_0_xyzzzz_xzzz[k] = -g_z_0_yzzzz_xzzz[k] * ab_x + g_z_0_yzzzz_xxzzz[k];

                g_z_0_xyzzzz_yyyy[k] = -g_z_0_yzzzz_yyyy[k] * ab_x + g_z_0_yzzzz_xyyyy[k];

                g_z_0_xyzzzz_yyyz[k] = -g_z_0_yzzzz_yyyz[k] * ab_x + g_z_0_yzzzz_xyyyz[k];

                g_z_0_xyzzzz_yyzz[k] = -g_z_0_yzzzz_yyzz[k] * ab_x + g_z_0_yzzzz_xyyzz[k];

                g_z_0_xyzzzz_yzzz[k] = -g_z_0_yzzzz_yzzz[k] * ab_x + g_z_0_yzzzz_xyzzz[k];

                g_z_0_xyzzzz_zzzz[k] = -g_z_0_yzzzz_zzzz[k] * ab_x + g_z_0_yzzzz_xzzzz[k];
            }

            /// Set up 1140-1155 components of targeted buffer : cbuffer.data(

            auto g_z_0_xzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 1140 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 1141 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 1142 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 1143 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 1144 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 1145 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 1146 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 1147 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 1148 * ccomps * dcomps);

            auto g_z_0_xzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 1149 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 1150 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 1151 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 1152 * ccomps * dcomps);

            auto g_z_0_xzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 1153 * ccomps * dcomps);

            auto g_z_0_xzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 1154 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_xzzzzz_xxxx, g_z_0_xzzzzz_xxxy, g_z_0_xzzzzz_xxxz, g_z_0_xzzzzz_xxyy, g_z_0_xzzzzz_xxyz, g_z_0_xzzzzz_xxzz, g_z_0_xzzzzz_xyyy, g_z_0_xzzzzz_xyyz, g_z_0_xzzzzz_xyzz, g_z_0_xzzzzz_xzzz, g_z_0_xzzzzz_yyyy, g_z_0_xzzzzz_yyyz, g_z_0_xzzzzz_yyzz, g_z_0_xzzzzz_yzzz, g_z_0_xzzzzz_zzzz, g_z_0_zzzzz_xxxx, g_z_0_zzzzz_xxxxx, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_yyyy, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_xzzzzz_xxxx[k] = -g_z_0_zzzzz_xxxx[k] * ab_x + g_z_0_zzzzz_xxxxx[k];

                g_z_0_xzzzzz_xxxy[k] = -g_z_0_zzzzz_xxxy[k] * ab_x + g_z_0_zzzzz_xxxxy[k];

                g_z_0_xzzzzz_xxxz[k] = -g_z_0_zzzzz_xxxz[k] * ab_x + g_z_0_zzzzz_xxxxz[k];

                g_z_0_xzzzzz_xxyy[k] = -g_z_0_zzzzz_xxyy[k] * ab_x + g_z_0_zzzzz_xxxyy[k];

                g_z_0_xzzzzz_xxyz[k] = -g_z_0_zzzzz_xxyz[k] * ab_x + g_z_0_zzzzz_xxxyz[k];

                g_z_0_xzzzzz_xxzz[k] = -g_z_0_zzzzz_xxzz[k] * ab_x + g_z_0_zzzzz_xxxzz[k];

                g_z_0_xzzzzz_xyyy[k] = -g_z_0_zzzzz_xyyy[k] * ab_x + g_z_0_zzzzz_xxyyy[k];

                g_z_0_xzzzzz_xyyz[k] = -g_z_0_zzzzz_xyyz[k] * ab_x + g_z_0_zzzzz_xxyyz[k];

                g_z_0_xzzzzz_xyzz[k] = -g_z_0_zzzzz_xyzz[k] * ab_x + g_z_0_zzzzz_xxyzz[k];

                g_z_0_xzzzzz_xzzz[k] = -g_z_0_zzzzz_xzzz[k] * ab_x + g_z_0_zzzzz_xxzzz[k];

                g_z_0_xzzzzz_yyyy[k] = -g_z_0_zzzzz_yyyy[k] * ab_x + g_z_0_zzzzz_xyyyy[k];

                g_z_0_xzzzzz_yyyz[k] = -g_z_0_zzzzz_yyyz[k] * ab_x + g_z_0_zzzzz_xyyyz[k];

                g_z_0_xzzzzz_yyzz[k] = -g_z_0_zzzzz_yyzz[k] * ab_x + g_z_0_zzzzz_xyyzz[k];

                g_z_0_xzzzzz_yzzz[k] = -g_z_0_zzzzz_yzzz[k] * ab_x + g_z_0_zzzzz_xyzzz[k];

                g_z_0_xzzzzz_zzzz[k] = -g_z_0_zzzzz_zzzz[k] * ab_x + g_z_0_zzzzz_xzzzz[k];
            }

            /// Set up 1155-1170 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyy_xxxx = cbuffer.data(ig_geom_10_off + 1155 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxy = cbuffer.data(ig_geom_10_off + 1156 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxxz = cbuffer.data(ig_geom_10_off + 1157 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxyy = cbuffer.data(ig_geom_10_off + 1158 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxyz = cbuffer.data(ig_geom_10_off + 1159 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xxzz = cbuffer.data(ig_geom_10_off + 1160 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyyy = cbuffer.data(ig_geom_10_off + 1161 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyyz = cbuffer.data(ig_geom_10_off + 1162 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xyzz = cbuffer.data(ig_geom_10_off + 1163 * ccomps * dcomps);

            auto g_z_0_yyyyyy_xzzz = cbuffer.data(ig_geom_10_off + 1164 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyyy = cbuffer.data(ig_geom_10_off + 1165 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyyz = cbuffer.data(ig_geom_10_off + 1166 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yyzz = cbuffer.data(ig_geom_10_off + 1167 * ccomps * dcomps);

            auto g_z_0_yyyyyy_yzzz = cbuffer.data(ig_geom_10_off + 1168 * ccomps * dcomps);

            auto g_z_0_yyyyyy_zzzz = cbuffer.data(ig_geom_10_off + 1169 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyy_xxxx, g_z_0_yyyyy_xxxxy, g_z_0_yyyyy_xxxy, g_z_0_yyyyy_xxxyy, g_z_0_yyyyy_xxxyz, g_z_0_yyyyy_xxxz, g_z_0_yyyyy_xxyy, g_z_0_yyyyy_xxyyy, g_z_0_yyyyy_xxyyz, g_z_0_yyyyy_xxyz, g_z_0_yyyyy_xxyzz, g_z_0_yyyyy_xxzz, g_z_0_yyyyy_xyyy, g_z_0_yyyyy_xyyyy, g_z_0_yyyyy_xyyyz, g_z_0_yyyyy_xyyz, g_z_0_yyyyy_xyyzz, g_z_0_yyyyy_xyzz, g_z_0_yyyyy_xyzzz, g_z_0_yyyyy_xzzz, g_z_0_yyyyy_yyyy, g_z_0_yyyyy_yyyyy, g_z_0_yyyyy_yyyyz, g_z_0_yyyyy_yyyz, g_z_0_yyyyy_yyyzz, g_z_0_yyyyy_yyzz, g_z_0_yyyyy_yyzzz, g_z_0_yyyyy_yzzz, g_z_0_yyyyy_yzzzz, g_z_0_yyyyy_zzzz, g_z_0_yyyyyy_xxxx, g_z_0_yyyyyy_xxxy, g_z_0_yyyyyy_xxxz, g_z_0_yyyyyy_xxyy, g_z_0_yyyyyy_xxyz, g_z_0_yyyyyy_xxzz, g_z_0_yyyyyy_xyyy, g_z_0_yyyyyy_xyyz, g_z_0_yyyyyy_xyzz, g_z_0_yyyyyy_xzzz, g_z_0_yyyyyy_yyyy, g_z_0_yyyyyy_yyyz, g_z_0_yyyyyy_yyzz, g_z_0_yyyyyy_yzzz, g_z_0_yyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyy_xxxx[k] = -g_z_0_yyyyy_xxxx[k] * ab_y + g_z_0_yyyyy_xxxxy[k];

                g_z_0_yyyyyy_xxxy[k] = -g_z_0_yyyyy_xxxy[k] * ab_y + g_z_0_yyyyy_xxxyy[k];

                g_z_0_yyyyyy_xxxz[k] = -g_z_0_yyyyy_xxxz[k] * ab_y + g_z_0_yyyyy_xxxyz[k];

                g_z_0_yyyyyy_xxyy[k] = -g_z_0_yyyyy_xxyy[k] * ab_y + g_z_0_yyyyy_xxyyy[k];

                g_z_0_yyyyyy_xxyz[k] = -g_z_0_yyyyy_xxyz[k] * ab_y + g_z_0_yyyyy_xxyyz[k];

                g_z_0_yyyyyy_xxzz[k] = -g_z_0_yyyyy_xxzz[k] * ab_y + g_z_0_yyyyy_xxyzz[k];

                g_z_0_yyyyyy_xyyy[k] = -g_z_0_yyyyy_xyyy[k] * ab_y + g_z_0_yyyyy_xyyyy[k];

                g_z_0_yyyyyy_xyyz[k] = -g_z_0_yyyyy_xyyz[k] * ab_y + g_z_0_yyyyy_xyyyz[k];

                g_z_0_yyyyyy_xyzz[k] = -g_z_0_yyyyy_xyzz[k] * ab_y + g_z_0_yyyyy_xyyzz[k];

                g_z_0_yyyyyy_xzzz[k] = -g_z_0_yyyyy_xzzz[k] * ab_y + g_z_0_yyyyy_xyzzz[k];

                g_z_0_yyyyyy_yyyy[k] = -g_z_0_yyyyy_yyyy[k] * ab_y + g_z_0_yyyyy_yyyyy[k];

                g_z_0_yyyyyy_yyyz[k] = -g_z_0_yyyyy_yyyz[k] * ab_y + g_z_0_yyyyy_yyyyz[k];

                g_z_0_yyyyyy_yyzz[k] = -g_z_0_yyyyy_yyzz[k] * ab_y + g_z_0_yyyyy_yyyzz[k];

                g_z_0_yyyyyy_yzzz[k] = -g_z_0_yyyyy_yzzz[k] * ab_y + g_z_0_yyyyy_yyzzz[k];

                g_z_0_yyyyyy_zzzz[k] = -g_z_0_yyyyy_zzzz[k] * ab_y + g_z_0_yyyyy_yzzzz[k];
            }

            /// Set up 1170-1185 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyyz_xxxx = cbuffer.data(ig_geom_10_off + 1170 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxy = cbuffer.data(ig_geom_10_off + 1171 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxxz = cbuffer.data(ig_geom_10_off + 1172 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxyy = cbuffer.data(ig_geom_10_off + 1173 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxyz = cbuffer.data(ig_geom_10_off + 1174 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xxzz = cbuffer.data(ig_geom_10_off + 1175 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyyy = cbuffer.data(ig_geom_10_off + 1176 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyyz = cbuffer.data(ig_geom_10_off + 1177 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xyzz = cbuffer.data(ig_geom_10_off + 1178 * ccomps * dcomps);

            auto g_z_0_yyyyyz_xzzz = cbuffer.data(ig_geom_10_off + 1179 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyyy = cbuffer.data(ig_geom_10_off + 1180 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyyz = cbuffer.data(ig_geom_10_off + 1181 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yyzz = cbuffer.data(ig_geom_10_off + 1182 * ccomps * dcomps);

            auto g_z_0_yyyyyz_yzzz = cbuffer.data(ig_geom_10_off + 1183 * ccomps * dcomps);

            auto g_z_0_yyyyyz_zzzz = cbuffer.data(ig_geom_10_off + 1184 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyyz_xxxx, g_z_0_yyyyyz_xxxy, g_z_0_yyyyyz_xxxz, g_z_0_yyyyyz_xxyy, g_z_0_yyyyyz_xxyz, g_z_0_yyyyyz_xxzz, g_z_0_yyyyyz_xyyy, g_z_0_yyyyyz_xyyz, g_z_0_yyyyyz_xyzz, g_z_0_yyyyyz_xzzz, g_z_0_yyyyyz_yyyy, g_z_0_yyyyyz_yyyz, g_z_0_yyyyyz_yyzz, g_z_0_yyyyyz_yzzz, g_z_0_yyyyyz_zzzz, g_z_0_yyyyz_xxxx, g_z_0_yyyyz_xxxxy, g_z_0_yyyyz_xxxy, g_z_0_yyyyz_xxxyy, g_z_0_yyyyz_xxxyz, g_z_0_yyyyz_xxxz, g_z_0_yyyyz_xxyy, g_z_0_yyyyz_xxyyy, g_z_0_yyyyz_xxyyz, g_z_0_yyyyz_xxyz, g_z_0_yyyyz_xxyzz, g_z_0_yyyyz_xxzz, g_z_0_yyyyz_xyyy, g_z_0_yyyyz_xyyyy, g_z_0_yyyyz_xyyyz, g_z_0_yyyyz_xyyz, g_z_0_yyyyz_xyyzz, g_z_0_yyyyz_xyzz, g_z_0_yyyyz_xyzzz, g_z_0_yyyyz_xzzz, g_z_0_yyyyz_yyyy, g_z_0_yyyyz_yyyyy, g_z_0_yyyyz_yyyyz, g_z_0_yyyyz_yyyz, g_z_0_yyyyz_yyyzz, g_z_0_yyyyz_yyzz, g_z_0_yyyyz_yyzzz, g_z_0_yyyyz_yzzz, g_z_0_yyyyz_yzzzz, g_z_0_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyyz_xxxx[k] = -g_z_0_yyyyz_xxxx[k] * ab_y + g_z_0_yyyyz_xxxxy[k];

                g_z_0_yyyyyz_xxxy[k] = -g_z_0_yyyyz_xxxy[k] * ab_y + g_z_0_yyyyz_xxxyy[k];

                g_z_0_yyyyyz_xxxz[k] = -g_z_0_yyyyz_xxxz[k] * ab_y + g_z_0_yyyyz_xxxyz[k];

                g_z_0_yyyyyz_xxyy[k] = -g_z_0_yyyyz_xxyy[k] * ab_y + g_z_0_yyyyz_xxyyy[k];

                g_z_0_yyyyyz_xxyz[k] = -g_z_0_yyyyz_xxyz[k] * ab_y + g_z_0_yyyyz_xxyyz[k];

                g_z_0_yyyyyz_xxzz[k] = -g_z_0_yyyyz_xxzz[k] * ab_y + g_z_0_yyyyz_xxyzz[k];

                g_z_0_yyyyyz_xyyy[k] = -g_z_0_yyyyz_xyyy[k] * ab_y + g_z_0_yyyyz_xyyyy[k];

                g_z_0_yyyyyz_xyyz[k] = -g_z_0_yyyyz_xyyz[k] * ab_y + g_z_0_yyyyz_xyyyz[k];

                g_z_0_yyyyyz_xyzz[k] = -g_z_0_yyyyz_xyzz[k] * ab_y + g_z_0_yyyyz_xyyzz[k];

                g_z_0_yyyyyz_xzzz[k] = -g_z_0_yyyyz_xzzz[k] * ab_y + g_z_0_yyyyz_xyzzz[k];

                g_z_0_yyyyyz_yyyy[k] = -g_z_0_yyyyz_yyyy[k] * ab_y + g_z_0_yyyyz_yyyyy[k];

                g_z_0_yyyyyz_yyyz[k] = -g_z_0_yyyyz_yyyz[k] * ab_y + g_z_0_yyyyz_yyyyz[k];

                g_z_0_yyyyyz_yyzz[k] = -g_z_0_yyyyz_yyzz[k] * ab_y + g_z_0_yyyyz_yyyzz[k];

                g_z_0_yyyyyz_yzzz[k] = -g_z_0_yyyyz_yzzz[k] * ab_y + g_z_0_yyyyz_yyzzz[k];

                g_z_0_yyyyyz_zzzz[k] = -g_z_0_yyyyz_zzzz[k] * ab_y + g_z_0_yyyyz_yzzzz[k];
            }

            /// Set up 1185-1200 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyyzz_xxxx = cbuffer.data(ig_geom_10_off + 1185 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxy = cbuffer.data(ig_geom_10_off + 1186 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxxz = cbuffer.data(ig_geom_10_off + 1187 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxyy = cbuffer.data(ig_geom_10_off + 1188 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxyz = cbuffer.data(ig_geom_10_off + 1189 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xxzz = cbuffer.data(ig_geom_10_off + 1190 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyyy = cbuffer.data(ig_geom_10_off + 1191 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyyz = cbuffer.data(ig_geom_10_off + 1192 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xyzz = cbuffer.data(ig_geom_10_off + 1193 * ccomps * dcomps);

            auto g_z_0_yyyyzz_xzzz = cbuffer.data(ig_geom_10_off + 1194 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyyy = cbuffer.data(ig_geom_10_off + 1195 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyyz = cbuffer.data(ig_geom_10_off + 1196 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yyzz = cbuffer.data(ig_geom_10_off + 1197 * ccomps * dcomps);

            auto g_z_0_yyyyzz_yzzz = cbuffer.data(ig_geom_10_off + 1198 * ccomps * dcomps);

            auto g_z_0_yyyyzz_zzzz = cbuffer.data(ig_geom_10_off + 1199 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyyzz_xxxx, g_z_0_yyyyzz_xxxy, g_z_0_yyyyzz_xxxz, g_z_0_yyyyzz_xxyy, g_z_0_yyyyzz_xxyz, g_z_0_yyyyzz_xxzz, g_z_0_yyyyzz_xyyy, g_z_0_yyyyzz_xyyz, g_z_0_yyyyzz_xyzz, g_z_0_yyyyzz_xzzz, g_z_0_yyyyzz_yyyy, g_z_0_yyyyzz_yyyz, g_z_0_yyyyzz_yyzz, g_z_0_yyyyzz_yzzz, g_z_0_yyyyzz_zzzz, g_z_0_yyyzz_xxxx, g_z_0_yyyzz_xxxxy, g_z_0_yyyzz_xxxy, g_z_0_yyyzz_xxxyy, g_z_0_yyyzz_xxxyz, g_z_0_yyyzz_xxxz, g_z_0_yyyzz_xxyy, g_z_0_yyyzz_xxyyy, g_z_0_yyyzz_xxyyz, g_z_0_yyyzz_xxyz, g_z_0_yyyzz_xxyzz, g_z_0_yyyzz_xxzz, g_z_0_yyyzz_xyyy, g_z_0_yyyzz_xyyyy, g_z_0_yyyzz_xyyyz, g_z_0_yyyzz_xyyz, g_z_0_yyyzz_xyyzz, g_z_0_yyyzz_xyzz, g_z_0_yyyzz_xyzzz, g_z_0_yyyzz_xzzz, g_z_0_yyyzz_yyyy, g_z_0_yyyzz_yyyyy, g_z_0_yyyzz_yyyyz, g_z_0_yyyzz_yyyz, g_z_0_yyyzz_yyyzz, g_z_0_yyyzz_yyzz, g_z_0_yyyzz_yyzzz, g_z_0_yyyzz_yzzz, g_z_0_yyyzz_yzzzz, g_z_0_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyyzz_xxxx[k] = -g_z_0_yyyzz_xxxx[k] * ab_y + g_z_0_yyyzz_xxxxy[k];

                g_z_0_yyyyzz_xxxy[k] = -g_z_0_yyyzz_xxxy[k] * ab_y + g_z_0_yyyzz_xxxyy[k];

                g_z_0_yyyyzz_xxxz[k] = -g_z_0_yyyzz_xxxz[k] * ab_y + g_z_0_yyyzz_xxxyz[k];

                g_z_0_yyyyzz_xxyy[k] = -g_z_0_yyyzz_xxyy[k] * ab_y + g_z_0_yyyzz_xxyyy[k];

                g_z_0_yyyyzz_xxyz[k] = -g_z_0_yyyzz_xxyz[k] * ab_y + g_z_0_yyyzz_xxyyz[k];

                g_z_0_yyyyzz_xxzz[k] = -g_z_0_yyyzz_xxzz[k] * ab_y + g_z_0_yyyzz_xxyzz[k];

                g_z_0_yyyyzz_xyyy[k] = -g_z_0_yyyzz_xyyy[k] * ab_y + g_z_0_yyyzz_xyyyy[k];

                g_z_0_yyyyzz_xyyz[k] = -g_z_0_yyyzz_xyyz[k] * ab_y + g_z_0_yyyzz_xyyyz[k];

                g_z_0_yyyyzz_xyzz[k] = -g_z_0_yyyzz_xyzz[k] * ab_y + g_z_0_yyyzz_xyyzz[k];

                g_z_0_yyyyzz_xzzz[k] = -g_z_0_yyyzz_xzzz[k] * ab_y + g_z_0_yyyzz_xyzzz[k];

                g_z_0_yyyyzz_yyyy[k] = -g_z_0_yyyzz_yyyy[k] * ab_y + g_z_0_yyyzz_yyyyy[k];

                g_z_0_yyyyzz_yyyz[k] = -g_z_0_yyyzz_yyyz[k] * ab_y + g_z_0_yyyzz_yyyyz[k];

                g_z_0_yyyyzz_yyzz[k] = -g_z_0_yyyzz_yyzz[k] * ab_y + g_z_0_yyyzz_yyyzz[k];

                g_z_0_yyyyzz_yzzz[k] = -g_z_0_yyyzz_yzzz[k] * ab_y + g_z_0_yyyzz_yyzzz[k];

                g_z_0_yyyyzz_zzzz[k] = -g_z_0_yyyzz_zzzz[k] * ab_y + g_z_0_yyyzz_yzzzz[k];
            }

            /// Set up 1200-1215 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyyzzz_xxxx = cbuffer.data(ig_geom_10_off + 1200 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxy = cbuffer.data(ig_geom_10_off + 1201 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxxz = cbuffer.data(ig_geom_10_off + 1202 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxyy = cbuffer.data(ig_geom_10_off + 1203 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxyz = cbuffer.data(ig_geom_10_off + 1204 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xxzz = cbuffer.data(ig_geom_10_off + 1205 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyyy = cbuffer.data(ig_geom_10_off + 1206 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyyz = cbuffer.data(ig_geom_10_off + 1207 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xyzz = cbuffer.data(ig_geom_10_off + 1208 * ccomps * dcomps);

            auto g_z_0_yyyzzz_xzzz = cbuffer.data(ig_geom_10_off + 1209 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyyy = cbuffer.data(ig_geom_10_off + 1210 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyyz = cbuffer.data(ig_geom_10_off + 1211 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yyzz = cbuffer.data(ig_geom_10_off + 1212 * ccomps * dcomps);

            auto g_z_0_yyyzzz_yzzz = cbuffer.data(ig_geom_10_off + 1213 * ccomps * dcomps);

            auto g_z_0_yyyzzz_zzzz = cbuffer.data(ig_geom_10_off + 1214 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyyzzz_xxxx, g_z_0_yyyzzz_xxxy, g_z_0_yyyzzz_xxxz, g_z_0_yyyzzz_xxyy, g_z_0_yyyzzz_xxyz, g_z_0_yyyzzz_xxzz, g_z_0_yyyzzz_xyyy, g_z_0_yyyzzz_xyyz, g_z_0_yyyzzz_xyzz, g_z_0_yyyzzz_xzzz, g_z_0_yyyzzz_yyyy, g_z_0_yyyzzz_yyyz, g_z_0_yyyzzz_yyzz, g_z_0_yyyzzz_yzzz, g_z_0_yyyzzz_zzzz, g_z_0_yyzzz_xxxx, g_z_0_yyzzz_xxxxy, g_z_0_yyzzz_xxxy, g_z_0_yyzzz_xxxyy, g_z_0_yyzzz_xxxyz, g_z_0_yyzzz_xxxz, g_z_0_yyzzz_xxyy, g_z_0_yyzzz_xxyyy, g_z_0_yyzzz_xxyyz, g_z_0_yyzzz_xxyz, g_z_0_yyzzz_xxyzz, g_z_0_yyzzz_xxzz, g_z_0_yyzzz_xyyy, g_z_0_yyzzz_xyyyy, g_z_0_yyzzz_xyyyz, g_z_0_yyzzz_xyyz, g_z_0_yyzzz_xyyzz, g_z_0_yyzzz_xyzz, g_z_0_yyzzz_xyzzz, g_z_0_yyzzz_xzzz, g_z_0_yyzzz_yyyy, g_z_0_yyzzz_yyyyy, g_z_0_yyzzz_yyyyz, g_z_0_yyzzz_yyyz, g_z_0_yyzzz_yyyzz, g_z_0_yyzzz_yyzz, g_z_0_yyzzz_yyzzz, g_z_0_yyzzz_yzzz, g_z_0_yyzzz_yzzzz, g_z_0_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyyzzz_xxxx[k] = -g_z_0_yyzzz_xxxx[k] * ab_y + g_z_0_yyzzz_xxxxy[k];

                g_z_0_yyyzzz_xxxy[k] = -g_z_0_yyzzz_xxxy[k] * ab_y + g_z_0_yyzzz_xxxyy[k];

                g_z_0_yyyzzz_xxxz[k] = -g_z_0_yyzzz_xxxz[k] * ab_y + g_z_0_yyzzz_xxxyz[k];

                g_z_0_yyyzzz_xxyy[k] = -g_z_0_yyzzz_xxyy[k] * ab_y + g_z_0_yyzzz_xxyyy[k];

                g_z_0_yyyzzz_xxyz[k] = -g_z_0_yyzzz_xxyz[k] * ab_y + g_z_0_yyzzz_xxyyz[k];

                g_z_0_yyyzzz_xxzz[k] = -g_z_0_yyzzz_xxzz[k] * ab_y + g_z_0_yyzzz_xxyzz[k];

                g_z_0_yyyzzz_xyyy[k] = -g_z_0_yyzzz_xyyy[k] * ab_y + g_z_0_yyzzz_xyyyy[k];

                g_z_0_yyyzzz_xyyz[k] = -g_z_0_yyzzz_xyyz[k] * ab_y + g_z_0_yyzzz_xyyyz[k];

                g_z_0_yyyzzz_xyzz[k] = -g_z_0_yyzzz_xyzz[k] * ab_y + g_z_0_yyzzz_xyyzz[k];

                g_z_0_yyyzzz_xzzz[k] = -g_z_0_yyzzz_xzzz[k] * ab_y + g_z_0_yyzzz_xyzzz[k];

                g_z_0_yyyzzz_yyyy[k] = -g_z_0_yyzzz_yyyy[k] * ab_y + g_z_0_yyzzz_yyyyy[k];

                g_z_0_yyyzzz_yyyz[k] = -g_z_0_yyzzz_yyyz[k] * ab_y + g_z_0_yyzzz_yyyyz[k];

                g_z_0_yyyzzz_yyzz[k] = -g_z_0_yyzzz_yyzz[k] * ab_y + g_z_0_yyzzz_yyyzz[k];

                g_z_0_yyyzzz_yzzz[k] = -g_z_0_yyzzz_yzzz[k] * ab_y + g_z_0_yyzzz_yyzzz[k];

                g_z_0_yyyzzz_zzzz[k] = -g_z_0_yyzzz_zzzz[k] * ab_y + g_z_0_yyzzz_yzzzz[k];
            }

            /// Set up 1215-1230 components of targeted buffer : cbuffer.data(

            auto g_z_0_yyzzzz_xxxx = cbuffer.data(ig_geom_10_off + 1215 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxy = cbuffer.data(ig_geom_10_off + 1216 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxxz = cbuffer.data(ig_geom_10_off + 1217 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxyy = cbuffer.data(ig_geom_10_off + 1218 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxyz = cbuffer.data(ig_geom_10_off + 1219 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xxzz = cbuffer.data(ig_geom_10_off + 1220 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyyy = cbuffer.data(ig_geom_10_off + 1221 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyyz = cbuffer.data(ig_geom_10_off + 1222 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xyzz = cbuffer.data(ig_geom_10_off + 1223 * ccomps * dcomps);

            auto g_z_0_yyzzzz_xzzz = cbuffer.data(ig_geom_10_off + 1224 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyyy = cbuffer.data(ig_geom_10_off + 1225 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyyz = cbuffer.data(ig_geom_10_off + 1226 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yyzz = cbuffer.data(ig_geom_10_off + 1227 * ccomps * dcomps);

            auto g_z_0_yyzzzz_yzzz = cbuffer.data(ig_geom_10_off + 1228 * ccomps * dcomps);

            auto g_z_0_yyzzzz_zzzz = cbuffer.data(ig_geom_10_off + 1229 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yyzzzz_xxxx, g_z_0_yyzzzz_xxxy, g_z_0_yyzzzz_xxxz, g_z_0_yyzzzz_xxyy, g_z_0_yyzzzz_xxyz, g_z_0_yyzzzz_xxzz, g_z_0_yyzzzz_xyyy, g_z_0_yyzzzz_xyyz, g_z_0_yyzzzz_xyzz, g_z_0_yyzzzz_xzzz, g_z_0_yyzzzz_yyyy, g_z_0_yyzzzz_yyyz, g_z_0_yyzzzz_yyzz, g_z_0_yyzzzz_yzzz, g_z_0_yyzzzz_zzzz, g_z_0_yzzzz_xxxx, g_z_0_yzzzz_xxxxy, g_z_0_yzzzz_xxxy, g_z_0_yzzzz_xxxyy, g_z_0_yzzzz_xxxyz, g_z_0_yzzzz_xxxz, g_z_0_yzzzz_xxyy, g_z_0_yzzzz_xxyyy, g_z_0_yzzzz_xxyyz, g_z_0_yzzzz_xxyz, g_z_0_yzzzz_xxyzz, g_z_0_yzzzz_xxzz, g_z_0_yzzzz_xyyy, g_z_0_yzzzz_xyyyy, g_z_0_yzzzz_xyyyz, g_z_0_yzzzz_xyyz, g_z_0_yzzzz_xyyzz, g_z_0_yzzzz_xyzz, g_z_0_yzzzz_xyzzz, g_z_0_yzzzz_xzzz, g_z_0_yzzzz_yyyy, g_z_0_yzzzz_yyyyy, g_z_0_yzzzz_yyyyz, g_z_0_yzzzz_yyyz, g_z_0_yzzzz_yyyzz, g_z_0_yzzzz_yyzz, g_z_0_yzzzz_yyzzz, g_z_0_yzzzz_yzzz, g_z_0_yzzzz_yzzzz, g_z_0_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yyzzzz_xxxx[k] = -g_z_0_yzzzz_xxxx[k] * ab_y + g_z_0_yzzzz_xxxxy[k];

                g_z_0_yyzzzz_xxxy[k] = -g_z_0_yzzzz_xxxy[k] * ab_y + g_z_0_yzzzz_xxxyy[k];

                g_z_0_yyzzzz_xxxz[k] = -g_z_0_yzzzz_xxxz[k] * ab_y + g_z_0_yzzzz_xxxyz[k];

                g_z_0_yyzzzz_xxyy[k] = -g_z_0_yzzzz_xxyy[k] * ab_y + g_z_0_yzzzz_xxyyy[k];

                g_z_0_yyzzzz_xxyz[k] = -g_z_0_yzzzz_xxyz[k] * ab_y + g_z_0_yzzzz_xxyyz[k];

                g_z_0_yyzzzz_xxzz[k] = -g_z_0_yzzzz_xxzz[k] * ab_y + g_z_0_yzzzz_xxyzz[k];

                g_z_0_yyzzzz_xyyy[k] = -g_z_0_yzzzz_xyyy[k] * ab_y + g_z_0_yzzzz_xyyyy[k];

                g_z_0_yyzzzz_xyyz[k] = -g_z_0_yzzzz_xyyz[k] * ab_y + g_z_0_yzzzz_xyyyz[k];

                g_z_0_yyzzzz_xyzz[k] = -g_z_0_yzzzz_xyzz[k] * ab_y + g_z_0_yzzzz_xyyzz[k];

                g_z_0_yyzzzz_xzzz[k] = -g_z_0_yzzzz_xzzz[k] * ab_y + g_z_0_yzzzz_xyzzz[k];

                g_z_0_yyzzzz_yyyy[k] = -g_z_0_yzzzz_yyyy[k] * ab_y + g_z_0_yzzzz_yyyyy[k];

                g_z_0_yyzzzz_yyyz[k] = -g_z_0_yzzzz_yyyz[k] * ab_y + g_z_0_yzzzz_yyyyz[k];

                g_z_0_yyzzzz_yyzz[k] = -g_z_0_yzzzz_yyzz[k] * ab_y + g_z_0_yzzzz_yyyzz[k];

                g_z_0_yyzzzz_yzzz[k] = -g_z_0_yzzzz_yzzz[k] * ab_y + g_z_0_yzzzz_yyzzz[k];

                g_z_0_yyzzzz_zzzz[k] = -g_z_0_yzzzz_zzzz[k] * ab_y + g_z_0_yzzzz_yzzzz[k];
            }

            /// Set up 1230-1245 components of targeted buffer : cbuffer.data(

            auto g_z_0_yzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 1230 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 1231 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 1232 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 1233 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 1234 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 1235 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 1236 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 1237 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 1238 * ccomps * dcomps);

            auto g_z_0_yzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 1239 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 1240 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 1241 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 1242 * ccomps * dcomps);

            auto g_z_0_yzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 1243 * ccomps * dcomps);

            auto g_z_0_yzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 1244 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_yzzzzz_xxxx, g_z_0_yzzzzz_xxxy, g_z_0_yzzzzz_xxxz, g_z_0_yzzzzz_xxyy, g_z_0_yzzzzz_xxyz, g_z_0_yzzzzz_xxzz, g_z_0_yzzzzz_xyyy, g_z_0_yzzzzz_xyyz, g_z_0_yzzzzz_xyzz, g_z_0_yzzzzz_xzzz, g_z_0_yzzzzz_yyyy, g_z_0_yzzzzz_yyyz, g_z_0_yzzzzz_yyzz, g_z_0_yzzzzz_yzzz, g_z_0_yzzzzz_zzzz, g_z_0_zzzzz_xxxx, g_z_0_zzzzz_xxxxy, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxxyy, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyyy, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyyy, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_yyyy, g_z_0_zzzzz_yyyyy, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_yzzzzz_xxxx[k] = -g_z_0_zzzzz_xxxx[k] * ab_y + g_z_0_zzzzz_xxxxy[k];

                g_z_0_yzzzzz_xxxy[k] = -g_z_0_zzzzz_xxxy[k] * ab_y + g_z_0_zzzzz_xxxyy[k];

                g_z_0_yzzzzz_xxxz[k] = -g_z_0_zzzzz_xxxz[k] * ab_y + g_z_0_zzzzz_xxxyz[k];

                g_z_0_yzzzzz_xxyy[k] = -g_z_0_zzzzz_xxyy[k] * ab_y + g_z_0_zzzzz_xxyyy[k];

                g_z_0_yzzzzz_xxyz[k] = -g_z_0_zzzzz_xxyz[k] * ab_y + g_z_0_zzzzz_xxyyz[k];

                g_z_0_yzzzzz_xxzz[k] = -g_z_0_zzzzz_xxzz[k] * ab_y + g_z_0_zzzzz_xxyzz[k];

                g_z_0_yzzzzz_xyyy[k] = -g_z_0_zzzzz_xyyy[k] * ab_y + g_z_0_zzzzz_xyyyy[k];

                g_z_0_yzzzzz_xyyz[k] = -g_z_0_zzzzz_xyyz[k] * ab_y + g_z_0_zzzzz_xyyyz[k];

                g_z_0_yzzzzz_xyzz[k] = -g_z_0_zzzzz_xyzz[k] * ab_y + g_z_0_zzzzz_xyyzz[k];

                g_z_0_yzzzzz_xzzz[k] = -g_z_0_zzzzz_xzzz[k] * ab_y + g_z_0_zzzzz_xyzzz[k];

                g_z_0_yzzzzz_yyyy[k] = -g_z_0_zzzzz_yyyy[k] * ab_y + g_z_0_zzzzz_yyyyy[k];

                g_z_0_yzzzzz_yyyz[k] = -g_z_0_zzzzz_yyyz[k] * ab_y + g_z_0_zzzzz_yyyyz[k];

                g_z_0_yzzzzz_yyzz[k] = -g_z_0_zzzzz_yyzz[k] * ab_y + g_z_0_zzzzz_yyyzz[k];

                g_z_0_yzzzzz_yzzz[k] = -g_z_0_zzzzz_yzzz[k] * ab_y + g_z_0_zzzzz_yyzzz[k];

                g_z_0_yzzzzz_zzzz[k] = -g_z_0_zzzzz_zzzz[k] * ab_y + g_z_0_zzzzz_yzzzz[k];
            }

            /// Set up 1245-1260 components of targeted buffer : cbuffer.data(

            auto g_z_0_zzzzzz_xxxx = cbuffer.data(ig_geom_10_off + 1245 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxy = cbuffer.data(ig_geom_10_off + 1246 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxxz = cbuffer.data(ig_geom_10_off + 1247 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxyy = cbuffer.data(ig_geom_10_off + 1248 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxyz = cbuffer.data(ig_geom_10_off + 1249 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xxzz = cbuffer.data(ig_geom_10_off + 1250 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyyy = cbuffer.data(ig_geom_10_off + 1251 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyyz = cbuffer.data(ig_geom_10_off + 1252 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xyzz = cbuffer.data(ig_geom_10_off + 1253 * ccomps * dcomps);

            auto g_z_0_zzzzzz_xzzz = cbuffer.data(ig_geom_10_off + 1254 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyyy = cbuffer.data(ig_geom_10_off + 1255 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyyz = cbuffer.data(ig_geom_10_off + 1256 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yyzz = cbuffer.data(ig_geom_10_off + 1257 * ccomps * dcomps);

            auto g_z_0_zzzzzz_yzzz = cbuffer.data(ig_geom_10_off + 1258 * ccomps * dcomps);

            auto g_z_0_zzzzzz_zzzz = cbuffer.data(ig_geom_10_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_z_0_zzzzz_xxxx, g_z_0_zzzzz_xxxxz, g_z_0_zzzzz_xxxy, g_z_0_zzzzz_xxxyz, g_z_0_zzzzz_xxxz, g_z_0_zzzzz_xxxzz, g_z_0_zzzzz_xxyy, g_z_0_zzzzz_xxyyz, g_z_0_zzzzz_xxyz, g_z_0_zzzzz_xxyzz, g_z_0_zzzzz_xxzz, g_z_0_zzzzz_xxzzz, g_z_0_zzzzz_xyyy, g_z_0_zzzzz_xyyyz, g_z_0_zzzzz_xyyz, g_z_0_zzzzz_xyyzz, g_z_0_zzzzz_xyzz, g_z_0_zzzzz_xyzzz, g_z_0_zzzzz_xzzz, g_z_0_zzzzz_xzzzz, g_z_0_zzzzz_yyyy, g_z_0_zzzzz_yyyyz, g_z_0_zzzzz_yyyz, g_z_0_zzzzz_yyyzz, g_z_0_zzzzz_yyzz, g_z_0_zzzzz_yyzzz, g_z_0_zzzzz_yzzz, g_z_0_zzzzz_yzzzz, g_z_0_zzzzz_zzzz, g_z_0_zzzzz_zzzzz, g_z_0_zzzzzz_xxxx, g_z_0_zzzzzz_xxxy, g_z_0_zzzzzz_xxxz, g_z_0_zzzzzz_xxyy, g_z_0_zzzzzz_xxyz, g_z_0_zzzzzz_xxzz, g_z_0_zzzzzz_xyyy, g_z_0_zzzzzz_xyyz, g_z_0_zzzzzz_xyzz, g_z_0_zzzzzz_xzzz, g_z_0_zzzzzz_yyyy, g_z_0_zzzzzz_yyyz, g_z_0_zzzzzz_yyzz, g_z_0_zzzzzz_yzzz, g_z_0_zzzzzz_zzzz, g_zzzzz_xxxx, g_zzzzz_xxxy, g_zzzzz_xxxz, g_zzzzz_xxyy, g_zzzzz_xxyz, g_zzzzz_xxzz, g_zzzzz_xyyy, g_zzzzz_xyyz, g_zzzzz_xyzz, g_zzzzz_xzzz, g_zzzzz_yyyy, g_zzzzz_yyyz, g_zzzzz_yyzz, g_zzzzz_yzzz, g_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_z_0_zzzzzz_xxxx[k] = -g_zzzzz_xxxx[k] - g_z_0_zzzzz_xxxx[k] * ab_z + g_z_0_zzzzz_xxxxz[k];

                g_z_0_zzzzzz_xxxy[k] = -g_zzzzz_xxxy[k] - g_z_0_zzzzz_xxxy[k] * ab_z + g_z_0_zzzzz_xxxyz[k];

                g_z_0_zzzzzz_xxxz[k] = -g_zzzzz_xxxz[k] - g_z_0_zzzzz_xxxz[k] * ab_z + g_z_0_zzzzz_xxxzz[k];

                g_z_0_zzzzzz_xxyy[k] = -g_zzzzz_xxyy[k] - g_z_0_zzzzz_xxyy[k] * ab_z + g_z_0_zzzzz_xxyyz[k];

                g_z_0_zzzzzz_xxyz[k] = -g_zzzzz_xxyz[k] - g_z_0_zzzzz_xxyz[k] * ab_z + g_z_0_zzzzz_xxyzz[k];

                g_z_0_zzzzzz_xxzz[k] = -g_zzzzz_xxzz[k] - g_z_0_zzzzz_xxzz[k] * ab_z + g_z_0_zzzzz_xxzzz[k];

                g_z_0_zzzzzz_xyyy[k] = -g_zzzzz_xyyy[k] - g_z_0_zzzzz_xyyy[k] * ab_z + g_z_0_zzzzz_xyyyz[k];

                g_z_0_zzzzzz_xyyz[k] = -g_zzzzz_xyyz[k] - g_z_0_zzzzz_xyyz[k] * ab_z + g_z_0_zzzzz_xyyzz[k];

                g_z_0_zzzzzz_xyzz[k] = -g_zzzzz_xyzz[k] - g_z_0_zzzzz_xyzz[k] * ab_z + g_z_0_zzzzz_xyzzz[k];

                g_z_0_zzzzzz_xzzz[k] = -g_zzzzz_xzzz[k] - g_z_0_zzzzz_xzzz[k] * ab_z + g_z_0_zzzzz_xzzzz[k];

                g_z_0_zzzzzz_yyyy[k] = -g_zzzzz_yyyy[k] - g_z_0_zzzzz_yyyy[k] * ab_z + g_z_0_zzzzz_yyyyz[k];

                g_z_0_zzzzzz_yyyz[k] = -g_zzzzz_yyyz[k] - g_z_0_zzzzz_yyyz[k] * ab_z + g_z_0_zzzzz_yyyzz[k];

                g_z_0_zzzzzz_yyzz[k] = -g_zzzzz_yyzz[k] - g_z_0_zzzzz_yyzz[k] * ab_z + g_z_0_zzzzz_yyzzz[k];

                g_z_0_zzzzzz_yzzz[k] = -g_zzzzz_yzzz[k] - g_z_0_zzzzz_yzzz[k] * ab_z + g_z_0_zzzzz_yzzzz[k];

                g_z_0_zzzzzz_zzzz[k] = -g_zzzzz_zzzz[k] - g_z_0_zzzzz_zzzz[k] * ab_z + g_z_0_zzzzz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

