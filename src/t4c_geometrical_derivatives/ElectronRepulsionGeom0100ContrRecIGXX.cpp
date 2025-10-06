#include "ElectronRepulsionGeom0100ContrRecIGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_igxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_igxx,
                                            const size_t idx_hgxx,
                                            const size_t idx_geom_01_hgxx,
                                            const size_t idx_geom_01_hhxx,
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

            const auto hg_geom_01_off = idx_geom_01_hgxx + i * dcomps + j;

            auto g_0_x_xxxxx_xxxx = cbuffer.data(hg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxy = cbuffer.data(hg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxxz = cbuffer.data(hg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxyy = cbuffer.data(hg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxyz = cbuffer.data(hg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxx_xxzz = cbuffer.data(hg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyyy = cbuffer.data(hg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyyz = cbuffer.data(hg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxx_xyzz = cbuffer.data(hg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxx_xzzz = cbuffer.data(hg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyyy = cbuffer.data(hg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyyz = cbuffer.data(hg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxx_yyzz = cbuffer.data(hg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxx_yzzz = cbuffer.data(hg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxx_zzzz = cbuffer.data(hg_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxx = cbuffer.data(hg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxy = cbuffer.data(hg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxxz = cbuffer.data(hg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxyy = cbuffer.data(hg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxyz = cbuffer.data(hg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxy_xxzz = cbuffer.data(hg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyyy = cbuffer.data(hg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyyz = cbuffer.data(hg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxy_xyzz = cbuffer.data(hg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxy_xzzz = cbuffer.data(hg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyyy = cbuffer.data(hg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyyz = cbuffer.data(hg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxy_yyzz = cbuffer.data(hg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxy_yzzz = cbuffer.data(hg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxy_zzzz = cbuffer.data(hg_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxx = cbuffer.data(hg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxy = cbuffer.data(hg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxxz = cbuffer.data(hg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxyy = cbuffer.data(hg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxyz = cbuffer.data(hg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxz_xxzz = cbuffer.data(hg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyyy = cbuffer.data(hg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyyz = cbuffer.data(hg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxxz_xyzz = cbuffer.data(hg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxxz_xzzz = cbuffer.data(hg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyyy = cbuffer.data(hg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyyz = cbuffer.data(hg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxxz_yyzz = cbuffer.data(hg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxxz_yzzz = cbuffer.data(hg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxxz_zzzz = cbuffer.data(hg_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxx = cbuffer.data(hg_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxy = cbuffer.data(hg_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxxz = cbuffer.data(hg_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxyy = cbuffer.data(hg_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxyz = cbuffer.data(hg_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxyy_xxzz = cbuffer.data(hg_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyyy = cbuffer.data(hg_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyyz = cbuffer.data(hg_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxyy_xyzz = cbuffer.data(hg_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxyy_xzzz = cbuffer.data(hg_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyyy = cbuffer.data(hg_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyyz = cbuffer.data(hg_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxyy_yyzz = cbuffer.data(hg_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxyy_yzzz = cbuffer.data(hg_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxyy_zzzz = cbuffer.data(hg_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxx = cbuffer.data(hg_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxy = cbuffer.data(hg_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxxz = cbuffer.data(hg_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxyy = cbuffer.data(hg_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxyz = cbuffer.data(hg_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxyz_xxzz = cbuffer.data(hg_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyyy = cbuffer.data(hg_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyyz = cbuffer.data(hg_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxyz_xyzz = cbuffer.data(hg_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxyz_xzzz = cbuffer.data(hg_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyyy = cbuffer.data(hg_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyyz = cbuffer.data(hg_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxyz_yyzz = cbuffer.data(hg_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxyz_yzzz = cbuffer.data(hg_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxyz_zzzz = cbuffer.data(hg_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxx = cbuffer.data(hg_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxy = cbuffer.data(hg_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxxz = cbuffer.data(hg_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxyy = cbuffer.data(hg_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxyz = cbuffer.data(hg_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxzz_xxzz = cbuffer.data(hg_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyyy = cbuffer.data(hg_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyyz = cbuffer.data(hg_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxzz_xyzz = cbuffer.data(hg_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxxzz_xzzz = cbuffer.data(hg_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyyy = cbuffer.data(hg_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyyz = cbuffer.data(hg_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxzz_yyzz = cbuffer.data(hg_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxzz_yzzz = cbuffer.data(hg_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxzz_zzzz = cbuffer.data(hg_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxx = cbuffer.data(hg_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxy = cbuffer.data(hg_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxxz = cbuffer.data(hg_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxyy = cbuffer.data(hg_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxyz = cbuffer.data(hg_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxyyy_xxzz = cbuffer.data(hg_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyyy = cbuffer.data(hg_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyyz = cbuffer.data(hg_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxyyy_xyzz = cbuffer.data(hg_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxyyy_xzzz = cbuffer.data(hg_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyyy = cbuffer.data(hg_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyyz = cbuffer.data(hg_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxyyy_yyzz = cbuffer.data(hg_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxyyy_yzzz = cbuffer.data(hg_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxyyy_zzzz = cbuffer.data(hg_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxx = cbuffer.data(hg_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxy = cbuffer.data(hg_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxxz = cbuffer.data(hg_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxyy = cbuffer.data(hg_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxyz = cbuffer.data(hg_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxyyz_xxzz = cbuffer.data(hg_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyyy = cbuffer.data(hg_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyyz = cbuffer.data(hg_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxyyz_xyzz = cbuffer.data(hg_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxyyz_xzzz = cbuffer.data(hg_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyyy = cbuffer.data(hg_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyyz = cbuffer.data(hg_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxyyz_yyzz = cbuffer.data(hg_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxyyz_yzzz = cbuffer.data(hg_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxyyz_zzzz = cbuffer.data(hg_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxx = cbuffer.data(hg_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxy = cbuffer.data(hg_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxxz = cbuffer.data(hg_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxyy = cbuffer.data(hg_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxyz = cbuffer.data(hg_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxyzz_xxzz = cbuffer.data(hg_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyyy = cbuffer.data(hg_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyyz = cbuffer.data(hg_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxyzz_xyzz = cbuffer.data(hg_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxyzz_xzzz = cbuffer.data(hg_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyyy = cbuffer.data(hg_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyyz = cbuffer.data(hg_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxyzz_yyzz = cbuffer.data(hg_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxyzz_yzzz = cbuffer.data(hg_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxyzz_zzzz = cbuffer.data(hg_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxx = cbuffer.data(hg_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxy = cbuffer.data(hg_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxxz = cbuffer.data(hg_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxyy = cbuffer.data(hg_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxyz = cbuffer.data(hg_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxzzz_xxzz = cbuffer.data(hg_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyyy = cbuffer.data(hg_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyyz = cbuffer.data(hg_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxzzz_xyzz = cbuffer.data(hg_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxzzz_xzzz = cbuffer.data(hg_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyyy = cbuffer.data(hg_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyyz = cbuffer.data(hg_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxzzz_yyzz = cbuffer.data(hg_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxzzz_yzzz = cbuffer.data(hg_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxzzz_zzzz = cbuffer.data(hg_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxx = cbuffer.data(hg_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxy = cbuffer.data(hg_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxxz = cbuffer.data(hg_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxyy = cbuffer.data(hg_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxyz = cbuffer.data(hg_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xyyyy_xxzz = cbuffer.data(hg_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyyy = cbuffer.data(hg_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyyz = cbuffer.data(hg_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xyyyy_xyzz = cbuffer.data(hg_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xyyyy_xzzz = cbuffer.data(hg_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyyy = cbuffer.data(hg_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyyz = cbuffer.data(hg_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xyyyy_yyzz = cbuffer.data(hg_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xyyyy_yzzz = cbuffer.data(hg_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xyyyy_zzzz = cbuffer.data(hg_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxx = cbuffer.data(hg_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxy = cbuffer.data(hg_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxxz = cbuffer.data(hg_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxyy = cbuffer.data(hg_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxyz = cbuffer.data(hg_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xyyyz_xxzz = cbuffer.data(hg_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyyy = cbuffer.data(hg_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyyz = cbuffer.data(hg_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xyyyz_xyzz = cbuffer.data(hg_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xyyyz_xzzz = cbuffer.data(hg_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyyy = cbuffer.data(hg_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyyz = cbuffer.data(hg_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xyyyz_yyzz = cbuffer.data(hg_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xyyyz_yzzz = cbuffer.data(hg_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xyyyz_zzzz = cbuffer.data(hg_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxx = cbuffer.data(hg_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxy = cbuffer.data(hg_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxxz = cbuffer.data(hg_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxyy = cbuffer.data(hg_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxyz = cbuffer.data(hg_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xyyzz_xxzz = cbuffer.data(hg_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyyy = cbuffer.data(hg_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyyz = cbuffer.data(hg_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xyyzz_xyzz = cbuffer.data(hg_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xyyzz_xzzz = cbuffer.data(hg_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyyy = cbuffer.data(hg_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyyz = cbuffer.data(hg_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xyyzz_yyzz = cbuffer.data(hg_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xyyzz_yzzz = cbuffer.data(hg_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xyyzz_zzzz = cbuffer.data(hg_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxx = cbuffer.data(hg_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxy = cbuffer.data(hg_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxxz = cbuffer.data(hg_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxyy = cbuffer.data(hg_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxyz = cbuffer.data(hg_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xyzzz_xxzz = cbuffer.data(hg_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyyy = cbuffer.data(hg_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyyz = cbuffer.data(hg_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xyzzz_xyzz = cbuffer.data(hg_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xyzzz_xzzz = cbuffer.data(hg_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyyy = cbuffer.data(hg_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyyz = cbuffer.data(hg_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xyzzz_yyzz = cbuffer.data(hg_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xyzzz_yzzz = cbuffer.data(hg_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xyzzz_zzzz = cbuffer.data(hg_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxx = cbuffer.data(hg_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxy = cbuffer.data(hg_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxxz = cbuffer.data(hg_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxyy = cbuffer.data(hg_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxyz = cbuffer.data(hg_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xzzzz_xxzz = cbuffer.data(hg_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyyy = cbuffer.data(hg_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyyz = cbuffer.data(hg_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xzzzz_xyzz = cbuffer.data(hg_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xzzzz_xzzz = cbuffer.data(hg_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyyy = cbuffer.data(hg_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyyz = cbuffer.data(hg_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xzzzz_yyzz = cbuffer.data(hg_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xzzzz_yzzz = cbuffer.data(hg_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xzzzz_zzzz = cbuffer.data(hg_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxx = cbuffer.data(hg_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxy = cbuffer.data(hg_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxxz = cbuffer.data(hg_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxyy = cbuffer.data(hg_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxyz = cbuffer.data(hg_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_yyyyy_xxzz = cbuffer.data(hg_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyyy = cbuffer.data(hg_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyyz = cbuffer.data(hg_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_yyyyy_xyzz = cbuffer.data(hg_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_yyyyy_xzzz = cbuffer.data(hg_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyyy = cbuffer.data(hg_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyyz = cbuffer.data(hg_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_yyyyy_yyzz = cbuffer.data(hg_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_yyyyy_yzzz = cbuffer.data(hg_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_yyyyy_zzzz = cbuffer.data(hg_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxx = cbuffer.data(hg_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxy = cbuffer.data(hg_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxxz = cbuffer.data(hg_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxyy = cbuffer.data(hg_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxyz = cbuffer.data(hg_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_yyyyz_xxzz = cbuffer.data(hg_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyyy = cbuffer.data(hg_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyyz = cbuffer.data(hg_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_yyyyz_xyzz = cbuffer.data(hg_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_yyyyz_xzzz = cbuffer.data(hg_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyyy = cbuffer.data(hg_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyyz = cbuffer.data(hg_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_yyyyz_yyzz = cbuffer.data(hg_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_yyyyz_yzzz = cbuffer.data(hg_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_yyyyz_zzzz = cbuffer.data(hg_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxx = cbuffer.data(hg_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxy = cbuffer.data(hg_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxxz = cbuffer.data(hg_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxyy = cbuffer.data(hg_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxyz = cbuffer.data(hg_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_yyyzz_xxzz = cbuffer.data(hg_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyyy = cbuffer.data(hg_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyyz = cbuffer.data(hg_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_yyyzz_xyzz = cbuffer.data(hg_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_yyyzz_xzzz = cbuffer.data(hg_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyyy = cbuffer.data(hg_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyyz = cbuffer.data(hg_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_yyyzz_yyzz = cbuffer.data(hg_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_yyyzz_yzzz = cbuffer.data(hg_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_yyyzz_zzzz = cbuffer.data(hg_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxx = cbuffer.data(hg_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxy = cbuffer.data(hg_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxxz = cbuffer.data(hg_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxyy = cbuffer.data(hg_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxyz = cbuffer.data(hg_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_yyzzz_xxzz = cbuffer.data(hg_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyyy = cbuffer.data(hg_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyyz = cbuffer.data(hg_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_yyzzz_xyzz = cbuffer.data(hg_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_yyzzz_xzzz = cbuffer.data(hg_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyyy = cbuffer.data(hg_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyyz = cbuffer.data(hg_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_yyzzz_yyzz = cbuffer.data(hg_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_yyzzz_yzzz = cbuffer.data(hg_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_yyzzz_zzzz = cbuffer.data(hg_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxx = cbuffer.data(hg_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxy = cbuffer.data(hg_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxxz = cbuffer.data(hg_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxyy = cbuffer.data(hg_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxyz = cbuffer.data(hg_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_yzzzz_xxzz = cbuffer.data(hg_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyyy = cbuffer.data(hg_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyyz = cbuffer.data(hg_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_yzzzz_xyzz = cbuffer.data(hg_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_yzzzz_xzzz = cbuffer.data(hg_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyyy = cbuffer.data(hg_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyyz = cbuffer.data(hg_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_yzzzz_yyzz = cbuffer.data(hg_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_yzzzz_yzzz = cbuffer.data(hg_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_yzzzz_zzzz = cbuffer.data(hg_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxx = cbuffer.data(hg_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxy = cbuffer.data(hg_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxxz = cbuffer.data(hg_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxyy = cbuffer.data(hg_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxyz = cbuffer.data(hg_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_zzzzz_xxzz = cbuffer.data(hg_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyyy = cbuffer.data(hg_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyyz = cbuffer.data(hg_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_zzzzz_xyzz = cbuffer.data(hg_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_zzzzz_xzzz = cbuffer.data(hg_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyyy = cbuffer.data(hg_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyyz = cbuffer.data(hg_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_zzzzz_yyzz = cbuffer.data(hg_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_zzzzz_yzzz = cbuffer.data(hg_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_zzzzz_zzzz = cbuffer.data(hg_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxx = cbuffer.data(hg_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxy = cbuffer.data(hg_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxxz = cbuffer.data(hg_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxyy = cbuffer.data(hg_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxyz = cbuffer.data(hg_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_y_xxxxx_xxzz = cbuffer.data(hg_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyyy = cbuffer.data(hg_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyyz = cbuffer.data(hg_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_y_xxxxx_xyzz = cbuffer.data(hg_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_y_xxxxx_xzzz = cbuffer.data(hg_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyyy = cbuffer.data(hg_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyyz = cbuffer.data(hg_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_y_xxxxx_yyzz = cbuffer.data(hg_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_y_xxxxx_yzzz = cbuffer.data(hg_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_y_xxxxx_zzzz = cbuffer.data(hg_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxx = cbuffer.data(hg_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxy = cbuffer.data(hg_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxxz = cbuffer.data(hg_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxyy = cbuffer.data(hg_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxyz = cbuffer.data(hg_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_y_xxxxy_xxzz = cbuffer.data(hg_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyyy = cbuffer.data(hg_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyyz = cbuffer.data(hg_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_y_xxxxy_xyzz = cbuffer.data(hg_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_y_xxxxy_xzzz = cbuffer.data(hg_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyyy = cbuffer.data(hg_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyyz = cbuffer.data(hg_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_y_xxxxy_yyzz = cbuffer.data(hg_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_y_xxxxy_yzzz = cbuffer.data(hg_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_y_xxxxy_zzzz = cbuffer.data(hg_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxx = cbuffer.data(hg_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxy = cbuffer.data(hg_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxxz = cbuffer.data(hg_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxyy = cbuffer.data(hg_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxyz = cbuffer.data(hg_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_y_xxxxz_xxzz = cbuffer.data(hg_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyyy = cbuffer.data(hg_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyyz = cbuffer.data(hg_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_y_xxxxz_xyzz = cbuffer.data(hg_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_y_xxxxz_xzzz = cbuffer.data(hg_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyyy = cbuffer.data(hg_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyyz = cbuffer.data(hg_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_y_xxxxz_yyzz = cbuffer.data(hg_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_y_xxxxz_yzzz = cbuffer.data(hg_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_y_xxxxz_zzzz = cbuffer.data(hg_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxx = cbuffer.data(hg_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxy = cbuffer.data(hg_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxxz = cbuffer.data(hg_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxyy = cbuffer.data(hg_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxyz = cbuffer.data(hg_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_y_xxxyy_xxzz = cbuffer.data(hg_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyyy = cbuffer.data(hg_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyyz = cbuffer.data(hg_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_y_xxxyy_xyzz = cbuffer.data(hg_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_y_xxxyy_xzzz = cbuffer.data(hg_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyyy = cbuffer.data(hg_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyyz = cbuffer.data(hg_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_y_xxxyy_yyzz = cbuffer.data(hg_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_y_xxxyy_yzzz = cbuffer.data(hg_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_y_xxxyy_zzzz = cbuffer.data(hg_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxx = cbuffer.data(hg_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxy = cbuffer.data(hg_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxxz = cbuffer.data(hg_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxyy = cbuffer.data(hg_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxyz = cbuffer.data(hg_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_y_xxxyz_xxzz = cbuffer.data(hg_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyyy = cbuffer.data(hg_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyyz = cbuffer.data(hg_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_y_xxxyz_xyzz = cbuffer.data(hg_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_y_xxxyz_xzzz = cbuffer.data(hg_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyyy = cbuffer.data(hg_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyyz = cbuffer.data(hg_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_y_xxxyz_yyzz = cbuffer.data(hg_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_y_xxxyz_yzzz = cbuffer.data(hg_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_y_xxxyz_zzzz = cbuffer.data(hg_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxx = cbuffer.data(hg_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxy = cbuffer.data(hg_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxxz = cbuffer.data(hg_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxyy = cbuffer.data(hg_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxyz = cbuffer.data(hg_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_y_xxxzz_xxzz = cbuffer.data(hg_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyyy = cbuffer.data(hg_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyyz = cbuffer.data(hg_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_y_xxxzz_xyzz = cbuffer.data(hg_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_y_xxxzz_xzzz = cbuffer.data(hg_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyyy = cbuffer.data(hg_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyyz = cbuffer.data(hg_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_y_xxxzz_yyzz = cbuffer.data(hg_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_y_xxxzz_yzzz = cbuffer.data(hg_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_y_xxxzz_zzzz = cbuffer.data(hg_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxx = cbuffer.data(hg_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxy = cbuffer.data(hg_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxxz = cbuffer.data(hg_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxyy = cbuffer.data(hg_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxyz = cbuffer.data(hg_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_y_xxyyy_xxzz = cbuffer.data(hg_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyyy = cbuffer.data(hg_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyyz = cbuffer.data(hg_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_y_xxyyy_xyzz = cbuffer.data(hg_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_y_xxyyy_xzzz = cbuffer.data(hg_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyyy = cbuffer.data(hg_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyyz = cbuffer.data(hg_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_y_xxyyy_yyzz = cbuffer.data(hg_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_y_xxyyy_yzzz = cbuffer.data(hg_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_y_xxyyy_zzzz = cbuffer.data(hg_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxx = cbuffer.data(hg_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxy = cbuffer.data(hg_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxxz = cbuffer.data(hg_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxyy = cbuffer.data(hg_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxyz = cbuffer.data(hg_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_xxyyz_xxzz = cbuffer.data(hg_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyyy = cbuffer.data(hg_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyyz = cbuffer.data(hg_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_xxyyz_xyzz = cbuffer.data(hg_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_xxyyz_xzzz = cbuffer.data(hg_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyyy = cbuffer.data(hg_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyyz = cbuffer.data(hg_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_xxyyz_yyzz = cbuffer.data(hg_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_xxyyz_yzzz = cbuffer.data(hg_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_xxyyz_zzzz = cbuffer.data(hg_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxx = cbuffer.data(hg_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxy = cbuffer.data(hg_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxxz = cbuffer.data(hg_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxyy = cbuffer.data(hg_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxyz = cbuffer.data(hg_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_y_xxyzz_xxzz = cbuffer.data(hg_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyyy = cbuffer.data(hg_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyyz = cbuffer.data(hg_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_xxyzz_xyzz = cbuffer.data(hg_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_xxyzz_xzzz = cbuffer.data(hg_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyyy = cbuffer.data(hg_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyyz = cbuffer.data(hg_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_xxyzz_yyzz = cbuffer.data(hg_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_xxyzz_yzzz = cbuffer.data(hg_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_xxyzz_zzzz = cbuffer.data(hg_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxx = cbuffer.data(hg_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxy = cbuffer.data(hg_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxxz = cbuffer.data(hg_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxyy = cbuffer.data(hg_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxyz = cbuffer.data(hg_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_xxzzz_xxzz = cbuffer.data(hg_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyyy = cbuffer.data(hg_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyyz = cbuffer.data(hg_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_xxzzz_xyzz = cbuffer.data(hg_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_xxzzz_xzzz = cbuffer.data(hg_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyyy = cbuffer.data(hg_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyyz = cbuffer.data(hg_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_xxzzz_yyzz = cbuffer.data(hg_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_xxzzz_yzzz = cbuffer.data(hg_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_xxzzz_zzzz = cbuffer.data(hg_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxx = cbuffer.data(hg_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxy = cbuffer.data(hg_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxxz = cbuffer.data(hg_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxyy = cbuffer.data(hg_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxyz = cbuffer.data(hg_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_xyyyy_xxzz = cbuffer.data(hg_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyyy = cbuffer.data(hg_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyyz = cbuffer.data(hg_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_xyyyy_xyzz = cbuffer.data(hg_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_xyyyy_xzzz = cbuffer.data(hg_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyyy = cbuffer.data(hg_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyyz = cbuffer.data(hg_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_xyyyy_yyzz = cbuffer.data(hg_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_xyyyy_yzzz = cbuffer.data(hg_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_xyyyy_zzzz = cbuffer.data(hg_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxx = cbuffer.data(hg_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxy = cbuffer.data(hg_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxxz = cbuffer.data(hg_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxyy = cbuffer.data(hg_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxyz = cbuffer.data(hg_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_xyyyz_xxzz = cbuffer.data(hg_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyyy = cbuffer.data(hg_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyyz = cbuffer.data(hg_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_xyyyz_xyzz = cbuffer.data(hg_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_xyyyz_xzzz = cbuffer.data(hg_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyyy = cbuffer.data(hg_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyyz = cbuffer.data(hg_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_xyyyz_yyzz = cbuffer.data(hg_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_xyyyz_yzzz = cbuffer.data(hg_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_xyyyz_zzzz = cbuffer.data(hg_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxx = cbuffer.data(hg_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxy = cbuffer.data(hg_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxxz = cbuffer.data(hg_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxyy = cbuffer.data(hg_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxyz = cbuffer.data(hg_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_xyyzz_xxzz = cbuffer.data(hg_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyyy = cbuffer.data(hg_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyyz = cbuffer.data(hg_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_xyyzz_xyzz = cbuffer.data(hg_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_xyyzz_xzzz = cbuffer.data(hg_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyyy = cbuffer.data(hg_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyyz = cbuffer.data(hg_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_xyyzz_yyzz = cbuffer.data(hg_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_xyyzz_yzzz = cbuffer.data(hg_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_xyyzz_zzzz = cbuffer.data(hg_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxx = cbuffer.data(hg_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxy = cbuffer.data(hg_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxxz = cbuffer.data(hg_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxyy = cbuffer.data(hg_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxyz = cbuffer.data(hg_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_xyzzz_xxzz = cbuffer.data(hg_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyyy = cbuffer.data(hg_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyyz = cbuffer.data(hg_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_xyzzz_xyzz = cbuffer.data(hg_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_xyzzz_xzzz = cbuffer.data(hg_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyyy = cbuffer.data(hg_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyyz = cbuffer.data(hg_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_xyzzz_yyzz = cbuffer.data(hg_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_xyzzz_yzzz = cbuffer.data(hg_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_xyzzz_zzzz = cbuffer.data(hg_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxx = cbuffer.data(hg_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxy = cbuffer.data(hg_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxxz = cbuffer.data(hg_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxyy = cbuffer.data(hg_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxyz = cbuffer.data(hg_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_xzzzz_xxzz = cbuffer.data(hg_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyyy = cbuffer.data(hg_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyyz = cbuffer.data(hg_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_xzzzz_xyzz = cbuffer.data(hg_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_xzzzz_xzzz = cbuffer.data(hg_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyyy = cbuffer.data(hg_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyyz = cbuffer.data(hg_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_xzzzz_yyzz = cbuffer.data(hg_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_xzzzz_yzzz = cbuffer.data(hg_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_xzzzz_zzzz = cbuffer.data(hg_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxx = cbuffer.data(hg_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxy = cbuffer.data(hg_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxxz = cbuffer.data(hg_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxyy = cbuffer.data(hg_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxyz = cbuffer.data(hg_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_yyyyy_xxzz = cbuffer.data(hg_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyyy = cbuffer.data(hg_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyyz = cbuffer.data(hg_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_yyyyy_xyzz = cbuffer.data(hg_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_yyyyy_xzzz = cbuffer.data(hg_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyyy = cbuffer.data(hg_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyyz = cbuffer.data(hg_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_yyyyy_yyzz = cbuffer.data(hg_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_yyyyy_yzzz = cbuffer.data(hg_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_yyyyy_zzzz = cbuffer.data(hg_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxx = cbuffer.data(hg_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxy = cbuffer.data(hg_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxxz = cbuffer.data(hg_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxyy = cbuffer.data(hg_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxyz = cbuffer.data(hg_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_y_yyyyz_xxzz = cbuffer.data(hg_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyyy = cbuffer.data(hg_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyyz = cbuffer.data(hg_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_y_yyyyz_xyzz = cbuffer.data(hg_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_y_yyyyz_xzzz = cbuffer.data(hg_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyyy = cbuffer.data(hg_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyyz = cbuffer.data(hg_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_y_yyyyz_yyzz = cbuffer.data(hg_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_y_yyyyz_yzzz = cbuffer.data(hg_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_y_yyyyz_zzzz = cbuffer.data(hg_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxx = cbuffer.data(hg_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxy = cbuffer.data(hg_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxxz = cbuffer.data(hg_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxyy = cbuffer.data(hg_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxyz = cbuffer.data(hg_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_y_yyyzz_xxzz = cbuffer.data(hg_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyyy = cbuffer.data(hg_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyyz = cbuffer.data(hg_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_y_yyyzz_xyzz = cbuffer.data(hg_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_y_yyyzz_xzzz = cbuffer.data(hg_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyyy = cbuffer.data(hg_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyyz = cbuffer.data(hg_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_y_yyyzz_yyzz = cbuffer.data(hg_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_y_yyyzz_yzzz = cbuffer.data(hg_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_y_yyyzz_zzzz = cbuffer.data(hg_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxx = cbuffer.data(hg_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxy = cbuffer.data(hg_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxxz = cbuffer.data(hg_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxyy = cbuffer.data(hg_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxyz = cbuffer.data(hg_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_y_yyzzz_xxzz = cbuffer.data(hg_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyyy = cbuffer.data(hg_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyyz = cbuffer.data(hg_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_yyzzz_xyzz = cbuffer.data(hg_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_yyzzz_xzzz = cbuffer.data(hg_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyyy = cbuffer.data(hg_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyyz = cbuffer.data(hg_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_yyzzz_yyzz = cbuffer.data(hg_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_yyzzz_yzzz = cbuffer.data(hg_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_yyzzz_zzzz = cbuffer.data(hg_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxx = cbuffer.data(hg_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxy = cbuffer.data(hg_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxxz = cbuffer.data(hg_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxyy = cbuffer.data(hg_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxyz = cbuffer.data(hg_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_yzzzz_xxzz = cbuffer.data(hg_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyyy = cbuffer.data(hg_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyyz = cbuffer.data(hg_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_yzzzz_xyzz = cbuffer.data(hg_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_y_yzzzz_xzzz = cbuffer.data(hg_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyyy = cbuffer.data(hg_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyyz = cbuffer.data(hg_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_yzzzz_yyzz = cbuffer.data(hg_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_yzzzz_yzzz = cbuffer.data(hg_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_yzzzz_zzzz = cbuffer.data(hg_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxx = cbuffer.data(hg_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxy = cbuffer.data(hg_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxxz = cbuffer.data(hg_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxyy = cbuffer.data(hg_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxyz = cbuffer.data(hg_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_y_zzzzz_xxzz = cbuffer.data(hg_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyyy = cbuffer.data(hg_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyyz = cbuffer.data(hg_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_zzzzz_xyzz = cbuffer.data(hg_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_zzzzz_xzzz = cbuffer.data(hg_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyyy = cbuffer.data(hg_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyyz = cbuffer.data(hg_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_zzzzz_yyzz = cbuffer.data(hg_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_zzzzz_yzzz = cbuffer.data(hg_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_zzzzz_zzzz = cbuffer.data(hg_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxx = cbuffer.data(hg_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxy = cbuffer.data(hg_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxxz = cbuffer.data(hg_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxyy = cbuffer.data(hg_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxyz = cbuffer.data(hg_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_z_xxxxx_xxzz = cbuffer.data(hg_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyyy = cbuffer.data(hg_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyyz = cbuffer.data(hg_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_z_xxxxx_xyzz = cbuffer.data(hg_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_z_xxxxx_xzzz = cbuffer.data(hg_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyyy = cbuffer.data(hg_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyyz = cbuffer.data(hg_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_z_xxxxx_yyzz = cbuffer.data(hg_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_z_xxxxx_yzzz = cbuffer.data(hg_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_z_xxxxx_zzzz = cbuffer.data(hg_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxx = cbuffer.data(hg_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxy = cbuffer.data(hg_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxxz = cbuffer.data(hg_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxyy = cbuffer.data(hg_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxyz = cbuffer.data(hg_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_z_xxxxy_xxzz = cbuffer.data(hg_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyyy = cbuffer.data(hg_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyyz = cbuffer.data(hg_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_z_xxxxy_xyzz = cbuffer.data(hg_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_z_xxxxy_xzzz = cbuffer.data(hg_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyyy = cbuffer.data(hg_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyyz = cbuffer.data(hg_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_z_xxxxy_yyzz = cbuffer.data(hg_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_z_xxxxy_yzzz = cbuffer.data(hg_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_z_xxxxy_zzzz = cbuffer.data(hg_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxx = cbuffer.data(hg_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxy = cbuffer.data(hg_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxxz = cbuffer.data(hg_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxyy = cbuffer.data(hg_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxyz = cbuffer.data(hg_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_z_xxxxz_xxzz = cbuffer.data(hg_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyyy = cbuffer.data(hg_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyyz = cbuffer.data(hg_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_z_xxxxz_xyzz = cbuffer.data(hg_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_z_xxxxz_xzzz = cbuffer.data(hg_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyyy = cbuffer.data(hg_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyyz = cbuffer.data(hg_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_z_xxxxz_yyzz = cbuffer.data(hg_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_z_xxxxz_yzzz = cbuffer.data(hg_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_z_xxxxz_zzzz = cbuffer.data(hg_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxx = cbuffer.data(hg_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxy = cbuffer.data(hg_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxxz = cbuffer.data(hg_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxyy = cbuffer.data(hg_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxyz = cbuffer.data(hg_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_z_xxxyy_xxzz = cbuffer.data(hg_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyyy = cbuffer.data(hg_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyyz = cbuffer.data(hg_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_z_xxxyy_xyzz = cbuffer.data(hg_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_z_xxxyy_xzzz = cbuffer.data(hg_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyyy = cbuffer.data(hg_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyyz = cbuffer.data(hg_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_z_xxxyy_yyzz = cbuffer.data(hg_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_z_xxxyy_yzzz = cbuffer.data(hg_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_z_xxxyy_zzzz = cbuffer.data(hg_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxx = cbuffer.data(hg_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxy = cbuffer.data(hg_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxxz = cbuffer.data(hg_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxyy = cbuffer.data(hg_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxyz = cbuffer.data(hg_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_z_xxxyz_xxzz = cbuffer.data(hg_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyyy = cbuffer.data(hg_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyyz = cbuffer.data(hg_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_z_xxxyz_xyzz = cbuffer.data(hg_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_z_xxxyz_xzzz = cbuffer.data(hg_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyyy = cbuffer.data(hg_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyyz = cbuffer.data(hg_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_z_xxxyz_yyzz = cbuffer.data(hg_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_z_xxxyz_yzzz = cbuffer.data(hg_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_z_xxxyz_zzzz = cbuffer.data(hg_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxx = cbuffer.data(hg_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxy = cbuffer.data(hg_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxxz = cbuffer.data(hg_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxyy = cbuffer.data(hg_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxyz = cbuffer.data(hg_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_z_xxxzz_xxzz = cbuffer.data(hg_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyyy = cbuffer.data(hg_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyyz = cbuffer.data(hg_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_z_xxxzz_xyzz = cbuffer.data(hg_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_z_xxxzz_xzzz = cbuffer.data(hg_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyyy = cbuffer.data(hg_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyyz = cbuffer.data(hg_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_z_xxxzz_yyzz = cbuffer.data(hg_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_z_xxxzz_yzzz = cbuffer.data(hg_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_z_xxxzz_zzzz = cbuffer.data(hg_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxx = cbuffer.data(hg_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxy = cbuffer.data(hg_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxxz = cbuffer.data(hg_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxyy = cbuffer.data(hg_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxyz = cbuffer.data(hg_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_z_xxyyy_xxzz = cbuffer.data(hg_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyyy = cbuffer.data(hg_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyyz = cbuffer.data(hg_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_z_xxyyy_xyzz = cbuffer.data(hg_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_z_xxyyy_xzzz = cbuffer.data(hg_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyyy = cbuffer.data(hg_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyyz = cbuffer.data(hg_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_z_xxyyy_yyzz = cbuffer.data(hg_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_z_xxyyy_yzzz = cbuffer.data(hg_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_z_xxyyy_zzzz = cbuffer.data(hg_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxx = cbuffer.data(hg_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxy = cbuffer.data(hg_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxxz = cbuffer.data(hg_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxyy = cbuffer.data(hg_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxyz = cbuffer.data(hg_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_z_xxyyz_xxzz = cbuffer.data(hg_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyyy = cbuffer.data(hg_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyyz = cbuffer.data(hg_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_z_xxyyz_xyzz = cbuffer.data(hg_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_z_xxyyz_xzzz = cbuffer.data(hg_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyyy = cbuffer.data(hg_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyyz = cbuffer.data(hg_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_z_xxyyz_yyzz = cbuffer.data(hg_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_z_xxyyz_yzzz = cbuffer.data(hg_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_z_xxyyz_zzzz = cbuffer.data(hg_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxx = cbuffer.data(hg_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxy = cbuffer.data(hg_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxxz = cbuffer.data(hg_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxyy = cbuffer.data(hg_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxyz = cbuffer.data(hg_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_z_xxyzz_xxzz = cbuffer.data(hg_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyyy = cbuffer.data(hg_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyyz = cbuffer.data(hg_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_z_xxyzz_xyzz = cbuffer.data(hg_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_z_xxyzz_xzzz = cbuffer.data(hg_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyyy = cbuffer.data(hg_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyyz = cbuffer.data(hg_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_z_xxyzz_yyzz = cbuffer.data(hg_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_z_xxyzz_yzzz = cbuffer.data(hg_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_z_xxyzz_zzzz = cbuffer.data(hg_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxx = cbuffer.data(hg_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxy = cbuffer.data(hg_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxxz = cbuffer.data(hg_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxyy = cbuffer.data(hg_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxyz = cbuffer.data(hg_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_z_xxzzz_xxzz = cbuffer.data(hg_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyyy = cbuffer.data(hg_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyyz = cbuffer.data(hg_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_z_xxzzz_xyzz = cbuffer.data(hg_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_z_xxzzz_xzzz = cbuffer.data(hg_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyyy = cbuffer.data(hg_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyyz = cbuffer.data(hg_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_z_xxzzz_yyzz = cbuffer.data(hg_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_z_xxzzz_yzzz = cbuffer.data(hg_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_z_xxzzz_zzzz = cbuffer.data(hg_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxx = cbuffer.data(hg_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxy = cbuffer.data(hg_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxxz = cbuffer.data(hg_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxyy = cbuffer.data(hg_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxyz = cbuffer.data(hg_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_z_xyyyy_xxzz = cbuffer.data(hg_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyyy = cbuffer.data(hg_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyyz = cbuffer.data(hg_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_z_xyyyy_xyzz = cbuffer.data(hg_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_z_xyyyy_xzzz = cbuffer.data(hg_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyyy = cbuffer.data(hg_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyyz = cbuffer.data(hg_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_z_xyyyy_yyzz = cbuffer.data(hg_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_z_xyyyy_yzzz = cbuffer.data(hg_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_z_xyyyy_zzzz = cbuffer.data(hg_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxx = cbuffer.data(hg_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxy = cbuffer.data(hg_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxxz = cbuffer.data(hg_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxyy = cbuffer.data(hg_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxyz = cbuffer.data(hg_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_z_xyyyz_xxzz = cbuffer.data(hg_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyyy = cbuffer.data(hg_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyyz = cbuffer.data(hg_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_z_xyyyz_xyzz = cbuffer.data(hg_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_z_xyyyz_xzzz = cbuffer.data(hg_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyyy = cbuffer.data(hg_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyyz = cbuffer.data(hg_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_z_xyyyz_yyzz = cbuffer.data(hg_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_z_xyyyz_yzzz = cbuffer.data(hg_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_z_xyyyz_zzzz = cbuffer.data(hg_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxx = cbuffer.data(hg_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxy = cbuffer.data(hg_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxxz = cbuffer.data(hg_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxyy = cbuffer.data(hg_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxyz = cbuffer.data(hg_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_z_xyyzz_xxzz = cbuffer.data(hg_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyyy = cbuffer.data(hg_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyyz = cbuffer.data(hg_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_z_xyyzz_xyzz = cbuffer.data(hg_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_z_xyyzz_xzzz = cbuffer.data(hg_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyyy = cbuffer.data(hg_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyyz = cbuffer.data(hg_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_z_xyyzz_yyzz = cbuffer.data(hg_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_z_xyyzz_yzzz = cbuffer.data(hg_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_z_xyyzz_zzzz = cbuffer.data(hg_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxx = cbuffer.data(hg_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxy = cbuffer.data(hg_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxxz = cbuffer.data(hg_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxyy = cbuffer.data(hg_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxyz = cbuffer.data(hg_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_z_xyzzz_xxzz = cbuffer.data(hg_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyyy = cbuffer.data(hg_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyyz = cbuffer.data(hg_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_z_xyzzz_xyzz = cbuffer.data(hg_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_z_xyzzz_xzzz = cbuffer.data(hg_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyyy = cbuffer.data(hg_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyyz = cbuffer.data(hg_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_z_xyzzz_yyzz = cbuffer.data(hg_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_z_xyzzz_yzzz = cbuffer.data(hg_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_z_xyzzz_zzzz = cbuffer.data(hg_geom_01_off + 839 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxx = cbuffer.data(hg_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxy = cbuffer.data(hg_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxxz = cbuffer.data(hg_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxyy = cbuffer.data(hg_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxyz = cbuffer.data(hg_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_z_xzzzz_xxzz = cbuffer.data(hg_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyyy = cbuffer.data(hg_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyyz = cbuffer.data(hg_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_z_xzzzz_xyzz = cbuffer.data(hg_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_z_xzzzz_xzzz = cbuffer.data(hg_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyyy = cbuffer.data(hg_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyyz = cbuffer.data(hg_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_z_xzzzz_yyzz = cbuffer.data(hg_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_z_xzzzz_yzzz = cbuffer.data(hg_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_z_xzzzz_zzzz = cbuffer.data(hg_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxx = cbuffer.data(hg_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxy = cbuffer.data(hg_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxxz = cbuffer.data(hg_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxyy = cbuffer.data(hg_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxyz = cbuffer.data(hg_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_z_yyyyy_xxzz = cbuffer.data(hg_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyyy = cbuffer.data(hg_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyyz = cbuffer.data(hg_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_z_yyyyy_xyzz = cbuffer.data(hg_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_z_yyyyy_xzzz = cbuffer.data(hg_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyyy = cbuffer.data(hg_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyyz = cbuffer.data(hg_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_z_yyyyy_yyzz = cbuffer.data(hg_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_z_yyyyy_yzzz = cbuffer.data(hg_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_z_yyyyy_zzzz = cbuffer.data(hg_geom_01_off + 869 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxx = cbuffer.data(hg_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxy = cbuffer.data(hg_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxxz = cbuffer.data(hg_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxyy = cbuffer.data(hg_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxyz = cbuffer.data(hg_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_z_yyyyz_xxzz = cbuffer.data(hg_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyyy = cbuffer.data(hg_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyyz = cbuffer.data(hg_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_z_yyyyz_xyzz = cbuffer.data(hg_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_z_yyyyz_xzzz = cbuffer.data(hg_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyyy = cbuffer.data(hg_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyyz = cbuffer.data(hg_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_z_yyyyz_yyzz = cbuffer.data(hg_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_z_yyyyz_yzzz = cbuffer.data(hg_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_z_yyyyz_zzzz = cbuffer.data(hg_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxx = cbuffer.data(hg_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxy = cbuffer.data(hg_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxxz = cbuffer.data(hg_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxyy = cbuffer.data(hg_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxyz = cbuffer.data(hg_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_z_yyyzz_xxzz = cbuffer.data(hg_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyyy = cbuffer.data(hg_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyyz = cbuffer.data(hg_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_z_yyyzz_xyzz = cbuffer.data(hg_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_z_yyyzz_xzzz = cbuffer.data(hg_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyyy = cbuffer.data(hg_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyyz = cbuffer.data(hg_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_z_yyyzz_yyzz = cbuffer.data(hg_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_z_yyyzz_yzzz = cbuffer.data(hg_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_z_yyyzz_zzzz = cbuffer.data(hg_geom_01_off + 899 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxx = cbuffer.data(hg_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxy = cbuffer.data(hg_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxxz = cbuffer.data(hg_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxyy = cbuffer.data(hg_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxyz = cbuffer.data(hg_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_z_yyzzz_xxzz = cbuffer.data(hg_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyyy = cbuffer.data(hg_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyyz = cbuffer.data(hg_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_z_yyzzz_xyzz = cbuffer.data(hg_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_z_yyzzz_xzzz = cbuffer.data(hg_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyyy = cbuffer.data(hg_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyyz = cbuffer.data(hg_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_z_yyzzz_yyzz = cbuffer.data(hg_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_z_yyzzz_yzzz = cbuffer.data(hg_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_z_yyzzz_zzzz = cbuffer.data(hg_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxx = cbuffer.data(hg_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxy = cbuffer.data(hg_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxxz = cbuffer.data(hg_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxyy = cbuffer.data(hg_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxyz = cbuffer.data(hg_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_z_yzzzz_xxzz = cbuffer.data(hg_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyyy = cbuffer.data(hg_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyyz = cbuffer.data(hg_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_z_yzzzz_xyzz = cbuffer.data(hg_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_z_yzzzz_xzzz = cbuffer.data(hg_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyyy = cbuffer.data(hg_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyyz = cbuffer.data(hg_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_z_yzzzz_yyzz = cbuffer.data(hg_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_z_yzzzz_yzzz = cbuffer.data(hg_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_z_yzzzz_zzzz = cbuffer.data(hg_geom_01_off + 929 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxx = cbuffer.data(hg_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxy = cbuffer.data(hg_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxxz = cbuffer.data(hg_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxyy = cbuffer.data(hg_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxyz = cbuffer.data(hg_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_z_zzzzz_xxzz = cbuffer.data(hg_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyyy = cbuffer.data(hg_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyyz = cbuffer.data(hg_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_z_zzzzz_xyzz = cbuffer.data(hg_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_z_zzzzz_xzzz = cbuffer.data(hg_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyyy = cbuffer.data(hg_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyyz = cbuffer.data(hg_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_z_zzzzz_yyzz = cbuffer.data(hg_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_z_zzzzz_yzzz = cbuffer.data(hg_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_z_zzzzz_zzzz = cbuffer.data(hg_geom_01_off + 944 * ccomps * dcomps);

            /// Set up components of auxilary buffer : HHSS

            const auto hh_geom_01_off = idx_geom_01_hhxx + i * dcomps + j;

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

            /// set up bra offset for contr_buffer_igxx

            const auto ig_geom_01_off = idx_geom_01_igxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxx_xxxx = cbuffer.data(ig_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxxy = cbuffer.data(ig_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxxz = cbuffer.data(ig_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxyy = cbuffer.data(ig_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxyz = cbuffer.data(ig_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxzz = cbuffer.data(ig_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyyy = cbuffer.data(ig_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyyz = cbuffer.data(ig_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyzz = cbuffer.data(ig_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xzzz = cbuffer.data(ig_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyyy = cbuffer.data(ig_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyyz = cbuffer.data(ig_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyzz = cbuffer.data(ig_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yzzz = cbuffer.data(ig_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxx_zzzz = cbuffer.data(ig_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxx_xxxx, g_0_x_xxxxx_xxxxx, g_0_x_xxxxx_xxxxy, g_0_x_xxxxx_xxxxz, g_0_x_xxxxx_xxxy, g_0_x_xxxxx_xxxyy, g_0_x_xxxxx_xxxyz, g_0_x_xxxxx_xxxz, g_0_x_xxxxx_xxxzz, g_0_x_xxxxx_xxyy, g_0_x_xxxxx_xxyyy, g_0_x_xxxxx_xxyyz, g_0_x_xxxxx_xxyz, g_0_x_xxxxx_xxyzz, g_0_x_xxxxx_xxzz, g_0_x_xxxxx_xxzzz, g_0_x_xxxxx_xyyy, g_0_x_xxxxx_xyyyy, g_0_x_xxxxx_xyyyz, g_0_x_xxxxx_xyyz, g_0_x_xxxxx_xyyzz, g_0_x_xxxxx_xyzz, g_0_x_xxxxx_xyzzz, g_0_x_xxxxx_xzzz, g_0_x_xxxxx_xzzzz, g_0_x_xxxxx_yyyy, g_0_x_xxxxx_yyyz, g_0_x_xxxxx_yyzz, g_0_x_xxxxx_yzzz, g_0_x_xxxxx_zzzz, g_0_x_xxxxxx_xxxx, g_0_x_xxxxxx_xxxy, g_0_x_xxxxxx_xxxz, g_0_x_xxxxxx_xxyy, g_0_x_xxxxxx_xxyz, g_0_x_xxxxxx_xxzz, g_0_x_xxxxxx_xyyy, g_0_x_xxxxxx_xyyz, g_0_x_xxxxxx_xyzz, g_0_x_xxxxxx_xzzz, g_0_x_xxxxxx_yyyy, g_0_x_xxxxxx_yyyz, g_0_x_xxxxxx_yyzz, g_0_x_xxxxxx_yzzz, g_0_x_xxxxxx_zzzz, g_xxxxx_xxxx, g_xxxxx_xxxy, g_xxxxx_xxxz, g_xxxxx_xxyy, g_xxxxx_xxyz, g_xxxxx_xxzz, g_xxxxx_xyyy, g_xxxxx_xyyz, g_xxxxx_xyzz, g_xxxxx_xzzz, g_xxxxx_yyyy, g_xxxxx_yyyz, g_xxxxx_yyzz, g_xxxxx_yzzz, g_xxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxx_xxxx[k] = g_xxxxx_xxxx[k] - g_0_x_xxxxx_xxxx[k] * ab_x + g_0_x_xxxxx_xxxxx[k];

                g_0_x_xxxxxx_xxxy[k] = g_xxxxx_xxxy[k] - g_0_x_xxxxx_xxxy[k] * ab_x + g_0_x_xxxxx_xxxxy[k];

                g_0_x_xxxxxx_xxxz[k] = g_xxxxx_xxxz[k] - g_0_x_xxxxx_xxxz[k] * ab_x + g_0_x_xxxxx_xxxxz[k];

                g_0_x_xxxxxx_xxyy[k] = g_xxxxx_xxyy[k] - g_0_x_xxxxx_xxyy[k] * ab_x + g_0_x_xxxxx_xxxyy[k];

                g_0_x_xxxxxx_xxyz[k] = g_xxxxx_xxyz[k] - g_0_x_xxxxx_xxyz[k] * ab_x + g_0_x_xxxxx_xxxyz[k];

                g_0_x_xxxxxx_xxzz[k] = g_xxxxx_xxzz[k] - g_0_x_xxxxx_xxzz[k] * ab_x + g_0_x_xxxxx_xxxzz[k];

                g_0_x_xxxxxx_xyyy[k] = g_xxxxx_xyyy[k] - g_0_x_xxxxx_xyyy[k] * ab_x + g_0_x_xxxxx_xxyyy[k];

                g_0_x_xxxxxx_xyyz[k] = g_xxxxx_xyyz[k] - g_0_x_xxxxx_xyyz[k] * ab_x + g_0_x_xxxxx_xxyyz[k];

                g_0_x_xxxxxx_xyzz[k] = g_xxxxx_xyzz[k] - g_0_x_xxxxx_xyzz[k] * ab_x + g_0_x_xxxxx_xxyzz[k];

                g_0_x_xxxxxx_xzzz[k] = g_xxxxx_xzzz[k] - g_0_x_xxxxx_xzzz[k] * ab_x + g_0_x_xxxxx_xxzzz[k];

                g_0_x_xxxxxx_yyyy[k] = g_xxxxx_yyyy[k] - g_0_x_xxxxx_yyyy[k] * ab_x + g_0_x_xxxxx_xyyyy[k];

                g_0_x_xxxxxx_yyyz[k] = g_xxxxx_yyyz[k] - g_0_x_xxxxx_yyyz[k] * ab_x + g_0_x_xxxxx_xyyyz[k];

                g_0_x_xxxxxx_yyzz[k] = g_xxxxx_yyzz[k] - g_0_x_xxxxx_yyzz[k] * ab_x + g_0_x_xxxxx_xyyzz[k];

                g_0_x_xxxxxx_yzzz[k] = g_xxxxx_yzzz[k] - g_0_x_xxxxx_yzzz[k] * ab_x + g_0_x_xxxxx_xyzzz[k];

                g_0_x_xxxxxx_zzzz[k] = g_xxxxx_zzzz[k] - g_0_x_xxxxx_zzzz[k] * ab_x + g_0_x_xxxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxy_xxxx = cbuffer.data(ig_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxxy = cbuffer.data(ig_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxxz = cbuffer.data(ig_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxyy = cbuffer.data(ig_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxyz = cbuffer.data(ig_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxzz = cbuffer.data(ig_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyyy = cbuffer.data(ig_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyyz = cbuffer.data(ig_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyzz = cbuffer.data(ig_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xzzz = cbuffer.data(ig_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyyy = cbuffer.data(ig_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyyz = cbuffer.data(ig_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyzz = cbuffer.data(ig_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yzzz = cbuffer.data(ig_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxxy_zzzz = cbuffer.data(ig_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxx_xxxx, g_0_x_xxxxx_xxxxy, g_0_x_xxxxx_xxxy, g_0_x_xxxxx_xxxyy, g_0_x_xxxxx_xxxyz, g_0_x_xxxxx_xxxz, g_0_x_xxxxx_xxyy, g_0_x_xxxxx_xxyyy, g_0_x_xxxxx_xxyyz, g_0_x_xxxxx_xxyz, g_0_x_xxxxx_xxyzz, g_0_x_xxxxx_xxzz, g_0_x_xxxxx_xyyy, g_0_x_xxxxx_xyyyy, g_0_x_xxxxx_xyyyz, g_0_x_xxxxx_xyyz, g_0_x_xxxxx_xyyzz, g_0_x_xxxxx_xyzz, g_0_x_xxxxx_xyzzz, g_0_x_xxxxx_xzzz, g_0_x_xxxxx_yyyy, g_0_x_xxxxx_yyyyy, g_0_x_xxxxx_yyyyz, g_0_x_xxxxx_yyyz, g_0_x_xxxxx_yyyzz, g_0_x_xxxxx_yyzz, g_0_x_xxxxx_yyzzz, g_0_x_xxxxx_yzzz, g_0_x_xxxxx_yzzzz, g_0_x_xxxxx_zzzz, g_0_x_xxxxxy_xxxx, g_0_x_xxxxxy_xxxy, g_0_x_xxxxxy_xxxz, g_0_x_xxxxxy_xxyy, g_0_x_xxxxxy_xxyz, g_0_x_xxxxxy_xxzz, g_0_x_xxxxxy_xyyy, g_0_x_xxxxxy_xyyz, g_0_x_xxxxxy_xyzz, g_0_x_xxxxxy_xzzz, g_0_x_xxxxxy_yyyy, g_0_x_xxxxxy_yyyz, g_0_x_xxxxxy_yyzz, g_0_x_xxxxxy_yzzz, g_0_x_xxxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxy_xxxx[k] = -g_0_x_xxxxx_xxxx[k] * ab_y + g_0_x_xxxxx_xxxxy[k];

                g_0_x_xxxxxy_xxxy[k] = -g_0_x_xxxxx_xxxy[k] * ab_y + g_0_x_xxxxx_xxxyy[k];

                g_0_x_xxxxxy_xxxz[k] = -g_0_x_xxxxx_xxxz[k] * ab_y + g_0_x_xxxxx_xxxyz[k];

                g_0_x_xxxxxy_xxyy[k] = -g_0_x_xxxxx_xxyy[k] * ab_y + g_0_x_xxxxx_xxyyy[k];

                g_0_x_xxxxxy_xxyz[k] = -g_0_x_xxxxx_xxyz[k] * ab_y + g_0_x_xxxxx_xxyyz[k];

                g_0_x_xxxxxy_xxzz[k] = -g_0_x_xxxxx_xxzz[k] * ab_y + g_0_x_xxxxx_xxyzz[k];

                g_0_x_xxxxxy_xyyy[k] = -g_0_x_xxxxx_xyyy[k] * ab_y + g_0_x_xxxxx_xyyyy[k];

                g_0_x_xxxxxy_xyyz[k] = -g_0_x_xxxxx_xyyz[k] * ab_y + g_0_x_xxxxx_xyyyz[k];

                g_0_x_xxxxxy_xyzz[k] = -g_0_x_xxxxx_xyzz[k] * ab_y + g_0_x_xxxxx_xyyzz[k];

                g_0_x_xxxxxy_xzzz[k] = -g_0_x_xxxxx_xzzz[k] * ab_y + g_0_x_xxxxx_xyzzz[k];

                g_0_x_xxxxxy_yyyy[k] = -g_0_x_xxxxx_yyyy[k] * ab_y + g_0_x_xxxxx_yyyyy[k];

                g_0_x_xxxxxy_yyyz[k] = -g_0_x_xxxxx_yyyz[k] * ab_y + g_0_x_xxxxx_yyyyz[k];

                g_0_x_xxxxxy_yyzz[k] = -g_0_x_xxxxx_yyzz[k] * ab_y + g_0_x_xxxxx_yyyzz[k];

                g_0_x_xxxxxy_yzzz[k] = -g_0_x_xxxxx_yzzz[k] * ab_y + g_0_x_xxxxx_yyzzz[k];

                g_0_x_xxxxxy_zzzz[k] = -g_0_x_xxxxx_zzzz[k] * ab_y + g_0_x_xxxxx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxz_xxxx = cbuffer.data(ig_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxxy = cbuffer.data(ig_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxxz = cbuffer.data(ig_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxyy = cbuffer.data(ig_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxyz = cbuffer.data(ig_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxzz = cbuffer.data(ig_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyyy = cbuffer.data(ig_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyyz = cbuffer.data(ig_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyzz = cbuffer.data(ig_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xzzz = cbuffer.data(ig_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyyy = cbuffer.data(ig_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyyz = cbuffer.data(ig_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyzz = cbuffer.data(ig_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yzzz = cbuffer.data(ig_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxxxz_zzzz = cbuffer.data(ig_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxx_xxxx, g_0_x_xxxxx_xxxxz, g_0_x_xxxxx_xxxy, g_0_x_xxxxx_xxxyz, g_0_x_xxxxx_xxxz, g_0_x_xxxxx_xxxzz, g_0_x_xxxxx_xxyy, g_0_x_xxxxx_xxyyz, g_0_x_xxxxx_xxyz, g_0_x_xxxxx_xxyzz, g_0_x_xxxxx_xxzz, g_0_x_xxxxx_xxzzz, g_0_x_xxxxx_xyyy, g_0_x_xxxxx_xyyyz, g_0_x_xxxxx_xyyz, g_0_x_xxxxx_xyyzz, g_0_x_xxxxx_xyzz, g_0_x_xxxxx_xyzzz, g_0_x_xxxxx_xzzz, g_0_x_xxxxx_xzzzz, g_0_x_xxxxx_yyyy, g_0_x_xxxxx_yyyyz, g_0_x_xxxxx_yyyz, g_0_x_xxxxx_yyyzz, g_0_x_xxxxx_yyzz, g_0_x_xxxxx_yyzzz, g_0_x_xxxxx_yzzz, g_0_x_xxxxx_yzzzz, g_0_x_xxxxx_zzzz, g_0_x_xxxxx_zzzzz, g_0_x_xxxxxz_xxxx, g_0_x_xxxxxz_xxxy, g_0_x_xxxxxz_xxxz, g_0_x_xxxxxz_xxyy, g_0_x_xxxxxz_xxyz, g_0_x_xxxxxz_xxzz, g_0_x_xxxxxz_xyyy, g_0_x_xxxxxz_xyyz, g_0_x_xxxxxz_xyzz, g_0_x_xxxxxz_xzzz, g_0_x_xxxxxz_yyyy, g_0_x_xxxxxz_yyyz, g_0_x_xxxxxz_yyzz, g_0_x_xxxxxz_yzzz, g_0_x_xxxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxz_xxxx[k] = -g_0_x_xxxxx_xxxx[k] * ab_z + g_0_x_xxxxx_xxxxz[k];

                g_0_x_xxxxxz_xxxy[k] = -g_0_x_xxxxx_xxxy[k] * ab_z + g_0_x_xxxxx_xxxyz[k];

                g_0_x_xxxxxz_xxxz[k] = -g_0_x_xxxxx_xxxz[k] * ab_z + g_0_x_xxxxx_xxxzz[k];

                g_0_x_xxxxxz_xxyy[k] = -g_0_x_xxxxx_xxyy[k] * ab_z + g_0_x_xxxxx_xxyyz[k];

                g_0_x_xxxxxz_xxyz[k] = -g_0_x_xxxxx_xxyz[k] * ab_z + g_0_x_xxxxx_xxyzz[k];

                g_0_x_xxxxxz_xxzz[k] = -g_0_x_xxxxx_xxzz[k] * ab_z + g_0_x_xxxxx_xxzzz[k];

                g_0_x_xxxxxz_xyyy[k] = -g_0_x_xxxxx_xyyy[k] * ab_z + g_0_x_xxxxx_xyyyz[k];

                g_0_x_xxxxxz_xyyz[k] = -g_0_x_xxxxx_xyyz[k] * ab_z + g_0_x_xxxxx_xyyzz[k];

                g_0_x_xxxxxz_xyzz[k] = -g_0_x_xxxxx_xyzz[k] * ab_z + g_0_x_xxxxx_xyzzz[k];

                g_0_x_xxxxxz_xzzz[k] = -g_0_x_xxxxx_xzzz[k] * ab_z + g_0_x_xxxxx_xzzzz[k];

                g_0_x_xxxxxz_yyyy[k] = -g_0_x_xxxxx_yyyy[k] * ab_z + g_0_x_xxxxx_yyyyz[k];

                g_0_x_xxxxxz_yyyz[k] = -g_0_x_xxxxx_yyyz[k] * ab_z + g_0_x_xxxxx_yyyzz[k];

                g_0_x_xxxxxz_yyzz[k] = -g_0_x_xxxxx_yyzz[k] * ab_z + g_0_x_xxxxx_yyzzz[k];

                g_0_x_xxxxxz_yzzz[k] = -g_0_x_xxxxx_yzzz[k] * ab_z + g_0_x_xxxxx_yzzzz[k];

                g_0_x_xxxxxz_zzzz[k] = -g_0_x_xxxxx_zzzz[k] * ab_z + g_0_x_xxxxx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyy_xxxx = cbuffer.data(ig_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxxy = cbuffer.data(ig_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxxz = cbuffer.data(ig_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxyy = cbuffer.data(ig_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxyz = cbuffer.data(ig_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxzz = cbuffer.data(ig_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyyy = cbuffer.data(ig_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyyz = cbuffer.data(ig_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyzz = cbuffer.data(ig_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xzzz = cbuffer.data(ig_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyyy = cbuffer.data(ig_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyyz = cbuffer.data(ig_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyzz = cbuffer.data(ig_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yzzz = cbuffer.data(ig_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxxyy_zzzz = cbuffer.data(ig_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxy_xxxx, g_0_x_xxxxy_xxxxy, g_0_x_xxxxy_xxxy, g_0_x_xxxxy_xxxyy, g_0_x_xxxxy_xxxyz, g_0_x_xxxxy_xxxz, g_0_x_xxxxy_xxyy, g_0_x_xxxxy_xxyyy, g_0_x_xxxxy_xxyyz, g_0_x_xxxxy_xxyz, g_0_x_xxxxy_xxyzz, g_0_x_xxxxy_xxzz, g_0_x_xxxxy_xyyy, g_0_x_xxxxy_xyyyy, g_0_x_xxxxy_xyyyz, g_0_x_xxxxy_xyyz, g_0_x_xxxxy_xyyzz, g_0_x_xxxxy_xyzz, g_0_x_xxxxy_xyzzz, g_0_x_xxxxy_xzzz, g_0_x_xxxxy_yyyy, g_0_x_xxxxy_yyyyy, g_0_x_xxxxy_yyyyz, g_0_x_xxxxy_yyyz, g_0_x_xxxxy_yyyzz, g_0_x_xxxxy_yyzz, g_0_x_xxxxy_yyzzz, g_0_x_xxxxy_yzzz, g_0_x_xxxxy_yzzzz, g_0_x_xxxxy_zzzz, g_0_x_xxxxyy_xxxx, g_0_x_xxxxyy_xxxy, g_0_x_xxxxyy_xxxz, g_0_x_xxxxyy_xxyy, g_0_x_xxxxyy_xxyz, g_0_x_xxxxyy_xxzz, g_0_x_xxxxyy_xyyy, g_0_x_xxxxyy_xyyz, g_0_x_xxxxyy_xyzz, g_0_x_xxxxyy_xzzz, g_0_x_xxxxyy_yyyy, g_0_x_xxxxyy_yyyz, g_0_x_xxxxyy_yyzz, g_0_x_xxxxyy_yzzz, g_0_x_xxxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyy_xxxx[k] = -g_0_x_xxxxy_xxxx[k] * ab_y + g_0_x_xxxxy_xxxxy[k];

                g_0_x_xxxxyy_xxxy[k] = -g_0_x_xxxxy_xxxy[k] * ab_y + g_0_x_xxxxy_xxxyy[k];

                g_0_x_xxxxyy_xxxz[k] = -g_0_x_xxxxy_xxxz[k] * ab_y + g_0_x_xxxxy_xxxyz[k];

                g_0_x_xxxxyy_xxyy[k] = -g_0_x_xxxxy_xxyy[k] * ab_y + g_0_x_xxxxy_xxyyy[k];

                g_0_x_xxxxyy_xxyz[k] = -g_0_x_xxxxy_xxyz[k] * ab_y + g_0_x_xxxxy_xxyyz[k];

                g_0_x_xxxxyy_xxzz[k] = -g_0_x_xxxxy_xxzz[k] * ab_y + g_0_x_xxxxy_xxyzz[k];

                g_0_x_xxxxyy_xyyy[k] = -g_0_x_xxxxy_xyyy[k] * ab_y + g_0_x_xxxxy_xyyyy[k];

                g_0_x_xxxxyy_xyyz[k] = -g_0_x_xxxxy_xyyz[k] * ab_y + g_0_x_xxxxy_xyyyz[k];

                g_0_x_xxxxyy_xyzz[k] = -g_0_x_xxxxy_xyzz[k] * ab_y + g_0_x_xxxxy_xyyzz[k];

                g_0_x_xxxxyy_xzzz[k] = -g_0_x_xxxxy_xzzz[k] * ab_y + g_0_x_xxxxy_xyzzz[k];

                g_0_x_xxxxyy_yyyy[k] = -g_0_x_xxxxy_yyyy[k] * ab_y + g_0_x_xxxxy_yyyyy[k];

                g_0_x_xxxxyy_yyyz[k] = -g_0_x_xxxxy_yyyz[k] * ab_y + g_0_x_xxxxy_yyyyz[k];

                g_0_x_xxxxyy_yyzz[k] = -g_0_x_xxxxy_yyzz[k] * ab_y + g_0_x_xxxxy_yyyzz[k];

                g_0_x_xxxxyy_yzzz[k] = -g_0_x_xxxxy_yzzz[k] * ab_y + g_0_x_xxxxy_yyzzz[k];

                g_0_x_xxxxyy_zzzz[k] = -g_0_x_xxxxy_zzzz[k] * ab_y + g_0_x_xxxxy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyz_xxxx = cbuffer.data(ig_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxxy = cbuffer.data(ig_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxxz = cbuffer.data(ig_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxyy = cbuffer.data(ig_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxyz = cbuffer.data(ig_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxzz = cbuffer.data(ig_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyyy = cbuffer.data(ig_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyyz = cbuffer.data(ig_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyzz = cbuffer.data(ig_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xzzz = cbuffer.data(ig_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyyy = cbuffer.data(ig_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyyz = cbuffer.data(ig_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyzz = cbuffer.data(ig_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yzzz = cbuffer.data(ig_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxxyz_zzzz = cbuffer.data(ig_geom_01_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyz_xxxx, g_0_x_xxxxyz_xxxy, g_0_x_xxxxyz_xxxz, g_0_x_xxxxyz_xxyy, g_0_x_xxxxyz_xxyz, g_0_x_xxxxyz_xxzz, g_0_x_xxxxyz_xyyy, g_0_x_xxxxyz_xyyz, g_0_x_xxxxyz_xyzz, g_0_x_xxxxyz_xzzz, g_0_x_xxxxyz_yyyy, g_0_x_xxxxyz_yyyz, g_0_x_xxxxyz_yyzz, g_0_x_xxxxyz_yzzz, g_0_x_xxxxyz_zzzz, g_0_x_xxxxz_xxxx, g_0_x_xxxxz_xxxxy, g_0_x_xxxxz_xxxy, g_0_x_xxxxz_xxxyy, g_0_x_xxxxz_xxxyz, g_0_x_xxxxz_xxxz, g_0_x_xxxxz_xxyy, g_0_x_xxxxz_xxyyy, g_0_x_xxxxz_xxyyz, g_0_x_xxxxz_xxyz, g_0_x_xxxxz_xxyzz, g_0_x_xxxxz_xxzz, g_0_x_xxxxz_xyyy, g_0_x_xxxxz_xyyyy, g_0_x_xxxxz_xyyyz, g_0_x_xxxxz_xyyz, g_0_x_xxxxz_xyyzz, g_0_x_xxxxz_xyzz, g_0_x_xxxxz_xyzzz, g_0_x_xxxxz_xzzz, g_0_x_xxxxz_yyyy, g_0_x_xxxxz_yyyyy, g_0_x_xxxxz_yyyyz, g_0_x_xxxxz_yyyz, g_0_x_xxxxz_yyyzz, g_0_x_xxxxz_yyzz, g_0_x_xxxxz_yyzzz, g_0_x_xxxxz_yzzz, g_0_x_xxxxz_yzzzz, g_0_x_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyz_xxxx[k] = -g_0_x_xxxxz_xxxx[k] * ab_y + g_0_x_xxxxz_xxxxy[k];

                g_0_x_xxxxyz_xxxy[k] = -g_0_x_xxxxz_xxxy[k] * ab_y + g_0_x_xxxxz_xxxyy[k];

                g_0_x_xxxxyz_xxxz[k] = -g_0_x_xxxxz_xxxz[k] * ab_y + g_0_x_xxxxz_xxxyz[k];

                g_0_x_xxxxyz_xxyy[k] = -g_0_x_xxxxz_xxyy[k] * ab_y + g_0_x_xxxxz_xxyyy[k];

                g_0_x_xxxxyz_xxyz[k] = -g_0_x_xxxxz_xxyz[k] * ab_y + g_0_x_xxxxz_xxyyz[k];

                g_0_x_xxxxyz_xxzz[k] = -g_0_x_xxxxz_xxzz[k] * ab_y + g_0_x_xxxxz_xxyzz[k];

                g_0_x_xxxxyz_xyyy[k] = -g_0_x_xxxxz_xyyy[k] * ab_y + g_0_x_xxxxz_xyyyy[k];

                g_0_x_xxxxyz_xyyz[k] = -g_0_x_xxxxz_xyyz[k] * ab_y + g_0_x_xxxxz_xyyyz[k];

                g_0_x_xxxxyz_xyzz[k] = -g_0_x_xxxxz_xyzz[k] * ab_y + g_0_x_xxxxz_xyyzz[k];

                g_0_x_xxxxyz_xzzz[k] = -g_0_x_xxxxz_xzzz[k] * ab_y + g_0_x_xxxxz_xyzzz[k];

                g_0_x_xxxxyz_yyyy[k] = -g_0_x_xxxxz_yyyy[k] * ab_y + g_0_x_xxxxz_yyyyy[k];

                g_0_x_xxxxyz_yyyz[k] = -g_0_x_xxxxz_yyyz[k] * ab_y + g_0_x_xxxxz_yyyyz[k];

                g_0_x_xxxxyz_yyzz[k] = -g_0_x_xxxxz_yyzz[k] * ab_y + g_0_x_xxxxz_yyyzz[k];

                g_0_x_xxxxyz_yzzz[k] = -g_0_x_xxxxz_yzzz[k] * ab_y + g_0_x_xxxxz_yyzzz[k];

                g_0_x_xxxxyz_zzzz[k] = -g_0_x_xxxxz_zzzz[k] * ab_y + g_0_x_xxxxz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxzz_xxxx = cbuffer.data(ig_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxxy = cbuffer.data(ig_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxxz = cbuffer.data(ig_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxyy = cbuffer.data(ig_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxyz = cbuffer.data(ig_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxzz = cbuffer.data(ig_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyyy = cbuffer.data(ig_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyyz = cbuffer.data(ig_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyzz = cbuffer.data(ig_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xzzz = cbuffer.data(ig_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyyy = cbuffer.data(ig_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyyz = cbuffer.data(ig_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyzz = cbuffer.data(ig_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yzzz = cbuffer.data(ig_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxxzz_zzzz = cbuffer.data(ig_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxz_xxxx, g_0_x_xxxxz_xxxxz, g_0_x_xxxxz_xxxy, g_0_x_xxxxz_xxxyz, g_0_x_xxxxz_xxxz, g_0_x_xxxxz_xxxzz, g_0_x_xxxxz_xxyy, g_0_x_xxxxz_xxyyz, g_0_x_xxxxz_xxyz, g_0_x_xxxxz_xxyzz, g_0_x_xxxxz_xxzz, g_0_x_xxxxz_xxzzz, g_0_x_xxxxz_xyyy, g_0_x_xxxxz_xyyyz, g_0_x_xxxxz_xyyz, g_0_x_xxxxz_xyyzz, g_0_x_xxxxz_xyzz, g_0_x_xxxxz_xyzzz, g_0_x_xxxxz_xzzz, g_0_x_xxxxz_xzzzz, g_0_x_xxxxz_yyyy, g_0_x_xxxxz_yyyyz, g_0_x_xxxxz_yyyz, g_0_x_xxxxz_yyyzz, g_0_x_xxxxz_yyzz, g_0_x_xxxxz_yyzzz, g_0_x_xxxxz_yzzz, g_0_x_xxxxz_yzzzz, g_0_x_xxxxz_zzzz, g_0_x_xxxxz_zzzzz, g_0_x_xxxxzz_xxxx, g_0_x_xxxxzz_xxxy, g_0_x_xxxxzz_xxxz, g_0_x_xxxxzz_xxyy, g_0_x_xxxxzz_xxyz, g_0_x_xxxxzz_xxzz, g_0_x_xxxxzz_xyyy, g_0_x_xxxxzz_xyyz, g_0_x_xxxxzz_xyzz, g_0_x_xxxxzz_xzzz, g_0_x_xxxxzz_yyyy, g_0_x_xxxxzz_yyyz, g_0_x_xxxxzz_yyzz, g_0_x_xxxxzz_yzzz, g_0_x_xxxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxzz_xxxx[k] = -g_0_x_xxxxz_xxxx[k] * ab_z + g_0_x_xxxxz_xxxxz[k];

                g_0_x_xxxxzz_xxxy[k] = -g_0_x_xxxxz_xxxy[k] * ab_z + g_0_x_xxxxz_xxxyz[k];

                g_0_x_xxxxzz_xxxz[k] = -g_0_x_xxxxz_xxxz[k] * ab_z + g_0_x_xxxxz_xxxzz[k];

                g_0_x_xxxxzz_xxyy[k] = -g_0_x_xxxxz_xxyy[k] * ab_z + g_0_x_xxxxz_xxyyz[k];

                g_0_x_xxxxzz_xxyz[k] = -g_0_x_xxxxz_xxyz[k] * ab_z + g_0_x_xxxxz_xxyzz[k];

                g_0_x_xxxxzz_xxzz[k] = -g_0_x_xxxxz_xxzz[k] * ab_z + g_0_x_xxxxz_xxzzz[k];

                g_0_x_xxxxzz_xyyy[k] = -g_0_x_xxxxz_xyyy[k] * ab_z + g_0_x_xxxxz_xyyyz[k];

                g_0_x_xxxxzz_xyyz[k] = -g_0_x_xxxxz_xyyz[k] * ab_z + g_0_x_xxxxz_xyyzz[k];

                g_0_x_xxxxzz_xyzz[k] = -g_0_x_xxxxz_xyzz[k] * ab_z + g_0_x_xxxxz_xyzzz[k];

                g_0_x_xxxxzz_xzzz[k] = -g_0_x_xxxxz_xzzz[k] * ab_z + g_0_x_xxxxz_xzzzz[k];

                g_0_x_xxxxzz_yyyy[k] = -g_0_x_xxxxz_yyyy[k] * ab_z + g_0_x_xxxxz_yyyyz[k];

                g_0_x_xxxxzz_yyyz[k] = -g_0_x_xxxxz_yyyz[k] * ab_z + g_0_x_xxxxz_yyyzz[k];

                g_0_x_xxxxzz_yyzz[k] = -g_0_x_xxxxz_yyzz[k] * ab_z + g_0_x_xxxxz_yyzzz[k];

                g_0_x_xxxxzz_yzzz[k] = -g_0_x_xxxxz_yzzz[k] * ab_z + g_0_x_xxxxz_yzzzz[k];

                g_0_x_xxxxzz_zzzz[k] = -g_0_x_xxxxz_zzzz[k] * ab_z + g_0_x_xxxxz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyy_xxxx = cbuffer.data(ig_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxxy = cbuffer.data(ig_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxxz = cbuffer.data(ig_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxyy = cbuffer.data(ig_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxyz = cbuffer.data(ig_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxzz = cbuffer.data(ig_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyyy = cbuffer.data(ig_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyyz = cbuffer.data(ig_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyzz = cbuffer.data(ig_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xzzz = cbuffer.data(ig_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyyy = cbuffer.data(ig_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyyz = cbuffer.data(ig_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyzz = cbuffer.data(ig_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yzzz = cbuffer.data(ig_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxxyyy_zzzz = cbuffer.data(ig_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyy_xxxx, g_0_x_xxxyy_xxxxy, g_0_x_xxxyy_xxxy, g_0_x_xxxyy_xxxyy, g_0_x_xxxyy_xxxyz, g_0_x_xxxyy_xxxz, g_0_x_xxxyy_xxyy, g_0_x_xxxyy_xxyyy, g_0_x_xxxyy_xxyyz, g_0_x_xxxyy_xxyz, g_0_x_xxxyy_xxyzz, g_0_x_xxxyy_xxzz, g_0_x_xxxyy_xyyy, g_0_x_xxxyy_xyyyy, g_0_x_xxxyy_xyyyz, g_0_x_xxxyy_xyyz, g_0_x_xxxyy_xyyzz, g_0_x_xxxyy_xyzz, g_0_x_xxxyy_xyzzz, g_0_x_xxxyy_xzzz, g_0_x_xxxyy_yyyy, g_0_x_xxxyy_yyyyy, g_0_x_xxxyy_yyyyz, g_0_x_xxxyy_yyyz, g_0_x_xxxyy_yyyzz, g_0_x_xxxyy_yyzz, g_0_x_xxxyy_yyzzz, g_0_x_xxxyy_yzzz, g_0_x_xxxyy_yzzzz, g_0_x_xxxyy_zzzz, g_0_x_xxxyyy_xxxx, g_0_x_xxxyyy_xxxy, g_0_x_xxxyyy_xxxz, g_0_x_xxxyyy_xxyy, g_0_x_xxxyyy_xxyz, g_0_x_xxxyyy_xxzz, g_0_x_xxxyyy_xyyy, g_0_x_xxxyyy_xyyz, g_0_x_xxxyyy_xyzz, g_0_x_xxxyyy_xzzz, g_0_x_xxxyyy_yyyy, g_0_x_xxxyyy_yyyz, g_0_x_xxxyyy_yyzz, g_0_x_xxxyyy_yzzz, g_0_x_xxxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyy_xxxx[k] = -g_0_x_xxxyy_xxxx[k] * ab_y + g_0_x_xxxyy_xxxxy[k];

                g_0_x_xxxyyy_xxxy[k] = -g_0_x_xxxyy_xxxy[k] * ab_y + g_0_x_xxxyy_xxxyy[k];

                g_0_x_xxxyyy_xxxz[k] = -g_0_x_xxxyy_xxxz[k] * ab_y + g_0_x_xxxyy_xxxyz[k];

                g_0_x_xxxyyy_xxyy[k] = -g_0_x_xxxyy_xxyy[k] * ab_y + g_0_x_xxxyy_xxyyy[k];

                g_0_x_xxxyyy_xxyz[k] = -g_0_x_xxxyy_xxyz[k] * ab_y + g_0_x_xxxyy_xxyyz[k];

                g_0_x_xxxyyy_xxzz[k] = -g_0_x_xxxyy_xxzz[k] * ab_y + g_0_x_xxxyy_xxyzz[k];

                g_0_x_xxxyyy_xyyy[k] = -g_0_x_xxxyy_xyyy[k] * ab_y + g_0_x_xxxyy_xyyyy[k];

                g_0_x_xxxyyy_xyyz[k] = -g_0_x_xxxyy_xyyz[k] * ab_y + g_0_x_xxxyy_xyyyz[k];

                g_0_x_xxxyyy_xyzz[k] = -g_0_x_xxxyy_xyzz[k] * ab_y + g_0_x_xxxyy_xyyzz[k];

                g_0_x_xxxyyy_xzzz[k] = -g_0_x_xxxyy_xzzz[k] * ab_y + g_0_x_xxxyy_xyzzz[k];

                g_0_x_xxxyyy_yyyy[k] = -g_0_x_xxxyy_yyyy[k] * ab_y + g_0_x_xxxyy_yyyyy[k];

                g_0_x_xxxyyy_yyyz[k] = -g_0_x_xxxyy_yyyz[k] * ab_y + g_0_x_xxxyy_yyyyz[k];

                g_0_x_xxxyyy_yyzz[k] = -g_0_x_xxxyy_yyzz[k] * ab_y + g_0_x_xxxyy_yyyzz[k];

                g_0_x_xxxyyy_yzzz[k] = -g_0_x_xxxyy_yzzz[k] * ab_y + g_0_x_xxxyy_yyzzz[k];

                g_0_x_xxxyyy_zzzz[k] = -g_0_x_xxxyy_zzzz[k] * ab_y + g_0_x_xxxyy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyz_xxxx = cbuffer.data(ig_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxxy = cbuffer.data(ig_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxxz = cbuffer.data(ig_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxyy = cbuffer.data(ig_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxyz = cbuffer.data(ig_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxzz = cbuffer.data(ig_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyyy = cbuffer.data(ig_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyyz = cbuffer.data(ig_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyzz = cbuffer.data(ig_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xzzz = cbuffer.data(ig_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyyy = cbuffer.data(ig_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyyz = cbuffer.data(ig_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyzz = cbuffer.data(ig_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yzzz = cbuffer.data(ig_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxxyyz_zzzz = cbuffer.data(ig_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyz_xxxx, g_0_x_xxxyyz_xxxy, g_0_x_xxxyyz_xxxz, g_0_x_xxxyyz_xxyy, g_0_x_xxxyyz_xxyz, g_0_x_xxxyyz_xxzz, g_0_x_xxxyyz_xyyy, g_0_x_xxxyyz_xyyz, g_0_x_xxxyyz_xyzz, g_0_x_xxxyyz_xzzz, g_0_x_xxxyyz_yyyy, g_0_x_xxxyyz_yyyz, g_0_x_xxxyyz_yyzz, g_0_x_xxxyyz_yzzz, g_0_x_xxxyyz_zzzz, g_0_x_xxxyz_xxxx, g_0_x_xxxyz_xxxxy, g_0_x_xxxyz_xxxy, g_0_x_xxxyz_xxxyy, g_0_x_xxxyz_xxxyz, g_0_x_xxxyz_xxxz, g_0_x_xxxyz_xxyy, g_0_x_xxxyz_xxyyy, g_0_x_xxxyz_xxyyz, g_0_x_xxxyz_xxyz, g_0_x_xxxyz_xxyzz, g_0_x_xxxyz_xxzz, g_0_x_xxxyz_xyyy, g_0_x_xxxyz_xyyyy, g_0_x_xxxyz_xyyyz, g_0_x_xxxyz_xyyz, g_0_x_xxxyz_xyyzz, g_0_x_xxxyz_xyzz, g_0_x_xxxyz_xyzzz, g_0_x_xxxyz_xzzz, g_0_x_xxxyz_yyyy, g_0_x_xxxyz_yyyyy, g_0_x_xxxyz_yyyyz, g_0_x_xxxyz_yyyz, g_0_x_xxxyz_yyyzz, g_0_x_xxxyz_yyzz, g_0_x_xxxyz_yyzzz, g_0_x_xxxyz_yzzz, g_0_x_xxxyz_yzzzz, g_0_x_xxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyz_xxxx[k] = -g_0_x_xxxyz_xxxx[k] * ab_y + g_0_x_xxxyz_xxxxy[k];

                g_0_x_xxxyyz_xxxy[k] = -g_0_x_xxxyz_xxxy[k] * ab_y + g_0_x_xxxyz_xxxyy[k];

                g_0_x_xxxyyz_xxxz[k] = -g_0_x_xxxyz_xxxz[k] * ab_y + g_0_x_xxxyz_xxxyz[k];

                g_0_x_xxxyyz_xxyy[k] = -g_0_x_xxxyz_xxyy[k] * ab_y + g_0_x_xxxyz_xxyyy[k];

                g_0_x_xxxyyz_xxyz[k] = -g_0_x_xxxyz_xxyz[k] * ab_y + g_0_x_xxxyz_xxyyz[k];

                g_0_x_xxxyyz_xxzz[k] = -g_0_x_xxxyz_xxzz[k] * ab_y + g_0_x_xxxyz_xxyzz[k];

                g_0_x_xxxyyz_xyyy[k] = -g_0_x_xxxyz_xyyy[k] * ab_y + g_0_x_xxxyz_xyyyy[k];

                g_0_x_xxxyyz_xyyz[k] = -g_0_x_xxxyz_xyyz[k] * ab_y + g_0_x_xxxyz_xyyyz[k];

                g_0_x_xxxyyz_xyzz[k] = -g_0_x_xxxyz_xyzz[k] * ab_y + g_0_x_xxxyz_xyyzz[k];

                g_0_x_xxxyyz_xzzz[k] = -g_0_x_xxxyz_xzzz[k] * ab_y + g_0_x_xxxyz_xyzzz[k];

                g_0_x_xxxyyz_yyyy[k] = -g_0_x_xxxyz_yyyy[k] * ab_y + g_0_x_xxxyz_yyyyy[k];

                g_0_x_xxxyyz_yyyz[k] = -g_0_x_xxxyz_yyyz[k] * ab_y + g_0_x_xxxyz_yyyyz[k];

                g_0_x_xxxyyz_yyzz[k] = -g_0_x_xxxyz_yyzz[k] * ab_y + g_0_x_xxxyz_yyyzz[k];

                g_0_x_xxxyyz_yzzz[k] = -g_0_x_xxxyz_yzzz[k] * ab_y + g_0_x_xxxyz_yyzzz[k];

                g_0_x_xxxyyz_zzzz[k] = -g_0_x_xxxyz_zzzz[k] * ab_y + g_0_x_xxxyz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyzz_xxxx = cbuffer.data(ig_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxxy = cbuffer.data(ig_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxxz = cbuffer.data(ig_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxyy = cbuffer.data(ig_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxyz = cbuffer.data(ig_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxzz = cbuffer.data(ig_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyyy = cbuffer.data(ig_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyyz = cbuffer.data(ig_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyzz = cbuffer.data(ig_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xzzz = cbuffer.data(ig_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyyy = cbuffer.data(ig_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyyz = cbuffer.data(ig_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyzz = cbuffer.data(ig_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yzzz = cbuffer.data(ig_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxxyzz_zzzz = cbuffer.data(ig_geom_01_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyzz_xxxx, g_0_x_xxxyzz_xxxy, g_0_x_xxxyzz_xxxz, g_0_x_xxxyzz_xxyy, g_0_x_xxxyzz_xxyz, g_0_x_xxxyzz_xxzz, g_0_x_xxxyzz_xyyy, g_0_x_xxxyzz_xyyz, g_0_x_xxxyzz_xyzz, g_0_x_xxxyzz_xzzz, g_0_x_xxxyzz_yyyy, g_0_x_xxxyzz_yyyz, g_0_x_xxxyzz_yyzz, g_0_x_xxxyzz_yzzz, g_0_x_xxxyzz_zzzz, g_0_x_xxxzz_xxxx, g_0_x_xxxzz_xxxxy, g_0_x_xxxzz_xxxy, g_0_x_xxxzz_xxxyy, g_0_x_xxxzz_xxxyz, g_0_x_xxxzz_xxxz, g_0_x_xxxzz_xxyy, g_0_x_xxxzz_xxyyy, g_0_x_xxxzz_xxyyz, g_0_x_xxxzz_xxyz, g_0_x_xxxzz_xxyzz, g_0_x_xxxzz_xxzz, g_0_x_xxxzz_xyyy, g_0_x_xxxzz_xyyyy, g_0_x_xxxzz_xyyyz, g_0_x_xxxzz_xyyz, g_0_x_xxxzz_xyyzz, g_0_x_xxxzz_xyzz, g_0_x_xxxzz_xyzzz, g_0_x_xxxzz_xzzz, g_0_x_xxxzz_yyyy, g_0_x_xxxzz_yyyyy, g_0_x_xxxzz_yyyyz, g_0_x_xxxzz_yyyz, g_0_x_xxxzz_yyyzz, g_0_x_xxxzz_yyzz, g_0_x_xxxzz_yyzzz, g_0_x_xxxzz_yzzz, g_0_x_xxxzz_yzzzz, g_0_x_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyzz_xxxx[k] = -g_0_x_xxxzz_xxxx[k] * ab_y + g_0_x_xxxzz_xxxxy[k];

                g_0_x_xxxyzz_xxxy[k] = -g_0_x_xxxzz_xxxy[k] * ab_y + g_0_x_xxxzz_xxxyy[k];

                g_0_x_xxxyzz_xxxz[k] = -g_0_x_xxxzz_xxxz[k] * ab_y + g_0_x_xxxzz_xxxyz[k];

                g_0_x_xxxyzz_xxyy[k] = -g_0_x_xxxzz_xxyy[k] * ab_y + g_0_x_xxxzz_xxyyy[k];

                g_0_x_xxxyzz_xxyz[k] = -g_0_x_xxxzz_xxyz[k] * ab_y + g_0_x_xxxzz_xxyyz[k];

                g_0_x_xxxyzz_xxzz[k] = -g_0_x_xxxzz_xxzz[k] * ab_y + g_0_x_xxxzz_xxyzz[k];

                g_0_x_xxxyzz_xyyy[k] = -g_0_x_xxxzz_xyyy[k] * ab_y + g_0_x_xxxzz_xyyyy[k];

                g_0_x_xxxyzz_xyyz[k] = -g_0_x_xxxzz_xyyz[k] * ab_y + g_0_x_xxxzz_xyyyz[k];

                g_0_x_xxxyzz_xyzz[k] = -g_0_x_xxxzz_xyzz[k] * ab_y + g_0_x_xxxzz_xyyzz[k];

                g_0_x_xxxyzz_xzzz[k] = -g_0_x_xxxzz_xzzz[k] * ab_y + g_0_x_xxxzz_xyzzz[k];

                g_0_x_xxxyzz_yyyy[k] = -g_0_x_xxxzz_yyyy[k] * ab_y + g_0_x_xxxzz_yyyyy[k];

                g_0_x_xxxyzz_yyyz[k] = -g_0_x_xxxzz_yyyz[k] * ab_y + g_0_x_xxxzz_yyyyz[k];

                g_0_x_xxxyzz_yyzz[k] = -g_0_x_xxxzz_yyzz[k] * ab_y + g_0_x_xxxzz_yyyzz[k];

                g_0_x_xxxyzz_yzzz[k] = -g_0_x_xxxzz_yzzz[k] * ab_y + g_0_x_xxxzz_yyzzz[k];

                g_0_x_xxxyzz_zzzz[k] = -g_0_x_xxxzz_zzzz[k] * ab_y + g_0_x_xxxzz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzzz_xxxx = cbuffer.data(ig_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxxy = cbuffer.data(ig_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxxz = cbuffer.data(ig_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxyy = cbuffer.data(ig_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxyz = cbuffer.data(ig_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxzz = cbuffer.data(ig_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyyy = cbuffer.data(ig_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyyz = cbuffer.data(ig_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyzz = cbuffer.data(ig_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xzzz = cbuffer.data(ig_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyyy = cbuffer.data(ig_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyyz = cbuffer.data(ig_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyzz = cbuffer.data(ig_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yzzz = cbuffer.data(ig_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxxzzz_zzzz = cbuffer.data(ig_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxzz_xxxx, g_0_x_xxxzz_xxxxz, g_0_x_xxxzz_xxxy, g_0_x_xxxzz_xxxyz, g_0_x_xxxzz_xxxz, g_0_x_xxxzz_xxxzz, g_0_x_xxxzz_xxyy, g_0_x_xxxzz_xxyyz, g_0_x_xxxzz_xxyz, g_0_x_xxxzz_xxyzz, g_0_x_xxxzz_xxzz, g_0_x_xxxzz_xxzzz, g_0_x_xxxzz_xyyy, g_0_x_xxxzz_xyyyz, g_0_x_xxxzz_xyyz, g_0_x_xxxzz_xyyzz, g_0_x_xxxzz_xyzz, g_0_x_xxxzz_xyzzz, g_0_x_xxxzz_xzzz, g_0_x_xxxzz_xzzzz, g_0_x_xxxzz_yyyy, g_0_x_xxxzz_yyyyz, g_0_x_xxxzz_yyyz, g_0_x_xxxzz_yyyzz, g_0_x_xxxzz_yyzz, g_0_x_xxxzz_yyzzz, g_0_x_xxxzz_yzzz, g_0_x_xxxzz_yzzzz, g_0_x_xxxzz_zzzz, g_0_x_xxxzz_zzzzz, g_0_x_xxxzzz_xxxx, g_0_x_xxxzzz_xxxy, g_0_x_xxxzzz_xxxz, g_0_x_xxxzzz_xxyy, g_0_x_xxxzzz_xxyz, g_0_x_xxxzzz_xxzz, g_0_x_xxxzzz_xyyy, g_0_x_xxxzzz_xyyz, g_0_x_xxxzzz_xyzz, g_0_x_xxxzzz_xzzz, g_0_x_xxxzzz_yyyy, g_0_x_xxxzzz_yyyz, g_0_x_xxxzzz_yyzz, g_0_x_xxxzzz_yzzz, g_0_x_xxxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzzz_xxxx[k] = -g_0_x_xxxzz_xxxx[k] * ab_z + g_0_x_xxxzz_xxxxz[k];

                g_0_x_xxxzzz_xxxy[k] = -g_0_x_xxxzz_xxxy[k] * ab_z + g_0_x_xxxzz_xxxyz[k];

                g_0_x_xxxzzz_xxxz[k] = -g_0_x_xxxzz_xxxz[k] * ab_z + g_0_x_xxxzz_xxxzz[k];

                g_0_x_xxxzzz_xxyy[k] = -g_0_x_xxxzz_xxyy[k] * ab_z + g_0_x_xxxzz_xxyyz[k];

                g_0_x_xxxzzz_xxyz[k] = -g_0_x_xxxzz_xxyz[k] * ab_z + g_0_x_xxxzz_xxyzz[k];

                g_0_x_xxxzzz_xxzz[k] = -g_0_x_xxxzz_xxzz[k] * ab_z + g_0_x_xxxzz_xxzzz[k];

                g_0_x_xxxzzz_xyyy[k] = -g_0_x_xxxzz_xyyy[k] * ab_z + g_0_x_xxxzz_xyyyz[k];

                g_0_x_xxxzzz_xyyz[k] = -g_0_x_xxxzz_xyyz[k] * ab_z + g_0_x_xxxzz_xyyzz[k];

                g_0_x_xxxzzz_xyzz[k] = -g_0_x_xxxzz_xyzz[k] * ab_z + g_0_x_xxxzz_xyzzz[k];

                g_0_x_xxxzzz_xzzz[k] = -g_0_x_xxxzz_xzzz[k] * ab_z + g_0_x_xxxzz_xzzzz[k];

                g_0_x_xxxzzz_yyyy[k] = -g_0_x_xxxzz_yyyy[k] * ab_z + g_0_x_xxxzz_yyyyz[k];

                g_0_x_xxxzzz_yyyz[k] = -g_0_x_xxxzz_yyyz[k] * ab_z + g_0_x_xxxzz_yyyzz[k];

                g_0_x_xxxzzz_yyzz[k] = -g_0_x_xxxzz_yyzz[k] * ab_z + g_0_x_xxxzz_yyzzz[k];

                g_0_x_xxxzzz_yzzz[k] = -g_0_x_xxxzz_yzzz[k] * ab_z + g_0_x_xxxzz_yzzzz[k];

                g_0_x_xxxzzz_zzzz[k] = -g_0_x_xxxzz_zzzz[k] * ab_z + g_0_x_xxxzz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyy_xxxx = cbuffer.data(ig_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxxy = cbuffer.data(ig_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxxz = cbuffer.data(ig_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxyy = cbuffer.data(ig_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxyz = cbuffer.data(ig_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxzz = cbuffer.data(ig_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyyy = cbuffer.data(ig_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyyz = cbuffer.data(ig_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyzz = cbuffer.data(ig_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xzzz = cbuffer.data(ig_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyyy = cbuffer.data(ig_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyyz = cbuffer.data(ig_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyzz = cbuffer.data(ig_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yzzz = cbuffer.data(ig_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xxyyyy_zzzz = cbuffer.data(ig_geom_01_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyy_xxxx, g_0_x_xxyyy_xxxxy, g_0_x_xxyyy_xxxy, g_0_x_xxyyy_xxxyy, g_0_x_xxyyy_xxxyz, g_0_x_xxyyy_xxxz, g_0_x_xxyyy_xxyy, g_0_x_xxyyy_xxyyy, g_0_x_xxyyy_xxyyz, g_0_x_xxyyy_xxyz, g_0_x_xxyyy_xxyzz, g_0_x_xxyyy_xxzz, g_0_x_xxyyy_xyyy, g_0_x_xxyyy_xyyyy, g_0_x_xxyyy_xyyyz, g_0_x_xxyyy_xyyz, g_0_x_xxyyy_xyyzz, g_0_x_xxyyy_xyzz, g_0_x_xxyyy_xyzzz, g_0_x_xxyyy_xzzz, g_0_x_xxyyy_yyyy, g_0_x_xxyyy_yyyyy, g_0_x_xxyyy_yyyyz, g_0_x_xxyyy_yyyz, g_0_x_xxyyy_yyyzz, g_0_x_xxyyy_yyzz, g_0_x_xxyyy_yyzzz, g_0_x_xxyyy_yzzz, g_0_x_xxyyy_yzzzz, g_0_x_xxyyy_zzzz, g_0_x_xxyyyy_xxxx, g_0_x_xxyyyy_xxxy, g_0_x_xxyyyy_xxxz, g_0_x_xxyyyy_xxyy, g_0_x_xxyyyy_xxyz, g_0_x_xxyyyy_xxzz, g_0_x_xxyyyy_xyyy, g_0_x_xxyyyy_xyyz, g_0_x_xxyyyy_xyzz, g_0_x_xxyyyy_xzzz, g_0_x_xxyyyy_yyyy, g_0_x_xxyyyy_yyyz, g_0_x_xxyyyy_yyzz, g_0_x_xxyyyy_yzzz, g_0_x_xxyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyy_xxxx[k] = -g_0_x_xxyyy_xxxx[k] * ab_y + g_0_x_xxyyy_xxxxy[k];

                g_0_x_xxyyyy_xxxy[k] = -g_0_x_xxyyy_xxxy[k] * ab_y + g_0_x_xxyyy_xxxyy[k];

                g_0_x_xxyyyy_xxxz[k] = -g_0_x_xxyyy_xxxz[k] * ab_y + g_0_x_xxyyy_xxxyz[k];

                g_0_x_xxyyyy_xxyy[k] = -g_0_x_xxyyy_xxyy[k] * ab_y + g_0_x_xxyyy_xxyyy[k];

                g_0_x_xxyyyy_xxyz[k] = -g_0_x_xxyyy_xxyz[k] * ab_y + g_0_x_xxyyy_xxyyz[k];

                g_0_x_xxyyyy_xxzz[k] = -g_0_x_xxyyy_xxzz[k] * ab_y + g_0_x_xxyyy_xxyzz[k];

                g_0_x_xxyyyy_xyyy[k] = -g_0_x_xxyyy_xyyy[k] * ab_y + g_0_x_xxyyy_xyyyy[k];

                g_0_x_xxyyyy_xyyz[k] = -g_0_x_xxyyy_xyyz[k] * ab_y + g_0_x_xxyyy_xyyyz[k];

                g_0_x_xxyyyy_xyzz[k] = -g_0_x_xxyyy_xyzz[k] * ab_y + g_0_x_xxyyy_xyyzz[k];

                g_0_x_xxyyyy_xzzz[k] = -g_0_x_xxyyy_xzzz[k] * ab_y + g_0_x_xxyyy_xyzzz[k];

                g_0_x_xxyyyy_yyyy[k] = -g_0_x_xxyyy_yyyy[k] * ab_y + g_0_x_xxyyy_yyyyy[k];

                g_0_x_xxyyyy_yyyz[k] = -g_0_x_xxyyy_yyyz[k] * ab_y + g_0_x_xxyyy_yyyyz[k];

                g_0_x_xxyyyy_yyzz[k] = -g_0_x_xxyyy_yyzz[k] * ab_y + g_0_x_xxyyy_yyyzz[k];

                g_0_x_xxyyyy_yzzz[k] = -g_0_x_xxyyy_yzzz[k] * ab_y + g_0_x_xxyyy_yyzzz[k];

                g_0_x_xxyyyy_zzzz[k] = -g_0_x_xxyyy_zzzz[k] * ab_y + g_0_x_xxyyy_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyz_xxxx = cbuffer.data(ig_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxxy = cbuffer.data(ig_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxxz = cbuffer.data(ig_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxyy = cbuffer.data(ig_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxyz = cbuffer.data(ig_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxzz = cbuffer.data(ig_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyyy = cbuffer.data(ig_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyyz = cbuffer.data(ig_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyzz = cbuffer.data(ig_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xzzz = cbuffer.data(ig_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyyy = cbuffer.data(ig_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyyz = cbuffer.data(ig_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyzz = cbuffer.data(ig_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yzzz = cbuffer.data(ig_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xxyyyz_zzzz = cbuffer.data(ig_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyz_xxxx, g_0_x_xxyyyz_xxxy, g_0_x_xxyyyz_xxxz, g_0_x_xxyyyz_xxyy, g_0_x_xxyyyz_xxyz, g_0_x_xxyyyz_xxzz, g_0_x_xxyyyz_xyyy, g_0_x_xxyyyz_xyyz, g_0_x_xxyyyz_xyzz, g_0_x_xxyyyz_xzzz, g_0_x_xxyyyz_yyyy, g_0_x_xxyyyz_yyyz, g_0_x_xxyyyz_yyzz, g_0_x_xxyyyz_yzzz, g_0_x_xxyyyz_zzzz, g_0_x_xxyyz_xxxx, g_0_x_xxyyz_xxxxy, g_0_x_xxyyz_xxxy, g_0_x_xxyyz_xxxyy, g_0_x_xxyyz_xxxyz, g_0_x_xxyyz_xxxz, g_0_x_xxyyz_xxyy, g_0_x_xxyyz_xxyyy, g_0_x_xxyyz_xxyyz, g_0_x_xxyyz_xxyz, g_0_x_xxyyz_xxyzz, g_0_x_xxyyz_xxzz, g_0_x_xxyyz_xyyy, g_0_x_xxyyz_xyyyy, g_0_x_xxyyz_xyyyz, g_0_x_xxyyz_xyyz, g_0_x_xxyyz_xyyzz, g_0_x_xxyyz_xyzz, g_0_x_xxyyz_xyzzz, g_0_x_xxyyz_xzzz, g_0_x_xxyyz_yyyy, g_0_x_xxyyz_yyyyy, g_0_x_xxyyz_yyyyz, g_0_x_xxyyz_yyyz, g_0_x_xxyyz_yyyzz, g_0_x_xxyyz_yyzz, g_0_x_xxyyz_yyzzz, g_0_x_xxyyz_yzzz, g_0_x_xxyyz_yzzzz, g_0_x_xxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyz_xxxx[k] = -g_0_x_xxyyz_xxxx[k] * ab_y + g_0_x_xxyyz_xxxxy[k];

                g_0_x_xxyyyz_xxxy[k] = -g_0_x_xxyyz_xxxy[k] * ab_y + g_0_x_xxyyz_xxxyy[k];

                g_0_x_xxyyyz_xxxz[k] = -g_0_x_xxyyz_xxxz[k] * ab_y + g_0_x_xxyyz_xxxyz[k];

                g_0_x_xxyyyz_xxyy[k] = -g_0_x_xxyyz_xxyy[k] * ab_y + g_0_x_xxyyz_xxyyy[k];

                g_0_x_xxyyyz_xxyz[k] = -g_0_x_xxyyz_xxyz[k] * ab_y + g_0_x_xxyyz_xxyyz[k];

                g_0_x_xxyyyz_xxzz[k] = -g_0_x_xxyyz_xxzz[k] * ab_y + g_0_x_xxyyz_xxyzz[k];

                g_0_x_xxyyyz_xyyy[k] = -g_0_x_xxyyz_xyyy[k] * ab_y + g_0_x_xxyyz_xyyyy[k];

                g_0_x_xxyyyz_xyyz[k] = -g_0_x_xxyyz_xyyz[k] * ab_y + g_0_x_xxyyz_xyyyz[k];

                g_0_x_xxyyyz_xyzz[k] = -g_0_x_xxyyz_xyzz[k] * ab_y + g_0_x_xxyyz_xyyzz[k];

                g_0_x_xxyyyz_xzzz[k] = -g_0_x_xxyyz_xzzz[k] * ab_y + g_0_x_xxyyz_xyzzz[k];

                g_0_x_xxyyyz_yyyy[k] = -g_0_x_xxyyz_yyyy[k] * ab_y + g_0_x_xxyyz_yyyyy[k];

                g_0_x_xxyyyz_yyyz[k] = -g_0_x_xxyyz_yyyz[k] * ab_y + g_0_x_xxyyz_yyyyz[k];

                g_0_x_xxyyyz_yyzz[k] = -g_0_x_xxyyz_yyzz[k] * ab_y + g_0_x_xxyyz_yyyzz[k];

                g_0_x_xxyyyz_yzzz[k] = -g_0_x_xxyyz_yzzz[k] * ab_y + g_0_x_xxyyz_yyzzz[k];

                g_0_x_xxyyyz_zzzz[k] = -g_0_x_xxyyz_zzzz[k] * ab_y + g_0_x_xxyyz_yzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyzz_xxxx = cbuffer.data(ig_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxxy = cbuffer.data(ig_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxxz = cbuffer.data(ig_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxyy = cbuffer.data(ig_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxyz = cbuffer.data(ig_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxzz = cbuffer.data(ig_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyyy = cbuffer.data(ig_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyyz = cbuffer.data(ig_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyzz = cbuffer.data(ig_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xzzz = cbuffer.data(ig_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyyy = cbuffer.data(ig_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyyz = cbuffer.data(ig_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyzz = cbuffer.data(ig_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yzzz = cbuffer.data(ig_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xxyyzz_zzzz = cbuffer.data(ig_geom_01_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyzz_xxxx, g_0_x_xxyyzz_xxxy, g_0_x_xxyyzz_xxxz, g_0_x_xxyyzz_xxyy, g_0_x_xxyyzz_xxyz, g_0_x_xxyyzz_xxzz, g_0_x_xxyyzz_xyyy, g_0_x_xxyyzz_xyyz, g_0_x_xxyyzz_xyzz, g_0_x_xxyyzz_xzzz, g_0_x_xxyyzz_yyyy, g_0_x_xxyyzz_yyyz, g_0_x_xxyyzz_yyzz, g_0_x_xxyyzz_yzzz, g_0_x_xxyyzz_zzzz, g_0_x_xxyzz_xxxx, g_0_x_xxyzz_xxxxy, g_0_x_xxyzz_xxxy, g_0_x_xxyzz_xxxyy, g_0_x_xxyzz_xxxyz, g_0_x_xxyzz_xxxz, g_0_x_xxyzz_xxyy, g_0_x_xxyzz_xxyyy, g_0_x_xxyzz_xxyyz, g_0_x_xxyzz_xxyz, g_0_x_xxyzz_xxyzz, g_0_x_xxyzz_xxzz, g_0_x_xxyzz_xyyy, g_0_x_xxyzz_xyyyy, g_0_x_xxyzz_xyyyz, g_0_x_xxyzz_xyyz, g_0_x_xxyzz_xyyzz, g_0_x_xxyzz_xyzz, g_0_x_xxyzz_xyzzz, g_0_x_xxyzz_xzzz, g_0_x_xxyzz_yyyy, g_0_x_xxyzz_yyyyy, g_0_x_xxyzz_yyyyz, g_0_x_xxyzz_yyyz, g_0_x_xxyzz_yyyzz, g_0_x_xxyzz_yyzz, g_0_x_xxyzz_yyzzz, g_0_x_xxyzz_yzzz, g_0_x_xxyzz_yzzzz, g_0_x_xxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyzz_xxxx[k] = -g_0_x_xxyzz_xxxx[k] * ab_y + g_0_x_xxyzz_xxxxy[k];

                g_0_x_xxyyzz_xxxy[k] = -g_0_x_xxyzz_xxxy[k] * ab_y + g_0_x_xxyzz_xxxyy[k];

                g_0_x_xxyyzz_xxxz[k] = -g_0_x_xxyzz_xxxz[k] * ab_y + g_0_x_xxyzz_xxxyz[k];

                g_0_x_xxyyzz_xxyy[k] = -g_0_x_xxyzz_xxyy[k] * ab_y + g_0_x_xxyzz_xxyyy[k];

                g_0_x_xxyyzz_xxyz[k] = -g_0_x_xxyzz_xxyz[k] * ab_y + g_0_x_xxyzz_xxyyz[k];

                g_0_x_xxyyzz_xxzz[k] = -g_0_x_xxyzz_xxzz[k] * ab_y + g_0_x_xxyzz_xxyzz[k];

                g_0_x_xxyyzz_xyyy[k] = -g_0_x_xxyzz_xyyy[k] * ab_y + g_0_x_xxyzz_xyyyy[k];

                g_0_x_xxyyzz_xyyz[k] = -g_0_x_xxyzz_xyyz[k] * ab_y + g_0_x_xxyzz_xyyyz[k];

                g_0_x_xxyyzz_xyzz[k] = -g_0_x_xxyzz_xyzz[k] * ab_y + g_0_x_xxyzz_xyyzz[k];

                g_0_x_xxyyzz_xzzz[k] = -g_0_x_xxyzz_xzzz[k] * ab_y + g_0_x_xxyzz_xyzzz[k];

                g_0_x_xxyyzz_yyyy[k] = -g_0_x_xxyzz_yyyy[k] * ab_y + g_0_x_xxyzz_yyyyy[k];

                g_0_x_xxyyzz_yyyz[k] = -g_0_x_xxyzz_yyyz[k] * ab_y + g_0_x_xxyzz_yyyyz[k];

                g_0_x_xxyyzz_yyzz[k] = -g_0_x_xxyzz_yyzz[k] * ab_y + g_0_x_xxyzz_yyyzz[k];

                g_0_x_xxyyzz_yzzz[k] = -g_0_x_xxyzz_yzzz[k] * ab_y + g_0_x_xxyzz_yyzzz[k];

                g_0_x_xxyyzz_zzzz[k] = -g_0_x_xxyzz_zzzz[k] * ab_y + g_0_x_xxyzz_yzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzzz_xxxx = cbuffer.data(ig_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxxy = cbuffer.data(ig_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxxz = cbuffer.data(ig_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxyy = cbuffer.data(ig_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxyz = cbuffer.data(ig_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxzz = cbuffer.data(ig_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyyy = cbuffer.data(ig_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyyz = cbuffer.data(ig_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyzz = cbuffer.data(ig_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xzzz = cbuffer.data(ig_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyyy = cbuffer.data(ig_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyyz = cbuffer.data(ig_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyzz = cbuffer.data(ig_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yzzz = cbuffer.data(ig_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xxyzzz_zzzz = cbuffer.data(ig_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzzz_xxxx, g_0_x_xxyzzz_xxxy, g_0_x_xxyzzz_xxxz, g_0_x_xxyzzz_xxyy, g_0_x_xxyzzz_xxyz, g_0_x_xxyzzz_xxzz, g_0_x_xxyzzz_xyyy, g_0_x_xxyzzz_xyyz, g_0_x_xxyzzz_xyzz, g_0_x_xxyzzz_xzzz, g_0_x_xxyzzz_yyyy, g_0_x_xxyzzz_yyyz, g_0_x_xxyzzz_yyzz, g_0_x_xxyzzz_yzzz, g_0_x_xxyzzz_zzzz, g_0_x_xxzzz_xxxx, g_0_x_xxzzz_xxxxy, g_0_x_xxzzz_xxxy, g_0_x_xxzzz_xxxyy, g_0_x_xxzzz_xxxyz, g_0_x_xxzzz_xxxz, g_0_x_xxzzz_xxyy, g_0_x_xxzzz_xxyyy, g_0_x_xxzzz_xxyyz, g_0_x_xxzzz_xxyz, g_0_x_xxzzz_xxyzz, g_0_x_xxzzz_xxzz, g_0_x_xxzzz_xyyy, g_0_x_xxzzz_xyyyy, g_0_x_xxzzz_xyyyz, g_0_x_xxzzz_xyyz, g_0_x_xxzzz_xyyzz, g_0_x_xxzzz_xyzz, g_0_x_xxzzz_xyzzz, g_0_x_xxzzz_xzzz, g_0_x_xxzzz_yyyy, g_0_x_xxzzz_yyyyy, g_0_x_xxzzz_yyyyz, g_0_x_xxzzz_yyyz, g_0_x_xxzzz_yyyzz, g_0_x_xxzzz_yyzz, g_0_x_xxzzz_yyzzz, g_0_x_xxzzz_yzzz, g_0_x_xxzzz_yzzzz, g_0_x_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzzz_xxxx[k] = -g_0_x_xxzzz_xxxx[k] * ab_y + g_0_x_xxzzz_xxxxy[k];

                g_0_x_xxyzzz_xxxy[k] = -g_0_x_xxzzz_xxxy[k] * ab_y + g_0_x_xxzzz_xxxyy[k];

                g_0_x_xxyzzz_xxxz[k] = -g_0_x_xxzzz_xxxz[k] * ab_y + g_0_x_xxzzz_xxxyz[k];

                g_0_x_xxyzzz_xxyy[k] = -g_0_x_xxzzz_xxyy[k] * ab_y + g_0_x_xxzzz_xxyyy[k];

                g_0_x_xxyzzz_xxyz[k] = -g_0_x_xxzzz_xxyz[k] * ab_y + g_0_x_xxzzz_xxyyz[k];

                g_0_x_xxyzzz_xxzz[k] = -g_0_x_xxzzz_xxzz[k] * ab_y + g_0_x_xxzzz_xxyzz[k];

                g_0_x_xxyzzz_xyyy[k] = -g_0_x_xxzzz_xyyy[k] * ab_y + g_0_x_xxzzz_xyyyy[k];

                g_0_x_xxyzzz_xyyz[k] = -g_0_x_xxzzz_xyyz[k] * ab_y + g_0_x_xxzzz_xyyyz[k];

                g_0_x_xxyzzz_xyzz[k] = -g_0_x_xxzzz_xyzz[k] * ab_y + g_0_x_xxzzz_xyyzz[k];

                g_0_x_xxyzzz_xzzz[k] = -g_0_x_xxzzz_xzzz[k] * ab_y + g_0_x_xxzzz_xyzzz[k];

                g_0_x_xxyzzz_yyyy[k] = -g_0_x_xxzzz_yyyy[k] * ab_y + g_0_x_xxzzz_yyyyy[k];

                g_0_x_xxyzzz_yyyz[k] = -g_0_x_xxzzz_yyyz[k] * ab_y + g_0_x_xxzzz_yyyyz[k];

                g_0_x_xxyzzz_yyzz[k] = -g_0_x_xxzzz_yyzz[k] * ab_y + g_0_x_xxzzz_yyyzz[k];

                g_0_x_xxyzzz_yzzz[k] = -g_0_x_xxzzz_yzzz[k] * ab_y + g_0_x_xxzzz_yyzzz[k];

                g_0_x_xxyzzz_zzzz[k] = -g_0_x_xxzzz_zzzz[k] * ab_y + g_0_x_xxzzz_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzzz_xxxx = cbuffer.data(ig_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxxy = cbuffer.data(ig_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxxz = cbuffer.data(ig_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxyy = cbuffer.data(ig_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxyz = cbuffer.data(ig_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxzz = cbuffer.data(ig_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyyy = cbuffer.data(ig_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyyz = cbuffer.data(ig_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyzz = cbuffer.data(ig_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xzzz = cbuffer.data(ig_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyyy = cbuffer.data(ig_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyyz = cbuffer.data(ig_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyzz = cbuffer.data(ig_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yzzz = cbuffer.data(ig_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xxzzzz_zzzz = cbuffer.data(ig_geom_01_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzzz_xxxx, g_0_x_xxzzz_xxxxz, g_0_x_xxzzz_xxxy, g_0_x_xxzzz_xxxyz, g_0_x_xxzzz_xxxz, g_0_x_xxzzz_xxxzz, g_0_x_xxzzz_xxyy, g_0_x_xxzzz_xxyyz, g_0_x_xxzzz_xxyz, g_0_x_xxzzz_xxyzz, g_0_x_xxzzz_xxzz, g_0_x_xxzzz_xxzzz, g_0_x_xxzzz_xyyy, g_0_x_xxzzz_xyyyz, g_0_x_xxzzz_xyyz, g_0_x_xxzzz_xyyzz, g_0_x_xxzzz_xyzz, g_0_x_xxzzz_xyzzz, g_0_x_xxzzz_xzzz, g_0_x_xxzzz_xzzzz, g_0_x_xxzzz_yyyy, g_0_x_xxzzz_yyyyz, g_0_x_xxzzz_yyyz, g_0_x_xxzzz_yyyzz, g_0_x_xxzzz_yyzz, g_0_x_xxzzz_yyzzz, g_0_x_xxzzz_yzzz, g_0_x_xxzzz_yzzzz, g_0_x_xxzzz_zzzz, g_0_x_xxzzz_zzzzz, g_0_x_xxzzzz_xxxx, g_0_x_xxzzzz_xxxy, g_0_x_xxzzzz_xxxz, g_0_x_xxzzzz_xxyy, g_0_x_xxzzzz_xxyz, g_0_x_xxzzzz_xxzz, g_0_x_xxzzzz_xyyy, g_0_x_xxzzzz_xyyz, g_0_x_xxzzzz_xyzz, g_0_x_xxzzzz_xzzz, g_0_x_xxzzzz_yyyy, g_0_x_xxzzzz_yyyz, g_0_x_xxzzzz_yyzz, g_0_x_xxzzzz_yzzz, g_0_x_xxzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzzz_xxxx[k] = -g_0_x_xxzzz_xxxx[k] * ab_z + g_0_x_xxzzz_xxxxz[k];

                g_0_x_xxzzzz_xxxy[k] = -g_0_x_xxzzz_xxxy[k] * ab_z + g_0_x_xxzzz_xxxyz[k];

                g_0_x_xxzzzz_xxxz[k] = -g_0_x_xxzzz_xxxz[k] * ab_z + g_0_x_xxzzz_xxxzz[k];

                g_0_x_xxzzzz_xxyy[k] = -g_0_x_xxzzz_xxyy[k] * ab_z + g_0_x_xxzzz_xxyyz[k];

                g_0_x_xxzzzz_xxyz[k] = -g_0_x_xxzzz_xxyz[k] * ab_z + g_0_x_xxzzz_xxyzz[k];

                g_0_x_xxzzzz_xxzz[k] = -g_0_x_xxzzz_xxzz[k] * ab_z + g_0_x_xxzzz_xxzzz[k];

                g_0_x_xxzzzz_xyyy[k] = -g_0_x_xxzzz_xyyy[k] * ab_z + g_0_x_xxzzz_xyyyz[k];

                g_0_x_xxzzzz_xyyz[k] = -g_0_x_xxzzz_xyyz[k] * ab_z + g_0_x_xxzzz_xyyzz[k];

                g_0_x_xxzzzz_xyzz[k] = -g_0_x_xxzzz_xyzz[k] * ab_z + g_0_x_xxzzz_xyzzz[k];

                g_0_x_xxzzzz_xzzz[k] = -g_0_x_xxzzz_xzzz[k] * ab_z + g_0_x_xxzzz_xzzzz[k];

                g_0_x_xxzzzz_yyyy[k] = -g_0_x_xxzzz_yyyy[k] * ab_z + g_0_x_xxzzz_yyyyz[k];

                g_0_x_xxzzzz_yyyz[k] = -g_0_x_xxzzz_yyyz[k] * ab_z + g_0_x_xxzzz_yyyzz[k];

                g_0_x_xxzzzz_yyzz[k] = -g_0_x_xxzzz_yyzz[k] * ab_z + g_0_x_xxzzz_yyzzz[k];

                g_0_x_xxzzzz_yzzz[k] = -g_0_x_xxzzz_yzzz[k] * ab_z + g_0_x_xxzzz_yzzzz[k];

                g_0_x_xxzzzz_zzzz[k] = -g_0_x_xxzzz_zzzz[k] * ab_z + g_0_x_xxzzz_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyy_xxxx = cbuffer.data(ig_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxxy = cbuffer.data(ig_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxxz = cbuffer.data(ig_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxyy = cbuffer.data(ig_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxyz = cbuffer.data(ig_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxzz = cbuffer.data(ig_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyyy = cbuffer.data(ig_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyyz = cbuffer.data(ig_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyzz = cbuffer.data(ig_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xzzz = cbuffer.data(ig_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyyy = cbuffer.data(ig_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyyz = cbuffer.data(ig_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyzz = cbuffer.data(ig_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yzzz = cbuffer.data(ig_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_xyyyyy_zzzz = cbuffer.data(ig_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyy_xxxx, g_0_x_xyyyy_xxxxy, g_0_x_xyyyy_xxxy, g_0_x_xyyyy_xxxyy, g_0_x_xyyyy_xxxyz, g_0_x_xyyyy_xxxz, g_0_x_xyyyy_xxyy, g_0_x_xyyyy_xxyyy, g_0_x_xyyyy_xxyyz, g_0_x_xyyyy_xxyz, g_0_x_xyyyy_xxyzz, g_0_x_xyyyy_xxzz, g_0_x_xyyyy_xyyy, g_0_x_xyyyy_xyyyy, g_0_x_xyyyy_xyyyz, g_0_x_xyyyy_xyyz, g_0_x_xyyyy_xyyzz, g_0_x_xyyyy_xyzz, g_0_x_xyyyy_xyzzz, g_0_x_xyyyy_xzzz, g_0_x_xyyyy_yyyy, g_0_x_xyyyy_yyyyy, g_0_x_xyyyy_yyyyz, g_0_x_xyyyy_yyyz, g_0_x_xyyyy_yyyzz, g_0_x_xyyyy_yyzz, g_0_x_xyyyy_yyzzz, g_0_x_xyyyy_yzzz, g_0_x_xyyyy_yzzzz, g_0_x_xyyyy_zzzz, g_0_x_xyyyyy_xxxx, g_0_x_xyyyyy_xxxy, g_0_x_xyyyyy_xxxz, g_0_x_xyyyyy_xxyy, g_0_x_xyyyyy_xxyz, g_0_x_xyyyyy_xxzz, g_0_x_xyyyyy_xyyy, g_0_x_xyyyyy_xyyz, g_0_x_xyyyyy_xyzz, g_0_x_xyyyyy_xzzz, g_0_x_xyyyyy_yyyy, g_0_x_xyyyyy_yyyz, g_0_x_xyyyyy_yyzz, g_0_x_xyyyyy_yzzz, g_0_x_xyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyy_xxxx[k] = -g_0_x_xyyyy_xxxx[k] * ab_y + g_0_x_xyyyy_xxxxy[k];

                g_0_x_xyyyyy_xxxy[k] = -g_0_x_xyyyy_xxxy[k] * ab_y + g_0_x_xyyyy_xxxyy[k];

                g_0_x_xyyyyy_xxxz[k] = -g_0_x_xyyyy_xxxz[k] * ab_y + g_0_x_xyyyy_xxxyz[k];

                g_0_x_xyyyyy_xxyy[k] = -g_0_x_xyyyy_xxyy[k] * ab_y + g_0_x_xyyyy_xxyyy[k];

                g_0_x_xyyyyy_xxyz[k] = -g_0_x_xyyyy_xxyz[k] * ab_y + g_0_x_xyyyy_xxyyz[k];

                g_0_x_xyyyyy_xxzz[k] = -g_0_x_xyyyy_xxzz[k] * ab_y + g_0_x_xyyyy_xxyzz[k];

                g_0_x_xyyyyy_xyyy[k] = -g_0_x_xyyyy_xyyy[k] * ab_y + g_0_x_xyyyy_xyyyy[k];

                g_0_x_xyyyyy_xyyz[k] = -g_0_x_xyyyy_xyyz[k] * ab_y + g_0_x_xyyyy_xyyyz[k];

                g_0_x_xyyyyy_xyzz[k] = -g_0_x_xyyyy_xyzz[k] * ab_y + g_0_x_xyyyy_xyyzz[k];

                g_0_x_xyyyyy_xzzz[k] = -g_0_x_xyyyy_xzzz[k] * ab_y + g_0_x_xyyyy_xyzzz[k];

                g_0_x_xyyyyy_yyyy[k] = -g_0_x_xyyyy_yyyy[k] * ab_y + g_0_x_xyyyy_yyyyy[k];

                g_0_x_xyyyyy_yyyz[k] = -g_0_x_xyyyy_yyyz[k] * ab_y + g_0_x_xyyyy_yyyyz[k];

                g_0_x_xyyyyy_yyzz[k] = -g_0_x_xyyyy_yyzz[k] * ab_y + g_0_x_xyyyy_yyyzz[k];

                g_0_x_xyyyyy_yzzz[k] = -g_0_x_xyyyy_yzzz[k] * ab_y + g_0_x_xyyyy_yyzzz[k];

                g_0_x_xyyyyy_zzzz[k] = -g_0_x_xyyyy_zzzz[k] * ab_y + g_0_x_xyyyy_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyz_xxxx = cbuffer.data(ig_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxxy = cbuffer.data(ig_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxxz = cbuffer.data(ig_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxyy = cbuffer.data(ig_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxyz = cbuffer.data(ig_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxzz = cbuffer.data(ig_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyyy = cbuffer.data(ig_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyyz = cbuffer.data(ig_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyzz = cbuffer.data(ig_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xzzz = cbuffer.data(ig_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyyy = cbuffer.data(ig_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyyz = cbuffer.data(ig_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyzz = cbuffer.data(ig_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yzzz = cbuffer.data(ig_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_xyyyyz_zzzz = cbuffer.data(ig_geom_01_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyz_xxxx, g_0_x_xyyyyz_xxxy, g_0_x_xyyyyz_xxxz, g_0_x_xyyyyz_xxyy, g_0_x_xyyyyz_xxyz, g_0_x_xyyyyz_xxzz, g_0_x_xyyyyz_xyyy, g_0_x_xyyyyz_xyyz, g_0_x_xyyyyz_xyzz, g_0_x_xyyyyz_xzzz, g_0_x_xyyyyz_yyyy, g_0_x_xyyyyz_yyyz, g_0_x_xyyyyz_yyzz, g_0_x_xyyyyz_yzzz, g_0_x_xyyyyz_zzzz, g_0_x_xyyyz_xxxx, g_0_x_xyyyz_xxxxy, g_0_x_xyyyz_xxxy, g_0_x_xyyyz_xxxyy, g_0_x_xyyyz_xxxyz, g_0_x_xyyyz_xxxz, g_0_x_xyyyz_xxyy, g_0_x_xyyyz_xxyyy, g_0_x_xyyyz_xxyyz, g_0_x_xyyyz_xxyz, g_0_x_xyyyz_xxyzz, g_0_x_xyyyz_xxzz, g_0_x_xyyyz_xyyy, g_0_x_xyyyz_xyyyy, g_0_x_xyyyz_xyyyz, g_0_x_xyyyz_xyyz, g_0_x_xyyyz_xyyzz, g_0_x_xyyyz_xyzz, g_0_x_xyyyz_xyzzz, g_0_x_xyyyz_xzzz, g_0_x_xyyyz_yyyy, g_0_x_xyyyz_yyyyy, g_0_x_xyyyz_yyyyz, g_0_x_xyyyz_yyyz, g_0_x_xyyyz_yyyzz, g_0_x_xyyyz_yyzz, g_0_x_xyyyz_yyzzz, g_0_x_xyyyz_yzzz, g_0_x_xyyyz_yzzzz, g_0_x_xyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyz_xxxx[k] = -g_0_x_xyyyz_xxxx[k] * ab_y + g_0_x_xyyyz_xxxxy[k];

                g_0_x_xyyyyz_xxxy[k] = -g_0_x_xyyyz_xxxy[k] * ab_y + g_0_x_xyyyz_xxxyy[k];

                g_0_x_xyyyyz_xxxz[k] = -g_0_x_xyyyz_xxxz[k] * ab_y + g_0_x_xyyyz_xxxyz[k];

                g_0_x_xyyyyz_xxyy[k] = -g_0_x_xyyyz_xxyy[k] * ab_y + g_0_x_xyyyz_xxyyy[k];

                g_0_x_xyyyyz_xxyz[k] = -g_0_x_xyyyz_xxyz[k] * ab_y + g_0_x_xyyyz_xxyyz[k];

                g_0_x_xyyyyz_xxzz[k] = -g_0_x_xyyyz_xxzz[k] * ab_y + g_0_x_xyyyz_xxyzz[k];

                g_0_x_xyyyyz_xyyy[k] = -g_0_x_xyyyz_xyyy[k] * ab_y + g_0_x_xyyyz_xyyyy[k];

                g_0_x_xyyyyz_xyyz[k] = -g_0_x_xyyyz_xyyz[k] * ab_y + g_0_x_xyyyz_xyyyz[k];

                g_0_x_xyyyyz_xyzz[k] = -g_0_x_xyyyz_xyzz[k] * ab_y + g_0_x_xyyyz_xyyzz[k];

                g_0_x_xyyyyz_xzzz[k] = -g_0_x_xyyyz_xzzz[k] * ab_y + g_0_x_xyyyz_xyzzz[k];

                g_0_x_xyyyyz_yyyy[k] = -g_0_x_xyyyz_yyyy[k] * ab_y + g_0_x_xyyyz_yyyyy[k];

                g_0_x_xyyyyz_yyyz[k] = -g_0_x_xyyyz_yyyz[k] * ab_y + g_0_x_xyyyz_yyyyz[k];

                g_0_x_xyyyyz_yyzz[k] = -g_0_x_xyyyz_yyzz[k] * ab_y + g_0_x_xyyyz_yyyzz[k];

                g_0_x_xyyyyz_yzzz[k] = -g_0_x_xyyyz_yzzz[k] * ab_y + g_0_x_xyyyz_yyzzz[k];

                g_0_x_xyyyyz_zzzz[k] = -g_0_x_xyyyz_zzzz[k] * ab_y + g_0_x_xyyyz_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyzz_xxxx = cbuffer.data(ig_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxxy = cbuffer.data(ig_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxxz = cbuffer.data(ig_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxyy = cbuffer.data(ig_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxyz = cbuffer.data(ig_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxzz = cbuffer.data(ig_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyyy = cbuffer.data(ig_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyyz = cbuffer.data(ig_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyzz = cbuffer.data(ig_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xzzz = cbuffer.data(ig_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyyy = cbuffer.data(ig_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyyz = cbuffer.data(ig_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyzz = cbuffer.data(ig_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yzzz = cbuffer.data(ig_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_xyyyzz_zzzz = cbuffer.data(ig_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyzz_xxxx, g_0_x_xyyyzz_xxxy, g_0_x_xyyyzz_xxxz, g_0_x_xyyyzz_xxyy, g_0_x_xyyyzz_xxyz, g_0_x_xyyyzz_xxzz, g_0_x_xyyyzz_xyyy, g_0_x_xyyyzz_xyyz, g_0_x_xyyyzz_xyzz, g_0_x_xyyyzz_xzzz, g_0_x_xyyyzz_yyyy, g_0_x_xyyyzz_yyyz, g_0_x_xyyyzz_yyzz, g_0_x_xyyyzz_yzzz, g_0_x_xyyyzz_zzzz, g_0_x_xyyzz_xxxx, g_0_x_xyyzz_xxxxy, g_0_x_xyyzz_xxxy, g_0_x_xyyzz_xxxyy, g_0_x_xyyzz_xxxyz, g_0_x_xyyzz_xxxz, g_0_x_xyyzz_xxyy, g_0_x_xyyzz_xxyyy, g_0_x_xyyzz_xxyyz, g_0_x_xyyzz_xxyz, g_0_x_xyyzz_xxyzz, g_0_x_xyyzz_xxzz, g_0_x_xyyzz_xyyy, g_0_x_xyyzz_xyyyy, g_0_x_xyyzz_xyyyz, g_0_x_xyyzz_xyyz, g_0_x_xyyzz_xyyzz, g_0_x_xyyzz_xyzz, g_0_x_xyyzz_xyzzz, g_0_x_xyyzz_xzzz, g_0_x_xyyzz_yyyy, g_0_x_xyyzz_yyyyy, g_0_x_xyyzz_yyyyz, g_0_x_xyyzz_yyyz, g_0_x_xyyzz_yyyzz, g_0_x_xyyzz_yyzz, g_0_x_xyyzz_yyzzz, g_0_x_xyyzz_yzzz, g_0_x_xyyzz_yzzzz, g_0_x_xyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyzz_xxxx[k] = -g_0_x_xyyzz_xxxx[k] * ab_y + g_0_x_xyyzz_xxxxy[k];

                g_0_x_xyyyzz_xxxy[k] = -g_0_x_xyyzz_xxxy[k] * ab_y + g_0_x_xyyzz_xxxyy[k];

                g_0_x_xyyyzz_xxxz[k] = -g_0_x_xyyzz_xxxz[k] * ab_y + g_0_x_xyyzz_xxxyz[k];

                g_0_x_xyyyzz_xxyy[k] = -g_0_x_xyyzz_xxyy[k] * ab_y + g_0_x_xyyzz_xxyyy[k];

                g_0_x_xyyyzz_xxyz[k] = -g_0_x_xyyzz_xxyz[k] * ab_y + g_0_x_xyyzz_xxyyz[k];

                g_0_x_xyyyzz_xxzz[k] = -g_0_x_xyyzz_xxzz[k] * ab_y + g_0_x_xyyzz_xxyzz[k];

                g_0_x_xyyyzz_xyyy[k] = -g_0_x_xyyzz_xyyy[k] * ab_y + g_0_x_xyyzz_xyyyy[k];

                g_0_x_xyyyzz_xyyz[k] = -g_0_x_xyyzz_xyyz[k] * ab_y + g_0_x_xyyzz_xyyyz[k];

                g_0_x_xyyyzz_xyzz[k] = -g_0_x_xyyzz_xyzz[k] * ab_y + g_0_x_xyyzz_xyyzz[k];

                g_0_x_xyyyzz_xzzz[k] = -g_0_x_xyyzz_xzzz[k] * ab_y + g_0_x_xyyzz_xyzzz[k];

                g_0_x_xyyyzz_yyyy[k] = -g_0_x_xyyzz_yyyy[k] * ab_y + g_0_x_xyyzz_yyyyy[k];

                g_0_x_xyyyzz_yyyz[k] = -g_0_x_xyyzz_yyyz[k] * ab_y + g_0_x_xyyzz_yyyyz[k];

                g_0_x_xyyyzz_yyzz[k] = -g_0_x_xyyzz_yyzz[k] * ab_y + g_0_x_xyyzz_yyyzz[k];

                g_0_x_xyyyzz_yzzz[k] = -g_0_x_xyyzz_yzzz[k] * ab_y + g_0_x_xyyzz_yyzzz[k];

                g_0_x_xyyyzz_zzzz[k] = -g_0_x_xyyzz_zzzz[k] * ab_y + g_0_x_xyyzz_yzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzzz_xxxx = cbuffer.data(ig_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxxy = cbuffer.data(ig_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxxz = cbuffer.data(ig_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxyy = cbuffer.data(ig_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxyz = cbuffer.data(ig_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxzz = cbuffer.data(ig_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyyy = cbuffer.data(ig_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyyz = cbuffer.data(ig_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyzz = cbuffer.data(ig_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xzzz = cbuffer.data(ig_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyyy = cbuffer.data(ig_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyyz = cbuffer.data(ig_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyzz = cbuffer.data(ig_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yzzz = cbuffer.data(ig_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_xyyzzz_zzzz = cbuffer.data(ig_geom_01_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzzz_xxxx, g_0_x_xyyzzz_xxxy, g_0_x_xyyzzz_xxxz, g_0_x_xyyzzz_xxyy, g_0_x_xyyzzz_xxyz, g_0_x_xyyzzz_xxzz, g_0_x_xyyzzz_xyyy, g_0_x_xyyzzz_xyyz, g_0_x_xyyzzz_xyzz, g_0_x_xyyzzz_xzzz, g_0_x_xyyzzz_yyyy, g_0_x_xyyzzz_yyyz, g_0_x_xyyzzz_yyzz, g_0_x_xyyzzz_yzzz, g_0_x_xyyzzz_zzzz, g_0_x_xyzzz_xxxx, g_0_x_xyzzz_xxxxy, g_0_x_xyzzz_xxxy, g_0_x_xyzzz_xxxyy, g_0_x_xyzzz_xxxyz, g_0_x_xyzzz_xxxz, g_0_x_xyzzz_xxyy, g_0_x_xyzzz_xxyyy, g_0_x_xyzzz_xxyyz, g_0_x_xyzzz_xxyz, g_0_x_xyzzz_xxyzz, g_0_x_xyzzz_xxzz, g_0_x_xyzzz_xyyy, g_0_x_xyzzz_xyyyy, g_0_x_xyzzz_xyyyz, g_0_x_xyzzz_xyyz, g_0_x_xyzzz_xyyzz, g_0_x_xyzzz_xyzz, g_0_x_xyzzz_xyzzz, g_0_x_xyzzz_xzzz, g_0_x_xyzzz_yyyy, g_0_x_xyzzz_yyyyy, g_0_x_xyzzz_yyyyz, g_0_x_xyzzz_yyyz, g_0_x_xyzzz_yyyzz, g_0_x_xyzzz_yyzz, g_0_x_xyzzz_yyzzz, g_0_x_xyzzz_yzzz, g_0_x_xyzzz_yzzzz, g_0_x_xyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzzz_xxxx[k] = -g_0_x_xyzzz_xxxx[k] * ab_y + g_0_x_xyzzz_xxxxy[k];

                g_0_x_xyyzzz_xxxy[k] = -g_0_x_xyzzz_xxxy[k] * ab_y + g_0_x_xyzzz_xxxyy[k];

                g_0_x_xyyzzz_xxxz[k] = -g_0_x_xyzzz_xxxz[k] * ab_y + g_0_x_xyzzz_xxxyz[k];

                g_0_x_xyyzzz_xxyy[k] = -g_0_x_xyzzz_xxyy[k] * ab_y + g_0_x_xyzzz_xxyyy[k];

                g_0_x_xyyzzz_xxyz[k] = -g_0_x_xyzzz_xxyz[k] * ab_y + g_0_x_xyzzz_xxyyz[k];

                g_0_x_xyyzzz_xxzz[k] = -g_0_x_xyzzz_xxzz[k] * ab_y + g_0_x_xyzzz_xxyzz[k];

                g_0_x_xyyzzz_xyyy[k] = -g_0_x_xyzzz_xyyy[k] * ab_y + g_0_x_xyzzz_xyyyy[k];

                g_0_x_xyyzzz_xyyz[k] = -g_0_x_xyzzz_xyyz[k] * ab_y + g_0_x_xyzzz_xyyyz[k];

                g_0_x_xyyzzz_xyzz[k] = -g_0_x_xyzzz_xyzz[k] * ab_y + g_0_x_xyzzz_xyyzz[k];

                g_0_x_xyyzzz_xzzz[k] = -g_0_x_xyzzz_xzzz[k] * ab_y + g_0_x_xyzzz_xyzzz[k];

                g_0_x_xyyzzz_yyyy[k] = -g_0_x_xyzzz_yyyy[k] * ab_y + g_0_x_xyzzz_yyyyy[k];

                g_0_x_xyyzzz_yyyz[k] = -g_0_x_xyzzz_yyyz[k] * ab_y + g_0_x_xyzzz_yyyyz[k];

                g_0_x_xyyzzz_yyzz[k] = -g_0_x_xyzzz_yyzz[k] * ab_y + g_0_x_xyzzz_yyyzz[k];

                g_0_x_xyyzzz_yzzz[k] = -g_0_x_xyzzz_yzzz[k] * ab_y + g_0_x_xyzzz_yyzzz[k];

                g_0_x_xyyzzz_zzzz[k] = -g_0_x_xyzzz_zzzz[k] * ab_y + g_0_x_xyzzz_yzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzzz_xxxx = cbuffer.data(ig_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxxy = cbuffer.data(ig_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxxz = cbuffer.data(ig_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxyy = cbuffer.data(ig_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxyz = cbuffer.data(ig_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxzz = cbuffer.data(ig_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyyy = cbuffer.data(ig_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyyz = cbuffer.data(ig_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyzz = cbuffer.data(ig_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xzzz = cbuffer.data(ig_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyyy = cbuffer.data(ig_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyyz = cbuffer.data(ig_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyzz = cbuffer.data(ig_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yzzz = cbuffer.data(ig_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_xyzzzz_zzzz = cbuffer.data(ig_geom_01_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzzz_xxxx, g_0_x_xyzzzz_xxxy, g_0_x_xyzzzz_xxxz, g_0_x_xyzzzz_xxyy, g_0_x_xyzzzz_xxyz, g_0_x_xyzzzz_xxzz, g_0_x_xyzzzz_xyyy, g_0_x_xyzzzz_xyyz, g_0_x_xyzzzz_xyzz, g_0_x_xyzzzz_xzzz, g_0_x_xyzzzz_yyyy, g_0_x_xyzzzz_yyyz, g_0_x_xyzzzz_yyzz, g_0_x_xyzzzz_yzzz, g_0_x_xyzzzz_zzzz, g_0_x_xzzzz_xxxx, g_0_x_xzzzz_xxxxy, g_0_x_xzzzz_xxxy, g_0_x_xzzzz_xxxyy, g_0_x_xzzzz_xxxyz, g_0_x_xzzzz_xxxz, g_0_x_xzzzz_xxyy, g_0_x_xzzzz_xxyyy, g_0_x_xzzzz_xxyyz, g_0_x_xzzzz_xxyz, g_0_x_xzzzz_xxyzz, g_0_x_xzzzz_xxzz, g_0_x_xzzzz_xyyy, g_0_x_xzzzz_xyyyy, g_0_x_xzzzz_xyyyz, g_0_x_xzzzz_xyyz, g_0_x_xzzzz_xyyzz, g_0_x_xzzzz_xyzz, g_0_x_xzzzz_xyzzz, g_0_x_xzzzz_xzzz, g_0_x_xzzzz_yyyy, g_0_x_xzzzz_yyyyy, g_0_x_xzzzz_yyyyz, g_0_x_xzzzz_yyyz, g_0_x_xzzzz_yyyzz, g_0_x_xzzzz_yyzz, g_0_x_xzzzz_yyzzz, g_0_x_xzzzz_yzzz, g_0_x_xzzzz_yzzzz, g_0_x_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzzz_xxxx[k] = -g_0_x_xzzzz_xxxx[k] * ab_y + g_0_x_xzzzz_xxxxy[k];

                g_0_x_xyzzzz_xxxy[k] = -g_0_x_xzzzz_xxxy[k] * ab_y + g_0_x_xzzzz_xxxyy[k];

                g_0_x_xyzzzz_xxxz[k] = -g_0_x_xzzzz_xxxz[k] * ab_y + g_0_x_xzzzz_xxxyz[k];

                g_0_x_xyzzzz_xxyy[k] = -g_0_x_xzzzz_xxyy[k] * ab_y + g_0_x_xzzzz_xxyyy[k];

                g_0_x_xyzzzz_xxyz[k] = -g_0_x_xzzzz_xxyz[k] * ab_y + g_0_x_xzzzz_xxyyz[k];

                g_0_x_xyzzzz_xxzz[k] = -g_0_x_xzzzz_xxzz[k] * ab_y + g_0_x_xzzzz_xxyzz[k];

                g_0_x_xyzzzz_xyyy[k] = -g_0_x_xzzzz_xyyy[k] * ab_y + g_0_x_xzzzz_xyyyy[k];

                g_0_x_xyzzzz_xyyz[k] = -g_0_x_xzzzz_xyyz[k] * ab_y + g_0_x_xzzzz_xyyyz[k];

                g_0_x_xyzzzz_xyzz[k] = -g_0_x_xzzzz_xyzz[k] * ab_y + g_0_x_xzzzz_xyyzz[k];

                g_0_x_xyzzzz_xzzz[k] = -g_0_x_xzzzz_xzzz[k] * ab_y + g_0_x_xzzzz_xyzzz[k];

                g_0_x_xyzzzz_yyyy[k] = -g_0_x_xzzzz_yyyy[k] * ab_y + g_0_x_xzzzz_yyyyy[k];

                g_0_x_xyzzzz_yyyz[k] = -g_0_x_xzzzz_yyyz[k] * ab_y + g_0_x_xzzzz_yyyyz[k];

                g_0_x_xyzzzz_yyzz[k] = -g_0_x_xzzzz_yyzz[k] * ab_y + g_0_x_xzzzz_yyyzz[k];

                g_0_x_xyzzzz_yzzz[k] = -g_0_x_xzzzz_yzzz[k] * ab_y + g_0_x_xzzzz_yyzzz[k];

                g_0_x_xyzzzz_zzzz[k] = -g_0_x_xzzzz_zzzz[k] * ab_y + g_0_x_xzzzz_yzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzzz_xxxx = cbuffer.data(ig_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxxy = cbuffer.data(ig_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxxz = cbuffer.data(ig_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxyy = cbuffer.data(ig_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxyz = cbuffer.data(ig_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxzz = cbuffer.data(ig_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyyy = cbuffer.data(ig_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyyz = cbuffer.data(ig_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyzz = cbuffer.data(ig_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xzzz = cbuffer.data(ig_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyyy = cbuffer.data(ig_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyyz = cbuffer.data(ig_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyzz = cbuffer.data(ig_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yzzz = cbuffer.data(ig_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_xzzzzz_zzzz = cbuffer.data(ig_geom_01_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzzz_xxxx, g_0_x_xzzzz_xxxxz, g_0_x_xzzzz_xxxy, g_0_x_xzzzz_xxxyz, g_0_x_xzzzz_xxxz, g_0_x_xzzzz_xxxzz, g_0_x_xzzzz_xxyy, g_0_x_xzzzz_xxyyz, g_0_x_xzzzz_xxyz, g_0_x_xzzzz_xxyzz, g_0_x_xzzzz_xxzz, g_0_x_xzzzz_xxzzz, g_0_x_xzzzz_xyyy, g_0_x_xzzzz_xyyyz, g_0_x_xzzzz_xyyz, g_0_x_xzzzz_xyyzz, g_0_x_xzzzz_xyzz, g_0_x_xzzzz_xyzzz, g_0_x_xzzzz_xzzz, g_0_x_xzzzz_xzzzz, g_0_x_xzzzz_yyyy, g_0_x_xzzzz_yyyyz, g_0_x_xzzzz_yyyz, g_0_x_xzzzz_yyyzz, g_0_x_xzzzz_yyzz, g_0_x_xzzzz_yyzzz, g_0_x_xzzzz_yzzz, g_0_x_xzzzz_yzzzz, g_0_x_xzzzz_zzzz, g_0_x_xzzzz_zzzzz, g_0_x_xzzzzz_xxxx, g_0_x_xzzzzz_xxxy, g_0_x_xzzzzz_xxxz, g_0_x_xzzzzz_xxyy, g_0_x_xzzzzz_xxyz, g_0_x_xzzzzz_xxzz, g_0_x_xzzzzz_xyyy, g_0_x_xzzzzz_xyyz, g_0_x_xzzzzz_xyzz, g_0_x_xzzzzz_xzzz, g_0_x_xzzzzz_yyyy, g_0_x_xzzzzz_yyyz, g_0_x_xzzzzz_yyzz, g_0_x_xzzzzz_yzzz, g_0_x_xzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzzz_xxxx[k] = -g_0_x_xzzzz_xxxx[k] * ab_z + g_0_x_xzzzz_xxxxz[k];

                g_0_x_xzzzzz_xxxy[k] = -g_0_x_xzzzz_xxxy[k] * ab_z + g_0_x_xzzzz_xxxyz[k];

                g_0_x_xzzzzz_xxxz[k] = -g_0_x_xzzzz_xxxz[k] * ab_z + g_0_x_xzzzz_xxxzz[k];

                g_0_x_xzzzzz_xxyy[k] = -g_0_x_xzzzz_xxyy[k] * ab_z + g_0_x_xzzzz_xxyyz[k];

                g_0_x_xzzzzz_xxyz[k] = -g_0_x_xzzzz_xxyz[k] * ab_z + g_0_x_xzzzz_xxyzz[k];

                g_0_x_xzzzzz_xxzz[k] = -g_0_x_xzzzz_xxzz[k] * ab_z + g_0_x_xzzzz_xxzzz[k];

                g_0_x_xzzzzz_xyyy[k] = -g_0_x_xzzzz_xyyy[k] * ab_z + g_0_x_xzzzz_xyyyz[k];

                g_0_x_xzzzzz_xyyz[k] = -g_0_x_xzzzz_xyyz[k] * ab_z + g_0_x_xzzzz_xyyzz[k];

                g_0_x_xzzzzz_xyzz[k] = -g_0_x_xzzzz_xyzz[k] * ab_z + g_0_x_xzzzz_xyzzz[k];

                g_0_x_xzzzzz_xzzz[k] = -g_0_x_xzzzz_xzzz[k] * ab_z + g_0_x_xzzzz_xzzzz[k];

                g_0_x_xzzzzz_yyyy[k] = -g_0_x_xzzzz_yyyy[k] * ab_z + g_0_x_xzzzz_yyyyz[k];

                g_0_x_xzzzzz_yyyz[k] = -g_0_x_xzzzz_yyyz[k] * ab_z + g_0_x_xzzzz_yyyzz[k];

                g_0_x_xzzzzz_yyzz[k] = -g_0_x_xzzzz_yyzz[k] * ab_z + g_0_x_xzzzz_yyzzz[k];

                g_0_x_xzzzzz_yzzz[k] = -g_0_x_xzzzz_yzzz[k] * ab_z + g_0_x_xzzzz_yzzzz[k];

                g_0_x_xzzzzz_zzzz[k] = -g_0_x_xzzzz_zzzz[k] * ab_z + g_0_x_xzzzz_zzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyy_xxxx = cbuffer.data(ig_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxxy = cbuffer.data(ig_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxxz = cbuffer.data(ig_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxyy = cbuffer.data(ig_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxyz = cbuffer.data(ig_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxzz = cbuffer.data(ig_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyyy = cbuffer.data(ig_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyyz = cbuffer.data(ig_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyzz = cbuffer.data(ig_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xzzz = cbuffer.data(ig_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyyy = cbuffer.data(ig_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyyz = cbuffer.data(ig_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyzz = cbuffer.data(ig_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yzzz = cbuffer.data(ig_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_yyyyyy_zzzz = cbuffer.data(ig_geom_01_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyy_xxxx, g_0_x_yyyyy_xxxxy, g_0_x_yyyyy_xxxy, g_0_x_yyyyy_xxxyy, g_0_x_yyyyy_xxxyz, g_0_x_yyyyy_xxxz, g_0_x_yyyyy_xxyy, g_0_x_yyyyy_xxyyy, g_0_x_yyyyy_xxyyz, g_0_x_yyyyy_xxyz, g_0_x_yyyyy_xxyzz, g_0_x_yyyyy_xxzz, g_0_x_yyyyy_xyyy, g_0_x_yyyyy_xyyyy, g_0_x_yyyyy_xyyyz, g_0_x_yyyyy_xyyz, g_0_x_yyyyy_xyyzz, g_0_x_yyyyy_xyzz, g_0_x_yyyyy_xyzzz, g_0_x_yyyyy_xzzz, g_0_x_yyyyy_yyyy, g_0_x_yyyyy_yyyyy, g_0_x_yyyyy_yyyyz, g_0_x_yyyyy_yyyz, g_0_x_yyyyy_yyyzz, g_0_x_yyyyy_yyzz, g_0_x_yyyyy_yyzzz, g_0_x_yyyyy_yzzz, g_0_x_yyyyy_yzzzz, g_0_x_yyyyy_zzzz, g_0_x_yyyyyy_xxxx, g_0_x_yyyyyy_xxxy, g_0_x_yyyyyy_xxxz, g_0_x_yyyyyy_xxyy, g_0_x_yyyyyy_xxyz, g_0_x_yyyyyy_xxzz, g_0_x_yyyyyy_xyyy, g_0_x_yyyyyy_xyyz, g_0_x_yyyyyy_xyzz, g_0_x_yyyyyy_xzzz, g_0_x_yyyyyy_yyyy, g_0_x_yyyyyy_yyyz, g_0_x_yyyyyy_yyzz, g_0_x_yyyyyy_yzzz, g_0_x_yyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyy_xxxx[k] = -g_0_x_yyyyy_xxxx[k] * ab_y + g_0_x_yyyyy_xxxxy[k];

                g_0_x_yyyyyy_xxxy[k] = -g_0_x_yyyyy_xxxy[k] * ab_y + g_0_x_yyyyy_xxxyy[k];

                g_0_x_yyyyyy_xxxz[k] = -g_0_x_yyyyy_xxxz[k] * ab_y + g_0_x_yyyyy_xxxyz[k];

                g_0_x_yyyyyy_xxyy[k] = -g_0_x_yyyyy_xxyy[k] * ab_y + g_0_x_yyyyy_xxyyy[k];

                g_0_x_yyyyyy_xxyz[k] = -g_0_x_yyyyy_xxyz[k] * ab_y + g_0_x_yyyyy_xxyyz[k];

                g_0_x_yyyyyy_xxzz[k] = -g_0_x_yyyyy_xxzz[k] * ab_y + g_0_x_yyyyy_xxyzz[k];

                g_0_x_yyyyyy_xyyy[k] = -g_0_x_yyyyy_xyyy[k] * ab_y + g_0_x_yyyyy_xyyyy[k];

                g_0_x_yyyyyy_xyyz[k] = -g_0_x_yyyyy_xyyz[k] * ab_y + g_0_x_yyyyy_xyyyz[k];

                g_0_x_yyyyyy_xyzz[k] = -g_0_x_yyyyy_xyzz[k] * ab_y + g_0_x_yyyyy_xyyzz[k];

                g_0_x_yyyyyy_xzzz[k] = -g_0_x_yyyyy_xzzz[k] * ab_y + g_0_x_yyyyy_xyzzz[k];

                g_0_x_yyyyyy_yyyy[k] = -g_0_x_yyyyy_yyyy[k] * ab_y + g_0_x_yyyyy_yyyyy[k];

                g_0_x_yyyyyy_yyyz[k] = -g_0_x_yyyyy_yyyz[k] * ab_y + g_0_x_yyyyy_yyyyz[k];

                g_0_x_yyyyyy_yyzz[k] = -g_0_x_yyyyy_yyzz[k] * ab_y + g_0_x_yyyyy_yyyzz[k];

                g_0_x_yyyyyy_yzzz[k] = -g_0_x_yyyyy_yzzz[k] * ab_y + g_0_x_yyyyy_yyzzz[k];

                g_0_x_yyyyyy_zzzz[k] = -g_0_x_yyyyy_zzzz[k] * ab_y + g_0_x_yyyyy_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyz_xxxx = cbuffer.data(ig_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxxy = cbuffer.data(ig_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxxz = cbuffer.data(ig_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxyy = cbuffer.data(ig_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxyz = cbuffer.data(ig_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxzz = cbuffer.data(ig_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyyy = cbuffer.data(ig_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyyz = cbuffer.data(ig_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyzz = cbuffer.data(ig_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xzzz = cbuffer.data(ig_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyyy = cbuffer.data(ig_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyyz = cbuffer.data(ig_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyzz = cbuffer.data(ig_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yzzz = cbuffer.data(ig_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_yyyyyz_zzzz = cbuffer.data(ig_geom_01_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyz_xxxx, g_0_x_yyyyyz_xxxy, g_0_x_yyyyyz_xxxz, g_0_x_yyyyyz_xxyy, g_0_x_yyyyyz_xxyz, g_0_x_yyyyyz_xxzz, g_0_x_yyyyyz_xyyy, g_0_x_yyyyyz_xyyz, g_0_x_yyyyyz_xyzz, g_0_x_yyyyyz_xzzz, g_0_x_yyyyyz_yyyy, g_0_x_yyyyyz_yyyz, g_0_x_yyyyyz_yyzz, g_0_x_yyyyyz_yzzz, g_0_x_yyyyyz_zzzz, g_0_x_yyyyz_xxxx, g_0_x_yyyyz_xxxxy, g_0_x_yyyyz_xxxy, g_0_x_yyyyz_xxxyy, g_0_x_yyyyz_xxxyz, g_0_x_yyyyz_xxxz, g_0_x_yyyyz_xxyy, g_0_x_yyyyz_xxyyy, g_0_x_yyyyz_xxyyz, g_0_x_yyyyz_xxyz, g_0_x_yyyyz_xxyzz, g_0_x_yyyyz_xxzz, g_0_x_yyyyz_xyyy, g_0_x_yyyyz_xyyyy, g_0_x_yyyyz_xyyyz, g_0_x_yyyyz_xyyz, g_0_x_yyyyz_xyyzz, g_0_x_yyyyz_xyzz, g_0_x_yyyyz_xyzzz, g_0_x_yyyyz_xzzz, g_0_x_yyyyz_yyyy, g_0_x_yyyyz_yyyyy, g_0_x_yyyyz_yyyyz, g_0_x_yyyyz_yyyz, g_0_x_yyyyz_yyyzz, g_0_x_yyyyz_yyzz, g_0_x_yyyyz_yyzzz, g_0_x_yyyyz_yzzz, g_0_x_yyyyz_yzzzz, g_0_x_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyz_xxxx[k] = -g_0_x_yyyyz_xxxx[k] * ab_y + g_0_x_yyyyz_xxxxy[k];

                g_0_x_yyyyyz_xxxy[k] = -g_0_x_yyyyz_xxxy[k] * ab_y + g_0_x_yyyyz_xxxyy[k];

                g_0_x_yyyyyz_xxxz[k] = -g_0_x_yyyyz_xxxz[k] * ab_y + g_0_x_yyyyz_xxxyz[k];

                g_0_x_yyyyyz_xxyy[k] = -g_0_x_yyyyz_xxyy[k] * ab_y + g_0_x_yyyyz_xxyyy[k];

                g_0_x_yyyyyz_xxyz[k] = -g_0_x_yyyyz_xxyz[k] * ab_y + g_0_x_yyyyz_xxyyz[k];

                g_0_x_yyyyyz_xxzz[k] = -g_0_x_yyyyz_xxzz[k] * ab_y + g_0_x_yyyyz_xxyzz[k];

                g_0_x_yyyyyz_xyyy[k] = -g_0_x_yyyyz_xyyy[k] * ab_y + g_0_x_yyyyz_xyyyy[k];

                g_0_x_yyyyyz_xyyz[k] = -g_0_x_yyyyz_xyyz[k] * ab_y + g_0_x_yyyyz_xyyyz[k];

                g_0_x_yyyyyz_xyzz[k] = -g_0_x_yyyyz_xyzz[k] * ab_y + g_0_x_yyyyz_xyyzz[k];

                g_0_x_yyyyyz_xzzz[k] = -g_0_x_yyyyz_xzzz[k] * ab_y + g_0_x_yyyyz_xyzzz[k];

                g_0_x_yyyyyz_yyyy[k] = -g_0_x_yyyyz_yyyy[k] * ab_y + g_0_x_yyyyz_yyyyy[k];

                g_0_x_yyyyyz_yyyz[k] = -g_0_x_yyyyz_yyyz[k] * ab_y + g_0_x_yyyyz_yyyyz[k];

                g_0_x_yyyyyz_yyzz[k] = -g_0_x_yyyyz_yyzz[k] * ab_y + g_0_x_yyyyz_yyyzz[k];

                g_0_x_yyyyyz_yzzz[k] = -g_0_x_yyyyz_yzzz[k] * ab_y + g_0_x_yyyyz_yyzzz[k];

                g_0_x_yyyyyz_zzzz[k] = -g_0_x_yyyyz_zzzz[k] * ab_y + g_0_x_yyyyz_yzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyzz_xxxx = cbuffer.data(ig_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxxy = cbuffer.data(ig_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxxz = cbuffer.data(ig_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxyy = cbuffer.data(ig_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxyz = cbuffer.data(ig_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxzz = cbuffer.data(ig_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyyy = cbuffer.data(ig_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyyz = cbuffer.data(ig_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyzz = cbuffer.data(ig_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xzzz = cbuffer.data(ig_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyyy = cbuffer.data(ig_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyyz = cbuffer.data(ig_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyzz = cbuffer.data(ig_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yzzz = cbuffer.data(ig_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_yyyyzz_zzzz = cbuffer.data(ig_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyzz_xxxx, g_0_x_yyyyzz_xxxy, g_0_x_yyyyzz_xxxz, g_0_x_yyyyzz_xxyy, g_0_x_yyyyzz_xxyz, g_0_x_yyyyzz_xxzz, g_0_x_yyyyzz_xyyy, g_0_x_yyyyzz_xyyz, g_0_x_yyyyzz_xyzz, g_0_x_yyyyzz_xzzz, g_0_x_yyyyzz_yyyy, g_0_x_yyyyzz_yyyz, g_0_x_yyyyzz_yyzz, g_0_x_yyyyzz_yzzz, g_0_x_yyyyzz_zzzz, g_0_x_yyyzz_xxxx, g_0_x_yyyzz_xxxxy, g_0_x_yyyzz_xxxy, g_0_x_yyyzz_xxxyy, g_0_x_yyyzz_xxxyz, g_0_x_yyyzz_xxxz, g_0_x_yyyzz_xxyy, g_0_x_yyyzz_xxyyy, g_0_x_yyyzz_xxyyz, g_0_x_yyyzz_xxyz, g_0_x_yyyzz_xxyzz, g_0_x_yyyzz_xxzz, g_0_x_yyyzz_xyyy, g_0_x_yyyzz_xyyyy, g_0_x_yyyzz_xyyyz, g_0_x_yyyzz_xyyz, g_0_x_yyyzz_xyyzz, g_0_x_yyyzz_xyzz, g_0_x_yyyzz_xyzzz, g_0_x_yyyzz_xzzz, g_0_x_yyyzz_yyyy, g_0_x_yyyzz_yyyyy, g_0_x_yyyzz_yyyyz, g_0_x_yyyzz_yyyz, g_0_x_yyyzz_yyyzz, g_0_x_yyyzz_yyzz, g_0_x_yyyzz_yyzzz, g_0_x_yyyzz_yzzz, g_0_x_yyyzz_yzzzz, g_0_x_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyzz_xxxx[k] = -g_0_x_yyyzz_xxxx[k] * ab_y + g_0_x_yyyzz_xxxxy[k];

                g_0_x_yyyyzz_xxxy[k] = -g_0_x_yyyzz_xxxy[k] * ab_y + g_0_x_yyyzz_xxxyy[k];

                g_0_x_yyyyzz_xxxz[k] = -g_0_x_yyyzz_xxxz[k] * ab_y + g_0_x_yyyzz_xxxyz[k];

                g_0_x_yyyyzz_xxyy[k] = -g_0_x_yyyzz_xxyy[k] * ab_y + g_0_x_yyyzz_xxyyy[k];

                g_0_x_yyyyzz_xxyz[k] = -g_0_x_yyyzz_xxyz[k] * ab_y + g_0_x_yyyzz_xxyyz[k];

                g_0_x_yyyyzz_xxzz[k] = -g_0_x_yyyzz_xxzz[k] * ab_y + g_0_x_yyyzz_xxyzz[k];

                g_0_x_yyyyzz_xyyy[k] = -g_0_x_yyyzz_xyyy[k] * ab_y + g_0_x_yyyzz_xyyyy[k];

                g_0_x_yyyyzz_xyyz[k] = -g_0_x_yyyzz_xyyz[k] * ab_y + g_0_x_yyyzz_xyyyz[k];

                g_0_x_yyyyzz_xyzz[k] = -g_0_x_yyyzz_xyzz[k] * ab_y + g_0_x_yyyzz_xyyzz[k];

                g_0_x_yyyyzz_xzzz[k] = -g_0_x_yyyzz_xzzz[k] * ab_y + g_0_x_yyyzz_xyzzz[k];

                g_0_x_yyyyzz_yyyy[k] = -g_0_x_yyyzz_yyyy[k] * ab_y + g_0_x_yyyzz_yyyyy[k];

                g_0_x_yyyyzz_yyyz[k] = -g_0_x_yyyzz_yyyz[k] * ab_y + g_0_x_yyyzz_yyyyz[k];

                g_0_x_yyyyzz_yyzz[k] = -g_0_x_yyyzz_yyzz[k] * ab_y + g_0_x_yyyzz_yyyzz[k];

                g_0_x_yyyyzz_yzzz[k] = -g_0_x_yyyzz_yzzz[k] * ab_y + g_0_x_yyyzz_yyzzz[k];

                g_0_x_yyyyzz_zzzz[k] = -g_0_x_yyyzz_zzzz[k] * ab_y + g_0_x_yyyzz_yzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzzz_xxxx = cbuffer.data(ig_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxxy = cbuffer.data(ig_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxxz = cbuffer.data(ig_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxyy = cbuffer.data(ig_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxyz = cbuffer.data(ig_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxzz = cbuffer.data(ig_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyyy = cbuffer.data(ig_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyyz = cbuffer.data(ig_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyzz = cbuffer.data(ig_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xzzz = cbuffer.data(ig_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyyy = cbuffer.data(ig_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyyz = cbuffer.data(ig_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyzz = cbuffer.data(ig_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yzzz = cbuffer.data(ig_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_x_yyyzzz_zzzz = cbuffer.data(ig_geom_01_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzzz_xxxx, g_0_x_yyyzzz_xxxy, g_0_x_yyyzzz_xxxz, g_0_x_yyyzzz_xxyy, g_0_x_yyyzzz_xxyz, g_0_x_yyyzzz_xxzz, g_0_x_yyyzzz_xyyy, g_0_x_yyyzzz_xyyz, g_0_x_yyyzzz_xyzz, g_0_x_yyyzzz_xzzz, g_0_x_yyyzzz_yyyy, g_0_x_yyyzzz_yyyz, g_0_x_yyyzzz_yyzz, g_0_x_yyyzzz_yzzz, g_0_x_yyyzzz_zzzz, g_0_x_yyzzz_xxxx, g_0_x_yyzzz_xxxxy, g_0_x_yyzzz_xxxy, g_0_x_yyzzz_xxxyy, g_0_x_yyzzz_xxxyz, g_0_x_yyzzz_xxxz, g_0_x_yyzzz_xxyy, g_0_x_yyzzz_xxyyy, g_0_x_yyzzz_xxyyz, g_0_x_yyzzz_xxyz, g_0_x_yyzzz_xxyzz, g_0_x_yyzzz_xxzz, g_0_x_yyzzz_xyyy, g_0_x_yyzzz_xyyyy, g_0_x_yyzzz_xyyyz, g_0_x_yyzzz_xyyz, g_0_x_yyzzz_xyyzz, g_0_x_yyzzz_xyzz, g_0_x_yyzzz_xyzzz, g_0_x_yyzzz_xzzz, g_0_x_yyzzz_yyyy, g_0_x_yyzzz_yyyyy, g_0_x_yyzzz_yyyyz, g_0_x_yyzzz_yyyz, g_0_x_yyzzz_yyyzz, g_0_x_yyzzz_yyzz, g_0_x_yyzzz_yyzzz, g_0_x_yyzzz_yzzz, g_0_x_yyzzz_yzzzz, g_0_x_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzzz_xxxx[k] = -g_0_x_yyzzz_xxxx[k] * ab_y + g_0_x_yyzzz_xxxxy[k];

                g_0_x_yyyzzz_xxxy[k] = -g_0_x_yyzzz_xxxy[k] * ab_y + g_0_x_yyzzz_xxxyy[k];

                g_0_x_yyyzzz_xxxz[k] = -g_0_x_yyzzz_xxxz[k] * ab_y + g_0_x_yyzzz_xxxyz[k];

                g_0_x_yyyzzz_xxyy[k] = -g_0_x_yyzzz_xxyy[k] * ab_y + g_0_x_yyzzz_xxyyy[k];

                g_0_x_yyyzzz_xxyz[k] = -g_0_x_yyzzz_xxyz[k] * ab_y + g_0_x_yyzzz_xxyyz[k];

                g_0_x_yyyzzz_xxzz[k] = -g_0_x_yyzzz_xxzz[k] * ab_y + g_0_x_yyzzz_xxyzz[k];

                g_0_x_yyyzzz_xyyy[k] = -g_0_x_yyzzz_xyyy[k] * ab_y + g_0_x_yyzzz_xyyyy[k];

                g_0_x_yyyzzz_xyyz[k] = -g_0_x_yyzzz_xyyz[k] * ab_y + g_0_x_yyzzz_xyyyz[k];

                g_0_x_yyyzzz_xyzz[k] = -g_0_x_yyzzz_xyzz[k] * ab_y + g_0_x_yyzzz_xyyzz[k];

                g_0_x_yyyzzz_xzzz[k] = -g_0_x_yyzzz_xzzz[k] * ab_y + g_0_x_yyzzz_xyzzz[k];

                g_0_x_yyyzzz_yyyy[k] = -g_0_x_yyzzz_yyyy[k] * ab_y + g_0_x_yyzzz_yyyyy[k];

                g_0_x_yyyzzz_yyyz[k] = -g_0_x_yyzzz_yyyz[k] * ab_y + g_0_x_yyzzz_yyyyz[k];

                g_0_x_yyyzzz_yyzz[k] = -g_0_x_yyzzz_yyzz[k] * ab_y + g_0_x_yyzzz_yyyzz[k];

                g_0_x_yyyzzz_yzzz[k] = -g_0_x_yyzzz_yzzz[k] * ab_y + g_0_x_yyzzz_yyzzz[k];

                g_0_x_yyyzzz_zzzz[k] = -g_0_x_yyzzz_zzzz[k] * ab_y + g_0_x_yyzzz_yzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzzz_xxxx = cbuffer.data(ig_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxxy = cbuffer.data(ig_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxxz = cbuffer.data(ig_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxyy = cbuffer.data(ig_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxyz = cbuffer.data(ig_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxzz = cbuffer.data(ig_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyyy = cbuffer.data(ig_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyyz = cbuffer.data(ig_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyzz = cbuffer.data(ig_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xzzz = cbuffer.data(ig_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyyy = cbuffer.data(ig_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyyz = cbuffer.data(ig_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyzz = cbuffer.data(ig_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yzzz = cbuffer.data(ig_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_x_yyzzzz_zzzz = cbuffer.data(ig_geom_01_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzzz_xxxx, g_0_x_yyzzzz_xxxy, g_0_x_yyzzzz_xxxz, g_0_x_yyzzzz_xxyy, g_0_x_yyzzzz_xxyz, g_0_x_yyzzzz_xxzz, g_0_x_yyzzzz_xyyy, g_0_x_yyzzzz_xyyz, g_0_x_yyzzzz_xyzz, g_0_x_yyzzzz_xzzz, g_0_x_yyzzzz_yyyy, g_0_x_yyzzzz_yyyz, g_0_x_yyzzzz_yyzz, g_0_x_yyzzzz_yzzz, g_0_x_yyzzzz_zzzz, g_0_x_yzzzz_xxxx, g_0_x_yzzzz_xxxxy, g_0_x_yzzzz_xxxy, g_0_x_yzzzz_xxxyy, g_0_x_yzzzz_xxxyz, g_0_x_yzzzz_xxxz, g_0_x_yzzzz_xxyy, g_0_x_yzzzz_xxyyy, g_0_x_yzzzz_xxyyz, g_0_x_yzzzz_xxyz, g_0_x_yzzzz_xxyzz, g_0_x_yzzzz_xxzz, g_0_x_yzzzz_xyyy, g_0_x_yzzzz_xyyyy, g_0_x_yzzzz_xyyyz, g_0_x_yzzzz_xyyz, g_0_x_yzzzz_xyyzz, g_0_x_yzzzz_xyzz, g_0_x_yzzzz_xyzzz, g_0_x_yzzzz_xzzz, g_0_x_yzzzz_yyyy, g_0_x_yzzzz_yyyyy, g_0_x_yzzzz_yyyyz, g_0_x_yzzzz_yyyz, g_0_x_yzzzz_yyyzz, g_0_x_yzzzz_yyzz, g_0_x_yzzzz_yyzzz, g_0_x_yzzzz_yzzz, g_0_x_yzzzz_yzzzz, g_0_x_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzzz_xxxx[k] = -g_0_x_yzzzz_xxxx[k] * ab_y + g_0_x_yzzzz_xxxxy[k];

                g_0_x_yyzzzz_xxxy[k] = -g_0_x_yzzzz_xxxy[k] * ab_y + g_0_x_yzzzz_xxxyy[k];

                g_0_x_yyzzzz_xxxz[k] = -g_0_x_yzzzz_xxxz[k] * ab_y + g_0_x_yzzzz_xxxyz[k];

                g_0_x_yyzzzz_xxyy[k] = -g_0_x_yzzzz_xxyy[k] * ab_y + g_0_x_yzzzz_xxyyy[k];

                g_0_x_yyzzzz_xxyz[k] = -g_0_x_yzzzz_xxyz[k] * ab_y + g_0_x_yzzzz_xxyyz[k];

                g_0_x_yyzzzz_xxzz[k] = -g_0_x_yzzzz_xxzz[k] * ab_y + g_0_x_yzzzz_xxyzz[k];

                g_0_x_yyzzzz_xyyy[k] = -g_0_x_yzzzz_xyyy[k] * ab_y + g_0_x_yzzzz_xyyyy[k];

                g_0_x_yyzzzz_xyyz[k] = -g_0_x_yzzzz_xyyz[k] * ab_y + g_0_x_yzzzz_xyyyz[k];

                g_0_x_yyzzzz_xyzz[k] = -g_0_x_yzzzz_xyzz[k] * ab_y + g_0_x_yzzzz_xyyzz[k];

                g_0_x_yyzzzz_xzzz[k] = -g_0_x_yzzzz_xzzz[k] * ab_y + g_0_x_yzzzz_xyzzz[k];

                g_0_x_yyzzzz_yyyy[k] = -g_0_x_yzzzz_yyyy[k] * ab_y + g_0_x_yzzzz_yyyyy[k];

                g_0_x_yyzzzz_yyyz[k] = -g_0_x_yzzzz_yyyz[k] * ab_y + g_0_x_yzzzz_yyyyz[k];

                g_0_x_yyzzzz_yyzz[k] = -g_0_x_yzzzz_yyzz[k] * ab_y + g_0_x_yzzzz_yyyzz[k];

                g_0_x_yyzzzz_yzzz[k] = -g_0_x_yzzzz_yzzz[k] * ab_y + g_0_x_yzzzz_yyzzz[k];

                g_0_x_yyzzzz_zzzz[k] = -g_0_x_yzzzz_zzzz[k] * ab_y + g_0_x_yzzzz_yzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzzz_xxxx = cbuffer.data(ig_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxxy = cbuffer.data(ig_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxxz = cbuffer.data(ig_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxyy = cbuffer.data(ig_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxyz = cbuffer.data(ig_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxzz = cbuffer.data(ig_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyyy = cbuffer.data(ig_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyyz = cbuffer.data(ig_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyzz = cbuffer.data(ig_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xzzz = cbuffer.data(ig_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyyy = cbuffer.data(ig_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyyz = cbuffer.data(ig_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyzz = cbuffer.data(ig_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yzzz = cbuffer.data(ig_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_x_yzzzzz_zzzz = cbuffer.data(ig_geom_01_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzzz_xxxx, g_0_x_yzzzzz_xxxy, g_0_x_yzzzzz_xxxz, g_0_x_yzzzzz_xxyy, g_0_x_yzzzzz_xxyz, g_0_x_yzzzzz_xxzz, g_0_x_yzzzzz_xyyy, g_0_x_yzzzzz_xyyz, g_0_x_yzzzzz_xyzz, g_0_x_yzzzzz_xzzz, g_0_x_yzzzzz_yyyy, g_0_x_yzzzzz_yyyz, g_0_x_yzzzzz_yyzz, g_0_x_yzzzzz_yzzz, g_0_x_yzzzzz_zzzz, g_0_x_zzzzz_xxxx, g_0_x_zzzzz_xxxxy, g_0_x_zzzzz_xxxy, g_0_x_zzzzz_xxxyy, g_0_x_zzzzz_xxxyz, g_0_x_zzzzz_xxxz, g_0_x_zzzzz_xxyy, g_0_x_zzzzz_xxyyy, g_0_x_zzzzz_xxyyz, g_0_x_zzzzz_xxyz, g_0_x_zzzzz_xxyzz, g_0_x_zzzzz_xxzz, g_0_x_zzzzz_xyyy, g_0_x_zzzzz_xyyyy, g_0_x_zzzzz_xyyyz, g_0_x_zzzzz_xyyz, g_0_x_zzzzz_xyyzz, g_0_x_zzzzz_xyzz, g_0_x_zzzzz_xyzzz, g_0_x_zzzzz_xzzz, g_0_x_zzzzz_yyyy, g_0_x_zzzzz_yyyyy, g_0_x_zzzzz_yyyyz, g_0_x_zzzzz_yyyz, g_0_x_zzzzz_yyyzz, g_0_x_zzzzz_yyzz, g_0_x_zzzzz_yyzzz, g_0_x_zzzzz_yzzz, g_0_x_zzzzz_yzzzz, g_0_x_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzzz_xxxx[k] = -g_0_x_zzzzz_xxxx[k] * ab_y + g_0_x_zzzzz_xxxxy[k];

                g_0_x_yzzzzz_xxxy[k] = -g_0_x_zzzzz_xxxy[k] * ab_y + g_0_x_zzzzz_xxxyy[k];

                g_0_x_yzzzzz_xxxz[k] = -g_0_x_zzzzz_xxxz[k] * ab_y + g_0_x_zzzzz_xxxyz[k];

                g_0_x_yzzzzz_xxyy[k] = -g_0_x_zzzzz_xxyy[k] * ab_y + g_0_x_zzzzz_xxyyy[k];

                g_0_x_yzzzzz_xxyz[k] = -g_0_x_zzzzz_xxyz[k] * ab_y + g_0_x_zzzzz_xxyyz[k];

                g_0_x_yzzzzz_xxzz[k] = -g_0_x_zzzzz_xxzz[k] * ab_y + g_0_x_zzzzz_xxyzz[k];

                g_0_x_yzzzzz_xyyy[k] = -g_0_x_zzzzz_xyyy[k] * ab_y + g_0_x_zzzzz_xyyyy[k];

                g_0_x_yzzzzz_xyyz[k] = -g_0_x_zzzzz_xyyz[k] * ab_y + g_0_x_zzzzz_xyyyz[k];

                g_0_x_yzzzzz_xyzz[k] = -g_0_x_zzzzz_xyzz[k] * ab_y + g_0_x_zzzzz_xyyzz[k];

                g_0_x_yzzzzz_xzzz[k] = -g_0_x_zzzzz_xzzz[k] * ab_y + g_0_x_zzzzz_xyzzz[k];

                g_0_x_yzzzzz_yyyy[k] = -g_0_x_zzzzz_yyyy[k] * ab_y + g_0_x_zzzzz_yyyyy[k];

                g_0_x_yzzzzz_yyyz[k] = -g_0_x_zzzzz_yyyz[k] * ab_y + g_0_x_zzzzz_yyyyz[k];

                g_0_x_yzzzzz_yyzz[k] = -g_0_x_zzzzz_yyzz[k] * ab_y + g_0_x_zzzzz_yyyzz[k];

                g_0_x_yzzzzz_yzzz[k] = -g_0_x_zzzzz_yzzz[k] * ab_y + g_0_x_zzzzz_yyzzz[k];

                g_0_x_yzzzzz_zzzz[k] = -g_0_x_zzzzz_zzzz[k] * ab_y + g_0_x_zzzzz_yzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzzz_xxxx = cbuffer.data(ig_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxxy = cbuffer.data(ig_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxxz = cbuffer.data(ig_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxyy = cbuffer.data(ig_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxyz = cbuffer.data(ig_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxzz = cbuffer.data(ig_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyyy = cbuffer.data(ig_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyyz = cbuffer.data(ig_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyzz = cbuffer.data(ig_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xzzz = cbuffer.data(ig_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyyy = cbuffer.data(ig_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyyz = cbuffer.data(ig_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyzz = cbuffer.data(ig_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yzzz = cbuffer.data(ig_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_x_zzzzzz_zzzz = cbuffer.data(ig_geom_01_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzzz_xxxx, g_0_x_zzzzz_xxxxz, g_0_x_zzzzz_xxxy, g_0_x_zzzzz_xxxyz, g_0_x_zzzzz_xxxz, g_0_x_zzzzz_xxxzz, g_0_x_zzzzz_xxyy, g_0_x_zzzzz_xxyyz, g_0_x_zzzzz_xxyz, g_0_x_zzzzz_xxyzz, g_0_x_zzzzz_xxzz, g_0_x_zzzzz_xxzzz, g_0_x_zzzzz_xyyy, g_0_x_zzzzz_xyyyz, g_0_x_zzzzz_xyyz, g_0_x_zzzzz_xyyzz, g_0_x_zzzzz_xyzz, g_0_x_zzzzz_xyzzz, g_0_x_zzzzz_xzzz, g_0_x_zzzzz_xzzzz, g_0_x_zzzzz_yyyy, g_0_x_zzzzz_yyyyz, g_0_x_zzzzz_yyyz, g_0_x_zzzzz_yyyzz, g_0_x_zzzzz_yyzz, g_0_x_zzzzz_yyzzz, g_0_x_zzzzz_yzzz, g_0_x_zzzzz_yzzzz, g_0_x_zzzzz_zzzz, g_0_x_zzzzz_zzzzz, g_0_x_zzzzzz_xxxx, g_0_x_zzzzzz_xxxy, g_0_x_zzzzzz_xxxz, g_0_x_zzzzzz_xxyy, g_0_x_zzzzzz_xxyz, g_0_x_zzzzzz_xxzz, g_0_x_zzzzzz_xyyy, g_0_x_zzzzzz_xyyz, g_0_x_zzzzzz_xyzz, g_0_x_zzzzzz_xzzz, g_0_x_zzzzzz_yyyy, g_0_x_zzzzzz_yyyz, g_0_x_zzzzzz_yyzz, g_0_x_zzzzzz_yzzz, g_0_x_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzzz_xxxx[k] = -g_0_x_zzzzz_xxxx[k] * ab_z + g_0_x_zzzzz_xxxxz[k];

                g_0_x_zzzzzz_xxxy[k] = -g_0_x_zzzzz_xxxy[k] * ab_z + g_0_x_zzzzz_xxxyz[k];

                g_0_x_zzzzzz_xxxz[k] = -g_0_x_zzzzz_xxxz[k] * ab_z + g_0_x_zzzzz_xxxzz[k];

                g_0_x_zzzzzz_xxyy[k] = -g_0_x_zzzzz_xxyy[k] * ab_z + g_0_x_zzzzz_xxyyz[k];

                g_0_x_zzzzzz_xxyz[k] = -g_0_x_zzzzz_xxyz[k] * ab_z + g_0_x_zzzzz_xxyzz[k];

                g_0_x_zzzzzz_xxzz[k] = -g_0_x_zzzzz_xxzz[k] * ab_z + g_0_x_zzzzz_xxzzz[k];

                g_0_x_zzzzzz_xyyy[k] = -g_0_x_zzzzz_xyyy[k] * ab_z + g_0_x_zzzzz_xyyyz[k];

                g_0_x_zzzzzz_xyyz[k] = -g_0_x_zzzzz_xyyz[k] * ab_z + g_0_x_zzzzz_xyyzz[k];

                g_0_x_zzzzzz_xyzz[k] = -g_0_x_zzzzz_xyzz[k] * ab_z + g_0_x_zzzzz_xyzzz[k];

                g_0_x_zzzzzz_xzzz[k] = -g_0_x_zzzzz_xzzz[k] * ab_z + g_0_x_zzzzz_xzzzz[k];

                g_0_x_zzzzzz_yyyy[k] = -g_0_x_zzzzz_yyyy[k] * ab_z + g_0_x_zzzzz_yyyyz[k];

                g_0_x_zzzzzz_yyyz[k] = -g_0_x_zzzzz_yyyz[k] * ab_z + g_0_x_zzzzz_yyyzz[k];

                g_0_x_zzzzzz_yyzz[k] = -g_0_x_zzzzz_yyzz[k] * ab_z + g_0_x_zzzzz_yyzzz[k];

                g_0_x_zzzzzz_yzzz[k] = -g_0_x_zzzzz_yzzz[k] * ab_z + g_0_x_zzzzz_yzzzz[k];

                g_0_x_zzzzzz_zzzz[k] = -g_0_x_zzzzz_zzzz[k] * ab_z + g_0_x_zzzzz_zzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxx_xxxx = cbuffer.data(ig_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxxy = cbuffer.data(ig_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxxz = cbuffer.data(ig_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxyy = cbuffer.data(ig_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxyz = cbuffer.data(ig_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxzz = cbuffer.data(ig_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyyy = cbuffer.data(ig_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyyz = cbuffer.data(ig_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyzz = cbuffer.data(ig_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xzzz = cbuffer.data(ig_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyyy = cbuffer.data(ig_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyyz = cbuffer.data(ig_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyzz = cbuffer.data(ig_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yzzz = cbuffer.data(ig_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_y_xxxxxx_zzzz = cbuffer.data(ig_geom_01_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxx_xxxx, g_0_y_xxxxx_xxxxx, g_0_y_xxxxx_xxxxy, g_0_y_xxxxx_xxxxz, g_0_y_xxxxx_xxxy, g_0_y_xxxxx_xxxyy, g_0_y_xxxxx_xxxyz, g_0_y_xxxxx_xxxz, g_0_y_xxxxx_xxxzz, g_0_y_xxxxx_xxyy, g_0_y_xxxxx_xxyyy, g_0_y_xxxxx_xxyyz, g_0_y_xxxxx_xxyz, g_0_y_xxxxx_xxyzz, g_0_y_xxxxx_xxzz, g_0_y_xxxxx_xxzzz, g_0_y_xxxxx_xyyy, g_0_y_xxxxx_xyyyy, g_0_y_xxxxx_xyyyz, g_0_y_xxxxx_xyyz, g_0_y_xxxxx_xyyzz, g_0_y_xxxxx_xyzz, g_0_y_xxxxx_xyzzz, g_0_y_xxxxx_xzzz, g_0_y_xxxxx_xzzzz, g_0_y_xxxxx_yyyy, g_0_y_xxxxx_yyyz, g_0_y_xxxxx_yyzz, g_0_y_xxxxx_yzzz, g_0_y_xxxxx_zzzz, g_0_y_xxxxxx_xxxx, g_0_y_xxxxxx_xxxy, g_0_y_xxxxxx_xxxz, g_0_y_xxxxxx_xxyy, g_0_y_xxxxxx_xxyz, g_0_y_xxxxxx_xxzz, g_0_y_xxxxxx_xyyy, g_0_y_xxxxxx_xyyz, g_0_y_xxxxxx_xyzz, g_0_y_xxxxxx_xzzz, g_0_y_xxxxxx_yyyy, g_0_y_xxxxxx_yyyz, g_0_y_xxxxxx_yyzz, g_0_y_xxxxxx_yzzz, g_0_y_xxxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxx_xxxx[k] = -g_0_y_xxxxx_xxxx[k] * ab_x + g_0_y_xxxxx_xxxxx[k];

                g_0_y_xxxxxx_xxxy[k] = -g_0_y_xxxxx_xxxy[k] * ab_x + g_0_y_xxxxx_xxxxy[k];

                g_0_y_xxxxxx_xxxz[k] = -g_0_y_xxxxx_xxxz[k] * ab_x + g_0_y_xxxxx_xxxxz[k];

                g_0_y_xxxxxx_xxyy[k] = -g_0_y_xxxxx_xxyy[k] * ab_x + g_0_y_xxxxx_xxxyy[k];

                g_0_y_xxxxxx_xxyz[k] = -g_0_y_xxxxx_xxyz[k] * ab_x + g_0_y_xxxxx_xxxyz[k];

                g_0_y_xxxxxx_xxzz[k] = -g_0_y_xxxxx_xxzz[k] * ab_x + g_0_y_xxxxx_xxxzz[k];

                g_0_y_xxxxxx_xyyy[k] = -g_0_y_xxxxx_xyyy[k] * ab_x + g_0_y_xxxxx_xxyyy[k];

                g_0_y_xxxxxx_xyyz[k] = -g_0_y_xxxxx_xyyz[k] * ab_x + g_0_y_xxxxx_xxyyz[k];

                g_0_y_xxxxxx_xyzz[k] = -g_0_y_xxxxx_xyzz[k] * ab_x + g_0_y_xxxxx_xxyzz[k];

                g_0_y_xxxxxx_xzzz[k] = -g_0_y_xxxxx_xzzz[k] * ab_x + g_0_y_xxxxx_xxzzz[k];

                g_0_y_xxxxxx_yyyy[k] = -g_0_y_xxxxx_yyyy[k] * ab_x + g_0_y_xxxxx_xyyyy[k];

                g_0_y_xxxxxx_yyyz[k] = -g_0_y_xxxxx_yyyz[k] * ab_x + g_0_y_xxxxx_xyyyz[k];

                g_0_y_xxxxxx_yyzz[k] = -g_0_y_xxxxx_yyzz[k] * ab_x + g_0_y_xxxxx_xyyzz[k];

                g_0_y_xxxxxx_yzzz[k] = -g_0_y_xxxxx_yzzz[k] * ab_x + g_0_y_xxxxx_xyzzz[k];

                g_0_y_xxxxxx_zzzz[k] = -g_0_y_xxxxx_zzzz[k] * ab_x + g_0_y_xxxxx_xzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxy_xxxx = cbuffer.data(ig_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxxy = cbuffer.data(ig_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxxz = cbuffer.data(ig_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxyy = cbuffer.data(ig_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxyz = cbuffer.data(ig_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxzz = cbuffer.data(ig_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyyy = cbuffer.data(ig_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyyz = cbuffer.data(ig_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyzz = cbuffer.data(ig_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xzzz = cbuffer.data(ig_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyyy = cbuffer.data(ig_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyyz = cbuffer.data(ig_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyzz = cbuffer.data(ig_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yzzz = cbuffer.data(ig_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_y_xxxxxy_zzzz = cbuffer.data(ig_geom_01_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxy_xxxx, g_0_y_xxxxxy_xxxy, g_0_y_xxxxxy_xxxz, g_0_y_xxxxxy_xxyy, g_0_y_xxxxxy_xxyz, g_0_y_xxxxxy_xxzz, g_0_y_xxxxxy_xyyy, g_0_y_xxxxxy_xyyz, g_0_y_xxxxxy_xyzz, g_0_y_xxxxxy_xzzz, g_0_y_xxxxxy_yyyy, g_0_y_xxxxxy_yyyz, g_0_y_xxxxxy_yyzz, g_0_y_xxxxxy_yzzz, g_0_y_xxxxxy_zzzz, g_0_y_xxxxy_xxxx, g_0_y_xxxxy_xxxxx, g_0_y_xxxxy_xxxxy, g_0_y_xxxxy_xxxxz, g_0_y_xxxxy_xxxy, g_0_y_xxxxy_xxxyy, g_0_y_xxxxy_xxxyz, g_0_y_xxxxy_xxxz, g_0_y_xxxxy_xxxzz, g_0_y_xxxxy_xxyy, g_0_y_xxxxy_xxyyy, g_0_y_xxxxy_xxyyz, g_0_y_xxxxy_xxyz, g_0_y_xxxxy_xxyzz, g_0_y_xxxxy_xxzz, g_0_y_xxxxy_xxzzz, g_0_y_xxxxy_xyyy, g_0_y_xxxxy_xyyyy, g_0_y_xxxxy_xyyyz, g_0_y_xxxxy_xyyz, g_0_y_xxxxy_xyyzz, g_0_y_xxxxy_xyzz, g_0_y_xxxxy_xyzzz, g_0_y_xxxxy_xzzz, g_0_y_xxxxy_xzzzz, g_0_y_xxxxy_yyyy, g_0_y_xxxxy_yyyz, g_0_y_xxxxy_yyzz, g_0_y_xxxxy_yzzz, g_0_y_xxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxy_xxxx[k] = -g_0_y_xxxxy_xxxx[k] * ab_x + g_0_y_xxxxy_xxxxx[k];

                g_0_y_xxxxxy_xxxy[k] = -g_0_y_xxxxy_xxxy[k] * ab_x + g_0_y_xxxxy_xxxxy[k];

                g_0_y_xxxxxy_xxxz[k] = -g_0_y_xxxxy_xxxz[k] * ab_x + g_0_y_xxxxy_xxxxz[k];

                g_0_y_xxxxxy_xxyy[k] = -g_0_y_xxxxy_xxyy[k] * ab_x + g_0_y_xxxxy_xxxyy[k];

                g_0_y_xxxxxy_xxyz[k] = -g_0_y_xxxxy_xxyz[k] * ab_x + g_0_y_xxxxy_xxxyz[k];

                g_0_y_xxxxxy_xxzz[k] = -g_0_y_xxxxy_xxzz[k] * ab_x + g_0_y_xxxxy_xxxzz[k];

                g_0_y_xxxxxy_xyyy[k] = -g_0_y_xxxxy_xyyy[k] * ab_x + g_0_y_xxxxy_xxyyy[k];

                g_0_y_xxxxxy_xyyz[k] = -g_0_y_xxxxy_xyyz[k] * ab_x + g_0_y_xxxxy_xxyyz[k];

                g_0_y_xxxxxy_xyzz[k] = -g_0_y_xxxxy_xyzz[k] * ab_x + g_0_y_xxxxy_xxyzz[k];

                g_0_y_xxxxxy_xzzz[k] = -g_0_y_xxxxy_xzzz[k] * ab_x + g_0_y_xxxxy_xxzzz[k];

                g_0_y_xxxxxy_yyyy[k] = -g_0_y_xxxxy_yyyy[k] * ab_x + g_0_y_xxxxy_xyyyy[k];

                g_0_y_xxxxxy_yyyz[k] = -g_0_y_xxxxy_yyyz[k] * ab_x + g_0_y_xxxxy_xyyyz[k];

                g_0_y_xxxxxy_yyzz[k] = -g_0_y_xxxxy_yyzz[k] * ab_x + g_0_y_xxxxy_xyyzz[k];

                g_0_y_xxxxxy_yzzz[k] = -g_0_y_xxxxy_yzzz[k] * ab_x + g_0_y_xxxxy_xyzzz[k];

                g_0_y_xxxxxy_zzzz[k] = -g_0_y_xxxxy_zzzz[k] * ab_x + g_0_y_xxxxy_xzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxz_xxxx = cbuffer.data(ig_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxxy = cbuffer.data(ig_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxxz = cbuffer.data(ig_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxyy = cbuffer.data(ig_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxyz = cbuffer.data(ig_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxzz = cbuffer.data(ig_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyyy = cbuffer.data(ig_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyyz = cbuffer.data(ig_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyzz = cbuffer.data(ig_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xzzz = cbuffer.data(ig_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyyy = cbuffer.data(ig_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyyz = cbuffer.data(ig_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyzz = cbuffer.data(ig_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yzzz = cbuffer.data(ig_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_y_xxxxxz_zzzz = cbuffer.data(ig_geom_01_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxz_xxxx, g_0_y_xxxxxz_xxxy, g_0_y_xxxxxz_xxxz, g_0_y_xxxxxz_xxyy, g_0_y_xxxxxz_xxyz, g_0_y_xxxxxz_xxzz, g_0_y_xxxxxz_xyyy, g_0_y_xxxxxz_xyyz, g_0_y_xxxxxz_xyzz, g_0_y_xxxxxz_xzzz, g_0_y_xxxxxz_yyyy, g_0_y_xxxxxz_yyyz, g_0_y_xxxxxz_yyzz, g_0_y_xxxxxz_yzzz, g_0_y_xxxxxz_zzzz, g_0_y_xxxxz_xxxx, g_0_y_xxxxz_xxxxx, g_0_y_xxxxz_xxxxy, g_0_y_xxxxz_xxxxz, g_0_y_xxxxz_xxxy, g_0_y_xxxxz_xxxyy, g_0_y_xxxxz_xxxyz, g_0_y_xxxxz_xxxz, g_0_y_xxxxz_xxxzz, g_0_y_xxxxz_xxyy, g_0_y_xxxxz_xxyyy, g_0_y_xxxxz_xxyyz, g_0_y_xxxxz_xxyz, g_0_y_xxxxz_xxyzz, g_0_y_xxxxz_xxzz, g_0_y_xxxxz_xxzzz, g_0_y_xxxxz_xyyy, g_0_y_xxxxz_xyyyy, g_0_y_xxxxz_xyyyz, g_0_y_xxxxz_xyyz, g_0_y_xxxxz_xyyzz, g_0_y_xxxxz_xyzz, g_0_y_xxxxz_xyzzz, g_0_y_xxxxz_xzzz, g_0_y_xxxxz_xzzzz, g_0_y_xxxxz_yyyy, g_0_y_xxxxz_yyyz, g_0_y_xxxxz_yyzz, g_0_y_xxxxz_yzzz, g_0_y_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxz_xxxx[k] = -g_0_y_xxxxz_xxxx[k] * ab_x + g_0_y_xxxxz_xxxxx[k];

                g_0_y_xxxxxz_xxxy[k] = -g_0_y_xxxxz_xxxy[k] * ab_x + g_0_y_xxxxz_xxxxy[k];

                g_0_y_xxxxxz_xxxz[k] = -g_0_y_xxxxz_xxxz[k] * ab_x + g_0_y_xxxxz_xxxxz[k];

                g_0_y_xxxxxz_xxyy[k] = -g_0_y_xxxxz_xxyy[k] * ab_x + g_0_y_xxxxz_xxxyy[k];

                g_0_y_xxxxxz_xxyz[k] = -g_0_y_xxxxz_xxyz[k] * ab_x + g_0_y_xxxxz_xxxyz[k];

                g_0_y_xxxxxz_xxzz[k] = -g_0_y_xxxxz_xxzz[k] * ab_x + g_0_y_xxxxz_xxxzz[k];

                g_0_y_xxxxxz_xyyy[k] = -g_0_y_xxxxz_xyyy[k] * ab_x + g_0_y_xxxxz_xxyyy[k];

                g_0_y_xxxxxz_xyyz[k] = -g_0_y_xxxxz_xyyz[k] * ab_x + g_0_y_xxxxz_xxyyz[k];

                g_0_y_xxxxxz_xyzz[k] = -g_0_y_xxxxz_xyzz[k] * ab_x + g_0_y_xxxxz_xxyzz[k];

                g_0_y_xxxxxz_xzzz[k] = -g_0_y_xxxxz_xzzz[k] * ab_x + g_0_y_xxxxz_xxzzz[k];

                g_0_y_xxxxxz_yyyy[k] = -g_0_y_xxxxz_yyyy[k] * ab_x + g_0_y_xxxxz_xyyyy[k];

                g_0_y_xxxxxz_yyyz[k] = -g_0_y_xxxxz_yyyz[k] * ab_x + g_0_y_xxxxz_xyyyz[k];

                g_0_y_xxxxxz_yyzz[k] = -g_0_y_xxxxz_yyzz[k] * ab_x + g_0_y_xxxxz_xyyzz[k];

                g_0_y_xxxxxz_yzzz[k] = -g_0_y_xxxxz_yzzz[k] * ab_x + g_0_y_xxxxz_xyzzz[k];

                g_0_y_xxxxxz_zzzz[k] = -g_0_y_xxxxz_zzzz[k] * ab_x + g_0_y_xxxxz_xzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyy_xxxx = cbuffer.data(ig_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxxy = cbuffer.data(ig_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxxz = cbuffer.data(ig_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxyy = cbuffer.data(ig_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxyz = cbuffer.data(ig_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxzz = cbuffer.data(ig_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyyy = cbuffer.data(ig_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyyz = cbuffer.data(ig_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyzz = cbuffer.data(ig_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xzzz = cbuffer.data(ig_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyyy = cbuffer.data(ig_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyyz = cbuffer.data(ig_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyzz = cbuffer.data(ig_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yzzz = cbuffer.data(ig_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_y_xxxxyy_zzzz = cbuffer.data(ig_geom_01_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyy_xxxx, g_0_y_xxxxyy_xxxy, g_0_y_xxxxyy_xxxz, g_0_y_xxxxyy_xxyy, g_0_y_xxxxyy_xxyz, g_0_y_xxxxyy_xxzz, g_0_y_xxxxyy_xyyy, g_0_y_xxxxyy_xyyz, g_0_y_xxxxyy_xyzz, g_0_y_xxxxyy_xzzz, g_0_y_xxxxyy_yyyy, g_0_y_xxxxyy_yyyz, g_0_y_xxxxyy_yyzz, g_0_y_xxxxyy_yzzz, g_0_y_xxxxyy_zzzz, g_0_y_xxxyy_xxxx, g_0_y_xxxyy_xxxxx, g_0_y_xxxyy_xxxxy, g_0_y_xxxyy_xxxxz, g_0_y_xxxyy_xxxy, g_0_y_xxxyy_xxxyy, g_0_y_xxxyy_xxxyz, g_0_y_xxxyy_xxxz, g_0_y_xxxyy_xxxzz, g_0_y_xxxyy_xxyy, g_0_y_xxxyy_xxyyy, g_0_y_xxxyy_xxyyz, g_0_y_xxxyy_xxyz, g_0_y_xxxyy_xxyzz, g_0_y_xxxyy_xxzz, g_0_y_xxxyy_xxzzz, g_0_y_xxxyy_xyyy, g_0_y_xxxyy_xyyyy, g_0_y_xxxyy_xyyyz, g_0_y_xxxyy_xyyz, g_0_y_xxxyy_xyyzz, g_0_y_xxxyy_xyzz, g_0_y_xxxyy_xyzzz, g_0_y_xxxyy_xzzz, g_0_y_xxxyy_xzzzz, g_0_y_xxxyy_yyyy, g_0_y_xxxyy_yyyz, g_0_y_xxxyy_yyzz, g_0_y_xxxyy_yzzz, g_0_y_xxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyy_xxxx[k] = -g_0_y_xxxyy_xxxx[k] * ab_x + g_0_y_xxxyy_xxxxx[k];

                g_0_y_xxxxyy_xxxy[k] = -g_0_y_xxxyy_xxxy[k] * ab_x + g_0_y_xxxyy_xxxxy[k];

                g_0_y_xxxxyy_xxxz[k] = -g_0_y_xxxyy_xxxz[k] * ab_x + g_0_y_xxxyy_xxxxz[k];

                g_0_y_xxxxyy_xxyy[k] = -g_0_y_xxxyy_xxyy[k] * ab_x + g_0_y_xxxyy_xxxyy[k];

                g_0_y_xxxxyy_xxyz[k] = -g_0_y_xxxyy_xxyz[k] * ab_x + g_0_y_xxxyy_xxxyz[k];

                g_0_y_xxxxyy_xxzz[k] = -g_0_y_xxxyy_xxzz[k] * ab_x + g_0_y_xxxyy_xxxzz[k];

                g_0_y_xxxxyy_xyyy[k] = -g_0_y_xxxyy_xyyy[k] * ab_x + g_0_y_xxxyy_xxyyy[k];

                g_0_y_xxxxyy_xyyz[k] = -g_0_y_xxxyy_xyyz[k] * ab_x + g_0_y_xxxyy_xxyyz[k];

                g_0_y_xxxxyy_xyzz[k] = -g_0_y_xxxyy_xyzz[k] * ab_x + g_0_y_xxxyy_xxyzz[k];

                g_0_y_xxxxyy_xzzz[k] = -g_0_y_xxxyy_xzzz[k] * ab_x + g_0_y_xxxyy_xxzzz[k];

                g_0_y_xxxxyy_yyyy[k] = -g_0_y_xxxyy_yyyy[k] * ab_x + g_0_y_xxxyy_xyyyy[k];

                g_0_y_xxxxyy_yyyz[k] = -g_0_y_xxxyy_yyyz[k] * ab_x + g_0_y_xxxyy_xyyyz[k];

                g_0_y_xxxxyy_yyzz[k] = -g_0_y_xxxyy_yyzz[k] * ab_x + g_0_y_xxxyy_xyyzz[k];

                g_0_y_xxxxyy_yzzz[k] = -g_0_y_xxxyy_yzzz[k] * ab_x + g_0_y_xxxyy_xyzzz[k];

                g_0_y_xxxxyy_zzzz[k] = -g_0_y_xxxyy_zzzz[k] * ab_x + g_0_y_xxxyy_xzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyz_xxxx = cbuffer.data(ig_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxxy = cbuffer.data(ig_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxxz = cbuffer.data(ig_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxyy = cbuffer.data(ig_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxyz = cbuffer.data(ig_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxzz = cbuffer.data(ig_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyyy = cbuffer.data(ig_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyyz = cbuffer.data(ig_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyzz = cbuffer.data(ig_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xzzz = cbuffer.data(ig_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyyy = cbuffer.data(ig_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyyz = cbuffer.data(ig_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyzz = cbuffer.data(ig_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yzzz = cbuffer.data(ig_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_y_xxxxyz_zzzz = cbuffer.data(ig_geom_01_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyz_xxxx, g_0_y_xxxxyz_xxxy, g_0_y_xxxxyz_xxxz, g_0_y_xxxxyz_xxyy, g_0_y_xxxxyz_xxyz, g_0_y_xxxxyz_xxzz, g_0_y_xxxxyz_xyyy, g_0_y_xxxxyz_xyyz, g_0_y_xxxxyz_xyzz, g_0_y_xxxxyz_xzzz, g_0_y_xxxxyz_yyyy, g_0_y_xxxxyz_yyyz, g_0_y_xxxxyz_yyzz, g_0_y_xxxxyz_yzzz, g_0_y_xxxxyz_zzzz, g_0_y_xxxyz_xxxx, g_0_y_xxxyz_xxxxx, g_0_y_xxxyz_xxxxy, g_0_y_xxxyz_xxxxz, g_0_y_xxxyz_xxxy, g_0_y_xxxyz_xxxyy, g_0_y_xxxyz_xxxyz, g_0_y_xxxyz_xxxz, g_0_y_xxxyz_xxxzz, g_0_y_xxxyz_xxyy, g_0_y_xxxyz_xxyyy, g_0_y_xxxyz_xxyyz, g_0_y_xxxyz_xxyz, g_0_y_xxxyz_xxyzz, g_0_y_xxxyz_xxzz, g_0_y_xxxyz_xxzzz, g_0_y_xxxyz_xyyy, g_0_y_xxxyz_xyyyy, g_0_y_xxxyz_xyyyz, g_0_y_xxxyz_xyyz, g_0_y_xxxyz_xyyzz, g_0_y_xxxyz_xyzz, g_0_y_xxxyz_xyzzz, g_0_y_xxxyz_xzzz, g_0_y_xxxyz_xzzzz, g_0_y_xxxyz_yyyy, g_0_y_xxxyz_yyyz, g_0_y_xxxyz_yyzz, g_0_y_xxxyz_yzzz, g_0_y_xxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyz_xxxx[k] = -g_0_y_xxxyz_xxxx[k] * ab_x + g_0_y_xxxyz_xxxxx[k];

                g_0_y_xxxxyz_xxxy[k] = -g_0_y_xxxyz_xxxy[k] * ab_x + g_0_y_xxxyz_xxxxy[k];

                g_0_y_xxxxyz_xxxz[k] = -g_0_y_xxxyz_xxxz[k] * ab_x + g_0_y_xxxyz_xxxxz[k];

                g_0_y_xxxxyz_xxyy[k] = -g_0_y_xxxyz_xxyy[k] * ab_x + g_0_y_xxxyz_xxxyy[k];

                g_0_y_xxxxyz_xxyz[k] = -g_0_y_xxxyz_xxyz[k] * ab_x + g_0_y_xxxyz_xxxyz[k];

                g_0_y_xxxxyz_xxzz[k] = -g_0_y_xxxyz_xxzz[k] * ab_x + g_0_y_xxxyz_xxxzz[k];

                g_0_y_xxxxyz_xyyy[k] = -g_0_y_xxxyz_xyyy[k] * ab_x + g_0_y_xxxyz_xxyyy[k];

                g_0_y_xxxxyz_xyyz[k] = -g_0_y_xxxyz_xyyz[k] * ab_x + g_0_y_xxxyz_xxyyz[k];

                g_0_y_xxxxyz_xyzz[k] = -g_0_y_xxxyz_xyzz[k] * ab_x + g_0_y_xxxyz_xxyzz[k];

                g_0_y_xxxxyz_xzzz[k] = -g_0_y_xxxyz_xzzz[k] * ab_x + g_0_y_xxxyz_xxzzz[k];

                g_0_y_xxxxyz_yyyy[k] = -g_0_y_xxxyz_yyyy[k] * ab_x + g_0_y_xxxyz_xyyyy[k];

                g_0_y_xxxxyz_yyyz[k] = -g_0_y_xxxyz_yyyz[k] * ab_x + g_0_y_xxxyz_xyyyz[k];

                g_0_y_xxxxyz_yyzz[k] = -g_0_y_xxxyz_yyzz[k] * ab_x + g_0_y_xxxyz_xyyzz[k];

                g_0_y_xxxxyz_yzzz[k] = -g_0_y_xxxyz_yzzz[k] * ab_x + g_0_y_xxxyz_xyzzz[k];

                g_0_y_xxxxyz_zzzz[k] = -g_0_y_xxxyz_zzzz[k] * ab_x + g_0_y_xxxyz_xzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxzz_xxxx = cbuffer.data(ig_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxxy = cbuffer.data(ig_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxxz = cbuffer.data(ig_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxyy = cbuffer.data(ig_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxyz = cbuffer.data(ig_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxzz = cbuffer.data(ig_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyyy = cbuffer.data(ig_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyyz = cbuffer.data(ig_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyzz = cbuffer.data(ig_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xzzz = cbuffer.data(ig_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyyy = cbuffer.data(ig_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyyz = cbuffer.data(ig_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyzz = cbuffer.data(ig_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yzzz = cbuffer.data(ig_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_y_xxxxzz_zzzz = cbuffer.data(ig_geom_01_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxzz_xxxx, g_0_y_xxxxzz_xxxy, g_0_y_xxxxzz_xxxz, g_0_y_xxxxzz_xxyy, g_0_y_xxxxzz_xxyz, g_0_y_xxxxzz_xxzz, g_0_y_xxxxzz_xyyy, g_0_y_xxxxzz_xyyz, g_0_y_xxxxzz_xyzz, g_0_y_xxxxzz_xzzz, g_0_y_xxxxzz_yyyy, g_0_y_xxxxzz_yyyz, g_0_y_xxxxzz_yyzz, g_0_y_xxxxzz_yzzz, g_0_y_xxxxzz_zzzz, g_0_y_xxxzz_xxxx, g_0_y_xxxzz_xxxxx, g_0_y_xxxzz_xxxxy, g_0_y_xxxzz_xxxxz, g_0_y_xxxzz_xxxy, g_0_y_xxxzz_xxxyy, g_0_y_xxxzz_xxxyz, g_0_y_xxxzz_xxxz, g_0_y_xxxzz_xxxzz, g_0_y_xxxzz_xxyy, g_0_y_xxxzz_xxyyy, g_0_y_xxxzz_xxyyz, g_0_y_xxxzz_xxyz, g_0_y_xxxzz_xxyzz, g_0_y_xxxzz_xxzz, g_0_y_xxxzz_xxzzz, g_0_y_xxxzz_xyyy, g_0_y_xxxzz_xyyyy, g_0_y_xxxzz_xyyyz, g_0_y_xxxzz_xyyz, g_0_y_xxxzz_xyyzz, g_0_y_xxxzz_xyzz, g_0_y_xxxzz_xyzzz, g_0_y_xxxzz_xzzz, g_0_y_xxxzz_xzzzz, g_0_y_xxxzz_yyyy, g_0_y_xxxzz_yyyz, g_0_y_xxxzz_yyzz, g_0_y_xxxzz_yzzz, g_0_y_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxzz_xxxx[k] = -g_0_y_xxxzz_xxxx[k] * ab_x + g_0_y_xxxzz_xxxxx[k];

                g_0_y_xxxxzz_xxxy[k] = -g_0_y_xxxzz_xxxy[k] * ab_x + g_0_y_xxxzz_xxxxy[k];

                g_0_y_xxxxzz_xxxz[k] = -g_0_y_xxxzz_xxxz[k] * ab_x + g_0_y_xxxzz_xxxxz[k];

                g_0_y_xxxxzz_xxyy[k] = -g_0_y_xxxzz_xxyy[k] * ab_x + g_0_y_xxxzz_xxxyy[k];

                g_0_y_xxxxzz_xxyz[k] = -g_0_y_xxxzz_xxyz[k] * ab_x + g_0_y_xxxzz_xxxyz[k];

                g_0_y_xxxxzz_xxzz[k] = -g_0_y_xxxzz_xxzz[k] * ab_x + g_0_y_xxxzz_xxxzz[k];

                g_0_y_xxxxzz_xyyy[k] = -g_0_y_xxxzz_xyyy[k] * ab_x + g_0_y_xxxzz_xxyyy[k];

                g_0_y_xxxxzz_xyyz[k] = -g_0_y_xxxzz_xyyz[k] * ab_x + g_0_y_xxxzz_xxyyz[k];

                g_0_y_xxxxzz_xyzz[k] = -g_0_y_xxxzz_xyzz[k] * ab_x + g_0_y_xxxzz_xxyzz[k];

                g_0_y_xxxxzz_xzzz[k] = -g_0_y_xxxzz_xzzz[k] * ab_x + g_0_y_xxxzz_xxzzz[k];

                g_0_y_xxxxzz_yyyy[k] = -g_0_y_xxxzz_yyyy[k] * ab_x + g_0_y_xxxzz_xyyyy[k];

                g_0_y_xxxxzz_yyyz[k] = -g_0_y_xxxzz_yyyz[k] * ab_x + g_0_y_xxxzz_xyyyz[k];

                g_0_y_xxxxzz_yyzz[k] = -g_0_y_xxxzz_yyzz[k] * ab_x + g_0_y_xxxzz_xyyzz[k];

                g_0_y_xxxxzz_yzzz[k] = -g_0_y_xxxzz_yzzz[k] * ab_x + g_0_y_xxxzz_xyzzz[k];

                g_0_y_xxxxzz_zzzz[k] = -g_0_y_xxxzz_zzzz[k] * ab_x + g_0_y_xxxzz_xzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyy_xxxx = cbuffer.data(ig_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxxy = cbuffer.data(ig_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxxz = cbuffer.data(ig_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxyy = cbuffer.data(ig_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxyz = cbuffer.data(ig_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxzz = cbuffer.data(ig_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyyy = cbuffer.data(ig_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyyz = cbuffer.data(ig_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyzz = cbuffer.data(ig_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xzzz = cbuffer.data(ig_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyyy = cbuffer.data(ig_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyyz = cbuffer.data(ig_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyzz = cbuffer.data(ig_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yzzz = cbuffer.data(ig_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_y_xxxyyy_zzzz = cbuffer.data(ig_geom_01_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyy_xxxx, g_0_y_xxxyyy_xxxy, g_0_y_xxxyyy_xxxz, g_0_y_xxxyyy_xxyy, g_0_y_xxxyyy_xxyz, g_0_y_xxxyyy_xxzz, g_0_y_xxxyyy_xyyy, g_0_y_xxxyyy_xyyz, g_0_y_xxxyyy_xyzz, g_0_y_xxxyyy_xzzz, g_0_y_xxxyyy_yyyy, g_0_y_xxxyyy_yyyz, g_0_y_xxxyyy_yyzz, g_0_y_xxxyyy_yzzz, g_0_y_xxxyyy_zzzz, g_0_y_xxyyy_xxxx, g_0_y_xxyyy_xxxxx, g_0_y_xxyyy_xxxxy, g_0_y_xxyyy_xxxxz, g_0_y_xxyyy_xxxy, g_0_y_xxyyy_xxxyy, g_0_y_xxyyy_xxxyz, g_0_y_xxyyy_xxxz, g_0_y_xxyyy_xxxzz, g_0_y_xxyyy_xxyy, g_0_y_xxyyy_xxyyy, g_0_y_xxyyy_xxyyz, g_0_y_xxyyy_xxyz, g_0_y_xxyyy_xxyzz, g_0_y_xxyyy_xxzz, g_0_y_xxyyy_xxzzz, g_0_y_xxyyy_xyyy, g_0_y_xxyyy_xyyyy, g_0_y_xxyyy_xyyyz, g_0_y_xxyyy_xyyz, g_0_y_xxyyy_xyyzz, g_0_y_xxyyy_xyzz, g_0_y_xxyyy_xyzzz, g_0_y_xxyyy_xzzz, g_0_y_xxyyy_xzzzz, g_0_y_xxyyy_yyyy, g_0_y_xxyyy_yyyz, g_0_y_xxyyy_yyzz, g_0_y_xxyyy_yzzz, g_0_y_xxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyy_xxxx[k] = -g_0_y_xxyyy_xxxx[k] * ab_x + g_0_y_xxyyy_xxxxx[k];

                g_0_y_xxxyyy_xxxy[k] = -g_0_y_xxyyy_xxxy[k] * ab_x + g_0_y_xxyyy_xxxxy[k];

                g_0_y_xxxyyy_xxxz[k] = -g_0_y_xxyyy_xxxz[k] * ab_x + g_0_y_xxyyy_xxxxz[k];

                g_0_y_xxxyyy_xxyy[k] = -g_0_y_xxyyy_xxyy[k] * ab_x + g_0_y_xxyyy_xxxyy[k];

                g_0_y_xxxyyy_xxyz[k] = -g_0_y_xxyyy_xxyz[k] * ab_x + g_0_y_xxyyy_xxxyz[k];

                g_0_y_xxxyyy_xxzz[k] = -g_0_y_xxyyy_xxzz[k] * ab_x + g_0_y_xxyyy_xxxzz[k];

                g_0_y_xxxyyy_xyyy[k] = -g_0_y_xxyyy_xyyy[k] * ab_x + g_0_y_xxyyy_xxyyy[k];

                g_0_y_xxxyyy_xyyz[k] = -g_0_y_xxyyy_xyyz[k] * ab_x + g_0_y_xxyyy_xxyyz[k];

                g_0_y_xxxyyy_xyzz[k] = -g_0_y_xxyyy_xyzz[k] * ab_x + g_0_y_xxyyy_xxyzz[k];

                g_0_y_xxxyyy_xzzz[k] = -g_0_y_xxyyy_xzzz[k] * ab_x + g_0_y_xxyyy_xxzzz[k];

                g_0_y_xxxyyy_yyyy[k] = -g_0_y_xxyyy_yyyy[k] * ab_x + g_0_y_xxyyy_xyyyy[k];

                g_0_y_xxxyyy_yyyz[k] = -g_0_y_xxyyy_yyyz[k] * ab_x + g_0_y_xxyyy_xyyyz[k];

                g_0_y_xxxyyy_yyzz[k] = -g_0_y_xxyyy_yyzz[k] * ab_x + g_0_y_xxyyy_xyyzz[k];

                g_0_y_xxxyyy_yzzz[k] = -g_0_y_xxyyy_yzzz[k] * ab_x + g_0_y_xxyyy_xyzzz[k];

                g_0_y_xxxyyy_zzzz[k] = -g_0_y_xxyyy_zzzz[k] * ab_x + g_0_y_xxyyy_xzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyz_xxxx = cbuffer.data(ig_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxxy = cbuffer.data(ig_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxxz = cbuffer.data(ig_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxyy = cbuffer.data(ig_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxyz = cbuffer.data(ig_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxzz = cbuffer.data(ig_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyyy = cbuffer.data(ig_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyyz = cbuffer.data(ig_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyzz = cbuffer.data(ig_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xzzz = cbuffer.data(ig_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyyy = cbuffer.data(ig_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyyz = cbuffer.data(ig_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyzz = cbuffer.data(ig_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yzzz = cbuffer.data(ig_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_y_xxxyyz_zzzz = cbuffer.data(ig_geom_01_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyz_xxxx, g_0_y_xxxyyz_xxxy, g_0_y_xxxyyz_xxxz, g_0_y_xxxyyz_xxyy, g_0_y_xxxyyz_xxyz, g_0_y_xxxyyz_xxzz, g_0_y_xxxyyz_xyyy, g_0_y_xxxyyz_xyyz, g_0_y_xxxyyz_xyzz, g_0_y_xxxyyz_xzzz, g_0_y_xxxyyz_yyyy, g_0_y_xxxyyz_yyyz, g_0_y_xxxyyz_yyzz, g_0_y_xxxyyz_yzzz, g_0_y_xxxyyz_zzzz, g_0_y_xxyyz_xxxx, g_0_y_xxyyz_xxxxx, g_0_y_xxyyz_xxxxy, g_0_y_xxyyz_xxxxz, g_0_y_xxyyz_xxxy, g_0_y_xxyyz_xxxyy, g_0_y_xxyyz_xxxyz, g_0_y_xxyyz_xxxz, g_0_y_xxyyz_xxxzz, g_0_y_xxyyz_xxyy, g_0_y_xxyyz_xxyyy, g_0_y_xxyyz_xxyyz, g_0_y_xxyyz_xxyz, g_0_y_xxyyz_xxyzz, g_0_y_xxyyz_xxzz, g_0_y_xxyyz_xxzzz, g_0_y_xxyyz_xyyy, g_0_y_xxyyz_xyyyy, g_0_y_xxyyz_xyyyz, g_0_y_xxyyz_xyyz, g_0_y_xxyyz_xyyzz, g_0_y_xxyyz_xyzz, g_0_y_xxyyz_xyzzz, g_0_y_xxyyz_xzzz, g_0_y_xxyyz_xzzzz, g_0_y_xxyyz_yyyy, g_0_y_xxyyz_yyyz, g_0_y_xxyyz_yyzz, g_0_y_xxyyz_yzzz, g_0_y_xxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyz_xxxx[k] = -g_0_y_xxyyz_xxxx[k] * ab_x + g_0_y_xxyyz_xxxxx[k];

                g_0_y_xxxyyz_xxxy[k] = -g_0_y_xxyyz_xxxy[k] * ab_x + g_0_y_xxyyz_xxxxy[k];

                g_0_y_xxxyyz_xxxz[k] = -g_0_y_xxyyz_xxxz[k] * ab_x + g_0_y_xxyyz_xxxxz[k];

                g_0_y_xxxyyz_xxyy[k] = -g_0_y_xxyyz_xxyy[k] * ab_x + g_0_y_xxyyz_xxxyy[k];

                g_0_y_xxxyyz_xxyz[k] = -g_0_y_xxyyz_xxyz[k] * ab_x + g_0_y_xxyyz_xxxyz[k];

                g_0_y_xxxyyz_xxzz[k] = -g_0_y_xxyyz_xxzz[k] * ab_x + g_0_y_xxyyz_xxxzz[k];

                g_0_y_xxxyyz_xyyy[k] = -g_0_y_xxyyz_xyyy[k] * ab_x + g_0_y_xxyyz_xxyyy[k];

                g_0_y_xxxyyz_xyyz[k] = -g_0_y_xxyyz_xyyz[k] * ab_x + g_0_y_xxyyz_xxyyz[k];

                g_0_y_xxxyyz_xyzz[k] = -g_0_y_xxyyz_xyzz[k] * ab_x + g_0_y_xxyyz_xxyzz[k];

                g_0_y_xxxyyz_xzzz[k] = -g_0_y_xxyyz_xzzz[k] * ab_x + g_0_y_xxyyz_xxzzz[k];

                g_0_y_xxxyyz_yyyy[k] = -g_0_y_xxyyz_yyyy[k] * ab_x + g_0_y_xxyyz_xyyyy[k];

                g_0_y_xxxyyz_yyyz[k] = -g_0_y_xxyyz_yyyz[k] * ab_x + g_0_y_xxyyz_xyyyz[k];

                g_0_y_xxxyyz_yyzz[k] = -g_0_y_xxyyz_yyzz[k] * ab_x + g_0_y_xxyyz_xyyzz[k];

                g_0_y_xxxyyz_yzzz[k] = -g_0_y_xxyyz_yzzz[k] * ab_x + g_0_y_xxyyz_xyzzz[k];

                g_0_y_xxxyyz_zzzz[k] = -g_0_y_xxyyz_zzzz[k] * ab_x + g_0_y_xxyyz_xzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyzz_xxxx = cbuffer.data(ig_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxxy = cbuffer.data(ig_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxxz = cbuffer.data(ig_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxyy = cbuffer.data(ig_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxyz = cbuffer.data(ig_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxzz = cbuffer.data(ig_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyyy = cbuffer.data(ig_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyyz = cbuffer.data(ig_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyzz = cbuffer.data(ig_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xzzz = cbuffer.data(ig_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyyy = cbuffer.data(ig_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyyz = cbuffer.data(ig_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyzz = cbuffer.data(ig_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yzzz = cbuffer.data(ig_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_xxxyzz_zzzz = cbuffer.data(ig_geom_01_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyzz_xxxx, g_0_y_xxxyzz_xxxy, g_0_y_xxxyzz_xxxz, g_0_y_xxxyzz_xxyy, g_0_y_xxxyzz_xxyz, g_0_y_xxxyzz_xxzz, g_0_y_xxxyzz_xyyy, g_0_y_xxxyzz_xyyz, g_0_y_xxxyzz_xyzz, g_0_y_xxxyzz_xzzz, g_0_y_xxxyzz_yyyy, g_0_y_xxxyzz_yyyz, g_0_y_xxxyzz_yyzz, g_0_y_xxxyzz_yzzz, g_0_y_xxxyzz_zzzz, g_0_y_xxyzz_xxxx, g_0_y_xxyzz_xxxxx, g_0_y_xxyzz_xxxxy, g_0_y_xxyzz_xxxxz, g_0_y_xxyzz_xxxy, g_0_y_xxyzz_xxxyy, g_0_y_xxyzz_xxxyz, g_0_y_xxyzz_xxxz, g_0_y_xxyzz_xxxzz, g_0_y_xxyzz_xxyy, g_0_y_xxyzz_xxyyy, g_0_y_xxyzz_xxyyz, g_0_y_xxyzz_xxyz, g_0_y_xxyzz_xxyzz, g_0_y_xxyzz_xxzz, g_0_y_xxyzz_xxzzz, g_0_y_xxyzz_xyyy, g_0_y_xxyzz_xyyyy, g_0_y_xxyzz_xyyyz, g_0_y_xxyzz_xyyz, g_0_y_xxyzz_xyyzz, g_0_y_xxyzz_xyzz, g_0_y_xxyzz_xyzzz, g_0_y_xxyzz_xzzz, g_0_y_xxyzz_xzzzz, g_0_y_xxyzz_yyyy, g_0_y_xxyzz_yyyz, g_0_y_xxyzz_yyzz, g_0_y_xxyzz_yzzz, g_0_y_xxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyzz_xxxx[k] = -g_0_y_xxyzz_xxxx[k] * ab_x + g_0_y_xxyzz_xxxxx[k];

                g_0_y_xxxyzz_xxxy[k] = -g_0_y_xxyzz_xxxy[k] * ab_x + g_0_y_xxyzz_xxxxy[k];

                g_0_y_xxxyzz_xxxz[k] = -g_0_y_xxyzz_xxxz[k] * ab_x + g_0_y_xxyzz_xxxxz[k];

                g_0_y_xxxyzz_xxyy[k] = -g_0_y_xxyzz_xxyy[k] * ab_x + g_0_y_xxyzz_xxxyy[k];

                g_0_y_xxxyzz_xxyz[k] = -g_0_y_xxyzz_xxyz[k] * ab_x + g_0_y_xxyzz_xxxyz[k];

                g_0_y_xxxyzz_xxzz[k] = -g_0_y_xxyzz_xxzz[k] * ab_x + g_0_y_xxyzz_xxxzz[k];

                g_0_y_xxxyzz_xyyy[k] = -g_0_y_xxyzz_xyyy[k] * ab_x + g_0_y_xxyzz_xxyyy[k];

                g_0_y_xxxyzz_xyyz[k] = -g_0_y_xxyzz_xyyz[k] * ab_x + g_0_y_xxyzz_xxyyz[k];

                g_0_y_xxxyzz_xyzz[k] = -g_0_y_xxyzz_xyzz[k] * ab_x + g_0_y_xxyzz_xxyzz[k];

                g_0_y_xxxyzz_xzzz[k] = -g_0_y_xxyzz_xzzz[k] * ab_x + g_0_y_xxyzz_xxzzz[k];

                g_0_y_xxxyzz_yyyy[k] = -g_0_y_xxyzz_yyyy[k] * ab_x + g_0_y_xxyzz_xyyyy[k];

                g_0_y_xxxyzz_yyyz[k] = -g_0_y_xxyzz_yyyz[k] * ab_x + g_0_y_xxyzz_xyyyz[k];

                g_0_y_xxxyzz_yyzz[k] = -g_0_y_xxyzz_yyzz[k] * ab_x + g_0_y_xxyzz_xyyzz[k];

                g_0_y_xxxyzz_yzzz[k] = -g_0_y_xxyzz_yzzz[k] * ab_x + g_0_y_xxyzz_xyzzz[k];

                g_0_y_xxxyzz_zzzz[k] = -g_0_y_xxyzz_zzzz[k] * ab_x + g_0_y_xxyzz_xzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzzz_xxxx = cbuffer.data(ig_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxxy = cbuffer.data(ig_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxxz = cbuffer.data(ig_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxyy = cbuffer.data(ig_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxyz = cbuffer.data(ig_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxzz = cbuffer.data(ig_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyyy = cbuffer.data(ig_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyyz = cbuffer.data(ig_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyzz = cbuffer.data(ig_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xzzz = cbuffer.data(ig_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyyy = cbuffer.data(ig_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyyz = cbuffer.data(ig_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyzz = cbuffer.data(ig_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yzzz = cbuffer.data(ig_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_y_xxxzzz_zzzz = cbuffer.data(ig_geom_01_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzzz_xxxx, g_0_y_xxxzzz_xxxy, g_0_y_xxxzzz_xxxz, g_0_y_xxxzzz_xxyy, g_0_y_xxxzzz_xxyz, g_0_y_xxxzzz_xxzz, g_0_y_xxxzzz_xyyy, g_0_y_xxxzzz_xyyz, g_0_y_xxxzzz_xyzz, g_0_y_xxxzzz_xzzz, g_0_y_xxxzzz_yyyy, g_0_y_xxxzzz_yyyz, g_0_y_xxxzzz_yyzz, g_0_y_xxxzzz_yzzz, g_0_y_xxxzzz_zzzz, g_0_y_xxzzz_xxxx, g_0_y_xxzzz_xxxxx, g_0_y_xxzzz_xxxxy, g_0_y_xxzzz_xxxxz, g_0_y_xxzzz_xxxy, g_0_y_xxzzz_xxxyy, g_0_y_xxzzz_xxxyz, g_0_y_xxzzz_xxxz, g_0_y_xxzzz_xxxzz, g_0_y_xxzzz_xxyy, g_0_y_xxzzz_xxyyy, g_0_y_xxzzz_xxyyz, g_0_y_xxzzz_xxyz, g_0_y_xxzzz_xxyzz, g_0_y_xxzzz_xxzz, g_0_y_xxzzz_xxzzz, g_0_y_xxzzz_xyyy, g_0_y_xxzzz_xyyyy, g_0_y_xxzzz_xyyyz, g_0_y_xxzzz_xyyz, g_0_y_xxzzz_xyyzz, g_0_y_xxzzz_xyzz, g_0_y_xxzzz_xyzzz, g_0_y_xxzzz_xzzz, g_0_y_xxzzz_xzzzz, g_0_y_xxzzz_yyyy, g_0_y_xxzzz_yyyz, g_0_y_xxzzz_yyzz, g_0_y_xxzzz_yzzz, g_0_y_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzzz_xxxx[k] = -g_0_y_xxzzz_xxxx[k] * ab_x + g_0_y_xxzzz_xxxxx[k];

                g_0_y_xxxzzz_xxxy[k] = -g_0_y_xxzzz_xxxy[k] * ab_x + g_0_y_xxzzz_xxxxy[k];

                g_0_y_xxxzzz_xxxz[k] = -g_0_y_xxzzz_xxxz[k] * ab_x + g_0_y_xxzzz_xxxxz[k];

                g_0_y_xxxzzz_xxyy[k] = -g_0_y_xxzzz_xxyy[k] * ab_x + g_0_y_xxzzz_xxxyy[k];

                g_0_y_xxxzzz_xxyz[k] = -g_0_y_xxzzz_xxyz[k] * ab_x + g_0_y_xxzzz_xxxyz[k];

                g_0_y_xxxzzz_xxzz[k] = -g_0_y_xxzzz_xxzz[k] * ab_x + g_0_y_xxzzz_xxxzz[k];

                g_0_y_xxxzzz_xyyy[k] = -g_0_y_xxzzz_xyyy[k] * ab_x + g_0_y_xxzzz_xxyyy[k];

                g_0_y_xxxzzz_xyyz[k] = -g_0_y_xxzzz_xyyz[k] * ab_x + g_0_y_xxzzz_xxyyz[k];

                g_0_y_xxxzzz_xyzz[k] = -g_0_y_xxzzz_xyzz[k] * ab_x + g_0_y_xxzzz_xxyzz[k];

                g_0_y_xxxzzz_xzzz[k] = -g_0_y_xxzzz_xzzz[k] * ab_x + g_0_y_xxzzz_xxzzz[k];

                g_0_y_xxxzzz_yyyy[k] = -g_0_y_xxzzz_yyyy[k] * ab_x + g_0_y_xxzzz_xyyyy[k];

                g_0_y_xxxzzz_yyyz[k] = -g_0_y_xxzzz_yyyz[k] * ab_x + g_0_y_xxzzz_xyyyz[k];

                g_0_y_xxxzzz_yyzz[k] = -g_0_y_xxzzz_yyzz[k] * ab_x + g_0_y_xxzzz_xyyzz[k];

                g_0_y_xxxzzz_yzzz[k] = -g_0_y_xxzzz_yzzz[k] * ab_x + g_0_y_xxzzz_xyzzz[k];

                g_0_y_xxxzzz_zzzz[k] = -g_0_y_xxzzz_zzzz[k] * ab_x + g_0_y_xxzzz_xzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyy_xxxx = cbuffer.data(ig_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxxy = cbuffer.data(ig_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxxz = cbuffer.data(ig_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxyy = cbuffer.data(ig_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxyz = cbuffer.data(ig_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxzz = cbuffer.data(ig_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyyy = cbuffer.data(ig_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyyz = cbuffer.data(ig_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyzz = cbuffer.data(ig_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xzzz = cbuffer.data(ig_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyyy = cbuffer.data(ig_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyyz = cbuffer.data(ig_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyzz = cbuffer.data(ig_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yzzz = cbuffer.data(ig_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_y_xxyyyy_zzzz = cbuffer.data(ig_geom_01_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyy_xxxx, g_0_y_xxyyyy_xxxy, g_0_y_xxyyyy_xxxz, g_0_y_xxyyyy_xxyy, g_0_y_xxyyyy_xxyz, g_0_y_xxyyyy_xxzz, g_0_y_xxyyyy_xyyy, g_0_y_xxyyyy_xyyz, g_0_y_xxyyyy_xyzz, g_0_y_xxyyyy_xzzz, g_0_y_xxyyyy_yyyy, g_0_y_xxyyyy_yyyz, g_0_y_xxyyyy_yyzz, g_0_y_xxyyyy_yzzz, g_0_y_xxyyyy_zzzz, g_0_y_xyyyy_xxxx, g_0_y_xyyyy_xxxxx, g_0_y_xyyyy_xxxxy, g_0_y_xyyyy_xxxxz, g_0_y_xyyyy_xxxy, g_0_y_xyyyy_xxxyy, g_0_y_xyyyy_xxxyz, g_0_y_xyyyy_xxxz, g_0_y_xyyyy_xxxzz, g_0_y_xyyyy_xxyy, g_0_y_xyyyy_xxyyy, g_0_y_xyyyy_xxyyz, g_0_y_xyyyy_xxyz, g_0_y_xyyyy_xxyzz, g_0_y_xyyyy_xxzz, g_0_y_xyyyy_xxzzz, g_0_y_xyyyy_xyyy, g_0_y_xyyyy_xyyyy, g_0_y_xyyyy_xyyyz, g_0_y_xyyyy_xyyz, g_0_y_xyyyy_xyyzz, g_0_y_xyyyy_xyzz, g_0_y_xyyyy_xyzzz, g_0_y_xyyyy_xzzz, g_0_y_xyyyy_xzzzz, g_0_y_xyyyy_yyyy, g_0_y_xyyyy_yyyz, g_0_y_xyyyy_yyzz, g_0_y_xyyyy_yzzz, g_0_y_xyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyy_xxxx[k] = -g_0_y_xyyyy_xxxx[k] * ab_x + g_0_y_xyyyy_xxxxx[k];

                g_0_y_xxyyyy_xxxy[k] = -g_0_y_xyyyy_xxxy[k] * ab_x + g_0_y_xyyyy_xxxxy[k];

                g_0_y_xxyyyy_xxxz[k] = -g_0_y_xyyyy_xxxz[k] * ab_x + g_0_y_xyyyy_xxxxz[k];

                g_0_y_xxyyyy_xxyy[k] = -g_0_y_xyyyy_xxyy[k] * ab_x + g_0_y_xyyyy_xxxyy[k];

                g_0_y_xxyyyy_xxyz[k] = -g_0_y_xyyyy_xxyz[k] * ab_x + g_0_y_xyyyy_xxxyz[k];

                g_0_y_xxyyyy_xxzz[k] = -g_0_y_xyyyy_xxzz[k] * ab_x + g_0_y_xyyyy_xxxzz[k];

                g_0_y_xxyyyy_xyyy[k] = -g_0_y_xyyyy_xyyy[k] * ab_x + g_0_y_xyyyy_xxyyy[k];

                g_0_y_xxyyyy_xyyz[k] = -g_0_y_xyyyy_xyyz[k] * ab_x + g_0_y_xyyyy_xxyyz[k];

                g_0_y_xxyyyy_xyzz[k] = -g_0_y_xyyyy_xyzz[k] * ab_x + g_0_y_xyyyy_xxyzz[k];

                g_0_y_xxyyyy_xzzz[k] = -g_0_y_xyyyy_xzzz[k] * ab_x + g_0_y_xyyyy_xxzzz[k];

                g_0_y_xxyyyy_yyyy[k] = -g_0_y_xyyyy_yyyy[k] * ab_x + g_0_y_xyyyy_xyyyy[k];

                g_0_y_xxyyyy_yyyz[k] = -g_0_y_xyyyy_yyyz[k] * ab_x + g_0_y_xyyyy_xyyyz[k];

                g_0_y_xxyyyy_yyzz[k] = -g_0_y_xyyyy_yyzz[k] * ab_x + g_0_y_xyyyy_xyyzz[k];

                g_0_y_xxyyyy_yzzz[k] = -g_0_y_xyyyy_yzzz[k] * ab_x + g_0_y_xyyyy_xyzzz[k];

                g_0_y_xxyyyy_zzzz[k] = -g_0_y_xyyyy_zzzz[k] * ab_x + g_0_y_xyyyy_xzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyz_xxxx = cbuffer.data(ig_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxxy = cbuffer.data(ig_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxxz = cbuffer.data(ig_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxyy = cbuffer.data(ig_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxyz = cbuffer.data(ig_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxzz = cbuffer.data(ig_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyyy = cbuffer.data(ig_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyyz = cbuffer.data(ig_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyzz = cbuffer.data(ig_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xzzz = cbuffer.data(ig_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyyy = cbuffer.data(ig_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyyz = cbuffer.data(ig_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyzz = cbuffer.data(ig_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yzzz = cbuffer.data(ig_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_xxyyyz_zzzz = cbuffer.data(ig_geom_01_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyz_xxxx, g_0_y_xxyyyz_xxxy, g_0_y_xxyyyz_xxxz, g_0_y_xxyyyz_xxyy, g_0_y_xxyyyz_xxyz, g_0_y_xxyyyz_xxzz, g_0_y_xxyyyz_xyyy, g_0_y_xxyyyz_xyyz, g_0_y_xxyyyz_xyzz, g_0_y_xxyyyz_xzzz, g_0_y_xxyyyz_yyyy, g_0_y_xxyyyz_yyyz, g_0_y_xxyyyz_yyzz, g_0_y_xxyyyz_yzzz, g_0_y_xxyyyz_zzzz, g_0_y_xyyyz_xxxx, g_0_y_xyyyz_xxxxx, g_0_y_xyyyz_xxxxy, g_0_y_xyyyz_xxxxz, g_0_y_xyyyz_xxxy, g_0_y_xyyyz_xxxyy, g_0_y_xyyyz_xxxyz, g_0_y_xyyyz_xxxz, g_0_y_xyyyz_xxxzz, g_0_y_xyyyz_xxyy, g_0_y_xyyyz_xxyyy, g_0_y_xyyyz_xxyyz, g_0_y_xyyyz_xxyz, g_0_y_xyyyz_xxyzz, g_0_y_xyyyz_xxzz, g_0_y_xyyyz_xxzzz, g_0_y_xyyyz_xyyy, g_0_y_xyyyz_xyyyy, g_0_y_xyyyz_xyyyz, g_0_y_xyyyz_xyyz, g_0_y_xyyyz_xyyzz, g_0_y_xyyyz_xyzz, g_0_y_xyyyz_xyzzz, g_0_y_xyyyz_xzzz, g_0_y_xyyyz_xzzzz, g_0_y_xyyyz_yyyy, g_0_y_xyyyz_yyyz, g_0_y_xyyyz_yyzz, g_0_y_xyyyz_yzzz, g_0_y_xyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyz_xxxx[k] = -g_0_y_xyyyz_xxxx[k] * ab_x + g_0_y_xyyyz_xxxxx[k];

                g_0_y_xxyyyz_xxxy[k] = -g_0_y_xyyyz_xxxy[k] * ab_x + g_0_y_xyyyz_xxxxy[k];

                g_0_y_xxyyyz_xxxz[k] = -g_0_y_xyyyz_xxxz[k] * ab_x + g_0_y_xyyyz_xxxxz[k];

                g_0_y_xxyyyz_xxyy[k] = -g_0_y_xyyyz_xxyy[k] * ab_x + g_0_y_xyyyz_xxxyy[k];

                g_0_y_xxyyyz_xxyz[k] = -g_0_y_xyyyz_xxyz[k] * ab_x + g_0_y_xyyyz_xxxyz[k];

                g_0_y_xxyyyz_xxzz[k] = -g_0_y_xyyyz_xxzz[k] * ab_x + g_0_y_xyyyz_xxxzz[k];

                g_0_y_xxyyyz_xyyy[k] = -g_0_y_xyyyz_xyyy[k] * ab_x + g_0_y_xyyyz_xxyyy[k];

                g_0_y_xxyyyz_xyyz[k] = -g_0_y_xyyyz_xyyz[k] * ab_x + g_0_y_xyyyz_xxyyz[k];

                g_0_y_xxyyyz_xyzz[k] = -g_0_y_xyyyz_xyzz[k] * ab_x + g_0_y_xyyyz_xxyzz[k];

                g_0_y_xxyyyz_xzzz[k] = -g_0_y_xyyyz_xzzz[k] * ab_x + g_0_y_xyyyz_xxzzz[k];

                g_0_y_xxyyyz_yyyy[k] = -g_0_y_xyyyz_yyyy[k] * ab_x + g_0_y_xyyyz_xyyyy[k];

                g_0_y_xxyyyz_yyyz[k] = -g_0_y_xyyyz_yyyz[k] * ab_x + g_0_y_xyyyz_xyyyz[k];

                g_0_y_xxyyyz_yyzz[k] = -g_0_y_xyyyz_yyzz[k] * ab_x + g_0_y_xyyyz_xyyzz[k];

                g_0_y_xxyyyz_yzzz[k] = -g_0_y_xyyyz_yzzz[k] * ab_x + g_0_y_xyyyz_xyzzz[k];

                g_0_y_xxyyyz_zzzz[k] = -g_0_y_xyyyz_zzzz[k] * ab_x + g_0_y_xyyyz_xzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyzz_xxxx = cbuffer.data(ig_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxxy = cbuffer.data(ig_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxxz = cbuffer.data(ig_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxyy = cbuffer.data(ig_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxyz = cbuffer.data(ig_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxzz = cbuffer.data(ig_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyyy = cbuffer.data(ig_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyyz = cbuffer.data(ig_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyzz = cbuffer.data(ig_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xzzz = cbuffer.data(ig_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyyy = cbuffer.data(ig_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyyz = cbuffer.data(ig_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyzz = cbuffer.data(ig_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yzzz = cbuffer.data(ig_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_xxyyzz_zzzz = cbuffer.data(ig_geom_01_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyzz_xxxx, g_0_y_xxyyzz_xxxy, g_0_y_xxyyzz_xxxz, g_0_y_xxyyzz_xxyy, g_0_y_xxyyzz_xxyz, g_0_y_xxyyzz_xxzz, g_0_y_xxyyzz_xyyy, g_0_y_xxyyzz_xyyz, g_0_y_xxyyzz_xyzz, g_0_y_xxyyzz_xzzz, g_0_y_xxyyzz_yyyy, g_0_y_xxyyzz_yyyz, g_0_y_xxyyzz_yyzz, g_0_y_xxyyzz_yzzz, g_0_y_xxyyzz_zzzz, g_0_y_xyyzz_xxxx, g_0_y_xyyzz_xxxxx, g_0_y_xyyzz_xxxxy, g_0_y_xyyzz_xxxxz, g_0_y_xyyzz_xxxy, g_0_y_xyyzz_xxxyy, g_0_y_xyyzz_xxxyz, g_0_y_xyyzz_xxxz, g_0_y_xyyzz_xxxzz, g_0_y_xyyzz_xxyy, g_0_y_xyyzz_xxyyy, g_0_y_xyyzz_xxyyz, g_0_y_xyyzz_xxyz, g_0_y_xyyzz_xxyzz, g_0_y_xyyzz_xxzz, g_0_y_xyyzz_xxzzz, g_0_y_xyyzz_xyyy, g_0_y_xyyzz_xyyyy, g_0_y_xyyzz_xyyyz, g_0_y_xyyzz_xyyz, g_0_y_xyyzz_xyyzz, g_0_y_xyyzz_xyzz, g_0_y_xyyzz_xyzzz, g_0_y_xyyzz_xzzz, g_0_y_xyyzz_xzzzz, g_0_y_xyyzz_yyyy, g_0_y_xyyzz_yyyz, g_0_y_xyyzz_yyzz, g_0_y_xyyzz_yzzz, g_0_y_xyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyzz_xxxx[k] = -g_0_y_xyyzz_xxxx[k] * ab_x + g_0_y_xyyzz_xxxxx[k];

                g_0_y_xxyyzz_xxxy[k] = -g_0_y_xyyzz_xxxy[k] * ab_x + g_0_y_xyyzz_xxxxy[k];

                g_0_y_xxyyzz_xxxz[k] = -g_0_y_xyyzz_xxxz[k] * ab_x + g_0_y_xyyzz_xxxxz[k];

                g_0_y_xxyyzz_xxyy[k] = -g_0_y_xyyzz_xxyy[k] * ab_x + g_0_y_xyyzz_xxxyy[k];

                g_0_y_xxyyzz_xxyz[k] = -g_0_y_xyyzz_xxyz[k] * ab_x + g_0_y_xyyzz_xxxyz[k];

                g_0_y_xxyyzz_xxzz[k] = -g_0_y_xyyzz_xxzz[k] * ab_x + g_0_y_xyyzz_xxxzz[k];

                g_0_y_xxyyzz_xyyy[k] = -g_0_y_xyyzz_xyyy[k] * ab_x + g_0_y_xyyzz_xxyyy[k];

                g_0_y_xxyyzz_xyyz[k] = -g_0_y_xyyzz_xyyz[k] * ab_x + g_0_y_xyyzz_xxyyz[k];

                g_0_y_xxyyzz_xyzz[k] = -g_0_y_xyyzz_xyzz[k] * ab_x + g_0_y_xyyzz_xxyzz[k];

                g_0_y_xxyyzz_xzzz[k] = -g_0_y_xyyzz_xzzz[k] * ab_x + g_0_y_xyyzz_xxzzz[k];

                g_0_y_xxyyzz_yyyy[k] = -g_0_y_xyyzz_yyyy[k] * ab_x + g_0_y_xyyzz_xyyyy[k];

                g_0_y_xxyyzz_yyyz[k] = -g_0_y_xyyzz_yyyz[k] * ab_x + g_0_y_xyyzz_xyyyz[k];

                g_0_y_xxyyzz_yyzz[k] = -g_0_y_xyyzz_yyzz[k] * ab_x + g_0_y_xyyzz_xyyzz[k];

                g_0_y_xxyyzz_yzzz[k] = -g_0_y_xyyzz_yzzz[k] * ab_x + g_0_y_xyyzz_xyzzz[k];

                g_0_y_xxyyzz_zzzz[k] = -g_0_y_xyyzz_zzzz[k] * ab_x + g_0_y_xyyzz_xzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzzz_xxxx = cbuffer.data(ig_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxxy = cbuffer.data(ig_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxxz = cbuffer.data(ig_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxyy = cbuffer.data(ig_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxyz = cbuffer.data(ig_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxzz = cbuffer.data(ig_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyyy = cbuffer.data(ig_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyyz = cbuffer.data(ig_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyzz = cbuffer.data(ig_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xzzz = cbuffer.data(ig_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyyy = cbuffer.data(ig_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyyz = cbuffer.data(ig_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyzz = cbuffer.data(ig_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yzzz = cbuffer.data(ig_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_xxyzzz_zzzz = cbuffer.data(ig_geom_01_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzzz_xxxx, g_0_y_xxyzzz_xxxy, g_0_y_xxyzzz_xxxz, g_0_y_xxyzzz_xxyy, g_0_y_xxyzzz_xxyz, g_0_y_xxyzzz_xxzz, g_0_y_xxyzzz_xyyy, g_0_y_xxyzzz_xyyz, g_0_y_xxyzzz_xyzz, g_0_y_xxyzzz_xzzz, g_0_y_xxyzzz_yyyy, g_0_y_xxyzzz_yyyz, g_0_y_xxyzzz_yyzz, g_0_y_xxyzzz_yzzz, g_0_y_xxyzzz_zzzz, g_0_y_xyzzz_xxxx, g_0_y_xyzzz_xxxxx, g_0_y_xyzzz_xxxxy, g_0_y_xyzzz_xxxxz, g_0_y_xyzzz_xxxy, g_0_y_xyzzz_xxxyy, g_0_y_xyzzz_xxxyz, g_0_y_xyzzz_xxxz, g_0_y_xyzzz_xxxzz, g_0_y_xyzzz_xxyy, g_0_y_xyzzz_xxyyy, g_0_y_xyzzz_xxyyz, g_0_y_xyzzz_xxyz, g_0_y_xyzzz_xxyzz, g_0_y_xyzzz_xxzz, g_0_y_xyzzz_xxzzz, g_0_y_xyzzz_xyyy, g_0_y_xyzzz_xyyyy, g_0_y_xyzzz_xyyyz, g_0_y_xyzzz_xyyz, g_0_y_xyzzz_xyyzz, g_0_y_xyzzz_xyzz, g_0_y_xyzzz_xyzzz, g_0_y_xyzzz_xzzz, g_0_y_xyzzz_xzzzz, g_0_y_xyzzz_yyyy, g_0_y_xyzzz_yyyz, g_0_y_xyzzz_yyzz, g_0_y_xyzzz_yzzz, g_0_y_xyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzzz_xxxx[k] = -g_0_y_xyzzz_xxxx[k] * ab_x + g_0_y_xyzzz_xxxxx[k];

                g_0_y_xxyzzz_xxxy[k] = -g_0_y_xyzzz_xxxy[k] * ab_x + g_0_y_xyzzz_xxxxy[k];

                g_0_y_xxyzzz_xxxz[k] = -g_0_y_xyzzz_xxxz[k] * ab_x + g_0_y_xyzzz_xxxxz[k];

                g_0_y_xxyzzz_xxyy[k] = -g_0_y_xyzzz_xxyy[k] * ab_x + g_0_y_xyzzz_xxxyy[k];

                g_0_y_xxyzzz_xxyz[k] = -g_0_y_xyzzz_xxyz[k] * ab_x + g_0_y_xyzzz_xxxyz[k];

                g_0_y_xxyzzz_xxzz[k] = -g_0_y_xyzzz_xxzz[k] * ab_x + g_0_y_xyzzz_xxxzz[k];

                g_0_y_xxyzzz_xyyy[k] = -g_0_y_xyzzz_xyyy[k] * ab_x + g_0_y_xyzzz_xxyyy[k];

                g_0_y_xxyzzz_xyyz[k] = -g_0_y_xyzzz_xyyz[k] * ab_x + g_0_y_xyzzz_xxyyz[k];

                g_0_y_xxyzzz_xyzz[k] = -g_0_y_xyzzz_xyzz[k] * ab_x + g_0_y_xyzzz_xxyzz[k];

                g_0_y_xxyzzz_xzzz[k] = -g_0_y_xyzzz_xzzz[k] * ab_x + g_0_y_xyzzz_xxzzz[k];

                g_0_y_xxyzzz_yyyy[k] = -g_0_y_xyzzz_yyyy[k] * ab_x + g_0_y_xyzzz_xyyyy[k];

                g_0_y_xxyzzz_yyyz[k] = -g_0_y_xyzzz_yyyz[k] * ab_x + g_0_y_xyzzz_xyyyz[k];

                g_0_y_xxyzzz_yyzz[k] = -g_0_y_xyzzz_yyzz[k] * ab_x + g_0_y_xyzzz_xyyzz[k];

                g_0_y_xxyzzz_yzzz[k] = -g_0_y_xyzzz_yzzz[k] * ab_x + g_0_y_xyzzz_xyzzz[k];

                g_0_y_xxyzzz_zzzz[k] = -g_0_y_xyzzz_zzzz[k] * ab_x + g_0_y_xyzzz_xzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzzz_xxxx = cbuffer.data(ig_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxxy = cbuffer.data(ig_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxxz = cbuffer.data(ig_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxyy = cbuffer.data(ig_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxyz = cbuffer.data(ig_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxzz = cbuffer.data(ig_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyyy = cbuffer.data(ig_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyyz = cbuffer.data(ig_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyzz = cbuffer.data(ig_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xzzz = cbuffer.data(ig_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyyy = cbuffer.data(ig_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyyz = cbuffer.data(ig_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyzz = cbuffer.data(ig_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yzzz = cbuffer.data(ig_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_y_xxzzzz_zzzz = cbuffer.data(ig_geom_01_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzzz_xxxx, g_0_y_xxzzzz_xxxy, g_0_y_xxzzzz_xxxz, g_0_y_xxzzzz_xxyy, g_0_y_xxzzzz_xxyz, g_0_y_xxzzzz_xxzz, g_0_y_xxzzzz_xyyy, g_0_y_xxzzzz_xyyz, g_0_y_xxzzzz_xyzz, g_0_y_xxzzzz_xzzz, g_0_y_xxzzzz_yyyy, g_0_y_xxzzzz_yyyz, g_0_y_xxzzzz_yyzz, g_0_y_xxzzzz_yzzz, g_0_y_xxzzzz_zzzz, g_0_y_xzzzz_xxxx, g_0_y_xzzzz_xxxxx, g_0_y_xzzzz_xxxxy, g_0_y_xzzzz_xxxxz, g_0_y_xzzzz_xxxy, g_0_y_xzzzz_xxxyy, g_0_y_xzzzz_xxxyz, g_0_y_xzzzz_xxxz, g_0_y_xzzzz_xxxzz, g_0_y_xzzzz_xxyy, g_0_y_xzzzz_xxyyy, g_0_y_xzzzz_xxyyz, g_0_y_xzzzz_xxyz, g_0_y_xzzzz_xxyzz, g_0_y_xzzzz_xxzz, g_0_y_xzzzz_xxzzz, g_0_y_xzzzz_xyyy, g_0_y_xzzzz_xyyyy, g_0_y_xzzzz_xyyyz, g_0_y_xzzzz_xyyz, g_0_y_xzzzz_xyyzz, g_0_y_xzzzz_xyzz, g_0_y_xzzzz_xyzzz, g_0_y_xzzzz_xzzz, g_0_y_xzzzz_xzzzz, g_0_y_xzzzz_yyyy, g_0_y_xzzzz_yyyz, g_0_y_xzzzz_yyzz, g_0_y_xzzzz_yzzz, g_0_y_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzzz_xxxx[k] = -g_0_y_xzzzz_xxxx[k] * ab_x + g_0_y_xzzzz_xxxxx[k];

                g_0_y_xxzzzz_xxxy[k] = -g_0_y_xzzzz_xxxy[k] * ab_x + g_0_y_xzzzz_xxxxy[k];

                g_0_y_xxzzzz_xxxz[k] = -g_0_y_xzzzz_xxxz[k] * ab_x + g_0_y_xzzzz_xxxxz[k];

                g_0_y_xxzzzz_xxyy[k] = -g_0_y_xzzzz_xxyy[k] * ab_x + g_0_y_xzzzz_xxxyy[k];

                g_0_y_xxzzzz_xxyz[k] = -g_0_y_xzzzz_xxyz[k] * ab_x + g_0_y_xzzzz_xxxyz[k];

                g_0_y_xxzzzz_xxzz[k] = -g_0_y_xzzzz_xxzz[k] * ab_x + g_0_y_xzzzz_xxxzz[k];

                g_0_y_xxzzzz_xyyy[k] = -g_0_y_xzzzz_xyyy[k] * ab_x + g_0_y_xzzzz_xxyyy[k];

                g_0_y_xxzzzz_xyyz[k] = -g_0_y_xzzzz_xyyz[k] * ab_x + g_0_y_xzzzz_xxyyz[k];

                g_0_y_xxzzzz_xyzz[k] = -g_0_y_xzzzz_xyzz[k] * ab_x + g_0_y_xzzzz_xxyzz[k];

                g_0_y_xxzzzz_xzzz[k] = -g_0_y_xzzzz_xzzz[k] * ab_x + g_0_y_xzzzz_xxzzz[k];

                g_0_y_xxzzzz_yyyy[k] = -g_0_y_xzzzz_yyyy[k] * ab_x + g_0_y_xzzzz_xyyyy[k];

                g_0_y_xxzzzz_yyyz[k] = -g_0_y_xzzzz_yyyz[k] * ab_x + g_0_y_xzzzz_xyyyz[k];

                g_0_y_xxzzzz_yyzz[k] = -g_0_y_xzzzz_yyzz[k] * ab_x + g_0_y_xzzzz_xyyzz[k];

                g_0_y_xxzzzz_yzzz[k] = -g_0_y_xzzzz_yzzz[k] * ab_x + g_0_y_xzzzz_xyzzz[k];

                g_0_y_xxzzzz_zzzz[k] = -g_0_y_xzzzz_zzzz[k] * ab_x + g_0_y_xzzzz_xzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyy_xxxx = cbuffer.data(ig_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxxy = cbuffer.data(ig_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxxz = cbuffer.data(ig_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxyy = cbuffer.data(ig_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxyz = cbuffer.data(ig_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxzz = cbuffer.data(ig_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyyy = cbuffer.data(ig_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyyz = cbuffer.data(ig_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyzz = cbuffer.data(ig_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xzzz = cbuffer.data(ig_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyyy = cbuffer.data(ig_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyyz = cbuffer.data(ig_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyzz = cbuffer.data(ig_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yzzz = cbuffer.data(ig_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_y_xyyyyy_zzzz = cbuffer.data(ig_geom_01_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyy_xxxx, g_0_y_xyyyyy_xxxy, g_0_y_xyyyyy_xxxz, g_0_y_xyyyyy_xxyy, g_0_y_xyyyyy_xxyz, g_0_y_xyyyyy_xxzz, g_0_y_xyyyyy_xyyy, g_0_y_xyyyyy_xyyz, g_0_y_xyyyyy_xyzz, g_0_y_xyyyyy_xzzz, g_0_y_xyyyyy_yyyy, g_0_y_xyyyyy_yyyz, g_0_y_xyyyyy_yyzz, g_0_y_xyyyyy_yzzz, g_0_y_xyyyyy_zzzz, g_0_y_yyyyy_xxxx, g_0_y_yyyyy_xxxxx, g_0_y_yyyyy_xxxxy, g_0_y_yyyyy_xxxxz, g_0_y_yyyyy_xxxy, g_0_y_yyyyy_xxxyy, g_0_y_yyyyy_xxxyz, g_0_y_yyyyy_xxxz, g_0_y_yyyyy_xxxzz, g_0_y_yyyyy_xxyy, g_0_y_yyyyy_xxyyy, g_0_y_yyyyy_xxyyz, g_0_y_yyyyy_xxyz, g_0_y_yyyyy_xxyzz, g_0_y_yyyyy_xxzz, g_0_y_yyyyy_xxzzz, g_0_y_yyyyy_xyyy, g_0_y_yyyyy_xyyyy, g_0_y_yyyyy_xyyyz, g_0_y_yyyyy_xyyz, g_0_y_yyyyy_xyyzz, g_0_y_yyyyy_xyzz, g_0_y_yyyyy_xyzzz, g_0_y_yyyyy_xzzz, g_0_y_yyyyy_xzzzz, g_0_y_yyyyy_yyyy, g_0_y_yyyyy_yyyz, g_0_y_yyyyy_yyzz, g_0_y_yyyyy_yzzz, g_0_y_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyy_xxxx[k] = -g_0_y_yyyyy_xxxx[k] * ab_x + g_0_y_yyyyy_xxxxx[k];

                g_0_y_xyyyyy_xxxy[k] = -g_0_y_yyyyy_xxxy[k] * ab_x + g_0_y_yyyyy_xxxxy[k];

                g_0_y_xyyyyy_xxxz[k] = -g_0_y_yyyyy_xxxz[k] * ab_x + g_0_y_yyyyy_xxxxz[k];

                g_0_y_xyyyyy_xxyy[k] = -g_0_y_yyyyy_xxyy[k] * ab_x + g_0_y_yyyyy_xxxyy[k];

                g_0_y_xyyyyy_xxyz[k] = -g_0_y_yyyyy_xxyz[k] * ab_x + g_0_y_yyyyy_xxxyz[k];

                g_0_y_xyyyyy_xxzz[k] = -g_0_y_yyyyy_xxzz[k] * ab_x + g_0_y_yyyyy_xxxzz[k];

                g_0_y_xyyyyy_xyyy[k] = -g_0_y_yyyyy_xyyy[k] * ab_x + g_0_y_yyyyy_xxyyy[k];

                g_0_y_xyyyyy_xyyz[k] = -g_0_y_yyyyy_xyyz[k] * ab_x + g_0_y_yyyyy_xxyyz[k];

                g_0_y_xyyyyy_xyzz[k] = -g_0_y_yyyyy_xyzz[k] * ab_x + g_0_y_yyyyy_xxyzz[k];

                g_0_y_xyyyyy_xzzz[k] = -g_0_y_yyyyy_xzzz[k] * ab_x + g_0_y_yyyyy_xxzzz[k];

                g_0_y_xyyyyy_yyyy[k] = -g_0_y_yyyyy_yyyy[k] * ab_x + g_0_y_yyyyy_xyyyy[k];

                g_0_y_xyyyyy_yyyz[k] = -g_0_y_yyyyy_yyyz[k] * ab_x + g_0_y_yyyyy_xyyyz[k];

                g_0_y_xyyyyy_yyzz[k] = -g_0_y_yyyyy_yyzz[k] * ab_x + g_0_y_yyyyy_xyyzz[k];

                g_0_y_xyyyyy_yzzz[k] = -g_0_y_yyyyy_yzzz[k] * ab_x + g_0_y_yyyyy_xyzzz[k];

                g_0_y_xyyyyy_zzzz[k] = -g_0_y_yyyyy_zzzz[k] * ab_x + g_0_y_yyyyy_xzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyz_xxxx = cbuffer.data(ig_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxxy = cbuffer.data(ig_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxxz = cbuffer.data(ig_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxyy = cbuffer.data(ig_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxyz = cbuffer.data(ig_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxzz = cbuffer.data(ig_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyyy = cbuffer.data(ig_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyyz = cbuffer.data(ig_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyzz = cbuffer.data(ig_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xzzz = cbuffer.data(ig_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyyy = cbuffer.data(ig_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyyz = cbuffer.data(ig_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyzz = cbuffer.data(ig_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yzzz = cbuffer.data(ig_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_y_xyyyyz_zzzz = cbuffer.data(ig_geom_01_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyz_xxxx, g_0_y_xyyyyz_xxxy, g_0_y_xyyyyz_xxxz, g_0_y_xyyyyz_xxyy, g_0_y_xyyyyz_xxyz, g_0_y_xyyyyz_xxzz, g_0_y_xyyyyz_xyyy, g_0_y_xyyyyz_xyyz, g_0_y_xyyyyz_xyzz, g_0_y_xyyyyz_xzzz, g_0_y_xyyyyz_yyyy, g_0_y_xyyyyz_yyyz, g_0_y_xyyyyz_yyzz, g_0_y_xyyyyz_yzzz, g_0_y_xyyyyz_zzzz, g_0_y_yyyyz_xxxx, g_0_y_yyyyz_xxxxx, g_0_y_yyyyz_xxxxy, g_0_y_yyyyz_xxxxz, g_0_y_yyyyz_xxxy, g_0_y_yyyyz_xxxyy, g_0_y_yyyyz_xxxyz, g_0_y_yyyyz_xxxz, g_0_y_yyyyz_xxxzz, g_0_y_yyyyz_xxyy, g_0_y_yyyyz_xxyyy, g_0_y_yyyyz_xxyyz, g_0_y_yyyyz_xxyz, g_0_y_yyyyz_xxyzz, g_0_y_yyyyz_xxzz, g_0_y_yyyyz_xxzzz, g_0_y_yyyyz_xyyy, g_0_y_yyyyz_xyyyy, g_0_y_yyyyz_xyyyz, g_0_y_yyyyz_xyyz, g_0_y_yyyyz_xyyzz, g_0_y_yyyyz_xyzz, g_0_y_yyyyz_xyzzz, g_0_y_yyyyz_xzzz, g_0_y_yyyyz_xzzzz, g_0_y_yyyyz_yyyy, g_0_y_yyyyz_yyyz, g_0_y_yyyyz_yyzz, g_0_y_yyyyz_yzzz, g_0_y_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyz_xxxx[k] = -g_0_y_yyyyz_xxxx[k] * ab_x + g_0_y_yyyyz_xxxxx[k];

                g_0_y_xyyyyz_xxxy[k] = -g_0_y_yyyyz_xxxy[k] * ab_x + g_0_y_yyyyz_xxxxy[k];

                g_0_y_xyyyyz_xxxz[k] = -g_0_y_yyyyz_xxxz[k] * ab_x + g_0_y_yyyyz_xxxxz[k];

                g_0_y_xyyyyz_xxyy[k] = -g_0_y_yyyyz_xxyy[k] * ab_x + g_0_y_yyyyz_xxxyy[k];

                g_0_y_xyyyyz_xxyz[k] = -g_0_y_yyyyz_xxyz[k] * ab_x + g_0_y_yyyyz_xxxyz[k];

                g_0_y_xyyyyz_xxzz[k] = -g_0_y_yyyyz_xxzz[k] * ab_x + g_0_y_yyyyz_xxxzz[k];

                g_0_y_xyyyyz_xyyy[k] = -g_0_y_yyyyz_xyyy[k] * ab_x + g_0_y_yyyyz_xxyyy[k];

                g_0_y_xyyyyz_xyyz[k] = -g_0_y_yyyyz_xyyz[k] * ab_x + g_0_y_yyyyz_xxyyz[k];

                g_0_y_xyyyyz_xyzz[k] = -g_0_y_yyyyz_xyzz[k] * ab_x + g_0_y_yyyyz_xxyzz[k];

                g_0_y_xyyyyz_xzzz[k] = -g_0_y_yyyyz_xzzz[k] * ab_x + g_0_y_yyyyz_xxzzz[k];

                g_0_y_xyyyyz_yyyy[k] = -g_0_y_yyyyz_yyyy[k] * ab_x + g_0_y_yyyyz_xyyyy[k];

                g_0_y_xyyyyz_yyyz[k] = -g_0_y_yyyyz_yyyz[k] * ab_x + g_0_y_yyyyz_xyyyz[k];

                g_0_y_xyyyyz_yyzz[k] = -g_0_y_yyyyz_yyzz[k] * ab_x + g_0_y_yyyyz_xyyzz[k];

                g_0_y_xyyyyz_yzzz[k] = -g_0_y_yyyyz_yzzz[k] * ab_x + g_0_y_yyyyz_xyzzz[k];

                g_0_y_xyyyyz_zzzz[k] = -g_0_y_yyyyz_zzzz[k] * ab_x + g_0_y_yyyyz_xzzzz[k];
            }

            /// Set up 675-690 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyzz_xxxx = cbuffer.data(ig_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxxy = cbuffer.data(ig_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxxz = cbuffer.data(ig_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxyy = cbuffer.data(ig_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxyz = cbuffer.data(ig_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxzz = cbuffer.data(ig_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyyy = cbuffer.data(ig_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyyz = cbuffer.data(ig_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyzz = cbuffer.data(ig_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xzzz = cbuffer.data(ig_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyyy = cbuffer.data(ig_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyyz = cbuffer.data(ig_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyzz = cbuffer.data(ig_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yzzz = cbuffer.data(ig_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_y_xyyyzz_zzzz = cbuffer.data(ig_geom_01_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyzz_xxxx, g_0_y_xyyyzz_xxxy, g_0_y_xyyyzz_xxxz, g_0_y_xyyyzz_xxyy, g_0_y_xyyyzz_xxyz, g_0_y_xyyyzz_xxzz, g_0_y_xyyyzz_xyyy, g_0_y_xyyyzz_xyyz, g_0_y_xyyyzz_xyzz, g_0_y_xyyyzz_xzzz, g_0_y_xyyyzz_yyyy, g_0_y_xyyyzz_yyyz, g_0_y_xyyyzz_yyzz, g_0_y_xyyyzz_yzzz, g_0_y_xyyyzz_zzzz, g_0_y_yyyzz_xxxx, g_0_y_yyyzz_xxxxx, g_0_y_yyyzz_xxxxy, g_0_y_yyyzz_xxxxz, g_0_y_yyyzz_xxxy, g_0_y_yyyzz_xxxyy, g_0_y_yyyzz_xxxyz, g_0_y_yyyzz_xxxz, g_0_y_yyyzz_xxxzz, g_0_y_yyyzz_xxyy, g_0_y_yyyzz_xxyyy, g_0_y_yyyzz_xxyyz, g_0_y_yyyzz_xxyz, g_0_y_yyyzz_xxyzz, g_0_y_yyyzz_xxzz, g_0_y_yyyzz_xxzzz, g_0_y_yyyzz_xyyy, g_0_y_yyyzz_xyyyy, g_0_y_yyyzz_xyyyz, g_0_y_yyyzz_xyyz, g_0_y_yyyzz_xyyzz, g_0_y_yyyzz_xyzz, g_0_y_yyyzz_xyzzz, g_0_y_yyyzz_xzzz, g_0_y_yyyzz_xzzzz, g_0_y_yyyzz_yyyy, g_0_y_yyyzz_yyyz, g_0_y_yyyzz_yyzz, g_0_y_yyyzz_yzzz, g_0_y_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyzz_xxxx[k] = -g_0_y_yyyzz_xxxx[k] * ab_x + g_0_y_yyyzz_xxxxx[k];

                g_0_y_xyyyzz_xxxy[k] = -g_0_y_yyyzz_xxxy[k] * ab_x + g_0_y_yyyzz_xxxxy[k];

                g_0_y_xyyyzz_xxxz[k] = -g_0_y_yyyzz_xxxz[k] * ab_x + g_0_y_yyyzz_xxxxz[k];

                g_0_y_xyyyzz_xxyy[k] = -g_0_y_yyyzz_xxyy[k] * ab_x + g_0_y_yyyzz_xxxyy[k];

                g_0_y_xyyyzz_xxyz[k] = -g_0_y_yyyzz_xxyz[k] * ab_x + g_0_y_yyyzz_xxxyz[k];

                g_0_y_xyyyzz_xxzz[k] = -g_0_y_yyyzz_xxzz[k] * ab_x + g_0_y_yyyzz_xxxzz[k];

                g_0_y_xyyyzz_xyyy[k] = -g_0_y_yyyzz_xyyy[k] * ab_x + g_0_y_yyyzz_xxyyy[k];

                g_0_y_xyyyzz_xyyz[k] = -g_0_y_yyyzz_xyyz[k] * ab_x + g_0_y_yyyzz_xxyyz[k];

                g_0_y_xyyyzz_xyzz[k] = -g_0_y_yyyzz_xyzz[k] * ab_x + g_0_y_yyyzz_xxyzz[k];

                g_0_y_xyyyzz_xzzz[k] = -g_0_y_yyyzz_xzzz[k] * ab_x + g_0_y_yyyzz_xxzzz[k];

                g_0_y_xyyyzz_yyyy[k] = -g_0_y_yyyzz_yyyy[k] * ab_x + g_0_y_yyyzz_xyyyy[k];

                g_0_y_xyyyzz_yyyz[k] = -g_0_y_yyyzz_yyyz[k] * ab_x + g_0_y_yyyzz_xyyyz[k];

                g_0_y_xyyyzz_yyzz[k] = -g_0_y_yyyzz_yyzz[k] * ab_x + g_0_y_yyyzz_xyyzz[k];

                g_0_y_xyyyzz_yzzz[k] = -g_0_y_yyyzz_yzzz[k] * ab_x + g_0_y_yyyzz_xyzzz[k];

                g_0_y_xyyyzz_zzzz[k] = -g_0_y_yyyzz_zzzz[k] * ab_x + g_0_y_yyyzz_xzzzz[k];
            }

            /// Set up 690-705 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzzz_xxxx = cbuffer.data(ig_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxxy = cbuffer.data(ig_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxxz = cbuffer.data(ig_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxyy = cbuffer.data(ig_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxyz = cbuffer.data(ig_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxzz = cbuffer.data(ig_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyyy = cbuffer.data(ig_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyyz = cbuffer.data(ig_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyzz = cbuffer.data(ig_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xzzz = cbuffer.data(ig_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyyy = cbuffer.data(ig_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyyz = cbuffer.data(ig_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyzz = cbuffer.data(ig_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yzzz = cbuffer.data(ig_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_y_xyyzzz_zzzz = cbuffer.data(ig_geom_01_off + 704 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzzz_xxxx, g_0_y_xyyzzz_xxxy, g_0_y_xyyzzz_xxxz, g_0_y_xyyzzz_xxyy, g_0_y_xyyzzz_xxyz, g_0_y_xyyzzz_xxzz, g_0_y_xyyzzz_xyyy, g_0_y_xyyzzz_xyyz, g_0_y_xyyzzz_xyzz, g_0_y_xyyzzz_xzzz, g_0_y_xyyzzz_yyyy, g_0_y_xyyzzz_yyyz, g_0_y_xyyzzz_yyzz, g_0_y_xyyzzz_yzzz, g_0_y_xyyzzz_zzzz, g_0_y_yyzzz_xxxx, g_0_y_yyzzz_xxxxx, g_0_y_yyzzz_xxxxy, g_0_y_yyzzz_xxxxz, g_0_y_yyzzz_xxxy, g_0_y_yyzzz_xxxyy, g_0_y_yyzzz_xxxyz, g_0_y_yyzzz_xxxz, g_0_y_yyzzz_xxxzz, g_0_y_yyzzz_xxyy, g_0_y_yyzzz_xxyyy, g_0_y_yyzzz_xxyyz, g_0_y_yyzzz_xxyz, g_0_y_yyzzz_xxyzz, g_0_y_yyzzz_xxzz, g_0_y_yyzzz_xxzzz, g_0_y_yyzzz_xyyy, g_0_y_yyzzz_xyyyy, g_0_y_yyzzz_xyyyz, g_0_y_yyzzz_xyyz, g_0_y_yyzzz_xyyzz, g_0_y_yyzzz_xyzz, g_0_y_yyzzz_xyzzz, g_0_y_yyzzz_xzzz, g_0_y_yyzzz_xzzzz, g_0_y_yyzzz_yyyy, g_0_y_yyzzz_yyyz, g_0_y_yyzzz_yyzz, g_0_y_yyzzz_yzzz, g_0_y_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzzz_xxxx[k] = -g_0_y_yyzzz_xxxx[k] * ab_x + g_0_y_yyzzz_xxxxx[k];

                g_0_y_xyyzzz_xxxy[k] = -g_0_y_yyzzz_xxxy[k] * ab_x + g_0_y_yyzzz_xxxxy[k];

                g_0_y_xyyzzz_xxxz[k] = -g_0_y_yyzzz_xxxz[k] * ab_x + g_0_y_yyzzz_xxxxz[k];

                g_0_y_xyyzzz_xxyy[k] = -g_0_y_yyzzz_xxyy[k] * ab_x + g_0_y_yyzzz_xxxyy[k];

                g_0_y_xyyzzz_xxyz[k] = -g_0_y_yyzzz_xxyz[k] * ab_x + g_0_y_yyzzz_xxxyz[k];

                g_0_y_xyyzzz_xxzz[k] = -g_0_y_yyzzz_xxzz[k] * ab_x + g_0_y_yyzzz_xxxzz[k];

                g_0_y_xyyzzz_xyyy[k] = -g_0_y_yyzzz_xyyy[k] * ab_x + g_0_y_yyzzz_xxyyy[k];

                g_0_y_xyyzzz_xyyz[k] = -g_0_y_yyzzz_xyyz[k] * ab_x + g_0_y_yyzzz_xxyyz[k];

                g_0_y_xyyzzz_xyzz[k] = -g_0_y_yyzzz_xyzz[k] * ab_x + g_0_y_yyzzz_xxyzz[k];

                g_0_y_xyyzzz_xzzz[k] = -g_0_y_yyzzz_xzzz[k] * ab_x + g_0_y_yyzzz_xxzzz[k];

                g_0_y_xyyzzz_yyyy[k] = -g_0_y_yyzzz_yyyy[k] * ab_x + g_0_y_yyzzz_xyyyy[k];

                g_0_y_xyyzzz_yyyz[k] = -g_0_y_yyzzz_yyyz[k] * ab_x + g_0_y_yyzzz_xyyyz[k];

                g_0_y_xyyzzz_yyzz[k] = -g_0_y_yyzzz_yyzz[k] * ab_x + g_0_y_yyzzz_xyyzz[k];

                g_0_y_xyyzzz_yzzz[k] = -g_0_y_yyzzz_yzzz[k] * ab_x + g_0_y_yyzzz_xyzzz[k];

                g_0_y_xyyzzz_zzzz[k] = -g_0_y_yyzzz_zzzz[k] * ab_x + g_0_y_yyzzz_xzzzz[k];
            }

            /// Set up 705-720 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzzz_xxxx = cbuffer.data(ig_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxxy = cbuffer.data(ig_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxxz = cbuffer.data(ig_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxyy = cbuffer.data(ig_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxyz = cbuffer.data(ig_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxzz = cbuffer.data(ig_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyyy = cbuffer.data(ig_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyyz = cbuffer.data(ig_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyzz = cbuffer.data(ig_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xzzz = cbuffer.data(ig_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyyy = cbuffer.data(ig_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyyz = cbuffer.data(ig_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyzz = cbuffer.data(ig_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yzzz = cbuffer.data(ig_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_y_xyzzzz_zzzz = cbuffer.data(ig_geom_01_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzzz_xxxx, g_0_y_xyzzzz_xxxy, g_0_y_xyzzzz_xxxz, g_0_y_xyzzzz_xxyy, g_0_y_xyzzzz_xxyz, g_0_y_xyzzzz_xxzz, g_0_y_xyzzzz_xyyy, g_0_y_xyzzzz_xyyz, g_0_y_xyzzzz_xyzz, g_0_y_xyzzzz_xzzz, g_0_y_xyzzzz_yyyy, g_0_y_xyzzzz_yyyz, g_0_y_xyzzzz_yyzz, g_0_y_xyzzzz_yzzz, g_0_y_xyzzzz_zzzz, g_0_y_yzzzz_xxxx, g_0_y_yzzzz_xxxxx, g_0_y_yzzzz_xxxxy, g_0_y_yzzzz_xxxxz, g_0_y_yzzzz_xxxy, g_0_y_yzzzz_xxxyy, g_0_y_yzzzz_xxxyz, g_0_y_yzzzz_xxxz, g_0_y_yzzzz_xxxzz, g_0_y_yzzzz_xxyy, g_0_y_yzzzz_xxyyy, g_0_y_yzzzz_xxyyz, g_0_y_yzzzz_xxyz, g_0_y_yzzzz_xxyzz, g_0_y_yzzzz_xxzz, g_0_y_yzzzz_xxzzz, g_0_y_yzzzz_xyyy, g_0_y_yzzzz_xyyyy, g_0_y_yzzzz_xyyyz, g_0_y_yzzzz_xyyz, g_0_y_yzzzz_xyyzz, g_0_y_yzzzz_xyzz, g_0_y_yzzzz_xyzzz, g_0_y_yzzzz_xzzz, g_0_y_yzzzz_xzzzz, g_0_y_yzzzz_yyyy, g_0_y_yzzzz_yyyz, g_0_y_yzzzz_yyzz, g_0_y_yzzzz_yzzz, g_0_y_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzzz_xxxx[k] = -g_0_y_yzzzz_xxxx[k] * ab_x + g_0_y_yzzzz_xxxxx[k];

                g_0_y_xyzzzz_xxxy[k] = -g_0_y_yzzzz_xxxy[k] * ab_x + g_0_y_yzzzz_xxxxy[k];

                g_0_y_xyzzzz_xxxz[k] = -g_0_y_yzzzz_xxxz[k] * ab_x + g_0_y_yzzzz_xxxxz[k];

                g_0_y_xyzzzz_xxyy[k] = -g_0_y_yzzzz_xxyy[k] * ab_x + g_0_y_yzzzz_xxxyy[k];

                g_0_y_xyzzzz_xxyz[k] = -g_0_y_yzzzz_xxyz[k] * ab_x + g_0_y_yzzzz_xxxyz[k];

                g_0_y_xyzzzz_xxzz[k] = -g_0_y_yzzzz_xxzz[k] * ab_x + g_0_y_yzzzz_xxxzz[k];

                g_0_y_xyzzzz_xyyy[k] = -g_0_y_yzzzz_xyyy[k] * ab_x + g_0_y_yzzzz_xxyyy[k];

                g_0_y_xyzzzz_xyyz[k] = -g_0_y_yzzzz_xyyz[k] * ab_x + g_0_y_yzzzz_xxyyz[k];

                g_0_y_xyzzzz_xyzz[k] = -g_0_y_yzzzz_xyzz[k] * ab_x + g_0_y_yzzzz_xxyzz[k];

                g_0_y_xyzzzz_xzzz[k] = -g_0_y_yzzzz_xzzz[k] * ab_x + g_0_y_yzzzz_xxzzz[k];

                g_0_y_xyzzzz_yyyy[k] = -g_0_y_yzzzz_yyyy[k] * ab_x + g_0_y_yzzzz_xyyyy[k];

                g_0_y_xyzzzz_yyyz[k] = -g_0_y_yzzzz_yyyz[k] * ab_x + g_0_y_yzzzz_xyyyz[k];

                g_0_y_xyzzzz_yyzz[k] = -g_0_y_yzzzz_yyzz[k] * ab_x + g_0_y_yzzzz_xyyzz[k];

                g_0_y_xyzzzz_yzzz[k] = -g_0_y_yzzzz_yzzz[k] * ab_x + g_0_y_yzzzz_xyzzz[k];

                g_0_y_xyzzzz_zzzz[k] = -g_0_y_yzzzz_zzzz[k] * ab_x + g_0_y_yzzzz_xzzzz[k];
            }

            /// Set up 720-735 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzzz_xxxx = cbuffer.data(ig_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxxy = cbuffer.data(ig_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxxz = cbuffer.data(ig_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxyy = cbuffer.data(ig_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxyz = cbuffer.data(ig_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxzz = cbuffer.data(ig_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyyy = cbuffer.data(ig_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyyz = cbuffer.data(ig_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyzz = cbuffer.data(ig_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xzzz = cbuffer.data(ig_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyyy = cbuffer.data(ig_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyyz = cbuffer.data(ig_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyzz = cbuffer.data(ig_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yzzz = cbuffer.data(ig_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_y_xzzzzz_zzzz = cbuffer.data(ig_geom_01_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzzz_xxxx, g_0_y_xzzzzz_xxxy, g_0_y_xzzzzz_xxxz, g_0_y_xzzzzz_xxyy, g_0_y_xzzzzz_xxyz, g_0_y_xzzzzz_xxzz, g_0_y_xzzzzz_xyyy, g_0_y_xzzzzz_xyyz, g_0_y_xzzzzz_xyzz, g_0_y_xzzzzz_xzzz, g_0_y_xzzzzz_yyyy, g_0_y_xzzzzz_yyyz, g_0_y_xzzzzz_yyzz, g_0_y_xzzzzz_yzzz, g_0_y_xzzzzz_zzzz, g_0_y_zzzzz_xxxx, g_0_y_zzzzz_xxxxx, g_0_y_zzzzz_xxxxy, g_0_y_zzzzz_xxxxz, g_0_y_zzzzz_xxxy, g_0_y_zzzzz_xxxyy, g_0_y_zzzzz_xxxyz, g_0_y_zzzzz_xxxz, g_0_y_zzzzz_xxxzz, g_0_y_zzzzz_xxyy, g_0_y_zzzzz_xxyyy, g_0_y_zzzzz_xxyyz, g_0_y_zzzzz_xxyz, g_0_y_zzzzz_xxyzz, g_0_y_zzzzz_xxzz, g_0_y_zzzzz_xxzzz, g_0_y_zzzzz_xyyy, g_0_y_zzzzz_xyyyy, g_0_y_zzzzz_xyyyz, g_0_y_zzzzz_xyyz, g_0_y_zzzzz_xyyzz, g_0_y_zzzzz_xyzz, g_0_y_zzzzz_xyzzz, g_0_y_zzzzz_xzzz, g_0_y_zzzzz_xzzzz, g_0_y_zzzzz_yyyy, g_0_y_zzzzz_yyyz, g_0_y_zzzzz_yyzz, g_0_y_zzzzz_yzzz, g_0_y_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzzz_xxxx[k] = -g_0_y_zzzzz_xxxx[k] * ab_x + g_0_y_zzzzz_xxxxx[k];

                g_0_y_xzzzzz_xxxy[k] = -g_0_y_zzzzz_xxxy[k] * ab_x + g_0_y_zzzzz_xxxxy[k];

                g_0_y_xzzzzz_xxxz[k] = -g_0_y_zzzzz_xxxz[k] * ab_x + g_0_y_zzzzz_xxxxz[k];

                g_0_y_xzzzzz_xxyy[k] = -g_0_y_zzzzz_xxyy[k] * ab_x + g_0_y_zzzzz_xxxyy[k];

                g_0_y_xzzzzz_xxyz[k] = -g_0_y_zzzzz_xxyz[k] * ab_x + g_0_y_zzzzz_xxxyz[k];

                g_0_y_xzzzzz_xxzz[k] = -g_0_y_zzzzz_xxzz[k] * ab_x + g_0_y_zzzzz_xxxzz[k];

                g_0_y_xzzzzz_xyyy[k] = -g_0_y_zzzzz_xyyy[k] * ab_x + g_0_y_zzzzz_xxyyy[k];

                g_0_y_xzzzzz_xyyz[k] = -g_0_y_zzzzz_xyyz[k] * ab_x + g_0_y_zzzzz_xxyyz[k];

                g_0_y_xzzzzz_xyzz[k] = -g_0_y_zzzzz_xyzz[k] * ab_x + g_0_y_zzzzz_xxyzz[k];

                g_0_y_xzzzzz_xzzz[k] = -g_0_y_zzzzz_xzzz[k] * ab_x + g_0_y_zzzzz_xxzzz[k];

                g_0_y_xzzzzz_yyyy[k] = -g_0_y_zzzzz_yyyy[k] * ab_x + g_0_y_zzzzz_xyyyy[k];

                g_0_y_xzzzzz_yyyz[k] = -g_0_y_zzzzz_yyyz[k] * ab_x + g_0_y_zzzzz_xyyyz[k];

                g_0_y_xzzzzz_yyzz[k] = -g_0_y_zzzzz_yyzz[k] * ab_x + g_0_y_zzzzz_xyyzz[k];

                g_0_y_xzzzzz_yzzz[k] = -g_0_y_zzzzz_yzzz[k] * ab_x + g_0_y_zzzzz_xyzzz[k];

                g_0_y_xzzzzz_zzzz[k] = -g_0_y_zzzzz_zzzz[k] * ab_x + g_0_y_zzzzz_xzzzz[k];
            }

            /// Set up 735-750 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyy_xxxx = cbuffer.data(ig_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxxy = cbuffer.data(ig_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxxz = cbuffer.data(ig_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxyy = cbuffer.data(ig_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxyz = cbuffer.data(ig_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxzz = cbuffer.data(ig_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyyy = cbuffer.data(ig_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyyz = cbuffer.data(ig_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyzz = cbuffer.data(ig_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xzzz = cbuffer.data(ig_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyyy = cbuffer.data(ig_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyyz = cbuffer.data(ig_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyzz = cbuffer.data(ig_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yzzz = cbuffer.data(ig_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_y_yyyyyy_zzzz = cbuffer.data(ig_geom_01_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyy_xxxx, g_0_y_yyyyy_xxxxy, g_0_y_yyyyy_xxxy, g_0_y_yyyyy_xxxyy, g_0_y_yyyyy_xxxyz, g_0_y_yyyyy_xxxz, g_0_y_yyyyy_xxyy, g_0_y_yyyyy_xxyyy, g_0_y_yyyyy_xxyyz, g_0_y_yyyyy_xxyz, g_0_y_yyyyy_xxyzz, g_0_y_yyyyy_xxzz, g_0_y_yyyyy_xyyy, g_0_y_yyyyy_xyyyy, g_0_y_yyyyy_xyyyz, g_0_y_yyyyy_xyyz, g_0_y_yyyyy_xyyzz, g_0_y_yyyyy_xyzz, g_0_y_yyyyy_xyzzz, g_0_y_yyyyy_xzzz, g_0_y_yyyyy_yyyy, g_0_y_yyyyy_yyyyy, g_0_y_yyyyy_yyyyz, g_0_y_yyyyy_yyyz, g_0_y_yyyyy_yyyzz, g_0_y_yyyyy_yyzz, g_0_y_yyyyy_yyzzz, g_0_y_yyyyy_yzzz, g_0_y_yyyyy_yzzzz, g_0_y_yyyyy_zzzz, g_0_y_yyyyyy_xxxx, g_0_y_yyyyyy_xxxy, g_0_y_yyyyyy_xxxz, g_0_y_yyyyyy_xxyy, g_0_y_yyyyyy_xxyz, g_0_y_yyyyyy_xxzz, g_0_y_yyyyyy_xyyy, g_0_y_yyyyyy_xyyz, g_0_y_yyyyyy_xyzz, g_0_y_yyyyyy_xzzz, g_0_y_yyyyyy_yyyy, g_0_y_yyyyyy_yyyz, g_0_y_yyyyyy_yyzz, g_0_y_yyyyyy_yzzz, g_0_y_yyyyyy_zzzz, g_yyyyy_xxxx, g_yyyyy_xxxy, g_yyyyy_xxxz, g_yyyyy_xxyy, g_yyyyy_xxyz, g_yyyyy_xxzz, g_yyyyy_xyyy, g_yyyyy_xyyz, g_yyyyy_xyzz, g_yyyyy_xzzz, g_yyyyy_yyyy, g_yyyyy_yyyz, g_yyyyy_yyzz, g_yyyyy_yzzz, g_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyy_xxxx[k] = g_yyyyy_xxxx[k] - g_0_y_yyyyy_xxxx[k] * ab_y + g_0_y_yyyyy_xxxxy[k];

                g_0_y_yyyyyy_xxxy[k] = g_yyyyy_xxxy[k] - g_0_y_yyyyy_xxxy[k] * ab_y + g_0_y_yyyyy_xxxyy[k];

                g_0_y_yyyyyy_xxxz[k] = g_yyyyy_xxxz[k] - g_0_y_yyyyy_xxxz[k] * ab_y + g_0_y_yyyyy_xxxyz[k];

                g_0_y_yyyyyy_xxyy[k] = g_yyyyy_xxyy[k] - g_0_y_yyyyy_xxyy[k] * ab_y + g_0_y_yyyyy_xxyyy[k];

                g_0_y_yyyyyy_xxyz[k] = g_yyyyy_xxyz[k] - g_0_y_yyyyy_xxyz[k] * ab_y + g_0_y_yyyyy_xxyyz[k];

                g_0_y_yyyyyy_xxzz[k] = g_yyyyy_xxzz[k] - g_0_y_yyyyy_xxzz[k] * ab_y + g_0_y_yyyyy_xxyzz[k];

                g_0_y_yyyyyy_xyyy[k] = g_yyyyy_xyyy[k] - g_0_y_yyyyy_xyyy[k] * ab_y + g_0_y_yyyyy_xyyyy[k];

                g_0_y_yyyyyy_xyyz[k] = g_yyyyy_xyyz[k] - g_0_y_yyyyy_xyyz[k] * ab_y + g_0_y_yyyyy_xyyyz[k];

                g_0_y_yyyyyy_xyzz[k] = g_yyyyy_xyzz[k] - g_0_y_yyyyy_xyzz[k] * ab_y + g_0_y_yyyyy_xyyzz[k];

                g_0_y_yyyyyy_xzzz[k] = g_yyyyy_xzzz[k] - g_0_y_yyyyy_xzzz[k] * ab_y + g_0_y_yyyyy_xyzzz[k];

                g_0_y_yyyyyy_yyyy[k] = g_yyyyy_yyyy[k] - g_0_y_yyyyy_yyyy[k] * ab_y + g_0_y_yyyyy_yyyyy[k];

                g_0_y_yyyyyy_yyyz[k] = g_yyyyy_yyyz[k] - g_0_y_yyyyy_yyyz[k] * ab_y + g_0_y_yyyyy_yyyyz[k];

                g_0_y_yyyyyy_yyzz[k] = g_yyyyy_yyzz[k] - g_0_y_yyyyy_yyzz[k] * ab_y + g_0_y_yyyyy_yyyzz[k];

                g_0_y_yyyyyy_yzzz[k] = g_yyyyy_yzzz[k] - g_0_y_yyyyy_yzzz[k] * ab_y + g_0_y_yyyyy_yyzzz[k];

                g_0_y_yyyyyy_zzzz[k] = g_yyyyy_zzzz[k] - g_0_y_yyyyy_zzzz[k] * ab_y + g_0_y_yyyyy_yzzzz[k];
            }

            /// Set up 750-765 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyz_xxxx = cbuffer.data(ig_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxxy = cbuffer.data(ig_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxxz = cbuffer.data(ig_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxyy = cbuffer.data(ig_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxyz = cbuffer.data(ig_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxzz = cbuffer.data(ig_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyyy = cbuffer.data(ig_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyyz = cbuffer.data(ig_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyzz = cbuffer.data(ig_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xzzz = cbuffer.data(ig_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyyy = cbuffer.data(ig_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyyz = cbuffer.data(ig_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyzz = cbuffer.data(ig_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yzzz = cbuffer.data(ig_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_y_yyyyyz_zzzz = cbuffer.data(ig_geom_01_off + 764 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyy_xxxx, g_0_y_yyyyy_xxxxz, g_0_y_yyyyy_xxxy, g_0_y_yyyyy_xxxyz, g_0_y_yyyyy_xxxz, g_0_y_yyyyy_xxxzz, g_0_y_yyyyy_xxyy, g_0_y_yyyyy_xxyyz, g_0_y_yyyyy_xxyz, g_0_y_yyyyy_xxyzz, g_0_y_yyyyy_xxzz, g_0_y_yyyyy_xxzzz, g_0_y_yyyyy_xyyy, g_0_y_yyyyy_xyyyz, g_0_y_yyyyy_xyyz, g_0_y_yyyyy_xyyzz, g_0_y_yyyyy_xyzz, g_0_y_yyyyy_xyzzz, g_0_y_yyyyy_xzzz, g_0_y_yyyyy_xzzzz, g_0_y_yyyyy_yyyy, g_0_y_yyyyy_yyyyz, g_0_y_yyyyy_yyyz, g_0_y_yyyyy_yyyzz, g_0_y_yyyyy_yyzz, g_0_y_yyyyy_yyzzz, g_0_y_yyyyy_yzzz, g_0_y_yyyyy_yzzzz, g_0_y_yyyyy_zzzz, g_0_y_yyyyy_zzzzz, g_0_y_yyyyyz_xxxx, g_0_y_yyyyyz_xxxy, g_0_y_yyyyyz_xxxz, g_0_y_yyyyyz_xxyy, g_0_y_yyyyyz_xxyz, g_0_y_yyyyyz_xxzz, g_0_y_yyyyyz_xyyy, g_0_y_yyyyyz_xyyz, g_0_y_yyyyyz_xyzz, g_0_y_yyyyyz_xzzz, g_0_y_yyyyyz_yyyy, g_0_y_yyyyyz_yyyz, g_0_y_yyyyyz_yyzz, g_0_y_yyyyyz_yzzz, g_0_y_yyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyz_xxxx[k] = -g_0_y_yyyyy_xxxx[k] * ab_z + g_0_y_yyyyy_xxxxz[k];

                g_0_y_yyyyyz_xxxy[k] = -g_0_y_yyyyy_xxxy[k] * ab_z + g_0_y_yyyyy_xxxyz[k];

                g_0_y_yyyyyz_xxxz[k] = -g_0_y_yyyyy_xxxz[k] * ab_z + g_0_y_yyyyy_xxxzz[k];

                g_0_y_yyyyyz_xxyy[k] = -g_0_y_yyyyy_xxyy[k] * ab_z + g_0_y_yyyyy_xxyyz[k];

                g_0_y_yyyyyz_xxyz[k] = -g_0_y_yyyyy_xxyz[k] * ab_z + g_0_y_yyyyy_xxyzz[k];

                g_0_y_yyyyyz_xxzz[k] = -g_0_y_yyyyy_xxzz[k] * ab_z + g_0_y_yyyyy_xxzzz[k];

                g_0_y_yyyyyz_xyyy[k] = -g_0_y_yyyyy_xyyy[k] * ab_z + g_0_y_yyyyy_xyyyz[k];

                g_0_y_yyyyyz_xyyz[k] = -g_0_y_yyyyy_xyyz[k] * ab_z + g_0_y_yyyyy_xyyzz[k];

                g_0_y_yyyyyz_xyzz[k] = -g_0_y_yyyyy_xyzz[k] * ab_z + g_0_y_yyyyy_xyzzz[k];

                g_0_y_yyyyyz_xzzz[k] = -g_0_y_yyyyy_xzzz[k] * ab_z + g_0_y_yyyyy_xzzzz[k];

                g_0_y_yyyyyz_yyyy[k] = -g_0_y_yyyyy_yyyy[k] * ab_z + g_0_y_yyyyy_yyyyz[k];

                g_0_y_yyyyyz_yyyz[k] = -g_0_y_yyyyy_yyyz[k] * ab_z + g_0_y_yyyyy_yyyzz[k];

                g_0_y_yyyyyz_yyzz[k] = -g_0_y_yyyyy_yyzz[k] * ab_z + g_0_y_yyyyy_yyzzz[k];

                g_0_y_yyyyyz_yzzz[k] = -g_0_y_yyyyy_yzzz[k] * ab_z + g_0_y_yyyyy_yzzzz[k];

                g_0_y_yyyyyz_zzzz[k] = -g_0_y_yyyyy_zzzz[k] * ab_z + g_0_y_yyyyy_zzzzz[k];
            }

            /// Set up 765-780 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyzz_xxxx = cbuffer.data(ig_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxxy = cbuffer.data(ig_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxxz = cbuffer.data(ig_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxyy = cbuffer.data(ig_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxyz = cbuffer.data(ig_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxzz = cbuffer.data(ig_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyyy = cbuffer.data(ig_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyyz = cbuffer.data(ig_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyzz = cbuffer.data(ig_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xzzz = cbuffer.data(ig_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyyy = cbuffer.data(ig_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyyz = cbuffer.data(ig_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyzz = cbuffer.data(ig_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yzzz = cbuffer.data(ig_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_y_yyyyzz_zzzz = cbuffer.data(ig_geom_01_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyz_xxxx, g_0_y_yyyyz_xxxxz, g_0_y_yyyyz_xxxy, g_0_y_yyyyz_xxxyz, g_0_y_yyyyz_xxxz, g_0_y_yyyyz_xxxzz, g_0_y_yyyyz_xxyy, g_0_y_yyyyz_xxyyz, g_0_y_yyyyz_xxyz, g_0_y_yyyyz_xxyzz, g_0_y_yyyyz_xxzz, g_0_y_yyyyz_xxzzz, g_0_y_yyyyz_xyyy, g_0_y_yyyyz_xyyyz, g_0_y_yyyyz_xyyz, g_0_y_yyyyz_xyyzz, g_0_y_yyyyz_xyzz, g_0_y_yyyyz_xyzzz, g_0_y_yyyyz_xzzz, g_0_y_yyyyz_xzzzz, g_0_y_yyyyz_yyyy, g_0_y_yyyyz_yyyyz, g_0_y_yyyyz_yyyz, g_0_y_yyyyz_yyyzz, g_0_y_yyyyz_yyzz, g_0_y_yyyyz_yyzzz, g_0_y_yyyyz_yzzz, g_0_y_yyyyz_yzzzz, g_0_y_yyyyz_zzzz, g_0_y_yyyyz_zzzzz, g_0_y_yyyyzz_xxxx, g_0_y_yyyyzz_xxxy, g_0_y_yyyyzz_xxxz, g_0_y_yyyyzz_xxyy, g_0_y_yyyyzz_xxyz, g_0_y_yyyyzz_xxzz, g_0_y_yyyyzz_xyyy, g_0_y_yyyyzz_xyyz, g_0_y_yyyyzz_xyzz, g_0_y_yyyyzz_xzzz, g_0_y_yyyyzz_yyyy, g_0_y_yyyyzz_yyyz, g_0_y_yyyyzz_yyzz, g_0_y_yyyyzz_yzzz, g_0_y_yyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyzz_xxxx[k] = -g_0_y_yyyyz_xxxx[k] * ab_z + g_0_y_yyyyz_xxxxz[k];

                g_0_y_yyyyzz_xxxy[k] = -g_0_y_yyyyz_xxxy[k] * ab_z + g_0_y_yyyyz_xxxyz[k];

                g_0_y_yyyyzz_xxxz[k] = -g_0_y_yyyyz_xxxz[k] * ab_z + g_0_y_yyyyz_xxxzz[k];

                g_0_y_yyyyzz_xxyy[k] = -g_0_y_yyyyz_xxyy[k] * ab_z + g_0_y_yyyyz_xxyyz[k];

                g_0_y_yyyyzz_xxyz[k] = -g_0_y_yyyyz_xxyz[k] * ab_z + g_0_y_yyyyz_xxyzz[k];

                g_0_y_yyyyzz_xxzz[k] = -g_0_y_yyyyz_xxzz[k] * ab_z + g_0_y_yyyyz_xxzzz[k];

                g_0_y_yyyyzz_xyyy[k] = -g_0_y_yyyyz_xyyy[k] * ab_z + g_0_y_yyyyz_xyyyz[k];

                g_0_y_yyyyzz_xyyz[k] = -g_0_y_yyyyz_xyyz[k] * ab_z + g_0_y_yyyyz_xyyzz[k];

                g_0_y_yyyyzz_xyzz[k] = -g_0_y_yyyyz_xyzz[k] * ab_z + g_0_y_yyyyz_xyzzz[k];

                g_0_y_yyyyzz_xzzz[k] = -g_0_y_yyyyz_xzzz[k] * ab_z + g_0_y_yyyyz_xzzzz[k];

                g_0_y_yyyyzz_yyyy[k] = -g_0_y_yyyyz_yyyy[k] * ab_z + g_0_y_yyyyz_yyyyz[k];

                g_0_y_yyyyzz_yyyz[k] = -g_0_y_yyyyz_yyyz[k] * ab_z + g_0_y_yyyyz_yyyzz[k];

                g_0_y_yyyyzz_yyzz[k] = -g_0_y_yyyyz_yyzz[k] * ab_z + g_0_y_yyyyz_yyzzz[k];

                g_0_y_yyyyzz_yzzz[k] = -g_0_y_yyyyz_yzzz[k] * ab_z + g_0_y_yyyyz_yzzzz[k];

                g_0_y_yyyyzz_zzzz[k] = -g_0_y_yyyyz_zzzz[k] * ab_z + g_0_y_yyyyz_zzzzz[k];
            }

            /// Set up 780-795 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzzz_xxxx = cbuffer.data(ig_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxxy = cbuffer.data(ig_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxxz = cbuffer.data(ig_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxyy = cbuffer.data(ig_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxyz = cbuffer.data(ig_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxzz = cbuffer.data(ig_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyyy = cbuffer.data(ig_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyyz = cbuffer.data(ig_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyzz = cbuffer.data(ig_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xzzz = cbuffer.data(ig_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyyy = cbuffer.data(ig_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyyz = cbuffer.data(ig_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyzz = cbuffer.data(ig_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yzzz = cbuffer.data(ig_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_y_yyyzzz_zzzz = cbuffer.data(ig_geom_01_off + 794 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyzz_xxxx, g_0_y_yyyzz_xxxxz, g_0_y_yyyzz_xxxy, g_0_y_yyyzz_xxxyz, g_0_y_yyyzz_xxxz, g_0_y_yyyzz_xxxzz, g_0_y_yyyzz_xxyy, g_0_y_yyyzz_xxyyz, g_0_y_yyyzz_xxyz, g_0_y_yyyzz_xxyzz, g_0_y_yyyzz_xxzz, g_0_y_yyyzz_xxzzz, g_0_y_yyyzz_xyyy, g_0_y_yyyzz_xyyyz, g_0_y_yyyzz_xyyz, g_0_y_yyyzz_xyyzz, g_0_y_yyyzz_xyzz, g_0_y_yyyzz_xyzzz, g_0_y_yyyzz_xzzz, g_0_y_yyyzz_xzzzz, g_0_y_yyyzz_yyyy, g_0_y_yyyzz_yyyyz, g_0_y_yyyzz_yyyz, g_0_y_yyyzz_yyyzz, g_0_y_yyyzz_yyzz, g_0_y_yyyzz_yyzzz, g_0_y_yyyzz_yzzz, g_0_y_yyyzz_yzzzz, g_0_y_yyyzz_zzzz, g_0_y_yyyzz_zzzzz, g_0_y_yyyzzz_xxxx, g_0_y_yyyzzz_xxxy, g_0_y_yyyzzz_xxxz, g_0_y_yyyzzz_xxyy, g_0_y_yyyzzz_xxyz, g_0_y_yyyzzz_xxzz, g_0_y_yyyzzz_xyyy, g_0_y_yyyzzz_xyyz, g_0_y_yyyzzz_xyzz, g_0_y_yyyzzz_xzzz, g_0_y_yyyzzz_yyyy, g_0_y_yyyzzz_yyyz, g_0_y_yyyzzz_yyzz, g_0_y_yyyzzz_yzzz, g_0_y_yyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzzz_xxxx[k] = -g_0_y_yyyzz_xxxx[k] * ab_z + g_0_y_yyyzz_xxxxz[k];

                g_0_y_yyyzzz_xxxy[k] = -g_0_y_yyyzz_xxxy[k] * ab_z + g_0_y_yyyzz_xxxyz[k];

                g_0_y_yyyzzz_xxxz[k] = -g_0_y_yyyzz_xxxz[k] * ab_z + g_0_y_yyyzz_xxxzz[k];

                g_0_y_yyyzzz_xxyy[k] = -g_0_y_yyyzz_xxyy[k] * ab_z + g_0_y_yyyzz_xxyyz[k];

                g_0_y_yyyzzz_xxyz[k] = -g_0_y_yyyzz_xxyz[k] * ab_z + g_0_y_yyyzz_xxyzz[k];

                g_0_y_yyyzzz_xxzz[k] = -g_0_y_yyyzz_xxzz[k] * ab_z + g_0_y_yyyzz_xxzzz[k];

                g_0_y_yyyzzz_xyyy[k] = -g_0_y_yyyzz_xyyy[k] * ab_z + g_0_y_yyyzz_xyyyz[k];

                g_0_y_yyyzzz_xyyz[k] = -g_0_y_yyyzz_xyyz[k] * ab_z + g_0_y_yyyzz_xyyzz[k];

                g_0_y_yyyzzz_xyzz[k] = -g_0_y_yyyzz_xyzz[k] * ab_z + g_0_y_yyyzz_xyzzz[k];

                g_0_y_yyyzzz_xzzz[k] = -g_0_y_yyyzz_xzzz[k] * ab_z + g_0_y_yyyzz_xzzzz[k];

                g_0_y_yyyzzz_yyyy[k] = -g_0_y_yyyzz_yyyy[k] * ab_z + g_0_y_yyyzz_yyyyz[k];

                g_0_y_yyyzzz_yyyz[k] = -g_0_y_yyyzz_yyyz[k] * ab_z + g_0_y_yyyzz_yyyzz[k];

                g_0_y_yyyzzz_yyzz[k] = -g_0_y_yyyzz_yyzz[k] * ab_z + g_0_y_yyyzz_yyzzz[k];

                g_0_y_yyyzzz_yzzz[k] = -g_0_y_yyyzz_yzzz[k] * ab_z + g_0_y_yyyzz_yzzzz[k];

                g_0_y_yyyzzz_zzzz[k] = -g_0_y_yyyzz_zzzz[k] * ab_z + g_0_y_yyyzz_zzzzz[k];
            }

            /// Set up 795-810 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzzz_xxxx = cbuffer.data(ig_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxxy = cbuffer.data(ig_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxxz = cbuffer.data(ig_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxyy = cbuffer.data(ig_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxyz = cbuffer.data(ig_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxzz = cbuffer.data(ig_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyyy = cbuffer.data(ig_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyyz = cbuffer.data(ig_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyzz = cbuffer.data(ig_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xzzz = cbuffer.data(ig_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyyy = cbuffer.data(ig_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyyz = cbuffer.data(ig_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyzz = cbuffer.data(ig_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yzzz = cbuffer.data(ig_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_y_yyzzzz_zzzz = cbuffer.data(ig_geom_01_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzzz_xxxx, g_0_y_yyzzz_xxxxz, g_0_y_yyzzz_xxxy, g_0_y_yyzzz_xxxyz, g_0_y_yyzzz_xxxz, g_0_y_yyzzz_xxxzz, g_0_y_yyzzz_xxyy, g_0_y_yyzzz_xxyyz, g_0_y_yyzzz_xxyz, g_0_y_yyzzz_xxyzz, g_0_y_yyzzz_xxzz, g_0_y_yyzzz_xxzzz, g_0_y_yyzzz_xyyy, g_0_y_yyzzz_xyyyz, g_0_y_yyzzz_xyyz, g_0_y_yyzzz_xyyzz, g_0_y_yyzzz_xyzz, g_0_y_yyzzz_xyzzz, g_0_y_yyzzz_xzzz, g_0_y_yyzzz_xzzzz, g_0_y_yyzzz_yyyy, g_0_y_yyzzz_yyyyz, g_0_y_yyzzz_yyyz, g_0_y_yyzzz_yyyzz, g_0_y_yyzzz_yyzz, g_0_y_yyzzz_yyzzz, g_0_y_yyzzz_yzzz, g_0_y_yyzzz_yzzzz, g_0_y_yyzzz_zzzz, g_0_y_yyzzz_zzzzz, g_0_y_yyzzzz_xxxx, g_0_y_yyzzzz_xxxy, g_0_y_yyzzzz_xxxz, g_0_y_yyzzzz_xxyy, g_0_y_yyzzzz_xxyz, g_0_y_yyzzzz_xxzz, g_0_y_yyzzzz_xyyy, g_0_y_yyzzzz_xyyz, g_0_y_yyzzzz_xyzz, g_0_y_yyzzzz_xzzz, g_0_y_yyzzzz_yyyy, g_0_y_yyzzzz_yyyz, g_0_y_yyzzzz_yyzz, g_0_y_yyzzzz_yzzz, g_0_y_yyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzzz_xxxx[k] = -g_0_y_yyzzz_xxxx[k] * ab_z + g_0_y_yyzzz_xxxxz[k];

                g_0_y_yyzzzz_xxxy[k] = -g_0_y_yyzzz_xxxy[k] * ab_z + g_0_y_yyzzz_xxxyz[k];

                g_0_y_yyzzzz_xxxz[k] = -g_0_y_yyzzz_xxxz[k] * ab_z + g_0_y_yyzzz_xxxzz[k];

                g_0_y_yyzzzz_xxyy[k] = -g_0_y_yyzzz_xxyy[k] * ab_z + g_0_y_yyzzz_xxyyz[k];

                g_0_y_yyzzzz_xxyz[k] = -g_0_y_yyzzz_xxyz[k] * ab_z + g_0_y_yyzzz_xxyzz[k];

                g_0_y_yyzzzz_xxzz[k] = -g_0_y_yyzzz_xxzz[k] * ab_z + g_0_y_yyzzz_xxzzz[k];

                g_0_y_yyzzzz_xyyy[k] = -g_0_y_yyzzz_xyyy[k] * ab_z + g_0_y_yyzzz_xyyyz[k];

                g_0_y_yyzzzz_xyyz[k] = -g_0_y_yyzzz_xyyz[k] * ab_z + g_0_y_yyzzz_xyyzz[k];

                g_0_y_yyzzzz_xyzz[k] = -g_0_y_yyzzz_xyzz[k] * ab_z + g_0_y_yyzzz_xyzzz[k];

                g_0_y_yyzzzz_xzzz[k] = -g_0_y_yyzzz_xzzz[k] * ab_z + g_0_y_yyzzz_xzzzz[k];

                g_0_y_yyzzzz_yyyy[k] = -g_0_y_yyzzz_yyyy[k] * ab_z + g_0_y_yyzzz_yyyyz[k];

                g_0_y_yyzzzz_yyyz[k] = -g_0_y_yyzzz_yyyz[k] * ab_z + g_0_y_yyzzz_yyyzz[k];

                g_0_y_yyzzzz_yyzz[k] = -g_0_y_yyzzz_yyzz[k] * ab_z + g_0_y_yyzzz_yyzzz[k];

                g_0_y_yyzzzz_yzzz[k] = -g_0_y_yyzzz_yzzz[k] * ab_z + g_0_y_yyzzz_yzzzz[k];

                g_0_y_yyzzzz_zzzz[k] = -g_0_y_yyzzz_zzzz[k] * ab_z + g_0_y_yyzzz_zzzzz[k];
            }

            /// Set up 810-825 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzzz_xxxx = cbuffer.data(ig_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxxy = cbuffer.data(ig_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxxz = cbuffer.data(ig_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxyy = cbuffer.data(ig_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxyz = cbuffer.data(ig_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxzz = cbuffer.data(ig_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyyy = cbuffer.data(ig_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyyz = cbuffer.data(ig_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyzz = cbuffer.data(ig_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xzzz = cbuffer.data(ig_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyyy = cbuffer.data(ig_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyyz = cbuffer.data(ig_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyzz = cbuffer.data(ig_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yzzz = cbuffer.data(ig_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_y_yzzzzz_zzzz = cbuffer.data(ig_geom_01_off + 824 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzzz_xxxx, g_0_y_yzzzz_xxxxz, g_0_y_yzzzz_xxxy, g_0_y_yzzzz_xxxyz, g_0_y_yzzzz_xxxz, g_0_y_yzzzz_xxxzz, g_0_y_yzzzz_xxyy, g_0_y_yzzzz_xxyyz, g_0_y_yzzzz_xxyz, g_0_y_yzzzz_xxyzz, g_0_y_yzzzz_xxzz, g_0_y_yzzzz_xxzzz, g_0_y_yzzzz_xyyy, g_0_y_yzzzz_xyyyz, g_0_y_yzzzz_xyyz, g_0_y_yzzzz_xyyzz, g_0_y_yzzzz_xyzz, g_0_y_yzzzz_xyzzz, g_0_y_yzzzz_xzzz, g_0_y_yzzzz_xzzzz, g_0_y_yzzzz_yyyy, g_0_y_yzzzz_yyyyz, g_0_y_yzzzz_yyyz, g_0_y_yzzzz_yyyzz, g_0_y_yzzzz_yyzz, g_0_y_yzzzz_yyzzz, g_0_y_yzzzz_yzzz, g_0_y_yzzzz_yzzzz, g_0_y_yzzzz_zzzz, g_0_y_yzzzz_zzzzz, g_0_y_yzzzzz_xxxx, g_0_y_yzzzzz_xxxy, g_0_y_yzzzzz_xxxz, g_0_y_yzzzzz_xxyy, g_0_y_yzzzzz_xxyz, g_0_y_yzzzzz_xxzz, g_0_y_yzzzzz_xyyy, g_0_y_yzzzzz_xyyz, g_0_y_yzzzzz_xyzz, g_0_y_yzzzzz_xzzz, g_0_y_yzzzzz_yyyy, g_0_y_yzzzzz_yyyz, g_0_y_yzzzzz_yyzz, g_0_y_yzzzzz_yzzz, g_0_y_yzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzzz_xxxx[k] = -g_0_y_yzzzz_xxxx[k] * ab_z + g_0_y_yzzzz_xxxxz[k];

                g_0_y_yzzzzz_xxxy[k] = -g_0_y_yzzzz_xxxy[k] * ab_z + g_0_y_yzzzz_xxxyz[k];

                g_0_y_yzzzzz_xxxz[k] = -g_0_y_yzzzz_xxxz[k] * ab_z + g_0_y_yzzzz_xxxzz[k];

                g_0_y_yzzzzz_xxyy[k] = -g_0_y_yzzzz_xxyy[k] * ab_z + g_0_y_yzzzz_xxyyz[k];

                g_0_y_yzzzzz_xxyz[k] = -g_0_y_yzzzz_xxyz[k] * ab_z + g_0_y_yzzzz_xxyzz[k];

                g_0_y_yzzzzz_xxzz[k] = -g_0_y_yzzzz_xxzz[k] * ab_z + g_0_y_yzzzz_xxzzz[k];

                g_0_y_yzzzzz_xyyy[k] = -g_0_y_yzzzz_xyyy[k] * ab_z + g_0_y_yzzzz_xyyyz[k];

                g_0_y_yzzzzz_xyyz[k] = -g_0_y_yzzzz_xyyz[k] * ab_z + g_0_y_yzzzz_xyyzz[k];

                g_0_y_yzzzzz_xyzz[k] = -g_0_y_yzzzz_xyzz[k] * ab_z + g_0_y_yzzzz_xyzzz[k];

                g_0_y_yzzzzz_xzzz[k] = -g_0_y_yzzzz_xzzz[k] * ab_z + g_0_y_yzzzz_xzzzz[k];

                g_0_y_yzzzzz_yyyy[k] = -g_0_y_yzzzz_yyyy[k] * ab_z + g_0_y_yzzzz_yyyyz[k];

                g_0_y_yzzzzz_yyyz[k] = -g_0_y_yzzzz_yyyz[k] * ab_z + g_0_y_yzzzz_yyyzz[k];

                g_0_y_yzzzzz_yyzz[k] = -g_0_y_yzzzz_yyzz[k] * ab_z + g_0_y_yzzzz_yyzzz[k];

                g_0_y_yzzzzz_yzzz[k] = -g_0_y_yzzzz_yzzz[k] * ab_z + g_0_y_yzzzz_yzzzz[k];

                g_0_y_yzzzzz_zzzz[k] = -g_0_y_yzzzz_zzzz[k] * ab_z + g_0_y_yzzzz_zzzzz[k];
            }

            /// Set up 825-840 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzzz_xxxx = cbuffer.data(ig_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxxy = cbuffer.data(ig_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxxz = cbuffer.data(ig_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxyy = cbuffer.data(ig_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxyz = cbuffer.data(ig_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxzz = cbuffer.data(ig_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyyy = cbuffer.data(ig_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyyz = cbuffer.data(ig_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyzz = cbuffer.data(ig_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xzzz = cbuffer.data(ig_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyyy = cbuffer.data(ig_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyyz = cbuffer.data(ig_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyzz = cbuffer.data(ig_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yzzz = cbuffer.data(ig_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_y_zzzzzz_zzzz = cbuffer.data(ig_geom_01_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzzz_xxxx, g_0_y_zzzzz_xxxxz, g_0_y_zzzzz_xxxy, g_0_y_zzzzz_xxxyz, g_0_y_zzzzz_xxxz, g_0_y_zzzzz_xxxzz, g_0_y_zzzzz_xxyy, g_0_y_zzzzz_xxyyz, g_0_y_zzzzz_xxyz, g_0_y_zzzzz_xxyzz, g_0_y_zzzzz_xxzz, g_0_y_zzzzz_xxzzz, g_0_y_zzzzz_xyyy, g_0_y_zzzzz_xyyyz, g_0_y_zzzzz_xyyz, g_0_y_zzzzz_xyyzz, g_0_y_zzzzz_xyzz, g_0_y_zzzzz_xyzzz, g_0_y_zzzzz_xzzz, g_0_y_zzzzz_xzzzz, g_0_y_zzzzz_yyyy, g_0_y_zzzzz_yyyyz, g_0_y_zzzzz_yyyz, g_0_y_zzzzz_yyyzz, g_0_y_zzzzz_yyzz, g_0_y_zzzzz_yyzzz, g_0_y_zzzzz_yzzz, g_0_y_zzzzz_yzzzz, g_0_y_zzzzz_zzzz, g_0_y_zzzzz_zzzzz, g_0_y_zzzzzz_xxxx, g_0_y_zzzzzz_xxxy, g_0_y_zzzzzz_xxxz, g_0_y_zzzzzz_xxyy, g_0_y_zzzzzz_xxyz, g_0_y_zzzzzz_xxzz, g_0_y_zzzzzz_xyyy, g_0_y_zzzzzz_xyyz, g_0_y_zzzzzz_xyzz, g_0_y_zzzzzz_xzzz, g_0_y_zzzzzz_yyyy, g_0_y_zzzzzz_yyyz, g_0_y_zzzzzz_yyzz, g_0_y_zzzzzz_yzzz, g_0_y_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzzz_xxxx[k] = -g_0_y_zzzzz_xxxx[k] * ab_z + g_0_y_zzzzz_xxxxz[k];

                g_0_y_zzzzzz_xxxy[k] = -g_0_y_zzzzz_xxxy[k] * ab_z + g_0_y_zzzzz_xxxyz[k];

                g_0_y_zzzzzz_xxxz[k] = -g_0_y_zzzzz_xxxz[k] * ab_z + g_0_y_zzzzz_xxxzz[k];

                g_0_y_zzzzzz_xxyy[k] = -g_0_y_zzzzz_xxyy[k] * ab_z + g_0_y_zzzzz_xxyyz[k];

                g_0_y_zzzzzz_xxyz[k] = -g_0_y_zzzzz_xxyz[k] * ab_z + g_0_y_zzzzz_xxyzz[k];

                g_0_y_zzzzzz_xxzz[k] = -g_0_y_zzzzz_xxzz[k] * ab_z + g_0_y_zzzzz_xxzzz[k];

                g_0_y_zzzzzz_xyyy[k] = -g_0_y_zzzzz_xyyy[k] * ab_z + g_0_y_zzzzz_xyyyz[k];

                g_0_y_zzzzzz_xyyz[k] = -g_0_y_zzzzz_xyyz[k] * ab_z + g_0_y_zzzzz_xyyzz[k];

                g_0_y_zzzzzz_xyzz[k] = -g_0_y_zzzzz_xyzz[k] * ab_z + g_0_y_zzzzz_xyzzz[k];

                g_0_y_zzzzzz_xzzz[k] = -g_0_y_zzzzz_xzzz[k] * ab_z + g_0_y_zzzzz_xzzzz[k];

                g_0_y_zzzzzz_yyyy[k] = -g_0_y_zzzzz_yyyy[k] * ab_z + g_0_y_zzzzz_yyyyz[k];

                g_0_y_zzzzzz_yyyz[k] = -g_0_y_zzzzz_yyyz[k] * ab_z + g_0_y_zzzzz_yyyzz[k];

                g_0_y_zzzzzz_yyzz[k] = -g_0_y_zzzzz_yyzz[k] * ab_z + g_0_y_zzzzz_yyzzz[k];

                g_0_y_zzzzzz_yzzz[k] = -g_0_y_zzzzz_yzzz[k] * ab_z + g_0_y_zzzzz_yzzzz[k];

                g_0_y_zzzzzz_zzzz[k] = -g_0_y_zzzzz_zzzz[k] * ab_z + g_0_y_zzzzz_zzzzz[k];
            }

            /// Set up 840-855 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxx_xxxx = cbuffer.data(ig_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxxy = cbuffer.data(ig_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxxz = cbuffer.data(ig_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxyy = cbuffer.data(ig_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxyz = cbuffer.data(ig_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxzz = cbuffer.data(ig_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyyy = cbuffer.data(ig_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyyz = cbuffer.data(ig_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyzz = cbuffer.data(ig_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xzzz = cbuffer.data(ig_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyyy = cbuffer.data(ig_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyyz = cbuffer.data(ig_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyzz = cbuffer.data(ig_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yzzz = cbuffer.data(ig_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_z_xxxxxx_zzzz = cbuffer.data(ig_geom_01_off + 854 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxx_xxxx, g_0_z_xxxxx_xxxxx, g_0_z_xxxxx_xxxxy, g_0_z_xxxxx_xxxxz, g_0_z_xxxxx_xxxy, g_0_z_xxxxx_xxxyy, g_0_z_xxxxx_xxxyz, g_0_z_xxxxx_xxxz, g_0_z_xxxxx_xxxzz, g_0_z_xxxxx_xxyy, g_0_z_xxxxx_xxyyy, g_0_z_xxxxx_xxyyz, g_0_z_xxxxx_xxyz, g_0_z_xxxxx_xxyzz, g_0_z_xxxxx_xxzz, g_0_z_xxxxx_xxzzz, g_0_z_xxxxx_xyyy, g_0_z_xxxxx_xyyyy, g_0_z_xxxxx_xyyyz, g_0_z_xxxxx_xyyz, g_0_z_xxxxx_xyyzz, g_0_z_xxxxx_xyzz, g_0_z_xxxxx_xyzzz, g_0_z_xxxxx_xzzz, g_0_z_xxxxx_xzzzz, g_0_z_xxxxx_yyyy, g_0_z_xxxxx_yyyz, g_0_z_xxxxx_yyzz, g_0_z_xxxxx_yzzz, g_0_z_xxxxx_zzzz, g_0_z_xxxxxx_xxxx, g_0_z_xxxxxx_xxxy, g_0_z_xxxxxx_xxxz, g_0_z_xxxxxx_xxyy, g_0_z_xxxxxx_xxyz, g_0_z_xxxxxx_xxzz, g_0_z_xxxxxx_xyyy, g_0_z_xxxxxx_xyyz, g_0_z_xxxxxx_xyzz, g_0_z_xxxxxx_xzzz, g_0_z_xxxxxx_yyyy, g_0_z_xxxxxx_yyyz, g_0_z_xxxxxx_yyzz, g_0_z_xxxxxx_yzzz, g_0_z_xxxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxx_xxxx[k] = -g_0_z_xxxxx_xxxx[k] * ab_x + g_0_z_xxxxx_xxxxx[k];

                g_0_z_xxxxxx_xxxy[k] = -g_0_z_xxxxx_xxxy[k] * ab_x + g_0_z_xxxxx_xxxxy[k];

                g_0_z_xxxxxx_xxxz[k] = -g_0_z_xxxxx_xxxz[k] * ab_x + g_0_z_xxxxx_xxxxz[k];

                g_0_z_xxxxxx_xxyy[k] = -g_0_z_xxxxx_xxyy[k] * ab_x + g_0_z_xxxxx_xxxyy[k];

                g_0_z_xxxxxx_xxyz[k] = -g_0_z_xxxxx_xxyz[k] * ab_x + g_0_z_xxxxx_xxxyz[k];

                g_0_z_xxxxxx_xxzz[k] = -g_0_z_xxxxx_xxzz[k] * ab_x + g_0_z_xxxxx_xxxzz[k];

                g_0_z_xxxxxx_xyyy[k] = -g_0_z_xxxxx_xyyy[k] * ab_x + g_0_z_xxxxx_xxyyy[k];

                g_0_z_xxxxxx_xyyz[k] = -g_0_z_xxxxx_xyyz[k] * ab_x + g_0_z_xxxxx_xxyyz[k];

                g_0_z_xxxxxx_xyzz[k] = -g_0_z_xxxxx_xyzz[k] * ab_x + g_0_z_xxxxx_xxyzz[k];

                g_0_z_xxxxxx_xzzz[k] = -g_0_z_xxxxx_xzzz[k] * ab_x + g_0_z_xxxxx_xxzzz[k];

                g_0_z_xxxxxx_yyyy[k] = -g_0_z_xxxxx_yyyy[k] * ab_x + g_0_z_xxxxx_xyyyy[k];

                g_0_z_xxxxxx_yyyz[k] = -g_0_z_xxxxx_yyyz[k] * ab_x + g_0_z_xxxxx_xyyyz[k];

                g_0_z_xxxxxx_yyzz[k] = -g_0_z_xxxxx_yyzz[k] * ab_x + g_0_z_xxxxx_xyyzz[k];

                g_0_z_xxxxxx_yzzz[k] = -g_0_z_xxxxx_yzzz[k] * ab_x + g_0_z_xxxxx_xyzzz[k];

                g_0_z_xxxxxx_zzzz[k] = -g_0_z_xxxxx_zzzz[k] * ab_x + g_0_z_xxxxx_xzzzz[k];
            }

            /// Set up 855-870 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxy_xxxx = cbuffer.data(ig_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxxy = cbuffer.data(ig_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxxz = cbuffer.data(ig_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxyy = cbuffer.data(ig_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxyz = cbuffer.data(ig_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxzz = cbuffer.data(ig_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyyy = cbuffer.data(ig_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyyz = cbuffer.data(ig_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyzz = cbuffer.data(ig_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xzzz = cbuffer.data(ig_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyyy = cbuffer.data(ig_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyyz = cbuffer.data(ig_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyzz = cbuffer.data(ig_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yzzz = cbuffer.data(ig_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_z_xxxxxy_zzzz = cbuffer.data(ig_geom_01_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxy_xxxx, g_0_z_xxxxxy_xxxy, g_0_z_xxxxxy_xxxz, g_0_z_xxxxxy_xxyy, g_0_z_xxxxxy_xxyz, g_0_z_xxxxxy_xxzz, g_0_z_xxxxxy_xyyy, g_0_z_xxxxxy_xyyz, g_0_z_xxxxxy_xyzz, g_0_z_xxxxxy_xzzz, g_0_z_xxxxxy_yyyy, g_0_z_xxxxxy_yyyz, g_0_z_xxxxxy_yyzz, g_0_z_xxxxxy_yzzz, g_0_z_xxxxxy_zzzz, g_0_z_xxxxy_xxxx, g_0_z_xxxxy_xxxxx, g_0_z_xxxxy_xxxxy, g_0_z_xxxxy_xxxxz, g_0_z_xxxxy_xxxy, g_0_z_xxxxy_xxxyy, g_0_z_xxxxy_xxxyz, g_0_z_xxxxy_xxxz, g_0_z_xxxxy_xxxzz, g_0_z_xxxxy_xxyy, g_0_z_xxxxy_xxyyy, g_0_z_xxxxy_xxyyz, g_0_z_xxxxy_xxyz, g_0_z_xxxxy_xxyzz, g_0_z_xxxxy_xxzz, g_0_z_xxxxy_xxzzz, g_0_z_xxxxy_xyyy, g_0_z_xxxxy_xyyyy, g_0_z_xxxxy_xyyyz, g_0_z_xxxxy_xyyz, g_0_z_xxxxy_xyyzz, g_0_z_xxxxy_xyzz, g_0_z_xxxxy_xyzzz, g_0_z_xxxxy_xzzz, g_0_z_xxxxy_xzzzz, g_0_z_xxxxy_yyyy, g_0_z_xxxxy_yyyz, g_0_z_xxxxy_yyzz, g_0_z_xxxxy_yzzz, g_0_z_xxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxy_xxxx[k] = -g_0_z_xxxxy_xxxx[k] * ab_x + g_0_z_xxxxy_xxxxx[k];

                g_0_z_xxxxxy_xxxy[k] = -g_0_z_xxxxy_xxxy[k] * ab_x + g_0_z_xxxxy_xxxxy[k];

                g_0_z_xxxxxy_xxxz[k] = -g_0_z_xxxxy_xxxz[k] * ab_x + g_0_z_xxxxy_xxxxz[k];

                g_0_z_xxxxxy_xxyy[k] = -g_0_z_xxxxy_xxyy[k] * ab_x + g_0_z_xxxxy_xxxyy[k];

                g_0_z_xxxxxy_xxyz[k] = -g_0_z_xxxxy_xxyz[k] * ab_x + g_0_z_xxxxy_xxxyz[k];

                g_0_z_xxxxxy_xxzz[k] = -g_0_z_xxxxy_xxzz[k] * ab_x + g_0_z_xxxxy_xxxzz[k];

                g_0_z_xxxxxy_xyyy[k] = -g_0_z_xxxxy_xyyy[k] * ab_x + g_0_z_xxxxy_xxyyy[k];

                g_0_z_xxxxxy_xyyz[k] = -g_0_z_xxxxy_xyyz[k] * ab_x + g_0_z_xxxxy_xxyyz[k];

                g_0_z_xxxxxy_xyzz[k] = -g_0_z_xxxxy_xyzz[k] * ab_x + g_0_z_xxxxy_xxyzz[k];

                g_0_z_xxxxxy_xzzz[k] = -g_0_z_xxxxy_xzzz[k] * ab_x + g_0_z_xxxxy_xxzzz[k];

                g_0_z_xxxxxy_yyyy[k] = -g_0_z_xxxxy_yyyy[k] * ab_x + g_0_z_xxxxy_xyyyy[k];

                g_0_z_xxxxxy_yyyz[k] = -g_0_z_xxxxy_yyyz[k] * ab_x + g_0_z_xxxxy_xyyyz[k];

                g_0_z_xxxxxy_yyzz[k] = -g_0_z_xxxxy_yyzz[k] * ab_x + g_0_z_xxxxy_xyyzz[k];

                g_0_z_xxxxxy_yzzz[k] = -g_0_z_xxxxy_yzzz[k] * ab_x + g_0_z_xxxxy_xyzzz[k];

                g_0_z_xxxxxy_zzzz[k] = -g_0_z_xxxxy_zzzz[k] * ab_x + g_0_z_xxxxy_xzzzz[k];
            }

            /// Set up 870-885 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxz_xxxx = cbuffer.data(ig_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxxy = cbuffer.data(ig_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxxz = cbuffer.data(ig_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxyy = cbuffer.data(ig_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxyz = cbuffer.data(ig_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxzz = cbuffer.data(ig_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyyy = cbuffer.data(ig_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyyz = cbuffer.data(ig_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyzz = cbuffer.data(ig_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xzzz = cbuffer.data(ig_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyyy = cbuffer.data(ig_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyyz = cbuffer.data(ig_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyzz = cbuffer.data(ig_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yzzz = cbuffer.data(ig_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_z_xxxxxz_zzzz = cbuffer.data(ig_geom_01_off + 884 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxz_xxxx, g_0_z_xxxxxz_xxxy, g_0_z_xxxxxz_xxxz, g_0_z_xxxxxz_xxyy, g_0_z_xxxxxz_xxyz, g_0_z_xxxxxz_xxzz, g_0_z_xxxxxz_xyyy, g_0_z_xxxxxz_xyyz, g_0_z_xxxxxz_xyzz, g_0_z_xxxxxz_xzzz, g_0_z_xxxxxz_yyyy, g_0_z_xxxxxz_yyyz, g_0_z_xxxxxz_yyzz, g_0_z_xxxxxz_yzzz, g_0_z_xxxxxz_zzzz, g_0_z_xxxxz_xxxx, g_0_z_xxxxz_xxxxx, g_0_z_xxxxz_xxxxy, g_0_z_xxxxz_xxxxz, g_0_z_xxxxz_xxxy, g_0_z_xxxxz_xxxyy, g_0_z_xxxxz_xxxyz, g_0_z_xxxxz_xxxz, g_0_z_xxxxz_xxxzz, g_0_z_xxxxz_xxyy, g_0_z_xxxxz_xxyyy, g_0_z_xxxxz_xxyyz, g_0_z_xxxxz_xxyz, g_0_z_xxxxz_xxyzz, g_0_z_xxxxz_xxzz, g_0_z_xxxxz_xxzzz, g_0_z_xxxxz_xyyy, g_0_z_xxxxz_xyyyy, g_0_z_xxxxz_xyyyz, g_0_z_xxxxz_xyyz, g_0_z_xxxxz_xyyzz, g_0_z_xxxxz_xyzz, g_0_z_xxxxz_xyzzz, g_0_z_xxxxz_xzzz, g_0_z_xxxxz_xzzzz, g_0_z_xxxxz_yyyy, g_0_z_xxxxz_yyyz, g_0_z_xxxxz_yyzz, g_0_z_xxxxz_yzzz, g_0_z_xxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxz_xxxx[k] = -g_0_z_xxxxz_xxxx[k] * ab_x + g_0_z_xxxxz_xxxxx[k];

                g_0_z_xxxxxz_xxxy[k] = -g_0_z_xxxxz_xxxy[k] * ab_x + g_0_z_xxxxz_xxxxy[k];

                g_0_z_xxxxxz_xxxz[k] = -g_0_z_xxxxz_xxxz[k] * ab_x + g_0_z_xxxxz_xxxxz[k];

                g_0_z_xxxxxz_xxyy[k] = -g_0_z_xxxxz_xxyy[k] * ab_x + g_0_z_xxxxz_xxxyy[k];

                g_0_z_xxxxxz_xxyz[k] = -g_0_z_xxxxz_xxyz[k] * ab_x + g_0_z_xxxxz_xxxyz[k];

                g_0_z_xxxxxz_xxzz[k] = -g_0_z_xxxxz_xxzz[k] * ab_x + g_0_z_xxxxz_xxxzz[k];

                g_0_z_xxxxxz_xyyy[k] = -g_0_z_xxxxz_xyyy[k] * ab_x + g_0_z_xxxxz_xxyyy[k];

                g_0_z_xxxxxz_xyyz[k] = -g_0_z_xxxxz_xyyz[k] * ab_x + g_0_z_xxxxz_xxyyz[k];

                g_0_z_xxxxxz_xyzz[k] = -g_0_z_xxxxz_xyzz[k] * ab_x + g_0_z_xxxxz_xxyzz[k];

                g_0_z_xxxxxz_xzzz[k] = -g_0_z_xxxxz_xzzz[k] * ab_x + g_0_z_xxxxz_xxzzz[k];

                g_0_z_xxxxxz_yyyy[k] = -g_0_z_xxxxz_yyyy[k] * ab_x + g_0_z_xxxxz_xyyyy[k];

                g_0_z_xxxxxz_yyyz[k] = -g_0_z_xxxxz_yyyz[k] * ab_x + g_0_z_xxxxz_xyyyz[k];

                g_0_z_xxxxxz_yyzz[k] = -g_0_z_xxxxz_yyzz[k] * ab_x + g_0_z_xxxxz_xyyzz[k];

                g_0_z_xxxxxz_yzzz[k] = -g_0_z_xxxxz_yzzz[k] * ab_x + g_0_z_xxxxz_xyzzz[k];

                g_0_z_xxxxxz_zzzz[k] = -g_0_z_xxxxz_zzzz[k] * ab_x + g_0_z_xxxxz_xzzzz[k];
            }

            /// Set up 885-900 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyy_xxxx = cbuffer.data(ig_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxxy = cbuffer.data(ig_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxxz = cbuffer.data(ig_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxyy = cbuffer.data(ig_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxyz = cbuffer.data(ig_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxzz = cbuffer.data(ig_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyyy = cbuffer.data(ig_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyyz = cbuffer.data(ig_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyzz = cbuffer.data(ig_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xzzz = cbuffer.data(ig_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyyy = cbuffer.data(ig_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyyz = cbuffer.data(ig_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyzz = cbuffer.data(ig_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yzzz = cbuffer.data(ig_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_z_xxxxyy_zzzz = cbuffer.data(ig_geom_01_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyy_xxxx, g_0_z_xxxxyy_xxxy, g_0_z_xxxxyy_xxxz, g_0_z_xxxxyy_xxyy, g_0_z_xxxxyy_xxyz, g_0_z_xxxxyy_xxzz, g_0_z_xxxxyy_xyyy, g_0_z_xxxxyy_xyyz, g_0_z_xxxxyy_xyzz, g_0_z_xxxxyy_xzzz, g_0_z_xxxxyy_yyyy, g_0_z_xxxxyy_yyyz, g_0_z_xxxxyy_yyzz, g_0_z_xxxxyy_yzzz, g_0_z_xxxxyy_zzzz, g_0_z_xxxyy_xxxx, g_0_z_xxxyy_xxxxx, g_0_z_xxxyy_xxxxy, g_0_z_xxxyy_xxxxz, g_0_z_xxxyy_xxxy, g_0_z_xxxyy_xxxyy, g_0_z_xxxyy_xxxyz, g_0_z_xxxyy_xxxz, g_0_z_xxxyy_xxxzz, g_0_z_xxxyy_xxyy, g_0_z_xxxyy_xxyyy, g_0_z_xxxyy_xxyyz, g_0_z_xxxyy_xxyz, g_0_z_xxxyy_xxyzz, g_0_z_xxxyy_xxzz, g_0_z_xxxyy_xxzzz, g_0_z_xxxyy_xyyy, g_0_z_xxxyy_xyyyy, g_0_z_xxxyy_xyyyz, g_0_z_xxxyy_xyyz, g_0_z_xxxyy_xyyzz, g_0_z_xxxyy_xyzz, g_0_z_xxxyy_xyzzz, g_0_z_xxxyy_xzzz, g_0_z_xxxyy_xzzzz, g_0_z_xxxyy_yyyy, g_0_z_xxxyy_yyyz, g_0_z_xxxyy_yyzz, g_0_z_xxxyy_yzzz, g_0_z_xxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyy_xxxx[k] = -g_0_z_xxxyy_xxxx[k] * ab_x + g_0_z_xxxyy_xxxxx[k];

                g_0_z_xxxxyy_xxxy[k] = -g_0_z_xxxyy_xxxy[k] * ab_x + g_0_z_xxxyy_xxxxy[k];

                g_0_z_xxxxyy_xxxz[k] = -g_0_z_xxxyy_xxxz[k] * ab_x + g_0_z_xxxyy_xxxxz[k];

                g_0_z_xxxxyy_xxyy[k] = -g_0_z_xxxyy_xxyy[k] * ab_x + g_0_z_xxxyy_xxxyy[k];

                g_0_z_xxxxyy_xxyz[k] = -g_0_z_xxxyy_xxyz[k] * ab_x + g_0_z_xxxyy_xxxyz[k];

                g_0_z_xxxxyy_xxzz[k] = -g_0_z_xxxyy_xxzz[k] * ab_x + g_0_z_xxxyy_xxxzz[k];

                g_0_z_xxxxyy_xyyy[k] = -g_0_z_xxxyy_xyyy[k] * ab_x + g_0_z_xxxyy_xxyyy[k];

                g_0_z_xxxxyy_xyyz[k] = -g_0_z_xxxyy_xyyz[k] * ab_x + g_0_z_xxxyy_xxyyz[k];

                g_0_z_xxxxyy_xyzz[k] = -g_0_z_xxxyy_xyzz[k] * ab_x + g_0_z_xxxyy_xxyzz[k];

                g_0_z_xxxxyy_xzzz[k] = -g_0_z_xxxyy_xzzz[k] * ab_x + g_0_z_xxxyy_xxzzz[k];

                g_0_z_xxxxyy_yyyy[k] = -g_0_z_xxxyy_yyyy[k] * ab_x + g_0_z_xxxyy_xyyyy[k];

                g_0_z_xxxxyy_yyyz[k] = -g_0_z_xxxyy_yyyz[k] * ab_x + g_0_z_xxxyy_xyyyz[k];

                g_0_z_xxxxyy_yyzz[k] = -g_0_z_xxxyy_yyzz[k] * ab_x + g_0_z_xxxyy_xyyzz[k];

                g_0_z_xxxxyy_yzzz[k] = -g_0_z_xxxyy_yzzz[k] * ab_x + g_0_z_xxxyy_xyzzz[k];

                g_0_z_xxxxyy_zzzz[k] = -g_0_z_xxxyy_zzzz[k] * ab_x + g_0_z_xxxyy_xzzzz[k];
            }

            /// Set up 900-915 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyz_xxxx = cbuffer.data(ig_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxxy = cbuffer.data(ig_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxxz = cbuffer.data(ig_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxyy = cbuffer.data(ig_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxyz = cbuffer.data(ig_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxzz = cbuffer.data(ig_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyyy = cbuffer.data(ig_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyyz = cbuffer.data(ig_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyzz = cbuffer.data(ig_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xzzz = cbuffer.data(ig_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyyy = cbuffer.data(ig_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyyz = cbuffer.data(ig_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyzz = cbuffer.data(ig_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yzzz = cbuffer.data(ig_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_z_xxxxyz_zzzz = cbuffer.data(ig_geom_01_off + 914 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyz_xxxx, g_0_z_xxxxyz_xxxy, g_0_z_xxxxyz_xxxz, g_0_z_xxxxyz_xxyy, g_0_z_xxxxyz_xxyz, g_0_z_xxxxyz_xxzz, g_0_z_xxxxyz_xyyy, g_0_z_xxxxyz_xyyz, g_0_z_xxxxyz_xyzz, g_0_z_xxxxyz_xzzz, g_0_z_xxxxyz_yyyy, g_0_z_xxxxyz_yyyz, g_0_z_xxxxyz_yyzz, g_0_z_xxxxyz_yzzz, g_0_z_xxxxyz_zzzz, g_0_z_xxxyz_xxxx, g_0_z_xxxyz_xxxxx, g_0_z_xxxyz_xxxxy, g_0_z_xxxyz_xxxxz, g_0_z_xxxyz_xxxy, g_0_z_xxxyz_xxxyy, g_0_z_xxxyz_xxxyz, g_0_z_xxxyz_xxxz, g_0_z_xxxyz_xxxzz, g_0_z_xxxyz_xxyy, g_0_z_xxxyz_xxyyy, g_0_z_xxxyz_xxyyz, g_0_z_xxxyz_xxyz, g_0_z_xxxyz_xxyzz, g_0_z_xxxyz_xxzz, g_0_z_xxxyz_xxzzz, g_0_z_xxxyz_xyyy, g_0_z_xxxyz_xyyyy, g_0_z_xxxyz_xyyyz, g_0_z_xxxyz_xyyz, g_0_z_xxxyz_xyyzz, g_0_z_xxxyz_xyzz, g_0_z_xxxyz_xyzzz, g_0_z_xxxyz_xzzz, g_0_z_xxxyz_xzzzz, g_0_z_xxxyz_yyyy, g_0_z_xxxyz_yyyz, g_0_z_xxxyz_yyzz, g_0_z_xxxyz_yzzz, g_0_z_xxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyz_xxxx[k] = -g_0_z_xxxyz_xxxx[k] * ab_x + g_0_z_xxxyz_xxxxx[k];

                g_0_z_xxxxyz_xxxy[k] = -g_0_z_xxxyz_xxxy[k] * ab_x + g_0_z_xxxyz_xxxxy[k];

                g_0_z_xxxxyz_xxxz[k] = -g_0_z_xxxyz_xxxz[k] * ab_x + g_0_z_xxxyz_xxxxz[k];

                g_0_z_xxxxyz_xxyy[k] = -g_0_z_xxxyz_xxyy[k] * ab_x + g_0_z_xxxyz_xxxyy[k];

                g_0_z_xxxxyz_xxyz[k] = -g_0_z_xxxyz_xxyz[k] * ab_x + g_0_z_xxxyz_xxxyz[k];

                g_0_z_xxxxyz_xxzz[k] = -g_0_z_xxxyz_xxzz[k] * ab_x + g_0_z_xxxyz_xxxzz[k];

                g_0_z_xxxxyz_xyyy[k] = -g_0_z_xxxyz_xyyy[k] * ab_x + g_0_z_xxxyz_xxyyy[k];

                g_0_z_xxxxyz_xyyz[k] = -g_0_z_xxxyz_xyyz[k] * ab_x + g_0_z_xxxyz_xxyyz[k];

                g_0_z_xxxxyz_xyzz[k] = -g_0_z_xxxyz_xyzz[k] * ab_x + g_0_z_xxxyz_xxyzz[k];

                g_0_z_xxxxyz_xzzz[k] = -g_0_z_xxxyz_xzzz[k] * ab_x + g_0_z_xxxyz_xxzzz[k];

                g_0_z_xxxxyz_yyyy[k] = -g_0_z_xxxyz_yyyy[k] * ab_x + g_0_z_xxxyz_xyyyy[k];

                g_0_z_xxxxyz_yyyz[k] = -g_0_z_xxxyz_yyyz[k] * ab_x + g_0_z_xxxyz_xyyyz[k];

                g_0_z_xxxxyz_yyzz[k] = -g_0_z_xxxyz_yyzz[k] * ab_x + g_0_z_xxxyz_xyyzz[k];

                g_0_z_xxxxyz_yzzz[k] = -g_0_z_xxxyz_yzzz[k] * ab_x + g_0_z_xxxyz_xyzzz[k];

                g_0_z_xxxxyz_zzzz[k] = -g_0_z_xxxyz_zzzz[k] * ab_x + g_0_z_xxxyz_xzzzz[k];
            }

            /// Set up 915-930 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxzz_xxxx = cbuffer.data(ig_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxxy = cbuffer.data(ig_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxxz = cbuffer.data(ig_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxyy = cbuffer.data(ig_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxyz = cbuffer.data(ig_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxzz = cbuffer.data(ig_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyyy = cbuffer.data(ig_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyyz = cbuffer.data(ig_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyzz = cbuffer.data(ig_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xzzz = cbuffer.data(ig_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyyy = cbuffer.data(ig_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyyz = cbuffer.data(ig_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyzz = cbuffer.data(ig_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yzzz = cbuffer.data(ig_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_z_xxxxzz_zzzz = cbuffer.data(ig_geom_01_off + 929 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxzz_xxxx, g_0_z_xxxxzz_xxxy, g_0_z_xxxxzz_xxxz, g_0_z_xxxxzz_xxyy, g_0_z_xxxxzz_xxyz, g_0_z_xxxxzz_xxzz, g_0_z_xxxxzz_xyyy, g_0_z_xxxxzz_xyyz, g_0_z_xxxxzz_xyzz, g_0_z_xxxxzz_xzzz, g_0_z_xxxxzz_yyyy, g_0_z_xxxxzz_yyyz, g_0_z_xxxxzz_yyzz, g_0_z_xxxxzz_yzzz, g_0_z_xxxxzz_zzzz, g_0_z_xxxzz_xxxx, g_0_z_xxxzz_xxxxx, g_0_z_xxxzz_xxxxy, g_0_z_xxxzz_xxxxz, g_0_z_xxxzz_xxxy, g_0_z_xxxzz_xxxyy, g_0_z_xxxzz_xxxyz, g_0_z_xxxzz_xxxz, g_0_z_xxxzz_xxxzz, g_0_z_xxxzz_xxyy, g_0_z_xxxzz_xxyyy, g_0_z_xxxzz_xxyyz, g_0_z_xxxzz_xxyz, g_0_z_xxxzz_xxyzz, g_0_z_xxxzz_xxzz, g_0_z_xxxzz_xxzzz, g_0_z_xxxzz_xyyy, g_0_z_xxxzz_xyyyy, g_0_z_xxxzz_xyyyz, g_0_z_xxxzz_xyyz, g_0_z_xxxzz_xyyzz, g_0_z_xxxzz_xyzz, g_0_z_xxxzz_xyzzz, g_0_z_xxxzz_xzzz, g_0_z_xxxzz_xzzzz, g_0_z_xxxzz_yyyy, g_0_z_xxxzz_yyyz, g_0_z_xxxzz_yyzz, g_0_z_xxxzz_yzzz, g_0_z_xxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxzz_xxxx[k] = -g_0_z_xxxzz_xxxx[k] * ab_x + g_0_z_xxxzz_xxxxx[k];

                g_0_z_xxxxzz_xxxy[k] = -g_0_z_xxxzz_xxxy[k] * ab_x + g_0_z_xxxzz_xxxxy[k];

                g_0_z_xxxxzz_xxxz[k] = -g_0_z_xxxzz_xxxz[k] * ab_x + g_0_z_xxxzz_xxxxz[k];

                g_0_z_xxxxzz_xxyy[k] = -g_0_z_xxxzz_xxyy[k] * ab_x + g_0_z_xxxzz_xxxyy[k];

                g_0_z_xxxxzz_xxyz[k] = -g_0_z_xxxzz_xxyz[k] * ab_x + g_0_z_xxxzz_xxxyz[k];

                g_0_z_xxxxzz_xxzz[k] = -g_0_z_xxxzz_xxzz[k] * ab_x + g_0_z_xxxzz_xxxzz[k];

                g_0_z_xxxxzz_xyyy[k] = -g_0_z_xxxzz_xyyy[k] * ab_x + g_0_z_xxxzz_xxyyy[k];

                g_0_z_xxxxzz_xyyz[k] = -g_0_z_xxxzz_xyyz[k] * ab_x + g_0_z_xxxzz_xxyyz[k];

                g_0_z_xxxxzz_xyzz[k] = -g_0_z_xxxzz_xyzz[k] * ab_x + g_0_z_xxxzz_xxyzz[k];

                g_0_z_xxxxzz_xzzz[k] = -g_0_z_xxxzz_xzzz[k] * ab_x + g_0_z_xxxzz_xxzzz[k];

                g_0_z_xxxxzz_yyyy[k] = -g_0_z_xxxzz_yyyy[k] * ab_x + g_0_z_xxxzz_xyyyy[k];

                g_0_z_xxxxzz_yyyz[k] = -g_0_z_xxxzz_yyyz[k] * ab_x + g_0_z_xxxzz_xyyyz[k];

                g_0_z_xxxxzz_yyzz[k] = -g_0_z_xxxzz_yyzz[k] * ab_x + g_0_z_xxxzz_xyyzz[k];

                g_0_z_xxxxzz_yzzz[k] = -g_0_z_xxxzz_yzzz[k] * ab_x + g_0_z_xxxzz_xyzzz[k];

                g_0_z_xxxxzz_zzzz[k] = -g_0_z_xxxzz_zzzz[k] * ab_x + g_0_z_xxxzz_xzzzz[k];
            }

            /// Set up 930-945 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyy_xxxx = cbuffer.data(ig_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxxy = cbuffer.data(ig_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxxz = cbuffer.data(ig_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxyy = cbuffer.data(ig_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxyz = cbuffer.data(ig_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxzz = cbuffer.data(ig_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyyy = cbuffer.data(ig_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyyz = cbuffer.data(ig_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyzz = cbuffer.data(ig_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xzzz = cbuffer.data(ig_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyyy = cbuffer.data(ig_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyyz = cbuffer.data(ig_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyzz = cbuffer.data(ig_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yzzz = cbuffer.data(ig_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_z_xxxyyy_zzzz = cbuffer.data(ig_geom_01_off + 944 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyy_xxxx, g_0_z_xxxyyy_xxxy, g_0_z_xxxyyy_xxxz, g_0_z_xxxyyy_xxyy, g_0_z_xxxyyy_xxyz, g_0_z_xxxyyy_xxzz, g_0_z_xxxyyy_xyyy, g_0_z_xxxyyy_xyyz, g_0_z_xxxyyy_xyzz, g_0_z_xxxyyy_xzzz, g_0_z_xxxyyy_yyyy, g_0_z_xxxyyy_yyyz, g_0_z_xxxyyy_yyzz, g_0_z_xxxyyy_yzzz, g_0_z_xxxyyy_zzzz, g_0_z_xxyyy_xxxx, g_0_z_xxyyy_xxxxx, g_0_z_xxyyy_xxxxy, g_0_z_xxyyy_xxxxz, g_0_z_xxyyy_xxxy, g_0_z_xxyyy_xxxyy, g_0_z_xxyyy_xxxyz, g_0_z_xxyyy_xxxz, g_0_z_xxyyy_xxxzz, g_0_z_xxyyy_xxyy, g_0_z_xxyyy_xxyyy, g_0_z_xxyyy_xxyyz, g_0_z_xxyyy_xxyz, g_0_z_xxyyy_xxyzz, g_0_z_xxyyy_xxzz, g_0_z_xxyyy_xxzzz, g_0_z_xxyyy_xyyy, g_0_z_xxyyy_xyyyy, g_0_z_xxyyy_xyyyz, g_0_z_xxyyy_xyyz, g_0_z_xxyyy_xyyzz, g_0_z_xxyyy_xyzz, g_0_z_xxyyy_xyzzz, g_0_z_xxyyy_xzzz, g_0_z_xxyyy_xzzzz, g_0_z_xxyyy_yyyy, g_0_z_xxyyy_yyyz, g_0_z_xxyyy_yyzz, g_0_z_xxyyy_yzzz, g_0_z_xxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyy_xxxx[k] = -g_0_z_xxyyy_xxxx[k] * ab_x + g_0_z_xxyyy_xxxxx[k];

                g_0_z_xxxyyy_xxxy[k] = -g_0_z_xxyyy_xxxy[k] * ab_x + g_0_z_xxyyy_xxxxy[k];

                g_0_z_xxxyyy_xxxz[k] = -g_0_z_xxyyy_xxxz[k] * ab_x + g_0_z_xxyyy_xxxxz[k];

                g_0_z_xxxyyy_xxyy[k] = -g_0_z_xxyyy_xxyy[k] * ab_x + g_0_z_xxyyy_xxxyy[k];

                g_0_z_xxxyyy_xxyz[k] = -g_0_z_xxyyy_xxyz[k] * ab_x + g_0_z_xxyyy_xxxyz[k];

                g_0_z_xxxyyy_xxzz[k] = -g_0_z_xxyyy_xxzz[k] * ab_x + g_0_z_xxyyy_xxxzz[k];

                g_0_z_xxxyyy_xyyy[k] = -g_0_z_xxyyy_xyyy[k] * ab_x + g_0_z_xxyyy_xxyyy[k];

                g_0_z_xxxyyy_xyyz[k] = -g_0_z_xxyyy_xyyz[k] * ab_x + g_0_z_xxyyy_xxyyz[k];

                g_0_z_xxxyyy_xyzz[k] = -g_0_z_xxyyy_xyzz[k] * ab_x + g_0_z_xxyyy_xxyzz[k];

                g_0_z_xxxyyy_xzzz[k] = -g_0_z_xxyyy_xzzz[k] * ab_x + g_0_z_xxyyy_xxzzz[k];

                g_0_z_xxxyyy_yyyy[k] = -g_0_z_xxyyy_yyyy[k] * ab_x + g_0_z_xxyyy_xyyyy[k];

                g_0_z_xxxyyy_yyyz[k] = -g_0_z_xxyyy_yyyz[k] * ab_x + g_0_z_xxyyy_xyyyz[k];

                g_0_z_xxxyyy_yyzz[k] = -g_0_z_xxyyy_yyzz[k] * ab_x + g_0_z_xxyyy_xyyzz[k];

                g_0_z_xxxyyy_yzzz[k] = -g_0_z_xxyyy_yzzz[k] * ab_x + g_0_z_xxyyy_xyzzz[k];

                g_0_z_xxxyyy_zzzz[k] = -g_0_z_xxyyy_zzzz[k] * ab_x + g_0_z_xxyyy_xzzzz[k];
            }

            /// Set up 945-960 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyz_xxxx = cbuffer.data(ig_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxxy = cbuffer.data(ig_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxxz = cbuffer.data(ig_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxyy = cbuffer.data(ig_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxyz = cbuffer.data(ig_geom_01_off + 949 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxzz = cbuffer.data(ig_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyyy = cbuffer.data(ig_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyyz = cbuffer.data(ig_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyzz = cbuffer.data(ig_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xzzz = cbuffer.data(ig_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyyy = cbuffer.data(ig_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyyz = cbuffer.data(ig_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyzz = cbuffer.data(ig_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yzzz = cbuffer.data(ig_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_z_xxxyyz_zzzz = cbuffer.data(ig_geom_01_off + 959 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyz_xxxx, g_0_z_xxxyyz_xxxy, g_0_z_xxxyyz_xxxz, g_0_z_xxxyyz_xxyy, g_0_z_xxxyyz_xxyz, g_0_z_xxxyyz_xxzz, g_0_z_xxxyyz_xyyy, g_0_z_xxxyyz_xyyz, g_0_z_xxxyyz_xyzz, g_0_z_xxxyyz_xzzz, g_0_z_xxxyyz_yyyy, g_0_z_xxxyyz_yyyz, g_0_z_xxxyyz_yyzz, g_0_z_xxxyyz_yzzz, g_0_z_xxxyyz_zzzz, g_0_z_xxyyz_xxxx, g_0_z_xxyyz_xxxxx, g_0_z_xxyyz_xxxxy, g_0_z_xxyyz_xxxxz, g_0_z_xxyyz_xxxy, g_0_z_xxyyz_xxxyy, g_0_z_xxyyz_xxxyz, g_0_z_xxyyz_xxxz, g_0_z_xxyyz_xxxzz, g_0_z_xxyyz_xxyy, g_0_z_xxyyz_xxyyy, g_0_z_xxyyz_xxyyz, g_0_z_xxyyz_xxyz, g_0_z_xxyyz_xxyzz, g_0_z_xxyyz_xxzz, g_0_z_xxyyz_xxzzz, g_0_z_xxyyz_xyyy, g_0_z_xxyyz_xyyyy, g_0_z_xxyyz_xyyyz, g_0_z_xxyyz_xyyz, g_0_z_xxyyz_xyyzz, g_0_z_xxyyz_xyzz, g_0_z_xxyyz_xyzzz, g_0_z_xxyyz_xzzz, g_0_z_xxyyz_xzzzz, g_0_z_xxyyz_yyyy, g_0_z_xxyyz_yyyz, g_0_z_xxyyz_yyzz, g_0_z_xxyyz_yzzz, g_0_z_xxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyz_xxxx[k] = -g_0_z_xxyyz_xxxx[k] * ab_x + g_0_z_xxyyz_xxxxx[k];

                g_0_z_xxxyyz_xxxy[k] = -g_0_z_xxyyz_xxxy[k] * ab_x + g_0_z_xxyyz_xxxxy[k];

                g_0_z_xxxyyz_xxxz[k] = -g_0_z_xxyyz_xxxz[k] * ab_x + g_0_z_xxyyz_xxxxz[k];

                g_0_z_xxxyyz_xxyy[k] = -g_0_z_xxyyz_xxyy[k] * ab_x + g_0_z_xxyyz_xxxyy[k];

                g_0_z_xxxyyz_xxyz[k] = -g_0_z_xxyyz_xxyz[k] * ab_x + g_0_z_xxyyz_xxxyz[k];

                g_0_z_xxxyyz_xxzz[k] = -g_0_z_xxyyz_xxzz[k] * ab_x + g_0_z_xxyyz_xxxzz[k];

                g_0_z_xxxyyz_xyyy[k] = -g_0_z_xxyyz_xyyy[k] * ab_x + g_0_z_xxyyz_xxyyy[k];

                g_0_z_xxxyyz_xyyz[k] = -g_0_z_xxyyz_xyyz[k] * ab_x + g_0_z_xxyyz_xxyyz[k];

                g_0_z_xxxyyz_xyzz[k] = -g_0_z_xxyyz_xyzz[k] * ab_x + g_0_z_xxyyz_xxyzz[k];

                g_0_z_xxxyyz_xzzz[k] = -g_0_z_xxyyz_xzzz[k] * ab_x + g_0_z_xxyyz_xxzzz[k];

                g_0_z_xxxyyz_yyyy[k] = -g_0_z_xxyyz_yyyy[k] * ab_x + g_0_z_xxyyz_xyyyy[k];

                g_0_z_xxxyyz_yyyz[k] = -g_0_z_xxyyz_yyyz[k] * ab_x + g_0_z_xxyyz_xyyyz[k];

                g_0_z_xxxyyz_yyzz[k] = -g_0_z_xxyyz_yyzz[k] * ab_x + g_0_z_xxyyz_xyyzz[k];

                g_0_z_xxxyyz_yzzz[k] = -g_0_z_xxyyz_yzzz[k] * ab_x + g_0_z_xxyyz_xyzzz[k];

                g_0_z_xxxyyz_zzzz[k] = -g_0_z_xxyyz_zzzz[k] * ab_x + g_0_z_xxyyz_xzzzz[k];
            }

            /// Set up 960-975 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyzz_xxxx = cbuffer.data(ig_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxxy = cbuffer.data(ig_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxxz = cbuffer.data(ig_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxyy = cbuffer.data(ig_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxyz = cbuffer.data(ig_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxzz = cbuffer.data(ig_geom_01_off + 965 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyyy = cbuffer.data(ig_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyyz = cbuffer.data(ig_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyzz = cbuffer.data(ig_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xzzz = cbuffer.data(ig_geom_01_off + 969 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyyy = cbuffer.data(ig_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyyz = cbuffer.data(ig_geom_01_off + 971 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyzz = cbuffer.data(ig_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yzzz = cbuffer.data(ig_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_z_xxxyzz_zzzz = cbuffer.data(ig_geom_01_off + 974 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyzz_xxxx, g_0_z_xxxyzz_xxxy, g_0_z_xxxyzz_xxxz, g_0_z_xxxyzz_xxyy, g_0_z_xxxyzz_xxyz, g_0_z_xxxyzz_xxzz, g_0_z_xxxyzz_xyyy, g_0_z_xxxyzz_xyyz, g_0_z_xxxyzz_xyzz, g_0_z_xxxyzz_xzzz, g_0_z_xxxyzz_yyyy, g_0_z_xxxyzz_yyyz, g_0_z_xxxyzz_yyzz, g_0_z_xxxyzz_yzzz, g_0_z_xxxyzz_zzzz, g_0_z_xxyzz_xxxx, g_0_z_xxyzz_xxxxx, g_0_z_xxyzz_xxxxy, g_0_z_xxyzz_xxxxz, g_0_z_xxyzz_xxxy, g_0_z_xxyzz_xxxyy, g_0_z_xxyzz_xxxyz, g_0_z_xxyzz_xxxz, g_0_z_xxyzz_xxxzz, g_0_z_xxyzz_xxyy, g_0_z_xxyzz_xxyyy, g_0_z_xxyzz_xxyyz, g_0_z_xxyzz_xxyz, g_0_z_xxyzz_xxyzz, g_0_z_xxyzz_xxzz, g_0_z_xxyzz_xxzzz, g_0_z_xxyzz_xyyy, g_0_z_xxyzz_xyyyy, g_0_z_xxyzz_xyyyz, g_0_z_xxyzz_xyyz, g_0_z_xxyzz_xyyzz, g_0_z_xxyzz_xyzz, g_0_z_xxyzz_xyzzz, g_0_z_xxyzz_xzzz, g_0_z_xxyzz_xzzzz, g_0_z_xxyzz_yyyy, g_0_z_xxyzz_yyyz, g_0_z_xxyzz_yyzz, g_0_z_xxyzz_yzzz, g_0_z_xxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyzz_xxxx[k] = -g_0_z_xxyzz_xxxx[k] * ab_x + g_0_z_xxyzz_xxxxx[k];

                g_0_z_xxxyzz_xxxy[k] = -g_0_z_xxyzz_xxxy[k] * ab_x + g_0_z_xxyzz_xxxxy[k];

                g_0_z_xxxyzz_xxxz[k] = -g_0_z_xxyzz_xxxz[k] * ab_x + g_0_z_xxyzz_xxxxz[k];

                g_0_z_xxxyzz_xxyy[k] = -g_0_z_xxyzz_xxyy[k] * ab_x + g_0_z_xxyzz_xxxyy[k];

                g_0_z_xxxyzz_xxyz[k] = -g_0_z_xxyzz_xxyz[k] * ab_x + g_0_z_xxyzz_xxxyz[k];

                g_0_z_xxxyzz_xxzz[k] = -g_0_z_xxyzz_xxzz[k] * ab_x + g_0_z_xxyzz_xxxzz[k];

                g_0_z_xxxyzz_xyyy[k] = -g_0_z_xxyzz_xyyy[k] * ab_x + g_0_z_xxyzz_xxyyy[k];

                g_0_z_xxxyzz_xyyz[k] = -g_0_z_xxyzz_xyyz[k] * ab_x + g_0_z_xxyzz_xxyyz[k];

                g_0_z_xxxyzz_xyzz[k] = -g_0_z_xxyzz_xyzz[k] * ab_x + g_0_z_xxyzz_xxyzz[k];

                g_0_z_xxxyzz_xzzz[k] = -g_0_z_xxyzz_xzzz[k] * ab_x + g_0_z_xxyzz_xxzzz[k];

                g_0_z_xxxyzz_yyyy[k] = -g_0_z_xxyzz_yyyy[k] * ab_x + g_0_z_xxyzz_xyyyy[k];

                g_0_z_xxxyzz_yyyz[k] = -g_0_z_xxyzz_yyyz[k] * ab_x + g_0_z_xxyzz_xyyyz[k];

                g_0_z_xxxyzz_yyzz[k] = -g_0_z_xxyzz_yyzz[k] * ab_x + g_0_z_xxyzz_xyyzz[k];

                g_0_z_xxxyzz_yzzz[k] = -g_0_z_xxyzz_yzzz[k] * ab_x + g_0_z_xxyzz_xyzzz[k];

                g_0_z_xxxyzz_zzzz[k] = -g_0_z_xxyzz_zzzz[k] * ab_x + g_0_z_xxyzz_xzzzz[k];
            }

            /// Set up 975-990 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzzz_xxxx = cbuffer.data(ig_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxxy = cbuffer.data(ig_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxxz = cbuffer.data(ig_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxyy = cbuffer.data(ig_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxyz = cbuffer.data(ig_geom_01_off + 979 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxzz = cbuffer.data(ig_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyyy = cbuffer.data(ig_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyyz = cbuffer.data(ig_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyzz = cbuffer.data(ig_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xzzz = cbuffer.data(ig_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyyy = cbuffer.data(ig_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyyz = cbuffer.data(ig_geom_01_off + 986 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyzz = cbuffer.data(ig_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yzzz = cbuffer.data(ig_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_z_xxxzzz_zzzz = cbuffer.data(ig_geom_01_off + 989 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzzz_xxxx, g_0_z_xxxzzz_xxxy, g_0_z_xxxzzz_xxxz, g_0_z_xxxzzz_xxyy, g_0_z_xxxzzz_xxyz, g_0_z_xxxzzz_xxzz, g_0_z_xxxzzz_xyyy, g_0_z_xxxzzz_xyyz, g_0_z_xxxzzz_xyzz, g_0_z_xxxzzz_xzzz, g_0_z_xxxzzz_yyyy, g_0_z_xxxzzz_yyyz, g_0_z_xxxzzz_yyzz, g_0_z_xxxzzz_yzzz, g_0_z_xxxzzz_zzzz, g_0_z_xxzzz_xxxx, g_0_z_xxzzz_xxxxx, g_0_z_xxzzz_xxxxy, g_0_z_xxzzz_xxxxz, g_0_z_xxzzz_xxxy, g_0_z_xxzzz_xxxyy, g_0_z_xxzzz_xxxyz, g_0_z_xxzzz_xxxz, g_0_z_xxzzz_xxxzz, g_0_z_xxzzz_xxyy, g_0_z_xxzzz_xxyyy, g_0_z_xxzzz_xxyyz, g_0_z_xxzzz_xxyz, g_0_z_xxzzz_xxyzz, g_0_z_xxzzz_xxzz, g_0_z_xxzzz_xxzzz, g_0_z_xxzzz_xyyy, g_0_z_xxzzz_xyyyy, g_0_z_xxzzz_xyyyz, g_0_z_xxzzz_xyyz, g_0_z_xxzzz_xyyzz, g_0_z_xxzzz_xyzz, g_0_z_xxzzz_xyzzz, g_0_z_xxzzz_xzzz, g_0_z_xxzzz_xzzzz, g_0_z_xxzzz_yyyy, g_0_z_xxzzz_yyyz, g_0_z_xxzzz_yyzz, g_0_z_xxzzz_yzzz, g_0_z_xxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzzz_xxxx[k] = -g_0_z_xxzzz_xxxx[k] * ab_x + g_0_z_xxzzz_xxxxx[k];

                g_0_z_xxxzzz_xxxy[k] = -g_0_z_xxzzz_xxxy[k] * ab_x + g_0_z_xxzzz_xxxxy[k];

                g_0_z_xxxzzz_xxxz[k] = -g_0_z_xxzzz_xxxz[k] * ab_x + g_0_z_xxzzz_xxxxz[k];

                g_0_z_xxxzzz_xxyy[k] = -g_0_z_xxzzz_xxyy[k] * ab_x + g_0_z_xxzzz_xxxyy[k];

                g_0_z_xxxzzz_xxyz[k] = -g_0_z_xxzzz_xxyz[k] * ab_x + g_0_z_xxzzz_xxxyz[k];

                g_0_z_xxxzzz_xxzz[k] = -g_0_z_xxzzz_xxzz[k] * ab_x + g_0_z_xxzzz_xxxzz[k];

                g_0_z_xxxzzz_xyyy[k] = -g_0_z_xxzzz_xyyy[k] * ab_x + g_0_z_xxzzz_xxyyy[k];

                g_0_z_xxxzzz_xyyz[k] = -g_0_z_xxzzz_xyyz[k] * ab_x + g_0_z_xxzzz_xxyyz[k];

                g_0_z_xxxzzz_xyzz[k] = -g_0_z_xxzzz_xyzz[k] * ab_x + g_0_z_xxzzz_xxyzz[k];

                g_0_z_xxxzzz_xzzz[k] = -g_0_z_xxzzz_xzzz[k] * ab_x + g_0_z_xxzzz_xxzzz[k];

                g_0_z_xxxzzz_yyyy[k] = -g_0_z_xxzzz_yyyy[k] * ab_x + g_0_z_xxzzz_xyyyy[k];

                g_0_z_xxxzzz_yyyz[k] = -g_0_z_xxzzz_yyyz[k] * ab_x + g_0_z_xxzzz_xyyyz[k];

                g_0_z_xxxzzz_yyzz[k] = -g_0_z_xxzzz_yyzz[k] * ab_x + g_0_z_xxzzz_xyyzz[k];

                g_0_z_xxxzzz_yzzz[k] = -g_0_z_xxzzz_yzzz[k] * ab_x + g_0_z_xxzzz_xyzzz[k];

                g_0_z_xxxzzz_zzzz[k] = -g_0_z_xxzzz_zzzz[k] * ab_x + g_0_z_xxzzz_xzzzz[k];
            }

            /// Set up 990-1005 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyy_xxxx = cbuffer.data(ig_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxxy = cbuffer.data(ig_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxxz = cbuffer.data(ig_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxyy = cbuffer.data(ig_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxyz = cbuffer.data(ig_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxzz = cbuffer.data(ig_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyyy = cbuffer.data(ig_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyyz = cbuffer.data(ig_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyzz = cbuffer.data(ig_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xzzz = cbuffer.data(ig_geom_01_off + 999 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyyy = cbuffer.data(ig_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyyz = cbuffer.data(ig_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyzz = cbuffer.data(ig_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yzzz = cbuffer.data(ig_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_z_xxyyyy_zzzz = cbuffer.data(ig_geom_01_off + 1004 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyy_xxxx, g_0_z_xxyyyy_xxxy, g_0_z_xxyyyy_xxxz, g_0_z_xxyyyy_xxyy, g_0_z_xxyyyy_xxyz, g_0_z_xxyyyy_xxzz, g_0_z_xxyyyy_xyyy, g_0_z_xxyyyy_xyyz, g_0_z_xxyyyy_xyzz, g_0_z_xxyyyy_xzzz, g_0_z_xxyyyy_yyyy, g_0_z_xxyyyy_yyyz, g_0_z_xxyyyy_yyzz, g_0_z_xxyyyy_yzzz, g_0_z_xxyyyy_zzzz, g_0_z_xyyyy_xxxx, g_0_z_xyyyy_xxxxx, g_0_z_xyyyy_xxxxy, g_0_z_xyyyy_xxxxz, g_0_z_xyyyy_xxxy, g_0_z_xyyyy_xxxyy, g_0_z_xyyyy_xxxyz, g_0_z_xyyyy_xxxz, g_0_z_xyyyy_xxxzz, g_0_z_xyyyy_xxyy, g_0_z_xyyyy_xxyyy, g_0_z_xyyyy_xxyyz, g_0_z_xyyyy_xxyz, g_0_z_xyyyy_xxyzz, g_0_z_xyyyy_xxzz, g_0_z_xyyyy_xxzzz, g_0_z_xyyyy_xyyy, g_0_z_xyyyy_xyyyy, g_0_z_xyyyy_xyyyz, g_0_z_xyyyy_xyyz, g_0_z_xyyyy_xyyzz, g_0_z_xyyyy_xyzz, g_0_z_xyyyy_xyzzz, g_0_z_xyyyy_xzzz, g_0_z_xyyyy_xzzzz, g_0_z_xyyyy_yyyy, g_0_z_xyyyy_yyyz, g_0_z_xyyyy_yyzz, g_0_z_xyyyy_yzzz, g_0_z_xyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyy_xxxx[k] = -g_0_z_xyyyy_xxxx[k] * ab_x + g_0_z_xyyyy_xxxxx[k];

                g_0_z_xxyyyy_xxxy[k] = -g_0_z_xyyyy_xxxy[k] * ab_x + g_0_z_xyyyy_xxxxy[k];

                g_0_z_xxyyyy_xxxz[k] = -g_0_z_xyyyy_xxxz[k] * ab_x + g_0_z_xyyyy_xxxxz[k];

                g_0_z_xxyyyy_xxyy[k] = -g_0_z_xyyyy_xxyy[k] * ab_x + g_0_z_xyyyy_xxxyy[k];

                g_0_z_xxyyyy_xxyz[k] = -g_0_z_xyyyy_xxyz[k] * ab_x + g_0_z_xyyyy_xxxyz[k];

                g_0_z_xxyyyy_xxzz[k] = -g_0_z_xyyyy_xxzz[k] * ab_x + g_0_z_xyyyy_xxxzz[k];

                g_0_z_xxyyyy_xyyy[k] = -g_0_z_xyyyy_xyyy[k] * ab_x + g_0_z_xyyyy_xxyyy[k];

                g_0_z_xxyyyy_xyyz[k] = -g_0_z_xyyyy_xyyz[k] * ab_x + g_0_z_xyyyy_xxyyz[k];

                g_0_z_xxyyyy_xyzz[k] = -g_0_z_xyyyy_xyzz[k] * ab_x + g_0_z_xyyyy_xxyzz[k];

                g_0_z_xxyyyy_xzzz[k] = -g_0_z_xyyyy_xzzz[k] * ab_x + g_0_z_xyyyy_xxzzz[k];

                g_0_z_xxyyyy_yyyy[k] = -g_0_z_xyyyy_yyyy[k] * ab_x + g_0_z_xyyyy_xyyyy[k];

                g_0_z_xxyyyy_yyyz[k] = -g_0_z_xyyyy_yyyz[k] * ab_x + g_0_z_xyyyy_xyyyz[k];

                g_0_z_xxyyyy_yyzz[k] = -g_0_z_xyyyy_yyzz[k] * ab_x + g_0_z_xyyyy_xyyzz[k];

                g_0_z_xxyyyy_yzzz[k] = -g_0_z_xyyyy_yzzz[k] * ab_x + g_0_z_xyyyy_xyzzz[k];

                g_0_z_xxyyyy_zzzz[k] = -g_0_z_xyyyy_zzzz[k] * ab_x + g_0_z_xyyyy_xzzzz[k];
            }

            /// Set up 1005-1020 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyz_xxxx = cbuffer.data(ig_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxxy = cbuffer.data(ig_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxxz = cbuffer.data(ig_geom_01_off + 1007 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxyy = cbuffer.data(ig_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxyz = cbuffer.data(ig_geom_01_off + 1009 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxzz = cbuffer.data(ig_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyyy = cbuffer.data(ig_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyyz = cbuffer.data(ig_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyzz = cbuffer.data(ig_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xzzz = cbuffer.data(ig_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyyy = cbuffer.data(ig_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyyz = cbuffer.data(ig_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyzz = cbuffer.data(ig_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yzzz = cbuffer.data(ig_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_z_xxyyyz_zzzz = cbuffer.data(ig_geom_01_off + 1019 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyz_xxxx, g_0_z_xxyyyz_xxxy, g_0_z_xxyyyz_xxxz, g_0_z_xxyyyz_xxyy, g_0_z_xxyyyz_xxyz, g_0_z_xxyyyz_xxzz, g_0_z_xxyyyz_xyyy, g_0_z_xxyyyz_xyyz, g_0_z_xxyyyz_xyzz, g_0_z_xxyyyz_xzzz, g_0_z_xxyyyz_yyyy, g_0_z_xxyyyz_yyyz, g_0_z_xxyyyz_yyzz, g_0_z_xxyyyz_yzzz, g_0_z_xxyyyz_zzzz, g_0_z_xyyyz_xxxx, g_0_z_xyyyz_xxxxx, g_0_z_xyyyz_xxxxy, g_0_z_xyyyz_xxxxz, g_0_z_xyyyz_xxxy, g_0_z_xyyyz_xxxyy, g_0_z_xyyyz_xxxyz, g_0_z_xyyyz_xxxz, g_0_z_xyyyz_xxxzz, g_0_z_xyyyz_xxyy, g_0_z_xyyyz_xxyyy, g_0_z_xyyyz_xxyyz, g_0_z_xyyyz_xxyz, g_0_z_xyyyz_xxyzz, g_0_z_xyyyz_xxzz, g_0_z_xyyyz_xxzzz, g_0_z_xyyyz_xyyy, g_0_z_xyyyz_xyyyy, g_0_z_xyyyz_xyyyz, g_0_z_xyyyz_xyyz, g_0_z_xyyyz_xyyzz, g_0_z_xyyyz_xyzz, g_0_z_xyyyz_xyzzz, g_0_z_xyyyz_xzzz, g_0_z_xyyyz_xzzzz, g_0_z_xyyyz_yyyy, g_0_z_xyyyz_yyyz, g_0_z_xyyyz_yyzz, g_0_z_xyyyz_yzzz, g_0_z_xyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyz_xxxx[k] = -g_0_z_xyyyz_xxxx[k] * ab_x + g_0_z_xyyyz_xxxxx[k];

                g_0_z_xxyyyz_xxxy[k] = -g_0_z_xyyyz_xxxy[k] * ab_x + g_0_z_xyyyz_xxxxy[k];

                g_0_z_xxyyyz_xxxz[k] = -g_0_z_xyyyz_xxxz[k] * ab_x + g_0_z_xyyyz_xxxxz[k];

                g_0_z_xxyyyz_xxyy[k] = -g_0_z_xyyyz_xxyy[k] * ab_x + g_0_z_xyyyz_xxxyy[k];

                g_0_z_xxyyyz_xxyz[k] = -g_0_z_xyyyz_xxyz[k] * ab_x + g_0_z_xyyyz_xxxyz[k];

                g_0_z_xxyyyz_xxzz[k] = -g_0_z_xyyyz_xxzz[k] * ab_x + g_0_z_xyyyz_xxxzz[k];

                g_0_z_xxyyyz_xyyy[k] = -g_0_z_xyyyz_xyyy[k] * ab_x + g_0_z_xyyyz_xxyyy[k];

                g_0_z_xxyyyz_xyyz[k] = -g_0_z_xyyyz_xyyz[k] * ab_x + g_0_z_xyyyz_xxyyz[k];

                g_0_z_xxyyyz_xyzz[k] = -g_0_z_xyyyz_xyzz[k] * ab_x + g_0_z_xyyyz_xxyzz[k];

                g_0_z_xxyyyz_xzzz[k] = -g_0_z_xyyyz_xzzz[k] * ab_x + g_0_z_xyyyz_xxzzz[k];

                g_0_z_xxyyyz_yyyy[k] = -g_0_z_xyyyz_yyyy[k] * ab_x + g_0_z_xyyyz_xyyyy[k];

                g_0_z_xxyyyz_yyyz[k] = -g_0_z_xyyyz_yyyz[k] * ab_x + g_0_z_xyyyz_xyyyz[k];

                g_0_z_xxyyyz_yyzz[k] = -g_0_z_xyyyz_yyzz[k] * ab_x + g_0_z_xyyyz_xyyzz[k];

                g_0_z_xxyyyz_yzzz[k] = -g_0_z_xyyyz_yzzz[k] * ab_x + g_0_z_xyyyz_xyzzz[k];

                g_0_z_xxyyyz_zzzz[k] = -g_0_z_xyyyz_zzzz[k] * ab_x + g_0_z_xyyyz_xzzzz[k];
            }

            /// Set up 1020-1035 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyzz_xxxx = cbuffer.data(ig_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxxy = cbuffer.data(ig_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxxz = cbuffer.data(ig_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxyy = cbuffer.data(ig_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxyz = cbuffer.data(ig_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxzz = cbuffer.data(ig_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyyy = cbuffer.data(ig_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyyz = cbuffer.data(ig_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyzz = cbuffer.data(ig_geom_01_off + 1028 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xzzz = cbuffer.data(ig_geom_01_off + 1029 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyyy = cbuffer.data(ig_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyyz = cbuffer.data(ig_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyzz = cbuffer.data(ig_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yzzz = cbuffer.data(ig_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_z_xxyyzz_zzzz = cbuffer.data(ig_geom_01_off + 1034 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyzz_xxxx, g_0_z_xxyyzz_xxxy, g_0_z_xxyyzz_xxxz, g_0_z_xxyyzz_xxyy, g_0_z_xxyyzz_xxyz, g_0_z_xxyyzz_xxzz, g_0_z_xxyyzz_xyyy, g_0_z_xxyyzz_xyyz, g_0_z_xxyyzz_xyzz, g_0_z_xxyyzz_xzzz, g_0_z_xxyyzz_yyyy, g_0_z_xxyyzz_yyyz, g_0_z_xxyyzz_yyzz, g_0_z_xxyyzz_yzzz, g_0_z_xxyyzz_zzzz, g_0_z_xyyzz_xxxx, g_0_z_xyyzz_xxxxx, g_0_z_xyyzz_xxxxy, g_0_z_xyyzz_xxxxz, g_0_z_xyyzz_xxxy, g_0_z_xyyzz_xxxyy, g_0_z_xyyzz_xxxyz, g_0_z_xyyzz_xxxz, g_0_z_xyyzz_xxxzz, g_0_z_xyyzz_xxyy, g_0_z_xyyzz_xxyyy, g_0_z_xyyzz_xxyyz, g_0_z_xyyzz_xxyz, g_0_z_xyyzz_xxyzz, g_0_z_xyyzz_xxzz, g_0_z_xyyzz_xxzzz, g_0_z_xyyzz_xyyy, g_0_z_xyyzz_xyyyy, g_0_z_xyyzz_xyyyz, g_0_z_xyyzz_xyyz, g_0_z_xyyzz_xyyzz, g_0_z_xyyzz_xyzz, g_0_z_xyyzz_xyzzz, g_0_z_xyyzz_xzzz, g_0_z_xyyzz_xzzzz, g_0_z_xyyzz_yyyy, g_0_z_xyyzz_yyyz, g_0_z_xyyzz_yyzz, g_0_z_xyyzz_yzzz, g_0_z_xyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyzz_xxxx[k] = -g_0_z_xyyzz_xxxx[k] * ab_x + g_0_z_xyyzz_xxxxx[k];

                g_0_z_xxyyzz_xxxy[k] = -g_0_z_xyyzz_xxxy[k] * ab_x + g_0_z_xyyzz_xxxxy[k];

                g_0_z_xxyyzz_xxxz[k] = -g_0_z_xyyzz_xxxz[k] * ab_x + g_0_z_xyyzz_xxxxz[k];

                g_0_z_xxyyzz_xxyy[k] = -g_0_z_xyyzz_xxyy[k] * ab_x + g_0_z_xyyzz_xxxyy[k];

                g_0_z_xxyyzz_xxyz[k] = -g_0_z_xyyzz_xxyz[k] * ab_x + g_0_z_xyyzz_xxxyz[k];

                g_0_z_xxyyzz_xxzz[k] = -g_0_z_xyyzz_xxzz[k] * ab_x + g_0_z_xyyzz_xxxzz[k];

                g_0_z_xxyyzz_xyyy[k] = -g_0_z_xyyzz_xyyy[k] * ab_x + g_0_z_xyyzz_xxyyy[k];

                g_0_z_xxyyzz_xyyz[k] = -g_0_z_xyyzz_xyyz[k] * ab_x + g_0_z_xyyzz_xxyyz[k];

                g_0_z_xxyyzz_xyzz[k] = -g_0_z_xyyzz_xyzz[k] * ab_x + g_0_z_xyyzz_xxyzz[k];

                g_0_z_xxyyzz_xzzz[k] = -g_0_z_xyyzz_xzzz[k] * ab_x + g_0_z_xyyzz_xxzzz[k];

                g_0_z_xxyyzz_yyyy[k] = -g_0_z_xyyzz_yyyy[k] * ab_x + g_0_z_xyyzz_xyyyy[k];

                g_0_z_xxyyzz_yyyz[k] = -g_0_z_xyyzz_yyyz[k] * ab_x + g_0_z_xyyzz_xyyyz[k];

                g_0_z_xxyyzz_yyzz[k] = -g_0_z_xyyzz_yyzz[k] * ab_x + g_0_z_xyyzz_xyyzz[k];

                g_0_z_xxyyzz_yzzz[k] = -g_0_z_xyyzz_yzzz[k] * ab_x + g_0_z_xyyzz_xyzzz[k];

                g_0_z_xxyyzz_zzzz[k] = -g_0_z_xyyzz_zzzz[k] * ab_x + g_0_z_xyyzz_xzzzz[k];
            }

            /// Set up 1035-1050 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzzz_xxxx = cbuffer.data(ig_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxxy = cbuffer.data(ig_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxxz = cbuffer.data(ig_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxyy = cbuffer.data(ig_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxyz = cbuffer.data(ig_geom_01_off + 1039 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxzz = cbuffer.data(ig_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyyy = cbuffer.data(ig_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyyz = cbuffer.data(ig_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyzz = cbuffer.data(ig_geom_01_off + 1043 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xzzz = cbuffer.data(ig_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyyy = cbuffer.data(ig_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyyz = cbuffer.data(ig_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyzz = cbuffer.data(ig_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yzzz = cbuffer.data(ig_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_z_xxyzzz_zzzz = cbuffer.data(ig_geom_01_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzzz_xxxx, g_0_z_xxyzzz_xxxy, g_0_z_xxyzzz_xxxz, g_0_z_xxyzzz_xxyy, g_0_z_xxyzzz_xxyz, g_0_z_xxyzzz_xxzz, g_0_z_xxyzzz_xyyy, g_0_z_xxyzzz_xyyz, g_0_z_xxyzzz_xyzz, g_0_z_xxyzzz_xzzz, g_0_z_xxyzzz_yyyy, g_0_z_xxyzzz_yyyz, g_0_z_xxyzzz_yyzz, g_0_z_xxyzzz_yzzz, g_0_z_xxyzzz_zzzz, g_0_z_xyzzz_xxxx, g_0_z_xyzzz_xxxxx, g_0_z_xyzzz_xxxxy, g_0_z_xyzzz_xxxxz, g_0_z_xyzzz_xxxy, g_0_z_xyzzz_xxxyy, g_0_z_xyzzz_xxxyz, g_0_z_xyzzz_xxxz, g_0_z_xyzzz_xxxzz, g_0_z_xyzzz_xxyy, g_0_z_xyzzz_xxyyy, g_0_z_xyzzz_xxyyz, g_0_z_xyzzz_xxyz, g_0_z_xyzzz_xxyzz, g_0_z_xyzzz_xxzz, g_0_z_xyzzz_xxzzz, g_0_z_xyzzz_xyyy, g_0_z_xyzzz_xyyyy, g_0_z_xyzzz_xyyyz, g_0_z_xyzzz_xyyz, g_0_z_xyzzz_xyyzz, g_0_z_xyzzz_xyzz, g_0_z_xyzzz_xyzzz, g_0_z_xyzzz_xzzz, g_0_z_xyzzz_xzzzz, g_0_z_xyzzz_yyyy, g_0_z_xyzzz_yyyz, g_0_z_xyzzz_yyzz, g_0_z_xyzzz_yzzz, g_0_z_xyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzzz_xxxx[k] = -g_0_z_xyzzz_xxxx[k] * ab_x + g_0_z_xyzzz_xxxxx[k];

                g_0_z_xxyzzz_xxxy[k] = -g_0_z_xyzzz_xxxy[k] * ab_x + g_0_z_xyzzz_xxxxy[k];

                g_0_z_xxyzzz_xxxz[k] = -g_0_z_xyzzz_xxxz[k] * ab_x + g_0_z_xyzzz_xxxxz[k];

                g_0_z_xxyzzz_xxyy[k] = -g_0_z_xyzzz_xxyy[k] * ab_x + g_0_z_xyzzz_xxxyy[k];

                g_0_z_xxyzzz_xxyz[k] = -g_0_z_xyzzz_xxyz[k] * ab_x + g_0_z_xyzzz_xxxyz[k];

                g_0_z_xxyzzz_xxzz[k] = -g_0_z_xyzzz_xxzz[k] * ab_x + g_0_z_xyzzz_xxxzz[k];

                g_0_z_xxyzzz_xyyy[k] = -g_0_z_xyzzz_xyyy[k] * ab_x + g_0_z_xyzzz_xxyyy[k];

                g_0_z_xxyzzz_xyyz[k] = -g_0_z_xyzzz_xyyz[k] * ab_x + g_0_z_xyzzz_xxyyz[k];

                g_0_z_xxyzzz_xyzz[k] = -g_0_z_xyzzz_xyzz[k] * ab_x + g_0_z_xyzzz_xxyzz[k];

                g_0_z_xxyzzz_xzzz[k] = -g_0_z_xyzzz_xzzz[k] * ab_x + g_0_z_xyzzz_xxzzz[k];

                g_0_z_xxyzzz_yyyy[k] = -g_0_z_xyzzz_yyyy[k] * ab_x + g_0_z_xyzzz_xyyyy[k];

                g_0_z_xxyzzz_yyyz[k] = -g_0_z_xyzzz_yyyz[k] * ab_x + g_0_z_xyzzz_xyyyz[k];

                g_0_z_xxyzzz_yyzz[k] = -g_0_z_xyzzz_yyzz[k] * ab_x + g_0_z_xyzzz_xyyzz[k];

                g_0_z_xxyzzz_yzzz[k] = -g_0_z_xyzzz_yzzz[k] * ab_x + g_0_z_xyzzz_xyzzz[k];

                g_0_z_xxyzzz_zzzz[k] = -g_0_z_xyzzz_zzzz[k] * ab_x + g_0_z_xyzzz_xzzzz[k];
            }

            /// Set up 1050-1065 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzzz_xxxx = cbuffer.data(ig_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxxy = cbuffer.data(ig_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxxz = cbuffer.data(ig_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxyy = cbuffer.data(ig_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxyz = cbuffer.data(ig_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxzz = cbuffer.data(ig_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyyy = cbuffer.data(ig_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyyz = cbuffer.data(ig_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyzz = cbuffer.data(ig_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xzzz = cbuffer.data(ig_geom_01_off + 1059 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyyy = cbuffer.data(ig_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyyz = cbuffer.data(ig_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyzz = cbuffer.data(ig_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yzzz = cbuffer.data(ig_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_z_xxzzzz_zzzz = cbuffer.data(ig_geom_01_off + 1064 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzzz_xxxx, g_0_z_xxzzzz_xxxy, g_0_z_xxzzzz_xxxz, g_0_z_xxzzzz_xxyy, g_0_z_xxzzzz_xxyz, g_0_z_xxzzzz_xxzz, g_0_z_xxzzzz_xyyy, g_0_z_xxzzzz_xyyz, g_0_z_xxzzzz_xyzz, g_0_z_xxzzzz_xzzz, g_0_z_xxzzzz_yyyy, g_0_z_xxzzzz_yyyz, g_0_z_xxzzzz_yyzz, g_0_z_xxzzzz_yzzz, g_0_z_xxzzzz_zzzz, g_0_z_xzzzz_xxxx, g_0_z_xzzzz_xxxxx, g_0_z_xzzzz_xxxxy, g_0_z_xzzzz_xxxxz, g_0_z_xzzzz_xxxy, g_0_z_xzzzz_xxxyy, g_0_z_xzzzz_xxxyz, g_0_z_xzzzz_xxxz, g_0_z_xzzzz_xxxzz, g_0_z_xzzzz_xxyy, g_0_z_xzzzz_xxyyy, g_0_z_xzzzz_xxyyz, g_0_z_xzzzz_xxyz, g_0_z_xzzzz_xxyzz, g_0_z_xzzzz_xxzz, g_0_z_xzzzz_xxzzz, g_0_z_xzzzz_xyyy, g_0_z_xzzzz_xyyyy, g_0_z_xzzzz_xyyyz, g_0_z_xzzzz_xyyz, g_0_z_xzzzz_xyyzz, g_0_z_xzzzz_xyzz, g_0_z_xzzzz_xyzzz, g_0_z_xzzzz_xzzz, g_0_z_xzzzz_xzzzz, g_0_z_xzzzz_yyyy, g_0_z_xzzzz_yyyz, g_0_z_xzzzz_yyzz, g_0_z_xzzzz_yzzz, g_0_z_xzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzzz_xxxx[k] = -g_0_z_xzzzz_xxxx[k] * ab_x + g_0_z_xzzzz_xxxxx[k];

                g_0_z_xxzzzz_xxxy[k] = -g_0_z_xzzzz_xxxy[k] * ab_x + g_0_z_xzzzz_xxxxy[k];

                g_0_z_xxzzzz_xxxz[k] = -g_0_z_xzzzz_xxxz[k] * ab_x + g_0_z_xzzzz_xxxxz[k];

                g_0_z_xxzzzz_xxyy[k] = -g_0_z_xzzzz_xxyy[k] * ab_x + g_0_z_xzzzz_xxxyy[k];

                g_0_z_xxzzzz_xxyz[k] = -g_0_z_xzzzz_xxyz[k] * ab_x + g_0_z_xzzzz_xxxyz[k];

                g_0_z_xxzzzz_xxzz[k] = -g_0_z_xzzzz_xxzz[k] * ab_x + g_0_z_xzzzz_xxxzz[k];

                g_0_z_xxzzzz_xyyy[k] = -g_0_z_xzzzz_xyyy[k] * ab_x + g_0_z_xzzzz_xxyyy[k];

                g_0_z_xxzzzz_xyyz[k] = -g_0_z_xzzzz_xyyz[k] * ab_x + g_0_z_xzzzz_xxyyz[k];

                g_0_z_xxzzzz_xyzz[k] = -g_0_z_xzzzz_xyzz[k] * ab_x + g_0_z_xzzzz_xxyzz[k];

                g_0_z_xxzzzz_xzzz[k] = -g_0_z_xzzzz_xzzz[k] * ab_x + g_0_z_xzzzz_xxzzz[k];

                g_0_z_xxzzzz_yyyy[k] = -g_0_z_xzzzz_yyyy[k] * ab_x + g_0_z_xzzzz_xyyyy[k];

                g_0_z_xxzzzz_yyyz[k] = -g_0_z_xzzzz_yyyz[k] * ab_x + g_0_z_xzzzz_xyyyz[k];

                g_0_z_xxzzzz_yyzz[k] = -g_0_z_xzzzz_yyzz[k] * ab_x + g_0_z_xzzzz_xyyzz[k];

                g_0_z_xxzzzz_yzzz[k] = -g_0_z_xzzzz_yzzz[k] * ab_x + g_0_z_xzzzz_xyzzz[k];

                g_0_z_xxzzzz_zzzz[k] = -g_0_z_xzzzz_zzzz[k] * ab_x + g_0_z_xzzzz_xzzzz[k];
            }

            /// Set up 1065-1080 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyy_xxxx = cbuffer.data(ig_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxxy = cbuffer.data(ig_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxxz = cbuffer.data(ig_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxyy = cbuffer.data(ig_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxyz = cbuffer.data(ig_geom_01_off + 1069 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxzz = cbuffer.data(ig_geom_01_off + 1070 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyyy = cbuffer.data(ig_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyyz = cbuffer.data(ig_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyzz = cbuffer.data(ig_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xzzz = cbuffer.data(ig_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyyy = cbuffer.data(ig_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyyz = cbuffer.data(ig_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyzz = cbuffer.data(ig_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yzzz = cbuffer.data(ig_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_z_xyyyyy_zzzz = cbuffer.data(ig_geom_01_off + 1079 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyy_xxxx, g_0_z_xyyyyy_xxxy, g_0_z_xyyyyy_xxxz, g_0_z_xyyyyy_xxyy, g_0_z_xyyyyy_xxyz, g_0_z_xyyyyy_xxzz, g_0_z_xyyyyy_xyyy, g_0_z_xyyyyy_xyyz, g_0_z_xyyyyy_xyzz, g_0_z_xyyyyy_xzzz, g_0_z_xyyyyy_yyyy, g_0_z_xyyyyy_yyyz, g_0_z_xyyyyy_yyzz, g_0_z_xyyyyy_yzzz, g_0_z_xyyyyy_zzzz, g_0_z_yyyyy_xxxx, g_0_z_yyyyy_xxxxx, g_0_z_yyyyy_xxxxy, g_0_z_yyyyy_xxxxz, g_0_z_yyyyy_xxxy, g_0_z_yyyyy_xxxyy, g_0_z_yyyyy_xxxyz, g_0_z_yyyyy_xxxz, g_0_z_yyyyy_xxxzz, g_0_z_yyyyy_xxyy, g_0_z_yyyyy_xxyyy, g_0_z_yyyyy_xxyyz, g_0_z_yyyyy_xxyz, g_0_z_yyyyy_xxyzz, g_0_z_yyyyy_xxzz, g_0_z_yyyyy_xxzzz, g_0_z_yyyyy_xyyy, g_0_z_yyyyy_xyyyy, g_0_z_yyyyy_xyyyz, g_0_z_yyyyy_xyyz, g_0_z_yyyyy_xyyzz, g_0_z_yyyyy_xyzz, g_0_z_yyyyy_xyzzz, g_0_z_yyyyy_xzzz, g_0_z_yyyyy_xzzzz, g_0_z_yyyyy_yyyy, g_0_z_yyyyy_yyyz, g_0_z_yyyyy_yyzz, g_0_z_yyyyy_yzzz, g_0_z_yyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyy_xxxx[k] = -g_0_z_yyyyy_xxxx[k] * ab_x + g_0_z_yyyyy_xxxxx[k];

                g_0_z_xyyyyy_xxxy[k] = -g_0_z_yyyyy_xxxy[k] * ab_x + g_0_z_yyyyy_xxxxy[k];

                g_0_z_xyyyyy_xxxz[k] = -g_0_z_yyyyy_xxxz[k] * ab_x + g_0_z_yyyyy_xxxxz[k];

                g_0_z_xyyyyy_xxyy[k] = -g_0_z_yyyyy_xxyy[k] * ab_x + g_0_z_yyyyy_xxxyy[k];

                g_0_z_xyyyyy_xxyz[k] = -g_0_z_yyyyy_xxyz[k] * ab_x + g_0_z_yyyyy_xxxyz[k];

                g_0_z_xyyyyy_xxzz[k] = -g_0_z_yyyyy_xxzz[k] * ab_x + g_0_z_yyyyy_xxxzz[k];

                g_0_z_xyyyyy_xyyy[k] = -g_0_z_yyyyy_xyyy[k] * ab_x + g_0_z_yyyyy_xxyyy[k];

                g_0_z_xyyyyy_xyyz[k] = -g_0_z_yyyyy_xyyz[k] * ab_x + g_0_z_yyyyy_xxyyz[k];

                g_0_z_xyyyyy_xyzz[k] = -g_0_z_yyyyy_xyzz[k] * ab_x + g_0_z_yyyyy_xxyzz[k];

                g_0_z_xyyyyy_xzzz[k] = -g_0_z_yyyyy_xzzz[k] * ab_x + g_0_z_yyyyy_xxzzz[k];

                g_0_z_xyyyyy_yyyy[k] = -g_0_z_yyyyy_yyyy[k] * ab_x + g_0_z_yyyyy_xyyyy[k];

                g_0_z_xyyyyy_yyyz[k] = -g_0_z_yyyyy_yyyz[k] * ab_x + g_0_z_yyyyy_xyyyz[k];

                g_0_z_xyyyyy_yyzz[k] = -g_0_z_yyyyy_yyzz[k] * ab_x + g_0_z_yyyyy_xyyzz[k];

                g_0_z_xyyyyy_yzzz[k] = -g_0_z_yyyyy_yzzz[k] * ab_x + g_0_z_yyyyy_xyzzz[k];

                g_0_z_xyyyyy_zzzz[k] = -g_0_z_yyyyy_zzzz[k] * ab_x + g_0_z_yyyyy_xzzzz[k];
            }

            /// Set up 1080-1095 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyz_xxxx = cbuffer.data(ig_geom_01_off + 1080 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxxy = cbuffer.data(ig_geom_01_off + 1081 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxxz = cbuffer.data(ig_geom_01_off + 1082 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxyy = cbuffer.data(ig_geom_01_off + 1083 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxyz = cbuffer.data(ig_geom_01_off + 1084 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxzz = cbuffer.data(ig_geom_01_off + 1085 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyyy = cbuffer.data(ig_geom_01_off + 1086 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyyz = cbuffer.data(ig_geom_01_off + 1087 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyzz = cbuffer.data(ig_geom_01_off + 1088 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xzzz = cbuffer.data(ig_geom_01_off + 1089 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyyy = cbuffer.data(ig_geom_01_off + 1090 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyyz = cbuffer.data(ig_geom_01_off + 1091 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyzz = cbuffer.data(ig_geom_01_off + 1092 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yzzz = cbuffer.data(ig_geom_01_off + 1093 * ccomps * dcomps);

            auto g_0_z_xyyyyz_zzzz = cbuffer.data(ig_geom_01_off + 1094 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyz_xxxx, g_0_z_xyyyyz_xxxy, g_0_z_xyyyyz_xxxz, g_0_z_xyyyyz_xxyy, g_0_z_xyyyyz_xxyz, g_0_z_xyyyyz_xxzz, g_0_z_xyyyyz_xyyy, g_0_z_xyyyyz_xyyz, g_0_z_xyyyyz_xyzz, g_0_z_xyyyyz_xzzz, g_0_z_xyyyyz_yyyy, g_0_z_xyyyyz_yyyz, g_0_z_xyyyyz_yyzz, g_0_z_xyyyyz_yzzz, g_0_z_xyyyyz_zzzz, g_0_z_yyyyz_xxxx, g_0_z_yyyyz_xxxxx, g_0_z_yyyyz_xxxxy, g_0_z_yyyyz_xxxxz, g_0_z_yyyyz_xxxy, g_0_z_yyyyz_xxxyy, g_0_z_yyyyz_xxxyz, g_0_z_yyyyz_xxxz, g_0_z_yyyyz_xxxzz, g_0_z_yyyyz_xxyy, g_0_z_yyyyz_xxyyy, g_0_z_yyyyz_xxyyz, g_0_z_yyyyz_xxyz, g_0_z_yyyyz_xxyzz, g_0_z_yyyyz_xxzz, g_0_z_yyyyz_xxzzz, g_0_z_yyyyz_xyyy, g_0_z_yyyyz_xyyyy, g_0_z_yyyyz_xyyyz, g_0_z_yyyyz_xyyz, g_0_z_yyyyz_xyyzz, g_0_z_yyyyz_xyzz, g_0_z_yyyyz_xyzzz, g_0_z_yyyyz_xzzz, g_0_z_yyyyz_xzzzz, g_0_z_yyyyz_yyyy, g_0_z_yyyyz_yyyz, g_0_z_yyyyz_yyzz, g_0_z_yyyyz_yzzz, g_0_z_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyz_xxxx[k] = -g_0_z_yyyyz_xxxx[k] * ab_x + g_0_z_yyyyz_xxxxx[k];

                g_0_z_xyyyyz_xxxy[k] = -g_0_z_yyyyz_xxxy[k] * ab_x + g_0_z_yyyyz_xxxxy[k];

                g_0_z_xyyyyz_xxxz[k] = -g_0_z_yyyyz_xxxz[k] * ab_x + g_0_z_yyyyz_xxxxz[k];

                g_0_z_xyyyyz_xxyy[k] = -g_0_z_yyyyz_xxyy[k] * ab_x + g_0_z_yyyyz_xxxyy[k];

                g_0_z_xyyyyz_xxyz[k] = -g_0_z_yyyyz_xxyz[k] * ab_x + g_0_z_yyyyz_xxxyz[k];

                g_0_z_xyyyyz_xxzz[k] = -g_0_z_yyyyz_xxzz[k] * ab_x + g_0_z_yyyyz_xxxzz[k];

                g_0_z_xyyyyz_xyyy[k] = -g_0_z_yyyyz_xyyy[k] * ab_x + g_0_z_yyyyz_xxyyy[k];

                g_0_z_xyyyyz_xyyz[k] = -g_0_z_yyyyz_xyyz[k] * ab_x + g_0_z_yyyyz_xxyyz[k];

                g_0_z_xyyyyz_xyzz[k] = -g_0_z_yyyyz_xyzz[k] * ab_x + g_0_z_yyyyz_xxyzz[k];

                g_0_z_xyyyyz_xzzz[k] = -g_0_z_yyyyz_xzzz[k] * ab_x + g_0_z_yyyyz_xxzzz[k];

                g_0_z_xyyyyz_yyyy[k] = -g_0_z_yyyyz_yyyy[k] * ab_x + g_0_z_yyyyz_xyyyy[k];

                g_0_z_xyyyyz_yyyz[k] = -g_0_z_yyyyz_yyyz[k] * ab_x + g_0_z_yyyyz_xyyyz[k];

                g_0_z_xyyyyz_yyzz[k] = -g_0_z_yyyyz_yyzz[k] * ab_x + g_0_z_yyyyz_xyyzz[k];

                g_0_z_xyyyyz_yzzz[k] = -g_0_z_yyyyz_yzzz[k] * ab_x + g_0_z_yyyyz_xyzzz[k];

                g_0_z_xyyyyz_zzzz[k] = -g_0_z_yyyyz_zzzz[k] * ab_x + g_0_z_yyyyz_xzzzz[k];
            }

            /// Set up 1095-1110 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyzz_xxxx = cbuffer.data(ig_geom_01_off + 1095 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxxy = cbuffer.data(ig_geom_01_off + 1096 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxxz = cbuffer.data(ig_geom_01_off + 1097 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxyy = cbuffer.data(ig_geom_01_off + 1098 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxyz = cbuffer.data(ig_geom_01_off + 1099 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxzz = cbuffer.data(ig_geom_01_off + 1100 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyyy = cbuffer.data(ig_geom_01_off + 1101 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyyz = cbuffer.data(ig_geom_01_off + 1102 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyzz = cbuffer.data(ig_geom_01_off + 1103 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xzzz = cbuffer.data(ig_geom_01_off + 1104 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyyy = cbuffer.data(ig_geom_01_off + 1105 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyyz = cbuffer.data(ig_geom_01_off + 1106 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyzz = cbuffer.data(ig_geom_01_off + 1107 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yzzz = cbuffer.data(ig_geom_01_off + 1108 * ccomps * dcomps);

            auto g_0_z_xyyyzz_zzzz = cbuffer.data(ig_geom_01_off + 1109 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyzz_xxxx, g_0_z_xyyyzz_xxxy, g_0_z_xyyyzz_xxxz, g_0_z_xyyyzz_xxyy, g_0_z_xyyyzz_xxyz, g_0_z_xyyyzz_xxzz, g_0_z_xyyyzz_xyyy, g_0_z_xyyyzz_xyyz, g_0_z_xyyyzz_xyzz, g_0_z_xyyyzz_xzzz, g_0_z_xyyyzz_yyyy, g_0_z_xyyyzz_yyyz, g_0_z_xyyyzz_yyzz, g_0_z_xyyyzz_yzzz, g_0_z_xyyyzz_zzzz, g_0_z_yyyzz_xxxx, g_0_z_yyyzz_xxxxx, g_0_z_yyyzz_xxxxy, g_0_z_yyyzz_xxxxz, g_0_z_yyyzz_xxxy, g_0_z_yyyzz_xxxyy, g_0_z_yyyzz_xxxyz, g_0_z_yyyzz_xxxz, g_0_z_yyyzz_xxxzz, g_0_z_yyyzz_xxyy, g_0_z_yyyzz_xxyyy, g_0_z_yyyzz_xxyyz, g_0_z_yyyzz_xxyz, g_0_z_yyyzz_xxyzz, g_0_z_yyyzz_xxzz, g_0_z_yyyzz_xxzzz, g_0_z_yyyzz_xyyy, g_0_z_yyyzz_xyyyy, g_0_z_yyyzz_xyyyz, g_0_z_yyyzz_xyyz, g_0_z_yyyzz_xyyzz, g_0_z_yyyzz_xyzz, g_0_z_yyyzz_xyzzz, g_0_z_yyyzz_xzzz, g_0_z_yyyzz_xzzzz, g_0_z_yyyzz_yyyy, g_0_z_yyyzz_yyyz, g_0_z_yyyzz_yyzz, g_0_z_yyyzz_yzzz, g_0_z_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyzz_xxxx[k] = -g_0_z_yyyzz_xxxx[k] * ab_x + g_0_z_yyyzz_xxxxx[k];

                g_0_z_xyyyzz_xxxy[k] = -g_0_z_yyyzz_xxxy[k] * ab_x + g_0_z_yyyzz_xxxxy[k];

                g_0_z_xyyyzz_xxxz[k] = -g_0_z_yyyzz_xxxz[k] * ab_x + g_0_z_yyyzz_xxxxz[k];

                g_0_z_xyyyzz_xxyy[k] = -g_0_z_yyyzz_xxyy[k] * ab_x + g_0_z_yyyzz_xxxyy[k];

                g_0_z_xyyyzz_xxyz[k] = -g_0_z_yyyzz_xxyz[k] * ab_x + g_0_z_yyyzz_xxxyz[k];

                g_0_z_xyyyzz_xxzz[k] = -g_0_z_yyyzz_xxzz[k] * ab_x + g_0_z_yyyzz_xxxzz[k];

                g_0_z_xyyyzz_xyyy[k] = -g_0_z_yyyzz_xyyy[k] * ab_x + g_0_z_yyyzz_xxyyy[k];

                g_0_z_xyyyzz_xyyz[k] = -g_0_z_yyyzz_xyyz[k] * ab_x + g_0_z_yyyzz_xxyyz[k];

                g_0_z_xyyyzz_xyzz[k] = -g_0_z_yyyzz_xyzz[k] * ab_x + g_0_z_yyyzz_xxyzz[k];

                g_0_z_xyyyzz_xzzz[k] = -g_0_z_yyyzz_xzzz[k] * ab_x + g_0_z_yyyzz_xxzzz[k];

                g_0_z_xyyyzz_yyyy[k] = -g_0_z_yyyzz_yyyy[k] * ab_x + g_0_z_yyyzz_xyyyy[k];

                g_0_z_xyyyzz_yyyz[k] = -g_0_z_yyyzz_yyyz[k] * ab_x + g_0_z_yyyzz_xyyyz[k];

                g_0_z_xyyyzz_yyzz[k] = -g_0_z_yyyzz_yyzz[k] * ab_x + g_0_z_yyyzz_xyyzz[k];

                g_0_z_xyyyzz_yzzz[k] = -g_0_z_yyyzz_yzzz[k] * ab_x + g_0_z_yyyzz_xyzzz[k];

                g_0_z_xyyyzz_zzzz[k] = -g_0_z_yyyzz_zzzz[k] * ab_x + g_0_z_yyyzz_xzzzz[k];
            }

            /// Set up 1110-1125 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzzz_xxxx = cbuffer.data(ig_geom_01_off + 1110 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxxy = cbuffer.data(ig_geom_01_off + 1111 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxxz = cbuffer.data(ig_geom_01_off + 1112 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxyy = cbuffer.data(ig_geom_01_off + 1113 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxyz = cbuffer.data(ig_geom_01_off + 1114 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxzz = cbuffer.data(ig_geom_01_off + 1115 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyyy = cbuffer.data(ig_geom_01_off + 1116 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyyz = cbuffer.data(ig_geom_01_off + 1117 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyzz = cbuffer.data(ig_geom_01_off + 1118 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xzzz = cbuffer.data(ig_geom_01_off + 1119 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyyy = cbuffer.data(ig_geom_01_off + 1120 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyyz = cbuffer.data(ig_geom_01_off + 1121 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyzz = cbuffer.data(ig_geom_01_off + 1122 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yzzz = cbuffer.data(ig_geom_01_off + 1123 * ccomps * dcomps);

            auto g_0_z_xyyzzz_zzzz = cbuffer.data(ig_geom_01_off + 1124 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzzz_xxxx, g_0_z_xyyzzz_xxxy, g_0_z_xyyzzz_xxxz, g_0_z_xyyzzz_xxyy, g_0_z_xyyzzz_xxyz, g_0_z_xyyzzz_xxzz, g_0_z_xyyzzz_xyyy, g_0_z_xyyzzz_xyyz, g_0_z_xyyzzz_xyzz, g_0_z_xyyzzz_xzzz, g_0_z_xyyzzz_yyyy, g_0_z_xyyzzz_yyyz, g_0_z_xyyzzz_yyzz, g_0_z_xyyzzz_yzzz, g_0_z_xyyzzz_zzzz, g_0_z_yyzzz_xxxx, g_0_z_yyzzz_xxxxx, g_0_z_yyzzz_xxxxy, g_0_z_yyzzz_xxxxz, g_0_z_yyzzz_xxxy, g_0_z_yyzzz_xxxyy, g_0_z_yyzzz_xxxyz, g_0_z_yyzzz_xxxz, g_0_z_yyzzz_xxxzz, g_0_z_yyzzz_xxyy, g_0_z_yyzzz_xxyyy, g_0_z_yyzzz_xxyyz, g_0_z_yyzzz_xxyz, g_0_z_yyzzz_xxyzz, g_0_z_yyzzz_xxzz, g_0_z_yyzzz_xxzzz, g_0_z_yyzzz_xyyy, g_0_z_yyzzz_xyyyy, g_0_z_yyzzz_xyyyz, g_0_z_yyzzz_xyyz, g_0_z_yyzzz_xyyzz, g_0_z_yyzzz_xyzz, g_0_z_yyzzz_xyzzz, g_0_z_yyzzz_xzzz, g_0_z_yyzzz_xzzzz, g_0_z_yyzzz_yyyy, g_0_z_yyzzz_yyyz, g_0_z_yyzzz_yyzz, g_0_z_yyzzz_yzzz, g_0_z_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzzz_xxxx[k] = -g_0_z_yyzzz_xxxx[k] * ab_x + g_0_z_yyzzz_xxxxx[k];

                g_0_z_xyyzzz_xxxy[k] = -g_0_z_yyzzz_xxxy[k] * ab_x + g_0_z_yyzzz_xxxxy[k];

                g_0_z_xyyzzz_xxxz[k] = -g_0_z_yyzzz_xxxz[k] * ab_x + g_0_z_yyzzz_xxxxz[k];

                g_0_z_xyyzzz_xxyy[k] = -g_0_z_yyzzz_xxyy[k] * ab_x + g_0_z_yyzzz_xxxyy[k];

                g_0_z_xyyzzz_xxyz[k] = -g_0_z_yyzzz_xxyz[k] * ab_x + g_0_z_yyzzz_xxxyz[k];

                g_0_z_xyyzzz_xxzz[k] = -g_0_z_yyzzz_xxzz[k] * ab_x + g_0_z_yyzzz_xxxzz[k];

                g_0_z_xyyzzz_xyyy[k] = -g_0_z_yyzzz_xyyy[k] * ab_x + g_0_z_yyzzz_xxyyy[k];

                g_0_z_xyyzzz_xyyz[k] = -g_0_z_yyzzz_xyyz[k] * ab_x + g_0_z_yyzzz_xxyyz[k];

                g_0_z_xyyzzz_xyzz[k] = -g_0_z_yyzzz_xyzz[k] * ab_x + g_0_z_yyzzz_xxyzz[k];

                g_0_z_xyyzzz_xzzz[k] = -g_0_z_yyzzz_xzzz[k] * ab_x + g_0_z_yyzzz_xxzzz[k];

                g_0_z_xyyzzz_yyyy[k] = -g_0_z_yyzzz_yyyy[k] * ab_x + g_0_z_yyzzz_xyyyy[k];

                g_0_z_xyyzzz_yyyz[k] = -g_0_z_yyzzz_yyyz[k] * ab_x + g_0_z_yyzzz_xyyyz[k];

                g_0_z_xyyzzz_yyzz[k] = -g_0_z_yyzzz_yyzz[k] * ab_x + g_0_z_yyzzz_xyyzz[k];

                g_0_z_xyyzzz_yzzz[k] = -g_0_z_yyzzz_yzzz[k] * ab_x + g_0_z_yyzzz_xyzzz[k];

                g_0_z_xyyzzz_zzzz[k] = -g_0_z_yyzzz_zzzz[k] * ab_x + g_0_z_yyzzz_xzzzz[k];
            }

            /// Set up 1125-1140 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzzz_xxxx = cbuffer.data(ig_geom_01_off + 1125 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxxy = cbuffer.data(ig_geom_01_off + 1126 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxxz = cbuffer.data(ig_geom_01_off + 1127 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxyy = cbuffer.data(ig_geom_01_off + 1128 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxyz = cbuffer.data(ig_geom_01_off + 1129 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxzz = cbuffer.data(ig_geom_01_off + 1130 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyyy = cbuffer.data(ig_geom_01_off + 1131 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyyz = cbuffer.data(ig_geom_01_off + 1132 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyzz = cbuffer.data(ig_geom_01_off + 1133 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xzzz = cbuffer.data(ig_geom_01_off + 1134 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyyy = cbuffer.data(ig_geom_01_off + 1135 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyyz = cbuffer.data(ig_geom_01_off + 1136 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyzz = cbuffer.data(ig_geom_01_off + 1137 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yzzz = cbuffer.data(ig_geom_01_off + 1138 * ccomps * dcomps);

            auto g_0_z_xyzzzz_zzzz = cbuffer.data(ig_geom_01_off + 1139 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzzz_xxxx, g_0_z_xyzzzz_xxxy, g_0_z_xyzzzz_xxxz, g_0_z_xyzzzz_xxyy, g_0_z_xyzzzz_xxyz, g_0_z_xyzzzz_xxzz, g_0_z_xyzzzz_xyyy, g_0_z_xyzzzz_xyyz, g_0_z_xyzzzz_xyzz, g_0_z_xyzzzz_xzzz, g_0_z_xyzzzz_yyyy, g_0_z_xyzzzz_yyyz, g_0_z_xyzzzz_yyzz, g_0_z_xyzzzz_yzzz, g_0_z_xyzzzz_zzzz, g_0_z_yzzzz_xxxx, g_0_z_yzzzz_xxxxx, g_0_z_yzzzz_xxxxy, g_0_z_yzzzz_xxxxz, g_0_z_yzzzz_xxxy, g_0_z_yzzzz_xxxyy, g_0_z_yzzzz_xxxyz, g_0_z_yzzzz_xxxz, g_0_z_yzzzz_xxxzz, g_0_z_yzzzz_xxyy, g_0_z_yzzzz_xxyyy, g_0_z_yzzzz_xxyyz, g_0_z_yzzzz_xxyz, g_0_z_yzzzz_xxyzz, g_0_z_yzzzz_xxzz, g_0_z_yzzzz_xxzzz, g_0_z_yzzzz_xyyy, g_0_z_yzzzz_xyyyy, g_0_z_yzzzz_xyyyz, g_0_z_yzzzz_xyyz, g_0_z_yzzzz_xyyzz, g_0_z_yzzzz_xyzz, g_0_z_yzzzz_xyzzz, g_0_z_yzzzz_xzzz, g_0_z_yzzzz_xzzzz, g_0_z_yzzzz_yyyy, g_0_z_yzzzz_yyyz, g_0_z_yzzzz_yyzz, g_0_z_yzzzz_yzzz, g_0_z_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzzz_xxxx[k] = -g_0_z_yzzzz_xxxx[k] * ab_x + g_0_z_yzzzz_xxxxx[k];

                g_0_z_xyzzzz_xxxy[k] = -g_0_z_yzzzz_xxxy[k] * ab_x + g_0_z_yzzzz_xxxxy[k];

                g_0_z_xyzzzz_xxxz[k] = -g_0_z_yzzzz_xxxz[k] * ab_x + g_0_z_yzzzz_xxxxz[k];

                g_0_z_xyzzzz_xxyy[k] = -g_0_z_yzzzz_xxyy[k] * ab_x + g_0_z_yzzzz_xxxyy[k];

                g_0_z_xyzzzz_xxyz[k] = -g_0_z_yzzzz_xxyz[k] * ab_x + g_0_z_yzzzz_xxxyz[k];

                g_0_z_xyzzzz_xxzz[k] = -g_0_z_yzzzz_xxzz[k] * ab_x + g_0_z_yzzzz_xxxzz[k];

                g_0_z_xyzzzz_xyyy[k] = -g_0_z_yzzzz_xyyy[k] * ab_x + g_0_z_yzzzz_xxyyy[k];

                g_0_z_xyzzzz_xyyz[k] = -g_0_z_yzzzz_xyyz[k] * ab_x + g_0_z_yzzzz_xxyyz[k];

                g_0_z_xyzzzz_xyzz[k] = -g_0_z_yzzzz_xyzz[k] * ab_x + g_0_z_yzzzz_xxyzz[k];

                g_0_z_xyzzzz_xzzz[k] = -g_0_z_yzzzz_xzzz[k] * ab_x + g_0_z_yzzzz_xxzzz[k];

                g_0_z_xyzzzz_yyyy[k] = -g_0_z_yzzzz_yyyy[k] * ab_x + g_0_z_yzzzz_xyyyy[k];

                g_0_z_xyzzzz_yyyz[k] = -g_0_z_yzzzz_yyyz[k] * ab_x + g_0_z_yzzzz_xyyyz[k];

                g_0_z_xyzzzz_yyzz[k] = -g_0_z_yzzzz_yyzz[k] * ab_x + g_0_z_yzzzz_xyyzz[k];

                g_0_z_xyzzzz_yzzz[k] = -g_0_z_yzzzz_yzzz[k] * ab_x + g_0_z_yzzzz_xyzzz[k];

                g_0_z_xyzzzz_zzzz[k] = -g_0_z_yzzzz_zzzz[k] * ab_x + g_0_z_yzzzz_xzzzz[k];
            }

            /// Set up 1140-1155 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzzz_xxxx = cbuffer.data(ig_geom_01_off + 1140 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxxy = cbuffer.data(ig_geom_01_off + 1141 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxxz = cbuffer.data(ig_geom_01_off + 1142 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxyy = cbuffer.data(ig_geom_01_off + 1143 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxyz = cbuffer.data(ig_geom_01_off + 1144 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxzz = cbuffer.data(ig_geom_01_off + 1145 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyyy = cbuffer.data(ig_geom_01_off + 1146 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyyz = cbuffer.data(ig_geom_01_off + 1147 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyzz = cbuffer.data(ig_geom_01_off + 1148 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xzzz = cbuffer.data(ig_geom_01_off + 1149 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyyy = cbuffer.data(ig_geom_01_off + 1150 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyyz = cbuffer.data(ig_geom_01_off + 1151 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyzz = cbuffer.data(ig_geom_01_off + 1152 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yzzz = cbuffer.data(ig_geom_01_off + 1153 * ccomps * dcomps);

            auto g_0_z_xzzzzz_zzzz = cbuffer.data(ig_geom_01_off + 1154 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzzz_xxxx, g_0_z_xzzzzz_xxxy, g_0_z_xzzzzz_xxxz, g_0_z_xzzzzz_xxyy, g_0_z_xzzzzz_xxyz, g_0_z_xzzzzz_xxzz, g_0_z_xzzzzz_xyyy, g_0_z_xzzzzz_xyyz, g_0_z_xzzzzz_xyzz, g_0_z_xzzzzz_xzzz, g_0_z_xzzzzz_yyyy, g_0_z_xzzzzz_yyyz, g_0_z_xzzzzz_yyzz, g_0_z_xzzzzz_yzzz, g_0_z_xzzzzz_zzzz, g_0_z_zzzzz_xxxx, g_0_z_zzzzz_xxxxx, g_0_z_zzzzz_xxxxy, g_0_z_zzzzz_xxxxz, g_0_z_zzzzz_xxxy, g_0_z_zzzzz_xxxyy, g_0_z_zzzzz_xxxyz, g_0_z_zzzzz_xxxz, g_0_z_zzzzz_xxxzz, g_0_z_zzzzz_xxyy, g_0_z_zzzzz_xxyyy, g_0_z_zzzzz_xxyyz, g_0_z_zzzzz_xxyz, g_0_z_zzzzz_xxyzz, g_0_z_zzzzz_xxzz, g_0_z_zzzzz_xxzzz, g_0_z_zzzzz_xyyy, g_0_z_zzzzz_xyyyy, g_0_z_zzzzz_xyyyz, g_0_z_zzzzz_xyyz, g_0_z_zzzzz_xyyzz, g_0_z_zzzzz_xyzz, g_0_z_zzzzz_xyzzz, g_0_z_zzzzz_xzzz, g_0_z_zzzzz_xzzzz, g_0_z_zzzzz_yyyy, g_0_z_zzzzz_yyyz, g_0_z_zzzzz_yyzz, g_0_z_zzzzz_yzzz, g_0_z_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzzz_xxxx[k] = -g_0_z_zzzzz_xxxx[k] * ab_x + g_0_z_zzzzz_xxxxx[k];

                g_0_z_xzzzzz_xxxy[k] = -g_0_z_zzzzz_xxxy[k] * ab_x + g_0_z_zzzzz_xxxxy[k];

                g_0_z_xzzzzz_xxxz[k] = -g_0_z_zzzzz_xxxz[k] * ab_x + g_0_z_zzzzz_xxxxz[k];

                g_0_z_xzzzzz_xxyy[k] = -g_0_z_zzzzz_xxyy[k] * ab_x + g_0_z_zzzzz_xxxyy[k];

                g_0_z_xzzzzz_xxyz[k] = -g_0_z_zzzzz_xxyz[k] * ab_x + g_0_z_zzzzz_xxxyz[k];

                g_0_z_xzzzzz_xxzz[k] = -g_0_z_zzzzz_xxzz[k] * ab_x + g_0_z_zzzzz_xxxzz[k];

                g_0_z_xzzzzz_xyyy[k] = -g_0_z_zzzzz_xyyy[k] * ab_x + g_0_z_zzzzz_xxyyy[k];

                g_0_z_xzzzzz_xyyz[k] = -g_0_z_zzzzz_xyyz[k] * ab_x + g_0_z_zzzzz_xxyyz[k];

                g_0_z_xzzzzz_xyzz[k] = -g_0_z_zzzzz_xyzz[k] * ab_x + g_0_z_zzzzz_xxyzz[k];

                g_0_z_xzzzzz_xzzz[k] = -g_0_z_zzzzz_xzzz[k] * ab_x + g_0_z_zzzzz_xxzzz[k];

                g_0_z_xzzzzz_yyyy[k] = -g_0_z_zzzzz_yyyy[k] * ab_x + g_0_z_zzzzz_xyyyy[k];

                g_0_z_xzzzzz_yyyz[k] = -g_0_z_zzzzz_yyyz[k] * ab_x + g_0_z_zzzzz_xyyyz[k];

                g_0_z_xzzzzz_yyzz[k] = -g_0_z_zzzzz_yyzz[k] * ab_x + g_0_z_zzzzz_xyyzz[k];

                g_0_z_xzzzzz_yzzz[k] = -g_0_z_zzzzz_yzzz[k] * ab_x + g_0_z_zzzzz_xyzzz[k];

                g_0_z_xzzzzz_zzzz[k] = -g_0_z_zzzzz_zzzz[k] * ab_x + g_0_z_zzzzz_xzzzz[k];
            }

            /// Set up 1155-1170 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyy_xxxx = cbuffer.data(ig_geom_01_off + 1155 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxxy = cbuffer.data(ig_geom_01_off + 1156 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxxz = cbuffer.data(ig_geom_01_off + 1157 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxyy = cbuffer.data(ig_geom_01_off + 1158 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxyz = cbuffer.data(ig_geom_01_off + 1159 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxzz = cbuffer.data(ig_geom_01_off + 1160 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyyy = cbuffer.data(ig_geom_01_off + 1161 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyyz = cbuffer.data(ig_geom_01_off + 1162 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyzz = cbuffer.data(ig_geom_01_off + 1163 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xzzz = cbuffer.data(ig_geom_01_off + 1164 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyyy = cbuffer.data(ig_geom_01_off + 1165 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyyz = cbuffer.data(ig_geom_01_off + 1166 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyzz = cbuffer.data(ig_geom_01_off + 1167 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yzzz = cbuffer.data(ig_geom_01_off + 1168 * ccomps * dcomps);

            auto g_0_z_yyyyyy_zzzz = cbuffer.data(ig_geom_01_off + 1169 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyy_xxxx, g_0_z_yyyyy_xxxxy, g_0_z_yyyyy_xxxy, g_0_z_yyyyy_xxxyy, g_0_z_yyyyy_xxxyz, g_0_z_yyyyy_xxxz, g_0_z_yyyyy_xxyy, g_0_z_yyyyy_xxyyy, g_0_z_yyyyy_xxyyz, g_0_z_yyyyy_xxyz, g_0_z_yyyyy_xxyzz, g_0_z_yyyyy_xxzz, g_0_z_yyyyy_xyyy, g_0_z_yyyyy_xyyyy, g_0_z_yyyyy_xyyyz, g_0_z_yyyyy_xyyz, g_0_z_yyyyy_xyyzz, g_0_z_yyyyy_xyzz, g_0_z_yyyyy_xyzzz, g_0_z_yyyyy_xzzz, g_0_z_yyyyy_yyyy, g_0_z_yyyyy_yyyyy, g_0_z_yyyyy_yyyyz, g_0_z_yyyyy_yyyz, g_0_z_yyyyy_yyyzz, g_0_z_yyyyy_yyzz, g_0_z_yyyyy_yyzzz, g_0_z_yyyyy_yzzz, g_0_z_yyyyy_yzzzz, g_0_z_yyyyy_zzzz, g_0_z_yyyyyy_xxxx, g_0_z_yyyyyy_xxxy, g_0_z_yyyyyy_xxxz, g_0_z_yyyyyy_xxyy, g_0_z_yyyyyy_xxyz, g_0_z_yyyyyy_xxzz, g_0_z_yyyyyy_xyyy, g_0_z_yyyyyy_xyyz, g_0_z_yyyyyy_xyzz, g_0_z_yyyyyy_xzzz, g_0_z_yyyyyy_yyyy, g_0_z_yyyyyy_yyyz, g_0_z_yyyyyy_yyzz, g_0_z_yyyyyy_yzzz, g_0_z_yyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyy_xxxx[k] = -g_0_z_yyyyy_xxxx[k] * ab_y + g_0_z_yyyyy_xxxxy[k];

                g_0_z_yyyyyy_xxxy[k] = -g_0_z_yyyyy_xxxy[k] * ab_y + g_0_z_yyyyy_xxxyy[k];

                g_0_z_yyyyyy_xxxz[k] = -g_0_z_yyyyy_xxxz[k] * ab_y + g_0_z_yyyyy_xxxyz[k];

                g_0_z_yyyyyy_xxyy[k] = -g_0_z_yyyyy_xxyy[k] * ab_y + g_0_z_yyyyy_xxyyy[k];

                g_0_z_yyyyyy_xxyz[k] = -g_0_z_yyyyy_xxyz[k] * ab_y + g_0_z_yyyyy_xxyyz[k];

                g_0_z_yyyyyy_xxzz[k] = -g_0_z_yyyyy_xxzz[k] * ab_y + g_0_z_yyyyy_xxyzz[k];

                g_0_z_yyyyyy_xyyy[k] = -g_0_z_yyyyy_xyyy[k] * ab_y + g_0_z_yyyyy_xyyyy[k];

                g_0_z_yyyyyy_xyyz[k] = -g_0_z_yyyyy_xyyz[k] * ab_y + g_0_z_yyyyy_xyyyz[k];

                g_0_z_yyyyyy_xyzz[k] = -g_0_z_yyyyy_xyzz[k] * ab_y + g_0_z_yyyyy_xyyzz[k];

                g_0_z_yyyyyy_xzzz[k] = -g_0_z_yyyyy_xzzz[k] * ab_y + g_0_z_yyyyy_xyzzz[k];

                g_0_z_yyyyyy_yyyy[k] = -g_0_z_yyyyy_yyyy[k] * ab_y + g_0_z_yyyyy_yyyyy[k];

                g_0_z_yyyyyy_yyyz[k] = -g_0_z_yyyyy_yyyz[k] * ab_y + g_0_z_yyyyy_yyyyz[k];

                g_0_z_yyyyyy_yyzz[k] = -g_0_z_yyyyy_yyzz[k] * ab_y + g_0_z_yyyyy_yyyzz[k];

                g_0_z_yyyyyy_yzzz[k] = -g_0_z_yyyyy_yzzz[k] * ab_y + g_0_z_yyyyy_yyzzz[k];

                g_0_z_yyyyyy_zzzz[k] = -g_0_z_yyyyy_zzzz[k] * ab_y + g_0_z_yyyyy_yzzzz[k];
            }

            /// Set up 1170-1185 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyz_xxxx = cbuffer.data(ig_geom_01_off + 1170 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxxy = cbuffer.data(ig_geom_01_off + 1171 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxxz = cbuffer.data(ig_geom_01_off + 1172 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxyy = cbuffer.data(ig_geom_01_off + 1173 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxyz = cbuffer.data(ig_geom_01_off + 1174 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxzz = cbuffer.data(ig_geom_01_off + 1175 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyyy = cbuffer.data(ig_geom_01_off + 1176 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyyz = cbuffer.data(ig_geom_01_off + 1177 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyzz = cbuffer.data(ig_geom_01_off + 1178 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xzzz = cbuffer.data(ig_geom_01_off + 1179 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyyy = cbuffer.data(ig_geom_01_off + 1180 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyyz = cbuffer.data(ig_geom_01_off + 1181 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyzz = cbuffer.data(ig_geom_01_off + 1182 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yzzz = cbuffer.data(ig_geom_01_off + 1183 * ccomps * dcomps);

            auto g_0_z_yyyyyz_zzzz = cbuffer.data(ig_geom_01_off + 1184 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyz_xxxx, g_0_z_yyyyyz_xxxy, g_0_z_yyyyyz_xxxz, g_0_z_yyyyyz_xxyy, g_0_z_yyyyyz_xxyz, g_0_z_yyyyyz_xxzz, g_0_z_yyyyyz_xyyy, g_0_z_yyyyyz_xyyz, g_0_z_yyyyyz_xyzz, g_0_z_yyyyyz_xzzz, g_0_z_yyyyyz_yyyy, g_0_z_yyyyyz_yyyz, g_0_z_yyyyyz_yyzz, g_0_z_yyyyyz_yzzz, g_0_z_yyyyyz_zzzz, g_0_z_yyyyz_xxxx, g_0_z_yyyyz_xxxxy, g_0_z_yyyyz_xxxy, g_0_z_yyyyz_xxxyy, g_0_z_yyyyz_xxxyz, g_0_z_yyyyz_xxxz, g_0_z_yyyyz_xxyy, g_0_z_yyyyz_xxyyy, g_0_z_yyyyz_xxyyz, g_0_z_yyyyz_xxyz, g_0_z_yyyyz_xxyzz, g_0_z_yyyyz_xxzz, g_0_z_yyyyz_xyyy, g_0_z_yyyyz_xyyyy, g_0_z_yyyyz_xyyyz, g_0_z_yyyyz_xyyz, g_0_z_yyyyz_xyyzz, g_0_z_yyyyz_xyzz, g_0_z_yyyyz_xyzzz, g_0_z_yyyyz_xzzz, g_0_z_yyyyz_yyyy, g_0_z_yyyyz_yyyyy, g_0_z_yyyyz_yyyyz, g_0_z_yyyyz_yyyz, g_0_z_yyyyz_yyyzz, g_0_z_yyyyz_yyzz, g_0_z_yyyyz_yyzzz, g_0_z_yyyyz_yzzz, g_0_z_yyyyz_yzzzz, g_0_z_yyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyz_xxxx[k] = -g_0_z_yyyyz_xxxx[k] * ab_y + g_0_z_yyyyz_xxxxy[k];

                g_0_z_yyyyyz_xxxy[k] = -g_0_z_yyyyz_xxxy[k] * ab_y + g_0_z_yyyyz_xxxyy[k];

                g_0_z_yyyyyz_xxxz[k] = -g_0_z_yyyyz_xxxz[k] * ab_y + g_0_z_yyyyz_xxxyz[k];

                g_0_z_yyyyyz_xxyy[k] = -g_0_z_yyyyz_xxyy[k] * ab_y + g_0_z_yyyyz_xxyyy[k];

                g_0_z_yyyyyz_xxyz[k] = -g_0_z_yyyyz_xxyz[k] * ab_y + g_0_z_yyyyz_xxyyz[k];

                g_0_z_yyyyyz_xxzz[k] = -g_0_z_yyyyz_xxzz[k] * ab_y + g_0_z_yyyyz_xxyzz[k];

                g_0_z_yyyyyz_xyyy[k] = -g_0_z_yyyyz_xyyy[k] * ab_y + g_0_z_yyyyz_xyyyy[k];

                g_0_z_yyyyyz_xyyz[k] = -g_0_z_yyyyz_xyyz[k] * ab_y + g_0_z_yyyyz_xyyyz[k];

                g_0_z_yyyyyz_xyzz[k] = -g_0_z_yyyyz_xyzz[k] * ab_y + g_0_z_yyyyz_xyyzz[k];

                g_0_z_yyyyyz_xzzz[k] = -g_0_z_yyyyz_xzzz[k] * ab_y + g_0_z_yyyyz_xyzzz[k];

                g_0_z_yyyyyz_yyyy[k] = -g_0_z_yyyyz_yyyy[k] * ab_y + g_0_z_yyyyz_yyyyy[k];

                g_0_z_yyyyyz_yyyz[k] = -g_0_z_yyyyz_yyyz[k] * ab_y + g_0_z_yyyyz_yyyyz[k];

                g_0_z_yyyyyz_yyzz[k] = -g_0_z_yyyyz_yyzz[k] * ab_y + g_0_z_yyyyz_yyyzz[k];

                g_0_z_yyyyyz_yzzz[k] = -g_0_z_yyyyz_yzzz[k] * ab_y + g_0_z_yyyyz_yyzzz[k];

                g_0_z_yyyyyz_zzzz[k] = -g_0_z_yyyyz_zzzz[k] * ab_y + g_0_z_yyyyz_yzzzz[k];
            }

            /// Set up 1185-1200 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyzz_xxxx = cbuffer.data(ig_geom_01_off + 1185 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxxy = cbuffer.data(ig_geom_01_off + 1186 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxxz = cbuffer.data(ig_geom_01_off + 1187 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxyy = cbuffer.data(ig_geom_01_off + 1188 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxyz = cbuffer.data(ig_geom_01_off + 1189 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxzz = cbuffer.data(ig_geom_01_off + 1190 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyyy = cbuffer.data(ig_geom_01_off + 1191 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyyz = cbuffer.data(ig_geom_01_off + 1192 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyzz = cbuffer.data(ig_geom_01_off + 1193 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xzzz = cbuffer.data(ig_geom_01_off + 1194 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyyy = cbuffer.data(ig_geom_01_off + 1195 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyyz = cbuffer.data(ig_geom_01_off + 1196 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyzz = cbuffer.data(ig_geom_01_off + 1197 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yzzz = cbuffer.data(ig_geom_01_off + 1198 * ccomps * dcomps);

            auto g_0_z_yyyyzz_zzzz = cbuffer.data(ig_geom_01_off + 1199 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyzz_xxxx, g_0_z_yyyyzz_xxxy, g_0_z_yyyyzz_xxxz, g_0_z_yyyyzz_xxyy, g_0_z_yyyyzz_xxyz, g_0_z_yyyyzz_xxzz, g_0_z_yyyyzz_xyyy, g_0_z_yyyyzz_xyyz, g_0_z_yyyyzz_xyzz, g_0_z_yyyyzz_xzzz, g_0_z_yyyyzz_yyyy, g_0_z_yyyyzz_yyyz, g_0_z_yyyyzz_yyzz, g_0_z_yyyyzz_yzzz, g_0_z_yyyyzz_zzzz, g_0_z_yyyzz_xxxx, g_0_z_yyyzz_xxxxy, g_0_z_yyyzz_xxxy, g_0_z_yyyzz_xxxyy, g_0_z_yyyzz_xxxyz, g_0_z_yyyzz_xxxz, g_0_z_yyyzz_xxyy, g_0_z_yyyzz_xxyyy, g_0_z_yyyzz_xxyyz, g_0_z_yyyzz_xxyz, g_0_z_yyyzz_xxyzz, g_0_z_yyyzz_xxzz, g_0_z_yyyzz_xyyy, g_0_z_yyyzz_xyyyy, g_0_z_yyyzz_xyyyz, g_0_z_yyyzz_xyyz, g_0_z_yyyzz_xyyzz, g_0_z_yyyzz_xyzz, g_0_z_yyyzz_xyzzz, g_0_z_yyyzz_xzzz, g_0_z_yyyzz_yyyy, g_0_z_yyyzz_yyyyy, g_0_z_yyyzz_yyyyz, g_0_z_yyyzz_yyyz, g_0_z_yyyzz_yyyzz, g_0_z_yyyzz_yyzz, g_0_z_yyyzz_yyzzz, g_0_z_yyyzz_yzzz, g_0_z_yyyzz_yzzzz, g_0_z_yyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyzz_xxxx[k] = -g_0_z_yyyzz_xxxx[k] * ab_y + g_0_z_yyyzz_xxxxy[k];

                g_0_z_yyyyzz_xxxy[k] = -g_0_z_yyyzz_xxxy[k] * ab_y + g_0_z_yyyzz_xxxyy[k];

                g_0_z_yyyyzz_xxxz[k] = -g_0_z_yyyzz_xxxz[k] * ab_y + g_0_z_yyyzz_xxxyz[k];

                g_0_z_yyyyzz_xxyy[k] = -g_0_z_yyyzz_xxyy[k] * ab_y + g_0_z_yyyzz_xxyyy[k];

                g_0_z_yyyyzz_xxyz[k] = -g_0_z_yyyzz_xxyz[k] * ab_y + g_0_z_yyyzz_xxyyz[k];

                g_0_z_yyyyzz_xxzz[k] = -g_0_z_yyyzz_xxzz[k] * ab_y + g_0_z_yyyzz_xxyzz[k];

                g_0_z_yyyyzz_xyyy[k] = -g_0_z_yyyzz_xyyy[k] * ab_y + g_0_z_yyyzz_xyyyy[k];

                g_0_z_yyyyzz_xyyz[k] = -g_0_z_yyyzz_xyyz[k] * ab_y + g_0_z_yyyzz_xyyyz[k];

                g_0_z_yyyyzz_xyzz[k] = -g_0_z_yyyzz_xyzz[k] * ab_y + g_0_z_yyyzz_xyyzz[k];

                g_0_z_yyyyzz_xzzz[k] = -g_0_z_yyyzz_xzzz[k] * ab_y + g_0_z_yyyzz_xyzzz[k];

                g_0_z_yyyyzz_yyyy[k] = -g_0_z_yyyzz_yyyy[k] * ab_y + g_0_z_yyyzz_yyyyy[k];

                g_0_z_yyyyzz_yyyz[k] = -g_0_z_yyyzz_yyyz[k] * ab_y + g_0_z_yyyzz_yyyyz[k];

                g_0_z_yyyyzz_yyzz[k] = -g_0_z_yyyzz_yyzz[k] * ab_y + g_0_z_yyyzz_yyyzz[k];

                g_0_z_yyyyzz_yzzz[k] = -g_0_z_yyyzz_yzzz[k] * ab_y + g_0_z_yyyzz_yyzzz[k];

                g_0_z_yyyyzz_zzzz[k] = -g_0_z_yyyzz_zzzz[k] * ab_y + g_0_z_yyyzz_yzzzz[k];
            }

            /// Set up 1200-1215 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzzz_xxxx = cbuffer.data(ig_geom_01_off + 1200 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxxy = cbuffer.data(ig_geom_01_off + 1201 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxxz = cbuffer.data(ig_geom_01_off + 1202 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxyy = cbuffer.data(ig_geom_01_off + 1203 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxyz = cbuffer.data(ig_geom_01_off + 1204 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxzz = cbuffer.data(ig_geom_01_off + 1205 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyyy = cbuffer.data(ig_geom_01_off + 1206 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyyz = cbuffer.data(ig_geom_01_off + 1207 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyzz = cbuffer.data(ig_geom_01_off + 1208 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xzzz = cbuffer.data(ig_geom_01_off + 1209 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyyy = cbuffer.data(ig_geom_01_off + 1210 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyyz = cbuffer.data(ig_geom_01_off + 1211 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyzz = cbuffer.data(ig_geom_01_off + 1212 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yzzz = cbuffer.data(ig_geom_01_off + 1213 * ccomps * dcomps);

            auto g_0_z_yyyzzz_zzzz = cbuffer.data(ig_geom_01_off + 1214 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzzz_xxxx, g_0_z_yyyzzz_xxxy, g_0_z_yyyzzz_xxxz, g_0_z_yyyzzz_xxyy, g_0_z_yyyzzz_xxyz, g_0_z_yyyzzz_xxzz, g_0_z_yyyzzz_xyyy, g_0_z_yyyzzz_xyyz, g_0_z_yyyzzz_xyzz, g_0_z_yyyzzz_xzzz, g_0_z_yyyzzz_yyyy, g_0_z_yyyzzz_yyyz, g_0_z_yyyzzz_yyzz, g_0_z_yyyzzz_yzzz, g_0_z_yyyzzz_zzzz, g_0_z_yyzzz_xxxx, g_0_z_yyzzz_xxxxy, g_0_z_yyzzz_xxxy, g_0_z_yyzzz_xxxyy, g_0_z_yyzzz_xxxyz, g_0_z_yyzzz_xxxz, g_0_z_yyzzz_xxyy, g_0_z_yyzzz_xxyyy, g_0_z_yyzzz_xxyyz, g_0_z_yyzzz_xxyz, g_0_z_yyzzz_xxyzz, g_0_z_yyzzz_xxzz, g_0_z_yyzzz_xyyy, g_0_z_yyzzz_xyyyy, g_0_z_yyzzz_xyyyz, g_0_z_yyzzz_xyyz, g_0_z_yyzzz_xyyzz, g_0_z_yyzzz_xyzz, g_0_z_yyzzz_xyzzz, g_0_z_yyzzz_xzzz, g_0_z_yyzzz_yyyy, g_0_z_yyzzz_yyyyy, g_0_z_yyzzz_yyyyz, g_0_z_yyzzz_yyyz, g_0_z_yyzzz_yyyzz, g_0_z_yyzzz_yyzz, g_0_z_yyzzz_yyzzz, g_0_z_yyzzz_yzzz, g_0_z_yyzzz_yzzzz, g_0_z_yyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzzz_xxxx[k] = -g_0_z_yyzzz_xxxx[k] * ab_y + g_0_z_yyzzz_xxxxy[k];

                g_0_z_yyyzzz_xxxy[k] = -g_0_z_yyzzz_xxxy[k] * ab_y + g_0_z_yyzzz_xxxyy[k];

                g_0_z_yyyzzz_xxxz[k] = -g_0_z_yyzzz_xxxz[k] * ab_y + g_0_z_yyzzz_xxxyz[k];

                g_0_z_yyyzzz_xxyy[k] = -g_0_z_yyzzz_xxyy[k] * ab_y + g_0_z_yyzzz_xxyyy[k];

                g_0_z_yyyzzz_xxyz[k] = -g_0_z_yyzzz_xxyz[k] * ab_y + g_0_z_yyzzz_xxyyz[k];

                g_0_z_yyyzzz_xxzz[k] = -g_0_z_yyzzz_xxzz[k] * ab_y + g_0_z_yyzzz_xxyzz[k];

                g_0_z_yyyzzz_xyyy[k] = -g_0_z_yyzzz_xyyy[k] * ab_y + g_0_z_yyzzz_xyyyy[k];

                g_0_z_yyyzzz_xyyz[k] = -g_0_z_yyzzz_xyyz[k] * ab_y + g_0_z_yyzzz_xyyyz[k];

                g_0_z_yyyzzz_xyzz[k] = -g_0_z_yyzzz_xyzz[k] * ab_y + g_0_z_yyzzz_xyyzz[k];

                g_0_z_yyyzzz_xzzz[k] = -g_0_z_yyzzz_xzzz[k] * ab_y + g_0_z_yyzzz_xyzzz[k];

                g_0_z_yyyzzz_yyyy[k] = -g_0_z_yyzzz_yyyy[k] * ab_y + g_0_z_yyzzz_yyyyy[k];

                g_0_z_yyyzzz_yyyz[k] = -g_0_z_yyzzz_yyyz[k] * ab_y + g_0_z_yyzzz_yyyyz[k];

                g_0_z_yyyzzz_yyzz[k] = -g_0_z_yyzzz_yyzz[k] * ab_y + g_0_z_yyzzz_yyyzz[k];

                g_0_z_yyyzzz_yzzz[k] = -g_0_z_yyzzz_yzzz[k] * ab_y + g_0_z_yyzzz_yyzzz[k];

                g_0_z_yyyzzz_zzzz[k] = -g_0_z_yyzzz_zzzz[k] * ab_y + g_0_z_yyzzz_yzzzz[k];
            }

            /// Set up 1215-1230 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzzz_xxxx = cbuffer.data(ig_geom_01_off + 1215 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxxy = cbuffer.data(ig_geom_01_off + 1216 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxxz = cbuffer.data(ig_geom_01_off + 1217 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxyy = cbuffer.data(ig_geom_01_off + 1218 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxyz = cbuffer.data(ig_geom_01_off + 1219 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxzz = cbuffer.data(ig_geom_01_off + 1220 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyyy = cbuffer.data(ig_geom_01_off + 1221 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyyz = cbuffer.data(ig_geom_01_off + 1222 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyzz = cbuffer.data(ig_geom_01_off + 1223 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xzzz = cbuffer.data(ig_geom_01_off + 1224 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyyy = cbuffer.data(ig_geom_01_off + 1225 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyyz = cbuffer.data(ig_geom_01_off + 1226 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyzz = cbuffer.data(ig_geom_01_off + 1227 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yzzz = cbuffer.data(ig_geom_01_off + 1228 * ccomps * dcomps);

            auto g_0_z_yyzzzz_zzzz = cbuffer.data(ig_geom_01_off + 1229 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzzz_xxxx, g_0_z_yyzzzz_xxxy, g_0_z_yyzzzz_xxxz, g_0_z_yyzzzz_xxyy, g_0_z_yyzzzz_xxyz, g_0_z_yyzzzz_xxzz, g_0_z_yyzzzz_xyyy, g_0_z_yyzzzz_xyyz, g_0_z_yyzzzz_xyzz, g_0_z_yyzzzz_xzzz, g_0_z_yyzzzz_yyyy, g_0_z_yyzzzz_yyyz, g_0_z_yyzzzz_yyzz, g_0_z_yyzzzz_yzzz, g_0_z_yyzzzz_zzzz, g_0_z_yzzzz_xxxx, g_0_z_yzzzz_xxxxy, g_0_z_yzzzz_xxxy, g_0_z_yzzzz_xxxyy, g_0_z_yzzzz_xxxyz, g_0_z_yzzzz_xxxz, g_0_z_yzzzz_xxyy, g_0_z_yzzzz_xxyyy, g_0_z_yzzzz_xxyyz, g_0_z_yzzzz_xxyz, g_0_z_yzzzz_xxyzz, g_0_z_yzzzz_xxzz, g_0_z_yzzzz_xyyy, g_0_z_yzzzz_xyyyy, g_0_z_yzzzz_xyyyz, g_0_z_yzzzz_xyyz, g_0_z_yzzzz_xyyzz, g_0_z_yzzzz_xyzz, g_0_z_yzzzz_xyzzz, g_0_z_yzzzz_xzzz, g_0_z_yzzzz_yyyy, g_0_z_yzzzz_yyyyy, g_0_z_yzzzz_yyyyz, g_0_z_yzzzz_yyyz, g_0_z_yzzzz_yyyzz, g_0_z_yzzzz_yyzz, g_0_z_yzzzz_yyzzz, g_0_z_yzzzz_yzzz, g_0_z_yzzzz_yzzzz, g_0_z_yzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzzz_xxxx[k] = -g_0_z_yzzzz_xxxx[k] * ab_y + g_0_z_yzzzz_xxxxy[k];

                g_0_z_yyzzzz_xxxy[k] = -g_0_z_yzzzz_xxxy[k] * ab_y + g_0_z_yzzzz_xxxyy[k];

                g_0_z_yyzzzz_xxxz[k] = -g_0_z_yzzzz_xxxz[k] * ab_y + g_0_z_yzzzz_xxxyz[k];

                g_0_z_yyzzzz_xxyy[k] = -g_0_z_yzzzz_xxyy[k] * ab_y + g_0_z_yzzzz_xxyyy[k];

                g_0_z_yyzzzz_xxyz[k] = -g_0_z_yzzzz_xxyz[k] * ab_y + g_0_z_yzzzz_xxyyz[k];

                g_0_z_yyzzzz_xxzz[k] = -g_0_z_yzzzz_xxzz[k] * ab_y + g_0_z_yzzzz_xxyzz[k];

                g_0_z_yyzzzz_xyyy[k] = -g_0_z_yzzzz_xyyy[k] * ab_y + g_0_z_yzzzz_xyyyy[k];

                g_0_z_yyzzzz_xyyz[k] = -g_0_z_yzzzz_xyyz[k] * ab_y + g_0_z_yzzzz_xyyyz[k];

                g_0_z_yyzzzz_xyzz[k] = -g_0_z_yzzzz_xyzz[k] * ab_y + g_0_z_yzzzz_xyyzz[k];

                g_0_z_yyzzzz_xzzz[k] = -g_0_z_yzzzz_xzzz[k] * ab_y + g_0_z_yzzzz_xyzzz[k];

                g_0_z_yyzzzz_yyyy[k] = -g_0_z_yzzzz_yyyy[k] * ab_y + g_0_z_yzzzz_yyyyy[k];

                g_0_z_yyzzzz_yyyz[k] = -g_0_z_yzzzz_yyyz[k] * ab_y + g_0_z_yzzzz_yyyyz[k];

                g_0_z_yyzzzz_yyzz[k] = -g_0_z_yzzzz_yyzz[k] * ab_y + g_0_z_yzzzz_yyyzz[k];

                g_0_z_yyzzzz_yzzz[k] = -g_0_z_yzzzz_yzzz[k] * ab_y + g_0_z_yzzzz_yyzzz[k];

                g_0_z_yyzzzz_zzzz[k] = -g_0_z_yzzzz_zzzz[k] * ab_y + g_0_z_yzzzz_yzzzz[k];
            }

            /// Set up 1230-1245 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzzz_xxxx = cbuffer.data(ig_geom_01_off + 1230 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxxy = cbuffer.data(ig_geom_01_off + 1231 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxxz = cbuffer.data(ig_geom_01_off + 1232 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxyy = cbuffer.data(ig_geom_01_off + 1233 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxyz = cbuffer.data(ig_geom_01_off + 1234 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxzz = cbuffer.data(ig_geom_01_off + 1235 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyyy = cbuffer.data(ig_geom_01_off + 1236 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyyz = cbuffer.data(ig_geom_01_off + 1237 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyzz = cbuffer.data(ig_geom_01_off + 1238 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xzzz = cbuffer.data(ig_geom_01_off + 1239 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyyy = cbuffer.data(ig_geom_01_off + 1240 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyyz = cbuffer.data(ig_geom_01_off + 1241 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyzz = cbuffer.data(ig_geom_01_off + 1242 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yzzz = cbuffer.data(ig_geom_01_off + 1243 * ccomps * dcomps);

            auto g_0_z_yzzzzz_zzzz = cbuffer.data(ig_geom_01_off + 1244 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzzz_xxxx, g_0_z_yzzzzz_xxxy, g_0_z_yzzzzz_xxxz, g_0_z_yzzzzz_xxyy, g_0_z_yzzzzz_xxyz, g_0_z_yzzzzz_xxzz, g_0_z_yzzzzz_xyyy, g_0_z_yzzzzz_xyyz, g_0_z_yzzzzz_xyzz, g_0_z_yzzzzz_xzzz, g_0_z_yzzzzz_yyyy, g_0_z_yzzzzz_yyyz, g_0_z_yzzzzz_yyzz, g_0_z_yzzzzz_yzzz, g_0_z_yzzzzz_zzzz, g_0_z_zzzzz_xxxx, g_0_z_zzzzz_xxxxy, g_0_z_zzzzz_xxxy, g_0_z_zzzzz_xxxyy, g_0_z_zzzzz_xxxyz, g_0_z_zzzzz_xxxz, g_0_z_zzzzz_xxyy, g_0_z_zzzzz_xxyyy, g_0_z_zzzzz_xxyyz, g_0_z_zzzzz_xxyz, g_0_z_zzzzz_xxyzz, g_0_z_zzzzz_xxzz, g_0_z_zzzzz_xyyy, g_0_z_zzzzz_xyyyy, g_0_z_zzzzz_xyyyz, g_0_z_zzzzz_xyyz, g_0_z_zzzzz_xyyzz, g_0_z_zzzzz_xyzz, g_0_z_zzzzz_xyzzz, g_0_z_zzzzz_xzzz, g_0_z_zzzzz_yyyy, g_0_z_zzzzz_yyyyy, g_0_z_zzzzz_yyyyz, g_0_z_zzzzz_yyyz, g_0_z_zzzzz_yyyzz, g_0_z_zzzzz_yyzz, g_0_z_zzzzz_yyzzz, g_0_z_zzzzz_yzzz, g_0_z_zzzzz_yzzzz, g_0_z_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzzz_xxxx[k] = -g_0_z_zzzzz_xxxx[k] * ab_y + g_0_z_zzzzz_xxxxy[k];

                g_0_z_yzzzzz_xxxy[k] = -g_0_z_zzzzz_xxxy[k] * ab_y + g_0_z_zzzzz_xxxyy[k];

                g_0_z_yzzzzz_xxxz[k] = -g_0_z_zzzzz_xxxz[k] * ab_y + g_0_z_zzzzz_xxxyz[k];

                g_0_z_yzzzzz_xxyy[k] = -g_0_z_zzzzz_xxyy[k] * ab_y + g_0_z_zzzzz_xxyyy[k];

                g_0_z_yzzzzz_xxyz[k] = -g_0_z_zzzzz_xxyz[k] * ab_y + g_0_z_zzzzz_xxyyz[k];

                g_0_z_yzzzzz_xxzz[k] = -g_0_z_zzzzz_xxzz[k] * ab_y + g_0_z_zzzzz_xxyzz[k];

                g_0_z_yzzzzz_xyyy[k] = -g_0_z_zzzzz_xyyy[k] * ab_y + g_0_z_zzzzz_xyyyy[k];

                g_0_z_yzzzzz_xyyz[k] = -g_0_z_zzzzz_xyyz[k] * ab_y + g_0_z_zzzzz_xyyyz[k];

                g_0_z_yzzzzz_xyzz[k] = -g_0_z_zzzzz_xyzz[k] * ab_y + g_0_z_zzzzz_xyyzz[k];

                g_0_z_yzzzzz_xzzz[k] = -g_0_z_zzzzz_xzzz[k] * ab_y + g_0_z_zzzzz_xyzzz[k];

                g_0_z_yzzzzz_yyyy[k] = -g_0_z_zzzzz_yyyy[k] * ab_y + g_0_z_zzzzz_yyyyy[k];

                g_0_z_yzzzzz_yyyz[k] = -g_0_z_zzzzz_yyyz[k] * ab_y + g_0_z_zzzzz_yyyyz[k];

                g_0_z_yzzzzz_yyzz[k] = -g_0_z_zzzzz_yyzz[k] * ab_y + g_0_z_zzzzz_yyyzz[k];

                g_0_z_yzzzzz_yzzz[k] = -g_0_z_zzzzz_yzzz[k] * ab_y + g_0_z_zzzzz_yyzzz[k];

                g_0_z_yzzzzz_zzzz[k] = -g_0_z_zzzzz_zzzz[k] * ab_y + g_0_z_zzzzz_yzzzz[k];
            }

            /// Set up 1245-1260 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzzz_xxxx = cbuffer.data(ig_geom_01_off + 1245 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxxy = cbuffer.data(ig_geom_01_off + 1246 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxxz = cbuffer.data(ig_geom_01_off + 1247 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxyy = cbuffer.data(ig_geom_01_off + 1248 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxyz = cbuffer.data(ig_geom_01_off + 1249 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxzz = cbuffer.data(ig_geom_01_off + 1250 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyyy = cbuffer.data(ig_geom_01_off + 1251 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyyz = cbuffer.data(ig_geom_01_off + 1252 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyzz = cbuffer.data(ig_geom_01_off + 1253 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xzzz = cbuffer.data(ig_geom_01_off + 1254 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyyy = cbuffer.data(ig_geom_01_off + 1255 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyyz = cbuffer.data(ig_geom_01_off + 1256 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyzz = cbuffer.data(ig_geom_01_off + 1257 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yzzz = cbuffer.data(ig_geom_01_off + 1258 * ccomps * dcomps);

            auto g_0_z_zzzzzz_zzzz = cbuffer.data(ig_geom_01_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzzz_xxxx, g_0_z_zzzzz_xxxxz, g_0_z_zzzzz_xxxy, g_0_z_zzzzz_xxxyz, g_0_z_zzzzz_xxxz, g_0_z_zzzzz_xxxzz, g_0_z_zzzzz_xxyy, g_0_z_zzzzz_xxyyz, g_0_z_zzzzz_xxyz, g_0_z_zzzzz_xxyzz, g_0_z_zzzzz_xxzz, g_0_z_zzzzz_xxzzz, g_0_z_zzzzz_xyyy, g_0_z_zzzzz_xyyyz, g_0_z_zzzzz_xyyz, g_0_z_zzzzz_xyyzz, g_0_z_zzzzz_xyzz, g_0_z_zzzzz_xyzzz, g_0_z_zzzzz_xzzz, g_0_z_zzzzz_xzzzz, g_0_z_zzzzz_yyyy, g_0_z_zzzzz_yyyyz, g_0_z_zzzzz_yyyz, g_0_z_zzzzz_yyyzz, g_0_z_zzzzz_yyzz, g_0_z_zzzzz_yyzzz, g_0_z_zzzzz_yzzz, g_0_z_zzzzz_yzzzz, g_0_z_zzzzz_zzzz, g_0_z_zzzzz_zzzzz, g_0_z_zzzzzz_xxxx, g_0_z_zzzzzz_xxxy, g_0_z_zzzzzz_xxxz, g_0_z_zzzzzz_xxyy, g_0_z_zzzzzz_xxyz, g_0_z_zzzzzz_xxzz, g_0_z_zzzzzz_xyyy, g_0_z_zzzzzz_xyyz, g_0_z_zzzzzz_xyzz, g_0_z_zzzzzz_xzzz, g_0_z_zzzzzz_yyyy, g_0_z_zzzzzz_yyyz, g_0_z_zzzzzz_yyzz, g_0_z_zzzzzz_yzzz, g_0_z_zzzzzz_zzzz, g_zzzzz_xxxx, g_zzzzz_xxxy, g_zzzzz_xxxz, g_zzzzz_xxyy, g_zzzzz_xxyz, g_zzzzz_xxzz, g_zzzzz_xyyy, g_zzzzz_xyyz, g_zzzzz_xyzz, g_zzzzz_xzzz, g_zzzzz_yyyy, g_zzzzz_yyyz, g_zzzzz_yyzz, g_zzzzz_yzzz, g_zzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzzz_xxxx[k] = g_zzzzz_xxxx[k] - g_0_z_zzzzz_xxxx[k] * ab_z + g_0_z_zzzzz_xxxxz[k];

                g_0_z_zzzzzz_xxxy[k] = g_zzzzz_xxxy[k] - g_0_z_zzzzz_xxxy[k] * ab_z + g_0_z_zzzzz_xxxyz[k];

                g_0_z_zzzzzz_xxxz[k] = g_zzzzz_xxxz[k] - g_0_z_zzzzz_xxxz[k] * ab_z + g_0_z_zzzzz_xxxzz[k];

                g_0_z_zzzzzz_xxyy[k] = g_zzzzz_xxyy[k] - g_0_z_zzzzz_xxyy[k] * ab_z + g_0_z_zzzzz_xxyyz[k];

                g_0_z_zzzzzz_xxyz[k] = g_zzzzz_xxyz[k] - g_0_z_zzzzz_xxyz[k] * ab_z + g_0_z_zzzzz_xxyzz[k];

                g_0_z_zzzzzz_xxzz[k] = g_zzzzz_xxzz[k] - g_0_z_zzzzz_xxzz[k] * ab_z + g_0_z_zzzzz_xxzzz[k];

                g_0_z_zzzzzz_xyyy[k] = g_zzzzz_xyyy[k] - g_0_z_zzzzz_xyyy[k] * ab_z + g_0_z_zzzzz_xyyyz[k];

                g_0_z_zzzzzz_xyyz[k] = g_zzzzz_xyyz[k] - g_0_z_zzzzz_xyyz[k] * ab_z + g_0_z_zzzzz_xyyzz[k];

                g_0_z_zzzzzz_xyzz[k] = g_zzzzz_xyzz[k] - g_0_z_zzzzz_xyzz[k] * ab_z + g_0_z_zzzzz_xyzzz[k];

                g_0_z_zzzzzz_xzzz[k] = g_zzzzz_xzzz[k] - g_0_z_zzzzz_xzzz[k] * ab_z + g_0_z_zzzzz_xzzzz[k];

                g_0_z_zzzzzz_yyyy[k] = g_zzzzz_yyyy[k] - g_0_z_zzzzz_yyyy[k] * ab_z + g_0_z_zzzzz_yyyyz[k];

                g_0_z_zzzzzz_yyyz[k] = g_zzzzz_yyyz[k] - g_0_z_zzzzz_yyyz[k] * ab_z + g_0_z_zzzzz_yyyzz[k];

                g_0_z_zzzzzz_yyzz[k] = g_zzzzz_yyzz[k] - g_0_z_zzzzz_yyzz[k] * ab_z + g_0_z_zzzzz_yyzzz[k];

                g_0_z_zzzzzz_yzzz[k] = g_zzzzz_yzzz[k] - g_0_z_zzzzz_yzzz[k] * ab_z + g_0_z_zzzzz_yzzzz[k];

                g_0_z_zzzzzz_zzzz[k] = g_zzzzz_zzzz[k] - g_0_z_zzzzz_zzzz[k] * ab_z + g_0_z_zzzzz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

