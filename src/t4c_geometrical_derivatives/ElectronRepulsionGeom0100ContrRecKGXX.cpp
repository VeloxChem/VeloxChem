#include "ElectronRepulsionGeom0100ContrRecKGXX.hpp"

#include "TensorComponents.hpp"

namespace erirec { // erirec namespace

auto
comp_bra_geom01_hrr_electron_repulsion_kgxx(CSimdArray<double>& cbuffer,
                                            const size_t idx_geom_01_kgxx,
                                            const size_t idx_igxx,
                                            const size_t idx_geom_01_igxx,
                                            const size_t idx_geom_01_ihxx,
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
            /// Set up components of auxilary buffer : IGSS

            const auto ig_off = idx_igxx + i * dcomps + j;

            auto g_xxxxxx_xxxx = cbuffer.data(ig_off + 0 * ccomps * dcomps);

            auto g_xxxxxx_xxxy = cbuffer.data(ig_off + 1 * ccomps * dcomps);

            auto g_xxxxxx_xxxz = cbuffer.data(ig_off + 2 * ccomps * dcomps);

            auto g_xxxxxx_xxyy = cbuffer.data(ig_off + 3 * ccomps * dcomps);

            auto g_xxxxxx_xxyz = cbuffer.data(ig_off + 4 * ccomps * dcomps);

            auto g_xxxxxx_xxzz = cbuffer.data(ig_off + 5 * ccomps * dcomps);

            auto g_xxxxxx_xyyy = cbuffer.data(ig_off + 6 * ccomps * dcomps);

            auto g_xxxxxx_xyyz = cbuffer.data(ig_off + 7 * ccomps * dcomps);

            auto g_xxxxxx_xyzz = cbuffer.data(ig_off + 8 * ccomps * dcomps);

            auto g_xxxxxx_xzzz = cbuffer.data(ig_off + 9 * ccomps * dcomps);

            auto g_xxxxxx_yyyy = cbuffer.data(ig_off + 10 * ccomps * dcomps);

            auto g_xxxxxx_yyyz = cbuffer.data(ig_off + 11 * ccomps * dcomps);

            auto g_xxxxxx_yyzz = cbuffer.data(ig_off + 12 * ccomps * dcomps);

            auto g_xxxxxx_yzzz = cbuffer.data(ig_off + 13 * ccomps * dcomps);

            auto g_xxxxxx_zzzz = cbuffer.data(ig_off + 14 * ccomps * dcomps);

            auto g_xxxxxy_xxxx = cbuffer.data(ig_off + 15 * ccomps * dcomps);

            auto g_xxxxxy_xxxy = cbuffer.data(ig_off + 16 * ccomps * dcomps);

            auto g_xxxxxy_xxxz = cbuffer.data(ig_off + 17 * ccomps * dcomps);

            auto g_xxxxxy_xxyy = cbuffer.data(ig_off + 18 * ccomps * dcomps);

            auto g_xxxxxy_xxyz = cbuffer.data(ig_off + 19 * ccomps * dcomps);

            auto g_xxxxxy_xxzz = cbuffer.data(ig_off + 20 * ccomps * dcomps);

            auto g_xxxxxy_xyyy = cbuffer.data(ig_off + 21 * ccomps * dcomps);

            auto g_xxxxxy_xyyz = cbuffer.data(ig_off + 22 * ccomps * dcomps);

            auto g_xxxxxy_xyzz = cbuffer.data(ig_off + 23 * ccomps * dcomps);

            auto g_xxxxxy_xzzz = cbuffer.data(ig_off + 24 * ccomps * dcomps);

            auto g_xxxxxy_yyyy = cbuffer.data(ig_off + 25 * ccomps * dcomps);

            auto g_xxxxxy_yyyz = cbuffer.data(ig_off + 26 * ccomps * dcomps);

            auto g_xxxxxy_yyzz = cbuffer.data(ig_off + 27 * ccomps * dcomps);

            auto g_xxxxxy_yzzz = cbuffer.data(ig_off + 28 * ccomps * dcomps);

            auto g_xxxxxy_zzzz = cbuffer.data(ig_off + 29 * ccomps * dcomps);

            auto g_xxxxxz_xxxx = cbuffer.data(ig_off + 30 * ccomps * dcomps);

            auto g_xxxxxz_xxxy = cbuffer.data(ig_off + 31 * ccomps * dcomps);

            auto g_xxxxxz_xxxz = cbuffer.data(ig_off + 32 * ccomps * dcomps);

            auto g_xxxxxz_xxyy = cbuffer.data(ig_off + 33 * ccomps * dcomps);

            auto g_xxxxxz_xxyz = cbuffer.data(ig_off + 34 * ccomps * dcomps);

            auto g_xxxxxz_xxzz = cbuffer.data(ig_off + 35 * ccomps * dcomps);

            auto g_xxxxxz_xyyy = cbuffer.data(ig_off + 36 * ccomps * dcomps);

            auto g_xxxxxz_xyyz = cbuffer.data(ig_off + 37 * ccomps * dcomps);

            auto g_xxxxxz_xyzz = cbuffer.data(ig_off + 38 * ccomps * dcomps);

            auto g_xxxxxz_xzzz = cbuffer.data(ig_off + 39 * ccomps * dcomps);

            auto g_xxxxxz_yyyy = cbuffer.data(ig_off + 40 * ccomps * dcomps);

            auto g_xxxxxz_yyyz = cbuffer.data(ig_off + 41 * ccomps * dcomps);

            auto g_xxxxxz_yyzz = cbuffer.data(ig_off + 42 * ccomps * dcomps);

            auto g_xxxxxz_yzzz = cbuffer.data(ig_off + 43 * ccomps * dcomps);

            auto g_xxxxxz_zzzz = cbuffer.data(ig_off + 44 * ccomps * dcomps);

            auto g_xxxxyy_xxxx = cbuffer.data(ig_off + 45 * ccomps * dcomps);

            auto g_xxxxyy_xxxy = cbuffer.data(ig_off + 46 * ccomps * dcomps);

            auto g_xxxxyy_xxxz = cbuffer.data(ig_off + 47 * ccomps * dcomps);

            auto g_xxxxyy_xxyy = cbuffer.data(ig_off + 48 * ccomps * dcomps);

            auto g_xxxxyy_xxyz = cbuffer.data(ig_off + 49 * ccomps * dcomps);

            auto g_xxxxyy_xxzz = cbuffer.data(ig_off + 50 * ccomps * dcomps);

            auto g_xxxxyy_xyyy = cbuffer.data(ig_off + 51 * ccomps * dcomps);

            auto g_xxxxyy_xyyz = cbuffer.data(ig_off + 52 * ccomps * dcomps);

            auto g_xxxxyy_xyzz = cbuffer.data(ig_off + 53 * ccomps * dcomps);

            auto g_xxxxyy_xzzz = cbuffer.data(ig_off + 54 * ccomps * dcomps);

            auto g_xxxxyy_yyyy = cbuffer.data(ig_off + 55 * ccomps * dcomps);

            auto g_xxxxyy_yyyz = cbuffer.data(ig_off + 56 * ccomps * dcomps);

            auto g_xxxxyy_yyzz = cbuffer.data(ig_off + 57 * ccomps * dcomps);

            auto g_xxxxyy_yzzz = cbuffer.data(ig_off + 58 * ccomps * dcomps);

            auto g_xxxxyy_zzzz = cbuffer.data(ig_off + 59 * ccomps * dcomps);

            auto g_xxxxyz_xxxx = cbuffer.data(ig_off + 60 * ccomps * dcomps);

            auto g_xxxxyz_xxxy = cbuffer.data(ig_off + 61 * ccomps * dcomps);

            auto g_xxxxyz_xxxz = cbuffer.data(ig_off + 62 * ccomps * dcomps);

            auto g_xxxxyz_xxyy = cbuffer.data(ig_off + 63 * ccomps * dcomps);

            auto g_xxxxyz_xxyz = cbuffer.data(ig_off + 64 * ccomps * dcomps);

            auto g_xxxxyz_xxzz = cbuffer.data(ig_off + 65 * ccomps * dcomps);

            auto g_xxxxyz_xyyy = cbuffer.data(ig_off + 66 * ccomps * dcomps);

            auto g_xxxxyz_xyyz = cbuffer.data(ig_off + 67 * ccomps * dcomps);

            auto g_xxxxyz_xyzz = cbuffer.data(ig_off + 68 * ccomps * dcomps);

            auto g_xxxxyz_xzzz = cbuffer.data(ig_off + 69 * ccomps * dcomps);

            auto g_xxxxyz_yyyy = cbuffer.data(ig_off + 70 * ccomps * dcomps);

            auto g_xxxxyz_yyyz = cbuffer.data(ig_off + 71 * ccomps * dcomps);

            auto g_xxxxyz_yyzz = cbuffer.data(ig_off + 72 * ccomps * dcomps);

            auto g_xxxxyz_yzzz = cbuffer.data(ig_off + 73 * ccomps * dcomps);

            auto g_xxxxyz_zzzz = cbuffer.data(ig_off + 74 * ccomps * dcomps);

            auto g_xxxxzz_xxxx = cbuffer.data(ig_off + 75 * ccomps * dcomps);

            auto g_xxxxzz_xxxy = cbuffer.data(ig_off + 76 * ccomps * dcomps);

            auto g_xxxxzz_xxxz = cbuffer.data(ig_off + 77 * ccomps * dcomps);

            auto g_xxxxzz_xxyy = cbuffer.data(ig_off + 78 * ccomps * dcomps);

            auto g_xxxxzz_xxyz = cbuffer.data(ig_off + 79 * ccomps * dcomps);

            auto g_xxxxzz_xxzz = cbuffer.data(ig_off + 80 * ccomps * dcomps);

            auto g_xxxxzz_xyyy = cbuffer.data(ig_off + 81 * ccomps * dcomps);

            auto g_xxxxzz_xyyz = cbuffer.data(ig_off + 82 * ccomps * dcomps);

            auto g_xxxxzz_xyzz = cbuffer.data(ig_off + 83 * ccomps * dcomps);

            auto g_xxxxzz_xzzz = cbuffer.data(ig_off + 84 * ccomps * dcomps);

            auto g_xxxxzz_yyyy = cbuffer.data(ig_off + 85 * ccomps * dcomps);

            auto g_xxxxzz_yyyz = cbuffer.data(ig_off + 86 * ccomps * dcomps);

            auto g_xxxxzz_yyzz = cbuffer.data(ig_off + 87 * ccomps * dcomps);

            auto g_xxxxzz_yzzz = cbuffer.data(ig_off + 88 * ccomps * dcomps);

            auto g_xxxxzz_zzzz = cbuffer.data(ig_off + 89 * ccomps * dcomps);

            auto g_xxxyyy_xxxx = cbuffer.data(ig_off + 90 * ccomps * dcomps);

            auto g_xxxyyy_xxxy = cbuffer.data(ig_off + 91 * ccomps * dcomps);

            auto g_xxxyyy_xxxz = cbuffer.data(ig_off + 92 * ccomps * dcomps);

            auto g_xxxyyy_xxyy = cbuffer.data(ig_off + 93 * ccomps * dcomps);

            auto g_xxxyyy_xxyz = cbuffer.data(ig_off + 94 * ccomps * dcomps);

            auto g_xxxyyy_xxzz = cbuffer.data(ig_off + 95 * ccomps * dcomps);

            auto g_xxxyyy_xyyy = cbuffer.data(ig_off + 96 * ccomps * dcomps);

            auto g_xxxyyy_xyyz = cbuffer.data(ig_off + 97 * ccomps * dcomps);

            auto g_xxxyyy_xyzz = cbuffer.data(ig_off + 98 * ccomps * dcomps);

            auto g_xxxyyy_xzzz = cbuffer.data(ig_off + 99 * ccomps * dcomps);

            auto g_xxxyyy_yyyy = cbuffer.data(ig_off + 100 * ccomps * dcomps);

            auto g_xxxyyy_yyyz = cbuffer.data(ig_off + 101 * ccomps * dcomps);

            auto g_xxxyyy_yyzz = cbuffer.data(ig_off + 102 * ccomps * dcomps);

            auto g_xxxyyy_yzzz = cbuffer.data(ig_off + 103 * ccomps * dcomps);

            auto g_xxxyyy_zzzz = cbuffer.data(ig_off + 104 * ccomps * dcomps);

            auto g_xxxyyz_xxxx = cbuffer.data(ig_off + 105 * ccomps * dcomps);

            auto g_xxxyyz_xxxy = cbuffer.data(ig_off + 106 * ccomps * dcomps);

            auto g_xxxyyz_xxxz = cbuffer.data(ig_off + 107 * ccomps * dcomps);

            auto g_xxxyyz_xxyy = cbuffer.data(ig_off + 108 * ccomps * dcomps);

            auto g_xxxyyz_xxyz = cbuffer.data(ig_off + 109 * ccomps * dcomps);

            auto g_xxxyyz_xxzz = cbuffer.data(ig_off + 110 * ccomps * dcomps);

            auto g_xxxyyz_xyyy = cbuffer.data(ig_off + 111 * ccomps * dcomps);

            auto g_xxxyyz_xyyz = cbuffer.data(ig_off + 112 * ccomps * dcomps);

            auto g_xxxyyz_xyzz = cbuffer.data(ig_off + 113 * ccomps * dcomps);

            auto g_xxxyyz_xzzz = cbuffer.data(ig_off + 114 * ccomps * dcomps);

            auto g_xxxyyz_yyyy = cbuffer.data(ig_off + 115 * ccomps * dcomps);

            auto g_xxxyyz_yyyz = cbuffer.data(ig_off + 116 * ccomps * dcomps);

            auto g_xxxyyz_yyzz = cbuffer.data(ig_off + 117 * ccomps * dcomps);

            auto g_xxxyyz_yzzz = cbuffer.data(ig_off + 118 * ccomps * dcomps);

            auto g_xxxyyz_zzzz = cbuffer.data(ig_off + 119 * ccomps * dcomps);

            auto g_xxxyzz_xxxx = cbuffer.data(ig_off + 120 * ccomps * dcomps);

            auto g_xxxyzz_xxxy = cbuffer.data(ig_off + 121 * ccomps * dcomps);

            auto g_xxxyzz_xxxz = cbuffer.data(ig_off + 122 * ccomps * dcomps);

            auto g_xxxyzz_xxyy = cbuffer.data(ig_off + 123 * ccomps * dcomps);

            auto g_xxxyzz_xxyz = cbuffer.data(ig_off + 124 * ccomps * dcomps);

            auto g_xxxyzz_xxzz = cbuffer.data(ig_off + 125 * ccomps * dcomps);

            auto g_xxxyzz_xyyy = cbuffer.data(ig_off + 126 * ccomps * dcomps);

            auto g_xxxyzz_xyyz = cbuffer.data(ig_off + 127 * ccomps * dcomps);

            auto g_xxxyzz_xyzz = cbuffer.data(ig_off + 128 * ccomps * dcomps);

            auto g_xxxyzz_xzzz = cbuffer.data(ig_off + 129 * ccomps * dcomps);

            auto g_xxxyzz_yyyy = cbuffer.data(ig_off + 130 * ccomps * dcomps);

            auto g_xxxyzz_yyyz = cbuffer.data(ig_off + 131 * ccomps * dcomps);

            auto g_xxxyzz_yyzz = cbuffer.data(ig_off + 132 * ccomps * dcomps);

            auto g_xxxyzz_yzzz = cbuffer.data(ig_off + 133 * ccomps * dcomps);

            auto g_xxxyzz_zzzz = cbuffer.data(ig_off + 134 * ccomps * dcomps);

            auto g_xxxzzz_xxxx = cbuffer.data(ig_off + 135 * ccomps * dcomps);

            auto g_xxxzzz_xxxy = cbuffer.data(ig_off + 136 * ccomps * dcomps);

            auto g_xxxzzz_xxxz = cbuffer.data(ig_off + 137 * ccomps * dcomps);

            auto g_xxxzzz_xxyy = cbuffer.data(ig_off + 138 * ccomps * dcomps);

            auto g_xxxzzz_xxyz = cbuffer.data(ig_off + 139 * ccomps * dcomps);

            auto g_xxxzzz_xxzz = cbuffer.data(ig_off + 140 * ccomps * dcomps);

            auto g_xxxzzz_xyyy = cbuffer.data(ig_off + 141 * ccomps * dcomps);

            auto g_xxxzzz_xyyz = cbuffer.data(ig_off + 142 * ccomps * dcomps);

            auto g_xxxzzz_xyzz = cbuffer.data(ig_off + 143 * ccomps * dcomps);

            auto g_xxxzzz_xzzz = cbuffer.data(ig_off + 144 * ccomps * dcomps);

            auto g_xxxzzz_yyyy = cbuffer.data(ig_off + 145 * ccomps * dcomps);

            auto g_xxxzzz_yyyz = cbuffer.data(ig_off + 146 * ccomps * dcomps);

            auto g_xxxzzz_yyzz = cbuffer.data(ig_off + 147 * ccomps * dcomps);

            auto g_xxxzzz_yzzz = cbuffer.data(ig_off + 148 * ccomps * dcomps);

            auto g_xxxzzz_zzzz = cbuffer.data(ig_off + 149 * ccomps * dcomps);

            auto g_xxyyyy_xxxx = cbuffer.data(ig_off + 150 * ccomps * dcomps);

            auto g_xxyyyy_xxxy = cbuffer.data(ig_off + 151 * ccomps * dcomps);

            auto g_xxyyyy_xxxz = cbuffer.data(ig_off + 152 * ccomps * dcomps);

            auto g_xxyyyy_xxyy = cbuffer.data(ig_off + 153 * ccomps * dcomps);

            auto g_xxyyyy_xxyz = cbuffer.data(ig_off + 154 * ccomps * dcomps);

            auto g_xxyyyy_xxzz = cbuffer.data(ig_off + 155 * ccomps * dcomps);

            auto g_xxyyyy_xyyy = cbuffer.data(ig_off + 156 * ccomps * dcomps);

            auto g_xxyyyy_xyyz = cbuffer.data(ig_off + 157 * ccomps * dcomps);

            auto g_xxyyyy_xyzz = cbuffer.data(ig_off + 158 * ccomps * dcomps);

            auto g_xxyyyy_xzzz = cbuffer.data(ig_off + 159 * ccomps * dcomps);

            auto g_xxyyyy_yyyy = cbuffer.data(ig_off + 160 * ccomps * dcomps);

            auto g_xxyyyy_yyyz = cbuffer.data(ig_off + 161 * ccomps * dcomps);

            auto g_xxyyyy_yyzz = cbuffer.data(ig_off + 162 * ccomps * dcomps);

            auto g_xxyyyy_yzzz = cbuffer.data(ig_off + 163 * ccomps * dcomps);

            auto g_xxyyyy_zzzz = cbuffer.data(ig_off + 164 * ccomps * dcomps);

            auto g_xxyyyz_xxxx = cbuffer.data(ig_off + 165 * ccomps * dcomps);

            auto g_xxyyyz_xxxy = cbuffer.data(ig_off + 166 * ccomps * dcomps);

            auto g_xxyyyz_xxxz = cbuffer.data(ig_off + 167 * ccomps * dcomps);

            auto g_xxyyyz_xxyy = cbuffer.data(ig_off + 168 * ccomps * dcomps);

            auto g_xxyyyz_xxyz = cbuffer.data(ig_off + 169 * ccomps * dcomps);

            auto g_xxyyyz_xxzz = cbuffer.data(ig_off + 170 * ccomps * dcomps);

            auto g_xxyyyz_xyyy = cbuffer.data(ig_off + 171 * ccomps * dcomps);

            auto g_xxyyyz_xyyz = cbuffer.data(ig_off + 172 * ccomps * dcomps);

            auto g_xxyyyz_xyzz = cbuffer.data(ig_off + 173 * ccomps * dcomps);

            auto g_xxyyyz_xzzz = cbuffer.data(ig_off + 174 * ccomps * dcomps);

            auto g_xxyyyz_yyyy = cbuffer.data(ig_off + 175 * ccomps * dcomps);

            auto g_xxyyyz_yyyz = cbuffer.data(ig_off + 176 * ccomps * dcomps);

            auto g_xxyyyz_yyzz = cbuffer.data(ig_off + 177 * ccomps * dcomps);

            auto g_xxyyyz_yzzz = cbuffer.data(ig_off + 178 * ccomps * dcomps);

            auto g_xxyyyz_zzzz = cbuffer.data(ig_off + 179 * ccomps * dcomps);

            auto g_xxyyzz_xxxx = cbuffer.data(ig_off + 180 * ccomps * dcomps);

            auto g_xxyyzz_xxxy = cbuffer.data(ig_off + 181 * ccomps * dcomps);

            auto g_xxyyzz_xxxz = cbuffer.data(ig_off + 182 * ccomps * dcomps);

            auto g_xxyyzz_xxyy = cbuffer.data(ig_off + 183 * ccomps * dcomps);

            auto g_xxyyzz_xxyz = cbuffer.data(ig_off + 184 * ccomps * dcomps);

            auto g_xxyyzz_xxzz = cbuffer.data(ig_off + 185 * ccomps * dcomps);

            auto g_xxyyzz_xyyy = cbuffer.data(ig_off + 186 * ccomps * dcomps);

            auto g_xxyyzz_xyyz = cbuffer.data(ig_off + 187 * ccomps * dcomps);

            auto g_xxyyzz_xyzz = cbuffer.data(ig_off + 188 * ccomps * dcomps);

            auto g_xxyyzz_xzzz = cbuffer.data(ig_off + 189 * ccomps * dcomps);

            auto g_xxyyzz_yyyy = cbuffer.data(ig_off + 190 * ccomps * dcomps);

            auto g_xxyyzz_yyyz = cbuffer.data(ig_off + 191 * ccomps * dcomps);

            auto g_xxyyzz_yyzz = cbuffer.data(ig_off + 192 * ccomps * dcomps);

            auto g_xxyyzz_yzzz = cbuffer.data(ig_off + 193 * ccomps * dcomps);

            auto g_xxyyzz_zzzz = cbuffer.data(ig_off + 194 * ccomps * dcomps);

            auto g_xxyzzz_xxxx = cbuffer.data(ig_off + 195 * ccomps * dcomps);

            auto g_xxyzzz_xxxy = cbuffer.data(ig_off + 196 * ccomps * dcomps);

            auto g_xxyzzz_xxxz = cbuffer.data(ig_off + 197 * ccomps * dcomps);

            auto g_xxyzzz_xxyy = cbuffer.data(ig_off + 198 * ccomps * dcomps);

            auto g_xxyzzz_xxyz = cbuffer.data(ig_off + 199 * ccomps * dcomps);

            auto g_xxyzzz_xxzz = cbuffer.data(ig_off + 200 * ccomps * dcomps);

            auto g_xxyzzz_xyyy = cbuffer.data(ig_off + 201 * ccomps * dcomps);

            auto g_xxyzzz_xyyz = cbuffer.data(ig_off + 202 * ccomps * dcomps);

            auto g_xxyzzz_xyzz = cbuffer.data(ig_off + 203 * ccomps * dcomps);

            auto g_xxyzzz_xzzz = cbuffer.data(ig_off + 204 * ccomps * dcomps);

            auto g_xxyzzz_yyyy = cbuffer.data(ig_off + 205 * ccomps * dcomps);

            auto g_xxyzzz_yyyz = cbuffer.data(ig_off + 206 * ccomps * dcomps);

            auto g_xxyzzz_yyzz = cbuffer.data(ig_off + 207 * ccomps * dcomps);

            auto g_xxyzzz_yzzz = cbuffer.data(ig_off + 208 * ccomps * dcomps);

            auto g_xxyzzz_zzzz = cbuffer.data(ig_off + 209 * ccomps * dcomps);

            auto g_xxzzzz_xxxx = cbuffer.data(ig_off + 210 * ccomps * dcomps);

            auto g_xxzzzz_xxxy = cbuffer.data(ig_off + 211 * ccomps * dcomps);

            auto g_xxzzzz_xxxz = cbuffer.data(ig_off + 212 * ccomps * dcomps);

            auto g_xxzzzz_xxyy = cbuffer.data(ig_off + 213 * ccomps * dcomps);

            auto g_xxzzzz_xxyz = cbuffer.data(ig_off + 214 * ccomps * dcomps);

            auto g_xxzzzz_xxzz = cbuffer.data(ig_off + 215 * ccomps * dcomps);

            auto g_xxzzzz_xyyy = cbuffer.data(ig_off + 216 * ccomps * dcomps);

            auto g_xxzzzz_xyyz = cbuffer.data(ig_off + 217 * ccomps * dcomps);

            auto g_xxzzzz_xyzz = cbuffer.data(ig_off + 218 * ccomps * dcomps);

            auto g_xxzzzz_xzzz = cbuffer.data(ig_off + 219 * ccomps * dcomps);

            auto g_xxzzzz_yyyy = cbuffer.data(ig_off + 220 * ccomps * dcomps);

            auto g_xxzzzz_yyyz = cbuffer.data(ig_off + 221 * ccomps * dcomps);

            auto g_xxzzzz_yyzz = cbuffer.data(ig_off + 222 * ccomps * dcomps);

            auto g_xxzzzz_yzzz = cbuffer.data(ig_off + 223 * ccomps * dcomps);

            auto g_xxzzzz_zzzz = cbuffer.data(ig_off + 224 * ccomps * dcomps);

            auto g_xyyyyy_xxxx = cbuffer.data(ig_off + 225 * ccomps * dcomps);

            auto g_xyyyyy_xxxy = cbuffer.data(ig_off + 226 * ccomps * dcomps);

            auto g_xyyyyy_xxxz = cbuffer.data(ig_off + 227 * ccomps * dcomps);

            auto g_xyyyyy_xxyy = cbuffer.data(ig_off + 228 * ccomps * dcomps);

            auto g_xyyyyy_xxyz = cbuffer.data(ig_off + 229 * ccomps * dcomps);

            auto g_xyyyyy_xxzz = cbuffer.data(ig_off + 230 * ccomps * dcomps);

            auto g_xyyyyy_xyyy = cbuffer.data(ig_off + 231 * ccomps * dcomps);

            auto g_xyyyyy_xyyz = cbuffer.data(ig_off + 232 * ccomps * dcomps);

            auto g_xyyyyy_xyzz = cbuffer.data(ig_off + 233 * ccomps * dcomps);

            auto g_xyyyyy_xzzz = cbuffer.data(ig_off + 234 * ccomps * dcomps);

            auto g_xyyyyy_yyyy = cbuffer.data(ig_off + 235 * ccomps * dcomps);

            auto g_xyyyyy_yyyz = cbuffer.data(ig_off + 236 * ccomps * dcomps);

            auto g_xyyyyy_yyzz = cbuffer.data(ig_off + 237 * ccomps * dcomps);

            auto g_xyyyyy_yzzz = cbuffer.data(ig_off + 238 * ccomps * dcomps);

            auto g_xyyyyy_zzzz = cbuffer.data(ig_off + 239 * ccomps * dcomps);

            auto g_xyyyyz_xxxx = cbuffer.data(ig_off + 240 * ccomps * dcomps);

            auto g_xyyyyz_xxxy = cbuffer.data(ig_off + 241 * ccomps * dcomps);

            auto g_xyyyyz_xxxz = cbuffer.data(ig_off + 242 * ccomps * dcomps);

            auto g_xyyyyz_xxyy = cbuffer.data(ig_off + 243 * ccomps * dcomps);

            auto g_xyyyyz_xxyz = cbuffer.data(ig_off + 244 * ccomps * dcomps);

            auto g_xyyyyz_xxzz = cbuffer.data(ig_off + 245 * ccomps * dcomps);

            auto g_xyyyyz_xyyy = cbuffer.data(ig_off + 246 * ccomps * dcomps);

            auto g_xyyyyz_xyyz = cbuffer.data(ig_off + 247 * ccomps * dcomps);

            auto g_xyyyyz_xyzz = cbuffer.data(ig_off + 248 * ccomps * dcomps);

            auto g_xyyyyz_xzzz = cbuffer.data(ig_off + 249 * ccomps * dcomps);

            auto g_xyyyyz_yyyy = cbuffer.data(ig_off + 250 * ccomps * dcomps);

            auto g_xyyyyz_yyyz = cbuffer.data(ig_off + 251 * ccomps * dcomps);

            auto g_xyyyyz_yyzz = cbuffer.data(ig_off + 252 * ccomps * dcomps);

            auto g_xyyyyz_yzzz = cbuffer.data(ig_off + 253 * ccomps * dcomps);

            auto g_xyyyyz_zzzz = cbuffer.data(ig_off + 254 * ccomps * dcomps);

            auto g_xyyyzz_xxxx = cbuffer.data(ig_off + 255 * ccomps * dcomps);

            auto g_xyyyzz_xxxy = cbuffer.data(ig_off + 256 * ccomps * dcomps);

            auto g_xyyyzz_xxxz = cbuffer.data(ig_off + 257 * ccomps * dcomps);

            auto g_xyyyzz_xxyy = cbuffer.data(ig_off + 258 * ccomps * dcomps);

            auto g_xyyyzz_xxyz = cbuffer.data(ig_off + 259 * ccomps * dcomps);

            auto g_xyyyzz_xxzz = cbuffer.data(ig_off + 260 * ccomps * dcomps);

            auto g_xyyyzz_xyyy = cbuffer.data(ig_off + 261 * ccomps * dcomps);

            auto g_xyyyzz_xyyz = cbuffer.data(ig_off + 262 * ccomps * dcomps);

            auto g_xyyyzz_xyzz = cbuffer.data(ig_off + 263 * ccomps * dcomps);

            auto g_xyyyzz_xzzz = cbuffer.data(ig_off + 264 * ccomps * dcomps);

            auto g_xyyyzz_yyyy = cbuffer.data(ig_off + 265 * ccomps * dcomps);

            auto g_xyyyzz_yyyz = cbuffer.data(ig_off + 266 * ccomps * dcomps);

            auto g_xyyyzz_yyzz = cbuffer.data(ig_off + 267 * ccomps * dcomps);

            auto g_xyyyzz_yzzz = cbuffer.data(ig_off + 268 * ccomps * dcomps);

            auto g_xyyyzz_zzzz = cbuffer.data(ig_off + 269 * ccomps * dcomps);

            auto g_xyyzzz_xxxx = cbuffer.data(ig_off + 270 * ccomps * dcomps);

            auto g_xyyzzz_xxxy = cbuffer.data(ig_off + 271 * ccomps * dcomps);

            auto g_xyyzzz_xxxz = cbuffer.data(ig_off + 272 * ccomps * dcomps);

            auto g_xyyzzz_xxyy = cbuffer.data(ig_off + 273 * ccomps * dcomps);

            auto g_xyyzzz_xxyz = cbuffer.data(ig_off + 274 * ccomps * dcomps);

            auto g_xyyzzz_xxzz = cbuffer.data(ig_off + 275 * ccomps * dcomps);

            auto g_xyyzzz_xyyy = cbuffer.data(ig_off + 276 * ccomps * dcomps);

            auto g_xyyzzz_xyyz = cbuffer.data(ig_off + 277 * ccomps * dcomps);

            auto g_xyyzzz_xyzz = cbuffer.data(ig_off + 278 * ccomps * dcomps);

            auto g_xyyzzz_xzzz = cbuffer.data(ig_off + 279 * ccomps * dcomps);

            auto g_xyyzzz_yyyy = cbuffer.data(ig_off + 280 * ccomps * dcomps);

            auto g_xyyzzz_yyyz = cbuffer.data(ig_off + 281 * ccomps * dcomps);

            auto g_xyyzzz_yyzz = cbuffer.data(ig_off + 282 * ccomps * dcomps);

            auto g_xyyzzz_yzzz = cbuffer.data(ig_off + 283 * ccomps * dcomps);

            auto g_xyyzzz_zzzz = cbuffer.data(ig_off + 284 * ccomps * dcomps);

            auto g_xyzzzz_xxxx = cbuffer.data(ig_off + 285 * ccomps * dcomps);

            auto g_xyzzzz_xxxy = cbuffer.data(ig_off + 286 * ccomps * dcomps);

            auto g_xyzzzz_xxxz = cbuffer.data(ig_off + 287 * ccomps * dcomps);

            auto g_xyzzzz_xxyy = cbuffer.data(ig_off + 288 * ccomps * dcomps);

            auto g_xyzzzz_xxyz = cbuffer.data(ig_off + 289 * ccomps * dcomps);

            auto g_xyzzzz_xxzz = cbuffer.data(ig_off + 290 * ccomps * dcomps);

            auto g_xyzzzz_xyyy = cbuffer.data(ig_off + 291 * ccomps * dcomps);

            auto g_xyzzzz_xyyz = cbuffer.data(ig_off + 292 * ccomps * dcomps);

            auto g_xyzzzz_xyzz = cbuffer.data(ig_off + 293 * ccomps * dcomps);

            auto g_xyzzzz_xzzz = cbuffer.data(ig_off + 294 * ccomps * dcomps);

            auto g_xyzzzz_yyyy = cbuffer.data(ig_off + 295 * ccomps * dcomps);

            auto g_xyzzzz_yyyz = cbuffer.data(ig_off + 296 * ccomps * dcomps);

            auto g_xyzzzz_yyzz = cbuffer.data(ig_off + 297 * ccomps * dcomps);

            auto g_xyzzzz_yzzz = cbuffer.data(ig_off + 298 * ccomps * dcomps);

            auto g_xyzzzz_zzzz = cbuffer.data(ig_off + 299 * ccomps * dcomps);

            auto g_xzzzzz_xxxx = cbuffer.data(ig_off + 300 * ccomps * dcomps);

            auto g_xzzzzz_xxxy = cbuffer.data(ig_off + 301 * ccomps * dcomps);

            auto g_xzzzzz_xxxz = cbuffer.data(ig_off + 302 * ccomps * dcomps);

            auto g_xzzzzz_xxyy = cbuffer.data(ig_off + 303 * ccomps * dcomps);

            auto g_xzzzzz_xxyz = cbuffer.data(ig_off + 304 * ccomps * dcomps);

            auto g_xzzzzz_xxzz = cbuffer.data(ig_off + 305 * ccomps * dcomps);

            auto g_xzzzzz_xyyy = cbuffer.data(ig_off + 306 * ccomps * dcomps);

            auto g_xzzzzz_xyyz = cbuffer.data(ig_off + 307 * ccomps * dcomps);

            auto g_xzzzzz_xyzz = cbuffer.data(ig_off + 308 * ccomps * dcomps);

            auto g_xzzzzz_xzzz = cbuffer.data(ig_off + 309 * ccomps * dcomps);

            auto g_xzzzzz_yyyy = cbuffer.data(ig_off + 310 * ccomps * dcomps);

            auto g_xzzzzz_yyyz = cbuffer.data(ig_off + 311 * ccomps * dcomps);

            auto g_xzzzzz_yyzz = cbuffer.data(ig_off + 312 * ccomps * dcomps);

            auto g_xzzzzz_yzzz = cbuffer.data(ig_off + 313 * ccomps * dcomps);

            auto g_xzzzzz_zzzz = cbuffer.data(ig_off + 314 * ccomps * dcomps);

            auto g_yyyyyy_xxxx = cbuffer.data(ig_off + 315 * ccomps * dcomps);

            auto g_yyyyyy_xxxy = cbuffer.data(ig_off + 316 * ccomps * dcomps);

            auto g_yyyyyy_xxxz = cbuffer.data(ig_off + 317 * ccomps * dcomps);

            auto g_yyyyyy_xxyy = cbuffer.data(ig_off + 318 * ccomps * dcomps);

            auto g_yyyyyy_xxyz = cbuffer.data(ig_off + 319 * ccomps * dcomps);

            auto g_yyyyyy_xxzz = cbuffer.data(ig_off + 320 * ccomps * dcomps);

            auto g_yyyyyy_xyyy = cbuffer.data(ig_off + 321 * ccomps * dcomps);

            auto g_yyyyyy_xyyz = cbuffer.data(ig_off + 322 * ccomps * dcomps);

            auto g_yyyyyy_xyzz = cbuffer.data(ig_off + 323 * ccomps * dcomps);

            auto g_yyyyyy_xzzz = cbuffer.data(ig_off + 324 * ccomps * dcomps);

            auto g_yyyyyy_yyyy = cbuffer.data(ig_off + 325 * ccomps * dcomps);

            auto g_yyyyyy_yyyz = cbuffer.data(ig_off + 326 * ccomps * dcomps);

            auto g_yyyyyy_yyzz = cbuffer.data(ig_off + 327 * ccomps * dcomps);

            auto g_yyyyyy_yzzz = cbuffer.data(ig_off + 328 * ccomps * dcomps);

            auto g_yyyyyy_zzzz = cbuffer.data(ig_off + 329 * ccomps * dcomps);

            auto g_yyyyyz_xxxx = cbuffer.data(ig_off + 330 * ccomps * dcomps);

            auto g_yyyyyz_xxxy = cbuffer.data(ig_off + 331 * ccomps * dcomps);

            auto g_yyyyyz_xxxz = cbuffer.data(ig_off + 332 * ccomps * dcomps);

            auto g_yyyyyz_xxyy = cbuffer.data(ig_off + 333 * ccomps * dcomps);

            auto g_yyyyyz_xxyz = cbuffer.data(ig_off + 334 * ccomps * dcomps);

            auto g_yyyyyz_xxzz = cbuffer.data(ig_off + 335 * ccomps * dcomps);

            auto g_yyyyyz_xyyy = cbuffer.data(ig_off + 336 * ccomps * dcomps);

            auto g_yyyyyz_xyyz = cbuffer.data(ig_off + 337 * ccomps * dcomps);

            auto g_yyyyyz_xyzz = cbuffer.data(ig_off + 338 * ccomps * dcomps);

            auto g_yyyyyz_xzzz = cbuffer.data(ig_off + 339 * ccomps * dcomps);

            auto g_yyyyyz_yyyy = cbuffer.data(ig_off + 340 * ccomps * dcomps);

            auto g_yyyyyz_yyyz = cbuffer.data(ig_off + 341 * ccomps * dcomps);

            auto g_yyyyyz_yyzz = cbuffer.data(ig_off + 342 * ccomps * dcomps);

            auto g_yyyyyz_yzzz = cbuffer.data(ig_off + 343 * ccomps * dcomps);

            auto g_yyyyyz_zzzz = cbuffer.data(ig_off + 344 * ccomps * dcomps);

            auto g_yyyyzz_xxxx = cbuffer.data(ig_off + 345 * ccomps * dcomps);

            auto g_yyyyzz_xxxy = cbuffer.data(ig_off + 346 * ccomps * dcomps);

            auto g_yyyyzz_xxxz = cbuffer.data(ig_off + 347 * ccomps * dcomps);

            auto g_yyyyzz_xxyy = cbuffer.data(ig_off + 348 * ccomps * dcomps);

            auto g_yyyyzz_xxyz = cbuffer.data(ig_off + 349 * ccomps * dcomps);

            auto g_yyyyzz_xxzz = cbuffer.data(ig_off + 350 * ccomps * dcomps);

            auto g_yyyyzz_xyyy = cbuffer.data(ig_off + 351 * ccomps * dcomps);

            auto g_yyyyzz_xyyz = cbuffer.data(ig_off + 352 * ccomps * dcomps);

            auto g_yyyyzz_xyzz = cbuffer.data(ig_off + 353 * ccomps * dcomps);

            auto g_yyyyzz_xzzz = cbuffer.data(ig_off + 354 * ccomps * dcomps);

            auto g_yyyyzz_yyyy = cbuffer.data(ig_off + 355 * ccomps * dcomps);

            auto g_yyyyzz_yyyz = cbuffer.data(ig_off + 356 * ccomps * dcomps);

            auto g_yyyyzz_yyzz = cbuffer.data(ig_off + 357 * ccomps * dcomps);

            auto g_yyyyzz_yzzz = cbuffer.data(ig_off + 358 * ccomps * dcomps);

            auto g_yyyyzz_zzzz = cbuffer.data(ig_off + 359 * ccomps * dcomps);

            auto g_yyyzzz_xxxx = cbuffer.data(ig_off + 360 * ccomps * dcomps);

            auto g_yyyzzz_xxxy = cbuffer.data(ig_off + 361 * ccomps * dcomps);

            auto g_yyyzzz_xxxz = cbuffer.data(ig_off + 362 * ccomps * dcomps);

            auto g_yyyzzz_xxyy = cbuffer.data(ig_off + 363 * ccomps * dcomps);

            auto g_yyyzzz_xxyz = cbuffer.data(ig_off + 364 * ccomps * dcomps);

            auto g_yyyzzz_xxzz = cbuffer.data(ig_off + 365 * ccomps * dcomps);

            auto g_yyyzzz_xyyy = cbuffer.data(ig_off + 366 * ccomps * dcomps);

            auto g_yyyzzz_xyyz = cbuffer.data(ig_off + 367 * ccomps * dcomps);

            auto g_yyyzzz_xyzz = cbuffer.data(ig_off + 368 * ccomps * dcomps);

            auto g_yyyzzz_xzzz = cbuffer.data(ig_off + 369 * ccomps * dcomps);

            auto g_yyyzzz_yyyy = cbuffer.data(ig_off + 370 * ccomps * dcomps);

            auto g_yyyzzz_yyyz = cbuffer.data(ig_off + 371 * ccomps * dcomps);

            auto g_yyyzzz_yyzz = cbuffer.data(ig_off + 372 * ccomps * dcomps);

            auto g_yyyzzz_yzzz = cbuffer.data(ig_off + 373 * ccomps * dcomps);

            auto g_yyyzzz_zzzz = cbuffer.data(ig_off + 374 * ccomps * dcomps);

            auto g_yyzzzz_xxxx = cbuffer.data(ig_off + 375 * ccomps * dcomps);

            auto g_yyzzzz_xxxy = cbuffer.data(ig_off + 376 * ccomps * dcomps);

            auto g_yyzzzz_xxxz = cbuffer.data(ig_off + 377 * ccomps * dcomps);

            auto g_yyzzzz_xxyy = cbuffer.data(ig_off + 378 * ccomps * dcomps);

            auto g_yyzzzz_xxyz = cbuffer.data(ig_off + 379 * ccomps * dcomps);

            auto g_yyzzzz_xxzz = cbuffer.data(ig_off + 380 * ccomps * dcomps);

            auto g_yyzzzz_xyyy = cbuffer.data(ig_off + 381 * ccomps * dcomps);

            auto g_yyzzzz_xyyz = cbuffer.data(ig_off + 382 * ccomps * dcomps);

            auto g_yyzzzz_xyzz = cbuffer.data(ig_off + 383 * ccomps * dcomps);

            auto g_yyzzzz_xzzz = cbuffer.data(ig_off + 384 * ccomps * dcomps);

            auto g_yyzzzz_yyyy = cbuffer.data(ig_off + 385 * ccomps * dcomps);

            auto g_yyzzzz_yyyz = cbuffer.data(ig_off + 386 * ccomps * dcomps);

            auto g_yyzzzz_yyzz = cbuffer.data(ig_off + 387 * ccomps * dcomps);

            auto g_yyzzzz_yzzz = cbuffer.data(ig_off + 388 * ccomps * dcomps);

            auto g_yyzzzz_zzzz = cbuffer.data(ig_off + 389 * ccomps * dcomps);

            auto g_yzzzzz_xxxx = cbuffer.data(ig_off + 390 * ccomps * dcomps);

            auto g_yzzzzz_xxxy = cbuffer.data(ig_off + 391 * ccomps * dcomps);

            auto g_yzzzzz_xxxz = cbuffer.data(ig_off + 392 * ccomps * dcomps);

            auto g_yzzzzz_xxyy = cbuffer.data(ig_off + 393 * ccomps * dcomps);

            auto g_yzzzzz_xxyz = cbuffer.data(ig_off + 394 * ccomps * dcomps);

            auto g_yzzzzz_xxzz = cbuffer.data(ig_off + 395 * ccomps * dcomps);

            auto g_yzzzzz_xyyy = cbuffer.data(ig_off + 396 * ccomps * dcomps);

            auto g_yzzzzz_xyyz = cbuffer.data(ig_off + 397 * ccomps * dcomps);

            auto g_yzzzzz_xyzz = cbuffer.data(ig_off + 398 * ccomps * dcomps);

            auto g_yzzzzz_xzzz = cbuffer.data(ig_off + 399 * ccomps * dcomps);

            auto g_yzzzzz_yyyy = cbuffer.data(ig_off + 400 * ccomps * dcomps);

            auto g_yzzzzz_yyyz = cbuffer.data(ig_off + 401 * ccomps * dcomps);

            auto g_yzzzzz_yyzz = cbuffer.data(ig_off + 402 * ccomps * dcomps);

            auto g_yzzzzz_yzzz = cbuffer.data(ig_off + 403 * ccomps * dcomps);

            auto g_yzzzzz_zzzz = cbuffer.data(ig_off + 404 * ccomps * dcomps);

            auto g_zzzzzz_xxxx = cbuffer.data(ig_off + 405 * ccomps * dcomps);

            auto g_zzzzzz_xxxy = cbuffer.data(ig_off + 406 * ccomps * dcomps);

            auto g_zzzzzz_xxxz = cbuffer.data(ig_off + 407 * ccomps * dcomps);

            auto g_zzzzzz_xxyy = cbuffer.data(ig_off + 408 * ccomps * dcomps);

            auto g_zzzzzz_xxyz = cbuffer.data(ig_off + 409 * ccomps * dcomps);

            auto g_zzzzzz_xxzz = cbuffer.data(ig_off + 410 * ccomps * dcomps);

            auto g_zzzzzz_xyyy = cbuffer.data(ig_off + 411 * ccomps * dcomps);

            auto g_zzzzzz_xyyz = cbuffer.data(ig_off + 412 * ccomps * dcomps);

            auto g_zzzzzz_xyzz = cbuffer.data(ig_off + 413 * ccomps * dcomps);

            auto g_zzzzzz_xzzz = cbuffer.data(ig_off + 414 * ccomps * dcomps);

            auto g_zzzzzz_yyyy = cbuffer.data(ig_off + 415 * ccomps * dcomps);

            auto g_zzzzzz_yyyz = cbuffer.data(ig_off + 416 * ccomps * dcomps);

            auto g_zzzzzz_yyzz = cbuffer.data(ig_off + 417 * ccomps * dcomps);

            auto g_zzzzzz_yzzz = cbuffer.data(ig_off + 418 * ccomps * dcomps);

            auto g_zzzzzz_zzzz = cbuffer.data(ig_off + 419 * ccomps * dcomps);

            /// Set up components of auxilary buffer : IGSS

            const auto ig_geom_01_off = idx_geom_01_igxx + i * dcomps + j;

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

            /// Set up components of auxilary buffer : IHSS

            const auto ih_geom_01_off = idx_geom_01_ihxx + i * dcomps + j;

            auto g_0_x_xxxxxx_xxxxx = cbuffer.data(ih_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxxxy = cbuffer.data(ih_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxxxz = cbuffer.data(ih_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxxyy = cbuffer.data(ih_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxxyz = cbuffer.data(ih_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxxzz = cbuffer.data(ih_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxyyy = cbuffer.data(ih_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxyyz = cbuffer.data(ih_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxyzz = cbuffer.data(ih_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xxzzz = cbuffer.data(ih_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyyyy = cbuffer.data(ih_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyyyz = cbuffer.data(ih_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyyzz = cbuffer.data(ih_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xyzzz = cbuffer.data(ih_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxx_xzzzz = cbuffer.data(ih_geom_01_off + 14 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyyyy = cbuffer.data(ih_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyyyz = cbuffer.data(ih_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyyzz = cbuffer.data(ih_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yyzzz = cbuffer.data(ih_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxxx_yzzzz = cbuffer.data(ih_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxxx_zzzzz = cbuffer.data(ih_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxxxx = cbuffer.data(ih_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxxxy = cbuffer.data(ih_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxxxz = cbuffer.data(ih_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxxyy = cbuffer.data(ih_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxxyz = cbuffer.data(ih_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxxzz = cbuffer.data(ih_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxyyy = cbuffer.data(ih_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxyyz = cbuffer.data(ih_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxyzz = cbuffer.data(ih_geom_01_off + 29 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xxzzz = cbuffer.data(ih_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyyyy = cbuffer.data(ih_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyyyz = cbuffer.data(ih_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyyzz = cbuffer.data(ih_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xyzzz = cbuffer.data(ih_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxxy_xzzzz = cbuffer.data(ih_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyyyy = cbuffer.data(ih_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyyyz = cbuffer.data(ih_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyyzz = cbuffer.data(ih_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yyzzz = cbuffer.data(ih_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxxxy_yzzzz = cbuffer.data(ih_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxxxy_zzzzz = cbuffer.data(ih_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxxxx = cbuffer.data(ih_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxxxy = cbuffer.data(ih_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxxxz = cbuffer.data(ih_geom_01_off + 44 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxxyy = cbuffer.data(ih_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxxyz = cbuffer.data(ih_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxxzz = cbuffer.data(ih_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxyyy = cbuffer.data(ih_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxyyz = cbuffer.data(ih_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxyzz = cbuffer.data(ih_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xxzzz = cbuffer.data(ih_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyyyy = cbuffer.data(ih_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyyyz = cbuffer.data(ih_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyyzz = cbuffer.data(ih_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xyzzz = cbuffer.data(ih_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxxxz_xzzzz = cbuffer.data(ih_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyyyy = cbuffer.data(ih_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyyyz = cbuffer.data(ih_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyyzz = cbuffer.data(ih_geom_01_off + 59 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yyzzz = cbuffer.data(ih_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxxxz_yzzzz = cbuffer.data(ih_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxxxz_zzzzz = cbuffer.data(ih_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxxxx = cbuffer.data(ih_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxxxy = cbuffer.data(ih_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxxxz = cbuffer.data(ih_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxxyy = cbuffer.data(ih_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxxyz = cbuffer.data(ih_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxxzz = cbuffer.data(ih_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxyyy = cbuffer.data(ih_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxyyz = cbuffer.data(ih_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxyzz = cbuffer.data(ih_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xxzzz = cbuffer.data(ih_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyyyy = cbuffer.data(ih_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyyyz = cbuffer.data(ih_geom_01_off + 74 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyyzz = cbuffer.data(ih_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xyzzz = cbuffer.data(ih_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxxyy_xzzzz = cbuffer.data(ih_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyyyy = cbuffer.data(ih_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyyyz = cbuffer.data(ih_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyyzz = cbuffer.data(ih_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yyzzz = cbuffer.data(ih_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxxyy_yzzzz = cbuffer.data(ih_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxxyy_zzzzz = cbuffer.data(ih_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxxxx = cbuffer.data(ih_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxxxy = cbuffer.data(ih_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxxxz = cbuffer.data(ih_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxxyy = cbuffer.data(ih_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxxyz = cbuffer.data(ih_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxxzz = cbuffer.data(ih_geom_01_off + 89 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxyyy = cbuffer.data(ih_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxyyz = cbuffer.data(ih_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxyzz = cbuffer.data(ih_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xxzzz = cbuffer.data(ih_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyyyy = cbuffer.data(ih_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyyyz = cbuffer.data(ih_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyyzz = cbuffer.data(ih_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xyzzz = cbuffer.data(ih_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxxxyz_xzzzz = cbuffer.data(ih_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyyyy = cbuffer.data(ih_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyyyz = cbuffer.data(ih_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyyzz = cbuffer.data(ih_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yyzzz = cbuffer.data(ih_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxxxyz_yzzzz = cbuffer.data(ih_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxxxyz_zzzzz = cbuffer.data(ih_geom_01_off + 104 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxxxx = cbuffer.data(ih_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxxxy = cbuffer.data(ih_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxxxz = cbuffer.data(ih_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxxyy = cbuffer.data(ih_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxxyz = cbuffer.data(ih_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxxzz = cbuffer.data(ih_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxyyy = cbuffer.data(ih_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxyyz = cbuffer.data(ih_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxyzz = cbuffer.data(ih_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xxzzz = cbuffer.data(ih_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyyyy = cbuffer.data(ih_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyyyz = cbuffer.data(ih_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyyzz = cbuffer.data(ih_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xyzzz = cbuffer.data(ih_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxxxzz_xzzzz = cbuffer.data(ih_geom_01_off + 119 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyyyy = cbuffer.data(ih_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyyyz = cbuffer.data(ih_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyyzz = cbuffer.data(ih_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yyzzz = cbuffer.data(ih_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxxxzz_yzzzz = cbuffer.data(ih_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxxxzz_zzzzz = cbuffer.data(ih_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxxxx = cbuffer.data(ih_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxxxy = cbuffer.data(ih_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxxxz = cbuffer.data(ih_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxxyy = cbuffer.data(ih_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxxyz = cbuffer.data(ih_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxxzz = cbuffer.data(ih_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxyyy = cbuffer.data(ih_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxyyz = cbuffer.data(ih_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxyzz = cbuffer.data(ih_geom_01_off + 134 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xxzzz = cbuffer.data(ih_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyyyy = cbuffer.data(ih_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyyyz = cbuffer.data(ih_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyyzz = cbuffer.data(ih_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xyzzz = cbuffer.data(ih_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxxyyy_xzzzz = cbuffer.data(ih_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyyyy = cbuffer.data(ih_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyyyz = cbuffer.data(ih_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyyzz = cbuffer.data(ih_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yyzzz = cbuffer.data(ih_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxxyyy_yzzzz = cbuffer.data(ih_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxxyyy_zzzzz = cbuffer.data(ih_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxxxx = cbuffer.data(ih_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxxxy = cbuffer.data(ih_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxxxz = cbuffer.data(ih_geom_01_off + 149 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxxyy = cbuffer.data(ih_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxxyz = cbuffer.data(ih_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxxzz = cbuffer.data(ih_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxyyy = cbuffer.data(ih_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxyyz = cbuffer.data(ih_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxyzz = cbuffer.data(ih_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xxzzz = cbuffer.data(ih_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyyyy = cbuffer.data(ih_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyyyz = cbuffer.data(ih_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyyzz = cbuffer.data(ih_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xyzzz = cbuffer.data(ih_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xxxyyz_xzzzz = cbuffer.data(ih_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyyyy = cbuffer.data(ih_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyyyz = cbuffer.data(ih_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyyzz = cbuffer.data(ih_geom_01_off + 164 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yyzzz = cbuffer.data(ih_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xxxyyz_yzzzz = cbuffer.data(ih_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xxxyyz_zzzzz = cbuffer.data(ih_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxxxx = cbuffer.data(ih_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxxxy = cbuffer.data(ih_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxxxz = cbuffer.data(ih_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxxyy = cbuffer.data(ih_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxxyz = cbuffer.data(ih_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxxzz = cbuffer.data(ih_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxyyy = cbuffer.data(ih_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxyyz = cbuffer.data(ih_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxyzz = cbuffer.data(ih_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xxzzz = cbuffer.data(ih_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyyyy = cbuffer.data(ih_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyyyz = cbuffer.data(ih_geom_01_off + 179 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyyzz = cbuffer.data(ih_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xyzzz = cbuffer.data(ih_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xxxyzz_xzzzz = cbuffer.data(ih_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyyyy = cbuffer.data(ih_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyyyz = cbuffer.data(ih_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyyzz = cbuffer.data(ih_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yyzzz = cbuffer.data(ih_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xxxyzz_yzzzz = cbuffer.data(ih_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xxxyzz_zzzzz = cbuffer.data(ih_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxxxx = cbuffer.data(ih_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxxxy = cbuffer.data(ih_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxxxz = cbuffer.data(ih_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxxyy = cbuffer.data(ih_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxxyz = cbuffer.data(ih_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxxzz = cbuffer.data(ih_geom_01_off + 194 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxyyy = cbuffer.data(ih_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxyyz = cbuffer.data(ih_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxyzz = cbuffer.data(ih_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xxzzz = cbuffer.data(ih_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyyyy = cbuffer.data(ih_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyyyz = cbuffer.data(ih_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyyzz = cbuffer.data(ih_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xyzzz = cbuffer.data(ih_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xxxzzz_xzzzz = cbuffer.data(ih_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyyyy = cbuffer.data(ih_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyyyz = cbuffer.data(ih_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyyzz = cbuffer.data(ih_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yyzzz = cbuffer.data(ih_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xxxzzz_yzzzz = cbuffer.data(ih_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xxxzzz_zzzzz = cbuffer.data(ih_geom_01_off + 209 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxxxx = cbuffer.data(ih_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxxxy = cbuffer.data(ih_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxxxz = cbuffer.data(ih_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxxyy = cbuffer.data(ih_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxxyz = cbuffer.data(ih_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxxzz = cbuffer.data(ih_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxyyy = cbuffer.data(ih_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxyyz = cbuffer.data(ih_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxyzz = cbuffer.data(ih_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xxzzz = cbuffer.data(ih_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyyyy = cbuffer.data(ih_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyyyz = cbuffer.data(ih_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyyzz = cbuffer.data(ih_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xyzzz = cbuffer.data(ih_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xxyyyy_xzzzz = cbuffer.data(ih_geom_01_off + 224 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyyyy = cbuffer.data(ih_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyyyz = cbuffer.data(ih_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyyzz = cbuffer.data(ih_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yyzzz = cbuffer.data(ih_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_xxyyyy_yzzzz = cbuffer.data(ih_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_xxyyyy_zzzzz = cbuffer.data(ih_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxxxx = cbuffer.data(ih_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxxxy = cbuffer.data(ih_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxxxz = cbuffer.data(ih_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxxyy = cbuffer.data(ih_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxxyz = cbuffer.data(ih_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxxzz = cbuffer.data(ih_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxyyy = cbuffer.data(ih_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxyyz = cbuffer.data(ih_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxyzz = cbuffer.data(ih_geom_01_off + 239 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xxzzz = cbuffer.data(ih_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyyyy = cbuffer.data(ih_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyyyz = cbuffer.data(ih_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyyzz = cbuffer.data(ih_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xyzzz = cbuffer.data(ih_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_xxyyyz_xzzzz = cbuffer.data(ih_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyyyy = cbuffer.data(ih_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyyyz = cbuffer.data(ih_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyyzz = cbuffer.data(ih_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yyzzz = cbuffer.data(ih_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_xxyyyz_yzzzz = cbuffer.data(ih_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_xxyyyz_zzzzz = cbuffer.data(ih_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxxxx = cbuffer.data(ih_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxxxy = cbuffer.data(ih_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxxxz = cbuffer.data(ih_geom_01_off + 254 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxxyy = cbuffer.data(ih_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxxyz = cbuffer.data(ih_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxxzz = cbuffer.data(ih_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxyyy = cbuffer.data(ih_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxyyz = cbuffer.data(ih_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxyzz = cbuffer.data(ih_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xxzzz = cbuffer.data(ih_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyyyy = cbuffer.data(ih_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyyyz = cbuffer.data(ih_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyyzz = cbuffer.data(ih_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xyzzz = cbuffer.data(ih_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_xxyyzz_xzzzz = cbuffer.data(ih_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyyyy = cbuffer.data(ih_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyyyz = cbuffer.data(ih_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyyzz = cbuffer.data(ih_geom_01_off + 269 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yyzzz = cbuffer.data(ih_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_xxyyzz_yzzzz = cbuffer.data(ih_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_xxyyzz_zzzzz = cbuffer.data(ih_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxxxx = cbuffer.data(ih_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxxxy = cbuffer.data(ih_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxxxz = cbuffer.data(ih_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxxyy = cbuffer.data(ih_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxxyz = cbuffer.data(ih_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxxzz = cbuffer.data(ih_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxyyy = cbuffer.data(ih_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxyyz = cbuffer.data(ih_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxyzz = cbuffer.data(ih_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xxzzz = cbuffer.data(ih_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyyyy = cbuffer.data(ih_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyyyz = cbuffer.data(ih_geom_01_off + 284 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyyzz = cbuffer.data(ih_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xyzzz = cbuffer.data(ih_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_xxyzzz_xzzzz = cbuffer.data(ih_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyyyy = cbuffer.data(ih_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyyyz = cbuffer.data(ih_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyyzz = cbuffer.data(ih_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yyzzz = cbuffer.data(ih_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_xxyzzz_yzzzz = cbuffer.data(ih_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_xxyzzz_zzzzz = cbuffer.data(ih_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 299 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_xxzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_xxzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_xxzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 314 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxxxx = cbuffer.data(ih_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxxxy = cbuffer.data(ih_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxxxz = cbuffer.data(ih_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxxyy = cbuffer.data(ih_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxxyz = cbuffer.data(ih_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxxzz = cbuffer.data(ih_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxyyy = cbuffer.data(ih_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxyyz = cbuffer.data(ih_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxyzz = cbuffer.data(ih_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xxzzz = cbuffer.data(ih_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyyyy = cbuffer.data(ih_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyyyz = cbuffer.data(ih_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyyzz = cbuffer.data(ih_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xyzzz = cbuffer.data(ih_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_xyyyyy_xzzzz = cbuffer.data(ih_geom_01_off + 329 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyyyy = cbuffer.data(ih_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyyyz = cbuffer.data(ih_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyyzz = cbuffer.data(ih_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yyzzz = cbuffer.data(ih_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_xyyyyy_yzzzz = cbuffer.data(ih_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_xyyyyy_zzzzz = cbuffer.data(ih_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxxxx = cbuffer.data(ih_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxxxy = cbuffer.data(ih_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxxxz = cbuffer.data(ih_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxxyy = cbuffer.data(ih_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxxyz = cbuffer.data(ih_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxxzz = cbuffer.data(ih_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxyyy = cbuffer.data(ih_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxyyz = cbuffer.data(ih_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxyzz = cbuffer.data(ih_geom_01_off + 344 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xxzzz = cbuffer.data(ih_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyyyy = cbuffer.data(ih_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyyyz = cbuffer.data(ih_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyyzz = cbuffer.data(ih_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xyzzz = cbuffer.data(ih_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_x_xyyyyz_xzzzz = cbuffer.data(ih_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyyyy = cbuffer.data(ih_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyyyz = cbuffer.data(ih_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyyzz = cbuffer.data(ih_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yyzzz = cbuffer.data(ih_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_xyyyyz_yzzzz = cbuffer.data(ih_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_xyyyyz_zzzzz = cbuffer.data(ih_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxxxx = cbuffer.data(ih_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxxxy = cbuffer.data(ih_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxxxz = cbuffer.data(ih_geom_01_off + 359 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxxyy = cbuffer.data(ih_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxxyz = cbuffer.data(ih_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxxzz = cbuffer.data(ih_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxyyy = cbuffer.data(ih_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxyyz = cbuffer.data(ih_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxyzz = cbuffer.data(ih_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xxzzz = cbuffer.data(ih_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyyyy = cbuffer.data(ih_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyyyz = cbuffer.data(ih_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyyzz = cbuffer.data(ih_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xyzzz = cbuffer.data(ih_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_x_xyyyzz_xzzzz = cbuffer.data(ih_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyyyy = cbuffer.data(ih_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyyyz = cbuffer.data(ih_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyyzz = cbuffer.data(ih_geom_01_off + 374 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yyzzz = cbuffer.data(ih_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_x_xyyyzz_yzzzz = cbuffer.data(ih_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_x_xyyyzz_zzzzz = cbuffer.data(ih_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxxxx = cbuffer.data(ih_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxxxy = cbuffer.data(ih_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxxxz = cbuffer.data(ih_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxxyy = cbuffer.data(ih_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxxyz = cbuffer.data(ih_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxxzz = cbuffer.data(ih_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxyyy = cbuffer.data(ih_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxyyz = cbuffer.data(ih_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxyzz = cbuffer.data(ih_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xxzzz = cbuffer.data(ih_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyyyy = cbuffer.data(ih_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyyyz = cbuffer.data(ih_geom_01_off + 389 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyyzz = cbuffer.data(ih_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xyzzz = cbuffer.data(ih_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_x_xyyzzz_xzzzz = cbuffer.data(ih_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyyyy = cbuffer.data(ih_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyyyz = cbuffer.data(ih_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyyzz = cbuffer.data(ih_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yyzzz = cbuffer.data(ih_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_x_xyyzzz_yzzzz = cbuffer.data(ih_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_x_xyyzzz_zzzzz = cbuffer.data(ih_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 404 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_x_xyzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_x_xyzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_x_xyzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 419 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_x_xzzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 434 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_x_xzzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_x_xzzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxxxx = cbuffer.data(ih_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxxxy = cbuffer.data(ih_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxxxz = cbuffer.data(ih_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxxyy = cbuffer.data(ih_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxxyz = cbuffer.data(ih_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxxzz = cbuffer.data(ih_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxyyy = cbuffer.data(ih_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxyyz = cbuffer.data(ih_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxyzz = cbuffer.data(ih_geom_01_off + 449 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xxzzz = cbuffer.data(ih_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyyyy = cbuffer.data(ih_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyyyz = cbuffer.data(ih_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyyzz = cbuffer.data(ih_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xyzzz = cbuffer.data(ih_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_x_yyyyyy_xzzzz = cbuffer.data(ih_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyyyy = cbuffer.data(ih_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyyyz = cbuffer.data(ih_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyyzz = cbuffer.data(ih_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yyzzz = cbuffer.data(ih_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_x_yyyyyy_yzzzz = cbuffer.data(ih_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_x_yyyyyy_zzzzz = cbuffer.data(ih_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxxxx = cbuffer.data(ih_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxxxy = cbuffer.data(ih_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxxxz = cbuffer.data(ih_geom_01_off + 464 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxxyy = cbuffer.data(ih_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxxyz = cbuffer.data(ih_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxxzz = cbuffer.data(ih_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxyyy = cbuffer.data(ih_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxyyz = cbuffer.data(ih_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxyzz = cbuffer.data(ih_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xxzzz = cbuffer.data(ih_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyyyy = cbuffer.data(ih_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyyyz = cbuffer.data(ih_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyyzz = cbuffer.data(ih_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xyzzz = cbuffer.data(ih_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_x_yyyyyz_xzzzz = cbuffer.data(ih_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyyyy = cbuffer.data(ih_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyyyz = cbuffer.data(ih_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyyzz = cbuffer.data(ih_geom_01_off + 479 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yyzzz = cbuffer.data(ih_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_x_yyyyyz_yzzzz = cbuffer.data(ih_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_x_yyyyyz_zzzzz = cbuffer.data(ih_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxxxx = cbuffer.data(ih_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxxxy = cbuffer.data(ih_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxxxz = cbuffer.data(ih_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxxyy = cbuffer.data(ih_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxxyz = cbuffer.data(ih_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxxzz = cbuffer.data(ih_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxyyy = cbuffer.data(ih_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxyyz = cbuffer.data(ih_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxyzz = cbuffer.data(ih_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xxzzz = cbuffer.data(ih_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyyyy = cbuffer.data(ih_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyyyz = cbuffer.data(ih_geom_01_off + 494 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyyzz = cbuffer.data(ih_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xyzzz = cbuffer.data(ih_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_x_yyyyzz_xzzzz = cbuffer.data(ih_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyyyy = cbuffer.data(ih_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyyyz = cbuffer.data(ih_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyyzz = cbuffer.data(ih_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yyzzz = cbuffer.data(ih_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_x_yyyyzz_yzzzz = cbuffer.data(ih_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_x_yyyyzz_zzzzz = cbuffer.data(ih_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxxxx = cbuffer.data(ih_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxxxy = cbuffer.data(ih_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxxxz = cbuffer.data(ih_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxxyy = cbuffer.data(ih_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxxyz = cbuffer.data(ih_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxxzz = cbuffer.data(ih_geom_01_off + 509 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxyyy = cbuffer.data(ih_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxyyz = cbuffer.data(ih_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxyzz = cbuffer.data(ih_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xxzzz = cbuffer.data(ih_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyyyy = cbuffer.data(ih_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyyyz = cbuffer.data(ih_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyyzz = cbuffer.data(ih_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xyzzz = cbuffer.data(ih_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_x_yyyzzz_xzzzz = cbuffer.data(ih_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyyyy = cbuffer.data(ih_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyyyz = cbuffer.data(ih_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyyzz = cbuffer.data(ih_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yyzzz = cbuffer.data(ih_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_x_yyyzzz_yzzzz = cbuffer.data(ih_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_x_yyyzzz_zzzzz = cbuffer.data(ih_geom_01_off + 524 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_x_yyzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 539 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_x_yyzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_x_yyzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 554 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_x_yzzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_x_yzzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_x_yzzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 569 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_x_zzzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 584 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_x_zzzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_x_zzzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxxxx = cbuffer.data(ih_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxxxy = cbuffer.data(ih_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxxxz = cbuffer.data(ih_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxxyy = cbuffer.data(ih_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxxyz = cbuffer.data(ih_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxxzz = cbuffer.data(ih_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxyyy = cbuffer.data(ih_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxyyz = cbuffer.data(ih_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxyzz = cbuffer.data(ih_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xxzzz = cbuffer.data(ih_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyyyy = cbuffer.data(ih_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyyyz = cbuffer.data(ih_geom_01_off + 599 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyyzz = cbuffer.data(ih_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xyzzz = cbuffer.data(ih_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_xxxxxx_xzzzz = cbuffer.data(ih_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyyyy = cbuffer.data(ih_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyyyz = cbuffer.data(ih_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyyzz = cbuffer.data(ih_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yyzzz = cbuffer.data(ih_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_xxxxxx_yzzzz = cbuffer.data(ih_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_xxxxxx_zzzzz = cbuffer.data(ih_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxxxx = cbuffer.data(ih_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxxxy = cbuffer.data(ih_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxxxz = cbuffer.data(ih_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxxyy = cbuffer.data(ih_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxxyz = cbuffer.data(ih_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxxzz = cbuffer.data(ih_geom_01_off + 614 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxyyy = cbuffer.data(ih_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxyyz = cbuffer.data(ih_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxyzz = cbuffer.data(ih_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xxzzz = cbuffer.data(ih_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyyyy = cbuffer.data(ih_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyyyz = cbuffer.data(ih_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyyzz = cbuffer.data(ih_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xyzzz = cbuffer.data(ih_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_xxxxxy_xzzzz = cbuffer.data(ih_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyyyy = cbuffer.data(ih_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyyyz = cbuffer.data(ih_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyyzz = cbuffer.data(ih_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yyzzz = cbuffer.data(ih_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_xxxxxy_yzzzz = cbuffer.data(ih_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_xxxxxy_zzzzz = cbuffer.data(ih_geom_01_off + 629 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxxxx = cbuffer.data(ih_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxxxy = cbuffer.data(ih_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxxxz = cbuffer.data(ih_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxxyy = cbuffer.data(ih_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxxyz = cbuffer.data(ih_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxxzz = cbuffer.data(ih_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxyyy = cbuffer.data(ih_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxyyz = cbuffer.data(ih_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxyzz = cbuffer.data(ih_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xxzzz = cbuffer.data(ih_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyyyy = cbuffer.data(ih_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyyyz = cbuffer.data(ih_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyyzz = cbuffer.data(ih_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xyzzz = cbuffer.data(ih_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_y_xxxxxz_xzzzz = cbuffer.data(ih_geom_01_off + 644 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyyyy = cbuffer.data(ih_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyyyz = cbuffer.data(ih_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyyzz = cbuffer.data(ih_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yyzzz = cbuffer.data(ih_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_y_xxxxxz_yzzzz = cbuffer.data(ih_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_y_xxxxxz_zzzzz = cbuffer.data(ih_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxxxx = cbuffer.data(ih_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxxxy = cbuffer.data(ih_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxxxz = cbuffer.data(ih_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxxyy = cbuffer.data(ih_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxxyz = cbuffer.data(ih_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxxzz = cbuffer.data(ih_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxyyy = cbuffer.data(ih_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxyyz = cbuffer.data(ih_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxyzz = cbuffer.data(ih_geom_01_off + 659 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xxzzz = cbuffer.data(ih_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyyyy = cbuffer.data(ih_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyyyz = cbuffer.data(ih_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyyzz = cbuffer.data(ih_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xyzzz = cbuffer.data(ih_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_y_xxxxyy_xzzzz = cbuffer.data(ih_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyyyy = cbuffer.data(ih_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyyyz = cbuffer.data(ih_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyyzz = cbuffer.data(ih_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yyzzz = cbuffer.data(ih_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_y_xxxxyy_yzzzz = cbuffer.data(ih_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_y_xxxxyy_zzzzz = cbuffer.data(ih_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxxxx = cbuffer.data(ih_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxxxy = cbuffer.data(ih_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxxxz = cbuffer.data(ih_geom_01_off + 674 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxxyy = cbuffer.data(ih_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxxyz = cbuffer.data(ih_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxxzz = cbuffer.data(ih_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxyyy = cbuffer.data(ih_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxyyz = cbuffer.data(ih_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxyzz = cbuffer.data(ih_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xxzzz = cbuffer.data(ih_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyyyy = cbuffer.data(ih_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyyyz = cbuffer.data(ih_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyyzz = cbuffer.data(ih_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xyzzz = cbuffer.data(ih_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_y_xxxxyz_xzzzz = cbuffer.data(ih_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyyyy = cbuffer.data(ih_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyyyz = cbuffer.data(ih_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyyzz = cbuffer.data(ih_geom_01_off + 689 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yyzzz = cbuffer.data(ih_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_y_xxxxyz_yzzzz = cbuffer.data(ih_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_y_xxxxyz_zzzzz = cbuffer.data(ih_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxxxx = cbuffer.data(ih_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxxxy = cbuffer.data(ih_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxxxz = cbuffer.data(ih_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxxyy = cbuffer.data(ih_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxxyz = cbuffer.data(ih_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxxzz = cbuffer.data(ih_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxyyy = cbuffer.data(ih_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxyyz = cbuffer.data(ih_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxyzz = cbuffer.data(ih_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xxzzz = cbuffer.data(ih_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyyyy = cbuffer.data(ih_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyyyz = cbuffer.data(ih_geom_01_off + 704 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyyzz = cbuffer.data(ih_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xyzzz = cbuffer.data(ih_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_y_xxxxzz_xzzzz = cbuffer.data(ih_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyyyy = cbuffer.data(ih_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyyyz = cbuffer.data(ih_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyyzz = cbuffer.data(ih_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yyzzz = cbuffer.data(ih_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_y_xxxxzz_yzzzz = cbuffer.data(ih_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_y_xxxxzz_zzzzz = cbuffer.data(ih_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxxxx = cbuffer.data(ih_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxxxy = cbuffer.data(ih_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxxxz = cbuffer.data(ih_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxxyy = cbuffer.data(ih_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxxyz = cbuffer.data(ih_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxxzz = cbuffer.data(ih_geom_01_off + 719 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxyyy = cbuffer.data(ih_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxyyz = cbuffer.data(ih_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxyzz = cbuffer.data(ih_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xxzzz = cbuffer.data(ih_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyyyy = cbuffer.data(ih_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyyyz = cbuffer.data(ih_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyyzz = cbuffer.data(ih_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xyzzz = cbuffer.data(ih_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_y_xxxyyy_xzzzz = cbuffer.data(ih_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyyyy = cbuffer.data(ih_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyyyz = cbuffer.data(ih_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyyzz = cbuffer.data(ih_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yyzzz = cbuffer.data(ih_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_y_xxxyyy_yzzzz = cbuffer.data(ih_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_y_xxxyyy_zzzzz = cbuffer.data(ih_geom_01_off + 734 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxxxx = cbuffer.data(ih_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxxxy = cbuffer.data(ih_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxxxz = cbuffer.data(ih_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxxyy = cbuffer.data(ih_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxxyz = cbuffer.data(ih_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxxzz = cbuffer.data(ih_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxyyy = cbuffer.data(ih_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxyyz = cbuffer.data(ih_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxyzz = cbuffer.data(ih_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xxzzz = cbuffer.data(ih_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyyyy = cbuffer.data(ih_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyyyz = cbuffer.data(ih_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyyzz = cbuffer.data(ih_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xyzzz = cbuffer.data(ih_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_y_xxxyyz_xzzzz = cbuffer.data(ih_geom_01_off + 749 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyyyy = cbuffer.data(ih_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyyyz = cbuffer.data(ih_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyyzz = cbuffer.data(ih_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yyzzz = cbuffer.data(ih_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_y_xxxyyz_yzzzz = cbuffer.data(ih_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_y_xxxyyz_zzzzz = cbuffer.data(ih_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxxxx = cbuffer.data(ih_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxxxy = cbuffer.data(ih_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxxxz = cbuffer.data(ih_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxxyy = cbuffer.data(ih_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxxyz = cbuffer.data(ih_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxxzz = cbuffer.data(ih_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxyyy = cbuffer.data(ih_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxyyz = cbuffer.data(ih_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxyzz = cbuffer.data(ih_geom_01_off + 764 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xxzzz = cbuffer.data(ih_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyyyy = cbuffer.data(ih_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyyyz = cbuffer.data(ih_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyyzz = cbuffer.data(ih_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xyzzz = cbuffer.data(ih_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_y_xxxyzz_xzzzz = cbuffer.data(ih_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyyyy = cbuffer.data(ih_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyyyz = cbuffer.data(ih_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyyzz = cbuffer.data(ih_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yyzzz = cbuffer.data(ih_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_y_xxxyzz_yzzzz = cbuffer.data(ih_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_y_xxxyzz_zzzzz = cbuffer.data(ih_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxxxx = cbuffer.data(ih_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxxxy = cbuffer.data(ih_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxxxz = cbuffer.data(ih_geom_01_off + 779 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxxyy = cbuffer.data(ih_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxxyz = cbuffer.data(ih_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxxzz = cbuffer.data(ih_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxyyy = cbuffer.data(ih_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxyyz = cbuffer.data(ih_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxyzz = cbuffer.data(ih_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xxzzz = cbuffer.data(ih_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyyyy = cbuffer.data(ih_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyyyz = cbuffer.data(ih_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyyzz = cbuffer.data(ih_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xyzzz = cbuffer.data(ih_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_y_xxxzzz_xzzzz = cbuffer.data(ih_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyyyy = cbuffer.data(ih_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyyyz = cbuffer.data(ih_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyyzz = cbuffer.data(ih_geom_01_off + 794 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yyzzz = cbuffer.data(ih_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_y_xxxzzz_yzzzz = cbuffer.data(ih_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_y_xxxzzz_zzzzz = cbuffer.data(ih_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxxxx = cbuffer.data(ih_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxxxy = cbuffer.data(ih_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxxxz = cbuffer.data(ih_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxxyy = cbuffer.data(ih_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxxyz = cbuffer.data(ih_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxxzz = cbuffer.data(ih_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxyyy = cbuffer.data(ih_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxyyz = cbuffer.data(ih_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxyzz = cbuffer.data(ih_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xxzzz = cbuffer.data(ih_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyyyy = cbuffer.data(ih_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyyyz = cbuffer.data(ih_geom_01_off + 809 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyyzz = cbuffer.data(ih_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xyzzz = cbuffer.data(ih_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_y_xxyyyy_xzzzz = cbuffer.data(ih_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyyyy = cbuffer.data(ih_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyyyz = cbuffer.data(ih_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyyzz = cbuffer.data(ih_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yyzzz = cbuffer.data(ih_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_y_xxyyyy_yzzzz = cbuffer.data(ih_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_y_xxyyyy_zzzzz = cbuffer.data(ih_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxxxx = cbuffer.data(ih_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxxxy = cbuffer.data(ih_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxxxz = cbuffer.data(ih_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxxyy = cbuffer.data(ih_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxxyz = cbuffer.data(ih_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxxzz = cbuffer.data(ih_geom_01_off + 824 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxyyy = cbuffer.data(ih_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxyyz = cbuffer.data(ih_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxyzz = cbuffer.data(ih_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xxzzz = cbuffer.data(ih_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyyyy = cbuffer.data(ih_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyyyz = cbuffer.data(ih_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyyzz = cbuffer.data(ih_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xyzzz = cbuffer.data(ih_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_y_xxyyyz_xzzzz = cbuffer.data(ih_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyyyy = cbuffer.data(ih_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyyyz = cbuffer.data(ih_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyyzz = cbuffer.data(ih_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yyzzz = cbuffer.data(ih_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_y_xxyyyz_yzzzz = cbuffer.data(ih_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_y_xxyyyz_zzzzz = cbuffer.data(ih_geom_01_off + 839 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxxxx = cbuffer.data(ih_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxxxy = cbuffer.data(ih_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxxxz = cbuffer.data(ih_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxxyy = cbuffer.data(ih_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxxyz = cbuffer.data(ih_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxxzz = cbuffer.data(ih_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxyyy = cbuffer.data(ih_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxyyz = cbuffer.data(ih_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxyzz = cbuffer.data(ih_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xxzzz = cbuffer.data(ih_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyyyy = cbuffer.data(ih_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyyyz = cbuffer.data(ih_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyyzz = cbuffer.data(ih_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xyzzz = cbuffer.data(ih_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_y_xxyyzz_xzzzz = cbuffer.data(ih_geom_01_off + 854 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyyyy = cbuffer.data(ih_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyyyz = cbuffer.data(ih_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyyzz = cbuffer.data(ih_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yyzzz = cbuffer.data(ih_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_y_xxyyzz_yzzzz = cbuffer.data(ih_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_y_xxyyzz_zzzzz = cbuffer.data(ih_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxxxx = cbuffer.data(ih_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxxxy = cbuffer.data(ih_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxxxz = cbuffer.data(ih_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxxyy = cbuffer.data(ih_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxxyz = cbuffer.data(ih_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxxzz = cbuffer.data(ih_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxyyy = cbuffer.data(ih_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxyyz = cbuffer.data(ih_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxyzz = cbuffer.data(ih_geom_01_off + 869 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xxzzz = cbuffer.data(ih_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyyyy = cbuffer.data(ih_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyyyz = cbuffer.data(ih_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyyzz = cbuffer.data(ih_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xyzzz = cbuffer.data(ih_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_y_xxyzzz_xzzzz = cbuffer.data(ih_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyyyy = cbuffer.data(ih_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyyyz = cbuffer.data(ih_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyyzz = cbuffer.data(ih_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yyzzz = cbuffer.data(ih_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_y_xxyzzz_yzzzz = cbuffer.data(ih_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_y_xxyzzz_zzzzz = cbuffer.data(ih_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 884 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_y_xxzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 899 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_y_xxzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_y_xxzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxxxx = cbuffer.data(ih_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxxxy = cbuffer.data(ih_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxxxz = cbuffer.data(ih_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxxyy = cbuffer.data(ih_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxxyz = cbuffer.data(ih_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxxzz = cbuffer.data(ih_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxyyy = cbuffer.data(ih_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxyyz = cbuffer.data(ih_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxyzz = cbuffer.data(ih_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xxzzz = cbuffer.data(ih_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyyyy = cbuffer.data(ih_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyyyz = cbuffer.data(ih_geom_01_off + 914 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyyzz = cbuffer.data(ih_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xyzzz = cbuffer.data(ih_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_y_xyyyyy_xzzzz = cbuffer.data(ih_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyyyy = cbuffer.data(ih_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyyyz = cbuffer.data(ih_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyyzz = cbuffer.data(ih_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yyzzz = cbuffer.data(ih_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_y_xyyyyy_yzzzz = cbuffer.data(ih_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_y_xyyyyy_zzzzz = cbuffer.data(ih_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxxxx = cbuffer.data(ih_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxxxy = cbuffer.data(ih_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxxxz = cbuffer.data(ih_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxxyy = cbuffer.data(ih_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxxyz = cbuffer.data(ih_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxxzz = cbuffer.data(ih_geom_01_off + 929 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxyyy = cbuffer.data(ih_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxyyz = cbuffer.data(ih_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxyzz = cbuffer.data(ih_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xxzzz = cbuffer.data(ih_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyyyy = cbuffer.data(ih_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyyyz = cbuffer.data(ih_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyyzz = cbuffer.data(ih_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xyzzz = cbuffer.data(ih_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_y_xyyyyz_xzzzz = cbuffer.data(ih_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyyyy = cbuffer.data(ih_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyyyz = cbuffer.data(ih_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyyzz = cbuffer.data(ih_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yyzzz = cbuffer.data(ih_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_y_xyyyyz_yzzzz = cbuffer.data(ih_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_y_xyyyyz_zzzzz = cbuffer.data(ih_geom_01_off + 944 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxxxx = cbuffer.data(ih_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxxxy = cbuffer.data(ih_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxxxz = cbuffer.data(ih_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxxyy = cbuffer.data(ih_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxxyz = cbuffer.data(ih_geom_01_off + 949 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxxzz = cbuffer.data(ih_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxyyy = cbuffer.data(ih_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxyyz = cbuffer.data(ih_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxyzz = cbuffer.data(ih_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xxzzz = cbuffer.data(ih_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyyyy = cbuffer.data(ih_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyyyz = cbuffer.data(ih_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyyzz = cbuffer.data(ih_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xyzzz = cbuffer.data(ih_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_y_xyyyzz_xzzzz = cbuffer.data(ih_geom_01_off + 959 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyyyy = cbuffer.data(ih_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyyyz = cbuffer.data(ih_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyyzz = cbuffer.data(ih_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yyzzz = cbuffer.data(ih_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_y_xyyyzz_yzzzz = cbuffer.data(ih_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_y_xyyyzz_zzzzz = cbuffer.data(ih_geom_01_off + 965 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxxxx = cbuffer.data(ih_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxxxy = cbuffer.data(ih_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxxxz = cbuffer.data(ih_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxxyy = cbuffer.data(ih_geom_01_off + 969 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxxyz = cbuffer.data(ih_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxxzz = cbuffer.data(ih_geom_01_off + 971 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxyyy = cbuffer.data(ih_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxyyz = cbuffer.data(ih_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxyzz = cbuffer.data(ih_geom_01_off + 974 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xxzzz = cbuffer.data(ih_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyyyy = cbuffer.data(ih_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyyyz = cbuffer.data(ih_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyyzz = cbuffer.data(ih_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xyzzz = cbuffer.data(ih_geom_01_off + 979 * ccomps * dcomps);

            auto g_0_y_xyyzzz_xzzzz = cbuffer.data(ih_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyyyy = cbuffer.data(ih_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyyyz = cbuffer.data(ih_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyyzz = cbuffer.data(ih_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yyzzz = cbuffer.data(ih_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_y_xyyzzz_yzzzz = cbuffer.data(ih_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_y_xyyzzz_zzzzz = cbuffer.data(ih_geom_01_off + 986 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 989 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 999 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_y_xyzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1004 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_y_xyzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_y_xyzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1007 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1009 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1019 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_y_xzzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_y_xzzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_y_xzzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1028 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxxxx = cbuffer.data(ih_geom_01_off + 1029 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxxxy = cbuffer.data(ih_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxxxz = cbuffer.data(ih_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxxyy = cbuffer.data(ih_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxxyz = cbuffer.data(ih_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxxzz = cbuffer.data(ih_geom_01_off + 1034 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxyyy = cbuffer.data(ih_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxyyz = cbuffer.data(ih_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxyzz = cbuffer.data(ih_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xxzzz = cbuffer.data(ih_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyyyy = cbuffer.data(ih_geom_01_off + 1039 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyyyz = cbuffer.data(ih_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyyzz = cbuffer.data(ih_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xyzzz = cbuffer.data(ih_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_y_yyyyyy_xzzzz = cbuffer.data(ih_geom_01_off + 1043 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyyyy = cbuffer.data(ih_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyyyz = cbuffer.data(ih_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyyzz = cbuffer.data(ih_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yyzzz = cbuffer.data(ih_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_y_yyyyyy_yzzzz = cbuffer.data(ih_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_y_yyyyyy_zzzzz = cbuffer.data(ih_geom_01_off + 1049 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxxxx = cbuffer.data(ih_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxxxy = cbuffer.data(ih_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxxxz = cbuffer.data(ih_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxxyy = cbuffer.data(ih_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxxyz = cbuffer.data(ih_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxxzz = cbuffer.data(ih_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxyyy = cbuffer.data(ih_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxyyz = cbuffer.data(ih_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxyzz = cbuffer.data(ih_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xxzzz = cbuffer.data(ih_geom_01_off + 1059 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyyyy = cbuffer.data(ih_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyyyz = cbuffer.data(ih_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyyzz = cbuffer.data(ih_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xyzzz = cbuffer.data(ih_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_y_yyyyyz_xzzzz = cbuffer.data(ih_geom_01_off + 1064 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyyyy = cbuffer.data(ih_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyyyz = cbuffer.data(ih_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyyzz = cbuffer.data(ih_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yyzzz = cbuffer.data(ih_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_y_yyyyyz_yzzzz = cbuffer.data(ih_geom_01_off + 1069 * ccomps * dcomps);

            auto g_0_y_yyyyyz_zzzzz = cbuffer.data(ih_geom_01_off + 1070 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxxxx = cbuffer.data(ih_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxxxy = cbuffer.data(ih_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxxxz = cbuffer.data(ih_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxxyy = cbuffer.data(ih_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxxyz = cbuffer.data(ih_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxxzz = cbuffer.data(ih_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxyyy = cbuffer.data(ih_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxyyz = cbuffer.data(ih_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxyzz = cbuffer.data(ih_geom_01_off + 1079 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xxzzz = cbuffer.data(ih_geom_01_off + 1080 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyyyy = cbuffer.data(ih_geom_01_off + 1081 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyyyz = cbuffer.data(ih_geom_01_off + 1082 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyyzz = cbuffer.data(ih_geom_01_off + 1083 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xyzzz = cbuffer.data(ih_geom_01_off + 1084 * ccomps * dcomps);

            auto g_0_y_yyyyzz_xzzzz = cbuffer.data(ih_geom_01_off + 1085 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyyyy = cbuffer.data(ih_geom_01_off + 1086 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyyyz = cbuffer.data(ih_geom_01_off + 1087 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyyzz = cbuffer.data(ih_geom_01_off + 1088 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yyzzz = cbuffer.data(ih_geom_01_off + 1089 * ccomps * dcomps);

            auto g_0_y_yyyyzz_yzzzz = cbuffer.data(ih_geom_01_off + 1090 * ccomps * dcomps);

            auto g_0_y_yyyyzz_zzzzz = cbuffer.data(ih_geom_01_off + 1091 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1092 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1093 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1094 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1095 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1096 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1097 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1098 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1099 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1100 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1101 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1102 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1103 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1104 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1105 * ccomps * dcomps);

            auto g_0_y_yyyzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1106 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1107 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1108 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1109 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1110 * ccomps * dcomps);

            auto g_0_y_yyyzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1111 * ccomps * dcomps);

            auto g_0_y_yyyzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1112 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1113 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1114 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1115 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1116 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1117 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1118 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1119 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1120 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1121 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1122 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1123 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1124 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1125 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1126 * ccomps * dcomps);

            auto g_0_y_yyzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1127 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1128 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1129 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1130 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1131 * ccomps * dcomps);

            auto g_0_y_yyzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1132 * ccomps * dcomps);

            auto g_0_y_yyzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1133 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1134 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1135 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1136 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1137 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1138 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1139 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1140 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1141 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1142 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1143 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1144 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1145 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1146 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1147 * ccomps * dcomps);

            auto g_0_y_yzzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1148 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1149 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1150 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1151 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1152 * ccomps * dcomps);

            auto g_0_y_yzzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1153 * ccomps * dcomps);

            auto g_0_y_yzzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1154 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1155 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1156 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1157 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1158 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1159 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1160 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1161 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1162 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1163 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1164 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1165 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1166 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1167 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1168 * ccomps * dcomps);

            auto g_0_y_zzzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1169 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1170 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1171 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1172 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1173 * ccomps * dcomps);

            auto g_0_y_zzzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1174 * ccomps * dcomps);

            auto g_0_y_zzzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1175 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxxxx = cbuffer.data(ih_geom_01_off + 1176 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxxxy = cbuffer.data(ih_geom_01_off + 1177 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxxxz = cbuffer.data(ih_geom_01_off + 1178 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxxyy = cbuffer.data(ih_geom_01_off + 1179 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxxyz = cbuffer.data(ih_geom_01_off + 1180 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxxzz = cbuffer.data(ih_geom_01_off + 1181 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxyyy = cbuffer.data(ih_geom_01_off + 1182 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxyyz = cbuffer.data(ih_geom_01_off + 1183 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxyzz = cbuffer.data(ih_geom_01_off + 1184 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xxzzz = cbuffer.data(ih_geom_01_off + 1185 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyyyy = cbuffer.data(ih_geom_01_off + 1186 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyyyz = cbuffer.data(ih_geom_01_off + 1187 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyyzz = cbuffer.data(ih_geom_01_off + 1188 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xyzzz = cbuffer.data(ih_geom_01_off + 1189 * ccomps * dcomps);

            auto g_0_z_xxxxxx_xzzzz = cbuffer.data(ih_geom_01_off + 1190 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyyyy = cbuffer.data(ih_geom_01_off + 1191 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyyyz = cbuffer.data(ih_geom_01_off + 1192 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyyzz = cbuffer.data(ih_geom_01_off + 1193 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yyzzz = cbuffer.data(ih_geom_01_off + 1194 * ccomps * dcomps);

            auto g_0_z_xxxxxx_yzzzz = cbuffer.data(ih_geom_01_off + 1195 * ccomps * dcomps);

            auto g_0_z_xxxxxx_zzzzz = cbuffer.data(ih_geom_01_off + 1196 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxxxx = cbuffer.data(ih_geom_01_off + 1197 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxxxy = cbuffer.data(ih_geom_01_off + 1198 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxxxz = cbuffer.data(ih_geom_01_off + 1199 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxxyy = cbuffer.data(ih_geom_01_off + 1200 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxxyz = cbuffer.data(ih_geom_01_off + 1201 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxxzz = cbuffer.data(ih_geom_01_off + 1202 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxyyy = cbuffer.data(ih_geom_01_off + 1203 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxyyz = cbuffer.data(ih_geom_01_off + 1204 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxyzz = cbuffer.data(ih_geom_01_off + 1205 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xxzzz = cbuffer.data(ih_geom_01_off + 1206 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyyyy = cbuffer.data(ih_geom_01_off + 1207 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyyyz = cbuffer.data(ih_geom_01_off + 1208 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyyzz = cbuffer.data(ih_geom_01_off + 1209 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xyzzz = cbuffer.data(ih_geom_01_off + 1210 * ccomps * dcomps);

            auto g_0_z_xxxxxy_xzzzz = cbuffer.data(ih_geom_01_off + 1211 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyyyy = cbuffer.data(ih_geom_01_off + 1212 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyyyz = cbuffer.data(ih_geom_01_off + 1213 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyyzz = cbuffer.data(ih_geom_01_off + 1214 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yyzzz = cbuffer.data(ih_geom_01_off + 1215 * ccomps * dcomps);

            auto g_0_z_xxxxxy_yzzzz = cbuffer.data(ih_geom_01_off + 1216 * ccomps * dcomps);

            auto g_0_z_xxxxxy_zzzzz = cbuffer.data(ih_geom_01_off + 1217 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxxxx = cbuffer.data(ih_geom_01_off + 1218 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxxxy = cbuffer.data(ih_geom_01_off + 1219 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxxxz = cbuffer.data(ih_geom_01_off + 1220 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxxyy = cbuffer.data(ih_geom_01_off + 1221 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxxyz = cbuffer.data(ih_geom_01_off + 1222 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxxzz = cbuffer.data(ih_geom_01_off + 1223 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxyyy = cbuffer.data(ih_geom_01_off + 1224 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxyyz = cbuffer.data(ih_geom_01_off + 1225 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxyzz = cbuffer.data(ih_geom_01_off + 1226 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xxzzz = cbuffer.data(ih_geom_01_off + 1227 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyyyy = cbuffer.data(ih_geom_01_off + 1228 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyyyz = cbuffer.data(ih_geom_01_off + 1229 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyyzz = cbuffer.data(ih_geom_01_off + 1230 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xyzzz = cbuffer.data(ih_geom_01_off + 1231 * ccomps * dcomps);

            auto g_0_z_xxxxxz_xzzzz = cbuffer.data(ih_geom_01_off + 1232 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyyyy = cbuffer.data(ih_geom_01_off + 1233 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyyyz = cbuffer.data(ih_geom_01_off + 1234 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyyzz = cbuffer.data(ih_geom_01_off + 1235 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yyzzz = cbuffer.data(ih_geom_01_off + 1236 * ccomps * dcomps);

            auto g_0_z_xxxxxz_yzzzz = cbuffer.data(ih_geom_01_off + 1237 * ccomps * dcomps);

            auto g_0_z_xxxxxz_zzzzz = cbuffer.data(ih_geom_01_off + 1238 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxxxx = cbuffer.data(ih_geom_01_off + 1239 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxxxy = cbuffer.data(ih_geom_01_off + 1240 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxxxz = cbuffer.data(ih_geom_01_off + 1241 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxxyy = cbuffer.data(ih_geom_01_off + 1242 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxxyz = cbuffer.data(ih_geom_01_off + 1243 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxxzz = cbuffer.data(ih_geom_01_off + 1244 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxyyy = cbuffer.data(ih_geom_01_off + 1245 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxyyz = cbuffer.data(ih_geom_01_off + 1246 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxyzz = cbuffer.data(ih_geom_01_off + 1247 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xxzzz = cbuffer.data(ih_geom_01_off + 1248 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyyyy = cbuffer.data(ih_geom_01_off + 1249 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyyyz = cbuffer.data(ih_geom_01_off + 1250 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyyzz = cbuffer.data(ih_geom_01_off + 1251 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xyzzz = cbuffer.data(ih_geom_01_off + 1252 * ccomps * dcomps);

            auto g_0_z_xxxxyy_xzzzz = cbuffer.data(ih_geom_01_off + 1253 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyyyy = cbuffer.data(ih_geom_01_off + 1254 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyyyz = cbuffer.data(ih_geom_01_off + 1255 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyyzz = cbuffer.data(ih_geom_01_off + 1256 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yyzzz = cbuffer.data(ih_geom_01_off + 1257 * ccomps * dcomps);

            auto g_0_z_xxxxyy_yzzzz = cbuffer.data(ih_geom_01_off + 1258 * ccomps * dcomps);

            auto g_0_z_xxxxyy_zzzzz = cbuffer.data(ih_geom_01_off + 1259 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxxxx = cbuffer.data(ih_geom_01_off + 1260 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxxxy = cbuffer.data(ih_geom_01_off + 1261 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxxxz = cbuffer.data(ih_geom_01_off + 1262 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxxyy = cbuffer.data(ih_geom_01_off + 1263 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxxyz = cbuffer.data(ih_geom_01_off + 1264 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxxzz = cbuffer.data(ih_geom_01_off + 1265 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxyyy = cbuffer.data(ih_geom_01_off + 1266 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxyyz = cbuffer.data(ih_geom_01_off + 1267 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxyzz = cbuffer.data(ih_geom_01_off + 1268 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xxzzz = cbuffer.data(ih_geom_01_off + 1269 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyyyy = cbuffer.data(ih_geom_01_off + 1270 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyyyz = cbuffer.data(ih_geom_01_off + 1271 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyyzz = cbuffer.data(ih_geom_01_off + 1272 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xyzzz = cbuffer.data(ih_geom_01_off + 1273 * ccomps * dcomps);

            auto g_0_z_xxxxyz_xzzzz = cbuffer.data(ih_geom_01_off + 1274 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyyyy = cbuffer.data(ih_geom_01_off + 1275 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyyyz = cbuffer.data(ih_geom_01_off + 1276 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyyzz = cbuffer.data(ih_geom_01_off + 1277 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yyzzz = cbuffer.data(ih_geom_01_off + 1278 * ccomps * dcomps);

            auto g_0_z_xxxxyz_yzzzz = cbuffer.data(ih_geom_01_off + 1279 * ccomps * dcomps);

            auto g_0_z_xxxxyz_zzzzz = cbuffer.data(ih_geom_01_off + 1280 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxxxx = cbuffer.data(ih_geom_01_off + 1281 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxxxy = cbuffer.data(ih_geom_01_off + 1282 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxxxz = cbuffer.data(ih_geom_01_off + 1283 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxxyy = cbuffer.data(ih_geom_01_off + 1284 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxxyz = cbuffer.data(ih_geom_01_off + 1285 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxxzz = cbuffer.data(ih_geom_01_off + 1286 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxyyy = cbuffer.data(ih_geom_01_off + 1287 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxyyz = cbuffer.data(ih_geom_01_off + 1288 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxyzz = cbuffer.data(ih_geom_01_off + 1289 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xxzzz = cbuffer.data(ih_geom_01_off + 1290 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyyyy = cbuffer.data(ih_geom_01_off + 1291 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyyyz = cbuffer.data(ih_geom_01_off + 1292 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyyzz = cbuffer.data(ih_geom_01_off + 1293 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xyzzz = cbuffer.data(ih_geom_01_off + 1294 * ccomps * dcomps);

            auto g_0_z_xxxxzz_xzzzz = cbuffer.data(ih_geom_01_off + 1295 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyyyy = cbuffer.data(ih_geom_01_off + 1296 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyyyz = cbuffer.data(ih_geom_01_off + 1297 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyyzz = cbuffer.data(ih_geom_01_off + 1298 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yyzzz = cbuffer.data(ih_geom_01_off + 1299 * ccomps * dcomps);

            auto g_0_z_xxxxzz_yzzzz = cbuffer.data(ih_geom_01_off + 1300 * ccomps * dcomps);

            auto g_0_z_xxxxzz_zzzzz = cbuffer.data(ih_geom_01_off + 1301 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxxxx = cbuffer.data(ih_geom_01_off + 1302 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxxxy = cbuffer.data(ih_geom_01_off + 1303 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxxxz = cbuffer.data(ih_geom_01_off + 1304 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxxyy = cbuffer.data(ih_geom_01_off + 1305 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxxyz = cbuffer.data(ih_geom_01_off + 1306 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxxzz = cbuffer.data(ih_geom_01_off + 1307 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxyyy = cbuffer.data(ih_geom_01_off + 1308 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxyyz = cbuffer.data(ih_geom_01_off + 1309 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxyzz = cbuffer.data(ih_geom_01_off + 1310 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xxzzz = cbuffer.data(ih_geom_01_off + 1311 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyyyy = cbuffer.data(ih_geom_01_off + 1312 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyyyz = cbuffer.data(ih_geom_01_off + 1313 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyyzz = cbuffer.data(ih_geom_01_off + 1314 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xyzzz = cbuffer.data(ih_geom_01_off + 1315 * ccomps * dcomps);

            auto g_0_z_xxxyyy_xzzzz = cbuffer.data(ih_geom_01_off + 1316 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyyyy = cbuffer.data(ih_geom_01_off + 1317 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyyyz = cbuffer.data(ih_geom_01_off + 1318 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyyzz = cbuffer.data(ih_geom_01_off + 1319 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yyzzz = cbuffer.data(ih_geom_01_off + 1320 * ccomps * dcomps);

            auto g_0_z_xxxyyy_yzzzz = cbuffer.data(ih_geom_01_off + 1321 * ccomps * dcomps);

            auto g_0_z_xxxyyy_zzzzz = cbuffer.data(ih_geom_01_off + 1322 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxxxx = cbuffer.data(ih_geom_01_off + 1323 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxxxy = cbuffer.data(ih_geom_01_off + 1324 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxxxz = cbuffer.data(ih_geom_01_off + 1325 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxxyy = cbuffer.data(ih_geom_01_off + 1326 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxxyz = cbuffer.data(ih_geom_01_off + 1327 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxxzz = cbuffer.data(ih_geom_01_off + 1328 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxyyy = cbuffer.data(ih_geom_01_off + 1329 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxyyz = cbuffer.data(ih_geom_01_off + 1330 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxyzz = cbuffer.data(ih_geom_01_off + 1331 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xxzzz = cbuffer.data(ih_geom_01_off + 1332 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyyyy = cbuffer.data(ih_geom_01_off + 1333 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyyyz = cbuffer.data(ih_geom_01_off + 1334 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyyzz = cbuffer.data(ih_geom_01_off + 1335 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xyzzz = cbuffer.data(ih_geom_01_off + 1336 * ccomps * dcomps);

            auto g_0_z_xxxyyz_xzzzz = cbuffer.data(ih_geom_01_off + 1337 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyyyy = cbuffer.data(ih_geom_01_off + 1338 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyyyz = cbuffer.data(ih_geom_01_off + 1339 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyyzz = cbuffer.data(ih_geom_01_off + 1340 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yyzzz = cbuffer.data(ih_geom_01_off + 1341 * ccomps * dcomps);

            auto g_0_z_xxxyyz_yzzzz = cbuffer.data(ih_geom_01_off + 1342 * ccomps * dcomps);

            auto g_0_z_xxxyyz_zzzzz = cbuffer.data(ih_geom_01_off + 1343 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxxxx = cbuffer.data(ih_geom_01_off + 1344 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxxxy = cbuffer.data(ih_geom_01_off + 1345 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxxxz = cbuffer.data(ih_geom_01_off + 1346 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxxyy = cbuffer.data(ih_geom_01_off + 1347 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxxyz = cbuffer.data(ih_geom_01_off + 1348 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxxzz = cbuffer.data(ih_geom_01_off + 1349 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxyyy = cbuffer.data(ih_geom_01_off + 1350 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxyyz = cbuffer.data(ih_geom_01_off + 1351 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxyzz = cbuffer.data(ih_geom_01_off + 1352 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xxzzz = cbuffer.data(ih_geom_01_off + 1353 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyyyy = cbuffer.data(ih_geom_01_off + 1354 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyyyz = cbuffer.data(ih_geom_01_off + 1355 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyyzz = cbuffer.data(ih_geom_01_off + 1356 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xyzzz = cbuffer.data(ih_geom_01_off + 1357 * ccomps * dcomps);

            auto g_0_z_xxxyzz_xzzzz = cbuffer.data(ih_geom_01_off + 1358 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyyyy = cbuffer.data(ih_geom_01_off + 1359 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyyyz = cbuffer.data(ih_geom_01_off + 1360 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyyzz = cbuffer.data(ih_geom_01_off + 1361 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yyzzz = cbuffer.data(ih_geom_01_off + 1362 * ccomps * dcomps);

            auto g_0_z_xxxyzz_yzzzz = cbuffer.data(ih_geom_01_off + 1363 * ccomps * dcomps);

            auto g_0_z_xxxyzz_zzzzz = cbuffer.data(ih_geom_01_off + 1364 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1365 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1366 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1367 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1368 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1369 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1370 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1371 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1372 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1373 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1374 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1375 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1376 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1377 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1378 * ccomps * dcomps);

            auto g_0_z_xxxzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1379 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1380 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1381 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1382 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1383 * ccomps * dcomps);

            auto g_0_z_xxxzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1384 * ccomps * dcomps);

            auto g_0_z_xxxzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1385 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxxxx = cbuffer.data(ih_geom_01_off + 1386 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxxxy = cbuffer.data(ih_geom_01_off + 1387 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxxxz = cbuffer.data(ih_geom_01_off + 1388 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxxyy = cbuffer.data(ih_geom_01_off + 1389 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxxyz = cbuffer.data(ih_geom_01_off + 1390 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxxzz = cbuffer.data(ih_geom_01_off + 1391 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxyyy = cbuffer.data(ih_geom_01_off + 1392 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxyyz = cbuffer.data(ih_geom_01_off + 1393 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxyzz = cbuffer.data(ih_geom_01_off + 1394 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xxzzz = cbuffer.data(ih_geom_01_off + 1395 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyyyy = cbuffer.data(ih_geom_01_off + 1396 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyyyz = cbuffer.data(ih_geom_01_off + 1397 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyyzz = cbuffer.data(ih_geom_01_off + 1398 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xyzzz = cbuffer.data(ih_geom_01_off + 1399 * ccomps * dcomps);

            auto g_0_z_xxyyyy_xzzzz = cbuffer.data(ih_geom_01_off + 1400 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyyyy = cbuffer.data(ih_geom_01_off + 1401 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyyyz = cbuffer.data(ih_geom_01_off + 1402 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyyzz = cbuffer.data(ih_geom_01_off + 1403 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yyzzz = cbuffer.data(ih_geom_01_off + 1404 * ccomps * dcomps);

            auto g_0_z_xxyyyy_yzzzz = cbuffer.data(ih_geom_01_off + 1405 * ccomps * dcomps);

            auto g_0_z_xxyyyy_zzzzz = cbuffer.data(ih_geom_01_off + 1406 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxxxx = cbuffer.data(ih_geom_01_off + 1407 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxxxy = cbuffer.data(ih_geom_01_off + 1408 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxxxz = cbuffer.data(ih_geom_01_off + 1409 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxxyy = cbuffer.data(ih_geom_01_off + 1410 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxxyz = cbuffer.data(ih_geom_01_off + 1411 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxxzz = cbuffer.data(ih_geom_01_off + 1412 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxyyy = cbuffer.data(ih_geom_01_off + 1413 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxyyz = cbuffer.data(ih_geom_01_off + 1414 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxyzz = cbuffer.data(ih_geom_01_off + 1415 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xxzzz = cbuffer.data(ih_geom_01_off + 1416 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyyyy = cbuffer.data(ih_geom_01_off + 1417 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyyyz = cbuffer.data(ih_geom_01_off + 1418 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyyzz = cbuffer.data(ih_geom_01_off + 1419 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xyzzz = cbuffer.data(ih_geom_01_off + 1420 * ccomps * dcomps);

            auto g_0_z_xxyyyz_xzzzz = cbuffer.data(ih_geom_01_off + 1421 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyyyy = cbuffer.data(ih_geom_01_off + 1422 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyyyz = cbuffer.data(ih_geom_01_off + 1423 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyyzz = cbuffer.data(ih_geom_01_off + 1424 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yyzzz = cbuffer.data(ih_geom_01_off + 1425 * ccomps * dcomps);

            auto g_0_z_xxyyyz_yzzzz = cbuffer.data(ih_geom_01_off + 1426 * ccomps * dcomps);

            auto g_0_z_xxyyyz_zzzzz = cbuffer.data(ih_geom_01_off + 1427 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxxxx = cbuffer.data(ih_geom_01_off + 1428 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxxxy = cbuffer.data(ih_geom_01_off + 1429 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxxxz = cbuffer.data(ih_geom_01_off + 1430 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxxyy = cbuffer.data(ih_geom_01_off + 1431 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxxyz = cbuffer.data(ih_geom_01_off + 1432 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxxzz = cbuffer.data(ih_geom_01_off + 1433 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxyyy = cbuffer.data(ih_geom_01_off + 1434 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxyyz = cbuffer.data(ih_geom_01_off + 1435 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxyzz = cbuffer.data(ih_geom_01_off + 1436 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xxzzz = cbuffer.data(ih_geom_01_off + 1437 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyyyy = cbuffer.data(ih_geom_01_off + 1438 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyyyz = cbuffer.data(ih_geom_01_off + 1439 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyyzz = cbuffer.data(ih_geom_01_off + 1440 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xyzzz = cbuffer.data(ih_geom_01_off + 1441 * ccomps * dcomps);

            auto g_0_z_xxyyzz_xzzzz = cbuffer.data(ih_geom_01_off + 1442 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyyyy = cbuffer.data(ih_geom_01_off + 1443 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyyyz = cbuffer.data(ih_geom_01_off + 1444 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyyzz = cbuffer.data(ih_geom_01_off + 1445 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yyzzz = cbuffer.data(ih_geom_01_off + 1446 * ccomps * dcomps);

            auto g_0_z_xxyyzz_yzzzz = cbuffer.data(ih_geom_01_off + 1447 * ccomps * dcomps);

            auto g_0_z_xxyyzz_zzzzz = cbuffer.data(ih_geom_01_off + 1448 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1449 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1450 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1451 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1452 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1453 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1454 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1455 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1456 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1457 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1458 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1459 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1460 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1461 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1462 * ccomps * dcomps);

            auto g_0_z_xxyzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1463 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1464 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1465 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1466 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1467 * ccomps * dcomps);

            auto g_0_z_xxyzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1468 * ccomps * dcomps);

            auto g_0_z_xxyzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1469 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1470 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1471 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1472 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1473 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1474 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1475 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1476 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1477 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1478 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1479 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1480 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1481 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1482 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1483 * ccomps * dcomps);

            auto g_0_z_xxzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1484 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1485 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1486 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1487 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1488 * ccomps * dcomps);

            auto g_0_z_xxzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1489 * ccomps * dcomps);

            auto g_0_z_xxzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1490 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxxxx = cbuffer.data(ih_geom_01_off + 1491 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxxxy = cbuffer.data(ih_geom_01_off + 1492 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxxxz = cbuffer.data(ih_geom_01_off + 1493 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxxyy = cbuffer.data(ih_geom_01_off + 1494 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxxyz = cbuffer.data(ih_geom_01_off + 1495 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxxzz = cbuffer.data(ih_geom_01_off + 1496 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxyyy = cbuffer.data(ih_geom_01_off + 1497 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxyyz = cbuffer.data(ih_geom_01_off + 1498 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxyzz = cbuffer.data(ih_geom_01_off + 1499 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xxzzz = cbuffer.data(ih_geom_01_off + 1500 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyyyy = cbuffer.data(ih_geom_01_off + 1501 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyyyz = cbuffer.data(ih_geom_01_off + 1502 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyyzz = cbuffer.data(ih_geom_01_off + 1503 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xyzzz = cbuffer.data(ih_geom_01_off + 1504 * ccomps * dcomps);

            auto g_0_z_xyyyyy_xzzzz = cbuffer.data(ih_geom_01_off + 1505 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyyyy = cbuffer.data(ih_geom_01_off + 1506 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyyyz = cbuffer.data(ih_geom_01_off + 1507 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyyzz = cbuffer.data(ih_geom_01_off + 1508 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yyzzz = cbuffer.data(ih_geom_01_off + 1509 * ccomps * dcomps);

            auto g_0_z_xyyyyy_yzzzz = cbuffer.data(ih_geom_01_off + 1510 * ccomps * dcomps);

            auto g_0_z_xyyyyy_zzzzz = cbuffer.data(ih_geom_01_off + 1511 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxxxx = cbuffer.data(ih_geom_01_off + 1512 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxxxy = cbuffer.data(ih_geom_01_off + 1513 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxxxz = cbuffer.data(ih_geom_01_off + 1514 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxxyy = cbuffer.data(ih_geom_01_off + 1515 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxxyz = cbuffer.data(ih_geom_01_off + 1516 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxxzz = cbuffer.data(ih_geom_01_off + 1517 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxyyy = cbuffer.data(ih_geom_01_off + 1518 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxyyz = cbuffer.data(ih_geom_01_off + 1519 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxyzz = cbuffer.data(ih_geom_01_off + 1520 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xxzzz = cbuffer.data(ih_geom_01_off + 1521 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyyyy = cbuffer.data(ih_geom_01_off + 1522 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyyyz = cbuffer.data(ih_geom_01_off + 1523 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyyzz = cbuffer.data(ih_geom_01_off + 1524 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xyzzz = cbuffer.data(ih_geom_01_off + 1525 * ccomps * dcomps);

            auto g_0_z_xyyyyz_xzzzz = cbuffer.data(ih_geom_01_off + 1526 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyyyy = cbuffer.data(ih_geom_01_off + 1527 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyyyz = cbuffer.data(ih_geom_01_off + 1528 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyyzz = cbuffer.data(ih_geom_01_off + 1529 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yyzzz = cbuffer.data(ih_geom_01_off + 1530 * ccomps * dcomps);

            auto g_0_z_xyyyyz_yzzzz = cbuffer.data(ih_geom_01_off + 1531 * ccomps * dcomps);

            auto g_0_z_xyyyyz_zzzzz = cbuffer.data(ih_geom_01_off + 1532 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxxxx = cbuffer.data(ih_geom_01_off + 1533 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxxxy = cbuffer.data(ih_geom_01_off + 1534 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxxxz = cbuffer.data(ih_geom_01_off + 1535 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxxyy = cbuffer.data(ih_geom_01_off + 1536 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxxyz = cbuffer.data(ih_geom_01_off + 1537 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxxzz = cbuffer.data(ih_geom_01_off + 1538 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxyyy = cbuffer.data(ih_geom_01_off + 1539 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxyyz = cbuffer.data(ih_geom_01_off + 1540 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxyzz = cbuffer.data(ih_geom_01_off + 1541 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xxzzz = cbuffer.data(ih_geom_01_off + 1542 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyyyy = cbuffer.data(ih_geom_01_off + 1543 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyyyz = cbuffer.data(ih_geom_01_off + 1544 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyyzz = cbuffer.data(ih_geom_01_off + 1545 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xyzzz = cbuffer.data(ih_geom_01_off + 1546 * ccomps * dcomps);

            auto g_0_z_xyyyzz_xzzzz = cbuffer.data(ih_geom_01_off + 1547 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyyyy = cbuffer.data(ih_geom_01_off + 1548 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyyyz = cbuffer.data(ih_geom_01_off + 1549 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyyzz = cbuffer.data(ih_geom_01_off + 1550 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yyzzz = cbuffer.data(ih_geom_01_off + 1551 * ccomps * dcomps);

            auto g_0_z_xyyyzz_yzzzz = cbuffer.data(ih_geom_01_off + 1552 * ccomps * dcomps);

            auto g_0_z_xyyyzz_zzzzz = cbuffer.data(ih_geom_01_off + 1553 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1554 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1555 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1556 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1557 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1558 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1559 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1560 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1561 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1562 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1563 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1564 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1565 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1566 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1567 * ccomps * dcomps);

            auto g_0_z_xyyzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1568 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1569 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1570 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1571 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1572 * ccomps * dcomps);

            auto g_0_z_xyyzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1573 * ccomps * dcomps);

            auto g_0_z_xyyzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1574 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1575 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1576 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1577 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1578 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1579 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1580 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1581 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1582 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1583 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1584 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1585 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1586 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1587 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1588 * ccomps * dcomps);

            auto g_0_z_xyzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1589 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1590 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1591 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1592 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1593 * ccomps * dcomps);

            auto g_0_z_xyzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1594 * ccomps * dcomps);

            auto g_0_z_xyzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1595 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1596 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1597 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1598 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1599 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1600 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1601 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1602 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1603 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1604 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1605 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1606 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1607 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1608 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1609 * ccomps * dcomps);

            auto g_0_z_xzzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1610 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1611 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1612 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1613 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1614 * ccomps * dcomps);

            auto g_0_z_xzzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1615 * ccomps * dcomps);

            auto g_0_z_xzzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1616 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxxxx = cbuffer.data(ih_geom_01_off + 1617 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxxxy = cbuffer.data(ih_geom_01_off + 1618 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxxxz = cbuffer.data(ih_geom_01_off + 1619 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxxyy = cbuffer.data(ih_geom_01_off + 1620 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxxyz = cbuffer.data(ih_geom_01_off + 1621 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxxzz = cbuffer.data(ih_geom_01_off + 1622 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxyyy = cbuffer.data(ih_geom_01_off + 1623 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxyyz = cbuffer.data(ih_geom_01_off + 1624 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxyzz = cbuffer.data(ih_geom_01_off + 1625 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xxzzz = cbuffer.data(ih_geom_01_off + 1626 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyyyy = cbuffer.data(ih_geom_01_off + 1627 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyyyz = cbuffer.data(ih_geom_01_off + 1628 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyyzz = cbuffer.data(ih_geom_01_off + 1629 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xyzzz = cbuffer.data(ih_geom_01_off + 1630 * ccomps * dcomps);

            auto g_0_z_yyyyyy_xzzzz = cbuffer.data(ih_geom_01_off + 1631 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyyyy = cbuffer.data(ih_geom_01_off + 1632 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyyyz = cbuffer.data(ih_geom_01_off + 1633 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyyzz = cbuffer.data(ih_geom_01_off + 1634 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yyzzz = cbuffer.data(ih_geom_01_off + 1635 * ccomps * dcomps);

            auto g_0_z_yyyyyy_yzzzz = cbuffer.data(ih_geom_01_off + 1636 * ccomps * dcomps);

            auto g_0_z_yyyyyy_zzzzz = cbuffer.data(ih_geom_01_off + 1637 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxxxx = cbuffer.data(ih_geom_01_off + 1638 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxxxy = cbuffer.data(ih_geom_01_off + 1639 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxxxz = cbuffer.data(ih_geom_01_off + 1640 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxxyy = cbuffer.data(ih_geom_01_off + 1641 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxxyz = cbuffer.data(ih_geom_01_off + 1642 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxxzz = cbuffer.data(ih_geom_01_off + 1643 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxyyy = cbuffer.data(ih_geom_01_off + 1644 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxyyz = cbuffer.data(ih_geom_01_off + 1645 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxyzz = cbuffer.data(ih_geom_01_off + 1646 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xxzzz = cbuffer.data(ih_geom_01_off + 1647 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyyyy = cbuffer.data(ih_geom_01_off + 1648 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyyyz = cbuffer.data(ih_geom_01_off + 1649 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyyzz = cbuffer.data(ih_geom_01_off + 1650 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xyzzz = cbuffer.data(ih_geom_01_off + 1651 * ccomps * dcomps);

            auto g_0_z_yyyyyz_xzzzz = cbuffer.data(ih_geom_01_off + 1652 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyyyy = cbuffer.data(ih_geom_01_off + 1653 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyyyz = cbuffer.data(ih_geom_01_off + 1654 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyyzz = cbuffer.data(ih_geom_01_off + 1655 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yyzzz = cbuffer.data(ih_geom_01_off + 1656 * ccomps * dcomps);

            auto g_0_z_yyyyyz_yzzzz = cbuffer.data(ih_geom_01_off + 1657 * ccomps * dcomps);

            auto g_0_z_yyyyyz_zzzzz = cbuffer.data(ih_geom_01_off + 1658 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxxxx = cbuffer.data(ih_geom_01_off + 1659 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxxxy = cbuffer.data(ih_geom_01_off + 1660 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxxxz = cbuffer.data(ih_geom_01_off + 1661 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxxyy = cbuffer.data(ih_geom_01_off + 1662 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxxyz = cbuffer.data(ih_geom_01_off + 1663 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxxzz = cbuffer.data(ih_geom_01_off + 1664 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxyyy = cbuffer.data(ih_geom_01_off + 1665 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxyyz = cbuffer.data(ih_geom_01_off + 1666 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxyzz = cbuffer.data(ih_geom_01_off + 1667 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xxzzz = cbuffer.data(ih_geom_01_off + 1668 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyyyy = cbuffer.data(ih_geom_01_off + 1669 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyyyz = cbuffer.data(ih_geom_01_off + 1670 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyyzz = cbuffer.data(ih_geom_01_off + 1671 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xyzzz = cbuffer.data(ih_geom_01_off + 1672 * ccomps * dcomps);

            auto g_0_z_yyyyzz_xzzzz = cbuffer.data(ih_geom_01_off + 1673 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyyyy = cbuffer.data(ih_geom_01_off + 1674 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyyyz = cbuffer.data(ih_geom_01_off + 1675 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyyzz = cbuffer.data(ih_geom_01_off + 1676 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yyzzz = cbuffer.data(ih_geom_01_off + 1677 * ccomps * dcomps);

            auto g_0_z_yyyyzz_yzzzz = cbuffer.data(ih_geom_01_off + 1678 * ccomps * dcomps);

            auto g_0_z_yyyyzz_zzzzz = cbuffer.data(ih_geom_01_off + 1679 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1680 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1681 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1682 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1683 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1684 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1685 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1686 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1687 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1688 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1689 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1690 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1691 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1692 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1693 * ccomps * dcomps);

            auto g_0_z_yyyzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1694 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1695 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1696 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1697 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1698 * ccomps * dcomps);

            auto g_0_z_yyyzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1699 * ccomps * dcomps);

            auto g_0_z_yyyzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1700 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1701 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1702 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1703 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1704 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1705 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1706 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1707 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1708 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1709 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1710 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1711 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1712 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1713 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1714 * ccomps * dcomps);

            auto g_0_z_yyzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1715 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1716 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1717 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1718 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1719 * ccomps * dcomps);

            auto g_0_z_yyzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1720 * ccomps * dcomps);

            auto g_0_z_yyzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1721 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1722 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1723 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1724 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1725 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1726 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1727 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1728 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1729 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1730 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1731 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1732 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1733 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1734 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1735 * ccomps * dcomps);

            auto g_0_z_yzzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1736 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1737 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1738 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1739 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1740 * ccomps * dcomps);

            auto g_0_z_yzzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1741 * ccomps * dcomps);

            auto g_0_z_yzzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1742 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxxxx = cbuffer.data(ih_geom_01_off + 1743 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxxxy = cbuffer.data(ih_geom_01_off + 1744 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxxxz = cbuffer.data(ih_geom_01_off + 1745 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxxyy = cbuffer.data(ih_geom_01_off + 1746 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxxyz = cbuffer.data(ih_geom_01_off + 1747 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxxzz = cbuffer.data(ih_geom_01_off + 1748 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxyyy = cbuffer.data(ih_geom_01_off + 1749 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxyyz = cbuffer.data(ih_geom_01_off + 1750 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxyzz = cbuffer.data(ih_geom_01_off + 1751 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xxzzz = cbuffer.data(ih_geom_01_off + 1752 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyyyy = cbuffer.data(ih_geom_01_off + 1753 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyyyz = cbuffer.data(ih_geom_01_off + 1754 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyyzz = cbuffer.data(ih_geom_01_off + 1755 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xyzzz = cbuffer.data(ih_geom_01_off + 1756 * ccomps * dcomps);

            auto g_0_z_zzzzzz_xzzzz = cbuffer.data(ih_geom_01_off + 1757 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyyyy = cbuffer.data(ih_geom_01_off + 1758 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyyyz = cbuffer.data(ih_geom_01_off + 1759 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyyzz = cbuffer.data(ih_geom_01_off + 1760 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yyzzz = cbuffer.data(ih_geom_01_off + 1761 * ccomps * dcomps);

            auto g_0_z_zzzzzz_yzzzz = cbuffer.data(ih_geom_01_off + 1762 * ccomps * dcomps);

            auto g_0_z_zzzzzz_zzzzz = cbuffer.data(ih_geom_01_off + 1763 * ccomps * dcomps);

            /// set up bra offset for contr_buffer_kgxx

            const auto kg_geom_01_off = idx_geom_01_kgxx + i * dcomps + j;

            /// Set up 0-15 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxx_xxxx = cbuffer.data(kg_geom_01_off + 0 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xxxy = cbuffer.data(kg_geom_01_off + 1 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xxxz = cbuffer.data(kg_geom_01_off + 2 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xxyy = cbuffer.data(kg_geom_01_off + 3 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xxyz = cbuffer.data(kg_geom_01_off + 4 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xxzz = cbuffer.data(kg_geom_01_off + 5 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xyyy = cbuffer.data(kg_geom_01_off + 6 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xyyz = cbuffer.data(kg_geom_01_off + 7 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xyzz = cbuffer.data(kg_geom_01_off + 8 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_xzzz = cbuffer.data(kg_geom_01_off + 9 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_yyyy = cbuffer.data(kg_geom_01_off + 10 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_yyyz = cbuffer.data(kg_geom_01_off + 11 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_yyzz = cbuffer.data(kg_geom_01_off + 12 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_yzzz = cbuffer.data(kg_geom_01_off + 13 * ccomps * dcomps);

            auto g_0_x_xxxxxxx_zzzz = cbuffer.data(kg_geom_01_off + 14 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_xxxx, g_0_x_xxxxxx_xxxxx, g_0_x_xxxxxx_xxxxy, g_0_x_xxxxxx_xxxxz, g_0_x_xxxxxx_xxxy, g_0_x_xxxxxx_xxxyy, g_0_x_xxxxxx_xxxyz, g_0_x_xxxxxx_xxxz, g_0_x_xxxxxx_xxxzz, g_0_x_xxxxxx_xxyy, g_0_x_xxxxxx_xxyyy, g_0_x_xxxxxx_xxyyz, g_0_x_xxxxxx_xxyz, g_0_x_xxxxxx_xxyzz, g_0_x_xxxxxx_xxzz, g_0_x_xxxxxx_xxzzz, g_0_x_xxxxxx_xyyy, g_0_x_xxxxxx_xyyyy, g_0_x_xxxxxx_xyyyz, g_0_x_xxxxxx_xyyz, g_0_x_xxxxxx_xyyzz, g_0_x_xxxxxx_xyzz, g_0_x_xxxxxx_xyzzz, g_0_x_xxxxxx_xzzz, g_0_x_xxxxxx_xzzzz, g_0_x_xxxxxx_yyyy, g_0_x_xxxxxx_yyyz, g_0_x_xxxxxx_yyzz, g_0_x_xxxxxx_yzzz, g_0_x_xxxxxx_zzzz, g_0_x_xxxxxxx_xxxx, g_0_x_xxxxxxx_xxxy, g_0_x_xxxxxxx_xxxz, g_0_x_xxxxxxx_xxyy, g_0_x_xxxxxxx_xxyz, g_0_x_xxxxxxx_xxzz, g_0_x_xxxxxxx_xyyy, g_0_x_xxxxxxx_xyyz, g_0_x_xxxxxxx_xyzz, g_0_x_xxxxxxx_xzzz, g_0_x_xxxxxxx_yyyy, g_0_x_xxxxxxx_yyyz, g_0_x_xxxxxxx_yyzz, g_0_x_xxxxxxx_yzzz, g_0_x_xxxxxxx_zzzz, g_xxxxxx_xxxx, g_xxxxxx_xxxy, g_xxxxxx_xxxz, g_xxxxxx_xxyy, g_xxxxxx_xxyz, g_xxxxxx_xxzz, g_xxxxxx_xyyy, g_xxxxxx_xyyz, g_xxxxxx_xyzz, g_xxxxxx_xzzz, g_xxxxxx_yyyy, g_xxxxxx_yyyz, g_xxxxxx_yyzz, g_xxxxxx_yzzz, g_xxxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxx_xxxx[k] = g_xxxxxx_xxxx[k] - g_0_x_xxxxxx_xxxx[k] * ab_x + g_0_x_xxxxxx_xxxxx[k];

                g_0_x_xxxxxxx_xxxy[k] = g_xxxxxx_xxxy[k] - g_0_x_xxxxxx_xxxy[k] * ab_x + g_0_x_xxxxxx_xxxxy[k];

                g_0_x_xxxxxxx_xxxz[k] = g_xxxxxx_xxxz[k] - g_0_x_xxxxxx_xxxz[k] * ab_x + g_0_x_xxxxxx_xxxxz[k];

                g_0_x_xxxxxxx_xxyy[k] = g_xxxxxx_xxyy[k] - g_0_x_xxxxxx_xxyy[k] * ab_x + g_0_x_xxxxxx_xxxyy[k];

                g_0_x_xxxxxxx_xxyz[k] = g_xxxxxx_xxyz[k] - g_0_x_xxxxxx_xxyz[k] * ab_x + g_0_x_xxxxxx_xxxyz[k];

                g_0_x_xxxxxxx_xxzz[k] = g_xxxxxx_xxzz[k] - g_0_x_xxxxxx_xxzz[k] * ab_x + g_0_x_xxxxxx_xxxzz[k];

                g_0_x_xxxxxxx_xyyy[k] = g_xxxxxx_xyyy[k] - g_0_x_xxxxxx_xyyy[k] * ab_x + g_0_x_xxxxxx_xxyyy[k];

                g_0_x_xxxxxxx_xyyz[k] = g_xxxxxx_xyyz[k] - g_0_x_xxxxxx_xyyz[k] * ab_x + g_0_x_xxxxxx_xxyyz[k];

                g_0_x_xxxxxxx_xyzz[k] = g_xxxxxx_xyzz[k] - g_0_x_xxxxxx_xyzz[k] * ab_x + g_0_x_xxxxxx_xxyzz[k];

                g_0_x_xxxxxxx_xzzz[k] = g_xxxxxx_xzzz[k] - g_0_x_xxxxxx_xzzz[k] * ab_x + g_0_x_xxxxxx_xxzzz[k];

                g_0_x_xxxxxxx_yyyy[k] = g_xxxxxx_yyyy[k] - g_0_x_xxxxxx_yyyy[k] * ab_x + g_0_x_xxxxxx_xyyyy[k];

                g_0_x_xxxxxxx_yyyz[k] = g_xxxxxx_yyyz[k] - g_0_x_xxxxxx_yyyz[k] * ab_x + g_0_x_xxxxxx_xyyyz[k];

                g_0_x_xxxxxxx_yyzz[k] = g_xxxxxx_yyzz[k] - g_0_x_xxxxxx_yyzz[k] * ab_x + g_0_x_xxxxxx_xyyzz[k];

                g_0_x_xxxxxxx_yzzz[k] = g_xxxxxx_yzzz[k] - g_0_x_xxxxxx_yzzz[k] * ab_x + g_0_x_xxxxxx_xyzzz[k];

                g_0_x_xxxxxxx_zzzz[k] = g_xxxxxx_zzzz[k] - g_0_x_xxxxxx_zzzz[k] * ab_x + g_0_x_xxxxxx_xzzzz[k];
            }

            /// Set up 15-30 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxy_xxxx = cbuffer.data(kg_geom_01_off + 15 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xxxy = cbuffer.data(kg_geom_01_off + 16 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xxxz = cbuffer.data(kg_geom_01_off + 17 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xxyy = cbuffer.data(kg_geom_01_off + 18 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xxyz = cbuffer.data(kg_geom_01_off + 19 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xxzz = cbuffer.data(kg_geom_01_off + 20 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xyyy = cbuffer.data(kg_geom_01_off + 21 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xyyz = cbuffer.data(kg_geom_01_off + 22 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xyzz = cbuffer.data(kg_geom_01_off + 23 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_xzzz = cbuffer.data(kg_geom_01_off + 24 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_yyyy = cbuffer.data(kg_geom_01_off + 25 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_yyyz = cbuffer.data(kg_geom_01_off + 26 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_yyzz = cbuffer.data(kg_geom_01_off + 27 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_yzzz = cbuffer.data(kg_geom_01_off + 28 * ccomps * dcomps);

            auto g_0_x_xxxxxxy_zzzz = cbuffer.data(kg_geom_01_off + 29 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_xxxx, g_0_x_xxxxxx_xxxxy, g_0_x_xxxxxx_xxxy, g_0_x_xxxxxx_xxxyy, g_0_x_xxxxxx_xxxyz, g_0_x_xxxxxx_xxxz, g_0_x_xxxxxx_xxyy, g_0_x_xxxxxx_xxyyy, g_0_x_xxxxxx_xxyyz, g_0_x_xxxxxx_xxyz, g_0_x_xxxxxx_xxyzz, g_0_x_xxxxxx_xxzz, g_0_x_xxxxxx_xyyy, g_0_x_xxxxxx_xyyyy, g_0_x_xxxxxx_xyyyz, g_0_x_xxxxxx_xyyz, g_0_x_xxxxxx_xyyzz, g_0_x_xxxxxx_xyzz, g_0_x_xxxxxx_xyzzz, g_0_x_xxxxxx_xzzz, g_0_x_xxxxxx_yyyy, g_0_x_xxxxxx_yyyyy, g_0_x_xxxxxx_yyyyz, g_0_x_xxxxxx_yyyz, g_0_x_xxxxxx_yyyzz, g_0_x_xxxxxx_yyzz, g_0_x_xxxxxx_yyzzz, g_0_x_xxxxxx_yzzz, g_0_x_xxxxxx_yzzzz, g_0_x_xxxxxx_zzzz, g_0_x_xxxxxxy_xxxx, g_0_x_xxxxxxy_xxxy, g_0_x_xxxxxxy_xxxz, g_0_x_xxxxxxy_xxyy, g_0_x_xxxxxxy_xxyz, g_0_x_xxxxxxy_xxzz, g_0_x_xxxxxxy_xyyy, g_0_x_xxxxxxy_xyyz, g_0_x_xxxxxxy_xyzz, g_0_x_xxxxxxy_xzzz, g_0_x_xxxxxxy_yyyy, g_0_x_xxxxxxy_yyyz, g_0_x_xxxxxxy_yyzz, g_0_x_xxxxxxy_yzzz, g_0_x_xxxxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxy_xxxx[k] = -g_0_x_xxxxxx_xxxx[k] * ab_y + g_0_x_xxxxxx_xxxxy[k];

                g_0_x_xxxxxxy_xxxy[k] = -g_0_x_xxxxxx_xxxy[k] * ab_y + g_0_x_xxxxxx_xxxyy[k];

                g_0_x_xxxxxxy_xxxz[k] = -g_0_x_xxxxxx_xxxz[k] * ab_y + g_0_x_xxxxxx_xxxyz[k];

                g_0_x_xxxxxxy_xxyy[k] = -g_0_x_xxxxxx_xxyy[k] * ab_y + g_0_x_xxxxxx_xxyyy[k];

                g_0_x_xxxxxxy_xxyz[k] = -g_0_x_xxxxxx_xxyz[k] * ab_y + g_0_x_xxxxxx_xxyyz[k];

                g_0_x_xxxxxxy_xxzz[k] = -g_0_x_xxxxxx_xxzz[k] * ab_y + g_0_x_xxxxxx_xxyzz[k];

                g_0_x_xxxxxxy_xyyy[k] = -g_0_x_xxxxxx_xyyy[k] * ab_y + g_0_x_xxxxxx_xyyyy[k];

                g_0_x_xxxxxxy_xyyz[k] = -g_0_x_xxxxxx_xyyz[k] * ab_y + g_0_x_xxxxxx_xyyyz[k];

                g_0_x_xxxxxxy_xyzz[k] = -g_0_x_xxxxxx_xyzz[k] * ab_y + g_0_x_xxxxxx_xyyzz[k];

                g_0_x_xxxxxxy_xzzz[k] = -g_0_x_xxxxxx_xzzz[k] * ab_y + g_0_x_xxxxxx_xyzzz[k];

                g_0_x_xxxxxxy_yyyy[k] = -g_0_x_xxxxxx_yyyy[k] * ab_y + g_0_x_xxxxxx_yyyyy[k];

                g_0_x_xxxxxxy_yyyz[k] = -g_0_x_xxxxxx_yyyz[k] * ab_y + g_0_x_xxxxxx_yyyyz[k];

                g_0_x_xxxxxxy_yyzz[k] = -g_0_x_xxxxxx_yyzz[k] * ab_y + g_0_x_xxxxxx_yyyzz[k];

                g_0_x_xxxxxxy_yzzz[k] = -g_0_x_xxxxxx_yzzz[k] * ab_y + g_0_x_xxxxxx_yyzzz[k];

                g_0_x_xxxxxxy_zzzz[k] = -g_0_x_xxxxxx_zzzz[k] * ab_y + g_0_x_xxxxxx_yzzzz[k];
            }

            /// Set up 30-45 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxxz_xxxx = cbuffer.data(kg_geom_01_off + 30 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xxxy = cbuffer.data(kg_geom_01_off + 31 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xxxz = cbuffer.data(kg_geom_01_off + 32 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xxyy = cbuffer.data(kg_geom_01_off + 33 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xxyz = cbuffer.data(kg_geom_01_off + 34 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xxzz = cbuffer.data(kg_geom_01_off + 35 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xyyy = cbuffer.data(kg_geom_01_off + 36 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xyyz = cbuffer.data(kg_geom_01_off + 37 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xyzz = cbuffer.data(kg_geom_01_off + 38 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_xzzz = cbuffer.data(kg_geom_01_off + 39 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_yyyy = cbuffer.data(kg_geom_01_off + 40 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_yyyz = cbuffer.data(kg_geom_01_off + 41 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_yyzz = cbuffer.data(kg_geom_01_off + 42 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_yzzz = cbuffer.data(kg_geom_01_off + 43 * ccomps * dcomps);

            auto g_0_x_xxxxxxz_zzzz = cbuffer.data(kg_geom_01_off + 44 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxx_xxxx, g_0_x_xxxxxx_xxxxz, g_0_x_xxxxxx_xxxy, g_0_x_xxxxxx_xxxyz, g_0_x_xxxxxx_xxxz, g_0_x_xxxxxx_xxxzz, g_0_x_xxxxxx_xxyy, g_0_x_xxxxxx_xxyyz, g_0_x_xxxxxx_xxyz, g_0_x_xxxxxx_xxyzz, g_0_x_xxxxxx_xxzz, g_0_x_xxxxxx_xxzzz, g_0_x_xxxxxx_xyyy, g_0_x_xxxxxx_xyyyz, g_0_x_xxxxxx_xyyz, g_0_x_xxxxxx_xyyzz, g_0_x_xxxxxx_xyzz, g_0_x_xxxxxx_xyzzz, g_0_x_xxxxxx_xzzz, g_0_x_xxxxxx_xzzzz, g_0_x_xxxxxx_yyyy, g_0_x_xxxxxx_yyyyz, g_0_x_xxxxxx_yyyz, g_0_x_xxxxxx_yyyzz, g_0_x_xxxxxx_yyzz, g_0_x_xxxxxx_yyzzz, g_0_x_xxxxxx_yzzz, g_0_x_xxxxxx_yzzzz, g_0_x_xxxxxx_zzzz, g_0_x_xxxxxx_zzzzz, g_0_x_xxxxxxz_xxxx, g_0_x_xxxxxxz_xxxy, g_0_x_xxxxxxz_xxxz, g_0_x_xxxxxxz_xxyy, g_0_x_xxxxxxz_xxyz, g_0_x_xxxxxxz_xxzz, g_0_x_xxxxxxz_xyyy, g_0_x_xxxxxxz_xyyz, g_0_x_xxxxxxz_xyzz, g_0_x_xxxxxxz_xzzz, g_0_x_xxxxxxz_yyyy, g_0_x_xxxxxxz_yyyz, g_0_x_xxxxxxz_yyzz, g_0_x_xxxxxxz_yzzz, g_0_x_xxxxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxxz_xxxx[k] = -g_0_x_xxxxxx_xxxx[k] * ab_z + g_0_x_xxxxxx_xxxxz[k];

                g_0_x_xxxxxxz_xxxy[k] = -g_0_x_xxxxxx_xxxy[k] * ab_z + g_0_x_xxxxxx_xxxyz[k];

                g_0_x_xxxxxxz_xxxz[k] = -g_0_x_xxxxxx_xxxz[k] * ab_z + g_0_x_xxxxxx_xxxzz[k];

                g_0_x_xxxxxxz_xxyy[k] = -g_0_x_xxxxxx_xxyy[k] * ab_z + g_0_x_xxxxxx_xxyyz[k];

                g_0_x_xxxxxxz_xxyz[k] = -g_0_x_xxxxxx_xxyz[k] * ab_z + g_0_x_xxxxxx_xxyzz[k];

                g_0_x_xxxxxxz_xxzz[k] = -g_0_x_xxxxxx_xxzz[k] * ab_z + g_0_x_xxxxxx_xxzzz[k];

                g_0_x_xxxxxxz_xyyy[k] = -g_0_x_xxxxxx_xyyy[k] * ab_z + g_0_x_xxxxxx_xyyyz[k];

                g_0_x_xxxxxxz_xyyz[k] = -g_0_x_xxxxxx_xyyz[k] * ab_z + g_0_x_xxxxxx_xyyzz[k];

                g_0_x_xxxxxxz_xyzz[k] = -g_0_x_xxxxxx_xyzz[k] * ab_z + g_0_x_xxxxxx_xyzzz[k];

                g_0_x_xxxxxxz_xzzz[k] = -g_0_x_xxxxxx_xzzz[k] * ab_z + g_0_x_xxxxxx_xzzzz[k];

                g_0_x_xxxxxxz_yyyy[k] = -g_0_x_xxxxxx_yyyy[k] * ab_z + g_0_x_xxxxxx_yyyyz[k];

                g_0_x_xxxxxxz_yyyz[k] = -g_0_x_xxxxxx_yyyz[k] * ab_z + g_0_x_xxxxxx_yyyzz[k];

                g_0_x_xxxxxxz_yyzz[k] = -g_0_x_xxxxxx_yyzz[k] * ab_z + g_0_x_xxxxxx_yyzzz[k];

                g_0_x_xxxxxxz_yzzz[k] = -g_0_x_xxxxxx_yzzz[k] * ab_z + g_0_x_xxxxxx_yzzzz[k];

                g_0_x_xxxxxxz_zzzz[k] = -g_0_x_xxxxxx_zzzz[k] * ab_z + g_0_x_xxxxxx_zzzzz[k];
            }

            /// Set up 45-60 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxyy_xxxx = cbuffer.data(kg_geom_01_off + 45 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xxxy = cbuffer.data(kg_geom_01_off + 46 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xxxz = cbuffer.data(kg_geom_01_off + 47 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xxyy = cbuffer.data(kg_geom_01_off + 48 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xxyz = cbuffer.data(kg_geom_01_off + 49 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xxzz = cbuffer.data(kg_geom_01_off + 50 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xyyy = cbuffer.data(kg_geom_01_off + 51 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xyyz = cbuffer.data(kg_geom_01_off + 52 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xyzz = cbuffer.data(kg_geom_01_off + 53 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_xzzz = cbuffer.data(kg_geom_01_off + 54 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_yyyy = cbuffer.data(kg_geom_01_off + 55 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_yyyz = cbuffer.data(kg_geom_01_off + 56 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_yyzz = cbuffer.data(kg_geom_01_off + 57 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_yzzz = cbuffer.data(kg_geom_01_off + 58 * ccomps * dcomps);

            auto g_0_x_xxxxxyy_zzzz = cbuffer.data(kg_geom_01_off + 59 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxy_xxxx, g_0_x_xxxxxy_xxxxy, g_0_x_xxxxxy_xxxy, g_0_x_xxxxxy_xxxyy, g_0_x_xxxxxy_xxxyz, g_0_x_xxxxxy_xxxz, g_0_x_xxxxxy_xxyy, g_0_x_xxxxxy_xxyyy, g_0_x_xxxxxy_xxyyz, g_0_x_xxxxxy_xxyz, g_0_x_xxxxxy_xxyzz, g_0_x_xxxxxy_xxzz, g_0_x_xxxxxy_xyyy, g_0_x_xxxxxy_xyyyy, g_0_x_xxxxxy_xyyyz, g_0_x_xxxxxy_xyyz, g_0_x_xxxxxy_xyyzz, g_0_x_xxxxxy_xyzz, g_0_x_xxxxxy_xyzzz, g_0_x_xxxxxy_xzzz, g_0_x_xxxxxy_yyyy, g_0_x_xxxxxy_yyyyy, g_0_x_xxxxxy_yyyyz, g_0_x_xxxxxy_yyyz, g_0_x_xxxxxy_yyyzz, g_0_x_xxxxxy_yyzz, g_0_x_xxxxxy_yyzzz, g_0_x_xxxxxy_yzzz, g_0_x_xxxxxy_yzzzz, g_0_x_xxxxxy_zzzz, g_0_x_xxxxxyy_xxxx, g_0_x_xxxxxyy_xxxy, g_0_x_xxxxxyy_xxxz, g_0_x_xxxxxyy_xxyy, g_0_x_xxxxxyy_xxyz, g_0_x_xxxxxyy_xxzz, g_0_x_xxxxxyy_xyyy, g_0_x_xxxxxyy_xyyz, g_0_x_xxxxxyy_xyzz, g_0_x_xxxxxyy_xzzz, g_0_x_xxxxxyy_yyyy, g_0_x_xxxxxyy_yyyz, g_0_x_xxxxxyy_yyzz, g_0_x_xxxxxyy_yzzz, g_0_x_xxxxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxyy_xxxx[k] = -g_0_x_xxxxxy_xxxx[k] * ab_y + g_0_x_xxxxxy_xxxxy[k];

                g_0_x_xxxxxyy_xxxy[k] = -g_0_x_xxxxxy_xxxy[k] * ab_y + g_0_x_xxxxxy_xxxyy[k];

                g_0_x_xxxxxyy_xxxz[k] = -g_0_x_xxxxxy_xxxz[k] * ab_y + g_0_x_xxxxxy_xxxyz[k];

                g_0_x_xxxxxyy_xxyy[k] = -g_0_x_xxxxxy_xxyy[k] * ab_y + g_0_x_xxxxxy_xxyyy[k];

                g_0_x_xxxxxyy_xxyz[k] = -g_0_x_xxxxxy_xxyz[k] * ab_y + g_0_x_xxxxxy_xxyyz[k];

                g_0_x_xxxxxyy_xxzz[k] = -g_0_x_xxxxxy_xxzz[k] * ab_y + g_0_x_xxxxxy_xxyzz[k];

                g_0_x_xxxxxyy_xyyy[k] = -g_0_x_xxxxxy_xyyy[k] * ab_y + g_0_x_xxxxxy_xyyyy[k];

                g_0_x_xxxxxyy_xyyz[k] = -g_0_x_xxxxxy_xyyz[k] * ab_y + g_0_x_xxxxxy_xyyyz[k];

                g_0_x_xxxxxyy_xyzz[k] = -g_0_x_xxxxxy_xyzz[k] * ab_y + g_0_x_xxxxxy_xyyzz[k];

                g_0_x_xxxxxyy_xzzz[k] = -g_0_x_xxxxxy_xzzz[k] * ab_y + g_0_x_xxxxxy_xyzzz[k];

                g_0_x_xxxxxyy_yyyy[k] = -g_0_x_xxxxxy_yyyy[k] * ab_y + g_0_x_xxxxxy_yyyyy[k];

                g_0_x_xxxxxyy_yyyz[k] = -g_0_x_xxxxxy_yyyz[k] * ab_y + g_0_x_xxxxxy_yyyyz[k];

                g_0_x_xxxxxyy_yyzz[k] = -g_0_x_xxxxxy_yyzz[k] * ab_y + g_0_x_xxxxxy_yyyzz[k];

                g_0_x_xxxxxyy_yzzz[k] = -g_0_x_xxxxxy_yzzz[k] * ab_y + g_0_x_xxxxxy_yyzzz[k];

                g_0_x_xxxxxyy_zzzz[k] = -g_0_x_xxxxxy_zzzz[k] * ab_y + g_0_x_xxxxxy_yzzzz[k];
            }

            /// Set up 60-75 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxyz_xxxx = cbuffer.data(kg_geom_01_off + 60 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xxxy = cbuffer.data(kg_geom_01_off + 61 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xxxz = cbuffer.data(kg_geom_01_off + 62 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xxyy = cbuffer.data(kg_geom_01_off + 63 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xxyz = cbuffer.data(kg_geom_01_off + 64 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xxzz = cbuffer.data(kg_geom_01_off + 65 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xyyy = cbuffer.data(kg_geom_01_off + 66 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xyyz = cbuffer.data(kg_geom_01_off + 67 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xyzz = cbuffer.data(kg_geom_01_off + 68 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_xzzz = cbuffer.data(kg_geom_01_off + 69 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_yyyy = cbuffer.data(kg_geom_01_off + 70 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_yyyz = cbuffer.data(kg_geom_01_off + 71 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_yyzz = cbuffer.data(kg_geom_01_off + 72 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_yzzz = cbuffer.data(kg_geom_01_off + 73 * ccomps * dcomps);

            auto g_0_x_xxxxxyz_zzzz = cbuffer.data(kg_geom_01_off + 74 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxyz_xxxx, g_0_x_xxxxxyz_xxxy, g_0_x_xxxxxyz_xxxz, g_0_x_xxxxxyz_xxyy, g_0_x_xxxxxyz_xxyz, g_0_x_xxxxxyz_xxzz, g_0_x_xxxxxyz_xyyy, g_0_x_xxxxxyz_xyyz, g_0_x_xxxxxyz_xyzz, g_0_x_xxxxxyz_xzzz, g_0_x_xxxxxyz_yyyy, g_0_x_xxxxxyz_yyyz, g_0_x_xxxxxyz_yyzz, g_0_x_xxxxxyz_yzzz, g_0_x_xxxxxyz_zzzz, g_0_x_xxxxxz_xxxx, g_0_x_xxxxxz_xxxxy, g_0_x_xxxxxz_xxxy, g_0_x_xxxxxz_xxxyy, g_0_x_xxxxxz_xxxyz, g_0_x_xxxxxz_xxxz, g_0_x_xxxxxz_xxyy, g_0_x_xxxxxz_xxyyy, g_0_x_xxxxxz_xxyyz, g_0_x_xxxxxz_xxyz, g_0_x_xxxxxz_xxyzz, g_0_x_xxxxxz_xxzz, g_0_x_xxxxxz_xyyy, g_0_x_xxxxxz_xyyyy, g_0_x_xxxxxz_xyyyz, g_0_x_xxxxxz_xyyz, g_0_x_xxxxxz_xyyzz, g_0_x_xxxxxz_xyzz, g_0_x_xxxxxz_xyzzz, g_0_x_xxxxxz_xzzz, g_0_x_xxxxxz_yyyy, g_0_x_xxxxxz_yyyyy, g_0_x_xxxxxz_yyyyz, g_0_x_xxxxxz_yyyz, g_0_x_xxxxxz_yyyzz, g_0_x_xxxxxz_yyzz, g_0_x_xxxxxz_yyzzz, g_0_x_xxxxxz_yzzz, g_0_x_xxxxxz_yzzzz, g_0_x_xxxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxyz_xxxx[k] = -g_0_x_xxxxxz_xxxx[k] * ab_y + g_0_x_xxxxxz_xxxxy[k];

                g_0_x_xxxxxyz_xxxy[k] = -g_0_x_xxxxxz_xxxy[k] * ab_y + g_0_x_xxxxxz_xxxyy[k];

                g_0_x_xxxxxyz_xxxz[k] = -g_0_x_xxxxxz_xxxz[k] * ab_y + g_0_x_xxxxxz_xxxyz[k];

                g_0_x_xxxxxyz_xxyy[k] = -g_0_x_xxxxxz_xxyy[k] * ab_y + g_0_x_xxxxxz_xxyyy[k];

                g_0_x_xxxxxyz_xxyz[k] = -g_0_x_xxxxxz_xxyz[k] * ab_y + g_0_x_xxxxxz_xxyyz[k];

                g_0_x_xxxxxyz_xxzz[k] = -g_0_x_xxxxxz_xxzz[k] * ab_y + g_0_x_xxxxxz_xxyzz[k];

                g_0_x_xxxxxyz_xyyy[k] = -g_0_x_xxxxxz_xyyy[k] * ab_y + g_0_x_xxxxxz_xyyyy[k];

                g_0_x_xxxxxyz_xyyz[k] = -g_0_x_xxxxxz_xyyz[k] * ab_y + g_0_x_xxxxxz_xyyyz[k];

                g_0_x_xxxxxyz_xyzz[k] = -g_0_x_xxxxxz_xyzz[k] * ab_y + g_0_x_xxxxxz_xyyzz[k];

                g_0_x_xxxxxyz_xzzz[k] = -g_0_x_xxxxxz_xzzz[k] * ab_y + g_0_x_xxxxxz_xyzzz[k];

                g_0_x_xxxxxyz_yyyy[k] = -g_0_x_xxxxxz_yyyy[k] * ab_y + g_0_x_xxxxxz_yyyyy[k];

                g_0_x_xxxxxyz_yyyz[k] = -g_0_x_xxxxxz_yyyz[k] * ab_y + g_0_x_xxxxxz_yyyyz[k];

                g_0_x_xxxxxyz_yyzz[k] = -g_0_x_xxxxxz_yyzz[k] * ab_y + g_0_x_xxxxxz_yyyzz[k];

                g_0_x_xxxxxyz_yzzz[k] = -g_0_x_xxxxxz_yzzz[k] * ab_y + g_0_x_xxxxxz_yyzzz[k];

                g_0_x_xxxxxyz_zzzz[k] = -g_0_x_xxxxxz_zzzz[k] * ab_y + g_0_x_xxxxxz_yzzzz[k];
            }

            /// Set up 75-90 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxxzz_xxxx = cbuffer.data(kg_geom_01_off + 75 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xxxy = cbuffer.data(kg_geom_01_off + 76 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xxxz = cbuffer.data(kg_geom_01_off + 77 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xxyy = cbuffer.data(kg_geom_01_off + 78 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xxyz = cbuffer.data(kg_geom_01_off + 79 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xxzz = cbuffer.data(kg_geom_01_off + 80 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xyyy = cbuffer.data(kg_geom_01_off + 81 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xyyz = cbuffer.data(kg_geom_01_off + 82 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xyzz = cbuffer.data(kg_geom_01_off + 83 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_xzzz = cbuffer.data(kg_geom_01_off + 84 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_yyyy = cbuffer.data(kg_geom_01_off + 85 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_yyyz = cbuffer.data(kg_geom_01_off + 86 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_yyzz = cbuffer.data(kg_geom_01_off + 87 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_yzzz = cbuffer.data(kg_geom_01_off + 88 * ccomps * dcomps);

            auto g_0_x_xxxxxzz_zzzz = cbuffer.data(kg_geom_01_off + 89 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxxz_xxxx, g_0_x_xxxxxz_xxxxz, g_0_x_xxxxxz_xxxy, g_0_x_xxxxxz_xxxyz, g_0_x_xxxxxz_xxxz, g_0_x_xxxxxz_xxxzz, g_0_x_xxxxxz_xxyy, g_0_x_xxxxxz_xxyyz, g_0_x_xxxxxz_xxyz, g_0_x_xxxxxz_xxyzz, g_0_x_xxxxxz_xxzz, g_0_x_xxxxxz_xxzzz, g_0_x_xxxxxz_xyyy, g_0_x_xxxxxz_xyyyz, g_0_x_xxxxxz_xyyz, g_0_x_xxxxxz_xyyzz, g_0_x_xxxxxz_xyzz, g_0_x_xxxxxz_xyzzz, g_0_x_xxxxxz_xzzz, g_0_x_xxxxxz_xzzzz, g_0_x_xxxxxz_yyyy, g_0_x_xxxxxz_yyyyz, g_0_x_xxxxxz_yyyz, g_0_x_xxxxxz_yyyzz, g_0_x_xxxxxz_yyzz, g_0_x_xxxxxz_yyzzz, g_0_x_xxxxxz_yzzz, g_0_x_xxxxxz_yzzzz, g_0_x_xxxxxz_zzzz, g_0_x_xxxxxz_zzzzz, g_0_x_xxxxxzz_xxxx, g_0_x_xxxxxzz_xxxy, g_0_x_xxxxxzz_xxxz, g_0_x_xxxxxzz_xxyy, g_0_x_xxxxxzz_xxyz, g_0_x_xxxxxzz_xxzz, g_0_x_xxxxxzz_xyyy, g_0_x_xxxxxzz_xyyz, g_0_x_xxxxxzz_xyzz, g_0_x_xxxxxzz_xzzz, g_0_x_xxxxxzz_yyyy, g_0_x_xxxxxzz_yyyz, g_0_x_xxxxxzz_yyzz, g_0_x_xxxxxzz_yzzz, g_0_x_xxxxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxxzz_xxxx[k] = -g_0_x_xxxxxz_xxxx[k] * ab_z + g_0_x_xxxxxz_xxxxz[k];

                g_0_x_xxxxxzz_xxxy[k] = -g_0_x_xxxxxz_xxxy[k] * ab_z + g_0_x_xxxxxz_xxxyz[k];

                g_0_x_xxxxxzz_xxxz[k] = -g_0_x_xxxxxz_xxxz[k] * ab_z + g_0_x_xxxxxz_xxxzz[k];

                g_0_x_xxxxxzz_xxyy[k] = -g_0_x_xxxxxz_xxyy[k] * ab_z + g_0_x_xxxxxz_xxyyz[k];

                g_0_x_xxxxxzz_xxyz[k] = -g_0_x_xxxxxz_xxyz[k] * ab_z + g_0_x_xxxxxz_xxyzz[k];

                g_0_x_xxxxxzz_xxzz[k] = -g_0_x_xxxxxz_xxzz[k] * ab_z + g_0_x_xxxxxz_xxzzz[k];

                g_0_x_xxxxxzz_xyyy[k] = -g_0_x_xxxxxz_xyyy[k] * ab_z + g_0_x_xxxxxz_xyyyz[k];

                g_0_x_xxxxxzz_xyyz[k] = -g_0_x_xxxxxz_xyyz[k] * ab_z + g_0_x_xxxxxz_xyyzz[k];

                g_0_x_xxxxxzz_xyzz[k] = -g_0_x_xxxxxz_xyzz[k] * ab_z + g_0_x_xxxxxz_xyzzz[k];

                g_0_x_xxxxxzz_xzzz[k] = -g_0_x_xxxxxz_xzzz[k] * ab_z + g_0_x_xxxxxz_xzzzz[k];

                g_0_x_xxxxxzz_yyyy[k] = -g_0_x_xxxxxz_yyyy[k] * ab_z + g_0_x_xxxxxz_yyyyz[k];

                g_0_x_xxxxxzz_yyyz[k] = -g_0_x_xxxxxz_yyyz[k] * ab_z + g_0_x_xxxxxz_yyyzz[k];

                g_0_x_xxxxxzz_yyzz[k] = -g_0_x_xxxxxz_yyzz[k] * ab_z + g_0_x_xxxxxz_yyzzz[k];

                g_0_x_xxxxxzz_yzzz[k] = -g_0_x_xxxxxz_yzzz[k] * ab_z + g_0_x_xxxxxz_yzzzz[k];

                g_0_x_xxxxxzz_zzzz[k] = -g_0_x_xxxxxz_zzzz[k] * ab_z + g_0_x_xxxxxz_zzzzz[k];
            }

            /// Set up 90-105 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyyy_xxxx = cbuffer.data(kg_geom_01_off + 90 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xxxy = cbuffer.data(kg_geom_01_off + 91 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xxxz = cbuffer.data(kg_geom_01_off + 92 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xxyy = cbuffer.data(kg_geom_01_off + 93 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xxyz = cbuffer.data(kg_geom_01_off + 94 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xxzz = cbuffer.data(kg_geom_01_off + 95 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xyyy = cbuffer.data(kg_geom_01_off + 96 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xyyz = cbuffer.data(kg_geom_01_off + 97 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xyzz = cbuffer.data(kg_geom_01_off + 98 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_xzzz = cbuffer.data(kg_geom_01_off + 99 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_yyyy = cbuffer.data(kg_geom_01_off + 100 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_yyyz = cbuffer.data(kg_geom_01_off + 101 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_yyzz = cbuffer.data(kg_geom_01_off + 102 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_yzzz = cbuffer.data(kg_geom_01_off + 103 * ccomps * dcomps);

            auto g_0_x_xxxxyyy_zzzz = cbuffer.data(kg_geom_01_off + 104 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyy_xxxx, g_0_x_xxxxyy_xxxxy, g_0_x_xxxxyy_xxxy, g_0_x_xxxxyy_xxxyy, g_0_x_xxxxyy_xxxyz, g_0_x_xxxxyy_xxxz, g_0_x_xxxxyy_xxyy, g_0_x_xxxxyy_xxyyy, g_0_x_xxxxyy_xxyyz, g_0_x_xxxxyy_xxyz, g_0_x_xxxxyy_xxyzz, g_0_x_xxxxyy_xxzz, g_0_x_xxxxyy_xyyy, g_0_x_xxxxyy_xyyyy, g_0_x_xxxxyy_xyyyz, g_0_x_xxxxyy_xyyz, g_0_x_xxxxyy_xyyzz, g_0_x_xxxxyy_xyzz, g_0_x_xxxxyy_xyzzz, g_0_x_xxxxyy_xzzz, g_0_x_xxxxyy_yyyy, g_0_x_xxxxyy_yyyyy, g_0_x_xxxxyy_yyyyz, g_0_x_xxxxyy_yyyz, g_0_x_xxxxyy_yyyzz, g_0_x_xxxxyy_yyzz, g_0_x_xxxxyy_yyzzz, g_0_x_xxxxyy_yzzz, g_0_x_xxxxyy_yzzzz, g_0_x_xxxxyy_zzzz, g_0_x_xxxxyyy_xxxx, g_0_x_xxxxyyy_xxxy, g_0_x_xxxxyyy_xxxz, g_0_x_xxxxyyy_xxyy, g_0_x_xxxxyyy_xxyz, g_0_x_xxxxyyy_xxzz, g_0_x_xxxxyyy_xyyy, g_0_x_xxxxyyy_xyyz, g_0_x_xxxxyyy_xyzz, g_0_x_xxxxyyy_xzzz, g_0_x_xxxxyyy_yyyy, g_0_x_xxxxyyy_yyyz, g_0_x_xxxxyyy_yyzz, g_0_x_xxxxyyy_yzzz, g_0_x_xxxxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyyy_xxxx[k] = -g_0_x_xxxxyy_xxxx[k] * ab_y + g_0_x_xxxxyy_xxxxy[k];

                g_0_x_xxxxyyy_xxxy[k] = -g_0_x_xxxxyy_xxxy[k] * ab_y + g_0_x_xxxxyy_xxxyy[k];

                g_0_x_xxxxyyy_xxxz[k] = -g_0_x_xxxxyy_xxxz[k] * ab_y + g_0_x_xxxxyy_xxxyz[k];

                g_0_x_xxxxyyy_xxyy[k] = -g_0_x_xxxxyy_xxyy[k] * ab_y + g_0_x_xxxxyy_xxyyy[k];

                g_0_x_xxxxyyy_xxyz[k] = -g_0_x_xxxxyy_xxyz[k] * ab_y + g_0_x_xxxxyy_xxyyz[k];

                g_0_x_xxxxyyy_xxzz[k] = -g_0_x_xxxxyy_xxzz[k] * ab_y + g_0_x_xxxxyy_xxyzz[k];

                g_0_x_xxxxyyy_xyyy[k] = -g_0_x_xxxxyy_xyyy[k] * ab_y + g_0_x_xxxxyy_xyyyy[k];

                g_0_x_xxxxyyy_xyyz[k] = -g_0_x_xxxxyy_xyyz[k] * ab_y + g_0_x_xxxxyy_xyyyz[k];

                g_0_x_xxxxyyy_xyzz[k] = -g_0_x_xxxxyy_xyzz[k] * ab_y + g_0_x_xxxxyy_xyyzz[k];

                g_0_x_xxxxyyy_xzzz[k] = -g_0_x_xxxxyy_xzzz[k] * ab_y + g_0_x_xxxxyy_xyzzz[k];

                g_0_x_xxxxyyy_yyyy[k] = -g_0_x_xxxxyy_yyyy[k] * ab_y + g_0_x_xxxxyy_yyyyy[k];

                g_0_x_xxxxyyy_yyyz[k] = -g_0_x_xxxxyy_yyyz[k] * ab_y + g_0_x_xxxxyy_yyyyz[k];

                g_0_x_xxxxyyy_yyzz[k] = -g_0_x_xxxxyy_yyzz[k] * ab_y + g_0_x_xxxxyy_yyyzz[k];

                g_0_x_xxxxyyy_yzzz[k] = -g_0_x_xxxxyy_yzzz[k] * ab_y + g_0_x_xxxxyy_yyzzz[k];

                g_0_x_xxxxyyy_zzzz[k] = -g_0_x_xxxxyy_zzzz[k] * ab_y + g_0_x_xxxxyy_yzzzz[k];
            }

            /// Set up 105-120 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyyz_xxxx = cbuffer.data(kg_geom_01_off + 105 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xxxy = cbuffer.data(kg_geom_01_off + 106 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xxxz = cbuffer.data(kg_geom_01_off + 107 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xxyy = cbuffer.data(kg_geom_01_off + 108 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xxyz = cbuffer.data(kg_geom_01_off + 109 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xxzz = cbuffer.data(kg_geom_01_off + 110 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xyyy = cbuffer.data(kg_geom_01_off + 111 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xyyz = cbuffer.data(kg_geom_01_off + 112 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xyzz = cbuffer.data(kg_geom_01_off + 113 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_xzzz = cbuffer.data(kg_geom_01_off + 114 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_yyyy = cbuffer.data(kg_geom_01_off + 115 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_yyyz = cbuffer.data(kg_geom_01_off + 116 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_yyzz = cbuffer.data(kg_geom_01_off + 117 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_yzzz = cbuffer.data(kg_geom_01_off + 118 * ccomps * dcomps);

            auto g_0_x_xxxxyyz_zzzz = cbuffer.data(kg_geom_01_off + 119 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyyz_xxxx, g_0_x_xxxxyyz_xxxy, g_0_x_xxxxyyz_xxxz, g_0_x_xxxxyyz_xxyy, g_0_x_xxxxyyz_xxyz, g_0_x_xxxxyyz_xxzz, g_0_x_xxxxyyz_xyyy, g_0_x_xxxxyyz_xyyz, g_0_x_xxxxyyz_xyzz, g_0_x_xxxxyyz_xzzz, g_0_x_xxxxyyz_yyyy, g_0_x_xxxxyyz_yyyz, g_0_x_xxxxyyz_yyzz, g_0_x_xxxxyyz_yzzz, g_0_x_xxxxyyz_zzzz, g_0_x_xxxxyz_xxxx, g_0_x_xxxxyz_xxxxy, g_0_x_xxxxyz_xxxy, g_0_x_xxxxyz_xxxyy, g_0_x_xxxxyz_xxxyz, g_0_x_xxxxyz_xxxz, g_0_x_xxxxyz_xxyy, g_0_x_xxxxyz_xxyyy, g_0_x_xxxxyz_xxyyz, g_0_x_xxxxyz_xxyz, g_0_x_xxxxyz_xxyzz, g_0_x_xxxxyz_xxzz, g_0_x_xxxxyz_xyyy, g_0_x_xxxxyz_xyyyy, g_0_x_xxxxyz_xyyyz, g_0_x_xxxxyz_xyyz, g_0_x_xxxxyz_xyyzz, g_0_x_xxxxyz_xyzz, g_0_x_xxxxyz_xyzzz, g_0_x_xxxxyz_xzzz, g_0_x_xxxxyz_yyyy, g_0_x_xxxxyz_yyyyy, g_0_x_xxxxyz_yyyyz, g_0_x_xxxxyz_yyyz, g_0_x_xxxxyz_yyyzz, g_0_x_xxxxyz_yyzz, g_0_x_xxxxyz_yyzzz, g_0_x_xxxxyz_yzzz, g_0_x_xxxxyz_yzzzz, g_0_x_xxxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyyz_xxxx[k] = -g_0_x_xxxxyz_xxxx[k] * ab_y + g_0_x_xxxxyz_xxxxy[k];

                g_0_x_xxxxyyz_xxxy[k] = -g_0_x_xxxxyz_xxxy[k] * ab_y + g_0_x_xxxxyz_xxxyy[k];

                g_0_x_xxxxyyz_xxxz[k] = -g_0_x_xxxxyz_xxxz[k] * ab_y + g_0_x_xxxxyz_xxxyz[k];

                g_0_x_xxxxyyz_xxyy[k] = -g_0_x_xxxxyz_xxyy[k] * ab_y + g_0_x_xxxxyz_xxyyy[k];

                g_0_x_xxxxyyz_xxyz[k] = -g_0_x_xxxxyz_xxyz[k] * ab_y + g_0_x_xxxxyz_xxyyz[k];

                g_0_x_xxxxyyz_xxzz[k] = -g_0_x_xxxxyz_xxzz[k] * ab_y + g_0_x_xxxxyz_xxyzz[k];

                g_0_x_xxxxyyz_xyyy[k] = -g_0_x_xxxxyz_xyyy[k] * ab_y + g_0_x_xxxxyz_xyyyy[k];

                g_0_x_xxxxyyz_xyyz[k] = -g_0_x_xxxxyz_xyyz[k] * ab_y + g_0_x_xxxxyz_xyyyz[k];

                g_0_x_xxxxyyz_xyzz[k] = -g_0_x_xxxxyz_xyzz[k] * ab_y + g_0_x_xxxxyz_xyyzz[k];

                g_0_x_xxxxyyz_xzzz[k] = -g_0_x_xxxxyz_xzzz[k] * ab_y + g_0_x_xxxxyz_xyzzz[k];

                g_0_x_xxxxyyz_yyyy[k] = -g_0_x_xxxxyz_yyyy[k] * ab_y + g_0_x_xxxxyz_yyyyy[k];

                g_0_x_xxxxyyz_yyyz[k] = -g_0_x_xxxxyz_yyyz[k] * ab_y + g_0_x_xxxxyz_yyyyz[k];

                g_0_x_xxxxyyz_yyzz[k] = -g_0_x_xxxxyz_yyzz[k] * ab_y + g_0_x_xxxxyz_yyyzz[k];

                g_0_x_xxxxyyz_yzzz[k] = -g_0_x_xxxxyz_yzzz[k] * ab_y + g_0_x_xxxxyz_yyzzz[k];

                g_0_x_xxxxyyz_zzzz[k] = -g_0_x_xxxxyz_zzzz[k] * ab_y + g_0_x_xxxxyz_yzzzz[k];
            }

            /// Set up 120-135 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxyzz_xxxx = cbuffer.data(kg_geom_01_off + 120 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xxxy = cbuffer.data(kg_geom_01_off + 121 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xxxz = cbuffer.data(kg_geom_01_off + 122 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xxyy = cbuffer.data(kg_geom_01_off + 123 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xxyz = cbuffer.data(kg_geom_01_off + 124 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xxzz = cbuffer.data(kg_geom_01_off + 125 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xyyy = cbuffer.data(kg_geom_01_off + 126 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xyyz = cbuffer.data(kg_geom_01_off + 127 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xyzz = cbuffer.data(kg_geom_01_off + 128 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_xzzz = cbuffer.data(kg_geom_01_off + 129 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_yyyy = cbuffer.data(kg_geom_01_off + 130 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_yyyz = cbuffer.data(kg_geom_01_off + 131 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_yyzz = cbuffer.data(kg_geom_01_off + 132 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_yzzz = cbuffer.data(kg_geom_01_off + 133 * ccomps * dcomps);

            auto g_0_x_xxxxyzz_zzzz = cbuffer.data(kg_geom_01_off + 134 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxyzz_xxxx, g_0_x_xxxxyzz_xxxy, g_0_x_xxxxyzz_xxxz, g_0_x_xxxxyzz_xxyy, g_0_x_xxxxyzz_xxyz, g_0_x_xxxxyzz_xxzz, g_0_x_xxxxyzz_xyyy, g_0_x_xxxxyzz_xyyz, g_0_x_xxxxyzz_xyzz, g_0_x_xxxxyzz_xzzz, g_0_x_xxxxyzz_yyyy, g_0_x_xxxxyzz_yyyz, g_0_x_xxxxyzz_yyzz, g_0_x_xxxxyzz_yzzz, g_0_x_xxxxyzz_zzzz, g_0_x_xxxxzz_xxxx, g_0_x_xxxxzz_xxxxy, g_0_x_xxxxzz_xxxy, g_0_x_xxxxzz_xxxyy, g_0_x_xxxxzz_xxxyz, g_0_x_xxxxzz_xxxz, g_0_x_xxxxzz_xxyy, g_0_x_xxxxzz_xxyyy, g_0_x_xxxxzz_xxyyz, g_0_x_xxxxzz_xxyz, g_0_x_xxxxzz_xxyzz, g_0_x_xxxxzz_xxzz, g_0_x_xxxxzz_xyyy, g_0_x_xxxxzz_xyyyy, g_0_x_xxxxzz_xyyyz, g_0_x_xxxxzz_xyyz, g_0_x_xxxxzz_xyyzz, g_0_x_xxxxzz_xyzz, g_0_x_xxxxzz_xyzzz, g_0_x_xxxxzz_xzzz, g_0_x_xxxxzz_yyyy, g_0_x_xxxxzz_yyyyy, g_0_x_xxxxzz_yyyyz, g_0_x_xxxxzz_yyyz, g_0_x_xxxxzz_yyyzz, g_0_x_xxxxzz_yyzz, g_0_x_xxxxzz_yyzzz, g_0_x_xxxxzz_yzzz, g_0_x_xxxxzz_yzzzz, g_0_x_xxxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxyzz_xxxx[k] = -g_0_x_xxxxzz_xxxx[k] * ab_y + g_0_x_xxxxzz_xxxxy[k];

                g_0_x_xxxxyzz_xxxy[k] = -g_0_x_xxxxzz_xxxy[k] * ab_y + g_0_x_xxxxzz_xxxyy[k];

                g_0_x_xxxxyzz_xxxz[k] = -g_0_x_xxxxzz_xxxz[k] * ab_y + g_0_x_xxxxzz_xxxyz[k];

                g_0_x_xxxxyzz_xxyy[k] = -g_0_x_xxxxzz_xxyy[k] * ab_y + g_0_x_xxxxzz_xxyyy[k];

                g_0_x_xxxxyzz_xxyz[k] = -g_0_x_xxxxzz_xxyz[k] * ab_y + g_0_x_xxxxzz_xxyyz[k];

                g_0_x_xxxxyzz_xxzz[k] = -g_0_x_xxxxzz_xxzz[k] * ab_y + g_0_x_xxxxzz_xxyzz[k];

                g_0_x_xxxxyzz_xyyy[k] = -g_0_x_xxxxzz_xyyy[k] * ab_y + g_0_x_xxxxzz_xyyyy[k];

                g_0_x_xxxxyzz_xyyz[k] = -g_0_x_xxxxzz_xyyz[k] * ab_y + g_0_x_xxxxzz_xyyyz[k];

                g_0_x_xxxxyzz_xyzz[k] = -g_0_x_xxxxzz_xyzz[k] * ab_y + g_0_x_xxxxzz_xyyzz[k];

                g_0_x_xxxxyzz_xzzz[k] = -g_0_x_xxxxzz_xzzz[k] * ab_y + g_0_x_xxxxzz_xyzzz[k];

                g_0_x_xxxxyzz_yyyy[k] = -g_0_x_xxxxzz_yyyy[k] * ab_y + g_0_x_xxxxzz_yyyyy[k];

                g_0_x_xxxxyzz_yyyz[k] = -g_0_x_xxxxzz_yyyz[k] * ab_y + g_0_x_xxxxzz_yyyyz[k];

                g_0_x_xxxxyzz_yyzz[k] = -g_0_x_xxxxzz_yyzz[k] * ab_y + g_0_x_xxxxzz_yyyzz[k];

                g_0_x_xxxxyzz_yzzz[k] = -g_0_x_xxxxzz_yzzz[k] * ab_y + g_0_x_xxxxzz_yyzzz[k];

                g_0_x_xxxxyzz_zzzz[k] = -g_0_x_xxxxzz_zzzz[k] * ab_y + g_0_x_xxxxzz_yzzzz[k];
            }

            /// Set up 135-150 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxxzzz_xxxx = cbuffer.data(kg_geom_01_off + 135 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xxxy = cbuffer.data(kg_geom_01_off + 136 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xxxz = cbuffer.data(kg_geom_01_off + 137 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xxyy = cbuffer.data(kg_geom_01_off + 138 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xxyz = cbuffer.data(kg_geom_01_off + 139 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xxzz = cbuffer.data(kg_geom_01_off + 140 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xyyy = cbuffer.data(kg_geom_01_off + 141 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xyyz = cbuffer.data(kg_geom_01_off + 142 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xyzz = cbuffer.data(kg_geom_01_off + 143 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_xzzz = cbuffer.data(kg_geom_01_off + 144 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_yyyy = cbuffer.data(kg_geom_01_off + 145 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_yyyz = cbuffer.data(kg_geom_01_off + 146 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_yyzz = cbuffer.data(kg_geom_01_off + 147 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_yzzz = cbuffer.data(kg_geom_01_off + 148 * ccomps * dcomps);

            auto g_0_x_xxxxzzz_zzzz = cbuffer.data(kg_geom_01_off + 149 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxxzz_xxxx, g_0_x_xxxxzz_xxxxz, g_0_x_xxxxzz_xxxy, g_0_x_xxxxzz_xxxyz, g_0_x_xxxxzz_xxxz, g_0_x_xxxxzz_xxxzz, g_0_x_xxxxzz_xxyy, g_0_x_xxxxzz_xxyyz, g_0_x_xxxxzz_xxyz, g_0_x_xxxxzz_xxyzz, g_0_x_xxxxzz_xxzz, g_0_x_xxxxzz_xxzzz, g_0_x_xxxxzz_xyyy, g_0_x_xxxxzz_xyyyz, g_0_x_xxxxzz_xyyz, g_0_x_xxxxzz_xyyzz, g_0_x_xxxxzz_xyzz, g_0_x_xxxxzz_xyzzz, g_0_x_xxxxzz_xzzz, g_0_x_xxxxzz_xzzzz, g_0_x_xxxxzz_yyyy, g_0_x_xxxxzz_yyyyz, g_0_x_xxxxzz_yyyz, g_0_x_xxxxzz_yyyzz, g_0_x_xxxxzz_yyzz, g_0_x_xxxxzz_yyzzz, g_0_x_xxxxzz_yzzz, g_0_x_xxxxzz_yzzzz, g_0_x_xxxxzz_zzzz, g_0_x_xxxxzz_zzzzz, g_0_x_xxxxzzz_xxxx, g_0_x_xxxxzzz_xxxy, g_0_x_xxxxzzz_xxxz, g_0_x_xxxxzzz_xxyy, g_0_x_xxxxzzz_xxyz, g_0_x_xxxxzzz_xxzz, g_0_x_xxxxzzz_xyyy, g_0_x_xxxxzzz_xyyz, g_0_x_xxxxzzz_xyzz, g_0_x_xxxxzzz_xzzz, g_0_x_xxxxzzz_yyyy, g_0_x_xxxxzzz_yyyz, g_0_x_xxxxzzz_yyzz, g_0_x_xxxxzzz_yzzz, g_0_x_xxxxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxxzzz_xxxx[k] = -g_0_x_xxxxzz_xxxx[k] * ab_z + g_0_x_xxxxzz_xxxxz[k];

                g_0_x_xxxxzzz_xxxy[k] = -g_0_x_xxxxzz_xxxy[k] * ab_z + g_0_x_xxxxzz_xxxyz[k];

                g_0_x_xxxxzzz_xxxz[k] = -g_0_x_xxxxzz_xxxz[k] * ab_z + g_0_x_xxxxzz_xxxzz[k];

                g_0_x_xxxxzzz_xxyy[k] = -g_0_x_xxxxzz_xxyy[k] * ab_z + g_0_x_xxxxzz_xxyyz[k];

                g_0_x_xxxxzzz_xxyz[k] = -g_0_x_xxxxzz_xxyz[k] * ab_z + g_0_x_xxxxzz_xxyzz[k];

                g_0_x_xxxxzzz_xxzz[k] = -g_0_x_xxxxzz_xxzz[k] * ab_z + g_0_x_xxxxzz_xxzzz[k];

                g_0_x_xxxxzzz_xyyy[k] = -g_0_x_xxxxzz_xyyy[k] * ab_z + g_0_x_xxxxzz_xyyyz[k];

                g_0_x_xxxxzzz_xyyz[k] = -g_0_x_xxxxzz_xyyz[k] * ab_z + g_0_x_xxxxzz_xyyzz[k];

                g_0_x_xxxxzzz_xyzz[k] = -g_0_x_xxxxzz_xyzz[k] * ab_z + g_0_x_xxxxzz_xyzzz[k];

                g_0_x_xxxxzzz_xzzz[k] = -g_0_x_xxxxzz_xzzz[k] * ab_z + g_0_x_xxxxzz_xzzzz[k];

                g_0_x_xxxxzzz_yyyy[k] = -g_0_x_xxxxzz_yyyy[k] * ab_z + g_0_x_xxxxzz_yyyyz[k];

                g_0_x_xxxxzzz_yyyz[k] = -g_0_x_xxxxzz_yyyz[k] * ab_z + g_0_x_xxxxzz_yyyzz[k];

                g_0_x_xxxxzzz_yyzz[k] = -g_0_x_xxxxzz_yyzz[k] * ab_z + g_0_x_xxxxzz_yyzzz[k];

                g_0_x_xxxxzzz_yzzz[k] = -g_0_x_xxxxzz_yzzz[k] * ab_z + g_0_x_xxxxzz_yzzzz[k];

                g_0_x_xxxxzzz_zzzz[k] = -g_0_x_xxxxzz_zzzz[k] * ab_z + g_0_x_xxxxzz_zzzzz[k];
            }

            /// Set up 150-165 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyyy_xxxx = cbuffer.data(kg_geom_01_off + 150 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xxxy = cbuffer.data(kg_geom_01_off + 151 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xxxz = cbuffer.data(kg_geom_01_off + 152 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xxyy = cbuffer.data(kg_geom_01_off + 153 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xxyz = cbuffer.data(kg_geom_01_off + 154 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xxzz = cbuffer.data(kg_geom_01_off + 155 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xyyy = cbuffer.data(kg_geom_01_off + 156 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xyyz = cbuffer.data(kg_geom_01_off + 157 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xyzz = cbuffer.data(kg_geom_01_off + 158 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_xzzz = cbuffer.data(kg_geom_01_off + 159 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_yyyy = cbuffer.data(kg_geom_01_off + 160 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_yyyz = cbuffer.data(kg_geom_01_off + 161 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_yyzz = cbuffer.data(kg_geom_01_off + 162 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_yzzz = cbuffer.data(kg_geom_01_off + 163 * ccomps * dcomps);

            auto g_0_x_xxxyyyy_zzzz = cbuffer.data(kg_geom_01_off + 164 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyy_xxxx, g_0_x_xxxyyy_xxxxy, g_0_x_xxxyyy_xxxy, g_0_x_xxxyyy_xxxyy, g_0_x_xxxyyy_xxxyz, g_0_x_xxxyyy_xxxz, g_0_x_xxxyyy_xxyy, g_0_x_xxxyyy_xxyyy, g_0_x_xxxyyy_xxyyz, g_0_x_xxxyyy_xxyz, g_0_x_xxxyyy_xxyzz, g_0_x_xxxyyy_xxzz, g_0_x_xxxyyy_xyyy, g_0_x_xxxyyy_xyyyy, g_0_x_xxxyyy_xyyyz, g_0_x_xxxyyy_xyyz, g_0_x_xxxyyy_xyyzz, g_0_x_xxxyyy_xyzz, g_0_x_xxxyyy_xyzzz, g_0_x_xxxyyy_xzzz, g_0_x_xxxyyy_yyyy, g_0_x_xxxyyy_yyyyy, g_0_x_xxxyyy_yyyyz, g_0_x_xxxyyy_yyyz, g_0_x_xxxyyy_yyyzz, g_0_x_xxxyyy_yyzz, g_0_x_xxxyyy_yyzzz, g_0_x_xxxyyy_yzzz, g_0_x_xxxyyy_yzzzz, g_0_x_xxxyyy_zzzz, g_0_x_xxxyyyy_xxxx, g_0_x_xxxyyyy_xxxy, g_0_x_xxxyyyy_xxxz, g_0_x_xxxyyyy_xxyy, g_0_x_xxxyyyy_xxyz, g_0_x_xxxyyyy_xxzz, g_0_x_xxxyyyy_xyyy, g_0_x_xxxyyyy_xyyz, g_0_x_xxxyyyy_xyzz, g_0_x_xxxyyyy_xzzz, g_0_x_xxxyyyy_yyyy, g_0_x_xxxyyyy_yyyz, g_0_x_xxxyyyy_yyzz, g_0_x_xxxyyyy_yzzz, g_0_x_xxxyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyyy_xxxx[k] = -g_0_x_xxxyyy_xxxx[k] * ab_y + g_0_x_xxxyyy_xxxxy[k];

                g_0_x_xxxyyyy_xxxy[k] = -g_0_x_xxxyyy_xxxy[k] * ab_y + g_0_x_xxxyyy_xxxyy[k];

                g_0_x_xxxyyyy_xxxz[k] = -g_0_x_xxxyyy_xxxz[k] * ab_y + g_0_x_xxxyyy_xxxyz[k];

                g_0_x_xxxyyyy_xxyy[k] = -g_0_x_xxxyyy_xxyy[k] * ab_y + g_0_x_xxxyyy_xxyyy[k];

                g_0_x_xxxyyyy_xxyz[k] = -g_0_x_xxxyyy_xxyz[k] * ab_y + g_0_x_xxxyyy_xxyyz[k];

                g_0_x_xxxyyyy_xxzz[k] = -g_0_x_xxxyyy_xxzz[k] * ab_y + g_0_x_xxxyyy_xxyzz[k];

                g_0_x_xxxyyyy_xyyy[k] = -g_0_x_xxxyyy_xyyy[k] * ab_y + g_0_x_xxxyyy_xyyyy[k];

                g_0_x_xxxyyyy_xyyz[k] = -g_0_x_xxxyyy_xyyz[k] * ab_y + g_0_x_xxxyyy_xyyyz[k];

                g_0_x_xxxyyyy_xyzz[k] = -g_0_x_xxxyyy_xyzz[k] * ab_y + g_0_x_xxxyyy_xyyzz[k];

                g_0_x_xxxyyyy_xzzz[k] = -g_0_x_xxxyyy_xzzz[k] * ab_y + g_0_x_xxxyyy_xyzzz[k];

                g_0_x_xxxyyyy_yyyy[k] = -g_0_x_xxxyyy_yyyy[k] * ab_y + g_0_x_xxxyyy_yyyyy[k];

                g_0_x_xxxyyyy_yyyz[k] = -g_0_x_xxxyyy_yyyz[k] * ab_y + g_0_x_xxxyyy_yyyyz[k];

                g_0_x_xxxyyyy_yyzz[k] = -g_0_x_xxxyyy_yyzz[k] * ab_y + g_0_x_xxxyyy_yyyzz[k];

                g_0_x_xxxyyyy_yzzz[k] = -g_0_x_xxxyyy_yzzz[k] * ab_y + g_0_x_xxxyyy_yyzzz[k];

                g_0_x_xxxyyyy_zzzz[k] = -g_0_x_xxxyyy_zzzz[k] * ab_y + g_0_x_xxxyyy_yzzzz[k];
            }

            /// Set up 165-180 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyyz_xxxx = cbuffer.data(kg_geom_01_off + 165 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xxxy = cbuffer.data(kg_geom_01_off + 166 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xxxz = cbuffer.data(kg_geom_01_off + 167 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xxyy = cbuffer.data(kg_geom_01_off + 168 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xxyz = cbuffer.data(kg_geom_01_off + 169 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xxzz = cbuffer.data(kg_geom_01_off + 170 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xyyy = cbuffer.data(kg_geom_01_off + 171 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xyyz = cbuffer.data(kg_geom_01_off + 172 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xyzz = cbuffer.data(kg_geom_01_off + 173 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_xzzz = cbuffer.data(kg_geom_01_off + 174 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_yyyy = cbuffer.data(kg_geom_01_off + 175 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_yyyz = cbuffer.data(kg_geom_01_off + 176 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_yyzz = cbuffer.data(kg_geom_01_off + 177 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_yzzz = cbuffer.data(kg_geom_01_off + 178 * ccomps * dcomps);

            auto g_0_x_xxxyyyz_zzzz = cbuffer.data(kg_geom_01_off + 179 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyyz_xxxx, g_0_x_xxxyyyz_xxxy, g_0_x_xxxyyyz_xxxz, g_0_x_xxxyyyz_xxyy, g_0_x_xxxyyyz_xxyz, g_0_x_xxxyyyz_xxzz, g_0_x_xxxyyyz_xyyy, g_0_x_xxxyyyz_xyyz, g_0_x_xxxyyyz_xyzz, g_0_x_xxxyyyz_xzzz, g_0_x_xxxyyyz_yyyy, g_0_x_xxxyyyz_yyyz, g_0_x_xxxyyyz_yyzz, g_0_x_xxxyyyz_yzzz, g_0_x_xxxyyyz_zzzz, g_0_x_xxxyyz_xxxx, g_0_x_xxxyyz_xxxxy, g_0_x_xxxyyz_xxxy, g_0_x_xxxyyz_xxxyy, g_0_x_xxxyyz_xxxyz, g_0_x_xxxyyz_xxxz, g_0_x_xxxyyz_xxyy, g_0_x_xxxyyz_xxyyy, g_0_x_xxxyyz_xxyyz, g_0_x_xxxyyz_xxyz, g_0_x_xxxyyz_xxyzz, g_0_x_xxxyyz_xxzz, g_0_x_xxxyyz_xyyy, g_0_x_xxxyyz_xyyyy, g_0_x_xxxyyz_xyyyz, g_0_x_xxxyyz_xyyz, g_0_x_xxxyyz_xyyzz, g_0_x_xxxyyz_xyzz, g_0_x_xxxyyz_xyzzz, g_0_x_xxxyyz_xzzz, g_0_x_xxxyyz_yyyy, g_0_x_xxxyyz_yyyyy, g_0_x_xxxyyz_yyyyz, g_0_x_xxxyyz_yyyz, g_0_x_xxxyyz_yyyzz, g_0_x_xxxyyz_yyzz, g_0_x_xxxyyz_yyzzz, g_0_x_xxxyyz_yzzz, g_0_x_xxxyyz_yzzzz, g_0_x_xxxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyyz_xxxx[k] = -g_0_x_xxxyyz_xxxx[k] * ab_y + g_0_x_xxxyyz_xxxxy[k];

                g_0_x_xxxyyyz_xxxy[k] = -g_0_x_xxxyyz_xxxy[k] * ab_y + g_0_x_xxxyyz_xxxyy[k];

                g_0_x_xxxyyyz_xxxz[k] = -g_0_x_xxxyyz_xxxz[k] * ab_y + g_0_x_xxxyyz_xxxyz[k];

                g_0_x_xxxyyyz_xxyy[k] = -g_0_x_xxxyyz_xxyy[k] * ab_y + g_0_x_xxxyyz_xxyyy[k];

                g_0_x_xxxyyyz_xxyz[k] = -g_0_x_xxxyyz_xxyz[k] * ab_y + g_0_x_xxxyyz_xxyyz[k];

                g_0_x_xxxyyyz_xxzz[k] = -g_0_x_xxxyyz_xxzz[k] * ab_y + g_0_x_xxxyyz_xxyzz[k];

                g_0_x_xxxyyyz_xyyy[k] = -g_0_x_xxxyyz_xyyy[k] * ab_y + g_0_x_xxxyyz_xyyyy[k];

                g_0_x_xxxyyyz_xyyz[k] = -g_0_x_xxxyyz_xyyz[k] * ab_y + g_0_x_xxxyyz_xyyyz[k];

                g_0_x_xxxyyyz_xyzz[k] = -g_0_x_xxxyyz_xyzz[k] * ab_y + g_0_x_xxxyyz_xyyzz[k];

                g_0_x_xxxyyyz_xzzz[k] = -g_0_x_xxxyyz_xzzz[k] * ab_y + g_0_x_xxxyyz_xyzzz[k];

                g_0_x_xxxyyyz_yyyy[k] = -g_0_x_xxxyyz_yyyy[k] * ab_y + g_0_x_xxxyyz_yyyyy[k];

                g_0_x_xxxyyyz_yyyz[k] = -g_0_x_xxxyyz_yyyz[k] * ab_y + g_0_x_xxxyyz_yyyyz[k];

                g_0_x_xxxyyyz_yyzz[k] = -g_0_x_xxxyyz_yyzz[k] * ab_y + g_0_x_xxxyyz_yyyzz[k];

                g_0_x_xxxyyyz_yzzz[k] = -g_0_x_xxxyyz_yzzz[k] * ab_y + g_0_x_xxxyyz_yyzzz[k];

                g_0_x_xxxyyyz_zzzz[k] = -g_0_x_xxxyyz_zzzz[k] * ab_y + g_0_x_xxxyyz_yzzzz[k];
            }

            /// Set up 180-195 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyyzz_xxxx = cbuffer.data(kg_geom_01_off + 180 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xxxy = cbuffer.data(kg_geom_01_off + 181 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xxxz = cbuffer.data(kg_geom_01_off + 182 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xxyy = cbuffer.data(kg_geom_01_off + 183 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xxyz = cbuffer.data(kg_geom_01_off + 184 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xxzz = cbuffer.data(kg_geom_01_off + 185 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xyyy = cbuffer.data(kg_geom_01_off + 186 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xyyz = cbuffer.data(kg_geom_01_off + 187 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xyzz = cbuffer.data(kg_geom_01_off + 188 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_xzzz = cbuffer.data(kg_geom_01_off + 189 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_yyyy = cbuffer.data(kg_geom_01_off + 190 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_yyyz = cbuffer.data(kg_geom_01_off + 191 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_yyzz = cbuffer.data(kg_geom_01_off + 192 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_yzzz = cbuffer.data(kg_geom_01_off + 193 * ccomps * dcomps);

            auto g_0_x_xxxyyzz_zzzz = cbuffer.data(kg_geom_01_off + 194 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyyzz_xxxx, g_0_x_xxxyyzz_xxxy, g_0_x_xxxyyzz_xxxz, g_0_x_xxxyyzz_xxyy, g_0_x_xxxyyzz_xxyz, g_0_x_xxxyyzz_xxzz, g_0_x_xxxyyzz_xyyy, g_0_x_xxxyyzz_xyyz, g_0_x_xxxyyzz_xyzz, g_0_x_xxxyyzz_xzzz, g_0_x_xxxyyzz_yyyy, g_0_x_xxxyyzz_yyyz, g_0_x_xxxyyzz_yyzz, g_0_x_xxxyyzz_yzzz, g_0_x_xxxyyzz_zzzz, g_0_x_xxxyzz_xxxx, g_0_x_xxxyzz_xxxxy, g_0_x_xxxyzz_xxxy, g_0_x_xxxyzz_xxxyy, g_0_x_xxxyzz_xxxyz, g_0_x_xxxyzz_xxxz, g_0_x_xxxyzz_xxyy, g_0_x_xxxyzz_xxyyy, g_0_x_xxxyzz_xxyyz, g_0_x_xxxyzz_xxyz, g_0_x_xxxyzz_xxyzz, g_0_x_xxxyzz_xxzz, g_0_x_xxxyzz_xyyy, g_0_x_xxxyzz_xyyyy, g_0_x_xxxyzz_xyyyz, g_0_x_xxxyzz_xyyz, g_0_x_xxxyzz_xyyzz, g_0_x_xxxyzz_xyzz, g_0_x_xxxyzz_xyzzz, g_0_x_xxxyzz_xzzz, g_0_x_xxxyzz_yyyy, g_0_x_xxxyzz_yyyyy, g_0_x_xxxyzz_yyyyz, g_0_x_xxxyzz_yyyz, g_0_x_xxxyzz_yyyzz, g_0_x_xxxyzz_yyzz, g_0_x_xxxyzz_yyzzz, g_0_x_xxxyzz_yzzz, g_0_x_xxxyzz_yzzzz, g_0_x_xxxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyyzz_xxxx[k] = -g_0_x_xxxyzz_xxxx[k] * ab_y + g_0_x_xxxyzz_xxxxy[k];

                g_0_x_xxxyyzz_xxxy[k] = -g_0_x_xxxyzz_xxxy[k] * ab_y + g_0_x_xxxyzz_xxxyy[k];

                g_0_x_xxxyyzz_xxxz[k] = -g_0_x_xxxyzz_xxxz[k] * ab_y + g_0_x_xxxyzz_xxxyz[k];

                g_0_x_xxxyyzz_xxyy[k] = -g_0_x_xxxyzz_xxyy[k] * ab_y + g_0_x_xxxyzz_xxyyy[k];

                g_0_x_xxxyyzz_xxyz[k] = -g_0_x_xxxyzz_xxyz[k] * ab_y + g_0_x_xxxyzz_xxyyz[k];

                g_0_x_xxxyyzz_xxzz[k] = -g_0_x_xxxyzz_xxzz[k] * ab_y + g_0_x_xxxyzz_xxyzz[k];

                g_0_x_xxxyyzz_xyyy[k] = -g_0_x_xxxyzz_xyyy[k] * ab_y + g_0_x_xxxyzz_xyyyy[k];

                g_0_x_xxxyyzz_xyyz[k] = -g_0_x_xxxyzz_xyyz[k] * ab_y + g_0_x_xxxyzz_xyyyz[k];

                g_0_x_xxxyyzz_xyzz[k] = -g_0_x_xxxyzz_xyzz[k] * ab_y + g_0_x_xxxyzz_xyyzz[k];

                g_0_x_xxxyyzz_xzzz[k] = -g_0_x_xxxyzz_xzzz[k] * ab_y + g_0_x_xxxyzz_xyzzz[k];

                g_0_x_xxxyyzz_yyyy[k] = -g_0_x_xxxyzz_yyyy[k] * ab_y + g_0_x_xxxyzz_yyyyy[k];

                g_0_x_xxxyyzz_yyyz[k] = -g_0_x_xxxyzz_yyyz[k] * ab_y + g_0_x_xxxyzz_yyyyz[k];

                g_0_x_xxxyyzz_yyzz[k] = -g_0_x_xxxyzz_yyzz[k] * ab_y + g_0_x_xxxyzz_yyyzz[k];

                g_0_x_xxxyyzz_yzzz[k] = -g_0_x_xxxyzz_yzzz[k] * ab_y + g_0_x_xxxyzz_yyzzz[k];

                g_0_x_xxxyyzz_zzzz[k] = -g_0_x_xxxyzz_zzzz[k] * ab_y + g_0_x_xxxyzz_yzzzz[k];
            }

            /// Set up 195-210 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxyzzz_xxxx = cbuffer.data(kg_geom_01_off + 195 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xxxy = cbuffer.data(kg_geom_01_off + 196 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xxxz = cbuffer.data(kg_geom_01_off + 197 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xxyy = cbuffer.data(kg_geom_01_off + 198 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xxyz = cbuffer.data(kg_geom_01_off + 199 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xxzz = cbuffer.data(kg_geom_01_off + 200 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xyyy = cbuffer.data(kg_geom_01_off + 201 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xyyz = cbuffer.data(kg_geom_01_off + 202 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xyzz = cbuffer.data(kg_geom_01_off + 203 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_xzzz = cbuffer.data(kg_geom_01_off + 204 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_yyyy = cbuffer.data(kg_geom_01_off + 205 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_yyyz = cbuffer.data(kg_geom_01_off + 206 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_yyzz = cbuffer.data(kg_geom_01_off + 207 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_yzzz = cbuffer.data(kg_geom_01_off + 208 * ccomps * dcomps);

            auto g_0_x_xxxyzzz_zzzz = cbuffer.data(kg_geom_01_off + 209 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxyzzz_xxxx, g_0_x_xxxyzzz_xxxy, g_0_x_xxxyzzz_xxxz, g_0_x_xxxyzzz_xxyy, g_0_x_xxxyzzz_xxyz, g_0_x_xxxyzzz_xxzz, g_0_x_xxxyzzz_xyyy, g_0_x_xxxyzzz_xyyz, g_0_x_xxxyzzz_xyzz, g_0_x_xxxyzzz_xzzz, g_0_x_xxxyzzz_yyyy, g_0_x_xxxyzzz_yyyz, g_0_x_xxxyzzz_yyzz, g_0_x_xxxyzzz_yzzz, g_0_x_xxxyzzz_zzzz, g_0_x_xxxzzz_xxxx, g_0_x_xxxzzz_xxxxy, g_0_x_xxxzzz_xxxy, g_0_x_xxxzzz_xxxyy, g_0_x_xxxzzz_xxxyz, g_0_x_xxxzzz_xxxz, g_0_x_xxxzzz_xxyy, g_0_x_xxxzzz_xxyyy, g_0_x_xxxzzz_xxyyz, g_0_x_xxxzzz_xxyz, g_0_x_xxxzzz_xxyzz, g_0_x_xxxzzz_xxzz, g_0_x_xxxzzz_xyyy, g_0_x_xxxzzz_xyyyy, g_0_x_xxxzzz_xyyyz, g_0_x_xxxzzz_xyyz, g_0_x_xxxzzz_xyyzz, g_0_x_xxxzzz_xyzz, g_0_x_xxxzzz_xyzzz, g_0_x_xxxzzz_xzzz, g_0_x_xxxzzz_yyyy, g_0_x_xxxzzz_yyyyy, g_0_x_xxxzzz_yyyyz, g_0_x_xxxzzz_yyyz, g_0_x_xxxzzz_yyyzz, g_0_x_xxxzzz_yyzz, g_0_x_xxxzzz_yyzzz, g_0_x_xxxzzz_yzzz, g_0_x_xxxzzz_yzzzz, g_0_x_xxxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxyzzz_xxxx[k] = -g_0_x_xxxzzz_xxxx[k] * ab_y + g_0_x_xxxzzz_xxxxy[k];

                g_0_x_xxxyzzz_xxxy[k] = -g_0_x_xxxzzz_xxxy[k] * ab_y + g_0_x_xxxzzz_xxxyy[k];

                g_0_x_xxxyzzz_xxxz[k] = -g_0_x_xxxzzz_xxxz[k] * ab_y + g_0_x_xxxzzz_xxxyz[k];

                g_0_x_xxxyzzz_xxyy[k] = -g_0_x_xxxzzz_xxyy[k] * ab_y + g_0_x_xxxzzz_xxyyy[k];

                g_0_x_xxxyzzz_xxyz[k] = -g_0_x_xxxzzz_xxyz[k] * ab_y + g_0_x_xxxzzz_xxyyz[k];

                g_0_x_xxxyzzz_xxzz[k] = -g_0_x_xxxzzz_xxzz[k] * ab_y + g_0_x_xxxzzz_xxyzz[k];

                g_0_x_xxxyzzz_xyyy[k] = -g_0_x_xxxzzz_xyyy[k] * ab_y + g_0_x_xxxzzz_xyyyy[k];

                g_0_x_xxxyzzz_xyyz[k] = -g_0_x_xxxzzz_xyyz[k] * ab_y + g_0_x_xxxzzz_xyyyz[k];

                g_0_x_xxxyzzz_xyzz[k] = -g_0_x_xxxzzz_xyzz[k] * ab_y + g_0_x_xxxzzz_xyyzz[k];

                g_0_x_xxxyzzz_xzzz[k] = -g_0_x_xxxzzz_xzzz[k] * ab_y + g_0_x_xxxzzz_xyzzz[k];

                g_0_x_xxxyzzz_yyyy[k] = -g_0_x_xxxzzz_yyyy[k] * ab_y + g_0_x_xxxzzz_yyyyy[k];

                g_0_x_xxxyzzz_yyyz[k] = -g_0_x_xxxzzz_yyyz[k] * ab_y + g_0_x_xxxzzz_yyyyz[k];

                g_0_x_xxxyzzz_yyzz[k] = -g_0_x_xxxzzz_yyzz[k] * ab_y + g_0_x_xxxzzz_yyyzz[k];

                g_0_x_xxxyzzz_yzzz[k] = -g_0_x_xxxzzz_yzzz[k] * ab_y + g_0_x_xxxzzz_yyzzz[k];

                g_0_x_xxxyzzz_zzzz[k] = -g_0_x_xxxzzz_zzzz[k] * ab_y + g_0_x_xxxzzz_yzzzz[k];
            }

            /// Set up 210-225 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxxzzzz_xxxx = cbuffer.data(kg_geom_01_off + 210 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xxxy = cbuffer.data(kg_geom_01_off + 211 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xxxz = cbuffer.data(kg_geom_01_off + 212 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xxyy = cbuffer.data(kg_geom_01_off + 213 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xxyz = cbuffer.data(kg_geom_01_off + 214 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xxzz = cbuffer.data(kg_geom_01_off + 215 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xyyy = cbuffer.data(kg_geom_01_off + 216 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xyyz = cbuffer.data(kg_geom_01_off + 217 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xyzz = cbuffer.data(kg_geom_01_off + 218 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_xzzz = cbuffer.data(kg_geom_01_off + 219 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_yyyy = cbuffer.data(kg_geom_01_off + 220 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_yyyz = cbuffer.data(kg_geom_01_off + 221 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_yyzz = cbuffer.data(kg_geom_01_off + 222 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_yzzz = cbuffer.data(kg_geom_01_off + 223 * ccomps * dcomps);

            auto g_0_x_xxxzzzz_zzzz = cbuffer.data(kg_geom_01_off + 224 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxxzzz_xxxx, g_0_x_xxxzzz_xxxxz, g_0_x_xxxzzz_xxxy, g_0_x_xxxzzz_xxxyz, g_0_x_xxxzzz_xxxz, g_0_x_xxxzzz_xxxzz, g_0_x_xxxzzz_xxyy, g_0_x_xxxzzz_xxyyz, g_0_x_xxxzzz_xxyz, g_0_x_xxxzzz_xxyzz, g_0_x_xxxzzz_xxzz, g_0_x_xxxzzz_xxzzz, g_0_x_xxxzzz_xyyy, g_0_x_xxxzzz_xyyyz, g_0_x_xxxzzz_xyyz, g_0_x_xxxzzz_xyyzz, g_0_x_xxxzzz_xyzz, g_0_x_xxxzzz_xyzzz, g_0_x_xxxzzz_xzzz, g_0_x_xxxzzz_xzzzz, g_0_x_xxxzzz_yyyy, g_0_x_xxxzzz_yyyyz, g_0_x_xxxzzz_yyyz, g_0_x_xxxzzz_yyyzz, g_0_x_xxxzzz_yyzz, g_0_x_xxxzzz_yyzzz, g_0_x_xxxzzz_yzzz, g_0_x_xxxzzz_yzzzz, g_0_x_xxxzzz_zzzz, g_0_x_xxxzzz_zzzzz, g_0_x_xxxzzzz_xxxx, g_0_x_xxxzzzz_xxxy, g_0_x_xxxzzzz_xxxz, g_0_x_xxxzzzz_xxyy, g_0_x_xxxzzzz_xxyz, g_0_x_xxxzzzz_xxzz, g_0_x_xxxzzzz_xyyy, g_0_x_xxxzzzz_xyyz, g_0_x_xxxzzzz_xyzz, g_0_x_xxxzzzz_xzzz, g_0_x_xxxzzzz_yyyy, g_0_x_xxxzzzz_yyyz, g_0_x_xxxzzzz_yyzz, g_0_x_xxxzzzz_yzzz, g_0_x_xxxzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxxzzzz_xxxx[k] = -g_0_x_xxxzzz_xxxx[k] * ab_z + g_0_x_xxxzzz_xxxxz[k];

                g_0_x_xxxzzzz_xxxy[k] = -g_0_x_xxxzzz_xxxy[k] * ab_z + g_0_x_xxxzzz_xxxyz[k];

                g_0_x_xxxzzzz_xxxz[k] = -g_0_x_xxxzzz_xxxz[k] * ab_z + g_0_x_xxxzzz_xxxzz[k];

                g_0_x_xxxzzzz_xxyy[k] = -g_0_x_xxxzzz_xxyy[k] * ab_z + g_0_x_xxxzzz_xxyyz[k];

                g_0_x_xxxzzzz_xxyz[k] = -g_0_x_xxxzzz_xxyz[k] * ab_z + g_0_x_xxxzzz_xxyzz[k];

                g_0_x_xxxzzzz_xxzz[k] = -g_0_x_xxxzzz_xxzz[k] * ab_z + g_0_x_xxxzzz_xxzzz[k];

                g_0_x_xxxzzzz_xyyy[k] = -g_0_x_xxxzzz_xyyy[k] * ab_z + g_0_x_xxxzzz_xyyyz[k];

                g_0_x_xxxzzzz_xyyz[k] = -g_0_x_xxxzzz_xyyz[k] * ab_z + g_0_x_xxxzzz_xyyzz[k];

                g_0_x_xxxzzzz_xyzz[k] = -g_0_x_xxxzzz_xyzz[k] * ab_z + g_0_x_xxxzzz_xyzzz[k];

                g_0_x_xxxzzzz_xzzz[k] = -g_0_x_xxxzzz_xzzz[k] * ab_z + g_0_x_xxxzzz_xzzzz[k];

                g_0_x_xxxzzzz_yyyy[k] = -g_0_x_xxxzzz_yyyy[k] * ab_z + g_0_x_xxxzzz_yyyyz[k];

                g_0_x_xxxzzzz_yyyz[k] = -g_0_x_xxxzzz_yyyz[k] * ab_z + g_0_x_xxxzzz_yyyzz[k];

                g_0_x_xxxzzzz_yyzz[k] = -g_0_x_xxxzzz_yyzz[k] * ab_z + g_0_x_xxxzzz_yyzzz[k];

                g_0_x_xxxzzzz_yzzz[k] = -g_0_x_xxxzzz_yzzz[k] * ab_z + g_0_x_xxxzzz_yzzzz[k];

                g_0_x_xxxzzzz_zzzz[k] = -g_0_x_xxxzzz_zzzz[k] * ab_z + g_0_x_xxxzzz_zzzzz[k];
            }

            /// Set up 225-240 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyyy_xxxx = cbuffer.data(kg_geom_01_off + 225 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xxxy = cbuffer.data(kg_geom_01_off + 226 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xxxz = cbuffer.data(kg_geom_01_off + 227 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xxyy = cbuffer.data(kg_geom_01_off + 228 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xxyz = cbuffer.data(kg_geom_01_off + 229 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xxzz = cbuffer.data(kg_geom_01_off + 230 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xyyy = cbuffer.data(kg_geom_01_off + 231 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xyyz = cbuffer.data(kg_geom_01_off + 232 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xyzz = cbuffer.data(kg_geom_01_off + 233 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_xzzz = cbuffer.data(kg_geom_01_off + 234 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_yyyy = cbuffer.data(kg_geom_01_off + 235 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_yyyz = cbuffer.data(kg_geom_01_off + 236 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_yyzz = cbuffer.data(kg_geom_01_off + 237 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_yzzz = cbuffer.data(kg_geom_01_off + 238 * ccomps * dcomps);

            auto g_0_x_xxyyyyy_zzzz = cbuffer.data(kg_geom_01_off + 239 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyy_xxxx, g_0_x_xxyyyy_xxxxy, g_0_x_xxyyyy_xxxy, g_0_x_xxyyyy_xxxyy, g_0_x_xxyyyy_xxxyz, g_0_x_xxyyyy_xxxz, g_0_x_xxyyyy_xxyy, g_0_x_xxyyyy_xxyyy, g_0_x_xxyyyy_xxyyz, g_0_x_xxyyyy_xxyz, g_0_x_xxyyyy_xxyzz, g_0_x_xxyyyy_xxzz, g_0_x_xxyyyy_xyyy, g_0_x_xxyyyy_xyyyy, g_0_x_xxyyyy_xyyyz, g_0_x_xxyyyy_xyyz, g_0_x_xxyyyy_xyyzz, g_0_x_xxyyyy_xyzz, g_0_x_xxyyyy_xyzzz, g_0_x_xxyyyy_xzzz, g_0_x_xxyyyy_yyyy, g_0_x_xxyyyy_yyyyy, g_0_x_xxyyyy_yyyyz, g_0_x_xxyyyy_yyyz, g_0_x_xxyyyy_yyyzz, g_0_x_xxyyyy_yyzz, g_0_x_xxyyyy_yyzzz, g_0_x_xxyyyy_yzzz, g_0_x_xxyyyy_yzzzz, g_0_x_xxyyyy_zzzz, g_0_x_xxyyyyy_xxxx, g_0_x_xxyyyyy_xxxy, g_0_x_xxyyyyy_xxxz, g_0_x_xxyyyyy_xxyy, g_0_x_xxyyyyy_xxyz, g_0_x_xxyyyyy_xxzz, g_0_x_xxyyyyy_xyyy, g_0_x_xxyyyyy_xyyz, g_0_x_xxyyyyy_xyzz, g_0_x_xxyyyyy_xzzz, g_0_x_xxyyyyy_yyyy, g_0_x_xxyyyyy_yyyz, g_0_x_xxyyyyy_yyzz, g_0_x_xxyyyyy_yzzz, g_0_x_xxyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyyy_xxxx[k] = -g_0_x_xxyyyy_xxxx[k] * ab_y + g_0_x_xxyyyy_xxxxy[k];

                g_0_x_xxyyyyy_xxxy[k] = -g_0_x_xxyyyy_xxxy[k] * ab_y + g_0_x_xxyyyy_xxxyy[k];

                g_0_x_xxyyyyy_xxxz[k] = -g_0_x_xxyyyy_xxxz[k] * ab_y + g_0_x_xxyyyy_xxxyz[k];

                g_0_x_xxyyyyy_xxyy[k] = -g_0_x_xxyyyy_xxyy[k] * ab_y + g_0_x_xxyyyy_xxyyy[k];

                g_0_x_xxyyyyy_xxyz[k] = -g_0_x_xxyyyy_xxyz[k] * ab_y + g_0_x_xxyyyy_xxyyz[k];

                g_0_x_xxyyyyy_xxzz[k] = -g_0_x_xxyyyy_xxzz[k] * ab_y + g_0_x_xxyyyy_xxyzz[k];

                g_0_x_xxyyyyy_xyyy[k] = -g_0_x_xxyyyy_xyyy[k] * ab_y + g_0_x_xxyyyy_xyyyy[k];

                g_0_x_xxyyyyy_xyyz[k] = -g_0_x_xxyyyy_xyyz[k] * ab_y + g_0_x_xxyyyy_xyyyz[k];

                g_0_x_xxyyyyy_xyzz[k] = -g_0_x_xxyyyy_xyzz[k] * ab_y + g_0_x_xxyyyy_xyyzz[k];

                g_0_x_xxyyyyy_xzzz[k] = -g_0_x_xxyyyy_xzzz[k] * ab_y + g_0_x_xxyyyy_xyzzz[k];

                g_0_x_xxyyyyy_yyyy[k] = -g_0_x_xxyyyy_yyyy[k] * ab_y + g_0_x_xxyyyy_yyyyy[k];

                g_0_x_xxyyyyy_yyyz[k] = -g_0_x_xxyyyy_yyyz[k] * ab_y + g_0_x_xxyyyy_yyyyz[k];

                g_0_x_xxyyyyy_yyzz[k] = -g_0_x_xxyyyy_yyzz[k] * ab_y + g_0_x_xxyyyy_yyyzz[k];

                g_0_x_xxyyyyy_yzzz[k] = -g_0_x_xxyyyy_yzzz[k] * ab_y + g_0_x_xxyyyy_yyzzz[k];

                g_0_x_xxyyyyy_zzzz[k] = -g_0_x_xxyyyy_zzzz[k] * ab_y + g_0_x_xxyyyy_yzzzz[k];
            }

            /// Set up 240-255 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyyz_xxxx = cbuffer.data(kg_geom_01_off + 240 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xxxy = cbuffer.data(kg_geom_01_off + 241 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xxxz = cbuffer.data(kg_geom_01_off + 242 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xxyy = cbuffer.data(kg_geom_01_off + 243 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xxyz = cbuffer.data(kg_geom_01_off + 244 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xxzz = cbuffer.data(kg_geom_01_off + 245 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xyyy = cbuffer.data(kg_geom_01_off + 246 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xyyz = cbuffer.data(kg_geom_01_off + 247 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xyzz = cbuffer.data(kg_geom_01_off + 248 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_xzzz = cbuffer.data(kg_geom_01_off + 249 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_yyyy = cbuffer.data(kg_geom_01_off + 250 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_yyyz = cbuffer.data(kg_geom_01_off + 251 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_yyzz = cbuffer.data(kg_geom_01_off + 252 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_yzzz = cbuffer.data(kg_geom_01_off + 253 * ccomps * dcomps);

            auto g_0_x_xxyyyyz_zzzz = cbuffer.data(kg_geom_01_off + 254 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyyz_xxxx, g_0_x_xxyyyyz_xxxy, g_0_x_xxyyyyz_xxxz, g_0_x_xxyyyyz_xxyy, g_0_x_xxyyyyz_xxyz, g_0_x_xxyyyyz_xxzz, g_0_x_xxyyyyz_xyyy, g_0_x_xxyyyyz_xyyz, g_0_x_xxyyyyz_xyzz, g_0_x_xxyyyyz_xzzz, g_0_x_xxyyyyz_yyyy, g_0_x_xxyyyyz_yyyz, g_0_x_xxyyyyz_yyzz, g_0_x_xxyyyyz_yzzz, g_0_x_xxyyyyz_zzzz, g_0_x_xxyyyz_xxxx, g_0_x_xxyyyz_xxxxy, g_0_x_xxyyyz_xxxy, g_0_x_xxyyyz_xxxyy, g_0_x_xxyyyz_xxxyz, g_0_x_xxyyyz_xxxz, g_0_x_xxyyyz_xxyy, g_0_x_xxyyyz_xxyyy, g_0_x_xxyyyz_xxyyz, g_0_x_xxyyyz_xxyz, g_0_x_xxyyyz_xxyzz, g_0_x_xxyyyz_xxzz, g_0_x_xxyyyz_xyyy, g_0_x_xxyyyz_xyyyy, g_0_x_xxyyyz_xyyyz, g_0_x_xxyyyz_xyyz, g_0_x_xxyyyz_xyyzz, g_0_x_xxyyyz_xyzz, g_0_x_xxyyyz_xyzzz, g_0_x_xxyyyz_xzzz, g_0_x_xxyyyz_yyyy, g_0_x_xxyyyz_yyyyy, g_0_x_xxyyyz_yyyyz, g_0_x_xxyyyz_yyyz, g_0_x_xxyyyz_yyyzz, g_0_x_xxyyyz_yyzz, g_0_x_xxyyyz_yyzzz, g_0_x_xxyyyz_yzzz, g_0_x_xxyyyz_yzzzz, g_0_x_xxyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyyz_xxxx[k] = -g_0_x_xxyyyz_xxxx[k] * ab_y + g_0_x_xxyyyz_xxxxy[k];

                g_0_x_xxyyyyz_xxxy[k] = -g_0_x_xxyyyz_xxxy[k] * ab_y + g_0_x_xxyyyz_xxxyy[k];

                g_0_x_xxyyyyz_xxxz[k] = -g_0_x_xxyyyz_xxxz[k] * ab_y + g_0_x_xxyyyz_xxxyz[k];

                g_0_x_xxyyyyz_xxyy[k] = -g_0_x_xxyyyz_xxyy[k] * ab_y + g_0_x_xxyyyz_xxyyy[k];

                g_0_x_xxyyyyz_xxyz[k] = -g_0_x_xxyyyz_xxyz[k] * ab_y + g_0_x_xxyyyz_xxyyz[k];

                g_0_x_xxyyyyz_xxzz[k] = -g_0_x_xxyyyz_xxzz[k] * ab_y + g_0_x_xxyyyz_xxyzz[k];

                g_0_x_xxyyyyz_xyyy[k] = -g_0_x_xxyyyz_xyyy[k] * ab_y + g_0_x_xxyyyz_xyyyy[k];

                g_0_x_xxyyyyz_xyyz[k] = -g_0_x_xxyyyz_xyyz[k] * ab_y + g_0_x_xxyyyz_xyyyz[k];

                g_0_x_xxyyyyz_xyzz[k] = -g_0_x_xxyyyz_xyzz[k] * ab_y + g_0_x_xxyyyz_xyyzz[k];

                g_0_x_xxyyyyz_xzzz[k] = -g_0_x_xxyyyz_xzzz[k] * ab_y + g_0_x_xxyyyz_xyzzz[k];

                g_0_x_xxyyyyz_yyyy[k] = -g_0_x_xxyyyz_yyyy[k] * ab_y + g_0_x_xxyyyz_yyyyy[k];

                g_0_x_xxyyyyz_yyyz[k] = -g_0_x_xxyyyz_yyyz[k] * ab_y + g_0_x_xxyyyz_yyyyz[k];

                g_0_x_xxyyyyz_yyzz[k] = -g_0_x_xxyyyz_yyzz[k] * ab_y + g_0_x_xxyyyz_yyyzz[k];

                g_0_x_xxyyyyz_yzzz[k] = -g_0_x_xxyyyz_yzzz[k] * ab_y + g_0_x_xxyyyz_yyzzz[k];

                g_0_x_xxyyyyz_zzzz[k] = -g_0_x_xxyyyz_zzzz[k] * ab_y + g_0_x_xxyyyz_yzzzz[k];
            }

            /// Set up 255-270 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyyzz_xxxx = cbuffer.data(kg_geom_01_off + 255 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xxxy = cbuffer.data(kg_geom_01_off + 256 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xxxz = cbuffer.data(kg_geom_01_off + 257 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xxyy = cbuffer.data(kg_geom_01_off + 258 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xxyz = cbuffer.data(kg_geom_01_off + 259 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xxzz = cbuffer.data(kg_geom_01_off + 260 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xyyy = cbuffer.data(kg_geom_01_off + 261 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xyyz = cbuffer.data(kg_geom_01_off + 262 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xyzz = cbuffer.data(kg_geom_01_off + 263 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_xzzz = cbuffer.data(kg_geom_01_off + 264 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_yyyy = cbuffer.data(kg_geom_01_off + 265 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_yyyz = cbuffer.data(kg_geom_01_off + 266 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_yyzz = cbuffer.data(kg_geom_01_off + 267 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_yzzz = cbuffer.data(kg_geom_01_off + 268 * ccomps * dcomps);

            auto g_0_x_xxyyyzz_zzzz = cbuffer.data(kg_geom_01_off + 269 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyyzz_xxxx, g_0_x_xxyyyzz_xxxy, g_0_x_xxyyyzz_xxxz, g_0_x_xxyyyzz_xxyy, g_0_x_xxyyyzz_xxyz, g_0_x_xxyyyzz_xxzz, g_0_x_xxyyyzz_xyyy, g_0_x_xxyyyzz_xyyz, g_0_x_xxyyyzz_xyzz, g_0_x_xxyyyzz_xzzz, g_0_x_xxyyyzz_yyyy, g_0_x_xxyyyzz_yyyz, g_0_x_xxyyyzz_yyzz, g_0_x_xxyyyzz_yzzz, g_0_x_xxyyyzz_zzzz, g_0_x_xxyyzz_xxxx, g_0_x_xxyyzz_xxxxy, g_0_x_xxyyzz_xxxy, g_0_x_xxyyzz_xxxyy, g_0_x_xxyyzz_xxxyz, g_0_x_xxyyzz_xxxz, g_0_x_xxyyzz_xxyy, g_0_x_xxyyzz_xxyyy, g_0_x_xxyyzz_xxyyz, g_0_x_xxyyzz_xxyz, g_0_x_xxyyzz_xxyzz, g_0_x_xxyyzz_xxzz, g_0_x_xxyyzz_xyyy, g_0_x_xxyyzz_xyyyy, g_0_x_xxyyzz_xyyyz, g_0_x_xxyyzz_xyyz, g_0_x_xxyyzz_xyyzz, g_0_x_xxyyzz_xyzz, g_0_x_xxyyzz_xyzzz, g_0_x_xxyyzz_xzzz, g_0_x_xxyyzz_yyyy, g_0_x_xxyyzz_yyyyy, g_0_x_xxyyzz_yyyyz, g_0_x_xxyyzz_yyyz, g_0_x_xxyyzz_yyyzz, g_0_x_xxyyzz_yyzz, g_0_x_xxyyzz_yyzzz, g_0_x_xxyyzz_yzzz, g_0_x_xxyyzz_yzzzz, g_0_x_xxyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyyzz_xxxx[k] = -g_0_x_xxyyzz_xxxx[k] * ab_y + g_0_x_xxyyzz_xxxxy[k];

                g_0_x_xxyyyzz_xxxy[k] = -g_0_x_xxyyzz_xxxy[k] * ab_y + g_0_x_xxyyzz_xxxyy[k];

                g_0_x_xxyyyzz_xxxz[k] = -g_0_x_xxyyzz_xxxz[k] * ab_y + g_0_x_xxyyzz_xxxyz[k];

                g_0_x_xxyyyzz_xxyy[k] = -g_0_x_xxyyzz_xxyy[k] * ab_y + g_0_x_xxyyzz_xxyyy[k];

                g_0_x_xxyyyzz_xxyz[k] = -g_0_x_xxyyzz_xxyz[k] * ab_y + g_0_x_xxyyzz_xxyyz[k];

                g_0_x_xxyyyzz_xxzz[k] = -g_0_x_xxyyzz_xxzz[k] * ab_y + g_0_x_xxyyzz_xxyzz[k];

                g_0_x_xxyyyzz_xyyy[k] = -g_0_x_xxyyzz_xyyy[k] * ab_y + g_0_x_xxyyzz_xyyyy[k];

                g_0_x_xxyyyzz_xyyz[k] = -g_0_x_xxyyzz_xyyz[k] * ab_y + g_0_x_xxyyzz_xyyyz[k];

                g_0_x_xxyyyzz_xyzz[k] = -g_0_x_xxyyzz_xyzz[k] * ab_y + g_0_x_xxyyzz_xyyzz[k];

                g_0_x_xxyyyzz_xzzz[k] = -g_0_x_xxyyzz_xzzz[k] * ab_y + g_0_x_xxyyzz_xyzzz[k];

                g_0_x_xxyyyzz_yyyy[k] = -g_0_x_xxyyzz_yyyy[k] * ab_y + g_0_x_xxyyzz_yyyyy[k];

                g_0_x_xxyyyzz_yyyz[k] = -g_0_x_xxyyzz_yyyz[k] * ab_y + g_0_x_xxyyzz_yyyyz[k];

                g_0_x_xxyyyzz_yyzz[k] = -g_0_x_xxyyzz_yyzz[k] * ab_y + g_0_x_xxyyzz_yyyzz[k];

                g_0_x_xxyyyzz_yzzz[k] = -g_0_x_xxyyzz_yzzz[k] * ab_y + g_0_x_xxyyzz_yyzzz[k];

                g_0_x_xxyyyzz_zzzz[k] = -g_0_x_xxyyzz_zzzz[k] * ab_y + g_0_x_xxyyzz_yzzzz[k];
            }

            /// Set up 270-285 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyyzzz_xxxx = cbuffer.data(kg_geom_01_off + 270 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xxxy = cbuffer.data(kg_geom_01_off + 271 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xxxz = cbuffer.data(kg_geom_01_off + 272 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xxyy = cbuffer.data(kg_geom_01_off + 273 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xxyz = cbuffer.data(kg_geom_01_off + 274 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xxzz = cbuffer.data(kg_geom_01_off + 275 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xyyy = cbuffer.data(kg_geom_01_off + 276 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xyyz = cbuffer.data(kg_geom_01_off + 277 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xyzz = cbuffer.data(kg_geom_01_off + 278 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_xzzz = cbuffer.data(kg_geom_01_off + 279 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_yyyy = cbuffer.data(kg_geom_01_off + 280 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_yyyz = cbuffer.data(kg_geom_01_off + 281 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_yyzz = cbuffer.data(kg_geom_01_off + 282 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_yzzz = cbuffer.data(kg_geom_01_off + 283 * ccomps * dcomps);

            auto g_0_x_xxyyzzz_zzzz = cbuffer.data(kg_geom_01_off + 284 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyyzzz_xxxx, g_0_x_xxyyzzz_xxxy, g_0_x_xxyyzzz_xxxz, g_0_x_xxyyzzz_xxyy, g_0_x_xxyyzzz_xxyz, g_0_x_xxyyzzz_xxzz, g_0_x_xxyyzzz_xyyy, g_0_x_xxyyzzz_xyyz, g_0_x_xxyyzzz_xyzz, g_0_x_xxyyzzz_xzzz, g_0_x_xxyyzzz_yyyy, g_0_x_xxyyzzz_yyyz, g_0_x_xxyyzzz_yyzz, g_0_x_xxyyzzz_yzzz, g_0_x_xxyyzzz_zzzz, g_0_x_xxyzzz_xxxx, g_0_x_xxyzzz_xxxxy, g_0_x_xxyzzz_xxxy, g_0_x_xxyzzz_xxxyy, g_0_x_xxyzzz_xxxyz, g_0_x_xxyzzz_xxxz, g_0_x_xxyzzz_xxyy, g_0_x_xxyzzz_xxyyy, g_0_x_xxyzzz_xxyyz, g_0_x_xxyzzz_xxyz, g_0_x_xxyzzz_xxyzz, g_0_x_xxyzzz_xxzz, g_0_x_xxyzzz_xyyy, g_0_x_xxyzzz_xyyyy, g_0_x_xxyzzz_xyyyz, g_0_x_xxyzzz_xyyz, g_0_x_xxyzzz_xyyzz, g_0_x_xxyzzz_xyzz, g_0_x_xxyzzz_xyzzz, g_0_x_xxyzzz_xzzz, g_0_x_xxyzzz_yyyy, g_0_x_xxyzzz_yyyyy, g_0_x_xxyzzz_yyyyz, g_0_x_xxyzzz_yyyz, g_0_x_xxyzzz_yyyzz, g_0_x_xxyzzz_yyzz, g_0_x_xxyzzz_yyzzz, g_0_x_xxyzzz_yzzz, g_0_x_xxyzzz_yzzzz, g_0_x_xxyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyyzzz_xxxx[k] = -g_0_x_xxyzzz_xxxx[k] * ab_y + g_0_x_xxyzzz_xxxxy[k];

                g_0_x_xxyyzzz_xxxy[k] = -g_0_x_xxyzzz_xxxy[k] * ab_y + g_0_x_xxyzzz_xxxyy[k];

                g_0_x_xxyyzzz_xxxz[k] = -g_0_x_xxyzzz_xxxz[k] * ab_y + g_0_x_xxyzzz_xxxyz[k];

                g_0_x_xxyyzzz_xxyy[k] = -g_0_x_xxyzzz_xxyy[k] * ab_y + g_0_x_xxyzzz_xxyyy[k];

                g_0_x_xxyyzzz_xxyz[k] = -g_0_x_xxyzzz_xxyz[k] * ab_y + g_0_x_xxyzzz_xxyyz[k];

                g_0_x_xxyyzzz_xxzz[k] = -g_0_x_xxyzzz_xxzz[k] * ab_y + g_0_x_xxyzzz_xxyzz[k];

                g_0_x_xxyyzzz_xyyy[k] = -g_0_x_xxyzzz_xyyy[k] * ab_y + g_0_x_xxyzzz_xyyyy[k];

                g_0_x_xxyyzzz_xyyz[k] = -g_0_x_xxyzzz_xyyz[k] * ab_y + g_0_x_xxyzzz_xyyyz[k];

                g_0_x_xxyyzzz_xyzz[k] = -g_0_x_xxyzzz_xyzz[k] * ab_y + g_0_x_xxyzzz_xyyzz[k];

                g_0_x_xxyyzzz_xzzz[k] = -g_0_x_xxyzzz_xzzz[k] * ab_y + g_0_x_xxyzzz_xyzzz[k];

                g_0_x_xxyyzzz_yyyy[k] = -g_0_x_xxyzzz_yyyy[k] * ab_y + g_0_x_xxyzzz_yyyyy[k];

                g_0_x_xxyyzzz_yyyz[k] = -g_0_x_xxyzzz_yyyz[k] * ab_y + g_0_x_xxyzzz_yyyyz[k];

                g_0_x_xxyyzzz_yyzz[k] = -g_0_x_xxyzzz_yyzz[k] * ab_y + g_0_x_xxyzzz_yyyzz[k];

                g_0_x_xxyyzzz_yzzz[k] = -g_0_x_xxyzzz_yzzz[k] * ab_y + g_0_x_xxyzzz_yyzzz[k];

                g_0_x_xxyyzzz_zzzz[k] = -g_0_x_xxyzzz_zzzz[k] * ab_y + g_0_x_xxyzzz_yzzzz[k];
            }

            /// Set up 285-300 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxyzzzz_xxxx = cbuffer.data(kg_geom_01_off + 285 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xxxy = cbuffer.data(kg_geom_01_off + 286 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xxxz = cbuffer.data(kg_geom_01_off + 287 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xxyy = cbuffer.data(kg_geom_01_off + 288 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xxyz = cbuffer.data(kg_geom_01_off + 289 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xxzz = cbuffer.data(kg_geom_01_off + 290 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xyyy = cbuffer.data(kg_geom_01_off + 291 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xyyz = cbuffer.data(kg_geom_01_off + 292 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xyzz = cbuffer.data(kg_geom_01_off + 293 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_xzzz = cbuffer.data(kg_geom_01_off + 294 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_yyyy = cbuffer.data(kg_geom_01_off + 295 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_yyyz = cbuffer.data(kg_geom_01_off + 296 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_yyzz = cbuffer.data(kg_geom_01_off + 297 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_yzzz = cbuffer.data(kg_geom_01_off + 298 * ccomps * dcomps);

            auto g_0_x_xxyzzzz_zzzz = cbuffer.data(kg_geom_01_off + 299 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxyzzzz_xxxx, g_0_x_xxyzzzz_xxxy, g_0_x_xxyzzzz_xxxz, g_0_x_xxyzzzz_xxyy, g_0_x_xxyzzzz_xxyz, g_0_x_xxyzzzz_xxzz, g_0_x_xxyzzzz_xyyy, g_0_x_xxyzzzz_xyyz, g_0_x_xxyzzzz_xyzz, g_0_x_xxyzzzz_xzzz, g_0_x_xxyzzzz_yyyy, g_0_x_xxyzzzz_yyyz, g_0_x_xxyzzzz_yyzz, g_0_x_xxyzzzz_yzzz, g_0_x_xxyzzzz_zzzz, g_0_x_xxzzzz_xxxx, g_0_x_xxzzzz_xxxxy, g_0_x_xxzzzz_xxxy, g_0_x_xxzzzz_xxxyy, g_0_x_xxzzzz_xxxyz, g_0_x_xxzzzz_xxxz, g_0_x_xxzzzz_xxyy, g_0_x_xxzzzz_xxyyy, g_0_x_xxzzzz_xxyyz, g_0_x_xxzzzz_xxyz, g_0_x_xxzzzz_xxyzz, g_0_x_xxzzzz_xxzz, g_0_x_xxzzzz_xyyy, g_0_x_xxzzzz_xyyyy, g_0_x_xxzzzz_xyyyz, g_0_x_xxzzzz_xyyz, g_0_x_xxzzzz_xyyzz, g_0_x_xxzzzz_xyzz, g_0_x_xxzzzz_xyzzz, g_0_x_xxzzzz_xzzz, g_0_x_xxzzzz_yyyy, g_0_x_xxzzzz_yyyyy, g_0_x_xxzzzz_yyyyz, g_0_x_xxzzzz_yyyz, g_0_x_xxzzzz_yyyzz, g_0_x_xxzzzz_yyzz, g_0_x_xxzzzz_yyzzz, g_0_x_xxzzzz_yzzz, g_0_x_xxzzzz_yzzzz, g_0_x_xxzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxyzzzz_xxxx[k] = -g_0_x_xxzzzz_xxxx[k] * ab_y + g_0_x_xxzzzz_xxxxy[k];

                g_0_x_xxyzzzz_xxxy[k] = -g_0_x_xxzzzz_xxxy[k] * ab_y + g_0_x_xxzzzz_xxxyy[k];

                g_0_x_xxyzzzz_xxxz[k] = -g_0_x_xxzzzz_xxxz[k] * ab_y + g_0_x_xxzzzz_xxxyz[k];

                g_0_x_xxyzzzz_xxyy[k] = -g_0_x_xxzzzz_xxyy[k] * ab_y + g_0_x_xxzzzz_xxyyy[k];

                g_0_x_xxyzzzz_xxyz[k] = -g_0_x_xxzzzz_xxyz[k] * ab_y + g_0_x_xxzzzz_xxyyz[k];

                g_0_x_xxyzzzz_xxzz[k] = -g_0_x_xxzzzz_xxzz[k] * ab_y + g_0_x_xxzzzz_xxyzz[k];

                g_0_x_xxyzzzz_xyyy[k] = -g_0_x_xxzzzz_xyyy[k] * ab_y + g_0_x_xxzzzz_xyyyy[k];

                g_0_x_xxyzzzz_xyyz[k] = -g_0_x_xxzzzz_xyyz[k] * ab_y + g_0_x_xxzzzz_xyyyz[k];

                g_0_x_xxyzzzz_xyzz[k] = -g_0_x_xxzzzz_xyzz[k] * ab_y + g_0_x_xxzzzz_xyyzz[k];

                g_0_x_xxyzzzz_xzzz[k] = -g_0_x_xxzzzz_xzzz[k] * ab_y + g_0_x_xxzzzz_xyzzz[k];

                g_0_x_xxyzzzz_yyyy[k] = -g_0_x_xxzzzz_yyyy[k] * ab_y + g_0_x_xxzzzz_yyyyy[k];

                g_0_x_xxyzzzz_yyyz[k] = -g_0_x_xxzzzz_yyyz[k] * ab_y + g_0_x_xxzzzz_yyyyz[k];

                g_0_x_xxyzzzz_yyzz[k] = -g_0_x_xxzzzz_yyzz[k] * ab_y + g_0_x_xxzzzz_yyyzz[k];

                g_0_x_xxyzzzz_yzzz[k] = -g_0_x_xxzzzz_yzzz[k] * ab_y + g_0_x_xxzzzz_yyzzz[k];

                g_0_x_xxyzzzz_zzzz[k] = -g_0_x_xxzzzz_zzzz[k] * ab_y + g_0_x_xxzzzz_yzzzz[k];
            }

            /// Set up 300-315 components of targeted buffer : cbuffer.data(

            auto g_0_x_xxzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 300 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 301 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 302 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 303 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 304 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 305 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 306 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 307 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 308 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 309 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 310 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 311 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 312 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 313 * ccomps * dcomps);

            auto g_0_x_xxzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 314 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xxzzzz_xxxx, g_0_x_xxzzzz_xxxxz, g_0_x_xxzzzz_xxxy, g_0_x_xxzzzz_xxxyz, g_0_x_xxzzzz_xxxz, g_0_x_xxzzzz_xxxzz, g_0_x_xxzzzz_xxyy, g_0_x_xxzzzz_xxyyz, g_0_x_xxzzzz_xxyz, g_0_x_xxzzzz_xxyzz, g_0_x_xxzzzz_xxzz, g_0_x_xxzzzz_xxzzz, g_0_x_xxzzzz_xyyy, g_0_x_xxzzzz_xyyyz, g_0_x_xxzzzz_xyyz, g_0_x_xxzzzz_xyyzz, g_0_x_xxzzzz_xyzz, g_0_x_xxzzzz_xyzzz, g_0_x_xxzzzz_xzzz, g_0_x_xxzzzz_xzzzz, g_0_x_xxzzzz_yyyy, g_0_x_xxzzzz_yyyyz, g_0_x_xxzzzz_yyyz, g_0_x_xxzzzz_yyyzz, g_0_x_xxzzzz_yyzz, g_0_x_xxzzzz_yyzzz, g_0_x_xxzzzz_yzzz, g_0_x_xxzzzz_yzzzz, g_0_x_xxzzzz_zzzz, g_0_x_xxzzzz_zzzzz, g_0_x_xxzzzzz_xxxx, g_0_x_xxzzzzz_xxxy, g_0_x_xxzzzzz_xxxz, g_0_x_xxzzzzz_xxyy, g_0_x_xxzzzzz_xxyz, g_0_x_xxzzzzz_xxzz, g_0_x_xxzzzzz_xyyy, g_0_x_xxzzzzz_xyyz, g_0_x_xxzzzzz_xyzz, g_0_x_xxzzzzz_xzzz, g_0_x_xxzzzzz_yyyy, g_0_x_xxzzzzz_yyyz, g_0_x_xxzzzzz_yyzz, g_0_x_xxzzzzz_yzzz, g_0_x_xxzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xxzzzzz_xxxx[k] = -g_0_x_xxzzzz_xxxx[k] * ab_z + g_0_x_xxzzzz_xxxxz[k];

                g_0_x_xxzzzzz_xxxy[k] = -g_0_x_xxzzzz_xxxy[k] * ab_z + g_0_x_xxzzzz_xxxyz[k];

                g_0_x_xxzzzzz_xxxz[k] = -g_0_x_xxzzzz_xxxz[k] * ab_z + g_0_x_xxzzzz_xxxzz[k];

                g_0_x_xxzzzzz_xxyy[k] = -g_0_x_xxzzzz_xxyy[k] * ab_z + g_0_x_xxzzzz_xxyyz[k];

                g_0_x_xxzzzzz_xxyz[k] = -g_0_x_xxzzzz_xxyz[k] * ab_z + g_0_x_xxzzzz_xxyzz[k];

                g_0_x_xxzzzzz_xxzz[k] = -g_0_x_xxzzzz_xxzz[k] * ab_z + g_0_x_xxzzzz_xxzzz[k];

                g_0_x_xxzzzzz_xyyy[k] = -g_0_x_xxzzzz_xyyy[k] * ab_z + g_0_x_xxzzzz_xyyyz[k];

                g_0_x_xxzzzzz_xyyz[k] = -g_0_x_xxzzzz_xyyz[k] * ab_z + g_0_x_xxzzzz_xyyzz[k];

                g_0_x_xxzzzzz_xyzz[k] = -g_0_x_xxzzzz_xyzz[k] * ab_z + g_0_x_xxzzzz_xyzzz[k];

                g_0_x_xxzzzzz_xzzz[k] = -g_0_x_xxzzzz_xzzz[k] * ab_z + g_0_x_xxzzzz_xzzzz[k];

                g_0_x_xxzzzzz_yyyy[k] = -g_0_x_xxzzzz_yyyy[k] * ab_z + g_0_x_xxzzzz_yyyyz[k];

                g_0_x_xxzzzzz_yyyz[k] = -g_0_x_xxzzzz_yyyz[k] * ab_z + g_0_x_xxzzzz_yyyzz[k];

                g_0_x_xxzzzzz_yyzz[k] = -g_0_x_xxzzzz_yyzz[k] * ab_z + g_0_x_xxzzzz_yyzzz[k];

                g_0_x_xxzzzzz_yzzz[k] = -g_0_x_xxzzzz_yzzz[k] * ab_z + g_0_x_xxzzzz_yzzzz[k];

                g_0_x_xxzzzzz_zzzz[k] = -g_0_x_xxzzzz_zzzz[k] * ab_z + g_0_x_xxzzzz_zzzzz[k];
            }

            /// Set up 315-330 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyyy_xxxx = cbuffer.data(kg_geom_01_off + 315 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xxxy = cbuffer.data(kg_geom_01_off + 316 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xxxz = cbuffer.data(kg_geom_01_off + 317 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xxyy = cbuffer.data(kg_geom_01_off + 318 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xxyz = cbuffer.data(kg_geom_01_off + 319 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xxzz = cbuffer.data(kg_geom_01_off + 320 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xyyy = cbuffer.data(kg_geom_01_off + 321 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xyyz = cbuffer.data(kg_geom_01_off + 322 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xyzz = cbuffer.data(kg_geom_01_off + 323 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_xzzz = cbuffer.data(kg_geom_01_off + 324 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_yyyy = cbuffer.data(kg_geom_01_off + 325 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_yyyz = cbuffer.data(kg_geom_01_off + 326 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_yyzz = cbuffer.data(kg_geom_01_off + 327 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_yzzz = cbuffer.data(kg_geom_01_off + 328 * ccomps * dcomps);

            auto g_0_x_xyyyyyy_zzzz = cbuffer.data(kg_geom_01_off + 329 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyy_xxxx, g_0_x_xyyyyy_xxxxy, g_0_x_xyyyyy_xxxy, g_0_x_xyyyyy_xxxyy, g_0_x_xyyyyy_xxxyz, g_0_x_xyyyyy_xxxz, g_0_x_xyyyyy_xxyy, g_0_x_xyyyyy_xxyyy, g_0_x_xyyyyy_xxyyz, g_0_x_xyyyyy_xxyz, g_0_x_xyyyyy_xxyzz, g_0_x_xyyyyy_xxzz, g_0_x_xyyyyy_xyyy, g_0_x_xyyyyy_xyyyy, g_0_x_xyyyyy_xyyyz, g_0_x_xyyyyy_xyyz, g_0_x_xyyyyy_xyyzz, g_0_x_xyyyyy_xyzz, g_0_x_xyyyyy_xyzzz, g_0_x_xyyyyy_xzzz, g_0_x_xyyyyy_yyyy, g_0_x_xyyyyy_yyyyy, g_0_x_xyyyyy_yyyyz, g_0_x_xyyyyy_yyyz, g_0_x_xyyyyy_yyyzz, g_0_x_xyyyyy_yyzz, g_0_x_xyyyyy_yyzzz, g_0_x_xyyyyy_yzzz, g_0_x_xyyyyy_yzzzz, g_0_x_xyyyyy_zzzz, g_0_x_xyyyyyy_xxxx, g_0_x_xyyyyyy_xxxy, g_0_x_xyyyyyy_xxxz, g_0_x_xyyyyyy_xxyy, g_0_x_xyyyyyy_xxyz, g_0_x_xyyyyyy_xxzz, g_0_x_xyyyyyy_xyyy, g_0_x_xyyyyyy_xyyz, g_0_x_xyyyyyy_xyzz, g_0_x_xyyyyyy_xzzz, g_0_x_xyyyyyy_yyyy, g_0_x_xyyyyyy_yyyz, g_0_x_xyyyyyy_yyzz, g_0_x_xyyyyyy_yzzz, g_0_x_xyyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyyy_xxxx[k] = -g_0_x_xyyyyy_xxxx[k] * ab_y + g_0_x_xyyyyy_xxxxy[k];

                g_0_x_xyyyyyy_xxxy[k] = -g_0_x_xyyyyy_xxxy[k] * ab_y + g_0_x_xyyyyy_xxxyy[k];

                g_0_x_xyyyyyy_xxxz[k] = -g_0_x_xyyyyy_xxxz[k] * ab_y + g_0_x_xyyyyy_xxxyz[k];

                g_0_x_xyyyyyy_xxyy[k] = -g_0_x_xyyyyy_xxyy[k] * ab_y + g_0_x_xyyyyy_xxyyy[k];

                g_0_x_xyyyyyy_xxyz[k] = -g_0_x_xyyyyy_xxyz[k] * ab_y + g_0_x_xyyyyy_xxyyz[k];

                g_0_x_xyyyyyy_xxzz[k] = -g_0_x_xyyyyy_xxzz[k] * ab_y + g_0_x_xyyyyy_xxyzz[k];

                g_0_x_xyyyyyy_xyyy[k] = -g_0_x_xyyyyy_xyyy[k] * ab_y + g_0_x_xyyyyy_xyyyy[k];

                g_0_x_xyyyyyy_xyyz[k] = -g_0_x_xyyyyy_xyyz[k] * ab_y + g_0_x_xyyyyy_xyyyz[k];

                g_0_x_xyyyyyy_xyzz[k] = -g_0_x_xyyyyy_xyzz[k] * ab_y + g_0_x_xyyyyy_xyyzz[k];

                g_0_x_xyyyyyy_xzzz[k] = -g_0_x_xyyyyy_xzzz[k] * ab_y + g_0_x_xyyyyy_xyzzz[k];

                g_0_x_xyyyyyy_yyyy[k] = -g_0_x_xyyyyy_yyyy[k] * ab_y + g_0_x_xyyyyy_yyyyy[k];

                g_0_x_xyyyyyy_yyyz[k] = -g_0_x_xyyyyy_yyyz[k] * ab_y + g_0_x_xyyyyy_yyyyz[k];

                g_0_x_xyyyyyy_yyzz[k] = -g_0_x_xyyyyy_yyzz[k] * ab_y + g_0_x_xyyyyy_yyyzz[k];

                g_0_x_xyyyyyy_yzzz[k] = -g_0_x_xyyyyy_yzzz[k] * ab_y + g_0_x_xyyyyy_yyzzz[k];

                g_0_x_xyyyyyy_zzzz[k] = -g_0_x_xyyyyy_zzzz[k] * ab_y + g_0_x_xyyyyy_yzzzz[k];
            }

            /// Set up 330-345 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyyz_xxxx = cbuffer.data(kg_geom_01_off + 330 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xxxy = cbuffer.data(kg_geom_01_off + 331 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xxxz = cbuffer.data(kg_geom_01_off + 332 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xxyy = cbuffer.data(kg_geom_01_off + 333 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xxyz = cbuffer.data(kg_geom_01_off + 334 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xxzz = cbuffer.data(kg_geom_01_off + 335 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xyyy = cbuffer.data(kg_geom_01_off + 336 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xyyz = cbuffer.data(kg_geom_01_off + 337 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xyzz = cbuffer.data(kg_geom_01_off + 338 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_xzzz = cbuffer.data(kg_geom_01_off + 339 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_yyyy = cbuffer.data(kg_geom_01_off + 340 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_yyyz = cbuffer.data(kg_geom_01_off + 341 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_yyzz = cbuffer.data(kg_geom_01_off + 342 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_yzzz = cbuffer.data(kg_geom_01_off + 343 * ccomps * dcomps);

            auto g_0_x_xyyyyyz_zzzz = cbuffer.data(kg_geom_01_off + 344 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyyz_xxxx, g_0_x_xyyyyyz_xxxy, g_0_x_xyyyyyz_xxxz, g_0_x_xyyyyyz_xxyy, g_0_x_xyyyyyz_xxyz, g_0_x_xyyyyyz_xxzz, g_0_x_xyyyyyz_xyyy, g_0_x_xyyyyyz_xyyz, g_0_x_xyyyyyz_xyzz, g_0_x_xyyyyyz_xzzz, g_0_x_xyyyyyz_yyyy, g_0_x_xyyyyyz_yyyz, g_0_x_xyyyyyz_yyzz, g_0_x_xyyyyyz_yzzz, g_0_x_xyyyyyz_zzzz, g_0_x_xyyyyz_xxxx, g_0_x_xyyyyz_xxxxy, g_0_x_xyyyyz_xxxy, g_0_x_xyyyyz_xxxyy, g_0_x_xyyyyz_xxxyz, g_0_x_xyyyyz_xxxz, g_0_x_xyyyyz_xxyy, g_0_x_xyyyyz_xxyyy, g_0_x_xyyyyz_xxyyz, g_0_x_xyyyyz_xxyz, g_0_x_xyyyyz_xxyzz, g_0_x_xyyyyz_xxzz, g_0_x_xyyyyz_xyyy, g_0_x_xyyyyz_xyyyy, g_0_x_xyyyyz_xyyyz, g_0_x_xyyyyz_xyyz, g_0_x_xyyyyz_xyyzz, g_0_x_xyyyyz_xyzz, g_0_x_xyyyyz_xyzzz, g_0_x_xyyyyz_xzzz, g_0_x_xyyyyz_yyyy, g_0_x_xyyyyz_yyyyy, g_0_x_xyyyyz_yyyyz, g_0_x_xyyyyz_yyyz, g_0_x_xyyyyz_yyyzz, g_0_x_xyyyyz_yyzz, g_0_x_xyyyyz_yyzzz, g_0_x_xyyyyz_yzzz, g_0_x_xyyyyz_yzzzz, g_0_x_xyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyyz_xxxx[k] = -g_0_x_xyyyyz_xxxx[k] * ab_y + g_0_x_xyyyyz_xxxxy[k];

                g_0_x_xyyyyyz_xxxy[k] = -g_0_x_xyyyyz_xxxy[k] * ab_y + g_0_x_xyyyyz_xxxyy[k];

                g_0_x_xyyyyyz_xxxz[k] = -g_0_x_xyyyyz_xxxz[k] * ab_y + g_0_x_xyyyyz_xxxyz[k];

                g_0_x_xyyyyyz_xxyy[k] = -g_0_x_xyyyyz_xxyy[k] * ab_y + g_0_x_xyyyyz_xxyyy[k];

                g_0_x_xyyyyyz_xxyz[k] = -g_0_x_xyyyyz_xxyz[k] * ab_y + g_0_x_xyyyyz_xxyyz[k];

                g_0_x_xyyyyyz_xxzz[k] = -g_0_x_xyyyyz_xxzz[k] * ab_y + g_0_x_xyyyyz_xxyzz[k];

                g_0_x_xyyyyyz_xyyy[k] = -g_0_x_xyyyyz_xyyy[k] * ab_y + g_0_x_xyyyyz_xyyyy[k];

                g_0_x_xyyyyyz_xyyz[k] = -g_0_x_xyyyyz_xyyz[k] * ab_y + g_0_x_xyyyyz_xyyyz[k];

                g_0_x_xyyyyyz_xyzz[k] = -g_0_x_xyyyyz_xyzz[k] * ab_y + g_0_x_xyyyyz_xyyzz[k];

                g_0_x_xyyyyyz_xzzz[k] = -g_0_x_xyyyyz_xzzz[k] * ab_y + g_0_x_xyyyyz_xyzzz[k];

                g_0_x_xyyyyyz_yyyy[k] = -g_0_x_xyyyyz_yyyy[k] * ab_y + g_0_x_xyyyyz_yyyyy[k];

                g_0_x_xyyyyyz_yyyz[k] = -g_0_x_xyyyyz_yyyz[k] * ab_y + g_0_x_xyyyyz_yyyyz[k];

                g_0_x_xyyyyyz_yyzz[k] = -g_0_x_xyyyyz_yyzz[k] * ab_y + g_0_x_xyyyyz_yyyzz[k];

                g_0_x_xyyyyyz_yzzz[k] = -g_0_x_xyyyyz_yzzz[k] * ab_y + g_0_x_xyyyyz_yyzzz[k];

                g_0_x_xyyyyyz_zzzz[k] = -g_0_x_xyyyyz_zzzz[k] * ab_y + g_0_x_xyyyyz_yzzzz[k];
            }

            /// Set up 345-360 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyyzz_xxxx = cbuffer.data(kg_geom_01_off + 345 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xxxy = cbuffer.data(kg_geom_01_off + 346 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xxxz = cbuffer.data(kg_geom_01_off + 347 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xxyy = cbuffer.data(kg_geom_01_off + 348 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xxyz = cbuffer.data(kg_geom_01_off + 349 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xxzz = cbuffer.data(kg_geom_01_off + 350 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xyyy = cbuffer.data(kg_geom_01_off + 351 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xyyz = cbuffer.data(kg_geom_01_off + 352 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xyzz = cbuffer.data(kg_geom_01_off + 353 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_xzzz = cbuffer.data(kg_geom_01_off + 354 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_yyyy = cbuffer.data(kg_geom_01_off + 355 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_yyyz = cbuffer.data(kg_geom_01_off + 356 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_yyzz = cbuffer.data(kg_geom_01_off + 357 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_yzzz = cbuffer.data(kg_geom_01_off + 358 * ccomps * dcomps);

            auto g_0_x_xyyyyzz_zzzz = cbuffer.data(kg_geom_01_off + 359 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyyzz_xxxx, g_0_x_xyyyyzz_xxxy, g_0_x_xyyyyzz_xxxz, g_0_x_xyyyyzz_xxyy, g_0_x_xyyyyzz_xxyz, g_0_x_xyyyyzz_xxzz, g_0_x_xyyyyzz_xyyy, g_0_x_xyyyyzz_xyyz, g_0_x_xyyyyzz_xyzz, g_0_x_xyyyyzz_xzzz, g_0_x_xyyyyzz_yyyy, g_0_x_xyyyyzz_yyyz, g_0_x_xyyyyzz_yyzz, g_0_x_xyyyyzz_yzzz, g_0_x_xyyyyzz_zzzz, g_0_x_xyyyzz_xxxx, g_0_x_xyyyzz_xxxxy, g_0_x_xyyyzz_xxxy, g_0_x_xyyyzz_xxxyy, g_0_x_xyyyzz_xxxyz, g_0_x_xyyyzz_xxxz, g_0_x_xyyyzz_xxyy, g_0_x_xyyyzz_xxyyy, g_0_x_xyyyzz_xxyyz, g_0_x_xyyyzz_xxyz, g_0_x_xyyyzz_xxyzz, g_0_x_xyyyzz_xxzz, g_0_x_xyyyzz_xyyy, g_0_x_xyyyzz_xyyyy, g_0_x_xyyyzz_xyyyz, g_0_x_xyyyzz_xyyz, g_0_x_xyyyzz_xyyzz, g_0_x_xyyyzz_xyzz, g_0_x_xyyyzz_xyzzz, g_0_x_xyyyzz_xzzz, g_0_x_xyyyzz_yyyy, g_0_x_xyyyzz_yyyyy, g_0_x_xyyyzz_yyyyz, g_0_x_xyyyzz_yyyz, g_0_x_xyyyzz_yyyzz, g_0_x_xyyyzz_yyzz, g_0_x_xyyyzz_yyzzz, g_0_x_xyyyzz_yzzz, g_0_x_xyyyzz_yzzzz, g_0_x_xyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyyzz_xxxx[k] = -g_0_x_xyyyzz_xxxx[k] * ab_y + g_0_x_xyyyzz_xxxxy[k];

                g_0_x_xyyyyzz_xxxy[k] = -g_0_x_xyyyzz_xxxy[k] * ab_y + g_0_x_xyyyzz_xxxyy[k];

                g_0_x_xyyyyzz_xxxz[k] = -g_0_x_xyyyzz_xxxz[k] * ab_y + g_0_x_xyyyzz_xxxyz[k];

                g_0_x_xyyyyzz_xxyy[k] = -g_0_x_xyyyzz_xxyy[k] * ab_y + g_0_x_xyyyzz_xxyyy[k];

                g_0_x_xyyyyzz_xxyz[k] = -g_0_x_xyyyzz_xxyz[k] * ab_y + g_0_x_xyyyzz_xxyyz[k];

                g_0_x_xyyyyzz_xxzz[k] = -g_0_x_xyyyzz_xxzz[k] * ab_y + g_0_x_xyyyzz_xxyzz[k];

                g_0_x_xyyyyzz_xyyy[k] = -g_0_x_xyyyzz_xyyy[k] * ab_y + g_0_x_xyyyzz_xyyyy[k];

                g_0_x_xyyyyzz_xyyz[k] = -g_0_x_xyyyzz_xyyz[k] * ab_y + g_0_x_xyyyzz_xyyyz[k];

                g_0_x_xyyyyzz_xyzz[k] = -g_0_x_xyyyzz_xyzz[k] * ab_y + g_0_x_xyyyzz_xyyzz[k];

                g_0_x_xyyyyzz_xzzz[k] = -g_0_x_xyyyzz_xzzz[k] * ab_y + g_0_x_xyyyzz_xyzzz[k];

                g_0_x_xyyyyzz_yyyy[k] = -g_0_x_xyyyzz_yyyy[k] * ab_y + g_0_x_xyyyzz_yyyyy[k];

                g_0_x_xyyyyzz_yyyz[k] = -g_0_x_xyyyzz_yyyz[k] * ab_y + g_0_x_xyyyzz_yyyyz[k];

                g_0_x_xyyyyzz_yyzz[k] = -g_0_x_xyyyzz_yyzz[k] * ab_y + g_0_x_xyyyzz_yyyzz[k];

                g_0_x_xyyyyzz_yzzz[k] = -g_0_x_xyyyzz_yzzz[k] * ab_y + g_0_x_xyyyzz_yyzzz[k];

                g_0_x_xyyyyzz_zzzz[k] = -g_0_x_xyyyzz_zzzz[k] * ab_y + g_0_x_xyyyzz_yzzzz[k];
            }

            /// Set up 360-375 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyyzzz_xxxx = cbuffer.data(kg_geom_01_off + 360 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xxxy = cbuffer.data(kg_geom_01_off + 361 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xxxz = cbuffer.data(kg_geom_01_off + 362 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xxyy = cbuffer.data(kg_geom_01_off + 363 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xxyz = cbuffer.data(kg_geom_01_off + 364 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xxzz = cbuffer.data(kg_geom_01_off + 365 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xyyy = cbuffer.data(kg_geom_01_off + 366 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xyyz = cbuffer.data(kg_geom_01_off + 367 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xyzz = cbuffer.data(kg_geom_01_off + 368 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_xzzz = cbuffer.data(kg_geom_01_off + 369 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_yyyy = cbuffer.data(kg_geom_01_off + 370 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_yyyz = cbuffer.data(kg_geom_01_off + 371 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_yyzz = cbuffer.data(kg_geom_01_off + 372 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_yzzz = cbuffer.data(kg_geom_01_off + 373 * ccomps * dcomps);

            auto g_0_x_xyyyzzz_zzzz = cbuffer.data(kg_geom_01_off + 374 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyyzzz_xxxx, g_0_x_xyyyzzz_xxxy, g_0_x_xyyyzzz_xxxz, g_0_x_xyyyzzz_xxyy, g_0_x_xyyyzzz_xxyz, g_0_x_xyyyzzz_xxzz, g_0_x_xyyyzzz_xyyy, g_0_x_xyyyzzz_xyyz, g_0_x_xyyyzzz_xyzz, g_0_x_xyyyzzz_xzzz, g_0_x_xyyyzzz_yyyy, g_0_x_xyyyzzz_yyyz, g_0_x_xyyyzzz_yyzz, g_0_x_xyyyzzz_yzzz, g_0_x_xyyyzzz_zzzz, g_0_x_xyyzzz_xxxx, g_0_x_xyyzzz_xxxxy, g_0_x_xyyzzz_xxxy, g_0_x_xyyzzz_xxxyy, g_0_x_xyyzzz_xxxyz, g_0_x_xyyzzz_xxxz, g_0_x_xyyzzz_xxyy, g_0_x_xyyzzz_xxyyy, g_0_x_xyyzzz_xxyyz, g_0_x_xyyzzz_xxyz, g_0_x_xyyzzz_xxyzz, g_0_x_xyyzzz_xxzz, g_0_x_xyyzzz_xyyy, g_0_x_xyyzzz_xyyyy, g_0_x_xyyzzz_xyyyz, g_0_x_xyyzzz_xyyz, g_0_x_xyyzzz_xyyzz, g_0_x_xyyzzz_xyzz, g_0_x_xyyzzz_xyzzz, g_0_x_xyyzzz_xzzz, g_0_x_xyyzzz_yyyy, g_0_x_xyyzzz_yyyyy, g_0_x_xyyzzz_yyyyz, g_0_x_xyyzzz_yyyz, g_0_x_xyyzzz_yyyzz, g_0_x_xyyzzz_yyzz, g_0_x_xyyzzz_yyzzz, g_0_x_xyyzzz_yzzz, g_0_x_xyyzzz_yzzzz, g_0_x_xyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyyzzz_xxxx[k] = -g_0_x_xyyzzz_xxxx[k] * ab_y + g_0_x_xyyzzz_xxxxy[k];

                g_0_x_xyyyzzz_xxxy[k] = -g_0_x_xyyzzz_xxxy[k] * ab_y + g_0_x_xyyzzz_xxxyy[k];

                g_0_x_xyyyzzz_xxxz[k] = -g_0_x_xyyzzz_xxxz[k] * ab_y + g_0_x_xyyzzz_xxxyz[k];

                g_0_x_xyyyzzz_xxyy[k] = -g_0_x_xyyzzz_xxyy[k] * ab_y + g_0_x_xyyzzz_xxyyy[k];

                g_0_x_xyyyzzz_xxyz[k] = -g_0_x_xyyzzz_xxyz[k] * ab_y + g_0_x_xyyzzz_xxyyz[k];

                g_0_x_xyyyzzz_xxzz[k] = -g_0_x_xyyzzz_xxzz[k] * ab_y + g_0_x_xyyzzz_xxyzz[k];

                g_0_x_xyyyzzz_xyyy[k] = -g_0_x_xyyzzz_xyyy[k] * ab_y + g_0_x_xyyzzz_xyyyy[k];

                g_0_x_xyyyzzz_xyyz[k] = -g_0_x_xyyzzz_xyyz[k] * ab_y + g_0_x_xyyzzz_xyyyz[k];

                g_0_x_xyyyzzz_xyzz[k] = -g_0_x_xyyzzz_xyzz[k] * ab_y + g_0_x_xyyzzz_xyyzz[k];

                g_0_x_xyyyzzz_xzzz[k] = -g_0_x_xyyzzz_xzzz[k] * ab_y + g_0_x_xyyzzz_xyzzz[k];

                g_0_x_xyyyzzz_yyyy[k] = -g_0_x_xyyzzz_yyyy[k] * ab_y + g_0_x_xyyzzz_yyyyy[k];

                g_0_x_xyyyzzz_yyyz[k] = -g_0_x_xyyzzz_yyyz[k] * ab_y + g_0_x_xyyzzz_yyyyz[k];

                g_0_x_xyyyzzz_yyzz[k] = -g_0_x_xyyzzz_yyzz[k] * ab_y + g_0_x_xyyzzz_yyyzz[k];

                g_0_x_xyyyzzz_yzzz[k] = -g_0_x_xyyzzz_yzzz[k] * ab_y + g_0_x_xyyzzz_yyzzz[k];

                g_0_x_xyyyzzz_zzzz[k] = -g_0_x_xyyzzz_zzzz[k] * ab_y + g_0_x_xyyzzz_yzzzz[k];
            }

            /// Set up 375-390 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyyzzzz_xxxx = cbuffer.data(kg_geom_01_off + 375 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xxxy = cbuffer.data(kg_geom_01_off + 376 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xxxz = cbuffer.data(kg_geom_01_off + 377 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xxyy = cbuffer.data(kg_geom_01_off + 378 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xxyz = cbuffer.data(kg_geom_01_off + 379 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xxzz = cbuffer.data(kg_geom_01_off + 380 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xyyy = cbuffer.data(kg_geom_01_off + 381 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xyyz = cbuffer.data(kg_geom_01_off + 382 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xyzz = cbuffer.data(kg_geom_01_off + 383 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_xzzz = cbuffer.data(kg_geom_01_off + 384 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_yyyy = cbuffer.data(kg_geom_01_off + 385 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_yyyz = cbuffer.data(kg_geom_01_off + 386 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_yyzz = cbuffer.data(kg_geom_01_off + 387 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_yzzz = cbuffer.data(kg_geom_01_off + 388 * ccomps * dcomps);

            auto g_0_x_xyyzzzz_zzzz = cbuffer.data(kg_geom_01_off + 389 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyyzzzz_xxxx, g_0_x_xyyzzzz_xxxy, g_0_x_xyyzzzz_xxxz, g_0_x_xyyzzzz_xxyy, g_0_x_xyyzzzz_xxyz, g_0_x_xyyzzzz_xxzz, g_0_x_xyyzzzz_xyyy, g_0_x_xyyzzzz_xyyz, g_0_x_xyyzzzz_xyzz, g_0_x_xyyzzzz_xzzz, g_0_x_xyyzzzz_yyyy, g_0_x_xyyzzzz_yyyz, g_0_x_xyyzzzz_yyzz, g_0_x_xyyzzzz_yzzz, g_0_x_xyyzzzz_zzzz, g_0_x_xyzzzz_xxxx, g_0_x_xyzzzz_xxxxy, g_0_x_xyzzzz_xxxy, g_0_x_xyzzzz_xxxyy, g_0_x_xyzzzz_xxxyz, g_0_x_xyzzzz_xxxz, g_0_x_xyzzzz_xxyy, g_0_x_xyzzzz_xxyyy, g_0_x_xyzzzz_xxyyz, g_0_x_xyzzzz_xxyz, g_0_x_xyzzzz_xxyzz, g_0_x_xyzzzz_xxzz, g_0_x_xyzzzz_xyyy, g_0_x_xyzzzz_xyyyy, g_0_x_xyzzzz_xyyyz, g_0_x_xyzzzz_xyyz, g_0_x_xyzzzz_xyyzz, g_0_x_xyzzzz_xyzz, g_0_x_xyzzzz_xyzzz, g_0_x_xyzzzz_xzzz, g_0_x_xyzzzz_yyyy, g_0_x_xyzzzz_yyyyy, g_0_x_xyzzzz_yyyyz, g_0_x_xyzzzz_yyyz, g_0_x_xyzzzz_yyyzz, g_0_x_xyzzzz_yyzz, g_0_x_xyzzzz_yyzzz, g_0_x_xyzzzz_yzzz, g_0_x_xyzzzz_yzzzz, g_0_x_xyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyyzzzz_xxxx[k] = -g_0_x_xyzzzz_xxxx[k] * ab_y + g_0_x_xyzzzz_xxxxy[k];

                g_0_x_xyyzzzz_xxxy[k] = -g_0_x_xyzzzz_xxxy[k] * ab_y + g_0_x_xyzzzz_xxxyy[k];

                g_0_x_xyyzzzz_xxxz[k] = -g_0_x_xyzzzz_xxxz[k] * ab_y + g_0_x_xyzzzz_xxxyz[k];

                g_0_x_xyyzzzz_xxyy[k] = -g_0_x_xyzzzz_xxyy[k] * ab_y + g_0_x_xyzzzz_xxyyy[k];

                g_0_x_xyyzzzz_xxyz[k] = -g_0_x_xyzzzz_xxyz[k] * ab_y + g_0_x_xyzzzz_xxyyz[k];

                g_0_x_xyyzzzz_xxzz[k] = -g_0_x_xyzzzz_xxzz[k] * ab_y + g_0_x_xyzzzz_xxyzz[k];

                g_0_x_xyyzzzz_xyyy[k] = -g_0_x_xyzzzz_xyyy[k] * ab_y + g_0_x_xyzzzz_xyyyy[k];

                g_0_x_xyyzzzz_xyyz[k] = -g_0_x_xyzzzz_xyyz[k] * ab_y + g_0_x_xyzzzz_xyyyz[k];

                g_0_x_xyyzzzz_xyzz[k] = -g_0_x_xyzzzz_xyzz[k] * ab_y + g_0_x_xyzzzz_xyyzz[k];

                g_0_x_xyyzzzz_xzzz[k] = -g_0_x_xyzzzz_xzzz[k] * ab_y + g_0_x_xyzzzz_xyzzz[k];

                g_0_x_xyyzzzz_yyyy[k] = -g_0_x_xyzzzz_yyyy[k] * ab_y + g_0_x_xyzzzz_yyyyy[k];

                g_0_x_xyyzzzz_yyyz[k] = -g_0_x_xyzzzz_yyyz[k] * ab_y + g_0_x_xyzzzz_yyyyz[k];

                g_0_x_xyyzzzz_yyzz[k] = -g_0_x_xyzzzz_yyzz[k] * ab_y + g_0_x_xyzzzz_yyyzz[k];

                g_0_x_xyyzzzz_yzzz[k] = -g_0_x_xyzzzz_yzzz[k] * ab_y + g_0_x_xyzzzz_yyzzz[k];

                g_0_x_xyyzzzz_zzzz[k] = -g_0_x_xyzzzz_zzzz[k] * ab_y + g_0_x_xyzzzz_yzzzz[k];
            }

            /// Set up 390-405 components of targeted buffer : cbuffer.data(

            auto g_0_x_xyzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 390 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 391 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 392 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 393 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 394 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 395 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 396 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 397 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 398 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 399 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 400 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 401 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 402 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 403 * ccomps * dcomps);

            auto g_0_x_xyzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 404 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xyzzzzz_xxxx, g_0_x_xyzzzzz_xxxy, g_0_x_xyzzzzz_xxxz, g_0_x_xyzzzzz_xxyy, g_0_x_xyzzzzz_xxyz, g_0_x_xyzzzzz_xxzz, g_0_x_xyzzzzz_xyyy, g_0_x_xyzzzzz_xyyz, g_0_x_xyzzzzz_xyzz, g_0_x_xyzzzzz_xzzz, g_0_x_xyzzzzz_yyyy, g_0_x_xyzzzzz_yyyz, g_0_x_xyzzzzz_yyzz, g_0_x_xyzzzzz_yzzz, g_0_x_xyzzzzz_zzzz, g_0_x_xzzzzz_xxxx, g_0_x_xzzzzz_xxxxy, g_0_x_xzzzzz_xxxy, g_0_x_xzzzzz_xxxyy, g_0_x_xzzzzz_xxxyz, g_0_x_xzzzzz_xxxz, g_0_x_xzzzzz_xxyy, g_0_x_xzzzzz_xxyyy, g_0_x_xzzzzz_xxyyz, g_0_x_xzzzzz_xxyz, g_0_x_xzzzzz_xxyzz, g_0_x_xzzzzz_xxzz, g_0_x_xzzzzz_xyyy, g_0_x_xzzzzz_xyyyy, g_0_x_xzzzzz_xyyyz, g_0_x_xzzzzz_xyyz, g_0_x_xzzzzz_xyyzz, g_0_x_xzzzzz_xyzz, g_0_x_xzzzzz_xyzzz, g_0_x_xzzzzz_xzzz, g_0_x_xzzzzz_yyyy, g_0_x_xzzzzz_yyyyy, g_0_x_xzzzzz_yyyyz, g_0_x_xzzzzz_yyyz, g_0_x_xzzzzz_yyyzz, g_0_x_xzzzzz_yyzz, g_0_x_xzzzzz_yyzzz, g_0_x_xzzzzz_yzzz, g_0_x_xzzzzz_yzzzz, g_0_x_xzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xyzzzzz_xxxx[k] = -g_0_x_xzzzzz_xxxx[k] * ab_y + g_0_x_xzzzzz_xxxxy[k];

                g_0_x_xyzzzzz_xxxy[k] = -g_0_x_xzzzzz_xxxy[k] * ab_y + g_0_x_xzzzzz_xxxyy[k];

                g_0_x_xyzzzzz_xxxz[k] = -g_0_x_xzzzzz_xxxz[k] * ab_y + g_0_x_xzzzzz_xxxyz[k];

                g_0_x_xyzzzzz_xxyy[k] = -g_0_x_xzzzzz_xxyy[k] * ab_y + g_0_x_xzzzzz_xxyyy[k];

                g_0_x_xyzzzzz_xxyz[k] = -g_0_x_xzzzzz_xxyz[k] * ab_y + g_0_x_xzzzzz_xxyyz[k];

                g_0_x_xyzzzzz_xxzz[k] = -g_0_x_xzzzzz_xxzz[k] * ab_y + g_0_x_xzzzzz_xxyzz[k];

                g_0_x_xyzzzzz_xyyy[k] = -g_0_x_xzzzzz_xyyy[k] * ab_y + g_0_x_xzzzzz_xyyyy[k];

                g_0_x_xyzzzzz_xyyz[k] = -g_0_x_xzzzzz_xyyz[k] * ab_y + g_0_x_xzzzzz_xyyyz[k];

                g_0_x_xyzzzzz_xyzz[k] = -g_0_x_xzzzzz_xyzz[k] * ab_y + g_0_x_xzzzzz_xyyzz[k];

                g_0_x_xyzzzzz_xzzz[k] = -g_0_x_xzzzzz_xzzz[k] * ab_y + g_0_x_xzzzzz_xyzzz[k];

                g_0_x_xyzzzzz_yyyy[k] = -g_0_x_xzzzzz_yyyy[k] * ab_y + g_0_x_xzzzzz_yyyyy[k];

                g_0_x_xyzzzzz_yyyz[k] = -g_0_x_xzzzzz_yyyz[k] * ab_y + g_0_x_xzzzzz_yyyyz[k];

                g_0_x_xyzzzzz_yyzz[k] = -g_0_x_xzzzzz_yyzz[k] * ab_y + g_0_x_xzzzzz_yyyzz[k];

                g_0_x_xyzzzzz_yzzz[k] = -g_0_x_xzzzzz_yzzz[k] * ab_y + g_0_x_xzzzzz_yyzzz[k];

                g_0_x_xyzzzzz_zzzz[k] = -g_0_x_xzzzzz_zzzz[k] * ab_y + g_0_x_xzzzzz_yzzzz[k];
            }

            /// Set up 405-420 components of targeted buffer : cbuffer.data(

            auto g_0_x_xzzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 405 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 406 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 407 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 408 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 409 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 410 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 411 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 412 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 413 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 414 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 415 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 416 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 417 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 418 * ccomps * dcomps);

            auto g_0_x_xzzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 419 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_xzzzzz_xxxx, g_0_x_xzzzzz_xxxxz, g_0_x_xzzzzz_xxxy, g_0_x_xzzzzz_xxxyz, g_0_x_xzzzzz_xxxz, g_0_x_xzzzzz_xxxzz, g_0_x_xzzzzz_xxyy, g_0_x_xzzzzz_xxyyz, g_0_x_xzzzzz_xxyz, g_0_x_xzzzzz_xxyzz, g_0_x_xzzzzz_xxzz, g_0_x_xzzzzz_xxzzz, g_0_x_xzzzzz_xyyy, g_0_x_xzzzzz_xyyyz, g_0_x_xzzzzz_xyyz, g_0_x_xzzzzz_xyyzz, g_0_x_xzzzzz_xyzz, g_0_x_xzzzzz_xyzzz, g_0_x_xzzzzz_xzzz, g_0_x_xzzzzz_xzzzz, g_0_x_xzzzzz_yyyy, g_0_x_xzzzzz_yyyyz, g_0_x_xzzzzz_yyyz, g_0_x_xzzzzz_yyyzz, g_0_x_xzzzzz_yyzz, g_0_x_xzzzzz_yyzzz, g_0_x_xzzzzz_yzzz, g_0_x_xzzzzz_yzzzz, g_0_x_xzzzzz_zzzz, g_0_x_xzzzzz_zzzzz, g_0_x_xzzzzzz_xxxx, g_0_x_xzzzzzz_xxxy, g_0_x_xzzzzzz_xxxz, g_0_x_xzzzzzz_xxyy, g_0_x_xzzzzzz_xxyz, g_0_x_xzzzzzz_xxzz, g_0_x_xzzzzzz_xyyy, g_0_x_xzzzzzz_xyyz, g_0_x_xzzzzzz_xyzz, g_0_x_xzzzzzz_xzzz, g_0_x_xzzzzzz_yyyy, g_0_x_xzzzzzz_yyyz, g_0_x_xzzzzzz_yyzz, g_0_x_xzzzzzz_yzzz, g_0_x_xzzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_xzzzzzz_xxxx[k] = -g_0_x_xzzzzz_xxxx[k] * ab_z + g_0_x_xzzzzz_xxxxz[k];

                g_0_x_xzzzzzz_xxxy[k] = -g_0_x_xzzzzz_xxxy[k] * ab_z + g_0_x_xzzzzz_xxxyz[k];

                g_0_x_xzzzzzz_xxxz[k] = -g_0_x_xzzzzz_xxxz[k] * ab_z + g_0_x_xzzzzz_xxxzz[k];

                g_0_x_xzzzzzz_xxyy[k] = -g_0_x_xzzzzz_xxyy[k] * ab_z + g_0_x_xzzzzz_xxyyz[k];

                g_0_x_xzzzzzz_xxyz[k] = -g_0_x_xzzzzz_xxyz[k] * ab_z + g_0_x_xzzzzz_xxyzz[k];

                g_0_x_xzzzzzz_xxzz[k] = -g_0_x_xzzzzz_xxzz[k] * ab_z + g_0_x_xzzzzz_xxzzz[k];

                g_0_x_xzzzzzz_xyyy[k] = -g_0_x_xzzzzz_xyyy[k] * ab_z + g_0_x_xzzzzz_xyyyz[k];

                g_0_x_xzzzzzz_xyyz[k] = -g_0_x_xzzzzz_xyyz[k] * ab_z + g_0_x_xzzzzz_xyyzz[k];

                g_0_x_xzzzzzz_xyzz[k] = -g_0_x_xzzzzz_xyzz[k] * ab_z + g_0_x_xzzzzz_xyzzz[k];

                g_0_x_xzzzzzz_xzzz[k] = -g_0_x_xzzzzz_xzzz[k] * ab_z + g_0_x_xzzzzz_xzzzz[k];

                g_0_x_xzzzzzz_yyyy[k] = -g_0_x_xzzzzz_yyyy[k] * ab_z + g_0_x_xzzzzz_yyyyz[k];

                g_0_x_xzzzzzz_yyyz[k] = -g_0_x_xzzzzz_yyyz[k] * ab_z + g_0_x_xzzzzz_yyyzz[k];

                g_0_x_xzzzzzz_yyzz[k] = -g_0_x_xzzzzz_yyzz[k] * ab_z + g_0_x_xzzzzz_yyzzz[k];

                g_0_x_xzzzzzz_yzzz[k] = -g_0_x_xzzzzz_yzzz[k] * ab_z + g_0_x_xzzzzz_yzzzz[k];

                g_0_x_xzzzzzz_zzzz[k] = -g_0_x_xzzzzz_zzzz[k] * ab_z + g_0_x_xzzzzz_zzzzz[k];
            }

            /// Set up 420-435 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyyy_xxxx = cbuffer.data(kg_geom_01_off + 420 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xxxy = cbuffer.data(kg_geom_01_off + 421 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xxxz = cbuffer.data(kg_geom_01_off + 422 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xxyy = cbuffer.data(kg_geom_01_off + 423 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xxyz = cbuffer.data(kg_geom_01_off + 424 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xxzz = cbuffer.data(kg_geom_01_off + 425 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xyyy = cbuffer.data(kg_geom_01_off + 426 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xyyz = cbuffer.data(kg_geom_01_off + 427 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xyzz = cbuffer.data(kg_geom_01_off + 428 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_xzzz = cbuffer.data(kg_geom_01_off + 429 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_yyyy = cbuffer.data(kg_geom_01_off + 430 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_yyyz = cbuffer.data(kg_geom_01_off + 431 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_yyzz = cbuffer.data(kg_geom_01_off + 432 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_yzzz = cbuffer.data(kg_geom_01_off + 433 * ccomps * dcomps);

            auto g_0_x_yyyyyyy_zzzz = cbuffer.data(kg_geom_01_off + 434 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyy_xxxx, g_0_x_yyyyyy_xxxxy, g_0_x_yyyyyy_xxxy, g_0_x_yyyyyy_xxxyy, g_0_x_yyyyyy_xxxyz, g_0_x_yyyyyy_xxxz, g_0_x_yyyyyy_xxyy, g_0_x_yyyyyy_xxyyy, g_0_x_yyyyyy_xxyyz, g_0_x_yyyyyy_xxyz, g_0_x_yyyyyy_xxyzz, g_0_x_yyyyyy_xxzz, g_0_x_yyyyyy_xyyy, g_0_x_yyyyyy_xyyyy, g_0_x_yyyyyy_xyyyz, g_0_x_yyyyyy_xyyz, g_0_x_yyyyyy_xyyzz, g_0_x_yyyyyy_xyzz, g_0_x_yyyyyy_xyzzz, g_0_x_yyyyyy_xzzz, g_0_x_yyyyyy_yyyy, g_0_x_yyyyyy_yyyyy, g_0_x_yyyyyy_yyyyz, g_0_x_yyyyyy_yyyz, g_0_x_yyyyyy_yyyzz, g_0_x_yyyyyy_yyzz, g_0_x_yyyyyy_yyzzz, g_0_x_yyyyyy_yzzz, g_0_x_yyyyyy_yzzzz, g_0_x_yyyyyy_zzzz, g_0_x_yyyyyyy_xxxx, g_0_x_yyyyyyy_xxxy, g_0_x_yyyyyyy_xxxz, g_0_x_yyyyyyy_xxyy, g_0_x_yyyyyyy_xxyz, g_0_x_yyyyyyy_xxzz, g_0_x_yyyyyyy_xyyy, g_0_x_yyyyyyy_xyyz, g_0_x_yyyyyyy_xyzz, g_0_x_yyyyyyy_xzzz, g_0_x_yyyyyyy_yyyy, g_0_x_yyyyyyy_yyyz, g_0_x_yyyyyyy_yyzz, g_0_x_yyyyyyy_yzzz, g_0_x_yyyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyyy_xxxx[k] = -g_0_x_yyyyyy_xxxx[k] * ab_y + g_0_x_yyyyyy_xxxxy[k];

                g_0_x_yyyyyyy_xxxy[k] = -g_0_x_yyyyyy_xxxy[k] * ab_y + g_0_x_yyyyyy_xxxyy[k];

                g_0_x_yyyyyyy_xxxz[k] = -g_0_x_yyyyyy_xxxz[k] * ab_y + g_0_x_yyyyyy_xxxyz[k];

                g_0_x_yyyyyyy_xxyy[k] = -g_0_x_yyyyyy_xxyy[k] * ab_y + g_0_x_yyyyyy_xxyyy[k];

                g_0_x_yyyyyyy_xxyz[k] = -g_0_x_yyyyyy_xxyz[k] * ab_y + g_0_x_yyyyyy_xxyyz[k];

                g_0_x_yyyyyyy_xxzz[k] = -g_0_x_yyyyyy_xxzz[k] * ab_y + g_0_x_yyyyyy_xxyzz[k];

                g_0_x_yyyyyyy_xyyy[k] = -g_0_x_yyyyyy_xyyy[k] * ab_y + g_0_x_yyyyyy_xyyyy[k];

                g_0_x_yyyyyyy_xyyz[k] = -g_0_x_yyyyyy_xyyz[k] * ab_y + g_0_x_yyyyyy_xyyyz[k];

                g_0_x_yyyyyyy_xyzz[k] = -g_0_x_yyyyyy_xyzz[k] * ab_y + g_0_x_yyyyyy_xyyzz[k];

                g_0_x_yyyyyyy_xzzz[k] = -g_0_x_yyyyyy_xzzz[k] * ab_y + g_0_x_yyyyyy_xyzzz[k];

                g_0_x_yyyyyyy_yyyy[k] = -g_0_x_yyyyyy_yyyy[k] * ab_y + g_0_x_yyyyyy_yyyyy[k];

                g_0_x_yyyyyyy_yyyz[k] = -g_0_x_yyyyyy_yyyz[k] * ab_y + g_0_x_yyyyyy_yyyyz[k];

                g_0_x_yyyyyyy_yyzz[k] = -g_0_x_yyyyyy_yyzz[k] * ab_y + g_0_x_yyyyyy_yyyzz[k];

                g_0_x_yyyyyyy_yzzz[k] = -g_0_x_yyyyyy_yzzz[k] * ab_y + g_0_x_yyyyyy_yyzzz[k];

                g_0_x_yyyyyyy_zzzz[k] = -g_0_x_yyyyyy_zzzz[k] * ab_y + g_0_x_yyyyyy_yzzzz[k];
            }

            /// Set up 435-450 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyyz_xxxx = cbuffer.data(kg_geom_01_off + 435 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xxxy = cbuffer.data(kg_geom_01_off + 436 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xxxz = cbuffer.data(kg_geom_01_off + 437 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xxyy = cbuffer.data(kg_geom_01_off + 438 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xxyz = cbuffer.data(kg_geom_01_off + 439 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xxzz = cbuffer.data(kg_geom_01_off + 440 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xyyy = cbuffer.data(kg_geom_01_off + 441 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xyyz = cbuffer.data(kg_geom_01_off + 442 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xyzz = cbuffer.data(kg_geom_01_off + 443 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_xzzz = cbuffer.data(kg_geom_01_off + 444 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_yyyy = cbuffer.data(kg_geom_01_off + 445 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_yyyz = cbuffer.data(kg_geom_01_off + 446 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_yyzz = cbuffer.data(kg_geom_01_off + 447 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_yzzz = cbuffer.data(kg_geom_01_off + 448 * ccomps * dcomps);

            auto g_0_x_yyyyyyz_zzzz = cbuffer.data(kg_geom_01_off + 449 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyyz_xxxx, g_0_x_yyyyyyz_xxxy, g_0_x_yyyyyyz_xxxz, g_0_x_yyyyyyz_xxyy, g_0_x_yyyyyyz_xxyz, g_0_x_yyyyyyz_xxzz, g_0_x_yyyyyyz_xyyy, g_0_x_yyyyyyz_xyyz, g_0_x_yyyyyyz_xyzz, g_0_x_yyyyyyz_xzzz, g_0_x_yyyyyyz_yyyy, g_0_x_yyyyyyz_yyyz, g_0_x_yyyyyyz_yyzz, g_0_x_yyyyyyz_yzzz, g_0_x_yyyyyyz_zzzz, g_0_x_yyyyyz_xxxx, g_0_x_yyyyyz_xxxxy, g_0_x_yyyyyz_xxxy, g_0_x_yyyyyz_xxxyy, g_0_x_yyyyyz_xxxyz, g_0_x_yyyyyz_xxxz, g_0_x_yyyyyz_xxyy, g_0_x_yyyyyz_xxyyy, g_0_x_yyyyyz_xxyyz, g_0_x_yyyyyz_xxyz, g_0_x_yyyyyz_xxyzz, g_0_x_yyyyyz_xxzz, g_0_x_yyyyyz_xyyy, g_0_x_yyyyyz_xyyyy, g_0_x_yyyyyz_xyyyz, g_0_x_yyyyyz_xyyz, g_0_x_yyyyyz_xyyzz, g_0_x_yyyyyz_xyzz, g_0_x_yyyyyz_xyzzz, g_0_x_yyyyyz_xzzz, g_0_x_yyyyyz_yyyy, g_0_x_yyyyyz_yyyyy, g_0_x_yyyyyz_yyyyz, g_0_x_yyyyyz_yyyz, g_0_x_yyyyyz_yyyzz, g_0_x_yyyyyz_yyzz, g_0_x_yyyyyz_yyzzz, g_0_x_yyyyyz_yzzz, g_0_x_yyyyyz_yzzzz, g_0_x_yyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyyz_xxxx[k] = -g_0_x_yyyyyz_xxxx[k] * ab_y + g_0_x_yyyyyz_xxxxy[k];

                g_0_x_yyyyyyz_xxxy[k] = -g_0_x_yyyyyz_xxxy[k] * ab_y + g_0_x_yyyyyz_xxxyy[k];

                g_0_x_yyyyyyz_xxxz[k] = -g_0_x_yyyyyz_xxxz[k] * ab_y + g_0_x_yyyyyz_xxxyz[k];

                g_0_x_yyyyyyz_xxyy[k] = -g_0_x_yyyyyz_xxyy[k] * ab_y + g_0_x_yyyyyz_xxyyy[k];

                g_0_x_yyyyyyz_xxyz[k] = -g_0_x_yyyyyz_xxyz[k] * ab_y + g_0_x_yyyyyz_xxyyz[k];

                g_0_x_yyyyyyz_xxzz[k] = -g_0_x_yyyyyz_xxzz[k] * ab_y + g_0_x_yyyyyz_xxyzz[k];

                g_0_x_yyyyyyz_xyyy[k] = -g_0_x_yyyyyz_xyyy[k] * ab_y + g_0_x_yyyyyz_xyyyy[k];

                g_0_x_yyyyyyz_xyyz[k] = -g_0_x_yyyyyz_xyyz[k] * ab_y + g_0_x_yyyyyz_xyyyz[k];

                g_0_x_yyyyyyz_xyzz[k] = -g_0_x_yyyyyz_xyzz[k] * ab_y + g_0_x_yyyyyz_xyyzz[k];

                g_0_x_yyyyyyz_xzzz[k] = -g_0_x_yyyyyz_xzzz[k] * ab_y + g_0_x_yyyyyz_xyzzz[k];

                g_0_x_yyyyyyz_yyyy[k] = -g_0_x_yyyyyz_yyyy[k] * ab_y + g_0_x_yyyyyz_yyyyy[k];

                g_0_x_yyyyyyz_yyyz[k] = -g_0_x_yyyyyz_yyyz[k] * ab_y + g_0_x_yyyyyz_yyyyz[k];

                g_0_x_yyyyyyz_yyzz[k] = -g_0_x_yyyyyz_yyzz[k] * ab_y + g_0_x_yyyyyz_yyyzz[k];

                g_0_x_yyyyyyz_yzzz[k] = -g_0_x_yyyyyz_yzzz[k] * ab_y + g_0_x_yyyyyz_yyzzz[k];

                g_0_x_yyyyyyz_zzzz[k] = -g_0_x_yyyyyz_zzzz[k] * ab_y + g_0_x_yyyyyz_yzzzz[k];
            }

            /// Set up 450-465 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyyzz_xxxx = cbuffer.data(kg_geom_01_off + 450 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xxxy = cbuffer.data(kg_geom_01_off + 451 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xxxz = cbuffer.data(kg_geom_01_off + 452 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xxyy = cbuffer.data(kg_geom_01_off + 453 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xxyz = cbuffer.data(kg_geom_01_off + 454 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xxzz = cbuffer.data(kg_geom_01_off + 455 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xyyy = cbuffer.data(kg_geom_01_off + 456 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xyyz = cbuffer.data(kg_geom_01_off + 457 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xyzz = cbuffer.data(kg_geom_01_off + 458 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_xzzz = cbuffer.data(kg_geom_01_off + 459 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_yyyy = cbuffer.data(kg_geom_01_off + 460 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_yyyz = cbuffer.data(kg_geom_01_off + 461 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_yyzz = cbuffer.data(kg_geom_01_off + 462 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_yzzz = cbuffer.data(kg_geom_01_off + 463 * ccomps * dcomps);

            auto g_0_x_yyyyyzz_zzzz = cbuffer.data(kg_geom_01_off + 464 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyyzz_xxxx, g_0_x_yyyyyzz_xxxy, g_0_x_yyyyyzz_xxxz, g_0_x_yyyyyzz_xxyy, g_0_x_yyyyyzz_xxyz, g_0_x_yyyyyzz_xxzz, g_0_x_yyyyyzz_xyyy, g_0_x_yyyyyzz_xyyz, g_0_x_yyyyyzz_xyzz, g_0_x_yyyyyzz_xzzz, g_0_x_yyyyyzz_yyyy, g_0_x_yyyyyzz_yyyz, g_0_x_yyyyyzz_yyzz, g_0_x_yyyyyzz_yzzz, g_0_x_yyyyyzz_zzzz, g_0_x_yyyyzz_xxxx, g_0_x_yyyyzz_xxxxy, g_0_x_yyyyzz_xxxy, g_0_x_yyyyzz_xxxyy, g_0_x_yyyyzz_xxxyz, g_0_x_yyyyzz_xxxz, g_0_x_yyyyzz_xxyy, g_0_x_yyyyzz_xxyyy, g_0_x_yyyyzz_xxyyz, g_0_x_yyyyzz_xxyz, g_0_x_yyyyzz_xxyzz, g_0_x_yyyyzz_xxzz, g_0_x_yyyyzz_xyyy, g_0_x_yyyyzz_xyyyy, g_0_x_yyyyzz_xyyyz, g_0_x_yyyyzz_xyyz, g_0_x_yyyyzz_xyyzz, g_0_x_yyyyzz_xyzz, g_0_x_yyyyzz_xyzzz, g_0_x_yyyyzz_xzzz, g_0_x_yyyyzz_yyyy, g_0_x_yyyyzz_yyyyy, g_0_x_yyyyzz_yyyyz, g_0_x_yyyyzz_yyyz, g_0_x_yyyyzz_yyyzz, g_0_x_yyyyzz_yyzz, g_0_x_yyyyzz_yyzzz, g_0_x_yyyyzz_yzzz, g_0_x_yyyyzz_yzzzz, g_0_x_yyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyyzz_xxxx[k] = -g_0_x_yyyyzz_xxxx[k] * ab_y + g_0_x_yyyyzz_xxxxy[k];

                g_0_x_yyyyyzz_xxxy[k] = -g_0_x_yyyyzz_xxxy[k] * ab_y + g_0_x_yyyyzz_xxxyy[k];

                g_0_x_yyyyyzz_xxxz[k] = -g_0_x_yyyyzz_xxxz[k] * ab_y + g_0_x_yyyyzz_xxxyz[k];

                g_0_x_yyyyyzz_xxyy[k] = -g_0_x_yyyyzz_xxyy[k] * ab_y + g_0_x_yyyyzz_xxyyy[k];

                g_0_x_yyyyyzz_xxyz[k] = -g_0_x_yyyyzz_xxyz[k] * ab_y + g_0_x_yyyyzz_xxyyz[k];

                g_0_x_yyyyyzz_xxzz[k] = -g_0_x_yyyyzz_xxzz[k] * ab_y + g_0_x_yyyyzz_xxyzz[k];

                g_0_x_yyyyyzz_xyyy[k] = -g_0_x_yyyyzz_xyyy[k] * ab_y + g_0_x_yyyyzz_xyyyy[k];

                g_0_x_yyyyyzz_xyyz[k] = -g_0_x_yyyyzz_xyyz[k] * ab_y + g_0_x_yyyyzz_xyyyz[k];

                g_0_x_yyyyyzz_xyzz[k] = -g_0_x_yyyyzz_xyzz[k] * ab_y + g_0_x_yyyyzz_xyyzz[k];

                g_0_x_yyyyyzz_xzzz[k] = -g_0_x_yyyyzz_xzzz[k] * ab_y + g_0_x_yyyyzz_xyzzz[k];

                g_0_x_yyyyyzz_yyyy[k] = -g_0_x_yyyyzz_yyyy[k] * ab_y + g_0_x_yyyyzz_yyyyy[k];

                g_0_x_yyyyyzz_yyyz[k] = -g_0_x_yyyyzz_yyyz[k] * ab_y + g_0_x_yyyyzz_yyyyz[k];

                g_0_x_yyyyyzz_yyzz[k] = -g_0_x_yyyyzz_yyzz[k] * ab_y + g_0_x_yyyyzz_yyyzz[k];

                g_0_x_yyyyyzz_yzzz[k] = -g_0_x_yyyyzz_yzzz[k] * ab_y + g_0_x_yyyyzz_yyzzz[k];

                g_0_x_yyyyyzz_zzzz[k] = -g_0_x_yyyyzz_zzzz[k] * ab_y + g_0_x_yyyyzz_yzzzz[k];
            }

            /// Set up 465-480 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyyzzz_xxxx = cbuffer.data(kg_geom_01_off + 465 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xxxy = cbuffer.data(kg_geom_01_off + 466 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xxxz = cbuffer.data(kg_geom_01_off + 467 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xxyy = cbuffer.data(kg_geom_01_off + 468 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xxyz = cbuffer.data(kg_geom_01_off + 469 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xxzz = cbuffer.data(kg_geom_01_off + 470 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xyyy = cbuffer.data(kg_geom_01_off + 471 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xyyz = cbuffer.data(kg_geom_01_off + 472 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xyzz = cbuffer.data(kg_geom_01_off + 473 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_xzzz = cbuffer.data(kg_geom_01_off + 474 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_yyyy = cbuffer.data(kg_geom_01_off + 475 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_yyyz = cbuffer.data(kg_geom_01_off + 476 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_yyzz = cbuffer.data(kg_geom_01_off + 477 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_yzzz = cbuffer.data(kg_geom_01_off + 478 * ccomps * dcomps);

            auto g_0_x_yyyyzzz_zzzz = cbuffer.data(kg_geom_01_off + 479 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyyzzz_xxxx, g_0_x_yyyyzzz_xxxy, g_0_x_yyyyzzz_xxxz, g_0_x_yyyyzzz_xxyy, g_0_x_yyyyzzz_xxyz, g_0_x_yyyyzzz_xxzz, g_0_x_yyyyzzz_xyyy, g_0_x_yyyyzzz_xyyz, g_0_x_yyyyzzz_xyzz, g_0_x_yyyyzzz_xzzz, g_0_x_yyyyzzz_yyyy, g_0_x_yyyyzzz_yyyz, g_0_x_yyyyzzz_yyzz, g_0_x_yyyyzzz_yzzz, g_0_x_yyyyzzz_zzzz, g_0_x_yyyzzz_xxxx, g_0_x_yyyzzz_xxxxy, g_0_x_yyyzzz_xxxy, g_0_x_yyyzzz_xxxyy, g_0_x_yyyzzz_xxxyz, g_0_x_yyyzzz_xxxz, g_0_x_yyyzzz_xxyy, g_0_x_yyyzzz_xxyyy, g_0_x_yyyzzz_xxyyz, g_0_x_yyyzzz_xxyz, g_0_x_yyyzzz_xxyzz, g_0_x_yyyzzz_xxzz, g_0_x_yyyzzz_xyyy, g_0_x_yyyzzz_xyyyy, g_0_x_yyyzzz_xyyyz, g_0_x_yyyzzz_xyyz, g_0_x_yyyzzz_xyyzz, g_0_x_yyyzzz_xyzz, g_0_x_yyyzzz_xyzzz, g_0_x_yyyzzz_xzzz, g_0_x_yyyzzz_yyyy, g_0_x_yyyzzz_yyyyy, g_0_x_yyyzzz_yyyyz, g_0_x_yyyzzz_yyyz, g_0_x_yyyzzz_yyyzz, g_0_x_yyyzzz_yyzz, g_0_x_yyyzzz_yyzzz, g_0_x_yyyzzz_yzzz, g_0_x_yyyzzz_yzzzz, g_0_x_yyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyyzzz_xxxx[k] = -g_0_x_yyyzzz_xxxx[k] * ab_y + g_0_x_yyyzzz_xxxxy[k];

                g_0_x_yyyyzzz_xxxy[k] = -g_0_x_yyyzzz_xxxy[k] * ab_y + g_0_x_yyyzzz_xxxyy[k];

                g_0_x_yyyyzzz_xxxz[k] = -g_0_x_yyyzzz_xxxz[k] * ab_y + g_0_x_yyyzzz_xxxyz[k];

                g_0_x_yyyyzzz_xxyy[k] = -g_0_x_yyyzzz_xxyy[k] * ab_y + g_0_x_yyyzzz_xxyyy[k];

                g_0_x_yyyyzzz_xxyz[k] = -g_0_x_yyyzzz_xxyz[k] * ab_y + g_0_x_yyyzzz_xxyyz[k];

                g_0_x_yyyyzzz_xxzz[k] = -g_0_x_yyyzzz_xxzz[k] * ab_y + g_0_x_yyyzzz_xxyzz[k];

                g_0_x_yyyyzzz_xyyy[k] = -g_0_x_yyyzzz_xyyy[k] * ab_y + g_0_x_yyyzzz_xyyyy[k];

                g_0_x_yyyyzzz_xyyz[k] = -g_0_x_yyyzzz_xyyz[k] * ab_y + g_0_x_yyyzzz_xyyyz[k];

                g_0_x_yyyyzzz_xyzz[k] = -g_0_x_yyyzzz_xyzz[k] * ab_y + g_0_x_yyyzzz_xyyzz[k];

                g_0_x_yyyyzzz_xzzz[k] = -g_0_x_yyyzzz_xzzz[k] * ab_y + g_0_x_yyyzzz_xyzzz[k];

                g_0_x_yyyyzzz_yyyy[k] = -g_0_x_yyyzzz_yyyy[k] * ab_y + g_0_x_yyyzzz_yyyyy[k];

                g_0_x_yyyyzzz_yyyz[k] = -g_0_x_yyyzzz_yyyz[k] * ab_y + g_0_x_yyyzzz_yyyyz[k];

                g_0_x_yyyyzzz_yyzz[k] = -g_0_x_yyyzzz_yyzz[k] * ab_y + g_0_x_yyyzzz_yyyzz[k];

                g_0_x_yyyyzzz_yzzz[k] = -g_0_x_yyyzzz_yzzz[k] * ab_y + g_0_x_yyyzzz_yyzzz[k];

                g_0_x_yyyyzzz_zzzz[k] = -g_0_x_yyyzzz_zzzz[k] * ab_y + g_0_x_yyyzzz_yzzzz[k];
            }

            /// Set up 480-495 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyyzzzz_xxxx = cbuffer.data(kg_geom_01_off + 480 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xxxy = cbuffer.data(kg_geom_01_off + 481 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xxxz = cbuffer.data(kg_geom_01_off + 482 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xxyy = cbuffer.data(kg_geom_01_off + 483 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xxyz = cbuffer.data(kg_geom_01_off + 484 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xxzz = cbuffer.data(kg_geom_01_off + 485 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xyyy = cbuffer.data(kg_geom_01_off + 486 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xyyz = cbuffer.data(kg_geom_01_off + 487 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xyzz = cbuffer.data(kg_geom_01_off + 488 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_xzzz = cbuffer.data(kg_geom_01_off + 489 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_yyyy = cbuffer.data(kg_geom_01_off + 490 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_yyyz = cbuffer.data(kg_geom_01_off + 491 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_yyzz = cbuffer.data(kg_geom_01_off + 492 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_yzzz = cbuffer.data(kg_geom_01_off + 493 * ccomps * dcomps);

            auto g_0_x_yyyzzzz_zzzz = cbuffer.data(kg_geom_01_off + 494 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyyzzzz_xxxx, g_0_x_yyyzzzz_xxxy, g_0_x_yyyzzzz_xxxz, g_0_x_yyyzzzz_xxyy, g_0_x_yyyzzzz_xxyz, g_0_x_yyyzzzz_xxzz, g_0_x_yyyzzzz_xyyy, g_0_x_yyyzzzz_xyyz, g_0_x_yyyzzzz_xyzz, g_0_x_yyyzzzz_xzzz, g_0_x_yyyzzzz_yyyy, g_0_x_yyyzzzz_yyyz, g_0_x_yyyzzzz_yyzz, g_0_x_yyyzzzz_yzzz, g_0_x_yyyzzzz_zzzz, g_0_x_yyzzzz_xxxx, g_0_x_yyzzzz_xxxxy, g_0_x_yyzzzz_xxxy, g_0_x_yyzzzz_xxxyy, g_0_x_yyzzzz_xxxyz, g_0_x_yyzzzz_xxxz, g_0_x_yyzzzz_xxyy, g_0_x_yyzzzz_xxyyy, g_0_x_yyzzzz_xxyyz, g_0_x_yyzzzz_xxyz, g_0_x_yyzzzz_xxyzz, g_0_x_yyzzzz_xxzz, g_0_x_yyzzzz_xyyy, g_0_x_yyzzzz_xyyyy, g_0_x_yyzzzz_xyyyz, g_0_x_yyzzzz_xyyz, g_0_x_yyzzzz_xyyzz, g_0_x_yyzzzz_xyzz, g_0_x_yyzzzz_xyzzz, g_0_x_yyzzzz_xzzz, g_0_x_yyzzzz_yyyy, g_0_x_yyzzzz_yyyyy, g_0_x_yyzzzz_yyyyz, g_0_x_yyzzzz_yyyz, g_0_x_yyzzzz_yyyzz, g_0_x_yyzzzz_yyzz, g_0_x_yyzzzz_yyzzz, g_0_x_yyzzzz_yzzz, g_0_x_yyzzzz_yzzzz, g_0_x_yyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyyzzzz_xxxx[k] = -g_0_x_yyzzzz_xxxx[k] * ab_y + g_0_x_yyzzzz_xxxxy[k];

                g_0_x_yyyzzzz_xxxy[k] = -g_0_x_yyzzzz_xxxy[k] * ab_y + g_0_x_yyzzzz_xxxyy[k];

                g_0_x_yyyzzzz_xxxz[k] = -g_0_x_yyzzzz_xxxz[k] * ab_y + g_0_x_yyzzzz_xxxyz[k];

                g_0_x_yyyzzzz_xxyy[k] = -g_0_x_yyzzzz_xxyy[k] * ab_y + g_0_x_yyzzzz_xxyyy[k];

                g_0_x_yyyzzzz_xxyz[k] = -g_0_x_yyzzzz_xxyz[k] * ab_y + g_0_x_yyzzzz_xxyyz[k];

                g_0_x_yyyzzzz_xxzz[k] = -g_0_x_yyzzzz_xxzz[k] * ab_y + g_0_x_yyzzzz_xxyzz[k];

                g_0_x_yyyzzzz_xyyy[k] = -g_0_x_yyzzzz_xyyy[k] * ab_y + g_0_x_yyzzzz_xyyyy[k];

                g_0_x_yyyzzzz_xyyz[k] = -g_0_x_yyzzzz_xyyz[k] * ab_y + g_0_x_yyzzzz_xyyyz[k];

                g_0_x_yyyzzzz_xyzz[k] = -g_0_x_yyzzzz_xyzz[k] * ab_y + g_0_x_yyzzzz_xyyzz[k];

                g_0_x_yyyzzzz_xzzz[k] = -g_0_x_yyzzzz_xzzz[k] * ab_y + g_0_x_yyzzzz_xyzzz[k];

                g_0_x_yyyzzzz_yyyy[k] = -g_0_x_yyzzzz_yyyy[k] * ab_y + g_0_x_yyzzzz_yyyyy[k];

                g_0_x_yyyzzzz_yyyz[k] = -g_0_x_yyzzzz_yyyz[k] * ab_y + g_0_x_yyzzzz_yyyyz[k];

                g_0_x_yyyzzzz_yyzz[k] = -g_0_x_yyzzzz_yyzz[k] * ab_y + g_0_x_yyzzzz_yyyzz[k];

                g_0_x_yyyzzzz_yzzz[k] = -g_0_x_yyzzzz_yzzz[k] * ab_y + g_0_x_yyzzzz_yyzzz[k];

                g_0_x_yyyzzzz_zzzz[k] = -g_0_x_yyzzzz_zzzz[k] * ab_y + g_0_x_yyzzzz_yzzzz[k];
            }

            /// Set up 495-510 components of targeted buffer : cbuffer.data(

            auto g_0_x_yyzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 495 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 496 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 497 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 498 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 499 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 500 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 501 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 502 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 503 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 504 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 505 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 506 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 507 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 508 * ccomps * dcomps);

            auto g_0_x_yyzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 509 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yyzzzzz_xxxx, g_0_x_yyzzzzz_xxxy, g_0_x_yyzzzzz_xxxz, g_0_x_yyzzzzz_xxyy, g_0_x_yyzzzzz_xxyz, g_0_x_yyzzzzz_xxzz, g_0_x_yyzzzzz_xyyy, g_0_x_yyzzzzz_xyyz, g_0_x_yyzzzzz_xyzz, g_0_x_yyzzzzz_xzzz, g_0_x_yyzzzzz_yyyy, g_0_x_yyzzzzz_yyyz, g_0_x_yyzzzzz_yyzz, g_0_x_yyzzzzz_yzzz, g_0_x_yyzzzzz_zzzz, g_0_x_yzzzzz_xxxx, g_0_x_yzzzzz_xxxxy, g_0_x_yzzzzz_xxxy, g_0_x_yzzzzz_xxxyy, g_0_x_yzzzzz_xxxyz, g_0_x_yzzzzz_xxxz, g_0_x_yzzzzz_xxyy, g_0_x_yzzzzz_xxyyy, g_0_x_yzzzzz_xxyyz, g_0_x_yzzzzz_xxyz, g_0_x_yzzzzz_xxyzz, g_0_x_yzzzzz_xxzz, g_0_x_yzzzzz_xyyy, g_0_x_yzzzzz_xyyyy, g_0_x_yzzzzz_xyyyz, g_0_x_yzzzzz_xyyz, g_0_x_yzzzzz_xyyzz, g_0_x_yzzzzz_xyzz, g_0_x_yzzzzz_xyzzz, g_0_x_yzzzzz_xzzz, g_0_x_yzzzzz_yyyy, g_0_x_yzzzzz_yyyyy, g_0_x_yzzzzz_yyyyz, g_0_x_yzzzzz_yyyz, g_0_x_yzzzzz_yyyzz, g_0_x_yzzzzz_yyzz, g_0_x_yzzzzz_yyzzz, g_0_x_yzzzzz_yzzz, g_0_x_yzzzzz_yzzzz, g_0_x_yzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yyzzzzz_xxxx[k] = -g_0_x_yzzzzz_xxxx[k] * ab_y + g_0_x_yzzzzz_xxxxy[k];

                g_0_x_yyzzzzz_xxxy[k] = -g_0_x_yzzzzz_xxxy[k] * ab_y + g_0_x_yzzzzz_xxxyy[k];

                g_0_x_yyzzzzz_xxxz[k] = -g_0_x_yzzzzz_xxxz[k] * ab_y + g_0_x_yzzzzz_xxxyz[k];

                g_0_x_yyzzzzz_xxyy[k] = -g_0_x_yzzzzz_xxyy[k] * ab_y + g_0_x_yzzzzz_xxyyy[k];

                g_0_x_yyzzzzz_xxyz[k] = -g_0_x_yzzzzz_xxyz[k] * ab_y + g_0_x_yzzzzz_xxyyz[k];

                g_0_x_yyzzzzz_xxzz[k] = -g_0_x_yzzzzz_xxzz[k] * ab_y + g_0_x_yzzzzz_xxyzz[k];

                g_0_x_yyzzzzz_xyyy[k] = -g_0_x_yzzzzz_xyyy[k] * ab_y + g_0_x_yzzzzz_xyyyy[k];

                g_0_x_yyzzzzz_xyyz[k] = -g_0_x_yzzzzz_xyyz[k] * ab_y + g_0_x_yzzzzz_xyyyz[k];

                g_0_x_yyzzzzz_xyzz[k] = -g_0_x_yzzzzz_xyzz[k] * ab_y + g_0_x_yzzzzz_xyyzz[k];

                g_0_x_yyzzzzz_xzzz[k] = -g_0_x_yzzzzz_xzzz[k] * ab_y + g_0_x_yzzzzz_xyzzz[k];

                g_0_x_yyzzzzz_yyyy[k] = -g_0_x_yzzzzz_yyyy[k] * ab_y + g_0_x_yzzzzz_yyyyy[k];

                g_0_x_yyzzzzz_yyyz[k] = -g_0_x_yzzzzz_yyyz[k] * ab_y + g_0_x_yzzzzz_yyyyz[k];

                g_0_x_yyzzzzz_yyzz[k] = -g_0_x_yzzzzz_yyzz[k] * ab_y + g_0_x_yzzzzz_yyyzz[k];

                g_0_x_yyzzzzz_yzzz[k] = -g_0_x_yzzzzz_yzzz[k] * ab_y + g_0_x_yzzzzz_yyzzz[k];

                g_0_x_yyzzzzz_zzzz[k] = -g_0_x_yzzzzz_zzzz[k] * ab_y + g_0_x_yzzzzz_yzzzz[k];
            }

            /// Set up 510-525 components of targeted buffer : cbuffer.data(

            auto g_0_x_yzzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 510 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 511 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 512 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 513 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 514 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 515 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 516 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 517 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 518 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 519 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 520 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 521 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 522 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 523 * ccomps * dcomps);

            auto g_0_x_yzzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 524 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_yzzzzzz_xxxx, g_0_x_yzzzzzz_xxxy, g_0_x_yzzzzzz_xxxz, g_0_x_yzzzzzz_xxyy, g_0_x_yzzzzzz_xxyz, g_0_x_yzzzzzz_xxzz, g_0_x_yzzzzzz_xyyy, g_0_x_yzzzzzz_xyyz, g_0_x_yzzzzzz_xyzz, g_0_x_yzzzzzz_xzzz, g_0_x_yzzzzzz_yyyy, g_0_x_yzzzzzz_yyyz, g_0_x_yzzzzzz_yyzz, g_0_x_yzzzzzz_yzzz, g_0_x_yzzzzzz_zzzz, g_0_x_zzzzzz_xxxx, g_0_x_zzzzzz_xxxxy, g_0_x_zzzzzz_xxxy, g_0_x_zzzzzz_xxxyy, g_0_x_zzzzzz_xxxyz, g_0_x_zzzzzz_xxxz, g_0_x_zzzzzz_xxyy, g_0_x_zzzzzz_xxyyy, g_0_x_zzzzzz_xxyyz, g_0_x_zzzzzz_xxyz, g_0_x_zzzzzz_xxyzz, g_0_x_zzzzzz_xxzz, g_0_x_zzzzzz_xyyy, g_0_x_zzzzzz_xyyyy, g_0_x_zzzzzz_xyyyz, g_0_x_zzzzzz_xyyz, g_0_x_zzzzzz_xyyzz, g_0_x_zzzzzz_xyzz, g_0_x_zzzzzz_xyzzz, g_0_x_zzzzzz_xzzz, g_0_x_zzzzzz_yyyy, g_0_x_zzzzzz_yyyyy, g_0_x_zzzzzz_yyyyz, g_0_x_zzzzzz_yyyz, g_0_x_zzzzzz_yyyzz, g_0_x_zzzzzz_yyzz, g_0_x_zzzzzz_yyzzz, g_0_x_zzzzzz_yzzz, g_0_x_zzzzzz_yzzzz, g_0_x_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_yzzzzzz_xxxx[k] = -g_0_x_zzzzzz_xxxx[k] * ab_y + g_0_x_zzzzzz_xxxxy[k];

                g_0_x_yzzzzzz_xxxy[k] = -g_0_x_zzzzzz_xxxy[k] * ab_y + g_0_x_zzzzzz_xxxyy[k];

                g_0_x_yzzzzzz_xxxz[k] = -g_0_x_zzzzzz_xxxz[k] * ab_y + g_0_x_zzzzzz_xxxyz[k];

                g_0_x_yzzzzzz_xxyy[k] = -g_0_x_zzzzzz_xxyy[k] * ab_y + g_0_x_zzzzzz_xxyyy[k];

                g_0_x_yzzzzzz_xxyz[k] = -g_0_x_zzzzzz_xxyz[k] * ab_y + g_0_x_zzzzzz_xxyyz[k];

                g_0_x_yzzzzzz_xxzz[k] = -g_0_x_zzzzzz_xxzz[k] * ab_y + g_0_x_zzzzzz_xxyzz[k];

                g_0_x_yzzzzzz_xyyy[k] = -g_0_x_zzzzzz_xyyy[k] * ab_y + g_0_x_zzzzzz_xyyyy[k];

                g_0_x_yzzzzzz_xyyz[k] = -g_0_x_zzzzzz_xyyz[k] * ab_y + g_0_x_zzzzzz_xyyyz[k];

                g_0_x_yzzzzzz_xyzz[k] = -g_0_x_zzzzzz_xyzz[k] * ab_y + g_0_x_zzzzzz_xyyzz[k];

                g_0_x_yzzzzzz_xzzz[k] = -g_0_x_zzzzzz_xzzz[k] * ab_y + g_0_x_zzzzzz_xyzzz[k];

                g_0_x_yzzzzzz_yyyy[k] = -g_0_x_zzzzzz_yyyy[k] * ab_y + g_0_x_zzzzzz_yyyyy[k];

                g_0_x_yzzzzzz_yyyz[k] = -g_0_x_zzzzzz_yyyz[k] * ab_y + g_0_x_zzzzzz_yyyyz[k];

                g_0_x_yzzzzzz_yyzz[k] = -g_0_x_zzzzzz_yyzz[k] * ab_y + g_0_x_zzzzzz_yyyzz[k];

                g_0_x_yzzzzzz_yzzz[k] = -g_0_x_zzzzzz_yzzz[k] * ab_y + g_0_x_zzzzzz_yyzzz[k];

                g_0_x_yzzzzzz_zzzz[k] = -g_0_x_zzzzzz_zzzz[k] * ab_y + g_0_x_zzzzzz_yzzzz[k];
            }

            /// Set up 525-540 components of targeted buffer : cbuffer.data(

            auto g_0_x_zzzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 525 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 526 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 527 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 528 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 529 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 530 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 531 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 532 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 533 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 534 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 535 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 536 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 537 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 538 * ccomps * dcomps);

            auto g_0_x_zzzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 539 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_x_zzzzzz_xxxx, g_0_x_zzzzzz_xxxxz, g_0_x_zzzzzz_xxxy, g_0_x_zzzzzz_xxxyz, g_0_x_zzzzzz_xxxz, g_0_x_zzzzzz_xxxzz, g_0_x_zzzzzz_xxyy, g_0_x_zzzzzz_xxyyz, g_0_x_zzzzzz_xxyz, g_0_x_zzzzzz_xxyzz, g_0_x_zzzzzz_xxzz, g_0_x_zzzzzz_xxzzz, g_0_x_zzzzzz_xyyy, g_0_x_zzzzzz_xyyyz, g_0_x_zzzzzz_xyyz, g_0_x_zzzzzz_xyyzz, g_0_x_zzzzzz_xyzz, g_0_x_zzzzzz_xyzzz, g_0_x_zzzzzz_xzzz, g_0_x_zzzzzz_xzzzz, g_0_x_zzzzzz_yyyy, g_0_x_zzzzzz_yyyyz, g_0_x_zzzzzz_yyyz, g_0_x_zzzzzz_yyyzz, g_0_x_zzzzzz_yyzz, g_0_x_zzzzzz_yyzzz, g_0_x_zzzzzz_yzzz, g_0_x_zzzzzz_yzzzz, g_0_x_zzzzzz_zzzz, g_0_x_zzzzzz_zzzzz, g_0_x_zzzzzzz_xxxx, g_0_x_zzzzzzz_xxxy, g_0_x_zzzzzzz_xxxz, g_0_x_zzzzzzz_xxyy, g_0_x_zzzzzzz_xxyz, g_0_x_zzzzzzz_xxzz, g_0_x_zzzzzzz_xyyy, g_0_x_zzzzzzz_xyyz, g_0_x_zzzzzzz_xyzz, g_0_x_zzzzzzz_xzzz, g_0_x_zzzzzzz_yyyy, g_0_x_zzzzzzz_yyyz, g_0_x_zzzzzzz_yyzz, g_0_x_zzzzzzz_yzzz, g_0_x_zzzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_x_zzzzzzz_xxxx[k] = -g_0_x_zzzzzz_xxxx[k] * ab_z + g_0_x_zzzzzz_xxxxz[k];

                g_0_x_zzzzzzz_xxxy[k] = -g_0_x_zzzzzz_xxxy[k] * ab_z + g_0_x_zzzzzz_xxxyz[k];

                g_0_x_zzzzzzz_xxxz[k] = -g_0_x_zzzzzz_xxxz[k] * ab_z + g_0_x_zzzzzz_xxxzz[k];

                g_0_x_zzzzzzz_xxyy[k] = -g_0_x_zzzzzz_xxyy[k] * ab_z + g_0_x_zzzzzz_xxyyz[k];

                g_0_x_zzzzzzz_xxyz[k] = -g_0_x_zzzzzz_xxyz[k] * ab_z + g_0_x_zzzzzz_xxyzz[k];

                g_0_x_zzzzzzz_xxzz[k] = -g_0_x_zzzzzz_xxzz[k] * ab_z + g_0_x_zzzzzz_xxzzz[k];

                g_0_x_zzzzzzz_xyyy[k] = -g_0_x_zzzzzz_xyyy[k] * ab_z + g_0_x_zzzzzz_xyyyz[k];

                g_0_x_zzzzzzz_xyyz[k] = -g_0_x_zzzzzz_xyyz[k] * ab_z + g_0_x_zzzzzz_xyyzz[k];

                g_0_x_zzzzzzz_xyzz[k] = -g_0_x_zzzzzz_xyzz[k] * ab_z + g_0_x_zzzzzz_xyzzz[k];

                g_0_x_zzzzzzz_xzzz[k] = -g_0_x_zzzzzz_xzzz[k] * ab_z + g_0_x_zzzzzz_xzzzz[k];

                g_0_x_zzzzzzz_yyyy[k] = -g_0_x_zzzzzz_yyyy[k] * ab_z + g_0_x_zzzzzz_yyyyz[k];

                g_0_x_zzzzzzz_yyyz[k] = -g_0_x_zzzzzz_yyyz[k] * ab_z + g_0_x_zzzzzz_yyyzz[k];

                g_0_x_zzzzzzz_yyzz[k] = -g_0_x_zzzzzz_yyzz[k] * ab_z + g_0_x_zzzzzz_yyzzz[k];

                g_0_x_zzzzzzz_yzzz[k] = -g_0_x_zzzzzz_yzzz[k] * ab_z + g_0_x_zzzzzz_yzzzz[k];

                g_0_x_zzzzzzz_zzzz[k] = -g_0_x_zzzzzz_zzzz[k] * ab_z + g_0_x_zzzzzz_zzzzz[k];
            }

            /// Set up 540-555 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxx_xxxx = cbuffer.data(kg_geom_01_off + 540 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xxxy = cbuffer.data(kg_geom_01_off + 541 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xxxz = cbuffer.data(kg_geom_01_off + 542 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xxyy = cbuffer.data(kg_geom_01_off + 543 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xxyz = cbuffer.data(kg_geom_01_off + 544 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xxzz = cbuffer.data(kg_geom_01_off + 545 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xyyy = cbuffer.data(kg_geom_01_off + 546 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xyyz = cbuffer.data(kg_geom_01_off + 547 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xyzz = cbuffer.data(kg_geom_01_off + 548 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_xzzz = cbuffer.data(kg_geom_01_off + 549 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_yyyy = cbuffer.data(kg_geom_01_off + 550 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_yyyz = cbuffer.data(kg_geom_01_off + 551 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_yyzz = cbuffer.data(kg_geom_01_off + 552 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_yzzz = cbuffer.data(kg_geom_01_off + 553 * ccomps * dcomps);

            auto g_0_y_xxxxxxx_zzzz = cbuffer.data(kg_geom_01_off + 554 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxx_xxxx, g_0_y_xxxxxx_xxxxx, g_0_y_xxxxxx_xxxxy, g_0_y_xxxxxx_xxxxz, g_0_y_xxxxxx_xxxy, g_0_y_xxxxxx_xxxyy, g_0_y_xxxxxx_xxxyz, g_0_y_xxxxxx_xxxz, g_0_y_xxxxxx_xxxzz, g_0_y_xxxxxx_xxyy, g_0_y_xxxxxx_xxyyy, g_0_y_xxxxxx_xxyyz, g_0_y_xxxxxx_xxyz, g_0_y_xxxxxx_xxyzz, g_0_y_xxxxxx_xxzz, g_0_y_xxxxxx_xxzzz, g_0_y_xxxxxx_xyyy, g_0_y_xxxxxx_xyyyy, g_0_y_xxxxxx_xyyyz, g_0_y_xxxxxx_xyyz, g_0_y_xxxxxx_xyyzz, g_0_y_xxxxxx_xyzz, g_0_y_xxxxxx_xyzzz, g_0_y_xxxxxx_xzzz, g_0_y_xxxxxx_xzzzz, g_0_y_xxxxxx_yyyy, g_0_y_xxxxxx_yyyz, g_0_y_xxxxxx_yyzz, g_0_y_xxxxxx_yzzz, g_0_y_xxxxxx_zzzz, g_0_y_xxxxxxx_xxxx, g_0_y_xxxxxxx_xxxy, g_0_y_xxxxxxx_xxxz, g_0_y_xxxxxxx_xxyy, g_0_y_xxxxxxx_xxyz, g_0_y_xxxxxxx_xxzz, g_0_y_xxxxxxx_xyyy, g_0_y_xxxxxxx_xyyz, g_0_y_xxxxxxx_xyzz, g_0_y_xxxxxxx_xzzz, g_0_y_xxxxxxx_yyyy, g_0_y_xxxxxxx_yyyz, g_0_y_xxxxxxx_yyzz, g_0_y_xxxxxxx_yzzz, g_0_y_xxxxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxx_xxxx[k] = -g_0_y_xxxxxx_xxxx[k] * ab_x + g_0_y_xxxxxx_xxxxx[k];

                g_0_y_xxxxxxx_xxxy[k] = -g_0_y_xxxxxx_xxxy[k] * ab_x + g_0_y_xxxxxx_xxxxy[k];

                g_0_y_xxxxxxx_xxxz[k] = -g_0_y_xxxxxx_xxxz[k] * ab_x + g_0_y_xxxxxx_xxxxz[k];

                g_0_y_xxxxxxx_xxyy[k] = -g_0_y_xxxxxx_xxyy[k] * ab_x + g_0_y_xxxxxx_xxxyy[k];

                g_0_y_xxxxxxx_xxyz[k] = -g_0_y_xxxxxx_xxyz[k] * ab_x + g_0_y_xxxxxx_xxxyz[k];

                g_0_y_xxxxxxx_xxzz[k] = -g_0_y_xxxxxx_xxzz[k] * ab_x + g_0_y_xxxxxx_xxxzz[k];

                g_0_y_xxxxxxx_xyyy[k] = -g_0_y_xxxxxx_xyyy[k] * ab_x + g_0_y_xxxxxx_xxyyy[k];

                g_0_y_xxxxxxx_xyyz[k] = -g_0_y_xxxxxx_xyyz[k] * ab_x + g_0_y_xxxxxx_xxyyz[k];

                g_0_y_xxxxxxx_xyzz[k] = -g_0_y_xxxxxx_xyzz[k] * ab_x + g_0_y_xxxxxx_xxyzz[k];

                g_0_y_xxxxxxx_xzzz[k] = -g_0_y_xxxxxx_xzzz[k] * ab_x + g_0_y_xxxxxx_xxzzz[k];

                g_0_y_xxxxxxx_yyyy[k] = -g_0_y_xxxxxx_yyyy[k] * ab_x + g_0_y_xxxxxx_xyyyy[k];

                g_0_y_xxxxxxx_yyyz[k] = -g_0_y_xxxxxx_yyyz[k] * ab_x + g_0_y_xxxxxx_xyyyz[k];

                g_0_y_xxxxxxx_yyzz[k] = -g_0_y_xxxxxx_yyzz[k] * ab_x + g_0_y_xxxxxx_xyyzz[k];

                g_0_y_xxxxxxx_yzzz[k] = -g_0_y_xxxxxx_yzzz[k] * ab_x + g_0_y_xxxxxx_xyzzz[k];

                g_0_y_xxxxxxx_zzzz[k] = -g_0_y_xxxxxx_zzzz[k] * ab_x + g_0_y_xxxxxx_xzzzz[k];
            }

            /// Set up 555-570 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxy_xxxx = cbuffer.data(kg_geom_01_off + 555 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xxxy = cbuffer.data(kg_geom_01_off + 556 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xxxz = cbuffer.data(kg_geom_01_off + 557 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xxyy = cbuffer.data(kg_geom_01_off + 558 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xxyz = cbuffer.data(kg_geom_01_off + 559 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xxzz = cbuffer.data(kg_geom_01_off + 560 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xyyy = cbuffer.data(kg_geom_01_off + 561 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xyyz = cbuffer.data(kg_geom_01_off + 562 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xyzz = cbuffer.data(kg_geom_01_off + 563 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_xzzz = cbuffer.data(kg_geom_01_off + 564 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_yyyy = cbuffer.data(kg_geom_01_off + 565 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_yyyz = cbuffer.data(kg_geom_01_off + 566 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_yyzz = cbuffer.data(kg_geom_01_off + 567 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_yzzz = cbuffer.data(kg_geom_01_off + 568 * ccomps * dcomps);

            auto g_0_y_xxxxxxy_zzzz = cbuffer.data(kg_geom_01_off + 569 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxxy_xxxx, g_0_y_xxxxxxy_xxxy, g_0_y_xxxxxxy_xxxz, g_0_y_xxxxxxy_xxyy, g_0_y_xxxxxxy_xxyz, g_0_y_xxxxxxy_xxzz, g_0_y_xxxxxxy_xyyy, g_0_y_xxxxxxy_xyyz, g_0_y_xxxxxxy_xyzz, g_0_y_xxxxxxy_xzzz, g_0_y_xxxxxxy_yyyy, g_0_y_xxxxxxy_yyyz, g_0_y_xxxxxxy_yyzz, g_0_y_xxxxxxy_yzzz, g_0_y_xxxxxxy_zzzz, g_0_y_xxxxxy_xxxx, g_0_y_xxxxxy_xxxxx, g_0_y_xxxxxy_xxxxy, g_0_y_xxxxxy_xxxxz, g_0_y_xxxxxy_xxxy, g_0_y_xxxxxy_xxxyy, g_0_y_xxxxxy_xxxyz, g_0_y_xxxxxy_xxxz, g_0_y_xxxxxy_xxxzz, g_0_y_xxxxxy_xxyy, g_0_y_xxxxxy_xxyyy, g_0_y_xxxxxy_xxyyz, g_0_y_xxxxxy_xxyz, g_0_y_xxxxxy_xxyzz, g_0_y_xxxxxy_xxzz, g_0_y_xxxxxy_xxzzz, g_0_y_xxxxxy_xyyy, g_0_y_xxxxxy_xyyyy, g_0_y_xxxxxy_xyyyz, g_0_y_xxxxxy_xyyz, g_0_y_xxxxxy_xyyzz, g_0_y_xxxxxy_xyzz, g_0_y_xxxxxy_xyzzz, g_0_y_xxxxxy_xzzz, g_0_y_xxxxxy_xzzzz, g_0_y_xxxxxy_yyyy, g_0_y_xxxxxy_yyyz, g_0_y_xxxxxy_yyzz, g_0_y_xxxxxy_yzzz, g_0_y_xxxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxy_xxxx[k] = -g_0_y_xxxxxy_xxxx[k] * ab_x + g_0_y_xxxxxy_xxxxx[k];

                g_0_y_xxxxxxy_xxxy[k] = -g_0_y_xxxxxy_xxxy[k] * ab_x + g_0_y_xxxxxy_xxxxy[k];

                g_0_y_xxxxxxy_xxxz[k] = -g_0_y_xxxxxy_xxxz[k] * ab_x + g_0_y_xxxxxy_xxxxz[k];

                g_0_y_xxxxxxy_xxyy[k] = -g_0_y_xxxxxy_xxyy[k] * ab_x + g_0_y_xxxxxy_xxxyy[k];

                g_0_y_xxxxxxy_xxyz[k] = -g_0_y_xxxxxy_xxyz[k] * ab_x + g_0_y_xxxxxy_xxxyz[k];

                g_0_y_xxxxxxy_xxzz[k] = -g_0_y_xxxxxy_xxzz[k] * ab_x + g_0_y_xxxxxy_xxxzz[k];

                g_0_y_xxxxxxy_xyyy[k] = -g_0_y_xxxxxy_xyyy[k] * ab_x + g_0_y_xxxxxy_xxyyy[k];

                g_0_y_xxxxxxy_xyyz[k] = -g_0_y_xxxxxy_xyyz[k] * ab_x + g_0_y_xxxxxy_xxyyz[k];

                g_0_y_xxxxxxy_xyzz[k] = -g_0_y_xxxxxy_xyzz[k] * ab_x + g_0_y_xxxxxy_xxyzz[k];

                g_0_y_xxxxxxy_xzzz[k] = -g_0_y_xxxxxy_xzzz[k] * ab_x + g_0_y_xxxxxy_xxzzz[k];

                g_0_y_xxxxxxy_yyyy[k] = -g_0_y_xxxxxy_yyyy[k] * ab_x + g_0_y_xxxxxy_xyyyy[k];

                g_0_y_xxxxxxy_yyyz[k] = -g_0_y_xxxxxy_yyyz[k] * ab_x + g_0_y_xxxxxy_xyyyz[k];

                g_0_y_xxxxxxy_yyzz[k] = -g_0_y_xxxxxy_yyzz[k] * ab_x + g_0_y_xxxxxy_xyyzz[k];

                g_0_y_xxxxxxy_yzzz[k] = -g_0_y_xxxxxy_yzzz[k] * ab_x + g_0_y_xxxxxy_xyzzz[k];

                g_0_y_xxxxxxy_zzzz[k] = -g_0_y_xxxxxy_zzzz[k] * ab_x + g_0_y_xxxxxy_xzzzz[k];
            }

            /// Set up 570-585 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxxz_xxxx = cbuffer.data(kg_geom_01_off + 570 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xxxy = cbuffer.data(kg_geom_01_off + 571 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xxxz = cbuffer.data(kg_geom_01_off + 572 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xxyy = cbuffer.data(kg_geom_01_off + 573 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xxyz = cbuffer.data(kg_geom_01_off + 574 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xxzz = cbuffer.data(kg_geom_01_off + 575 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xyyy = cbuffer.data(kg_geom_01_off + 576 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xyyz = cbuffer.data(kg_geom_01_off + 577 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xyzz = cbuffer.data(kg_geom_01_off + 578 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_xzzz = cbuffer.data(kg_geom_01_off + 579 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_yyyy = cbuffer.data(kg_geom_01_off + 580 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_yyyz = cbuffer.data(kg_geom_01_off + 581 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_yyzz = cbuffer.data(kg_geom_01_off + 582 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_yzzz = cbuffer.data(kg_geom_01_off + 583 * ccomps * dcomps);

            auto g_0_y_xxxxxxz_zzzz = cbuffer.data(kg_geom_01_off + 584 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxxz_xxxx, g_0_y_xxxxxxz_xxxy, g_0_y_xxxxxxz_xxxz, g_0_y_xxxxxxz_xxyy, g_0_y_xxxxxxz_xxyz, g_0_y_xxxxxxz_xxzz, g_0_y_xxxxxxz_xyyy, g_0_y_xxxxxxz_xyyz, g_0_y_xxxxxxz_xyzz, g_0_y_xxxxxxz_xzzz, g_0_y_xxxxxxz_yyyy, g_0_y_xxxxxxz_yyyz, g_0_y_xxxxxxz_yyzz, g_0_y_xxxxxxz_yzzz, g_0_y_xxxxxxz_zzzz, g_0_y_xxxxxz_xxxx, g_0_y_xxxxxz_xxxxx, g_0_y_xxxxxz_xxxxy, g_0_y_xxxxxz_xxxxz, g_0_y_xxxxxz_xxxy, g_0_y_xxxxxz_xxxyy, g_0_y_xxxxxz_xxxyz, g_0_y_xxxxxz_xxxz, g_0_y_xxxxxz_xxxzz, g_0_y_xxxxxz_xxyy, g_0_y_xxxxxz_xxyyy, g_0_y_xxxxxz_xxyyz, g_0_y_xxxxxz_xxyz, g_0_y_xxxxxz_xxyzz, g_0_y_xxxxxz_xxzz, g_0_y_xxxxxz_xxzzz, g_0_y_xxxxxz_xyyy, g_0_y_xxxxxz_xyyyy, g_0_y_xxxxxz_xyyyz, g_0_y_xxxxxz_xyyz, g_0_y_xxxxxz_xyyzz, g_0_y_xxxxxz_xyzz, g_0_y_xxxxxz_xyzzz, g_0_y_xxxxxz_xzzz, g_0_y_xxxxxz_xzzzz, g_0_y_xxxxxz_yyyy, g_0_y_xxxxxz_yyyz, g_0_y_xxxxxz_yyzz, g_0_y_xxxxxz_yzzz, g_0_y_xxxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxxz_xxxx[k] = -g_0_y_xxxxxz_xxxx[k] * ab_x + g_0_y_xxxxxz_xxxxx[k];

                g_0_y_xxxxxxz_xxxy[k] = -g_0_y_xxxxxz_xxxy[k] * ab_x + g_0_y_xxxxxz_xxxxy[k];

                g_0_y_xxxxxxz_xxxz[k] = -g_0_y_xxxxxz_xxxz[k] * ab_x + g_0_y_xxxxxz_xxxxz[k];

                g_0_y_xxxxxxz_xxyy[k] = -g_0_y_xxxxxz_xxyy[k] * ab_x + g_0_y_xxxxxz_xxxyy[k];

                g_0_y_xxxxxxz_xxyz[k] = -g_0_y_xxxxxz_xxyz[k] * ab_x + g_0_y_xxxxxz_xxxyz[k];

                g_0_y_xxxxxxz_xxzz[k] = -g_0_y_xxxxxz_xxzz[k] * ab_x + g_0_y_xxxxxz_xxxzz[k];

                g_0_y_xxxxxxz_xyyy[k] = -g_0_y_xxxxxz_xyyy[k] * ab_x + g_0_y_xxxxxz_xxyyy[k];

                g_0_y_xxxxxxz_xyyz[k] = -g_0_y_xxxxxz_xyyz[k] * ab_x + g_0_y_xxxxxz_xxyyz[k];

                g_0_y_xxxxxxz_xyzz[k] = -g_0_y_xxxxxz_xyzz[k] * ab_x + g_0_y_xxxxxz_xxyzz[k];

                g_0_y_xxxxxxz_xzzz[k] = -g_0_y_xxxxxz_xzzz[k] * ab_x + g_0_y_xxxxxz_xxzzz[k];

                g_0_y_xxxxxxz_yyyy[k] = -g_0_y_xxxxxz_yyyy[k] * ab_x + g_0_y_xxxxxz_xyyyy[k];

                g_0_y_xxxxxxz_yyyz[k] = -g_0_y_xxxxxz_yyyz[k] * ab_x + g_0_y_xxxxxz_xyyyz[k];

                g_0_y_xxxxxxz_yyzz[k] = -g_0_y_xxxxxz_yyzz[k] * ab_x + g_0_y_xxxxxz_xyyzz[k];

                g_0_y_xxxxxxz_yzzz[k] = -g_0_y_xxxxxz_yzzz[k] * ab_x + g_0_y_xxxxxz_xyzzz[k];

                g_0_y_xxxxxxz_zzzz[k] = -g_0_y_xxxxxz_zzzz[k] * ab_x + g_0_y_xxxxxz_xzzzz[k];
            }

            /// Set up 585-600 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxyy_xxxx = cbuffer.data(kg_geom_01_off + 585 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xxxy = cbuffer.data(kg_geom_01_off + 586 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xxxz = cbuffer.data(kg_geom_01_off + 587 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xxyy = cbuffer.data(kg_geom_01_off + 588 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xxyz = cbuffer.data(kg_geom_01_off + 589 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xxzz = cbuffer.data(kg_geom_01_off + 590 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xyyy = cbuffer.data(kg_geom_01_off + 591 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xyyz = cbuffer.data(kg_geom_01_off + 592 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xyzz = cbuffer.data(kg_geom_01_off + 593 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_xzzz = cbuffer.data(kg_geom_01_off + 594 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_yyyy = cbuffer.data(kg_geom_01_off + 595 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_yyyz = cbuffer.data(kg_geom_01_off + 596 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_yyzz = cbuffer.data(kg_geom_01_off + 597 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_yzzz = cbuffer.data(kg_geom_01_off + 598 * ccomps * dcomps);

            auto g_0_y_xxxxxyy_zzzz = cbuffer.data(kg_geom_01_off + 599 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxyy_xxxx, g_0_y_xxxxxyy_xxxy, g_0_y_xxxxxyy_xxxz, g_0_y_xxxxxyy_xxyy, g_0_y_xxxxxyy_xxyz, g_0_y_xxxxxyy_xxzz, g_0_y_xxxxxyy_xyyy, g_0_y_xxxxxyy_xyyz, g_0_y_xxxxxyy_xyzz, g_0_y_xxxxxyy_xzzz, g_0_y_xxxxxyy_yyyy, g_0_y_xxxxxyy_yyyz, g_0_y_xxxxxyy_yyzz, g_0_y_xxxxxyy_yzzz, g_0_y_xxxxxyy_zzzz, g_0_y_xxxxyy_xxxx, g_0_y_xxxxyy_xxxxx, g_0_y_xxxxyy_xxxxy, g_0_y_xxxxyy_xxxxz, g_0_y_xxxxyy_xxxy, g_0_y_xxxxyy_xxxyy, g_0_y_xxxxyy_xxxyz, g_0_y_xxxxyy_xxxz, g_0_y_xxxxyy_xxxzz, g_0_y_xxxxyy_xxyy, g_0_y_xxxxyy_xxyyy, g_0_y_xxxxyy_xxyyz, g_0_y_xxxxyy_xxyz, g_0_y_xxxxyy_xxyzz, g_0_y_xxxxyy_xxzz, g_0_y_xxxxyy_xxzzz, g_0_y_xxxxyy_xyyy, g_0_y_xxxxyy_xyyyy, g_0_y_xxxxyy_xyyyz, g_0_y_xxxxyy_xyyz, g_0_y_xxxxyy_xyyzz, g_0_y_xxxxyy_xyzz, g_0_y_xxxxyy_xyzzz, g_0_y_xxxxyy_xzzz, g_0_y_xxxxyy_xzzzz, g_0_y_xxxxyy_yyyy, g_0_y_xxxxyy_yyyz, g_0_y_xxxxyy_yyzz, g_0_y_xxxxyy_yzzz, g_0_y_xxxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxyy_xxxx[k] = -g_0_y_xxxxyy_xxxx[k] * ab_x + g_0_y_xxxxyy_xxxxx[k];

                g_0_y_xxxxxyy_xxxy[k] = -g_0_y_xxxxyy_xxxy[k] * ab_x + g_0_y_xxxxyy_xxxxy[k];

                g_0_y_xxxxxyy_xxxz[k] = -g_0_y_xxxxyy_xxxz[k] * ab_x + g_0_y_xxxxyy_xxxxz[k];

                g_0_y_xxxxxyy_xxyy[k] = -g_0_y_xxxxyy_xxyy[k] * ab_x + g_0_y_xxxxyy_xxxyy[k];

                g_0_y_xxxxxyy_xxyz[k] = -g_0_y_xxxxyy_xxyz[k] * ab_x + g_0_y_xxxxyy_xxxyz[k];

                g_0_y_xxxxxyy_xxzz[k] = -g_0_y_xxxxyy_xxzz[k] * ab_x + g_0_y_xxxxyy_xxxzz[k];

                g_0_y_xxxxxyy_xyyy[k] = -g_0_y_xxxxyy_xyyy[k] * ab_x + g_0_y_xxxxyy_xxyyy[k];

                g_0_y_xxxxxyy_xyyz[k] = -g_0_y_xxxxyy_xyyz[k] * ab_x + g_0_y_xxxxyy_xxyyz[k];

                g_0_y_xxxxxyy_xyzz[k] = -g_0_y_xxxxyy_xyzz[k] * ab_x + g_0_y_xxxxyy_xxyzz[k];

                g_0_y_xxxxxyy_xzzz[k] = -g_0_y_xxxxyy_xzzz[k] * ab_x + g_0_y_xxxxyy_xxzzz[k];

                g_0_y_xxxxxyy_yyyy[k] = -g_0_y_xxxxyy_yyyy[k] * ab_x + g_0_y_xxxxyy_xyyyy[k];

                g_0_y_xxxxxyy_yyyz[k] = -g_0_y_xxxxyy_yyyz[k] * ab_x + g_0_y_xxxxyy_xyyyz[k];

                g_0_y_xxxxxyy_yyzz[k] = -g_0_y_xxxxyy_yyzz[k] * ab_x + g_0_y_xxxxyy_xyyzz[k];

                g_0_y_xxxxxyy_yzzz[k] = -g_0_y_xxxxyy_yzzz[k] * ab_x + g_0_y_xxxxyy_xyzzz[k];

                g_0_y_xxxxxyy_zzzz[k] = -g_0_y_xxxxyy_zzzz[k] * ab_x + g_0_y_xxxxyy_xzzzz[k];
            }

            /// Set up 600-615 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxyz_xxxx = cbuffer.data(kg_geom_01_off + 600 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xxxy = cbuffer.data(kg_geom_01_off + 601 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xxxz = cbuffer.data(kg_geom_01_off + 602 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xxyy = cbuffer.data(kg_geom_01_off + 603 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xxyz = cbuffer.data(kg_geom_01_off + 604 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xxzz = cbuffer.data(kg_geom_01_off + 605 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xyyy = cbuffer.data(kg_geom_01_off + 606 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xyyz = cbuffer.data(kg_geom_01_off + 607 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xyzz = cbuffer.data(kg_geom_01_off + 608 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_xzzz = cbuffer.data(kg_geom_01_off + 609 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_yyyy = cbuffer.data(kg_geom_01_off + 610 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_yyyz = cbuffer.data(kg_geom_01_off + 611 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_yyzz = cbuffer.data(kg_geom_01_off + 612 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_yzzz = cbuffer.data(kg_geom_01_off + 613 * ccomps * dcomps);

            auto g_0_y_xxxxxyz_zzzz = cbuffer.data(kg_geom_01_off + 614 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxyz_xxxx, g_0_y_xxxxxyz_xxxy, g_0_y_xxxxxyz_xxxz, g_0_y_xxxxxyz_xxyy, g_0_y_xxxxxyz_xxyz, g_0_y_xxxxxyz_xxzz, g_0_y_xxxxxyz_xyyy, g_0_y_xxxxxyz_xyyz, g_0_y_xxxxxyz_xyzz, g_0_y_xxxxxyz_xzzz, g_0_y_xxxxxyz_yyyy, g_0_y_xxxxxyz_yyyz, g_0_y_xxxxxyz_yyzz, g_0_y_xxxxxyz_yzzz, g_0_y_xxxxxyz_zzzz, g_0_y_xxxxyz_xxxx, g_0_y_xxxxyz_xxxxx, g_0_y_xxxxyz_xxxxy, g_0_y_xxxxyz_xxxxz, g_0_y_xxxxyz_xxxy, g_0_y_xxxxyz_xxxyy, g_0_y_xxxxyz_xxxyz, g_0_y_xxxxyz_xxxz, g_0_y_xxxxyz_xxxzz, g_0_y_xxxxyz_xxyy, g_0_y_xxxxyz_xxyyy, g_0_y_xxxxyz_xxyyz, g_0_y_xxxxyz_xxyz, g_0_y_xxxxyz_xxyzz, g_0_y_xxxxyz_xxzz, g_0_y_xxxxyz_xxzzz, g_0_y_xxxxyz_xyyy, g_0_y_xxxxyz_xyyyy, g_0_y_xxxxyz_xyyyz, g_0_y_xxxxyz_xyyz, g_0_y_xxxxyz_xyyzz, g_0_y_xxxxyz_xyzz, g_0_y_xxxxyz_xyzzz, g_0_y_xxxxyz_xzzz, g_0_y_xxxxyz_xzzzz, g_0_y_xxxxyz_yyyy, g_0_y_xxxxyz_yyyz, g_0_y_xxxxyz_yyzz, g_0_y_xxxxyz_yzzz, g_0_y_xxxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxyz_xxxx[k] = -g_0_y_xxxxyz_xxxx[k] * ab_x + g_0_y_xxxxyz_xxxxx[k];

                g_0_y_xxxxxyz_xxxy[k] = -g_0_y_xxxxyz_xxxy[k] * ab_x + g_0_y_xxxxyz_xxxxy[k];

                g_0_y_xxxxxyz_xxxz[k] = -g_0_y_xxxxyz_xxxz[k] * ab_x + g_0_y_xxxxyz_xxxxz[k];

                g_0_y_xxxxxyz_xxyy[k] = -g_0_y_xxxxyz_xxyy[k] * ab_x + g_0_y_xxxxyz_xxxyy[k];

                g_0_y_xxxxxyz_xxyz[k] = -g_0_y_xxxxyz_xxyz[k] * ab_x + g_0_y_xxxxyz_xxxyz[k];

                g_0_y_xxxxxyz_xxzz[k] = -g_0_y_xxxxyz_xxzz[k] * ab_x + g_0_y_xxxxyz_xxxzz[k];

                g_0_y_xxxxxyz_xyyy[k] = -g_0_y_xxxxyz_xyyy[k] * ab_x + g_0_y_xxxxyz_xxyyy[k];

                g_0_y_xxxxxyz_xyyz[k] = -g_0_y_xxxxyz_xyyz[k] * ab_x + g_0_y_xxxxyz_xxyyz[k];

                g_0_y_xxxxxyz_xyzz[k] = -g_0_y_xxxxyz_xyzz[k] * ab_x + g_0_y_xxxxyz_xxyzz[k];

                g_0_y_xxxxxyz_xzzz[k] = -g_0_y_xxxxyz_xzzz[k] * ab_x + g_0_y_xxxxyz_xxzzz[k];

                g_0_y_xxxxxyz_yyyy[k] = -g_0_y_xxxxyz_yyyy[k] * ab_x + g_0_y_xxxxyz_xyyyy[k];

                g_0_y_xxxxxyz_yyyz[k] = -g_0_y_xxxxyz_yyyz[k] * ab_x + g_0_y_xxxxyz_xyyyz[k];

                g_0_y_xxxxxyz_yyzz[k] = -g_0_y_xxxxyz_yyzz[k] * ab_x + g_0_y_xxxxyz_xyyzz[k];

                g_0_y_xxxxxyz_yzzz[k] = -g_0_y_xxxxyz_yzzz[k] * ab_x + g_0_y_xxxxyz_xyzzz[k];

                g_0_y_xxxxxyz_zzzz[k] = -g_0_y_xxxxyz_zzzz[k] * ab_x + g_0_y_xxxxyz_xzzzz[k];
            }

            /// Set up 615-630 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxxzz_xxxx = cbuffer.data(kg_geom_01_off + 615 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xxxy = cbuffer.data(kg_geom_01_off + 616 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xxxz = cbuffer.data(kg_geom_01_off + 617 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xxyy = cbuffer.data(kg_geom_01_off + 618 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xxyz = cbuffer.data(kg_geom_01_off + 619 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xxzz = cbuffer.data(kg_geom_01_off + 620 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xyyy = cbuffer.data(kg_geom_01_off + 621 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xyyz = cbuffer.data(kg_geom_01_off + 622 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xyzz = cbuffer.data(kg_geom_01_off + 623 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_xzzz = cbuffer.data(kg_geom_01_off + 624 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_yyyy = cbuffer.data(kg_geom_01_off + 625 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_yyyz = cbuffer.data(kg_geom_01_off + 626 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_yyzz = cbuffer.data(kg_geom_01_off + 627 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_yzzz = cbuffer.data(kg_geom_01_off + 628 * ccomps * dcomps);

            auto g_0_y_xxxxxzz_zzzz = cbuffer.data(kg_geom_01_off + 629 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxxzz_xxxx, g_0_y_xxxxxzz_xxxy, g_0_y_xxxxxzz_xxxz, g_0_y_xxxxxzz_xxyy, g_0_y_xxxxxzz_xxyz, g_0_y_xxxxxzz_xxzz, g_0_y_xxxxxzz_xyyy, g_0_y_xxxxxzz_xyyz, g_0_y_xxxxxzz_xyzz, g_0_y_xxxxxzz_xzzz, g_0_y_xxxxxzz_yyyy, g_0_y_xxxxxzz_yyyz, g_0_y_xxxxxzz_yyzz, g_0_y_xxxxxzz_yzzz, g_0_y_xxxxxzz_zzzz, g_0_y_xxxxzz_xxxx, g_0_y_xxxxzz_xxxxx, g_0_y_xxxxzz_xxxxy, g_0_y_xxxxzz_xxxxz, g_0_y_xxxxzz_xxxy, g_0_y_xxxxzz_xxxyy, g_0_y_xxxxzz_xxxyz, g_0_y_xxxxzz_xxxz, g_0_y_xxxxzz_xxxzz, g_0_y_xxxxzz_xxyy, g_0_y_xxxxzz_xxyyy, g_0_y_xxxxzz_xxyyz, g_0_y_xxxxzz_xxyz, g_0_y_xxxxzz_xxyzz, g_0_y_xxxxzz_xxzz, g_0_y_xxxxzz_xxzzz, g_0_y_xxxxzz_xyyy, g_0_y_xxxxzz_xyyyy, g_0_y_xxxxzz_xyyyz, g_0_y_xxxxzz_xyyz, g_0_y_xxxxzz_xyyzz, g_0_y_xxxxzz_xyzz, g_0_y_xxxxzz_xyzzz, g_0_y_xxxxzz_xzzz, g_0_y_xxxxzz_xzzzz, g_0_y_xxxxzz_yyyy, g_0_y_xxxxzz_yyyz, g_0_y_xxxxzz_yyzz, g_0_y_xxxxzz_yzzz, g_0_y_xxxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxxzz_xxxx[k] = -g_0_y_xxxxzz_xxxx[k] * ab_x + g_0_y_xxxxzz_xxxxx[k];

                g_0_y_xxxxxzz_xxxy[k] = -g_0_y_xxxxzz_xxxy[k] * ab_x + g_0_y_xxxxzz_xxxxy[k];

                g_0_y_xxxxxzz_xxxz[k] = -g_0_y_xxxxzz_xxxz[k] * ab_x + g_0_y_xxxxzz_xxxxz[k];

                g_0_y_xxxxxzz_xxyy[k] = -g_0_y_xxxxzz_xxyy[k] * ab_x + g_0_y_xxxxzz_xxxyy[k];

                g_0_y_xxxxxzz_xxyz[k] = -g_0_y_xxxxzz_xxyz[k] * ab_x + g_0_y_xxxxzz_xxxyz[k];

                g_0_y_xxxxxzz_xxzz[k] = -g_0_y_xxxxzz_xxzz[k] * ab_x + g_0_y_xxxxzz_xxxzz[k];

                g_0_y_xxxxxzz_xyyy[k] = -g_0_y_xxxxzz_xyyy[k] * ab_x + g_0_y_xxxxzz_xxyyy[k];

                g_0_y_xxxxxzz_xyyz[k] = -g_0_y_xxxxzz_xyyz[k] * ab_x + g_0_y_xxxxzz_xxyyz[k];

                g_0_y_xxxxxzz_xyzz[k] = -g_0_y_xxxxzz_xyzz[k] * ab_x + g_0_y_xxxxzz_xxyzz[k];

                g_0_y_xxxxxzz_xzzz[k] = -g_0_y_xxxxzz_xzzz[k] * ab_x + g_0_y_xxxxzz_xxzzz[k];

                g_0_y_xxxxxzz_yyyy[k] = -g_0_y_xxxxzz_yyyy[k] * ab_x + g_0_y_xxxxzz_xyyyy[k];

                g_0_y_xxxxxzz_yyyz[k] = -g_0_y_xxxxzz_yyyz[k] * ab_x + g_0_y_xxxxzz_xyyyz[k];

                g_0_y_xxxxxzz_yyzz[k] = -g_0_y_xxxxzz_yyzz[k] * ab_x + g_0_y_xxxxzz_xyyzz[k];

                g_0_y_xxxxxzz_yzzz[k] = -g_0_y_xxxxzz_yzzz[k] * ab_x + g_0_y_xxxxzz_xyzzz[k];

                g_0_y_xxxxxzz_zzzz[k] = -g_0_y_xxxxzz_zzzz[k] * ab_x + g_0_y_xxxxzz_xzzzz[k];
            }

            /// Set up 630-645 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyyy_xxxx = cbuffer.data(kg_geom_01_off + 630 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xxxy = cbuffer.data(kg_geom_01_off + 631 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xxxz = cbuffer.data(kg_geom_01_off + 632 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xxyy = cbuffer.data(kg_geom_01_off + 633 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xxyz = cbuffer.data(kg_geom_01_off + 634 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xxzz = cbuffer.data(kg_geom_01_off + 635 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xyyy = cbuffer.data(kg_geom_01_off + 636 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xyyz = cbuffer.data(kg_geom_01_off + 637 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xyzz = cbuffer.data(kg_geom_01_off + 638 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_xzzz = cbuffer.data(kg_geom_01_off + 639 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_yyyy = cbuffer.data(kg_geom_01_off + 640 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_yyyz = cbuffer.data(kg_geom_01_off + 641 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_yyzz = cbuffer.data(kg_geom_01_off + 642 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_yzzz = cbuffer.data(kg_geom_01_off + 643 * ccomps * dcomps);

            auto g_0_y_xxxxyyy_zzzz = cbuffer.data(kg_geom_01_off + 644 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyyy_xxxx, g_0_y_xxxxyyy_xxxy, g_0_y_xxxxyyy_xxxz, g_0_y_xxxxyyy_xxyy, g_0_y_xxxxyyy_xxyz, g_0_y_xxxxyyy_xxzz, g_0_y_xxxxyyy_xyyy, g_0_y_xxxxyyy_xyyz, g_0_y_xxxxyyy_xyzz, g_0_y_xxxxyyy_xzzz, g_0_y_xxxxyyy_yyyy, g_0_y_xxxxyyy_yyyz, g_0_y_xxxxyyy_yyzz, g_0_y_xxxxyyy_yzzz, g_0_y_xxxxyyy_zzzz, g_0_y_xxxyyy_xxxx, g_0_y_xxxyyy_xxxxx, g_0_y_xxxyyy_xxxxy, g_0_y_xxxyyy_xxxxz, g_0_y_xxxyyy_xxxy, g_0_y_xxxyyy_xxxyy, g_0_y_xxxyyy_xxxyz, g_0_y_xxxyyy_xxxz, g_0_y_xxxyyy_xxxzz, g_0_y_xxxyyy_xxyy, g_0_y_xxxyyy_xxyyy, g_0_y_xxxyyy_xxyyz, g_0_y_xxxyyy_xxyz, g_0_y_xxxyyy_xxyzz, g_0_y_xxxyyy_xxzz, g_0_y_xxxyyy_xxzzz, g_0_y_xxxyyy_xyyy, g_0_y_xxxyyy_xyyyy, g_0_y_xxxyyy_xyyyz, g_0_y_xxxyyy_xyyz, g_0_y_xxxyyy_xyyzz, g_0_y_xxxyyy_xyzz, g_0_y_xxxyyy_xyzzz, g_0_y_xxxyyy_xzzz, g_0_y_xxxyyy_xzzzz, g_0_y_xxxyyy_yyyy, g_0_y_xxxyyy_yyyz, g_0_y_xxxyyy_yyzz, g_0_y_xxxyyy_yzzz, g_0_y_xxxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyyy_xxxx[k] = -g_0_y_xxxyyy_xxxx[k] * ab_x + g_0_y_xxxyyy_xxxxx[k];

                g_0_y_xxxxyyy_xxxy[k] = -g_0_y_xxxyyy_xxxy[k] * ab_x + g_0_y_xxxyyy_xxxxy[k];

                g_0_y_xxxxyyy_xxxz[k] = -g_0_y_xxxyyy_xxxz[k] * ab_x + g_0_y_xxxyyy_xxxxz[k];

                g_0_y_xxxxyyy_xxyy[k] = -g_0_y_xxxyyy_xxyy[k] * ab_x + g_0_y_xxxyyy_xxxyy[k];

                g_0_y_xxxxyyy_xxyz[k] = -g_0_y_xxxyyy_xxyz[k] * ab_x + g_0_y_xxxyyy_xxxyz[k];

                g_0_y_xxxxyyy_xxzz[k] = -g_0_y_xxxyyy_xxzz[k] * ab_x + g_0_y_xxxyyy_xxxzz[k];

                g_0_y_xxxxyyy_xyyy[k] = -g_0_y_xxxyyy_xyyy[k] * ab_x + g_0_y_xxxyyy_xxyyy[k];

                g_0_y_xxxxyyy_xyyz[k] = -g_0_y_xxxyyy_xyyz[k] * ab_x + g_0_y_xxxyyy_xxyyz[k];

                g_0_y_xxxxyyy_xyzz[k] = -g_0_y_xxxyyy_xyzz[k] * ab_x + g_0_y_xxxyyy_xxyzz[k];

                g_0_y_xxxxyyy_xzzz[k] = -g_0_y_xxxyyy_xzzz[k] * ab_x + g_0_y_xxxyyy_xxzzz[k];

                g_0_y_xxxxyyy_yyyy[k] = -g_0_y_xxxyyy_yyyy[k] * ab_x + g_0_y_xxxyyy_xyyyy[k];

                g_0_y_xxxxyyy_yyyz[k] = -g_0_y_xxxyyy_yyyz[k] * ab_x + g_0_y_xxxyyy_xyyyz[k];

                g_0_y_xxxxyyy_yyzz[k] = -g_0_y_xxxyyy_yyzz[k] * ab_x + g_0_y_xxxyyy_xyyzz[k];

                g_0_y_xxxxyyy_yzzz[k] = -g_0_y_xxxyyy_yzzz[k] * ab_x + g_0_y_xxxyyy_xyzzz[k];

                g_0_y_xxxxyyy_zzzz[k] = -g_0_y_xxxyyy_zzzz[k] * ab_x + g_0_y_xxxyyy_xzzzz[k];
            }

            /// Set up 645-660 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyyz_xxxx = cbuffer.data(kg_geom_01_off + 645 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xxxy = cbuffer.data(kg_geom_01_off + 646 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xxxz = cbuffer.data(kg_geom_01_off + 647 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xxyy = cbuffer.data(kg_geom_01_off + 648 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xxyz = cbuffer.data(kg_geom_01_off + 649 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xxzz = cbuffer.data(kg_geom_01_off + 650 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xyyy = cbuffer.data(kg_geom_01_off + 651 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xyyz = cbuffer.data(kg_geom_01_off + 652 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xyzz = cbuffer.data(kg_geom_01_off + 653 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_xzzz = cbuffer.data(kg_geom_01_off + 654 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_yyyy = cbuffer.data(kg_geom_01_off + 655 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_yyyz = cbuffer.data(kg_geom_01_off + 656 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_yyzz = cbuffer.data(kg_geom_01_off + 657 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_yzzz = cbuffer.data(kg_geom_01_off + 658 * ccomps * dcomps);

            auto g_0_y_xxxxyyz_zzzz = cbuffer.data(kg_geom_01_off + 659 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyyz_xxxx, g_0_y_xxxxyyz_xxxy, g_0_y_xxxxyyz_xxxz, g_0_y_xxxxyyz_xxyy, g_0_y_xxxxyyz_xxyz, g_0_y_xxxxyyz_xxzz, g_0_y_xxxxyyz_xyyy, g_0_y_xxxxyyz_xyyz, g_0_y_xxxxyyz_xyzz, g_0_y_xxxxyyz_xzzz, g_0_y_xxxxyyz_yyyy, g_0_y_xxxxyyz_yyyz, g_0_y_xxxxyyz_yyzz, g_0_y_xxxxyyz_yzzz, g_0_y_xxxxyyz_zzzz, g_0_y_xxxyyz_xxxx, g_0_y_xxxyyz_xxxxx, g_0_y_xxxyyz_xxxxy, g_0_y_xxxyyz_xxxxz, g_0_y_xxxyyz_xxxy, g_0_y_xxxyyz_xxxyy, g_0_y_xxxyyz_xxxyz, g_0_y_xxxyyz_xxxz, g_0_y_xxxyyz_xxxzz, g_0_y_xxxyyz_xxyy, g_0_y_xxxyyz_xxyyy, g_0_y_xxxyyz_xxyyz, g_0_y_xxxyyz_xxyz, g_0_y_xxxyyz_xxyzz, g_0_y_xxxyyz_xxzz, g_0_y_xxxyyz_xxzzz, g_0_y_xxxyyz_xyyy, g_0_y_xxxyyz_xyyyy, g_0_y_xxxyyz_xyyyz, g_0_y_xxxyyz_xyyz, g_0_y_xxxyyz_xyyzz, g_0_y_xxxyyz_xyzz, g_0_y_xxxyyz_xyzzz, g_0_y_xxxyyz_xzzz, g_0_y_xxxyyz_xzzzz, g_0_y_xxxyyz_yyyy, g_0_y_xxxyyz_yyyz, g_0_y_xxxyyz_yyzz, g_0_y_xxxyyz_yzzz, g_0_y_xxxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyyz_xxxx[k] = -g_0_y_xxxyyz_xxxx[k] * ab_x + g_0_y_xxxyyz_xxxxx[k];

                g_0_y_xxxxyyz_xxxy[k] = -g_0_y_xxxyyz_xxxy[k] * ab_x + g_0_y_xxxyyz_xxxxy[k];

                g_0_y_xxxxyyz_xxxz[k] = -g_0_y_xxxyyz_xxxz[k] * ab_x + g_0_y_xxxyyz_xxxxz[k];

                g_0_y_xxxxyyz_xxyy[k] = -g_0_y_xxxyyz_xxyy[k] * ab_x + g_0_y_xxxyyz_xxxyy[k];

                g_0_y_xxxxyyz_xxyz[k] = -g_0_y_xxxyyz_xxyz[k] * ab_x + g_0_y_xxxyyz_xxxyz[k];

                g_0_y_xxxxyyz_xxzz[k] = -g_0_y_xxxyyz_xxzz[k] * ab_x + g_0_y_xxxyyz_xxxzz[k];

                g_0_y_xxxxyyz_xyyy[k] = -g_0_y_xxxyyz_xyyy[k] * ab_x + g_0_y_xxxyyz_xxyyy[k];

                g_0_y_xxxxyyz_xyyz[k] = -g_0_y_xxxyyz_xyyz[k] * ab_x + g_0_y_xxxyyz_xxyyz[k];

                g_0_y_xxxxyyz_xyzz[k] = -g_0_y_xxxyyz_xyzz[k] * ab_x + g_0_y_xxxyyz_xxyzz[k];

                g_0_y_xxxxyyz_xzzz[k] = -g_0_y_xxxyyz_xzzz[k] * ab_x + g_0_y_xxxyyz_xxzzz[k];

                g_0_y_xxxxyyz_yyyy[k] = -g_0_y_xxxyyz_yyyy[k] * ab_x + g_0_y_xxxyyz_xyyyy[k];

                g_0_y_xxxxyyz_yyyz[k] = -g_0_y_xxxyyz_yyyz[k] * ab_x + g_0_y_xxxyyz_xyyyz[k];

                g_0_y_xxxxyyz_yyzz[k] = -g_0_y_xxxyyz_yyzz[k] * ab_x + g_0_y_xxxyyz_xyyzz[k];

                g_0_y_xxxxyyz_yzzz[k] = -g_0_y_xxxyyz_yzzz[k] * ab_x + g_0_y_xxxyyz_xyzzz[k];

                g_0_y_xxxxyyz_zzzz[k] = -g_0_y_xxxyyz_zzzz[k] * ab_x + g_0_y_xxxyyz_xzzzz[k];
            }

            /// Set up 660-675 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxyzz_xxxx = cbuffer.data(kg_geom_01_off + 660 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xxxy = cbuffer.data(kg_geom_01_off + 661 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xxxz = cbuffer.data(kg_geom_01_off + 662 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xxyy = cbuffer.data(kg_geom_01_off + 663 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xxyz = cbuffer.data(kg_geom_01_off + 664 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xxzz = cbuffer.data(kg_geom_01_off + 665 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xyyy = cbuffer.data(kg_geom_01_off + 666 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xyyz = cbuffer.data(kg_geom_01_off + 667 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xyzz = cbuffer.data(kg_geom_01_off + 668 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_xzzz = cbuffer.data(kg_geom_01_off + 669 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_yyyy = cbuffer.data(kg_geom_01_off + 670 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_yyyz = cbuffer.data(kg_geom_01_off + 671 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_yyzz = cbuffer.data(kg_geom_01_off + 672 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_yzzz = cbuffer.data(kg_geom_01_off + 673 * ccomps * dcomps);

            auto g_0_y_xxxxyzz_zzzz = cbuffer.data(kg_geom_01_off + 674 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxyzz_xxxx, g_0_y_xxxxyzz_xxxy, g_0_y_xxxxyzz_xxxz, g_0_y_xxxxyzz_xxyy, g_0_y_xxxxyzz_xxyz, g_0_y_xxxxyzz_xxzz, g_0_y_xxxxyzz_xyyy, g_0_y_xxxxyzz_xyyz, g_0_y_xxxxyzz_xyzz, g_0_y_xxxxyzz_xzzz, g_0_y_xxxxyzz_yyyy, g_0_y_xxxxyzz_yyyz, g_0_y_xxxxyzz_yyzz, g_0_y_xxxxyzz_yzzz, g_0_y_xxxxyzz_zzzz, g_0_y_xxxyzz_xxxx, g_0_y_xxxyzz_xxxxx, g_0_y_xxxyzz_xxxxy, g_0_y_xxxyzz_xxxxz, g_0_y_xxxyzz_xxxy, g_0_y_xxxyzz_xxxyy, g_0_y_xxxyzz_xxxyz, g_0_y_xxxyzz_xxxz, g_0_y_xxxyzz_xxxzz, g_0_y_xxxyzz_xxyy, g_0_y_xxxyzz_xxyyy, g_0_y_xxxyzz_xxyyz, g_0_y_xxxyzz_xxyz, g_0_y_xxxyzz_xxyzz, g_0_y_xxxyzz_xxzz, g_0_y_xxxyzz_xxzzz, g_0_y_xxxyzz_xyyy, g_0_y_xxxyzz_xyyyy, g_0_y_xxxyzz_xyyyz, g_0_y_xxxyzz_xyyz, g_0_y_xxxyzz_xyyzz, g_0_y_xxxyzz_xyzz, g_0_y_xxxyzz_xyzzz, g_0_y_xxxyzz_xzzz, g_0_y_xxxyzz_xzzzz, g_0_y_xxxyzz_yyyy, g_0_y_xxxyzz_yyyz, g_0_y_xxxyzz_yyzz, g_0_y_xxxyzz_yzzz, g_0_y_xxxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxyzz_xxxx[k] = -g_0_y_xxxyzz_xxxx[k] * ab_x + g_0_y_xxxyzz_xxxxx[k];

                g_0_y_xxxxyzz_xxxy[k] = -g_0_y_xxxyzz_xxxy[k] * ab_x + g_0_y_xxxyzz_xxxxy[k];

                g_0_y_xxxxyzz_xxxz[k] = -g_0_y_xxxyzz_xxxz[k] * ab_x + g_0_y_xxxyzz_xxxxz[k];

                g_0_y_xxxxyzz_xxyy[k] = -g_0_y_xxxyzz_xxyy[k] * ab_x + g_0_y_xxxyzz_xxxyy[k];

                g_0_y_xxxxyzz_xxyz[k] = -g_0_y_xxxyzz_xxyz[k] * ab_x + g_0_y_xxxyzz_xxxyz[k];

                g_0_y_xxxxyzz_xxzz[k] = -g_0_y_xxxyzz_xxzz[k] * ab_x + g_0_y_xxxyzz_xxxzz[k];

                g_0_y_xxxxyzz_xyyy[k] = -g_0_y_xxxyzz_xyyy[k] * ab_x + g_0_y_xxxyzz_xxyyy[k];

                g_0_y_xxxxyzz_xyyz[k] = -g_0_y_xxxyzz_xyyz[k] * ab_x + g_0_y_xxxyzz_xxyyz[k];

                g_0_y_xxxxyzz_xyzz[k] = -g_0_y_xxxyzz_xyzz[k] * ab_x + g_0_y_xxxyzz_xxyzz[k];

                g_0_y_xxxxyzz_xzzz[k] = -g_0_y_xxxyzz_xzzz[k] * ab_x + g_0_y_xxxyzz_xxzzz[k];

                g_0_y_xxxxyzz_yyyy[k] = -g_0_y_xxxyzz_yyyy[k] * ab_x + g_0_y_xxxyzz_xyyyy[k];

                g_0_y_xxxxyzz_yyyz[k] = -g_0_y_xxxyzz_yyyz[k] * ab_x + g_0_y_xxxyzz_xyyyz[k];

                g_0_y_xxxxyzz_yyzz[k] = -g_0_y_xxxyzz_yyzz[k] * ab_x + g_0_y_xxxyzz_xyyzz[k];

                g_0_y_xxxxyzz_yzzz[k] = -g_0_y_xxxyzz_yzzz[k] * ab_x + g_0_y_xxxyzz_xyzzz[k];

                g_0_y_xxxxyzz_zzzz[k] = -g_0_y_xxxyzz_zzzz[k] * ab_x + g_0_y_xxxyzz_xzzzz[k];
            }

            /// Set up 675-690 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxxzzz_xxxx = cbuffer.data(kg_geom_01_off + 675 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xxxy = cbuffer.data(kg_geom_01_off + 676 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xxxz = cbuffer.data(kg_geom_01_off + 677 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xxyy = cbuffer.data(kg_geom_01_off + 678 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xxyz = cbuffer.data(kg_geom_01_off + 679 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xxzz = cbuffer.data(kg_geom_01_off + 680 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xyyy = cbuffer.data(kg_geom_01_off + 681 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xyyz = cbuffer.data(kg_geom_01_off + 682 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xyzz = cbuffer.data(kg_geom_01_off + 683 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_xzzz = cbuffer.data(kg_geom_01_off + 684 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_yyyy = cbuffer.data(kg_geom_01_off + 685 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_yyyz = cbuffer.data(kg_geom_01_off + 686 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_yyzz = cbuffer.data(kg_geom_01_off + 687 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_yzzz = cbuffer.data(kg_geom_01_off + 688 * ccomps * dcomps);

            auto g_0_y_xxxxzzz_zzzz = cbuffer.data(kg_geom_01_off + 689 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxxzzz_xxxx, g_0_y_xxxxzzz_xxxy, g_0_y_xxxxzzz_xxxz, g_0_y_xxxxzzz_xxyy, g_0_y_xxxxzzz_xxyz, g_0_y_xxxxzzz_xxzz, g_0_y_xxxxzzz_xyyy, g_0_y_xxxxzzz_xyyz, g_0_y_xxxxzzz_xyzz, g_0_y_xxxxzzz_xzzz, g_0_y_xxxxzzz_yyyy, g_0_y_xxxxzzz_yyyz, g_0_y_xxxxzzz_yyzz, g_0_y_xxxxzzz_yzzz, g_0_y_xxxxzzz_zzzz, g_0_y_xxxzzz_xxxx, g_0_y_xxxzzz_xxxxx, g_0_y_xxxzzz_xxxxy, g_0_y_xxxzzz_xxxxz, g_0_y_xxxzzz_xxxy, g_0_y_xxxzzz_xxxyy, g_0_y_xxxzzz_xxxyz, g_0_y_xxxzzz_xxxz, g_0_y_xxxzzz_xxxzz, g_0_y_xxxzzz_xxyy, g_0_y_xxxzzz_xxyyy, g_0_y_xxxzzz_xxyyz, g_0_y_xxxzzz_xxyz, g_0_y_xxxzzz_xxyzz, g_0_y_xxxzzz_xxzz, g_0_y_xxxzzz_xxzzz, g_0_y_xxxzzz_xyyy, g_0_y_xxxzzz_xyyyy, g_0_y_xxxzzz_xyyyz, g_0_y_xxxzzz_xyyz, g_0_y_xxxzzz_xyyzz, g_0_y_xxxzzz_xyzz, g_0_y_xxxzzz_xyzzz, g_0_y_xxxzzz_xzzz, g_0_y_xxxzzz_xzzzz, g_0_y_xxxzzz_yyyy, g_0_y_xxxzzz_yyyz, g_0_y_xxxzzz_yyzz, g_0_y_xxxzzz_yzzz, g_0_y_xxxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxxzzz_xxxx[k] = -g_0_y_xxxzzz_xxxx[k] * ab_x + g_0_y_xxxzzz_xxxxx[k];

                g_0_y_xxxxzzz_xxxy[k] = -g_0_y_xxxzzz_xxxy[k] * ab_x + g_0_y_xxxzzz_xxxxy[k];

                g_0_y_xxxxzzz_xxxz[k] = -g_0_y_xxxzzz_xxxz[k] * ab_x + g_0_y_xxxzzz_xxxxz[k];

                g_0_y_xxxxzzz_xxyy[k] = -g_0_y_xxxzzz_xxyy[k] * ab_x + g_0_y_xxxzzz_xxxyy[k];

                g_0_y_xxxxzzz_xxyz[k] = -g_0_y_xxxzzz_xxyz[k] * ab_x + g_0_y_xxxzzz_xxxyz[k];

                g_0_y_xxxxzzz_xxzz[k] = -g_0_y_xxxzzz_xxzz[k] * ab_x + g_0_y_xxxzzz_xxxzz[k];

                g_0_y_xxxxzzz_xyyy[k] = -g_0_y_xxxzzz_xyyy[k] * ab_x + g_0_y_xxxzzz_xxyyy[k];

                g_0_y_xxxxzzz_xyyz[k] = -g_0_y_xxxzzz_xyyz[k] * ab_x + g_0_y_xxxzzz_xxyyz[k];

                g_0_y_xxxxzzz_xyzz[k] = -g_0_y_xxxzzz_xyzz[k] * ab_x + g_0_y_xxxzzz_xxyzz[k];

                g_0_y_xxxxzzz_xzzz[k] = -g_0_y_xxxzzz_xzzz[k] * ab_x + g_0_y_xxxzzz_xxzzz[k];

                g_0_y_xxxxzzz_yyyy[k] = -g_0_y_xxxzzz_yyyy[k] * ab_x + g_0_y_xxxzzz_xyyyy[k];

                g_0_y_xxxxzzz_yyyz[k] = -g_0_y_xxxzzz_yyyz[k] * ab_x + g_0_y_xxxzzz_xyyyz[k];

                g_0_y_xxxxzzz_yyzz[k] = -g_0_y_xxxzzz_yyzz[k] * ab_x + g_0_y_xxxzzz_xyyzz[k];

                g_0_y_xxxxzzz_yzzz[k] = -g_0_y_xxxzzz_yzzz[k] * ab_x + g_0_y_xxxzzz_xyzzz[k];

                g_0_y_xxxxzzz_zzzz[k] = -g_0_y_xxxzzz_zzzz[k] * ab_x + g_0_y_xxxzzz_xzzzz[k];
            }

            /// Set up 690-705 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyyy_xxxx = cbuffer.data(kg_geom_01_off + 690 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xxxy = cbuffer.data(kg_geom_01_off + 691 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xxxz = cbuffer.data(kg_geom_01_off + 692 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xxyy = cbuffer.data(kg_geom_01_off + 693 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xxyz = cbuffer.data(kg_geom_01_off + 694 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xxzz = cbuffer.data(kg_geom_01_off + 695 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xyyy = cbuffer.data(kg_geom_01_off + 696 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xyyz = cbuffer.data(kg_geom_01_off + 697 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xyzz = cbuffer.data(kg_geom_01_off + 698 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_xzzz = cbuffer.data(kg_geom_01_off + 699 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_yyyy = cbuffer.data(kg_geom_01_off + 700 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_yyyz = cbuffer.data(kg_geom_01_off + 701 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_yyzz = cbuffer.data(kg_geom_01_off + 702 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_yzzz = cbuffer.data(kg_geom_01_off + 703 * ccomps * dcomps);

            auto g_0_y_xxxyyyy_zzzz = cbuffer.data(kg_geom_01_off + 704 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyyy_xxxx, g_0_y_xxxyyyy_xxxy, g_0_y_xxxyyyy_xxxz, g_0_y_xxxyyyy_xxyy, g_0_y_xxxyyyy_xxyz, g_0_y_xxxyyyy_xxzz, g_0_y_xxxyyyy_xyyy, g_0_y_xxxyyyy_xyyz, g_0_y_xxxyyyy_xyzz, g_0_y_xxxyyyy_xzzz, g_0_y_xxxyyyy_yyyy, g_0_y_xxxyyyy_yyyz, g_0_y_xxxyyyy_yyzz, g_0_y_xxxyyyy_yzzz, g_0_y_xxxyyyy_zzzz, g_0_y_xxyyyy_xxxx, g_0_y_xxyyyy_xxxxx, g_0_y_xxyyyy_xxxxy, g_0_y_xxyyyy_xxxxz, g_0_y_xxyyyy_xxxy, g_0_y_xxyyyy_xxxyy, g_0_y_xxyyyy_xxxyz, g_0_y_xxyyyy_xxxz, g_0_y_xxyyyy_xxxzz, g_0_y_xxyyyy_xxyy, g_0_y_xxyyyy_xxyyy, g_0_y_xxyyyy_xxyyz, g_0_y_xxyyyy_xxyz, g_0_y_xxyyyy_xxyzz, g_0_y_xxyyyy_xxzz, g_0_y_xxyyyy_xxzzz, g_0_y_xxyyyy_xyyy, g_0_y_xxyyyy_xyyyy, g_0_y_xxyyyy_xyyyz, g_0_y_xxyyyy_xyyz, g_0_y_xxyyyy_xyyzz, g_0_y_xxyyyy_xyzz, g_0_y_xxyyyy_xyzzz, g_0_y_xxyyyy_xzzz, g_0_y_xxyyyy_xzzzz, g_0_y_xxyyyy_yyyy, g_0_y_xxyyyy_yyyz, g_0_y_xxyyyy_yyzz, g_0_y_xxyyyy_yzzz, g_0_y_xxyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyyy_xxxx[k] = -g_0_y_xxyyyy_xxxx[k] * ab_x + g_0_y_xxyyyy_xxxxx[k];

                g_0_y_xxxyyyy_xxxy[k] = -g_0_y_xxyyyy_xxxy[k] * ab_x + g_0_y_xxyyyy_xxxxy[k];

                g_0_y_xxxyyyy_xxxz[k] = -g_0_y_xxyyyy_xxxz[k] * ab_x + g_0_y_xxyyyy_xxxxz[k];

                g_0_y_xxxyyyy_xxyy[k] = -g_0_y_xxyyyy_xxyy[k] * ab_x + g_0_y_xxyyyy_xxxyy[k];

                g_0_y_xxxyyyy_xxyz[k] = -g_0_y_xxyyyy_xxyz[k] * ab_x + g_0_y_xxyyyy_xxxyz[k];

                g_0_y_xxxyyyy_xxzz[k] = -g_0_y_xxyyyy_xxzz[k] * ab_x + g_0_y_xxyyyy_xxxzz[k];

                g_0_y_xxxyyyy_xyyy[k] = -g_0_y_xxyyyy_xyyy[k] * ab_x + g_0_y_xxyyyy_xxyyy[k];

                g_0_y_xxxyyyy_xyyz[k] = -g_0_y_xxyyyy_xyyz[k] * ab_x + g_0_y_xxyyyy_xxyyz[k];

                g_0_y_xxxyyyy_xyzz[k] = -g_0_y_xxyyyy_xyzz[k] * ab_x + g_0_y_xxyyyy_xxyzz[k];

                g_0_y_xxxyyyy_xzzz[k] = -g_0_y_xxyyyy_xzzz[k] * ab_x + g_0_y_xxyyyy_xxzzz[k];

                g_0_y_xxxyyyy_yyyy[k] = -g_0_y_xxyyyy_yyyy[k] * ab_x + g_0_y_xxyyyy_xyyyy[k];

                g_0_y_xxxyyyy_yyyz[k] = -g_0_y_xxyyyy_yyyz[k] * ab_x + g_0_y_xxyyyy_xyyyz[k];

                g_0_y_xxxyyyy_yyzz[k] = -g_0_y_xxyyyy_yyzz[k] * ab_x + g_0_y_xxyyyy_xyyzz[k];

                g_0_y_xxxyyyy_yzzz[k] = -g_0_y_xxyyyy_yzzz[k] * ab_x + g_0_y_xxyyyy_xyzzz[k];

                g_0_y_xxxyyyy_zzzz[k] = -g_0_y_xxyyyy_zzzz[k] * ab_x + g_0_y_xxyyyy_xzzzz[k];
            }

            /// Set up 705-720 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyyz_xxxx = cbuffer.data(kg_geom_01_off + 705 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xxxy = cbuffer.data(kg_geom_01_off + 706 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xxxz = cbuffer.data(kg_geom_01_off + 707 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xxyy = cbuffer.data(kg_geom_01_off + 708 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xxyz = cbuffer.data(kg_geom_01_off + 709 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xxzz = cbuffer.data(kg_geom_01_off + 710 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xyyy = cbuffer.data(kg_geom_01_off + 711 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xyyz = cbuffer.data(kg_geom_01_off + 712 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xyzz = cbuffer.data(kg_geom_01_off + 713 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_xzzz = cbuffer.data(kg_geom_01_off + 714 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_yyyy = cbuffer.data(kg_geom_01_off + 715 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_yyyz = cbuffer.data(kg_geom_01_off + 716 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_yyzz = cbuffer.data(kg_geom_01_off + 717 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_yzzz = cbuffer.data(kg_geom_01_off + 718 * ccomps * dcomps);

            auto g_0_y_xxxyyyz_zzzz = cbuffer.data(kg_geom_01_off + 719 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyyz_xxxx, g_0_y_xxxyyyz_xxxy, g_0_y_xxxyyyz_xxxz, g_0_y_xxxyyyz_xxyy, g_0_y_xxxyyyz_xxyz, g_0_y_xxxyyyz_xxzz, g_0_y_xxxyyyz_xyyy, g_0_y_xxxyyyz_xyyz, g_0_y_xxxyyyz_xyzz, g_0_y_xxxyyyz_xzzz, g_0_y_xxxyyyz_yyyy, g_0_y_xxxyyyz_yyyz, g_0_y_xxxyyyz_yyzz, g_0_y_xxxyyyz_yzzz, g_0_y_xxxyyyz_zzzz, g_0_y_xxyyyz_xxxx, g_0_y_xxyyyz_xxxxx, g_0_y_xxyyyz_xxxxy, g_0_y_xxyyyz_xxxxz, g_0_y_xxyyyz_xxxy, g_0_y_xxyyyz_xxxyy, g_0_y_xxyyyz_xxxyz, g_0_y_xxyyyz_xxxz, g_0_y_xxyyyz_xxxzz, g_0_y_xxyyyz_xxyy, g_0_y_xxyyyz_xxyyy, g_0_y_xxyyyz_xxyyz, g_0_y_xxyyyz_xxyz, g_0_y_xxyyyz_xxyzz, g_0_y_xxyyyz_xxzz, g_0_y_xxyyyz_xxzzz, g_0_y_xxyyyz_xyyy, g_0_y_xxyyyz_xyyyy, g_0_y_xxyyyz_xyyyz, g_0_y_xxyyyz_xyyz, g_0_y_xxyyyz_xyyzz, g_0_y_xxyyyz_xyzz, g_0_y_xxyyyz_xyzzz, g_0_y_xxyyyz_xzzz, g_0_y_xxyyyz_xzzzz, g_0_y_xxyyyz_yyyy, g_0_y_xxyyyz_yyyz, g_0_y_xxyyyz_yyzz, g_0_y_xxyyyz_yzzz, g_0_y_xxyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyyz_xxxx[k] = -g_0_y_xxyyyz_xxxx[k] * ab_x + g_0_y_xxyyyz_xxxxx[k];

                g_0_y_xxxyyyz_xxxy[k] = -g_0_y_xxyyyz_xxxy[k] * ab_x + g_0_y_xxyyyz_xxxxy[k];

                g_0_y_xxxyyyz_xxxz[k] = -g_0_y_xxyyyz_xxxz[k] * ab_x + g_0_y_xxyyyz_xxxxz[k];

                g_0_y_xxxyyyz_xxyy[k] = -g_0_y_xxyyyz_xxyy[k] * ab_x + g_0_y_xxyyyz_xxxyy[k];

                g_0_y_xxxyyyz_xxyz[k] = -g_0_y_xxyyyz_xxyz[k] * ab_x + g_0_y_xxyyyz_xxxyz[k];

                g_0_y_xxxyyyz_xxzz[k] = -g_0_y_xxyyyz_xxzz[k] * ab_x + g_0_y_xxyyyz_xxxzz[k];

                g_0_y_xxxyyyz_xyyy[k] = -g_0_y_xxyyyz_xyyy[k] * ab_x + g_0_y_xxyyyz_xxyyy[k];

                g_0_y_xxxyyyz_xyyz[k] = -g_0_y_xxyyyz_xyyz[k] * ab_x + g_0_y_xxyyyz_xxyyz[k];

                g_0_y_xxxyyyz_xyzz[k] = -g_0_y_xxyyyz_xyzz[k] * ab_x + g_0_y_xxyyyz_xxyzz[k];

                g_0_y_xxxyyyz_xzzz[k] = -g_0_y_xxyyyz_xzzz[k] * ab_x + g_0_y_xxyyyz_xxzzz[k];

                g_0_y_xxxyyyz_yyyy[k] = -g_0_y_xxyyyz_yyyy[k] * ab_x + g_0_y_xxyyyz_xyyyy[k];

                g_0_y_xxxyyyz_yyyz[k] = -g_0_y_xxyyyz_yyyz[k] * ab_x + g_0_y_xxyyyz_xyyyz[k];

                g_0_y_xxxyyyz_yyzz[k] = -g_0_y_xxyyyz_yyzz[k] * ab_x + g_0_y_xxyyyz_xyyzz[k];

                g_0_y_xxxyyyz_yzzz[k] = -g_0_y_xxyyyz_yzzz[k] * ab_x + g_0_y_xxyyyz_xyzzz[k];

                g_0_y_xxxyyyz_zzzz[k] = -g_0_y_xxyyyz_zzzz[k] * ab_x + g_0_y_xxyyyz_xzzzz[k];
            }

            /// Set up 720-735 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyyzz_xxxx = cbuffer.data(kg_geom_01_off + 720 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xxxy = cbuffer.data(kg_geom_01_off + 721 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xxxz = cbuffer.data(kg_geom_01_off + 722 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xxyy = cbuffer.data(kg_geom_01_off + 723 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xxyz = cbuffer.data(kg_geom_01_off + 724 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xxzz = cbuffer.data(kg_geom_01_off + 725 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xyyy = cbuffer.data(kg_geom_01_off + 726 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xyyz = cbuffer.data(kg_geom_01_off + 727 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xyzz = cbuffer.data(kg_geom_01_off + 728 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_xzzz = cbuffer.data(kg_geom_01_off + 729 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_yyyy = cbuffer.data(kg_geom_01_off + 730 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_yyyz = cbuffer.data(kg_geom_01_off + 731 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_yyzz = cbuffer.data(kg_geom_01_off + 732 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_yzzz = cbuffer.data(kg_geom_01_off + 733 * ccomps * dcomps);

            auto g_0_y_xxxyyzz_zzzz = cbuffer.data(kg_geom_01_off + 734 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyyzz_xxxx, g_0_y_xxxyyzz_xxxy, g_0_y_xxxyyzz_xxxz, g_0_y_xxxyyzz_xxyy, g_0_y_xxxyyzz_xxyz, g_0_y_xxxyyzz_xxzz, g_0_y_xxxyyzz_xyyy, g_0_y_xxxyyzz_xyyz, g_0_y_xxxyyzz_xyzz, g_0_y_xxxyyzz_xzzz, g_0_y_xxxyyzz_yyyy, g_0_y_xxxyyzz_yyyz, g_0_y_xxxyyzz_yyzz, g_0_y_xxxyyzz_yzzz, g_0_y_xxxyyzz_zzzz, g_0_y_xxyyzz_xxxx, g_0_y_xxyyzz_xxxxx, g_0_y_xxyyzz_xxxxy, g_0_y_xxyyzz_xxxxz, g_0_y_xxyyzz_xxxy, g_0_y_xxyyzz_xxxyy, g_0_y_xxyyzz_xxxyz, g_0_y_xxyyzz_xxxz, g_0_y_xxyyzz_xxxzz, g_0_y_xxyyzz_xxyy, g_0_y_xxyyzz_xxyyy, g_0_y_xxyyzz_xxyyz, g_0_y_xxyyzz_xxyz, g_0_y_xxyyzz_xxyzz, g_0_y_xxyyzz_xxzz, g_0_y_xxyyzz_xxzzz, g_0_y_xxyyzz_xyyy, g_0_y_xxyyzz_xyyyy, g_0_y_xxyyzz_xyyyz, g_0_y_xxyyzz_xyyz, g_0_y_xxyyzz_xyyzz, g_0_y_xxyyzz_xyzz, g_0_y_xxyyzz_xyzzz, g_0_y_xxyyzz_xzzz, g_0_y_xxyyzz_xzzzz, g_0_y_xxyyzz_yyyy, g_0_y_xxyyzz_yyyz, g_0_y_xxyyzz_yyzz, g_0_y_xxyyzz_yzzz, g_0_y_xxyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyyzz_xxxx[k] = -g_0_y_xxyyzz_xxxx[k] * ab_x + g_0_y_xxyyzz_xxxxx[k];

                g_0_y_xxxyyzz_xxxy[k] = -g_0_y_xxyyzz_xxxy[k] * ab_x + g_0_y_xxyyzz_xxxxy[k];

                g_0_y_xxxyyzz_xxxz[k] = -g_0_y_xxyyzz_xxxz[k] * ab_x + g_0_y_xxyyzz_xxxxz[k];

                g_0_y_xxxyyzz_xxyy[k] = -g_0_y_xxyyzz_xxyy[k] * ab_x + g_0_y_xxyyzz_xxxyy[k];

                g_0_y_xxxyyzz_xxyz[k] = -g_0_y_xxyyzz_xxyz[k] * ab_x + g_0_y_xxyyzz_xxxyz[k];

                g_0_y_xxxyyzz_xxzz[k] = -g_0_y_xxyyzz_xxzz[k] * ab_x + g_0_y_xxyyzz_xxxzz[k];

                g_0_y_xxxyyzz_xyyy[k] = -g_0_y_xxyyzz_xyyy[k] * ab_x + g_0_y_xxyyzz_xxyyy[k];

                g_0_y_xxxyyzz_xyyz[k] = -g_0_y_xxyyzz_xyyz[k] * ab_x + g_0_y_xxyyzz_xxyyz[k];

                g_0_y_xxxyyzz_xyzz[k] = -g_0_y_xxyyzz_xyzz[k] * ab_x + g_0_y_xxyyzz_xxyzz[k];

                g_0_y_xxxyyzz_xzzz[k] = -g_0_y_xxyyzz_xzzz[k] * ab_x + g_0_y_xxyyzz_xxzzz[k];

                g_0_y_xxxyyzz_yyyy[k] = -g_0_y_xxyyzz_yyyy[k] * ab_x + g_0_y_xxyyzz_xyyyy[k];

                g_0_y_xxxyyzz_yyyz[k] = -g_0_y_xxyyzz_yyyz[k] * ab_x + g_0_y_xxyyzz_xyyyz[k];

                g_0_y_xxxyyzz_yyzz[k] = -g_0_y_xxyyzz_yyzz[k] * ab_x + g_0_y_xxyyzz_xyyzz[k];

                g_0_y_xxxyyzz_yzzz[k] = -g_0_y_xxyyzz_yzzz[k] * ab_x + g_0_y_xxyyzz_xyzzz[k];

                g_0_y_xxxyyzz_zzzz[k] = -g_0_y_xxyyzz_zzzz[k] * ab_x + g_0_y_xxyyzz_xzzzz[k];
            }

            /// Set up 735-750 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxyzzz_xxxx = cbuffer.data(kg_geom_01_off + 735 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xxxy = cbuffer.data(kg_geom_01_off + 736 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xxxz = cbuffer.data(kg_geom_01_off + 737 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xxyy = cbuffer.data(kg_geom_01_off + 738 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xxyz = cbuffer.data(kg_geom_01_off + 739 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xxzz = cbuffer.data(kg_geom_01_off + 740 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xyyy = cbuffer.data(kg_geom_01_off + 741 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xyyz = cbuffer.data(kg_geom_01_off + 742 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xyzz = cbuffer.data(kg_geom_01_off + 743 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_xzzz = cbuffer.data(kg_geom_01_off + 744 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_yyyy = cbuffer.data(kg_geom_01_off + 745 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_yyyz = cbuffer.data(kg_geom_01_off + 746 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_yyzz = cbuffer.data(kg_geom_01_off + 747 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_yzzz = cbuffer.data(kg_geom_01_off + 748 * ccomps * dcomps);

            auto g_0_y_xxxyzzz_zzzz = cbuffer.data(kg_geom_01_off + 749 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxyzzz_xxxx, g_0_y_xxxyzzz_xxxy, g_0_y_xxxyzzz_xxxz, g_0_y_xxxyzzz_xxyy, g_0_y_xxxyzzz_xxyz, g_0_y_xxxyzzz_xxzz, g_0_y_xxxyzzz_xyyy, g_0_y_xxxyzzz_xyyz, g_0_y_xxxyzzz_xyzz, g_0_y_xxxyzzz_xzzz, g_0_y_xxxyzzz_yyyy, g_0_y_xxxyzzz_yyyz, g_0_y_xxxyzzz_yyzz, g_0_y_xxxyzzz_yzzz, g_0_y_xxxyzzz_zzzz, g_0_y_xxyzzz_xxxx, g_0_y_xxyzzz_xxxxx, g_0_y_xxyzzz_xxxxy, g_0_y_xxyzzz_xxxxz, g_0_y_xxyzzz_xxxy, g_0_y_xxyzzz_xxxyy, g_0_y_xxyzzz_xxxyz, g_0_y_xxyzzz_xxxz, g_0_y_xxyzzz_xxxzz, g_0_y_xxyzzz_xxyy, g_0_y_xxyzzz_xxyyy, g_0_y_xxyzzz_xxyyz, g_0_y_xxyzzz_xxyz, g_0_y_xxyzzz_xxyzz, g_0_y_xxyzzz_xxzz, g_0_y_xxyzzz_xxzzz, g_0_y_xxyzzz_xyyy, g_0_y_xxyzzz_xyyyy, g_0_y_xxyzzz_xyyyz, g_0_y_xxyzzz_xyyz, g_0_y_xxyzzz_xyyzz, g_0_y_xxyzzz_xyzz, g_0_y_xxyzzz_xyzzz, g_0_y_xxyzzz_xzzz, g_0_y_xxyzzz_xzzzz, g_0_y_xxyzzz_yyyy, g_0_y_xxyzzz_yyyz, g_0_y_xxyzzz_yyzz, g_0_y_xxyzzz_yzzz, g_0_y_xxyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxyzzz_xxxx[k] = -g_0_y_xxyzzz_xxxx[k] * ab_x + g_0_y_xxyzzz_xxxxx[k];

                g_0_y_xxxyzzz_xxxy[k] = -g_0_y_xxyzzz_xxxy[k] * ab_x + g_0_y_xxyzzz_xxxxy[k];

                g_0_y_xxxyzzz_xxxz[k] = -g_0_y_xxyzzz_xxxz[k] * ab_x + g_0_y_xxyzzz_xxxxz[k];

                g_0_y_xxxyzzz_xxyy[k] = -g_0_y_xxyzzz_xxyy[k] * ab_x + g_0_y_xxyzzz_xxxyy[k];

                g_0_y_xxxyzzz_xxyz[k] = -g_0_y_xxyzzz_xxyz[k] * ab_x + g_0_y_xxyzzz_xxxyz[k];

                g_0_y_xxxyzzz_xxzz[k] = -g_0_y_xxyzzz_xxzz[k] * ab_x + g_0_y_xxyzzz_xxxzz[k];

                g_0_y_xxxyzzz_xyyy[k] = -g_0_y_xxyzzz_xyyy[k] * ab_x + g_0_y_xxyzzz_xxyyy[k];

                g_0_y_xxxyzzz_xyyz[k] = -g_0_y_xxyzzz_xyyz[k] * ab_x + g_0_y_xxyzzz_xxyyz[k];

                g_0_y_xxxyzzz_xyzz[k] = -g_0_y_xxyzzz_xyzz[k] * ab_x + g_0_y_xxyzzz_xxyzz[k];

                g_0_y_xxxyzzz_xzzz[k] = -g_0_y_xxyzzz_xzzz[k] * ab_x + g_0_y_xxyzzz_xxzzz[k];

                g_0_y_xxxyzzz_yyyy[k] = -g_0_y_xxyzzz_yyyy[k] * ab_x + g_0_y_xxyzzz_xyyyy[k];

                g_0_y_xxxyzzz_yyyz[k] = -g_0_y_xxyzzz_yyyz[k] * ab_x + g_0_y_xxyzzz_xyyyz[k];

                g_0_y_xxxyzzz_yyzz[k] = -g_0_y_xxyzzz_yyzz[k] * ab_x + g_0_y_xxyzzz_xyyzz[k];

                g_0_y_xxxyzzz_yzzz[k] = -g_0_y_xxyzzz_yzzz[k] * ab_x + g_0_y_xxyzzz_xyzzz[k];

                g_0_y_xxxyzzz_zzzz[k] = -g_0_y_xxyzzz_zzzz[k] * ab_x + g_0_y_xxyzzz_xzzzz[k];
            }

            /// Set up 750-765 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxxzzzz_xxxx = cbuffer.data(kg_geom_01_off + 750 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xxxy = cbuffer.data(kg_geom_01_off + 751 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xxxz = cbuffer.data(kg_geom_01_off + 752 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xxyy = cbuffer.data(kg_geom_01_off + 753 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xxyz = cbuffer.data(kg_geom_01_off + 754 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xxzz = cbuffer.data(kg_geom_01_off + 755 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xyyy = cbuffer.data(kg_geom_01_off + 756 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xyyz = cbuffer.data(kg_geom_01_off + 757 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xyzz = cbuffer.data(kg_geom_01_off + 758 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_xzzz = cbuffer.data(kg_geom_01_off + 759 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_yyyy = cbuffer.data(kg_geom_01_off + 760 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_yyyz = cbuffer.data(kg_geom_01_off + 761 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_yyzz = cbuffer.data(kg_geom_01_off + 762 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_yzzz = cbuffer.data(kg_geom_01_off + 763 * ccomps * dcomps);

            auto g_0_y_xxxzzzz_zzzz = cbuffer.data(kg_geom_01_off + 764 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxxzzzz_xxxx, g_0_y_xxxzzzz_xxxy, g_0_y_xxxzzzz_xxxz, g_0_y_xxxzzzz_xxyy, g_0_y_xxxzzzz_xxyz, g_0_y_xxxzzzz_xxzz, g_0_y_xxxzzzz_xyyy, g_0_y_xxxzzzz_xyyz, g_0_y_xxxzzzz_xyzz, g_0_y_xxxzzzz_xzzz, g_0_y_xxxzzzz_yyyy, g_0_y_xxxzzzz_yyyz, g_0_y_xxxzzzz_yyzz, g_0_y_xxxzzzz_yzzz, g_0_y_xxxzzzz_zzzz, g_0_y_xxzzzz_xxxx, g_0_y_xxzzzz_xxxxx, g_0_y_xxzzzz_xxxxy, g_0_y_xxzzzz_xxxxz, g_0_y_xxzzzz_xxxy, g_0_y_xxzzzz_xxxyy, g_0_y_xxzzzz_xxxyz, g_0_y_xxzzzz_xxxz, g_0_y_xxzzzz_xxxzz, g_0_y_xxzzzz_xxyy, g_0_y_xxzzzz_xxyyy, g_0_y_xxzzzz_xxyyz, g_0_y_xxzzzz_xxyz, g_0_y_xxzzzz_xxyzz, g_0_y_xxzzzz_xxzz, g_0_y_xxzzzz_xxzzz, g_0_y_xxzzzz_xyyy, g_0_y_xxzzzz_xyyyy, g_0_y_xxzzzz_xyyyz, g_0_y_xxzzzz_xyyz, g_0_y_xxzzzz_xyyzz, g_0_y_xxzzzz_xyzz, g_0_y_xxzzzz_xyzzz, g_0_y_xxzzzz_xzzz, g_0_y_xxzzzz_xzzzz, g_0_y_xxzzzz_yyyy, g_0_y_xxzzzz_yyyz, g_0_y_xxzzzz_yyzz, g_0_y_xxzzzz_yzzz, g_0_y_xxzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxxzzzz_xxxx[k] = -g_0_y_xxzzzz_xxxx[k] * ab_x + g_0_y_xxzzzz_xxxxx[k];

                g_0_y_xxxzzzz_xxxy[k] = -g_0_y_xxzzzz_xxxy[k] * ab_x + g_0_y_xxzzzz_xxxxy[k];

                g_0_y_xxxzzzz_xxxz[k] = -g_0_y_xxzzzz_xxxz[k] * ab_x + g_0_y_xxzzzz_xxxxz[k];

                g_0_y_xxxzzzz_xxyy[k] = -g_0_y_xxzzzz_xxyy[k] * ab_x + g_0_y_xxzzzz_xxxyy[k];

                g_0_y_xxxzzzz_xxyz[k] = -g_0_y_xxzzzz_xxyz[k] * ab_x + g_0_y_xxzzzz_xxxyz[k];

                g_0_y_xxxzzzz_xxzz[k] = -g_0_y_xxzzzz_xxzz[k] * ab_x + g_0_y_xxzzzz_xxxzz[k];

                g_0_y_xxxzzzz_xyyy[k] = -g_0_y_xxzzzz_xyyy[k] * ab_x + g_0_y_xxzzzz_xxyyy[k];

                g_0_y_xxxzzzz_xyyz[k] = -g_0_y_xxzzzz_xyyz[k] * ab_x + g_0_y_xxzzzz_xxyyz[k];

                g_0_y_xxxzzzz_xyzz[k] = -g_0_y_xxzzzz_xyzz[k] * ab_x + g_0_y_xxzzzz_xxyzz[k];

                g_0_y_xxxzzzz_xzzz[k] = -g_0_y_xxzzzz_xzzz[k] * ab_x + g_0_y_xxzzzz_xxzzz[k];

                g_0_y_xxxzzzz_yyyy[k] = -g_0_y_xxzzzz_yyyy[k] * ab_x + g_0_y_xxzzzz_xyyyy[k];

                g_0_y_xxxzzzz_yyyz[k] = -g_0_y_xxzzzz_yyyz[k] * ab_x + g_0_y_xxzzzz_xyyyz[k];

                g_0_y_xxxzzzz_yyzz[k] = -g_0_y_xxzzzz_yyzz[k] * ab_x + g_0_y_xxzzzz_xyyzz[k];

                g_0_y_xxxzzzz_yzzz[k] = -g_0_y_xxzzzz_yzzz[k] * ab_x + g_0_y_xxzzzz_xyzzz[k];

                g_0_y_xxxzzzz_zzzz[k] = -g_0_y_xxzzzz_zzzz[k] * ab_x + g_0_y_xxzzzz_xzzzz[k];
            }

            /// Set up 765-780 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyyy_xxxx = cbuffer.data(kg_geom_01_off + 765 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xxxy = cbuffer.data(kg_geom_01_off + 766 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xxxz = cbuffer.data(kg_geom_01_off + 767 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xxyy = cbuffer.data(kg_geom_01_off + 768 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xxyz = cbuffer.data(kg_geom_01_off + 769 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xxzz = cbuffer.data(kg_geom_01_off + 770 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xyyy = cbuffer.data(kg_geom_01_off + 771 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xyyz = cbuffer.data(kg_geom_01_off + 772 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xyzz = cbuffer.data(kg_geom_01_off + 773 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_xzzz = cbuffer.data(kg_geom_01_off + 774 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_yyyy = cbuffer.data(kg_geom_01_off + 775 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_yyyz = cbuffer.data(kg_geom_01_off + 776 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_yyzz = cbuffer.data(kg_geom_01_off + 777 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_yzzz = cbuffer.data(kg_geom_01_off + 778 * ccomps * dcomps);

            auto g_0_y_xxyyyyy_zzzz = cbuffer.data(kg_geom_01_off + 779 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyyy_xxxx, g_0_y_xxyyyyy_xxxy, g_0_y_xxyyyyy_xxxz, g_0_y_xxyyyyy_xxyy, g_0_y_xxyyyyy_xxyz, g_0_y_xxyyyyy_xxzz, g_0_y_xxyyyyy_xyyy, g_0_y_xxyyyyy_xyyz, g_0_y_xxyyyyy_xyzz, g_0_y_xxyyyyy_xzzz, g_0_y_xxyyyyy_yyyy, g_0_y_xxyyyyy_yyyz, g_0_y_xxyyyyy_yyzz, g_0_y_xxyyyyy_yzzz, g_0_y_xxyyyyy_zzzz, g_0_y_xyyyyy_xxxx, g_0_y_xyyyyy_xxxxx, g_0_y_xyyyyy_xxxxy, g_0_y_xyyyyy_xxxxz, g_0_y_xyyyyy_xxxy, g_0_y_xyyyyy_xxxyy, g_0_y_xyyyyy_xxxyz, g_0_y_xyyyyy_xxxz, g_0_y_xyyyyy_xxxzz, g_0_y_xyyyyy_xxyy, g_0_y_xyyyyy_xxyyy, g_0_y_xyyyyy_xxyyz, g_0_y_xyyyyy_xxyz, g_0_y_xyyyyy_xxyzz, g_0_y_xyyyyy_xxzz, g_0_y_xyyyyy_xxzzz, g_0_y_xyyyyy_xyyy, g_0_y_xyyyyy_xyyyy, g_0_y_xyyyyy_xyyyz, g_0_y_xyyyyy_xyyz, g_0_y_xyyyyy_xyyzz, g_0_y_xyyyyy_xyzz, g_0_y_xyyyyy_xyzzz, g_0_y_xyyyyy_xzzz, g_0_y_xyyyyy_xzzzz, g_0_y_xyyyyy_yyyy, g_0_y_xyyyyy_yyyz, g_0_y_xyyyyy_yyzz, g_0_y_xyyyyy_yzzz, g_0_y_xyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyyy_xxxx[k] = -g_0_y_xyyyyy_xxxx[k] * ab_x + g_0_y_xyyyyy_xxxxx[k];

                g_0_y_xxyyyyy_xxxy[k] = -g_0_y_xyyyyy_xxxy[k] * ab_x + g_0_y_xyyyyy_xxxxy[k];

                g_0_y_xxyyyyy_xxxz[k] = -g_0_y_xyyyyy_xxxz[k] * ab_x + g_0_y_xyyyyy_xxxxz[k];

                g_0_y_xxyyyyy_xxyy[k] = -g_0_y_xyyyyy_xxyy[k] * ab_x + g_0_y_xyyyyy_xxxyy[k];

                g_0_y_xxyyyyy_xxyz[k] = -g_0_y_xyyyyy_xxyz[k] * ab_x + g_0_y_xyyyyy_xxxyz[k];

                g_0_y_xxyyyyy_xxzz[k] = -g_0_y_xyyyyy_xxzz[k] * ab_x + g_0_y_xyyyyy_xxxzz[k];

                g_0_y_xxyyyyy_xyyy[k] = -g_0_y_xyyyyy_xyyy[k] * ab_x + g_0_y_xyyyyy_xxyyy[k];

                g_0_y_xxyyyyy_xyyz[k] = -g_0_y_xyyyyy_xyyz[k] * ab_x + g_0_y_xyyyyy_xxyyz[k];

                g_0_y_xxyyyyy_xyzz[k] = -g_0_y_xyyyyy_xyzz[k] * ab_x + g_0_y_xyyyyy_xxyzz[k];

                g_0_y_xxyyyyy_xzzz[k] = -g_0_y_xyyyyy_xzzz[k] * ab_x + g_0_y_xyyyyy_xxzzz[k];

                g_0_y_xxyyyyy_yyyy[k] = -g_0_y_xyyyyy_yyyy[k] * ab_x + g_0_y_xyyyyy_xyyyy[k];

                g_0_y_xxyyyyy_yyyz[k] = -g_0_y_xyyyyy_yyyz[k] * ab_x + g_0_y_xyyyyy_xyyyz[k];

                g_0_y_xxyyyyy_yyzz[k] = -g_0_y_xyyyyy_yyzz[k] * ab_x + g_0_y_xyyyyy_xyyzz[k];

                g_0_y_xxyyyyy_yzzz[k] = -g_0_y_xyyyyy_yzzz[k] * ab_x + g_0_y_xyyyyy_xyzzz[k];

                g_0_y_xxyyyyy_zzzz[k] = -g_0_y_xyyyyy_zzzz[k] * ab_x + g_0_y_xyyyyy_xzzzz[k];
            }

            /// Set up 780-795 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyyz_xxxx = cbuffer.data(kg_geom_01_off + 780 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xxxy = cbuffer.data(kg_geom_01_off + 781 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xxxz = cbuffer.data(kg_geom_01_off + 782 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xxyy = cbuffer.data(kg_geom_01_off + 783 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xxyz = cbuffer.data(kg_geom_01_off + 784 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xxzz = cbuffer.data(kg_geom_01_off + 785 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xyyy = cbuffer.data(kg_geom_01_off + 786 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xyyz = cbuffer.data(kg_geom_01_off + 787 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xyzz = cbuffer.data(kg_geom_01_off + 788 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_xzzz = cbuffer.data(kg_geom_01_off + 789 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_yyyy = cbuffer.data(kg_geom_01_off + 790 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_yyyz = cbuffer.data(kg_geom_01_off + 791 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_yyzz = cbuffer.data(kg_geom_01_off + 792 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_yzzz = cbuffer.data(kg_geom_01_off + 793 * ccomps * dcomps);

            auto g_0_y_xxyyyyz_zzzz = cbuffer.data(kg_geom_01_off + 794 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyyz_xxxx, g_0_y_xxyyyyz_xxxy, g_0_y_xxyyyyz_xxxz, g_0_y_xxyyyyz_xxyy, g_0_y_xxyyyyz_xxyz, g_0_y_xxyyyyz_xxzz, g_0_y_xxyyyyz_xyyy, g_0_y_xxyyyyz_xyyz, g_0_y_xxyyyyz_xyzz, g_0_y_xxyyyyz_xzzz, g_0_y_xxyyyyz_yyyy, g_0_y_xxyyyyz_yyyz, g_0_y_xxyyyyz_yyzz, g_0_y_xxyyyyz_yzzz, g_0_y_xxyyyyz_zzzz, g_0_y_xyyyyz_xxxx, g_0_y_xyyyyz_xxxxx, g_0_y_xyyyyz_xxxxy, g_0_y_xyyyyz_xxxxz, g_0_y_xyyyyz_xxxy, g_0_y_xyyyyz_xxxyy, g_0_y_xyyyyz_xxxyz, g_0_y_xyyyyz_xxxz, g_0_y_xyyyyz_xxxzz, g_0_y_xyyyyz_xxyy, g_0_y_xyyyyz_xxyyy, g_0_y_xyyyyz_xxyyz, g_0_y_xyyyyz_xxyz, g_0_y_xyyyyz_xxyzz, g_0_y_xyyyyz_xxzz, g_0_y_xyyyyz_xxzzz, g_0_y_xyyyyz_xyyy, g_0_y_xyyyyz_xyyyy, g_0_y_xyyyyz_xyyyz, g_0_y_xyyyyz_xyyz, g_0_y_xyyyyz_xyyzz, g_0_y_xyyyyz_xyzz, g_0_y_xyyyyz_xyzzz, g_0_y_xyyyyz_xzzz, g_0_y_xyyyyz_xzzzz, g_0_y_xyyyyz_yyyy, g_0_y_xyyyyz_yyyz, g_0_y_xyyyyz_yyzz, g_0_y_xyyyyz_yzzz, g_0_y_xyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyyz_xxxx[k] = -g_0_y_xyyyyz_xxxx[k] * ab_x + g_0_y_xyyyyz_xxxxx[k];

                g_0_y_xxyyyyz_xxxy[k] = -g_0_y_xyyyyz_xxxy[k] * ab_x + g_0_y_xyyyyz_xxxxy[k];

                g_0_y_xxyyyyz_xxxz[k] = -g_0_y_xyyyyz_xxxz[k] * ab_x + g_0_y_xyyyyz_xxxxz[k];

                g_0_y_xxyyyyz_xxyy[k] = -g_0_y_xyyyyz_xxyy[k] * ab_x + g_0_y_xyyyyz_xxxyy[k];

                g_0_y_xxyyyyz_xxyz[k] = -g_0_y_xyyyyz_xxyz[k] * ab_x + g_0_y_xyyyyz_xxxyz[k];

                g_0_y_xxyyyyz_xxzz[k] = -g_0_y_xyyyyz_xxzz[k] * ab_x + g_0_y_xyyyyz_xxxzz[k];

                g_0_y_xxyyyyz_xyyy[k] = -g_0_y_xyyyyz_xyyy[k] * ab_x + g_0_y_xyyyyz_xxyyy[k];

                g_0_y_xxyyyyz_xyyz[k] = -g_0_y_xyyyyz_xyyz[k] * ab_x + g_0_y_xyyyyz_xxyyz[k];

                g_0_y_xxyyyyz_xyzz[k] = -g_0_y_xyyyyz_xyzz[k] * ab_x + g_0_y_xyyyyz_xxyzz[k];

                g_0_y_xxyyyyz_xzzz[k] = -g_0_y_xyyyyz_xzzz[k] * ab_x + g_0_y_xyyyyz_xxzzz[k];

                g_0_y_xxyyyyz_yyyy[k] = -g_0_y_xyyyyz_yyyy[k] * ab_x + g_0_y_xyyyyz_xyyyy[k];

                g_0_y_xxyyyyz_yyyz[k] = -g_0_y_xyyyyz_yyyz[k] * ab_x + g_0_y_xyyyyz_xyyyz[k];

                g_0_y_xxyyyyz_yyzz[k] = -g_0_y_xyyyyz_yyzz[k] * ab_x + g_0_y_xyyyyz_xyyzz[k];

                g_0_y_xxyyyyz_yzzz[k] = -g_0_y_xyyyyz_yzzz[k] * ab_x + g_0_y_xyyyyz_xyzzz[k];

                g_0_y_xxyyyyz_zzzz[k] = -g_0_y_xyyyyz_zzzz[k] * ab_x + g_0_y_xyyyyz_xzzzz[k];
            }

            /// Set up 795-810 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyyzz_xxxx = cbuffer.data(kg_geom_01_off + 795 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xxxy = cbuffer.data(kg_geom_01_off + 796 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xxxz = cbuffer.data(kg_geom_01_off + 797 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xxyy = cbuffer.data(kg_geom_01_off + 798 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xxyz = cbuffer.data(kg_geom_01_off + 799 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xxzz = cbuffer.data(kg_geom_01_off + 800 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xyyy = cbuffer.data(kg_geom_01_off + 801 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xyyz = cbuffer.data(kg_geom_01_off + 802 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xyzz = cbuffer.data(kg_geom_01_off + 803 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_xzzz = cbuffer.data(kg_geom_01_off + 804 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_yyyy = cbuffer.data(kg_geom_01_off + 805 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_yyyz = cbuffer.data(kg_geom_01_off + 806 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_yyzz = cbuffer.data(kg_geom_01_off + 807 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_yzzz = cbuffer.data(kg_geom_01_off + 808 * ccomps * dcomps);

            auto g_0_y_xxyyyzz_zzzz = cbuffer.data(kg_geom_01_off + 809 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyyzz_xxxx, g_0_y_xxyyyzz_xxxy, g_0_y_xxyyyzz_xxxz, g_0_y_xxyyyzz_xxyy, g_0_y_xxyyyzz_xxyz, g_0_y_xxyyyzz_xxzz, g_0_y_xxyyyzz_xyyy, g_0_y_xxyyyzz_xyyz, g_0_y_xxyyyzz_xyzz, g_0_y_xxyyyzz_xzzz, g_0_y_xxyyyzz_yyyy, g_0_y_xxyyyzz_yyyz, g_0_y_xxyyyzz_yyzz, g_0_y_xxyyyzz_yzzz, g_0_y_xxyyyzz_zzzz, g_0_y_xyyyzz_xxxx, g_0_y_xyyyzz_xxxxx, g_0_y_xyyyzz_xxxxy, g_0_y_xyyyzz_xxxxz, g_0_y_xyyyzz_xxxy, g_0_y_xyyyzz_xxxyy, g_0_y_xyyyzz_xxxyz, g_0_y_xyyyzz_xxxz, g_0_y_xyyyzz_xxxzz, g_0_y_xyyyzz_xxyy, g_0_y_xyyyzz_xxyyy, g_0_y_xyyyzz_xxyyz, g_0_y_xyyyzz_xxyz, g_0_y_xyyyzz_xxyzz, g_0_y_xyyyzz_xxzz, g_0_y_xyyyzz_xxzzz, g_0_y_xyyyzz_xyyy, g_0_y_xyyyzz_xyyyy, g_0_y_xyyyzz_xyyyz, g_0_y_xyyyzz_xyyz, g_0_y_xyyyzz_xyyzz, g_0_y_xyyyzz_xyzz, g_0_y_xyyyzz_xyzzz, g_0_y_xyyyzz_xzzz, g_0_y_xyyyzz_xzzzz, g_0_y_xyyyzz_yyyy, g_0_y_xyyyzz_yyyz, g_0_y_xyyyzz_yyzz, g_0_y_xyyyzz_yzzz, g_0_y_xyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyyzz_xxxx[k] = -g_0_y_xyyyzz_xxxx[k] * ab_x + g_0_y_xyyyzz_xxxxx[k];

                g_0_y_xxyyyzz_xxxy[k] = -g_0_y_xyyyzz_xxxy[k] * ab_x + g_0_y_xyyyzz_xxxxy[k];

                g_0_y_xxyyyzz_xxxz[k] = -g_0_y_xyyyzz_xxxz[k] * ab_x + g_0_y_xyyyzz_xxxxz[k];

                g_0_y_xxyyyzz_xxyy[k] = -g_0_y_xyyyzz_xxyy[k] * ab_x + g_0_y_xyyyzz_xxxyy[k];

                g_0_y_xxyyyzz_xxyz[k] = -g_0_y_xyyyzz_xxyz[k] * ab_x + g_0_y_xyyyzz_xxxyz[k];

                g_0_y_xxyyyzz_xxzz[k] = -g_0_y_xyyyzz_xxzz[k] * ab_x + g_0_y_xyyyzz_xxxzz[k];

                g_0_y_xxyyyzz_xyyy[k] = -g_0_y_xyyyzz_xyyy[k] * ab_x + g_0_y_xyyyzz_xxyyy[k];

                g_0_y_xxyyyzz_xyyz[k] = -g_0_y_xyyyzz_xyyz[k] * ab_x + g_0_y_xyyyzz_xxyyz[k];

                g_0_y_xxyyyzz_xyzz[k] = -g_0_y_xyyyzz_xyzz[k] * ab_x + g_0_y_xyyyzz_xxyzz[k];

                g_0_y_xxyyyzz_xzzz[k] = -g_0_y_xyyyzz_xzzz[k] * ab_x + g_0_y_xyyyzz_xxzzz[k];

                g_0_y_xxyyyzz_yyyy[k] = -g_0_y_xyyyzz_yyyy[k] * ab_x + g_0_y_xyyyzz_xyyyy[k];

                g_0_y_xxyyyzz_yyyz[k] = -g_0_y_xyyyzz_yyyz[k] * ab_x + g_0_y_xyyyzz_xyyyz[k];

                g_0_y_xxyyyzz_yyzz[k] = -g_0_y_xyyyzz_yyzz[k] * ab_x + g_0_y_xyyyzz_xyyzz[k];

                g_0_y_xxyyyzz_yzzz[k] = -g_0_y_xyyyzz_yzzz[k] * ab_x + g_0_y_xyyyzz_xyzzz[k];

                g_0_y_xxyyyzz_zzzz[k] = -g_0_y_xyyyzz_zzzz[k] * ab_x + g_0_y_xyyyzz_xzzzz[k];
            }

            /// Set up 810-825 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyyzzz_xxxx = cbuffer.data(kg_geom_01_off + 810 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xxxy = cbuffer.data(kg_geom_01_off + 811 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xxxz = cbuffer.data(kg_geom_01_off + 812 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xxyy = cbuffer.data(kg_geom_01_off + 813 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xxyz = cbuffer.data(kg_geom_01_off + 814 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xxzz = cbuffer.data(kg_geom_01_off + 815 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xyyy = cbuffer.data(kg_geom_01_off + 816 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xyyz = cbuffer.data(kg_geom_01_off + 817 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xyzz = cbuffer.data(kg_geom_01_off + 818 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_xzzz = cbuffer.data(kg_geom_01_off + 819 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_yyyy = cbuffer.data(kg_geom_01_off + 820 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_yyyz = cbuffer.data(kg_geom_01_off + 821 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_yyzz = cbuffer.data(kg_geom_01_off + 822 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_yzzz = cbuffer.data(kg_geom_01_off + 823 * ccomps * dcomps);

            auto g_0_y_xxyyzzz_zzzz = cbuffer.data(kg_geom_01_off + 824 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyyzzz_xxxx, g_0_y_xxyyzzz_xxxy, g_0_y_xxyyzzz_xxxz, g_0_y_xxyyzzz_xxyy, g_0_y_xxyyzzz_xxyz, g_0_y_xxyyzzz_xxzz, g_0_y_xxyyzzz_xyyy, g_0_y_xxyyzzz_xyyz, g_0_y_xxyyzzz_xyzz, g_0_y_xxyyzzz_xzzz, g_0_y_xxyyzzz_yyyy, g_0_y_xxyyzzz_yyyz, g_0_y_xxyyzzz_yyzz, g_0_y_xxyyzzz_yzzz, g_0_y_xxyyzzz_zzzz, g_0_y_xyyzzz_xxxx, g_0_y_xyyzzz_xxxxx, g_0_y_xyyzzz_xxxxy, g_0_y_xyyzzz_xxxxz, g_0_y_xyyzzz_xxxy, g_0_y_xyyzzz_xxxyy, g_0_y_xyyzzz_xxxyz, g_0_y_xyyzzz_xxxz, g_0_y_xyyzzz_xxxzz, g_0_y_xyyzzz_xxyy, g_0_y_xyyzzz_xxyyy, g_0_y_xyyzzz_xxyyz, g_0_y_xyyzzz_xxyz, g_0_y_xyyzzz_xxyzz, g_0_y_xyyzzz_xxzz, g_0_y_xyyzzz_xxzzz, g_0_y_xyyzzz_xyyy, g_0_y_xyyzzz_xyyyy, g_0_y_xyyzzz_xyyyz, g_0_y_xyyzzz_xyyz, g_0_y_xyyzzz_xyyzz, g_0_y_xyyzzz_xyzz, g_0_y_xyyzzz_xyzzz, g_0_y_xyyzzz_xzzz, g_0_y_xyyzzz_xzzzz, g_0_y_xyyzzz_yyyy, g_0_y_xyyzzz_yyyz, g_0_y_xyyzzz_yyzz, g_0_y_xyyzzz_yzzz, g_0_y_xyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyyzzz_xxxx[k] = -g_0_y_xyyzzz_xxxx[k] * ab_x + g_0_y_xyyzzz_xxxxx[k];

                g_0_y_xxyyzzz_xxxy[k] = -g_0_y_xyyzzz_xxxy[k] * ab_x + g_0_y_xyyzzz_xxxxy[k];

                g_0_y_xxyyzzz_xxxz[k] = -g_0_y_xyyzzz_xxxz[k] * ab_x + g_0_y_xyyzzz_xxxxz[k];

                g_0_y_xxyyzzz_xxyy[k] = -g_0_y_xyyzzz_xxyy[k] * ab_x + g_0_y_xyyzzz_xxxyy[k];

                g_0_y_xxyyzzz_xxyz[k] = -g_0_y_xyyzzz_xxyz[k] * ab_x + g_0_y_xyyzzz_xxxyz[k];

                g_0_y_xxyyzzz_xxzz[k] = -g_0_y_xyyzzz_xxzz[k] * ab_x + g_0_y_xyyzzz_xxxzz[k];

                g_0_y_xxyyzzz_xyyy[k] = -g_0_y_xyyzzz_xyyy[k] * ab_x + g_0_y_xyyzzz_xxyyy[k];

                g_0_y_xxyyzzz_xyyz[k] = -g_0_y_xyyzzz_xyyz[k] * ab_x + g_0_y_xyyzzz_xxyyz[k];

                g_0_y_xxyyzzz_xyzz[k] = -g_0_y_xyyzzz_xyzz[k] * ab_x + g_0_y_xyyzzz_xxyzz[k];

                g_0_y_xxyyzzz_xzzz[k] = -g_0_y_xyyzzz_xzzz[k] * ab_x + g_0_y_xyyzzz_xxzzz[k];

                g_0_y_xxyyzzz_yyyy[k] = -g_0_y_xyyzzz_yyyy[k] * ab_x + g_0_y_xyyzzz_xyyyy[k];

                g_0_y_xxyyzzz_yyyz[k] = -g_0_y_xyyzzz_yyyz[k] * ab_x + g_0_y_xyyzzz_xyyyz[k];

                g_0_y_xxyyzzz_yyzz[k] = -g_0_y_xyyzzz_yyzz[k] * ab_x + g_0_y_xyyzzz_xyyzz[k];

                g_0_y_xxyyzzz_yzzz[k] = -g_0_y_xyyzzz_yzzz[k] * ab_x + g_0_y_xyyzzz_xyzzz[k];

                g_0_y_xxyyzzz_zzzz[k] = -g_0_y_xyyzzz_zzzz[k] * ab_x + g_0_y_xyyzzz_xzzzz[k];
            }

            /// Set up 825-840 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxyzzzz_xxxx = cbuffer.data(kg_geom_01_off + 825 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xxxy = cbuffer.data(kg_geom_01_off + 826 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xxxz = cbuffer.data(kg_geom_01_off + 827 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xxyy = cbuffer.data(kg_geom_01_off + 828 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xxyz = cbuffer.data(kg_geom_01_off + 829 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xxzz = cbuffer.data(kg_geom_01_off + 830 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xyyy = cbuffer.data(kg_geom_01_off + 831 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xyyz = cbuffer.data(kg_geom_01_off + 832 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xyzz = cbuffer.data(kg_geom_01_off + 833 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_xzzz = cbuffer.data(kg_geom_01_off + 834 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_yyyy = cbuffer.data(kg_geom_01_off + 835 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_yyyz = cbuffer.data(kg_geom_01_off + 836 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_yyzz = cbuffer.data(kg_geom_01_off + 837 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_yzzz = cbuffer.data(kg_geom_01_off + 838 * ccomps * dcomps);

            auto g_0_y_xxyzzzz_zzzz = cbuffer.data(kg_geom_01_off + 839 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxyzzzz_xxxx, g_0_y_xxyzzzz_xxxy, g_0_y_xxyzzzz_xxxz, g_0_y_xxyzzzz_xxyy, g_0_y_xxyzzzz_xxyz, g_0_y_xxyzzzz_xxzz, g_0_y_xxyzzzz_xyyy, g_0_y_xxyzzzz_xyyz, g_0_y_xxyzzzz_xyzz, g_0_y_xxyzzzz_xzzz, g_0_y_xxyzzzz_yyyy, g_0_y_xxyzzzz_yyyz, g_0_y_xxyzzzz_yyzz, g_0_y_xxyzzzz_yzzz, g_0_y_xxyzzzz_zzzz, g_0_y_xyzzzz_xxxx, g_0_y_xyzzzz_xxxxx, g_0_y_xyzzzz_xxxxy, g_0_y_xyzzzz_xxxxz, g_0_y_xyzzzz_xxxy, g_0_y_xyzzzz_xxxyy, g_0_y_xyzzzz_xxxyz, g_0_y_xyzzzz_xxxz, g_0_y_xyzzzz_xxxzz, g_0_y_xyzzzz_xxyy, g_0_y_xyzzzz_xxyyy, g_0_y_xyzzzz_xxyyz, g_0_y_xyzzzz_xxyz, g_0_y_xyzzzz_xxyzz, g_0_y_xyzzzz_xxzz, g_0_y_xyzzzz_xxzzz, g_0_y_xyzzzz_xyyy, g_0_y_xyzzzz_xyyyy, g_0_y_xyzzzz_xyyyz, g_0_y_xyzzzz_xyyz, g_0_y_xyzzzz_xyyzz, g_0_y_xyzzzz_xyzz, g_0_y_xyzzzz_xyzzz, g_0_y_xyzzzz_xzzz, g_0_y_xyzzzz_xzzzz, g_0_y_xyzzzz_yyyy, g_0_y_xyzzzz_yyyz, g_0_y_xyzzzz_yyzz, g_0_y_xyzzzz_yzzz, g_0_y_xyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxyzzzz_xxxx[k] = -g_0_y_xyzzzz_xxxx[k] * ab_x + g_0_y_xyzzzz_xxxxx[k];

                g_0_y_xxyzzzz_xxxy[k] = -g_0_y_xyzzzz_xxxy[k] * ab_x + g_0_y_xyzzzz_xxxxy[k];

                g_0_y_xxyzzzz_xxxz[k] = -g_0_y_xyzzzz_xxxz[k] * ab_x + g_0_y_xyzzzz_xxxxz[k];

                g_0_y_xxyzzzz_xxyy[k] = -g_0_y_xyzzzz_xxyy[k] * ab_x + g_0_y_xyzzzz_xxxyy[k];

                g_0_y_xxyzzzz_xxyz[k] = -g_0_y_xyzzzz_xxyz[k] * ab_x + g_0_y_xyzzzz_xxxyz[k];

                g_0_y_xxyzzzz_xxzz[k] = -g_0_y_xyzzzz_xxzz[k] * ab_x + g_0_y_xyzzzz_xxxzz[k];

                g_0_y_xxyzzzz_xyyy[k] = -g_0_y_xyzzzz_xyyy[k] * ab_x + g_0_y_xyzzzz_xxyyy[k];

                g_0_y_xxyzzzz_xyyz[k] = -g_0_y_xyzzzz_xyyz[k] * ab_x + g_0_y_xyzzzz_xxyyz[k];

                g_0_y_xxyzzzz_xyzz[k] = -g_0_y_xyzzzz_xyzz[k] * ab_x + g_0_y_xyzzzz_xxyzz[k];

                g_0_y_xxyzzzz_xzzz[k] = -g_0_y_xyzzzz_xzzz[k] * ab_x + g_0_y_xyzzzz_xxzzz[k];

                g_0_y_xxyzzzz_yyyy[k] = -g_0_y_xyzzzz_yyyy[k] * ab_x + g_0_y_xyzzzz_xyyyy[k];

                g_0_y_xxyzzzz_yyyz[k] = -g_0_y_xyzzzz_yyyz[k] * ab_x + g_0_y_xyzzzz_xyyyz[k];

                g_0_y_xxyzzzz_yyzz[k] = -g_0_y_xyzzzz_yyzz[k] * ab_x + g_0_y_xyzzzz_xyyzz[k];

                g_0_y_xxyzzzz_yzzz[k] = -g_0_y_xyzzzz_yzzz[k] * ab_x + g_0_y_xyzzzz_xyzzz[k];

                g_0_y_xxyzzzz_zzzz[k] = -g_0_y_xyzzzz_zzzz[k] * ab_x + g_0_y_xyzzzz_xzzzz[k];
            }

            /// Set up 840-855 components of targeted buffer : cbuffer.data(

            auto g_0_y_xxzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 840 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 841 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 842 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 843 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 844 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 845 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 846 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 847 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 848 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 849 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 850 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 851 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 852 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 853 * ccomps * dcomps);

            auto g_0_y_xxzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 854 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xxzzzzz_xxxx, g_0_y_xxzzzzz_xxxy, g_0_y_xxzzzzz_xxxz, g_0_y_xxzzzzz_xxyy, g_0_y_xxzzzzz_xxyz, g_0_y_xxzzzzz_xxzz, g_0_y_xxzzzzz_xyyy, g_0_y_xxzzzzz_xyyz, g_0_y_xxzzzzz_xyzz, g_0_y_xxzzzzz_xzzz, g_0_y_xxzzzzz_yyyy, g_0_y_xxzzzzz_yyyz, g_0_y_xxzzzzz_yyzz, g_0_y_xxzzzzz_yzzz, g_0_y_xxzzzzz_zzzz, g_0_y_xzzzzz_xxxx, g_0_y_xzzzzz_xxxxx, g_0_y_xzzzzz_xxxxy, g_0_y_xzzzzz_xxxxz, g_0_y_xzzzzz_xxxy, g_0_y_xzzzzz_xxxyy, g_0_y_xzzzzz_xxxyz, g_0_y_xzzzzz_xxxz, g_0_y_xzzzzz_xxxzz, g_0_y_xzzzzz_xxyy, g_0_y_xzzzzz_xxyyy, g_0_y_xzzzzz_xxyyz, g_0_y_xzzzzz_xxyz, g_0_y_xzzzzz_xxyzz, g_0_y_xzzzzz_xxzz, g_0_y_xzzzzz_xxzzz, g_0_y_xzzzzz_xyyy, g_0_y_xzzzzz_xyyyy, g_0_y_xzzzzz_xyyyz, g_0_y_xzzzzz_xyyz, g_0_y_xzzzzz_xyyzz, g_0_y_xzzzzz_xyzz, g_0_y_xzzzzz_xyzzz, g_0_y_xzzzzz_xzzz, g_0_y_xzzzzz_xzzzz, g_0_y_xzzzzz_yyyy, g_0_y_xzzzzz_yyyz, g_0_y_xzzzzz_yyzz, g_0_y_xzzzzz_yzzz, g_0_y_xzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xxzzzzz_xxxx[k] = -g_0_y_xzzzzz_xxxx[k] * ab_x + g_0_y_xzzzzz_xxxxx[k];

                g_0_y_xxzzzzz_xxxy[k] = -g_0_y_xzzzzz_xxxy[k] * ab_x + g_0_y_xzzzzz_xxxxy[k];

                g_0_y_xxzzzzz_xxxz[k] = -g_0_y_xzzzzz_xxxz[k] * ab_x + g_0_y_xzzzzz_xxxxz[k];

                g_0_y_xxzzzzz_xxyy[k] = -g_0_y_xzzzzz_xxyy[k] * ab_x + g_0_y_xzzzzz_xxxyy[k];

                g_0_y_xxzzzzz_xxyz[k] = -g_0_y_xzzzzz_xxyz[k] * ab_x + g_0_y_xzzzzz_xxxyz[k];

                g_0_y_xxzzzzz_xxzz[k] = -g_0_y_xzzzzz_xxzz[k] * ab_x + g_0_y_xzzzzz_xxxzz[k];

                g_0_y_xxzzzzz_xyyy[k] = -g_0_y_xzzzzz_xyyy[k] * ab_x + g_0_y_xzzzzz_xxyyy[k];

                g_0_y_xxzzzzz_xyyz[k] = -g_0_y_xzzzzz_xyyz[k] * ab_x + g_0_y_xzzzzz_xxyyz[k];

                g_0_y_xxzzzzz_xyzz[k] = -g_0_y_xzzzzz_xyzz[k] * ab_x + g_0_y_xzzzzz_xxyzz[k];

                g_0_y_xxzzzzz_xzzz[k] = -g_0_y_xzzzzz_xzzz[k] * ab_x + g_0_y_xzzzzz_xxzzz[k];

                g_0_y_xxzzzzz_yyyy[k] = -g_0_y_xzzzzz_yyyy[k] * ab_x + g_0_y_xzzzzz_xyyyy[k];

                g_0_y_xxzzzzz_yyyz[k] = -g_0_y_xzzzzz_yyyz[k] * ab_x + g_0_y_xzzzzz_xyyyz[k];

                g_0_y_xxzzzzz_yyzz[k] = -g_0_y_xzzzzz_yyzz[k] * ab_x + g_0_y_xzzzzz_xyyzz[k];

                g_0_y_xxzzzzz_yzzz[k] = -g_0_y_xzzzzz_yzzz[k] * ab_x + g_0_y_xzzzzz_xyzzz[k];

                g_0_y_xxzzzzz_zzzz[k] = -g_0_y_xzzzzz_zzzz[k] * ab_x + g_0_y_xzzzzz_xzzzz[k];
            }

            /// Set up 855-870 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyyy_xxxx = cbuffer.data(kg_geom_01_off + 855 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xxxy = cbuffer.data(kg_geom_01_off + 856 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xxxz = cbuffer.data(kg_geom_01_off + 857 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xxyy = cbuffer.data(kg_geom_01_off + 858 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xxyz = cbuffer.data(kg_geom_01_off + 859 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xxzz = cbuffer.data(kg_geom_01_off + 860 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xyyy = cbuffer.data(kg_geom_01_off + 861 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xyyz = cbuffer.data(kg_geom_01_off + 862 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xyzz = cbuffer.data(kg_geom_01_off + 863 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_xzzz = cbuffer.data(kg_geom_01_off + 864 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_yyyy = cbuffer.data(kg_geom_01_off + 865 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_yyyz = cbuffer.data(kg_geom_01_off + 866 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_yyzz = cbuffer.data(kg_geom_01_off + 867 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_yzzz = cbuffer.data(kg_geom_01_off + 868 * ccomps * dcomps);

            auto g_0_y_xyyyyyy_zzzz = cbuffer.data(kg_geom_01_off + 869 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyyy_xxxx, g_0_y_xyyyyyy_xxxy, g_0_y_xyyyyyy_xxxz, g_0_y_xyyyyyy_xxyy, g_0_y_xyyyyyy_xxyz, g_0_y_xyyyyyy_xxzz, g_0_y_xyyyyyy_xyyy, g_0_y_xyyyyyy_xyyz, g_0_y_xyyyyyy_xyzz, g_0_y_xyyyyyy_xzzz, g_0_y_xyyyyyy_yyyy, g_0_y_xyyyyyy_yyyz, g_0_y_xyyyyyy_yyzz, g_0_y_xyyyyyy_yzzz, g_0_y_xyyyyyy_zzzz, g_0_y_yyyyyy_xxxx, g_0_y_yyyyyy_xxxxx, g_0_y_yyyyyy_xxxxy, g_0_y_yyyyyy_xxxxz, g_0_y_yyyyyy_xxxy, g_0_y_yyyyyy_xxxyy, g_0_y_yyyyyy_xxxyz, g_0_y_yyyyyy_xxxz, g_0_y_yyyyyy_xxxzz, g_0_y_yyyyyy_xxyy, g_0_y_yyyyyy_xxyyy, g_0_y_yyyyyy_xxyyz, g_0_y_yyyyyy_xxyz, g_0_y_yyyyyy_xxyzz, g_0_y_yyyyyy_xxzz, g_0_y_yyyyyy_xxzzz, g_0_y_yyyyyy_xyyy, g_0_y_yyyyyy_xyyyy, g_0_y_yyyyyy_xyyyz, g_0_y_yyyyyy_xyyz, g_0_y_yyyyyy_xyyzz, g_0_y_yyyyyy_xyzz, g_0_y_yyyyyy_xyzzz, g_0_y_yyyyyy_xzzz, g_0_y_yyyyyy_xzzzz, g_0_y_yyyyyy_yyyy, g_0_y_yyyyyy_yyyz, g_0_y_yyyyyy_yyzz, g_0_y_yyyyyy_yzzz, g_0_y_yyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyyy_xxxx[k] = -g_0_y_yyyyyy_xxxx[k] * ab_x + g_0_y_yyyyyy_xxxxx[k];

                g_0_y_xyyyyyy_xxxy[k] = -g_0_y_yyyyyy_xxxy[k] * ab_x + g_0_y_yyyyyy_xxxxy[k];

                g_0_y_xyyyyyy_xxxz[k] = -g_0_y_yyyyyy_xxxz[k] * ab_x + g_0_y_yyyyyy_xxxxz[k];

                g_0_y_xyyyyyy_xxyy[k] = -g_0_y_yyyyyy_xxyy[k] * ab_x + g_0_y_yyyyyy_xxxyy[k];

                g_0_y_xyyyyyy_xxyz[k] = -g_0_y_yyyyyy_xxyz[k] * ab_x + g_0_y_yyyyyy_xxxyz[k];

                g_0_y_xyyyyyy_xxzz[k] = -g_0_y_yyyyyy_xxzz[k] * ab_x + g_0_y_yyyyyy_xxxzz[k];

                g_0_y_xyyyyyy_xyyy[k] = -g_0_y_yyyyyy_xyyy[k] * ab_x + g_0_y_yyyyyy_xxyyy[k];

                g_0_y_xyyyyyy_xyyz[k] = -g_0_y_yyyyyy_xyyz[k] * ab_x + g_0_y_yyyyyy_xxyyz[k];

                g_0_y_xyyyyyy_xyzz[k] = -g_0_y_yyyyyy_xyzz[k] * ab_x + g_0_y_yyyyyy_xxyzz[k];

                g_0_y_xyyyyyy_xzzz[k] = -g_0_y_yyyyyy_xzzz[k] * ab_x + g_0_y_yyyyyy_xxzzz[k];

                g_0_y_xyyyyyy_yyyy[k] = -g_0_y_yyyyyy_yyyy[k] * ab_x + g_0_y_yyyyyy_xyyyy[k];

                g_0_y_xyyyyyy_yyyz[k] = -g_0_y_yyyyyy_yyyz[k] * ab_x + g_0_y_yyyyyy_xyyyz[k];

                g_0_y_xyyyyyy_yyzz[k] = -g_0_y_yyyyyy_yyzz[k] * ab_x + g_0_y_yyyyyy_xyyzz[k];

                g_0_y_xyyyyyy_yzzz[k] = -g_0_y_yyyyyy_yzzz[k] * ab_x + g_0_y_yyyyyy_xyzzz[k];

                g_0_y_xyyyyyy_zzzz[k] = -g_0_y_yyyyyy_zzzz[k] * ab_x + g_0_y_yyyyyy_xzzzz[k];
            }

            /// Set up 870-885 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyyz_xxxx = cbuffer.data(kg_geom_01_off + 870 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xxxy = cbuffer.data(kg_geom_01_off + 871 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xxxz = cbuffer.data(kg_geom_01_off + 872 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xxyy = cbuffer.data(kg_geom_01_off + 873 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xxyz = cbuffer.data(kg_geom_01_off + 874 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xxzz = cbuffer.data(kg_geom_01_off + 875 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xyyy = cbuffer.data(kg_geom_01_off + 876 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xyyz = cbuffer.data(kg_geom_01_off + 877 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xyzz = cbuffer.data(kg_geom_01_off + 878 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_xzzz = cbuffer.data(kg_geom_01_off + 879 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_yyyy = cbuffer.data(kg_geom_01_off + 880 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_yyyz = cbuffer.data(kg_geom_01_off + 881 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_yyzz = cbuffer.data(kg_geom_01_off + 882 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_yzzz = cbuffer.data(kg_geom_01_off + 883 * ccomps * dcomps);

            auto g_0_y_xyyyyyz_zzzz = cbuffer.data(kg_geom_01_off + 884 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyyz_xxxx, g_0_y_xyyyyyz_xxxy, g_0_y_xyyyyyz_xxxz, g_0_y_xyyyyyz_xxyy, g_0_y_xyyyyyz_xxyz, g_0_y_xyyyyyz_xxzz, g_0_y_xyyyyyz_xyyy, g_0_y_xyyyyyz_xyyz, g_0_y_xyyyyyz_xyzz, g_0_y_xyyyyyz_xzzz, g_0_y_xyyyyyz_yyyy, g_0_y_xyyyyyz_yyyz, g_0_y_xyyyyyz_yyzz, g_0_y_xyyyyyz_yzzz, g_0_y_xyyyyyz_zzzz, g_0_y_yyyyyz_xxxx, g_0_y_yyyyyz_xxxxx, g_0_y_yyyyyz_xxxxy, g_0_y_yyyyyz_xxxxz, g_0_y_yyyyyz_xxxy, g_0_y_yyyyyz_xxxyy, g_0_y_yyyyyz_xxxyz, g_0_y_yyyyyz_xxxz, g_0_y_yyyyyz_xxxzz, g_0_y_yyyyyz_xxyy, g_0_y_yyyyyz_xxyyy, g_0_y_yyyyyz_xxyyz, g_0_y_yyyyyz_xxyz, g_0_y_yyyyyz_xxyzz, g_0_y_yyyyyz_xxzz, g_0_y_yyyyyz_xxzzz, g_0_y_yyyyyz_xyyy, g_0_y_yyyyyz_xyyyy, g_0_y_yyyyyz_xyyyz, g_0_y_yyyyyz_xyyz, g_0_y_yyyyyz_xyyzz, g_0_y_yyyyyz_xyzz, g_0_y_yyyyyz_xyzzz, g_0_y_yyyyyz_xzzz, g_0_y_yyyyyz_xzzzz, g_0_y_yyyyyz_yyyy, g_0_y_yyyyyz_yyyz, g_0_y_yyyyyz_yyzz, g_0_y_yyyyyz_yzzz, g_0_y_yyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyyz_xxxx[k] = -g_0_y_yyyyyz_xxxx[k] * ab_x + g_0_y_yyyyyz_xxxxx[k];

                g_0_y_xyyyyyz_xxxy[k] = -g_0_y_yyyyyz_xxxy[k] * ab_x + g_0_y_yyyyyz_xxxxy[k];

                g_0_y_xyyyyyz_xxxz[k] = -g_0_y_yyyyyz_xxxz[k] * ab_x + g_0_y_yyyyyz_xxxxz[k];

                g_0_y_xyyyyyz_xxyy[k] = -g_0_y_yyyyyz_xxyy[k] * ab_x + g_0_y_yyyyyz_xxxyy[k];

                g_0_y_xyyyyyz_xxyz[k] = -g_0_y_yyyyyz_xxyz[k] * ab_x + g_0_y_yyyyyz_xxxyz[k];

                g_0_y_xyyyyyz_xxzz[k] = -g_0_y_yyyyyz_xxzz[k] * ab_x + g_0_y_yyyyyz_xxxzz[k];

                g_0_y_xyyyyyz_xyyy[k] = -g_0_y_yyyyyz_xyyy[k] * ab_x + g_0_y_yyyyyz_xxyyy[k];

                g_0_y_xyyyyyz_xyyz[k] = -g_0_y_yyyyyz_xyyz[k] * ab_x + g_0_y_yyyyyz_xxyyz[k];

                g_0_y_xyyyyyz_xyzz[k] = -g_0_y_yyyyyz_xyzz[k] * ab_x + g_0_y_yyyyyz_xxyzz[k];

                g_0_y_xyyyyyz_xzzz[k] = -g_0_y_yyyyyz_xzzz[k] * ab_x + g_0_y_yyyyyz_xxzzz[k];

                g_0_y_xyyyyyz_yyyy[k] = -g_0_y_yyyyyz_yyyy[k] * ab_x + g_0_y_yyyyyz_xyyyy[k];

                g_0_y_xyyyyyz_yyyz[k] = -g_0_y_yyyyyz_yyyz[k] * ab_x + g_0_y_yyyyyz_xyyyz[k];

                g_0_y_xyyyyyz_yyzz[k] = -g_0_y_yyyyyz_yyzz[k] * ab_x + g_0_y_yyyyyz_xyyzz[k];

                g_0_y_xyyyyyz_yzzz[k] = -g_0_y_yyyyyz_yzzz[k] * ab_x + g_0_y_yyyyyz_xyzzz[k];

                g_0_y_xyyyyyz_zzzz[k] = -g_0_y_yyyyyz_zzzz[k] * ab_x + g_0_y_yyyyyz_xzzzz[k];
            }

            /// Set up 885-900 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyyzz_xxxx = cbuffer.data(kg_geom_01_off + 885 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xxxy = cbuffer.data(kg_geom_01_off + 886 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xxxz = cbuffer.data(kg_geom_01_off + 887 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xxyy = cbuffer.data(kg_geom_01_off + 888 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xxyz = cbuffer.data(kg_geom_01_off + 889 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xxzz = cbuffer.data(kg_geom_01_off + 890 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xyyy = cbuffer.data(kg_geom_01_off + 891 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xyyz = cbuffer.data(kg_geom_01_off + 892 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xyzz = cbuffer.data(kg_geom_01_off + 893 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_xzzz = cbuffer.data(kg_geom_01_off + 894 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_yyyy = cbuffer.data(kg_geom_01_off + 895 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_yyyz = cbuffer.data(kg_geom_01_off + 896 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_yyzz = cbuffer.data(kg_geom_01_off + 897 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_yzzz = cbuffer.data(kg_geom_01_off + 898 * ccomps * dcomps);

            auto g_0_y_xyyyyzz_zzzz = cbuffer.data(kg_geom_01_off + 899 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyyzz_xxxx, g_0_y_xyyyyzz_xxxy, g_0_y_xyyyyzz_xxxz, g_0_y_xyyyyzz_xxyy, g_0_y_xyyyyzz_xxyz, g_0_y_xyyyyzz_xxzz, g_0_y_xyyyyzz_xyyy, g_0_y_xyyyyzz_xyyz, g_0_y_xyyyyzz_xyzz, g_0_y_xyyyyzz_xzzz, g_0_y_xyyyyzz_yyyy, g_0_y_xyyyyzz_yyyz, g_0_y_xyyyyzz_yyzz, g_0_y_xyyyyzz_yzzz, g_0_y_xyyyyzz_zzzz, g_0_y_yyyyzz_xxxx, g_0_y_yyyyzz_xxxxx, g_0_y_yyyyzz_xxxxy, g_0_y_yyyyzz_xxxxz, g_0_y_yyyyzz_xxxy, g_0_y_yyyyzz_xxxyy, g_0_y_yyyyzz_xxxyz, g_0_y_yyyyzz_xxxz, g_0_y_yyyyzz_xxxzz, g_0_y_yyyyzz_xxyy, g_0_y_yyyyzz_xxyyy, g_0_y_yyyyzz_xxyyz, g_0_y_yyyyzz_xxyz, g_0_y_yyyyzz_xxyzz, g_0_y_yyyyzz_xxzz, g_0_y_yyyyzz_xxzzz, g_0_y_yyyyzz_xyyy, g_0_y_yyyyzz_xyyyy, g_0_y_yyyyzz_xyyyz, g_0_y_yyyyzz_xyyz, g_0_y_yyyyzz_xyyzz, g_0_y_yyyyzz_xyzz, g_0_y_yyyyzz_xyzzz, g_0_y_yyyyzz_xzzz, g_0_y_yyyyzz_xzzzz, g_0_y_yyyyzz_yyyy, g_0_y_yyyyzz_yyyz, g_0_y_yyyyzz_yyzz, g_0_y_yyyyzz_yzzz, g_0_y_yyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyyzz_xxxx[k] = -g_0_y_yyyyzz_xxxx[k] * ab_x + g_0_y_yyyyzz_xxxxx[k];

                g_0_y_xyyyyzz_xxxy[k] = -g_0_y_yyyyzz_xxxy[k] * ab_x + g_0_y_yyyyzz_xxxxy[k];

                g_0_y_xyyyyzz_xxxz[k] = -g_0_y_yyyyzz_xxxz[k] * ab_x + g_0_y_yyyyzz_xxxxz[k];

                g_0_y_xyyyyzz_xxyy[k] = -g_0_y_yyyyzz_xxyy[k] * ab_x + g_0_y_yyyyzz_xxxyy[k];

                g_0_y_xyyyyzz_xxyz[k] = -g_0_y_yyyyzz_xxyz[k] * ab_x + g_0_y_yyyyzz_xxxyz[k];

                g_0_y_xyyyyzz_xxzz[k] = -g_0_y_yyyyzz_xxzz[k] * ab_x + g_0_y_yyyyzz_xxxzz[k];

                g_0_y_xyyyyzz_xyyy[k] = -g_0_y_yyyyzz_xyyy[k] * ab_x + g_0_y_yyyyzz_xxyyy[k];

                g_0_y_xyyyyzz_xyyz[k] = -g_0_y_yyyyzz_xyyz[k] * ab_x + g_0_y_yyyyzz_xxyyz[k];

                g_0_y_xyyyyzz_xyzz[k] = -g_0_y_yyyyzz_xyzz[k] * ab_x + g_0_y_yyyyzz_xxyzz[k];

                g_0_y_xyyyyzz_xzzz[k] = -g_0_y_yyyyzz_xzzz[k] * ab_x + g_0_y_yyyyzz_xxzzz[k];

                g_0_y_xyyyyzz_yyyy[k] = -g_0_y_yyyyzz_yyyy[k] * ab_x + g_0_y_yyyyzz_xyyyy[k];

                g_0_y_xyyyyzz_yyyz[k] = -g_0_y_yyyyzz_yyyz[k] * ab_x + g_0_y_yyyyzz_xyyyz[k];

                g_0_y_xyyyyzz_yyzz[k] = -g_0_y_yyyyzz_yyzz[k] * ab_x + g_0_y_yyyyzz_xyyzz[k];

                g_0_y_xyyyyzz_yzzz[k] = -g_0_y_yyyyzz_yzzz[k] * ab_x + g_0_y_yyyyzz_xyzzz[k];

                g_0_y_xyyyyzz_zzzz[k] = -g_0_y_yyyyzz_zzzz[k] * ab_x + g_0_y_yyyyzz_xzzzz[k];
            }

            /// Set up 900-915 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyyzzz_xxxx = cbuffer.data(kg_geom_01_off + 900 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xxxy = cbuffer.data(kg_geom_01_off + 901 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xxxz = cbuffer.data(kg_geom_01_off + 902 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xxyy = cbuffer.data(kg_geom_01_off + 903 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xxyz = cbuffer.data(kg_geom_01_off + 904 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xxzz = cbuffer.data(kg_geom_01_off + 905 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xyyy = cbuffer.data(kg_geom_01_off + 906 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xyyz = cbuffer.data(kg_geom_01_off + 907 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xyzz = cbuffer.data(kg_geom_01_off + 908 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_xzzz = cbuffer.data(kg_geom_01_off + 909 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_yyyy = cbuffer.data(kg_geom_01_off + 910 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_yyyz = cbuffer.data(kg_geom_01_off + 911 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_yyzz = cbuffer.data(kg_geom_01_off + 912 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_yzzz = cbuffer.data(kg_geom_01_off + 913 * ccomps * dcomps);

            auto g_0_y_xyyyzzz_zzzz = cbuffer.data(kg_geom_01_off + 914 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyyzzz_xxxx, g_0_y_xyyyzzz_xxxy, g_0_y_xyyyzzz_xxxz, g_0_y_xyyyzzz_xxyy, g_0_y_xyyyzzz_xxyz, g_0_y_xyyyzzz_xxzz, g_0_y_xyyyzzz_xyyy, g_0_y_xyyyzzz_xyyz, g_0_y_xyyyzzz_xyzz, g_0_y_xyyyzzz_xzzz, g_0_y_xyyyzzz_yyyy, g_0_y_xyyyzzz_yyyz, g_0_y_xyyyzzz_yyzz, g_0_y_xyyyzzz_yzzz, g_0_y_xyyyzzz_zzzz, g_0_y_yyyzzz_xxxx, g_0_y_yyyzzz_xxxxx, g_0_y_yyyzzz_xxxxy, g_0_y_yyyzzz_xxxxz, g_0_y_yyyzzz_xxxy, g_0_y_yyyzzz_xxxyy, g_0_y_yyyzzz_xxxyz, g_0_y_yyyzzz_xxxz, g_0_y_yyyzzz_xxxzz, g_0_y_yyyzzz_xxyy, g_0_y_yyyzzz_xxyyy, g_0_y_yyyzzz_xxyyz, g_0_y_yyyzzz_xxyz, g_0_y_yyyzzz_xxyzz, g_0_y_yyyzzz_xxzz, g_0_y_yyyzzz_xxzzz, g_0_y_yyyzzz_xyyy, g_0_y_yyyzzz_xyyyy, g_0_y_yyyzzz_xyyyz, g_0_y_yyyzzz_xyyz, g_0_y_yyyzzz_xyyzz, g_0_y_yyyzzz_xyzz, g_0_y_yyyzzz_xyzzz, g_0_y_yyyzzz_xzzz, g_0_y_yyyzzz_xzzzz, g_0_y_yyyzzz_yyyy, g_0_y_yyyzzz_yyyz, g_0_y_yyyzzz_yyzz, g_0_y_yyyzzz_yzzz, g_0_y_yyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyyzzz_xxxx[k] = -g_0_y_yyyzzz_xxxx[k] * ab_x + g_0_y_yyyzzz_xxxxx[k];

                g_0_y_xyyyzzz_xxxy[k] = -g_0_y_yyyzzz_xxxy[k] * ab_x + g_0_y_yyyzzz_xxxxy[k];

                g_0_y_xyyyzzz_xxxz[k] = -g_0_y_yyyzzz_xxxz[k] * ab_x + g_0_y_yyyzzz_xxxxz[k];

                g_0_y_xyyyzzz_xxyy[k] = -g_0_y_yyyzzz_xxyy[k] * ab_x + g_0_y_yyyzzz_xxxyy[k];

                g_0_y_xyyyzzz_xxyz[k] = -g_0_y_yyyzzz_xxyz[k] * ab_x + g_0_y_yyyzzz_xxxyz[k];

                g_0_y_xyyyzzz_xxzz[k] = -g_0_y_yyyzzz_xxzz[k] * ab_x + g_0_y_yyyzzz_xxxzz[k];

                g_0_y_xyyyzzz_xyyy[k] = -g_0_y_yyyzzz_xyyy[k] * ab_x + g_0_y_yyyzzz_xxyyy[k];

                g_0_y_xyyyzzz_xyyz[k] = -g_0_y_yyyzzz_xyyz[k] * ab_x + g_0_y_yyyzzz_xxyyz[k];

                g_0_y_xyyyzzz_xyzz[k] = -g_0_y_yyyzzz_xyzz[k] * ab_x + g_0_y_yyyzzz_xxyzz[k];

                g_0_y_xyyyzzz_xzzz[k] = -g_0_y_yyyzzz_xzzz[k] * ab_x + g_0_y_yyyzzz_xxzzz[k];

                g_0_y_xyyyzzz_yyyy[k] = -g_0_y_yyyzzz_yyyy[k] * ab_x + g_0_y_yyyzzz_xyyyy[k];

                g_0_y_xyyyzzz_yyyz[k] = -g_0_y_yyyzzz_yyyz[k] * ab_x + g_0_y_yyyzzz_xyyyz[k];

                g_0_y_xyyyzzz_yyzz[k] = -g_0_y_yyyzzz_yyzz[k] * ab_x + g_0_y_yyyzzz_xyyzz[k];

                g_0_y_xyyyzzz_yzzz[k] = -g_0_y_yyyzzz_yzzz[k] * ab_x + g_0_y_yyyzzz_xyzzz[k];

                g_0_y_xyyyzzz_zzzz[k] = -g_0_y_yyyzzz_zzzz[k] * ab_x + g_0_y_yyyzzz_xzzzz[k];
            }

            /// Set up 915-930 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyyzzzz_xxxx = cbuffer.data(kg_geom_01_off + 915 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xxxy = cbuffer.data(kg_geom_01_off + 916 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xxxz = cbuffer.data(kg_geom_01_off + 917 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xxyy = cbuffer.data(kg_geom_01_off + 918 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xxyz = cbuffer.data(kg_geom_01_off + 919 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xxzz = cbuffer.data(kg_geom_01_off + 920 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xyyy = cbuffer.data(kg_geom_01_off + 921 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xyyz = cbuffer.data(kg_geom_01_off + 922 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xyzz = cbuffer.data(kg_geom_01_off + 923 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_xzzz = cbuffer.data(kg_geom_01_off + 924 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_yyyy = cbuffer.data(kg_geom_01_off + 925 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_yyyz = cbuffer.data(kg_geom_01_off + 926 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_yyzz = cbuffer.data(kg_geom_01_off + 927 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_yzzz = cbuffer.data(kg_geom_01_off + 928 * ccomps * dcomps);

            auto g_0_y_xyyzzzz_zzzz = cbuffer.data(kg_geom_01_off + 929 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyyzzzz_xxxx, g_0_y_xyyzzzz_xxxy, g_0_y_xyyzzzz_xxxz, g_0_y_xyyzzzz_xxyy, g_0_y_xyyzzzz_xxyz, g_0_y_xyyzzzz_xxzz, g_0_y_xyyzzzz_xyyy, g_0_y_xyyzzzz_xyyz, g_0_y_xyyzzzz_xyzz, g_0_y_xyyzzzz_xzzz, g_0_y_xyyzzzz_yyyy, g_0_y_xyyzzzz_yyyz, g_0_y_xyyzzzz_yyzz, g_0_y_xyyzzzz_yzzz, g_0_y_xyyzzzz_zzzz, g_0_y_yyzzzz_xxxx, g_0_y_yyzzzz_xxxxx, g_0_y_yyzzzz_xxxxy, g_0_y_yyzzzz_xxxxz, g_0_y_yyzzzz_xxxy, g_0_y_yyzzzz_xxxyy, g_0_y_yyzzzz_xxxyz, g_0_y_yyzzzz_xxxz, g_0_y_yyzzzz_xxxzz, g_0_y_yyzzzz_xxyy, g_0_y_yyzzzz_xxyyy, g_0_y_yyzzzz_xxyyz, g_0_y_yyzzzz_xxyz, g_0_y_yyzzzz_xxyzz, g_0_y_yyzzzz_xxzz, g_0_y_yyzzzz_xxzzz, g_0_y_yyzzzz_xyyy, g_0_y_yyzzzz_xyyyy, g_0_y_yyzzzz_xyyyz, g_0_y_yyzzzz_xyyz, g_0_y_yyzzzz_xyyzz, g_0_y_yyzzzz_xyzz, g_0_y_yyzzzz_xyzzz, g_0_y_yyzzzz_xzzz, g_0_y_yyzzzz_xzzzz, g_0_y_yyzzzz_yyyy, g_0_y_yyzzzz_yyyz, g_0_y_yyzzzz_yyzz, g_0_y_yyzzzz_yzzz, g_0_y_yyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyyzzzz_xxxx[k] = -g_0_y_yyzzzz_xxxx[k] * ab_x + g_0_y_yyzzzz_xxxxx[k];

                g_0_y_xyyzzzz_xxxy[k] = -g_0_y_yyzzzz_xxxy[k] * ab_x + g_0_y_yyzzzz_xxxxy[k];

                g_0_y_xyyzzzz_xxxz[k] = -g_0_y_yyzzzz_xxxz[k] * ab_x + g_0_y_yyzzzz_xxxxz[k];

                g_0_y_xyyzzzz_xxyy[k] = -g_0_y_yyzzzz_xxyy[k] * ab_x + g_0_y_yyzzzz_xxxyy[k];

                g_0_y_xyyzzzz_xxyz[k] = -g_0_y_yyzzzz_xxyz[k] * ab_x + g_0_y_yyzzzz_xxxyz[k];

                g_0_y_xyyzzzz_xxzz[k] = -g_0_y_yyzzzz_xxzz[k] * ab_x + g_0_y_yyzzzz_xxxzz[k];

                g_0_y_xyyzzzz_xyyy[k] = -g_0_y_yyzzzz_xyyy[k] * ab_x + g_0_y_yyzzzz_xxyyy[k];

                g_0_y_xyyzzzz_xyyz[k] = -g_0_y_yyzzzz_xyyz[k] * ab_x + g_0_y_yyzzzz_xxyyz[k];

                g_0_y_xyyzzzz_xyzz[k] = -g_0_y_yyzzzz_xyzz[k] * ab_x + g_0_y_yyzzzz_xxyzz[k];

                g_0_y_xyyzzzz_xzzz[k] = -g_0_y_yyzzzz_xzzz[k] * ab_x + g_0_y_yyzzzz_xxzzz[k];

                g_0_y_xyyzzzz_yyyy[k] = -g_0_y_yyzzzz_yyyy[k] * ab_x + g_0_y_yyzzzz_xyyyy[k];

                g_0_y_xyyzzzz_yyyz[k] = -g_0_y_yyzzzz_yyyz[k] * ab_x + g_0_y_yyzzzz_xyyyz[k];

                g_0_y_xyyzzzz_yyzz[k] = -g_0_y_yyzzzz_yyzz[k] * ab_x + g_0_y_yyzzzz_xyyzz[k];

                g_0_y_xyyzzzz_yzzz[k] = -g_0_y_yyzzzz_yzzz[k] * ab_x + g_0_y_yyzzzz_xyzzz[k];

                g_0_y_xyyzzzz_zzzz[k] = -g_0_y_yyzzzz_zzzz[k] * ab_x + g_0_y_yyzzzz_xzzzz[k];
            }

            /// Set up 930-945 components of targeted buffer : cbuffer.data(

            auto g_0_y_xyzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 930 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 931 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 932 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 933 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 934 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 935 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 936 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 937 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 938 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 939 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 940 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 941 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 942 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 943 * ccomps * dcomps);

            auto g_0_y_xyzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 944 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xyzzzzz_xxxx, g_0_y_xyzzzzz_xxxy, g_0_y_xyzzzzz_xxxz, g_0_y_xyzzzzz_xxyy, g_0_y_xyzzzzz_xxyz, g_0_y_xyzzzzz_xxzz, g_0_y_xyzzzzz_xyyy, g_0_y_xyzzzzz_xyyz, g_0_y_xyzzzzz_xyzz, g_0_y_xyzzzzz_xzzz, g_0_y_xyzzzzz_yyyy, g_0_y_xyzzzzz_yyyz, g_0_y_xyzzzzz_yyzz, g_0_y_xyzzzzz_yzzz, g_0_y_xyzzzzz_zzzz, g_0_y_yzzzzz_xxxx, g_0_y_yzzzzz_xxxxx, g_0_y_yzzzzz_xxxxy, g_0_y_yzzzzz_xxxxz, g_0_y_yzzzzz_xxxy, g_0_y_yzzzzz_xxxyy, g_0_y_yzzzzz_xxxyz, g_0_y_yzzzzz_xxxz, g_0_y_yzzzzz_xxxzz, g_0_y_yzzzzz_xxyy, g_0_y_yzzzzz_xxyyy, g_0_y_yzzzzz_xxyyz, g_0_y_yzzzzz_xxyz, g_0_y_yzzzzz_xxyzz, g_0_y_yzzzzz_xxzz, g_0_y_yzzzzz_xxzzz, g_0_y_yzzzzz_xyyy, g_0_y_yzzzzz_xyyyy, g_0_y_yzzzzz_xyyyz, g_0_y_yzzzzz_xyyz, g_0_y_yzzzzz_xyyzz, g_0_y_yzzzzz_xyzz, g_0_y_yzzzzz_xyzzz, g_0_y_yzzzzz_xzzz, g_0_y_yzzzzz_xzzzz, g_0_y_yzzzzz_yyyy, g_0_y_yzzzzz_yyyz, g_0_y_yzzzzz_yyzz, g_0_y_yzzzzz_yzzz, g_0_y_yzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xyzzzzz_xxxx[k] = -g_0_y_yzzzzz_xxxx[k] * ab_x + g_0_y_yzzzzz_xxxxx[k];

                g_0_y_xyzzzzz_xxxy[k] = -g_0_y_yzzzzz_xxxy[k] * ab_x + g_0_y_yzzzzz_xxxxy[k];

                g_0_y_xyzzzzz_xxxz[k] = -g_0_y_yzzzzz_xxxz[k] * ab_x + g_0_y_yzzzzz_xxxxz[k];

                g_0_y_xyzzzzz_xxyy[k] = -g_0_y_yzzzzz_xxyy[k] * ab_x + g_0_y_yzzzzz_xxxyy[k];

                g_0_y_xyzzzzz_xxyz[k] = -g_0_y_yzzzzz_xxyz[k] * ab_x + g_0_y_yzzzzz_xxxyz[k];

                g_0_y_xyzzzzz_xxzz[k] = -g_0_y_yzzzzz_xxzz[k] * ab_x + g_0_y_yzzzzz_xxxzz[k];

                g_0_y_xyzzzzz_xyyy[k] = -g_0_y_yzzzzz_xyyy[k] * ab_x + g_0_y_yzzzzz_xxyyy[k];

                g_0_y_xyzzzzz_xyyz[k] = -g_0_y_yzzzzz_xyyz[k] * ab_x + g_0_y_yzzzzz_xxyyz[k];

                g_0_y_xyzzzzz_xyzz[k] = -g_0_y_yzzzzz_xyzz[k] * ab_x + g_0_y_yzzzzz_xxyzz[k];

                g_0_y_xyzzzzz_xzzz[k] = -g_0_y_yzzzzz_xzzz[k] * ab_x + g_0_y_yzzzzz_xxzzz[k];

                g_0_y_xyzzzzz_yyyy[k] = -g_0_y_yzzzzz_yyyy[k] * ab_x + g_0_y_yzzzzz_xyyyy[k];

                g_0_y_xyzzzzz_yyyz[k] = -g_0_y_yzzzzz_yyyz[k] * ab_x + g_0_y_yzzzzz_xyyyz[k];

                g_0_y_xyzzzzz_yyzz[k] = -g_0_y_yzzzzz_yyzz[k] * ab_x + g_0_y_yzzzzz_xyyzz[k];

                g_0_y_xyzzzzz_yzzz[k] = -g_0_y_yzzzzz_yzzz[k] * ab_x + g_0_y_yzzzzz_xyzzz[k];

                g_0_y_xyzzzzz_zzzz[k] = -g_0_y_yzzzzz_zzzz[k] * ab_x + g_0_y_yzzzzz_xzzzz[k];
            }

            /// Set up 945-960 components of targeted buffer : cbuffer.data(

            auto g_0_y_xzzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 945 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 946 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 947 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 948 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 949 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 950 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 951 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 952 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 953 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 954 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 955 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 956 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 957 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 958 * ccomps * dcomps);

            auto g_0_y_xzzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 959 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_xzzzzzz_xxxx, g_0_y_xzzzzzz_xxxy, g_0_y_xzzzzzz_xxxz, g_0_y_xzzzzzz_xxyy, g_0_y_xzzzzzz_xxyz, g_0_y_xzzzzzz_xxzz, g_0_y_xzzzzzz_xyyy, g_0_y_xzzzzzz_xyyz, g_0_y_xzzzzzz_xyzz, g_0_y_xzzzzzz_xzzz, g_0_y_xzzzzzz_yyyy, g_0_y_xzzzzzz_yyyz, g_0_y_xzzzzzz_yyzz, g_0_y_xzzzzzz_yzzz, g_0_y_xzzzzzz_zzzz, g_0_y_zzzzzz_xxxx, g_0_y_zzzzzz_xxxxx, g_0_y_zzzzzz_xxxxy, g_0_y_zzzzzz_xxxxz, g_0_y_zzzzzz_xxxy, g_0_y_zzzzzz_xxxyy, g_0_y_zzzzzz_xxxyz, g_0_y_zzzzzz_xxxz, g_0_y_zzzzzz_xxxzz, g_0_y_zzzzzz_xxyy, g_0_y_zzzzzz_xxyyy, g_0_y_zzzzzz_xxyyz, g_0_y_zzzzzz_xxyz, g_0_y_zzzzzz_xxyzz, g_0_y_zzzzzz_xxzz, g_0_y_zzzzzz_xxzzz, g_0_y_zzzzzz_xyyy, g_0_y_zzzzzz_xyyyy, g_0_y_zzzzzz_xyyyz, g_0_y_zzzzzz_xyyz, g_0_y_zzzzzz_xyyzz, g_0_y_zzzzzz_xyzz, g_0_y_zzzzzz_xyzzz, g_0_y_zzzzzz_xzzz, g_0_y_zzzzzz_xzzzz, g_0_y_zzzzzz_yyyy, g_0_y_zzzzzz_yyyz, g_0_y_zzzzzz_yyzz, g_0_y_zzzzzz_yzzz, g_0_y_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_xzzzzzz_xxxx[k] = -g_0_y_zzzzzz_xxxx[k] * ab_x + g_0_y_zzzzzz_xxxxx[k];

                g_0_y_xzzzzzz_xxxy[k] = -g_0_y_zzzzzz_xxxy[k] * ab_x + g_0_y_zzzzzz_xxxxy[k];

                g_0_y_xzzzzzz_xxxz[k] = -g_0_y_zzzzzz_xxxz[k] * ab_x + g_0_y_zzzzzz_xxxxz[k];

                g_0_y_xzzzzzz_xxyy[k] = -g_0_y_zzzzzz_xxyy[k] * ab_x + g_0_y_zzzzzz_xxxyy[k];

                g_0_y_xzzzzzz_xxyz[k] = -g_0_y_zzzzzz_xxyz[k] * ab_x + g_0_y_zzzzzz_xxxyz[k];

                g_0_y_xzzzzzz_xxzz[k] = -g_0_y_zzzzzz_xxzz[k] * ab_x + g_0_y_zzzzzz_xxxzz[k];

                g_0_y_xzzzzzz_xyyy[k] = -g_0_y_zzzzzz_xyyy[k] * ab_x + g_0_y_zzzzzz_xxyyy[k];

                g_0_y_xzzzzzz_xyyz[k] = -g_0_y_zzzzzz_xyyz[k] * ab_x + g_0_y_zzzzzz_xxyyz[k];

                g_0_y_xzzzzzz_xyzz[k] = -g_0_y_zzzzzz_xyzz[k] * ab_x + g_0_y_zzzzzz_xxyzz[k];

                g_0_y_xzzzzzz_xzzz[k] = -g_0_y_zzzzzz_xzzz[k] * ab_x + g_0_y_zzzzzz_xxzzz[k];

                g_0_y_xzzzzzz_yyyy[k] = -g_0_y_zzzzzz_yyyy[k] * ab_x + g_0_y_zzzzzz_xyyyy[k];

                g_0_y_xzzzzzz_yyyz[k] = -g_0_y_zzzzzz_yyyz[k] * ab_x + g_0_y_zzzzzz_xyyyz[k];

                g_0_y_xzzzzzz_yyzz[k] = -g_0_y_zzzzzz_yyzz[k] * ab_x + g_0_y_zzzzzz_xyyzz[k];

                g_0_y_xzzzzzz_yzzz[k] = -g_0_y_zzzzzz_yzzz[k] * ab_x + g_0_y_zzzzzz_xyzzz[k];

                g_0_y_xzzzzzz_zzzz[k] = -g_0_y_zzzzzz_zzzz[k] * ab_x + g_0_y_zzzzzz_xzzzz[k];
            }

            /// Set up 960-975 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyyy_xxxx = cbuffer.data(kg_geom_01_off + 960 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xxxy = cbuffer.data(kg_geom_01_off + 961 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xxxz = cbuffer.data(kg_geom_01_off + 962 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xxyy = cbuffer.data(kg_geom_01_off + 963 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xxyz = cbuffer.data(kg_geom_01_off + 964 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xxzz = cbuffer.data(kg_geom_01_off + 965 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xyyy = cbuffer.data(kg_geom_01_off + 966 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xyyz = cbuffer.data(kg_geom_01_off + 967 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xyzz = cbuffer.data(kg_geom_01_off + 968 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_xzzz = cbuffer.data(kg_geom_01_off + 969 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_yyyy = cbuffer.data(kg_geom_01_off + 970 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_yyyz = cbuffer.data(kg_geom_01_off + 971 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_yyzz = cbuffer.data(kg_geom_01_off + 972 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_yzzz = cbuffer.data(kg_geom_01_off + 973 * ccomps * dcomps);

            auto g_0_y_yyyyyyy_zzzz = cbuffer.data(kg_geom_01_off + 974 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyy_xxxx, g_0_y_yyyyyy_xxxxy, g_0_y_yyyyyy_xxxy, g_0_y_yyyyyy_xxxyy, g_0_y_yyyyyy_xxxyz, g_0_y_yyyyyy_xxxz, g_0_y_yyyyyy_xxyy, g_0_y_yyyyyy_xxyyy, g_0_y_yyyyyy_xxyyz, g_0_y_yyyyyy_xxyz, g_0_y_yyyyyy_xxyzz, g_0_y_yyyyyy_xxzz, g_0_y_yyyyyy_xyyy, g_0_y_yyyyyy_xyyyy, g_0_y_yyyyyy_xyyyz, g_0_y_yyyyyy_xyyz, g_0_y_yyyyyy_xyyzz, g_0_y_yyyyyy_xyzz, g_0_y_yyyyyy_xyzzz, g_0_y_yyyyyy_xzzz, g_0_y_yyyyyy_yyyy, g_0_y_yyyyyy_yyyyy, g_0_y_yyyyyy_yyyyz, g_0_y_yyyyyy_yyyz, g_0_y_yyyyyy_yyyzz, g_0_y_yyyyyy_yyzz, g_0_y_yyyyyy_yyzzz, g_0_y_yyyyyy_yzzz, g_0_y_yyyyyy_yzzzz, g_0_y_yyyyyy_zzzz, g_0_y_yyyyyyy_xxxx, g_0_y_yyyyyyy_xxxy, g_0_y_yyyyyyy_xxxz, g_0_y_yyyyyyy_xxyy, g_0_y_yyyyyyy_xxyz, g_0_y_yyyyyyy_xxzz, g_0_y_yyyyyyy_xyyy, g_0_y_yyyyyyy_xyyz, g_0_y_yyyyyyy_xyzz, g_0_y_yyyyyyy_xzzz, g_0_y_yyyyyyy_yyyy, g_0_y_yyyyyyy_yyyz, g_0_y_yyyyyyy_yyzz, g_0_y_yyyyyyy_yzzz, g_0_y_yyyyyyy_zzzz, g_yyyyyy_xxxx, g_yyyyyy_xxxy, g_yyyyyy_xxxz, g_yyyyyy_xxyy, g_yyyyyy_xxyz, g_yyyyyy_xxzz, g_yyyyyy_xyyy, g_yyyyyy_xyyz, g_yyyyyy_xyzz, g_yyyyyy_xzzz, g_yyyyyy_yyyy, g_yyyyyy_yyyz, g_yyyyyy_yyzz, g_yyyyyy_yzzz, g_yyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyyy_xxxx[k] = g_yyyyyy_xxxx[k] - g_0_y_yyyyyy_xxxx[k] * ab_y + g_0_y_yyyyyy_xxxxy[k];

                g_0_y_yyyyyyy_xxxy[k] = g_yyyyyy_xxxy[k] - g_0_y_yyyyyy_xxxy[k] * ab_y + g_0_y_yyyyyy_xxxyy[k];

                g_0_y_yyyyyyy_xxxz[k] = g_yyyyyy_xxxz[k] - g_0_y_yyyyyy_xxxz[k] * ab_y + g_0_y_yyyyyy_xxxyz[k];

                g_0_y_yyyyyyy_xxyy[k] = g_yyyyyy_xxyy[k] - g_0_y_yyyyyy_xxyy[k] * ab_y + g_0_y_yyyyyy_xxyyy[k];

                g_0_y_yyyyyyy_xxyz[k] = g_yyyyyy_xxyz[k] - g_0_y_yyyyyy_xxyz[k] * ab_y + g_0_y_yyyyyy_xxyyz[k];

                g_0_y_yyyyyyy_xxzz[k] = g_yyyyyy_xxzz[k] - g_0_y_yyyyyy_xxzz[k] * ab_y + g_0_y_yyyyyy_xxyzz[k];

                g_0_y_yyyyyyy_xyyy[k] = g_yyyyyy_xyyy[k] - g_0_y_yyyyyy_xyyy[k] * ab_y + g_0_y_yyyyyy_xyyyy[k];

                g_0_y_yyyyyyy_xyyz[k] = g_yyyyyy_xyyz[k] - g_0_y_yyyyyy_xyyz[k] * ab_y + g_0_y_yyyyyy_xyyyz[k];

                g_0_y_yyyyyyy_xyzz[k] = g_yyyyyy_xyzz[k] - g_0_y_yyyyyy_xyzz[k] * ab_y + g_0_y_yyyyyy_xyyzz[k];

                g_0_y_yyyyyyy_xzzz[k] = g_yyyyyy_xzzz[k] - g_0_y_yyyyyy_xzzz[k] * ab_y + g_0_y_yyyyyy_xyzzz[k];

                g_0_y_yyyyyyy_yyyy[k] = g_yyyyyy_yyyy[k] - g_0_y_yyyyyy_yyyy[k] * ab_y + g_0_y_yyyyyy_yyyyy[k];

                g_0_y_yyyyyyy_yyyz[k] = g_yyyyyy_yyyz[k] - g_0_y_yyyyyy_yyyz[k] * ab_y + g_0_y_yyyyyy_yyyyz[k];

                g_0_y_yyyyyyy_yyzz[k] = g_yyyyyy_yyzz[k] - g_0_y_yyyyyy_yyzz[k] * ab_y + g_0_y_yyyyyy_yyyzz[k];

                g_0_y_yyyyyyy_yzzz[k] = g_yyyyyy_yzzz[k] - g_0_y_yyyyyy_yzzz[k] * ab_y + g_0_y_yyyyyy_yyzzz[k];

                g_0_y_yyyyyyy_zzzz[k] = g_yyyyyy_zzzz[k] - g_0_y_yyyyyy_zzzz[k] * ab_y + g_0_y_yyyyyy_yzzzz[k];
            }

            /// Set up 975-990 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyyz_xxxx = cbuffer.data(kg_geom_01_off + 975 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xxxy = cbuffer.data(kg_geom_01_off + 976 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xxxz = cbuffer.data(kg_geom_01_off + 977 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xxyy = cbuffer.data(kg_geom_01_off + 978 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xxyz = cbuffer.data(kg_geom_01_off + 979 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xxzz = cbuffer.data(kg_geom_01_off + 980 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xyyy = cbuffer.data(kg_geom_01_off + 981 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xyyz = cbuffer.data(kg_geom_01_off + 982 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xyzz = cbuffer.data(kg_geom_01_off + 983 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_xzzz = cbuffer.data(kg_geom_01_off + 984 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_yyyy = cbuffer.data(kg_geom_01_off + 985 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_yyyz = cbuffer.data(kg_geom_01_off + 986 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_yyzz = cbuffer.data(kg_geom_01_off + 987 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_yzzz = cbuffer.data(kg_geom_01_off + 988 * ccomps * dcomps);

            auto g_0_y_yyyyyyz_zzzz = cbuffer.data(kg_geom_01_off + 989 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyy_xxxx, g_0_y_yyyyyy_xxxxz, g_0_y_yyyyyy_xxxy, g_0_y_yyyyyy_xxxyz, g_0_y_yyyyyy_xxxz, g_0_y_yyyyyy_xxxzz, g_0_y_yyyyyy_xxyy, g_0_y_yyyyyy_xxyyz, g_0_y_yyyyyy_xxyz, g_0_y_yyyyyy_xxyzz, g_0_y_yyyyyy_xxzz, g_0_y_yyyyyy_xxzzz, g_0_y_yyyyyy_xyyy, g_0_y_yyyyyy_xyyyz, g_0_y_yyyyyy_xyyz, g_0_y_yyyyyy_xyyzz, g_0_y_yyyyyy_xyzz, g_0_y_yyyyyy_xyzzz, g_0_y_yyyyyy_xzzz, g_0_y_yyyyyy_xzzzz, g_0_y_yyyyyy_yyyy, g_0_y_yyyyyy_yyyyz, g_0_y_yyyyyy_yyyz, g_0_y_yyyyyy_yyyzz, g_0_y_yyyyyy_yyzz, g_0_y_yyyyyy_yyzzz, g_0_y_yyyyyy_yzzz, g_0_y_yyyyyy_yzzzz, g_0_y_yyyyyy_zzzz, g_0_y_yyyyyy_zzzzz, g_0_y_yyyyyyz_xxxx, g_0_y_yyyyyyz_xxxy, g_0_y_yyyyyyz_xxxz, g_0_y_yyyyyyz_xxyy, g_0_y_yyyyyyz_xxyz, g_0_y_yyyyyyz_xxzz, g_0_y_yyyyyyz_xyyy, g_0_y_yyyyyyz_xyyz, g_0_y_yyyyyyz_xyzz, g_0_y_yyyyyyz_xzzz, g_0_y_yyyyyyz_yyyy, g_0_y_yyyyyyz_yyyz, g_0_y_yyyyyyz_yyzz, g_0_y_yyyyyyz_yzzz, g_0_y_yyyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyyz_xxxx[k] = -g_0_y_yyyyyy_xxxx[k] * ab_z + g_0_y_yyyyyy_xxxxz[k];

                g_0_y_yyyyyyz_xxxy[k] = -g_0_y_yyyyyy_xxxy[k] * ab_z + g_0_y_yyyyyy_xxxyz[k];

                g_0_y_yyyyyyz_xxxz[k] = -g_0_y_yyyyyy_xxxz[k] * ab_z + g_0_y_yyyyyy_xxxzz[k];

                g_0_y_yyyyyyz_xxyy[k] = -g_0_y_yyyyyy_xxyy[k] * ab_z + g_0_y_yyyyyy_xxyyz[k];

                g_0_y_yyyyyyz_xxyz[k] = -g_0_y_yyyyyy_xxyz[k] * ab_z + g_0_y_yyyyyy_xxyzz[k];

                g_0_y_yyyyyyz_xxzz[k] = -g_0_y_yyyyyy_xxzz[k] * ab_z + g_0_y_yyyyyy_xxzzz[k];

                g_0_y_yyyyyyz_xyyy[k] = -g_0_y_yyyyyy_xyyy[k] * ab_z + g_0_y_yyyyyy_xyyyz[k];

                g_0_y_yyyyyyz_xyyz[k] = -g_0_y_yyyyyy_xyyz[k] * ab_z + g_0_y_yyyyyy_xyyzz[k];

                g_0_y_yyyyyyz_xyzz[k] = -g_0_y_yyyyyy_xyzz[k] * ab_z + g_0_y_yyyyyy_xyzzz[k];

                g_0_y_yyyyyyz_xzzz[k] = -g_0_y_yyyyyy_xzzz[k] * ab_z + g_0_y_yyyyyy_xzzzz[k];

                g_0_y_yyyyyyz_yyyy[k] = -g_0_y_yyyyyy_yyyy[k] * ab_z + g_0_y_yyyyyy_yyyyz[k];

                g_0_y_yyyyyyz_yyyz[k] = -g_0_y_yyyyyy_yyyz[k] * ab_z + g_0_y_yyyyyy_yyyzz[k];

                g_0_y_yyyyyyz_yyzz[k] = -g_0_y_yyyyyy_yyzz[k] * ab_z + g_0_y_yyyyyy_yyzzz[k];

                g_0_y_yyyyyyz_yzzz[k] = -g_0_y_yyyyyy_yzzz[k] * ab_z + g_0_y_yyyyyy_yzzzz[k];

                g_0_y_yyyyyyz_zzzz[k] = -g_0_y_yyyyyy_zzzz[k] * ab_z + g_0_y_yyyyyy_zzzzz[k];
            }

            /// Set up 990-1005 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyyzz_xxxx = cbuffer.data(kg_geom_01_off + 990 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xxxy = cbuffer.data(kg_geom_01_off + 991 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xxxz = cbuffer.data(kg_geom_01_off + 992 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xxyy = cbuffer.data(kg_geom_01_off + 993 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xxyz = cbuffer.data(kg_geom_01_off + 994 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xxzz = cbuffer.data(kg_geom_01_off + 995 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xyyy = cbuffer.data(kg_geom_01_off + 996 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xyyz = cbuffer.data(kg_geom_01_off + 997 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xyzz = cbuffer.data(kg_geom_01_off + 998 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_xzzz = cbuffer.data(kg_geom_01_off + 999 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_yyyy = cbuffer.data(kg_geom_01_off + 1000 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_yyyz = cbuffer.data(kg_geom_01_off + 1001 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_yyzz = cbuffer.data(kg_geom_01_off + 1002 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_yzzz = cbuffer.data(kg_geom_01_off + 1003 * ccomps * dcomps);

            auto g_0_y_yyyyyzz_zzzz = cbuffer.data(kg_geom_01_off + 1004 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyyz_xxxx, g_0_y_yyyyyz_xxxxz, g_0_y_yyyyyz_xxxy, g_0_y_yyyyyz_xxxyz, g_0_y_yyyyyz_xxxz, g_0_y_yyyyyz_xxxzz, g_0_y_yyyyyz_xxyy, g_0_y_yyyyyz_xxyyz, g_0_y_yyyyyz_xxyz, g_0_y_yyyyyz_xxyzz, g_0_y_yyyyyz_xxzz, g_0_y_yyyyyz_xxzzz, g_0_y_yyyyyz_xyyy, g_0_y_yyyyyz_xyyyz, g_0_y_yyyyyz_xyyz, g_0_y_yyyyyz_xyyzz, g_0_y_yyyyyz_xyzz, g_0_y_yyyyyz_xyzzz, g_0_y_yyyyyz_xzzz, g_0_y_yyyyyz_xzzzz, g_0_y_yyyyyz_yyyy, g_0_y_yyyyyz_yyyyz, g_0_y_yyyyyz_yyyz, g_0_y_yyyyyz_yyyzz, g_0_y_yyyyyz_yyzz, g_0_y_yyyyyz_yyzzz, g_0_y_yyyyyz_yzzz, g_0_y_yyyyyz_yzzzz, g_0_y_yyyyyz_zzzz, g_0_y_yyyyyz_zzzzz, g_0_y_yyyyyzz_xxxx, g_0_y_yyyyyzz_xxxy, g_0_y_yyyyyzz_xxxz, g_0_y_yyyyyzz_xxyy, g_0_y_yyyyyzz_xxyz, g_0_y_yyyyyzz_xxzz, g_0_y_yyyyyzz_xyyy, g_0_y_yyyyyzz_xyyz, g_0_y_yyyyyzz_xyzz, g_0_y_yyyyyzz_xzzz, g_0_y_yyyyyzz_yyyy, g_0_y_yyyyyzz_yyyz, g_0_y_yyyyyzz_yyzz, g_0_y_yyyyyzz_yzzz, g_0_y_yyyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyyzz_xxxx[k] = -g_0_y_yyyyyz_xxxx[k] * ab_z + g_0_y_yyyyyz_xxxxz[k];

                g_0_y_yyyyyzz_xxxy[k] = -g_0_y_yyyyyz_xxxy[k] * ab_z + g_0_y_yyyyyz_xxxyz[k];

                g_0_y_yyyyyzz_xxxz[k] = -g_0_y_yyyyyz_xxxz[k] * ab_z + g_0_y_yyyyyz_xxxzz[k];

                g_0_y_yyyyyzz_xxyy[k] = -g_0_y_yyyyyz_xxyy[k] * ab_z + g_0_y_yyyyyz_xxyyz[k];

                g_0_y_yyyyyzz_xxyz[k] = -g_0_y_yyyyyz_xxyz[k] * ab_z + g_0_y_yyyyyz_xxyzz[k];

                g_0_y_yyyyyzz_xxzz[k] = -g_0_y_yyyyyz_xxzz[k] * ab_z + g_0_y_yyyyyz_xxzzz[k];

                g_0_y_yyyyyzz_xyyy[k] = -g_0_y_yyyyyz_xyyy[k] * ab_z + g_0_y_yyyyyz_xyyyz[k];

                g_0_y_yyyyyzz_xyyz[k] = -g_0_y_yyyyyz_xyyz[k] * ab_z + g_0_y_yyyyyz_xyyzz[k];

                g_0_y_yyyyyzz_xyzz[k] = -g_0_y_yyyyyz_xyzz[k] * ab_z + g_0_y_yyyyyz_xyzzz[k];

                g_0_y_yyyyyzz_xzzz[k] = -g_0_y_yyyyyz_xzzz[k] * ab_z + g_0_y_yyyyyz_xzzzz[k];

                g_0_y_yyyyyzz_yyyy[k] = -g_0_y_yyyyyz_yyyy[k] * ab_z + g_0_y_yyyyyz_yyyyz[k];

                g_0_y_yyyyyzz_yyyz[k] = -g_0_y_yyyyyz_yyyz[k] * ab_z + g_0_y_yyyyyz_yyyzz[k];

                g_0_y_yyyyyzz_yyzz[k] = -g_0_y_yyyyyz_yyzz[k] * ab_z + g_0_y_yyyyyz_yyzzz[k];

                g_0_y_yyyyyzz_yzzz[k] = -g_0_y_yyyyyz_yzzz[k] * ab_z + g_0_y_yyyyyz_yzzzz[k];

                g_0_y_yyyyyzz_zzzz[k] = -g_0_y_yyyyyz_zzzz[k] * ab_z + g_0_y_yyyyyz_zzzzz[k];
            }

            /// Set up 1005-1020 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyyzzz_xxxx = cbuffer.data(kg_geom_01_off + 1005 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xxxy = cbuffer.data(kg_geom_01_off + 1006 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xxxz = cbuffer.data(kg_geom_01_off + 1007 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xxyy = cbuffer.data(kg_geom_01_off + 1008 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xxyz = cbuffer.data(kg_geom_01_off + 1009 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xxzz = cbuffer.data(kg_geom_01_off + 1010 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xyyy = cbuffer.data(kg_geom_01_off + 1011 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xyyz = cbuffer.data(kg_geom_01_off + 1012 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xyzz = cbuffer.data(kg_geom_01_off + 1013 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_xzzz = cbuffer.data(kg_geom_01_off + 1014 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_yyyy = cbuffer.data(kg_geom_01_off + 1015 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_yyyz = cbuffer.data(kg_geom_01_off + 1016 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_yyzz = cbuffer.data(kg_geom_01_off + 1017 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_yzzz = cbuffer.data(kg_geom_01_off + 1018 * ccomps * dcomps);

            auto g_0_y_yyyyzzz_zzzz = cbuffer.data(kg_geom_01_off + 1019 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyyzz_xxxx, g_0_y_yyyyzz_xxxxz, g_0_y_yyyyzz_xxxy, g_0_y_yyyyzz_xxxyz, g_0_y_yyyyzz_xxxz, g_0_y_yyyyzz_xxxzz, g_0_y_yyyyzz_xxyy, g_0_y_yyyyzz_xxyyz, g_0_y_yyyyzz_xxyz, g_0_y_yyyyzz_xxyzz, g_0_y_yyyyzz_xxzz, g_0_y_yyyyzz_xxzzz, g_0_y_yyyyzz_xyyy, g_0_y_yyyyzz_xyyyz, g_0_y_yyyyzz_xyyz, g_0_y_yyyyzz_xyyzz, g_0_y_yyyyzz_xyzz, g_0_y_yyyyzz_xyzzz, g_0_y_yyyyzz_xzzz, g_0_y_yyyyzz_xzzzz, g_0_y_yyyyzz_yyyy, g_0_y_yyyyzz_yyyyz, g_0_y_yyyyzz_yyyz, g_0_y_yyyyzz_yyyzz, g_0_y_yyyyzz_yyzz, g_0_y_yyyyzz_yyzzz, g_0_y_yyyyzz_yzzz, g_0_y_yyyyzz_yzzzz, g_0_y_yyyyzz_zzzz, g_0_y_yyyyzz_zzzzz, g_0_y_yyyyzzz_xxxx, g_0_y_yyyyzzz_xxxy, g_0_y_yyyyzzz_xxxz, g_0_y_yyyyzzz_xxyy, g_0_y_yyyyzzz_xxyz, g_0_y_yyyyzzz_xxzz, g_0_y_yyyyzzz_xyyy, g_0_y_yyyyzzz_xyyz, g_0_y_yyyyzzz_xyzz, g_0_y_yyyyzzz_xzzz, g_0_y_yyyyzzz_yyyy, g_0_y_yyyyzzz_yyyz, g_0_y_yyyyzzz_yyzz, g_0_y_yyyyzzz_yzzz, g_0_y_yyyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyyzzz_xxxx[k] = -g_0_y_yyyyzz_xxxx[k] * ab_z + g_0_y_yyyyzz_xxxxz[k];

                g_0_y_yyyyzzz_xxxy[k] = -g_0_y_yyyyzz_xxxy[k] * ab_z + g_0_y_yyyyzz_xxxyz[k];

                g_0_y_yyyyzzz_xxxz[k] = -g_0_y_yyyyzz_xxxz[k] * ab_z + g_0_y_yyyyzz_xxxzz[k];

                g_0_y_yyyyzzz_xxyy[k] = -g_0_y_yyyyzz_xxyy[k] * ab_z + g_0_y_yyyyzz_xxyyz[k];

                g_0_y_yyyyzzz_xxyz[k] = -g_0_y_yyyyzz_xxyz[k] * ab_z + g_0_y_yyyyzz_xxyzz[k];

                g_0_y_yyyyzzz_xxzz[k] = -g_0_y_yyyyzz_xxzz[k] * ab_z + g_0_y_yyyyzz_xxzzz[k];

                g_0_y_yyyyzzz_xyyy[k] = -g_0_y_yyyyzz_xyyy[k] * ab_z + g_0_y_yyyyzz_xyyyz[k];

                g_0_y_yyyyzzz_xyyz[k] = -g_0_y_yyyyzz_xyyz[k] * ab_z + g_0_y_yyyyzz_xyyzz[k];

                g_0_y_yyyyzzz_xyzz[k] = -g_0_y_yyyyzz_xyzz[k] * ab_z + g_0_y_yyyyzz_xyzzz[k];

                g_0_y_yyyyzzz_xzzz[k] = -g_0_y_yyyyzz_xzzz[k] * ab_z + g_0_y_yyyyzz_xzzzz[k];

                g_0_y_yyyyzzz_yyyy[k] = -g_0_y_yyyyzz_yyyy[k] * ab_z + g_0_y_yyyyzz_yyyyz[k];

                g_0_y_yyyyzzz_yyyz[k] = -g_0_y_yyyyzz_yyyz[k] * ab_z + g_0_y_yyyyzz_yyyzz[k];

                g_0_y_yyyyzzz_yyzz[k] = -g_0_y_yyyyzz_yyzz[k] * ab_z + g_0_y_yyyyzz_yyzzz[k];

                g_0_y_yyyyzzz_yzzz[k] = -g_0_y_yyyyzz_yzzz[k] * ab_z + g_0_y_yyyyzz_yzzzz[k];

                g_0_y_yyyyzzz_zzzz[k] = -g_0_y_yyyyzz_zzzz[k] * ab_z + g_0_y_yyyyzz_zzzzz[k];
            }

            /// Set up 1020-1035 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyyzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1020 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1021 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1022 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1023 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1024 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1025 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1026 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1027 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1028 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1029 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1030 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1031 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1032 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1033 * ccomps * dcomps);

            auto g_0_y_yyyzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1034 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyyzzz_xxxx, g_0_y_yyyzzz_xxxxz, g_0_y_yyyzzz_xxxy, g_0_y_yyyzzz_xxxyz, g_0_y_yyyzzz_xxxz, g_0_y_yyyzzz_xxxzz, g_0_y_yyyzzz_xxyy, g_0_y_yyyzzz_xxyyz, g_0_y_yyyzzz_xxyz, g_0_y_yyyzzz_xxyzz, g_0_y_yyyzzz_xxzz, g_0_y_yyyzzz_xxzzz, g_0_y_yyyzzz_xyyy, g_0_y_yyyzzz_xyyyz, g_0_y_yyyzzz_xyyz, g_0_y_yyyzzz_xyyzz, g_0_y_yyyzzz_xyzz, g_0_y_yyyzzz_xyzzz, g_0_y_yyyzzz_xzzz, g_0_y_yyyzzz_xzzzz, g_0_y_yyyzzz_yyyy, g_0_y_yyyzzz_yyyyz, g_0_y_yyyzzz_yyyz, g_0_y_yyyzzz_yyyzz, g_0_y_yyyzzz_yyzz, g_0_y_yyyzzz_yyzzz, g_0_y_yyyzzz_yzzz, g_0_y_yyyzzz_yzzzz, g_0_y_yyyzzz_zzzz, g_0_y_yyyzzz_zzzzz, g_0_y_yyyzzzz_xxxx, g_0_y_yyyzzzz_xxxy, g_0_y_yyyzzzz_xxxz, g_0_y_yyyzzzz_xxyy, g_0_y_yyyzzzz_xxyz, g_0_y_yyyzzzz_xxzz, g_0_y_yyyzzzz_xyyy, g_0_y_yyyzzzz_xyyz, g_0_y_yyyzzzz_xyzz, g_0_y_yyyzzzz_xzzz, g_0_y_yyyzzzz_yyyy, g_0_y_yyyzzzz_yyyz, g_0_y_yyyzzzz_yyzz, g_0_y_yyyzzzz_yzzz, g_0_y_yyyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyyzzzz_xxxx[k] = -g_0_y_yyyzzz_xxxx[k] * ab_z + g_0_y_yyyzzz_xxxxz[k];

                g_0_y_yyyzzzz_xxxy[k] = -g_0_y_yyyzzz_xxxy[k] * ab_z + g_0_y_yyyzzz_xxxyz[k];

                g_0_y_yyyzzzz_xxxz[k] = -g_0_y_yyyzzz_xxxz[k] * ab_z + g_0_y_yyyzzz_xxxzz[k];

                g_0_y_yyyzzzz_xxyy[k] = -g_0_y_yyyzzz_xxyy[k] * ab_z + g_0_y_yyyzzz_xxyyz[k];

                g_0_y_yyyzzzz_xxyz[k] = -g_0_y_yyyzzz_xxyz[k] * ab_z + g_0_y_yyyzzz_xxyzz[k];

                g_0_y_yyyzzzz_xxzz[k] = -g_0_y_yyyzzz_xxzz[k] * ab_z + g_0_y_yyyzzz_xxzzz[k];

                g_0_y_yyyzzzz_xyyy[k] = -g_0_y_yyyzzz_xyyy[k] * ab_z + g_0_y_yyyzzz_xyyyz[k];

                g_0_y_yyyzzzz_xyyz[k] = -g_0_y_yyyzzz_xyyz[k] * ab_z + g_0_y_yyyzzz_xyyzz[k];

                g_0_y_yyyzzzz_xyzz[k] = -g_0_y_yyyzzz_xyzz[k] * ab_z + g_0_y_yyyzzz_xyzzz[k];

                g_0_y_yyyzzzz_xzzz[k] = -g_0_y_yyyzzz_xzzz[k] * ab_z + g_0_y_yyyzzz_xzzzz[k];

                g_0_y_yyyzzzz_yyyy[k] = -g_0_y_yyyzzz_yyyy[k] * ab_z + g_0_y_yyyzzz_yyyyz[k];

                g_0_y_yyyzzzz_yyyz[k] = -g_0_y_yyyzzz_yyyz[k] * ab_z + g_0_y_yyyzzz_yyyzz[k];

                g_0_y_yyyzzzz_yyzz[k] = -g_0_y_yyyzzz_yyzz[k] * ab_z + g_0_y_yyyzzz_yyzzz[k];

                g_0_y_yyyzzzz_yzzz[k] = -g_0_y_yyyzzz_yzzz[k] * ab_z + g_0_y_yyyzzz_yzzzz[k];

                g_0_y_yyyzzzz_zzzz[k] = -g_0_y_yyyzzz_zzzz[k] * ab_z + g_0_y_yyyzzz_zzzzz[k];
            }

            /// Set up 1035-1050 components of targeted buffer : cbuffer.data(

            auto g_0_y_yyzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1035 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1036 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1037 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1038 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1039 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1040 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1041 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1042 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1043 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1044 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1045 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1046 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1047 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1048 * ccomps * dcomps);

            auto g_0_y_yyzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1049 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yyzzzz_xxxx, g_0_y_yyzzzz_xxxxz, g_0_y_yyzzzz_xxxy, g_0_y_yyzzzz_xxxyz, g_0_y_yyzzzz_xxxz, g_0_y_yyzzzz_xxxzz, g_0_y_yyzzzz_xxyy, g_0_y_yyzzzz_xxyyz, g_0_y_yyzzzz_xxyz, g_0_y_yyzzzz_xxyzz, g_0_y_yyzzzz_xxzz, g_0_y_yyzzzz_xxzzz, g_0_y_yyzzzz_xyyy, g_0_y_yyzzzz_xyyyz, g_0_y_yyzzzz_xyyz, g_0_y_yyzzzz_xyyzz, g_0_y_yyzzzz_xyzz, g_0_y_yyzzzz_xyzzz, g_0_y_yyzzzz_xzzz, g_0_y_yyzzzz_xzzzz, g_0_y_yyzzzz_yyyy, g_0_y_yyzzzz_yyyyz, g_0_y_yyzzzz_yyyz, g_0_y_yyzzzz_yyyzz, g_0_y_yyzzzz_yyzz, g_0_y_yyzzzz_yyzzz, g_0_y_yyzzzz_yzzz, g_0_y_yyzzzz_yzzzz, g_0_y_yyzzzz_zzzz, g_0_y_yyzzzz_zzzzz, g_0_y_yyzzzzz_xxxx, g_0_y_yyzzzzz_xxxy, g_0_y_yyzzzzz_xxxz, g_0_y_yyzzzzz_xxyy, g_0_y_yyzzzzz_xxyz, g_0_y_yyzzzzz_xxzz, g_0_y_yyzzzzz_xyyy, g_0_y_yyzzzzz_xyyz, g_0_y_yyzzzzz_xyzz, g_0_y_yyzzzzz_xzzz, g_0_y_yyzzzzz_yyyy, g_0_y_yyzzzzz_yyyz, g_0_y_yyzzzzz_yyzz, g_0_y_yyzzzzz_yzzz, g_0_y_yyzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yyzzzzz_xxxx[k] = -g_0_y_yyzzzz_xxxx[k] * ab_z + g_0_y_yyzzzz_xxxxz[k];

                g_0_y_yyzzzzz_xxxy[k] = -g_0_y_yyzzzz_xxxy[k] * ab_z + g_0_y_yyzzzz_xxxyz[k];

                g_0_y_yyzzzzz_xxxz[k] = -g_0_y_yyzzzz_xxxz[k] * ab_z + g_0_y_yyzzzz_xxxzz[k];

                g_0_y_yyzzzzz_xxyy[k] = -g_0_y_yyzzzz_xxyy[k] * ab_z + g_0_y_yyzzzz_xxyyz[k];

                g_0_y_yyzzzzz_xxyz[k] = -g_0_y_yyzzzz_xxyz[k] * ab_z + g_0_y_yyzzzz_xxyzz[k];

                g_0_y_yyzzzzz_xxzz[k] = -g_0_y_yyzzzz_xxzz[k] * ab_z + g_0_y_yyzzzz_xxzzz[k];

                g_0_y_yyzzzzz_xyyy[k] = -g_0_y_yyzzzz_xyyy[k] * ab_z + g_0_y_yyzzzz_xyyyz[k];

                g_0_y_yyzzzzz_xyyz[k] = -g_0_y_yyzzzz_xyyz[k] * ab_z + g_0_y_yyzzzz_xyyzz[k];

                g_0_y_yyzzzzz_xyzz[k] = -g_0_y_yyzzzz_xyzz[k] * ab_z + g_0_y_yyzzzz_xyzzz[k];

                g_0_y_yyzzzzz_xzzz[k] = -g_0_y_yyzzzz_xzzz[k] * ab_z + g_0_y_yyzzzz_xzzzz[k];

                g_0_y_yyzzzzz_yyyy[k] = -g_0_y_yyzzzz_yyyy[k] * ab_z + g_0_y_yyzzzz_yyyyz[k];

                g_0_y_yyzzzzz_yyyz[k] = -g_0_y_yyzzzz_yyyz[k] * ab_z + g_0_y_yyzzzz_yyyzz[k];

                g_0_y_yyzzzzz_yyzz[k] = -g_0_y_yyzzzz_yyzz[k] * ab_z + g_0_y_yyzzzz_yyzzz[k];

                g_0_y_yyzzzzz_yzzz[k] = -g_0_y_yyzzzz_yzzz[k] * ab_z + g_0_y_yyzzzz_yzzzz[k];

                g_0_y_yyzzzzz_zzzz[k] = -g_0_y_yyzzzz_zzzz[k] * ab_z + g_0_y_yyzzzz_zzzzz[k];
            }

            /// Set up 1050-1065 components of targeted buffer : cbuffer.data(

            auto g_0_y_yzzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1050 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1051 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1052 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1053 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1054 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1055 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1056 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1057 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1058 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1059 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1060 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1061 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1062 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1063 * ccomps * dcomps);

            auto g_0_y_yzzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1064 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_yzzzzz_xxxx, g_0_y_yzzzzz_xxxxz, g_0_y_yzzzzz_xxxy, g_0_y_yzzzzz_xxxyz, g_0_y_yzzzzz_xxxz, g_0_y_yzzzzz_xxxzz, g_0_y_yzzzzz_xxyy, g_0_y_yzzzzz_xxyyz, g_0_y_yzzzzz_xxyz, g_0_y_yzzzzz_xxyzz, g_0_y_yzzzzz_xxzz, g_0_y_yzzzzz_xxzzz, g_0_y_yzzzzz_xyyy, g_0_y_yzzzzz_xyyyz, g_0_y_yzzzzz_xyyz, g_0_y_yzzzzz_xyyzz, g_0_y_yzzzzz_xyzz, g_0_y_yzzzzz_xyzzz, g_0_y_yzzzzz_xzzz, g_0_y_yzzzzz_xzzzz, g_0_y_yzzzzz_yyyy, g_0_y_yzzzzz_yyyyz, g_0_y_yzzzzz_yyyz, g_0_y_yzzzzz_yyyzz, g_0_y_yzzzzz_yyzz, g_0_y_yzzzzz_yyzzz, g_0_y_yzzzzz_yzzz, g_0_y_yzzzzz_yzzzz, g_0_y_yzzzzz_zzzz, g_0_y_yzzzzz_zzzzz, g_0_y_yzzzzzz_xxxx, g_0_y_yzzzzzz_xxxy, g_0_y_yzzzzzz_xxxz, g_0_y_yzzzzzz_xxyy, g_0_y_yzzzzzz_xxyz, g_0_y_yzzzzzz_xxzz, g_0_y_yzzzzzz_xyyy, g_0_y_yzzzzzz_xyyz, g_0_y_yzzzzzz_xyzz, g_0_y_yzzzzzz_xzzz, g_0_y_yzzzzzz_yyyy, g_0_y_yzzzzzz_yyyz, g_0_y_yzzzzzz_yyzz, g_0_y_yzzzzzz_yzzz, g_0_y_yzzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_yzzzzzz_xxxx[k] = -g_0_y_yzzzzz_xxxx[k] * ab_z + g_0_y_yzzzzz_xxxxz[k];

                g_0_y_yzzzzzz_xxxy[k] = -g_0_y_yzzzzz_xxxy[k] * ab_z + g_0_y_yzzzzz_xxxyz[k];

                g_0_y_yzzzzzz_xxxz[k] = -g_0_y_yzzzzz_xxxz[k] * ab_z + g_0_y_yzzzzz_xxxzz[k];

                g_0_y_yzzzzzz_xxyy[k] = -g_0_y_yzzzzz_xxyy[k] * ab_z + g_0_y_yzzzzz_xxyyz[k];

                g_0_y_yzzzzzz_xxyz[k] = -g_0_y_yzzzzz_xxyz[k] * ab_z + g_0_y_yzzzzz_xxyzz[k];

                g_0_y_yzzzzzz_xxzz[k] = -g_0_y_yzzzzz_xxzz[k] * ab_z + g_0_y_yzzzzz_xxzzz[k];

                g_0_y_yzzzzzz_xyyy[k] = -g_0_y_yzzzzz_xyyy[k] * ab_z + g_0_y_yzzzzz_xyyyz[k];

                g_0_y_yzzzzzz_xyyz[k] = -g_0_y_yzzzzz_xyyz[k] * ab_z + g_0_y_yzzzzz_xyyzz[k];

                g_0_y_yzzzzzz_xyzz[k] = -g_0_y_yzzzzz_xyzz[k] * ab_z + g_0_y_yzzzzz_xyzzz[k];

                g_0_y_yzzzzzz_xzzz[k] = -g_0_y_yzzzzz_xzzz[k] * ab_z + g_0_y_yzzzzz_xzzzz[k];

                g_0_y_yzzzzzz_yyyy[k] = -g_0_y_yzzzzz_yyyy[k] * ab_z + g_0_y_yzzzzz_yyyyz[k];

                g_0_y_yzzzzzz_yyyz[k] = -g_0_y_yzzzzz_yyyz[k] * ab_z + g_0_y_yzzzzz_yyyzz[k];

                g_0_y_yzzzzzz_yyzz[k] = -g_0_y_yzzzzz_yyzz[k] * ab_z + g_0_y_yzzzzz_yyzzz[k];

                g_0_y_yzzzzzz_yzzz[k] = -g_0_y_yzzzzz_yzzz[k] * ab_z + g_0_y_yzzzzz_yzzzz[k];

                g_0_y_yzzzzzz_zzzz[k] = -g_0_y_yzzzzz_zzzz[k] * ab_z + g_0_y_yzzzzz_zzzzz[k];
            }

            /// Set up 1065-1080 components of targeted buffer : cbuffer.data(

            auto g_0_y_zzzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1065 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1066 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1067 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1068 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1069 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1070 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1071 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1072 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1073 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1074 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1075 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1076 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1077 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1078 * ccomps * dcomps);

            auto g_0_y_zzzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1079 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_y_zzzzzz_xxxx, g_0_y_zzzzzz_xxxxz, g_0_y_zzzzzz_xxxy, g_0_y_zzzzzz_xxxyz, g_0_y_zzzzzz_xxxz, g_0_y_zzzzzz_xxxzz, g_0_y_zzzzzz_xxyy, g_0_y_zzzzzz_xxyyz, g_0_y_zzzzzz_xxyz, g_0_y_zzzzzz_xxyzz, g_0_y_zzzzzz_xxzz, g_0_y_zzzzzz_xxzzz, g_0_y_zzzzzz_xyyy, g_0_y_zzzzzz_xyyyz, g_0_y_zzzzzz_xyyz, g_0_y_zzzzzz_xyyzz, g_0_y_zzzzzz_xyzz, g_0_y_zzzzzz_xyzzz, g_0_y_zzzzzz_xzzz, g_0_y_zzzzzz_xzzzz, g_0_y_zzzzzz_yyyy, g_0_y_zzzzzz_yyyyz, g_0_y_zzzzzz_yyyz, g_0_y_zzzzzz_yyyzz, g_0_y_zzzzzz_yyzz, g_0_y_zzzzzz_yyzzz, g_0_y_zzzzzz_yzzz, g_0_y_zzzzzz_yzzzz, g_0_y_zzzzzz_zzzz, g_0_y_zzzzzz_zzzzz, g_0_y_zzzzzzz_xxxx, g_0_y_zzzzzzz_xxxy, g_0_y_zzzzzzz_xxxz, g_0_y_zzzzzzz_xxyy, g_0_y_zzzzzzz_xxyz, g_0_y_zzzzzzz_xxzz, g_0_y_zzzzzzz_xyyy, g_0_y_zzzzzzz_xyyz, g_0_y_zzzzzzz_xyzz, g_0_y_zzzzzzz_xzzz, g_0_y_zzzzzzz_yyyy, g_0_y_zzzzzzz_yyyz, g_0_y_zzzzzzz_yyzz, g_0_y_zzzzzzz_yzzz, g_0_y_zzzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_y_zzzzzzz_xxxx[k] = -g_0_y_zzzzzz_xxxx[k] * ab_z + g_0_y_zzzzzz_xxxxz[k];

                g_0_y_zzzzzzz_xxxy[k] = -g_0_y_zzzzzz_xxxy[k] * ab_z + g_0_y_zzzzzz_xxxyz[k];

                g_0_y_zzzzzzz_xxxz[k] = -g_0_y_zzzzzz_xxxz[k] * ab_z + g_0_y_zzzzzz_xxxzz[k];

                g_0_y_zzzzzzz_xxyy[k] = -g_0_y_zzzzzz_xxyy[k] * ab_z + g_0_y_zzzzzz_xxyyz[k];

                g_0_y_zzzzzzz_xxyz[k] = -g_0_y_zzzzzz_xxyz[k] * ab_z + g_0_y_zzzzzz_xxyzz[k];

                g_0_y_zzzzzzz_xxzz[k] = -g_0_y_zzzzzz_xxzz[k] * ab_z + g_0_y_zzzzzz_xxzzz[k];

                g_0_y_zzzzzzz_xyyy[k] = -g_0_y_zzzzzz_xyyy[k] * ab_z + g_0_y_zzzzzz_xyyyz[k];

                g_0_y_zzzzzzz_xyyz[k] = -g_0_y_zzzzzz_xyyz[k] * ab_z + g_0_y_zzzzzz_xyyzz[k];

                g_0_y_zzzzzzz_xyzz[k] = -g_0_y_zzzzzz_xyzz[k] * ab_z + g_0_y_zzzzzz_xyzzz[k];

                g_0_y_zzzzzzz_xzzz[k] = -g_0_y_zzzzzz_xzzz[k] * ab_z + g_0_y_zzzzzz_xzzzz[k];

                g_0_y_zzzzzzz_yyyy[k] = -g_0_y_zzzzzz_yyyy[k] * ab_z + g_0_y_zzzzzz_yyyyz[k];

                g_0_y_zzzzzzz_yyyz[k] = -g_0_y_zzzzzz_yyyz[k] * ab_z + g_0_y_zzzzzz_yyyzz[k];

                g_0_y_zzzzzzz_yyzz[k] = -g_0_y_zzzzzz_yyzz[k] * ab_z + g_0_y_zzzzzz_yyzzz[k];

                g_0_y_zzzzzzz_yzzz[k] = -g_0_y_zzzzzz_yzzz[k] * ab_z + g_0_y_zzzzzz_yzzzz[k];

                g_0_y_zzzzzzz_zzzz[k] = -g_0_y_zzzzzz_zzzz[k] * ab_z + g_0_y_zzzzzz_zzzzz[k];
            }

            /// Set up 1080-1095 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxx_xxxx = cbuffer.data(kg_geom_01_off + 1080 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xxxy = cbuffer.data(kg_geom_01_off + 1081 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xxxz = cbuffer.data(kg_geom_01_off + 1082 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xxyy = cbuffer.data(kg_geom_01_off + 1083 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xxyz = cbuffer.data(kg_geom_01_off + 1084 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xxzz = cbuffer.data(kg_geom_01_off + 1085 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xyyy = cbuffer.data(kg_geom_01_off + 1086 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xyyz = cbuffer.data(kg_geom_01_off + 1087 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xyzz = cbuffer.data(kg_geom_01_off + 1088 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_xzzz = cbuffer.data(kg_geom_01_off + 1089 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_yyyy = cbuffer.data(kg_geom_01_off + 1090 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_yyyz = cbuffer.data(kg_geom_01_off + 1091 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_yyzz = cbuffer.data(kg_geom_01_off + 1092 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_yzzz = cbuffer.data(kg_geom_01_off + 1093 * ccomps * dcomps);

            auto g_0_z_xxxxxxx_zzzz = cbuffer.data(kg_geom_01_off + 1094 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxx_xxxx, g_0_z_xxxxxx_xxxxx, g_0_z_xxxxxx_xxxxy, g_0_z_xxxxxx_xxxxz, g_0_z_xxxxxx_xxxy, g_0_z_xxxxxx_xxxyy, g_0_z_xxxxxx_xxxyz, g_0_z_xxxxxx_xxxz, g_0_z_xxxxxx_xxxzz, g_0_z_xxxxxx_xxyy, g_0_z_xxxxxx_xxyyy, g_0_z_xxxxxx_xxyyz, g_0_z_xxxxxx_xxyz, g_0_z_xxxxxx_xxyzz, g_0_z_xxxxxx_xxzz, g_0_z_xxxxxx_xxzzz, g_0_z_xxxxxx_xyyy, g_0_z_xxxxxx_xyyyy, g_0_z_xxxxxx_xyyyz, g_0_z_xxxxxx_xyyz, g_0_z_xxxxxx_xyyzz, g_0_z_xxxxxx_xyzz, g_0_z_xxxxxx_xyzzz, g_0_z_xxxxxx_xzzz, g_0_z_xxxxxx_xzzzz, g_0_z_xxxxxx_yyyy, g_0_z_xxxxxx_yyyz, g_0_z_xxxxxx_yyzz, g_0_z_xxxxxx_yzzz, g_0_z_xxxxxx_zzzz, g_0_z_xxxxxxx_xxxx, g_0_z_xxxxxxx_xxxy, g_0_z_xxxxxxx_xxxz, g_0_z_xxxxxxx_xxyy, g_0_z_xxxxxxx_xxyz, g_0_z_xxxxxxx_xxzz, g_0_z_xxxxxxx_xyyy, g_0_z_xxxxxxx_xyyz, g_0_z_xxxxxxx_xyzz, g_0_z_xxxxxxx_xzzz, g_0_z_xxxxxxx_yyyy, g_0_z_xxxxxxx_yyyz, g_0_z_xxxxxxx_yyzz, g_0_z_xxxxxxx_yzzz, g_0_z_xxxxxxx_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxx_xxxx[k] = -g_0_z_xxxxxx_xxxx[k] * ab_x + g_0_z_xxxxxx_xxxxx[k];

                g_0_z_xxxxxxx_xxxy[k] = -g_0_z_xxxxxx_xxxy[k] * ab_x + g_0_z_xxxxxx_xxxxy[k];

                g_0_z_xxxxxxx_xxxz[k] = -g_0_z_xxxxxx_xxxz[k] * ab_x + g_0_z_xxxxxx_xxxxz[k];

                g_0_z_xxxxxxx_xxyy[k] = -g_0_z_xxxxxx_xxyy[k] * ab_x + g_0_z_xxxxxx_xxxyy[k];

                g_0_z_xxxxxxx_xxyz[k] = -g_0_z_xxxxxx_xxyz[k] * ab_x + g_0_z_xxxxxx_xxxyz[k];

                g_0_z_xxxxxxx_xxzz[k] = -g_0_z_xxxxxx_xxzz[k] * ab_x + g_0_z_xxxxxx_xxxzz[k];

                g_0_z_xxxxxxx_xyyy[k] = -g_0_z_xxxxxx_xyyy[k] * ab_x + g_0_z_xxxxxx_xxyyy[k];

                g_0_z_xxxxxxx_xyyz[k] = -g_0_z_xxxxxx_xyyz[k] * ab_x + g_0_z_xxxxxx_xxyyz[k];

                g_0_z_xxxxxxx_xyzz[k] = -g_0_z_xxxxxx_xyzz[k] * ab_x + g_0_z_xxxxxx_xxyzz[k];

                g_0_z_xxxxxxx_xzzz[k] = -g_0_z_xxxxxx_xzzz[k] * ab_x + g_0_z_xxxxxx_xxzzz[k];

                g_0_z_xxxxxxx_yyyy[k] = -g_0_z_xxxxxx_yyyy[k] * ab_x + g_0_z_xxxxxx_xyyyy[k];

                g_0_z_xxxxxxx_yyyz[k] = -g_0_z_xxxxxx_yyyz[k] * ab_x + g_0_z_xxxxxx_xyyyz[k];

                g_0_z_xxxxxxx_yyzz[k] = -g_0_z_xxxxxx_yyzz[k] * ab_x + g_0_z_xxxxxx_xyyzz[k];

                g_0_z_xxxxxxx_yzzz[k] = -g_0_z_xxxxxx_yzzz[k] * ab_x + g_0_z_xxxxxx_xyzzz[k];

                g_0_z_xxxxxxx_zzzz[k] = -g_0_z_xxxxxx_zzzz[k] * ab_x + g_0_z_xxxxxx_xzzzz[k];
            }

            /// Set up 1095-1110 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxy_xxxx = cbuffer.data(kg_geom_01_off + 1095 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xxxy = cbuffer.data(kg_geom_01_off + 1096 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xxxz = cbuffer.data(kg_geom_01_off + 1097 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xxyy = cbuffer.data(kg_geom_01_off + 1098 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xxyz = cbuffer.data(kg_geom_01_off + 1099 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xxzz = cbuffer.data(kg_geom_01_off + 1100 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xyyy = cbuffer.data(kg_geom_01_off + 1101 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xyyz = cbuffer.data(kg_geom_01_off + 1102 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xyzz = cbuffer.data(kg_geom_01_off + 1103 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_xzzz = cbuffer.data(kg_geom_01_off + 1104 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_yyyy = cbuffer.data(kg_geom_01_off + 1105 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_yyyz = cbuffer.data(kg_geom_01_off + 1106 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_yyzz = cbuffer.data(kg_geom_01_off + 1107 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_yzzz = cbuffer.data(kg_geom_01_off + 1108 * ccomps * dcomps);

            auto g_0_z_xxxxxxy_zzzz = cbuffer.data(kg_geom_01_off + 1109 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxxy_xxxx, g_0_z_xxxxxxy_xxxy, g_0_z_xxxxxxy_xxxz, g_0_z_xxxxxxy_xxyy, g_0_z_xxxxxxy_xxyz, g_0_z_xxxxxxy_xxzz, g_0_z_xxxxxxy_xyyy, g_0_z_xxxxxxy_xyyz, g_0_z_xxxxxxy_xyzz, g_0_z_xxxxxxy_xzzz, g_0_z_xxxxxxy_yyyy, g_0_z_xxxxxxy_yyyz, g_0_z_xxxxxxy_yyzz, g_0_z_xxxxxxy_yzzz, g_0_z_xxxxxxy_zzzz, g_0_z_xxxxxy_xxxx, g_0_z_xxxxxy_xxxxx, g_0_z_xxxxxy_xxxxy, g_0_z_xxxxxy_xxxxz, g_0_z_xxxxxy_xxxy, g_0_z_xxxxxy_xxxyy, g_0_z_xxxxxy_xxxyz, g_0_z_xxxxxy_xxxz, g_0_z_xxxxxy_xxxzz, g_0_z_xxxxxy_xxyy, g_0_z_xxxxxy_xxyyy, g_0_z_xxxxxy_xxyyz, g_0_z_xxxxxy_xxyz, g_0_z_xxxxxy_xxyzz, g_0_z_xxxxxy_xxzz, g_0_z_xxxxxy_xxzzz, g_0_z_xxxxxy_xyyy, g_0_z_xxxxxy_xyyyy, g_0_z_xxxxxy_xyyyz, g_0_z_xxxxxy_xyyz, g_0_z_xxxxxy_xyyzz, g_0_z_xxxxxy_xyzz, g_0_z_xxxxxy_xyzzz, g_0_z_xxxxxy_xzzz, g_0_z_xxxxxy_xzzzz, g_0_z_xxxxxy_yyyy, g_0_z_xxxxxy_yyyz, g_0_z_xxxxxy_yyzz, g_0_z_xxxxxy_yzzz, g_0_z_xxxxxy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxy_xxxx[k] = -g_0_z_xxxxxy_xxxx[k] * ab_x + g_0_z_xxxxxy_xxxxx[k];

                g_0_z_xxxxxxy_xxxy[k] = -g_0_z_xxxxxy_xxxy[k] * ab_x + g_0_z_xxxxxy_xxxxy[k];

                g_0_z_xxxxxxy_xxxz[k] = -g_0_z_xxxxxy_xxxz[k] * ab_x + g_0_z_xxxxxy_xxxxz[k];

                g_0_z_xxxxxxy_xxyy[k] = -g_0_z_xxxxxy_xxyy[k] * ab_x + g_0_z_xxxxxy_xxxyy[k];

                g_0_z_xxxxxxy_xxyz[k] = -g_0_z_xxxxxy_xxyz[k] * ab_x + g_0_z_xxxxxy_xxxyz[k];

                g_0_z_xxxxxxy_xxzz[k] = -g_0_z_xxxxxy_xxzz[k] * ab_x + g_0_z_xxxxxy_xxxzz[k];

                g_0_z_xxxxxxy_xyyy[k] = -g_0_z_xxxxxy_xyyy[k] * ab_x + g_0_z_xxxxxy_xxyyy[k];

                g_0_z_xxxxxxy_xyyz[k] = -g_0_z_xxxxxy_xyyz[k] * ab_x + g_0_z_xxxxxy_xxyyz[k];

                g_0_z_xxxxxxy_xyzz[k] = -g_0_z_xxxxxy_xyzz[k] * ab_x + g_0_z_xxxxxy_xxyzz[k];

                g_0_z_xxxxxxy_xzzz[k] = -g_0_z_xxxxxy_xzzz[k] * ab_x + g_0_z_xxxxxy_xxzzz[k];

                g_0_z_xxxxxxy_yyyy[k] = -g_0_z_xxxxxy_yyyy[k] * ab_x + g_0_z_xxxxxy_xyyyy[k];

                g_0_z_xxxxxxy_yyyz[k] = -g_0_z_xxxxxy_yyyz[k] * ab_x + g_0_z_xxxxxy_xyyyz[k];

                g_0_z_xxxxxxy_yyzz[k] = -g_0_z_xxxxxy_yyzz[k] * ab_x + g_0_z_xxxxxy_xyyzz[k];

                g_0_z_xxxxxxy_yzzz[k] = -g_0_z_xxxxxy_yzzz[k] * ab_x + g_0_z_xxxxxy_xyzzz[k];

                g_0_z_xxxxxxy_zzzz[k] = -g_0_z_xxxxxy_zzzz[k] * ab_x + g_0_z_xxxxxy_xzzzz[k];
            }

            /// Set up 1110-1125 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxxz_xxxx = cbuffer.data(kg_geom_01_off + 1110 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xxxy = cbuffer.data(kg_geom_01_off + 1111 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xxxz = cbuffer.data(kg_geom_01_off + 1112 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xxyy = cbuffer.data(kg_geom_01_off + 1113 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xxyz = cbuffer.data(kg_geom_01_off + 1114 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xxzz = cbuffer.data(kg_geom_01_off + 1115 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xyyy = cbuffer.data(kg_geom_01_off + 1116 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xyyz = cbuffer.data(kg_geom_01_off + 1117 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xyzz = cbuffer.data(kg_geom_01_off + 1118 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_xzzz = cbuffer.data(kg_geom_01_off + 1119 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_yyyy = cbuffer.data(kg_geom_01_off + 1120 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_yyyz = cbuffer.data(kg_geom_01_off + 1121 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_yyzz = cbuffer.data(kg_geom_01_off + 1122 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_yzzz = cbuffer.data(kg_geom_01_off + 1123 * ccomps * dcomps);

            auto g_0_z_xxxxxxz_zzzz = cbuffer.data(kg_geom_01_off + 1124 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxxz_xxxx, g_0_z_xxxxxxz_xxxy, g_0_z_xxxxxxz_xxxz, g_0_z_xxxxxxz_xxyy, g_0_z_xxxxxxz_xxyz, g_0_z_xxxxxxz_xxzz, g_0_z_xxxxxxz_xyyy, g_0_z_xxxxxxz_xyyz, g_0_z_xxxxxxz_xyzz, g_0_z_xxxxxxz_xzzz, g_0_z_xxxxxxz_yyyy, g_0_z_xxxxxxz_yyyz, g_0_z_xxxxxxz_yyzz, g_0_z_xxxxxxz_yzzz, g_0_z_xxxxxxz_zzzz, g_0_z_xxxxxz_xxxx, g_0_z_xxxxxz_xxxxx, g_0_z_xxxxxz_xxxxy, g_0_z_xxxxxz_xxxxz, g_0_z_xxxxxz_xxxy, g_0_z_xxxxxz_xxxyy, g_0_z_xxxxxz_xxxyz, g_0_z_xxxxxz_xxxz, g_0_z_xxxxxz_xxxzz, g_0_z_xxxxxz_xxyy, g_0_z_xxxxxz_xxyyy, g_0_z_xxxxxz_xxyyz, g_0_z_xxxxxz_xxyz, g_0_z_xxxxxz_xxyzz, g_0_z_xxxxxz_xxzz, g_0_z_xxxxxz_xxzzz, g_0_z_xxxxxz_xyyy, g_0_z_xxxxxz_xyyyy, g_0_z_xxxxxz_xyyyz, g_0_z_xxxxxz_xyyz, g_0_z_xxxxxz_xyyzz, g_0_z_xxxxxz_xyzz, g_0_z_xxxxxz_xyzzz, g_0_z_xxxxxz_xzzz, g_0_z_xxxxxz_xzzzz, g_0_z_xxxxxz_yyyy, g_0_z_xxxxxz_yyyz, g_0_z_xxxxxz_yyzz, g_0_z_xxxxxz_yzzz, g_0_z_xxxxxz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxxz_xxxx[k] = -g_0_z_xxxxxz_xxxx[k] * ab_x + g_0_z_xxxxxz_xxxxx[k];

                g_0_z_xxxxxxz_xxxy[k] = -g_0_z_xxxxxz_xxxy[k] * ab_x + g_0_z_xxxxxz_xxxxy[k];

                g_0_z_xxxxxxz_xxxz[k] = -g_0_z_xxxxxz_xxxz[k] * ab_x + g_0_z_xxxxxz_xxxxz[k];

                g_0_z_xxxxxxz_xxyy[k] = -g_0_z_xxxxxz_xxyy[k] * ab_x + g_0_z_xxxxxz_xxxyy[k];

                g_0_z_xxxxxxz_xxyz[k] = -g_0_z_xxxxxz_xxyz[k] * ab_x + g_0_z_xxxxxz_xxxyz[k];

                g_0_z_xxxxxxz_xxzz[k] = -g_0_z_xxxxxz_xxzz[k] * ab_x + g_0_z_xxxxxz_xxxzz[k];

                g_0_z_xxxxxxz_xyyy[k] = -g_0_z_xxxxxz_xyyy[k] * ab_x + g_0_z_xxxxxz_xxyyy[k];

                g_0_z_xxxxxxz_xyyz[k] = -g_0_z_xxxxxz_xyyz[k] * ab_x + g_0_z_xxxxxz_xxyyz[k];

                g_0_z_xxxxxxz_xyzz[k] = -g_0_z_xxxxxz_xyzz[k] * ab_x + g_0_z_xxxxxz_xxyzz[k];

                g_0_z_xxxxxxz_xzzz[k] = -g_0_z_xxxxxz_xzzz[k] * ab_x + g_0_z_xxxxxz_xxzzz[k];

                g_0_z_xxxxxxz_yyyy[k] = -g_0_z_xxxxxz_yyyy[k] * ab_x + g_0_z_xxxxxz_xyyyy[k];

                g_0_z_xxxxxxz_yyyz[k] = -g_0_z_xxxxxz_yyyz[k] * ab_x + g_0_z_xxxxxz_xyyyz[k];

                g_0_z_xxxxxxz_yyzz[k] = -g_0_z_xxxxxz_yyzz[k] * ab_x + g_0_z_xxxxxz_xyyzz[k];

                g_0_z_xxxxxxz_yzzz[k] = -g_0_z_xxxxxz_yzzz[k] * ab_x + g_0_z_xxxxxz_xyzzz[k];

                g_0_z_xxxxxxz_zzzz[k] = -g_0_z_xxxxxz_zzzz[k] * ab_x + g_0_z_xxxxxz_xzzzz[k];
            }

            /// Set up 1125-1140 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxyy_xxxx = cbuffer.data(kg_geom_01_off + 1125 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xxxy = cbuffer.data(kg_geom_01_off + 1126 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xxxz = cbuffer.data(kg_geom_01_off + 1127 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xxyy = cbuffer.data(kg_geom_01_off + 1128 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xxyz = cbuffer.data(kg_geom_01_off + 1129 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xxzz = cbuffer.data(kg_geom_01_off + 1130 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xyyy = cbuffer.data(kg_geom_01_off + 1131 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xyyz = cbuffer.data(kg_geom_01_off + 1132 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xyzz = cbuffer.data(kg_geom_01_off + 1133 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_xzzz = cbuffer.data(kg_geom_01_off + 1134 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_yyyy = cbuffer.data(kg_geom_01_off + 1135 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_yyyz = cbuffer.data(kg_geom_01_off + 1136 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_yyzz = cbuffer.data(kg_geom_01_off + 1137 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_yzzz = cbuffer.data(kg_geom_01_off + 1138 * ccomps * dcomps);

            auto g_0_z_xxxxxyy_zzzz = cbuffer.data(kg_geom_01_off + 1139 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxyy_xxxx, g_0_z_xxxxxyy_xxxy, g_0_z_xxxxxyy_xxxz, g_0_z_xxxxxyy_xxyy, g_0_z_xxxxxyy_xxyz, g_0_z_xxxxxyy_xxzz, g_0_z_xxxxxyy_xyyy, g_0_z_xxxxxyy_xyyz, g_0_z_xxxxxyy_xyzz, g_0_z_xxxxxyy_xzzz, g_0_z_xxxxxyy_yyyy, g_0_z_xxxxxyy_yyyz, g_0_z_xxxxxyy_yyzz, g_0_z_xxxxxyy_yzzz, g_0_z_xxxxxyy_zzzz, g_0_z_xxxxyy_xxxx, g_0_z_xxxxyy_xxxxx, g_0_z_xxxxyy_xxxxy, g_0_z_xxxxyy_xxxxz, g_0_z_xxxxyy_xxxy, g_0_z_xxxxyy_xxxyy, g_0_z_xxxxyy_xxxyz, g_0_z_xxxxyy_xxxz, g_0_z_xxxxyy_xxxzz, g_0_z_xxxxyy_xxyy, g_0_z_xxxxyy_xxyyy, g_0_z_xxxxyy_xxyyz, g_0_z_xxxxyy_xxyz, g_0_z_xxxxyy_xxyzz, g_0_z_xxxxyy_xxzz, g_0_z_xxxxyy_xxzzz, g_0_z_xxxxyy_xyyy, g_0_z_xxxxyy_xyyyy, g_0_z_xxxxyy_xyyyz, g_0_z_xxxxyy_xyyz, g_0_z_xxxxyy_xyyzz, g_0_z_xxxxyy_xyzz, g_0_z_xxxxyy_xyzzz, g_0_z_xxxxyy_xzzz, g_0_z_xxxxyy_xzzzz, g_0_z_xxxxyy_yyyy, g_0_z_xxxxyy_yyyz, g_0_z_xxxxyy_yyzz, g_0_z_xxxxyy_yzzz, g_0_z_xxxxyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxyy_xxxx[k] = -g_0_z_xxxxyy_xxxx[k] * ab_x + g_0_z_xxxxyy_xxxxx[k];

                g_0_z_xxxxxyy_xxxy[k] = -g_0_z_xxxxyy_xxxy[k] * ab_x + g_0_z_xxxxyy_xxxxy[k];

                g_0_z_xxxxxyy_xxxz[k] = -g_0_z_xxxxyy_xxxz[k] * ab_x + g_0_z_xxxxyy_xxxxz[k];

                g_0_z_xxxxxyy_xxyy[k] = -g_0_z_xxxxyy_xxyy[k] * ab_x + g_0_z_xxxxyy_xxxyy[k];

                g_0_z_xxxxxyy_xxyz[k] = -g_0_z_xxxxyy_xxyz[k] * ab_x + g_0_z_xxxxyy_xxxyz[k];

                g_0_z_xxxxxyy_xxzz[k] = -g_0_z_xxxxyy_xxzz[k] * ab_x + g_0_z_xxxxyy_xxxzz[k];

                g_0_z_xxxxxyy_xyyy[k] = -g_0_z_xxxxyy_xyyy[k] * ab_x + g_0_z_xxxxyy_xxyyy[k];

                g_0_z_xxxxxyy_xyyz[k] = -g_0_z_xxxxyy_xyyz[k] * ab_x + g_0_z_xxxxyy_xxyyz[k];

                g_0_z_xxxxxyy_xyzz[k] = -g_0_z_xxxxyy_xyzz[k] * ab_x + g_0_z_xxxxyy_xxyzz[k];

                g_0_z_xxxxxyy_xzzz[k] = -g_0_z_xxxxyy_xzzz[k] * ab_x + g_0_z_xxxxyy_xxzzz[k];

                g_0_z_xxxxxyy_yyyy[k] = -g_0_z_xxxxyy_yyyy[k] * ab_x + g_0_z_xxxxyy_xyyyy[k];

                g_0_z_xxxxxyy_yyyz[k] = -g_0_z_xxxxyy_yyyz[k] * ab_x + g_0_z_xxxxyy_xyyyz[k];

                g_0_z_xxxxxyy_yyzz[k] = -g_0_z_xxxxyy_yyzz[k] * ab_x + g_0_z_xxxxyy_xyyzz[k];

                g_0_z_xxxxxyy_yzzz[k] = -g_0_z_xxxxyy_yzzz[k] * ab_x + g_0_z_xxxxyy_xyzzz[k];

                g_0_z_xxxxxyy_zzzz[k] = -g_0_z_xxxxyy_zzzz[k] * ab_x + g_0_z_xxxxyy_xzzzz[k];
            }

            /// Set up 1140-1155 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxyz_xxxx = cbuffer.data(kg_geom_01_off + 1140 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xxxy = cbuffer.data(kg_geom_01_off + 1141 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xxxz = cbuffer.data(kg_geom_01_off + 1142 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xxyy = cbuffer.data(kg_geom_01_off + 1143 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xxyz = cbuffer.data(kg_geom_01_off + 1144 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xxzz = cbuffer.data(kg_geom_01_off + 1145 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xyyy = cbuffer.data(kg_geom_01_off + 1146 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xyyz = cbuffer.data(kg_geom_01_off + 1147 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xyzz = cbuffer.data(kg_geom_01_off + 1148 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_xzzz = cbuffer.data(kg_geom_01_off + 1149 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_yyyy = cbuffer.data(kg_geom_01_off + 1150 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_yyyz = cbuffer.data(kg_geom_01_off + 1151 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_yyzz = cbuffer.data(kg_geom_01_off + 1152 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_yzzz = cbuffer.data(kg_geom_01_off + 1153 * ccomps * dcomps);

            auto g_0_z_xxxxxyz_zzzz = cbuffer.data(kg_geom_01_off + 1154 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxyz_xxxx, g_0_z_xxxxxyz_xxxy, g_0_z_xxxxxyz_xxxz, g_0_z_xxxxxyz_xxyy, g_0_z_xxxxxyz_xxyz, g_0_z_xxxxxyz_xxzz, g_0_z_xxxxxyz_xyyy, g_0_z_xxxxxyz_xyyz, g_0_z_xxxxxyz_xyzz, g_0_z_xxxxxyz_xzzz, g_0_z_xxxxxyz_yyyy, g_0_z_xxxxxyz_yyyz, g_0_z_xxxxxyz_yyzz, g_0_z_xxxxxyz_yzzz, g_0_z_xxxxxyz_zzzz, g_0_z_xxxxyz_xxxx, g_0_z_xxxxyz_xxxxx, g_0_z_xxxxyz_xxxxy, g_0_z_xxxxyz_xxxxz, g_0_z_xxxxyz_xxxy, g_0_z_xxxxyz_xxxyy, g_0_z_xxxxyz_xxxyz, g_0_z_xxxxyz_xxxz, g_0_z_xxxxyz_xxxzz, g_0_z_xxxxyz_xxyy, g_0_z_xxxxyz_xxyyy, g_0_z_xxxxyz_xxyyz, g_0_z_xxxxyz_xxyz, g_0_z_xxxxyz_xxyzz, g_0_z_xxxxyz_xxzz, g_0_z_xxxxyz_xxzzz, g_0_z_xxxxyz_xyyy, g_0_z_xxxxyz_xyyyy, g_0_z_xxxxyz_xyyyz, g_0_z_xxxxyz_xyyz, g_0_z_xxxxyz_xyyzz, g_0_z_xxxxyz_xyzz, g_0_z_xxxxyz_xyzzz, g_0_z_xxxxyz_xzzz, g_0_z_xxxxyz_xzzzz, g_0_z_xxxxyz_yyyy, g_0_z_xxxxyz_yyyz, g_0_z_xxxxyz_yyzz, g_0_z_xxxxyz_yzzz, g_0_z_xxxxyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxyz_xxxx[k] = -g_0_z_xxxxyz_xxxx[k] * ab_x + g_0_z_xxxxyz_xxxxx[k];

                g_0_z_xxxxxyz_xxxy[k] = -g_0_z_xxxxyz_xxxy[k] * ab_x + g_0_z_xxxxyz_xxxxy[k];

                g_0_z_xxxxxyz_xxxz[k] = -g_0_z_xxxxyz_xxxz[k] * ab_x + g_0_z_xxxxyz_xxxxz[k];

                g_0_z_xxxxxyz_xxyy[k] = -g_0_z_xxxxyz_xxyy[k] * ab_x + g_0_z_xxxxyz_xxxyy[k];

                g_0_z_xxxxxyz_xxyz[k] = -g_0_z_xxxxyz_xxyz[k] * ab_x + g_0_z_xxxxyz_xxxyz[k];

                g_0_z_xxxxxyz_xxzz[k] = -g_0_z_xxxxyz_xxzz[k] * ab_x + g_0_z_xxxxyz_xxxzz[k];

                g_0_z_xxxxxyz_xyyy[k] = -g_0_z_xxxxyz_xyyy[k] * ab_x + g_0_z_xxxxyz_xxyyy[k];

                g_0_z_xxxxxyz_xyyz[k] = -g_0_z_xxxxyz_xyyz[k] * ab_x + g_0_z_xxxxyz_xxyyz[k];

                g_0_z_xxxxxyz_xyzz[k] = -g_0_z_xxxxyz_xyzz[k] * ab_x + g_0_z_xxxxyz_xxyzz[k];

                g_0_z_xxxxxyz_xzzz[k] = -g_0_z_xxxxyz_xzzz[k] * ab_x + g_0_z_xxxxyz_xxzzz[k];

                g_0_z_xxxxxyz_yyyy[k] = -g_0_z_xxxxyz_yyyy[k] * ab_x + g_0_z_xxxxyz_xyyyy[k];

                g_0_z_xxxxxyz_yyyz[k] = -g_0_z_xxxxyz_yyyz[k] * ab_x + g_0_z_xxxxyz_xyyyz[k];

                g_0_z_xxxxxyz_yyzz[k] = -g_0_z_xxxxyz_yyzz[k] * ab_x + g_0_z_xxxxyz_xyyzz[k];

                g_0_z_xxxxxyz_yzzz[k] = -g_0_z_xxxxyz_yzzz[k] * ab_x + g_0_z_xxxxyz_xyzzz[k];

                g_0_z_xxxxxyz_zzzz[k] = -g_0_z_xxxxyz_zzzz[k] * ab_x + g_0_z_xxxxyz_xzzzz[k];
            }

            /// Set up 1155-1170 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxxzz_xxxx = cbuffer.data(kg_geom_01_off + 1155 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xxxy = cbuffer.data(kg_geom_01_off + 1156 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xxxz = cbuffer.data(kg_geom_01_off + 1157 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xxyy = cbuffer.data(kg_geom_01_off + 1158 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xxyz = cbuffer.data(kg_geom_01_off + 1159 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xxzz = cbuffer.data(kg_geom_01_off + 1160 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xyyy = cbuffer.data(kg_geom_01_off + 1161 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xyyz = cbuffer.data(kg_geom_01_off + 1162 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xyzz = cbuffer.data(kg_geom_01_off + 1163 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_xzzz = cbuffer.data(kg_geom_01_off + 1164 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_yyyy = cbuffer.data(kg_geom_01_off + 1165 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_yyyz = cbuffer.data(kg_geom_01_off + 1166 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_yyzz = cbuffer.data(kg_geom_01_off + 1167 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_yzzz = cbuffer.data(kg_geom_01_off + 1168 * ccomps * dcomps);

            auto g_0_z_xxxxxzz_zzzz = cbuffer.data(kg_geom_01_off + 1169 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxxzz_xxxx, g_0_z_xxxxxzz_xxxy, g_0_z_xxxxxzz_xxxz, g_0_z_xxxxxzz_xxyy, g_0_z_xxxxxzz_xxyz, g_0_z_xxxxxzz_xxzz, g_0_z_xxxxxzz_xyyy, g_0_z_xxxxxzz_xyyz, g_0_z_xxxxxzz_xyzz, g_0_z_xxxxxzz_xzzz, g_0_z_xxxxxzz_yyyy, g_0_z_xxxxxzz_yyyz, g_0_z_xxxxxzz_yyzz, g_0_z_xxxxxzz_yzzz, g_0_z_xxxxxzz_zzzz, g_0_z_xxxxzz_xxxx, g_0_z_xxxxzz_xxxxx, g_0_z_xxxxzz_xxxxy, g_0_z_xxxxzz_xxxxz, g_0_z_xxxxzz_xxxy, g_0_z_xxxxzz_xxxyy, g_0_z_xxxxzz_xxxyz, g_0_z_xxxxzz_xxxz, g_0_z_xxxxzz_xxxzz, g_0_z_xxxxzz_xxyy, g_0_z_xxxxzz_xxyyy, g_0_z_xxxxzz_xxyyz, g_0_z_xxxxzz_xxyz, g_0_z_xxxxzz_xxyzz, g_0_z_xxxxzz_xxzz, g_0_z_xxxxzz_xxzzz, g_0_z_xxxxzz_xyyy, g_0_z_xxxxzz_xyyyy, g_0_z_xxxxzz_xyyyz, g_0_z_xxxxzz_xyyz, g_0_z_xxxxzz_xyyzz, g_0_z_xxxxzz_xyzz, g_0_z_xxxxzz_xyzzz, g_0_z_xxxxzz_xzzz, g_0_z_xxxxzz_xzzzz, g_0_z_xxxxzz_yyyy, g_0_z_xxxxzz_yyyz, g_0_z_xxxxzz_yyzz, g_0_z_xxxxzz_yzzz, g_0_z_xxxxzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxxzz_xxxx[k] = -g_0_z_xxxxzz_xxxx[k] * ab_x + g_0_z_xxxxzz_xxxxx[k];

                g_0_z_xxxxxzz_xxxy[k] = -g_0_z_xxxxzz_xxxy[k] * ab_x + g_0_z_xxxxzz_xxxxy[k];

                g_0_z_xxxxxzz_xxxz[k] = -g_0_z_xxxxzz_xxxz[k] * ab_x + g_0_z_xxxxzz_xxxxz[k];

                g_0_z_xxxxxzz_xxyy[k] = -g_0_z_xxxxzz_xxyy[k] * ab_x + g_0_z_xxxxzz_xxxyy[k];

                g_0_z_xxxxxzz_xxyz[k] = -g_0_z_xxxxzz_xxyz[k] * ab_x + g_0_z_xxxxzz_xxxyz[k];

                g_0_z_xxxxxzz_xxzz[k] = -g_0_z_xxxxzz_xxzz[k] * ab_x + g_0_z_xxxxzz_xxxzz[k];

                g_0_z_xxxxxzz_xyyy[k] = -g_0_z_xxxxzz_xyyy[k] * ab_x + g_0_z_xxxxzz_xxyyy[k];

                g_0_z_xxxxxzz_xyyz[k] = -g_0_z_xxxxzz_xyyz[k] * ab_x + g_0_z_xxxxzz_xxyyz[k];

                g_0_z_xxxxxzz_xyzz[k] = -g_0_z_xxxxzz_xyzz[k] * ab_x + g_0_z_xxxxzz_xxyzz[k];

                g_0_z_xxxxxzz_xzzz[k] = -g_0_z_xxxxzz_xzzz[k] * ab_x + g_0_z_xxxxzz_xxzzz[k];

                g_0_z_xxxxxzz_yyyy[k] = -g_0_z_xxxxzz_yyyy[k] * ab_x + g_0_z_xxxxzz_xyyyy[k];

                g_0_z_xxxxxzz_yyyz[k] = -g_0_z_xxxxzz_yyyz[k] * ab_x + g_0_z_xxxxzz_xyyyz[k];

                g_0_z_xxxxxzz_yyzz[k] = -g_0_z_xxxxzz_yyzz[k] * ab_x + g_0_z_xxxxzz_xyyzz[k];

                g_0_z_xxxxxzz_yzzz[k] = -g_0_z_xxxxzz_yzzz[k] * ab_x + g_0_z_xxxxzz_xyzzz[k];

                g_0_z_xxxxxzz_zzzz[k] = -g_0_z_xxxxzz_zzzz[k] * ab_x + g_0_z_xxxxzz_xzzzz[k];
            }

            /// Set up 1170-1185 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyyy_xxxx = cbuffer.data(kg_geom_01_off + 1170 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xxxy = cbuffer.data(kg_geom_01_off + 1171 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xxxz = cbuffer.data(kg_geom_01_off + 1172 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xxyy = cbuffer.data(kg_geom_01_off + 1173 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xxyz = cbuffer.data(kg_geom_01_off + 1174 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xxzz = cbuffer.data(kg_geom_01_off + 1175 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xyyy = cbuffer.data(kg_geom_01_off + 1176 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xyyz = cbuffer.data(kg_geom_01_off + 1177 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xyzz = cbuffer.data(kg_geom_01_off + 1178 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_xzzz = cbuffer.data(kg_geom_01_off + 1179 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_yyyy = cbuffer.data(kg_geom_01_off + 1180 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_yyyz = cbuffer.data(kg_geom_01_off + 1181 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_yyzz = cbuffer.data(kg_geom_01_off + 1182 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_yzzz = cbuffer.data(kg_geom_01_off + 1183 * ccomps * dcomps);

            auto g_0_z_xxxxyyy_zzzz = cbuffer.data(kg_geom_01_off + 1184 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyyy_xxxx, g_0_z_xxxxyyy_xxxy, g_0_z_xxxxyyy_xxxz, g_0_z_xxxxyyy_xxyy, g_0_z_xxxxyyy_xxyz, g_0_z_xxxxyyy_xxzz, g_0_z_xxxxyyy_xyyy, g_0_z_xxxxyyy_xyyz, g_0_z_xxxxyyy_xyzz, g_0_z_xxxxyyy_xzzz, g_0_z_xxxxyyy_yyyy, g_0_z_xxxxyyy_yyyz, g_0_z_xxxxyyy_yyzz, g_0_z_xxxxyyy_yzzz, g_0_z_xxxxyyy_zzzz, g_0_z_xxxyyy_xxxx, g_0_z_xxxyyy_xxxxx, g_0_z_xxxyyy_xxxxy, g_0_z_xxxyyy_xxxxz, g_0_z_xxxyyy_xxxy, g_0_z_xxxyyy_xxxyy, g_0_z_xxxyyy_xxxyz, g_0_z_xxxyyy_xxxz, g_0_z_xxxyyy_xxxzz, g_0_z_xxxyyy_xxyy, g_0_z_xxxyyy_xxyyy, g_0_z_xxxyyy_xxyyz, g_0_z_xxxyyy_xxyz, g_0_z_xxxyyy_xxyzz, g_0_z_xxxyyy_xxzz, g_0_z_xxxyyy_xxzzz, g_0_z_xxxyyy_xyyy, g_0_z_xxxyyy_xyyyy, g_0_z_xxxyyy_xyyyz, g_0_z_xxxyyy_xyyz, g_0_z_xxxyyy_xyyzz, g_0_z_xxxyyy_xyzz, g_0_z_xxxyyy_xyzzz, g_0_z_xxxyyy_xzzz, g_0_z_xxxyyy_xzzzz, g_0_z_xxxyyy_yyyy, g_0_z_xxxyyy_yyyz, g_0_z_xxxyyy_yyzz, g_0_z_xxxyyy_yzzz, g_0_z_xxxyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyyy_xxxx[k] = -g_0_z_xxxyyy_xxxx[k] * ab_x + g_0_z_xxxyyy_xxxxx[k];

                g_0_z_xxxxyyy_xxxy[k] = -g_0_z_xxxyyy_xxxy[k] * ab_x + g_0_z_xxxyyy_xxxxy[k];

                g_0_z_xxxxyyy_xxxz[k] = -g_0_z_xxxyyy_xxxz[k] * ab_x + g_0_z_xxxyyy_xxxxz[k];

                g_0_z_xxxxyyy_xxyy[k] = -g_0_z_xxxyyy_xxyy[k] * ab_x + g_0_z_xxxyyy_xxxyy[k];

                g_0_z_xxxxyyy_xxyz[k] = -g_0_z_xxxyyy_xxyz[k] * ab_x + g_0_z_xxxyyy_xxxyz[k];

                g_0_z_xxxxyyy_xxzz[k] = -g_0_z_xxxyyy_xxzz[k] * ab_x + g_0_z_xxxyyy_xxxzz[k];

                g_0_z_xxxxyyy_xyyy[k] = -g_0_z_xxxyyy_xyyy[k] * ab_x + g_0_z_xxxyyy_xxyyy[k];

                g_0_z_xxxxyyy_xyyz[k] = -g_0_z_xxxyyy_xyyz[k] * ab_x + g_0_z_xxxyyy_xxyyz[k];

                g_0_z_xxxxyyy_xyzz[k] = -g_0_z_xxxyyy_xyzz[k] * ab_x + g_0_z_xxxyyy_xxyzz[k];

                g_0_z_xxxxyyy_xzzz[k] = -g_0_z_xxxyyy_xzzz[k] * ab_x + g_0_z_xxxyyy_xxzzz[k];

                g_0_z_xxxxyyy_yyyy[k] = -g_0_z_xxxyyy_yyyy[k] * ab_x + g_0_z_xxxyyy_xyyyy[k];

                g_0_z_xxxxyyy_yyyz[k] = -g_0_z_xxxyyy_yyyz[k] * ab_x + g_0_z_xxxyyy_xyyyz[k];

                g_0_z_xxxxyyy_yyzz[k] = -g_0_z_xxxyyy_yyzz[k] * ab_x + g_0_z_xxxyyy_xyyzz[k];

                g_0_z_xxxxyyy_yzzz[k] = -g_0_z_xxxyyy_yzzz[k] * ab_x + g_0_z_xxxyyy_xyzzz[k];

                g_0_z_xxxxyyy_zzzz[k] = -g_0_z_xxxyyy_zzzz[k] * ab_x + g_0_z_xxxyyy_xzzzz[k];
            }

            /// Set up 1185-1200 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyyz_xxxx = cbuffer.data(kg_geom_01_off + 1185 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xxxy = cbuffer.data(kg_geom_01_off + 1186 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xxxz = cbuffer.data(kg_geom_01_off + 1187 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xxyy = cbuffer.data(kg_geom_01_off + 1188 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xxyz = cbuffer.data(kg_geom_01_off + 1189 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xxzz = cbuffer.data(kg_geom_01_off + 1190 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xyyy = cbuffer.data(kg_geom_01_off + 1191 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xyyz = cbuffer.data(kg_geom_01_off + 1192 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xyzz = cbuffer.data(kg_geom_01_off + 1193 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_xzzz = cbuffer.data(kg_geom_01_off + 1194 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_yyyy = cbuffer.data(kg_geom_01_off + 1195 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_yyyz = cbuffer.data(kg_geom_01_off + 1196 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_yyzz = cbuffer.data(kg_geom_01_off + 1197 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_yzzz = cbuffer.data(kg_geom_01_off + 1198 * ccomps * dcomps);

            auto g_0_z_xxxxyyz_zzzz = cbuffer.data(kg_geom_01_off + 1199 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyyz_xxxx, g_0_z_xxxxyyz_xxxy, g_0_z_xxxxyyz_xxxz, g_0_z_xxxxyyz_xxyy, g_0_z_xxxxyyz_xxyz, g_0_z_xxxxyyz_xxzz, g_0_z_xxxxyyz_xyyy, g_0_z_xxxxyyz_xyyz, g_0_z_xxxxyyz_xyzz, g_0_z_xxxxyyz_xzzz, g_0_z_xxxxyyz_yyyy, g_0_z_xxxxyyz_yyyz, g_0_z_xxxxyyz_yyzz, g_0_z_xxxxyyz_yzzz, g_0_z_xxxxyyz_zzzz, g_0_z_xxxyyz_xxxx, g_0_z_xxxyyz_xxxxx, g_0_z_xxxyyz_xxxxy, g_0_z_xxxyyz_xxxxz, g_0_z_xxxyyz_xxxy, g_0_z_xxxyyz_xxxyy, g_0_z_xxxyyz_xxxyz, g_0_z_xxxyyz_xxxz, g_0_z_xxxyyz_xxxzz, g_0_z_xxxyyz_xxyy, g_0_z_xxxyyz_xxyyy, g_0_z_xxxyyz_xxyyz, g_0_z_xxxyyz_xxyz, g_0_z_xxxyyz_xxyzz, g_0_z_xxxyyz_xxzz, g_0_z_xxxyyz_xxzzz, g_0_z_xxxyyz_xyyy, g_0_z_xxxyyz_xyyyy, g_0_z_xxxyyz_xyyyz, g_0_z_xxxyyz_xyyz, g_0_z_xxxyyz_xyyzz, g_0_z_xxxyyz_xyzz, g_0_z_xxxyyz_xyzzz, g_0_z_xxxyyz_xzzz, g_0_z_xxxyyz_xzzzz, g_0_z_xxxyyz_yyyy, g_0_z_xxxyyz_yyyz, g_0_z_xxxyyz_yyzz, g_0_z_xxxyyz_yzzz, g_0_z_xxxyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyyz_xxxx[k] = -g_0_z_xxxyyz_xxxx[k] * ab_x + g_0_z_xxxyyz_xxxxx[k];

                g_0_z_xxxxyyz_xxxy[k] = -g_0_z_xxxyyz_xxxy[k] * ab_x + g_0_z_xxxyyz_xxxxy[k];

                g_0_z_xxxxyyz_xxxz[k] = -g_0_z_xxxyyz_xxxz[k] * ab_x + g_0_z_xxxyyz_xxxxz[k];

                g_0_z_xxxxyyz_xxyy[k] = -g_0_z_xxxyyz_xxyy[k] * ab_x + g_0_z_xxxyyz_xxxyy[k];

                g_0_z_xxxxyyz_xxyz[k] = -g_0_z_xxxyyz_xxyz[k] * ab_x + g_0_z_xxxyyz_xxxyz[k];

                g_0_z_xxxxyyz_xxzz[k] = -g_0_z_xxxyyz_xxzz[k] * ab_x + g_0_z_xxxyyz_xxxzz[k];

                g_0_z_xxxxyyz_xyyy[k] = -g_0_z_xxxyyz_xyyy[k] * ab_x + g_0_z_xxxyyz_xxyyy[k];

                g_0_z_xxxxyyz_xyyz[k] = -g_0_z_xxxyyz_xyyz[k] * ab_x + g_0_z_xxxyyz_xxyyz[k];

                g_0_z_xxxxyyz_xyzz[k] = -g_0_z_xxxyyz_xyzz[k] * ab_x + g_0_z_xxxyyz_xxyzz[k];

                g_0_z_xxxxyyz_xzzz[k] = -g_0_z_xxxyyz_xzzz[k] * ab_x + g_0_z_xxxyyz_xxzzz[k];

                g_0_z_xxxxyyz_yyyy[k] = -g_0_z_xxxyyz_yyyy[k] * ab_x + g_0_z_xxxyyz_xyyyy[k];

                g_0_z_xxxxyyz_yyyz[k] = -g_0_z_xxxyyz_yyyz[k] * ab_x + g_0_z_xxxyyz_xyyyz[k];

                g_0_z_xxxxyyz_yyzz[k] = -g_0_z_xxxyyz_yyzz[k] * ab_x + g_0_z_xxxyyz_xyyzz[k];

                g_0_z_xxxxyyz_yzzz[k] = -g_0_z_xxxyyz_yzzz[k] * ab_x + g_0_z_xxxyyz_xyzzz[k];

                g_0_z_xxxxyyz_zzzz[k] = -g_0_z_xxxyyz_zzzz[k] * ab_x + g_0_z_xxxyyz_xzzzz[k];
            }

            /// Set up 1200-1215 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxyzz_xxxx = cbuffer.data(kg_geom_01_off + 1200 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xxxy = cbuffer.data(kg_geom_01_off + 1201 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xxxz = cbuffer.data(kg_geom_01_off + 1202 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xxyy = cbuffer.data(kg_geom_01_off + 1203 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xxyz = cbuffer.data(kg_geom_01_off + 1204 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xxzz = cbuffer.data(kg_geom_01_off + 1205 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xyyy = cbuffer.data(kg_geom_01_off + 1206 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xyyz = cbuffer.data(kg_geom_01_off + 1207 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xyzz = cbuffer.data(kg_geom_01_off + 1208 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_xzzz = cbuffer.data(kg_geom_01_off + 1209 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_yyyy = cbuffer.data(kg_geom_01_off + 1210 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_yyyz = cbuffer.data(kg_geom_01_off + 1211 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_yyzz = cbuffer.data(kg_geom_01_off + 1212 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_yzzz = cbuffer.data(kg_geom_01_off + 1213 * ccomps * dcomps);

            auto g_0_z_xxxxyzz_zzzz = cbuffer.data(kg_geom_01_off + 1214 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxyzz_xxxx, g_0_z_xxxxyzz_xxxy, g_0_z_xxxxyzz_xxxz, g_0_z_xxxxyzz_xxyy, g_0_z_xxxxyzz_xxyz, g_0_z_xxxxyzz_xxzz, g_0_z_xxxxyzz_xyyy, g_0_z_xxxxyzz_xyyz, g_0_z_xxxxyzz_xyzz, g_0_z_xxxxyzz_xzzz, g_0_z_xxxxyzz_yyyy, g_0_z_xxxxyzz_yyyz, g_0_z_xxxxyzz_yyzz, g_0_z_xxxxyzz_yzzz, g_0_z_xxxxyzz_zzzz, g_0_z_xxxyzz_xxxx, g_0_z_xxxyzz_xxxxx, g_0_z_xxxyzz_xxxxy, g_0_z_xxxyzz_xxxxz, g_0_z_xxxyzz_xxxy, g_0_z_xxxyzz_xxxyy, g_0_z_xxxyzz_xxxyz, g_0_z_xxxyzz_xxxz, g_0_z_xxxyzz_xxxzz, g_0_z_xxxyzz_xxyy, g_0_z_xxxyzz_xxyyy, g_0_z_xxxyzz_xxyyz, g_0_z_xxxyzz_xxyz, g_0_z_xxxyzz_xxyzz, g_0_z_xxxyzz_xxzz, g_0_z_xxxyzz_xxzzz, g_0_z_xxxyzz_xyyy, g_0_z_xxxyzz_xyyyy, g_0_z_xxxyzz_xyyyz, g_0_z_xxxyzz_xyyz, g_0_z_xxxyzz_xyyzz, g_0_z_xxxyzz_xyzz, g_0_z_xxxyzz_xyzzz, g_0_z_xxxyzz_xzzz, g_0_z_xxxyzz_xzzzz, g_0_z_xxxyzz_yyyy, g_0_z_xxxyzz_yyyz, g_0_z_xxxyzz_yyzz, g_0_z_xxxyzz_yzzz, g_0_z_xxxyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxyzz_xxxx[k] = -g_0_z_xxxyzz_xxxx[k] * ab_x + g_0_z_xxxyzz_xxxxx[k];

                g_0_z_xxxxyzz_xxxy[k] = -g_0_z_xxxyzz_xxxy[k] * ab_x + g_0_z_xxxyzz_xxxxy[k];

                g_0_z_xxxxyzz_xxxz[k] = -g_0_z_xxxyzz_xxxz[k] * ab_x + g_0_z_xxxyzz_xxxxz[k];

                g_0_z_xxxxyzz_xxyy[k] = -g_0_z_xxxyzz_xxyy[k] * ab_x + g_0_z_xxxyzz_xxxyy[k];

                g_0_z_xxxxyzz_xxyz[k] = -g_0_z_xxxyzz_xxyz[k] * ab_x + g_0_z_xxxyzz_xxxyz[k];

                g_0_z_xxxxyzz_xxzz[k] = -g_0_z_xxxyzz_xxzz[k] * ab_x + g_0_z_xxxyzz_xxxzz[k];

                g_0_z_xxxxyzz_xyyy[k] = -g_0_z_xxxyzz_xyyy[k] * ab_x + g_0_z_xxxyzz_xxyyy[k];

                g_0_z_xxxxyzz_xyyz[k] = -g_0_z_xxxyzz_xyyz[k] * ab_x + g_0_z_xxxyzz_xxyyz[k];

                g_0_z_xxxxyzz_xyzz[k] = -g_0_z_xxxyzz_xyzz[k] * ab_x + g_0_z_xxxyzz_xxyzz[k];

                g_0_z_xxxxyzz_xzzz[k] = -g_0_z_xxxyzz_xzzz[k] * ab_x + g_0_z_xxxyzz_xxzzz[k];

                g_0_z_xxxxyzz_yyyy[k] = -g_0_z_xxxyzz_yyyy[k] * ab_x + g_0_z_xxxyzz_xyyyy[k];

                g_0_z_xxxxyzz_yyyz[k] = -g_0_z_xxxyzz_yyyz[k] * ab_x + g_0_z_xxxyzz_xyyyz[k];

                g_0_z_xxxxyzz_yyzz[k] = -g_0_z_xxxyzz_yyzz[k] * ab_x + g_0_z_xxxyzz_xyyzz[k];

                g_0_z_xxxxyzz_yzzz[k] = -g_0_z_xxxyzz_yzzz[k] * ab_x + g_0_z_xxxyzz_xyzzz[k];

                g_0_z_xxxxyzz_zzzz[k] = -g_0_z_xxxyzz_zzzz[k] * ab_x + g_0_z_xxxyzz_xzzzz[k];
            }

            /// Set up 1215-1230 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxxzzz_xxxx = cbuffer.data(kg_geom_01_off + 1215 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xxxy = cbuffer.data(kg_geom_01_off + 1216 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xxxz = cbuffer.data(kg_geom_01_off + 1217 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xxyy = cbuffer.data(kg_geom_01_off + 1218 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xxyz = cbuffer.data(kg_geom_01_off + 1219 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xxzz = cbuffer.data(kg_geom_01_off + 1220 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xyyy = cbuffer.data(kg_geom_01_off + 1221 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xyyz = cbuffer.data(kg_geom_01_off + 1222 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xyzz = cbuffer.data(kg_geom_01_off + 1223 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_xzzz = cbuffer.data(kg_geom_01_off + 1224 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_yyyy = cbuffer.data(kg_geom_01_off + 1225 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_yyyz = cbuffer.data(kg_geom_01_off + 1226 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_yyzz = cbuffer.data(kg_geom_01_off + 1227 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_yzzz = cbuffer.data(kg_geom_01_off + 1228 * ccomps * dcomps);

            auto g_0_z_xxxxzzz_zzzz = cbuffer.data(kg_geom_01_off + 1229 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxxzzz_xxxx, g_0_z_xxxxzzz_xxxy, g_0_z_xxxxzzz_xxxz, g_0_z_xxxxzzz_xxyy, g_0_z_xxxxzzz_xxyz, g_0_z_xxxxzzz_xxzz, g_0_z_xxxxzzz_xyyy, g_0_z_xxxxzzz_xyyz, g_0_z_xxxxzzz_xyzz, g_0_z_xxxxzzz_xzzz, g_0_z_xxxxzzz_yyyy, g_0_z_xxxxzzz_yyyz, g_0_z_xxxxzzz_yyzz, g_0_z_xxxxzzz_yzzz, g_0_z_xxxxzzz_zzzz, g_0_z_xxxzzz_xxxx, g_0_z_xxxzzz_xxxxx, g_0_z_xxxzzz_xxxxy, g_0_z_xxxzzz_xxxxz, g_0_z_xxxzzz_xxxy, g_0_z_xxxzzz_xxxyy, g_0_z_xxxzzz_xxxyz, g_0_z_xxxzzz_xxxz, g_0_z_xxxzzz_xxxzz, g_0_z_xxxzzz_xxyy, g_0_z_xxxzzz_xxyyy, g_0_z_xxxzzz_xxyyz, g_0_z_xxxzzz_xxyz, g_0_z_xxxzzz_xxyzz, g_0_z_xxxzzz_xxzz, g_0_z_xxxzzz_xxzzz, g_0_z_xxxzzz_xyyy, g_0_z_xxxzzz_xyyyy, g_0_z_xxxzzz_xyyyz, g_0_z_xxxzzz_xyyz, g_0_z_xxxzzz_xyyzz, g_0_z_xxxzzz_xyzz, g_0_z_xxxzzz_xyzzz, g_0_z_xxxzzz_xzzz, g_0_z_xxxzzz_xzzzz, g_0_z_xxxzzz_yyyy, g_0_z_xxxzzz_yyyz, g_0_z_xxxzzz_yyzz, g_0_z_xxxzzz_yzzz, g_0_z_xxxzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxxzzz_xxxx[k] = -g_0_z_xxxzzz_xxxx[k] * ab_x + g_0_z_xxxzzz_xxxxx[k];

                g_0_z_xxxxzzz_xxxy[k] = -g_0_z_xxxzzz_xxxy[k] * ab_x + g_0_z_xxxzzz_xxxxy[k];

                g_0_z_xxxxzzz_xxxz[k] = -g_0_z_xxxzzz_xxxz[k] * ab_x + g_0_z_xxxzzz_xxxxz[k];

                g_0_z_xxxxzzz_xxyy[k] = -g_0_z_xxxzzz_xxyy[k] * ab_x + g_0_z_xxxzzz_xxxyy[k];

                g_0_z_xxxxzzz_xxyz[k] = -g_0_z_xxxzzz_xxyz[k] * ab_x + g_0_z_xxxzzz_xxxyz[k];

                g_0_z_xxxxzzz_xxzz[k] = -g_0_z_xxxzzz_xxzz[k] * ab_x + g_0_z_xxxzzz_xxxzz[k];

                g_0_z_xxxxzzz_xyyy[k] = -g_0_z_xxxzzz_xyyy[k] * ab_x + g_0_z_xxxzzz_xxyyy[k];

                g_0_z_xxxxzzz_xyyz[k] = -g_0_z_xxxzzz_xyyz[k] * ab_x + g_0_z_xxxzzz_xxyyz[k];

                g_0_z_xxxxzzz_xyzz[k] = -g_0_z_xxxzzz_xyzz[k] * ab_x + g_0_z_xxxzzz_xxyzz[k];

                g_0_z_xxxxzzz_xzzz[k] = -g_0_z_xxxzzz_xzzz[k] * ab_x + g_0_z_xxxzzz_xxzzz[k];

                g_0_z_xxxxzzz_yyyy[k] = -g_0_z_xxxzzz_yyyy[k] * ab_x + g_0_z_xxxzzz_xyyyy[k];

                g_0_z_xxxxzzz_yyyz[k] = -g_0_z_xxxzzz_yyyz[k] * ab_x + g_0_z_xxxzzz_xyyyz[k];

                g_0_z_xxxxzzz_yyzz[k] = -g_0_z_xxxzzz_yyzz[k] * ab_x + g_0_z_xxxzzz_xyyzz[k];

                g_0_z_xxxxzzz_yzzz[k] = -g_0_z_xxxzzz_yzzz[k] * ab_x + g_0_z_xxxzzz_xyzzz[k];

                g_0_z_xxxxzzz_zzzz[k] = -g_0_z_xxxzzz_zzzz[k] * ab_x + g_0_z_xxxzzz_xzzzz[k];
            }

            /// Set up 1230-1245 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyyy_xxxx = cbuffer.data(kg_geom_01_off + 1230 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xxxy = cbuffer.data(kg_geom_01_off + 1231 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xxxz = cbuffer.data(kg_geom_01_off + 1232 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xxyy = cbuffer.data(kg_geom_01_off + 1233 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xxyz = cbuffer.data(kg_geom_01_off + 1234 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xxzz = cbuffer.data(kg_geom_01_off + 1235 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xyyy = cbuffer.data(kg_geom_01_off + 1236 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xyyz = cbuffer.data(kg_geom_01_off + 1237 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xyzz = cbuffer.data(kg_geom_01_off + 1238 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_xzzz = cbuffer.data(kg_geom_01_off + 1239 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_yyyy = cbuffer.data(kg_geom_01_off + 1240 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_yyyz = cbuffer.data(kg_geom_01_off + 1241 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_yyzz = cbuffer.data(kg_geom_01_off + 1242 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_yzzz = cbuffer.data(kg_geom_01_off + 1243 * ccomps * dcomps);

            auto g_0_z_xxxyyyy_zzzz = cbuffer.data(kg_geom_01_off + 1244 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyyy_xxxx, g_0_z_xxxyyyy_xxxy, g_0_z_xxxyyyy_xxxz, g_0_z_xxxyyyy_xxyy, g_0_z_xxxyyyy_xxyz, g_0_z_xxxyyyy_xxzz, g_0_z_xxxyyyy_xyyy, g_0_z_xxxyyyy_xyyz, g_0_z_xxxyyyy_xyzz, g_0_z_xxxyyyy_xzzz, g_0_z_xxxyyyy_yyyy, g_0_z_xxxyyyy_yyyz, g_0_z_xxxyyyy_yyzz, g_0_z_xxxyyyy_yzzz, g_0_z_xxxyyyy_zzzz, g_0_z_xxyyyy_xxxx, g_0_z_xxyyyy_xxxxx, g_0_z_xxyyyy_xxxxy, g_0_z_xxyyyy_xxxxz, g_0_z_xxyyyy_xxxy, g_0_z_xxyyyy_xxxyy, g_0_z_xxyyyy_xxxyz, g_0_z_xxyyyy_xxxz, g_0_z_xxyyyy_xxxzz, g_0_z_xxyyyy_xxyy, g_0_z_xxyyyy_xxyyy, g_0_z_xxyyyy_xxyyz, g_0_z_xxyyyy_xxyz, g_0_z_xxyyyy_xxyzz, g_0_z_xxyyyy_xxzz, g_0_z_xxyyyy_xxzzz, g_0_z_xxyyyy_xyyy, g_0_z_xxyyyy_xyyyy, g_0_z_xxyyyy_xyyyz, g_0_z_xxyyyy_xyyz, g_0_z_xxyyyy_xyyzz, g_0_z_xxyyyy_xyzz, g_0_z_xxyyyy_xyzzz, g_0_z_xxyyyy_xzzz, g_0_z_xxyyyy_xzzzz, g_0_z_xxyyyy_yyyy, g_0_z_xxyyyy_yyyz, g_0_z_xxyyyy_yyzz, g_0_z_xxyyyy_yzzz, g_0_z_xxyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyyy_xxxx[k] = -g_0_z_xxyyyy_xxxx[k] * ab_x + g_0_z_xxyyyy_xxxxx[k];

                g_0_z_xxxyyyy_xxxy[k] = -g_0_z_xxyyyy_xxxy[k] * ab_x + g_0_z_xxyyyy_xxxxy[k];

                g_0_z_xxxyyyy_xxxz[k] = -g_0_z_xxyyyy_xxxz[k] * ab_x + g_0_z_xxyyyy_xxxxz[k];

                g_0_z_xxxyyyy_xxyy[k] = -g_0_z_xxyyyy_xxyy[k] * ab_x + g_0_z_xxyyyy_xxxyy[k];

                g_0_z_xxxyyyy_xxyz[k] = -g_0_z_xxyyyy_xxyz[k] * ab_x + g_0_z_xxyyyy_xxxyz[k];

                g_0_z_xxxyyyy_xxzz[k] = -g_0_z_xxyyyy_xxzz[k] * ab_x + g_0_z_xxyyyy_xxxzz[k];

                g_0_z_xxxyyyy_xyyy[k] = -g_0_z_xxyyyy_xyyy[k] * ab_x + g_0_z_xxyyyy_xxyyy[k];

                g_0_z_xxxyyyy_xyyz[k] = -g_0_z_xxyyyy_xyyz[k] * ab_x + g_0_z_xxyyyy_xxyyz[k];

                g_0_z_xxxyyyy_xyzz[k] = -g_0_z_xxyyyy_xyzz[k] * ab_x + g_0_z_xxyyyy_xxyzz[k];

                g_0_z_xxxyyyy_xzzz[k] = -g_0_z_xxyyyy_xzzz[k] * ab_x + g_0_z_xxyyyy_xxzzz[k];

                g_0_z_xxxyyyy_yyyy[k] = -g_0_z_xxyyyy_yyyy[k] * ab_x + g_0_z_xxyyyy_xyyyy[k];

                g_0_z_xxxyyyy_yyyz[k] = -g_0_z_xxyyyy_yyyz[k] * ab_x + g_0_z_xxyyyy_xyyyz[k];

                g_0_z_xxxyyyy_yyzz[k] = -g_0_z_xxyyyy_yyzz[k] * ab_x + g_0_z_xxyyyy_xyyzz[k];

                g_0_z_xxxyyyy_yzzz[k] = -g_0_z_xxyyyy_yzzz[k] * ab_x + g_0_z_xxyyyy_xyzzz[k];

                g_0_z_xxxyyyy_zzzz[k] = -g_0_z_xxyyyy_zzzz[k] * ab_x + g_0_z_xxyyyy_xzzzz[k];
            }

            /// Set up 1245-1260 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyyz_xxxx = cbuffer.data(kg_geom_01_off + 1245 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xxxy = cbuffer.data(kg_geom_01_off + 1246 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xxxz = cbuffer.data(kg_geom_01_off + 1247 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xxyy = cbuffer.data(kg_geom_01_off + 1248 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xxyz = cbuffer.data(kg_geom_01_off + 1249 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xxzz = cbuffer.data(kg_geom_01_off + 1250 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xyyy = cbuffer.data(kg_geom_01_off + 1251 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xyyz = cbuffer.data(kg_geom_01_off + 1252 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xyzz = cbuffer.data(kg_geom_01_off + 1253 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_xzzz = cbuffer.data(kg_geom_01_off + 1254 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_yyyy = cbuffer.data(kg_geom_01_off + 1255 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_yyyz = cbuffer.data(kg_geom_01_off + 1256 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_yyzz = cbuffer.data(kg_geom_01_off + 1257 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_yzzz = cbuffer.data(kg_geom_01_off + 1258 * ccomps * dcomps);

            auto g_0_z_xxxyyyz_zzzz = cbuffer.data(kg_geom_01_off + 1259 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyyz_xxxx, g_0_z_xxxyyyz_xxxy, g_0_z_xxxyyyz_xxxz, g_0_z_xxxyyyz_xxyy, g_0_z_xxxyyyz_xxyz, g_0_z_xxxyyyz_xxzz, g_0_z_xxxyyyz_xyyy, g_0_z_xxxyyyz_xyyz, g_0_z_xxxyyyz_xyzz, g_0_z_xxxyyyz_xzzz, g_0_z_xxxyyyz_yyyy, g_0_z_xxxyyyz_yyyz, g_0_z_xxxyyyz_yyzz, g_0_z_xxxyyyz_yzzz, g_0_z_xxxyyyz_zzzz, g_0_z_xxyyyz_xxxx, g_0_z_xxyyyz_xxxxx, g_0_z_xxyyyz_xxxxy, g_0_z_xxyyyz_xxxxz, g_0_z_xxyyyz_xxxy, g_0_z_xxyyyz_xxxyy, g_0_z_xxyyyz_xxxyz, g_0_z_xxyyyz_xxxz, g_0_z_xxyyyz_xxxzz, g_0_z_xxyyyz_xxyy, g_0_z_xxyyyz_xxyyy, g_0_z_xxyyyz_xxyyz, g_0_z_xxyyyz_xxyz, g_0_z_xxyyyz_xxyzz, g_0_z_xxyyyz_xxzz, g_0_z_xxyyyz_xxzzz, g_0_z_xxyyyz_xyyy, g_0_z_xxyyyz_xyyyy, g_0_z_xxyyyz_xyyyz, g_0_z_xxyyyz_xyyz, g_0_z_xxyyyz_xyyzz, g_0_z_xxyyyz_xyzz, g_0_z_xxyyyz_xyzzz, g_0_z_xxyyyz_xzzz, g_0_z_xxyyyz_xzzzz, g_0_z_xxyyyz_yyyy, g_0_z_xxyyyz_yyyz, g_0_z_xxyyyz_yyzz, g_0_z_xxyyyz_yzzz, g_0_z_xxyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyyz_xxxx[k] = -g_0_z_xxyyyz_xxxx[k] * ab_x + g_0_z_xxyyyz_xxxxx[k];

                g_0_z_xxxyyyz_xxxy[k] = -g_0_z_xxyyyz_xxxy[k] * ab_x + g_0_z_xxyyyz_xxxxy[k];

                g_0_z_xxxyyyz_xxxz[k] = -g_0_z_xxyyyz_xxxz[k] * ab_x + g_0_z_xxyyyz_xxxxz[k];

                g_0_z_xxxyyyz_xxyy[k] = -g_0_z_xxyyyz_xxyy[k] * ab_x + g_0_z_xxyyyz_xxxyy[k];

                g_0_z_xxxyyyz_xxyz[k] = -g_0_z_xxyyyz_xxyz[k] * ab_x + g_0_z_xxyyyz_xxxyz[k];

                g_0_z_xxxyyyz_xxzz[k] = -g_0_z_xxyyyz_xxzz[k] * ab_x + g_0_z_xxyyyz_xxxzz[k];

                g_0_z_xxxyyyz_xyyy[k] = -g_0_z_xxyyyz_xyyy[k] * ab_x + g_0_z_xxyyyz_xxyyy[k];

                g_0_z_xxxyyyz_xyyz[k] = -g_0_z_xxyyyz_xyyz[k] * ab_x + g_0_z_xxyyyz_xxyyz[k];

                g_0_z_xxxyyyz_xyzz[k] = -g_0_z_xxyyyz_xyzz[k] * ab_x + g_0_z_xxyyyz_xxyzz[k];

                g_0_z_xxxyyyz_xzzz[k] = -g_0_z_xxyyyz_xzzz[k] * ab_x + g_0_z_xxyyyz_xxzzz[k];

                g_0_z_xxxyyyz_yyyy[k] = -g_0_z_xxyyyz_yyyy[k] * ab_x + g_0_z_xxyyyz_xyyyy[k];

                g_0_z_xxxyyyz_yyyz[k] = -g_0_z_xxyyyz_yyyz[k] * ab_x + g_0_z_xxyyyz_xyyyz[k];

                g_0_z_xxxyyyz_yyzz[k] = -g_0_z_xxyyyz_yyzz[k] * ab_x + g_0_z_xxyyyz_xyyzz[k];

                g_0_z_xxxyyyz_yzzz[k] = -g_0_z_xxyyyz_yzzz[k] * ab_x + g_0_z_xxyyyz_xyzzz[k];

                g_0_z_xxxyyyz_zzzz[k] = -g_0_z_xxyyyz_zzzz[k] * ab_x + g_0_z_xxyyyz_xzzzz[k];
            }

            /// Set up 1260-1275 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyyzz_xxxx = cbuffer.data(kg_geom_01_off + 1260 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xxxy = cbuffer.data(kg_geom_01_off + 1261 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xxxz = cbuffer.data(kg_geom_01_off + 1262 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xxyy = cbuffer.data(kg_geom_01_off + 1263 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xxyz = cbuffer.data(kg_geom_01_off + 1264 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xxzz = cbuffer.data(kg_geom_01_off + 1265 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xyyy = cbuffer.data(kg_geom_01_off + 1266 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xyyz = cbuffer.data(kg_geom_01_off + 1267 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xyzz = cbuffer.data(kg_geom_01_off + 1268 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_xzzz = cbuffer.data(kg_geom_01_off + 1269 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_yyyy = cbuffer.data(kg_geom_01_off + 1270 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_yyyz = cbuffer.data(kg_geom_01_off + 1271 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_yyzz = cbuffer.data(kg_geom_01_off + 1272 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_yzzz = cbuffer.data(kg_geom_01_off + 1273 * ccomps * dcomps);

            auto g_0_z_xxxyyzz_zzzz = cbuffer.data(kg_geom_01_off + 1274 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyyzz_xxxx, g_0_z_xxxyyzz_xxxy, g_0_z_xxxyyzz_xxxz, g_0_z_xxxyyzz_xxyy, g_0_z_xxxyyzz_xxyz, g_0_z_xxxyyzz_xxzz, g_0_z_xxxyyzz_xyyy, g_0_z_xxxyyzz_xyyz, g_0_z_xxxyyzz_xyzz, g_0_z_xxxyyzz_xzzz, g_0_z_xxxyyzz_yyyy, g_0_z_xxxyyzz_yyyz, g_0_z_xxxyyzz_yyzz, g_0_z_xxxyyzz_yzzz, g_0_z_xxxyyzz_zzzz, g_0_z_xxyyzz_xxxx, g_0_z_xxyyzz_xxxxx, g_0_z_xxyyzz_xxxxy, g_0_z_xxyyzz_xxxxz, g_0_z_xxyyzz_xxxy, g_0_z_xxyyzz_xxxyy, g_0_z_xxyyzz_xxxyz, g_0_z_xxyyzz_xxxz, g_0_z_xxyyzz_xxxzz, g_0_z_xxyyzz_xxyy, g_0_z_xxyyzz_xxyyy, g_0_z_xxyyzz_xxyyz, g_0_z_xxyyzz_xxyz, g_0_z_xxyyzz_xxyzz, g_0_z_xxyyzz_xxzz, g_0_z_xxyyzz_xxzzz, g_0_z_xxyyzz_xyyy, g_0_z_xxyyzz_xyyyy, g_0_z_xxyyzz_xyyyz, g_0_z_xxyyzz_xyyz, g_0_z_xxyyzz_xyyzz, g_0_z_xxyyzz_xyzz, g_0_z_xxyyzz_xyzzz, g_0_z_xxyyzz_xzzz, g_0_z_xxyyzz_xzzzz, g_0_z_xxyyzz_yyyy, g_0_z_xxyyzz_yyyz, g_0_z_xxyyzz_yyzz, g_0_z_xxyyzz_yzzz, g_0_z_xxyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyyzz_xxxx[k] = -g_0_z_xxyyzz_xxxx[k] * ab_x + g_0_z_xxyyzz_xxxxx[k];

                g_0_z_xxxyyzz_xxxy[k] = -g_0_z_xxyyzz_xxxy[k] * ab_x + g_0_z_xxyyzz_xxxxy[k];

                g_0_z_xxxyyzz_xxxz[k] = -g_0_z_xxyyzz_xxxz[k] * ab_x + g_0_z_xxyyzz_xxxxz[k];

                g_0_z_xxxyyzz_xxyy[k] = -g_0_z_xxyyzz_xxyy[k] * ab_x + g_0_z_xxyyzz_xxxyy[k];

                g_0_z_xxxyyzz_xxyz[k] = -g_0_z_xxyyzz_xxyz[k] * ab_x + g_0_z_xxyyzz_xxxyz[k];

                g_0_z_xxxyyzz_xxzz[k] = -g_0_z_xxyyzz_xxzz[k] * ab_x + g_0_z_xxyyzz_xxxzz[k];

                g_0_z_xxxyyzz_xyyy[k] = -g_0_z_xxyyzz_xyyy[k] * ab_x + g_0_z_xxyyzz_xxyyy[k];

                g_0_z_xxxyyzz_xyyz[k] = -g_0_z_xxyyzz_xyyz[k] * ab_x + g_0_z_xxyyzz_xxyyz[k];

                g_0_z_xxxyyzz_xyzz[k] = -g_0_z_xxyyzz_xyzz[k] * ab_x + g_0_z_xxyyzz_xxyzz[k];

                g_0_z_xxxyyzz_xzzz[k] = -g_0_z_xxyyzz_xzzz[k] * ab_x + g_0_z_xxyyzz_xxzzz[k];

                g_0_z_xxxyyzz_yyyy[k] = -g_0_z_xxyyzz_yyyy[k] * ab_x + g_0_z_xxyyzz_xyyyy[k];

                g_0_z_xxxyyzz_yyyz[k] = -g_0_z_xxyyzz_yyyz[k] * ab_x + g_0_z_xxyyzz_xyyyz[k];

                g_0_z_xxxyyzz_yyzz[k] = -g_0_z_xxyyzz_yyzz[k] * ab_x + g_0_z_xxyyzz_xyyzz[k];

                g_0_z_xxxyyzz_yzzz[k] = -g_0_z_xxyyzz_yzzz[k] * ab_x + g_0_z_xxyyzz_xyzzz[k];

                g_0_z_xxxyyzz_zzzz[k] = -g_0_z_xxyyzz_zzzz[k] * ab_x + g_0_z_xxyyzz_xzzzz[k];
            }

            /// Set up 1275-1290 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxyzzz_xxxx = cbuffer.data(kg_geom_01_off + 1275 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xxxy = cbuffer.data(kg_geom_01_off + 1276 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xxxz = cbuffer.data(kg_geom_01_off + 1277 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xxyy = cbuffer.data(kg_geom_01_off + 1278 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xxyz = cbuffer.data(kg_geom_01_off + 1279 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xxzz = cbuffer.data(kg_geom_01_off + 1280 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xyyy = cbuffer.data(kg_geom_01_off + 1281 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xyyz = cbuffer.data(kg_geom_01_off + 1282 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xyzz = cbuffer.data(kg_geom_01_off + 1283 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_xzzz = cbuffer.data(kg_geom_01_off + 1284 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_yyyy = cbuffer.data(kg_geom_01_off + 1285 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_yyyz = cbuffer.data(kg_geom_01_off + 1286 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_yyzz = cbuffer.data(kg_geom_01_off + 1287 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_yzzz = cbuffer.data(kg_geom_01_off + 1288 * ccomps * dcomps);

            auto g_0_z_xxxyzzz_zzzz = cbuffer.data(kg_geom_01_off + 1289 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxyzzz_xxxx, g_0_z_xxxyzzz_xxxy, g_0_z_xxxyzzz_xxxz, g_0_z_xxxyzzz_xxyy, g_0_z_xxxyzzz_xxyz, g_0_z_xxxyzzz_xxzz, g_0_z_xxxyzzz_xyyy, g_0_z_xxxyzzz_xyyz, g_0_z_xxxyzzz_xyzz, g_0_z_xxxyzzz_xzzz, g_0_z_xxxyzzz_yyyy, g_0_z_xxxyzzz_yyyz, g_0_z_xxxyzzz_yyzz, g_0_z_xxxyzzz_yzzz, g_0_z_xxxyzzz_zzzz, g_0_z_xxyzzz_xxxx, g_0_z_xxyzzz_xxxxx, g_0_z_xxyzzz_xxxxy, g_0_z_xxyzzz_xxxxz, g_0_z_xxyzzz_xxxy, g_0_z_xxyzzz_xxxyy, g_0_z_xxyzzz_xxxyz, g_0_z_xxyzzz_xxxz, g_0_z_xxyzzz_xxxzz, g_0_z_xxyzzz_xxyy, g_0_z_xxyzzz_xxyyy, g_0_z_xxyzzz_xxyyz, g_0_z_xxyzzz_xxyz, g_0_z_xxyzzz_xxyzz, g_0_z_xxyzzz_xxzz, g_0_z_xxyzzz_xxzzz, g_0_z_xxyzzz_xyyy, g_0_z_xxyzzz_xyyyy, g_0_z_xxyzzz_xyyyz, g_0_z_xxyzzz_xyyz, g_0_z_xxyzzz_xyyzz, g_0_z_xxyzzz_xyzz, g_0_z_xxyzzz_xyzzz, g_0_z_xxyzzz_xzzz, g_0_z_xxyzzz_xzzzz, g_0_z_xxyzzz_yyyy, g_0_z_xxyzzz_yyyz, g_0_z_xxyzzz_yyzz, g_0_z_xxyzzz_yzzz, g_0_z_xxyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxyzzz_xxxx[k] = -g_0_z_xxyzzz_xxxx[k] * ab_x + g_0_z_xxyzzz_xxxxx[k];

                g_0_z_xxxyzzz_xxxy[k] = -g_0_z_xxyzzz_xxxy[k] * ab_x + g_0_z_xxyzzz_xxxxy[k];

                g_0_z_xxxyzzz_xxxz[k] = -g_0_z_xxyzzz_xxxz[k] * ab_x + g_0_z_xxyzzz_xxxxz[k];

                g_0_z_xxxyzzz_xxyy[k] = -g_0_z_xxyzzz_xxyy[k] * ab_x + g_0_z_xxyzzz_xxxyy[k];

                g_0_z_xxxyzzz_xxyz[k] = -g_0_z_xxyzzz_xxyz[k] * ab_x + g_0_z_xxyzzz_xxxyz[k];

                g_0_z_xxxyzzz_xxzz[k] = -g_0_z_xxyzzz_xxzz[k] * ab_x + g_0_z_xxyzzz_xxxzz[k];

                g_0_z_xxxyzzz_xyyy[k] = -g_0_z_xxyzzz_xyyy[k] * ab_x + g_0_z_xxyzzz_xxyyy[k];

                g_0_z_xxxyzzz_xyyz[k] = -g_0_z_xxyzzz_xyyz[k] * ab_x + g_0_z_xxyzzz_xxyyz[k];

                g_0_z_xxxyzzz_xyzz[k] = -g_0_z_xxyzzz_xyzz[k] * ab_x + g_0_z_xxyzzz_xxyzz[k];

                g_0_z_xxxyzzz_xzzz[k] = -g_0_z_xxyzzz_xzzz[k] * ab_x + g_0_z_xxyzzz_xxzzz[k];

                g_0_z_xxxyzzz_yyyy[k] = -g_0_z_xxyzzz_yyyy[k] * ab_x + g_0_z_xxyzzz_xyyyy[k];

                g_0_z_xxxyzzz_yyyz[k] = -g_0_z_xxyzzz_yyyz[k] * ab_x + g_0_z_xxyzzz_xyyyz[k];

                g_0_z_xxxyzzz_yyzz[k] = -g_0_z_xxyzzz_yyzz[k] * ab_x + g_0_z_xxyzzz_xyyzz[k];

                g_0_z_xxxyzzz_yzzz[k] = -g_0_z_xxyzzz_yzzz[k] * ab_x + g_0_z_xxyzzz_xyzzz[k];

                g_0_z_xxxyzzz_zzzz[k] = -g_0_z_xxyzzz_zzzz[k] * ab_x + g_0_z_xxyzzz_xzzzz[k];
            }

            /// Set up 1290-1305 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxxzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1290 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1291 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1292 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1293 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1294 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1295 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1296 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1297 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1298 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1299 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1300 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1301 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1302 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1303 * ccomps * dcomps);

            auto g_0_z_xxxzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1304 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxxzzzz_xxxx, g_0_z_xxxzzzz_xxxy, g_0_z_xxxzzzz_xxxz, g_0_z_xxxzzzz_xxyy, g_0_z_xxxzzzz_xxyz, g_0_z_xxxzzzz_xxzz, g_0_z_xxxzzzz_xyyy, g_0_z_xxxzzzz_xyyz, g_0_z_xxxzzzz_xyzz, g_0_z_xxxzzzz_xzzz, g_0_z_xxxzzzz_yyyy, g_0_z_xxxzzzz_yyyz, g_0_z_xxxzzzz_yyzz, g_0_z_xxxzzzz_yzzz, g_0_z_xxxzzzz_zzzz, g_0_z_xxzzzz_xxxx, g_0_z_xxzzzz_xxxxx, g_0_z_xxzzzz_xxxxy, g_0_z_xxzzzz_xxxxz, g_0_z_xxzzzz_xxxy, g_0_z_xxzzzz_xxxyy, g_0_z_xxzzzz_xxxyz, g_0_z_xxzzzz_xxxz, g_0_z_xxzzzz_xxxzz, g_0_z_xxzzzz_xxyy, g_0_z_xxzzzz_xxyyy, g_0_z_xxzzzz_xxyyz, g_0_z_xxzzzz_xxyz, g_0_z_xxzzzz_xxyzz, g_0_z_xxzzzz_xxzz, g_0_z_xxzzzz_xxzzz, g_0_z_xxzzzz_xyyy, g_0_z_xxzzzz_xyyyy, g_0_z_xxzzzz_xyyyz, g_0_z_xxzzzz_xyyz, g_0_z_xxzzzz_xyyzz, g_0_z_xxzzzz_xyzz, g_0_z_xxzzzz_xyzzz, g_0_z_xxzzzz_xzzz, g_0_z_xxzzzz_xzzzz, g_0_z_xxzzzz_yyyy, g_0_z_xxzzzz_yyyz, g_0_z_xxzzzz_yyzz, g_0_z_xxzzzz_yzzz, g_0_z_xxzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxxzzzz_xxxx[k] = -g_0_z_xxzzzz_xxxx[k] * ab_x + g_0_z_xxzzzz_xxxxx[k];

                g_0_z_xxxzzzz_xxxy[k] = -g_0_z_xxzzzz_xxxy[k] * ab_x + g_0_z_xxzzzz_xxxxy[k];

                g_0_z_xxxzzzz_xxxz[k] = -g_0_z_xxzzzz_xxxz[k] * ab_x + g_0_z_xxzzzz_xxxxz[k];

                g_0_z_xxxzzzz_xxyy[k] = -g_0_z_xxzzzz_xxyy[k] * ab_x + g_0_z_xxzzzz_xxxyy[k];

                g_0_z_xxxzzzz_xxyz[k] = -g_0_z_xxzzzz_xxyz[k] * ab_x + g_0_z_xxzzzz_xxxyz[k];

                g_0_z_xxxzzzz_xxzz[k] = -g_0_z_xxzzzz_xxzz[k] * ab_x + g_0_z_xxzzzz_xxxzz[k];

                g_0_z_xxxzzzz_xyyy[k] = -g_0_z_xxzzzz_xyyy[k] * ab_x + g_0_z_xxzzzz_xxyyy[k];

                g_0_z_xxxzzzz_xyyz[k] = -g_0_z_xxzzzz_xyyz[k] * ab_x + g_0_z_xxzzzz_xxyyz[k];

                g_0_z_xxxzzzz_xyzz[k] = -g_0_z_xxzzzz_xyzz[k] * ab_x + g_0_z_xxzzzz_xxyzz[k];

                g_0_z_xxxzzzz_xzzz[k] = -g_0_z_xxzzzz_xzzz[k] * ab_x + g_0_z_xxzzzz_xxzzz[k];

                g_0_z_xxxzzzz_yyyy[k] = -g_0_z_xxzzzz_yyyy[k] * ab_x + g_0_z_xxzzzz_xyyyy[k];

                g_0_z_xxxzzzz_yyyz[k] = -g_0_z_xxzzzz_yyyz[k] * ab_x + g_0_z_xxzzzz_xyyyz[k];

                g_0_z_xxxzzzz_yyzz[k] = -g_0_z_xxzzzz_yyzz[k] * ab_x + g_0_z_xxzzzz_xyyzz[k];

                g_0_z_xxxzzzz_yzzz[k] = -g_0_z_xxzzzz_yzzz[k] * ab_x + g_0_z_xxzzzz_xyzzz[k];

                g_0_z_xxxzzzz_zzzz[k] = -g_0_z_xxzzzz_zzzz[k] * ab_x + g_0_z_xxzzzz_xzzzz[k];
            }

            /// Set up 1305-1320 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyyy_xxxx = cbuffer.data(kg_geom_01_off + 1305 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xxxy = cbuffer.data(kg_geom_01_off + 1306 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xxxz = cbuffer.data(kg_geom_01_off + 1307 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xxyy = cbuffer.data(kg_geom_01_off + 1308 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xxyz = cbuffer.data(kg_geom_01_off + 1309 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xxzz = cbuffer.data(kg_geom_01_off + 1310 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xyyy = cbuffer.data(kg_geom_01_off + 1311 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xyyz = cbuffer.data(kg_geom_01_off + 1312 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xyzz = cbuffer.data(kg_geom_01_off + 1313 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_xzzz = cbuffer.data(kg_geom_01_off + 1314 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_yyyy = cbuffer.data(kg_geom_01_off + 1315 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_yyyz = cbuffer.data(kg_geom_01_off + 1316 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_yyzz = cbuffer.data(kg_geom_01_off + 1317 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_yzzz = cbuffer.data(kg_geom_01_off + 1318 * ccomps * dcomps);

            auto g_0_z_xxyyyyy_zzzz = cbuffer.data(kg_geom_01_off + 1319 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyyy_xxxx, g_0_z_xxyyyyy_xxxy, g_0_z_xxyyyyy_xxxz, g_0_z_xxyyyyy_xxyy, g_0_z_xxyyyyy_xxyz, g_0_z_xxyyyyy_xxzz, g_0_z_xxyyyyy_xyyy, g_0_z_xxyyyyy_xyyz, g_0_z_xxyyyyy_xyzz, g_0_z_xxyyyyy_xzzz, g_0_z_xxyyyyy_yyyy, g_0_z_xxyyyyy_yyyz, g_0_z_xxyyyyy_yyzz, g_0_z_xxyyyyy_yzzz, g_0_z_xxyyyyy_zzzz, g_0_z_xyyyyy_xxxx, g_0_z_xyyyyy_xxxxx, g_0_z_xyyyyy_xxxxy, g_0_z_xyyyyy_xxxxz, g_0_z_xyyyyy_xxxy, g_0_z_xyyyyy_xxxyy, g_0_z_xyyyyy_xxxyz, g_0_z_xyyyyy_xxxz, g_0_z_xyyyyy_xxxzz, g_0_z_xyyyyy_xxyy, g_0_z_xyyyyy_xxyyy, g_0_z_xyyyyy_xxyyz, g_0_z_xyyyyy_xxyz, g_0_z_xyyyyy_xxyzz, g_0_z_xyyyyy_xxzz, g_0_z_xyyyyy_xxzzz, g_0_z_xyyyyy_xyyy, g_0_z_xyyyyy_xyyyy, g_0_z_xyyyyy_xyyyz, g_0_z_xyyyyy_xyyz, g_0_z_xyyyyy_xyyzz, g_0_z_xyyyyy_xyzz, g_0_z_xyyyyy_xyzzz, g_0_z_xyyyyy_xzzz, g_0_z_xyyyyy_xzzzz, g_0_z_xyyyyy_yyyy, g_0_z_xyyyyy_yyyz, g_0_z_xyyyyy_yyzz, g_0_z_xyyyyy_yzzz, g_0_z_xyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyyy_xxxx[k] = -g_0_z_xyyyyy_xxxx[k] * ab_x + g_0_z_xyyyyy_xxxxx[k];

                g_0_z_xxyyyyy_xxxy[k] = -g_0_z_xyyyyy_xxxy[k] * ab_x + g_0_z_xyyyyy_xxxxy[k];

                g_0_z_xxyyyyy_xxxz[k] = -g_0_z_xyyyyy_xxxz[k] * ab_x + g_0_z_xyyyyy_xxxxz[k];

                g_0_z_xxyyyyy_xxyy[k] = -g_0_z_xyyyyy_xxyy[k] * ab_x + g_0_z_xyyyyy_xxxyy[k];

                g_0_z_xxyyyyy_xxyz[k] = -g_0_z_xyyyyy_xxyz[k] * ab_x + g_0_z_xyyyyy_xxxyz[k];

                g_0_z_xxyyyyy_xxzz[k] = -g_0_z_xyyyyy_xxzz[k] * ab_x + g_0_z_xyyyyy_xxxzz[k];

                g_0_z_xxyyyyy_xyyy[k] = -g_0_z_xyyyyy_xyyy[k] * ab_x + g_0_z_xyyyyy_xxyyy[k];

                g_0_z_xxyyyyy_xyyz[k] = -g_0_z_xyyyyy_xyyz[k] * ab_x + g_0_z_xyyyyy_xxyyz[k];

                g_0_z_xxyyyyy_xyzz[k] = -g_0_z_xyyyyy_xyzz[k] * ab_x + g_0_z_xyyyyy_xxyzz[k];

                g_0_z_xxyyyyy_xzzz[k] = -g_0_z_xyyyyy_xzzz[k] * ab_x + g_0_z_xyyyyy_xxzzz[k];

                g_0_z_xxyyyyy_yyyy[k] = -g_0_z_xyyyyy_yyyy[k] * ab_x + g_0_z_xyyyyy_xyyyy[k];

                g_0_z_xxyyyyy_yyyz[k] = -g_0_z_xyyyyy_yyyz[k] * ab_x + g_0_z_xyyyyy_xyyyz[k];

                g_0_z_xxyyyyy_yyzz[k] = -g_0_z_xyyyyy_yyzz[k] * ab_x + g_0_z_xyyyyy_xyyzz[k];

                g_0_z_xxyyyyy_yzzz[k] = -g_0_z_xyyyyy_yzzz[k] * ab_x + g_0_z_xyyyyy_xyzzz[k];

                g_0_z_xxyyyyy_zzzz[k] = -g_0_z_xyyyyy_zzzz[k] * ab_x + g_0_z_xyyyyy_xzzzz[k];
            }

            /// Set up 1320-1335 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyyz_xxxx = cbuffer.data(kg_geom_01_off + 1320 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xxxy = cbuffer.data(kg_geom_01_off + 1321 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xxxz = cbuffer.data(kg_geom_01_off + 1322 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xxyy = cbuffer.data(kg_geom_01_off + 1323 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xxyz = cbuffer.data(kg_geom_01_off + 1324 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xxzz = cbuffer.data(kg_geom_01_off + 1325 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xyyy = cbuffer.data(kg_geom_01_off + 1326 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xyyz = cbuffer.data(kg_geom_01_off + 1327 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xyzz = cbuffer.data(kg_geom_01_off + 1328 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_xzzz = cbuffer.data(kg_geom_01_off + 1329 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_yyyy = cbuffer.data(kg_geom_01_off + 1330 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_yyyz = cbuffer.data(kg_geom_01_off + 1331 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_yyzz = cbuffer.data(kg_geom_01_off + 1332 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_yzzz = cbuffer.data(kg_geom_01_off + 1333 * ccomps * dcomps);

            auto g_0_z_xxyyyyz_zzzz = cbuffer.data(kg_geom_01_off + 1334 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyyz_xxxx, g_0_z_xxyyyyz_xxxy, g_0_z_xxyyyyz_xxxz, g_0_z_xxyyyyz_xxyy, g_0_z_xxyyyyz_xxyz, g_0_z_xxyyyyz_xxzz, g_0_z_xxyyyyz_xyyy, g_0_z_xxyyyyz_xyyz, g_0_z_xxyyyyz_xyzz, g_0_z_xxyyyyz_xzzz, g_0_z_xxyyyyz_yyyy, g_0_z_xxyyyyz_yyyz, g_0_z_xxyyyyz_yyzz, g_0_z_xxyyyyz_yzzz, g_0_z_xxyyyyz_zzzz, g_0_z_xyyyyz_xxxx, g_0_z_xyyyyz_xxxxx, g_0_z_xyyyyz_xxxxy, g_0_z_xyyyyz_xxxxz, g_0_z_xyyyyz_xxxy, g_0_z_xyyyyz_xxxyy, g_0_z_xyyyyz_xxxyz, g_0_z_xyyyyz_xxxz, g_0_z_xyyyyz_xxxzz, g_0_z_xyyyyz_xxyy, g_0_z_xyyyyz_xxyyy, g_0_z_xyyyyz_xxyyz, g_0_z_xyyyyz_xxyz, g_0_z_xyyyyz_xxyzz, g_0_z_xyyyyz_xxzz, g_0_z_xyyyyz_xxzzz, g_0_z_xyyyyz_xyyy, g_0_z_xyyyyz_xyyyy, g_0_z_xyyyyz_xyyyz, g_0_z_xyyyyz_xyyz, g_0_z_xyyyyz_xyyzz, g_0_z_xyyyyz_xyzz, g_0_z_xyyyyz_xyzzz, g_0_z_xyyyyz_xzzz, g_0_z_xyyyyz_xzzzz, g_0_z_xyyyyz_yyyy, g_0_z_xyyyyz_yyyz, g_0_z_xyyyyz_yyzz, g_0_z_xyyyyz_yzzz, g_0_z_xyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyyz_xxxx[k] = -g_0_z_xyyyyz_xxxx[k] * ab_x + g_0_z_xyyyyz_xxxxx[k];

                g_0_z_xxyyyyz_xxxy[k] = -g_0_z_xyyyyz_xxxy[k] * ab_x + g_0_z_xyyyyz_xxxxy[k];

                g_0_z_xxyyyyz_xxxz[k] = -g_0_z_xyyyyz_xxxz[k] * ab_x + g_0_z_xyyyyz_xxxxz[k];

                g_0_z_xxyyyyz_xxyy[k] = -g_0_z_xyyyyz_xxyy[k] * ab_x + g_0_z_xyyyyz_xxxyy[k];

                g_0_z_xxyyyyz_xxyz[k] = -g_0_z_xyyyyz_xxyz[k] * ab_x + g_0_z_xyyyyz_xxxyz[k];

                g_0_z_xxyyyyz_xxzz[k] = -g_0_z_xyyyyz_xxzz[k] * ab_x + g_0_z_xyyyyz_xxxzz[k];

                g_0_z_xxyyyyz_xyyy[k] = -g_0_z_xyyyyz_xyyy[k] * ab_x + g_0_z_xyyyyz_xxyyy[k];

                g_0_z_xxyyyyz_xyyz[k] = -g_0_z_xyyyyz_xyyz[k] * ab_x + g_0_z_xyyyyz_xxyyz[k];

                g_0_z_xxyyyyz_xyzz[k] = -g_0_z_xyyyyz_xyzz[k] * ab_x + g_0_z_xyyyyz_xxyzz[k];

                g_0_z_xxyyyyz_xzzz[k] = -g_0_z_xyyyyz_xzzz[k] * ab_x + g_0_z_xyyyyz_xxzzz[k];

                g_0_z_xxyyyyz_yyyy[k] = -g_0_z_xyyyyz_yyyy[k] * ab_x + g_0_z_xyyyyz_xyyyy[k];

                g_0_z_xxyyyyz_yyyz[k] = -g_0_z_xyyyyz_yyyz[k] * ab_x + g_0_z_xyyyyz_xyyyz[k];

                g_0_z_xxyyyyz_yyzz[k] = -g_0_z_xyyyyz_yyzz[k] * ab_x + g_0_z_xyyyyz_xyyzz[k];

                g_0_z_xxyyyyz_yzzz[k] = -g_0_z_xyyyyz_yzzz[k] * ab_x + g_0_z_xyyyyz_xyzzz[k];

                g_0_z_xxyyyyz_zzzz[k] = -g_0_z_xyyyyz_zzzz[k] * ab_x + g_0_z_xyyyyz_xzzzz[k];
            }

            /// Set up 1335-1350 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyyzz_xxxx = cbuffer.data(kg_geom_01_off + 1335 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xxxy = cbuffer.data(kg_geom_01_off + 1336 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xxxz = cbuffer.data(kg_geom_01_off + 1337 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xxyy = cbuffer.data(kg_geom_01_off + 1338 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xxyz = cbuffer.data(kg_geom_01_off + 1339 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xxzz = cbuffer.data(kg_geom_01_off + 1340 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xyyy = cbuffer.data(kg_geom_01_off + 1341 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xyyz = cbuffer.data(kg_geom_01_off + 1342 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xyzz = cbuffer.data(kg_geom_01_off + 1343 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_xzzz = cbuffer.data(kg_geom_01_off + 1344 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_yyyy = cbuffer.data(kg_geom_01_off + 1345 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_yyyz = cbuffer.data(kg_geom_01_off + 1346 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_yyzz = cbuffer.data(kg_geom_01_off + 1347 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_yzzz = cbuffer.data(kg_geom_01_off + 1348 * ccomps * dcomps);

            auto g_0_z_xxyyyzz_zzzz = cbuffer.data(kg_geom_01_off + 1349 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyyzz_xxxx, g_0_z_xxyyyzz_xxxy, g_0_z_xxyyyzz_xxxz, g_0_z_xxyyyzz_xxyy, g_0_z_xxyyyzz_xxyz, g_0_z_xxyyyzz_xxzz, g_0_z_xxyyyzz_xyyy, g_0_z_xxyyyzz_xyyz, g_0_z_xxyyyzz_xyzz, g_0_z_xxyyyzz_xzzz, g_0_z_xxyyyzz_yyyy, g_0_z_xxyyyzz_yyyz, g_0_z_xxyyyzz_yyzz, g_0_z_xxyyyzz_yzzz, g_0_z_xxyyyzz_zzzz, g_0_z_xyyyzz_xxxx, g_0_z_xyyyzz_xxxxx, g_0_z_xyyyzz_xxxxy, g_0_z_xyyyzz_xxxxz, g_0_z_xyyyzz_xxxy, g_0_z_xyyyzz_xxxyy, g_0_z_xyyyzz_xxxyz, g_0_z_xyyyzz_xxxz, g_0_z_xyyyzz_xxxzz, g_0_z_xyyyzz_xxyy, g_0_z_xyyyzz_xxyyy, g_0_z_xyyyzz_xxyyz, g_0_z_xyyyzz_xxyz, g_0_z_xyyyzz_xxyzz, g_0_z_xyyyzz_xxzz, g_0_z_xyyyzz_xxzzz, g_0_z_xyyyzz_xyyy, g_0_z_xyyyzz_xyyyy, g_0_z_xyyyzz_xyyyz, g_0_z_xyyyzz_xyyz, g_0_z_xyyyzz_xyyzz, g_0_z_xyyyzz_xyzz, g_0_z_xyyyzz_xyzzz, g_0_z_xyyyzz_xzzz, g_0_z_xyyyzz_xzzzz, g_0_z_xyyyzz_yyyy, g_0_z_xyyyzz_yyyz, g_0_z_xyyyzz_yyzz, g_0_z_xyyyzz_yzzz, g_0_z_xyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyyzz_xxxx[k] = -g_0_z_xyyyzz_xxxx[k] * ab_x + g_0_z_xyyyzz_xxxxx[k];

                g_0_z_xxyyyzz_xxxy[k] = -g_0_z_xyyyzz_xxxy[k] * ab_x + g_0_z_xyyyzz_xxxxy[k];

                g_0_z_xxyyyzz_xxxz[k] = -g_0_z_xyyyzz_xxxz[k] * ab_x + g_0_z_xyyyzz_xxxxz[k];

                g_0_z_xxyyyzz_xxyy[k] = -g_0_z_xyyyzz_xxyy[k] * ab_x + g_0_z_xyyyzz_xxxyy[k];

                g_0_z_xxyyyzz_xxyz[k] = -g_0_z_xyyyzz_xxyz[k] * ab_x + g_0_z_xyyyzz_xxxyz[k];

                g_0_z_xxyyyzz_xxzz[k] = -g_0_z_xyyyzz_xxzz[k] * ab_x + g_0_z_xyyyzz_xxxzz[k];

                g_0_z_xxyyyzz_xyyy[k] = -g_0_z_xyyyzz_xyyy[k] * ab_x + g_0_z_xyyyzz_xxyyy[k];

                g_0_z_xxyyyzz_xyyz[k] = -g_0_z_xyyyzz_xyyz[k] * ab_x + g_0_z_xyyyzz_xxyyz[k];

                g_0_z_xxyyyzz_xyzz[k] = -g_0_z_xyyyzz_xyzz[k] * ab_x + g_0_z_xyyyzz_xxyzz[k];

                g_0_z_xxyyyzz_xzzz[k] = -g_0_z_xyyyzz_xzzz[k] * ab_x + g_0_z_xyyyzz_xxzzz[k];

                g_0_z_xxyyyzz_yyyy[k] = -g_0_z_xyyyzz_yyyy[k] * ab_x + g_0_z_xyyyzz_xyyyy[k];

                g_0_z_xxyyyzz_yyyz[k] = -g_0_z_xyyyzz_yyyz[k] * ab_x + g_0_z_xyyyzz_xyyyz[k];

                g_0_z_xxyyyzz_yyzz[k] = -g_0_z_xyyyzz_yyzz[k] * ab_x + g_0_z_xyyyzz_xyyzz[k];

                g_0_z_xxyyyzz_yzzz[k] = -g_0_z_xyyyzz_yzzz[k] * ab_x + g_0_z_xyyyzz_xyzzz[k];

                g_0_z_xxyyyzz_zzzz[k] = -g_0_z_xyyyzz_zzzz[k] * ab_x + g_0_z_xyyyzz_xzzzz[k];
            }

            /// Set up 1350-1365 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyyzzz_xxxx = cbuffer.data(kg_geom_01_off + 1350 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xxxy = cbuffer.data(kg_geom_01_off + 1351 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xxxz = cbuffer.data(kg_geom_01_off + 1352 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xxyy = cbuffer.data(kg_geom_01_off + 1353 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xxyz = cbuffer.data(kg_geom_01_off + 1354 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xxzz = cbuffer.data(kg_geom_01_off + 1355 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xyyy = cbuffer.data(kg_geom_01_off + 1356 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xyyz = cbuffer.data(kg_geom_01_off + 1357 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xyzz = cbuffer.data(kg_geom_01_off + 1358 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_xzzz = cbuffer.data(kg_geom_01_off + 1359 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_yyyy = cbuffer.data(kg_geom_01_off + 1360 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_yyyz = cbuffer.data(kg_geom_01_off + 1361 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_yyzz = cbuffer.data(kg_geom_01_off + 1362 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_yzzz = cbuffer.data(kg_geom_01_off + 1363 * ccomps * dcomps);

            auto g_0_z_xxyyzzz_zzzz = cbuffer.data(kg_geom_01_off + 1364 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyyzzz_xxxx, g_0_z_xxyyzzz_xxxy, g_0_z_xxyyzzz_xxxz, g_0_z_xxyyzzz_xxyy, g_0_z_xxyyzzz_xxyz, g_0_z_xxyyzzz_xxzz, g_0_z_xxyyzzz_xyyy, g_0_z_xxyyzzz_xyyz, g_0_z_xxyyzzz_xyzz, g_0_z_xxyyzzz_xzzz, g_0_z_xxyyzzz_yyyy, g_0_z_xxyyzzz_yyyz, g_0_z_xxyyzzz_yyzz, g_0_z_xxyyzzz_yzzz, g_0_z_xxyyzzz_zzzz, g_0_z_xyyzzz_xxxx, g_0_z_xyyzzz_xxxxx, g_0_z_xyyzzz_xxxxy, g_0_z_xyyzzz_xxxxz, g_0_z_xyyzzz_xxxy, g_0_z_xyyzzz_xxxyy, g_0_z_xyyzzz_xxxyz, g_0_z_xyyzzz_xxxz, g_0_z_xyyzzz_xxxzz, g_0_z_xyyzzz_xxyy, g_0_z_xyyzzz_xxyyy, g_0_z_xyyzzz_xxyyz, g_0_z_xyyzzz_xxyz, g_0_z_xyyzzz_xxyzz, g_0_z_xyyzzz_xxzz, g_0_z_xyyzzz_xxzzz, g_0_z_xyyzzz_xyyy, g_0_z_xyyzzz_xyyyy, g_0_z_xyyzzz_xyyyz, g_0_z_xyyzzz_xyyz, g_0_z_xyyzzz_xyyzz, g_0_z_xyyzzz_xyzz, g_0_z_xyyzzz_xyzzz, g_0_z_xyyzzz_xzzz, g_0_z_xyyzzz_xzzzz, g_0_z_xyyzzz_yyyy, g_0_z_xyyzzz_yyyz, g_0_z_xyyzzz_yyzz, g_0_z_xyyzzz_yzzz, g_0_z_xyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyyzzz_xxxx[k] = -g_0_z_xyyzzz_xxxx[k] * ab_x + g_0_z_xyyzzz_xxxxx[k];

                g_0_z_xxyyzzz_xxxy[k] = -g_0_z_xyyzzz_xxxy[k] * ab_x + g_0_z_xyyzzz_xxxxy[k];

                g_0_z_xxyyzzz_xxxz[k] = -g_0_z_xyyzzz_xxxz[k] * ab_x + g_0_z_xyyzzz_xxxxz[k];

                g_0_z_xxyyzzz_xxyy[k] = -g_0_z_xyyzzz_xxyy[k] * ab_x + g_0_z_xyyzzz_xxxyy[k];

                g_0_z_xxyyzzz_xxyz[k] = -g_0_z_xyyzzz_xxyz[k] * ab_x + g_0_z_xyyzzz_xxxyz[k];

                g_0_z_xxyyzzz_xxzz[k] = -g_0_z_xyyzzz_xxzz[k] * ab_x + g_0_z_xyyzzz_xxxzz[k];

                g_0_z_xxyyzzz_xyyy[k] = -g_0_z_xyyzzz_xyyy[k] * ab_x + g_0_z_xyyzzz_xxyyy[k];

                g_0_z_xxyyzzz_xyyz[k] = -g_0_z_xyyzzz_xyyz[k] * ab_x + g_0_z_xyyzzz_xxyyz[k];

                g_0_z_xxyyzzz_xyzz[k] = -g_0_z_xyyzzz_xyzz[k] * ab_x + g_0_z_xyyzzz_xxyzz[k];

                g_0_z_xxyyzzz_xzzz[k] = -g_0_z_xyyzzz_xzzz[k] * ab_x + g_0_z_xyyzzz_xxzzz[k];

                g_0_z_xxyyzzz_yyyy[k] = -g_0_z_xyyzzz_yyyy[k] * ab_x + g_0_z_xyyzzz_xyyyy[k];

                g_0_z_xxyyzzz_yyyz[k] = -g_0_z_xyyzzz_yyyz[k] * ab_x + g_0_z_xyyzzz_xyyyz[k];

                g_0_z_xxyyzzz_yyzz[k] = -g_0_z_xyyzzz_yyzz[k] * ab_x + g_0_z_xyyzzz_xyyzz[k];

                g_0_z_xxyyzzz_yzzz[k] = -g_0_z_xyyzzz_yzzz[k] * ab_x + g_0_z_xyyzzz_xyzzz[k];

                g_0_z_xxyyzzz_zzzz[k] = -g_0_z_xyyzzz_zzzz[k] * ab_x + g_0_z_xyyzzz_xzzzz[k];
            }

            /// Set up 1365-1380 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxyzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1365 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1366 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1367 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1368 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1369 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1370 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1371 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1372 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1373 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1374 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1375 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1376 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1377 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1378 * ccomps * dcomps);

            auto g_0_z_xxyzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1379 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxyzzzz_xxxx, g_0_z_xxyzzzz_xxxy, g_0_z_xxyzzzz_xxxz, g_0_z_xxyzzzz_xxyy, g_0_z_xxyzzzz_xxyz, g_0_z_xxyzzzz_xxzz, g_0_z_xxyzzzz_xyyy, g_0_z_xxyzzzz_xyyz, g_0_z_xxyzzzz_xyzz, g_0_z_xxyzzzz_xzzz, g_0_z_xxyzzzz_yyyy, g_0_z_xxyzzzz_yyyz, g_0_z_xxyzzzz_yyzz, g_0_z_xxyzzzz_yzzz, g_0_z_xxyzzzz_zzzz, g_0_z_xyzzzz_xxxx, g_0_z_xyzzzz_xxxxx, g_0_z_xyzzzz_xxxxy, g_0_z_xyzzzz_xxxxz, g_0_z_xyzzzz_xxxy, g_0_z_xyzzzz_xxxyy, g_0_z_xyzzzz_xxxyz, g_0_z_xyzzzz_xxxz, g_0_z_xyzzzz_xxxzz, g_0_z_xyzzzz_xxyy, g_0_z_xyzzzz_xxyyy, g_0_z_xyzzzz_xxyyz, g_0_z_xyzzzz_xxyz, g_0_z_xyzzzz_xxyzz, g_0_z_xyzzzz_xxzz, g_0_z_xyzzzz_xxzzz, g_0_z_xyzzzz_xyyy, g_0_z_xyzzzz_xyyyy, g_0_z_xyzzzz_xyyyz, g_0_z_xyzzzz_xyyz, g_0_z_xyzzzz_xyyzz, g_0_z_xyzzzz_xyzz, g_0_z_xyzzzz_xyzzz, g_0_z_xyzzzz_xzzz, g_0_z_xyzzzz_xzzzz, g_0_z_xyzzzz_yyyy, g_0_z_xyzzzz_yyyz, g_0_z_xyzzzz_yyzz, g_0_z_xyzzzz_yzzz, g_0_z_xyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxyzzzz_xxxx[k] = -g_0_z_xyzzzz_xxxx[k] * ab_x + g_0_z_xyzzzz_xxxxx[k];

                g_0_z_xxyzzzz_xxxy[k] = -g_0_z_xyzzzz_xxxy[k] * ab_x + g_0_z_xyzzzz_xxxxy[k];

                g_0_z_xxyzzzz_xxxz[k] = -g_0_z_xyzzzz_xxxz[k] * ab_x + g_0_z_xyzzzz_xxxxz[k];

                g_0_z_xxyzzzz_xxyy[k] = -g_0_z_xyzzzz_xxyy[k] * ab_x + g_0_z_xyzzzz_xxxyy[k];

                g_0_z_xxyzzzz_xxyz[k] = -g_0_z_xyzzzz_xxyz[k] * ab_x + g_0_z_xyzzzz_xxxyz[k];

                g_0_z_xxyzzzz_xxzz[k] = -g_0_z_xyzzzz_xxzz[k] * ab_x + g_0_z_xyzzzz_xxxzz[k];

                g_0_z_xxyzzzz_xyyy[k] = -g_0_z_xyzzzz_xyyy[k] * ab_x + g_0_z_xyzzzz_xxyyy[k];

                g_0_z_xxyzzzz_xyyz[k] = -g_0_z_xyzzzz_xyyz[k] * ab_x + g_0_z_xyzzzz_xxyyz[k];

                g_0_z_xxyzzzz_xyzz[k] = -g_0_z_xyzzzz_xyzz[k] * ab_x + g_0_z_xyzzzz_xxyzz[k];

                g_0_z_xxyzzzz_xzzz[k] = -g_0_z_xyzzzz_xzzz[k] * ab_x + g_0_z_xyzzzz_xxzzz[k];

                g_0_z_xxyzzzz_yyyy[k] = -g_0_z_xyzzzz_yyyy[k] * ab_x + g_0_z_xyzzzz_xyyyy[k];

                g_0_z_xxyzzzz_yyyz[k] = -g_0_z_xyzzzz_yyyz[k] * ab_x + g_0_z_xyzzzz_xyyyz[k];

                g_0_z_xxyzzzz_yyzz[k] = -g_0_z_xyzzzz_yyzz[k] * ab_x + g_0_z_xyzzzz_xyyzz[k];

                g_0_z_xxyzzzz_yzzz[k] = -g_0_z_xyzzzz_yzzz[k] * ab_x + g_0_z_xyzzzz_xyzzz[k];

                g_0_z_xxyzzzz_zzzz[k] = -g_0_z_xyzzzz_zzzz[k] * ab_x + g_0_z_xyzzzz_xzzzz[k];
            }

            /// Set up 1380-1395 components of targeted buffer : cbuffer.data(

            auto g_0_z_xxzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1380 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1381 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1382 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1383 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1384 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1385 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1386 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1387 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1388 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1389 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1390 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1391 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1392 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1393 * ccomps * dcomps);

            auto g_0_z_xxzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1394 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xxzzzzz_xxxx, g_0_z_xxzzzzz_xxxy, g_0_z_xxzzzzz_xxxz, g_0_z_xxzzzzz_xxyy, g_0_z_xxzzzzz_xxyz, g_0_z_xxzzzzz_xxzz, g_0_z_xxzzzzz_xyyy, g_0_z_xxzzzzz_xyyz, g_0_z_xxzzzzz_xyzz, g_0_z_xxzzzzz_xzzz, g_0_z_xxzzzzz_yyyy, g_0_z_xxzzzzz_yyyz, g_0_z_xxzzzzz_yyzz, g_0_z_xxzzzzz_yzzz, g_0_z_xxzzzzz_zzzz, g_0_z_xzzzzz_xxxx, g_0_z_xzzzzz_xxxxx, g_0_z_xzzzzz_xxxxy, g_0_z_xzzzzz_xxxxz, g_0_z_xzzzzz_xxxy, g_0_z_xzzzzz_xxxyy, g_0_z_xzzzzz_xxxyz, g_0_z_xzzzzz_xxxz, g_0_z_xzzzzz_xxxzz, g_0_z_xzzzzz_xxyy, g_0_z_xzzzzz_xxyyy, g_0_z_xzzzzz_xxyyz, g_0_z_xzzzzz_xxyz, g_0_z_xzzzzz_xxyzz, g_0_z_xzzzzz_xxzz, g_0_z_xzzzzz_xxzzz, g_0_z_xzzzzz_xyyy, g_0_z_xzzzzz_xyyyy, g_0_z_xzzzzz_xyyyz, g_0_z_xzzzzz_xyyz, g_0_z_xzzzzz_xyyzz, g_0_z_xzzzzz_xyzz, g_0_z_xzzzzz_xyzzz, g_0_z_xzzzzz_xzzz, g_0_z_xzzzzz_xzzzz, g_0_z_xzzzzz_yyyy, g_0_z_xzzzzz_yyyz, g_0_z_xzzzzz_yyzz, g_0_z_xzzzzz_yzzz, g_0_z_xzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xxzzzzz_xxxx[k] = -g_0_z_xzzzzz_xxxx[k] * ab_x + g_0_z_xzzzzz_xxxxx[k];

                g_0_z_xxzzzzz_xxxy[k] = -g_0_z_xzzzzz_xxxy[k] * ab_x + g_0_z_xzzzzz_xxxxy[k];

                g_0_z_xxzzzzz_xxxz[k] = -g_0_z_xzzzzz_xxxz[k] * ab_x + g_0_z_xzzzzz_xxxxz[k];

                g_0_z_xxzzzzz_xxyy[k] = -g_0_z_xzzzzz_xxyy[k] * ab_x + g_0_z_xzzzzz_xxxyy[k];

                g_0_z_xxzzzzz_xxyz[k] = -g_0_z_xzzzzz_xxyz[k] * ab_x + g_0_z_xzzzzz_xxxyz[k];

                g_0_z_xxzzzzz_xxzz[k] = -g_0_z_xzzzzz_xxzz[k] * ab_x + g_0_z_xzzzzz_xxxzz[k];

                g_0_z_xxzzzzz_xyyy[k] = -g_0_z_xzzzzz_xyyy[k] * ab_x + g_0_z_xzzzzz_xxyyy[k];

                g_0_z_xxzzzzz_xyyz[k] = -g_0_z_xzzzzz_xyyz[k] * ab_x + g_0_z_xzzzzz_xxyyz[k];

                g_0_z_xxzzzzz_xyzz[k] = -g_0_z_xzzzzz_xyzz[k] * ab_x + g_0_z_xzzzzz_xxyzz[k];

                g_0_z_xxzzzzz_xzzz[k] = -g_0_z_xzzzzz_xzzz[k] * ab_x + g_0_z_xzzzzz_xxzzz[k];

                g_0_z_xxzzzzz_yyyy[k] = -g_0_z_xzzzzz_yyyy[k] * ab_x + g_0_z_xzzzzz_xyyyy[k];

                g_0_z_xxzzzzz_yyyz[k] = -g_0_z_xzzzzz_yyyz[k] * ab_x + g_0_z_xzzzzz_xyyyz[k];

                g_0_z_xxzzzzz_yyzz[k] = -g_0_z_xzzzzz_yyzz[k] * ab_x + g_0_z_xzzzzz_xyyzz[k];

                g_0_z_xxzzzzz_yzzz[k] = -g_0_z_xzzzzz_yzzz[k] * ab_x + g_0_z_xzzzzz_xyzzz[k];

                g_0_z_xxzzzzz_zzzz[k] = -g_0_z_xzzzzz_zzzz[k] * ab_x + g_0_z_xzzzzz_xzzzz[k];
            }

            /// Set up 1395-1410 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyyy_xxxx = cbuffer.data(kg_geom_01_off + 1395 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xxxy = cbuffer.data(kg_geom_01_off + 1396 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xxxz = cbuffer.data(kg_geom_01_off + 1397 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xxyy = cbuffer.data(kg_geom_01_off + 1398 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xxyz = cbuffer.data(kg_geom_01_off + 1399 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xxzz = cbuffer.data(kg_geom_01_off + 1400 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xyyy = cbuffer.data(kg_geom_01_off + 1401 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xyyz = cbuffer.data(kg_geom_01_off + 1402 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xyzz = cbuffer.data(kg_geom_01_off + 1403 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_xzzz = cbuffer.data(kg_geom_01_off + 1404 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_yyyy = cbuffer.data(kg_geom_01_off + 1405 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_yyyz = cbuffer.data(kg_geom_01_off + 1406 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_yyzz = cbuffer.data(kg_geom_01_off + 1407 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_yzzz = cbuffer.data(kg_geom_01_off + 1408 * ccomps * dcomps);

            auto g_0_z_xyyyyyy_zzzz = cbuffer.data(kg_geom_01_off + 1409 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyyy_xxxx, g_0_z_xyyyyyy_xxxy, g_0_z_xyyyyyy_xxxz, g_0_z_xyyyyyy_xxyy, g_0_z_xyyyyyy_xxyz, g_0_z_xyyyyyy_xxzz, g_0_z_xyyyyyy_xyyy, g_0_z_xyyyyyy_xyyz, g_0_z_xyyyyyy_xyzz, g_0_z_xyyyyyy_xzzz, g_0_z_xyyyyyy_yyyy, g_0_z_xyyyyyy_yyyz, g_0_z_xyyyyyy_yyzz, g_0_z_xyyyyyy_yzzz, g_0_z_xyyyyyy_zzzz, g_0_z_yyyyyy_xxxx, g_0_z_yyyyyy_xxxxx, g_0_z_yyyyyy_xxxxy, g_0_z_yyyyyy_xxxxz, g_0_z_yyyyyy_xxxy, g_0_z_yyyyyy_xxxyy, g_0_z_yyyyyy_xxxyz, g_0_z_yyyyyy_xxxz, g_0_z_yyyyyy_xxxzz, g_0_z_yyyyyy_xxyy, g_0_z_yyyyyy_xxyyy, g_0_z_yyyyyy_xxyyz, g_0_z_yyyyyy_xxyz, g_0_z_yyyyyy_xxyzz, g_0_z_yyyyyy_xxzz, g_0_z_yyyyyy_xxzzz, g_0_z_yyyyyy_xyyy, g_0_z_yyyyyy_xyyyy, g_0_z_yyyyyy_xyyyz, g_0_z_yyyyyy_xyyz, g_0_z_yyyyyy_xyyzz, g_0_z_yyyyyy_xyzz, g_0_z_yyyyyy_xyzzz, g_0_z_yyyyyy_xzzz, g_0_z_yyyyyy_xzzzz, g_0_z_yyyyyy_yyyy, g_0_z_yyyyyy_yyyz, g_0_z_yyyyyy_yyzz, g_0_z_yyyyyy_yzzz, g_0_z_yyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyyy_xxxx[k] = -g_0_z_yyyyyy_xxxx[k] * ab_x + g_0_z_yyyyyy_xxxxx[k];

                g_0_z_xyyyyyy_xxxy[k] = -g_0_z_yyyyyy_xxxy[k] * ab_x + g_0_z_yyyyyy_xxxxy[k];

                g_0_z_xyyyyyy_xxxz[k] = -g_0_z_yyyyyy_xxxz[k] * ab_x + g_0_z_yyyyyy_xxxxz[k];

                g_0_z_xyyyyyy_xxyy[k] = -g_0_z_yyyyyy_xxyy[k] * ab_x + g_0_z_yyyyyy_xxxyy[k];

                g_0_z_xyyyyyy_xxyz[k] = -g_0_z_yyyyyy_xxyz[k] * ab_x + g_0_z_yyyyyy_xxxyz[k];

                g_0_z_xyyyyyy_xxzz[k] = -g_0_z_yyyyyy_xxzz[k] * ab_x + g_0_z_yyyyyy_xxxzz[k];

                g_0_z_xyyyyyy_xyyy[k] = -g_0_z_yyyyyy_xyyy[k] * ab_x + g_0_z_yyyyyy_xxyyy[k];

                g_0_z_xyyyyyy_xyyz[k] = -g_0_z_yyyyyy_xyyz[k] * ab_x + g_0_z_yyyyyy_xxyyz[k];

                g_0_z_xyyyyyy_xyzz[k] = -g_0_z_yyyyyy_xyzz[k] * ab_x + g_0_z_yyyyyy_xxyzz[k];

                g_0_z_xyyyyyy_xzzz[k] = -g_0_z_yyyyyy_xzzz[k] * ab_x + g_0_z_yyyyyy_xxzzz[k];

                g_0_z_xyyyyyy_yyyy[k] = -g_0_z_yyyyyy_yyyy[k] * ab_x + g_0_z_yyyyyy_xyyyy[k];

                g_0_z_xyyyyyy_yyyz[k] = -g_0_z_yyyyyy_yyyz[k] * ab_x + g_0_z_yyyyyy_xyyyz[k];

                g_0_z_xyyyyyy_yyzz[k] = -g_0_z_yyyyyy_yyzz[k] * ab_x + g_0_z_yyyyyy_xyyzz[k];

                g_0_z_xyyyyyy_yzzz[k] = -g_0_z_yyyyyy_yzzz[k] * ab_x + g_0_z_yyyyyy_xyzzz[k];

                g_0_z_xyyyyyy_zzzz[k] = -g_0_z_yyyyyy_zzzz[k] * ab_x + g_0_z_yyyyyy_xzzzz[k];
            }

            /// Set up 1410-1425 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyyz_xxxx = cbuffer.data(kg_geom_01_off + 1410 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xxxy = cbuffer.data(kg_geom_01_off + 1411 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xxxz = cbuffer.data(kg_geom_01_off + 1412 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xxyy = cbuffer.data(kg_geom_01_off + 1413 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xxyz = cbuffer.data(kg_geom_01_off + 1414 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xxzz = cbuffer.data(kg_geom_01_off + 1415 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xyyy = cbuffer.data(kg_geom_01_off + 1416 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xyyz = cbuffer.data(kg_geom_01_off + 1417 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xyzz = cbuffer.data(kg_geom_01_off + 1418 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_xzzz = cbuffer.data(kg_geom_01_off + 1419 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_yyyy = cbuffer.data(kg_geom_01_off + 1420 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_yyyz = cbuffer.data(kg_geom_01_off + 1421 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_yyzz = cbuffer.data(kg_geom_01_off + 1422 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_yzzz = cbuffer.data(kg_geom_01_off + 1423 * ccomps * dcomps);

            auto g_0_z_xyyyyyz_zzzz = cbuffer.data(kg_geom_01_off + 1424 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyyz_xxxx, g_0_z_xyyyyyz_xxxy, g_0_z_xyyyyyz_xxxz, g_0_z_xyyyyyz_xxyy, g_0_z_xyyyyyz_xxyz, g_0_z_xyyyyyz_xxzz, g_0_z_xyyyyyz_xyyy, g_0_z_xyyyyyz_xyyz, g_0_z_xyyyyyz_xyzz, g_0_z_xyyyyyz_xzzz, g_0_z_xyyyyyz_yyyy, g_0_z_xyyyyyz_yyyz, g_0_z_xyyyyyz_yyzz, g_0_z_xyyyyyz_yzzz, g_0_z_xyyyyyz_zzzz, g_0_z_yyyyyz_xxxx, g_0_z_yyyyyz_xxxxx, g_0_z_yyyyyz_xxxxy, g_0_z_yyyyyz_xxxxz, g_0_z_yyyyyz_xxxy, g_0_z_yyyyyz_xxxyy, g_0_z_yyyyyz_xxxyz, g_0_z_yyyyyz_xxxz, g_0_z_yyyyyz_xxxzz, g_0_z_yyyyyz_xxyy, g_0_z_yyyyyz_xxyyy, g_0_z_yyyyyz_xxyyz, g_0_z_yyyyyz_xxyz, g_0_z_yyyyyz_xxyzz, g_0_z_yyyyyz_xxzz, g_0_z_yyyyyz_xxzzz, g_0_z_yyyyyz_xyyy, g_0_z_yyyyyz_xyyyy, g_0_z_yyyyyz_xyyyz, g_0_z_yyyyyz_xyyz, g_0_z_yyyyyz_xyyzz, g_0_z_yyyyyz_xyzz, g_0_z_yyyyyz_xyzzz, g_0_z_yyyyyz_xzzz, g_0_z_yyyyyz_xzzzz, g_0_z_yyyyyz_yyyy, g_0_z_yyyyyz_yyyz, g_0_z_yyyyyz_yyzz, g_0_z_yyyyyz_yzzz, g_0_z_yyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyyz_xxxx[k] = -g_0_z_yyyyyz_xxxx[k] * ab_x + g_0_z_yyyyyz_xxxxx[k];

                g_0_z_xyyyyyz_xxxy[k] = -g_0_z_yyyyyz_xxxy[k] * ab_x + g_0_z_yyyyyz_xxxxy[k];

                g_0_z_xyyyyyz_xxxz[k] = -g_0_z_yyyyyz_xxxz[k] * ab_x + g_0_z_yyyyyz_xxxxz[k];

                g_0_z_xyyyyyz_xxyy[k] = -g_0_z_yyyyyz_xxyy[k] * ab_x + g_0_z_yyyyyz_xxxyy[k];

                g_0_z_xyyyyyz_xxyz[k] = -g_0_z_yyyyyz_xxyz[k] * ab_x + g_0_z_yyyyyz_xxxyz[k];

                g_0_z_xyyyyyz_xxzz[k] = -g_0_z_yyyyyz_xxzz[k] * ab_x + g_0_z_yyyyyz_xxxzz[k];

                g_0_z_xyyyyyz_xyyy[k] = -g_0_z_yyyyyz_xyyy[k] * ab_x + g_0_z_yyyyyz_xxyyy[k];

                g_0_z_xyyyyyz_xyyz[k] = -g_0_z_yyyyyz_xyyz[k] * ab_x + g_0_z_yyyyyz_xxyyz[k];

                g_0_z_xyyyyyz_xyzz[k] = -g_0_z_yyyyyz_xyzz[k] * ab_x + g_0_z_yyyyyz_xxyzz[k];

                g_0_z_xyyyyyz_xzzz[k] = -g_0_z_yyyyyz_xzzz[k] * ab_x + g_0_z_yyyyyz_xxzzz[k];

                g_0_z_xyyyyyz_yyyy[k] = -g_0_z_yyyyyz_yyyy[k] * ab_x + g_0_z_yyyyyz_xyyyy[k];

                g_0_z_xyyyyyz_yyyz[k] = -g_0_z_yyyyyz_yyyz[k] * ab_x + g_0_z_yyyyyz_xyyyz[k];

                g_0_z_xyyyyyz_yyzz[k] = -g_0_z_yyyyyz_yyzz[k] * ab_x + g_0_z_yyyyyz_xyyzz[k];

                g_0_z_xyyyyyz_yzzz[k] = -g_0_z_yyyyyz_yzzz[k] * ab_x + g_0_z_yyyyyz_xyzzz[k];

                g_0_z_xyyyyyz_zzzz[k] = -g_0_z_yyyyyz_zzzz[k] * ab_x + g_0_z_yyyyyz_xzzzz[k];
            }

            /// Set up 1425-1440 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyyzz_xxxx = cbuffer.data(kg_geom_01_off + 1425 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xxxy = cbuffer.data(kg_geom_01_off + 1426 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xxxz = cbuffer.data(kg_geom_01_off + 1427 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xxyy = cbuffer.data(kg_geom_01_off + 1428 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xxyz = cbuffer.data(kg_geom_01_off + 1429 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xxzz = cbuffer.data(kg_geom_01_off + 1430 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xyyy = cbuffer.data(kg_geom_01_off + 1431 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xyyz = cbuffer.data(kg_geom_01_off + 1432 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xyzz = cbuffer.data(kg_geom_01_off + 1433 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_xzzz = cbuffer.data(kg_geom_01_off + 1434 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_yyyy = cbuffer.data(kg_geom_01_off + 1435 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_yyyz = cbuffer.data(kg_geom_01_off + 1436 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_yyzz = cbuffer.data(kg_geom_01_off + 1437 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_yzzz = cbuffer.data(kg_geom_01_off + 1438 * ccomps * dcomps);

            auto g_0_z_xyyyyzz_zzzz = cbuffer.data(kg_geom_01_off + 1439 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyyzz_xxxx, g_0_z_xyyyyzz_xxxy, g_0_z_xyyyyzz_xxxz, g_0_z_xyyyyzz_xxyy, g_0_z_xyyyyzz_xxyz, g_0_z_xyyyyzz_xxzz, g_0_z_xyyyyzz_xyyy, g_0_z_xyyyyzz_xyyz, g_0_z_xyyyyzz_xyzz, g_0_z_xyyyyzz_xzzz, g_0_z_xyyyyzz_yyyy, g_0_z_xyyyyzz_yyyz, g_0_z_xyyyyzz_yyzz, g_0_z_xyyyyzz_yzzz, g_0_z_xyyyyzz_zzzz, g_0_z_yyyyzz_xxxx, g_0_z_yyyyzz_xxxxx, g_0_z_yyyyzz_xxxxy, g_0_z_yyyyzz_xxxxz, g_0_z_yyyyzz_xxxy, g_0_z_yyyyzz_xxxyy, g_0_z_yyyyzz_xxxyz, g_0_z_yyyyzz_xxxz, g_0_z_yyyyzz_xxxzz, g_0_z_yyyyzz_xxyy, g_0_z_yyyyzz_xxyyy, g_0_z_yyyyzz_xxyyz, g_0_z_yyyyzz_xxyz, g_0_z_yyyyzz_xxyzz, g_0_z_yyyyzz_xxzz, g_0_z_yyyyzz_xxzzz, g_0_z_yyyyzz_xyyy, g_0_z_yyyyzz_xyyyy, g_0_z_yyyyzz_xyyyz, g_0_z_yyyyzz_xyyz, g_0_z_yyyyzz_xyyzz, g_0_z_yyyyzz_xyzz, g_0_z_yyyyzz_xyzzz, g_0_z_yyyyzz_xzzz, g_0_z_yyyyzz_xzzzz, g_0_z_yyyyzz_yyyy, g_0_z_yyyyzz_yyyz, g_0_z_yyyyzz_yyzz, g_0_z_yyyyzz_yzzz, g_0_z_yyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyyzz_xxxx[k] = -g_0_z_yyyyzz_xxxx[k] * ab_x + g_0_z_yyyyzz_xxxxx[k];

                g_0_z_xyyyyzz_xxxy[k] = -g_0_z_yyyyzz_xxxy[k] * ab_x + g_0_z_yyyyzz_xxxxy[k];

                g_0_z_xyyyyzz_xxxz[k] = -g_0_z_yyyyzz_xxxz[k] * ab_x + g_0_z_yyyyzz_xxxxz[k];

                g_0_z_xyyyyzz_xxyy[k] = -g_0_z_yyyyzz_xxyy[k] * ab_x + g_0_z_yyyyzz_xxxyy[k];

                g_0_z_xyyyyzz_xxyz[k] = -g_0_z_yyyyzz_xxyz[k] * ab_x + g_0_z_yyyyzz_xxxyz[k];

                g_0_z_xyyyyzz_xxzz[k] = -g_0_z_yyyyzz_xxzz[k] * ab_x + g_0_z_yyyyzz_xxxzz[k];

                g_0_z_xyyyyzz_xyyy[k] = -g_0_z_yyyyzz_xyyy[k] * ab_x + g_0_z_yyyyzz_xxyyy[k];

                g_0_z_xyyyyzz_xyyz[k] = -g_0_z_yyyyzz_xyyz[k] * ab_x + g_0_z_yyyyzz_xxyyz[k];

                g_0_z_xyyyyzz_xyzz[k] = -g_0_z_yyyyzz_xyzz[k] * ab_x + g_0_z_yyyyzz_xxyzz[k];

                g_0_z_xyyyyzz_xzzz[k] = -g_0_z_yyyyzz_xzzz[k] * ab_x + g_0_z_yyyyzz_xxzzz[k];

                g_0_z_xyyyyzz_yyyy[k] = -g_0_z_yyyyzz_yyyy[k] * ab_x + g_0_z_yyyyzz_xyyyy[k];

                g_0_z_xyyyyzz_yyyz[k] = -g_0_z_yyyyzz_yyyz[k] * ab_x + g_0_z_yyyyzz_xyyyz[k];

                g_0_z_xyyyyzz_yyzz[k] = -g_0_z_yyyyzz_yyzz[k] * ab_x + g_0_z_yyyyzz_xyyzz[k];

                g_0_z_xyyyyzz_yzzz[k] = -g_0_z_yyyyzz_yzzz[k] * ab_x + g_0_z_yyyyzz_xyzzz[k];

                g_0_z_xyyyyzz_zzzz[k] = -g_0_z_yyyyzz_zzzz[k] * ab_x + g_0_z_yyyyzz_xzzzz[k];
            }

            /// Set up 1440-1455 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyyzzz_xxxx = cbuffer.data(kg_geom_01_off + 1440 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xxxy = cbuffer.data(kg_geom_01_off + 1441 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xxxz = cbuffer.data(kg_geom_01_off + 1442 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xxyy = cbuffer.data(kg_geom_01_off + 1443 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xxyz = cbuffer.data(kg_geom_01_off + 1444 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xxzz = cbuffer.data(kg_geom_01_off + 1445 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xyyy = cbuffer.data(kg_geom_01_off + 1446 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xyyz = cbuffer.data(kg_geom_01_off + 1447 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xyzz = cbuffer.data(kg_geom_01_off + 1448 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_xzzz = cbuffer.data(kg_geom_01_off + 1449 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_yyyy = cbuffer.data(kg_geom_01_off + 1450 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_yyyz = cbuffer.data(kg_geom_01_off + 1451 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_yyzz = cbuffer.data(kg_geom_01_off + 1452 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_yzzz = cbuffer.data(kg_geom_01_off + 1453 * ccomps * dcomps);

            auto g_0_z_xyyyzzz_zzzz = cbuffer.data(kg_geom_01_off + 1454 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyyzzz_xxxx, g_0_z_xyyyzzz_xxxy, g_0_z_xyyyzzz_xxxz, g_0_z_xyyyzzz_xxyy, g_0_z_xyyyzzz_xxyz, g_0_z_xyyyzzz_xxzz, g_0_z_xyyyzzz_xyyy, g_0_z_xyyyzzz_xyyz, g_0_z_xyyyzzz_xyzz, g_0_z_xyyyzzz_xzzz, g_0_z_xyyyzzz_yyyy, g_0_z_xyyyzzz_yyyz, g_0_z_xyyyzzz_yyzz, g_0_z_xyyyzzz_yzzz, g_0_z_xyyyzzz_zzzz, g_0_z_yyyzzz_xxxx, g_0_z_yyyzzz_xxxxx, g_0_z_yyyzzz_xxxxy, g_0_z_yyyzzz_xxxxz, g_0_z_yyyzzz_xxxy, g_0_z_yyyzzz_xxxyy, g_0_z_yyyzzz_xxxyz, g_0_z_yyyzzz_xxxz, g_0_z_yyyzzz_xxxzz, g_0_z_yyyzzz_xxyy, g_0_z_yyyzzz_xxyyy, g_0_z_yyyzzz_xxyyz, g_0_z_yyyzzz_xxyz, g_0_z_yyyzzz_xxyzz, g_0_z_yyyzzz_xxzz, g_0_z_yyyzzz_xxzzz, g_0_z_yyyzzz_xyyy, g_0_z_yyyzzz_xyyyy, g_0_z_yyyzzz_xyyyz, g_0_z_yyyzzz_xyyz, g_0_z_yyyzzz_xyyzz, g_0_z_yyyzzz_xyzz, g_0_z_yyyzzz_xyzzz, g_0_z_yyyzzz_xzzz, g_0_z_yyyzzz_xzzzz, g_0_z_yyyzzz_yyyy, g_0_z_yyyzzz_yyyz, g_0_z_yyyzzz_yyzz, g_0_z_yyyzzz_yzzz, g_0_z_yyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyyzzz_xxxx[k] = -g_0_z_yyyzzz_xxxx[k] * ab_x + g_0_z_yyyzzz_xxxxx[k];

                g_0_z_xyyyzzz_xxxy[k] = -g_0_z_yyyzzz_xxxy[k] * ab_x + g_0_z_yyyzzz_xxxxy[k];

                g_0_z_xyyyzzz_xxxz[k] = -g_0_z_yyyzzz_xxxz[k] * ab_x + g_0_z_yyyzzz_xxxxz[k];

                g_0_z_xyyyzzz_xxyy[k] = -g_0_z_yyyzzz_xxyy[k] * ab_x + g_0_z_yyyzzz_xxxyy[k];

                g_0_z_xyyyzzz_xxyz[k] = -g_0_z_yyyzzz_xxyz[k] * ab_x + g_0_z_yyyzzz_xxxyz[k];

                g_0_z_xyyyzzz_xxzz[k] = -g_0_z_yyyzzz_xxzz[k] * ab_x + g_0_z_yyyzzz_xxxzz[k];

                g_0_z_xyyyzzz_xyyy[k] = -g_0_z_yyyzzz_xyyy[k] * ab_x + g_0_z_yyyzzz_xxyyy[k];

                g_0_z_xyyyzzz_xyyz[k] = -g_0_z_yyyzzz_xyyz[k] * ab_x + g_0_z_yyyzzz_xxyyz[k];

                g_0_z_xyyyzzz_xyzz[k] = -g_0_z_yyyzzz_xyzz[k] * ab_x + g_0_z_yyyzzz_xxyzz[k];

                g_0_z_xyyyzzz_xzzz[k] = -g_0_z_yyyzzz_xzzz[k] * ab_x + g_0_z_yyyzzz_xxzzz[k];

                g_0_z_xyyyzzz_yyyy[k] = -g_0_z_yyyzzz_yyyy[k] * ab_x + g_0_z_yyyzzz_xyyyy[k];

                g_0_z_xyyyzzz_yyyz[k] = -g_0_z_yyyzzz_yyyz[k] * ab_x + g_0_z_yyyzzz_xyyyz[k];

                g_0_z_xyyyzzz_yyzz[k] = -g_0_z_yyyzzz_yyzz[k] * ab_x + g_0_z_yyyzzz_xyyzz[k];

                g_0_z_xyyyzzz_yzzz[k] = -g_0_z_yyyzzz_yzzz[k] * ab_x + g_0_z_yyyzzz_xyzzz[k];

                g_0_z_xyyyzzz_zzzz[k] = -g_0_z_yyyzzz_zzzz[k] * ab_x + g_0_z_yyyzzz_xzzzz[k];
            }

            /// Set up 1455-1470 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyyzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1455 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1456 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1457 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1458 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1459 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1460 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1461 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1462 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1463 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1464 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1465 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1466 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1467 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1468 * ccomps * dcomps);

            auto g_0_z_xyyzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1469 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyyzzzz_xxxx, g_0_z_xyyzzzz_xxxy, g_0_z_xyyzzzz_xxxz, g_0_z_xyyzzzz_xxyy, g_0_z_xyyzzzz_xxyz, g_0_z_xyyzzzz_xxzz, g_0_z_xyyzzzz_xyyy, g_0_z_xyyzzzz_xyyz, g_0_z_xyyzzzz_xyzz, g_0_z_xyyzzzz_xzzz, g_0_z_xyyzzzz_yyyy, g_0_z_xyyzzzz_yyyz, g_0_z_xyyzzzz_yyzz, g_0_z_xyyzzzz_yzzz, g_0_z_xyyzzzz_zzzz, g_0_z_yyzzzz_xxxx, g_0_z_yyzzzz_xxxxx, g_0_z_yyzzzz_xxxxy, g_0_z_yyzzzz_xxxxz, g_0_z_yyzzzz_xxxy, g_0_z_yyzzzz_xxxyy, g_0_z_yyzzzz_xxxyz, g_0_z_yyzzzz_xxxz, g_0_z_yyzzzz_xxxzz, g_0_z_yyzzzz_xxyy, g_0_z_yyzzzz_xxyyy, g_0_z_yyzzzz_xxyyz, g_0_z_yyzzzz_xxyz, g_0_z_yyzzzz_xxyzz, g_0_z_yyzzzz_xxzz, g_0_z_yyzzzz_xxzzz, g_0_z_yyzzzz_xyyy, g_0_z_yyzzzz_xyyyy, g_0_z_yyzzzz_xyyyz, g_0_z_yyzzzz_xyyz, g_0_z_yyzzzz_xyyzz, g_0_z_yyzzzz_xyzz, g_0_z_yyzzzz_xyzzz, g_0_z_yyzzzz_xzzz, g_0_z_yyzzzz_xzzzz, g_0_z_yyzzzz_yyyy, g_0_z_yyzzzz_yyyz, g_0_z_yyzzzz_yyzz, g_0_z_yyzzzz_yzzz, g_0_z_yyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyyzzzz_xxxx[k] = -g_0_z_yyzzzz_xxxx[k] * ab_x + g_0_z_yyzzzz_xxxxx[k];

                g_0_z_xyyzzzz_xxxy[k] = -g_0_z_yyzzzz_xxxy[k] * ab_x + g_0_z_yyzzzz_xxxxy[k];

                g_0_z_xyyzzzz_xxxz[k] = -g_0_z_yyzzzz_xxxz[k] * ab_x + g_0_z_yyzzzz_xxxxz[k];

                g_0_z_xyyzzzz_xxyy[k] = -g_0_z_yyzzzz_xxyy[k] * ab_x + g_0_z_yyzzzz_xxxyy[k];

                g_0_z_xyyzzzz_xxyz[k] = -g_0_z_yyzzzz_xxyz[k] * ab_x + g_0_z_yyzzzz_xxxyz[k];

                g_0_z_xyyzzzz_xxzz[k] = -g_0_z_yyzzzz_xxzz[k] * ab_x + g_0_z_yyzzzz_xxxzz[k];

                g_0_z_xyyzzzz_xyyy[k] = -g_0_z_yyzzzz_xyyy[k] * ab_x + g_0_z_yyzzzz_xxyyy[k];

                g_0_z_xyyzzzz_xyyz[k] = -g_0_z_yyzzzz_xyyz[k] * ab_x + g_0_z_yyzzzz_xxyyz[k];

                g_0_z_xyyzzzz_xyzz[k] = -g_0_z_yyzzzz_xyzz[k] * ab_x + g_0_z_yyzzzz_xxyzz[k];

                g_0_z_xyyzzzz_xzzz[k] = -g_0_z_yyzzzz_xzzz[k] * ab_x + g_0_z_yyzzzz_xxzzz[k];

                g_0_z_xyyzzzz_yyyy[k] = -g_0_z_yyzzzz_yyyy[k] * ab_x + g_0_z_yyzzzz_xyyyy[k];

                g_0_z_xyyzzzz_yyyz[k] = -g_0_z_yyzzzz_yyyz[k] * ab_x + g_0_z_yyzzzz_xyyyz[k];

                g_0_z_xyyzzzz_yyzz[k] = -g_0_z_yyzzzz_yyzz[k] * ab_x + g_0_z_yyzzzz_xyyzz[k];

                g_0_z_xyyzzzz_yzzz[k] = -g_0_z_yyzzzz_yzzz[k] * ab_x + g_0_z_yyzzzz_xyzzz[k];

                g_0_z_xyyzzzz_zzzz[k] = -g_0_z_yyzzzz_zzzz[k] * ab_x + g_0_z_yyzzzz_xzzzz[k];
            }

            /// Set up 1470-1485 components of targeted buffer : cbuffer.data(

            auto g_0_z_xyzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1470 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1471 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1472 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1473 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1474 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1475 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1476 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1477 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1478 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1479 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1480 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1481 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1482 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1483 * ccomps * dcomps);

            auto g_0_z_xyzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1484 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xyzzzzz_xxxx, g_0_z_xyzzzzz_xxxy, g_0_z_xyzzzzz_xxxz, g_0_z_xyzzzzz_xxyy, g_0_z_xyzzzzz_xxyz, g_0_z_xyzzzzz_xxzz, g_0_z_xyzzzzz_xyyy, g_0_z_xyzzzzz_xyyz, g_0_z_xyzzzzz_xyzz, g_0_z_xyzzzzz_xzzz, g_0_z_xyzzzzz_yyyy, g_0_z_xyzzzzz_yyyz, g_0_z_xyzzzzz_yyzz, g_0_z_xyzzzzz_yzzz, g_0_z_xyzzzzz_zzzz, g_0_z_yzzzzz_xxxx, g_0_z_yzzzzz_xxxxx, g_0_z_yzzzzz_xxxxy, g_0_z_yzzzzz_xxxxz, g_0_z_yzzzzz_xxxy, g_0_z_yzzzzz_xxxyy, g_0_z_yzzzzz_xxxyz, g_0_z_yzzzzz_xxxz, g_0_z_yzzzzz_xxxzz, g_0_z_yzzzzz_xxyy, g_0_z_yzzzzz_xxyyy, g_0_z_yzzzzz_xxyyz, g_0_z_yzzzzz_xxyz, g_0_z_yzzzzz_xxyzz, g_0_z_yzzzzz_xxzz, g_0_z_yzzzzz_xxzzz, g_0_z_yzzzzz_xyyy, g_0_z_yzzzzz_xyyyy, g_0_z_yzzzzz_xyyyz, g_0_z_yzzzzz_xyyz, g_0_z_yzzzzz_xyyzz, g_0_z_yzzzzz_xyzz, g_0_z_yzzzzz_xyzzz, g_0_z_yzzzzz_xzzz, g_0_z_yzzzzz_xzzzz, g_0_z_yzzzzz_yyyy, g_0_z_yzzzzz_yyyz, g_0_z_yzzzzz_yyzz, g_0_z_yzzzzz_yzzz, g_0_z_yzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xyzzzzz_xxxx[k] = -g_0_z_yzzzzz_xxxx[k] * ab_x + g_0_z_yzzzzz_xxxxx[k];

                g_0_z_xyzzzzz_xxxy[k] = -g_0_z_yzzzzz_xxxy[k] * ab_x + g_0_z_yzzzzz_xxxxy[k];

                g_0_z_xyzzzzz_xxxz[k] = -g_0_z_yzzzzz_xxxz[k] * ab_x + g_0_z_yzzzzz_xxxxz[k];

                g_0_z_xyzzzzz_xxyy[k] = -g_0_z_yzzzzz_xxyy[k] * ab_x + g_0_z_yzzzzz_xxxyy[k];

                g_0_z_xyzzzzz_xxyz[k] = -g_0_z_yzzzzz_xxyz[k] * ab_x + g_0_z_yzzzzz_xxxyz[k];

                g_0_z_xyzzzzz_xxzz[k] = -g_0_z_yzzzzz_xxzz[k] * ab_x + g_0_z_yzzzzz_xxxzz[k];

                g_0_z_xyzzzzz_xyyy[k] = -g_0_z_yzzzzz_xyyy[k] * ab_x + g_0_z_yzzzzz_xxyyy[k];

                g_0_z_xyzzzzz_xyyz[k] = -g_0_z_yzzzzz_xyyz[k] * ab_x + g_0_z_yzzzzz_xxyyz[k];

                g_0_z_xyzzzzz_xyzz[k] = -g_0_z_yzzzzz_xyzz[k] * ab_x + g_0_z_yzzzzz_xxyzz[k];

                g_0_z_xyzzzzz_xzzz[k] = -g_0_z_yzzzzz_xzzz[k] * ab_x + g_0_z_yzzzzz_xxzzz[k];

                g_0_z_xyzzzzz_yyyy[k] = -g_0_z_yzzzzz_yyyy[k] * ab_x + g_0_z_yzzzzz_xyyyy[k];

                g_0_z_xyzzzzz_yyyz[k] = -g_0_z_yzzzzz_yyyz[k] * ab_x + g_0_z_yzzzzz_xyyyz[k];

                g_0_z_xyzzzzz_yyzz[k] = -g_0_z_yzzzzz_yyzz[k] * ab_x + g_0_z_yzzzzz_xyyzz[k];

                g_0_z_xyzzzzz_yzzz[k] = -g_0_z_yzzzzz_yzzz[k] * ab_x + g_0_z_yzzzzz_xyzzz[k];

                g_0_z_xyzzzzz_zzzz[k] = -g_0_z_yzzzzz_zzzz[k] * ab_x + g_0_z_yzzzzz_xzzzz[k];
            }

            /// Set up 1485-1500 components of targeted buffer : cbuffer.data(

            auto g_0_z_xzzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1485 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1486 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1487 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1488 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1489 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1490 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1491 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1492 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1493 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1494 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1495 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1496 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1497 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1498 * ccomps * dcomps);

            auto g_0_z_xzzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1499 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_xzzzzzz_xxxx, g_0_z_xzzzzzz_xxxy, g_0_z_xzzzzzz_xxxz, g_0_z_xzzzzzz_xxyy, g_0_z_xzzzzzz_xxyz, g_0_z_xzzzzzz_xxzz, g_0_z_xzzzzzz_xyyy, g_0_z_xzzzzzz_xyyz, g_0_z_xzzzzzz_xyzz, g_0_z_xzzzzzz_xzzz, g_0_z_xzzzzzz_yyyy, g_0_z_xzzzzzz_yyyz, g_0_z_xzzzzzz_yyzz, g_0_z_xzzzzzz_yzzz, g_0_z_xzzzzzz_zzzz, g_0_z_zzzzzz_xxxx, g_0_z_zzzzzz_xxxxx, g_0_z_zzzzzz_xxxxy, g_0_z_zzzzzz_xxxxz, g_0_z_zzzzzz_xxxy, g_0_z_zzzzzz_xxxyy, g_0_z_zzzzzz_xxxyz, g_0_z_zzzzzz_xxxz, g_0_z_zzzzzz_xxxzz, g_0_z_zzzzzz_xxyy, g_0_z_zzzzzz_xxyyy, g_0_z_zzzzzz_xxyyz, g_0_z_zzzzzz_xxyz, g_0_z_zzzzzz_xxyzz, g_0_z_zzzzzz_xxzz, g_0_z_zzzzzz_xxzzz, g_0_z_zzzzzz_xyyy, g_0_z_zzzzzz_xyyyy, g_0_z_zzzzzz_xyyyz, g_0_z_zzzzzz_xyyz, g_0_z_zzzzzz_xyyzz, g_0_z_zzzzzz_xyzz, g_0_z_zzzzzz_xyzzz, g_0_z_zzzzzz_xzzz, g_0_z_zzzzzz_xzzzz, g_0_z_zzzzzz_yyyy, g_0_z_zzzzzz_yyyz, g_0_z_zzzzzz_yyzz, g_0_z_zzzzzz_yzzz, g_0_z_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_xzzzzzz_xxxx[k] = -g_0_z_zzzzzz_xxxx[k] * ab_x + g_0_z_zzzzzz_xxxxx[k];

                g_0_z_xzzzzzz_xxxy[k] = -g_0_z_zzzzzz_xxxy[k] * ab_x + g_0_z_zzzzzz_xxxxy[k];

                g_0_z_xzzzzzz_xxxz[k] = -g_0_z_zzzzzz_xxxz[k] * ab_x + g_0_z_zzzzzz_xxxxz[k];

                g_0_z_xzzzzzz_xxyy[k] = -g_0_z_zzzzzz_xxyy[k] * ab_x + g_0_z_zzzzzz_xxxyy[k];

                g_0_z_xzzzzzz_xxyz[k] = -g_0_z_zzzzzz_xxyz[k] * ab_x + g_0_z_zzzzzz_xxxyz[k];

                g_0_z_xzzzzzz_xxzz[k] = -g_0_z_zzzzzz_xxzz[k] * ab_x + g_0_z_zzzzzz_xxxzz[k];

                g_0_z_xzzzzzz_xyyy[k] = -g_0_z_zzzzzz_xyyy[k] * ab_x + g_0_z_zzzzzz_xxyyy[k];

                g_0_z_xzzzzzz_xyyz[k] = -g_0_z_zzzzzz_xyyz[k] * ab_x + g_0_z_zzzzzz_xxyyz[k];

                g_0_z_xzzzzzz_xyzz[k] = -g_0_z_zzzzzz_xyzz[k] * ab_x + g_0_z_zzzzzz_xxyzz[k];

                g_0_z_xzzzzzz_xzzz[k] = -g_0_z_zzzzzz_xzzz[k] * ab_x + g_0_z_zzzzzz_xxzzz[k];

                g_0_z_xzzzzzz_yyyy[k] = -g_0_z_zzzzzz_yyyy[k] * ab_x + g_0_z_zzzzzz_xyyyy[k];

                g_0_z_xzzzzzz_yyyz[k] = -g_0_z_zzzzzz_yyyz[k] * ab_x + g_0_z_zzzzzz_xyyyz[k];

                g_0_z_xzzzzzz_yyzz[k] = -g_0_z_zzzzzz_yyzz[k] * ab_x + g_0_z_zzzzzz_xyyzz[k];

                g_0_z_xzzzzzz_yzzz[k] = -g_0_z_zzzzzz_yzzz[k] * ab_x + g_0_z_zzzzzz_xyzzz[k];

                g_0_z_xzzzzzz_zzzz[k] = -g_0_z_zzzzzz_zzzz[k] * ab_x + g_0_z_zzzzzz_xzzzz[k];
            }

            /// Set up 1500-1515 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyyy_xxxx = cbuffer.data(kg_geom_01_off + 1500 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xxxy = cbuffer.data(kg_geom_01_off + 1501 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xxxz = cbuffer.data(kg_geom_01_off + 1502 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xxyy = cbuffer.data(kg_geom_01_off + 1503 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xxyz = cbuffer.data(kg_geom_01_off + 1504 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xxzz = cbuffer.data(kg_geom_01_off + 1505 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xyyy = cbuffer.data(kg_geom_01_off + 1506 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xyyz = cbuffer.data(kg_geom_01_off + 1507 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xyzz = cbuffer.data(kg_geom_01_off + 1508 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_xzzz = cbuffer.data(kg_geom_01_off + 1509 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_yyyy = cbuffer.data(kg_geom_01_off + 1510 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_yyyz = cbuffer.data(kg_geom_01_off + 1511 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_yyzz = cbuffer.data(kg_geom_01_off + 1512 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_yzzz = cbuffer.data(kg_geom_01_off + 1513 * ccomps * dcomps);

            auto g_0_z_yyyyyyy_zzzz = cbuffer.data(kg_geom_01_off + 1514 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyy_xxxx, g_0_z_yyyyyy_xxxxy, g_0_z_yyyyyy_xxxy, g_0_z_yyyyyy_xxxyy, g_0_z_yyyyyy_xxxyz, g_0_z_yyyyyy_xxxz, g_0_z_yyyyyy_xxyy, g_0_z_yyyyyy_xxyyy, g_0_z_yyyyyy_xxyyz, g_0_z_yyyyyy_xxyz, g_0_z_yyyyyy_xxyzz, g_0_z_yyyyyy_xxzz, g_0_z_yyyyyy_xyyy, g_0_z_yyyyyy_xyyyy, g_0_z_yyyyyy_xyyyz, g_0_z_yyyyyy_xyyz, g_0_z_yyyyyy_xyyzz, g_0_z_yyyyyy_xyzz, g_0_z_yyyyyy_xyzzz, g_0_z_yyyyyy_xzzz, g_0_z_yyyyyy_yyyy, g_0_z_yyyyyy_yyyyy, g_0_z_yyyyyy_yyyyz, g_0_z_yyyyyy_yyyz, g_0_z_yyyyyy_yyyzz, g_0_z_yyyyyy_yyzz, g_0_z_yyyyyy_yyzzz, g_0_z_yyyyyy_yzzz, g_0_z_yyyyyy_yzzzz, g_0_z_yyyyyy_zzzz, g_0_z_yyyyyyy_xxxx, g_0_z_yyyyyyy_xxxy, g_0_z_yyyyyyy_xxxz, g_0_z_yyyyyyy_xxyy, g_0_z_yyyyyyy_xxyz, g_0_z_yyyyyyy_xxzz, g_0_z_yyyyyyy_xyyy, g_0_z_yyyyyyy_xyyz, g_0_z_yyyyyyy_xyzz, g_0_z_yyyyyyy_xzzz, g_0_z_yyyyyyy_yyyy, g_0_z_yyyyyyy_yyyz, g_0_z_yyyyyyy_yyzz, g_0_z_yyyyyyy_yzzz, g_0_z_yyyyyyy_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyyy_xxxx[k] = -g_0_z_yyyyyy_xxxx[k] * ab_y + g_0_z_yyyyyy_xxxxy[k];

                g_0_z_yyyyyyy_xxxy[k] = -g_0_z_yyyyyy_xxxy[k] * ab_y + g_0_z_yyyyyy_xxxyy[k];

                g_0_z_yyyyyyy_xxxz[k] = -g_0_z_yyyyyy_xxxz[k] * ab_y + g_0_z_yyyyyy_xxxyz[k];

                g_0_z_yyyyyyy_xxyy[k] = -g_0_z_yyyyyy_xxyy[k] * ab_y + g_0_z_yyyyyy_xxyyy[k];

                g_0_z_yyyyyyy_xxyz[k] = -g_0_z_yyyyyy_xxyz[k] * ab_y + g_0_z_yyyyyy_xxyyz[k];

                g_0_z_yyyyyyy_xxzz[k] = -g_0_z_yyyyyy_xxzz[k] * ab_y + g_0_z_yyyyyy_xxyzz[k];

                g_0_z_yyyyyyy_xyyy[k] = -g_0_z_yyyyyy_xyyy[k] * ab_y + g_0_z_yyyyyy_xyyyy[k];

                g_0_z_yyyyyyy_xyyz[k] = -g_0_z_yyyyyy_xyyz[k] * ab_y + g_0_z_yyyyyy_xyyyz[k];

                g_0_z_yyyyyyy_xyzz[k] = -g_0_z_yyyyyy_xyzz[k] * ab_y + g_0_z_yyyyyy_xyyzz[k];

                g_0_z_yyyyyyy_xzzz[k] = -g_0_z_yyyyyy_xzzz[k] * ab_y + g_0_z_yyyyyy_xyzzz[k];

                g_0_z_yyyyyyy_yyyy[k] = -g_0_z_yyyyyy_yyyy[k] * ab_y + g_0_z_yyyyyy_yyyyy[k];

                g_0_z_yyyyyyy_yyyz[k] = -g_0_z_yyyyyy_yyyz[k] * ab_y + g_0_z_yyyyyy_yyyyz[k];

                g_0_z_yyyyyyy_yyzz[k] = -g_0_z_yyyyyy_yyzz[k] * ab_y + g_0_z_yyyyyy_yyyzz[k];

                g_0_z_yyyyyyy_yzzz[k] = -g_0_z_yyyyyy_yzzz[k] * ab_y + g_0_z_yyyyyy_yyzzz[k];

                g_0_z_yyyyyyy_zzzz[k] = -g_0_z_yyyyyy_zzzz[k] * ab_y + g_0_z_yyyyyy_yzzzz[k];
            }

            /// Set up 1515-1530 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyyz_xxxx = cbuffer.data(kg_geom_01_off + 1515 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xxxy = cbuffer.data(kg_geom_01_off + 1516 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xxxz = cbuffer.data(kg_geom_01_off + 1517 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xxyy = cbuffer.data(kg_geom_01_off + 1518 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xxyz = cbuffer.data(kg_geom_01_off + 1519 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xxzz = cbuffer.data(kg_geom_01_off + 1520 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xyyy = cbuffer.data(kg_geom_01_off + 1521 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xyyz = cbuffer.data(kg_geom_01_off + 1522 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xyzz = cbuffer.data(kg_geom_01_off + 1523 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_xzzz = cbuffer.data(kg_geom_01_off + 1524 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_yyyy = cbuffer.data(kg_geom_01_off + 1525 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_yyyz = cbuffer.data(kg_geom_01_off + 1526 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_yyzz = cbuffer.data(kg_geom_01_off + 1527 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_yzzz = cbuffer.data(kg_geom_01_off + 1528 * ccomps * dcomps);

            auto g_0_z_yyyyyyz_zzzz = cbuffer.data(kg_geom_01_off + 1529 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyyz_xxxx, g_0_z_yyyyyyz_xxxy, g_0_z_yyyyyyz_xxxz, g_0_z_yyyyyyz_xxyy, g_0_z_yyyyyyz_xxyz, g_0_z_yyyyyyz_xxzz, g_0_z_yyyyyyz_xyyy, g_0_z_yyyyyyz_xyyz, g_0_z_yyyyyyz_xyzz, g_0_z_yyyyyyz_xzzz, g_0_z_yyyyyyz_yyyy, g_0_z_yyyyyyz_yyyz, g_0_z_yyyyyyz_yyzz, g_0_z_yyyyyyz_yzzz, g_0_z_yyyyyyz_zzzz, g_0_z_yyyyyz_xxxx, g_0_z_yyyyyz_xxxxy, g_0_z_yyyyyz_xxxy, g_0_z_yyyyyz_xxxyy, g_0_z_yyyyyz_xxxyz, g_0_z_yyyyyz_xxxz, g_0_z_yyyyyz_xxyy, g_0_z_yyyyyz_xxyyy, g_0_z_yyyyyz_xxyyz, g_0_z_yyyyyz_xxyz, g_0_z_yyyyyz_xxyzz, g_0_z_yyyyyz_xxzz, g_0_z_yyyyyz_xyyy, g_0_z_yyyyyz_xyyyy, g_0_z_yyyyyz_xyyyz, g_0_z_yyyyyz_xyyz, g_0_z_yyyyyz_xyyzz, g_0_z_yyyyyz_xyzz, g_0_z_yyyyyz_xyzzz, g_0_z_yyyyyz_xzzz, g_0_z_yyyyyz_yyyy, g_0_z_yyyyyz_yyyyy, g_0_z_yyyyyz_yyyyz, g_0_z_yyyyyz_yyyz, g_0_z_yyyyyz_yyyzz, g_0_z_yyyyyz_yyzz, g_0_z_yyyyyz_yyzzz, g_0_z_yyyyyz_yzzz, g_0_z_yyyyyz_yzzzz, g_0_z_yyyyyz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyyz_xxxx[k] = -g_0_z_yyyyyz_xxxx[k] * ab_y + g_0_z_yyyyyz_xxxxy[k];

                g_0_z_yyyyyyz_xxxy[k] = -g_0_z_yyyyyz_xxxy[k] * ab_y + g_0_z_yyyyyz_xxxyy[k];

                g_0_z_yyyyyyz_xxxz[k] = -g_0_z_yyyyyz_xxxz[k] * ab_y + g_0_z_yyyyyz_xxxyz[k];

                g_0_z_yyyyyyz_xxyy[k] = -g_0_z_yyyyyz_xxyy[k] * ab_y + g_0_z_yyyyyz_xxyyy[k];

                g_0_z_yyyyyyz_xxyz[k] = -g_0_z_yyyyyz_xxyz[k] * ab_y + g_0_z_yyyyyz_xxyyz[k];

                g_0_z_yyyyyyz_xxzz[k] = -g_0_z_yyyyyz_xxzz[k] * ab_y + g_0_z_yyyyyz_xxyzz[k];

                g_0_z_yyyyyyz_xyyy[k] = -g_0_z_yyyyyz_xyyy[k] * ab_y + g_0_z_yyyyyz_xyyyy[k];

                g_0_z_yyyyyyz_xyyz[k] = -g_0_z_yyyyyz_xyyz[k] * ab_y + g_0_z_yyyyyz_xyyyz[k];

                g_0_z_yyyyyyz_xyzz[k] = -g_0_z_yyyyyz_xyzz[k] * ab_y + g_0_z_yyyyyz_xyyzz[k];

                g_0_z_yyyyyyz_xzzz[k] = -g_0_z_yyyyyz_xzzz[k] * ab_y + g_0_z_yyyyyz_xyzzz[k];

                g_0_z_yyyyyyz_yyyy[k] = -g_0_z_yyyyyz_yyyy[k] * ab_y + g_0_z_yyyyyz_yyyyy[k];

                g_0_z_yyyyyyz_yyyz[k] = -g_0_z_yyyyyz_yyyz[k] * ab_y + g_0_z_yyyyyz_yyyyz[k];

                g_0_z_yyyyyyz_yyzz[k] = -g_0_z_yyyyyz_yyzz[k] * ab_y + g_0_z_yyyyyz_yyyzz[k];

                g_0_z_yyyyyyz_yzzz[k] = -g_0_z_yyyyyz_yzzz[k] * ab_y + g_0_z_yyyyyz_yyzzz[k];

                g_0_z_yyyyyyz_zzzz[k] = -g_0_z_yyyyyz_zzzz[k] * ab_y + g_0_z_yyyyyz_yzzzz[k];
            }

            /// Set up 1530-1545 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyyzz_xxxx = cbuffer.data(kg_geom_01_off + 1530 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xxxy = cbuffer.data(kg_geom_01_off + 1531 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xxxz = cbuffer.data(kg_geom_01_off + 1532 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xxyy = cbuffer.data(kg_geom_01_off + 1533 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xxyz = cbuffer.data(kg_geom_01_off + 1534 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xxzz = cbuffer.data(kg_geom_01_off + 1535 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xyyy = cbuffer.data(kg_geom_01_off + 1536 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xyyz = cbuffer.data(kg_geom_01_off + 1537 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xyzz = cbuffer.data(kg_geom_01_off + 1538 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_xzzz = cbuffer.data(kg_geom_01_off + 1539 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_yyyy = cbuffer.data(kg_geom_01_off + 1540 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_yyyz = cbuffer.data(kg_geom_01_off + 1541 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_yyzz = cbuffer.data(kg_geom_01_off + 1542 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_yzzz = cbuffer.data(kg_geom_01_off + 1543 * ccomps * dcomps);

            auto g_0_z_yyyyyzz_zzzz = cbuffer.data(kg_geom_01_off + 1544 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyyzz_xxxx, g_0_z_yyyyyzz_xxxy, g_0_z_yyyyyzz_xxxz, g_0_z_yyyyyzz_xxyy, g_0_z_yyyyyzz_xxyz, g_0_z_yyyyyzz_xxzz, g_0_z_yyyyyzz_xyyy, g_0_z_yyyyyzz_xyyz, g_0_z_yyyyyzz_xyzz, g_0_z_yyyyyzz_xzzz, g_0_z_yyyyyzz_yyyy, g_0_z_yyyyyzz_yyyz, g_0_z_yyyyyzz_yyzz, g_0_z_yyyyyzz_yzzz, g_0_z_yyyyyzz_zzzz, g_0_z_yyyyzz_xxxx, g_0_z_yyyyzz_xxxxy, g_0_z_yyyyzz_xxxy, g_0_z_yyyyzz_xxxyy, g_0_z_yyyyzz_xxxyz, g_0_z_yyyyzz_xxxz, g_0_z_yyyyzz_xxyy, g_0_z_yyyyzz_xxyyy, g_0_z_yyyyzz_xxyyz, g_0_z_yyyyzz_xxyz, g_0_z_yyyyzz_xxyzz, g_0_z_yyyyzz_xxzz, g_0_z_yyyyzz_xyyy, g_0_z_yyyyzz_xyyyy, g_0_z_yyyyzz_xyyyz, g_0_z_yyyyzz_xyyz, g_0_z_yyyyzz_xyyzz, g_0_z_yyyyzz_xyzz, g_0_z_yyyyzz_xyzzz, g_0_z_yyyyzz_xzzz, g_0_z_yyyyzz_yyyy, g_0_z_yyyyzz_yyyyy, g_0_z_yyyyzz_yyyyz, g_0_z_yyyyzz_yyyz, g_0_z_yyyyzz_yyyzz, g_0_z_yyyyzz_yyzz, g_0_z_yyyyzz_yyzzz, g_0_z_yyyyzz_yzzz, g_0_z_yyyyzz_yzzzz, g_0_z_yyyyzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyyzz_xxxx[k] = -g_0_z_yyyyzz_xxxx[k] * ab_y + g_0_z_yyyyzz_xxxxy[k];

                g_0_z_yyyyyzz_xxxy[k] = -g_0_z_yyyyzz_xxxy[k] * ab_y + g_0_z_yyyyzz_xxxyy[k];

                g_0_z_yyyyyzz_xxxz[k] = -g_0_z_yyyyzz_xxxz[k] * ab_y + g_0_z_yyyyzz_xxxyz[k];

                g_0_z_yyyyyzz_xxyy[k] = -g_0_z_yyyyzz_xxyy[k] * ab_y + g_0_z_yyyyzz_xxyyy[k];

                g_0_z_yyyyyzz_xxyz[k] = -g_0_z_yyyyzz_xxyz[k] * ab_y + g_0_z_yyyyzz_xxyyz[k];

                g_0_z_yyyyyzz_xxzz[k] = -g_0_z_yyyyzz_xxzz[k] * ab_y + g_0_z_yyyyzz_xxyzz[k];

                g_0_z_yyyyyzz_xyyy[k] = -g_0_z_yyyyzz_xyyy[k] * ab_y + g_0_z_yyyyzz_xyyyy[k];

                g_0_z_yyyyyzz_xyyz[k] = -g_0_z_yyyyzz_xyyz[k] * ab_y + g_0_z_yyyyzz_xyyyz[k];

                g_0_z_yyyyyzz_xyzz[k] = -g_0_z_yyyyzz_xyzz[k] * ab_y + g_0_z_yyyyzz_xyyzz[k];

                g_0_z_yyyyyzz_xzzz[k] = -g_0_z_yyyyzz_xzzz[k] * ab_y + g_0_z_yyyyzz_xyzzz[k];

                g_0_z_yyyyyzz_yyyy[k] = -g_0_z_yyyyzz_yyyy[k] * ab_y + g_0_z_yyyyzz_yyyyy[k];

                g_0_z_yyyyyzz_yyyz[k] = -g_0_z_yyyyzz_yyyz[k] * ab_y + g_0_z_yyyyzz_yyyyz[k];

                g_0_z_yyyyyzz_yyzz[k] = -g_0_z_yyyyzz_yyzz[k] * ab_y + g_0_z_yyyyzz_yyyzz[k];

                g_0_z_yyyyyzz_yzzz[k] = -g_0_z_yyyyzz_yzzz[k] * ab_y + g_0_z_yyyyzz_yyzzz[k];

                g_0_z_yyyyyzz_zzzz[k] = -g_0_z_yyyyzz_zzzz[k] * ab_y + g_0_z_yyyyzz_yzzzz[k];
            }

            /// Set up 1545-1560 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyyzzz_xxxx = cbuffer.data(kg_geom_01_off + 1545 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xxxy = cbuffer.data(kg_geom_01_off + 1546 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xxxz = cbuffer.data(kg_geom_01_off + 1547 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xxyy = cbuffer.data(kg_geom_01_off + 1548 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xxyz = cbuffer.data(kg_geom_01_off + 1549 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xxzz = cbuffer.data(kg_geom_01_off + 1550 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xyyy = cbuffer.data(kg_geom_01_off + 1551 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xyyz = cbuffer.data(kg_geom_01_off + 1552 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xyzz = cbuffer.data(kg_geom_01_off + 1553 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_xzzz = cbuffer.data(kg_geom_01_off + 1554 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_yyyy = cbuffer.data(kg_geom_01_off + 1555 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_yyyz = cbuffer.data(kg_geom_01_off + 1556 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_yyzz = cbuffer.data(kg_geom_01_off + 1557 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_yzzz = cbuffer.data(kg_geom_01_off + 1558 * ccomps * dcomps);

            auto g_0_z_yyyyzzz_zzzz = cbuffer.data(kg_geom_01_off + 1559 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyyzzz_xxxx, g_0_z_yyyyzzz_xxxy, g_0_z_yyyyzzz_xxxz, g_0_z_yyyyzzz_xxyy, g_0_z_yyyyzzz_xxyz, g_0_z_yyyyzzz_xxzz, g_0_z_yyyyzzz_xyyy, g_0_z_yyyyzzz_xyyz, g_0_z_yyyyzzz_xyzz, g_0_z_yyyyzzz_xzzz, g_0_z_yyyyzzz_yyyy, g_0_z_yyyyzzz_yyyz, g_0_z_yyyyzzz_yyzz, g_0_z_yyyyzzz_yzzz, g_0_z_yyyyzzz_zzzz, g_0_z_yyyzzz_xxxx, g_0_z_yyyzzz_xxxxy, g_0_z_yyyzzz_xxxy, g_0_z_yyyzzz_xxxyy, g_0_z_yyyzzz_xxxyz, g_0_z_yyyzzz_xxxz, g_0_z_yyyzzz_xxyy, g_0_z_yyyzzz_xxyyy, g_0_z_yyyzzz_xxyyz, g_0_z_yyyzzz_xxyz, g_0_z_yyyzzz_xxyzz, g_0_z_yyyzzz_xxzz, g_0_z_yyyzzz_xyyy, g_0_z_yyyzzz_xyyyy, g_0_z_yyyzzz_xyyyz, g_0_z_yyyzzz_xyyz, g_0_z_yyyzzz_xyyzz, g_0_z_yyyzzz_xyzz, g_0_z_yyyzzz_xyzzz, g_0_z_yyyzzz_xzzz, g_0_z_yyyzzz_yyyy, g_0_z_yyyzzz_yyyyy, g_0_z_yyyzzz_yyyyz, g_0_z_yyyzzz_yyyz, g_0_z_yyyzzz_yyyzz, g_0_z_yyyzzz_yyzz, g_0_z_yyyzzz_yyzzz, g_0_z_yyyzzz_yzzz, g_0_z_yyyzzz_yzzzz, g_0_z_yyyzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyyzzz_xxxx[k] = -g_0_z_yyyzzz_xxxx[k] * ab_y + g_0_z_yyyzzz_xxxxy[k];

                g_0_z_yyyyzzz_xxxy[k] = -g_0_z_yyyzzz_xxxy[k] * ab_y + g_0_z_yyyzzz_xxxyy[k];

                g_0_z_yyyyzzz_xxxz[k] = -g_0_z_yyyzzz_xxxz[k] * ab_y + g_0_z_yyyzzz_xxxyz[k];

                g_0_z_yyyyzzz_xxyy[k] = -g_0_z_yyyzzz_xxyy[k] * ab_y + g_0_z_yyyzzz_xxyyy[k];

                g_0_z_yyyyzzz_xxyz[k] = -g_0_z_yyyzzz_xxyz[k] * ab_y + g_0_z_yyyzzz_xxyyz[k];

                g_0_z_yyyyzzz_xxzz[k] = -g_0_z_yyyzzz_xxzz[k] * ab_y + g_0_z_yyyzzz_xxyzz[k];

                g_0_z_yyyyzzz_xyyy[k] = -g_0_z_yyyzzz_xyyy[k] * ab_y + g_0_z_yyyzzz_xyyyy[k];

                g_0_z_yyyyzzz_xyyz[k] = -g_0_z_yyyzzz_xyyz[k] * ab_y + g_0_z_yyyzzz_xyyyz[k];

                g_0_z_yyyyzzz_xyzz[k] = -g_0_z_yyyzzz_xyzz[k] * ab_y + g_0_z_yyyzzz_xyyzz[k];

                g_0_z_yyyyzzz_xzzz[k] = -g_0_z_yyyzzz_xzzz[k] * ab_y + g_0_z_yyyzzz_xyzzz[k];

                g_0_z_yyyyzzz_yyyy[k] = -g_0_z_yyyzzz_yyyy[k] * ab_y + g_0_z_yyyzzz_yyyyy[k];

                g_0_z_yyyyzzz_yyyz[k] = -g_0_z_yyyzzz_yyyz[k] * ab_y + g_0_z_yyyzzz_yyyyz[k];

                g_0_z_yyyyzzz_yyzz[k] = -g_0_z_yyyzzz_yyzz[k] * ab_y + g_0_z_yyyzzz_yyyzz[k];

                g_0_z_yyyyzzz_yzzz[k] = -g_0_z_yyyzzz_yzzz[k] * ab_y + g_0_z_yyyzzz_yyzzz[k];

                g_0_z_yyyyzzz_zzzz[k] = -g_0_z_yyyzzz_zzzz[k] * ab_y + g_0_z_yyyzzz_yzzzz[k];
            }

            /// Set up 1560-1575 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyyzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1560 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1561 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1562 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1563 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1564 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1565 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1566 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1567 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1568 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1569 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1570 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1571 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1572 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1573 * ccomps * dcomps);

            auto g_0_z_yyyzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1574 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyyzzzz_xxxx, g_0_z_yyyzzzz_xxxy, g_0_z_yyyzzzz_xxxz, g_0_z_yyyzzzz_xxyy, g_0_z_yyyzzzz_xxyz, g_0_z_yyyzzzz_xxzz, g_0_z_yyyzzzz_xyyy, g_0_z_yyyzzzz_xyyz, g_0_z_yyyzzzz_xyzz, g_0_z_yyyzzzz_xzzz, g_0_z_yyyzzzz_yyyy, g_0_z_yyyzzzz_yyyz, g_0_z_yyyzzzz_yyzz, g_0_z_yyyzzzz_yzzz, g_0_z_yyyzzzz_zzzz, g_0_z_yyzzzz_xxxx, g_0_z_yyzzzz_xxxxy, g_0_z_yyzzzz_xxxy, g_0_z_yyzzzz_xxxyy, g_0_z_yyzzzz_xxxyz, g_0_z_yyzzzz_xxxz, g_0_z_yyzzzz_xxyy, g_0_z_yyzzzz_xxyyy, g_0_z_yyzzzz_xxyyz, g_0_z_yyzzzz_xxyz, g_0_z_yyzzzz_xxyzz, g_0_z_yyzzzz_xxzz, g_0_z_yyzzzz_xyyy, g_0_z_yyzzzz_xyyyy, g_0_z_yyzzzz_xyyyz, g_0_z_yyzzzz_xyyz, g_0_z_yyzzzz_xyyzz, g_0_z_yyzzzz_xyzz, g_0_z_yyzzzz_xyzzz, g_0_z_yyzzzz_xzzz, g_0_z_yyzzzz_yyyy, g_0_z_yyzzzz_yyyyy, g_0_z_yyzzzz_yyyyz, g_0_z_yyzzzz_yyyz, g_0_z_yyzzzz_yyyzz, g_0_z_yyzzzz_yyzz, g_0_z_yyzzzz_yyzzz, g_0_z_yyzzzz_yzzz, g_0_z_yyzzzz_yzzzz, g_0_z_yyzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyyzzzz_xxxx[k] = -g_0_z_yyzzzz_xxxx[k] * ab_y + g_0_z_yyzzzz_xxxxy[k];

                g_0_z_yyyzzzz_xxxy[k] = -g_0_z_yyzzzz_xxxy[k] * ab_y + g_0_z_yyzzzz_xxxyy[k];

                g_0_z_yyyzzzz_xxxz[k] = -g_0_z_yyzzzz_xxxz[k] * ab_y + g_0_z_yyzzzz_xxxyz[k];

                g_0_z_yyyzzzz_xxyy[k] = -g_0_z_yyzzzz_xxyy[k] * ab_y + g_0_z_yyzzzz_xxyyy[k];

                g_0_z_yyyzzzz_xxyz[k] = -g_0_z_yyzzzz_xxyz[k] * ab_y + g_0_z_yyzzzz_xxyyz[k];

                g_0_z_yyyzzzz_xxzz[k] = -g_0_z_yyzzzz_xxzz[k] * ab_y + g_0_z_yyzzzz_xxyzz[k];

                g_0_z_yyyzzzz_xyyy[k] = -g_0_z_yyzzzz_xyyy[k] * ab_y + g_0_z_yyzzzz_xyyyy[k];

                g_0_z_yyyzzzz_xyyz[k] = -g_0_z_yyzzzz_xyyz[k] * ab_y + g_0_z_yyzzzz_xyyyz[k];

                g_0_z_yyyzzzz_xyzz[k] = -g_0_z_yyzzzz_xyzz[k] * ab_y + g_0_z_yyzzzz_xyyzz[k];

                g_0_z_yyyzzzz_xzzz[k] = -g_0_z_yyzzzz_xzzz[k] * ab_y + g_0_z_yyzzzz_xyzzz[k];

                g_0_z_yyyzzzz_yyyy[k] = -g_0_z_yyzzzz_yyyy[k] * ab_y + g_0_z_yyzzzz_yyyyy[k];

                g_0_z_yyyzzzz_yyyz[k] = -g_0_z_yyzzzz_yyyz[k] * ab_y + g_0_z_yyzzzz_yyyyz[k];

                g_0_z_yyyzzzz_yyzz[k] = -g_0_z_yyzzzz_yyzz[k] * ab_y + g_0_z_yyzzzz_yyyzz[k];

                g_0_z_yyyzzzz_yzzz[k] = -g_0_z_yyzzzz_yzzz[k] * ab_y + g_0_z_yyzzzz_yyzzz[k];

                g_0_z_yyyzzzz_zzzz[k] = -g_0_z_yyzzzz_zzzz[k] * ab_y + g_0_z_yyzzzz_yzzzz[k];
            }

            /// Set up 1575-1590 components of targeted buffer : cbuffer.data(

            auto g_0_z_yyzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1575 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1576 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1577 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1578 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1579 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1580 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1581 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1582 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1583 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1584 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1585 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1586 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1587 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1588 * ccomps * dcomps);

            auto g_0_z_yyzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1589 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yyzzzzz_xxxx, g_0_z_yyzzzzz_xxxy, g_0_z_yyzzzzz_xxxz, g_0_z_yyzzzzz_xxyy, g_0_z_yyzzzzz_xxyz, g_0_z_yyzzzzz_xxzz, g_0_z_yyzzzzz_xyyy, g_0_z_yyzzzzz_xyyz, g_0_z_yyzzzzz_xyzz, g_0_z_yyzzzzz_xzzz, g_0_z_yyzzzzz_yyyy, g_0_z_yyzzzzz_yyyz, g_0_z_yyzzzzz_yyzz, g_0_z_yyzzzzz_yzzz, g_0_z_yyzzzzz_zzzz, g_0_z_yzzzzz_xxxx, g_0_z_yzzzzz_xxxxy, g_0_z_yzzzzz_xxxy, g_0_z_yzzzzz_xxxyy, g_0_z_yzzzzz_xxxyz, g_0_z_yzzzzz_xxxz, g_0_z_yzzzzz_xxyy, g_0_z_yzzzzz_xxyyy, g_0_z_yzzzzz_xxyyz, g_0_z_yzzzzz_xxyz, g_0_z_yzzzzz_xxyzz, g_0_z_yzzzzz_xxzz, g_0_z_yzzzzz_xyyy, g_0_z_yzzzzz_xyyyy, g_0_z_yzzzzz_xyyyz, g_0_z_yzzzzz_xyyz, g_0_z_yzzzzz_xyyzz, g_0_z_yzzzzz_xyzz, g_0_z_yzzzzz_xyzzz, g_0_z_yzzzzz_xzzz, g_0_z_yzzzzz_yyyy, g_0_z_yzzzzz_yyyyy, g_0_z_yzzzzz_yyyyz, g_0_z_yzzzzz_yyyz, g_0_z_yzzzzz_yyyzz, g_0_z_yzzzzz_yyzz, g_0_z_yzzzzz_yyzzz, g_0_z_yzzzzz_yzzz, g_0_z_yzzzzz_yzzzz, g_0_z_yzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yyzzzzz_xxxx[k] = -g_0_z_yzzzzz_xxxx[k] * ab_y + g_0_z_yzzzzz_xxxxy[k];

                g_0_z_yyzzzzz_xxxy[k] = -g_0_z_yzzzzz_xxxy[k] * ab_y + g_0_z_yzzzzz_xxxyy[k];

                g_0_z_yyzzzzz_xxxz[k] = -g_0_z_yzzzzz_xxxz[k] * ab_y + g_0_z_yzzzzz_xxxyz[k];

                g_0_z_yyzzzzz_xxyy[k] = -g_0_z_yzzzzz_xxyy[k] * ab_y + g_0_z_yzzzzz_xxyyy[k];

                g_0_z_yyzzzzz_xxyz[k] = -g_0_z_yzzzzz_xxyz[k] * ab_y + g_0_z_yzzzzz_xxyyz[k];

                g_0_z_yyzzzzz_xxzz[k] = -g_0_z_yzzzzz_xxzz[k] * ab_y + g_0_z_yzzzzz_xxyzz[k];

                g_0_z_yyzzzzz_xyyy[k] = -g_0_z_yzzzzz_xyyy[k] * ab_y + g_0_z_yzzzzz_xyyyy[k];

                g_0_z_yyzzzzz_xyyz[k] = -g_0_z_yzzzzz_xyyz[k] * ab_y + g_0_z_yzzzzz_xyyyz[k];

                g_0_z_yyzzzzz_xyzz[k] = -g_0_z_yzzzzz_xyzz[k] * ab_y + g_0_z_yzzzzz_xyyzz[k];

                g_0_z_yyzzzzz_xzzz[k] = -g_0_z_yzzzzz_xzzz[k] * ab_y + g_0_z_yzzzzz_xyzzz[k];

                g_0_z_yyzzzzz_yyyy[k] = -g_0_z_yzzzzz_yyyy[k] * ab_y + g_0_z_yzzzzz_yyyyy[k];

                g_0_z_yyzzzzz_yyyz[k] = -g_0_z_yzzzzz_yyyz[k] * ab_y + g_0_z_yzzzzz_yyyyz[k];

                g_0_z_yyzzzzz_yyzz[k] = -g_0_z_yzzzzz_yyzz[k] * ab_y + g_0_z_yzzzzz_yyyzz[k];

                g_0_z_yyzzzzz_yzzz[k] = -g_0_z_yzzzzz_yzzz[k] * ab_y + g_0_z_yzzzzz_yyzzz[k];

                g_0_z_yyzzzzz_zzzz[k] = -g_0_z_yzzzzz_zzzz[k] * ab_y + g_0_z_yzzzzz_yzzzz[k];
            }

            /// Set up 1590-1605 components of targeted buffer : cbuffer.data(

            auto g_0_z_yzzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1590 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1591 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1592 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1593 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1594 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1595 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1596 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1597 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1598 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1599 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1600 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1601 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1602 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1603 * ccomps * dcomps);

            auto g_0_z_yzzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1604 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_yzzzzzz_xxxx, g_0_z_yzzzzzz_xxxy, g_0_z_yzzzzzz_xxxz, g_0_z_yzzzzzz_xxyy, g_0_z_yzzzzzz_xxyz, g_0_z_yzzzzzz_xxzz, g_0_z_yzzzzzz_xyyy, g_0_z_yzzzzzz_xyyz, g_0_z_yzzzzzz_xyzz, g_0_z_yzzzzzz_xzzz, g_0_z_yzzzzzz_yyyy, g_0_z_yzzzzzz_yyyz, g_0_z_yzzzzzz_yyzz, g_0_z_yzzzzzz_yzzz, g_0_z_yzzzzzz_zzzz, g_0_z_zzzzzz_xxxx, g_0_z_zzzzzz_xxxxy, g_0_z_zzzzzz_xxxy, g_0_z_zzzzzz_xxxyy, g_0_z_zzzzzz_xxxyz, g_0_z_zzzzzz_xxxz, g_0_z_zzzzzz_xxyy, g_0_z_zzzzzz_xxyyy, g_0_z_zzzzzz_xxyyz, g_0_z_zzzzzz_xxyz, g_0_z_zzzzzz_xxyzz, g_0_z_zzzzzz_xxzz, g_0_z_zzzzzz_xyyy, g_0_z_zzzzzz_xyyyy, g_0_z_zzzzzz_xyyyz, g_0_z_zzzzzz_xyyz, g_0_z_zzzzzz_xyyzz, g_0_z_zzzzzz_xyzz, g_0_z_zzzzzz_xyzzz, g_0_z_zzzzzz_xzzz, g_0_z_zzzzzz_yyyy, g_0_z_zzzzzz_yyyyy, g_0_z_zzzzzz_yyyyz, g_0_z_zzzzzz_yyyz, g_0_z_zzzzzz_yyyzz, g_0_z_zzzzzz_yyzz, g_0_z_zzzzzz_yyzzz, g_0_z_zzzzzz_yzzz, g_0_z_zzzzzz_yzzzz, g_0_z_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_yzzzzzz_xxxx[k] = -g_0_z_zzzzzz_xxxx[k] * ab_y + g_0_z_zzzzzz_xxxxy[k];

                g_0_z_yzzzzzz_xxxy[k] = -g_0_z_zzzzzz_xxxy[k] * ab_y + g_0_z_zzzzzz_xxxyy[k];

                g_0_z_yzzzzzz_xxxz[k] = -g_0_z_zzzzzz_xxxz[k] * ab_y + g_0_z_zzzzzz_xxxyz[k];

                g_0_z_yzzzzzz_xxyy[k] = -g_0_z_zzzzzz_xxyy[k] * ab_y + g_0_z_zzzzzz_xxyyy[k];

                g_0_z_yzzzzzz_xxyz[k] = -g_0_z_zzzzzz_xxyz[k] * ab_y + g_0_z_zzzzzz_xxyyz[k];

                g_0_z_yzzzzzz_xxzz[k] = -g_0_z_zzzzzz_xxzz[k] * ab_y + g_0_z_zzzzzz_xxyzz[k];

                g_0_z_yzzzzzz_xyyy[k] = -g_0_z_zzzzzz_xyyy[k] * ab_y + g_0_z_zzzzzz_xyyyy[k];

                g_0_z_yzzzzzz_xyyz[k] = -g_0_z_zzzzzz_xyyz[k] * ab_y + g_0_z_zzzzzz_xyyyz[k];

                g_0_z_yzzzzzz_xyzz[k] = -g_0_z_zzzzzz_xyzz[k] * ab_y + g_0_z_zzzzzz_xyyzz[k];

                g_0_z_yzzzzzz_xzzz[k] = -g_0_z_zzzzzz_xzzz[k] * ab_y + g_0_z_zzzzzz_xyzzz[k];

                g_0_z_yzzzzzz_yyyy[k] = -g_0_z_zzzzzz_yyyy[k] * ab_y + g_0_z_zzzzzz_yyyyy[k];

                g_0_z_yzzzzzz_yyyz[k] = -g_0_z_zzzzzz_yyyz[k] * ab_y + g_0_z_zzzzzz_yyyyz[k];

                g_0_z_yzzzzzz_yyzz[k] = -g_0_z_zzzzzz_yyzz[k] * ab_y + g_0_z_zzzzzz_yyyzz[k];

                g_0_z_yzzzzzz_yzzz[k] = -g_0_z_zzzzzz_yzzz[k] * ab_y + g_0_z_zzzzzz_yyzzz[k];

                g_0_z_yzzzzzz_zzzz[k] = -g_0_z_zzzzzz_zzzz[k] * ab_y + g_0_z_zzzzzz_yzzzz[k];
            }

            /// Set up 1605-1620 components of targeted buffer : cbuffer.data(

            auto g_0_z_zzzzzzz_xxxx = cbuffer.data(kg_geom_01_off + 1605 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xxxy = cbuffer.data(kg_geom_01_off + 1606 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xxxz = cbuffer.data(kg_geom_01_off + 1607 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xxyy = cbuffer.data(kg_geom_01_off + 1608 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xxyz = cbuffer.data(kg_geom_01_off + 1609 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xxzz = cbuffer.data(kg_geom_01_off + 1610 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xyyy = cbuffer.data(kg_geom_01_off + 1611 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xyyz = cbuffer.data(kg_geom_01_off + 1612 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xyzz = cbuffer.data(kg_geom_01_off + 1613 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_xzzz = cbuffer.data(kg_geom_01_off + 1614 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_yyyy = cbuffer.data(kg_geom_01_off + 1615 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_yyyz = cbuffer.data(kg_geom_01_off + 1616 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_yyzz = cbuffer.data(kg_geom_01_off + 1617 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_yzzz = cbuffer.data(kg_geom_01_off + 1618 * ccomps * dcomps);

            auto g_0_z_zzzzzzz_zzzz = cbuffer.data(kg_geom_01_off + 1619 * ccomps * dcomps);

            #pragma omp simd aligned(g_0_z_zzzzzz_xxxx, g_0_z_zzzzzz_xxxxz, g_0_z_zzzzzz_xxxy, g_0_z_zzzzzz_xxxyz, g_0_z_zzzzzz_xxxz, g_0_z_zzzzzz_xxxzz, g_0_z_zzzzzz_xxyy, g_0_z_zzzzzz_xxyyz, g_0_z_zzzzzz_xxyz, g_0_z_zzzzzz_xxyzz, g_0_z_zzzzzz_xxzz, g_0_z_zzzzzz_xxzzz, g_0_z_zzzzzz_xyyy, g_0_z_zzzzzz_xyyyz, g_0_z_zzzzzz_xyyz, g_0_z_zzzzzz_xyyzz, g_0_z_zzzzzz_xyzz, g_0_z_zzzzzz_xyzzz, g_0_z_zzzzzz_xzzz, g_0_z_zzzzzz_xzzzz, g_0_z_zzzzzz_yyyy, g_0_z_zzzzzz_yyyyz, g_0_z_zzzzzz_yyyz, g_0_z_zzzzzz_yyyzz, g_0_z_zzzzzz_yyzz, g_0_z_zzzzzz_yyzzz, g_0_z_zzzzzz_yzzz, g_0_z_zzzzzz_yzzzz, g_0_z_zzzzzz_zzzz, g_0_z_zzzzzz_zzzzz, g_0_z_zzzzzzz_xxxx, g_0_z_zzzzzzz_xxxy, g_0_z_zzzzzzz_xxxz, g_0_z_zzzzzzz_xxyy, g_0_z_zzzzzzz_xxyz, g_0_z_zzzzzzz_xxzz, g_0_z_zzzzzzz_xyyy, g_0_z_zzzzzzz_xyyz, g_0_z_zzzzzzz_xyzz, g_0_z_zzzzzzz_xzzz, g_0_z_zzzzzzz_yyyy, g_0_z_zzzzzzz_yyyz, g_0_z_zzzzzzz_yyzz, g_0_z_zzzzzzz_yzzz, g_0_z_zzzzzzz_zzzz, g_zzzzzz_xxxx, g_zzzzzz_xxxy, g_zzzzzz_xxxz, g_zzzzzz_xxyy, g_zzzzzz_xxyz, g_zzzzzz_xxzz, g_zzzzzz_xyyy, g_zzzzzz_xyyz, g_zzzzzz_xyzz, g_zzzzzz_xzzz, g_zzzzzz_yyyy, g_zzzzzz_yyyz, g_zzzzzz_yyzz, g_zzzzzz_yzzz, g_zzzzzz_zzzz  : 64)
            for (size_t k = 0; k < nelems; k++)
            {
                g_0_z_zzzzzzz_xxxx[k] = g_zzzzzz_xxxx[k] - g_0_z_zzzzzz_xxxx[k] * ab_z + g_0_z_zzzzzz_xxxxz[k];

                g_0_z_zzzzzzz_xxxy[k] = g_zzzzzz_xxxy[k] - g_0_z_zzzzzz_xxxy[k] * ab_z + g_0_z_zzzzzz_xxxyz[k];

                g_0_z_zzzzzzz_xxxz[k] = g_zzzzzz_xxxz[k] - g_0_z_zzzzzz_xxxz[k] * ab_z + g_0_z_zzzzzz_xxxzz[k];

                g_0_z_zzzzzzz_xxyy[k] = g_zzzzzz_xxyy[k] - g_0_z_zzzzzz_xxyy[k] * ab_z + g_0_z_zzzzzz_xxyyz[k];

                g_0_z_zzzzzzz_xxyz[k] = g_zzzzzz_xxyz[k] - g_0_z_zzzzzz_xxyz[k] * ab_z + g_0_z_zzzzzz_xxyzz[k];

                g_0_z_zzzzzzz_xxzz[k] = g_zzzzzz_xxzz[k] - g_0_z_zzzzzz_xxzz[k] * ab_z + g_0_z_zzzzzz_xxzzz[k];

                g_0_z_zzzzzzz_xyyy[k] = g_zzzzzz_xyyy[k] - g_0_z_zzzzzz_xyyy[k] * ab_z + g_0_z_zzzzzz_xyyyz[k];

                g_0_z_zzzzzzz_xyyz[k] = g_zzzzzz_xyyz[k] - g_0_z_zzzzzz_xyyz[k] * ab_z + g_0_z_zzzzzz_xyyzz[k];

                g_0_z_zzzzzzz_xyzz[k] = g_zzzzzz_xyzz[k] - g_0_z_zzzzzz_xyzz[k] * ab_z + g_0_z_zzzzzz_xyzzz[k];

                g_0_z_zzzzzzz_xzzz[k] = g_zzzzzz_xzzz[k] - g_0_z_zzzzzz_xzzz[k] * ab_z + g_0_z_zzzzzz_xzzzz[k];

                g_0_z_zzzzzzz_yyyy[k] = g_zzzzzz_yyyy[k] - g_0_z_zzzzzz_yyyy[k] * ab_z + g_0_z_zzzzzz_yyyyz[k];

                g_0_z_zzzzzzz_yyyz[k] = g_zzzzzz_yyyz[k] - g_0_z_zzzzzz_yyyz[k] * ab_z + g_0_z_zzzzzz_yyyzz[k];

                g_0_z_zzzzzzz_yyzz[k] = g_zzzzzz_yyzz[k] - g_0_z_zzzzzz_yyzz[k] * ab_z + g_0_z_zzzzzz_yyzzz[k];

                g_0_z_zzzzzzz_yzzz[k] = g_zzzzzz_yzzz[k] - g_0_z_zzzzzz_yzzz[k] * ab_z + g_0_z_zzzzzz_yzzzz[k];

                g_0_z_zzzzzzz_zzzz[k] = g_zzzzzz_zzzz[k] - g_0_z_zzzzzz_zzzz[k] * ab_z + g_0_z_zzzzzz_zzzzz[k];
            }
        }
    }
}

} // erirec namespace

